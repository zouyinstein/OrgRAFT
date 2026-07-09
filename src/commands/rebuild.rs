use std::cmp::Ordering;
use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet};
use std::fs::{self, File};
use std::io::{self, Write};
use std::path::{Path, PathBuf};
use std::process::{self, Command};
use std::time::{Instant, SystemTime, UNIX_EPOCH};

use crate::commands::shared::{
    print_contract, resolve_gfa_editor_cli, run_gfa_editor_image, CommandContract, GfaImageExport,
};
use crate::error::OrgraftError;

const HELP: &str = r#"orgraft rebuild

Rebuild a compact verified GFA from an edited GFA and polished FASTA.

Usage:
  orgraft rebuild --organelle NAME --subgraph ID --edited-gfa FILE --polished-fasta FILE [options]

Inputs:
  --organelle NAME              organelle name for this rebuild run [mito]
  --subgraph ID                 subgraph/ring id [subgraph_001]
  --edited-gfa FILE             edited/check graph GFA
  --polished-fasta FILE         final polished FASTA
  --soft-paths FILE             tool paths file [soft_paths.txt]

Outputs:
  --out-dir DIR                 rebuild root output directory [results/rebuild]
  --force                       replace existing output directory

Additional Parameters:
  --threads N                   minimap2 threads for node-to-polished projection [4]
  --image-reference-fasta FILE  reference FASTA for graph colouring; enables PDF/SVG export

Layout: OUT/SUBGRAPH/rebuild_SUBGRAPH* plus OUT/logs/*.tsv (pdf/svg need reference)
"#;

const DEFAULT_SOFT_PATHS: &str = "soft_paths.txt";
const DEFAULT_THREADS: usize = 4;
const DEFAULT_ORGANELLE: &str = "mito";
const DEFAULT_SUBGRAPH: &str = "subgraph_001";
const DEFAULT_OUT_DIR: &str = "results/rebuild";

pub fn run(args: &[String]) -> Result<(), OrgraftError> {
    if args.is_empty() {
        println!("{HELP}");
        return Ok(());
    }
    if args.iter().any(|arg| arg == "-h" || arg == "--help") {
        println!("{HELP}");
        return Ok(());
    }
    if args.iter().any(|arg| arg == "--contract") {
        print_contract(&contract());
        return Ok(());
    }
    let options = RebuildOptions::from_args(args)?;
    run_rebuild(&options)
}

fn run_rebuild(options: &RebuildOptions) -> Result<(), OrgraftError> {
    let started = Instant::now();
    let edited_gfa = canonicalize_existing(&options.edited_gfa, "--edited-gfa")?;
    let polished_fasta = canonicalize_existing(&options.polished_fasta, "--polished-fasta")?;
    let image_reference_fasta = options
        .image_reference_fasta
        .as_ref()
        .map(|path| canonicalize_existing(path, "--image-reference-fasta"))
        .transpose()?;
    let merged_template = options
        .merged_gfa_template
        .as_ref()
        .map(|path| canonicalize_existing(path, "--merged-gfa-template"))
        .transpose()?;

    if options.out_dir.exists() {
        if options.force {
            remove_path_if_exists(&options.out_dir)?;
        } else {
            return Err(OrgraftError::InvalidArgument(format!(
                "{} already exists; use --force to replace it",
                options.out_dir.display()
            )));
        }
    }

    let paths = OutputPaths::new(&options.out_dir, &options.subgraph);
    paths.create(options.keep_debug)?;

    let raw_gfa = Gfa::read(&edited_gfa)?;
    let polished_records = read_fasta(&polished_fasta)?;
    if polished_records.len() != 1 {
        return Err(OrgraftError::InvalidArgument(format!(
            "{} must contain exactly one FASTA record, found {}",
            polished_fasta.display(),
            polished_records.len()
        )));
    }
    let polished_record = &polished_records[0];
    let (merged_gfa, merge_mode) = if let Some(template) = &merged_template {
        (Gfa::read(template)?, "template".to_string())
    } else {
        merge_unambiguous_gfa(&raw_gfa)
    };

    if options.keep_debug {
        merged_gfa.write(&paths.debug_dir.join("merged_skeleton.gfa"))?;
    }

    let soft_paths = read_soft_paths(&options.soft_paths)?;
    let minimap2 = options
        .minimap2
        .clone()
        .or_else(|| soft_paths.get("minimap2").cloned())
        .unwrap_or_else(|| PathBuf::from("minimap2"));
    let blastn = options
        .blastn
        .clone()
        .or_else(|| soft_paths.get("blastn").cloned())
        .unwrap_or_else(|| PathBuf::from("blastn"));
    let mapping = run_node_projection(
        &merged_gfa,
        polished_record,
        &minimap2,
        options.threads,
        &paths,
        options.keep_debug,
    )?;

    let node_classes = classify_nodes(&merged_gfa, &mapping);
    let sequence_plan = infer_verified_sequences_from_repeat_cores(
        &merged_gfa,
        polished_record,
        &polished_fasta,
        &mapping,
        &node_classes,
        &blastn,
        &paths,
        options.keep_debug,
    )?;
    let sequence_projection_method = sequence_plan.method.clone();
    let mapping = sequence_plan.mapping;
    let mut verified_gfa = make_verified_gfa_from_sequences(
        &merged_gfa,
        &sequence_plan.sequences,
        &mapping,
        &node_classes,
        &sequence_projection_method,
    );
    let validate_data_dir = discover_validate_data_dir(&polished_fasta);
    let coverage_rows =
        compute_node_remapped_coverage(&mapping, &node_classes, validate_data_dir.as_deref())?;
    annotate_gfa_with_remapped_coverage(&mut verified_gfa, &coverage_rows);
    let repeat_path_rows = repeat_path_support_rows(
        &verified_gfa,
        &mapping,
        &node_classes,
        polished_record,
        validate_data_dir.as_deref(),
    )?;
    attach_repeat_path_support_paths(&mut verified_gfa, &repeat_path_rows);
    write_single_fasta(
        &paths.verified_fasta,
        &polished_record.header,
        &polished_record.sequence,
    )?;
    verified_gfa.write(&paths.verified_gfa)?;
    write_node_fasta(&paths.verified_nodes_fasta, &verified_gfa)?;
    let (gfa_editor_cli, image_exports) =
        if let Some(image_reference_fasta) = image_reference_fasta.as_deref() {
            match resolve_gfa_editor_cli(&soft_paths) {
                Ok(gfa_editor_cli) => {
                    let image_exports = export_gfa_reference_images(
                        &paths,
                        image_reference_fasta,
                        &gfa_editor_cli,
                        &soft_paths,
                    );
                    for row in image_exports.iter().filter(|row| row.status != "written") {
                        eprintln!(
                            "Warning: optional GFA_Editor {} export {} for {}; see {}",
                            row.format,
                            row.status.replace('_', " "),
                            paths.verified_gfa.display(),
                            paths.run_report.display()
                        );
                    }
                    (Some(gfa_editor_cli), image_exports)
                }
                Err(error) => {
                    eprintln!("Warning: optional GFA_Editor image export skipped: {error}");
                    (
                        None,
                        skipped_gfa_reference_images(&paths, &error.to_string()),
                    )
                }
            }
        } else {
            (None, Vec::new())
        };
    let consistency = CoordinateConsistency::new(
        polished_record,
        &raw_gfa,
        &merged_gfa,
        &verified_gfa,
        &mapping,
    );

    let raw_stats = graph_stats(&raw_gfa);
    let merged_stats = graph_stats(&merged_gfa);
    let verified_stats = graph_stats(&verified_gfa);
    let status =
        if consistency.unmapped_segments == 0 && consistency.linear_tiling_status() == "PASS" {
            "PASS"
        } else {
            "WARN"
        };
    let summary = vec![
        ("status", status.to_string()),
        ("organelle", options.organelle.clone()),
        ("subgraph", options.subgraph.clone()),
        ("edited_gfa", edited_gfa.display().to_string()),
        ("polished_fasta", polished_fasta.display().to_string()),
        ("merge_mode", merge_mode.clone()),
        ("repeat_resolution", "disabled".to_string()),
        (
            "sequence_projection_method",
            sequence_projection_method.clone(),
        ),
        ("raw_segments", raw_stats.segments.to_string()),
        ("raw_links", raw_stats.links.to_string()),
        ("raw_bp", raw_stats.bp.to_string()),
        ("merged_segments", merged_stats.segments.to_string()),
        ("merged_links", merged_stats.links.to_string()),
        ("merged_bp", merged_stats.bp.to_string()),
        ("verified_segments", verified_stats.segments.to_string()),
        ("verified_links", verified_stats.links.to_string()),
        ("verified_bp", verified_stats.bp.to_string()),
        ("mapped_segments", consistency.mapped_segments.to_string()),
        (
            "unmapped_segments",
            consistency.unmapped_segments.to_string(),
        ),
        ("polished_length", consistency.polished_length.to_string()),
        ("covered_bases", consistency.covered_bases.to_string()),
        ("gap_bases", consistency.gap_bases.to_string()),
        (
            "multi_covered_bases",
            consistency.multi_covered_bases.to_string(),
        ),
        (
            "linear_tiling",
            consistency.linear_tiling_status().to_string(),
        ),
        (
            "coverage_fraction",
            format!("{:.6}", consistency.coverage_fraction),
        ),
        (
            "runtime_seconds",
            format!("{:.6}", started.elapsed().as_secs_f64()),
        ),
    ];
    write_rebuild_extract(
        &paths.extract_report,
        &options.subgraph,
        &verified_gfa,
        &mapping,
        &node_classes,
        &sequence_plan.node_source_rows,
        &coverage_rows,
    )?;
    write_run_report(
        &paths.run_report,
        &options.out_dir,
        &paths,
        options,
        &edited_gfa,
        &polished_fasta,
        image_reference_fasta.as_deref(),
        validate_data_dir.as_deref(),
        &merge_mode,
        &sequence_projection_method,
        &minimap2,
        &blastn,
        gfa_editor_cli.as_deref(),
        &mapping,
        &image_exports,
        &started,
    )?;
    write_result_stats(
        &paths.result_stats,
        &options.subgraph,
        &summary,
        &raw_stats,
        &merged_stats,
        &verified_stats,
        &consistency,
        &coverage_rows,
        &repeat_path_rows,
    )?;

    println!("Wrote {}", paths.verified_gfa.display());
    println!("Wrote {}", paths.run_report.display());
    Ok(())
}

fn contract() -> CommandContract {
    CommandContract {
        command: "rebuild",
        origin: "high-quality graph generation after orgraft resolve and polish",
        purpose: "rebuild final verified graph and compact reports",
        inputs: &[
            "--organelle NAME",
            "--subgraph ID",
            "--edited-gfa FILE",
            "--polished-fasta FILE",
            "optional --image-reference-fasta FILE",
            "optional --merged-gfa-template FILE",
            "soft_paths.txt containing minimap2, blastn, and optional gfa_editor_cli",
        ],
        outputs: &[
            "subgraph_001/rebuild_subgraph_001.gfa",
            "subgraph_001/rebuild_subgraph_001.fasta",
            "subgraph_001/rebuild_subgraph_001_nodes.fasta",
            "subgraph_001/rebuild_subgraph_001.pdf",
            "subgraph_001/rebuild_subgraph_001.svg",
            "logs/rebuild_subgraph_001_extract.tsv",
            "logs/rebuild_subgraph_001_run_report.tsv",
            "logs/rebuild_subgraph_001_result_stats.tsv",
        ],
        notes: &[
            "rebuild logic is Rust; optional PDF/SVG export calls the GFA_Editor CLI",
            "auto merge follows the conservative non-repeat linear compaction used by resolve/verified_gfa",
            "auto repeat resolution is intentionally disabled",
            "minimap2 projects merged graph nodes; blastn self-BLAST finds maximal exact repeat cores",
        ],
    }
}

#[derive(Debug, Clone)]
struct RebuildOptions {
    edited_gfa: PathBuf,
    polished_fasta: PathBuf,
    out_dir: PathBuf,
    soft_paths: PathBuf,
    image_reference_fasta: Option<PathBuf>,
    merged_gfa_template: Option<PathBuf>,
    minimap2: Option<PathBuf>,
    blastn: Option<PathBuf>,
    organelle: String,
    subgraph: String,
    threads: usize,
    force: bool,
    keep_debug: bool,
}

impl RebuildOptions {
    fn from_args(args: &[String]) -> Result<Self, OrgraftError> {
        let mut edited_gfa = None;
        let mut polished_fasta = None;
        let mut out_dir = PathBuf::from(DEFAULT_OUT_DIR);
        let mut soft_paths = PathBuf::from(DEFAULT_SOFT_PATHS);
        let mut image_reference_fasta = None;
        let mut merged_gfa_template = None;
        let mut minimap2 = None;
        let mut blastn = None;
        let mut organelle = DEFAULT_ORGANELLE.to_string();
        let mut subgraph = DEFAULT_SUBGRAPH.to_string();
        let mut threads = DEFAULT_THREADS;
        let mut force = false;
        let mut keep_debug = false;

        let mut index = 0usize;
        while index < args.len() {
            let arg = &args[index];
            match arg.as_str() {
                "--edited-gfa" | "--checked-draft-gfa" | "--draft-gfa" | "--raw-gfa" | "--gfa" => {
                    edited_gfa = Some(PathBuf::from(required_value(args, &mut index, arg)?));
                }
                "--polished-fasta" | "--verified-fasta" | "--final-fasta" | "--fasta" => {
                    polished_fasta = Some(PathBuf::from(required_value(args, &mut index, arg)?));
                }
                "--out-dir" | "--output-dir" => {
                    out_dir = PathBuf::from(required_value(args, &mut index, arg)?);
                }
                "--soft-paths" => {
                    soft_paths = PathBuf::from(required_value(args, &mut index, arg)?);
                }
                "--image-reference-fasta" | "--reference-fasta" => {
                    image_reference_fasta =
                        Some(PathBuf::from(required_value(args, &mut index, arg)?));
                }
                "--merged-gfa-template" | "--merged-template" => {
                    merged_gfa_template =
                        Some(PathBuf::from(required_value(args, &mut index, arg)?));
                }
                "--blastn" => {
                    blastn = Some(PathBuf::from(required_value(args, &mut index, arg)?));
                }
                "--minimap2" => {
                    minimap2 = Some(PathBuf::from(required_value(args, &mut index, arg)?));
                }
                "--organelle" => {
                    organelle = required_value(args, &mut index, arg)?.to_string();
                }
                "--subgraph" | "--ring" => {
                    subgraph = parse_subgraph_id(required_value(args, &mut index, arg)?)?;
                }
                "--threads" => {
                    threads = parse_usize(required_value(args, &mut index, arg)?, arg)?;
                    if threads == 0 {
                        return Err(OrgraftError::InvalidArgument(
                            "--threads must be at least 1".to_string(),
                        ));
                    }
                }
                "--force" => force = true,
                "--keep-debug" => keep_debug = true,
                other => {
                    return Err(OrgraftError::InvalidArgument(format!(
                        "unknown orgraft rebuild option `{other}`"
                    )));
                }
            }
            index += 1;
        }

        Ok(Self {
            edited_gfa: edited_gfa.ok_or_else(|| {
                OrgraftError::InvalidArgument("missing --edited-gfa FILE".to_string())
            })?,
            polished_fasta: polished_fasta.ok_or_else(|| {
                OrgraftError::InvalidArgument("missing --polished-fasta FILE".to_string())
            })?,
            out_dir,
            soft_paths,
            image_reference_fasta,
            merged_gfa_template,
            minimap2,
            blastn,
            organelle,
            subgraph,
            threads,
            force,
            keep_debug,
        })
    }
}

#[derive(Debug, Clone)]
struct OutputPaths {
    subgraph_dir: PathBuf,
    logs_dir: PathBuf,
    debug_dir: PathBuf,
    verified_gfa: PathBuf,
    verified_fasta: PathBuf,
    verified_nodes_fasta: PathBuf,
    verified_pdf: PathBuf,
    verified_svg: PathBuf,
    extract_report: PathBuf,
    run_report: PathBuf,
    result_stats: PathBuf,
}

impl OutputPaths {
    fn new(root: &Path, subgraph: &str) -> Self {
        let stem = format!("rebuild_{subgraph}");
        let subgraph_dir = root.join(subgraph);
        let logs_dir = root.join("logs");
        Self {
            subgraph_dir: subgraph_dir.clone(),
            logs_dir: logs_dir.clone(),
            debug_dir: root.join("debug"),
            verified_gfa: subgraph_dir.join(format!("{stem}.gfa")),
            verified_fasta: subgraph_dir.join(format!("{stem}.fasta")),
            verified_nodes_fasta: subgraph_dir.join(format!("{stem}_nodes.fasta")),
            verified_pdf: subgraph_dir.join(format!("{stem}.pdf")),
            verified_svg: subgraph_dir.join(format!("{stem}.svg")),
            extract_report: logs_dir.join(format!("{stem}_extract.tsv")),
            run_report: logs_dir.join(format!("{stem}_run_report.tsv")),
            result_stats: logs_dir.join(format!("{stem}_result_stats.tsv")),
        }
    }

    fn create(&self, keep_debug: bool) -> Result<(), OrgraftError> {
        fs::create_dir_all(&self.subgraph_dir)?;
        fs::create_dir_all(&self.logs_dir)?;
        if keep_debug {
            fs::create_dir_all(&self.debug_dir)?;
        }
        Ok(())
    }
}

fn relative_output_path(root: &Path, path: &Path) -> String {
    path.strip_prefix(root)
        .unwrap_or(path)
        .to_string_lossy()
        .replace(std::path::MAIN_SEPARATOR, "/")
}

#[derive(Debug, Clone)]
struct FastaRecord {
    header: String,
    id: String,
    sequence: String,
}

#[derive(Debug, Clone)]
struct Segment {
    name: String,
    sequence: String,
    tags: Vec<String>,
}

#[derive(Debug, Clone)]
struct Link {
    from_name: String,
    from_orient: char,
    to_name: String,
    to_orient: char,
    overlap: String,
    tags: Vec<String>,
}

#[derive(Debug, Clone)]
struct Gfa {
    headers: Vec<Vec<String>>,
    segments: BTreeMap<String, Segment>,
    order: Vec<String>,
    links: Vec<Link>,
    other_lines: Vec<Vec<String>>,
}

impl Gfa {
    fn read(path: &Path) -> Result<Self, OrgraftError> {
        let text = fs::read_to_string(path).map_err(|error| {
            OrgraftError::InvalidArgument(format!("cannot read {}: {error}", path.display()))
        })?;
        let mut headers = Vec::new();
        let mut segments = BTreeMap::new();
        let mut order = Vec::new();
        let mut links = Vec::new();
        let mut other_lines = Vec::new();
        for (line_index, line) in text.lines().enumerate() {
            if line.trim().is_empty() {
                continue;
            }
            let fields: Vec<String> = line.split('\t').map(str::to_string).collect();
            match fields.first().map(String::as_str) {
                Some("H") => headers.push(fields),
                Some("S") => {
                    if fields.len() < 3 {
                        return Err(OrgraftError::InvalidArgument(format!(
                            "{}:{} malformed S line",
                            path.display(),
                            line_index + 1
                        )));
                    }
                    let name = fields[1].clone();
                    order.push(name.clone());
                    segments.insert(
                        name.clone(),
                        Segment {
                            name,
                            sequence: fields[2].clone(),
                            tags: fields[3..].to_vec(),
                        },
                    );
                }
                Some("L") => {
                    if fields.len() < 6 {
                        return Err(OrgraftError::InvalidArgument(format!(
                            "{}:{} malformed L line",
                            path.display(),
                            line_index + 1
                        )));
                    }
                    links.push(Link {
                        from_name: fields[1].clone(),
                        from_orient: parse_orient(&fields[2], path, line_index + 1)?,
                        to_name: fields[3].clone(),
                        to_orient: parse_orient(&fields[4], path, line_index + 1)?,
                        overlap: fields[5].clone(),
                        tags: fields[6..].to_vec(),
                    });
                }
                Some(_) => other_lines.push(fields),
                None => {}
            }
        }
        if segments.is_empty() {
            return Err(OrgraftError::InvalidArgument(format!(
                "{} contains no segments",
                path.display()
            )));
        }
        Ok(Self {
            headers,
            segments,
            order,
            links,
            other_lines,
        })
    }

    fn write(&self, path: &Path) -> Result<(), OrgraftError> {
        let mut out = File::create(path)?;
        for fields in &self.headers {
            writeln!(out, "{}", fields.join("\t"))?;
        }
        for name in &self.order {
            let segment = &self.segments[name];
            let mut fields = vec![
                "S".to_string(),
                segment.name.clone(),
                segment.sequence.clone(),
            ];
            fields.extend(segment.tags.clone());
            writeln!(out, "{}", fields.join("\t"))?;
        }
        for link in &self.links {
            let mut fields = vec![
                "L".to_string(),
                link.from_name.clone(),
                link.from_orient.to_string(),
                link.to_name.clone(),
                link.to_orient.to_string(),
                link.overlap.clone(),
            ];
            fields.extend(link.tags.clone());
            writeln!(out, "{}", fields.join("\t"))?;
        }
        for fields in &self.other_lines {
            writeln!(out, "{}", fields.join("\t"))?;
        }
        Ok(())
    }

    fn degrees(&self) -> BTreeMap<String, usize> {
        let mut degrees = BTreeMap::new();
        for name in &self.order {
            degrees.insert(name.clone(), 0usize);
        }
        for link in &self.links {
            *degrees.entry(link.from_name.clone()).or_insert(0) += 1;
            *degrees.entry(link.to_name.clone()).or_insert(0) += 1;
        }
        degrees
    }
}

#[derive(Debug, Clone)]
struct PathInfo {
    path: Vec<String>,
    orientations: Vec<char>,
    index: usize,
}

fn merge_unambiguous_gfa(raw: &Gfa) -> (Gfa, String) {
    let protected = repeat_like_nodes(raw);
    let mut adjacency: BTreeMap<String, Vec<String>> = BTreeMap::new();
    for name in &raw.order {
        adjacency.insert(name.clone(), Vec::new());
    }
    for link in &raw.links {
        if protected.contains(&link.from_name) || protected.contains(&link.to_name) {
            continue;
        }
        adjacency
            .entry(link.from_name.clone())
            .or_default()
            .push(link.to_name.clone());
        adjacency
            .entry(link.to_name.clone())
            .or_default()
            .push(link.from_name.clone());
    }

    let mut visited = HashSet::new();
    let mut components = Vec::new();
    for name in &raw.order {
        if visited.contains(name) || protected.contains(name) {
            continue;
        }
        let mut stack = vec![name.clone()];
        let mut component = Vec::new();
        visited.insert(name.clone());
        while let Some(current) = stack.pop() {
            component.push(current.clone());
            for next in adjacency.get(&current).cloned().unwrap_or_default() {
                if visited.insert(next.clone()) {
                    stack.push(next);
                }
            }
        }
        components.push(component);
    }

    let mut path_records: Vec<(Vec<String>, Vec<char>)> = Vec::new();
    let mut used = HashSet::new();
    for component in components {
        if let Some(order) = ordered_linear_component(&component, &adjacency) {
            let orientations = orient_path(raw, &order);
            used.extend(order.iter().cloned());
            path_records.push((order, orientations));
        } else {
            let mut names = component;
            names.sort_by(|a, b| natural_cmp(a, b));
            for name in names {
                used.insert(name.clone());
                path_records.push((vec![name], vec!['+']));
            }
        }
    }
    for name in &raw.order {
        if used.insert(name.clone()) {
            path_records.push((vec![name.clone()], vec!['+']));
        }
    }

    let mut old_to_new = HashMap::new();
    let mut old_to_path_info = HashMap::new();
    let mut segments = BTreeMap::new();
    let mut order = Vec::new();
    let mut mergeable_path_count = 0usize;
    for (path, orientations) in &path_records {
        if path.len() > 1 {
            mergeable_path_count += 1;
        }
        let name = path.join("_");
        let mut sequence = String::new();
        for (index, (old_name, orient)) in path.iter().zip(orientations.iter()).enumerate() {
            let mut part = raw.segments[old_name].sequence.clone();
            if *orient == '-' {
                part = reverse_complement(&part);
            }
            sequence.push_str(&part);
            old_to_new.insert(old_name.clone(), name.clone());
            old_to_path_info.insert(
                old_name.clone(),
                PathInfo {
                    path: path.clone(),
                    orientations: orientations.clone(),
                    index,
                },
            );
        }
        let tags = vec![
            format!("LN:i:{}", sequence.len()),
            if path.len() > 1 {
                "SC:Z:linear_compaction"
            } else {
                "SC:Z:preserved_node"
            }
            .to_string(),
            "RR:Z:disabled".to_string(),
        ];
        order.push(name.clone());
        segments.insert(
            name.clone(),
            Segment {
                name,
                sequence,
                tags,
            },
        );
    }

    let mut links = Vec::new();
    let mut seen = HashSet::new();
    for link in &raw.links {
        if let Some(converted) = convert_link(link, &old_to_new, &old_to_path_info) {
            let key = format!(
                "{}{}:{}{}:{}",
                converted.from_name,
                converted.from_orient,
                converted.to_name,
                converted.to_orient,
                converted.overlap
            );
            if seen.insert(key) {
                links.push(converted);
            }
        }
    }

    let mode = if mergeable_path_count == 0 {
        "input_already_merged_or_no_linear_compaction"
    } else {
        "auto_non_repeat_linear_compaction"
    };
    (
        Gfa {
            headers: raw.headers.clone(),
            segments,
            order,
            links,
            other_lines: raw.other_lines.clone(),
        },
        mode.to_string(),
    )
}

fn repeat_like_nodes(gfa: &Gfa) -> BTreeSet<String> {
    let degrees = gfa.degrees();
    let branch: Vec<usize> = gfa
        .order
        .iter()
        .filter_map(|name| degrees.get(name).copied())
        .filter(|degree| *degree >= 3)
        .collect();
    if branch.is_empty() {
        return BTreeSet::new();
    }
    let branch_floor = *branch.iter().min().unwrap();
    let has_higher = branch.iter().any(|degree| *degree > branch_floor);
    gfa.order
        .iter()
        .filter(|name| {
            let degree = degrees.get(*name).copied().unwrap_or(0);
            if has_higher {
                degree > branch_floor
            } else {
                degree >= 3
            }
        })
        .cloned()
        .collect()
}

fn ordered_linear_component(
    component: &[String],
    adjacency: &BTreeMap<String, Vec<String>>,
) -> Option<Vec<String>> {
    let component_set: HashSet<String> = component.iter().cloned().collect();
    let mut sub_degrees = BTreeMap::new();
    for name in component {
        let degree = adjacency
            .get(name)
            .unwrap_or(&Vec::new())
            .iter()
            .filter(|next| component_set.contains(*next))
            .count();
        if degree > 2 {
            return None;
        }
        sub_degrees.insert(name.clone(), degree);
    }
    let mut endpoints: Vec<String> = sub_degrees
        .iter()
        .filter(|(_, degree)| **degree <= 1)
        .map(|(name, _)| name.clone())
        .collect();
    if component.len() > 1 && endpoints.len() != 2 {
        return None;
    }
    endpoints.sort_by(|a, b| natural_cmp(a, b));
    let mut order = Vec::new();
    let mut previous: Option<String> = None;
    let mut current = endpoints
        .first()
        .cloned()
        .or_else(|| component.first().cloned())?;
    loop {
        order.push(current.clone());
        let mut candidates: Vec<String> = adjacency
            .get(&current)
            .unwrap_or(&Vec::new())
            .iter()
            .filter(|next| component_set.contains(*next))
            .filter(|next| previous.as_ref() != Some(*next))
            .cloned()
            .collect();
        candidates.sort_by(|a, b| natural_cmp(a, b));
        let Some(next) = candidates.first().cloned() else {
            break;
        };
        previous = Some(current);
        current = next;
    }
    if order.len() == component.len() {
        Some(order)
    } else {
        None
    }
}

fn orient_path(gfa: &Gfa, order: &[String]) -> Vec<char> {
    for start_orient in ['+', '-'] {
        let mut orientations = vec![start_orient];
        let mut ok = true;
        for pair in order.windows(2) {
            let current = &pair[0];
            let next = &pair[1];
            let Some(orient) = next_orient(gfa, current, *orientations.last().unwrap(), next)
            else {
                ok = false;
                break;
            };
            orientations.push(orient);
        }
        if ok {
            return orientations;
        }
    }
    vec!['+'; order.len()]
}

fn next_orient(gfa: &Gfa, current: &str, current_orient: char, next: &str) -> Option<char> {
    for link in &gfa.links {
        if link.from_name == current && link.to_name == next && link.from_orient == current_orient {
            return Some(link.to_orient);
        }
        if link.to_name == current
            && link.from_name == next
            && flip_orient(link.to_orient) == current_orient
        {
            return Some(flip_orient(link.from_orient));
        }
    }
    None
}

fn convert_link(
    link: &Link,
    old_to_new: &HashMap<String, String>,
    old_to_path_info: &HashMap<String, PathInfo>,
) -> Option<Link> {
    let merged_from = old_to_new.get(&link.from_name)?;
    let merged_to = old_to_new.get(&link.to_name)?;
    if merged_from == merged_to {
        return None;
    }
    let from_endpoint = old_from_endpoint(link.from_orient);
    let to_endpoint = old_to_endpoint(link.to_orient);
    let merged_from_endpoint =
        map_old_endpoint_to_merged(&link.from_name, from_endpoint, old_to_path_info)?;
    let merged_to_endpoint =
        map_old_endpoint_to_merged(&link.to_name, to_endpoint, old_to_path_info)?;
    Some(Link {
        from_name: merged_from.clone(),
        from_orient: orient_for_from_endpoint(merged_from_endpoint),
        to_name: merged_to.clone(),
        to_orient: orient_for_to_endpoint(merged_to_endpoint),
        overlap: link.overlap.clone(),
        tags: link.tags.clone(),
    })
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum SegmentEndpoint {
    Left,
    Right,
}

fn old_from_endpoint(orient: char) -> SegmentEndpoint {
    match orient {
        '+' => SegmentEndpoint::Right,
        '-' => SegmentEndpoint::Left,
        _ => SegmentEndpoint::Right,
    }
}

fn old_to_endpoint(orient: char) -> SegmentEndpoint {
    match orient {
        '+' => SegmentEndpoint::Left,
        '-' => SegmentEndpoint::Right,
        _ => SegmentEndpoint::Left,
    }
}

fn orient_for_from_endpoint(endpoint: SegmentEndpoint) -> char {
    match endpoint {
        SegmentEndpoint::Right => '+',
        SegmentEndpoint::Left => '-',
    }
}

fn orient_for_to_endpoint(endpoint: SegmentEndpoint) -> char {
    match endpoint {
        SegmentEndpoint::Left => '+',
        SegmentEndpoint::Right => '-',
    }
}

fn map_old_endpoint_to_merged(
    old_name: &str,
    endpoint: SegmentEndpoint,
    old_to_path_info: &HashMap<String, PathInfo>,
) -> Option<SegmentEndpoint> {
    let info = old_to_path_info.get(old_name)?;
    let path_orient = info.orientations[info.index];
    if info.path.len() == 1 {
        return Some(match (path_orient, endpoint) {
            ('+', SegmentEndpoint::Left) | ('-', SegmentEndpoint::Right) => SegmentEndpoint::Left,
            ('+', SegmentEndpoint::Right) | ('-', SegmentEndpoint::Left) => SegmentEndpoint::Right,
            _ => endpoint,
        });
    }
    let first_endpoint = match path_orient {
        '+' => SegmentEndpoint::Left,
        '-' => SegmentEndpoint::Right,
        _ => SegmentEndpoint::Left,
    };
    if info.index == 0 && endpoint == first_endpoint {
        return Some(SegmentEndpoint::Left);
    }
    let last_endpoint = match path_orient {
        '+' => SegmentEndpoint::Right,
        '-' => SegmentEndpoint::Left,
        _ => SegmentEndpoint::Right,
    };
    if info.index == info.path.len() - 1 && endpoint == last_endpoint {
        return Some(SegmentEndpoint::Right);
    }
    None
}

#[derive(Debug, Clone)]
struct PafHit {
    query_name: String,
    query_length: usize,
    query_start: usize,
    query_end: usize,
    strand: char,
    target_start: usize,
    target_end: usize,
    matches: usize,
    block_length: usize,
    mapq: usize,
    primary: bool,
}

impl PafHit {
    fn identity(&self) -> f64 {
        if self.block_length == 0 {
            0.0
        } else {
            self.matches as f64 / self.block_length as f64
        }
    }
    fn aligned_fraction(&self) -> f64 {
        if self.query_length == 0 {
            0.0
        } else {
            (self.query_end.saturating_sub(self.query_start)) as f64 / self.query_length as f64
        }
    }
    fn folded(&self, reference_len: usize) -> FoldedInterval {
        FoldedInterval::new(self.target_start, self.target_end, reference_len)
    }
}

#[derive(Debug, Clone)]
struct FoldedInterval {
    start: usize,
    end: usize,
    wraps: bool,
}

impl FoldedInterval {
    fn new(start_zero: usize, end_zero: usize, reference_len: usize) -> Self {
        let span = end_zero.saturating_sub(start_zero);
        let start = (start_zero % reference_len) + 1;
        let end = ((end_zero.saturating_sub(1)) % reference_len) + 1;
        Self {
            start,
            end,
            wraps: span > 0 && start > end,
        }
    }
    fn parts_text(&self, reference_len: usize) -> String {
        if self.wraps {
            format!("{}-{},1-{}", self.start, reference_len, self.end)
        } else {
            format!("{}-{}", self.start, self.end)
        }
    }
}

#[derive(Debug, Clone)]
struct NodeMapping {
    selected: Option<PafHit>,
    accepted_hits: Vec<PafHit>,
    all_hits: Vec<PafHit>,
}

#[derive(Debug, Clone)]
struct Mapping {
    by_node: BTreeMap<String, NodeMapping>,
    reference_len: usize,
    command: String,
    stderr: String,
}

fn run_node_projection(
    gfa: &Gfa,
    record: &FastaRecord,
    minimap2: &Path,
    threads: usize,
    paths: &OutputPaths,
    keep_debug: bool,
) -> Result<Mapping, OrgraftError> {
    let temp_dir = if keep_debug {
        paths.debug_dir.clone()
    } else {
        temp_work_dir("orgraft-rebuild")
    };
    fs::create_dir_all(&temp_dir)?;
    let query_fasta = temp_dir.join("merged_nodes.fasta");
    let target_fasta = temp_dir.join("polished_doubled.fasta");
    write_node_fasta_from_order(&query_fasta, gfa)?;
    write_single_fasta(
        &target_fasta,
        &format!("{}_doubled", record.id),
        &format!("{}{}", record.sequence, record.sequence),
    )?;

    let command_text = format!(
        "{} -x asm5 -c --eqx -N 20 -t {} {} {}",
        minimap2.display(),
        threads,
        target_fasta.display(),
        query_fasta.display()
    );
    let output = Command::new(minimap2)
        .args([
            "-x",
            "asm5",
            "-c",
            "--eqx",
            "-N",
            "20",
            "-t",
            &threads.to_string(),
        ])
        .arg(&target_fasta)
        .arg(&query_fasta)
        .output()
        .map_err(|error| {
            OrgraftError::InvalidArgument(format!("failed to run minimap2: {error}"))
        })?;
    let stdout = String::from_utf8_lossy(&output.stdout).to_string();
    let stderr = String::from_utf8_lossy(&output.stderr).to_string();
    if keep_debug {
        fs::write(paths.debug_dir.join("minimap2.paf"), &stdout)?;
        fs::write(paths.debug_dir.join("minimap2.stderr.log"), &stderr)?;
    }
    if !output.status.success() {
        if !keep_debug {
            let _ = fs::remove_dir_all(&temp_dir);
        }
        return Err(OrgraftError::InvalidArgument(format!(
            "minimap2 failed: {stderr}"
        )));
    }

    let mut raw_hits: BTreeMap<String, Vec<PafHit>> = BTreeMap::new();
    for line in stdout.lines() {
        if line.trim().is_empty() {
            continue;
        }
        let hit = parse_paf(line)?;
        raw_hits
            .entry(hit.query_name.clone())
            .or_default()
            .push(hit);
    }
    let mut by_node = BTreeMap::new();
    for name in &gfa.order {
        let mut hits = raw_hits.remove(name).unwrap_or_default();
        hits.retain(|hit| hit.identity() >= 0.97 && hit.aligned_fraction() >= 0.05);
        sort_hits(&mut hits);
        let accepted_hits = dedupe_hits(
            hits.iter()
                .filter(|hit| hit.identity() >= 0.99 && hit.aligned_fraction() >= 0.75)
                .cloned()
                .collect(),
            record.sequence.len(),
        );
        by_node.insert(
            name.clone(),
            NodeMapping {
                selected: hits.first().cloned(),
                accepted_hits,
                all_hits: hits,
            },
        );
    }
    if !keep_debug {
        let _ = fs::remove_dir_all(&temp_dir);
    }
    Ok(Mapping {
        by_node,
        reference_len: record.sequence.len(),
        command: command_text,
        stderr,
    })
}

fn classify_nodes(gfa: &Gfa, mapping: &Mapping) -> BTreeMap<String, String> {
    let repeats = repeat_like_nodes(gfa);
    let mut classes = BTreeMap::new();
    for name in &gfa.order {
        let mapped_repeat = mapping
            .by_node
            .get(name)
            .map(|node| node.accepted_hits.len() > 1)
            .unwrap_or(false);
        let class = if repeats.contains(name) || mapped_repeat {
            "repeat_node"
        } else {
            "single_copy_node"
        };
        classes.insert(name.clone(), class.to_string());
    }
    classes
}

#[derive(Debug, Clone)]
struct VerifiedNodeSourceRow {
    node: String,
    node_class: String,
    source: String,
    length: usize,
    parts: String,
    notes: String,
}

#[derive(Debug, Clone)]
struct SequencePlan {
    sequences: BTreeMap<String, String>,
    mapping: Mapping,
    node_source_rows: Vec<VerifiedNodeSourceRow>,
    method: String,
}

fn make_verified_gfa_from_sequences(
    merged: &Gfa,
    sequences: &BTreeMap<String, String>,
    mapping: &Mapping,
    node_classes: &BTreeMap<String, String>,
    method: &str,
) -> Gfa {
    let mut segments = BTreeMap::new();
    for name in &merged.order {
        let segment = &merged.segments[name];
        let node_mapping = &mapping.by_node[name];
        let sequence = sequences
            .get(name)
            .cloned()
            .unwrap_or_else(|| segment.sequence.clone());
        let interval_text = if node_mapping.accepted_hits.is_empty() {
            ".".to_string()
        } else {
            node_mapping
                .accepted_hits
                .iter()
                .map(|hit| {
                    let interval = hit.folded(mapping.reference_len);
                    format!(
                        "{}:{}",
                        interval.parts_text(mapping.reference_len),
                        hit.strand
                    )
                })
                .collect::<Vec<_>>()
                .join(",")
        };
        let mut tags = replace_tags(
            &segment.tags,
            &[
                format!("VC:Z:{method}"),
                "VZ:Z:verified_fasta".to_string(),
                format!("VS:Z:{interval_text}"),
                format!("NC:Z:{}", node_classes[name]),
                format!("LN:i:{}", sequence.len()),
                format!("OL:i:{}", segment.sequence.len()),
                format!("VL:i:{}", sequence.len()),
            ],
        );
        if node_classes[name] == "repeat_node" {
            tags = replace_tag(&tags, "RT:Z:unresolved".to_string());
            tags = replace_tag(&tags, "PM:Z:edge_rc_min_proxy".to_string());
        }
        segments.insert(
            name.clone(),
            Segment {
                name: name.clone(),
                sequence,
                tags,
            },
        );
    }
    Gfa {
        headers: merged.headers.clone(),
        segments,
        order: merged.order.clone(),
        links: merged.links.clone(),
        other_lines: merged.other_lines.clone(),
    }
}

#[derive(Debug, Clone)]
struct BlastRow {
    pident: f64,
    length: usize,
    mismatch: usize,
    gapopen: usize,
    qstart: usize,
    qend: usize,
    sstart: usize,
    send: usize,
    bitscore: f64,
}

#[derive(Debug, Clone)]
struct RepeatCopy {
    start: usize,
    end: usize,
    strand: char,
}

#[derive(Debug, Clone)]
struct RepeatCorePair {
    length: usize,
    source_blast_length: usize,
    source_pident: f64,
    source_gapopen: usize,
    source_mismatch: usize,
    boundary_extension_left: usize,
    boundary_extension_right: usize,
    copies: Vec<RepeatCopy>,
}

#[derive(Debug, Clone)]
struct RepeatGap {
    left_repeat_node: String,
    left_copy_index: usize,
    right_repeat_node: String,
    right_copy_index: usize,
    parts: Vec<(usize, usize)>,
}

fn infer_verified_sequences_from_repeat_cores(
    merged: &Gfa,
    record: &FastaRecord,
    verified_fasta: &Path,
    mapping: &Mapping,
    node_classes: &BTreeMap<String, String>,
    blastn: &Path,
    paths: &OutputPaths,
    keep_debug: bool,
) -> Result<SequencePlan, OrgraftError> {
    let repeat_count = node_classes
        .values()
        .filter(|class| class.as_str() == "repeat_node")
        .count();
    if repeat_count == 0 {
        return Ok(fallback_sequence_plan(
            merged,
            record,
            mapping,
            node_classes,
            "minimap2_fallback",
        ));
    }

    let core_pairs = infer_repeat_core_pairs(verified_fasta, record, blastn, repeat_count)?;
    if keep_debug {
        write_self_blast_repeat_core_debug(
            &paths.debug_dir.join("self_blast_repeat_cores.tsv"),
            &core_pairs,
        )?;
    }
    if core_pairs.len() < repeat_count {
        return Ok(fallback_sequence_plan(
            merged,
            record,
            mapping,
            node_classes,
            "minimap2_fallback",
        ));
    }
    let repeat_assignments =
        assign_repeat_nodes_to_core_pairs(merged, node_classes, mapping, &core_pairs);
    if repeat_assignments.len() < repeat_count {
        return Ok(fallback_sequence_plan(
            merged,
            record,
            mapping,
            node_classes,
            "minimap2_fallback",
        ));
    }

    let gaps = build_repeat_core_gaps(&repeat_assignments, record.sequence.len());
    let single_assignments = assign_single_copy_nodes_to_gaps(merged, node_classes, mapping, &gaps);
    let mut sequences = BTreeMap::new();
    let mut by_node = BTreeMap::new();
    let mut rows = Vec::new();

    for name in &merged.order {
        if node_classes[name] == "repeat_node" {
            let Some(pair) = repeat_assignments.get(name) else {
                continue;
            };
            let (copy, strand) = choose_repeat_copy_for_node(name, pair, mapping);
            let sequence = record_interval_sequence(record, copy.start, copy.end, strand);
            let synthetic_hits =
                synthetic_hits_for_repeat_copies(name, &sequence, record, &pair.copies);
            sequences.insert(name.clone(), sequence.clone());
            by_node.insert(
                name.clone(),
                NodeMapping {
                    selected: synthetic_hits.first().cloned(),
                    accepted_hits: synthetic_hits.clone(),
                    all_hits: synthetic_hits,
                },
            );
            rows.push(VerifiedNodeSourceRow {
                node: name.clone(),
                node_class: "repeat_node".to_string(),
                source: "verified_fasta_self_blast_repeat_core".to_string(),
                length: sequence.len(),
                parts: format!("{}-{}:{}", copy.start, copy.end, strand),
                notes: format!(
                    "core_copies={};blast_length={};blast_pident={:.3};boundary_extension_left={};boundary_extension_right={}",
                    pair.copies
                        .iter()
                        .map(|item| format!("{}-{}:{}", item.start, item.end, item.strand))
                        .collect::<Vec<_>>()
                        .join(","),
                    pair.source_blast_length,
                    pair.source_pident,
                    pair.boundary_extension_left,
                    pair.boundary_extension_right
                ),
            });
        }
    }

    for name in &merged.order {
        if node_classes[name] != "single_copy_node" {
            continue;
        }
        let Some(gap) = single_assignments.get(name) else {
            continue;
        };
        let strand = best_orientation_for_parts(name, &gap.parts, mapping);
        let sequence = verified_sequence_from_parts(record, &gap.parts, strand);
        let synthetic_hits =
            synthetic_hits_for_parts(name, sequence.len(), record, &gap.parts, strand);
        sequences.insert(name.clone(), sequence.clone());
        by_node.insert(
            name.clone(),
            NodeMapping {
                selected: synthetic_hits.first().cloned(),
                accepted_hits: synthetic_hits.clone(),
                all_hits: synthetic_hits,
            },
        );
        rows.push(VerifiedNodeSourceRow {
            node: name.clone(),
            node_class: "single_copy_node".to_string(),
            source: "verified_fasta_between_repeat_cores".to_string(),
            length: sequence.len(),
            parts: format!(
                "{}:{}",
                gap.parts
                    .iter()
                    .map(|(start, end)| format!("{start}-{end}"))
                    .collect::<Vec<_>>()
                    .join(","),
                strand
            ),
            notes: format!(
                "left={}:copy{};right={}:copy{}",
                gap.left_repeat_node,
                gap.left_copy_index,
                gap.right_repeat_node,
                gap.right_copy_index
            ),
        });
    }

    for name in &merged.order {
        if sequences.contains_key(name) {
            continue;
        }
        let segment = &merged.segments[name];
        let original = mapping.by_node.get(name);
        let hit = original.and_then(|node| node.selected.clone());
        let (sequence, accepted_hits, source, parts, notes) = if let Some(hit) = hit {
            let sequence = extract_verified_sequence(record, &hit);
            (
                sequence,
                vec![hit.clone()],
                "verified_fasta_minimap2_fallback".to_string(),
                format!("{}-{}:{}", hit.target_start + 1, hit.target_end, hit.strand),
                "fallback".to_string(),
            )
        } else {
            (
                segment.sequence.clone(),
                Vec::new(),
                "unmapped_original_sequence".to_string(),
                ".".to_string(),
                "fallback_no_alignment".to_string(),
            )
        };
        sequences.insert(name.clone(), sequence.clone());
        by_node.insert(
            name.clone(),
            NodeMapping {
                selected: accepted_hits.first().cloned(),
                accepted_hits: accepted_hits.clone(),
                all_hits: accepted_hits,
            },
        );
        rows.push(VerifiedNodeSourceRow {
            node: name.clone(),
            node_class: node_classes[name].clone(),
            source,
            length: sequence.len(),
            parts,
            notes,
        });
    }

    Ok(SequencePlan {
        sequences,
        mapping: Mapping {
            by_node,
            reference_len: record.sequence.len(),
            command: mapping.command.clone(),
            stderr: mapping.stderr.clone(),
        },
        node_source_rows: rows,
        method: "verified_fasta_self_blast_repeat_core".to_string(),
    })
}

fn fallback_sequence_plan(
    merged: &Gfa,
    record: &FastaRecord,
    mapping: &Mapping,
    node_classes: &BTreeMap<String, String>,
    method: &str,
) -> SequencePlan {
    let mut sequences = BTreeMap::new();
    let mut rows = Vec::new();
    for name in &merged.order {
        let segment = &merged.segments[name];
        let node = &mapping.by_node[name];
        if let Some(hit) = &node.selected {
            let sequence = extract_verified_sequence(record, hit);
            sequences.insert(name.clone(), sequence.clone());
            rows.push(VerifiedNodeSourceRow {
                node: name.clone(),
                node_class: node_classes[name].clone(),
                source: "polished_fasta_minimap2".to_string(),
                length: sequence.len(),
                parts: format!("{}-{}:{}", hit.target_start + 1, hit.target_end, hit.strand),
                notes: format!(
                    "identity={:.6};aligned_fraction={:.6}",
                    hit.identity(),
                    hit.aligned_fraction()
                ),
            });
        } else {
            sequences.insert(name.clone(), segment.sequence.clone());
            rows.push(VerifiedNodeSourceRow {
                node: name.clone(),
                node_class: node_classes[name].clone(),
                source: "unmapped_original_sequence".to_string(),
                length: segment.sequence.len(),
                parts: ".".to_string(),
                notes: "selected_hit=missing".to_string(),
            });
        }
    }
    SequencePlan {
        sequences,
        mapping: mapping.clone(),
        node_source_rows: rows,
        method: method.to_string(),
    }
}

fn infer_repeat_core_pairs(
    verified_fasta: &Path,
    record: &FastaRecord,
    blastn: &Path,
    repeat_count: usize,
) -> Result<Vec<RepeatCorePair>, OrgraftError> {
    let mut rows = run_self_blastn_repeat_hits(verified_fasta, blastn)?;
    rows.retain(|row| row.length >= 1000 && row.pident >= 95.0);
    rows.sort_by(|a, b| {
        b.length
            .cmp(&a.length)
            .then_with(|| b.pident.partial_cmp(&a.pident).unwrap_or(Ordering::Equal))
            .then_with(|| {
                b.bitscore
                    .partial_cmp(&a.bitscore)
                    .unwrap_or(Ordering::Equal)
            })
    });

    let mut selected = Vec::new();
    let mut seen = HashSet::new();
    let mut used_intervals: Vec<(usize, usize)> = Vec::new();
    for row in rows {
        let key = core_pair_key(&row);
        if !seen.insert(key) {
            continue;
        }
        let (q_low, q_high) = sorted_pair(row.qstart, row.qend);
        let (s_low, s_high) = sorted_pair(row.sstart, row.send);
        if used_intervals.iter().any(|(old_start, old_end)| {
            interval_overlap(q_low, q_high, *old_start, *old_end) > 100
                || interval_overlap(s_low, s_high, *old_start, *old_end) > 100
        }) {
            continue;
        }
        let pair = repeat_core_pair_from_blast_row(&row, record);
        if pair.length == 0 {
            continue;
        }
        used_intervals.push((pair.copies[0].start, pair.copies[0].end));
        used_intervals.push((pair.copies[1].start, pair.copies[1].end));
        selected.push(pair);
        if selected.len() >= repeat_count {
            break;
        }
    }
    Ok(selected)
}

fn run_self_blastn_repeat_hits(
    verified_fasta: &Path,
    blastn: &Path,
) -> Result<Vec<BlastRow>, OrgraftError> {
    let output = Command::new(blastn)
        .arg("-query")
        .arg(verified_fasta)
        .arg("-subject")
        .arg(verified_fasta)
        .arg("-outfmt")
        .arg("6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore")
        .arg("-dust")
        .arg("no")
        .arg("-soft_masking")
        .arg("false")
        .output()
        .map_err(|error| {
            OrgraftError::InvalidArgument(format!("failed to run blastn self alignment: {error}"))
        })?;
    let stdout = String::from_utf8_lossy(&output.stdout);
    let stderr = String::from_utf8_lossy(&output.stderr);
    if !output.status.success() {
        return Err(OrgraftError::InvalidArgument(format!(
            "blastn self alignment failed: {stderr}"
        )));
    }
    let mut rows = Vec::new();
    for line in stdout.lines() {
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 12 {
            continue;
        }
        let row = BlastRow {
            pident: fields[2].parse::<f64>().unwrap_or(0.0),
            length: fields[3].parse::<usize>().unwrap_or(0),
            mismatch: fields[4].parse::<usize>().unwrap_or(0),
            gapopen: fields[5].parse::<usize>().unwrap_or(0),
            qstart: fields[6].parse::<usize>().unwrap_or(0),
            qend: fields[7].parse::<usize>().unwrap_or(0),
            sstart: fields[8].parse::<usize>().unwrap_or(0),
            send: fields[9].parse::<usize>().unwrap_or(0),
            bitscore: fields[11].parse::<f64>().unwrap_or(0.0),
        };
        let (q_low, q_high) = sorted_pair(row.qstart, row.qend);
        let (s_low, s_high) = sorted_pair(row.sstart, row.send);
        if q_low == s_low && q_high == s_high {
            continue;
        }
        rows.push(row);
    }
    Ok(rows)
}

fn core_pair_key(row: &BlastRow) -> String {
    let (q_low, q_high) = sorted_pair(row.qstart, row.qend);
    let (s_low, s_high) = sorted_pair(row.sstart, row.send);
    let mut pairs = [(q_low, q_high), (s_low, s_high)];
    pairs.sort();
    format!(
        "{}-{}:{}-{}",
        pairs[0].0, pairs[0].1, pairs[1].0, pairs[1].1
    )
}

fn repeat_core_pair_from_blast_row(row: &BlastRow, record: &FastaRecord) -> RepeatCorePair {
    let (q_seq, q_strand, q_low, q_high) = blast_interval_sequence(record, row.qstart, row.qend);
    let (s_seq, s_strand, s_low, s_high) = blast_interval_sequence(record, row.sstart, row.send);
    let (q_offset, s_offset, core_length) = longest_exact_common_substring(&q_seq, &s_seq);
    let (q_start, q_end, q_core_strand) =
        map_oriented_offset_to_record(q_low, q_high, q_strand, q_offset, core_length);
    let (s_start, s_end, s_core_strand) =
        map_oriented_offset_to_record(s_low, s_high, s_strand, s_offset, core_length);
    let mut pair = RepeatCorePair {
        length: core_length,
        source_blast_length: row.length,
        source_pident: row.pident,
        source_gapopen: row.gapopen,
        source_mismatch: row.mismatch,
        boundary_extension_left: 0,
        boundary_extension_right: 0,
        copies: vec![
            RepeatCopy {
                start: q_start,
                end: q_end,
                strand: q_core_strand,
            },
            RepeatCopy {
                start: s_start,
                end: s_end,
                strand: s_core_strand,
            },
        ],
    };
    extend_repeat_core_pair_to_max_exact(&mut pair, record);
    pair
}

fn blast_interval_sequence(
    record: &FastaRecord,
    start: usize,
    end: usize,
) -> (String, char, usize, usize) {
    if start <= end {
        (
            record_interval_sequence(record, start, end, '+'),
            '+',
            start,
            end,
        )
    } else {
        (
            record_interval_sequence(record, end, start, '-'),
            '-',
            end,
            start,
        )
    }
}

fn longest_exact_common_substring(seq_a: &str, seq_b: &str) -> (usize, usize, usize) {
    let a = seq_a.as_bytes();
    let b = seq_b.as_bytes();
    let mut previous = vec![0usize; b.len() + 1];
    let mut best_start_a = 0usize;
    let mut best_start_b = 0usize;
    let mut best_len = 0usize;
    for (index_a, base_a) in a.iter().enumerate() {
        let mut current = vec![0usize; b.len() + 1];
        for (index_b, base_b) in b.iter().enumerate() {
            if base_a.eq_ignore_ascii_case(base_b) {
                current[index_b + 1] = previous[index_b] + 1;
                if current[index_b + 1] > best_len {
                    best_len = current[index_b + 1];
                    best_start_a = index_a + 1 - best_len;
                    best_start_b = index_b + 1 - best_len;
                }
            }
        }
        previous = current;
    }
    (best_start_a, best_start_b, best_len)
}

fn map_oriented_offset_to_record(
    start: usize,
    end: usize,
    strand: char,
    offset: usize,
    length: usize,
) -> (usize, usize, char) {
    if strand == '+' {
        let core_start = start + offset;
        (core_start, core_start + length.saturating_sub(1), '+')
    } else {
        let core_end = end.saturating_sub(offset);
        (core_end + 1 - length, core_end, '-')
    }
}

fn extend_repeat_core_pair_to_max_exact(pair: &mut RepeatCorePair, record: &FastaRecord) {
    while try_extend_repeat_core_pair(pair, record, "left") {
        pair.boundary_extension_left += 1;
    }
    while try_extend_repeat_core_pair(pair, record, "right") {
        pair.boundary_extension_right += 1;
    }
}

fn try_extend_repeat_core_pair(
    pair: &mut RepeatCorePair,
    record: &FastaRecord,
    side: &str,
) -> bool {
    let reference_len = record.sequence.len();
    let Some(pos_a) = repeat_copy_extension_position(&pair.copies[0], side, reference_len) else {
        return false;
    };
    let Some(pos_b) = repeat_copy_extension_position(&pair.copies[1], side, reference_len) else {
        return false;
    };
    let base_a = oriented_base_at(record, pos_a, pair.copies[0].strand);
    let base_b = oriented_base_at(record, pos_b, pair.copies[1].strand);
    if !base_a.eq_ignore_ascii_case(&base_b) {
        return false;
    }
    let copy_a = extended_repeat_copy(&pair.copies[0], pos_a, side);
    let copy_b = extended_repeat_copy(&pair.copies[1], pos_b, side);
    if interval_overlap(copy_a.start, copy_a.end, copy_b.start, copy_b.end) > 0 {
        return false;
    }
    pair.copies = vec![copy_a, copy_b];
    pair.length += 1;
    true
}

fn repeat_copy_extension_position(
    copy: &RepeatCopy,
    side: &str,
    reference_len: usize,
) -> Option<usize> {
    if side == "left" {
        if copy.strand == '+' {
            (copy.start > 1).then_some(copy.start - 1)
        } else {
            (copy.end < reference_len).then_some(copy.end + 1)
        }
    } else if copy.strand == '+' {
        (copy.end < reference_len).then_some(copy.end + 1)
    } else {
        (copy.start > 1).then_some(copy.start - 1)
    }
}

fn extended_repeat_copy(copy: &RepeatCopy, position: usize, side: &str) -> RepeatCopy {
    let mut updated = copy.clone();
    if side == "left" {
        if updated.strand == '+' {
            updated.start = position;
        } else {
            updated.end = position;
        }
    } else if updated.strand == '+' {
        updated.end = position;
    } else {
        updated.start = position;
    }
    updated
}

fn oriented_base_at(record: &FastaRecord, position: usize, strand: char) -> char {
    let base = record
        .sequence
        .as_bytes()
        .get(position.saturating_sub(1))
        .copied()
        .unwrap_or(b'N') as char;
    if strand == '-' {
        complement_base(base)
    } else {
        base
    }
}

fn assign_repeat_nodes_to_core_pairs(
    merged: &Gfa,
    node_classes: &BTreeMap<String, String>,
    mapping: &Mapping,
    core_pairs: &[RepeatCorePair],
) -> BTreeMap<String, RepeatCorePair> {
    let mut assignments = BTreeMap::new();
    let mut used = HashSet::new();
    for name in &merged.order {
        if node_classes[name] != "repeat_node" {
            continue;
        }
        let mut best_index = None;
        let mut best_score = 0usize;
        for (index, pair) in core_pairs.iter().enumerate() {
            if used.contains(&index) {
                continue;
            }
            let score: usize = mapping.by_node[name]
                .all_hits
                .iter()
                .map(|hit| {
                    pair.copies
                        .iter()
                        .map(|copy| hit_interval_overlap(hit, copy))
                        .sum::<usize>()
                })
                .sum();
            if best_index.is_none() || score > best_score {
                best_index = Some(index);
                best_score = score;
            }
        }
        if let Some(index) = best_index {
            used.insert(index);
            assignments.insert(name.clone(), core_pairs[index].clone());
        }
    }
    assignments
}

fn build_repeat_core_gaps(
    repeat_assignments: &BTreeMap<String, RepeatCorePair>,
    reference_len: usize,
) -> Vec<RepeatGap> {
    let mut copies = Vec::new();
    for (node, pair) in repeat_assignments {
        for (copy_index, copy) in pair.copies.iter().enumerate() {
            copies.push((node.clone(), copy_index + 1, copy.clone()));
        }
    }
    copies.sort_by_key(|(_, _, copy)| copy.start);
    let mut gaps = Vec::new();
    if copies.is_empty() {
        return gaps;
    }
    for index in 0..copies.len() {
        let (left_node, left_copy_index, left_copy) = &copies[index];
        let (right_node, right_copy_index, right_copy) = &copies[(index + 1) % copies.len()];
        let parts = repeat_gap_parts(left_copy, right_copy, reference_len);
        gaps.push(RepeatGap {
            left_repeat_node: left_node.clone(),
            left_copy_index: *left_copy_index,
            right_repeat_node: right_node.clone(),
            right_copy_index: *right_copy_index,
            parts,
        });
    }
    gaps
}

fn repeat_gap_parts(
    left_copy: &RepeatCopy,
    right_copy: &RepeatCopy,
    reference_len: usize,
) -> Vec<(usize, usize)> {
    let start = left_copy.end + 1;
    let end = right_copy.start.saturating_sub(1);
    if start <= end {
        vec![(start, end)]
    } else {
        let mut parts = Vec::new();
        if start <= reference_len {
            parts.push((start, reference_len));
        }
        if end >= 1 {
            parts.push((1, end));
        }
        parts
    }
}

fn assign_single_copy_nodes_to_gaps(
    merged: &Gfa,
    node_classes: &BTreeMap<String, String>,
    mapping: &Mapping,
    gaps: &[RepeatGap],
) -> BTreeMap<String, RepeatGap> {
    let mut assignments = BTreeMap::new();
    let mut used = HashSet::new();
    for name in &merged.order {
        if node_classes[name] != "single_copy_node" {
            continue;
        }
        let mut best_index = None;
        let mut best_score = 0usize;
        for (index, gap) in gaps.iter().enumerate() {
            if used.contains(&index) {
                continue;
            }
            let score: usize = mapping.by_node[name]
                .all_hits
                .iter()
                .map(|hit| hit_overlap_with_parts(hit, &gap.parts))
                .sum();
            if best_index.is_none() || score > best_score {
                best_index = Some(index);
                best_score = score;
            }
        }
        if let Some(index) = best_index {
            used.insert(index);
            assignments.insert(name.clone(), gaps[index].clone());
        }
    }
    assignments
}

fn choose_repeat_copy_for_node(
    name: &str,
    pair: &RepeatCorePair,
    mapping: &Mapping,
) -> (RepeatCopy, char) {
    let mut best_copy = pair.copies[0].clone();
    let mut best_strand = pair.copies[0].strand;
    let mut best_score = 0usize;
    for hit in &mapping.by_node[name].all_hits {
        for copy in &pair.copies {
            let score = hit_interval_overlap(hit, copy);
            if score > best_score {
                best_score = score;
                best_copy = copy.clone();
                best_strand = hit.strand;
            }
        }
    }
    (best_copy, best_strand)
}

fn best_orientation_for_parts(name: &str, parts: &[(usize, usize)], mapping: &Mapping) -> char {
    let mut best_strand = '+';
    let mut best_score = 0usize;
    for hit in &mapping.by_node[name].all_hits {
        let score = hit_overlap_with_parts(hit, parts);
        if score > best_score {
            best_score = score;
            best_strand = hit.strand;
        }
    }
    best_strand
}

fn hit_interval_overlap(hit: &PafHit, copy: &RepeatCopy) -> usize {
    interval_overlap(hit.target_start + 1, hit.target_end, copy.start, copy.end)
}

fn hit_overlap_with_parts(hit: &PafHit, parts: &[(usize, usize)]) -> usize {
    parts
        .iter()
        .map(|(start, end)| interval_overlap(hit.target_start + 1, hit.target_end, *start, *end))
        .sum()
}

fn synthetic_hits_for_repeat_copies(
    name: &str,
    node_sequence: &str,
    record: &FastaRecord,
    copies: &[RepeatCopy],
) -> Vec<PafHit> {
    let mut hits = Vec::new();
    for copy in copies {
        let plus_sequence = record_interval_sequence(record, copy.start, copy.end, '+');
        let minus_sequence = reverse_complement(&plus_sequence);
        let strand = if node_sequence == plus_sequence {
            '+'
        } else if node_sequence == minus_sequence {
            '-'
        } else {
            copy.strand
        };
        hits.push(make_synthetic_hit(
            name,
            node_sequence.len(),
            1,
            node_sequence.len(),
            copy.start,
            copy.end,
            strand,
        ));
    }
    hits
}

fn synthetic_hits_for_parts(
    name: &str,
    node_len: usize,
    record: &FastaRecord,
    parts: &[(usize, usize)],
    strand: char,
) -> Vec<PafHit> {
    let mut hits = Vec::new();
    let mut query_pos = 1usize;
    let mut ordered_parts = parts.to_vec();
    if strand == '-' {
        ordered_parts.reverse();
    }
    for (start, end) in ordered_parts {
        let part_len = end - start + 1;
        hits.push(make_synthetic_hit(
            name,
            node_len,
            query_pos,
            query_pos + part_len - 1,
            start,
            end,
            strand,
        ));
        query_pos += part_len;
    }
    // Keep reference id in use so future refactoring does not accidentally change semantics.
    let _ = &record.id;
    hits
}

fn make_synthetic_hit(
    name: &str,
    node_len: usize,
    query_start: usize,
    query_end: usize,
    target_start: usize,
    target_end: usize,
    strand: char,
) -> PafHit {
    PafHit {
        query_name: name.to_string(),
        query_length: node_len,
        query_start: query_start - 1,
        query_end,
        strand,
        target_start: target_start - 1,
        target_end,
        matches: query_end - query_start + 1,
        block_length: query_end - query_start + 1,
        mapq: 60,
        primary: true,
    }
}

fn verified_sequence_from_parts(
    record: &FastaRecord,
    parts: &[(usize, usize)],
    strand: char,
) -> String {
    let mut sequence = String::new();
    for (start, end) in parts {
        sequence.push_str(&record_interval_sequence(record, *start, *end, '+'));
    }
    if strand == '-' {
        reverse_complement(&sequence)
    } else {
        sequence
    }
}

fn extract_verified_sequence(record: &FastaRecord, hit: &PafHit) -> String {
    let mut sequence = extract_circular(&record.sequence, hit.target_start, hit.target_end);
    if hit.strand == '-' {
        sequence = reverse_complement(&sequence);
    }
    sequence
}

fn record_interval_sequence(
    record: &FastaRecord,
    start: usize,
    end: usize,
    strand: char,
) -> String {
    let sequence = record
        .sequence
        .chars()
        .skip(start.saturating_sub(1))
        .take(end.saturating_sub(start) + 1)
        .collect::<String>();
    if strand == '-' {
        reverse_complement(&sequence)
    } else {
        sequence
    }
}

fn interval_overlap(start_a: usize, end_a: usize, start_b: usize, end_b: usize) -> usize {
    let start = start_a.max(start_b);
    let end = end_a.min(end_b);
    if end < start {
        0
    } else {
        end - start + 1
    }
}

fn sorted_pair(a: usize, b: usize) -> (usize, usize) {
    if a <= b {
        (a, b)
    } else {
        (b, a)
    }
}

fn write_self_blast_repeat_core_debug(
    path: &Path,
    pairs: &[RepeatCorePair],
) -> Result<(), OrgraftError> {
    let mut out = File::create(path)?;
    writeln!(out, "pair_index\tlength\tsource_blast_length\tsource_pident\tsource_mismatch\tsource_gapopen\tcopy_index\tstart\tend\tstrand\tboundary_extension_left\tboundary_extension_right")?;
    for (pair_index, pair) in pairs.iter().enumerate() {
        for (copy_index, copy) in pair.copies.iter().enumerate() {
            writeln!(
                out,
                "{}\t{}\t{}\t{:.3}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                pair_index + 1,
                pair.length,
                pair.source_blast_length,
                pair.source_pident,
                pair.source_mismatch,
                pair.source_gapopen,
                copy_index + 1,
                copy.start,
                copy.end,
                copy.strand,
                pair.boundary_extension_left,
                pair.boundary_extension_right
            )?;
        }
    }
    Ok(())
}

fn write_rebuild_extract(
    path: &Path,
    subgraph: &str,
    gfa: &Gfa,
    mapping: &Mapping,
    node_classes: &BTreeMap<String, String>,
    source_rows: &[VerifiedNodeSourceRow],
    coverage_rows: &[CoverageRow],
) -> Result<(), OrgraftError> {
    let degrees = gfa.degrees();
    let source_by_node: BTreeMap<String, VerifiedNodeSourceRow> = source_rows
        .iter()
        .map(|row| (row.node.clone(), row.clone()))
        .collect();
    let coverage_by_node: BTreeMap<String, CoverageRow> = coverage_rows
        .iter()
        .map(|row| (row.node.clone(), row.clone()))
        .collect();
    let mut rows: Vec<ExtractRow> = Vec::new();
    for name in &gfa.order {
        let segment = &gfa.segments[name];
        let Some(node) = mapping.by_node.get(name) else {
            rows.push(ExtractRow::unmapped(
                subgraph,
                name,
                segment.sequence.len(),
                node_classes[name].as_str(),
                degrees.get(name).copied().unwrap_or(0),
            ));
            continue;
        };
        for (copy_index, hit) in node.accepted_hits.iter().enumerate() {
            let interval = hit.folded(mapping.reference_len);
            let source = source_by_node.get(name);
            let coverage = coverage_by_node.get(name);
            let (node_start, node_end) = if hit.strand == '-' {
                (hit.query_end, hit.query_start + 1)
            } else {
                (hit.query_start + 1, hit.query_end)
            };
            let mut push_interval = |linear_start: usize, linear_end: usize, notes: &str| {
                rows.push(ExtractRow {
                    subgraph: subgraph.to_string(),
                    row_type: "node_copy".to_string(),
                    linear_start: linear_start.to_string(),
                    linear_end: linear_end.to_string(),
                    node: name.clone(),
                    node_start,
                    node_end,
                    strand: hit.strand,
                    node_length: segment.sequence.len(),
                    node_class: node_classes[name].clone(),
                    copy_index: copy_index + 1,
                    source_kind: source
                        .map(|row| row.source.clone())
                        .unwrap_or_else(|| "polished_fasta_minimap2".to_string()),
                    source_node_class: source
                        .map(|row| row.node_class.clone())
                        .unwrap_or_else(|| ".".to_string()),
                    source_length: source
                        .map(|row| row.length.to_string())
                        .unwrap_or_else(|| ".".to_string()),
                    source_parts: source
                        .map(|row| row.parts.clone())
                        .unwrap_or_else(|| ".".to_string()),
                    source_notes: source
                        .map(|row| row.notes.clone())
                        .unwrap_or_else(|| ".".to_string()),
                    depth_mean: optional_f64(coverage.and_then(|row| row.depth_mean)),
                    coverage_intervals: coverage
                        .map(|row| row.intervals.clone())
                        .unwrap_or_else(|| ".".to_string()),
                    coverage_source: coverage
                        .map(|row| row.coverage_source.clone())
                        .unwrap_or_else(|| ".".to_string()),
                    degree: degrees.get(name).copied().unwrap_or(0),
                    identity: format!("{:.6}", hit.identity()),
                    aligned_fraction: format!("{:.6}", hit.aligned_fraction()),
                    mapq: hit.mapq.to_string(),
                    notes: notes.to_string(),
                });
            };
            if interval.wraps {
                push_interval(
                    interval.start,
                    mapping.reference_len,
                    "wrapped_circular_interval",
                );
                push_interval(1, interval.end, "wrapped_circular_interval");
            } else {
                push_interval(interval.start, interval.end, ".");
            }
        }
        if node.accepted_hits.is_empty() {
            rows.push(ExtractRow::unmapped(
                subgraph,
                name,
                segment.sequence.len(),
                node_classes[name].as_str(),
                degrees.get(name).copied().unwrap_or(0),
            ));
        }
    }
    rows.sort_by(|a, b| {
        parse_sortable_usize(&a.linear_start)
            .cmp(&parse_sortable_usize(&b.linear_start))
            .then_with(|| {
                parse_sortable_usize(&a.linear_end).cmp(&parse_sortable_usize(&b.linear_end))
            })
            .then_with(|| natural_cmp(&a.node, &b.node))
            .then_with(|| a.copy_index.cmp(&b.copy_index))
    });
    let mut out = File::create(path)?;
    writeln!(out, "subgraph\trow_type\tlinear_start\tlinear_end\tnode\tnode_start\tnode_end\tstrand\tnode_length\tnode_class\tcopy_index\tsource_kind\tsource_node_class\tsource_length\tsource_parts\tsource_notes\tdepth_mean\tcoverage_intervals\tcoverage_source\tdegree\tidentity\taligned_fraction\tmapq\tnotes")?;
    for row in rows {
        writeln!(
            out,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            row.subgraph,
            row.row_type,
            row.linear_start,
            row.linear_end,
            row.node,
            row.node_start,
            row.node_end,
            row.strand,
            row.node_length,
            row.node_class,
            row.copy_index,
            row.source_kind,
            row.source_node_class,
            row.source_length,
            row.source_parts,
            row.source_notes,
            row.depth_mean,
            row.coverage_intervals,
            row.coverage_source,
            row.degree,
            row.identity,
            row.aligned_fraction,
            row.mapq,
            row.notes
        )?;
    }
    Ok(())
}

#[derive(Debug)]
struct ExtractRow {
    subgraph: String,
    row_type: String,
    linear_start: String,
    linear_end: String,
    node: String,
    node_start: usize,
    node_end: usize,
    strand: char,
    node_length: usize,
    node_class: String,
    copy_index: usize,
    source_kind: String,
    source_node_class: String,
    source_length: String,
    source_parts: String,
    source_notes: String,
    depth_mean: String,
    coverage_intervals: String,
    coverage_source: String,
    degree: usize,
    identity: String,
    aligned_fraction: String,
    mapq: String,
    notes: String,
}

impl ExtractRow {
    fn unmapped(
        subgraph: &str,
        node: &str,
        node_length: usize,
        node_class: &str,
        degree: usize,
    ) -> Self {
        Self {
            subgraph: subgraph.to_string(),
            row_type: "node_unmapped".to_string(),
            linear_start: ".".to_string(),
            linear_end: ".".to_string(),
            node: node.to_string(),
            node_start: 0,
            node_end: 0,
            strand: '.',
            node_length,
            node_class: node_class.to_string(),
            copy_index: 0,
            source_kind: ".".to_string(),
            source_node_class: ".".to_string(),
            source_length: ".".to_string(),
            source_parts: ".".to_string(),
            source_notes: "no_alignment".to_string(),
            depth_mean: ".".to_string(),
            coverage_intervals: ".".to_string(),
            coverage_source: ".".to_string(),
            degree,
            identity: ".".to_string(),
            aligned_fraction: ".".to_string(),
            mapq: ".".to_string(),
            notes: "no_alignment".to_string(),
        }
    }
}

fn parse_sortable_usize(value: &str) -> usize {
    value.parse::<usize>().unwrap_or(usize::MAX)
}

fn export_gfa_reference_images(
    paths: &OutputPaths,
    image_reference_fasta: &Path,
    gfa_editor_cli: &Path,
    soft_paths: &HashMap<String, PathBuf>,
) -> Vec<GfaImageExport> {
    let mut rows = Vec::new();
    for (format, output_path) in [
        ("pdf", paths.verified_pdf.as_path()),
        ("svg", paths.verified_svg.as_path()),
    ] {
        rows.push(run_gfa_editor_image(
            gfa_editor_cli,
            soft_paths,
            &paths.verified_gfa,
            output_path,
            image_reference_fasta,
            format,
        ));
    }
    rows
}

fn skipped_gfa_reference_images(paths: &OutputPaths, reason: &str) -> Vec<GfaImageExport> {
    [
        ("pdf", paths.verified_pdf.as_path()),
        ("svg", paths.verified_svg.as_path()),
    ]
    .into_iter()
    .map(|(format, output_path)| GfaImageExport {
        format: format.to_string(),
        output: output_path.to_path_buf(),
        command: ".".to_string(),
        status: "skipped_missing_gfa_editor_cli".to_string(),
        stdout: String::new(),
        stderr: reason.to_string(),
    })
    .collect()
}

fn write_run_report(
    path: &Path,
    out_dir: &Path,
    paths: &OutputPaths,
    options: &RebuildOptions,
    edited_gfa: &Path,
    polished_fasta: &Path,
    image_reference_fasta: Option<&Path>,
    validate_data_dir: Option<&Path>,
    merge_mode: &str,
    sequence_projection_method: &str,
    minimap2: &Path,
    blastn: &Path,
    gfa_editor_cli: Option<&Path>,
    mapping: &Mapping,
    image_exports: &[GfaImageExport],
    started: &Instant,
) -> Result<(), OrgraftError> {
    let mut out = File::create(path)?;
    writeln!(out, "section\tkey\tvalue")?;
    let mut rows = vec![
        ("run", "status", "completed".to_string()),
        ("run", "organelle", options.organelle.clone()),
        ("run", "subgraph", options.subgraph.clone()),
        ("run", "threads", options.threads.to_string()),
        (
            "run",
            "runtime_seconds",
            format!("{:.6}", started.elapsed().as_secs_f64()),
        ),
        ("input", "edited_gfa", edited_gfa.display().to_string()),
        (
            "input",
            "polished_fasta",
            polished_fasta.display().to_string(),
        ),
        (
            "input",
            "image_reference_fasta",
            image_reference_fasta
                .map(|path| path.display().to_string())
                .unwrap_or_else(|| ".".to_string()),
        ),
        (
            "input",
            "validate_data_dir",
            validate_data_dir
                .map(|path| path.display().to_string())
                .unwrap_or_else(|| ".".to_string()),
        ),
        (
            "output",
            "rebuild_gfa",
            relative_output_path(out_dir, &paths.verified_gfa),
        ),
        (
            "output",
            "rebuild_fasta",
            relative_output_path(out_dir, &paths.verified_fasta),
        ),
        (
            "output",
            "rebuild_nodes_fasta",
            relative_output_path(out_dir, &paths.verified_nodes_fasta),
        ),
        (
            "output",
            "extract_table",
            relative_output_path(out_dir, &paths.extract_report),
        ),
        (
            "output",
            "run_report",
            relative_output_path(out_dir, &paths.run_report),
        ),
        (
            "output",
            "result_stats",
            relative_output_path(out_dir, &paths.result_stats),
        ),
        (
            "file_description",
            "rebuild_fasta",
            "complete verified/polished linear sequence for the subgraph".to_string(),
        ),
        (
            "file_description",
            "rebuild_nodes_fasta",
            "FASTA records extracted from each S node in the rebuilt GFA".to_string(),
        ),
        ("method", "merge_mode", merge_mode.to_string()),
        ("method", "repeat_resolution", "disabled".to_string()),
        (
            "method",
            "sequence_projection",
            sequence_projection_method.to_string(),
        ),
        ("tool", "minimap2", minimap2.display().to_string()),
        ("tool", "blastn", blastn.display().to_string()),
        (
            "tool",
            "gfa_editor_cli",
            gfa_editor_cli
                .map(|path| path.display().to_string())
                .unwrap_or_else(|| {
                    if image_reference_fasta.is_some() {
                        "skipped_missing_gfa_editor_cli".to_string()
                    } else {
                        "skipped_no_image_reference".to_string()
                    }
                }),
        ),
        ("tool", "minimap2_command", mapping.command.clone()),
        (
            "tool",
            "minimap2_stderr",
            mapping.stderr.replace('\n', "\\n"),
        ),
    ];
    for row in image_exports.iter().filter(|row| row.status == "written") {
        let key = match row.format.as_str() {
            "pdf" => "rebuild_pdf",
            "svg" => "rebuild_svg",
            other => other,
        };
        rows.push(("output", key, relative_output_path(out_dir, &row.output)));
    }
    for (section, key, value) in rows {
        writeln!(out, "{section}\t{key}\t{}", value.replace('\n', "\\n"))?;
    }
    for row in image_exports {
        writeln!(
            out,
            "image\t{}_output\t{}",
            row.format,
            relative_output_path(out_dir, &row.output)
        )?;
        writeln!(out, "image\t{}_status\t{}", row.format, row.status)?;
        writeln!(
            out,
            "image\t{}_command\t{}",
            row.format,
            row.command.replace('\n', "\\n")
        )?;
        writeln!(
            out,
            "image\t{}_stdout\t{}",
            row.format,
            row.stdout.replace('\n', "\\n")
        )?;
        writeln!(
            out,
            "image\t{}_stderr\t{}",
            row.format,
            row.stderr.replace('\n', "\\n")
        )?;
    }
    Ok(())
}

fn write_result_stats(
    path: &Path,
    subgraph: &str,
    summary: &[(&str, String)],
    raw_stats: &GraphStats,
    merged_stats: &GraphStats,
    verified_stats: &GraphStats,
    consistency: &CoordinateConsistency,
    coverage_rows: &[CoverageRow],
    repeat_rows: &[RepeatPathSupportRow],
) -> Result<(), OrgraftError> {
    let mut out = File::create(path)?;
    writeln!(out, "section\tsubgraph\titem\tmetric\tvalue\textra")?;
    for (key, value) in summary {
        writeln!(
            out,
            "summary\t{subgraph}\t.\t{key}\t{}\t.",
            value.replace('\n', "\\n")
        )?;
    }
    for (name, stats) in [
        ("raw", raw_stats),
        ("merged", merged_stats),
        ("verified", verified_stats),
    ] {
        writeln!(
            out,
            "graph\t{subgraph}\t{name}\tsegments\t{}\t.",
            stats.segments
        )?;
        writeln!(out, "graph\t{subgraph}\t{name}\tlinks\t{}\t.", stats.links)?;
        writeln!(out, "graph\t{subgraph}\t{name}\tbp\t{}\t.", stats.bp)?;
    }
    writeln!(
        out,
        "consistency\t{subgraph}\tpolished\tlength\t{}\thash64={}",
        consistency.polished_length, consistency.polished_hash64
    )?;
    writeln!(
        out,
        "consistency\t{subgraph}\tsegments\tmapped\t{}\tunmapped={}",
        consistency.mapped_segments, consistency.unmapped_segments
    )?;
    writeln!(
        out,
        "consistency\t{subgraph}\tcoverage\tcovered_bases\t{}\tmulti_covered_bases={};gap_bases={};coverage_fraction={:.6};linear_tiling={}",
        consistency.covered_bases,
        consistency.multi_covered_bases,
        consistency.gap_bases,
        consistency.coverage_fraction,
        consistency.linear_tiling_status()
    )?;
    writeln!(
        out,
        "consistency\t{subgraph}\tgraph_segments\tverified\t{}\traw={};merged={}",
        consistency.verified_segments, consistency.raw_segments, consistency.merged_segments
    )?;
    for row in coverage_rows {
        writeln!(
            out,
            "node_depth\t{subgraph}\t{}\tdepth_mean\t{}\tnode_class={};covered_bases={};intervals={};source={};coverage_path={}",
            row.node,
            optional_f64(row.depth_mean),
            row.node_class,
            optional_usize(row.covered_bases),
            row.intervals,
            row.coverage_source,
            row.coverage_path
        )?;
    }
    for row in repeat_rows {
        writeln!(
            out,
            "repeat_path\t{subgraph}\t{}\tpath_support\t{}\tstatus={};method={};left={};right={};ratio={};left_edge={};right_edge={};left_support={};right_support={}",
            row.repeat_node,
            optional_f64(row.path_support),
            row.repeat_status,
            row.support_method,
            row.left_endpoint,
            row.right_endpoint,
            optional_f64(row.path_ratio),
            optional_usize(row.left_edge_index),
            optional_usize(row.right_edge_index),
            optional_f64(row.left_support),
            optional_f64(row.right_support)
        )?;
    }
    Ok(())
}

#[derive(Debug, Clone)]
struct CoverageRow {
    node: String,
    node_class: String,
    coverage_source: String,
    depth_mean: Option<f64>,
    depth_bases: Option<f64>,
    covered_bases: Option<usize>,
    interval_count: usize,
    intervals: String,
    coverage_path: String,
}

fn compute_node_remapped_coverage(
    mapping: &Mapping,
    node_classes: &BTreeMap<String, String>,
    validate_data_dir: Option<&Path>,
) -> Result<Vec<CoverageRow>, OrgraftError> {
    let coverage_path = validate_data_dir
        .map(|dir| dir.join("sv_coverage.tsv"))
        .unwrap_or_default();
    let depths = read_remapped_depths(&coverage_path)?;
    let mut rows = Vec::new();
    for (name, node) in &mapping.by_node {
        let intervals: Vec<(usize, usize)> = node
            .accepted_hits
            .iter()
            .map(|hit| (hit.target_start + 1, hit.target_end))
            .collect();
        let merged = merge_intervals(&intervals);
        if let Some((depth_mean, depth_bases, covered_bases)) =
            coverage_for_intervals(&depths, &merged)
        {
            rows.push(CoverageRow {
                node: name.clone(),
                node_class: node_classes[name].clone(),
                coverage_source: "remapped_FL_to_verified_fasta".to_string(),
                depth_mean: Some(depth_mean),
                depth_bases: Some(depth_bases),
                covered_bases: Some(covered_bases),
                interval_count: merged.len(),
                intervals: intervals_text(&merged),
                coverage_path: coverage_path.display().to_string(),
            });
        } else {
            rows.push(CoverageRow {
                node: name.clone(),
                node_class: node_classes[name].clone(),
                coverage_source: "remapped_FL_to_verified_fasta_unavailable".to_string(),
                depth_mean: None,
                depth_bases: None,
                covered_bases: None,
                interval_count: 0,
                intervals: ".".to_string(),
                coverage_path: coverage_path.display().to_string(),
            });
        }
    }
    Ok(rows)
}

fn read_remapped_depths(path: &Path) -> Result<Vec<f64>, OrgraftError> {
    if !path.exists() {
        return Ok(Vec::new());
    }
    let text = fs::read_to_string(path).map_err(|error| {
        OrgraftError::InvalidArgument(format!("cannot read {}: {error}", path.display()))
    })?;
    let mut lines = text.lines();
    let Some(header) = lines.next() else {
        return Ok(Vec::new());
    };
    let columns: Vec<&str> = header.split('\t').collect();
    let depth_index = columns
        .iter()
        .position(|column| *column == "fl_depth")
        .or_else(|| columns.iter().position(|column| *column == "depth"))
        .unwrap_or_else(|| columns.len().saturating_sub(1));
    let mut depths = vec![0.0];
    for line in lines {
        let fields: Vec<&str> = line.split('\t').collect();
        let value = fields
            .get(depth_index)
            .and_then(|value| value.parse::<f64>().ok())
            .unwrap_or(0.0);
        depths.push(value);
    }
    Ok(depths)
}

fn merge_intervals(intervals: &[(usize, usize)]) -> Vec<(usize, usize)> {
    let mut sorted: Vec<(usize, usize)> = intervals
        .iter()
        .copied()
        .filter(|(start, end)| start <= end)
        .collect();
    sorted.sort();
    let mut merged: Vec<(usize, usize)> = Vec::new();
    for (start, end) in sorted {
        if let Some(last) = merged.last_mut() {
            if start <= last.1 + 1 {
                last.1 = last.1.max(end);
                continue;
            }
        }
        merged.push((start, end));
    }
    merged
}

fn coverage_for_intervals(
    depths: &[f64],
    intervals: &[(usize, usize)],
) -> Option<(f64, f64, usize)> {
    if depths.is_empty() || intervals.is_empty() {
        return None;
    }
    let max_position = depths.len().saturating_sub(1);
    let mut depth_bases = 0.0;
    let mut covered_bases = 0usize;
    for (start, end) in intervals {
        let start = (*start).max(1);
        let end = (*end).min(max_position);
        if end < start {
            continue;
        }
        for value in &depths[start..=end] {
            depth_bases += *value;
        }
        covered_bases += end - start + 1;
    }
    if covered_bases == 0 {
        None
    } else {
        Some((
            depth_bases / covered_bases as f64,
            depth_bases,
            covered_bases,
        ))
    }
}

fn intervals_text(intervals: &[(usize, usize)]) -> String {
    if intervals.is_empty() {
        ".".to_string()
    } else {
        intervals
            .iter()
            .map(|(start, end)| format!("{start}-{end}"))
            .collect::<Vec<_>>()
            .join(",")
    }
}

fn annotate_gfa_with_remapped_coverage(gfa: &mut Gfa, rows: &[CoverageRow]) {
    let by_node: BTreeMap<&str, &CoverageRow> =
        rows.iter().map(|row| (row.node.as_str(), row)).collect();
    for name in &gfa.order {
        let Some(row) = by_node.get(name.as_str()) else {
            continue;
        };
        let Some(segment) = gfa.segments.get_mut(name) else {
            continue;
        };
        let mut tags = segment.tags.clone();
        if let (Some(depth_mean), Some(depth_bases)) = (row.depth_mean, row.depth_bases) {
            tags = replace_tags(
                &tags,
                &[
                    format!("DP:f:{depth_mean:.6}"),
                    "CM:Z:remapped_FL_to_verified_fasta".to_string(),
                    format!("RI:i:{}", row.interval_count),
                    format!("RB:f:{depth_bases:.3}"),
                ],
            );
        } else {
            tags = replace_tag(
                &tags,
                "CM:Z:remapped_FL_to_verified_fasta_unavailable".to_string(),
            );
        }
        segment.tags = tags;
    }
}

#[derive(Debug, Clone)]
struct RepeatPathSupportRow {
    repeat_node: String,
    repeat_status: String,
    support_method: String,
    left_edge_index: Option<usize>,
    right_edge_index: Option<usize>,
    left_endpoint: String,
    right_endpoint: String,
    left_support: Option<f64>,
    right_support: Option<f64>,
    path_support: Option<f64>,
    path_ratio: Option<f64>,
}

fn repeat_path_support_rows(
    gfa: &Gfa,
    _mapping: &Mapping,
    node_classes: &BTreeMap<String, String>,
    _record: &FastaRecord,
    _validate_data_dir: Option<&Path>,
) -> Result<Vec<RepeatPathSupportRow>, OrgraftError> {
    let mut rows = Vec::new();
    for (node, class) in node_classes {
        if class != "repeat_node" {
            continue;
        }
        let mut left = Vec::new();
        let mut right = Vec::new();
        for (index, link) in gfa.links.iter().enumerate() {
            if let Some(side) = link_endpoint_side(link, node) {
                let item = (
                    index + 1,
                    other_link_endpoint(link, node),
                    link_support(link),
                );
                if side == '-' {
                    left.push(item)
                } else {
                    right.push(item)
                }
            }
        }
        let status = if left.len() == 2 && right.len() == 2 {
            "unresolved"
        } else if left.len() <= 1 && right.len() <= 1 {
            "resolved"
        } else {
            "ambiguous"
        };
        let combos: Vec<_> = left
            .iter()
            .flat_map(|l| right.iter().map(move |r| (l, r, l.2.min(r.2))))
            .collect();
        let total: f64 = combos.iter().map(|(_, _, support)| *support).sum();
        if combos.is_empty() {
            rows.push(RepeatPathSupportRow {
                repeat_node: node.clone(),
                repeat_status: status.to_string(),
                support_method: "edge_rc_min_proxy".to_string(),
                left_edge_index: None,
                right_edge_index: None,
                left_endpoint: ".".to_string(),
                right_endpoint: ".".to_string(),
                left_support: None,
                right_support: None,
                path_support: None,
                path_ratio: None,
            });
        }
        for (l, r, support) in combos {
            let ratio = if total > 0.0 { support / total } else { 0.0 };
            rows.push(RepeatPathSupportRow {
                repeat_node: node.clone(),
                repeat_status: status.to_string(),
                support_method: "edge_rc_min_proxy".to_string(),
                left_edge_index: Some(l.0),
                right_edge_index: Some(r.0),
                left_endpoint: l.1.clone(),
                right_endpoint: r.1.clone(),
                left_support: Some(l.2),
                right_support: Some(r.2),
                path_support: Some(support),
                path_ratio: Some(ratio),
            });
        }
    }
    Ok(rows)
}

fn attach_repeat_path_support_paths(gfa: &mut Gfa, rows: &[RepeatPathSupportRow]) {
    gfa.other_lines
        .retain(|fields| fields.first().map(String::as_str) != Some("P"));
    let mut counters: BTreeMap<&str, usize> = BTreeMap::new();
    for row in rows {
        if row.left_edge_index.is_none() || row.right_edge_index.is_none() {
            continue;
        }
        let counter = counters.entry(&row.repeat_node).or_default();
        *counter += 1;
        let path_index = format!("p{counter}");
        let path = format!(
            "{},{},{}",
            row.left_endpoint,
            oriented_node(&row.repeat_node, '+'),
            row.right_endpoint
        );
        gfa.other_lines.push(vec![
            "P".to_string(),
            format!("repeat_{}_{}", row.repeat_node, path_index),
            path,
            "*,*".to_string(),
            "PT:Z:repeat_path_support".to_string(),
            format!("RN:Z:{}", row.repeat_node),
            format!("PI:Z:{path_index}"),
            format!("RS:Z:{}", row.repeat_status),
            format!("PM:Z:{}", row.support_method),
            format!("RC:f:{:.3}", row.path_support.unwrap_or(0.0)),
            format!("PR:f:{:.6}", row.path_ratio.unwrap_or(0.0)),
            format!("LE:Z:{}", row.left_endpoint),
            format!("RE:Z:{}", row.right_endpoint),
        ]);
    }
}

fn oriented_node(name: &str, orient: char) -> String {
    format!("{name}{orient}")
}

fn optional_usize(value: Option<usize>) -> String {
    value
        .map(|value| value.to_string())
        .unwrap_or_else(|| ".".to_string())
}

fn optional_f64(value: Option<f64>) -> String {
    value
        .map(|value| format!("{value:.6}"))
        .unwrap_or_else(|| ".".to_string())
}

#[derive(Debug)]
struct CoordinateConsistency {
    polished_length: usize,
    polished_hash64: String,
    mapped_segments: usize,
    unmapped_segments: usize,
    covered_bases: usize,
    multi_covered_bases: usize,
    gap_bases: usize,
    coverage_fraction: f64,
    raw_segments: usize,
    merged_segments: usize,
    verified_segments: usize,
}

impl CoordinateConsistency {
    fn new(
        record: &FastaRecord,
        raw: &Gfa,
        merged: &Gfa,
        verified: &Gfa,
        mapping: &Mapping,
    ) -> Self {
        let mut coverage = vec![0usize; record.sequence.len()];
        let mut mapped = 0usize;
        let mut unmapped = 0usize;
        for name in &merged.order {
            let node_mapping = &mapping.by_node[name];
            let hits = if node_mapping.accepted_hits.is_empty() {
                node_mapping.selected.iter().collect::<Vec<_>>()
            } else {
                node_mapping.accepted_hits.iter().collect::<Vec<_>>()
            };
            if hits.is_empty() {
                unmapped += 1;
            } else {
                mapped += 1;
            }
            for hit in hits {
                for pos in hit.target_start..hit.target_end {
                    coverage[pos % record.sequence.len()] += 1;
                }
            }
        }
        let covered = coverage.iter().filter(|value| **value > 0).count();
        let multi = coverage.iter().filter(|value| **value > 1).count();
        let gap_bases = record.sequence.len().saturating_sub(covered);
        Self {
            polished_length: record.sequence.len(),
            polished_hash64: stable_hash64(&record.sequence),
            mapped_segments: mapped,
            unmapped_segments: unmapped,
            covered_bases: covered,
            multi_covered_bases: multi,
            gap_bases,
            coverage_fraction: covered as f64 / record.sequence.len() as f64,
            raw_segments: raw.order.len(),
            merged_segments: merged.order.len(),
            verified_segments: verified.order.len(),
        }
    }

    fn linear_tiling_status(&self) -> &'static str {
        if self.gap_bases == 0 && self.multi_covered_bases == 0 {
            "PASS"
        } else {
            "WARN"
        }
    }
}

#[derive(Debug)]
struct GraphStats {
    segments: usize,
    links: usize,
    bp: usize,
}

fn graph_stats(gfa: &Gfa) -> GraphStats {
    GraphStats {
        segments: gfa.order.len(),
        links: gfa.links.len(),
        bp: gfa
            .order
            .iter()
            .map(|name| gfa.segments[name].sequence.len())
            .sum(),
    }
}

fn parse_paf(line: &str) -> Result<PafHit, OrgraftError> {
    let fields: Vec<&str> = line.split('\t').collect();
    if fields.len() < 12 {
        return Err(OrgraftError::InvalidArgument(format!(
            "malformed PAF line: {line}"
        )));
    }
    Ok(PafHit {
        query_name: fields[0].to_string(),
        query_length: parse_usize_value(fields[1], "PAF query length")?,
        query_start: parse_usize_value(fields[2], "PAF query start")?,
        query_end: parse_usize_value(fields[3], "PAF query end")?,
        strand: parse_strand(fields[4])?,
        target_start: parse_usize_value(fields[7], "PAF target start")?,
        target_end: parse_usize_value(fields[8], "PAF target end")?,
        matches: parse_usize_value(fields[9], "PAF matches")?,
        block_length: parse_usize_value(fields[10], "PAF block length")?,
        mapq: parse_usize_value(fields[11], "PAF mapq")?,
        primary: fields.iter().any(|field| *field == "tp:A:P"),
    })
}

fn sort_hits(hits: &mut [PafHit]) {
    hits.sort_by(|a, b| {
        b.query_end
            .saturating_sub(b.query_start)
            .cmp(&a.query_end.saturating_sub(a.query_start))
            .then_with(|| b.matches.cmp(&a.matches))
            .then_with(|| {
                b.identity()
                    .partial_cmp(&a.identity())
                    .unwrap_or(Ordering::Equal)
            })
            .then_with(|| b.primary.cmp(&a.primary))
            .then_with(|| b.mapq.cmp(&a.mapq))
            .then_with(|| a.target_start.cmp(&b.target_start))
    });
}

fn dedupe_hits(mut hits: Vec<PafHit>, reference_len: usize) -> Vec<PafHit> {
    sort_hits(&mut hits);
    let mut out: Vec<PafHit> = Vec::new();
    for hit in hits {
        if out.iter().any(|old| {
            let a = old.folded(reference_len);
            let b = hit.folded(reference_len);
            a.start == b.start && a.end == b.end && a.wraps == b.wraps && old.strand == hit.strand
        }) {
            continue;
        }
        out.push(hit);
    }
    out
}

fn read_fasta(path: &Path) -> Result<Vec<FastaRecord>, OrgraftError> {
    let text = fs::read_to_string(path).map_err(|error| {
        OrgraftError::InvalidArgument(format!("cannot read {}: {error}", path.display()))
    })?;
    let mut records = Vec::new();
    let mut header: Option<String> = None;
    let mut sequence = String::new();
    for line in text.lines() {
        if let Some(next_header) = line.strip_prefix('>') {
            if let Some(old_header) = header.replace(next_header.trim().to_string()) {
                records.push(FastaRecord {
                    id: fasta_id(&old_header),
                    header: old_header,
                    sequence: sequence.clone(),
                });
                sequence.clear();
            }
        } else if !line.trim().is_empty() {
            sequence.push_str(&line.trim().to_ascii_uppercase());
        }
    }
    if let Some(old_header) = header {
        records.push(FastaRecord {
            id: fasta_id(&old_header),
            header: old_header,
            sequence,
        });
    }
    if records.is_empty() {
        return Err(OrgraftError::InvalidArgument(format!(
            "{} contains no FASTA records",
            path.display()
        )));
    }
    Ok(records)
}

fn write_single_fasta(path: &Path, header: &str, sequence: &str) -> Result<(), OrgraftError> {
    let mut out = File::create(path)?;
    write_fasta_record(&mut out, header, sequence)
}

fn write_node_fasta(path: &Path, gfa: &Gfa) -> Result<(), OrgraftError> {
    write_node_fasta_from_order(path, gfa)
}

fn write_node_fasta_from_order(path: &Path, gfa: &Gfa) -> Result<(), OrgraftError> {
    let mut out = File::create(path)?;
    for name in &gfa.order {
        let segment = &gfa.segments[name];
        write_fasta_record(&mut out, &segment.name, &segment.sequence)?;
    }
    Ok(())
}

fn write_fasta_record<W: Write>(
    out: &mut W,
    header: &str,
    sequence: &str,
) -> Result<(), OrgraftError> {
    writeln!(out, ">{header}")?;
    for chunk in sequence.as_bytes().chunks(80) {
        writeln!(out, "{}", String::from_utf8_lossy(chunk))?;
    }
    Ok(())
}

fn extract_circular(sequence: &str, start: usize, end: usize) -> String {
    let bytes = sequence.as_bytes();
    let len = bytes.len();
    let mut out = Vec::with_capacity(end.saturating_sub(start));
    for pos in start..end {
        out.push(bytes[pos % len]);
    }
    String::from_utf8_lossy(&out).to_string()
}

fn reverse_complement(sequence: &str) -> String {
    sequence
        .bytes()
        .rev()
        .map(|base| match base {
            b'A' | b'a' => 'T',
            b'C' | b'c' => 'G',
            b'G' | b'g' => 'C',
            b'T' | b't' => 'A',
            b'N' | b'n' => 'N',
            other => other as char,
        })
        .collect()
}

fn complement_base(base: char) -> char {
    match base {
        'A' | 'a' => 'T',
        'C' | 'c' => 'G',
        'G' | 'g' => 'C',
        'T' | 't' => 'A',
        'N' | 'n' => 'N',
        other => other,
    }
}

fn flip_orient(orient: char) -> char {
    if orient == '+' {
        '-'
    } else {
        '+'
    }
}

fn parse_orient(value: &str, path: &Path, line: usize) -> Result<char, OrgraftError> {
    match value {
        "+" => Ok('+'),
        "-" => Ok('-'),
        _ => Err(OrgraftError::InvalidArgument(format!(
            "{}:{} invalid orientation `{value}`",
            path.display(),
            line
        ))),
    }
}

fn parse_strand(value: &str) -> Result<char, OrgraftError> {
    match value {
        "+" => Ok('+'),
        "-" => Ok('-'),
        _ => Err(OrgraftError::InvalidArgument(format!(
            "invalid strand `{value}`"
        ))),
    }
}

fn required_value<'a>(
    args: &'a [String],
    index: &mut usize,
    option: &str,
) -> Result<&'a str, OrgraftError> {
    *index += 1;
    args.get(*index)
        .map(String::as_str)
        .ok_or_else(|| OrgraftError::InvalidArgument(format!("{option} requires a value")))
}

fn parse_usize(value: &str, option: &str) -> Result<usize, OrgraftError> {
    value.parse::<usize>().map_err(|_| {
        OrgraftError::InvalidArgument(format!("{option} must be an integer, got `{value}`"))
    })
}

fn parse_subgraph_id(value: &str) -> Result<String, OrgraftError> {
    if value.is_empty() {
        return Err(OrgraftError::InvalidArgument(
            "--subgraph must not be empty".to_string(),
        ));
    }
    if !value
        .chars()
        .all(|ch| ch.is_ascii_alphanumeric() || ch == '_' || ch == '-')
    {
        return Err(OrgraftError::InvalidArgument(format!(
            "--subgraph may only contain ASCII letters, digits, `_`, and `-`, got `{value}`"
        )));
    }
    Ok(value.to_string())
}

fn parse_usize_value(value: &str, label: &str) -> Result<usize, OrgraftError> {
    value.parse::<usize>().map_err(|_| {
        OrgraftError::InvalidArgument(format!("{label} must be an integer, got `{value}`"))
    })
}

fn canonicalize_existing(path: &Path, label: &str) -> Result<PathBuf, OrgraftError> {
    fs::canonicalize(path).map_err(|error| {
        OrgraftError::InvalidArgument(format!("cannot read {label} {}: {error}", path.display()))
    })
}

fn remove_path_if_exists(path: &Path) -> Result<(), OrgraftError> {
    match fs::symlink_metadata(path) {
        Ok(metadata) if metadata.is_dir() => fs::remove_dir_all(path).map_err(OrgraftError::Io),
        Ok(_) => fs::remove_file(path).map_err(OrgraftError::Io),
        Err(error) if error.kind() == io::ErrorKind::NotFound => Ok(()),
        Err(error) => Err(OrgraftError::Io(error)),
    }
}

fn read_soft_paths(path: &Path) -> Result<HashMap<String, PathBuf>, OrgraftError> {
    if !path.exists() {
        return Ok(HashMap::new());
    }
    let text = fs::read_to_string(path).map_err(|error| {
        OrgraftError::InvalidArgument(format!("cannot read {}: {error}", path.display()))
    })?;
    let mut map = HashMap::new();
    for line in text.lines() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let fields: Vec<&str> = line.split_whitespace().collect();
        if fields.len() >= 2 {
            map.insert(fields[0].to_string(), PathBuf::from(fields[1]));
        }
    }
    Ok(map)
}

fn discover_validate_data_dir(polished_fasta: &Path) -> Option<PathBuf> {
    let file_dir = polished_fasta.parent()?;
    let round_dir = match file_dir.file_name().and_then(|name| name.to_str()) {
        Some("01.inputs") | Some("02.polish") => file_dir.parent()?,
        _ => return None,
    };
    let candidate = round_dir.join("03.validate/01.data");
    if candidate.is_dir() {
        Some(candidate)
    } else {
        None
    }
}

fn fasta_id(header: &str) -> String {
    header
        .split_whitespace()
        .next()
        .unwrap_or(header)
        .to_string()
}

fn replace_tag(tags: &[String], new_tag: String) -> Vec<String> {
    let tag_name = new_tag.split(':').next().unwrap_or("");
    let prefix = format!("{tag_name}:");
    let mut out: Vec<String> = tags
        .iter()
        .filter(|tag| !tag.starts_with(&prefix))
        .cloned()
        .collect();
    out.push(new_tag);
    out
}

fn replace_tags(tags: &[String], new_tags: &[String]) -> Vec<String> {
    let mut out = tags.to_vec();
    for tag in new_tags {
        out = replace_tag(&out, tag.clone());
    }
    out
}

fn tag_value(tags: &[String], tag_name: &str) -> Option<String> {
    let prefix = format!("{tag_name}:");
    tags.iter().find_map(|tag| {
        if !tag.starts_with(&prefix) {
            return None;
        }
        tag.splitn(3, ':').nth(2).map(str::to_string)
    })
}

fn numeric_tag(tags: &[String], tag_name: &str) -> Option<f64> {
    tag_value(tags, tag_name).and_then(|value| value.parse::<f64>().ok())
}

fn link_support(link: &Link) -> f64 {
    for tag_name in ["RC", "SK", "PA"] {
        if let Some(value) = numeric_tag(&link.tags, tag_name) {
            return value;
        }
    }
    0.0
}

fn link_endpoint_side(link: &Link, node: &str) -> Option<char> {
    if link.from_name == node {
        Some(if link.from_orient == '-' { '-' } else { '+' })
    } else if link.to_name == node {
        Some(if link.to_orient == '-' { '+' } else { '-' })
    } else {
        None
    }
}

fn other_link_endpoint(link: &Link, node: &str) -> String {
    if link.from_name == node {
        format!("{}{}", link.to_name, link.to_orient)
    } else if link.to_name == node {
        format!("{}{}", link.from_name, link.from_orient)
    } else {
        ".".to_string()
    }
}

fn natural_cmp(left: &str, right: &str) -> Ordering {
    natural_key(left).cmp(&natural_key(right))
}

fn natural_key(value: &str) -> Vec<NaturalPart> {
    let mut parts = Vec::new();
    let mut current = String::new();
    let mut in_digits = None;
    for ch in value.chars() {
        let digit = ch.is_ascii_digit();
        match in_digits {
            Some(state) if state == digit => current.push(ch),
            Some(state) => {
                parts.push(if state {
                    NaturalPart::Number(current.parse::<u64>().unwrap_or(0))
                } else {
                    NaturalPart::Text(current.clone())
                });
                current.clear();
                current.push(ch);
                in_digits = Some(digit);
            }
            None => {
                current.push(ch);
                in_digits = Some(digit);
            }
        }
    }
    if let Some(state) = in_digits {
        parts.push(if state {
            NaturalPart::Number(current.parse::<u64>().unwrap_or(0))
        } else {
            NaturalPart::Text(current)
        });
    }
    parts
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
enum NaturalPart {
    Text(String),
    Number(u64),
}

fn stable_hash64(sequence: &str) -> String {
    let mut hash: u64 = 0xcbf29ce484222325;
    for byte in sequence.as_bytes() {
        hash ^= *byte as u64;
        hash = hash.wrapping_mul(0x100000001b3);
    }
    format!("{hash:016x}")
}

fn temp_work_dir(label: &str) -> PathBuf {
    let stamp = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map(|value| value.as_nanos())
        .unwrap_or(0);
    std::env::temp_dir().join(format!("{label}-{}-{stamp}", process::id()))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_options() {
        let args = vec![
            "--edited-gfa".to_string(),
            "graph.gfa".to_string(),
            "--polished-fasta".to_string(),
            "final.fa".to_string(),
            "--image-reference-fasta".to_string(),
            "bait.fa".to_string(),
            "--out-dir".to_string(),
            "out".to_string(),
        ];
        let options = RebuildOptions::from_args(&args).unwrap();
        assert_eq!(options.edited_gfa, PathBuf::from("graph.gfa"));
        assert_eq!(options.polished_fasta, PathBuf::from("final.fa"));
        assert_eq!(
            options.image_reference_fasta,
            Some(PathBuf::from("bait.fa"))
        );
        assert_eq!(options.out_dir, PathBuf::from("out"));
        assert_eq!(options.subgraph, "subgraph_001");
    }

    #[test]
    fn rebuild_out_dir_defaults_to_results_rebuild() {
        let args = vec![
            "--edited-gfa".to_string(),
            "graph.gfa".to_string(),
            "--polished-fasta".to_string(),
            "final.fa".to_string(),
        ];
        let options = RebuildOptions::from_args(&args).unwrap();
        assert_eq!(options.out_dir, PathBuf::from("results/rebuild"));
        assert_eq!(options.image_reference_fasta, None);
    }

    #[test]
    fn image_reference_fasta_is_optional() {
        let args = vec![
            "--edited-gfa".to_string(),
            "graph.gfa".to_string(),
            "--polished-fasta".to_string(),
            "final.fa".to_string(),
        ];
        let options = RebuildOptions::from_args(&args).unwrap();
        assert_eq!(options.image_reference_fasta, None);
    }

    #[test]
    fn discover_validate_data_dir_uses_workflow_round_directory() {
        let base = std::env::temp_dir().join(format!(
            "orgraft-rebuild-round-suffix-{}",
            std::process::id()
        ));
        let round_dir = base.join("04.polish/mito/subgraph_001/round_2");
        let data_dir = round_dir.join("03.validate/01.data");
        fs::create_dir_all(&data_dir).unwrap();
        let polished = round_dir.join("01.inputs/linear_subgraph.round_2.fasta");
        fs::create_dir_all(polished.parent().unwrap()).unwrap();
        fs::write(&polished, ">x\nACGT\n").unwrap();

        assert_eq!(discover_validate_data_dir(&polished), Some(data_dir));
        let _ = fs::remove_dir_all(base);
    }

    #[test]
    fn gfa_image_export_failure_is_recorded_without_error() {
        let row = run_gfa_editor_image(
            Path::new("/definitely/missing/gfa_editor_cli"),
            &HashMap::new(),
            Path::new("input.gfa"),
            Path::new("output.pdf"),
            Path::new("reference.fa"),
            "pdf",
        );

        assert_eq!(row.format, "pdf");
        assert_eq!(row.status, "failed_to_run");
        assert!(row.stderr.contains("failed to run GFA_Editor image export"));
    }

    #[test]
    fn skipped_gfa_reference_images_records_both_formats() {
        let paths = OutputPaths::new(Path::new("out"), "subgraph_001");
        let rows = skipped_gfa_reference_images(&paths, "missing cli");

        assert_eq!(rows.len(), 2);
        assert_eq!(rows[0].format, "pdf");
        assert_eq!(rows[0].status, "skipped_missing_gfa_editor_cli");
        assert_eq!(rows[0].stderr, "missing cli");
        assert_eq!(rows[1].format, "svg");
        assert_eq!(rows[1].status, "skipped_missing_gfa_editor_cli");
    }

    #[test]
    fn help_places_subgraph_in_inputs_and_uses_layout_summary() {
        assert!(HELP.contains(
            "Usage:\n  orgraft rebuild --organelle NAME --subgraph ID --edited-gfa FILE --polished-fasta FILE [options]"
        ));
        assert!(!HELP.contains(
            "Usage:\n  orgraft rebuild --organelle NAME --subgraph ID --edited-gfa FILE --polished-fasta FILE --out-dir DIR"
        ));
        assert!(HELP.contains(
            "  --out-dir DIR                 rebuild root output directory [results/rebuild]\n  --force                       replace existing output directory"
        ));
        assert!(HELP.contains(
            "Inputs:\n  --organelle NAME              organelle name for this rebuild run [mito]\n  --subgraph ID                 subgraph/ring id [subgraph_001]"
        ));
        assert!(HELP.contains(
            "Additional Parameters:\n  --threads N                   minimap2 threads for node-to-polished projection [4]\n  --image-reference-fasta FILE  reference FASTA for graph colouring; enables PDF/SVG export"
        ));
        assert!(HELP.contains(
            "Layout: OUT/SUBGRAPH/rebuild_SUBGRAPH* plus OUT/logs/*.tsv (pdf/svg need reference)"
        ));
        assert!(!HELP.contains("DIR/logs/manifest.tsv"));
        assert!(!HELP.contains("--genome"));
    }

    #[test]
    fn extracts_circular_sequence() {
        assert_eq!(extract_circular("ABCDEF", 4, 8), "EFAB");
    }

    #[test]
    fn consistency_uses_all_accepted_node_copies_for_linear_tiling() {
        let record = test_record("ACGTACGTAA");
        let gfa = test_gfa(&[("repeat", "AC"), ("single", "GTGTAA")]);
        let mut by_node = BTreeMap::new();
        let repeat_hits = vec![
            make_synthetic_hit("repeat", 2, 1, 2, 1, 2, '+'),
            make_synthetic_hit("repeat", 2, 1, 2, 5, 6, '+'),
        ];
        let single_hits = vec![
            make_synthetic_hit("single", 6, 1, 2, 3, 4, '+'),
            make_synthetic_hit("single", 6, 3, 6, 7, 10, '+'),
        ];
        by_node.insert("repeat".to_string(), test_mapping(repeat_hits));
        by_node.insert("single".to_string(), test_mapping(single_hits));
        let mapping = Mapping {
            by_node,
            reference_len: record.sequence.len(),
            command: String::new(),
            stderr: String::new(),
        };

        let consistency = CoordinateConsistency::new(&record, &gfa, &gfa, &gfa, &mapping);

        assert_eq!(consistency.covered_bases, 10);
        assert_eq!(consistency.gap_bases, 0);
        assert_eq!(consistency.multi_covered_bases, 0);
        assert_eq!(consistency.linear_tiling_status(), "PASS");
        assert_eq!(consistency.coverage_fraction, 1.0);
    }

    #[test]
    fn consistency_reports_gaps_and_overlaps() {
        let record = test_record("ACGTAC");
        let gfa = test_gfa(&[("left", "ACG"), ("right", "GT")]);
        let mut by_node = BTreeMap::new();
        by_node.insert(
            "left".to_string(),
            test_mapping(vec![make_synthetic_hit("left", 3, 1, 3, 1, 3, '+')]),
        );
        by_node.insert(
            "right".to_string(),
            test_mapping(vec![make_synthetic_hit("right", 2, 1, 2, 3, 4, '+')]),
        );
        let mapping = Mapping {
            by_node,
            reference_len: record.sequence.len(),
            command: String::new(),
            stderr: String::new(),
        };

        let consistency = CoordinateConsistency::new(&record, &gfa, &gfa, &gfa, &mapping);

        assert_eq!(consistency.covered_bases, 4);
        assert_eq!(consistency.gap_bases, 2);
        assert_eq!(consistency.multi_covered_bases, 1);
        assert_eq!(consistency.linear_tiling_status(), "WARN");
    }

    #[test]
    fn output_paths_use_subgraph_id() {
        let paths = OutputPaths::new(Path::new("out"), "subgraph_007");
        assert_eq!(
            paths.verified_gfa,
            PathBuf::from("out/subgraph_007/rebuild_subgraph_007.gfa")
        );
        assert_eq!(
            paths.verified_fasta,
            PathBuf::from("out/subgraph_007/rebuild_subgraph_007.fasta")
        );
        assert_eq!(
            paths.verified_nodes_fasta,
            PathBuf::from("out/subgraph_007/rebuild_subgraph_007_nodes.fasta")
        );
        assert_eq!(
            paths.verified_pdf,
            PathBuf::from("out/subgraph_007/rebuild_subgraph_007.pdf")
        );
        assert_eq!(
            paths.verified_svg,
            PathBuf::from("out/subgraph_007/rebuild_subgraph_007.svg")
        );
    }

    #[test]
    fn rejects_unsafe_subgraph_id() {
        assert!(parse_subgraph_id("../subgraph_001").is_err());
        assert!(parse_subgraph_id("subgraph/001").is_err());
        assert!(parse_subgraph_id("subgraph_001").is_ok());
    }

    fn test_record(sequence: &str) -> FastaRecord {
        FastaRecord {
            header: "test".to_string(),
            id: "test".to_string(),
            sequence: sequence.to_string(),
        }
    }

    fn test_gfa(nodes: &[(&str, &str)]) -> Gfa {
        let mut segments = BTreeMap::new();
        let mut order = Vec::new();
        for (name, sequence) in nodes {
            order.push((*name).to_string());
            segments.insert(
                (*name).to_string(),
                Segment {
                    name: (*name).to_string(),
                    sequence: (*sequence).to_string(),
                    tags: Vec::new(),
                },
            );
        }
        Gfa {
            headers: Vec::new(),
            segments,
            order,
            links: Vec::new(),
            other_lines: Vec::new(),
        }
    }

    fn test_mapping(hits: Vec<PafHit>) -> NodeMapping {
        NodeMapping {
            selected: hits.first().cloned(),
            accepted_hits: hits.clone(),
            all_hits: hits,
        }
    }
}
