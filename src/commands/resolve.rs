use std::cmp::Ordering;
use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet};
use std::fs::{self, File};
use std::io::{self, BufReader, Write};
use std::path::{Path, PathBuf};
use std::process::{self, Command, Stdio};
use std::time::{Instant, SystemTime, UNIX_EPOCH};

use crate::commands::shared::{print_contract, CommandContract};
use crate::error::OrgraftError;
use crate::topology::{analyze_gfa, TopologyReport};

const HELP: &str = r#"orgraft resolve

Resolve a checked draft GFA into reference-oriented graph and FASTA products.

Usage:
  orgraft resolve --checked-draft-gfa FILE --reference FILE [options]

Inputs:
  --checked-draft-gfa FILE       checked draft GFA
  --reference FILE               linear reference FASTA for automatic rotation
  --pre-rotated-reference FILE   overrides --reference and skips internal rotation
  --organelle NAME               output subdirectory name
  --soft-paths FILE              tool paths file [soft_paths.txt]

Outputs:
  --out-dir DIR                  resolve output directory [resolve_gfa]
  --force                        replace an existing output directory

Additional Parameters:
  --gfa-editor-mode MODE         rust|cli repeat candidate selection engine [rust]
  --max-states N                 auto-repeat search state limit [5000]
  --max-candidates N             auto-repeat final candidate limit [100]

Output layout: DIR[/NAME]/{logs,graph,fasta}/.
"#;

const MIN_NON_REPEAT_LEN: usize = 5000;
const DEFAULT_AUTO_REPEAT_MAX_STATES: usize = 5000;
const DEFAULT_AUTO_REPEAT_MAX_CANDIDATES: usize = 100;

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

    let options = ResolveOptions::from_args(args)?;
    run_resolve(&options)
}

fn run_resolve(options: &ResolveOptions) -> Result<(), OrgraftError> {
    let started = Instant::now();
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

    fs::create_dir_all(&options.out_dir)?;
    let paths = OutputPaths::new(&options.out_dir);
    fs::create_dir_all(&paths.logs_dir)?;
    fs::create_dir_all(&paths.graph_dir)?;
    fs::create_dir_all(&paths.fasta_dir)?;

    let checked_gfa = canonicalize_existing(&options.checked_draft_gfa)?;
    let reference_input = options.reference.canonicalize()?;

    let raw_graph = GfaGraph::read(&checked_gfa)?;
    let topology = build_topology_audit(&checked_gfa, &raw_graph)?;

    let soft_paths = read_soft_paths(&options.soft_paths)?;
    let blastn = require_tool(&soft_paths, "blastn")?;
    let gfa_editor_cli = resolve_gfa_editor_cli(&soft_paths, options.gfa_editor_mode)?;
    let reference_state =
        prepare_reference(&reference_input, &paths.reference_rotated_fasta, &blastn)?;

    let rust_prepare = prepare_graph_in_rust(&raw_graph, &paths.merged_unresolved_gfa)?;
    let merged_graph = GfaGraph::read(&paths.merged_unresolved_gfa)?;
    let split_component_gfas = write_component_gfas(&merged_graph, &paths.graph_dir)?;
    let draft_subgraph_fasta = temp_file_path("orgraft-resolve-subgraphs", "fasta");
    let repeat_resolution = build_reference_ordered_subgraph_fasta(
        &reference_state.rotated_fasta,
        &merged_graph,
        &draft_subgraph_fasta,
        &blastn,
        options.max_states,
        options.max_candidates,
        options.gfa_editor_mode,
        gfa_editor_cli.as_deref(),
    )?;
    let alignment = align_fasta_to_reference_coordinates(
        &reference_state.rotated_fasta,
        &draft_subgraph_fasta,
        &repeat_resolution,
        &paths.final_fasta,
        &blastn,
    )?;
    let _ = fs::remove_file(&draft_subgraph_fasta);

    write_id_map(&paths.id_map, &repeat_resolution)?;
    write_resolve_report(
        &paths.report,
        options,
        &checked_gfa,
        &paths,
        &topology,
        &reference_state,
        &rust_prepare,
        &split_component_gfas,
        &repeat_resolution,
        &alignment,
        started.elapsed().as_secs_f64(),
    )?;
    write_resolve_details(
        &paths.details,
        options,
        &checked_gfa,
        &paths,
        &topology,
        &reference_state,
        &rust_prepare,
        &split_component_gfas,
        &repeat_resolution,
        &alignment,
        started.elapsed().as_secs_f64(),
    )?;
    println!("Wrote {}", paths.report.display());
    println!("Wrote {}", paths.final_fasta.display());
    Ok(())
}

fn contract() -> CommandContract {
    CommandContract {
        command: "resolve",
        origin: "checked draft graph resolution workflow",
        purpose: "resolve checked draft GFA files into reference-oriented FASTA, graph reports, and run details",
        inputs: &[
            "--checked-draft-gfa FILE",
            "--reference FILE for automatic reference rotation, or --pre-rotated-reference FILE to skip internal rotation",
            "--organelle NAME optional output subdirectory name; no algorithm effect",
            "soft_paths.txt containing blastn and optional gfa_editor_cli",
        ],
        outputs: &[
            "logs/resolve_report.md",
            "logs/id_map.tsv",
            "logs/resolve_details.tsv",
            "graph/merged_unresolved.gfa",
            "graph/merged_unresolved_subgraph_*.gfa split component GFAs",
            "fasta/rotated_reference.fasta",
            "fasta/resolved_subgraphs.fasta",
        ],
        notes: &[
            "Rust owns GFA/FASTA parsing, topology/component reports, conservative unambiguous compaction, reference-oriented FASTA assembly, BLAST checks, and coordinate re-orientation",
            "gfa_editor_cli is optional for upstream auto-repeat candidate generation/selection in cli mode",
        ],
    }
}

#[derive(Debug, Clone)]
struct ResolveOptions {
    checked_draft_gfa: PathBuf,
    reference: ReferenceInput,
    soft_paths: PathBuf,
    out_dir: PathBuf,
    organelle: Option<String>,
    force: bool,
    gfa_editor_mode: GfaEditorMode,
    max_states: usize,
    max_candidates: usize,
}

#[derive(Debug, Clone)]
enum ReferenceInput {
    AutoRotate(PathBuf),
    PreRotated(PathBuf),
}

impl ReferenceInput {
    fn canonicalize(&self) -> Result<Self, OrgraftError> {
        match self {
            Self::AutoRotate(path) => Ok(Self::AutoRotate(canonicalize_existing(path)?)),
            Self::PreRotated(path) => Ok(Self::PreRotated(canonicalize_existing(path)?)),
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum GfaEditorMode {
    Rust,
    Cli,
}

impl GfaEditorMode {
    fn parse(value: &str) -> Result<Self, OrgraftError> {
        match value {
            "rust" => Ok(Self::Rust),
            "cli" | "gfa-editor" | "gfa_editor" => Ok(Self::Cli),
            other => Err(OrgraftError::InvalidArgument(format!(
                "unknown --gfa-editor-mode `{other}`; expected rust or cli"
            ))),
        }
    }

    fn as_str(self) -> &'static str {
        match self {
            Self::Rust => "rust",
            Self::Cli => "cli",
        }
    }
}

#[derive(Debug, Clone)]
struct OutputPaths {
    logs_dir: PathBuf,
    graph_dir: PathBuf,
    fasta_dir: PathBuf,
    report: PathBuf,
    id_map: PathBuf,
    details: PathBuf,
    final_fasta: PathBuf,
    merged_unresolved_gfa: PathBuf,
    reference_rotated_fasta: PathBuf,
}

impl OutputPaths {
    fn new(out_dir: &Path) -> Self {
        let logs_dir = out_dir.join("logs");
        let graph_dir = out_dir.join("graph");
        let fasta_dir = out_dir.join("fasta");
        Self {
            logs_dir: logs_dir.clone(),
            graph_dir: graph_dir.clone(),
            fasta_dir: fasta_dir.clone(),
            report: logs_dir.join("resolve_report.md"),
            id_map: logs_dir.join("id_map.tsv"),
            details: logs_dir.join("resolve_details.tsv"),
            final_fasta: fasta_dir.join("resolved_subgraphs.fasta"),
            merged_unresolved_gfa: graph_dir.join("merged_unresolved.gfa"),
            reference_rotated_fasta: fasta_dir.join("rotated_reference.fasta"),
        }
    }
}

impl ResolveOptions {
    fn from_args(args: &[String]) -> Result<Self, OrgraftError> {
        let mut checked_draft_gfa = None;
        let mut reference = None;
        let mut pre_rotated_reference = None;
        let mut soft_paths = PathBuf::from("soft_paths.txt");
        let mut out_dir: Option<PathBuf> = None;
        let mut organelle: Option<String> = None;
        let mut force = false;
        let mut gfa_editor_mode = GfaEditorMode::Rust;
        let mut max_states = DEFAULT_AUTO_REPEAT_MAX_STATES;
        let mut max_candidates = DEFAULT_AUTO_REPEAT_MAX_CANDIDATES;

        let mut index = 0usize;
        while index < args.len() {
            let arg = &args[index];
            match arg.as_str() {
                "--checked-draft-gfa" | "--gfa" | "--input" => {
                    checked_draft_gfa = Some(PathBuf::from(required_value(args, &mut index, arg)?));
                }
                "--reference" => {
                    reference = Some(PathBuf::from(required_value(args, &mut index, arg)?));
                }
                "--pre-rotated-reference" => {
                    pre_rotated_reference =
                        Some(PathBuf::from(required_value(args, &mut index, arg)?));
                }
                "--soft-paths" => {
                    soft_paths = PathBuf::from(required_value(args, &mut index, arg)?);
                }
                "--out-dir" => {
                    out_dir = Some(PathBuf::from(required_value(args, &mut index, arg)?));
                }
                "--organelle" => {
                    organelle = Some(parse_output_label(required_value(args, &mut index, arg)?)?);
                }
                "--force" => {
                    force = true;
                }
                "--gfa-editor-mode" => {
                    gfa_editor_mode = GfaEditorMode::parse(required_value(args, &mut index, arg)?)?;
                }
                "--use-gfa-editor-cli" => {
                    gfa_editor_mode = GfaEditorMode::Cli;
                }
                "--no-gfa-editor-cli" => {
                    gfa_editor_mode = GfaEditorMode::Rust;
                }
                "--max-states" => {
                    max_states = parse_usize(required_value(args, &mut index, arg)?, arg)?;
                }
                "--max-candidates" => {
                    max_candidates = parse_usize(required_value(args, &mut index, arg)?, arg)?;
                }
                other => {
                    return Err(OrgraftError::InvalidArgument(format!(
                        "unknown orgraft resolve option `{other}`"
                    )));
                }
            }
            index += 1;
        }

        let checked_draft_gfa = checked_draft_gfa.ok_or_else(|| {
            OrgraftError::InvalidArgument("missing --checked-draft-gfa FILE".to_string())
        })?;
        let reference = match (reference, pre_rotated_reference) {
            (Some(reference), None) => ReferenceInput::AutoRotate(reference),
            (None, Some(pre_rotated_reference)) => {
                ReferenceInput::PreRotated(pre_rotated_reference)
            }
            (Some(_), Some(pre_rotated_reference)) => {
                ReferenceInput::PreRotated(pre_rotated_reference)
            }
            (None, None) => {
                return Err(OrgraftError::InvalidArgument(
                    "missing --reference FILE (or --pre-rotated-reference FILE)".to_string(),
                ));
            }
        };
        let mut out_dir = out_dir.unwrap_or_else(|| PathBuf::from("resolve_gfa"));
        if let Some(label) = &organelle {
            out_dir = out_dir.join(label);
        }
        Ok(Self {
            checked_draft_gfa,
            reference,
            soft_paths,
            out_dir,
            organelle,
            force,
            gfa_editor_mode,
            max_states: max_states.max(1),
            max_candidates: max_candidates.max(1),
        })
    }
}

fn required_value<'a>(
    args: &'a [String],
    index: &mut usize,
    name: &str,
) -> Result<&'a str, OrgraftError> {
    *index += 1;
    args.get(*index)
        .map(String::as_str)
        .ok_or_else(|| OrgraftError::InvalidArgument(format!("missing value for {name}")))
}

fn parse_usize(value: &str, name: &str) -> Result<usize, OrgraftError> {
    value.parse::<usize>().map_err(|_| {
        OrgraftError::InvalidArgument(format!("{name} expects a positive integer, got `{value}`"))
    })
}

fn parse_output_label(value: &str) -> Result<String, OrgraftError> {
    let label = value.trim();
    let is_simple_name = !label.is_empty()
        && label != "."
        && label != ".."
        && label
            .chars()
            .all(|ch| ch.is_ascii_alphanumeric() || ch == '-' || ch == '_' || ch == '.');
    if is_simple_name {
        Ok(label.to_string())
    } else {
        Err(OrgraftError::InvalidArgument(format!(
            "--organelle expects a simple output label, got `{value}`"
        )))
    }
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
struct GfaGraph {
    headers: Vec<Vec<String>>,
    segment_order: Vec<String>,
    segments: BTreeMap<String, Segment>,
    links: Vec<Link>,
    other_lines: Vec<Vec<String>>,
}

impl GfaGraph {
    fn read(path: &Path) -> Result<Self, OrgraftError> {
        let text = fs::read_to_string(path).map_err(|error| {
            OrgraftError::InvalidArgument(format!("cannot read {}: {error}", path.display()))
        })?;
        Self::parse(&text, path)
    }

    fn parse(text: &str, path: &Path) -> Result<Self, OrgraftError> {
        let mut headers = Vec::new();
        let mut segment_order = Vec::new();
        let mut segments = BTreeMap::new();
        let mut links = Vec::new();
        let mut other_lines = Vec::new();

        for (index, line) in text.lines().enumerate() {
            if line.trim().is_empty() || line.starts_with('#') {
                continue;
            }
            let fields = line.split('\t').map(str::to_string).collect::<Vec<_>>();
            match fields.first().map(String::as_str) {
                Some("H") => headers.push(fields),
                Some("S") => {
                    if fields.len() < 3 {
                        return Err(OrgraftError::InvalidArgument(format!(
                            "{}:{} malformed GFA segment line",
                            path.display(),
                            index + 1
                        )));
                    }
                    if fields[2] == "*" {
                        return Err(OrgraftError::InvalidArgument(format!(
                            "{}:{} segment `{}` has `*` sequence; resolve needs concrete sequences",
                            path.display(),
                            index + 1,
                            fields[1]
                        )));
                    }
                    let segment = Segment {
                        name: fields[1].clone(),
                        sequence: fields[2].clone(),
                        tags: fields[3..].to_vec(),
                    };
                    if !segments.contains_key(&segment.name) {
                        segment_order.push(segment.name.clone());
                    }
                    segments.insert(segment.name.clone(), segment);
                }
                Some("L") => {
                    if fields.len() < 6 {
                        return Err(OrgraftError::InvalidArgument(format!(
                            "{}:{} malformed GFA link line",
                            path.display(),
                            index + 1
                        )));
                    }
                    links.push(Link {
                        from_name: fields[1].clone(),
                        from_orient: parse_orient(&fields[2], path, index + 1)?,
                        to_name: fields[3].clone(),
                        to_orient: parse_orient(&fields[4], path, index + 1)?,
                        overlap: fields[5].clone(),
                        tags: fields[6..].to_vec(),
                    });
                }
                _ => other_lines.push(fields),
            }
        }

        if segments.is_empty() {
            return Err(OrgraftError::InvalidArgument(format!(
                "{} contains no GFA segments",
                path.display()
            )));
        }

        Ok(Self {
            headers,
            segment_order,
            segments,
            links,
            other_lines,
        })
    }

    fn write(&self, path: &Path) -> Result<(), OrgraftError> {
        if let Some(parent) = path.parent() {
            fs::create_dir_all(parent)?;
        }
        let mut out = File::create(path)?;
        for fields in &self.headers {
            writeln!(out, "{}", fields.join("\t"))?;
        }
        for segment in self.ordered_segments() {
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

    fn write_fasta(&self, path: &Path) -> Result<(), OrgraftError> {
        if let Some(parent) = path.parent() {
            fs::create_dir_all(parent)?;
        }
        let mut out = File::create(path)?;
        for segment in self.ordered_segments() {
            write_fasta_record(&mut out, &segment.name, &segment.sequence)?;
        }
        Ok(())
    }

    fn ordered_segments(&self) -> Vec<&Segment> {
        let mut seen = HashSet::new();
        let mut ordered = Vec::new();
        for name in &self.segment_order {
            if let Some(segment) = self.segments.get(name) {
                seen.insert(name.clone());
                ordered.push(segment);
            }
        }
        for (name, segment) in &self.segments {
            if seen.insert(name.clone()) {
                ordered.push(segment);
            }
        }
        ordered
    }

    fn ordered_segment_names(&self) -> Vec<String> {
        self.ordered_segments()
            .into_iter()
            .map(|segment| segment.name.clone())
            .collect()
    }

    fn subgraph(&self, node_ids: &[String]) -> Self {
        let node_set = node_ids.iter().cloned().collect::<HashSet<_>>();
        let mut segment_order = self
            .segment_order
            .iter()
            .filter(|name| node_set.contains(*name))
            .cloned()
            .collect::<Vec<_>>();
        let segments = self
            .segments
            .iter()
            .filter(|(name, _)| node_set.contains(*name))
            .map(|(name, segment)| (name.clone(), segment.clone()))
            .collect::<BTreeMap<_, _>>();
        let links = self
            .links
            .iter()
            .filter(|link| node_set.contains(&link.from_name) && node_set.contains(&link.to_name))
            .cloned()
            .collect::<Vec<_>>();
        for name in segments.keys() {
            if !segment_order.contains(name) {
                segment_order.push(name.clone());
            }
        }
        Self {
            headers: self.headers.clone(),
            segment_order,
            segments,
            links,
            other_lines: Vec::new(),
        }
    }

    fn component_report(&self) -> ComponentReport {
        let mut names = self.segments.keys().cloned().collect::<Vec<_>>();
        names.sort_by(|left, right| natural_cmp(left, right));
        let index_by_name = names
            .iter()
            .enumerate()
            .map(|(index, name)| (name.clone(), index))
            .collect::<HashMap<_, _>>();
        let mut dsu = Dsu::new(names.len());
        for link in &self.links {
            if let (Some(left), Some(right)) = (
                index_by_name.get(&link.from_name),
                index_by_name.get(&link.to_name),
            ) {
                dsu.union(*left, *right);
            }
        }

        let mut by_root: BTreeMap<usize, Vec<String>> = BTreeMap::new();
        for (index, name) in names.iter().enumerate() {
            by_root
                .entry(dsu.find(index))
                .or_default()
                .push(name.clone());
        }

        let mut components = by_root
            .into_values()
            .map(|mut node_ids| {
                node_ids.sort_by(|left, right| natural_cmp(left, right));
                let node_set = node_ids.iter().cloned().collect::<HashSet<_>>();
                let link_count = self
                    .links
                    .iter()
                    .filter(|link| {
                        node_set.contains(&link.from_name) && node_set.contains(&link.to_name)
                    })
                    .count();
                let total_bp = node_ids
                    .iter()
                    .filter_map(|name| self.segments.get(name))
                    .map(|segment| segment.sequence.len())
                    .sum::<usize>();
                Component {
                    id: String::new(),
                    node_ids,
                    link_count,
                    total_bp,
                }
            })
            .collect::<Vec<_>>();

        components.sort_by(|left, right| {
            let left_element_count = left.node_ids.len() + left.link_count;
            let right_element_count = right.node_ids.len() + right.link_count;
            right_element_count
                .cmp(&left_element_count)
                .then_with(|| right.total_bp.cmp(&left.total_bp))
                .then_with(|| natural_cmp(&left.node_ids[0], &right.node_ids[0]))
        });
        for (index, component) in components.iter_mut().enumerate() {
            component.id = format!("subgraph_{:03}", index + 1);
        }
        let kind = if components.len() <= 1 {
            "single-subgraph"
        } else {
            "multi-subgraph"
        }
        .to_string();
        ComponentReport { kind, components }
    }
}

fn parse_orient(value: &str, path: &Path, line: usize) -> Result<char, OrgraftError> {
    match value {
        "+" => Ok('+'),
        "-" => Ok('-'),
        _ => Err(OrgraftError::InvalidArgument(format!(
            "{}:{line} invalid GFA orientation `{value}`",
            path.display()
        ))),
    }
}

#[derive(Debug, Clone)]
struct ComponentReport {
    kind: String,
    components: Vec<Component>,
}

#[derive(Debug, Clone)]
struct Component {
    id: String,
    node_ids: Vec<String>,
    link_count: usize,
    total_bp: usize,
}

#[derive(Debug, Clone)]
struct TopologyAudit {
    report: TopologyReport,
    components: ComponentReport,
}

fn build_topology_audit(input_gfa: &Path, graph: &GfaGraph) -> Result<TopologyAudit, OrgraftError> {
    let report = analyze_gfa(BufReader::new(File::open(input_gfa)?))?;
    let components = graph.component_report();
    Ok(TopologyAudit { report, components })
}

#[derive(Debug, Clone)]
struct ReferenceState {
    unrotated_fasta: Option<PathBuf>,
    pre_rotated_input_fasta: Option<PathBuf>,
    rotated_fasta: PathBuf,
    records: Vec<ReferenceRecordState>,
}

#[derive(Debug, Clone)]
struct ReferenceRecordState {
    alias: String,
    unrotated_id: String,
    rotated_id: String,
    reference_length: usize,
    orientation: char,
    rotation: usize,
    declared_rotation: Option<usize>,
    status: String,
    rollback_note: String,
}

fn prepare_reference(
    reference: &ReferenceInput,
    rotated_path: &Path,
    blastn: &Path,
) -> Result<ReferenceState, OrgraftError> {
    if let Some(parent) = rotated_path.parent() {
        fs::create_dir_all(parent)?;
    }
    let mut rotated_out = File::create(rotated_path)?;
    let mut states = Vec::new();

    match reference {
        ReferenceInput::AutoRotate(reference) => {
            let unrotated_records = read_fasta(reference)?;
            for (index, unrotated_record) in unrotated_records.iter().enumerate() {
                let alias = format!("reference_{:03}", index + 1);
                let single_reference = temp_file_path("orgraft-resolve-reference", "fasta");
                write_single_fasta(
                    &single_reference,
                    &unrotated_record.id,
                    &unrotated_record.sequence,
                )?;
                let selected_rotation =
                    choose_reference_rotation(unrotated_record, &single_reference, blastn)?;
                let _ = fs::remove_file(&single_reference);
                let rotated_record = FastaRecord {
                    id: format!("{} [rotation={}]", unrotated_record.id, selected_rotation),
                    sequence: rotate_sequence(
                        &unrotated_record.sequence,
                        selected_rotation as isize,
                    ),
                };
                let (orientation, detected_rotation, detected_status) =
                    detect_reference_rotation(&unrotated_record.sequence, &rotated_record.sequence)
                        .unwrap_or(('+', 0, "reference_rotation_unverified".to_string()));
                let status = if detected_status == "reference_rotation_unverified" {
                    detected_status
                } else {
                    "rotated_by_blastn_non_repeat_region".to_string()
                };
                write_fasta_record(&mut rotated_out, &alias, &rotated_record.sequence)?;
                states.push(ReferenceRecordState {
                    alias,
                    unrotated_id: unrotated_record.id.clone(),
                    rotated_id: rotated_record.id,
                    reference_length: rotated_record.sequence.len(),
                    orientation,
                    rotation: detected_rotation,
                    declared_rotation: Some(selected_rotation),
                    status,
                    rollback_note: rollback_note(
                        orientation,
                        detected_rotation,
                        unrotated_record.sequence.len(),
                    ),
                });
            }
            let state = ReferenceState {
                unrotated_fasta: Some(reference.to_path_buf()),
                pre_rotated_input_fasta: None,
                rotated_fasta: rotated_path.to_path_buf(),
                records: states,
            };
            Ok(state)
        }
        ReferenceInput::PreRotated(pre_rotated_reference) => {
            let rotated_records = read_fasta(pre_rotated_reference)?;
            for (index, rotated_record) in rotated_records.iter().enumerate() {
                let alias = format!("reference_{:03}", index + 1);
                write_fasta_record(&mut rotated_out, &alias, &rotated_record.sequence)?;
                states.push(ReferenceRecordState {
                    alias,
                    unrotated_id: ".".to_string(),
                    rotated_id: rotated_record.id.clone(),
                    reference_length: rotated_record.sequence.len(),
                    orientation: '+',
                    rotation: 0,
                    declared_rotation: parse_declared_rotation(&rotated_record.id),
                    status: "provided_pre_rotated_reference_unverified".to_string(),
                    rollback_note:
                        "pre-rotated reference was provided; original linear origin is not available"
                            .to_string(),
                });
            }
            let state = ReferenceState {
                unrotated_fasta: None,
                pre_rotated_input_fasta: Some(pre_rotated_reference.to_path_buf()),
                rotated_fasta: rotated_path.to_path_buf(),
                records: states,
            };
            Ok(state)
        }
    }
}

fn detect_reference_rotation(original: &str, candidate: &str) -> Option<(char, usize, String)> {
    if original.len() != candidate.len() {
        return None;
    }
    let original_upper = original.to_ascii_uppercase();
    let candidate_upper = candidate.to_ascii_uppercase();
    let doubled = format!("{original_upper}{original_upper}");
    if let Some(index) = doubled.find(&candidate_upper) {
        if index < original_upper.len() {
            return Some(('+', index, "reference_rotation_verified_plus".to_string()));
        }
    }

    let rc = reverse_complement(&original_upper);
    let doubled_rc = format!("{rc}{rc}");
    doubled_rc.find(&candidate_upper).and_then(|index| {
        (index < rc.len()).then_some((
            '-',
            index,
            "reference_rotation_verified_reverse_complement".to_string(),
        ))
    })
}

fn parse_declared_rotation(header: &str) -> Option<usize> {
    let marker = "rotation=";
    let start = header.find(marker)? + marker.len();
    let digits = header[start..]
        .chars()
        .take_while(|ch| ch.is_ascii_digit())
        .collect::<String>();
    digits.parse::<usize>().ok()
}

fn choose_reference_rotation(
    record: &FastaRecord,
    reference_fasta: &Path,
    blastn: &Path,
) -> Result<usize, OrgraftError> {
    let self_blast = temp_file_path("orgraft-resolve-reference-self", "tsv");
    run_blastn(
        blastn,
        reference_fasta,
        reference_fasta,
        &self_blast,
        BlastMode::Standard,
    )?;
    let lines = read_nonempty_lines(&self_blast)?;
    let _ = fs::remove_file(&self_blast);
    let repeat_mask = build_repeat_mask_from_blastn(&lines, record.sequence.len());
    let non_repeat_intervals = intervals_from_mask(&repeat_mask.mask, false);
    let circular_blocks = circular_non_repeat_blocks(&non_repeat_intervals, record.sequence.len());

    let Some(best) = circular_blocks.first() else {
        return Err(OrgraftError::InvalidArgument(format!(
            "no non-repeat region found while rotating reference"
        )));
    };
    if best.length <= MIN_NON_REPEAT_LEN {
        return Ok(0);
    }
    Ok(midpoint_on_circle(best, record.sequence.len()))
}

#[derive(Debug)]
struct RepeatMask {
    mask: Vec<bool>,
}

fn build_repeat_mask_from_blastn(lines: &[String], ref_len: usize) -> RepeatMask {
    let mut diff = vec![0i32; ref_len + 1];
    for line in lines {
        let fields = line.split('\t').collect::<Vec<_>>();
        if fields.len() < 12 {
            continue;
        }
        if is_full_self_hit(&fields, ref_len) {
            continue;
        }
        let (q_start, q_end) = sorted_pair_usize(fields[6], fields[7]);
        let start = q_start.saturating_sub(1).min(ref_len);
        let end = q_end.min(ref_len);
        if start >= end {
            continue;
        }
        diff[start] += 1;
        diff[end] -= 1;
    }

    let mut current = 0i32;
    let mask = diff
        .iter()
        .take(ref_len)
        .map(|value| {
            current += *value;
            current > 0
        })
        .collect::<Vec<_>>();
    RepeatMask { mask }
}

fn is_full_self_hit(fields: &[&str], ref_len: usize) -> bool {
    let (q_start, q_end) = sorted_pair_usize(fields[6], fields[7]);
    let (s_start, s_end) = sorted_pair_usize(fields[8], fields[9]);
    let aln_len = fields[3].parse::<usize>().unwrap_or(0);
    let pident = fields[2].parse::<f64>().unwrap_or(0.0);
    q_start == 1
        && q_end == ref_len
        && s_start == 1
        && s_end == ref_len
        && aln_len >= ref_len
        && pident >= 99.999
}

fn sorted_pair_usize(left: &str, right: &str) -> (usize, usize) {
    let left = left.parse::<usize>().unwrap_or(0);
    let right = right.parse::<usize>().unwrap_or(0);
    if left <= right {
        (left, right)
    } else {
        (right, left)
    }
}

#[derive(Debug, Clone)]
struct Interval {
    start: usize,
    end: usize,
    length: usize,
}

fn intervals_from_mask(mask: &[bool], want_repeat: bool) -> Vec<Interval> {
    let mut intervals = Vec::new();
    let mut start = None;
    for (index, value) in mask.iter().enumerate() {
        let selected = if want_repeat { *value } else { !*value };
        match (selected, start) {
            (true, None) => start = Some(index),
            (false, Some(start_index)) => {
                intervals.push(Interval {
                    start: start_index,
                    end: index - 1,
                    length: index - start_index,
                });
                start = None;
            }
            _ => {}
        }
    }
    if let Some(start_index) = start {
        intervals.push(Interval {
            start: start_index,
            end: mask.len() - 1,
            length: mask.len() - start_index,
        });
    }
    intervals
}

fn circular_non_repeat_blocks(non_repeat: &[Interval], ref_len: usize) -> Vec<Interval> {
    if non_repeat.is_empty() {
        return Vec::new();
    }
    let mut blocks = non_repeat.to_vec();
    if blocks.len() > 1
        && blocks.first().is_some_and(|item| item.start == 0)
        && blocks.last().is_some_and(|item| item.end == ref_len - 1)
    {
        let first = blocks.remove(0);
        let last = blocks.pop().expect("last interval exists");
        blocks.push(Interval {
            start: last.start,
            end: first.end,
            length: last.length + first.length,
        });
    }
    blocks.sort_by(|left, right| right.length.cmp(&left.length));
    blocks
}

fn midpoint_on_circle(block: &Interval, ref_len: usize) -> usize {
    if block.start <= block.end {
        (block.start + block.end) / 2
    } else {
        ((block.start + block.end + ref_len) / 2) % ref_len
    }
}

fn rollback_note(orientation: char, rotation: usize, length: usize) -> String {
    match orientation {
        '+' => format!("rotate aligned output right by {rotation} bp on length {length} to recover the original linear origin"),
        '-' => format!("rotate right by {rotation} bp on reverse-complement coordinates, then reverse-complement to recover the original linear origin"),
        _ => "unknown orientation; original reference copy is retained".to_string(),
    }
}

#[derive(Debug, Clone)]
struct RustPrepare {
    input_segments: usize,
    input_links: usize,
    merged_segments: usize,
    merged_links: usize,
    protected_nodes: Vec<String>,
    merge_mode: String,
    merge_action: String,
    merged_raw_gfa: PathBuf,
}

fn prepare_graph_in_rust(
    raw_graph: &GfaGraph,
    merged_raw_gfa: &Path,
) -> Result<RustPrepare, OrgraftError> {
    let (merged, merge_mode, protected_nodes) = merge_unambiguous_gfa(raw_graph);
    merged.write(&merged_raw_gfa)?;

    let merge_action = if merged.segments.len() < raw_graph.segments.len() {
        "auto_checked_and_merged"
    } else {
        "auto_checked_no_merge_needed"
    }
    .to_string();

    let prepare = RustPrepare {
        input_segments: raw_graph.segments.len(),
        input_links: raw_graph.links.len(),
        merged_segments: merged.segments.len(),
        merged_links: merged.links.len(),
        protected_nodes,
        merge_mode,
        merge_action,
        merged_raw_gfa: merged_raw_gfa.to_path_buf(),
    };
    Ok(prepare)
}

fn write_component_gfas(graph: &GfaGraph, graph_dir: &Path) -> Result<Vec<PathBuf>, OrgraftError> {
    fs::create_dir_all(graph_dir)?;
    let mut paths = Vec::new();
    for component in graph.component_report().components {
        let path = graph_dir.join(subgraph_gfa_filename(&component.id));
        graph.subgraph(&component.node_ids).write(&path)?;
        paths.push(path);
    }
    Ok(paths)
}

fn subgraph_gfa_filename(subgraph_id: &str) -> String {
    let suffix = subgraph_id.strip_prefix("subgraph_").unwrap_or(subgraph_id);
    format!("merged_unresolved_subgraph_{suffix}.gfa")
}

fn merge_unambiguous_gfa(raw: &GfaGraph) -> (GfaGraph, String, Vec<String>) {
    let protected_nodes = repeat_like_nodes(raw);
    let protected_set = protected_nodes.iter().cloned().collect::<HashSet<_>>();
    let mut adjacency: HashMap<String, Vec<String>> = HashMap::new();
    for name in raw.segments.keys() {
        adjacency.entry(name.clone()).or_default();
    }
    for link in &raw.links {
        if protected_set.contains(&link.from_name) || protected_set.contains(&link.to_name) {
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
    let mut names = raw.segments.keys().cloned().collect::<Vec<_>>();
    names.sort_by(|left, right| natural_cmp(left, right));
    for name in &names {
        if visited.contains(name) || protected_set.contains(name) {
            continue;
        }
        let mut stack = vec![name.clone()];
        let mut component = Vec::new();
        visited.insert(name.clone());
        while let Some(current) = stack.pop() {
            component.push(current.clone());
            let mut neighbors = adjacency.get(&current).cloned().unwrap_or_default();
            neighbors.sort_by(|left, right| natural_cmp(left, right));
            for next in neighbors {
                if visited.insert(next.clone()) {
                    stack.push(next);
                }
            }
        }
        components.push(component);
    }

    let links_by_pair = links_by_pair(raw);
    let mut path_records: Vec<(Vec<String>, Vec<char>)> = Vec::new();
    let mut used_nodes = HashSet::new();
    for component in components {
        if let Some(order) = ordered_linear_component(&component, &adjacency) {
            let orientations = orient_path(&order, &links_by_pair);
            used_nodes.extend(order.iter().cloned());
            path_records.push((order, orientations));
        } else {
            let mut component = component;
            component.sort_by(|left, right| natural_cmp(left, right));
            for name in component {
                used_nodes.insert(name.clone());
                path_records.push((vec![name], vec!['+']));
            }
        }
    }

    for name in &names {
        if used_nodes.insert(name.clone()) {
            path_records.push((vec![name.clone()], vec!['+']));
        }
    }

    let mut old_to_new = HashMap::new();
    let mut old_to_info = HashMap::new();
    let mut merged_segments = BTreeMap::new();
    let mut merged_segment_order = Vec::new();
    let mut mergeable_path_count = 0usize;
    for (order, orientations) in &path_records {
        if order.len() > 1 {
            mergeable_path_count += 1;
        }
        let name = order.join("_");
        let mut sequence = String::new();
        for (index, (old_name, orient)) in order.iter().zip(orientations).enumerate() {
            let mut part = raw
                .segments
                .get(old_name)
                .map(|segment| segment.sequence.clone())
                .unwrap_or_default();
            if *orient == '-' {
                part = reverse_complement(&part);
            }
            sequence.push_str(&part);
            old_to_new.insert(old_name.clone(), name.clone());
            old_to_info.insert(
                old_name.clone(),
                PathInfo {
                    path: order.clone(),
                    orientations: orientations.clone(),
                    index,
                },
            );
        }
        let tags = merged_segment_tags(raw, order, sequence.len());
        merged_segment_order.push(name.clone());
        merged_segments.insert(
            name.clone(),
            Segment {
                name,
                sequence,
                tags,
            },
        );
    }

    let mut seen_links = HashSet::new();
    let mut merged_links = Vec::new();
    for link in &raw.links {
        let Some(converted) = convert_link(link, &old_to_new, &old_to_info) else {
            continue;
        };
        let key = format!(
            "{}\t{}\t{}\t{}\t{}",
            converted.from_name,
            converted.from_orient,
            converted.to_name,
            converted.to_orient,
            converted.overlap
        );
        if seen_links.insert(key) {
            merged_links.push(converted);
        }
    }

    let merge_mode = if mergeable_path_count > 0 {
        "auto_non_repeat_linear_compaction"
    } else {
        "input_already_merged_or_no_linear_compaction"
    }
    .to_string();
    (
        GfaGraph {
            headers: raw.headers.clone(),
            segment_order: merged_segment_order,
            segments: merged_segments,
            links: merged_links,
            other_lines: raw.other_lines.clone(),
        },
        merge_mode,
        protected_nodes,
    )
}

fn merged_segment_tags(raw: &GfaGraph, order: &[String], length: usize) -> Vec<String> {
    let mut tags = Vec::new();
    tags.push(format!("LN:i:{length}"));
    tags.extend(coverage_tags_for_path(raw, order, length));
    tags.push(if order.len() > 1 {
        "SC:Z:linear_compaction".to_string()
    } else {
        "SC:Z:preserved_node".to_string()
    });
    tags.push("RR:Z:unresolved".to_string());
    tags
}

fn coverage_tags_for_path(raw: &GfaGraph, order: &[String], total_len: usize) -> Vec<String> {
    let mut depth_bases = 0.0f64;
    let mut depth_len = 0usize;
    let mut ab_sum = 0i64;
    let mut ac_sum = 0i64;
    let mut has_ab = false;
    let mut has_ac = false;

    for name in order {
        let Some(segment) = raw.segments.get(name) else {
            continue;
        };
        let segment_len = segment.sequence.len();
        if let Some(depth) = segment_depth(segment) {
            depth_bases += depth * segment_len as f64;
            depth_len += segment_len;
        }
        if let Some(value) = integer_tag(&segment.tags, "AB") {
            ab_sum += value;
            has_ab = true;
        }
        if let Some(value) = integer_tag(&segment.tags, "AC") {
            ac_sum += value;
            has_ac = true;
        }
    }

    let mut tags = Vec::new();
    if depth_len > 0 {
        tags.push(format!("DP:f:{:.6}", depth_bases / depth_len as f64));
    }
    tags.push("CM:Z:raw_node_DP_length_weighted".to_string());
    tags.push(format!("RL:i:{total_len}"));
    if depth_len > 0 {
        tags.push(format!("DB:f:{depth_bases:.3}"));
    }
    if has_ab {
        tags.push(format!("AB:i:{ab_sum}"));
    }
    if has_ac {
        tags.push(format!("AC:i:{ac_sum}"));
    }
    tags
}

fn segment_depth(segment: &Segment) -> Option<f64> {
    numeric_tag(
        &segment.tags,
        &["DP", "dp", "rd", "RD", "cov", "COV", "KC", "RC"],
    )
    .or_else(|| {
        let aligned_bases = integer_tag(&segment.tags, "AB")?;
        (!segment.sequence.is_empty())
            .then_some(aligned_bases as f64 / segment.sequence.len() as f64)
    })
}

fn numeric_tag(tags: &[String], names: &[&str]) -> Option<f64> {
    for name in names {
        for tag in tags {
            let Some((tag_name, value)) = parse_typed_tag(tag) else {
                continue;
            };
            if tag_name == *name {
                if let Ok(parsed) = value.parse::<f64>() {
                    return Some(parsed);
                }
            }
        }
    }
    None
}

fn integer_tag(tags: &[String], name: &str) -> Option<i64> {
    for tag in tags {
        let Some((tag_name, value)) = parse_typed_tag(tag) else {
            continue;
        };
        if tag_name == name {
            if let Ok(parsed) = value.parse::<i64>() {
                return Some(parsed);
            }
        }
    }
    None
}

fn parse_typed_tag(tag: &str) -> Option<(&str, &str)> {
    let mut parts = tag.splitn(3, ':');
    let tag_name = parts.next()?;
    let _tag_type = parts.next()?;
    let value = parts.next()?;
    Some((tag_name, value))
}

#[derive(Debug, Clone)]
struct PathInfo {
    path: Vec<String>,
    orientations: Vec<char>,
    index: usize,
}

fn repeat_like_nodes(raw: &GfaGraph) -> Vec<String> {
    let mut degrees: HashMap<String, usize> = raw
        .segments
        .keys()
        .map(|name| (name.clone(), 0usize))
        .collect();
    for link in &raw.links {
        *degrees.entry(link.from_name.clone()).or_default() += 1;
        *degrees.entry(link.to_name.clone()).or_default() += 1;
    }
    let branch_degrees = raw
        .segments
        .keys()
        .filter_map(|name| degrees.get(name).copied())
        .filter(|degree| *degree >= 3)
        .collect::<Vec<_>>();
    if branch_degrees.is_empty() {
        return Vec::new();
    }
    let branch_floor = *branch_degrees.iter().min().unwrap_or(&3);
    let has_higher_order_branch = raw
        .segments
        .keys()
        .any(|name| degrees.get(name).copied().unwrap_or(0) > branch_floor);
    let mut nodes = raw
        .segments
        .keys()
        .filter(|name| {
            let degree = degrees.get(*name).copied().unwrap_or(0);
            if has_higher_order_branch {
                degree > branch_floor
            } else {
                degree >= 3
            }
        })
        .cloned()
        .collect::<Vec<_>>();
    nodes.sort_by(|left, right| natural_cmp(left, right));
    nodes
}

fn ordered_linear_component(
    component: &[String],
    adjacency: &HashMap<String, Vec<String>>,
) -> Option<Vec<String>> {
    let component_set = component.iter().cloned().collect::<HashSet<_>>();
    let mut sub_degrees = HashMap::new();
    for name in component {
        let degree = adjacency
            .get(name)
            .into_iter()
            .flatten()
            .filter(|neighbor| component_set.contains(*neighbor))
            .count();
        if degree > 2 {
            return None;
        }
        sub_degrees.insert(name.clone(), degree);
    }
    let mut endpoints = sub_degrees
        .iter()
        .filter(|(_, degree)| **degree <= 1)
        .map(|(name, _)| name.clone())
        .collect::<Vec<_>>();
    endpoints.sort_by(|left, right| natural_cmp(left, right));
    if component.len() > 1 && endpoints.len() != 2 {
        return None;
    }
    let start = endpoints.first().cloned().or_else(|| {
        component
            .iter()
            .min_by(|left, right| natural_cmp(left, right))
            .cloned()
    })?;
    let mut order = Vec::new();
    let mut previous: Option<String> = None;
    let mut current = Some(start);
    while let Some(name) = current {
        order.push(name.clone());
        let mut candidates = adjacency
            .get(&name)
            .into_iter()
            .flatten()
            .filter(|neighbor| component_set.contains(*neighbor))
            .filter(|neighbor| previous.as_ref() != Some(*neighbor))
            .cloned()
            .collect::<Vec<_>>();
        candidates.sort_by(|left, right| natural_cmp(left, right));
        previous = Some(name);
        current = candidates.first().cloned();
    }
    (order.len() == component.len()).then_some(order)
}

fn links_by_pair(raw: &GfaGraph) -> HashMap<BTreeSet<String>, Vec<Link>> {
    let mut by_pair: HashMap<BTreeSet<String>, Vec<Link>> = HashMap::new();
    for link in &raw.links {
        let pair = BTreeSet::from([link.from_name.clone(), link.to_name.clone()]);
        by_pair.entry(pair).or_default().push(link.clone());
    }
    by_pair
}

fn orient_path(
    order: &[String],
    links_by_pair: &HashMap<BTreeSet<String>, Vec<Link>>,
) -> Vec<char> {
    for start_orient in ['+', '-'] {
        let mut orientations = vec![start_orient];
        let mut ok = true;
        for window in order.windows(2) {
            let Some(orient) = next_orient(
                &window[0],
                *orientations.last().unwrap(),
                &window[1],
                links_by_pair,
            ) else {
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

fn next_orient(
    current: &str,
    current_orient: char,
    next: &str,
    links_by_pair: &HashMap<BTreeSet<String>, Vec<Link>>,
) -> Option<char> {
    let pair = BTreeSet::from([current.to_string(), next.to_string()]);
    for link in links_by_pair.get(&pair)? {
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
    old_to_info: &HashMap<String, PathInfo>,
) -> Option<Link> {
    let merged_from = old_to_new.get(&link.from_name)?;
    let merged_to = old_to_new.get(&link.to_name)?;
    if merged_from == merged_to {
        return None;
    }
    let from_endpoint = old_from_endpoint(link.from_orient);
    let to_endpoint = old_to_endpoint(link.to_orient);
    let merged_from_endpoint =
        map_old_endpoint_to_merged(&link.from_name, from_endpoint, old_to_info)?;
    let merged_to_endpoint = map_old_endpoint_to_merged(&link.to_name, to_endpoint, old_to_info)?;
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
    old_to_info: &HashMap<String, PathInfo>,
) -> Option<SegmentEndpoint> {
    let info = old_to_info.get(old_name)?;
    let path_orient = *info.orientations.get(info.index)?;
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

fn flip_orient(orient: char) -> char {
    match orient {
        '+' => '-',
        '-' => '+',
        other => other,
    }
}

#[derive(Debug, Clone)]
struct RepeatResolutionResult {
    component_kind: String,
    subgraphs: Vec<SubgraphResolution>,
}

#[derive(Debug, Clone)]
struct SubgraphResolution {
    subgraph_id: String,
    reference_alias: String,
    unresolved_node_count: usize,
    unresolved_total_bp: usize,
    node_count: usize,
    total_bp: usize,
    resolution_engine: String,
    ready_repeat_nodes: Vec<String>,
    candidate_count: usize,
    selected_candidate: String,
    selected_circular: bool,
    selected_order: Vec<String>,
    score_method: String,
    score_value: Option<f64>,
    score_orientation: Option<char>,
    length_delta: Option<usize>,
    continuous_bp: Option<usize>,
    ordered_nodes: Vec<String>,
    oriented_nodes: Vec<String>,
    missing_reference_hits: Vec<String>,
}

#[derive(Debug, Clone)]
struct NodePlacement {
    node_id: String,
    reference_alias: Option<String>,
    reference_start: Option<usize>,
    orientation: char,
    sequence: String,
}

#[derive(Debug, Clone)]
struct AutoRepeatCandidate {
    id: String,
    graph: GfaGraph,
    order: Vec<RepeatStep>,
    signature: String,
    circular: bool,
    merged_order_count: usize,
}

#[derive(Debug, Clone)]
struct RepeatStep {
    node_id: String,
    duplicate_id: String,
    strategy: char,
}

#[derive(Debug, Clone)]
struct MergedSequence {
    sequence: String,
    ordered_nodes: Vec<String>,
    oriented_nodes: Vec<String>,
}

#[derive(Debug, Clone)]
struct CandidateScore {
    candidate_id: String,
    candidate_index: usize,
    reference_record: String,
    score: f64,
    method: String,
    orientation: char,
    length_delta: usize,
    continuous_bp: usize,
    continuous_fraction: f64,
    diagonal_fraction: f64,
}

#[derive(Debug, Clone)]
struct SequenceScore {
    score: f64,
    method: String,
    orientation: char,
    continuous_bp: usize,
    continuous_fraction: f64,
    diagonal_fraction: f64,
}

#[derive(Debug, Clone)]
struct ReferenceKmerIndex {
    record_id: String,
    sequence: String,
    indexes: Vec<KmerIndex>,
}

#[derive(Debug, Clone)]
struct KmerIndex {
    kmer: usize,
    index: HashMap<String, Vec<usize>>,
}

#[derive(Debug, Clone)]
struct CliRepeatSummary {
    candidate_count: usize,
    selected_candidate: String,
    selected_order: Vec<String>,
    reference_alias: String,
    score_method: String,
    score_value: Option<f64>,
    score_orientation: Option<char>,
    length_delta: Option<usize>,
    continuous_bp: Option<usize>,
}

fn auto_repeat_ready_node_ids(graph: &GfaGraph) -> Vec<String> {
    let mut ready = graph
        .segments
        .keys()
        .filter(|node_id| {
            let (counts, has_self_loop) = node_side_counts(graph, node_id);
            !has_self_loop
                && counts.get(&'-').copied().unwrap_or(0) == 2
                && counts.get(&'+').copied().unwrap_or(0) == 2
        })
        .cloned()
        .collect::<Vec<_>>();
    ready.sort();
    ready
}

fn node_side_counts(graph: &GfaGraph, node_id: &str) -> (HashMap<char, usize>, bool) {
    let mut counts = HashMap::from([('-', 0usize), ('+', 0usize)]);
    let mut has_self_loop = false;
    for link in unique_valid_links(graph) {
        if link.from_name == node_id && link.to_name == node_id {
            has_self_loop = true;
            continue;
        }
        if link.from_name == node_id {
            *counts
                .entry(endpoint_side(link.from_orient, "source"))
                .or_default() += 1;
        }
        if link.to_name == node_id {
            *counts
                .entry(endpoint_side(link.to_orient, "target"))
                .or_default() += 1;
        }
    }
    (counts, has_self_loop)
}

fn endpoint_side(orient: char, role: &str) -> char {
    match role {
        "target" => {
            if orient == '-' {
                '+'
            } else {
                '-'
            }
        }
        _ => {
            if orient == '-' {
                '-'
            } else {
                '+'
            }
        }
    }
}

fn link_endpoint_side(link: &Link, node_id: &str) -> Option<char> {
    let is_source = link.from_name == node_id;
    let is_target = link.to_name == node_id;
    if is_source && is_target {
        return None;
    }
    if is_source {
        return Some(endpoint_side(link.from_orient, "source"));
    }
    if is_target {
        return Some(endpoint_side(link.to_orient, "target"));
    }
    None
}

fn unique_valid_links(graph: &GfaGraph) -> Vec<Link> {
    let mut seen = HashSet::new();
    let mut unique = Vec::new();
    for link in &graph.links {
        if !graph.segments.contains_key(&link.from_name)
            || !graph.segments.contains_key(&link.to_name)
        {
            continue;
        }
        let key = graph_link_key(link);
        if seen.insert(key) {
            unique.push(link.clone());
        }
    }
    unique
}

fn graph_link_key(link: &Link) -> String {
    format!(
        "{}|{}|{}",
        canonical_link_key(link),
        link.overlap,
        link.tags.join("|")
    )
}

fn canonical_link_key(link: &Link) -> String {
    let mut endpoints = [
        (
            link.from_name.clone(),
            endpoint_side(link.from_orient, "source"),
        ),
        (
            link.to_name.clone(),
            endpoint_side(link.to_orient, "target"),
        ),
    ];
    endpoints.sort();
    format!(
        "{}{}--{}{}",
        endpoints[0].0, endpoints[0].1, endpoints[1].0, endpoints[1].1
    )
}

fn graph_topology_signature(graph: &GfaGraph) -> String {
    let nodes = graph
        .segments
        .iter()
        .map(|(name, segment)| {
            format!(
                "{name}:{}:{}",
                segment.sequence.len(),
                segment.tags.join(",")
            )
        })
        .collect::<Vec<_>>()
        .join(";");
    let mut links = graph.links.iter().map(graph_link_key).collect::<Vec<_>>();
    links.sort();
    format!("nodes=[{nodes}] links=[{}]", links.join(";"))
}

fn graph_is_connected(graph: &GfaGraph) -> bool {
    graph.component_report().components.len() == 1
}

fn graph_is_circular_subgraph(graph: &GfaGraph) -> bool {
    if !graph_is_connected(graph) || unique_valid_links(graph).len() != graph.segments.len() {
        return false;
    }
    graph.segments.keys().all(|node_id| {
        let (counts, has_self_loop) = node_side_counts(graph, node_id);
        !has_self_loop
            && counts.get(&'-').copied().unwrap_or(0) == 1
            && counts.get(&'+').copied().unwrap_or(0) == 1
    })
}

fn build_auto_repeat_resolution_candidates(
    graph: &GfaGraph,
    max_states: usize,
    max_candidates: usize,
) -> Result<(Vec<AutoRepeatCandidate>, Option<String>), OrgraftError> {
    let mut search_graph = graph.clone();
    deduplicate_links(&mut search_graph);

    if !graph_is_connected(&search_graph) {
        return Err(OrgraftError::InvalidArgument(
            "auto repeat resolution requires a connected subgraph".to_string(),
        ));
    }
    let targets = auto_repeat_ready_node_ids(&search_graph);
    if targets.is_empty() {
        return Ok((
            Vec::new(),
            Some("No 2-in/2-out repeat nodes were found in the selected subgraph.".to_string()),
        ));
    }

    #[derive(Clone)]
    struct State {
        graph: GfaGraph,
        remaining: Vec<String>,
        order: Vec<RepeatStep>,
    }

    let mut states = vec![State {
        graph: search_graph,
        remaining: targets,
        order: Vec::new(),
    }];
    let mut visited = HashSet::new();
    let mut final_by_signature: HashMap<String, usize> = HashMap::new();
    let mut finals: Vec<AutoRepeatCandidate> = Vec::new();
    let mut explored_state_count = 0usize;
    let mut truncated = false;

    while let Some(state) = states.pop() {
        let mut remaining = state.remaining.clone();
        remaining.sort();
        let state_key = format!(
            "{}|{}",
            remaining.join(","),
            graph_topology_signature(&state.graph)
        );
        if !visited.insert(state_key) {
            continue;
        }
        explored_state_count += 1;
        if explored_state_count > max_states {
            truncated = true;
            break;
        }

        if remaining.is_empty() {
            let final_signature = graph_topology_signature(&state.graph);
            if let Some(index) = final_by_signature.get(&final_signature).copied() {
                finals[index].merged_order_count += 1;
                continue;
            }
            let candidate = AutoRepeatCandidate {
                id: String::new(),
                graph: state.graph.clone(),
                order: state.order.clone(),
                signature: final_signature.clone(),
                circular: graph_is_circular_subgraph(&state.graph),
                merged_order_count: 1,
            };
            final_by_signature.insert(final_signature, finals.len());
            finals.push(candidate);
            if finals.len() >= max_candidates {
                truncated = !states.is_empty();
                break;
            }
            continue;
        }

        for node_id in &remaining {
            let (counts, has_self_loop) = node_side_counts(&state.graph, node_id);
            if has_self_loop
                || counts.get(&'-').copied().unwrap_or(0) != 2
                || counts.get(&'+').copied().unwrap_or(0) != 2
            {
                continue;
            }
            for strategy in ['A', 'B'] {
                let mut next_graph = state.graph.clone();
                let Ok(duplicate_id) = duplicate_node(&mut next_graph, node_id) else {
                    continue;
                };
                if repeat_resolve_node(&mut next_graph, node_id, &duplicate_id, strategy).is_err() {
                    continue;
                }
                deduplicate_links(&mut next_graph);
                if !graph_is_connected(&next_graph) {
                    continue;
                }
                let next_remaining = remaining
                    .iter()
                    .filter(|candidate| *candidate != node_id)
                    .cloned()
                    .collect::<Vec<_>>();
                let mut next_order = state.order.clone();
                next_order.push(RepeatStep {
                    node_id: node_id.clone(),
                    duplicate_id,
                    strategy,
                });
                states.push(State {
                    graph: next_graph,
                    remaining: next_remaining,
                    order: next_order,
                });
            }
        }
    }

    let mut candidates = order_auto_repeat_candidates(finals);
    for (index, candidate) in candidates.iter_mut().enumerate() {
        candidate.id = format!("auto_repeat_{:03}", index + 1);
    }
    let warning = if truncated {
        Some(format!(
            "Search stopped after {explored_state_count} states. Showing {} unique candidate results.",
            candidates.len()
        ))
    } else if candidates.is_empty() {
        Some("No connected result could resolve all 2-in/2-out repeat nodes.".to_string())
    } else {
        None
    };
    Ok((candidates, warning))
}

fn duplicate_node(graph: &mut GfaGraph, node_id: &str) -> Result<String, OrgraftError> {
    let segment = graph
        .segments
        .get(node_id)
        .cloned()
        .ok_or_else(|| OrgraftError::InvalidArgument(format!("node not found: {node_id}")))?;
    let duplicate_id = next_duplicate_id(graph, node_id);
    let mut duplicate = segment;
    duplicate.name = duplicate_id.clone();
    graph.segments.insert(duplicate_id.clone(), duplicate);
    graph.segment_order.push(duplicate_id.clone());

    let incident_links = graph
        .links
        .iter()
        .filter(|link| link.from_name == node_id || link.to_name == node_id)
        .cloned()
        .collect::<Vec<_>>();
    let mut seen = graph
        .links
        .iter()
        .map(canonical_link_key)
        .collect::<HashSet<_>>();
    for mut link in incident_links {
        if link.from_name == node_id {
            link.from_name = duplicate_id.clone();
        }
        if link.to_name == node_id {
            link.to_name = duplicate_id.clone();
        }
        if seen.insert(canonical_link_key(&link)) {
            graph.links.push(link);
        }
    }
    Ok(duplicate_id)
}

fn next_duplicate_id(graph: &GfaGraph, source_id: &str) -> String {
    let mut index = 1usize;
    loop {
        let candidate = format!("{source_id}_copy{index}");
        if !graph.segments.contains_key(&candidate) {
            return candidate;
        }
        index += 1;
    }
}

fn repeat_resolve_node(
    graph: &mut GfaGraph,
    node_id: &str,
    duplicate_id: &str,
    strategy: char,
) -> Result<(), OrgraftError> {
    let source_groups = incident_link_indices_by_side(graph, node_id);
    let duplicate_groups = incident_link_indices_by_side(graph, duplicate_id);
    let source_minus = source_groups.get(&'-').cloned().unwrap_or_default();
    let source_plus = source_groups.get(&'+').cloned().unwrap_or_default();
    let duplicate_minus = duplicate_groups.get(&'-').cloned().unwrap_or_default();
    let duplicate_plus = duplicate_groups.get(&'+').cloned().unwrap_or_default();
    if source_minus.len() != 2
        || source_plus.len() != 2
        || duplicate_minus.len() != 2
        || duplicate_plus.len() != 2
    {
        return Err(OrgraftError::InvalidArgument(format!(
            "repeat resolution expects {node_id} and {duplicate_id} to each have two links on both ends"
        )));
    }

    let copied_by_source = match_duplicate_links(
        graph,
        &[source_minus.clone(), source_plus.clone()].concat(),
        &[duplicate_minus, duplicate_plus].concat(),
        node_id,
        duplicate_id,
    )?;
    let source_keep = match strategy {
        'A' => vec![source_minus[0], source_plus[0]],
        'B' => vec![source_minus[0], source_plus[1]],
        _ => {
            return Err(OrgraftError::InvalidArgument(
                "repeat resolution strategy must be A or B".to_string(),
            ))
        }
    };
    let duplicate_keep = match strategy {
        'A' => vec![
            *copied_by_source.get(&source_minus[1]).ok_or_else(|| {
                OrgraftError::InvalidArgument("missing duplicated minus link".to_string())
            })?,
            *copied_by_source.get(&source_plus[1]).ok_or_else(|| {
                OrgraftError::InvalidArgument("missing duplicated plus link".to_string())
            })?,
        ],
        'B' => vec![
            *copied_by_source.get(&source_minus[1]).ok_or_else(|| {
                OrgraftError::InvalidArgument("missing duplicated minus link".to_string())
            })?,
            *copied_by_source.get(&source_plus[0]).ok_or_else(|| {
                OrgraftError::InvalidArgument("missing duplicated plus link".to_string())
            })?,
        ],
        _ => unreachable!(),
    };

    let keep_ids = source_keep
        .into_iter()
        .chain(duplicate_keep)
        .collect::<HashSet<_>>();
    let candidate_ids = source_groups
        .values()
        .chain(duplicate_groups.values())
        .flatten()
        .copied()
        .collect::<HashSet<_>>();
    graph.links = graph
        .links
        .iter()
        .enumerate()
        .filter(|(index, _)| !candidate_ids.contains(index) || keep_ids.contains(index))
        .map(|(_, link)| link.clone())
        .collect();
    Ok(())
}

fn incident_link_indices_by_side(graph: &GfaGraph, node_id: &str) -> HashMap<char, Vec<usize>> {
    let mut groups = HashMap::from([('-', Vec::new()), ('+', Vec::new())]);
    for (index, link) in graph.links.iter().enumerate() {
        if let Some(side) = link_endpoint_side(link, node_id) {
            groups.entry(side).or_default().push(index);
        }
    }
    groups
}

fn match_duplicate_links(
    graph: &GfaGraph,
    source_indices: &[usize],
    duplicate_indices: &[usize],
    node_id: &str,
    duplicate_id: &str,
) -> Result<HashMap<usize, usize>, OrgraftError> {
    let mut duplicate_by_signature: HashMap<String, Vec<usize>> = HashMap::new();
    for index in duplicate_indices {
        let signature = repeat_link_signature(&graph.links[*index], duplicate_id);
        duplicate_by_signature
            .entry(signature)
            .or_default()
            .push(*index);
    }

    let mut copied_by_source = HashMap::new();
    for index in source_indices {
        let signature = repeat_link_signature(&graph.links[*index], node_id);
        let candidates = duplicate_by_signature.get_mut(&signature).ok_or_else(|| {
            OrgraftError::InvalidArgument(
                "repeat resolution could not match duplicated links".to_string(),
            )
        })?;
        if candidates.is_empty() {
            return Err(OrgraftError::InvalidArgument(
                "repeat resolution duplicate link signature was exhausted".to_string(),
            ));
        }
        copied_by_source.insert(*index, candidates.remove(0));
    }
    Ok(copied_by_source)
}

fn repeat_link_signature(link: &Link, node_id: &str) -> String {
    let source = if link.from_name == node_id {
        "__SELF__"
    } else {
        &link.from_name
    };
    let target = if link.to_name == node_id {
        "__SELF__"
    } else {
        &link.to_name
    };
    format!(
        "{}\t{}\t{}\t{}\t{}\t{}",
        source,
        link.from_orient,
        target,
        link.to_orient,
        link.overlap,
        link.tags.join("|")
    )
}

fn deduplicate_links(graph: &mut GfaGraph) -> usize {
    let mut seen = HashSet::new();
    let before = graph.links.len();
    graph.links = graph
        .links
        .iter()
        .filter(|link| seen.insert(canonical_link_key(link)))
        .cloned()
        .collect();
    before - graph.links.len()
}

fn order_auto_repeat_candidates(candidates: Vec<AutoRepeatCandidate>) -> Vec<AutoRepeatCandidate> {
    let mut indexed = candidates.into_iter().enumerate().collect::<Vec<_>>();
    indexed.sort_by(|(left_index, left), (right_index, right)| {
        auto_repeat_candidate_key(left, *left_index)
            .cmp(&auto_repeat_candidate_key(right, *right_index))
    });
    indexed
        .into_iter()
        .map(|(_, candidate)| candidate)
        .collect()
}

fn auto_repeat_candidate_key(
    candidate: &AutoRepeatCandidate,
    original_index: usize,
) -> (u8, usize, String, String, usize) {
    match gfa_editor_merge_all_sequence(&candidate.graph) {
        Ok(merged) => (
            0,
            merged.len(),
            head_to_tail_sequence_feature(&merged),
            candidate.signature.clone(),
            original_index,
        ),
        Err(_) => (
            1,
            candidate.graph.segments.len(),
            String::new(),
            candidate.signature.clone(),
            original_index,
        ),
    }
}

fn gfa_editor_merge_all_sequence(graph: &GfaGraph) -> Result<String, OrgraftError> {
    if graph.segments.is_empty() {
        return Err(OrgraftError::InvalidArgument(
            "cannot merge an empty graph".to_string(),
        ));
    }
    if graph.segments.len() == 1 {
        let segment = graph
            .ordered_segments()
            .first()
            .copied()
            .or_else(|| graph.segments.values().next())
            .unwrap();
        return Ok(segment.sequence.clone());
    }

    let mut merged_graph = graph.clone();
    let node_ids = merged_graph.ordered_segment_names();
    let (path_node_ids, retained_cycle_link) =
        gfa_editor_selected_merge_path(&merged_graph, &node_ids)?;
    if let Some(retained_index) = retained_cycle_link {
        if retained_index < merged_graph.links.len() {
            merged_graph.links.remove(retained_index);
        }
    }

    let mut current_node_id = path_node_ids
        .first()
        .cloned()
        .ok_or_else(|| OrgraftError::InvalidArgument("empty merge path".to_string()))?;
    for next_node_id in path_node_ids.iter().skip(1) {
        current_node_id = gfa_editor_merge_unique_link_between(
            &mut merged_graph,
            &current_node_id,
            next_node_id,
        )?;
    }
    merged_graph
        .segments
        .get(&current_node_id)
        .map(|segment| segment.sequence.clone())
        .ok_or_else(|| {
            OrgraftError::InvalidArgument(
                "GFA_Editor-style merge did not produce a segment".to_string(),
            )
        })
}

fn gfa_editor_selected_merge_path(
    graph: &GfaGraph,
    node_ids: &[String],
) -> Result<(Vec<String>, Option<usize>), OrgraftError> {
    let selected = node_ids.iter().cloned().collect::<HashSet<_>>();
    let internal_links = graph
        .links
        .iter()
        .enumerate()
        .filter(|(_, link)| {
            selected.contains(&link.from_name)
                && selected.contains(&link.to_name)
                && link.from_name != link.to_name
        })
        .map(|(index, _)| index)
        .collect::<Vec<_>>();
    let mut adjacency = node_ids
        .iter()
        .map(|node_id| (node_id.clone(), Vec::<usize>::new()))
        .collect::<HashMap<_, _>>();
    let mut pair_counts: HashMap<(String, String), usize> = HashMap::new();
    for link_index in &internal_links {
        let link = &graph.links[*link_index];
        let pair = sorted_pair(&link.from_name, &link.to_name);
        *pair_counts.entry(pair).or_default() += 1;
        adjacency
            .entry(link.from_name.clone())
            .or_default()
            .push(*link_index);
        adjacency
            .entry(link.to_name.clone())
            .or_default()
            .push(*link_index);
    }

    if selected.len() == 2 && internal_links.len() == 2 {
        return Ok((node_ids.to_vec(), Some(internal_links[1])));
    }
    if pair_counts.values().any(|count| *count != 1) {
        return Err(OrgraftError::InvalidArgument(
            "selected contigs must have exactly one link between each connected pair".to_string(),
        ));
    }
    if internal_links.len() == selected.len() {
        return gfa_editor_selected_merge_cycle_path(graph, node_ids, &adjacency);
    }
    if internal_links.len() != selected.len().saturating_sub(1) {
        return Err(OrgraftError::InvalidArgument(
            "selected contigs must form a single path or simple cycle".to_string(),
        ));
    }

    let endpoints = node_ids
        .iter()
        .filter(|node_id| adjacency.get(*node_id).map(Vec::len).unwrap_or(0) == 1)
        .cloned()
        .collect::<Vec<_>>();
    let middle_count = node_ids
        .iter()
        .filter(|node_id| adjacency.get(*node_id).map(Vec::len).unwrap_or(0) == 2)
        .count();
    let valid_degrees = if selected.len() == 2 {
        endpoints.len() == 2
    } else {
        endpoints.len() == 2 && middle_count == selected.len().saturating_sub(2)
    };
    if !valid_degrees {
        return Err(OrgraftError::InvalidArgument(
            "selected contigs must form one unbranched head-to-tail path".to_string(),
        ));
    }

    let start_node_id = endpoints
        .iter()
        .min_by_key(|node_id| {
            node_ids
                .iter()
                .position(|candidate| candidate == *node_id)
                .unwrap_or(usize::MAX)
        })
        .cloned()
        .unwrap();
    let mut path_node_ids = vec![start_node_id.clone()];
    let mut seen = HashSet::from([start_node_id.clone()]);
    let mut previous_node_id: Option<String> = None;
    let mut current_node_id = start_node_id;
    while path_node_ids.len() < selected.len() {
        let candidates = adjacency
            .get(&current_node_id)
            .cloned()
            .unwrap_or_default()
            .into_iter()
            .filter(|link_index| {
                previous_node_id.as_ref()
                    != Some(&other_link_node(
                        &graph.links[*link_index],
                        &current_node_id,
                    ))
            })
            .collect::<Vec<_>>();
        if candidates.len() != 1 {
            return Err(OrgraftError::InvalidArgument(
                "selected contigs must form one unbranched head-to-tail path".to_string(),
            ));
        }
        let next_node_id = other_link_node(&graph.links[candidates[0]], &current_node_id);
        if !seen.insert(next_node_id.clone()) {
            return Err(OrgraftError::InvalidArgument(
                "selected contigs must form one unbranched head-to-tail path".to_string(),
            ));
        }
        path_node_ids.push(next_node_id.clone());
        previous_node_id = Some(current_node_id);
        current_node_id = next_node_id;
    }
    Ok((path_node_ids, None))
}

fn gfa_editor_selected_merge_cycle_path(
    graph: &GfaGraph,
    node_ids: &[String],
    adjacency: &HashMap<String, Vec<usize>>,
) -> Result<(Vec<String>, Option<usize>), OrgraftError> {
    if node_ids.len() < 3
        || node_ids
            .iter()
            .any(|node_id| adjacency.get(node_id).map(Vec::len).unwrap_or(0) != 2)
    {
        return Err(OrgraftError::InvalidArgument(
            "selected contigs must form one simple cycle".to_string(),
        ));
    }
    let start_node_id = node_ids[0].clone();
    let start_links = adjacency.get(&start_node_id).cloned().unwrap_or_default();
    let first_link = start_links
        .iter()
        .copied()
        .find(|link_index| {
            node_ids.len() > 1
                && other_link_node(&graph.links[*link_index], &start_node_id) == node_ids[1]
        })
        .or_else(|| start_links.first().copied())
        .ok_or_else(|| {
            OrgraftError::InvalidArgument("selected contigs must form one simple cycle".to_string())
        })?;
    let mut path_node_ids = vec![start_node_id.clone()];
    let mut seen = HashSet::from([start_node_id.clone()]);
    let mut previous_link: Option<usize> = None;
    let mut current_node_id = start_node_id;

    while path_node_ids.len() < node_ids.len() {
        let mut candidates = adjacency
            .get(&current_node_id)
            .cloned()
            .unwrap_or_default()
            .into_iter()
            .filter(|link_index| Some(*link_index) != previous_link)
            .collect::<Vec<_>>();
        if current_node_id == node_ids[0] {
            candidates = vec![first_link];
        }
        if candidates.len() != 1 {
            return Err(OrgraftError::InvalidArgument(
                "selected contigs must form one simple cycle".to_string(),
            ));
        }
        let link_index = candidates[0];
        let next_node_id = other_link_node(&graph.links[link_index], &current_node_id);
        if !seen.insert(next_node_id.clone()) {
            return Err(OrgraftError::InvalidArgument(
                "selected contigs must form one simple cycle".to_string(),
            ));
        }
        path_node_ids.push(next_node_id.clone());
        previous_link = Some(link_index);
        current_node_id = next_node_id;
    }

    let retained_candidates = adjacency
        .get(&current_node_id)
        .cloned()
        .unwrap_or_default()
        .into_iter()
        .filter(|link_index| Some(*link_index) != previous_link)
        .filter(|link_index| {
            other_link_node(&graph.links[*link_index], &current_node_id) == node_ids[0]
        })
        .collect::<Vec<_>>();
    if retained_candidates.len() != 1 || seen != node_ids.iter().cloned().collect::<HashSet<_>>() {
        return Err(OrgraftError::InvalidArgument(
            "selected contigs must form one simple cycle".to_string(),
        ));
    }
    Ok((path_node_ids, Some(retained_candidates[0])))
}

fn gfa_editor_merge_unique_link_between(
    graph: &mut GfaGraph,
    first_node_id: &str,
    second_node_id: &str,
) -> Result<String, OrgraftError> {
    let matches = graph
        .links
        .iter()
        .enumerate()
        .filter(|(_, link)| {
            link.from_name != link.to_name
                && ((link.from_name == first_node_id && link.to_name == second_node_id)
                    || (link.from_name == second_node_id && link.to_name == first_node_id))
        })
        .map(|(index, _)| index)
        .collect::<Vec<_>>();
    if matches.len() != 1 {
        return Err(OrgraftError::InvalidArgument(
            "merge path no longer has exactly one link between adjacent contigs".to_string(),
        ));
    }
    gfa_editor_merge_link(graph, matches[0])
}

fn gfa_editor_merge_link(graph: &mut GfaGraph, link_index: usize) -> Result<String, OrgraftError> {
    let merge_link_record = graph
        .links
        .get(link_index)
        .cloned()
        .ok_or_else(|| OrgraftError::InvalidArgument("link not found".to_string()))?;
    let source_id = merge_link_record.from_name.clone();
    let target_id = merge_link_record.to_name.clone();
    if source_id == target_id {
        return Err(OrgraftError::InvalidArgument(
            "cannot merge a self-link".to_string(),
        ));
    }
    let source_segment = graph.segments.get(&source_id).cloned().ok_or_else(|| {
        OrgraftError::InvalidArgument(format!("cannot merge missing segment `{source_id}`"))
    })?;
    let target_segment = graph.segments.get(&target_id).cloned().ok_or_else(|| {
        OrgraftError::InvalidArgument(format!("cannot merge missing segment `{target_id}`"))
    })?;
    let overlap = overlap_length_from_cigar(&merge_link_record.overlap)
        .min(source_segment.sequence.len())
        .min(target_segment.sequence.len());
    let new_id = next_merged_id(graph, &source_id, &target_id);
    let source_sequence = if merge_link_record.from_orient == '-' {
        reverse_complement(&source_segment.sequence)
    } else {
        source_segment.sequence.clone()
    };
    let target_sequence = if merge_link_record.to_orient == '-' {
        reverse_complement(&target_segment.sequence)
    } else {
        target_segment.sequence.clone()
    };
    let mut sequence = source_sequence;
    sequence.push_str(&target_sequence[overlap.min(target_sequence.len())..]);

    let mut rewired_links = Vec::new();
    for (index, link) in graph.links.clone().into_iter().enumerate() {
        if index == link_index {
            continue;
        }
        let touches_source = link.from_name == source_id || link.to_name == source_id;
        let touches_target = link.from_name == target_id || link.to_name == target_id;
        if touches_source && touches_target {
            return Err(OrgraftError::InvalidArgument(
                "cannot merge nodes with links that would become self-links".to_string(),
            ));
        }
        let mut duplicate = link;
        if touches_source {
            replace_link_endpoint(&mut duplicate, &source_id, &new_id, '-');
        }
        if touches_target {
            replace_link_endpoint(&mut duplicate, &target_id, &new_id, '+');
        }
        rewired_links.push(duplicate);
    }

    let mut new_segment_order = Vec::new();
    let mut inserted = false;
    for segment_id in &graph.segment_order {
        if segment_id == &source_id {
            new_segment_order.push(new_id.clone());
            inserted = true;
            continue;
        }
        if segment_id == &target_id {
            if !inserted {
                new_segment_order.push(new_id.clone());
                inserted = true;
            }
            continue;
        }
        if graph.segments.contains_key(segment_id) {
            new_segment_order.push(segment_id.clone());
        }
    }
    if !inserted {
        new_segment_order.push(new_id.clone());
    }

    graph.segments.remove(&source_id);
    graph.segments.remove(&target_id);
    graph.segments.insert(
        new_id.clone(),
        Segment {
            name: new_id.clone(),
            sequence,
            tags: Vec::new(),
        },
    );
    graph.segment_order = new_segment_order;
    graph.links = rewired_links;
    Ok(new_id)
}

fn next_merged_id(graph: &GfaGraph, source_id: &str, target_id: &str) -> String {
    let base = format!("{source_id}_{target_id}");
    if !graph.segments.contains_key(&base) {
        return base;
    }
    let mut index = 1usize;
    loop {
        let candidate = format!("{base}_merge{index}");
        if !graph.segments.contains_key(&candidate) {
            return candidate;
        }
        index += 1;
    }
}

fn replace_link_endpoint(link: &mut Link, old_id: &str, new_id: &str, new_side: char) {
    if link.from_name == old_id {
        link.from_name = new_id.to_string();
        link.from_orient = orient_for_endpoint_side(new_side, "source");
    }
    if link.to_name == old_id {
        link.to_name = new_id.to_string();
        link.to_orient = orient_for_endpoint_side(new_side, "target");
    }
}

fn orient_for_endpoint_side(side: char, role: &str) -> char {
    match role {
        "target" => {
            if side == '+' {
                '-'
            } else {
                '+'
            }
        }
        _ => {
            if side == '+' {
                '+'
            } else {
                '-'
            }
        }
    }
}

fn sorted_pair(left: &str, right: &str) -> (String, String) {
    if left <= right {
        (left.to_string(), right.to_string())
    } else {
        (right.to_string(), left.to_string())
    }
}

fn other_link_node(link: &Link, node_id: &str) -> String {
    if link.from_name == node_id {
        return link.to_name.clone();
    }
    if link.to_name == node_id {
        return link.from_name.clone();
    }
    node_id.to_string()
}

fn head_to_tail_sequence_feature(sequence: &str) -> String {
    let sequence_length = sequence.len();
    let kmer_size = if sequence_length >= 31 {
        31
    } else {
        sequence_length.max(1)
    };
    let anchor_count = 256.min((sequence_length.saturating_sub(kmer_size) + 1).max(1));
    let mut anchors = Vec::new();
    if anchor_count == 1 {
        anchors.push(sequence.get(..kmer_size).unwrap_or(sequence).to_string());
    } else {
        let last_position = sequence_length - kmer_size;
        for index in 0..anchor_count {
            let position = (index * last_position) / (anchor_count - 1);
            anchors.push(sequence[position..position + kmer_size].to_string());
        }
    }
    let head = if sequence_length > 256 {
        &sequence[..256]
    } else {
        sequence
    };
    let tail = if sequence_length > 256 {
        &sequence[sequence_length - 256..]
    } else {
        sequence
    };
    format!(
        "{kmer_size}\u{1f}{}\u{1f}{head}\u{1f}{tail}\u{1f}{}",
        anchors.join("\u{1e}"),
        sha256_hex(sequence.as_bytes())
    )
}

fn sha256_hex(input: &[u8]) -> String {
    const K: [u32; 64] = [
        0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5, 0x3956c25b, 0x59f111f1, 0x923f82a4,
        0xab1c5ed5, 0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3, 0x72be5d74, 0x80deb1fe,
        0x9bdc06a7, 0xc19bf174, 0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc, 0x2de92c6f,
        0x4a7484aa, 0x5cb0a9dc, 0x76f988da, 0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7,
        0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967, 0x27b70a85, 0x2e1b2138, 0x4d2c6dfc,
        0x53380d13, 0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85, 0xa2bfe8a1, 0xa81a664b,
        0xc24b8b70, 0xc76c51a3, 0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070, 0x19a4c116,
        0x1e376c08, 0x2748774c, 0x34b0bcb5, 0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3,
        0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208, 0x90befffa, 0xa4506ceb, 0xbef9a3f7,
        0xc67178f2,
    ];
    let mut h = [
        0x6a09e667u32,
        0xbb67ae85,
        0x3c6ef372,
        0xa54ff53a,
        0x510e527f,
        0x9b05688c,
        0x1f83d9ab,
        0x5be0cd19,
    ];
    let bit_len = (input.len() as u64).wrapping_mul(8);
    let mut message = input.to_vec();
    message.push(0x80);
    while message.len() % 64 != 56 {
        message.push(0);
    }
    message.extend_from_slice(&bit_len.to_be_bytes());

    for chunk in message.chunks_exact(64) {
        let mut w = [0u32; 64];
        for (index, word) in w.iter_mut().take(16).enumerate() {
            let offset = index * 4;
            *word = u32::from_be_bytes([
                chunk[offset],
                chunk[offset + 1],
                chunk[offset + 2],
                chunk[offset + 3],
            ]);
        }
        for index in 16..64 {
            let s0 = w[index - 15].rotate_right(7)
                ^ w[index - 15].rotate_right(18)
                ^ (w[index - 15] >> 3);
            let s1 = w[index - 2].rotate_right(17)
                ^ w[index - 2].rotate_right(19)
                ^ (w[index - 2] >> 10);
            w[index] = w[index - 16]
                .wrapping_add(s0)
                .wrapping_add(w[index - 7])
                .wrapping_add(s1);
        }

        let mut a = h[0];
        let mut b = h[1];
        let mut c = h[2];
        let mut d = h[3];
        let mut e = h[4];
        let mut f = h[5];
        let mut g = h[6];
        let mut hh = h[7];

        for index in 0..64 {
            let s1 = e.rotate_right(6) ^ e.rotate_right(11) ^ e.rotate_right(25);
            let ch = (e & f) ^ ((!e) & g);
            let temp1 = hh
                .wrapping_add(s1)
                .wrapping_add(ch)
                .wrapping_add(K[index])
                .wrapping_add(w[index]);
            let s0 = a.rotate_right(2) ^ a.rotate_right(13) ^ a.rotate_right(22);
            let maj = (a & b) ^ (a & c) ^ (b & c);
            let temp2 = s0.wrapping_add(maj);

            hh = g;
            g = f;
            f = e;
            e = d.wrapping_add(temp1);
            d = c;
            c = b;
            b = a;
            a = temp1.wrapping_add(temp2);
        }

        h[0] = h[0].wrapping_add(a);
        h[1] = h[1].wrapping_add(b);
        h[2] = h[2].wrapping_add(c);
        h[3] = h[3].wrapping_add(d);
        h[4] = h[4].wrapping_add(e);
        h[5] = h[5].wrapping_add(f);
        h[6] = h[6].wrapping_add(g);
        h[7] = h[7].wrapping_add(hh);
    }

    h.iter()
        .map(|word| format!("{word:08x}"))
        .collect::<String>()
}

fn merge_graph_to_sequence(graph: &GfaGraph) -> Result<MergedSequence, OrgraftError> {
    if graph.segments.is_empty() {
        return Err(OrgraftError::InvalidArgument(
            "cannot merge an empty graph".to_string(),
        ));
    }
    if graph.segments.len() == 1 {
        let segment = graph
            .ordered_segments()
            .first()
            .copied()
            .or_else(|| graph.segments.values().next())
            .unwrap();
        return Ok(MergedSequence {
            sequence: segment.sequence.clone(),
            ordered_nodes: vec![segment.name.clone()],
            oriented_nodes: vec![format!("{}+", segment.name)],
        });
    }

    let order = graph_traversal_order(graph)?;
    let links_by_pair = links_by_pair(graph);
    let orientations = orient_path(&order, &links_by_pair);
    let mut sequence = String::new();
    let mut oriented_nodes = Vec::new();
    for (index, (node_id, orient)) in order.iter().zip(orientations.iter()).enumerate() {
        let segment = graph.segments.get(node_id).ok_or_else(|| {
            OrgraftError::InvalidArgument(format!(
                "merge traversal references missing node `{node_id}`"
            ))
        })?;
        let mut part = if *orient == '-' {
            reverse_complement(&segment.sequence)
        } else {
            segment.sequence.clone()
        };
        if index > 0 {
            let overlap = overlap_between(graph, &order[index - 1], node_id).unwrap_or(0);
            let trim = overlap.min(part.len());
            part = part[trim..].to_string();
        }
        sequence.push_str(&part);
        oriented_nodes.push(format!("{node_id}{orient}"));
    }
    Ok(MergedSequence {
        sequence,
        ordered_nodes: order,
        oriented_nodes,
    })
}

fn graph_traversal_order(graph: &GfaGraph) -> Result<Vec<String>, OrgraftError> {
    let ordered_names = graph.ordered_segment_names();
    let order_index = ordered_names
        .iter()
        .enumerate()
        .map(|(index, name)| (name.clone(), index))
        .collect::<HashMap<_, _>>();
    let mut adjacency: HashMap<String, Vec<String>> = graph
        .segments
        .keys()
        .map(|name| (name.clone(), Vec::new()))
        .collect();
    for link in unique_valid_links(graph) {
        if link.from_name == link.to_name {
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
    for neighbors in adjacency.values_mut() {
        neighbors.sort_by(|left, right| {
            order_index
                .get(left)
                .copied()
                .unwrap_or(usize::MAX)
                .cmp(&order_index.get(right).copied().unwrap_or(usize::MAX))
                .then_with(|| natural_cmp(left, right))
        });
    }

    let mut endpoints = adjacency
        .iter()
        .filter(|(_, neighbors)| neighbors.len() <= 1)
        .map(|(node_id, _)| node_id.clone())
        .collect::<Vec<_>>();
    endpoints.sort_by(|left, right| {
        order_index
            .get(left)
            .copied()
            .unwrap_or(usize::MAX)
            .cmp(&order_index.get(right).copied().unwrap_or(usize::MAX))
            .then_with(|| natural_cmp(left, right))
    });
    let start = endpoints.first().cloned().unwrap_or_else(|| {
        ordered_names
            .first()
            .cloned()
            .or_else(|| graph.segments.keys().next().cloned())
            .unwrap()
    });
    let preferred_second = ordered_names.get(1).cloned();

    let mut order = Vec::new();
    let mut seen = HashSet::new();
    let mut previous: Option<String> = None;
    let mut current = start;
    loop {
        if !seen.insert(current.clone()) {
            break;
        }
        order.push(current.clone());
        let candidates = adjacency
            .get(&current)
            .cloned()
            .unwrap_or_default()
            .into_iter()
            .filter(|neighbor| previous.as_ref() != Some(neighbor))
            .filter(|neighbor| !seen.contains(neighbor))
            .collect::<Vec<_>>();
        if candidates.is_empty() {
            break;
        }
        previous = Some(current);
        current = if order.len() == 1 {
            preferred_second
                .as_ref()
                .filter(|preferred| candidates.contains(preferred))
                .cloned()
                .unwrap_or_else(|| candidates[0].clone())
        } else {
            candidates[0].clone()
        };
    }
    if order.len() != graph.segments.len() {
        return Err(OrgraftError::InvalidArgument(
            "auto-repeat candidate graph could not be traversed as one path or cycle".to_string(),
        ));
    }
    Ok(order)
}

fn overlap_between(graph: &GfaGraph, left: &str, right: &str) -> Option<usize> {
    graph
        .links
        .iter()
        .find(|link| {
            (link.from_name == left && link.to_name == right)
                || (link.from_name == right && link.to_name == left)
        })
        .map(|link| overlap_length_from_cigar(&link.overlap))
}

fn overlap_length_from_cigar(cigar: &str) -> usize {
    if cigar.is_empty() || cigar == "*" {
        return 0;
    }
    let mut total = 0usize;
    let mut digits = String::new();
    for ch in cigar.chars() {
        if ch.is_ascii_digit() {
            digits.push(ch);
        } else {
            if matches!(ch, 'M' | '=' | 'X') {
                total += digits.parse::<usize>().unwrap_or(0);
            }
            digits.clear();
        }
    }
    total
}

fn build_reference_indexes(records: &[FastaRecord]) -> Vec<ReferenceKmerIndex> {
    records
        .iter()
        .map(|record| {
            let sequence = record.sequence.to_ascii_uppercase();
            let indexes = [31usize, 21, 15, 11]
                .into_iter()
                .filter(|kmer| sequence.len() >= *kmer)
                .map(|kmer| KmerIndex {
                    kmer,
                    index: build_circular_kmer_index(&sequence, kmer),
                })
                .collect::<Vec<_>>();
            ReferenceKmerIndex {
                record_id: record.id.clone(),
                sequence,
                indexes,
            }
        })
        .collect()
}

fn score_candidates_against_references(
    candidates: &[AutoRepeatCandidate],
    references: &[ReferenceKmerIndex],
) -> Result<Vec<CandidateScore>, OrgraftError> {
    let mut scores = Vec::new();
    for (candidate_index, candidate) in candidates.iter().enumerate() {
        let candidate_sequence =
            gfa_editor_merge_all_sequence(&candidate.graph)?.to_ascii_uppercase();
        let mut best_for_candidate: Option<CandidateScore> = None;
        for reference in references {
            let sequence_score = score_sequence_arrangement(
                &candidate_sequence,
                &reference.sequence,
                &reference.indexes,
            );
            let score = CandidateScore {
                candidate_id: candidate.id.clone(),
                candidate_index: candidate_index + 1,
                reference_record: reference.record_id.clone(),
                score: sequence_score.score,
                method: sequence_score.method,
                orientation: sequence_score.orientation,
                length_delta: candidate_sequence.len().abs_diff(reference.sequence.len()),
                continuous_bp: sequence_score.continuous_bp,
                continuous_fraction: sequence_score.continuous_fraction,
                diagonal_fraction: sequence_score.diagonal_fraction,
            };
            let replace = best_for_candidate
                .as_ref()
                .map(|old| candidate_score_order(&score, old) == Ordering::Greater)
                .unwrap_or(true);
            if replace {
                best_for_candidate = Some(score);
            }
        }
        if let Some(score) = best_for_candidate {
            scores.push(score);
        }
    }
    Ok(scores)
}

fn best_candidate_score(scores: &[CandidateScore]) -> Option<CandidateScore> {
    scores
        .iter()
        .max_by(|left, right| candidate_score_order(left, right))
        .cloned()
}

fn candidate_score_order(left: &CandidateScore, right: &CandidateScore) -> Ordering {
    left.score
        .partial_cmp(&right.score)
        .unwrap_or(Ordering::Equal)
        .then_with(|| {
            left.continuous_fraction
                .partial_cmp(&right.continuous_fraction)
                .unwrap_or(Ordering::Equal)
        })
        .then_with(|| {
            left.diagonal_fraction
                .partial_cmp(&right.diagonal_fraction)
                .unwrap_or(Ordering::Equal)
        })
        .then_with(|| right.length_delta.cmp(&left.length_delta))
        .then_with(|| right.candidate_index.cmp(&left.candidate_index))
}

fn score_sequence_arrangement(
    candidate_sequence: &str,
    reference_sequence: &str,
    reference_indexes: &[KmerIndex],
) -> SequenceScore {
    if let Some(score) = exact_circular_sequence_score(candidate_sequence, reference_sequence) {
        return score;
    }
    let mut best = SequenceScore {
        score: 0.0,
        method: "sequence-global-kmer-chain".to_string(),
        orientation: '+',
        continuous_bp: 0,
        continuous_fraction: 0.0,
        diagonal_fraction: 0.0,
    };
    for reference_index in reference_indexes {
        if candidate_sequence.len() < reference_index.kmer || reference_index.index.is_empty() {
            continue;
        }
        for (orientation, oriented_sequence) in [
            ('+', candidate_sequence.to_string()),
            ('-', reverse_complement(candidate_sequence)),
        ] {
            let score = global_kmer_chain_score(
                &oriented_sequence,
                reference_sequence,
                &reference_index.index,
                reference_index.kmer,
            );
            let replace = score.score > best.score
                || (score.score == best.score
                    && score.continuous_fraction > best.continuous_fraction)
                || (score.score == best.score
                    && score.continuous_fraction == best.continuous_fraction
                    && orientation == '+'
                    && best.orientation != '+');
            if replace {
                best = SequenceScore {
                    orientation,
                    ..score
                };
            }
        }
        if best.score > 0.0 {
            break;
        }
    }
    best
}

fn exact_circular_sequence_score(
    candidate_sequence: &str,
    reference_sequence: &str,
) -> Option<SequenceScore> {
    if candidate_sequence.len() != reference_sequence.len() || candidate_sequence.is_empty() {
        return None;
    }
    let doubled_reference = format!("{reference_sequence}{reference_sequence}");
    if doubled_reference.contains(candidate_sequence) {
        return Some(SequenceScore {
            score: 1.0,
            method: "sequence-exact-circular".to_string(),
            orientation: '+',
            continuous_bp: candidate_sequence.len(),
            continuous_fraction: 1.0,
            diagonal_fraction: 1.0,
        });
    }
    let reverse_candidate = reverse_complement(candidate_sequence);
    if doubled_reference.contains(&reverse_candidate) {
        return Some(SequenceScore {
            score: 1.0,
            method: "sequence-exact-circular".to_string(),
            orientation: '-',
            continuous_bp: candidate_sequence.len(),
            continuous_fraction: 1.0,
            diagonal_fraction: 1.0,
        });
    }
    None
}

fn build_circular_kmer_index(sequence: &str, kmer_size: usize) -> HashMap<String, Vec<usize>> {
    let mut index: HashMap<String, Vec<usize>> = HashMap::new();
    if sequence.len() < kmer_size {
        return index;
    }
    let stride = (sequence.len() / 500_000).max(1);
    for position in (0..sequence.len()).step_by(stride) {
        let kmer = circular_kmer(sequence, position, kmer_size);
        if !valid_dna_kmer(&kmer) {
            continue;
        }
        let positions = index.entry(kmer).or_default();
        if positions.len() <= 50 {
            positions.push(position);
        }
    }
    index
}

fn circular_kmer(sequence: &str, position: usize, kmer_size: usize) -> String {
    let end = position + kmer_size;
    if end <= sequence.len() {
        sequence[position..end].to_string()
    } else {
        format!(
            "{}{}",
            &sequence[position..],
            &sequence[..end - sequence.len()]
        )
    }
}

fn global_kmer_chain_score(
    candidate_sequence: &str,
    reference_sequence: &str,
    reference_index: &HashMap<String, Vec<usize>>,
    kmer_size: usize,
) -> SequenceScore {
    let reference_length = reference_sequence.len();
    let candidate_limit = candidate_sequence.len().saturating_sub(kmer_size) + 1;
    let stride = (candidate_sequence.len() / 5000).max(1);
    let bin_size = (reference_length / 5000).max(25);
    let mut sampled = 0usize;
    let mut diagonal_bins: HashMap<usize, usize> = HashMap::new();
    let reference_copies = (candidate_sequence.len() / reference_length.max(1) + 2).max(2);
    let mut anchors = Vec::new();

    for query_position in (0..candidate_limit).step_by(stride) {
        let kmer = &candidate_sequence[query_position..query_position + kmer_size];
        if !valid_dna_kmer(kmer) {
            continue;
        }
        sampled += 1;
        let Some(reference_positions) = reference_index.get(kmer) else {
            continue;
        };
        if reference_positions.is_empty() || reference_positions.len() > 50 {
            continue;
        }
        for reference_position in reference_positions {
            let diagonal = (reference_position + reference_length
                - (query_position % reference_length))
                % reference_length;
            let diagonal_bin = diagonal / bin_size;
            *diagonal_bins.entry(diagonal_bin).or_default() += 1;
            for copy_index in 0..reference_copies {
                anchors.push((
                    query_position,
                    reference_position + (copy_index * reference_length),
                    *reference_position,
                ));
            }
        }
    }

    let (_chain_kmers, chain_bp) =
        longest_global_kmer_chain(&anchors, kmer_size, stride, candidate_sequence.len());
    let best_diagonal_count = diagonal_bins
        .iter()
        .map(|(_, count)| *count)
        .max()
        .unwrap_or(0);
    let denominator = candidate_sequence
        .len()
        .max(reference_sequence.len())
        .max(1);
    let continuous_fraction = chain_bp as f64 / denominator as f64;
    let diagonal_fraction = if sampled == 0 {
        0.0
    } else {
        best_diagonal_count as f64 / sampled as f64
    };
    let score = if sampled == 0 {
        0.0
    } else {
        _chain_kmers as f64 / sampled as f64
    };
    SequenceScore {
        score,
        method: format!("sequence-global-kmer-chain-{kmer_size}"),
        orientation: '+',
        continuous_bp: chain_bp,
        continuous_fraction,
        diagonal_fraction,
    }
}

#[derive(Debug, Clone)]
struct GlobalChainState {
    count: usize,
    query_start: usize,
    query_end: usize,
    reference_start: usize,
    reference_end: usize,
}

fn longest_global_kmer_chain(
    anchors: &[(usize, usize, usize)],
    kmer_size: usize,
    stride: usize,
    candidate_len: usize,
) -> (usize, usize) {
    if anchors.is_empty() {
        return (0, 0);
    }
    let mut sorted_anchors = anchors.to_vec();
    sorted_anchors.sort_unstable();
    sorted_anchors.dedup();
    let mut reference_coordinates = sorted_anchors
        .iter()
        .map(|(_, reference_position, _)| *reference_position)
        .collect::<Vec<_>>();
    reference_coordinates.sort_unstable();
    reference_coordinates.dedup();
    let coordinate_rank = reference_coordinates
        .iter()
        .enumerate()
        .map(|(index, coordinate)| (*coordinate, index + 1))
        .collect::<HashMap<_, _>>();
    let mut states: Vec<GlobalChainState> = Vec::new();
    let mut tree: Vec<Option<usize>> = vec![None; reference_coordinates.len() + 1];

    fn state_key(
        states: &[GlobalChainState],
        index: Option<usize>,
        kmer_size: usize,
    ) -> (usize, usize, usize, usize, usize) {
        let Some(index) = index else {
            return (0, 0, 0, 0, 0);
        };
        let state = &states[index];
        (
            state.count,
            state.query_end - state.query_start + kmer_size,
            state.reference_end - state.reference_start + kmer_size,
            usize::MAX - state.query_start,
            usize::MAX - state.reference_start,
        )
    }

    fn better_state(
        states: &[GlobalChainState],
        left: Option<usize>,
        right: Option<usize>,
        kmer_size: usize,
    ) -> Option<usize> {
        if state_key(states, left, kmer_size) >= state_key(states, right, kmer_size) {
            left
        } else {
            right
        }
    }

    fn update_tree(
        tree: &mut [Option<usize>],
        states: &[GlobalChainState],
        mut rank: usize,
        state_index: usize,
        kmer_size: usize,
    ) {
        while rank < tree.len() {
            tree[rank] = better_state(states, Some(state_index), tree[rank], kmer_size);
            rank += rank & rank.wrapping_neg();
        }
    }

    fn query_tree(
        tree: &[Option<usize>],
        states: &[GlobalChainState],
        mut rank: usize,
        kmer_size: usize,
    ) -> Option<usize> {
        let mut best = None;
        while rank > 0 {
            best = better_state(states, best, tree[rank], kmer_size);
            rank -= rank & rank.wrapping_neg();
        }
        best
    }

    let mut index = 0usize;
    let mut best_state = None;
    while index < sorted_anchors.len() {
        let query_position = sorted_anchors[index].0;
        let mut group_updates = Vec::new();
        while index < sorted_anchors.len() && sorted_anchors[index].0 == query_position {
            let (query_pos, reference_pos, _) = sorted_anchors[index];
            let rank = coordinate_rank[&reference_pos];
            let predecessor = query_tree(&tree, &states, rank.saturating_sub(1), kmer_size);
            let state = if let Some(predecessor_index) = predecessor {
                let predecessor_state = &states[predecessor_index];
                GlobalChainState {
                    count: predecessor_state.count + 1,
                    query_start: predecessor_state.query_start,
                    query_end: query_pos,
                    reference_start: predecessor_state.reference_start,
                    reference_end: reference_pos,
                }
            } else {
                GlobalChainState {
                    count: 1,
                    query_start: query_pos,
                    query_end: query_pos,
                    reference_start: reference_pos,
                    reference_end: reference_pos,
                }
            };
            let state_index = states.len();
            states.push(state);
            group_updates.push((rank, state_index));
            best_state = better_state(&states, Some(state_index), best_state, kmer_size);
            index += 1;
        }
        for (rank, state_index) in group_updates {
            update_tree(&mut tree, &states, rank, state_index, kmer_size);
        }
    }
    let Some(best_index) = best_state else {
        return (0, 0);
    };
    let best = &states[best_index];
    let chain_bp = if best.count == 0 {
        0
    } else {
        ((best.count - 1) * stride + kmer_size).min(candidate_len)
    };
    (best.count, chain_bp)
}

fn valid_dna_kmer(kmer: &str) -> bool {
    !kmer.is_empty()
        && kmer
            .as_bytes()
            .iter()
            .all(|base| matches!(base, b'A' | b'C' | b'G' | b'T'))
}

fn resolve_component_with_gfa_editor_cli(
    component: &Component,
    component_graph: &GfaGraph,
    reference_fasta: &Path,
    gfa_editor_cli: &Path,
    max_states: usize,
    max_candidates: usize,
    ready_repeat_nodes: Vec<String>,
) -> Result<(FastaRecord, SubgraphResolution), OrgraftError> {
    let input_gfa = temp_file_path("orgraft-resolve-gfa-editor-input", "gfa");
    let resolved_gfa = temp_file_path("orgraft-resolve-gfa-editor-resolved", "gfa");
    let summary_json = temp_file_path("orgraft-resolve-gfa-editor-summary", "json");
    let history_json = temp_file_path("orgraft-resolve-gfa-editor-history", "json");
    component_graph.write(&input_gfa)?;

    let mut command = Command::new(gfa_editor_cli);
    command
        .arg("auto-repeat")
        .arg(&input_gfa)
        .arg(&resolved_gfa)
        .arg("--reference-fasta")
        .arg(reference_fasta)
        .arg("--summary-json")
        .arg(&summary_json)
        .arg("--history-json")
        .arg(&history_json)
        .arg("--max-states")
        .arg(max_states.to_string())
        .arg("--max-candidates")
        .arg(max_candidates.to_string());
    let completed = command.output()?;
    if !completed.status.success() {
        let _ = fs::remove_file(&input_gfa);
        let _ = fs::remove_file(&resolved_gfa);
        let _ = fs::remove_file(&summary_json);
        let _ = fs::remove_file(&history_json);
        return Err(OrgraftError::InvalidArgument(format!(
            "gfa_editor_cli auto-repeat failed with exit code {:?}: {}",
            completed.status.code(),
            String::from_utf8_lossy(&completed.stderr)
        )));
    }

    let summary = parse_gfa_editor_summary(&summary_json)?;
    let resolved_graph = GfaGraph::read(&resolved_gfa)?;
    let merged = merge_graph_to_sequence(&resolved_graph)?;
    let record = FastaRecord {
        id: component.id.clone(),
        sequence: merged.sequence.clone(),
    };
    let resolution = SubgraphResolution {
        subgraph_id: component.id.clone(),
        reference_alias: summary.reference_alias,
        unresolved_node_count: component.node_ids.len(),
        unresolved_total_bp: component.total_bp,
        node_count: resolved_graph.segments.len(),
        total_bp: merged.sequence.len(),
        resolution_engine: "gfa_editor_cli".to_string(),
        ready_repeat_nodes,
        candidate_count: summary.candidate_count,
        selected_candidate: summary.selected_candidate,
        selected_circular: graph_is_circular_subgraph(&resolved_graph),
        selected_order: summary.selected_order,
        score_method: summary.score_method,
        score_value: summary.score_value,
        score_orientation: summary.score_orientation,
        length_delta: summary.length_delta,
        continuous_bp: summary.continuous_bp,
        ordered_nodes: merged.ordered_nodes,
        oriented_nodes: merged.oriented_nodes,
        missing_reference_hits: Vec::new(),
    };

    let _ = fs::remove_file(&input_gfa);
    let _ = fs::remove_file(&resolved_gfa);
    let _ = fs::remove_file(&summary_json);
    let _ = fs::remove_file(&history_json);
    Ok((record, resolution))
}

fn parse_gfa_editor_summary(path: &Path) -> Result<CliRepeatSummary, OrgraftError> {
    let text = fs::read_to_string(path)?;
    let selected = json_block_after(&text, "\"selected\"", '{', '}').unwrap_or("");
    let best = json_block_after(&text, "\"best\"", '{', '}').unwrap_or("");
    let selected_order_block = json_block_after(selected, "\"order\"", '[', ']').unwrap_or("");
    let selected_order = json_objects(selected_order_block)
        .into_iter()
        .filter_map(|item| {
            let node_id = json_string_value(item, "nodeId")?;
            let strategy = json_string_value(item, "strategy")?;
            let duplicate_id = json_string_value(item, "duplicateId")?;
            Some(format!("{node_id}:{strategy}->{duplicate_id}"))
        })
        .collect::<Vec<_>>();
    Ok(CliRepeatSummary {
        candidate_count: json_usize_value(&text, "candidate_count").unwrap_or(0),
        selected_candidate: json_string_value(selected, "id").unwrap_or_else(|| ".".to_string()),
        selected_order,
        reference_alias: json_string_value(best, "referenceRecord")
            .unwrap_or_else(|| ".".to_string()),
        score_method: json_string_value(best, "method").unwrap_or_else(|| ".".to_string()),
        score_value: json_f64_value(best, "score"),
        score_orientation: json_string_value(best, "orientation")
            .and_then(|value| value.chars().next()),
        length_delta: json_usize_value(best, "lengthDelta"),
        continuous_bp: json_usize_value(best, "continuousBp"),
    })
}

fn json_block_after<'a>(text: &'a str, key: &str, open: char, close: char) -> Option<&'a str> {
    let key_start = text.find(key)?;
    let relative_open = text[key_start..].find(open)?;
    let start = key_start + relative_open;
    json_balanced_body_from(text, start, open, close)
}

fn json_balanced_body_from(text: &str, start: usize, open: char, close: char) -> Option<&str> {
    let mut depth = 0isize;
    for (offset, ch) in text[start..].char_indices() {
        if ch == open {
            depth += 1;
        } else if ch == close {
            depth -= 1;
            if depth == 0 {
                return Some(&text[start + 1..start + offset]);
            }
        }
    }
    None
}

fn json_objects(array_body: &str) -> Vec<&str> {
    let mut objects = Vec::new();
    let mut index = 0usize;
    while let Some(relative) = array_body[index..].find('{') {
        let start = index + relative;
        let Some(body) = json_balanced_body_from(array_body, start, '{', '}') else {
            break;
        };
        objects.push(body);
        index = start + body.len() + 2;
    }
    objects
}

fn json_string_value(text: &str, key: &str) -> Option<String> {
    let raw = json_raw_value(text, key)?;
    let value = raw.strip_prefix('"')?;
    let end = value.find('"')?;
    Some(value[..end].to_string())
}

fn json_usize_value(text: &str, key: &str) -> Option<usize> {
    json_raw_value(text, key)?.parse::<usize>().ok()
}

fn json_f64_value(text: &str, key: &str) -> Option<f64> {
    json_raw_value(text, key)?.parse::<f64>().ok()
}

fn json_raw_value<'a>(text: &'a str, key: &str) -> Option<&'a str> {
    let quoted = format!("\"{key}\"");
    let key_start = text.find(&quoted)?;
    let colon = text[key_start..].find(':')?;
    let value_start = key_start + colon + 1;
    let value = text[value_start..].trim_start();
    if value.starts_with('"') {
        let end = value[1..].find('"')? + 2;
        Some(&value[..end])
    } else {
        let end = value
            .find(|ch: char| ch == ',' || ch == '}' || ch == ']' || ch.is_whitespace())
            .unwrap_or(value.len());
        Some(&value[..end])
    }
}

fn build_reference_ordered_subgraph_fasta(
    reference_fasta: &Path,
    graph: &GfaGraph,
    draft_fasta: &Path,
    blastn: &Path,
    max_states: usize,
    max_candidates: usize,
    gfa_editor_mode: GfaEditorMode,
    gfa_editor_cli: Option<&Path>,
) -> Result<RepeatResolutionResult, OrgraftError> {
    let component_report = graph.component_report();
    let reference_records = read_fasta(reference_fasta)?;
    let reference_indexes = build_reference_indexes(&reference_records);
    let segment_fasta = temp_file_path("orgraft-resolve-merged-segments", "fasta");
    let segment_blast = temp_file_path("orgraft-resolve-segment-placement", "tsv");
    graph.write_fasta(&segment_fasta)?;
    run_blastn(
        blastn,
        reference_fasta,
        &segment_fasta,
        &segment_blast,
        BlastMode::WithLengths,
    )?;
    let hits = parse_blast_hits(&segment_blast)?;
    let _ = fs::remove_file(&segment_fasta);
    let _ = fs::remove_file(&segment_blast);

    let best_overall = best_hits_by_subject(&hits);
    let mut hits_by_subject_reference: HashMap<(String, String), BlastHit> = HashMap::new();
    for hit in &hits {
        let key = (hit.subject_id.clone(), hit.query_id.clone());
        let replace = hits_by_subject_reference
            .get(&key)
            .map(|old| blast_hit_order(hit, old) == Ordering::Greater)
            .unwrap_or(true);
        if replace {
            hits_by_subject_reference.insert(key, hit.clone());
        }
    }

    let mut records = Vec::new();
    let mut resolutions = Vec::new();
    for component in &component_report.components {
        let reference_alias = choose_component_reference(component, &best_overall);
        let component_graph = graph.subgraph(&component.node_ids);
        let ready_repeat_nodes = auto_repeat_ready_node_ids(&component_graph);
        if !ready_repeat_nodes.is_empty() {
            if matches!(gfa_editor_mode, GfaEditorMode::Cli) {
                if let Some(cli) = gfa_editor_cli {
                    let (record, resolution) = resolve_component_with_gfa_editor_cli(
                        component,
                        &component_graph,
                        reference_fasta,
                        cli,
                        max_states,
                        max_candidates,
                        ready_repeat_nodes.clone(),
                    )?;
                    records.push(record);
                    resolutions.push(resolution);
                    continue;
                }
            }
            let (candidates, _warning) = build_auto_repeat_resolution_candidates(
                &component_graph,
                max_states,
                max_candidates,
            )?;
            if !candidates.is_empty() {
                let scores = score_candidates_against_references(&candidates, &reference_indexes)?;
                if let Some(best_score) = best_candidate_score(&scores) {
                    let selected = candidates
                        .iter()
                        .find(|candidate| candidate.id == best_score.candidate_id)
                        .ok_or_else(|| {
                            OrgraftError::InvalidArgument(format!(
                                "selected auto-repeat candidate `{}` is missing",
                                best_score.candidate_id
                            ))
                        })?;
                    let merged = merge_graph_to_sequence(&selected.graph)?;
                    records.push(FastaRecord {
                        id: component.id.clone(),
                        sequence: merged.sequence.clone(),
                    });
                    resolutions.push(SubgraphResolution {
                        subgraph_id: component.id.clone(),
                        reference_alias: best_score.reference_record.clone(),
                        unresolved_node_count: component.node_ids.len(),
                        unresolved_total_bp: component.total_bp,
                        node_count: selected.graph.segments.len(),
                        total_bp: merged.sequence.len(),
                        resolution_engine: "rust".to_string(),
                        ready_repeat_nodes,
                        candidate_count: candidates.len(),
                        selected_candidate: selected.id.clone(),
                        selected_circular: selected.circular,
                        selected_order: selected
                            .order
                            .iter()
                            .map(|step| {
                                format!("{}:{}->{}", step.node_id, step.strategy, step.duplicate_id)
                            })
                            .collect(),
                        score_method: best_score.method.clone(),
                        score_value: Some(best_score.score),
                        score_orientation: Some(best_score.orientation),
                        length_delta: Some(best_score.length_delta),
                        continuous_bp: Some(best_score.continuous_bp),
                        ordered_nodes: merged.ordered_nodes,
                        oriented_nodes: merged.oriented_nodes,
                        missing_reference_hits: Vec::new(),
                    });
                    continue;
                }
            }
        }

        let mut placements = Vec::new();
        for node_id in &component.node_ids {
            let segment = graph.segments.get(node_id).ok_or_else(|| {
                OrgraftError::InvalidArgument(format!(
                    "component refers to missing node `{node_id}`"
                ))
            })?;
            let reference_hit = reference_alias.as_ref().and_then(|alias| {
                hits_by_subject_reference
                    .get(&(node_id.clone(), alias.clone()))
                    .cloned()
            });
            let best = reference_hit.or_else(|| best_overall.get(node_id).cloned());
            let orientation = best.as_ref().map(BlastHit::subject_strand).unwrap_or('+');
            let sequence = if orientation == '-' {
                reverse_complement(&segment.sequence)
            } else {
                segment.sequence.clone()
            };
            placements.push(NodePlacement {
                node_id: node_id.clone(),
                reference_alias: best.as_ref().map(|hit| hit.query_id.clone()),
                reference_start: best.as_ref().map(|hit| hit.query_start),
                orientation,
                sequence,
            });
        }
        placements.sort_by(|left, right| {
            match (left.reference_start, right.reference_start) {
                (Some(a), Some(b)) => a.cmp(&b),
                (Some(_), None) => Ordering::Less,
                (None, Some(_)) => Ordering::Greater,
                (None, None) => natural_cmp(&left.node_id, &right.node_id),
            }
            .then_with(|| natural_cmp(&left.node_id, &right.node_id))
        });
        let sequence = placements
            .iter()
            .map(|placement| placement.sequence.as_str())
            .collect::<String>();
        records.push(FastaRecord {
            id: component.id.clone(),
            sequence,
        });
        let ordered_nodes = placements
            .iter()
            .map(|placement| placement.node_id.clone())
            .collect::<Vec<_>>();
        let oriented_nodes = placements
            .iter()
            .map(|placement| {
                let reference = placement
                    .reference_alias
                    .clone()
                    .unwrap_or_else(|| ".".to_string());
                let start = placement
                    .reference_start
                    .map(|value| value.to_string())
                    .unwrap_or_else(|| ".".to_string());
                format!(
                    "{}{}@{}:{}",
                    placement.node_id, placement.orientation, reference, start
                )
            })
            .collect::<Vec<_>>();
        let missing_reference_hits = placements
            .iter()
            .filter(|placement| placement.reference_start.is_none())
            .map(|placement| placement.node_id.clone())
            .collect::<Vec<_>>();
        resolutions.push(SubgraphResolution {
            subgraph_id: component.id.clone(),
            reference_alias: reference_alias.unwrap_or_else(|| ".".to_string()),
            unresolved_node_count: component.node_ids.len(),
            unresolved_total_bp: component.total_bp,
            node_count: component.node_ids.len(),
            total_bp: component.total_bp,
            resolution_engine: "rust".to_string(),
            ready_repeat_nodes,
            candidate_count: 0,
            selected_candidate: ".".to_string(),
            selected_circular: false,
            selected_order: Vec::new(),
            score_method: "reference-placement-fallback".to_string(),
            score_value: None,
            score_orientation: None,
            length_delta: None,
            continuous_bp: None,
            ordered_nodes,
            oriented_nodes,
            missing_reference_hits,
        });
    }

    write_fasta_records(draft_fasta, &records)?;
    let result = RepeatResolutionResult {
        component_kind: component_report.kind,
        subgraphs: resolutions,
    };
    Ok(result)
}

fn choose_component_reference(
    component: &Component,
    best_hits: &HashMap<String, BlastHit>,
) -> Option<String> {
    let mut scores: HashMap<String, f64> = HashMap::new();
    for node_id in &component.node_ids {
        if let Some(hit) = best_hits.get(node_id) {
            *scores.entry(hit.query_id.clone()).or_default() += hit.bitscore;
        }
    }
    scores
        .into_iter()
        .max_by(|left, right| {
            left.1
                .partial_cmp(&right.1)
                .unwrap_or(Ordering::Equal)
                .then_with(|| right.0.cmp(&left.0))
        })
        .map(|(reference_alias, _)| reference_alias)
}

fn write_id_map(path: &Path, result: &RepeatResolutionResult) -> Result<(), OrgraftError> {
    if let Some(parent) = path.parent() {
        fs::create_dir_all(parent)?;
    }
    let mut out = File::create(path)?;
    writeln!(
        out,
        "rotated_reference_id\tresolved_subgraph_id\tsubgraph_filename"
    )?;
    for subgraph in &result.subgraphs {
        writeln!(
            out,
            "{}\t{}\t{}",
            subgraph.reference_alias,
            subgraph.subgraph_id,
            subgraph_gfa_filename(&subgraph.subgraph_id)
        )?;
    }
    Ok(())
}

#[derive(Debug, Clone)]
struct AlignmentResult {
    aligned_fasta: PathBuf,
    records: Vec<RecordAlignment>,
}

#[derive(Debug, Clone)]
struct RecordAlignment {
    record_id: String,
    reference_alias: String,
    input_length: usize,
    orientation: char,
    reverse_complemented: bool,
    rotation_step: isize,
    best_pident: f64,
    best_length: usize,
    query_start: usize,
    subject_start: usize,
    subject_end: usize,
}

fn align_fasta_to_reference_coordinates(
    reference_fasta: &Path,
    selected_fasta: &Path,
    repeat_resolution: &RepeatResolutionResult,
    aligned_fasta: &Path,
    blastn: &Path,
) -> Result<AlignmentResult, OrgraftError> {
    let reference_by_subgraph = repeat_resolution
        .subgraphs
        .iter()
        .map(|item| (item.subgraph_id.clone(), item.reference_alias.clone()))
        .collect::<HashMap<_, _>>();

    let initial_blast = temp_file_path("orgraft-resolve-subgraph-initial", "tsv");
    run_blastn(
        blastn,
        reference_fasta,
        selected_fasta,
        &initial_blast,
        BlastMode::WithLengths,
    )?;
    let input_records = read_fasta(selected_fasta)?;
    let initial_hits = parse_blast_hits(&initial_blast)?;
    let _ = fs::remove_file(&initial_blast);
    let best_by_subject = best_hits_by_subject_and_reference(&initial_hits, &reference_by_subgraph);

    let oriented_fasta = temp_file_path("orgraft-resolve-oriented", "fasta");
    let mut oriented_records = Vec::new();
    for record in &input_records {
        let best = best_by_subject.get(&record.id);
        let reverse = best.is_some_and(|hit| hit.subject_strand() == '-');
        let sequence = if reverse {
            reverse_complement(&record.sequence)
        } else {
            record.sequence.clone()
        };
        oriented_records.push(FastaRecord {
            id: record.id.clone(),
            sequence,
        });
    }
    write_fasta_records(&oriented_fasta, &oriented_records)?;

    let oriented_blast = temp_file_path("orgraft-resolve-subgraph-oriented", "tsv");
    run_blastn(
        blastn,
        reference_fasta,
        &oriented_fasta,
        &oriented_blast,
        BlastMode::WithLengths,
    )?;
    let oriented_hits = parse_blast_hits(&oriented_blast)?;
    let _ = fs::remove_file(&oriented_fasta);
    let _ = fs::remove_file(&oriented_blast);
    let oriented_best = best_hits_by_subject_and_reference(&oriented_hits, &reference_by_subgraph);

    let mut aligned_records = Vec::new();
    let mut alignment_rows = Vec::new();
    for record in &oriented_records {
        let original = input_records
            .iter()
            .find(|item| item.id == record.id)
            .expect("oriented record comes from input");
        let was_rc = best_by_subject
            .get(&record.id)
            .is_some_and(|hit| hit.subject_strand() == '-');
        let best = oriented_best.get(&record.id);
        let rotation_step = best
            .map(|hit| hit.subject_start as isize - hit.query_start as isize)
            .unwrap_or(0);
        let reference_alias = reference_by_subgraph
            .get(&record.id)
            .cloned()
            .or_else(|| best.map(|hit| hit.query_id.clone()))
            .unwrap_or_else(|| ".".to_string());
        let aligned_sequence = rotate_sequence(&record.sequence, rotation_step);
        aligned_records.push(FastaRecord {
            id: format!(
                "{} [reference={};orientation={};rotation={}]",
                record.id,
                reference_alias,
                if was_rc { "-" } else { "+" },
                rotation_step
            ),
            sequence: aligned_sequence,
        });
        let row = RecordAlignment {
            record_id: original.id.clone(),
            reference_alias,
            input_length: original.sequence.len(),
            orientation: if was_rc { '-' } else { '+' },
            reverse_complemented: was_rc,
            rotation_step,
            best_pident: best.map(|hit| hit.pident).unwrap_or(0.0),
            best_length: best.map(|hit| hit.length).unwrap_or(0),
            query_start: best.map(|hit| hit.query_start).unwrap_or(0),
            subject_start: best.map(|hit| hit.subject_start).unwrap_or(0),
            subject_end: best.map(|hit| hit.subject_end).unwrap_or(0),
        };
        alignment_rows.push(row);
    }
    write_fasta_records(aligned_fasta, &aligned_records)?;
    Ok(AlignmentResult {
        aligned_fasta: aligned_fasta.to_path_buf(),
        records: alignment_rows,
    })
}

fn best_hits_by_subject(hits: &[BlastHit]) -> HashMap<String, BlastHit> {
    let mut best = HashMap::new();
    for hit in hits {
        let replace = best
            .get(&hit.subject_id)
            .map(|old| blast_hit_order(hit, old) == Ordering::Greater)
            .unwrap_or(true);
        if replace {
            best.insert(hit.subject_id.clone(), hit.clone());
        }
    }
    best
}

fn best_hits_by_subject_and_reference(
    hits: &[BlastHit],
    reference_by_subject: &HashMap<String, String>,
) -> HashMap<String, BlastHit> {
    let mut best = HashMap::new();
    for hit in hits {
        if let Some(reference_alias) = reference_by_subject.get(&hit.subject_id) {
            if reference_alias != "." && &hit.query_id != reference_alias {
                continue;
            }
        }
        let replace = best
            .get(&hit.subject_id)
            .map(|old| blast_hit_order(hit, old) == Ordering::Greater)
            .unwrap_or(true);
        if replace {
            best.insert(hit.subject_id.clone(), hit.clone());
        }
    }
    if best.len() < reference_by_subject.len() {
        for hit in hits {
            if best.contains_key(&hit.subject_id) {
                continue;
            }
            let replace = best
                .get(&hit.subject_id)
                .map(|old| blast_hit_order(hit, old) == Ordering::Greater)
                .unwrap_or(true);
            if replace {
                best.insert(hit.subject_id.clone(), hit.clone());
            }
        }
    }
    best
}

#[derive(Debug, Clone)]
struct BlastHit {
    query_id: String,
    subject_id: String,
    pident: f64,
    length: usize,
    query_start: usize,
    subject_start: usize,
    subject_end: usize,
    bitscore: f64,
}

impl BlastHit {
    fn subject_strand(&self) -> char {
        if self.subject_start <= self.subject_end {
            '+'
        } else {
            '-'
        }
    }
}

#[derive(Debug, Clone, Copy)]
enum BlastMode {
    Standard,
    WithLengths,
}

fn run_blastn(
    blastn: &Path,
    query: &Path,
    subject: &Path,
    output: &Path,
    mode: BlastMode,
) -> Result<(), OrgraftError> {
    let outfmt = match mode {
        BlastMode::Standard => "6",
        BlastMode::WithLengths => {
            "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen"
        }
    };
    let output_handle = File::create(output)?;
    let completed = Command::new(blastn)
        .args([
            "-query",
            &query.display().to_string(),
            "-subject",
            &subject.display().to_string(),
            "-outfmt",
            outfmt,
        ])
        .stdout(Stdio::from(output_handle))
        .stderr(Stdio::piped())
        .status()?;
    if !completed.success() {
        return Err(OrgraftError::InvalidArgument(format!(
            "blastn failed for query {} subject {} with exit code {:?}",
            query.display(),
            subject.display(),
            completed.code()
        )));
    }
    Ok(())
}

fn parse_blast_hits(path: &Path) -> Result<Vec<BlastHit>, OrgraftError> {
    let lines = read_nonempty_lines(path)?;
    let mut hits = Vec::new();
    for (index, line) in lines.iter().enumerate() {
        let fields = line.split('\t').collect::<Vec<_>>();
        if fields.len() < 12 {
            continue;
        }
        hits.push(BlastHit {
            query_id: fields[0].to_string(),
            subject_id: fields[1].to_string(),
            pident: parse_f64_field(fields[2], path, index + 1)?,
            length: parse_usize_field(fields[3], path, index + 1)?,
            query_start: parse_usize_field(fields[6], path, index + 1)?,
            subject_start: parse_usize_field(fields[8], path, index + 1)?,
            subject_end: parse_usize_field(fields[9], path, index + 1)?,
            bitscore: parse_f64_field(fields[11], path, index + 1)?,
        });
    }
    Ok(hits)
}

fn parse_usize_field(value: &str, path: &Path, line: usize) -> Result<usize, OrgraftError> {
    value.parse::<usize>().map_err(|_| {
        OrgraftError::InvalidArgument(format!(
            "{}:{line} expected integer BLAST field, got `{value}`",
            path.display()
        ))
    })
}

fn parse_f64_field(value: &str, path: &Path, line: usize) -> Result<f64, OrgraftError> {
    value.parse::<f64>().map_err(|_| {
        OrgraftError::InvalidArgument(format!(
            "{}:{line} expected numeric BLAST field, got `{value}`",
            path.display()
        ))
    })
}

fn blast_hit_order(left: &BlastHit, right: &BlastHit) -> Ordering {
    left.bitscore
        .partial_cmp(&right.bitscore)
        .unwrap_or(Ordering::Equal)
        .then_with(|| left.length.cmp(&right.length))
        .then_with(|| {
            left.pident
                .partial_cmp(&right.pident)
                .unwrap_or(Ordering::Equal)
        })
}

fn write_resolve_report(
    path: &Path,
    options: &ResolveOptions,
    checked_gfa: &Path,
    paths: &OutputPaths,
    topology: &TopologyAudit,
    reference: &ReferenceState,
    rust_prepare: &RustPrepare,
    split_component_gfas: &[PathBuf],
    repeat_resolution: &RepeatResolutionResult,
    alignment: &AlignmentResult,
    elapsed_seconds: f64,
) -> Result<(), OrgraftError> {
    let mut out = File::create(path)?;
    writeln!(out, "# orgraft resolve report")?;
    writeln!(out)?;
    writeln!(out, "## Inputs")?;
    writeln!(out)?;
    writeln!(out, "- checked draft GFA: `{}`", checked_gfa.display())?;
    if let Some(input) = &reference.unrotated_fasta {
        writeln!(out, "- unrotated reference: `{}`", input.display())?;
    }
    writeln!(
        out,
        "- normalized rotated reference: `{}`",
        reference.rotated_fasta.display()
    )?;
    if let Some(input) = &reference.pre_rotated_input_fasta {
        writeln!(
            out,
            "- provided pre-rotated reference: `{}`",
            input.display()
        )?;
    }
    writeln!(out)?;
    writeln!(out, "## Outputs")?;
    writeln!(out)?;
    writeln!(out, "- report: `{}`", paths.report.display())?;
    writeln!(out, "- id map: `{}`", paths.id_map.display())?;
    writeln!(out, "- details: `{}`", paths.details.display())?;
    writeln!(
        out,
        "- merged unresolved GFA: `{}`",
        paths.merged_unresolved_gfa.display()
    )?;
    if !split_component_gfas.is_empty() {
        writeln!(
            out,
            "- split component GFAs: `{}`",
            split_component_gfas
                .iter()
                .map(|path| path.display().to_string())
                .collect::<Vec<_>>()
                .join("`, `")
        )?;
    }
    writeln!(
        out,
        "- resolved subgraphs FASTA: `{}`",
        paths.final_fasta.display()
    )?;
    writeln!(out, "- logs directory: `{}`", paths.logs_dir.display())?;
    writeln!(out, "- graph directory: `{}`", paths.graph_dir.display())?;
    writeln!(out, "- FASTA directory: `{}`", paths.fasta_dir.display())?;
    writeln!(out)?;
    writeln!(out, "## Reference")?;
    writeln!(out)?;
    writeln!(out, "- reference records: {}", reference.records.len())?;
    writeln!(out, "- details file: `{}`", paths.details.display())?;
    for record in &reference.records {
        writeln!(
            out,
            "- `{}` <= `{}`: status `{}`, orientation `{}`, rotation `{}`",
            record.alias, record.unrotated_id, record.status, record.orientation, record.rotation
        )?;
    }
    writeln!(out)?;
    writeln!(out, "## Topology")?;
    writeln!(out)?;
    writeln!(
        out,
        "- checked draft topology: `{}`",
        topology.components.kind
    )?;
    writeln!(
        out,
        "- checked draft graph: {} nodes, {} links",
        topology.report.node_count, topology.report.link_count
    )?;
    writeln!(out, "- subgraphs: {}", topology.components.components.len())?;
    for component in &topology.components.components {
        writeln!(
            out,
            "- `{}`: {} nodes, {} links, {} bp",
            component.id,
            component.node_ids.len(),
            component.link_count,
            component.total_bp
        )?;
    }
    writeln!(out, "- details file: `{}`", paths.details.display())?;
    writeln!(out)?;
    writeln!(out, "## Resolution")?;
    writeln!(out)?;
    writeln!(
        out,
        "- repeat resolution engine mode: `{}`",
        options.gfa_editor_mode.as_str()
    )?;
    writeln!(
        out,
        "- auto merge: `{}` (`{}`), segments {} -> {}, links {} -> {}",
        rust_prepare.merge_action,
        rust_prepare.merge_mode,
        rust_prepare.input_segments,
        rust_prepare.merged_segments,
        rust_prepare.input_links,
        rust_prepare.merged_links
    )?;
    writeln!(
        out,
        "- merged unresolved GFA keeps repeat resolution out of the graph file: `{}`",
        rust_prepare.merged_raw_gfa.display()
    )?;
    if rust_prepare.protected_nodes.is_empty() {
        writeln!(out, "- Protected repeat-like nodes: none")?;
    } else {
        writeln!(
            out,
            "- Protected repeat-like nodes: `{}`",
            rust_prepare.protected_nodes.join(", ")
        )?;
    }
    writeln!(
        out,
        "- auto repeat resolution: 2-in/2-out repeat candidate search with reference-scored selection"
    )?;
    for subgraph in &repeat_resolution.subgraphs {
        if subgraph.candidate_count > 0 {
            writeln!(
                out,
                "- `{}` -> `{}`: engine `{}` selected `{}` from {} candidates; unresolved {} nodes/{} bp, resolved {} nodes/{} bp; score `{}` {} orientation `{}` length_delta {} global_chain_bp {}; ready repeats `{}`; order `{}`",
                subgraph.subgraph_id,
                subgraph.reference_alias,
                subgraph.resolution_engine,
                subgraph.selected_candidate,
                subgraph.candidate_count,
                subgraph.unresolved_node_count,
                subgraph.unresolved_total_bp,
                subgraph.node_count,
                subgraph.total_bp,
                subgraph.score_method,
                subgraph
                    .score_value
                    .map(|value| format!("{value:.6}"))
                    .unwrap_or_else(|| ".".to_string()),
                subgraph
                    .score_orientation
                    .map(|value| value.to_string())
                    .unwrap_or_else(|| ".".to_string()),
                subgraph
                    .length_delta
                    .map(|value| value.to_string())
                    .unwrap_or_else(|| ".".to_string()),
                subgraph
                    .continuous_bp
                    .map(|value| value.to_string())
                    .unwrap_or_else(|| ".".to_string()),
                subgraph.ready_repeat_nodes.join(", "),
                subgraph.selected_order.join(", ")
            )?;
        } else {
            writeln!(
                out,
                "- `{}` -> `{}`: no repeat candidate selected; reference-ordered fallback, {} nodes, {} bp",
                subgraph.subgraph_id,
                subgraph.reference_alias,
                subgraph.node_count,
                subgraph.total_bp
            )?;
        }
    }
    writeln!(out, "- details file: `{}`", paths.details.display())?;
    writeln!(out)?;
    writeln!(out, "## Coordinate Alignment")?;
    writeln!(out)?;
    writeln!(
        out,
        "- aligned FASTA: `{}`",
        alignment.aligned_fasta.display()
    )?;
    writeln!(out, "- details file: `{}`", paths.details.display())?;
    for row in &alignment.records {
        writeln!(
            out,
            "- `{}` vs `{}`: orientation `{}`, rotation `{}`, best {:.3}% over {} bp",
            row.record_id,
            row.reference_alias,
            row.orientation,
            row.rotation_step,
            row.best_pident,
            row.best_length
        )?;
    }
    writeln!(
        out,
        "- Multi-subgraph cross-merge validation is still deferred; subgraphs are resolved independently against their assigned references."
    )?;
    writeln!(out)?;
    writeln!(out, "Elapsed seconds: {:.3}", elapsed_seconds)?;
    Ok(())
}

fn write_resolve_details(
    path: &Path,
    options: &ResolveOptions,
    checked_gfa: &Path,
    paths: &OutputPaths,
    topology: &TopologyAudit,
    reference: &ReferenceState,
    rust_prepare: &RustPrepare,
    split_component_gfas: &[PathBuf],
    repeat_resolution: &RepeatResolutionResult,
    alignment: &AlignmentResult,
    elapsed_seconds: f64,
) -> Result<(), OrgraftError> {
    let mut out = File::create(path)?;
    writeln!(out, "section\tid\tkey\tvalue\tnotes")?;
    write_detail(&mut out, "run", ".", "command", "orgraft resolve", "")?;
    write_detail(
        &mut out,
        "run",
        ".",
        "checked_draft_gfa",
        checked_gfa.display(),
        "",
    )?;
    write_detail(
        &mut out,
        "run",
        ".",
        "out_dir",
        options.out_dir.display(),
        "",
    )?;
    write_detail(
        &mut out,
        "run",
        ".",
        "organelle",
        options.organelle.as_deref().unwrap_or("."),
        "output subdirectory name; no algorithm effect",
    )?;
    write_detail(
        &mut out,
        "run",
        ".",
        "gfa_editor_mode",
        options.gfa_editor_mode.as_str(),
        "",
    )?;
    write_detail(&mut out, "run", ".", "report", paths.report.display(), "")?;
    write_detail(&mut out, "run", ".", "id_map", paths.id_map.display(), "")?;
    write_detail(&mut out, "run", ".", "details", paths.details.display(), "")?;
    write_detail(
        &mut out,
        "run",
        ".",
        "merged_unresolved_gfa",
        paths.merged_unresolved_gfa.display(),
        "",
    )?;
    write_detail(
        &mut out,
        "run",
        ".",
        "split_component_gfas",
        if split_component_gfas.is_empty() {
            ".".to_string()
        } else {
            split_component_gfas
                .iter()
                .map(|path| path.display().to_string())
                .collect::<Vec<_>>()
                .join(",")
        },
        "",
    )?;
    write_detail(
        &mut out,
        "run",
        ".",
        "resolved_fasta",
        paths.final_fasta.display(),
        "",
    )?;
    write_detail(
        &mut out,
        "run",
        ".",
        "logs_dir",
        paths.logs_dir.display(),
        "",
    )?;
    write_detail(
        &mut out,
        "run",
        ".",
        "graph_dir",
        paths.graph_dir.display(),
        "",
    )?;
    write_detail(
        &mut out,
        "run",
        ".",
        "fasta_dir",
        paths.fasta_dir.display(),
        "",
    )?;
    write_detail(
        &mut out,
        "run",
        ".",
        "elapsed_seconds",
        format!("{elapsed_seconds:.3}"),
        "",
    )?;
    write_detail(
        &mut out,
        "topology",
        ".",
        "kind",
        &topology.components.kind,
        "checked draft GFA connected-component classification",
    )?;
    write_detail(
        &mut out,
        "resolution",
        ".",
        "component_kind",
        &repeat_resolution.component_kind,
        "",
    )?;
    write_detail(
        &mut out,
        "topology",
        ".",
        "node_count",
        topology.report.node_count,
        "number of segment records plus linked implicit nodes",
    )?;
    write_detail(
        &mut out,
        "topology",
        ".",
        "link_count",
        topology.report.link_count,
        "number of GFA L records",
    )?;
    let mut class_counts: BTreeMap<&str, usize> = BTreeMap::new();
    for node in &topology.report.nodes {
        *class_counts.entry(node.taxon.code).or_default() += 1;
    }
    for (class_code, count) in class_counts {
        write_detail(
            &mut out,
            "topology",
            ".",
            format!("class:{class_code}"),
            count,
            "node endpoint-degree class count",
        )?;
    }
    for component in &topology.components.components {
        write_detail(
            &mut out,
            "component",
            &component.id,
            "node_count",
            component.node_ids.len(),
            &topology.components.kind,
        )?;
        write_detail(
            &mut out,
            "component",
            &component.id,
            "link_count",
            component.link_count,
            &topology.components.kind,
        )?;
        write_detail(
            &mut out,
            "component",
            &component.id,
            "total_bp",
            component.total_bp,
            &topology.components.kind,
        )?;
        write_detail(
            &mut out,
            "component",
            &component.id,
            "node_ids",
            component.node_ids.join(","),
            &topology.components.kind,
        )?;
    }
    for node in &topology.report.nodes {
        write_detail(
            &mut out,
            "node",
            &node.node_id,
            "endpoint_class",
            node.taxon.code,
            format!(
                "left={};right={};self_links={};{}",
                node.degrees.left, node.degrees.right, node.degrees.self_links, node.taxon.name
            ),
        )?;
    }
    write_detail(
        &mut out,
        "reference",
        ".",
        "record_count",
        reference.records.len(),
        "",
    )?;
    for record in &reference.records {
        write_detail(
            &mut out,
            "reference",
            &record.alias,
            "unrotated_id",
            &record.unrotated_id,
            "",
        )?;
        write_detail(
            &mut out,
            "reference",
            &record.alias,
            "rotated_id",
            &record.rotated_id,
            "",
        )?;
        write_detail(
            &mut out,
            "reference",
            &record.alias,
            "reference_length",
            record.reference_length,
            "",
        )?;
        write_detail(
            &mut out,
            "reference",
            &record.alias,
            "status",
            &record.status,
            "",
        )?;
        write_detail(
            &mut out,
            "reference",
            &record.alias,
            "orientation",
            record.orientation,
            "",
        )?;
        write_detail(
            &mut out,
            "reference",
            &record.alias,
            "sequence_rotation_0based",
            record.rotation,
            "",
        )?;
        write_detail(
            &mut out,
            "reference",
            &record.alias,
            "declared_rotation_0based",
            record
                .declared_rotation
                .map(|value| value.to_string())
                .unwrap_or_else(|| ".".to_string()),
            "",
        )?;
        write_detail(
            &mut out,
            "reference",
            &record.alias,
            "rollback",
            &record.rollback_note,
            "",
        )?;
    }
    write_detail(
        &mut out,
        "reference",
        ".",
        "rotated_fasta",
        reference.rotated_fasta.display(),
        "",
    )?;
    write_detail(
        &mut out,
        "reference",
        ".",
        "unrotated_fasta",
        reference
            .unrotated_fasta
            .as_ref()
            .map(|path| path.display().to_string())
            .unwrap_or_else(|| ".".to_string()),
        "",
    )?;
    write_detail(
        &mut out,
        "prepare",
        ".",
        "merge_action",
        &rust_prepare.merge_action,
        "",
    )?;
    write_detail(
        &mut out,
        "prepare",
        ".",
        "merge_mode",
        &rust_prepare.merge_mode,
        "",
    )?;
    write_detail(
        &mut out,
        "prepare",
        ".",
        "segments",
        format!(
            "{}->{}",
            rust_prepare.input_segments, rust_prepare.merged_segments
        ),
        "",
    )?;
    write_detail(
        &mut out,
        "prepare",
        ".",
        "links",
        format!(
            "{}->{}",
            rust_prepare.input_links, rust_prepare.merged_links
        ),
        "",
    )?;
    write_detail(
        &mut out,
        "prepare",
        ".",
        "protected_repeat_like_nodes",
        if rust_prepare.protected_nodes.is_empty() {
            ".".to_string()
        } else {
            rust_prepare.protected_nodes.join(",")
        },
        "",
    )?;
    for subgraph in &repeat_resolution.subgraphs {
        let action = if subgraph.candidate_count > 0 {
            "auto_repeat_resolution"
        } else {
            "rust_reference_ordered_fallback"
        };
        write_detail(
            &mut out,
            "resolution",
            &subgraph.subgraph_id,
            "action",
            action,
            "",
        )?;
        write_detail(
            &mut out,
            "resolution",
            &subgraph.subgraph_id,
            "reference_alias",
            &subgraph.reference_alias,
            "",
        )?;
        write_detail(
            &mut out,
            "resolution",
            &subgraph.subgraph_id,
            "engine",
            &subgraph.resolution_engine,
            "",
        )?;
        write_detail(
            &mut out,
            "resolution",
            &subgraph.subgraph_id,
            "unresolved_nodes",
            subgraph.unresolved_node_count,
            "",
        )?;
        write_detail(
            &mut out,
            "resolution",
            &subgraph.subgraph_id,
            "unresolved_bp",
            subgraph.unresolved_total_bp,
            "",
        )?;
        write_detail(
            &mut out,
            "resolution",
            &subgraph.subgraph_id,
            "resolved_nodes",
            subgraph.node_count,
            "",
        )?;
        write_detail(
            &mut out,
            "resolution",
            &subgraph.subgraph_id,
            "resolved_bp",
            subgraph.total_bp,
            "",
        )?;
        write_detail(
            &mut out,
            "resolution",
            &subgraph.subgraph_id,
            "ready_repeat_nodes",
            if subgraph.ready_repeat_nodes.is_empty() {
                ".".to_string()
            } else {
                subgraph.ready_repeat_nodes.join(",")
            },
            "",
        )?;
        write_detail(
            &mut out,
            "resolution",
            &subgraph.subgraph_id,
            "candidate_count",
            subgraph.candidate_count,
            "",
        )?;
        write_detail(
            &mut out,
            "resolution",
            &subgraph.subgraph_id,
            "selected_candidate",
            &subgraph.selected_candidate,
            "",
        )?;
        write_detail(
            &mut out,
            "resolution",
            &subgraph.subgraph_id,
            "selected_circular",
            subgraph.selected_circular,
            "",
        )?;
        write_detail(
            &mut out,
            "resolution",
            &subgraph.subgraph_id,
            "selected_order",
            if subgraph.selected_order.is_empty() {
                ".".to_string()
            } else {
                subgraph.selected_order.join(",")
            },
            "",
        )?;
        write_detail(
            &mut out,
            "resolution",
            &subgraph.subgraph_id,
            "score_method",
            &subgraph.score_method,
            "",
        )?;
        write_detail(
            &mut out,
            "resolution",
            &subgraph.subgraph_id,
            "score",
            subgraph
                .score_value
                .map(|value| format!("{value:.6}"))
                .unwrap_or_else(|| ".".to_string()),
            "",
        )?;
        write_detail(
            &mut out,
            "resolution",
            &subgraph.subgraph_id,
            "score_orientation",
            subgraph
                .score_orientation
                .map(|value| value.to_string())
                .unwrap_or_else(|| ".".to_string()),
            "",
        )?;
        write_detail(
            &mut out,
            "resolution",
            &subgraph.subgraph_id,
            "length_delta",
            subgraph
                .length_delta
                .map(|value| value.to_string())
                .unwrap_or_else(|| ".".to_string()),
            "",
        )?;
        write_detail(
            &mut out,
            "resolution",
            &subgraph.subgraph_id,
            "global_chain_bp",
            subgraph
                .continuous_bp
                .map(|value| value.to_string())
                .unwrap_or_else(|| ".".to_string()),
            "longest global k-mer chain span used by reference-scored repeat resolution",
        )?;
        write_detail(
            &mut out,
            "resolution",
            &subgraph.subgraph_id,
            "continuous_bp",
            subgraph
                .continuous_bp
                .map(|value| value.to_string())
                .unwrap_or_else(|| ".".to_string()),
            "compatibility alias for global_chain_bp",
        )?;
        write_detail(
            &mut out,
            "resolution",
            &subgraph.subgraph_id,
            "ordered_nodes",
            subgraph.ordered_nodes.join(","),
            "",
        )?;
        write_detail(
            &mut out,
            "resolution",
            &subgraph.subgraph_id,
            "oriented_nodes",
            subgraph.oriented_nodes.join(","),
            "",
        )?;
        write_detail(
            &mut out,
            "resolution",
            &subgraph.subgraph_id,
            "missing_reference_hits",
            if subgraph.missing_reference_hits.is_empty() {
                ".".to_string()
            } else {
                subgraph.missing_reference_hits.join(",")
            },
            "",
        )?;
    }
    for row in &alignment.records {
        write_detail(
            &mut out,
            "alignment",
            &row.record_id,
            "reference_alias",
            &row.reference_alias,
            "",
        )?;
        write_detail(
            &mut out,
            "alignment",
            &row.record_id,
            "input_length",
            row.input_length,
            "",
        )?;
        write_detail(
            &mut out,
            "alignment",
            &row.record_id,
            "orientation",
            row.orientation,
            "",
        )?;
        write_detail(
            &mut out,
            "alignment",
            &row.record_id,
            "reverse_complemented",
            row.reverse_complemented,
            "",
        )?;
        write_detail(
            &mut out,
            "alignment",
            &row.record_id,
            "rotation_step",
            row.rotation_step,
            "",
        )?;
        write_detail(
            &mut out,
            "alignment",
            &row.record_id,
            "best_pident",
            format!("{:.6}", row.best_pident),
            "",
        )?;
        write_detail(
            &mut out,
            "alignment",
            &row.record_id,
            "best_length",
            row.best_length,
            "",
        )?;
        write_detail(
            &mut out,
            "alignment",
            &row.record_id,
            "query_start",
            row.query_start,
            "",
        )?;
        write_detail(
            &mut out,
            "alignment",
            &row.record_id,
            "subject_start",
            row.subject_start,
            "",
        )?;
        write_detail(
            &mut out,
            "alignment",
            &row.record_id,
            "subject_end",
            row.subject_end,
            "",
        )?;
    }
    Ok(())
}

fn write_detail(
    out: &mut File,
    section: impl std::fmt::Display,
    id: impl std::fmt::Display,
    key: impl std::fmt::Display,
    value: impl std::fmt::Display,
    notes: impl std::fmt::Display,
) -> Result<(), OrgraftError> {
    writeln!(
        out,
        "{}\t{}\t{}\t{}\t{}",
        tsv_cell(section),
        tsv_cell(id),
        tsv_cell(key),
        tsv_cell(value),
        tsv_cell(notes)
    )?;
    Ok(())
}

fn tsv_cell(value: impl std::fmt::Display) -> String {
    value
        .to_string()
        .replace(['\t', '\n', '\r'], " ")
        .trim()
        .to_string()
}

#[derive(Debug, Clone)]
struct FastaRecord {
    id: String,
    sequence: String,
}

fn read_fasta(path: &Path) -> Result<Vec<FastaRecord>, OrgraftError> {
    let text = fs::read_to_string(path).map_err(|error| {
        OrgraftError::InvalidArgument(format!("cannot read {}: {error}", path.display()))
    })?;
    let mut records = Vec::new();
    let mut current_id: Option<String> = None;
    let mut current_seq = String::new();
    for line in text.lines() {
        if let Some(header) = line.strip_prefix('>') {
            if let Some(id) = current_id.replace(header.trim().to_string()) {
                records.push(FastaRecord {
                    id,
                    sequence: current_seq.clone(),
                });
                current_seq.clear();
            }
        } else if !line.trim().is_empty() {
            current_seq.push_str(line.trim());
        }
    }
    if let Some(id) = current_id {
        records.push(FastaRecord {
            id,
            sequence: current_seq,
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

fn write_single_fasta(path: &Path, id: &str, sequence: &str) -> Result<(), OrgraftError> {
    let mut out = File::create(path)?;
    write_fasta_record(&mut out, id, sequence)
}

fn write_fasta_records(path: &Path, records: &[FastaRecord]) -> Result<(), OrgraftError> {
    let mut out = File::create(path)?;
    for record in records {
        write_fasta_record(&mut out, &record.id, &record.sequence)?;
    }
    Ok(())
}

fn write_fasta_record<W: Write>(out: &mut W, id: &str, sequence: &str) -> Result<(), OrgraftError> {
    writeln!(out, ">{id}")?;
    for chunk in sequence.as_bytes().chunks(80) {
        writeln!(out, "{}", String::from_utf8_lossy(chunk))?;
    }
    Ok(())
}

fn rotate_sequence(sequence: &str, step: isize) -> String {
    if sequence.is_empty() {
        return String::new();
    }
    let len = sequence.len() as isize;
    let normalized = ((step % len) + len) % len;
    let split = normalized as usize;
    format!("{}{}", &sequence[split..], &sequence[..split])
}

fn reverse_complement(sequence: &str) -> String {
    sequence.chars().rev().map(complement_base).collect()
}

fn complement_base(base: char) -> char {
    match base {
        'A' => 'T',
        'C' => 'G',
        'G' => 'C',
        'T' => 'A',
        'U' => 'A',
        'R' => 'Y',
        'Y' => 'R',
        'K' => 'M',
        'M' => 'K',
        'S' => 'S',
        'W' => 'W',
        'B' => 'V',
        'D' => 'H',
        'H' => 'D',
        'V' => 'B',
        'N' => 'N',
        'a' => 't',
        'c' => 'g',
        'g' => 'c',
        't' => 'a',
        'u' => 'a',
        'r' => 'y',
        'y' => 'r',
        'k' => 'm',
        'm' => 'k',
        's' => 's',
        'w' => 'w',
        'b' => 'v',
        'd' => 'h',
        'h' => 'd',
        'v' => 'b',
        'n' => 'n',
        other => other,
    }
}

fn read_soft_paths(path: &Path) -> Result<HashMap<String, PathBuf>, OrgraftError> {
    let text = fs::read_to_string(path).map_err(|error| {
        OrgraftError::InvalidArgument(format!("cannot read {}: {error}", path.display()))
    })?;
    let mut tools = HashMap::new();
    for (index, line) in text.lines().enumerate() {
        let line = strip_comment(line).trim();
        if line.is_empty() {
            continue;
        }
        let (name, value) = split_tool_line(line).ok_or_else(|| {
            OrgraftError::InvalidArgument(format!(
                "{}:{} expected software_name<TAB>absolute_path_to_executable",
                path.display(),
                index + 1
            ))
        })?;
        tools.insert(name.to_string(), PathBuf::from(value));
    }
    Ok(tools)
}

fn require_tool(
    soft_paths: &HashMap<String, PathBuf>,
    name: &str,
) -> Result<PathBuf, OrgraftError> {
    let path = soft_paths.get(name).ok_or_else(|| {
        OrgraftError::InvalidArgument(format!("soft_paths.txt is missing `{name}`"))
    })?;
    if !path.exists() {
        return Err(OrgraftError::InvalidArgument(format!(
            "`{name}` not found at {}",
            path.display()
        )));
    }
    Ok(path.clone())
}

fn resolve_gfa_editor_cli(
    soft_paths: &HashMap<String, PathBuf>,
    mode: GfaEditorMode,
) -> Result<Option<PathBuf>, OrgraftError> {
    match mode {
        GfaEditorMode::Rust => Ok(None),
        GfaEditorMode::Cli => require_tool(soft_paths, "gfa_editor_cli").map(Some),
    }
}

fn split_tool_line(line: &str) -> Option<(&str, &str)> {
    line.split_once('\t')
        .or_else(|| line.split_once(char::is_whitespace))
        .map(|(name, path)| (name.trim(), path.trim()))
        .filter(|(name, path)| !name.is_empty() && !path.is_empty())
}

fn strip_comment(line: &str) -> &str {
    line.split_once('#').map(|(value, _)| value).unwrap_or(line)
}

fn canonicalize_existing(path: &Path) -> Result<PathBuf, OrgraftError> {
    fs::canonicalize(path).map_err(|error| {
        OrgraftError::InvalidArgument(format!("cannot read {}: {error}", path.display()))
    })
}

fn temp_file_path(label: &str, extension: &str) -> PathBuf {
    let stamp = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map(|value| value.as_nanos())
        .unwrap_or(0);
    std::env::temp_dir().join(format!("{label}-{}-{stamp}.{extension}", process::id()))
}

fn remove_path_if_exists(path: &Path) -> Result<(), OrgraftError> {
    match fs::symlink_metadata(path) {
        Ok(metadata) if metadata.is_dir() => fs::remove_dir_all(path).map_err(OrgraftError::Io),
        Ok(_) => fs::remove_file(path).map_err(OrgraftError::Io),
        Err(error) if error.kind() == io::ErrorKind::NotFound => Ok(()),
        Err(error) => Err(OrgraftError::Io(error)),
    }
}

fn read_nonempty_lines(path: &Path) -> Result<Vec<String>, OrgraftError> {
    let text = fs::read_to_string(path)?;
    Ok(text
        .lines()
        .map(str::trim)
        .filter(|line| !line.is_empty())
        .map(str::to_string)
        .collect())
}

#[derive(Debug, Clone)]
struct Dsu {
    parent: Vec<usize>,
    rank: Vec<u8>,
}

impl Dsu {
    fn new(size: usize) -> Self {
        Self {
            parent: (0..size).collect(),
            rank: vec![0; size],
        }
    }

    fn find(&mut self, value: usize) -> usize {
        if self.parent[value] != value {
            self.parent[value] = self.find(self.parent[value]);
        }
        self.parent[value]
    }

    fn union(&mut self, left: usize, right: usize) {
        let left_root = self.find(left);
        let right_root = self.find(right);
        if left_root == right_root {
            return;
        }
        if self.rank[left_root] < self.rank[right_root] {
            self.parent[left_root] = right_root;
        } else if self.rank[left_root] > self.rank[right_root] {
            self.parent[right_root] = left_root;
        } else {
            self.parent[right_root] = left_root;
            self.rank[left_root] += 1;
        }
    }
}

fn natural_cmp(left: &str, right: &str) -> Ordering {
    let mut left_parts = NaturalParts::new(left);
    let mut right_parts = NaturalParts::new(right);
    loop {
        match (left_parts.next(), right_parts.next()) {
            (Some(NaturalPart::Number(a)), Some(NaturalPart::Number(b))) => match a.cmp(&b) {
                Ordering::Equal => {}
                ordering => return ordering,
            },
            (Some(NaturalPart::Text(a)), Some(NaturalPart::Text(b))) => match a.cmp(b) {
                Ordering::Equal => {}
                ordering => return ordering,
            },
            (Some(NaturalPart::Number(_)), Some(NaturalPart::Text(_))) => return Ordering::Less,
            (Some(NaturalPart::Text(_)), Some(NaturalPart::Number(_))) => return Ordering::Greater,
            (Some(_), None) => return Ordering::Greater,
            (None, Some(_)) => return Ordering::Less,
            (None, None) => return left.cmp(right),
        }
    }
}

enum NaturalPart<'a> {
    Number(u64),
    Text(&'a str),
}

struct NaturalParts<'a> {
    value: &'a str,
    offset: usize,
}

impl<'a> NaturalParts<'a> {
    fn new(value: &'a str) -> Self {
        Self { value, offset: 0 }
    }

    fn next(&mut self) -> Option<NaturalPart<'a>> {
        if self.offset >= self.value.len() {
            return None;
        }
        let bytes = self.value.as_bytes();
        let start = self.offset;
        let digit = bytes[start].is_ascii_digit();
        self.offset += 1;
        while self.offset < bytes.len() && bytes[self.offset].is_ascii_digit() == digit {
            self.offset += 1;
        }
        let part = &self.value[start..self.offset];
        if digit {
            Some(NaturalPart::Number(part.parse::<u64>().unwrap_or(0)))
        } else {
            Some(NaturalPart::Text(part))
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rotates_sequence_left_by_offset() {
        assert_eq!(rotate_sequence("ABCDE", 2), "CDEAB");
        assert_eq!(rotate_sequence("ABCDE", -1), "EABCD");
    }

    #[test]
    fn detects_plus_rotation() {
        let detected = detect_reference_rotation("ABCDE", "CDEAB").unwrap();
        assert_eq!(detected.0, '+');
        assert_eq!(detected.1, 2);
    }

    #[test]
    fn parses_declared_rotation_from_header() {
        assert_eq!(
            parse_declared_rotation("mito_1 [rotation=293434]"),
            Some(293434)
        );
        assert_eq!(parse_declared_rotation("mito_1"), None);
    }

    #[test]
    fn natural_sort_places_utg2_before_utg10() {
        assert_eq!(natural_cmp("utg2", "utg10"), Ordering::Less);
    }

    fn resolve_args(extra: &[&str]) -> Vec<String> {
        ["--checked-draft-gfa", "draft.gfa"]
            .into_iter()
            .chain(extra.iter().copied())
            .map(str::to_string)
            .collect()
    }

    #[test]
    fn resolve_reference_uses_auto_rotation_mode() {
        let options =
            ResolveOptions::from_args(&resolve_args(&["--reference", "reference.fa"])).unwrap();

        match &options.reference {
            ReferenceInput::AutoRotate(path) => assert_eq!(path, &PathBuf::from("reference.fa")),
            ReferenceInput::PreRotated(_) => panic!("expected auto-rotated reference input"),
        }
    }

    #[test]
    fn resolve_accepts_pre_rotated_reference_without_auto_rotation() {
        let options =
            ResolveOptions::from_args(&resolve_args(&["--pre-rotated-reference", "rotated.fa"]))
                .unwrap();

        match &options.reference {
            ReferenceInput::PreRotated(path) => assert_eq!(path, &PathBuf::from("rotated.fa")),
            ReferenceInput::AutoRotate(_) => panic!("expected pre-rotated reference input"),
        }
    }

    #[test]
    fn resolve_default_out_dir_is_resolve_gfa() {
        let options =
            ResolveOptions::from_args(&resolve_args(&["--pre-rotated-reference", "rotated.fa"]))
                .unwrap();

        assert_eq!(options.out_dir, PathBuf::from("resolve_gfa"));
    }

    #[test]
    fn resolve_pre_rotated_reference_overrides_reference() {
        let options = ResolveOptions::from_args(&resolve_args(&[
            "--reference",
            "unrotated.fa",
            "--pre-rotated-reference",
            "rotated.fa",
        ]))
        .unwrap();

        match &options.reference {
            ReferenceInput::PreRotated(path) => assert_eq!(path, &PathBuf::from("rotated.fa")),
            ReferenceInput::AutoRotate(_) => panic!("expected pre-rotated reference to override"),
        }
    }

    #[test]
    fn resolve_organelle_appends_output_label_only() {
        let options = ResolveOptions::from_args(&resolve_args(&[
            "--pre-rotated-reference",
            "rotated.fa",
            "--out-dir",
            "resolve_gfa",
            "--organelle",
            "plastid",
        ]))
        .unwrap();

        assert_eq!(options.organelle.as_deref(), Some("plastid"));
        assert_eq!(
            options.out_dir,
            PathBuf::from("resolve_gfa").join("plastid")
        );
    }

    #[test]
    fn resolve_rejects_unsafe_organelle_output_label() {
        let error = ResolveOptions::from_args(&resolve_args(&[
            "--pre-rotated-reference",
            "rotated.fa",
            "--organelle",
            "../mito",
        ]))
        .unwrap_err();

        assert_eq!(
            error.to_string(),
            "--organelle expects a simple output label, got `../mito`"
        );
    }

    #[test]
    fn resolve_rejects_removed_candidate_option() {
        let candidate_error = ResolveOptions::from_args(&resolve_args(&[
            "--pre-rotated-reference",
            "rotated.fa",
            "--candidate",
            "1",
        ]))
        .unwrap_err();

        assert_eq!(
            candidate_error.to_string(),
            "unknown orgraft resolve option `--candidate`"
        );
    }

    #[test]
    fn resolve_rejects_removed_gfa_editor_auto_mode() {
        let error = ResolveOptions::from_args(&resolve_args(&[
            "--reference",
            "reference.fa",
            "--gfa-editor-mode",
            "auto",
        ]))
        .unwrap_err();

        assert_eq!(
            error.to_string(),
            "unknown --gfa-editor-mode `auto`; expected rust or cli"
        );
    }

    #[test]
    fn resolve_requires_exactly_one_reference_input() {
        let error = ResolveOptions::from_args(&resolve_args(&[])).unwrap_err();

        assert_eq!(
            error.to_string(),
            "missing --reference FILE (or --pre-rotated-reference FILE)"
        );
    }

    #[test]
    fn rust_auto_repeat_candidates_match_plastid_cross_link_topology() {
        let raw = GfaGraph::parse(
            concat!(
                "S\tutg0\tAAAA\n",
                "S\tutg1\tCCCC\n",
                "S\tutg2\tGGGG\n",
                "L\tutg0\t+\tutg2\t+\t0M\tRC:i:687\tJL:Z:terminal\n",
                "L\tutg0\t-\tutg2\t+\t0M\tRC:i:666\tJL:Z:terminal\n",
                "L\tutg1\t+\tutg2\t-\t0M\tRC:i:838\tJL:Z:terminal\n",
                "L\tutg1\t-\tutg2\t-\t0M\tRC:i:638\tJL:Z:terminal\n",
                "L\tutg2\t+\tutg1\t+\t0M\tRC:i:638\tJL:Z:terminal\n",
                "L\tutg2\t+\tutg1\t-\t0M\tRC:i:838\tJL:Z:terminal\n",
                "L\tutg2\t-\tutg0\t+\t0M\tRC:i:666\tJL:Z:terminal\n",
                "L\tutg2\t-\tutg0\t-\t0M\tRC:i:687\tJL:Z:terminal\n",
            ),
            Path::new("test.gfa"),
        )
        .unwrap();

        assert_eq!(auto_repeat_ready_node_ids(&raw), vec!["utg2".to_string()]);

        let (candidates, warning) =
            build_auto_repeat_resolution_candidates(&raw, 5000, 100).unwrap();

        assert_eq!(warning, None);
        assert_eq!(candidates.len(), 2);
        assert!(candidates.iter().all(|candidate| candidate.circular));
    }

    #[test]
    fn sha256_hex_matches_standard_digest() {
        assert_eq!(
            sha256_hex(b"abc"),
            "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad"
        );
    }

    #[test]
    fn merged_linear_path_carries_length_weighted_depth_tags() {
        let raw = GfaGraph::parse(
            "H\tVN:Z:1.0\nS\ta\tAAA\tLN:i:3\tDP:f:2.0\tAB:i:6\tAC:i:1\nS\tb\tCC\tLN:i:2\tDP:f:5.0\tAB:i:10\tAC:i:2\nL\ta\t+\tb\t+\t0M\n",
            Path::new("test.gfa"),
        )
        .unwrap();

        let (merged, _, _) = merge_unambiguous_gfa(&raw);
        let segment = merged.segments.get("a_b").unwrap();

        assert_eq!(segment.sequence, "AAACC");
        assert!(segment.tags.contains(&"LN:i:5".to_string()));
        assert!(segment.tags.contains(&"DP:f:3.200000".to_string()));
        assert!(segment
            .tags
            .contains(&"CM:Z:raw_node_DP_length_weighted".to_string()));
        assert!(segment.tags.contains(&"RL:i:5".to_string()));
        assert!(segment.tags.contains(&"DB:f:16.000".to_string()));
        assert!(segment.tags.contains(&"AB:i:16".to_string()));
        assert!(segment.tags.contains(&"AC:i:3".to_string()));
        assert!(segment.tags.contains(&"SC:Z:linear_compaction".to_string()));
    }

    #[test]
    fn preserved_node_carries_existing_depth_tags() {
        let raw = GfaGraph::parse(
            "S\tx\tAAAA\tLN:i:4\tDP:f:7.5\tAB:i:30\tAC:i:6\n",
            Path::new("test.gfa"),
        )
        .unwrap();

        let (merged, _, _) = merge_unambiguous_gfa(&raw);
        let segment = merged.segments.get("x").unwrap();

        assert!(segment.tags.contains(&"DP:f:7.500000".to_string()));
        assert!(segment.tags.contains(&"AB:i:30".to_string()));
        assert!(segment.tags.contains(&"AC:i:6".to_string()));
        assert!(segment.tags.contains(&"SC:Z:preserved_node".to_string()));
    }
}
