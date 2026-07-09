use std::collections::{BTreeSet, HashMap, HashSet};
use std::env;
use std::ffi::OsString;
use std::fs::{self, File};
use std::io::{self, BufRead, BufReader, Write};
use std::path::{Path, PathBuf};
use std::time::Instant;

use crate::commands::asm_core;
use crate::commands::shared::{
    print_contract, resolve_gfa_editor_cli, run_gfa_editor_image, CommandContract, GfaImageExport,
};
use crate::error::OrgraftError;

// Command shell for `orgraft asm`: parse user-facing options, prepare the
// top-level output layout, call the assembly core, and publish the final graph.
const HELP: &str = r#"orgraft asm

Conservative draft graph assembly from recruited organelle reads.

Usage:
  orgraft asm --reads FILE --organelle mito|plastid [options]

Inputs:
  --reads FILE              recruited organelle reads
  --organelle NAME          mito or plastid
  --soft-paths FILE         software path table for minimap2 [soft_paths.txt]

Outputs:
  --out-dir DIR             draft assembly directory [results/draft_asm]
  --force                   replace existing step outputs for this organelle

Additional Parameters:
  --profile NAME            low | standard | high [standard]
  --profile-help            show profile presets and advanced assembly parameters
  --threads N               threads passed to the assembly core [8]
  --image-reference-fasta FILE
                            reference FASTA for graph colouring; exports graph.pdf/svg

Layout: OUT/ORGANELLE/{01.input_reads,02.anchor_graph_core,03.finalize_graph,logs}

Finalize note:
  02.anchor_graph_core keeps algorithm/debug graphs unchanged. 03.finalize_graph
  publishes the checkpoint-facing graph and removes reverse-complement duplicate
  L records so one physical connection is counted once by topology checks.
"#;

const PROFILE_HELP: &str = r#"orgraft asm --profile-help

Profile presets and advanced assembly parameters.

Profiles are presets on one 01-08 workflow frame. Shared steps keep the
same file names; unused or inapplicable steps write skipped reports.

  --profile standard | --min-graph-coverage 18 --branch-ratio 0.30 --tip-len 3000 --link-support 20
    normal-depth baseline, about 200-450x; mito uses skeleton-link,
    plastid uses direct-anchor.
  --profile low      | --min-graph-coverage 12 --branch-ratio 0.30 --tip-len 3000 --link-support 20
    low-depth or corrected reads; mito additionally resolves
    04.low_depth_bridge_rescue.
  --profile high     | --subsets ~= 300x * default_genome_size / input_read_bases
    high-depth reads; adds read-depth normalization before the same organelle
    workflow. Genome-size presets: mito 500kb, plastid 150kb.

Workflow frame:
  01.anchor_walk_support       common anchor/read-walk support
  02.unitig_graph              common filtered unitig graph
  03.read_junction_graph       skeleton-link rescue evidence; skipped for direct-anchor
  04.low_depth_bridge_rescue   mito low extra step; otherwise skipped
  05.skeleton_link_evidence    mito skeleton-link evidence; skipped for direct-anchor
  06.repeat_aware_resolution   --stable extra step; otherwise skipped
  07.linked_graph              common selected graph handoff
  08.workflow_summary          reports exact active values and skipped steps

03.finalize_graph publishes the selected graph for checkpoint use. Intermediate
graphs keep raw bidirectional/read-orientation evidence; by default the finalize
copy removes reverse-complement duplicate L records such as A+->B- and B+->A-.
Use `--finalize-dedup-rc-links off` to publish the selected graph unchanged.

Each run writes exact active values to 08.workflow_summary/profile_parameters.tsv

Advanced parameters:
  --stable                  for unstable mitogenomes; enable repeat-aware resolution
  --min-graph-coverage N    override shared anchor/edge coverage floor
  --min-link-ratio FLOAT    optional weak-link ratio filter
  --subsets LIST            manual read subset percent(s) for high-depth reruns
  --finalize-dedup-rc-links on|off  deduplicate RC L records in finalize graph [on]
  --image-reference-fasta FILE       export graph.pdf/svg from finalize graph
  --keep-debug-files        keep .full.* graph companions and subset FASTA files
"#;

const DEFAULT_THREADS: usize = 8;

pub fn run(args: &[String]) -> Result<(), OrgraftError> {
    if args.is_empty() {
        println!("{HELP}");
        return Ok(());
    }

    if args.iter().any(|arg| arg == "-h" || arg == "--help") {
        println!("{HELP}");
        return Ok(());
    }

    if args.iter().any(|arg| arg == "--profile-help") {
        println!("{PROFILE_HELP}");
        return Ok(());
    }

    if args.iter().any(|arg| arg == "--contract") {
        print_contract(&contract());
        return Ok(());
    }

    let options = AsmOptions::from_args(args)?;
    run_asm(&options)
}

fn run_asm(options: &AsmOptions) -> Result<(), OrgraftError> {
    fs::create_dir_all(&options.out_dir)?;
    let summary = run_one_input(options, &options.input)?;
    println!("Wrote {}", summary.output_gfa.display());
    Ok(())
}

fn run_one_input(options: &AsmOptions, input: &AsmInput) -> Result<RunSummary, OrgraftError> {
    let started = Instant::now();
    let label = input.organelle.as_str();
    let profile = input.profile;

    let organelle_dir = options.organelle_dir(input.organelle);
    let input_step = organelle_dir.join("01.input_reads");
    let work_dir = organelle_dir.join("02.anchor_graph_core");
    let finalize_step = organelle_dir.join("03.finalize_graph");
    let logs_dir = organelle_dir.join("logs");
    let output_gfa = finalize_step.join("graph.gfa");

    if !options.force && (work_dir.exists() || output_gfa.exists()) {
        return Err(OrgraftError::InvalidArgument(format!(
            "{} already exists; use --force to replace outputs for {label}",
            work_dir.display()
        )));
    }
    if options.force {
        remove_path_if_exists(&input_step)?;
        remove_path_if_exists(&work_dir)?;
        remove_path_if_exists(&finalize_step)?;
        remove_path_if_exists(&logs_dir)?;
    }

    fs::create_dir_all(&input_step)?;
    fs::create_dir_all(&work_dir)?;
    fs::create_dir_all(&finalize_step)?;
    fs::create_dir_all(&logs_dir)?;
    write_algorithm_notes(&logs_dir.join("algorithm.md"))?;

    let source_reads = canonicalize_existing(&input.reads)?;
    let image_reference_fasta = options
        .image_reference_fasta
        .as_ref()
        .map(|path| canonicalize_existing(path))
        .transpose()?;
    let reads_for_command = prepare_input_link(&source_reads, &input_step)?;
    write_input_manifest(
        &input_step.join("manifest.tsv"),
        input,
        &source_reads,
        &reads_for_command,
    )?;

    let core_request = asm_core::DraftAssemblyRequest {
        organelle: input.organelle.to_core(),
        data_mode: profile.to_core_data_mode(),
        auto_read_subset: profile.auto_read_subset() && options.read_subsets.is_none(),
        repeat_aware_resolution: input.repeat_aware_resolution,
        reads: vec![reads_for_command.clone()],
        out_dir: work_dir.clone(),
        threads: options.threads,
        min_graph_coverage: options.min_graph_coverage,
        min_branch_ratio: options.min_branch_ratio,
        min_tip_len: options.min_tip_len,
        min_link_support: options.min_link_support,
        min_link_ratio: options.min_link_ratio,
        read_subsets: options.read_subsets.clone(),
        keep_debug_files: options.keep_debug_files,
    };
    write_core_request_manifest(&logs_dir.join("asm_core.request.tsv"), &core_request)?;

    let path_guard = PathEnvGuard::set_augmented(&options.soft_paths);
    let core_result = asm_core::run_draft_assembly(core_request);
    drop(path_guard);

    let status_log = logs_dir.join("asm_core.status.log");
    match core_result {
        Ok(()) => fs::write(&status_log, "status\tok\n")?,
        Err(error) => {
            fs::write(&status_log, format!("status\terror\nmessage\t{error}\n"))?;
            return Err(OrgraftError::InvalidArgument(format!(
                "assembly core failed for {label}; see {}",
                status_log.display()
            )));
        }
    }
    write_core_output_manifest(
        &work_dir.join("08.workflow_summary/manifest.tsv"),
        &work_dir,
    )?;

    let selected_graph = find_selected_graph(&work_dir).ok_or_else(|| {
        OrgraftError::InvalidArgument(format!(
            "assembly core did not produce graph.gfa under {}",
            work_dir.display()
        ))
    })?;
    let finalize_stats = copy_finalize_graph(
        &selected_graph,
        &output_gfa,
        options.finalize_dedup_rc_links,
    )?;
    let image_exports = export_finalize_graph_images(
        &output_gfa,
        &finalize_step,
        image_reference_fasta.as_deref(),
        &options.soft_paths,
    );
    for row in image_exports.iter().filter(|row| row.status != "written") {
        if image_reference_fasta.is_some() {
            eprintln!(
                "Warning: optional GFA_Editor {} export {} for {}; see {}",
                row.format,
                row.status,
                row.output.display(),
                finalize_step.join("manifest.tsv").display()
            );
        }
    }
    write_finalize_manifest(
        &finalize_step.join("manifest.tsv"),
        &selected_graph,
        &output_gfa,
        finalize_stats,
        options.finalize_dedup_rc_links,
        image_reference_fasta.as_deref(),
        &image_exports,
    )?;

    let summary = RunSummary {
        organelle: input.organelle,
        profile,
        repeat_aware_resolution: input.repeat_aware_resolution,
        reads: source_reads,
        work_dir,
        selected_graph,
        output_gfa,
        elapsed_seconds: started.elapsed().as_secs_f64(),
    };
    write_run_manifest(&logs_dir.join("run.tsv"), &[summary.clone()])?;
    Ok(summary)
}

fn contract() -> CommandContract {
    CommandContract {
        command: "asm",
        origin: "refactored internal anchor graph construction logic",
        purpose: "build conservative plant organelle draft graphs from recruited reads",
        inputs: &[
            "--reads FILE",
            "--organelle mito|plastid",
            "minimap2 executable for skeleton-link workflows",
        ],
        outputs: &[
            "results/draft_asm/<organelle>/01.input_reads",
            "results/draft_asm/<organelle>/02.anchor_graph_core",
            "results/draft_asm/<organelle>/03.finalize_graph/graph.gfa",
            "optional results/draft_asm/<organelle>/03.finalize_graph/graph.pdf",
            "optional results/draft_asm/<organelle>/03.finalize_graph/graph.svg",
            "results/draft_asm/<organelle>/logs/*.log",
        ],
        notes: &[
            "OrgRAFT runs one explicit profile for one recruited read input",
            "mito defaults to the standard skeleton-link workflow",
            "plastid defaults to the standard direct-anchor workflow",
            "graph PDF/SVG export runs only when --image-reference-fasta is provided",
            "polishing is intentionally not part of this draft assembly command",
        ],
    }
}

#[derive(Debug, Clone)]
struct AsmOptions {
    input: AsmInput,
    out_dir: PathBuf,
    soft_paths: HashMap<String, PathBuf>,
    threads: usize,
    force: bool,
    min_graph_coverage: Option<u32>,
    min_branch_ratio: Option<f64>,
    min_tip_len: Option<usize>,
    min_link_support: Option<u32>,
    min_link_ratio: Option<f64>,
    read_subsets: Option<Vec<u16>>,
    finalize_dedup_rc_links: bool,
    image_reference_fasta: Option<PathBuf>,
    keep_debug_files: bool,
}

impl AsmOptions {
    fn from_args(args: &[String]) -> Result<Self, OrgraftError> {
        let mut out_dir = PathBuf::from("results/draft_asm");
        let mut soft_paths_file = PathBuf::from("soft_paths.txt");
        let mut threads = DEFAULT_THREADS;
        let mut force = false;
        let mut min_graph_coverage = None;
        let mut min_branch_ratio = None;
        let mut min_tip_len = None;
        let mut min_link_support = None;
        let mut min_link_ratio = None;
        let mut read_subsets = None;
        let mut finalize_dedup_rc_links = true;
        let mut image_reference_fasta = None;
        let mut keep_debug_files = false;
        let mut stable = false;
        let mut profile = AsmProfile::Standard;
        let mut reads: Option<PathBuf> = None;
        let mut organelle: Option<Organelle> = None;

        let mut index = 0;
        while index < args.len() {
            match args[index].as_str() {
                "--reads" | "-i" => {
                    let name = args[index].clone();
                    if reads.is_some() {
                        return Err(OrgraftError::InvalidArgument(
                            "--reads may be provided only once".to_string(),
                        ));
                    }
                    reads = Some(PathBuf::from(value_after(args, &mut index, &name)?));
                }
                "--organelle" => {
                    organelle = Some(parse_organelle(value_after(
                        args,
                        &mut index,
                        "--organelle",
                    )?)?);
                }
                "--out-dir" => {
                    out_dir = PathBuf::from(value_after(args, &mut index, "--out-dir")?);
                }
                "--soft-paths" => {
                    soft_paths_file = PathBuf::from(value_after(args, &mut index, "--soft-paths")?);
                }
                "--threads" | "-t" => {
                    let name = args[index].clone();
                    threads = parse_usize(value_after(args, &mut index, &name)?, &name)?;
                    if threads == 0 {
                        return Err(OrgraftError::InvalidArgument(
                            "--threads must be greater than 0".to_string(),
                        ));
                    }
                }
                "--min-graph-coverage" => {
                    let value = parse_u32(
                        value_after(args, &mut index, "--min-graph-coverage")?,
                        "--min-graph-coverage",
                    )?;
                    if value == 0 {
                        return Err(OrgraftError::InvalidArgument(
                            "--min-graph-coverage must be greater than 0".to_string(),
                        ));
                    }
                    min_graph_coverage = Some(value);
                }
                "--branch-ratio" | "--min-branch-ratio" => {
                    let name = args[index].clone();
                    let value = parse_f64(value_after(args, &mut index, &name)?, &name)?;
                    if !(0.0..=1.0).contains(&value) {
                        return Err(OrgraftError::InvalidArgument(format!(
                            "{name} must be between 0 and 1"
                        )));
                    }
                    min_branch_ratio = Some(value);
                }
                "--tip-len" | "--min-tip-len" => {
                    let name = args[index].clone();
                    let value = parse_usize(value_after(args, &mut index, &name)?, &name)?;
                    min_tip_len = Some(value);
                }
                "--link-support" | "--min-link-support" => {
                    let name = args[index].clone();
                    let value = parse_u32(value_after(args, &mut index, &name)?, &name)?;
                    if value == 0 {
                        return Err(OrgraftError::InvalidArgument(format!(
                            "{name} must be greater than 0"
                        )));
                    }
                    min_link_support = Some(value);
                }
                "--min-link-ratio" => {
                    let value = parse_f64(
                        value_after(args, &mut index, "--min-link-ratio")?,
                        "--min-link-ratio",
                    )?;
                    if !(0.0..=1.0).contains(&value) {
                        return Err(OrgraftError::InvalidArgument(
                            "--min-link-ratio must be between 0 and 1".to_string(),
                        ));
                    }
                    min_link_ratio = Some(value);
                }
                value
                    if value.starts_with("--subsets=") || value.starts_with("--read-subsets=") =>
                {
                    let (name, subset_spec) = value.split_once('=').expect("matched assignment");
                    read_subsets = Some(parse_subset_list(subset_spec, name)?);
                }
                "--subsets" | "--read-subsets" => {
                    let name = args[index].clone();
                    read_subsets = Some(parse_subset_list(
                        value_after(args, &mut index, &name)?,
                        &name,
                    )?);
                }
                "--profile" => {
                    profile = parse_profile(value_after(args, &mut index, "--profile")?)?;
                }
                "--finalize-dedup-rc-links" => {
                    finalize_dedup_rc_links = parse_on_off(
                        value_after(args, &mut index, "--finalize-dedup-rc-links")?,
                        "--finalize-dedup-rc-links",
                    )?;
                }
                "--image-reference-fasta" => {
                    image_reference_fasta = Some(PathBuf::from(value_after(
                        args,
                        &mut index,
                        "--image-reference-fasta",
                    )?));
                }
                "--stable" => stable = true,
                "--keep-debug-files" => keep_debug_files = true,
                "--force" => force = true,
                other => {
                    return Err(OrgraftError::InvalidArgument(format!(
                        "unknown asm option `{other}`"
                    )));
                }
            }
            index += 1;
        }

        let reads = reads
            .ok_or_else(|| OrgraftError::InvalidArgument("missing --reads FILE".to_string()))?;
        let organelle = organelle.ok_or_else(|| {
            OrgraftError::InvalidArgument("missing --organelle mito|plastid".to_string())
        })?;
        if stable {
            if organelle != Organelle::Mito {
                return Err(OrgraftError::InvalidArgument(
                    "--stable is only available with --organelle mito".to_string(),
                ));
            }
        }

        let soft_paths = read_soft_paths_optional(&soft_paths_file)?;

        Ok(Self {
            input: AsmInput {
                organelle,
                reads,
                profile,
                repeat_aware_resolution: stable,
            },
            out_dir,
            soft_paths,
            threads,
            force,
            min_graph_coverage,
            min_branch_ratio,
            min_tip_len,
            min_link_support,
            min_link_ratio,
            read_subsets,
            finalize_dedup_rc_links,
            image_reference_fasta,
            keep_debug_files,
        })
    }

    fn organelle_dir(&self, organelle: Organelle) -> PathBuf {
        self.out_dir.join(organelle.as_str())
    }
}

#[derive(Debug, Clone)]
struct AsmInput {
    organelle: Organelle,
    reads: PathBuf,
    profile: AsmProfile,
    repeat_aware_resolution: bool,
}

#[derive(Debug, Clone)]
struct RunSummary {
    organelle: Organelle,
    profile: AsmProfile,
    repeat_aware_resolution: bool,
    reads: PathBuf,
    work_dir: PathBuf,
    selected_graph: PathBuf,
    output_gfa: PathBuf,
    elapsed_seconds: f64,
}

#[derive(Debug, Clone, Copy, Eq, PartialEq)]
enum Organelle {
    Mito,
    Plastid,
}

impl Organelle {
    fn as_str(self) -> &'static str {
        match self {
            Self::Mito => "mito",
            Self::Plastid => "plastid",
        }
    }

    fn to_core(self) -> asm_core::DraftOrganelle {
        match self {
            Self::Mito => asm_core::DraftOrganelle::Mito,
            Self::Plastid => asm_core::DraftOrganelle::Plastid,
        }
    }
}

#[derive(Debug, Clone, Copy, Eq, PartialEq)]
enum AsmProfile {
    Low,
    Standard,
    High,
}

impl AsmProfile {
    fn as_str(self) -> &'static str {
        match self {
            Self::Low => "low",
            Self::Standard => "standard",
            Self::High => "high",
        }
    }

    fn to_core_data_mode(self) -> asm_core::DraftDataMode {
        match self {
            Self::Low => asm_core::DraftDataMode::Low,
            Self::Standard | Self::High => asm_core::DraftDataMode::Standard,
        }
    }

    fn auto_read_subset(self) -> bool {
        self == Self::High
    }
}

fn value_after<'a>(
    args: &'a [String],
    index: &mut usize,
    name: &str,
) -> Result<&'a str, OrgraftError> {
    *index += 1;
    let Some(value) = args.get(*index) else {
        return Err(OrgraftError::InvalidArgument(format!(
            "missing value for {name}"
        )));
    };
    if value.starts_with("--") {
        return Err(OrgraftError::InvalidArgument(format!(
            "missing value for {name}"
        )));
    }
    Ok(value)
}

fn parse_usize(value: &str, name: &str) -> Result<usize, OrgraftError> {
    value.parse::<usize>().map_err(|_| {
        OrgraftError::InvalidArgument(format!("invalid {name}: expected a positive integer"))
    })
}

fn parse_u32(value: &str, name: &str) -> Result<u32, OrgraftError> {
    value.parse::<u32>().map_err(|_| {
        OrgraftError::InvalidArgument(format!("invalid {name}: expected a positive integer"))
    })
}

fn parse_f64(value: &str, name: &str) -> Result<f64, OrgraftError> {
    value
        .parse::<f64>()
        .map_err(|_| OrgraftError::InvalidArgument(format!("invalid {name}: expected a number")))
}

fn parse_on_off(value: &str, name: &str) -> Result<bool, OrgraftError> {
    match value {
        "on" | "true" => Ok(true),
        "off" | "false" => Ok(false),
        _ => Err(OrgraftError::InvalidArgument(format!(
            "{name} expects on/off"
        ))),
    }
}

fn parse_subset_list(value: &str, name: &str) -> Result<Vec<u16>, OrgraftError> {
    let mut subsets = Vec::new();
    for raw_part in value.split(',') {
        let part = raw_part.trim();
        if part.is_empty() {
            return Err(OrgraftError::InvalidArgument(format!(
                "invalid {name}: expected comma-separated percentages"
            )));
        }
        let percent = part.parse::<f64>().map_err(|_| {
            OrgraftError::InvalidArgument(format!("invalid {name}: expected percentages"))
        })?;
        if !(0.0..=100.0).contains(&percent) || percent == 0.0 {
            return Err(OrgraftError::InvalidArgument(format!(
                "{name} values must be between 1 and 100"
            )));
        }
        let basis_points = (percent * 100.0).round();
        if (basis_points / 100.0 - percent).abs() > 1e-9 {
            return Err(OrgraftError::InvalidArgument(format!(
                "{name} values support at most two decimal places"
            )));
        }
        subsets.push(basis_points as u16);
    }
    subsets.sort_unstable();
    subsets.dedup();
    Ok(subsets)
}

fn parse_profile(value: &str) -> Result<AsmProfile, OrgraftError> {
    match value {
        "low" => Ok(AsmProfile::Low),
        "standard" => Ok(AsmProfile::Standard),
        "high" => Ok(AsmProfile::High),
        other => Err(OrgraftError::InvalidArgument(format!(
            "unknown asm profile `{other}`; expected low, standard, or high"
        ))),
    }
}

fn parse_organelle(value: &str) -> Result<Organelle, OrgraftError> {
    match value {
        "mito" | "mitochondria" | "mitochondrion" | "mt" => Ok(Organelle::Mito),
        "plastid" | "plasti" | "chloroplast" | "cp" => Ok(Organelle::Plastid),
        other => Err(OrgraftError::InvalidArgument(format!(
            "unknown organelle `{other}`; expected mito or plastid"
        ))),
    }
}

fn canonicalize_existing(path: &Path) -> Result<PathBuf, OrgraftError> {
    fs::canonicalize(path).map_err(|error| {
        OrgraftError::InvalidArgument(format!("cannot read {}: {error}", path.display()))
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

fn prepare_input_link(source: &Path, input_step: &Path) -> Result<PathBuf, OrgraftError> {
    let link_name = if source
        .file_name()
        .and_then(|name| name.to_str())
        .is_some_and(|name| name.ends_with(".gz"))
    {
        "reads.fastq.gz"
    } else {
        "reads.fastq"
    };
    let link = input_step.join(link_name);
    remove_path_if_exists(&link)?;

    #[cfg(unix)]
    {
        std::os::unix::fs::symlink(source, &link)?;
        Ok(link)
    }

    #[cfg(not(unix))]
    {
        let note = input_step.join("reads.symlink.unavailable.txt");
        fs::write(
            &note,
            format!("Using original reads path: {}\n", source.display()),
        )?;
        Ok(source.to_path_buf())
    }
}

fn write_input_manifest(
    path: &Path,
    input: &AsmInput,
    source_reads: &Path,
    command_reads: &Path,
) -> Result<(), OrgraftError> {
    let mut out = File::create(path)?;
    writeln!(out, "key\tvalue")?;
    writeln!(out, "organelle\t{}", input.organelle.as_str())?;
    writeln!(out, "profile\t{}", input.profile.as_str())?;
    writeln!(
        out,
        "repeat_aware_resolution\t{}",
        input.repeat_aware_resolution
    )?;
    writeln!(out, "source_reads\t{}", source_reads.display())?;
    writeln!(out, "command_reads\t{}", command_reads.display())?;
    Ok(())
}

fn write_core_request_manifest(
    path: &Path,
    request: &asm_core::DraftAssemblyRequest,
) -> Result<(), OrgraftError> {
    let mut out = File::create(path)?;
    writeln!(out, "key\tvalue")?;
    writeln!(out, "organelle\t{}", core_organelle_name(request.organelle))?;
    writeln!(out, "data_mode\t{}", core_data_mode_name(request.data_mode))?;
    writeln!(out, "auto_read_subset\t{}", request.auto_read_subset)?;
    writeln!(out, "out_dir\t{}", request.out_dir.display())?;
    writeln!(out, "threads\t{}", request.threads)?;
    if let Some(min_graph_coverage) = request.min_graph_coverage {
        writeln!(out, "min_graph_coverage\t{min_graph_coverage}")?;
    }
    if let Some(min_branch_ratio) = request.min_branch_ratio {
        writeln!(out, "min_branch_ratio\t{min_branch_ratio}")?;
    }
    if let Some(min_tip_len) = request.min_tip_len {
        writeln!(out, "min_tip_len\t{min_tip_len}")?;
    }
    if let Some(min_link_support) = request.min_link_support {
        writeln!(out, "min_link_support\t{min_link_support}")?;
    }
    if let Some(min_link_ratio) = request.min_link_ratio {
        writeln!(out, "min_link_ratio\t{min_link_ratio}")?;
    }
    writeln!(
        out,
        "repeat_aware_resolution\t{}",
        request.repeat_aware_resolution
    )?;
    if let Some(read_subsets) = &request.read_subsets {
        writeln!(out, "subsets\t{}", format_subset_list(read_subsets))?;
    }
    writeln!(out, "keep_debug_files\t{}", request.keep_debug_files)?;
    for (index, read) in request.reads.iter().enumerate() {
        writeln!(out, "reads_{}\t{}", index + 1, read.display())?;
    }
    Ok(())
}

fn format_subset_list(subsets: &[u16]) -> String {
    subsets
        .iter()
        .map(|basis_points| {
            let whole = basis_points / 100;
            let decimal = basis_points % 100;
            if decimal == 0 {
                whole.to_string()
            } else if decimal % 10 == 0 {
                format!("{}.{}", whole, decimal / 10)
            } else {
                format!("{}.{:02}", whole, decimal)
            }
        })
        .collect::<Vec<_>>()
        .join(",")
}

fn core_organelle_name(organelle: asm_core::DraftOrganelle) -> &'static str {
    match organelle {
        asm_core::DraftOrganelle::Mito => "mito",
        asm_core::DraftOrganelle::Plastid => "plastid",
    }
}

fn core_data_mode_name(data_mode: asm_core::DraftDataMode) -> &'static str {
    match data_mode {
        asm_core::DraftDataMode::Low => "low",
        asm_core::DraftDataMode::Standard => "standard",
    }
}

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
struct FinalizeGraphStats {
    input_links: usize,
    output_links: usize,
    rc_duplicate_links_removed: usize,
}

fn copy_finalize_graph(
    selected_graph: &Path,
    output_gfa: &Path,
    dedup_rc_links: bool,
) -> Result<FinalizeGraphStats, OrgraftError> {
    let input = BufReader::new(File::open(selected_graph)?);
    let mut output = File::create(output_gfa)?;
    let mut seen_links = HashSet::new();
    let mut stats = FinalizeGraphStats::default();

    for (line_index, line_result) in input.lines().enumerate() {
        let line_number = line_index + 1;
        let line = line_result?;
        if let Some(key) = canonical_finalize_link_key(&line, line_number)? {
            stats.input_links += 1;
            if dedup_rc_links && !seen_links.insert(key) {
                stats.rc_duplicate_links_removed += 1;
                continue;
            }
            stats.output_links += 1;
        }
        writeln!(output, "{line}")?;
    }

    Ok(stats)
}

fn export_finalize_graph_images(
    output_gfa: &Path,
    finalize_step: &Path,
    image_reference_fasta: Option<&Path>,
    soft_paths: &HashMap<String, PathBuf>,
) -> Vec<GfaImageExport> {
    let outputs = [
        ("pdf", finalize_step.join("graph.pdf")),
        ("svg", finalize_step.join("graph.svg")),
    ];
    let Some(image_reference_fasta) = image_reference_fasta else {
        return skipped_finalize_graph_images(&outputs, "skipped_no_image_reference");
    };
    let gfa_editor_cli = match resolve_gfa_editor_cli(soft_paths) {
        Ok(path) => path,
        Err(error) => return skipped_finalize_graph_images(&outputs, &error),
    };
    outputs
        .iter()
        .map(|(format, output_path)| {
            run_gfa_editor_image(
                &gfa_editor_cli,
                soft_paths,
                output_gfa,
                output_path,
                image_reference_fasta,
                format,
            )
        })
        .collect()
}

fn skipped_finalize_graph_images(outputs: &[(&str, PathBuf)], reason: &str) -> Vec<GfaImageExport> {
    let status = if reason == "skipped_no_image_reference" {
        "skipped_no_image_reference"
    } else {
        "skipped_missing_gfa_editor_cli"
    };
    outputs
        .iter()
        .map(|(format, output_path)| GfaImageExport {
            format: (*format).to_string(),
            output: output_path.clone(),
            command: ".".to_string(),
            status: status.to_string(),
            stdout: String::new(),
            stderr: reason.to_string(),
        })
        .collect()
}

fn canonical_finalize_link_key(
    line: &str,
    line_number: usize,
) -> Result<Option<(String, char, String, char)>, OrgraftError> {
    let fields: Vec<&str> = line.split('\t').collect();
    if fields.first().copied() != Some("L") {
        return Ok(None);
    }
    if fields.len() < 5 {
        return Err(OrgraftError::InvalidArgument(format!(
            "GFA line {line_number}: link record is missing required fields"
        )));
    }

    let from = fields[1];
    let from_orient = parse_gfa_orientation(fields[2], line_number)?;
    let to = fields[3];
    let to_orient = parse_gfa_orientation(fields[4], line_number)?;
    let forward = (from.to_string(), from_orient, to.to_string(), to_orient);
    let reverse = (
        to.to_string(),
        flip_gfa_orientation(to_orient),
        from.to_string(),
        flip_gfa_orientation(from_orient),
    );
    Ok(Some(if forward <= reverse { forward } else { reverse }))
}

fn parse_gfa_orientation(value: &str, line_number: usize) -> Result<char, OrgraftError> {
    match value {
        "+" => Ok('+'),
        "-" => Ok('-'),
        _ => Err(OrgraftError::InvalidArgument(format!(
            "GFA line {line_number}: invalid orientation `{value}`"
        ))),
    }
}

fn flip_gfa_orientation(orientation: char) -> char {
    if orientation == '+' {
        '-'
    } else {
        '+'
    }
}

fn write_finalize_manifest(
    path: &Path,
    selected_graph: &Path,
    output_gfa: &Path,
    stats: FinalizeGraphStats,
    dedup_rc_links: bool,
    image_reference_fasta: Option<&Path>,
    image_exports: &[GfaImageExport],
) -> Result<(), OrgraftError> {
    let mut out = File::create(path)?;
    writeln!(out, "key\tvalue")?;
    writeln!(out, "selected_graph\t{}", selected_graph.display())?;
    writeln!(out, "output_gfa\t{}", output_gfa.display())?;
    writeln!(
        out,
        "finalize_dedup_rc_links\t{}",
        if dedup_rc_links { "on" } else { "off" }
    )?;
    writeln!(out, "input_links\t{}", stats.input_links)?;
    writeln!(out, "output_links\t{}", stats.output_links)?;
    writeln!(
        out,
        "rc_duplicate_links_removed\t{}",
        stats.rc_duplicate_links_removed
    )?;
    writeln!(
        out,
        "image_reference_fasta\t{}",
        image_reference_fasta
            .map(|path| path.display().to_string())
            .unwrap_or_else(|| ".".to_string())
    )?;
    for row in image_exports {
        writeln!(out, "image_{}_output\t{}", row.format, row.output.display())?;
        writeln!(out, "image_{}_status\t{}", row.format, row.status)?;
        writeln!(
            out,
            "image_{}_command\t{}",
            row.format,
            row.command.replace('\n', "\\n")
        )?;
        writeln!(
            out,
            "image_{}_stdout\t{}",
            row.format,
            row.stdout.replace('\n', "\\n")
        )?;
        writeln!(
            out,
            "image_{}_stderr\t{}",
            row.format,
            row.stderr.replace('\n', "\\n")
        )?;
    }
    Ok(())
}

fn write_core_output_manifest(path: &Path, work_dir: &Path) -> Result<(), OrgraftError> {
    if let Some(parent) = path.parent() {
        fs::create_dir_all(parent)?;
    }
    let expected = [
        (
            "01.read_subset_selection/summary.tsv",
            "01.read_subset_selection",
            "deterministic read-subset selection summary for high profiles",
        ),
        (
            "01.anchor_walk_support/anchors.tsv",
            "01.anchor_walk_support",
            "anchor nodes retained before graph compression",
        ),
        (
            "01.anchor_walk_support/edges.tsv",
            "01.anchor_walk_support",
            "directed anchor-to-anchor edge support",
        ),
        (
            "02.unitig_graph/unitigs.fasta",
            "02.unitig_graph",
            "compressed non-branching graph paths",
        ),
        (
            "02.unitig_graph/graph.gfa",
            "02.unitig_graph",
            "direct-anchor graph or conservative skeleton seed graph",
        ),
        (
            "02.unitig_graph/graph.full.gfa",
            "02.unitig_graph",
            "less-filtered companion graph for debugging",
        ),
        (
            "03.read_junction_graph/graph.gfa",
            "03.read_junction_graph",
            "companion graph with read-junction links enabled",
        ),
        (
            "03.read_junction_graph/report.txt",
            "03.read_junction_graph",
            "read-junction graph report or skipped-step report",
        ),
        (
            "04.low_depth_bridge_rescue/report.txt",
            "04.low_depth_bridge_rescue",
            "low-depth bridge rescue status and evidence handoff report",
        ),
        (
            "04.low_depth_bridge_rescue/link_selection.report.txt",
            "04.low_depth_bridge_rescue",
            "low-depth bridge rescue link selection report",
        ),
        (
            "04.low_depth_bridge_rescue/selected_links.tsv",
            "04.low_depth_bridge_rescue",
            "bridge links selected during low-depth rescue",
        ),
        (
            "04.low_depth_bridge_rescue/pruned_links.tsv",
            "04.low_depth_bridge_rescue",
            "links pruned after low-depth bridge rescue",
        ),
        (
            "04.low_depth_bridge_rescue/focused_read_ids.txt",
            "04.low_depth_bridge_rescue",
            "reads used for focused low-depth bridge evaluation",
        ),
        (
            "05.skeleton_link_evidence/skeleton_segments.fasta",
            "05.skeleton_link_evidence",
            "unitig sequences used as skeleton remapping targets",
        ),
        (
            "05.skeleton_link_evidence/read_alignments.paf",
            "05.skeleton_link_evidence",
            "read-to-skeleton minimap2 alignments",
        ),
        (
            "05.skeleton_link_evidence/links.tsv",
            "05.skeleton_link_evidence",
            "skeleton link support table",
        ),
        (
            "05.skeleton_link_evidence/depth.tsv",
            "05.skeleton_link_evidence",
            "read depth on skeleton unitigs",
        ),
        (
            "05.skeleton_link_evidence/report.txt",
            "05.skeleton_link_evidence",
            "skeleton evidence report or skipped-step report",
        ),
        (
            "05.skeleton_link_evidence/graph.gfa",
            "05.skeleton_link_evidence",
            "linked graph produced directly from skeleton evidence",
        ),
        (
            "05.skeleton_link_evidence/linking.report.txt",
            "05.skeleton_link_evidence",
            "direct skeleton-linking report",
        ),
        (
            "06.repeat_aware_resolution/report.txt",
            "06.repeat_aware_resolution",
            "repeat-aware resolution, skeleton-link summary, and topology repair report",
        ),
        (
            "06.repeat_aware_resolution/node_repairs.tsv",
            "06.repeat_aware_resolution",
            "compact node-level summary of repeat-aware repairs and selected links",
        ),
        (
            "06.repeat_aware_resolution/links.tsv",
            "06.repeat_aware_resolution",
            "links retained after repeat-aware resolution",
        ),
        (
            "06.repeat_aware_resolution/depth.tsv",
            "06.repeat_aware_resolution",
            "depth retained after repeat-aware resolution",
        ),
        (
            "06.repeat_aware_resolution/graph.gfa",
            "06.repeat_aware_resolution",
            "repeat-aware resolved graph before final handoff",
        ),
        (
            "07.linked_graph/graph.gfa",
            "07.linked_graph",
            "linked graph selected from skeleton link evidence",
        ),
        (
            "07.linked_graph/links.tsv",
            "07.linked_graph",
            "links retained in the linked graph",
        ),
        (
            "07.linked_graph/depth.tsv",
            "07.linked_graph",
            "depth retained in the linked graph",
        ),
        (
            "07.linked_graph/report.txt",
            "07.linked_graph",
            "linked graph selection report",
        ),
        (
            "08.workflow_summary/report.txt",
            "08.workflow_summary",
            "anchor graph workflow report",
        ),
        (
            "08.workflow_summary/profile_parameters.tsv",
            "08.workflow_summary",
            "standard baseline versus active profile parameter report",
        ),
    ];

    let mut out = File::create(path)?;
    writeln!(out, "file\tstage\tdescription")?;
    for (relative, stage, description) in expected {
        if work_dir.join(relative).exists() {
            writeln!(out, "{relative}\t{stage}\t{description}")?;
        }
    }
    Ok(())
}

fn find_selected_graph(work_dir: &Path) -> Option<PathBuf> {
    for relative in [
        "07.linked_graph/graph.gfa",
        "06.repeat_aware_resolution/graph.gfa",
        "02.unitig_graph/graph.gfa",
        "graph.gfa",
    ] {
        let graph = work_dir.join(relative);
        if graph.exists() {
            return Some(graph);
        }
    }

    let mut subset_graphs = Vec::new();
    for entry in fs::read_dir(work_dir).ok()? {
        let entry = entry.ok()?;
        let name = entry.file_name();
        let name = name.to_string_lossy();
        let Some(percent) = name.strip_prefix("read_subset_") else {
            continue;
        };
        let mut graph = entry.path().join("07.linked_graph/graph.gfa");
        if !graph.exists() {
            graph = entry.path().join("06.repeat_aware_resolution/graph.gfa");
        }
        if !graph.exists() {
            graph = entry.path().join("02.unitig_graph/graph.gfa");
        }
        if !graph.exists() {
            graph = entry.path().join("graph.gfa");
        }
        if !graph.exists() {
            continue;
        }
        let rank = percent.parse::<usize>().unwrap_or(0);
        subset_graphs.push((rank, graph));
    }
    subset_graphs.sort_by(|a, b| a.0.cmp(&b.0));
    subset_graphs.pop().map(|(_, graph)| graph)
}

fn write_run_manifest(path: &Path, summaries: &[RunSummary]) -> Result<(), OrgraftError> {
    let mut out = File::create(path)?;
    writeln!(
        out,
        "organelle\tprofile\trepeat_aware_resolution\treads\twork_dir\tselected_graph\toutput_gfa\telapsed_seconds"
    )?;
    for summary in summaries {
        writeln!(
            out,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.3}",
            summary.organelle.as_str(),
            summary.profile.as_str(),
            summary.repeat_aware_resolution,
            summary.reads.display(),
            summary.work_dir.display(),
            summary.selected_graph.display(),
            summary.output_gfa.display(),
            summary.elapsed_seconds,
        )?;
    }
    Ok(())
}

fn write_algorithm_notes(path: &Path) -> Result<(), OrgraftError> {
    let text = r#"# OrgRAFT draft assembly algorithm notes

This command keeps one readable pipeline order around the validated
anchor graph assembly core.

## 01.input_reads

Use the recruited organelle FASTQ/FASTQ.GZ from `orgraft recruit`. OrgRAFT
records the original path and creates a local `reads.fastq.gz` symlink when the
platform supports symlinks.

## 02.anchor_graph_core

Run one organelle profile, not a batch of modes. `standard` is the baseline
parameter set; other profiles are expressed as standard plus the smallest
profile-specific differences:

- low: standard plus a lower shared graph coverage floor for small
  corrected-read datasets; mito low also records low-depth bridge rescue as
  Step 04.
- standard: organelle defaults; plastid uses the direct-anchor workflow, mito
  uses the skeleton-link workflow.
- `--stable`: standard plus repeat-aware resolution for unstable mitogenomes.
- high: standard plus deterministic read-subset candidates, selected after
  assembly by topology first and posterior `read_bases / graph_bases` near
  300x. Plastid high wraps the direct-anchor workflow; mito high wraps the
  skeleton-link workflow.

The assembly core then:

1. `01.anchor_walk_support`: read FASTQ/FASTA, find ordered syncmer
   anchors, convert reads into anchor walks, and write anchor/edge support.
2. `02.unitig_graph`: filter low-support anchors/edges, compress non-branching
   graph paths into unitigs, and write the primary graph. Full debug graph
   companions are kept only with `--keep-debug-files`.
3. `03.read_junction_graph`: for skeleton-link workflows, write the same graph
   with read-junction links enabled as rescue evidence.
4. `04.low_depth_bridge_rescue`: for mito low workflows, record the bridge
   rescue handoff from read-junction evidence into skeleton-link selection.
5. `05.skeleton_link_evidence`: for skeleton-link workflows, remap reads to
   skeleton unitigs and collect depth/link evidence.
6. `06.repeat_aware_resolution`: for repeat-aware mito workflows, resolve
   candidate links, copy choices, repeat expansions, and topology repairs.
7. `07.linked_graph`: write the unified linked-graph handoff, either copied
   from the direct-anchor graph or selected from skeleton evidence.
8. `08.workflow_summary`: write `report.txt` and `manifest.tsv` for every
   workflow; skipped steps are recorded in the report.

## 03.finalize_graph

Select `02.anchor_graph_core/07.linked_graph/graph.gfa` when present, otherwise
select `02.anchor_graph_core/02.unitig_graph/graph.gfa`. If a high profile
produced read-subset directories, select the topology-first candidate closest
to the target read depth. Publish only this selected graph into
`03.finalize_graph/graph.gfa`.

The assembly core may keep the same physical connection in both read
orientations, for example a link and its reverse-complement counterpart. The
finalize copy removes those reverse-complement duplicate L records by default
so topology checks count one physical connection once. Use
`--finalize-dedup-rc-links off` to publish the selected graph unchanged. Source
graphs under `02.anchor_graph_core` remain unchanged for debugging and
mito-stable analysis.
"#;
    fs::write(path, text)?;
    Ok(())
}

fn read_soft_paths_optional(path: &Path) -> Result<HashMap<String, PathBuf>, OrgraftError> {
    if !path.exists() {
        return Ok(HashMap::new());
    }
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

fn split_tool_line(line: &str) -> Option<(&str, &str)> {
    line.split_once('\t')
        .or_else(|| line.split_once(char::is_whitespace))
        .map(|(name, path)| (name.trim(), path.trim()))
        .filter(|(name, path)| !name.is_empty() && !path.is_empty())
}

fn strip_comment(line: &str) -> &str {
    line.split_once('#').map(|(value, _)| value).unwrap_or(line)
}

fn augmented_path(soft_paths: &HashMap<String, PathBuf>) -> String {
    let mut dirs = BTreeSet::new();
    for path in soft_paths.values() {
        if let Some(parent) = path.parent() {
            dirs.insert(parent.display().to_string());
        }
    }
    if let Some(path) = env::var_os("PATH") {
        for dir in env::split_paths(&path) {
            dirs.insert(dir.display().to_string());
        }
    }
    env::join_paths(dirs.iter().map(PathBuf::from))
        .ok()
        .and_then(|path| path.into_string().ok())
        .unwrap_or_default()
}

struct PathEnvGuard {
    previous: Option<OsString>,
}

impl PathEnvGuard {
    fn set_augmented(soft_paths: &HashMap<String, PathBuf>) -> Self {
        let previous = env::var_os("PATH");
        unsafe {
            env::set_var("PATH", augmented_path(soft_paths));
        }
        Self { previous }
    }
}

impl Drop for PathEnvGuard {
    fn drop(&mut self) {
        unsafe {
            match &self.previous {
                Some(value) => env::set_var("PATH", value),
                None => env::remove_var("PATH"),
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::time::{SystemTime, UNIX_EPOCH};

    #[test]
    fn empty_args_show_help_without_running() {
        assert!(run(&[]).is_ok());
    }

    #[test]
    fn profile_help_succeeds_without_required_inputs() {
        assert!(run(&["--profile-help".to_string()]).is_ok());
    }

    #[test]
    fn main_help_points_to_profile_help() {
        assert!(HELP.contains("--profile-help"));
        assert!(HELP.contains("--threads N               threads passed to the assembly core [8]"));
        assert!(HELP.contains("Outputs:"));
        assert!(HELP.contains(
            "Layout: OUT/ORGANELLE/{01.input_reads,02.anchor_graph_core,03.finalize_graph,logs}"
        ));
        assert!(!HELP.contains("<out-dir>/<organelle>/03.finalize_graph/graph.gfa"));
        assert!(!HELP.contains("Workflow frame:"));
        assert!(!HELP.contains("--stable                  for unstable mitogenomes"));
        assert!(!HELP.contains("Profiles are presets on one 01-08 workflow frame"));
    }

    #[test]
    fn profile_help_contains_profile_details() {
        assert!(PROFILE_HELP.contains("Profiles are presets on one 01-08 workflow frame"));
        assert!(PROFILE_HELP.contains("Workflow frame:"));
        assert!(PROFILE_HELP.contains("03.finalize_graph publishes the selected graph"));
        assert!(PROFILE_HELP.contains("--finalize-dedup-rc-links on|off"));
        assert!(PROFILE_HELP.contains("--image-reference-fasta FILE"));
        assert!(PROFILE_HELP.contains("--stable                  for unstable mitogenomes"));
        assert!(PROFILE_HELP.contains("--keep-debug-files"));
    }

    #[test]
    fn finalize_graph_removes_reverse_complement_duplicate_links() {
        let dir = test_dir("asm_finalize_rc_duplicates");
        fs::create_dir_all(&dir).unwrap();
        let selected = dir.join("selected.gfa");
        let finalized = dir.join("graph.gfa");
        fs::write(
            &selected,
            concat!(
                "H\tVN:Z:1.0\n",
                "S\tutg0\tAAAA\n",
                "S\tutg1\tCCCC\n",
                "S\tutg2\tGGGG\n",
                "L\tutg0\t+\tutg2\t+\t0M\tRC:i:838\n",
                "L\tutg0\t+\tutg2\t-\t0M\tRC:i:638\n",
                "L\tutg0\t-\tutg1\t+\t0M\tRC:i:672\n",
                "L\tutg0\t-\tutg1\t-\t0M\tRC:i:626\n",
                "L\tutg1\t+\tutg0\t+\t0M\tRC:i:626\n",
                "L\tutg1\t-\tutg0\t+\t0M\tRC:i:672\n",
                "L\tutg2\t+\tutg0\t-\t0M\tRC:i:638\n",
                "L\tutg2\t-\tutg0\t-\t0M\tRC:i:838\n",
            ),
        )
        .unwrap();

        let stats = copy_finalize_graph(&selected, &finalized, true).unwrap();
        let text = fs::read_to_string(&finalized).unwrap();
        let link_count = text.lines().filter(|line| line.starts_with("L\t")).count();

        assert_eq!(
            stats,
            FinalizeGraphStats {
                input_links: 8,
                output_links: 4,
                rc_duplicate_links_removed: 4,
            }
        );
        assert_eq!(link_count, 4);
        assert!(text.contains("L\tutg0\t+\tutg2\t+\t0M\tRC:i:838\n"));
        assert!(!text.contains("L\tutg2\t-\tutg0\t-\t0M\tRC:i:838\n"));
    }

    #[test]
    fn finalize_graph_can_keep_reverse_complement_duplicate_links() {
        let dir = test_dir("asm_finalize_keep_rc_duplicates");
        fs::create_dir_all(&dir).unwrap();
        let selected = dir.join("selected.gfa");
        let finalized = dir.join("graph.gfa");
        fs::write(
            &selected,
            concat!(
                "H\tVN:Z:1.0\n",
                "S\tutg0\tAAAA\n",
                "S\tutg1\tCCCC\n",
                "L\tutg0\t-\tutg1\t-\t0M\tRC:i:626\n",
                "L\tutg1\t+\tutg0\t+\t0M\tRC:i:626\n",
            ),
        )
        .unwrap();

        let stats = copy_finalize_graph(&selected, &finalized, false).unwrap();
        let text = fs::read_to_string(&finalized).unwrap();
        let link_count = text.lines().filter(|line| line.starts_with("L\t")).count();

        assert_eq!(
            stats,
            FinalizeGraphStats {
                input_links: 2,
                output_links: 2,
                rc_duplicate_links_removed: 0,
            }
        );
        assert_eq!(link_count, 2);
    }

    #[test]
    fn args_without_reads_fail_explicitly() {
        let err = AsmOptions::from_args(&["--threads".to_string(), "2".to_string()]).unwrap_err();
        assert!(err.to_string().contains("missing --reads FILE"));
    }

    #[test]
    fn reads_and_organelle_default_to_standard_profile() {
        let options = AsmOptions::from_args(&[
            "--reads".to_string(),
            "mito.fastq.gz".to_string(),
            "--organelle".to_string(),
            "mito".to_string(),
        ])
        .unwrap();

        assert_eq!(options.input.organelle, Organelle::Mito);
        assert_eq!(options.input.profile, AsmProfile::Standard);
        assert_eq!(options.threads, DEFAULT_THREADS);
        assert!(options.finalize_dedup_rc_links);
        assert_eq!(options.image_reference_fasta, None);
    }

    #[test]
    fn parses_shared_graph_coverage_override() {
        let options = AsmOptions::from_args(&[
            "--reads".to_string(),
            "mito.fastq.gz".to_string(),
            "--organelle".to_string(),
            "mito".to_string(),
            "--min-graph-coverage".to_string(),
            "12".to_string(),
            "--min-link-ratio".to_string(),
            "0.2".to_string(),
            "--subsets=3,5,10".to_string(),
        ])
        .unwrap();

        assert_eq!(options.min_graph_coverage, Some(12));
        assert_eq!(options.min_link_ratio, Some(0.2));
        assert_eq!(options.read_subsets, Some(vec![300, 500, 1000]));
    }

    #[test]
    fn parses_finalize_dedup_rc_links_off() {
        let options = AsmOptions::from_args(&[
            "--reads".to_string(),
            "plastid.fastq.gz".to_string(),
            "--organelle".to_string(),
            "plastid".to_string(),
            "--finalize-dedup-rc-links".to_string(),
            "off".to_string(),
        ])
        .unwrap();

        assert!(!options.finalize_dedup_rc_links);
    }

    #[test]
    fn stable_is_mito_additional_step() {
        let options = AsmOptions::from_args(&[
            "--reads".to_string(),
            "mito.fastq.gz".to_string(),
            "--organelle".to_string(),
            "mito".to_string(),
            "--stable".to_string(),
        ])
        .unwrap();

        assert_eq!(options.input.profile, AsmProfile::Standard);
        assert!(options.input.repeat_aware_resolution);
    }

    #[test]
    fn parses_keep_debug_files() {
        let options = AsmOptions::from_args(&[
            "--reads".to_string(),
            "mito.fastq.gz".to_string(),
            "--organelle".to_string(),
            "mito".to_string(),
            "--keep-debug-files".to_string(),
        ])
        .unwrap();

        assert!(options.keep_debug_files);
    }

    #[test]
    fn parses_image_reference_fasta() {
        let options = AsmOptions::from_args(&[
            "--reads".to_string(),
            "mito.fastq.gz".to_string(),
            "--organelle".to_string(),
            "mito".to_string(),
            "--image-reference-fasta".to_string(),
            "refs/mito.fa".to_string(),
        ])
        .unwrap();

        assert_eq!(
            options.image_reference_fasta,
            Some(PathBuf::from("refs/mito.fa"))
        );
    }

    #[test]
    fn stable_rejects_plastid() {
        let err = AsmOptions::from_args(&[
            "--reads".to_string(),
            "plastid.fastq.gz".to_string(),
            "--organelle".to_string(),
            "plastid".to_string(),
            "--stable".to_string(),
        ])
        .unwrap_err();

        assert!(err
            .to_string()
            .contains("--stable is only available with --organelle mito"));
    }

    fn test_dir(name: &str) -> PathBuf {
        let nanos = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        std::env::temp_dir().join(format!("orgraft_{name}_{nanos}"))
    }
}
