use std::collections::{BTreeMap, BTreeSet};
use std::fmt::Write as FmtWrite;
use std::fs::{self, File};
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};
use std::process::Command;

use crate::commands;
use crate::error::OrgraftError;
use crate::sv_repair::{
    detect_repeat_pairing_mismatch, repair_sv_subgroup, select_sv_subgroup_spec,
    write_repeat_pairing_mismatch_report, RepeatPairingMismatch, SvRepairRequest,
};
use crate::topology::{analyze_gfa, nodes_tsv, summary_tsv, TopologyReport};

const WORKFLOW_HARD_MAX_ROUNDS: usize = 10;

const HELP: &str = r#"orgraft workflow

Coordinate OrgRAFT project folders, manual checkpoints, and validation rounds.

Usage:
  orgraft workflow <command> [options]

Commands:
  template            print workflow TOML template
  init                write workflow TOML template

Template/init options:
    --sample NAME       project sample [inferred from output folder]
    --results-dir DIR   output root [results_workflow]
    --soft-paths FILE   tool paths [soft_paths.txt]

Workflow commands:
  plan                write runnable commands from workflow config
  run-script          generate command script, then execute it with bash
  run                 execute automatic workflow inside orgraft
  runtime-summary     write results_dir/runtime_summary.md

Checkpoint/test commands:
  checkpoint1         check topology and write checked_draft.gfa when simple
  checkpoint2         localize/check SV and prepare the next correction round
  postrepair-run      run repair and replay an explicitly configured human decision
  postrepair-check    verify the repair handoff and print its validated FASTA
  refine-bundle       verify rebuild inputs and prepare one organelle component bundle
  refine-check        report refine result and stop when review is required
  correct             apply pos/ref/alt table to one FASTA record
  test-correction     smoke test FASTA correction with explicit inputs
  test-fake-validate  simulate validate failure from swapped polish_aln FASTA

Common options:
  --config FILE       workflow config [orgraft.workflow.toml]
  --case NAME         select one workflow.case section
  --round N           checkpoint2 validation/correction round [1]
  --source-only       with postrepair-check, validate a declared read-only source round
  --sv-subgroup SPEC  manually repair one group_name:old_index at checkpoint2
  --force             overwrite generated checkpoint/correction outputs; with plan,
                      regenerate the two entry scripts without running them

The workflow layer is intentionally thin: resolve/polish/rebuild own their core
algorithms; workflow owns config parsing, project layout, run-command generation,
manual checkpoint status, max_rounds ordinary validation/correction rounds, and
extra SV-repair rounds up to a hard total of 10.
"#;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct WorkflowStep {
    pub command: &'static str,
    pub responsibility: &'static str,
}

pub const DEFAULT_WORKFLOW: &[WorkflowStep] = &[
    WorkflowStep {
        command: "workflow_config",
        responsibility: "define sample, output directory, software paths, and command inputs",
    },
    WorkflowStep {
        command: "recruit",
        responsibility: "enrich organelle HiFi reads",
    },
    WorkflowStep {
        command: "asm",
        responsibility: "build a conservative draft graph",
    },
    WorkflowStep {
        command: "resolve",
        responsibility: "prepare checked draft GFA for high-quality graph generation",
    },
    WorkflowStep {
        command: "polish",
        responsibility: "polish linearized subgraph FASTA and evaluate variants",
    },
    WorkflowStep {
        command: "postrepair",
        responsibility: "review residual errors from the completed validation round",
    },
    WorkflowStep {
        command: "rebuild",
        responsibility: "rebuild final verified graph and compact reports",
    },
    WorkflowStep {
        command: "refine",
        responsibility: "refine repeat representation and component connectivity",
    },
];

#[derive(Debug, Clone)]
struct WorkflowConfig {
    config_path: PathBuf,
    sample: String,
    results_dir: PathBuf,
    soft_paths: PathBuf,
    mode: String,
    max_rounds: usize,
    threads: usize,
    force: bool,
    auto_sv_correction: bool,
    auto_snv_indel_correction: bool,
    recruit: WorkflowRecruitConfig,
    asm: WorkflowAsmConfig,
    resolve: WorkflowResolveConfig,
    polish: WorkflowPolishConfig,
    postrepair: Option<WorkflowPostrepairConfig>,
    refine: Option<WorkflowRefineConfig>,
    rebuild: WorkflowRebuildConfig,
    topology_simple_allowed_classes: BTreeSet<String>,
    cases: Vec<WorkflowCase>,
}

#[derive(Debug, Clone)]
struct WorkflowCase {
    enabled: bool,
    name: String,
    sample: String,
    organelle: String,
    subgraph: String,
    draft_graph: PathBuf,
    unitig_graph: Option<PathBuf>,
    checked_draft_gfa: PathBuf,
    reference: Option<PathBuf>,
    pre_rotated_reference: Option<PathBuf>,
    reads: PathBuf,
    asm_reads: Option<PathBuf>,
    resolve_out_dir: PathBuf,
    polish_out_dir: PathBuf,
    postrepair_out_dir: Option<PathBuf>,
    postrepair_source_dir: Option<PathBuf>,
    resolve_details: Option<PathBuf>,
    first_polish_report: Option<PathBuf>,
    rebuild_out_dir: PathBuf,
    rebuild_edited_gfa: Option<PathBuf>,
    rebuild_polished_fasta: Option<PathBuf>,
    image_reference_fasta: Option<PathBuf>,
    workflow_dir: PathBuf,
    linearized_fasta: Option<PathBuf>,
    polish_reference: Option<PathBuf>,
    sv_correction_subgroup: Option<String>,
    refine_downgrade_representation: commands::refine::DowngradeRepresentation,
}

#[derive(Debug, Clone)]
struct WorkflowRecruitConfig {
    enabled: bool,
    reads: Option<PathBuf>,
    out_dir: PathBuf,
    threads: Option<usize>,
    baits: Vec<(String, PathBuf)>,
    prefix: Option<String>,
    bait_format: Option<String>,
    gfa_split: Option<String>,
    rename_bait: bool,
    write_id_map: bool,
    split_output: Option<String>,
    gzip_output: Option<bool>,
    minimap2: Option<String>,
    align_mode: Option<String>,
    platform: Option<String>,
    preset: Option<String>,
    min_mapq: Option<u8>,
    min_aln_len: Option<u64>,
    sam: Option<PathBuf>,
    max_reads: Vec<String>,
    random_seed: Option<u64>,
    write_sampled_ids: bool,
    read_stats: Option<String>,
    write_read_classification: bool,
    write_bait_partitions: bool,
    gzip_tool: Option<String>,
    mode: Option<String>,
    iterations: Option<usize>,
    extra_args: Vec<String>,
}

#[derive(Debug, Clone)]
struct WorkflowAsmConfig {
    enabled: bool,
    out_dir: PathBuf,
    threads: Option<usize>,
    profile: Option<String>,
    stable: bool,
    min_graph_coverage: Option<u32>,
    min_branch_ratio: Option<f64>,
    min_tip_len: Option<usize>,
    min_link_support: Option<u32>,
    min_link_ratio: Option<f64>,
    subsets: Option<String>,
    image_reference_fasta: Option<PathBuf>,
    keep_debug_files: bool,
    extra_args: Vec<String>,
}

#[derive(Debug, Clone)]
struct WorkflowResolveConfig {
    out_dir: PathBuf,
    gfa_editor_mode: Option<String>,
    max_states: Option<usize>,
    max_candidates: Option<usize>,
    extra_args: Vec<String>,
}

#[derive(Debug, Clone)]
struct WorkflowPolishConfig {
    out_dir: PathBuf,
    homorepair: bool,
    homorepair_frq: f64,
    threads: Option<usize>,
    per_read_variant_calls: Option<bool>,
    snv_indel_overlap_policy: Option<String>,
    plot_range: Option<String>,
    plot_dpi: Option<usize>,
    plot_output_format: Option<String>,
    coverage_plot_rasterize: Option<bool>,
    snv_indel_plot_rasterize: Option<bool>,
    sv_plot_highlight_subgroups: Option<String>,
    sv_plot_highlight_read_ids: Option<PathBuf>,
    sv_plot_highlight_min_fraction: Option<f64>,
    sv_plot_highlight_min_reads: Option<usize>,
    snv_indel_plot_low_confidence: Option<String>,
    snv_indel_plot_low_min_reads: Option<usize>,
    snv_indel_plot_low_min_fraction: Option<f64>,
    snv_indel_plot_high_risk_fraction: Option<f64>,
    extra_args: Vec<String>,
}

#[derive(Debug, Clone)]
struct WorkflowPostrepairConfig {
    enabled: bool,
    out_dir: PathBuf,
    threads: Option<usize>,
    frequency: f64,
    review_files: BTreeMap<String, PathBuf>,
}

#[derive(Debug, Clone)]
struct WorkflowRefineConfig {
    enabled: bool,
    out_dir: PathBuf,
    threads: Option<usize>,
    bundle: Option<PathBuf>,
    organelle_reads: Option<PathBuf>,
    adopt: BTreeMap<String, PathBuf>,
    intra_measurement: BTreeMap<String, commands::refine::IntraMeasurement>,
    duplicate_template_full_span: BTreeMap<String, bool>,
    linear_chain_flanks: BTreeMap<String, bool>,
    graph_full_family: BTreeMap<String, bool>,
    graph_measurement: BTreeMap<String, commands::refine::GraphMeasurement>,
}

#[derive(Debug, Clone)]
struct WorkflowRebuildConfig {
    enabled: bool,
    out_dir: PathBuf,
    threads: Option<usize>,
    edited_gfa: Option<PathBuf>,
    polished_fasta: Option<PathBuf>,
    image_reference_fasta: Option<PathBuf>,
    merged_gfa_template: Option<PathBuf>,
    minimap2: Option<PathBuf>,
    blastn: Option<PathBuf>,
    keep_debug: bool,
    extra_args: Vec<String>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
enum Checkpoint1Status {
    Checked,
    ManualRequired,
}

#[derive(Debug, Clone, PartialEq, Eq)]
enum Checkpoint2Status {
    Complete,
    NextRoundReady,
    RebuildReady,
    ManualRequired,
}

#[derive(Debug, Clone)]
struct VariantEdit {
    pos: usize,
    reference: String,
    alternate: String,
}

#[derive(Debug, Clone)]
struct RepeatPairingMarker {
    state: String,
    repeat_node: String,
    source_round: usize,
    mismatch_report: PathBuf,
    rebuild_fasta: Option<PathBuf>,
    validation_round: Option<usize>,
}

pub fn run(args: &[String]) -> Result<(), OrgraftError> {
    if args.is_empty() || args.iter().any(|arg| arg == "-h" || arg == "--help") {
        println!("{HELP}");
        return Ok(());
    }

    match args.first().map(String::as_str) {
        Some("template") => {
            let options = TemplateOptions::from_args(&args[1..])?;
            print!("{}", workflow_config_template(&options));
            Ok(())
        }
        Some("init") => init_workflow_config(&args[1..]),
        Some("plan") => run_plan(&args[1..]),
        Some("run-script") | Some("script") => run_generated_script(&args[1..]),
        Some("runtime-summary") | Some("summary") => run_runtime_summary(&args[1..]),
        Some("checkpoint1") | Some("check-topology") => run_checkpoint1(&args[1..]),
        Some("checkpoint2") | Some("check-validation") => run_checkpoint2(&args[1..]),
        Some("postrepair-check") | Some("homorepair-check") => run_homorepair_check(&args[1..]),
        Some("postrepair-run") => run_configured_postrepair(&args[1..]),
        Some("refine-bundle") => run_refine_bundle(&args[1..]),
        Some("refine-check") => run_refine_check(&args[1..]),
        Some("correct") => run_correct(&args[1..]),
        Some("test-correction") => run_test_correction(&args[1..]),
        Some("test-fake-validate") => run_test_fake_validate(&args[1..]),
        Some("run") => run_automatic_workflow(&args[1..]),
        Some(other) => Err(OrgraftError::InvalidArgument(format!(
            "unknown workflow action `{other}`"
        ))),
        None => {
            println!("{HELP}");
            Ok(())
        }
    }
}

#[derive(Debug, Clone)]
struct TemplateOptions {
    sample: String,
    results_dir: String,
    soft_paths: String,
}

impl TemplateOptions {
    fn from_args(args: &[String]) -> Result<Self, OrgraftError> {
        if option_value(args, "--command-mode")?.is_some() {
            return Err(OrgraftError::InvalidArgument(
                "workflow templates no longer support --command-mode; scripts are always expanded"
                    .to_string(),
            ));
        }
        let sample = option_value(args, "--sample")?.map(ToString::to_string)
            .unwrap_or_else(|| {
                let output = option_value(args, "--out").ok().flatten().map(PathBuf::from);
                let result = option_value(args, "--results-dir").ok().flatten().map(PathBuf::from);
                output.as_deref().and_then(Path::parent).filter(|p| !p.as_os_str().is_empty())
                    .or(result.as_deref()).and_then(Path::file_name).and_then(|s| s.to_str())
                    .filter(|s| *s != ".").map(ToString::to_string)
                    .unwrap_or_else(|| "sample".into())
            });
        Ok(Self {
            sample,
            results_dir: option_value(args, "--results-dir")?.unwrap_or(".").to_string(),
            soft_paths: option_value(args, "--soft-paths")?
                .unwrap_or("soft_paths.txt")
                .to_string(),
        })
    }
}

fn init_workflow_config(args: &[String]) -> Result<(), OrgraftError> {
    let options = TemplateOptions::from_args(args)?;
    let out = option_value(args, "--out")?.unwrap_or("orgraft.workflow.toml");
    let force = has_flag(args, "--force");
    let out_path = Path::new(out);

    if out_path.exists() && !force {
        return Err(OrgraftError::InvalidArgument(format!(
            "{} already exists; pass --force to overwrite",
            out_path.display()
        )));
    }

    write_with_parent(out_path, &workflow_config_template(&options))?;
    println!("Wrote {}", out_path.display());
    Ok(())
}

fn workflow_config_template(options: &TemplateOptions) -> String {
    let sample = toml_string(&options.sample);
    let results_dir = toml_string(&options.results_dir);
    let soft_paths = toml_string(&options.soft_paths);

    format!(
        r#"# OrgRAFT workflow configuration.
# This file drives both stepwise batch mode and automatic checkpoint mode.

[project]
sample = {sample}
results_dir = {results_dir}

[software]
soft_paths = {soft_paths}

[workflow]
# stepwise: generate/check files, then let an external batch runner call commands.
# automatic: `orgraft workflow run` may execute resolve/polish rounds directly.
mode = "stepwise"
max_rounds = 3
# An accepted SV repair adds one actual validation round without consuming the
# max_rounds budget. Total actual rounds are always capped at 10.
auto_sv_correction = true
# Fallback when a command-specific threads value is omitted.
threads = 64
force = false
# When true, checkpoint2 applies high-confidence SNV/InDel rows to make the next
# polish input. Set false when `sv_snv_indel_summary.tsv status=pass` should be
# accepted as the round result and correction is tested separately.
auto_snv_indel_correction = true

# Default policy: higher-order complex nodes and self-associated branches require
# manual graph editing at checkpoint 1.
topology_simple_allowed_classes = "0-0,0-1/1-0,1-1,1-2/2-1,2-2"

[commands.recruit]
enabled = true
raw_reads = "/path/to/raw_hifi.fastq.gz"
baits = "mito=/path/to/mito.fasta,plastid=/path/to/plastid.fasta"
out_dir = "${{results_dir}}/01.recruit"
threads = 16
# platform = "HiFi"
# bait_format = "auto"
# gzip_output = true
# max_reads = "all,20000"
# extra_args = "--write-read-classification --write-sampled-ids"

[commands.asm]
enabled = true
out_dir = "${{results_dir}}/02.draft_asm"
threads = 8
# profile = "standard"
# stable = true
# min_graph_coverage = 18
# min_branch_ratio = 0.30
# min_tip_len = 3000
# min_link_support = 20
# min_link_ratio = 0.05
# subsets = "3,5,10"
# image_reference_fasta = "/path/to/reference.fasta"  # workflow default: matching case reference/bait
# keep_debug_files = true

[commands.resolve]
out_dir = "${{results_dir}}/03.resolve_gfa"
# gfa_editor_mode = "rust"
# max_states = 5000
# max_candidates = 100
# extra_args = ""

[commands.polish]
out_dir = "${{results_dir}}/04.polish"
threads = 64
# per_read_variant_calls = true
# snv_indel_overlap_policy = "mark-overlap"
# plot_range = "1-50000"
# plot_dpi = 300
# plot_output_format = "png"
# coverage_plot_rasterize = true
# snv_indel_plot_rasterize = true
# sv_plot_highlight_subgroups = "subgraph_001:0"
# sv_plot_highlight_read_ids = "/path/to/read_ids.txt"
# sv_plot_highlight_min_fraction = 0.005
# sv_plot_highlight_min_reads = 10
# snv_indel_plot_low_confidence = "non-high"
# snv_indel_plot_low_min_reads = 3
# snv_indel_plot_low_min_fraction = 0
# snv_indel_plot_high_risk_fraction = 0.5
# extra_args = ""

[commands.postrepair]
enabled = true
out_dir = "${{results_dir}}/05.postrepair"
frequency = 0.2
threads = 16
# Source evidence remains in 04.polish. Changed sequences receive a new
# validate-only round in 05.postrepair; unresolved decisions stop the handoff.

[commands.rebuild]
enabled = true
out_dir = "${{results_dir}}/06.rebuild"
threads = 16
# The generated workflow script passes the final checkpoint2 polished FASTA to
# rebuild automatically, together with the case reads used for real repeat-path
# support. Set edited_gfa/polished_fasta only for standalone rebuild overrides.
# image_reference_fasta enables optional PDF/SVG export through GFA_Editor; export
# failures are recorded in the rebuild run report and do not fail core GFA/FASTA output.
# merged_gfa_template = "${{results_dir}}/03.resolve_gfa/mito/graph/merged_unresolved.gfa"
# keep_debug = true

[commands.refine]
enabled = true
out_dir = "${{results_dir}}/07.refine"
threads = 16
# For split component read pools, set organelle_reads to the full organelle pool.
# organelle_reads = "${{results_dir}}/01.recruit/${{organelle}}.fastq.gz"
# Without this value, all enabled cases for the organelle must use the same reads.
# An optional sample_bundle supplies a reviewed multi-component input.
# Otherwise workflow prepares a bundle from the matching rebuild output.
# sample_bundle = "/path/to/components.tsv"

# Use enabled = false to keep a case in this config while excluding it from the
# generated master workflow script. Explicit --case NAME can still generate/run it.

[workflow.case.mito_subgraph_001]
enabled = true
name = "mito_subgraph_001"
organelle = "mito"
subgraph = "subgraph_001"

# Numbered output handoff paths.
# 01.recruit writes reads -> 02.draft_asm builds draft graph -> workflow
# records checkpoint files -> 03.resolve_gfa writes resolved FASTA ->
# 04.polish validates -> 05.postrepair reviews residuals -> 06.rebuild writes
# graph/sequence products -> 07.refine reviews repeat representation.
workflow_dir = "${{results_dir}}/workflow/${{organelle}}/${{subgraph}}"

# Checkpoint 1 consumes the draft graph and writes checked_draft_gfa only when
# topology is simple and all GFA links reference declared segment records.
draft_graph = "${{results_dir}}/02.draft_asm/${{organelle}}/03.finalize_graph/graph.gfa"
# Optional override used only by SV check to project breakpoints back to the
# original unitig graph; the default is derived from commands.asm.out_dir.
# unitig_graph = "${{results_dir}}/02.draft_asm/${{organelle}}/02.anchor_graph_core/02.unitig_graph/graph.gfa"
checked_draft_gfa = "${{results_dir}}/workflow/${{organelle}}/${{subgraph}}/checkpoint_1/checked_draft.gfa"

# Resolve uses the matching FASTA from commands.recruit.baits when reference is omitted.
# reference = "/path/to/mito.fasta"
# resolve_out_dir defaults to commands.resolve.out_dir.
# resolve_out_dir = "${{results_dir}}/03.resolve_gfa"

# To start directly at 05 from an existing native round, keep recruit/asm disabled.
# postrepair_source_dir = "/old-run/04.polish/${{organelle}}/${{subgraph}}/round_2"
# resolve_details = "/old-run/03.resolve_gfa/${{organelle}}/logs/resolve_details.tsv"
# first_polish_report = "/old-run/04.polish/${{organelle}}/${{subgraph}}/round_1/logs/report.tsv"

# Explicit polish handoff inputs from 01.recruit and 03.resolve_gfa.
reads = "${{results_dir}}/01.recruit/${{organelle}}.fastq.gz"
# linearized_fasta defaults to commands.resolve.out_dir/${{organelle}}/fasta/resolved_subgraphs.fasta.
# linearized_fasta = "${{results_dir}}/03.resolve_gfa/${{organelle}}/fasta/resolved_subgraphs.fasta"
# polish_reference defaults to commands.resolve.out_dir/${{organelle}}/fasta/rotated_reference.fasta.
# polish_reference = "${{results_dir}}/03.resolve_gfa/${{organelle}}/fasta/rotated_reference.fasta"
# polish_out_dir defaults to commands.polish.out_dir.
# polish_out_dir = "${{results_dir}}/04.polish"
# Optional manual override for checkpoint2 SV repair. The automatic mode only
# selects `possible_reference_sv_error` rows whose local type_1 support is weak.
# sv_correction_subgroup = "type_3_subtype_rep_rep_NA:4"

rebuild_out_dir = "${{results_dir}}/06.rebuild/${{organelle}}"
rebuild_edited_gfa = "${{results_dir}}/03.resolve_gfa/${{organelle}}/graph/merged_unresolved_${{subgraph}}.gfa"
# checkpoint2 selects the final polished FASTA automatically; set this only for manual override.
# rebuild_polished_fasta = "${{results_dir}}/workflow/${{organelle}}/${{subgraph}}/checkpoint_2/round_1/polish_aln_v2.fasta"
image_reference_fasta = "${{results_dir}}/03.resolve_gfa/${{organelle}}/fasta/rotated_reference.fasta"

[workflow.case.plastid_subgraph_001]
enabled = true
name = "plastid_subgraph_001"
organelle = "plastid"
subgraph = "subgraph_001"

workflow_dir = "${{results_dir}}/workflow/${{organelle}}/${{subgraph}}"
draft_graph = "${{results_dir}}/02.draft_asm/${{organelle}}/03.finalize_graph/graph.gfa"
checked_draft_gfa = "${{results_dir}}/workflow/${{organelle}}/${{subgraph}}/checkpoint_1/checked_draft.gfa"
reads = "${{results_dir}}/01.recruit/${{organelle}}.fastq.gz"
# linearized_fasta = "${{results_dir}}/03.resolve_gfa/${{organelle}}/fasta/resolved_subgraphs.fasta"
# polish_reference = "${{results_dir}}/03.resolve_gfa/${{organelle}}/fasta/rotated_reference.fasta"
rebuild_out_dir = "${{results_dir}}/06.rebuild/${{organelle}}"
rebuild_edited_gfa = "${{results_dir}}/03.resolve_gfa/${{organelle}}/graph/merged_unresolved_${{subgraph}}.gfa"
image_reference_fasta = "${{results_dir}}/03.resolve_gfa/${{organelle}}/fasta/rotated_reference.fasta"
"#
    )
}

fn run_plan(args: &[String]) -> Result<(), OrgraftError> {
    let config = load_config_from_args(args)?;
    let out = write_plan_from_args(&config, args)?;
    println!("Wrote {}", out.display());
    Ok(())
}

fn run_generated_script(args: &[String]) -> Result<(), OrgraftError> {
    let config = load_config_from_args(args)?;
    let out = write_plan_from_args(&config, args)?;
    println!("Wrote {}", out.display());

    let mut command = Command::new("bash");
    command.arg(&out);
    if std::env::var_os("ORGRAFT_BIN").is_none() {
        if let Ok(current_exe) = std::env::current_exe() {
            command.env("ORGRAFT_BIN", current_exe);
        }
    }
    let status = command.status()?;
    if status.success() {
        Ok(())
    } else {
        Err(OrgraftError::InvalidArgument(format!(
            "{} exited with {status}",
            out.display()
        )))
    }
}

fn run_runtime_summary(args: &[String]) -> Result<(), OrgraftError> {
    let config = load_config_from_args(args)?;
    let out = option_value(args, "--out")?
        .map(PathBuf::from)
        .unwrap_or_else(|| config.results_dir.join("runtime_summary.md"));
    let force = config.force || has_flag(args, "--force");
    if enabled_cases(&config).iter().any(|case| case.postrepair_source_dir.is_some())
        && !workflow_future_path(&out)?.starts_with(&workflow_future_path(&config.results_dir)?) {
        return Err(OrgraftError::InvalidArgument("runtime summary must remain under results_dir".into()));
    }
    if out.exists() && !force {
        return Err(OrgraftError::InvalidArgument(format!(
            "{} already exists; pass --force to overwrite",
            out.display()
        )));
    }

    let summary = workflow_runtime_summary_markdown(&config)?;
    write_with_parent(&out, &summary)?;
    println!("Wrote {}", out.display());
    Ok(())
}

fn write_plan_from_args(config: &WorkflowConfig, args: &[String]) -> Result<PathBuf, OrgraftError> {
    declared_postrepair_reviews(config)?;
    let selected = option_value(args, "--case")?.is_some();
    let cases = if selected { vec![select_case(config, args)?] } else { enabled_cases(config) };
    let sample_name: String = config.sample.chars().map(|c| if c.is_ascii_alphanumeric() || "-_.".contains(c) { c } else { '_' }).collect();
    let out = option_value(args, "--out")?.map(PathBuf::from)
        .unwrap_or_else(|| config.results_dir.join(format!("run_{sample_name}.sh")));
    let manual = out.parent().unwrap_or(Path::new(".")).join(format!("run_{sample_name}_manual.sh"));
    if cases.iter().any(|case| case.postrepair_source_dir.is_some()) {
        let root = workflow_future_path(&config.results_dir)?;
        if !workflow_future_path(&out)?.starts_with(&root) || !workflow_future_path(&manual)?.starts_with(&root) {
            return Err(OrgraftError::InvalidArgument("postrepair source workflow scripts must be written under results_dir".into()));
        }
    }
    if out == manual {
        return Err(OrgraftError::InvalidArgument("master and manual script paths must differ".into()));
    }
    let script = if selected || cases.len() == 1 {
        let case = cases.first().ok_or_else(|| OrgraftError::InvalidArgument("no enabled workflow cases".into()))?;
        let mut script = render_plan_script(config, case, None, false, true, None)?;
        write_refine_stage(&mut script, config, &case.organelle)?;
        script
    } else {
        render_all_cases_plan_script(config)?
    };
    let manual_script = render_manual_plan_script(config, &cases)?;
    // Check both destinations before creating either entrypoint. A retry may
    // reuse identical scripts, but cannot replace a previous plan silently.
    let replace_scripts = has_flag(args, "--force");
    check_plan_destination(&out, &script, replace_scripts)?;
    check_plan_destination(&manual, &manual_script, replace_scripts)?;
    write_new_plan_script(&out, &script, replace_scripts)?;
    write_new_plan_script(&manual, &manual_script, replace_scripts)?;
    Ok(out)
}

fn run_checkpoint1(args: &[String]) -> Result<(), OrgraftError> {
    let config = load_config_from_args(args)?;
    let case = select_case(&config, args)?;
    let force = config.force || has_flag(args, "--force");
    let status = checkpoint1_impl(&config, case, force)?;
    match status {
        Checkpoint1Status::Checked => {
            println!("checkpoint1 checked: {}", case.checked_draft_gfa.display());
        }
        Checkpoint1Status::ManualRequired => {
            println!(
                "checkpoint1 manual_required: edit {} and then provide checked_draft_gfa",
                case.workflow_dir
                    .join("checkpoint_1")
                    .join("manual_edit_required.gfa")
                    .display()
            );
        }
    }
    Ok(())
}

fn run_checkpoint2(args: &[String]) -> Result<(), OrgraftError> {
    let config = load_config_from_args(args)?;
    let case = select_case(&config, args)?;
    let round = parse_round(args)?;
    let force = config.force || has_flag(args, "--force");
    let run_next = has_flag(args, "--run-next");
    let manual_sv_subgroup = option_value(args, "--sv-subgroup")?;
    let status = checkpoint2_impl(&config, case, round, force, run_next, manual_sv_subgroup)?;
    match status {
        Checkpoint2Status::Complete => println!("checkpoint2 complete at round {round}"),
        Checkpoint2Status::NextRoundReady => println!(
            "checkpoint2 next_round_ready: {}",
            corrected_fasta_path(case, round).display()
        ),
        Checkpoint2Status::RebuildReady => println!(
            "checkpoint2 rebuild_ready: constrain repeat during rebuild, then validate {}",
            rebuild_fasta_path(case).display()
        ),
        Checkpoint2Status::ManualRequired => {
            println!("checkpoint2 manual_required at round {round}");
        }
    }
    Ok(())
}

fn run_homorepair_check(args: &[String]) -> Result<(), OrgraftError> {
    let config = load_config_from_args(args)?;
    let case = select_case(&config, args)?;
    if has_flag(args, "--source-only") {
        if has_flag(args, "--rebuild-args") { return Err(OrgraftError::InvalidArgument("--source-only cannot be combined with --rebuild-args".into())); }
        let round = postrepair_source_round(case)?.ok_or_else(|| OrgraftError::InvalidArgument("--source-only requires postrepair_source_dir".into()))?;
        if option_value(args, "--round")?.is_some() && parse_round(args)? != round {
            return Err(OrgraftError::InvalidArgument("source check round differs from postrepair_source_dir".into()));
        }
        println!("{}", repair_source_fasta(case, round).display());
    } else if args.iter().any(|arg| arg == "--rebuild-args") {
        let polished = option_value(args, "--polished-fasta")?.ok_or_else(||
            OrgraftError::InvalidArgument("--rebuild-args requires --polished-fasta FILE".into()))?;
        for arg in postrepair_rebuild_args(&config, case, Path::new(polished))? { println!("{arg}"); }
    } else {
        let round = if option_value(args, "--round")?.is_some() { parse_round(args)? } else { postrepair_source_round(case)?.unwrap_or(1) };
        println!("{}", selected_repair_fasta(&config, case, round)?.display());
    }
    Ok(())
}

// Only the explicit stages require the recorded resolve-path source contract.
// The legacy polish homorepair setting does not change the rebuild source mode.
fn rebuild_uses_recorded_source(config: &WorkflowConfig) -> bool {
    config.postrepair.as_ref().is_some_and(|repair| repair.enabled)
        || config.refine.as_ref().is_some_and(|refine| refine.enabled)
}

fn rebuild_source_gfa(config: &WorkflowConfig, case: &WorkflowCase) -> PathBuf {
    case.rebuild_edited_gfa.as_ref().or(config.rebuild.edited_gfa.as_ref()).cloned()
        .unwrap_or_else(|| if rebuild_uses_recorded_source(config) {
            case.resolve_out_dir.join(&case.organelle).join("graph").join(format!("merged_unresolved_{}.gfa", case.subgraph))
        } else { case.checked_draft_gfa.clone() })
}

fn resolve_details_path(case: &WorkflowCase) -> PathBuf {
    case.resolve_details.clone().unwrap_or_else(|| case.resolve_out_dir.join(&case.organelle).join("logs/resolve_details.tsv"))
}

// A supplied ordinary round is a read-only input, not a checkpoint in this run.
fn postrepair_source_round(case: &WorkflowCase) -> Result<Option<usize>, OrgraftError> {
    let Some(source) = &case.postrepair_source_dir else {
        return Ok(None);
    };
    let round = source
        .file_name()
        .and_then(|name| name.to_str())
        .and_then(|name| name.strip_prefix("round_"))
        .and_then(|number| number.parse::<usize>().ok())
        .filter(|number| *number > 0 && *number <= WORKFLOW_HARD_MAX_ROUNDS)
        .ok_or_else(|| {
            OrgraftError::InvalidArgument(
                "postrepair_source_dir must name an ordinary round_N directory (1–10)".into(),
            )
        })?;
    if source.file_name().and_then(|name| name.to_str()) != Some(format!("round_{round}").as_str())
    {
        return Err(OrgraftError::InvalidArgument(
            "postrepair_source_dir round name is not canonical".into(),
        ));
    }
    Ok(Some(round))
}

fn repair_source_dir(case: &WorkflowCase, round: usize) -> PathBuf {
    case.postrepair_source_dir
        .clone()
        .unwrap_or_else(|| polish_round_dir(case, round))
}

fn repair_source_fasta(case: &WorkflowCase, round: usize) -> PathBuf {
    let source = repair_source_dir(case, round);
    if round == 1 {
        source.join("02.polish/polished_aln.fasta")
    } else {
        source.join(format!("01.inputs/linear_subgraph.round_{round}.fasta"))
    }
}

fn first_polish_report_path(case: &WorkflowCase) -> Option<PathBuf> {
    case.first_polish_report.clone().or_else(|| {
        case.postrepair_source_dir
            .as_ref()
            .and_then(|source| source.parent())
            .map(|component| component.join("round_1/logs/report.tsv"))
    })
}

// Canonicalize existing prefixes before interpreting '..', including symlinks.
fn workflow_future_path(path: &Path) -> Result<PathBuf, OrgraftError> {
    use std::path::Component;
    let absolute = if path.is_absolute() {
        path.to_path_buf()
    } else {
        std::env::current_dir()?.join(path)
    };
    let mut resolved = PathBuf::new();
    for part in absolute.components() {
        if !resolved.as_os_str().is_empty() {
            match fs::metadata(&resolved) {
                Ok(metadata) if !metadata.is_dir() => {
                    return Err(OrgraftError::InvalidArgument(
                        "workflow output traverses a non-directory".into(),
                    ))
                }
                Ok(_) => {}
                Err(cause) if cause.kind() == std::io::ErrorKind::NotFound => {}
                Err(cause) => return Err(cause.into()),
            }
        }
        match part {
            Component::Prefix(_) | Component::RootDir => resolved.push(part.as_os_str()),
            Component::CurDir => {}
            Component::ParentDir => {
                resolved.pop();
            }
            Component::Normal(name) => {
                resolved.push(name);
                match fs::symlink_metadata(&resolved) {
                    Ok(_) => resolved = fs::canonicalize(&resolved)?,
                    Err(cause) if cause.kind() == std::io::ErrorKind::NotFound => {}
                    Err(cause) => return Err(cause.into()),
                }
            }
        }
    }
    Ok(resolved)
}

fn validate_postrepair_sources(config: &WorkflowConfig) -> Result<(), OrgraftError> {
    let sources: Vec<_> = enabled_cases(config)
        .into_iter()
        .filter(|case| case.postrepair_source_dir.is_some())
        .collect();
    if sources.is_empty() {
        return Ok(());
    }
    let invalid = |message: &str| {
        OrgraftError::InvalidArgument(format!("postrepair source workflow: {message}"))
    };
    if config.recruit.enabled || config.asm.enabled {
        return Err(invalid(
            "postrepair_source_dir requires recruit and asm to be disabled",
        ));
    }
    if !config
        .postrepair
        .as_ref()
        .is_some_and(|repair| repair.enabled)
    {
        return Err(invalid(
            "postrepair_source_dir requires commands.postrepair.enabled=true",
        ));
    }
    for argument in config
        .polish
        .extra_args
        .iter()
        .chain(&config.rebuild.extra_args)
    {
        let flag = argument.split('=').next().unwrap_or(argument);
        if matches!(
            flag,
            "--out-dir"
                | "--output-dir"
                | "--linear-subgraph"
                | "--subgraph-fasta"
                | "--rotated-reference"
                | "--subgraph-reads"
                | "--subgraph-id"
                | "--workflow-round"
                | "--checked-draft-gfa"
                | "--draft-gfa"
                | "--raw-gfa"
                | "--gfa"
                | "--verified-fasta"
                | "--final-fasta"
                | "--fasta"
                | "--draft"
                | "--reference"
                | "--reads"
                | "--organelle"
                | "--subgraph"
                | "--source-dir"
                | "--checked-gfa"
                | "--resolve-details"
                | "--ordinary-source-dir"
                | "--first-polish-report"
                | "--edited-gfa"
                | "--polished-fasta"
                | "--joint-ir-trial-gfa"
                | "--joint-ir-evidence"
                | "--validate-round"
                | "--validate-only"
        ) {
            return Err(invalid(
                "extra_args cannot override native source or output bindings",
            ));
        }
    }
    let results = workflow_future_path(&config.results_dir)?;
    for case in &sources {
        if case.postrepair_out_dir.is_none() {
            return Err(invalid("source case requires a separate postrepair output"));
        }
        let round = postrepair_source_round(case)?.unwrap();
        let source = fs::canonicalize(case.postrepair_source_dir.as_ref().unwrap())?;
        let component = source
            .parent()
            .ok_or_else(|| invalid("source component is absent"))?;
        let organelle = component
            .parent()
            .ok_or_else(|| invalid("source organelle is absent"))?;
        let polish_root = organelle
            .parent()
            .ok_or_else(|| invalid("source polish root is absent"))?;
        if source.file_name() != case.postrepair_source_dir.as_ref().unwrap().file_name()
            || component.file_name().and_then(|name| name.to_str()) != Some(case.subgraph.as_str())
            || organelle.file_name().and_then(|name| name.to_str()) != Some(case.organelle.as_str())
        {
            return Err(invalid(
                "source directory identity differs from its workflow case",
            ));
        }
        if results.starts_with(polish_root) || polish_root.starts_with(&results) {
            return Err(invalid(
                "results_dir must be separate from the source polish root",
            ));
        }
        let report = homorepair_rows(&source.join("logs/report.tsv"))?;
        let value = |section: &str, name: &str, metric: &str| -> Result<String, OrgraftError> {
            let values: Vec<_> = report
                .iter()
                .filter(|row| {
                    row.get("section").map(String::as_str) == Some(section)
                        && row.get("name").map(String::as_str) == Some(name)
                        && row.get("metric").map(String::as_str) == Some(metric)
                })
                .filter_map(|row| row.get("value"))
                .collect();
            if values.len() != 1 {
                return Err(invalid(&format!(
                    "source report requires one {section}/{name}/{metric}"
                )));
            }
            Ok(values[0].clone())
        };
        for (metric, expected) in [
            ("organelle", case.organelle.clone()),
            ("subgraph", case.subgraph.clone()),
            ("validate_round", round.to_string()),
            ("per_read_variant_calls", "true".into()),
        ] {
            if value("global", "run", metric)? != expected {
                return Err(invalid(&format!(
                    "source report {metric} disagrees with the required input"
                )));
            }
        }
        for stage in ["sv-eval", "snv-indel-eval", "plot"] {
            if value("stage", &format!("{stage}:round_{round}"), "status")? != "ok" {
                return Err(invalid("source native validation is incomplete"));
            }
        }
        let fasta = repair_source_fasta(case, round);
        if fs::canonicalize(value("output", "validation_fasta", "path")?)?
            != fs::canonicalize(&fasta)?
            || read_single_fasta(&fasta)?.1.is_empty()
        {
            return Err(invalid(
                "source report does not bind its native validation FASTA",
            ));
        }
        let data = source.join("03.validate/01.data");
        let reports = source.join("03.validate/03.reports");
        if read_metric_value(&reports.join("sv_snv_indel_summary.tsv"), "status")?.as_deref()
            != Some("pass")
            || ordinary_validation_edit_count(&reports.join("snv_indel_high.tsv"))? != 0
        {
            return Err(invalid(
                "source still requires ordinary SV or small-variant correction",
            ));
        }
        // Native single-round SV sections use round_1 even when validate_round is N.
        if value("sv_eval", "round_1", "sv_support_status")? != "pass" {
            return Err(invalid("source native SV report is not pass"));
        }
        for file in [
            source.join("logs/report.tsv"),
            fasta,
            reports.join("sv_snv_indel_summary.tsv"),
            reports.join("snv_indel_high.tsv"),
            data.join("snv_indel_calls.tsv"),
            data.join("snv_indel_variants_combined.tsv"),
            source.join("01.inputs/subgraph_reads.fastq.gz"),
            source.join("01.inputs/rotated_reference.fasta"),
            case.reads.clone(),
            case.checked_draft_gfa.clone(),
            rebuild_source_gfa(config, case),
            resolve_details_path(case),
            first_polish_report_path(case).unwrap(),
        ] {
            if !file.is_file() {
                return Err(invalid(&format!(
                    "required native input is missing: {}",
                    file.display()
                )));
            }
            if fs::canonicalize(&file)?.starts_with(&results) {
                return Err(invalid(
                    "read-only source evidence must be outside the new results_dir",
                ));
            }
        }
        if let Some(reference) = &case.polish_reference {
            if fs::read(reference)? != fs::read(source.join("01.inputs/rotated_reference.fasta"))? {
                return Err(invalid(
                    "polish_reference differs from the source round reference",
                ));
            }
        }
    }
    // Every possible write for these configurations stays under the new root.
    // This also prevents an ordinary sibling case from writing a source round.
    for case in enabled_cases(config) {
        let mut outputs = vec![case.workflow_dir.clone(), case.rebuild_out_dir.clone()];
        if let Some(root) = &case.postrepair_out_dir {
            outputs.push(root.clone());
        }
        if case.postrepair_source_dir.is_none() {
            outputs.extend([
                case.resolve_out_dir.clone(),
                case.polish_out_dir.clone(),
                case.checked_draft_gfa.clone(),
            ]);
        }
        if config.refine.as_ref().is_some_and(|refine| refine.enabled) {
            outputs.extend([
                refine_out_dir(config, &case.organelle),
                refine_bundle_path(config, &case.organelle),
            ]);
        }
        for output in outputs {
            if !workflow_future_path(&output)?.starts_with(&results) {
                return Err(invalid(&format!(
                    "output must remain under results_dir: {}",
                    output.display()
                )));
            }
        }
    }
    Ok(())
}

fn reject_source_checkpoint(case: &WorkflowCase) -> Result<(), OrgraftError> {
    if case.postrepair_source_dir.is_some() {
        return Err(OrgraftError::InvalidArgument("this case starts at postrepair; its source round is read-only and no new checkpoint 1 or 2 may be written".into()));
    }
    Ok(())
}

fn write_postrepair_source_stages(
    script: &mut String,
    config: &WorkflowConfig,
    case: &WorkflowCase,
) -> Result<(), OrgraftError> {
    let round = postrepair_source_round(case)?
        .ok_or_else(|| OrgraftError::InvalidArgument("postrepair_source_dir is absent".into()))?;
    writeln!(
        script,
        "# Start from read-only native validation: {}",
        repair_source_dir(case, round).display()
    )
    .unwrap();
    writeln!(
        script,
        "{} workflow postrepair-check --config {} --case {} --source-only > /dev/null",
        orgraft_shell_token(),
        shell_quote(&config.config_path),
        shell_quote_str(&case.name)
    )
    .unwrap();
    write_homorepair_after_checkpoint(script, config, case, round)?;
    if config.rebuild.enabled {
        write_stage_header(
            script,
            "",
            "06.rebuild",
            "rebuild from the strictly verified postrepair handoff",
        );
        writeln!(
            script,
            "{}",
            rebuild_command_with_polished_arg(
                config,
                case,
                config.force,
                "\"${final_polished}\"",
                None
            )
        )
        .unwrap();
    }
    Ok(())
}

fn refine_bundle_path(config: &WorkflowConfig, organelle: &str) -> PathBuf {
    config.refine.as_ref().and_then(|r| r.bundle.as_ref())
        .map(|path| PathBuf::from(path.to_string_lossy().replace("${organelle}", organelle)))
        .unwrap_or_else(|| config.results_dir.join("workflow").join(organelle).join("refine.inputs.tsv"))
}

fn refine_out_dir(config: &WorkflowConfig, organelle: &str) -> PathBuf {
    config.refine.as_ref().map(|r| r.out_dir.join(organelle))
        .unwrap_or_else(|| config.results_dir.join("07.refine").join(organelle))
}

fn refine_reads(config: &WorkflowConfig, organelle: &str) -> Result<PathBuf, OrgraftError> {
    if let Some(path) = config.refine.as_ref().and_then(|r| r.organelle_reads.as_ref()) {
        return Ok(PathBuf::from(path.to_string_lossy().replace("${organelle}", organelle)));
    }
    let reads: BTreeSet<_> = enabled_cases(config).into_iter()
        .filter(|case| case.organelle == organelle).map(|case| &case.reads).collect();
    if reads.len() == 1 {
        return Ok((*reads.first().unwrap()).clone());
    }
    Err(OrgraftError::InvalidArgument(format!(
        "refine requires commands.refine.organelle_reads for {organelle}: enabled components do not share one reads file; set it to the full organelle read pool"
    )))
}

fn refine_downgrade_representation(
    config: &WorkflowConfig,
    organelle: &str,
) -> Result<commands::refine::DowngradeRepresentation, OrgraftError> {
    let values: BTreeSet<_> = enabled_cases(config)
        .into_iter()
        .filter(|case| case.organelle == organelle)
        .map(|case| case.refine_downgrade_representation)
        .collect();
    if values.len() != 1 {
        return Err(OrgraftError::InvalidArgument(format!(
            "refine requires one consistent refine_downgrade_representation across enabled {organelle} cases"
        )));
    }
    Ok(*values.first().unwrap())
}

fn refine_intra_measurement(config: &WorkflowConfig, organelle: &str) -> commands::refine::IntraMeasurement {
    config.refine.as_ref().and_then(|r| r.intra_measurement.get(organelle)).copied()
        .unwrap_or(commands::refine::IntraMeasurement::ExactCore)
}

fn refine_duplicate_template_full_span(config: &WorkflowConfig, organelle: &str) -> bool {
    config
        .refine
        .as_ref()
        .and_then(|r| r.duplicate_template_full_span.get(organelle))
        .copied()
        .unwrap_or(false)
}

fn refine_linear_chain_flanks(config: &WorkflowConfig, organelle: &str) -> bool {
    config
        .refine
        .as_ref()
        .and_then(|r| r.linear_chain_flanks.get(organelle))
        .copied()
        .unwrap_or(false)
}

fn refine_graph_measurement(config: &WorkflowConfig, organelle: &str) -> commands::refine::GraphMeasurement {
    config.refine.as_ref().and_then(|r|r.graph_measurement.get(organelle)).copied()
        .unwrap_or(commands::refine::GraphMeasurement::Default)
}

fn refine_graph_full_family(config: &WorkflowConfig, organelle: &str) -> bool {
    config
        .refine
        .as_ref()
        .and_then(|r| r.graph_full_family.get(organelle))
        .copied()
        .unwrap_or(false)
}

fn refine_args(config: &WorkflowConfig, organelle: &str) -> Result<Vec<String>, OrgraftError> {
    let mut args = vec![
        "--sample-bundle".into(),
        refine_bundle_path(config, organelle).display().to_string(),
        "--reads".into(),
        refine_reads(config, organelle)?.display().to_string(),
        "--sample".into(),
        config.sample.clone(),
        "--organelle".into(),
        organelle.into(),
        "--out-dir".into(),
        refine_out_dir(config, organelle).display().to_string(),
        "--downgrade-representation".into(),
        refine_downgrade_representation(config, organelle)?
            .as_str()
            .into(),
        "--auto".into(),
        "--threads".into(),
        config
            .refine
            .as_ref()
            .and_then(|r| r.threads)
            .unwrap_or(config.threads)
            .to_string(),
        "--soft-paths".into(),
        config.soft_paths.display().to_string(),
    ];
    if config.refine.as_ref().is_some_and(|r| r.intra_measurement.contains_key(organelle)) {
        args.extend(["--intra-measurement".into(), refine_intra_measurement(config, organelle).as_str().into()]);
    }
    if refine_duplicate_template_full_span(config, organelle) {
        args.push("--duplicate-template-full-span".into());
    }
    if refine_linear_chain_flanks(config, organelle) {
        args.push("--linear-chain-flanks".into());
    }
    let graph_measurement=refine_graph_measurement(config,organelle);
    if graph_measurement != commands::refine::GraphMeasurement::Default {
        args.extend(["--graph-measurement".into(),graph_measurement.as_str().into()]);
    }
    if refine_graph_full_family(config, organelle) {
        args.push("--graph-full-family".into());
    }
    if let Some(path) = config.refine.as_ref().and_then(|r| r.adopt.get(organelle)) {
        if !path.is_file() {
            return Err(OrgraftError::InvalidArgument(format!(
                "commands.refine.adopt.{organelle} is not a readable file: {}",
                path.display()
            )));
        }
        fs::File::open(path)?;
        args.extend(["--adopt".into(), path.display().to_string()]);
    }
    Ok(args)
}

fn write_refine_stage(script: &mut String, config: &WorkflowConfig, organelle: &str) -> Result<(), OrgraftError> {
    if !config.refine.as_ref().is_some_and(|r| r.enabled) { return Ok(()); }
    write_stage_header(script, "", "07.refine", "review all rebuilt components together for this organelle");
    writeln!(script, "{} workflow refine-bundle --config {} --organelle {}", orgraft_shell_token(), shell_quote(&config.config_path), shell_quote_str(organelle)).unwrap();
    writeln!(script, "{} refine {}", orgraft_shell_token(), refine_args(config, organelle)?.iter().map(|arg| shell_quote_str(arg)).collect::<Vec<_>>().join(" ")).unwrap();
    writeln!(script, "{} workflow refine-check --config {} --organelle {}", orgraft_shell_token(), shell_quote(&config.config_path), shell_quote_str(organelle)).unwrap();
    Ok(())
}

fn refine_organelle(config: &WorkflowConfig, args: &[String]) -> Result<String, OrgraftError> {
    if let Some(organelle) = option_value(args, "--organelle")? {
        if enabled_cases(config).iter().any(|case| case.organelle == organelle) { return Ok(organelle.into()); }
        return Err(OrgraftError::InvalidArgument(format!("no enabled cases for organelle {organelle}")));
    }
    Ok(select_case(config, args)?.organelle.clone())
}

fn run_refine_bundle(args: &[String]) -> Result<(), OrgraftError> {
    let config = load_config_from_args(args)?;
    let organelle = refine_organelle(&config, args)?;
    let path = prepare_refine_bundle(&config, &organelle)?;
    println!("{}", path.display());
    Ok(())
}

fn refine_validation_reference(validation: &Path) -> Result<PathBuf, OrgraftError> {
    let round = validation.parent().ok_or_else(|| {
        OrgraftError::InvalidArgument("refine validation has no owning round".into())
    })?;
    Ok(round.join("01.inputs/rotated_reference.fasta"))
}

fn prepare_refine_bundle(config: &WorkflowConfig, organelle: &str) -> Result<PathBuf, OrgraftError> {
    let path = refine_bundle_path(config, organelle);
    let cases: Vec<_> = enabled_cases(config).into_iter().filter(|c| c.organelle == organelle).collect();
    if cases.is_empty() { return Err(OrgraftError::InvalidArgument(format!("refine pending: no enabled {organelle} components"))); }
    refine_reads(config, organelle)?;
    refine_downgrade_representation(config, organelle)?;
    // First check the whole organelle. Do not write a partial bundle if another
    // component still needs a checkpoint, postrepair, or rebuild decision.
    let mut inputs = Vec::new();
    for case in &cases {
        let round = match postrepair_source_round(case)? { Some(round) => Some(round), None => completed_round(case)? }.ok_or_else(|| OrgraftError::InvalidArgument(format!("refine pending: {} has no complete checkpoint 2", case.name)))?;
        let gfa = case.rebuild_out_dir.join(&case.subgraph).join(format!("rebuild_{}.gfa", case.subgraph));
        let fasta = rebuild_fasta_path(case);
        if !gfa.is_file() || !fasta.is_file() {
            return Err(OrgraftError::InvalidArgument(format!("refine pending: {} rebuild is not complete", case.name)));
        }
        let mut validation = repair_source_dir(case, round).join("03.validate");
        if postrepair_enabled(config) {
            selected_repair_fasta(config, case, round)?;
            if case.postrepair_out_dir.is_some() {
                let component = postrepair_component_for_source(case, round)?;
                validation = commands::postrepair::selected_validation_dir(&component)?;
            } else {
                let outcome = read_homorepair_run_outcome(&homorepair_outcome_path(case))?;
                if let Some(final_dir) = outcome.get("final_validation_dir") { validation = PathBuf::from(final_dir).join("03.validate"); }
            }
        }
        inputs.push((*case, gfa, fasta, validation));
    }
    if config.refine.as_ref().is_some_and(|r| r.bundle.is_some()) {
        if !path.is_file() { return Err(OrgraftError::InvalidArgument(format!("refine sample_bundle does not exist: {}", path.display()))); }
        return Ok(path);
    }
    let mut text = String::from("component_id\tgfa\tcomponent_fasta\tvalidate_data_dir\tcomponent_reads\trotated_reference\tsource_contract\n");
    for (case, gfa, fasta, validation) in inputs {
        let report = rebuild_report_path(case);
        let stats = rebuild_result_stats_path(case);
        let source = rebuild_source_gfa(config, case);
        // The final validation owns the component reference, including after a repair.
        let reference = refine_validation_reference(&validation)?;
        let contract = case.workflow_dir.join("refine_source.tsv");
        if let Some(parent) = contract.parent() { fs::create_dir_all(parent)?; }
        commands::refine::prepare_rebuild_component(&case.subgraph, &gfa, &fasta, &report, &stats,
            &source, &resolve_details_path(case), &validation, &case.reads, &reference, &contract)?;
        let fields = [case.subgraph.clone(), gfa.display().to_string(), fasta.display().to_string(),
            validation.display().to_string(), case.reads.display().to_string(),
            reference.display().to_string(), contract.display().to_string()];
        if fields.iter().any(|field| field.contains(['\t','\n','\r'])) {
            return Err(OrgraftError::InvalidArgument("refine bundle fields must not contain tabs or newlines".into()));
        }
        writeln!(text, "{}", fields.join("\t")).unwrap();
    }
    if path.exists() {
        if fs::read_to_string(&path)? != text { return Err(OrgraftError::InvalidArgument(format!("refine bundle inputs changed; review existing {} before replacement", path.display()))); }
    } else { write_with_parent(&path, &text)?; }
    Ok(path)
}

fn check_refine_status(config: &WorkflowConfig, organelle: &str) -> Result<(), OrgraftError> {
    let record = refine_out_dir(config, organelle).join("auto_complete.tsv");
    let status = read_metric_value(&record, "status")?.unwrap_or_else(|| "missing".into());
    if !matches!(status.as_str(), "confirmed" | "no_change") {
        return Err(OrgraftError::InvalidArgument(format!("refine {status}: inspect {}; no final acceptance was inferred", record.display())));
    }
    commands::refine::verify_intra_measurement(&refine_out_dir(config, organelle),
        refine_intra_measurement(config, organelle))?;
    commands::refine::verify_duplicate_template_full_span(
        &refine_out_dir(config, organelle),
        refine_duplicate_template_full_span(config, organelle),
    )?;
    commands::refine::verify_linear_chain_flanks(&refine_out_dir(config, organelle),
        refine_linear_chain_flanks(config, organelle))?;
    commands::refine::verify_graph_measurement(&refine_out_dir(config,organelle),refine_graph_measurement(config,organelle))?;
    commands::refine::verify_graph_full_family(
        &refine_out_dir(config, organelle),
        refine_graph_full_family(config, organelle),
    )?;
    println!("refine {organelle}: {status}");
    Ok(())
}

fn run_refine_check(args: &[String]) -> Result<(), OrgraftError> {
    let config = load_config_from_args(args)?;
    let organelle = refine_organelle(&config, args)?;
    check_refine_status(&config, &organelle)
}

fn run_refine_for_organelle(config: &WorkflowConfig, organelle: &str) -> Result<(), OrgraftError> {
    if !config.refine.as_ref().is_some_and(|r| r.enabled) { return Ok(()); }
    prepare_refine_bundle(config, organelle)?;
    commands::refine::run(&refine_args(config, organelle)?)?;
    check_refine_status(config, organelle)
}

fn run_correct(args: &[String]) -> Result<(), OrgraftError> {
    let input = option_value(args, "--input-fasta")?
        .or_else(|| option_value(args, "--old-fasta").ok().flatten())
        .ok_or_else(|| OrgraftError::InvalidArgument("missing --input-fasta FILE".to_string()))?;
    let variants = option_value(args, "--pos-ref-alt")?
        .ok_or_else(|| OrgraftError::InvalidArgument("missing --pos-ref-alt FILE".to_string()))?;
    let out = option_value(args, "--out")?.ok_or_else(|| {
        OrgraftError::InvalidArgument("missing --out FILE for corrected FASTA".to_string())
    })?;

    let edits = read_pos_ref_alt(Path::new(variants))?;
    let summary = apply_edits_to_fasta(Path::new(input), Path::new(out), &edits)?;
    println!("Wrote {out}");
    println!("{summary}");
    Ok(())
}

fn run_test_correction(args: &[String]) -> Result<(), OrgraftError> {
    let input = option_value(args, "--input-fasta")?
        .ok_or_else(|| OrgraftError::InvalidArgument("missing --input-fasta FILE".to_string()))?;
    let variants = option_value(args, "--pos-ref-alt")?
        .ok_or_else(|| OrgraftError::InvalidArgument("missing --pos-ref-alt FILE".to_string()))?;
    let out_dir = option_value(args, "--out-dir")?
        .map(PathBuf::from)
        .unwrap_or_else(|| PathBuf::from("results_workflow/correction_test"));
    let force = has_flag(args, "--force");
    if out_dir.exists() && !force {
        return Err(OrgraftError::InvalidArgument(format!(
            "{} already exists; pass --force to replace test outputs",
            out_dir.display()
        )));
    }
    if out_dir.exists() {
        fs::remove_dir_all(&out_dir)?;
    }
    fs::create_dir_all(&out_dir)?;

    let corrected = out_dir.join("corrected.fasta");
    let edits = read_pos_ref_alt(Path::new(variants))?;
    let summary = apply_edits_to_fasta(Path::new(input), &corrected, &edits)?;
    fs::write(
        out_dir.join("summary.tsv"),
        format!(
            "metric\tvalue\ninput_fasta\t{}\npos_ref_alt\t{}\ncorrected_fasta\t{}\n{}\n",
            input,
            variants,
            corrected.display(),
            correction_summary_tsv(&summary)
        ),
    )?;
    println!("Wrote {}", corrected.display());
    Ok(())
}

fn run_test_fake_validate(args: &[String]) -> Result<(), OrgraftError> {
    let input = option_value(args, "--input-fasta")?
        .or_else(|| option_value(args, "--error-fasta").ok().flatten())
        .ok_or_else(|| {
            OrgraftError::InvalidArgument(
                "missing --input-fasta FILE for fake validate smoke test".to_string(),
            )
        })?;
    let variants = option_value(args, "--pos-ref-alt")?.ok_or_else(|| {
        OrgraftError::InvalidArgument(
            "missing --pos-ref-alt FILE for fake validate smoke test".to_string(),
        )
    })?;
    let out_dir = option_value(args, "--out-dir")?
        .map(PathBuf::from)
        .unwrap_or_else(|| PathBuf::from("results_workflow/fake_validate"));
    let organelle = option_value(args, "--organelle")?.unwrap_or("mito");
    let subgraph = option_value(args, "--subgraph")?.unwrap_or("subgraph_001");
    let force = has_flag(args, "--force");

    if out_dir.exists() && !force {
        return Err(OrgraftError::InvalidArgument(format!(
            "{} already exists; pass --force to replace fake validate outputs",
            out_dir.display()
        )));
    }
    if out_dir.exists() {
        fs::remove_dir_all(&out_dir)?;
    }
    fs::create_dir_all(&out_dir)?;

    let config_path = out_dir.join("orgraft.fake-validate.toml");
    let config = WorkflowConfig {
        config_path: config_path.clone(),
        sample: "fake_validate".to_string(),
        results_dir: out_dir.clone(),
        soft_paths: PathBuf::from("soft_paths.txt"),
        mode: "stepwise".to_string(),
        max_rounds: 3,
        threads: 1,
        force: true,
        auto_sv_correction: false,
        auto_snv_indel_correction: true,
        recruit: WorkflowRecruitConfig {
            enabled: false,
            reads: None,
            out_dir: out_dir.join("recruit"),
            threads: None,
            baits: Vec::new(),
            prefix: None,
            bait_format: None,
            gfa_split: None,
            rename_bait: false,
            write_id_map: false,
            split_output: None,
            gzip_output: None,
            minimap2: None,
            align_mode: None,
            platform: None,
            preset: None,
            min_mapq: None,
            min_aln_len: None,
            sam: None,
            max_reads: Vec::new(),
            random_seed: None,
            write_sampled_ids: false,
            read_stats: None,
            write_read_classification: false,
            write_bait_partitions: false,
            gzip_tool: None,
            mode: None,
            iterations: None,
            extra_args: Vec::new(),
        },
        asm: WorkflowAsmConfig {
            enabled: false,
            out_dir: out_dir.join("draft_asm"),
            threads: None,
            profile: None,
            stable: false,
            min_graph_coverage: None,
            min_branch_ratio: None,
            min_tip_len: None,
            min_link_support: None,
            min_link_ratio: None,
            subsets: None,
            image_reference_fasta: None,
            keep_debug_files: false,
            extra_args: Vec::new(),
        },
        resolve: WorkflowResolveConfig {
            out_dir: out_dir.join("resolve_gfa"),
            gfa_editor_mode: None,
            max_states: None,
            max_candidates: None,
            extra_args: Vec::new(),
        },
        polish: WorkflowPolishConfig {
            out_dir: out_dir.join("polish"),
            homorepair: false,
            homorepair_frq: 0.2,
            threads: None,
            per_read_variant_calls: None,
            snv_indel_overlap_policy: None,
            plot_range: None,
            plot_dpi: None,
            plot_output_format: None,
            coverage_plot_rasterize: None,
            snv_indel_plot_rasterize: None,
            sv_plot_highlight_subgroups: None,
            sv_plot_highlight_read_ids: None,
            sv_plot_highlight_min_fraction: None,
            sv_plot_highlight_min_reads: None,
            snv_indel_plot_low_confidence: None,
            snv_indel_plot_low_min_reads: None,
            snv_indel_plot_low_min_fraction: None,
            snv_indel_plot_high_risk_fraction: None,
            extra_args: Vec::new(),
        },
        postrepair: None,
        refine: None,
        rebuild: WorkflowRebuildConfig {
            enabled: false,
            out_dir: out_dir.join("rebuild"),
            threads: None,
            edited_gfa: None,
            polished_fasta: None,
            image_reference_fasta: None,
            merged_gfa_template: None,
            minimap2: None,
            blastn: None,
            keep_debug: false,
            extra_args: Vec::new(),
        },
        topology_simple_allowed_classes: default_simple_classes(),
        cases: vec![WorkflowCase {
            enabled: true,
            name: format!("{organelle}_{subgraph}"),
            sample: "fake_validate".to_string(),
            organelle: organelle.to_string(),
            subgraph: subgraph.to_string(),
            draft_graph: out_dir.join("draft.gfa"),
            unitig_graph: None,
            checked_draft_gfa: out_dir
                .join("workflow")
                .join(organelle)
                .join(subgraph)
                .join("checkpoint_1/checked_draft.gfa"),
            reference: None,
            pre_rotated_reference: None,
            reads: out_dir.join("reads.fastq.gz"),
            asm_reads: None,
            resolve_out_dir: out_dir.join("resolve_gfa"),
            polish_out_dir: out_dir.join("polish"),
            postrepair_out_dir: None,
            postrepair_source_dir: None,
            resolve_details: None,
            first_polish_report: None,
            rebuild_out_dir: out_dir.join("rebuild"),
            rebuild_edited_gfa: None,
            rebuild_polished_fasta: None,
            image_reference_fasta: None,
            workflow_dir: out_dir.join("workflow").join(organelle).join(subgraph),
            linearized_fasta: None,
            polish_reference: None,
            sv_correction_subgroup: None,
            refine_downgrade_representation: commands::refine::DowngradeRepresentation::Recorded,
        }],
    };
    write_with_parent(
        &config_path,
        "# Generated fake validate config for workflow checkpoint2 smoke testing.\n",
    )?;
    let case = &config.cases[0];

    let edits = read_pos_ref_alt(Path::new(variants))?;
    if edits.is_empty() {
        return Err(OrgraftError::InvalidArgument(
            "--pos-ref-alt must contain at least one fake validation error".to_string(),
        ));
    }

    let round1_polished = polished_aln_path(case, 1);
    copy_with_parent(Path::new(input), &round1_polished, true)?;
    write_fake_sv_pass_summary(&sv_summary_path(case, 1))?;
    write_fake_snv_high_table(&snv_indel_high_path(case, 1), &edits)?;

    let round1_status = checkpoint2_impl(&config, case, 1, true, false, None)?;
    if round1_status != Checkpoint2Status::NextRoundReady {
        return Err(OrgraftError::InvalidArgument(format!(
            "fake round_1 expected next_round_ready, got {round1_status:?}"
        )));
    }

    let round2_polished = polished_aln_path(case, 2);
    copy_with_parent(&corrected_fasta_path(case, 1), &round2_polished, true)?;
    write_fake_sv_pass_summary(&sv_summary_path(case, 2))?;
    write_fake_snv_high_table(&snv_indel_high_path(case, 2), &[])?;

    let round2_status = checkpoint2_impl(&config, case, 2, true, false, None)?;
    if round2_status != Checkpoint2Status::Complete {
        return Err(OrgraftError::InvalidArgument(format!(
            "fake round_2 expected complete, got {round2_status:?}"
        )));
    }

    write_fake_validate_summary(&out_dir.join("summary.tsv"), input, variants, case)?;
    println!("Wrote {}", out_dir.join("summary.tsv").display());
    Ok(())
}

fn run_automatic_workflow(args: &[String]) -> Result<(), OrgraftError> {
    let config = load_config_from_args(args)?;
    let case = select_case(&config, args)?;
    let force = config.force || has_flag(args, "--force");

    if let Some(round) = postrepair_source_round(case)? {
        let selected = run_homorepair_after_checkpoint(&config, case, round)?;
        if config.rebuild.enabled {
            commands::rebuild::run(&rebuild_args_with_polished(&config, case, force, Some(&selected), None)?)?;
        }
        run_refine_for_organelle(&config, &case.organelle)?;
        return Ok(());
    }

    if config.recruit.enabled {
        commands::recruit::run(&recruit_args(&config)?)?;
    }
    if config.asm.enabled {
        commands::asm::run(&asm_args(&config, case))?;
    }

    if checkpoint1_impl(&config, case, force)? != Checkpoint1Status::Checked {
        return Ok(());
    }

    run_resolve_for_case(&config, case, force)?;
    run_polish_round(&config, case, 1, &round_draft_path(case, 1), force)?;

    for round in 1..=WORKFLOW_HARD_MAX_ROUNDS {
        match checkpoint2_impl(&config, case, round, force, false, None)? {
            Checkpoint2Status::Complete => {
                let final_fasta = if postrepair_enabled(&config) {
                    run_homorepair_after_checkpoint(&config, case, round)?
                } else {
                    polished_aln_path(case, round)
                };
                if config.rebuild.enabled {
                    let marker = read_repeat_pairing_marker(case)?;
                    let repeat_node = marker
                        .as_ref()
                        .filter(|marker| marker.state == "validated")
                        .map(|marker| marker.repeat_node.as_str());
                    let rebuild_fasta = rebuild_fasta_path(case);
                    let already_current = repeat_node.is_some()
                        && rebuild_fasta.is_file()
                        && files_equal(&final_fasta, &rebuild_fasta)?;
                    if !already_current {
                        commands::rebuild::run(&rebuild_args_with_polished(
                            &config,
                            case,
                            force || rebuild_fasta.exists(),
                            Some(&final_fasta),
                            repeat_node,
                        )?)?;
                    }
                }
                run_refine_for_organelle(&config, &case.organelle)?;
                return Ok(());
            }
            Checkpoint2Status::ManualRequired => return Ok(()),
            Checkpoint2Status::NextRoundReady => {
                let next_round = round + 1;
                if next_round > WORKFLOW_HARD_MAX_ROUNDS {
                    return Ok(());
                }
                let draft = corrected_fasta_path(case, round);
                run_polish_round(&config, case, next_round, &draft, force)?;
            }
            Checkpoint2Status::RebuildReady => {
                let marker = read_repeat_pairing_marker(case)?.ok_or_else(|| {
                    OrgraftError::InvalidArgument(
                        "checkpoint2 returned rebuild_ready without a repeat pairing marker"
                            .to_string(),
                    )
                })?;
                run_repeat_rebuild_then_validate(
                    &config,
                    case,
                    round,
                    &polished_aln_path(case, round),
                    &marker,
                    force,
                )?;
            }
        }
    }

    Ok(())
}

fn checkpoint1_impl(
    config: &WorkflowConfig,
    case: &WorkflowCase,
    force: bool,
) -> Result<Checkpoint1Status, OrgraftError> {
    reject_source_checkpoint(case)?;
    let checkpoint_dir = case.workflow_dir.join("checkpoint_1");
    fs::create_dir_all(&checkpoint_dir)?;

    let reader = BufReader::new(File::open(&case.draft_graph)?);
    let report = analyze_gfa(reader)?;
    let missing_segments = missing_link_segments(&case.draft_graph)?;
    let complex_nodes = complex_nodes(&report, &config.topology_simple_allowed_classes);
    let is_checked = missing_segments.is_empty() && complex_nodes.is_empty();

    fs::write(
        checkpoint_dir.join("topology_nodes.tsv"),
        nodes_tsv(&report),
    )?;
    fs::write(
        checkpoint_dir.join("topology_summary.tsv"),
        summary_tsv(&report),
    )?;

    if is_checked {
        copy_with_parent(&case.draft_graph, &case.checked_draft_gfa, force)?;
        write_checkpoint1_status(
            &checkpoint_dir.join("checkpoint_1.status.tsv"),
            "checked",
            "topology is simple and all link endpoints reference declared segments",
            case,
            &report,
            &missing_segments,
            &complex_nodes,
        )?;
        Ok(Checkpoint1Status::Checked)
    } else {
        let manual_gfa = checkpoint_dir.join("manual_edit_required.gfa");
        copy_with_parent(&case.draft_graph, &manual_gfa, true)?;
        write_checkpoint1_status(
            &checkpoint_dir.join("checkpoint_1.status.tsv"),
            "manual_required",
            "complex topology or inconsistent GFA links require manual graph editing",
            case,
            &report,
            &missing_segments,
            &complex_nodes,
        )?;
        Ok(Checkpoint1Status::ManualRequired)
    }
}

fn checkpoint2_impl(
    config: &WorkflowConfig,
    case: &WorkflowCase,
    round: usize,
    force: bool,
    run_next: bool,
    manual_sv_subgroup: Option<&str>,
) -> Result<Checkpoint2Status, OrgraftError> {
    reject_source_checkpoint(case)?;
    if round == 0 {
        return Err(OrgraftError::InvalidArgument(
            "--round must be greater than 0".to_string(),
        ));
    }
    let checkpoint_dir = case
        .workflow_dir
        .join("checkpoint_2")
        .join(format!("round_{round}"));
    if force && checkpoint_dir.exists() {
        fs::remove_dir_all(&checkpoint_dir)?;
    }
    fs::create_dir_all(&checkpoint_dir)?;

    let summary_path = sv_summary_path(case, round);
    let high_path = snv_indel_high_path(case, round);
    let polished_aln = polished_aln_path(case, round);
    let sv_status = read_metric_value(&summary_path, "status")?;
    let requested_sv_subgroup = manual_sv_subgroup.or(case.sv_correction_subgroup.as_deref());
    let configured_sv_subgroup = match requested_sv_subgroup {
        Some(spec) if sv_subgroup_already_corrected(case, round, spec)? => None,
        other => other,
    };

    // Record a dominant alternate repeat pairing, but defer topology selection until
    // rebuild has measured complete flank-repeat-flank support for all four paths.
    if round == 1 && force {
        remove_path_if_exists(&repeat_pairing_marker_path(case))?;
    }
    let mut repeat_marker = read_repeat_pairing_marker(case)?;
    let repeat_rebuild_applied = match repeat_marker.as_ref() {
        Some(marker) if marker.state == "rebuild_ready" => {
            rebuild_has_repeat_resolution(case, &marker.repeat_node)?
        }
        _ => false,
    };

    if config.auto_sv_correction
        && configured_sv_subgroup.is_none()
        && (repeat_marker.is_none() || repeat_rebuild_applied)
    {
        let unitig_graph = discover_unitig_graph(config, case);
        if let Some(mismatch) = detect_repeat_pairing_mismatch(
            &sv_subgroup_stats_path(case, round),
            &sv_read_index_path(case, round),
            &sv_read_evidence_path(case, round),
            unitig_graph.as_deref(),
        )? {
            let repair_dir = checkpoint_dir.join("repeat_pairing_repair");
            fs::create_dir_all(&repair_dir)?;
            let mismatch_report = repair_dir.join("repeat_pairing_mismatch.tsv");
            write_repeat_pairing_mismatch_report(&mismatch_report, &mismatch)?;

            if repeat_rebuild_applied
                && repeat_marker
                    .as_ref()
                    .is_some_and(|marker| marker.repeat_node == mismatch.repeat_node)
            {
                let status_path = checkpoint_dir.join("checkpoint_2.status.tsv");
                write_checkpoint2_status(
                    &status_path,
                    "manual_required",
                    "repeat pairing mismatch remains after rebuild-constrained FASTA validation",
                    round,
                    &summary_path,
                    &high_path,
                    None,
                    None,
                )?;
                append_repeat_pairing_metrics(
                    &status_path,
                    &mismatch,
                    &mismatch_report,
                    "repeat_pairing_rebuild_validation_failed",
                )?;
                if let Some(marker) = repeat_marker.as_mut() {
                    marker.state = "failed".to_string();
                    marker.validation_round = Some(round);
                    write_repeat_pairing_marker(case, marker)?;
                }
                write_next_round_script(
                    &checkpoint_dir.join("next_round.sh"),
                    config,
                    case,
                    round + 1,
                    None,
                    Some("repeat pairing remains inconsistent after rebuild constraint"),
                )?;
                return Ok(Checkpoint2Status::ManualRequired);
            }

            let marker = RepeatPairingMarker {
                state: "pending".to_string(),
                repeat_node: mismatch.repeat_node.clone(),
                source_round: round,
                mismatch_report: mismatch_report.clone(),
                rebuild_fasta: None,
                validation_round: None,
            };
            write_repeat_pairing_marker(case, &marker)?;
            repeat_marker = Some(marker);
        }
    }

    if configured_sv_subgroup.is_some()
        || (config.auto_sv_correction && sv_status.as_deref() != Some("pass"))
    {
        let sv_high_path = sv_high_subgroups_path(case, round);
        let selected_sv_spec = if sv_high_path.exists() {
            select_sv_subgroup_spec(&sv_high_path, configured_sv_subgroup)?
        } else {
            None
        };
        if let Some(selected_spec) = selected_sv_spec {
            if round >= WORKFLOW_HARD_MAX_ROUNDS {
                write_checkpoint2_status(
                    &checkpoint_dir.join("checkpoint_2.status.tsv"),
                    "manual_required",
                    "SV correction is still required, but the hard total of 10 validation rounds has been reached",
                    round,
                    &summary_path,
                    &high_path,
                    None,
                    None,
                )?;
                append_checkpoint2_metrics(
                    &checkpoint_dir.join("checkpoint_2.status.tsv"),
                    &[
                        ("correction_kind", "sv".to_string()),
                        ("sv_subgroup", selected_spec),
                        ("hard_max_rounds", WORKFLOW_HARD_MAX_ROUNDS.to_string()),
                    ],
                )?;
                return Ok(Checkpoint2Status::ManualRequired);
            }

            let corrected = corrected_fasta_path(case, round);
            if corrected.exists() && !force {
                return Err(OrgraftError::InvalidArgument(format!(
                    "{} already exists; pass --force to overwrite",
                    corrected.display()
                )));
            }
            let prior_fastas = checkpoint_history_fastas(case, round);
            let repair_dir = checkpoint_dir.join("sv_repair");
            let unitig_graph = discover_unitig_graph(config, case);
            let repair = repair_sv_subgroup(&SvRepairRequest {
                reference: &polished_aln,
                reads: &case.reads,
                soft_paths: &config.soft_paths,
                summary: &summary_path,
                high_subgroups: &sv_high_path,
                segments: &snv_indel_segments_path(case, round),
                read_index: &sv_read_index_path(case, round),
                unitig_graph: unitig_graph.as_deref(),
                output_dir: &repair_dir,
                manual_subgroup: Some(&selected_spec),
                threads: polish_threads(config).min(8),
                prior_fastas: &prior_fastas,
                max_candidate_evaluations: None,
            })?;

            let Some(repair) = repair else {
                let repair_report = repair_dir.join("sv_repair.tsv");
                write_checkpoint2_status(
                    &checkpoint_dir.join("checkpoint_2.status.tsv"),
                    "manual_required",
                    "SV subgroup is abnormal, but no candidate both fixed it and improved global reference support",
                    round,
                    &summary_path,
                    &high_path,
                    None,
                    None,
                )?;
                append_checkpoint2_metrics(
                    &checkpoint_dir.join("checkpoint_2.status.tsv"),
                    &[
                        ("correction_kind", "sv".to_string()),
                        ("sv_subgroup", selected_spec),
                        ("sv_repair_report", repair_report.display().to_string()),
                        (
                            "sv_graph_localization",
                            repair_dir
                                .join("sv_graph_localization.tsv")
                                .display()
                                .to_string(),
                        ),
                    ],
                )?;
                write_next_round_script(
                    &checkpoint_dir.join("next_round.sh"),
                    config,
                    case,
                    round + 1,
                    None,
                    Some("SV candidate search requires manual review"),
                )?;
                return Ok(Checkpoint2Status::ManualRequired);
            };

            copy_with_parent(&repair.corrected_fasta, &corrected, force)?;
            let sv_corrections = sv_correction_count(case, round)? + 1;
            let effective_max_round =
                (config.max_rounds + sv_corrections).min(WORKFLOW_HARD_MAX_ROUNDS);
            write_checkpoint2_status(
                &checkpoint_dir.join("checkpoint_2.status.tsv"),
                "next_round_ready",
                "one abnormal SV subgroup was corrected and validated; continue with all reads without consuming the ordinary max_rounds budget",
                round,
                &summary_path,
                &high_path,
                None,
                Some(&corrected),
            )?;
            append_checkpoint2_metrics(
                &checkpoint_dir.join("checkpoint_2.status.tsv"),
                &[
                    ("correction_kind", "sv".to_string()),
                    ("sv_subgroup", repair.subgroup_spec.clone()),
                    (
                        "sv_repair_report",
                        repair.repair_report.display().to_string(),
                    ),
                    (
                        "sv_candidate_scores",
                        repair.score_table.display().to_string(),
                    ),
                    (
                        "sv_graph_localization",
                        repair.graph_localization.report.display().to_string(),
                    ),
                    (
                        "sv_graph_problem_scope",
                        repair.graph_localization.problem_scope.clone(),
                    ),
                    (
                        "sv_graph_suspect_segments",
                        if repair.graph_localization.suspect_segments.is_empty() {
                            ".".to_string()
                        } else {
                            repair.graph_localization.suspect_segments.join(",")
                        },
                    ),
                    ("sv_candidate_count", repair.candidate_count.to_string()),
                    (
                        "sv_evaluated_candidates",
                        repair.evaluated_candidates.to_string(),
                    ),
                    ("sv_target_reads", repair.target_reads.to_string()),
                    (
                        "sv_target_type1_reads",
                        repair.target_type1_reads.to_string(),
                    ),
                    (
                        "sv_after_low_green_window_fraction",
                        format!("{:.6}", repair.evaluation.low_green_window_fraction),
                    ),
                    ("sv_correction_count", sv_corrections.to_string()),
                    ("ordinary_max_rounds", config.max_rounds.to_string()),
                    ("effective_max_round", effective_max_round.to_string()),
                    ("hard_max_rounds", WORKFLOW_HARD_MAX_ROUNDS.to_string()),
                ],
            )?;
            write_next_round_script(
                &checkpoint_dir.join("next_round.sh"),
                config,
                case,
                round + 1,
                Some(&corrected),
                None,
            )?;
            if run_next {
                run_polish_round(config, case, round + 1, &corrected, force)?;
            }
            return Ok(Checkpoint2Status::NextRoundReady);
        }
    }

    if sv_status.as_deref() != Some("pass") {
        write_checkpoint2_status(
            &checkpoint_dir.join("checkpoint_2.status.tsv"),
            "manual_required",
            "SV support did not pass; manually inspect/correct the graph or polished sequence",
            round,
            &summary_path,
            &high_path,
            None,
            None,
        )?;
        write_next_round_script(
            &checkpoint_dir.join("next_round.sh"),
            config,
            case,
            round + 1,
            None,
            Some("SV not pass; do not run automatic SNV/InDel correction"),
        )?;
        return Ok(Checkpoint2Status::ManualRequired);
    }

    let edits = read_high_variant_edits(&high_path)?;
    let repeat_can_advance = edits.is_empty() || !config.auto_snv_indel_correction;
    if repeat_can_advance && rebuild_uses_recorded_source(config) {
        if let Some(marker) = repeat_marker.as_ref().filter(|marker| {
            marker.state == "pending"
                || (marker.state == "rebuild_ready" && !repeat_rebuild_applied)
        }) {
            // Keep the mismatch and marker intact: the native repeat resolver
            // cannot consume a recorded resolve-path contract. This also stops
            // an unapplied rebuild_ready marker left by an earlier workflow.
            let reason = "repeat pairing correction is incompatible with the recorded resolve-path source required by enabled postrepair or refine; inspect the source without running automatic repeat rebuild";
            let status_path = checkpoint_dir.join("checkpoint_2.status.tsv");
            write_checkpoint2_status(
                &status_path,
                "manual_required",
                reason,
                round,
                &summary_path,
                &high_path,
                None,
                None,
            )?;
            append_repeat_marker_metrics(
                &status_path,
                marker,
                "repeat_pairing_recorded_source_blocked",
            )?;
            write_next_round_script(
                &checkpoint_dir.join("next_round.sh"),
                config,
                case,
                round + 1,
                None,
                Some(reason),
            )?;
            return Ok(Checkpoint2Status::ManualRequired);
        }
    }
    if repeat_can_advance
        && repeat_marker
            .as_ref()
            .is_some_and(|marker| marker.state == "pending"
                || (marker.state == "rebuild_ready" && !repeat_rebuild_applied))
    {
        let marker = repeat_marker.as_mut().unwrap();
        let status_path = checkpoint_dir.join("checkpoint_2.status.tsv");
        if !config.rebuild.enabled || round >= WORKFLOW_HARD_MAX_ROUNDS {
            write_checkpoint2_status(
                &status_path,
                "manual_required",
                "repeat pairing correction requires enabled rebuild and one remaining validation round",
                round,
                &summary_path,
                &high_path,
                None,
                None,
            )?;
            append_repeat_marker_metrics(&status_path, marker, "repeat_pairing_rebuild_blocked")?;
            return Ok(Checkpoint2Status::ManualRequired);
        }

        marker.state = "rebuild_ready".to_string();
        marker.rebuild_fasta = Some(rebuild_fasta_path(case));
        marker.validation_round = Some(round + 1);
        write_repeat_pairing_marker(case, marker)?;
        write_checkpoint2_status(
            &status_path,
            "rebuild_ready",
            "repeat pairing mismatch is deferred to rebuild; constrain the triggered repeat with complete spanning-read support, then run one validate-only round",
            round,
            &summary_path,
            &high_path,
            None,
            Some(&rebuild_fasta_path(case)),
        )?;
        append_repeat_marker_metrics(&status_path, marker, "repeat_pairing_rebuild")?;
        append_checkpoint2_metrics(
            &status_path,
            &[
                ("correction_kind", "sv".to_string()),
                (
                    "sv_correction_count",
                    (sv_correction_count(case, round)? + 1).to_string(),
                ),
                ("ordinary_max_rounds", config.max_rounds.to_string()),
                ("hard_max_rounds", WORKFLOW_HARD_MAX_ROUNDS.to_string()),
            ],
        )?;
        write_repeat_rebuild_round_script(
            &checkpoint_dir.join("next_round.sh"),
            config,
            case,
            round,
            &polished_aln,
            marker,
        )?;
        if run_next {
            run_repeat_rebuild_then_validate(config, case, round, &polished_aln, marker, force)?;
        }
        return Ok(Checkpoint2Status::RebuildReady);
    }

    if edits.is_empty() {
        if repeat_rebuild_applied {
            if let Some(marker) = repeat_marker.as_mut() {
                marker.state = "validated".to_string();
                marker.validation_round = Some(round);
                write_repeat_pairing_marker(case, marker)?;
            }
        }
        let sv_corrections = sv_correction_count(case, round)?;
        let ordinary_round = round.saturating_sub(sv_corrections);
        write_checkpoint2_status(
            &checkpoint_dir.join("checkpoint_2.status.tsv"),
            "complete",
            "SV passed and no correction-worthy SNV/InDel rows remain",
            round,
            &summary_path,
            &high_path,
            None,
            None,
        )?;
        append_checkpoint2_metrics(
            &checkpoint_dir.join("checkpoint_2.status.tsv"),
            &[
                ("sv_correction_count", sv_corrections.to_string()),
                ("ordinary_round", ordinary_round.to_string()),
                ("ordinary_max_rounds", config.max_rounds.to_string()),
                ("hard_max_rounds", WORKFLOW_HARD_MAX_ROUNDS.to_string()),
            ],
        )?;
        return Ok(Checkpoint2Status::Complete);
    }

    if !config.auto_snv_indel_correction {
        if repeat_rebuild_applied {
            if let Some(marker) = repeat_marker.as_mut() {
                marker.state = "validated".to_string();
                marker.validation_round = Some(round);
                write_repeat_pairing_marker(case, marker)?;
            }
        }
        write_checkpoint2_status(
            &checkpoint_dir.join("checkpoint_2.status.tsv"),
            "complete",
            "SV passed; high-confidence SNV/InDel rows remain, but automatic correction is disabled for this workflow",
            round,
            &summary_path,
            &high_path,
            None,
            None,
        )?;
        return Ok(Checkpoint2Status::Complete);
    }

    let sv_corrections = sv_correction_count(case, round)?;
    let ordinary_round = round.saturating_sub(sv_corrections);
    if ordinary_round >= config.max_rounds || round >= WORKFLOW_HARD_MAX_ROUNDS {
        let pos_ref_alt = checkpoint_dir.join("pos_ref_alt.txt");
        write_pos_ref_alt(&pos_ref_alt, &edits)?;
        write_checkpoint2_status(
            &checkpoint_dir.join("checkpoint_2.status.tsv"),
            "manual_required",
            "SNV/InDel still fails after the ordinary max_rounds budget or the hard total of 10 rounds",
            round,
            &summary_path,
            &high_path,
            Some(&pos_ref_alt),
            None,
        )?;
        write_next_round_script(
            &checkpoint_dir.join("next_round.sh"),
            config,
            case,
            round + 1,
            None,
            Some("ordinary max_rounds or hard total round limit reached"),
        )?;
        append_checkpoint2_metrics(
            &checkpoint_dir.join("checkpoint_2.status.tsv"),
            &[
                ("correction_kind", "snv_indel".to_string()),
                ("sv_correction_count", sv_corrections.to_string()),
                ("ordinary_round", ordinary_round.to_string()),
                ("ordinary_max_rounds", config.max_rounds.to_string()),
                ("hard_max_rounds", WORKFLOW_HARD_MAX_ROUNDS.to_string()),
            ],
        )?;
        return Ok(Checkpoint2Status::ManualRequired);
    }

    let pos_ref_alt = checkpoint_dir.join("pos_ref_alt.txt");
    let corrected = corrected_fasta_path(case, round);
    if corrected.exists() && !force {
        return Err(OrgraftError::InvalidArgument(format!(
            "{} already exists; pass --force to overwrite",
            corrected.display()
        )));
    }
    write_pos_ref_alt(&pos_ref_alt, &edits)?;
    apply_edits_to_fasta(&polished_aln, &corrected, &edits)?;
    write_checkpoint2_status(
        &checkpoint_dir.join("checkpoint_2.status.tsv"),
        "next_round_ready",
        "SV passed; high-confidence SNV/InDel corrections were applied to the next round input",
        round,
        &summary_path,
        &high_path,
        Some(&pos_ref_alt),
        Some(&corrected),
    )?;
    append_checkpoint2_metrics(
        &checkpoint_dir.join("checkpoint_2.status.tsv"),
        &[
            ("correction_kind", "snv_indel".to_string()),
            ("sv_correction_count", sv_corrections.to_string()),
            ("ordinary_round", ordinary_round.to_string()),
            ("ordinary_max_rounds", config.max_rounds.to_string()),
            ("hard_max_rounds", WORKFLOW_HARD_MAX_ROUNDS.to_string()),
        ],
    )?;
    write_next_round_script(
        &checkpoint_dir.join("next_round.sh"),
        config,
        case,
        round + 1,
        Some(&corrected),
        None,
    )?;

    if run_next {
        run_polish_round(config, case, round + 1, &corrected, force)?;
    }

    Ok(Checkpoint2Status::NextRoundReady)
}

#[cfg(test)]
fn write_plan_script(
    path: &Path,
    config: &WorkflowConfig,
    case: &WorkflowCase,
    manual_message: Option<&str>,
    comment_checkpoint1: bool,
    include_recruit: bool,
) -> Result<(), OrgraftError> {
    let mut script = render_plan_script(config, case, manual_message, comment_checkpoint1, include_recruit, None)?;
    write_refine_stage(&mut script, config, &case.organelle)?;
    write_executable_text(path, &script)
}

fn render_plan_script(
    config: &WorkflowConfig,
    case: &WorkflowCase,
    manual_message: Option<&str>,
    comment_checkpoint1: bool,
    include_recruit: bool,
    shared_resolve: Option<&WorkflowCase>,
) -> Result<String, OrgraftError> {
    let mut script = String::new();
    writeln!(script, "#!/usr/bin/env bash").unwrap();
    writeln!(script, "set -euo pipefail").unwrap();
    writeln!(script).unwrap();
    writeln!(
        script,
        "# Generated by orgraft workflow for case {} ({}/{})",
        case.name, case.organelle, case.subgraph
    )
    .unwrap();
    writeln!(script, "# sample: {}", case.sample).unwrap();
    writeln!(script, "# mode: {}", config.mode).unwrap();
    writeln!(script, "# project_sample: {}", config.sample).unwrap();
    writeln!(script, "# results_dir: {}", config.results_dir.display()).unwrap();
    if let Some(message) = manual_message {
        writeln!(script, "# manual_required: {message}").unwrap();
    }
    writeln!(script).unwrap();
    write_script_runtime_header(&mut script);
    if case.postrepair_source_dir.is_some() {
        write_postrepair_source_stages(&mut script, config, case)?;
        return Ok(script);
    }

    let checkpoint1_status_path = case
        .workflow_dir
        .join("checkpoint_1")
        .join("checkpoint_1.status.tsv");
    writeln!(
        script,
        "checkpoint1_status_file={}",
        shell_quote(&checkpoint1_status_path)
    )
    .unwrap();
    writeln!(
        script,
        "checkpoint1_checked_gfa={}",
        shell_quote(&case.checked_draft_gfa)
    )
    .unwrap();
    writeln!(script, "resume_after_checkpoint1=0").unwrap();
    writeln!(
        script,
        "if [[ -s \"$checkpoint1_status_file\" && -s \"$checkpoint1_checked_gfa\" ]]; then"
    )
    .unwrap();
    writeln!(
        script,
        "  checkpoint1_status=\"$(status_value \"$checkpoint1_status_file\" status || true)\""
    )
    .unwrap();
    writeln!(
        script,
        "  if [[ \"$checkpoint1_status\" == \"checked\" ]]; then"
    )
    .unwrap();
    writeln!(script, "    resume_after_checkpoint1=1").unwrap();
    writeln!(
        script,
        "    echo \"checkpoint1 already checked; skip recruit/asm/checkpoint1 and resume at resolve\""
    )
    .unwrap();
    writeln!(script, "  fi").unwrap();
    writeln!(script, "fi").unwrap();
    writeln!(script).unwrap();
    writeln!(
        script,
        "if [[ \"$resume_after_checkpoint1\" != \"1\" ]]; then"
    )
    .unwrap();

    write_stage_header(
        &mut script,
        "  ",
        "01.recruit",
        "select organelle reads from raw HiFi reads and bait references",
    );
    if include_recruit && config.recruit.enabled {
        writeln!(script, "  {}", recruit_command(config)?).unwrap();
    } else if config.recruit.enabled {
        writeln!(
            script,
            "  # recruit handled once by the master workflow script; using reads from {}",
            case.reads.display()
        )
        .unwrap();
    } else {
        writeln!(
            script,
            "  # recruit disabled; using reads from {}",
            case.reads.display()
        )
        .unwrap();
    }

    write_stage_header(
        &mut script,
        "  ",
        "02.draft_asm",
        "assemble selected reads into the conservative draft graph",
    );
    if config.asm.enabled && shared_resolve.is_none() {
        write_checkpoint_draft_backup(&mut script, case);
        writeln!(script, "  {}", asm_command(config, case)).unwrap();
        write_checkpoint_draft_restore(&mut script, case);
    } else {
        writeln!(
            script,
            "  # asm disabled; using draft graph from {}",
            case.draft_graph.display()
        )
        .unwrap();
    }

    let checkpoint1 = format!(
        "{} workflow checkpoint1 --config {} --case {}{}",
        orgraft_shell_token(),
        shell_quote(&config.config_path),
        shell_quote_str(&case.name),
        if config.force { " --force" } else { "" }
    );
    write_stage_header(
        &mut script,
        "  ",
        "checkpoint1",
        "check graph topology and materialize checked_draft.gfa when simple",
    );
    if comment_checkpoint1 {
        writeln!(script, "  # already-run: {checkpoint1}").unwrap();
    } else {
        writeln!(script, "  {checkpoint1}").unwrap();
    }

    writeln!(
        script,
        "  checkpoint1_status=\"$(status_value \"$checkpoint1_status_file\" status)\""
    )
    .unwrap();
    writeln!(
        script,
        "  if [[ \"$checkpoint1_status\" != \"checked\" ]]; then"
    )
    .unwrap();
    writeln!(
        script,
        "    echo \"checkpoint1 status is ${{checkpoint1_status:-unknown}}; stop for manual graph editing\""
    )
    .unwrap();
    writeln!(script, "    exit 0").unwrap();
    writeln!(script, "  fi").unwrap();
    writeln!(script, "fi").unwrap();
    writeln!(script).unwrap();

    write_stage_header(
        &mut script,
        "",
        "03.resolve_gfa",
        "linearize the checked graph and prepare rotated reference handoff files",
    );
    write_case_resolve_stage(&mut script, config, case, shared_resolve);
    writeln!(script).unwrap();
    writeln!(script, "final_polished=\"\"").unwrap();
    writeln!(script, "next_round_draft=\"\"").unwrap();
    writeln!(script, "repeat_resolution_node=\"\"").unwrap();
    writeln!(script, "provisional_rebuild_fasta=\"\"").unwrap();
    writeln!(script).unwrap();
    write_stage_header(
        &mut script,
        "",
        "04.polish_checkpoint2",
        "polish, localize/validate SV evidence, and advance correction rounds",
    );

    for round in 1..=WORKFLOW_HARD_MAX_ROUNDS {
        let draft = round_draft_path(case, round);
        writeln!(script, "if [[ -z \"$final_polished\" ]]; then").unwrap();
        writeln!(script, "  # 04.polish round {round}").unwrap();
        if round == 1 {
            writeln!(script, "  round_draft={}", shell_quote(&draft)).unwrap();
        } else {
            writeln!(script, "  if [[ -n \"$next_round_draft\" ]]; then").unwrap();
            writeln!(script, "    round_draft=\"$next_round_draft\"").unwrap();
            writeln!(script, "    next_round_draft=\"\"").unwrap();
            writeln!(script, "  else").unwrap();
            writeln!(script, "    round_draft={}", shell_quote(&draft)).unwrap();
            writeln!(script, "  fi").unwrap();
        }
        writeln!(
            script,
            "  {}",
            polish_command_with_draft_arg(config, case, round, "\"${round_draft}\"", config.force,)
        )
        .unwrap();
        writeln!(script, "  # checkpoint2 round {round}").unwrap();
        writeln!(
            script,
            "  {} workflow checkpoint2 --config {} --case {} --round {}{}",
            orgraft_shell_token(),
            shell_quote(&config.config_path),
            shell_quote_str(&case.name),
            round,
            if config.force { " --force" } else { "" }
        )
        .unwrap();
        let checkpoint2_status_path = case
            .workflow_dir
            .join("checkpoint_2")
            .join(format!("round_{round}"))
            .join("checkpoint_2.status.tsv");
        writeln!(
            script,
            "  checkpoint2_status=\"$(status_value {} status)\"",
            shell_quote(&checkpoint2_status_path)
        )
        .unwrap();
        writeln!(script, "  case \"$checkpoint2_status\" in").unwrap();
        writeln!(script, "    complete)").unwrap();
        let validated_fasta = polished_aln_path(case, round);
        writeln!(
            script,
            "      final_polished={}",
            shell_quote(&validated_fasta)
        )
        .unwrap();
        if postrepair_enabled(config) {
            write_homorepair_after_checkpoint(&mut script, config, case, round)?;
        }
        writeln!(
            script,
            "      echo \"workflow complete at checkpoint2 round {round}\""
        )
        .unwrap();
        writeln!(script, "      ;;").unwrap();
        writeln!(script, "    next_round_ready)").unwrap();
        if round < WORKFLOW_HARD_MAX_ROUNDS {
            writeln!(
                script,
                "      echo \"checkpoint2 round {round} ready; continuing to round {}\"",
                round + 1
            )
            .unwrap();
            writeln!(script, "      ;;").unwrap();
        } else {
            writeln!(
                script,
                "      echo \"checkpoint2 round {round} requested another round, but the hard total of {} rounds has been reached\"",
                WORKFLOW_HARD_MAX_ROUNDS
            )
            .unwrap();
            writeln!(script, "      exit 0").unwrap();
            writeln!(script, "      ;;").unwrap();
        }
        writeln!(script, "    rebuild_ready)").unwrap();
        if config.rebuild.enabled && round < WORKFLOW_HARD_MAX_ROUNDS {
            writeln!(
                script,
                "      repeat_resolution_node=\"$(status_value {} repeat_node)\"",
                shell_quote(&checkpoint2_status_path)
            )
            .unwrap();
            writeln!(
                script,
                "      if [[ -z \"$repeat_resolution_node\" ]]; then"
            )
            .unwrap();
            writeln!(
                script,
                "        echo \"checkpoint2 rebuild_ready did not provide repeat_node\" >&2"
            )
            .unwrap();
            writeln!(script, "        exit 1").unwrap();
            writeln!(script, "      fi").unwrap();
            writeln!(
                script,
                "      {}",
                rebuild_command_with_polished_arg(
                    config,
                    case,
                    true,
                    &shell_quote(&validated_fasta),
                    Some("\"${repeat_resolution_node}\""),
                )
            )
            .unwrap();
            writeln!(
                script,
                "      provisional_rebuild_fasta={}",
                shell_quote(&rebuild_fasta_path(case))
            )
            .unwrap();
            writeln!(
                script,
                "      if [[ ! -s \"$provisional_rebuild_fasta\" ]]; then"
            )
            .unwrap();
            writeln!(
                script,
                "        echo \"repeat-constrained rebuild FASTA was not created\" >&2"
            )
            .unwrap();
            writeln!(script, "        exit 1").unwrap();
            writeln!(script, "      fi").unwrap();
            writeln!(
                script,
                "      next_round_draft=\"$provisional_rebuild_fasta\""
            )
            .unwrap();
            writeln!(
                script,
                "      echo \"checkpoint2 round {round} rebuilt repeat $repeat_resolution_node; validating in round {}\"",
                round + 1
            )
            .unwrap();
            writeln!(script, "      ;;").unwrap();
        } else {
            writeln!(
                script,
                "      echo \"checkpoint2 round {round} requested repeat rebuild, but rebuild is disabled or no validation round remains\" >&2"
            )
            .unwrap();
            writeln!(script, "      exit 1").unwrap();
            writeln!(script, "      ;;").unwrap();
        }
        writeln!(script, "    manual_required)").unwrap();
        writeln!(
            script,
            "      echo \"checkpoint2 round {round} requires manual review\""
        )
        .unwrap();
        writeln!(script, "      exit 0").unwrap();
        writeln!(script, "      ;;").unwrap();
        writeln!(script, "    *)").unwrap();
        writeln!(
            script,
            "      echo \"unknown checkpoint2 status '${{checkpoint2_status:-empty}}' at round {round}\" >&2"
        )
        .unwrap();
        writeln!(script, "      exit 1").unwrap();
        writeln!(script, "      ;;").unwrap();
        writeln!(script, "  esac").unwrap();
        writeln!(script, "fi").unwrap();
        writeln!(script).unwrap();
    }

    writeln!(script).unwrap();
    writeln!(script, "if [[ -z \"$final_polished\" ]]; then").unwrap();
    writeln!(
        script,
        "  echo \"workflow stopped before a complete final polished FASTA was selected\""
    )
    .unwrap();
    writeln!(script, "  exit 0").unwrap();
    writeln!(script, "fi").unwrap();
    writeln!(script).unwrap();
    writeln!(
        script,
        "# Rebuild after checkpoint2 complete, using the final verified graph/FASTA:"
    )
    .unwrap();
    write_stage_header(
        &mut script,
        "",
        "06.rebuild",
        "rebuild the final verified graph and compact reports",
    );
    if config.rebuild.enabled {
        writeln!(script, "if [[ -n \"$repeat_resolution_node\" ]]; then").unwrap();
        writeln!(
            script,
            "  if [[ -s \"$provisional_rebuild_fasta\" ]] && cmp -s \"$final_polished\" \"$provisional_rebuild_fasta\"; then"
        )
        .unwrap();
        writeln!(
            script,
            "    echo \"repeat-constrained rebuild already matches the final validated FASTA; skip duplicate rebuild\""
        )
        .unwrap();
        writeln!(script, "  else").unwrap();
        writeln!(
            script,
            "    {}",
            rebuild_command_with_polished_arg(
                config,
                case,
                true,
                "\"${final_polished}\"",
                Some("\"${repeat_resolution_node}\""),
            )
        )
        .unwrap();
        writeln!(script, "  fi").unwrap();
        writeln!(script, "else").unwrap();
        writeln!(
            script,
            "  {}",
            rebuild_command_with_polished_arg(
                config,
                case,
                config.force,
                "\"${final_polished}\"",
                None,
            )
        )
        .unwrap();
        writeln!(script, "fi").unwrap();
    } else {
        for line in rebuild_command_with_polished(
            config, case, config.force, "FINAL_POLISHED_FASTA",
        ).lines() {
            writeln!(script, "# {line}").unwrap();
        }
    }

    Ok(script)
}

fn write_stage_header(script: &mut String, indent: &str, stage: &str, summary: &str) {
    writeln!(script, "{indent}# {stage}: {summary}").unwrap();
}

fn checkpoint_draft_backup_path(case: &WorkflowCase) -> Option<PathBuf> {
    let file_name = case.draft_graph.file_name()?.to_str()?;
    if file_name != "graph.edited.gfa" {
        return None;
    }
    Some(
        case.workflow_dir
            .join("checkpoint_1")
            .join("draft_graph.input_backup.gfa"),
    )
}

fn write_checkpoint_draft_backup(script: &mut String, case: &WorkflowCase) {
    let Some(backup) = checkpoint_draft_backup_path(case) else {
        return;
    };
    writeln!(
        script,
        "checkpoint1_draft_graph={}",
        shell_quote(&case.draft_graph)
    )
    .unwrap();
    writeln!(script, "checkpoint1_draft_backup={}", shell_quote(&backup)).unwrap();
    writeln!(script, "if [[ -f \"$checkpoint1_draft_graph\" ]]; then").unwrap();
    writeln!(
        script,
        "  mkdir -p \"$(dirname \"$checkpoint1_draft_backup\")\""
    )
    .unwrap();
    writeln!(
        script,
        "  cp \"$checkpoint1_draft_graph\" \"$checkpoint1_draft_backup\""
    )
    .unwrap();
    writeln!(script, "fi").unwrap();
}

fn write_checkpoint_draft_restore(script: &mut String, case: &WorkflowCase) {
    if checkpoint_draft_backup_path(case).is_none() {
        return;
    }
    writeln!(script, "if [[ -f \"$checkpoint1_draft_backup\" ]]; then").unwrap();
    writeln!(
        script,
        "  mkdir -p \"$(dirname \"$checkpoint1_draft_graph\")\""
    )
    .unwrap();
    writeln!(
        script,
        "  cp \"$checkpoint1_draft_backup\" \"$checkpoint1_draft_graph\""
    )
    .unwrap();
    writeln!(script, "fi").unwrap();
}

fn shared_case_source<'a>(config: &WorkflowConfig, case: &WorkflowCase, previous: &[&'a WorkflowCase]) -> Result<Option<&'a WorkflowCase>, OrgraftError> {
    if case.postrepair_source_dir.is_some() { return Ok(None); }
    for prior in previous {
        if prior.postrepair_source_dir.is_none() && prior.organelle == case.organelle && prior.resolve_out_dir == case.resolve_out_dir {
            if prior.sample != case.sample || prior.draft_graph != case.draft_graph
                || prior.reference != case.reference || prior.pre_rotated_reference != case.pre_rotated_reference
                || (config.asm.enabled && (prior.asm_reads != case.asm_reads || prior.reads != case.reads)) {
                return Err(OrgraftError::InvalidArgument(format!("cases {} and {} use different inputs for the same resolve output; give them separate resolve_out_dir values", prior.name, case.name)));
            }
            return Ok(Some(prior));
        }
    }
    Ok(None)
}

fn write_case_resolve_stage(script: &mut String, config: &WorkflowConfig, case: &WorkflowCase, shared: Option<&WorkflowCase>) {
    if let Some(prior) = shared {
        writeln!(script, "# Reuse the shared resolve only after both checked inputs agree.").unwrap();
        writeln!(script, "cmp -s {} {} || {{ echo 'Checked component inputs differ; review the shared resolve configuration' >&2; exit 1; }}", shell_quote(&prior.checked_draft_gfa), shell_quote(&case.checked_draft_gfa)).unwrap();
        writeln!(script, "[[ -s {} && -s {} ]] || {{ echo 'Shared resolve is incomplete' >&2; exit 1; }}", shell_quote(&resolve_details_path(case)), shell_quote(&rebuild_source_gfa(config, case))).unwrap();
    } else {
        writeln!(script, "{}", resolve_command(config, case, config.force)).unwrap();
    }
}

fn disabled_plan_script(config: &WorkflowConfig) -> String {
    let message = format!("No runnable workflow cases: every configured case has enabled=false in {}. Review the disabled cases before running this route.", config.config_path.display());
    format!("#!/usr/bin/env bash\nset -euo pipefail\nprintf '%s\\n' {} >&2\nexit 1\n", shell_quote_str(&message))
}

fn render_all_cases_plan_script(config: &WorkflowConfig) -> Result<String, OrgraftError> {
    let enabled = enabled_cases(config);
    if enabled.is_empty() {
        return Ok(disabled_plan_script(config));
    }
    let mut script = String::from("#!/usr/bin/env bash\nset -euo pipefail\n\n# Generated by orgraft workflow for all configured cases\n");
    writeln!(script, "# sample: {}", config.sample).unwrap();
    write_script_runtime_header(&mut script);
    if enabled.iter().any(|case| case.postrepair_source_dir.is_none()) {
        writeln!(script, "# 01.recruit is shared across enabled cases and runs at most once here.").unwrap();
    }
    write_master_recruit_stage(&mut script, config, &enabled)?;
    let mut completed_sources: Vec<&WorkflowCase> = Vec::new();
    for case in &enabled {
        let shared = shared_case_source(config, case, &completed_sources)?;
        writeln!(script, "\n# Case {} ({}/{})\n(", case.name, case.organelle, case.subgraph).unwrap();
        // A subshell isolates case variables and preserves each manual checkpoint stop.
        // Every command is written here; no generated child script is required.
        script.push_str(&render_plan_script(config, case, None, false, false, shared)?);
        writeln!(script, ")").unwrap();
        completed_sources.push(case);
    }
    let organelles: BTreeSet<_> = enabled.iter().map(|case| case.organelle.as_str()).collect();
    for organelle in organelles { write_refine_stage(&mut script, config, organelle)?; }
    Ok(script)
}

fn completed_round(case: &WorkflowCase) -> Result<Option<usize>, OrgraftError> {
    for round in (1..=WORKFLOW_HARD_MAX_ROUNDS).rev() {
        let status = case.workflow_dir.join("checkpoint_2").join(format!("round_{round}")).join("checkpoint_2.status.tsv");
        if status.is_file() && read_metric_value(&status, "status")?.as_deref() == Some("complete") {
            return Ok(Some(round));
        }
    }
    Ok(None)
}

fn render_manual_plan_script(config: &WorkflowConfig, cases: &[&WorkflowCase]) -> Result<String, OrgraftError> {
    if cases.is_empty() { return Ok(disabled_plan_script(config)); }
    let mut script = String::from("#!/usr/bin/env bash\nset -euo pipefail\n# Each OrgRAFT invocation below is a standalone command.\n# This script shows one validation round per case. If a checkpoint requests\n# another round or manual review, use run_<sample>.sh to resume after review.\n");
    write_script_runtime_header(&mut script);
    write_master_recruit_stage(&mut script, config, cases)?;
    let mut completed_sources: Vec<&WorkflowCase> = Vec::new();
    for case in cases {
        if case.postrepair_source_dir.is_some() {
            writeln!(script, "\n# {}: read-only ordinary validation input.", case.name).unwrap();
            write_postrepair_source_stages(&mut script, config, case)?;
            continue;
        }
        let shared = shared_case_source(config, case, &completed_sources)?;
        let round = completed_round(case)?.unwrap_or(1);
        writeln!(script, "\n# {}: inspect each checkpoint before continuing.", case.name).unwrap();
        if config.asm.enabled && shared.is_none() { writeln!(script, "{}", asm_command(config, case)).unwrap(); }
        writeln!(script, "{} workflow checkpoint1 --config {} --case {}{}", orgraft_shell_token(), shell_quote(&config.config_path), shell_quote_str(&case.name), if config.force { " --force" } else { "" }).unwrap();
        let check1 = case.workflow_dir.join("checkpoint_1/checkpoint_1.status.tsv");
        writeln!(script, "[[ $(status_value {} status) == checked ]] || {{ echo 'Stop for checkpoint 1 review'; exit 0; }}", shell_quote(&check1)).unwrap();
        write_case_resolve_stage(&mut script, config, case, shared);
        writeln!(script, "{}", polish_command_with_draft_arg(config, case, round, &shell_quote(&round_draft_path(case, round)), config.force)).unwrap();
        writeln!(script, "{} workflow checkpoint2 --config {} --case {} --round {round}{}", orgraft_shell_token(), shell_quote(&config.config_path), shell_quote_str(&case.name), if config.force { " --force" } else { "" }).unwrap();
        let check2 = case.workflow_dir.join("checkpoint_2").join(format!("round_{round}/checkpoint_2.status.tsv"));
        writeln!(script, "[[ $(status_value {} status) == complete ]] || {{ echo 'Stop for checkpoint 2 review or the next round'; exit 0; }}", shell_quote(&check2)).unwrap();
        writeln!(script, "final_polished={}", shell_quote(&polished_aln_path(case, round))).unwrap();
        if postrepair_enabled(config) { write_homorepair_after_checkpoint(&mut script, config, case, round)?; }
        if config.rebuild.enabled { writeln!(script, "{}", rebuild_command_with_polished_arg(config, case, config.force, "\"${final_polished}\"", None)).unwrap(); }
        completed_sources.push(case);
    }
    let organelles: BTreeSet<_> = cases.iter().map(|case| case.organelle.as_str()).collect();
    for organelle in organelles { write_refine_stage(&mut script, config, organelle)?; }
    Ok(script)
}

fn write_master_recruit_stage(
    script: &mut String,
    config: &WorkflowConfig,
    cases: &[&WorkflowCase],
) -> Result<(), OrgraftError> {
    if cases.iter().all(|case| case.postrepair_source_dir.is_some()) { return Ok(()); }
    write_stage_header(
        script,
        "",
        "01.recruit",
        "select organelle reads once for all enabled workflow cases",
    );
    if config.recruit.enabled {
        let command = recruit_command(config)?;
        if config.force {
            writeln!(script, "{command}").unwrap();
        } else {
            writeln!(script, "recruit_ready=1").unwrap();
            for case in cases {
                writeln!(
                    script,
                    "[[ -s {} ]] || recruit_ready=0",
                    shell_quote(&case.reads)
                )
                .unwrap();
            }
            writeln!(script, "if [[ \"$recruit_ready\" == \"1\" ]]; then").unwrap();
            writeln!(
                script,
                "  echo \"recruit outputs already exist; skip global recruit\""
            )
            .unwrap();
            writeln!(script, "else").unwrap();
            writeln!(script, "  {command}").unwrap();
            writeln!(script, "fi").unwrap();
        }
    } else {
        writeln!(
            script,
            "# recruit disabled; each case uses its configured reads path"
        )
        .unwrap();
    }
    writeln!(script).unwrap();
    Ok(())
}

fn write_next_round_script(
    path: &Path,
    config: &WorkflowConfig,
    case: &WorkflowCase,
    round: usize,
    draft: Option<&Path>,
    manual_message: Option<&str>,
) -> Result<(), OrgraftError> {
    let mut script = String::new();
    writeln!(script, "#!/usr/bin/env bash").unwrap();
    writeln!(script, "set -euo pipefail").unwrap();
    writeln!(script).unwrap();
    write_script_runtime_header(&mut script);
    if let Some(message) = manual_message {
        writeln!(script, "# manual_required: {message}").unwrap();
        writeln!(script, "printf '%s\\n' {} >&2", shell_quote_str(message)).unwrap();
        writeln!(script, "exit 1").unwrap();
    } else if let Some(draft) = draft {
        writeln!(
            script,
            "{}",
            polish_command(config, case, round, draft, config.force)
        )
        .unwrap();
        writeln!(
            script,
            "{} workflow checkpoint2 --config {} --case {} --round {}{}",
            orgraft_shell_token(),
            shell_quote(&config.config_path),
            shell_quote_str(&case.name),
            round,
            if config.force { " --force" } else { "" }
        )
        .unwrap();
    }
    write_executable_text(path, &script)
}

fn write_repeat_rebuild_round_script(
    path: &Path,
    config: &WorkflowConfig,
    case: &WorkflowCase,
    source_round: usize,
    polished_fasta: &Path,
    marker: &RepeatPairingMarker,
) -> Result<(), OrgraftError> {
    let mut script = String::new();
    writeln!(script, "#!/usr/bin/env bash").unwrap();
    writeln!(script, "set -euo pipefail").unwrap();
    writeln!(script).unwrap();
    write_script_runtime_header(&mut script);
    writeln!(
        script,
        "{}",
        rebuild_command_with_polished_and_repeat(
            config,
            case,
            config.force,
            &polished_fasta.display().to_string(),
            Some(&marker.repeat_node),
        )
    )
    .unwrap();
    writeln!(
        script,
        "{}",
        polish_command(
            config,
            case,
            source_round + 1,
            &rebuild_fasta_path(case),
            config.force,
        )
    )
    .unwrap();
    writeln!(
        script,
        "{} workflow checkpoint2 --config {} --case {} --round {}{}",
        orgraft_shell_token(),
        shell_quote(&config.config_path),
        shell_quote_str(&case.name),
        source_round + 1,
        if config.force { " --force" } else { "" }
    )
    .unwrap();
    write_executable_text(path, &script)
}

fn run_repeat_rebuild_then_validate(
    config: &WorkflowConfig,
    case: &WorkflowCase,
    source_round: usize,
    polished_fasta: &Path,
    marker: &RepeatPairingMarker,
    force: bool,
) -> Result<(), OrgraftError> {
    let rebuild_fasta = rebuild_fasta_path(case);
    let already_built = !force
        && rebuild_fasta.is_file()
        && rebuild_has_repeat_resolution(case, &marker.repeat_node)?;
    if !already_built {
        commands::rebuild::run(&rebuild_args_with_polished(
            config,
            case,
            force || case.rebuild_out_dir.exists(),
            Some(polished_fasta),
            Some(&marker.repeat_node),
        )?)?;
    }
    if !rebuild_has_repeat_resolution(case, &marker.repeat_node)? {
        return Err(OrgraftError::InvalidArgument(format!(
            "rebuild did not record constrained resolution for repeat `{}`",
            marker.repeat_node
        )));
    }
    let mut applied = marker.clone();
    applied.state = "rebuild_ready".to_string();
    applied.rebuild_fasta = Some(rebuild_fasta.clone());
    applied.validation_round = Some(source_round + 1);
    write_repeat_pairing_marker(case, &applied)?;
    run_polish_round(config, case, source_round + 1, &rebuild_fasta, force)
}

fn write_script_runtime_header(script: &mut String) {
    writeln!(script, "if [[ -z \"${{ORGRAFT_BIN:-}}\" ]]; then").unwrap();
    writeln!(script, "  ORGRAFT_BIN={}", default_orgraft_bin()).unwrap();
    writeln!(script, "fi").unwrap();
    writeln!(script).unwrap();
    writeln!(script, "status_value() {{").unwrap();
    writeln!(script, "  local file=\"$1\"").unwrap();
    writeln!(script, "  local key=\"$2\"").unwrap();
    writeln!(
        script,
        "  awk -F '\\t' -v key=\"$key\" '$1 == key {{ print $2; exit }}' \"$file\""
    )
    .unwrap();
    writeln!(script, "}}").unwrap();
    writeln!(script).unwrap();
}

fn orgraft_shell_token() -> &'static str {
    "\"${ORGRAFT_BIN}\""
}

fn recruit_threads(config: &WorkflowConfig) -> usize {
    config.recruit.threads.unwrap_or(config.threads)
}

fn asm_threads(config: &WorkflowConfig) -> usize {
    config.asm.threads.unwrap_or(config.threads)
}

fn asm_image_reference_fasta<'a>(
    config: &'a WorkflowConfig,
    case: &'a WorkflowCase,
) -> Option<&'a Path> {
    config
        .asm
        .image_reference_fasta
        .as_deref()
        .or(case.reference.as_deref())
        .or(case.pre_rotated_reference.as_deref())
}

fn polish_threads(config: &WorkflowConfig) -> usize {
    config.polish.threads.unwrap_or(config.threads)
}

fn recruit_args(config: &WorkflowConfig) -> Result<Vec<String>, OrgraftError> {
    let recruit = &config.recruit;
    let reads = recruit.reads.as_ref().ok_or_else(|| {
        OrgraftError::InvalidArgument(
            "commands.recruit.enabled=true requires commands.recruit.raw_reads".to_string(),
        )
    })?;
    if recruit.baits.is_empty() {
        return Err(OrgraftError::InvalidArgument(
            "commands.recruit.enabled=true requires commands.recruit.baits".to_string(),
        ));
    }

    let mut args = vec![
        "--reads".to_string(),
        reads.display().to_string(),
        "--out-dir".to_string(),
        recruit.out_dir.display().to_string(),
        "--threads".to_string(),
        recruit_threads(config).to_string(),
    ];
    for (label, path) in &recruit.baits {
        match label.as_str() {
            "mito" => {
                args.push("--mito".to_string());
                args.push(path.display().to_string());
            }
            "plastid" | "plasti" => {
                args.push("--plastid".to_string());
                args.push(path.display().to_string());
            }
            _ => {
                args.push("--bait".to_string());
                args.push(format!("{label}={}", path.display()));
            }
        }
    }
    push_raw_string_option(&mut args, "--prefix", recruit.prefix.as_deref());
    push_raw_string_option(&mut args, "--bait-format", recruit.bait_format.as_deref());
    push_raw_string_option(&mut args, "--gfa-split", recruit.gfa_split.as_deref());
    push_raw_flag(&mut args, "--rename-bait", recruit.rename_bait);
    push_raw_flag(&mut args, "--write-id-map", recruit.write_id_map);
    push_raw_string_option(&mut args, "--split-output", recruit.split_output.as_deref());
    if let Some(gzip_output) = recruit.gzip_output {
        push_raw_string_option(
            &mut args,
            "--gzip-output",
            Some(if gzip_output { "on" } else { "off" }),
        );
    }
    push_raw_string_option(&mut args, "--minimap2", recruit.minimap2.as_deref());
    push_raw_string_option(&mut args, "--align-mode", recruit.align_mode.as_deref());
    push_raw_string_option(&mut args, "--platform", recruit.platform.as_deref());
    push_raw_string_option(&mut args, "--preset", recruit.preset.as_deref());
    push_raw_display_option(&mut args, "--min-mapq", recruit.min_mapq);
    push_raw_display_option(&mut args, "--min-aln-len", recruit.min_aln_len);
    push_raw_path_option(&mut args, "--sam", recruit.sam.as_deref());
    for max_reads in &recruit.max_reads {
        push_raw_string_option(&mut args, "--max-reads", Some(max_reads));
    }
    push_raw_display_option(&mut args, "--random-seed", recruit.random_seed);
    push_raw_flag(&mut args, "--write-sampled-ids", recruit.write_sampled_ids);
    push_raw_string_option(&mut args, "--read-stats", recruit.read_stats.as_deref());
    push_raw_flag(
        &mut args,
        "--write-read-classification",
        recruit.write_read_classification,
    );
    push_raw_flag(
        &mut args,
        "--write-bait-partitions",
        recruit.write_bait_partitions,
    );
    push_raw_string_option(&mut args, "--gzip-tool", recruit.gzip_tool.as_deref());
    push_raw_string_option(&mut args, "--mode", recruit.mode.as_deref());
    push_raw_display_option(&mut args, "--iterations", recruit.iterations);
    args.extend(recruit.extra_args.iter().cloned());
    Ok(args)
}

fn asm_args(config: &WorkflowConfig, case: &WorkflowCase) -> Vec<String> {
    let reads = case.asm_reads.as_ref().cloned().unwrap_or_else(|| {
        if config.recruit.enabled {
            config
                .recruit
                .out_dir
                .join(format!("{}.fastq.gz", case.organelle))
        } else {
            case.reads.clone()
        }
    });
    let asm = &config.asm;
    let mut args = vec![
        "--reads".to_string(),
        reads.display().to_string(),
        "--organelle".to_string(),
        case.organelle.clone(),
        "--soft-paths".to_string(),
        config.soft_paths.display().to_string(),
        "--out-dir".to_string(),
        asm.out_dir.display().to_string(),
        "--threads".to_string(),
        asm_threads(config).to_string(),
    ];
    push_raw_string_option(&mut args, "--profile", asm.profile.as_deref());
    push_raw_flag(&mut args, "--stable", asm.stable);
    push_raw_display_option(&mut args, "--min-graph-coverage", asm.min_graph_coverage);
    push_raw_display_option(&mut args, "--branch-ratio", asm.min_branch_ratio);
    push_raw_display_option(&mut args, "--tip-len", asm.min_tip_len);
    push_raw_display_option(&mut args, "--link-support", asm.min_link_support);
    push_raw_display_option(&mut args, "--min-link-ratio", asm.min_link_ratio);
    push_raw_string_option(&mut args, "--subsets", asm.subsets.as_deref());
    push_raw_path_option(
        &mut args,
        "--image-reference-fasta",
        asm_image_reference_fasta(config, case),
    );
    push_raw_flag(&mut args, "--keep-debug-files", asm.keep_debug_files);
    push_raw_flag(&mut args, "--force", config.force);
    args.extend(asm.extra_args.iter().cloned());
    args
}

fn rebuild_args_with_polished(
    config: &WorkflowConfig,
    case: &WorkflowCase,
    force: bool,
    polished_override: Option<&Path>,
    repeat_resolution_node: Option<&str>,
) -> Result<Vec<String>, OrgraftError> {
    let rebuild = &config.rebuild;
    let polished = polished_override
        .or(case.rebuild_polished_fasta.as_deref())
        .or(rebuild.polished_fasta.as_deref())
        .ok_or_else(|| {
            OrgraftError::InvalidArgument(
                "commands.rebuild.enabled=true requires final polished FASTA from checkpoint2 or rebuild_polished_fasta/commands.rebuild.polished_fasta".to_string(),
            )
        })?;
    let edited_gfa = rebuild_source_gfa(config, case);
    let default_reference = polish_reference_path(case);
    let image_reference_fasta = case
        .image_reference_fasta
        .as_ref()
        .or(rebuild.image_reference_fasta.as_ref())
        .or_else(|| config.refine.as_ref().filter(|r| r.enabled).map(|_| &default_reference));

    let mut args = vec![
        "--organelle".to_string(),
        case.organelle.clone(),
        "--subgraph".to_string(),
        case.subgraph.clone(),
        "--edited-gfa".to_string(),
        edited_gfa.display().to_string(),
        "--polished-fasta".to_string(),
        polished.display().to_string(),
        "--reads".to_string(),
        case.reads.display().to_string(),
        "--soft-paths".to_string(),
        config.soft_paths.display().to_string(),
        "--out-dir".to_string(),
        case.rebuild_out_dir.display().to_string(),
    ];
    if rebuild_uses_recorded_source(config) {
        args.extend(["--checked-gfa".into(), case.checked_draft_gfa.display().to_string(),
            "--resolve-details".into(), resolve_details_path(case).display().to_string()]);
    }
    if repeat_resolution_node.is_none() {
        args.extend(postrepair_rebuild_args(config, case, polished)?);
    }
    push_raw_display_option(&mut args, "--threads", rebuild.threads);
    push_raw_path_option(
        &mut args,
        "--image-reference-fasta",
        image_reference_fasta.map(PathBuf::as_path),
    );
    push_raw_path_option(
        &mut args,
        "--merged-gfa-template",
        rebuild.merged_gfa_template.as_deref(),
    );
    push_raw_path_option(&mut args, "--minimap2", rebuild.minimap2.as_deref());
    push_raw_path_option(&mut args, "--blastn", rebuild.blastn.as_deref());
    push_raw_string_option(
        &mut args,
        "--repeat-resolution-node",
        repeat_resolution_node,
    );
    push_raw_flag(&mut args, "--keep-debug", rebuild.keep_debug);
    push_raw_flag(&mut args, "--force", force);
    args.extend(rebuild.extra_args.iter().cloned());
    Ok(args)
}

fn recruit_command(config: &WorkflowConfig) -> Result<String, OrgraftError> {
    let recruit = &config.recruit;
    let reads = recruit.reads.as_ref().ok_or_else(|| {
        OrgraftError::InvalidArgument(
            "commands.recruit.enabled=true requires commands.recruit.raw_reads".to_string(),
        )
    })?;
    if recruit.baits.is_empty() {
        return Err(OrgraftError::InvalidArgument(
            "commands.recruit.enabled=true requires commands.recruit.baits".to_string(),
        ));
    }

    let mut args = vec![
        orgraft_shell_token().to_string(),
        "recruit".to_string(),
        "--reads".to_string(),
        shell_quote(reads),
        "--out-dir".to_string(),
        shell_quote(&recruit.out_dir),
        "--threads".to_string(),
        recruit_threads(config).to_string(),
    ];
    for (label, path) in &recruit.baits {
        match label.as_str() {
            "mito" => {
                args.push("--mito".to_string());
                args.push(shell_quote(path));
            }
            "plastid" | "plasti" => {
                args.push("--plastid".to_string());
                args.push(shell_quote(path));
            }
            _ => {
                args.push("--bait".to_string());
                args.push(shell_quote_str(&format!("{label}={}", path.display())));
            }
        }
    }
    push_string_option(&mut args, "--prefix", recruit.prefix.as_deref());
    push_string_option(&mut args, "--bait-format", recruit.bait_format.as_deref());
    push_string_option(&mut args, "--gfa-split", recruit.gfa_split.as_deref());
    push_flag(&mut args, "--rename-bait", recruit.rename_bait);
    push_flag(&mut args, "--write-id-map", recruit.write_id_map);
    push_string_option(&mut args, "--split-output", recruit.split_output.as_deref());
    if let Some(gzip_output) = recruit.gzip_output {
        push_string_option(
            &mut args,
            "--gzip-output",
            Some(if gzip_output { "on" } else { "off" }),
        );
    }
    push_string_option(&mut args, "--minimap2", recruit.minimap2.as_deref());
    push_string_option(&mut args, "--align-mode", recruit.align_mode.as_deref());
    push_string_option(&mut args, "--platform", recruit.platform.as_deref());
    push_string_option(&mut args, "--preset", recruit.preset.as_deref());
    push_display_option(&mut args, "--min-mapq", recruit.min_mapq);
    push_display_option(&mut args, "--min-aln-len", recruit.min_aln_len);
    push_path_option(&mut args, "--sam", recruit.sam.as_deref());
    for max_reads in &recruit.max_reads {
        push_string_option(&mut args, "--max-reads", Some(max_reads));
    }
    push_display_option(&mut args, "--random-seed", recruit.random_seed);
    push_flag(&mut args, "--write-sampled-ids", recruit.write_sampled_ids);
    push_string_option(&mut args, "--read-stats", recruit.read_stats.as_deref());
    push_flag(
        &mut args,
        "--write-read-classification",
        recruit.write_read_classification,
    );
    push_flag(
        &mut args,
        "--write-bait-partitions",
        recruit.write_bait_partitions,
    );
    push_string_option(&mut args, "--gzip-tool", recruit.gzip_tool.as_deref());
    push_string_option(&mut args, "--mode", recruit.mode.as_deref());
    push_display_option(&mut args, "--iterations", recruit.iterations);
    args.extend(recruit.extra_args.iter().cloned());
    Ok(args.join(" "))
}

fn asm_command(config: &WorkflowConfig, case: &WorkflowCase) -> String {
    let reads = case.asm_reads.as_ref().cloned().unwrap_or_else(|| {
        if config.recruit.enabled {
            config
                .recruit
                .out_dir
                .join(format!("{}.fastq.gz", case.organelle))
        } else {
            case.reads.clone()
        }
    });
    let asm = &config.asm;
    let mut args = vec![
        orgraft_shell_token().to_string(),
        "asm".to_string(),
        "--reads".to_string(),
        shell_quote(&reads),
        "--organelle".to_string(),
        shell_quote_str(&case.organelle),
        "--soft-paths".to_string(),
        shell_quote(&config.soft_paths),
        "--out-dir".to_string(),
        shell_quote(&asm.out_dir),
        "--threads".to_string(),
        asm_threads(config).to_string(),
    ];
    push_string_option(&mut args, "--profile", asm.profile.as_deref());
    push_flag(&mut args, "--stable", asm.stable);
    push_display_option(&mut args, "--min-graph-coverage", asm.min_graph_coverage);
    push_display_option(&mut args, "--branch-ratio", asm.min_branch_ratio);
    push_display_option(&mut args, "--tip-len", asm.min_tip_len);
    push_display_option(&mut args, "--link-support", asm.min_link_support);
    push_display_option(&mut args, "--min-link-ratio", asm.min_link_ratio);
    push_string_option(&mut args, "--subsets", asm.subsets.as_deref());
    push_path_option(
        &mut args,
        "--image-reference-fasta",
        asm_image_reference_fasta(config, case),
    );
    push_flag(&mut args, "--keep-debug-files", asm.keep_debug_files);
    push_flag(&mut args, "--force", config.force);
    args.extend(asm.extra_args.iter().cloned());
    args.join(" ")
}

fn rebuild_command_with_polished(
    config: &WorkflowConfig,
    case: &WorkflowCase,
    force: bool,
    polished_fasta: &str,
) -> String {
    rebuild_command_with_polished_arg(config, case, force, &shell_quote_str(polished_fasta), None)
}

fn rebuild_command_with_polished_and_repeat(
    config: &WorkflowConfig,
    case: &WorkflowCase,
    force: bool,
    polished_fasta: &str,
    repeat_resolution_node: Option<&str>,
) -> String {
    rebuild_command_with_polished_arg(
        config,
        case,
        force,
        &shell_quote_str(polished_fasta),
        repeat_resolution_node.map(shell_quote_str).as_deref(),
    )
}

fn rebuild_command_with_polished_arg(
    config: &WorkflowConfig,
    case: &WorkflowCase,
    force: bool,
    polished_fasta_arg: &str,
    repeat_resolution_node_arg: Option<&str>,
) -> String {
    let rebuild = &config.rebuild;
    let edited_gfa = rebuild_source_gfa(config, case);
    let out_dir = if case.rebuild_out_dir.as_os_str().is_empty() {
        &rebuild.out_dir
    } else {
        &case.rebuild_out_dir
    };
    let default_reference = polish_reference_path(case);
    let image_reference_fasta = case
        .image_reference_fasta
        .as_ref()
        .or(rebuild.image_reference_fasta.as_ref())
        .or_else(|| config.refine.as_ref().filter(|r| r.enabled).map(|_| &default_reference));

    let mut args = vec![
        orgraft_shell_token().to_string(),
        "rebuild".to_string(),
        "--organelle".to_string(),
        shell_quote_str(&case.organelle),
        "--subgraph".to_string(),
        shell_quote_str(&case.subgraph),
        "--edited-gfa".to_string(),
        shell_quote(&edited_gfa),
        "--polished-fasta".to_string(),
        polished_fasta_arg.to_string(),
        "--reads".to_string(),
        shell_quote(&case.reads),
        "--soft-paths".to_string(),
        shell_quote(&config.soft_paths),
        "--out-dir".to_string(),
        shell_quote(out_dir),
    ];
    if rebuild_uses_recorded_source(config) {
        args.extend(["--checked-gfa".into(), shell_quote(&case.checked_draft_gfa),
            "--resolve-details".into(), shell_quote(&resolve_details_path(case))]);
    }
    push_display_option(&mut args, "--threads", rebuild.threads);
    push_path_option(
        &mut args,
        "--image-reference-fasta",
        image_reference_fasta.map(PathBuf::as_path),
    );
    push_path_option(
        &mut args,
        "--merged-gfa-template",
        rebuild.merged_gfa_template.as_deref(),
    );
    push_path_option(&mut args, "--minimap2", rebuild.minimap2.as_deref());
    push_path_option(&mut args, "--blastn", rebuild.blastn.as_deref());
    if let Some(value) = repeat_resolution_node_arg {
        args.push("--repeat-resolution-node".to_string());
        args.push(value.to_string());
    }
    push_flag(&mut args, "--keep-debug", rebuild.keep_debug);
    push_flag(&mut args, "--force", force);
    args.extend(rebuild.extra_args.iter().cloned());
    if postrepair_enabled(config) && case.postrepair_out_dir.is_some() && repeat_resolution_node_arg.is_none() {
        args.push("${joint_ir_args[@]+\"${joint_ir_args[@]}\"}".into());
        format!("(\n  joint_ir_args=()\n  repair_args=\"$({} workflow postrepair-check --config {} --case {} --rebuild-args --polished-fasta {})\" || exit $?\n  if [[ -n \"$repair_args\" ]]; then\n    while IFS= read -r arg; do joint_ir_args+=(\"$arg\"); done <<< \"$repair_args\"\n  fi\n  {}\n)",
            orgraft_shell_token(), shell_quote(&config.config_path), shell_quote_str(&case.name),
            polished_fasta_arg, args.join(" "))
    } else {
        args.join(" ")
    }
}

fn append_resolve_raw_args(args: &mut Vec<String>, resolve: &WorkflowResolveConfig) {
    push_raw_string_option(
        args,
        "--gfa-editor-mode",
        resolve.gfa_editor_mode.as_deref(),
    );
    push_raw_display_option(args, "--max-states", resolve.max_states);
    push_raw_display_option(args, "--max-candidates", resolve.max_candidates);
    args.extend(resolve.extra_args.iter().cloned());
}

fn append_resolve_shell_args(args: &mut Vec<String>, resolve: &WorkflowResolveConfig) {
    push_string_option(
        args,
        "--gfa-editor-mode",
        resolve.gfa_editor_mode.as_deref(),
    );
    push_display_option(args, "--max-states", resolve.max_states);
    push_display_option(args, "--max-candidates", resolve.max_candidates);
    args.extend(resolve.extra_args.iter().cloned());
}

fn append_polish_raw_args(args: &mut Vec<String>, polish: &WorkflowPolishConfig) {
    push_raw_string_option(
        args,
        "--per-read-variant-calls",
        polish.per_read_variant_calls.map(on_off),
    );
    push_raw_string_option(
        args,
        "--snv-indel-overlap-policy",
        polish.snv_indel_overlap_policy.as_deref(),
    );
    push_raw_string_option(args, "--plot-range", polish.plot_range.as_deref());
    push_raw_display_option(args, "--plot-dpi", polish.plot_dpi);
    push_raw_string_option(
        args,
        "--plot-output-format",
        polish.plot_output_format.as_deref(),
    );
    push_raw_string_option(
        args,
        "--coverage-plot-rasterize",
        polish.coverage_plot_rasterize.map(on_off),
    );
    push_raw_string_option(
        args,
        "--snv-indel-plot-rasterize",
        polish.snv_indel_plot_rasterize.map(on_off),
    );
    push_raw_string_option(
        args,
        "--sv-plot-highlight-subgroups",
        polish.sv_plot_highlight_subgroups.as_deref(),
    );
    push_raw_path_option(
        args,
        "--sv-plot-highlight-read-ids",
        polish.sv_plot_highlight_read_ids.as_deref(),
    );
    push_raw_display_option(
        args,
        "--sv-plot-highlight-min-fraction",
        polish.sv_plot_highlight_min_fraction,
    );
    push_raw_display_option(
        args,
        "--sv-plot-highlight-min-reads",
        polish.sv_plot_highlight_min_reads,
    );
    push_raw_string_option(
        args,
        "--snv-indel-plot-low-confidence",
        polish.snv_indel_plot_low_confidence.as_deref(),
    );
    push_raw_display_option(
        args,
        "--snv-indel-plot-low-min-reads",
        polish.snv_indel_plot_low_min_reads,
    );
    push_raw_display_option(
        args,
        "--snv-indel-plot-low-min-fraction",
        polish.snv_indel_plot_low_min_fraction,
    );
    push_raw_display_option(
        args,
        "--snv-indel-plot-high-risk-fraction",
        polish.snv_indel_plot_high_risk_fraction,
    );
    args.extend(polish.extra_args.iter().cloned());
}

fn append_polish_shell_args(args: &mut Vec<String>, polish: &WorkflowPolishConfig) {
    push_string_option(
        args,
        "--per-read-variant-calls",
        polish.per_read_variant_calls.map(on_off),
    );
    push_string_option(
        args,
        "--snv-indel-overlap-policy",
        polish.snv_indel_overlap_policy.as_deref(),
    );
    push_string_option(args, "--plot-range", polish.plot_range.as_deref());
    push_display_option(args, "--plot-dpi", polish.plot_dpi);
    push_string_option(
        args,
        "--plot-output-format",
        polish.plot_output_format.as_deref(),
    );
    push_string_option(
        args,
        "--coverage-plot-rasterize",
        polish.coverage_plot_rasterize.map(on_off),
    );
    push_string_option(
        args,
        "--snv-indel-plot-rasterize",
        polish.snv_indel_plot_rasterize.map(on_off),
    );
    push_string_option(
        args,
        "--sv-plot-highlight-subgroups",
        polish.sv_plot_highlight_subgroups.as_deref(),
    );
    push_path_option(
        args,
        "--sv-plot-highlight-read-ids",
        polish.sv_plot_highlight_read_ids.as_deref(),
    );
    push_display_option(
        args,
        "--sv-plot-highlight-min-fraction",
        polish.sv_plot_highlight_min_fraction,
    );
    push_display_option(
        args,
        "--sv-plot-highlight-min-reads",
        polish.sv_plot_highlight_min_reads,
    );
    push_string_option(
        args,
        "--snv-indel-plot-low-confidence",
        polish.snv_indel_plot_low_confidence.as_deref(),
    );
    push_display_option(
        args,
        "--snv-indel-plot-low-min-reads",
        polish.snv_indel_plot_low_min_reads,
    );
    push_display_option(
        args,
        "--snv-indel-plot-low-min-fraction",
        polish.snv_indel_plot_low_min_fraction,
    );
    push_display_option(
        args,
        "--snv-indel-plot-high-risk-fraction",
        polish.snv_indel_plot_high_risk_fraction,
    );
    args.extend(polish.extra_args.iter().cloned());
}

fn on_off(value: bool) -> &'static str {
    if value {
        "on"
    } else {
        "off"
    }
}

fn run_resolve_for_case(
    config: &WorkflowConfig,
    case: &WorkflowCase,
    force: bool,
) -> Result<(), OrgraftError> {
    commands::resolve::run(&resolve_args(config, case, force)?)
}

fn resolve_args(
    config: &WorkflowConfig,
    case: &WorkflowCase,
    force: bool,
) -> Result<Vec<String>, OrgraftError> {
    let mut args = vec![
        "--checked-draft-gfa".to_string(),
        case.checked_draft_gfa.display().to_string(),
        "--soft-paths".to_string(),
        config.soft_paths.display().to_string(),
        "--out-dir".to_string(),
        case.resolve_out_dir.display().to_string(),
        "--organelle".to_string(),
        case.organelle.clone(),
    ];
    if let Some(pre_rotated_reference) = &case.pre_rotated_reference {
        args.push("--pre-rotated-reference".to_string());
        args.push(pre_rotated_reference.display().to_string());
    } else if let Some(reference) = &case.reference {
        args.push("--reference".to_string());
        args.push(reference.display().to_string());
    } else {
        return Err(OrgraftError::InvalidArgument(format!(
            "workflow case `{}` needs reference or pre_rotated_reference",
            case.name
        )));
    }
    if force {
        args.push("--force".to_string());
    }
    append_resolve_raw_args(&mut args, &config.resolve);
    Ok(args)
}

fn run_polish_round(
    config: &WorkflowConfig,
    case: &WorkflowCase,
    round: usize,
    draft: &Path,
    force: bool,
) -> Result<(), OrgraftError> {
    let reference = polish_reference_path(case);
    let mut args = vec![
        "--organelle".to_string(),
        case.organelle.clone(),
        "--subgraph".to_string(),
        case.subgraph.clone(),
        "--draft".to_string(),
        draft.display().to_string(),
        "--reference".to_string(),
        reference.display().to_string(),
        "--reads".to_string(),
        case.reads.display().to_string(),
        "--soft-paths".to_string(),
        config.soft_paths.display().to_string(),
        "--out-dir".to_string(),
        polish_run_out_dir(case, round).display().to_string(),
        "--validate-round".to_string(),
        round.to_string(),
        "--threads".to_string(),
        polish_threads(config).to_string(),
        "--max-rounds".to_string(),
        config.max_rounds.to_string(),
    ];
    if force {
        args.push("--force".to_string());
    }
    append_polish_raw_args(&mut args, &config.polish);
    commands::polish::run(&args)
}

fn homorepair_out_dir(case: &WorkflowCase) -> PathBuf {
    let base = case.postrepair_out_dir.as_ref().unwrap_or(&case.polish_out_dir)
        .join(&case.organelle).join(&case.subgraph);
    base.join("round_hr")
}

fn homorepair_outcome_path(case: &WorkflowCase) -> PathBuf {
    homorepair_out_dir(case).join("logs/run.tsv")
}

// Historical layouts remain visible only in summaries. Current handoff always
// requires run.tsv and the current data paths from the same HR directory.
fn homorepair_existing_outcome_path(case: &WorkflowCase) -> PathBuf {
    let current = homorepair_outcome_path(case);
    if current.exists() { return current; }
    let old_single = homorepair_out_dir(case).join("outcome.tsv");
    if old_single.exists() { return old_single; }
    if homorepair_out_dir(case).exists() { return current; }
    let base = homorepair_out_dir(case).parent().unwrap().to_path_buf();
    for name in ["round_hr2", "round_hr1"] {
        let old = base.join(name).join("outcome.tsv");
        if old.exists() { return old; }
    }
    current
}

fn read_homorepair_run_outcome(path: &Path) -> Result<BTreeMap<String, String>, OrgraftError> {
    let mut lines = BufReader::new(File::open(path)?).lines();
    let invalid = |message: &str| OrgraftError::InvalidArgument(format!(
        "malformed homorepair run table {}: {message}", path.display()));
    if lines.next().transpose()?.as_deref() != Some("section\tname\tmetric\tvalue") {
        return Err(invalid("expected section/name/metric/value header"));
    }
    let mut values = BTreeMap::new();
    for line in lines {
        let line = line?;
        if line.trim().is_empty() { continue; }
        let fields: Vec<_> = line.split('\t').collect();
        if fields.len() != 4 { return Err(invalid("expected four fields")); }
        if fields[0] == "outcome" && fields[1] == "homorepair" {
            if fields[2].is_empty() || values.insert(fields[2].to_string(), fields[3].to_string()).is_some() {
                return Err(invalid("empty or duplicate outcome metric"));
            }
        }
    }
    Ok(values)
}

fn homorepair_probe_records(path: &Path, count: usize) -> Result<BTreeMap<String, BTreeMap<String, String>>, OrgraftError> {
    let invalid = |message: &str| OrgraftError::InvalidArgument(format!(
        "malformed homorepair probe metadata {}: {message}", path.display()));
    let mut records = BTreeMap::<String, BTreeMap<String, String>>::new();
    for row in homorepair_rows(path)? {
        if hr_field(&row, "section")? != "probe" { continue; }
        let name = hr_field(&row, "name")?;
        let metric = hr_field(&row, "metric")?;
        let value = hr_field(&row, "value")?;
        if metric.is_empty() || records.entry(name.to_string()).or_default()
            .insert(metric.to_string(), value.to_string()).is_some() {
            return Err(invalid("empty or duplicate probe metric"));
        }
    }
    if records.len() != count { return Err(invalid("probe records disagree with native_probe_count")); }
    for number in 1..=count {
        let name = format!("probe_{number}");
        let record = records.get(&name).ok_or_else(|| invalid("probe IDs are not consecutive"))?;
        if hr_field(record, "evidence_round")? != (number + 1).to_string()
            || hr_field(record, "evidence_dir")? != format!("03.validate/01.data/{name}")
            || hr_field(record, "probe_fasta")? != format!("03.validate/01.data/{name}.fasta")
            || hr_field(record, "SV_status")? != "pass"
            || hr_field(record, "candidate_count")?.parse::<usize>().ok().filter(|n| *n > 0).is_none() {
            return Err(invalid("probe identity, candidate count or SV status is invalid"));
        }
    }
    Ok(records)
}

const HOMOREPAIR_CANDIDATE_HEADER: &str = "source_rows\ttract_start\ttract_end\tmotif\treference_copies\tselected_copies\treview_status\treason\tevidence_id\tprobe_start\tprobe_end\tprobe_copies\tdepth\ttotal_ALT\tmain_ALT\tremainder\tpool_reads\tuncertain_reads\tcopy_counts\tuncertain_reasons\tcontext_changed\tedit_pos\tedit_ref\tedit_alt";

fn homorepair_candidate_rows(path: &Path) -> Result<Vec<BTreeMap<String, String>>, OrgraftError> {
    if BufReader::new(File::open(path)?).lines().next().transpose()?.as_deref()
        != Some(HOMOREPAIR_CANDIDATE_HEADER) {
        return Err(OrgraftError::InvalidArgument(format!(
            "homorepair candidates table {} does not use the current 24-column schema", path.display())));
    }
    homorepair_rows(path)
}

fn read_homorepair_summary_outcome(path: &Path) -> Result<BTreeMap<String, String>, OrgraftError> {
    if !path.exists() { return Ok(BTreeMap::new()); }
    if path.file_name().is_some_and(|name| name == "run.tsv") {
        return read_homorepair_run_outcome(path);
    }
    let mut values = BTreeMap::new();
    for row in homorepair_rows(path)? {
        let metric = hr_field(&row, "metric")?;
        let value = hr_field(&row, "value")?;
        if values.insert(metric.to_string(), value.to_string()).is_some() {
            return Err(OrgraftError::InvalidArgument(format!(
                "duplicate historical homorepair metric in {}", path.display())));
        }
    }
    Ok(values)
}

fn homorepair_args(config: &WorkflowConfig, case: &WorkflowCase, round: usize) -> Vec<String> {
    let reference = if case.postrepair_out_dir.is_some() {
        // The ordinary round has already selected this component's reference record.
        repair_source_dir(case, round).join("01.inputs/rotated_reference.fasta")
    } else {
        polish_reference_path(case)
    };
    let mut args = vec![
        "--organelle".into(), case.organelle.clone(),
        "--subgraph".into(), case.subgraph.clone(),
        "--draft".into(), repair_source_fasta(case, round).display().to_string(),
        "--reference".into(), reference.display().to_string(),
        "--reads".into(), case.reads.display().to_string(),
        "--soft-paths".into(), config.soft_paths.display().to_string(),
        "--out-dir".into(), case.postrepair_out_dir.as_ref().unwrap_or(&case.polish_out_dir).display().to_string(),
        "--homorepair-frq".into(), config.postrepair.as_ref().map(|p| p.frequency).unwrap_or(config.polish.homorepair_frq).to_string(),
        "--threads".into(), config.postrepair.as_ref().and_then(|p| p.threads).unwrap_or_else(|| polish_threads(config)).to_string(),
        "--max-rounds".into(), config.max_rounds.to_string(),
    ];
    if case.postrepair_out_dir.is_some() {
        args.extend(["--source-dir".into(), repair_source_dir(case, round).display().to_string(),
            "--checked-gfa".into(), case.checked_draft_gfa.display().to_string(),
            "--resolve-details".into(), resolve_details_path(case).display().to_string()]);
    } else {
        args.extend(["--homorepair-source-round".into(), round.to_string()]);
    }
    if let Some(report) = first_polish_report_path(case) {
        args.extend(["--first-polish-report".into(), report.display().to_string()]);
    }
    append_polish_raw_args(&mut args, &config.polish);
    args
}

const POSTREPAIR_REVIEW_HEADER: &str = "component_id\tdecision\tevent_type\tevent\treason";

#[derive(Debug, Clone)]
struct DeclaredPostrepairReview {
    file: PathBuf,
    file_sha256: String,
    decision: Option<ManualKeepDecision>,
    arguments: Vec<String>,
}

fn declared_postrepair_reviews(
    config: &WorkflowConfig,
) -> Result<BTreeMap<(String, String), DeclaredPostrepairReview>, OrgraftError> {
    let mut result = BTreeMap::new();
    let Some(postrepair) = &config.postrepair else {
        return Ok(result);
    };
    for (organelle, path) in &postrepair.review_files {
        if !postrepair.enabled || !path.is_file() {
            return Err(OrgraftError::InvalidArgument(format!(
                "postrepair review requires an enabled stage and readable file: {}",
                path.display()
            )));
        }
        let bytes = fs::read(path)?;
        let text = std::str::from_utf8(&bytes).map_err(|_| {
            OrgraftError::InvalidArgument("postrepair review TSV is not UTF-8".into())
        })?;
        let mut lines = text.lines();
        if lines.next() != Some(POSTREPAIR_REVIEW_HEADER) {
            return Err(OrgraftError::InvalidArgument(format!(
                "postrepair review requires exactly {POSTREPAIR_REVIEW_HEADER}"
            )));
        }
        let mut groups: BTreeMap<String, (String, String, Vec<(String, String)>)> = BTreeMap::new();
        for line in lines.filter(|line| !line.is_empty()) {
            let fields = line.split('\t').collect::<Vec<_>>();
            if fields.len() != 5
                || fields
                    .iter()
                    .any(|field| field.is_empty() || field.contains('\r'))
            {
                return Err(OrgraftError::InvalidArgument(
                    "postrepair review has an empty or malformed row".into(),
                ));
            }
            let cases = config
                .cases
                .iter()
                .filter(|case| {
                    case.enabled && case.organelle == *organelle && case.subgraph == fields[0]
                })
                .collect::<Vec<_>>();
            if cases.len() != 1 || cases[0].postrepair_out_dir.is_none() {
                return Err(OrgraftError::InvalidArgument("postrepair review component is unknown, disabled, ambiguous or lacks a standalone output".into()));
            }
            let group = groups
                .entry(fields[0].into())
                .or_insert_with(|| (fields[1].into(), fields[4].into(), Vec::new()));
            if group.0 != fields[1] || group.1 != fields[4] {
                return Err(OrgraftError::InvalidArgument(
                    "postrepair review component has conflicting decisions or reasons".into(),
                ));
            }
            let event = if fields[2] == "evidence_file" {
                config_relative_path(config.config_path.parent().unwrap(), Path::new(fields[3]))
                    .display()
                    .to_string()
            } else {
                fields[3].into()
            };
            group.2.push((fields[2].into(), event));
        }
        if groups.is_empty() {
            return Err(OrgraftError::InvalidArgument(
                "postrepair review TSV has no declared events".into(),
            ));
        }
        for (component, (mode, reason, events)) in groups {
            let (decision, mut arguments) =
                commands::postrepair::declared_review_arguments(&mode, &events, &reason)?;
            if decision.is_none() {
                fs::File::open(&config.soft_paths)?;
                arguments.extend([
                    "--soft-paths".into(),
                    config.soft_paths.display().to_string(),
                    "--threads".into(),
                    postrepair
                        .threads
                        .unwrap_or_else(|| polish_threads(config))
                        .to_string(),
                ]);
            }
            result.insert(
                (organelle.clone(), component),
                DeclaredPostrepairReview {
                    file: path.clone(),
                    file_sha256: commands::resolve::sha256_hex(&bytes),
                    decision,
                    arguments,
                },
            );
        }
    }
    Ok(result)
}

fn declared_postrepair_review(
    config: &WorkflowConfig,
    case: &WorkflowCase,
) -> Result<Option<DeclaredPostrepairReview>, OrgraftError> {
    Ok(declared_postrepair_reviews(config)?
        .remove(&(case.organelle.clone(), case.subgraph.clone())))
}

fn run_configured_postrepair(args: &[String]) -> Result<(), OrgraftError> {
    let config = load_config_from_args(args)?;
    let case = select_case(&config, args)?;
    let round = if option_value(args, "--round")?.is_some() {
        parse_round(args)?
    } else {
        postrepair_source_round(case)?.unwrap_or(1)
    };
    if declared_postrepair_review(&config, case)?.is_none() {
        return Err(OrgraftError::InvalidArgument(
            "postrepair-run requires an explicit decision for this component".into(),
        ));
    }
    run_homorepair_after_checkpoint(&config, case, round)?;
    Ok(())
}

fn finish_declared_postrepair(
    config: &WorkflowConfig,
    case: &WorkflowCase,
    round: usize,
    declaration: &DeclaredPostrepairReview,
    automatic: Result<(), OrgraftError>,
) -> Result<PathBuf, OrgraftError> {
    let component = homorepair_out_dir(case).parent().unwrap().to_path_buf();
    if let Err(cause) = automatic {
        // Only the producer's explicit terminal review is eligible. I/O, tool
        // failures and partial outputs keep their original error.
        let terminal = format!(
            "postrepair requires review; inspect the selected sequence and evidence in {}",
            component.display()
        );
        if !matches!(&cause, OrgraftError::InvalidArgument(message) if message == &terminal) {
            return Err(cause);
        }
    }
    postrepair_component_for_source(case, round)?;
    if commands::resolve::sha256_hex(&fs::read(&declaration.file)?) != declaration.file_sha256 {
        return Err(OrgraftError::InvalidArgument(
            "postrepair human declaration changed during repair".into(),
        ));
    }
    if fs::symlink_metadata(component.join("manual_review.tsv")).is_ok() {
        commands::postrepair::verify_declared_review(
            &component,
            declaration.decision.as_ref(),
            &declaration.arguments,
        )?;
    } else {
        let mut args = vec![
            "review".into(),
            "--result-dir".into(),
            component.display().to_string(),
        ];
        args.extend(declaration.arguments.clone());
        // Native review checks the full terminal state and evidence. Keep
        // decisions preserve DNA; explicit correction uses native validation.
        commands::postrepair::run(&args)?;
    }
    selected_repair_fasta(config, case, round)
}

fn postrepair_enabled(config: &WorkflowConfig) -> bool {
    config.postrepair.as_ref().map(|p| p.enabled).unwrap_or(config.polish.homorepair)
}

fn run_homorepair_after_checkpoint(config: &WorkflowConfig, case: &WorkflowCase, round: usize) -> Result<PathBuf, OrgraftError> {
    if let Some(declaration) = declared_postrepair_review(config, case)? {
        if postrepair_source_round(case)?.is_some_and(|expected| expected != round) {
            return Err(OrgraftError::InvalidArgument("postrepair run round differs from postrepair_source_dir".into()));
        }
        let component = homorepair_out_dir(case).parent().unwrap().to_path_buf();
        let automatic = if fs::symlink_metadata(component.join("postrepair.tsv")).is_ok()
            || fs::symlink_metadata(component.join("manual_review.tsv")).is_ok() { Ok(()) }
            else { commands::postrepair::run(&homorepair_args(config, case, round)) };
        return finish_declared_postrepair(config, case, round, &declaration, automatic);
    }
    let manual = homorepair_out_dir(case).parent().unwrap().join("manual_review.tsv");
    if case.postrepair_out_dir.is_some() && fs::symlink_metadata(&manual).is_ok() {
        return selected_repair_fasta(config, case, round);
    }
    commands::postrepair::run(&homorepair_args(config, case, round))?;
    selected_repair_fasta(config, case, round)
}

fn postrepair_component_for_source(case: &WorkflowCase, round: usize) -> Result<PathBuf, OrgraftError> {
    if postrepair_source_round(case)?.is_some_and(|expected| expected != round) {
        return Err(OrgraftError::InvalidArgument("postrepair check round differs from postrepair_source_dir".into()));
    }
    let root = case.postrepair_out_dir.as_ref().ok_or_else(|| OrgraftError::InvalidArgument("postrepair output is not configured".into()))?;
    let component = root.join(&case.organelle).join(&case.subgraph);
    let recorded = read_metric_value(&component.join("postrepair.tsv"), "source_round")?
        .ok_or_else(|| OrgraftError::InvalidArgument("postrepair summary lacks its source round".into()))?;
    if fs::canonicalize(recorded)? != fs::canonicalize(repair_source_dir(case, round))? {
        return Err(OrgraftError::InvalidArgument("postrepair source does not match this workflow case and round".into()));
    }
    Ok(component)
}

/// Share the verified postrepair handoff between Rust runs and generated scripts.
/// Arguments are emitted as separate lines for a Bash array, never shell code.
fn postrepair_rebuild_args(
    config: &WorkflowConfig,
    case: &WorkflowCase,
    polished: &Path,
) -> Result<Vec<String>, OrgraftError> {
    if !postrepair_enabled(config) || case.postrepair_out_dir.is_none() {
        return Ok(Vec::new());
    }
    let invalid = |message: &str| OrgraftError::InvalidArgument(format!("postrepair rebuild handoff: {message}"));
    let component = case.postrepair_out_dir.as_ref().unwrap().join(&case.organelle).join(&case.subgraph);
    if let Some(declaration) = declared_postrepair_review(config, case)? {
        commands::postrepair::verify_declared_review(&component, declaration.decision.as_ref(), &declaration.arguments)?;
    }
    let summary = component.join("postrepair.tsv");
    let source = read_metric_value(&summary, "source_round")?
        .ok_or_else(|| invalid("source round is absent"))?;
    let round = Path::new(&source).file_name().and_then(|name| name.to_str())
        .and_then(|name| name.strip_prefix("round_")).and_then(|number| number.parse::<usize>().ok())
        .filter(|number| *number > 0 && *number <= WORKFLOW_HARD_MAX_ROUNDS)
        .ok_or_else(|| invalid("source round is invalid"))?;
    let validated = selected_repair_fasta(config, case, round)?;
    if fs::canonicalize(polished)? != fs::canonicalize(&validated)? {
        return Err(invalid("polished FASTA is not the selected native validation input"));
    }
    if read_metric_value(&summary, "selection_method")?.as_deref() != Some("joint_ir") {
        if let Some(source) = &case.postrepair_source_dir {
            if commands::rebuild::ordinary_source_checkpoint(&case.checked_draft_gfa)?.is_some_and(|checkpoint|checkpoint.join("repeat_pairing.status.tsv").is_file()) {
                return Ok(vec!["--ordinary-source-dir".into(), fs::canonicalize(source)?.display().to_string()]);
            }
        }
        return Ok(Vec::new());
    }
    let source_fasta = read_metric_value(&summary, "source_fasta")?
        .ok_or_else(|| invalid("source FASTA is absent"))?;
    if fs::canonicalize(source_fasta)? != fs::canonicalize(repair_source_fasta(case, round))? {
        return Err(invalid("source FASTA is outside this case and round"));
    }
    let component = fs::canonicalize(&component)?;
    let proposal = component.join("joint_ir/proposal");
    let candidate = fs::canonicalize(proposal.join("candidate.gfa"))?;
    let evidence = fs::canonicalize(proposal.join("joint_ir.tsv"))?;
    let selected = fs::canonicalize(proposal.join("candidate.fasta"))?;
    let validation = read_metric_value(&summary, "final_validation_dir")?
        .ok_or_else(|| invalid("final validation directory is absent"))?;
    if [&candidate, &evidence, &selected, &validated, &fs::canonicalize(validation)?]
        .iter().any(|path| !path.starts_with(&component)) {
        return Err(invalid("joint IR files are outside this case's postrepair output"));
    }
    let digest = |path: &Path| -> Result<String, OrgraftError> {
        Ok(commands::resolve::sha256_hex(&fs::read(path)?))
    };
    let mut evidence_bindings = Vec::new();
    for row in homorepair_rows(&summary)? {
        if let Some(path) = hr_field(&row, "metric")?.strip_prefix("evidence_sha256:") {
            if fs::canonicalize(path)? == evidence {
                evidence_bindings.push(hr_field(&row, "value")?.to_string());
            }
        }
    }
    if evidence_bindings != vec![digest(&evidence)?] {
        return Err(invalid("joint IR evidence is not uniquely bound by the postrepair summary"));
    }
    let mut values = BTreeMap::new();
    for row in homorepair_rows(&evidence)? {
        let key = hr_field(&row, "key")?;
        if values.insert(key.to_string(), hr_field(&row, "value")?.to_string()).is_some() {
            return Err(invalid("joint IR evidence has a duplicate key"));
        }
    }
    for (key, expected) in [("status", "selected_for_ordinary_validation"),
        ("R2_path_links", "verified_zero_overlap"), ("opposite_pair_count", "1")] {
        if values.get(key).map(String::as_str) != Some(expected) {
            return Err(invalid("joint IR proposal was not selected with its original checked path"));
        }
    }
    for (path_key, hash_key, expected) in [
        ("checked_gfa_path", "checked_gfa_sha256", case.checked_draft_gfa.clone()),
        ("resolve_details_path", "resolve_details_sha256", resolve_details_path(case)),
    ] {
        let recorded = values.get(path_key).ok_or_else(|| invalid("joint IR source path is absent"))?;
        if fs::canonicalize(recorded)? != fs::canonicalize(&expected)?
            || values.get(hash_key) != Some(&digest(&expected)?) {
            return Err(invalid("joint IR proposal uses a different checked graph or resolve record"));
        }
    }
    for (key, path) in [("candidate_gfa_sha256", candidate.as_path()),
        ("candidate_fasta_sha256", selected.as_path()), ("candidate_fasta_sha256", validated.as_path())] {
        if values.get(key) != Some(&digest(path)?) {
            return Err(invalid("joint IR candidate differs from its bound validation input"));
        }
    }
    let args = vec!["--joint-ir-trial-gfa".into(), candidate.display().to_string(),
        "--joint-ir-evidence".into(), evidence.display().to_string()];
    if args.iter().any(|arg: &String| arg.contains(['\n', '\r'])) {
        return Err(invalid("argument contains a line break"));
    }
    Ok(args)
}

fn selected_repair_fasta(config: &WorkflowConfig, case: &WorkflowCase, round: usize) -> Result<PathBuf, OrgraftError> {
    if case.postrepair_out_dir.is_none() {
        return homorepair_selected_fasta(config, case, round);
    }
    let component = postrepair_component_for_source(case, round)?;
    if let Some(declaration) = declared_postrepair_review(config, case)? {
        commands::postrepair::verify_declared_review(&component, declaration.decision.as_ref(), &declaration.arguments)?;
    }
    // Postrepair owns replay of both ordinary and separately labelled repair
    // evidence. Its public handoff must reject unresolved or unbound decisions.
    let selected = commands::postrepair::selected_fasta(&component)?;
    if read_metric_value(&component.join("postrepair.tsv"), "selection_method")?.as_deref() == Some("ordinary") {
        let frequency = config.postrepair.as_ref().map(|p| p.frequency).unwrap_or(config.polish.homorepair_frq);
        let replayed = commands::postrepair::verify_ordinary_selection(&component, &repair_source_dir(case, round), frequency)?;
        if read_single_fasta(&replayed)?.1 != read_single_fasta(&selected)?.1 {
            return Err(OrgraftError::InvalidArgument("ordinary postrepair evidence selects a different sequence".into()));
        }
    }
    let validation = commands::postrepair::selected_validation_dir(&component)?;
    let report = validation.parent().ok_or_else(|| OrgraftError::InvalidArgument("postrepair validation has no round directory".into()))?.join("logs/report.tsv");
    let rows = homorepair_rows(&report)?;
    let paths: Vec<_> = rows.iter().filter(|row|
        row.get("section").map(String::as_str) == Some("output")
        && row.get("name").map(String::as_str) == Some("validation_fasta")
        && row.get("metric").map(String::as_str) == Some("path"))
        .filter_map(|row| row.get("value")).collect();
    if paths.len() != 1 {
        return Err(OrgraftError::InvalidArgument("postrepair final validation must identify one FASTA".into()));
    }
    let validated = fs::canonicalize(paths[0])?;
    if read_single_fasta(&validated)?.1 != read_single_fasta(&selected)?.1 {
        return Err(OrgraftError::InvalidArgument("postrepair selected sequence differs from its final validation input".into()));
    }
    // Rebuild discovers validation beside 01.inputs/02.polish. Give it the
    // actually validated input rather than a candidate in a repair report folder.
    Ok(validated)
}

fn replay_circle_end_edit(
    data: &Path,
    outcome: &BTreeMap<String, String>,
    reconstructed: &mut String,
    source_round: &Path,
    original: &str,
    ordinary: &[BTreeMap<String, String>],
    prior_edits: usize,
    min_alt_fraction: f64,
) -> Result<usize, OrgraftError> {
    replay_circle_end_edit_with_keep(data, outcome, reconstructed, source_round, original,
        ordinary, prior_edits, min_alt_fraction, None)
}

fn replay_circle_end_edit_with_keep(
    data: &Path,
    outcome: &BTreeMap<String, String>,
    reconstructed: &mut String,
    source_round: &Path,
    original: &str,
    ordinary: &[BTreeMap<String, String>],
    prior_edits: usize,
    min_alt_fraction: f64,
    kept: Option<&ManualKeepRepeatExtra>,
) -> Result<usize, OrgraftError> {
    let invalid = |message: &str| OrgraftError::InvalidArgument(format!(
        "homorepair circle-end {message}; rebuild stopped"));
    let parse_outcome = |key: &str| -> Result<usize, OrgraftError> {
        outcome.get(key).map(String::as_str).unwrap_or("0").parse::<usize>()
            .map_err(|_| invalid(&format!("{key} is invalid")))
    };
    let candidate_count = parse_outcome("circle_end_candidate_count")?;
    let adopted_count = parse_outcome("circle_end_adopted_count")?;
    let review_count = parse_outcome("circle_end_review_count")?;
    let reason = outcome.get("circle_end_reason").map(String::as_str).unwrap_or("disabled");
    if adopted_count > candidate_count || review_count != usize::from(kept.is_some()) {
        return Err(invalid("candidates are unresolved or internally inconsistent"));
    }
    if reason == "disabled" {
        if candidate_count != 0 || adopted_count != 0 || kept.is_some() {
            return Err(invalid("disabled stage has candidate decisions"));
        }
        return Ok(0);
    }
    let root = data.join("circle_end");
    let candidate_path = root.join("candidates.tsv");
    const HEADER: &str = "source_rows\tsource_length\tcut_0based\trotated_start_0based\trotated_end_0based\toriginal_start_0based\toriginal_end_0based\tmotif\treference_copies\tselected_copies\tdepth\ttotal_ALT\tmain_ALT\tremainder\tuncertain_reads\tstatus\treason";
    if BufReader::new(File::open(&candidate_path)?).lines().next().transpose()?.as_deref() != Some(HEADER) {
        return Err(invalid("candidate table schema is invalid"));
    }
    let rows = homorepair_rows(&candidate_path)?;
    if rows.len() != candidate_count || candidate_count > 1 {
        return Err(invalid("candidate table does not match one boundary decision"));
    }
    if candidate_count == 0 {
        if reason != "no_boundary_candidate" || kept.is_some() { return Err(invalid("empty candidate reason is invalid")); }
        return Ok(0);
    }
    let row = &rows[0];
    let manual = kept.is_some();
    let ordinary_skip = hr_field(row, "status")? == "skipped"
        && hr_field(row, "reason")? == "ordinary_probe_confirms_original_copies";
    if reason != hr_field(row, "reason")?
        || (ordinary_skip && adopted_count != 0)
        || (manual && (adopted_count != 0 || hr_field(row, "status")? != "review_required"
            || reason != "fixed_read_pool_strict_gate_failed" || hr_field(row, "selected_copies")? != "."))
        || (!ordinary_skip && !manual && (adopted_count != 1 || hr_field(row, "status")? != "adopted"
            || reason != "unique_main_ALT_above_fixed_read_pool_remainder")) {
        return Err(invalid("candidate has no strict adoption or ordinary ownership decision"));
    }
    if hr_field(row, "source_rows")?.parse::<usize>().ok().filter(|n| *n > 0).is_none() {
        return Err(invalid("candidate source row is invalid"));
    }
    let number = |key: &str| hr_field(row, key)?.parse::<usize>()
        .map_err(|_| invalid(&format!("candidate {key} is invalid")));
    let source = reconstructed.clone();
    let length = number("source_length")?;
    let cut = number("cut_0based")?;
    let start = number("rotated_start_0based")?;
    let end = number("rotated_end_0based")?;
    let mapped_start = number("original_start_0based")?;
    let mapped_end = number("original_end_0based")?;
    let copies = number("reference_copies")?;
    let selected_copies = if manual { copies } else { number("selected_copies")? };
    let motif = hr_field(row, "motif")?;
    if length != source.len() || cut == 0 || cut >= length || start == 0 || start >= end
        || end > length || motif.is_empty() || !motif.is_ascii()
        || copies == 0 || motif.len().checked_mul(copies) != Some(end - start) {
        return Err(invalid("source length, cut, or repeat tract is invalid"));
    }
    let boundary = length - cut;
    if start > boundary.saturating_add(120) || end < boundary.saturating_sub(120) {
        return Err(invalid("candidate is outside the original-origin window"));
    }
    if start < boundary && end > boundary {
        return Err(invalid("repeat tract crosses the original linear origin"));
    }
    let expected_map = if end <= boundary {
        (cut.checked_add(start).ok_or_else(|| invalid("mapped coordinate overflow"))?,
            cut.checked_add(end).ok_or_else(|| invalid("mapped coordinate overflow"))?)
    } else { (start - boundary, end - boundary) };
    if (mapped_start, mapped_end) != expected_map {
        return Err(invalid("original coordinates disagree with the rotated tract"));
    }
    let rotated = format!("{}{}", &source[cut..], &source[..cut]);
    let tract = rotated.get(start..end).ok_or_else(|| invalid("rotated tract is outside the source"))?;
    if tract.as_bytes().chunks(motif.len()).any(|unit| unit != motif.as_bytes())
        || source.get(mapped_start..mapped_end) != Some(tract) {
        return Err(invalid("repeat tract differs from the source sequence"));
    }
    if let Some(kept) = kept {
        if kept.start != mapped_start + 1 || kept.end != mapped_end || kept.motif != motif || kept.copies != copies {
            return Err(invalid("manual circle event differs from its current sequence"));
        }
        commands::postrepair::verify_current_circle_source(data, source_round, &source, row, min_alt_fraction)?;
    }
    if ordinary_skip {
        if prior_edits != 0 || source != original || selected_copies != copies {
            return Err(invalid("ordinary ownership requires the unchanged complete source"));
        }
        for key in ["depth", "total_ALT", "main_ALT", "remainder", "uncertain_reads"] {
            if hr_field(row, key)? != "." { return Err(invalid("ordinary skip invents circle read support")); }
        }
        let class_path = root.join("read_classes.tsv");
        if fs::read_to_string(&class_path)? != "source_rows\tread_id\tclass\texact_copies\n" {
            return Err(invalid("ordinary skip has independent circle read classifications"));
        }
        let matching: Vec<_> = ordinary.iter().filter(|candidate| {
            hr_field(candidate, "tract_start").ok().and_then(|n| n.parse::<usize>().ok()) == Some(mapped_start + 1)
                && hr_field(candidate, "tract_end").ok().and_then(|n| n.parse::<usize>().ok()) == Some(mapped_end)
                && hr_field(candidate, "motif").ok() == Some(motif)
                && hr_field(candidate, "reference_copies").ok().and_then(|n| n.parse::<usize>().ok()) == Some(copies)
        }).collect();
        if matching.len() != 1 { return Err(invalid("ordinary ownership is missing or ambiguous")); }
        for other in ordinary.iter().filter(|candidate| !std::ptr::eq(*candidate, matching[0])) {
            let start = hr_field(other, "tract_start")?.parse::<usize>()
                .map_err(|_| invalid("ordinary tract start is invalid"))?.checked_sub(1)
                .ok_or_else(|| invalid("ordinary tract start is zero"))?;
            let end = hr_field(other, "tract_end")?.parse::<usize>()
                .map_err(|_| invalid("ordinary tract end is invalid"))?;
            if start.saturating_sub(1) < mapped_end && mapped_start.saturating_sub(1) < end {
                return Err(invalid("ordinary ownership conflicts with another candidate"));
            }
        }
        commands::postrepair::verify_circle_ordinary_owner(
            source_round, data, original, row, matching[0], min_alt_fraction)?;
        return Ok(0);
    }
    let depth = number("depth")?;
    let total = number("total_ALT")?;
    let main = number("main_ALT")?;
    let remainder = number("remainder")?;
    let uncertain = number("uncertain_reads")?;
    let class_path = root.join("read_classes.tsv");
    if BufReader::new(File::open(&class_path)?).lines().next().transpose()?.as_deref()
        != Some("source_rows\tread_id\tclass\texact_copies") {
        return Err(invalid("read class table schema is invalid"));
    }
    let mut seen = std::collections::HashSet::new();
    let mut exact = BTreeMap::<usize, usize>::new();
    let mut uncertain_rows = 0usize;
    for class in homorepair_rows(&class_path)? {
        if hr_field(&class, "source_rows")? != hr_field(row, "source_rows")? {
            return Err(invalid("read class source does not match the candidate"));
        }
        let id = hr_field(&class, "read_id")?;
        if id.is_empty() || !seen.insert(id.to_string()) {
            return Err(invalid("read pool has an empty or duplicate identifier"));
        }
        match hr_field(&class, "class")? {
            "source" | "ALT" => {
                let copy = hr_field(&class, "exact_copies")?.parse::<usize>()
                    .map_err(|_| invalid("exact copy class is invalid"))?;
                if (copy == copies) != (hr_field(&class, "class")? == "source") {
                    return Err(invalid("read source/ALT class disagrees with exact copies"));
                }
                *exact.entry(copy).or_default() += 1;
            }
            "uncertain" if hr_field(&class, "exact_copies")? == "." => uncertain_rows += 1,
            _ => return Err(invalid("read class is invalid")),
        }
    }
    let alt_counts: Vec<_> = exact.iter().filter(|(copy, _)| **copy != copies)
        .map(|(_, count)| *count).collect();
    let maximum = alt_counts.iter().copied().max().unwrap_or(0);
    if depth == 0 || seen.len() != depth || uncertain_rows != uncertain
        || exact.values().sum::<usize>() + uncertain != depth
        || alt_counts.iter().sum::<usize>() != total || total > depth
        || remainder != depth - total || main != maximum {
        return Err(invalid("read classes do not reproduce the recorded fixed pool"));
    }
    if manual {
        if maximum > remainder && alt_counts.iter().filter(|n| **n == maximum).count() == 1 {
            return Err(invalid("manual circle record claims failure for a strict accepted ALT"));
        }
        return Ok(0);
    }
    if exact.get(&selected_copies).copied() != Some(main)
        || alt_counts.iter().filter(|count| **count == main).count() != 1
        || main <= remainder || selected_copies == copies {
        return Err(invalid("read classes do not reproduce the strict fixed-pool decision"));
    }
    let replacement = motif.repeat(selected_copies);
    reconstructed.replace_range(mapped_start..mapped_end, &replacement);
    let mut rotated_changed = rotated;
    rotated_changed.replace_range(start..end, &replacement);
    let changed_boundary = if end <= boundary {
        boundary.checked_add(replacement.len()).and_then(|n| n.checked_sub(end - start))
            .ok_or_else(|| invalid("changed circular boundary overflows"))?
    } else { boundary };
    if *reconstructed != format!("{}{}", &rotated_changed[changed_boundary..],
        &rotated_changed[..changed_boundary]) {
        return Err(invalid("original-coordinate edit does not replay from the rotated model"));
    }
    Ok(1)
}

fn replay_repeat_extra(
    data: &Path,
    outcome: &BTreeMap<String, String>,
    original: &str,
    ordinary: &[BTreeMap<String, String>],
    reconstructed: &mut String,
) -> Result<(usize, usize), OrgraftError> {
    replay_repeat_extra_with_keep(data, outcome, original, ordinary, reconstructed, &[])
}

fn replay_repeat_extra_with_keep(
    data: &Path,
    outcome: &BTreeMap<String, String>,
    original: &str,
    ordinary: &[BTreeMap<String, String>],
    reconstructed: &mut String,
    keep: &[ManualKeepRepeatExtra],
) -> Result<(usize, usize), OrgraftError> {
    replay_repeat_extra_with_deferred(data, outcome, original, ordinary, reconstructed, keep, &[])
}

fn replay_repeat_extra_with_deferred(
    data: &Path, outcome: &BTreeMap<String, String>, original: &str,
    ordinary: &[BTreeMap<String, String>], reconstructed: &mut String,
    keep: &[ManualKeepRepeatExtra], deferred: &[commands::postrepair::DeferredRepeat],
) -> Result<(usize, usize), OrgraftError> {
    let invalid = |message: &str| OrgraftError::InvalidArgument(format!(
        "homorepair supplemental repeat {message}; rebuild stopped"));
    let metric = |key: &str| -> Result<usize, OrgraftError> {
        outcome.get(key).map(String::as_str).unwrap_or("0").parse::<usize>()
            .map_err(|_| invalid(&format!("{key} is invalid")))
    };
    let candidate_count = metric("repeat_extra_candidate_count")?;
    let adopted_count = metric("repeat_extra_adopted_count")?;
    let review_count = metric("repeat_extra_review_count")?;
    let requested = keep.iter().cloned().collect::<std::collections::BTreeSet<_>>();
    let deferred_set = deferred.iter().map(|site| ManualKeepRepeatExtra {
        start:site.start, end:site.end, motif:site.motif.clone(), copies:site.copies,
    }).collect::<BTreeSet<_>>();
    if (!keep.is_empty() && !deferred.is_empty()) || deferred_set.len() != deferred.len()
        || review_count != keep.len() + deferred.len() || requested.len() != keep.len() || adopted_count > candidate_count
        || (!deferred.is_empty() && adopted_count != 0)
        || (!keep.is_empty() && (adopted_count != 0 || reconstructed != original)) {
        return Err(invalid("candidates are unresolved or internally inconsistent"));
    }
    let mut kept = std::collections::BTreeSet::new();
    let root = data.join("repeat_extra");
    let path = root.join("candidates.tsv");
    if !path.exists() {
        if candidate_count == 0 && adopted_count == 0 && keep.is_empty() { return Ok((0,0)); }
        return Err(invalid("candidate table is missing"));
    }
    let mut checked = 0usize;
    let mut applied = 0usize;
    let mut previous_end = 0usize;
    for row in homorepair_rows(&path)? {
        let number = |key| hr_field(&row,key)?.parse::<usize>()
            .map_err(|_| invalid(&format!("candidate {key} is invalid")));
        let start = number("source_start")?.checked_sub(1)
            .ok_or_else(|| invalid("source start is zero"))?;
        let end = number("source_end")?;
        let copies = number("original_copies")?;
        let motif = hr_field(&row,"motif")?;
        if motif.len()!=1 || !matches!(motif.as_bytes()[0],b'A'|b'C'|b'G'|b'T')
            || copies<=10 || end.checked_sub(start)!=Some(copies)
            || original.get(start..end)!=Some(motif.repeat(copies).as_str())
            || start<previous_end || hr_field(&row,"entry")?!="long_homopolymer_residual" {
            return Err(invalid("source tract or row order disagrees with the original sequence"));
        }
        previous_end=end;
        let owned = ordinary.iter().any(|r| {
            hr_field(r,"tract_start").ok()==Some((start+1).to_string().as_str())
                && hr_field(r,"tract_end").ok()==Some(end.to_string().as_str())
        });
        let status=hr_field(&row,"status")?;
        let reason=hr_field(&row,"reason")?;
        if status=="skipped" {
            if !owned || reason!="site_already_owned_by_old_repeat"
                || hr_field(&row,"selected_copies")?!="."
                || ["depth","total_ALT","main_ALT","remainder","untyped"].iter().any(|k| hr_field(&row,k).ok()!=Some("0")) {
                return Err(invalid("ordinary-owned skip is not supported by its original review"));
            }
            continue;
        }
        if owned { return Err(invalid("ordinary-owned tract was processed again")); }
        checked+=1;
        if status == "review_required" {
            let event = ManualKeepRepeatExtra { start: start + 1, end, motif: motif.into(), copies };
            let is_deferred = deferred_set.contains(&event);
            if (!requested.contains(&event) && !is_deferred) || !kept.insert(event)
                || reason != "source_context_not_unique_on_selected"
                || hr_field(&row, "selected_copies")? != "."
                || hr_field(&row, "copy_counts")? != "."
                || ["depth", "total_ALT", "main_ALT", "remainder", "untyped"].iter()
                    .any(|key| hr_field(&row, key).ok() != Some("0")) {
                return Err(invalid("manual decision does not match the complete unresolved context-ambiguity set"));
            }
            let right = end.checked_add(30).ok_or_else(|| invalid("manual site flank overflows"))?;
            let context = if start < 30 || right > original.len() {
                None
            } else {
                Some(original.get(start - 30..right).ok_or_else(|| invalid("manual source context is invalid"))?)
            };
            // The producer also records this unmeasured reason when a linear
            // endpoint lacks a complete flank. Do not truncate or wrap that flank.
            // Deferred rows retain their separate complete-context requirement.
            if is_deferred && context.is_none() {
                return Err(invalid("deferred site lacks its complete source context"));
            }
            if !is_deferred && context.is_some_and(|context|
                original.as_bytes().windows(context.len()).filter(|window| *window == context.as_bytes()).take(2).count() != 2) {
                return Err(invalid("manual site does not have multiple complete source contexts"));
            }
            continue;
        }
        if status=="unchanged" {
            if reason!="no_unique_dominant_alternative_in_complete_reads"
                || hr_field(&row,"selected_copies")?!="." {
                return Err(invalid("unchanged decision records an edit or unrecognized reason"));
            }
            continue;
        }
        if status!="support_verified" || reason!="unique_main_ALT_above_fixed_pool_remainder" {
            return Err(invalid("contains a candidate without a confirmed fixed-pool decision"));
        }
        let selected=number("selected_copies")?;
        let depth=number("depth")?;
        let total=number("total_ALT")?;
        let main=number("main_ALT")?;
        let remainder=number("remainder")?;
        let unknown=number("untyped")?;
        if selected==0 || selected==copies || depth==0 || total>depth || remainder!=depth-total
            || main<=remainder || main>total || unknown>remainder {
            return Err(invalid("original strict support formula is not met"));
        }
        let mut histogram=BTreeMap::<usize,usize>::new();
        for pair in hr_field(&row,"copy_counts")?.split(',') {
            let (copy,count)=pair.split_once(':').ok_or_else(|| invalid("copy histogram is malformed"))?;
            let copy=copy.parse::<usize>().map_err(|_| invalid("copy histogram number is invalid"))?;
            let count=count.parse::<usize>().map_err(|_| invalid("copy histogram count is invalid"))?;
            if copy==0 || count==0 || histogram.insert(copy,count).is_some() {
                return Err(invalid("copy histogram contains a zero or duplicate entry"));
            }
        }
        let mut ids=std::collections::HashSet::new();
        let mut actual=BTreeMap::<usize,usize>::new();
        let mut actual_unknown=0usize;
        for class in homorepair_rows(&root.join(format!("site_{}_read_classes.tsv",start+1)))? {
            let id=hr_field(&class,"read_id")?;
            if id.is_empty() || !ids.insert(id.to_string()) { return Err(invalid("fixed read pool has an empty or repeated ID")); }
            let copy=hr_field(&class,"copy_count")?;
            if copy=="unknown" { actual_unknown+=1; }
            else {
                let copy=copy.parse::<usize>().map_err(|_| invalid("read copy class is invalid"))?;
                if copy==0 { return Err(invalid("read copy class is zero")); }
                *actual.entry(copy).or_default()+=1;
            }
        }
        let actual_total=actual.iter().filter(|(c,_)| **c!=copies).map(|(_,n)|*n).sum::<usize>();
        let actual_main=actual.iter().filter(|(c,_)| **c!=copies).map(|(_,n)|*n).max().unwrap_or(0);
        if ids.len()!=depth || actual_unknown!=unknown || actual!=histogram || actual_total!=total
            || actual_main!=main || actual.get(&selected).copied()!=Some(main)
            || actual.iter().filter(|(c,n)| **c!=copies && **n==main).count()!=1 {
            return Err(invalid("fixed read classes do not reproduce the unique strict winner"));
        }
        let context_start=start.checked_sub(30).ok_or_else(|| invalid("source lacks thirty-base left flank"))?;
        let context_end=end.checked_add(30).ok_or_else(|| invalid("source flank overflow"))?;
        let context=original.get(context_start..context_end).ok_or_else(|| invalid("source lacks thirty-base right flank"))?;
        let matches=reconstructed.match_indices(context).map(|(p,_)|p).collect::<Vec<_>>();
        if matches.len()!=1 { return Err(invalid("source context is not unique on the preceding selected sequence")); }
        let mapped_start=matches[0]+30;
        let mapped_end=mapped_start+copies;
        let source_model=root.join(format!("source_{applied}"));
        let alternate_model=root.join(format!("site_{}",start+1));
        let (source_id,source_sequence)=read_single_fasta(&source_model.join("candidate.fasta"))?;
        let expected_alt=format!("{}{}{}",&reconstructed[..mapped_start],motif.repeat(selected),&reconstructed[mapped_end..]);
        let (alt_id,alt_sequence)=read_single_fasta(&alternate_model.join("candidate.fasta"))?;
        if source_sequence!=*reconstructed || alt_sequence!=expected_alt
            || source_id.split_whitespace().next()!=alt_id.split_whitespace().next() {
            return Err(invalid("source or ALT model differs from the selected background and candidate edit"));
        }
        let target_id=source_id.split_whitespace().next().ok_or_else(||invalid("source model lacks its record ID"))?;
        let mut fixed_ids=std::collections::HashSet::new();
        for model in [&source_model,&alternate_model] {
            for evidence in homorepair_rows(&model.join("round_2/03.validate/01.data/sv_read_evidence.tsv"))? {
                if hr_field(&evidence,"target_id")?!=target_id || hr_field(&evidence,"alignment_role")?=="no_alignment" { continue; }
                let s=hr_field(&evidence,"target_start")?.parse::<usize>().map_err(|_|invalid("native pool start is invalid"))?;
                let e=hr_field(&evidence,"target_end")?.parse::<usize>().map_err(|_|invalid("native pool end is invalid"))?;
                let anchor=mapped_start.checked_sub(1).ok_or_else(||invalid("candidate lacks an anchor"))?;
                if s==0 || s-1>anchor || anchor>=e { continue; }
                if !matches!(hr_field(&evidence,"alignment_role")?,"primary"|"supplementary"|"secondary") {
                    return Err(invalid("native fixed pool has an unrecognized locus-hit role"));
                }
                fixed_ids.insert(hr_field(&evidence,"read_id")?.to_string());
            }
        }
        if fixed_ids!=ids { return Err(invalid("classified IDs differ from the source/ALT native fixed union")); }
        reconstructed.replace_range(mapped_start..mapped_end,&motif.repeat(selected));
        applied+=1;
    }
    let expected_kept = requested.union(&deferred_set).cloned().collect::<BTreeSet<_>>();
    if checked!=candidate_count || applied!=adopted_count || kept != expected_kept {
        return Err(invalid("candidate, adopted, or manual event count disagrees with the evidence table"));
    }
    Ok((candidate_count,adopted_count))
}

fn long_indel_event_shape(source_ref: &str, source_alt: &str) -> bool {
    if !source_ref.is_ascii() || !source_alt.is_ascii() || source_ref.is_empty() || source_alt.is_empty() {
        return false;
    }
    let prefix = source_ref.bytes().zip(source_alt.bytes()).take_while(|(a, b)| a == b).count();
    let reference = &source_ref[prefix..];
    let alternate = &source_alt[prefix..];
    let suffix = reference.bytes().rev().zip(alternate.bytes().rev()).take_while(|(a, b)| a == b).count();
    matches!((reference.len() - suffix, alternate.len() - suffix), (0, 100..) | (100.., 0))
}

// The native 20 bp review can stop before typing any reads. Replay that
// template failure on the same source, provisional ALT and selected ALT DNA.
fn unmeasured_flank20_template_failure_matches(
    reason: &str, source: &str, source_start: usize, source_ref: &str,
    source_alt: &str, selected: &str, complete_alt: &str,
) -> bool {
    let Some(left_start) = source_start.checked_sub(20) else { return false; };
    let Some(source_end) = source_start.checked_add(source_ref.len()) else { return false; };
    let Some(right_end) = source_end.checked_add(20) else { return false; };
    let (Some(left), Some(right)) = (source.get(left_start..source_start), source.get(source_end..right_end)) else {
        return false;
    };
    let unique = |sequence: &str, template: &str| {
        !template.is_empty() && sequence.as_bytes().windows(template.len())
            .filter(|window| *window == template.as_bytes()).take(2).count() == 1
    };
    let source_unique = unique(source, &format!("{left}{source_ref}{right}"));
    match reason {
        "20bp_source_full_template_not_unique" => !source_unique,
        "20bp_highest_ALT_full_template_not_unique" if source_unique => {
            let temporary = format!("{}{}{}", &source[..source_start], source_alt, &source[source_end..]);
            !unique(&temporary, &format!("{left}{source_alt}{right}"))
                || !unique(selected, &format!("{left}{complete_alt}{right}"))
        }
        _ => false,
    }
}

fn short_complete_allele_shape(source_ref: &str, source_alt: &str, adopted: bool) -> bool {
    let Some(inserted_len) = source_alt.len().checked_sub(source_ref.len()) else {
        return false;
    };
    source_ref.is_ascii()
        && source_alt.is_ascii()
        && inserted_len > 0
        && (!adopted || (12..=99).contains(&inserted_len))
        && (0..=source_ref.len()).any(|split| {
            source_alt.starts_with(&source_ref[..split])
                && source_alt.ends_with(&source_ref[split..])
        })
}

/// Replay the existing ordinary repair evidence without trusting the stage summary.
/// The separate postrepair module also uses this before publishing its handoff.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub(crate) struct ManualKeepEvent {
    pub(crate) pos: usize,
    pub(crate) reference: String,
    pub(crate) alternative: String,
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub(crate) struct ManualKeepRepeatExtra {
    pub(crate) start: usize,
    pub(crate) end: usize,
    pub(crate) motif: String,
    pub(crate) copies: usize,
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub(crate) enum ManualCurrentEvent {
    Indel {
        step: usize,
        event: ManualKeepEvent,
    },
    Circle(ManualKeepRepeatExtra),
    Junction {
        left: String,
        right: String,
        base: String,
        current: String,
    },
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) enum ManualKeepDecision {
    Current(Vec<ManualCurrentEvent>),
    Indel(ManualKeepEvent),
    RepeatExtra(Vec<ManualKeepRepeatExtra>),
}

impl ManualKeepDecision {
    pub(crate) fn len(&self) -> usize {
        match self {
            Self::Indel(_) => 1,
            Self::RepeatExtra(events) => events.len(),
            Self::Current(events) => events.len(),
        }
    }
    pub(crate) fn current(&self) -> &[ManualCurrentEvent] {
        match self {
            Self::Current(events) => events,
            _ => &[],
        }
    }
    pub(crate) fn ordinary_len(&self) -> usize {
        match self {
            Self::Current(events) => events
                .iter()
                .filter(|event| !matches!(event, ManualCurrentEvent::Junction { .. }))
                .count(),
            _ => self.len(),
        }
    }
    fn indel(&self) -> Option<&ManualKeepEvent> {
        match self {
            Self::Indel(event) => Some(event),
            _ => None,
        }
    }
    fn repeat_extra(&self) -> &[ManualKeepRepeatExtra] {
        match self {
            Self::RepeatExtra(events) => events,
            _ => &[],
        }
    }
}

pub(crate) fn verify_ordinary_postrepair(component: &Path, source: &Path, min_alt_fraction: f64) -> Result<PathBuf, OrgraftError> {
    verify_ordinary_postrepair_with_keep(component, source, min_alt_fraction, None)
}

pub(crate) fn verify_ordinary_postrepair_with_keep(
    component: &Path,
    source: &Path,
    min_alt_fraction: f64,
    keep: Option<&ManualKeepDecision>,
) -> Result<PathBuf, OrgraftError> {
    verify_ordinary_postrepair_pass_with_keep(component, source, min_alt_fraction, keep, 1, false)
}

pub(crate) fn verify_ordinary_postrepair_pass(
    component: &Path,
    source: &Path,
    min_alt_fraction: f64,
    pass: usize,
    defer_contexts: bool,
) -> Result<PathBuf, OrgraftError> {
    verify_ordinary_postrepair_pass_with_keep(
        component,
        source,
        min_alt_fraction,
        None,
        pass,
        defer_contexts,
    )
}

fn verify_ordinary_postrepair_pass_with_keep(
    component: &Path,
    source: &Path,
    min_alt_fraction: f64,
    keep: Option<&ManualKeepDecision>,
    pass: usize,
    defer_contexts: bool,
) -> Result<PathBuf, OrgraftError> {
    let invalid =
        |message: &str| OrgraftError::InvalidArgument(format!("ordinary postrepair {message}"));
    if pass == 0 || (keep.is_some() && (pass != 1 || defer_contexts)) {
        return Err(invalid("invalid continuation or mixed manual decision"));
    }
    if !min_alt_fraction.is_finite() || !(0.0..=1.0).contains(&min_alt_fraction) {
        return Err(invalid("ALT fraction must be between 0 and 1"));
    }
    let source = fs::canonicalize(source)?;
    let component = fs::canonicalize(component)?;
    let round = source
        .file_name()
        .and_then(|name| name.to_str())
        .and_then(|name| name.strip_prefix("round_"))
        .and_then(|value| value.parse::<usize>().ok())
        .filter(|round| *round > 0)
        .ok_or_else(|| invalid("source must be an ordinary round_N directory"))?;
    let source_component = source
        .parent()
        .ok_or_else(|| invalid("source has no component directory"))?;
    let source_organelle = source_component
        .parent()
        .ok_or_else(|| invalid("source has no organelle directory"))?;
    let polish_root = source_organelle
        .parent()
        .ok_or_else(|| invalid("source has no polish root"))?;
    let organelle = source_organelle
        .file_name()
        .and_then(|name| name.to_str())
        .ok_or_else(|| invalid("organelle name is invalid"))?;
    let subgraph = source_component
        .file_name()
        .and_then(|name| name.to_str())
        .ok_or_else(|| invalid("component name is invalid"))?;
    let output_organelle = component
        .parent()
        .ok_or_else(|| invalid("output has no organelle directory"))?;
    let output_root = output_organelle
        .parent()
        .ok_or_else(|| invalid("output has no postrepair root"))?;
    if component.file_name() != source_component.file_name()
        || output_organelle.file_name() != source_organelle.file_name()
        || source.file_name().and_then(|name| name.to_str())
            != Some(format!("round_{round}").as_str())
    {
        return Err(invalid("output identity differs from its source round"));
    }
    let section = BTreeMap::from([
        ("organelle".to_string(), organelle.to_string()),
        ("subgraph".to_string(), subgraph.to_string()),
    ]);
    // This read-only check uses only the case identity and the two stage roots.
    let postrepair_root = output_root.to_path_buf();
    let case = WorkflowCase::from_section(
        &section,
        "postrepair",
        "postrepair",
        &output_root.display().to_string(),
        &polish_root.display().to_string(),
        polish_root,
        polish_root,
        output_root,
        Some(&postrepair_root),
    )?;
    let deferred = if defer_contexts {
        commands::postrepair::ordinary_continuation_sites(
            &component,
            &source,
            pass,
            min_alt_fraction,
        )?
    } else {
        Vec::new()
    };
    verify_ordinary_homorepair_at(
        &case,
        round,
        min_alt_fraction,
        keep,
        &commands::postrepair::ordinary_pass_dir(&component, pass),
        &deferred,
    )
}

fn homorepair_selected_fasta(config: &WorkflowConfig, case: &WorkflowCase, round: usize) -> Result<PathBuf, OrgraftError> {
    verify_ordinary_homorepair(case, round, config.polish.homorepair_frq)
}

fn verify_ordinary_homorepair(case: &WorkflowCase, round: usize, min_alt_fraction: f64) -> Result<PathBuf, OrgraftError> {
    verify_ordinary_homorepair_inner(case, round, min_alt_fraction, None)
}

fn verify_ordinary_homorepair_inner(case: &WorkflowCase, round: usize, min_alt_fraction: f64, keep: Option<&ManualKeepDecision>) -> Result<PathBuf, OrgraftError> {
    verify_ordinary_homorepair_at(case, round, min_alt_fraction, keep, &homorepair_out_dir(case), &[])
}

fn verify_ordinary_homorepair_at(
    case: &WorkflowCase, round: usize, min_alt_fraction: f64, keep: Option<&ManualKeepDecision>,
    hr: &Path, deferred: &[commands::postrepair::DeferredRepeat],
) -> Result<PathBuf, OrgraftError> {
    let run = hr.join("logs/run.tsv");
    let outcome = read_homorepair_run_outcome(&run)?;
    let invalid = |message: &str| OrgraftError::InvalidArgument(format!("homorepair {message}; rebuild stopped"));
    let value = |key: &str| -> Result<String, OrgraftError> {
        outcome.get(key).cloned().ok_or_else(|| invalid(&format!("outcome is missing {key}")))
    };
    let kept_ordinary = keep.map_or(0, ManualKeepDecision::ordinary_len);
    let expected_status = if kept_ordinary > 0 || !deferred.is_empty() { "review_required" } else { "complete" };
    let expected_unresolved = (kept_ordinary + deferred.len()).to_string();
    if value("status")? != expected_status || value("unresolved_count")? != expected_unresolved
        || value("stage_count")? != "1" || value("source_round")? != round.to_string() {
        return Err(invalid("is incomplete or does not match the current source round"));
    }
    let probes = value("native_probe_count")?.parse::<usize>().map_err(|_| invalid("probe count is invalid"))?;
    if probes > 2 { return Err(invalid("probe count exceeds the two internal probes")); }
    let probe_records = homorepair_probe_records(&run, probes)?;
    let local_batches = value("local_batch_count")?.parse::<usize>().map_err(|_| invalid("local batch count is invalid"))?;
    let source = fs::canonicalize(polished_aln_path(case, round))?;
    let original = read_single_fasta(&source)?.1;
    let snapshot = read_single_fasta(&hr.join("01.inputs/linear_subgraph.round_hr.fasta"))?.1;
    if original.is_empty() || original != snapshot.to_ascii_uppercase()
        || fs::canonicalize(Path::new(&value("input_fasta")?))? != source {
        return Err(invalid("input sequence does not match its source snapshot"));
    }
    let data = hr.join("03.validate/01.data");
    let reports = hr.join("03.validate/03.reports");
    let selected = PathBuf::from(value("selected_fasta")?);
    let selected_sequence = read_single_fasta(&selected)?.1;
    let reviews = homorepair_candidate_rows(&data.join("candidates.tsv"))?;
    let adopted = value("adopted_edit_count")?.parse::<usize>().map_err(|_| invalid("edit count is invalid"))?;
    let expected_selected = if adopted == 0 { source.clone() } else { reports.join("polished.fasta") };
    if fs::canonicalize(&selected)? != fs::canonicalize(&expected_selected)? {
        return Err(invalid("selected FASTA is not the recorded source or revised sequence"));
    }
    let indel_keep = keep.and_then(ManualKeepDecision::indel);
    let extra_keep = keep.map(ManualKeepDecision::repeat_extra).unwrap_or(&[]);
    let manual_reason = if indel_keep.is_some() { "nonrepeat_indel_candidates_unconfirmed" }
        else { "supplemental_repeat_candidates_unconfirmed" };
    if keep.is_some_and(|decision| !matches!(decision, ManualKeepDecision::Current(_))) && (keep.is_some_and(|decision| decision.len() == 0)
        || adopted != 0 || selected_sequence != original
        || value("reason")? != manual_reason
        || ["pending_count", "uncertain_count", "indel_nested_repeat_edits"].iter()
            .any(|key| outcome.get(*key).map(String::as_str).unwrap_or("0") != "0")) {
        return Err(invalid("manual keep-source requires an unchanged source and only its explicitly listed unresolved events"));
    }
    let current = keep.map(ManualKeepDecision::current).unwrap_or(&[]);
    if !current.is_empty() && (["pending_count", "uncertain_count"].iter()
        .any(|key| outcome.get(*key).map(String::as_str).unwrap_or("0") != "0")
        || !matches!(value("reason")?.as_str(), "nonrepeat_indel_candidates_unconfirmed" | "circle_end_candidates_unconfirmed") && kept_ordinary > 0) {
        return Err(invalid("keep-current has unrelated ordinary reviews"));
    }
    let current_indels = current.iter().filter_map(|event| match event {
        ManualCurrentEvent::Indel { step, event } => Some((*step, event)), _ => None,
    }).collect::<Vec<_>>();
    let current_circle = current.iter().filter_map(|event| match event {
        ManualCurrentEvent::Circle(site) => Some(site), _ => None,
    }).collect::<Vec<_>>();
    if current_circle.len() > 1 { return Err(invalid("keep-current has multiple circle decisions")); }
    let mut edits = Vec::new();
    let mut source_rows = std::collections::HashSet::new();
    for row in &reviews {
        let status = hr_field(row, "review_status")?;
        if !matches!(status, "support" | "suggest_modify" | "support_verified") {
            return Err(invalid("contains an unconfirmed candidate or screened-out row"));
        }
        let number = |key| hr_field(row, key)?.parse::<usize>().map_err(|_| invalid(&format!("candidate {key} is invalid")));
        for identifier in hr_field(row, "source_rows")?.split(',') {
            let identifier = identifier.parse::<usize>().map_err(|_| invalid("candidate source row is invalid"))?;
            if identifier == 0 || !source_rows.insert(identifier) {
                return Err(invalid("candidate source rows are zero or duplicated"));
            }
        }
        if probes == 0 { return Err(invalid("candidate decisions have no internal probe evidence")); }
        let evidence_id = hr_field(row, "evidence_id")?;
        let reference_copies = number("reference_copies")?;
        let selected_copies = number("selected_copies")?;
        let motif = hr_field(row, "motif")?;
        let start = number("tract_start")?.checked_sub(1).ok_or_else(|| invalid("tract start is zero"))?;
        let end = number("tract_end")?;
        let tract = original.get(start..end).ok_or_else(|| invalid("tract is outside the original sequence"))?;
        if motif.is_empty() || !motif.is_ascii()
            || reference_copies == 0 || motif.len().checked_mul(reference_copies) != Some(tract.len())
            || tract.as_bytes().chunks(motif.len()).any(|unit| unit != motif.as_bytes()) {
            return Err(invalid("tract motif or reference copy number disagrees with the source"));
        }
        let probe_start = number("probe_start")?.checked_sub(1).ok_or_else(|| invalid("probe start is zero"))?;
        let probe_copies = number("probe_copies")?;
        let probe_length = motif.len().checked_mul(probe_copies).ok_or_else(|| invalid("probe length overflow"))?;
        if number("probe_end")?.checked_sub(probe_start) != Some(probe_length) {
            return Err(invalid("probe coordinates and copy number disagree"));
        }
        let depth = number("depth")?;
        let total = number("total_ALT")?;
        let main = number("main_ALT")?;
        let remainder = number("remainder")?;
        if status == "support_verified" {
            let batch = evidence_id.strip_prefix("batch_").and_then(|n| n.parse::<usize>().ok())
                .filter(|n| *n > 0 && *n <= local_batches)
                .ok_or_else(|| invalid("local candidate has an invalid batch ID"))?;
            if evidence_id != format!("batch_{batch:03}") {
                return Err(invalid("local candidate batch ID is not canonical"));
            }
            // Local support chooses the unique highest exact-copy count, including
            // the probe's own copy number. Unclear reads never become REF votes.
            let mut counts = BTreeMap::new();
            for entry in hr_field(row, "copy_counts")?.split(',') {
                let (copy, count) = entry.split_once(':').ok_or_else(|| invalid("local copy histogram is invalid"))?;
                let copy = copy.parse::<usize>().map_err(|_| invalid("local copy number is invalid"))?;
                let count = count.parse::<usize>().map_err(|_| invalid("local copy count is invalid"))?;
                if counts.insert(copy, count).is_some() { return Err(invalid("local copy histogram contains duplicate copies")); }
            }
            let exact_depth = counts.values().try_fold(0usize, |sum, n| sum.checked_add(*n))
                .ok_or_else(|| invalid("local depth overflow"))?;
            let maximum = counts.values().copied().max().unwrap_or(0);
            let selected_count = counts.get(&selected_copies).copied().unwrap_or(0);
            let probe_count = counts.get(&probe_copies).copied().unwrap_or(0);
            let main_count = counts.iter().filter(|(copy, _)| **copy != probe_copies)
                .map(|(_, count)| *count).max().unwrap_or(0);
            let mut reasons = std::collections::HashSet::new();
            let mut uncertain_count = 0usize;
            for entry in hr_field(row, "uncertain_reasons")?.split(';').filter(|entry| !entry.is_empty()) {
                let (reason, count) = entry.rsplit_once(':').ok_or_else(|| invalid("local uncertain reasons are invalid"))?;
                let count = count.parse::<usize>().map_err(|_| invalid("local uncertain count is invalid"))?;
                if reason.is_empty() || !reasons.insert(reason) { return Err(invalid("local uncertain reasons are empty or duplicated")); }
                uncertain_count = uncertain_count.checked_add(count).ok_or_else(|| invalid("local uncertain count overflow"))?;
            }
            if depth == 0 || depth != exact_depth || remainder != probe_count
                || total != depth - probe_count || main != main_count
                || depth.checked_add(number("uncertain_reads")?) != Some(number("pool_reads")?)
                || uncertain_count != number("uncertain_reads")?
                || !matches!(hr_field(row, "context_changed")?, "0" | "1")
                || selected_count != maximum || counts.values().filter(|n| **n == maximum).count() != 1 {
                return Err(invalid("local support does not identify one supported copy in its retained read pool"));
            }
        } else {
            if depth == 0 || total > depth || remainder != depth-total || main <= remainder || main > total
                || !probe_records.contains_key(evidence_id) {
                return Err(invalid("review does not satisfy the internal ALT support bound"));
            }
            for key in ["pool_reads", "uncertain_reads", "copy_counts", "uncertain_reasons", "context_changed"] {
                if hr_field(row, key)? != "." { return Err(invalid("native candidate contains local-only evidence")); }
            }
            if (selected_copies == reference_copies) != (status == "support") {
                return Err(invalid("review copy number and decision disagree"));
            }
        }
        let edit_pos = hr_field(row, "edit_pos")?;
        let edit_ref = hr_field(row, "edit_ref")?;
        let edit_alt = hr_field(row, "edit_alt")?;
        if selected_copies == reference_copies {
            if [edit_pos, edit_ref, edit_alt] != [".", ".", "."] {
                return Err(invalid("unchanged candidate has an applied edit"));
            }
        } else {
            let anchor = start.checked_sub(1).ok_or_else(|| invalid("applied edit lacks a source anchor"))?;
            let replacement_len = motif.len().checked_mul(selected_copies)
                .filter(|n| *n <= selected_sequence.len()).ok_or_else(|| invalid("selected repeat length is invalid"))?;
            let replacement = motif.repeat(selected_copies);
            if number("edit_pos")? != start || original.get(anchor..end) != Some(edit_ref)
                || edit_alt.len().checked_sub(1) != Some(replacement_len)
                || edit_alt.get(..1) != original.get(anchor..start)
                || edit_alt.get(1..) != Some(replacement.as_str()) {
                return Err(invalid("applied edit does not match its original anchored tract and selected copies"));
            }
            edits.push((anchor, end, edit_alt.to_string()));
        }
    }
    // Replay anchored alleles independently against the original sequence.
    // The source candidate screen also includes anchors in overlap checks.
    edits.sort_unstable_by_key(|edit| edit.0);
    let mut reconstructed = String::new();
    let mut cursor = 0usize;
    for (start, end, replacement) in &edits {
        if *start < cursor { return Err(invalid("applied anchored repeat intervals overlap")); }
        reconstructed.push_str(&original[cursor..*start]);
        reconstructed.push_str(replacement);
        cursor = *end;
    }
    reconstructed.push_str(&original[cursor..]);
    let (extra_candidates,extra_adopted) = if !deferred.is_empty() {
        replay_repeat_extra_with_deferred(&data,&outcome,&original,&reviews,&mut reconstructed,&[],deferred)?
    } else if extra_keep.is_empty() {
        replay_repeat_extra(&data,&outcome,&original,&reviews,&mut reconstructed)?
    } else {
        replay_repeat_extra_with_keep(&data,&outcome,&original,&reviews,&mut reconstructed,extra_keep)?
    };
    let indel_candidates = outcome.get("indel_candidate_count").map(String::as_str).unwrap_or("0")
        .parse::<usize>().map_err(|_| invalid("nonrepeat InDel candidate count is invalid"))?;
    let indel_adopted = outcome.get("indel_adopted_count").map(String::as_str).unwrap_or("0")
        .parse::<usize>().map_err(|_| invalid("nonrepeat InDel adopted count is invalid"))?;
    let indel_reviews = outcome.get("indel_review_count").map(String::as_str).unwrap_or("0")
        .parse::<usize>().map_err(|_| invalid("nonrepeat InDel review count is invalid"))?;
    if (indel_keep.is_none() && indel_reviews != current_indels.len()) || indel_adopted > indel_candidates || indel_reviews > indel_candidates
        || (indel_keep.is_some() && (indel_reviews != 1 || indel_candidates != 1 || indel_adopted != 0)) {
        return Err(invalid("nonrepeat InDel candidates are unresolved or internally inconsistent"));
    }
    let mut kept_event_verified = !extra_keep.is_empty() || !current.is_empty();
    let indel_base = reconstructed.clone();
    let mut indel_frames = BTreeMap::from([(0usize, indel_base)]);
    if indel_candidates > 0 {
        let rows = homorepair_rows(&data.join("indel_candidates.tsv"))?;
        if rows.len() != indel_candidates {
            return Err(invalid("nonrepeat InDel candidate table has the wrong row count"));
        }
        let flank20_path = data.join("indel_flank20_candidates.tsv");
        let flank20_rows = if flank20_path.exists() { homorepair_rows(&flank20_path)? } else { Vec::new() };
        let complete_path = data.join("indel_short_candidates.tsv");
        let complete_rows = if complete_path.exists() { homorepair_rows(&complete_path)? } else { Vec::new() };
        let event_path = data.join("indel_event_candidates.tsv");
        let event_rows = if event_path.exists() { homorepair_rows(&event_path)? } else { Vec::new() };
        let row_key = |row: &BTreeMap<String, String>| -> Result<String, OrgraftError> {
            Ok(["step", "source_row", "source_pos", "source_ref", "source_alt"]
                .iter().map(|field| hr_field(row, field)).collect::<Result<Vec<_>, _>>()?.join("\t"))
        };
        let current_key = |row: &BTreeMap<String, String>| -> bool {
            current_indels.iter().any(|(step, event)| hr_field(row, "step").ok() == Some(step.to_string().as_str())
                && hr_field(row, "source_pos").ok() == Some(event.pos.to_string().as_str())
                && hr_field(row, "source_ref").ok() == Some(event.reference.as_str())
                && hr_field(row, "source_alt").ok() == Some(event.alternative.as_str()))
        };
        let mut kept_rows = BTreeMap::<usize, Vec<&BTreeMap<String, String>>>::new();
        let mut flank20_by_key = BTreeMap::new();
        for row in &flank20_rows {
            let status = hr_field(row, "status")?;
            let reason = hr_field(row, "reason")?;
            let accepted_20 = matches!(status, "adopted" | "not_selected")
                && matches!(reason, "unique_main_ALT_above_fixed_20bp_pool_remainder" | "selected_after_fixed_20bp_review");
            let failed_20 = status == "review_required" && matches!(reason,
                "20bp_fixed_pool_gate_failed" | "20bp_source_full_template_not_unique"
                    | "20bp_highest_ALT_full_template_not_unique");
            if !(accepted_20 || failed_20)
                || flank20_by_key.insert(row_key(row)?, row).is_some() {
                return Err(invalid("20 bp InDel review has an unresolved or duplicate decision"));
            }
        }
        let mut complete_by_key = BTreeMap::new();
        for row in &complete_rows {
            let status = hr_field(row, "status")?;
            let reason = hr_field(row, "reason")?;
            let source_ref = hr_field(row, "source_ref")?;
            let source_alt = hr_field(row, "source_alt")?;
            if !((status == "adopted" && reason == "unique_main_ALT_above_fixed_pool_remainder")
                || (status == "not_selected" && reason == "selected_after_short_complete_allele")
                || (current_key(row) && status == "review_required" && reason == "short_complete_allele_gate_failed"))
                || hr_field(row, "nested_repeat_edits")? != "0"
                || !short_complete_allele_shape(source_ref, source_alt, status == "adopted")
                || complete_by_key.insert(row_key(row)?, row).is_some() {
                return Err(invalid("short complete-allele review is unresolved or duplicated"));
            }
        }
        let junction_path = data.join("indel_junction_candidates.tsv");
        let junction_rows = if junction_path.exists() { homorepair_rows(&junction_path)? } else { Vec::new() };
        let mut junction_by_key = BTreeMap::new();
        for row in &junction_rows {
            if hr_field(row, "status")? != "adopted"
                || hr_field(row, "reason")? != "unique_main_junction_above_fixed_pool_remainder"
                || hr_field(row, "nested_repeat_edits")? != "0"
                || junction_by_key.insert(row_key(row)?, row).is_some() {
                return Err(invalid("exact-junction decision is not one unambiguous unmodified-model adoption"));
            }
        }
        let mut event_by_key = BTreeMap::new();
        for row in &event_rows {
            let status = hr_field(row, "status")?;
            let reason = hr_field(row, "reason")?;
            if !((status == "adopted" && reason == "unique_main_event_above_fixed_pool_remainder")
                || (status == "not_selected" && reason == "selected_after_paired_CIGAR_event")
                || (status == "review_required" && reason == "paired_event_strict_gate_failed"
                    && (junction_by_key.contains_key(&row_key(row)?) || current_key(row))))
                || (status == "not_selected" && ["edit_pos", "edit_ref", "edit_alt"].iter().any(|field| hr_field(row, field).ok() != Some(".")))
                || event_by_key.insert(row_key(row)?, row).is_some() {
                return Err(invalid("long InDel event review is unresolved, duplicated, or has an invalid event shape"));
            }
        }
        #[derive(Clone, Copy, PartialEq, Eq)]
        enum IndelFrame { Flank60, Flank20, CompleteShort, LongEvent, ExactJunction }
        let mut selected_rows = Vec::new();
        let mut unmeasured_flank20 = BTreeMap::new();
        for row in &rows {
            let key = row_key(row)?;
            let matching = flank20_by_key.remove(&key);
            let complete = complete_by_key.remove(&key);
            let event = event_by_key.remove(&key);
            let junction = junction_by_key.remove(&key);
            if current_key(row) {
                if junction.is_some() || hr_field(row, "status")? != "review_required"
                    || hr_field(row, "reason")? != "highest_source_ALT_failed_fixed_pool_gate"
                    || matching.is_none() || [matching, complete, event].into_iter().flatten()
                        .any(|r| hr_field(r, "status").ok() != Some("review_required")) {
                    return Err(invalid("keep-current does not identify an unadopted residual candidate"));
                }
                let step = hr_field(row, "step")?.parse::<usize>().map_err(|_| invalid("invalid kept step"))?;
                kept_rows.entry(step).or_default().push(row);
                continue;
            }
            let event = if let Some(junction) = junction {
                let failed = event.ok_or_else(|| invalid("junction adoption lacks its original failed CIGAR gate"))?;
                if hr_field(failed, "status")? != "review_required"
                    || hr_field(failed, "reason")? != "paired_event_strict_gate_failed"
                    || ["edit_pos", "edit_ref", "edit_alt"].iter().any(|f| hr_field(failed, f).ok() != Some(".")) {
                    return Err(invalid("junction adoption changed the original CIGAR decision"));
                }
                for field in ["source_rows", "source_support", "source_depth", "source_total_ALT", "trial_candidates",
                              "tested_probes", "probe_pos", "probe_delta", "entry_event_support", "nested_repeat_edits"] {
                    if hr_field(failed, field)? != hr_field(junction, field)? {
                        return Err(invalid("junction adoption changed its original source or model"));
                    }
                }
                Some(junction)
            } else { event };
            if complete.is_some() && event.is_some() {
                return Err(invalid("one InDel candidate has both short and long fallback decisions"));
            }
            match hr_field(row, "status")? {
                "adopted" | "not_selected" if matching.is_none() && complete.is_none() && event.is_none() => {
                    if hr_field(row, "status")? == "adopted" {
                        selected_rows.push((row, IndelFrame::Flank60));
                    }
                }
                "review_required" if hr_field(row, "reason")? == "highest_source_ALT_failed_fixed_pool_gate" => {
                    let flank20 = matching.ok_or_else(|| invalid("60 bp InDel failure lacks its 20 bp decision"))?;
                    for field in ["source_rows", "source_support", "source_depth", "source_total_ALT", "trial_candidates", "tested_probes", "probe_pos", "probe_delta"] {
                        if hr_field(row, field)? != hr_field(flank20, field)? {
                            return Err(invalid("20 bp InDel decision does not match the 60 bp source and fixed read pool"));
                        }
                    }
                    if hr_field(flank20, "depth")? == "." {
                        let allowed = (keep.is_none() || keep.is_some_and(|k| matches!(k, ManualKeepDecision::Current(_)))) && complete.is_none()
                            && event.is_some_and(|event| hr_field(event, "status").ok() == Some("adopted"))
                            && hr_field(flank20, "status")? == "review_required"
                            && matches!(hr_field(flank20, "reason")?,
                                "20bp_source_full_template_not_unique" | "20bp_highest_ALT_full_template_not_unique")
                            && ["depth", "total_ALT", "main_ALT", "remainder", "uncertain_reads", "untried_exact_reads",
                                "edit_pos", "edit_ref", "edit_alt"].iter().all(|field| hr_field(flank20, field).ok() == Some("."))
                            && hr_field(flank20, "nested_repeat_edits")? == "0";
                        let depth60 = hr_field(row, "depth")?.parse::<usize>()
                            .map_err(|_| invalid("unmeasured 20 bp review lacks a numeric 60 bp fixed pool"))?;
                        if !allowed || depth60 == 0 {
                            return Err(invalid("unmeasured 20 bp review is not an unedited template failure followed by one long-event adoption"));
                        }
                        unmeasured_flank20.insert(key.clone(), (row, flank20));
                    } else if hr_field(row, "depth")? != hr_field(flank20, "depth")? {
                        return Err(invalid("20 bp InDel decision does not match the 60 bp source and fixed read pool"));
                    }
                    match hr_field(flank20, "status")? {
                        "adopted" | "not_selected" if complete.is_none() && event.is_none() => {
                            if hr_field(flank20, "status")? == "adopted" {
                                selected_rows.push((flank20, IndelFrame::Flank20));
                            }
                        }
                        "review_required" => {
                            if let Some(kept) = indel_keep.filter(|_| complete.is_none() && event.is_none()) {
                                verify_manual_kept_indel(&data, row, flank20, kept, &original, min_alt_fraction)?;
                                kept_event_verified = true;
                                continue;
                            }
                            let (fallback, frame) = match (complete, event) {
                                (Some(complete), None) => (complete, IndelFrame::CompleteShort),
                                (None, Some(event)) => (event, if junction.is_some() { IndelFrame::ExactJunction } else { IndelFrame::LongEvent }),
                                _ => return Err(invalid("failed 20 bp InDel review lacks a complete-allele or long-event decision")),
                            };
                            for field in ["source_rows", "source_support", "source_depth", "source_total_ALT", "trial_candidates", "tested_probes", "probe_pos", "probe_delta"] {
                                if hr_field(flank20, field)? != hr_field(fallback, field)? {
                                    return Err(invalid("InDel fallback decision changed its source or tested probes"));
                                }
                            }
                            if matches!(frame, IndelFrame::LongEvent | IndelFrame::ExactJunction) {
                                for source in [row, flank20] {
                                    if hr_field(source, "entry_event_support")? != hr_field(fallback, "entry_event_support")? {
                                        return Err(invalid("long InDel event decision changed its source entry support"));
                                    }
                                }
                            } else if hr_field(flank20, "depth")? != "."
                                && hr_field(flank20, "depth")? != hr_field(fallback, "depth")? {
                                return Err(invalid("short complete-allele decision changed the 20 bp fixed read pool"));
                            }
                            for failed in [row, flank20] {
                                if ["edit_pos", "edit_ref", "edit_alt"].iter()
                                    .any(|field| hr_field(failed, field).ok() != Some(".")) {
                                    return Err(invalid("failed InDel gate contains an applied edit"));
                                }
                            }
                            if hr_field(fallback, "status")? == "adopted" {
                                selected_rows.push((fallback, frame));
                            }
                        }
                        _ => return Err(invalid("20 bp InDel decision disagrees with the short complete-allele review")),
                    }
                }
                _ => return Err(invalid("nonrepeat InDel candidate table contains an unresolved decision")),
            }
        }
        if !flank20_by_key.is_empty() || !complete_by_key.is_empty() || !event_by_key.is_empty() || !junction_by_key.is_empty() {
            return Err(invalid("InDel fallback review has no matching 60 bp source row"));
        }
        if kept_rows.values().map(Vec::len).sum::<usize>() != current_indels.len() {
            return Err(invalid("keep-current event set is missing, duplicated, or contains an adopted event"));
        }
        if selected_rows.len() != indel_adopted {
            return Err(invalid("nonrepeat InDel adopted status disagrees with the run outcome"));
        }
        selected_rows.sort_by_key(|(row, _)| hr_field(row, "step").ok().and_then(|step| step.parse::<usize>().ok()).unwrap_or(0));
        let classes60 = if selected_rows.iter().any(|(_, frame)| matches!(frame, IndelFrame::Flank60 | IndelFrame::CompleteShort | IndelFrame::LongEvent | IndelFrame::ExactJunction)) { homorepair_rows(&data.join("indel/read_classes.tsv"))? } else { Vec::new() };
        let flank20_classes_path = data.join("indel/flank20_read_classes.tsv");
        let classes20 = if selected_rows.iter().any(|(_, frame)| *frame == IndelFrame::Flank20)
            || (selected_rows.iter().any(|(_, frame)| matches!(frame, IndelFrame::CompleteShort | IndelFrame::LongEvent | IndelFrame::ExactJunction)) && flank20_classes_path.exists()) {
            homorepair_rows(&flank20_classes_path)?
        } else { Vec::new() };
        if !unmeasured_flank20.is_empty() && classes20.iter().any(|row|
            hr_field(row, "step").ok().and_then(|step| step.parse::<usize>().ok()).is_none_or(|step| step == 0)) {
            return Err(invalid("unmeasured 20 bp review has a read class without a valid step"));
        }
        let classes_complete = if selected_rows.iter().any(|(_, frame)| *frame == IndelFrame::CompleteShort) { homorepair_rows(&data.join("indel/short_read_classes.tsv"))? } else { Vec::new() };
        let classes_event = if selected_rows.iter().any(|(_, frame)| *frame == IndelFrame::LongEvent) { homorepair_rows(&data.join("indel/event_read_classes.tsv"))? } else { Vec::new() };
        let mut previous_step = 0;
        for (row, frame) in selected_rows {
            let parse = |field| hr_field(row, field)?.parse::<usize>()
                .map_err(|_| invalid(&format!("nonrepeat InDel {field} is invalid")));
            let step = parse("step")?;
            if step <= previous_step { return Err(invalid("nonrepeat InDel adopted steps are duplicate or out of order")); }
            let first_adoption = previous_step == 0;
            previous_step = step;
            let depth = parse("depth")?;
            let total = parse("total_ALT")?;
            let main = parse("main_ALT")?;
            let remainder = parse("remainder")?;
            let uncertain = parse("uncertain_reads")?;
            if depth == 0 || total > depth || remainder != depth - total
                || main <= remainder || main > total || uncertain > remainder {
                return Err(invalid("nonrepeat InDel strict support bound is not met"));
            }
            if matches!(frame, IndelFrame::CompleteShort | IndelFrame::LongEvent | IndelFrame::ExactJunction)
                && main as f64 / depth as f64 <= min_alt_fraction {
                return Err(invalid("InDel fallback support did not exceed the configured entry fraction"));
            }
            if matches!(frame, IndelFrame::LongEvent | IndelFrame::ExactJunction) {
                let entry = parse("entry_event_support")?;
                let source_depth = parse("source_depth")?;
                if entry < 3 || entry > source_depth || source_depth == 0
                    || entry as f64 / source_depth as f64 <= min_alt_fraction {
                    return Err(invalid("long InDel event did not satisfy its source entry bound"));
                }
            }
            let source_pos = parse("source_pos")?;
            let source_ref = hr_field(row, "source_ref")?;
            let source_start = source_pos.checked_sub(1).ok_or_else(|| invalid("nonrepeat InDel source position is zero"))?;
            let source_end = source_start.checked_add(source_ref.len()).ok_or_else(|| invalid("nonrepeat InDel source interval overflows"))?;
            if reconstructed.get(source_start..source_end) != Some(source_ref) {
                return Err(invalid("nonrepeat InDel source coordinate does not match the current sequence"));
            }
            let pos = parse("edit_pos")?;
            let edit_ref = hr_field(row, "edit_ref")?;
            let edit_alt = hr_field(row, "edit_alt")?;
            let edit_start = pos.checked_sub(1).ok_or_else(|| invalid("nonrepeat InDel edit has no anchor"))?;
            let edit_end = edit_start.checked_add(edit_ref.len()).ok_or_else(|| invalid("nonrepeat InDel edit coordinate overflow"))?;
            if edit_ref.is_empty() || edit_alt.is_empty()
                || reconstructed.get(edit_start..edit_end) != Some(edit_ref) {
                return Err(invalid("nonrepeat InDel edit differs from repeat-selected sequence"));
            }
            let old_len = reconstructed.len();
            let unchanged_prefix = reconstructed[..source_start].to_string();
            let unchanged_suffix = reconstructed[source_end..].to_string();
            reconstructed.replace_range(edit_start..edit_end, edit_alt);
            let complete_len = reconstructed.len().checked_add(source_ref.len())
                .and_then(|n| n.checked_sub(old_len))
                .ok_or_else(|| invalid("nonrepeat InDel complete ALT length is invalid"))?;
            let complete_end = source_start.checked_add(complete_len)
                .ok_or_else(|| invalid("nonrepeat InDel complete ALT interval overflows"))?;
            let complete_alt = reconstructed.get(source_start..complete_end)
                .ok_or_else(|| invalid("nonrepeat InDel complete ALT is outside the revised sequence"))?;
            if complete_alt.is_empty() || complete_alt == source_ref
                || reconstructed.get(..source_start) != Some(unchanged_prefix.as_str())
                || reconstructed.get(complete_end..) != Some(unchanged_suffix.as_str()) {
                return Err(invalid("nonrepeat InDel edit changes sequence outside its source allele"));
            }
            if matches!(frame, IndelFrame::LongEvent | IndelFrame::ExactJunction) && !long_indel_event_shape(source_ref, complete_alt) {
                return Err(invalid("long InDel event does not reconstruct a pure indel of at least 100 bp"));
            }
            let source_alt = hr_field(row, "source_alt")?;
            if let Some((_, failed20)) = unmeasured_flank20.get(&row_key(row)?) {
                let before_edit = format!("{unchanged_prefix}{source_ref}{unchanged_suffix}");
                if !unmeasured_flank20_template_failure_matches(hr_field(failed20, "reason")?,
                    &before_edit, source_start, source_ref, source_alt, &reconstructed, complete_alt) {
                    return Err(invalid("20 bp template failure is not reproduced by the source and ALT DNA"));
                }
            }
            let nested_edits = parse("nested_repeat_edits")?;
            if nested_edits > 1 || (nested_edits == 0 && complete_alt != source_alt) {
                return Err(invalid("nonrepeat InDel complete ALT disagrees with its recorded source or nested edit"));
            }
            let mut seen_ids = std::collections::HashSet::new();
            let mut class_counts = BTreeMap::<String, usize>::new();
            let mut uncertain_rows = 0usize;
            let step_text = step.to_string();
            let junction_classes;
            if frame == IndelFrame::ExactJunction {
                if nested_edits != 0 || complete_alt.len() <= source_ref.len() {
                    return Err(invalid("junction adoption is not the original single insertion"));
                }
                let initial_variants = polish_round_dir(case, round)
                    .join("03.validate/01.data/snv_indel_variants_combined.tsv");
                let initial = (first_adoption && edits.is_empty() && extra_adopted == 0)
                    .then_some(initial_variants.as_path());
                let target = read_single_fasta(&source)?.0;
                let target = target.split_whitespace().next().ok_or_else(|| invalid("source FASTA ID is empty"))?;
                let before_edit = format!("{unchanged_prefix}{source_ref}{unchanged_suffix}");
                let classes_path = commands::postrepair::junction::verify_automatic(
                    &data, step, target, &before_edit, &reconstructed,
                    &polish_round_dir(case, round).join("01.inputs/subgraph_reads.fastq.gz"), initial,
                )?;
                junction_classes = homorepair_rows(&classes_path)?;
            } else { junction_classes = Vec::new(); }
            let classes = match frame {
                IndelFrame::Flank60 => &classes60,
                IndelFrame::Flank20 => &classes20,
                IndelFrame::CompleteShort => &classes_complete,
                IndelFrame::LongEvent => &classes_event,
                IndelFrame::ExactJunction => &junction_classes,
            };
            let mut exact_by_id = BTreeMap::new();
            if matches!(frame, IndelFrame::LongEvent | IndelFrame::ExactJunction) {
                for exact in classes60.iter().filter(|class| hr_field(class, "step").ok() == Some(step_text.as_str())) {
                    let id = hr_field(exact, "read_id")?;
                    let allele = hr_field(exact, "exact_allele")?;
                    let label = hr_field(exact, "class")?;
                    if id.is_empty() || exact_by_id.insert(id, (label, allele)).is_some()
                        || !matches!((label, allele), ("uncertain", ".")) && !(
                            (label == "source" && allele == source_ref)
                            || (label == "ALT" && allele != source_ref && allele != "." && !allele.is_empty())) {
                        return Err(invalid("long InDel event has an invalid original exact read class"));
                    }
                }
            }
            for class in classes.iter().filter(|class| hr_field(class, "step").ok() == Some(step_text.as_str())) {
                let id = hr_field(class, "read_id")?;
                if id.is_empty() || !seen_ids.insert(id) { return Err(invalid("nonrepeat InDel read pool has an empty or duplicate identifier")); }
                if matches!(frame, IndelFrame::LongEvent | IndelFrame::ExactJunction) {
                    let exact = exact_by_id.get(id).copied();
                    match hr_field(class, "event_class")? {
                        "source" if !matches!(exact, Some(("ALT", _))) => {
                            *class_counts.entry(source_ref.into()).or_default() += 1;
                        }
                        "selected_ALT_event" if !matches!(exact, Some(("source", _))) => {
                            if let Some(("ALT", allele)) = exact {
                                if allele != complete_alt && !commands::postrepair::same_long_indel_event(source_ref, complete_alt, allele) {
                                    return Err(invalid("long InDel event conflicts with a different exact ALT"));
                                }
                            }
                            *class_counts.entry(complete_alt.into()).or_default() += 1;
                        }
                        "other_exact_ALT" => {
                            let Some(("ALT", allele)) = exact else {
                                return Err(invalid("long InDel competitor has no original exact ALT evidence"));
                            };
                            if allele == complete_alt || commands::postrepair::same_long_indel_event(source_ref, complete_alt, allele) {
                                return Err(invalid("long InDel event family was counted as a different exact ALT"));
                            }
                            *class_counts.entry(allele.into()).or_default() += 1;
                        }
                        "uncertain" => uncertain_rows += 1,
                        _ => return Err(invalid("long InDel event read class is invalid or conflicts with exact evidence")),
                    }
                    continue;
                }
                match hr_field(class, "class")? {
                    "source" | "ALT" => {
                        let allele = hr_field(class, "exact_allele")?;
                        *class_counts.entry(allele.into()).or_default() += 1;
                    }
                    "uncertain" if hr_field(class, "exact_allele")? == "." => uncertain_rows += 1,
                    _ => return Err(invalid("nonrepeat InDel read class is invalid")),
                }
            }
            if matches!(frame, IndelFrame::CompleteShort | IndelFrame::LongEvent | IndelFrame::ExactJunction) {
                for (label, fixed_classes) in [("60 bp", &classes60), ("20 bp", &classes20)] {
                    let mut fixed_ids = std::collections::HashSet::new();
                    for class in fixed_classes.iter().filter(|class| hr_field(class, "step").ok() == Some(step_text.as_str())) {
                        let id = hr_field(class, "read_id")?;
                        if !fixed_ids.insert(id) {
                            return Err(invalid(&format!("{label} fixed read pool has a duplicate identifier")));
                        }
                    }
                    if let Some((sixty, _)) = unmeasured_flank20.get(&row_key(row)?) {
                        if (label == "20 bp" && !fixed_ids.is_empty())
                            || (label == "60 bp" && hr_field(sixty, "depth")?.parse::<usize>().ok() != Some(fixed_ids.len())) {
                            return Err(invalid("unmeasured 20 bp review has typed reads or changed the original 60 bp pool"));
                        }
                    }
                    if (label == "60 bp" && fixed_ids.is_empty())
                        || (frame == IndelFrame::CompleteShort && !fixed_ids.is_empty() && fixed_ids != seen_ids)
                        || (matches!(frame, IndelFrame::LongEvent | IndelFrame::ExactJunction) && !fixed_ids.is_subset(&seen_ids)) {
                        return Err(invalid(&format!("InDel fallback review changed or lost the {label} fixed read pool")));
                    }
                }
            }
            let alt_counts: Vec<_> = class_counts.iter().filter(|(allele, _)| allele.as_str() != source_ref)
                .map(|(_, count)| *count).collect();
            if seen_ids.len() != depth || uncertain_rows != uncertain
                || alt_counts.iter().sum::<usize>() != total
                || class_counts.get(complete_alt).copied() != Some(main)
                || alt_counts.iter().filter(|count| **count == main).count() != 1
                || (matches!(frame, IndelFrame::LongEvent | IndelFrame::ExactJunction) && alt_counts.iter().any(|count| *count > main)) {
                return Err(invalid("nonrepeat InDel read classes do not reproduce the strict support bound"));
            }
            indel_frames.insert(step, reconstructed.clone());
        }
        for (step, kept) in kept_rows {
            let (previous, sequence) = indel_frames.range(..step).next_back()
                .ok_or_else(|| invalid("kept event has no current source frame"))?;
            commands::postrepair::verify_current_indel_group(&data, &polish_round_dir(case, round),
                step, previous + 1, sequence, edits.len() + extra_adopted == 0 && *previous == 0,
                min_alt_fraction, &kept.into_iter().cloned().collect::<Vec<_>>())?;
        }
    }
    let circle_adopted = replay_circle_end_edit_with_keep(&data, &outcome, &mut reconstructed,
        &polish_round_dir(case, round), &original, &reviews,
        edits.len() + extra_adopted + indel_adopted, min_alt_fraction, current_circle.first().copied())?;
    if edits.len() + extra_adopted + indel_adopted + circle_adopted != adopted || reconstructed != selected_sequence {
        return Err(invalid("selected sequence and applied candidate edits disagree"));
    }
    if keep.is_some() && !kept_event_verified {
        return Err(invalid("manual keep-source event was not the unique unresolved fixed-pool decision"));
    }
    let candidates = reviews.len();
    let support_status = value("support_status")?;
    let expected_support = if kept_ordinary > 0 || !deferred.is_empty() { "review_required" } else { "support_verified" };
    let circle_candidates = outcome.get("circle_end_candidate_count").map(String::as_str).unwrap_or("0")
        .parse::<usize>().map_err(|_| invalid("circle-end candidate count is invalid"))?;
    if support_status == "unchanged" {
        // The replay above permits zero edits only for explicitly unchanged
        // supplemental rows; no ordinary probe or other candidate may share this status.
        if kept_ordinary > 0 || candidates != 0 || extra_candidates == 0
            || indel_candidates != 0 || circle_candidates != 0
            || adopted != 0 || probes != 0 || local_batches != 0
            || selected_sequence != original
            || value("reason")? != "supplemental_repeat_candidates_unchanged"
            || value("failed_SV_frames")? != "."
            || ["pending_count", "uncertain_count", "indel_nested_repeat_edits", "indel_probe_count",
                "repeat_extra_adopted_count", "repeat_extra_review_count",
                "circle_end_adopted_count", "circle_end_review_count"].iter()
                .any(|key| outcome.get(*key).map(String::as_str) != Some("0"))
            || homorepair_rows(&data.join("repeat_extra/candidates.tsv"))?.iter()
                .any(|row| hr_field(row, "status").ok() != Some("unchanged"))
            || read_metric_value(&sv_summary_path(case, round), "status")?.as_deref() != Some("pass")
            || ordinary_validation_edit_count(&snv_indel_high_path(case, round))? != 0 {
            return Err(invalid("unchanged supplemental candidates are not the sole completed source-preserving decisions"));
        }
    } else if (candidates == 0 && extra_candidates == 0 && indel_candidates == 0 && circle_candidates == 0
            && (support_status != "no_candidates" || adopted != 0 || probes != 0))
        || (candidates == 0 && (extra_candidates > 0 || indel_candidates > 0 || circle_candidates > 0)
            && support_status != expected_support)
        || (candidates > 0 && support_status != expected_support) {
        return Err(invalid("support status conflicts with the candidate decisions"));
    }
    let final_status = value("final_validation_status")?;
    if adopted == 0 {
        if final_status != "skipped_unchanged" || reports.join("polished.fasta").exists()
            || ["final_validation_round", "final_validation_dir", "final_validation_fasta", "validation_round"]
                .iter().any(|key| outcome.contains_key(*key)) {
            return Err(invalid("unchanged source has unexpected revised output or final validation"));
        }
        return Ok(selected);
    }
    let final_round = round.checked_add(1).ok_or_else(|| invalid("round overflow"))?;
    let mut validation_case = case.clone();
    if let Some(out_dir) = &case.postrepair_out_dir {
        validation_case.polish_out_dir = out_dir.clone();
    }
    let case = &validation_case;
    if final_status != "passed" || value("final_validation_round")? != final_round.to_string()
        || value("validation_round")? != final_round.to_string()
        || fs::canonicalize(Path::new(&value("final_validation_dir")?))? != fs::canonicalize(polish_round_dir(case, final_round))? {
        return Err(invalid("changed sequence has not passed its next ordinary validation round"));
    }
    let final_fasta = polished_aln_path(case, final_round);
    if fs::canonicalize(Path::new(&value("final_validation_fasta")?))? != fs::canonicalize(&final_fasta)?
        || read_single_fasta(&final_fasta)?.1 != selected_sequence {
        return Err(invalid("ordinary validation input does not match the selected sequence"));
    }
    let report = homorepair_rows(&polish_report_path(case, final_round))?;
    for stage in ["sv-eval", "snv-indel-eval", "plot"] {
        let name = format!("{stage}:round_{final_round}");
        if !report.iter().any(|r| hr_field(r,"section").ok() == Some("stage")
            && hr_field(r,"name").ok() == Some(name.as_str())
            && hr_field(r,"metric").ok() == Some("status") && hr_field(r,"value").ok() == Some("ok")) {
            return Err(invalid("ordinary validation or plot stage is incomplete"));
        }
    }
    if read_metric_value(&sv_summary_path(case, final_round), "status")?.as_deref() != Some("pass")
        || ordinary_validation_edit_count(&snv_indel_high_path(case, final_round))? != 0 {
        return Err(invalid("selected sequence has unresolved SV or requires ordinary small-variant correction"));
    }
    Ok(final_fasta)
}

// Recount the failed fixed pools without turning a human decision into ALT support.
fn verify_manual_kept_indel(
    data: &Path, row: &BTreeMap<String, String>, flank20: &BTreeMap<String, String>,
    event: &ManualKeepEvent, original: &str, fraction: f64,
) -> Result<(), OrgraftError> {
    let invalid = |message: &str| OrgraftError::InvalidArgument(format!("manual keep-source {message}"));
    let number = |record: &BTreeMap<String, String>, key: &str| -> Result<usize, OrgraftError> {
        hr_field(record, key)?.parse().map_err(|_| invalid(&format!("invalid {key}")))
    };
    let start = event.pos.checked_sub(1).ok_or_else(|| invalid("position is zero"))?;
    let end = start.checked_add(event.reference.len()).ok_or_else(|| invalid("source interval overflows"))?;
    if event.reference.is_empty() || event.alternative.is_empty()
        || event.reference.len() == event.alternative.len()
        || !event.reference.bytes().chain(event.alternative.bytes()).all(|b| matches!(b,b'A'|b'C'|b'G'|b'T'|b'N'))
        || original.get(start..end) != Some(event.reference.as_str()) {
        return Err(invalid("event REF/ALT does not identify an InDel on the unchanged source"));
    }
    for (record, reason) in [(row, "highest_source_ALT_failed_fixed_pool_gate"), (flank20, "20bp_fixed_pool_gate_failed")] {
        if number(record,"step")? != 1 || number(record,"source_row")? == 0
            || number(record,"source_pos")? != event.pos || hr_field(record,"source_ref")? != event.reference
            || hr_field(record,"source_alt")? != event.alternative
            || hr_field(record,"status")? != "review_required" || hr_field(record,"reason")? != reason
            || number(record,"nested_repeat_edits")? != 0
            || number(record,"trial_candidates")? != 1 || number(record,"tested_probes")? != 1
            || ["edit_pos","edit_ref","edit_alt"].iter().any(|key| hr_field(record,key).ok()!=Some(".")) {
            return Err(invalid("requires one exact unresolved event with no applied or nested edit"));
        }
        let support=number(record,"source_support")?;
        let depth=number(record,"source_depth")?;
        let total=number(record,"source_total_ALT")?;
        if support<3 || depth==0 || support>total || total>depth || support as f64/depth as f64<=fraction {
            return Err(invalid("event does not reproduce its source entry bound"));
        }
    }
    let mut pools=Vec::new();
    for (record,relative) in [(row,"indel/read_classes.tsv"),(flank20,"indel/flank20_read_classes.tsv")] {
        let mut ids=BTreeSet::new();
        let mut counts=BTreeMap::<String,usize>::new();
        let mut unknown=0usize;
        for class in homorepair_rows(&data.join(relative))? {
            let id=hr_field(&class,"read_id")?;
            let allele=hr_field(&class,"exact_allele")?;
            if hr_field(&class,"step")?!="1" || id.is_empty() || !ids.insert(id.to_string()) {
                return Err(invalid("fixed read pool contains a foreign step or duplicate/empty ID"));
            }
            match hr_field(&class,"class")? {
                "source" if allele==event.reference => {*counts.entry(allele.into()).or_default()+=1;},
                "ALT" if allele!=event.reference && !allele.is_empty()
                    && allele.bytes().all(|b|matches!(b,b'A'|b'C'|b'G'|b'T'|b'N')) => {*counts.entry(allele.into()).or_default()+=1;},
                "uncertain" if allele=="." => unknown+=1,
                _ => return Err(invalid("fixed read class disagrees with its complete allele")),
            }
        }
        let depth=number(record,"depth")?;
        let total=counts.iter().filter(|(allele,_)|**allele!=event.reference).map(|(_,n)|*n).sum::<usize>();
        let main=counts.get(&event.alternative).copied().unwrap_or(0);
        if depth==0 || ids.len()!=depth || total>depth || number(record,"total_ALT")?!=total
            || number(record,"main_ALT")?!=main || number(record,"remainder")?!=depth-total
            || number(record,"uncertain_reads")?!=unknown || unknown>depth-total
            || main>depth-total || counts.iter().any(|(allele,n)| *allele!=event.reference && *n>main) {
            return Err(invalid("read classes do not reproduce the recorded failed fixed-pool bound"));
        }
        pools.push(ids);
    }
    if pools[0]!=pools[1] {return Err(invalid("20 bp review changed the original fixed read pool"));}
    Ok(())
}

fn homorepair_rows(path: &Path) -> Result<Vec<BTreeMap<String,String>>, OrgraftError> {
    let mut lines = BufReader::new(File::open(path)?).lines();
    let bad = || OrgraftError::InvalidArgument(format!("malformed homorepair table {}", path.display()));
    let header = lines.next().transpose()?.ok_or_else(bad)?;
    let keys: Vec<_> = header.split('\t').collect();
    if keys.is_empty() || keys.iter().collect::<std::collections::HashSet<_>>().len() != keys.len() { return Err(bad()); }
    let mut rows = Vec::new();
    for line in lines {
        let line = line?;
        if line.trim().is_empty() { continue; }
        let fields: Vec<_> = line.split('\t').collect();
        if fields.len() != keys.len() { return Err(bad()); }
        rows.push(keys.iter().zip(fields).map(|(k,v)| (k.to_string(),v.to_string())).collect());
    }
    Ok(rows)
}
fn hr_field<'a>(row: &'a BTreeMap<String,String>, name: &str) -> Result<&'a str, OrgraftError> {
    row.get(name).map(String::as_str).ok_or_else(|| OrgraftError::InvalidArgument(format!("homorepair table is missing {name}")))
}

fn write_homorepair_after_checkpoint(script: &mut String, config: &WorkflowConfig, case: &WorkflowCase, round: usize) -> Result<(), OrgraftError> {
    let args = homorepair_args(config, case, round).iter().map(|arg| shell_quote_str(arg)).collect::<Vec<_>>().join(" ");
    writeln!(script, "      # 05.postrepair: review residuals from the validated ordinary round.").unwrap();
    if declared_postrepair_review(config, case)?.is_some() {
        writeln!(script, "      {} workflow postrepair-run --config {} --case {} --round {round}",
            orgraft_shell_token(), shell_quote(&config.config_path), shell_quote_str(&case.name)).unwrap();
    } else if case.postrepair_out_dir.is_some() {
        let manual = homorepair_out_dir(case).parent().unwrap().join("manual_review.tsv");
        let quoted = shell_quote(&manual);
        writeln!(script, "      if [[ ! -e {quoted} && ! -L {quoted} ]]; then").unwrap();
        writeln!(script, "        {} postrepair {}", orgraft_shell_token(), args).unwrap();
        writeln!(script, "      fi").unwrap();
    } else {
        writeln!(script, "      {} postrepair {}", orgraft_shell_token(), args).unwrap();
    }
    let check = if case.postrepair_out_dir.is_some() { "postrepair-check" } else { "homorepair-check" };
    writeln!(script, "      final_polished=\"$({} workflow {check} --config {} --case {} --round {round})\"",
        orgraft_shell_token(), shell_quote(&config.config_path), shell_quote_str(&case.name)).unwrap();
    Ok(())
}

fn load_config_from_args(args: &[String]) -> Result<WorkflowConfig, OrgraftError> {
    let config_path = option_value(args, "--config")?.unwrap_or("orgraft.workflow.toml");
    WorkflowConfig::from_path(Path::new(config_path))
}

impl WorkflowConfig {
    fn from_path(path: &Path) -> Result<Self, OrgraftError> {
        let config_path = fs::canonicalize(path)?;
        let base = config_path.parent().ok_or_else(|| OrgraftError::InvalidArgument("configuration has no parent directory".into()))?;
        let raw = read_toml_like(&config_path)?;
        let project = raw.section("project");
        let software = raw.section("software");
        let workflow = raw.section("workflow");

        let sample = project
            .and_then(|section| section.get("sample"))
            .cloned()
            .filter(|sample| !sample.trim().is_empty())
            .unwrap_or_else(|| base.file_name().and_then(|n| n.to_str()).unwrap_or("sample").to_string());
        let results_dir_string = project
            .and_then(|section| section.get("results_dir"))
            .cloned()
            .unwrap_or_else(|| "results_workflow".to_string());
        let results_dir_string = config_relative_path(base, Path::new(&results_dir_string.replace("${sample}", &sample))).display().to_string();
        let source_results_dir_string = project
            .and_then(|section| section.get("source_results_dir"))
            .cloned()
            .unwrap_or_else(|| results_dir_string.clone());
        let source_results_dir_string = config_relative_path(base, Path::new(&source_results_dir_string.replace("${sample}", &sample))).display().to_string();
        let results_dir = PathBuf::from(expand_template_value(
            &results_dir_string,
            &sample,
            &results_dir_string,
            &source_results_dir_string,
            "",
            "",
        ));
        let soft_paths = PathBuf::from(expand_template_value(
            software
                .and_then(|section| section.get("soft_paths"))
                .map(String::as_str)
                .unwrap_or("soft_paths.txt"),
            &sample,
            &results_dir_string,
            &source_results_dir_string,
            "",
            "",
        ));
        let mode = workflow
            .and_then(|section| section.get("mode"))
            .cloned()
            .unwrap_or_else(|| "stepwise".to_string());
        let max_rounds = parse_config_usize(workflow, "max_rounds", 3)?.max(1);
        let threads = parse_config_usize(workflow, "threads", 64)?.max(1);
        let force = parse_config_bool(workflow, "force", false)?;
        let auto_sv_correction = parse_config_bool(workflow, "auto_sv_correction", true)?;
        let auto_snv_indel_correction =
            parse_config_bool(workflow, "auto_snv_indel_correction", true)?;
        let topology_simple_allowed_classes = workflow
            .and_then(|section| section.get("topology_simple_allowed_classes"))
            .map(|value| split_csv_set(value))
            .unwrap_or_else(default_simple_classes);
        let recruit = WorkflowRecruitConfig::from_section(
            raw.section("commands.recruit"),
            &sample,
            &results_dir_string,
            &source_results_dir_string,
        )?;
        let asm = WorkflowAsmConfig::from_section(
            raw.section("commands.asm"),
            &sample,
            &results_dir_string,
            &source_results_dir_string,
        )?;
        let resolve = WorkflowResolveConfig::from_section(
            raw.section("commands.resolve"),
            &sample,
            &results_dir_string,
            &source_results_dir_string,
        )?;
        let polish = WorkflowPolishConfig::from_section(
            raw.section("commands.polish"),
            &sample,
            &results_dir_string,
            &source_results_dir_string,
        )?;
        if raw.section("commands.postrepair.review").is_some() && raw.section("commands.postrepair").is_none() {
            return Err(OrgraftError::InvalidArgument("commands.postrepair.review requires commands.postrepair".into()));
        }
        let postrepair = raw.section("commands.postrepair").map(|section| {
            let frequency = parse_optional_f64(section, "frequency")?.unwrap_or(0.2);
            if !frequency.is_finite() || !(0.0..=1.0).contains(&frequency) {
                return Err(OrgraftError::InvalidArgument("commands.postrepair.frequency must be between 0 and 1".into()));
            }
            let expand = |value: &str| expand_template_value(value, &sample, &results_dir_string, &source_results_dir_string, "", "");
            Ok::<_, OrgraftError>(WorkflowPostrepairConfig {
                enabled: parse_section_bool(section, "enabled", true)?,
                out_dir: expanded_path_or(section, "out_dir", "${results_dir}/05.postrepair", &expand),
                threads: parse_optional_usize(section, "threads")?, frequency,
                review_files: raw.section("commands.postrepair.review").into_iter().flat_map(|section| section.iter())
                    .map(|(organelle, value)| {
                        if !matches!(organelle.as_str(), "mito" | "plastid") || value.trim().is_empty() {
                            return Err(OrgraftError::InvalidArgument("commands.postrepair.review requires a nonempty mito or plastid file".into()));
                        }
                        Ok((organelle.clone(), config_relative_path(base, Path::new(&expand(value)))))
                    }).collect::<Result<BTreeMap<_, _>, OrgraftError>>()?,
            })
        }).transpose()?;
        if (raw.section("commands.refine.graph_full_family").is_some() || raw.section("commands.refine.graph_measurement").is_some())
            && raw.section("commands.refine").is_none()
        {
            return Err(OrgraftError::InvalidArgument(
                "commands.refine.graph_full_family requires commands.refine".into(),
            ));
        }
        let refine = raw.section("commands.refine").map(|section| {
            let expand = |value: &str| expand_template_value(value, &sample, &results_dir_string, &source_results_dir_string, "", "");
            Ok::<_, OrgraftError>(WorkflowRefineConfig {
                enabled: parse_section_bool(section, "enabled", true)?,
                out_dir: expanded_path_or(section, "out_dir", "${results_dir}/07.refine", &expand),
                threads: parse_optional_usize(section, "threads")?,
                bundle: optional_expanded_path(section, "sample_bundle", &expand),
                organelle_reads: section.get("organelle_reads").map(|value| {
                    if value.trim().is_empty() {
                        return Err(OrgraftError::InvalidArgument("commands.refine.organelle_reads must not be empty".into()));
                    }
                    Ok(PathBuf::from(expand_template_value(value, &sample, &results_dir_string, &source_results_dir_string, "${organelle}", "")))
                }).transpose()?,
                intra_measurement: raw.section("commands.refine.intra_measurement").into_iter().flat_map(|section| section.iter())
                    .map(|(organelle, value)| {
                        if !matches!(organelle.as_str(), "mito" | "plastid") {
                            return Err(OrgraftError::InvalidArgument("commands.refine.intra_measurement requires mito or plastid".into()));
                        }
                        Ok((organelle.clone(), commands::refine::IntraMeasurement::parse(value)?))
                    }).collect::<Result<BTreeMap<_, _>, OrgraftError>>()?,
                duplicate_template_full_span: raw
                    .section("commands.refine.duplicate_template_full_span")
                    .into_iter()
                    .flat_map(|section| section.iter())
                    .map(|(organelle, value)| {
                        if !matches!(organelle.as_str(), "mito" | "plastid") {
                            return Err(OrgraftError::InvalidArgument(
                                "commands.refine.duplicate_template_full_span requires mito or plastid".into(),
                            ));
                        }
                        Ok((
                            organelle.clone(),
                            parse_bool_value(value, "commands.refine.duplicate_template_full_span")?,
                        ))
                    })
                    .collect::<Result<BTreeMap<_, _>, OrgraftError>>()?,
                linear_chain_flanks: raw
                    .section("commands.refine.linear_chain_flanks")
                    .into_iter()
                    .flat_map(|section| section.iter())
                    .map(|(organelle, value)| {
                        if !matches!(organelle.as_str(), "mito" | "plastid") {
                            return Err(OrgraftError::InvalidArgument(
                                "commands.refine.linear_chain_flanks requires mito or plastid".into(),
                            ));
                        }
                        Ok((
                            organelle.clone(),
                            parse_bool_value(value, "commands.refine.linear_chain_flanks")?,
                        ))
                    })
                    .collect::<Result<BTreeMap<_, _>, OrgraftError>>()?,
                graph_measurement: raw.section("commands.refine.graph_measurement").into_iter().flat_map(|section|section.iter())
                    .map(|(organelle,value)|{
                        if !matches!(organelle.as_str(),"mito"|"plastid") {return Err(OrgraftError::InvalidArgument("commands.refine.graph_measurement requires mito or plastid".into()));}
                        Ok((organelle.clone(),commands::refine::GraphMeasurement::parse(value)?))
                    }).collect::<Result<BTreeMap<_,_>,OrgraftError>>()?,
                graph_full_family: raw
                    .section("commands.refine.graph_full_family")
                    .into_iter()
                    .flat_map(|section| section.iter())
                    .map(|(organelle, value)| {
                        if !matches!(organelle.as_str(), "mito" | "plastid") {
                            return Err(OrgraftError::InvalidArgument(
                                "commands.refine.graph_full_family requires mito or plastid".into(),
                            ));
                        }
                        Ok((
                            organelle.clone(),
                            parse_bool_value(value, "commands.refine.graph_full_family")?,
                        ))
                    })
                    .collect::<Result<BTreeMap<_, _>, OrgraftError>>()?,
                adopt: raw.section("commands.refine.adopt").into_iter().flat_map(|section| section.iter())
                    .map(|(organelle, value)| {
                        if !matches!(organelle.as_str(), "mito" | "plastid") || value.trim().is_empty() {
                            return Err(OrgraftError::InvalidArgument("commands.refine.adopt requires a nonempty mito or plastid file".into()));
                        }
                        Ok((organelle.clone(), config_relative_path(base, Path::new(&expand(value)))))
                    }).collect::<Result<BTreeMap<_, _>, OrgraftError>>()?,
            })
        }).transpose()?;
        if let Some(refine) = &refine {
            for (organelle,mode) in &refine.graph_measurement {
                if *mode != commands::refine::GraphMeasurement::Default && (refine.graph_full_family.get(organelle)==Some(&true)
                    || refine.duplicate_template_full_span.get(organelle)==Some(&true) || refine.linear_chain_flanks.get(organelle)==Some(&true)) {
                    return Err(OrgraftError::InvalidArgument("explicit graph measurement conflicts with another graph template mode".into()));
                }
            }
            for (organelle, enabled) in &refine.graph_full_family {
                if *enabled
                    && (refine.duplicate_template_full_span.get(organelle) == Some(&true)
                        || refine.linear_chain_flanks.get(organelle) == Some(&true))
                {
                    return Err(OrgraftError::InvalidArgument(format!(
                        "graph_full_family cannot be combined with duplicate_template_full_span or linear_chain_flanks for {organelle}"
                    )));
                }
            }
        }
        if postrepair.as_ref().is_some_and(|p| p.enabled) && polish.per_read_variant_calls == Some(false) {
            return Err(OrgraftError::InvalidArgument("postrepair requires per_read_variant_calls".into()));
        }
        let rebuild = WorkflowRebuildConfig::from_section(
            raw.section("commands.rebuild"),
            &sample,
            &results_dir_string,
            &source_results_dir_string,
        )?;

        if refine.as_ref().is_some_and(|r| r.enabled) && !rebuild.enabled {
            return Err(OrgraftError::InvalidArgument("commands.refine requires commands.rebuild.enabled=true".into()));
        }
        let mut cases = Vec::new();
        for (section_name, section) in &raw.sections {
            if section_name == "workflow.case" || section_name.starts_with("workflow.case.") {
                let default_name = section_name
                    .strip_prefix("workflow.case.")
                    .unwrap_or("default");
                cases.push(WorkflowCase::from_section(
                    section,
                    default_name,
                    &sample,
                    &results_dir_string,
                    &source_results_dir_string,
                    &resolve.out_dir,
                    &polish.out_dir,
                    &rebuild.out_dir,
                    postrepair.as_ref().filter(|p| p.enabled).map(|p| &p.out_dir),
                )?);
            }
        }
        if cases.is_empty() {
            cases.push(WorkflowCase::from_section(
                &BTreeMap::new(),
                "default",
                &sample,
                &results_dir_string,
                &source_results_dir_string,
                &resolve.out_dir,
                &polish.out_dir,
                &rebuild.out_dir,
                postrepair.as_ref().filter(|p| p.enabled).map(|p| &p.out_dir),
            )?);
        }
        if refine.as_ref().is_some_and(|r| r.enabled) {
            let mut identities = BTreeSet::new();
            for case in cases.iter().filter(|c| c.enabled) {
                if case.sample != sample || !identities.insert((case.organelle.clone(), case.subgraph.clone())) {
                    return Err(OrgraftError::InvalidArgument("refine requires unique components from one project sample; use a separate config for each sample".into()));
                }
            }
        }
        for case in &mut cases {
            if case.reference.is_none() && case.pre_rotated_reference.is_none() {
                case.reference = recruit_bait_for_organelle(&recruit, &case.organelle);
            }
        }

        let mut config = Self {
            config_path: config_path.clone(),
            sample,
            results_dir,
            soft_paths,
            mode,
            max_rounds,
            threads,
            force,
            auto_sv_correction,
            auto_snv_indel_correction,
            recruit,
            asm,
            resolve,
            polish,
            postrepair,
            refine,
            rebuild,
            topology_simple_allowed_classes,
            cases,
        };
        config.resolve_relative_paths(base);
        validate_postrepair_sources(&config)?;
        declared_postrepair_reviews(&config)?;
        Ok(config)
    }

    fn resolve_relative_paths(&mut self, base: &Path) {
        let resolve = |path: &mut PathBuf| *path = config_relative_path(base, path);
        let optional = |path: &mut Option<PathBuf>| { if let Some(path) = path { resolve(path); } };
        resolve(&mut self.results_dir); resolve(&mut self.soft_paths);
        resolve(&mut self.recruit.out_dir); optional(&mut self.recruit.reads); optional(&mut self.recruit.sam);
        for (_, path) in &mut self.recruit.baits { resolve(path); }
        resolve(&mut self.asm.out_dir); optional(&mut self.asm.image_reference_fasta);
        resolve(&mut self.resolve.out_dir); resolve(&mut self.polish.out_dir);
        optional(&mut self.polish.sv_plot_highlight_read_ids);
        if let Some(postrepair) = &mut self.postrepair { resolve(&mut postrepair.out_dir); }
        if let Some(refine) = &mut self.refine { resolve(&mut refine.out_dir); optional(&mut refine.bundle); optional(&mut refine.organelle_reads); }
        resolve(&mut self.rebuild.out_dir);
        optional(&mut self.rebuild.edited_gfa); optional(&mut self.rebuild.polished_fasta);
        optional(&mut self.rebuild.image_reference_fasta); optional(&mut self.rebuild.merged_gfa_template);
        optional(&mut self.rebuild.minimap2); optional(&mut self.rebuild.blastn);
        for case in &mut self.cases {
            resolve(&mut case.draft_graph); optional(&mut case.unitig_graph); resolve(&mut case.checked_draft_gfa);
            optional(&mut case.reference); optional(&mut case.pre_rotated_reference); resolve(&mut case.reads);
            optional(&mut case.asm_reads); resolve(&mut case.resolve_out_dir); resolve(&mut case.polish_out_dir);
            optional(&mut case.postrepair_out_dir); resolve(&mut case.rebuild_out_dir);
            optional(&mut case.postrepair_source_dir); optional(&mut case.resolve_details); optional(&mut case.first_polish_report);
            optional(&mut case.rebuild_edited_gfa); optional(&mut case.rebuild_polished_fasta);
            optional(&mut case.image_reference_fasta); resolve(&mut case.workflow_dir);
            optional(&mut case.linearized_fasta); optional(&mut case.polish_reference);
        }
    }
}

fn config_relative_path(base: &Path, path: &Path) -> PathBuf {
    let path = if path.is_absolute() { path.to_path_buf() } else { base.join(path) };
    path.components().collect()
}

impl WorkflowCase {
    fn from_section(
        section: &BTreeMap<String, String>,
        default_name: &str,
        sample: &str,
        results_dir: &str,
        source_results_dir: &str,
        default_resolve_out_dir: &Path,
        default_polish_out_dir: &Path,
        default_rebuild_out_dir: &Path,
        default_postrepair_out_dir: Option<&PathBuf>,
    ) -> Result<Self, OrgraftError> {
        let case_sample = section
            .get("sample")
            .cloned()
            .unwrap_or_else(|| sample.to_string());
        let organelle = section
            .get("organelle")
            .cloned()
            .unwrap_or_else(|| "mito".to_string());
        let subgraph = section
            .get("subgraph")
            .cloned()
            .unwrap_or_else(|| "subgraph_001".to_string());
        let name = section
            .get("name")
            .cloned()
            .unwrap_or_else(|| default_name.to_string());
        let enabled = parse_section_bool(section, "enabled", true)?;

        let expand = |value: &str| {
            expand_template_value(
                value,
                &case_sample,
                results_dir,
                source_results_dir,
                &organelle,
                &subgraph,
            )
        };
        let path_or = |key: &str, default: String| {
            section
                .get(key)
                .map(|value| PathBuf::from(expand(value)))
                .unwrap_or_else(|| PathBuf::from(expand(&default)))
        };
        let optional_path = |key: &str| section.get(key).map(|value| PathBuf::from(expand(value)));

        let workflow_dir = path_or(
            "workflow_dir",
            "${results_dir}/workflow/${organelle}/${subgraph}".to_string(),
        );
        let draft_graph = path_or(
            "draft_graph",
            "${results_dir}/02.draft_asm/${organelle}/03.finalize_graph/graph.gfa".to_string(),
        );
        let unitig_graph = optional_path("unitig_graph");
        let checked_draft_gfa = path_or(
            "checked_draft_gfa",
            "${results_dir}/workflow/${organelle}/${subgraph}/checkpoint_1/checked_draft.gfa"
                .to_string(),
        );
        let reads = path_or(
            "reads",
            "${results_dir}/01.recruit/${organelle}.fastq.gz".to_string(),
        );
        let asm_reads = optional_path("asm_reads");
        let resolve_out_dir = path_or(
            "resolve_out_dir",
            default_resolve_out_dir.display().to_string(),
        );
        let polish_out_dir = path_or(
            "polish_out_dir",
            default_polish_out_dir.display().to_string(),
        );
        let postrepair_out_dir = optional_path("postrepair_out_dir").or_else(|| default_postrepair_out_dir.cloned());
        let postrepair_source_dir = optional_path("postrepair_source_dir");
        let resolve_details = optional_path("resolve_details");
        let first_polish_report = optional_path("first_polish_report");
        let rebuild_out_dir = path_or("rebuild_out_dir", default_rebuild_out_dir.join(&organelle).display().to_string());
        let rebuild_edited_gfa = optional_path("rebuild_edited_gfa");
        let rebuild_polished_fasta = optional_path("rebuild_polished_fasta");
        let image_reference_fasta = optional_path("image_reference_fasta");
        let reference = optional_path("reference");
        let pre_rotated_reference = optional_path("pre_rotated_reference");
        let linearized_fasta = optional_path("linearized_fasta");
        let polish_reference = optional_path("polish_reference");
        let sv_correction_subgroup = section.get("sv_correction_subgroup").cloned();
        let refine_downgrade_representation = commands::refine::DowngradeRepresentation::parse(
            section.get("refine_downgrade_representation").map(String::as_str).unwrap_or("recorded"),
        )?;

        Ok(Self {
            enabled,
            name,
            sample: case_sample,
            organelle,
            subgraph,
            draft_graph,
            unitig_graph,
            checked_draft_gfa,
            reference,
            pre_rotated_reference,
            reads,
            asm_reads,
            resolve_out_dir,
            polish_out_dir,
            postrepair_out_dir,
            postrepair_source_dir,
            resolve_details,
            first_polish_report,
            rebuild_out_dir,
            rebuild_edited_gfa,
            rebuild_polished_fasta,
            image_reference_fasta,
            workflow_dir,
            linearized_fasta,
            polish_reference,
            sv_correction_subgroup,
            refine_downgrade_representation,
        })
    }
}

impl WorkflowRecruitConfig {
    fn from_section(
        section: Option<&BTreeMap<String, String>>,
        sample: &str,
        results_dir: &str,
        source_results_dir: &str,
    ) -> Result<Self, OrgraftError> {
        let section = section.cloned().unwrap_or_default();
        let expand = |value: &str| {
            expand_template_value(value, sample, results_dir, source_results_dir, "", "")
        };
        Ok(Self {
            enabled: parse_section_bool(&section, "enabled", false)?,
            reads: optional_expanded_path(&section, "raw_reads", &expand)
                .or_else(|| optional_expanded_path(&section, "reads", &expand)),
            out_dir: expanded_path_or(&section, "out_dir", "${results_dir}/01.recruit", &expand),
            threads: parse_optional_usize(&section, "threads")?,
            baits: parse_baits(section.get("baits"), &expand),
            prefix: section.get("prefix").cloned(),
            bait_format: section.get("bait_format").cloned(),
            gfa_split: section.get("gfa_split").cloned(),
            rename_bait: parse_section_bool(&section, "rename_bait", false)?,
            write_id_map: parse_section_bool(&section, "write_id_map", false)?,
            split_output: section.get("split_output").cloned(),
            gzip_output: section
                .get("gzip_output")
                .map(|value| parse_bool_value(value, "commands.recruit.gzip_output"))
                .transpose()?,
            minimap2: section.get("minimap2").cloned(),
            align_mode: section.get("align_mode").cloned(),
            platform: section.get("platform").cloned(),
            preset: section.get("preset").cloned(),
            min_mapq: parse_optional_u8(&section, "min_mapq")?,
            min_aln_len: parse_optional_u64(&section, "min_aln_len")?,
            sam: optional_expanded_path(&section, "sam", &expand),
            max_reads: split_semicolon(section.get("max_reads")),
            random_seed: parse_optional_u64(&section, "random_seed")?,
            write_sampled_ids: parse_section_bool(&section, "write_sampled_ids", false)?,
            read_stats: section.get("read_stats").cloned(),
            write_read_classification: parse_section_bool(
                &section,
                "write_read_classification",
                false,
            )?,
            write_bait_partitions: parse_section_bool(&section, "write_bait_partitions", false)?,
            gzip_tool: section.get("gzip_tool").cloned(),
            mode: section.get("mode").cloned(),
            iterations: parse_optional_usize(&section, "iterations")?,
            extra_args: split_extra_args(section.get("extra_args")),
        })
    }
}

impl WorkflowAsmConfig {
    fn from_section(
        section: Option<&BTreeMap<String, String>>,
        sample: &str,
        results_dir: &str,
        source_results_dir: &str,
    ) -> Result<Self, OrgraftError> {
        let section = section.cloned().unwrap_or_default();
        let expand = |value: &str| {
            expand_template_value(value, sample, results_dir, source_results_dir, "", "")
        };
        Ok(Self {
            enabled: parse_section_bool(&section, "enabled", false)?,
            out_dir: expanded_path_or(&section, "out_dir", "${results_dir}/02.draft_asm", &expand),
            threads: parse_optional_usize(&section, "threads")?,
            profile: section.get("profile").cloned(),
            stable: parse_section_bool(&section, "stable", false)?,
            min_graph_coverage: parse_optional_u32(&section, "min_graph_coverage")?,
            min_branch_ratio: parse_optional_f64(&section, "min_branch_ratio")?,
            min_tip_len: parse_optional_usize(&section, "min_tip_len")?,
            min_link_support: parse_optional_u32(&section, "min_link_support")?,
            min_link_ratio: parse_optional_f64(&section, "min_link_ratio")?,
            subsets: section.get("subsets").cloned(),
            image_reference_fasta: optional_expanded_path(
                &section,
                "image_reference_fasta",
                &expand,
            ),
            keep_debug_files: parse_section_bool(&section, "keep_debug_files", false)?,
            extra_args: split_extra_args(section.get("extra_args")),
        })
    }
}

impl WorkflowResolveConfig {
    fn from_section(
        section: Option<&BTreeMap<String, String>>,
        sample: &str,
        results_dir: &str,
        source_results_dir: &str,
    ) -> Result<Self, OrgraftError> {
        let section = section.cloned().unwrap_or_default();
        let expand = |value: &str| {
            expand_template_value(value, sample, results_dir, source_results_dir, "", "")
        };
        Ok(Self {
            out_dir: expanded_path_or(
                &section,
                "out_dir",
                "${results_dir}/03.resolve_gfa",
                &expand,
            ),
            gfa_editor_mode: section.get("gfa_editor_mode").cloned(),
            max_states: parse_optional_usize(&section, "max_states")?,
            max_candidates: parse_optional_usize(&section, "max_candidates")?,
            extra_args: split_extra_args(section.get("extra_args")),
        })
    }
}

impl WorkflowPolishConfig {
    fn from_section(
        section: Option<&BTreeMap<String, String>>,
        sample: &str,
        results_dir: &str,
        source_results_dir: &str,
    ) -> Result<Self, OrgraftError> {
        let section = section.cloned().unwrap_or_default();
        let expand = |value: &str| {
            expand_template_value(value, sample, results_dir, source_results_dir, "", "")
        };
        let homorepair = match section.get("homorepair").map(String::as_str).unwrap_or("0") {
            "0" => false,
            "1" => true,
            _ => return Err(OrgraftError::InvalidArgument("commands.polish.homorepair expects 0 or 1".into())),
        };
        let homorepair_frq = parse_optional_f64(&section, "homorepair_frq")?.unwrap_or(0.2);
        if !homorepair_frq.is_finite() || !(0.0..=1.0).contains(&homorepair_frq) {
            return Err(OrgraftError::InvalidArgument("commands.polish.homorepair_frq expects a finite number from 0 to 1".into()));
        }
        let extra_args = split_extra_args(section.get("extra_args"));
        if extra_args.iter().any(|arg| arg == "-hr" || arg.starts_with("-hr=") || arg.starts_with("--homorepair")) {
            return Err(OrgraftError::InvalidArgument("configure homorepair in commands.polish, not extra_args; it runs only after the final checkpoint".into()));
        }
        if homorepair && parse_optional_bool(&section, "per_read_variant_calls")? == Some(false) {
            return Err(OrgraftError::InvalidArgument("homorepair requires per_read_variant_calls".into()));
        }
        Ok(Self {
            homorepair,
            homorepair_frq,
            out_dir: expanded_path_or(&section, "out_dir", "${results_dir}/04.polish", &expand),
            threads: parse_optional_usize(&section, "threads")?,
            per_read_variant_calls: parse_optional_bool(&section, "per_read_variant_calls")?,
            snv_indel_overlap_policy: section.get("snv_indel_overlap_policy").cloned(),
            plot_range: section.get("plot_range").cloned(),
            plot_dpi: parse_optional_usize(&section, "plot_dpi")?,
            plot_output_format: section.get("plot_output_format").cloned(),
            coverage_plot_rasterize: parse_optional_bool(&section, "coverage_plot_rasterize")?,
            snv_indel_plot_rasterize: parse_optional_bool(&section, "snv_indel_plot_rasterize")?,
            sv_plot_highlight_subgroups: section.get("sv_plot_highlight_subgroups").cloned(),
            sv_plot_highlight_read_ids: optional_expanded_path(
                &section,
                "sv_plot_highlight_read_ids",
                &expand,
            ),
            sv_plot_highlight_min_fraction: parse_optional_f64(
                &section,
                "sv_plot_highlight_min_fraction",
            )?,
            sv_plot_highlight_min_reads: parse_optional_usize(
                &section,
                "sv_plot_highlight_min_reads",
            )?,
            snv_indel_plot_low_confidence: section.get("snv_indel_plot_low_confidence").cloned(),
            snv_indel_plot_low_min_reads: parse_optional_usize(
                &section,
                "snv_indel_plot_low_min_reads",
            )?,
            snv_indel_plot_low_min_fraction: parse_optional_f64(
                &section,
                "snv_indel_plot_low_min_fraction",
            )?,
            snv_indel_plot_high_risk_fraction: parse_optional_f64(
                &section,
                "snv_indel_plot_high_risk_fraction",
            )?,
            extra_args,
        })
    }
}

impl WorkflowRebuildConfig {
    fn from_section(
        section: Option<&BTreeMap<String, String>>,
        sample: &str,
        results_dir: &str,
        source_results_dir: &str,
    ) -> Result<Self, OrgraftError> {
        let section = section.cloned().unwrap_or_default();
        let expand = |value: &str| {
            expand_template_value(value, sample, results_dir, source_results_dir, "", "")
        };
        Ok(Self {
            enabled: parse_section_bool(&section, "enabled", false)?,
            out_dir: expanded_path_or(&section, "out_dir", "${results_dir}/05.rebuild", &expand),
            threads: parse_optional_usize(&section, "threads")?,
            edited_gfa: optional_expanded_path(&section, "edited_gfa", &expand),
            polished_fasta: optional_expanded_path(&section, "polished_fasta", &expand),
            image_reference_fasta: optional_expanded_path(
                &section,
                "image_reference_fasta",
                &expand,
            ),
            merged_gfa_template: optional_expanded_path(&section, "merged_gfa_template", &expand),
            minimap2: optional_expanded_path(&section, "minimap2", &expand),
            blastn: optional_expanded_path(&section, "blastn", &expand),
            keep_debug: parse_section_bool(&section, "keep_debug", false)?,
            extra_args: split_extra_args(section.get("extra_args")),
        })
    }
}

#[derive(Debug, Default)]
struct RawConfig {
    sections: BTreeMap<String, BTreeMap<String, String>>,
}

impl RawConfig {
    fn section(&self, name: &str) -> Option<&BTreeMap<String, String>> {
        self.sections.get(name)
    }
}

fn read_toml_like(path: &Path) -> Result<RawConfig, OrgraftError> {
    let file = File::open(path)?;
    let mut current = String::new();
    let mut raw = RawConfig::default();

    for (line_index, line_result) in BufReader::new(file).lines().enumerate() {
        let line_number = line_index + 1;
        let line = strip_inline_comment(&line_result?);
        let line = line.trim();
        if line.is_empty() {
            continue;
        }
        if line.starts_with('[') && line.ends_with(']') {
            current = line
                .trim_start_matches('[')
                .trim_end_matches(']')
                .trim()
                .to_string();
            raw.sections.entry(current.clone()).or_default();
            continue;
        }
        let Some((key, value)) = line.split_once('=') else {
            return Err(OrgraftError::InvalidArgument(format!(
                "{}:{line_number}: expected key = value",
                path.display()
            )));
        };
        let key = key.trim();
        if key.is_empty() {
            return Err(OrgraftError::InvalidArgument(format!(
                "{}:{line_number}: empty config key",
                path.display()
            )));
        }
        let value = parse_scalar_value(value.trim());
        raw.sections
            .entry(current.clone())
            .or_default()
            .insert(key.to_string(), value);
    }

    Ok(raw)
}

fn strip_inline_comment(line: &str) -> String {
    let mut in_string = false;
    let mut escaped = false;
    let mut output = String::new();
    for character in line.chars() {
        if escaped {
            output.push(character);
            escaped = false;
            continue;
        }
        match character {
            '\\' if in_string => {
                output.push(character);
                escaped = true;
            }
            '"' => {
                in_string = !in_string;
                output.push(character);
            }
            '#' if !in_string => break,
            _ => output.push(character),
        }
    }
    output
}

fn parse_scalar_value(value: &str) -> String {
    let trimmed = value.trim();
    if trimmed.starts_with('"') && trimmed.ends_with('"') && trimmed.len() >= 2 {
        unescape_toml_string(&trimmed[1..trimmed.len() - 1])
    } else {
        trimmed.to_string()
    }
}

fn unescape_toml_string(value: &str) -> String {
    let mut output = String::new();
    let mut escaped = false;
    for character in value.chars() {
        if escaped {
            match character {
                'n' => output.push('\n'),
                'r' => output.push('\r'),
                't' => output.push('\t'),
                '\\' => output.push('\\'),
                '"' => output.push('"'),
                other => output.push(other),
            }
            escaped = false;
        } else if character == '\\' {
            escaped = true;
        } else {
            output.push(character);
        }
    }
    output
}

fn parse_config_usize(
    section: Option<&BTreeMap<String, String>>,
    key: &str,
    default: usize,
) -> Result<usize, OrgraftError> {
    let Some(value) = section.and_then(|section| section.get(key)) else {
        return Ok(default);
    };
    value.parse::<usize>().map_err(|_| {
        OrgraftError::InvalidArgument(format!("workflow.{key} expects a positive integer"))
    })
}

fn parse_config_bool(
    section: Option<&BTreeMap<String, String>>,
    key: &str,
    default: bool,
) -> Result<bool, OrgraftError> {
    let Some(value) = section.and_then(|section| section.get(key)) else {
        return Ok(default);
    };
    match value {
        value if value.eq_ignore_ascii_case("true") => Ok(true),
        value if value.eq_ignore_ascii_case("false") => Ok(false),
        _ => Err(OrgraftError::InvalidArgument(format!(
            "workflow.{key} expects true or false"
        ))),
    }
}

fn parse_section_bool(
    section: &BTreeMap<String, String>,
    key: &str,
    default: bool,
) -> Result<bool, OrgraftError> {
    let Some(value) = section.get(key) else {
        return Ok(default);
    };
    parse_bool_value(value, key)
}

fn parse_bool_value(value: &str, key: &str) -> Result<bool, OrgraftError> {
    match value {
        value if value.eq_ignore_ascii_case("true") => Ok(true),
        value if value.eq_ignore_ascii_case("false") => Ok(false),
        value if value.eq_ignore_ascii_case("on") => Ok(true),
        value if value.eq_ignore_ascii_case("off") => Ok(false),
        _ => Err(OrgraftError::InvalidArgument(format!(
            "{key} expects true/false or on/off"
        ))),
    }
}

fn parse_optional_bool(
    section: &BTreeMap<String, String>,
    key: &str,
) -> Result<Option<bool>, OrgraftError> {
    section
        .get(key)
        .map(|value| parse_bool_value(value, key))
        .transpose()
}

fn parse_optional_usize(
    section: &BTreeMap<String, String>,
    key: &str,
) -> Result<Option<usize>, OrgraftError> {
    section
        .get(key)
        .map(|value| {
            value.parse::<usize>().map_err(|_| {
                OrgraftError::InvalidArgument(format!("{key} expects a positive integer"))
            })
        })
        .transpose()
}

fn parse_optional_u32(
    section: &BTreeMap<String, String>,
    key: &str,
) -> Result<Option<u32>, OrgraftError> {
    section
        .get(key)
        .map(|value| {
            value
                .parse::<u32>()
                .map_err(|_| OrgraftError::InvalidArgument(format!("{key} expects an integer")))
        })
        .transpose()
}

fn parse_optional_u8(
    section: &BTreeMap<String, String>,
    key: &str,
) -> Result<Option<u8>, OrgraftError> {
    section
        .get(key)
        .map(|value| {
            value
                .parse::<u8>()
                .map_err(|_| OrgraftError::InvalidArgument(format!("{key} expects 0-255")))
        })
        .transpose()
}

fn parse_optional_u64(
    section: &BTreeMap<String, String>,
    key: &str,
) -> Result<Option<u64>, OrgraftError> {
    section
        .get(key)
        .map(|value| {
            value
                .parse::<u64>()
                .map_err(|_| OrgraftError::InvalidArgument(format!("{key} expects an integer")))
        })
        .transpose()
}

fn parse_optional_f64(
    section: &BTreeMap<String, String>,
    key: &str,
) -> Result<Option<f64>, OrgraftError> {
    section
        .get(key)
        .map(|value| {
            value
                .parse::<f64>()
                .map_err(|_| OrgraftError::InvalidArgument(format!("{key} expects a number")))
        })
        .transpose()
}

fn expanded_path_or<F>(
    section: &BTreeMap<String, String>,
    key: &str,
    default: &str,
    expand: F,
) -> PathBuf
where
    F: Fn(&str) -> String,
{
    section
        .get(key)
        .map(|value| PathBuf::from(expand(value)))
        .unwrap_or_else(|| PathBuf::from(expand(default)))
}

fn optional_expanded_path<F>(
    section: &BTreeMap<String, String>,
    key: &str,
    expand: F,
) -> Option<PathBuf>
where
    F: Fn(&str) -> String,
{
    section.get(key).map(|value| PathBuf::from(expand(value)))
}

fn parse_baits<F>(value: Option<&String>, expand: F) -> Vec<(String, PathBuf)>
where
    F: Fn(&str) -> String,
{
    value
        .map(|value| {
            value
                .split(',')
                .filter_map(|entry| {
                    let (label, path) = entry.split_once('=')?;
                    let label = label.trim();
                    let path = path.trim();
                    if label.is_empty() || path.is_empty() {
                        None
                    } else {
                        Some((label.to_string(), PathBuf::from(expand(path))))
                    }
                })
                .collect()
        })
        .unwrap_or_default()
}

fn recruit_bait_for_organelle(recruit: &WorkflowRecruitConfig, organelle: &str) -> Option<PathBuf> {
    let wanted = organelle.trim();
    recruit
        .baits
        .iter()
        .find(|(label, _)| bait_label_matches_organelle(label, wanted))
        .map(|(_, path)| path.clone())
}

fn bait_label_matches_organelle(label: &str, organelle: &str) -> bool {
    let label = label.trim();
    let organelle = organelle.trim();
    if label == organelle {
        return true;
    }
    matches!(
        (label, organelle),
        ("plasti", "plastid") | ("plastid", "plasti")
    )
}

fn split_semicolon(value: Option<&String>) -> Vec<String> {
    value
        .map(|value| {
            value
                .split(';')
                .map(str::trim)
                .filter(|entry| !entry.is_empty())
                .map(ToString::to_string)
                .collect()
        })
        .unwrap_or_default()
}

fn split_extra_args(value: Option<&String>) -> Vec<String> {
    value
        .map(|value| value.split_whitespace().map(ToString::to_string).collect())
        .unwrap_or_default()
}

fn expand_template_value(
    value: &str,
    sample: &str,
    results_dir: &str,
    source_results_dir: &str,
    organelle: &str,
    subgraph: &str,
) -> String {
    value
        .replace("${sample}", sample)
        .replace("${results_dir}", results_dir)
        .replace("${source_results_dir}", source_results_dir)
        .replace("${organelle}", organelle)
        .replace("${subgraph}", subgraph)
}

fn select_case<'a>(
    config: &'a WorkflowConfig,
    args: &[String],
) -> Result<&'a WorkflowCase, OrgraftError> {
    if let Some(name) = option_value(args, "--case")? {
        let case = config.cases.iter().find(|case| case.name == name).ok_or_else(|| {
            OrgraftError::InvalidArgument(format!("workflow case `{name}` was not found"))
        })?;
        if !case.enabled {
            let mut selected_config = config.clone();
            selected_config.cases.iter_mut().find(|selected| selected.name == name).unwrap().enabled = true;
            validate_postrepair_sources(&selected_config)?;
        }
        return Ok(case);
    }
    let enabled = enabled_cases(config);
    if enabled.len() == 1 {
        Ok(enabled[0])
    } else if enabled.is_empty() {
        Err(OrgraftError::InvalidArgument(
            "no enabled workflow cases are configured".to_string(),
        ))
    } else {
        Err(OrgraftError::InvalidArgument(
            "multiple enabled workflow cases are configured; pass --case NAME".to_string(),
        ))
    }
}

fn enabled_cases(config: &WorkflowConfig) -> Vec<&WorkflowCase> {
    config.cases.iter().filter(|case| case.enabled).collect()
}

fn parse_round(args: &[String]) -> Result<usize, OrgraftError> {
    option_value(args, "--round")?
        .unwrap_or("1")
        .parse::<usize>()
        .map_err(|_| {
            OrgraftError::InvalidArgument("--round expects a positive integer".to_string())
        })
}

fn missing_link_segments(path: &Path) -> Result<BTreeSet<String>, OrgraftError> {
    let mut segments = BTreeSet::new();
    let mut linked = BTreeSet::new();
    for (line_index, line_result) in BufReader::new(File::open(path)?).lines().enumerate() {
        let line_number = line_index + 1;
        let line = line_result?;
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        match fields.first().copied() {
            Some("S") => {
                let segment_id = fields.get(1).ok_or_else(|| {
                    OrgraftError::InvalidArgument(format!(
                        "GFA line {line_number}: segment record is missing an id"
                    ))
                })?;
                segments.insert((*segment_id).to_string());
            }
            Some("L") => {
                let from = fields.get(1).ok_or_else(|| {
                    OrgraftError::InvalidArgument(format!(
                        "GFA line {line_number}: link record is missing from segment"
                    ))
                })?;
                let to = fields.get(3).ok_or_else(|| {
                    OrgraftError::InvalidArgument(format!(
                        "GFA line {line_number}: link record is missing to segment"
                    ))
                })?;
                linked.insert((*from).to_string());
                linked.insert((*to).to_string());
            }
            _ => {}
        }
    }
    Ok(linked.difference(&segments).cloned().collect())
}

fn complex_nodes(report: &TopologyReport, allowed: &BTreeSet<String>) -> Vec<String> {
    report
        .nodes
        .iter()
        .filter(|node| !allowed.contains(node.taxon.code))
        .map(|node| format!("{}:{}", node.node_id, node.taxon.code))
        .collect()
}

fn default_simple_classes() -> BTreeSet<String> {
    split_csv_set("0-0,0-1/1-0,1-1,1-2/2-1,2-2")
}

fn split_csv_set(value: &str) -> BTreeSet<String> {
    value
        .split(',')
        .map(str::trim)
        .filter(|value| !value.is_empty())
        .map(ToString::to_string)
        .collect()
}

#[derive(Debug, Clone)]
struct RuntimeStageRow {
    scope: String,
    stage: String,
    seconds: Option<f64>,
    source: String,
}

fn workflow_runtime_summary_markdown(config: &WorkflowConfig) -> Result<String, OrgraftError> {
    let mut out = String::new();
    writeln!(out, "# OrgRAFT workflow runtime summary").unwrap();
    writeln!(out).unwrap();
    writeln!(out, "Generated by `orgraft workflow runtime-summary`.").unwrap();
    writeln!(out).unwrap();
    writeln!(out, "- Config: `{}`", config.config_path.display()).unwrap();
    writeln!(out, "- Results dir: `{}`", config.results_dir.display()).unwrap();
    writeln!(out, "- Sample: `{}`", config.sample).unwrap();
    writeln!(out).unwrap();

    write_runtime_case_statuses(&mut out, config)?;
    write_runtime_repair_results(&mut out, config)?;
    write_runtime_homorepair_statuses(&mut out, config)?;
    write_runtime_stage_table(&mut out, config)?;
    write_runtime_finalize_table(&mut out, config)?;
    write_runtime_fake_validate(&mut out, config)?;

    writeln!(out, "## Notes").unwrap();
    writeln!(
        out,
        "- Timing rows are parsed from existing logs/reports; missing outputs are shown as `pending`."
    )
    .unwrap();
    writeln!(
        out,
        "- Recruit timing uses the latest `minimap2.stderr.log` real time when present."
    )
    .unwrap();
    writeln!(
        out,
        "- `03.finalize_graph` link cleanup is reported from each asm finalize manifest."
    )
    .unwrap();

    Ok(out)
}

fn write_runtime_case_statuses(
    out: &mut String,
    config: &WorkflowConfig,
) -> Result<(), OrgraftError> {
    writeln!(out, "## Case Status").unwrap();
    writeln!(
        out,
        "| case | organelle | subgraph | checkpoint1 | checkpoint2 | rebuild |"
    )
    .unwrap();
    writeln!(out, "| --- | --- | --- | --- | --- | --- |").unwrap();
    for case in &config.cases {
        let checkpoint1 = if case.postrepair_source_dir.is_some() { "not_run (source round input)".into() } else { optional_metric_value(
            &case
                .workflow_dir
                .join("checkpoint_1")
                .join("checkpoint_1.status.tsv"),
            "status",
        )?
        .unwrap_or_else(|| "pending".to_string()) };
        let checkpoint2 = if case.postrepair_source_dir.is_some() { "not_run (source round input)".into() } else { latest_checkpoint2_status(config, case)? };
        let rebuild = rebuild_status(case)?.unwrap_or_else(|| "pending".to_string());
        writeln!(
            out,
            "| `{}` | {} | {} | {} | {} | {} |",
            case.name, case.organelle, case.subgraph, checkpoint1, checkpoint2, rebuild
        )
        .unwrap();
    }
    writeln!(out).unwrap();
    Ok(())
}

fn write_runtime_repair_results(out: &mut String, config: &WorkflowConfig) -> Result<(), OrgraftError> {
    if config.postrepair.as_ref().is_some_and(|repair| repair.enabled) {
        writeln!(out, "## Sequence repair\n").unwrap();
        writeln!(out, "Counts describe the initial review; any additional repair method and its evidence are linked separately.\n").unwrap();
        writeln!(out, "| case | status | repeat candidates | InDel candidates | initial edits | unresolved | selection |\n| --- | --- | ---: | ---: | ---: | ---: | --- |").unwrap();
        for case in config.cases.iter().filter(|case| case.enabled) {
            let Some(root) = &case.postrepair_out_dir else { continue; };
            let component = root.join(&case.organelle).join(&case.subgraph);
            let summary = component.join("postrepair.tsv");
            let run = homorepair_outcome_path(case);
            let values = read_homorepair_summary_outcome(&run)?;
            let metric = |key: &str| values.get(key).cloned().unwrap_or_else(|| "not recorded".into());
            let mut candidates = "not recorded".to_string();
            if run.is_file() {
                for row in homorepair_rows(&run)? {
                    if row.get("section").map(String::as_str) == Some("homorepair")
                        && row.get("metric").map(String::as_str) == Some("candidate_count") {
                        candidates = hr_field(&row, "value")?.to_string();
                    }
                }
            }
            let status = optional_metric_value(&summary, "status")?.unwrap_or_else(|| "pending".into());
            let unresolved = optional_metric_value(&summary, "unresolved_count")?.unwrap_or_else(|| metric("unresolved_count"));
            let method = optional_metric_value(&summary, "selection_method")?.unwrap_or_else(|| "pending".into());
            let selected = optional_metric_value(&summary, "selected_fasta")?;
            let mut selection = if let Some(path) = selected {
                format!("[sequence](<{}>); [{}](<{}>)", path, method, summary.display())
            } else { "pending".into() };
            if fs::symlink_metadata(component.join("manual_review.tsv")).is_ok() {
                let label = commands::postrepair::manual_review_label(&component)
                    .unwrap_or("manual decision needs verification");
                write!(selection, "; [{}](<{}>)", label, component.join("manual_review.tsv").display()).unwrap();
            }
            writeln!(out, "| {} | {} | {} | {} | {} | {} | {} |", case.name, status,
                candidates, metric("indel_candidate_count"), metric("adopted_edit_count"), unresolved, selection).unwrap();
        }
        writeln!(out).unwrap();
    }
    if config.refine.as_ref().is_some_and(|refine| refine.enabled) {
        writeln!(out, "## Repeat refinement\n\n| organelle | status | retained result |\n| --- | --- | --- |").unwrap();
        let organelles: BTreeSet<_> = config.cases.iter().filter(|case| case.enabled).map(|case| case.organelle.as_str()).collect();
        for organelle in organelles {
            let complete = refine_out_dir(config, organelle).join("auto_complete.tsv");
            let status = optional_metric_value(&complete, "status")?.unwrap_or_else(|| "pending".into());
            let retained = optional_metric_value(&complete, "retained_bundle")?;
            let retained = retained.map(|path| format!("[components](<{path}>)")).unwrap_or_else(|| "not recorded".into());
            writeln!(out, "| {} | {} | {} |", organelle, status, retained).unwrap();
        }
        writeln!(out).unwrap();
    }
    Ok(())
}

fn write_runtime_homorepair_statuses(out: &mut String, config: &WorkflowConfig) -> Result<(), OrgraftError> {
    if !config.polish.homorepair { return Ok(()); }
    writeln!(out, "## Homorepair\n").unwrap();
    writeln!(out, "| case | status | repeat support | final validation | source round | selected FASTA | outcome |\n| --- | --- | --- | --- | --- | --- | --- |").unwrap();
    for case in &config.cases {
        let path = homorepair_existing_outcome_path(case);
        let values = read_homorepair_summary_outcome(&path)?;
        let status = values.get("status").cloned().unwrap_or_else(|| "pending".into());
        let source_round = values.get("source_round").cloned().unwrap_or_else(|| ".".into());
        let selected = values.get("selected_fasta").cloned().unwrap_or_else(|| ".".into());
        let support = values.get("support_status").cloned().unwrap_or_else(|| ".".into());
        let final_validation = values.get("final_validation_status").cloned().unwrap_or_else(|| {
            if path.exists() { "not_recorded".into() } else { "pending".into() }
        });
        writeln!(out, "| {} | {} | {} | {} | {} | `{}` | `{}` |", case.name, status, support, final_validation, source_round, selected, path.display()).unwrap();
    }
    writeln!(out).unwrap();
    Ok(())
}

fn write_runtime_stage_table(
    out: &mut String,
    config: &WorkflowConfig,
) -> Result<(), OrgraftError> {
    let rows = runtime_stage_rows(config)?;
    let logged_total = rows.iter().filter_map(|row| row.seconds).sum::<f64>();

    writeln!(out, "## Timings").unwrap();
    writeln!(out, "| scope | stage | seconds | source |").unwrap();
    writeln!(out, "| --- | --- | ---: | --- |").unwrap();
    for row in &rows {
        writeln!(
            out,
            "| {} | {} | {} | `{}` |",
            row.scope,
            row.stage,
            format_optional_seconds(row.seconds),
            row.source
        )
        .unwrap();
    }
    writeln!(
        out,
        "| workflow | parsed logged total | {:.3} | parsed rows with numeric seconds |",
        logged_total
    )
    .unwrap();
    writeln!(out).unwrap();
    Ok(())
}

fn runtime_stage_rows(config: &WorkflowConfig) -> Result<Vec<RuntimeStageRow>, OrgraftError> {
    let mut rows = Vec::new();
    let recruit_log = config
        .recruit
        .out_dir
        .join("logs")
        .join("minimap2.stderr.log");
    rows.push(RuntimeStageRow {
        scope: "recruit".to_string(),
        stage: "latest minimap2 alignment".to_string(),
        seconds: read_minimap2_real_time_seconds(&recruit_log)?,
        source: recruit_log.display().to_string(),
    });

    for case in &config.cases {
        if case.postrepair_source_dir.is_none() {
        let asm_run = config
            .asm
            .out_dir
            .join(&case.organelle)
            .join("logs")
            .join("run.tsv");
        rows.push(RuntimeStageRow {
            scope: case.name.clone(),
            stage: "asm".to_string(),
            seconds: read_asm_elapsed_seconds(&asm_run)?,
            source: asm_run.display().to_string(),
        });

        let resolve_details = case
            .resolve_out_dir
            .join(&case.organelle)
            .join("logs")
            .join("resolve_details.tsv");
        rows.push(RuntimeStageRow {
            scope: case.name.clone(),
            stage: "resolve".to_string(),
            seconds: optional_section_metric_f64(&resolve_details, "run", "elapsed_seconds")?,
            source: resolve_details.display().to_string(),
        });

        }
        let homorepair_outcome = homorepair_existing_outcome_path(case);
        let current_hr = homorepair_outcome.starts_with(homorepair_out_dir(case));
        let hr_values = if postrepair_enabled(config) {
            read_homorepair_summary_outcome(&homorepair_outcome)?
        } else { BTreeMap::new() };
        let hr_validation_round = if postrepair_enabled(config) {
            let field = if current_hr { "final_validation_round" } else { "validation_round" };
            hr_values.get(field).and_then(|s| s.parse::<usize>().ok())
        } else { None };
        let hr_source_round = hr_values.get("source_round").and_then(|s| s.parse::<usize>().ok());
        for round in 1..=WORKFLOW_HARD_MAX_ROUNDS {
            // Homorepair's elapsed time includes its final ordinary validation.
            // Older layouts also included the numbered temporary probe rounds.
            let included = if current_hr { hr_validation_round == Some(round) } else {
                hr_source_round.zip(hr_validation_round).is_some_and(|(first,last)| round > first && round <= last)
            };
            if included { continue; }
            let polish_report = polish_report_path(case, round);
            if case.postrepair_source_dir.is_none() && polish_report.exists() {
                rows.push(RuntimeStageRow {
                    scope: case.name.clone(),
                    stage: format!("polish round {round} stage sum"),
                    seconds: sum_polish_stage_seconds(&polish_report)?,
                    source: polish_report.display().to_string(),
                });
            }
        }

        if postrepair_enabled(config) {
            rows.push(RuntimeStageRow {
                scope: case.name.clone(),
                stage: "polish homorepair".into(),
                seconds: hr_values.get("elapsed_seconds").and_then(|s| s.parse().ok()),
                source: homorepair_outcome.display().to_string(),
            });
        }
        let rebuild_report = rebuild_report_path(case);
        rows.push(RuntimeStageRow {
            scope: case.name.clone(),
            stage: "rebuild".to_string(),
            seconds: optional_section_metric_f64(&rebuild_report, "run", "runtime_seconds")?,
            source: rebuild_report.display().to_string(),
        });
    }

    Ok(rows)
}

fn write_runtime_finalize_table(
    out: &mut String,
    config: &WorkflowConfig,
) -> Result<(), OrgraftError> {
    writeln!(out, "## Finalize Link Cleanup").unwrap();
    writeln!(
        out,
        "| organelle | dedup_rc_links | input_links | output_links | removed | manifest |"
    )
    .unwrap();
    writeln!(out, "| --- | --- | ---: | ---: | ---: | --- |").unwrap();
    let mut seen_organelles = BTreeSet::new();
    for case in &config.cases {
        if !seen_organelles.insert(case.organelle.clone()) {
            continue;
        }
        let manifest = config
            .asm
            .out_dir
            .join(&case.organelle)
            .join("03.finalize_graph")
            .join("manifest.tsv");
        let dedup = optional_metric_value(&manifest, "finalize_dedup_rc_links")?
            .unwrap_or_else(|| "pending".to_string());
        let input = optional_metric_value(&manifest, "input_links")?.unwrap_or_default();
        let output = optional_metric_value(&manifest, "output_links")?.unwrap_or_default();
        let removed =
            optional_metric_value(&manifest, "rc_duplicate_links_removed")?.unwrap_or_default();
        writeln!(
            out,
            "| {} | {} | {} | {} | {} | `{}` |",
            case.organelle,
            dedup,
            table_value_or_pending(&input),
            table_value_or_pending(&output),
            table_value_or_pending(&removed),
            manifest.display()
        )
        .unwrap();
    }
    writeln!(out).unwrap();
    Ok(())
}

fn write_runtime_fake_validate(
    out: &mut String,
    config: &WorkflowConfig,
) -> Result<(), OrgraftError> {
    let summary = config.results_dir.join("fake_validate").join("summary.tsv");
    writeln!(out, "## Fake Validate").unwrap();
    if !summary.exists() {
        writeln!(out, "`fake_validate` output is pending.").unwrap();
        writeln!(out).unwrap();
        return Ok(());
    }

    let round_2_status_path = optional_metric_value(&summary, "round_2_status")?;
    let round_2_status = round_2_status_path
        .as_deref()
        .map(Path::new)
        .map(|path| optional_metric_value(path, "status"))
        .transpose()?
        .flatten()
        .unwrap_or_else(|| "pending".to_string());
    let round_1_status = optional_metric_value(
        &config
            .results_dir
            .join("fake_validate")
            .join("workflow/mito/subgraph_001/checkpoint_2/round_1/checkpoint_2.status.tsv"),
        "status",
    )?
    .unwrap_or_else(|| "pending".to_string());

    writeln!(out, "| step | status | source |").unwrap();
    writeln!(out, "| --- | --- | --- |").unwrap();
    writeln!(
        out,
        "| round 1 checkpoint2 | {} | `{}` |",
        round_1_status,
        summary.display()
    )
    .unwrap();
    writeln!(
        out,
        "| round 2 checkpoint2 | {} | `{}` |",
        round_2_status,
        summary.display()
    )
    .unwrap();
    writeln!(out).unwrap();
    Ok(())
}

fn latest_checkpoint2_status(
    _config: &WorkflowConfig,
    case: &WorkflowCase,
) -> Result<String, OrgraftError> {
    let mut statuses = Vec::new();
    for round in 1..=WORKFLOW_HARD_MAX_ROUNDS {
        let path = case
            .workflow_dir
            .join("checkpoint_2")
            .join(format!("round_{round}"))
            .join("checkpoint_2.status.tsv");
        if let Some(status) = optional_metric_value(&path, "status")? {
            statuses.push(format!("round_{round}:{status}"));
        }
    }
    if statuses.is_empty() {
        Ok("pending".to_string())
    } else {
        Ok(statuses.join(", "))
    }
}

fn rebuild_status(case: &WorkflowCase) -> Result<Option<String>, OrgraftError> {
    optional_section_metric_string(&rebuild_report_path(case), "run", "status")
}

fn rebuild_report_path(case: &WorkflowCase) -> PathBuf {
    case.rebuild_out_dir
        .join("logs")
        .join(format!("rebuild_{}_run_report.tsv", case.subgraph))
}

fn optional_metric_value(path: &Path, metric: &str) -> Result<Option<String>, OrgraftError> {
    if !path.exists() {
        return Ok(None);
    }
    let file = File::open(path)?;
    for (line_index, line_result) in BufReader::new(file).lines().enumerate() {
        let line = line_result?;
        if line_index == 0 && line.starts_with("metric\t") {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() >= 2 && fields[0] == metric {
            return Ok(Some(fields[1].to_string()));
        }
    }
    Ok(None)
}

fn optional_section_metric_string(
    path: &Path,
    section: &str,
    metric: &str,
) -> Result<Option<String>, OrgraftError> {
    if !path.exists() {
        return Ok(None);
    }
    let file = File::open(path)?;
    for line_result in BufReader::new(file).lines() {
        let line = line_result?;
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() >= 3 && fields[0] == section && fields[1] == metric {
            return Ok(Some(fields[2].to_string()));
        }
        if fields.len() >= 4 && fields[0] == section && fields[2] == metric {
            return Ok(Some(fields[3].to_string()));
        }
    }
    Ok(None)
}

fn optional_section_metric_f64(
    path: &Path,
    section: &str,
    metric: &str,
) -> Result<Option<f64>, OrgraftError> {
    optional_section_metric_string(path, section, metric)?
        .map(|value| parse_runtime_f64(&value, path, metric))
        .transpose()
}

fn sum_polish_stage_seconds(path: &Path) -> Result<Option<f64>, OrgraftError> {
    if !path.exists() {
        return Ok(None);
    }
    let mut sum = 0.0;
    let mut found = false;
    let file = File::open(path)?;
    for line_result in BufReader::new(file).lines() {
        let line = line_result?;
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() >= 4 && fields[0] == "stage" && fields[2] == "elapsed_seconds" {
            sum += parse_runtime_f64(fields[3], path, "elapsed_seconds")?;
            found = true;
        }
    }
    Ok(found.then_some(sum))
}

fn read_asm_elapsed_seconds(path: &Path) -> Result<Option<f64>, OrgraftError> {
    if !path.exists() {
        return Ok(None);
    }
    let file = File::open(path)?;
    let mut lines = BufReader::new(file).lines();
    let Some(header) = lines.next().transpose()? else {
        return Ok(None);
    };
    let headers: Vec<&str> = header.split('\t').collect();
    let Some(index) = headers
        .iter()
        .position(|header| *header == "elapsed_seconds")
    else {
        return Ok(None);
    };
    for line_result in lines {
        let line = line_result?;
        let fields: Vec<&str> = line.split('\t').collect();
        if let Some(value) = fields.get(index) {
            return Ok(Some(parse_runtime_f64(value, path, "elapsed_seconds")?));
        }
    }
    Ok(None)
}

fn read_minimap2_real_time_seconds(path: &Path) -> Result<Option<f64>, OrgraftError> {
    if !path.exists() {
        return Ok(None);
    }
    let file = File::open(path)?;
    for line_result in BufReader::new(file).lines() {
        let line = line_result?;
        if let Some((_, rest)) = line.split_once("Real time:") {
            let value = rest.trim().split_whitespace().next().unwrap_or_default();
            return Ok(Some(parse_runtime_f64(value, path, "Real time")?));
        }
    }
    Ok(None)
}

fn parse_runtime_f64(value: &str, path: &Path, metric: &str) -> Result<f64, OrgraftError> {
    value.parse::<f64>().map_err(|_| {
        OrgraftError::InvalidArgument(format!(
            "{} has invalid numeric value for {metric}: `{value}`",
            path.display()
        ))
    })
}

fn format_optional_seconds(seconds: Option<f64>) -> String {
    seconds
        .map(|seconds| format!("{seconds:.3}"))
        .unwrap_or_else(|| "pending".to_string())
}

fn table_value_or_pending(value: &str) -> &str {
    if value.is_empty() {
        "pending"
    } else {
        value
    }
}

fn write_checkpoint1_status(
    path: &Path,
    status: &str,
    message: &str,
    case: &WorkflowCase,
    report: &TopologyReport,
    missing_segments: &BTreeSet<String>,
    complex_nodes: &[String],
) -> Result<(), OrgraftError> {
    let mut out = String::from("metric\tvalue\n");
    push_metric(&mut out, "status", status);
    push_metric(&mut out, "message", message);
    push_metric(
        &mut out,
        "draft_graph",
        &case.draft_graph.display().to_string(),
    );
    push_metric(
        &mut out,
        "checked_draft_gfa",
        &case.checked_draft_gfa.display().to_string(),
    );
    push_metric(&mut out, "node_count", &report.node_count.to_string());
    push_metric(&mut out, "link_count", &report.link_count.to_string());
    push_metric(
        &mut out,
        "missing_link_segment_count",
        &missing_segments.len().to_string(),
    );
    push_metric(
        &mut out,
        "missing_link_segments",
        &missing_segments
            .iter()
            .cloned()
            .collect::<Vec<_>>()
            .join(","),
    );
    push_metric(
        &mut out,
        "complex_node_count",
        &complex_nodes.len().to_string(),
    );
    push_metric(&mut out, "complex_nodes", &complex_nodes.join(","));
    write_with_parent(path, &out)
}

fn write_checkpoint2_status(
    path: &Path,
    status: &str,
    message: &str,
    round: usize,
    summary_path: &Path,
    high_path: &Path,
    pos_ref_alt: Option<&Path>,
    corrected_fasta: Option<&Path>,
) -> Result<(), OrgraftError> {
    let mut out = String::from("metric\tvalue\n");
    push_metric(&mut out, "status", status);
    push_metric(&mut out, "message", message);
    push_metric(&mut out, "round", &round.to_string());
    push_metric(&mut out, "sv_summary", &summary_path.display().to_string());
    push_metric(&mut out, "snv_indel_high", &high_path.display().to_string());
    if let Some(path) = pos_ref_alt {
        push_metric(&mut out, "pos_ref_alt", &path.display().to_string());
    }
    if let Some(path) = corrected_fasta {
        push_metric(&mut out, "corrected_fasta", &path.display().to_string());
    }
    write_with_parent(path, &out)
}

fn append_checkpoint2_metrics(path: &Path, metrics: &[(&str, String)]) -> Result<(), OrgraftError> {
    let mut output = fs::read_to_string(path)?;
    for (key, value) in metrics {
        push_metric(&mut output, key, value);
    }
    fs::write(path, output)?;
    Ok(())
}

fn repeat_pairing_marker_path(case: &WorkflowCase) -> PathBuf {
    case.workflow_dir
        .join("checkpoint_2")
        .join("repeat_pairing.status.tsv")
}

fn write_repeat_pairing_marker(
    case: &WorkflowCase,
    marker: &RepeatPairingMarker,
) -> Result<(), OrgraftError> {
    let mut output = String::from("metric\tvalue\n");
    push_metric(&mut output, "state", &marker.state);
    push_metric(&mut output, "repeat_node", &marker.repeat_node);
    push_metric(
        &mut output,
        "source_round",
        &marker.source_round.to_string(),
    );
    push_metric(
        &mut output,
        "mismatch_report",
        &marker.mismatch_report.display().to_string(),
    );
    push_metric(
        &mut output,
        "rebuild_fasta",
        &marker
            .rebuild_fasta
            .as_ref()
            .map(|path| path.display().to_string())
            .unwrap_or_else(|| ".".to_string()),
    );
    push_metric(
        &mut output,
        "validation_round",
        &marker
            .validation_round
            .map(|round| round.to_string())
            .unwrap_or_else(|| ".".to_string()),
    );
    write_with_parent(&repeat_pairing_marker_path(case), &output)
}

fn read_repeat_pairing_marker(
    case: &WorkflowCase,
) -> Result<Option<RepeatPairingMarker>, OrgraftError> {
    let path = repeat_pairing_marker_path(case);
    if !path.is_file() {
        return Ok(None);
    }
    let required = |metric: &str| -> Result<String, OrgraftError> {
        read_metric_value(&path, metric)?.ok_or_else(|| {
            OrgraftError::InvalidArgument(format!("{} is missing `{metric}`", path.display()))
        })
    };
    let optional_path = |metric: &str| -> Result<Option<PathBuf>, OrgraftError> {
        Ok(match required(metric)?.as_str() {
            "." | "" => None,
            value => Some(PathBuf::from(value)),
        })
    };
    let optional_round = |metric: &str| -> Result<Option<usize>, OrgraftError> {
        Ok(match required(metric)?.as_str() {
            "." | "" => None,
            value => Some(value.parse::<usize>().map_err(|_| {
                OrgraftError::InvalidArgument(format!(
                    "{} has invalid `{metric}` value `{value}`",
                    path.display()
                ))
            })?),
        })
    };
    Ok(Some(RepeatPairingMarker {
        state: required("state")?,
        repeat_node: required("repeat_node")?,
        source_round: required("source_round")?.parse::<usize>().map_err(|_| {
            OrgraftError::InvalidArgument(format!("{} has invalid source_round", path.display()))
        })?,
        mismatch_report: PathBuf::from(required("mismatch_report")?),
        rebuild_fasta: optional_path("rebuild_fasta")?,
        validation_round: optional_round("validation_round")?,
    }))
}

fn append_repeat_pairing_metrics(
    status_path: &Path,
    mismatch: &RepeatPairingMismatch,
    mismatch_report: &Path,
    mode: &str,
) -> Result<(), OrgraftError> {
    append_checkpoint2_metrics(
        status_path,
        &[
            ("sv_correction_mode", mode.to_string()),
            ("repeat_node", mismatch.repeat_node.clone()),
            ("repeat_length", mismatch.repeat_length.to_string()),
            ("sv_subgroups", mismatch.subgroup_specs.join(",")),
            (
                "alternate_path_reads",
                mismatch
                    .alternate_path_reads
                    .iter()
                    .map(ToString::to_string)
                    .collect::<Vec<_>>()
                    .join(","),
            ),
            (
                "reference_spanning_reads",
                mismatch
                    .reference_spanning_reads
                    .iter()
                    .map(ToString::to_string)
                    .collect::<Vec<_>>()
                    .join(","),
            ),
            (
                "repeat_pairing_dominance_ratio",
                format!("{:.6}", mismatch.dominance_ratio),
            ),
            (
                "repeat_pairing_mismatch_report",
                mismatch_report.display().to_string(),
            ),
        ],
    )
}

fn append_repeat_marker_metrics(
    status_path: &Path,
    marker: &RepeatPairingMarker,
    mode: &str,
) -> Result<(), OrgraftError> {
    append_checkpoint2_metrics(
        status_path,
        &[
            ("sv_correction_mode", mode.to_string()),
            ("repeat_node", marker.repeat_node.clone()),
            (
                "repeat_pairing_marker",
                repeat_pairing_marker_path_for_display(status_path, marker),
            ),
            (
                "repeat_pairing_mismatch_report",
                marker.mismatch_report.display().to_string(),
            ),
            (
                "repeat_pairing_source_round",
                marker.source_round.to_string(),
            ),
            (
                "repeat_pairing_validation_round",
                marker
                    .validation_round
                    .map(|round| round.to_string())
                    .unwrap_or_else(|| ".".to_string()),
            ),
        ],
    )
}

fn repeat_pairing_marker_path_for_display(
    status_path: &Path,
    marker: &RepeatPairingMarker,
) -> String {
    status_path
        .parent()
        .and_then(Path::parent)
        .map(|path| path.join("repeat_pairing.status.tsv"))
        .unwrap_or_else(|| marker.mismatch_report.clone())
        .display()
        .to_string()
}

fn rebuild_fasta_path(case: &WorkflowCase) -> PathBuf {
    case.rebuild_out_dir
        .join(&case.subgraph)
        .join(format!("rebuild_{}.fasta", case.subgraph))
}

fn rebuild_result_stats_path(case: &WorkflowCase) -> PathBuf {
    case.rebuild_out_dir
        .join("logs")
        .join(format!("rebuild_{}_result_stats.tsv", case.subgraph))
}

fn rebuild_has_repeat_resolution(
    case: &WorkflowCase,
    repeat_node: &str,
) -> Result<bool, OrgraftError> {
    let path = rebuild_result_stats_path(case);
    if !path.is_file() {
        return Ok(false);
    }
    for line_result in BufReader::new(File::open(path)?).lines() {
        let line = line_result?;
        let fields = line.split('\t').collect::<Vec<_>>();
        if fields.len() >= 5 && fields[0] == "summary" && fields[3] == "repeat_resolution_node" {
            return Ok(fields[4] == repeat_node);
        }
    }
    Ok(false)
}

fn files_equal(left: &Path, right: &Path) -> Result<bool, OrgraftError> {
    if fs::metadata(left)?.len() != fs::metadata(right)?.len() {
        return Ok(false);
    }
    Ok(fs::read(left)? == fs::read(right)?)
}

fn remove_path_if_exists(path: &Path) -> Result<(), OrgraftError> {
    if path.exists() {
        fs::remove_file(path)?;
    }
    Ok(())
}

fn sv_correction_count(case: &WorkflowCase, before_round: usize) -> Result<usize, OrgraftError> {
    let mut count = 0usize;
    for round in 1..before_round {
        let status_path = case
            .workflow_dir
            .join("checkpoint_2")
            .join(format!("round_{round}"))
            .join("checkpoint_2.status.tsv");
        if read_metric_value(&status_path, "correction_kind")?.as_deref() == Some("sv") {
            count += 1;
        }
    }
    Ok(count)
}

fn sv_subgroup_already_corrected(
    case: &WorkflowCase,
    before_round: usize,
    requested: &str,
) -> Result<bool, OrgraftError> {
    let requested = canonical_sv_subgroup_spec(requested);
    for round in 1..before_round {
        let status_path = case
            .workflow_dir
            .join("checkpoint_2")
            .join(format!("round_{round}"))
            .join("checkpoint_2.status.tsv");
        if read_metric_value(&status_path, "correction_kind")?.as_deref() != Some("sv") {
            continue;
        }
        if read_metric_value(&status_path, "sv_subgroup")?
            .as_deref()
            .is_some_and(|value| canonical_sv_subgroup_spec(value) == requested)
        {
            return Ok(true);
        }
    }
    Ok(false)
}

fn canonical_sv_subgroup_spec(value: &str) -> String {
    let Some((group, index)) = value.rsplit_once(':') else {
        return value.to_string();
    };
    format!(
        "{}:{}",
        group.replace("_subtype_", "_").trim_end_matches("_NA"),
        index
    )
}

fn checkpoint_history_fastas(case: &WorkflowCase, round: usize) -> Vec<PathBuf> {
    let mut paths = Vec::new();
    for previous_round in 1..round {
        paths.push(polished_aln_path(case, previous_round));
        paths.push(corrected_fasta_path(case, previous_round));
    }
    paths.push(polished_aln_path(case, round));
    paths
}

fn read_metric_value(path: &Path, metric: &str) -> Result<Option<String>, OrgraftError> {
    let file = File::open(path).map_err(|error| {
        OrgraftError::InvalidArgument(format!(
            "could not read {}; run polish for this round first ({error})",
            path.display()
        ))
    })?;
    for (line_index, line_result) in BufReader::new(file).lines().enumerate() {
        let line = line_result?;
        if line_index == 0 && line.starts_with("metric\t") {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() >= 2 && fields[0] == metric {
            return Ok(Some(fields[1].to_string()));
        }
    }
    Ok(None)
}

fn write_fake_sv_pass_summary(path: &Path) -> Result<(), OrgraftError> {
    write_with_parent(path, "metric\tvalue\nstatus\tpass\n")
}

fn write_fake_snv_high_table(path: &Path, edits: &[VariantEdit]) -> Result<(), OrgraftError> {
    let mut out = String::from("pos\tref\talt\tfixed_ref\n");
    for edit in edits {
        writeln!(
            out,
            "{}\t{}\t{}\tNo",
            edit.pos, edit.reference, edit.alternate
        )
        .unwrap();
    }
    write_with_parent(path, &out)
}

fn write_fake_validate_summary(
    path: &Path,
    input: &str,
    variants: &str,
    case: &WorkflowCase,
) -> Result<(), OrgraftError> {
    let mut out = String::from("metric\tvalue\n");
    push_metric(&mut out, "input_fasta", input);
    push_metric(&mut out, "pos_ref_alt", variants);
    push_metric(
        &mut out,
        "round_1_swapped_polish_aln",
        &polished_aln_path(case, 1).display().to_string(),
    );
    push_metric(
        &mut out,
        "round_1_corrected_fasta",
        &corrected_fasta_path(case, 1).display().to_string(),
    );
    push_metric(
        &mut out,
        "round_2_polish_aln",
        &polished_aln_path(case, 2).display().to_string(),
    );
    push_metric(
        &mut out,
        "round_2_status",
        &case
            .workflow_dir
            .join("checkpoint_2/round_2/checkpoint_2.status.tsv")
            .display()
            .to_string(),
    );
    write_with_parent(path, &out)
}

// Use checkpoint2's ordinary correction rule, rather than counting every high
// variant row: a high-frequency row can still favor REF over every single ALT.
pub(crate) fn ordinary_validation_edit_count(path: &Path) -> Result<usize, OrgraftError> {
    Ok(read_high_variant_edits(path)?.len())
}

fn read_high_variant_edits(path: &Path) -> Result<Vec<VariantEdit>, OrgraftError> {
    let file = File::open(path).map_err(|error| {
        OrgraftError::InvalidArgument(format!(
            "could not read {}; run polish for this round first ({error})",
            path.display()
        ))
    })?;
    let mut lines = BufReader::new(file).lines();
    let Some(header) = lines.next() else {
        return Ok(Vec::new());
    };
    let header = header?;
    let columns: Vec<String> = header.split('\t').map(ToString::to_string).collect();
    let index = |name: &str| columns.iter().position(|column| column == name);
    let pos_index = index("pos").ok_or_else(|| {
        OrgraftError::InvalidArgument(format!("{} is missing `pos` column", path.display()))
    })?;
    let ref_index = index("ref").ok_or_else(|| {
        OrgraftError::InvalidArgument(format!("{} is missing `ref` column", path.display()))
    })?;
    let alt_index = index("alt").ok_or_else(|| {
        OrgraftError::InvalidArgument(format!("{} is missing `alt` column", path.display()))
    })?;
    let counts_index = index("counts");
    let total_count_index = index("total_count");
    let depth_index = index("depth");

    let mut edits = Vec::new();
    for line_result in lines {
        let line = line_result?;
        if line.trim().is_empty() {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        if reference_support_is_top(&fields, counts_index, total_count_index, depth_index) {
            continue;
        }
        let pos = fields
            .get(pos_index)
            .ok_or_else(|| {
                OrgraftError::InvalidArgument(format!(
                    "{} contains a row without pos",
                    path.display()
                ))
            })?
            .parse::<usize>()
            .map_err(|_| {
                OrgraftError::InvalidArgument(format!(
                    "{} contains a non-integer pos",
                    path.display()
                ))
            })?;
        let reference = fields.get(ref_index).copied().unwrap_or("").to_string();
        let alternate = choose_alternate(
            fields.get(alt_index).copied().unwrap_or(""),
            counts_index.and_then(|idx| fields.get(idx).copied()),
        );
        if !alternate.is_empty() {
            edits.push(VariantEdit {
                pos,
                reference,
                alternate,
            });
        }
    }
    edits.sort_by_key(|edit| edit.pos);
    Ok(edits)
}

fn reference_support_is_top(
    fields: &[&str],
    counts_index: Option<usize>,
    total_count_index: Option<usize>,
    depth_index: Option<usize>,
) -> bool {
    let Some(counts) = counts_index
        .and_then(|idx| fields.get(idx))
        .map(|value| parse_count_list(value))
    else {
        return false;
    };
    let Some(total_count) =
        total_count_index.and_then(|idx| fields.get(idx)?.parse::<usize>().ok())
    else {
        return false;
    };
    let Some(depth) = depth_index.and_then(|idx| fields.get(idx)?.parse::<usize>().ok()) else {
        return false;
    };
    let Some(max_alt_count) = counts.iter().copied().max() else {
        return false;
    };
    depth.saturating_sub(total_count) >= max_alt_count
}

fn choose_alternate(alt: &str, counts: Option<&str>) -> String {
    let alts: Vec<&str> = alt.split('#').collect();
    if alts.len() <= 1 {
        return alt.to_string();
    }
    let Some(counts) = counts else {
        return alts[0].to_string();
    };
    let parsed_counts = parse_count_list(counts);
    if parsed_counts.len() != alts.len() {
        return alts[0].to_string();
    }
    let best_index = parsed_counts
        .iter()
        .enumerate()
        .max_by_key(|(_, count)| *count)
        .map(|(index, _)| index)
        .unwrap_or(0);
    alts[best_index].to_string()
}

fn parse_count_list(counts: &str) -> Vec<usize> {
    counts
        .split('#')
        .map(|value| value.parse::<usize>().unwrap_or(0))
        .collect()
}

fn read_pos_ref_alt(path: &Path) -> Result<Vec<VariantEdit>, OrgraftError> {
    let file = File::open(path)?;
    let mut edits = Vec::new();
    for (line_index, line_result) in BufReader::new(file).lines().enumerate() {
        let line_number = line_index + 1;
        let line = line_result?;
        if line.trim().is_empty() || line.starts_with('#') {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 3 {
            return Err(OrgraftError::InvalidArgument(format!(
                "{}:{line_number}: expected pos<TAB>ref<TAB>alt",
                path.display()
            )));
        }
        let pos = fields[0].parse::<usize>().map_err(|_| {
            OrgraftError::InvalidArgument(format!(
                "{}:{line_number}: pos must be a positive integer",
                path.display()
            ))
        })?;
        if pos == 0 {
            return Err(OrgraftError::InvalidArgument(format!(
                "{}:{line_number}: pos must be 1-based and greater than 0",
                path.display()
            )));
        }
        edits.push(VariantEdit {
            pos,
            reference: fields[1].to_string(),
            alternate: fields[2].to_string(),
        });
    }
    edits.sort_by_key(|edit| edit.pos);
    Ok(edits)
}

fn write_pos_ref_alt(path: &Path, edits: &[VariantEdit]) -> Result<(), OrgraftError> {
    let mut out = String::new();
    for edit in edits {
        writeln!(out, "{}\t{}\t{}", edit.pos, edit.reference, edit.alternate).unwrap();
    }
    write_with_parent(path, &out)
}

fn apply_edits_to_fasta(
    input: &Path,
    output: &Path,
    edits: &[VariantEdit],
) -> Result<String, OrgraftError> {
    let (header, sequence) = read_single_fasta(input)?;
    let corrected = apply_edits_to_sequence(&sequence, edits)?;
    write_fasta(output, &header, &corrected)?;
    Ok(format!(
        "old_length={} new_length={} edit_count={}",
        sequence.len(),
        corrected.len(),
        edits.len()
    ))
}

fn correction_summary_tsv(summary: &str) -> String {
    summary
        .split_whitespace()
        .filter_map(|entry| entry.split_once('='))
        .map(|(key, value)| format!("{key}\t{value}"))
        .collect::<Vec<_>>()
        .join("\n")
}

fn read_single_fasta(path: &Path) -> Result<(String, String), OrgraftError> {
    let file = File::open(path)?;
    let mut header = None;
    let mut sequence = String::new();
    for line_result in BufReader::new(file).lines() {
        let line = line_result?;
        if let Some(rest) = line.strip_prefix('>') {
            if header.is_some() {
                return Err(OrgraftError::InvalidArgument(format!(
                    "{} contains more than one FASTA record",
                    path.display()
                )));
            }
            header = Some(rest.to_string());
        } else {
            sequence.push_str(line.trim());
        }
    }
    let header = header.ok_or_else(|| {
        OrgraftError::InvalidArgument(format!("{} contains no FASTA record", path.display()))
    })?;
    Ok((header, sequence))
}

fn apply_edits_to_sequence(sequence: &str, edits: &[VariantEdit]) -> Result<String, OrgraftError> {
    if edits.is_empty() {
        return Ok(sequence.to_string());
    }

    let mut sorted = edits.to_vec();
    sorted.sort_by_key(|edit| edit.pos);
    let mut groups: Vec<Vec<VariantEdit>> = Vec::new();
    let mut current: Vec<VariantEdit> = Vec::new();
    let mut current_end = 0usize;
    for edit in sorted {
        let start = edit.pos.checked_sub(1).ok_or_else(|| {
            OrgraftError::InvalidArgument("variant positions must be 1-based".to_string())
        })?;
        let end = start + edit.reference.len();
        if current.is_empty() || start < current_end {
            current_end = current_end.max(end);
            current.push(edit);
        } else {
            groups.push(current);
            current = vec![edit];
            current_end = end;
        }
    }
    if !current.is_empty() {
        groups.push(current);
    }

    let mut output = String::new();
    let mut cursor = 0usize;
    for group in groups {
        let start = group[0].pos - 1;
        let end = group
            .iter()
            .map(|edit| edit.pos - 1 + edit.reference.len())
            .max()
            .unwrap_or(start);
        if end > sequence.len() {
            return Err(OrgraftError::InvalidArgument(format!(
                "variant at position {} extends past FASTA length {}",
                group[0].pos,
                sequence.len()
            )));
        }
        if start > sequence.len() {
            return Err(OrgraftError::InvalidArgument(format!(
                "variant at position {} is past FASTA length {}",
                group[0].pos,
                sequence.len()
            )));
        }
        output.push_str(&sequence[cursor..start]);
        output.push_str(&replacement_for_group(sequence, &group)?);
        cursor = end;
    }
    output.push_str(&sequence[cursor..]);
    Ok(output)
}

fn replacement_for_group(sequence: &str, group: &[VariantEdit]) -> Result<String, OrgraftError> {
    if group.len() == 1 {
        return Ok(group[0].alternate.clone());
    }

    let mut replacement = String::new();
    for (index, edit) in group.iter().enumerate() {
        let start = edit.pos - 1;
        let mut ref_len = edit.reference.len();
        let mut alt = edit.alternate.clone();
        if let Some(next) = group.get(index + 1) {
            let next_start = next.pos - 1;
            let delta = start + ref_len;
            if delta > next_start {
                let overlap = delta - next_start;
                ref_len = ref_len.saturating_sub(overlap);
                let keep = alt.len().saturating_sub(overlap);
                alt.truncate(keep);
            }
            replacement.push_str(&alt);
            let gap_start = start + ref_len;
            if gap_start < next_start {
                replacement.push_str(sequence.get(gap_start..next_start).ok_or_else(|| {
                    OrgraftError::InvalidArgument(
                        "overlapping variant coordinates are outside FASTA bounds".to_string(),
                    )
                })?);
            }
        } else {
            replacement.push_str(&alt);
        }
    }
    Ok(replacement)
}

fn write_fasta(path: &Path, header: &str, sequence: &str) -> Result<(), OrgraftError> {
    let mut out = String::new();
    writeln!(out, ">{header}").unwrap();
    for chunk in sequence.as_bytes().chunks(80) {
        out.push_str(std::str::from_utf8(chunk).map_err(|_| {
            OrgraftError::InvalidArgument("FASTA sequence is not valid UTF-8".to_string())
        })?);
        out.push('\n');
    }
    write_with_parent(path, &out)
}

fn discover_unitig_graph(config: &WorkflowConfig, case: &WorkflowCase) -> Option<PathBuf> {
    if let Some(path) = &case.unitig_graph {
        return Some(path.clone());
    }
    let asm_path = config
        .asm
        .out_dir
        .join(&case.organelle)
        .join("02.anchor_graph_core/02.unitig_graph/graph.gfa");
    let draft_relative = case
        .draft_graph
        .parent()
        .and_then(Path::parent)
        .map(|organelle_dir| organelle_dir.join("02.anchor_graph_core/02.unitig_graph/graph.gfa"));
    if asm_path.exists() {
        Some(asm_path)
    } else if draft_relative.as_ref().is_some_and(|path| path.exists()) {
        draft_relative
    } else {
        Some(asm_path)
    }
}

fn round_draft_path(case: &WorkflowCase, round: usize) -> PathBuf {
    if round == 1 {
        case.linearized_fasta
            .clone()
            .unwrap_or_else(|| resolved_subgraphs_path(case))
    } else {
        corrected_fasta_path(case, round - 1)
    }
}

fn resolved_subgraphs_path(case: &WorkflowCase) -> PathBuf {
    case.resolve_out_dir
        .join(&case.organelle)
        .join("fasta")
        .join("resolved_subgraphs.fasta")
}

fn polish_reference_path(case: &WorkflowCase) -> PathBuf {
    if let Some(source) = &case.postrepair_source_dir { return source.join("01.inputs/rotated_reference.fasta"); }
    case.polish_reference.clone().unwrap_or_else(|| {
        case.resolve_out_dir
            .join(&case.organelle)
            .join("fasta")
            .join("rotated_reference.fasta")
    })
}

fn polish_round_dir(case: &WorkflowCase, round: usize) -> PathBuf {
    polish_run_out_dir(case, round)
        .join(&case.organelle)
        .join(&case.subgraph)
        .join(format!("round_{round}"))
}

fn polish_run_out_dir(case: &WorkflowCase, _round: usize) -> PathBuf {
    case.polish_out_dir.clone()
}

fn sv_summary_path(case: &WorkflowCase, round: usize) -> PathBuf {
    polish_round_dir(case, round).join("03.validate/03.reports/sv_snv_indel_summary.tsv")
}

fn snv_indel_high_path(case: &WorkflowCase, round: usize) -> PathBuf {
    polish_round_dir(case, round).join("03.validate/03.reports/snv_indel_high.tsv")
}

fn sv_high_subgroups_path(case: &WorkflowCase, round: usize) -> PathBuf {
    polish_round_dir(case, round).join("03.validate/03.reports/sv_high_subgroups.tsv")
}

fn sv_subgroup_stats_path(case: &WorkflowCase, round: usize) -> PathBuf {
    polish_round_dir(case, round).join("03.validate/03.reports/sv_subgroup_stats.tsv")
}

fn snv_indel_segments_path(case: &WorkflowCase, round: usize) -> PathBuf {
    polish_round_dir(case, round).join("03.validate/01.data/snv_indel_segments.tsv")
}

fn sv_read_index_path(case: &WorkflowCase, round: usize) -> PathBuf {
    polish_round_dir(case, round).join("03.validate/01.data/sv_read_index.tsv")
}

fn sv_read_evidence_path(case: &WorkflowCase, round: usize) -> PathBuf {
    polish_round_dir(case, round).join("03.validate/01.data/sv_read_evidence.tsv")
}

fn polished_aln_path(case: &WorkflowCase, round: usize) -> PathBuf {
    if round == 1 {
        polish_round_dir(case, round).join("02.polish/polished_aln.fasta")
    } else {
        polish_round_input_path(case, round)
    }
}

fn polish_round_input_path(case: &WorkflowCase, round: usize) -> PathBuf {
    polish_round_dir(case, round)
        .join("01.inputs")
        .join(format!("linear_subgraph.round_{round}.fasta"))
}

fn polish_report_path(case: &WorkflowCase, round: usize) -> PathBuf {
    polish_round_dir(case, round).join("logs/report.tsv")
}

fn corrected_fasta_path(case: &WorkflowCase, round: usize) -> PathBuf {
    case.workflow_dir
        .join("checkpoint_2")
        .join(format!("round_{round}"))
        .join(format!("polish_aln_v{}.fasta", round + 1))
}

fn resolve_command(config: &WorkflowConfig, case: &WorkflowCase, force: bool) -> String {
    let mut args = vec![
        orgraft_shell_token().to_string(),
        "resolve".to_string(),
        "--checked-draft-gfa".to_string(),
        shell_quote(&case.checked_draft_gfa),
        "--soft-paths".to_string(),
        shell_quote(&config.soft_paths),
        "--out-dir".to_string(),
        shell_quote(&case.resolve_out_dir),
        "--organelle".to_string(),
        shell_quote_str(&case.organelle),
    ];
    if let Some(pre_rotated_reference) = &case.pre_rotated_reference {
        args.push("--pre-rotated-reference".to_string());
        args.push(shell_quote(pre_rotated_reference));
    } else if let Some(reference) = &case.reference {
        args.push("--reference".to_string());
        args.push(shell_quote(reference));
    } else {
        args.push("--reference".to_string());
        args.push("FILE".to_string());
    }
    if force {
        args.push("--force".to_string());
    }
    append_resolve_shell_args(&mut args, &config.resolve);
    args.join(" ")
}

fn polish_command(
    config: &WorkflowConfig,
    case: &WorkflowCase,
    round: usize,
    draft: &Path,
    force: bool,
) -> String {
    polish_command_with_draft_arg(config, case, round, &shell_quote(draft), force)
}

fn polish_command_with_draft_arg(
    config: &WorkflowConfig,
    case: &WorkflowCase,
    round: usize,
    draft_arg: &str,
    force: bool,
) -> String {
    let mut args = vec![
        orgraft_shell_token().to_string(),
        "polish".to_string(),
        "--organelle".to_string(),
        shell_quote_str(&case.organelle),
        "--subgraph".to_string(),
        shell_quote_str(&case.subgraph),
        "--draft".to_string(),
        draft_arg.to_string(),
        "--reference".to_string(),
        shell_quote(&polish_reference_path(case)),
        "--reads".to_string(),
        shell_quote(&case.reads),
        "--soft-paths".to_string(),
        shell_quote(&config.soft_paths),
        "--out-dir".to_string(),
        shell_quote(&polish_run_out_dir(case, round)),
        "--validate-round".to_string(),
        round.to_string(),
        "--threads".to_string(),
        polish_threads(config).to_string(),
        "--max-rounds".to_string(),
        config.max_rounds.to_string(),
    ];
    if force {
        args.push("--force".to_string());
    }
    append_polish_shell_args(&mut args, &config.polish);
    args.join(" ")
}

fn check_plan_destination(path: &Path, content: &str, replace: bool) -> Result<(), OrgraftError> {
    match fs::symlink_metadata(path) {
        Ok(metadata) if metadata.is_file() && (replace || fs::read(path)? == content.as_bytes()) => Ok(()),
        Ok(_) => Err(OrgraftError::InvalidArgument(format!(
            "workflow plan differs from existing {}; review it and use plan --force to regenerate both scripts, or choose a new --out directory", path.display()
        ))),
        Err(error) if error.kind() == std::io::ErrorKind::NotFound => Ok(()),
        Err(error) => Err(error.into()),
    }
}

fn write_new_plan_script(path: &Path, content: &str, replace: bool) -> Result<(), OrgraftError> {
    use std::io::Write;
    if let Some(parent) = path.parent() { fs::create_dir_all(parent)?; }
    if replace && path.exists() {
        check_plan_destination(path, content, true)?;
        if fs::read(path)? == content.as_bytes() { return Ok(()); }
        let name = path.file_name().ok_or_else(|| OrgraftError::InvalidArgument("script path has no filename".into()))?;
        let temporary = path.with_file_name(format!(".{}.orgraft-{}", name.to_string_lossy(), std::process::id()));
        write_new_plan_script(&temporary, content, false)?;
        if let Err(error) = fs::rename(&temporary, path) {
            let _ = fs::remove_file(&temporary);
            return Err(error.into());
        }
        return Ok(());
    }
    let mut file = match fs::OpenOptions::new().write(true).create_new(true).open(path) {
        Ok(file) => file,
        Err(error) if error.kind() == std::io::ErrorKind::AlreadyExists => return check_plan_destination(path, content, false),
        Err(error) => return Err(error.into()),
    };
    file.write_all(content.as_bytes())?;
    #[cfg(unix)]
    {
        use std::os::unix::fs::PermissionsExt;
        file.set_permissions(fs::Permissions::from_mode(0o755))?;
    }
    Ok(())
}

fn write_executable_text(path: &Path, content: &str) -> Result<(), OrgraftError> {
    write_with_parent(path, content)?;
    #[cfg(unix)]
    {
        use std::os::unix::fs::PermissionsExt;
        let mut permissions = fs::metadata(path)?.permissions();
        permissions.set_mode(0o755);
        fs::set_permissions(path, permissions)?;
    }
    Ok(())
}

fn copy_with_parent(from: &Path, to: &Path, force: bool) -> Result<(), OrgraftError> {
    if to.exists() && !force {
        return Err(OrgraftError::InvalidArgument(format!(
            "{} already exists; pass --force to overwrite",
            to.display()
        )));
    }
    if let Some(parent) = to.parent() {
        fs::create_dir_all(parent)?;
    }
    fs::copy(from, to)?;
    Ok(())
}

fn write_with_parent(path: &Path, content: &str) -> Result<(), OrgraftError> {
    if let Some(parent) = path.parent() {
        fs::create_dir_all(parent)?;
    }
    fs::write(path, content)?;
    Ok(())
}

fn push_metric(out: &mut String, metric: &str, value: &str) {
    writeln!(out, "{metric}\t{}", value.replace('\n', " ")).unwrap();
}

fn option_value<'a>(args: &'a [String], name: &str) -> Result<Option<&'a str>, OrgraftError> {
    for (index, arg) in args.iter().enumerate() {
        if arg == name {
            let Some(value) = args.get(index + 1) else {
                return Err(OrgraftError::InvalidArgument(format!(
                    "missing value for {name}"
                )));
            };
            if value.starts_with("--") {
                return Err(OrgraftError::InvalidArgument(format!(
                    "missing value for {name}"
                )));
            }
            return Ok(Some(value));
        }
    }
    Ok(None)
}

fn has_flag(args: &[String], name: &str) -> bool {
    args.iter().any(|arg| arg == name)
}

fn push_raw_string_option(args: &mut Vec<String>, name: &str, value: Option<&str>) {
    if let Some(value) = value {
        args.push(name.to_string());
        args.push(value.to_string());
    }
}

fn push_raw_path_option(args: &mut Vec<String>, name: &str, value: Option<&Path>) {
    if let Some(value) = value {
        args.push(name.to_string());
        args.push(value.display().to_string());
    }
}

fn push_raw_display_option<T>(args: &mut Vec<String>, name: &str, value: Option<T>)
where
    T: std::fmt::Display,
{
    if let Some(value) = value {
        args.push(name.to_string());
        args.push(value.to_string());
    }
}

fn push_raw_flag(args: &mut Vec<String>, name: &str, enabled: bool) {
    if enabled {
        args.push(name.to_string());
    }
}

fn push_string_option(args: &mut Vec<String>, name: &str, value: Option<&str>) {
    if let Some(value) = value {
        args.push(name.to_string());
        args.push(shell_quote_str(value));
    }
}

fn push_path_option(args: &mut Vec<String>, name: &str, value: Option<&Path>) {
    if let Some(value) = value {
        args.push(name.to_string());
        args.push(shell_quote(value));
    }
}

fn push_display_option<T>(args: &mut Vec<String>, name: &str, value: Option<T>)
where
    T: std::fmt::Display,
{
    if let Some(value) = value {
        args.push(name.to_string());
        args.push(value.to_string());
    }
}

fn push_flag(args: &mut Vec<String>, name: &str, enabled: bool) {
    if enabled {
        args.push(name.to_string());
    }
}

fn shell_quote(path: &Path) -> String {
    shell_quote_str(&path.display().to_string())
}

fn default_orgraft_bin() -> String {
    std::env::current_exe()
        .ok()
        .map(|path| shell_quote(&path))
        .unwrap_or_else(|| "orgraft".to_string())
}

fn shell_quote_str(value: &str) -> String {
    if value
        .chars()
        .all(|ch| ch.is_ascii_alphanumeric() || matches!(ch, '/' | '.' | '_' | '-' | ':' | '+'))
    {
        return value.to_string();
    }
    format!("'{}'", value.replace('\'', "'\\''"))
}

fn toml_string(value: &str) -> String {
    let mut escaped = String::from("\"");
    for character in value.chars() {
        match character {
            '\\' => escaped.push_str("\\\\"),
            '"' => escaped.push_str("\\\""),
            '\n' => escaped.push_str("\\n"),
            '\r' => escaped.push_str("\\r"),
            '\t' => escaped.push_str("\\t"),
            _ => escaped.push(character),
        }
    }
    escaped.push('"');
    escaped
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::time::{SystemTime, UNIX_EPOCH};

    #[test]
    fn workflow_sample_name_comes_from_the_config_folder() {
        let options = TemplateOptions::from_args(&["--out".into(), "S42/orgraft.workflow.toml".into()]).unwrap();
        assert_eq!(options.sample, "S42");
        assert_eq!(options.results_dir, ".");
    }

    #[test]
    fn explicit_postrepair_separates_source_from_output_and_refine_waits_for_all_cases() {
        let dir = test_dir("workflow_explicit_postrepair");
        fs::create_dir_all(&dir).unwrap();
        let path = dir.join("orgraft.workflow.toml");
        fs::write(&path, "[project]\nsample = \"S42\"\nresults_dir = \".\"\n[commands.postrepair]\nenabled = true\n[commands.rebuild]\nenabled = true\nout_dir = \"${results_dir}/06.rebuild\"\n[commands.refine]\nenabled = true\n[workflow.case.first]\norganelle = \"mito\"\nsubgraph = \"subgraph_001\"\n[workflow.case.second]\norganelle = \"mito\"\nsubgraph = \"subgraph_002\"\n").unwrap();
        let config = WorkflowConfig::from_path(&path).unwrap();
        let case = &config.cases[0];
        let args = homorepair_args(&config, case, 2);
        assert!(args.windows(2).any(|p| p[0] == "--source-dir" && p[1].ends_with("04.polish/mito/subgraph_001/round_2")));
        assert!(args.windows(2).any(|p| p[0] == "--out-dir" && p[1].ends_with("05.postrepair")));
        assert!(!args.iter().any(|p| p == "--homorepair-source-round"));
        assert!(config.soft_paths.is_absolute());
        let err = prepare_refine_bundle(&config, "mito").unwrap_err().to_string();
        assert!(err.contains("refine pending"));
        assert!(!refine_bundle_path(&config, "mito").exists());
        let script = write_plan_from_args(&config, &[]).unwrap();
        let text = fs::read_to_string(&script).unwrap();
        assert!(script.ends_with("run_S42.sh"));
        assert!(dir.join("run_S42_manual.sh").exists());
        assert_eq!(text.matches(" refine --sample-bundle ").count(), 1);
        assert_eq!(text.matches(" resolve --checked-draft-gfa ").count(), 1);
        assert!(text.contains("workflow postrepair-check"));
        assert!(text.contains("Checked component inputs differ"));
        assert!(!text.contains("bash "));
        assert!(Command::new("bash").arg("-n").arg(&script).status().unwrap().success());
        assert!(Command::new("bash").arg("-n").arg(dir.join("run_S42_manual.sh")).status().unwrap().success());
    }

    // These small native reports exercise only the 06 -> 07 handoff, without
    // running mapping, repair or graph construction.
    fn refine_reference_fixture(dir: &Path) -> WorkflowConfig {
        fs::create_dir_all(dir).unwrap();
        let path = dir.join("workflow.toml");
        fs::write(&path, "[project]\nsample = \"fixture\"\nresults_dir = \".\"\n[commands.rebuild]\nenabled = true\n[commands.refine]\nenabled = true\n[workflow.case.first]\norganelle = \"mito\"\nsubgraph = \"subgraph_001\"\n[workflow.case.second]\norganelle = \"mito\"\nsubgraph = \"subgraph_002\"\n").unwrap();
        let config = WorkflowConfig::from_path(&path).unwrap();
        write_with_parent(
            &polish_reference_path(&config.cases[0]),
            ">reference_001\nACGTA\n>reference_002\nTTGCC\n",
        )
        .unwrap();
        for (index, case) in config.cases.iter().enumerate() {
            let id = &case.subgraph;
            let dna = if index == 0 { "ACGTA" } else { "TTGCC" };
            let round = polish_round_dir(case, 2);
            let validation = round.join("03.validate");
            fs::create_dir_all(validation.join("01.data")).unwrap();
            let reference = round.join("01.inputs/rotated_reference.fasta");
            write_with_parent(&reference, &format!(">reference_{:03}\n{dna}\n", index + 1)).unwrap();
            write_with_parent(&case.reads, "@read\nACGT\n+\nIIII\n").unwrap();
            fs::copy(&case.reads, round.join("01.inputs/subgraph_reads.fastq.gz")).unwrap();
            let source = rebuild_source_gfa(&config, case);
            let graph = format!("S\tn{index}\t{dna}\n");
            write_with_parent(&source, &graph).unwrap();
            write_with_parent(&case.checked_draft_gfa, &graph).unwrap();
            let details = resolve_details_path(case);
            let mut text = if details.exists() {
                fs::read_to_string(&details).unwrap()
            } else {
                "section\tid\tkey\tvalue\tnotes\n".into()
            };
            text.push_str(&format!("resolution\t{id}\tselected_circular\tfalse\nresolution\t{id}\tselected_order\t.\nresolution\t{id}\toriented_nodes\tn{index}+\nresolution\t{id}\tscore_method\tfixture\nalignment\t{id}\treverse_complemented\tfalse\nalignment\t{id}\trotation_step\t0\n"));
            write_with_parent(&details, &text).unwrap();
            let gfa = case
                .rebuild_out_dir
                .join(id)
                .join(format!("rebuild_{id}.gfa"));
            let fasta = rebuild_fasta_path(case);
            write_with_parent(&gfa, &graph).unwrap();
            write_with_parent(&fasta, &format!(">{id}\n{dna}\n")).unwrap();
            let report = rebuild_report_path(case);
            let stats = rebuild_result_stats_path(case);
            let logs = report.parent().unwrap();
            let support = logs.join(format!("rebuild_{id}_repeat_path_support.tsv"));
            let extract = logs.join(format!("rebuild_{id}_extract.tsv"));
            write_with_parent(&support, "repeat_node\tspanning_reads\n").unwrap();
            write_with_parent(
                &extract,
                &format!("subgraph\trow_type\tnode\n{id}\tnode_copy\tn{index}\n"),
            )
            .unwrap();
            let mut stats_text = "section\tsubgraph\tmetric\tvalue\n".to_string();
            for (key, value) in [
                ("status", "PASS"),
                ("unmapped_segments", "0"),
                ("gap_bases", "0"),
                ("linear_tiling", "PASS"),
                ("topology_contract", "PASS_checked_compaction"),
                ("recorded_path_replay", "PASS_exact_authoritative_sequence"),
                ("all_edge_sequence_overlaps", "PASS"),
                ("shared_node_sequences", "PASS_all_recorded_occurrences"),
                ("merge_mode", "checked_compaction_and_recorded_resolve_path"),
                ("recorded_path_occurrences", "1"),
                ("source_alignment_reverse_complemented", "false"),
                ("target_frame_offset", "0"),
            ] {
                writeln!(stats_text, "summary\t{id}\t{key}\t{value}").unwrap();
            }
            writeln!(
                stats_text,
                "summary\t{id}\tedited_gfa\t{}\nsummary\t{id}\tpolished_fasta\t{}",
                source.display(),
                fasta.display()
            )
            .unwrap();
            write_with_parent(&stats, &stats_text).unwrap();
            let mut report_text =
                format!("section\tkey\tvalue\nrun\tstatus\tcompleted\nrun\tsubgraph\t{id}\n");
            for (key, path) in [
                ("edited_gfa", &source),
                ("checked_gfa", &case.checked_draft_gfa),
                ("resolve_details", &details),
                ("polished_fasta", &fasta),
                ("reads", &case.reads),
                ("validate_data_dir", &validation.join("01.data")),
            ] {
                writeln!(report_text, "input\t{key}\t{}", path.display()).unwrap();
            }
            for (key, path) in [
                ("rebuild_gfa", &gfa),
                ("rebuild_fasta", &fasta),
                ("result_stats", &stats),
                ("run_report", &report),
                ("repeat_path_support", &support),
                ("extract_table", &extract),
            ] {
                writeln!(report_text, "output\t{key}\t{}", path.display()).unwrap();
            }
            write_with_parent(&report, &report_text).unwrap();
            write_with_parent(&round.join("logs/report.tsv"), &format!("section\tname\tmetric\tvalue\nsv_eval\tround_1\tsv_support_status\tpass\noutput\tround_1\tdir\t{}\noutput\tvalidation_fasta\tpath\t{}\ninput\treads\tpath\t{}\noutput\treference\tpath\t{}\n", validation.display(), fasta.display(), case.reads.display(), reference.display())).unwrap();
            write_with_parent(
                &case
                    .workflow_dir
                    .join("checkpoint_2/round_2/checkpoint_2.status.tsv"),
                "metric\tvalue\nstatus\tcomplete\n",
            )
            .unwrap();
        }
        config
    }

    #[test]
    fn refine_reference_uses_each_validated_record_and_preserves_the_shared_reference() {
        let dir = test_dir("refine component references");
        let config = refine_reference_fixture(&dir);
        let shared = polish_reference_path(&config.cases[0]);
        let before = fs::read(&shared).unwrap();
        let bundle = prepare_refine_bundle(&config, "mito").unwrap();
        let rows = homorepair_rows(&bundle).unwrap();
        assert_eq!(rows.len(), 2);
        for (row, case) in rows.iter().zip(&config.cases) {
            let reference = polish_round_dir(case, 2).join("01.inputs/rotated_reference.fasta");
            assert_eq!(row["rotated_reference"], reference.display().to_string());
            let contract = fs::read_to_string(case.workflow_dir.join("refine_source.tsv")).unwrap();
            assert!(contract.contains(&format!("\trotated_reference\t{}\n", reference.display())));
            assert!(contract.contains(&format!(
                "\trotated_reference_sha256\t{}\n",
                commands::resolve::sha256_hex(&fs::read(&reference).unwrap())
            )));
        }
        assert_eq!(before, fs::read(&shared).unwrap());
        fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn refine_reference_rejects_wrong_component_record_and_tampered_validation_report() {
        for fault in ["reference", "report"] {
            let dir = test_dir(&format!("refine wrong {fault}"));
            let config = refine_reference_fixture(&dir);
            let round = polish_round_dir(&config.cases[0], 2);
            let other = polish_round_dir(&config.cases[1], 2).join("01.inputs/rotated_reference.fasta");
            let reference = round.join("01.inputs/rotated_reference.fasta");
            let report = round.join("logs/report.tsv");
            if fault == "reference" {
                let preserved = round.join("reference_before.fasta");
                fs::copy(&reference, &preserved).unwrap();
                fs::write(
                    &report,
                    fs::read_to_string(&report)
                        .unwrap()
                        .replace(reference.to_str().unwrap(), preserved.to_str().unwrap()),
                )
                .unwrap();
                fs::copy(&other, &reference).unwrap();
            } else {
                fs::write(
                    &report,
                    fs::read_to_string(&report)
                        .unwrap()
                        .replace(reference.to_str().unwrap(), other.to_str().unwrap()),
                )
                .unwrap();
            }
            let error = prepare_refine_bundle(&config, "mito")
                .unwrap_err()
                .to_string();
            assert!(
                error.contains("component reference differs"),
                "{fault}: {error}"
            );
            assert!(!refine_bundle_path(&config, "mito").exists());
            assert!(!config.cases[0]
                .workflow_dir
                .join("refine_source.tsv")
                .exists());
            fs::remove_dir_all(dir).unwrap();
        }
    }

    #[test]
    fn standalone_postrepair_reference_follows_its_source_round_in_every_entrypoint() {
        let dir = test_dir("postrepair_source_reference");
        let mut config = homorepair_test_config(&dir, true);
        let legacy_reference = polish_reference_path(&config.cases[0]);
        write_with_parent(
            &legacy_reference,
            ">reference_001\nACGT\n>reference_002\nTGCA\n",
        )
        .unwrap();
        let original = fs::read(&legacy_reference).unwrap();
        for round in [1, 2, 3, 10] {
            let args = homorepair_args(&config, &config.cases[0], round);
            assert_eq!(
                option_value(&args, "--reference").unwrap(),
                legacy_reference.to_str()
            );
        }
        let output = dir.join("05.postrepair");
        config.postrepair = Some(WorkflowPostrepairConfig {
            enabled: true,
            out_dir: output.clone(),
            threads: None,
            frequency: 0.2,
            review_files: BTreeMap::new(),
        });
        config.cases[0].postrepair_out_dir = Some(output);
        for supplied in [false, true] {
            config.cases[0].postrepair_source_dir = supplied.then(|| dir.join("prior/round_3"));
            for round in if supplied { vec![3] } else { vec![1, 2, 3, 10] } {
                let case = &config.cases[0];
                let source = repair_source_dir(case, round);
                let reference = source.join("01.inputs/rotated_reference.fasta");
                write_with_parent(&reference, ">reference_001\nACGT\n").unwrap();
                let args = homorepair_args(&config, case, round);
                assert_eq!(
                    option_value(&args, "--source-dir").unwrap(),
                    source.to_str()
                );
                assert_eq!(
                    option_value(&args, "--reference").unwrap(),
                    reference.to_str()
                );
                assert_eq!(
                    args.iter().filter(|value| *value == "--reference").count(),
                    1
                );
                assert!(!args.contains(&legacy_reference.display().to_string()));
            }
            let destination = dir.join(if supplied { "supplied" } else { "current" }).join(format!("run_{}.sh", config.sample));
            let master = write_plan_from_args(
                &config,
                &["--out".into(), destination.display().to_string()],
            )
            .unwrap();
            let manual = master.with_file_name(format!("run_{}_manual.sh", config.sample));
            for path in [&master, &manual] {
                let text = fs::read_to_string(path).unwrap();
                let lines: Vec<_> = text
                    .lines()
                    .filter(|line| line.contains(" postrepair --organelle "))
                    .collect();
                assert!(!lines.is_empty());
                for line in lines {
                    let args: Vec<String> = line.split_whitespace().map(String::from).collect();
                    let source = option_value(&args, "--source-dir").unwrap().unwrap();
                    let reference = option_value(&args, "--reference").unwrap().unwrap();
                    assert_eq!(
                        PathBuf::from(reference),
                        Path::new(source).join("01.inputs/rotated_reference.fasta")
                    );
                    if supplied {
                        assert_eq!(PathBuf::from(source), dir.join("prior/round_3"));
                    }
                }
                assert!(Command::new("bash")
                    .arg("-n")
                    .arg(path)
                    .status()
                    .unwrap()
                    .success());
            }
        }
        assert_eq!(fs::read(&legacy_reference).unwrap(), original);
        fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn disabled_rebuild_comments_every_line_of_postrepair_handoff() {
        let dir = test_dir("workflow_disabled_rebuild_postrepair");
        fs::create_dir_all(&dir).unwrap();
        let path = dir.join("orgraft.workflow.toml");
        fs::write(&path, "[project]\nsample = \"S42\"\nresults_dir = \".\"\n[commands.postrepair]\nenabled = true\n[commands.rebuild]\nenabled = false\n[workflow.case.mito]\norganelle = \"mito\"\nsubgraph = \"subgraph_001\"\n").unwrap();
        let config = WorkflowConfig::from_path(&path).unwrap();
        let master = write_plan_from_args(&config, &[]).unwrap();
        let command = rebuild_command_with_polished(
            &config, &config.cases[0], false, "FINAL_POLISHED_FASTA",
        );
        assert!(command.lines().count() > 1);
        let master_text = fs::read_to_string(&master).unwrap();
        for line in command.lines() {
            assert!(master_text.lines().any(|actual| actual == format!("# {line}")), "uncommented rebuild example: {line}");
        }
        for script in [master, dir.join("run_S42_manual.sh")] {
            let text = fs::read_to_string(&script).unwrap();
            assert!(!text.lines().any(|line| line.trim_start().starts_with("joint_ir_args=")));
            assert!(Command::new("bash").arg("-n").arg(&script).status().unwrap().success());
        }
        fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn disabled_route_plans_exit_with_a_reason_without_science() {
        let dir = test_dir("workflow_disabled_route");
        fs::create_dir_all(&dir).unwrap();
        let config_path = dir.join("orgraft.workflow.toml");
        fs::write(&config_path, "[project]\nsample = \"S42\"\nresults_dir = \"outputs\"\n[commands.recruit]\nenabled = true\n[workflow.case.mito]\nenabled = false\norganelle = \"mito\"\n").unwrap();
        let config = WorkflowConfig::from_path(&config_path).unwrap();
        let master = write_plan_from_args(&config, &[]).unwrap();
        let manual = config.results_dir.join("run_S42_manual.sh");
        for script in [&master, &manual] {
            let text = fs::read_to_string(script).unwrap();
            assert!(!text.contains("ORGRAFT_BIN"));
            assert!(!text.contains(" recruit "));
            let result = Command::new("bash").arg(script).output().unwrap();
            assert_eq!(result.status.code(), Some(1));
            assert!(String::from_utf8_lossy(&result.stderr).contains("every configured case has enabled=false"));
        }
        assert_eq!(master, config.results_dir.join("run_S42.sh"));
        fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn explicit_plan_output_keeps_both_scripts_together_and_refuses_changed_files() {
        let dir = test_dir("workflow_explicit_plan_output");
        fs::create_dir_all(&dir).unwrap();
        let config_path = dir.join("orgraft.workflow.toml");
        fs::write(&config_path, "[project]\nsample = \"S42\"\nresults_dir = \"outputs\"\n[workflow.case.mito]\norganelle = \"mito\"\n").unwrap();
        let config = WorkflowConfig::from_path(&config_path).unwrap();
        let master = dir.join("route/run_S42.sh");
        let manual = dir.join("route/run_S42_manual.sh");
        let args = vec!["--out".into(), master.display().to_string()];
        assert_eq!(write_plan_from_args(&config, &args).unwrap(), master);
        let original = fs::read(&master).unwrap();
        let original_manual = fs::read(&manual).unwrap();
        write_plan_from_args(&config, &args).unwrap();
        assert_eq!(fs::read(&master).unwrap(), original);
        assert!(!config.results_dir.join("run_S42_manual.sh").exists());
        fs::write(&manual, "manual changes\n").unwrap();
        let error = write_plan_from_args(&config, &args).unwrap_err().to_string();
        assert!(error.contains("differs from existing"));
        assert_eq!(fs::read(&master).unwrap(), original);
        assert_eq!(fs::read_to_string(&manual).unwrap(), "manual changes\n");
        fs::write(&manual, original_manual).unwrap();
        fs::write(&master, "master changes\n").unwrap();
        assert!(write_plan_from_args(&config, &args).is_err());
        assert_eq!(fs::read_to_string(&master).unwrap(), "master changes\n");
        let sentinel = config.results_dir.join("accepted.fasta");
        fs::create_dir_all(&config.results_dir).unwrap();
        fs::write(&sentinel, ">accepted\nACGT\n").unwrap();
        let mut force_args = args.clone();
        force_args.push("--force".into());
        write_plan_from_args(&config, &force_args).unwrap();
        assert_eq!(fs::read(&master).unwrap(), original);
        assert_eq!(fs::read_to_string(&sentinel).unwrap(), ">accepted\nACGT\n");
        assert!(!fs::read_to_string(&master).unwrap().contains("--force"));
        #[cfg(unix)] {
            fs::remove_file(&manual).unwrap();
            std::os::unix::fs::symlink(&sentinel, &manual).unwrap();
            assert!(write_plan_from_args(&config, &force_args).is_err());
            assert_eq!(fs::read_to_string(&sentinel).unwrap(), ">accepted\nACGT\n");
        }
        fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn refine_representation_is_consistent_per_organelle_in_all_entrypoints() {
        let dir = test_dir("workflow_refine_representation");
        fs::create_dir_all(&dir).unwrap();
        let path = dir.join("orgraft.workflow.toml");
        fs::write(
            &path,
            r#"[project]
    sample = "S42"
    results_dir = "outputs"
    [commands.rebuild]
    enabled = true
    [commands.refine]
    enabled = true
    organelle_reads = "reads/${organelle}.fastq.gz"
    [workflow.case.first]
    organelle = "mito"
    subgraph = "subgraph_001"
    refine_downgrade_representation = "linear-compaction"
    [workflow.case.second]
    organelle = "mito"
    subgraph = "subgraph_002"
    refine_downgrade_representation = "linear-compaction"
    [workflow.case.disabled]
    enabled = false
    organelle = "mito"
    subgraph = "subgraph_003"
    [workflow.case.plastid]
    organelle = "plastid"
    subgraph = "subgraph_001"
    "#,
        )
        .unwrap();
        let config = WorkflowConfig::from_path(&path).unwrap();
        for (organelle, expected) in [("mito", "linear-compaction"), ("plastid", "recorded")] {
            let args = refine_args(&config, organelle).unwrap();
            assert_eq!(
                args.iter()
                    .filter(|v| *v == "--downgrade-representation")
                    .count(),
                1
            );
            assert!(args
                .windows(2)
                .any(|v| v[0] == "--downgrade-representation" && v[1] == expected));
            let out = refine_out_dir(&config, organelle);
            write_with_parent(
                &out.join("auto_complete.tsv"),
                "key\tvalue\nstatus\tno_change\n",
            )
            .unwrap();
            write_with_parent(
                &out.join("auto_inputs.tsv"),
                &format!("key\tvalue\nmeta:downgrade_representation\t{expected}\n"),
            )
            .unwrap();
            check_refine_status(&config, organelle).unwrap();
        }
        let master = write_plan_from_args(&config, &[]).unwrap();
        let manual = master.with_file_name("run_S42_manual.sh");
        for script in [&master, &manual] {
            let text = fs::read_to_string(script).unwrap();
            let commands: Vec<_> = text
                .lines()
                .filter(|line| line.contains(" refine --sample-bundle "))
                .collect();
            assert_eq!(commands.len(), 2);
            assert!(commands
                .iter()
                .any(|line| line.contains("--downgrade-representation linear-compaction")));
            assert!(commands
                .iter()
                .any(|line| line.contains("--downgrade-representation recorded")));
            assert!(Command::new("bash")
                .arg("-n")
                .arg(script)
                .status()
                .unwrap()
                .success());
        }
        let mut inconsistent = config.clone();
        inconsistent
            .cases
            .iter_mut()
            .find(|c| c.name == "second")
            .unwrap()
            .refine_downgrade_representation = commands::refine::DowngradeRepresentation::Recorded;
        assert!(refine_args(&inconsistent, "mito").is_err());
        assert!(prepare_refine_bundle(&inconsistent, "mito")
            .unwrap_err()
            .to_string()
            .contains("consistent refine_downgrade_representation"));
        fs::write(
            &path,
            fs::read_to_string(&path)
                .unwrap()
                .replace("linear-compaction", "guess"),
        )
        .unwrap();
        assert!(WorkflowConfig::from_path(&path).is_err());
        fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn refine_intra_measurement_is_explicit_per_organelle_and_shared_by_entrypoints() {
        let dir = test_dir("refine full copy");
        fs::create_dir_all(&dir).unwrap();
        let path = dir.join("orgraft.workflow.toml");
        let text = "[project]\nsample=\"measurement_fixture\"\nresults_dir=\"outputs\"\n[commands.rebuild]\nenabled=true\n[commands.refine]\nenabled=true\n[commands.refine.intra_measurement]\nmito=\"full-copy\"\n[workflow.case.mito]\norganelle=\"mito\"\n[workflow.case.plastid]\norganelle=\"plastid\"\n";
        fs::write(&path, text).unwrap();
        let config = WorkflowConfig::from_path(&path).unwrap();
        let args = refine_args(&config, "mito").unwrap();
        assert_eq!(args.iter().filter(|arg| *arg == "--intra-measurement").count(), 1);
        assert!(args.windows(2).any(|pair| pair == ["--intra-measurement", "full-copy"]));
        assert!(!refine_args(&config, "plastid").unwrap().iter().any(|arg| arg == "--intra-measurement"));
        let main = write_plan_from_args(&config, &[]).unwrap();
        let manual = main.with_file_name("run_measurement_fixture_manual.sh");
        for file in [&main, &manual] {
            let script = fs::read_to_string(file).unwrap();
            let commands = script.lines().filter(|line| line.contains(" refine --sample-bundle ")).collect::<Vec<_>>();
            assert_eq!(commands.len(), 2);
            let mito = commands.iter().find(|line| line.contains("--organelle mito")).unwrap();
            let plastid = commands.iter().find(|line| line.contains("--organelle plastid")).unwrap();
            assert!(mito.contains("--intra-measurement full-copy"));
            assert!(!plastid.contains("--intra-measurement"));
            assert!(Command::new("bash").arg("-n").arg(file).status().unwrap().success());
        }
        fs::write(&path, text.replace("full-copy", "best-count")).unwrap();
        assert!(WorkflowConfig::from_path(&path).is_err());
        fs::write(&path, text.replace("mito=\"full-copy\"", "nucleus=\"full-copy\"")).unwrap();
        assert!(WorkflowConfig::from_path(&path).is_err());
        fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn refine_duplicate_template_full_span_is_explicit_per_organelle_in_all_entrypoints() {
        let dir = test_dir("refine duplicate full span");
        fs::create_dir_all(&dir).unwrap();
        let path = dir.join("orgraft.workflow.toml");
        let text = "[project]\nsample=\"fixture\"\nresults_dir=\"outputs\"\n[commands.rebuild]\nenabled=true\n[commands.refine]\nenabled=true\n[commands.refine.duplicate_template_full_span]\nmito=true\nplastid=false\n[workflow.case.mito]\norganelle=\"mito\"\n[workflow.case.plastid]\norganelle=\"plastid\"\n";
        fs::write(&path, text).unwrap();
        let config = WorkflowConfig::from_path(&path).unwrap();
        assert_eq!(
            refine_args(&config, "mito")
                .unwrap()
                .iter()
                .filter(|a| *a == "--duplicate-template-full-span")
                .count(),
            1
        );
        assert!(!refine_args(&config, "plastid")
            .unwrap()
            .iter()
            .any(|a| a == "--duplicate-template-full-span"));
        let main = write_plan_from_args(&config, &[]).unwrap();
        let manual = main.with_file_name("run_fixture_manual.sh");
        for file in [&main, &manual] {
            let script = fs::read_to_string(file).unwrap();
            let commands = script
                .lines()
                .filter(|line| line.contains(" refine --sample-bundle "))
                .collect::<Vec<_>>();
            assert_eq!(commands.len(), 2);
            assert!(commands
                .iter()
                .find(|line| line.contains("--organelle mito"))
                .unwrap()
                .contains("--duplicate-template-full-span"));
            assert!(!commands
                .iter()
                .find(|line| line.contains("--organelle plastid"))
                .unwrap()
                .contains("--duplicate-template-full-span"));
            assert!(Command::new("bash")
                .arg("-n")
                .arg(file)
                .status()
                .unwrap()
                .success());
        }
        for bad in [
            text.replace("mito=true", "nucleus=true"),
            text.replace("mito=true", "mito=auto"),
        ] {
            fs::write(&path, bad).unwrap();
            assert!(WorkflowConfig::from_path(&path).is_err());
        }
        fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn refine_linear_chain_flanks_is_explicit_per_organelle_in_all_entrypoints() {
        let dir = test_dir("refine linear chain flanks");
        fs::create_dir_all(&dir).unwrap();
        let path = dir.join("orgraft.workflow.toml");
        let text = "[project]\nsample=\"fixture\"\nresults_dir=\"outputs\"\n[commands.rebuild]\nenabled=true\n[commands.refine]\nenabled=true\n[commands.refine.linear_chain_flanks]\nmito=true\nplastid=false\n[workflow.case.mito]\norganelle=\"mito\"\n[workflow.case.plastid]\norganelle=\"plastid\"\n";
        fs::write(&path, text).unwrap();
        let config = WorkflowConfig::from_path(&path).unwrap();
        assert_eq!(
            refine_args(&config, "mito")
                .unwrap()
                .iter()
                .filter(|a| *a == "--linear-chain-flanks")
                .count(),
            1
        );
        assert!(!refine_args(&config, "plastid")
            .unwrap()
            .iter()
            .any(|a| a == "--linear-chain-flanks"));
        let main = write_plan_from_args(&config, &[]).unwrap();
        let manual = main.with_file_name("run_fixture_manual.sh");
        for file in [&main, &manual] {
            let script = fs::read_to_string(file).unwrap();
            let commands = script
                .lines()
                .filter(|line| line.contains(" refine --sample-bundle "))
                .collect::<Vec<_>>();
            assert_eq!(commands.len(), 2);
            assert!(commands
                .iter()
                .find(|line| line.contains("--organelle mito"))
                .unwrap()
                .contains("--linear-chain-flanks"));
            assert!(!commands
                .iter()
                .find(|line| line.contains("--organelle plastid"))
                .unwrap()
                .contains("--linear-chain-flanks"));
            assert!(Command::new("bash")
                .arg("-n")
                .arg(file)
                .status()
                .unwrap()
                .success());
        }
        for bad in [
            text.replace("mito=true", "nucleus=true"),
            text.replace("mito=true", "mito=auto"),
        ] {
            fs::write(&path, bad).unwrap();
            assert!(WorkflowConfig::from_path(&path).is_err());
        }
        fs::remove_dir_all(dir).unwrap();
    }


    fn check_graph_measurement_workflow_all_entrypoints(mode: &str) {
        let dir = test_dir("refine graph full family");
        fs::create_dir_all(&dir).unwrap();
        let path = dir.join("orgraft.workflow.toml");
        let base = "[project]\nsample=\"fixture\"\nresults_dir=\"outputs\"\n[commands.rebuild]\nenabled=true\n[commands.refine]\nenabled=true\n[workflow.case.mito]\norganelle=\"mito\"\n[workflow.case.plastid]\norganelle=\"plastid\"\n";
        fs::write(&path, base).unwrap();
        let default = WorkflowConfig::from_path(&path).unwrap();
        let render = |config: &WorkflowConfig| {
            [
                render_all_cases_plan_script(config).unwrap(),
                render_manual_plan_script(config, &enabled_cases(config)).unwrap(),
            ]
        };
        let original = render(&default);
        let explicit = format!("{base}[commands.refine.graph_measurement]\nmito=\"default\"\nplastid=\"{mode}\"\n");
        fs::write(&path, &explicit).unwrap();
        let full = WorkflowConfig::from_path(&path).unwrap();
        assert_eq!(refine_args(&full, "mito").unwrap(), refine_args(&default, "mito").unwrap());
        let mut args = refine_args(&full, "plastid").unwrap();
        assert_eq!(args.iter().filter(|a| *a == "--graph-measurement").count(), 1);
        args.retain(|a| a != "--graph-measurement" && a != mode);
        assert_eq!(args, refine_args(&default, "plastid").unwrap());
        // These are the actual master/manual renderers; workflow run shares refine_args.
        for (index, (script, old)) in render(&full).into_iter().zip(original).enumerate() {
            assert_eq!(script.matches("--graph-measurement").count(), 1);
            assert_eq!(script.replace(&format!(" --graph-measurement {mode}"), ""), old);
            let line = script.lines().find(|line| line.contains("--graph-measurement")).unwrap();
            assert!(line.contains(" refine --sample-bundle ") && line.contains("--organelle plastid"));
            let file = dir.join(format!("entry_{index}.sh"));
            fs::write(&file, script).unwrap();
            assert!(Command::new("bash").arg("-n").arg(file).status().unwrap().success());
        }
        fs::write(&path, explicit.replace(mode, "default")).unwrap();
        let off = WorkflowConfig::from_path(&path).unwrap();
        assert_eq!(render(&off), render(&default));
        for bad in [
            explicit.replace("plastid=", "nucleus="),
            explicit.replace(mode, "unknown"),
            explicit.replace("[commands.refine]\nenabled=true\n", ""),
        ] {
            fs::write(&path, bad).unwrap();
            assert!(WorkflowConfig::from_path(&path).is_err());
        }
        fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn maximal_graph_measurement_workflow_all_entrypoints() { check_graph_measurement_workflow_all_entrypoints("maximal-exact-span"); }
    #[test]
    fn complete_family_alleles_workflow_all_entrypoints() { check_graph_measurement_workflow_all_entrypoints("complete-family-alleles"); }
    #[test]
    fn multi_context_workflow_all_entrypoints() { check_graph_measurement_workflow_all_entrypoints("multi-context"); }
    #[test]
    fn combined_context_workflow_all_entrypoints() { check_graph_measurement_workflow_all_entrypoints("multi-context-maximal-exact-span"); }

    #[test]
    fn refine_graph_full_family_is_explicit_per_organelle_in_all_entrypoints() {
        let dir = test_dir("refine graph full family");
        fs::create_dir_all(&dir).unwrap();
        let path = dir.join("orgraft.workflow.toml");
        let base = "[project]\nsample=\"fixture\"\nresults_dir=\"outputs\"\n[commands.rebuild]\nenabled=true\n[commands.refine]\nenabled=true\n[workflow.case.mito]\norganelle=\"mito\"\n[workflow.case.plastid]\norganelle=\"plastid\"\n";
        fs::write(&path, base).unwrap();
        let default = WorkflowConfig::from_path(&path).unwrap();
        let render = |config: &WorkflowConfig| {
            [
                render_all_cases_plan_script(config).unwrap(),
                render_manual_plan_script(config, &enabled_cases(config)).unwrap(),
            ]
        };
        let original = render(&default);
        let explicit = format!("{base}[commands.refine.graph_full_family]\nmito=false\nplastid=true\n");
        fs::write(&path, &explicit).unwrap();
        let full = WorkflowConfig::from_path(&path).unwrap();
        assert_eq!(refine_args(&full, "mito").unwrap(), refine_args(&default, "mito").unwrap());
        let mut args = refine_args(&full, "plastid").unwrap();
        assert_eq!(args.iter().filter(|a| *a == "--graph-full-family").count(), 1);
        args.retain(|a| a != "--graph-full-family");
        assert_eq!(args, refine_args(&default, "plastid").unwrap());
        // These are the actual master/manual renderers; workflow run shares refine_args.
        for (index, (script, old)) in render(&full).into_iter().zip(original).enumerate() {
            assert_eq!(script.matches("--graph-full-family").count(), 1);
            assert_eq!(script.replace(" --graph-full-family", ""), old);
            let line = script.lines().find(|line| line.contains("--graph-full-family")).unwrap();
            assert!(line.contains(" refine --sample-bundle ") && line.contains("--organelle plastid"));
            let file = dir.join(format!("entry_{index}.sh"));
            fs::write(&file, script).unwrap();
            assert!(Command::new("bash").arg("-n").arg(file).status().unwrap().success());
        }
        fs::write(&path, explicit.replace("plastid=true", "plastid=false")).unwrap();
        let off = WorkflowConfig::from_path(&path).unwrap();
        assert_eq!(render(&off), render(&default));
        for bad in [
            explicit.replace("plastid=true", "nucleus=true"),
            explicit.replace("plastid=true", "plastid=auto"),
            explicit.replace("[commands.refine]\nenabled=true\n", ""),
        ] {
            fs::write(&path, bad).unwrap();
            assert!(WorkflowConfig::from_path(&path).is_err());
        }
        fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn refine_graph_full_family_rejects_conflicting_methods_for_same_organelle() {
        let dir = test_dir("refine graph full family conflicts");
        fs::create_dir_all(&dir).unwrap();
        let path = dir.join("orgraft.workflow.toml");
        let base = "[project]\nsample=\"fixture\"\nresults_dir=\"outputs\"\n[commands.rebuild]\nenabled=true\n[commands.refine]\nenabled=true\n[commands.refine.graph_full_family]\nplastid=true\n[workflow.case.mito]\norganelle=\"mito\"\n[workflow.case.plastid]\norganelle=\"plastid\"\n";
        for method in ["duplicate_template_full_span", "linear_chain_flanks"] {
            fs::write(&path, format!("{base}[commands.refine.{method}]\nplastid=true\n")).unwrap();
            assert!(WorkflowConfig::from_path(&path).is_err());
            assert!(!dir.join("outputs").exists());
            fs::write(&path, format!("{base}[commands.refine.{method}]\nmito=true\n")).unwrap();
            let separate = WorkflowConfig::from_path(&path).unwrap();
            assert!(refine_args(&separate, "plastid").is_ok());
            assert!(refine_args(&separate, "mito").is_ok());
        }
        fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn refine_graph_full_family_check_rejects_old_or_unbound_execution() {
        let dir = test_dir("refine graph full family completion");
        fs::create_dir_all(&dir).unwrap();
        let path = dir.join("orgraft.workflow.toml");
        let base = "[project]\nsample=\"fixture\"\nresults_dir=\"outputs\"\n[commands.rebuild]\nenabled=true\n[commands.refine]\nenabled=true\n[workflow.case.plastid]\norganelle=\"plastid\"\n";
        fs::write(&path, base).unwrap();
        let old = WorkflowConfig::from_path(&path).unwrap();
        let out = refine_out_dir(&old, "plastid");
        write_with_parent(&out.join("auto_complete.tsv"), "key\tvalue\nstatus\tno_change\n").unwrap();
        check_refine_status(&old, "plastid").unwrap();
        fs::write(&path, format!("{base}[commands.refine.graph_full_family]\nplastid=true\n")).unwrap();
        let full = WorkflowConfig::from_path(&path).unwrap();
        assert!(check_refine_status(&full, "plastid").is_err());
        fs::write(out.join("auto_inputs.tsv"), "key\tvalue\nmeta:graph_full_family\tfalse\n").unwrap();
        assert!(check_refine_status(&full, "plastid").is_err());
        let input = "key\tvalue\nmeta:graph_full_family\ttrue\n";
        fs::write(out.join("auto_inputs.tsv"), input).unwrap();
        assert!(check_refine_status(&full, "plastid").is_err());
        assert!(check_refine_status(&old, "plastid").is_err());
        fs::write(out.join("auto_complete.tsv"), "key\tvalue\nschema\trefine_auto_complete_v2\nstatus\tno_change\ninput_digest\tstale\n").unwrap();
        assert!(check_refine_status(&full, "plastid").is_err());
        fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn refine_intra_measurement_check_keeps_legacy_default_and_binds_explicit_execution() {
        let dir = test_dir("refine measurement check");
        fs::create_dir_all(&dir).unwrap();
        let path = dir.join("orgraft.workflow.toml");
        let base = "[project]\nsample=\"fixture\"\nresults_dir=\"outputs\"\n[commands.rebuild]\nenabled=true\n[commands.refine]\nenabled=true\n[workflow.case.mito]\norganelle=\"mito\"\n";
        fs::write(&path, base).unwrap();
        let old = WorkflowConfig::from_path(&path).unwrap();
        let out = refine_out_dir(&old, "mito");
        write_with_parent(&out.join("auto_complete.tsv"), "key\tvalue\nstatus\tno_change\n").unwrap();
        check_refine_status(&old, "mito").unwrap();
        // A historical file without the new field still means exact-core.
        fs::write(out.join("auto_inputs.tsv"), "key\tvalue\nmeta:sample\tfixture\n").unwrap();
        check_refine_status(&old, "mito").unwrap();
        fs::write(&path, format!("{base}[commands.refine.intra_measurement]\nmito=\"full-copy\"\n")).unwrap();
        let full = WorkflowConfig::from_path(&path).unwrap();
        assert!(check_refine_status(&full, "mito").is_err());
        let input = "key\tvalue\nmeta:intra_measurement\tfull-copy\n";
        fs::write(out.join("auto_inputs.tsv"), input).unwrap();
        assert!(check_refine_status(&full, "mito").is_err());
        let complete = format!("key\tvalue\ninput_digest\t{}\nschema\trefine_auto_complete_v2\nstatus\tno_change\n",
            commands::resolve::sha256_hex(input.as_bytes()));
        fs::write(out.join("auto_complete.tsv"), &complete).unwrap();
        check_refine_status(&full, "mito").unwrap();
        assert!(check_refine_status(&old, "mito").is_err());
        fs::write(out.join("auto_inputs.tsv"), input.replace("full-copy", "exact-core")).unwrap();
        assert!(check_refine_status(&full, "mito").is_err());
        assert!(check_refine_status(&old, "mito").is_err()); // stale completion digest
        fs::write(out.join("auto_inputs.tsv"), format!("{input}meta:intra_measurement\tfull-copy\n")).unwrap();
        assert!(check_refine_status(&full, "mito").is_err());
        fs::write(out.join("auto_inputs.tsv"), input).unwrap();
        fs::write(out.join("auto_complete.tsv"), complete.replace("no_change", "review_required")).unwrap();
        assert!(check_refine_status(&full, "mito").is_err());
        fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn refine_adopt_file_is_organelle_specific_and_shared_by_native_entrypoints() {
        let dir=test_dir("refine adopt with spaces");
        fs::create_dir_all(&dir).unwrap();
        let path=dir.join("orgraft.workflow.toml");
        let text="[project]\nsample=\"fixture\"\nresults_dir=\"outputs\"\n[commands.rebuild]\nenabled=true\n[commands.refine]\nenabled=true\n[commands.refine.adopt]\nmito=\"inputs/mito.refine.tsv\"\n[workflow.case.mito]\norganelle=\"mito\"\n[workflow.case.plastid]\norganelle=\"plastid\"\n";
        fs::write(&path,text).unwrap();
        let config=WorkflowConfig::from_path(&path).unwrap();
        assert!(refine_args(&config,"mito").is_err());
        let decision=dir.join("inputs/mito.refine.tsv");
        write_with_parent(&decision,"repeat_id\taction\tnote\nintra_repeat_0001\tupgrade\tExplicit current review\n").unwrap();
        let expected=fs::canonicalize(&decision).unwrap();
        let args=refine_args(&config,"mito").unwrap();
        assert!(args.windows(2).any(|p|p[0]=="--adopt"&&Path::new(&p[1])==expected));
        assert!(!refine_args(&config,"plastid").unwrap().iter().any(|v|v=="--adopt"));
        let master=write_plan_from_args(&config,&[]).unwrap();
        for script in [&master,&master.with_file_name("run_fixture_manual.sh")] {
            let text=fs::read_to_string(script).unwrap();
            let calls=text.lines().filter(|l|l.contains(" refine --sample-bundle ")).collect::<Vec<_>>();
            assert_eq!(calls.len(),2);
            assert!(calls.iter().find(|l|l.contains("--organelle mito ")).unwrap().contains(&format!("--adopt {}",shell_quote(&expected))));
            assert!(!calls.iter().find(|l|l.contains("--organelle plastid ")).unwrap().contains("--adopt"));
            assert!(Command::new("bash").arg("-n").arg(script).status().unwrap().success());
        }
        fs::remove_file(&decision).unwrap();
        assert!(refine_args(&config,"mito").is_err());
        fs::write(&path,text.replace("mito=\"inputs/mito.refine.tsv\"","unknown=\"x\"")).unwrap();
        assert!(WorkflowConfig::from_path(&path).is_err());
        fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn refine_requires_full_organelle_reads_for_split_components() {
        let dir = test_dir("workflow_refine_full_reads");
        fs::create_dir_all(&dir).unwrap();
        let config_path = dir.join("orgraft.workflow.toml");
        let text = "[project]\nsample = \"S42\"\nresults_dir = \"outputs\"\nsource_results_dir = \"source\"\n[commands.rebuild]\nenabled = true\n[commands.refine]\nenabled = true\n[workflow.case.first]\norganelle = \"mito\"\nsubgraph = \"subgraph_001\"\nreads = \"split/one.fastq.gz\"\n[workflow.case.second]\norganelle = \"mito\"\nsubgraph = \"subgraph_002\"\nreads = \"split/two.fastq.gz\"\n";
        fs::write(&config_path, text).unwrap();
        let mut config = WorkflowConfig::from_path(&config_path).unwrap();
        let error = write_plan_from_args(&config, &[]).unwrap_err().to_string();
        assert!(error.contains("commands.refine.organelle_reads"), "{error}");
        assert!(!config.results_dir.join("run_S42.sh").exists());
        assert!(prepare_refine_bundle(&config, "mito").unwrap_err().to_string().contains("full organelle read pool"));
        config.cases[1].reads = config.cases[0].reads.clone();
        assert_eq!(refine_reads(&config, "mito").unwrap(), config.cases[0].reads);
        fs::write(&config_path, text.replace("[commands.refine]\n", "[commands.refine]\norganelle_reads = \"${source_results_dir}/01.recruit/${organelle}/reads.fastq.gz\"\n")).unwrap();
        let config = WorkflowConfig::from_path(&config_path).unwrap();
        let full_reads = fs::canonicalize(&dir).unwrap().join("source/01.recruit/mito/reads.fastq.gz");
        assert_eq!(refine_reads(&config, "mito").unwrap(), full_reads);
        assert!(refine_args(&config, "mito").unwrap().windows(2).any(|pair| pair[0] == "--reads" && pair[1] == full_reads.display().to_string()));
        let master = write_plan_from_args(&config, &[]).unwrap();
        let text = fs::read_to_string(&master).unwrap();
        assert_eq!(text.matches(" resolve --checked-draft-gfa ").count(), 1);
        assert_eq!(text.matches(" refine --sample-bundle ").count(), 1);
        assert!(text.contains("split/one.fastq.gz") && text.contains("split/two.fastq.gz"));
        assert!(Command::new("bash").arg("-n").arg(&master).status().unwrap().success());
        let mut assemble = config.clone();
        assemble.asm.enabled = true;
        assert!(shared_case_source(&assemble, &assemble.cases[1], &[&assemble.cases[0]]).is_err());
        let mut changed_graph = config.clone();
        changed_graph.cases[1].draft_graph = dir.join("different.gfa");
        assert!(shared_case_source(&changed_graph, &changed_graph.cases[1], &[&changed_graph.cases[0]]).is_err());
        fs::remove_dir_all(dir).unwrap();
    }

    fn joint_ir_handoff_fixture(dir: &Path) -> (WorkflowConfig, PathBuf) {
        let mut config = homorepair_test_config(dir, false);
        let out = dir.join("05.postrepair");
        config.postrepair = Some(WorkflowPostrepairConfig { enabled: true, out_dir: out.clone(), threads: None, frequency: 0.2, review_files: BTreeMap::new() });
        config.cases[0].postrepair_out_dir = Some(out);
        let case = &config.cases[0];
        let source = polished_aln_path(case, 2);
        write_with_parent(&source, ">subgraph_001\nACGT\n").unwrap();
        write_with_parent(&case.checked_draft_gfa, "S\tnode\tACGT\n").unwrap();
        write_with_parent(&resolve_details_path(case), "section\tid\tkey\tvalue\n").unwrap();
        let component = homorepair_out_dir(case).parent().unwrap().to_path_buf();
        let joint = component.join("joint_ir");
        let selected = joint.join("proposal/candidate.fasta");
        let trial = joint.join("proposal/candidate.gfa");
        let evidence = joint.join("proposal/joint_ir.tsv");
        write_with_parent(&selected, ">subgraph_001\nAACGTT\n").unwrap();
        write_with_parent(&trial, "S\tnode\tAACGTT\n").unwrap();
        write_with_parent(&joint.join("source.fasta"), ">subgraph_001\nACGT\n").unwrap();
        write_with_parent(&joint.join("event_candidates.tsv"), "source_pos\tsource_ref\tsource_alt\n1\tA\tAA\n4\tT\tTT\n").unwrap();
        write_with_parent(&joint.join("event_read_classes.tsv"), "step\tread_id\tevent_class\n1\ta\tselected_ALT_event\n1\tb\tselected_ALT_event\n1\te\tuncertain\n1\tf\tuncertain\n2\tc\tselected_ALT_event\n2\td\tselected_ALT_event\n2\te\tuncertain\n2\tf\tuncertain\n").unwrap();
        let round = component.join("round_3");
        let validated = round.join("01.inputs/linear_subgraph.round_3.fasta");
        write_with_parent(&validated, &fs::read_to_string(&selected).unwrap()).unwrap();
        let report = round.join("logs/report.tsv");
        write_with_parent(&report, &format!("section\tname\tmetric\tvalue\nglobal\trun\torganelle\t{}\nglobal\trun\tsubgraph\t{}\noutput\tvalidation_fasta\tpath\t{}\nstage\tsv-eval:round_3\tstatus\tok\nstage\tsnv-indel-eval:round_3\tstatus\tok\nstage\tplot:round_3\tstatus\tok\n",case.organelle,case.subgraph,validated.display())).unwrap();
        write_fake_sv_pass_summary(&round.join("03.validate/03.reports/sv_snv_indel_summary.tsv")).unwrap();
        write_with_parent(&round.join("03.validate/03.reports/snv_indel_high.tsv"), "pos\tref\talt\tcounts\tdepth\ttype\n").unwrap();
        let digest = |path: &Path| commands::resolve::sha256_hex(&fs::read(path).unwrap());
        write_with_parent(&evidence, &format!("key\tvalue\nstatus\tselected_for_ordinary_validation\nR2_path_links\tverified_zero_overlap\nopposite_pair_count\t1\nchecked_gfa_path\t{}\nchecked_gfa_sha256\t{}\nresolve_details_path\t{}\nresolve_details_sha256\t{}\ncandidate_fasta_sha256\t{}\ncandidate_gfa_sha256\t{}\n", case.checked_draft_gfa.display(), digest(&case.checked_draft_gfa), resolve_details_path(case).display(), digest(&resolve_details_path(case)), digest(&selected), digest(&trial))).unwrap();
        let mut summary = format!("metric\tvalue\nstatus\tcomplete\nunresolved_count\t0\nsource_round\t{}\nsource_fasta\t{}\nsource_sequence_sha256\t{}\nselected_fasta\t{}\nselected_sequence_sha256\t{}\nselection_method\tjoint_ir\nfinal_validation_dir\t{}\n",polish_round_dir(case,2).display(),source.display(),commands::resolve::sha256_hex(b"ACGT"),selected.display(),commands::resolve::sha256_hex(b"AACGTT"),round.join("03.validate").display());
        for path in [&evidence, &report, &joint.join("event_candidates.tsv"), &joint.join("event_read_classes.tsv")] {
            writeln!(summary, "evidence_sha256:{}\t{}", path.display(), digest(path)).unwrap();
        }
        write_with_parent(&component.join("postrepair.tsv"), &summary).unwrap();
        (config, validated)
    }

    #[test]
    fn joint_ir_workflow_handoff_uses_only_this_cases_validated_input() {
        let dir = test_dir("joint IR handoff");
        let (config, validated) = joint_ir_handoff_fixture(&dir);
        let case = &config.cases[0];
        let args = postrepair_rebuild_args(&config, case, &validated).unwrap();
        assert_eq!(args.len(), 4);
        assert_eq!(args[0], "--joint-ir-trial-gfa");
        assert_eq!(args[2], "--joint-ir-evidence");
        assert!(args[1].ends_with("joint_ir/proposal/candidate.gfa"));
        let raw = rebuild_args_with_polished(&config, case, false, Some(&validated), None).unwrap();
        assert!(raw.windows(4).any(|window| window == args));
        let same_sequence = dir.join("another_validation.fasta");
        fs::copy(&validated, &same_sequence).unwrap();
        assert!(postrepair_rebuild_args(&config, case, &same_sequence).unwrap_err().to_string().contains("not the selected native validation input"));
        fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn joint_ir_workflow_handoff_rejects_wrong_source_selection_and_missing_files() {
        for fault in ["source_round", "source_fasta", "selection", "missing_trial", "evidence_drift", "checked_path"] {
            let dir = test_dir(&format!("joint handoff {fault}"));
            let (config, validated) = joint_ir_handoff_fixture(&dir);
            let case = &config.cases[0];
            let component = homorepair_out_dir(case).parent().unwrap().to_path_buf();
            let summary = component.join("postrepair.tsv");
            let evidence = component.join("joint_ir/proposal/joint_ir.tsv");
            match fault {
                "source_round" => {
                    let other = dir.join("other/round_2"); fs::create_dir_all(&other).unwrap();
                    let text = fs::read_to_string(&summary).unwrap().replace(&format!("source_round\t{}",polish_round_dir(case,2).display()), &format!("source_round\t{}",other.display()));
                    fs::write(&summary,text).unwrap();
                },
                "source_fasta" => {
                    let other = dir.join("other_source.fasta"); fs::copy(polished_aln_path(case,2),&other).unwrap();
                    let text=fs::read_to_string(&summary).unwrap().replace(&format!("source_fasta\t{}",polished_aln_path(case,2).display()),&format!("source_fasta\t{}",other.display()));fs::write(&summary,text).unwrap();
                },
                "selection" => { let text=fs::read_to_string(&summary).unwrap().replace("selection_method\tjoint_ir","selection_method\tunknown");fs::write(&summary,text).unwrap(); },
                "missing_trial" => fs::remove_file(component.join("joint_ir/proposal/candidate.gfa")).unwrap(),
                "evidence_drift" => { let mut text=fs::read_to_string(&evidence).unwrap();text.push_str("unexpected\trow\n");fs::write(&evidence,text).unwrap(); },
                "checked_path" => {
                    let other=dir.join("other_checked.gfa");fs::copy(&case.checked_draft_gfa,&other).unwrap();
                    let old=commands::resolve::sha256_hex(&fs::read(&evidence).unwrap());
                    let text=fs::read_to_string(&evidence).unwrap().replace(&format!("checked_gfa_path\t{}",case.checked_draft_gfa.display()),&format!("checked_gfa_path\t{}",other.display()));fs::write(&evidence,text).unwrap();
                    let new=commands::resolve::sha256_hex(&fs::read(&evidence).unwrap());
                    let text=fs::read_to_string(&summary).unwrap().replace(&old,&new);fs::write(&summary,text).unwrap();
                },
                _ => unreachable!(),
            }
            assert!(postrepair_rebuild_args(&config,case,&validated).is_err(),"{fault}");
            fs::remove_dir_all(dir).unwrap();
        }
    }

    #[test]
    fn joint_ir_generated_command_defers_handoff_and_preserves_failure_status() {
        let dir = test_dir("joint IR shell with spaces");
        let mut config = homorepair_test_config(&dir, false);
        config.postrepair = Some(WorkflowPostrepairConfig { enabled: true, out_dir: dir.join("05.postrepair"), threads: None, frequency: 0.2, review_files: BTreeMap::new() });
        config.cases[0].postrepair_out_dir=Some(dir.join("05.postrepair"));
        let command = rebuild_command_with_polished_arg(&config,&config.cases[0],false,"\"$final_polished\"",None);
        assert!(command.contains("--rebuild-args --polished-fasta"));
        assert!(!command.contains("eval"));
        assert!(!command.contains("< <("));
        let shell=format!(r#"set -euo pipefail
ORGRAFT_BIN=orgraft_stub
final_polished='/validation path/input.fasta'
orgraft_stub() {{
  if [[ $1 == workflow ]]; then
    if [[ $CHECK_STATUS != 0 ]]; then return "$CHECK_STATUS"; fi
    if [[ $JOINT == 1 ]]; then printf '%s\n' --joint-ir-trial-gfa '/proposal path/candidate.gfa' --joint-ir-evidence '/proposal path/joint_ir.tsv'; fi
    return 0
  fi
  printf '%s\n' "$@"
}}
{command}
"#);
        for (joint,code) in [("1","0"),("0","0"),("1","7")] {
            let result=Command::new("bash").arg("-c").arg(&shell).env("JOINT",joint).env("CHECK_STATUS",code).output().unwrap();
            assert_eq!(result.status.code().unwrap().to_string(),code,"{}",String::from_utf8_lossy(&result.stderr));
            let output=String::from_utf8_lossy(&result.stdout);
            if code!="0" {assert!(output.is_empty());}
            else if joint=="1" {assert!(output.contains("--joint-ir-trial-gfa\n/proposal path/candidate.gfa\n--joint-ir-evidence\n/proposal path/joint_ir.tsv\n"));}
            else {assert!(!output.contains("--joint-ir-"));assert!(output.starts_with("rebuild\n"));}
        }
        fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn explicit_postrepair_rejects_another_source_round_and_missing_native_handoff() {
        let dir = test_dir("workflow_postrepair_source_binding");
        let mut config = homorepair_test_config(&dir, false);
        let out = dir.join("05.postrepair");
        config.postrepair = Some(WorkflowPostrepairConfig { enabled: true, out_dir: out.clone(), threads: None, frequency: 0.2, review_files: BTreeMap::new() });
        config.cases[0].postrepair_out_dir = Some(out);
        let case = &config.cases[0];
        let source = polish_round_dir(case, 2);
        let wrong_source = polish_round_dir(case, 1);
        fs::create_dir_all(&source).unwrap();
        fs::create_dir_all(&wrong_source).unwrap();
        let summary = homorepair_out_dir(case).parent().unwrap().join("postrepair.tsv");
        write_with_parent(&summary, &format!("metric\tvalue\nstatus\tcomplete\nunresolved_count\t0\nsource_round\t{}\n", wrong_source.display())).unwrap();
        assert!(selected_repair_fasta(&config, case, 2).unwrap_err().to_string().contains("does not match"));
        fs::write(&summary, format!("metric\tvalue\nstatus\tcomplete\nunresolved_count\t0\nsource_round\t{}\n", source.display())).unwrap();
        assert!(selected_repair_fasta(&config, case, 2).is_err());
        fs::remove_dir_all(dir).unwrap();
    }

    fn postrepair_source_fixture(dir: &Path, round: usize) -> WorkflowConfig {
        let source = dir.join(format!("old/04.polish/mito/subgraph_001/round_{round}"));
        let fasta = if round == 1 {
            source.join("02.polish/polished_aln.fasta")
        } else {
            source.join(format!("01.inputs/linear_subgraph.round_{round}.fasta"))
        };
        write_with_parent(&fasta, ">subgraph_001\nACGTAC\n").unwrap();
        write_with_parent(
            &source.join("01.inputs/rotated_reference.fasta"),
            ">reference\nACGTAC\n",
        )
        .unwrap();
        write_with_parent(
            &source.join("01.inputs/subgraph_reads.fastq.gz"),
            "fixture reads\n",
        )
        .unwrap();
        write_with_parent(
            &source.join("03.validate/01.data/snv_indel_calls.tsv"),
            "read_id\tpos\tref\talt\n",
        )
        .unwrap();
        write_with_parent(
            &source.join("03.validate/01.data/snv_indel_variants_combined.tsv"),
            "pos\tref\talt\n",
        )
        .unwrap();
        write_fake_sv_pass_summary(&source.join("03.validate/03.reports/sv_snv_indel_summary.tsv"))
            .unwrap();
        write_with_parent(
            &source.join("03.validate/03.reports/snv_indel_high.tsv"),
            "pos\tref\talt\n",
        )
        .unwrap();
        let report = format!("section\tname\tmetric\tvalue\nglobal\trun\torganelle\tmito\nglobal\trun\tsubgraph\tsubgraph_001\nglobal\trun\tvalidate_round\t{round}\nglobal\trun\tper_read_variant_calls\ttrue\noutput\tvalidation_fasta\tpath\t{}\nsv_eval\tround_1\tsv_support_status\tpass\nstage\tsv-eval:round_{round}\tstatus\tok\nstage\tsnv-indel-eval:round_{round}\tstatus\tok\nstage\tplot:round_{round}\tstatus\tok\n", fasta.display());
        write_with_parent(&source.join("logs/report.tsv"), &report).unwrap();
        if round > 1 {
            write_with_parent(
                &source.parent().unwrap().join("round_1/logs/report.tsv"),
                "original first polish provenance\n",
            )
            .unwrap();
        }
        let resolve = dir.join("old/03.resolve_gfa/mito");
        write_with_parent(
            &resolve.join("logs/resolve_details.tsv"),
            "section\tid\tkey\tvalue\n",
        )
        .unwrap();
        write_with_parent(
            &resolve.join("graph/merged_unresolved_subgraph_001.gfa"),
            "S\tn1\tACGTAC\n",
        )
        .unwrap();
        write_with_parent(&dir.join("old/checked.gfa"), "S\tn1\tACGTAC\n").unwrap();
        let path = dir.join("orgraft.workflow.toml");
        write_with_parent(
            &path,
            &format!(
                r#"[project]
sample = "S"
results_dir = "runs"
[commands.recruit]
enabled = false
[commands.asm]
enabled = false
[commands.postrepair]
enabled = true
[commands.rebuild]
enabled = true
[commands.refine]
enabled = true
[workflow.case.mito]
organelle = "mito"
subgraph = "subgraph_001"
postrepair_source_dir = "{}"
reads = "{}"
checked_draft_gfa = "old/checked.gfa"
resolve_out_dir = "old/03.resolve_gfa"
"#,
                source.display(),
                source.join("01.inputs/subgraph_reads.fastq.gz").display()
            ),
        )
        .unwrap();
        WorkflowConfig::from_path(&path).unwrap()
    }

    fn source_tree_bytes(path: &Path) -> BTreeMap<PathBuf, Vec<u8>> {
        fn visit(root: &Path, path: &Path, output: &mut BTreeMap<PathBuf, Vec<u8>>) {
            for entry in fs::read_dir(path).unwrap() {
                let entry = entry.unwrap();
                if entry.file_type().unwrap().is_dir() {
                    visit(root, &entry.path(), output);
                } else {
                    output.insert(
                        entry.path().strip_prefix(root).unwrap().to_path_buf(),
                        fs::read(entry.path()).unwrap(),
                    );
                }
            }
        }
        let mut files = BTreeMap::new();
        visit(path, path, &mut files);
        files
    }

    #[test]
    fn postrepair_source_plan_and_manual_skip_old_stages_and_leave_source_unchanged() {
        let dir = test_dir("postrepair_source_entrypoints");
        let config = postrepair_source_fixture(&dir, 3);
        let case = &config.cases[0];
        let before = source_tree_bytes(&dir.join("old"));
        let master = write_plan_from_args(&config, &[]).unwrap();
        let manual = config.results_dir.join("run_S_manual.sh");
        for file in [&master, &manual] {
            let text = fs::read_to_string(file).unwrap();
            for forbidden in [
                " workflow checkpoint1 ",
                " workflow checkpoint2 ",
                " resolve --",
                " polish --",
                " recruit --",
                " asm --",
            ] {
                assert!(!text.contains(forbidden), "{forbidden}");
            }
            assert!(
                text.contains("--source-dir")
                    && text.contains("round_3")
                    && text.contains("--first-polish-report")
            );
            assert_eq!(text.matches(" postrepair --").count(), 1);
            assert_eq!(text.matches(" rebuild --").count(), 1);
            assert_eq!(text.matches(" refine --").count(), 1);
            assert!(Command::new("bash")
                .arg("-n")
                .arg(file)
                .status()
                .unwrap()
                .success());
        }
        assert_eq!(postrepair_source_round(case).unwrap(), Some(3));
        assert_eq!(completed_round(case).unwrap(), None);
        assert!(checkpoint1_impl(&config, case, true)
            .unwrap_err()
            .to_string()
            .contains("read-only"));
        assert!(checkpoint2_impl(&config, case, 3, true, true, None)
            .unwrap_err()
            .to_string()
            .contains("read-only"));
        assert!(!case.workflow_dir.exists());
        let mut summary = String::new();
        write_runtime_case_statuses(&mut summary, &config).unwrap();
        assert!(summary.contains("not_run (source round input)"));
        assert_eq!(source_tree_bytes(&dir.join("old")), before);
        fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn postrepair_source_scripts_stop_on_review_and_only_run_verified_downstream_commands() {
        let dir = test_dir("postrepair_source_script_flow");
        let config = postrepair_source_fixture(&dir, 2);
        let case = &config.cases[0];
        let scripts = [
            render_plan_script(&config, case, None, false, true, None).unwrap(),
            render_manual_plan_script(&config, &[case]).unwrap(),
        ];
        for (index, text) in scripts.iter().enumerate() {
            let text = text.replace("#!/usr/bin/env bash", "");
            for code in [0, 7] {
                let log = dir.join(format!("commands_{index}_{code}.txt"));
                let setup = format!(
                    r#"set -euo pipefail
orgraft_stub() {{
  printf '%s\n' "$1 ${{2:-}}" >> {}
  if [[ $1 == postrepair ]]; then return {code}; fi
  if [[ $1 == workflow && $2 == postrepair-check && " $* " != *" --rebuild-args "* ]]; then
    printf '%s\n' {}
  fi
  return 0
}}
export -f orgraft_stub
export ORGRAFT_BIN=orgraft_stub
"#,
                    shell_quote(&log),
                    shell_quote(&repair_source_fasta(case, 2))
                );
                let result = Command::new("bash")
                    .arg("-c")
                    .arg(setup + &text)
                    .output()
                    .unwrap();
                assert_eq!(
                    result.status.code(),
                    Some(code),
                    "{}",
                    String::from_utf8_lossy(&result.stderr)
                );
                let commands = fs::read_to_string(log).unwrap();
                assert!(commands.starts_with("workflow postrepair-check\npostrepair --organelle\n"));
                if code != 0 {
                    assert_eq!(commands.lines().count(), 2);
                } else {
                    assert!(
                        commands.contains("workflow postrepair-check")
                            && commands.contains("rebuild --organelle")
                    );
                }
            }
        }
        fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn postrepair_source_rejects_incomplete_native_evidence_and_input_overrides() {
        for fault in [
            "sv",
            "ordinary",
            "per_read",
            "round",
            "identity",
            "stage",
            "validation_fasta",
            "calls",
            "reference",
            "extra_args",
            "recruit",
        ] {
            let dir = test_dir(&format!("postrepair_source_bad_{fault}"));
            let mut config = postrepair_source_fixture(&dir, 2);
            let case = &config.cases[0];
            let source = repair_source_dir(case, 2);
            let report = source.join("logs/report.tsv");
            match fault {
                "sv" => fs::write(
                    source.join("03.validate/03.reports/sv_snv_indel_summary.tsv"),
                    "metric\tvalue\nstatus\tfail\n",
                )
                .unwrap(),
                "ordinary" => fs::write(
                    source.join("03.validate/03.reports/snv_indel_high.tsv"),
                    "pos\tref\talt\n1\tA\tT\n",
                )
                .unwrap(),
                "per_read" => fs::write(
                    &report,
                    fs::read_to_string(&report).unwrap().replace(
                        "per_read_variant_calls\ttrue",
                        "per_read_variant_calls\tfalse",
                    ),
                )
                .unwrap(),
                "round" => fs::write(
                    &report,
                    fs::read_to_string(&report)
                        .unwrap()
                        .replace("validate_round\t2", "validate_round\t1"),
                )
                .unwrap(),
                "identity" => fs::write(
                    &report,
                    fs::read_to_string(&report)
                        .unwrap()
                        .replace("subgraph\tsubgraph_001", "subgraph\tsubgraph_002"),
                )
                .unwrap(),
                "stage" => fs::write(
                    &report,
                    fs::read_to_string(&report).unwrap().replace(
                        "snv-indel-eval:round_2\tstatus\tok",
                        "snv-indel-eval:round_2\tstatus\tfailed",
                    ),
                )
                .unwrap(),
                "validation_fasta" => fs::write(
                    &report,
                    fs::read_to_string(&report)
                        .unwrap()
                        .replace("linear_subgraph.round_2.fasta", "rotated_reference.fasta"),
                )
                .unwrap(),
                "calls" => {
                    fs::remove_file(source.join("03.validate/01.data/snv_indel_calls.tsv")).unwrap()
                }
                "reference" => {
                    let path = dir.join("wrong_reference.fasta");
                    fs::write(&path, ">wrong\nAAA\n").unwrap();
                    config.cases[0].polish_reference = Some(path);
                }
                "extra_args" => {
                    config.polish.extra_args =
                        vec!["--out-dir".into(), source.display().to_string()]
                }
                "recruit" => config.recruit.enabled = true,
                _ => unreachable!(),
            }
            assert!(validate_postrepair_sources(&config).is_err(), "{fault}");
            if fault == "sv" {
                let args = vec!["--config".into(), config.config_path.display().to_string()];
                assert!(run_automatic_workflow(&args).is_err());
                let mut source_check = args;
                source_check.push("--source-only".into());
                assert!(run_homorepair_check(&source_check).is_err());
            }
            assert!(!config.results_dir.exists());
            fs::remove_dir_all(dir).unwrap();
        }
    }

    #[test]
    fn postrepair_source_rejects_output_aliases_and_explicit_checkpoint_execution() {
        let dir = test_dir("postrepair_source_paths");
        let mut config = postrepair_source_fixture(&dir, 1);
        let before = source_tree_bytes(&dir.join("old"));
        let args = vec![
            "--config".into(),
            config.config_path.display().to_string(),
            "--force".into(),
        ];
        assert!(run_checkpoint1(&args).is_err());
        assert!(run_checkpoint2(&args).is_err());
        assert!(write_plan_from_args(
            &config,
            &["--out".into(), dir.join("old/run.sh").display().to_string()]
        )
        .is_err());
        config.cases[0].workflow_dir = repair_source_dir(&config.cases[0], 1).join("workflow");
        assert!(validate_postrepair_sources(&config).is_err());
        config.cases[0].workflow_dir = config.results_dir.join("workflow");
        #[cfg(unix)]
        {
            let child = repair_source_dir(&config.cases[0], 1).join("child");
            fs::create_dir(&child).unwrap();
            let alias = dir.join("alias");
            std::os::unix::fs::symlink(&child, &alias).unwrap();
            config.results_dir = alias.join("../new");
            assert!(validate_postrepair_sources(&config).is_err());
            fs::remove_file(&alias).unwrap();
            std::os::unix::fs::symlink(dir.join("missing"), &alias).unwrap();
            assert!(workflow_future_path(&alias.join("../new")).is_err());
            fs::remove_dir(child).unwrap();
        }
        assert_eq!(source_tree_bytes(&dir.join("old")), before);
        fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn postrepair_source_rejects_binding_aliases_and_source_file_links_into_results() {
        let dir = test_dir("postrepair_source_binding_aliases");
        let config = postrepair_source_fixture(&dir, 2);
        for alias in [
            "--linear-subgraph",
            "--subgraph-fasta",
            "--rotated-reference",
            "--subgraph-reads",
            "--subgraph-id",
            "--workflow-round",
            "--output-dir",
            "--checked-draft-gfa",
            "--draft-gfa",
            "--raw-gfa",
            "--gfa",
            "--verified-fasta",
            "--final-fasta",
            "--fasta",
        ] {
            let mut changed = config.clone();
            changed.polish.extra_args = vec![alias.into(), "different-input".into()];
            assert!(
                validate_postrepair_sources(&changed)
                    .unwrap_err()
                    .to_string()
                    .contains("extra_args"),
                "{alias}"
            );
            changed.polish.extra_args.clear();
            changed.rebuild.extra_args = vec![format!("{alias}=different-input")];
            assert!(
                validate_postrepair_sources(&changed)
                    .unwrap_err()
                    .to_string()
                    .contains("extra_args"),
                "{alias}"
            );
        }
        #[cfg(unix)]
        {
            let source = repair_source_dir(&config.cases[0], 2);
            let target = config.results_dir.join("copied_evidence");
            fs::create_dir_all(&config.results_dir).unwrap();
            for path in [
                source.join("logs/report.tsv"),
                repair_source_fasta(&config.cases[0], 2),
                source.join("03.validate/03.reports/sv_snv_indel_summary.tsv"),
                source.join("03.validate/03.reports/snv_indel_high.tsv"),
            ] {
                let original = fs::read(&path).unwrap();
                fs::write(&target, &original).unwrap();
                fs::remove_file(&path).unwrap();
                std::os::unix::fs::symlink(&target, &path).unwrap();
                let message = validate_postrepair_sources(&config)
                    .unwrap_err()
                    .to_string();
                assert!(
                    message.contains("read-only source evidence"),
                    "{}: {message}",
                    path.display()
                );
                fs::remove_file(&path).unwrap();
                fs::write(&path, original).unwrap();
            }
            validate_postrepair_sources(&config).unwrap();
        }
        fs::remove_dir_all(dir).unwrap();
    }
    #[test]
    fn postrepair_source_handoff_uses_declared_round_without_checkpoint_fabrication() {
        let dir = test_dir("postrepair_source_handoff");
        let (mut config, validated) = joint_ir_handoff_fixture(&dir);
        let case = &mut config.cases[0];
        let source = polish_round_dir(case, 2);
        case.postrepair_source_dir = Some(source.clone());
        case.polish_out_dir = dir.join("new/04.polish");
        assert_eq!(
            selected_repair_fasta(&config, &config.cases[0], 2).unwrap(),
            fs::canonicalize(&validated).unwrap()
        );
        assert!(selected_repair_fasta(&config, &config.cases[0], 1)
            .unwrap_err()
            .to_string()
            .contains("round differs"));
        assert_eq!(
            postrepair_rebuild_args(&config, &config.cases[0], &validated)
                .unwrap()
                .len(),
            4
        );
        assert_eq!(completed_round(&config.cases[0]).unwrap(), None);
        let error = prepare_refine_bundle(&config, "plastid")
            .unwrap_err()
            .to_string();
        assert!(error.contains("rebuild is not complete"), "{error}");
        assert!(!config.cases[0].workflow_dir.join("checkpoint_2").exists());
        fs::remove_dir_all(dir).unwrap();
    }
    #[test]
    fn workflow_starts_with_workflow_config() {
        assert_eq!(DEFAULT_WORKFLOW.first().unwrap().command, "workflow_config");
    }

    #[test]
    fn workflow_ends_with_refine() {
        assert_eq!(DEFAULT_WORKFLOW.last().unwrap().command, "refine");
    }

    #[test]
    fn template_uses_results_workflow_and_checkpoint_sections() {
        let template = workflow_config_template(&TemplateOptions {
            sample: "sample_001".to_string(),
            results_dir: "results_workflow".to_string(),
            soft_paths: "soft_paths.txt".to_string(),
        });
        assert!(template.contains("results_workflow"));
        assert!(!template.contains("source_results_dir"));
        assert!(template.contains("[workflow]"));
        assert!(template.contains("[commands.recruit]"));
        assert!(template.contains("enabled = true"));
        assert!(template.contains("raw_reads = \"/path/to/raw_hifi.fastq.gz\""));
        assert!(template
            .contains("baits = \"mito=/path/to/mito.fasta,plastid=/path/to/plastid.fasta\""));
        assert!(template.contains("[commands.asm]"));
        assert!(template.contains("[commands.resolve]"));
        assert!(template.contains("out_dir = \"${results_dir}/03.resolve_gfa\""));
        assert!(template.contains("[commands.polish]"));
        assert!(template.contains("out_dir = \"${results_dir}/04.polish\""));
        assert!(template.contains("threads = 64"));
        assert!(template.contains("# per_read_variant_calls = true"));
        assert!(template.contains("# sv_plot_highlight_min_reads = 10"));
        assert!(template.contains("[commands.rebuild]"));
        assert!(template.contains("[workflow.case.mito_subgraph_001]"));
        assert!(template.contains("[workflow.case.plastid_subgraph_001]"));
        assert!(template.contains("organelle = \"mito\""));
        assert!(template.contains("organelle = \"plastid\""));
        assert!(template.contains("max_rounds = 3"));
        assert!(template.contains("auto_sv_correction = true"));
        assert!(template.contains("capped at 10"));
        assert!(template.contains("sv_correction_subgroup = \"type_3_subtype_rep_rep_NA:4\""));
        assert!(template.contains("unitig_graph = \"${results_dir}/02.draft_asm/${organelle}/02.anchor_graph_core/02.unitig_graph/graph.gfa\""));
        assert!(!template.contains("command_mode"));
        assert!(template.contains("auto_snv_indel_correction = true"));
        assert!(template.contains(
            "draft_graph = \"${results_dir}/02.draft_asm/${organelle}/03.finalize_graph/graph.gfa\""
        ));
        assert!(template
            .contains("workflow_dir = \"${results_dir}/workflow/${organelle}/${subgraph}\""));
        assert!(template.contains("resolve_out_dir = \"${results_dir}/03.resolve_gfa\""));
        assert!(template.contains("polish_out_dir = \"${results_dir}/04.polish\""));
        assert!(template.contains("rebuild_out_dir = \"${results_dir}/06.rebuild/${organelle}\""));
        assert!(template.contains("checkpoint_1/checked_draft.gfa"));
    }

    #[test]
    fn template_options_reject_removed_command_mode() {
        let mode_args = vec![
            "--command-mode".to_string(),
            "classic".to_string(),
            "--sample".to_string(),
            "S1".to_string(),
        ];
        let err = TemplateOptions::from_args(&mode_args).unwrap_err();
        assert!(err.to_string().contains("no longer support --command-mode"));

        let default_options = TemplateOptions::from_args(&[]).unwrap();
        let default_template = workflow_config_template(&default_options);
        assert!(default_template.contains("sample = \"sample\""));
        assert!(!default_template.contains("command_mode"));

        assert!(!HELP.contains("--command-mode MODE"));
        assert!(!HELP.contains("classic|detailed|concise"));
    }

    #[test]
    fn config_parser_expands_case_paths() {
        let dir = test_dir("workflow_config_parser");
        fs::create_dir_all(&dir).unwrap();
        let config_path = dir.join("orgraft.workflow.toml");
        fs::write(
            &config_path,
            r#"
[project]
sample = "S1"
results_dir = "results_workflow"

[software]
soft_paths = "soft_paths.txt"

[workflow]
max_rounds = 2

[workflow.case]
name = "mito_subgraph_001"
organelle = "mito"
subgraph = "subgraph_001"
draft_graph = "${results_dir}/draft/${organelle}.gfa"
unitig_graph = "${results_dir}/draft/${organelle}.unitig.gfa"
sv_correction_subgroup = "type_3_rep_rep:4"
"#,
        )
        .unwrap();

        let dir = fs::canonicalize(&dir).unwrap();
        let config = WorkflowConfig::from_path(&config_path).unwrap();
        assert_eq!(config.sample, "S1");
        assert_eq!(config.max_rounds, 2);
        assert!(config.auto_sv_correction);
        assert_eq!(
            config.cases[0].sv_correction_subgroup.as_deref(),
            Some("type_3_rep_rep:4")
        );
        assert_eq!(
            config.cases[0].draft_graph,
            dir.join("results_workflow/draft/mito.gfa")
        );
        assert_eq!(
            config.cases[0].unitig_graph,
            Some(dir.join("results_workflow/draft/mito.unitig.gfa"))
        );
        assert_eq!(
            config.cases[0].checked_draft_gfa,
            dir.join("results_workflow/workflow/mito/subgraph_001/checkpoint_1/checked_draft.gfa")
        );
    }

    #[test]
    fn plan_script_gates_rounds_on_checkpoint_status() {
        let dir = test_dir("workflow_plan_script_status_gates");
        fs::create_dir_all(&dir).unwrap();
        let config_path = dir.join("orgraft.workflow.toml");
        fs::write(
            &config_path,
            format!(
                r#"
[project]
sample = "S1"
results_dir = "{}/results_workflow"

[workflow]
max_rounds = 2

[commands.rebuild]
enabled = true

[workflow.case]
name = "mito_subgraph_001"
organelle = "mito"
subgraph = "subgraph_001"
reference = "mito.fa"
"#,
                dir.display()
            ),
        )
        .unwrap();

        let config = WorkflowConfig::from_path(&config_path).unwrap();
        let case = &config.cases[0];
        assert_eq!(
            polish_round_dir(case, 1),
            dir.join("results_workflow")
                .join("04.polish")
                .join("mito")
                .join("subgraph_001")
                .join("round_1")
        );
        assert_eq!(
            polish_round_dir(case, 2),
            dir.join("results_workflow")
                .join("04.polish")
                .join("mito")
                .join("subgraph_001")
                .join("round_2")
        );
        assert_eq!(
            sv_summary_path(case, 2),
            dir.join("results_workflow")
                .join("04.polish")
                .join("mito")
                .join("subgraph_001")
                .join("round_2")
                .join("03.validate")
                .join("03.reports")
                .join("sv_snv_indel_summary.tsv")
        );
        assert_eq!(
            polished_aln_path(case, 2),
            dir.join("results_workflow")
                .join("04.polish")
                .join("mito")
                .join("subgraph_001")
                .join("round_2")
                .join("01.inputs")
                .join("linear_subgraph.round_2.fasta")
        );
        let script_path = dir.join("workflow.commands.sh");
        write_plan_script(&script_path, &config, case, None, false, true).unwrap();

        let script = fs::read_to_string(script_path).unwrap();
        assert!(script.contains("ORGRAFT_BIN="));
        assert!(script.contains("status_value()"));
        assert!(script.contains("checkpoint1_status_file="));
        assert!(script.contains("checkpoint1_checked_gfa="));
        assert!(script.contains("resume_after_checkpoint1=0"));
        assert!(script.contains("checkpoint1 already checked; skip recruit/asm/checkpoint1"));
        assert!(script.contains("if [[ \"$resume_after_checkpoint1\" != \"1\" ]]; then"));
        assert!(script.contains("checkpoint1_status="));
        assert!(script.contains("case \"$checkpoint2_status\" in"));
        assert!(script.contains("next_round_ready)"));
        assert!(script.contains("rebuild_ready)"));
        assert!(script.contains("manual_required)"));
        assert!(script.contains("final_polished="));
        assert!(script.contains("next_round_draft="));
        assert!(script.contains("--repeat-resolution-node \"${repeat_resolution_node}\""));
        assert!(!script.contains("next_round_repolish"));
        assert!(script.contains("\"${ORGRAFT_BIN}\" polish"));
        assert!(!script.contains("/04.polish/round_1/"));
        assert!(!script.contains("/04.polish/round_2/"));
        assert!(script.contains("--validate-round 2"));
        assert!(script.contains("--validate-round 10"));
        assert!(script.contains("workflow complete at checkpoint2 round 1"));
    }

    #[test]
    fn plan_script_expands_workflow_stages() {
        let dir = test_dir("workflow_plan_script_expanded");
        fs::create_dir_all(&dir).unwrap();
        let config_path = dir.join("orgraft.workflow.toml");
        fs::write(
            &config_path,
            format!(
                r#"
[project]
sample = "S1"
results_dir = "{}/results_workflow"

[commands.recruit]
enabled = true
raw_reads = "reads.fastq"
baits = "mito=mito.fa"

[commands.asm]
enabled = true

[commands.resolve]
gfa_editor_mode = "cli"
max_states = 123
max_candidates = 7

[commands.polish]
threads = 12
per_read_variant_calls = false
snv_indel_overlap_policy = "mask-both"
plot_output_format = "both"
sv_plot_highlight_min_reads = 25
snv_indel_plot_high_risk_fraction = 0.7

[workflow.case]
name = "mito_subgraph_001"
organelle = "mito"
subgraph = "subgraph_001"
reference = "mito.fa"
"#,
                dir.display()
            ),
        )
        .unwrap();
        let config = WorkflowConfig::from_path(&config_path).unwrap();
        let case = config.cases[0].clone();

        let script_path = dir.join("workflow.commands.sh");
        write_plan_script(&script_path, &config, &case, None, false, true).unwrap();
        let script = fs::read_to_string(&script_path).unwrap();
        assert!(script.contains("# 01.recruit:"));
        assert!(script.contains("# 02.draft_asm:"));
        assert!(script.contains("# 03.resolve_gfa:"));
        assert!(script.contains("# 04.polish_checkpoint2:"));
        assert!(script.contains("# 06.rebuild:"));
        assert!(script.contains("\"${ORGRAFT_BIN}\" recruit"));
        assert!(script.contains("\"${ORGRAFT_BIN}\" asm"));
        assert!(script.contains(&format!("--image-reference-fasta {}", shell_quote(config.cases[0].reference.as_ref().unwrap()))));
        assert!(script.contains("\"${ORGRAFT_BIN}\" resolve"));
        assert!(script.contains("\"${ORGRAFT_BIN}\" polish"));
        assert!(script.contains("--gfa-editor-mode cli"));
        assert!(script.contains("--max-states 123"));
        assert!(script.contains("--max-candidates 7"));
        assert!(script.contains("--threads 12"));
        assert!(script.contains("--per-read-variant-calls off"));
        assert!(script.contains("--snv-indel-overlap-policy mask-both"));
        assert!(script.contains("--plot-output-format both"));
        assert!(script.contains("--sv-plot-highlight-min-reads 25"));
        assert!(script.contains("--snv-indel-plot-high-risk-fraction 0.7"));
        assert!(!script.contains("\"${ORGRAFT_BIN}\" workflow run"));
        assert!(!script.contains("common args: --platform"));
    }

    #[test]
    fn case_reference_defaults_to_matching_recruit_bait() {
        let dir = test_dir("workflow_case_reference_from_bait");
        fs::create_dir_all(&dir).unwrap();
        let config_path = dir.join("orgraft.workflow.toml");
        fs::write(
            &config_path,
            r#"
[project]
sample = "S1"
results_dir = "results_workflow"

[commands.recruit]
enabled = true
raw_reads = "reads.fastq"
baits = "mito=refs/mito.fa,plastid=refs/plastid.fa"

[workflow.case]
name = "plastid_subgraph_001"
organelle = "plastid"
subgraph = "subgraph_001"
"#,
        )
        .unwrap();

        let dir = fs::canonicalize(&dir).unwrap();
        let config = WorkflowConfig::from_path(&config_path).unwrap();
        assert_eq!(
            config.cases[0].reference,
            Some(dir.join("refs/plastid.fa"))
        );
        assert_eq!(
            config.cases[0].draft_graph,
            dir.join("results_workflow/02.draft_asm/plastid/03.finalize_graph/graph.gfa")
        );
        assert_eq!(
            config.cases[0].reads,
            dir.join("results_workflow/01.recruit/plastid.fastq.gz")
        );
    }

    #[test]
    fn plan_without_case_writes_master_script_for_multiple_cases() {
        let dir = test_dir("workflow_plan_multiple_cases");
        fs::create_dir_all(&dir).unwrap();
        let config_path = dir.join("orgraft.workflow.toml");
        fs::write(
            &config_path,
            format!(
                r#"
[project]
sample = "S1"
results_dir = "{}/results_workflow"

[commands.recruit]
enabled = true
raw_reads = "reads.fastq"
baits = "mito=mito.fa,plastid=plastid.fa"

[commands.asm]
enabled = true

[workflow.case.mito_subgraph_001]
name = "mito_subgraph_001"
organelle = "mito"
subgraph = "subgraph_001"
reference = "mito.fa"

[workflow.case.plastid_subgraph_001]
name = "plastid_subgraph_001"
organelle = "plastid"
subgraph = "subgraph_001"
reference = "plastid.fa"
"#,
                dir.display()
            ),
        )
        .unwrap();
        let config = WorkflowConfig::from_path(&config_path).unwrap();
        let master = dir.join("all.commands.sh");
        write_new_plan_script(&master, &render_all_cases_plan_script(&config).unwrap(), false).unwrap();

        let script = fs::read_to_string(&master).unwrap();
        assert!(script.contains("Generated by orgraft workflow for all configured cases"));
        assert!(script.contains("01.recruit is shared across enabled cases"));
        assert_eq!(script.matches("\"${ORGRAFT_BIN}\" recruit").count(), 1);
        assert_eq!(script.matches("\"${ORGRAFT_BIN}\" asm").count(), 2);
        assert!(script.contains("recruit_ready=1"));
        assert!(script.contains("recruit outputs already exist; skip global recruit"));
        assert!(!script.contains("bash "));
        assert!(script.contains("--organelle mito"));
        assert!(script.contains("--organelle plastid"));
        for case in &config.cases {
            assert!(!case.workflow_dir.join("workflow.commands.sh").exists());
        }
        assert!(Command::new("bash").arg("-n").arg(&master).status().unwrap().success());

    }

    #[test]
    fn master_plan_skips_disabled_cases() {
        let dir = test_dir("workflow_plan_disabled_cases");
        fs::create_dir_all(&dir).unwrap();
        let config_path = dir.join("orgraft.workflow.toml");
        fs::write(
            &config_path,
            format!(
                r#"
[project]
sample = "S1"
results_dir = "{}/results_workflow"

[workflow.case.mito_subgraph_001]
enabled = false
name = "mito_subgraph_001"
organelle = "mito"
subgraph = "subgraph_001"
reference = "mito.fa"

[workflow.case.plastid_subgraph_001]
enabled = true
name = "plastid_subgraph_001"
organelle = "plastid"
subgraph = "subgraph_001"
reference = "plastid.fa"
"#,
                dir.display()
            ),
        )
        .unwrap();
        let config = WorkflowConfig::from_path(&config_path).unwrap();
        let master = dir.join("all.commands.sh");
        write_new_plan_script(&master, &render_all_cases_plan_script(&config).unwrap(), false).unwrap();

        let script = fs::read_to_string(&master).unwrap();
        assert!(!script.contains("mito_subgraph_001"));
        assert!(script.contains("plastid_subgraph_001"));
        assert!(!config.cases[0]
            .workflow_dir
            .join("workflow.commands.sh")
            .exists());
        assert!(!config.cases[1].workflow_dir.join("workflow.commands.sh").exists());

        let args = vec![
            "--case".to_string(),
            "mito_subgraph_001".to_string(),
            "--out".to_string(),
            dir.join("mito.explicit.sh").display().to_string(),
        ];
        let explicit = write_plan_from_args(&config, &args).unwrap();
        let explicit_script = fs::read_to_string(explicit).unwrap();
        assert!(explicit_script.contains("--organelle mito"));
    }

    #[test]
    fn plan_preserves_graph_edited_checkpoint_when_asm_runs() {
        let dir = test_dir("workflow_plan_preserve_graph_edited");
        fs::create_dir_all(&dir).unwrap();
        let config_path = dir.join("orgraft.workflow.toml");
        fs::write(
            &config_path,
            format!(
                r#"
[project]
sample = "S1"
results_dir = "{}/results_workflow"

[commands.recruit]
enabled = false

[commands.asm]
enabled = true

[workflow.case]
name = "mito_subgraph_001"
organelle = "mito"
subgraph = "subgraph_001"
draft_graph = "${{results_dir}}/draft_asm/${{organelle}}/03.finalize_graph/graph.edited.gfa"
reference = "mito.fa"
"#,
                dir.display()
            ),
        )
        .unwrap();
        let config = WorkflowConfig::from_path(&config_path).unwrap();
        let case = &config.cases[0];
        let script_path = dir.join("workflow.commands.sh");
        write_plan_script(&script_path, &config, case, None, false, true).unwrap();

        let script = fs::read_to_string(script_path).unwrap();
        assert!(script.contains("checkpoint1_draft_graph="));
        assert!(script.contains("checkpoint1_draft_backup="));
        assert!(script.contains("draft_graph.input_backup.gfa"));
        assert!(script.contains("cp \"$checkpoint1_draft_graph\" \"$checkpoint1_draft_backup\""));
        assert!(script.contains("cp \"$checkpoint1_draft_backup\" \"$checkpoint1_draft_graph\""));
    }

    #[test]
    fn checkpoint1_auto_checks_simple_graph() {
        let dir = test_dir("workflow_checkpoint1_simple");
        fs::create_dir_all(&dir).unwrap();
        let graph = dir.join("graph.gfa");
        fs::write(&graph, "S\tA\tACGT\nS\tB\tACGT\nL\tA\t+\tB\t+\t0M\n").unwrap();
        let config_path = dir.join("workflow.toml");
        fs::write(
            &config_path,
            format!(
                r#"
[project]
results_dir = "{}"
[workflow.case]
draft_graph = "{}"
checked_draft_gfa = "{}/checked.gfa"
"#,
                dir.display(),
                graph.display(),
                dir.display()
            ),
        )
        .unwrap();
        let config = WorkflowConfig::from_path(&config_path).unwrap();
        let case = &config.cases[0];
        let status = checkpoint1_impl(&config, case, true).unwrap();
        assert_eq!(status, Checkpoint1Status::Checked);
        assert!(case.checked_draft_gfa.exists());
        assert!(case
            .workflow_dir
            .join("checkpoint_1/topology_summary.tsv")
            .exists());
    }

    #[test]
    fn checkpoint1_blocks_missing_link_segment() {
        let dir = test_dir("workflow_checkpoint1_missing_segment");
        fs::create_dir_all(&dir).unwrap();
        let graph = dir.join("graph.gfa");
        fs::write(&graph, "S\tA\tACGT\nL\tA\t+\tB\t+\t0M\n").unwrap();
        let config_path = dir.join("workflow.toml");
        fs::write(
            &config_path,
            format!(
                r#"
[project]
results_dir = "{}"
[workflow.case]
draft_graph = "{}"
checked_draft_gfa = "{}/checked.gfa"
"#,
                dir.display(),
                graph.display(),
                dir.display()
            ),
        )
        .unwrap();
        let config = WorkflowConfig::from_path(&config_path).unwrap();
        let case = &config.cases[0];
        let status = checkpoint1_impl(&config, case, true).unwrap();
        assert_eq!(status, Checkpoint1Status::ManualRequired);
        assert!(!case.checked_draft_gfa.exists());
        assert!(case
            .workflow_dir
            .join("checkpoint_1/manual_edit_required.gfa")
            .exists());
    }

    #[test]
    fn pos_ref_alt_correction_preserves_header() {
        let dir = test_dir("workflow_correction");
        fs::create_dir_all(&dir).unwrap();
        let input = dir.join("input.fasta");
        let output = dir.join("output.fasta");
        fs::write(&input, ">subgraph_001\nAACCGGTT\n").unwrap();
        let edits = vec![
            VariantEdit {
                pos: 3,
                reference: "C".to_string(),
                alternate: "T".to_string(),
            },
            VariantEdit {
                pos: 7,
                reference: "T".to_string(),
                alternate: "TA".to_string(),
            },
        ];

        let summary = apply_edits_to_fasta(&input, &output, &edits).unwrap();
        let corrected = fs::read_to_string(output).unwrap();
        assert!(corrected.starts_with(">subgraph_001\n"));
        assert!(corrected.contains("AATCGGTAT\n"));
        assert!(summary.contains("edit_count=2"));
    }

    #[test]
    fn high_variant_reader_selects_major_multi_allelic_alt() {
        let dir = test_dir("workflow_high_variant_reader");
        fs::create_dir_all(&dir).unwrap();
        let high = dir.join("snv_indel_high.tsv");
        fs::write(
            &high,
            "pos\tref\talt\tcounts\tfixed_ref\n5\tA\tC#G#T\t1#9#2\tNo\n",
        )
        .unwrap();
        let edits = read_high_variant_edits(&high).unwrap();
        assert_eq!(edits.len(), 1);
        assert_eq!(edits[0].alternate, "G");
    }

    #[test]
    fn high_variant_reader_skips_when_reference_support_is_top() {
        let dir = test_dir("workflow_high_variant_reader_ref_top");
        fs::create_dir_all(&dir).unwrap();
        let high = dir.join("snv_indel_high.tsv");
        fs::write(
            &high,
            "pos\tref\talt\ttype\ttotal_count\tdepth\tcounts\tfixed_ref\n311439\tGCCCCCCCCCCCCC\tGCCCCCCCCCCC#GCCCCCCCCCCCC#GCCCCCCCCCCCCCC\tInDel,homopolymer\t92\t156\t5#61#20\tNo\n",
        )
        .unwrap();

        let edits = read_high_variant_edits(&high).unwrap();
        assert!(edits.is_empty());
    }

    #[test]
    fn high_variant_reader_keeps_alt_when_alt_support_exceeds_reference() {
        let dir = test_dir("workflow_high_variant_reader_alt_top");
        fs::create_dir_all(&dir).unwrap();
        let high = dir.join("snv_indel_high.tsv");
        fs::write(
            &high,
            "pos\tref\talt\ttype\ttotal_count\tdepth\tcounts\tfixed_ref\n5\tA\tC#G#T\tSNV\t92\t156\t5#70#17\tYes\n",
        )
        .unwrap();

        let edits = read_high_variant_edits(&high).unwrap();
        assert_eq!(edits.len(), 1);
        assert_eq!(edits[0].alternate, "G");
    }

    #[test]
    fn checkpoint2_can_accept_pass_summary_without_auto_correction() {
        let dir = test_dir("workflow_checkpoint2_no_auto_correction");
        fs::create_dir_all(&dir).unwrap();
        let config_path = dir.join("workflow.toml");
        fs::write(
            &config_path,
            format!(
                r#"
[project]
results_dir = "{}"

[workflow]
max_rounds = 2
auto_snv_indel_correction = false

[workflow.case]
name = "mito_subgraph_001"
"#,
                dir.display()
            ),
        )
        .unwrap();
        let config = WorkflowConfig::from_path(&config_path).unwrap();
        let case = &config.cases[0];
        write_fake_sv_pass_summary(&sv_summary_path(case, 1)).unwrap();
        write_with_parent(
            &snv_indel_high_path(case, 1),
            "pos\tref\talt\tfixed_ref\n5\tA\tG\tNo\n",
        )
        .unwrap();

        let status = checkpoint2_impl(&config, case, 1, true, false, None).unwrap();
        assert_eq!(status, Checkpoint2Status::Complete);
        assert!(!corrected_fasta_path(case, 1).exists());
    }

    #[test]
    fn checkpoint2_completes_when_snv_indel_reference_support_is_top() {
        let dir = test_dir("workflow_checkpoint2_ref_top");
        fs::create_dir_all(&dir).unwrap();
        let config_path = dir.join("workflow.toml");
        fs::write(
            &config_path,
            format!(
                r#"
[project]
results_dir = "{}"

[workflow]
max_rounds = 2
auto_snv_indel_correction = true

[workflow.case]
name = "mito_subgraph_001"
"#,
                dir.display()
            ),
        )
        .unwrap();
        let config = WorkflowConfig::from_path(&config_path).unwrap();
        let case = &config.cases[0];
        write_fake_sv_pass_summary(&sv_summary_path(case, 1)).unwrap();
        write_with_parent(
            &snv_indel_high_path(case, 1),
            "pos\tref\talt\ttype\ttotal_count\tdepth\tcounts\tfixed_ref\n311439\tGCCCCCCCCCCCCC\tGCCCCCCCCCCC#GCCCCCCCCCCCC#GCCCCCCCCCCCCCC\tInDel,homopolymer\t92\t156\t5#61#20\tNo\n",
        )
        .unwrap();

        let status = checkpoint2_impl(&config, case, 1, true, false, None).unwrap();
        assert_eq!(status, Checkpoint2Status::Complete);
        assert!(!corrected_fasta_path(case, 1).exists());
    }


    fn recorded_source_test_config(dir: &Path, stages: Option<(bool, bool)>) -> WorkflowConfig {
        fs::create_dir_all(dir).unwrap();
        let path = dir.join("workflow.toml");
        let stage_sections = stages.map(|(repair, refine)| format!(
            "[commands.postrepair]\nenabled = {repair}\n[commands.refine]\nenabled = {refine}\n"
        )).unwrap_or_default();
        fs::write(&path, format!(
            "[project]\nresults_dir = {}\n[workflow]\nauto_sv_correction = false\n[commands.rebuild]\nenabled = true\nthreads = 16\n{stage_sections}[workflow.case]\nname = \"mito_subgraph_001\"\nreference = \"reference.fasta\"\n",
            toml_string(&dir.display().to_string())
        )).unwrap();
        let mut config = WorkflowConfig::from_path(&path).unwrap();
        // Isolate command construction from the separate postrepair handoff.
        config.cases[0].postrepair_out_dir = None;
        config.cases[0].image_reference_fasta = Some(dir.join("image.fasta"));
        config
    }

    #[test]
    fn recorded_source_enabled_controls_gfa_and_both_commands() {
        let dir = test_dir("workflow_recorded_source_modes");
        let mut config = recorded_source_test_config(&dir, Some((false, false)));
        let repair = config.postrepair.clone().unwrap();
        let refine = config.refine.clone().unwrap();
        // None/None with legacy homorepair must remain the original legacy mode.
        config.polish.homorepair = true;
        let polished = dir.join("polished.fasta");
        let mut legacy_args = None;
        let mut legacy_script = None;
        for stages in [None, Some((false, false)), Some((true, false)), Some((false, true))] {
            config.postrepair = stages.map(|(enabled, _)| WorkflowPostrepairConfig { enabled, ..repair.clone() });
            config.refine = stages.map(|(_, enabled)| WorkflowRefineConfig { enabled, ..refine.clone() });
            let active = stages.is_some_and(|(a, b)| a || b);
            let case = &config.cases[0];
            let expected_source = if active {
                case.resolve_out_dir.join(&case.organelle).join("graph")
                    .join(format!("merged_unresolved_{}.gfa", case.subgraph))
            } else { case.checked_draft_gfa.clone() };
            assert_eq!(rebuild_source_gfa(&config, case), expected_source);
            let args = rebuild_args_with_polished(&config, case, false, Some(&polished), None).unwrap();
            let script = rebuild_command_with_polished_arg(&config, case, false, "POLISHED", None);
            for flag in ["--checked-gfa", "--resolve-details"] {
                assert_eq!(args.iter().filter(|a| a.as_str() == flag).count(), usize::from(active));
                assert_eq!(script.matches(flag).count(), usize::from(active));
            }
            let mut normalized = Vec::new();
            let mut index = 0;
            while index < args.len() {
                if matches!(args[index].as_str(), "--checked-gfa" | "--resolve-details") { index += 2; continue; }
                normalized.push(args[index].clone());
                if args[index] == "--edited-gfa" {
                    normalized.push(case.checked_draft_gfa.display().to_string()); index += 2;
                } else { index += 1; }
            }
            let normalized_script = script.replace(
                &format!(" --checked-gfa {} --resolve-details {}", shell_quote(&case.checked_draft_gfa), shell_quote(&resolve_details_path(case))), ""
            ).replace(&shell_quote(&expected_source), &shell_quote(&case.checked_draft_gfa));
            if let Some(expected) = &legacy_args { assert_eq!(&normalized, expected); }
            else { legacy_args = Some(normalized); }
            if let Some(expected) = &legacy_script { assert_eq!(&normalized_script, expected); }
            else { legacy_script = Some(normalized_script); }
        }
        fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn recorded_source_checkpoint_stops_pending_and_unapplied_ready() {
        for stages in [(true, false), (false, true)] {
            for state in ["pending", "rebuild_ready"] {
                let dir = test_dir("workflow_recorded_source_repeat_stop");
                let config = recorded_source_test_config(&dir, Some(stages));
                let case = &config.cases[0];
                write_fake_sv_pass_summary(&sv_summary_path(case, 1)).unwrap();
                write_fake_snv_high_table(&snv_indel_high_path(case, 1), &[]).unwrap();
                let mismatch = dir.join("mismatch.tsv");
                fs::write(&mismatch, "original mismatch evidence\n").unwrap();
                let marker = RepeatPairingMarker {
                    state: state.into(), repeat_node: "repeat1".into(), source_round: 1,
                    mismatch_report: mismatch.clone(), rebuild_fasta: None, validation_round: None,
                };
                write_repeat_pairing_marker(case, &marker).unwrap();
                let marker_path = repeat_pairing_marker_path(case);
                let marker_before = fs::read(&marker_path).unwrap();
                let checkpoint = case.workflow_dir.join("checkpoint_2/round_1");
                let next = checkpoint.join("next_round.sh");
                write_with_parent(&next, "#!/bin/bash\norgraft rebuild --repeat-resolution-node repeat1\norgraft polish\n").unwrap();
                // run_next=true cannot invoke a native operation on this branch.
                assert_eq!(checkpoint2_impl(&config, case, 1, false, true, None).unwrap(), Checkpoint2Status::ManualRequired);
                assert_eq!(fs::read(&marker_path).unwrap(), marker_before);
                assert_eq!(fs::read_to_string(&mismatch).unwrap(), "original mismatch evidence\n");
                assert_eq!(read_metric_value(&checkpoint.join("checkpoint_2.status.tsv"), "status").unwrap().as_deref(), Some("manual_required"));
                let script = fs::read_to_string(&next).unwrap();
                assert!(script.contains("manual_required"));
                assert!(script.contains("exit 1"));
                assert!(!script.contains("--repeat-resolution-node"));
                assert!(!script.contains("\"${ORGRAFT_BIN}\" polish"));
                assert!(!script.contains("\"${ORGRAFT_BIN}\" rebuild"));
                assert!(!rebuild_result_stats_path(case).exists());
                assert!(!polish_round_dir(case, 2).exists());
                assert!(Command::new("bash").arg("-n").arg(&next).status().unwrap().success());
                assert!(!Command::new("bash").arg(&next).output().unwrap().status.success());
                fs::remove_dir_all(dir).unwrap();
            }
        }
    }

    #[test]
    fn recorded_source_preserves_legacy_repeat_and_completed_validation() {
        let dir = test_dir("workflow_recorded_source_legacy_repeat");
        let mut config = recorded_source_test_config(&dir, Some((false, false)));
        let case = &config.cases[0];
        write_fake_sv_pass_summary(&sv_summary_path(case, 1)).unwrap();
        write_fake_snv_high_table(&snv_indel_high_path(case, 1), &[]).unwrap();
        let mismatch = dir.join("mismatch.tsv");
        fs::write(&mismatch, "original mismatch evidence\n").unwrap();
        write_repeat_pairing_marker(case, &RepeatPairingMarker {
            state: "pending".into(), repeat_node: "repeat1".into(), source_round: 1,
            mismatch_report: mismatch, rebuild_fasta: None, validation_round: None,
        }).unwrap();
        assert_eq!(checkpoint2_impl(&config, case, 1, false, false, None).unwrap(), Checkpoint2Status::RebuildReady);
        let next = case.workflow_dir.join("checkpoint_2/round_1/next_round.sh");
        let script = fs::read_to_string(next).unwrap();
        assert!(script.contains("--repeat-resolution-node"));
        assert!(!script.contains("--checked-gfa"));
        assert!(!script.contains("--resolve-details"));
        assert_eq!(read_repeat_pairing_marker(case).unwrap().unwrap().state, "rebuild_ready");
        // An old ready marker without an actual rebuild must reissue the same
        // legacy operation, not silently become complete.
        assert!(!rebuild_has_repeat_resolution(case, "repeat1").unwrap());
        assert_eq!(checkpoint2_impl(&config, case, 1, false, false, None).unwrap(), Checkpoint2Status::RebuildReady);
        assert_eq!(read_repeat_pairing_marker(case).unwrap().unwrap().state, "rebuild_ready");
        assert!(!polish_round_dir(case, 2).exists());
        // A completed native repeat rebuild must retain the original validation path.
        write_with_parent(&rebuild_result_stats_path(case),
            "section\tid\tname\tkey\tvalue\nsummary\t.\t.\trepeat_resolution_node\trepeat1\n").unwrap();
        config.refine.as_mut().unwrap().enabled = true;
        let case = &config.cases[0];
        write_fake_sv_pass_summary(&sv_summary_path(case, 2)).unwrap();
        write_fake_snv_high_table(&snv_indel_high_path(case, 2), &[]).unwrap();
        assert_eq!(checkpoint2_impl(&config, case, 2, false, false, None).unwrap(), Checkpoint2Status::Complete);
        assert_eq!(read_repeat_pairing_marker(case).unwrap().unwrap().state, "validated");
        fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn repeat_rebuild_validation_counts_as_an_sv_round() {
        let dir = test_dir("workflow_repeat_rebuild_budget");
        fs::create_dir_all(&dir).unwrap();
        let config_path = dir.join("workflow.toml");
        fs::write(
            &config_path,
            format!(
                "[project]\nresults_dir = \"{}\"\n\n[workflow]\nmax_rounds = 3\n\n[workflow.case]\nname = \"mito_subgraph_001\"\n",
                dir.display()
            ),
        )
        .unwrap();
        let config = WorkflowConfig::from_path(&config_path).unwrap();
        let case = &config.cases[0];
        let round1_status = case
            .workflow_dir
            .join("checkpoint_2/round_1/checkpoint_2.status.tsv");
        write_with_parent(
            &round1_status,
            "metric\tvalue\nstatus\trebuild_ready\ncorrection_kind\tsv\nsv_correction_mode\trepeat_pairing_rebuild\nrepeat_node\tutg5\n",
        )
        .unwrap();

        assert_eq!(sv_correction_count(case, 2).unwrap(), 1);
        assert_eq!(
            2usize.saturating_sub(sv_correction_count(case, 2).unwrap()),
            1
        );
        fs::remove_dir_all(dir).unwrap();
    }


    fn homorepair_test_config(dir: &Path, enabled: bool) -> WorkflowConfig {
        fs::create_dir_all(dir).unwrap();
        let path = dir.join("workflow.toml");
        fs::write(&path, format!("[project]\nresults_dir = {}\n[commands.polish]\nhomorepair = {}\n[workflow.case]\nname = \"plastid_subgraph_001\"\norganelle = \"plastid\"\n",
            toml_string(&dir.display().to_string()), u8::from(enabled))).unwrap();
        WorkflowConfig::from_path(&path).unwrap()
    }

    #[test]
    fn homorepair_workflow_config_is_opt_in_and_cannot_enter_ordinary_extra_args() {
        let empty = WorkflowPolishConfig::from_section(None, "S", "results", "results").unwrap();
        assert!(!empty.homorepair);
        assert_eq!(empty.homorepair_frq, 0.2);
        for (key, value) in [("homorepair", "true"), ("homorepair", "2"),
            ("extra_args", "--homorepair 1"), ("extra_args", "-hr=1"),
            ("extra_args", "--homorepair-source old/round_2"),
            ("extra_args", "--homorepair-source-round 3"),
            ("extra_args", "--homorepair-source-round=10"),
            ("extra_args", "--homorepair-source-dir old/round_3"),
            ("extra_args", "--homorepair-source-dir=old/round_10"),
            ("extra_args", "--homorepair-frq 0.25"),
            ("extra_args", "--homorepair-frq=0.25")] {
            let section = BTreeMap::from([(key.to_string(), value.to_string())]);
            assert!(WorkflowPolishConfig::from_section(Some(&section), "S", "results", "results").is_err(), "{key}={value}");
        }
    }

    #[test]
    fn homorepair_frequency_is_validated_and_only_passed_to_repeat_review() {
        for value in ["-0.1", "1.01", "NaN", "inf", "-inf", "ten"] {
            let section = BTreeMap::from([("homorepair_frq".into(), value.into())]);
            assert!(WorkflowPolishConfig::from_section(Some(&section), "S", "results", "results").is_err(), "{value}");
        }
        let dir = test_dir("homorepair frequency");
        let mut config = homorepair_test_config(&dir, true);
        for (value, expected) in [("0", 0.0), ("0.25", 0.25), ("1", 1.0)] {
            let section = BTreeMap::from([("homorepair_frq".into(), value.into())]);
            config.polish = WorkflowPolishConfig::from_section(Some(&section), "S", "results", "results").unwrap();
            assert_eq!(config.polish.homorepair_frq, expected);
            assert!(!config.polish.homorepair);
            config.polish.homorepair = true;
            let case = &config.cases[0];
            let args = homorepair_args(&config, case, 3);
            assert_eq!(args.windows(2).find(|pair| pair[0] == "--homorepair-frq").unwrap()[1], expected.to_string());
            assert_eq!(args.iter().filter(|arg| *arg == "--homorepair-frq").count(), 1);
            assert!(!polish_command(&config, case, 3, &polished_aln_path(case, 3), false).contains("--homorepair"));
        }
        fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn homorepair_plan_is_called_only_after_complete_and_without_force() {
        let dir = test_dir("homorepair_plan");
        let mut config = homorepair_test_config(&dir, true);
        config.force = true;
        config.rebuild.enabled = false;
        let case = &config.cases[0];
        let args = homorepair_args(&config, case, 2);
        assert_eq!(args.windows(2).find(|pair| pair[0] == "--homorepair-frq").unwrap()[1], "0.2");
        for round in [1, 2, 3, 10, 21] {
            let args = homorepair_args(&config, case, round);
            let draft = args.windows(2).find(|pair| pair[0] == "--draft").unwrap()[1].clone();
            assert_eq!(PathBuf::from(draft), polished_aln_path(case, round));
            let source_round = args.windows(2).find(|pair| pair[0] == "--homorepair-source-round").unwrap()[1].clone();
            assert_eq!(source_round, round.to_string());
            assert_eq!(args.iter().filter(|arg| *arg == "--homorepair-source-round").count(), 1);
            for removed in ["--validate-round", "--homorepair", "-hr", "--homorepair-source",
                "--homorepair-variants", "--homorepair-out-dir"] {
                assert!(!args.iter().any(|arg| arg == removed), "unexpected {removed}");
            }
        }
        assert!(!args.iter().any(|arg| arg == "--force"));
        assert!(!args.iter().any(|arg| arg == "--homorepair-max-rounds"));
        let ordinary = polish_command(&config, case, 2, &polished_aln_path(case, 2), true);
        assert!(!ordinary.contains("--homorepair"));
        assert!(ordinary.contains("--validate-round 2"));
        let script_path = dir.join("workflow.sh");
        write_plan_script(&script_path, &config, case, None, false, true).unwrap();
        let script = fs::read_to_string(&script_path).unwrap();
        let mut source_rounds = Vec::new();
        for branch in script.split("    complete)").skip(1) {
            let (complete, rest) = branch.split_once("    next_round_ready)").unwrap();
            assert_eq!(complete.matches("--homorepair-source-round ").count(), 1);
            assert_eq!(complete.matches("--homorepair-frq 0.2").count(), 1);
            let source_round = complete.split_once("--homorepair-source-round ").unwrap().1
                .split_whitespace().next().unwrap().parse::<usize>().unwrap();
            source_rounds.push(source_round);
            assert!(!rest.split("  esac").next().unwrap().contains("--homorepair-source-round "));
            assert!(complete.find("final_polished=").unwrap() < complete.find("--homorepair-source-round ").unwrap());
            assert!(complete.contains("workflow homorepair-check"));
            assert!(!complete.contains("homorepair_status"));
        }
        assert_eq!(source_rounds, (1..=WORKFLOW_HARD_MAX_ROUNDS).collect::<Vec<_>>());
        assert!(Command::new("bash").arg("-n").arg(&script_path).status().unwrap().success());
        config.polish.homorepair = false;
        // A stale optional result cannot affect an ordinary disabled run.
        write_with_parent(&homorepair_outcome_path(&config.cases[0]), "metric\tvalue\nstatus\tcomplete\nelapsed_seconds\t123\n").unwrap();
        let mut summary = String::new();
        write_runtime_homorepair_statuses(&mut summary, &config).unwrap();
        assert!(summary.is_empty());
        assert!(!runtime_stage_rows(&config).unwrap().iter().any(|row| row.stage == "polish homorepair"));
        write_plan_script(&script_path, &config, &config.cases[0], None, false, true).unwrap();
        assert!(!fs::read_to_string(&script_path).unwrap().contains("--homorepair"));
        fs::remove_dir_all(dir).unwrap();
    }

    fn homorepair_test_candidates(case: &WorkflowCase, rows: &[Vec<String>]) {
        let mut text = format!("{HOMOREPAIR_CANDIDATE_HEADER}\n");
        for row in rows { writeln!(text, "{}", row.join("\t")).unwrap(); }
        write_with_parent(&homorepair_out_dir(case).join("03.validate/01.data/candidates.tsv"), &text).unwrap();
    }

    fn homorepair_test_candidate(changed: bool) -> Vec<String> {
        vec!["1", "4", "4", "T", "1", if changed { "2" } else { "1" },
            if changed { "suggest_modify" } else { "support" }, "ALT_bound", "probe_1",
            "4", "3", "0", "100", "80", "60", "20", ".", ".", ".", ".", ".",
            if changed { "3" } else { "." }, if changed { "GT" } else { "." }, if changed { "GTT" } else { "." }]
            .into_iter().map(str::to_string).collect()
    }

    fn homorepair_test_set_candidate(row: &mut [String], key: &str, value: &str) {
        let index = HOMOREPAIR_CANDIDATE_HEADER.split('\t').position(|name| name == key).unwrap();
        row[index] = value.to_string();
    }

    fn homorepair_test_evidence(case: &WorkflowCase, changed: bool) -> (PathBuf, PathBuf) {
        let input = polished_aln_path(case, 2);
        let hr = homorepair_out_dir(case);
        write_with_parent(&input, ">source\nACGT\n").unwrap();
        write_with_parent(&hr.join("01.inputs/linear_subgraph.round_hr.fasta"), ">snapshot\nACGT\n").unwrap();
        let reports = hr.join("03.validate/03.reports");
        let selected = if changed { reports.join("polished.fasta") } else { input.clone() };
        if changed { write_with_parent(&selected, ">selected\nACGTT\n").unwrap(); }
        homorepair_test_candidates(case, &[homorepair_test_candidate(changed)]);
        let final_fasta = if changed { polished_aln_path(case, 3) } else { input.clone() };
        let mut outcome = format!("metric\tvalue\nstatus\tcomplete\nsupport_status\tsupport_verified\nsource_round\t2\ninput_fasta\t{}\nselected_fasta\t{}\nvalidation_status\t{}\nstage_count\t1\nnative_probe_count\t1\nadopted_edit_count\t{}\nunresolved_count\t0\nfinal_validation_status\t{}\n",
            input.display(),selected.display(),if changed {"completed"} else {"skipped_unchanged"},u8::from(changed),if changed {"passed"} else {"skipped_unchanged"});
        if changed {
            write_with_parent(&final_fasta, ">validated\nACGTT\n").unwrap();
            write_fake_sv_pass_summary(&sv_summary_path(case,3)).unwrap();
            write_with_parent(&snv_indel_high_path(case,3), "pos\tref\talt\ttotal_count\tdepth\tcounts\n").unwrap();
            write_with_parent(&polish_report_path(case,3), "section\tname\tmetric\tvalue\nstage\tpolish:round_3\tstatus\tskipped\nstage\tsv-eval:round_3\tstatus\tok\nstage\tsnv-indel-eval:round_3\tstatus\tok\nstage\tplot:round_3\tstatus\tok\n").unwrap();
            writeln!(outcome,"validation_round\t3\nfinal_validation_round\t3\nfinal_validation_dir\t{}\nfinal_validation_fasta\t{}",polish_round_dir(case,3).display(),final_fasta.display()).unwrap();
        }
        let mut run = format!("section\tname\tmetric\tvalue\n{}", outcome.lines().skip(1)
            .map(|line| format!("outcome\thomorepair\t{line}\n")).collect::<String>());
        run.push_str(concat!(
            "outcome\thomorepair\tlocal_batch_count\t0\n",
            "probe\tprobe_1\tevidence_round\t2\n",
            "probe\tprobe_1\tevidence_dir\t03.validate/01.data/probe_1\n",
            "probe\tprobe_1\tprobe_fasta\t03.validate/01.data/probe_1.fasta\n",
            "probe\tprobe_1\tcandidate_count\t1\n",
            "probe\tprobe_1\tSV_status\tpass\n"));
        write_with_parent(&homorepair_outcome_path(case), &run).unwrap();
        (input,final_fasta)
    }

    #[test]
    fn homorepair_handoff_uses_source_when_unchanged_and_final_round_when_changed() {
        for changed in [false,true] {
            let dir = test_dir(if changed {"homorepair changed handoff"} else {"homorepair unchanged handoff"});
            let config = homorepair_test_config(&dir,true);
            let case = &config.cases[0];
            let (source,final_fasta) = homorepair_test_evidence(case,changed);
            assert_eq!(homorepair_selected_fasta(&config,case,2).unwrap(),final_fasta);
            assert!(homorepair_out_dir(case).ends_with("round_hr"));
            if !changed {
                assert_eq!(final_fasta,source);
                assert!(!homorepair_out_dir(case).join("03.validate/03.reports/polished.fasta").exists());
                assert!(!polish_round_dir(case,3).exists());
            }
            fs::remove_dir_all(dir).unwrap();
        }
    }

    #[test]
    fn ordinary_postrepair_wrapper_replays_separate_outputs_and_rejects_changed_evidence() {
        for changed in [false, true] {
            let dir = test_dir(if changed { "postrepair separate changed" } else { "postrepair separate unchanged" });
            let mut config = homorepair_test_config(&dir, true);
            config.cases[0].postrepair_out_dir = Some(dir.join("05.postrepair"));
            let case = &config.cases[0];
            let (_, mut expected) = homorepair_test_evidence(case, changed);
            let component = homorepair_out_dir(case).parent().unwrap().to_path_buf();
            if changed {
                let old_round = polish_round_dir(case, 3);
                let new_round = component.join("round_3");
                fs::rename(&old_round, &new_round).unwrap();
                let outcome = homorepair_outcome_path(case);
                let text = fs::read_to_string(&outcome).unwrap().replace(&old_round.display().to_string(), &new_round.display().to_string());
                fs::write(outcome, text).unwrap();
                expected = new_round.join("01.inputs/linear_subgraph.round_3.fasta");
            }
            let source = polish_round_dir(case, 2);
            assert_eq!(fs::canonicalize(verify_ordinary_postrepair(&component, &source, 0.2).unwrap()).unwrap(), fs::canonicalize(expected).unwrap());
            assert!(verify_ordinary_postrepair(&component, &source, f64::NAN).is_err());
            let wrong_component = component.parent().unwrap().join("subgraph_002");
            fs::create_dir_all(&wrong_component).unwrap();
            assert!(verify_ordinary_postrepair(&wrong_component, &source, 0.2).unwrap_err().to_string().contains("identity"));
            let mut row = homorepair_test_candidate(changed);
            homorepair_test_set_candidate(&mut row, "main_ALT", "20");
            homorepair_test_candidates(case, &[row]);
            assert!(verify_ordinary_postrepair(&component, &source, 0.2).unwrap_err().to_string().contains("support bound"));
            fs::remove_dir_all(dir).unwrap();
        }
    }

    fn declare_review_fixture(config: &mut WorkflowConfig, rows: &str) -> PathBuf {
        let path = config
            .config_path
            .parent()
            .unwrap()
            .join("inputs/human's review.tsv");
        write_with_parent(&path, &format!("{POSTREPAIR_REVIEW_HEADER}\n{rows}")).unwrap();
        config
            .postrepair
            .as_mut()
            .unwrap()
            .review_files
            .insert(config.cases[0].organelle.clone(), path.clone());
        path
    }

    #[test]
    fn declared_postrepair_file_drives_native_keep_and_preserves_receipt_and_automatic_state() {
        let dir = test_dir("declared_review_keep");
        let (mut config, component, _) = manual_keep_fixture(&dir);
        let file = declare_review_fixture(
            &mut config,
            "subgraph_001\tkeep_source\tindel\t3:GT:GTT\tHuman's explicit choice; retain this event\n",
        );
        let summary = fs::read(component.join("postrepair.tsv")).unwrap();
        let run = fs::read(homorepair_outcome_path(&config.cases[0])).unwrap();
        let source = polished_aln_path(&config.cases[0], 2);
        let before = fs::read(&source).unwrap();
        let declaration = declared_postrepair_review(&config, &config.cases[0])
            .unwrap()
            .unwrap();
        let terminal = OrgraftError::InvalidArgument(format!(
            "postrepair requires review; inspect the selected sequence and evidence in {}",
            component.display()
        ));
        assert_eq!(
            finish_declared_postrepair(&config, &config.cases[0], 2, &declaration, Err(terminal))
                .unwrap(),
            fs::canonicalize(&source).unwrap()
        );
        let receipt = fs::read(component.join("manual_review.tsv")).unwrap();
        assert_eq!(
            run_homorepair_after_checkpoint(&config, &config.cases[0], 2).unwrap(),
            fs::canonicalize(&source).unwrap()
        );
        assert_eq!(
            fs::read(component.join("manual_review.tsv")).unwrap(),
            receipt
        );
        assert_eq!(fs::read(component.join("postrepair.tsv")).unwrap(), summary);
        assert_eq!(
            fs::read(homorepair_outcome_path(&config.cases[0])).unwrap(),
            run
        );
        assert_eq!(fs::read(&source).unwrap(), before);
        let text = fs::read_to_string(&file).unwrap();
        for changed in [
            text.replace("Human's explicit choice", "Another choice"),
            text.replace("3:GT:GTT", "3:GT:GTA"),
        ] {
            fs::write(&file, changed).unwrap();
            assert!(selected_repair_fasta(&config, &config.cases[0], 2).is_err());
            assert!(run_homorepair_after_checkpoint(&config, &config.cases[0], 2).is_err());
            assert!(postrepair_rebuild_args(&config, &config.cases[0], &source).is_err());
            assert_eq!(
                fs::read(component.join("manual_review.tsv")).unwrap(),
                receipt
            );
        }
        fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn declared_postrepair_does_not_swallow_errors_or_accept_incomplete_and_wrong_source_outputs() {
        for fault in [
            "io",
            "ordinary_error",
            "incomplete",
            "wrong_round",
            "source",
            "declaration_changed",
        ] {
            let dir = test_dir(&format!("declared_review_{fault}"));
            let (mut config, component, _) = manual_keep_fixture(&dir);
            let file = declare_review_fixture(
                &mut config,
                "subgraph_001\tkeep_source\tindel\t3:GT:GTT\tExplicit keep\n",
            );
            let declaration = declared_postrepair_review(&config, &config.cases[0])
                .unwrap()
                .unwrap();
            let automatic = match fault {
                "io" => Err(OrgraftError::Io(std::io::Error::other(
                    "failed to write report",
                ))),
                "ordinary_error" => Err(OrgraftError::InvalidArgument("alignment failed".into())),
                _ => Ok(()),
            };
            match fault {
                "incomplete" => {
                    fs::remove_file(homorepair_outcome_path(&config.cases[0])).unwrap();
                }
                "source" => {
                    fs::write(
                        polished_aln_path(&config.cases[0], 2),
                        ">subgraph_001\nAAAA\n",
                    )
                    .unwrap();
                }
                "declaration_changed" => {
                    fs::write(
                        &file,
                        fs::read_to_string(&file)
                            .unwrap()
                            .replace("Explicit keep", "Changed decision"),
                    )
                    .unwrap();
                }
                _ => {}
            }
            let result = finish_declared_postrepair(
                &config,
                &config.cases[0],
                if fault == "wrong_round" { 1 } else { 2 },
                &declaration,
                automatic,
            );
            let message = result.unwrap_err().to_string();
            if fault == "io" {
                assert_eq!(message, "failed to write report");
            }
            if fault == "ordinary_error" {
                assert_eq!(message, "alignment failed");
            }
            assert!(!component.join("manual_review.tsv").exists(), "{fault}");
            fs::remove_dir_all(dir).unwrap();
        }
    }

    #[test]
    fn declared_postrepair_current_replays_accepted_edits_and_correction_matches_full_declaration() {
        let dir = test_dir("declared_review_current");
        let (mut config, component, args) = manual_current_fixture(&dir, true);
        let events = args
            .windows(2)
            .filter(|pair| pair[0] == "--retain-indel")
            .map(|pair| {
                format!(
                    "subgraph_001\tkeep_current\tindel\t{}\tExplicit current keep\n",
                    pair[1]
                )
            })
            .collect::<String>();
        declare_review_fixture(&mut config, &events);
        let before = fs::read(component.join("postrepair.tsv")).unwrap();
        let selected = run_homorepair_after_checkpoint(&config, &config.cases[0], 2).unwrap();
        assert_ne!(
            read_single_fasta(&selected).unwrap().1,
            read_single_fasta(&polished_aln_path(&config.cases[0], 2))
                .unwrap()
                .1
        );
        assert_eq!(fs::read(component.join("postrepair.tsv")).unwrap(), before);
        fs::remove_dir_all(dir).unwrap();

        let dir = test_dir("declared_review_correction");
        let (mut config, component, old_args) = manual_correction_fixture(&dir);
        config.soft_paths = dir.join("soft_paths.txt");
        config.postrepair.as_mut().unwrap().threads = Some(1);
        let mut rows = String::from("subgraph_001\tcorrect_snv\tcorrect_snv\t1:A:C\tApproved correction\nsubgraph_001\tcorrect_snv\tretain_indel\t3:GT:GTT\tApproved correction\n");
        for pair in old_args.windows(2).filter(|pair| pair[0] == "--evidence") {
            let relative = Path::new(&pair[1]).strip_prefix(&dir).unwrap();
            writeln!(
                rows,
                "subgraph_001\tcorrect_snv\tevidence_file\t{}\tApproved correction",
                relative.display()
            )
            .unwrap();
        }
        let file = declare_review_fixture(&mut config, &rows);
        let declaration = declared_postrepair_review(&config, &config.cases[0])
            .unwrap()
            .unwrap();
        let mut args = vec!["--result-dir".into(), component.display().to_string()];
        args.extend(declaration.arguments.clone());
        commands::postrepair::manual_correction_fixture(&args, manual_correction_native_fixture)
            .unwrap();
        let selected = run_homorepair_after_checkpoint(&config, &config.cases[0], 2).unwrap();
        assert_eq!(read_single_fasta(&selected).unwrap().1, "CCGT");
        let receipt = fs::read(component.join("manual_review.tsv")).unwrap();
        let original = fs::read_to_string(&file).unwrap();
        for changed in [
            original.replace("1:A:C", "1:A:G"),
            original.replace("Approved correction", "Different reason"),
            original
                .lines()
                .filter(|line| !line.contains("attachment_1.txt"))
                .collect::<Vec<_>>()
                .join("\n"),
        ] {
            fs::write(&file, changed).unwrap();
            assert!(selected_repair_fasta(&config, &config.cases[0], 2).is_err());
        }
        fs::write(&file, original).unwrap();
        config.postrepair.as_mut().unwrap().threads = Some(2);
        assert!(selected_repair_fasta(&config, &config.cases[0], 2).is_err());
        assert_eq!(
            fs::read(component.join("manual_review.tsv")).unwrap(),
            receipt
        );
        fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn declared_postrepair_configuration_and_scripts_share_native_syntax_and_reject_conflicts() {
        let dir = test_dir("declared review ' paths");
        fs::create_dir_all(&dir).unwrap();
        let path = dir.join("workflow.toml");
        let text = "[project]\nsample=\"fixture\"\nresults_dir=\"outputs\"\n[commands.postrepair]\nenabled=true\n[commands.postrepair.review]\nmito=\"inputs/human's review.tsv\"\n[workflow.case.one]\norganelle=\"mito\"\nsubgraph=\"subgraph_001\"\n[workflow.case.two]\norganelle=\"plastid\"\nsubgraph=\"subgraph_001\"\n";
        fs::write(&path, text).unwrap();
        assert!(WorkflowConfig::from_path(&path).is_err());
        let file = dir.join("inputs/human's review.tsv");
        let row = "subgraph_001\tkeep_source\trepeat_extra\t31:41:T:11\tHuman's exact decision\n";
        write_with_parent(&file, &format!("{POSTREPAIR_REVIEW_HEADER}\n{row}")).unwrap();
        let config = WorkflowConfig::from_path(&path).unwrap();
        let review = declared_postrepair_review(&config, &config.cases[0])
            .unwrap()
            .unwrap();
        assert!(review
            .arguments
            .windows(2)
            .any(|pair| pair == ["--repeat-extra", "31:41:T:11"]));
        assert!(declared_postrepair_review(&config, &config.cases[1])
            .unwrap()
            .is_none());
        let main = write_plan_from_args(&config, &[]).unwrap();
        for script in [&main, &main.with_file_name("run_fixture_manual.sh")] {
            let content = fs::read_to_string(script).unwrap();
            assert!(content.contains(&format!(
                "workflow postrepair-run --config {} --case one",
                shell_quote(&config.config_path)
            )));
            let declared_calls = content
                .lines()
                .filter(|line| line.contains("workflow postrepair-run"))
                .collect::<Vec<_>>();
            assert!(!declared_calls.is_empty());
            assert!(declared_calls
                .iter()
                .all(|line| line.contains("--case one ") && !line.contains("||")));
            assert!(content
                .lines()
                .any(|line| line.contains(" postrepair ") && line.contains("--organelle plastid")));
            assert!(Command::new("bash")
                .arg("-n")
                .arg(script)
                .status()
                .unwrap()
                .success());
        }
        for bad in [row.repeat(2), row.replace("31:41:T:11", "031:41:T:11") + row,
            row.replace("subgraph_001", "unknown"), row.replace("keep_source", "guess"),
            row.replace("repeat_extra", "circle"), row.to_string() + &row.replace("exact decision", "other decision"),
            "subgraph_001\tkeep_current\tjunction\ta+:b-:C:GC\tKeep\nsubgraph_001\tkeep_current\tcircle\t1:2:A:2\tKeep\n".into(),
            "subgraph_001\tcorrect_snv\tcorrect_snv\t1:A:C\tKeep\n".into()] {
            fs::write(&file, format!("{POSTREPAIR_REVIEW_HEADER}\n{bad}")).unwrap();
            assert!(WorkflowConfig::from_path(&path).is_err(), "{bad}");
        }
        fs::write(&file, format!("{POSTREPAIR_REVIEW_HEADER}\n{row}")).unwrap();
        fs::write(
            &path,
            text.replace(
                "[workflow.case.one]\n",
                "[workflow.case.one]\nenabled=false\n",
            ),
        )
        .unwrap();
        assert!(WorkflowConfig::from_path(&path).is_err());
        fs::remove_dir_all(dir).unwrap();
    }

    fn manual_keep_fixture(dir: &Path) -> (WorkflowConfig, PathBuf, Vec<String>) {
        let mut config=homorepair_test_config(dir,false);
        let out=dir.join("05.postrepair");
        config.postrepair=Some(WorkflowPostrepairConfig{enabled:true,out_dir:out.clone(),threads:None,frequency:0.2,review_files:BTreeMap::new()});
        config.cases[0].postrepair_out_dir=Some(out);
        let case=&config.cases[0];
        let (source,_)=homorepair_test_evidence(case,false);
        let component=homorepair_out_dir(case).parent().unwrap().to_path_buf();
        let run=homorepair_outcome_path(case);
        let text=fs::read_to_string(&run).unwrap().replace("\tstatus\tcomplete\n","\tstatus\treview_required\n")
            .replace("\tsupport_status\tsupport_verified\n","\tsupport_status\treview_required\n")
            .replace("\tunresolved_count\t0\n","\tunresolved_count\t1\n");
        fs::write(&run,format!("{text}outcome\thomorepair\treason\tnonrepeat_indel_candidates_unconfirmed\noutcome\thomorepair\tindel_candidate_count\t1\noutcome\thomorepair\tindel_adopted_count\t0\noutcome\thomorepair\tindel_review_count\t1\n")).unwrap();
        let data=homorepair_out_dir(case).join("03.validate/01.data");
        let header="step\tsource_row\tsource_rows\tsource_pos\tsource_ref\tsource_alt\tsource_support\tsource_depth\tsource_total_ALT\tprobe_pos\tprobe_delta\tdepth\ttotal_ALT\tmain_ALT\tremainder\tuncertain_reads\tstatus\treason\tnested_repeat_edits\tedit_pos\tedit_ref\tedit_alt\ttrial_candidates\ttested_probes\n";
        for (name,main,remainder,reason) in [("indel_candidates.tsv",4,6,"highest_source_ALT_failed_fixed_pool_gate"),("indel_flank20_candidates.tsv",5,5,"20bp_fixed_pool_gate_failed")] {
            write_with_parent(&data.join(name),&format!("{header}1\t1\t1\t3\tGT\tGTT\t7\t10\t7\t3\t1\t10\t{main}\t{main}\t{remainder}\t1\treview_required\t{reason}\t0\t.\t.\t.\t1\t1\n")).unwrap();
        }
        for (name,alts) in [("read_classes.tsv",4),("flank20_read_classes.tsv",5)] {
            let mut text=String::from("step\tread_id\texact_allele\tclass\n");
            for i in 0..10 {let (allele,class)=if i<alts {("GTT","ALT")} else if i==9 {(".","uncertain")}else{("GT","source")};writeln!(text,"1\tr{i}\t{allele}\t{class}").unwrap();}
            write_with_parent(&data.join("indel").join(name),&text).unwrap();
        }
        let round=polish_round_dir(case,2);
        write_with_parent(&round.join("03.validate/01.data/snv_indel_variants_combined.tsv"),"pos\tref\talt\tcounts\tdepth\ttype\n3\tGT\tGTT\t7\t10\tInDel\n").unwrap();
        write_fake_sv_pass_summary(&sv_summary_path(case,2)).unwrap();
        write_with_parent(&snv_indel_high_path(case,2),"pos\tref\talt\tcounts\tdepth\ttype\n").unwrap();
        write_with_parent(&polish_report_path(case,2),&format!("section\tname\tmetric\tvalue\nglobal\trun\torganelle\t{}\nglobal\trun\tsubgraph\t{}\noutput\tvalidation_fasta\tpath\t{}\nstage\tsv-eval:round_2\tstatus\tok\nstage\tsnv-indel-eval:round_2\tstatus\tok\nstage\tplot:round_2\tstatus\tok\n",case.organelle,case.subgraph,source.display())).unwrap();
        let sha=commands::resolve::sha256_hex(b"ACGT");
        write_with_parent(&component.join("postrepair.tsv"),&format!("metric\tvalue\nstatus\treview_required\nunresolved_count\t1\nselection_method\tordinary\nsource_round\t{}\nsource_fasta\t{}\nselected_fasta\t{}\nsource_sequence_sha256\t{sha}\nselected_sequence_sha256\t{sha}\nfinal_validation_dir\t{}\nmin_alt_fraction\t0.2\n",round.display(),source.display(),source.display(),round.join("03.validate").display())).unwrap();
        let args=vec!["review".into(),"--result-dir".into(),component.display().to_string(),"--keep-source".into(),"--pos".into(),"3".into(),"--ref".into(),"GT".into(),"--alt".into(),"GTT".into(),"--reason".into(),"Current human review: retain the source allele".into()];
        (config,component,args)
    }

    fn manual_current_fixture(dir: &Path, changed: bool) -> (WorkflowConfig, PathBuf, Vec<String>) {
        let (config, component, _) = manual_keep_fixture(dir);
        let case = &config.cases[0];
        let source_round = polish_round_dir(case, 2);
        let source = polished_aln_path(case, 2);
        let hr = homorepair_out_dir(case);
        let data = hr.join("03.validate/01.data");
        let dna = |mut seed: u64, n: usize| -> String {
            (0..n)
                .map(|_| {
                    seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1);
                    b"ACGT"[(seed >> 32) as usize % 4] as char
                })
                .collect()
        };
        let original = format!("ACGT{}CA{}", dna(31, 594), dna(71, 145));
        let current = if changed {
            format!("ACGTT{}", &original[4..])
        } else {
            original.clone()
        };
        let pos = 599 + usize::from(changed);
        let fa = |seq: &str| format!(">{}\n{seq}\n", case.subgraph);
        write_with_parent(&source, &fa(&original)).unwrap();
        write_with_parent(
            &hr.join("01.inputs/linear_subgraph.round_hr.fasta"),
            &fa(&original),
        )
        .unwrap();
        let selected = if changed {
            hr.join("03.validate/03.reports/polished.fasta")
        } else {
            source.clone()
        };
        write_with_parent(&selected, &fa(&current)).unwrap();
        let reads = source_round.join("01.inputs/subgraph_reads.fastq.gz");
        write_with_parent(&reads, "@r0\nACGT\n+\nIIII\n").unwrap();
        let variants = |p: usize| {
            format!("pos\tref\talt\tcounts\tdepth\ttype\tread_ids\n{p}\tCA\tCAAACAAAGAGA#CATGCG\t7#1\t10\tInDel\tr0,r1,r2,r3,r4,r5,r6#r7\n")
        };
        write_with_parent(
            &source_round.join("03.validate/01.data/snv_indel_variants_combined.tsv"),
            &variants(599),
        )
        .unwrap();
        homorepair_test_candidates(case, &[homorepair_test_candidate(changed)]);
        let final_round = if changed {
            component.join("round_3")
        } else {
            source_round.clone()
        };
        let final_fa = if changed {
            final_round.join("01.inputs/linear_subgraph.round_3.fasta")
        } else {
            source.clone()
        };
        write_with_parent(&final_fa, &fa(&current)).unwrap();
        write_fake_sv_pass_summary(
            &final_round.join("03.validate/03.reports/sv_snv_indel_summary.tsv"),
        )
        .unwrap();
        write_with_parent(
            &final_round.join("03.validate/03.reports/snv_indel_high.tsv"),
            "pos\tref\talt\tcounts\tdepth\ttype\n",
        )
        .unwrap();
        write_with_parent(
            &final_round.join("03.validate/01.data/snv_indel_variants_combined.tsv"),
            &variants(pos),
        )
        .unwrap();
        let report = |round: usize, fa: &Path| {
            format!("section\tname\tmetric\tvalue\nglobal\trun\torganelle\t{}\nglobal\trun\tsubgraph\t{}\nglobal\trun\tvalidate_round\t{round}\noutput\tvalidation_fasta\tpath\t{}\nstage\tsv-eval:round_{round}\tstatus\tok\nstage\tsnv-indel-eval:round_{round}\tstatus\tok\nstage\tplot:round_{round}\tstatus\tok\n",case.organelle,case.subgraph,fa.display())
        };
        write_with_parent(&source_round.join("logs/report.tsv"), &report(2, &source)).unwrap();
        write_with_parent(
            &final_round.join("logs/report.tsv"),
            &report(if changed { 3 } else { 2 }, &final_fa),
        )
        .unwrap();
        let mut run = fs::read_to_string(homorepair_outcome_path(case)).unwrap();
        for (k, v) in [
            ("adopted_edit_count", usize::from(changed).to_string()),
            ("unresolved_count", "2".into()),
            ("indel_candidate_count", "2".into()),
            ("indel_review_count", "2".into()),
            ("selected_fasta", selected.display().to_string()),
            (
                "final_validation_status",
                if changed {
                    "passed"
                } else {
                    "skipped_unchanged"
                }
                .into(),
            ),
            (
                "validation_status",
                if changed {
                    "completed"
                } else {
                    "skipped_unchanged"
                }
                .into(),
            ),
        ] {
            run = run
                .lines()
                .map(|l| {
                    let prefix = format!("outcome\thomorepair\t{k}\t");
                    if l.starts_with(&prefix) {
                        format!("{prefix}{v}")
                    } else {
                        l.into()
                    }
                })
                .collect::<Vec<_>>()
                .join("\n")
                + "\n";
        }
        if changed {
            run.push_str(&format!("outcome\thomorepair\tvalidation_round\t3\noutcome\thomorepair\tfinal_validation_round\t3\noutcome\thomorepair\tfinal_validation_dir\t{}\noutcome\thomorepair\tfinal_validation_fasta\t{}\n",final_round.display(),final_fa.display()));
        }
        fs::write(homorepair_outcome_path(case), run).unwrap();
        let header="step\tsource_row\tsource_rows\tsource_pos\tsource_ref\tsource_alt\tsource_support\tsource_depth\tsource_total_ALT\tprobe_pos\tprobe_delta\tdepth\ttotal_ALT\tmain_ALT\tremainder\tuncertain_reads\tstatus\treason\tnested_repeat_edits\tedit_pos\tedit_ref\tedit_alt\ttrial_candidates\ttested_probes\tuntried_exact_reads\tentry_event_support\n";
        for (file, total, main, reason) in [
            (
                "indel_candidates.tsv",
                4,
                3,
                "highest_source_ALT_failed_fixed_pool_gate",
            ),
            (
                "indel_flank20_candidates.tsv",
                5,
                4,
                "20bp_fixed_pool_gate_failed",
            ),
        ] {
            let mut text = header.to_string();
            for (alt, support) in [("CAAACAAAGAGA", 7), ("CATGCG", 1)] {
                text.push_str(&format!("1\t1\t1\t{pos}\tCA\t{alt}\t{support}\t10\t8\t{pos}\t{}\t10\t{total}\t{main}\t{}\t1\treview_required\t{reason}\t0\t.\t.\t.\t2\t1\t1\t{support}\n",alt.len()-2,10-total));
            }
            fs::write(data.join(file), text).unwrap();
        }
        for (file, main) in [("read_classes.tsv", 3), ("flank20_read_classes.tsv", 4)] {
            let mut text = "step\tread_id\texact_allele\tclass\n".to_string();
            for i in 0..10 {
                let (a, c) = if i < main {
                    ("CAAACAAAGAGA", "ALT")
                } else if i == main {
                    ("CATGCG", "ALT")
                } else if i == 9 {
                    (".", "uncertain")
                } else {
                    ("CA", "source")
                };
                text.push_str(&format!("1\tr{i}\t{a}\t{c}\n"));
            }
            fs::write(data.join("indel").join(file), text).unwrap();
        }
        let model = format!("{}CAAACAAAGAGA{}", &current[..pos - 1], &current[pos + 1..]);
        for (name, seq) in [("source", current.as_str()), ("model_001", model.as_str())] {
            let root = data.join("indel/step_001").join(name);
            write_with_parent(&root.join("candidate.fasta"), &fa(seq)).unwrap();
            let native = root.join("round_2");
            let input = native.join("01.inputs/linear_subgraph.round_2.fasta");
            write_with_parent(&input, &fa(seq)).unwrap();
            write_with_parent(&native.join("01.inputs/rotated_reference.fasta"), &fa(seq)).unwrap();
            std::os::unix::fs::symlink(&reads, native.join("01.inputs/subgraph_reads.fastq.gz"))
                .unwrap();
            write_with_parent(&native.join("logs/report.tsv"), &report(2, &input)).unwrap();
            write_fake_sv_pass_summary(&native.join("03.validate/03.reports/sv_snv_indel_summary.tsv"))
                .unwrap();
            if name == "source" {
                write_with_parent(
                    &native.join("03.validate/01.data/snv_indel_variants_combined.tsv"),
                    &variants(pos),
                )
                .unwrap();
            }
        }
        let sha = |s: &str| commands::resolve::sha256_hex(s.as_bytes());
        fs::write(component.join("postrepair.tsv"),format!("metric\tvalue\nstatus\treview_required\nunresolved_count\t2\nselection_method\tordinary\nsource_round\t{}\nsource_fasta\t{}\nselected_fasta\t{}\nsource_sequence_sha256\t{}\nselected_sequence_sha256\t{}\nfinal_validation_dir\t{}\nmin_alt_fraction\t0.2\n",source_round.display(),source.display(),selected.display(),sha(&original),sha(&current),final_round.join("03.validate").display())).unwrap();
        let args = vec![
            "review".into(),
            "--result-dir".into(),
            component.display().to_string(),
            "--keep-current".into(),
            "--retain-indel".into(),
            format!("1:{pos}:CA:CAAACAAAGAGA"),
            "--retain-indel".into(),
            format!("1:{pos}:CA:CATGCG"),
            "--reason".into(),
            "Current human review: retain the current selected sequence at both unresolved ALTs".into(),
        ];
        (config, component, args)
    }

    #[test]
    fn manual_current_replays_a_group_and_preserves_accepted_edits() {
        for changed in [false, true] {
            let dir = test_dir(if changed {
                "current changed"
            } else {
                "current unchanged"
            });
            let (config, component, args) = manual_current_fixture(&dir, changed);
            let summary = fs::read(component.join("postrepair.tsv")).unwrap();
            let run = fs::read(homorepair_outcome_path(&config.cases[0])).unwrap();
            commands::postrepair::run(&args).unwrap();
            assert_eq!(fs::read(component.join("postrepair.tsv")).unwrap(), summary);
            assert_eq!(
                fs::read(homorepair_outcome_path(&config.cases[0])).unwrap(),
                run
            );
            let selected = commands::postrepair::selected_fasta(&component).unwrap();
            assert_eq!(
                read_single_fasta(&selected).unwrap().0,
                config.cases[0].subgraph
            );
            assert!(selected_repair_fasta(&config, &config.cases[0], 2).is_ok());
            let table = component.join("round_hr/03.validate/01.data/indel/read_classes.tsv");
            let original = fs::read(&table).unwrap();
            fs::write(&table, b"changed").unwrap();
            assert!(commands::postrepair::selected_fasta(&component).is_err());
            fs::write(table, original).unwrap();
            assert!(commands::postrepair::run(&args).is_err());
            fs::remove_dir_all(dir).unwrap();
        }
    }

    #[test]
    fn manual_current_rejects_incomplete_sets_changed_frames_and_false_pool_claims() {
        for fault in [
            "missing",
            "duplicate",
            "frame",
            "counts",
            "group",
            "source_ids",
            "record",
            "adopted",
            "review",
            "empty_stage",
        ] {
            let dir = test_dir(fault);
            let (config, component, mut args) = manual_current_fixture(&dir, true);
            let data = component.join("round_hr/03.validate/01.data");
            match fault {
                "missing" => {
                    args.drain(6..8);
                }
                "duplicate" => {
                    args[7] = args[5].clone();
                }
                "frame" => {
                    fs::write(
                        data.join("indel/step_001/source/candidate.fasta"),
                        ">subgraph_001\nACGT\n",
                    )
                    .unwrap();
                }
                "counts" | "group" => {
                    let p = data.join("indel_candidates.tsv");
                    let t = fs::read_to_string(&p).unwrap();
                    fs::write(
                        p,
                        if fault == "counts" {
                            t.replace("\t10\t4\t3\t6\t1\t", "\t10\t4\t4\t6\t1\t")
                        } else {
                            t.replacen("\t10\t4\t3\t6\t1\t", "\t10\t5\t3\t5\t1\t", 1)
                        },
                    )
                    .unwrap();
                }
                "source_ids" => {
                    let p=data.join("indel/step_001/source/round_2/03.validate/01.data/snv_indel_variants_combined.tsv");
                    let t = fs::read_to_string(&p).unwrap();
                    fs::write(p, t.replace("r0,r1", "foreign,r1")).unwrap();
                }
                "record" => {
                    let p = component.join("round_3/01.inputs/linear_subgraph.round_3.fasta");
                    let t = fs::read_to_string(&p).unwrap();
                    fs::write(p, t.replace(">subgraph_001", ">wrong")).unwrap();
                }
                "adopted" => {
                    let p = homorepair_outcome_path(&config.cases[0]);
                    let t = fs::read_to_string(&p).unwrap();
                    fs::write(
                        p,
                        t.replace("adopted_edit_count\t1", "adopted_edit_count\t2"),
                    )
                    .unwrap();
                }
                "empty_stage" => {
                    let p = homorepair_outcome_path(&config.cases[0]);
                    let t = fs::read_to_string(&p).unwrap();
                    fs::write(
                        p,
                        t.replace("indel_candidate_count\t2", "indel_candidate_count\t0"),
                    )
                    .unwrap();
                }
                "review" => {
                    let p = homorepair_outcome_path(&config.cases[0]);
                    let t = fs::read_to_string(&p).unwrap();
                    fs::write(p, t + "outcome\thomorepair\trepeat_extra_review_count\t1\n").unwrap();
                }
                _ => unreachable!(),
            }
            assert!(
                commands::postrepair::run(&args).is_err(),
                "accepted {fault}"
            );
            assert!(!component.join("manual_review.tsv").exists());
            fs::remove_dir_all(dir).unwrap();
        }
    }

    #[test]
    fn manual_current_replays_an_accepted_indel_before_the_retained_step() {
        let dir = test_dir("manual current preceding InDel");
        let (config, component, mut args) = manual_current_fixture(&dir, true);
        let case = &config.cases[0];
        let data = component.join("round_hr/03.validate/01.data");
        let selected = component.join("round_hr/03.validate/03.reports/polished.fasta");
        let prior = read_single_fasta(&selected).unwrap().1;
        let pos = 700usize;
        let reference = &prior[pos - 1..pos];
        let alt = format!("{reference}TGCG");
        let current = format!("{}{alt}{}", &prior[..pos - 1], &prior[pos..]);
        for p in [
            &selected,
            &component.join("round_3/01.inputs/linear_subgraph.round_3.fasta"),
        ] {
            fs::write(p, format!(">{}\n{current}\n", case.subgraph)).unwrap();
        }
        let old = data.join("indel/step_001");
        let new = data.join("indel/step_002");
        fs::rename(&old, &new).unwrap();
        // The existing failed model retains its exact nominated edit, in the
        // independently validated source created after the preceding adoption.
        for name in ["source", "model_001"] {
            let root = new.join(name);
            let previous = read_single_fasta(&root.join("candidate.fasta")).unwrap().1;
            let offset = if name == "model_001" {
                "CAAACAAAGAGA".len() - "CA".len()
            } else {
                0
            };
            let start = pos - 1 + offset;
            assert_eq!(&previous[start..start + 1], reference);
            let next = format!("{}{alt}{}", &previous[..start], &previous[start + 1..]);
            for relative in [
                "candidate.fasta",
                "round_2/01.inputs/linear_subgraph.round_2.fasta",
                "round_2/01.inputs/rotated_reference.fasta",
            ] {
                fs::write(root.join(relative), format!(">{}\n{next}\n", case.subgraph)).unwrap();
            }
            let report = root.join("round_2/logs/report.tsv");
            let text = fs::read_to_string(&report)
                .unwrap()
                .replace(&old.display().to_string(), &new.display().to_string());
            fs::write(report, text).unwrap();
        }
        for file in [
            "indel_candidates.tsv",
            "indel_flank20_candidates.tsv",
            "indel/read_classes.tsv",
            "indel/flank20_read_classes.tsv",
        ] {
            let path = data.join(file);
            let text = fs::read_to_string(&path).unwrap();
            let text = text
                .lines()
                .map(|line| {
                    line.strip_prefix("1\t")
                        .map(|rest| format!("2\t{rest}"))
                        .unwrap_or_else(|| line.into())
                })
                .collect::<Vec<_>>()
                .join("\n")
                + "\n";
            fs::write(path, text).unwrap();
        }
        let candidates = data.join("indel_candidates.tsv");
        let text = fs::read_to_string(&candidates).unwrap();
        let (first, rest) = text.split_once('\n').unwrap();
        fs::write(candidates,format!("{first}\n1\t2\t2\t{pos}\t{reference}\t{alt}\t8\t10\t8\t{pos}\t4\t10\t8\t8\t2\t0\tadopted\tunique_main_ALT_above_fixed_pool_remainder\t0\t{pos}\t{reference}\t{alt}\t1\t1\t0\t8\n{rest}")).unwrap();
        let classes = data.join("indel/read_classes.tsv");
        let text = fs::read_to_string(&classes).unwrap();
        let (h, r) = text.split_once('\n').unwrap();
        let mut next = format!("{h}\n");
        for i in 0..10 {
            next.push_str(&format!(
                "1\ta{i}\t{}\t{}\n",
                if i < 8 { alt.as_str() } else { reference },
                if i < 8 { "ALT" } else { "source" }
            ));
        }
        next.push_str(r);
        fs::write(classes, next).unwrap();
        let run = homorepair_outcome_path(case);
        let text = fs::read_to_string(&run)
            .unwrap()
            .replace("adopted_edit_count\t1", "adopted_edit_count\t2")
            .replace("indel_candidate_count\t2", "indel_candidate_count\t3")
            .replace("indel_adopted_count\t0", "indel_adopted_count\t1");
        fs::write(run, text).unwrap();
        let summary = component.join("postrepair.tsv");
        let text = fs::read_to_string(&summary).unwrap().replace(
            &commands::resolve::sha256_hex(prior.as_bytes()),
            &commands::resolve::sha256_hex(current.as_bytes()),
        );
        fs::write(&summary, text).unwrap();
        args[5] = args[5].replacen("1:", "2:", 1);
        args[7] = args[7].replacen("1:", "2:", 1);
        commands::postrepair::run(&args).unwrap();
        assert_eq!(
            read_single_fasta(&commands::postrepair::selected_fasta(&component).unwrap())
                .unwrap()
                .1,
            current
        );
        let table = data.join("indel_candidates.tsv");
        let text = fs::read_to_string(&table).unwrap();
        fs::write(
            table,
            text.replace("\t10\t8\t8\t2\t0\tadopted", "\t10\t8\t7\t2\t0\tadopted"),
        )
        .unwrap();
        assert!(commands::postrepair::selected_fasta(&component).is_err());
        fs::remove_dir_all(dir).unwrap();
    }

    fn manual_current_provenance(
        config: &WorkflowConfig,
        component: &Path,
        resolved_sequence: &str,
        graph: &str,
    ) {
        let case = &config.cases[0];
        let source = polish_round_dir(case, 2);
        let first = polish_round_dir(case, 1);
        let resolved = config.config_path.parent().unwrap().join("resolved.fasta");
        write_with_parent(
            &resolved,
            &format!(">{}\n{resolved_sequence}\n", case.subgraph),
        )
        .unwrap();
        let checked = config.config_path.parent().unwrap().join("checked.gfa");
        fs::write(&checked, graph).unwrap();
        let details = config
            .config_path
            .parent()
            .unwrap()
            .join("resolve_details.tsv");
        fs::write(&details,format!("section\tid\tkey\tvalue\nrun\t.\tresolved_fasta\t{}\nrun\t.\tchecked_draft_gfa\t{}\nresolution\t{}\tselected_circular\ttrue\n",resolved.display(),checked.display(),case.subgraph)).unwrap();
        let first_fa = first.join("01.inputs/linear_subgraph.round_1.fasta");
        write_with_parent(
            &first_fa,
            &format!(">{}\n{resolved_sequence}\n", case.subgraph),
        )
        .unwrap();
        let source_fa = polished_aln_path(case, 2);
        for (round, number, fa) in [(&first, 1, &first_fa), (&source, 2, &source_fa)] {
            let old = fs::read_to_string(round.join("logs/report.tsv"))
                .unwrap_or_else(|_| "section\tname\tmetric\tvalue\n".into());
            let new = old
                .lines()
                .filter(|l| {
                    !l.starts_with("global\t")
                        && !l.starts_with("input\t")
                        && !l.starts_with("output\t")
                })
                .collect::<Vec<_>>()
                .join("\n")
                + "\n";
            write_with_parent(&round.join("logs/report.tsv"),&format!("{new}global\trun\torganelle\t{}\nglobal\trun\tsubgraph\t{}\nglobal\trun\tvalidate_round\t{number}\ninput\tdraft\tpath\t{}\ninput\tdraft\trecord_id\t{}\noutput\tsubgraph\tdir\t{}\noutput\tworkflow_round\tdir\t{}\noutput\tvalidation_fasta\tpath\t{}\n",case.organelle,case.subgraph,resolved.display(),case.subgraph,source.parent().unwrap().display(),round.display(),fa.display())).unwrap();
        }
        fs::copy(&source_fa, source.join("01.inputs/rotated_reference.fasta")).unwrap();
        write_with_parent(
            &component.join("round_hr/03.validate/01.data/circle_end/decision.tsv"),
            &format!(
                "metric\tvalue\nresolve_details\t{}\nfirst_polish_report\t{}\n",
                details.display(),
                first.join("logs/report.tsv").display()
            ),
        )
        .unwrap();
    }

    fn manual_current_remove_indels(config: &WorkflowConfig, component: &Path, review: bool) {
        let data = component.join("round_hr/03.validate/01.data");
        for name in ["indel_candidates.tsv", "indel_flank20_candidates.tsv"] {
            fs::remove_file(data.join(name)).unwrap();
        }
        fs::remove_dir_all(data.join("indel")).unwrap();
        let path = homorepair_outcome_path(&config.cases[0]);
        let mut text = fs::read_to_string(&path).unwrap();
        for key in ["indel_candidate_count", "indel_review_count"] {
            text = text.replace(&format!("{key}\t2"), &format!("{key}\t0"));
        }
        text = text.replace(
            "unresolved_count\t2",
            if review {
                "unresolved_count\t1"
            } else {
                "unresolved_count\t0"
            },
        );
        text = text.replace(
            "reason\tnonrepeat_indel_candidates_unconfirmed",
            if review {
                "reason\tcircle_end_candidates_unconfirmed"
            } else {
                "reason\tselected_edits_passed_ordinary_validation"
            },
        );
        if !review {
            text = text
                .replace("status\treview_required", "status\tcomplete")
                .replace(
                    "support_status\tcomplete",
                    "support_status\tsupport_verified",
                );
        }
        fs::write(&path, text).unwrap();
        let summary = component.join("postrepair.tsv");
        let text = fs::read_to_string(&summary)
            .unwrap()
            .replace("unresolved_count\t2", "unresolved_count\t1");
        fs::write(summary, text).unwrap();
    }

    #[test]
    fn manual_current_circle_replays_origin_native_source_and_failed_classes() {
        for fault in [
            "valid",
            "missing",
            "copies",
            "false_pass",
            "class",
            "origin",
            "disabled",
            "empty_stage",
            "reference",
            "reads",
        ] {
            let dir = test_dir(&format!("manual circle {fault}"));
            let (config, component, _) = manual_current_fixture(&dir, false);
            manual_current_remove_indels(&config, &component, true);
            let case = &config.cases[0];
            let source = polished_aln_path(case, 2);
            let hr = homorepair_out_dir(case);
            let data = hr.join("03.validate/01.data");
            let mut dna = read_single_fasta(&source).unwrap().1;
            dna.replace_range(63..78, &format!("G{}C", "T".repeat(13)));
            let fa = format!(">{}\n{dna}\n", case.subgraph);
            fs::write(&source, &fa).unwrap();
            fs::write(hr.join("01.inputs/linear_subgraph.round_hr.fasta"), &fa).unwrap();
            let p = component.join("postrepair.tsv");
            let old = fs::read_to_string(&p).unwrap();
            let hash = commands::resolve::sha256_hex(dna.as_bytes());
            let text = old
                .lines()
                .map(|l| {
                    if l.starts_with("source_sequence_sha256\t")
                        || l.starts_with("selected_sequence_sha256\t")
                    {
                        format!("{}\t{hash}", l.split('\t').next().unwrap())
                    } else {
                        l.into()
                    }
                })
                .collect::<Vec<_>>()
                .join("\n")
                + "\n";
            fs::write(&p, text).unwrap();
            manual_current_provenance(
                &config,
                &component,
                &dna,
                &format!("S\tn\t{dna}\nL\tn\t+\tn\t+\t0M\n"),
            );
            let cut = dna.len() / 2;
            let rotated = format!("{}{}", &dna[cut..], &dna[..cut]);
            let start = dna.len() - cut + 64;
            let end = start + 13;
            let circle = data.join("circle_end");
            fs::write(
                circle.join("rotated.fasta"),
                format!(">{}\n{rotated}\n", case.subgraph),
            )
            .unwrap();
            write_with_parent(
                &circle.join("rotated_validate/01.inputs/linear_subgraph.round_2.fasta"),
                &format!(">{}\n{rotated}\n", case.subgraph),
            )
            .unwrap();
            fs::copy(
                circle.join("rotated_validate/01.inputs/linear_subgraph.round_2.fasta"),
                circle.join("rotated_validate/01.inputs/rotated_reference.fasta"),
            )
            .unwrap();
            std::os::unix::fs::symlink(
                polish_round_dir(case, 2).join("01.inputs/subgraph_reads.fastq.gz"),
                circle.join("rotated_validate/01.inputs/subgraph_reads.fastq.gz"),
            )
            .unwrap();
            write_with_parent(&circle.join("rotated_validate/logs/report.tsv"),&format!("section\tname\tmetric\tvalue\nglobal\trun\torganelle\t{}\nglobal\trun\tsubgraph\t{}\nglobal\trun\tvalidate_round\t2\noutput\tvalidation_fasta\tpath\t{}\n",case.organelle,case.subgraph,circle.join("rotated_validate/01.inputs/linear_subgraph.round_2.fasta").display())).unwrap();
            write_with_parent(
                &circle.join("rotated_validate/03.validate/03.reports/sv_snv_indel_summary.tsv"),
                "metric\tvalue\nstatus\tpass\n",
            )
            .unwrap();
            write_with_parent(
                &circle.join("rotated_validate/03.validate/01.data/snv_indel_variants_combined.tsv"),
                &format!(
                    "pos\tref\talt\tcounts\tdepth\ttype\n{start}\tG{}\tG{}\t20\t40\thomopolymer\n",
                    "T".repeat(13),
                    "T".repeat(12)
                ),
            )
            .unwrap();
            let header="source_rows\tsource_length\tcut_0based\trotated_start_0based\trotated_end_0based\toriginal_start_0based\toriginal_end_0based\tmotif\treference_copies\tselected_copies\tdepth\ttotal_ALT\tmain_ALT\tremainder\tuncertain_reads\tstatus\treason\n";
            fs::write(circle.join("candidates.tsv"),format!("{header}1\t{}\t{cut}\t{start}\t{end}\t64\t77\tT\t13\t.\t40\t15\t15\t25\t1\treview_required\tfixed_read_pool_strict_gate_failed\n",dna.len())).unwrap();
            let mut calls = "source_rows\tread_id\tclass\texact_copies\n".to_string();
            for i in 0..40 {
                let (c, n) = if i < 15 {
                    ("ALT", "12")
                } else if i == 39 {
                    ("uncertain", ".")
                } else {
                    ("source", "13")
                };
                calls.push_str(&format!("1\tr{i}\t{c}\t{n}\n"));
            }
            fs::write(circle.join("read_classes.tsv"), calls).unwrap();
            let run = homorepair_outcome_path(case);
            let text=fs::read_to_string(&run).unwrap()+"outcome\thomorepair\tcircle_end_candidate_count\t1\noutcome\thomorepair\tcircle_end_adopted_count\t0\noutcome\thomorepair\tcircle_end_review_count\t1\noutcome\thomorepair\tcircle_end_reason\tfixed_read_pool_strict_gate_failed\n";
            fs::write(run, text).unwrap();
            let mut args = vec![
                "review".into(),
                "--result-dir".into(),
                component.display().to_string(),
                "--keep-current".into(),
                "--retain-circle".into(),
                "65:77:T:13".into(),
                "--reason".into(),
                "Current review retains this exact circle tract".into(),
            ];
            match fault {
                "disabled" | "empty_stage" => {
                    let p = homorepair_outcome_path(case);
                    let t = fs::read_to_string(&p)
                        .unwrap()
                        .replace(
                            "circle_end_candidate_count\t1",
                            "circle_end_candidate_count\t0",
                        )
                        .replace(
                            "circle_end_reason\tfixed_read_pool_strict_gate_failed",
                            if fault == "disabled" {
                                "circle_end_reason\tdisabled"
                            } else {
                                "circle_end_reason\tno_boundary_candidate"
                            },
                        );
                    fs::write(p, t).unwrap();
                    fs::write(circle.join("candidates.tsv"), header).unwrap();
                }
                "reference" => {
                    fs::write(
                        circle.join("rotated_validate/01.inputs/rotated_reference.fasta"),
                        ">wrong\nACGT\n",
                    )
                    .unwrap();
                }
                "reads" => {
                    let p = circle.join("rotated_validate/01.inputs/subgraph_reads.fastq.gz");
                    fs::remove_file(&p).unwrap();
                    fs::write(p, ">foreign\nACTG\n").unwrap();
                }
                "missing" => {
                    args.drain(4..6);
                }
                "copies" => args[5] = "65:77:T:12".into(),
                "false_pass" => {
                    let p = circle.join("candidates.tsv");
                    let t = fs::read_to_string(&p).unwrap();
                    fs::write(p, t.replace("\t15\t15\t25\t1\t", "\t30\t30\t10\t1\t")).unwrap();
                }
                "class" => {
                    let p = circle.join("read_classes.tsv");
                    let t = fs::read_to_string(&p).unwrap();
                    fs::write(p, t.replacen("\tALT\t12", "\tsource\t12", 1)).unwrap();
                }
                "origin" => {
                    let p = dir.join("resolve_details.tsv");
                    let t = fs::read_to_string(&p).unwrap();
                    fs::write(
                        p,
                        t.replace("selected_circular\ttrue", "selected_circular\tfalse"),
                    )
                    .unwrap();
                }
                _ => {}
            }
            let result = commands::postrepair::run(&args);
            if fault == "valid" {
                result.unwrap();
                assert_eq!(
                    read_single_fasta(&commands::postrepair::selected_fasta(&component).unwrap())
                        .unwrap()
                        .1,
                    dna
                );
            } else {
                assert!(result.is_err(), "{fault}");
                assert!(!component.join("manual_review.tsv").exists());
            }
            fs::remove_dir_all(dir).unwrap();
        }
    }

    #[test]
    fn manual_current_junction_preserves_other_supported_events_and_graph_direction() {
        for fault in [
            "valid",
            "current",
            "base",
            "direction",
            "count",
            "extra_review",
            "supported",
            "graph",
        ] {
            let dir = test_dir(&format!("manual junction {fault}"));
            let (config, component, _) = manual_current_fixture(&dir, true);
            manual_current_remove_indels(&config, &component, false);
            let case = &config.cases[0];
            let source = read_single_fasta(&polished_aln_path(case, 2)).unwrap().1;
            let left = &source[..598];
            let right = &source[600..];
            let zero = format!("{left}{right}");
            // The second candidate is a unique already-present one-base junction.
            let extra_left = &source[100..200];
            let extra_right = &source[201..301];
            let extra_base = &source[200..201];
            let graph=format!("S\ta\t{left}\nS\tb\t{right}\nS\tc\t{extra_left}\nS\td\t{extra_right}\nL\ta\t+\tb\t+\t0M\nL\tc\t+\td\t+\t0M\n");
            manual_current_provenance(&config, &component, &zero, &graph);
            let junction = component.join("junction_1bp");
            fs::create_dir(&junction).unwrap();
            let mut table =
                "left\tright\tflank\tzero\tA\tC\tG\tT\tselected_base\tfixed_pool\tstatus\n".to_string();
            for (l, r, b, status) in [
                ("a+", "b+", "A", "requires_graph_review"),
                ("c+", "d+", extra_base, "already_supported"),
            ] {
                for f in [20, 30, 40, 60, 80] {
                    let counts = ["A", "C", "G", "T"].map(|a| if a == b { "20" } else { "0" });
                    table.push_str(&format!(
                        "{l}\t{r}\t{f}\t0\t{}\t{b}\t20\t{status}\n",
                        counts.join("\t")
                    ));
                }
            }
            fs::write(junction.join("junction.tsv"), table).unwrap();
            let tokens = [("a+", "b+", "A"), ("c+", "d+", extra_base)]
                .iter()
                .flat_map(|(l, r, b)| [20, 30, 40, 60, 80].map(|f| format!("{l}>{r}:{f}:{b}")))
                .collect::<Vec<_>>()
                .join(",");
            let mut calls = "read_id\tmodel_hits\n".to_string();
            for i in 0..20 {
                calls.push_str(&format!("r{i}\t{tokens}\n"));
            }
            fs::write(junction.join("read_support.tsv"), calls).unwrap();
            let mut args = vec![
                "review".into(),
                "--result-dir".into(),
                component.display().to_string(),
                "--keep-current".into(),
                "--retain-junction".into(),
                "a+:b+:A:CA".into(),
                "--reason".into(),
                "Current review retains the current graph-direction CA, whose support was not measured"
                    .into(),
            ];
            match fault {
                "current" => args[5] = "a+:b+:A:AC".into(),
                "base" => args[5] = "a+:b+:C:CA".into(),
                "direction" => args[5] = "b-:a-:A:CA".into(),
                "count" => {
                    let p = junction.join("junction.tsv");
                    let t = fs::read_to_string(&p).unwrap();
                    fs::write(p, t.replacen("\t20\t0\t0\t0\tA", "\t21\t0\t0\t0\tA", 1)).unwrap();
                }
                "extra_review" => {
                    let p = junction.join("junction.tsv");
                    let t = fs::read_to_string(&p).unwrap();
                    fs::write(p, t.replace("already_supported", "requires_graph_review")).unwrap();
                }
                "supported" => {
                    let p = component.join("round_hr/03.validate/03.reports/polished.fasta");
                    let mut t = read_single_fasta(&p).unwrap().1;
                    t.replace_range(201..202, if extra_base == "A" { "C" } else { "A" });
                    fs::write(p, format!(">{}\n{t}\n", case.subgraph)).unwrap();
                }
                "graph" => {
                    let p = dir.join("checked.gfa");
                    let t = fs::read_to_string(&p).unwrap();
                    fs::write(p, t.replace("L\ta\t+\tb\t+\t0M", "L\ta\t+\tb\t+\t1M")).unwrap();
                }
                _ => {}
            }
            let result = commands::postrepair::run(&args);
            if fault == "valid" {
                result.unwrap();
                assert!(commands::postrepair::selected_fasta(&component).is_ok());
            } else {
                assert!(result.is_err(), "{fault}");
                assert!(!component.join("manual_review.tsv").exists());
            }
            fs::remove_dir_all(dir).unwrap();
        }
    }

    fn manual_correction_fixture(dir: &Path) -> (WorkflowConfig, PathBuf, Vec<String>) {
        let (config, component, _) = manual_keep_fixture(dir);
        let case = &config.cases[0];
        let round = polish_round_dir(case, 2);
        let fasta = polished_aln_path(case, 2);
        fs::write(&fasta, format!(">{}\nACGT\n", case.subgraph)).unwrap();
        let reads = round.join("01.inputs/subgraph_reads.fastq.gz");
        fs::write(
            &reads,
            &[
                31, 139, 8, 0, 0, 0, 0, 0, 2, 255, 179, 43, 74, 77, 76, 49, 228, 114, 116, 118, 15,
                1, 97, 46, 0, 160, 157, 223, 66, 16, 0, 0, 0,
            ],
        )
        .unwrap();
        let reference = round.join("01.inputs/rotated_reference.fasta");
        fs::write(&reference, ">reference_001\nACGT\n").unwrap();
        let report = polish_report_path(case, 2);
        let mut text = fs::read_to_string(&report).unwrap();
        // Record actual source-cohort and reference paths, not test defaults.
        writeln!(
            text,
            "input\treads\tpath\t{}\noutput\treference\tpath\t{}",
            reads.display(),
            reference.display()
        )
        .unwrap();
        fs::write(&report, text).unwrap();
        let evidence_dir = dir.join("evidence");
        fs::create_dir(&evidence_dir).unwrap();
        let mut attachments = Vec::new();
        for index in 1..=4 {
            let path = evidence_dir.join(format!("attachment_{index}.txt"));
            fs::write(
                &path,
                format!("Supplementary evidence {index}; no automatic approval is inferred\n"),
            )
            .unwrap();
            attachments.push(path);
        }
        let tools = dir.join("soft_paths.txt");
        fs::write(&tools, "fixture; no external process is run\n").unwrap();
        let mut args = vec![
            "--result-dir".into(),
            component.display().to_string(),
            "--correct-snv".into(),
            "1:A:C".into(),
            "--retain-indel".into(),
            "3:GT:GTT".into(),
            "--reason".into(),
            "Explicit current correction; supplementary evidence is not recomputed".into(),
            "--soft-paths".into(),
            tools.display().to_string(),
            "--threads".into(),
            "1".into(),
        ];
        for attachment in attachments {
            args.extend(["--evidence".into(), attachment.display().to_string()]);
        }
        (config, component, args)
    }

    fn manual_correction_native_fixture(
        manual: &Path,
        proposal: &Path,
        number: usize,
    ) -> Result<(), OrgraftError> {
        let component = manual.parent().unwrap();
        let source = PathBuf::from(
            read_metric_value(&component.join("postrepair.tsv"), "source_round")?.unwrap(),
        );
        let source_report = source.join("logs/report.tsv");
        let rows = homorepair_rows(&source_report)?;
        let identity = |metric: &str| {
            rows.iter()
                .find(|r| {
                    r.get("section").map(String::as_str) == Some("global")
                        && r.get("metric").map(String::as_str) == Some(metric)
                })
                .unwrap()["value"]
                .clone()
        };
        let organelle = identity("organelle");
        let name = identity("subgraph");
        let round = manual.join(format!("round_{number}"));
        let validated = round.join(format!("01.inputs/linear_subgraph.round_{number}.fasta"));
        write_with_parent(&validated, &fs::read_to_string(proposal)?)?;
        let reference = fs::canonicalize(source.join("01.inputs/rotated_reference.fasta"))?;
        let reads = fs::canonicalize(source.join("01.inputs/subgraph_reads.fastq.gz"))?;
        fs::copy(&reference, round.join("01.inputs/rotated_reference.fasta"))?;
        fs::copy(&reads, round.join("01.inputs/subgraph_reads.fastq.gz"))?;
        let tools = component.ancestors().nth(3).unwrap().join("soft_paths.txt");
        let report = format!("section\tname\tmetric\tvalue\nglobal\trun\torganelle\t{organelle}\nglobal\trun\tsubgraph\t{name}\nglobal\trun\tvalidate_round\t{number}\nglobal\trun\tthreads\t1\ninput\tdraft\tpath\t{}\ninput\treference\tpath\t{}\ninput\treads\tpath\t{}\ninput\tsoft_paths\tpath\t{}\noutput\tvalidation_fasta\tpath\t{}\nstage\tsv-eval:round_{number}\tstatus\tok\nstage\tsnv-indel-eval:round_{number}\tstatus\tok\nstage\tplot:round_{number}\tstatus\tok\n", proposal.display(), reference.display(), reads.display(), tools.display(), validated.display());
        write_with_parent(&round.join("logs/report.tsv"), &report)?;
        write_fake_sv_pass_summary(&round.join("03.validate/03.reports/sv_snv_indel_summary.tsv"))?;
        write_with_parent(
            &round.join("03.validate/03.reports/snv_indel_high.tsv"),
            "pos\tref\talt\tcounts\tdepth\ttype\n",
        )?;
        write_with_parent(
            &round.join("03.validate/01.data/snv_indel_variants_combined.tsv"),
            "pos\tref\talt\tcounts\tdepth\ttype\n3\tGT\tGTT\t7\t10\tInDel\n",
        )?;
        Ok(())
    }

    #[test]
    fn manual_correction_handoff_replays_the_declared_edit_and_preserves_automatic_review() {
        for attach in [false, true] {
            let dir = test_dir(if attach {
                "manual_correction_handoff_evidence"
            } else {
                "manual_correction_handoff_no_evidence"
            });
            let (config, component, mut args) = manual_correction_fixture(&dir);
            if !attach {
                args = args
                    .chunks(2)
                    .filter(|pair| pair[0] != "--evidence")
                    .flatten()
                    .cloned()
                    .collect();
            }
            let original_summary = fs::read(component.join("postrepair.tsv")).unwrap();
            let original_run = fs::read(component.join("round_hr/logs/run.tsv")).unwrap();
            assert!(commands::postrepair::selected_fasta(&component).is_err());
            commands::postrepair::manual_correction_fixture(
                &args,
                manual_correction_native_fixture,
            )
            .unwrap();
            let selected = commands::postrepair::selected_fasta(&component).unwrap();
            assert_eq!(read_single_fasta(&selected).unwrap().1, "CCGT");
            assert_eq!(
                selected,
                selected_repair_fasta(&config, &config.cases[0], 2).unwrap()
            );
            assert_eq!(
                commands::postrepair::selected_validation_dir(&component).unwrap(),
                fs::canonicalize(component.join("manual/round_3/03.validate")).unwrap()
            );
            assert_eq!(
                commands::postrepair::verify_ordinary_selection(
                    &component,
                    &polish_round_dir(&config.cases[0], 2),
                    0.2
                )
                .unwrap(),
                selected
            );
            assert!(commands::postrepair::verify_ordinary_selection(
                &component,
                &polish_round_dir(&config.cases[0], 1),
                0.2
            )
            .is_err());
            assert!(commands::postrepair::manual_review_label(&component)
                .unwrap()
                .contains("manual_corrected"));
            assert_eq!(
                fs::read(component.join("postrepair.tsv")).unwrap(),
                original_summary
            );
            assert_eq!(
                fs::read(component.join("round_hr/logs/run.tsv")).unwrap(),
                original_run
            );
            let receipt = fs::read_to_string(component.join("manual_review.tsv")).unwrap();
            assert!(receipt.contains("supplementary_evidence_recomputed\tfalse"));
            assert!(receipt.contains("decision\texplicit_manual_correction"));
            assert!(receipt.contains("retain_indel\t3:GT:GTT"));
            assert!(commands::postrepair::manual_correction_fixture(
                &args,
                manual_correction_native_fixture
            )
            .is_err());
            fs::remove_dir_all(dir).unwrap();
        }
    }

    #[test]
    fn manual_correction_rejects_bad_declarations_other_reviews_and_failed_validation() {
        for fault in [
            "wrong_ref",
            "overlap",
            "extra_edit",
            "wrong_retain",
            "duplicate_evidence",
            "cohort",
            "missing_evidence",
            "other_review",
            "continuation",
            "native_failed",
            "partial",
        ] {
            let dir = test_dir(&format!("manual_correction_reject_{fault}"));
            let (config, component, mut args) = manual_correction_fixture(&dir);
            let replace = |args: &mut Vec<String>, flag: &str, value: &str| {
                let index = args.iter().position(|x| x == flag).unwrap();
                args[index + 1] = value.into();
            };
            match fault {
                "wrong_ref" => replace(&mut args, "--correct-snv", "1:G:C"),
                "overlap" => replace(&mut args, "--correct-snv", "3:G:C"),
                "extra_edit" => replace(&mut args, "--correct-snv", "1:AC:CC"),
                "wrong_retain" => replace(&mut args, "--retain-indel", "3:GT:GTA"),
                "duplicate_evidence" => args.extend([
                    "--evidence".into(),
                    dir.join("evidence/attachment_1.txt").display().to_string(),
                ]),
                "cohort" => {
                    let report = polish_report_path(&config.cases[0], 2);
                    let reads = polish_round_dir(&config.cases[0], 2)
                        .join("01.inputs/subgraph_reads.fastq.gz");
                    let other = dir.join("other.fastq.gz");
                    fs::copy(&reads, &other).unwrap();
                    let text = fs::read_to_string(&report)
                        .unwrap()
                        .replace(&reads.display().to_string(), &other.display().to_string());
                    fs::write(report, text).unwrap();
                }
                "missing_evidence" => replace(
                    &mut args,
                    "--evidence",
                    &dir.join("missing.txt").display().to_string(),
                ),
                "other_review" => {
                    let run = component.join("round_hr/logs/run.tsv");
                    let text = fs::read_to_string(&run)
                        .unwrap()
                        .replace("\tunresolved_count\t1\n", "\tunresolved_count\t2\n");
                    fs::write(run, text).unwrap();
                }
                "continuation" => {
                    let path = component.join("postrepair.tsv");
                    let mut text = fs::read_to_string(&path).unwrap();
                    text.push_str("ordinary_pass_count\t2\nordinary_pass_limit\t3\n");
                    fs::write(path, text).unwrap();
                }
                "partial" => fs::create_dir(component.join("manual")).unwrap(),
                _ => {}
            }
            let result = if fault == "native_failed" {
                commands::postrepair::manual_correction_fixture(&args, |_, _, _| {
                    Err(OrgraftError::InvalidArgument(
                        "fixture native validation failed".into(),
                    ))
                })
            } else {
                commands::postrepair::manual_correction_fixture(
                    &args,
                    manual_correction_native_fixture,
                )
            };
            assert!(result.is_err(), "{fault}");
            assert!(!component.join("manual_review.tsv").exists(), "{fault}");
            fs::remove_dir_all(dir).unwrap();
        }
    }

    #[test]
    fn manual_correction_receipt_rejects_tampered_outputs_bindings_and_symbolic_paths() {
        for fault in [
            "candidate",
            "validated",
            "reference",
            "reads",
            "report",
            "high",
            "history",
            "source",
            "tools",
            "ordinary",
            "receipt",
            "symlink",
        ] {
            let dir = test_dir(&format!("manual_correction_changed_{fault}"));
            let (config, component, args) = manual_correction_fixture(&dir);
            commands::postrepair::manual_correction_fixture(
                &args,
                manual_correction_native_fixture,
            )
            .unwrap();
            let round = component.join("manual/round_3");
            let path = match fault {
                "candidate" => component.join("manual/candidate.fasta"),
                "validated" => round.join("01.inputs/linear_subgraph.round_3.fasta"),
                "reference" => round.join("01.inputs/rotated_reference.fasta"),
                "reads" => round.join("01.inputs/subgraph_reads.fastq.gz"),
                "report" => round.join("logs/report.tsv"),
                "high" => round.join("03.validate/03.reports/snv_indel_high.tsv"),
                "history" => dir.join("evidence/attachment_1.txt"),
                "source" => polished_aln_path(&config.cases[0], 2),
                "tools" => dir.join("soft_paths.txt"),
                "ordinary" => component.join("round_hr/03.validate/01.data/indel/read_classes.tsv"),
                _ => component.join("manual_review.tsv"),
            };
            if fault == "symlink" {
                #[cfg(unix)]
                {
                    let moved = component.join("saved_receipt.tsv");
                    fs::rename(&path, &moved).unwrap();
                    std::os::unix::fs::symlink(&moved, &path).unwrap();
                }
                #[cfg(not(unix))]
                {
                    fs::write(&path, "bad receipt\n").unwrap();
                }
            } else {
                let mut file = fs::OpenOptions::new().append(true).open(&path).unwrap();
                std::io::Write::write_all(&mut file, b"changed\n").unwrap();
            }
            assert!(
                commands::postrepair::selected_fasta(&component).is_err(),
                "{fault}"
            );
            assert!(
                commands::postrepair::selected_validation_dir(&component).is_err(),
                "{fault}"
            );
            fs::remove_dir_all(dir).unwrap();
        }
    }

    fn manual_keep_repeat_extra_fixture(dir: &Path) -> (WorkflowConfig, PathBuf, Vec<String>) {
        let (config, component, _) = manual_keep_fixture(dir);
        let case = &config.cases[0];
        let block = format!("{}{}{}", "ACG".repeat(10), "T".repeat(11), "GCA".repeat(10));
        let sequence = format!("{block}ACGACG{block}");
        let source = polished_aln_path(case, 2);
        fs::write(&source, format!(">source\n{sequence}\n")).unwrap();
        fs::write(homorepair_out_dir(case).join("01.inputs/linear_subgraph.round_hr.fasta"), format!(">snapshot\n{sequence}\n")).unwrap();
        homorepair_test_candidates(case, &[]);
        let run = homorepair_outcome_path(case);
        let text = fs::read_to_string(&run).unwrap().lines().filter(|line| !line.starts_with("probe\t")).collect::<Vec<_>>().join("\n") + "\n";
        let text = text.replace("native_probe_count\t1", "native_probe_count\t0")
            .replace("unresolved_count\t1", "unresolved_count\t2")
            .replace("nonrepeat_indel_candidates_unconfirmed", "supplemental_repeat_candidates_unconfirmed")
            .replace("indel_candidate_count\t1", "indel_candidate_count\t0")
            .replace("indel_review_count\t1", "indel_review_count\t0");
        fs::write(&run, format!("{text}outcome\thomorepair\trepeat_extra_candidate_count\t2\noutcome\thomorepair\trepeat_extra_adopted_count\t0\noutcome\thomorepair\trepeat_extra_review_count\t2\n")).unwrap();
        let data = homorepair_out_dir(case).join("03.validate/01.data");
        for name in ["indel_candidates.tsv", "indel_flank20_candidates.tsv"] { fs::remove_file(data.join(name)).unwrap(); }
        fs::remove_dir_all(data.join("indel")).unwrap();
        let mut table = "source_start\tsource_end\tmotif\toriginal_copies\tselected_copies\tentry\tstatus\treason\tdepth\ttotal_ALT\tmain_ALT\tremainder\tuntyped\tcopy_counts\n".to_string();
        let mut variants = "pos\tref\talt\tcounts\tdepth\ttype\n".to_string();
        for start in [31, 108] {
            writeln!(table, "{start}\t{}\tT\t11\t.\tlong_homopolymer_residual\treview_required\tsource_context_not_unique_on_selected\t0\t0\t0\t0\t0\t.", start + 10).unwrap();
            writeln!(variants, "{}\tG{}\tG{}#G{}\t2#3\t100\thomopolymer", start - 1, "T".repeat(11), "T".repeat(10), "T".repeat(12)).unwrap();
        }
        write_with_parent(&data.join("repeat_extra/candidates.tsv"), &table).unwrap();
        fs::write(polish_round_dir(case, 2).join("03.validate/01.data/snv_indel_variants_combined.tsv"), variants).unwrap();
        let summary = component.join("postrepair.tsv");
        let text = fs::read_to_string(&summary).unwrap().replace("unresolved_count\t1", "unresolved_count\t2")
            .replace(&commands::resolve::sha256_hex(b"ACGT"), &commands::resolve::sha256_hex(sequence.as_bytes()));
        fs::write(summary, text).unwrap();
        let args = vec!["review".into(), "--result-dir".into(), component.display().to_string(), "--keep-source".into(),
            "--repeat-extra".into(), "31:41:T:11".into(), "--repeat-extra".into(), "108:118:T:11".into(),
            "--reason".into(), "Current human review retains both original copies; the shared context cannot resolve individual copies".into()];
        (config, component, args)
    }

    #[test]
    fn manual_keep_repeat_extra_preserves_automatic_review_and_requires_the_exact_set() {
        let dir = test_dir("manual repeat-extra exact set");
        let (config, component, mut args) = manual_keep_repeat_extra_fixture(&dir);
        let case = &config.cases[0];
        let source = polished_aln_path(case, 2);
        let paths = [component.join("postrepair.tsv"), homorepair_outcome_path(case), source.clone(),
            homorepair_out_dir(case).join("03.validate/01.data/repeat_extra/candidates.tsv")];
        let before = paths.iter().map(|path| fs::read(path).unwrap()).collect::<Vec<_>>();
        assert!(commands::postrepair::selected_fasta(&component).is_err());
        assert!(selected_repair_fasta(&config, case, 2).is_err());
        args.swap(5, 7); // The CLI canonicalizes order, but cannot omit or duplicate a site.
        commands::postrepair::run(&args).unwrap();
        assert_eq!(commands::postrepair::selected_fasta(&component).unwrap(), source);
        assert_eq!(selected_repair_fasta(&config, case, 2).unwrap(), fs::canonicalize(&source).unwrap());
        assert!(verify_ordinary_postrepair(&component, &polish_round_dir(case, 2), 0.2).is_err());
        assert!(commands::postrepair::run(&args).is_err());
        let receipt = fs::read_to_string(component.join("manual_review.tsv")).unwrap();
        assert!(receipt.contains("event_type\trepeat_extra\n") && receipt.contains("event_count\t2\n"));
        assert!(receipt.contains("site_1\t31:41:T:11\n") && receipt.contains("site_2\t108:118:T:11\n"));
        assert!(receipt.contains("/repeat_extra/candidates.tsv\t") && receipt.contains("/snv_indel_variants_combined.tsv\t"));
        for (path, bytes) in paths.iter().zip(before) { assert_eq!(fs::read(path).unwrap(), bytes); }
        fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn manual_keep_repeat_extra_linear_boundary_preserves_unmeasured_evidence() {
        for (fault, left, right, accepted) in [
            ("right_boundary", 30, 15, true),
            ("left_boundary", 15, 30, true),
            ("complete_unique", 30, 30, false),
            ("measured", 30, 15, false),
            ("wrong_reason", 30, 15, false),
            ("adopted", 30, 15, false),
            ("snapshot", 30, 15, false),
        ] {
            let dir = test_dir(&format!("manual repeat-extra boundary {fault}"));
            let (config, component, _) = manual_keep_repeat_extra_fixture(&dir);
            let case = &config.cases[0];
            let source = polished_aln_path(case, 2);
            let (_, previous) = read_single_fasta(&source).unwrap();
            let sequence = format!(
                "{}{}{}",
                "ACG".repeat(left / 3),
                "T".repeat(11),
                "GCA".repeat(right / 3)
            );
            fs::write(&source, format!(">source\n{sequence}\n")).unwrap();
            let snapshot = homorepair_out_dir(case).join("01.inputs/linear_subgraph.round_hr.fasta");
            fs::write(&snapshot, format!(">snapshot\n{sequence}\n")).unwrap();
            let run = homorepair_outcome_path(case);
            let run_text = fs::read_to_string(&run)
                .unwrap()
                .replace("unresolved_count\t2", "unresolved_count\t1")
                .replace(
                    "repeat_extra_candidate_count\t2",
                    "repeat_extra_candidate_count\t1",
                )
                .replace(
                    "repeat_extra_review_count\t2",
                    "repeat_extra_review_count\t1",
                );
            fs::write(
                &run,
                if fault == "adopted" {
                    run_text.replace("adopted_edit_count\t0", "adopted_edit_count\t1")
                } else {
                    run_text
                },
            )
            .unwrap();
            let summary = component.join("postrepair.tsv");
            let text = fs::read_to_string(&summary)
                .unwrap()
                .replace("unresolved_count\t2", "unresolved_count\t1")
                .replace(
                    &commands::resolve::sha256_hex(previous.as_bytes()),
                    &commands::resolve::sha256_hex(sequence.as_bytes()),
                );
            fs::write(&summary, text).unwrap();
            let table =
                homorepair_out_dir(case).join("03.validate/01.data/repeat_extra/candidates.tsv");
            let header = fs::read_to_string(&table)
                .unwrap()
                .lines()
                .next()
                .unwrap()
                .to_string();
            let reason = if fault == "wrong_reason" {
                "fixed_read_pool_strict_gate_failed"
            } else {
                "source_context_not_unique_on_selected"
            };
            let depth = if fault == "measured" { 1 } else { 0 };
            fs::write(&table, format!("{header}\n{}\t{}\tT\t11\t.\tlong_homopolymer_residual\treview_required\t{reason}\t{depth}\t0\t0\t0\t0\t.\n", left + 1, left + 11)).unwrap();
            fs::write(
                polish_round_dir(case, 2).join("03.validate/01.data/snv_indel_variants_combined.tsv"),
                format!(
                    "pos\tref\talt\tcounts\tdepth\ttype\n{left}\tG{}\tG{}#G{}\t2#3\t100\thomopolymer\n",
                    "T".repeat(11),
                    "T".repeat(10),
                    "T".repeat(12)
                ),
            )
            .unwrap();
            if fault == "snapshot" {
                fs::write(&snapshot, format!(">snapshot\nA{sequence}\n")).unwrap();
            }
            let args = vec![
                "review".into(),
                "--result-dir".into(),
                component.display().to_string(),
                "--keep-source".into(),
                "--repeat-extra".into(),
                format!("{}:{}:T:11", left + 1, left + 11),
                "--reason".into(),
                "Current human review retains the original tract with unavailable linear context"
                    .into(),
            ];
            let originals =
                [&source, &snapshot, &summary, &run, &table].map(|path| fs::read(path).unwrap());
            assert_eq!(
                commands::postrepair::run(&args).is_ok(),
                accepted,
                "{fault}"
            );
            assert_eq!(
                component.join("manual_review.tsv").exists(),
                accepted,
                "{fault}"
            );
            if accepted {
                assert_eq!(
                    commands::postrepair::selected_fasta(&component).unwrap(),
                    source
                );
                assert_eq!(
                    selected_repair_fasta(&config, case, 2).unwrap(),
                    fs::canonicalize(&source).unwrap()
                );
                assert!(
                    verify_ordinary_postrepair(&component, &polish_round_dir(case, 2), 0.2).is_err()
                );
            }
            for (path, bytes) in [&source, &snapshot, &summary, &run, &table]
                .into_iter()
                .zip(originals)
            {
                assert_eq!(fs::read(path).unwrap(), bytes, "{fault}");
            }
            fs::remove_dir_all(dir).unwrap();
        }
    }

    #[test]
    fn manual_keep_repeat_extra_rejects_incomplete_sets_other_reviews_and_false_contexts() {
        for fault in ["missing", "additional", "duplicate", "mixed", "copies", "source", "snapshot", "source_row", "source_count",
            "ordinary_review", "indel_review", "circle_review", "adopted", "extra_adopted", "selected_alt", "read_gate", "pool", "count", "unique_context", "source_sv", "source_high"] {
            let dir = test_dir(&format!("manual repeat-extra reject {fault}"));
            let (config, component, mut args) = manual_keep_repeat_extra_fixture(&dir);
            let case = &config.cases[0];
            let data = homorepair_out_dir(case).join("03.validate/01.data");
            let run = homorepair_outcome_path(case);
            let source = polished_aln_path(case, 2);
            let table = data.join("repeat_extra/candidates.tsv");
            let variants = polish_round_dir(case, 2).join("03.validate/01.data/snv_indel_variants_combined.tsv");
            let change = |path: &Path, a: &str, b: &str| {
                let text = fs::read_to_string(path).unwrap(); assert!(text.contains(a)); fs::write(path, text.replace(a, b)).unwrap();
            };
            match fault {
                "missing" => { args.drain(6..8); },
                "additional" => args.extend(["--repeat-extra".into(), "200:210:T:11".into()]),
                "duplicate" => args.extend(["--repeat-extra".into(), "031:041:T:011".into()]),
                "mixed" => args.extend(["--pos".into(), "30".into(), "--ref".into(), "GT".into(), "--alt".into(), "GTT".into()]),
                "copies" => args[5] = "31:41:T:12".into(),
                "source" => change(&source, ">source\nA", ">source\nC"),
                "snapshot" => change(&homorepair_out_dir(case).join("01.inputs/linear_subgraph.round_hr.fasta"), ">snapshot\nA", ">snapshot\nC"),
                "source_row" => change(&variants, "107\t", "106\t"),
                "source_count" => change(&variants, "2#3\t100", "2#300\t100"),
                "ordinary_review" => { let mut text=fs::read_to_string(&run).unwrap(); text.push_str("outcome\thomorepair\tpending_count\t1\n"); fs::write(&run, text).unwrap(); },
                "indel_review" => change(&run, "indel_review_count\t0", "indel_review_count\t1"),
                "circle_review" => { let mut text=fs::read_to_string(&run).unwrap(); text.push_str("outcome\thomorepair\tcircle_end_review_count\t1\n"); fs::write(&run, text).unwrap(); },
                "adopted" => change(&run, "adopted_edit_count\t0", "adopted_edit_count\t1"),
                "extra_adopted" => change(&run, "repeat_extra_adopted_count\t0", "repeat_extra_adopted_count\t1"),
                "selected_alt" => change(&table, "\t11\t.\tlong_", "\t11\t10\tlong_"),
                "read_gate" => change(&table, "source_context_not_unique_on_selected", "fixed_read_pool_strict_gate_failed"),
                "pool" => change(&table, "\t0\t0\t0\t0\t0\t.", "\t1\t0\t0\t1\t1\t."),
                "count" => change(&run, "repeat_extra_review_count\t2", "repeat_extra_review_count\t1"),
                "unique_context" => {
                    let (_, original) = read_single_fasta(&source).unwrap();
                    let mut changed = original.clone(); changed.replace_range(77..78, "C");
                    fs::write(&source, format!(">source\n{changed}\n")).unwrap();
                    fs::write(homorepair_out_dir(case).join("01.inputs/linear_subgraph.round_hr.fasta"), format!(">snapshot\n{changed}\n")).unwrap();
                    change(&component.join("postrepair.tsv"), &commands::resolve::sha256_hex(original.as_bytes()), &commands::resolve::sha256_hex(changed.as_bytes()));
                },
                "source_sv" => change(&sv_summary_path(case,2), "status\tpass", "status\tfail"),
                "source_high" => fs::write(snv_indel_high_path(case,2), "pos\tref\talt\tcounts\tdepth\ttype\n30\tGT\tGTT\t90\t100\tInDel\n").unwrap(),
                _ => unreachable!(),
            }
            assert!(commands::postrepair::run(&args).is_err(), "{fault}");
            assert!(!component.join("manual_review.tsv").exists(), "{fault}");
            fs::remove_dir_all(dir).unwrap();
        }
    }

    #[test]
    fn manual_keep_repeat_extra_receipt_rejects_changed_missing_and_new_evidence() {
        for fault in ["source", "variants", "extra", "summary", "missing", "new_classes", "event", "event_count"] {
            let dir = test_dir(&format!("manual repeat-extra binding {fault}"));
            let (config, component, args) = manual_keep_repeat_extra_fixture(&dir);
            let case = &config.cases[0];
            let data = homorepair_out_dir(case).join("03.validate/01.data");
            commands::postrepair::run(&args).unwrap();
            let path = match fault {
                "source" => polished_aln_path(case,2),
                "variants" => polish_round_dir(case,2).join("03.validate/01.data/snv_indel_variants_combined.tsv"),
                "extra" | "missing" => data.join("repeat_extra/candidates.tsv"),
                "summary" => component.join("postrepair.tsv"),
                "new_classes" => data.join("repeat_extra/site_31_read_classes.tsv"),
                _ => component.join("manual_review.tsv"),
            };
            if fault == "missing" { fs::remove_file(path).unwrap(); }
            else if fault == "new_classes" { fs::write(path, "read_id\tcopy_count\nr1\t11\n").unwrap(); }
            else {
                let text = fs::read_to_string(&path).unwrap();
                let text = match fault {
                    "event" => text.replace("31:41:T:11", "32:42:T:11"),
                    "event_count" => text.replace("event_count\t2", "event_count\t1"),
                    _ => text + "\n",
                };
                fs::write(path, text).unwrap();
            }
            assert!(commands::postrepair::selected_fasta(&component).is_err(), "{fault}");
            assert!(selected_repair_fasta(&config,case,2).is_err(), "{fault}");
            fs::remove_dir_all(dir).unwrap();
        }
    }

    #[test]
    fn manual_keep_handoff_preserves_automatic_failure_and_replays_current_evidence() {
        let dir=test_dir("manual keep-source handoff");
        let (config,component,args)=manual_keep_fixture(&dir);
        let case=&config.cases[0];
        let source=polished_aln_path(case,2);
        let summary=fs::read(component.join("postrepair.tsv")).unwrap();
        let run=fs::read(homorepair_outcome_path(case)).unwrap();
        let source_bytes=fs::read(&source).unwrap();
        assert!(commands::postrepair::selected_fasta(&component).is_err());
        assert!(selected_repair_fasta(&config,case,2).is_err());
        commands::postrepair::run(&args).unwrap();
        let receipt=component.join("manual_review.tsv");
        let receipt_bytes=fs::read(&receipt).unwrap();
        assert!(String::from_utf8_lossy(&receipt_bytes).contains("review_kind\tcurrent_manual_review"));
        assert!(String::from_utf8_lossy(&receipt_bytes).contains("decision\tkeep_source"));
        assert!(String::from_utf8_lossy(&receipt_bytes).contains("/indel/flank20_read_classes.tsv\t"));
        assert_eq!(commands::postrepair::selected_fasta(&component).unwrap(),source);
        assert_eq!(selected_repair_fasta(&config,case,2).unwrap(),fs::canonicalize(&source).unwrap());
        assert!(postrepair_rebuild_args(&config,case,&source).unwrap().is_empty());
        assert_eq!(run_homorepair_after_checkpoint(&config,case,2).unwrap(),fs::canonicalize(&source).unwrap());
        assert!(verify_ordinary_postrepair(&component,&polish_round_dir(case,2),0.2).is_err());
        assert!(commands::postrepair::run(&args).is_err());
        assert_eq!(fs::read(&receipt).unwrap(),receipt_bytes);
        assert_eq!(fs::read(component.join("postrepair.tsv")).unwrap(),summary);
        assert_eq!(fs::read(homorepair_outcome_path(case)).unwrap(),run);
        assert_eq!(fs::read(source).unwrap(),source_bytes);
        let mut script=String::from("set -euo pipefail\n");
        write_homorepair_after_checkpoint(&mut script,&config,case,2).unwrap();
        let shell=dir.join("resume.sh");fs::write(&shell,&script).unwrap();
        assert!(script.contains("manual_review.tsv")&&script.contains("workflow postrepair-check"));
        assert!(!script.contains("eval"));
        assert!(Command::new("bash").arg("-n").arg(&shell).status().unwrap().success());
        let class=homorepair_out_dir(case).join("03.validate/01.data/indel/read_classes.tsv");
        fs::write(&class,fs::read_to_string(&class).unwrap().replace("r0","changed_id")).unwrap();
        assert!(commands::postrepair::selected_fasta(&component).unwrap_err().to_string().contains("evidence changed"));
        fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn manual_keep_rejects_wrong_event_other_reviews_and_inconsistent_pools_before_writing() {
        for fault in ["pos","ref","alt","reason","adopted","extra_review","circle_review","multiple","source","ordinary","class_count","pool_ids","nonfailed_pool","nested","default","source_row","source_sv","overlap"] {
            let dir=test_dir(&format!("manual keep negative {fault}"));
            let (config,component,mut args)=manual_keep_fixture(&dir);
            let case=&config.cases[0];
            let run=homorepair_outcome_path(case);
            let data=homorepair_out_dir(case).join("03.validate/01.data");
            let change=|path:&Path,a:&str,b:&str|{let text=fs::read_to_string(path).unwrap();assert!(text.contains(a));fs::write(path,text.replace(a,b)).unwrap();};
            match fault {
                "pos"=>args[5]="2".into(),"ref"=>args[7]="GA".into(),"alt"=>args[9]="GTTT".into(),
                "reason"=>args[11]="first\nsecond".into(),
                "adopted"=>change(&run,"adopted_edit_count\t0","adopted_edit_count\t1"),
                "extra_review"=>{let mut text=fs::read_to_string(&run).unwrap();text.push_str("outcome\thomorepair\trepeat_extra_review_count\t1\n");fs::write(&run,text).unwrap();},
                "circle_review"=>{let mut text=fs::read_to_string(&run).unwrap();text.push_str("outcome\thomorepair\tcircle_end_review_count\t1\n");fs::write(&run,text).unwrap();},
                "multiple"=>change(&run,"indel_candidate_count\t1","indel_candidate_count\t2"),
                "source"=>change(&polished_aln_path(case,2),"ACGT","TCGT"),
                "ordinary"=>{let mut row=homorepair_test_candidate(false);homorepair_test_set_candidate(&mut row,"main_ALT","1");homorepair_test_candidates(case,&[row]);},
                "class_count"=>change(&data.join("indel/read_classes.tsv"),"r0\tGTT\tALT","r0\t.\tuncertain"),
                "pool_ids"=>change(&data.join("indel/flank20_read_classes.tsv"),"r0","other"),
                "nonfailed_pool"=>change(&data.join("indel_candidates.tsv"),"\t4\t4\t6\t1\t","\t7\t7\t3\t1\t"),
                "nested"=>change(&data.join("indel_candidates.tsv"),"fixed_pool_gate\t0","fixed_pool_gate\t1"),
                "default"=>{args.remove(3);},
                "source_row"=>change(&polish_round_dir(case,2).join("03.validate/01.data/snv_indel_variants_combined.tsv"),"\t7\t10\t","\t6\t10\t"),
                "source_sv"=>change(&sv_summary_path(case,2),"status\tpass","status\tfail"),
                "overlap"=>{
                    let nested=component.join("round_2");fs::create_dir_all(&nested).unwrap();
                    change(&component.join("postrepair.tsv"),&format!("source_round\t{}",polish_round_dir(case,2).display()),&format!("source_round\t{}",nested.display()));
                },
                _=>unreachable!(),
            }
            let summary=fs::read(component.join("postrepair.tsv")).unwrap();
            let run_bytes=fs::read(&run).unwrap();
            assert!(commands::postrepair::run(&args).is_err(),"{fault}");
            assert!(!component.join("manual_review.tsv").exists(),"{fault}");
            assert_eq!(fs::read(component.join("postrepair.tsv")).unwrap(),summary);
            assert_eq!(fs::read(&run).unwrap(),run_bytes);
            fs::remove_dir_all(dir).unwrap();
        }
    }

    #[cfg(unix)]
    #[test]
    fn manual_keep_rejects_symbolic_receipts_components_and_moved_decisions() {
        let dir=test_dir("manual keep path guards");
        let (config,component,mut args)=manual_keep_fixture(&dir);
        let receipt=component.join("manual_review.tsv");
        let external=dir.join("external.tsv");fs::write(&external,"sentinel\n").unwrap();
        std::os::unix::fs::symlink(&external,&receipt).unwrap();
        assert!(commands::postrepair::run(&args).is_err());
        assert!(commands::postrepair::selected_fasta(&component).is_err());
        assert_eq!(fs::read_to_string(&external).unwrap(),"sentinel\n");
        fs::remove_file(&receipt).unwrap();
        let alias=dir.join("component_alias");std::os::unix::fs::symlink(&component,&alias).unwrap();
        args[2]=alias.display().to_string();assert!(commands::postrepair::run(&args).is_err());
        args[2]=component.display().to_string();commands::postrepair::run(&args).unwrap();
        let text=fs::read_to_string(&receipt).unwrap().replace("pos\t3\n","pos\t2\n");
        fs::write(&receipt,text).unwrap();
        assert!(selected_repair_fasta(&config,&config.cases[0],2).is_err());
        fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn homorepair_handoff_replays_a_nonrepeat_indel_and_checks_fixed_read_pool() {
        let dir = test_dir("homorepair nonrepeat indel handoff");
        let config = homorepair_test_config(&dir, true);
        let case = &config.cases[0];
        let input = polished_aln_path(case, 2);
        let hr = homorepair_out_dir(case);
        let data = hr.join("03.validate/01.data");
        let selected = hr.join("03.validate/03.reports/polished.fasta");
        let final_fasta = polished_aln_path(case, 3);
        write_with_parent(&input, ">source\nACGT\n").unwrap();
        write_with_parent(&hr.join("01.inputs/linear_subgraph.round_hr.fasta"), ">snapshot\nACGT\n").unwrap();
        write_with_parent(&selected, ">selected\nACGTT\n").unwrap();
        write_with_parent(&final_fasta, ">validated\nACGTT\n").unwrap();
        homorepair_test_candidates(case, &[]);
        let candidate_path = data.join("indel_candidates.tsv");
        write_with_parent(&candidate_path, concat!(
            "step\tsource_row\tsource_pos\tsource_ref\tsource_alt\tsource_support\tsource_depth\tsource_total_ALT\tprobe_pos\tprobe_delta\tdepth\ttotal_ALT\tmain_ALT\tremainder\tuncertain_reads\tstatus\treason\tnested_repeat_edits\tedit_pos\tedit_ref\tedit_alt\n",
            "1\t1\t3\tGT\tGTT\t7\t10\t7\t3\t1\t10\t7\t7\t3\t1\tadopted\tunique_main_ALT_above_fixed_pool_remainder\t0\t3\tGT\tGTT\n"
        )).unwrap();
        let mut classes = String::from("step\tread_id\texact_allele\tclass\n");
        for i in 0..7 { writeln!(classes, "1\tr{i}\tGTT\tALT").unwrap(); }
        for i in 7..9 { writeln!(classes, "1\tr{i}\tGT\tsource").unwrap(); }
        classes.push_str("1\tr9\t.\tuncertain\n");
        write_with_parent(&data.join("indel/read_classes.tsv"), &classes).unwrap();
        write_fake_sv_pass_summary(&sv_summary_path(case,3)).unwrap();
        write_with_parent(&snv_indel_high_path(case,3), "pos\tref\talt\ttotal_count\tdepth\tcounts\n").unwrap();
        write_with_parent(&polish_report_path(case,3), "section\tname\tmetric\tvalue\nstage\tsv-eval:round_3\tstatus\tok\nstage\tsnv-indel-eval:round_3\tstatus\tok\nstage\tplot:round_3\tstatus\tok\n").unwrap();
        let outcome = format!("section\tname\tmetric\tvalue\noutcome\thomorepair\tstatus\tcomplete\noutcome\thomorepair\tsupport_status\tsupport_verified\noutcome\thomorepair\tsource_round\t2\noutcome\thomorepair\tinput_fasta\t{}\noutcome\thomorepair\tselected_fasta\t{}\noutcome\thomorepair\tstage_count\t1\noutcome\thomorepair\tnative_probe_count\t0\noutcome\thomorepair\tlocal_batch_count\t0\noutcome\thomorepair\tadopted_edit_count\t1\noutcome\thomorepair\tunresolved_count\t0\noutcome\thomorepair\tindel_candidate_count\t1\noutcome\thomorepair\tindel_adopted_count\t1\noutcome\thomorepair\tindel_review_count\t0\noutcome\thomorepair\tfinal_validation_status\tpassed\noutcome\thomorepair\tvalidation_round\t3\noutcome\thomorepair\tfinal_validation_round\t3\noutcome\thomorepair\tfinal_validation_dir\t{}\noutcome\thomorepair\tfinal_validation_fasta\t{}\n",
            input.display(), selected.display(), polish_round_dir(case,3).display(), final_fasta.display());
        write_with_parent(&homorepair_outcome_path(case), &outcome).unwrap();
        // With no nested repeat edit, the replayed complete allele must still
        // be the original source ALT.
        assert_eq!(homorepair_selected_fasta(&config, case, 2).unwrap(), final_fasta);
        let original_candidates = fs::read_to_string(&candidate_path).unwrap();
        let paired_header = concat!(
            "step\tsource_row\tsource_rows\tsource_pos\tsource_ref\tsource_alt\tsource_support\tsource_depth\tsource_total_ALT\tprobe_pos\tprobe_delta\tdepth\ttotal_ALT\tmain_ALT\tremainder\tuncertain_reads\tstatus\treason\tnested_repeat_edits\tedit_pos\tedit_ref\tedit_alt\ttrial_candidates\ttested_probes\n"
        );
        let sixty_failed = format!("{paired_header}1\t1\t1\t3\tGT\tGTT\t7\t10\t7\t3\t1\t10\t4\t4\t6\t1\treview_required\thighest_source_ALT_failed_fixed_pool_gate\t0\t.\t.\t.\t1\t1\n");
        let twenty_adopted = format!("{paired_header}1\t1\t1\t3\tGT\tGTT\t7\t10\t7\t3\t1\t10\t7\t7\t3\t1\tadopted\tunique_main_ALT_above_fixed_20bp_pool_remainder\t0\t3\tGT\tGTT\t1\t1\n");
        let short_path = data.join("indel_flank20_candidates.tsv");
        let short_classes = data.join("indel/flank20_read_classes.tsv");
        fs::write(&candidate_path, &sixty_failed).unwrap();
        write_with_parent(&short_path, &twenty_adopted).unwrap();
        write_with_parent(&short_classes, &classes).unwrap();
        assert_eq!(homorepair_selected_fasta(&config, case, 2).unwrap(), final_fasta);
        fs::write(&short_path, twenty_adopted.replace("\t10\t7\t7\t3\t1\tadopted", "\t11\t7\t7\t4\t1\tadopted")).unwrap();
        assert!(homorepair_selected_fasta(&config, case, 2).is_err());
        fs::write(&short_path, &twenty_adopted).unwrap();
        fs::remove_file(&short_classes).unwrap();
        assert!(homorepair_selected_fasta(&config, case, 2).is_err());
        fs::remove_file(&short_path).unwrap();
        fs::write(&candidate_path, &original_candidates).unwrap();
        let changed_source_alt = original_candidates.replace(
            "\tGT\tGTT\t7\t10\t7\t3\t1\t10",
            "\tGT\tGTTT\t7\t10\t7\t3\t2\t10",
        );
        fs::write(&candidate_path, &changed_source_alt).unwrap();
        assert!(homorepair_selected_fasta(&config, case, 2).is_err());
        // The source ALT can differ only when one nested repeat edit produced
        // the exact complete allele recorded by the read classes and replay.
        let nested_candidates = changed_source_alt.replace(
            "unique_main_ALT_above_fixed_pool_remainder\t0\t3\tGT\tGTT",
            "unique_main_ALT_above_fixed_pool_remainder\t1\t3\tGT\tGTT",
        );
        fs::write(&candidate_path, &nested_candidates).unwrap();
        assert_eq!(homorepair_selected_fasta(&config, case, 2).unwrap(), final_fasta);
        let offsite_edit = nested_candidates.replacen("\t3\tGT\tGTT\n", "\t2\tCGT\tTGTT\n", 1);
        assert_ne!(offsite_edit, nested_candidates);
        fs::write(&candidate_path, offsite_edit).unwrap();
        fs::write(&selected, ">offsite\nATGTT\n").unwrap();
        fs::write(&final_fasta, ">offsite\nATGTT\n").unwrap();
        assert!(homorepair_selected_fasta(&config, case, 2).is_err());
        fs::write(&candidate_path, &nested_candidates).unwrap();
        fs::write(&selected, ">selected\nACGTT\n").unwrap();
        fs::write(&final_fasta, ">validated\nACGTT\n").unwrap();
        let class_path = data.join("indel/read_classes.tsv");
        fs::write(&class_path, classes.replacen("GTT\tALT", "GT\tALT", 1)).unwrap();
        assert!(homorepair_selected_fasta(&config, case, 2).is_err());
        fs::write(&class_path, &classes).unwrap();
        fs::write(&selected, ">wrong\nACGTA\n").unwrap();
        fs::write(&final_fasta, ">wrong\nACGTA\n").unwrap();
        assert!(homorepair_selected_fasta(&config, case, 2).is_err());
        fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn short_complete_allele_shape_keeps_mixed_length_competitors() {
        let winner = "TTTTCAAACCATAAAAA";
        assert!(short_complete_allele_shape("T", winner, true));
        for other in ["TTA", "TTTA", "TTTTCAAACC", "TCACCCATAAAA", "TTCAAACCTAAA"] {
            assert!(short_complete_allele_shape("T", other, false));
            assert!(!short_complete_allele_shape("T", other, true));
        }
        assert!(!short_complete_allele_shape("T", "A", false));
        assert!(!short_complete_allele_shape("T", "AAAA", false));
        assert!(!short_complete_allele_shape("TT", "T", false));
    }

    #[test]
    fn homorepair_handoff_checks_short_complete_allele_and_original_read_pool() {
        let dir = test_dir("homorepair short complete allele handoff");
        let config = homorepair_test_config(&dir, true);
        let case = &config.cases[0];
        let input = polished_aln_path(case, 2);
        let hr = homorepair_out_dir(case);
        let data = hr.join("03.validate/01.data");
        let selected = hr.join("03.validate/03.reports/polished.fasta");
        let final_fasta = polished_aln_path(case, 3);
        let source_alt = "GTACGTTAGCTGCA";
        let selected_sequence = format!("AC{source_alt}");
        write_with_parent(&input, ">source\nACGT\n").unwrap();
        write_with_parent(&hr.join("01.inputs/linear_subgraph.round_hr.fasta"), ">snapshot\nACGT\n").unwrap();
        write_with_parent(&selected, &format!(">selected\n{selected_sequence}\n")).unwrap();
        write_with_parent(&final_fasta, &format!(">validated\n{selected_sequence}\n")).unwrap();
        homorepair_test_candidates(case, &[]);
        let header = "step\tsource_row\tsource_rows\tsource_pos\tsource_ref\tsource_alt\tsource_support\tsource_depth\tsource_total_ALT\tprobe_pos\tprobe_delta\tdepth\ttotal_ALT\tmain_ALT\tremainder\tuncertain_reads\tstatus\treason\tnested_repeat_edits\tedit_pos\tedit_ref\tedit_alt\ttrial_candidates\ttested_probes\n";
        let common = format!("1\t1\t1\t3\tGT\t{source_alt}\t7\t10\t7\t3\t12");
        let sixty = format!("{header}{common}\t10\t4\t4\t6\t1\treview_required\thighest_source_ALT_failed_fixed_pool_gate\t0\t.\t.\t.\t1\t1\n");
        let twenty = format!("{header}{common}\t10\t4\t4\t6\t1\treview_required\t20bp_fixed_pool_gate_failed\t0\t.\t.\t.\t1\t1\n");
        let short = format!("{header}{common}\t10\t7\t7\t3\t1\tadopted\tunique_main_ALT_above_fixed_pool_remainder\t0\t3\tGT\t{source_alt}\t1\t1\n");
        let sixty_path = data.join("indel_candidates.tsv");
        let twenty_path = data.join("indel_flank20_candidates.tsv");
        let short_path = data.join("indel_short_candidates.tsv");
        write_with_parent(&sixty_path, &sixty).unwrap();
        write_with_parent(&twenty_path, &twenty).unwrap();
        write_with_parent(&short_path, &short).unwrap();
        let class_header = "step\tread_id\texact_allele\tclass\n";
        let mut old_classes = class_header.to_string();
        let mut new_classes = class_header.to_string();
        for i in 0..10 {
            let (old_allele, old_class) = if i < 4 { (source_alt, "ALT") } else if i < 9 { ("GT", "source") } else { (".", "uncertain") };
            let (new_allele, new_class) = if i < 7 { (source_alt, "ALT") } else if i < 9 { ("GT", "source") } else { (".", "uncertain") };
            writeln!(old_classes, "1\tr{i}\t{old_allele}\t{old_class}").unwrap();
            writeln!(new_classes, "1\tr{i}\t{new_allele}\t{new_class}").unwrap();
        }
        write_with_parent(&data.join("indel/read_classes.tsv"), &old_classes).unwrap();
        write_with_parent(&data.join("indel/flank20_read_classes.tsv"), &old_classes).unwrap();
        let short_classes_path = data.join("indel/short_read_classes.tsv");
        write_with_parent(&short_classes_path, &new_classes).unwrap();
        write_fake_sv_pass_summary(&sv_summary_path(case,3)).unwrap();
        write_with_parent(&snv_indel_high_path(case,3), "pos\tref\talt\ttotal_count\tdepth\tcounts\n").unwrap();
        write_with_parent(&polish_report_path(case,3), "section\tname\tmetric\tvalue\nstage\tsv-eval:round_3\tstatus\tok\nstage\tsnv-indel-eval:round_3\tstatus\tok\nstage\tplot:round_3\tstatus\tok\n").unwrap();
        let outcome = format!("section\tname\tmetric\tvalue\noutcome\thomorepair\tstatus\tcomplete\noutcome\thomorepair\tsupport_status\tsupport_verified\noutcome\thomorepair\tsource_round\t2\noutcome\thomorepair\tinput_fasta\t{}\noutcome\thomorepair\tselected_fasta\t{}\noutcome\thomorepair\tstage_count\t1\noutcome\thomorepair\tnative_probe_count\t0\noutcome\thomorepair\tlocal_batch_count\t0\noutcome\thomorepair\tadopted_edit_count\t1\noutcome\thomorepair\tunresolved_count\t0\noutcome\thomorepair\tindel_candidate_count\t1\noutcome\thomorepair\tindel_adopted_count\t1\noutcome\thomorepair\tindel_review_count\t0\noutcome\thomorepair\tfinal_validation_status\tpassed\noutcome\thomorepair\tvalidation_round\t3\noutcome\thomorepair\tfinal_validation_round\t3\noutcome\thomorepair\tfinal_validation_dir\t{}\noutcome\thomorepair\tfinal_validation_fasta\t{}\n", input.display(), selected.display(), polish_round_dir(case,3).display(), final_fasta.display());
        write_with_parent(&homorepair_outcome_path(case), &outcome).unwrap();
        assert_eq!(homorepair_selected_fasta(&config, case, 2).unwrap(), final_fasta);
        write_with_parent(&short_path, &short.replacen("1\t1\t1\t3", "2\t1\t1\t3", 1)).unwrap();
        assert!(homorepair_selected_fasta(&config, case, 2).is_err());
        write_with_parent(&short_path, &short).unwrap();
        write_with_parent(&short_classes_path, &new_classes.replacen("1\tr0\t", "1\tnew_read\t", 1)).unwrap();
        assert!(homorepair_selected_fasta(&config, case, 2).is_err());
        write_with_parent(&short_classes_path, &new_classes).unwrap();
        write_with_parent(&short_path, &short.replacen("\t10\t7\t7\t3\t1\tadopted", "\t10\t7\t6\t3\t1\tadopted", 1)).unwrap();
        assert!(homorepair_selected_fasta(&config, case, 2).is_err());
        write_with_parent(&short_path, &short).unwrap();
        write_with_parent(&selected, ">wrong\nACGTACGTTAGCTGCT\n").unwrap();
        write_with_parent(&final_fasta, ">wrong\nACGTACGTTAGCTGCT\n").unwrap();
        assert!(homorepair_selected_fasta(&config, case, 2).is_err());
        fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn long_event_shape_requires_a_pure_indel_and_complete_length() {
        let inserted = format!("G{}T", "ACGT".repeat(25));
        assert!(long_indel_event_shape("GT", &inserted));
        assert!(long_indel_event_shape(&inserted, "GT"));
        assert!(!long_indel_event_shape("GT", &format!("G{}T", "A".repeat(99))));
        assert!(!long_indel_event_shape("GT", &format!("C{}T", "A".repeat(100))));
        assert!(!long_indel_event_shape("", &inserted));
    }

    #[test]
    fn homorepair_handoff_checks_long_event_union_and_separate_competing_alleles() {
        let dir = test_dir("homorepair long event handoff");
        let config = homorepair_test_config(&dir, true);
        let case = &config.cases[0];
        let input = polished_aln_path(case, 2);
        let hr = homorepair_out_dir(case);
        let data = hr.join("03.validate/01.data");
        let selected = hr.join("03.validate/03.reports/polished.fasta");
        let final_fasta = polished_aln_path(case, 3);
        let core = "ACGT".repeat(25);
        let source_alt = format!("GT{}", &core[..99]);
        let complete_alt = format!("GT{core}");
        let competitor = format!("GT{}", "C".repeat(100));
        let selected_sequence = format!("AC{complete_alt}");
        write_with_parent(&input, ">source\nACGT\n").unwrap();
        write_with_parent(&hr.join("01.inputs/linear_subgraph.round_hr.fasta"), ">snapshot\nACGT\n").unwrap();
        write_with_parent(&selected, &format!(">selected\n{selected_sequence}\n")).unwrap();
        write_with_parent(&final_fasta, &format!(">validated\n{selected_sequence}\n")).unwrap();
        homorepair_test_candidates(case, &[]);
        let header = "step\tsource_row\tsource_rows\tsource_pos\tsource_ref\tsource_alt\tsource_support\tsource_depth\tsource_total_ALT\tprobe_pos\tprobe_delta\tdepth\ttotal_ALT\tmain_ALT\tremainder\tuncertain_reads\tstatus\treason\tnested_repeat_edits\tedit_pos\tedit_ref\tedit_alt\ttrial_candidates\ttested_probes\tentry_event_support\n";
        let common = format!("1\t1\t1\t3\tGT\t{source_alt}\t4\t10\t7\t3\t99");
        let sixty = format!("{header}{common}\t10\t6\t4\t4\t1\treview_required\thighest_source_ALT_failed_fixed_pool_gate\t0\t.\t.\t.\t1\t1\t7\n");
        let twenty = format!("{header}{common}\t10\t6\t4\t4\t1\treview_required\t20bp_fixed_pool_gate_failed\t0\t.\t.\t.\t1\t1\t7\n");
        let event = format!("{header}{common}\t12\t8\t6\t4\t1\tadopted\tunique_main_event_above_fixed_pool_remainder\t1\t3\tGT\t{complete_alt}\t1\t1\t7\n");
        let sixty_path = data.join("indel_candidates.tsv");
        let twenty_path = data.join("indel_flank20_candidates.tsv");
        let event_path = data.join("indel_event_candidates.tsv");
        write_with_parent(&sixty_path, &sixty).unwrap();
        write_with_parent(&twenty_path, &twenty).unwrap();
        write_with_parent(&event_path, &event).unwrap();
        let mut exact = "step\tread_id\texact_allele\tclass\n".to_string();
        let mut event_classes = "step\tread_id\tevent_class\n".to_string();
        for i in 0..12 {
            if i < 10 {
                let (allele, label) = if i < 4 { (complete_alt.as_str(), "ALT") }
                    else if i < 6 { (competitor.as_str(), "ALT") }
                    else if i < 9 { ("GT", "source") } else { (".", "uncertain") };
                writeln!(exact, "1\tr{i}\t{allele}\t{label}").unwrap();
            }
            let label = if i < 4 || i >= 10 { "selected_ALT_event" }
                else if i < 6 { "other_exact_ALT" }
                else if i < 9 { "source" } else { "uncertain" };
            writeln!(event_classes, "1\tr{i}\t{label}").unwrap();
        }
        let exact_path = data.join("indel/read_classes.tsv");
        let classes_path = data.join("indel/event_read_classes.tsv");
        write_with_parent(&exact_path, &exact).unwrap();
        write_with_parent(&data.join("indel/flank20_read_classes.tsv"), &exact).unwrap();
        write_with_parent(&classes_path, &event_classes).unwrap();
        write_fake_sv_pass_summary(&sv_summary_path(case,3)).unwrap();
        write_with_parent(&snv_indel_high_path(case,3), "pos\tref\talt\ttotal_count\tdepth\tcounts\n").unwrap();
        write_with_parent(&polish_report_path(case,3), "section\tname\tmetric\tvalue\nstage\tsv-eval:round_3\tstatus\tok\nstage\tsnv-indel-eval:round_3\tstatus\tok\nstage\tplot:round_3\tstatus\tok\n").unwrap();
        let outcome = format!("section\tname\tmetric\tvalue\noutcome\thomorepair\tstatus\tcomplete\noutcome\thomorepair\tsupport_status\tsupport_verified\noutcome\thomorepair\tsource_round\t2\noutcome\thomorepair\tinput_fasta\t{}\noutcome\thomorepair\tselected_fasta\t{}\noutcome\thomorepair\tstage_count\t1\noutcome\thomorepair\tnative_probe_count\t0\noutcome\thomorepair\tlocal_batch_count\t0\noutcome\thomorepair\tadopted_edit_count\t1\noutcome\thomorepair\tunresolved_count\t0\noutcome\thomorepair\tindel_candidate_count\t1\noutcome\thomorepair\tindel_adopted_count\t1\noutcome\thomorepair\tindel_review_count\t0\noutcome\thomorepair\tfinal_validation_status\tpassed\noutcome\thomorepair\tvalidation_round\t3\noutcome\thomorepair\tfinal_validation_round\t3\noutcome\thomorepair\tfinal_validation_dir\t{}\noutcome\thomorepair\tfinal_validation_fasta\t{}\n", input.display(), selected.display(), polish_round_dir(case,3).display(), final_fasta.display());
        write_with_parent(&homorepair_outcome_path(case), &outcome).unwrap();
        assert_eq!(homorepair_selected_fasta(&config, case, 2).unwrap(), final_fasta);
        // The long frame may extend the read pool, but must retain every old ID.
        for faulty in [event_classes.replacen("1\tr0\t", "1\tnew_read\t", 1),
            event_classes.replacen("1\tr10\t", "1\tr0\t", 1),
            event_classes.replace("other_exact_ALT", "selected_ALT_event"),
            event_classes.replace("1\tr10\tselected_ALT_event", "1\tr10\tother_exact_ALT")] {
            write_with_parent(&classes_path, &faulty).unwrap();
            assert!(homorepair_selected_fasta(&config, case, 2).is_err());
        }
        write_with_parent(&classes_path, &event_classes).unwrap();
        for faulty in [event.replace("\t12\t8\t6\t4\t1\tadopted", "\t12\t8\t5\t4\t1\tadopted"),
            event.replace("\t1\t1\t7\n", "\t1\t1\t11\n"),
            event.replacen("\t4\t10\t7\t3", "\t5\t10\t7\t3", 1)] {
            write_with_parent(&event_path, &faulty).unwrap();
            assert!(homorepair_selected_fasta(&config, case, 2).is_err());
        }
        write_with_parent(&event_path, &event).unwrap();
        // An exact call in the selected event family is not a separate ALT vote.
        write_with_parent(&exact_path, &exact.replace(&competitor, &complete_alt)).unwrap();
        assert!(homorepair_selected_fasta(&config, case, 2).is_err());
        write_with_parent(&exact_path, &exact).unwrap();
        assert_eq!(homorepair_selected_fasta(&config, case, 2).unwrap(), final_fasta);
        fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn homorepair_handoff_replays_junction_without_erasing_failed_event_or_other_review() {
        let dir = test_dir("homorepair exact junction handoff");
        let mut config = homorepair_test_config(&dir, true);
        config.cases[0].reads = dir.join("reads.fastq");
        let case = &config.cases[0];
        let input = polished_aln_path(case, 2);
        let hr = homorepair_out_dir(case);
        let data = hr.join("03.validate/01.data");
        let variants =
            polish_round_dir(case, 2).join("03.validate/01.data/snv_indel_variants_combined.tsv");
        let fixture = commands::postrepair::junction::test_fixture(&data, &variants, &case.reads);
        let source_reads = polish_round_dir(case, 2).join("01.inputs/subgraph_reads.fastq.gz");
        fs::create_dir_all(source_reads.parent().unwrap()).unwrap();
        let compressed = Command::new("gzip")
            .arg("-c")
            .arg(&case.reads)
            .output()
            .unwrap();
        assert!(compressed.status.success());
        fs::write(&source_reads, &compressed.stdout).unwrap();

        let selected = hr.join("03.validate/03.reports/polished.fasta");
        let final_fasta = polished_aln_path(case, 3);
        write_with_parent(&input, &format!(">component\n{}\n", fixture.source)).unwrap();
        write_with_parent(
            &hr.join("01.inputs/linear_subgraph.round_hr.fasta"),
            &format!(">component\n{}\n", fixture.source),
        )
        .unwrap();
        for path in [&selected, &final_fasta] {
            write_with_parent(path, &format!(">component\n{}\n", fixture.candidate)).unwrap();
        }
        homorepair_test_candidates(case, &[]);
        let sixty = fixture.gate.replace(
            "paired_event_strict_gate_failed",
            "highest_source_ALT_failed_fixed_pool_gate",
        );
        let twenty = fixture.gate.replace(
            "paired_event_strict_gate_failed",
            "20bp_fixed_pool_gate_failed",
        );
        write_with_parent(&data.join("indel_candidates.tsv"), &sixty).unwrap();
        write_with_parent(&data.join("indel_flank20_candidates.tsv"), &twenty).unwrap();
        let junction = fixture.gate.replace("\t5\t2\t2\t3\t3\treview_required\tpaired_event_strict_gate_failed\t0\t.\t.\t.",
            &format!("\t5\t3\t3\t2\t2\tadopted\tunique_main_junction_above_fixed_pool_remainder\t0\t{}\t{}\t{}",fixture.site,fixture.reference,fixture.alternate));
        write_with_parent(&data.join("indel_junction_candidates.tsv"), &junction).unwrap();
        let mut exact = String::from("step\tread_id\texact_allele\tclass\n");
        for i in 0..5 {
            let (allele, label) = if i < 2 {
                (fixture.alternate.as_str(), "ALT")
            } else {
                (".", "uncertain")
            };
            writeln!(exact, "1\tr{i}\t{allele}\t{label}").unwrap();
        }
        write_with_parent(&data.join("indel/read_classes.tsv"), &exact).unwrap();
        write_with_parent(&data.join("indel/flank20_read_classes.tsv"), &exact).unwrap();
        write_fake_sv_pass_summary(&sv_summary_path(case, 3)).unwrap();
        write_with_parent(
            &snv_indel_high_path(case, 3),
            "pos\tref\talt\ttotal_count\tdepth\tcounts\n",
        )
        .unwrap();
        write_with_parent(&polish_report_path(case,3),"section\tname\tmetric\tvalue\nstage\tsv-eval:round_3\tstatus\tok\nstage\tsnv-indel-eval:round_3\tstatus\tok\nstage\tplot:round_3\tstatus\tok\n").unwrap();
        let outcome = format!("section\tname\tmetric\tvalue\noutcome\thomorepair\tstatus\tcomplete\noutcome\thomorepair\tsupport_status\tsupport_verified\noutcome\thomorepair\tsource_round\t2\noutcome\thomorepair\tinput_fasta\t{}\noutcome\thomorepair\tselected_fasta\t{}\noutcome\thomorepair\tstage_count\t1\noutcome\thomorepair\tnative_probe_count\t0\noutcome\thomorepair\tlocal_batch_count\t0\noutcome\thomorepair\tadopted_edit_count\t1\noutcome\thomorepair\tunresolved_count\t0\noutcome\thomorepair\tindel_candidate_count\t1\noutcome\thomorepair\tindel_adopted_count\t1\noutcome\thomorepair\tindel_review_count\t0\noutcome\thomorepair\tfinal_validation_status\tpassed\noutcome\thomorepair\tvalidation_round\t3\noutcome\thomorepair\tfinal_validation_round\t3\noutcome\thomorepair\tfinal_validation_dir\t{}\noutcome\thomorepair\tfinal_validation_fasta\t{}\n",input.display(),selected.display(),polish_round_dir(case,3).display(),final_fasta.display());
        write_with_parent(&homorepair_outcome_path(case), &outcome).unwrap();
        assert_eq!(
            homorepair_selected_fasta(&config, case, 2).unwrap(),
            final_fasta
        );
        // The public standalone handoff synthesizes a case without configured
        // reads. It must use the actual native source round, including copies.
        let report = polish_report_path(case, 3);
        let mut native_report = fs::read_to_string(&report).unwrap();
        native_report.push_str(&format!("global\trun\torganelle\t{}\nglobal\trun\tsubgraph\t{}\noutput\tvalidation_fasta\tpath\t{}\n",
            case.organelle,case.subgraph,final_fasta.display()));
        fs::write(&report, native_report).unwrap();
        let component = hr.parent().unwrap();
        let summary = format!("metric\tvalue\nstatus\tcomplete\nunresolved_count\t0\nselection_method\tordinary\nsource_round\t{}\nsource_fasta\t{}\nsource_sequence_sha256\t{}\nselected_fasta\t{}\nselected_sequence_sha256\t{}\nfinal_validation_dir\t{}\nmin_alt_fraction\t0.2\n",
            polish_round_dir(case,2).display(),input.display(),commands::resolve::sha256_hex(fixture.source.as_bytes()),
            final_fasta.display(),commands::resolve::sha256_hex(fixture.candidate.as_bytes()),
            polish_round_dir(case,3).join("03.validate").display());
        fs::write(component.join("postrepair.tsv"), summary).unwrap();
        assert_eq!(
            commands::postrepair::selected_fasta(component).unwrap(),
            final_fasta
        );
        assert_eq!(
            commands::postrepair::selected_validation_dir(component).unwrap(),
            polish_round_dir(case, 3).join("03.validate")
        );
        let bad = dir.join("different-cohort.fastq");
        fs::write(&bad, "@wrong\nACGT\n+\nIIII\n").unwrap();
        let wrong = Command::new("gzip").arg("-c").arg(&bad).output().unwrap();
        assert!(wrong.status.success());
        fs::write(&source_reads, wrong.stdout).unwrap();
        assert!(commands::postrepair::selected_fasta(component).is_err());
        fs::write(&source_reads, compressed.stdout).unwrap();
        assert!(commands::postrepair::selected_fasta(component).is_ok());

        assert_eq!(
            fs::read_to_string(data.join("indel_event_candidates.tsv")).unwrap(),
            fixture.gate
        );
        // An independent residual review still stops the unchanged ordinary gate.
        write_with_parent(
            &homorepair_outcome_path(case),
            &outcome.replace("indel_review_count\t0", "indel_review_count\t1"),
        )
        .unwrap();
        assert!(homorepair_selected_fasta(&config, case, 2).is_err());
        write_with_parent(&homorepair_outcome_path(case), &outcome).unwrap();
        for changed in [
            junction.replace("\t5\t3\t3\t2\t2\tadopted", "\t5\t4\t4\t1\t1\tadopted"),
            junction.replace(
                "unique_main_junction_above_fixed_pool_remainder",
                "unique_main_event_above_fixed_pool_remainder",
            ),
            junction.replace("\t123\t5", "\t124\t5"),
        ] {
            write_with_parent(&data.join("indel_junction_candidates.tsv"), &changed).unwrap();
            assert!(homorepair_selected_fasta(&config, case, 2).is_err());
        }
        write_with_parent(&data.join("indel_junction_candidates.tsv"), &junction).unwrap();
        let proof = fixture.output.join("junction_read_classes.tsv");
        fs::remove_file(&proof).unwrap();
        assert!(homorepair_selected_fasta(&config, case, 2).is_err());
        fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn homorepair_handoff_unmeasured_flank20_requires_long_event_and_real_template_failure() {
        let dir = test_dir("homorepair unmeasured twenty handoff");
        let config = homorepair_test_config(&dir, true);
        let case = &config.cases[0];
        let input = polished_aln_path(case, 2);
        let hr = homorepair_out_dir(case);
        let data = hr.join("03.validate/01.data");
        let selected = hr.join("03.validate/03.reports/polished.fasta");
        let final_fasta = polished_aln_path(case, 3);
        let left = "TTGCAACCTTGGAACCTTGA";
        let right = "CCTAGGATTCGGAATCCGTA";
        let core = "A".repeat(101);
        let source_alt = format!("GT{}", &core[..100]);
        let complete_alt = format!("GT{core}");
        let competitor = format!("GT{}", "C".repeat(101));
        // Source context is unique. The provisional ALT has a second full
        // 20 bp template elsewhere, before its one supported nested edit.
        let source_sequence = format!("{left}GT{right}ATATAT{left}{source_alt}{right}");
        let selected_sequence = format!("{left}{complete_alt}{right}ATATAT{left}{source_alt}{right}");
        let write_sequences = |source: &str, revised: &str| {
            write_with_parent(&input, &format!(">source\n{source}\n")).unwrap();
            write_with_parent(&hr.join("01.inputs/linear_subgraph.round_hr.fasta"), &format!(">snapshot\n{source}\n")).unwrap();
            write_with_parent(&selected, &format!(">selected\n{revised}\n")).unwrap();
            write_with_parent(&final_fasta, &format!(">validated\n{revised}\n")).unwrap();
        };
        write_sequences(&source_sequence, &selected_sequence);
        homorepair_test_candidates(case, &[]);
        let header = "step\tsource_row\tsource_rows\tsource_pos\tsource_ref\tsource_alt\tsource_support\tsource_depth\tsource_total_ALT\tprobe_pos\tprobe_delta\tdepth\ttotal_ALT\tmain_ALT\tremainder\tuncertain_reads\tstatus\treason\tnested_repeat_edits\tedit_pos\tedit_ref\tedit_alt\ttrial_candidates\ttested_probes\tuntried_exact_reads\tentry_event_support\n";
        let common = format!("1\t1\t1\t21\tGT\t{source_alt}\t4\t10\t7\t21\t100");
        let sixty = format!("{header}{common}\t10\t6\t4\t4\t1\treview_required\thighest_source_ALT_failed_fixed_pool_gate\t0\t.\t.\t.\t1\t1\t0\t7\n");
        let twenty = format!("{header}{common}\t.\t.\t.\t.\t.\treview_required\t20bp_highest_ALT_full_template_not_unique\t0\t.\t.\t.\t1\t1\t.\t7\n");
        let event = format!("{header}{common}\t12\t8\t6\t4\t1\tadopted\tunique_main_event_above_fixed_pool_remainder\t1\t21\tGT\t{complete_alt}\t1\t1\t0\t7\n");
        let sixty_path = data.join("indel_candidates.tsv");
        let twenty_path = data.join("indel_flank20_candidates.tsv");
        let event_path = data.join("indel_event_candidates.tsv");
        for (path, text) in [(&sixty_path, &sixty), (&twenty_path, &twenty), (&event_path, &event)] {
            write_with_parent(path, text).unwrap();
        }
        let mut exact = "step\tread_id\texact_allele\tclass\n".to_string();
        let mut event_classes = "step\tread_id\tevent_class\n".to_string();
        for i in 0..12 {
            if i < 10 {
                let (allele, label) = if i < 4 { (complete_alt.as_str(), "ALT") }
                    else if i < 6 { (competitor.as_str(), "ALT") }
                    else if i < 9 { ("GT", "source") } else { (".", "uncertain") };
                writeln!(exact, "1\tr{i}\t{allele}\t{label}").unwrap();
            }
            let label = if i < 4 || i >= 10 { "selected_ALT_event" }
                else if i < 6 { "other_exact_ALT" } else if i < 9 { "source" } else { "uncertain" };
            writeln!(event_classes, "1\tr{i}\t{label}").unwrap();
        }
        let exact_path = data.join("indel/read_classes.tsv");
        let classes_path = data.join("indel/event_read_classes.tsv");
        let twenty_classes = data.join("indel/flank20_read_classes.tsv");
        write_with_parent(&exact_path, &exact).unwrap();
        write_with_parent(&classes_path, &event_classes).unwrap();
        write_fake_sv_pass_summary(&sv_summary_path(case,3)).unwrap();
        write_with_parent(&snv_indel_high_path(case,3), "pos\tref\talt\ttotal_count\tdepth\tcounts\n").unwrap();
        write_with_parent(&polish_report_path(case,3), "section\tname\tmetric\tvalue\nstage\tsv-eval:round_3\tstatus\tok\nstage\tsnv-indel-eval:round_3\tstatus\tok\nstage\tplot:round_3\tstatus\tok\n").unwrap();
        let outcome = format!("section\tname\tmetric\tvalue\noutcome\thomorepair\tstatus\tcomplete\noutcome\thomorepair\tsupport_status\tsupport_verified\noutcome\thomorepair\tsource_round\t2\noutcome\thomorepair\tinput_fasta\t{}\noutcome\thomorepair\tselected_fasta\t{}\noutcome\thomorepair\tstage_count\t1\noutcome\thomorepair\tnative_probe_count\t0\noutcome\thomorepair\tlocal_batch_count\t0\noutcome\thomorepair\tadopted_edit_count\t1\noutcome\thomorepair\tunresolved_count\t0\noutcome\thomorepair\tindel_candidate_count\t1\noutcome\thomorepair\tindel_adopted_count\t1\noutcome\thomorepair\tindel_review_count\t0\noutcome\thomorepair\tfinal_validation_status\tpassed\noutcome\thomorepair\tvalidation_round\t3\noutcome\thomorepair\tfinal_validation_round\t3\noutcome\thomorepair\tfinal_validation_dir\t{}\noutcome\thomorepair\tfinal_validation_fasta\t{}\n", input.display(), selected.display(), polish_round_dir(case,3).display(), final_fasta.display());
        write_with_parent(&homorepair_outcome_path(case), &outcome).unwrap();
        assert_eq!(homorepair_selected_fasta(&config, case, 2).unwrap(), final_fasta);
        // Other steps may have a 20 bp pool; this failed step may not.
        write_with_parent(&twenty_classes, "step\tread_id\texact_allele\tclass\n2\tother_step\tGT\tsource\n").unwrap();
        assert_eq!(homorepair_selected_fasta(&config, case, 2).unwrap(), final_fasta);
        let changed = |table: &str, field: &str, value: &str| {
            let mut lines = table.lines();
            let header = lines.next().unwrap();
            let index = header.split('\t').position(|key| key == field).unwrap();
            let mut values = lines.next().unwrap().split('\t').map(str::to_string).collect::<Vec<_>>();
            values[index] = value.into();
            format!("{header}\n{}\n", values.join("\t"))
        };
        for (field, value) in [("depth", "11"), ("reason", "20bp_fixed_pool_gate_failed"),
            ("reason", "forged_template_not_unique"), ("reason", "20bp_source_full_template_not_unique"),
            ("main_ALT", "0"), ("untried_exact_reads", "0"), ("nested_repeat_edits", "1"),
            ("edit_pos", "21"), ("source_support", "5"), ("probe_pos", "22"), ("entry_event_support", "8")] {
            write_with_parent(&twenty_path, &changed(&twenty, field, value)).unwrap();
            assert!(homorepair_selected_fasta(&config, case, 2).is_err(), "{field}={value}");
        }
        let numeric_twenty = changed(&changed(&sixty, "reason", "20bp_fixed_pool_gate_failed"), "depth", "11");
        write_with_parent(&twenty_path, &numeric_twenty).unwrap();
        assert!(homorepair_selected_fasta(&config, case, 2).is_err());
        write_with_parent(&twenty_path, &twenty).unwrap();
        let keep = ManualKeepDecision::Indel(ManualKeepEvent { pos: 21, reference: "GT".into(), alternative: source_alt.clone() });
        assert!(verify_ordinary_homorepair_inner(case, 2, config.polish.homorepair_frq, Some(&keep)).is_err());
        fs::remove_file(&event_path).unwrap();
        assert!(homorepair_selected_fasta(&config, case, 2).is_err());
        let short_path = data.join("indel_short_candidates.tsv");
        let short = changed(&changed(&event, "reason", "unique_main_ALT_above_fixed_pool_remainder"), "nested_repeat_edits", "0");
        write_with_parent(&short_path, &short).unwrap();
        assert!(homorepair_selected_fasta(&config, case, 2).is_err());
        fs::remove_file(&short_path).unwrap();
        write_with_parent(&event_path, &event).unwrap();
        write_with_parent(&twenty_classes, "step\tread_id\texact_allele\tclass\n1\tr0\tGT\tsource\n").unwrap();
        assert!(homorepair_selected_fasta(&config, case, 2).is_err());
        write_with_parent(&twenty_classes, "step\tread_id\texact_allele\tclass\n.\tr0\tGT\tsource\n").unwrap();
        assert!(homorepair_selected_fasta(&config, case, 2).is_err());
        fs::remove_file(&twenty_classes).unwrap();
        for depth in [".", "0", "9"] {
            write_with_parent(&sixty_path, &changed(&sixty, "depth", depth)).unwrap();
            assert!(homorepair_selected_fasta(&config, case, 2).is_err(), "60 depth {depth}");
        }
        write_with_parent(&sixty_path, &sixty).unwrap();
        write_with_parent(&classes_path, &event_classes.replacen("1\tr0\t", "1\tnew_read\t", 1)).unwrap();
        assert!(homorepair_selected_fasta(&config, case, 2).is_err());
        write_with_parent(&classes_path, &event_classes).unwrap();
        write_with_parent(&event_path, &changed(&event, "main_ALT", "4")).unwrap();
        assert!(homorepair_selected_fasta(&config, case, 2).is_err());
        write_with_parent(&event_path, &event).unwrap();
        // Consistent but unique source/ALT DNA cannot substantiate a fake failure.
        write_sequences(&format!("{left}GT{right}ATATAT{left}{competitor}{right}"),
            &format!("{left}{complete_alt}{right}ATATAT{left}{competitor}{right}"));
        assert!(homorepair_selected_fasta(&config, case, 2).is_err());
        // The alternative native reason is valid only when the source template repeats.
        write_sequences(&format!("{left}GT{right}ATATAT{left}GT{right}"),
            &format!("{left}{complete_alt}{right}ATATAT{left}GT{right}"));
        write_with_parent(&twenty_path, &changed(&twenty, "reason", "20bp_source_full_template_not_unique")).unwrap();
        assert_eq!(homorepair_selected_fasta(&config, case, 2).unwrap(), final_fasta);
        write_sequences(&source_sequence, &selected_sequence);
        write_with_parent(&twenty_path, &twenty).unwrap();
        write_with_parent(&final_fasta, ">wrong\nACGT\n").unwrap();
        assert!(homorepair_selected_fasta(&config, case, 2).is_err());
        fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn homorepair_run_reader_selects_only_its_outcome_and_rejects_malformed_rows() {
        let dir = test_dir("homorepair consolidated run");
        let config = homorepair_test_config(&dir, true);
        let case = &config.cases[0];
        let (_, selected) = homorepair_test_evidence(case, true);
        let path = homorepair_outcome_path(case);
        let header = "section\tname\tmetric\tvalue\n";
        let unrelated = concat!(
            "global\trun\tstatus\tfailed\n",
            "outcome\tprobe\tstatus\treview_required\n",
            "homorepair\tround_hr\telapsed_seconds\t999\n",
            "batch\tresult\tjson\t{\"status\":\"failed\",\"note\":\"one\\ntwo\\tthree\"}\n");
        let original = fs::read_to_string(&path).unwrap();
        let run = original.replacen(header, &format!("{header}{unrelated}"), 1)
            + "outcome\thomorepair\telapsed_seconds\t42\n";
        fs::write(&path, &run).unwrap();
        let values = read_homorepair_run_outcome(&path).unwrap();
        assert_eq!(values.get("status").map(String::as_str), Some("complete"));
        assert_eq!(values.get("elapsed_seconds").map(String::as_str), Some("42"));
        assert!(!values.contains_key("json"));
        assert_eq!(homorepair_selected_fasta(&config, case, 2).unwrap(), selected);
        let timings = runtime_stage_rows(&config).unwrap();
        assert_eq!(timings.iter().find(|row| row.stage == "polish homorepair").unwrap().seconds, Some(42.0));
        assert!(!timings.iter().any(|row| row.stage == "polish round 3 stage sum"));
        for malformed in [
            format!("{run}outcome\thomorepair\tstatus\tcomplete\n"),
            format!("{run}outcome\thomorepair\t\tvalue\n"),
            format!("{run}batch\tresult\tjson\tbad\textra\n"),
            run.replacen(header, "metric\tvalue\n", 1),
        ] {
            fs::write(&path, malformed).unwrap();
            assert!(read_homorepair_run_outcome(&path).is_err());
            assert!(homorepair_selected_fasta(&config, case, 2).is_err());
        }
        fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn homorepair_current_handoff_never_falls_back_to_legacy_files() {
        let dir = test_dir("homorepair no mixed layouts");
        let config = homorepair_test_config(&dir, true);
        let case = &config.cases[0];
        homorepair_test_evidence(case, true);
        let hr = homorepair_out_dir(case);
        let run = homorepair_outcome_path(case);
        let data = hr.join("03.validate/01.data");
        let reports = hr.join("03.validate/03.reports");
        fs::rename(data.join("candidates.tsv"), reports.join("candidates.tsv")).unwrap();
        write_with_parent(&data.join("homorepair_review.tsv"), "review_status\nsupport_verified\n").unwrap();
        write_with_parent(&data.join("changes.tsv"), "old_pos\tref\talt\tstatus\n3\tGT\tGTT\tapplied\n").unwrap();
        assert!(homorepair_selected_fasta(&config, case, 2).is_err());
        let values = read_homorepair_run_outcome(&run).unwrap();
        let legacy = hr.join("outcome.tsv");
        let text = format!("metric\tvalue\n{}", values.iter()
            .map(|(metric, value)| format!("{metric}\t{value}\n")).collect::<String>());
        write_with_parent(&legacy, &text).unwrap();
        fs::remove_file(&run).unwrap();
        assert_eq!(homorepair_existing_outcome_path(case), legacy);
        let mut summary = String::new();
        write_runtime_homorepair_statuses(&mut summary, &config).unwrap();
        assert!(summary.contains("complete") && summary.contains("outcome.tsv"));
        assert!(homorepair_selected_fasta(&config, case, 2).is_err());
        fs::write(&run, "metric\tvalue\nstatus\tcomplete\n").unwrap();
        assert_eq!(homorepair_existing_outcome_path(case), run);
        assert!(write_runtime_homorepair_statuses(&mut String::new(), &config).is_err());
        assert!(homorepair_selected_fasta(&config, case, 2).is_err());
        fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn homorepair_handoff_rejects_unconfirmed_or_inconsistent_results() {
        for fault in ["status","unresolved","source_round","stage_count","probe_count","final_status","validation_round","snapshot","sequence","changes","bound","decision","sv","unrelated","final_input","final_dir","plot"] {
            let dir = test_dir(&format!("homorepair fault {fault}"));
            let config = homorepair_test_config(&dir,true);
            let case = &config.cases[0];
            let (_,final_fasta) = homorepair_test_evidence(case,true);
            let hr = homorepair_out_dir(case);
            let outcome = homorepair_outcome_path(case);
            let data = hr.join("03.validate/01.data");
            let reports = hr.join("03.validate/03.reports");
            match fault {
                "status" | "unresolved" | "source_round" | "stage_count" | "probe_count" | "final_status" | "validation_round" => {
                    let (before,after) = match fault {
                        "status" => ("status\tcomplete\n","status\treview_required\n"),
                        "unresolved" => ("unresolved_count\t0","unresolved_count\t1"),
                        "source_round" => ("source_round\t2","source_round\t1"),
                        "stage_count" => ("stage_count\t1","stage_count\t2"),
                        "probe_count" => ("native_probe_count\t1","native_probe_count\t3"),
                        "final_status" => ("final_validation_status\tpassed","final_validation_status\treview_required"),
                        _ => ("validation_round\t3","validation_round\t4"),
                    };
                    fs::write(&outcome,fs::read_to_string(&outcome).unwrap().replace(before,after)).unwrap();
                }
                "snapshot" => fs::write(hr.join("01.inputs/linear_subgraph.round_hr.fasta"),">bad\nACGG\n").unwrap(),
                "sequence" => fs::write(reports.join("polished.fasta"),">bad\nACGTA\n").unwrap(),
                "changes" => fs::write(data.join("candidates.tsv"),"status\nkept\n").unwrap(),
                "bound" | "decision" => {
                    let path=data.join("candidates.tsv");
                    let (before,after)=if fault=="bound" {("100\t80\t60\t20","100\t40\t30\t60")} else {("suggest_modify","support")};
                    fs::write(&path,fs::read_to_string(&path).unwrap().replace(before,after)).unwrap();
                }
                "sv" => write_with_parent(&sv_summary_path(case,3),"metric\tvalue\nstatus\tfail\n").unwrap(),
                "unrelated" => write_with_parent(&snv_indel_high_path(case,3),"pos\tref\talt\ttotal_count\tdepth\tcounts\n1\tA\tC\t60\t100\t60\n").unwrap(),
                "final_input" => fs::write(final_fasta,">wrong_input\nACGT\n").unwrap(),
                "final_dir" => fs::write(&outcome,fs::read_to_string(&outcome).unwrap().replace(&format!("final_validation_dir\t{}",polish_round_dir(case,3).display()),&format!("final_validation_dir\t{}",hr.display()))).unwrap(),
                "plot" => {
                    let path=polish_report_path(case,3);
                    fs::write(&path,fs::read_to_string(&path).unwrap().replace("plot:round_3\tstatus\tok","plot:round_3\tstatus\tskipped")).unwrap();
                }
                _=> unreachable!(),
            }
            assert!(homorepair_selected_fasta(&config,case,2).is_err(),"{fault}");
            fs::remove_dir_all(dir).unwrap();
        }
    }

    #[test]
    fn homorepair_final_validation_uses_ordinary_reference_support_rule() {
        let dir = test_dir("homorepair ordinary final correction rule");
        let config = homorepair_test_config(&dir, true);
        let case = &config.cases[0];
        let (_, final_fasta) = homorepair_test_evidence(case, true);
        let high = snv_indel_high_path(case, 3);
        // Combined ALT is 55%, but REF=45 exceeds either individual ALT (30,25).
        write_with_parent(&high, "pos\tref\talt\ttotal_count\tdepth\tcounts\n1\tA\tC#G\t55\t100\t30#25\n").unwrap();
        assert_eq!(ordinary_validation_edit_count(&high).unwrap(), 0);
        assert_eq!(homorepair_selected_fasta(&config, case, 2).unwrap(), final_fasta);
        // An actual ordinary correction remains a handoff blocker.
        write_with_parent(&high, "pos\tref\talt\ttotal_count\tdepth\tcounts\n1\tA\tC#G\t70\t100\t60#10\n").unwrap();
        assert_eq!(ordinary_validation_edit_count(&high).unwrap(), 1);
        assert!(homorepair_selected_fasta(&config, case, 2).is_err());
        // A valid-looking support split cannot compensate for incomplete evidence.
        write_with_parent(&high, "ref\talt\ttotal_count\tdepth\tcounts\nA\tC#G\t55\t100\t30#25\n").unwrap();
        assert!(ordinary_validation_edit_count(&high).is_err());
        assert!(homorepair_selected_fasta(&config, case, 2).is_err());
        fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn homorepair_local_support_uses_unique_copy_winner_and_final_validation() {
        let dir=test_dir("homorepair local handoff");
        let config=homorepair_test_config(&dir,true);
        let case=&config.cases[0];
        let (_,final_fasta)=homorepair_test_evidence(case,true);
        let outcome = homorepair_outcome_path(case);
        fs::write(&outcome, fs::read_to_string(&outcome).unwrap().replace("local_batch_count\t0", "local_batch_count\t1")).unwrap();
        let mut row = homorepair_test_candidate(true);
        for (key, value) in [("review_status", "support_verified"), ("evidence_id", "batch_001"),
            ("depth", "80"), ("total_ALT", "80"), ("main_ALT", "60"), ("remainder", "0"),
            ("pool_reads", "100"), ("uncertain_reads", "20"), ("copy_counts", "1:20,2:60"),
            ("uncertain_reasons", "ambiguous:20"), ("context_changed", "0")] {
            homorepair_test_set_candidate(&mut row, key, value);
        }
        homorepair_test_candidates(case, &[row.clone()]);
        assert_eq!(homorepair_selected_fasta(&config,case,2).unwrap(),final_fasta);
        for (key, value) in [("copy_counts", "1:40,2:40"), ("copy_counts", "1:60,2:20"),
            ("copy_counts", "1:20,2:30,2:30"), ("pool_reads", "101"), ("depth", "81"),
            ("total_ALT", "60"), ("main_ALT", "61"), ("remainder", "20"),
            ("uncertain_reasons", "ambiguous:19"), ("context_changed", "."),
            ("evidence_id", "probe_1"), ("evidence_id", "batch_002"), ("evidence_id", "batch_1")] {
            let mut bad = row.clone();
            homorepair_test_set_candidate(&mut bad, key, value);
            if value == "1:40,2:40" { homorepair_test_set_candidate(&mut bad, "main_ALT", "40"); }
            homorepair_test_candidates(case, &[bad]);
            assert!(homorepair_selected_fasta(&config,case,2).is_err(), "{key}={value}");
        }
        homorepair_test_candidates(case, &[row]);
        fs::write(&sv_summary_path(case,3),"metric\tvalue\nstatus\tfail\n").unwrap();
        assert!(homorepair_selected_fasta(&config,case,2).is_err());
        fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn homorepair_compact_candidates_require_original_anchored_edits() {
        let dir = test_dir("homorepair compact anchored edits");
        let config = homorepair_test_config(&dir, true);
        let case = &config.cases[0];
        let (_, selected) = homorepair_test_evidence(case, true);
        let row = homorepair_test_candidate(true);
        for (key, value) in [("edit_pos", "."), ("edit_ref", "."), ("edit_alt", "."),
            ("edit_ref", "GG"), ("edit_alt", "GTA"), ("selected_copies", "3"),
            ("reference_copies", "2"), ("source_rows", "1,1"),
            ("review_status", "unknown"), ("review_status", "review_required")] {
            let mut bad = row.clone();
            homorepair_test_set_candidate(&mut bad, key, value);
            homorepair_test_candidates(case, &[bad]);
            assert!(homorepair_selected_fasta(&config, case, 2).is_err(), "{key}={value}");
        }
        // This generic insertion makes the same FASTA, but is not the recorded
        // tract's original-reference anchored allele.
        let mut shifted_anchor = row.clone();
        for (key, value) in [("edit_pos", "4"), ("edit_ref", "T"), ("edit_alt", "TT")] {
            homorepair_test_set_candidate(&mut shifted_anchor, key, value);
        }
        homorepair_test_candidates(case, &[shifted_anchor]);
        assert!(homorepair_selected_fasta(&config, case, 2).is_err());
        homorepair_test_candidates(case, &[row.clone(), row.clone()]);
        assert!(homorepair_selected_fasta(&config, case, 2).is_err());
        homorepair_test_candidates(case, &[row]);
        assert_eq!(homorepair_selected_fasta(&config, case, 2).unwrap(), selected);
        fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn homorepair_compact_candidates_bind_native_probe_ids_to_run_metadata() {
        let dir = test_dir("homorepair compact probe identity");
        let config = homorepair_test_config(&dir, true);
        let case = &config.cases[0];
        let (_, selected) = homorepair_test_evidence(case, true);
        let mut row = homorepair_test_candidate(true);
        for id in ["probe_2", "probe_01", "batch_001", "."] {
            homorepair_test_set_candidate(&mut row, "evidence_id", id);
            homorepair_test_candidates(case, &[row.clone()]);
            assert!(homorepair_selected_fasta(&config, case, 2).is_err(), "{id}");
        }
        homorepair_test_set_candidate(&mut row, "evidence_id", "probe_1");
        homorepair_test_candidates(case, &[row]);
        let path = homorepair_outcome_path(case);
        let run = fs::read_to_string(&path).unwrap();
        for invalid in [
            run.replace("probe\tprobe_1\tevidence_round\t2", "probe\tprobe_1\tevidence_round\t3"),
            run.replace("probe\tprobe_1\tSV_status\tpass", "probe\tprobe_1\tSV_status\tnot_pass"),
            run.replace("probe\tprobe_1\tcandidate_count\t1", "probe\tprobe_1\tcandidate_count\t0"),
            run.replace("probe\tprobe_1\t", "probe\tprobe_2\t"),
            format!("{run}probe\tprobe_1\tSV_status\tpass\n"),
            run.lines().filter(|line| !line.starts_with("probe\t")).map(|line| format!("{line}\n")).collect(),
        ] {
            fs::write(&path, invalid).unwrap();
            assert!(homorepair_selected_fasta(&config, case, 2).is_err());
        }
        fs::write(&path, &run).unwrap();
        assert_eq!(homorepair_selected_fasta(&config, case, 2).unwrap(), selected);
        fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn homorepair_local_probe_copy_can_be_the_unique_winner_without_edits() {
        let dir = test_dir("homorepair local probe winner");
        let config = homorepair_test_config(&dir, true);
        let case = &config.cases[0];
        let (_, source) = homorepair_test_evidence(case, false);
        let path = homorepair_outcome_path(case);
        fs::write(&path, fs::read_to_string(&path).unwrap().replace("local_batch_count\t0", "local_batch_count\t1")).unwrap();
        let mut row = homorepair_test_candidate(false);
        for (key, value) in [("review_status", "support_verified"), ("evidence_id", "batch_001"),
            ("probe_end", "4"), ("probe_copies", "1"), ("depth", "80"), ("total_ALT", "20"),
            ("main_ALT", "20"), ("remainder", "60"), ("pool_reads", "80"), ("uncertain_reads", "0"),
            ("copy_counts", "1:60,2:20"), ("uncertain_reasons", ""), ("context_changed", "1")] {
            homorepair_test_set_candidate(&mut row, key, value);
        }
        homorepair_test_candidates(case, &[row]);
        assert_eq!(homorepair_selected_fasta(&config, case, 2).unwrap(), source);
        fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn homorepair_no_candidates_is_unchanged_but_not_site_confirmation() {
        let dir=test_dir("homorepair no candidates");
        let config=homorepair_test_config(&dir,true);
        let case=&config.cases[0];
        let (_,selected)=homorepair_test_evidence(case,false);
        homorepair_test_candidates(case, &[]);
        let outcome=homorepair_outcome_path(case);
        let text=fs::read_to_string(&outcome).unwrap().replace("support_status\tsupport_verified","support_status\tno_candidates").replace("native_probe_count\t1","native_probe_count\t0")
            .lines().filter(|line| !line.starts_with("probe\t")).map(|line| format!("{line}\n")).collect::<String>();
        fs::write(&outcome,text).unwrap();
        assert_eq!(homorepair_selected_fasta(&config,case,2).unwrap(),selected);
        for status in ["review_required", "unknown", "screened_out"] {
            let mut row = homorepair_test_candidate(false);
            homorepair_test_set_candidate(&mut row, "review_status", status);
            homorepair_test_set_candidate(&mut row, "selected_copies", ".");
            homorepair_test_candidates(case, &[row]);
            assert!(homorepair_selected_fasta(&config,case,2).is_err(), "{status}");
        }
        fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn supplemental_unchanged_handoff_replays_rows_and_rejects_mixed_or_changed_evidence() {
        let dir = test_dir("supplemental unchanged handoff");
        let mut config = homorepair_test_config(&dir, true);
        config.cases[0].polish_out_dir = dir.join("old/04.polish");
        config.cases[0].postrepair_out_dir = Some(dir.join("new/05.postrepair"));
        let case = &config.cases[0];
        let (source, _) = homorepair_test_evidence(case, false);
        let sequence = format!("ACG{}CGA", "T".repeat(11));
        let fasta = format!(">source\n{sequence}\n");
        fs::write(&source, &fasta).unwrap();
        fs::write(
            homorepair_out_dir(case).join("01.inputs/linear_subgraph.round_hr.fasta"),
            &fasta,
        )
        .unwrap();
        homorepair_test_candidates(case, &[]);
        let run = homorepair_outcome_path(case);
        let text = fs::read_to_string(&run)
            .unwrap()
            .replace(
                "support_status\tsupport_verified",
                "support_status\tunchanged",
            )
            .replace("native_probe_count\t1", "native_probe_count\t0")
            .lines()
            .filter(|line| !line.starts_with("probe\t"))
            .map(|line| format!("{line}\n"))
            .collect::<String>();
        let mut text = format!("{text}outcome\thomorepair\treason\tsupplemental_repeat_candidates_unchanged\noutcome\thomorepair\tfailed_SV_frames\t.\n");
        for key in [
            "pending_count",
            "uncertain_count",
            "indel_candidate_count",
            "indel_adopted_count",
            "indel_review_count",
            "indel_nested_repeat_edits",
            "indel_probe_count",
            "repeat_extra_adopted_count",
            "repeat_extra_review_count",
            "circle_end_candidate_count",
            "circle_end_adopted_count",
            "circle_end_review_count",
        ] {
            writeln!(text, "outcome\thomorepair\t{key}\t0").unwrap();
        }
        text.push_str("outcome\thomorepair\trepeat_extra_candidate_count\t1\n");
        fs::write(&run, &text).unwrap();
        let table =
            homorepair_out_dir(case).join("03.validate/01.data/repeat_extra/candidates.tsv");
        let rows = "source_start\tsource_end\tmotif\toriginal_copies\tselected_copies\tentry\tstatus\treason\tdepth\ttotal_ALT\tmain_ALT\tremainder\tuntyped\tcopy_counts\n4\t14\tT\t11\t.\tlong_homopolymer_residual\tunchanged\tno_unique_dominant_alternative_in_complete_reads\t100\t20\t20\t80\t0\t10:20,11:80\n";
        write_with_parent(&table, rows).unwrap();
        write_fake_sv_pass_summary(&sv_summary_path(case, 2)).unwrap();
        let sv = fs::read_to_string(sv_summary_path(case, 2)).unwrap();
        let high = "pos\tref\talt\tcounts\tdepth\ttype\n";
        write_with_parent(&snv_indel_high_path(case, 2), high).unwrap();
        assert_eq!(homorepair_selected_fasta(&config, case, 2).unwrap(), source);
        let component = homorepair_out_dir(case).parent().unwrap().to_path_buf();
        let source_round = polish_round_dir(case, 2);
        let mut external_case = case.clone();
        external_case.postrepair_source_dir = Some(source_round.clone());
        external_case.polish_out_dir = dir.join("new/04.polish");
        assert_eq!(
            verify_ordinary_postrepair(&component, &repair_source_dir(&external_case, 2), 0.2)
                .unwrap(),
            source
        );
        assert!(!external_case.polish_out_dir.exists());

        for fault in [
            "missing_table",
            "row_count",
            "review",
            "selected_alt",
            "reason",
            "source",
            "probe",
            "adopted",
            "ordinary",
            "indel",
            "circle",
            "pending",
            "uncertain",
            "failed_frame",
            "sv",
            "high",
        ] {
            fs::write(&run, &text).unwrap();
            fs::write(&table, rows).unwrap();
            fs::write(&source, &fasta).unwrap();
            homorepair_test_candidates(case, &[]);
            fs::write(sv_summary_path(case, 2), &sv).unwrap();
            fs::write(snv_indel_high_path(case, 2), high).unwrap();
            let change = |key: &str, replacement: &str| {
                fs::write(&run, text.replace(key, replacement)).unwrap();
            };
            match fault {
                "missing_table" => fs::remove_file(&table).unwrap(),
                "row_count" => change(
                    "repeat_extra_candidate_count\t1",
                    "repeat_extra_candidate_count\t2",
                ),
                "review" => {
                    fs::write(&table, rows.replace("\tunchanged\t", "\treview_required\t")).unwrap()
                }
                "selected_alt" => {
                    fs::write(&table, rows.replace("11\t.\tlong_", "11\t12\tlong_")).unwrap()
                }
                "reason" => fs::write(
                    &table,
                    rows.replace(
                        "no_unique_dominant_alternative_in_complete_reads",
                        "strict_fixed_pool_support_not_pass",
                    ),
                )
                .unwrap(),
                "source" => fs::write(&source, format!(">source\n{sequence}A\n")).unwrap(),
                "probe" => change("native_probe_count\t0", "native_probe_count\t1"),
                "adopted" => change("adopted_edit_count\t0", "adopted_edit_count\t1"),
                "ordinary" => homorepair_test_candidates(case, &[homorepair_test_candidate(false)]),
                "indel" => change("indel_candidate_count\t0", "indel_candidate_count\t1"),
                "circle" => change(
                    "circle_end_candidate_count\t0",
                    "circle_end_candidate_count\t1",
                ),
                "pending" => change("pending_count\t0", "pending_count\t1"),
                "uncertain" => change("uncertain_count\t0", "uncertain_count\t1"),
                "failed_frame" => change("failed_SV_frames\t.", "failed_SV_frames\tprobe_1"),
                "sv" => fs::write(sv_summary_path(case, 2), sv.replace("pass", "failed")).unwrap(),
                "high" => fs::write(
                    snv_indel_high_path(case, 2),
                    format!("{high}1\tA\tC\t80\t100\tSNV\n"),
                )
                .unwrap(),
                _ => unreachable!(),
            }
            assert!(
                homorepair_selected_fasta(&config, case, 2).is_err(),
                "{fault}"
            );
            assert!(
                verify_ordinary_postrepair(&component, &repair_source_dir(&external_case, 2), 0.2)
                    .is_err(),
                "{fault}"
            );
        }
        fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn homorepair_legacy_outcome_is_visible_but_never_a_current_handoff() {
        let dir=test_dir("homorepair historical summary");
        let config=homorepair_test_config(&dir,true);
        let case=&config.cases[0];
        let legacy=polish_round_dir(case,2).parent().unwrap().join("round_hr2/outcome.tsv");
        write_with_parent(&legacy,"metric\tvalue\nstatus\tcomplete\nsource_round\t2\nvalidation_round\t4\n").unwrap();
        assert_eq!(homorepair_existing_outcome_path(case),legacy);
        assert!(homorepair_selected_fasta(&config,case,2).is_err());
        fs::create_dir_all(homorepair_out_dir(case)).unwrap();
        assert_eq!(homorepair_existing_outcome_path(case),homorepair_outcome_path(case));
        fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn homorepair_legacy_single_hr_timing_uses_only_its_final_validation_round() {
        let dir = test_dir("homorepair single HR historical timing");
        let config = homorepair_test_config(&dir, true);
        let case = &config.cases[0];
        let legacy = homorepair_out_dir(case).join("outcome.tsv");
        write_with_parent(&legacy, "metric\tvalue\nstatus\tcomplete\nsource_round\t2\nvalidation_round\t4\nfinal_validation_round\t3\nelapsed_seconds\t7\n").unwrap();
        for round in [3, 4] {
            write_with_parent(&polish_report_path(case, round),
                "section\tname\tmetric\tvalue\nstage\tvalidation\telapsed_seconds\t5\n").unwrap();
        }
        let rows = runtime_stage_rows(&config).unwrap();
        assert!(!rows.iter().any(|row| row.stage == "polish round 3 stage sum"));
        assert!(rows.iter().any(|row| row.stage == "polish round 4 stage sum"));
        assert_eq!(rows.iter().find(|row| row.stage == "polish homorepair").unwrap().seconds, Some(7.0));
        assert!(homorepair_selected_fasta(&config, case, 2).is_err());
        fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn homorepair_shell_uses_shared_read_only_check_and_stops_on_failure() {
        let dir = test_dir("homorepair shell handoff");
        let config = homorepair_test_config(&dir, true);
        let case = &config.cases[0];
        let selected = polish_round_input_path(case, 3);
        let mut shell = String::from("set -euo pipefail\nORGRAFT_BIN=orgraft_stub\norgraft_stub() {\n  if [[ $1 == postrepair ]]; then return 0; fi\n");
        writeln!(shell, "  [[ $# == 8 && $1 == workflow && $2 == homorepair-check && $3 == --config && $4 == {} && $5 == --case && $6 == {} && $7 == --round && $8 == 2 ]] || return 98",
            shell_quote(&config.config_path), shell_quote_str(&case.name)).unwrap();
        writeln!(shell, "  if [[ $CHECK_STATUS != 0 ]]; then return \"$CHECK_STATUS\"; fi\n  printf '%s\\n' {}\n}}", shell_quote(&selected)).unwrap();
        write_homorepair_after_checkpoint(&mut shell, &config, case, 2).unwrap();
        writeln!(shell, "printf '%s\\n' \"$final_polished\"").unwrap();
        for code in [0, 7] {
            let output = Command::new("bash").arg("-c").arg(&shell).env("CHECK_STATUS", code.to_string()).output().unwrap();
            assert_eq!(output.status.code(), Some(code), "{}", String::from_utf8_lossy(&output.stderr));
            if code == 0 { assert_eq!(String::from_utf8_lossy(&output.stdout).trim(), selected.display().to_string()); }
            else { assert!(output.stdout.is_empty()); }
        }
        fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn runtime_summary_shows_recorded_repair_counts_and_retained_review_results() {
        let dir = test_dir("repair runtime summary");
        let mut config = homorepair_test_config(&dir, false);
        let root = dir.join("05.postrepair");
        config.postrepair = Some(WorkflowPostrepairConfig { enabled: true, out_dir: root.clone(), threads: None, frequency: 0.2, review_files: BTreeMap::new() });
        config.cases[0].postrepair_out_dir = Some(root);
        config.refine = Some(WorkflowRefineConfig { enabled: true, out_dir: dir.join("07.refine"), threads: None, bundle: None, organelle_reads: None, adopt: BTreeMap::new(), intra_measurement: BTreeMap::new(), duplicate_template_full_span: BTreeMap::new(), linear_chain_flanks: BTreeMap::new(), graph_full_family: BTreeMap::new(), graph_measurement: BTreeMap::new() });
        let case = &config.cases[0];
        let component = homorepair_out_dir(case).parent().unwrap().to_path_buf();
        write_with_parent(&refine_out_dir(&config, &case.organelle).join("auto_complete.tsv"), "key\tvalue\nstatus\tpartial_requires_review\nretained_bundle\t/preserved components.tsv\n").unwrap();
        write_with_parent(&component.join("postrepair.tsv"), "metric\tvalue\nstatus\treview_required\nunresolved_count\t2\nselection_method\tjoint_ir\nselected_fasta\t/preserved candidate.fasta\n").unwrap();
        write_with_parent(&homorepair_outcome_path(case), "section\tname\tmetric\tvalue\noutcome\thomorepair\tadopted_edit_count\t3\noutcome\thomorepair\tindel_candidate_count\t4\nhomorepair\tround_hr\tcandidate_count\t12\n").unwrap();
        let mut summary = String::new();
        write_runtime_repair_results(&mut summary,&config).unwrap();
        assert!(summary.contains("| review_required | 12 | 4 | 3 | 2 |"));
        assert!(summary.contains("[sequence](</preserved candidate.fasta>)"));
        assert!(summary.contains("| plastid | partial_requires_review | [components](</preserved components.tsv>) |"));
        assert!(!summary.contains("confirmed"));
        fs::remove_file(component.join("postrepair.tsv")).unwrap();
        summary.clear();write_runtime_repair_results(&mut summary,&config).unwrap();
        assert!(summary.contains("| pending | 12 | 4 | 3 | not recorded | pending |"));
        fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn homorepair_runtime_does_not_count_native_validation_twice() {
        let dir = test_dir("homorepair runtime");
        let mut config = homorepair_test_config(&dir, true);
        let case = &config.cases[0];
        homorepair_test_evidence(case, true);
        write_with_parent(&polish_report_path(case, 3), "section\tname\tmetric\tvalue\nstage\tvalidation\telapsed_seconds\t5\n").unwrap();
        let rows = runtime_stage_rows(&config).unwrap();
        assert!(!rows.iter().any(|r| r.stage == "polish round 3 stage sum"));
        assert!(rows.iter().any(|r| r.stage == "polish homorepair"));
        config.polish.homorepair = false;
        let rows = runtime_stage_rows(&config).unwrap();
        assert!(rows.iter().any(|r| r.stage == "polish round 3 stage sum"));
        assert!(!rows.iter().any(|r| r.stage == "polish homorepair"));
        fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn homorepair_enabled_does_not_change_standalone_checkpoint() {
        let dir = test_dir("homorepair_checkpoint");
        let config = homorepair_test_config(&dir, true);
        let case = &config.cases[0];
        write_fake_sv_pass_summary(&sv_summary_path(case, 1)).unwrap();
        write_with_parent(&snv_indel_high_path(case, 1), "pos\tref\talt\ttype\ttotal_count\tdepth\tcounts\tfixed_ref\n").unwrap();
        assert_eq!(checkpoint2_impl(&config, case, 1, true, false, None).unwrap(), Checkpoint2Status::Complete);
        assert!(!homorepair_out_dir(case).exists());
        fs::remove_dir_all(dir).unwrap();
    }

fn circle_ownership_native_fixture(
    dir: &Path,
) -> (
    PathBuf,
    PathBuf,
    String,
    Vec<BTreeMap<String, String>>,
    BTreeMap<String, String>,
) {
    let source = dir.join("source/round_1");
    let data = dir.join("round_hr/03.validate/01.data");
    let original = "CAAAAGCCCCGTGTGT".to_string();
    let rotated = format!("{}{}", &original[8..], &original[..8]);
    let variants = "pos\tref\talt\tcounts\tdepth\ttype\n";
    write_with_parent(
        &source.join("03.validate/01.data/snv_indel_variants_combined.tsv"),
        &format!("{variants}1\tCAAAA\tCAAA\t8\t10\thomopolymer\n"),
    )
    .unwrap();
    write_with_parent(&data.join("probe_1.fasta"), ">x\nCAAAGCCCCGTGTGT\n").unwrap();
    write_with_parent(
        &data.join("probe_1/01.inputs/linear_subgraph.round_2.fasta"),
        ">x\nCAAAGCCCCGTGTGT\n",
    )
    .unwrap();
    write_with_parent(
        &data.join("probe_1/03.validate/01.data/snv_indel_variants_combined.tsv"),
        &format!("{variants}1\tCAAA\tCAAAA\t8\t10\thomopolymer\n"),
    )
    .unwrap();
    write_with_parent(
        &data.join("probe_1/03.validate/03.reports/sv_snv_indel_summary.tsv"),
        "metric\tvalue\nstatus\tpass\n",
    )
    .unwrap();
    write_with_parent(
        &data.join("circle_end/rotated.fasta"),
        &format!(">x\n{rotated}\n"),
    )
    .unwrap();
    write_with_parent(
        &data.join("circle_end/rotated_validate/01.inputs/linear_subgraph.round_2.fasta"),
        &format!(">x\n{rotated}\n"),
    )
    .unwrap();
    write_with_parent(
        &data.join(
            "circle_end/rotated_validate/03.validate/01.data/snv_indel_variants_combined.tsv",
        ),
        &format!("{variants}9\tCAAAA\tCAAA\t8\t10\thomopolymer\n"),
    )
    .unwrap();
    write_with_parent(
        &data.join("circle_end/rotated_validate/03.validate/03.reports/sv_snv_indel_summary.tsv"),
        "metric\tvalue\nstatus\tpass\n",
    )
    .unwrap();
    write_with_parent(
        &data.join("circle_end/decision.tsv"),
        "metric\tvalue\nordinary_source_rows\t1\nordinary_evidence_id\tprobe_1\n",
    )
    .unwrap();
    write_with_parent(
        &data.join("circle_end/read_classes.tsv"),
        "source_rows\tread_id\tclass\texact_copies\n",
    )
    .unwrap();
    let header = "source_rows\tsource_length\tcut_0based\trotated_start_0based\trotated_end_0based\toriginal_start_0based\toriginal_end_0based\tmotif\treference_copies\tselected_copies\tdepth\ttotal_ALT\tmain_ALT\tremainder\tuncertain_reads\tstatus\treason\n";
    write_with_parent(&data.join("circle_end/candidates.tsv"), &format!("{header}1\t16\t8\t9\t13\t1\t5\tA\t4\t4\t.\t.\t.\t.\t.\tskipped\tordinary_probe_confirms_original_copies\n")).unwrap();
    let row = "1\t2\t5\tA\t4\t4\tsupport\tdominant_probe_ALT_restores_original_copies\tprobe_1\t2\t4\t3\t10\t8\t8\t2\t.\t.\t.\t.\t.\t.\t.\t.";
    let ordinary = vec![HOMOREPAIR_CANDIDATE_HEADER
        .split('\t')
        .zip(row.split('\t'))
        .map(|(k, v)| (k.to_string(), v.to_string()))
        .collect()];
    let outcome = BTreeMap::from([
        ("circle_end_candidate_count".into(), "1".into()),
        ("circle_end_adopted_count".into(), "0".into()),
        ("circle_end_review_count".into(), "0".into()),
        (
            "circle_end_reason".into(),
            "ordinary_probe_confirms_original_copies".into(),
        ),
    ]);
    (source, data, original, ordinary, outcome)
}

#[test]
fn circle_ordinary_handoff_replays_native_evidence_and_retains_source() {
    let dir = test_dir("circle ordinary replay");
    let (source, data, original, ordinary, outcome) = circle_ownership_native_fixture(&dir);
    let mut selected = original.clone();
    assert_eq!(
        replay_circle_end_edit(
            &data,
            &outcome,
            &mut selected,
            &source,
            &original,
            &ordinary,
            0,
            0.2
        )
        .unwrap(),
        0
    );
    assert_eq!(selected, original);
    let path = data.join("probe_1/03.validate/01.data/snv_indel_variants_combined.tsv");
    let saved = fs::read_to_string(&path).unwrap();
    fs::write(&path, saved.replace("\t8\t10\t", "\t2\t10\t")).unwrap();
    assert!(replay_circle_end_edit(
        &data,
        &outcome,
        &mut selected,
        &source,
        &original,
        &ordinary,
        0,
        0.2
    )
    .is_err());
    fs::remove_dir_all(dir).unwrap();
}

#[test]
fn circle_ordinary_handoff_rejects_changed_sequence_edits_or_saved_claims() {
    let dir = test_dir("circle ordinary invalid claims");
    let (source, data, original, ordinary, outcome) = circle_ownership_native_fixture(&dir);
    let mut selected = original.clone();
    assert!(replay_circle_end_edit(
        &data,
        &outcome,
        &mut selected,
        &source,
        &original,
        &ordinary,
        1,
        0.2
    )
    .is_err());
    selected.push('A');
    assert!(replay_circle_end_edit(
        &data,
        &outcome,
        &mut selected,
        &source,
        &original,
        &ordinary,
        0,
        0.2
    )
    .is_err());
    for (field, value) in [
        ("reason", "different_reason"),
        ("evidence_id", "batch_001"),
        ("selected_copies", "3"),
        ("review_status", "review_required"),
        ("main_ALT", "9"),
        ("tract_start", "3"),
    ] {
        let mut rows = ordinary.clone();
        rows[0].insert(field.into(), value.into());
        assert!(
            replay_circle_end_edit(
                &data,
                &outcome,
                &mut original.clone(),
                &source,
                &original,
                &rows,
                0,
                0.2
            )
            .is_err(),
            "{field}"
        );
    }
    let mut duplicate = ordinary.clone();
    duplicate.push(ordinary[0].clone());
    assert!(replay_circle_end_edit(
        &data,
        &outcome,
        &mut original.clone(),
        &source,
        &original,
        &duplicate,
        0,
        0.2
    )
    .is_err());
    fs::remove_dir_all(dir).unwrap();
}

#[test]
fn circle_ordinary_handoff_requires_real_source_probe_and_rotation_binding() {
    for kind in [
        "source",
        "probe",
        "SV",
        "rotated",
        "circle_source",
        "circle_classes",
        "link",
    ] {
        let dir = test_dir(&format!("circle ordinary binding {kind}"));
        let (source, data, original, ordinary, outcome) = circle_ownership_native_fixture(&dir);
        let (path,from,to)=match kind {
                "source" => (source.join("03.validate/01.data/snv_indel_variants_combined.tsv"), "\t8\t10\t", "\t1\t10\t"),
                "probe" => (data.join("probe_1.fasta"), "CAAAG", "CAAG"),
                "SV" => (data.join("probe_1/03.validate/03.reports/sv_snv_indel_summary.tsv"), "pass", "not_pass"),
                "rotated" => (data.join("circle_end/rotated.fasta"), "AAAA", "AAA"),
                "circle_source" => (data.join("circle_end/rotated_validate/03.validate/01.data/snv_indel_variants_combined.tsv"), "\t8\t10\t", "\t1\t10\t"),
                "circle_classes" => (data.join("circle_end/read_classes.tsv"), "exact_copies\n", "exact_copies\n1\tr1\tsource\t4\n"),
                _ => (data.join("circle_end/decision.tsv"), "probe_1", "probe_2"),
            };
        let old = fs::read_to_string(&path).unwrap();
        assert!(old.contains(from));
        fs::write(&path, old.replace(from, to)).unwrap();
        assert!(
            replay_circle_end_edit(
                &data,
                &outcome,
                &mut original.clone(),
                &source,
                &original,
                &ordinary,
                0,
                0.2
            )
            .is_err(),
            "{kind}"
        );
        fs::remove_dir_all(dir).unwrap();
    }
}

    #[test]
    fn circle_end_handoff_replays_fixed_pool_and_rejects_a_changed_vote() {
        let dir = test_dir("circle_end_handoff");
        let root = dir.join("circle_end");
        fs::create_dir_all(&root).unwrap();
        let candidate = root.join("candidates.tsv");
        let classes = root.join("read_classes.tsv");
        let header = "source_rows\tsource_length\tcut_0based\trotated_start_0based\trotated_end_0based\toriginal_start_0based\toriginal_end_0based\tmotif\treference_copies\tselected_copies\tdepth\ttotal_ALT\tmain_ALT\tremainder\tuncertain_reads\tstatus\treason\n";
        let row = "7\t16\t8\t8\t12\t0\t4\tA\t4\t3\t10\t6\t6\t4\t2\tadopted\tunique_main_ALT_above_fixed_read_pool_remainder\n";
        fs::write(&candidate, format!("{header}{row}")).unwrap();
        let mut votes = "source_rows\tread_id\tclass\texact_copies\n".to_string();
        for i in 0..6 { votes.push_str(&format!("7\ta{i}\tALT\t3\n")); }
        for i in 0..2 { votes.push_str(&format!("7\tr{i}\tsource\t4\n")); }
        for i in 0..2 { votes.push_str(&format!("7\tu{i}\tuncertain\t.\n")); }
        fs::write(&classes, &votes).unwrap();
        let outcome = BTreeMap::from([
            ("circle_end_candidate_count".into(), "1".into()),
            ("circle_end_adopted_count".into(), "1".into()),
            ("circle_end_review_count".into(), "0".into()),
            ("circle_end_reason".into(), "unique_main_ALT_above_fixed_read_pool_remainder".into()),
        ]);
        let mut source = "AAAACCCCGGGGTTTT".to_string();
        assert_eq!(replay_circle_end_edit(&dir, &outcome, &mut source, &dir, "AAAACCCCGGGGTTTT", &[], 0, 0.2).unwrap(), 1);
        assert_eq!(source, "AAACCCCGGGGTTTT");
        fs::write(&classes, votes.replacen("7\ta0\tALT\t3", "7\ta0\tuncertain\t.", 1)).unwrap();
        let mut source = "AAAACCCCGGGGTTTT".to_string();
        assert!(replay_circle_end_edit(&dir, &outcome, &mut source, &dir, "AAAACCCCGGGGTTTT", &[], 0, 0.2).is_err());
        fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn supplemental_repeat_handoff_replays_the_fixed_pool_before_later_edits() {
        let dir=test_dir("repeat extra handoff");
        let root=dir.join("repeat_extra");
        let original=format!("{}{}{}","G".repeat(30),"C".repeat(12),"A".repeat(30));
        let selected=format!("{}{}{}","G".repeat(30),"C".repeat(13),"A".repeat(30));
        let header="source_start\tsource_end\tmotif\toriginal_copies\tselected_copies\tentry\tstatus\treason\tdepth\ttotal_ALT\tmain_ALT\tremainder\tuntyped\tcopy_counts\n";
        let row="31\t42\tC\t12\t13\tlong_homopolymer_residual\tsupport_verified\tunique_main_ALT_above_fixed_pool_remainder\t10\t7\t7\t3\t1\t12:2,13:7\n";
        let candidate=root.join("candidates.tsv");
        write_with_parent(&candidate,&format!("{header}{row}")).unwrap();
        let mut classes="read_id\tcopy_count\n".to_string();
        let mut native="read_id\ttarget_id\ttarget_start\ttarget_end\talignment_role\n".to_string();
        for i in 0..10 {
            let copy=if i<7 {"13"}else if i<9 {"12"}else{"unknown"};
            writeln!(classes,"r{i}\t{copy}").unwrap();
            writeln!(native,"r{i}\tx\t1\t72\tprimary").unwrap();
        }
        let class_path=root.join("site_31_read_classes.tsv");
        write_with_parent(&class_path,&classes).unwrap();
        for (model,sequence) in [("source_0",&original),("site_31",&selected)] {
            write_with_parent(&root.join(model).join("candidate.fasta"),&format!(">x\n{sequence}\n")).unwrap();
            write_with_parent(&root.join(model).join("round_2/03.validate/01.data/sv_read_evidence.tsv"),&native).unwrap();
        }
        let outcome=BTreeMap::from([("repeat_extra_candidate_count".into(),"1".into()),("repeat_extra_adopted_count".into(),"1".into()),("repeat_extra_review_count".into(),"0".into())]);
        let mut replay=original.clone();
        assert_eq!(replay_repeat_extra(&dir,&outcome,&original,&[],&mut replay).unwrap(),(1,1));
        assert_eq!(replay,selected);
        fs::write(&class_path,classes.replace("r0\t13","outsider\t13")).unwrap();
        assert!(replay_repeat_extra(&dir,&outcome,&original,&[],&mut original.clone()).is_err());
        fs::write(&class_path,&classes).unwrap();
        fs::write(&candidate,format!("{header}{}",row.replace("\t7\t3\t1\t12:2", "\t3\t7\t1\t12:2"))).unwrap();
        assert!(replay_repeat_extra(&dir,&outcome,&original,&[],&mut original.clone()).is_err());
        fs::write(&candidate,format!("{header}{row}")).unwrap();
        let owned=vec![BTreeMap::from([("tract_start".into(),"31".into()),("tract_end".into(),"42".into())])];
        assert!(replay_repeat_extra(&dir,&outcome,&original,&owned,&mut original.clone()).is_err());
        fs::remove_dir_all(dir).unwrap();
    }

    fn test_dir(name: &str) -> PathBuf {
        let nanos = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        std::env::temp_dir().join(format!("orgraft_{name}_{nanos}"))
    }
}

// Source replay uses the same ordinary ALT selection as checkpoint2.
pub(crate) fn replay_ordinary_source_edits(high: &Path, saved: &Path, sequence: &str) -> Result<String, OrgraftError> {
    let expected=read_high_variant_edits(high)?;
    let actual=read_pos_ref_alt(saved)?;
    let tuples=|v: &[VariantEdit]| v.iter().map(|e|(e.pos,e.reference.clone(),e.alternate.clone())).collect::<Vec<_>>();
    if tuples(&expected)!=tuples(&actual) || expected.is_empty() { return Err(OrgraftError::InvalidArgument("ordinary source saved edits differ from native high-variant selection".into())); }
    for e in &actual {
        let start=e.pos.checked_sub(1).ok_or_else(||OrgraftError::InvalidArgument("ordinary source REF coordinate is zero".into()))?;
        if e.reference.is_empty() || sequence.get(start..start+e.reference.len()) != Some(e.reference.as_str()) { return Err(OrgraftError::InvalidArgument("ordinary source edit REF differs from complete input".into())); }
    }
    apply_edits_to_sequence(sequence,&actual)
}
