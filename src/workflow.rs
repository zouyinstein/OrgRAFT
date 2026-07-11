use std::collections::{BTreeMap, BTreeSet};
use std::fmt::Write as FmtWrite;
use std::fs::{self, File};
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};
use std::process::Command;

use crate::commands;
use crate::error::OrgraftError;
use crate::sv_repair::{repair_sv_subgroup, select_sv_subgroup_spec, SvRepairRequest};
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
    --sample NAME       project sample [sample_001]
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
  correct             apply pos/ref/alt table to one FASTA record
  test-correction     smoke test FASTA correction with explicit inputs
  test-fake-validate  simulate validate failure from swapped polish_aln FASTA

Common options:
  --config FILE       workflow config [orgraft.workflow.toml]
  --case NAME         select one workflow.case section
  --round N           checkpoint2 validation/correction round [1]
  --sv-subgroup SPEC  manually repair one group_name:old_index at checkpoint2
  --force             overwrite generated checkpoint/correction outputs

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
        command: "rebuild",
        responsibility: "rebuild final verified graph and compact reports",
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
    rebuild_out_dir: PathBuf,
    rebuild_edited_gfa: Option<PathBuf>,
    rebuild_polished_fasta: Option<PathBuf>,
    image_reference_fasta: Option<PathBuf>,
    workflow_dir: PathBuf,
    linearized_fasta: Option<PathBuf>,
    polish_reference: Option<PathBuf>,
    sv_correction_subgroup: Option<String>,
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
    ManualRequired,
}

#[derive(Debug, Clone)]
struct VariantEdit {
    pos: usize,
    reference: String,
    alternate: String,
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
        Ok(Self {
            sample: option_value(args, "--sample")?
                .unwrap_or("sample_001")
                .to_string(),
            results_dir: option_value(args, "--results-dir")?
                .unwrap_or("results_workflow")
                .to_string(),
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

    fs::write(out_path, workflow_config_template(&options))?;
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

[commands.rebuild]
enabled = true
out_dir = "${{results_dir}}/05.rebuild"
threads = 16
# The generated workflow script passes the final checkpoint2 polished FASTA to
# rebuild automatically. Set edited_gfa/polished_fasta only for standalone
# rebuild overrides.
# image_reference_fasta enables optional PDF/SVG export through GFA_Editor; export
# failures are recorded in the rebuild run report and do not fail core GFA/FASTA output.
# merged_gfa_template = "${{results_dir}}/03.resolve_gfa/mito/graph/merged_unresolved.gfa"
# keep_debug = true

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
# 04.polish writes validation evidence -> 05.rebuild writes final products.
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

rebuild_out_dir = "${{results_dir}}/05.rebuild/${{organelle}}"
rebuild_edited_gfa = "${{results_dir}}/workflow/${{organelle}}/${{subgraph}}/checkpoint_1/checked_draft.gfa"
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
rebuild_out_dir = "${{results_dir}}/05.rebuild/${{organelle}}"
rebuild_edited_gfa = "${{results_dir}}/workflow/${{organelle}}/${{subgraph}}/checkpoint_1/checked_draft.gfa"
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
    if option_value(args, "--case")?.is_some() || config.cases.len() == 1 {
        let case = select_case(config, args)?;
        let out = option_value(args, "--out")?
            .map(PathBuf::from)
            .unwrap_or_else(|| case.workflow_dir.join("workflow.commands.sh"));
        write_plan_script(&out, config, case, None, false, true)?;
        Ok(out)
    } else {
        let out = option_value(args, "--out")?
            .map(PathBuf::from)
            .unwrap_or_else(|| config.results_dir.join("workflow.commands.sh"));
        write_all_cases_plan_script(&out, config)?;
        Ok(out)
    }
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
        Checkpoint2Status::ManualRequired => {
            println!("checkpoint2 manual_required at round {round}");
        }
    }
    Ok(())
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
            rebuild_out_dir: out_dir.join("rebuild"),
            rebuild_edited_gfa: None,
            rebuild_polished_fasta: None,
            image_reference_fasta: None,
            workflow_dir: out_dir.join("workflow").join(organelle).join(subgraph),
            linearized_fasta: None,
            polish_reference: None,
            sv_correction_subgroup: None,
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
                if config.rebuild.enabled {
                    commands::rebuild::run(&rebuild_args_with_polished(
                        &config,
                        case,
                        force,
                        Some(&polished_aln_path(case, round)),
                    )?)?;
                }
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
        }
    }

    Ok(())
}

fn checkpoint1_impl(
    config: &WorkflowConfig,
    case: &WorkflowCase,
    force: bool,
) -> Result<Checkpoint1Status, OrgraftError> {
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
    if round == 0 {
        return Err(OrgraftError::InvalidArgument(
            "--round must be greater than 0".to_string(),
        ));
    }
    let checkpoint_dir = case
        .workflow_dir
        .join("checkpoint_2")
        .join(format!("round_{round}"));
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

    if config.auto_sv_correction || configured_sv_subgroup.is_some() {
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
    if edits.is_empty() {
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
        return Ok(Checkpoint2Status::Complete);
    }

    if !config.auto_snv_indel_correction {
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

fn write_plan_script(
    path: &Path,
    config: &WorkflowConfig,
    case: &WorkflowCase,
    manual_message: Option<&str>,
    comment_checkpoint1: bool,
    include_recruit: bool,
) -> Result<(), OrgraftError> {
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
    if config.asm.enabled {
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
    writeln!(script, "{}", resolve_command(config, case, config.force)).unwrap();
    writeln!(script).unwrap();
    writeln!(script, "final_polished=\"\"").unwrap();
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
        writeln!(
            script,
            "  {}",
            polish_command(config, case, round, &draft, config.force)
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
        writeln!(
            script,
            "      final_polished={}",
            shell_quote(&polished_aln_path(case, round))
        )
        .unwrap();
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
        "05.rebuild",
        "rebuild the final verified graph and compact reports",
    );
    if config.rebuild.enabled {
        writeln!(
            script,
            "{}",
            rebuild_command_with_polished_arg(config, case, config.force, "\"${final_polished}\"",)
        )
        .unwrap();
    } else {
        writeln!(
            script,
            "# {}",
            rebuild_command_with_polished(config, case, config.force, "FINAL_POLISHED_FASTA",)
        )
        .unwrap();
    }

    write_executable_text(path, &script)
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

fn write_all_cases_plan_script(path: &Path, config: &WorkflowConfig) -> Result<(), OrgraftError> {
    let enabled = enabled_cases(config);
    if enabled.is_empty() {
        return Err(OrgraftError::InvalidArgument(
            "no enabled workflow cases are configured".to_string(),
        ));
    }
    for case in &enabled {
        write_plan_script(
            &case.workflow_dir.join("workflow.commands.sh"),
            config,
            case,
            None,
            false,
            false,
        )?;
    }

    let mut script = String::new();
    writeln!(script, "#!/usr/bin/env bash").unwrap();
    writeln!(script, "set -euo pipefail").unwrap();
    writeln!(script).unwrap();
    writeln!(
        script,
        "# Generated by orgraft workflow for all configured cases"
    )
    .unwrap();
    writeln!(script, "# sample: {}", config.sample).unwrap();
    writeln!(script, "# mode: {}", config.mode).unwrap();
    writeln!(script, "# results_dir: {}", config.results_dir.display()).unwrap();
    writeln!(script).unwrap();
    write_script_runtime_header(&mut script);

    writeln!(
        script,
        "# 01.recruit is shared across enabled cases and runs at most once here."
    )
    .unwrap();
    writeln!(
        script,
        "# Each case script then expands case-local stages 02-05 in order."
    )
    .unwrap();
    writeln!(script).unwrap();
    write_master_recruit_stage(&mut script, config, &enabled)?;
    for case in enabled {
        let case_script = case.workflow_dir.join("workflow.commands.sh");
        writeln!(
            script,
            "echo {}",
            shell_quote_str(&format!(
                "==> running case {} ({}/{})",
                case.name, case.organelle, case.subgraph
            ))
        )
        .unwrap();
        writeln!(script, "bash {}", shell_quote(&case_script)).unwrap();
        writeln!(script).unwrap();
    }

    write_executable_text(path, &script)
}

fn write_master_recruit_stage(
    script: &mut String,
    config: &WorkflowConfig,
    cases: &[&WorkflowCase],
) -> Result<(), OrgraftError> {
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
    let edited_gfa = case
        .rebuild_edited_gfa
        .as_ref()
        .or(rebuild.edited_gfa.as_ref())
        .unwrap_or(&case.checked_draft_gfa);
    let image_reference_fasta = case
        .image_reference_fasta
        .as_ref()
        .or(rebuild.image_reference_fasta.as_ref());

    let mut args = vec![
        "--organelle".to_string(),
        case.organelle.clone(),
        "--subgraph".to_string(),
        case.subgraph.clone(),
        "--edited-gfa".to_string(),
        edited_gfa.display().to_string(),
        "--polished-fasta".to_string(),
        polished.display().to_string(),
        "--soft-paths".to_string(),
        config.soft_paths.display().to_string(),
        "--out-dir".to_string(),
        case.rebuild_out_dir.display().to_string(),
    ];
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
    rebuild_command_with_polished_arg(config, case, force, &shell_quote_str(polished_fasta))
}

fn rebuild_command_with_polished_arg(
    config: &WorkflowConfig,
    case: &WorkflowCase,
    force: bool,
    polished_fasta_arg: &str,
) -> String {
    let rebuild = &config.rebuild;
    let edited_gfa = case
        .rebuild_edited_gfa
        .as_ref()
        .or(rebuild.edited_gfa.as_ref())
        .unwrap_or(&case.checked_draft_gfa);
    let out_dir = if case.rebuild_out_dir.as_os_str().is_empty() {
        &rebuild.out_dir
    } else {
        &case.rebuild_out_dir
    };
    let image_reference_fasta = case
        .image_reference_fasta
        .as_ref()
        .or(rebuild.image_reference_fasta.as_ref());

    let mut args = vec![
        orgraft_shell_token().to_string(),
        "rebuild".to_string(),
        "--organelle".to_string(),
        shell_quote_str(&case.organelle),
        "--subgraph".to_string(),
        shell_quote_str(&case.subgraph),
        "--edited-gfa".to_string(),
        shell_quote(edited_gfa),
        "--polished-fasta".to_string(),
        polished_fasta_arg.to_string(),
        "--soft-paths".to_string(),
        shell_quote(&config.soft_paths),
        "--out-dir".to_string(),
        shell_quote(out_dir),
    ];
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
    push_flag(&mut args, "--keep-debug", rebuild.keep_debug);
    push_flag(&mut args, "--force", force);
    args.extend(rebuild.extra_args.iter().cloned());
    args.join(" ")
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

fn load_config_from_args(args: &[String]) -> Result<WorkflowConfig, OrgraftError> {
    let config_path = option_value(args, "--config")?.unwrap_or("orgraft.workflow.toml");
    WorkflowConfig::from_path(Path::new(config_path))
}

impl WorkflowConfig {
    fn from_path(path: &Path) -> Result<Self, OrgraftError> {
        let raw = read_toml_like(path)?;
        let project = raw.section("project");
        let software = raw.section("software");
        let workflow = raw.section("workflow");

        let sample = project
            .and_then(|section| section.get("sample"))
            .cloned()
            .unwrap_or_else(|| "sample_001".to_string());
        let results_dir_string = project
            .and_then(|section| section.get("results_dir"))
            .cloned()
            .unwrap_or_else(|| "results_workflow".to_string());
        let source_results_dir_string = project
            .and_then(|section| section.get("source_results_dir"))
            .cloned()
            .unwrap_or_else(|| results_dir_string.clone());
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
        let rebuild = WorkflowRebuildConfig::from_section(
            raw.section("commands.rebuild"),
            &sample,
            &results_dir_string,
            &source_results_dir_string,
        )?;

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
            )?);
        }
        for case in &mut cases {
            if case.reference.is_none() && case.pre_rotated_reference.is_none() {
                case.reference = recruit_bait_for_organelle(&recruit, &case.organelle);
            }
        }

        Ok(Self {
            config_path: path.to_path_buf(),
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
            rebuild,
            topology_simple_allowed_classes,
            cases,
        })
    }
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
        let rebuild_out_dir = path_or("rebuild_out_dir", "${results_dir}/05.rebuild".to_string());
        let rebuild_edited_gfa = optional_path("rebuild_edited_gfa");
        let rebuild_polished_fasta = optional_path("rebuild_polished_fasta");
        let image_reference_fasta = optional_path("image_reference_fasta");
        let reference = optional_path("reference");
        let pre_rotated_reference = optional_path("pre_rotated_reference");
        let linearized_fasta = optional_path("linearized_fasta");
        let polish_reference = optional_path("polish_reference");
        let sv_correction_subgroup = section.get("sv_correction_subgroup").cloned();

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
            rebuild_out_dir,
            rebuild_edited_gfa,
            rebuild_polished_fasta,
            image_reference_fasta,
            workflow_dir,
            linearized_fasta,
            polish_reference,
            sv_correction_subgroup,
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
        Ok(Self {
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
            extra_args: split_extra_args(section.get("extra_args")),
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
        return config
            .cases
            .iter()
            .find(|case| case.name == name)
            .ok_or_else(|| {
                OrgraftError::InvalidArgument(format!("workflow case `{name}` was not found"))
            });
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
        let checkpoint1 = optional_metric_value(
            &case
                .workflow_dir
                .join("checkpoint_1")
                .join("checkpoint_1.status.tsv"),
            "status",
        )?
        .unwrap_or_else(|| "pending".to_string());
        let checkpoint2 = latest_checkpoint2_status(config, case)?;
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

        for round in 1..=WORKFLOW_HARD_MAX_ROUNDS {
            let polish_report = polish_report_path(case, round);
            if polish_report.exists() {
                rows.push(RuntimeStageRow {
                    scope: case.name.clone(),
                    stage: format!("polish round {round} stage sum"),
                    seconds: sum_polish_stage_seconds(&polish_report)?,
                    source: polish_report.display().to_string(),
                });
            }
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

fn snv_indel_segments_path(case: &WorkflowCase, round: usize) -> PathBuf {
    polish_round_dir(case, round).join("03.validate/01.data/snv_indel_segments.tsv")
}

fn sv_read_index_path(case: &WorkflowCase, round: usize) -> PathBuf {
    polish_round_dir(case, round).join("03.validate/01.data/sv_read_index.tsv")
}

fn polished_aln_path(case: &WorkflowCase, round: usize) -> PathBuf {
    if round == 1 {
        polish_round_dir(case, round).join("02.polish/polished_aln.fasta")
    } else {
        polish_round_dir(case, round)
            .join("01.inputs")
            .join(format!("linear_subgraph.round_{round}.fasta"))
    }
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
    let mut args = vec![
        orgraft_shell_token().to_string(),
        "polish".to_string(),
        "--organelle".to_string(),
        shell_quote_str(&case.organelle),
        "--subgraph".to_string(),
        shell_quote_str(&case.subgraph),
        "--draft".to_string(),
        shell_quote(draft),
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
    fn workflow_starts_with_workflow_config() {
        assert_eq!(DEFAULT_WORKFLOW.first().unwrap().command, "workflow_config");
    }

    #[test]
    fn workflow_ends_with_rebuild() {
        assert_eq!(DEFAULT_WORKFLOW.last().unwrap().command, "rebuild");
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
        assert!(template.contains("rebuild_out_dir = \"${results_dir}/05.rebuild/${organelle}\""));
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
        assert!(default_template.contains("sample = \"sample_001\""));
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
            PathBuf::from("results_workflow/draft/mito.gfa")
        );
        assert_eq!(
            config.cases[0].unitig_graph,
            Some(PathBuf::from("results_workflow/draft/mito.unitig.gfa"))
        );
        assert_eq!(
            config.cases[0].checked_draft_gfa,
            PathBuf::from(
                "results_workflow/workflow/mito/subgraph_001/checkpoint_1/checked_draft.gfa"
            )
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
        assert!(script.contains("manual_required)"));
        assert!(script.contains("final_polished="));
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
        assert!(script.contains("# 05.rebuild:"));
        assert!(script.contains("\"${ORGRAFT_BIN}\" recruit"));
        assert!(script.contains("\"${ORGRAFT_BIN}\" asm"));
        assert!(script.contains("--image-reference-fasta mito.fa"));
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

        let config = WorkflowConfig::from_path(&config_path).unwrap();
        assert_eq!(
            config.cases[0].reference,
            Some(PathBuf::from("refs/plastid.fa"))
        );
        assert_eq!(
            config.cases[0].draft_graph,
            PathBuf::from("results_workflow/02.draft_asm/plastid/03.finalize_graph/graph.gfa")
        );
        assert_eq!(
            config.cases[0].reads,
            PathBuf::from("results_workflow/01.recruit/plastid.fastq.gz")
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
        write_all_cases_plan_script(&master, &config).unwrap();

        let script = fs::read_to_string(&master).unwrap();
        assert!(script.contains("Generated by orgraft workflow for all configured cases"));
        assert!(script.contains("01.recruit is shared across enabled cases"));
        assert_eq!(script.matches("\"${ORGRAFT_BIN}\" recruit").count(), 1);
        assert_eq!(script.matches("\"${ORGRAFT_BIN}\" asm").count(), 0);
        assert!(script.contains("recruit_ready=1"));
        assert!(script.contains("recruit outputs already exist; skip global recruit"));
        assert!(script.contains("bash "));
        assert!(script.contains("mito_subgraph_001"));
        assert!(script.contains("plastid_subgraph_001"));
        let mito_script_path = config.cases[0].workflow_dir.join("workflow.commands.sh");
        let plastid_script_path = config.cases[1].workflow_dir.join("workflow.commands.sh");
        assert!(mito_script_path.exists());
        assert!(plastid_script_path.exists());
        let mito_script = fs::read_to_string(&mito_script_path).unwrap();
        let plastid_script = fs::read_to_string(&plastid_script_path).unwrap();
        assert!(!mito_script.contains("\"${ORGRAFT_BIN}\" recruit"));
        assert!(mito_script.contains("recruit handled once by the master workflow script"));
        assert!(mito_script.contains("\"${ORGRAFT_BIN}\" asm --reads"));
        assert!(mito_script.contains("--image-reference-fasta mito.fa"));
        assert!(!plastid_script.contains("\"${ORGRAFT_BIN}\" recruit"));
        assert!(plastid_script.contains("recruit handled once by the master workflow script"));
        assert!(plastid_script.contains("\"${ORGRAFT_BIN}\" asm --reads"));
        assert!(plastid_script.contains("--image-reference-fasta plastid.fa"));
        assert!(mito_script.contains("--organelle mito"));
        assert!(plastid_script.contains("--organelle plastid"));
        assert!(config.cases[0]
            .workflow_dir
            .join("workflow.commands.sh")
            .exists());
        assert!(config.cases[1]
            .workflow_dir
            .join("workflow.commands.sh")
            .exists());
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
        write_all_cases_plan_script(&master, &config).unwrap();

        let script = fs::read_to_string(&master).unwrap();
        assert!(!script.contains("mito_subgraph_001"));
        assert!(script.contains("plastid_subgraph_001"));
        assert!(!config.cases[0]
            .workflow_dir
            .join("workflow.commands.sh")
            .exists());
        assert!(config.cases[1]
            .workflow_dir
            .join("workflow.commands.sh")
            .exists());

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

    fn test_dir(name: &str) -> PathBuf {
        let nanos = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        std::env::temp_dir().join(format!("orgraft_{name}_{nanos}"))
    }
}
