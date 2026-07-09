use std::cmp::{Ordering, Reverse};
use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet};
use std::fmt::Write as FmtWrite;
use std::fs::{self, File, OpenOptions};
use std::io::{BufRead, BufReader, BufWriter, Read, Write};
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};
use std::sync::{mpsc, Arc, Mutex};
use std::thread;
use std::time::{Instant, SystemTime, UNIX_EPOCH};

use crate::commands::shared::{print_contract, CommandContract};
use crate::error::OrgraftError;

const HELP: &str = r#"orgraft polish

Create a per-subgraph polish workspace and track polish/evaluation stages.

Usage:
  orgraft polish --organelle NAME --draft FILE --reference FILE --reads FILE [options]

Inputs:
  --organelle NAME                        organelle name for this polish run
  --subgraph ID                           subgraph/ring id [subgraph_001]
  --draft FILE                            linearized subgraph FASTA
  --reference FILE                        rotated linear reference FASTA
  --reads FILE                            reads already binned to this subgraph
  --soft-paths FILE                       tool paths file [soft_paths.txt]

Outputs:
  --out-dir DIR                           polish root output directory [results/polish]
  --force                                 replace this organelle/subgraph polish workspace

Additional Parameters:
  --threads N                             worker threads [8]
  --max-rounds N                          planned validation rounds [3]
  --validate-round N                      workflow validation round directory [1]
  --per-read-variant-calls on|off         run read-level SV/SNP-InDel evidence stages [on]
  --snv-indel-overlap-policy MODE         mark-overlap|mask-both|assign-downstream [mark-overlap]
  --plot-help                             show advanced plot tuning options

Layout: OUT/NAME/SUBGRAPH/round_N/{01.inputs,02.polish,03.validate,logs}
"#;

const PLOT_HELP: &str = r#"orgraft polish plot parameters

Common:
  Plots are attempted automatically when validation data exists; scripts are always written.
  --plot-range START-END                  restrict plots to a reference interval; full length when omitted
  --plot-dpi N                            plot raster DPI [300]
  --plot-output-format png|pdf|both       plot file format [png]

SV Coverage:
  --coverage-plot-rasterize on|off        rasterize dense coverage artists in PDF [on]
  --sv-plot-highlight-subgroups SPEC      comma-separated group:old_index specs for SV coverage highlight
  --sv-plot-highlight-read-ids FILE       read-id file for SV coverage highlight
  --sv-plot-highlight-min-fraction FLOAT  auto-highlight min FL fraction for non-reference subgroups [0.005]
  --sv-plot-highlight-min-reads N         auto-highlight min read count for non-reference subgroups [10]

SNV/InDel:
  --snv-indel-plot-rasterize on|off       rasterize dense SNV/InDel scatter artists in PDF [on]
  --snv-indel-plot-low-confidence TYPES   non-high, none, or comma-separated confidence labels [non-high]
  --snv-indel-plot-low-min-reads N        grey SNV/InDel points below this read count [3]
  --snv-indel-plot-low-min-fraction FLOAT grey SNV/InDel points below this frequency [0]
  --snv-indel-plot-high-risk-fraction F   orange/red highlight threshold for SNV/InDel points [0.5]

Generated scripts in round_N/03.validate/02.plots also provide their own --help for replotting.
"#;

const DEFAULT_SUBGRAPH: &str = "subgraph_001";
const DEFAULT_OUT_DIR: &str = "results/polish";
const DEFAULT_SOFT_PATHS: &str = "soft_paths.txt";
const SV_MINIMAP2_OPTIONS: &[&str] = &[
    "-x",
    "map-hifi",
    "-c",
    "-k",
    "11",
    "-w",
    "7",
    "--secondary=yes",
    "-N",
    "200",
    "-p",
    "0.01",
];
const TERMINAL_EXTENSION_WINDOW: usize = 8;
const TERMINAL_EXTENSION_MAX_GAP: usize = 2;
const TERMINAL_EXTENSION_MIN_ALIGNMENT_LENGTH: usize = 500;
const FL_PERCENT_TOTAL_THRESHOLD: f64 = 98.0;
const READ_GROUP_WATCH_FL_MIN: usize = 20;
const READ_SUBGROUP_WATCH_MIN: usize = 10;
const REFERENCE_SUPPORT_REP_MID_OLP_MIN: f64 = 1000.0;
const SV_SUPPORT_WINDOW_BP: usize = 1000;
const SV_SUPPORT_MIN_GREEN_FRACTION: f64 = 0.50;
const SV_SUPPORT_MAX_LOW_GREEN_WINDOW_FRACTION: f64 = 0.05;
const SV_SUPPORT_LOW_GREEN_FRACTION: f64 = 0.20;
const SV_SUPPORT_MIN_GREEN_DEPTH: f64 = 3.0;
const HIGH_SUBGROUP_MIN_FRACTION: f64 = 0.005;
const DEFAULT_HIGHLIGHT_MIN_FRACTION: f64 = 0.005;
const DEFAULT_HIGHLIGHT_MIN_READS: usize = READ_SUBGROUP_WATCH_MIN;
const DEFAULT_SNV_INDEL_PLOT_LOW_CONFIDENCE: &str = "non-high";
const DEFAULT_SNV_INDEL_PLOT_LOW_MIN_READS: usize = 3;
const DEFAULT_SNV_INDEL_PLOT_LOW_MIN_FRACTION: f64 = 0.0;
const DEFAULT_SNV_INDEL_PLOT_HIGH_RISK_FRACTION: f64 = 0.5;
const DEFAULT_PLOT_DPI: usize = 300;
const DEFAULT_COVERAGE_PLOT_RASTERIZE: bool = true;
const DEFAULT_SNV_INDEL_PLOT_RASTERIZE: bool = true;
const DEFAULT_PLOT_OUTPUT_FORMAT: PlotOutputFormat = PlotOutputFormat::Png;
const SNV_INDEL_SV_CONTEXT_FILTER: bool = true;
const BREAKPOINT_WINDOW_BP: usize = 1000;
const SNV_INDEL_SHORT_TERMINAL_SEGMENT_BP: usize = 1000;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum SnvIndelOverlapPolicy {
    AssignDownstream,
    MaskBoth,
    MarkOverlap,
}

impl SnvIndelOverlapPolicy {
    fn parse(value: &str, option: &str) -> Result<Self, OrgraftError> {
        match value {
            "assign-downstream" | "assign_downstream" | "downstream" => {
                Ok(Self::AssignDownstream)
            }
            "mask-both" | "mask_both" | "both" => Ok(Self::MaskBoth),
            "mark-overlap" | "mark_overlap" | "mark" => Ok(Self::MarkOverlap),
            _ => Err(OrgraftError::InvalidArgument(format!(
                "{option} must be `assign-downstream`, `mask-both`, or `mark-overlap`, got `{value}`"
            ))),
        }
    }

    fn as_str(self) -> &'static str {
        match self {
            Self::AssignDownstream => "assign_downstream",
            Self::MaskBoth => "mask_both",
            Self::MarkOverlap => "mark_overlap",
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum PlotOutputFormat {
    Png,
    Pdf,
    Both,
}

impl PlotOutputFormat {
    fn parse(value: &str, option: &str) -> Result<Self, OrgraftError> {
        match value {
            "png" => Ok(Self::Png),
            "pdf" => Ok(Self::Pdf),
            "both" | "all" | "png+pdf" | "pdf+png" => Ok(Self::Both),
            _ => Err(OrgraftError::InvalidArgument(format!(
                "{option} must be `png`, `pdf`, or `both`, got `{value}`"
            ))),
        }
    }

    fn as_str(self) -> &'static str {
        match self {
            Self::Png => "png",
            Self::Pdf => "pdf",
            Self::Both => "both",
        }
    }
}

pub fn run(args: &[String]) -> Result<(), OrgraftError> {
    if args.is_empty() {
        println!("{HELP}");
        return Ok(());
    }

    if args.iter().any(|arg| arg == "-h" || arg == "--help") {
        println!("{HELP}");
        return Ok(());
    }

    if args.iter().any(|arg| arg == "--plot-help") {
        println!("{PLOT_HELP}");
        return Ok(());
    }

    if args.iter().any(|arg| arg == "--contract") {
        print_contract(&contract());
        return Ok(());
    }

    let options = PolishOptions::from_args(args)?;
    run_polish_scaffold(&options)
}

fn run_polish_scaffold(options: &PolishOptions) -> Result<(), OrgraftError> {
    let started = Instant::now();
    let inputs = ResolvedInputs::from_options(options)?;
    let paths = PolishPaths::new(
        &options.out_dir,
        &options.organelle,
        &options.subgraph,
        options.validate_round,
    );

    if paths.round_dir.exists() {
        if options.force {
            fs::remove_dir_all(&paths.round_dir)?;
        } else {
            return Err(OrgraftError::InvalidArgument(format!(
                "{} already exists; use --force to replace this workflow round workspace",
                paths.round_dir.display()
            )));
        }
    }

    paths.create_dirs()?;
    if fs::symlink_metadata(&paths.input_reads).is_ok() {
        fs::remove_file(&paths.input_reads)?;
    }
    let extracted_draft = extract_fasta_record_by_id(
        &inputs.draft,
        &paths.input_draft,
        &options.subgraph,
        "--draft",
    )?;
    let reference_id = resolve_reference_id(&inputs, &options.subgraph, &extracted_draft)?;
    let extracted_reference = extract_fasta_record_by_id(
        &inputs.reference,
        &paths.input_reference,
        &reference_id,
        "--reference",
    )?;
    link_or_copy(&inputs.reads, &paths.input_reads)?;
    let prepare_seconds = started.elapsed().as_secs_f64();

    let mut command_records = Vec::new();
    let validation_round = format!("round_{}", options.validate_round);
    let mut stage_records = vec![StageRecord::ok(
        "prepare",
        "setup",
        prepare_seconds,
        "polish workspace initialized",
    )];

    let (round_reports, alignment_report) = if options.validate_round == 1 {
        let polish_started = Instant::now();
        let round_reports =
            run_minimap2_rust_polish(options, &inputs, &paths, &mut command_records)?;
        let polish_seconds = polish_started.elapsed().as_secs_f64();

        let align_started = Instant::now();
        let alignment_report = align_polished_to_reference(&inputs, &paths, &mut command_records)?;
        let align_seconds = align_started.elapsed().as_secs_f64();

        stage_records.push(StageRecord::ok(
            "polish",
            "setup",
            polish_seconds,
            "minimap2 pileup polish completed",
        ));
        stage_records.push(StageRecord::ok(
            "align-reference",
            "setup",
            align_seconds,
            "polished FASTA aligned to rotated reference coordinates",
        ));
        (round_reports, Some(alignment_report))
    } else {
        stage_records.push(StageRecord::skipped(
            "polish",
            validation_round.clone(),
            "validation-only workflow round; using corrected checkpoint FASTA directly",
        ));
        stage_records.push(StageRecord::skipped(
            "align-reference",
            validation_round.clone(),
            "validation-only workflow round; input is already in reference-aligned coordinates",
        ));
        (Vec::new(), None)
    };
    cleanup_input_sidecars(&paths)?;

    let sv_eval_started = Instant::now();
    let sv_eval_report = if options.per_read_variant_calls {
        let report = run_sv_eval_round1(options, &inputs, &paths, &mut command_records)?;
        Some(report)
    } else {
        None
    };
    let sv_eval_seconds = sv_eval_started.elapsed().as_secs_f64();

    let snv_indel_started = Instant::now();
    let snv_indel_report = if options.per_read_variant_calls {
        let sv_eval_report = sv_eval_report.as_ref().ok_or_else(|| {
            OrgraftError::InvalidArgument(
                "SNP/InDel evaluation requires SV read evidence from the same round".to_string(),
            )
        })?;
        let report = run_snv_indel_eval_round1(
            options,
            &inputs,
            &paths,
            sv_eval_report.alignment_records.as_slice(),
            &sv_eval_report.read_sequences,
            &mut command_records,
        )?;
        Some(report)
    } else {
        None
    };
    let snv_indel_seconds = snv_indel_started.elapsed().as_secs_f64();

    let plot_started = Instant::now();
    let plot_stage = if options.per_read_variant_calls {
        run_plot_script(
            options,
            &inputs,
            &paths,
            sv_eval_report.as_ref(),
            snv_indel_report.as_ref(),
            &mut command_records,
        )?;
        Some(StageRecord::ok(
            "plot",
            validation_round.clone(),
            plot_started.elapsed().as_secs_f64(),
            "SV support and SNV/InDel plots generated when Python was available",
        ))
    } else {
        Some(StageRecord::skipped(
            "plot",
            validation_round.clone(),
            "disabled because --per-read-variant-calls off",
        ))
    };

    if options.per_read_variant_calls {
        stage_records.push(StageRecord::ok(
            "sv-eval",
            validation_round.clone(),
            sv_eval_seconds,
            "read-level minimap2 SV evidence written",
        ));
        stage_records.push(StageRecord::ok(
            "snv-indel-eval",
            validation_round.clone(),
            snv_indel_seconds,
            "per-read SNP/InDel calls written",
        ));
        if let Some(record) = plot_stage {
            stage_records.push(record);
        }
    } else {
        stage_records.push(StageRecord::skipped(
            "sv-eval",
            validation_round.clone(),
            "disabled by --per-read-variant-calls off",
        ));
        stage_records.push(StageRecord::skipped(
            "snv-indel-eval",
            validation_round,
            "disabled by --per-read-variant-calls off",
        ));
        if let Some(record) = plot_stage {
            stage_records.push(record);
        }
    }
    write_report(
        &paths.report,
        options,
        &inputs,
        &paths,
        &extracted_draft,
        &extracted_reference,
        &stage_records,
        &round_reports,
        alignment_report.as_ref(),
        sv_eval_report.as_ref(),
        snv_indel_report.as_ref(),
        &command_records,
    )?;

    println!("Wrote {}", paths.round_dir.display());
    if options.validate_round == 1 {
        println!("Wrote {}", paths.polished_fasta.display());
        println!("Wrote {}", paths.aligned_fasta.display());
    }
    println!("Wrote {}", paths.report.display());
    Ok(())
}

fn contract() -> CommandContract {
    CommandContract {
        command: "polish",
        origin: "high-quality graph generation after orgraft resolve",
        purpose: "polish one linearized subgraph FASTA and prepare per-round SV/SNP-InDel evaluation outputs",
        inputs: &[
            "--organelle NAME",
            "--subgraph ID for one resolved graph/ring",
            "--draft FILE linearized subgraph FASTA",
            "--reference FILE rotated linear reference FASTA",
            "--reads FILE reads already binned to this subgraph",
            "soft_paths.txt containing minimap2, blastn, and optionally python for plotting",
        ],
        outputs: &[
            "ORGANELLE/SUBGRAPH/round_N/01.inputs with extracted draft/reference FASTA files and linked reads",
            "ORGANELLE/SUBGRAPH/round_1/02.polish/polished.fasta final polished sequence",
            "ORGANELLE/SUBGRAPH/round_1/02.polish/polished_aln.fasta reference-aligned polished sequence",
            "ORGANELLE/SUBGRAPH/round_N/logs/report.tsv for metadata, status, round metrics, and commands",
            "ORGANELLE/SUBGRAPH/round_N/logs/external.stderr.log for combined minimap2/blastn stderr",
            "ORGANELLE/SUBGRAPH/round_N/03.validate/{01.data,02.plots,03.reports}",
            "ORGANELLE/SUBGRAPH/round_N/03.validate/02.plots/plot_sv_support.py and plot_snv_indel.py for optional plotting",
            "ORGANELLE/SUBGRAPH/round_N/03.validate/01.data/snv_indel_calls.tsv with read-level SNP/InDel evidence",
        ],
        notes: &[
            "one command invocation owns one organelle subgraph; batch multi-subgraph runs can call this command repeatedly",
            "SV failures stop for manual judgement; SNP/InDel correction may continue through evaluation rounds",
            "Rust owns state, reporting, command logging, and compact tables before deeper algorithm rewrites",
            "polish uses minimap2-only linear pileup consensus",
            "normalized SV and SNP/InDel TSVs are the default validation outputs",
            "SNP/InDel evaluation uses the Rust CIGAR-diff caller by default; multi-segment FL read overlaps are marked as lower-confidence evidence unless an alternate overlap policy is selected",
        ],
    }
}

#[derive(Debug, Clone)]
struct PolishOptions {
    organelle: String,
    subgraph: String,
    draft: Option<PathBuf>,
    reference: Option<PathBuf>,
    reads: PathBuf,
    soft_paths: PathBuf,
    out_dir: PathBuf,
    threads: usize,
    max_rounds: usize,
    validate_round: usize,
    force: bool,
    per_read_variant_calls: bool,
    snv_indel_overlap_policy: SnvIndelOverlapPolicy,
    plot_range: Option<PlotRange>,
    plot_dpi: usize,
    plot_output_format: PlotOutputFormat,
    coverage_plot_rasterize: bool,
    snv_indel_plot_rasterize: bool,
    sv_plot_highlight_subgroups: Vec<String>,
    sv_plot_highlight_read_ids: Option<PathBuf>,
    sv_plot_highlight_min_fraction: f64,
    sv_plot_highlight_min_reads: usize,
    snv_indel_plot_low_confidence: String,
    snv_indel_plot_low_min_reads: usize,
    snv_indel_plot_low_min_fraction: f64,
    snv_indel_plot_high_risk_fraction: f64,
}

impl PolishOptions {
    fn from_args(args: &[String]) -> Result<Self, OrgraftError> {
        let mut organelle = None;
        let mut subgraph = DEFAULT_SUBGRAPH.to_string();
        let mut draft = None;
        let mut reference = None;
        let mut reads = None;
        let mut soft_paths = PathBuf::from(DEFAULT_SOFT_PATHS);
        let mut out_dir = PathBuf::from(DEFAULT_OUT_DIR);
        let mut threads = 8usize;
        let mut max_rounds = 3usize;
        let mut validate_round = 1usize;
        let mut force = false;
        let mut per_read_variant_calls = true;
        let mut snv_indel_overlap_policy = SnvIndelOverlapPolicy::MarkOverlap;
        let mut plot_range = None;
        let mut plot_dpi = DEFAULT_PLOT_DPI;
        let mut plot_output_format = DEFAULT_PLOT_OUTPUT_FORMAT;
        let mut coverage_plot_rasterize = DEFAULT_COVERAGE_PLOT_RASTERIZE;
        let mut snv_indel_plot_rasterize = DEFAULT_SNV_INDEL_PLOT_RASTERIZE;
        let mut sv_plot_highlight_subgroups = Vec::new();
        let mut sv_plot_highlight_read_ids = None;
        let mut sv_plot_highlight_min_fraction = DEFAULT_HIGHLIGHT_MIN_FRACTION;
        let mut sv_plot_highlight_min_reads = DEFAULT_HIGHLIGHT_MIN_READS;
        let mut snv_indel_plot_low_confidence = DEFAULT_SNV_INDEL_PLOT_LOW_CONFIDENCE.to_string();
        let mut snv_indel_plot_low_min_reads = DEFAULT_SNV_INDEL_PLOT_LOW_MIN_READS;
        let mut snv_indel_plot_low_min_fraction = DEFAULT_SNV_INDEL_PLOT_LOW_MIN_FRACTION;
        let mut snv_indel_plot_high_risk_fraction = DEFAULT_SNV_INDEL_PLOT_HIGH_RISK_FRACTION;

        let mut index = 0usize;
        while index < args.len() {
            let arg = &args[index];
            match arg.as_str() {
                "--organelle" => {
                    organelle = Some(parse_label(required_value(args, &mut index, arg)?, arg)?);
                }
                "--subgraph" | "--subgraph-id" => {
                    subgraph = parse_label(required_value(args, &mut index, arg)?, arg)?;
                }
                "--draft" | "--linear-subgraph" | "--subgraph-fasta" => {
                    draft = Some(PathBuf::from(required_value(args, &mut index, arg)?));
                }
                "--reference" | "--rotated-reference" => {
                    reference = Some(PathBuf::from(required_value(args, &mut index, arg)?));
                }
                "--reads" | "--subgraph-reads" => {
                    reads = Some(PathBuf::from(required_value(args, &mut index, arg)?));
                }
                "--soft-paths" => {
                    soft_paths = PathBuf::from(required_value(args, &mut index, arg)?);
                }
                "--out-dir" => {
                    out_dir = PathBuf::from(required_value(args, &mut index, arg)?);
                }
                "--threads" => {
                    threads = parse_usize(required_value(args, &mut index, arg)?, arg)?;
                }
                "--max-rounds" => {
                    max_rounds = parse_usize(required_value(args, &mut index, arg)?, arg)?;
                }
                "--validate-round" | "--workflow-round" => {
                    validate_round = parse_usize(required_value(args, &mut index, arg)?, arg)?;
                    if validate_round == 0 {
                        return Err(OrgraftError::InvalidArgument(
                            "--validate-round must be at least 1".to_string(),
                        ));
                    }
                }
                "--per-read-variant-calls" => {
                    per_read_variant_calls =
                        parse_on_off(required_value(args, &mut index, arg)?, arg)?;
                }
                "--snv-indel-overlap-policy"
                | "--snv_indel-overlap-policy"
                | "--snv_indel_overlap_policy"
                | "--variant-overlap-policy"
                | "--variant_overlap_policy" => {
                    snv_indel_overlap_policy =
                        SnvIndelOverlapPolicy::parse(required_value(args, &mut index, arg)?, arg)?;
                }
                "--plot-range" | "--plot-region" => {
                    plot_range = Some(parse_plot_range(
                        required_value(args, &mut index, arg)?,
                        arg,
                    )?);
                }
                "--plot-dpi" => {
                    plot_dpi = parse_usize(required_value(args, &mut index, arg)?, arg)?;
                    if plot_dpi == 0 {
                        return Err(OrgraftError::InvalidArgument(
                            "--plot-dpi must be at least 1".to_string(),
                        ));
                    }
                }
                "--plot-output-format" | "--plot-format" => {
                    plot_output_format =
                        PlotOutputFormat::parse(required_value(args, &mut index, arg)?, arg)?;
                }
                "--coverage-plot-rasterize" => {
                    coverage_plot_rasterize =
                        parse_on_off(required_value(args, &mut index, arg)?, arg)?;
                }
                "--snv-indel-plot-rasterize" => {
                    snv_indel_plot_rasterize =
                        parse_on_off(required_value(args, &mut index, arg)?, arg)?;
                }
                "--sv-plot-highlight-subgroups" => {
                    sv_plot_highlight_subgroups =
                        parse_comma_list(required_value(args, &mut index, arg)?);
                }
                "--sv-plot-highlight-read-ids"
                | "--sv-plot-highlight-read-id-file"
                | "--sv-plot-highlight-read_id_file" => {
                    sv_plot_highlight_read_ids =
                        Some(PathBuf::from(required_value(args, &mut index, arg)?));
                }
                "--sv-plot-highlight-min-fraction" => {
                    sv_plot_highlight_min_fraction =
                        parse_fraction(required_value(args, &mut index, arg)?, arg)?;
                }
                "--sv-plot-highlight-min-reads" | "--sv-plot-highlight-min-count" => {
                    sv_plot_highlight_min_reads =
                        parse_usize(required_value(args, &mut index, arg)?, arg)?;
                }
                "--snv-indel-plot-low-confidence"
                | "--snv-indel-plot-low-confidence-types"
                | "--variant-plot-low-confidence"
                | "--variant-plot-low-confidence-types" => {
                    snv_indel_plot_low_confidence =
                        required_value(args, &mut index, arg)?.to_string();
                    if snv_indel_plot_low_confidence.trim().is_empty() {
                        return Err(OrgraftError::InvalidArgument(format!(
                            "{arg} must be non-empty"
                        )));
                    }
                }
                "--snv-indel-plot-low-min-reads" | "--variant-plot-low-min-reads" => {
                    snv_indel_plot_low_min_reads =
                        parse_usize(required_value(args, &mut index, arg)?, arg)?;
                }
                "--snv-indel-plot-low-min-fraction" | "--variant-plot-low-min-fraction" => {
                    snv_indel_plot_low_min_fraction =
                        parse_fraction(required_value(args, &mut index, arg)?, arg)?;
                }
                "--snv-indel-plot-high-risk-fraction" | "--variant-plot-high-risk-fraction" => {
                    snv_indel_plot_high_risk_fraction =
                        parse_fraction(required_value(args, &mut index, arg)?, arg)?;
                }
                "--force" => force = true,
                other => {
                    return Err(OrgraftError::InvalidArgument(format!(
                        "unknown orgraft polish option `{other}`"
                    )));
                }
            }
            index += 1;
        }

        let organelle = organelle
            .ok_or_else(|| OrgraftError::InvalidArgument("missing --organelle NAME".to_string()))?;
        let reads = reads
            .ok_or_else(|| OrgraftError::InvalidArgument("missing --reads FILE".to_string()))?;

        Ok(Self {
            organelle,
            subgraph,
            draft,
            reference,
            reads,
            soft_paths,
            out_dir,
            threads: threads.max(1),
            max_rounds: max_rounds.max(1),
            validate_round,
            force,
            per_read_variant_calls,
            snv_indel_overlap_policy,
            plot_range,
            plot_dpi,
            plot_output_format,
            coverage_plot_rasterize,
            snv_indel_plot_rasterize,
            sv_plot_highlight_subgroups,
            sv_plot_highlight_read_ids,
            sv_plot_highlight_min_fraction,
            sv_plot_highlight_min_reads,
            snv_indel_plot_low_confidence,
            snv_indel_plot_low_min_reads,
            snv_indel_plot_low_min_fraction,
            snv_indel_plot_high_risk_fraction,
        })
    }
}

#[derive(Debug, Clone)]
struct PlotRange {
    start: usize,
    end: usize,
}

impl PlotRange {
    fn as_arg(&self) -> String {
        format!("{}-{}", self.start, self.end)
    }
}

#[derive(Debug, Clone)]
struct ResolvedInputs {
    draft: PathBuf,
    reference: PathBuf,
    reads: PathBuf,
    soft_paths: PathBuf,
}

impl ResolvedInputs {
    fn from_options(options: &PolishOptions) -> Result<Self, OrgraftError> {
        let draft = match &options.draft {
            Some(path) => path.clone(),
            None => {
                return Err(OrgraftError::InvalidArgument(
                    "missing --draft FILE".to_string(),
                ));
            }
        };
        let reference = match &options.reference {
            Some(path) => path.clone(),
            None => {
                return Err(OrgraftError::InvalidArgument(
                    "missing --reference FILE".to_string(),
                ));
            }
        };
        Ok(Self {
            draft: canonicalize_existing(&draft, "--draft")?,
            reference: canonicalize_existing(&reference, "--reference")?,
            reads: canonicalize_existing(&options.reads, "--reads")?,
            soft_paths: canonicalize_existing(&options.soft_paths, "--soft-paths")?,
        })
    }
}

#[derive(Debug, Clone)]
struct PolishPaths {
    validate_round: usize,
    subgraph_dir: PathBuf,
    round_dir: PathBuf,
    inputs_dir: PathBuf,
    polish_dir: PathBuf,
    eval_dir: PathBuf,
    round1_dir: PathBuf,
    round1_sv_dir: PathBuf,
    round1_sv_data_dir: PathBuf,
    round1_sv_reports_dir: PathBuf,
    round1_snv_indel_dir: PathBuf,
    round1_snv_indel_data_dir: PathBuf,
    round1_snv_indel_reports_dir: PathBuf,
    round1_snv_indel_plots_dir: PathBuf,
    round1_sv_plots_dir: PathBuf,
    logs_dir: PathBuf,
    report: PathBuf,
    external_stderr: PathBuf,
    input_draft: PathBuf,
    input_reference: PathBuf,
    input_reads: PathBuf,
    polished_fasta: PathBuf,
    aligned_fasta: PathBuf,
    round1_sv_whole_read_evidence: PathBuf,
    round1_sv_group_summary: PathBuf,
    round1_sv_subgroup_summary: PathBuf,
    round1_sv_group_ids: PathBuf,
    round1_sv_coverage: PathBuf,
    round1_sv_support_summary: PathBuf,
    round1_sv_high_subgroup_report: PathBuf,
    round1_plot_script: PathBuf,
    round1_snv_indel_per_variant_calls: PathBuf,
    round1_snv_indel_segments: PathBuf,
    round1_snv_indel_variant_type_annotations: PathBuf,
    round1_snv_indel_variant_type_annotations_combined: PathBuf,
    round1_snv_indel_variant_type_annotations_combined_high: PathBuf,
    round1_snv_indel_plot_points: PathBuf,
    round1_snv_indel_plot_script: PathBuf,
}

impl PolishPaths {
    fn new(out_dir: &Path, organelle: &str, subgraph: &str, validate_round: usize) -> Self {
        let subgraph_dir = out_dir.join(organelle).join(subgraph);
        let round_dir = subgraph_dir.join(format!("round_{validate_round}"));
        let inputs_dir = round_dir.join("01.inputs");
        let polish_dir = round_dir.join("02.polish");
        let eval_dir = round_dir.join("03.validate");
        let input_draft_name = if validate_round == 1 {
            "linear_subgraph.fasta".to_string()
        } else {
            format!("linear_subgraph.round_{validate_round}.fasta")
        };
        let round1_dir = eval_dir.clone();
        let round1_data_dir = round1_dir.join("01.data");
        let round1_plots_dir = round1_dir.join("02.plots");
        let round1_reports_dir = round1_dir.join("03.reports");
        let round1_sv_dir = round1_dir.clone();
        let round1_sv_data_dir = round1_data_dir.clone();
        let round1_sv_reports_dir = round1_reports_dir.clone();
        let round1_snv_indel_dir = round1_dir.clone();
        let round1_snv_indel_data_dir = round1_data_dir;
        let round1_snv_indel_reports_dir = round1_reports_dir;
        let round1_snv_indel_plots_dir = round1_plots_dir.clone();
        let round1_sv_plots_dir = round1_plots_dir;
        let logs_dir = round_dir.join("logs");
        Self {
            validate_round,
            round_dir,
            report: logs_dir.join("report.tsv"),
            external_stderr: logs_dir.join("external.stderr.log"),
            input_draft: inputs_dir.join(input_draft_name),
            input_reference: inputs_dir.join("rotated_reference.fasta"),
            input_reads: inputs_dir.join("subgraph_reads.fastq.gz"),
            polished_fasta: polish_dir.join("polished.fasta"),
            aligned_fasta: polish_dir.join("polished_aln.fasta"),
            round1_sv_whole_read_evidence: round1_sv_data_dir.join("sv_read_evidence.tsv"),
            round1_sv_group_summary: round1_sv_reports_dir.join("sv_group_stats.tsv"),
            round1_sv_subgroup_summary: round1_sv_reports_dir.join("sv_subgroup_stats.tsv"),
            round1_sv_group_ids: round1_sv_data_dir.join("sv_read_index.tsv"),
            round1_sv_coverage: round1_sv_data_dir.join("sv_coverage.tsv"),
            round1_sv_support_summary: round1_sv_reports_dir.join("sv_snv_indel_summary.tsv"),
            round1_sv_high_subgroup_report: round1_sv_reports_dir.join("sv_high_subgroups.tsv"),
            round1_plot_script: round1_sv_plots_dir.join("plot_sv_support.py"),
            round1_snv_indel_per_variant_calls: round1_snv_indel_data_dir
                .join("snv_indel_calls.tsv"),
            round1_snv_indel_segments: round1_snv_indel_data_dir.join("snv_indel_segments.tsv"),
            round1_snv_indel_variant_type_annotations: round1_snv_indel_data_dir
                .join("snv_indel_variants.tsv"),
            round1_snv_indel_variant_type_annotations_combined: round1_snv_indel_data_dir
                .join("snv_indel_variants_combined.tsv"),
            round1_snv_indel_variant_type_annotations_combined_high: round1_snv_indel_reports_dir
                .join("snv_indel_high.tsv"),
            round1_snv_indel_plot_points: round1_snv_indel_data_dir
                .join("snv_indel_plot_points.tsv"),
            round1_snv_indel_plot_script: round1_snv_indel_plots_dir.join("plot_snv_indel.py"),
            subgraph_dir,
            inputs_dir,
            polish_dir,
            eval_dir,
            round1_dir,
            round1_sv_dir,
            round1_sv_data_dir,
            round1_sv_reports_dir,
            round1_snv_indel_dir,
            round1_snv_indel_data_dir,
            round1_snv_indel_reports_dir,
            round1_snv_indel_plots_dir,
            round1_sv_plots_dir,
            logs_dir,
        }
    }

    fn create_dirs(&self) -> Result<(), OrgraftError> {
        for path in [
            &self.round_dir,
            &self.inputs_dir,
            &self.polish_dir,
            &self.eval_dir,
            &self.round1_dir,
            &self.round1_sv_dir,
            &self.round1_sv_data_dir,
            &self.round1_sv_reports_dir,
            &self.round1_snv_indel_dir,
            &self.round1_snv_indel_data_dir,
            &self.round1_snv_indel_reports_dir,
            &self.round1_snv_indel_plots_dir,
            &self.round1_sv_plots_dir,
            &self.logs_dir,
        ] {
            fs::create_dir_all(path)?;
        }
        Ok(())
    }

    fn polished_round_fasta(&self, round: usize) -> PathBuf {
        self.polish_dir
            .join(format!("polished_round_{round}.fasta"))
    }

    fn validation_fasta(&self) -> &Path {
        if self.validate_round == 1 {
            &self.aligned_fasta
        } else {
            &self.input_draft
        }
    }
}

#[derive(Debug)]
struct StageRecord {
    stage: String,
    round: String,
    status: &'static str,
    elapsed_seconds: Option<f64>,
    message: String,
}

impl StageRecord {
    fn ok(
        stage: impl Into<String>,
        round: impl Into<String>,
        elapsed_seconds: f64,
        message: impl Into<String>,
    ) -> Self {
        Self {
            stage: stage.into(),
            round: round.into(),
            status: "ok",
            elapsed_seconds: Some(elapsed_seconds),
            message: message.into(),
        }
    }

    fn skipped(
        stage: impl Into<String>,
        round: impl Into<String>,
        message: impl Into<String>,
    ) -> Self {
        Self {
            stage: stage.into(),
            round: round.into(),
            status: "skipped",
            elapsed_seconds: None,
            message: message.into(),
        }
    }
}

#[derive(Debug, Clone)]
struct CommandRecord {
    timestamp: String,
    stage: &'static str,
    round: String,
    status: &'static str,
    elapsed_seconds: f64,
    stdout: String,
    stderr: String,
    command: String,
}

fn run_minimap2_rust_polish(
    options: &PolishOptions,
    inputs: &ResolvedInputs,
    paths: &PolishPaths,
    commands: &mut Vec<CommandRecord>,
) -> Result<Vec<PileupRoundReport>, OrgraftError> {
    let soft_paths = read_soft_paths(&inputs.soft_paths)?;
    let minimap2 = require_tool(&soft_paths, "minimap2")?;
    if paths.external_stderr.exists() {
        fs::remove_file(&paths.external_stderr)?;
    }

    let mut current_reference = paths.input_draft.clone();
    let mut round_reports = Vec::new();
    for round in 1..=2usize {
        let polished_path = if round == 2 {
            paths.polished_fasta.clone()
        } else {
            paths.polished_round_fasta(round)
        };
        let report = run_rust_polish_round(
            round,
            &minimap2,
            options,
            &current_reference,
            &paths.input_reads,
            &polished_path,
            &paths.external_stderr,
            commands,
        )?;
        round_reports.push(report);
        if round < 2 {
            current_reference = polished_path;
        }
    }

    Ok(round_reports)
}

fn run_rust_polish_round(
    round: usize,
    minimap2: &Path,
    options: &PolishOptions,
    reference_path: &Path,
    reads_path: &Path,
    output_path: &Path,
    stderr_path: &Path,
    commands: &mut Vec<CommandRecord>,
) -> Result<PileupRoundReport, OrgraftError> {
    let (record_id, reference) = read_single_fasta_record(reference_path)?;
    let mut pileup = PileupState::new(reference.len());
    let mut command = Command::new(minimap2);
    command
        .arg("-ax")
        .arg("map-hifi")
        .arg("--secondary=no")
        .arg("-t")
        .arg(options.threads.to_string())
        .arg(reference_path)
        .arg(reads_path);
    let command_text = format!("{command:?}");
    let mut stderr_file = OpenOptions::new()
        .create(true)
        .append(true)
        .open(stderr_path)?;
    writeln!(stderr_file, "### round_{round} minimap2 stderr ###")?;
    let stderr_for_child = stderr_file.try_clone()?;
    let started = Instant::now();
    let mut child = command
        .stdout(Stdio::piped())
        .stderr(Stdio::from(stderr_for_child))
        .spawn()?;
    let stdout = child.stdout.take().ok_or_else(|| {
        OrgraftError::InvalidArgument("failed to capture minimap2 stdout".to_string())
    })?;
    let alignments_used_result = pileup.load_sam_reader(BufReader::new(stdout));
    let status = child.wait()?;
    let elapsed_seconds = started.elapsed().as_secs_f64();
    let status_text = if status.success() { "ok" } else { "failed" };
    writeln!(
        OpenOptions::new().append(true).open(stderr_path)?,
        "### round_{round} status={status_text} elapsed_seconds={elapsed_seconds:.3} ###\n"
    )?;
    commands.push(CommandRecord {
        timestamp: timestamp(),
        stage: "polish-align",
        round: format!("round_{round}"),
        status: status_text,
        elapsed_seconds,
        stdout: "stream:minimap2-sam".to_string(),
        stderr: display_path(stderr_path),
        command: command_text,
    });
    if !status.success() {
        return Err(OrgraftError::InvalidArgument(format!(
            "polish-align round_{round} failed; see {}",
            stderr_path.display()
        )));
    }

    let alignments_used = alignments_used_result?;
    let consensus = pileup.consensus(&reference);
    let report = PileupRoundReport {
        round,
        input_length: reference.len(),
        output_length: consensus.len(),
        alignments_used,
        substitutions: pileup.substitutions,
        deletions: pileup.deletions,
        inserted_bases: pileup.inserted_bases,
        low_coverage_bases: pileup.low_coverage_bases,
    };
    write_single_fasta(output_path, &record_id, &consensus)?;
    Ok(report)
}

#[derive(Debug, Clone)]
struct AlignmentReport {
    record_id: String,
    reference_id: String,
    input_length: usize,
    output_length: usize,
    orientation: char,
    reverse_complemented: bool,
    rotation_step: isize,
    best_pident: f64,
    best_length: usize,
    query_start: usize,
    subject_start: usize,
    subject_end: usize,
}

fn align_polished_to_reference(
    inputs: &ResolvedInputs,
    paths: &PolishPaths,
    commands: &mut Vec<CommandRecord>,
) -> Result<AlignmentReport, OrgraftError> {
    let soft_paths = read_soft_paths(&inputs.soft_paths)?;
    let blastn = require_tool(&soft_paths, "blastn")?;
    let (record_id, polished_sequence) = read_single_fasta_record(&paths.polished_fasta)?;
    let (reference_id, _) = read_single_fasta_record(&paths.input_reference)?;

    let initial_blast = temp_file_path("orgraft-polish-align-initial", "tsv");
    run_blastn(
        &blastn,
        &paths.input_reference,
        &paths.polished_fasta,
        &initial_blast,
        &paths.external_stderr,
        "align-reference-blastn",
        "initial",
        commands,
    )?;
    let initial_hits = parse_blast_hits(&initial_blast)?;
    let _ = fs::remove_file(&initial_blast);
    let initial_best = best_hit_for_subject(&initial_hits, &record_id);
    let reverse = initial_best
        .as_ref()
        .is_some_and(|hit| hit.subject_strand() == '-');
    let oriented_sequence = if reverse {
        reverse_complement(&polished_sequence)
    } else {
        polished_sequence.clone()
    };

    let oriented_fasta = temp_file_path("orgraft-polish-align-oriented", "fasta");
    write_single_fasta(&oriented_fasta, &record_id, &oriented_sequence)?;
    let oriented_blast = temp_file_path("orgraft-polish-align-oriented", "tsv");
    run_blastn(
        &blastn,
        &paths.input_reference,
        &oriented_fasta,
        &oriented_blast,
        &paths.external_stderr,
        "align-reference-blastn",
        "oriented",
        commands,
    )?;
    let oriented_hits = parse_blast_hits(&oriented_blast)?;
    let _ = fs::remove_file(&oriented_fasta);
    let _ = fs::remove_file(&oriented_blast);
    let best = best_hit_for_subject(&oriented_hits, &record_id);
    let rotation_step = best
        .as_ref()
        .map(|hit| hit.subject_start as isize - hit.query_start as isize)
        .unwrap_or(0);
    let aligned_sequence = rotate_sequence(&oriented_sequence, rotation_step);
    let orientation = if reverse { '-' } else { '+' };
    let aligned_id = format!(
        "{} [reference={};orientation={};rotation={}]",
        record_id, reference_id, orientation, rotation_step
    );
    write_single_fasta(&paths.aligned_fasta, &aligned_id, &aligned_sequence)?;

    Ok(AlignmentReport {
        record_id,
        reference_id,
        input_length: polished_sequence.len(),
        output_length: aligned_sequence.len(),
        orientation,
        reverse_complemented: reverse,
        rotation_step,
        best_pident: best.as_ref().map(|hit| hit.pident).unwrap_or(0.0),
        best_length: best.as_ref().map(|hit| hit.length).unwrap_or(0),
        query_start: best.as_ref().map(|hit| hit.query_start).unwrap_or(0),
        subject_start: best.as_ref().map(|hit| hit.subject_start).unwrap_or(0),
        subject_end: best.as_ref().map(|hit| hit.subject_end).unwrap_or(0),
    })
}

#[derive(Debug, Clone)]
struct SvEvalReport {
    read_count: usize,
    paf_alignments: usize,
    summary_rows: usize,
    no_alignment_reads: usize,
    whole_read_evidence_rows: usize,
    fl_reads: usize,
    partial_reads: usize,
    reference_support_reads: usize,
    read_group_count: usize,
    read_subgroup_count: usize,
    sv_support_status: String,
    minimap2_mode: &'static str,
    minimap2_workers: usize,
    alignment_records: Vec<AlignmentSummaryRecord>,
    whole_read_evidence_path: PathBuf,
    read_group_summary_path: PathBuf,
    read_subgroup_summary_path: PathBuf,
    read_group_ids_path: PathBuf,
    coverage_path: PathBuf,
    sv_support_summary_path: PathBuf,
    high_subgroup_report_path: PathBuf,
    plot_script_path: PathBuf,
    minimap2_options: String,
    auto_highlight_subgroups: Vec<String>,
    auto_sv_plot_highlight_min_fraction: f64,
    auto_sv_plot_highlight_min_reads: usize,
    read_sequences: HashMap<String, String>,
}

#[derive(Debug, Clone)]
struct SnvIndelReport {
    fl_read_count: usize,
    segment_count: usize,
    total_calls: usize,
    reads_with_calls: usize,
    failed_segments: usize,
    worker_count: usize,
    elapsed_seconds: f64,
    alignment_seconds: f64,
    sum_segment_seconds: f64,
    per_variant_calls_path: PathBuf,
    segments_path: PathBuf,
    variant_type_annotations_path: PathBuf,
    variant_type_annotations_combined_path: PathBuf,
    variant_type_annotations_combined_high_path: PathBuf,
    plot_points_path: PathBuf,
    plot_script_path: PathBuf,
    summary_path: PathBuf,
    reference_path: PathBuf,
    call_mode: &'static str,
    overlap_policy: &'static str,
    sv_context_filter: bool,
    minimap2_preset: &'static str,
    write_timings: SnvIndelWriteTimings,
}

#[derive(Debug, Clone, Default)]
struct SnvIndelWriteTimings {
    per_variant_calls_seconds: f64,
    per_variant_calls_bytes: usize,
    variant_segments_seconds: f64,
    variant_segments_bytes: usize,
    variant_type: VariantTypeAnnotationTimings,
    failed_segments_log_seconds: f64,
    total_seconds: f64,
}

#[derive(Debug, Clone, Default)]
struct VariantTypeAnnotationTimings {
    read_reference_seconds: f64,
    read_frequency_depth_seconds: f64,
    collect_source_rows_seconds: f64,
    write_raw_seconds: f64,
    write_raw_bytes: usize,
    combine_rows_seconds: f64,
    write_combined_seconds: f64,
    write_combined_bytes: usize,
    write_combined_high_seconds: f64,
    write_combined_high_bytes: usize,
    write_plot_points_seconds: f64,
    write_plot_points_bytes: usize,
    append_summary_seconds: f64,
    total_seconds: f64,
    source_rows: usize,
    combined_rows: usize,
}

#[derive(Debug, Clone, Default)]
struct TsvWriteTiming {
    total_seconds: f64,
    bytes: usize,
}

#[derive(Debug, Clone)]
struct ReadGroupReport {
    fl_reads: usize,
    partial_reads: usize,
    reference_support_reads: usize,
    read_group_count: usize,
    read_subgroup_count: usize,
    sv_support_status: String,
    auto_highlight_subgroups: Vec<String>,
}

#[derive(Debug, Clone)]
struct ReadRecord {
    id: String,
    sequence: String,
}

#[derive(Debug, Clone)]
struct PafAlignment {
    query_id: String,
    query_start: usize,
    query_end: usize,
    strand: char,
    target_id: String,
    target_start: usize,
    target_end: usize,
    matches: usize,
    block_len: usize,
    mapq: String,
    alignment_role: String,
    cigar: String,
}

#[derive(Debug, Clone)]
struct BlastLikeAlignment {
    query_id: String,
    subject_id: String,
    pident: f64,
    query_start: usize,
    query_end: usize,
    subject_start: usize,
    subject_end: usize,
}

fn run_sv_eval_round1(
    options: &PolishOptions,
    inputs: &ResolvedInputs,
    paths: &PolishPaths,
    commands: &mut Vec<CommandRecord>,
) -> Result<SvEvalReport, OrgraftError> {
    let soft_paths = read_soft_paths(&inputs.soft_paths)?;
    let minimap2 = require_tool(&soft_paths, "minimap2")?;
    let reads = read_sequence_records(&paths.input_reads)?;
    let reference_by_id = read_fasta_records_by_id(paths.validation_fasta())?;
    let reference_len = reference_by_id
        .values()
        .next()
        .map(|sequence| sequence.len())
        .unwrap_or(0);

    let (paf_by_read, minimap2_workers) = run_sv_minimap2_batch(
        &minimap2,
        options,
        paths.validation_fasta(),
        &paths.input_reads,
        &paths.external_stderr,
        commands,
    )?;

    let mut evidence_output = BufWriter::new(File::create(&paths.round1_sv_whole_read_evidence)?);
    write_whole_read_evidence_header(&mut evidence_output)?;
    let mut paf_alignments = 0usize;
    let mut summary_rows = 0usize;
    let mut no_alignment_reads = 0usize;
    let mut whole_read_evidence_rows = 0usize;
    let mut alignment_records = Vec::with_capacity(reads.len());
    for read in &reads {
        let paf_lines = paf_by_read.get(&read.id).map(Vec::as_slice).unwrap_or(&[]);
        paf_alignments += paf_lines.len();
        let adjusted = extend_paf_terminal_microindels(paf_lines, &read.sequence, &reference_by_id);
        let blast_like = paf_to_blast_like(&adjusted);
        let alignment_summary = build_sorted_alignment_summary(&blast_like, read.sequence.len());
        whole_read_evidence_rows += write_whole_read_evidence_rows(
            &mut evidence_output,
            read,
            paf_lines,
            alignment_summary
                .as_ref()
                .map(|(summary, _record)| summary.as_str()),
        )?;
        if let Some((_summary, record)) = alignment_summary {
            alignment_records.push(record);
            summary_rows += 1;
        } else {
            no_alignment_reads += 1;
        }
    }
    let read_group_report = write_read_group_reports(
        &alignment_records,
        paths,
        reference_len,
        options,
        &paf_by_read,
    )?;

    Ok(SvEvalReport {
        read_count: reads.len(),
        paf_alignments,
        summary_rows,
        no_alignment_reads,
        whole_read_evidence_rows,
        fl_reads: read_group_report.fl_reads,
        partial_reads: read_group_report.partial_reads,
        reference_support_reads: read_group_report.reference_support_reads,
        read_group_count: read_group_report.read_group_count,
        read_subgroup_count: read_group_report.read_subgroup_count,
        sv_support_status: read_group_report.sv_support_status,
        minimap2_mode: "batch",
        minimap2_workers,
        alignment_records,
        whole_read_evidence_path: paths.round1_sv_whole_read_evidence.clone(),
        read_group_summary_path: paths.round1_sv_group_summary.clone(),
        read_subgroup_summary_path: paths.round1_sv_subgroup_summary.clone(),
        read_group_ids_path: paths.round1_sv_group_ids.clone(),
        coverage_path: paths.round1_sv_coverage.clone(),
        sv_support_summary_path: paths.round1_sv_support_summary.clone(),
        high_subgroup_report_path: paths.round1_sv_high_subgroup_report.clone(),
        plot_script_path: paths.round1_plot_script.clone(),
        minimap2_options: SV_MINIMAP2_OPTIONS.join(" "),
        auto_highlight_subgroups: read_group_report.auto_highlight_subgroups,
        auto_sv_plot_highlight_min_fraction: options.sv_plot_highlight_min_fraction,
        auto_sv_plot_highlight_min_reads: options.sv_plot_highlight_min_reads,
        read_sequences: reads
            .into_iter()
            .map(|read| (read.id, read.sequence))
            .collect(),
    })
}

#[derive(Debug, Clone)]
struct ReadIndexMetadata {
    read_class: String,
    group_name: String,
    subgroup_old_index: String,
    subgroup_key: String,
}

#[derive(Debug, Clone)]
struct VariantSegment {
    read_id: String,
    segment_id: String,
    read_class: String,
    group_name: String,
    subgroup_old_index: String,
    subgroup_key: String,
    segment_index: usize,
    segment_count: usize,
    query_start: usize,
    query_end: usize,
    subject_start: isize,
    subject_end: isize,
    strand: char,
    sequence: String,
    trim_note: String,
    overlap_query_intervals: Vec<QueryInterval>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct QueryInterval {
    start: usize,
    end: usize,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum VariantConfidence {
    High,
    OverlapContext,
    ShortSegmentContext,
    ComplexContext,
    TerminalContext,
}

impl VariantConfidence {
    fn as_str(self) -> &'static str {
        match self {
            Self::High => "high",
            Self::OverlapContext => "overlap_context",
            Self::ShortSegmentContext => "short_segment_context",
            Self::ComplexContext => "complex_context",
            Self::TerminalContext => "terminal_context",
        }
    }

    fn reason(self) -> &'static str {
        match self {
            Self::High => ".",
            Self::OverlapContext => "multi_alignment_overlap",
            Self::ShortSegmentContext => "short_terminal_segment",
            Self::ComplexContext => "dense_or_long_indel_neighborhood",
            Self::TerminalContext => "near_segment_alignment_boundary",
        }
    }
}

#[derive(Debug, Clone)]
struct VariantCall {
    pos: usize,
    ref_allele: String,
    alt_allele: String,
    confidence: VariantConfidence,
    query_start: Option<usize>,
    query_end: Option<usize>,
}

impl VariantCall {
    fn high(pos: usize, ref_allele: String, alt_allele: String) -> Self {
        Self {
            pos,
            ref_allele,
            alt_allele,
            confidence: VariantConfidence::High,
            query_start: None,
            query_end: None,
        }
    }

    fn with_query_range(mut self, start: usize, end: usize) -> Self {
        self.query_start = Some(start);
        self.query_end = Some(end.max(start));
        self
    }

    fn with_metadata_from(mut self, other: &VariantCall) -> Self {
        self.confidence = other.confidence;
        self.query_start = other.query_start;
        self.query_end = other.query_end;
        self
    }
}

impl PartialEq for VariantCall {
    fn eq(&self, other: &Self) -> bool {
        self.pos == other.pos
            && self.ref_allele == other.ref_allele
            && self.alt_allele == other.alt_allele
    }
}

impl Eq for VariantCall {}

impl PartialOrd for VariantCall {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for VariantCall {
    fn cmp(&self, other: &Self) -> Ordering {
        self.pos
            .cmp(&other.pos)
            .then_with(|| self.ref_allele.cmp(&other.ref_allele))
            .then_with(|| self.alt_allele.cmp(&other.alt_allele))
    }
}

#[derive(Debug, Clone)]
struct SegmentCallResult {
    segment: VariantSegment,
    calls: Vec<VariantCall>,
    elapsed_seconds: f64,
    exit_status: i32,
    message: String,
    stderr: String,
}

#[derive(Debug, Clone)]
struct SplitSamRecord {
    flag: u16,
    mapq: u8,
    reference_name: String,
    position: usize,
    cigar: String,
    sequence: String,
}

#[derive(Debug, Clone, Copy)]
struct VariantStreamMetrics {
    alignment_seconds: f64,
    split_sam_records: usize,
}

#[derive(Debug, Clone)]
struct SnvIndelRuntimeSummary {
    fl_read_count: usize,
    segment_count: usize,
    total_calls: usize,
    reads_with_calls: usize,
    failed_segments: usize,
    workers: usize,
    elapsed_seconds: f64,
    shared_minimap2_stream_seconds: f64,
    split_sam_records: usize,
    sum_segment_seconds: f64,
    reference_path: String,
    call_mode: &'static str,
    overlap_policy: &'static str,
    sv_context_filter: bool,
    minimap2_preset: &'static str,
}

fn run_snv_indel_eval_round1(
    options: &PolishOptions,
    inputs: &ResolvedInputs,
    paths: &PolishPaths,
    alignment_records: &[AlignmentSummaryRecord],
    read_sequences: &HashMap<String, String>,
    commands: &mut Vec<CommandRecord>,
) -> Result<SnvIndelReport, OrgraftError> {
    let soft_paths = read_soft_paths(&inputs.soft_paths)?;
    let minimap2 = require_tool(&soft_paths, "minimap2")?;

    let read_index = read_index_from_alignment_records(alignment_records);
    let segments = build_variant_segments(
        alignment_records,
        &read_index,
        read_sequences,
        options.snv_indel_overlap_policy,
    )?;
    let segment_fasta = paths
        .logs_dir
        .join("variant_call_segments.round_1.tmp.fasta");
    write_variant_segment_fasta(&segment_fasta, &segments)?;

    let started = Instant::now();
    let (results, worker_count, stream_metrics) = run_custom_variant_caller_segments(
        &minimap2,
        paths.validation_fasta(),
        &segment_fasta,
        &segments,
        options.threads,
        SNV_INDEL_SV_CONTEXT_FILTER,
        &paths.external_stderr,
        commands,
    )?;
    match fs::remove_file(&segment_fasta) {
        Ok(()) => {}
        Err(error) if error.kind() == std::io::ErrorKind::NotFound => {}
        Err(error) => return Err(error.into()),
    }
    let elapsed_seconds = started.elapsed().as_secs_f64();

    let write_total_started = Instant::now();
    let mut write_timings = SnvIndelWriteTimings::default();
    let per_variant_timing = write_per_variant_calls(
        &paths.round1_snv_indel_per_variant_calls,
        &segments,
        &results,
    )?;
    write_timings.per_variant_calls_seconds = per_variant_timing.total_seconds;
    write_timings.per_variant_calls_bytes = per_variant_timing.bytes;
    let segment_timing =
        write_variant_segments_tsv(&paths.round1_snv_indel_segments, &segments, &results)?;
    write_timings.variant_segments_seconds = segment_timing.total_seconds;
    write_timings.variant_segments_bytes = segment_timing.bytes;
    let total_calls = results
        .iter()
        .map(|result| result.calls.len())
        .sum::<usize>();
    let failed_segments = results
        .iter()
        .filter(|result| result.exit_status != 0)
        .count();
    let reads_with_calls = reads_with_calls(&results);
    let sum_segment_seconds = results
        .iter()
        .map(|result| result.elapsed_seconds)
        .sum::<f64>();
    let fl_read_count = segments
        .iter()
        .map(|segment| segment.read_id.as_str())
        .collect::<HashSet<_>>()
        .len();
    let runtime_summary = SnvIndelRuntimeSummary {
        fl_read_count,
        segment_count: segments.len(),
        total_calls,
        reads_with_calls,
        failed_segments,
        workers: worker_count,
        elapsed_seconds,
        shared_minimap2_stream_seconds: stream_metrics.alignment_seconds,
        sum_segment_seconds,
        split_sam_records: stream_metrics.split_sam_records,
        reference_path: display_path(paths.validation_fasta()),
        call_mode: "rust",
        overlap_policy: options.snv_indel_overlap_policy.as_str(),
        sv_context_filter: SNV_INDEL_SV_CONTEXT_FILTER,
        minimap2_preset: "map-hifi",
    };
    write_timings.variant_type = write_variant_type_annotation_tables(
        &paths.round1_snv_indel_variant_type_annotations,
        &paths.round1_snv_indel_variant_type_annotations_combined,
        &paths.round1_snv_indel_variant_type_annotations_combined_high,
        &paths.round1_snv_indel_plot_points,
        &paths.round1_sv_support_summary,
        &runtime_summary,
        paths.validation_fasta(),
        &paths.round1_sv_coverage,
        &segments,
        &results,
    )?;
    let write_started = Instant::now();
    append_snv_failed_segments_to_log(
        &paths.external_stderr,
        &segments,
        &results,
        failed_segments,
    )?;
    write_timings.failed_segments_log_seconds = write_started.elapsed().as_secs_f64();
    write_timings.total_seconds = write_total_started.elapsed().as_secs_f64();
    append_snv_indel_write_timing_summary(&paths.round1_sv_support_summary, &write_timings)?;

    Ok(SnvIndelReport {
        fl_read_count,
        segment_count: segments.len(),
        total_calls,
        reads_with_calls,
        failed_segments,
        worker_count,
        elapsed_seconds,
        alignment_seconds: stream_metrics.alignment_seconds,
        sum_segment_seconds,
        per_variant_calls_path: paths.round1_snv_indel_per_variant_calls.clone(),
        segments_path: paths.round1_snv_indel_segments.clone(),
        variant_type_annotations_path: paths.round1_snv_indel_variant_type_annotations.clone(),
        variant_type_annotations_combined_path: paths
            .round1_snv_indel_variant_type_annotations_combined
            .clone(),
        variant_type_annotations_combined_high_path: paths
            .round1_snv_indel_variant_type_annotations_combined_high
            .clone(),
        plot_points_path: paths.round1_snv_indel_plot_points.clone(),
        plot_script_path: paths.round1_snv_indel_plot_script.clone(),
        summary_path: paths.round1_sv_support_summary.clone(),
        reference_path: paths.validation_fasta().to_path_buf(),
        call_mode: "rust",
        overlap_policy: options.snv_indel_overlap_policy.as_str(),
        sv_context_filter: SNV_INDEL_SV_CONTEXT_FILTER,
        minimap2_preset: "map-hifi",
        write_timings,
    })
}

fn read_index_from_alignment_records(
    records: &[AlignmentSummaryRecord],
) -> HashMap<String, ReadIndexMetadata> {
    let mut subgroups: BTreeMap<String, BTreeMap<SubgroupKey, Vec<usize>>> = BTreeMap::new();
    for (record_index, record) in records.iter().enumerate() {
        if let Some(key) = record.subgroup_key() {
            subgroups
                .entry(record.group_name())
                .or_default()
                .entry(key)
                .or_default()
                .push(record_index);
        }
    }

    let mut subgroup_old_index: BTreeMap<String, BTreeMap<SubgroupKey, usize>> = BTreeMap::new();
    for (group_name, group_subgroups) in &subgroups {
        let mut keys = group_subgroups.keys().cloned().collect::<Vec<_>>();
        keys.sort_by(subgroup_old_index_order);
        subgroup_old_index.insert(
            group_name.clone(),
            keys.into_iter()
                .enumerate()
                .map(|(index, key)| (key, index + 1))
                .collect(),
        );
    }

    records
        .iter()
        .map(|record| {
            let group_name = record.group_name();
            let (subgroup_old_index, subgroup_key) = record
                .subgroup_key()
                .map(|key| {
                    let old_index = subgroup_old_index
                        .get(&group_name)
                        .and_then(|indices| indices.get(&key))
                        .copied()
                        .unwrap_or(0);
                    (old_index.to_string(), key.label())
                })
                .unwrap_or_else(|| (".".to_string(), ".".to_string()));
            (
                record.read_id.clone(),
                ReadIndexMetadata {
                    read_class: if record.is_fl() { "FL" } else { "partial" }.to_string(),
                    group_name,
                    subgroup_old_index,
                    subgroup_key,
                },
            )
        })
        .collect()
}

fn tsv_column(columns: &[&str], name: &str, path: &Path) -> Result<usize, OrgraftError> {
    columns
        .iter()
        .position(|column| *column == name)
        .ok_or_else(|| {
            OrgraftError::InvalidArgument(format!("{} is missing {name} column", path.display()))
        })
}

fn build_variant_segments(
    records: &[AlignmentSummaryRecord],
    read_index: &HashMap<String, ReadIndexMetadata>,
    read_sequences: &HashMap<String, String>,
    overlap_policy: SnvIndelOverlapPolicy,
) -> Result<Vec<VariantSegment>, OrgraftError> {
    let mut segments = Vec::new();
    for record in records {
        let Some(metadata) = read_index.get(&record.read_id) else {
            continue;
        };
        if metadata.read_class != "FL" {
            continue;
        }
        let sequence = read_sequences.get(&record.read_id).ok_or_else(|| {
            OrgraftError::InvalidArgument(format!(
                "read {} is present in SV records but missing from input reads",
                record.read_id
            ))
        })?;
        segments.extend(variant_segments_for_record(
            record,
            metadata,
            sequence,
            overlap_policy,
        )?);
    }
    Ok(segments)
}

fn variant_segments_for_record(
    record: &AlignmentSummaryRecord,
    metadata: &ReadIndexMetadata,
    read_sequence: &str,
    overlap_policy: SnvIndelOverlapPolicy,
) -> Result<Vec<VariantSegment>, OrgraftError> {
    let mut intervals = record
        .alignments
        .iter()
        .enumerate()
        .map(|(index, alignment)| {
            let mut start = alignment.qs.max(1).min(read_sequence.len().max(1));
            let mut end = alignment.qe.max(start).min(read_sequence.len());
            if record.alignments.len() == 1 {
                start = 1;
                end = read_sequence.len();
            }
            (
                index,
                start,
                end,
                alignment,
                "none".to_string(),
                Vec::<QueryInterval>::new(),
            )
        })
        .collect::<Vec<_>>();

    for index in 0..intervals.len().saturating_sub(1) {
        let next_start = intervals[index + 1].1;
        let current_end = intervals[index].2;
        if current_end >= next_start {
            let overlap = QueryInterval {
                start: next_start,
                end: current_end,
            };
            match overlap_policy {
                SnvIndelOverlapPolicy::AssignDownstream => {
                    intervals[index].2 = next_start.saturating_sub(1);
                    intervals[index].4 = format!("trimmed_overlap_before_P{}", index + 2);
                }
                SnvIndelOverlapPolicy::MaskBoth => {
                    intervals[index].2 = next_start.saturating_sub(1);
                    intervals[index + 1].1 = current_end.saturating_add(1);
                    intervals[index].4 = format!("masked_overlap_before_P{}", index + 2);
                    intervals[index + 1].4 = format!("masked_overlap_after_P{}", index + 1);
                }
                SnvIndelOverlapPolicy::MarkOverlap => {
                    intervals[index].5.push(overlap);
                    intervals[index + 1].5.push(overlap);
                    intervals[index].4 = format!("marked_overlap_before_P{}", index + 2);
                    intervals[index + 1].4 = format!("marked_overlap_after_P{}", index + 1);
                }
            }
        }
    }

    let mut raw_segments = Vec::new();
    for (source_index, start, end, alignment, trim_note, overlap_query_intervals) in intervals {
        if start == 0 || end < start || start > read_sequence.len() {
            continue;
        }
        let end = end.min(read_sequence.len());
        let sequence = read_sequence[start - 1..end].to_string();
        if sequence.is_empty() {
            continue;
        }
        raw_segments.push((
            source_index,
            start,
            end,
            alignment.ss,
            alignment.se,
            alignment.strand,
            sequence,
            trim_note,
            overlap_query_intervals,
        ));
    }

    let segment_count = raw_segments.len();
    Ok(raw_segments
        .into_iter()
        .enumerate()
        .map(
            |(
                emitted_index,
                (
                    source_index,
                    start,
                    end,
                    subject_start,
                    subject_end,
                    strand,
                    sequence,
                    trim_note,
                    overlap_query_intervals,
                ),
            )| {
                let segment_index = emitted_index + 1;
                let segment_id = if segment_count == 1 {
                    record.read_id.clone()
                } else {
                    format!("{}_P{}", record.read_id, source_index + 1)
                };
                VariantSegment {
                    read_id: record.read_id.clone(),
                    segment_id,
                    read_class: metadata.read_class.clone(),
                    group_name: metadata.group_name.clone(),
                    subgroup_old_index: metadata.subgroup_old_index.clone(),
                    subgroup_key: metadata.subgroup_key.clone(),
                    segment_index,
                    segment_count,
                    query_start: start,
                    query_end: end,
                    subject_start,
                    subject_end,
                    strand,
                    sequence,
                    trim_note,
                    overlap_query_intervals,
                }
            },
        )
        .collect())
}

fn write_variant_segments_tsv(
    path: &Path,
    segments: &[VariantSegment],
    results: &[SegmentCallResult],
) -> Result<TsvWriteTiming, OrgraftError> {
    let build_started = Instant::now();
    let result_by_segment = result_by_segment_id(results);
    let mut buffer = String::with_capacity(segments.len().saturating_mul(512));
    buffer.push_str(
        "read_id\tsegment_id\tread_class\tgroup_name\tsubgroup_old_index\tsubgroup_key\tsegment_index\tsegment_count\tquery_start\tquery_end\tsegment_length\tsubject_start\tsubject_end\tstrand\ttrim_note\toverlap_query_intervals\tcall_count\tcall_elapsed_seconds\tcall_exit_status\tcall_message\n",
    );
    for segment in segments {
        let overlap_query_intervals = format_query_intervals(&segment.overlap_query_intervals);
        let result = result_by_segment.get(&segment.segment_id);
        let call_count = result
            .map(|result| result.calls.len().to_string())
            .unwrap_or_else(|| ".".to_string());
        let call_elapsed_seconds = result
            .map(|result| format!("{:.6}", result.elapsed_seconds))
            .unwrap_or_else(|| ".".to_string());
        let call_exit_status = result
            .map(|result| result.exit_status.to_string())
            .unwrap_or_else(|| ".".to_string());
        let call_message = result
            .map(|result| sanitize_tsv(&result.message))
            .unwrap_or_else(|| ".".to_string());
        writeln!(
            buffer,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            segment.read_id,
            segment.segment_id,
            segment.read_class,
            segment.group_name,
            segment.subgroup_old_index,
            segment.subgroup_key,
            segment.segment_index,
            segment.segment_count,
            segment.query_start,
            segment.query_end,
            segment.sequence.len(),
            segment.subject_start,
            segment.subject_end,
            segment.strand,
            segment.trim_note,
            overlap_query_intervals,
            call_count,
            call_elapsed_seconds,
            call_exit_status,
            call_message,
        )
        .expect("writing to String cannot fail");
    }
    write_built_tsv(path, buffer, build_started)
}

fn format_query_intervals(intervals: &[QueryInterval]) -> String {
    if intervals.is_empty() {
        ".".to_string()
    } else {
        intervals
            .iter()
            .map(|interval| format!("{}-{}", interval.start, interval.end))
            .collect::<Vec<_>>()
            .join(",")
    }
}

fn write_variant_segment_fasta(
    path: &Path,
    segments: &[VariantSegment],
) -> Result<(), OrgraftError> {
    let mut file = File::create(path)?;
    for segment in segments {
        writeln!(
            file,
            ">{} read_id={} group_name={} subgroup_old_index={} query={}-{}",
            segment.segment_id,
            segment.read_id,
            segment.group_name,
            segment.subgroup_old_index,
            segment.query_start,
            segment.query_end,
        )?;
        write_wrapped_sequence(&mut file, &segment.sequence)?;
    }
    Ok(())
}

fn run_custom_variant_caller_segments(
    minimap2: &Path,
    reference: &Path,
    segment_fasta: &Path,
    segments: &[VariantSegment],
    threads: usize,
    sv_context_filter: bool,
    stderr_path: &Path,
    commands: &mut Vec<CommandRecord>,
) -> Result<(Vec<SegmentCallResult>, usize, VariantStreamMetrics), OrgraftError> {
    writeln!(
        OpenOptions::new()
            .create(true)
            .append(true)
            .open(stderr_path)?,
        "### snv-indel-eval:round_1 custom CIGAR-diff caller stderr ###"
    )?;

    let reference_by_id = Arc::new(read_fasta_records_by_id(reference)?);
    let (_sam_header, sam_by_segment, stream_metrics) = run_minimap2_segment_sam_stream(
        minimap2,
        reference,
        segment_fasta,
        threads,
        stderr_path,
        commands,
    )?;
    let worker_count = threads.max(1).min(segments.len().max(1));
    let work_queue = Arc::new(Mutex::new(segments.to_vec().into_iter()));
    let sam_by_segment = Arc::new(sam_by_segment);
    let (sender, receiver) = mpsc::channel();
    let started = Instant::now();
    thread::scope(|scope| {
        for _ in 0..worker_count {
            let sender = sender.clone();
            let work_queue = Arc::clone(&work_queue);
            let sam_by_segment = Arc::clone(&sam_by_segment);
            let reference_by_id = Arc::clone(&reference_by_id);
            scope.spawn(move || {
                let mut chunk_results = Vec::new();
                loop {
                    let Some(segment) = work_queue
                        .lock()
                        .expect("segment queue lock poisoned")
                        .next()
                    else {
                        break;
                    };
                    chunk_results.push(run_custom_single_segment(
                        &reference_by_id,
                        &segment,
                        sam_by_segment.get(&segment.segment_id).map(Vec::as_slice),
                        sv_context_filter,
                    ));
                }
                let _ = sender.send(chunk_results);
            });
        }
    });
    drop(sender);

    let mut results = Vec::new();
    for worker_result in receiver {
        results.extend(worker_result);
    }
    let elapsed_seconds = started.elapsed().as_secs_f64();
    let failed_segments = results
        .iter()
        .filter(|result| result.exit_status != 0)
        .count();
    let total_calls = results
        .iter()
        .map(|result| result.calls.len())
        .sum::<usize>();
    writeln!(
        OpenOptions::new().append(true).open(stderr_path)?,
        "### snv-indel-eval:round_1 custom status=ok elapsed_seconds={elapsed_seconds:.3} workers={worker_count} segments={} total_calls={total_calls} failed_segments={failed_segments} sam_records={} ###\n",
        segments.len(),
        stream_metrics.split_sam_records,
    )?;
    commands.push(CommandRecord {
        timestamp: timestamp(),
        stage: "snv-indel-custom-cigar-diff",
        round: "round_1".to_string(),
        status: if failed_segments == 0 {
            "ok"
        } else {
            "partial"
        },
        elapsed_seconds,
        stdout: "memory:custom-cigar-diff-variants".to_string(),
        stderr: display_path(stderr_path),
        command: custom_variant_caller_command_template(reference, segments.len(), worker_count),
    });
    results.sort_by(|left, right| left.segment.segment_id.cmp(&right.segment.segment_id));
    Ok((results, worker_count, stream_metrics))
}

fn run_custom_single_segment(
    reference_by_id: &HashMap<String, String>,
    segment: &VariantSegment,
    sam_records: Option<&[SplitSamRecord]>,
    sv_context_filter: bool,
) -> SegmentCallResult {
    let started = Instant::now();
    let result = custom_calls_for_segment(reference_by_id, segment, sam_records, sv_context_filter);
    let elapsed_seconds = started.elapsed().as_secs_f64();
    match result {
        Ok(calls) => {
            let message = if sam_records.is_none() {
                "no_alignment"
            } else if calls.is_empty() {
                "no_calls"
            } else {
                "ok"
            }
            .to_string();
            SegmentCallResult {
                segment: segment.clone(),
                calls,
                elapsed_seconds,
                exit_status: 0,
                message,
                stderr: String::new(),
            }
        }
        Err(error) => SegmentCallResult {
            segment: segment.clone(),
            calls: Vec::new(),
            elapsed_seconds,
            exit_status: 1,
            message: "custom_caller_failed".to_string(),
            stderr: error.to_string(),
        },
    }
}

fn custom_calls_for_segment(
    reference_by_id: &HashMap<String, String>,
    segment: &VariantSegment,
    sam_records: Option<&[SplitSamRecord]>,
    sv_context_filter: bool,
) -> Result<Vec<VariantCall>, OrgraftError> {
    let Some(sam_records) = sam_records else {
        return Ok(Vec::new());
    };
    let mut calls = BTreeMap::<VariantCall, ()>::new();
    let selected_records = select_sam_records_for_segment(segment, sam_records, sv_context_filter)?;
    for record in selected_records {
        if record.reference_name == "*" || record.position == 0 || record.cigar == "*" {
            continue;
        }
        let reference = reference_by_id.get(&record.reference_name).ok_or_else(|| {
            OrgraftError::InvalidArgument(format!(
                "SAM reference `{}` is missing from polished reference FASTA",
                record.reference_name
            ))
        })?;
        for call in custom_calls_for_sam_record(record, reference)? {
            calls.insert(call, ());
        }
    }
    let mut calls = calls.keys().cloned().collect::<Vec<_>>();
    annotate_custom_call_confidence(&mut calls, segment);
    Ok(calls)
}

fn custom_calls_for_sam_record(
    record: &SplitSamRecord,
    reference: &str,
) -> Result<Vec<VariantCall>, OrgraftError> {
    let reference = reference.as_bytes();
    let query = record.sequence.as_bytes();
    let mut ref_index = record.position.saturating_sub(1);
    let mut query_index = 0usize;
    let mut calls = Vec::new();

    for (len, op) in parse_cigar(&record.cigar)? {
        match op {
            'M' | '=' | 'X' => {
                ensure_sam_bounds(reference, ref_index, len, "reference", record)?;
                ensure_sam_bounds(query, query_index, len, "query", record)?;
                for offset in 0..len {
                    let ref_base = reference[ref_index + offset].to_ascii_uppercase();
                    let query_base = query[query_index + offset].to_ascii_uppercase();
                    if ref_base != query_base
                        && base_index(ref_base).is_some()
                        && base_index(query_base).is_some()
                    {
                        calls.push(
                            VariantCall::high(
                                ref_index + offset + 1,
                                byte_to_base_string(ref_base),
                                byte_to_base_string(query_base),
                            )
                            .with_query_range(query_index + offset + 1, query_index + offset + 1),
                        );
                    }
                }
                ref_index += len;
                query_index += len;
            }
            'I' => {
                ensure_sam_bounds(query, query_index, len, "query", record)?;
                if let Some(call) = insertion_call(
                    reference,
                    ref_index,
                    &query[query_index..query_index + len],
                    query_index + 1,
                    query_index + len,
                ) {
                    calls.push(format_left_padded_indel(call, reference));
                }
                query_index += len;
            }
            'D' => {
                ensure_sam_bounds(reference, ref_index, len, "reference", record)?;
                if let Some(call) = deletion_call(
                    reference,
                    ref_index,
                    &reference[ref_index..ref_index + len],
                    query_index.max(1),
                ) {
                    calls.push(format_left_padded_indel(call, reference));
                }
                ref_index += len;
            }
            'N' => {
                ensure_sam_bounds(reference, ref_index, len, "reference", record)?;
                ref_index += len;
            }
            'S' => {
                ensure_sam_bounds(query, query_index, len, "query", record)?;
                query_index += len;
            }
            'H' | 'P' => {}
            _ => {
                return Err(OrgraftError::InvalidArgument(format!(
                    "unsupported CIGAR operation `{op}` in {}:{} {}",
                    record.reference_name, record.position, record.cigar
                )));
            }
        }
    }

    calls.sort();
    calls.dedup();
    Ok(calls)
}

fn ensure_sam_bounds(
    sequence: &[u8],
    start: usize,
    len: usize,
    label: &str,
    record: &SplitSamRecord,
) -> Result<(), OrgraftError> {
    if start
        .checked_add(len)
        .is_some_and(|end| end <= sequence.len())
    {
        Ok(())
    } else {
        Err(OrgraftError::InvalidArgument(format!(
            "SAM CIGAR {} runs past {label} length {} at {}:{} {}",
            len,
            sequence.len(),
            record.reference_name,
            record.position,
            record.cigar
        )))
    }
}

fn insertion_call(
    reference: &[u8],
    ref_index: usize,
    inserted: &[u8],
    query_start: usize,
    query_end: usize,
) -> Option<VariantCall> {
    if inserted.is_empty() || !inserted.iter().all(|base| base_index(*base).is_some()) {
        return None;
    }
    if ref_index > 0 {
        let anchor_index = ref_index - 1;
        let anchor = reference.get(anchor_index)?.to_ascii_uppercase();
        if base_index(anchor).is_none() {
            return None;
        }
        let mut alt = Vec::with_capacity(inserted.len() + 1);
        alt.push(anchor);
        alt.extend(inserted.iter().map(|base| base.to_ascii_uppercase()));
        Some(
            VariantCall::high(
                anchor_index + 1,
                byte_to_base_string(anchor),
                bases_to_string(&alt),
            )
            .with_query_range(query_start, query_end),
        )
    } else {
        let anchor = reference.first()?.to_ascii_uppercase();
        if base_index(anchor).is_none() {
            return None;
        }
        let mut alt = inserted
            .iter()
            .map(|base| base.to_ascii_uppercase())
            .collect::<Vec<_>>();
        alt.push(anchor);
        Some(
            VariantCall::high(1, byte_to_base_string(anchor), bases_to_string(&alt))
                .with_query_range(query_start, query_end),
        )
    }
}

fn deletion_call(
    reference: &[u8],
    ref_index: usize,
    deleted: &[u8],
    query_anchor: usize,
) -> Option<VariantCall> {
    if deleted.is_empty() || !deleted.iter().all(|base| base_index(*base).is_some()) {
        return None;
    }
    if ref_index > 0 {
        let anchor_index = ref_index - 1;
        let anchor = reference.get(anchor_index)?.to_ascii_uppercase();
        if base_index(anchor).is_none() {
            return None;
        }
        let mut ref_allele = Vec::with_capacity(deleted.len() + 1);
        ref_allele.push(anchor);
        ref_allele.extend(deleted.iter().map(|base| base.to_ascii_uppercase()));
        Some(
            VariantCall::high(
                anchor_index + 1,
                bases_to_string(&ref_allele),
                byte_to_base_string(anchor),
            )
            .with_query_range(query_anchor, query_anchor),
        )
    } else {
        let next_index = deleted.len();
        let anchor = reference.get(next_index)?.to_ascii_uppercase();
        if base_index(anchor).is_none() {
            return None;
        }
        let mut ref_allele = deleted
            .iter()
            .map(|base| base.to_ascii_uppercase())
            .collect::<Vec<_>>();
        ref_allele.push(anchor);
        Some(
            VariantCall::high(1, bases_to_string(&ref_allele), byte_to_base_string(anchor))
                .with_query_range(query_anchor, query_anchor),
        )
    }
}

fn format_left_padded_indel(call: VariantCall, reference: &[u8]) -> VariantCall {
    if call.ref_allele.len() == 1
        && call.alt_allele.len() > 1
        && call.alt_allele.starts_with(&call.ref_allele)
    {
        return extend_insertion_repeat_context(call, reference);
    }
    if call.alt_allele.len() == 1
        && call.ref_allele.len() > 1
        && call.ref_allele.starts_with(&call.alt_allele)
    {
        return extend_deletion_repeat_context(call, reference);
    }
    call
}

fn extend_insertion_repeat_context(call: VariantCall, reference: &[u8]) -> VariantCall {
    let inserted = call.alt_allele.as_bytes()[1..]
        .iter()
        .map(|base| base.to_ascii_uppercase())
        .collect::<Vec<_>>();
    if inserted.is_empty() || !inserted.iter().all(|base| base_index(*base).is_some()) {
        return call;
    }
    let context = repeated_right_context(reference, call.pos, &inserted);
    let anchor = reference[call.pos - 1].to_ascii_uppercase();
    let mut ref_allele = Vec::with_capacity(context.len() + 1);
    ref_allele.push(anchor);
    ref_allele.extend(context.iter().copied());
    let mut alt = Vec::with_capacity(inserted.len() + context.len() + 1);
    alt.push(anchor);
    alt.extend(inserted);
    alt.extend(context);
    VariantCall::high(
        call.pos,
        bases_to_string(&ref_allele),
        bases_to_string(&alt),
    )
    .with_metadata_from(&call)
}

fn extend_deletion_repeat_context(call: VariantCall, reference: &[u8]) -> VariantCall {
    let deleted = call.ref_allele.as_bytes()[1..]
        .iter()
        .map(|base| base.to_ascii_uppercase())
        .collect::<Vec<_>>();
    if deleted.is_empty() || !deleted.iter().all(|base| base_index(*base).is_some()) {
        return call;
    }
    let context = repeated_right_context(reference, call.pos + deleted.len(), &deleted);
    let anchor = reference[call.pos - 1].to_ascii_uppercase();
    let mut ref_allele = Vec::with_capacity(deleted.len() + context.len() + 1);
    ref_allele.push(anchor);
    ref_allele.extend(deleted);
    ref_allele.extend(context.iter().copied());
    let mut alt = Vec::with_capacity(context.len() + 1);
    alt.push(anchor);
    alt.extend(context);
    VariantCall::high(
        call.pos,
        bases_to_string(&ref_allele),
        bases_to_string(&alt),
    )
    .with_metadata_from(&call)
}

fn repeated_right_context(reference: &[u8], start: usize, motif: &[u8]) -> Vec<u8> {
    let mut context = Vec::new();
    if motif.is_empty() {
        return context;
    }
    let mut cursor = start;
    while cursor < reference.len() {
        let expected = motif[context.len() % motif.len()].to_ascii_uppercase();
        let observed = reference[cursor].to_ascii_uppercase();
        if observed != expected {
            break;
        }
        context.push(observed);
        cursor += 1;
    }
    context
}

fn annotate_custom_call_confidence(calls: &mut [VariantCall], segment: &VariantSegment) {
    let indel_indices = calls
        .iter()
        .enumerate()
        .filter_map(|(index, call)| (variant_call_type(call) == "InDel").then_some(index))
        .collect::<Vec<_>>();
    let segment_start = segment.subject_start.min(segment.subject_end).max(1) as usize;
    let segment_end = segment.subject_start.max(segment.subject_end).max(1) as usize;
    let short_terminal_segment = segment.segment_count > 1
        && (segment.segment_index == 1 || segment.segment_index == segment.segment_count)
        && segment.sequence.len() < SNV_INDEL_SHORT_TERMINAL_SEGMENT_BP;
    for index in 0..calls.len() {
        let call = &calls[index];
        let in_overlap = call_query_interval(call, segment)
            .map(|query_interval| {
                segment
                    .overlap_query_intervals
                    .iter()
                    .any(|overlap| query_intervals_overlap(query_interval, *overlap))
            })
            .unwrap_or(false);
        if in_overlap {
            calls[index].confidence = VariantConfidence::OverlapContext;
            continue;
        }
        if short_terminal_segment {
            calls[index].confidence = VariantConfidence::ShortSegmentContext;
            continue;
        }
        if variant_call_type(call) != "InDel" {
            continue;
        }
        let near_boundary = call.pos.abs_diff(segment_start) <= 30
            || call.pos.abs_diff(segment_end) <= 30
            || call.pos + call.ref_allele.len().max(call.alt_allele.len()) >= segment_end;
        let span = call.ref_allele.len().max(call.alt_allele.len());
        let delta = call.ref_allele.len().abs_diff(call.alt_allele.len());
        let dense_or_complex = span >= 20
            || delta >= 10
            || indel_indices.iter().any(|&other_index| {
                if other_index == index {
                    return false;
                }
                let other = &calls[other_index];
                let other_span = other.ref_allele.len().max(other.alt_allele.len());
                let distance = call.pos.abs_diff(other.pos);
                distance <= 20 || (distance <= 100 && (span >= 10 || other_span >= 10))
            });
        let confidence = if near_boundary {
            VariantConfidence::TerminalContext
        } else if dense_or_complex {
            VariantConfidence::ComplexContext
        } else {
            VariantConfidence::High
        };
        calls[index].confidence = confidence;
    }
}

fn call_query_interval(call: &VariantCall, segment: &VariantSegment) -> Option<QueryInterval> {
    let relative_start = call.query_start?;
    let relative_end = call.query_end?;
    Some(QueryInterval {
        start: segment.query_start + relative_start.saturating_sub(1),
        end: segment.query_start + relative_end.saturating_sub(1),
    })
}

fn query_intervals_overlap(left: QueryInterval, right: QueryInterval) -> bool {
    left.start <= right.end && right.start <= left.end
}

fn byte_to_base_string(base: u8) -> String {
    String::from_utf8(vec![base.to_ascii_uppercase()]).unwrap_or_else(|_| "N".to_string())
}

fn bases_to_string(bases: &[u8]) -> String {
    String::from_utf8(
        bases
            .iter()
            .map(|base| base.to_ascii_uppercase())
            .collect::<Vec<_>>(),
    )
    .unwrap_or_else(|_| "N".to_string())
}

fn run_minimap2_segment_sam_stream(
    minimap2: &Path,
    reference: &Path,
    segment_fasta: &Path,
    threads: usize,
    stderr_path: &Path,
    commands: &mut Vec<CommandRecord>,
) -> Result<
    (
        Vec<String>,
        HashMap<String, Vec<SplitSamRecord>>,
        VariantStreamMetrics,
    ),
    OrgraftError,
> {
    let mut stderr_file = OpenOptions::new()
        .create(true)
        .append(true)
        .open(stderr_path)?;
    writeln!(
        stderr_file,
        "### snv-indel-minimap2-stream:round_1 minimap2 stderr ###"
    )?;
    let stderr_for_child = stderr_file.try_clone()?;
    let mut command = Command::new(minimap2);
    command
        .arg("-t")
        .arg(threads.max(1).to_string())
        .arg("-ax")
        .arg("map-hifi")
        .arg(reference)
        .arg(segment_fasta);
    let command_text = format!("{command:?}");
    let started = Instant::now();
    let mut child = command
        .stdout(Stdio::piped())
        .stderr(Stdio::from(stderr_for_child))
        .spawn()?;
    let stdout = child.stdout.take().ok_or_else(|| {
        OrgraftError::InvalidArgument("failed to capture SNV minimap2 stdout".to_string())
    })?;
    let (header, by_segment, split_sam_records) =
        parse_minimap2_sam_stream(BufReader::new(stdout))?;
    let status = child.wait()?;
    let alignment_seconds = started.elapsed().as_secs_f64();
    let status_text = if status.success() { "ok" } else { "failed" };
    writeln!(
        OpenOptions::new().append(true).open(stderr_path)?,
        "### snv-indel-minimap2-stream:round_1 status={status_text} elapsed_seconds={alignment_seconds:.3} sam_records={split_sam_records} ###\n"
    )?;
    if !status.success() {
        return Err(OrgraftError::InvalidArgument(format!(
            "snv-indel minimap2 streaming round_1 failed; see {}",
            stderr_path.display()
        )));
    }
    commands.push(CommandRecord {
        timestamp: timestamp(),
        stage: "snv-indel-minimap2-stream",
        round: "round_1".to_string(),
        status: status_text,
        elapsed_seconds: alignment_seconds,
        stdout: "stream:minimap2-sam".to_string(),
        stderr: display_path(stderr_path),
        command: command_text,
    });
    Ok((
        header,
        by_segment,
        VariantStreamMetrics {
            alignment_seconds,
            split_sam_records,
        },
    ))
}

fn parse_minimap2_sam_stream<R: BufRead>(
    reader: R,
) -> Result<(Vec<String>, HashMap<String, Vec<SplitSamRecord>>, usize), OrgraftError> {
    let mut header = Vec::new();
    let mut by_segment: HashMap<String, Vec<SplitSamRecord>> = HashMap::new();
    let mut kept_records = 0usize;
    for line in reader.lines() {
        let line = line?;
        if line.starts_with('@') {
            header.push(line);
            continue;
        }
        if line.trim().is_empty() {
            continue;
        }
        let fields = line.split('\t').collect::<Vec<_>>();
        if fields.len() < 11 {
            continue;
        }
        let flag = parse_u16(fields[1], "SAM flag")?;
        if flag & 0x100 != 0 {
            continue;
        }
        let query_name = fields[0].to_string();
        let record = SplitSamRecord {
            flag,
            mapq: parse_u8(fields[4], "SAM MAPQ")?,
            reference_name: fields[2].to_string(),
            position: parse_usize_value(fields[3], "SAM POS")?,
            cigar: fields[5].to_string(),
            sequence: fields[9].to_ascii_uppercase(),
        };
        by_segment.entry(query_name).or_default().push(record);
        kept_records += 1;
    }
    for records in by_segment.values_mut() {
        records.sort_by(|left, right| {
            left.reference_name
                .cmp(&right.reference_name)
                .then_with(|| left.position.cmp(&right.position))
                .then_with(|| left.flag.cmp(&right.flag))
        });
    }
    Ok((header, by_segment, kept_records))
}

fn select_sam_records_for_segment<'a>(
    segment: &VariantSegment,
    sam_records: &'a [SplitSamRecord],
    sv_context_filter: bool,
) -> Result<Vec<&'a SplitSamRecord>, OrgraftError> {
    if !sv_context_filter {
        return Ok(sam_records.iter().collect());
    }
    let mut scored_records = Vec::new();
    let expected = QueryInterval {
        start: segment.subject_start.min(segment.subject_end).max(1) as usize,
        end: segment.subject_start.max(segment.subject_end).max(1) as usize,
    };
    for (index, record) in sam_records.iter().enumerate() {
        let Some(span) = sam_record_reference_span(record)? else {
            continue;
        };
        let distance = interval_distance(span, expected);
        let is_supplementary = usize::from(record.flag & 0x800 != 0);
        let edge_distance = span
            .start
            .abs_diff(expected.start)
            .min(span.end.abs_diff(expected.end));
        scored_records.push((
            distance,
            is_supplementary,
            Reverse(record.mapq),
            edge_distance,
            index,
            record,
        ));
    }

    if scored_records.is_empty() {
        return Ok(Vec::new());
    }
    if scored_records.len() == 1 {
        return Ok(vec![scored_records[0].5]);
    }

    scored_records.sort_by(|left, right| {
        left.0
            .cmp(&right.0)
            .then_with(|| left.1.cmp(&right.1))
            .then_with(|| left.2.cmp(&right.2))
            .then_with(|| left.3.cmp(&right.3))
            .then_with(|| left.4.cmp(&right.4))
    });
    Ok(vec![scored_records[0].5])
}

fn sam_record_reference_span(
    record: &SplitSamRecord,
) -> Result<Option<QueryInterval>, OrgraftError> {
    if record.reference_name == "*" || record.position == 0 || record.cigar == "*" {
        return Ok(None);
    }
    let mut reference_len = 0usize;
    for (len, op) in parse_cigar(&record.cigar)? {
        if matches!(op, 'M' | '=' | 'X' | 'D' | 'N') {
            reference_len += len;
        }
    }
    if reference_len == 0 {
        return Ok(None);
    }
    Ok(Some(QueryInterval {
        start: record.position,
        end: record.position + reference_len - 1,
    }))
}

fn interval_distance(left: QueryInterval, right: QueryInterval) -> usize {
    if query_intervals_overlap(left, right) {
        0
    } else if left.end < right.start {
        right.start - left.end
    } else {
        left.start - right.end
    }
}

fn write_per_variant_calls(
    path: &Path,
    segments: &[VariantSegment],
    results: &[SegmentCallResult],
) -> Result<TsvWriteTiming, OrgraftError> {
    let build_started = Instant::now();
    let result_by_segment = result_by_segment_id(results);
    let total_calls = results
        .iter()
        .map(|result| result.calls.len())
        .sum::<usize>();
    let mut buffer = String::with_capacity(total_calls.saturating_mul(160));
    buffer.push_str(
        "read_id\tsegment_id\tread_class\tgroup_name\tsubgroup_old_index\tsegment_index\tsegment_count\tpos\tref\talt\ttype\tconfidence\tconfidence_reason\n",
    );
    for segment in segments {
        let Some(result) = result_by_segment.get(&segment.segment_id) else {
            continue;
        };
        let mut calls = result.calls.clone();
        calls.sort();
        for call in calls {
            writeln!(
                buffer,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                segment.read_id,
                segment.segment_id,
                segment.read_class,
                segment.group_name,
                segment.subgroup_old_index,
                segment.segment_index,
                segment.segment_count,
                call.pos,
                call.ref_allele,
                call.alt_allele,
                variant_call_type(&call),
                call.confidence.as_str(),
                call.confidence.reason(),
            )
            .expect("writing to String cannot fail");
        }
    }
    write_built_tsv(path, buffer, build_started)
}

#[derive(Debug, Clone)]
struct LocalIndelAnnotation {
    indel_group: String,
    summary_anno: String,
    method: String,
    microhomology_size: Option<usize>,
    best_shift: isize,
}

#[derive(Debug, Clone)]
struct VariantTypeAnnotationSourceRow {
    pos: usize,
    ref_allele: String,
    alt_allele: String,
    variant_type: String,
    row_type: String,
    sample_info: String,
    counts: usize,
    id_list: String,
    call_count: usize,
    high_count: usize,
    overlap_context_count: usize,
    short_segment_context_count: usize,
    complex_context_count: usize,
    terminal_context_count: usize,
    primary_annotation: VariantTypeAnnotation,
    local_annotation: LocalIndelAnnotation,
}

impl VariantTypeAnnotationSourceRow {
    fn confidence_counts(&self) -> String {
        format!(
            "high={};overlap_context={};short_segment_context={};complex_context={};terminal_context={}",
            self.high_count,
            self.overlap_context_count,
            self.short_segment_context_count,
            self.complex_context_count,
            self.terminal_context_count,
        )
    }

    fn local_changed(&self) -> bool {
        self.variant_type == "InDel"
            && self.primary_annotation.indel_group != self.local_annotation.indel_group
    }
}

#[derive(Debug, Clone)]
struct CombinedVariantTypeAnnotationRow {
    compat: AnnotatedVariantCompatRow,
    source_variant_count: usize,
    call_count: usize,
    high_count: usize,
    overlap_context_count: usize,
    short_segment_context_count: usize,
    complex_context_count: usize,
    terminal_context_count: usize,
    primary_summary: String,
    local_changed_count: usize,
    same_pos_types: String,
}

impl CombinedVariantTypeAnnotationRow {
    fn confidence_counts(&self) -> String {
        format!(
            "high={};overlap_context={};short_segment_context={};complex_context={};terminal_context={}",
            self.high_count,
            self.overlap_context_count,
            self.short_segment_context_count,
            self.complex_context_count,
            self.terminal_context_count,
        )
    }
}

fn write_variant_type_annotation_tables(
    raw_path: &Path,
    combined_path: &Path,
    combined_high_path: &Path,
    plot_points_path: &Path,
    summary_path: &Path,
    runtime_summary: &SnvIndelRuntimeSummary,
    reference_path: &Path,
    coverage_path: &Path,
    segments: &[VariantSegment],
    results: &[SegmentCallResult],
) -> Result<VariantTypeAnnotationTimings, OrgraftError> {
    let total_started = Instant::now();
    let mut timings = VariantTypeAnnotationTimings::default();

    let started = Instant::now();
    let reference = first_sequence(reference_path)?;
    timings.read_reference_seconds = started.elapsed().as_secs_f64();

    let started = Instant::now();
    let frequency_depth = read_variant_frequency_depth(coverage_path)?;
    timings.read_frequency_depth_seconds = started.elapsed().as_secs_f64();

    let started = Instant::now();
    let (source_rows, compatible_type_by_pos) =
        collect_variant_type_annotation_rows(&reference, segments, results)?;
    timings.collect_source_rows_seconds = started.elapsed().as_secs_f64();
    timings.source_rows = source_rows.len();

    let raw_timing = write_raw_variant_type_annotations(
        raw_path,
        &source_rows,
        &frequency_depth,
        &compatible_type_by_pos,
    )?;
    timings.write_raw_seconds = raw_timing.total_seconds;
    timings.write_raw_bytes = raw_timing.bytes;

    let started = Instant::now();
    let combined_rows = combined_variant_type_annotation_rows(
        &source_rows,
        &frequency_depth,
        &compatible_type_by_pos,
    );
    timings.combine_rows_seconds = started.elapsed().as_secs_f64();
    timings.combined_rows = combined_rows.len();

    let combined_timing =
        write_combined_variant_type_annotations(combined_path, combined_rows.iter())?;
    timings.write_combined_seconds = combined_timing.total_seconds;
    timings.write_combined_bytes = combined_timing.bytes;

    let combined_high_timing = write_combined_variant_type_annotations(
        combined_high_path,
        combined_rows
            .iter()
            .filter(|row| row.compat.frequency >= 0.5),
    )?;
    timings.write_combined_high_seconds = combined_high_timing.total_seconds;
    timings.write_combined_high_bytes = combined_high_timing.bytes;

    let plot_points_timing = write_snv_indel_plot_points(plot_points_path, combined_rows.iter())?;
    timings.write_plot_points_seconds = plot_points_timing.total_seconds;
    timings.write_plot_points_bytes = plot_points_timing.bytes;

    let started = Instant::now();
    append_snv_indel_summary(summary_path, runtime_summary)?;
    timings.append_summary_seconds = started.elapsed().as_secs_f64();
    timings.total_seconds = total_started.elapsed().as_secs_f64();
    Ok(timings)
}

fn collect_variant_type_annotation_rows(
    reference: &str,
    segments: &[VariantSegment],
    results: &[SegmentCallResult],
) -> Result<
    (
        Vec<VariantTypeAnnotationSourceRow>,
        BTreeMap<usize, BTreeSet<String>>,
    ),
    OrgraftError,
> {
    let aggregates = aggregate_variant_calls_by_read(segments, results);
    let mut compatible_type_by_pos = BTreeMap::<usize, BTreeSet<String>>::new();
    let mut source_rows = Vec::new();
    for aggregate in aggregates {
        let call = &aggregate.call;
        let variant_type = variant_call_type(call);
        let annotation = annotate_variant_type(
            &reference,
            call.pos,
            &call.ref_allele,
            &call.alt_allele,
            variant_type,
        );
        let row_type = compatible_variant_row_type(variant_type, &annotation);
        compatible_type_by_pos
            .entry(call.pos)
            .or_default()
            .insert(row_type.clone());
        let local_annotation = annotate_local_indel_context(
            &reference,
            call.pos,
            &call.ref_allele,
            &call.alt_allele,
            variant_type,
            &annotation,
        );
        let single = AnnotatedSingleVariant {
            pos: call.pos,
            ref_allele: call.ref_allele.clone(),
            alt_allele: call.alt_allele.clone(),
            variant_type: variant_type.to_string(),
            id_list: aggregate.read_ids.join(","),
            counts: aggregate.read_ids.len(),
            type_annotation: annotation.clone(),
        };
        source_rows.push(VariantTypeAnnotationSourceRow {
            pos: single.pos,
            ref_allele: single.ref_allele.clone(),
            alt_allele: single.alt_allele.clone(),
            variant_type: single.variant_type.clone(),
            row_type,
            sample_info: variant_sample_info(&single),
            counts: single.counts,
            id_list: single.id_list,
            call_count: aggregate.call_count,
            high_count: aggregate.high_count,
            overlap_context_count: aggregate.overlap_context_count,
            short_segment_context_count: aggregate.short_segment_context_count,
            complex_context_count: aggregate.complex_context_count,
            terminal_context_count: aggregate.terminal_context_count,
            primary_annotation: annotation,
            local_annotation,
        });
    }

    source_rows.sort_by(compare_variant_type_annotation_rows);
    Ok((source_rows, compatible_type_by_pos))
}

fn write_raw_variant_type_annotations(
    path: &Path,
    source_rows: &[VariantTypeAnnotationSourceRow],
    frequency_depth: &[usize],
    compatible_type_by_pos: &BTreeMap<usize, BTreeSet<String>>,
) -> Result<TsvWriteTiming, OrgraftError> {
    let build_started = Instant::now();
    let mut buffer = String::with_capacity(source_rows.len().saturating_mul(320));
    let same_pos_labels = same_pos_types_by_pos(compatible_type_by_pos);
    buffer.push_str(
        "pos\tref\talt\ttype\tannotation_type\tread_count\tcall_count\tdepth\tfrequency\tconfidence_counts\tannotation_context\tannotation_group\tannotation_summary\tlocal_group\tlocal_summary\tlocal_method\tlocal_mh\tlocal_shift\tlocal_changed\tsame_pos_types\tread_ids\n",
    );
    for row in source_rows {
        let depth = frequency_depth.get(row.pos).copied().unwrap_or(0);
        let frequency = if depth == 0 {
            0.0
        } else {
            row.counts as f64 / depth as f64
        };
        write!(buffer, "{}", row.pos).expect("writing to String cannot fail");
        buffer.push('\t');
        push_sanitized_tsv(&mut buffer, &row.ref_allele);
        buffer.push('\t');
        push_sanitized_tsv(&mut buffer, &row.alt_allele);
        buffer.push('\t');
        push_sanitized_tsv(&mut buffer, &row.variant_type);
        buffer.push('\t');
        push_sanitized_tsv(&mut buffer, &row.row_type);
        write!(
            buffer,
            "\t{}\t{}\t{}\t{}",
            row.counts, row.call_count, depth, frequency
        )
        .expect("writing to String cannot fail");
        buffer.push('\t');
        push_sanitized_tsv(&mut buffer, &row.confidence_counts());
        buffer.push('\t');
        push_sanitized_tsv(&mut buffer, &row.primary_annotation.trinucleotide_context);
        buffer.push('\t');
        push_sanitized_tsv(&mut buffer, primary_group_label(row));
        buffer.push('\t');
        push_sanitized_tsv(&mut buffer, &row.primary_annotation.summary_anno);
        buffer.push('\t');
        push_sanitized_tsv(&mut buffer, &row.local_annotation.indel_group);
        buffer.push('\t');
        push_sanitized_tsv(&mut buffer, &row.local_annotation.summary_anno);
        buffer.push('\t');
        push_sanitized_tsv(&mut buffer, &row.local_annotation.method);
        buffer.push('\t');
        if let Some(value) = row.local_annotation.microhomology_size {
            write!(buffer, "{value}").expect("writing to String cannot fail");
        } else {
            buffer.push('.');
        }
        write!(
            buffer,
            "\t{}\t{}\t",
            row.local_annotation.best_shift,
            if row.local_changed() { "Yes" } else { "No" }
        )
        .expect("writing to String cannot fail");
        push_sanitized_tsv(
            &mut buffer,
            same_pos_labels
                .get(&row.pos)
                .map(String::as_str)
                .unwrap_or("."),
        );
        buffer.push('\t');
        push_sanitized_tsv(&mut buffer, &row.id_list);
        buffer.push('\n');
    }
    write_built_tsv(path, buffer, build_started)
}

fn combined_variant_type_annotation_rows(
    source_rows: &[VariantTypeAnnotationSourceRow],
    frequency_depth: &[usize],
    compatible_type_by_pos: &BTreeMap<usize, BTreeSet<String>>,
) -> Vec<CombinedVariantTypeAnnotationRow> {
    let same_pos_labels = same_pos_types_by_pos(compatible_type_by_pos);
    let mut by_pos = BTreeMap::<usize, Vec<VariantTypeAnnotationSourceRow>>::new();
    for row in source_rows.iter().cloned() {
        by_pos.entry(row.pos).or_default().push(row);
    }

    let mut combined_rows = Vec::new();
    for (_pos, rows) in by_pos {
        let mut snv_rows = Vec::new();
        let mut homopolymer_rows = Vec::new();
        let mut homodimer_rows = Vec::new();
        let mut tandem_rows = Vec::new();
        let mut other_rows = Vec::new();
        for row in rows {
            match row.row_type.as_str() {
                "SNV" => snv_rows.push(row),
                "InDel,homopolymer" => homopolymer_rows.push(row),
                "InDel,homodimer" => homodimer_rows.push(row),
                "InDel,tandem" => tandem_rows.push(row),
                _ => other_rows.push(row),
            }
        }
        push_combined_variant_type_annotation_group(
            &mut combined_rows,
            &snv_rows,
            Some("SNV"),
            &frequency_depth,
            &same_pos_labels,
        );
        push_combined_variant_type_annotation_group(
            &mut combined_rows,
            &homopolymer_rows,
            Some("InDel,homopolymer"),
            &frequency_depth,
            &same_pos_labels,
        );
        push_combined_variant_type_annotation_group(
            &mut combined_rows,
            &homodimer_rows,
            Some("InDel,homodimer"),
            &frequency_depth,
            &same_pos_labels,
        );
        push_combined_variant_type_annotation_group(
            &mut combined_rows,
            &tandem_rows,
            Some("InDel,tandem"),
            &frequency_depth,
            &same_pos_labels,
        );
        for row in other_rows {
            push_combined_variant_type_annotation_group(
                &mut combined_rows,
                &[row],
                None,
                &frequency_depth,
                &same_pos_labels,
            );
        }
    }
    combined_rows
}

fn push_combined_variant_type_annotation_group(
    combined_rows: &mut Vec<CombinedVariantTypeAnnotationRow>,
    source_rows: &[VariantTypeAnnotationSourceRow],
    row_type_override: Option<&str>,
    frequency_depth: &[usize],
    same_pos_labels: &BTreeMap<usize, String>,
) {
    if source_rows.is_empty() {
        return;
    }
    let compat_sources = source_rows
        .iter()
        .map(|row| VariantCombineSourceRow {
            pos: row.pos,
            ref_allele: row.ref_allele.clone(),
            row_type: row.row_type.clone(),
            alt_allele: row.alt_allele.clone(),
            sample_info: row.sample_info.clone(),
            counts: row.counts,
            id_list: row.id_list.clone(),
        })
        .collect::<Vec<_>>();
    let compat = finalize_compat_row(compat_sources, row_type_override, frequency_depth);
    combined_rows.push(CombinedVariantTypeAnnotationRow {
        same_pos_types: same_pos_labels
            .get(&compat.pos)
            .cloned()
            .unwrap_or_else(|| ".".to_string()),
        compat,
        source_variant_count: source_rows.len(),
        call_count: source_rows.iter().map(|row| row.call_count).sum::<usize>(),
        high_count: source_rows.iter().map(|row| row.high_count).sum::<usize>(),
        overlap_context_count: source_rows
            .iter()
            .map(|row| row.overlap_context_count)
            .sum::<usize>(),
        short_segment_context_count: source_rows
            .iter()
            .map(|row| row.short_segment_context_count)
            .sum::<usize>(),
        complex_context_count: source_rows
            .iter()
            .map(|row| row.complex_context_count)
            .sum::<usize>(),
        terminal_context_count: source_rows
            .iter()
            .map(|row| row.terminal_context_count)
            .sum::<usize>(),
        primary_summary: join_annotation_field(source_rows, |row| {
            row.primary_annotation.summary_anno.clone()
        }),
        local_changed_count: source_rows.iter().filter(|row| row.local_changed()).count(),
    });
}

fn write_combined_variant_type_annotations<'a, I>(
    path: &Path,
    rows: I,
) -> Result<TsvWriteTiming, OrgraftError>
where
    I: IntoIterator<Item = &'a CombinedVariantTypeAnnotationRow>,
{
    let build_started = Instant::now();
    let rows = rows.into_iter();
    let estimated_rows = rows.size_hint().1.unwrap_or(0);
    let mut buffer = String::with_capacity(estimated_rows.saturating_mul(300));
    buffer.push_str(
        "pos\tref\talt\ttype\ttotal_count\tdepth\tfrequency\tcombined_info\tmulti_allelic\tcounts\tread_ids\tfixed_ref\tsource_variant_count\tcall_count\tconfidence_counts\tannotation_summary\tlocal_changed_count\tsame_pos_types\n",
    );
    for row in rows {
        write!(buffer, "{}", row.compat.pos).expect("writing to String cannot fail");
        buffer.push('\t');
        push_sanitized_tsv(&mut buffer, &row.compat.ref_allele);
        buffer.push('\t');
        push_sanitized_tsv(&mut buffer, &row.compat.alt_allele);
        buffer.push('\t');
        push_sanitized_tsv(&mut buffer, &row.compat.row_type);
        write!(
            buffer,
            "\t{}\t{}\t{}",
            row.compat.total_count, row.compat.depth, row.compat.frequency
        )
        .expect("writing to String cannot fail");
        buffer.push('\t');
        push_sanitized_tsv(&mut buffer, &row.compat.combined_info);
        write!(buffer, "\t{}\t", row.compat.multi_allelic).expect("writing to String cannot fail");
        push_sanitized_tsv(&mut buffer, &row.compat.counts);
        buffer.push('\t');
        push_sanitized_tsv(&mut buffer, &row.compat.id_list);
        write!(
            buffer,
            "\t{}\t{}\t{}\t",
            row.compat.fixed_ref, row.source_variant_count, row.call_count
        )
        .expect("writing to String cannot fail");
        push_sanitized_tsv(&mut buffer, &row.confidence_counts());
        buffer.push('\t');
        push_sanitized_tsv(&mut buffer, &row.primary_summary);
        write!(buffer, "\t{}\t", row.local_changed_count).expect("writing to String cannot fail");
        push_sanitized_tsv(&mut buffer, &row.same_pos_types);
        buffer.push('\n');
    }
    write_built_tsv(path, buffer, build_started)
}

fn write_snv_indel_plot_points<'a, I>(path: &Path, rows: I) -> Result<TsvWriteTiming, OrgraftError>
where
    I: IntoIterator<Item = &'a CombinedVariantTypeAnnotationRow>,
{
    let build_started = Instant::now();
    let rows = rows.into_iter();
    let estimated_rows = rows.size_hint().1.unwrap_or(0);
    let mut buffer = String::with_capacity(estimated_rows.saturating_mul(240));
    buffer.push_str(
        "pos\tvariant_kind\ttype\tplot_type\ttotal_count\tdepth\tfrequency\tref_count\tref_frequency\tmax_alt_count\tmax_alt_frequency\tref_is_top\thigh_count\tnon_high_count\tconfidence_class\tconfidence_counts\tref\talt\tfixed_ref\tmulti_allelic\tcounts\n",
    );
    for row in rows {
        let variant_kind = if row.compat.row_type == "SNV" {
            "SNV"
        } else {
            "InDel"
        };
        let ref_count = row.compat.depth.saturating_sub(row.compat.total_count);
        let ref_frequency = fraction_or_zero(ref_count, row.compat.depth);
        let max_alt_count = max_count_field(&row.compat.counts).unwrap_or(row.compat.total_count);
        let max_alt_frequency = fraction_or_zero(max_alt_count, row.compat.depth);
        let ref_is_top = ref_count >= max_alt_count;
        let non_high_count = row.call_count.saturating_sub(row.high_count);
        let confidence_class = if non_high_count == 0 {
            "high"
        } else {
            "non_high"
        };
        write!(buffer, "{}\t{}\t", row.compat.pos, variant_kind)
            .expect("writing to String cannot fail");
        push_sanitized_tsv(&mut buffer, &row.compat.row_type);
        buffer.push('\t');
        push_sanitized_tsv(&mut buffer, &snv_indel_plot_type(row));
        write!(
            buffer,
            "\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t",
            row.compat.total_count,
            row.compat.depth,
            row.compat.frequency,
            ref_count,
            ref_frequency,
            max_alt_count,
            max_alt_frequency,
            ref_is_top,
            row.high_count,
            non_high_count,
            confidence_class,
        )
        .expect("writing to String cannot fail");
        push_sanitized_tsv(&mut buffer, &row.confidence_counts());
        buffer.push('\t');
        push_sanitized_tsv(&mut buffer, &row.compat.ref_allele);
        buffer.push('\t');
        push_sanitized_tsv(&mut buffer, &row.compat.alt_allele);
        write!(
            buffer,
            "\t{}\t{}\t",
            row.compat.fixed_ref, row.compat.multi_allelic
        )
        .expect("writing to String cannot fail");
        push_sanitized_tsv(&mut buffer, &row.compat.counts);
        buffer.push('\n');
    }
    write_built_tsv(path, buffer, build_started)
}

fn fraction_or_zero(count: usize, depth: usize) -> f64 {
    if depth == 0 {
        0.0
    } else {
        count as f64 / depth as f64
    }
}

fn max_count_field(counts: &str) -> Option<usize> {
    counts
        .split('#')
        .filter_map(|value| value.parse::<usize>().ok())
        .max()
}

fn snv_indel_plot_type(row: &CombinedVariantTypeAnnotationRow) -> String {
    if row.compat.row_type == "SNV" {
        return snv_plot_type(row);
    }
    indel_plot_type(row)
}

fn snv_plot_type(row: &CombinedVariantTypeAnnotationRow) -> String {
    let mut labels = BTreeSet::new();
    if row.compat.ref_allele.len() == 1 {
        for alt in row.compat.alt_allele.split('#') {
            if alt.len() == 1 {
                labels.insert(format!("{}>{}", row.compat.ref_allele, alt));
            }
        }
    }
    if labels.is_empty() {
        "SNV,other".to_string()
    } else {
        labels.into_iter().collect::<Vec<_>>().join(";")
    }
}

fn indel_plot_type(row: &CombinedVariantTypeAnnotationRow) -> String {
    match row.compat.row_type.as_str() {
        "InDel,homopolymer" => row
            .compat
            .combined_info
            .split(';')
            .next()
            .filter(|label| label.starts_with("poly-"))
            .unwrap_or("poly-other")
            .to_string(),
        "InDel,homodimer" => "tandem_size_2".to_string(),
        "InDel,tandem" => tandem_plot_type(&row.compat.combined_info),
        "InDel,MMEJ" => "InDel,MMEJ".to_string(),
        "InDel,NHEJ" => "InDel,NHEJ".to_string(),
        other => other.to_string(),
    }
}

fn tandem_plot_type(combined_info: &str) -> String {
    let unit_len = combined_info
        .split(';')
        .next()
        .and_then(|prefix| prefix.strip_prefix("tandem-"))
        .map(str::len)
        .unwrap_or(0);
    match unit_len {
        2..=6 => format!("tandem_size_{unit_len}"),
        _ => "tandem_size_other".to_string(),
    }
}

fn append_snv_indel_summary(
    path: &Path,
    runtime: &SnvIndelRuntimeSummary,
) -> Result<(), OrgraftError> {
    let mut file = OpenOptions::new().create(true).append(true).open(path)?;
    for (metric, value) in snv_runtime_summary_rows(runtime) {
        writeln!(
            file,
            "snv_indel_{}\t{}",
            sanitize_tsv(&metric),
            sanitize_tsv(&value),
        )?;
    }
    Ok(())
}

fn append_snv_indel_write_timing_summary(
    path: &Path,
    timings: &SnvIndelWriteTimings,
) -> Result<(), OrgraftError> {
    let mut file = OpenOptions::new().create(true).append(true).open(path)?;
    for (metric, value) in snv_write_timing_rows(timings) {
        writeln!(
            file,
            "snv_indel_{}\t{}",
            sanitize_tsv(&metric),
            sanitize_tsv(&value),
        )?;
    }
    Ok(())
}

fn snv_runtime_summary_rows(runtime: &SnvIndelRuntimeSummary) -> Vec<(String, String)> {
    vec![
        ("call_mode".to_string(), runtime.call_mode.to_string()),
        (
            "overlap_policy".to_string(),
            runtime.overlap_policy.to_string(),
        ),
        (
            "sv_context_filter".to_string(),
            runtime.sv_context_filter.to_string(),
        ),
        (
            "minimap2_preset".to_string(),
            runtime.minimap2_preset.to_string(),
        ),
        ("reference_path".to_string(), runtime.reference_path.clone()),
        (
            "fl_read_count".to_string(),
            runtime.fl_read_count.to_string(),
        ),
        (
            "segment_count".to_string(),
            runtime.segment_count.to_string(),
        ),
        ("total_calls".to_string(), runtime.total_calls.to_string()),
        (
            "reads_with_calls".to_string(),
            runtime.reads_with_calls.to_string(),
        ),
        (
            "failed_segments".to_string(),
            runtime.failed_segments.to_string(),
        ),
        ("workers".to_string(), runtime.workers.to_string()),
        (
            "elapsed_seconds".to_string(),
            format!("{:.6}", runtime.elapsed_seconds),
        ),
        (
            "shared_minimap2_stream_seconds".to_string(),
            format!("{:.6}", runtime.shared_minimap2_stream_seconds),
        ),
        (
            "split_sam_records".to_string(),
            runtime.split_sam_records.to_string(),
        ),
        (
            "sum_segment_seconds".to_string(),
            format!("{:.6}", runtime.sum_segment_seconds),
        ),
    ]
}

fn snv_write_timing_rows(timings: &SnvIndelWriteTimings) -> Vec<(String, String)> {
    let variant_type = &timings.variant_type;
    vec![
        (
            "write_per_variant_calls_seconds".to_string(),
            format!("{:.6}", timings.per_variant_calls_seconds),
        ),
        (
            "write_per_variant_calls_bytes".to_string(),
            timings.per_variant_calls_bytes.to_string(),
        ),
        (
            "write_variant_segments_seconds".to_string(),
            format!("{:.6}", timings.variant_segments_seconds),
        ),
        (
            "write_variant_segments_bytes".to_string(),
            timings.variant_segments_bytes.to_string(),
        ),
        (
            "write_variant_type_total_seconds".to_string(),
            format!("{:.6}", variant_type.total_seconds),
        ),
        (
            "write_variant_type_read_reference_seconds".to_string(),
            format!("{:.6}", variant_type.read_reference_seconds),
        ),
        (
            "write_variant_type_read_frequency_depth_seconds".to_string(),
            format!("{:.6}", variant_type.read_frequency_depth_seconds),
        ),
        (
            "write_variant_type_collect_source_rows_seconds".to_string(),
            format!("{:.6}", variant_type.collect_source_rows_seconds),
        ),
        (
            "write_variant_type_raw_tsv_seconds".to_string(),
            format!("{:.6}", variant_type.write_raw_seconds),
        ),
        (
            "write_variant_type_raw_tsv_bytes".to_string(),
            variant_type.write_raw_bytes.to_string(),
        ),
        (
            "write_variant_type_combine_rows_seconds".to_string(),
            format!("{:.6}", variant_type.combine_rows_seconds),
        ),
        (
            "write_variant_type_combined_tsv_seconds".to_string(),
            format!("{:.6}", variant_type.write_combined_seconds),
        ),
        (
            "write_variant_type_combined_tsv_bytes".to_string(),
            variant_type.write_combined_bytes.to_string(),
        ),
        (
            "write_variant_type_high_tsv_seconds".to_string(),
            format!("{:.6}", variant_type.write_combined_high_seconds),
        ),
        (
            "write_variant_type_high_tsv_bytes".to_string(),
            variant_type.write_combined_high_bytes.to_string(),
        ),
        (
            "write_variant_type_plot_points_tsv_seconds".to_string(),
            format!("{:.6}", variant_type.write_plot_points_seconds),
        ),
        (
            "write_variant_type_plot_points_tsv_bytes".to_string(),
            variant_type.write_plot_points_bytes.to_string(),
        ),
        (
            "write_variant_type_append_summary_seconds".to_string(),
            format!("{:.6}", variant_type.append_summary_seconds),
        ),
        (
            "write_variant_type_source_rows".to_string(),
            variant_type.source_rows.to_string(),
        ),
        (
            "write_variant_type_combined_rows".to_string(),
            variant_type.combined_rows.to_string(),
        ),
        (
            "write_failed_segments_log_seconds".to_string(),
            format!("{:.6}", timings.failed_segments_log_seconds),
        ),
        (
            "write_total_seconds".to_string(),
            format!("{:.6}", timings.total_seconds),
        ),
    ]
}

fn compare_variant_type_annotation_rows(
    left: &VariantTypeAnnotationSourceRow,
    right: &VariantTypeAnnotationSourceRow,
) -> Ordering {
    left.pos
        .cmp(&right.pos)
        .then_with(|| {
            variant_row_type_priority(&left.row_type)
                .cmp(&variant_row_type_priority(&right.row_type))
        })
        .then_with(|| left.ref_allele.cmp(&right.ref_allele))
        .then_with(|| left.alt_allele.cmp(&right.alt_allele))
}

fn variant_row_type_priority(row_type: &str) -> usize {
    match row_type {
        "SNV" => 0,
        "InDel,homopolymer" => 1,
        "InDel,homodimer" => 2,
        "InDel,tandem" => 3,
        "InDel,MMEJ" => 4,
        "InDel,NHEJ" => 5,
        _ => 6,
    }
}

fn primary_group_label(row: &VariantTypeAnnotationSourceRow) -> &str {
    if row.variant_type == "SNV" {
        &row.primary_annotation.snv_group
    } else {
        &row.primary_annotation.indel_group
    }
}

fn same_pos_types_label(
    compatible_type_by_pos: &BTreeMap<usize, BTreeSet<String>>,
    pos: usize,
) -> String {
    let Some(types) = compatible_type_by_pos.get(&pos) else {
        return ".".to_string();
    };
    let mut types = types.iter().cloned().collect::<Vec<_>>();
    types.sort_by(|left, right| {
        variant_row_type_priority(left)
            .cmp(&variant_row_type_priority(right))
            .then_with(|| left.cmp(right))
    });
    types.join(",")
}

fn join_annotation_field<F>(rows: &[VariantTypeAnnotationSourceRow], mut value: F) -> String
where
    F: FnMut(&VariantTypeAnnotationSourceRow) -> String,
{
    rows.iter()
        .map(|row| value(row))
        .collect::<Vec<_>>()
        .join("#")
}

fn compatible_variant_row_type(variant_type: &str, annotation: &VariantTypeAnnotation) -> String {
    if variant_type == "SNV" {
        "SNV".to_string()
    } else if variant_type == "InDel" {
        format!("InDel,{}", annotation.indel_group)
    } else {
        variant_type.to_string()
    }
}

fn annotate_local_indel_context(
    reference: &str,
    pos: usize,
    ref_allele: &str,
    alt_allele: &str,
    variant_type: &str,
    compatible: &VariantTypeAnnotation,
) -> LocalIndelAnnotation {
    let indel_size = alt_allele.len() as isize - ref_allele.len() as isize;
    if variant_type != "InDel" {
        return LocalIndelAnnotation {
            indel_group: "-".to_string(),
            summary_anno: "-".to_string(),
            method: "not_indel".to_string(),
            microhomology_size: None,
            best_shift: 0,
        };
    }

    if matches!(
        compatible.indel_group.as_str(),
        "homopolymer" | "homodimer" | "tandem"
    ) {
        return LocalIndelAnnotation {
            indel_group: compatible.indel_group.clone(),
            summary_anno: compatible.summary_anno.clone(),
            method: "compatible_repeat_slippage".to_string(),
            microhomology_size: None,
            best_shift: 0,
        };
    }

    if let Some(local) =
        local_shift_microhomology(reference, pos, ref_allele, alt_allele, indel_size)
    {
        return local;
    }

    LocalIndelAnnotation {
        indel_group: compatible.indel_group.clone(),
        summary_anno: compatible.summary_anno.clone(),
        method: "compatible_prefix_microhomology".to_string(),
        microhomology_size: microhomology_from_summary(&compatible.summary_anno),
        best_shift: 0,
    }
}

fn local_shift_microhomology(
    reference: &str,
    pos: usize,
    ref_allele: &str,
    alt_allele: &str,
    indel_size: isize,
) -> Option<LocalIndelAnnotation> {
    let prefix_len = common_prefix_len(ref_allele, alt_allele);
    let suffix_len = common_suffix_len_after_prefix(ref_allele, alt_allele, prefix_len);
    let ref_core_end = ref_allele.len().saturating_sub(suffix_len);
    let alt_core_end = alt_allele.len().saturating_sub(suffix_len);
    let ref_core = &ref_allele[prefix_len..ref_core_end];
    let alt_core = &alt_allele[prefix_len..alt_core_end];
    let event_start = pos.checked_sub(1)? + prefix_len;
    let window = 100usize;

    if !ref_core.is_empty() && alt_core.is_empty() {
        let deleted_len = ref_core.len();
        let window_start = event_start.saturating_sub(window);
        let window_end = (event_start + deleted_len + window).min(reference.len());
        let ref_window = &reference[window_start..window_end];
        let local_start = event_start - window_start;
        if local_start + deleted_len > ref_window.len() {
            return None;
        }
        let alt_window = format!(
            "{}{}",
            &ref_window[..local_start],
            &ref_window[local_start + deleted_len..]
        );
        let mut best_mh = 0usize;
        let mut best_shift = 0isize;
        for candidate_start in 0..=ref_window.len().saturating_sub(deleted_len) {
            if deletion_matches_alt_window(ref_window, &alt_window, candidate_start, deleted_len) {
                let global_start = window_start + candidate_start;
                let mh = suffix_match_len(
                    &reference.as_bytes()[..global_start + deleted_len],
                    &reference.as_bytes()[..global_start],
                    99,
                );
                let shift = candidate_start as isize - local_start as isize;
                if mh > best_mh || (mh == best_mh && shift.abs() < best_shift.abs()) {
                    best_mh = mh;
                    best_shift = shift;
                }
            }
        }
        return Some(local_microhomology_annotation(
            best_mh,
            indel_size,
            best_shift,
            "local_shift_deletion",
        ));
    }

    if ref_core.is_empty() && !alt_core.is_empty() {
        let inserted = alt_core;
        let window_start = event_start.saturating_sub(window);
        let window_end = (event_start + window).min(reference.len());
        let ref_window = &reference[window_start..window_end];
        let local_start = event_start - window_start;
        let alt_window = format!(
            "{}{}{}",
            &ref_window[..local_start],
            inserted,
            &ref_window[local_start..]
        );
        let mut best_mh = 0usize;
        let mut best_shift = 0isize;
        for candidate_start in 0..=ref_window.len() {
            if insertion_matches_alt_window(
                ref_window,
                &alt_window,
                candidate_start,
                inserted.len(),
            ) {
                let global_start = window_start + candidate_start;
                let mut alt_prefix = reference[..global_start].to_string();
                alt_prefix.push_str(inserted);
                let mh = suffix_match_len(
                    reference[..global_start].as_bytes(),
                    alt_prefix.as_bytes(),
                    99,
                );
                let shift = candidate_start as isize - local_start as isize;
                if mh > best_mh || (mh == best_mh && shift.abs() < best_shift.abs()) {
                    best_mh = mh;
                    best_shift = shift;
                }
            }
        }
        return Some(local_microhomology_annotation(
            best_mh,
            indel_size,
            best_shift,
            "local_shift_insertion",
        ));
    }

    None
}

fn local_microhomology_annotation(
    microhomology_size: usize,
    indel_size: isize,
    best_shift: isize,
    method: &str,
) -> LocalIndelAnnotation {
    let indel_group = if microhomology_size <= 1 {
        "NHEJ"
    } else {
        "MMEJ"
    };
    LocalIndelAnnotation {
        indel_group: indel_group.to_string(),
        summary_anno: format!("MH_size={microhomology_size};indel_size={indel_size}"),
        method: method.to_string(),
        microhomology_size: Some(microhomology_size),
        best_shift,
    }
}

fn deletion_matches_alt_window(
    ref_window: &str,
    alt_window: &str,
    deletion_start: usize,
    deletion_len: usize,
) -> bool {
    ref_window[..deletion_start] == alt_window[..deletion_start]
        && ref_window[deletion_start + deletion_len..] == alt_window[deletion_start..]
}

fn insertion_matches_alt_window(
    ref_window: &str,
    alt_window: &str,
    insertion_start: usize,
    insertion_len: usize,
) -> bool {
    ref_window[..insertion_start] == alt_window[..insertion_start]
        && ref_window[insertion_start..] == alt_window[insertion_start + insertion_len..]
}

fn microhomology_from_summary(summary: &str) -> Option<usize> {
    summary
        .split(';')
        .find_map(|field| field.strip_prefix("MH_size="))
        .and_then(|value| value.parse::<usize>().ok())
}

#[derive(Debug, Clone)]
struct VariantAggregate {
    call: VariantCall,
    read_ids: Vec<String>,
    call_count: usize,
    high_count: usize,
    overlap_context_count: usize,
    short_segment_context_count: usize,
    complex_context_count: usize,
    terminal_context_count: usize,
}

#[derive(Debug, Clone)]
struct VariantTypeAnnotation {
    annotation: &'static str,
    trinucleotide_context: String,
    snv_group: String,
    indel_group: String,
    summary_anno: String,
}

#[derive(Debug, Clone)]
struct AnnotatedSingleVariant {
    pos: usize,
    ref_allele: String,
    alt_allele: String,
    variant_type: String,
    id_list: String,
    counts: usize,
    type_annotation: VariantTypeAnnotation,
}

#[derive(Debug, Clone)]
struct VariantCombineSourceRow {
    pos: usize,
    ref_allele: String,
    row_type: String,
    alt_allele: String,
    sample_info: String,
    counts: usize,
    id_list: String,
}

#[derive(Debug, Clone)]
struct AnnotatedVariantCompatRow {
    pos: usize,
    ref_allele: String,
    alt_allele: String,
    row_type: String,
    total_count: usize,
    combined_info: String,
    multi_allelic: String,
    counts: String,
    id_list: String,
    fixed_ref: String,
    depth: usize,
    frequency: f64,
}

fn aggregate_variant_calls_by_read(
    segments: &[VariantSegment],
    results: &[SegmentCallResult],
) -> Vec<VariantAggregate> {
    let mut read_order = Vec::new();
    let mut segment_ids_by_read: HashMap<&str, Vec<&str>> = HashMap::new();
    let mut seen_reads = HashSet::new();
    for segment in segments {
        if seen_reads.insert(segment.read_id.as_str()) {
            read_order.push(segment.read_id.as_str());
        }
        segment_ids_by_read
            .entry(segment.read_id.as_str())
            .or_default()
            .push(segment.segment_id.as_str());
    }

    let result_by_segment = results
        .iter()
        .map(|result| (result.segment.segment_id.as_str(), result))
        .collect::<HashMap<_, _>>();
    let mut by_call = BTreeMap::<VariantCall, VariantAggregate>::new();
    for result in results {
        for call in &result.calls {
            let aggregate = by_call
                .entry(call.clone())
                .or_insert_with(|| VariantAggregate {
                    call: call.clone(),
                    read_ids: Vec::new(),
                    call_count: 0,
                    high_count: 0,
                    overlap_context_count: 0,
                    short_segment_context_count: 0,
                    complex_context_count: 0,
                    terminal_context_count: 0,
                });
            aggregate.call_count += 1;
            match call.confidence {
                VariantConfidence::High => aggregate.high_count += 1,
                VariantConfidence::OverlapContext => aggregate.overlap_context_count += 1,
                VariantConfidence::ShortSegmentContext => {
                    aggregate.short_segment_context_count += 1
                }
                VariantConfidence::ComplexContext => aggregate.complex_context_count += 1,
                VariantConfidence::TerminalContext => aggregate.terminal_context_count += 1,
            }
        }
    }
    for read_id in read_order {
        let mut calls_for_read = BTreeMap::<VariantCall, ()>::new();
        for segment_id in segment_ids_by_read
            .get(read_id)
            .map(Vec::as_slice)
            .unwrap_or(&[])
        {
            let Some(result) = result_by_segment.get(segment_id) else {
                continue;
            };
            for call in &result.calls {
                calls_for_read.insert(call.clone(), ());
            }
        }
        for call in calls_for_read.keys() {
            if let Some(aggregate) = by_call.get_mut(call) {
                aggregate.read_ids.push(read_id.to_string());
            }
        }
    }

    by_call.into_values().collect()
}

fn annotate_variant_type(
    reference: &str,
    pos: usize,
    ref_allele: &str,
    alt_allele: &str,
    variant_type: &str,
) -> VariantTypeAnnotation {
    if variant_type == "SNV" {
        return VariantTypeAnnotation {
            annotation: "Yes",
            trinucleotide_context: trinucleotide_context(reference, pos, ref_allele),
            snv_group: format!("{ref_allele}>{alt_allele}"),
            indel_group: "-".to_string(),
            summary_anno: "-".to_string(),
        };
    }

    let mut annotation = VariantTypeAnnotation {
        annotation: "No",
        trinucleotide_context: "-".to_string(),
        snv_group: "-".to_string(),
        indel_group: "-".to_string(),
        summary_anno: "-".to_string(),
    };
    if variant_type != "InDel" {
        return annotation;
    }

    if ref_allele.len() > 1 && alt_allele.len() > 1 {
        let ref_str = &ref_allele[1..];
        let alt_str = &alt_allele[1..];
        let ref_unit_base = &ref_allele[1..2];
        let alt_unit_base = &alt_allele[1..2];
        if ref_unit_base == alt_unit_base
            && ref_str
                .chars()
                .all(|base| base.to_string() == ref_unit_base)
            && alt_str
                .chars()
                .all(|base| base.to_string() == alt_unit_base)
        {
            annotation.annotation = "Yes";
            annotation.indel_group = "homopolymer".to_string();
            annotation.summary_anno = format!(
                "poly-{};ref_size={};alt_size={}",
                ref_unit_base,
                ref_str.len(),
                alt_str.len()
            );
        }
    }

    if ref_allele.len() > 2 && alt_allele.len() > 2 {
        let ref_str = &ref_allele[1..];
        let alt_str = &alt_allele[1..];
        let ref_unit = &ref_allele[1..3];
        let alt_unit = &alt_allele[1..3];
        if ref_unit == alt_unit
            && ref_unit.as_bytes()[0] != ref_unit.as_bytes()[1]
            && is_exact_repeat(ref_str, ref_unit)
            && is_exact_repeat(alt_str, alt_unit)
        {
            annotation.annotation = "Yes";
            annotation.indel_group = "homodimer".to_string();
            annotation.summary_anno = format!(
                "dimer-{};ref_size={};alt_size={}",
                ref_unit,
                ref_str.len() / 2,
                alt_str.len() / 2
            );
        }
    }

    if annotation.annotation == "No" && ref_allele.len() > 1 && alt_allele.len() > 1 {
        let ref_str = &ref_allele[1..];
        let alt_str = &alt_allele[1..];
        let unit_len = gcd_usize(ref_str.len(), alt_str.len());
        if unit_len > 0 {
            let ref_count = ref_str.len() / unit_len;
            let alt_count = alt_str.len() / unit_len;
            let ref_unit = &ref_str[..unit_len];
            if ref_allele.as_bytes()[0] == alt_allele.as_bytes()[0]
                && ref_str == ref_unit.repeat(ref_count)
                && alt_str == ref_unit.repeat(alt_count)
            {
                annotation.annotation = "Yes";
                annotation.indel_group = "tandem".to_string();
                annotation.summary_anno = format!(
                    "tandem-{};ref_size={};alt_size={}",
                    ref_unit, ref_count, alt_count
                );
            }
        }
    }

    if annotation.annotation == "No" {
        let mut leading_base_count = 2usize;
        while ref_prefix(ref_allele, leading_base_count)
            == ref_prefix(alt_allele, leading_base_count)
        {
            if ref_allele.len() > leading_base_count && alt_allele.len() > leading_base_count {
                let ref_str = &ref_allele[leading_base_count..];
                let alt_str = &alt_allele[leading_base_count..];
                let unit_len = gcd_usize(ref_str.len(), alt_str.len());
                if unit_len > 0 {
                    let ref_count = ref_str.len() / unit_len;
                    let alt_count = alt_str.len() / unit_len;
                    let ref_unit = &ref_str[..unit_len];
                    if ref_str == ref_unit.repeat(ref_count)
                        && alt_str == ref_unit.repeat(alt_count)
                    {
                        annotation.annotation = "Yes";
                        annotation.indel_group = "tandem".to_string();
                        annotation.summary_anno = format!(
                            "tandem-{};ref_size={};alt_size={}",
                            ref_unit, ref_count, alt_count
                        );
                        break;
                    }
                }
            }
            leading_base_count += 1;
            if leading_base_count == ref_allele.len() || leading_base_count == alt_allele.len() {
                break;
            }
        }
    }

    if annotation.annotation == "No" {
        let leading_base_count = common_prefix_len(ref_allele, alt_allele);
        if leading_base_count == ref_allele.len() {
            let ref_border_pos = pos + ref_allele.len() - 1;
            let ref_border_seq = reference[..ref_border_pos.min(reference.len())].to_string();
            let alt_border_seq = format!("{}{}", ref_border_seq, &alt_allele[ref_allele.len()..]);
            annotate_microhomology(
                &mut annotation,
                &ref_border_seq,
                &alt_border_seq,
                ref_allele,
                alt_allele,
            );
        } else if leading_base_count == alt_allele.len() {
            let prefix_end = pos.saturating_sub(1).min(reference.len());
            let ref_border_seq = format!("{}{}", &reference[..prefix_end], ref_allele);
            let alt_border_seq = format!("{}{}", &ref_border_seq[..prefix_end], alt_allele);
            annotate_microhomology(
                &mut annotation,
                &ref_border_seq,
                &alt_border_seq,
                ref_allele,
                alt_allele,
            );
        }
    }

    annotation
}

fn annotate_microhomology(
    annotation: &mut VariantTypeAnnotation,
    ref_border_seq: &str,
    alt_border_seq: &str,
    ref_allele: &str,
    alt_allele: &str,
) {
    let microhomology_size =
        suffix_match_len(ref_border_seq.as_bytes(), alt_border_seq.as_bytes(), 99);
    annotation.indel_group = if microhomology_size <= 1 {
        "NHEJ".to_string()
    } else {
        "MMEJ".to_string()
    };
    let indel_size = alt_allele.len() as isize - ref_allele.len() as isize;
    annotation.summary_anno = format!("MH_size={microhomology_size};indel_size={indel_size}");
    annotation.annotation = "Yes";
}

#[cfg(test)]
fn combine_annotated_variants(
    single_rows: Vec<AnnotatedSingleVariant>,
    frequency_depth: &[usize],
) -> Vec<AnnotatedVariantCompatRow> {
    let mut by_pos = BTreeMap::<usize, Vec<AnnotatedSingleVariant>>::new();
    for row in single_rows {
        by_pos.entry(row.pos).or_default().push(row);
    }

    let mut rows = Vec::new();
    for (pos, variants) in by_pos {
        let mut snv_rows = Vec::new();
        let mut homopolymer_rows = Vec::new();
        let mut homodimer_rows = Vec::new();
        let mut tandem_rows = Vec::new();
        let mut other_rows = Vec::new();
        for variant in variants {
            let sample_info = variant_sample_info(&variant);
            let source = VariantCombineSourceRow {
                pos,
                ref_allele: variant.ref_allele,
                row_type: variant.variant_type.clone(),
                alt_allele: variant.alt_allele,
                sample_info,
                counts: variant.counts,
                id_list: variant.id_list,
            };
            if variant.variant_type == "SNV" {
                snv_rows.push(source);
            } else if variant.variant_type == "InDel" {
                match variant.type_annotation.indel_group.as_str() {
                    "homopolymer" => {
                        let mut source = source;
                        source.row_type = "InDel,homopolymer".to_string();
                        homopolymer_rows.push(source);
                    }
                    "homodimer" => {
                        let mut source = source;
                        source.row_type = "InDel,homodimer".to_string();
                        homodimer_rows.push(source);
                    }
                    "tandem" => {
                        let mut source = source;
                        source.row_type = "InDel,tandem".to_string();
                        tandem_rows.push(source);
                    }
                    group => {
                        let mut source = source;
                        source.row_type = format!("InDel,{group}");
                        other_rows.push(source);
                    }
                }
            }
        }

        push_combined_variant_group(&mut rows, &snv_rows, "SNV", frequency_depth);
        push_combined_variant_group(
            &mut rows,
            &homopolymer_rows,
            "InDel,homopolymer",
            frequency_depth,
        );
        push_combined_variant_group(
            &mut rows,
            &homodimer_rows,
            "InDel,homodimer",
            frequency_depth,
        );
        push_combined_variant_group(&mut rows, &tandem_rows, "InDel,tandem", frequency_depth);
        for source in other_rows {
            rows.push(finalize_compat_row(vec![source], None, frequency_depth));
        }
    }
    rows
}

#[cfg(test)]
fn push_combined_variant_group(
    rows: &mut Vec<AnnotatedVariantCompatRow>,
    source_rows: &[VariantCombineSourceRow],
    row_type: &str,
    frequency_depth: &[usize],
) {
    if !source_rows.is_empty() {
        rows.push(finalize_compat_row(
            source_rows.to_vec(),
            Some(row_type),
            frequency_depth,
        ));
    }
}

fn finalize_compat_row(
    source_rows: Vec<VariantCombineSourceRow>,
    row_type_override: Option<&str>,
    frequency_depth: &[usize],
) -> AnnotatedVariantCompatRow {
    let pos = source_rows[0].pos;
    let row_type = row_type_override
        .map(ToString::to_string)
        .unwrap_or_else(|| source_rows[0].row_type.clone());
    let alt_allele = source_rows
        .iter()
        .map(|row| row.alt_allele.as_str())
        .collect::<Vec<_>>()
        .join("#");
    let sample_info = source_rows
        .iter()
        .map(|row| row.sample_info.as_str())
        .collect::<Vec<_>>()
        .join("#");
    let counts = source_rows
        .iter()
        .map(|row| row.counts.to_string())
        .collect::<Vec<_>>()
        .join("#");
    let id_list = source_rows
        .iter()
        .map(|row| row.id_list.as_str())
        .collect::<Vec<_>>()
        .join("#");
    let total_count = source_rows.iter().map(|row| row.counts).sum::<usize>();
    let multi_allelic = if source_rows.len() > 1 {
        "multi-allelic".to_string()
    } else {
        "di-allelic".to_string()
    };
    let combined_info = combined_variant_info(&row_type, &sample_info, &counts, total_count);
    let mut ref_allele = source_rows[0].ref_allele.clone();
    let fixed_ref = fix_combined_ref(&mut ref_allele, &sample_info, &multi_allelic);
    let depth = frequency_depth.get(pos).copied().unwrap_or(0);
    let frequency = if depth == 0 {
        0.0
    } else {
        total_count as f64 / depth as f64
    };

    AnnotatedVariantCompatRow {
        pos,
        ref_allele,
        alt_allele,
        row_type,
        total_count,
        combined_info,
        multi_allelic,
        counts,
        id_list,
        fixed_ref,
        depth,
        frequency,
    }
}

fn variant_sample_info(variant: &AnnotatedSingleVariant) -> String {
    format!(
        "{},{},{},{},{},{}|{}",
        variant.ref_allele,
        variant.alt_allele,
        variant.variant_type,
        variant.type_annotation.trinucleotide_context,
        variant.type_annotation.snv_group,
        variant.type_annotation.indel_group,
        variant.type_annotation.summary_anno
    )
}

fn combined_variant_info(
    row_type: &str,
    sample_info: &str,
    counts: &str,
    total_count: usize,
) -> String {
    if row_type == "SNV" {
        let sample_infos = sample_info.split('#').collect::<Vec<_>>();
        let count_values = counts.split('#').collect::<Vec<_>>();
        if sample_infos.len() > 1 {
            let mut context = "";
            let mut groups = Vec::new();
            for (index, info) in sample_infos.iter().enumerate() {
                let fields = info.split(',').collect::<Vec<_>>();
                if fields.len() >= 6 {
                    context = fields[3];
                    groups.push(format!(
                        "{}={}",
                        fields[4],
                        count_values.get(index).copied().unwrap_or("0")
                    ));
                }
            }
            return format!("{context},{}", groups.join(","));
        }
        let fields = sample_info.split(',').collect::<Vec<_>>();
        if fields.len() >= 6 {
            return format!("{},{}={}", fields[3], fields[4], total_count);
        }
        return "-".to_string();
    }

    if matches!(
        row_type,
        "InDel,homopolymer" | "InDel,homodimer" | "InDel,tandem"
    ) {
        let sample_infos = sample_info.split('#').collect::<Vec<_>>();
        let count_values = counts.split('#').collect::<Vec<_>>();
        let mut subtype = "";
        let mut ref_size_str = "";
        let mut change_counts = Vec::new();
        for (index, info) in sample_infos.iter().enumerate() {
            let Some((_left, right)) = info.split_once('|') else {
                continue;
            };
            let parts = right.split(';').collect::<Vec<_>>();
            if parts.len() != 3 {
                continue;
            }
            subtype = parts[0];
            ref_size_str = parts[1];
            let ref_size = parts[1]
                .split_once('=')
                .and_then(|(_, value)| value.parse::<isize>().ok())
                .unwrap_or(0);
            let alt_size = parts[2]
                .split_once('=')
                .and_then(|(_, value)| value.parse::<isize>().ok())
                .unwrap_or(0);
            change_counts.push(format!(
                "{}:{}",
                alt_size - ref_size,
                count_values.get(index).copied().unwrap_or("0")
            ));
        }
        if !subtype.is_empty() {
            return format!("{subtype};{ref_size_str};{}", change_counts.join(","));
        }
        return "-".to_string();
    }

    if row_type == "InDel,MMEJ" || row_type == "InDel,NHEJ" {
        if let Some((_left, right)) = sample_info.split_once('|') {
            return right.to_string();
        }
    }

    if row_type != "SNV" {
        "-".to_string()
    } else {
        String::new()
    }
}

fn fix_combined_ref(ref_allele: &mut String, sample_info: &str, multi_allelic: &str) -> String {
    if multi_allelic == "di-allelic" {
        let ref_new = sample_info.split(',').next().unwrap_or("");
        if ref_new != ref_allele {
            *ref_allele = ref_new.to_string();
            "Yes".to_string()
        } else {
            "No".to_string()
        }
    } else {
        let ref_values = sample_info
            .split('#')
            .filter_map(|info| info.split(',').next())
            .collect::<BTreeSet<_>>();
        if ref_values.len() == 1 {
            let ref_new = ref_values.iter().next().copied().unwrap_or("");
            if ref_new != ref_allele {
                *ref_allele = ref_new.to_string();
                "Yes".to_string()
            } else {
                "No".to_string()
            }
        } else {
            "No".to_string()
        }
    }
}

fn first_sequence(path: &Path) -> Result<String, OrgraftError> {
    read_sequence_records(path)?
        .into_iter()
        .next()
        .map(|record| record.sequence)
        .ok_or_else(|| {
            OrgraftError::InvalidArgument(format!(
                "{} does not contain a FASTA sequence",
                path.display()
            ))
        })
}

fn read_variant_frequency_depth(path: &Path) -> Result<Vec<usize>, OrgraftError> {
    let file = File::open(path)?;
    let mut lines = BufReader::new(file).lines();
    let header = lines
        .next()
        .ok_or_else(|| OrgraftError::InvalidArgument(format!("{} is empty", path.display())))??;
    let columns = header.split('\t').collect::<Vec<_>>();
    let position_index = tsv_column(&columns, "position", path)?;
    let depth_index = tsv_column(&columns, "fl_depth", path)?;
    let mut depths = vec![0usize];
    for line in lines {
        let line = line?;
        if line.trim().is_empty() {
            continue;
        }
        let fields = line.split('\t').collect::<Vec<_>>();
        if fields.len() <= position_index || fields.len() <= depth_index {
            continue;
        }
        let position = parse_usize_value(fields[position_index], "coverage position")?;
        let depth = parse_usize_value(fields[depth_index], "fl_depth")?;
        if position >= depths.len() {
            depths.resize(position + 1, 0);
        }
        depths[position] = depth;
    }
    Ok(depths)
}

fn trinucleotide_context(reference: &str, pos: usize, ref_allele: &str) -> String {
    if reference.is_empty() || pos == 0 {
        return "-".to_string();
    }
    let bytes = reference.as_bytes();
    let len = bytes.len();
    let left = if pos == 1 {
        bytes[len - 1]
    } else {
        bytes[pos - 2]
    };
    let right = if pos == len { bytes[0] } else { bytes[pos] };
    format!("{}{}{}", left as char, ref_allele, right as char)
}

fn is_exact_repeat(sequence: &str, unit: &str) -> bool {
    !unit.is_empty()
        && sequence.len() % unit.len() == 0
        && sequence == unit.repeat(sequence.len() / unit.len())
}

fn gcd_usize(mut left: usize, mut right: usize) -> usize {
    while right != 0 {
        let next = left % right;
        left = right;
        right = next;
    }
    left
}

fn ref_prefix(value: &str, len: usize) -> &str {
    &value[..len.min(value.len())]
}

fn common_prefix_len(left: &str, right: &str) -> usize {
    left.as_bytes()
        .iter()
        .zip(right.as_bytes())
        .take_while(|(left, right)| left == right)
        .count()
}

fn common_suffix_len_after_prefix(left: &str, right: &str, prefix_len: usize) -> usize {
    let left_tail = &left.as_bytes()[prefix_len.min(left.len())..];
    let right_tail = &right.as_bytes()[prefix_len.min(right.len())..];
    let limit = left_tail.len().min(right_tail.len());
    let mut count = 0usize;
    while count < limit
        && left_tail[left_tail.len() - 1 - count] == right_tail[right_tail.len() - 1 - count]
    {
        count += 1;
    }
    count
}

fn suffix_match_len(left: &[u8], right: &[u8], max_len: usize) -> usize {
    let limit = max_len.min(left.len()).min(right.len());
    let mut count = 0usize;
    while count < limit && left[left.len() - 1 - count] == right[right.len() - 1 - count] {
        count += 1;
    }
    count
}

fn append_snv_failed_segments_to_log(
    path: &Path,
    segments: &[VariantSegment],
    results: &[SegmentCallResult],
    failed_segments: usize,
) -> Result<(), OrgraftError> {
    let result_by_segment = result_by_segment_id(results);
    let mut file = OpenOptions::new().create(true).append(true).open(path)?;
    writeln!(
        file,
        "### snv-indel-failed-segments:round_1 count={failed_segments} ###"
    )?;
    if failed_segments == 0 {
        writeln!(file)?;
        return Ok(());
    }
    writeln!(file, "read_id\tsegment_id\texit_status\tmessage\tstderr")?;
    for segment in segments {
        let Some(result) = result_by_segment.get(&segment.segment_id) else {
            continue;
        };
        if result.exit_status == 0 {
            continue;
        }
        writeln!(
            file,
            "{}\t{}\t{}\t{}\t{}",
            segment.read_id,
            segment.segment_id,
            result.exit_status,
            result.message,
            result.stderr.replace(['\t', '\n', '\r'], " "),
        )?;
    }
    writeln!(file)?;
    Ok(())
}

fn result_by_segment_id(results: &[SegmentCallResult]) -> HashMap<String, &SegmentCallResult> {
    results
        .iter()
        .map(|result| (result.segment.segment_id.clone(), result))
        .collect()
}

fn reads_with_calls(results: &[SegmentCallResult]) -> usize {
    results
        .iter()
        .filter(|result| !result.calls.is_empty())
        .map(|result| result.segment.read_id.as_str())
        .collect::<HashSet<_>>()
        .len()
}

fn variant_call_type(call: &VariantCall) -> &'static str {
    if call.ref_allele.len() == 1 && call.alt_allele.len() == 1 {
        "SNV"
    } else {
        "InDel"
    }
}

fn custom_variant_caller_command_template(
    reference: &Path,
    segment_count: usize,
    worker_count: usize,
) -> String {
    format!(
        "custom Rust CIGAR-diff caller ({} segments, {} workers; reference={})",
        segment_count,
        worker_count,
        reference.display(),
    )
}

fn run_sv_minimap2_batch(
    minimap2: &Path,
    options: &PolishOptions,
    reference_path: &Path,
    reads_path: &Path,
    stderr_path: &Path,
    commands: &mut Vec<CommandRecord>,
) -> Result<(HashMap<String, Vec<PafAlignment>>, usize), OrgraftError> {
    let mut stderr_file = OpenOptions::new()
        .create(true)
        .append(true)
        .open(stderr_path)?;
    writeln!(stderr_file, "### sv-eval:round_1 batch minimap2 stderr ###")?;
    let stderr_for_child = stderr_file.try_clone()?;
    let mut command = Command::new(minimap2);
    command.arg("-t").arg(options.threads.to_string());
    for option in SV_MINIMAP2_OPTIONS {
        command.arg(option);
    }
    command.arg(reference_path).arg(reads_path);
    let command_text = format!("{command:?}");
    let started = Instant::now();
    let mut child = command
        .stdout(Stdio::piped())
        .stderr(Stdio::from(stderr_for_child))
        .spawn()?;
    let stdout = child.stdout.take().ok_or_else(|| {
        OrgraftError::InvalidArgument("failed to capture SV minimap2 stdout".to_string())
    })?;
    let paf_result = parse_paf_by_read(BufReader::new(stdout));
    let status = child.wait()?;
    let elapsed_seconds = started.elapsed().as_secs_f64();
    let status_text = if status.success() { "ok" } else { "failed" };
    writeln!(
        OpenOptions::new().append(true).open(stderr_path)?,
        "### sv-eval:round_1 batch status={status_text} elapsed_seconds={elapsed_seconds:.3} workers={} ###\n",
        options.threads
    )?;
    commands.push(CommandRecord {
        timestamp: timestamp(),
        stage: "sv-eval-minimap2-batch",
        round: "round_1".to_string(),
        status: status_text,
        elapsed_seconds,
        stdout: "stream:minimap2-paf".to_string(),
        stderr: display_path(stderr_path),
        command: command_text,
    });
    if !status.success() {
        return Err(OrgraftError::InvalidArgument(format!(
            "sv-eval round_1 batch minimap2 failed; see {}",
            stderr_path.display()
        )));
    }
    Ok((paf_result?, options.threads))
}

fn run_plot_script(
    options: &PolishOptions,
    inputs: &ResolvedInputs,
    paths: &PolishPaths,
    sv_eval_report: Option<&SvEvalReport>,
    snv_indel_report: Option<&SnvIndelReport>,
    commands: &mut Vec<CommandRecord>,
) -> Result<(), OrgraftError> {
    let soft_paths = read_soft_paths(&inputs.soft_paths)?;
    let python = soft_paths
        .get("python")
        .filter(|path| !path.is_absolute() || path.exists())
        .cloned()
        .unwrap_or_else(|| PathBuf::from("python3"));
    let mut sv_command = Command::new(&python);
    sv_command.arg(&paths.round1_plot_script);
    sv_command
        .arg("--plot-dpi")
        .arg(options.plot_dpi.to_string())
        .arg("--plot-output-format")
        .arg(options.plot_output_format.as_str())
        .arg("--coverage-plot-rasterize")
        .arg(if options.coverage_plot_rasterize {
            "on"
        } else {
            "off"
        });
    if let Some(range) = &options.plot_range {
        sv_command.arg("--plot-range").arg(range.as_arg());
    }
    let manual_highlight = !options.sv_plot_highlight_subgroups.is_empty()
        || options.sv_plot_highlight_read_ids.is_some();
    let highlight_subgroups = if options.sv_plot_highlight_subgroups.is_empty() && !manual_highlight
    {
        sv_eval_report
            .map(|report| report.auto_highlight_subgroups.clone())
            .unwrap_or_default()
    } else {
        options.sv_plot_highlight_subgroups.clone()
    };
    if !highlight_subgroups.is_empty() {
        sv_command
            .arg("--sv-plot-highlight-subgroups")
            .arg(highlight_subgroups.join(","));
    }
    if let Some(path) = &options.sv_plot_highlight_read_ids {
        sv_command.arg("--sv-plot-highlight-read-ids").arg(path);
    }
    run_plot_command(
        sv_command,
        "plot_sv",
        paths,
        &paths.round1_sv_plots_dir,
        commands,
    )?;
    if snv_indel_report.is_some() {
        let mut snv_command = Command::new(&python);
        snv_command
            .arg(&paths.round1_snv_indel_plot_script)
            .arg("--plot-dpi")
            .arg(options.plot_dpi.to_string())
            .arg("--plot-output-format")
            .arg(options.plot_output_format.as_str())
            .arg("--snv-indel-plot-rasterize")
            .arg(if options.snv_indel_plot_rasterize {
                "on"
            } else {
                "off"
            })
            .arg("--plot-range")
            .arg(
                options
                    .plot_range
                    .as_ref()
                    .map(PlotRange::as_arg)
                    .unwrap_or_else(|| ".".to_string()),
            )
            .arg("--low-confidence-types")
            .arg(&options.snv_indel_plot_low_confidence)
            .arg("--low-min-reads")
            .arg(options.snv_indel_plot_low_min_reads.to_string())
            .arg("--low-min-fraction")
            .arg(format!("{}", options.snv_indel_plot_low_min_fraction))
            .arg("--high-risk-fraction")
            .arg(format!("{}", options.snv_indel_plot_high_risk_fraction));
        run_plot_command(
            snv_command,
            "plot_snv_indel",
            paths,
            &paths.round1_snv_indel_plots_dir,
            commands,
        )?;
    }
    Ok(())
}

fn run_plot_command(
    mut command: Command,
    stage: &'static str,
    paths: &PolishPaths,
    output_dir: &Path,
    commands: &mut Vec<CommandRecord>,
) -> Result<(), OrgraftError> {
    let command_text = format!("{command:?}");
    let started = Instant::now();
    let output = command.output();
    let elapsed_seconds = started.elapsed().as_secs_f64();
    let (status_text, stdout, stderr) = match output {
        Ok(output) => {
            let status_text = if output.status.success() {
                "ok"
            } else {
                "skipped"
            };
            (
                status_text,
                String::from_utf8_lossy(&output.stdout).trim().to_string(),
                String::from_utf8_lossy(&output.stderr).trim().to_string(),
            )
        }
        Err(error) => (
            "skipped",
            String::new(),
            format!("plot command could not be started: {error}"),
        ),
    };
    let mut stderr_file = OpenOptions::new()
        .create(true)
        .append(true)
        .open(&paths.external_stderr)?;
    writeln!(
        stderr_file,
        "### {stage}:round_1 status={status_text} elapsed_seconds={elapsed_seconds:.3} ###"
    )?;
    if !stdout.is_empty() {
        writeln!(stderr_file, "stdout:\n{stdout}")?;
    }
    if !stderr.is_empty() {
        writeln!(stderr_file, "stderr:\n{stderr}")?;
    }
    writeln!(stderr_file)?;
    commands.push(CommandRecord {
        timestamp: timestamp(),
        stage,
        round: "round_1".to_string(),
        status: status_text,
        elapsed_seconds,
        stdout: display_path(output_dir),
        stderr: display_path(&paths.external_stderr),
        command: command_text,
    });
    Ok(())
}

fn parse_paf_by_read<R: BufRead>(
    reader: R,
) -> Result<HashMap<String, Vec<PafAlignment>>, OrgraftError> {
    let mut by_read: HashMap<String, Vec<PafAlignment>> = HashMap::new();
    for line in reader.lines() {
        let line = line?;
        if line.trim().is_empty() {
            continue;
        }
        let Some(alignment) = parse_paf_alignment(&line)? else {
            continue;
        };
        by_read
            .entry(alignment.query_id.clone())
            .or_default()
            .push(alignment);
    }
    Ok(by_read)
}

fn parse_paf_alignment(line: &str) -> Result<Option<PafAlignment>, OrgraftError> {
    let fields = line.split('\t').collect::<Vec<_>>();
    if fields.len() < 12 {
        return Ok(None);
    }
    Ok(Some(PafAlignment {
        query_id: normalize_read_id(fields[0]),
        query_start: parse_usize_value(fields[2], "PAF query start")?,
        query_end: parse_usize_value(fields[3], "PAF query end")?,
        strand: fields[4].chars().next().unwrap_or('+'),
        target_id: fields[5].to_string(),
        target_start: parse_usize_value(fields[7], "PAF target start")?,
        target_end: parse_usize_value(fields[8], "PAF target end")?,
        matches: parse_usize_value(fields[9], "PAF matches")?,
        block_len: parse_usize_value(fields[10], "PAF block length")?,
        mapq: fields[11].to_string(),
        alignment_role: paf_alignment_role(&fields),
        cigar: paf_tag_value(&fields, "cg").unwrap_or(".").to_string(),
    }))
}

fn paf_tag_value<'a>(fields: &'a [&str], tag: &str) -> Option<&'a str> {
    fields.iter().skip(12).find_map(|field| {
        let mut parts = field.splitn(3, ':');
        let current_tag = parts.next()?;
        let _kind = parts.next()?;
        let value = parts.next()?;
        (current_tag == tag).then_some(value)
    })
}

fn paf_alignment_role(fields: &[&str]) -> String {
    match paf_tag_value(fields, "tp").unwrap_or(".") {
        "P" => "primary",
        "S" => "secondary",
        "I" | "i" => "inversion",
        other => other,
    }
    .to_string()
}

fn write_whole_read_evidence_header<W: Write>(file: &mut W) -> Result<(), OrgraftError> {
    writeln!(
        file,
        "read_id\tread_length\tquery_start\tquery_end\ttarget_id\ttarget_start\ttarget_end\tstrand\tidentity\tmapq\talignment_role\tcigar\tcrosses_junction\tjunction_id\tcrosses_repeat_choice\trepeat_choice_id\talignment_summary"
    )?;
    Ok(())
}

fn write_whole_read_evidence_rows(
    file: &mut impl Write,
    read: &ReadRecord,
    paf_lines: &[PafAlignment],
    alignment_summary: Option<&str>,
) -> Result<usize, OrgraftError> {
    let summary_cell = alignment_summary
        .map(escape_tsv_cell)
        .unwrap_or_else(|| ".".to_string());
    if paf_lines.is_empty() {
        writeln!(
            file,
            "{}\t{}\t.\t.\t.\t.\t.\t.\t.\t.\tno_alignment\t.\tnot_evaluated\t.\tnot_evaluated\t.\t{}",
            read.id,
            read.sequence.len(),
            summary_cell,
        )?;
        return Ok(1);
    }

    for (index, alignment) in paf_lines.iter().enumerate() {
        let identity = (alignment.matches as f64) / (alignment.block_len.max(1) as f64) * 100.0;
        let summary_cell = if index == 0 {
            summary_cell.as_str()
        } else {
            "."
        };
        writeln!(
            file,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.3}\t{}\t{}\t{}\tnot_evaluated\t.\tnot_evaluated\t.\t{}",
            alignment.query_id,
            read.sequence.len(),
            alignment.query_start + 1,
            alignment.query_end,
            alignment.target_id,
            alignment.target_start + 1,
            alignment.target_end,
            alignment.strand,
            identity,
            alignment.mapq,
            alignment.alignment_role,
            alignment.cigar,
            summary_cell,
        )?;
    }
    Ok(paf_lines.len())
}

fn extend_paf_terminal_microindels(
    paf_lines: &[PafAlignment],
    read_sequence: &str,
    reference_by_id: &HashMap<String, String>,
) -> Vec<PafAlignment> {
    let read_sequence = read_sequence.to_ascii_uppercase();
    let mut adjusted = Vec::with_capacity(paf_lines.len());
    for alignment in paf_lines {
        let Some(target_sequence) = reference_by_id.get(&alignment.target_id) else {
            adjusted.push(alignment.clone());
            continue;
        };
        if alignment.query_end.saturating_sub(alignment.query_start)
            < TERMINAL_EXTENSION_MIN_ALIGNMENT_LENGTH
        {
            adjusted.push(alignment.clone());
            continue;
        }

        let query_right = slice_forward(
            &read_sequence,
            alignment.query_end,
            TERMINAL_EXTENSION_WINDOW,
        );
        let query_left = reverse_string(slice_before(
            &read_sequence,
            alignment.query_start,
            TERMINAL_EXTENSION_WINDOW,
        ));
        let (target_right, target_left) = if alignment.strand == '+' {
            (
                slice_forward(
                    target_sequence,
                    alignment.target_end,
                    TERMINAL_EXTENSION_WINDOW,
                )
                .to_string(),
                reverse_string(slice_before(
                    target_sequence,
                    alignment.target_start,
                    TERMINAL_EXTENSION_WINDOW,
                )),
            )
        } else {
            (
                reverse_complement(slice_before(
                    target_sequence,
                    alignment.target_start,
                    TERMINAL_EXTENSION_WINDOW,
                )),
                reverse_string(&reverse_complement(slice_forward(
                    target_sequence,
                    alignment.target_end,
                    TERMINAL_EXTENSION_WINDOW,
                ))),
            )
        };

        let (right_query, right_target, right_matches) =
            best_terminal_microindel_extension(&query_right, &target_right);
        let (left_query, left_target, left_matches) =
            best_terminal_microindel_extension(&query_left, &target_left);

        let mut next = alignment.clone();
        next.query_start = next.query_start.saturating_sub(left_query);
        next.query_end += right_query;
        if next.strand == '+' {
            next.target_start = next.target_start.saturating_sub(left_target);
            next.target_end += right_target;
        } else {
            next.target_start = next.target_start.saturating_sub(right_target);
            next.target_end += left_target;
        }
        next.matches += left_matches + right_matches;
        next.block_len += left_query.max(left_target) + right_query.max(right_target);
        next.block_len = next.block_len.max(1);
        adjusted.push(next);
    }
    adjusted
}

fn paf_to_blast_like(paf_lines: &[PafAlignment]) -> Vec<BlastLikeAlignment> {
    paf_lines
        .iter()
        .map(|alignment| {
            let block_len = alignment.block_len.max(1);
            let pident = (alignment.matches as f64) / (block_len as f64) * 100.0;
            let (subject_start, subject_end) = if alignment.strand == '-' {
                (alignment.target_end, alignment.target_start + 1)
            } else {
                (alignment.target_start + 1, alignment.target_end)
            };
            BlastLikeAlignment {
                query_id: alignment.query_id.clone(),
                subject_id: alignment.target_id.clone(),
                pident: round_to_three_decimals(pident),
                query_start: alignment.query_start + 1,
                query_end: alignment.query_end,
                subject_start,
                subject_end,
            }
        })
        .collect()
}

fn build_sorted_alignment_summary(
    blast_like: &[BlastLikeAlignment],
    read_len: usize,
) -> Option<(String, AlignmentSummaryRecord)> {
    let one_alignments = sorted_noncontained_alignments(blast_like);
    let first = one_alignments.first()?;
    let num_align = one_alignments.len();
    let mut part_lengths = Vec::with_capacity(num_align);
    let mut overlap_lengths = Vec::with_capacity(num_align);
    let mut summary_alignments = Vec::with_capacity(num_align);
    for alignment in &one_alignments {
        part_lengths.push(alignment.query_end.saturating_sub(alignment.query_start) + 1);
    }
    let mut subtype_parts = Vec::with_capacity(num_align);
    for index in 0..num_align {
        if index + 1 == num_align {
            overlap_lengths.push(None);
        } else {
            let value = (one_alignments[index].query_end as isize)
                - (one_alignments[index + 1].query_start as isize)
                + 1;
            overlap_lengths.push(Some(value));
        }
        let label = match overlap_lengths[index] {
            None => "NA",
            Some(0) => "ref",
            Some(value) if value > 0 => "rep",
            Some(_) => "ins",
        };
        subtype_parts.push(label);
    }
    let covered = part_lengths.iter().sum::<usize>() as isize
        - overlap_lengths.iter().flatten().sum::<isize>();
    let percent_total = if read_len == 0 {
        0.0
    } else {
        (covered as f64) / (read_len as f64) * 100.0
    };

    let subtype_csv = subtype_parts.join(",");
    let aln_type = format!("aln_type={num_align};{subtype_csv}");

    let mut row = format!(
        "{}\t{}\t{}\t{}\t{}",
        first.query_id,
        first.subject_id,
        read_len,
        aln_type,
        format_py_float(percent_total),
    );
    for (index, alignment) in one_alignments.iter().enumerate() {
        let strand = if alignment.subject_start > alignment.subject_end {
            '-'
        } else {
            '+'
        };
        let overlap = overlap_lengths[index]
            .map(|value| value.to_string())
            .unwrap_or_else(|| "NA".to_string());
        row.push('\t');
        row.push_str(&format!(
            "aln={};len={};olp={};idt={};strand={};qs={};qe={};ss={};se={};cn=1;c1=100.0,1,1",
            index + 1,
            part_lengths[index],
            overlap,
            format_py_float(alignment.pident),
            strand,
            alignment.query_start,
            alignment.query_end,
            alignment.subject_start,
            alignment.subject_end,
        ));
        summary_alignments.push(SummaryAlignment {
            olp: overlap_lengths[index],
            strand,
            qs: alignment.query_start,
            qe: alignment.query_end,
            ss: alignment.subject_start as isize,
            se: alignment.subject_end as isize,
            cn: 1,
        });
    }
    let record = AlignmentSummaryRecord {
        read_id: first.query_id.clone(),
        target_id: first.subject_id.clone(),
        read_len,
        num_align,
        subtype: subtype_csv.replace(',', "_"),
        percent_total,
        alignments: summary_alignments,
    };
    Some((row, record))
}

#[derive(Debug, Clone)]
struct AlignmentSummaryRecord {
    read_id: String,
    target_id: String,
    read_len: usize,
    num_align: usize,
    subtype: String,
    percent_total: f64,
    alignments: Vec<SummaryAlignment>,
}

impl AlignmentSummaryRecord {
    fn group_name(&self) -> String {
        format!("type_{}_subtype_{}", self.num_align, self.subtype)
    }

    fn is_fl(&self) -> bool {
        self.percent_total >= FL_PERCENT_TOTAL_THRESHOLD
    }

    fn is_multi(&self) -> bool {
        self.alignments.iter().any(|alignment| alignment.cn > 1)
    }

    fn subgroup_key(&self) -> Option<SubgroupKey> {
        if !self.is_fl() || self.alignments.len() < 2 {
            return None;
        }
        let mut pairs = Vec::with_capacity(self.alignments.len().saturating_sub(1));
        for index in 0..self.alignments.len() - 1 {
            pairs.push((self.alignments[index].se, self.alignments[index + 1].ss));
        }
        Some(SubgroupKey(pairs))
    }
}

#[derive(Debug, Clone)]
struct SummaryAlignment {
    olp: Option<isize>,
    strand: char,
    qs: usize,
    qe: usize,
    ss: isize,
    se: isize,
    cn: usize,
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
struct SubgroupKey(Vec<(isize, isize)>);

impl SubgroupKey {
    fn label(&self) -> String {
        self.0
            .iter()
            .enumerate()
            .map(|(index, (se, ss))| format!("se{}={},ss{}={}", index + 1, se, index + 2, ss))
            .collect::<Vec<_>>()
            .join(";")
    }
}

#[derive(Debug, Clone)]
struct ReadGroupStats {
    num_align: usize,
    subtype: String,
    total_reads: usize,
    fl_reads: usize,
    partial_reads: usize,
    fl_multi_reads: usize,
}

impl ReadGroupStats {
    fn new(record: &AlignmentSummaryRecord) -> Self {
        Self {
            num_align: record.num_align,
            subtype: record.subtype.clone(),
            total_reads: 0,
            fl_reads: 0,
            partial_reads: 0,
            fl_multi_reads: 0,
        }
    }

    fn add(&mut self, record: &AlignmentSummaryRecord) {
        self.total_reads += 1;
        if record.is_fl() {
            self.fl_reads += 1;
            if record.is_multi() {
                self.fl_multi_reads += 1;
            }
        } else {
            self.partial_reads += 1;
        }
    }
}

fn write_read_group_reports(
    records: &[AlignmentSummaryRecord],
    paths: &PolishPaths,
    reference_len: usize,
    options: &PolishOptions,
    paf_by_read: &HashMap<String, Vec<PafAlignment>>,
) -> Result<ReadGroupReport, OrgraftError> {
    let mut groups: BTreeMap<String, ReadGroupStats> = BTreeMap::new();
    let mut subgroups: BTreeMap<String, BTreeMap<SubgroupKey, Vec<usize>>> = BTreeMap::new();
    for (record_index, record) in records.iter().enumerate() {
        let group_name = record.group_name();
        groups
            .entry(group_name.clone())
            .or_insert_with(|| ReadGroupStats::new(record))
            .add(record);
        if let Some(key) = record.subgroup_key() {
            subgroups
                .entry(group_name)
                .or_default()
                .entry(key)
                .or_default()
                .push(record_index);
        }
    }

    let mut subgroup_old_index: BTreeMap<String, BTreeMap<SubgroupKey, usize>> = BTreeMap::new();
    for (group_name, group_subgroups) in &subgroups {
        let mut keys = group_subgroups.keys().cloned().collect::<Vec<_>>();
        keys.sort_by(subgroup_old_index_order);
        let mut indices = BTreeMap::new();
        for (index, key) in keys.into_iter().enumerate() {
            indices.insert(key, index + 1);
        }
        subgroup_old_index.insert(group_name.clone(), indices);
    }

    let sorted_groups = sorted_read_groups(&groups);
    write_read_group_summary(&paths.round1_sv_group_summary, &sorted_groups)?;
    let read_subgroup_count = write_read_subgroup_summary(
        &paths.round1_sv_subgroup_summary,
        &records,
        &groups,
        &subgroups,
        &subgroup_old_index,
    )?;
    write_read_group_ids(&paths.round1_sv_group_ids, &records, &subgroup_old_index)?;
    let support_report = write_sv_support_outputs(
        options,
        paths,
        reference_len,
        &records,
        &groups,
        &subgroups,
        &subgroup_old_index,
        paf_by_read,
    )?;

    Ok(ReadGroupReport {
        fl_reads: sorted_groups.iter().map(|(_, stats)| stats.fl_reads).sum(),
        partial_reads: sorted_groups
            .iter()
            .map(|(_, stats)| stats.partial_reads)
            .sum(),
        reference_support_reads: support_report.reference_support_reads,
        read_group_count: sorted_groups.len(),
        read_subgroup_count,
        sv_support_status: support_report.status,
        auto_highlight_subgroups: support_report.auto_highlight_subgroups,
    })
}

#[derive(Debug, Clone)]
struct SvSupportReport {
    reference_support_reads: usize,
    status: String,
    auto_highlight_subgroups: Vec<String>,
}

#[derive(Debug, Clone)]
struct CoverageData {
    fl: Vec<u32>,
    partial: Vec<u32>,
    reference_support: Vec<u32>,
}

#[derive(Debug, Clone, Copy)]
struct WindowDepthStats {
    mean_fl: f64,
    mean_partial: f64,
    mean_reference_support: f64,
    reference_fraction: f64,
}

fn write_sv_support_outputs(
    options: &PolishOptions,
    paths: &PolishPaths,
    reference_len: usize,
    records: &[AlignmentSummaryRecord],
    groups: &BTreeMap<String, ReadGroupStats>,
    subgroups: &BTreeMap<String, BTreeMap<SubgroupKey, Vec<usize>>>,
    subgroup_old_index: &BTreeMap<String, BTreeMap<SubgroupKey, usize>>,
    paf_by_read: &HashMap<String, Vec<PafAlignment>>,
) -> Result<SvSupportReport, OrgraftError> {
    let read_classes = read_class_by_id(records);
    let reference_support_ids =
        reference_support_read_ids(records, subgroups, subgroup_old_index, reference_len);
    let coverage = coverage_from_paf_by_read(
        paf_by_read,
        reference_len,
        &read_classes,
        &reference_support_ids,
    )?;
    write_coverage_tsv(&paths.round1_sv_coverage, &coverage)?;
    let total_fl = records.iter().filter(|record| record.is_fl()).count();
    let auto_highlight_subgroups = write_high_subgroup_report(
        &paths.round1_sv_high_subgroup_report,
        reference_len,
        total_fl,
        options.sv_plot_highlight_min_fraction,
        options.sv_plot_highlight_min_reads,
        records,
        groups,
        subgroups,
        subgroup_old_index,
        &reference_support_ids,
        &coverage,
    )?;
    let status = write_sv_support_summary(
        &paths.round1_sv_support_summary,
        reference_len,
        total_fl,
        records.len().saturating_sub(total_fl),
        reference_support_ids.len(),
        &coverage,
        options.sv_plot_highlight_min_fraction,
        options.sv_plot_highlight_min_reads,
        &auto_highlight_subgroups,
    )?;
    write_plot_script(paths)?;

    Ok(SvSupportReport {
        reference_support_reads: reference_support_ids.len(),
        status,
        auto_highlight_subgroups,
    })
}

fn read_class_by_id(records: &[AlignmentSummaryRecord]) -> HashMap<String, &'static str> {
    records
        .iter()
        .map(|record| {
            (
                record.read_id.clone(),
                if record.is_fl() { "FL" } else { "partial" },
            )
        })
        .collect()
}

fn reference_support_read_ids(
    records: &[AlignmentSummaryRecord],
    subgroups: &BTreeMap<String, BTreeMap<SubgroupKey, Vec<usize>>>,
    subgroup_old_index: &BTreeMap<String, BTreeMap<SubgroupKey, usize>>,
    reference_len: usize,
) -> HashSet<String> {
    let mut ids = HashSet::new();
    for record in records {
        if record.is_fl() && record.group_name() == "type_1_subtype_NA" {
            ids.insert(record.read_id.clone());
        }
    }
    for (group_name, group_subgroups) in subgroups {
        for (key, record_indices) in group_subgroups {
            let is_reference_subgroup = match group_name.as_str() {
                "type_2_subtype_ref_NA" => is_circular_terminal_subgroup(key, reference_len),
                "type_2_subtype_rep_NA" => {
                    subgroup_mid_olp(records, record_indices, 0)
                        >= REFERENCE_SUPPORT_REP_MID_OLP_MIN
                }
                _ => false,
            };
            if !is_reference_subgroup {
                continue;
            }
            if subgroup_old_index
                .get(group_name)
                .and_then(|indices| indices.get(key))
                .is_none()
            {
                continue;
            }
            for record_index in record_indices {
                ids.insert(records[*record_index].read_id.clone());
            }
        }
    }
    ids
}

fn is_circular_terminal_subgroup(key: &SubgroupKey, reference_len: usize) -> bool {
    matches!(
        key.0.as_slice(),
        [(1, ss2)] if *ss2 == reference_len as isize
    ) || matches!(
        key.0.as_slice(),
        [(se1, 1)] if *se1 == reference_len as isize
    )
}

fn subgroup_mid_olp(
    records: &[AlignmentSummaryRecord],
    record_indices: &[usize],
    align_index: usize,
) -> f64 {
    let values = record_indices
        .iter()
        .filter_map(|record_index| records[*record_index].alignments.get(align_index)?.olp)
        .collect::<Vec<_>>();
    if values.is_empty() {
        return 0.0;
    }
    let min = values.iter().min().copied().unwrap_or(0);
    let max = values.iter().max().copied().unwrap_or(0);
    (min as f64 + max as f64) / 2.0
}

fn coverage_from_paf_by_read(
    paf_by_read: &HashMap<String, Vec<PafAlignment>>,
    reference_len: usize,
    read_classes: &HashMap<String, &'static str>,
    reference_support_ids: &HashSet<String>,
) -> Result<CoverageData, OrgraftError> {
    let mut fl_diff = vec![0i32; reference_len + 1];
    let mut partial_diff = vec![0i32; reference_len + 1];
    let mut reference_support_diff = vec![0i32; reference_len + 1];
    for (read_id, alignments) in paf_by_read {
        let Some(read_class) = read_classes.get(read_id).copied() else {
            continue;
        };
        let is_reference_support = reference_support_ids.contains(read_id);
        for alignment in alignments {
            if alignment.alignment_role == "secondary" {
                continue;
            }
            match read_class {
                "FL" => {
                    add_paf_alignment_coverage(&mut fl_diff, alignment, reference_len)?;
                    if is_reference_support {
                        add_paf_alignment_coverage(
                            &mut reference_support_diff,
                            alignment,
                            reference_len,
                        )?;
                    }
                }
                "partial" => {
                    add_paf_alignment_coverage(&mut partial_diff, alignment, reference_len)?;
                }
                _ => {}
            }
        }
    }
    Ok(CoverageData {
        fl: diff_to_coverage(fl_diff),
        partial: diff_to_coverage(partial_diff),
        reference_support: diff_to_coverage(reference_support_diff),
    })
}

fn add_paf_alignment_coverage(
    diff: &mut [i32],
    alignment: &PafAlignment,
    reference_len: usize,
) -> Result<(), OrgraftError> {
    if alignment.cigar == "." {
        add_coverage_interval(
            diff,
            alignment.target_start + 1,
            alignment.target_end,
            reference_len,
        );
        return Ok(());
    }
    let mut reference_pos = alignment.target_start + 1;
    for (len, op) in parse_cigar(&alignment.cigar)? {
        match op {
            'M' | '=' | 'X' | 'D' | 'N' => {
                let end = reference_pos + len.saturating_sub(1);
                add_coverage_interval(diff, reference_pos, end, reference_len);
                reference_pos += len;
            }
            'I' | 'S' | 'H' | 'P' => {}
            _ => {}
        }
    }
    Ok(())
}

fn add_coverage_interval(diff: &mut [i32], start: usize, end: usize, reference_len: usize) {
    if reference_len == 0 {
        return;
    }
    let start = start.min(reference_len).max(1);
    let end = end.min(reference_len).max(start);
    diff[start - 1] += 1;
    diff[end] -= 1;
}

fn diff_to_coverage(diff: Vec<i32>) -> Vec<u32> {
    let mut coverage = Vec::with_capacity(diff.len().saturating_sub(1));
    let mut current = 0i32;
    for delta in diff.into_iter().take(coverage.capacity()) {
        current += delta;
        coverage.push(current.max(0) as u32);
    }
    coverage
}

fn write_coverage_tsv(path: &Path, coverage: &CoverageData) -> Result<(), OrgraftError> {
    let mut buffer = String::with_capacity(coverage.fl.len().saturating_mul(32));
    buffer.push_str(
        "position\tfl_depth\tpartial_depth\treference_support_depth\ttotal_depth\treference_support_fraction\n",
    );
    for index in 0..coverage.fl.len() {
        let fl = coverage.fl[index];
        let partial = coverage.partial[index];
        let reference_support = coverage.reference_support[index];
        let fraction = if fl == 0 {
            0.0
        } else {
            reference_support as f64 / fl as f64
        };
        writeln!(
            buffer,
            "{}\t{}\t{}\t{}\t{}\t{:.6}",
            index + 1,
            fl,
            partial,
            reference_support,
            fl + partial,
            fraction
        )
        .expect("writing to String cannot fail");
    }
    fs::write(path, buffer)?;
    Ok(())
}

fn write_sv_support_summary(
    path: &Path,
    reference_len: usize,
    fl_reads: usize,
    partial_reads: usize,
    reference_support_reads: usize,
    coverage: &CoverageData,
    auto_sv_plot_highlight_min_fraction: f64,
    auto_sv_plot_highlight_min_reads: usize,
    auto_highlight_subgroups: &[String],
) -> Result<String, OrgraftError> {
    let total_fl_depth = coverage.fl.iter().map(|value| *value as u64).sum::<u64>();
    let total_reference_depth = coverage
        .reference_support
        .iter()
        .map(|value| *value as u64)
        .sum::<u64>();
    let reference_area_fraction = if total_fl_depth == 0 {
        0.0
    } else {
        total_reference_depth as f64 / total_fl_depth as f64
    };
    let low_green_window_fraction = low_green_window_fraction(coverage);
    let green_median = median_u32(&coverage.reference_support);
    let fl_median = median_u32(&coverage.fl);
    let status = if reference_area_fraction >= SV_SUPPORT_MIN_GREEN_FRACTION
        && low_green_window_fraction <= SV_SUPPORT_MAX_LOW_GREEN_WINDOW_FRACTION
    {
        "pass"
    } else {
        "review"
    }
    .to_string();

    let mut file = File::create(path)?;
    writeln!(file, "metric\tvalue")?;
    writeln!(file, "status\t{status}")?;
    writeln!(file, "reference_length\t{reference_len}")?;
    writeln!(file, "fl_reads\t{fl_reads}")?;
    writeln!(file, "partial_reads\t{partial_reads}")?;
    writeln!(file, "reference_support_reads\t{reference_support_reads}")?;
    writeln!(
        file,
        "reference_support_read_fraction\t{:.6}",
        if fl_reads == 0 {
            0.0
        } else {
            reference_support_reads as f64 / fl_reads as f64
        }
    )?;
    writeln!(
        file,
        "reference_support_depth_area_fraction\t{reference_area_fraction:.6}"
    )?;
    writeln!(
        file,
        "low_green_window_fraction\t{low_green_window_fraction:.6}"
    )?;
    writeln!(file, "fl_depth_median\t{fl_median:.3}")?;
    writeln!(file, "reference_support_depth_median\t{green_median:.3}")?;
    writeln!(
        file,
        "auto_sv_plot_highlight_min_fraction\t{auto_sv_plot_highlight_min_fraction:.6}"
    )?;
    writeln!(
        file,
        "auto_sv_plot_highlight_min_reads\t{auto_sv_plot_highlight_min_reads}"
    )?;
    writeln!(
        file,
        "auto_highlight_subgroup_count\t{}",
        auto_highlight_subgroups.len()
    )?;
    writeln!(
        file,
        "auto_highlight_subgroups\t{}",
        if auto_highlight_subgroups.is_empty() {
            ".".to_string()
        } else {
            auto_highlight_subgroups.join(",")
        }
    )?;
    writeln!(
        file,
        "min_green_fraction_threshold\t{SV_SUPPORT_MIN_GREEN_FRACTION:.6}"
    )?;
    writeln!(
        file,
        "max_low_green_window_fraction_threshold\t{SV_SUPPORT_MAX_LOW_GREEN_WINDOW_FRACTION:.6}"
    )?;
    writeln!(file, "window_bp\t{SV_SUPPORT_WINDOW_BP}")?;
    Ok(status)
}

fn median_u32(values: &[u32]) -> f64 {
    if values.is_empty() {
        return 0.0;
    }
    let mut sorted = values.to_vec();
    sorted.sort_unstable();
    let mid = sorted.len() / 2;
    if sorted.len() % 2 == 0 {
        (sorted[mid - 1] as f64 + sorted[mid] as f64) / 2.0
    } else {
        sorted[mid] as f64
    }
}

fn low_green_window_fraction(coverage: &CoverageData) -> f64 {
    let mut total_windows = 0usize;
    let mut low_windows = 0usize;
    let mut start = 0usize;
    while start < coverage.fl.len() {
        let end = (start + SV_SUPPORT_WINDOW_BP).min(coverage.fl.len());
        let stats = depth_stats_range(coverage, start, end);
        if stats.mean_fl > 0.0 {
            total_windows += 1;
            if stats.mean_reference_support < SV_SUPPORT_MIN_GREEN_DEPTH
                || stats.reference_fraction < SV_SUPPORT_LOW_GREEN_FRACTION
            {
                low_windows += 1;
            }
        }
        start = end;
    }
    if total_windows == 0 {
        1.0
    } else {
        low_windows as f64 / total_windows as f64
    }
}

fn depth_stats_range(coverage: &CoverageData, start: usize, end: usize) -> WindowDepthStats {
    let end = end.min(coverage.fl.len()).max(start);
    let len = end.saturating_sub(start).max(1) as f64;
    let sum_fl = coverage.fl[start..end]
        .iter()
        .map(|value| *value as u64)
        .sum::<u64>();
    let sum_partial = coverage.partial[start..end]
        .iter()
        .map(|value| *value as u64)
        .sum::<u64>();
    let sum_reference = coverage.reference_support[start..end]
        .iter()
        .map(|value| *value as u64)
        .sum::<u64>();
    let mean_fl = sum_fl as f64 / len;
    let mean_reference_support = sum_reference as f64 / len;
    WindowDepthStats {
        mean_fl,
        mean_partial: sum_partial as f64 / len,
        mean_reference_support,
        reference_fraction: if mean_fl == 0.0 {
            0.0
        } else {
            mean_reference_support / mean_fl
        },
    }
}

fn depth_stats_around(
    coverage: &CoverageData,
    center_1_based: isize,
    window_bp: usize,
) -> WindowDepthStats {
    if coverage.fl.is_empty() {
        return WindowDepthStats {
            mean_fl: 0.0,
            mean_partial: 0.0,
            mean_reference_support: 0.0,
            reference_fraction: 0.0,
        };
    }
    let center = center_1_based.max(1) as usize;
    let start = center.saturating_sub(window_bp + 1);
    let end = (center + window_bp).min(coverage.fl.len());
    depth_stats_range(coverage, start, end)
}

fn write_high_subgroup_report(
    path: &Path,
    reference_len: usize,
    total_fl: usize,
    auto_sv_plot_highlight_min_fraction: f64,
    auto_sv_plot_highlight_min_reads: usize,
    records: &[AlignmentSummaryRecord],
    groups: &BTreeMap<String, ReadGroupStats>,
    subgroups: &BTreeMap<String, BTreeMap<SubgroupKey, Vec<usize>>>,
    subgroup_old_index: &BTreeMap<String, BTreeMap<SubgroupKey, usize>>,
    reference_support_ids: &HashSet<String>,
    coverage: &CoverageData,
) -> Result<Vec<String>, OrgraftError> {
    let mut file = File::create(path)?;
    writeln!(
        file,
        "group_name\told_index\tboundary_key\tsubgroup_reads\tsubgroup_fraction\tis_reference_support_subgroup\tauto_highlight_default\tmin_window_fl_depth\tmin_window_reference_support_depth\tmin_window_reference_support_fraction\tmax_window_partial_depth\tjudgement"
    )?;
    let mut auto_highlight_subgroups = Vec::new();
    for (group_name, group_subgroups) in subgroups {
        let Some(stats) = groups.get(group_name) else {
            continue;
        };
        for (key, record_indices) in group_subgroups {
            let subgroup_fraction = if total_fl == 0 {
                0.0
            } else {
                record_indices.len() as f64 / total_fl as f64
            };
            let old_index = subgroup_old_index
                .get(group_name)
                .and_then(|indices| indices.get(key))
                .copied()
                .unwrap_or(0);
            let is_reference_support_subgroup = record_indices
                .iter()
                .all(|index| reference_support_ids.contains(&records[*index].read_id))
                && !record_indices.is_empty();
            let auto_highlight_default = !is_reference_support_subgroup
                && record_indices.len() >= auto_sv_plot_highlight_min_reads
                && subgroup_fraction >= auto_sv_plot_highlight_min_fraction;
            if !auto_highlight_default
                && record_indices.len() < READ_SUBGROUP_WATCH_MIN
                && subgroup_fraction < HIGH_SUBGROUP_MIN_FRACTION
            {
                continue;
            }
            if auto_highlight_default && old_index > 0 {
                auto_highlight_subgroups.push(format!("{group_name}:{old_index}"));
            }
            let mut min_fl = f64::MAX;
            let mut min_reference = f64::MAX;
            let mut min_reference_fraction = f64::MAX;
            let mut max_partial: f64 = 0.0;
            for (se, ss) in &key.0 {
                for center in [*se, *ss] {
                    let window = depth_stats_around(coverage, center, BREAKPOINT_WINDOW_BP);
                    min_fl = min_fl.min(window.mean_fl);
                    min_reference = min_reference.min(window.mean_reference_support);
                    min_reference_fraction = min_reference_fraction.min(window.reference_fraction);
                    max_partial = max_partial.max(window.mean_partial);
                }
            }
            let judgement = high_subgroup_judgement(
                stats,
                is_reference_support_subgroup,
                min_reference,
                min_reference_fraction,
            );
            writeln!(
                file,
                "{}\t{}\t{}\t{}\t{:.6}\t{}\t{}\t{:.3}\t{:.3}\t{:.6}\t{:.3}\t{}",
                group_name,
                old_index,
                key.label(),
                record_indices.len(),
                subgroup_fraction,
                is_reference_support_subgroup,
                auto_highlight_default,
                if min_fl == f64::MAX { 0.0 } else { min_fl },
                if min_reference == f64::MAX {
                    0.0
                } else {
                    min_reference
                },
                if min_reference_fraction == f64::MAX {
                    0.0
                } else {
                    min_reference_fraction
                },
                max_partial,
                judgement,
            )?;
        }
    }
    let _ = reference_len;
    Ok(auto_highlight_subgroups)
}

fn high_subgroup_judgement(
    stats: &ReadGroupStats,
    is_reference_support_subgroup: bool,
    min_reference_depth: f64,
    min_reference_fraction: f64,
) -> &'static str {
    if is_reference_support_subgroup {
        "reference_support_configuration"
    } else if stats.num_align >= 3
        && (min_reference_depth < SV_SUPPORT_MIN_GREEN_DEPTH
            || min_reference_fraction < SV_SUPPORT_LOW_GREEN_FRACTION)
    {
        "possible_reference_sv_error"
    } else {
        "minor_recombination_or_alternative_configuration"
    }
}

fn write_plot_script(paths: &PolishPaths) -> Result<(), OrgraftError> {
    let mut file = File::create(&paths.round1_plot_script)?;
    file.write_all(PLOT_SV_SUPPORT_PY.as_bytes())?;
    let mut file = File::create(&paths.round1_snv_indel_plot_script)?;
    file.write_all(PLOT_SNV_INDEL_PY.as_bytes())?;
    Ok(())
}

const PLOT_SV_SUPPORT_PY: &str = r##"#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import os
import tempfile
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "orgraft_matplotlib"))


def read_coverage(path: Path):
    rows = []
    with path.open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            rows.append(
                {
                    "position": int(row["position"]),
                    "fl_depth": float(row["fl_depth"]),
                    "partial_depth": float(row["partial_depth"]),
                    "reference_support_depth": float(row["reference_support_depth"]),
                    "total_depth": float(row["total_depth"]),
                }
            )
    return rows


def read_bubble(path: Path):
    if not path.exists():
        return []
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_group_ids(path: Path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_read_ids_file(path: Path | None):
    if path is None:
        return set()
    ids = set()
    with path.open() as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            ids.add(line.split()[0])
    return ids


def normalize_specs(spec_text: str):
    specs = []
    for item in spec_text.split(","):
        item = item.strip()
        if not item:
            continue
        if ":" not in item:
            raise SystemExit(f"invalid subgroup spec {item!r}; expected group_name:old_index")
        group_name, old_index = item.rsplit(":", 1)
        specs.append((group_name, old_index))
    return specs


def ids_for_subgroup_specs(read_group_rows, specs):
    if not specs:
        return set()
    spec_set = {(group, str(old_index)) for group, old_index in specs}
    ids = set()
    for row in read_group_rows:
        key = (row["group_name"], row["subgroup_old_index"])
        if row["read_class"] == "FL" and key in spec_set:
            ids.add(row["read_id"])
    return ids


def read_whole_read_evidence(path: Path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def cigar_ref_blocks(cigar: str, target_start: int):
    if cigar == ".":
        return []
    blocks = []
    number = ""
    ref_pos = target_start
    for ch in cigar:
        if ch.isdigit():
            number += ch
            continue
        if not number:
            continue
        length = int(number)
        number = ""
        if ch in "M=XDN":
            blocks.append((ref_pos, ref_pos + length - 1))
            ref_pos += length
        elif ch in "ISHP":
            continue
    return blocks


def highlight_depth_from_evidence(evidence_rows, highlight_ids, reference_len: int):
    if not highlight_ids:
        return []
    diff = [0] * (reference_len + 1)
    for row in evidence_rows:
        if row["read_id"] not in highlight_ids:
            continue
        if row["alignment_role"] == "secondary" or row["target_start"] == ".":
            continue
        if row["cigar"] == ".":
            blocks = [(int(row["target_start"]), int(row["target_end"]))]
        else:
            blocks = cigar_ref_blocks(row["cigar"], int(row["target_start"]))
        for start, end in blocks:
            start = max(1, min(reference_len, start))
            end = max(start, min(reference_len, end))
            diff[start - 1] += 1
            diff[end] -= 1
    depth = []
    current = 0
    for delta in diff[:-1]:
        current += delta
        depth.append(max(current, 0))
    return depth


def parse_plot_range(text: str, reference_len: int):
    text = (text or ".").strip()
    if text in {"", "."}:
        return 1, reference_len, True
    sep = "-" if "-" in text else ":"
    if sep not in text:
        raise SystemExit(f"invalid --plot-range {text!r}; expected START-END")
    start_text, end_text = text.split(sep, 1)
    start = int(start_text.replace(",", ""))
    end = int(end_text.replace(",", ""))
    if start < 1 or end < start:
        raise SystemExit(f"invalid --plot-range {text!r}; expected 1-based START-END")
    return max(1, start), min(reference_len, end), start == 1 and end >= reference_len


def range_stem(prefix: str, start: int, end: int, is_full: bool):
    return prefix if is_full else f"{prefix}_{start}_{end}"


def place_legend_above(ax, ncol=3):
    handles, labels = ax.get_legend_handles_labels()
    if handles:
        ax.legend(
            handles,
            labels,
            loc="lower center",
            bbox_to_anchor=(0.5, 1.03),
            ncol=min(ncol, len(handles)),
            frameon=False,
            fontsize=7,
            borderaxespad=0,
            columnspacing=1.2,
            handletextpad=0.4,
        )


def slice_depth(depth, start: int, end: int):
    if not depth:
        return depth
    return depth[start - 1 : end]


def save_plot(fig, plots_dir: Path, stem: str, output_format: str):
    if output_format in ("pdf", "both"):
        fig.savefig(plots_dir / f"{stem}.pdf")
    if output_format in ("png", "both"):
        fig.savefig(plots_dir / f"{stem}.png")


def plot_coverage(
    rows,
    plots_dir: Path,
    plot_start: int,
    plot_end: int,
    is_full_range: bool,
    plot_dpi: int,
    output_format: str,
    coverage_rasterize: bool,
    highlight_depth=None,
    highlight_label="highlight",
):
    if not rows:
        return
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plot_rows = [row for row in rows if plot_start <= row["position"] <= plot_end]
    x = [row["position"] for row in plot_rows]
    fl = [row["fl_depth"] for row in plot_rows]
    partial = [row["partial_depth"] for row in plot_rows]
    total = [row["total_depth"] for row in plot_rows]
    reference = [row["reference_support_depth"] for row in plot_rows]
    highlight_depth = slice_depth(highlight_depth, plot_start, plot_end)

    fig, ax = plt.subplots(figsize=(12, 3), dpi=plot_dpi)
    ax.fill_between(x, total, color="#EAB13E", alpha=1.0, label="FL+partial", rasterized=coverage_rasterize)
    ax.fill_between(x, fl, color="#D1D1D1", alpha=1.0, label="FL", rasterized=coverage_rasterize)
    ax.plot(x, reference, color="#5CAB38", linewidth=1, label="reference_support", rasterized=coverage_rasterize)
    if highlight_depth:
        ax.fill_between(
            x,
            0,
            highlight_depth,
            color="#D62728",
            alpha=0.90,
            label=highlight_label,
            zorder=20,
            rasterized=coverage_rasterize,
        )
    ax.set_xlim(plot_start, plot_end)
    max_y = max(max(total) if total else 1, max(highlight_depth) if highlight_depth else 1)
    ax.set_ylim(0, max_y + 100)
    ax.grid(True, alpha=0.5)
    place_legend_above(ax, ncol=4)
    fig.subplots_adjust(left=0.07, right=0.99, bottom=0.14, top=0.78)
    prefix = "coverage_highlight" if highlight_depth else "coverage"
    stem = range_stem(prefix, plot_start, plot_end, is_full_range)
    save_plot(fig, plots_dir, stem, output_format)
    plt.close()


def plot_bubble(rows, output_prefix: str, plots_dir: Path, plot_dpi: int = 300, output_format: str = "png"):
    if not rows:
        return
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    se1 = [int(row["se1"]) for row in rows]
    ss2 = [int(row["ss2"]) for row in rows]
    size = [max(float(row["subgroup_count_norm_per_10000"]) * 10.0, 1.0) for row in rows]
    color = [row["color"] for row in rows]
    max_coord = max(max(se1), max(ss2))

    fig, ax = plt.subplots(figsize=(10, 10), dpi=plot_dpi)
    ax.scatter(se1, ss2, s=size, c=color, alpha=0.5, linewidths=0)
    ax.set_xlim(1, max_coord)
    ax.set_ylim(1, max_coord)
    ax.grid(True, alpha=0.5)
    fig.tight_layout()
    save_plot(fig, plots_dir, output_prefix, output_format)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser()
    default_plots_dir = Path(__file__).resolve().parent
    parser.add_argument("--round-dir", type=Path, default=None, help="round directory containing 01.data and 02.plots")
    parser.add_argument("--sv-dir", type=Path, default=None, help="legacy SV directory; data is read from SV_DIR/data")
    parser.add_argument("--plots-dir", type=Path, default=default_plots_dir, help="output directory for coverage and bubble plots")
    parser.add_argument("--plot-range", default=".", help="1-based START-END interval; omit for full reference")
    parser.add_argument("--plot-dpi", type=int, default=300, help="raster DPI for PNG output and rasterized PDF artists")
    parser.add_argument("--plot-output-format", choices=("png", "pdf", "both"), default="png", help="plot output format")
    parser.add_argument("--coverage-plot-rasterize", choices=("on", "off"), default="on", help="rasterize dense coverage artists in PDF")
    parser.add_argument(
        "--sv-plot-highlight-subgroups",
        default="",
        help="Comma-separated group_name:old_index specs, e.g. type_3_subtype_rep_rep_NA:3,type_3_subtype_rep_rep_NA:5",
    )
    parser.add_argument(
        "--sv-plot-highlight-read-ids",
        "--sv-plot-highlight-read-id-file",
        type=Path,
        default=None,
        help="File with one read id per line to draw as red bottom coverage blocks",
    )
    args = parser.parse_args()
    if args.sv_dir is not None:
        data_dir = args.sv_dir / "data"
    elif args.round_dir is not None:
        data_dir = args.round_dir / "01.data"
    else:
        data_dir = args.plots_dir.parent / "01.data"
    plots_dir = args.plots_dir
    plots_dir.mkdir(parents=True, exist_ok=True)

    coverage_rows = read_coverage(data_dir / "sv_coverage.tsv")
    plot_start, plot_end, is_full_range = parse_plot_range(args.plot_range, len(coverage_rows))
    coverage_rasterize = args.coverage_plot_rasterize == "on"
    plot_coverage(coverage_rows, plots_dir, plot_start, plot_end, is_full_range, args.plot_dpi, args.plot_output_format, coverage_rasterize)
    highlight_ids = ids_for_subgroup_specs(
        read_group_ids(data_dir / "sv_read_index.tsv"),
        normalize_specs(args.sv_plot_highlight_subgroups),
    )
    highlight_ids |= read_read_ids_file(args.sv_plot_highlight_read_ids)
    if highlight_ids:
        highlight_depth = highlight_depth_from_evidence(
            read_whole_read_evidence(data_dir / "sv_read_evidence.tsv"),
            highlight_ids,
            len(coverage_rows),
        )
        plot_coverage(
            coverage_rows,
            plots_dir,
            plot_start,
            plot_end,
            is_full_range,
            args.plot_dpi,
            args.plot_output_format,
            coverage_rasterize,
            highlight_depth,
            f"highlight n={len(highlight_ids)}",
        )


if __name__ == "__main__":
    main()
"##;

const PLOT_SNV_INDEL_PY: &str = r##"#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import math
import os
import tempfile
from collections import Counter
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "orgraft_matplotlib"))


SNV_ORDER = [
    "G>A",
    "G>C",
    "G>T",
    "A>G",
    "A>C",
    "A>T",
    "C>G",
    "C>A",
    "C>T",
    "T>G",
    "T>A",
    "T>C",
]
INDEL_ORDER = [
    "poly-A",
    "poly-T",
    "poly-C",
    "poly-G",
    "tandem_size_2",
    "tandem_size_3",
    "tandem_size_4",
    "tandem_size_5",
    "tandem_size_6",
    "tandem_size_other",
    "InDel,MMEJ",
    "InDel,NHEJ",
]


def read_points(path: Path):
    with path.open() as handle:
        rows = []
        for row in csv.DictReader(handle, delimiter="\t"):
            row["pos"] = int(row["pos"])
            row["total_count"] = int(row["total_count"])
            row["depth"] = int(row["depth"])
            row["frequency"] = float(row["frequency"])
            row["ref_count"] = int(row["ref_count"])
            row["ref_frequency"] = float(row["ref_frequency"])
            row["max_alt_count"] = int(row["max_alt_count"])
            row["max_alt_frequency"] = float(row["max_alt_frequency"])
            row["high_count"] = int(row["high_count"])
            row["non_high_count"] = int(row["non_high_count"])
            rows.append(row)
    return rows


def read_reference_len(data_dir: Path, rows):
    coverage = data_dir / "sv_coverage.tsv"
    if coverage.exists():
        last = 0
        with coverage.open() as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                last = int(row["position"])
        if last > 0:
            return last
    return max((row["pos"] for row in rows), default=0)


def parse_confidence_counts(text: str):
    counts = {}
    for item in text.split(";"):
        if "=" not in item:
            continue
        key, value = item.split("=", 1)
        try:
            counts[key] = int(value)
        except ValueError:
            counts[key] = 0
    return counts


def confidence_low(row, spec: str):
    spec = spec.strip().lower()
    if spec in {"", "none", "off", "false"}:
        return False
    counts = parse_confidence_counts(row["confidence_counts"])
    if spec in {"non-high", "non_high", "not-high", "not_high"}:
        return sum(value for key, value in counts.items() if key != "high") > 0
    selected = {item.strip() for item in spec.split(",") if item.strip()}
    return any(counts.get(item, 0) > 0 for item in selected)


def is_low_confidence(row, args):
    if confidence_low(row, args.low_confidence_types):
        return True
    if args.low_min_reads > 0 and row["total_count"] < args.low_min_reads:
        return True
    if args.low_min_fraction > 0 and row["frequency"] < args.low_min_fraction:
        return True
    return False


def risk_class(row, threshold: float):
    if row["frequency"] < threshold:
        return "none"
    if row["ref_count"] >= row["max_alt_count"]:
        return "orange"
    return "red"


def split_by_kind(rows, kind: str):
    return [row for row in rows if row["variant_kind"] == kind]


def parse_plot_range(text: str, reference_len: int):
    text = (text or ".").strip()
    if text in {"", "."}:
        return 1, reference_len, True
    sep = "-" if "-" in text else ":"
    if sep not in text:
        raise SystemExit(f"invalid --plot-range {text!r}; expected START-END")
    start_text, end_text = text.split(sep, 1)
    start = int(start_text.replace(",", ""))
    end = int(end_text.replace(",", ""))
    if start < 1 or end < start:
        raise SystemExit(f"invalid --plot-range {text!r}; expected 1-based START-END")
    return max(1, start), min(reference_len, end), start == 1 and end >= reference_len


def range_stem(prefix: str, start: int, end: int, is_full: bool):
    return prefix if is_full else f"{prefix}_{start}_{end}"


def place_legend_above(ax, ncol=4):
    handles, labels = ax.get_legend_handles_labels()
    if handles:
        ax.legend(
            handles,
            labels,
            loc="lower center",
            bbox_to_anchor=(0.5, 1.03),
            ncol=min(ncol, len(handles)),
            frameon=False,
            fontsize=7,
            borderaxespad=0,
            columnspacing=1.2,
            handletextpad=0.4,
        )


def y_limit(rows, threshold: float):
    max_freq = max([row["frequency"] for row in rows] + [threshold, 0.05])
    top = max(0.5, max_freq * 1.08)
    return min(1.0, top)


def save_plot(fig, plots_dir: Path, stem: str, output_format: str):
    if output_format in ("pdf", "both"):
        fig.savefig(plots_dir / f"{stem}.pdf")
    if output_format in ("png", "both"):
        fig.savefig(plots_dir / f"{stem}.png")


def plot_frequency(
    rows,
    kind: str,
    reference_len: int,
    plots_dir: Path,
    args,
    plot_start: int,
    plot_end: int,
    is_full_range: bool,
):
    kind_rows = [
        row for row in split_by_kind(rows, kind) if plot_start <= row["pos"] <= plot_end
    ]
    if not kind_rows:
        return
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    marker = "o" if kind == "SNV" else "^"
    low_rows = [row for row in kind_rows if is_low_confidence(row, args)]
    main_rows = [row for row in kind_rows if not is_low_confidence(row, args)]
    rasterize_points = args.snv_indel_plot_rasterize == "on"

    fig, ax = plt.subplots(figsize=(12, 3), dpi=args.plot_dpi)
    if main_rows:
        ax.scatter(
            [row["pos"] / 1000.0 for row in main_rows],
            [row["frequency"] for row in main_rows],
            s=8,
            marker=marker,
            c="#10AFCF",
            alpha=0.58,
            linewidths=0,
            rasterized=rasterize_points,
            label="high confidence",
            zorder=10,
        )
    if low_rows:
        ax.scatter(
            [row["pos"] / 1000.0 for row in low_rows],
            [row["frequency"] for row in low_rows],
            s=8,
            marker=marker,
            c="#CFCFCF",
            alpha=0.48,
            linewidths=0,
            rasterized=rasterize_points,
            label="low confidence",
            zorder=8,
        )
    orange_rows = [
        row for row in kind_rows if risk_class(row, args.high_risk_fraction) == "orange"
    ]
    red_rows = [row for row in kind_rows if risk_class(row, args.high_risk_fraction) == "red"]
    if orange_rows:
        ax.scatter(
            [row["pos"] / 1000.0 for row in orange_rows],
            [row["frequency"] for row in orange_rows],
            s=15,
            marker=marker,
            c="#F39C12",
            alpha=0.90,
            edgecolors="none",
            rasterized=rasterize_points,
            label="risk: ref still top",
            zorder=30,
        )
    if red_rows:
        ax.scatter(
            [row["pos"] / 1000.0 for row in red_rows],
            [row["frequency"] for row in red_rows],
            s=16,
            marker=marker,
            c="#E41A1C",
            alpha=0.92,
            edgecolors="none",
            rasterized=rasterize_points,
            label="risk: ref not top",
            zorder=35,
        )
    ax.axhline(0.5, color="#E41A1C", linestyle="--", linewidth=1.0, alpha=0.45, zorder=5)
    ax.set_xlim(plot_start / 1000.0, max(plot_end / 1000.0, plot_start / 1000.0 + 0.001))
    ax.set_ylim(0, y_limit(kind_rows, args.high_risk_fraction))
    ax.set_ylabel("Frequency")
    ax.set_xlabel("Position (kb)")
    ax.grid(True, alpha=0.45)
    place_legend_above(ax, ncol=4)
    fig.subplots_adjust(left=0.07, right=0.99, bottom=0.18, top=0.78)
    stem = f"{kind.lower()}_frequency"
    stem = range_stem(stem, plot_start, plot_end, is_full_range)
    save_plot(fig, plots_dir, stem, args.plot_output_format)
    plt.close(fig)


def type_labels(row):
    labels = [item.strip() for item in row["plot_type"].split(";") if item.strip()]
    return labels or [row["type"]]


def infer_organelle(plots_dir: Path):
    try:
        return plots_dir.resolve().parents[3].name
    except IndexError:
        return "current"


def plot_type_counts(rows, kind: str, plots_dir: Path, args):
    kind_rows = [row for row in split_by_kind(rows, kind) if not is_low_confidence(row, args)]
    if not kind_rows:
        return
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    counts = Counter()
    for row in kind_rows:
        for label in type_labels(row):
            counts[label] += 1
    base_order = SNV_ORDER if kind == "SNV" else INDEL_ORDER
    extras = sorted(label for label in counts if label not in base_order)
    labels = [label for label in base_order if counts.get(label, 0) > 0] + extras
    if not labels:
        return
    organelle = infer_organelle(plots_dir)
    color = "#FFB52E" if organelle == "mito" else "#05B83F"

    figure_width = max(5.0, min(12.0, len(labels) * 0.78))
    bar_width = 0.38 if len(labels) <= 6 else 0.46
    fig, ax = plt.subplots(figsize=(figure_width, 3), dpi=args.plot_dpi)
    xs = list(range(len(labels)))
    ax.bar(xs, [counts[label] for label in labels], color=color, width=bar_width, label=organelle)
    ax.set_xticks(xs)
    ax.set_xticklabels(labels, rotation=45, ha="right")
    ax.set_ylabel("SNV counts" if kind == "SNV" else "InDel counts")
    ax.grid(axis="y", alpha=0.5)
    place_legend_above(ax, ncol=1)
    fig.subplots_adjust(left=0.12, right=0.98, bottom=0.35, top=0.78)
    stem = "snv_type" if kind == "SNV" else "indel_type"
    save_plot(fig, plots_dir, stem, args.plot_output_format)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser()
    default_plots_dir = Path(__file__).resolve().parent
    parser.add_argument("--round-dir", type=Path, default=None, help="round directory containing 01.data and 02.plots")
    parser.add_argument("--data-dir", type=Path, default=None, help="directory containing snv_indel_plot_points.tsv")
    parser.add_argument("--plots-dir", type=Path, default=default_plots_dir, help="output directory for SNV/InDel plots")
    parser.add_argument("--points", type=Path, default=None, help="explicit snv_indel_plot_points.tsv path")
    parser.add_argument("--plot-range", default=".", help="1-based START-END interval; omit for full reference")
    parser.add_argument("--plot-dpi", type=int, default=300, help="raster DPI for PNG output and rasterized PDF artists")
    parser.add_argument("--plot-output-format", choices=("png", "pdf", "both"), default="png", help="plot output format")
    parser.add_argument("--snv-indel-plot-rasterize", choices=("on", "off"), default="on", help="rasterize dense SNV/InDel scatter artists in PDF")
    parser.add_argument("--low-confidence-types", default="non-high", help="grey point classes: non-high, none, or comma-separated confidence labels")
    parser.add_argument("--low-min-reads", type=int, default=3, help="grey points below this read count")
    parser.add_argument("--low-min-fraction", type=float, default=0.0, help="grey points below this variant frequency")
    parser.add_argument("--high-risk-fraction", type=float, default=0.5, help="orange/red threshold for high-risk points")
    args = parser.parse_args()

    if args.data_dir is not None:
        data_dir = args.data_dir
    elif args.round_dir is not None:
        data_dir = args.round_dir / "01.data"
    else:
        data_dir = args.plots_dir.parent / "01.data"
    plots_dir = args.plots_dir
    plots_dir.mkdir(parents=True, exist_ok=True)

    points_path = args.points or data_dir / "snv_indel_plot_points.tsv"
    if not points_path.exists():
        raise SystemExit(f"missing SNV/InDel plot points: {points_path}")
    rows = read_points(points_path)
    reference_len = read_reference_len(data_dir, rows)
    if reference_len <= 0 or not rows:
        return
    plot_start, plot_end, is_full_range = parse_plot_range(args.plot_range, reference_len)
    for kind in ("SNV", "InDel"):
        plot_frequency(
            rows,
            kind,
            reference_len,
            plots_dir,
            args,
            plot_start=plot_start,
            plot_end=plot_end,
            is_full_range=is_full_range,
        )
        plot_type_counts(rows, kind, plots_dir, args)


if __name__ == "__main__":
    main()
"##;

fn sorted_read_groups<'a>(
    groups: &'a BTreeMap<String, ReadGroupStats>,
) -> Vec<(&'a String, &'a ReadGroupStats)> {
    let mut sorted = groups.iter().collect::<Vec<_>>();
    sorted.sort_by(|left, right| {
        left.1
            .num_align
            .cmp(&right.1.num_align)
            .then_with(|| left.1.subtype.cmp(&right.1.subtype))
            .then_with(|| left.0.cmp(right.0))
    });
    sorted
}

fn write_read_group_summary(
    path: &Path,
    sorted_groups: &[(&String, &ReadGroupStats)],
) -> Result<(), OrgraftError> {
    let mut file = File::create(path)?;
    writeln!(
        file,
        "group_name\tnum_align\tsubtype\ttotal_reads\tfl_reads\tpartial_reads\tfl_multi_reads\tfl_fraction\tpartial_fraction\tpriority"
    )?;
    for (group_name, stats) in sorted_groups {
        let total = stats.total_reads.max(1) as f64;
        writeln!(
            file,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.6}\t{:.6}\t{}",
            group_name,
            stats.num_align,
            stats.subtype,
            stats.total_reads,
            stats.fl_reads,
            stats.partial_reads,
            stats.fl_multi_reads,
            stats.fl_reads as f64 / total,
            stats.partial_reads as f64 / total,
            read_group_priority(stats),
        )?;
    }
    Ok(())
}

fn write_read_subgroup_summary(
    path: &Path,
    records: &[AlignmentSummaryRecord],
    groups: &BTreeMap<String, ReadGroupStats>,
    subgroups: &BTreeMap<String, BTreeMap<SubgroupKey, Vec<usize>>>,
    subgroup_old_index: &BTreeMap<String, BTreeMap<SubgroupKey, usize>>,
) -> Result<usize, OrgraftError> {
    let mut file = File::create(path)?;
    writeln!(
        file,
        "group_name\tnum_align\tsubtype\told_index\tboundary_key\tstrand_str\tgroup_fl_reads\tsubgroup_reads\tsubgroup_multi_reads\tcount_rank\tpriority\tmin_olps\tmid_olps\tmax_olps"
    )?;
    let mut written = 0usize;
    for (group_name, group_subgroups) in subgroups {
        let Some(stats) = groups.get(group_name) else {
            continue;
        };
        let mut rows = group_subgroups
            .iter()
            .map(|(key, record_indices)| {
                let old_index = subgroup_old_index
                    .get(group_name)
                    .and_then(|indices| indices.get(key))
                    .copied()
                    .unwrap_or(0);
                (key, record_indices, old_index)
            })
            .collect::<Vec<_>>();
        let mut count_ranked = rows.clone();
        count_ranked.sort_by(|left, right| {
            right
                .1
                .len()
                .cmp(&left.1.len())
                .then_with(|| left.2.cmp(&right.2))
        });
        let count_rank_by_old_index = count_ranked
            .iter()
            .enumerate()
            .map(|(index, (_, _, old_index))| (*old_index, index + 1))
            .collect::<BTreeMap<_, _>>();

        rows.sort_by(|left, right| left.2.cmp(&right.2));
        for (key, record_indices, old_index) in rows {
            let subgroup_multi_reads = record_indices
                .iter()
                .filter(|index| records[**index].is_multi())
                .count();
            let (min_olps, mid_olps, max_olps) =
                subgroup_overlap_summaries(records, record_indices, stats.num_align);
            let priority = read_subgroup_priority(stats, record_indices.len());
            writeln!(
                file,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                group_name,
                stats.num_align,
                stats.subtype,
                old_index,
                key.label(),
                subgroup_strand_summary(records, record_indices, stats.num_align),
                stats.fl_reads,
                record_indices.len(),
                subgroup_multi_reads,
                count_rank_by_old_index
                    .get(&old_index)
                    .copied()
                    .unwrap_or(0),
                priority,
                min_olps,
                mid_olps,
                max_olps,
            )?;
            written += 1;
        }
    }
    Ok(written)
}

fn write_read_group_ids(
    path: &Path,
    records: &[AlignmentSummaryRecord],
    subgroup_old_index: &BTreeMap<String, BTreeMap<SubgroupKey, usize>>,
) -> Result<(), OrgraftError> {
    let mut file = File::create(path)?;
    writeln!(
        file,
        "read_id\ttarget_id\tread_length\tread_class\tgroup_name\tnum_align\tsubtype\tpercent_total\tis_multi\tsubgroup_old_index\tsubgroup_key"
    )?;
    for record in records {
        let group_name = record.group_name();
        let (subgroup_index, subgroup_key) = record
            .subgroup_key()
            .map(|key| {
                let index = subgroup_old_index
                    .get(&group_name)
                    .and_then(|indices| indices.get(&key))
                    .copied()
                    .unwrap_or(0);
                (index.to_string(), key.label())
            })
            .unwrap_or_else(|| (".".to_string(), ".".to_string()));
        writeln!(
            file,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            record.read_id,
            record.target_id,
            record.read_len,
            if record.is_fl() { "FL" } else { "partial" },
            group_name,
            record.num_align,
            record.subtype,
            format_py_float(record.percent_total),
            record.is_multi(),
            subgroup_index,
            subgroup_key,
        )?;
    }
    Ok(())
}

fn subgroup_overlap_summaries(
    records: &[AlignmentSummaryRecord],
    record_indices: &[usize],
    num_align: usize,
) -> (String, String, String) {
    let mut min_values = Vec::new();
    let mut mid_values = Vec::new();
    let mut max_values = Vec::new();
    for align_index in 0..num_align.saturating_sub(1) {
        let values = record_indices
            .iter()
            .filter_map(|record_index| records[*record_index].alignments[align_index].olp)
            .collect::<Vec<_>>();
        if values.is_empty() {
            min_values.push(".".to_string());
            mid_values.push(".".to_string());
            max_values.push(".".to_string());
            continue;
        }
        let min = values.iter().min().copied().unwrap_or(0);
        let max = values.iter().max().copied().unwrap_or(0);
        min_values.push(min.to_string());
        mid_values.push(format_py_float((min as f64 + max as f64) / 2.0));
        max_values.push(max.to_string());
    }
    (
        min_values.join(","),
        mid_values.join(","),
        max_values.join(","),
    )
}

fn subgroup_strand_summary(
    records: &[AlignmentSummaryRecord],
    record_indices: &[usize],
    num_align: usize,
) -> String {
    (0..num_align)
        .map(|align_index| {
            let mut counts = BTreeMap::new();
            for record_index in record_indices {
                let strand = records[*record_index].alignments[align_index].strand;
                *counts.entry(strand).or_insert(0usize) += 1;
            }
            counts
                .into_iter()
                .max_by(|left, right| left.1.cmp(&right.1).then_with(|| right.0.cmp(&left.0)))
                .map(|(strand, _)| strand.to_string())
                .unwrap_or_else(|| ".".to_string())
        })
        .collect::<Vec<_>>()
        .join(",")
}

fn subgroup_old_index_order(left: &SubgroupKey, right: &SubgroupKey) -> Ordering {
    left.0
        .first()
        .map(|pair| pair.0)
        .cmp(&right.0.first().map(|pair| pair.0))
        .then_with(|| left.cmp(right))
}

fn read_group_priority(stats: &ReadGroupStats) -> &'static str {
    if is_core_read_group(stats.num_align, &stats.subtype) {
        "core"
    } else if stats.num_align >= 3 && stats.fl_reads >= READ_GROUP_WATCH_FL_MIN {
        "watch"
    } else {
        "other"
    }
}

fn read_subgroup_priority(stats: &ReadGroupStats, subgroup_reads: usize) -> &'static str {
    if is_core_read_group(stats.num_align, &stats.subtype) {
        "core"
    } else if subgroup_reads >= READ_SUBGROUP_WATCH_MIN {
        "watch"
    } else {
        "other"
    }
}

fn is_core_read_group(num_align: usize, subtype: &str) -> bool {
    num_align == 1 && subtype == "NA"
        || num_align == 2 && matches!(subtype, "ins_NA" | "ref_NA" | "rep_NA")
}

#[cfg(test)]
fn parse_alignment_summary_record(
    line: &str,
    line_number: usize,
) -> Result<AlignmentSummaryRecord, OrgraftError> {
    let fields = line.split('\t').collect::<Vec<_>>();
    if fields.len() < 5 {
        return Err(OrgraftError::InvalidArgument(format!(
            "alignment_summary:{line_number} has fewer than 5 columns"
        )));
    }
    let (num_align, subtype) = parse_alignment_type(fields[3], line_number)?;
    let alignments = fields
        .iter()
        .skip(5)
        .map(|field| parse_summary_alignment(field, line_number))
        .collect::<Result<Vec<_>, _>>()?;
    if alignments.len() != num_align {
        return Err(OrgraftError::InvalidArgument(format!(
            "alignment_summary:{line_number} aln_type={num_align} but has {} alignment columns",
            alignments.len()
        )));
    }
    Ok(AlignmentSummaryRecord {
        read_id: fields[0].to_string(),
        target_id: fields[1].to_string(),
        read_len: parse_usize_value(fields[2], "alignment summary read length")?,
        num_align,
        subtype,
        percent_total: parse_f64_label(fields[4], "alignment summary percent_total")?,
        alignments,
    })
}

#[cfg(test)]
fn parse_alignment_type(value: &str, line_number: usize) -> Result<(usize, String), OrgraftError> {
    let Some(rest) = value.strip_prefix("aln_type=") else {
        return Err(OrgraftError::InvalidArgument(format!(
            "alignment_summary:{line_number} invalid aln_type field `{value}`"
        )));
    };
    let Some((num_align, subtype)) = rest.split_once(';') else {
        return Err(OrgraftError::InvalidArgument(format!(
            "alignment_summary:{line_number} invalid aln_type field `{value}`"
        )));
    };
    Ok((
        parse_usize_value(num_align, "alignment summary aln_type")?,
        subtype.replace(',', "_"),
    ))
}

#[cfg(test)]
fn parse_summary_alignment(
    value: &str,
    line_number: usize,
) -> Result<SummaryAlignment, OrgraftError> {
    let mut olp = None;
    let mut strand = None;
    let mut qs = None;
    let mut qe = None;
    let mut ss = None;
    let mut se = None;
    let mut cn = None;
    for part in value.split(';') {
        let Some((key, value)) = part.split_once('=') else {
            continue;
        };
        match key {
            "olp" => {
                if value != "NA" {
                    olp = Some(parse_isize_label(value, "alignment summary olp")?);
                }
            }
            "strand" => strand = value.chars().next(),
            "qs" => qs = Some(parse_usize_value(value, "alignment summary qs")?),
            "qe" => qe = Some(parse_usize_value(value, "alignment summary qe")?),
            "ss" => ss = Some(parse_isize_label(value, "alignment summary ss")?),
            "se" => se = Some(parse_isize_label(value, "alignment summary se")?),
            "cn" => cn = Some(parse_usize_value(value, "alignment summary cn")?),
            _ => {}
        }
    }
    Ok(SummaryAlignment {
        olp,
        strand: strand.ok_or_else(|| {
            OrgraftError::InvalidArgument(format!(
                "alignment_summary:{line_number} missing strand in `{value}`"
            ))
        })?,
        qs: qs.ok_or_else(|| {
            OrgraftError::InvalidArgument(format!(
                "alignment_summary:{line_number} missing qs in `{value}`"
            ))
        })?,
        qe: qe.ok_or_else(|| {
            OrgraftError::InvalidArgument(format!(
                "alignment_summary:{line_number} missing qe in `{value}`"
            ))
        })?,
        ss: ss.ok_or_else(|| {
            OrgraftError::InvalidArgument(format!(
                "alignment_summary:{line_number} missing ss in `{value}`"
            ))
        })?,
        se: se.ok_or_else(|| {
            OrgraftError::InvalidArgument(format!(
                "alignment_summary:{line_number} missing se in `{value}`"
            ))
        })?,
        cn: cn.unwrap_or(1),
    })
}

fn sorted_noncontained_alignments(alignments: &[BlastLikeAlignment]) -> Vec<BlastLikeAlignment> {
    let mut keep = vec![true; alignments.len()];
    loop {
        let mut removed = 0usize;
        for i in 0..alignments.len() {
            if !keep[i] {
                continue;
            }
            let i_start = alignments[i].query_start;
            let i_end = alignments[i].query_end;
            for j in 0..alignments.len() {
                if i == j || !keep[j] {
                    continue;
                }
                let j_start = alignments[j].query_start;
                let j_end = alignments[j].query_end;
                if i_start <= j_start && j_end <= i_end {
                    keep[j] = false;
                    removed += 1;
                } else if j_start <= i_start && i_end <= j_end {
                    keep[i] = false;
                    removed += 1;
                    break;
                }
            }
        }
        if removed == 0 {
            break;
        }
    }
    let mut kept = alignments
        .iter()
        .zip(keep)
        .filter_map(|(alignment, keep)| keep.then(|| alignment.clone()))
        .collect::<Vec<_>>();
    kept.sort_by_key(|alignment| alignment.query_start);
    kept
}

fn best_terminal_microindel_extension(query_seq: &str, target_seq: &str) -> (usize, usize, usize) {
    let query_seq = query_seq
        .as_bytes()
        .iter()
        .copied()
        .map(|base| base.to_ascii_uppercase())
        .collect::<Vec<_>>();
    let target_seq = target_seq
        .as_bytes()
        .iter()
        .copied()
        .map(|base| base.to_ascii_uppercase())
        .collect::<Vec<_>>();
    let mut best_query_len = 0usize;
    let mut best_target_len = 0usize;
    let mut best_key = (-1isize, -1isize, -1isize);
    for query_len in 3..=query_seq.len().min(TERMINAL_EXTENSION_WINDOW) {
        for target_len in 3..=target_seq.len().min(TERMINAL_EXTENSION_WINDOW) {
            let gap_len = query_len.abs_diff(target_len);
            if gap_len == 0 || gap_len > TERMINAL_EXTENSION_MAX_GAP {
                continue;
            }
            if !can_match_after_skipping(
                &query_seq[..query_len],
                &target_seq[..target_len],
                TERMINAL_EXTENSION_MAX_GAP,
            ) {
                continue;
            }
            let matches = query_len.min(target_len);
            let key = (
                matches as isize,
                -((query_len + target_len) as isize),
                -(gap_len as isize),
            );
            if key > best_key {
                best_key = key;
                best_query_len = query_len;
                best_target_len = target_len;
            }
        }
    }
    (
        best_query_len,
        best_target_len,
        best_query_len.min(best_target_len),
    )
}

fn can_match_after_skipping(sequence_a: &[u8], sequence_b: &[u8], max_skips: usize) -> bool {
    if sequence_a.len().abs_diff(sequence_b.len()) > max_skips {
        return false;
    }
    if sequence_a
        .iter()
        .chain(sequence_b.iter())
        .any(|base| !matches!(base, b'A' | b'C' | b'G' | b'T'))
    {
        return false;
    }
    if sequence_a.len() == sequence_b.len() {
        return sequence_a == sequence_b;
    }
    if sequence_a.len() > sequence_b.len() {
        return can_drop_bases_to_match(sequence_a, sequence_b, max_skips, 0, 0, 0);
    }
    can_drop_bases_to_match(sequence_b, sequence_a, max_skips, 0, 0, 0)
}

fn can_drop_bases_to_match(
    longer: &[u8],
    shorter: &[u8],
    max_skips: usize,
    long_index: usize,
    short_index: usize,
    skipped: usize,
) -> bool {
    if short_index == shorter.len() {
        return skipped + longer.len().saturating_sub(long_index) <= max_skips;
    }
    if long_index == longer.len() || skipped > max_skips {
        return false;
    }
    if longer[long_index] == shorter[short_index]
        && can_drop_bases_to_match(
            longer,
            shorter,
            max_skips,
            long_index + 1,
            short_index + 1,
            skipped,
        )
    {
        return true;
    }
    can_drop_bases_to_match(
        longer,
        shorter,
        max_skips,
        long_index + 1,
        short_index,
        skipped + 1,
    )
}

fn cleanup_input_sidecars(paths: &PolishPaths) -> Result<(), OrgraftError> {
    for source in [&paths.input_draft, &paths.input_reference] {
        let fai_path = PathBuf::from(format!("{}.fai", source.display()));
        if fai_path.exists() {
            fs::remove_file(fai_path)?;
        }
    }
    Ok(())
}

#[derive(Debug, Clone)]
struct PileupRoundReport {
    round: usize,
    input_length: usize,
    output_length: usize,
    alignments_used: usize,
    substitutions: usize,
    deletions: usize,
    inserted_bases: usize,
    low_coverage_bases: usize,
}

#[derive(Debug, Clone)]
struct PileupState {
    bases: Vec<[u32; 4]>,
    deletions_by_pos: Vec<u32>,
    coverage: Vec<u32>,
    insertions: Vec<HashMap<String, u32>>,
    substitutions: usize,
    deletions: usize,
    inserted_bases: usize,
    low_coverage_bases: usize,
}

impl PileupState {
    fn new(reference_len: usize) -> Self {
        Self {
            bases: vec![[0; 4]; reference_len],
            deletions_by_pos: vec![0; reference_len],
            coverage: vec![0; reference_len],
            insertions: vec![HashMap::new(); reference_len + 1],
            substitutions: 0,
            deletions: 0,
            inserted_bases: 0,
            low_coverage_bases: 0,
        }
    }

    fn load_sam_reader<R: BufRead>(&mut self, reader: R) -> Result<usize, OrgraftError> {
        let mut used = 0usize;
        for line in reader.lines() {
            let line = line?;
            if line.starts_with('@') || line.trim().is_empty() {
                continue;
            }
            if self.add_sam_alignment(&line)? {
                used += 1;
            }
        }
        Ok(used)
    }

    fn add_sam_alignment(&mut self, line: &str) -> Result<bool, OrgraftError> {
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 11 {
            return Ok(false);
        }
        let flag = parse_u16(fields[1], "SAM flag")?;
        if flag & 0x4 != 0 || flag & 0x100 != 0 || flag & 0x800 != 0 {
            return Ok(false);
        }
        let mapq = parse_u8(fields[4], "SAM MAPQ")?;
        if mapq < 20 {
            return Ok(false);
        }
        let pos = parse_usize_value(fields[3], "SAM POS")?;
        if pos == 0 {
            return Ok(false);
        }
        let cigar = fields[5];
        let read = fields[9].as_bytes();
        if cigar == "*" || read == b"*" {
            return Ok(false);
        }

        let mut ref_pos = pos - 1;
        let mut read_pos = 0usize;
        for (len, op) in parse_cigar(cigar)? {
            match op {
                'M' | '=' | 'X' => {
                    for offset in 0..len {
                        let Some(base) = read.get(read_pos + offset) else {
                            break;
                        };
                        let current_ref = ref_pos + offset;
                        if current_ref >= self.coverage.len() {
                            break;
                        }
                        self.coverage[current_ref] += 1;
                        if let Some(index) = base_index(*base) {
                            self.bases[current_ref][index] += 1;
                        }
                    }
                    read_pos += len;
                    ref_pos += len;
                }
                'I' => {
                    if ref_pos <= self.coverage.len() && read_pos + len <= read.len() {
                        let inserted = String::from_utf8_lossy(&read[read_pos..read_pos + len])
                            .to_ascii_uppercase();
                        if inserted.bytes().all(|base| base_index(base).is_some()) {
                            *self.insertions[ref_pos].entry(inserted).or_insert(0) += 1;
                        }
                    }
                    read_pos += len;
                }
                'D' | 'N' => {
                    for offset in 0..len {
                        let current_ref = ref_pos + offset;
                        if current_ref >= self.coverage.len() {
                            break;
                        }
                        self.coverage[current_ref] += 1;
                        self.deletions_by_pos[current_ref] += 1;
                    }
                    ref_pos += len;
                }
                'S' => read_pos += len,
                'H' | 'P' => {}
                other => {
                    return Err(OrgraftError::InvalidArgument(format!(
                        "unsupported CIGAR op `{other}` in {cigar}"
                    )));
                }
            }
        }
        Ok(true)
    }

    fn consensus(&mut self, reference: &str) -> String {
        const MIN_COVERAGE: u32 = 3;
        const BASE_FRACTION: f64 = 0.60;
        const DELETE_FRACTION: f64 = 0.60;
        const INSERT_FRACTION: f64 = 0.55;

        let mut out = String::with_capacity(reference.len());
        for (index, reference_base) in reference.bytes().enumerate() {
            self.push_consensus_insertion(index, INSERT_FRACTION, &mut out);
            let coverage = self.coverage[index];
            if coverage < MIN_COVERAGE {
                self.low_coverage_bases += 1;
                out.push(reference_base as char);
                continue;
            }

            let delete_count = self.deletions_by_pos[index];
            let (best_base_index, best_base_count) = best_base_count(self.bases[index]);
            if (delete_count as f64) / (coverage as f64) >= DELETE_FRACTION
                && delete_count > best_base_count
            {
                self.deletions += 1;
                continue;
            }

            if best_base_count > 0 && (best_base_count as f64) / (coverage as f64) >= BASE_FRACTION
            {
                let best_base = index_base(best_base_index);
                if reference_base.to_ascii_uppercase() != best_base {
                    self.substitutions += 1;
                }
                out.push(best_base as char);
            } else {
                out.push(reference_base as char);
            }
        }
        self.push_consensus_insertion(reference.len(), INSERT_FRACTION, &mut out);
        out
    }

    fn push_consensus_insertion(&mut self, index: usize, fraction: f64, out: &mut String) {
        let Some((inserted, support)) = self.insertions[index]
            .iter()
            .max_by_key(|(_, support)| *support)
        else {
            return;
        };
        let denominator = if index < self.coverage.len() {
            self.coverage[index].max(1)
        } else {
            self.coverage.last().copied().unwrap_or(1).max(1)
        };
        if *support >= 3 && (*support as f64) / (denominator as f64) >= fraction {
            self.inserted_bases += inserted.len();
            out.push_str(inserted);
        }
    }
}

#[derive(Debug, Clone)]
struct FastaRecordInfo {
    id: String,
    header: String,
}

fn write_report(
    path: &Path,
    options: &PolishOptions,
    inputs: &ResolvedInputs,
    paths: &PolishPaths,
    extracted_draft: &FastaRecordInfo,
    extracted_reference: &FastaRecordInfo,
    stage_records: &[StageRecord],
    round_reports: &[PileupRoundReport],
    alignment_report: Option<&AlignmentReport>,
    sv_eval_report: Option<&SvEvalReport>,
    snv_indel_report: Option<&SnvIndelReport>,
    command_records: &[CommandRecord],
) -> Result<(), OrgraftError> {
    let mut file = File::create(path)?;
    writeln!(file, "section\tname\tmetric\tvalue")?;
    write_report_row(
        &mut file,
        "global",
        "run",
        "orgraft_version",
        env!("CARGO_PKG_VERSION"),
    )?;
    write_report_row(&mut file, "global", "run", "organelle", &options.organelle)?;
    write_report_row(&mut file, "global", "run", "subgraph", &options.subgraph)?;
    write_report_row(
        &mut file,
        "global",
        "run",
        "threads",
        &options.threads.to_string(),
    )?;
    write_report_row(
        &mut file,
        "global",
        "run",
        "max_rounds",
        &options.max_rounds.to_string(),
    )?;
    write_report_row(
        &mut file,
        "global",
        "run",
        "validate_round",
        &options.validate_round.to_string(),
    )?;
    write_report_row(&mut file, "global", "run", "created_at", &timestamp())?;
    write_report_row(
        &mut file,
        "global",
        "run",
        "per_read_variant_calls",
        &options.per_read_variant_calls.to_string(),
    )?;
    write_report_row(
        &mut file,
        "global",
        "run",
        "snv_indel_overlap_policy",
        options.snv_indel_overlap_policy.as_str(),
    )?;
    write_report_row(&mut file, "global", "run", "plot_attempted", "true")?;
    if let Some(range) = &options.plot_range {
        write_report_row(&mut file, "global", "run", "plot_range", &range.as_arg())?;
    }
    write_report_row(
        &mut file,
        "global",
        "run",
        "plot_dpi",
        &options.plot_dpi.to_string(),
    )?;
    write_report_row(
        &mut file,
        "global",
        "run",
        "plot_output_format",
        options.plot_output_format.as_str(),
    )?;
    write_report_row(
        &mut file,
        "global",
        "run",
        "coverage_plot_rasterize",
        &options.coverage_plot_rasterize.to_string(),
    )?;
    write_report_row(
        &mut file,
        "global",
        "run",
        "snv_indel_plot_rasterize",
        &options.snv_indel_plot_rasterize.to_string(),
    )?;
    write_report_row(
        &mut file,
        "global",
        "run",
        "auto_sv_plot_highlight_min_fraction",
        &format!("{:.6}", options.sv_plot_highlight_min_fraction),
    )?;
    write_report_row(
        &mut file,
        "global",
        "run",
        "auto_sv_plot_highlight_min_reads",
        &options.sv_plot_highlight_min_reads.to_string(),
    )?;
    write_report_row(
        &mut file,
        "global",
        "run",
        "snv_indel_plot_low_confidence",
        &options.snv_indel_plot_low_confidence,
    )?;
    write_report_row(
        &mut file,
        "global",
        "run",
        "snv_indel_plot_low_min_reads",
        &options.snv_indel_plot_low_min_reads.to_string(),
    )?;
    write_report_row(
        &mut file,
        "global",
        "run",
        "snv_indel_plot_low_min_fraction",
        &format!("{:.6}", options.snv_indel_plot_low_min_fraction),
    )?;
    write_report_row(
        &mut file,
        "global",
        "run",
        "snv_indel_plot_high_risk_fraction",
        &format!("{:.6}", options.snv_indel_plot_high_risk_fraction),
    )?;
    if !options.sv_plot_highlight_subgroups.is_empty() {
        write_report_row(
            &mut file,
            "global",
            "run",
            "manual_highlight_subgroups",
            &options.sv_plot_highlight_subgroups.join(","),
        )?;
    }
    if let Some(path) = &options.sv_plot_highlight_read_ids {
        write_report_row(
            &mut file,
            "global",
            "run",
            "manual_highlight_read_ids",
            &display_path(path),
        )?;
    }
    write_report_row(
        &mut file,
        "input",
        "draft",
        "path",
        &display_path(&inputs.draft),
    )?;
    write_report_row(
        &mut file,
        "input",
        "draft",
        "record_id",
        &extracted_draft.id,
    )?;
    write_report_row(
        &mut file,
        "input",
        "reference",
        "path",
        &display_path(&inputs.reference),
    )?;
    write_report_row(
        &mut file,
        "input",
        "reference",
        "record_id",
        &extracted_reference.id,
    )?;
    write_report_row(
        &mut file,
        "input",
        "reads",
        "path",
        &display_path(&inputs.reads),
    )?;
    write_report_row(
        &mut file,
        "input",
        "soft_paths",
        "path",
        &display_path(&inputs.soft_paths),
    )?;
    write_report_row(
        &mut file,
        "output",
        "subgraph",
        "dir",
        &display_path(&paths.subgraph_dir),
    )?;
    write_report_row(
        &mut file,
        "output",
        "workflow_round",
        "dir",
        &display_path(&paths.round_dir),
    )?;
    write_report_row(
        &mut file,
        "output",
        "draft",
        "path",
        &display_path(&paths.input_draft),
    )?;
    write_report_row(
        &mut file,
        "output",
        "reference",
        "path",
        &display_path(&paths.input_reference),
    )?;
    write_report_row(
        &mut file,
        "output",
        "validation_fasta",
        "path",
        &display_path(paths.validation_fasta()),
    )?;
    if options.validate_round == 1 {
        write_report_row(
            &mut file,
            "output",
            "polished_round_1",
            "path",
            &display_path(&paths.polished_round_fasta(1)),
        )?;
        write_report_row(
            &mut file,
            "output",
            "polished",
            "path",
            &display_path(&paths.polished_fasta),
        )?;
        write_report_row(
            &mut file,
            "output",
            "polished_aligned",
            "path",
            &display_path(&paths.aligned_fasta),
        )?;
    }
    write_report_row(
        &mut file,
        "output",
        "logs",
        "dir",
        &display_path(&paths.logs_dir),
    )?;
    write_report_row(
        &mut file,
        "output",
        "external_stderr",
        "path",
        &display_path(&paths.external_stderr),
    )?;
    write_report_row(
        &mut file,
        "output",
        "sv_evidence_round_1",
        "path",
        &display_path(&paths.round1_sv_whole_read_evidence),
    )?;
    write_report_row(
        &mut file,
        "output",
        "sv_group_stats_round_1",
        "path",
        &display_path(&paths.round1_sv_group_summary),
    )?;
    write_report_row(
        &mut file,
        "output",
        "sv_subgroup_stats_round_1",
        "path",
        &display_path(&paths.round1_sv_subgroup_summary),
    )?;
    write_report_row(
        &mut file,
        "output",
        "sv_read_index_round_1",
        "path",
        &display_path(&paths.round1_sv_group_ids),
    )?;
    write_report_row(
        &mut file,
        "output",
        "sv_coverage_round_1",
        "path",
        &display_path(&paths.round1_sv_coverage),
    )?;
    write_report_row(
        &mut file,
        "output",
        "sv_snv_indel_summary_round_1",
        "path",
        &display_path(&paths.round1_sv_support_summary),
    )?;
    write_report_row(
        &mut file,
        "output",
        "sv_high_subgroups_round_1",
        "path",
        &display_path(&paths.round1_sv_high_subgroup_report),
    )?;
    write_report_row(
        &mut file,
        "output",
        "sv_plot_script_round_1",
        "path",
        &display_path(&paths.round1_plot_script),
    )?;
    write_report_row(
        &mut file,
        "output",
        "round_1_data",
        "dir",
        &display_path(&paths.round1_sv_data_dir),
    )?;
    write_report_row(
        &mut file,
        "output",
        "round_1_plots",
        "dir",
        &display_path(&paths.round1_sv_plots_dir),
    )?;
    write_report_row(
        &mut file,
        "output",
        "round_1_reports",
        "dir",
        &display_path(&paths.round1_sv_reports_dir),
    )?;
    write_report_row(
        &mut file,
        "output",
        "snv_indel_calls_round_1",
        "path",
        &display_path(&paths.round1_snv_indel_per_variant_calls),
    )?;
    write_report_row(
        &mut file,
        "output",
        "snv_indel_segments_round_1",
        "path",
        &display_path(&paths.round1_snv_indel_segments),
    )?;
    write_report_row(
        &mut file,
        "output",
        "snv_indel_variants_round_1",
        "path",
        &display_path(&paths.round1_snv_indel_variant_type_annotations),
    )?;
    write_report_row(
        &mut file,
        "output",
        "snv_indel_variants_combined_round_1",
        "path",
        &display_path(&paths.round1_snv_indel_variant_type_annotations_combined),
    )?;
    write_report_row(
        &mut file,
        "output",
        "snv_indel_high_round_1",
        "path",
        &display_path(&paths.round1_snv_indel_variant_type_annotations_combined_high),
    )?;
    write_report_row(
        &mut file,
        "output",
        "snv_indel_plot_points_round_1",
        "path",
        &display_path(&paths.round1_snv_indel_plot_points),
    )?;
    write_report_row(
        &mut file,
        "output",
        "snv_indel_plot_script_round_1",
        "path",
        &display_path(&paths.round1_snv_indel_plot_script),
    )?;
    write_report_row(
        &mut file,
        "output",
        "snv_indel_runtime_summary_appended_to",
        "path",
        &display_path(&paths.round1_sv_support_summary),
    )?;
    write_report_row(
        &mut file,
        "output",
        "round_1",
        "dir",
        &display_path(&paths.round1_dir),
    )?;
    for record in stage_records {
        let name = format!("{}:{}", record.stage, record.round);
        write_report_row(&mut file, "stage", &name, "status", record.status)?;
        if let Some(elapsed_seconds) = record.elapsed_seconds {
            write_report_row(
                &mut file,
                "stage",
                &name,
                "elapsed_seconds",
                &format!("{elapsed_seconds:.3}"),
            )?;
        }
        write_report_row(&mut file, "stage", &name, "message", &record.message)?;
    }
    for report in round_reports {
        let name = format!("round_{}", report.round);
        write_report_row(
            &mut file,
            "round",
            &name,
            "input_length",
            &report.input_length.to_string(),
        )?;
        write_report_row(
            &mut file,
            "round",
            &name,
            "output_length",
            &report.output_length.to_string(),
        )?;
        write_report_row(
            &mut file,
            "round",
            &name,
            "alignments_used",
            &report.alignments_used.to_string(),
        )?;
        write_report_row(
            &mut file,
            "round",
            &name,
            "substitutions",
            &report.substitutions.to_string(),
        )?;
        write_report_row(
            &mut file,
            "round",
            &name,
            "deletions",
            &report.deletions.to_string(),
        )?;
        write_report_row(
            &mut file,
            "round",
            &name,
            "inserted_bases",
            &report.inserted_bases.to_string(),
        )?;
        write_report_row(
            &mut file,
            "round",
            &name,
            "low_coverage_bases",
            &report.low_coverage_bases.to_string(),
        )?;
    }
    if let Some(alignment_report) = alignment_report {
        write_alignment_report_rows(&mut file, alignment_report)?;
    }
    if let Some(report) = sv_eval_report {
        write_sv_eval_report_rows(&mut file, report)?;
    }
    if let Some(report) = snv_indel_report {
        write_snv_indel_report_rows(&mut file, report)?;
    }
    for record in command_records {
        let name = format!("{}:{}", record.stage, record.round);
        write_report_row(
            &mut file,
            "command",
            &name,
            "timestamp",
            &record.timestamp.to_string(),
        )?;
        write_report_row(&mut file, "command", &name, "status", record.status)?;
        write_report_row(
            &mut file,
            "command",
            &name,
            "elapsed_seconds",
            &format!("{:.3}", record.elapsed_seconds),
        )?;
        write_report_row(&mut file, "command", &name, "stdout", &record.stdout)?;
        write_report_row(&mut file, "command", &name, "stderr", &record.stderr)?;
        write_report_row(&mut file, "command", &name, "command", &record.command)?;
    }
    Ok(())
}

fn write_sv_eval_report_rows(file: &mut File, report: &SvEvalReport) -> Result<(), OrgraftError> {
    let name = "round_1";
    write_report_row(
        file,
        "sv_eval",
        name,
        "read_count",
        &report.read_count.to_string(),
    )?;
    write_report_row(
        file,
        "sv_eval",
        name,
        "paf_alignments",
        &report.paf_alignments.to_string(),
    )?;
    write_report_row(
        file,
        "sv_eval",
        name,
        "summary_rows",
        &report.summary_rows.to_string(),
    )?;
    write_report_row(
        file,
        "sv_eval",
        name,
        "no_alignment_reads",
        &report.no_alignment_reads.to_string(),
    )?;
    write_report_row(
        file,
        "sv_eval",
        name,
        "whole_read_evidence_rows",
        &report.whole_read_evidence_rows.to_string(),
    )?;
    write_report_row(
        file,
        "sv_eval",
        name,
        "fl_reads",
        &report.fl_reads.to_string(),
    )?;
    write_report_row(
        file,
        "sv_eval",
        name,
        "partial_reads",
        &report.partial_reads.to_string(),
    )?;
    write_report_row(
        file,
        "sv_eval",
        name,
        "reference_support_reads",
        &report.reference_support_reads.to_string(),
    )?;
    write_report_row(
        file,
        "sv_eval",
        name,
        "read_group_count",
        &report.read_group_count.to_string(),
    )?;
    write_report_row(
        file,
        "sv_eval",
        name,
        "read_subgroup_count",
        &report.read_subgroup_count.to_string(),
    )?;
    write_report_row(
        file,
        "sv_eval",
        name,
        "sv_support_status",
        &report.sv_support_status,
    )?;
    write_report_row(
        file,
        "sv_eval",
        name,
        "reference_support_rule",
        "type_1_subtype_NA_FL + circular-terminal type_2_subtype_ref_NA_FL + type_2_subtype_rep_NA_FL subgroups with mid_olp_1 >= 1000",
    )?;
    write_report_row(
        file,
        "sv_eval",
        name,
        "sv_support_pass_rule",
        "reference_support_depth_area_fraction >= 0.50 and low_green_window_fraction <= 0.05",
    )?;
    write_report_row(
        file,
        "sv_eval",
        name,
        "high_subgroup_judgement_rule",
        "high non-reference-support type>=3 subgroup with local reference_support_depth < 3 or local reference_support_fraction < 0.20 is possible_reference_sv_error; otherwise minor_recombination_or_alternative_configuration",
    )?;
    write_report_row(
        file,
        "sv_eval",
        name,
        "reference_support_rep_mid_olp_min",
        &format!("{REFERENCE_SUPPORT_REP_MID_OLP_MIN:.3}"),
    )?;
    write_report_row(
        file,
        "sv_eval",
        name,
        "sv_support_window_bp",
        &SV_SUPPORT_WINDOW_BP.to_string(),
    )?;
    write_report_row(
        file,
        "sv_eval",
        name,
        "sv_support_min_green_fraction",
        &format!("{SV_SUPPORT_MIN_GREEN_FRACTION:.6}"),
    )?;
    write_report_row(
        file,
        "sv_eval",
        name,
        "sv_support_max_low_green_window_fraction",
        &format!("{SV_SUPPORT_MAX_LOW_GREEN_WINDOW_FRACTION:.6}"),
    )?;
    write_report_row(
        file,
        "sv_eval",
        name,
        "sv_support_low_green_fraction",
        &format!("{SV_SUPPORT_LOW_GREEN_FRACTION:.6}"),
    )?;
    write_report_row(
        file,
        "sv_eval",
        name,
        "sv_support_min_green_depth",
        &format!("{SV_SUPPORT_MIN_GREEN_DEPTH:.3}"),
    )?;
    write_report_row(
        file,
        "sv_eval",
        name,
        "high_subgroup_min_fraction",
        &format!("{HIGH_SUBGROUP_MIN_FRACTION:.6}"),
    )?;
    write_report_row(
        file,
        "sv_eval",
        name,
        "auto_sv_plot_highlight_min_fraction",
        &format!("{:.6}", report.auto_sv_plot_highlight_min_fraction),
    )?;
    write_report_row(
        file,
        "sv_eval",
        name,
        "auto_sv_plot_highlight_min_reads",
        &report.auto_sv_plot_highlight_min_reads.to_string(),
    )?;
    write_report_row(
        file,
        "sv_eval",
        name,
        "auto_highlight_subgroup_count",
        &report.auto_highlight_subgroups.len().to_string(),
    )?;
    write_report_row(
        file,
        "sv_eval",
        name,
        "auto_highlight_subgroups",
        &if report.auto_highlight_subgroups.is_empty() {
            ".".to_string()
        } else {
            report.auto_highlight_subgroups.join(",")
        },
    )?;
    write_report_row(
        file,
        "sv_eval",
        name,
        "breakpoint_window_bp",
        &BREAKPOINT_WINDOW_BP.to_string(),
    )?;
    write_report_row(file, "sv_eval", name, "minimap2_mode", report.minimap2_mode)?;
    write_report_row(
        file,
        "sv_eval",
        name,
        "minimap2_workers",
        &report.minimap2_workers.to_string(),
    )?;
    write_report_row(
        file,
        "sv_eval",
        name,
        "whole_read_evidence_path",
        &display_path(&report.whole_read_evidence_path),
    )?;
    write_report_row(
        file,
        "sv_eval",
        name,
        "read_group_stats_path",
        &display_path(&report.read_group_summary_path),
    )?;
    write_report_row(
        file,
        "sv_eval",
        name,
        "read_subgroup_stats_path",
        &display_path(&report.read_subgroup_summary_path),
    )?;
    write_report_row(
        file,
        "sv_eval",
        name,
        "read_group_ids_path",
        &display_path(&report.read_group_ids_path),
    )?;
    write_report_row(
        file,
        "sv_eval",
        name,
        "coverage_path",
        &display_path(&report.coverage_path),
    )?;
    write_report_row(
        file,
        "sv_eval",
        name,
        "sv_snv_indel_summary_path",
        &display_path(&report.sv_support_summary_path),
    )?;
    write_report_row(
        file,
        "sv_eval",
        name,
        "high_subgroup_report_path",
        &display_path(&report.high_subgroup_report_path),
    )?;
    write_report_row(
        file,
        "sv_eval",
        name,
        "plot_script_path",
        &display_path(&report.plot_script_path),
    )?;
    write_report_row(
        file,
        "sv_eval",
        name,
        "minimap2_options",
        &report.minimap2_options,
    )?;
    write_report_row(
        file,
        "sv_eval",
        name,
        "terminal_extension_window",
        &TERMINAL_EXTENSION_WINDOW.to_string(),
    )?;
    write_report_row(
        file,
        "sv_eval",
        name,
        "terminal_extension_max_gap",
        &TERMINAL_EXTENSION_MAX_GAP.to_string(),
    )?;
    write_report_row(
        file,
        "sv_eval",
        name,
        "terminal_extension_min_alignment_length",
        &TERMINAL_EXTENSION_MIN_ALIGNMENT_LENGTH.to_string(),
    )?;
    Ok(())
}

fn write_snv_indel_report_rows(
    file: &mut File,
    report: &SnvIndelReport,
) -> Result<(), OrgraftError> {
    let name = "round_1";
    write_report_row(file, "snv_indel", name, "call_mode", report.call_mode)?;
    write_report_row(
        file,
        "snv_indel",
        name,
        "overlap_policy",
        report.overlap_policy,
    )?;
    write_report_row(
        file,
        "snv_indel",
        name,
        "sv_context_filter",
        &report.sv_context_filter.to_string(),
    )?;
    write_report_row(
        file,
        "snv_indel",
        name,
        "minimap2_preset",
        report.minimap2_preset,
    )?;
    write_report_row(
        file,
        "snv_indel",
        name,
        "fl_read_count",
        &report.fl_read_count.to_string(),
    )?;
    write_report_row(
        file,
        "snv_indel",
        name,
        "segment_count",
        &report.segment_count.to_string(),
    )?;
    write_report_row(
        file,
        "snv_indel",
        name,
        "total_calls",
        &report.total_calls.to_string(),
    )?;
    write_report_row(
        file,
        "snv_indel",
        name,
        "reads_with_calls",
        &report.reads_with_calls.to_string(),
    )?;
    write_report_row(
        file,
        "snv_indel",
        name,
        "failed_segments",
        &report.failed_segments.to_string(),
    )?;
    write_report_row(
        file,
        "snv_indel",
        name,
        "workers",
        &report.worker_count.to_string(),
    )?;
    write_report_row(
        file,
        "snv_indel",
        name,
        "elapsed_seconds",
        &format!("{:.3}", report.elapsed_seconds),
    )?;
    write_report_row(
        file,
        "snv_indel",
        name,
        "shared_minimap2_stream_seconds",
        &format!("{:.3}", report.alignment_seconds),
    )?;
    write_report_row(
        file,
        "snv_indel",
        name,
        "sum_segment_seconds",
        &format!("{:.3}", report.sum_segment_seconds),
    )?;
    for (metric, value) in snv_write_timing_rows(&report.write_timings) {
        write_report_row(file, "snv_indel", name, &metric, &value)?;
    }
    write_report_row(
        file,
        "snv_indel",
        name,
        "reference_path",
        &display_path(&report.reference_path),
    )?;
    write_report_row(
        file,
        "snv_indel",
        name,
        "per_variant_calls_path",
        &display_path(&report.per_variant_calls_path),
    )?;
    write_report_row(
        file,
        "snv_indel",
        name,
        "segments_path",
        &display_path(&report.segments_path),
    )?;
    write_report_row(
        file,
        "snv_indel",
        name,
        "variants_path",
        &display_path(&report.variant_type_annotations_path),
    )?;
    write_report_row(
        file,
        "snv_indel",
        name,
        "variants_combined_path",
        &display_path(&report.variant_type_annotations_combined_path),
    )?;
    write_report_row(
        file,
        "snv_indel",
        name,
        "high_variants_path",
        &display_path(&report.variant_type_annotations_combined_high_path),
    )?;
    write_report_row(
        file,
        "snv_indel",
        name,
        "plot_points_path",
        &display_path(&report.plot_points_path),
    )?;
    write_report_row(
        file,
        "snv_indel",
        name,
        "plot_script_path",
        &display_path(&report.plot_script_path),
    )?;
    write_report_row(
        file,
        "snv_indel",
        name,
        "runtime_summary_appended_to",
        &display_path(&report.summary_path),
    )?;
    write_report_row(
        file,
        "snv_indel",
        name,
        "read_split_rule",
        "all FL reads from round_1/01.data/sv_read_index.tsv; multi-alignment reads are split by round_1 alignment_summary qs/qe into P1..Pn; overlap handling follows overlap_policy; when sv_context_filter=true, multi-mapping segment SAM records are selected by nearest round_1 subject interval",
    )?;
    write_report_row(
        file,
        "snv_indel",
        name,
        "per_variant_calls_format",
        "read_id, segment_id, read metadata, pos, ref, alt, type, confidence, confidence_reason; custom caller marks dense/long indel neighborhoods and segment-boundary indels",
    )?;
    Ok(())
}

fn write_alignment_report_rows(
    file: &mut File,
    report: &AlignmentReport,
) -> Result<(), OrgraftError> {
    let name = "polished_aln";
    write_report_row(file, "alignment", name, "record_id", &report.record_id)?;
    write_report_row(
        file,
        "alignment",
        name,
        "reference_id",
        &report.reference_id,
    )?;
    write_report_row(
        file,
        "alignment",
        name,
        "input_length",
        &report.input_length.to_string(),
    )?;
    write_report_row(
        file,
        "alignment",
        name,
        "output_length",
        &report.output_length.to_string(),
    )?;
    write_report_row(
        file,
        "alignment",
        name,
        "orientation",
        &report.orientation.to_string(),
    )?;
    write_report_row(
        file,
        "alignment",
        name,
        "reverse_complemented",
        &report.reverse_complemented.to_string(),
    )?;
    write_report_row(
        file,
        "alignment",
        name,
        "rotation_step",
        &report.rotation_step.to_string(),
    )?;
    write_report_row(
        file,
        "alignment",
        name,
        "best_pident",
        &format!("{:.3}", report.best_pident),
    )?;
    write_report_row(
        file,
        "alignment",
        name,
        "best_length",
        &report.best_length.to_string(),
    )?;
    write_report_row(
        file,
        "alignment",
        name,
        "query_start",
        &report.query_start.to_string(),
    )?;
    write_report_row(
        file,
        "alignment",
        name,
        "subject_start",
        &report.subject_start.to_string(),
    )?;
    write_report_row(
        file,
        "alignment",
        name,
        "subject_end",
        &report.subject_end.to_string(),
    )?;
    Ok(())
}

fn write_report_row(
    file: &mut File,
    section: &str,
    name: &str,
    metric: &str,
    value: &str,
) -> Result<(), OrgraftError> {
    writeln!(
        file,
        "{}\t{}\t{}\t{}",
        sanitize_tsv(section),
        sanitize_tsv(name),
        sanitize_tsv(metric),
        sanitize_tsv(value)
    )?;
    Ok(())
}

fn sanitize_tsv(value: &str) -> String {
    value.replace(['\t', '\n', '\r'], " ")
}

fn push_sanitized_tsv(buffer: &mut String, value: &str) {
    if value
        .bytes()
        .any(|byte| matches!(byte, b'\t' | b'\n' | b'\r'))
    {
        for ch in value.chars() {
            match ch {
                '\t' | '\n' | '\r' => buffer.push(' '),
                _ => buffer.push(ch),
            }
        }
    } else {
        buffer.push_str(value);
    }
}

fn same_pos_types_by_pos(
    compatible_type_by_pos: &BTreeMap<usize, BTreeSet<String>>,
) -> BTreeMap<usize, String> {
    compatible_type_by_pos
        .keys()
        .map(|pos| (*pos, same_pos_types_label(compatible_type_by_pos, *pos)))
        .collect()
}

fn write_built_tsv(
    path: &Path,
    buffer: String,
    started_at: Instant,
) -> Result<TsvWriteTiming, OrgraftError> {
    let bytes = buffer.len();
    fs::write(path, buffer)?;
    Ok(TsvWriteTiming {
        total_seconds: started_at.elapsed().as_secs_f64(),
        bytes,
    })
}

fn escape_tsv_cell(value: &str) -> String {
    value
        .replace('\\', "\\\\")
        .replace('\t', "\\t")
        .replace('\n', "\\n")
        .replace('\r', "\\r")
}

fn resolve_reference_id(
    inputs: &ResolvedInputs,
    subgraph: &str,
    draft: &FastaRecordInfo,
) -> Result<String, OrgraftError> {
    if let Some(reference_id) = reference_id_from_subgraph_header(&draft.header) {
        return Ok(reference_id);
    }
    let reference_ids = fasta_record_ids(&inputs.reference)?;
    match reference_ids.as_slice() {
        [only] => Ok(only.clone()),
        _ => Err(OrgraftError::InvalidArgument(format!(
            "could not infer reference record for `{subgraph}`; include reference=... in the subgraph FASTA header or pass a single-record reference FASTA"
        ))),
    }
}

fn reference_id_from_subgraph_header(header: &str) -> Option<String> {
    let start = header.find("reference=")? + "reference=".len();
    let rest = &header[start..];
    let end = rest
        .find(|ch: char| ch == ';' || ch == ']' || ch.is_whitespace())
        .unwrap_or(rest.len());
    let reference_id = &rest[..end];
    (!reference_id.is_empty()).then(|| reference_id.to_string())
}

fn extract_fasta_record_by_id(
    source: &Path,
    target: &Path,
    wanted_id: &str,
    label: &str,
) -> Result<FastaRecordInfo, OrgraftError> {
    let text = fs::read_to_string(source)?;
    let mut records = Vec::new();
    let mut current_header: Option<String> = None;
    let mut current_sequence = String::new();

    for line in text.lines() {
        if let Some(header) = line.strip_prefix('>') {
            if let Some(previous_header) = current_header.replace(header.to_string()) {
                records.push((previous_header, std::mem::take(&mut current_sequence)));
            }
        } else if current_header.is_some() {
            current_sequence.push_str(line.trim());
        }
    }
    if let Some(header) = current_header {
        records.push((header, current_sequence));
    }

    for (header, sequence) in records {
        let id = fasta_id(&header);
        if id == wanted_id {
            let mut file = File::create(target)?;
            writeln!(file, ">{header}")?;
            write_wrapped_sequence(&mut file, &sequence)?;
            return Ok(FastaRecordInfo { id, header });
        }
    }

    Err(OrgraftError::InvalidArgument(format!(
        "{} {} does not contain FASTA record `{wanted_id}`",
        label,
        source.display()
    )))
}

fn fasta_record_ids(path: &Path) -> Result<Vec<String>, OrgraftError> {
    let text = fs::read_to_string(path)?;
    Ok(text
        .lines()
        .filter_map(|line| line.strip_prefix('>'))
        .map(fasta_id)
        .collect())
}

fn fasta_id(header: &str) -> String {
    header
        .split_whitespace()
        .next()
        .unwrap_or(header)
        .to_string()
}

fn write_wrapped_sequence(file: &mut File, sequence: &str) -> Result<(), OrgraftError> {
    for chunk in sequence.as_bytes().chunks(80) {
        file.write_all(chunk)?;
        writeln!(file)?;
    }
    Ok(())
}

fn read_single_fasta_record(path: &Path) -> Result<(String, String), OrgraftError> {
    let text = fs::read_to_string(path)?;
    let mut current_header: Option<String> = None;
    let mut current_sequence = String::new();
    let mut records = Vec::new();
    for line in text.lines() {
        if let Some(header) = line.strip_prefix('>') {
            if let Some(previous_header) = current_header.replace(header.to_string()) {
                records.push((
                    fasta_id(&previous_header),
                    std::mem::take(&mut current_sequence),
                ));
            }
        } else if current_header.is_some() {
            current_sequence.push_str(&line.trim().to_ascii_uppercase());
        }
    }
    if let Some(header) = current_header {
        records.push((fasta_id(&header), current_sequence));
    }
    match records.len() {
        1 => Ok(records.remove(0)),
        0 => Err(OrgraftError::InvalidArgument(format!(
            "{} contains no FASTA records",
            path.display()
        ))),
        _ => Err(OrgraftError::InvalidArgument(format!(
            "{} contains multiple FASTA records; polish expects one linear sequence",
            path.display()
        ))),
    }
}

fn read_fasta_records_by_id(path: &Path) -> Result<HashMap<String, String>, OrgraftError> {
    let text = fs::read_to_string(path)?;
    let mut current_header: Option<String> = None;
    let mut current_sequence = String::new();
    let mut records = HashMap::new();
    for line in text.lines() {
        if let Some(header) = line.strip_prefix('>') {
            if let Some(previous_header) = current_header.replace(header.to_string()) {
                records.insert(
                    fasta_id(&previous_header),
                    std::mem::take(&mut current_sequence).to_ascii_uppercase(),
                );
            }
        } else if current_header.is_some() {
            current_sequence.push_str(line.trim());
        }
    }
    if let Some(header) = current_header {
        records.insert(fasta_id(&header), current_sequence.to_ascii_uppercase());
    }
    if records.is_empty() {
        return Err(OrgraftError::InvalidArgument(format!(
            "{} contains no FASTA records",
            path.display()
        )));
    }
    Ok(records)
}

fn read_sequence_records(path: &Path) -> Result<Vec<ReadRecord>, OrgraftError> {
    with_text_reader(path, |reader| {
        let mut first = String::new();
        if reader.read_line(&mut first)? == 0 {
            return Ok(Vec::new());
        }
        if first.starts_with('@') {
            read_fastq_records(reader, first)
        } else if first.starts_with('>') {
            read_fasta_records_from_reader(reader, first)
        } else {
            Err(OrgraftError::InvalidArgument(format!(
                "{} is not FASTQ or FASTA",
                path.display()
            )))
        }
    })
}

fn with_text_reader<F, T>(path: &Path, callback: F) -> Result<T, OrgraftError>
where
    F: FnOnce(&mut dyn BufRead) -> Result<T, OrgraftError>,
{
    if is_gzip_path(path) {
        let mut child = Command::new("gzip")
            .arg("-cd")
            .arg(path)
            .stdout(Stdio::piped())
            .spawn()
            .map_err(|error| {
                OrgraftError::InvalidArgument(format!(
                    "failed to spawn gzip for {}: {error}",
                    path.display()
                ))
            })?;
        let stdout = child.stdout.take().ok_or_else(|| {
            OrgraftError::InvalidArgument(format!(
                "failed to capture gzip stdout for {}",
                path.display()
            ))
        })?;
        let mut reader = BufReader::new(stdout);
        let result = callback(&mut reader);
        let status = child.wait()?;
        if !status.success() {
            return Err(OrgraftError::InvalidArgument(format!(
                "gzip failed while reading {}",
                path.display()
            )));
        }
        result
    } else {
        let file = File::open(path)?;
        let mut reader = BufReader::new(file);
        callback(&mut reader)
    }
}

fn read_fastq_records(
    reader: &mut dyn BufRead,
    first_header: String,
) -> Result<Vec<ReadRecord>, OrgraftError> {
    let mut records = Vec::new();
    let mut header = first_header;
    loop {
        let mut sequence = String::new();
        let mut plus = String::new();
        let mut quality = String::new();
        if !header.starts_with('@') {
            return Err(OrgraftError::InvalidArgument(format!(
                "invalid FASTQ header `{}`",
                header.trim_end()
            )));
        }
        if reader.read_line(&mut sequence)? == 0
            || reader.read_line(&mut plus)? == 0
            || reader.read_line(&mut quality)? == 0
        {
            return Err(OrgraftError::InvalidArgument(
                "truncated FASTQ record".to_string(),
            ));
        }
        if !plus.starts_with('+') {
            return Err(OrgraftError::InvalidArgument(format!(
                "invalid FASTQ plus line `{}`",
                plus.trim_end()
            )));
        }
        records.push(ReadRecord {
            id: read_header_id(&header[1..]),
            sequence: sequence.trim_end().to_ascii_uppercase(),
        });

        header.clear();
        if reader.read_line(&mut header)? == 0 {
            break;
        }
    }
    Ok(records)
}

fn read_fasta_records_from_reader(
    reader: &mut dyn BufRead,
    first_header: String,
) -> Result<Vec<ReadRecord>, OrgraftError> {
    let mut records = Vec::new();
    let mut current_id = read_header_id(first_header.trim_start_matches('>'));
    let mut current_sequence = String::new();
    for line in reader.lines() {
        let line = line?;
        if let Some(header) = line.strip_prefix('>') {
            records.push(ReadRecord {
                id: std::mem::take(&mut current_id),
                sequence: std::mem::take(&mut current_sequence).to_ascii_uppercase(),
            });
            current_id = read_header_id(header);
        } else {
            current_sequence.push_str(line.trim());
        }
    }
    records.push(ReadRecord {
        id: current_id,
        sequence: current_sequence.to_ascii_uppercase(),
    });
    Ok(records)
}

fn read_header_id(header: &str) -> String {
    normalize_read_id(header.split_whitespace().next().unwrap_or(header).trim())
}

fn normalize_read_id(id: &str) -> String {
    id.split_whitespace()
        .next()
        .unwrap_or(id)
        .trim()
        .replace('/', "_")
        .to_string()
}

fn write_single_fasta(path: &Path, id: &str, sequence: &str) -> Result<(), OrgraftError> {
    let mut file = File::create(path)?;
    writeln!(file, ">{id}")?;
    write_wrapped_sequence(&mut file, sequence)
}

fn parse_cigar(cigar: &str) -> Result<Vec<(usize, char)>, OrgraftError> {
    let mut result = Vec::new();
    let mut number = String::new();
    for ch in cigar.chars() {
        if ch.is_ascii_digit() {
            number.push(ch);
            continue;
        }
        if number.is_empty() {
            return Err(OrgraftError::InvalidArgument(format!(
                "invalid CIGAR `{cigar}`"
            )));
        }
        let len = parse_usize_value(&number, "CIGAR length")?;
        result.push((len, ch));
        number.clear();
    }
    if !number.is_empty() {
        return Err(OrgraftError::InvalidArgument(format!(
            "invalid trailing CIGAR length in `{cigar}`"
        )));
    }
    Ok(result)
}

#[derive(Debug, Clone)]
struct BlastHit {
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

fn run_blastn(
    blastn: &Path,
    query: &Path,
    subject: &Path,
    output: &Path,
    stderr_path: &Path,
    stage: &'static str,
    round: &str,
    commands: &mut Vec<CommandRecord>,
) -> Result<(), OrgraftError> {
    let mut stderr_file = OpenOptions::new()
        .create(true)
        .append(true)
        .open(stderr_path)?;
    writeln!(stderr_file, "### {stage}:{round} blastn stderr ###")?;
    let stderr_for_child = stderr_file.try_clone()?;
    let stdout_file = File::create(output)?;
    let outfmt = "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen";
    let mut command = Command::new(blastn);
    command
        .arg("-query")
        .arg(query)
        .arg("-subject")
        .arg(subject)
        .arg("-outfmt")
        .arg(outfmt);
    let command_text = format!("{command:?}");
    let started = Instant::now();
    let status = command
        .stdout(Stdio::from(stdout_file))
        .stderr(Stdio::from(stderr_for_child))
        .status()?;
    let elapsed_seconds = started.elapsed().as_secs_f64();
    let status_text = if status.success() { "ok" } else { "failed" };
    writeln!(
        OpenOptions::new().append(true).open(stderr_path)?,
        "### {stage}:{round} status={status_text} elapsed_seconds={elapsed_seconds:.3} ###\n"
    )?;
    commands.push(CommandRecord {
        timestamp: timestamp(),
        stage,
        round: round.to_string(),
        status: status_text,
        elapsed_seconds,
        stdout: display_path(output),
        stderr: display_path(stderr_path),
        command: command_text,
    });
    if status.success() {
        Ok(())
    } else {
        Err(OrgraftError::InvalidArgument(format!(
            "{stage} {round} failed; see {}",
            stderr_path.display()
        )))
    }
}

fn parse_blast_hits(path: &Path) -> Result<Vec<BlastHit>, OrgraftError> {
    let text = fs::read_to_string(path)?;
    let mut hits = Vec::new();
    for (index, line) in text
        .lines()
        .filter(|line| !line.trim().is_empty())
        .enumerate()
    {
        let fields = line.split('\t').collect::<Vec<_>>();
        if fields.len() < 12 {
            continue;
        }
        hits.push(BlastHit {
            subject_id: fields[1].to_string(),
            pident: parse_f64_value(fields[2], path, index + 1)?,
            length: parse_usize_value_at(fields[3], path, index + 1)?,
            query_start: parse_usize_value_at(fields[6], path, index + 1)?,
            subject_start: parse_usize_value_at(fields[8], path, index + 1)?,
            subject_end: parse_usize_value_at(fields[9], path, index + 1)?,
            bitscore: parse_f64_value(fields[11], path, index + 1)?,
        });
    }
    Ok(hits)
}

fn best_hit_for_subject(hits: &[BlastHit], subject_id: &str) -> Option<BlastHit> {
    hits.iter()
        .filter(|hit| hit.subject_id == subject_id)
        .max_by(|left, right| blast_hit_order(left, right))
        .cloned()
        .or_else(|| {
            hits.iter()
                .max_by(|left, right| blast_hit_order(left, right))
                .cloned()
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

fn rotate_sequence(sequence: &str, step: isize) -> String {
    if sequence.is_empty() {
        return String::new();
    }
    let len = sequence.len() as isize;
    let normalized = ((step % len) + len) % len;
    let split = normalized as usize;
    format!("{}{}", &sequence[split..], &sequence[..split])
}

fn slice_forward(sequence: &str, start: usize, len: usize) -> &str {
    let start = start.min(sequence.len());
    let end = (start + len).min(sequence.len());
    &sequence[start..end]
}

fn slice_before(sequence: &str, end: usize, len: usize) -> &str {
    let end = end.min(sequence.len());
    let start = end.saturating_sub(len);
    &sequence[start..end]
}

fn reverse_string(sequence: &str) -> String {
    sequence.chars().rev().collect()
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

fn base_index(base: u8) -> Option<usize> {
    match base.to_ascii_uppercase() {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' => Some(3),
        _ => None,
    }
}

fn index_base(index: usize) -> u8 {
    match index {
        0 => b'A',
        1 => b'C',
        2 => b'G',
        _ => b'T',
    }
}

fn best_base_count(counts: [u32; 4]) -> (usize, u32) {
    let mut best_index = 0usize;
    let mut best_count = counts[0];
    for (index, count) in counts.iter().copied().enumerate().skip(1) {
        if count > best_count {
            best_index = index;
            best_count = count;
        }
    }
    (best_index, best_count)
}

fn round_to_three_decimals(value: f64) -> f64 {
    (value * 1000.0).round() / 1000.0
}

fn format_py_float(value: f64) -> String {
    if value.is_finite() && (value.fract()).abs() < 1e-9 {
        format!("{value:.1}")
    } else {
        value.to_string()
    }
}

fn read_soft_paths(path: &Path) -> Result<HashMap<String, PathBuf>, OrgraftError> {
    let text = fs::read_to_string(path)?;
    let mut paths = HashMap::new();
    for line in text.lines() {
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        let mut fields = trimmed.split_whitespace();
        let Some(name) = fields.next() else {
            continue;
        };
        let Some(tool_path) = fields.next() else {
            continue;
        };
        paths.insert(name.to_string(), PathBuf::from(tool_path));
    }
    Ok(paths)
}

fn require_tool(paths: &HashMap<String, PathBuf>, name: &str) -> Result<PathBuf, OrgraftError> {
    paths.get(name).cloned().ok_or_else(|| {
        OrgraftError::InvalidArgument(format!("soft_paths is missing required tool `{name}`"))
    })
}

fn required_value<'a>(
    args: &'a [String],
    index: &mut usize,
    option: &str,
) -> Result<&'a str, OrgraftError> {
    *index += 1;
    args.get(*index)
        .map(|value| value.as_str())
        .ok_or_else(|| OrgraftError::InvalidArgument(format!("missing value for {option}")))
}

fn parse_label(value: &str, option: &str) -> Result<String, OrgraftError> {
    if value.is_empty()
        || value
            .chars()
            .any(|ch| !(ch.is_ascii_alphanumeric() || ch == '_' || ch == '-' || ch == '.'))
    {
        return Err(OrgraftError::InvalidArgument(format!(
            "invalid {option} `{value}`; use letters, numbers, '.', '_' or '-'"
        )));
    }
    Ok(value.to_string())
}

fn parse_usize(value: &str, option: &str) -> Result<usize, OrgraftError> {
    value.parse::<usize>().map_err(|_| {
        OrgraftError::InvalidArgument(format!("invalid value for {option}: `{value}`"))
    })
}

fn parse_fraction(value: &str, option: &str) -> Result<f64, OrgraftError> {
    let parsed = parse_f64_label(value, option)?;
    if (0.0..=1.0).contains(&parsed) {
        Ok(parsed)
    } else {
        Err(OrgraftError::InvalidArgument(format!(
            "invalid value for {option}: `{value}`; expected fraction between 0 and 1"
        )))
    }
}

fn parse_plot_range(value: &str, option: &str) -> Result<PlotRange, OrgraftError> {
    let Some((start, end)) = value.split_once('-').or_else(|| value.split_once(':')) else {
        return Err(OrgraftError::InvalidArgument(format!(
            "invalid value for {option}: `{value}`; expected START-END"
        )));
    };
    let start = parse_usize(start.trim(), option)?;
    let end = parse_usize(end.trim(), option)?;
    if start == 0 || end < start {
        return Err(OrgraftError::InvalidArgument(format!(
            "invalid value for {option}: `{value}`; expected 1-based START-END with START <= END"
        )));
    }
    Ok(PlotRange { start, end })
}

fn parse_on_off(value: &str, option: &str) -> Result<bool, OrgraftError> {
    match value {
        "on" => Ok(true),
        "off" => Ok(false),
        other => Err(OrgraftError::InvalidArgument(format!(
            "invalid value for {option}: `{other}`; expected on or off"
        ))),
    }
}

fn parse_comma_list(value: &str) -> Vec<String> {
    value
        .split(',')
        .map(str::trim)
        .filter(|item| !item.is_empty())
        .map(ToString::to_string)
        .collect()
}

fn parse_usize_value(value: &str, label: &str) -> Result<usize, OrgraftError> {
    value
        .parse::<usize>()
        .map_err(|_| OrgraftError::InvalidArgument(format!("invalid {label}: `{value}`")))
}

fn parse_usize_value_at(value: &str, path: &Path, line: usize) -> Result<usize, OrgraftError> {
    value.parse::<usize>().map_err(|_| {
        OrgraftError::InvalidArgument(format!(
            "{}:{line} expected integer BLAST field, got `{value}`",
            path.display()
        ))
    })
}

fn parse_f64_value(value: &str, path: &Path, line: usize) -> Result<f64, OrgraftError> {
    value.parse::<f64>().map_err(|_| {
        OrgraftError::InvalidArgument(format!(
            "{}:{line} expected numeric BLAST field, got `{value}`",
            path.display()
        ))
    })
}

fn parse_f64_label(value: &str, label: &str) -> Result<f64, OrgraftError> {
    value
        .parse::<f64>()
        .map_err(|_| OrgraftError::InvalidArgument(format!("invalid {label}: `{value}`")))
}

#[cfg(test)]
fn parse_isize_label(value: &str, label: &str) -> Result<isize, OrgraftError> {
    value
        .parse::<isize>()
        .map_err(|_| OrgraftError::InvalidArgument(format!("invalid {label}: `{value}`")))
}

fn parse_u16(value: &str, label: &str) -> Result<u16, OrgraftError> {
    value
        .parse::<u16>()
        .map_err(|_| OrgraftError::InvalidArgument(format!("invalid {label}: `{value}`")))
}

fn parse_u8(value: &str, label: &str) -> Result<u8, OrgraftError> {
    value
        .parse::<u8>()
        .map_err(|_| OrgraftError::InvalidArgument(format!("invalid {label}: `{value}`")))
}

fn canonicalize_existing(path: &Path, label: &str) -> Result<PathBuf, OrgraftError> {
    path.canonicalize().map_err(|error| {
        OrgraftError::InvalidArgument(format!(
            "{} {} could not be read: {}",
            label,
            path.display(),
            error
        ))
    })
}

fn display_path(path: &Path) -> String {
    path.display().to_string()
}

fn is_gzip_path(path: &Path) -> bool {
    let mut magic = [0u8; 2];
    File::open(path)
        .and_then(|mut file| file.read_exact(&mut magic))
        .is_ok_and(|_| magic == [0x1f, 0x8b])
}

fn timestamp() -> String {
    SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map(|duration| duration.as_secs().to_string())
        .unwrap_or_else(|_| "0".to_string())
}

fn temp_file_path(label: &str, extension: &str) -> PathBuf {
    let stamp = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map(|duration| duration.as_nanos())
        .unwrap_or(0);
    std::env::temp_dir().join(format!(
        "{label}-{}-{stamp}.{extension}",
        std::process::id()
    ))
}

#[cfg(unix)]
fn link_or_copy(source: &Path, target: &Path) -> Result<(), OrgraftError> {
    std::os::unix::fs::symlink(source, target)?;
    Ok(())
}

#[cfg(not(unix))]
fn link_or_copy(source: &Path, target: &Path) -> Result<(), OrgraftError> {
    fs::copy(source, target)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_explicit_subgraph_layout() {
        let args = vec![
            "--organelle".to_string(),
            "mito".to_string(),
            "--subgraph".to_string(),
            "subgraph_002".to_string(),
            "--draft".to_string(),
            "draft.fa".to_string(),
            "--reference".to_string(),
            "ref.fa".to_string(),
            "--reads".to_string(),
            "reads.fastq.gz".to_string(),
        ];
        let options = PolishOptions::from_args(&args).unwrap();
        let paths = PolishPaths::new(
            &options.out_dir,
            &options.organelle,
            &options.subgraph,
            options.validate_round,
        );
        assert_eq!(
            paths.subgraph_dir,
            PathBuf::from("results/polish/mito/subgraph_002")
        );
        assert_eq!(
            paths.round_dir,
            PathBuf::from("results/polish/mito/subgraph_002/round_1")
        );
        assert_eq!(
            paths.round1_snv_indel_dir,
            PathBuf::from("results/polish/mito/subgraph_002/round_1/03.validate")
        );
        assert_eq!(
            paths.round1_snv_indel_data_dir,
            PathBuf::from("results/polish/mito/subgraph_002/round_1/03.validate/01.data")
        );
        assert_eq!(
            paths.round1_snv_indel_reports_dir,
            PathBuf::from("results/polish/mito/subgraph_002/round_1/03.validate/03.reports")
        );
        assert_eq!(
            paths.round1_snv_indel_plots_dir,
            PathBuf::from("results/polish/mito/subgraph_002/round_1/03.validate/02.plots")
        );
        assert_eq!(
            paths.report,
            PathBuf::from("results/polish/mito/subgraph_002/round_1/logs/report.tsv")
        );
        assert_eq!(
            paths.external_stderr,
            PathBuf::from("results/polish/mito/subgraph_002/round_1/logs/external.stderr.log")
        );
        assert_eq!(
            paths.aligned_fasta,
            PathBuf::from("results/polish/mito/subgraph_002/round_1/02.polish/polished_aln.fasta")
        );
        assert_eq!(
            paths.polished_fasta,
            PathBuf::from("results/polish/mito/subgraph_002/round_1/02.polish/polished.fasta")
        );
        assert_eq!(
            paths.round1_sv_whole_read_evidence,
            PathBuf::from(
                "results/polish/mito/subgraph_002/round_1/03.validate/01.data/sv_read_evidence.tsv"
            )
        );
        assert_eq!(
            paths.round1_sv_group_summary,
            PathBuf::from(
                "results/polish/mito/subgraph_002/round_1/03.validate/03.reports/sv_group_stats.tsv"
            )
        );
        assert_eq!(
            paths.round1_sv_support_summary,
            PathBuf::from(
                "results/polish/mito/subgraph_002/round_1/03.validate/03.reports/sv_snv_indel_summary.tsv"
            )
        );
        assert_eq!(
            paths.round1_sv_high_subgroup_report,
            PathBuf::from(
                "results/polish/mito/subgraph_002/round_1/03.validate/03.reports/sv_high_subgroups.tsv"
            )
        );
        assert_eq!(
            paths.round1_sv_coverage,
            PathBuf::from(
                "results/polish/mito/subgraph_002/round_1/03.validate/01.data/sv_coverage.tsv"
            )
        );
        assert_eq!(
            paths.round1_plot_script,
            PathBuf::from(
                "results/polish/mito/subgraph_002/round_1/03.validate/02.plots/plot_sv_support.py"
            )
        );
        assert_eq!(
            paths.round1_snv_indel_variant_type_annotations,
            PathBuf::from(
                "results/polish/mito/subgraph_002/round_1/03.validate/01.data/snv_indel_variants.tsv"
            )
        );
        assert_eq!(
            paths.round1_snv_indel_variant_type_annotations_combined,
            PathBuf::from(
                "results/polish/mito/subgraph_002/round_1/03.validate/01.data/snv_indel_variants_combined.tsv"
            )
        );
        assert_eq!(
            paths.round1_snv_indel_variant_type_annotations_combined_high,
            PathBuf::from(
                "results/polish/mito/subgraph_002/round_1/03.validate/03.reports/snv_indel_high.tsv"
            )
        );
        assert_eq!(
            paths.round1_snv_indel_plot_points,
            PathBuf::from(
                "results/polish/mito/subgraph_002/round_1/03.validate/01.data/snv_indel_plot_points.tsv"
            )
        );
        assert_eq!(
            paths.round1_snv_indel_plot_script,
            PathBuf::from(
                "results/polish/mito/subgraph_002/round_1/03.validate/02.plots/plot_snv_indel.py"
            )
        );
        assert!(options.per_read_variant_calls);
        assert_eq!(
            options.snv_indel_overlap_policy,
            SnvIndelOverlapPolicy::MarkOverlap
        );
    }

    #[test]
    fn validate_round_changes_internal_validate_directory_only() {
        let args = vec![
            "--organelle".to_string(),
            "mito".to_string(),
            "--subgraph".to_string(),
            "subgraph_001".to_string(),
            "--draft".to_string(),
            "draft.fa".to_string(),
            "--reference".to_string(),
            "ref.fa".to_string(),
            "--reads".to_string(),
            "reads.fastq.gz".to_string(),
            "--out-dir".to_string(),
            "results_workflow/04.polish".to_string(),
            "--validate-round".to_string(),
            "2".to_string(),
        ];
        let options = PolishOptions::from_args(&args).unwrap();
        let paths = PolishPaths::new(
            &options.out_dir,
            &options.organelle,
            &options.subgraph,
            options.validate_round,
        );
        assert_eq!(options.validate_round, 2);
        assert_eq!(
            paths.subgraph_dir,
            PathBuf::from("results_workflow/04.polish/mito/subgraph_001")
        );
        assert_eq!(
            paths.round_dir,
            PathBuf::from("results_workflow/04.polish/mito/subgraph_001/round_2")
        );
        assert_eq!(
            paths.round1_snv_indel_dir,
            PathBuf::from("results_workflow/04.polish/mito/subgraph_001/round_2/03.validate")
        );
        assert_eq!(
            paths.round1_sv_support_summary,
            PathBuf::from(
                "results_workflow/04.polish/mito/subgraph_001/round_2/03.validate/03.reports/sv_snv_indel_summary.tsv"
            )
        );
        assert_eq!(
            paths.input_draft,
            PathBuf::from(
                "results_workflow/04.polish/mito/subgraph_001/round_2/01.inputs/linear_subgraph.round_2.fasta"
            )
        );
        assert_eq!(
            paths.input_reference,
            PathBuf::from(
                "results_workflow/04.polish/mito/subgraph_001/round_2/01.inputs/rotated_reference.fasta"
            )
        );
        assert_eq!(
            paths.input_reads,
            PathBuf::from(
                "results_workflow/04.polish/mito/subgraph_001/round_2/01.inputs/subgraph_reads.fastq.gz"
            )
        );
        assert_eq!(
            paths.validation_fasta(),
            Path::new("results_workflow/04.polish/mito/subgraph_001/round_2/01.inputs/linear_subgraph.round_2.fasta")
        );
        assert_eq!(
            paths.report,
            PathBuf::from("results_workflow/04.polish/mito/subgraph_001/round_2/logs/report.tsv")
        );
    }

    #[test]
    fn parses_alignment_summary_group_and_subgroup() {
        let line = "readA\tmito_adj_2_1\t1000\taln_type=2;rep,NA\t99.0\taln=1;len=600;olp=100;idt=100.0;strand=+;qs=1;qe=600;ss=1;se=600;cn=1;c1=100.0,1,1\taln=2;len=500;olp=NA;idt=99.9;strand=+;qs=501;qe=1000;ss=501;se=1000;cn=1;c1=100.0,1,1";
        let record = parse_alignment_summary_record(line, 1).unwrap();
        assert_eq!(record.group_name(), "type_2_subtype_rep_NA");
        assert!(record.is_fl());
        assert_eq!(record.subgroup_key().unwrap().label(), "se1=600,ss2=501");
        assert_eq!(record.alignments[0].qs, 1);
        assert_eq!(record.alignments[0].qe, 600);
    }

    #[test]
    fn variant_segments_split_multi_alignment_fl_read_and_trim_overlap() {
        let line = "readC\tmito\t1200\taln_type=3;rep,rep,NA\t100.0\taln=1;len=500;olp=101;idt=100.0;strand=+;qs=1;qe=500;ss=1;se=500;cn=1;c1=100.0,1,1\taln=2;len=500;olp=51;idt=100.0;strand=+;qs=400;qe=899;ss=400;se=899;cn=1;c1=100.0,1,1\taln=3;len=352;olp=NA;idt=100.0;strand=+;qs=849;qe=1200;ss=849;se=1200;cn=1;c1=100.0,1,1";
        let record = parse_alignment_summary_record(line, 1).unwrap();
        let metadata = ReadIndexMetadata {
            read_class: "FL".to_string(),
            group_name: record.group_name(),
            subgroup_old_index: "3".to_string(),
            subgroup_key: record.subgroup_key().unwrap().label(),
        };
        let read = ReadRecord {
            id: "readC".to_string(),
            sequence: "A".repeat(1200),
        };
        let segments = variant_segments_for_record(
            &record,
            &metadata,
            &read.sequence,
            SnvIndelOverlapPolicy::AssignDownstream,
        )
        .unwrap();
        assert_eq!(segments.len(), 3);
        assert_eq!(segments[0].segment_id, "readC_P1");
        assert_eq!(segments[0].query_start, 1);
        assert_eq!(segments[0].query_end, 399);
        assert_eq!(segments[0].trim_note, "trimmed_overlap_before_P2");
        assert_eq!(segments[1].segment_id, "readC_P2");
        assert_eq!(segments[1].query_start, 400);
        assert_eq!(segments[1].query_end, 848);
        assert_eq!(segments[1].trim_note, "trimmed_overlap_before_P3");
        assert_eq!(segments[2].segment_id, "readC_P3");
        assert_eq!(segments[2].query_start, 849);
        assert_eq!(segments[2].query_end, 1200);
    }

    #[test]
    fn variant_segments_mark_overlap_keeps_both_sides_with_overlap_metadata() {
        let line = "readC\tmito\t1200\taln_type=3;rep,rep,NA\t100.0\taln=1;len=500;olp=101;idt=100.0;strand=+;qs=1;qe=500;ss=1;se=500;cn=1;c1=100.0,1,1\taln=2;len=500;olp=51;idt=100.0;strand=+;qs=400;qe=899;ss=400;se=899;cn=1;c1=100.0,1,1\taln=3;len=352;olp=NA;idt=100.0;strand=+;qs=849;qe=1200;ss=849;se=1200;cn=1;c1=100.0,1,1";
        let record = parse_alignment_summary_record(line, 1).unwrap();
        let metadata = ReadIndexMetadata {
            read_class: "FL".to_string(),
            group_name: record.group_name(),
            subgroup_old_index: "3".to_string(),
            subgroup_key: record.subgroup_key().unwrap().label(),
        };
        let read = ReadRecord {
            id: "readC".to_string(),
            sequence: "A".repeat(1200),
        };
        let segments = variant_segments_for_record(
            &record,
            &metadata,
            &read.sequence,
            SnvIndelOverlapPolicy::MarkOverlap,
        )
        .unwrap();
        assert_eq!(segments[0].query_start, 1);
        assert_eq!(segments[0].query_end, 500);
        assert_eq!(segments[1].query_start, 400);
        assert_eq!(segments[1].query_end, 899);
        assert_eq!(segments[2].query_start, 849);
        assert_eq!(segments[2].query_end, 1200);
        assert_eq!(
            segments[0].overlap_query_intervals,
            vec![QueryInterval {
                start: 400,
                end: 500
            }]
        );
        assert_eq!(
            segments[1].overlap_query_intervals,
            vec![
                QueryInterval {
                    start: 400,
                    end: 500
                },
                QueryInterval {
                    start: 849,
                    end: 899
                }
            ]
        );
    }

    #[test]
    fn variant_segments_mask_both_removes_overlap_metadata() {
        let line = "readC\tmito\t1200\taln_type=3;rep,rep,NA\t100.0\taln=1;len=500;olp=101;idt=100.0;strand=+;qs=1;qe=500;ss=1;se=500;cn=1;c1=100.0,1,1\taln=2;len=500;olp=51;idt=100.0;strand=+;qs=400;qe=899;ss=400;se=899;cn=1;c1=100.0,1,1\taln=3;len=352;olp=NA;idt=100.0;strand=+;qs=849;qe=1200;ss=849;se=1200;cn=1;c1=100.0,1,1";
        let record = parse_alignment_summary_record(line, 1).unwrap();
        let metadata = ReadIndexMetadata {
            read_class: "FL".to_string(),
            group_name: record.group_name(),
            subgroup_old_index: "3".to_string(),
            subgroup_key: record.subgroup_key().unwrap().label(),
        };
        let read = ReadRecord {
            id: "readC".to_string(),
            sequence: "A".repeat(1200),
        };
        let segments = variant_segments_for_record(
            &record,
            &metadata,
            &read.sequence,
            SnvIndelOverlapPolicy::MaskBoth,
        )
        .unwrap();
        assert_eq!(segments[0].query_start, 1);
        assert_eq!(segments[0].query_end, 399);
        assert_eq!(segments[1].query_start, 501);
        assert_eq!(segments[1].query_end, 848);
        assert_eq!(segments[2].query_start, 900);
        assert_eq!(segments[2].query_end, 1200);
        assert!(segments
            .iter()
            .all(|segment| segment.overlap_query_intervals.is_empty()));
    }

    #[test]
    fn variant_type_annotation_matches_snv_context() {
        let annotation = annotate_variant_type("ACGT", 1, "A", "G", "SNV");
        assert_eq!(annotation.trinucleotide_context, "TAC");
        assert_eq!(annotation.snv_group, "A>G");
        assert_eq!(annotation.indel_group, "-");
        assert_eq!(annotation.summary_anno, "-");
    }

    #[test]
    fn variant_frequency_depth_uses_all_fl_depth() {
        let path = std::env::temp_dir().join(format!(
            "orgraft_coverage_depth_test_{}.tsv",
            std::process::id()
        ));
        std::fs::write(
            &path,
            "position\tfl_depth\tpartial_depth\treference_support_depth\ttotal_depth\treference_support_fraction\n1\t10\t0\t4\t10\t0.4\n2\t20\t0\t5\t20\t0.25\n",
        )
        .unwrap();

        let depths = read_variant_frequency_depth(&path).unwrap();
        std::fs::remove_file(&path).unwrap();

        assert_eq!(depths[1], 10);
        assert_eq!(depths[2], 20);
    }

    #[test]
    fn local_indel_annotation_detects_mmej_deletion_context() {
        let reference = "TTGATAATGATCCCC";
        let compatible = annotate_variant_type(reference, 1, "TTGATAATGAT", "TTGAT", "InDel");
        let local = annotate_local_indel_context(
            reference,
            1,
            "TTGATAATGAT",
            "TTGAT",
            "InDel",
            &compatible,
        );

        assert_eq!(local.indel_group, "MMEJ");
        assert_eq!(local.method, "local_shift_deletion");
        assert!(local.microhomology_size.unwrap_or(0) > 1);
        assert!(local.summary_anno.contains("MH_size="));
    }

    #[test]
    fn variant_type_annotation_combines_homopolymer_indels() {
        let rows = vec![
            AnnotatedSingleVariant {
                pos: 18,
                ref_allele: "GTT".to_string(),
                alt_allele: "GTTT".to_string(),
                variant_type: "InDel".to_string(),
                id_list: "read1".to_string(),
                counts: 1,
                type_annotation: annotate_variant_type("GTTTT", 18, "GTT", "GTTT", "InDel"),
            },
            AnnotatedSingleVariant {
                pos: 18,
                ref_allele: "GTT".to_string(),
                alt_allele: "GT".to_string(),
                variant_type: "InDel".to_string(),
                id_list: "read2".to_string(),
                counts: 1,
                type_annotation: annotate_variant_type("GTTTT", 18, "GTT", "GT", "InDel"),
            },
        ];
        let combined = combine_annotated_variants(rows, &[0; 32]);
        assert_eq!(combined.len(), 1);
        assert_eq!(combined[0].row_type, "InDel,homopolymer");
        assert_eq!(combined[0].alt_allele, "GTTT#GT");
        assert_eq!(combined[0].combined_info, "poly-T;ref_size=2;1:1,-1:1");
        assert_eq!(combined[0].multi_allelic, "multi-allelic");
    }

    #[test]
    fn snv_indel_segment_sam_selection_prefers_sv_expected_interval() {
        let segment = VariantSegment {
            read_id: "readC".to_string(),
            segment_id: "readC_P1".to_string(),
            read_class: "FL".to_string(),
            group_name: "type_2_subtype_rep_NA".to_string(),
            subgroup_old_index: "1".to_string(),
            subgroup_key: "se1=1500,ss2=2000".to_string(),
            segment_index: 1,
            segment_count: 2,
            query_start: 1,
            query_end: 100,
            subject_start: 1000,
            subject_end: 1100,
            strand: '+',
            sequence: "A".repeat(100),
            trim_note: "none".to_string(),
            overlap_query_intervals: Vec::new(),
        };
        let far_high_mapq = SplitSamRecord {
            flag: 0,
            mapq: 60,
            reference_name: "ref".to_string(),
            position: 5000,
            cigar: "100M".to_string(),
            sequence: "A".repeat(100),
        };
        let expected_low_mapq = SplitSamRecord {
            flag: 0,
            mapq: 10,
            reference_name: "ref".to_string(),
            position: 1001,
            cigar: "100M".to_string(),
            sequence: "A".repeat(100),
        };
        let records = vec![far_high_mapq, expected_low_mapq];
        let selected = select_sam_records_for_segment(&segment, &records, true).unwrap();
        assert_eq!(selected.len(), 1);
        assert_eq!(selected[0].position, 1001);
        let unfiltered = select_sam_records_for_segment(&segment, &records, false).unwrap();
        assert_eq!(unfiltered.len(), 2);
    }

    #[test]
    fn alignment_summary_uses_full_length_threshold() {
        let line = "readB\tmito_adj_2_1\t1000\taln_type=1;NA\t97.999\taln=1;len=970;olp=NA;idt=100.0;strand=+;qs=1;qe=970;ss=1;se=970;cn=1;c1=100.0,1,1";
        let record = parse_alignment_summary_record(line, 1).unwrap();
        assert_eq!(record.group_name(), "type_1_subtype_NA");
        assert!(!record.is_fl());
    }

    #[test]
    fn missing_reads_fails() {
        let args = vec!["--organelle".to_string(), "mito".to_string()];
        assert!(matches!(
            PolishOptions::from_args(&args),
            Err(OrgraftError::InvalidArgument(_))
        ));
    }

    #[test]
    fn invalid_label_fails() {
        let args = vec![
            "--organelle".to_string(),
            "../mito".to_string(),
            "--reads".to_string(),
            "reads.fastq.gz".to_string(),
        ];
        assert!(matches!(
            PolishOptions::from_args(&args),
            Err(OrgraftError::InvalidArgument(_))
        ));
    }

    #[test]
    fn per_read_variant_calls_can_be_disabled() {
        let args = vec![
            "--organelle".to_string(),
            "mito".to_string(),
            "--reads".to_string(),
            "reads.fastq.gz".to_string(),
            "--per-read-variant-calls".to_string(),
            "off".to_string(),
        ];
        let options = PolishOptions::from_args(&args).unwrap();
        assert!(!options.per_read_variant_calls);
    }

    #[test]
    fn per_read_variant_calls_rejects_unknown_value() {
        let args = vec![
            "--organelle".to_string(),
            "mito".to_string(),
            "--reads".to_string(),
            "reads.fastq.gz".to_string(),
            "--per-read-variant-calls".to_string(),
            "maybe".to_string(),
        ];
        assert!(matches!(
            PolishOptions::from_args(&args),
            Err(OrgraftError::InvalidArgument(_))
        ));
    }

    #[test]
    fn snv_indel_overlap_policy_can_use_mask_both() {
        let args = vec![
            "--organelle".to_string(),
            "mito".to_string(),
            "--reads".to_string(),
            "reads.fastq.gz".to_string(),
            "--snv-indel-overlap-policy".to_string(),
            "mask-both".to_string(),
        ];
        let options = PolishOptions::from_args(&args).unwrap();
        assert_eq!(
            options.snv_indel_overlap_policy,
            SnvIndelOverlapPolicy::MaskBoth
        );
    }

    #[test]
    fn plot_range_and_sv_highlights_are_parsed() {
        let args = vec![
            "--organelle".to_string(),
            "mito".to_string(),
            "--reads".to_string(),
            "reads.fastq.gz".to_string(),
            "--plot-range".to_string(),
            "100-5000".to_string(),
            "--plot-dpi".to_string(),
            "600".to_string(),
            "--plot-output-format".to_string(),
            "both".to_string(),
            "--coverage-plot-rasterize".to_string(),
            "off".to_string(),
            "--snv-indel-plot-rasterize".to_string(),
            "off".to_string(),
            "--sv-plot-highlight-subgroups".to_string(),
            "type_3_subtype_rep_rep_NA:3,type_3_subtype_rep_rep_NA:5".to_string(),
            "--sv-plot-highlight-read-ids".to_string(),
            "ids.txt".to_string(),
            "--sv-plot-highlight-min-fraction".to_string(),
            "0.01".to_string(),
            "--sv-plot-highlight-min-reads".to_string(),
            "25".to_string(),
        ];
        let options = PolishOptions::from_args(&args).unwrap();
        assert_eq!(options.plot_range.as_ref().unwrap().as_arg(), "100-5000");
        assert_eq!(options.plot_dpi, 600);
        assert_eq!(options.plot_output_format, PlotOutputFormat::Both);
        assert!(!options.coverage_plot_rasterize);
        assert!(!options.snv_indel_plot_rasterize);
        assert_eq!(
            options.sv_plot_highlight_subgroups,
            vec![
                "type_3_subtype_rep_rep_NA:3".to_string(),
                "type_3_subtype_rep_rep_NA:5".to_string()
            ]
        );
        assert_eq!(
            options.sv_plot_highlight_read_ids,
            Some(PathBuf::from("ids.txt"))
        );
        assert_eq!(options.sv_plot_highlight_min_fraction, 0.01);
        assert_eq!(options.sv_plot_highlight_min_reads, 25);
    }

    #[test]
    fn highlight_fraction_rejects_out_of_range_value() {
        let args = vec![
            "--organelle".to_string(),
            "mito".to_string(),
            "--reads".to_string(),
            "reads.fastq.gz".to_string(),
            "--sv-plot-highlight-min-fraction".to_string(),
            "1.2".to_string(),
        ];
        let err = PolishOptions::from_args(&args).unwrap_err();
        assert!(format!("{err}").contains("expected fraction between 0 and 1"));
    }

    #[test]
    fn sv_plot_highlight_read_id_file_alias_is_parsed() {
        let args = vec![
            "--organelle".to_string(),
            "mito".to_string(),
            "--reads".to_string(),
            "reads.fastq.gz".to_string(),
            "--sv-plot-highlight-read_id_file".to_string(),
            "ids.txt".to_string(),
        ];
        let options = PolishOptions::from_args(&args).unwrap();
        assert_eq!(
            options.sv_plot_highlight_read_ids,
            Some(PathBuf::from("ids.txt"))
        );
    }

    #[test]
    fn removed_polish_options_are_rejected() {
        for option in [
            "--sv-minimap2-mode",
            "--snv-indel-caller",
            "--legacy-sv-reports",
            "--id-map",
            "--dry-run",
            "--auto-plot",
            "--snv-indel-sv-context-filter",
            "--plot-rasterize",
            "--coverage-rasterize",
            "--highlight-min-fraction",
        ] {
            let args = vec![
                "--organelle".to_string(),
                "mito".to_string(),
                "--reads".to_string(),
                "reads.fastq.gz".to_string(),
                option.to_string(),
                "x".to_string(),
            ];
            assert!(matches!(
                PolishOptions::from_args(&args),
                Err(OrgraftError::InvalidArgument(_))
            ));
        }
    }
}
