use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet};
use std::fs::{self, File};
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};
use std::process::{Child, Command, Stdio};

use crate::commands::shared::{print_contract, CommandContract};
use crate::error::OrgraftError;

const HELP: &str = r#"orgraft recruit

Reference/seed-based read recruitment prototype.

Usage:
  orgraft recruit --reads reads.fastq[.gz] --mito mito.fa --plastid plastid.fa [options]

Inputs:
  --reads FILE              input reads; fastq, fq, fastq.gz, fq.gz only
  --mito FILE               mitochondrial FASTA bait
  --plastid FILE            plastid FASTA bait

Outputs:
  --out-dir DIR             output directory [results/recruit]
  --gzip-output MODE        on|off; write recruited FASTQ files as .fastq.gz [on]

Additional Parameters:
  --threads N               minimap2 and pigz threads [1]
  --platform NAME           HiFi, CLR, ONT, ultra-long [HiFi]
  --bait-format FORMAT      auto, fasta, gfa [auto]
  --max-reads NAME,N        cap selected reads for a bait label or partition
  --advanced-help           show bait, sampling, debug, and compression options

Layout: OUT/{*.fastq[.gz],logs}
"#;

const ADVANCED_HELP: &str = r#"orgraft recruit --advanced-help

Advanced recruit options.

Name/ID options:
  --bait LABEL=FILE             generic FASTA/GFA bait; repeatable
  --prefix NAME                 fallback label/prefix for unlabeled bait
  --rename-bait                 rewrite bait IDs as LABEL_1, LABEL_2, ...

Split options:
  --gfa-split MODE              all, components [all]
  --split-output MODE           label, partition, none [label]

Write options:
  --write-id-map                write bait_id_map.tsv
  --write-read-classification   write read_classification.tsv
  --write-sampled-ids           write sampled_read_ids.tsv
  --write-bait-partitions       write bait partition FASTA files

Other options:
  --random-seed N               deterministic sampling seed [42]
  --read-stats MODE             basic, full [basic]

"#;

pub fn run(args: &[String]) -> Result<(), OrgraftError> {
    if args.is_empty() || args.iter().any(|arg| arg == "-h" || arg == "--help") {
        println!("{HELP}");
        return Ok(());
    }

    if args.iter().any(|arg| arg == "--advanced-help") {
        println!("{ADVANCED_HELP}");
        return Ok(());
    }

    if args.iter().any(|arg| arg == "--contract") {
        print_contract(&contract());
        return Ok(());
    }

    let options = RecruitOptions::from_args(args)?;
    run_recruitment(&options)
}

fn run_recruitment(options: &RecruitOptions) -> Result<(), OrgraftError> {
    if options.mode == RecruitmentMode::Iterative || options.iterations != 1 {
        return Err(OrgraftError::InvalidArgument(
            "iterative recruitment is reserved in the CLI, but this prototype only runs one reference/seed-based round".to_string(),
        ));
    }

    validate_reads_path(&options.reads)?;
    fs::create_dir_all(&options.out_dir)?;
    fs::create_dir_all(logs_dir(options))?;

    let gzip = GzipRuntime::resolve(&options.gzip, options.threads)?;
    let bait = prepare_baits(options, &gzip)?;
    let alignment = match &options.sam {
        Some(path) => parse_sam(path, &bait.targets, options, &gzip)?,
        None => run_minimap2(options, &bait.fasta_path, &bait.targets, &gzip)?,
    };

    let selection = select_reads(&alignment.reads, &bait, options)?;
    if options.write_read_classification {
        write_read_classification(
            &options.out_dir.join("read_classification.tsv"),
            &alignment.reads,
            &selection.selected_all,
        )?;
    }

    let copy_stats = extract_fastq(options, &selection, &gzip)?;
    if options.write_sampled_ids {
        write_sampled_ids(&options.out_dir.join("sampled_read_ids.tsv"), &selection)?;
    }

    let summary = SummaryContext {
        options,
        gzip: &gzip,
        bait: &bait,
        alignment: &alignment,
        selection: &selection,
        copy_stats: &copy_stats,
    };
    write_summary(&summary_path(options), &summary)?;
    write_read_stats(
        &read_stats_path(options),
        &copy_stats.read_stats,
        &options.read_stats_mode,
    )?;

    for output in copy_stats.split_outputs.values() {
        println!("Wrote {}", output.display());
    }
    println!("Wrote {}", summary_path(options).display());
    println!("Wrote {}", read_stats_path(options).display());
    Ok(())
}

fn contract() -> CommandContract {
    CommandContract {
        command: "recruit",
        origin: "OrgRAFT read recruitment logic",
        purpose: "select organelle-enriched HiFi/CCS reads before graph assembly and validation",
        inputs: &[
            "raw HiFi/CCS reads in FASTQ/FASTQ.GZ",
            "reference/seed bait as FASTA or GFA",
            "optional precomputed SAM",
        ],
        outputs: &[
            "<label>.fastq[.gz]",
            "logs/recruitment_summary.tsv",
            "logs/read_stats.tsv",
        ],
        notes: &[
            "reference/seed-based recruitment is runnable through minimap2",
            "iterative recruitment flags are reserved but not implemented yet",
            "ambiguous reads are reported rather than silently dropped",
        ],
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
enum RecruitmentMode {
    Reference,
    Iterative,
}

#[derive(Debug, Clone, PartialEq, Eq)]
enum BaitFormat {
    Auto,
    Fasta,
    Gfa,
}

#[derive(Debug, Clone, PartialEq, Eq)]
enum GfaSplit {
    All,
    Components,
}

#[derive(Debug, Clone, PartialEq, Eq)]
enum SplitOutput {
    None,
    Label,
    Partition,
}

#[derive(Debug, Clone, PartialEq, Eq)]
enum GzipToolChoice {
    Auto,
    Pigz,
    Gzip,
}

#[derive(Debug, Clone, PartialEq, Eq)]
enum AlignMode {
    Sam,
    PafCigar,
}

#[derive(Debug, Clone, PartialEq, Eq)]
enum ReadStatsMode {
    Basic,
    Full,
}

#[derive(Debug, Clone)]
struct GzipConfig {
    choice: GzipToolChoice,
}

#[derive(Debug, Clone)]
struct BaitInput {
    label: Option<String>,
    path: PathBuf,
}

#[derive(Debug, Clone)]
struct RecruitOptions {
    reads: PathBuf,
    baits: Vec<BaitInput>,
    out_dir: PathBuf,
    prefix: Option<String>,
    bait_format: BaitFormat,
    gfa_split: GfaSplit,
    rename_bait: bool,
    write_id_map: bool,
    split_output: SplitOutput,
    gzip_output: bool,
    minimap2: String,
    align_mode: AlignMode,
    preset: String,
    min_mapq: u8,
    min_aln_len: u64,
    sam: Option<PathBuf>,
    max_reads: BTreeMap<String, usize>,
    random_seed: u64,
    write_sampled_ids: bool,
    read_stats_mode: ReadStatsMode,
    write_read_classification: bool,
    write_bait_partitions: bool,
    gzip: GzipConfig,
    mode: RecruitmentMode,
    iterations: usize,
    threads: usize,
}

impl RecruitOptions {
    fn from_args(args: &[String]) -> Result<Self, OrgraftError> {
        let mut options = Self {
            reads: PathBuf::new(),
            baits: Vec::new(),
            out_dir: PathBuf::from("results/recruit"),
            prefix: None,
            bait_format: BaitFormat::Auto,
            gfa_split: GfaSplit::All,
            rename_bait: false,
            write_id_map: false,
            split_output: SplitOutput::Label,
            gzip_output: true,
            minimap2: "minimap2".to_string(),
            align_mode: AlignMode::PafCigar,
            preset: "map-hifi".to_string(),
            min_mapq: 0,
            min_aln_len: 0,
            sam: None,
            max_reads: BTreeMap::new(),
            random_seed: 42,
            write_sampled_ids: false,
            read_stats_mode: ReadStatsMode::Basic,
            write_read_classification: false,
            write_bait_partitions: false,
            gzip: GzipConfig {
                choice: GzipToolChoice::Auto,
            },
            mode: RecruitmentMode::Reference,
            iterations: 1,
            threads: 1,
        };

        let mut index = 0;
        while index < args.len() {
            match args[index].as_str() {
                "--reads" => {
                    options.reads = PathBuf::from(value_after(args, &mut index, "--reads")?);
                }
                "--mito" => {
                    options.baits.push(BaitInput {
                        label: Some("mito".to_string()),
                        path: PathBuf::from(value_after(args, &mut index, "--mito")?),
                    });
                }
                "--plastid" | "--plasti" => {
                    let name = args[index].clone();
                    options.baits.push(BaitInput {
                        label: Some("plastid".to_string()),
                        path: PathBuf::from(value_after(args, &mut index, &name)?),
                    });
                }
                "--bait" => {
                    options
                        .baits
                        .push(parse_bait_input(value_after(args, &mut index, "--bait")?));
                }
                "--out-dir" => {
                    options.out_dir = PathBuf::from(value_after(args, &mut index, "--out-dir")?);
                }
                "--prefix" => {
                    options.prefix = Some(value_after(args, &mut index, "--prefix")?.to_string());
                }
                "--bait-format" => {
                    options.bait_format =
                        parse_bait_format(value_after(args, &mut index, "--bait-format")?)?;
                }
                "--gfa-split" => {
                    options.gfa_split =
                        parse_gfa_split(value_after(args, &mut index, "--gfa-split")?)?;
                }
                "--rename-bait" => options.rename_bait = true,
                "--write-id-map" => options.write_id_map = true,
                "--split-output" => {
                    options.split_output =
                        parse_split_output(value_after(args, &mut index, "--split-output")?)?;
                }
                "--gzip-output" => {
                    options.gzip_output =
                        parse_on_off(value_after(args, &mut index, "--gzip-output")?)?;
                }
                "--minimap2" => {
                    options.minimap2 = value_after(args, &mut index, "--minimap2")?.to_string();
                }
                "--align-mode" => {
                    options.align_mode =
                        parse_align_mode(value_after(args, &mut index, "--align-mode")?)?;
                }
                "--threads" => {
                    options.threads = parse_number(value_after(args, &mut index, "--threads")?)?;
                    if options.threads == 0 {
                        return Err(OrgraftError::InvalidArgument(
                            "--threads must be greater than 0".to_string(),
                        ));
                    }
                }
                "--platform" => {
                    options.preset =
                        platform_to_minimap2_preset(value_after(args, &mut index, "--platform")?)?
                            .to_string();
                }
                "--preset" => {
                    // Hidden compatibility escape hatch; the public CLI keeps
                    // platform as the single way to choose the minimap2 preset.
                    options.preset = value_after(args, &mut index, "--preset")?.to_string();
                }
                "--min-mapq" => {
                    // Hidden compatibility/debug option. The supported
                    // prototype behavior keeps both recruitment filters at 0.
                    options.min_mapq = parse_number(value_after(args, &mut index, "--min-mapq")?)?;
                }
                "--min-aln-len" => {
                    // Hidden compatibility/debug option. The supported
                    // prototype behavior keeps both recruitment filters at 0.
                    options.min_aln_len =
                        parse_number(value_after(args, &mut index, "--min-aln-len")?)?;
                }
                "--sam" => {
                    options.sam = Some(PathBuf::from(value_after(args, &mut index, "--sam")?));
                }
                "--max-reads" => {
                    let (name, max) =
                        parse_max_reads(value_after(args, &mut index, "--max-reads")?)?;
                    options.max_reads.insert(name, max);
                }
                "--random-seed" => {
                    options.random_seed =
                        parse_number(value_after(args, &mut index, "--random-seed")?)?;
                }
                "--write-sampled-ids" => options.write_sampled_ids = true,
                "--read-stats" => {
                    options.read_stats_mode =
                        parse_read_stats_mode(value_after(args, &mut index, "--read-stats")?)?;
                }
                "--write-read-classification" => options.write_read_classification = true,
                "--write-bait-partitions" => options.write_bait_partitions = true,
                "--gzip-tool" => {
                    options.gzip.choice =
                        parse_gzip_tool(value_after(args, &mut index, "--gzip-tool")?)?;
                }
                "--gzip-threads" => {
                    let _ = value_after(args, &mut index, "--gzip-threads")?;
                    return Err(OrgraftError::InvalidArgument(
                        "--gzip-threads has been removed; use --threads to control minimap2 and pigz threads".to_string(),
                    ));
                }
                "--mode" => {
                    options.mode = parse_mode(value_after(args, &mut index, "--mode")?)?;
                }
                "--iterations" | "--rounds" => {
                    let name = args[index].clone();
                    options.iterations = parse_number(value_after(args, &mut index, &name)?)?;
                }
                other => {
                    return Err(OrgraftError::InvalidArgument(format!(
                        "unknown recruit option `{other}`"
                    )));
                }
            }
            index += 1;
        }

        if options.reads.as_os_str().is_empty() {
            return Err(OrgraftError::InvalidArgument(
                "missing required --reads".to_string(),
            ));
        }
        if options.baits.is_empty() {
            return Err(OrgraftError::InvalidArgument(
                "missing required --bait".to_string(),
            ));
        }

        Ok(options)
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

fn parse_bait_input(value: &str) -> BaitInput {
    if let Some((label, path)) = value.split_once('=') {
        if is_label_like(label) && !path.is_empty() {
            return BaitInput {
                label: Some(label.to_string()),
                path: PathBuf::from(path),
            };
        }
    }

    BaitInput {
        label: None,
        path: PathBuf::from(value),
    }
}

fn is_label_like(value: &str) -> bool {
    !value.is_empty()
        && value
            .chars()
            .all(|character| character.is_ascii_alphanumeric() || matches!(character, '_' | '-'))
}

fn parse_bait_format(value: &str) -> Result<BaitFormat, OrgraftError> {
    match value {
        "auto" => Ok(BaitFormat::Auto),
        "fasta" | "fa" => Ok(BaitFormat::Fasta),
        "gfa" => Ok(BaitFormat::Gfa),
        other => Err(OrgraftError::InvalidArgument(format!(
            "unknown bait format `{other}`; expected auto, fasta, or gfa"
        ))),
    }
}

fn parse_gfa_split(value: &str) -> Result<GfaSplit, OrgraftError> {
    match value {
        "all" => Ok(GfaSplit::All),
        "components" | "component" | "connected" => Ok(GfaSplit::Components),
        other => Err(OrgraftError::InvalidArgument(format!(
            "unknown GFA split mode `{other}`; expected all or components"
        ))),
    }
}

fn parse_split_output(value: &str) -> Result<SplitOutput, OrgraftError> {
    match value {
        "none" | "all" => Ok(SplitOutput::None),
        "label" | "labels" | "organelle" => Ok(SplitOutput::Label),
        "partition" | "partitions" | "bait" => Ok(SplitOutput::Partition),
        other => Err(OrgraftError::InvalidArgument(format!(
            "unknown split-output mode `{other}`; expected none, label, or partition"
        ))),
    }
}

fn parse_gzip_tool(value: &str) -> Result<GzipToolChoice, OrgraftError> {
    match value {
        "auto" => Ok(GzipToolChoice::Auto),
        "pigz" => Ok(GzipToolChoice::Pigz),
        "gzip" => Ok(GzipToolChoice::Gzip),
        other => Err(OrgraftError::InvalidArgument(format!(
            "unknown gzip tool `{other}`; expected auto, pigz, or gzip"
        ))),
    }
}

fn parse_on_off(value: &str) -> Result<bool, OrgraftError> {
    match value {
        "on" => Ok(true),
        "off" => Ok(false),
        other => Err(OrgraftError::InvalidArgument(format!(
            "unknown on/off value `{other}`; expected on or off"
        ))),
    }
}

fn parse_align_mode(value: &str) -> Result<AlignMode, OrgraftError> {
    match value {
        "sam" => Ok(AlignMode::Sam),
        "paf-cigar" | "paf_cigar" | "pafc" | "paf-c" => Ok(AlignMode::PafCigar),
        other => Err(OrgraftError::InvalidArgument(format!(
            "unknown align mode `{other}`; expected sam or paf-cigar"
        ))),
    }
}

fn parse_read_stats_mode(value: &str) -> Result<ReadStatsMode, OrgraftError> {
    match value {
        "basic" => Ok(ReadStatsMode::Basic),
        "full" => Ok(ReadStatsMode::Full),
        other => Err(OrgraftError::InvalidArgument(format!(
            "unknown read-stats mode `{other}`; expected basic or full"
        ))),
    }
}

fn parse_mode(value: &str) -> Result<RecruitmentMode, OrgraftError> {
    match value {
        "reference" | "seed" | "seed-based" => Ok(RecruitmentMode::Reference),
        "iterative" => Ok(RecruitmentMode::Iterative),
        other => Err(OrgraftError::InvalidArgument(format!(
            "unknown recruitment mode `{other}`; expected reference or iterative"
        ))),
    }
}

fn platform_to_minimap2_preset(value: &str) -> Result<&'static str, OrgraftError> {
    match value {
        "HiFi" | "hifi" | "CCS" | "ccs" => Ok("map-hifi"),
        "CLR" | "clr" | "PB" | "pb" => Ok("map-pb"),
        "ONT" | "ont" | "nanopore" | "ultra-long" | "ultralong" => Ok("map-ont"),
        other => Err(OrgraftError::InvalidArgument(format!(
            "unknown read platform `{other}`; expected HiFi, CLR, ONT, or ultra-long"
        ))),
    }
}

fn parse_max_reads(value: &str) -> Result<(String, usize), OrgraftError> {
    let (name, count) = value
        .split_once(',')
        .ok_or_else(|| OrgraftError::InvalidArgument("--max-reads expects NAME,N".to_string()))?;
    let name = sanitize_name(name);
    if name.is_empty() {
        return Err(OrgraftError::InvalidArgument(
            "--max-reads NAME must not be empty".to_string(),
        ));
    }
    if name == "all" {
        return Err(OrgraftError::InvalidArgument(
            "--max-reads all,N is not supported; use a bait label such as mito,N or plastid,N"
                .to_string(),
        ));
    }
    let count = parse_number(count)?;
    Ok((name, count))
}

fn parse_number<T>(value: &str) -> Result<T, OrgraftError>
where
    T: std::str::FromStr,
{
    value
        .parse::<T>()
        .map_err(|_| OrgraftError::InvalidArgument(format!("expected a number, got `{value}`")))
}

fn validate_reads_path(path: &Path) -> Result<(), OrgraftError> {
    if !path.exists() {
        return Err(OrgraftError::InvalidArgument(format!(
            "reads file not found: {}",
            path.display()
        )));
    }
    if !is_fastq_path(path) {
        return Err(OrgraftError::InvalidArgument(format!(
            "unsupported reads format {}; expected fastq, fq, fastq.gz, or fq.gz",
            path.display()
        )));
    }
    Ok(())
}

fn logs_dir(options: &RecruitOptions) -> PathBuf {
    options.out_dir.join("logs")
}

fn summary_path(options: &RecruitOptions) -> PathBuf {
    logs_dir(options).join("recruitment_summary.tsv")
}

fn read_stats_path(options: &RecruitOptions) -> PathBuf {
    logs_dir(options).join("read_stats.tsv")
}

fn is_fastq_path(path: &Path) -> bool {
    let name = path
        .file_name()
        .and_then(|value| value.to_str())
        .unwrap_or("")
        .to_ascii_lowercase();
    matches!(
        name.as_str(),
        name if name.ends_with(".fastq")
            || name.ends_with(".fq")
            || name.ends_with(".fastq.gz")
            || name.ends_with(".fq.gz")
    )
}

#[derive(Debug, Clone)]
struct FastaRecord {
    id: String,
    seq: String,
}

#[derive(Debug, Clone)]
struct TargetInfo {
    label: String,
    partition: String,
    original_id: String,
}

#[derive(Debug)]
struct BaitPreparation {
    fasta_path: PathBuf,
    id_map_path: Option<PathBuf>,
    targets: HashMap<String, TargetInfo>,
    labels: BTreeSet<String>,
    partitions: BTreeMap<String, String>,
    record_count: usize,
}

fn prepare_baits(
    options: &RecruitOptions,
    gzip: &GzipRuntime,
) -> Result<BaitPreparation, OrgraftError> {
    let mut rows = Vec::new();
    let mut labels = BTreeSet::new();
    let mut partitions = BTreeMap::new();
    let mut used_ids = HashSet::new();
    let mut label_counts: BTreeMap<String, usize> = BTreeMap::new();

    for (bait_index, input) in options.baits.iter().enumerate() {
        if !input.path.exists() {
            return Err(OrgraftError::InvalidArgument(format!(
                "bait file not found: {}",
                input.path.display()
            )));
        }

        let label = bait_label(input, options, bait_index);
        labels.insert(label.clone());
        let format = resolve_bait_format(&input.path, &options.bait_format)?;

        match format {
            BaitFormat::Fasta => {
                let text = read_text(&input.path, gzip)?;
                let records = parse_fasta(&text, &input.path)?;
                let partition = label.clone();
                partitions.insert(partition.clone(), label.clone());
                for record in records {
                    let count = next_label_count(&mut label_counts, &label);
                    let new_id = bait_record_id(
                        options.rename_bait,
                        &label,
                        count,
                        &record.id,
                        &mut used_ids,
                    );
                    rows.push(PreparedBaitRecord {
                        new_id,
                        original_id: record.id,
                        label: label.clone(),
                        partition: partition.clone(),
                        seq: record.seq,
                        source: input.path.clone(),
                    });
                }
            }
            BaitFormat::Gfa => {
                let text = read_text(&input.path, gzip)?;
                let graph = parse_gfa(&text, &input.path)?;
                let components = match options.gfa_split {
                    GfaSplit::All => vec![graph
                        .segments
                        .iter()
                        .map(|segment| segment.id.clone())
                        .collect()],
                    GfaSplit::Components => graph.connected_components(),
                };
                for (component_index, component) in components.iter().enumerate() {
                    let partition = match options.gfa_split {
                        GfaSplit::All => label.clone(),
                        GfaSplit::Components => format!("{}_{}", label, component_index + 1),
                    };
                    partitions.insert(partition.clone(), label.clone());
                    for segment_id in component {
                        let Some(segment) = graph.segment(segment_id) else {
                            continue;
                        };
                        let count = next_label_count(&mut label_counts, &label);
                        let new_id = bait_record_id(
                            options.rename_bait,
                            &label,
                            count,
                            &segment.id,
                            &mut used_ids,
                        );
                        rows.push(PreparedBaitRecord {
                            new_id,
                            original_id: segment.id.clone(),
                            label: label.clone(),
                            partition: partition.clone(),
                            seq: segment.seq.clone(),
                            source: input.path.clone(),
                        });
                    }
                }
            }
            BaitFormat::Auto => unreachable!("auto format should be resolved"),
        }
    }

    if rows.is_empty() {
        return Err(OrgraftError::InvalidArgument(
            "bait did not contain any usable FASTA/GFA sequences".to_string(),
        ));
    }

    let fasta_path = logs_dir(options).join("bait.fasta");
    write_prepared_fasta(&fasta_path, &rows)?;
    if options.write_bait_partitions {
        write_partition_fastas(&options.out_dir.join("bait_partitions"), &rows)?;
    }

    let id_map_path = if options.write_id_map {
        let path = options.out_dir.join("bait_id_map.tsv");
        write_id_map(&path, &rows)?;
        Some(path)
    } else {
        None
    };

    let mut targets = HashMap::new();
    for row in &rows {
        targets.insert(
            row.new_id.clone(),
            TargetInfo {
                label: row.label.clone(),
                partition: row.partition.clone(),
                original_id: row.original_id.clone(),
            },
        );
    }

    Ok(BaitPreparation {
        fasta_path,
        id_map_path,
        targets,
        labels,
        partitions,
        record_count: rows.len(),
    })
}

#[derive(Debug, Clone)]
struct PreparedBaitRecord {
    new_id: String,
    original_id: String,
    label: String,
    partition: String,
    seq: String,
    source: PathBuf,
}

fn bait_label(input: &BaitInput, options: &RecruitOptions, index: usize) -> String {
    if let Some(label) = &input.label {
        return sanitize_name(label);
    }
    if let Some(prefix) = &options.prefix {
        return sanitize_name(prefix);
    }
    let stem = file_stem_without_gz(&input.path).unwrap_or_else(|| format!("bait{}", index + 1));
    sanitize_name(&stem)
}

fn resolve_bait_format(path: &Path, requested: &BaitFormat) -> Result<BaitFormat, OrgraftError> {
    if *requested != BaitFormat::Auto {
        return Ok(requested.clone());
    }

    let name = path
        .file_name()
        .and_then(|value| value.to_str())
        .unwrap_or("")
        .to_ascii_lowercase();
    let name = name.strip_suffix(".gz").unwrap_or(&name);
    if name.ends_with(".gfa") {
        Ok(BaitFormat::Gfa)
    } else if name.ends_with(".fa")
        || name.ends_with(".fasta")
        || name.ends_with(".fna")
        || name.ends_with(".fas")
    {
        Ok(BaitFormat::Fasta)
    } else {
        Err(OrgraftError::InvalidArgument(format!(
            "cannot infer bait format from {}; use --bait-format fasta or gfa",
            path.display()
        )))
    }
}

fn file_stem_without_gz(path: &Path) -> Option<String> {
    let name = path.file_name()?.to_str()?;
    let without_gz = name.strip_suffix(".gz").unwrap_or(name);
    let stem = without_gz
        .rsplit_once('.')
        .map(|(stem, _)| stem)
        .unwrap_or(without_gz);
    Some(stem.to_string())
}

fn next_label_count(counts: &mut BTreeMap<String, usize>, label: &str) -> usize {
    let counter = counts.entry(label.to_string()).or_insert(0);
    *counter += 1;
    *counter
}

fn bait_record_id(
    rename: bool,
    label: &str,
    count: usize,
    original_id: &str,
    used_ids: &mut HashSet<String>,
) -> String {
    let base = if rename {
        format!("{label}_{count}")
    } else {
        sanitize_name(original_id)
    };
    make_unique_id(&base, used_ids)
}

fn make_unique_id(base: &str, used_ids: &mut HashSet<String>) -> String {
    let base = if base.is_empty() { "seq" } else { base };
    if used_ids.insert(base.to_string()) {
        return base.to_string();
    }

    let mut index = 2;
    loop {
        let candidate = format!("{base}_{index}");
        if used_ids.insert(candidate.clone()) {
            return candidate;
        }
        index += 1;
    }
}

fn sanitize_name(value: &str) -> String {
    let mut result = String::new();
    for character in value.trim().chars() {
        if character.is_ascii_alphanumeric() || matches!(character, '_' | '-' | '.') {
            result.push(character);
        } else {
            result.push('_');
        }
    }
    result.trim_matches('_').to_string()
}

fn parse_fasta(text: &str, path: &Path) -> Result<Vec<FastaRecord>, OrgraftError> {
    let mut records = Vec::new();
    let mut current_id: Option<String> = None;
    let mut current_seq = String::new();

    for (index, line) in text.lines().enumerate() {
        if let Some(header) = line.strip_prefix('>') {
            if let Some(id) = current_id.take() {
                push_fasta_record(&mut records, id, &mut current_seq, path)?;
            }
            let id = header.split_whitespace().next().unwrap_or("").trim();
            if id.is_empty() {
                return Err(OrgraftError::InvalidArgument(format!(
                    "{}:{} FASTA header has no ID",
                    path.display(),
                    index + 1
                )));
            }
            current_id = Some(id.to_string());
        } else if !line.trim().is_empty() {
            if current_id.is_none() {
                return Err(OrgraftError::InvalidArgument(format!(
                    "{}:{} FASTA sequence appears before first header",
                    path.display(),
                    index + 1
                )));
            }
            current_seq.push_str(line.trim());
        }
    }

    if let Some(id) = current_id {
        push_fasta_record(&mut records, id, &mut current_seq, path)?;
    }

    if records.is_empty() {
        return Err(OrgraftError::InvalidArgument(format!(
            "{} contains no FASTA records",
            path.display()
        )));
    }
    Ok(records)
}

fn push_fasta_record(
    records: &mut Vec<FastaRecord>,
    id: String,
    seq: &mut String,
    path: &Path,
) -> Result<(), OrgraftError> {
    if seq.is_empty() {
        return Err(OrgraftError::InvalidArgument(format!(
            "{} FASTA record `{id}` has no sequence",
            path.display()
        )));
    }
    records.push(FastaRecord {
        id,
        seq: std::mem::take(seq),
    });
    Ok(())
}

#[derive(Debug, Clone)]
struct GfaSegment {
    id: String,
    seq: String,
}

#[derive(Debug, Clone)]
struct GfaGraph {
    segments: Vec<GfaSegment>,
    links: Vec<(String, String)>,
}

impl GfaGraph {
    fn segment(&self, id: &str) -> Option<&GfaSegment> {
        self.segments.iter().find(|segment| segment.id == id)
    }

    fn connected_components(&self) -> Vec<Vec<String>> {
        let mut index_by_id = HashMap::new();
        for (index, segment) in self.segments.iter().enumerate() {
            index_by_id.insert(segment.id.clone(), index);
        }

        let mut dsu = Dsu::new(self.segments.len());
        for (left, right) in &self.links {
            if let (Some(left_index), Some(right_index)) =
                (index_by_id.get(left), index_by_id.get(right))
            {
                dsu.union(*left_index, *right_index);
            }
        }

        let mut grouped: BTreeMap<usize, Vec<String>> = BTreeMap::new();
        for (index, segment) in self.segments.iter().enumerate() {
            grouped
                .entry(dsu.find(index))
                .or_default()
                .push(segment.id.clone());
        }

        let mut components: Vec<Vec<String>> = grouped
            .into_values()
            .map(|mut component| {
                component.sort();
                component
            })
            .collect();
        components.sort_by(|left, right| left.first().cmp(&right.first()));
        components
    }
}

fn parse_gfa(text: &str, path: &Path) -> Result<GfaGraph, OrgraftError> {
    let mut segments = Vec::new();
    let mut links = Vec::new();

    for (index, line) in text.lines().enumerate() {
        if line.trim().is_empty() || line.starts_with('#') {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        match fields.first().copied() {
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
                        "{}:{} GFA segment `{}` has `*` sequence and cannot be used as bait",
                        path.display(),
                        index + 1,
                        fields[1]
                    )));
                }
                segments.push(GfaSegment {
                    id: fields[1].to_string(),
                    seq: fields[2].to_string(),
                });
            }
            Some("L") => {
                if fields.len() >= 4 {
                    links.push((fields[1].to_string(), fields[3].to_string()));
                }
            }
            Some("J") => {
                if fields.len() >= 4 {
                    links.push((fields[1].to_string(), fields[3].to_string()));
                }
            }
            _ => {}
        }
    }

    if segments.is_empty() {
        return Err(OrgraftError::InvalidArgument(format!(
            "{} contains no GFA segment sequences",
            path.display()
        )));
    }

    Ok(GfaGraph { segments, links })
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

fn write_prepared_fasta(path: &Path, rows: &[PreparedBaitRecord]) -> Result<(), OrgraftError> {
    let mut writer = BufWriter::new(File::create(path)?);
    for row in rows {
        writeln!(writer, ">{}", row.new_id)?;
        write_wrapped_sequence(&mut writer, &row.seq)?;
    }
    Ok(())
}

fn write_partition_fastas(path: &Path, rows: &[PreparedBaitRecord]) -> Result<(), OrgraftError> {
    fs::create_dir_all(path)?;
    let mut by_partition: BTreeMap<&str, Vec<&PreparedBaitRecord>> = BTreeMap::new();
    for row in rows {
        by_partition.entry(&row.partition).or_default().push(row);
    }
    for (partition, partition_rows) in by_partition {
        let file_path = path.join(format!("{partition}.fasta"));
        let mut writer = BufWriter::new(File::create(file_path)?);
        for row in partition_rows {
            writeln!(writer, ">{}", row.new_id)?;
            write_wrapped_sequence(&mut writer, &row.seq)?;
        }
    }
    Ok(())
}

fn write_wrapped_sequence(writer: &mut dyn Write, seq: &str) -> Result<(), OrgraftError> {
    for chunk in seq.as_bytes().chunks(80) {
        writer.write_all(chunk)?;
        writer.write_all(b"\n")?;
    }
    Ok(())
}

fn write_id_map(path: &Path, rows: &[PreparedBaitRecord]) -> Result<(), OrgraftError> {
    let mut writer = BufWriter::new(File::create(path)?);
    writeln!(writer, "new_id\toriginal_id\tlabel\tpartition\tsource")?;
    for row in rows {
        writeln!(
            writer,
            "{}\t{}\t{}\t{}\t{}",
            tsv(&row.new_id),
            tsv(&row.original_id),
            tsv(&row.label),
            tsv(&row.partition),
            tsv(&row.source.display().to_string())
        )?;
    }
    Ok(())
}

#[derive(Debug, Clone)]
struct GzipRuntime {
    program: String,
    threads: Option<usize>,
}

impl GzipRuntime {
    fn resolve(config: &GzipConfig, threads: usize) -> Result<Self, OrgraftError> {
        let program = match config.choice {
            GzipToolChoice::Auto => {
                if command_available("pigz") {
                    "pigz"
                } else {
                    "gzip"
                }
            }
            GzipToolChoice::Pigz => "pigz",
            GzipToolChoice::Gzip => "gzip",
        };

        if !command_available(program) {
            return Err(OrgraftError::InvalidArgument(format!(
                "compression tool `{program}` was not found"
            )));
        }

        Ok(Self {
            program: program.to_string(),
            threads: (program == "pigz").then_some(threads),
        })
    }

    fn decompress_command(&self, path: &Path) -> Command {
        let mut command = Command::new(&self.program);
        if self.program == "pigz" {
            if let Some(threads) = self.threads {
                command.arg("-p").arg(threads.to_string());
            }
        }
        command.arg("-cd").arg(path);
        command
    }

    fn compress_command(&self) -> Command {
        let mut command = Command::new(&self.program);
        if self.program == "pigz" {
            if let Some(threads) = self.threads {
                command.arg("-p").arg(threads.to_string());
            }
        }
        command.arg("-c");
        command
    }
}

fn command_available(program: &str) -> bool {
    Command::new(program)
        .arg("--version")
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .status()
        .is_ok()
}

fn read_text(path: &Path, gzip: &GzipRuntime) -> Result<String, OrgraftError> {
    let mut text = String::new();
    with_text_reader(path, gzip, |reader| {
        reader.read_to_string(&mut text)?;
        Ok(())
    })?;
    Ok(text)
}

fn with_text_reader<T>(
    path: &Path,
    gzip: &GzipRuntime,
    read: impl FnOnce(&mut dyn BufRead) -> Result<T, OrgraftError>,
) -> Result<T, OrgraftError> {
    if is_gzip_path(path) {
        let mut child = gzip
            .decompress_command(path)
            .stdout(Stdio::piped())
            .spawn()
            .map_err(|error| {
                OrgraftError::InvalidArgument(format!(
                    "failed to start {} for {}: {error}",
                    gzip.program,
                    path.display()
                ))
            })?;
        let stdout = child.stdout.take().ok_or_else(|| {
            OrgraftError::InvalidArgument(format!("failed to capture {} stdout", gzip.program))
        })?;
        let mut reader = BufReader::new(stdout);
        let result = read(&mut reader);
        let status = child.wait()?;
        if !status.success() && result.is_ok() {
            return Err(OrgraftError::InvalidArgument(format!(
                "{} failed while reading {}",
                gzip.program,
                path.display()
            )));
        }
        result
    } else {
        let file = File::open(path)?;
        let mut reader = BufReader::new(file);
        read(&mut reader)
    }
}

fn is_gzip_path(path: &Path) -> bool {
    path.file_name()
        .and_then(|value| value.to_str())
        .map(|name| name.to_ascii_lowercase().ends_with(".gz"))
        .unwrap_or(false)
}

fn run_minimap2(
    options: &RecruitOptions,
    bait_fasta: &Path,
    targets: &HashMap<String, TargetInfo>,
    gzip: &GzipRuntime,
) -> Result<AlignmentParseResult, OrgraftError> {
    match options.align_mode {
        AlignMode::Sam => run_minimap2_sam(options, bait_fasta, targets, gzip),
        AlignMode::PafCigar => run_minimap2_paf_cigar(options, bait_fasta, targets, gzip),
    }
}

fn run_minimap2_sam(
    options: &RecruitOptions,
    bait_fasta: &Path,
    targets: &HashMap<String, TargetInfo>,
    gzip: &GzipRuntime,
) -> Result<AlignmentParseResult, OrgraftError> {
    let stderr_path = logs_dir(options).join("minimap2.stderr.log");
    let command_path = logs_dir(options).join("minimap2.command.txt");
    fs::write(
        &command_path,
        format!("{}\n", minimap2_command_preview(options, bait_fasta, gzip)),
    )?;
    let stderr = File::create(&stderr_path)?;

    let use_pigz_pipe = is_gzip_path(&options.reads) && gzip.program == "pigz";
    let mut minimap2 = Command::new(&options.minimap2);
    minimap2
        .arg("-t")
        .arg(options.threads.to_string())
        .arg("-a")
        .arg("-x")
        .arg(&options.preset)
        .arg(bait_fasta);

    let mut decompressor = if use_pigz_pipe {
        let mut decompressor = gzip
            .decompress_command(&options.reads)
            .stdout(Stdio::piped())
            .spawn()
            .map_err(|error| {
                OrgraftError::InvalidArgument(format!(
                    "failed to start pigz for {}: {error}",
                    options.reads.display()
                ))
            })?;
        let decompressed = decompressor.stdout.take().ok_or_else(|| {
            OrgraftError::InvalidArgument("failed to capture pigz stdout".to_string())
        })?;
        minimap2.arg("-");
        minimap2.stdin(Stdio::from(decompressed));
        Some(decompressor)
    } else {
        minimap2.arg(&options.reads);
        None
    };

    let mut child = minimap2
        .stdout(Stdio::piped())
        .stderr(Stdio::from(stderr))
        .spawn()
        .map_err(|error| {
            OrgraftError::InvalidArgument(format!(
                "failed to start minimap2 `{}`: {error}",
                options.minimap2
            ))
        })?;
    let stdout = child.stdout.take().ok_or_else(|| {
        OrgraftError::InvalidArgument("failed to capture minimap2 stdout".to_string())
    })?;
    let mut reader = BufReader::new(stdout);
    let result = parse_sam_reader(&mut reader, targets, options, "minimap2_stdout")?;
    let status = child.wait()?;
    if !status.success() {
        return Err(OrgraftError::InvalidArgument(format!(
            "minimap2 failed; see {}",
            stderr_path.display()
        )));
    }

    if let Some(mut decompressor) = decompressor.take() {
        let decompress_status = decompressor.wait()?;
        if !decompress_status.success() {
            return Err(OrgraftError::InvalidArgument(format!(
                "pigz failed while streaming {} to minimap2",
                options.reads.display()
            )));
        }
    }

    Ok(result)
}

fn run_minimap2_paf_cigar(
    options: &RecruitOptions,
    bait_fasta: &Path,
    targets: &HashMap<String, TargetInfo>,
    gzip: &GzipRuntime,
) -> Result<AlignmentParseResult, OrgraftError> {
    let stderr_path = logs_dir(options).join("minimap2.stderr.log");
    let command_path = logs_dir(options).join("minimap2.command.txt");
    fs::write(
        &command_path,
        format!("{}\n", minimap2_command_preview(options, bait_fasta, gzip)),
    )?;
    let stderr = File::create(&stderr_path)?;

    let use_pigz_pipe = is_gzip_path(&options.reads) && gzip.program == "pigz";
    let mut minimap2 = Command::new(&options.minimap2);
    minimap2
        .arg("-t")
        .arg(options.threads.to_string())
        .arg("-c")
        .arg("-x")
        .arg(&options.preset)
        .arg(bait_fasta);

    let mut decompressor = if use_pigz_pipe {
        let mut decompressor = gzip
            .decompress_command(&options.reads)
            .stdout(Stdio::piped())
            .spawn()
            .map_err(|error| {
                OrgraftError::InvalidArgument(format!(
                    "failed to start pigz for {}: {error}",
                    options.reads.display()
                ))
            })?;
        let decompressed = decompressor.stdout.take().ok_or_else(|| {
            OrgraftError::InvalidArgument("failed to capture pigz stdout".to_string())
        })?;
        minimap2.arg("-");
        minimap2.stdin(Stdio::from(decompressed));
        Some(decompressor)
    } else {
        minimap2.arg(&options.reads);
        None
    };

    let mut child = minimap2
        .stdout(Stdio::piped())
        .stderr(Stdio::from(stderr))
        .spawn()
        .map_err(|error| {
            OrgraftError::InvalidArgument(format!(
                "failed to start minimap2 `{}`: {error}",
                options.minimap2
            ))
        })?;
    let stdout = child.stdout.take().ok_or_else(|| {
        OrgraftError::InvalidArgument("failed to capture minimap2 stdout".to_string())
    })?;
    let mut reader = BufReader::new(stdout);
    let result = parse_paf_cigar_reader(&mut reader, targets, options, "minimap2_stdout")?;
    let status = child.wait()?;
    if !status.success() {
        return Err(OrgraftError::InvalidArgument(format!(
            "minimap2 failed; see {}",
            stderr_path.display()
        )));
    }

    if let Some(mut decompressor) = decompressor.take() {
        let decompress_status = decompressor.wait()?;
        if !decompress_status.success() {
            return Err(OrgraftError::InvalidArgument(format!(
                "pigz failed while streaming {} to minimap2",
                options.reads.display()
            )));
        }
    }

    Ok(result)
}

fn minimap2_command_preview(
    options: &RecruitOptions,
    bait_fasta: &Path,
    gzip: &GzipRuntime,
) -> String {
    let output_flag = match options.align_mode {
        AlignMode::Sam => "-a",
        AlignMode::PafCigar => "-c",
    };
    if is_gzip_path(&options.reads) && gzip.program == "pigz" {
        let mut gzip_parts = vec![gzip.program.clone()];
        if let Some(threads) = gzip.threads {
            gzip_parts.push("-p".to_string());
            gzip_parts.push(threads.to_string());
        }
        gzip_parts.push("-cd".to_string());
        gzip_parts.push(sh_quote(&options.reads.display().to_string()));
        format!(
            "{} | {} -t {} {} -x {} {} -",
            gzip_parts.join(" "),
            sh_quote(&options.minimap2),
            options.threads,
            output_flag,
            sh_quote(&options.preset),
            sh_quote(&bait_fasta.display().to_string())
        )
    } else {
        format!(
            "{} -t {} {} -x {} {} {}",
            sh_quote(&options.minimap2),
            options.threads,
            output_flag,
            sh_quote(&options.preset),
            sh_quote(&bait_fasta.display().to_string()),
            sh_quote(&options.reads.display().to_string())
        )
    }
}

#[derive(Debug)]
struct AlignmentRead {
    labels: BTreeSet<String>,
    partitions: BTreeSet<String>,
    targets: BTreeSet<String>,
    original_targets: BTreeSet<String>,
    best_mapq: u8,
    best_aln_len: u64,
    hit_count: usize,
}

#[derive(Debug)]
struct AlignmentParseResult {
    reads: BTreeMap<String, AlignmentRead>,
    source: String,
    total_alignments: usize,
    passing_alignments: usize,
    skipped_unmapped: usize,
    skipped_secondary: usize,
    skipped_supplementary: usize,
    unknown_targets: usize,
}

fn parse_sam(
    path: &Path,
    targets: &HashMap<String, TargetInfo>,
    options: &RecruitOptions,
    gzip: &GzipRuntime,
) -> Result<AlignmentParseResult, OrgraftError> {
    with_text_reader(path, gzip, |reader| {
        parse_sam_reader(reader, targets, options, &path.display().to_string())
    })
}

fn parse_sam_reader(
    reader: &mut dyn BufRead,
    targets: &HashMap<String, TargetInfo>,
    options: &RecruitOptions,
    source: &str,
) -> Result<AlignmentParseResult, OrgraftError> {
    let mut result = AlignmentParseResult {
        reads: BTreeMap::new(),
        source: source.to_string(),
        total_alignments: 0,
        passing_alignments: 0,
        skipped_unmapped: 0,
        skipped_secondary: 0,
        skipped_supplementary: 0,
        unknown_targets: 0,
    };

    for (index, line) in reader.lines().enumerate() {
        let line = line?;
        if line.trim().is_empty() || line.starts_with('@') {
            continue;
        }
        result.total_alignments += 1;
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 11 {
            return Err(OrgraftError::InvalidArgument(format!(
                "{source}:{} malformed SAM line; expected at least 11 columns",
                index + 1
            )));
        }

        let flag: u16 = fields[1].parse().map_err(|_| {
            OrgraftError::InvalidArgument(format!(
                "{source}:{} invalid SAM flag `{}`",
                index + 1,
                fields[1]
            ))
        })?;
        if flag & 0x4 != 0 {
            result.skipped_unmapped += 1;
            continue;
        }
        if flag & 0x100 != 0 {
            result.skipped_secondary += 1;
            continue;
        }
        if flag & 0x800 != 0 {
            result.skipped_supplementary += 1;
            continue;
        }

        let read_id = fields[0];
        let target_id = fields[2];
        let mapq: u8 = fields[4].parse().map_err(|_| {
            OrgraftError::InvalidArgument(format!(
                "{source}:{} invalid SAM mapq `{}`",
                index + 1,
                fields[4]
            ))
        })?;
        let aln_len = cigar_alignment_len(fields[5])?;
        if mapq < options.min_mapq || aln_len < options.min_aln_len {
            continue;
        }
        let Some(target) = targets.get(target_id) else {
            result.unknown_targets += 1;
            continue;
        };
        result.passing_alignments += 1;
        let read = result
            .reads
            .entry(read_id.to_string())
            .or_insert(AlignmentRead {
                labels: BTreeSet::new(),
                partitions: BTreeSet::new(),
                targets: BTreeSet::new(),
                original_targets: BTreeSet::new(),
                best_mapq: 0,
                best_aln_len: 0,
                hit_count: 0,
            });
        read.labels.insert(target.label.clone());
        read.partitions.insert(target.partition.clone());
        read.targets.insert(target_id.to_string());
        read.original_targets.insert(target.original_id.clone());
        read.best_mapq = read.best_mapq.max(mapq);
        read.best_aln_len = read.best_aln_len.max(aln_len);
        read.hit_count += 1;
    }

    Ok(result)
}

#[derive(Debug, Clone)]
struct BestPafHit {
    read_id: String,
    target_id: String,
    mapq: u8,
    aln_len: u64,
    matches: u64,
}

impl BestPafHit {
    fn is_better_than(&self, other: &Self) -> bool {
        (self.aln_len, self.matches, self.mapq) > (other.aln_len, other.matches, other.mapq)
    }
}

fn parse_paf_cigar_reader(
    reader: &mut dyn BufRead,
    targets: &HashMap<String, TargetInfo>,
    options: &RecruitOptions,
    source: &str,
) -> Result<AlignmentParseResult, OrgraftError> {
    let mut result = AlignmentParseResult {
        reads: BTreeMap::new(),
        source: source.to_string(),
        total_alignments: 0,
        passing_alignments: 0,
        skipped_unmapped: 0,
        skipped_secondary: 0,
        skipped_supplementary: 0,
        unknown_targets: 0,
    };
    let mut best_hits: BTreeMap<String, BestPafHit> = BTreeMap::new();

    for (index, line) in reader.lines().enumerate() {
        let line = line?;
        if line.trim().is_empty() || line.starts_with('#') {
            continue;
        }
        result.total_alignments += 1;
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 12 {
            return Err(OrgraftError::InvalidArgument(format!(
                "{source}:{} malformed PAF line; expected at least 12 columns",
                index + 1
            )));
        }

        let read_id = fields[0];
        let target_id = fields[5];
        let matches: u64 = fields[9].parse().map_err(|_| {
            OrgraftError::InvalidArgument(format!(
                "{source}:{} invalid PAF matching bases `{}`",
                index + 1,
                fields[9]
            ))
        })?;
        let aln_len: u64 = fields[10].parse().map_err(|_| {
            OrgraftError::InvalidArgument(format!(
                "{source}:{} invalid PAF alignment block length `{}`",
                index + 1,
                fields[10]
            ))
        })?;
        let mapq: u8 = fields[11].parse().map_err(|_| {
            OrgraftError::InvalidArgument(format!(
                "{source}:{} invalid PAF mapq `{}`",
                index + 1,
                fields[11]
            ))
        })?;
        if mapq < options.min_mapq || aln_len < options.min_aln_len {
            continue;
        }
        if !targets.contains_key(target_id) {
            result.unknown_targets += 1;
            continue;
        }

        let hit = BestPafHit {
            read_id: read_id.to_string(),
            target_id: target_id.to_string(),
            mapq,
            aln_len,
            matches,
        };
        match best_hits.get(read_id) {
            Some(previous) if !hit.is_better_than(previous) => {}
            _ => {
                best_hits.insert(read_id.to_string(), hit);
            }
        }
    }

    for hit in best_hits.into_values() {
        let Some(target) = targets.get(&hit.target_id) else {
            continue;
        };
        result.passing_alignments += 1;
        let mut read = AlignmentRead {
            labels: BTreeSet::new(),
            partitions: BTreeSet::new(),
            targets: BTreeSet::new(),
            original_targets: BTreeSet::new(),
            best_mapq: hit.mapq,
            best_aln_len: hit.aln_len,
            hit_count: 1,
        };
        read.labels.insert(target.label.clone());
        read.partitions.insert(target.partition.clone());
        read.targets.insert(hit.target_id);
        read.original_targets.insert(target.original_id.clone());
        result.reads.insert(hit.read_id, read);
    }

    Ok(result)
}

fn cigar_alignment_len(cigar: &str) -> Result<u64, OrgraftError> {
    if cigar == "*" {
        return Ok(0);
    }

    let mut length = 0u64;
    let mut number = String::new();
    for character in cigar.chars() {
        if character.is_ascii_digit() {
            number.push(character);
            continue;
        }
        if number.is_empty() {
            return Err(OrgraftError::InvalidArgument(format!(
                "invalid SAM CIGAR `{cigar}`"
            )));
        }
        let value: u64 = number
            .parse()
            .map_err(|_| OrgraftError::InvalidArgument(format!("invalid SAM CIGAR `{cigar}`")))?;
        if matches!(character, 'M' | '=' | 'X' | 'I' | 'D') {
            length += value;
        }
        number.clear();
    }
    if !number.is_empty() {
        return Err(OrgraftError::InvalidArgument(format!(
            "invalid SAM CIGAR `{cigar}`"
        )));
    }
    Ok(length)
}

#[derive(Debug)]
struct Selection {
    selected_all: BTreeSet<String>,
    selected_by_label: BTreeMap<String, BTreeSet<String>>,
    selected_by_partition: BTreeMap<String, BTreeSet<String>>,
    sampled_scopes: BTreeSet<String>,
}

fn select_reads(
    reads: &BTreeMap<String, AlignmentRead>,
    bait: &BaitPreparation,
    options: &RecruitOptions,
) -> Result<Selection, OrgraftError> {
    let mut by_label: BTreeMap<String, BTreeSet<String>> = bait
        .labels
        .iter()
        .map(|label| (label.clone(), BTreeSet::new()))
        .collect();
    let mut by_partition: BTreeMap<String, BTreeSet<String>> = bait
        .partitions
        .keys()
        .map(|partition| (partition.clone(), BTreeSet::new()))
        .collect();

    for (read_id, read) in reads {
        for label in &read.labels {
            by_label
                .entry(label.clone())
                .or_default()
                .insert(read_id.clone());
        }
        for partition in &read.partitions {
            by_partition
                .entry(partition.clone())
                .or_default()
                .insert(read_id.clone());
        }
    }

    let mut sampled_scopes = BTreeSet::new();
    let mut has_partition_limit = false;
    let mut has_label_limit = false;

    for (name, max_reads) in &options.max_reads {
        if let Some(values) = by_partition.get_mut(name) {
            *values = sample_read_ids(values, *max_reads, options.random_seed, name);
            sampled_scopes.insert(name.clone());
            has_partition_limit = true;
        } else if let Some(values) = by_label.get_mut(name) {
            *values = sample_read_ids(values, *max_reads, options.random_seed, name);
            sampled_scopes.insert(name.clone());
            has_label_limit = true;
        } else {
            return Err(OrgraftError::InvalidArgument(format!(
                "--max-reads references unknown label/partition `{name}`"
            )));
        }
    }

    let mut selected_all = BTreeSet::new();
    if has_partition_limit {
        for reads in by_partition.values() {
            selected_all.extend(reads.iter().cloned());
        }
    } else if has_label_limit {
        for reads in by_label.values() {
            selected_all.extend(reads.iter().cloned());
        }
    } else {
        selected_all.extend(reads.keys().cloned());
    }

    for values in by_label.values_mut() {
        *values = values.intersection(&selected_all).cloned().collect();
    }
    for values in by_partition.values_mut() {
        *values = values.intersection(&selected_all).cloned().collect();
    }

    Ok(Selection {
        selected_all,
        selected_by_label: by_label,
        selected_by_partition: by_partition,
        sampled_scopes,
    })
}

fn sample_read_ids(
    ids: &BTreeSet<String>,
    max_reads: usize,
    seed: u64,
    scope: &str,
) -> BTreeSet<String> {
    if ids.len() <= max_reads {
        return ids.clone();
    }

    let mut values: Vec<String> = ids.iter().cloned().collect();
    let mut state = seed ^ stable_hash(scope);
    shuffle(&mut values, &mut state);
    values.truncate(max_reads);
    values.into_iter().collect()
}

fn shuffle<T>(values: &mut [T], state: &mut u64) {
    for index in (1..values.len()).rev() {
        let offset = (next_random(state) as usize) % (index + 1);
        values.swap(index, offset);
    }
}

fn next_random(state: &mut u64) -> u64 {
    if *state == 0 {
        *state = 0x9e37_79b9_7f4a_7c15;
    }
    let mut value = *state;
    value ^= value >> 12;
    value ^= value << 25;
    value ^= value >> 27;
    *state = value;
    value.wrapping_mul(0x2545_f491_4f6c_dd1d)
}

fn stable_hash(value: &str) -> u64 {
    let mut hash = 0xcbf2_9ce4_8422_2325u64;
    for byte in value.as_bytes() {
        hash ^= u64::from(*byte);
        hash = hash.wrapping_mul(0x0000_0100_0000_01b3);
    }
    hash
}

#[derive(Debug)]
struct ReadStatsKey {
    scope: String,
    name: String,
}

impl ReadStatsKey {
    fn new(scope: &str, name: &str) -> Self {
        Self {
            scope: scope.to_string(),
            name: name.to_string(),
        }
    }
}

impl PartialEq for ReadStatsKey {
    fn eq(&self, other: &Self) -> bool {
        self.scope == other.scope && self.name == other.name
    }
}

impl Eq for ReadStatsKey {}

impl PartialOrd for ReadStatsKey {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for ReadStatsKey {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        (&self.scope, &self.name).cmp(&(&other.scope, &other.name))
    }
}

#[derive(Debug, Default)]
struct ReadStatsAccumulator {
    read_count: usize,
    bases: u64,
    min_len: Option<u64>,
    max_len: u64,
    lengths: Vec<u64>,
    mean_qual_sum: f64,
}

impl ReadStatsAccumulator {
    fn observe(&mut self, record: &FastqRecord, mode: &ReadStatsMode) {
        let length = fastq_line_len(&record.seq) as u64;
        self.read_count += 1;
        self.bases += length;
        self.min_len = Some(self.min_len.map_or(length, |current| current.min(length)));
        self.max_len = self.max_len.max(length);
        self.lengths.push(length);
        if *mode == ReadStatsMode::Full {
            self.mean_qual_sum += fastq_mean_quality(&record.qual);
        }
    }

    fn mean_len(&self) -> f64 {
        if self.read_count == 0 {
            0.0
        } else {
            self.bases as f64 / self.read_count as f64
        }
    }

    fn median_len(&self) -> f64 {
        if self.lengths.is_empty() {
            return 0.0;
        }
        let mut lengths = self.lengths.clone();
        lengths.sort_unstable();
        let middle = lengths.len() / 2;
        if lengths.len() % 2 == 0 {
            (lengths[middle - 1] as f64 + lengths[middle] as f64) / 2.0
        } else {
            lengths[middle] as f64
        }
    }

    fn n50(&self) -> u64 {
        if self.lengths.is_empty() {
            return 0;
        }
        let mut lengths = self.lengths.clone();
        lengths.sort_unstable_by(|left, right| right.cmp(left));
        let threshold = self.bases.div_ceil(2);
        let mut cumulative = 0u64;
        for length in lengths {
            cumulative += length;
            if cumulative >= threshold {
                return length;
            }
        }
        0
    }

    fn mean_qual(&self) -> f64 {
        if self.read_count == 0 {
            0.0
        } else {
            self.mean_qual_sum / self.read_count as f64
        }
    }
}

fn fastq_line_len(value: &str) -> usize {
    value.trim_end_matches(['\n', '\r']).len()
}

fn fastq_mean_quality(value: &str) -> f64 {
    let value = value.trim_end_matches(['\n', '\r']);
    if value.is_empty() {
        return 0.0;
    }
    let sum: u64 = value
        .as_bytes()
        .iter()
        .map(|byte| u64::from(byte.saturating_sub(33)))
        .sum();
    sum as f64 / value.len() as f64
}

#[derive(Debug)]
struct FastqCopyStats {
    total_reads: usize,
    selected_reads_seen: usize,
    split_counts: BTreeMap<String, usize>,
    split_outputs: BTreeMap<String, PathBuf>,
    read_stats: BTreeMap<ReadStatsKey, ReadStatsAccumulator>,
}

fn extract_fastq(
    options: &RecruitOptions,
    selection: &Selection,
    gzip: &GzipRuntime,
) -> Result<FastqCopyStats, OrgraftError> {
    let selected_all: HashSet<String> = selection.selected_all.iter().cloned().collect();

    let split_sets = split_sets(options, selection);
    let mut split_sinks = BTreeMap::new();
    let mut split_lookup = BTreeMap::new();
    for (name, ids) in split_sets {
        let path = options.out_dir.join(format!(
            "{}.fastq{}",
            sanitize_name(&name),
            if options.gzip_output { ".gz" } else { "" }
        ));
        split_lookup.insert(name.clone(), ids.into_iter().collect::<HashSet<_>>());
        split_sinks.insert(name, FastqSink::create(&path, options.gzip_output, gzip)?);
    }

    let mut stats = FastqCopyStats {
        total_reads: 0,
        selected_reads_seen: 0,
        split_counts: BTreeMap::new(),
        split_outputs: split_sinks
            .iter()
            .map(|(name, sink)| (name.clone(), sink.path.clone()))
            .collect(),
        read_stats: BTreeMap::new(),
    };
    let split_scope = split_output_scope(options).to_string();

    with_text_reader(&options.reads, gzip, |reader| {
        loop {
            let Some(record) = read_fastq_record(reader)? else {
                break;
            };
            stats.total_reads += 1;
            stats
                .read_stats
                .entry(ReadStatsKey::new("input", "all"))
                .or_default()
                .observe(&record, &options.read_stats_mode);
            if selected_all.contains(&record.id) {
                stats.selected_reads_seen += 1;
            }
            for (name, sink) in split_sinks.iter_mut() {
                if split_lookup
                    .get(name)
                    .map(|ids| ids.contains(&record.id))
                    .unwrap_or(false)
                {
                    sink.write_record(&record)?;
                    *stats.split_counts.entry(name.clone()).or_insert(0) += 1;
                    stats
                        .read_stats
                        .entry(ReadStatsKey::new(&split_scope, name))
                        .or_default()
                        .observe(&record, &options.read_stats_mode);
                }
            }
        }
        Ok(())
    })?;

    for (_, sink) in split_sinks {
        sink.finish()?;
    }

    Ok(stats)
}

fn split_output_scope(options: &RecruitOptions) -> &'static str {
    match options.split_output {
        SplitOutput::None => "output",
        SplitOutput::Label => "label",
        SplitOutput::Partition => "partition",
    }
}

fn split_sets(
    options: &RecruitOptions,
    selection: &Selection,
) -> BTreeMap<String, BTreeSet<String>> {
    match options.split_output {
        SplitOutput::None => BTreeMap::new(),
        SplitOutput::Label => selection.selected_by_label.clone(),
        SplitOutput::Partition => selection.selected_by_partition.clone(),
    }
}

#[derive(Debug)]
struct FastqRecord {
    id: String,
    header: String,
    seq: String,
    plus: String,
    qual: String,
}

fn read_fastq_record(reader: &mut dyn BufRead) -> Result<Option<FastqRecord>, OrgraftError> {
    let mut header = String::new();
    if reader.read_line(&mut header)? == 0 {
        return Ok(None);
    }
    let mut seq = String::new();
    let mut plus = String::new();
    let mut qual = String::new();
    if reader.read_line(&mut seq)? == 0
        || reader.read_line(&mut plus)? == 0
        || reader.read_line(&mut qual)? == 0
    {
        return Err(OrgraftError::InvalidArgument(
            "truncated FASTQ record".to_string(),
        ));
    }
    if !header.starts_with('@') {
        return Err(OrgraftError::InvalidArgument(format!(
            "FASTQ record header does not start with @: {}",
            header.trim_end()
        )));
    }
    if !plus.starts_with('+') {
        return Err(OrgraftError::InvalidArgument(format!(
            "FASTQ record plus line does not start with + for {}",
            header.trim_end()
        )));
    }
    let id = header[1..]
        .split_whitespace()
        .next()
        .unwrap_or("")
        .to_string();
    if id.is_empty() {
        return Err(OrgraftError::InvalidArgument(
            "FASTQ record has empty read ID".to_string(),
        ));
    }
    Ok(Some(FastqRecord {
        id,
        header,
        seq,
        plus,
        qual,
    }))
}

struct FastqSink {
    path: PathBuf,
    writer: Option<Box<dyn Write>>,
    child: Option<Child>,
}

impl FastqSink {
    fn create(path: &Path, gzip_output: bool, gzip: &GzipRuntime) -> Result<Self, OrgraftError> {
        if gzip_output {
            let output = File::create(path)?;
            let mut child = gzip
                .compress_command()
                .stdin(Stdio::piped())
                .stdout(Stdio::from(output))
                .spawn()
                .map_err(|error| {
                    OrgraftError::InvalidArgument(format!(
                        "failed to start {} for {}: {error}",
                        gzip.program,
                        path.display()
                    ))
                })?;
            let stdin = child.stdin.take().ok_or_else(|| {
                OrgraftError::InvalidArgument(format!("failed to open {} stdin", gzip.program))
            })?;
            Ok(Self {
                path: path.to_path_buf(),
                writer: Some(Box::new(BufWriter::new(stdin))),
                child: Some(child),
            })
        } else {
            Ok(Self {
                path: path.to_path_buf(),
                writer: Some(Box::new(BufWriter::new(File::create(path)?))),
                child: None,
            })
        }
    }

    fn write_record(&mut self, record: &FastqRecord) -> Result<(), OrgraftError> {
        let writer = self.writer.as_mut().ok_or_else(|| {
            OrgraftError::InvalidArgument(format!(
                "output {} is already closed",
                self.path.display()
            ))
        })?;
        writer.write_all(record.header.as_bytes())?;
        writer.write_all(record.seq.as_bytes())?;
        writer.write_all(record.plus.as_bytes())?;
        writer.write_all(record.qual.as_bytes())?;
        Ok(())
    }

    fn finish(mut self) -> Result<(), OrgraftError> {
        if let Some(mut writer) = self.writer.take() {
            writer.flush()?;
        }
        if let Some(mut child) = self.child.take() {
            let status = child.wait()?;
            if !status.success() {
                return Err(OrgraftError::InvalidArgument(format!(
                    "compression failed while writing {}",
                    self.path.display()
                )));
            }
        }
        Ok(())
    }
}

fn write_read_classification(
    path: &Path,
    reads: &BTreeMap<String, AlignmentRead>,
    selected: &BTreeSet<String>,
) -> Result<(), OrgraftError> {
    let mut writer = BufWriter::new(File::create(path)?);
    writeln!(
        writer,
        "read_id\tclassification\tlabels\tpartitions\ttarget_ids\toriginal_target_ids\tbest_mapq\tbest_aln_len\thit_count\tselected"
    )?;
    for (read_id, read) in reads {
        let classification = if read.labels.len() <= 1 {
            "organelle"
        } else {
            "ambiguous"
        };
        writeln!(
            writer,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            tsv(read_id),
            classification,
            join_set(&read.labels),
            join_set(&read.partitions),
            join_set(&read.targets),
            join_set(&read.original_targets),
            read.best_mapq,
            read.best_aln_len,
            read.hit_count,
            selected.contains(read_id)
        )?;
    }
    Ok(())
}

fn write_sampled_ids(path: &Path, selection: &Selection) -> Result<(), OrgraftError> {
    let mut writer = BufWriter::new(File::create(path)?);
    writeln!(writer, "scope\tname\tread_id")?;
    for read_id in &selection.selected_all {
        writeln!(writer, "all\tall\t{}", tsv(read_id))?;
    }
    for (label, reads) in &selection.selected_by_label {
        for read_id in reads {
            writeln!(writer, "label\t{}\t{}", tsv(label), tsv(read_id))?;
        }
    }
    for (partition, reads) in &selection.selected_by_partition {
        for read_id in reads {
            writeln!(writer, "partition\t{}\t{}", tsv(partition), tsv(read_id))?;
        }
    }
    Ok(())
}

fn write_read_stats(
    path: &Path,
    stats: &BTreeMap<ReadStatsKey, ReadStatsAccumulator>,
    mode: &ReadStatsMode,
) -> Result<(), OrgraftError> {
    let mut writer = BufWriter::new(File::create(path)?);
    match mode {
        ReadStatsMode::Basic => {
            writeln!(
                writer,
                "scope\tname\tread_count\tbases\tmin_len\tmean_len\tn50\tmax_len"
            )?;
        }
        ReadStatsMode::Full => {
            writeln!(
                writer,
                "scope\tname\tread_count\tbases\tmin_len\tmean_len\tmedian_len\tn50\tmax_len\tmean_qual"
            )?;
        }
    }
    for (key, value) in stats {
        match mode {
            ReadStatsMode::Basic => {
                writeln!(
                    writer,
                    "{}\t{}\t{}\t{}\t{}\t{:.2}\t{}\t{}",
                    tsv(&key.scope),
                    tsv(&key.name),
                    value.read_count,
                    value.bases,
                    value.min_len.unwrap_or(0),
                    value.mean_len(),
                    value.n50(),
                    value.max_len,
                )?;
            }
            ReadStatsMode::Full => {
                writeln!(
                    writer,
                    "{}\t{}\t{}\t{}\t{}\t{:.2}\t{:.2}\t{}\t{}\t{:.2}",
                    tsv(&key.scope),
                    tsv(&key.name),
                    value.read_count,
                    value.bases,
                    value.min_len.unwrap_or(0),
                    value.mean_len(),
                    value.median_len(),
                    value.n50(),
                    value.max_len,
                    value.mean_qual(),
                )?;
            }
        }
    }
    Ok(())
}

struct SummaryContext<'a> {
    options: &'a RecruitOptions,
    gzip: &'a GzipRuntime,
    bait: &'a BaitPreparation,
    alignment: &'a AlignmentParseResult,
    selection: &'a Selection,
    copy_stats: &'a FastqCopyStats,
}

fn write_summary(path: &Path, context: &SummaryContext<'_>) -> Result<(), OrgraftError> {
    let mut writer = BufWriter::new(File::create(path)?);
    writeln!(writer, "scope\tname\tmetric\tvalue")?;
    summary_row(&mut writer, "global", "run", "mode", "reference")?;
    summary_row(
        &mut writer,
        "global",
        "run",
        "reads_input",
        &context.options.reads.display().to_string(),
    )?;
    summary_row(
        &mut writer,
        "global",
        "run",
        "alignment_format",
        alignment_format_label(context.options),
    )?;
    summary_row(
        &mut writer,
        "global",
        "run",
        "alignment_source",
        &context.alignment.source,
    )?;
    summary_row(
        &mut writer,
        "global",
        "run",
        "minimap2_preset",
        &context.options.preset,
    )?;
    summary_row(
        &mut writer,
        "global",
        "run",
        "minimap2_threads",
        &context.options.threads.to_string(),
    )?;
    summary_row(
        &mut writer,
        "global",
        "run",
        "minimap2_command",
        &context
            .options
            .sam
            .as_ref()
            .map(|_| "not_run_precomputed_sam".to_string())
            .unwrap_or_else(|| {
                minimap2_command_preview(context.options, &context.bait.fasta_path, context.gzip)
            }),
    )?;
    summary_row(
        &mut writer,
        "global",
        "run",
        "bait_fasta",
        &context.bait.fasta_path.display().to_string(),
    )?;
    summary_row(
        &mut writer,
        "global",
        "run",
        "bait_id_map",
        context
            .bait
            .id_map_path
            .as_ref()
            .map(|path| path.display().to_string())
            .as_deref()
            .unwrap_or("not_requested"),
    )?;
    summary_row(
        &mut writer,
        "global",
        "run",
        "read_stats",
        &read_stats_path(context.options).display().to_string(),
    )?;
    summary_row(
        &mut writer,
        "global",
        "run",
        "read_stats_mode",
        read_stats_mode_label(&context.options.read_stats_mode),
    )?;
    summary_row(
        &mut writer,
        "global",
        "run",
        "gzip_tool",
        &context.gzip.program,
    )?;
    summary_row(
        &mut writer,
        "global",
        "run",
        "gzip_threads",
        &context
            .gzip
            .threads
            .map(|threads| threads.to_string())
            .unwrap_or_else(|| "not_applicable".to_string()),
    )?;
    summary_row(
        &mut writer,
        "global",
        "counts",
        "bait_records",
        &context.bait.record_count.to_string(),
    )?;
    summary_row(
        &mut writer,
        "global",
        "counts",
        "fastq_reads_scanned",
        &context.copy_stats.total_reads.to_string(),
    )?;
    summary_row(
        &mut writer,
        "global",
        "counts",
        "alignments_total",
        &context.alignment.total_alignments.to_string(),
    )?;
    summary_row(
        &mut writer,
        "global",
        "counts",
        "alignments_passing",
        &context.alignment.passing_alignments.to_string(),
    )?;
    summary_row(
        &mut writer,
        "global",
        "counts",
        "alignments_skipped_unmapped",
        &context.alignment.skipped_unmapped.to_string(),
    )?;
    summary_row(
        &mut writer,
        "global",
        "counts",
        "alignments_skipped_secondary",
        &context.alignment.skipped_secondary.to_string(),
    )?;
    summary_row(
        &mut writer,
        "global",
        "counts",
        "alignments_skipped_supplementary",
        &context.alignment.skipped_supplementary.to_string(),
    )?;
    summary_row(
        &mut writer,
        "global",
        "counts",
        "alignments_unknown_targets",
        &context.alignment.unknown_targets.to_string(),
    )?;
    summary_row(
        &mut writer,
        "global",
        "counts",
        "recruited_reads",
        &context.alignment.reads.len().to_string(),
    )?;
    summary_row(
        &mut writer,
        "global",
        "counts",
        "selected_reads",
        &context.selection.selected_all.len().to_string(),
    )?;
    summary_row(
        &mut writer,
        "global",
        "counts",
        "written_reads",
        &context.copy_stats.selected_reads_seen.to_string(),
    )?;
    summary_row(
        &mut writer,
        "global",
        "counts",
        "selected_reads_missing_from_fastq",
        &context
            .selection
            .selected_all
            .len()
            .saturating_sub(context.copy_stats.selected_reads_seen)
            .to_string(),
    )?;
    summary_row(
        &mut writer,
        "global",
        "sampling",
        "sampled_scopes",
        &join_set(&context.selection.sampled_scopes),
    )?;

    for (label, reads) in &context.selection.selected_by_label {
        let recruited = context
            .alignment
            .reads
            .values()
            .filter(|read| read.labels.contains(label))
            .count();
        summary_row(
            &mut writer,
            "label",
            label,
            "recruited_reads",
            &recruited.to_string(),
        )?;
        summary_row(
            &mut writer,
            "label",
            label,
            "selected_reads",
            &reads.len().to_string(),
        )?;
        if let Some(count) = context.copy_stats.split_counts.get(label) {
            summary_row(
                &mut writer,
                "label",
                label,
                "written_reads",
                &count.to_string(),
            )?;
        }
        if let Some(path) = context.copy_stats.split_outputs.get(label) {
            summary_row(
                &mut writer,
                "label",
                label,
                "output_fastq",
                &path.display().to_string(),
            )?;
        }
    }

    for (partition, reads) in &context.selection.selected_by_partition {
        let recruited = context
            .alignment
            .reads
            .values()
            .filter(|read| read.partitions.contains(partition))
            .count();
        summary_row(
            &mut writer,
            "partition",
            partition,
            "label",
            context
                .bait
                .partitions
                .get(partition)
                .map(String::as_str)
                .unwrap_or("unknown"),
        )?;
        summary_row(
            &mut writer,
            "partition",
            partition,
            "recruited_reads",
            &recruited.to_string(),
        )?;
        summary_row(
            &mut writer,
            "partition",
            partition,
            "selected_reads",
            &reads.len().to_string(),
        )?;
        if let Some(count) = context.copy_stats.split_counts.get(partition) {
            summary_row(
                &mut writer,
                "partition",
                partition,
                "written_reads",
                &count.to_string(),
            )?;
        }
        if let Some(path) = context.copy_stats.split_outputs.get(partition) {
            summary_row(
                &mut writer,
                "partition",
                partition,
                "output_fastq",
                &path.display().to_string(),
            )?;
        }
    }

    Ok(())
}

fn summary_row(
    writer: &mut dyn Write,
    scope: &str,
    name: &str,
    metric: &str,
    value: &str,
) -> Result<(), OrgraftError> {
    writeln!(
        writer,
        "{}\t{}\t{}\t{}",
        tsv(scope),
        tsv(name),
        tsv(metric),
        tsv(value)
    )?;
    Ok(())
}

fn tsv(value: &str) -> String {
    value.replace(['\t', '\n', '\r'], " ")
}

fn alignment_format_label(options: &RecruitOptions) -> &'static str {
    if options.sam.is_some() {
        "SAM"
    } else {
        match options.align_mode {
            AlignMode::Sam => "SAM",
            AlignMode::PafCigar => "PAF-CG",
        }
    }
}

fn read_stats_mode_label(mode: &ReadStatsMode) -> &'static str {
    match mode {
        ReadStatsMode::Basic => "basic",
        ReadStatsMode::Full => "full",
    }
}

fn sh_quote(value: &str) -> String {
    if value.chars().all(|character| {
        character.is_ascii_alphanumeric() || matches!(character, '/' | '.' | '_' | '-' | ':')
    }) {
        value.to_string()
    } else {
        format!("'{}'", value.replace('\'', "'\\''"))
    }
}

fn join_set(values: &BTreeSet<String>) -> String {
    if values.is_empty() {
        ".".to_string()
    } else {
        values.iter().cloned().collect::<Vec<_>>().join(",")
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_labeled_bait_argument() {
        let bait = parse_bait_input("mito=seed.fa");
        assert_eq!(bait.label.as_deref(), Some("mito"));
        assert_eq!(bait.path, PathBuf::from("seed.fa"));
    }

    #[test]
    fn parses_fasta_records() {
        let records = parse_fasta(">a desc\nAC\nGT\n>b\nTT\n", Path::new("bait.fa")).unwrap();
        assert_eq!(records.len(), 2);
        assert_eq!(records[0].id, "a");
        assert_eq!(records[0].seq, "ACGT");
        assert_eq!(records[1].id, "b");
    }

    #[test]
    fn splits_gfa_by_connected_components() {
        let graph = parse_gfa(
            "S\tA\tAC\nS\tB\tGT\nS\tC\tTA\nL\tA\t+\tB\t+\t0M\n",
            Path::new("toy.gfa"),
        )
        .unwrap();
        let components = graph.connected_components();
        assert_eq!(
            components,
            vec![
                vec!["A".to_string(), "B".to_string()],
                vec!["C".to_string()]
            ]
        );
    }

    #[test]
    fn samples_are_deterministic_and_capped() {
        let ids = ["r1", "r2", "r3", "r4"]
            .into_iter()
            .map(String::from)
            .collect::<BTreeSet<_>>();
        let first = sample_read_ids(&ids, 2, 42, "mito");
        let second = sample_read_ids(&ids, 2, 42, "mito");
        assert_eq!(first, second);
        assert_eq!(first.len(), 2);
    }

    #[test]
    fn command_options_keep_iterative_interface() {
        let args = vec![
            "--reads".to_string(),
            "reads.fastq".to_string(),
            "--bait".to_string(),
            "mito=seed.fa".to_string(),
            "--mode".to_string(),
            "iterative".to_string(),
            "--iterations".to_string(),
            "3".to_string(),
        ];
        let options = RecruitOptions::from_args(&args).unwrap();
        assert_eq!(options.mode, RecruitmentMode::Iterative);
        assert_eq!(options.iterations, 3);
    }

    #[test]
    fn default_align_mode_is_paf_cigar() {
        let args = vec![
            "--reads".to_string(),
            "reads.fastq".to_string(),
            "--bait".to_string(),
            "mito=seed.fa".to_string(),
        ];
        let options = RecruitOptions::from_args(&args).unwrap();
        assert_eq!(options.align_mode, AlignMode::PafCigar);
    }

    #[test]
    fn default_random_seed_is_42() {
        let options = RecruitOptions::from_args(&[
            "--reads".to_string(),
            "reads.fastq".to_string(),
            "--bait".to_string(),
            "mito=seed.fa".to_string(),
        ])
        .unwrap();
        assert_eq!(options.random_seed, 42);
    }

    #[test]
    fn default_read_stats_mode_is_basic() {
        let options = RecruitOptions::from_args(&[
            "--reads".to_string(),
            "reads.fastq".to_string(),
            "--bait".to_string(),
            "mito=seed.fa".to_string(),
        ])
        .unwrap();
        assert_eq!(options.read_stats_mode, ReadStatsMode::Basic);
        assert_eq!(parse_read_stats_mode("full").unwrap(), ReadStatsMode::Full);
    }

    #[test]
    fn gzip_output_defaults_on_and_accepts_on_off() {
        let default_options = RecruitOptions::from_args(&[
            "--reads".to_string(),
            "reads.fastq".to_string(),
            "--bait".to_string(),
            "mito=seed.fa".to_string(),
        ])
        .unwrap();
        assert!(default_options.gzip_output);

        let off_options = RecruitOptions::from_args(&[
            "--reads".to_string(),
            "reads.fastq".to_string(),
            "--bait".to_string(),
            "mito=seed.fa".to_string(),
            "--gzip-output".to_string(),
            "off".to_string(),
        ])
        .unwrap();
        assert!(!off_options.gzip_output);

        let on_options = RecruitOptions::from_args(&[
            "--reads".to_string(),
            "reads.fastq".to_string(),
            "--bait".to_string(),
            "mito=seed.fa".to_string(),
            "--gzip-output".to_string(),
            "on".to_string(),
        ])
        .unwrap();
        assert!(on_options.gzip_output);

        let error = RecruitOptions::from_args(&[
            "--reads".to_string(),
            "reads.fastq".to_string(),
            "--bait".to_string(),
            "mito=seed.fa".to_string(),
            "--gzip-output".to_string(),
            "maybe".to_string(),
        ])
        .unwrap_err();
        assert_eq!(
            error.to_string(),
            "unknown on/off value `maybe`; expected on or off"
        );
    }

    #[test]
    fn max_reads_all_is_not_supported() {
        let error = parse_max_reads("all,20000").unwrap_err();
        assert!(error
            .to_string()
            .contains("--max-reads all,N is not supported"));
    }

    #[test]
    fn gzip_threads_is_unified_into_threads() {
        assert!(HELP.contains("--threads N               minimap2 and pigz threads [1]"));
        assert!(HELP.contains("--advanced-help"));
        assert!(!HELP.contains("--gzip-threads"));

        let error = RecruitOptions::from_args(&[
            "--reads".to_string(),
            "reads.fastq".to_string(),
            "--bait".to_string(),
            "mito=seed.fa".to_string(),
            "--gzip-threads".to_string(),
            "8".to_string(),
        ])
        .unwrap_err();
        assert!(error
            .to_string()
            .contains("--gzip-threads has been removed; use --threads"));
    }

    #[test]
    fn advanced_help_succeeds_without_required_inputs() {
        assert!(run(&["--advanced-help".to_string()]).is_ok());
    }

    #[test]
    fn main_help_is_compact() {
        assert!(HELP.contains("Inputs:"));
        assert!(HELP.contains("Outputs:"));
        assert!(HELP.contains("Additional Parameters:"));
        assert!(HELP.contains("Layout: OUT/{*.fastq[.gz],logs}"));
        assert!(HELP.contains(
            "orgraft recruit --reads reads.fastq[.gz] --mito mito.fa --plastid plastid.fa [options]"
        ));
        assert!(HELP.contains("--mito FILE"));
        assert!(HELP.contains("--plastid FILE"));
        assert!(HELP.contains(
            "--gzip-output MODE        on|off; write recruited FASTQ files as .fastq.gz [on]"
        ));
        assert!(!HELP.contains("--bait LABEL=FILE"));
        assert!(HELP.contains("--platform NAME"));
        assert!(HELP.contains("--bait-format FORMAT"));
        assert!(!HELP.contains("--split-output MODE"));
        assert!(!HELP.contains("Additional bait inputs:"));
        assert!(!HELP.contains("Debug outputs:"));
        assert!(!HELP.contains("logs/minimap2.command.txt"));
        assert!(!HELP.contains("--random-seed N"));
        assert!(!HELP.contains("--read-stats MODE"));
        assert!(HELP.contains(
            "--max-reads NAME,N        cap selected reads for a bait label or partition"
        ));
        assert!(!HELP.contains("--preset PRESET"));
        assert!(!HELP.contains("--min-mapq N"));
        assert!(!HELP.contains("--min-aln-len N"));
    }

    #[test]
    fn advanced_help_exposes_sampling_debug_and_logs() {
        assert!(ADVANCED_HELP
            .contains("--random-seed N               deterministic sampling seed [42]"));
        assert!(ADVANCED_HELP.contains("--read-stats MODE             basic, full [basic]"));
        assert!(!ADVANCED_HELP.contains("--max-reads NAME,N"));
        assert!(ADVANCED_HELP.contains("Name/ID options:"));
        assert!(ADVANCED_HELP.contains("Split options:"));
        assert!(ADVANCED_HELP.contains("Write options:"));
        assert!(ADVANCED_HELP.contains("Other options:"));
        assert!(
            ADVANCED_HELP.contains("--split-output MODE           label, partition, none [label]")
        );
        assert!(!ADVANCED_HELP.contains("Bait inputs:"));
        assert!(!ADVANCED_HELP.contains("Sampling:"));
        assert!(!ADVANCED_HELP.contains("Debug outputs:"));
        assert!(ADVANCED_HELP
            .contains("--bait LABEL=FILE             generic FASTA/GFA bait; repeatable"));
        assert!(!ADVANCED_HELP.contains("--bait-format FORMAT"));
        assert!(!ADVANCED_HELP.contains("Recruitment:"));
        assert!(!ADVANCED_HELP.contains("--minimap2 PATH"));
        assert!(!ADVANCED_HELP.contains("--align-mode MODE"));
        assert!(!ADVANCED_HELP.contains("--platform NAME"));
        assert!(!ADVANCED_HELP.contains("--sam FILE"));
        assert!(
            ADVANCED_HELP.contains("--write-read-classification   write read_classification.tsv")
        );
        assert!(ADVANCED_HELP
            .contains("--write-bait-partitions       write bait partition FASTA files"));
        assert!(!ADVANCED_HELP.contains("Compression:"));
        assert!(!ADVANCED_HELP.contains("--gzip-tool"));
        assert!(!ADVANCED_HELP.contains("Reserved interface:"));
        assert!(!ADVANCED_HELP.contains("--mode reference|iterative"));
        assert!(!ADVANCED_HELP.contains("--iterations N"));
        assert!(!ADVANCED_HELP.contains("Outputs:"));
        assert!(!ADVANCED_HELP.contains("<label>.fastq"));
        assert!(!ADVANCED_HELP.contains("Logs:"));
        assert!(!ADVANCED_HELP.contains("logs/minimap2.command.txt"));
        assert!(!ADVANCED_HELP.contains("logs/recruitment_summary.tsv"));
        assert!(!ADVANCED_HELP.contains("cap selected reads for all"));
    }

    #[test]
    fn help_uses_logs_for_intermediates() {
        assert!(!ADVANCED_HELP.contains("logs/bait.fasta"));
        assert!(!ADVANCED_HELP.contains("logs/minimap2.command.txt"));
        assert!(!ADVANCED_HELP.contains("logs/minimap2.stderr.log"));
        assert!(!ADVANCED_HELP.contains("logs/recruitment_summary.tsv"));
        assert!(!ADVANCED_HELP.contains("logs/read_stats.tsv"));
        assert!(!ADVANCED_HELP.contains("work/bait.fasta"));
        assert!(!ADVANCED_HELP.contains("--report"));
        assert!(!ADVANCED_HELP.contains("recruitment_report.md"));
    }

    #[test]
    fn read_fastq_record_uses_first_header_token() {
        let text = b"@read1 comment\nAC\n+\n!!\n";
        let mut reader = BufReader::new(&text[..]);
        let record = read_fastq_record(&mut reader).unwrap().unwrap();
        assert_eq!(record.id, "read1");
    }

    #[test]
    fn run_with_existing_sam_writes_split_fastq_and_summary() {
        let root = unique_test_dir("orgraft_recruit_e2e");
        fs::create_dir_all(&root).unwrap();
        let reads = root.join("reads.fastq");
        let mito = root.join("mito.fa");
        let plastid = root.join("plastid.fa");
        let sam = root.join("hits.sam");
        let out_dir = root.join("out");

        fs::write(
            &reads,
            "@r1 comment\nACGTACGT\n+\n!!!!!!!!\n@r2\nTTTT\n+\n!!!!\n@r3\nGGGG\n+\n!!!!\n@r4\nCCCC\n+\n!!!!\n",
        )
        .unwrap();
        fs::write(&mito, ">old_mito\nACGTACGT\n").unwrap();
        fs::write(&plastid, ">old_plastid\nTTTT\n").unwrap();
        fs::write(
            &sam,
            "@SQ\tSN:mito_1\tLN:8\n@SQ\tSN:plastid_1\tLN:4\n\
r1\t0\tmito_1\t1\t60\t8M\t*\t0\t0\tACGTACGT\t!!!!!!!!\n\
r2\t0\tplastid_1\t1\t60\t4M\t*\t0\t0\tTTTT\t!!!!\n\
r3\t4\t*\t0\t0\t*\t*\t0\t0\tGGGG\t!!!!\n\
r4\t2048\tmito_1\t1\t60\t4M\t*\t0\t0\tCCCC\t!!!!\n",
        )
        .unwrap();

        let args = vec![
            "--reads".to_string(),
            reads.display().to_string(),
            "--bait".to_string(),
            format!("mito={}", mito.display()),
            "--bait".to_string(),
            format!("plastid={}", plastid.display()),
            "--rename-bait".to_string(),
            "--write-id-map".to_string(),
            "--split-output".to_string(),
            "label".to_string(),
            "--sam".to_string(),
            sam.display().to_string(),
            "--out-dir".to_string(),
            out_dir.display().to_string(),
            "--gzip-output".to_string(),
            "off".to_string(),
        ];

        run(&args).unwrap();

        assert!(!out_dir.join("organelle_reads.fastq").exists());
        assert!(fs::read_to_string(out_dir.join("mito.fastq"))
            .unwrap()
            .contains("@r1 comment\n"));
        assert!(fs::read_to_string(out_dir.join("plastid.fastq"))
            .unwrap()
            .contains("@r2\n"));
        assert!(!out_dir.join("recruitment_summary.tsv").exists());
        assert!(!out_dir.join("recruitment_report.md").exists());
        let summary =
            fs::read_to_string(out_dir.join("logs").join("recruitment_summary.tsv")).unwrap();
        assert!(summary.contains("global\tcounts\tselected_reads\t2\n"));
        assert!(summary.contains("global\trun\tminimap2_command\tnot_run_precomputed_sam\n"));
        assert!(summary.contains("global\trun\tread_stats\t"));
        assert!(summary.contains("global\trun\tread_stats_mode\tbasic\n"));
        assert!(summary.contains("global\tcounts\talignments_skipped_unmapped\t1\n"));
        assert!(summary.contains("global\tcounts\talignments_skipped_supplementary\t1\n"));
        let read_stats = fs::read_to_string(out_dir.join("logs").join("read_stats.tsv")).unwrap();
        assert!(read_stats
            .contains("scope\tname\tread_count\tbases\tmin_len\tmean_len\tn50\tmax_len\n"));
        assert!(!read_stats.contains("mean_qual"));
        assert!(!read_stats.contains("median_len"));
        assert!(read_stats.contains("input\tall\t4\t20\t4\t5.00\t4\t8\n"));
        assert!(read_stats.contains("label\tmito\t1\t8\t8\t8.00\t8\t8\n"));
        assert!(read_stats.contains("label\tplastid\t1\t4\t4\t4.00\t4\t4\n"));

        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn cigar_alignment_len_counts_alignment_ops() {
        assert_eq!(cigar_alignment_len("5S10M1I2D3=4X5H").unwrap(), 20);
        assert_eq!(cigar_alignment_len("*").unwrap(), 0);
    }

    #[test]
    fn fastq_mean_quality_uses_phred33() {
        assert_eq!(fastq_mean_quality("!!!!\n"), 0.0);
        assert_eq!(fastq_mean_quality("IIII\n"), 40.0);
    }

    #[test]
    fn full_read_stats_include_median_and_quality() {
        let mut stats = BTreeMap::new();
        let mut accumulator = ReadStatsAccumulator::default();
        accumulator.observe(
            &FastqRecord {
                id: "r1".to_string(),
                header: "@r1\n".to_string(),
                seq: "ACGT\n".to_string(),
                plus: "+\n".to_string(),
                qual: "IIII\n".to_string(),
            },
            &ReadStatsMode::Full,
        );
        stats.insert(ReadStatsKey::new("input", "all"), accumulator);

        let root = unique_test_dir("orgraft_read_stats_full");
        fs::create_dir_all(&root).unwrap();
        let path = root.join("read_stats.tsv");
        write_read_stats(&path, &stats, &ReadStatsMode::Full).unwrap();
        let text = fs::read_to_string(&path).unwrap();
        assert!(text.contains("median_len"));
        assert!(text.contains("mean_qual"));
        assert!(text.contains("input\tall\t1\t4\t4\t4.00\t4.00\t4\t4\t40.00\n"));
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn paf_cigar_reader_assigns_best_target_per_read() {
        let options = RecruitOptions::from_args(&[
            "--reads".to_string(),
            "reads.fastq".to_string(),
            "--bait".to_string(),
            "mito=mito.fa".to_string(),
        ])
        .unwrap();
        let mut targets = HashMap::new();
        targets.insert(
            "mito_1".to_string(),
            TargetInfo {
                label: "mito".to_string(),
                partition: "mito".to_string(),
                original_id: "old_mito".to_string(),
            },
        );
        targets.insert(
            "plastid_1".to_string(),
            TargetInfo {
                label: "plastid".to_string(),
                partition: "plastid".to_string(),
                original_id: "old_plastid".to_string(),
            },
        );

        let paf = b"r1\t100\t0\t90\t+\tmito_1\t1000\t0\t90\t80\t90\t60\tcg:Z:90M\n\
r1\t100\t0\t40\t+\tplastid_1\t1000\t0\t40\t40\t40\t60\tcg:Z:40M\n\
r2\t100\t0\t20\t+\tmito_1\t1000\t0\t20\t20\t20\t60\tcg:Z:20M\n\
r2\t100\t0\t80\t+\tplastid_1\t1000\t0\t80\t75\t80\t60\tcg:Z:80M\n";
        let mut reader = BufReader::new(&paf[..]);
        let result = parse_paf_cigar_reader(&mut reader, &targets, &options, "test.paf").unwrap();

        assert!(result.reads["r1"].labels.contains("mito"));
        assert!(!result.reads["r1"].labels.contains("plastid"));
        assert!(result.reads["r2"].labels.contains("plastid"));
        assert!(!result.reads["r2"].labels.contains("mito"));
    }

    fn unique_test_dir(prefix: &str) -> PathBuf {
        let nanos = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        std::env::temp_dir().join(format!("{prefix}_{}_{}", std::process::id(), nanos))
    }
}
