use std::cmp::Ordering;
use std::collections::{HashMap, HashSet, VecDeque};
use std::fs::{self, File};
use std::io::{self, BufRead, BufReader, BufWriter, Read, Write};
use std::path::{Path, PathBuf};
use std::process::{Child, ChildStdout, Command, Stdio};
use std::time::Instant;

// Assembly engine for `orgraft asm`: resolve profiles into Config, run the
// numbered graph-construction steps, and write algorithm-level intermediates.
#[derive(Debug, Clone, Copy, Eq, PartialEq)]
enum OrganelleProfile {
    Plastid,
    Mito,
}

impl OrganelleProfile {
    fn as_str(self) -> &'static str {
        match self {
            OrganelleProfile::Plastid => "plastid",
            OrganelleProfile::Mito => "mito",
        }
    }
}

#[derive(Debug, Clone, Copy, Eq, PartialEq)]
enum NumtInterference {
    Low,
    High,
}

impl NumtInterference {
    fn as_str(self) -> &'static str {
        match self {
            NumtInterference::Low => "low",
            NumtInterference::High => "high",
        }
    }
}

#[derive(Debug, Clone, Copy, Eq, PartialEq)]
enum DataMode {
    Standard,
    Low,
}

impl DataMode {
    fn as_str(self) -> &'static str {
        match self {
            DataMode::Standard => "standard",
            DataMode::Low => "low",
        }
    }
}

#[derive(Debug, Clone, Copy, Eq, PartialEq)]
enum MitoStableSelectionMode {
    Auto,
    #[cfg(test)]
    ForbidSelected,
    AllowSelected,
}

impl MitoStableSelectionMode {
    #[cfg(test)]
    fn from_str(value: &str) -> Option<Self> {
        match value {
            "auto" => Some(MitoStableSelectionMode::Auto),
            "forbid-selected" | "forbid_selected" | "forbid" => {
                Some(MitoStableSelectionMode::ForbidSelected)
            }
            "allow-selected" | "allow_selected" | "allow" => {
                Some(MitoStableSelectionMode::AllowSelected)
            }
            _ => None,
        }
    }

    fn as_str(self) -> &'static str {
        match self {
            MitoStableSelectionMode::Auto => "auto",
            #[cfg(test)]
            MitoStableSelectionMode::ForbidSelected => "forbid-selected",
            MitoStableSelectionMode::AllowSelected => "allow-selected",
        }
    }
}

#[derive(Debug, Default)]
struct OverrideFlags {
    k: bool,
    s: bool,
    numt_interference: bool,
    min_anchor_coverage: bool,
    min_edge_coverage: bool,
    min_branch_ratio: bool,
    max_edges_per_state: bool,
    dedup_kmer: bool,
    containment_ratio: bool,
    min_tip_len: bool,
    min_link_support: bool,
    min_link_ratio: bool,
    read_junction_links: bool,
    bidirectional_links: bool,
    read_subsets: bool,
    rounds: bool,
    skeleton_min_link_support: bool,
    skeleton_min_link_ratio: bool,
    skeleton_rescue_link_support: bool,
}

#[derive(Debug, Clone)]
struct Config {
    reads: Vec<PathBuf>,
    out_dir: PathBuf,
    organelle: Option<OrganelleProfile>,
    data_mode: DataMode,
    numt_interference: NumtInterference,
    k: usize,
    s: usize,
    syncmer_pos: Option<usize>,
    min_anchor_spacing: usize,
    min_anchor_coverage: u32,
    min_edge_coverage: u32,
    min_branch_ratio: f64,
    max_edges_per_state: usize,
    dedup_kmer: usize,
    containment_ratio: f64,
    min_unitig_len: usize,
    min_tip_len: usize,
    min_link_support: u32,
    min_link_ratio: f64,
    read_junction_links: bool,
    bidirectional_links: bool,
    junction_rescue_support: u32,
    min_read_len: usize,
    max_reads: Option<usize>,
    genome_size: Option<u64>,
    asm_coverage: Option<f64>,
    min_overlap: Option<usize>,
    iterations: Option<usize>,
    hifi_error_rate: f64,
    minimap_min_identity: f64,
    minimap_min_align_len: usize,
    paf_max_link_gap: isize,
    skeleton_gfa: Option<PathBuf>,
    skeleton_only: bool,
    skeleton_end_slop: usize,
    skeleton_min_link_support: u32,
    skeleton_min_link_ratio: f64,
    skeleton_rescue_gfa: Option<PathBuf>,
    skeleton_rescue_link_support: u32,
    rounds: usize,
    read_subsets: Vec<ReadSubset>,
    read_subsets_requested: bool,
    read_subsets_auto_recommended: bool,
    read_subset_default_genome_size: Option<u64>,
    read_subset_input_bases: Option<u64>,
    read_subset_target_depth: Option<f64>,
    read_subset_percent: ReadSubset,
    keep_debug_files: bool,
    threads: usize,
    run_minimap2: bool,
    mito_stable: bool,
    mito_stable_selection_mode: MitoStableSelectionMode,
    mito_stable_selected_nodes: HashSet<String>,
    mito_stable_allow_three_way: HashSet<String>,
    mito_stable_forbid_three_way: HashSet<String>,
    mito_stable_force_links: Vec<SkeletonLinkKey>,
    mito_stable_drop_links: Vec<SkeletonLinkKey>,
}

const READ_SUBSET_SCALE: u16 = 10_000;
const HIGH_PROFILE_TARGET_DEPTH: f64 = 300.0;
const HIGH_PROFILE_MITO_DEFAULT_GENOME_SIZE: u64 = 500_000;
const HIGH_PROFILE_PLASTID_DEFAULT_GENOME_SIZE: u64 = 150_000;
const DEFAULT_THREADS: usize = 8;

// OrgRAFT-facing API. The public command resolves user profiles before calling
// the core, so the engine receives algorithm-level switches.
#[derive(Debug, Clone)]
pub(super) struct DraftAssemblyRequest {
    pub(super) organelle: DraftOrganelle,
    pub(super) data_mode: DraftDataMode,
    pub(super) auto_read_subset: bool,
    pub(super) repeat_aware_resolution: bool,
    pub(super) reads: Vec<PathBuf>,
    pub(super) out_dir: PathBuf,
    pub(super) threads: usize,
    pub(super) min_graph_coverage: Option<u32>,
    pub(super) min_branch_ratio: Option<f64>,
    pub(super) min_tip_len: Option<usize>,
    pub(super) min_link_support: Option<u32>,
    pub(super) min_link_ratio: Option<f64>,
    pub(super) read_subsets: Option<Vec<u16>>,
    pub(super) keep_debug_files: bool,
}

#[derive(Debug, Clone, Copy, Eq, PartialEq)]
pub(super) enum DraftOrganelle {
    Mito,
    Plastid,
}

#[derive(Debug, Clone, Copy, Eq, PartialEq)]
pub(super) enum DraftDataMode {
    Standard,
    Low,
}

pub(super) fn run_draft_assembly(request: DraftAssemblyRequest) -> io::Result<()> {
    let started = Instant::now();
    let config = Config::from_draft_request(request)
        .map_err(|message| io::Error::new(io::ErrorKind::InvalidInput, message))?;
    run_algorithm_steps(config, started)
}

#[derive(Debug, Clone, Copy, Eq, PartialEq, Ord, PartialOrd)]
struct ReadSubset {
    basis_points: u16,
}

impl ReadSubset {
    fn full() -> Self {
        Self {
            basis_points: READ_SUBSET_SCALE,
        }
    }

    fn is_full(self) -> bool {
        self.basis_points >= READ_SUBSET_SCALE
    }

    fn dir_label(self) -> String {
        self.to_string().replace('.', "_")
    }

    fn includes_bucket(self, bucket: u16) -> bool {
        self.basis_points >= READ_SUBSET_SCALE || bucket < self.basis_points
    }
}

impl std::fmt::Display for ReadSubset {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let whole = self.basis_points / 100;
        let decimal = self.basis_points % 100;
        if decimal == 0 {
            write!(f, "{whole}")
        } else if decimal % 10 == 0 {
            write!(f, "{}.{}", whole, decimal / 10)
        } else {
            write!(f, "{}.{:02}", whole, decimal)
        }
    }
}

#[derive(Debug, Clone)]
struct AnchorNode {
    seq: String,
    coverage: u32,
}

#[derive(Debug, Clone, Copy, Eq, PartialEq, Ord, PartialOrd, Hash)]
struct EdgeKey {
    from: usize,
    to: usize,
}

#[derive(Debug, Clone)]
struct EdgeStats {
    coverage: u32,
    total_span: u64,
    min_span: usize,
    max_span: usize,
    suffix: String,
}

#[derive(Debug, Clone, Copy)]
struct AnchorHit {
    state: usize,
    pos: usize,
}

#[derive(Debug)]
struct ReadWalk {
    edges: Vec<EdgeKey>,
}

#[derive(Debug)]
struct Assembly {
    nodes: Vec<AnchorNode>,
    key_to_node: HashMap<String, usize>,
    edges: HashMap<EdgeKey, EdgeStats>,
    edge_junctions: HashMap<(EdgeKey, EdgeKey), u32>,
    read_walks: Vec<ReadWalk>,
    reads_seen: usize,
    bases_seen: u64,
    anchors_seen: u64,
}

#[derive(Debug, Clone)]
struct Unitig {
    id: usize,
    path_states: Vec<usize>,
    path_edges: Vec<EdgeKey>,
    sequence: String,
    coverage: f64,
}

struct AnchorGraphResult {
    assembly: Assembly,
    edge_junction_support: HashMap<EdgeKey, u32>,
    graph: FilteredGraph,
    compressed: CompressedGraph,
    junctions: HashMap<(usize, usize), u32>,
    link_support: HashMap<(UnitigEnd, UnitigEnd), u32>,
    full_graph: FilteredGraph,
    full_compressed: CompressedGraph,
    full_junctions: HashMap<(usize, usize), u32>,
    full_link_support: HashMap<(UnitigEnd, UnitigEnd), u32>,
}

struct ReadSubsetRun {
    percent: ReadSubset,
    config: Config,
    assembly: Assembly,
    target_bases: Option<u64>,
    read_id_path: Option<PathBuf>,
    read_id_writer: Option<BufWriter<File>>,
    subset_read_path: Option<PathBuf>,
    subset_read_writer: Option<BufWriter<File>>,
    selected_reads: usize,
    selected_bases: u64,
    started: Instant,
}

enum TextReader {
    Plain(BufReader<File>),
    Gzip {
        child: Child,
        reader: BufReader<ChildStdout>,
    },
}

impl Read for TextReader {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        match self {
            TextReader::Plain(reader) => reader.read(buf),
            TextReader::Gzip { reader, .. } => reader.read(buf),
        }
    }
}

impl BufRead for TextReader {
    fn fill_buf(&mut self) -> io::Result<&[u8]> {
        match self {
            TextReader::Plain(reader) => reader.fill_buf(),
            TextReader::Gzip { reader, .. } => reader.fill_buf(),
        }
    }

    fn consume(&mut self, amt: usize) {
        match self {
            TextReader::Plain(reader) => reader.consume(amt),
            TextReader::Gzip { reader, .. } => reader.consume(amt),
        }
    }
}

impl Drop for TextReader {
    fn drop(&mut self) {
        if let TextReader::Gzip { child, .. } = self {
            let _ = child.kill();
            let _ = child.wait();
        }
    }
}

// Workflow router: choose the direct anchor graph path, read-subset path,
// skeleton-only path, or skeleton-link workflow from resolved config.
fn run_algorithm_steps(mut config: Config, started: Instant) -> io::Result<()> {
    validate_config(&config)?;
    fs::create_dir_all(&config.out_dir)?;

    if config.read_subsets_requested {
        if config.read_subsets_auto_recommended {
            apply_high_profile_auto_subset(&mut config)?;
            validate_config(&config)?;
        }
        run_read_subset_workflows(config)?;
        return Ok(());
    }

    run_configured_workflow(&config, started)
}

fn run_configured_workflow(config: &Config, started: Instant) -> io::Result<()> {
    if config.mito_stable && config.skeleton_gfa.is_some() {
        run_mito_stable_linking(config, None)?;
        return Ok(());
    }
    if config.skeleton_only {
        if config.mito_stable {
            run_mito_stable_linking(config, None)?;
            return Ok(());
        }
        run_skeleton_linking(config)?;
        return Ok(());
    }

    if config.rounds >= 2 {
        run_two_round_skeleton_workflow(config, started)?;
        return Ok(());
    }

    run_anchor_graph_workflow(config, started)
}

fn is_mito_compact(config: &Config) -> bool {
    config.organelle == Some(OrganelleProfile::Mito) && config.data_mode == DataMode::Low
}

fn run_read_subset_workflows(config: Config) -> io::Result<()> {
    if read_subsets_need_materialized_reads(&config) {
        return run_materialized_read_subsets(config);
    }
    run_anchor_read_subsets(config)
}

fn read_subsets_need_materialized_reads(config: &Config) -> bool {
    config.rounds >= 2
        || config.skeleton_only
        || config.skeleton_gfa.is_some()
        || config.skeleton_rescue_gfa.is_some()
        || config.run_minimap2
}

fn read_subsets_from_basis_points(mut values: Vec<u16>) -> Result<Vec<ReadSubset>, String> {
    values.sort_unstable();
    values.dedup();
    if values.is_empty()
        || values
            .iter()
            .any(|value| *value == 0 || *value > READ_SUBSET_SCALE)
    {
        return Err("read subsets must be between 1 and 100".to_string());
    }
    Ok(values
        .into_iter()
        .map(|basis_points| ReadSubset { basis_points })
        .collect())
}

fn apply_high_profile_auto_subset(config: &mut Config) -> io::Result<()> {
    let genome_size = high_profile_default_genome_size(config).ok_or_else(|| {
        io::Error::new(
            io::ErrorKind::InvalidInput,
            "high profile requires --organelle mito or plastid",
        )
    })?;
    let input_bases = count_input_sequence_bases(config)?;
    if input_bases == 0 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "cannot recommend high-profile subset because input read bases are zero",
        ));
    }

    let subset = recommend_high_profile_subset(genome_size, input_bases);
    config.read_subsets = vec![subset];
    config.read_subset_default_genome_size = Some(genome_size);
    config.read_subset_input_bases = Some(input_bases);
    config.read_subset_target_depth = Some(HIGH_PROFILE_TARGET_DEPTH);
    Ok(())
}

fn high_profile_default_genome_size(config: &Config) -> Option<u64> {
    match config.organelle {
        Some(OrganelleProfile::Mito) => Some(HIGH_PROFILE_MITO_DEFAULT_GENOME_SIZE),
        Some(OrganelleProfile::Plastid) => Some(HIGH_PROFILE_PLASTID_DEFAULT_GENOME_SIZE),
        None => None,
    }
}

fn count_input_sequence_bases(config: &Config) -> io::Result<u64> {
    let mut bases = 0u64;
    for path in &config.reads {
        read_sequence_file(path, |_name, seq| {
            if seq.len() >= config.min_read_len {
                bases = bases.saturating_add(seq.len() as u64);
            }
            Ok(true)
        })?;
    }
    Ok(bases)
}

fn recommend_high_profile_subset(genome_size: u64, input_bases: u64) -> ReadSubset {
    let target_bases = genome_size as f64 * HIGH_PROFILE_TARGET_DEPTH;
    let basis_points = ((target_bases / input_bases as f64) * READ_SUBSET_SCALE as f64)
        .round()
        .clamp(1.0, READ_SUBSET_SCALE as f64) as u16;
    ReadSubset { basis_points }
}

fn write_read_subset_recommendation_summary<W: Write>(
    out: &mut W,
    config: &Config,
) -> io::Result<()> {
    if !config.read_subsets_auto_recommended {
        writeln!(out, "# subset_recommendation\tmanual")?;
        return Ok(());
    }
    writeln!(out, "# subset_recommendation\tauto")?;
    if let Some(genome_size) = config.read_subset_default_genome_size {
        writeln!(out, "# default_genome_size\t{genome_size}")?;
    }
    if let Some(input_bases) = config.read_subset_input_bases {
        writeln!(out, "# input_read_bases\t{input_bases}")?;
    }
    if let Some(target_depth) = config.read_subset_target_depth {
        writeln!(out, "# subset_target_depth\t{target_depth:.0}")?;
    }
    writeln!(
        out,
        "# recommendation_formula\tsubset_target_depth * default_genome_size / input_read_bases"
    )?;
    Ok(())
}

fn run_anchor_read_subsets(config: Config) -> io::Result<()> {
    fs::create_dir_all(&config.out_dir)?;

    let mut runs = new_read_subset_runs(&config, false)?;

    let mut reads_loaded = 0u64;
    'inputs: for path in &config.reads {
        read_sequence_file(path, |name, seq| {
            let read_index = reads_loaded;
            reads_loaded += 1;

            if seq.len() >= config.min_read_len {
                let bucket = read_subset_bucket(name, read_index);
                for run in &mut runs {
                    if run.percent.includes_bucket(bucket)
                        && !assembly_limit_reached(&run.config, &run.assembly, run.target_bases)
                    {
                        process_read(&run.config, &mut run.assembly, name, seq)?;
                        if let Some(writer) = run.read_id_writer.as_mut() {
                            writeln!(writer, "{name}")?;
                        }
                        run.selected_reads += 1;
                        run.selected_bases += seq.len() as u64;
                    }
                }
            }

            Ok(!all_read_subset_runs_finished(&runs))
        })?;

        if all_read_subset_runs_finished(&runs) {
            break 'inputs;
        }
    }

    let selection_dir = config.out_dir.join("01.read_subset_selection");
    fs::create_dir_all(&selection_dir)?;
    let summary_path = selection_dir.join("summary.tsv");
    let mut summary = File::create(&summary_path)?;
    write_read_subset_recommendation_summary(&mut summary, &config)?;
    writeln!(
        summary,
        "read_subset_percent\tout_dir\tread_ids_file\telapsed_seconds"
    )?;

    for mut run in runs {
        if let Some(writer) = run.read_id_writer.as_mut() {
            writer.flush()?;
        }
        let assembly = std::mem::replace(&mut run.assembly, new_assembly());
        let result = finish_anchor_graph(&run.config, assembly);
        let support_dir = anchor_support_dir(&run.config.out_dir);
        let graph_dir = unitig_graph_dir(&run.config.out_dir);
        let read_junction_dir = read_junction_graph_dir(&run.config.out_dir);
        let bridge_rescue_dir = low_depth_bridge_rescue_dir(&run.config.out_dir);
        let link_evidence_dir = skeleton_link_evidence_dir(&run.config.out_dir);
        let repeat_resolution_dir = repeat_aware_resolution_dir(&run.config.out_dir);
        let linked_graph_dir = linked_graph_dir(&run.config.out_dir);
        let workflow_dir = workflow_summary_dir(&run.config.out_dir);
        write_anchor_graph_outputs(
            &run.config,
            &result,
            run.started,
            Some(&support_dir),
            &graph_dir,
        )?;
        write_step_status_report(
            &read_junction_dir,
            "03.read_junction_graph",
            "skipped",
            "direct-anchor workflow does not build a read-junction companion graph",
            Some(&graph_dir),
            None,
        )?;
        write_step_status_report(
            &bridge_rescue_dir,
            "04.low_depth_bridge_rescue",
            "skipped",
            "low-depth bridge rescue is not required for direct-anchor workflow",
            Some(&graph_dir),
            None,
        )?;
        write_step_status_report(
            &link_evidence_dir,
            "05.skeleton_link_evidence",
            "skipped",
            "direct-anchor workflow does not remap reads to a skeleton graph",
            Some(&graph_dir),
            None,
        )?;
        write_step_status_report(
            &repeat_resolution_dir,
            "06.repeat_aware_resolution",
            "skipped",
            "repeat-aware resolution is not required for direct-anchor workflow",
            Some(&graph_dir),
            None,
        )?;
        copy_one_round_linked_graph(&graph_dir, &linked_graph_dir)?;
        write_one_round_report(
            &run.config,
            &support_dir,
            &graph_dir,
            &read_junction_dir,
            &bridge_rescue_dir,
            &link_evidence_dir,
            &repeat_resolution_dir,
            &linked_graph_dir,
            &workflow_dir,
            run.started,
        )?;
        writeln!(
            summary,
            "{}\t{}\t{}\t{:.3}",
            run.percent,
            run.config.out_dir.display(),
            run.read_id_path
                .as_ref()
                .map(|path| path.display().to_string())
                .unwrap_or_default(),
            run.started.elapsed().as_secs_f64()
        )?;
    }
    drop(summary);
    copy_selected_read_subset_output(&config, &summary_path)?;
    write_profile_parameter_report(
        &workflow_summary_dir(&config.out_dir).join("profile_parameters.tsv"),
        &config,
    )?;

    Ok(())
}

fn run_materialized_read_subsets(config: Config) -> io::Result<()> {
    fs::create_dir_all(&config.out_dir)?;
    let mut runs = new_read_subset_runs(&config, true)?;

    materialize_read_subsets(&config, &mut runs)?;
    for run in &mut runs {
        if let Some(writer) = run.read_id_writer.as_mut() {
            writer.flush()?;
        }
        if let Some(writer) = run.subset_read_writer.as_mut() {
            writer.flush()?;
        }
        run.read_id_writer = None;
        run.subset_read_writer = None;
    }

    let selection_dir = config.out_dir.join("01.read_subset_selection");
    fs::create_dir_all(&selection_dir)?;
    let summary_path = selection_dir.join("summary.tsv");
    let mut summary = File::create(&summary_path)?;
    write_read_subset_recommendation_summary(&mut summary, &config)?;
    writeln!(
        summary,
        "read_subset_percent\tout_dir\tread_ids_file\tsubset_reads_file\telapsed_seconds"
    )?;

    for run in &runs {
        run_configured_workflow(&run.config, run.started)?;
        let subset_read_status = match &run.subset_read_path {
            Some(path) if config.keep_debug_files => path.display().to_string(),
            Some(path) => match fs::remove_file(path) {
                Ok(()) => "deleted_after_run".to_string(),
                Err(error) if error.kind() == io::ErrorKind::NotFound => {
                    "deleted_after_run".to_string()
                }
                Err(error) => return Err(error),
            },
            None => String::new(),
        };
        writeln!(
            summary,
            "{}\t{}\t{}\t{}\t{:.3}",
            run.percent,
            run.config.out_dir.display(),
            run.read_id_path
                .as_ref()
                .map(|path| path.display().to_string())
                .unwrap_or_default(),
            subset_read_status,
            run.started.elapsed().as_secs_f64()
        )?;
    }
    drop(summary);
    copy_selected_read_subset_output(&config, &summary_path)?;
    write_profile_parameter_report(
        &workflow_summary_dir(&config.out_dir).join("profile_parameters.tsv"),
        &config,
    )?;

    Ok(())
}

fn new_read_subset_runs(
    config: &Config,
    materialize_reads: bool,
) -> io::Result<Vec<ReadSubsetRun>> {
    let mut runs = Vec::new();
    for subset in &config.read_subsets {
        let mut subset_config = config.clone();
        subset_config.out_dir = config
            .out_dir
            .join(format!("read_subset_{}", subset.dir_label()));
        subset_config.read_subsets_requested = false;
        subset_config.read_subsets = vec![*subset];
        subset_config.read_subset_percent = *subset;

        fs::create_dir_all(&subset_config.out_dir)?;
        let read_id_path = if subset.is_full() {
            None
        } else {
            Some(subset_config.out_dir.join("read_ids.txt"))
        };
        let read_id_writer = match &read_id_path {
            Some(path) => Some(BufWriter::new(File::create(path)?)),
            None => None,
        };
        let subset_read_path = if materialize_reads && !subset.is_full() {
            Some(subset_config.out_dir.join("reads.fasta"))
        } else {
            None
        };
        let subset_read_writer = match &subset_read_path {
            Some(path) => Some(BufWriter::new(File::create(path)?)),
            None => None,
        };
        if let Some(path) = &subset_read_path {
            subset_config.reads = vec![path.clone()];
        }

        runs.push(ReadSubsetRun {
            percent: *subset,
            target_bases: target_bases_for_config(&subset_config),
            config: subset_config,
            assembly: new_assembly(),
            read_id_path,
            read_id_writer,
            subset_read_path,
            subset_read_writer,
            selected_reads: 0,
            selected_bases: 0,
            started: Instant::now(),
        });
    }
    Ok(runs)
}

fn materialize_read_subsets(config: &Config, runs: &mut [ReadSubsetRun]) -> io::Result<()> {
    let mut reads_loaded = 0u64;
    'inputs: for path in &config.reads {
        read_sequence_file(path, |name, seq| {
            let read_index = reads_loaded;
            reads_loaded += 1;

            if seq.len() >= config.min_read_len {
                let bucket = read_subset_bucket(name, read_index);
                for run in runs.iter_mut() {
                    if run.percent.includes_bucket(bucket)
                        && !read_subset_selection_limit_reached(run)
                    {
                        if let Some(writer) = run.subset_read_writer.as_mut() {
                            write_fasta_record(writer, name, read_index, seq)?;
                        }
                        if let Some(writer) = run.read_id_writer.as_mut() {
                            writeln!(writer, "{}", read_output_name(name, read_index))?;
                        }
                        run.selected_reads += 1;
                        run.selected_bases += seq.len() as u64;
                    }
                }
            }

            Ok(!all_materialized_read_subset_runs_finished(runs))
        })?;

        if all_materialized_read_subset_runs_finished(runs) {
            break 'inputs;
        }
    }
    Ok(())
}

fn read_subset_selection_limit_reached(run: &ReadSubsetRun) -> bool {
    if let Some(max_reads) = run.config.max_reads {
        if run.selected_reads >= max_reads {
            return true;
        }
    }
    if let Some(target) = run.target_bases {
        if run.selected_bases >= target {
            return true;
        }
    }
    false
}

fn all_materialized_read_subset_runs_finished(runs: &[ReadSubsetRun]) -> bool {
    runs.iter().all(read_subset_selection_limit_reached)
}

fn write_fasta_record<W: Write>(
    writer: &mut W,
    name: &str,
    read_index: u64,
    seq: &str,
) -> io::Result<()> {
    writeln!(writer, ">{}", read_output_name(name, read_index))?;
    for chunk in seq.as_bytes().chunks(80) {
        writer.write_all(chunk)?;
        writer.write_all(b"\n")?;
    }
    Ok(())
}

fn read_output_name(name: &str, read_index: u64) -> String {
    if name.is_empty() {
        format!("read_{read_index}")
    } else {
        name.to_string()
    }
}

// Step 02 core: convert recruited reads into an anchor graph and compressed
// unitigs. The caller owns Step 01 input staging and Step 03 final graph copy.
fn run_anchor_graph_workflow(config: &Config, started: Instant) -> io::Result<()> {
    let result = compute_anchor_graph(config)?;
    let support_dir = anchor_support_dir(&config.out_dir);
    let graph_dir = unitig_graph_dir(&config.out_dir);
    let read_junction_dir = read_junction_graph_dir(&config.out_dir);
    let bridge_rescue_dir = low_depth_bridge_rescue_dir(&config.out_dir);
    let link_evidence_dir = skeleton_link_evidence_dir(&config.out_dir);
    let repeat_resolution_dir = repeat_aware_resolution_dir(&config.out_dir);
    let linked_graph_dir = linked_graph_dir(&config.out_dir);
    let summary_dir = workflow_summary_dir(&config.out_dir);
    write_anchor_graph_outputs(config, &result, started, Some(&support_dir), &graph_dir)?;

    if config.run_minimap2 {
        let mut graph_config = config.clone();
        graph_config.out_dir = graph_dir.clone();
        run_minimap2(&graph_config)?;
    }

    if config.skeleton_gfa.is_some() {
        run_skeleton_linking(config)?;
    }

    write_step_status_report(
        &read_junction_dir,
        "03.read_junction_graph",
        "skipped",
        "direct-anchor workflow does not build a read-junction companion graph",
        Some(&graph_dir),
        None,
    )?;
    write_step_status_report(
        &bridge_rescue_dir,
        "04.low_depth_bridge_rescue",
        "skipped",
        "low-depth bridge rescue is not required for direct-anchor workflow",
        Some(&graph_dir),
        None,
    )?;
    write_step_status_report(
        &link_evidence_dir,
        "05.skeleton_link_evidence",
        "skipped",
        "direct-anchor workflow does not remap reads to a skeleton graph",
        Some(&graph_dir),
        None,
    )?;
    write_step_status_report(
        &repeat_resolution_dir,
        "06.repeat_aware_resolution",
        "skipped",
        "repeat-aware resolution is not required for direct-anchor workflow",
        Some(&graph_dir),
        None,
    )?;
    copy_one_round_linked_graph(&graph_dir, &linked_graph_dir)?;

    write_one_round_report(
        config,
        &support_dir,
        &graph_dir,
        &read_junction_dir,
        &bridge_rescue_dir,
        &link_evidence_dir,
        &repeat_resolution_dir,
        &linked_graph_dir,
        &summary_dir,
        started,
    )?;
    Ok(())
}

fn compute_anchor_graph(config: &Config) -> io::Result<AnchorGraphResult> {
    let mut assembly = new_assembly();
    let target_bases = target_bases_for_config(config);

    'inputs: for path in &config.reads {
        read_sequence_file(path, |name, seq| {
            if seq.len() < config.min_read_len {
                return Ok(true);
            }
            if assembly_limit_reached(config, &assembly, target_bases) {
                return Ok(false);
            }

            process_read(&config, &mut assembly, name, seq)?;
            Ok(true)
        })?;

        if assembly_limit_reached(config, &assembly, target_bases) {
            break 'inputs;
        }
    }

    Ok(finish_anchor_graph(config, assembly))
}

fn new_assembly() -> Assembly {
    Assembly {
        nodes: Vec::new(),
        key_to_node: HashMap::new(),
        edges: HashMap::new(),
        edge_junctions: HashMap::new(),
        read_walks: Vec::new(),
        reads_seen: 0,
        bases_seen: 0,
        anchors_seen: 0,
    }
}

fn target_bases_for_config(config: &Config) -> Option<u64> {
    match (config.genome_size, config.asm_coverage) {
        (Some(genome_size), Some(coverage)) if coverage > 0.0 => {
            Some((genome_size as f64 * coverage).round() as u64)
        }
        _ => None,
    }
}

fn assembly_limit_reached(config: &Config, assembly: &Assembly, target_bases: Option<u64>) -> bool {
    if let Some(max_reads) = config.max_reads {
        if assembly.reads_seen >= max_reads {
            return true;
        }
    }
    if let Some(target) = target_bases {
        if assembly.bases_seen >= target {
            return true;
        }
    }
    false
}

fn all_read_subset_runs_finished(runs: &[ReadSubsetRun]) -> bool {
    runs.iter()
        .all(|run| assembly_limit_reached(&run.config, &run.assembly, run.target_bases))
}

fn finish_anchor_graph(config: &Config, assembly: Assembly) -> AnchorGraphResult {
    let edge_junction_support = edge_junction_support(&assembly);
    let graph = build_filtered_graph(&config, &assembly, &edge_junction_support);
    let compressed = compress_unitigs(&config, &assembly, &graph);
    let junctions = count_unitig_junctions(&assembly, &compressed.edge_to_unitig);
    let link_support = count_link_support(&assembly, &compressed.edge_to_placement);

    let full_config = full_graph_config(&config);
    let full_graph = build_filtered_graph(&full_config, &assembly, &edge_junction_support);
    let full_compressed = compress_unitigs(&full_config, &assembly, &full_graph);
    let full_junctions = count_unitig_junctions(&assembly, &full_compressed.edge_to_unitig);
    let full_link_support = count_link_support(&assembly, &full_compressed.edge_to_placement);

    AnchorGraphResult {
        assembly,
        edge_junction_support,
        graph,
        compressed,
        junctions,
        link_support,
        full_graph,
        full_compressed,
        full_junctions,
        full_link_support,
    }
}

fn read_subset_bucket(name: &str, read_index: u64) -> u16 {
    let hash = if name.is_empty() {
        mix64(read_index.wrapping_add(0x9e37_79b9_7f4a_7c15))
    } else {
        stable_hash_bytes_seed(name.as_bytes(), 0x9e37_79b9_7f4a_7c15, 0x1000_0000_01b3)
    };
    (mix64(hash) % READ_SUBSET_SCALE as u64) as u16
}

fn anchor_support_dir(root: &Path) -> PathBuf {
    root.join("01.anchor_walk_support")
}

fn unitig_graph_dir(root: &Path) -> PathBuf {
    root.join("02.unitig_graph")
}

fn read_junction_graph_dir(root: &Path) -> PathBuf {
    root.join("03.read_junction_graph")
}

fn low_depth_bridge_rescue_dir(root: &Path) -> PathBuf {
    root.join("04.low_depth_bridge_rescue")
}

fn skeleton_link_evidence_dir(root: &Path) -> PathBuf {
    root.join("05.skeleton_link_evidence")
}

fn repeat_aware_resolution_dir(root: &Path) -> PathBuf {
    root.join("06.repeat_aware_resolution")
}

fn linked_graph_dir(root: &Path) -> PathBuf {
    root.join("07.linked_graph")
}

fn workflow_summary_dir(root: &Path) -> PathBuf {
    root.join("08.workflow_summary")
}

fn write_anchor_graph_outputs(
    config: &Config,
    result: &AnchorGraphResult,
    started: Instant,
    support_dir: Option<&Path>,
    graph_dir: &Path,
) -> io::Result<()> {
    if let Some(support_dir) = support_dir {
        fs::create_dir_all(support_dir)?;
        let mut support_config = config.clone();
        support_config.out_dir = support_dir.to_path_buf();
        write_anchors(&support_config, &result.assembly, &result.graph)?;
        write_edges(
            &support_config,
            &result.graph,
            &result.edge_junction_support,
            &result.compressed.edge_to_unitig,
        )?;
    }

    fs::create_dir_all(graph_dir)?;
    let mut graph_config = config.clone();
    graph_config.out_dir = graph_dir.to_path_buf();
    write_unitigs(&graph_config, &result.compressed.unitigs, "unitigs.fasta")?;
    if graph_config.keep_debug_files {
        write_unitigs(
            &graph_config,
            &result.full_compressed.unitigs,
            "unitigs.full.fasta",
        )?;
    }
    write_gfa(
        &graph_config,
        &result.compressed.unitigs,
        &result.compressed.state_starts,
        &result.compressed.state_ends,
        &result.link_support,
        "graph.gfa",
    )?;
    if graph_config.keep_debug_files {
        write_gfa(
            &graph_config,
            &result.full_compressed.unitigs,
            &result.full_compressed.state_starts,
            &result.full_compressed.state_ends,
            &result.full_link_support,
            "graph.full.gfa",
        )?;
    }
    write_junctions(&graph_config, &result.junctions, "junctions.tsv")?;
    if graph_config.keep_debug_files {
        write_junctions(&graph_config, &result.full_junctions, "junctions.full.tsv")?;
    }
    write_report(
        &graph_config,
        &result.assembly,
        &result.graph,
        &result.full_graph,
        &result.compressed.unitigs,
        &result.full_compressed.unitigs,
        started.elapsed(),
    )?;
    Ok(())
}

fn run_two_round_skeleton_workflow(config: &Config, started: Instant) -> io::Result<()> {
    fs::create_dir_all(&config.out_dir)?;

    let support_dir = anchor_support_dir(&config.out_dir);
    let seed_graph_dir = unitig_graph_dir(&config.out_dir);
    let read_junction_dir = read_junction_graph_dir(&config.out_dir);
    let bridge_rescue_dir = low_depth_bridge_rescue_dir(&config.out_dir);
    let link_evidence_dir = skeleton_link_evidence_dir(&config.out_dir);
    let repeat_resolution_dir = repeat_aware_resolution_dir(&config.out_dir);
    let linked_graph_dir = linked_graph_dir(&config.out_dir);
    let summary_dir = workflow_summary_dir(&config.out_dir);

    let mut round1 = config.clone();
    round1.out_dir = seed_graph_dir.clone();
    round1.rounds = 1;
    round1.skeleton_gfa = None;
    round1.skeleton_rescue_gfa = None;
    round1.skeleton_only = false;
    round1.read_junction_links = false;
    round1.min_link_support = round1.min_link_support.min(10);
    round1.run_minimap2 = false;
    let round1_result = compute_anchor_graph(&round1)?;
    write_anchor_graph_outputs(
        &round1,
        &round1_result,
        started,
        Some(&support_dir),
        &seed_graph_dir,
    )?;

    let mut read_junction_graph_config = round1.clone();
    read_junction_graph_config.out_dir = read_junction_dir.clone();
    read_junction_graph_config.read_junction_links = true;
    write_anchor_graph_outputs(
        &read_junction_graph_config,
        &round1_result,
        started,
        None,
        &read_junction_dir,
    )?;

    let mut round2 = config.clone();
    round2.out_dir = link_evidence_dir.clone();
    round2.rounds = 1;
    round2.skeleton_only = true;
    round2.skeleton_gfa = Some(seed_graph_dir.join("graph.gfa"));
    round2.skeleton_rescue_gfa = Some(read_junction_dir.join("graph.gfa"));
    let final_graph_source_dir;
    if config.mito_stable {
        run_mito_stable_linking(&round2, Some(&repeat_resolution_dir))?;
        write_low_depth_bridge_rescue_report(
            config,
            &bridge_rescue_dir,
            "skipped",
            "low-depth bridge rescue is not required for repeat-aware mito workflow",
            &read_junction_dir,
            &link_evidence_dir,
            false,
        )?;
        write_step_status_report(
            &repeat_resolution_dir,
            "06.repeat_aware_resolution",
            "resolved",
            "repeat-aware mito workflow resolved candidate links and topology repairs",
            Some(&link_evidence_dir),
            Some(&repeat_resolution_dir),
        )?;
        final_graph_source_dir = repeat_resolution_dir.as_path();
    } else if is_mito_compact(config) {
        run_mito_compact_bridge_linking(&round2)?;
        move_low_depth_bridge_rescue_outputs(&link_evidence_dir, &bridge_rescue_dir)?;
        write_low_depth_bridge_rescue_report(
            config,
            &bridge_rescue_dir,
            "resolved",
            "low-depth mito workflow used read-junction rescue evidence during bridge selection",
            &read_junction_dir,
            &link_evidence_dir,
            true,
        )?;
        write_step_status_report(
            &repeat_resolution_dir,
            "06.repeat_aware_resolution",
            "skipped",
            "repeat-aware resolution is not required for low-depth mito workflow",
            Some(&link_evidence_dir),
            None,
        )?;
        final_graph_source_dir = link_evidence_dir.as_path();
    } else {
        run_skeleton_linking(&round2)?;
        write_low_depth_bridge_rescue_report(
            config,
            &bridge_rescue_dir,
            "skipped",
            "low-depth bridge rescue is only enabled for mito low data mode",
            &read_junction_dir,
            &link_evidence_dir,
            false,
        )?;
        write_step_status_report(
            &repeat_resolution_dir,
            "06.repeat_aware_resolution",
            "skipped",
            "repeat-aware resolution is not required for this skeleton-link workflow",
            Some(&link_evidence_dir),
            None,
        )?;
        final_graph_source_dir = link_evidence_dir.as_path();
    }

    copy_final_two_round_outputs(final_graph_source_dir, &linked_graph_dir)?;
    write_two_round_report(
        config,
        &support_dir,
        &seed_graph_dir,
        &read_junction_dir,
        &bridge_rescue_dir,
        &link_evidence_dir,
        &repeat_resolution_dir,
        if config.mito_stable {
            "resolved"
        } else {
            "skipped"
        },
        &linked_graph_dir,
        &summary_dir,
        started,
    )?;
    Ok(())
}

fn copy_final_two_round_outputs(evidence_dir: &Path, linked_graph_dir: &Path) -> io::Result<()> {
    fs::create_dir_all(linked_graph_dir)?;
    copy_if_exists(
        evidence_dir.join("graph.gfa"),
        linked_graph_dir.join("graph.gfa"),
    )?;
    copy_if_exists(
        evidence_dir.join("links.tsv"),
        linked_graph_dir.join("links.tsv"),
    )?;
    copy_if_exists(
        evidence_dir.join("depth.tsv"),
        linked_graph_dir.join("depth.tsv"),
    )?;
    copy_if_exists(
        evidence_dir.join("linking.report.txt"),
        linked_graph_dir.join("report.txt"),
    )?;
    Ok(())
}

fn copy_one_round_linked_graph(graph_dir: &Path, linked_graph_dir: &Path) -> io::Result<()> {
    fs::create_dir_all(linked_graph_dir)?;
    fs::copy(
        graph_dir.join("graph.gfa"),
        linked_graph_dir.join("graph.gfa"),
    )?;
    let mut out = File::create(linked_graph_dir.join("report.txt"))?;
    writeln!(out, "linked_graph_report")?;
    writeln!(out, "status\tselected")?;
    writeln!(out, "source_stage\t02.unitig_graph")?;
    writeln!(
        out,
        "source_graph\t{}",
        graph_dir.join("graph.gfa").display()
    )?;
    writeln!(
        out,
        "output_graph\t{}",
        linked_graph_dir.join("graph.gfa").display()
    )?;
    writeln!(out, "output_depth\tskipped")?;
    writeln!(out, "output_links\tskipped")?;
    Ok(())
}

const READ_SUBSET_TARGET_DEPTH: f64 = HIGH_PROFILE_TARGET_DEPTH;
const READ_SUBSET_MIN_STANDARD_DEPTH: f64 = 200.0;
const READ_SUBSET_MAX_STANDARD_DEPTH: f64 = 450.0;

#[derive(Debug)]
struct ReadSubsetCandidate {
    subset_name: String,
    subset_percent: String,
    graph: PathBuf,
    graph_dir: PathBuf,
    read_bases: Option<u64>,
    estimated_depth: Option<f64>,
    graph_stats: GfaTopologyStats,
    score: ReadSubsetScore,
}

#[derive(Debug, Clone, Copy, Eq, PartialEq, Ord, PartialOrd)]
struct ReadSubsetScore {
    components: usize,
    open_ends: usize,
    endpoint_overload: usize,
    branch_ends: usize,
    cycle_rank: usize,
    depth_penalty: u8,
    depth_distance_milli: u64,
    subset_basis_points: u16,
}

#[derive(Debug, Clone, Copy)]
struct GfaTopologyStats {
    segments: usize,
    links: usize,
    graph_bases: u64,
    components: usize,
    open_ends: usize,
    endpoint_overload: usize,
    branch_ends: usize,
    cycle_rank: usize,
}

fn copy_selected_read_subset_output(config: &Config, summary_path: &Path) -> io::Result<()> {
    let root_dir = &config.out_dir;
    let mut candidates = Vec::new();
    for entry in fs::read_dir(root_dir)? {
        let entry = entry?;
        if !entry.file_type()?.is_dir() {
            continue;
        }
        let name = entry.file_name().to_string_lossy().to_string();
        let Some(label) = name.strip_prefix("read_subset_") else {
            continue;
        };
        if let Some(graph) = selected_graph_in_work_dir(&entry.path()) {
            candidates.push(read_subset_candidate(&entry.path(), &name, label, graph)?);
        }
    }
    candidates.sort_by(read_subset_candidate_cmp);
    let Some(selected) = candidates.first() else {
        return Err(io::Error::new(
            io::ErrorKind::NotFound,
            format!("no read subset graph found under {}", root_dir.display()),
        ));
    };

    let linked_graph_dir = linked_graph_dir(root_dir);
    fs::create_dir_all(&linked_graph_dir)?;
    fs::copy(&selected.graph, linked_graph_dir.join("graph.gfa"))?;
    copy_if_exists(
        selected.graph_dir.join("links.tsv"),
        linked_graph_dir.join("links.tsv"),
    )?;
    copy_if_exists(
        selected.graph_dir.join("depth.tsv"),
        linked_graph_dir.join("depth.tsv"),
    )?;

    let summary_dir = workflow_summary_dir(root_dir);
    fs::create_dir_all(&summary_dir)?;
    copy_if_exists(
        root_dir
            .join(&selected.subset_name)
            .join("08.workflow_summary/profile_parameters.tsv"),
        summary_dir.join("profile_parameters.tsv"),
    )?;
    write_read_subset_selection_table(&summary_dir.join("subset_selection.tsv"), &candidates)?;

    let mut linked_report = File::create(linked_graph_dir.join("report.txt"))?;
    writeln!(linked_report, "linked_graph_report")?;
    writeln!(linked_report, "status\tselected")?;
    writeln!(
        linked_report,
        "selection_strategy\tposterior_depth_topology"
    )?;
    writeln!(
        linked_report,
        "depth_estimation\tread_bases / candidate_graph_bases after assembly"
    )?;
    writeln!(
        linked_report,
        "standard_depth_range\t{READ_SUBSET_MIN_STANDARD_DEPTH:.0}-{READ_SUBSET_MAX_STANDARD_DEPTH:.0}"
    )?;
    write_auto_subset_report_fields(&mut linked_report, config)?;
    writeln!(linked_report, "source_subset\t{}", selected.subset_name)?;
    writeln!(
        linked_report,
        "source_subset_percent\t{}",
        selected.subset_percent
    )?;
    writeln!(
        linked_report,
        "estimated_depth\t{}",
        format_optional_depth(selected.estimated_depth)
    )?;
    writeln!(
        linked_report,
        "topology\tcomponents={},open_ends={},endpoint_overload={},branch_ends={},cycle_rank={}",
        selected.graph_stats.components,
        selected.graph_stats.open_ends,
        selected.graph_stats.endpoint_overload,
        selected.graph_stats.branch_ends,
        selected.graph_stats.cycle_rank
    )?;
    writeln!(linked_report, "source_graph\t{}", selected.graph.display())?;
    writeln!(
        linked_report,
        "output_graph\t{}",
        linked_graph_dir.join("graph.gfa").display()
    )?;

    let mut workflow_report = File::create(summary_dir.join("report.txt"))?;
    writeln!(workflow_report, "orgraft_asm read-subset report")?;
    writeln!(
        workflow_report,
        "read_subset_summary\t{}",
        summary_path.display()
    )?;
    writeln!(
        workflow_report,
        "selection_strategy\tposterior_depth_topology"
    )?;
    writeln!(
        workflow_report,
        "depth_estimation\tread_bases / candidate_graph_bases after assembly"
    )?;
    writeln!(
        workflow_report,
        "standard_depth_range\t{READ_SUBSET_MIN_STANDARD_DEPTH:.0}-{READ_SUBSET_MAX_STANDARD_DEPTH:.0}"
    )?;
    write_auto_subset_report_fields(&mut workflow_report, config)?;
    writeln!(workflow_report, "selected_subset\t{}", selected.subset_name)?;
    writeln!(
        workflow_report,
        "selected_subset_percent\t{}",
        selected.subset_percent
    )?;
    writeln!(
        workflow_report,
        "estimated_depth\t{}",
        format_optional_depth(selected.estimated_depth)
    )?;
    writeln!(
        workflow_report,
        "subset_selection\t{}",
        summary_dir.join("subset_selection.tsv").display()
    )?;
    writeln!(
        workflow_report,
        "selected_graph\t{}",
        selected.graph.display()
    )?;
    writeln!(
        workflow_report,
        "linked_graph\t{}",
        linked_graph_dir.display()
    )?;
    Ok(())
}

fn read_subset_candidate(
    work_dir: &Path,
    subset_name: &str,
    subset_label: &str,
    graph: PathBuf,
) -> io::Result<ReadSubsetCandidate> {
    let graph_stats = gfa_topology_stats(&graph)?;
    let read_bases = read_report_u64(&work_dir.join("02.unitig_graph/report.txt"), "bases_seen")?;
    let estimated_depth = read_bases.and_then(|bases| {
        if graph_stats.graph_bases > 0 {
            Some(bases as f64 / graph_stats.graph_bases as f64)
        } else {
            None
        }
    });
    let subset_basis_points = read_subset_label_to_basis_points(subset_label);
    let depth_penalty = match estimated_depth {
        Some(depth)
            if (READ_SUBSET_MIN_STANDARD_DEPTH..=READ_SUBSET_MAX_STANDARD_DEPTH)
                .contains(&depth) =>
        {
            0
        }
        Some(_) => 1,
        None => 2,
    };
    let depth_distance_milli = estimated_depth
        .map(|depth| ((depth - READ_SUBSET_TARGET_DEPTH).abs() * 1000.0).round() as u64)
        .unwrap_or(u64::MAX);
    let score = ReadSubsetScore {
        components: graph_stats.components,
        open_ends: graph_stats.open_ends,
        endpoint_overload: graph_stats.endpoint_overload,
        branch_ends: graph_stats.branch_ends,
        cycle_rank: graph_stats.cycle_rank,
        depth_penalty,
        depth_distance_milli,
        subset_basis_points,
    };
    let graph_dir = graph
        .parent()
        .map(Path::to_path_buf)
        .unwrap_or_else(|| work_dir.to_path_buf());
    Ok(ReadSubsetCandidate {
        subset_name: subset_name.to_string(),
        subset_percent: subset_label.replace('_', "."),
        graph,
        graph_dir,
        read_bases,
        estimated_depth,
        graph_stats,
        score,
    })
}

fn write_auto_subset_report_fields<W: Write>(out: &mut W, config: &Config) -> io::Result<()> {
    if config.read_subsets_auto_recommended {
        writeln!(out, "subset_recommendation\tauto")?;
        if let Some(genome_size) = config.read_subset_default_genome_size {
            writeln!(out, "default_genome_size\t{genome_size}")?;
        }
        if let Some(input_bases) = config.read_subset_input_bases {
            writeln!(out, "input_read_bases\t{input_bases}")?;
        }
        if let Some(target_depth) = config.read_subset_target_depth {
            writeln!(out, "subset_target_depth\t{target_depth:.0}")?;
        }
        writeln!(
            out,
            "recommendation_formula\tsubset_target_depth * default_genome_size / input_read_bases"
        )?;
    } else {
        writeln!(out, "subset_recommendation\tmanual")?;
    }
    Ok(())
}

fn read_subset_candidate_cmp(left: &ReadSubsetCandidate, right: &ReadSubsetCandidate) -> Ordering {
    left.score
        .cmp(&right.score)
        .then_with(|| left.subset_name.cmp(&right.subset_name))
}

fn write_read_subset_selection_table(
    path: &Path,
    candidates: &[ReadSubsetCandidate],
) -> io::Result<()> {
    let mut out = File::create(path)?;
    writeln!(
        out,
        "rank\tsubset\tsubset_percent\tread_bases\testimated_depth\tsegments\tlinks\tgraph_bases\tcomponents\topen_ends\tendpoint_overload\tbranch_ends\tcycle_rank\tdepth_penalty\tdepth_distance"
    )?;
    for (rank, candidate) in candidates.iter().enumerate() {
        writeln!(
            out,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.3}",
            rank + 1,
            candidate.subset_name,
            candidate.subset_percent,
            candidate
                .read_bases
                .map(|value| value.to_string())
                .unwrap_or_else(|| "unknown".to_string()),
            format_optional_depth(candidate.estimated_depth),
            candidate.graph_stats.segments,
            candidate.graph_stats.links,
            candidate.graph_stats.graph_bases,
            candidate.graph_stats.components,
            candidate.graph_stats.open_ends,
            candidate.graph_stats.endpoint_overload,
            candidate.graph_stats.branch_ends,
            candidate.graph_stats.cycle_rank,
            candidate.score.depth_penalty,
            candidate.score.depth_distance_milli as f64 / 1000.0
        )?;
    }
    Ok(())
}

fn read_subset_label_to_basis_points(label: &str) -> u16 {
    label
        .replace('_', ".")
        .parse::<f64>()
        .ok()
        .map(|percent| {
            (percent * 100.0)
                .round()
                .clamp(1.0, READ_SUBSET_SCALE as f64) as u16
        })
        .unwrap_or(READ_SUBSET_SCALE)
}

fn format_optional_depth(depth: Option<f64>) -> String {
    depth
        .map(|value| format!("{value:.3}"))
        .unwrap_or_else(|| "unknown".to_string())
}

fn read_report_u64(path: &Path, key: &str) -> io::Result<Option<u64>> {
    if !path.exists() {
        return Ok(None);
    }
    let reader = BufReader::new(File::open(path)?);
    for line in reader.lines() {
        let line = line?;
        let mut fields = line.split('\t');
        if fields.next() == Some(key) {
            return Ok(fields.next().and_then(|value| value.parse::<u64>().ok()));
        }
    }
    Ok(None)
}

fn gfa_topology_stats(path: &Path) -> io::Result<GfaTopologyStats> {
    let reader = BufReader::new(File::open(path)?);
    let mut segments: HashMap<String, u64> = HashMap::new();
    let mut adjacency: HashMap<String, Vec<String>> = HashMap::new();
    let mut endpoint_degree: HashMap<(String, char), usize> = HashMap::new();
    let mut links = 0usize;
    let mut unique_links = HashSet::new();

    for line in reader.lines() {
        let line = line?;
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.is_empty() {
            continue;
        }
        match fields[0] {
            "S" if fields.len() >= 3 => {
                let len = fields[2].len() as u64;
                segments.insert(fields[1].to_string(), len);
                adjacency.entry(fields[1].to_string()).or_default();
            }
            "L" if fields.len() >= 5 => {
                let from = fields[1].to_string();
                let from_orient = fields[2].chars().next().unwrap_or('+');
                let to = fields[3].to_string();
                let to_orient = fields[4].chars().next().unwrap_or('+');
                links += 1;
                adjacency.entry(from.clone()).or_default().push(to.clone());
                adjacency.entry(to.clone()).or_default().push(from.clone());
                *endpoint_degree
                    .entry((from.clone(), from_orient))
                    .or_insert(0) += 1;
                *endpoint_degree.entry((to.clone(), to_orient)).or_insert(0) += 1;
                unique_links.insert(canonical_gfa_link_key(&from, from_orient, &to, to_orient));
            }
            _ => {}
        }
    }

    let components = graph_component_count(&segments, &adjacency);
    let mut open_ends = 0usize;
    let mut branch_ends = 0usize;
    let mut endpoint_overload = 0usize;
    for name in segments.keys() {
        for orient in ['+', '-'] {
            let degree = endpoint_degree
                .get(&(name.clone(), orient))
                .copied()
                .unwrap_or(0);
            if degree == 0 {
                open_ends += 1;
            } else if degree > 1 {
                branch_ends += 1;
                endpoint_overload += degree - 1;
            }
        }
    }
    let cycle_rank = unique_links
        .len()
        .saturating_sub(segments.len())
        .saturating_add(components);
    Ok(GfaTopologyStats {
        segments: segments.len(),
        links,
        graph_bases: segments.values().sum(),
        components,
        open_ends,
        endpoint_overload,
        branch_ends,
        cycle_rank,
    })
}

fn graph_component_count(
    segments: &HashMap<String, u64>,
    adjacency: &HashMap<String, Vec<String>>,
) -> usize {
    let mut seen = HashSet::new();
    let mut components = 0usize;
    for name in segments.keys() {
        if !seen.insert(name.clone()) {
            continue;
        }
        components += 1;
        let mut stack = vec![name.clone()];
        while let Some(node) = stack.pop() {
            if let Some(neighbors) = adjacency.get(&node) {
                for neighbor in neighbors {
                    if seen.insert(neighbor.clone()) {
                        stack.push(neighbor.clone());
                    }
                }
            }
        }
    }
    components
}

fn canonical_gfa_link_key(from: &str, from_orient: char, to: &str, to_orient: char) -> String {
    let forward = format!("{from}{from_orient}>{to}{to_orient}");
    let reverse = format!(
        "{}{}>{}{}",
        to,
        reverse_orient(to_orient),
        from,
        reverse_orient(from_orient)
    );
    if forward <= reverse {
        forward
    } else {
        reverse
    }
}

fn reverse_orient(orient: char) -> char {
    if orient == '+' {
        '-'
    } else {
        '+'
    }
}

fn selected_graph_in_work_dir(work_dir: &Path) -> Option<PathBuf> {
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
    None
}

fn copy_if_exists(from: PathBuf, to: PathBuf) -> io::Result<()> {
    if from.exists() {
        fs::copy(&from, &to)?;
    }
    Ok(())
}

fn move_if_exists(from: PathBuf, to: PathBuf) -> io::Result<()> {
    if from.exists() {
        if let Some(parent) = to.parent() {
            fs::create_dir_all(parent)?;
        }
        fs::copy(&from, &to)?;
        fs::remove_file(from)?;
    }
    Ok(())
}

fn move_low_depth_bridge_rescue_outputs(
    link_evidence_dir: &Path,
    bridge_rescue_dir: &Path,
) -> io::Result<()> {
    fs::create_dir_all(bridge_rescue_dir)?;
    for file_name in [
        "focused_read_ids.txt",
        "focused_reads.fasta",
        "link_selection.report.txt",
        "selected_links.tsv",
        "pruned_links.tsv",
    ] {
        move_if_exists(
            link_evidence_dir.join(file_name),
            bridge_rescue_dir.join(file_name),
        )?;
    }
    Ok(())
}

fn write_step_status_report(
    step_dir: &Path,
    step: &str,
    status: &str,
    reason: &str,
    input_dir: Option<&Path>,
    output_dir: Option<&Path>,
) -> io::Result<()> {
    fs::create_dir_all(step_dir)?;
    let mut out = File::create(step_dir.join("report.txt"))?;
    writeln!(out, "step_report")?;
    writeln!(out, "step\t{step}")?;
    writeln!(out, "status\t{status}")?;
    writeln!(out, "reason\t{reason}")?;
    if let Some(input_dir) = input_dir {
        writeln!(out, "input_dir\t{}", input_dir.display())?;
    } else {
        writeln!(out, "input_dir\tnone")?;
    }
    if let Some(output_dir) = output_dir {
        writeln!(out, "output_dir\t{}", output_dir.display())?;
    } else {
        writeln!(out, "output_dir\tnone")?;
    }
    Ok(())
}

fn write_low_depth_bridge_rescue_report(
    config: &Config,
    step_dir: &Path,
    status: &str,
    reason: &str,
    read_junction_dir: &Path,
    link_evidence_dir: &Path,
    wrote_selection_outputs: bool,
) -> io::Result<()> {
    fs::create_dir_all(step_dir)?;
    let mut out = File::create(step_dir.join("report.txt"))?;
    writeln!(out, "step_report")?;
    writeln!(out, "step\t04.low_depth_bridge_rescue")?;
    writeln!(out, "status\t{status}")?;
    writeln!(out, "reason\t{reason}")?;
    writeln!(out, "data_mode\t{}", config.data_mode.as_str())?;
    writeln!(
        out,
        "input_read_junction_graph\t{}",
        read_junction_dir.join("graph.gfa").display()
    )?;
    writeln!(
        out,
        "output_skeleton_link_evidence\t{}",
        link_evidence_dir.display()
    )?;
    writeln!(
        out,
        "skeleton_rescue_link_support\t{}",
        config.skeleton_rescue_link_support
    )?;
    if wrote_selection_outputs {
        writeln!(
            out,
            "link_selection_report\t{}",
            step_dir.join("link_selection.report.txt").display()
        )?;
        writeln!(
            out,
            "selected_links\t{}",
            step_dir.join("selected_links.tsv").display()
        )?;
        writeln!(
            out,
            "pruned_links\t{}",
            step_dir.join("pruned_links.tsv").display()
        )?;
        writeln!(
            out,
            "focused_read_ids\t{}",
            step_dir.join("focused_read_ids.txt").display()
        )?;
    } else {
        writeln!(out, "link_selection_report\tskipped")?;
        writeln!(out, "selected_links\tskipped")?;
        writeln!(out, "pruned_links\tskipped")?;
        writeln!(out, "focused_read_ids\tskipped")?;
    }
    Ok(())
}

fn write_one_round_report(
    config: &Config,
    support_dir: &Path,
    graph_dir: &Path,
    read_junction_dir: &Path,
    bridge_rescue_dir: &Path,
    link_evidence_dir: &Path,
    repeat_resolution_dir: &Path,
    linked_graph_dir: &Path,
    summary_dir: &Path,
    started: Instant,
) -> io::Result<()> {
    fs::create_dir_all(summary_dir)?;
    let mut out = File::create(summary_dir.join("report.txt"))?;
    writeln!(out, "orgraft_asm direct-anchor report")?;
    write_workflow_header(&mut out, config, started)?;
    writeln!(out, "graph_workflow\t{}", graph_workflow_label(config))?;
    writeln!(
        out,
        "workflow_steps\t01.anchor_walk_support,02.unitig_graph,03.read_junction_graph,04.low_depth_bridge_rescue,05.skeleton_link_evidence,06.repeat_aware_resolution,07.linked_graph,08.workflow_summary"
    )?;
    writeln!(out, "linking_strategy\tnone")?;
    writeln!(out, "anchor_walk_support\t{}", support_dir.display())?;
    writeln!(out, "unitig_graph\t{}", graph_dir.display())?;
    writeln!(
        out,
        "read_junction_graph\t{} (skipped)",
        read_junction_dir.display()
    )?;
    writeln!(
        out,
        "low_depth_bridge_rescue\t{} (skipped)",
        bridge_rescue_dir.display()
    )?;
    writeln!(
        out,
        "skeleton_link_evidence\t{} (skipped)",
        link_evidence_dir.display()
    )?;
    writeln!(
        out,
        "repeat_aware_resolution\t{} (skipped)",
        repeat_resolution_dir.display()
    )?;
    writeln!(out, "linked_graph\t{}", linked_graph_dir.display())?;
    writeln!(
        out,
        "output_graph\t{}",
        linked_graph_dir.join("graph.gfa").display()
    )?;
    writeln!(out, "output_depth\tskipped")?;
    writeln!(out, "output_links\tskipped")?;
    write_workflow_parameters(&mut out, config)?;
    write_profile_parameter_report(&summary_dir.join("profile_parameters.tsv"), config)?;
    Ok(())
}

fn write_two_round_report(
    config: &Config,
    support_dir: &Path,
    seed_graph_dir: &Path,
    read_junction_dir: &Path,
    bridge_rescue_dir: &Path,
    link_evidence_dir: &Path,
    repeat_resolution_dir: &Path,
    repeat_resolution_status: &str,
    linked_graph_dir: &Path,
    summary_dir: &Path,
    started: Instant,
) -> io::Result<()> {
    fs::create_dir_all(summary_dir)?;
    let mut out = File::create(summary_dir.join("report.txt"))?;
    writeln!(out, "orgraft_asm skeleton-link report")?;
    write_workflow_header(&mut out, config, started)?;
    writeln!(out, "graph_workflow\t{}", graph_workflow_label(config))?;
    writeln!(
        out,
        "workflow_steps\t01.anchor_walk_support,02.unitig_graph,03.read_junction_graph,04.low_depth_bridge_rescue,05.skeleton_link_evidence,06.repeat_aware_resolution,07.linked_graph,08.workflow_summary"
    )?;
    let linking_strategy = if config.mito_stable {
        "repeat_aware"
    } else if is_mito_compact(config) {
        "low_depth_bridge_rescue"
    } else {
        "skeleton_link"
    };
    writeln!(out, "linking_strategy\t{linking_strategy}")?;
    writeln!(out, "anchor_walk_support\t{}", support_dir.display())?;
    writeln!(out, "unitig_graph\t{}", seed_graph_dir.display())?;
    writeln!(out, "read_junction_graph\t{}", read_junction_dir.display())?;
    let bridge_rescue_status = if is_mito_compact(config) {
        "resolved"
    } else {
        "skipped"
    };
    writeln!(
        out,
        "low_depth_bridge_rescue\t{} ({bridge_rescue_status})",
        bridge_rescue_dir.display()
    )?;
    writeln!(
        out,
        "skeleton_link_evidence\t{}",
        link_evidence_dir.display()
    )?;
    writeln!(
        out,
        "repeat_aware_resolution\t{} ({repeat_resolution_status})",
        repeat_resolution_dir.display()
    )?;
    writeln!(out, "linked_graph\t{}", linked_graph_dir.display())?;
    writeln!(
        out,
        "output_graph\t{}",
        linked_graph_dir.join("graph.gfa").display()
    )?;
    writeln!(
        out,
        "output_depth\t{}",
        linked_graph_dir.join("depth.tsv").display()
    )?;
    writeln!(
        out,
        "output_links\t{}",
        linked_graph_dir.join("links.tsv").display()
    )?;
    write_workflow_parameters(&mut out, config)?;
    write_profile_parameter_report(&summary_dir.join("profile_parameters.tsv"), config)?;
    Ok(())
}

fn write_workflow_header<W: Write>(
    out: &mut W,
    config: &Config,
    started: Instant,
) -> io::Result<()> {
    writeln!(
        out,
        "elapsed_seconds\t{:.3}",
        started.elapsed().as_secs_f64()
    )?;
    writeln!(
        out,
        "organelle\t{}",
        config
            .organelle
            .map(|o| o.as_str())
            .unwrap_or("unspecified")
    )?;
    writeln!(out, "data_mode\t{}", config.data_mode.as_str())?;
    Ok(())
}

fn write_workflow_parameters<W: Write>(out: &mut W, config: &Config) -> io::Result<()> {
    writeln!(out, "k\t{}", config.k)?;
    writeln!(out, "s\t{}", config.s)?;
    writeln!(out, "min_anchor_coverage\t{}", config.min_anchor_coverage)?;
    writeln!(out, "min_edge_coverage\t{}", config.min_edge_coverage)?;
    writeln!(out, "min_branch_ratio\t{}", config.min_branch_ratio)?;
    writeln!(out, "max_edges_per_state\t{}", config.max_edges_per_state)?;
    writeln!(out, "dedup_kmer\t{}", config.dedup_kmer)?;
    writeln!(out, "containment_ratio\t{}", config.containment_ratio)?;
    writeln!(out, "min_tip_len\t{}", config.min_tip_len)?;
    writeln!(out, "min_link_support\t{}", config.min_link_support)?;
    writeln!(out, "min_link_ratio\t{}", config.min_link_ratio)?;
    writeln!(out, "read_subset_percent\t{}", config.read_subset_percent)?;
    writeln!(out, "skeleton_end_slop\t{}", config.skeleton_end_slop)?;
    writeln!(
        out,
        "skeleton_min_link_support\t{}",
        config.skeleton_min_link_support
    )?;
    writeln!(
        out,
        "skeleton_min_link_ratio\t{}",
        config.skeleton_min_link_ratio
    )?;
    writeln!(
        out,
        "skeleton_rescue_link_support\t{}",
        config.skeleton_rescue_link_support
    )?;
    writeln!(out, "minimap_min_identity\t{}", config.minimap_min_identity)?;
    writeln!(
        out,
        "minimap_min_align_len\t{}",
        config.minimap_min_align_len
    )?;
    writeln!(out, "paf_max_link_gap\t{}", config.paf_max_link_gap)?;
    writeln!(out, "keep_debug_files\t{}", config.keep_debug_files)?;
    Ok(())
}

fn graph_workflow_label(config: &Config) -> &'static str {
    if config.rounds >= 2 {
        if config.mito_stable {
            "repeat_aware_skeleton_link"
        } else {
            "skeleton_link"
        }
    } else {
        "direct_anchor"
    }
}

fn write_profile_parameter_report(path: &Path, config: &Config) -> io::Result<()> {
    let standard = standard_profile_baseline(config);
    let mut out = File::create(path)?;
    writeln!(
        out,
        "parameter\tstandard_value\tprofile_value\tstatus\tuser_option\tnote"
    )?;
    write_profile_parameter_row(
        &mut out,
        "graph_workflow",
        graph_workflow_label(&standard).to_string(),
        graph_workflow_label(config).to_string(),
        "",
        "step workflow selected by the profile",
    )?;
    write_profile_parameter_row(
        &mut out,
        "min_anchor_coverage",
        standard.min_anchor_coverage.to_string(),
        config.min_anchor_coverage.to_string(),
        "--min-graph-coverage",
        "anchor coverage floor",
    )?;
    write_profile_parameter_row(
        &mut out,
        "min_edge_coverage",
        standard.min_edge_coverage.to_string(),
        config.min_edge_coverage.to_string(),
        "--min-graph-coverage",
        "edge coverage floor",
    )?;
    write_profile_parameter_row(
        &mut out,
        "min_branch_ratio",
        standard.min_branch_ratio.to_string(),
        config.min_branch_ratio.to_string(),
        "",
        "branch pruning ratio",
    )?;
    write_profile_parameter_row(
        &mut out,
        "max_edges_per_state",
        standard.max_edges_per_state.to_string(),
        config.max_edges_per_state.to_string(),
        "",
        "maximum retained outgoing edges per anchor state",
    )?;
    write_profile_parameter_row(
        &mut out,
        "min_tip_len",
        standard.min_tip_len.to_string(),
        config.min_tip_len.to_string(),
        "",
        "short tip pruning length",
    )?;
    write_profile_parameter_row(
        &mut out,
        "min_link_support",
        standard.min_link_support.to_string(),
        config.min_link_support.to_string(),
        "",
        "GFA link support floor",
    )?;
    write_profile_parameter_row(
        &mut out,
        "min_link_ratio",
        standard.min_link_ratio.to_string(),
        config.min_link_ratio.to_string(),
        "--min-link-ratio",
        "optional weak-link ratio filter",
    )?;
    write_profile_parameter_row(
        &mut out,
        "skeleton_min_link_support",
        standard.skeleton_min_link_support.to_string(),
        config.skeleton_min_link_support.to_string(),
        "",
        "skeleton remap link support floor",
    )?;
    write_profile_parameter_row(
        &mut out,
        "skeleton_min_link_ratio",
        standard.skeleton_min_link_ratio.to_string(),
        config.skeleton_min_link_ratio.to_string(),
        "",
        "skeleton remap link ratio floor",
    )?;
    write_profile_parameter_row(
        &mut out,
        "skeleton_rescue_link_support",
        standard.skeleton_rescue_link_support.to_string(),
        config.skeleton_rescue_link_support.to_string(),
        "",
        "skeleton rescue support floor",
    )?;
    write_profile_parameter_row(
        &mut out,
        "repeat_aware_resolution",
        standard.mito_stable.to_string(),
        config.mito_stable.to_string(),
        "",
        "unstable mitogenome additional parameter enables repeat-aware resolution",
    )?;
    write_profile_parameter_row(
        &mut out,
        "read_subsets_requested",
        standard.read_subsets_requested.to_string(),
        config.read_subsets_requested.to_string(),
        "",
        "high profile enables deterministic subset candidates",
    )?;
    write_profile_parameter_row(
        &mut out,
        "read_subsets_auto_recommended",
        standard.read_subsets_auto_recommended.to_string(),
        config.read_subsets_auto_recommended.to_string(),
        "",
        "high profile estimates first-pass subset from input depth when --subsets is not supplied",
    )?;
    write_profile_parameter_row(
        &mut out,
        "read_subset_default_genome_size",
        format_optional_u64(standard.read_subset_default_genome_size),
        format_optional_u64(config.read_subset_default_genome_size),
        "",
        "default genome size for high-profile subset recommendation",
    )?;
    write_profile_parameter_row(
        &mut out,
        "read_subset_input_bases",
        format_optional_u64(standard.read_subset_input_bases),
        format_optional_u64(config.read_subset_input_bases),
        "",
        "input read bases counted before high-profile subset recommendation",
    )?;
    write_profile_parameter_row(
        &mut out,
        "read_subset_target_depth",
        format_optional_f64(standard.read_subset_target_depth),
        format_optional_f64(config.read_subset_target_depth),
        "",
        "standard target depth for high-profile subset recommendation",
    )?;
    write_profile_parameter_row(
        &mut out,
        "read_subsets",
        format_read_subsets(&standard.read_subsets),
        format_read_subsets(&config.read_subsets),
        "--subsets",
        "deterministic read subset percentages",
    )?;
    write_profile_parameter_row(
        &mut out,
        "keep_debug_files",
        standard.keep_debug_files.to_string(),
        config.keep_debug_files.to_string(),
        "--keep-debug-files",
        "retain full debug graphs and materialized subset FASTA files",
    )?;
    Ok(())
}

fn format_optional_u64(value: Option<u64>) -> String {
    value
        .map(|value| value.to_string())
        .unwrap_or_else(|| ".".to_string())
}

fn format_optional_f64(value: Option<f64>) -> String {
    value
        .map(|value| format!("{value:.3}"))
        .unwrap_or_else(|| ".".to_string())
}

fn write_profile_parameter_row<W: Write>(
    out: &mut W,
    parameter: &str,
    standard_value: String,
    profile_value: String,
    user_option: &str,
    note: &str,
) -> io::Result<()> {
    let status = if standard_value == profile_value {
        "same_as_standard"
    } else {
        "profile_modified"
    };
    let user_option = if user_option.is_empty() {
        "."
    } else {
        user_option
    };
    writeln!(
        out,
        "{parameter}\t{standard_value}\t{profile_value}\t{status}\t{user_option}\t{note}"
    )
}

fn standard_profile_baseline(config: &Config) -> Config {
    let mut standard = Config::base();
    standard.organelle = config.organelle;
    standard.data_mode = DataMode::Standard;
    standard.threads = config.threads;
    standard.reads = config.reads.clone();
    standard.out_dir = config.out_dir.clone();
    let overrides = OverrideFlags::default();
    apply_profiles(&mut standard, &overrides);
    standard
}

fn format_read_subsets(subsets: &[ReadSubset]) -> String {
    subsets
        .iter()
        .map(ToString::to_string)
        .collect::<Vec<_>>()
        .join(",")
}

fn full_graph_config(config: &Config) -> Config {
    let mut full = config.clone();
    full.min_branch_ratio = 0.0;
    full.max_edges_per_state = 0;
    full.dedup_kmer = 0;
    full.containment_ratio = 0.0;
    full.min_unitig_len = 0;
    full.min_tip_len = 0;
    full
}

fn validate_config(config: &Config) -> io::Result<()> {
    if config.reads.is_empty() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "at least one --reads/-i/--pacbio-hifi input is required",
        ));
    }
    if config.rounds == 0 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "--rounds must be >= 1",
        ));
    }
    if config.k == 0 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "--k must be > 0",
        ));
    }
    if config.s == 0 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "--s must be > 0",
        ));
    }
    if config.s > config.k {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "--s must be <= --k in syncmer mode",
        ));
    }
    let t = config.syncmer_pos.unwrap_or((config.k - config.s) / 2);
    if t > config.k - config.s {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "--syncmer-pos must be <= k - s",
        ));
    }
    if !(0.0..=1.0).contains(&config.min_branch_ratio) {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "--min-branch-ratio must be between 0 and 1",
        ));
    }
    if !(0.0..=1.0).contains(&config.containment_ratio) {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "--containment-ratio must be between 0 and 1",
        ));
    }
    if !(0.0..=1.0).contains(&config.min_link_ratio) {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "--min-link-ratio must be between 0 and 1",
        ));
    }
    if !(0.0..1.0).contains(&config.hifi_error_rate) {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "--hifi-error-rate must be >= 0 and < 1",
        ));
    }
    if !(0.0..=1.0).contains(&config.minimap_min_identity) {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "--minimap-min-identity must be between 0 and 1",
        ));
    }
    if config.skeleton_only && config.skeleton_gfa.is_none() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "--skeleton-only requires --skeleton-gfa",
        ));
    }
    if !(0.0..=1.0).contains(&config.skeleton_min_link_ratio) {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "--skeleton-min-link-ratio must be between 0 and 1",
        ));
    }
    if config.read_subsets.is_empty()
        || config
            .read_subsets
            .iter()
            .any(|subset| subset.basis_points == 0 || subset.basis_points > READ_SUBSET_SCALE)
    {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "--read-subsets values must be between 1 and 100",
        ));
    }
    if config.read_subset_percent.basis_points == 0
        || config.read_subset_percent.basis_points > READ_SUBSET_SCALE
    {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "internal read subset percent must be between 1 and 100",
        ));
    }
    Ok(())
}

impl Config {
    fn base() -> Self {
        Config {
            reads: Vec::new(),
            out_dir: PathBuf::from("result_graph"),
            organelle: None,
            data_mode: DataMode::Standard,
            numt_interference: NumtInterference::Low,
            k: 501,
            s: 31,
            syncmer_pos: None,
            min_anchor_spacing: 0,
            min_anchor_coverage: 2,
            min_edge_coverage: 2,
            min_branch_ratio: 0.0,
            max_edges_per_state: 0,
            dedup_kmer: 0,
            containment_ratio: 0.0,
            min_unitig_len: 0,
            min_tip_len: 0,
            min_link_support: 0,
            min_link_ratio: 0.0,
            read_junction_links: false,
            bidirectional_links: false,
            junction_rescue_support: 0,
            min_read_len: 0,
            max_reads: None,
            genome_size: None,
            asm_coverage: None,
            min_overlap: None,
            iterations: None,
            hifi_error_rate: 0.003,
            minimap_min_identity: 0.95,
            minimap_min_align_len: 500,
            paf_max_link_gap: 2_000,
            skeleton_gfa: None,
            skeleton_only: false,
            skeleton_end_slop: 1_500,
            skeleton_min_link_support: 10,
            skeleton_min_link_ratio: 0.20,
            skeleton_rescue_gfa: None,
            skeleton_rescue_link_support: 20,
            rounds: 0,
            read_subsets: vec![ReadSubset::full()],
            read_subsets_requested: false,
            read_subsets_auto_recommended: false,
            read_subset_default_genome_size: None,
            read_subset_input_bases: None,
            read_subset_target_depth: None,
            read_subset_percent: ReadSubset::full(),
            keep_debug_files: false,
            threads: DEFAULT_THREADS,
            run_minimap2: false,
            mito_stable: false,
            mito_stable_selection_mode: MitoStableSelectionMode::Auto,
            mito_stable_selected_nodes: HashSet::new(),
            mito_stable_allow_three_way: HashSet::new(),
            mito_stable_forbid_three_way: HashSet::new(),
            mito_stable_force_links: Vec::new(),
            mito_stable_drop_links: Vec::new(),
        }
    }

    fn from_draft_request(request: DraftAssemblyRequest) -> Result<Self, String> {
        if request.threads == 0 {
            return Err("threads must be greater than 0".to_string());
        }

        let organelle = request.organelle;
        let mut config = Config::base();
        config.reads = request.reads;
        config.out_dir = request.out_dir;
        config.threads = request.threads;
        config.keep_debug_files = request.keep_debug_files;
        config.organelle = Some(match organelle {
            DraftOrganelle::Mito => OrganelleProfile::Mito,
            DraftOrganelle::Plastid => OrganelleProfile::Plastid,
        });
        if request.repeat_aware_resolution {
            if organelle != DraftOrganelle::Mito {
                return Err("repeat-aware resolution is only available for mito".to_string());
            }
            config.mito_stable = true;
        }

        config.data_mode = match request.data_mode {
            DraftDataMode::Low => DataMode::Low,
            DraftDataMode::Standard => DataMode::Standard,
        };
        if request.auto_read_subset {
            enable_auto_read_subset(&mut config);
        }

        let mut overrides = OverrideFlags::default();
        if request.auto_read_subset {
            overrides.read_subsets = true;
        }
        if let Some(min_graph_coverage) = request.min_graph_coverage {
            config.min_anchor_coverage = min_graph_coverage;
            config.min_edge_coverage = min_graph_coverage;
            overrides.min_anchor_coverage = true;
            overrides.min_edge_coverage = true;
        }
        if let Some(min_branch_ratio) = request.min_branch_ratio {
            config.min_branch_ratio = min_branch_ratio;
            overrides.min_branch_ratio = true;
        }
        if let Some(min_tip_len) = request.min_tip_len {
            config.min_tip_len = min_tip_len;
            overrides.min_tip_len = true;
        }
        if let Some(min_link_support) = request.min_link_support {
            config.min_link_support = min_link_support;
            overrides.min_link_support = true;
        }
        if let Some(min_link_ratio) = request.min_link_ratio {
            config.min_link_ratio = min_link_ratio;
            overrides.min_link_ratio = true;
        }
        if let Some(read_subsets) = request.read_subsets {
            config.read_subsets = read_subsets_from_basis_points(read_subsets)?;
            config.read_subsets_requested = true;
            config.read_subsets_auto_recommended = false;
            overrides.read_subsets = true;
        }
        apply_mito_stable_selection_mode(&mut config)?;
        apply_profiles(&mut config, &overrides);
        Ok(config)
    }

    // Legacy standalone parser kept for regression tests.
    // OrgRAFT's command path enters through from_draft_request above.
    #[cfg(test)]
    fn from_args(args: Vec<String>) -> Result<Self, String> {
        let mut config = Config::base();

        let mut overrides = OverrideFlags::default();
        let mut i = 1;
        while i < args.len() {
            match args[i].as_str() {
                "-i" | "--reads" | "--pacbio-hifi" => {
                    i += 1;
                    config
                        .reads
                        .push(PathBuf::from(take_arg(&args, i, "input path")?));
                }
                "-o" | "--out-dir" => {
                    i += 1;
                    config.out_dir = PathBuf::from(take_arg(&args, i, "output directory")?);
                }
                "--organelle" => {
                    i += 1;
                    config.organelle = Some(match take_arg(&args, i, "organelle")?.as_str() {
                        "plastid" | "chloroplast" | "cp" => OrganelleProfile::Plastid,
                        "mito" | "mitochondria" | "mitochondrion" | "mt" => OrganelleProfile::Mito,
                        other => return Err(format!("unknown --organelle value: {other}")),
                    });
                }
                "--data-mode" => {
                    i += 1;
                    config.data_mode = match take_arg(&args, i, "data mode")?.as_str() {
                        "standard" => DataMode::Standard,
                        "low" => DataMode::Low,
                        other => return Err(format!("unknown --data-mode value: {other}")),
                    };
                }
                "--numt-interference" => {
                    i += 1;
                    config.numt_interference = match take_arg(&args, i, "NUMT interference")?
                        .as_str()
                    {
                        "low" | "small" | "minor" => NumtInterference::Low,
                        "high" | "large" | "strong" => NumtInterference::High,
                        other => return Err(format!("unknown --numt-interference value: {other}")),
                    };
                    overrides.numt_interference = true;
                }
                "--k" | "-k" => {
                    i += 1;
                    config.k = parse_usize(take_arg(&args, i, "k")?, "--k")?;
                    overrides.k = true;
                }
                "--s" => {
                    i += 1;
                    config.s = parse_usize(take_arg(&args, i, "s")?, "--s")?;
                    overrides.s = true;
                }
                "--syncmer-pos" => {
                    i += 1;
                    config.syncmer_pos = Some(parse_usize(
                        take_arg(&args, i, "syncmer position")?,
                        "--syncmer-pos",
                    )?);
                }
                "--min-anchor-spacing" => {
                    i += 1;
                    config.min_anchor_spacing = parse_usize(
                        take_arg(&args, i, "minimum anchor spacing")?,
                        "--min-anchor-spacing",
                    )?;
                }
                "--min-anchor-coverage" => {
                    i += 1;
                    config.min_anchor_coverage = parse_u32(
                        take_arg(&args, i, "min anchor coverage")?,
                        "--min-anchor-coverage",
                    )?;
                    overrides.min_anchor_coverage = true;
                }
                "--min-edge-coverage" => {
                    i += 1;
                    config.min_edge_coverage = parse_u32(
                        take_arg(&args, i, "min edge coverage")?,
                        "--min-edge-coverage",
                    )?;
                    overrides.min_edge_coverage = true;
                }
                "--min-branch-ratio" => {
                    i += 1;
                    config.min_branch_ratio = parse_f64(
                        take_arg(&args, i, "min branch ratio")?,
                        "--min-branch-ratio",
                    )?;
                    overrides.min_branch_ratio = true;
                }
                "--max-edges-per-state" => {
                    i += 1;
                    config.max_edges_per_state = parse_usize(
                        take_arg(&args, i, "max edges per state")?,
                        "--max-edges-per-state",
                    )?;
                    overrides.max_edges_per_state = true;
                }
                "--dedup-kmer" => {
                    i += 1;
                    config.dedup_kmer =
                        parse_usize(take_arg(&args, i, "dedup k-mer")?, "--dedup-kmer")?;
                    overrides.dedup_kmer = true;
                }
                "--containment-ratio" => {
                    i += 1;
                    config.containment_ratio = parse_f64(
                        take_arg(&args, i, "containment ratio")?,
                        "--containment-ratio",
                    )?;
                    overrides.containment_ratio = true;
                }
                "--min-unitig-len" => {
                    i += 1;
                    config.min_unitig_len =
                        parse_usize(take_arg(&args, i, "min unitig length")?, "--min-unitig-len")?;
                }
                "--min-tip-len" => {
                    i += 1;
                    config.min_tip_len =
                        parse_usize(take_arg(&args, i, "min tip length")?, "--min-tip-len")?;
                    overrides.min_tip_len = true;
                }
                "--min-link-support" => {
                    i += 1;
                    config.min_link_support = parse_u32(
                        take_arg(&args, i, "min link support")?,
                        "--min-link-support",
                    )?;
                    overrides.min_link_support = true;
                }
                "--min-link-ratio" => {
                    i += 1;
                    config.min_link_ratio =
                        parse_f64(take_arg(&args, i, "min link ratio")?, "--min-link-ratio")?;
                    overrides.min_link_ratio = true;
                }
                value
                    if value.starts_with("--read-subsets=") || value.starts_with("--subsets=") =>
                {
                    let (flag, subset_spec) = value
                        .split_once('=')
                        .expect("matched subset option with assignment");
                    config.read_subsets = parse_percent_list(subset_spec, flag)?;
                    config.read_subsets_requested = true;
                }
                value if value.starts_with("--mito-stable-mode=") => {
                    let (flag, mode) = value
                        .split_once('=')
                        .expect("matched mito-stable-mode option with assignment");
                    config.mito_stable_selection_mode =
                        parse_mito_stable_selection_mode(mode, flag)?;
                }
                value
                    if value.starts_with("--selected-nodes=")
                        || value.starts_with("--mito-stable-selected=") =>
                {
                    let (flag, names) = value
                        .split_once('=')
                        .expect("matched selected-nodes option with assignment");
                    config.mito_stable_selected_nodes = parse_name_set(names, flag)?;
                }
                value
                    if value.starts_with("--mito-stable-allow-three-way=")
                        || value.starts_with("--allow-three-way=") =>
                {
                    let (flag, names) = value
                        .split_once('=')
                        .expect("matched allow-three-way option with assignment");
                    let selected = parse_name_set(names, flag)?;
                    config.mito_stable_selection_mode = MitoStableSelectionMode::AllowSelected;
                    config.mito_stable_selected_nodes = selected.clone();
                    config.mito_stable_allow_three_way = selected;
                    config.mito_stable_forbid_three_way.clear();
                }
                value
                    if value.starts_with("--mito-stable-forbid-three-way=")
                        || value.starts_with("--forbid-three-way=") =>
                {
                    let (flag, names) = value
                        .split_once('=')
                        .expect("matched forbid-three-way option with assignment");
                    let selected = parse_name_set(names, flag)?;
                    config.mito_stable_selection_mode = MitoStableSelectionMode::ForbidSelected;
                    config.mito_stable_selected_nodes = selected.clone();
                    config.mito_stable_forbid_three_way = selected;
                    config.mito_stable_allow_three_way.clear();
                }
                value
                    if value.starts_with("--mito-stable-force-link=")
                        || value.starts_with("--force-link=") =>
                {
                    let (flag, links) = value
                        .split_once('=')
                        .expect("matched force-link option with assignment");
                    config.mito_stable_force_links = parse_skeleton_link_key_list(links, flag)?;
                }
                value
                    if value.starts_with("--mito-stable-drop-link=")
                        || value.starts_with("--drop-link=") =>
                {
                    let (flag, links) = value
                        .split_once('=')
                        .expect("matched drop-link option with assignment");
                    config.mito_stable_drop_links = parse_skeleton_link_key_list(links, flag)?;
                }
                "--read-subsets" | "--subsets" => {
                    let flag = args[i].clone();
                    i += 1;
                    config.read_subsets =
                        parse_percent_list(&take_arg(&args, i, "read subsets")?, &flag)?;
                    config.read_subsets_requested = true;
                }
                "--mito-stable-mode" => {
                    let flag = args[i].clone();
                    i += 1;
                    let mode = take_arg(&args, i, "mito stable mode")?;
                    config.mito_stable_selection_mode =
                        parse_mito_stable_selection_mode(&mode, &flag)?;
                }
                "--selected-nodes" | "--mito-stable-selected" => {
                    let flag = args[i].clone();
                    i += 1;
                    config.mito_stable_selected_nodes =
                        parse_name_set(&take_arg(&args, i, "selected nodes")?, &flag)?;
                }
                "--mito-stable-allow-three-way" | "--allow-three-way" => {
                    let flag = args[i].clone();
                    i += 1;
                    let selected =
                        parse_name_set(&take_arg(&args, i, "allowed three-way nodes")?, &flag)?;
                    config.mito_stable_selection_mode = MitoStableSelectionMode::AllowSelected;
                    config.mito_stable_selected_nodes = selected.clone();
                    config.mito_stable_allow_three_way = selected;
                    config.mito_stable_forbid_three_way.clear();
                }
                "--mito-stable-forbid-three-way" | "--forbid-three-way" => {
                    let flag = args[i].clone();
                    i += 1;
                    let selected =
                        parse_name_set(&take_arg(&args, i, "forbidden three-way nodes")?, &flag)?;
                    config.mito_stable_selection_mode = MitoStableSelectionMode::ForbidSelected;
                    config.mito_stable_selected_nodes = selected.clone();
                    config.mito_stable_forbid_three_way = selected;
                    config.mito_stable_allow_three_way.clear();
                }
                "--mito-stable-force-link" | "--force-link" => {
                    let flag = args[i].clone();
                    i += 1;
                    config.mito_stable_force_links =
                        parse_skeleton_link_key_list(&take_arg(&args, i, "forced links")?, &flag)?;
                }
                "--mito-stable-drop-link" | "--drop-link" => {
                    let flag = args[i].clone();
                    i += 1;
                    config.mito_stable_drop_links =
                        parse_skeleton_link_key_list(&take_arg(&args, i, "dropped links")?, &flag)?;
                }
                "--read-junction-links" => {
                    config.read_junction_links = true;
                    overrides.read_junction_links = true;
                }
                "--bidirectional-links" => {
                    config.bidirectional_links = true;
                    overrides.bidirectional_links = true;
                }
                "--junction-rescue-support" => {
                    i += 1;
                    config.junction_rescue_support = parse_u32(
                        take_arg(&args, i, "junction rescue support")?,
                        "--junction-rescue-support",
                    )?;
                }
                "--min-read-len" => {
                    i += 1;
                    config.min_read_len =
                        parse_usize(take_arg(&args, i, "min read length")?, "--min-read-len")?;
                }
                "--max-reads" => {
                    i += 1;
                    config.max_reads = Some(parse_usize(
                        take_arg(&args, i, "max reads")?,
                        "--max-reads",
                    )?);
                }
                "--genome-size" => {
                    i += 1;
                    config.genome_size = Some(parse_size(take_arg(&args, i, "genome size")?)?);
                }
                "--asm-coverage" => {
                    i += 1;
                    config.asm_coverage = Some(parse_f64(
                        take_arg(&args, i, "assembly coverage")?,
                        "--asm-coverage",
                    )?);
                }
                "--min-overlap" => {
                    i += 1;
                    config.min_overlap = Some(parse_usize(
                        take_arg(&args, i, "min overlap")?,
                        "--min-overlap",
                    )?);
                }
                "--iterations" => {
                    i += 1;
                    config.iterations = Some(parse_usize(
                        take_arg(&args, i, "iterations")?,
                        "--iterations",
                    )?);
                }
                "--hifi-error-rate" => {
                    i += 1;
                    config.hifi_error_rate =
                        parse_f64(take_arg(&args, i, "HiFi error rate")?, "--hifi-error-rate")?;
                }
                "--hifi-accuracy" => {
                    i += 1;
                    let accuracy =
                        parse_f64(take_arg(&args, i, "HiFi accuracy")?, "--hifi-accuracy")?;
                    config.hifi_error_rate = 1.0 - accuracy;
                }
                "--minimap-min-identity" => {
                    i += 1;
                    config.minimap_min_identity = parse_f64(
                        take_arg(&args, i, "minimap minimum identity")?,
                        "--minimap-min-identity",
                    )?;
                }
                "--minimap-min-align-len" => {
                    i += 1;
                    config.minimap_min_align_len = parse_usize(
                        take_arg(&args, i, "minimap minimum alignment length")?,
                        "--minimap-min-align-len",
                    )?;
                }
                "--paf-max-link-gap" => {
                    i += 1;
                    config.paf_max_link_gap = take_arg(&args, i, "PAF max link gap")?
                        .parse::<isize>()
                        .map_err(|_| "invalid --paf-max-link-gap".to_string())?;
                }
                "--skeleton-gfa" => {
                    i += 1;
                    config.skeleton_gfa = Some(PathBuf::from(take_arg(&args, i, "skeleton GFA")?));
                }
                "--skeleton-only" => {
                    config.skeleton_only = true;
                }
                "--skeleton-end-slop" => {
                    i += 1;
                    config.skeleton_end_slop = parse_usize(
                        take_arg(&args, i, "skeleton end slop")?,
                        "--skeleton-end-slop",
                    )?;
                }
                "--skeleton-min-link-support" => {
                    i += 1;
                    config.skeleton_min_link_support = parse_u32(
                        take_arg(&args, i, "skeleton minimum link support")?,
                        "--skeleton-min-link-support",
                    )?;
                    overrides.skeleton_min_link_support = true;
                }
                "--skeleton-min-link-ratio" => {
                    i += 1;
                    config.skeleton_min_link_ratio = parse_f64(
                        take_arg(&args, i, "skeleton minimum link ratio")?,
                        "--skeleton-min-link-ratio",
                    )?;
                    overrides.skeleton_min_link_ratio = true;
                }
                "--skeleton-rescue-gfa" => {
                    i += 1;
                    config.skeleton_rescue_gfa =
                        Some(PathBuf::from(take_arg(&args, i, "skeleton rescue GFA")?));
                }
                "--skeleton-rescue-link-support" => {
                    i += 1;
                    config.skeleton_rescue_link_support = parse_u32(
                        take_arg(&args, i, "skeleton rescue link support")?,
                        "--skeleton-rescue-link-support",
                    )?;
                    overrides.skeleton_rescue_link_support = true;
                }
                "--rounds" => {
                    i += 1;
                    config.rounds = parse_usize(take_arg(&args, i, "rounds")?, "--rounds")?;
                    overrides.rounds = true;
                }
                "-t" | "--threads" => {
                    i += 1;
                    config.threads = parse_usize(take_arg(&args, i, "threads")?, "--threads")?;
                }
                "--minimap2" => {
                    config.run_minimap2 = true;
                }
                "--keep-debug-files" => {
                    config.keep_debug_files = true;
                }
                "--mito-stable" => {
                    config.mito_stable = true;
                    config.organelle = Some(OrganelleProfile::Mito);
                    config.data_mode = DataMode::Standard;
                }
                "--help-advanced" => {
                    print_advanced_usage();
                    std::process::exit(0);
                }
                "-h" | "--help" => {
                    print_usage();
                    std::process::exit(0);
                }
                value if value.starts_with('-') => return Err(format!("unknown option: {value}")),
                value => {
                    config.reads.push(PathBuf::from(value));
                }
            }
            i += 1;
        }

        apply_mito_stable_selection_mode(&mut config)?;
        apply_profiles(&mut config, &overrides);
        Ok(config)
    }
}

fn apply_mito_stable_selection_mode(config: &mut Config) -> Result<(), String> {
    match config.mito_stable_selection_mode {
        MitoStableSelectionMode::Auto => {
            if !config.mito_stable_selected_nodes.is_empty() {
                return Err(
                    "selected nodes require mito-stable selection mode forbid-selected or allow-selected"
                        .to_string(),
                );
            }
            config.mito_stable_allow_three_way.clear();
            config.mito_stable_forbid_three_way.clear();
        }
        #[cfg(test)]
        MitoStableSelectionMode::ForbidSelected => {
            if config.mito_stable_selected_nodes.is_empty() {
                return Err(
                    "mito-stable selection mode forbid-selected requires selected nodes"
                        .to_string(),
                );
            }
            config.mito_stable_forbid_three_way = config.mito_stable_selected_nodes.clone();
            config.mito_stable_allow_three_way.clear();
        }
        MitoStableSelectionMode::AllowSelected => {
            if config.mito_stable_selected_nodes.is_empty() {
                return Err(
                    "mito-stable selection mode allow-selected requires selected nodes".to_string(),
                );
            }
            config.mito_stable_allow_three_way.clear();
            config.mito_stable_forbid_three_way.clear();
        }
    }
    Ok(())
}

fn enable_auto_read_subset(config: &mut Config) {
    config.read_subsets = vec![ReadSubset { basis_points: 500 }];
    config.read_subsets_requested = true;
    config.read_subsets_auto_recommended = true;
}

fn apply_profiles(config: &mut Config, overrides: &OverrideFlags) {
    if !overrides.rounds && config.rounds == 0 {
        config.rounds = match config.organelle {
            Some(OrganelleProfile::Mito) => 2,
            _ => 1,
        };
    }
    if config.organelle == Some(OrganelleProfile::Mito) && !overrides.numt_interference {
        config.numt_interference = NumtInterference::High;
    }

    if config.data_mode == DataMode::Low {
        apply_compact_profile(config, overrides);
        apply_mito_stable_profile(config, overrides);
        return;
    }

    match (config.organelle, config.numt_interference) {
        (Some(OrganelleProfile::Plastid), _) => {
            apply_common_organelle_defaults(config, overrides);
            if !overrides.min_anchor_coverage {
                config.min_anchor_coverage = 18;
            }
            if !overrides.min_edge_coverage {
                config.min_edge_coverage = 18;
            }
            if !overrides.min_branch_ratio && config.min_branch_ratio == 0.0 {
                config.min_branch_ratio = 0.30;
            }
            if !overrides.max_edges_per_state && config.max_edges_per_state == 0 {
                config.max_edges_per_state = 3;
            }
            if !overrides.min_tip_len {
                config.min_tip_len = 3000;
            }
            if !overrides.min_link_support {
                config.min_link_support = 20;
            }
            if !overrides.read_junction_links {
                config.read_junction_links = true;
            }
            if !overrides.bidirectional_links {
                config.bidirectional_links = true;
            }
        }
        (Some(OrganelleProfile::Mito), NumtInterference::Low) => {
            apply_common_organelle_defaults(config, overrides);
            if !overrides.min_anchor_coverage {
                config.min_anchor_coverage = 18;
            }
            if !overrides.min_edge_coverage {
                config.min_edge_coverage = 18;
            }
            if !overrides.min_branch_ratio && config.min_branch_ratio == 0.0 {
                config.min_branch_ratio = 0.25;
            }
            if !overrides.max_edges_per_state && config.max_edges_per_state == 0 {
                config.max_edges_per_state = 4;
            }
            if !overrides.min_tip_len {
                config.min_tip_len = 3000;
            }
            if !overrides.min_link_support {
                config.min_link_support = 20;
            }
            if !overrides.bidirectional_links {
                config.bidirectional_links = true;
            }
        }
        (Some(OrganelleProfile::Mito), NumtInterference::High) | (None, NumtInterference::High) => {
            apply_common_organelle_defaults(config, overrides);
            if !overrides.min_anchor_coverage {
                config.min_anchor_coverage = 18;
            }
            if !overrides.min_edge_coverage {
                config.min_edge_coverage = 18;
            }
            if !overrides.min_branch_ratio && config.min_branch_ratio == 0.0 {
                config.min_branch_ratio = 0.30;
            }
            if !overrides.max_edges_per_state && config.max_edges_per_state == 0 {
                config.max_edges_per_state = 3;
            }
            if !overrides.min_tip_len {
                config.min_tip_len = 3000;
            }
            if !overrides.min_link_support {
                config.min_link_support = 20;
            }
            if !overrides.bidirectional_links {
                config.bidirectional_links = true;
            }
        }
        (None, NumtInterference::Low) => {
            if config.rounds == 0 {
                config.rounds = 1;
            }
        }
    }
    apply_mito_stable_profile(config, overrides);
}

fn apply_mito_stable_profile(config: &mut Config, overrides: &OverrideFlags) {
    if !config.mito_stable {
        return;
    }
    if config.organelle != Some(OrganelleProfile::Mito) {
        return;
    }
    if !overrides.rounds {
        config.rounds = 2;
    }
    if !overrides.bidirectional_links {
        config.bidirectional_links = true;
    }
}

fn apply_compact_profile(config: &mut Config, overrides: &OverrideFlags) {
    apply_common_organelle_defaults(config, overrides);
    match config.organelle {
        Some(OrganelleProfile::Mito) => apply_compact_mito_profile(config, overrides),
        Some(OrganelleProfile::Plastid) | None => apply_compact_plastid_profile(config, overrides),
    }
}

fn apply_compact_plastid_profile(config: &mut Config, overrides: &OverrideFlags) {
    if !overrides.min_anchor_coverage {
        config.min_anchor_coverage = 12;
    }
    if !overrides.min_edge_coverage {
        config.min_edge_coverage = 12;
    }
    if !overrides.min_branch_ratio && config.min_branch_ratio == 0.0 {
        config.min_branch_ratio = 0.30;
    }
    if !overrides.max_edges_per_state && config.max_edges_per_state == 0 {
        config.max_edges_per_state = 3;
    }
    if !overrides.min_tip_len {
        config.min_tip_len = 3000;
    }
    if !overrides.min_link_support {
        config.min_link_support = 20;
    }
    if !overrides.min_link_ratio {
        config.min_link_ratio = 0.0;
    }
    if !overrides.read_junction_links {
        config.read_junction_links = true;
    }
    if !overrides.bidirectional_links {
        config.bidirectional_links = true;
    }
}

fn apply_compact_mito_profile(config: &mut Config, overrides: &OverrideFlags) {
    if !overrides.min_anchor_coverage {
        config.min_anchor_coverage = 12;
    }
    if !overrides.min_edge_coverage {
        config.min_edge_coverage = 12;
    }
    if !overrides.min_branch_ratio && config.min_branch_ratio == 0.0 {
        config.min_branch_ratio = 0.30;
    }
    if !overrides.max_edges_per_state && config.max_edges_per_state == 0 {
        config.max_edges_per_state = 3;
    }
    if !overrides.min_tip_len {
        config.min_tip_len = 3000;
    }
    if !overrides.min_link_support {
        config.min_link_support = 20;
    }
    if !overrides.min_link_ratio {
        config.min_link_ratio = 0.0;
    }
    if !overrides.read_junction_links {
        config.read_junction_links = true;
    }
    if !overrides.bidirectional_links {
        config.bidirectional_links = true;
    }
    if !overrides.skeleton_min_link_support {
        config.skeleton_min_link_support = 10;
    }
    if !overrides.skeleton_min_link_ratio {
        config.skeleton_min_link_ratio = 0.20;
    }
    if !overrides.skeleton_rescue_link_support {
        config.skeleton_rescue_link_support = 20;
    }
}

fn apply_common_organelle_defaults(config: &mut Config, overrides: &OverrideFlags) {
    if !overrides.k {
        config.k = 251;
    }
    if !overrides.s {
        config.s = 21;
    }
    if !overrides.dedup_kmer {
        config.dedup_kmer = 17;
    }
    if !overrides.containment_ratio {
        config.containment_ratio = 0.70;
    }
}

#[cfg(test)]
fn take_arg(args: &[String], index: usize, name: &str) -> Result<String, String> {
    args.get(index)
        .cloned()
        .ok_or_else(|| format!("missing value for {name}"))
}

#[cfg(test)]
fn parse_usize(value: String, flag: &str) -> Result<usize, String> {
    value
        .parse::<usize>()
        .map_err(|_| format!("invalid {flag}: {value}"))
}

#[cfg(test)]
fn parse_u32(value: String, flag: &str) -> Result<u32, String> {
    value
        .parse::<u32>()
        .map_err(|_| format!("invalid {flag}: {value}"))
}

#[cfg(test)]
fn parse_f64(value: String, flag: &str) -> Result<f64, String> {
    value
        .parse::<f64>()
        .map_err(|_| format!("invalid {flag}: {value}"))
}

#[cfg(test)]
fn parse_percent_list(value: &str, flag: &str) -> Result<Vec<ReadSubset>, String> {
    let mut percents = Vec::new();
    for raw_part in value.split(',') {
        let part = raw_part.trim();
        if part.is_empty() {
            return Err(format!("invalid {flag} value: {value}"));
        }
        percents.push(parse_read_subset_percent(part, flag)?);
    }
    percents.sort_unstable();
    percents.dedup();
    if percents.is_empty() {
        return Err(format!("{flag} requires at least one value"));
    }
    Ok(percents)
}

#[cfg(test)]
fn parse_name_set(value: &str, flag: &str) -> Result<HashSet<String>, String> {
    let mut names = HashSet::new();
    for raw_part in value.split(',') {
        let part = raw_part.trim();
        if part.is_empty() {
            return Err(format!("invalid {flag} value: {value}"));
        }
        if part.chars().any(char::is_whitespace) {
            return Err(format!("invalid {flag} node name: {part}"));
        }
        names.insert(part.to_string());
    }
    if names.is_empty() {
        return Err(format!("{flag} requires at least one node name"));
    }
    Ok(names)
}

#[cfg(test)]
fn parse_mito_stable_selection_mode(
    value: &str,
    flag: &str,
) -> Result<MitoStableSelectionMode, String> {
    let Some(mode) = MitoStableSelectionMode::from_str(value.trim()) else {
        return Err(format!(
            "unknown {flag} value: {value}; expected auto, forbid-selected, or allow-selected"
        ));
    };
    Ok(mode)
}

#[cfg(test)]
fn parse_skeleton_link_key_list(value: &str, flag: &str) -> Result<Vec<SkeletonLinkKey>, String> {
    let mut links = Vec::new();
    for raw_part in value.split(',') {
        let part = raw_part.trim();
        if part.is_empty() {
            return Err(format!("invalid {flag} value: {value}"));
        }
        links.push(parse_skeleton_link_key(part, flag)?);
    }
    if links.is_empty() {
        return Err(format!("{flag} requires at least one link"));
    }
    Ok(links)
}

#[cfg(test)]
fn parse_skeleton_link_key(value: &str, flag: &str) -> Result<SkeletonLinkKey, String> {
    let fields = value.split(':').collect::<Vec<_>>();
    if fields.len() != 4 {
        return Err(format!(
            "invalid {flag} link: {value}; expected from:from_orient:to:to_orient"
        ));
    }
    let from = fields[0].trim();
    let from_orient = parse_link_orient(fields[1].trim(), flag, value)?;
    let to = fields[2].trim();
    let to_orient = parse_link_orient(fields[3].trim(), flag, value)?;
    if from.is_empty() || to.is_empty() {
        return Err(format!("invalid {flag} link: {value}"));
    }
    Ok(SkeletonLinkKey {
        from: from.to_string(),
        from_orient,
        to: to.to_string(),
        to_orient,
    })
}

#[cfg(test)]
fn parse_link_orient(value: &str, flag: &str, link: &str) -> Result<char, String> {
    let mut chars = value.chars();
    let Some(orient) = chars.next() else {
        return Err(format!("invalid {flag} orientation in link: {link}"));
    };
    if chars.next().is_some() || !matches!(orient, '+' | '-') {
        return Err(format!("invalid {flag} orientation in link: {link}"));
    }
    Ok(orient)
}

#[cfg(test)]
fn parse_read_subset_percent(value: &str, flag: &str) -> Result<ReadSubset, String> {
    let (whole_part, decimal_part) = match value.split_once('.') {
        Some((whole, decimal)) => (whole, Some(decimal)),
        None => (value, None),
    };
    if whole_part.is_empty() || whole_part.starts_with('+') || whole_part.starts_with('-') {
        return Err(format!("invalid {flag} value: {value}"));
    }

    let whole = whole_part
        .parse::<u16>()
        .map_err(|_| format!("invalid {flag} value: {value}"))?;
    let decimal = match decimal_part {
        Some(decimal) => parse_percent_decimal(decimal, flag, value)?,
        None => 0,
    };
    let basis_points = whole
        .checked_mul(100)
        .and_then(|scaled| scaled.checked_add(decimal))
        .ok_or_else(|| format!("{flag} value must be > 0 and <= 100: {value}"))?;
    if basis_points == 0 || basis_points > READ_SUBSET_SCALE {
        return Err(format!("{flag} value must be > 0 and <= 100: {value}"));
    }
    Ok(ReadSubset { basis_points })
}

#[cfg(test)]
fn parse_percent_decimal(decimal: &str, flag: &str, value: &str) -> Result<u16, String> {
    if decimal.is_empty() || decimal.len() > 2 || !decimal.bytes().all(|b| b.is_ascii_digit()) {
        return Err(format!(
            "{flag} values support at most two decimal places: {value}"
        ));
    }
    let parsed = decimal
        .parse::<u16>()
        .map_err(|_| format!("invalid {flag} value: {value}"))?;
    Ok(if decimal.len() == 1 {
        parsed * 10
    } else {
        parsed
    })
}

#[cfg(test)]
fn parse_size(value: String) -> Result<u64, String> {
    let trimmed = value.trim();
    if trimmed.is_empty() {
        return Err("empty size".to_string());
    }
    let (number, multiplier) = match trimmed.as_bytes()[trimmed.len() - 1].to_ascii_lowercase() {
        b'k' => (&trimmed[..trimmed.len() - 1], 1_000_u64),
        b'm' => (&trimmed[..trimmed.len() - 1], 1_000_000_u64),
        b'g' => (&trimmed[..trimmed.len() - 1], 1_000_000_000_u64),
        _ => (trimmed, 1_u64),
    };
    let parsed = number
        .parse::<f64>()
        .map_err(|_| format!("invalid size: {value}"))?;
    Ok((parsed * multiplier as f64).round() as u64)
}

#[cfg(test)]
fn print_usage() {
    eprintln!(
        "Usage: simple_draft_asm --pacbio-hifi reads.fastq.gz -o result_graph [options]\n\
         \n\
         Common options:\n\
           --organelle plastid|mito     plastid defaults to 1 round; mito defaults to 2 rounds\n\
           --data-mode standard|low     low mode for low-depth or corrected-read datasets\n\
           --mito-stable                unstable mitogenomes; enable repeat-aware resolution\n\
           --rounds INT                 override profile rounds\n\
           --numt-interference low|high mito profile strictness, mito default: high\n\
           -i, --pacbio-hifi FILE       input reads; may be repeated\n\
           -o DIR                       output directory\n\
           -t INT                       threads\n\
           --help-advanced              show tuning and skeleton-linking options\n\
         \n\
         Profile defaults:\n\
           low:     corrected-read or low-depth mode\n\
           standard: organelle standard profile; plastid stays 1 round, mito stays 2 rounds\n\
           high:    high-data mode: standard profile plus subset candidate search\n\
           --mito-stable: add repeat-aware resolution for unstable mitogenomes\n\
         \n\
         Typical commands:\n\
           simple_draft_asm --organelle plastid -i data/plastid.fastq.gz -o result_plastid\n\
           simple_draft_asm --organelle mito -i data/mito.fastq.gz -o result_mito\n\
           simple_draft_asm --organelle mito --data-mode low -i data/mecat_mito.fasta.gz -o result_mito_low"
    );
}

#[cfg(test)]
fn print_advanced_usage() {
    eprintln!(
        "Advanced options supported by simple_draft_asm:\n\
         \n\
         Anchor selection:\n\
           --k INT\n\
           --s INT\n\
           --syncmer-pos INT\n\
           --min-anchor-spacing INT\n\
         \n\
         First-round graph filtering:\n\
           --min-anchor-coverage INT\n\
           --min-edge-coverage INT\n\
           --min-branch-ratio FLOAT\n\
           --max-edges-per-state INT\n\
           --junction-rescue-support INT\n\
           --dedup-kmer INT\n\
           --containment-ratio FLOAT\n\
           --min-unitig-len INT\n\
           --min-tip-len INT\n\
           --min-link-support INT\n\
           --min-link-ratio FLOAT\n\
           --read-subsets, --subsets LIST\n\
                                      fan out deterministic read subsets in one pass, e.g. 12.5,25,50,100\n\
           --read-junction-links\n\
           --bidirectional-links\n\
         \n\
         HiFi and minimap2 support:\n\
           --hifi-error-rate FLOAT\n\
           --hifi-accuracy FLOAT\n\
           --minimap2\n\
           --mito-stable\n\
           --mito-stable-mode auto|forbid-selected|allow-selected\n\
                                      mito-stable node selection mode; auto is default\n\
           --selected-nodes LIST\n\
                                      nodes used by forbid-selected or allow-selected mode\n\
           --allow-three-way LIST\n\
                                      alias for --mito-stable-mode allow-selected --selected-nodes LIST\n\
           --forbid-three-way LIST\n\
                                      alias for --mito-stable-mode forbid-selected --selected-nodes LIST\n\
           --force-link LINK[,LINK]\n\
                                      mito-stable: force GFA-style links, e.g. edge_8:-:edge_4:-\n\
           --drop-link LINK[,LINK]\n\
                                      mito-stable: remove GFA-style links before final validation\n\
           --minimap-min-identity FLOAT\n\
           --minimap-min-align-len INT\n\
           --paf-max-link-gap INT\n\
         \n\
         Skeleton second round:\n\
           --skeleton-gfa FILE\n\
           --skeleton-only\n\
           --skeleton-end-slop INT\n\
           --skeleton-min-link-support INT\n\
           --skeleton-min-link-ratio FLOAT\n\
           --skeleton-rescue-gfa FILE\n\
           --skeleton-rescue-link-support INT"
    );
}

// Step 01 parser support: read FASTQ/FASTA, optionally through gzip, and stream
// sequence records into the anchor graph builder.
fn read_sequence_file<F>(path: &Path, callback: F) -> io::Result<()>
where
    F: FnMut(&str, &str) -> io::Result<bool>,
{
    let mut reader = open_text(path)?;
    let mut line = String::new();
    if reader.read_line(&mut line)? == 0 {
        return Ok(());
    }
    trim_line(&mut line);

    if line.starts_with('@') {
        read_fastq(reader, line, callback)
    } else if line.starts_with('>') {
        read_fasta(reader, line, callback)
    } else {
        Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!("{} is not FASTA/FASTQ", path.display()),
        ))
    }
}

fn open_text(path: &Path) -> io::Result<TextReader> {
    if path.extension().and_then(|ext| ext.to_str()) == Some("gz") {
        let mut child = Command::new("gzip")
            .arg("-dc")
            .arg(path)
            .stdout(Stdio::piped())
            .spawn()
            .map_err(|err| {
                io::Error::new(
                    err.kind(),
                    format!("failed to spawn gzip for {}: {err}", path.display()),
                )
            })?;
        let stdout = child.stdout.take().ok_or_else(|| {
            io::Error::other(format!(
                "failed to capture gzip stdout for {}",
                path.display()
            ))
        })?;
        Ok(TextReader::Gzip {
            child,
            reader: BufReader::new(stdout),
        })
    } else {
        Ok(TextReader::Plain(BufReader::new(File::open(path)?)))
    }
}

fn read_fastq<F>(mut reader: TextReader, first_header: String, mut callback: F) -> io::Result<()>
where
    F: FnMut(&str, &str) -> io::Result<bool>,
{
    let mut header = first_header;
    let mut seq = String::new();
    let mut plus = String::new();
    let mut qual = String::new();
    loop {
        if !header.starts_with('@') {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("expected FASTQ header, got {header}"),
            ));
        }
        let name = header[1..].split_whitespace().next().unwrap_or("");

        seq.clear();
        plus.clear();
        qual.clear();
        if reader.read_line(&mut seq)? == 0 {
            return Err(io::Error::new(
                io::ErrorKind::UnexpectedEof,
                "truncated FASTQ sequence",
            ));
        }
        if reader.read_line(&mut plus)? == 0 {
            return Err(io::Error::new(
                io::ErrorKind::UnexpectedEof,
                "truncated FASTQ plus line",
            ));
        }
        if reader.read_line(&mut qual)? == 0 {
            return Err(io::Error::new(
                io::ErrorKind::UnexpectedEof,
                "truncated FASTQ quality",
            ));
        }
        trim_line(&mut seq);
        trim_line(&mut plus);
        trim_line(&mut qual);
        if !plus.starts_with('+') {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("expected FASTQ plus line for {name}"),
            ));
        }
        if !callback(name, &seq)? {
            break;
        }

        header.clear();
        if reader.read_line(&mut header)? == 0 {
            break;
        }
        trim_line(&mut header);
    }
    Ok(())
}

fn read_fasta<F>(mut reader: TextReader, first_header: String, mut callback: F) -> io::Result<()>
where
    F: FnMut(&str, &str) -> io::Result<bool>,
{
    let mut header = first_header;
    let mut seq = String::new();
    loop {
        let name = header[1..]
            .split_whitespace()
            .next()
            .unwrap_or("")
            .to_string();
        seq.clear();
        loop {
            let mut line = String::new();
            let bytes = reader.read_line(&mut line)?;
            if bytes == 0 {
                if !callback(&name, &seq)? {
                    return Ok(());
                }
                return Ok(());
            }
            trim_line(&mut line);
            if line.starts_with('>') {
                if !callback(&name, &seq)? {
                    return Ok(());
                }
                header = line;
                break;
            }
            seq.push_str(&line);
        }
    }
}

fn trim_line(line: &mut String) {
    while line.ends_with('\n') || line.ends_with('\r') {
        line.pop();
    }
}

// Step 02a anchor walk construction: select anchors on each read, register
// anchor states, and accumulate edge/junction/read-walk support.
fn process_read(
    config: &Config,
    assembly: &mut Assembly,
    _name: &str,
    seq: &str,
) -> io::Result<()> {
    assembly.reads_seen += 1;
    assembly.bases_seen += seq.len() as u64;
    let seq_bytes = seq.as_bytes();

    let raw_hits = select_syncmer_hits(seq_bytes, config.k, config.s, config.syncmer_pos);
    let raw_hits = thin_anchor_hits(raw_hits, config.min_anchor_spacing);

    assembly.anchors_seen += raw_hits.len() as u64;

    let mut hits: Vec<AnchorHit> = Vec::with_capacity(raw_hits.len());
    let mut last_state: Option<usize> = None;
    for pos in raw_hits {
        let Some((canonical, is_forward)) = canonical_kmer(seq_bytes, pos, config.k) else {
            continue;
        };
        let node_id = match assembly.key_to_node.get(&canonical) {
            Some(id) => *id,
            None => {
                let id = assembly.nodes.len();
                assembly.nodes.push(AnchorNode {
                    seq: canonical.clone(),
                    coverage: 0,
                });
                assembly.key_to_node.insert(canonical, id);
                id
            }
        };
        assembly.nodes[node_id].coverage += 1;
        let state = state_id(node_id, is_forward);
        if Some(state) != last_state {
            hits.push(AnchorHit { state, pos });
            last_state = Some(state);
        }
    }

    if hits.len() < 2 {
        assembly.read_walks.push(ReadWalk { edges: Vec::new() });
        return Ok(());
    }

    let mut read_edges = Vec::with_capacity(hits.len().saturating_sub(1));
    for pair in hits.windows(2) {
        let left = pair[0];
        let right = pair[1];
        if right.pos <= left.pos {
            continue;
        }
        let suffix_start = left.pos.saturating_add(config.k).min(seq_bytes.len());
        let suffix_end = right.pos.saturating_add(config.k).min(seq_bytes.len());
        if suffix_end < suffix_start {
            continue;
        }
        let span = right.pos.saturating_sub(left.pos);
        let edge_key = EdgeKey {
            from: left.state,
            to: right.state,
        };
        let suffix_len = suffix_end.saturating_sub(suffix_start);
        if let Some(entry) = assembly.edges.get_mut(&edge_key) {
            entry.coverage += 1;
            entry.total_span += span as u64;
            entry.min_span = entry.min_span.min(span);
            entry.max_span = entry.max_span.max(span);
            if suffix_len > entry.suffix.len() {
                entry.suffix = ascii_uppercase_string(&seq_bytes[suffix_start..suffix_end]);
            }
        } else {
            assembly.edges.insert(
                edge_key,
                EdgeStats {
                    coverage: 1,
                    total_span: span as u64,
                    min_span: span,
                    max_span: span,
                    suffix: ascii_uppercase_string(&seq_bytes[suffix_start..suffix_end]),
                },
            );
        }
        read_edges.push(edge_key);
    }

    for pair in read_edges.windows(2) {
        *assembly
            .edge_junctions
            .entry((pair[0], pair[1]))
            .or_insert(0) += 1;
    }

    assembly.read_walks.push(ReadWalk { edges: read_edges });
    Ok(())
}

fn select_syncmer_hits(seq: &[u8], k: usize, s: usize, syncmer_pos: Option<usize>) -> Vec<usize> {
    let target_offset = syncmer_pos.unwrap_or((k - s) / 2);
    let mut hits = Vec::new();

    for (segment_start, segment) in valid_segments(seq) {
        if segment.len() < k {
            continue;
        }
        let s_hashes = hash_smers(segment, s);
        let window = k - s + 1;
        let mut deque: VecDeque<usize> = VecDeque::new();

        for i in 0..s_hashes.len() {
            while let Some(&back) = deque.back() {
                if s_hashes[back] > s_hashes[i] || (s_hashes[back] == s_hashes[i] && back > i) {
                    deque.pop_back();
                } else {
                    break;
                }
            }
            deque.push_back(i);

            let k_start = (i + 1).saturating_sub(window);
            while let Some(&front) = deque.front() {
                if front < k_start {
                    deque.pop_front();
                } else {
                    break;
                }
            }

            if i + 1 >= window {
                let min_pos = *deque.front().expect("syncmer deque is not empty");
                if min_pos == k_start + target_offset {
                    hits.push(segment_start + k_start);
                }
            }
        }
    }

    hits
}

fn thin_anchor_hits(hits: Vec<usize>, min_spacing: usize) -> Vec<usize> {
    if min_spacing <= 1 || hits.len() < 2 {
        return hits;
    }

    let mut thinned = Vec::with_capacity(hits.len());
    let mut last_kept = None;
    for pos in hits {
        if last_kept.is_none_or(|last| pos >= last + min_spacing) {
            thinned.push(pos);
            last_kept = Some(pos);
        }
    }
    thinned
}

fn valid_segments(seq: &[u8]) -> Vec<(usize, &[u8])> {
    let mut segments = Vec::new();
    let mut start = None;
    for (i, &base) in seq.iter().enumerate() {
        if base_code(base).is_some() {
            if start.is_none() {
                start = Some(i);
            }
        } else if let Some(s) = start.take() {
            if i > s {
                segments.push((s, &seq[s..i]));
            }
        }
    }
    if let Some(s) = start {
        if seq.len() > s {
            segments.push((s, &seq[s..]));
        }
    }
    segments
}

fn hash_smers(segment: &[u8], s: usize) -> Vec<u64> {
    if s <= 32 {
        hash_smers_2bit(segment, s)
    } else {
        let mut hashes = Vec::with_capacity(segment.len().saturating_sub(s) + 1);
        for i in 0..=segment.len() - s {
            let forward = stable_hash_bytes(&segment[i..i + s]);
            let reverse = stable_hash_revcomp(&segment[i..i + s]);
            hashes.push(mix64(forward.min(reverse)));
        }
        hashes
    }
}

fn hash_smers_2bit(segment: &[u8], s: usize) -> Vec<u64> {
    let mask = if s == 32 {
        u64::MAX
    } else {
        (1_u64 << (2 * s)) - 1
    };
    let mut forward = 0_u64;
    let mut reverse = 0_u64;
    let mut hashes = Vec::with_capacity(segment.len().saturating_sub(s) + 1);

    for (i, &base) in segment.iter().enumerate() {
        let code = base_code(base).expect("valid segment contains only ACGT");
        forward = ((forward << 2) | code as u64) & mask;
        reverse = (reverse >> 2) | ((3 - code as u64) << (2 * (s - 1)));
        if i + 1 >= s {
            hashes.push(mix64(forward.min(reverse)));
        }
    }
    hashes
}

fn canonical_kmer(seq: &[u8], pos: usize, k: usize) -> Option<(String, bool)> {
    if pos + k > seq.len() {
        return None;
    }
    let slice = &seq[pos..pos + k];
    if !slice.iter().all(|&base| base_code(base).is_some()) {
        return None;
    }

    if forward_leq_revcomp(slice) {
        Some((ascii_uppercase_string(slice), true))
    } else {
        Some((revcomp_bytes_string(slice), false))
    }
}

fn forward_leq_revcomp(slice: &[u8]) -> bool {
    for i in 0..slice.len() {
        let left = slice[i].to_ascii_uppercase();
        let right = complement(slice[slice.len() - 1 - i].to_ascii_uppercase());
        if left != right {
            return left < right;
        }
    }
    true
}

fn ascii_uppercase_string(bytes: &[u8]) -> String {
    let mut out = String::with_capacity(bytes.len());
    for &byte in bytes {
        out.push(byte.to_ascii_uppercase() as char);
    }
    out
}

fn revcomp_bytes_string(bytes: &[u8]) -> String {
    let mut out = String::with_capacity(bytes.len());
    for &byte in bytes.iter().rev() {
        out.push(complement(byte.to_ascii_uppercase()) as char);
    }
    out
}

fn state_id(node_id: usize, is_forward: bool) -> usize {
    node_id * 2 + if is_forward { 0 } else { 1 }
}

fn node_id_from_state(state: usize) -> usize {
    state / 2
}

fn state_orientation(state: usize) -> char {
    if state % 2 == 0 {
        '+'
    } else {
        '-'
    }
}

fn oriented_anchor_seq(node: &AnchorNode, state: usize) -> String {
    if state % 2 == 0 {
        node.seq.clone()
    } else {
        revcomp_string(&node.seq)
    }
}

fn base_code(base: u8) -> Option<u8> {
    match base {
        b'A' | b'a' => Some(0),
        b'C' | b'c' => Some(1),
        b'G' | b'g' => Some(2),
        b'T' | b't' => Some(3),
        _ => None,
    }
}

fn complement(base: u8) -> u8 {
    match base {
        b'A' | b'a' => b'T',
        b'C' | b'c' => b'G',
        b'G' | b'g' => b'C',
        b'T' | b't' => b'A',
        _ => b'N',
    }
}

fn revcomp_string(seq: &str) -> String {
    let mut out = String::with_capacity(seq.len());
    for base in seq.bytes().rev() {
        out.push(complement(base) as char);
    }
    out
}

fn stable_hash_bytes(bytes: &[u8]) -> u64 {
    stable_hash_bytes_seed(bytes, 0xcbf2_9ce4_8422_2325, 0x1000_0000_01b3)
}

fn stable_hash_bytes_seed(bytes: &[u8], seed: u64, multiplier: u64) -> u64 {
    let mut hash = seed;
    for &byte in bytes {
        hash ^= byte.to_ascii_uppercase() as u64;
        hash = hash.wrapping_mul(multiplier);
    }
    hash
}

fn stable_hash_revcomp(bytes: &[u8]) -> u64 {
    stable_hash_revcomp_seed(bytes, 0xcbf2_9ce4_8422_2325, 0x1000_0000_01b3)
}

fn stable_hash_revcomp_seed(bytes: &[u8], seed: u64, multiplier: u64) -> u64 {
    let mut hash = seed;
    for &byte in bytes.iter().rev() {
        hash ^= complement(byte) as u64;
        hash = hash.wrapping_mul(multiplier);
    }
    hash
}

fn mix64(mut x: u64) -> u64 {
    x ^= x >> 33;
    x = x.wrapping_mul(0xff51_afd7_ed55_8ccd);
    x ^= x >> 33;
    x = x.wrapping_mul(0xc4ce_b9fe_1a85_ec53);
    x ^ (x >> 33)
}

// Step 02b graph filtering: keep supported anchors and edges, preserve selected
// branches, and prepare the directed graph for unitig compression.
struct FilteredGraph {
    edges: Vec<(EdgeKey, EdgeStats)>,
    outgoing: HashMap<usize, Vec<EdgeKey>>,
    incoming: HashMap<usize, Vec<EdgeKey>>,
}

fn build_filtered_graph(
    config: &Config,
    assembly: &Assembly,
    edge_junction_support: &HashMap<EdgeKey, u32>,
) -> FilteredGraph {
    let mut candidates = Vec::new();
    let mut outgoing_rank: HashMap<usize, Vec<(EdgeKey, u32)>> = HashMap::new();
    let mut incoming_rank: HashMap<usize, Vec<(EdgeKey, u32)>> = HashMap::new();
    let mut max_outgoing: HashMap<usize, u32> = HashMap::new();
    let mut max_incoming: HashMap<usize, u32> = HashMap::new();

    for (key, stats) in &assembly.edges {
        let from_node = node_id_from_state(key.from);
        let to_node = node_id_from_state(key.to);
        let raw_junction_support = edge_junction_support.get(key).copied().unwrap_or(0);
        let rescue = config.junction_rescue_support > 0
            && raw_junction_support >= config.junction_rescue_support;
        if assembly.nodes[from_node].coverage < config.min_anchor_coverage {
            continue;
        }
        if assembly.nodes[to_node].coverage < config.min_anchor_coverage {
            continue;
        }
        if !rescue && stats.coverage < config.min_edge_coverage {
            continue;
        }

        candidates.push((*key, stats.clone()));
        let rank_support = stats.coverage.max(raw_junction_support);
        outgoing_rank
            .entry(key.from)
            .or_default()
            .push((*key, rank_support));
        incoming_rank
            .entry(key.to)
            .or_default()
            .push((*key, rank_support));
        max_outgoing
            .entry(key.from)
            .and_modify(|max_cov| *max_cov = (*max_cov).max(stats.coverage))
            .or_insert(stats.coverage);
        max_incoming
            .entry(key.to)
            .and_modify(|max_cov| *max_cov = (*max_cov).max(stats.coverage))
            .or_insert(stats.coverage);
    }

    let outgoing_keep = rank_keep_set(&mut outgoing_rank, config.max_edges_per_state);
    let incoming_keep = rank_keep_set(&mut incoming_rank, config.max_edges_per_state);
    let mut edges = Vec::new();
    let mut outgoing: HashMap<usize, Vec<EdgeKey>> = HashMap::new();
    let mut incoming: HashMap<usize, Vec<EdgeKey>> = HashMap::new();

    for (key, stats) in candidates {
        let raw_junction_support = edge_junction_support.get(&key).copied().unwrap_or(0);
        let rescue = config.junction_rescue_support > 0
            && raw_junction_support >= config.junction_rescue_support;
        if config.min_branch_ratio > 0.0 {
            let out_best = max_outgoing
                .get(&key.from)
                .copied()
                .unwrap_or(stats.coverage);
            let in_best = max_incoming.get(&key.to).copied().unwrap_or(stats.coverage);
            let out_ratio = stats.coverage as f64 / out_best.max(1) as f64;
            let in_ratio = stats.coverage as f64 / in_best.max(1) as f64;
            if !rescue
                && (out_ratio < config.min_branch_ratio || in_ratio < config.min_branch_ratio)
            {
                continue;
            }
        }
        if !rescue && !outgoing_keep.is_empty() && !outgoing_keep.contains(&key) {
            continue;
        }
        if !rescue && !incoming_keep.is_empty() && !incoming_keep.contains(&key) {
            continue;
        }
        edges.push((key, stats.clone()));
        outgoing.entry(key.from).or_default().push(key);
        incoming.entry(key.to).or_default().push(key);
    }

    edges.sort_by_key(|(key, _)| (key.from, key.to));
    for values in outgoing.values_mut() {
        values.sort_by_key(|key| key.to);
    }
    for values in incoming.values_mut() {
        values.sort_by_key(|key| key.from);
    }

    FilteredGraph {
        edges,
        outgoing,
        incoming,
    }
}

fn edge_junction_support(assembly: &Assembly) -> HashMap<EdgeKey, u32> {
    let mut support = HashMap::new();
    for (&(left, right), &count) in &assembly.edge_junctions {
        support
            .entry(left)
            .and_modify(|max_count: &mut u32| *max_count = (*max_count).max(count))
            .or_insert(count);
        support
            .entry(right)
            .and_modify(|max_count: &mut u32| *max_count = (*max_count).max(count))
            .or_insert(count);
    }
    support
}

fn rank_keep_set(
    ranks: &mut HashMap<usize, Vec<(EdgeKey, u32)>>,
    max_edges_per_state: usize,
) -> HashSet<EdgeKey> {
    let mut keep = HashSet::new();
    if max_edges_per_state == 0 {
        return keep;
    }

    for edges in ranks.values_mut() {
        edges.sort_by(|(left_key, left_cov), (right_key, right_cov)| {
            right_cov
                .cmp(left_cov)
                .then_with(|| left_key.from.cmp(&right_key.from))
                .then_with(|| left_key.to.cmp(&right_key.to))
        });
        for &(edge, _) in edges.iter().take(max_edges_per_state) {
            keep.insert(edge);
        }
    }

    keep
}

// Step 02c unitig compression and GFA output helpers.
struct CompressedGraph {
    unitigs: Vec<Unitig>,
    edge_to_unitig: HashMap<EdgeKey, usize>,
    edge_to_placement: HashMap<EdgeKey, UnitigEnd>,
    state_starts: HashMap<usize, Vec<UnitigEnd>>,
    state_ends: HashMap<usize, Vec<UnitigEnd>>,
}

#[derive(Debug, Clone, Copy, Eq, PartialEq, Ord, PartialOrd, Hash)]
struct UnitigEnd {
    unitig_id: usize,
    orient: char,
}

#[derive(Debug, Clone, Copy)]
struct GfaLinkCandidate {
    left: UnitigEnd,
    right: UnitigEnd,
    support: u32,
    source: &'static str,
}

fn compress_unitigs(
    config: &Config,
    assembly: &Assembly,
    graph: &FilteredGraph,
) -> CompressedGraph {
    let mut visited: HashSet<EdgeKey> = HashSet::new();
    let mut directed_unitigs = Vec::new();

    let mut starts = Vec::new();
    for (edge, _) in &graph.edges {
        let indeg = graph.incoming.get(&edge.from).map_or(0, Vec::len);
        let outdeg = graph.outgoing.get(&edge.from).map_or(0, Vec::len);
        if indeg != 1 || outdeg != 1 {
            starts.push(*edge);
        }
    }

    for edge in starts {
        if visited.contains(&edge) {
            continue;
        }
        let unitig = trace_unitig(config, assembly, graph, edge, &mut visited);
        let unitig_id = directed_unitigs.len();
        directed_unitigs.push(Unitig {
            id: unitig_id,
            ..unitig
        });
    }

    for (edge, _) in &graph.edges {
        if visited.contains(edge) {
            continue;
        }
        let unitig = trace_unitig(config, assembly, graph, *edge, &mut visited);
        let unitig_id = directed_unitigs.len();
        directed_unitigs.push(Unitig {
            id: unitig_id,
            ..unitig
        });
    }

    let compressed = canonicalize_unitigs(directed_unitigs);
    let compressed = prune_redundant_unitigs(config, compressed);
    prune_short_tips(config, compressed)
}

fn canonicalize_unitigs(directed_unitigs: Vec<Unitig>) -> CompressedGraph {
    let mut path_index: HashMap<Vec<EdgeKey>, usize> = HashMap::new();
    for unitig in &directed_unitigs {
        path_index
            .entry(unitig.path_edges.clone())
            .or_insert(unitig.id);
    }

    let mut processed = vec![false; directed_unitigs.len()];
    let mut unitigs = Vec::new();
    let mut edge_to_unitig = HashMap::new();
    let mut edge_to_placement = HashMap::new();
    let mut state_starts: HashMap<usize, Vec<UnitigEnd>> = HashMap::new();
    let mut state_ends: HashMap<usize, Vec<UnitigEnd>> = HashMap::new();

    for old_id in 0..directed_unitigs.len() {
        if processed[old_id] {
            continue;
        }

        let old = &directed_unitigs[old_id];
        let rc_edges = reverse_complement_path(&old.path_edges);
        let rc_id = path_index.get(&rc_edges).copied();
        let use_old_forward = old.path_edges <= rc_edges;
        let segment_id = unitigs.len();
        let sequence = if use_old_forward {
            old.sequence.clone()
        } else {
            revcomp_string(&old.sequence)
        };
        let mut coverage_sum = old.coverage;
        let mut coverage_count = 1.0;

        processed[old_id] = true;
        let old_orient = if use_old_forward { '+' } else { '-' };
        register_unitig_projection(
            old,
            UnitigEnd {
                unitig_id: segment_id,
                orient: old_orient,
            },
            &mut edge_to_unitig,
            &mut edge_to_placement,
            &mut state_starts,
            &mut state_ends,
        );

        if let Some(pair_id) = rc_id {
            if pair_id != old_id {
                processed[pair_id] = true;
                let pair = &directed_unitigs[pair_id];
                coverage_sum += pair.coverage;
                coverage_count += 1.0;
                let pair_orient = if use_old_forward { '-' } else { '+' };
                register_unitig_projection(
                    pair,
                    UnitigEnd {
                        unitig_id: segment_id,
                        orient: pair_orient,
                    },
                    &mut edge_to_unitig,
                    &mut edge_to_placement,
                    &mut state_starts,
                    &mut state_ends,
                );
            }
        }

        unitigs.push(Unitig {
            id: segment_id,
            path_states: if use_old_forward {
                old.path_states.clone()
            } else {
                reverse_complement_states(&old.path_states)
            },
            path_edges: if use_old_forward {
                old.path_edges.clone()
            } else {
                rc_edges
            },
            sequence,
            coverage: coverage_sum / coverage_count,
        });
    }

    CompressedGraph {
        unitigs,
        edge_to_unitig,
        edge_to_placement,
        state_starts,
        state_ends,
    }
}

fn register_unitig_projection(
    unitig: &Unitig,
    placement: UnitigEnd,
    edge_to_unitig: &mut HashMap<EdgeKey, usize>,
    edge_to_placement: &mut HashMap<EdgeKey, UnitigEnd>,
    state_starts: &mut HashMap<usize, Vec<UnitigEnd>>,
    state_ends: &mut HashMap<usize, Vec<UnitigEnd>>,
) {
    for &edge in &unitig.path_edges {
        edge_to_unitig.insert(edge, placement.unitig_id);
        edge_to_placement.insert(edge, placement);
    }
    if let Some(&state) = unitig.path_states.first() {
        state_starts.entry(state).or_default().push(placement);
    }
    if let Some(&state) = unitig.path_states.last() {
        state_ends.entry(state).or_default().push(placement);
    }
}

fn reverse_complement_path(path_edges: &[EdgeKey]) -> Vec<EdgeKey> {
    path_edges
        .iter()
        .rev()
        .map(|edge| EdgeKey {
            from: reverse_state(edge.to),
            to: reverse_state(edge.from),
        })
        .collect()
}

fn reverse_complement_states(path_states: &[usize]) -> Vec<usize> {
    path_states
        .iter()
        .rev()
        .map(|&state| reverse_state(state))
        .collect()
}

fn reverse_state(state: usize) -> usize {
    state ^ 1
}

fn reverse_unitig_end(end: UnitigEnd) -> UnitigEnd {
    UnitigEnd {
        unitig_id: end.unitig_id,
        orient: if end.orient == '+' { '-' } else { '+' },
    }
}

fn prune_redundant_unitigs(config: &Config, graph: CompressedGraph) -> CompressedGraph {
    if config.min_unitig_len == 0
        && (config.dedup_kmer == 0 || config.containment_ratio <= 0.0 || graph.unitigs.len() < 2)
    {
        return graph;
    }

    let mut keep: Vec<bool> = graph
        .unitigs
        .iter()
        .map(|unitig| unitig.sequence.len() >= config.min_unitig_len)
        .collect();

    if config.dedup_kmer == 0 || config.containment_ratio <= 0.0 || graph.unitigs.len() < 2 {
        return rebuild_compressed_graph(graph, &keep);
    }

    let sketches: Vec<HashSet<u64>> = graph
        .unitigs
        .iter()
        .map(|unitig| sequence_kmer_set(&unitig.sequence, config.dedup_kmer))
        .collect();
    let mut order: Vec<usize> = (0..graph.unitigs.len()).collect();
    order.sort_by(|&left, &right| {
        graph.unitigs[right]
            .sequence
            .len()
            .cmp(&graph.unitigs[left].sequence.len())
            .then_with(|| left.cmp(&right))
    });

    for short_pos in 0..order.len() {
        let short_id = order[short_pos];
        if !keep[short_id] {
            continue;
        }
        if sketches[short_id].is_empty() {
            continue;
        }
        for &long_id in &order[..short_pos] {
            if !keep[long_id] || sketches[long_id].is_empty() {
                continue;
            }
            let shared = sketches[short_id]
                .iter()
                .filter(|hash| sketches[long_id].contains(hash))
                .count();
            let containment = shared as f64 / sketches[short_id].len() as f64;
            if containment >= config.containment_ratio {
                keep[short_id] = false;
                break;
            }
        }
    }

    rebuild_compressed_graph(graph, &keep)
}

fn prune_short_tips(config: &Config, graph: CompressedGraph) -> CompressedGraph {
    if config.min_tip_len == 0 || graph.unitigs.is_empty() {
        return graph;
    }

    let mut degrees = vec![0_usize; graph.unitigs.len()];
    let mut links = HashSet::new();
    for (state, end_ids) in &graph.state_ends {
        if let Some(start_ids) = graph.state_starts.get(state) {
            for &left in end_ids {
                for &right in start_ids {
                    if left.unitig_id == right.unitig_id && left.orient == right.orient {
                        continue;
                    }
                    if links.insert((left, right)) {
                        degrees[left.unitig_id] += 1;
                        degrees[right.unitig_id] += 1;
                    }
                }
            }
        }
    }

    let keep: Vec<bool> = graph
        .unitigs
        .iter()
        .map(|unitig| unitig.sequence.len() >= config.min_tip_len || degrees[unitig.id] > 1)
        .collect();
    rebuild_compressed_graph(graph, &keep)
}

fn sequence_kmer_set(sequence: &str, k: usize) -> HashSet<u64> {
    let seq = sequence.as_bytes();
    let mut hashes = HashSet::new();
    if k == 0 || seq.len() < k {
        return hashes;
    }
    for pos in 0..=seq.len() - k {
        if let Some(hash) = canonical_kmer_hash(seq, pos, k) {
            hashes.insert(hash);
        }
    }
    hashes
}

fn canonical_kmer_hash(seq: &[u8], pos: usize, k: usize) -> Option<u64> {
    if pos + k > seq.len() {
        return None;
    }
    let slice = &seq[pos..pos + k];
    if !slice.iter().all(|&base| base_code(base).is_some()) {
        return None;
    }
    let forward = stable_hash_bytes(slice);
    let reverse = stable_hash_revcomp(slice);
    Some(mix64(forward.min(reverse)))
}

fn rebuild_compressed_graph(graph: CompressedGraph, keep: &[bool]) -> CompressedGraph {
    let mut id_remap = HashMap::new();
    let mut unitigs = Vec::new();
    for unitig in graph.unitigs {
        if keep.get(unitig.id).copied().unwrap_or(false) {
            let new_id = unitigs.len();
            id_remap.insert(unitig.id, new_id);
            unitigs.push(Unitig {
                id: new_id,
                ..unitig
            });
        }
    }

    let edge_to_unitig = graph
        .edge_to_unitig
        .into_iter()
        .filter_map(|(edge, old_id)| id_remap.get(&old_id).map(|&new_id| (edge, new_id)))
        .collect();
    let edge_to_placement = graph
        .edge_to_placement
        .into_iter()
        .filter_map(|(edge, placement)| {
            id_remap.get(&placement.unitig_id).map(|&unitig_id| {
                (
                    edge,
                    UnitigEnd {
                        unitig_id,
                        orient: placement.orient,
                    },
                )
            })
        })
        .collect();
    let state_starts = remap_unitig_end_map(graph.state_starts, &id_remap);
    let state_ends = remap_unitig_end_map(graph.state_ends, &id_remap);

    CompressedGraph {
        unitigs,
        edge_to_unitig,
        edge_to_placement,
        state_starts,
        state_ends,
    }
}

fn remap_unitig_end_map(
    map: HashMap<usize, Vec<UnitigEnd>>,
    id_remap: &HashMap<usize, usize>,
) -> HashMap<usize, Vec<UnitigEnd>> {
    let mut remapped: HashMap<usize, Vec<UnitigEnd>> = HashMap::new();
    let mut seen = HashSet::new();
    for (state, placements) in map {
        for placement in placements {
            let Some(&unitig_id) = id_remap.get(&placement.unitig_id) else {
                continue;
            };
            let new_placement = UnitigEnd {
                unitig_id,
                orient: placement.orient,
            };
            if seen.insert((state, new_placement)) {
                remapped.entry(state).or_default().push(new_placement);
            }
        }
    }
    remapped
}

fn trace_unitig(
    config: &Config,
    assembly: &Assembly,
    graph: &FilteredGraph,
    first_edge: EdgeKey,
    visited: &mut HashSet<EdgeKey>,
) -> Unitig {
    let mut path_edges = Vec::new();
    let mut path_states = vec![first_edge.from];
    let mut sequence = oriented_anchor_seq(
        &assembly.nodes[node_id_from_state(first_edge.from)],
        first_edge.from,
    );
    let mut coverage_sum = 0_u64;
    let mut current = first_edge;

    loop {
        if visited.contains(&current) {
            break;
        }
        visited.insert(current);
        path_edges.push(current);
        path_states.push(current.to);

        if let Some(stats) = assembly.edges.get(&current) {
            sequence.push_str(&stats.suffix);
            coverage_sum += stats.coverage as u64;
        }

        let indeg = graph.incoming.get(&current.to).map_or(0, Vec::len);
        let out_edges = graph.outgoing.get(&current.to);
        let outdeg = out_edges.map_or(0, Vec::len);
        if indeg != 1 || outdeg != 1 {
            break;
        }
        let next = out_edges.expect("outdeg is 1")[0];
        if next == first_edge {
            break;
        }
        current = next;
    }

    let coverage = if path_edges.is_empty() {
        0.0
    } else {
        coverage_sum as f64 / path_edges.len() as f64
    };

    Unitig {
        id: usize::MAX,
        path_states,
        path_edges,
        sequence: if sequence.is_empty() {
            "N".repeat(config.k)
        } else {
            sequence
        },
        coverage,
    }
}

fn count_unitig_junctions(
    assembly: &Assembly,
    edge_to_unitig: &HashMap<EdgeKey, usize>,
) -> HashMap<(usize, usize), u32> {
    let mut counts = HashMap::new();
    for walk in &assembly.read_walks {
        let mut previous = None;
        for edge in &walk.edges {
            let Some(&unitig_id) = edge_to_unitig.get(edge) else {
                continue;
            };
            if previous != Some(unitig_id) {
                if let Some(prev_id) = previous {
                    *counts.entry((prev_id, unitig_id)).or_insert(0) += 1;
                }
                previous = Some(unitig_id);
            }
        }
    }
    counts
}

fn count_link_support(
    assembly: &Assembly,
    edge_to_placement: &HashMap<EdgeKey, UnitigEnd>,
) -> HashMap<(UnitigEnd, UnitigEnd), u32> {
    let mut counts = HashMap::new();
    for walk in &assembly.read_walks {
        let mut previous = None;
        for edge in &walk.edges {
            let Some(&placement) = edge_to_placement.get(edge) else {
                continue;
            };
            if previous != Some(placement) {
                if let Some(prev_placement) = previous {
                    *counts.entry((prev_placement, placement)).or_insert(0) += 1;
                }
                previous = Some(placement);
            }
        }
    }
    counts
}

fn write_anchors(config: &Config, assembly: &Assembly, graph: &FilteredGraph) -> io::Result<()> {
    let mut out = File::create(config.out_dir.join("anchors.tsv"))?;
    writeln!(
        out,
        "anchor_id\tcoverage\tforward_indegree\tforward_outdegree\treverse_indegree\treverse_outdegree\tsequence"
    )?;
    for (id, node) in assembly.nodes.iter().enumerate() {
        if node.coverage < config.min_anchor_coverage {
            continue;
        }
        let fwd = state_id(id, true);
        let rev = state_id(id, false);
        writeln!(
            out,
            "a{id}\t{}\t{}\t{}\t{}\t{}\t{}",
            node.coverage,
            graph.incoming.get(&fwd).map_or(0, Vec::len),
            graph.outgoing.get(&fwd).map_or(0, Vec::len),
            graph.incoming.get(&rev).map_or(0, Vec::len),
            graph.outgoing.get(&rev).map_or(0, Vec::len),
            node.seq
        )?;
    }
    Ok(())
}

fn write_edges(
    config: &Config,
    graph: &FilteredGraph,
    junction_support: &HashMap<EdgeKey, u32>,
    edge_to_unitig: &HashMap<EdgeKey, usize>,
) -> io::Result<()> {
    let mut out = File::create(config.out_dir.join("edges.tsv"))?;
    writeln!(
        out,
        "edge_id\tfrom_anchor\tfrom_orient\tto_anchor\tto_orient\tcoverage\tmax_raw_junction_support\tmean_span\tmin_span\tmax_span\tunitig"
    )?;
    for (i, (edge, stats)) in graph.edges.iter().enumerate() {
        let mean_span = stats.total_span as f64 / stats.coverage.max(1) as f64;
        let raw_junction_support = junction_support.get(edge).copied().unwrap_or(0);
        let unitig = edge_to_unitig
            .get(edge)
            .map(|id| format!("utg{id}"))
            .unwrap_or_else(|| ".".to_string());
        writeln!(
            out,
            "e{i}\ta{}\t{}\ta{}\t{}\t{}\t{}\t{mean_span:.2}\t{}\t{}\t{}",
            node_id_from_state(edge.from),
            state_orientation(edge.from),
            node_id_from_state(edge.to),
            state_orientation(edge.to),
            stats.coverage,
            raw_junction_support,
            stats.min_span,
            stats.max_span,
            unitig
        )?;
    }
    Ok(())
}

fn write_unitigs(config: &Config, unitigs: &[Unitig], filename: &str) -> io::Result<()> {
    let mut out = File::create(config.out_dir.join(filename))?;
    for unitig in unitigs {
        writeln!(
            out,
            ">utg{} len={} edges={} cov={:.2}",
            unitig.id,
            unitig.sequence.len(),
            unitig.path_edges.len(),
            unitig.coverage
        )?;
        write_wrapped(&mut out, &unitig.sequence, 80)?;
    }
    Ok(())
}

fn write_wrapped(out: &mut File, sequence: &str, width: usize) -> io::Result<()> {
    for chunk in sequence.as_bytes().chunks(width) {
        out.write_all(chunk)?;
        out.write_all(b"\n")?;
    }
    Ok(())
}

fn write_gfa(
    config: &Config,
    unitigs: &[Unitig],
    state_starts: &HashMap<usize, Vec<UnitigEnd>>,
    state_ends: &HashMap<usize, Vec<UnitigEnd>>,
    link_support: &HashMap<(UnitigEnd, UnitigEnd), u32>,
    filename: &str,
) -> io::Result<()> {
    let mut out = File::create(config.out_dir.join(filename))?;
    writeln!(out, "H\tVN:Z:1.0")?;
    for unitig in unitigs {
        writeln!(
            out,
            "S\tutg{}\t{}\tLN:i:{}\tRC:f:{:.2}",
            unitig.id,
            unitig.sequence,
            unitig.sequence.len(),
            unitig.coverage
        )?;
    }

    let candidates = collect_gfa_link_candidates(config, state_starts, state_ends, link_support);
    let (out_best, in_best) = gfa_link_best_support(&candidates);
    let mut seen_links = HashSet::new();
    let mut links = Vec::new();
    for candidate in candidates {
        let local_out_best = out_best.get(&candidate.left).copied().unwrap_or(0);
        let local_in_best = in_best.get(&candidate.right).copied().unwrap_or(0);
        if !keep_gfa_link(config, candidate.support, local_out_best, local_in_best) {
            continue;
        }
        add_gfa_link(
            &mut links,
            &mut seen_links,
            candidate.left,
            candidate.right,
            candidate.support,
            candidate.source,
            config.bidirectional_links,
        );
    }

    links.sort_by_key(|(left, right, _support, source)| (*left, *right, *source));
    for (left, right, support, source) in links {
        if config.read_junction_links {
            writeln!(
                out,
                "L\tutg{}\t{}\tutg{}\t{}\t0M\tRC:i:{}\tJL:Z:{}",
                left.unitig_id, left.orient, right.unitig_id, right.orient, support, source
            )?;
        } else {
            writeln!(
                out,
                "L\tutg{}\t{}\tutg{}\t{}\t0M\tRC:i:{}",
                left.unitig_id, left.orient, right.unitig_id, right.orient, support
            )?;
        }
    }
    Ok(())
}

fn collect_gfa_link_candidates(
    config: &Config,
    state_starts: &HashMap<usize, Vec<UnitigEnd>>,
    state_ends: &HashMap<usize, Vec<UnitigEnd>>,
    link_support: &HashMap<(UnitigEnd, UnitigEnd), u32>,
) -> Vec<GfaLinkCandidate> {
    let mut candidates = Vec::new();
    for (state, end_ids) in state_ends {
        if let Some(start_ids) = state_starts.get(state) {
            for &left in end_ids {
                for &right in start_ids {
                    if left.unitig_id == right.unitig_id && left.orient == right.orient {
                        continue;
                    }
                    candidates.push(GfaLinkCandidate {
                        left,
                        right,
                        support: link_support.get(&(left, right)).copied().unwrap_or(0),
                        source: "terminal",
                    });
                }
            }
        }
    }

    if config.read_junction_links {
        for (&(left, right), &support) in link_support {
            if left.unitig_id == right.unitig_id && left.orient == right.orient {
                continue;
            }
            candidates.push(GfaLinkCandidate {
                left,
                right,
                support,
                source: "read",
            });
        }
    }
    candidates
}

fn gfa_link_best_support(
    candidates: &[GfaLinkCandidate],
) -> (HashMap<UnitigEnd, u32>, HashMap<UnitigEnd, u32>) {
    let mut out_best: HashMap<UnitigEnd, u32> = HashMap::new();
    let mut in_best: HashMap<UnitigEnd, u32> = HashMap::new();
    for candidate in candidates {
        out_best
            .entry(candidate.left)
            .and_modify(|best| *best = (*best).max(candidate.support))
            .or_insert(candidate.support);
        in_best
            .entry(candidate.right)
            .and_modify(|best| *best = (*best).max(candidate.support))
            .or_insert(candidate.support);
    }
    (out_best, in_best)
}

fn keep_gfa_link(config: &Config, support: u32, out_best: u32, in_best: u32) -> bool {
    if support < config.min_link_support {
        return false;
    }
    if config.min_link_ratio <= 0.0 {
        return true;
    }
    if support == 0 {
        return false;
    }
    let out_ratio = support as f64 / out_best.max(1) as f64;
    let in_ratio = support as f64 / in_best.max(1) as f64;
    out_ratio.min(in_ratio) >= config.min_link_ratio
}

fn add_gfa_link(
    links: &mut Vec<(UnitigEnd, UnitigEnd, u32, &'static str)>,
    seen_links: &mut HashSet<(UnitigEnd, UnitigEnd)>,
    left: UnitigEnd,
    right: UnitigEnd,
    support: u32,
    source: &'static str,
    bidirectional: bool,
) {
    if seen_links.insert((left, right)) {
        links.push((left, right, support, source));
    }
    if bidirectional {
        let rc_left = reverse_unitig_end(right);
        let rc_right = reverse_unitig_end(left);
        if seen_links.insert((rc_left, rc_right)) {
            links.push((rc_left, rc_right, support, source));
        }
    }
}

fn write_junctions(
    config: &Config,
    junctions: &HashMap<(usize, usize), u32>,
    filename: &str,
) -> io::Result<()> {
    let mut out_totals: HashMap<usize, u32> = HashMap::new();
    for (&(from, _to), &count) in junctions {
        *out_totals.entry(from).or_insert(0) += count;
    }

    let mut rows: Vec<_> = junctions.iter().collect();
    rows.sort_by_key(|(&(from, to), _)| (from, to));

    let mut out = File::create(config.out_dir.join(filename))?;
    writeln!(out, "from_unitig\tto_unitig\tcount\tfrom_total\tfrequency")?;
    for (&(from, to), &count) in rows {
        let total = out_totals.get(&from).copied().unwrap_or(0);
        let freq = if total == 0 {
            0.0
        } else {
            count as f64 / total as f64
        };
        writeln!(out, "utg{from}\tutg{to}\t{count}\t{total}\t{freq:.6}")?;
    }
    Ok(())
}

fn write_report(
    config: &Config,
    assembly: &Assembly,
    graph: &FilteredGraph,
    full_graph: &FilteredGraph,
    unitigs: &[Unitig],
    full_unitigs: &[Unitig],
    elapsed: std::time::Duration,
) -> io::Result<()> {
    let mut out = File::create(config.out_dir.join("report.txt"))?;
    writeln!(out, "orgraft_asm report")?;
    writeln!(out, "elapsed_seconds\t{:.3}", elapsed.as_secs_f64())?;
    writeln!(out, "reads_seen\t{}", assembly.reads_seen)?;
    writeln!(out, "bases_seen\t{}", assembly.bases_seen)?;
    writeln!(out, "raw_anchor_hits\t{}", assembly.anchors_seen)?;
    writeln!(out, "unique_anchors\t{}", assembly.nodes.len())?;
    writeln!(out, "raw_edges\t{}", assembly.edges.len())?;
    writeln!(out, "raw_edge_junctions\t{}", assembly.edge_junctions.len())?;
    writeln!(out, "filtered_edges\t{}", graph.edges.len())?;
    writeln!(out, "full_edges\t{}", full_graph.edges.len())?;
    writeln!(out, "unitigs\t{}", unitigs.len())?;
    writeln!(
        out,
        "unitig_bases\t{}",
        unitigs
            .iter()
            .map(|unitig| unitig.sequence.len())
            .sum::<usize>()
    )?;
    writeln!(out, "full_unitigs\t{}", full_unitigs.len())?;
    writeln!(
        out,
        "full_unitig_bases\t{}",
        full_unitigs
            .iter()
            .map(|unitig| unitig.sequence.len())
            .sum::<usize>()
    )?;
    writeln!(
        out,
        "organelle\t{}",
        config
            .organelle
            .map(|o| o.as_str())
            .unwrap_or("unspecified")
    )?;
    writeln!(
        out,
        "numt_interference\t{}",
        config.numt_interference.as_str()
    )?;
    writeln!(out, "data_mode\t{}", config.data_mode.as_str())?;
    writeln!(out, "graph_workflow\t{}", graph_workflow_label(config))?;
    writeln!(out, "anchor_mode\tsyncmer")?;
    writeln!(out, "k\t{}", config.k)?;
    writeln!(out, "s\t{}", config.s)?;
    writeln!(
        out,
        "syncmer_pos\t{}",
        config.syncmer_pos.unwrap_or((config.k - config.s) / 2)
    )?;
    writeln!(out, "min_anchor_spacing\t{}", config.min_anchor_spacing)?;
    writeln!(out, "min_anchor_coverage\t{}", config.min_anchor_coverage)?;
    writeln!(out, "min_edge_coverage\t{}", config.min_edge_coverage)?;
    writeln!(out, "min_branch_ratio\t{}", config.min_branch_ratio)?;
    writeln!(out, "max_edges_per_state\t{}", config.max_edges_per_state)?;
    writeln!(out, "dedup_kmer\t{}", config.dedup_kmer)?;
    writeln!(out, "containment_ratio\t{}", config.containment_ratio)?;
    writeln!(out, "min_unitig_len\t{}", config.min_unitig_len)?;
    writeln!(out, "min_tip_len\t{}", config.min_tip_len)?;
    writeln!(out, "min_link_support\t{}", config.min_link_support)?;
    writeln!(out, "min_link_ratio\t{}", config.min_link_ratio)?;
    writeln!(out, "read_subset_percent\t{}", config.read_subset_percent)?;
    writeln!(out, "read_junction_links\t{}", config.read_junction_links)?;
    writeln!(out, "bidirectional_links\t{}", config.bidirectional_links)?;
    writeln!(
        out,
        "junction_rescue_support\t{}",
        config.junction_rescue_support
    )?;
    if let Some(genome_size) = config.genome_size {
        writeln!(out, "genome_size\t{genome_size}")?;
    }
    if let Some(coverage) = config.asm_coverage {
        writeln!(out, "asm_coverage\t{coverage}")?;
    }
    if let Some(min_overlap) = config.min_overlap {
        writeln!(out, "min_overlap\t{min_overlap}")?;
    }
    if let Some(iterations) = config.iterations {
        writeln!(out, "iterations\t{iterations}")?;
    }
    writeln!(out, "hifi_error_rate\t{:.5}", config.hifi_error_rate)?;
    writeln!(
        out,
        "expected_exact_anchor_survival\t{:.6}",
        expected_exact_survival(config.k, config.hifi_error_rate)
    )?;
    writeln!(
        out,
        "expected_exact_anchor_pair_survival\t{:.6}",
        expected_exact_survival(config.k, config.hifi_error_rate).powi(2)
    )?;
    writeln!(out, "minimap2\t{}", config.run_minimap2)?;
    writeln!(out, "minimap_min_identity\t{}", config.minimap_min_identity)?;
    writeln!(
        out,
        "minimap_min_align_len\t{}",
        config.minimap_min_align_len
    )?;
    writeln!(out, "paf_max_link_gap\t{}", config.paf_max_link_gap)?;
    writeln!(out, "threads\t{}", config.threads)?;
    writeln!(out, "keep_debug_files\t{}", config.keep_debug_files)?;
    writeln!(out, "graph_gfa\tgraph.gfa")?;
    if config.keep_debug_files {
        writeln!(out, "graph_full_gfa\tgraph.full.gfa")?;
    }
    writeln!(
        out,
        "support_note\tRC:i is observed exact-anchor/read-walk support; long-k HiFi anchor dropout can underestimate true support, so use minimap2 depth/link TSVs when --minimap2 is enabled."
    )?;
    Ok(())
}

fn expected_exact_survival(k: usize, error_rate: f64) -> f64 {
    (1.0 - error_rate).clamp(0.0, 1.0).powi(k as i32)
}

#[derive(Debug, Clone)]
struct PafAln {
    qname: String,
    qstart: usize,
    qend: usize,
    strand: char,
    tname: String,
    tlen: usize,
    tstart: usize,
    tend: usize,
    matches: usize,
    alen: usize,
}

impl PafAln {
    fn identity(&self) -> f64 {
        if self.alen == 0 {
            0.0
        } else {
            self.matches as f64 / self.alen as f64
        }
    }
}

// Step 02d skeleton linking: use the first-round anchor graph as a skeleton,
// align reads back to unitigs, and rescue graph links from read evidence.
#[derive(Debug, Clone)]
struct SkeletonSegment {
    name: String,
    sequence: String,
}

#[derive(Debug, Clone, Eq, PartialEq, Ord, PartialOrd, Hash)]
struct SkeletonLinkKey {
    from: String,
    from_orient: char,
    to: String,
    to_orient: char,
}

#[derive(Debug, Clone, Default)]
struct SkeletonLinkSupport {
    skeleton_support: u32,
    paf_support: u32,
    rescue_support: u32,
    out_ratio: f64,
    in_ratio: f64,
}

const MITO_COMPACT_LOCAL_BRIDGE_GAP_MULTIPLIER: usize = 20;
const MITO_COMPACT_MAX_ENDPOINT_DEGREE: usize = 2;
const MITO_COMPACT_SMALL_COMPONENT_MAX_SEGMENTS: usize = 3;
const MITO_COMPACT_SMALL_COMPONENT_MAX_BASES: usize = 30_000;

#[derive(Debug, Clone)]
struct SkeletonGraphStats {
    components: Vec<Vec<String>>,
    component_by_segment: HashMap<String, usize>,
    open_out: HashSet<(String, char)>,
    open_in: HashSet<(String, char)>,
    out_degree: HashMap<(String, char), usize>,
    in_degree: HashMap<(String, char), usize>,
    kept_links: usize,
}

#[derive(Debug, Clone, Copy, Eq, PartialEq, Ord, PartialOrd)]
struct MitoStableTopologyStats {
    components: usize,
    open_ends: usize,
    endpoint_overload: usize,
    unique_links: usize,
    branch_ends: usize,
    cycle_rank: usize,
}

impl MitoStableTopologyStats {
    fn score_key(
        self,
        selected_links: usize,
        support_sum: u64,
    ) -> (
        usize,
        usize,
        usize,
        usize,
        usize,
        usize,
        usize,
        std::cmp::Reverse<u64>,
    ) {
        (
            self.components,
            self.endpoint_overload,
            self.open_ends,
            self.cycle_rank,
            self.unique_links,
            self.branch_ends,
            selected_links,
            std::cmp::Reverse(support_sum),
        )
    }
}

#[derive(Debug, Clone, Copy, Eq, PartialEq)]
enum MitoCompactBridgeSource {
    LocalPaf,
    Rescue,
    CopyChoice,
}

impl MitoCompactBridgeSource {
    fn label(self) -> &'static str {
        match self {
            MitoCompactBridgeSource::LocalPaf => "local_paf",
            MitoCompactBridgeSource::Rescue => "read_walk",
            MitoCompactBridgeSource::CopyChoice => "copy_choice",
        }
    }

    fn priority(self) -> usize {
        match self {
            MitoCompactBridgeSource::LocalPaf => 1,
            MitoCompactBridgeSource::Rescue => 0,
            MitoCompactBridgeSource::CopyChoice => 2,
        }
    }
}

#[derive(Debug, Clone)]
struct MitoCompactBridgeCandidate {
    key: SkeletonLinkKey,
    support: u32,
    source: MitoCompactBridgeSource,
}

#[derive(Debug, Clone)]
struct MitoCompactSelectedBridge {
    candidate: MitoCompactBridgeCandidate,
    component_delta: usize,
    open_delta: usize,
}

#[derive(Debug, Default)]
struct MitoCompactBridgeReport {
    initial_components: Vec<usize>,
    initial_open_roles: usize,
    initial_kept_links: usize,
    focused_reads: usize,
    local_candidates: usize,
    rescue_candidates: usize,
    selected: Vec<MitoCompactSelectedBridge>,
    pruned_links: Vec<SkeletonLinkKey>,
    pruned_segments: Vec<String>,
    final_components: Vec<usize>,
    final_open_roles: usize,
    final_kept_links: usize,
}

#[derive(Debug, Clone)]
struct MitoStableSplitPoint {
    position: usize,
    support: u32,
}

#[derive(Debug, Clone)]
struct MitoStableSplitReportRow {
    segment: String,
    position: usize,
    support: u32,
}

#[derive(Debug, Clone)]
struct MitoStableSplitLayout {
    pieces: Vec<String>,
}

#[derive(Debug, Clone, Eq, PartialEq, Ord, PartialOrd)]
struct MitoStableRepeatExpansion {
    template: String,
    clone: String,
    removed_shortcut: SkeletonLinkKey,
    left_link: SkeletonLinkKey,
    right_link: SkeletonLinkKey,
    shortcut_overlap: usize,
    left_overlap: usize,
    right_overlap: usize,
    support: u32,
}

#[derive(Debug, Clone)]
struct MitoStableManualEdit {
    action: String,
    key: SkeletonLinkKey,
}

#[derive(Debug, Clone)]
struct MitoStableCopyChoiceRow {
    node_a: String,
    node_b: String,
    option: String,
    removed_links: Vec<SkeletonLinkKey>,
    added_link: Option<SkeletonLinkKey>,
    added_support: u32,
    added_support_source: String,
    components: usize,
    open_ends: usize,
    endpoint_overload: usize,
    branch_ends: usize,
    cycle_rank: usize,
    node_a_class: &'static str,
    node_b_class: &'static str,
    interpretation: String,
}

#[derive(Debug, Clone)]
struct MitoStableSelectedNodeDiagnostic {
    segment: String,
    left_degree: usize,
    right_degree: usize,
    class: &'static str,
    selected: bool,
    mode: String,
    incident_base_links: Vec<String>,
    incident_candidates: Vec<String>,
}

fn run_skeleton_linking(config: &Config) -> io::Result<()> {
    let skeleton_gfa = config
        .skeleton_gfa
        .as_ref()
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidInput, "--skeleton-gfa is required"))?;
    let (segments, skeleton_links) = read_skeleton_gfa(skeleton_gfa)?;
    if segments.is_empty() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!("{} has no S records", skeleton_gfa.display()),
        ));
    }

    fs::create_dir_all(&config.out_dir)?;
    let skeleton_fasta = config.out_dir.join("skeleton_segments.fasta");
    write_skeleton_fasta(&skeleton_fasta, &segments)?;

    let paf_path = config.out_dir.join("read_alignments.paf");
    let log_path = config.out_dir.join("read_alignments.minimap2.log");
    run_minimap2_to_target(config, &skeleton_fasta, &paf_path, &log_path)?;

    let (depth, paf_links) = summarize_skeleton_paf(config, &segments, &paf_path)?;
    write_skeleton_depth(&config.out_dir.join("depth.tsv"), &segments, &depth)?;
    let mut merged_links = merge_skeleton_links(config, skeleton_links, paf_links);
    if let Some(rescue_gfa) = &config.skeleton_rescue_gfa {
        let rescue_links = read_gfa_links(rescue_gfa)?;
        add_component_rescue_links(config, &segments, &mut merged_links, rescue_links);
    }
    write_skeleton_links(&config.out_dir.join("links.tsv"), &merged_links)?;
    write_skeleton_linked_gfa(
        config,
        &config.out_dir.join("graph.gfa"),
        &segments,
        &depth,
        &merged_links,
    )?;
    write_skeleton_report(config, skeleton_gfa, &segments, &merged_links)?;
    Ok(())
}

fn run_mito_compact_bridge_linking(config: &Config) -> io::Result<()> {
    let skeleton_gfa = config
        .skeleton_gfa
        .as_ref()
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidInput, "--skeleton-gfa is required"))?;
    let (segments, skeleton_links) = read_skeleton_gfa(skeleton_gfa)?;
    if segments.is_empty() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!("{} has no S records", skeleton_gfa.display()),
        ));
    }

    fs::create_dir_all(&config.out_dir)?;
    let skeleton_fasta = config.out_dir.join("skeleton_segments.fasta");
    write_skeleton_fasta(&skeleton_fasta, &segments)?;

    let paf_path = config.out_dir.join("read_alignments.paf");
    let log_path = config.out_dir.join("read_alignments.minimap2.log");
    run_minimap2_to_target(config, &skeleton_fasta, &paf_path, &log_path)?;

    let (depth, paf_links) = summarize_skeleton_paf(config, &segments, &paf_path)?;
    let paf_by_read = read_filtered_skeleton_paf_by_read(config, &segments, &paf_path)?;
    let mut merged_links = merge_skeleton_links(config, skeleton_links, paf_links);
    let rescue_links = match &config.skeleton_rescue_gfa {
        Some(rescue_gfa) => read_gfa_links(rescue_gfa)?,
        None => HashMap::new(),
    };

    let mut bridge_report = complete_mito_compact_links(
        config,
        &segments,
        &mut merged_links,
        rescue_links,
        &paf_by_read,
    )?;
    refresh_skeleton_link_ratios(config, &mut merged_links);
    let pruned_links = prune_mito_compact_secondary_links(config, &segments, &mut merged_links);
    bridge_report.pruned_links = pruned_links;
    refresh_skeleton_link_ratios(config, &mut merged_links);

    let final_segment_names =
        mito_compact_final_segment_names(config, &segments, &merged_links, &mut bridge_report);
    let final_segments: Vec<SkeletonSegment> = segments
        .iter()
        .filter(|segment| final_segment_names.contains(&segment.name))
        .cloned()
        .collect();
    let final_links = filter_skeleton_links_to_segments(&merged_links, &final_segment_names);
    let final_stats = skeleton_graph_stats(config, &final_segments, &final_links);
    bridge_report.final_components = final_stats.components.iter().map(Vec::len).collect();
    bridge_report.final_open_roles = final_stats.open_out.len() + final_stats.open_in.len();
    bridge_report.final_kept_links = final_stats.kept_links;

    write_skeleton_depth(&config.out_dir.join("depth.tsv"), &final_segments, &depth)?;
    write_skeleton_links(&config.out_dir.join("links.tsv"), &final_links)?;
    write_skeleton_linked_gfa(
        config,
        &config.out_dir.join("graph.gfa"),
        &final_segments,
        &depth,
        &final_links,
    )?;
    write_skeleton_report(config, skeleton_gfa, &final_segments, &final_links)?;
    write_mito_compact_bridge_report(config, &bridge_report)?;
    Ok(())
}

fn run_mito_stable_linking(config: &Config, resolution_dir: Option<&Path>) -> io::Result<()> {
    let skeleton_gfa = config
        .skeleton_gfa
        .as_ref()
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidInput, "--skeleton-gfa is required"))?;
    let (mut segments, mut skeleton_links) = read_skeleton_gfa(skeleton_gfa)?;
    if segments.is_empty() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!("{} has no S records", skeleton_gfa.display()),
        ));
    }

    fs::create_dir_all(&config.out_dir)?;
    let mut resolution_config = config.clone();
    if let Some(resolution_dir) = resolution_dir {
        resolution_config.out_dir = resolution_dir.to_path_buf();
    }
    fs::create_dir_all(&resolution_config.out_dir)?;

    let skeleton_fasta = config.out_dir.join("skeleton_segments.fasta");
    write_skeleton_fasta(&skeleton_fasta, &segments)?;

    let paf_path = config.out_dir.join("read_alignments.paf");
    let log_path = config.out_dir.join("read_alignments.minimap2.log");
    run_minimap2_to_target(config, &skeleton_fasta, &paf_path, &log_path)?;

    let presplit_paf_by_read = read_filtered_skeleton_paf_by_read(config, &segments, &paf_path)?;
    let split_points =
        detect_mito_stable_internal_split_points(config, &segments, &presplit_paf_by_read);
    let mut split_report = Vec::new();
    if !split_points.is_empty() {
        let presplit_fasta = resolution_config
            .out_dir
            .join("presplit_skeleton_segments.fasta");
        let presplit_paf = resolution_config
            .out_dir
            .join("read_alignments.presplit.paf");
        let presplit_log = resolution_config
            .out_dir
            .join("read_alignments.presplit.minimap2.log");
        fs::copy(&skeleton_fasta, presplit_fasta)?;
        fs::copy(&paf_path, presplit_paf)?;
        if log_path.exists() {
            fs::copy(&log_path, presplit_log)?;
        }
        let (refined_segments, refined_links, refined_report) =
            apply_mito_stable_internal_splits(config, segments, skeleton_links, split_points);
        segments = refined_segments;
        skeleton_links = refined_links;
        split_report = refined_report;
        write_skeleton_fasta(&skeleton_fasta, &segments)?;
        run_minimap2_to_target(config, &skeleton_fasta, &paf_path, &log_path)?;
    }

    let (mut depth, paf_links) = summarize_skeleton_paf(config, &segments, &paf_path)?;
    let paf_by_read = read_filtered_skeleton_paf_by_read(config, &segments, &paf_path)?;
    let rescue_links = match &config.skeleton_rescue_gfa {
        Some(rescue_gfa) => read_gfa_links(rescue_gfa)?,
        None => HashMap::new(),
    };

    let (mut base_links, mut candidates) =
        mito_stable_base_links_and_paf_candidates(skeleton_links, paf_links);
    refresh_skeleton_link_ratios(config, &mut base_links);
    write_skeleton_depth(&config.out_dir.join("depth.tsv"), &segments, &depth)?;
    write_skeleton_links(&config.out_dir.join("links.tsv"), &base_links)?;
    write_skeleton_linked_gfa(
        config,
        &config.out_dir.join("graph.gfa"),
        &segments,
        &depth,
        &base_links,
    )?;
    write_skeleton_report(config, skeleton_gfa, &segments, &base_links)?;
    let initial_stats = skeleton_graph_stats(config, &segments, &base_links);
    candidates
        .retain(|candidate| mito_compact_candidate_is_relevant(&initial_stats, &candidate.key));

    let focused_reads = mito_compact_focused_reads(&initial_stats, &paf_by_read);
    write_mito_compact_focused_reads(&resolution_config, &focused_reads)?;
    candidates.extend(mito_compact_local_paf_candidates(
        config,
        &initial_stats,
        &paf_by_read,
        &focused_reads,
    ));
    candidates.extend(mito_compact_rescue_candidates(&initial_stats, rescue_links));
    let all_candidates = dedupe_mito_bridge_candidates(candidates);
    let all_candidate_count = all_candidates.len();
    let copy_choice_rows =
        analyze_mito_stable_copy_choices(config, &segments, &base_links, &all_candidates);
    let selected_node_diagnostics =
        analyze_mito_stable_selected_nodes(config, &segments, &base_links, &all_candidates);
    let candidates = mito_stable_scoped_candidates(config, all_candidates.clone());

    let base_topology = mito_stable_topology_stats(config, &segments, &base_links);
    let copy_choice_bridges =
        select_mito_stable_copy_choice_bridges(config, &base_topology, &copy_choice_rows);
    let mut default_node_shape_violations = 0usize;
    let mut node_constraint_fallback = false;
    let selected = match config.mito_stable_selection_mode {
        MitoStableSelectionMode::Auto => {
            let mut selected =
                select_mito_stable_bridges(config, &segments, &base_links, &candidates, false);
            append_mito_stable_selected_bridges(&mut selected, copy_choice_bridges.clone());
            if mito_stable_has_explicit_three_way_constraints(config) {
                let mut trial_links = base_links.clone();
                for candidate in &selected {
                    insert_mito_compact_bridge(config, &mut trial_links, candidate);
                }
                refresh_skeleton_link_ratios(config, &mut trial_links);
                default_node_shape_violations =
                    mito_stable_node_shape_constraint_violations(config, &segments, &trial_links);
                if default_node_shape_violations > 0 {
                    node_constraint_fallback = true;
                    selected = select_mito_stable_bridges(
                        config,
                        &segments,
                        &base_links,
                        &candidates,
                        true,
                    );
                }
            }
            selected
        }
        #[cfg(test)]
        MitoStableSelectionMode::ForbidSelected => {
            default_node_shape_violations =
                mito_stable_node_shape_constraint_violations(config, &segments, &base_links);
            let mut selected =
                select_mito_stable_bridges(config, &segments, &base_links, &candidates, true);
            append_mito_stable_selected_bridges(&mut selected, copy_choice_bridges.clone());
            selected
        }
        MitoStableSelectionMode::AllowSelected => Vec::new(),
    };
    let mut final_links = base_links.clone();
    for candidate in &selected {
        insert_mito_compact_bridge(config, &mut final_links, candidate);
    }
    let manual_edits = apply_mito_stable_manual_edits(config, &mut final_links);
    refresh_skeleton_link_ratios(config, &mut final_links);
    let repeat_expansions = if config.mito_stable_selection_mode == MitoStableSelectionMode::Auto {
        let expansions = expand_mito_stable_repeat_shortcuts(
            config,
            &mut segments,
            &mut depth,
            &mut final_links,
            &candidates,
        );
        refresh_skeleton_link_ratios(config, &mut final_links);
        expansions
    } else {
        Vec::new()
    };
    let pruned_links = if config.mito_stable_selection_mode == MitoStableSelectionMode::Auto {
        prune_mito_stable_redundant_links(config, &segments, &mut final_links)
    } else {
        Vec::new()
    };
    refresh_skeleton_link_ratios(config, &mut final_links);
    let final_topology = mito_stable_topology_stats(config, &segments, &final_links);
    let final_node_shape_violations =
        mito_stable_node_shape_constraint_violations(config, &segments, &final_links);

    write_skeleton_depth(
        &resolution_config.out_dir.join("depth.tsv"),
        &segments,
        &depth,
    )?;
    write_skeleton_links(&resolution_config.out_dir.join("links.tsv"), &final_links)?;
    write_skeleton_linked_gfa(
        &resolution_config,
        &resolution_config.out_dir.join("graph.gfa"),
        &segments,
        &depth,
        &final_links,
    )?;
    write_skeleton_report(&resolution_config, skeleton_gfa, &segments, &final_links)?;
    write_mito_stable_bridge_report(
        &resolution_config,
        &base_topology,
        &final_topology,
        &candidates,
        all_candidate_count,
        &copy_choice_bridges,
        &selected,
        &pruned_links,
        &repeat_expansions,
        &manual_edits,
        default_node_shape_violations,
        node_constraint_fallback,
        final_node_shape_violations,
    )?;
    write_mito_stable_repeat_expansions(&resolution_config, &repeat_expansions)?;
    write_mito_stable_pruned_links(&resolution_config, &pruned_links)?;
    write_mito_stable_manual_edits(&resolution_config, &manual_edits)?;
    write_mito_stable_copy_choices(&resolution_config, &copy_choice_rows)?;
    write_mito_stable_selected_node_diagnostics(&resolution_config, &selected_node_diagnostics)?;
    write_mito_stable_node_degrees(&resolution_config, &segments, &final_links)?;
    write_mito_stable_node_repairs(
        &resolution_config,
        &segments,
        &base_links,
        &final_links,
        &selected,
        &manual_edits,
        &pruned_links,
        &repeat_expansions,
    )?;
    write_mito_stable_topology_scan(&resolution_config, &segments, &final_links)?;
    write_mito_stable_split_report(&resolution_config, &split_report)?;
    Ok(())
}

#[derive(Debug, Clone)]
struct MitoStablePafLocus {
    segment: String,
    coord: usize,
    tlen: usize,
}

fn detect_mito_stable_internal_split_points(
    config: &Config,
    segments: &[SkeletonSegment],
    paf_by_read: &HashMap<String, Vec<PafAln>>,
) -> HashMap<String, Vec<MitoStableSplitPoint>> {
    let segment_lengths: HashMap<String, usize> = segments
        .iter()
        .map(|segment| (segment.name.clone(), segment.sequence.len()))
        .collect();
    let mut events: HashMap<String, Vec<usize>> = HashMap::new();
    let max_gap = config.paf_max_link_gap;
    for alns in paf_by_read.values() {
        let chain = paf_non_overlapping_chain(alns.clone(), max_gap);
        for pair in chain.windows(2) {
            let left = &pair[0];
            let right = &pair[1];
            let gap = right.qstart as isize - left.qend as isize;
            if gap.abs() > max_gap {
                continue;
            }
            let left_locus = mito_stable_paf_exit_locus(left);
            let right_locus = mito_stable_paf_entry_locus(right);
            if mito_stable_locus_is_split_candidate(config, &segment_lengths, &left_locus) {
                events
                    .entry(left_locus.segment.clone())
                    .or_default()
                    .push(left_locus.coord);
            }
            if mito_stable_locus_is_split_candidate(config, &segment_lengths, &right_locus) {
                events
                    .entry(right_locus.segment.clone())
                    .or_default()
                    .push(right_locus.coord);
            }
        }
    }

    let min_support = mito_stable_internal_split_min_support(config);
    let cluster_slop = config.skeleton_end_slop.max(500);
    let min_piece_len = mito_stable_min_split_piece_len(config);
    let mut split_points = HashMap::new();
    for (segment, mut coords) in events {
        coords.sort_unstable();
        let Some(&segment_len) = segment_lengths.get(&segment) else {
            continue;
        };
        let mut clusters: Vec<Vec<usize>> = Vec::new();
        for coord in coords {
            if let Some(last_cluster) = clusters.last_mut() {
                let center =
                    last_cluster.iter().copied().sum::<usize>() / last_cluster.len().max(1);
                if coord.abs_diff(center) <= cluster_slop {
                    last_cluster.push(coord);
                    continue;
                }
            }
            clusters.push(vec![coord]);
        }

        let mut points = Vec::new();
        for cluster in clusters {
            let support = cluster.len() as u32;
            if support < min_support {
                continue;
            }
            let position = cluster.iter().copied().sum::<usize>() / cluster.len().max(1);
            if position < min_piece_len || segment_len.saturating_sub(position) < min_piece_len {
                continue;
            }
            points.push(MitoStableSplitPoint { position, support });
        }
        points.sort_by_key(|point| point.position);
        if !points.is_empty() {
            split_points.insert(segment, points);
        }
    }
    split_points
}

fn mito_stable_paf_exit_locus(aln: &PafAln) -> MitoStablePafLocus {
    let coord = if aln.strand == '+' {
        aln.tend
    } else {
        aln.tstart
    };
    MitoStablePafLocus {
        segment: aln.tname.clone(),
        coord,
        tlen: aln.tlen,
    }
}

fn mito_stable_paf_entry_locus(aln: &PafAln) -> MitoStablePafLocus {
    let coord = if aln.strand == '+' {
        aln.tstart
    } else {
        aln.tend
    };
    MitoStablePafLocus {
        segment: aln.tname.clone(),
        coord,
        tlen: aln.tlen,
    }
}

fn mito_stable_locus_is_split_candidate(
    config: &Config,
    segment_lengths: &HashMap<String, usize>,
    locus: &MitoStablePafLocus,
) -> bool {
    let segment_len = segment_lengths
        .get(&locus.segment)
        .copied()
        .unwrap_or(locus.tlen);
    locus.coord > config.skeleton_end_slop
        && segment_len.saturating_sub(locus.coord) > config.skeleton_end_slop
}

fn mito_stable_internal_split_min_support(config: &Config) -> u32 {
    config.skeleton_min_link_support.max(30)
}

fn mito_stable_min_split_piece_len(config: &Config) -> usize {
    config.skeleton_end_slop.min(1_000).max(500)
}

fn apply_mito_stable_internal_splits(
    config: &Config,
    segments: Vec<SkeletonSegment>,
    links: HashMap<SkeletonLinkKey, SkeletonLinkSupport>,
    split_points: HashMap<String, Vec<MitoStableSplitPoint>>,
) -> (
    Vec<SkeletonSegment>,
    HashMap<SkeletonLinkKey, SkeletonLinkSupport>,
    Vec<MitoStableSplitReportRow>,
) {
    let mut layouts: HashMap<String, MitoStableSplitLayout> = HashMap::new();
    let mut refined_segments = Vec::new();
    let mut report = Vec::new();

    for segment in segments {
        let mut breakpoints: Vec<_> = split_points
            .get(&segment.name)
            .map(|points| {
                points
                    .iter()
                    .map(|point| (point.position, point.support))
                    .collect()
            })
            .unwrap_or_default();
        breakpoints.sort_unstable_by_key(|(position, _)| *position);
        breakpoints.dedup_by_key(|(position, _)| *position);
        if breakpoints.is_empty() {
            layouts.insert(
                segment.name.clone(),
                MitoStableSplitLayout {
                    pieces: vec![segment.name.clone()],
                },
            );
            refined_segments.push(segment);
            continue;
        }

        let mut cuts = Vec::with_capacity(breakpoints.len() + 2);
        cuts.push(0usize);
        for (position, support) in breakpoints {
            if position == 0 || position >= segment.sequence.len() {
                continue;
            }
            cuts.push(position);
            report.push(MitoStableSplitReportRow {
                segment: segment.name.clone(),
                position,
                support,
            });
        }
        cuts.push(segment.sequence.len());
        cuts.sort_unstable();
        cuts.dedup();

        if cuts.len() <= 2 {
            layouts.insert(
                segment.name.clone(),
                MitoStableSplitLayout {
                    pieces: vec![segment.name.clone()],
                },
            );
            refined_segments.push(segment);
            continue;
        }

        let mut piece_names = Vec::new();
        for (index, window) in cuts.windows(2).enumerate() {
            let start = window[0];
            let end = window[1];
            if start >= end {
                continue;
            }
            let name = format!("{}_s{}", segment.name, index);
            piece_names.push(name.clone());
            refined_segments.push(SkeletonSegment {
                name,
                sequence: segment.sequence[start..end].to_string(),
            });
        }
        layouts.insert(
            segment.name,
            MitoStableSplitLayout {
                pieces: piece_names,
            },
        );
    }

    let mut refined_links: HashMap<SkeletonLinkKey, SkeletonLinkSupport> = HashMap::new();
    for (key, support) in links {
        let from = mito_stable_split_exit(&layouts, &key.from, key.from_orient);
        let to = mito_stable_split_entry(&layouts, &key.to, key.to_orient);
        let Some((from_name, from_orient)) = from else {
            continue;
        };
        let Some((to_name, to_orient)) = to else {
            continue;
        };
        let refined_key = SkeletonLinkKey {
            from: from_name,
            from_orient,
            to: to_name,
            to_orient,
        };
        merge_skeleton_link_support(&mut refined_links, refined_key, support);
    }

    let internal_support = config
        .min_link_support
        .max(config.skeleton_min_link_support)
        .max(1);
    for layout in layouts.values() {
        for pair in layout.pieces.windows(2) {
            let key = SkeletonLinkKey {
                from: pair[0].clone(),
                from_orient: '+',
                to: pair[1].clone(),
                to_orient: '+',
            };
            insert_mito_stable_internal_link(config, &mut refined_links, key, internal_support);
        }
    }

    report.sort_by(|left, right| {
        left.segment
            .cmp(&right.segment)
            .then_with(|| left.position.cmp(&right.position))
    });
    (refined_segments, refined_links, report)
}

fn mito_stable_split_exit(
    layouts: &HashMap<String, MitoStableSplitLayout>,
    segment: &str,
    orient: char,
) -> Option<(String, char)> {
    let layout = layouts.get(segment)?;
    if orient == '+' {
        layout.pieces.last().cloned().map(|piece| (piece, '+'))
    } else {
        layout.pieces.first().cloned().map(|piece| (piece, '-'))
    }
}

fn mito_stable_split_entry(
    layouts: &HashMap<String, MitoStableSplitLayout>,
    segment: &str,
    orient: char,
) -> Option<(String, char)> {
    let layout = layouts.get(segment)?;
    if orient == '+' {
        layout.pieces.first().cloned().map(|piece| (piece, '+'))
    } else {
        layout.pieces.last().cloned().map(|piece| (piece, '-'))
    }
}

fn merge_skeleton_link_support(
    links: &mut HashMap<SkeletonLinkKey, SkeletonLinkSupport>,
    key: SkeletonLinkKey,
    support: SkeletonLinkSupport,
) {
    let entry = links.entry(key).or_default();
    entry.skeleton_support = entry.skeleton_support.max(support.skeleton_support);
    entry.paf_support = entry.paf_support.max(support.paf_support);
    entry.rescue_support = entry.rescue_support.max(support.rescue_support);
    entry.out_ratio = entry.out_ratio.max(support.out_ratio);
    entry.in_ratio = entry.in_ratio.max(support.in_ratio);
}

fn insert_mito_stable_internal_link(
    config: &Config,
    links: &mut HashMap<SkeletonLinkKey, SkeletonLinkSupport>,
    key: SkeletonLinkKey,
    support: u32,
) {
    links.entry(key.clone()).or_default().skeleton_support = support;
    if config.bidirectional_links {
        let rc = reverse_skeleton_link_key(&key);
        links.entry(rc).or_default().skeleton_support = support;
    }
}

fn mito_stable_base_links_and_paf_candidates(
    mut base_links: HashMap<SkeletonLinkKey, SkeletonLinkSupport>,
    paf_links: HashMap<SkeletonLinkKey, u32>,
) -> (
    HashMap<SkeletonLinkKey, SkeletonLinkSupport>,
    Vec<MitoCompactBridgeCandidate>,
) {
    let mut candidates = Vec::new();
    for (key, support) in paf_links {
        if let Some(entry) = base_links.get_mut(&key) {
            entry.paf_support += support;
        } else {
            candidates.push(MitoCompactBridgeCandidate {
                key,
                support,
                source: MitoCompactBridgeSource::LocalPaf,
            });
        }
    }
    (base_links, candidates)
}

fn dedupe_mito_bridge_candidates(
    candidates: Vec<MitoCompactBridgeCandidate>,
) -> Vec<MitoCompactBridgeCandidate> {
    let mut by_key: HashMap<SkeletonLinkKey, MitoCompactBridgeCandidate> = HashMap::new();
    for candidate in candidates {
        let canonical = canonical_skeleton_link_key(&candidate.key);
        by_key
            .entry(canonical)
            .and_modify(|old| {
                let old_score = (old.source.priority(), old.support);
                let new_score = (candidate.source.priority(), candidate.support);
                if new_score > old_score {
                    *old = candidate.clone();
                }
            })
            .or_insert(candidate);
    }
    let mut candidates: Vec<_> = by_key.into_values().collect();
    candidates.sort_by(|left, right| {
        right
            .source
            .priority()
            .cmp(&left.source.priority())
            .then_with(|| right.support.cmp(&left.support))
            .then_with(|| left.key.cmp(&right.key))
    });
    candidates
}

fn canonical_skeleton_link_key(key: &SkeletonLinkKey) -> SkeletonLinkKey {
    let rc = reverse_skeleton_link_key(key);
    if *key <= rc {
        key.clone()
    } else {
        rc
    }
}

fn select_mito_stable_bridges(
    config: &Config,
    segments: &[SkeletonSegment],
    base_links: &HashMap<SkeletonLinkKey, SkeletonLinkSupport>,
    candidates: &[MitoCompactBridgeCandidate],
    respect_node_constraints: bool,
) -> Vec<MitoCompactBridgeCandidate> {
    if candidates.len() <= 20 {
        return select_mito_stable_bridges_exhaustive(
            config,
            segments,
            base_links,
            candidates,
            respect_node_constraints,
        );
    }
    select_mito_stable_bridges_beam(
        config,
        segments,
        base_links,
        candidates,
        respect_node_constraints,
    )
}

fn mito_stable_scoped_candidates(
    config: &Config,
    candidates: Vec<MitoCompactBridgeCandidate>,
) -> Vec<MitoCompactBridgeCandidate> {
    match config.mito_stable_selection_mode {
        MitoStableSelectionMode::Auto => candidates,
        #[cfg(test)]
        MitoStableSelectionMode::ForbidSelected => {
            scoped_mito_stable_candidates(config, candidates)
        }
        MitoStableSelectionMode::AllowSelected => scoped_mito_stable_candidates(config, candidates),
    }
}

fn scoped_mito_stable_candidates(
    config: &Config,
    candidates: Vec<MitoCompactBridgeCandidate>,
) -> Vec<MitoCompactBridgeCandidate> {
    candidates
        .into_iter()
        .filter(|candidate| {
            config
                .mito_stable_selected_nodes
                .contains(&candidate.key.from)
                || config
                    .mito_stable_selected_nodes
                    .contains(&candidate.key.to)
        })
        .collect()
}

fn append_mito_stable_selected_bridges(
    selected: &mut Vec<MitoCompactBridgeCandidate>,
    extra: Vec<MitoCompactBridgeCandidate>,
) {
    let mut seen_physical: HashSet<_> = selected
        .iter()
        .map(|candidate| canonical_mito_stable_physical_link(&candidate.key))
        .collect();
    for candidate in extra {
        if seen_physical.insert(canonical_mito_stable_physical_link(&candidate.key)) {
            selected.push(candidate);
        }
    }
}

fn select_mito_stable_bridges_exhaustive(
    config: &Config,
    segments: &[SkeletonSegment],
    base_links: &HashMap<SkeletonLinkKey, SkeletonLinkSupport>,
    candidates: &[MitoCompactBridgeCandidate],
    respect_node_constraints: bool,
) -> Vec<MitoCompactBridgeCandidate> {
    let mut best_mask = 0usize;
    let mut best_score = None;
    let mask_count = 1usize << candidates.len();
    for mask in 0..mask_count {
        let mut trial_links = base_links.clone();
        let mut selected_count = 0usize;
        let mut support_sum = 0u64;
        for (index, candidate) in candidates.iter().enumerate() {
            if (mask & (1usize << index)) == 0 {
                continue;
            }
            insert_mito_compact_bridge(config, &mut trial_links, candidate);
            selected_count += 1;
            support_sum += candidate.support as u64;
        }
        refresh_skeleton_link_ratios(config, &mut trial_links);
        let score = mito_stable_bridge_selection_score(
            config,
            segments,
            &trial_links,
            selected_count,
            support_sum,
            respect_node_constraints,
        );
        if best_score.as_ref().is_none_or(|best| score < *best) {
            best_score = Some(score);
            best_mask = mask;
        }
    }

    candidates
        .iter()
        .enumerate()
        .filter(|(index, _)| (best_mask & (1usize << index)) != 0)
        .map(|(_, candidate)| candidate.clone())
        .collect()
}

#[derive(Clone)]
struct MitoStableBeamState {
    links: HashMap<SkeletonLinkKey, SkeletonLinkSupport>,
    selected: Vec<MitoCompactBridgeCandidate>,
    support_sum: u64,
}

fn select_mito_stable_bridges_beam(
    config: &Config,
    segments: &[SkeletonSegment],
    base_links: &HashMap<SkeletonLinkKey, SkeletonLinkSupport>,
    candidates: &[MitoCompactBridgeCandidate],
    respect_node_constraints: bool,
) -> Vec<MitoCompactBridgeCandidate> {
    const BEAM_WIDTH: usize = 256;
    let mut states = vec![MitoStableBeamState {
        links: base_links.clone(),
        selected: Vec::new(),
        support_sum: 0,
    }];

    for candidate in candidates {
        let mut next = Vec::with_capacity(states.len() * 2);
        for state in states {
            next.push(state.clone());
            let mut with_candidate = state;
            insert_mito_compact_bridge(config, &mut with_candidate.links, candidate);
            with_candidate.selected.push(candidate.clone());
            with_candidate.support_sum += candidate.support as u64;
            next.push(with_candidate);
        }
        for state in &mut next {
            refresh_skeleton_link_ratios(config, &mut state.links);
        }
        next.sort_by_key(|state| {
            mito_stable_bridge_selection_score(
                config,
                segments,
                &state.links,
                state.selected.len(),
                state.support_sum,
                respect_node_constraints,
            )
        });
        next.truncate(BEAM_WIDTH);
        states = next;
    }

    states
        .into_iter()
        .min_by_key(|state| {
            mito_stable_bridge_selection_score(
                config,
                segments,
                &state.links,
                state.selected.len(),
                state.support_sum,
                respect_node_constraints,
            )
        })
        .map(|state| state.selected)
        .unwrap_or_default()
}

fn mito_stable_bridge_selection_score(
    config: &Config,
    segments: &[SkeletonSegment],
    links: &HashMap<SkeletonLinkKey, SkeletonLinkSupport>,
    selected_links: usize,
    support_sum: u64,
    respect_node_constraints: bool,
) -> (
    usize,
    usize,
    usize,
    usize,
    usize,
    usize,
    usize,
    usize,
    std::cmp::Reverse<u64>,
) {
    let constraint_violations = if respect_node_constraints {
        mito_stable_node_shape_constraint_violations(config, segments, links)
    } else {
        0
    };
    let topology = mito_stable_topology_stats(config, segments, links);
    let (
        components,
        endpoint_overload,
        open_ends,
        cycle_rank,
        unique_links,
        branch_ends,
        selected_links,
        support_sum,
    ) = topology.score_key(selected_links, support_sum);
    (
        constraint_violations,
        components,
        endpoint_overload,
        open_ends,
        cycle_rank,
        unique_links,
        branch_ends,
        selected_links,
        support_sum,
    )
}

fn mito_stable_topology_stats(
    config: &Config,
    segments: &[SkeletonSegment],
    links: &HashMap<SkeletonLinkKey, SkeletonLinkSupport>,
) -> MitoStableTopologyStats {
    let mut segment_names = HashSet::new();
    let mut endpoint_degree: HashMap<(String, char), usize> = HashMap::new();
    let mut adj: HashMap<String, HashSet<String>> = HashMap::new();
    let mut unique_links = HashSet::new();

    for segment in segments {
        segment_names.insert(segment.name.clone());
        adj.entry(segment.name.clone()).or_default();
    }

    for (key, support) in links {
        if !keep_skeleton_link(config, support) {
            continue;
        }
        if !segment_names.contains(&key.from) || !segment_names.contains(&key.to) {
            continue;
        }
        let left = skeleton_physical_from_end(&key.from, key.from_orient);
        let right = skeleton_physical_to_end(&key.to, key.to_orient);
        let edge = if left <= right {
            (left.clone(), right.clone())
        } else {
            (right.clone(), left.clone())
        };
        if unique_links.insert(edge) {
            *endpoint_degree.entry(left).or_insert(0) += 1;
            *endpoint_degree.entry(right).or_insert(0) += 1;
            adj.entry(key.from.clone())
                .or_default()
                .insert(key.to.clone());
            adj.entry(key.to.clone())
                .or_default()
                .insert(key.from.clone());
        }
    }

    let mut seen = HashSet::new();
    let mut components = 0usize;
    for segment in segments {
        if !seen.insert(segment.name.clone()) {
            continue;
        }
        components += 1;
        let mut stack = vec![segment.name.clone()];
        while let Some(node) = stack.pop() {
            if let Some(neighbors) = adj.get(&node) {
                for neighbor in neighbors {
                    if seen.insert(neighbor.clone()) {
                        stack.push(neighbor.clone());
                    }
                }
            }
        }
    }

    let mut open_ends = 0usize;
    let mut endpoint_overload = 0usize;
    let mut branch_ends = 0usize;
    for segment in segments {
        for end in ['L', 'R'] {
            let degree = endpoint_degree
                .get(&(segment.name.clone(), end))
                .copied()
                .unwrap_or(0);
            if degree == 0 {
                open_ends += 1;
            }
            if degree > 1 {
                branch_ends += 1;
            }
            endpoint_overload += degree.saturating_sub(2);
        }
    }
    let cycle_rank = unique_links
        .len()
        .saturating_add(components)
        .saturating_sub(segments.len());

    MitoStableTopologyStats {
        components,
        open_ends,
        endpoint_overload,
        unique_links: unique_links.len(),
        branch_ends,
        cycle_rank,
    }
}

fn mito_stable_physical_side_degrees(
    config: &Config,
    segments: &[SkeletonSegment],
    links: &HashMap<SkeletonLinkKey, SkeletonLinkSupport>,
) -> HashMap<String, (usize, usize)> {
    let segment_names: HashSet<_> = segments
        .iter()
        .map(|segment| segment.name.clone())
        .collect();
    let mut degrees: HashMap<String, (usize, usize)> = segments
        .iter()
        .map(|segment| (segment.name.clone(), (0, 0)))
        .collect();
    let mut seen = HashSet::new();
    for (key, support) in links {
        if !keep_skeleton_link(config, support) {
            continue;
        }
        if !segment_names.contains(&key.from) || !segment_names.contains(&key.to) {
            continue;
        }
        let left = skeleton_physical_from_end(&key.from, key.from_orient);
        let right = skeleton_physical_to_end(&key.to, key.to_orient);
        let edge = if left <= right {
            (left.clone(), right.clone())
        } else {
            (right.clone(), left.clone())
        };
        if !seen.insert(edge) {
            continue;
        }
        increment_mito_stable_side_degree(&mut degrees, &left);
        increment_mito_stable_side_degree(&mut degrees, &right);
    }
    degrees
}

fn increment_mito_stable_side_degree(
    degrees: &mut HashMap<String, (usize, usize)>,
    endpoint: &(String, char),
) {
    let entry = degrees.entry(endpoint.0.clone()).or_insert((0, 0));
    if endpoint.1 == 'L' {
        entry.0 += 1;
    } else {
        entry.1 += 1;
    }
}

fn mito_stable_three_way_nodes(
    config: &Config,
    segments: &[SkeletonSegment],
    links: &HashMap<SkeletonLinkKey, SkeletonLinkSupport>,
) -> HashSet<String> {
    mito_stable_physical_side_degrees(config, segments, links)
        .into_iter()
        .filter_map(|(segment, (left_degree, right_degree))| {
            if (left_degree == 1 && right_degree == 2) || (left_degree == 2 && right_degree == 1) {
                Some(segment)
            } else {
                None
            }
        })
        .collect()
}

fn mito_stable_introduces_new_three_way_nodes(
    config: &Config,
    before_segments: &[SkeletonSegment],
    before_links: &HashMap<SkeletonLinkKey, SkeletonLinkSupport>,
    after_segments: &[SkeletonSegment],
    after_links: &HashMap<SkeletonLinkKey, SkeletonLinkSupport>,
) -> bool {
    let before = mito_stable_three_way_nodes(config, before_segments, before_links);
    let after = mito_stable_three_way_nodes(config, after_segments, after_links);
    after.iter().any(|segment| !before.contains(segment))
}

fn mito_stable_has_explicit_three_way_constraints(config: &Config) -> bool {
    !config.mito_stable_allow_three_way.is_empty()
        || !config.mito_stable_forbid_three_way.is_empty()
}

fn mito_stable_node_constraints_ok(
    config: &Config,
    segments: &[SkeletonSegment],
    links: &HashMap<SkeletonLinkKey, SkeletonLinkSupport>,
) -> bool {
    mito_stable_three_way_constraint_violations(config, segments, links) == 0
}

fn mito_stable_three_way_constraint_violations(
    config: &Config,
    segments: &[SkeletonSegment],
    links: &HashMap<SkeletonLinkKey, SkeletonLinkSupport>,
) -> usize {
    if !mito_stable_has_explicit_three_way_constraints(config) {
        return 0;
    }
    let three_way = mito_stable_three_way_nodes(config, segments, links);
    let mut violations = HashSet::new();
    for segment in three_way {
        if config.mito_stable_forbid_three_way.contains(&segment)
            || (!config.mito_stable_allow_three_way.is_empty()
                && !config.mito_stable_allow_three_way.contains(&segment))
        {
            violations.insert(segment);
        }
    }
    violations.len()
}

fn mito_stable_node_shape_constraint_violations(
    config: &Config,
    segments: &[SkeletonSegment],
    links: &HashMap<SkeletonLinkKey, SkeletonLinkSupport>,
) -> usize {
    if !mito_stable_has_explicit_three_way_constraints(config) {
        return 0;
    }
    let mut violations = HashSet::new();
    for (segment, (left_degree, right_degree)) in
        mito_stable_physical_side_degrees(config, segments, links)
    {
        match mito_stable_node_degree_class(left_degree, right_degree) {
            "1-1" | "2-2" => {}
            "three_way" => {
                if config.mito_stable_forbid_three_way.contains(&segment)
                    || (!config.mito_stable_allow_three_way.is_empty()
                        && !config.mito_stable_allow_three_way.contains(&segment))
                {
                    violations.insert(segment);
                }
            }
            _ => {
                violations.insert(segment);
            }
        }
    }
    violations.len()
}

fn mito_stable_node_constraints_ok_for_transition(
    config: &Config,
    before_segments: &[SkeletonSegment],
    before_links: &HashMap<SkeletonLinkKey, SkeletonLinkSupport>,
    after_segments: &[SkeletonSegment],
    after_links: &HashMap<SkeletonLinkKey, SkeletonLinkSupport>,
) -> bool {
    if mito_stable_has_explicit_three_way_constraints(config) {
        return mito_stable_node_constraints_ok(config, after_segments, after_links);
    }
    !mito_stable_introduces_new_three_way_nodes(
        config,
        before_segments,
        before_links,
        after_segments,
        after_links,
    )
}

fn mito_stable_node_degree_class(left_degree: usize, right_degree: usize) -> &'static str {
    match (left_degree, right_degree) {
        (1, 1) => "1-1",
        (2, 2) => "2-2",
        (1, 2) | (2, 1) => "three_way",
        _ => "other",
    }
}

fn expand_mito_stable_repeat_shortcuts(
    config: &Config,
    segments: &mut Vec<SkeletonSegment>,
    depth: &mut HashMap<String, (usize, usize)>,
    links: &mut HashMap<SkeletonLinkKey, SkeletonLinkSupport>,
    candidates: &[MitoCompactBridgeCandidate],
) -> Vec<MitoStableRepeatExpansion> {
    const MAX_EXPANSIONS: usize = 16;
    let mut expansions = Vec::new();

    while expansions.len() < MAX_EXPANSIONS {
        let Some(expansion) =
            best_mito_stable_repeat_expansion(config, segments, links, candidates)
        else {
            break;
        };

        let Some(template_segment) = segments
            .iter()
            .find(|segment| segment.name == expansion.template)
            .cloned()
        else {
            break;
        };

        let mut left_support =
            skeleton_link_support_for_key(links, &expansion.left_link).unwrap_or_default();
        remove_mito_stable_physical_link_pair(links, &expansion.removed_shortcut);
        remove_mito_stable_physical_link_pair(links, &expansion.left_link);

        if left_support.skeleton_support == 0
            && left_support.paf_support == 0
            && left_support.rescue_support == 0
        {
            left_support.skeleton_support = expansion.support;
        }
        let mut right_support =
            mito_stable_available_link_support(config, links, candidates, &expansion.right_link)
                .unwrap_or_default();
        if right_support.skeleton_support == 0
            && right_support.paf_support == 0
            && right_support.rescue_support == 0
        {
            right_support.paf_support = expansion.support;
        }

        let clone_name = expansion.clone.clone();
        segments.push(SkeletonSegment {
            name: clone_name.clone(),
            sequence: template_segment.sequence,
        });
        if let Some(template_depth) = depth.get(&expansion.template).copied() {
            depth.insert(clone_name.clone(), template_depth);
        }

        let mut left_clone = expansion.left_link.clone();
        left_clone.to = clone_name.clone();
        let mut right_clone = expansion.right_link.clone();
        right_clone.from = clone_name;
        insert_mito_stable_supported_link_pair(config, links, left_clone, left_support);
        insert_mito_stable_supported_link_pair(config, links, right_clone, right_support);

        expansions.push(expansion);
        refresh_skeleton_link_ratios(config, links);
    }

    expansions
}

fn best_mito_stable_repeat_expansion(
    config: &Config,
    segments: &[SkeletonSegment],
    links: &HashMap<SkeletonLinkKey, SkeletonLinkSupport>,
    candidates: &[MitoCompactBridgeCandidate],
) -> Option<MitoStableRepeatExpansion> {
    const MIN_REPEAT_OVERLAP: usize = 100;
    const MAX_SHORTCUT_OVERLAP: usize = 20;

    let sequence_by_segment: HashMap<String, String> = segments
        .iter()
        .map(|segment| (segment.name.clone(), segment.sequence.clone()))
        .collect();
    let available = mito_stable_available_links(config, links, candidates);
    let mut seen_physical = HashSet::new();
    let mut best: Option<(
        std::cmp::Reverse<usize>,
        std::cmp::Reverse<u32>,
        u32,
        MitoStableRepeatExpansion,
    )> = None;
    let mut clone_index = 0usize;

    for (shortcut, shortcut_support) in links {
        if !keep_skeleton_link(config, shortcut_support) {
            continue;
        }
        if shortcut_support.skeleton_support > 0 || shortcut_support.rescue_support > 0 {
            continue;
        }
        if !seen_physical.insert(canonical_mito_stable_physical_link(shortcut)) {
            continue;
        }

        let shortcut_forms = [shortcut.clone(), reverse_skeleton_link_key(shortcut)];
        for shortcut_key in shortcut_forms {
            let shortcut_overlap =
                skeleton_link_exact_overlap(&sequence_by_segment, &shortcut_key, 2_000);
            if shortcut_overlap > MAX_SHORTCUT_OVERLAP {
                continue;
            }

            for (left_key, left_support) in &available {
                if left_key.from != shortcut_key.from
                    || left_key.from_orient != shortcut_key.from_orient
                    || left_key.to == shortcut_key.to
                {
                    continue;
                }
                if !skeleton_link_pair_is_kept(config, links, left_key) {
                    continue;
                }
                let left_overlap =
                    skeleton_link_exact_overlap(&sequence_by_segment, left_key, 2_000);
                if left_overlap < MIN_REPEAT_OVERLAP {
                    continue;
                }

                let right_key = SkeletonLinkKey {
                    from: left_key.to.clone(),
                    from_orient: left_key.to_orient,
                    to: shortcut_key.to.clone(),
                    to_orient: shortcut_key.to_orient,
                };
                let Some(right_support) = available.get(&right_key) else {
                    continue;
                };
                let right_overlap =
                    skeleton_link_exact_overlap(&sequence_by_segment, &right_key, 2_000);
                if right_overlap < MIN_REPEAT_OVERLAP {
                    continue;
                }
                if right_support.paf_support.max(right_support.rescue_support)
                    < config.skeleton_min_link_support
                    && right_support.skeleton_support < config.min_link_support
                {
                    continue;
                }

                let support = left_support
                    .skeleton_support
                    .max(left_support.paf_support)
                    .max(left_support.rescue_support)
                    .min(
                        right_support
                            .skeleton_support
                            .max(right_support.paf_support)
                            .max(right_support.rescue_support),
                    );
                let clone =
                    unique_mito_stable_clone_name(&sequence_by_segment, &left_key.to, clone_index);
                clone_index += 1;
                let expansion = MitoStableRepeatExpansion {
                    template: left_key.to.clone(),
                    clone,
                    removed_shortcut: shortcut_key.clone(),
                    left_link: left_key.clone(),
                    right_link: right_key,
                    shortcut_overlap,
                    left_overlap,
                    right_overlap,
                    support,
                };
                if !mito_stable_repeat_expansion_topology_ok(
                    config, segments, links, &available, &expansion,
                ) {
                    continue;
                }
                let score = (
                    std::cmp::Reverse(left_overlap.min(right_overlap)),
                    std::cmp::Reverse(support),
                    shortcut_support.paf_support,
                    expansion,
                );
                if best.as_ref().is_none_or(|old| score < *old) {
                    best = Some(score);
                }
            }
        }
    }

    best.map(|(_, _, _, expansion)| expansion)
}

fn mito_stable_repeat_expansion_topology_ok(
    config: &Config,
    segments: &[SkeletonSegment],
    links: &HashMap<SkeletonLinkKey, SkeletonLinkSupport>,
    available: &HashMap<SkeletonLinkKey, SkeletonLinkSupport>,
    expansion: &MitoStableRepeatExpansion,
) -> bool {
    let before = mito_stable_topology_stats(config, segments, links);
    let Some(template_segment) = segments
        .iter()
        .find(|segment| segment.name == expansion.template)
    else {
        return false;
    };
    let mut trial_segments = segments.to_vec();
    trial_segments.push(SkeletonSegment {
        name: expansion.clone.clone(),
        sequence: template_segment.sequence.clone(),
    });

    let mut trial_links = links.clone();
    let left_support =
        skeleton_link_support_for_key(&trial_links, &expansion.left_link).unwrap_or_default();
    let right_support = available
        .get(&expansion.right_link)
        .cloned()
        .unwrap_or_default();
    remove_mito_stable_physical_link_pair(&mut trial_links, &expansion.removed_shortcut);
    remove_mito_stable_physical_link_pair(&mut trial_links, &expansion.left_link);

    let mut left_clone = expansion.left_link.clone();
    left_clone.to = expansion.clone.clone();
    let mut right_clone = expansion.right_link.clone();
    right_clone.from = expansion.clone.clone();
    insert_mito_stable_supported_link_pair(config, &mut trial_links, left_clone, left_support);
    insert_mito_stable_supported_link_pair(config, &mut trial_links, right_clone, right_support);
    refresh_skeleton_link_ratios(config, &mut trial_links);

    let after = mito_stable_topology_stats(config, &trial_segments, &trial_links);
    after.components <= before.components
        && after.open_ends <= before.open_ends
        && after.endpoint_overload <= before.endpoint_overload
        && after.unique_links <= before.unique_links
        && (after.branch_ends < before.branch_ends || after.cycle_rank < before.cycle_rank)
        && mito_stable_node_constraints_ok_for_transition(
            config,
            segments,
            links,
            &trial_segments,
            &trial_links,
        )
}

fn unique_mito_stable_clone_name(
    sequence_by_segment: &HashMap<String, String>,
    template: &str,
    mut index: usize,
) -> String {
    loop {
        let name = format!("{}_copy{}", template, index);
        if !sequence_by_segment.contains_key(&name) {
            return name;
        }
        index += 1;
    }
}

fn mito_stable_available_links(
    config: &Config,
    links: &HashMap<SkeletonLinkKey, SkeletonLinkSupport>,
    candidates: &[MitoCompactBridgeCandidate],
) -> HashMap<SkeletonLinkKey, SkeletonLinkSupport> {
    let mut available = links.clone();
    for candidate in candidates {
        let mut support = SkeletonLinkSupport::default();
        match candidate.source {
            MitoCompactBridgeSource::LocalPaf | MitoCompactBridgeSource::CopyChoice => {
                support.paf_support = candidate.support;
            }
            MitoCompactBridgeSource::Rescue => support.rescue_support = candidate.support,
        }
        merge_skeleton_link_support(&mut available, candidate.key.clone(), support.clone());
        if config.bidirectional_links {
            merge_skeleton_link_support(
                &mut available,
                reverse_skeleton_link_key(&candidate.key),
                support,
            );
        }
    }
    available
}

fn mito_stable_available_link_support(
    config: &Config,
    links: &HashMap<SkeletonLinkKey, SkeletonLinkSupport>,
    candidates: &[MitoCompactBridgeCandidate],
    key: &SkeletonLinkKey,
) -> Option<SkeletonLinkSupport> {
    mito_stable_available_links(config, links, candidates)
        .get(key)
        .cloned()
}

fn skeleton_link_support_for_key(
    links: &HashMap<SkeletonLinkKey, SkeletonLinkSupport>,
    key: &SkeletonLinkKey,
) -> Option<SkeletonLinkSupport> {
    links
        .get(key)
        .cloned()
        .or_else(|| links.get(&reverse_skeleton_link_key(key)).cloned())
}

fn skeleton_link_pair_is_kept(
    config: &Config,
    links: &HashMap<SkeletonLinkKey, SkeletonLinkSupport>,
    key: &SkeletonLinkKey,
) -> bool {
    skeleton_link_support_for_key(links, key)
        .is_some_and(|support| keep_skeleton_link(config, &support))
}

fn insert_mito_stable_supported_link_pair(
    config: &Config,
    links: &mut HashMap<SkeletonLinkKey, SkeletonLinkSupport>,
    key: SkeletonLinkKey,
    support: SkeletonLinkSupport,
) {
    merge_skeleton_link_support(links, key.clone(), support.clone());
    if config.bidirectional_links {
        merge_skeleton_link_support(links, reverse_skeleton_link_key(&key), support);
    }
}

fn skeleton_link_exact_overlap(
    sequence_by_segment: &HashMap<String, String>,
    key: &SkeletonLinkKey,
    max_overlap: usize,
) -> usize {
    let Some(from_sequence) = sequence_by_segment.get(&key.from) else {
        return 0;
    };
    let Some(to_sequence) = sequence_by_segment.get(&key.to) else {
        return 0;
    };
    let from = oriented_skeleton_sequence(from_sequence, key.from_orient);
    let to = oriented_skeleton_sequence(to_sequence, key.to_orient);
    let max_len = max_overlap.min(from.len()).min(to.len());
    let mut best = 0usize;
    for len in 1..=max_len {
        if from[from.len() - len..] == to[..len] {
            best = len;
        }
    }
    best
}

fn oriented_skeleton_sequence(sequence: &str, orient: char) -> String {
    if orient == '+' {
        sequence.to_string()
    } else {
        revcomp_string(sequence)
    }
}

fn prune_mito_stable_redundant_links(
    config: &Config,
    segments: &[SkeletonSegment],
    links: &mut HashMap<SkeletonLinkKey, SkeletonLinkSupport>,
) -> Vec<SkeletonLinkKey> {
    let mut pruned = Vec::new();
    loop {
        let before = mito_stable_topology_stats(config, segments, links);
        let mut seen = HashSet::new();
        let mut best: Option<(MitoStableTopologyStats, bool, u32, SkeletonLinkKey)> = None;

        let keys: Vec<_> = links
            .iter()
            .filter(|(_, support)| keep_skeleton_link(config, support))
            .filter_map(|(key, support)| {
                let physical = canonical_mito_stable_physical_link(key);
                if !seen.insert(physical) {
                    return None;
                }
                let support_score = support
                    .skeleton_support
                    .max(support.paf_support)
                    .max(support.rescue_support);
                Some((key.clone(), support.skeleton_support > 0, support_score))
            })
            .collect();

        for (key, has_skeleton_support, support_score) in keys {
            let mut trial = links.clone();
            remove_mito_stable_physical_link_pair(&mut trial, &key);
            refresh_skeleton_link_ratios(config, &mut trial);
            let after = mito_stable_topology_stats(config, segments, &trial);
            if !mito_stable_redundant_prune_improves(before, after) {
                continue;
            }
            if !mito_stable_node_constraints_ok_for_transition(
                config, segments, links, segments, &trial,
            ) {
                continue;
            }
            let score = (after, has_skeleton_support, support_score, key.clone());
            if best.as_ref().is_none_or(|old| score < *old) {
                best = Some(score);
            }
        }

        let Some((_, _, _, key)) = best else {
            break;
        };
        remove_mito_stable_physical_link_pair(links, &key);
        pruned.push(key);
    }
    pruned
}

fn mito_stable_redundant_prune_improves(
    before: MitoStableTopologyStats,
    after: MitoStableTopologyStats,
) -> bool {
    after.components == before.components
        && after.open_ends <= before.open_ends
        && after.endpoint_overload <= before.endpoint_overload
        && after.unique_links < before.unique_links
        && after.branch_ends < before.branch_ends
        && after.cycle_rank < before.cycle_rank
}

fn canonical_mito_stable_physical_link(key: &SkeletonLinkKey) -> ((String, char), (String, char)) {
    let left = skeleton_physical_from_end(&key.from, key.from_orient);
    let right = skeleton_physical_to_end(&key.to, key.to_orient);
    canonical_physical_end_pair(left, right)
}

fn canonical_physical_end_pair(
    left: (String, char),
    right: (String, char),
) -> ((String, char), (String, char)) {
    if left <= right {
        (left, right)
    } else {
        (right, left)
    }
}

fn remove_mito_stable_physical_link_pair(
    links: &mut HashMap<SkeletonLinkKey, SkeletonLinkSupport>,
    key: &SkeletonLinkKey,
) {
    let physical = canonical_mito_stable_physical_link(key);
    links.retain(|candidate, _| canonical_mito_stable_physical_link(candidate) != physical);
}

fn apply_mito_stable_manual_edits(
    config: &Config,
    links: &mut HashMap<SkeletonLinkKey, SkeletonLinkSupport>,
) -> Vec<MitoStableManualEdit> {
    let mut edits = Vec::new();
    for key in &config.mito_stable_drop_links {
        remove_mito_stable_physical_link_pair(links, key);
        edits.push(MitoStableManualEdit {
            action: "drop".to_string(),
            key: key.clone(),
        });
    }
    for key in &config.mito_stable_force_links {
        insert_mito_stable_supported_link_pair(
            config,
            links,
            key.clone(),
            SkeletonLinkSupport {
                skeleton_support: config.min_link_support,
                paf_support: config.skeleton_min_link_support,
                rescue_support: config.skeleton_rescue_link_support,
                ..SkeletonLinkSupport::default()
            },
        );
        edits.push(MitoStableManualEdit {
            action: "force".to_string(),
            key: key.clone(),
        });
    }
    edits
}

fn analyze_mito_stable_copy_choices(
    config: &Config,
    segments: &[SkeletonSegment],
    base_links: &HashMap<SkeletonLinkKey, SkeletonLinkSupport>,
    candidates: &[MitoCompactBridgeCandidate],
) -> Vec<MitoStableCopyChoiceRow> {
    let degrees = mito_stable_physical_side_degrees(config, segments, base_links);
    let mut selected_three_way = if config.mito_stable_selected_nodes.is_empty() {
        segments
            .iter()
            .filter_map(|segment| {
                let (left, right) = degrees.get(&segment.name).copied().unwrap_or((0, 0));
                (mito_stable_node_degree_class(left, right) == "three_way")
                    .then(|| segment.name.clone())
            })
            .collect::<Vec<_>>()
    } else {
        config
            .mito_stable_selected_nodes
            .iter()
            .filter_map(|segment| {
                let (left, right) = degrees.get(segment).copied().unwrap_or((0, 0));
                (mito_stable_node_degree_class(left, right) == "three_way").then(|| segment.clone())
            })
            .collect::<Vec<_>>()
    };
    if selected_three_way.len() < 2 {
        return Vec::new();
    }
    selected_three_way.sort();

    let support_by_physical = mito_stable_candidate_support_by_physical(candidates);
    let physical_links = mito_stable_kept_physical_links(config, base_links);
    let mut rows = Vec::new();
    for left_index in 0..selected_three_way.len() {
        for right_index in (left_index + 1)..selected_three_way.len() {
            let node_a = &selected_three_way[left_index];
            let node_b = &selected_three_way[right_index];
            let Some((node_a_single, node_a_branch)) =
                mito_stable_three_way_single_and_branch_sides(&degrees, node_a)
            else {
                continue;
            };
            let Some((node_b_single, node_b_branch)) =
                mito_stable_three_way_single_and_branch_sides(&degrees, node_b)
            else {
                continue;
            };

            let double_link = skeleton_link_key_from_physical_ends(
                (node_a.clone(), node_a_single),
                (node_b.clone(), node_b_single),
            );
            let inferred_double_support = mito_stable_double_copy_bridge_chain_support(
                config,
                segments,
                base_links,
                &(node_a.clone(), node_a_single),
                &(node_b.clone(), node_b_single),
            );
            rows.push(mito_stable_copy_choice_row(
                config,
                segments,
                base_links,
                node_a,
                node_b,
                "double-copy-add-single-sides",
                &[],
                Some(double_link),
                &support_by_physical,
                inferred_double_support,
            ));

            let node_a_branch_links = mito_stable_links_touching_physical_side(
                &physical_links,
                &(node_a.clone(), node_a_branch),
            );
            let node_b_branch_links = mito_stable_links_touching_physical_side(
                &physical_links,
                &(node_b.clone(), node_b_branch),
            );
            for remove_a in &node_a_branch_links {
                for remove_b in &node_b_branch_links {
                    if canonical_physical_end_pair(remove_a.0.clone(), remove_a.1.clone())
                        == canonical_physical_end_pair(remove_b.0.clone(), remove_b.1.clone())
                    {
                        continue;
                    }
                    let Some(free_a) =
                        other_physical_end(remove_a, &(node_a.clone(), node_a_branch))
                    else {
                        continue;
                    };
                    let Some(free_b) =
                        other_physical_end(remove_b, &(node_b.clone(), node_b_branch))
                    else {
                        continue;
                    };
                    let removed = vec![
                        skeleton_link_key_from_physical_ends(
                            remove_a.0.clone(),
                            remove_a.1.clone(),
                        ),
                        skeleton_link_key_from_physical_ends(
                            remove_b.0.clone(),
                            remove_b.1.clone(),
                        ),
                    ];
                    let added = skeleton_link_key_from_physical_ends(free_a, free_b);
                    rows.push(mito_stable_copy_choice_row(
                        config,
                        segments,
                        base_links,
                        node_a,
                        node_b,
                        "single-copy-remove-branches-and-reconnect",
                        &removed,
                        Some(added),
                        &support_by_physical,
                        None,
                    ));
                }
            }
        }
    }
    rows.sort_by(|left, right| {
        left.node_a
            .cmp(&right.node_a)
            .then_with(|| left.node_b.cmp(&right.node_b))
            .then_with(|| left.option.cmp(&right.option))
            .then_with(|| right.added_support.cmp(&left.added_support))
            .then_with(|| left.added_link.cmp(&right.added_link))
    });
    rows
}

fn select_mito_stable_copy_choice_bridges(
    config: &Config,
    base_topology: &MitoStableTopologyStats,
    rows: &[MitoStableCopyChoiceRow],
) -> Vec<MitoCompactBridgeCandidate> {
    if config.mito_stable_selection_mode == MitoStableSelectionMode::AllowSelected {
        return Vec::new();
    }

    let min_double_support = config.skeleton_min_link_support.saturating_mul(3).max(30);
    let mut rows_by_pair: HashMap<(String, String), Vec<&MitoStableCopyChoiceRow>> = HashMap::new();
    for row in rows {
        rows_by_pair
            .entry((row.node_a.clone(), row.node_b.clone()))
            .or_default()
            .push(row);
    }

    let mut candidates = Vec::new();
    for ((node_a, node_b), pair_rows) in rows_by_pair {
        let Some(double_row) = pair_rows
            .iter()
            .copied()
            .find(|row| row.option == "double-copy-add-single-sides")
        else {
            continue;
        };
        let Some(added_link) = double_row.added_link.clone() else {
            continue;
        };
        if double_row.added_support < min_double_support {
            continue;
        }
        if double_row.node_a_class != "2-2" || double_row.node_b_class != "2-2" {
            continue;
        }
        if double_row.components > base_topology.components
            || double_row.open_ends > base_topology.open_ends
            || double_row.endpoint_overload > base_topology.endpoint_overload
        {
            continue;
        }

        let max_single_support = pair_rows
            .iter()
            .filter(|row| row.option == "single-copy-remove-branches-and-reconnect")
            .map(|row| row.added_support)
            .max()
            .unwrap_or(0);
        if max_single_support > config.skeleton_min_link_support
            && double_row.added_support < max_single_support.saturating_mul(3)
        {
            continue;
        }

        candidates.push((
            double_row.added_support,
            max_single_support,
            node_a,
            node_b,
            added_link,
        ));
    }

    candidates.sort_by(|left, right| {
        right
            .0
            .cmp(&left.0)
            .then_with(|| left.1.cmp(&right.1))
            .then_with(|| left.2.cmp(&right.2))
            .then_with(|| left.3.cmp(&right.3))
            .then_with(|| left.4.cmp(&right.4))
    });

    let mut used_nodes = HashSet::new();
    let mut seen_physical = HashSet::new();
    let mut selected = Vec::new();
    for (support, _max_single_support, node_a, node_b, key) in candidates {
        if used_nodes.contains(&node_a) || used_nodes.contains(&node_b) {
            continue;
        }
        if !seen_physical.insert(canonical_mito_stable_physical_link(&key)) {
            continue;
        }
        used_nodes.insert(node_a);
        used_nodes.insert(node_b);
        selected.push(MitoCompactBridgeCandidate {
            key,
            support,
            source: MitoCompactBridgeSource::CopyChoice,
        });
    }
    selected
}

fn analyze_mito_stable_selected_nodes(
    config: &Config,
    segments: &[SkeletonSegment],
    base_links: &HashMap<SkeletonLinkKey, SkeletonLinkSupport>,
    candidates: &[MitoCompactBridgeCandidate],
) -> Vec<MitoStableSelectedNodeDiagnostic> {
    let segment_names = segments
        .iter()
        .map(|segment| segment.name.clone())
        .collect::<HashSet<_>>();
    let degrees = mito_stable_physical_side_degrees(config, segments, base_links);
    let physical_links = mito_stable_kept_physical_links(config, base_links);
    let mut rows = Vec::new();
    let mut selected = config
        .mito_stable_selected_nodes
        .iter()
        .filter(|segment| segment_names.contains(*segment))
        .cloned()
        .collect::<Vec<_>>();
    selected.sort();
    for segment in selected {
        let (left_degree, right_degree) = degrees.get(&segment).copied().unwrap_or((0, 0));
        let mut incident_base_links = physical_links
            .iter()
            .filter(|(left, right)| left.0 == segment || right.0 == segment)
            .map(|(left, right)| format_physical_link(left, right))
            .collect::<Vec<_>>();
        incident_base_links.sort();

        let mut incident_candidates = candidates
            .iter()
            .filter(|candidate| candidate.key.from == segment || candidate.key.to == segment)
            .map(|candidate| {
                format!(
                    "{}:{}",
                    format_skeleton_link_key(&candidate.key),
                    candidate.support
                )
            })
            .collect::<Vec<_>>();
        incident_candidates.sort();

        rows.push(MitoStableSelectedNodeDiagnostic {
            segment,
            left_degree,
            right_degree,
            class: mito_stable_node_degree_class(left_degree, right_degree),
            selected: true,
            mode: config.mito_stable_selection_mode.as_str().to_string(),
            incident_base_links,
            incident_candidates,
        });
    }
    rows
}

fn mito_stable_copy_choice_row(
    config: &Config,
    segments: &[SkeletonSegment],
    base_links: &HashMap<SkeletonLinkKey, SkeletonLinkSupport>,
    node_a: &str,
    node_b: &str,
    option: &str,
    removed_links: &[SkeletonLinkKey],
    added_link: Option<SkeletonLinkKey>,
    support_by_physical: &HashMap<((String, char), (String, char)), u32>,
    inferred_added_support: Option<(u32, String)>,
) -> MitoStableCopyChoiceRow {
    let direct_support = added_link
        .as_ref()
        .and_then(|key| {
            support_by_physical
                .get(&canonical_mito_stable_physical_link(key))
                .copied()
        })
        .unwrap_or(0);
    let inferred_support = inferred_added_support
        .as_ref()
        .map(|(support, _)| *support)
        .unwrap_or(0);
    let added_support = direct_support.max(inferred_support);
    let added_support_source = match inferred_added_support {
        Some((support, source)) if support > direct_support => source,
        _ if direct_support > 0 => "read_walk".to_string(),
        _ => ".".to_string(),
    };

    let mut trial = base_links.clone();
    for key in removed_links {
        remove_mito_stable_physical_link_pair(&mut trial, key);
    }
    if let Some(key) = &added_link {
        insert_mito_stable_supported_link_pair(
            config,
            &mut trial,
            key.clone(),
            SkeletonLinkSupport {
                skeleton_support: config.min_link_support,
                paf_support: added_support,
                ..SkeletonLinkSupport::default()
            },
        );
    }
    refresh_skeleton_link_ratios(config, &mut trial);
    let topology = mito_stable_topology_stats(config, segments, &trial);
    let degrees = mito_stable_physical_side_degrees(config, segments, &trial);
    let (node_a_left, node_a_right) = degrees.get(node_a).copied().unwrap_or((0, 0));
    let (node_b_left, node_b_right) = degrees.get(node_b).copied().unwrap_or((0, 0));
    let node_a_class = mito_stable_node_degree_class(node_a_left, node_a_right);
    let node_b_class = mito_stable_node_degree_class(node_b_left, node_b_right);
    let interpretation = if option.starts_with("double-copy") {
        if added_support >= config.skeleton_min_link_support {
            "read_supported_double_copy_bridge"
        } else {
            "weak_double_copy_bridge"
        }
    } else if added_support <= 1 {
        "weak_single_copy_reconnect"
    } else if added_support >= config.skeleton_min_link_support {
        "read_supported_single_copy_reconnect"
    } else {
        "low_support_single_copy_reconnect"
    }
    .to_string();
    MitoStableCopyChoiceRow {
        node_a: node_a.to_string(),
        node_b: node_b.to_string(),
        option: option.to_string(),
        removed_links: removed_links.to_vec(),
        added_link,
        added_support,
        added_support_source,
        components: topology.components,
        open_ends: topology.open_ends,
        endpoint_overload: topology.endpoint_overload,
        branch_ends: topology.branch_ends,
        cycle_rank: topology.cycle_rank,
        node_a_class,
        node_b_class,
        interpretation,
    }
}

fn mito_stable_three_way_single_and_branch_sides(
    degrees: &HashMap<String, (usize, usize)>,
    segment: &str,
) -> Option<(char, char)> {
    let (left, right) = degrees.get(segment).copied().unwrap_or((0, 0));
    match (left, right) {
        (1, 2) => Some(('L', 'R')),
        (2, 1) => Some(('R', 'L')),
        _ => None,
    }
}

fn mito_stable_candidate_support_by_physical(
    candidates: &[MitoCompactBridgeCandidate],
) -> HashMap<((String, char), (String, char)), u32> {
    let mut support: HashMap<((String, char), (String, char)), u32> = HashMap::new();
    for candidate in candidates {
        support
            .entry(canonical_mito_stable_physical_link(&candidate.key))
            .and_modify(|value| *value = (*value).max(candidate.support))
            .or_insert(candidate.support);
    }
    support
}

fn mito_stable_double_copy_bridge_chain_support(
    config: &Config,
    segments: &[SkeletonSegment],
    links: &HashMap<SkeletonLinkKey, SkeletonLinkSupport>,
    left_end: &(String, char),
    right_end: &(String, char),
) -> Option<(u32, String)> {
    let degrees = mito_stable_physical_side_degrees(config, segments, links);
    let lengths: HashMap<String, usize> = segments
        .iter()
        .map(|segment| (segment.name.clone(), segment.sequence.len()))
        .collect();
    let physical_support = mito_stable_kept_physical_link_support_scores(config, links);
    let left_links =
        mito_stable_supported_links_touching_physical_side(&physical_support, left_end);
    let right_links =
        mito_stable_supported_links_touching_physical_side(&physical_support, right_end);
    let max_bridge_len = mito_stable_small_repeat_bridge_max_len(config);
    let mut best: Option<(u32, String, String)> = None;

    for (left_physical, left_support) in &left_links {
        let Some(left_bridge_end) = other_physical_end(left_physical, left_end) else {
            continue;
        };
        let bridge = left_bridge_end.0.clone();
        if bridge == left_end.0 || bridge == right_end.0 {
            continue;
        }
        let (bridge_left_degree, bridge_right_degree) =
            degrees.get(&bridge).copied().unwrap_or((0, 0));
        if mito_stable_node_degree_class(bridge_left_degree, bridge_right_degree) != "2-2" {
            continue;
        }
        if lengths.get(&bridge).copied().unwrap_or(usize::MAX) > max_bridge_len {
            continue;
        }

        for (right_physical, right_support) in &right_links {
            let Some(right_bridge_end) = other_physical_end(right_physical, right_end) else {
                continue;
            };
            if right_bridge_end.0 != bridge || right_bridge_end.1 == left_bridge_end.1 {
                continue;
            }
            let support = (*left_support).min(*right_support);
            let source = format!(
                "bridge_chain:{}:{}:{}",
                bridge,
                format_physical_link(&left_physical.0, &left_physical.1),
                format_physical_link(&right_physical.0, &right_physical.1)
            );
            let score = (support, bridge.clone(), source);
            if best.as_ref().is_none_or(|old| score > *old) {
                best = Some(score);
            }
        }
    }

    best.map(|(support, _bridge, source)| (support, source))
}

fn mito_stable_kept_physical_link_support_scores(
    config: &Config,
    links: &HashMap<SkeletonLinkKey, SkeletonLinkSupport>,
) -> HashMap<((String, char), (String, char)), u32> {
    let mut support_by_physical: HashMap<((String, char), (String, char)), u32> = HashMap::new();
    for (key, support) in links {
        if !keep_skeleton_link(config, support) {
            continue;
        }
        let score = mito_stable_link_support_score(support);
        support_by_physical
            .entry(canonical_mito_stable_physical_link(key))
            .and_modify(|old| *old = (*old).max(score))
            .or_insert(score);
    }
    support_by_physical
}

fn mito_stable_supported_links_touching_physical_side(
    physical_support: &HashMap<((String, char), (String, char)), u32>,
    endpoint: &(String, char),
) -> Vec<(((String, char), (String, char)), u32)> {
    physical_support
        .iter()
        .filter(|((left, right), _)| left == endpoint || right == endpoint)
        .map(|(physical, support)| (physical.clone(), *support))
        .collect()
}

fn mito_stable_link_support_score(support: &SkeletonLinkSupport) -> u32 {
    support
        .skeleton_support
        .max(support.paf_support)
        .max(support.rescue_support)
}

fn mito_stable_kept_physical_links(
    config: &Config,
    links: &HashMap<SkeletonLinkKey, SkeletonLinkSupport>,
) -> Vec<((String, char), (String, char))> {
    let mut seen = HashSet::new();
    let mut physical_links = Vec::new();
    for (key, support) in links {
        if !keep_skeleton_link(config, support) {
            continue;
        }
        let physical = canonical_mito_stable_physical_link(key);
        if seen.insert(physical.clone()) {
            physical_links.push(physical);
        }
    }
    physical_links.sort();
    physical_links
}

fn mito_stable_links_touching_physical_side(
    physical_links: &[((String, char), (String, char))],
    endpoint: &(String, char),
) -> Vec<((String, char), (String, char))> {
    physical_links
        .iter()
        .filter(|(left, right)| left == endpoint || right == endpoint)
        .cloned()
        .collect()
}

fn other_physical_end(
    physical_link: &((String, char), (String, char)),
    endpoint: &(String, char),
) -> Option<(String, char)> {
    if &physical_link.0 == endpoint {
        Some(physical_link.1.clone())
    } else if &physical_link.1 == endpoint {
        Some(physical_link.0.clone())
    } else {
        None
    }
}

fn skeleton_link_key_from_physical_ends(
    from_end: (String, char),
    to_end: (String, char),
) -> SkeletonLinkKey {
    SkeletonLinkKey {
        from: from_end.0,
        from_orient: if from_end.1 == 'R' { '+' } else { '-' },
        to: to_end.0,
        to_orient: if to_end.1 == 'L' { '+' } else { '-' },
    }
}

fn skeleton_physical_from_end(segment: &str, orient: char) -> (String, char) {
    (segment.to_string(), if orient == '+' { 'R' } else { 'L' })
}

fn skeleton_physical_to_end(segment: &str, orient: char) -> (String, char) {
    (segment.to_string(), if orient == '+' { 'L' } else { 'R' })
}

fn read_skeleton_gfa(
    path: &Path,
) -> io::Result<(
    Vec<SkeletonSegment>,
    HashMap<SkeletonLinkKey, SkeletonLinkSupport>,
)> {
    let reader = BufReader::new(File::open(path)?);
    let mut segments = Vec::new();
    let mut links: HashMap<SkeletonLinkKey, SkeletonLinkSupport> = HashMap::new();
    for line in reader.lines() {
        let line = line?;
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.is_empty() {
            continue;
        }
        match fields[0] {
            "S" if fields.len() >= 3 => segments.push(SkeletonSegment {
                name: fields[1].to_string(),
                sequence: fields[2].to_string(),
            }),
            "L" if fields.len() >= 5 => {
                let mut support = 0u32;
                for field in &fields[5..] {
                    if let Some(value) = field.strip_prefix("RC:i:") {
                        support = value.parse().unwrap_or(0);
                    }
                }
                let key = SkeletonLinkKey {
                    from: fields[1].to_string(),
                    from_orient: fields[2].chars().next().unwrap_or('+'),
                    to: fields[3].to_string(),
                    to_orient: fields[4].chars().next().unwrap_or('+'),
                };
                links.entry(key).or_default().skeleton_support = support;
            }
            _ => {}
        }
    }
    Ok((segments, links))
}

fn read_gfa_links(path: &Path) -> io::Result<HashMap<SkeletonLinkKey, u32>> {
    let reader = BufReader::new(File::open(path)?);
    let mut links: HashMap<SkeletonLinkKey, u32> = HashMap::new();
    for line in reader.lines() {
        let line = line?;
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 5 || fields[0] != "L" {
            continue;
        }
        let mut support = 0u32;
        for field in &fields[5..] {
            if let Some(value) = field.strip_prefix("RC:i:") {
                support = value.parse().unwrap_or(0);
            }
        }
        let key = SkeletonLinkKey {
            from: fields[1].to_string(),
            from_orient: fields[2].chars().next().unwrap_or('+'),
            to: fields[3].to_string(),
            to_orient: fields[4].chars().next().unwrap_or('+'),
        };
        links
            .entry(key)
            .and_modify(|old| *old = (*old).max(support))
            .or_insert(support);
    }
    Ok(links)
}

fn add_component_rescue_links(
    config: &Config,
    segments: &[SkeletonSegment],
    links: &mut HashMap<SkeletonLinkKey, SkeletonLinkSupport>,
    rescue_links: HashMap<SkeletonLinkKey, u32>,
) {
    let components = skeleton_components(config, segments, links);
    for (key, support) in rescue_links {
        if support < config.skeleton_rescue_link_support {
            continue;
        }
        let Some(&from_component) = components.get(&key.from) else {
            continue;
        };
        let Some(&to_component) = components.get(&key.to) else {
            continue;
        };
        if from_component == to_component {
            continue;
        }
        links.entry(key).or_default().rescue_support = support;
    }
}

fn skeleton_components(
    config: &Config,
    segments: &[SkeletonSegment],
    links: &HashMap<SkeletonLinkKey, SkeletonLinkSupport>,
) -> HashMap<String, usize> {
    let mut adj: HashMap<String, Vec<String>> = HashMap::new();
    for segment in segments {
        adj.entry(segment.name.clone()).or_default();
    }
    for (key, support) in links {
        if !keep_skeleton_link(config, support) {
            continue;
        }
        adj.entry(key.from.clone())
            .or_default()
            .push(key.to.clone());
        adj.entry(key.to.clone())
            .or_default()
            .push(key.from.clone());
    }
    let mut component = HashMap::new();
    let mut next_component = 0usize;
    for segment in segments {
        if component.contains_key(&segment.name) {
            continue;
        }
        next_component += 1;
        let mut stack = vec![segment.name.clone()];
        component.insert(segment.name.clone(), next_component);
        while let Some(node) = stack.pop() {
            if let Some(neighbors) = adj.get(&node) {
                for neighbor in neighbors {
                    if !component.contains_key(neighbor) {
                        component.insert(neighbor.clone(), next_component);
                        stack.push(neighbor.clone());
                    }
                }
            }
        }
    }
    component
}

fn read_filtered_skeleton_paf_by_read(
    config: &Config,
    segments: &[SkeletonSegment],
    paf_path: &Path,
) -> io::Result<HashMap<String, Vec<PafAln>>> {
    let lengths: HashSet<String> = segments
        .iter()
        .map(|segment| segment.name.clone())
        .collect();
    let mut by_read: HashMap<String, Vec<PafAln>> = HashMap::new();
    let reader = BufReader::new(File::open(paf_path)?);
    for line in reader.lines() {
        let line = line?;
        let Some(aln) = parse_paf_line(&line) else {
            continue;
        };
        if aln.identity() < config.minimap_min_identity || aln.alen < config.minimap_min_align_len {
            continue;
        }
        if !lengths.contains(&aln.tname) {
            continue;
        }
        by_read.entry(aln.qname.clone()).or_default().push(aln);
    }
    Ok(by_read)
}

fn complete_mito_compact_links(
    config: &Config,
    segments: &[SkeletonSegment],
    links: &mut HashMap<SkeletonLinkKey, SkeletonLinkSupport>,
    rescue_links: HashMap<SkeletonLinkKey, u32>,
    paf_by_read: &HashMap<String, Vec<PafAln>>,
) -> io::Result<MitoCompactBridgeReport> {
    let initial_stats = skeleton_graph_stats(config, segments, links);
    let mut report = MitoCompactBridgeReport {
        initial_components: initial_stats.components.iter().map(Vec::len).collect(),
        initial_open_roles: initial_stats.open_out.len() + initial_stats.open_in.len(),
        initial_kept_links: initial_stats.kept_links,
        ..MitoCompactBridgeReport::default()
    };

    let focused_reads = mito_compact_focused_reads(&initial_stats, paf_by_read);
    report.focused_reads = focused_reads.len();
    write_mito_compact_focused_reads(config, &focused_reads)?;

    let local_candidates =
        mito_compact_local_paf_candidates(config, &initial_stats, paf_by_read, &focused_reads);
    let rescue_candidates = mito_compact_rescue_candidates(&initial_stats, rescue_links);
    report.local_candidates = local_candidates.len();
    report.rescue_candidates = rescue_candidates.len();

    let mut candidates = Vec::with_capacity(local_candidates.len() + rescue_candidates.len());
    candidates.extend(local_candidates);
    candidates.extend(rescue_candidates);

    loop {
        let before = skeleton_graph_stats(config, segments, links);
        let mut best: Option<(usize, (usize, usize, usize, u32), SkeletonGraphStats)> = None;
        for (index, candidate) in candidates.iter().enumerate() {
            if candidate.support < mito_compact_candidate_min_support(config, candidate.source) {
                continue;
            }
            if skeleton_link_is_kept(config, links, &candidate.key) {
                continue;
            }
            let mut trial_links = links.clone();
            insert_mito_compact_bridge(config, &mut trial_links, candidate);
            refresh_skeleton_link_ratios(config, &mut trial_links);
            let after = skeleton_graph_stats(config, segments, &trial_links);
            let component_delta = before
                .components
                .len()
                .saturating_sub(after.components.len());
            let open_before = before.open_out.len() + before.open_in.len();
            let open_after = after.open_out.len() + after.open_in.len();
            let open_delta = open_before.saturating_sub(open_after);
            if component_delta == 0 && open_delta == 0 {
                continue;
            }
            if !mito_compact_candidate_degree_ok(&after, &candidate.key)
                && component_delta == 0
                && open_delta == 0
            {
                continue;
            }
            let score = (
                component_delta,
                open_delta,
                candidate.source.priority(),
                candidate.support,
            );
            if best
                .as_ref()
                .is_none_or(|(_, best_score, _)| score > *best_score)
            {
                best = Some((index, score, after));
            }
        }

        let Some((index, score, _after)) = best else {
            break;
        };
        let candidate = candidates.swap_remove(index);
        insert_mito_compact_bridge(config, links, &candidate);
        report.selected.push(MitoCompactSelectedBridge {
            candidate,
            component_delta: score.0,
            open_delta: score.1,
        });
    }

    Ok(report)
}

fn mito_compact_candidate_min_support(config: &Config, source: MitoCompactBridgeSource) -> u32 {
    match source {
        MitoCompactBridgeSource::LocalPaf | MitoCompactBridgeSource::CopyChoice => {
            config.skeleton_min_link_support
        }
        MitoCompactBridgeSource::Rescue => config.skeleton_rescue_link_support,
    }
}

fn mito_compact_focused_reads(
    stats: &SkeletonGraphStats,
    paf_by_read: &HashMap<String, Vec<PafAln>>,
) -> HashSet<String> {
    let mut focus_segments: HashSet<String> = HashSet::new();
    for (segment, _) in stats.open_out.iter().chain(stats.open_in.iter()) {
        focus_segments.insert(segment.clone());
    }
    for component in stats.components.iter().skip(1) {
        for segment in component {
            focus_segments.insert(segment.clone());
        }
    }

    let mut read_ids = HashSet::new();
    for (read, alns) in paf_by_read {
        if alns.iter().any(|aln| focus_segments.contains(&aln.tname)) {
            read_ids.insert(read.clone());
        }
    }
    read_ids
}

fn write_mito_compact_focused_reads(
    config: &Config,
    focused_reads: &HashSet<String>,
) -> io::Result<()> {
    let ids_path = config.out_dir.join("focused_read_ids.txt");
    let mut ids: Vec<_> = focused_reads.iter().collect();
    ids.sort();
    let mut id_out = File::create(ids_path)?;
    for read in &ids {
        writeln!(id_out, "{read}")?;
    }

    if !config.keep_debug_files {
        return Ok(());
    }

    let fasta_path = config.out_dir.join("focused_reads.fasta");
    let mut fasta = BufWriter::new(File::create(fasta_path)?);
    let mut read_index = 0u64;
    for path in &config.reads {
        read_sequence_file(path, |name, seq| {
            if focused_reads.contains(name) {
                write_fasta_record(&mut fasta, name, read_index, seq)?;
            }
            read_index += 1;
            Ok(true)
        })?;
    }
    fasta.flush()?;
    Ok(())
}

fn mito_compact_local_paf_candidates(
    config: &Config,
    stats: &SkeletonGraphStats,
    paf_by_read: &HashMap<String, Vec<PafAln>>,
    focused_reads: &HashSet<String>,
) -> Vec<MitoCompactBridgeCandidate> {
    let local_gap = mito_compact_local_bridge_gap(config);
    let mut support: HashMap<SkeletonLinkKey, u32> = HashMap::new();
    for read in focused_reads {
        let Some(alns) = paf_by_read.get(read) else {
            continue;
        };
        let chain = paf_non_overlapping_chain(alns.clone(), local_gap);
        let mut read_keys = HashSet::new();
        for left_index in 0..chain.len() {
            for right_index in left_index + 1..chain.len() {
                let left = &chain[left_index];
                let right = &chain[right_index];
                let gap = right.qstart as isize - left.qend as isize;
                if gap.abs() > local_gap {
                    continue;
                }
                let Some(from_orient) = paf_exit_orient(left, config.skeleton_end_slop) else {
                    continue;
                };
                let Some(to_orient) = paf_entry_orient(right, config.skeleton_end_slop) else {
                    continue;
                };
                if left.tname == right.tname && from_orient == to_orient {
                    continue;
                }
                let key = SkeletonLinkKey {
                    from: left.tname.clone(),
                    from_orient,
                    to: right.tname.clone(),
                    to_orient,
                };
                if mito_compact_candidate_is_relevant(stats, &key) {
                    read_keys.insert(key);
                }
            }
        }
        for key in read_keys {
            *support.entry(key).or_insert(0) += 1;
        }
    }

    support
        .into_iter()
        .map(|(key, support)| MitoCompactBridgeCandidate {
            key,
            support,
            source: MitoCompactBridgeSource::LocalPaf,
        })
        .collect()
}

fn mito_compact_rescue_candidates(
    stats: &SkeletonGraphStats,
    rescue_links: HashMap<SkeletonLinkKey, u32>,
) -> Vec<MitoCompactBridgeCandidate> {
    rescue_links
        .into_iter()
        .filter(|(key, _)| mito_compact_candidate_is_relevant(stats, key))
        .map(|(key, support)| MitoCompactBridgeCandidate {
            key,
            support,
            source: MitoCompactBridgeSource::Rescue,
        })
        .collect()
}

fn mito_compact_local_bridge_gap(config: &Config) -> isize {
    let local_gap = config
        .skeleton_end_slop
        .saturating_mul(MITO_COMPACT_LOCAL_BRIDGE_GAP_MULTIPLIER);
    config.paf_max_link_gap.max(local_gap as isize)
}

fn mito_compact_candidate_is_relevant(stats: &SkeletonGraphStats, key: &SkeletonLinkKey) -> bool {
    let from_component = stats.component_by_segment.get(&key.from);
    let to_component = stats.component_by_segment.get(&key.to);
    from_component
        .zip(to_component)
        .is_some_and(|(from, to)| from != to)
        || stats
            .open_out
            .contains(&(key.from.clone(), key.from_orient))
        || stats.open_in.contains(&(key.to.clone(), key.to_orient))
}

fn insert_mito_compact_bridge(
    config: &Config,
    links: &mut HashMap<SkeletonLinkKey, SkeletonLinkSupport>,
    candidate: &MitoCompactBridgeCandidate,
) {
    insert_mito_compact_bridge_key(links, &candidate.key, candidate.source, candidate.support);
    if config.bidirectional_links {
        let rc = reverse_skeleton_link_key(&candidate.key);
        insert_mito_compact_bridge_key(links, &rc, candidate.source, candidate.support);
    }
}

fn insert_mito_compact_bridge_key(
    links: &mut HashMap<SkeletonLinkKey, SkeletonLinkSupport>,
    key: &SkeletonLinkKey,
    source: MitoCompactBridgeSource,
    support: u32,
) {
    let entry = links.entry(key.clone()).or_default();
    match source {
        MitoCompactBridgeSource::LocalPaf | MitoCompactBridgeSource::CopyChoice => {
            entry.paf_support = entry.paf_support.max(support);
        }
        MitoCompactBridgeSource::Rescue => {
            entry.rescue_support = entry.rescue_support.max(support);
        }
    }
}

fn mito_compact_candidate_degree_ok(stats: &SkeletonGraphStats, key: &SkeletonLinkKey) -> bool {
    let rc = reverse_skeleton_link_key(key);
    skeleton_out_degree(stats, &key.from, key.from_orient) <= MITO_COMPACT_MAX_ENDPOINT_DEGREE
        && skeleton_in_degree(stats, &key.to, key.to_orient) <= MITO_COMPACT_MAX_ENDPOINT_DEGREE
        && skeleton_out_degree(stats, &rc.from, rc.from_orient) <= MITO_COMPACT_MAX_ENDPOINT_DEGREE
        && skeleton_in_degree(stats, &rc.to, rc.to_orient) <= MITO_COMPACT_MAX_ENDPOINT_DEGREE
}

fn skeleton_out_degree(stats: &SkeletonGraphStats, segment: &str, orient: char) -> usize {
    stats
        .out_degree
        .get(&(segment.to_string(), orient))
        .copied()
        .unwrap_or(0)
}

fn skeleton_in_degree(stats: &SkeletonGraphStats, segment: &str, orient: char) -> usize {
    stats
        .in_degree
        .get(&(segment.to_string(), orient))
        .copied()
        .unwrap_or(0)
}

fn skeleton_link_is_kept(
    config: &Config,
    links: &HashMap<SkeletonLinkKey, SkeletonLinkSupport>,
    key: &SkeletonLinkKey,
) -> bool {
    links
        .get(key)
        .is_some_and(|support| keep_skeleton_link(config, support))
}

fn skeleton_graph_stats(
    config: &Config,
    segments: &[SkeletonSegment],
    links: &HashMap<SkeletonLinkKey, SkeletonLinkSupport>,
) -> SkeletonGraphStats {
    let mut adj: HashMap<String, Vec<String>> = HashMap::new();
    let mut out_degree: HashMap<(String, char), usize> = HashMap::new();
    let mut in_degree: HashMap<(String, char), usize> = HashMap::new();
    let mut kept_links = 0usize;

    for segment in segments {
        adj.entry(segment.name.clone()).or_default();
    }
    for (key, support) in links {
        if !keep_skeleton_link(config, support) {
            continue;
        }
        kept_links += 1;
        adj.entry(key.from.clone())
            .or_default()
            .push(key.to.clone());
        adj.entry(key.to.clone())
            .or_default()
            .push(key.from.clone());
        *out_degree
            .entry((key.from.clone(), key.from_orient))
            .or_insert(0) += 1;
        *in_degree
            .entry((key.to.clone(), key.to_orient))
            .or_insert(0) += 1;
    }

    let mut components = Vec::new();
    let mut component_by_segment = HashMap::new();
    let mut seen = HashSet::new();
    for segment in segments {
        if seen.contains(&segment.name) {
            continue;
        }
        let component_id = components.len();
        let mut component = Vec::new();
        let mut stack = vec![segment.name.clone()];
        seen.insert(segment.name.clone());
        while let Some(node) = stack.pop() {
            component_by_segment.insert(node.clone(), component_id);
            component.push(node.clone());
            if let Some(neighbors) = adj.get(&node) {
                for neighbor in neighbors {
                    if seen.insert(neighbor.clone()) {
                        stack.push(neighbor.clone());
                    }
                }
            }
        }
        components.push(component);
    }

    components.sort_by(|left, right| {
        right
            .len()
            .cmp(&left.len())
            .then_with(|| left.first().cmp(&right.first()))
    });
    component_by_segment.clear();
    for (component_id, component) in components.iter().enumerate() {
        for segment in component {
            component_by_segment.insert(segment.clone(), component_id);
        }
    }

    let mut open_out = HashSet::new();
    let mut open_in = HashSet::new();
    for segment in segments {
        for orient in ['+', '-'] {
            let endpoint = (segment.name.clone(), orient);
            if out_degree.get(&endpoint).copied().unwrap_or(0) == 0 {
                open_out.insert(endpoint.clone());
            }
            if in_degree.get(&endpoint).copied().unwrap_or(0) == 0 {
                open_in.insert(endpoint);
            }
        }
    }

    SkeletonGraphStats {
        components,
        component_by_segment,
        open_out,
        open_in,
        out_degree,
        in_degree,
        kept_links,
    }
}

fn prune_mito_compact_secondary_links(
    config: &Config,
    segments: &[SkeletonSegment],
    links: &mut HashMap<SkeletonLinkKey, SkeletonLinkSupport>,
) -> Vec<SkeletonLinkKey> {
    let mut pruned = Vec::new();
    loop {
        let before = skeleton_graph_stats(config, segments, links);
        let before_largest = before.components.first().map_or(0, Vec::len);
        let before_open = before.open_out.len() + before.open_in.len();
        let overloaded = mito_compact_overloaded_endpoints(&before);
        if overloaded.is_empty() {
            break;
        }

        let mut candidates: Vec<_> = links
            .iter()
            .filter(|(key, support)| {
                keep_skeleton_link(config, support)
                    && support.skeleton_support == 0
                    && (overloaded.contains(&(key.from.clone(), key.from_orient, true))
                        || overloaded.contains(&(key.to.clone(), key.to_orient, false)))
            })
            .map(|(key, support)| {
                let source_penalty = if support.rescue_support > 0 && support.paf_support == 0 {
                    0
                } else {
                    1
                };
                let support_score = support
                    .paf_support
                    .max(support.rescue_support)
                    .max(support.skeleton_support);
                (key.clone(), source_penalty, support_score)
            })
            .collect();
        candidates.sort_by(|a, b| {
            a.1.cmp(&b.1)
                .then_with(|| a.2.cmp(&b.2))
                .then_with(|| a.0.cmp(&b.0))
        });

        let mut removed = None;
        for (key, _, _) in candidates {
            let mut trial = links.clone();
            remove_skeleton_link_pair(&mut trial, &key);
            let after = skeleton_graph_stats(config, segments, &trial);
            let after_largest = after.components.first().map_or(0, Vec::len);
            let after_open = after.open_out.len() + after.open_in.len();
            if after.components.len() <= before.components.len()
                && after_largest >= before_largest
                && after_open <= before_open
            {
                remove_skeleton_link_pair(links, &key);
                removed = Some(key);
                break;
            }
        }

        let Some(key) = removed else {
            break;
        };
        pruned.push(key);
    }
    pruned
}

fn mito_compact_overloaded_endpoints(stats: &SkeletonGraphStats) -> HashSet<(String, char, bool)> {
    let mut overloaded = HashSet::new();
    for ((segment, orient), degree) in &stats.out_degree {
        if *degree > MITO_COMPACT_MAX_ENDPOINT_DEGREE {
            overloaded.insert((segment.clone(), *orient, true));
        }
    }
    for ((segment, orient), degree) in &stats.in_degree {
        if *degree > MITO_COMPACT_MAX_ENDPOINT_DEGREE {
            overloaded.insert((segment.clone(), *orient, false));
        }
    }
    overloaded
}

fn remove_skeleton_link_pair(
    links: &mut HashMap<SkeletonLinkKey, SkeletonLinkSupport>,
    key: &SkeletonLinkKey,
) {
    links.remove(key);
    let rc = reverse_skeleton_link_key(key);
    links.remove(&rc);
}

fn mito_compact_final_segment_names(
    config: &Config,
    segments: &[SkeletonSegment],
    links: &HashMap<SkeletonLinkKey, SkeletonLinkSupport>,
    report: &mut MitoCompactBridgeReport,
) -> HashSet<String> {
    let stats = skeleton_graph_stats(config, segments, links);
    let lengths: HashMap<String, usize> = segments
        .iter()
        .map(|segment| (segment.name.clone(), segment.sequence.len()))
        .collect();
    let mut keep = HashSet::new();
    let Some(main_component) = stats.components.first() else {
        return keep;
    };
    for segment in main_component {
        keep.insert(segment.clone());
    }

    for component in stats.components.iter().skip(1) {
        let bases: usize = component
            .iter()
            .map(|segment| lengths.get(segment).copied().unwrap_or(0))
            .sum();
        let is_small = component.len() <= MITO_COMPACT_SMALL_COMPONENT_MAX_SEGMENTS
            || bases <= MITO_COMPACT_SMALL_COMPONENT_MAX_BASES;
        if is_small {
            report.pruned_segments.extend(component.iter().cloned());
        } else {
            keep.extend(component.iter().cloned());
        }
    }
    report.pruned_segments.sort();
    keep
}

fn filter_skeleton_links_to_segments(
    links: &HashMap<SkeletonLinkKey, SkeletonLinkSupport>,
    segment_names: &HashSet<String>,
) -> HashMap<SkeletonLinkKey, SkeletonLinkSupport> {
    links
        .iter()
        .filter(|(key, _)| segment_names.contains(&key.from) && segment_names.contains(&key.to))
        .map(|(key, support)| (key.clone(), support.clone()))
        .collect()
}

fn write_mito_compact_bridge_report(
    config: &Config,
    report: &MitoCompactBridgeReport,
) -> io::Result<()> {
    let mut out = File::create(config.out_dir.join("link_selection.report.txt"))?;
    writeln!(out, "link_selection_report")?;
    writeln!(out, "mode\tlocal_open_end_bridge")?;
    writeln!(
        out,
        "local_bridge_gap\t{}",
        mito_compact_local_bridge_gap(config)
    )?;
    writeln!(out, "initial_kept_links\t{}", report.initial_kept_links)?;
    writeln!(
        out,
        "initial_components\t{}",
        join_usize_list(&report.initial_components)
    )?;
    writeln!(out, "initial_open_roles\t{}", report.initial_open_roles)?;
    writeln!(out, "focused_reads\t{}", report.focused_reads)?;
    writeln!(out, "local_candidates\t{}", report.local_candidates)?;
    writeln!(out, "rescue_candidates\t{}", report.rescue_candidates)?;
    writeln!(out, "selected_bridges\t{}", report.selected.len())?;
    writeln!(out, "pruned_links\t{}", report.pruned_links.len())?;
    writeln!(
        out,
        "pruned_segments\t{}",
        join_string_list(&report.pruned_segments)
    )?;
    writeln!(out, "final_kept_links\t{}", report.final_kept_links)?;
    writeln!(
        out,
        "final_components\t{}",
        join_usize_list(&report.final_components)
    )?;
    writeln!(out, "final_open_roles\t{}", report.final_open_roles)?;

    let mut bridge_out = File::create(config.out_dir.join("selected_links.tsv"))?;
    writeln!(
        bridge_out,
        "from\tfrom_orient\tto\tto_orient\tsource\tsupport\tcomponent_delta\topen_delta"
    )?;
    for selected in &report.selected {
        let key = &selected.candidate.key;
        writeln!(
            bridge_out,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            key.from,
            key.from_orient,
            key.to,
            key.to_orient,
            selected.candidate.source.label(),
            selected.candidate.support,
            selected.component_delta,
            selected.open_delta
        )?;
    }

    let mut pruned_out = File::create(config.out_dir.join("pruned_links.tsv"))?;
    writeln!(pruned_out, "kind\tfrom\tfrom_orient\tto\tto_orient")?;
    for key in &report.pruned_links {
        writeln!(
            pruned_out,
            "link\t{}\t{}\t{}\t{}",
            key.from, key.from_orient, key.to, key.to_orient
        )?;
    }
    for segment in &report.pruned_segments {
        writeln!(pruned_out, "segment\t{segment}\t.\t.\t.")?;
    }
    Ok(())
}

fn write_mito_stable_bridge_report(
    config: &Config,
    initial: &MitoStableTopologyStats,
    final_stats: &MitoStableTopologyStats,
    candidates: &[MitoCompactBridgeCandidate],
    all_candidate_count: usize,
    auto_copy_choice_bridges: &[MitoCompactBridgeCandidate],
    selected: &[MitoCompactBridgeCandidate],
    pruned_links: &[SkeletonLinkKey],
    repeat_expansions: &[MitoStableRepeatExpansion],
    manual_edits: &[MitoStableManualEdit],
    default_node_shape_violations: usize,
    node_constraint_fallback: bool,
    final_node_shape_violations: usize,
) -> io::Result<()> {
    let mut out = File::create(config.out_dir.join("link_selection.report.txt"))?;
    writeln!(out, "link_selection_report")?;
    writeln!(out, "mode\tglobal_candidate_selection")?;
    writeln!(
        out,
        "selection_mode\t{}",
        config.mito_stable_selection_mode.as_str()
    )?;
    let mut selected_nodes = config
        .mito_stable_selected_nodes
        .iter()
        .cloned()
        .collect::<Vec<_>>();
    selected_nodes.sort();
    writeln!(
        out,
        "selected_nodes\t{}",
        if selected_nodes.is_empty() {
            ".".to_string()
        } else {
            join_string_list(&selected_nodes)
        }
    )?;
    writeln!(
        out,
        "node_constraints_explicit\t{}",
        mito_stable_has_explicit_three_way_constraints(config)
    )?;
    writeln!(
        out,
        "default_node_shape_violations\t{}",
        default_node_shape_violations
    )?;
    writeln!(
        out,
        "node_constraint_fallback\t{}",
        node_constraint_fallback
    )?;
    writeln!(
        out,
        "final_node_shape_violations\t{}",
        final_node_shape_violations
    )?;
    writeln!(out, "initial_components\t{}", initial.components)?;
    writeln!(out, "initial_open_ends\t{}", initial.open_ends)?;
    writeln!(
        out,
        "initial_endpoint_overload\t{}",
        initial.endpoint_overload
    )?;
    writeln!(out, "initial_unique_links\t{}", initial.unique_links)?;
    writeln!(out, "initial_branch_ends\t{}", initial.branch_ends)?;
    writeln!(out, "initial_cycle_rank\t{}", initial.cycle_rank)?;
    writeln!(out, "all_candidates\t{}", all_candidate_count)?;
    writeln!(out, "candidates\t{}", candidates.len())?;
    writeln!(
        out,
        "auto_copy_choice_bridges\t{}",
        auto_copy_choice_bridges.len()
    )?;
    writeln!(out, "selected_bridges\t{}", selected.len())?;
    writeln!(
        out,
        "manual_forced_links\t{}",
        manual_edits
            .iter()
            .filter(|edit| edit.action == "force")
            .count()
    )?;
    writeln!(
        out,
        "manual_dropped_links\t{}",
        manual_edits
            .iter()
            .filter(|edit| edit.action == "drop")
            .count()
    )?;
    writeln!(out, "repeat_expansions\t{}", repeat_expansions.len())?;
    writeln!(out, "pruned_redundant_links\t{}", pruned_links.len())?;
    writeln!(out, "final_components\t{}", final_stats.components)?;
    writeln!(out, "final_open_ends\t{}", final_stats.open_ends)?;
    writeln!(
        out,
        "final_endpoint_overload\t{}",
        final_stats.endpoint_overload
    )?;
    writeln!(out, "final_unique_links\t{}", final_stats.unique_links)?;
    writeln!(out, "final_branch_ends\t{}", final_stats.branch_ends)?;
    writeln!(out, "final_cycle_rank\t{}", final_stats.cycle_rank)?;

    let mut candidate_out = File::create(config.out_dir.join("link_candidates.tsv"))?;
    writeln!(
        candidate_out,
        "selected\tfrom\tfrom_orient\tto\tto_orient\tsource\tsupport"
    )?;
    let selected_keys: HashSet<_> = selected
        .iter()
        .map(|candidate| canonical_skeleton_link_key(&candidate.key))
        .collect();
    for candidate in candidates {
        let key = &candidate.key;
        let selected = selected_keys.contains(&canonical_skeleton_link_key(key));
        writeln!(
            candidate_out,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}",
            selected,
            key.from,
            key.from_orient,
            key.to,
            key.to_orient,
            candidate.source.label(),
            candidate.support
        )?;
    }

    let mut selected_out = File::create(config.out_dir.join("selected_links.tsv"))?;
    writeln!(
        selected_out,
        "from\tfrom_orient\tto\tto_orient\tsource\tsupport"
    )?;
    for candidate in selected {
        let key = &candidate.key;
        writeln!(
            selected_out,
            "{}\t{}\t{}\t{}\t{}\t{}",
            key.from,
            key.from_orient,
            key.to,
            key.to_orient,
            candidate.source.label(),
            candidate.support
        )?;
    }
    Ok(())
}

fn write_mito_stable_repeat_expansions(
    config: &Config,
    expansions: &[MitoStableRepeatExpansion],
) -> io::Result<()> {
    let mut out = File::create(config.out_dir.join("repeat_expansions.tsv"))?;
    writeln!(
        out,
        "template\tclone\tremoved_from\tremoved_from_orient\tremoved_to\tremoved_to_orient\tleft_from\tleft_from_orient\tleft_to\tleft_to_orient\tright_from\tright_from_orient\tright_to\tright_to_orient\tshortcut_overlap\tleft_overlap\tright_overlap\tsupport"
    )?;
    for expansion in expansions {
        writeln!(
            out,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            expansion.template,
            expansion.clone,
            expansion.removed_shortcut.from,
            expansion.removed_shortcut.from_orient,
            expansion.removed_shortcut.to,
            expansion.removed_shortcut.to_orient,
            expansion.left_link.from,
            expansion.left_link.from_orient,
            expansion.left_link.to,
            expansion.left_link.to_orient,
            expansion.right_link.from,
            expansion.right_link.from_orient,
            expansion.right_link.to,
            expansion.right_link.to_orient,
            expansion.shortcut_overlap,
            expansion.left_overlap,
            expansion.right_overlap,
            expansion.support
        )?;
    }
    Ok(())
}

fn write_mito_stable_pruned_links(
    config: &Config,
    pruned_links: &[SkeletonLinkKey],
) -> io::Result<()> {
    let mut out = File::create(config.out_dir.join("redundant_links_pruned.tsv"))?;
    writeln!(out, "from\tfrom_orient\tto\tto_orient")?;
    for key in pruned_links {
        writeln!(
            out,
            "{}\t{}\t{}\t{}",
            key.from, key.from_orient, key.to, key.to_orient
        )?;
    }
    Ok(())
}

fn write_mito_stable_manual_edits(
    config: &Config,
    manual_edits: &[MitoStableManualEdit],
) -> io::Result<()> {
    let mut out = File::create(config.out_dir.join("manual_link_edits.tsv"))?;
    writeln!(out, "action\tfrom\tfrom_orient\tto\tto_orient")?;
    for edit in manual_edits {
        writeln!(
            out,
            "{}\t{}\t{}\t{}\t{}",
            edit.action, edit.key.from, edit.key.from_orient, edit.key.to, edit.key.to_orient
        )?;
    }
    Ok(())
}

fn write_mito_stable_copy_choices(
    config: &Config,
    rows: &[MitoStableCopyChoiceRow],
) -> io::Result<()> {
    let mut out = File::create(config.out_dir.join("copy_choice_links.tsv"))?;
    writeln!(
        out,
        "node_a\tnode_b\toption\tremoved_links\tadded_link\tadded_support\tsupport_source\tcomponents\topen_ends\tendpoint_overload\tbranch_ends\tcycle_rank\tnode_a_class\tnode_b_class\tinterpretation"
    )?;
    for row in rows {
        let removed_links = if row.removed_links.is_empty() {
            ".".to_string()
        } else {
            row.removed_links
                .iter()
                .map(format_skeleton_link_key)
                .collect::<Vec<_>>()
                .join(",")
        };
        let added_link = row
            .added_link
            .as_ref()
            .map(format_skeleton_link_key)
            .unwrap_or_else(|| ".".to_string());
        writeln!(
            out,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            row.node_a,
            row.node_b,
            row.option,
            removed_links,
            added_link,
            row.added_support,
            row.added_support_source,
            row.components,
            row.open_ends,
            row.endpoint_overload,
            row.branch_ends,
            row.cycle_rank,
            row.node_a_class,
            row.node_b_class,
            row.interpretation
        )?;
    }
    Ok(())
}

fn write_mito_stable_selected_node_diagnostics(
    config: &Config,
    rows: &[MitoStableSelectedNodeDiagnostic],
) -> io::Result<()> {
    let mut out = File::create(config.out_dir.join("constrained_nodes.tsv"))?;
    writeln!(
        out,
        "segment\tselected\tmode\tleft_degree\tright_degree\tclass\tincident_base_links\tincident_candidates"
    )?;
    for row in rows {
        writeln!(
            out,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            row.segment,
            row.selected,
            row.mode,
            row.left_degree,
            row.right_degree,
            row.class,
            if row.incident_base_links.is_empty() {
                ".".to_string()
            } else {
                join_string_list(&row.incident_base_links)
            },
            if row.incident_candidates.is_empty() {
                ".".to_string()
            } else {
                join_string_list(&row.incident_candidates)
            }
        )?;
    }
    Ok(())
}

fn write_mito_stable_node_degrees(
    config: &Config,
    segments: &[SkeletonSegment],
    links: &HashMap<SkeletonLinkKey, SkeletonLinkSupport>,
) -> io::Result<()> {
    let degrees = mito_stable_physical_side_degrees(config, segments, links);
    let mut out = File::create(config.out_dir.join("node_degrees.tsv"))?;
    writeln!(out, "segment\tleft_degree\tright_degree\tclass")?;
    let mut rows: Vec<_> = degrees.into_iter().collect();
    rows.sort_by(|left, right| left.0.cmp(&right.0));
    for (segment, (left_degree, right_degree)) in rows {
        writeln!(
            out,
            "{}\t{}\t{}\t{}",
            segment,
            left_degree,
            right_degree,
            mito_stable_node_degree_class(left_degree, right_degree)
        )?;
    }
    Ok(())
}

fn write_mito_stable_node_repairs(
    config: &Config,
    segments: &[SkeletonSegment],
    before_links: &HashMap<SkeletonLinkKey, SkeletonLinkSupport>,
    after_links: &HashMap<SkeletonLinkKey, SkeletonLinkSupport>,
    selected: &[MitoCompactBridgeCandidate],
    manual_edits: &[MitoStableManualEdit],
    pruned_links: &[SkeletonLinkKey],
    repeat_expansions: &[MitoStableRepeatExpansion],
) -> io::Result<()> {
    let before = mito_stable_physical_side_degrees(config, segments, before_links);
    let after = mito_stable_physical_side_degrees(config, segments, after_links);
    let mut events: HashMap<String, Vec<String>> = HashMap::new();

    for candidate in selected {
        let event = format!(
            "selected_bridge:{}:{}:{}:{}:{}",
            candidate.key.from,
            candidate.key.from_orient,
            candidate.key.to,
            candidate.key.to_orient,
            candidate.support
        );
        record_mito_stable_link_event(&mut events, &candidate.key, event);
    }
    for edit in manual_edits {
        let event = format!(
            "manual_{}:{}",
            edit.action,
            format_skeleton_link_key(&edit.key)
        );
        record_mito_stable_link_event(&mut events, &edit.key, event);
    }
    for key in pruned_links {
        let event = format!("pruned:{}", format_skeleton_link_key(key));
        record_mito_stable_link_event(&mut events, key, event);
    }
    for expansion in repeat_expansions {
        events
            .entry(expansion.template.clone())
            .or_default()
            .push(format!(
                "repeat_expansion_template:{}",
                format_skeleton_link_key(&expansion.removed_shortcut)
            ));
        events
            .entry(expansion.clone.clone())
            .or_default()
            .push(format!(
                "repeat_expansion_clone:{}:{}",
                expansion.template, expansion.support
            ));
    }

    let mut out = File::create(config.out_dir.join("node_repairs.tsv"))?;
    writeln!(
        out,
        "segment\tbefore_left\tbefore_right\tbefore_class\tafter_left\tafter_right\tafter_class\trepair_events"
    )?;
    let mut sorted_segments = segments.to_vec();
    sorted_segments.sort_by(|left, right| left.name.cmp(&right.name));
    for segment in sorted_segments {
        let (before_left, before_right) = before.get(&segment.name).copied().unwrap_or((0, 0));
        let (after_left, after_right) = after.get(&segment.name).copied().unwrap_or((0, 0));
        let mut node_events = events.remove(&segment.name).unwrap_or_default();
        node_events.sort();
        node_events.dedup();
        if before_left == after_left && before_right == after_right && node_events.is_empty() {
            continue;
        }
        writeln!(
            out,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            segment.name,
            before_left,
            before_right,
            mito_stable_node_degree_class(before_left, before_right),
            after_left,
            after_right,
            mito_stable_node_degree_class(after_left, after_right),
            if node_events.is_empty() {
                ".".to_string()
            } else {
                join_string_list(&node_events)
            }
        )?;
    }
    Ok(())
}

fn record_mito_stable_link_event(
    events: &mut HashMap<String, Vec<String>>,
    key: &SkeletonLinkKey,
    event: String,
) {
    events
        .entry(key.from.clone())
        .or_default()
        .push(event.clone());
    events.entry(key.to.clone()).or_default().push(event);
}

fn format_skeleton_link_key(key: &SkeletonLinkKey) -> String {
    format!(
        "{}:{}:{}:{}",
        key.from, key.from_orient, key.to, key.to_orient
    )
}

fn format_physical_link(left: &(String, char), right: &(String, char)) -> String {
    format!("{}:{}--{}:{}", left.0, left.1, right.0, right.1)
}

fn write_mito_stable_topology_scan(
    config: &Config,
    segments: &[SkeletonSegment],
    links: &HashMap<SkeletonLinkKey, SkeletonLinkSupport>,
) -> io::Result<()> {
    let degrees = mito_stable_physical_side_degrees(config, segments, links);
    let adjacency = mito_stable_segment_adjacency(config, segments, links);
    let merge_evidence = mito_stable_three_way_merge_evidence(config, segments, links);
    let bridge_evidence = mito_stable_small_repeat_bridge_evidence(config, segments, links);
    let mut out = File::create(config.out_dir.join("topology_scan.tsv"))?;
    writeln!(
        out,
        "segment\tlength\tleft_degree\tright_degree\tclass\tauto_role\tevidence\tneighbors"
    )?;
    let mut sorted_segments = segments.to_vec();
    sorted_segments.sort_by(|left, right| left.name.cmp(&right.name));
    for segment in sorted_segments {
        let (left_degree, right_degree) = degrees.get(&segment.name).copied().unwrap_or((0, 0));
        let class = mito_stable_node_degree_class(left_degree, right_degree);
        let mut neighbors = adjacency
            .get(&segment.name)
            .cloned()
            .unwrap_or_default()
            .into_iter()
            .collect::<Vec<_>>();
        neighbors.sort();
        let (auto_role, evidence) = if let Some(evidence) = merge_evidence.get(&segment.name) {
            (
                "accepted_three_way_merge_candidate",
                join_string_list(evidence),
            )
        } else if let Some(evidence) = bridge_evidence.get(&segment.name) {
            ("small_2_2_repeat_bridge", join_string_list(evidence))
        } else {
            match class {
                "three_way" => ("unresolved_three_way", ".".to_string()),
                "2-2" => ("repeat_2_2", ".".to_string()),
                "1-1" => ("linear", ".".to_string()),
                _ => ("invalid_shape", ".".to_string()),
            }
        };
        writeln!(
            out,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            segment.name,
            segment.sequence.len(),
            left_degree,
            right_degree,
            class,
            auto_role,
            evidence,
            join_string_list(&neighbors)
        )?;
    }
    Ok(())
}

fn mito_stable_segment_adjacency(
    config: &Config,
    segments: &[SkeletonSegment],
    links: &HashMap<SkeletonLinkKey, SkeletonLinkSupport>,
) -> HashMap<String, HashSet<String>> {
    let segment_names: HashSet<_> = segments
        .iter()
        .map(|segment| segment.name.clone())
        .collect();
    let mut adjacency: HashMap<String, HashSet<String>> = segments
        .iter()
        .map(|segment| (segment.name.clone(), HashSet::new()))
        .collect();
    let mut seen = HashSet::new();
    for (key, support) in links {
        if !keep_skeleton_link(config, support) {
            continue;
        }
        if !segment_names.contains(&key.from) || !segment_names.contains(&key.to) {
            continue;
        }
        let physical = canonical_mito_stable_physical_link(key);
        if !seen.insert(physical) {
            continue;
        }
        adjacency
            .entry(key.from.clone())
            .or_default()
            .insert(key.to.clone());
        adjacency
            .entry(key.to.clone())
            .or_default()
            .insert(key.from.clone());
    }
    adjacency
}

fn mito_stable_three_way_merge_evidence(
    config: &Config,
    segments: &[SkeletonSegment],
    links: &HashMap<SkeletonLinkKey, SkeletonLinkSupport>,
) -> HashMap<String, Vec<String>> {
    let bridge_evidence = mito_stable_small_repeat_bridge_evidence(config, segments, links);
    let mut evidence: HashMap<String, Vec<String>> = HashMap::new();
    for (bridge, three_way_neighbors) in bridge_evidence {
        for three_way in &three_way_neighbors {
            let partners = three_way_neighbors
                .iter()
                .filter(|partner| *partner != three_way)
                .cloned()
                .collect::<Vec<_>>();
            if partners.is_empty() {
                continue;
            }
            evidence.entry(three_way.clone()).or_default().push(format!(
                "bridge={bridge};partners={}",
                join_string_list(&partners)
            ));
        }
    }
    for values in evidence.values_mut() {
        values.sort();
        values.dedup();
    }
    evidence
}

fn mito_stable_small_repeat_bridge_evidence(
    config: &Config,
    segments: &[SkeletonSegment],
    links: &HashMap<SkeletonLinkKey, SkeletonLinkSupport>,
) -> HashMap<String, Vec<String>> {
    let degrees = mito_stable_physical_side_degrees(config, segments, links);
    let adjacency = mito_stable_segment_adjacency(config, segments, links);
    let lengths: HashMap<String, usize> = segments
        .iter()
        .map(|segment| (segment.name.clone(), segment.sequence.len()))
        .collect();
    let max_bridge_len = mito_stable_small_repeat_bridge_max_len(config);
    let mut evidence = HashMap::new();
    for segment in segments {
        let (left_degree, right_degree) = degrees.get(&segment.name).copied().unwrap_or((0, 0));
        if mito_stable_node_degree_class(left_degree, right_degree) != "2-2" {
            continue;
        }
        if lengths.get(&segment.name).copied().unwrap_or(usize::MAX) > max_bridge_len {
            continue;
        }
        let mut three_way_neighbors = adjacency
            .get(&segment.name)
            .into_iter()
            .flat_map(|neighbors| neighbors.iter())
            .filter_map(|neighbor| {
                let (left_degree, right_degree) = degrees.get(neighbor).copied().unwrap_or((0, 0));
                (mito_stable_node_degree_class(left_degree, right_degree) == "three_way")
                    .then(|| neighbor.clone())
            })
            .collect::<Vec<_>>();
        three_way_neighbors.sort();
        three_way_neighbors.dedup();
        if three_way_neighbors.len() >= 2 {
            evidence.insert(segment.name.clone(), three_way_neighbors);
        }
    }
    evidence
}

fn mito_stable_small_repeat_bridge_max_len(config: &Config) -> usize {
    config.min_tip_len.max(5_000)
}

fn write_mito_stable_split_report(
    config: &Config,
    split_report: &[MitoStableSplitReportRow],
) -> io::Result<()> {
    let mut out = File::create(config.out_dir.join("segment_splits.tsv"))?;
    writeln!(out, "segment\tposition\tsupport")?;
    for row in split_report {
        writeln!(out, "{}\t{}\t{}", row.segment, row.position, row.support)?;
    }
    Ok(())
}

fn join_usize_list(values: &[usize]) -> String {
    values
        .iter()
        .map(|value| value.to_string())
        .collect::<Vec<_>>()
        .join(",")
}

fn join_string_list(values: &[String]) -> String {
    values.join(",")
}

fn write_skeleton_fasta(path: &Path, segments: &[SkeletonSegment]) -> io::Result<()> {
    let mut out = File::create(path)?;
    for segment in segments {
        writeln!(out, ">{}", segment.name)?;
        write_wrapped(&mut out, &segment.sequence, 80)?;
    }
    Ok(())
}

fn run_minimap2_to_target(
    config: &Config,
    target: &Path,
    paf_path: &Path,
    log_path: &Path,
) -> io::Result<()> {
    let paf = File::create(paf_path)?;
    let mut command = Command::new("minimap2");
    command
        .arg("-x")
        .arg("map-hifi")
        .arg("-t")
        .arg(config.threads.to_string())
        .arg(target);
    for read in &config.reads {
        command.arg(read);
    }
    let output = command
        .stdout(Stdio::from(paf))
        .stderr(Stdio::piped())
        .output();
    match output {
        Ok(output) if output.status.success() => {
            fs::write(log_path, &output.stderr)?;
            Ok(())
        }
        Ok(output) => {
            fs::write(log_path, &output.stderr)?;
            Err(io::Error::other(format!(
                "minimap2 failed; see {}",
                log_path.display()
            )))
        }
        Err(err) if err.kind() == io::ErrorKind::NotFound => Err(io::Error::new(
            io::ErrorKind::NotFound,
            "minimap2 not found in PATH",
        )),
        Err(err) => Err(err),
    }
}

fn summarize_skeleton_paf(
    config: &Config,
    segments: &[SkeletonSegment],
    paf_path: &Path,
) -> io::Result<(
    HashMap<String, (usize, usize)>,
    HashMap<SkeletonLinkKey, u32>,
)> {
    let lengths: HashMap<String, usize> = segments
        .iter()
        .map(|segment| (segment.name.clone(), segment.sequence.len()))
        .collect();
    let mut depth: HashMap<String, (usize, usize)> = HashMap::new();
    let mut by_read: HashMap<String, Vec<PafAln>> = HashMap::new();
    let reader = BufReader::new(File::open(paf_path)?);
    for line in reader.lines() {
        let line = line?;
        let Some(aln) = parse_paf_line(&line) else {
            continue;
        };
        if aln.identity() < config.minimap_min_identity || aln.alen < config.minimap_min_align_len {
            continue;
        }
        if !lengths.contains_key(&aln.tname) {
            continue;
        }
        let entry = depth.entry(aln.tname.clone()).or_insert((0, 0));
        entry.0 += aln.tend.saturating_sub(aln.tstart);
        entry.1 += 1;
        by_read.entry(aln.qname.clone()).or_default().push(aln);
    }

    let mut links: HashMap<SkeletonLinkKey, u32> = HashMap::new();
    for (_read, mut alns) in by_read {
        alns.sort_by(|a, b| {
            a.qstart
                .cmp(&b.qstart)
                .then_with(|| b.alen.cmp(&a.alen))
                .then_with(|| {
                    b.identity()
                        .partial_cmp(&a.identity())
                        .unwrap_or(std::cmp::Ordering::Equal)
                })
        });
        let chain = paf_non_overlapping_chain(alns, config.paf_max_link_gap);
        for pair in chain.windows(2) {
            let left = &pair[0];
            let right = &pair[1];
            let gap = right.qstart as isize - left.qend as isize;
            if gap.abs() > config.paf_max_link_gap {
                continue;
            }
            let Some(from_orient) = paf_exit_orient(left, config.skeleton_end_slop) else {
                continue;
            };
            let Some(to_orient) = paf_entry_orient(right, config.skeleton_end_slop) else {
                continue;
            };
            if left.tname == right.tname && from_orient == to_orient {
                continue;
            }
            let key = SkeletonLinkKey {
                from: left.tname.clone(),
                from_orient,
                to: right.tname.clone(),
                to_orient,
            };
            *links.entry(key).or_insert(0) += 1;
        }
    }
    Ok((depth, links))
}

fn paf_non_overlapping_chain(mut alns: Vec<PafAln>, max_gap: isize) -> Vec<PafAln> {
    let mut chain: Vec<PafAln> = Vec::new();
    alns.sort_by(|a, b| {
        a.qstart
            .cmp(&b.qstart)
            .then_with(|| b.alen.cmp(&a.alen))
            .then_with(|| {
                b.identity()
                    .partial_cmp(&a.identity())
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
    });
    for aln in alns {
        if let Some(last) = chain.last() {
            let overlap = last.qend as isize - aln.qstart as isize;
            if overlap > max_gap {
                continue;
            }
        }
        chain.push(aln);
    }
    chain
}

fn paf_exit_orient(aln: &PafAln, slop: usize) -> Option<char> {
    match aln.strand {
        '+' if aln.tlen.saturating_sub(aln.tend) <= slop => Some('+'),
        '-' if aln.tstart <= slop => Some('-'),
        _ => None,
    }
}

fn paf_entry_orient(aln: &PafAln, slop: usize) -> Option<char> {
    match aln.strand {
        '+' if aln.tstart <= slop => Some('+'),
        '-' if aln.tlen.saturating_sub(aln.tend) <= slop => Some('-'),
        _ => None,
    }
}

fn reverse_skeleton_link_key(key: &SkeletonLinkKey) -> SkeletonLinkKey {
    SkeletonLinkKey {
        from: key.to.clone(),
        from_orient: flip_orient(key.to_orient),
        to: key.from.clone(),
        to_orient: flip_orient(key.from_orient),
    }
}

fn merge_skeleton_links(
    config: &Config,
    mut links: HashMap<SkeletonLinkKey, SkeletonLinkSupport>,
    paf_links: HashMap<SkeletonLinkKey, u32>,
) -> HashMap<SkeletonLinkKey, SkeletonLinkSupport> {
    for (key, support) in paf_links {
        links.entry(key).or_default().paf_support += support;
    }
    refresh_skeleton_link_ratios(config, &mut links);
    add_bidirectional_skeleton_links(config, &mut links);
    links
}

fn refresh_skeleton_link_ratios(
    _config: &Config,
    links: &mut HashMap<SkeletonLinkKey, SkeletonLinkSupport>,
) {
    let mut out_max: HashMap<(String, char), u32> = HashMap::new();
    let mut in_max: HashMap<(String, char), u32> = HashMap::new();
    for (key, support) in links.iter() {
        if support.paf_support == 0 {
            continue;
        }
        out_max
            .entry((key.from.clone(), key.from_orient))
            .and_modify(|max| *max = (*max).max(support.paf_support))
            .or_insert(support.paf_support);
        in_max
            .entry((key.to.clone(), key.to_orient))
            .and_modify(|max| *max = (*max).max(support.paf_support))
            .or_insert(support.paf_support);
    }
    for (key, support) in links.iter_mut() {
        let out_best = out_max
            .get(&(key.from.clone(), key.from_orient))
            .copied()
            .unwrap_or(support.paf_support.max(1));
        let in_best = in_max
            .get(&(key.to.clone(), key.to_orient))
            .copied()
            .unwrap_or(support.paf_support.max(1));
        support.out_ratio = support.paf_support as f64 / out_best.max(1) as f64;
        support.in_ratio = support.paf_support as f64 / in_best.max(1) as f64;
    }
}

fn add_bidirectional_skeleton_links(
    config: &Config,
    links: &mut HashMap<SkeletonLinkKey, SkeletonLinkSupport>,
) {
    if config.bidirectional_links {
        let mut extra = Vec::new();
        for (key, support) in links.iter() {
            let rc = reverse_skeleton_link_key(key);
            if !links.contains_key(&rc) {
                extra.push((rc, support.clone()));
            }
        }
        for (key, support) in extra {
            links.insert(key, support);
        }
    }
}

fn write_skeleton_depth(
    path: &Path,
    segments: &[SkeletonSegment],
    depth: &HashMap<String, (usize, usize)>,
) -> io::Result<()> {
    let mut out = File::create(path)?;
    writeln!(
        out,
        "segment\tlength\taligned_bases\tmean_depth\talignments"
    )?;
    for segment in segments {
        let (bases, count) = depth.get(&segment.name).copied().unwrap_or((0, 0));
        let mean = if segment.sequence.is_empty() {
            0.0
        } else {
            bases as f64 / segment.sequence.len() as f64
        };
        writeln!(
            out,
            "{}\t{}\t{}\t{:.3}\t{}",
            segment.name,
            segment.sequence.len(),
            bases,
            mean,
            count
        )?;
    }
    Ok(())
}

fn write_skeleton_links(
    path: &Path,
    links: &HashMap<SkeletonLinkKey, SkeletonLinkSupport>,
) -> io::Result<()> {
    let mut out = File::create(path)?;
    writeln!(
        out,
        "from\tfrom_orient\tto\tto_orient\tskeleton_support\tpaf_support\trescue_support\tout_ratio\tin_ratio"
    )?;
    let mut rows: Vec<_> = links.iter().collect();
    rows.sort_by(|a, b| {
        b.1.paf_support
            .cmp(&a.1.paf_support)
            .then_with(|| b.1.skeleton_support.cmp(&a.1.skeleton_support))
            .then_with(|| b.1.rescue_support.cmp(&a.1.rescue_support))
            .then_with(|| a.0.cmp(b.0))
    });
    for (key, support) in rows {
        writeln!(
            out,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.4}\t{:.4}",
            key.from,
            key.from_orient,
            key.to,
            key.to_orient,
            support.skeleton_support,
            support.paf_support,
            support.rescue_support,
            support.out_ratio,
            support.in_ratio
        )?;
    }
    Ok(())
}

fn write_skeleton_linked_gfa(
    config: &Config,
    path: &Path,
    segments: &[SkeletonSegment],
    depth: &HashMap<String, (usize, usize)>,
    links: &HashMap<SkeletonLinkKey, SkeletonLinkSupport>,
) -> io::Result<()> {
    let mut out = File::create(path)?;
    writeln!(out, "H\tVN:Z:1.0\tPG:Z:orgraft_asm\tST:Z:skeleton_linked")?;
    for segment in segments {
        let (aligned_bases, alignments) = depth.get(&segment.name).copied().unwrap_or((0, 0));
        let mean_depth = if segment.sequence.is_empty() {
            0.0
        } else {
            aligned_bases as f64 / segment.sequence.len() as f64
        };
        writeln!(
            out,
            "S\t{}\t{}\tLN:i:{}\tDP:f:{:.3}\tAB:i:{}\tAC:i:{}",
            segment.name,
            segment.sequence,
            segment.sequence.len(),
            mean_depth,
            aligned_bases,
            alignments
        )?;
    }
    let mut rows: Vec<_> = links.iter().collect();
    rows.sort_by(|a, b| a.0.cmp(b.0));
    for (key, support) in rows {
        if !keep_skeleton_link(config, support) {
            continue;
        }
        let ratio = support.out_ratio.min(support.in_ratio);
        let rc = support
            .skeleton_support
            .max(support.paf_support)
            .max(support.rescue_support);
        let source = match (support.skeleton_support > 0, support.paf_support > 0) {
            (true, true) if support.rescue_support > 0 => "skeleton+paf+rescue",
            (true, true) => "skeleton+paf",
            (true, false) if support.rescue_support > 0 => "skeleton+rescue",
            (true, false) => "skeleton",
            (false, true) if support.rescue_support > 0 => "paf+rescue",
            (false, true) => "paf",
            (false, false) if support.rescue_support > 0 => "rescue",
            (false, false) => "none",
        };
        writeln!(
            out,
            "L\t{}\t{}\t{}\t{}\t0M\tRC:i:{}\tSK:i:{}\tPA:i:{}\tRS:i:{}\tLR:f:{:.4}\tSC:Z:{}",
            key.from,
            key.from_orient,
            key.to,
            key.to_orient,
            rc,
            support.skeleton_support,
            support.paf_support,
            support.rescue_support,
            ratio,
            source
        )?;
    }
    Ok(())
}

fn write_skeleton_report(
    config: &Config,
    skeleton_gfa: &Path,
    segments: &[SkeletonSegment],
    links: &HashMap<SkeletonLinkKey, SkeletonLinkSupport>,
) -> io::Result<()> {
    let mut out = File::create(config.out_dir.join("linking.report.txt"))?;
    let kept = links
        .values()
        .filter(|support| keep_skeleton_link(config, support))
        .count();
    writeln!(out, "skeleton_gfa\t{}", skeleton_gfa.display())?;
    writeln!(out, "segments\t{}", segments.len())?;
    writeln!(out, "links_total\t{}", links.len())?;
    writeln!(out, "links_kept\t{}", kept)?;
    writeln!(out, "skeleton_end_slop\t{}", config.skeleton_end_slop)?;
    writeln!(
        out,
        "skeleton_min_link_support\t{}",
        config.skeleton_min_link_support
    )?;
    writeln!(
        out,
        "skeleton_min_link_ratio\t{}",
        config.skeleton_min_link_ratio
    )?;
    if let Some(rescue_gfa) = &config.skeleton_rescue_gfa {
        writeln!(out, "skeleton_rescue_gfa\t{}", rescue_gfa.display())?;
    }
    writeln!(
        out,
        "skeleton_rescue_link_support\t{}",
        config.skeleton_rescue_link_support
    )?;
    writeln!(out, "minimap_min_identity\t{}", config.minimap_min_identity)?;
    writeln!(
        out,
        "minimap_min_align_len\t{}",
        config.minimap_min_align_len
    )?;
    writeln!(out, "paf_max_link_gap\t{}", config.paf_max_link_gap)?;
    writeln!(out, "output_gfa\tgraph.gfa")?;
    Ok(())
}

fn keep_skeleton_link(config: &Config, support: &SkeletonLinkSupport) -> bool {
    support.skeleton_support >= config.min_link_support
        || (support.paf_support >= config.skeleton_min_link_support
            && support.out_ratio.min(support.in_ratio) >= config.skeleton_min_link_ratio)
        || support.rescue_support >= config.skeleton_rescue_link_support
}

fn flip_orient(orient: char) -> char {
    if orient == '+' {
        '-'
    } else {
        '+'
    }
}

fn run_minimap2(config: &Config) -> io::Result<()> {
    let jobs = [
        (
            "dominant",
            config.out_dir.join("unitigs.fasta"),
            config.out_dir.join("read_to_unitigs.paf"),
            config.out_dir.join("read_to_unitigs.minimap2.log"),
            config.out_dir.join("depth.minimap2.tsv"),
            config.out_dir.join("junctions.minimap2.tsv"),
        ),
        (
            "full",
            config.out_dir.join("unitigs.full.fasta"),
            config.out_dir.join("read_to_unitigs.full.paf"),
            config.out_dir.join("read_to_unitigs.full.minimap2.log"),
            config.out_dir.join("depth.minimap2.full.tsv"),
            config.out_dir.join("junctions.minimap2.full.tsv"),
        ),
    ];

    for (label, unitigs, paf_path, log_path, depth_path, junction_path) in jobs {
        if fasta_is_empty(&unitigs)? {
            continue;
        }
        let paf = File::create(&paf_path)?;
        let mut command = Command::new("minimap2");
        command
            .arg("-x")
            .arg("map-hifi")
            .arg("-t")
            .arg(config.threads.to_string())
            .arg(&unitigs);
        for read in &config.reads {
            command.arg(read);
        }
        let output = command
            .stdout(Stdio::from(paf))
            .stderr(Stdio::piped())
            .output();
        match output {
            Ok(output) if output.status.success() => {
                fs::write(&log_path, &output.stderr)?;
                summarize_paf(config, &paf_path, &unitigs, &depth_path, &junction_path)?;
            }
            Ok(output) => {
                fs::write(&log_path, &output.stderr)?;
                return Err(io::Error::other(format!(
                    "minimap2 failed for {label}; see {}",
                    log_path.display()
                )));
            }
            Err(err) if err.kind() == io::ErrorKind::NotFound => {
                let message = "warning: --minimap2 requested, but minimap2 was not found in PATH\n";
                fs::write(&log_path, message)?;
                eprint!("{message}");
            }
            Err(err) => return Err(err),
        }
    }
    Ok(())
}

fn fasta_is_empty(path: &Path) -> io::Result<bool> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    for line in reader.lines() {
        if line?.starts_with('>') {
            return Ok(false);
        }
    }
    Ok(true)
}

fn summarize_paf(
    config: &Config,
    paf_path: &Path,
    target_fasta: &Path,
    depth_path: &Path,
    junction_path: &Path,
) -> io::Result<()> {
    let mut target_lengths = read_fasta_lengths(target_fasta)?;
    let mut depth_bases: HashMap<String, usize> = HashMap::new();
    let mut alignment_counts: HashMap<String, usize> = HashMap::new();
    let mut by_read: HashMap<String, Vec<PafAln>> = HashMap::new();

    let reader = BufReader::new(File::open(paf_path)?);
    for line in reader.lines() {
        let line = line?;
        let Some(aln) = parse_paf_line(&line) else {
            continue;
        };
        if aln.identity() < config.minimap_min_identity || aln.alen < config.minimap_min_align_len {
            continue;
        }
        target_lengths.entry(aln.tname.clone()).or_insert(aln.tlen);
        *depth_bases.entry(aln.tname.clone()).or_insert(0) += aln.tend.saturating_sub(aln.tstart);
        *alignment_counts.entry(aln.tname.clone()).or_insert(0) += 1;
        by_read.entry(aln.qname.clone()).or_default().push(aln);
    }

    let mut depth_out = File::create(depth_path)?;
    writeln!(
        depth_out,
        "unitig\tlength\taligned_bases\tmean_depth\talignments"
    )?;
    let mut names: Vec<_> = target_lengths.keys().cloned().collect();
    names.sort();
    for name in names {
        let len = target_lengths.get(&name).copied().unwrap_or(0);
        let bases = depth_bases.get(&name).copied().unwrap_or(0);
        let count = alignment_counts.get(&name).copied().unwrap_or(0);
        let depth = if len == 0 {
            0.0
        } else {
            bases as f64 / len as f64
        };
        writeln!(depth_out, "{name}\t{len}\t{bases}\t{depth:.3}\t{count}")?;
    }

    let mut links: HashMap<(String, char, String, char), usize> = HashMap::new();
    for (_read, mut alns) in by_read {
        alns.sort_by(|a, b| {
            a.qstart
                .cmp(&b.qstart)
                .then_with(|| b.alen.cmp(&a.alen))
                .then_with(|| {
                    b.identity()
                        .partial_cmp(&a.identity())
                        .unwrap_or(std::cmp::Ordering::Equal)
                })
        });
        let mut chain: Vec<PafAln> = Vec::new();
        for aln in alns {
            if let Some(last) = chain.last() {
                let overlap = last.qend as isize - aln.qstart as isize;
                if overlap > config.paf_max_link_gap {
                    continue;
                }
            }
            chain.push(aln);
        }
        for pair in chain.windows(2) {
            let left = &pair[0];
            let right = &pair[1];
            if left.tname == right.tname && left.strand == right.strand {
                continue;
            }
            let gap = right.qstart as isize - left.qend as isize;
            if gap.abs() <= config.paf_max_link_gap {
                *links
                    .entry((
                        left.tname.clone(),
                        left.strand,
                        right.tname.clone(),
                        right.strand,
                    ))
                    .or_insert(0) += 1;
            }
        }
    }

    let mut junction_out = File::create(junction_path)?;
    writeln!(
        junction_out,
        "from_unitig\tfrom_orient\tto_unitig\tto_orient\tcount"
    )?;
    let mut rows: Vec<_> = links.into_iter().collect();
    rows.sort_by(|a, b| b.1.cmp(&a.1).then_with(|| a.0.cmp(&b.0)));
    for ((from, from_orient, to, to_orient), count) in rows {
        writeln!(
            junction_out,
            "{from}\t{from_orient}\t{to}\t{to_orient}\t{count}"
        )?;
    }
    Ok(())
}

fn parse_paf_line(line: &str) -> Option<PafAln> {
    let fields: Vec<&str> = line.split('\t').collect();
    if fields.len() < 12 {
        return None;
    }
    Some(PafAln {
        qname: fields[0].to_string(),
        qstart: fields[2].parse().ok()?,
        qend: fields[3].parse().ok()?,
        strand: fields[4].chars().next()?,
        tname: fields[5].to_string(),
        tlen: fields[6].parse().ok()?,
        tstart: fields[7].parse().ok()?,
        tend: fields[8].parse().ok()?,
        matches: fields[9].parse().ok()?,
        alen: fields[10].parse().ok()?,
    })
}

fn read_fasta_lengths(path: &Path) -> io::Result<HashMap<String, usize>> {
    let reader = BufReader::new(File::open(path)?);
    let mut lengths = HashMap::new();
    let mut current = None;
    let mut len = 0usize;
    for line in reader.lines() {
        let line = line?;
        if let Some(header) = line.strip_prefix('>') {
            if let Some(name) = current.take() {
                lengths.insert(name, len);
            }
            current = Some(
                header
                    .split_whitespace()
                    .next()
                    .unwrap_or(header)
                    .to_string(),
            );
            len = 0;
        } else {
            len += line.trim().len();
        }
    }
    if let Some(name) = current {
        lengths.insert(name, len);
    }
    Ok(lengths)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_size_suffixes() {
        assert_eq!(parse_size("500k".to_string()).unwrap(), 500_000);
        assert_eq!(parse_size("1.5m".to_string()).unwrap(), 1_500_000);
    }

    #[test]
    fn recommends_high_profile_subset_from_subset_target_depth() {
        assert_eq!(
            recommend_high_profile_subset(500_000, 1_000_000_000),
            ReadSubset { basis_points: 1500 }
        );
        assert_eq!(
            recommend_high_profile_subset(150_000, 1_000_000_000),
            ReadSubset { basis_points: 450 }
        );
    }

    #[test]
    fn reverse_complements() {
        assert_eq!(revcomp_string("ACGTTA"), "TAACGT");
    }

    #[test]
    fn canonicalizes_kmers() {
        let (canon, forward) = canonical_kmer(b"TTTACG", 0, 6).unwrap();
        assert_eq!(canon, "CGTAAA");
        assert!(!forward);
    }

    #[test]
    fn syncmer_selects_inside_valid_segment() {
        let seq = b"ACGTACGTACGTACGTACGT";
        let hits = select_syncmer_hits(seq, 7, 3, Some(2));
        assert!(hits.iter().all(|&pos| pos + 7 <= seq.len()));
    }

    #[test]
    fn thins_anchor_hits_by_minimum_spacing() {
        assert_eq!(
            thin_anchor_hits(vec![2, 5, 13, 17, 28], 10),
            vec![2, 13, 28]
        );
    }

    #[test]
    fn parses_percent_lists() {
        assert_eq!(
            parse_percent_list("50,12.5,25,12.50", "--subsets").unwrap(),
            vec![
                ReadSubset { basis_points: 1250 },
                ReadSubset { basis_points: 2500 },
                ReadSubset { basis_points: 5000 },
            ]
        );
        assert!(parse_percent_list("0,50", "--read-subsets").is_err());
        assert!(parse_percent_list("25,", "--read-subsets").is_err());
        assert!(parse_percent_list("12.345", "--subsets").is_err());
    }

    #[test]
    fn parses_subsets_alias_arg() {
        let config = Config::from_args(
            [
                "simple_draft_asm",
                "--organelle",
                "plastid",
                "-i",
                "reads.fastq.gz",
                "--subsets=12.5,25",
            ]
            .iter()
            .map(|value| value.to_string())
            .collect(),
        )
        .unwrap();

        assert!(config.read_subsets_requested);
        assert_eq!(
            config.read_subsets,
            vec![
                ReadSubset { basis_points: 1250 },
                ReadSubset { basis_points: 2500 },
            ]
        );
    }

    #[test]
    fn draft_request_mito_stable_uses_two_round_defaults() {
        let config = Config::from_draft_request(DraftAssemblyRequest {
            organelle: DraftOrganelle::Mito,
            data_mode: DraftDataMode::Standard,
            auto_read_subset: false,
            repeat_aware_resolution: true,
            reads: vec![PathBuf::from("reads.fastq.gz")],
            out_dir: PathBuf::from("draft_asm/mito/02.anchor_graph_core"),
            threads: 4,
            min_graph_coverage: None,
            min_branch_ratio: None,
            min_tip_len: None,
            min_link_support: None,
            min_link_ratio: None,
            read_subsets: None,
            keep_debug_files: false,
        })
        .unwrap();

        assert_eq!(config.organelle, Some(OrganelleProfile::Mito));
        assert_eq!(config.data_mode, DataMode::Standard);
        assert_eq!(config.rounds, 2);
        assert!(config.mito_stable);
        assert_eq!(config.min_anchor_coverage, 18);
        assert_eq!(config.min_edge_coverage, 18);
        assert!((config.min_branch_ratio - 0.30).abs() < f64::EPSILON);
        assert_eq!(config.skeleton_min_link_support, 10);
        assert!((config.skeleton_min_link_ratio - 0.20).abs() < f64::EPSILON);
        assert_eq!(config.skeleton_rescue_link_support, 20);
        assert_eq!(config.threads, 4);
    }

    #[test]
    fn draft_request_plastid_standard_uses_one_round_defaults() {
        let config = Config::from_draft_request(DraftAssemblyRequest {
            organelle: DraftOrganelle::Plastid,
            data_mode: DraftDataMode::Standard,
            auto_read_subset: false,
            repeat_aware_resolution: false,
            reads: vec![PathBuf::from("reads.fastq.gz")],
            out_dir: PathBuf::from("draft_asm/plastid/02.anchor_graph_core"),
            threads: 2,
            min_graph_coverage: None,
            min_branch_ratio: None,
            min_tip_len: None,
            min_link_support: None,
            min_link_ratio: None,
            read_subsets: None,
            keep_debug_files: false,
        })
        .unwrap();

        assert_eq!(config.organelle, Some(OrganelleProfile::Plastid));
        assert_eq!(config.data_mode, DataMode::Standard);
        assert_eq!(config.rounds, 1);
        assert!(!config.mito_stable);
        assert_eq!(config.min_anchor_coverage, 18);
        assert_eq!(config.min_edge_coverage, 18);
        assert!((config.min_branch_ratio - 0.30).abs() < f64::EPSILON);
        assert_eq!(config.max_edges_per_state, 3);
        assert_eq!(config.min_tip_len, 3000);
        assert!(config.read_junction_links);
        assert!(config.bidirectional_links);
        assert_eq!(config.threads, 2);
    }

    #[test]
    fn draft_request_auto_subset_uses_standard_data_mode_without_link_ratio_filter() {
        let config = Config::from_draft_request(DraftAssemblyRequest {
            organelle: DraftOrganelle::Plastid,
            data_mode: DraftDataMode::Standard,
            auto_read_subset: true,
            repeat_aware_resolution: false,
            reads: vec![PathBuf::from("reads.fastq.gz")],
            out_dir: PathBuf::from("draft_asm/plastid/02.anchor_graph_core"),
            threads: 2,
            min_graph_coverage: None,
            min_branch_ratio: None,
            min_tip_len: None,
            min_link_support: None,
            min_link_ratio: None,
            read_subsets: None,
            keep_debug_files: false,
        })
        .unwrap();

        assert_eq!(config.organelle, Some(OrganelleProfile::Plastid));
        assert_eq!(config.data_mode, DataMode::Standard);
        assert_eq!(config.rounds, 1);
        assert_eq!(config.min_link_ratio, 0.0);
        assert!(config.read_subsets_requested);
        assert!(config.read_subsets_auto_recommended);
        assert_eq!(
            high_profile_default_genome_size(&config),
            Some(HIGH_PROFILE_PLASTID_DEFAULT_GENOME_SIZE)
        );
    }

    #[test]
    fn parses_mito_stable_three_way_node_constraints() {
        let config = Config::from_args(
            [
                "simple_draft_asm",
                "--organelle",
                "mito",
                "--mito-stable",
                "-i",
                "reads.fastq.gz",
                "--mito-stable-mode",
                "forbid-selected",
                "--selected-nodes=edge_12,edge_17",
            ]
            .iter()
            .map(|value| value.to_string())
            .collect(),
        )
        .unwrap();

        assert_eq!(
            config.mito_stable_selection_mode,
            MitoStableSelectionMode::ForbidSelected
        );
        assert_eq!(
            config.mito_stable_selected_nodes,
            HashSet::from(["edge_12".to_string(), "edge_17".to_string()])
        );
        assert!(config.mito_stable_allow_three_way.is_empty());
        assert_eq!(
            config.mito_stable_forbid_three_way,
            HashSet::from(["edge_12".to_string(), "edge_17".to_string()])
        );
    }

    #[test]
    fn parses_mito_stable_allow_selected_as_analysis_mode() {
        let config = Config::from_args(
            [
                "simple_draft_asm",
                "--organelle",
                "mito",
                "--mito-stable",
                "-i",
                "reads.fastq.gz",
                "--mito-stable-mode=allow-selected",
                "--selected-nodes",
                "edge_4,edge_8",
            ]
            .iter()
            .map(|value| value.to_string())
            .collect(),
        )
        .unwrap();

        assert_eq!(
            config.mito_stable_selection_mode,
            MitoStableSelectionMode::AllowSelected
        );
        assert_eq!(
            config.mito_stable_selected_nodes,
            HashSet::from(["edge_4".to_string(), "edge_8".to_string()])
        );
        assert!(config.mito_stable_allow_three_way.is_empty());
        assert!(config.mito_stable_forbid_three_way.is_empty());
    }

    #[test]
    fn mito_low_uses_low_depth_defaults() {
        let config = Config::from_args(
            [
                "simple_draft_asm",
                "--organelle",
                "mito",
                "--data-mode",
                "low",
                "-i",
                "reads.fasta.gz",
            ]
            .iter()
            .map(|value| value.to_string())
            .collect(),
        )
        .unwrap();

        assert_eq!(config.organelle, Some(OrganelleProfile::Mito));
        assert_eq!(config.data_mode, DataMode::Low);
        assert_eq!(config.rounds, 2);
        assert_eq!(config.skeleton_min_link_support, 10);
    }

    #[test]
    fn explicit_depth_search_subsets_are_parsed_without_link_ratio_filter() {
        let config = Config::from_args(
            [
                "simple_draft_asm",
                "--organelle",
                "mito",
                "--subsets=3,5,10",
                "-i",
                "reads.fastq.gz",
            ]
            .iter()
            .map(|value| value.to_string())
            .collect(),
        )
        .unwrap();

        assert_eq!(config.organelle, Some(OrganelleProfile::Mito));
        assert_eq!(config.data_mode, DataMode::Standard);
        assert_eq!(config.rounds, 2);
        assert_eq!(config.min_link_ratio, 0.0);
        assert!(config.read_subsets_requested);
        assert_eq!(
            config.read_subsets,
            vec![
                ReadSubset { basis_points: 300 },
                ReadSubset { basis_points: 500 },
                ReadSubset { basis_points: 1000 },
            ]
        );
    }

    #[test]
    fn plastid_profiles_keep_direct_anchor_round_logic() {
        let low = Config::from_args(
            [
                "simple_draft_asm",
                "--organelle",
                "plastid",
                "--data-mode",
                "low",
                "-i",
                "reads.fasta.gz",
            ]
            .iter()
            .map(|value| value.to_string())
            .collect(),
        )
        .unwrap();
        let standard = Config::from_args(
            [
                "simple_draft_asm",
                "--organelle",
                "plastid",
                "-i",
                "reads.fastq.gz",
            ]
            .iter()
            .map(|value| value.to_string())
            .collect(),
        )
        .unwrap();
        let high = Config::from_args(
            [
                "simple_draft_asm",
                "--organelle",
                "plastid",
                "--subsets=3,5,10",
                "-i",
                "reads.fastq.gz",
            ]
            .iter()
            .map(|value| value.to_string())
            .collect(),
        )
        .unwrap();

        assert_eq!(low.organelle, Some(OrganelleProfile::Plastid));
        assert_eq!(low.data_mode, DataMode::Low);
        assert_eq!(low.rounds, 1);
        assert_eq!(standard.organelle, Some(OrganelleProfile::Plastid));
        assert_eq!(standard.data_mode, DataMode::Standard);
        assert_eq!(standard.rounds, 1);
        assert_eq!(high.organelle, Some(OrganelleProfile::Plastid));
        assert_eq!(high.data_mode, DataMode::Standard);
        assert_eq!(high.rounds, 1);
        assert_eq!(high.min_link_ratio, 0.0);
        assert!(high.read_subsets_requested);
    }

    #[test]
    fn explicit_subset_and_link_ratio_are_parsed() {
        let config = Config::from_args(
            [
                "simple_draft_asm",
                "--organelle",
                "mito",
                "--subsets=25",
                "--min-link-ratio",
                "0.10",
                "-i",
                "reads.fastq.gz",
            ]
            .iter()
            .map(|value| value.to_string())
            .collect(),
        )
        .unwrap();

        assert_eq!(config.read_subsets, vec![ReadSubset { basis_points: 2500 }]);
        assert!((config.min_link_ratio - 0.10).abs() < f64::EPSILON);
    }

    #[test]
    fn mito_stable_uses_complex_repeat_defaults() {
        let config = Config::from_args(
            ["simple_draft_asm", "--mito-stable", "-i", "reads.fastq.gz"]
                .iter()
                .map(|value| value.to_string())
                .collect(),
        )
        .unwrap();

        assert_eq!(config.organelle, Some(OrganelleProfile::Mito));
        assert_eq!(config.data_mode, DataMode::Standard);
        assert!(config.mito_stable);
        assert_eq!(config.rounds, 2);
        assert_eq!(config.min_anchor_coverage, 18);
        assert_eq!(config.min_edge_coverage, 18);
        assert!((config.min_branch_ratio - 0.30).abs() < f64::EPSILON);
        assert_eq!(config.skeleton_min_link_support, 10);
        assert!((config.skeleton_min_link_ratio - 0.20).abs() < f64::EPSILON);
        assert_eq!(config.skeleton_rescue_link_support, 20);
        assert!(!config.read_subsets_requested);
    }

    #[test]
    fn low_plastid_profile_uses_low_coverage_defaults() {
        let config = Config::from_args(
            [
                "simple_draft_asm",
                "--organelle",
                "plastid",
                "--data-mode",
                "low",
                "-i",
                "reads.fasta.gz",
            ]
            .iter()
            .map(|value| value.to_string())
            .collect(),
        )
        .unwrap();

        assert_eq!(config.data_mode, DataMode::Low);
        assert_eq!(config.rounds, 1);
        assert_eq!(config.min_anchor_coverage, 12);
        assert_eq!(config.min_edge_coverage, 12);
        assert!((config.min_branch_ratio - 0.30).abs() < f64::EPSILON);
        assert_eq!(config.min_tip_len, 3000);
        assert_eq!(config.min_link_support, 20);
        assert!(config.read_junction_links);
    }

    #[test]
    fn low_mito_profile_uses_rescue_friendly_defaults() {
        let config = Config::from_args(
            [
                "simple_draft_asm",
                "--organelle",
                "mito",
                "--data-mode",
                "low",
                "-i",
                "reads.fasta.gz",
            ]
            .iter()
            .map(|value| value.to_string())
            .collect(),
        )
        .unwrap();

        assert_eq!(config.data_mode, DataMode::Low);
        assert_eq!(config.rounds, 2);
        assert_eq!(config.min_anchor_coverage, 12);
        assert_eq!(config.min_edge_coverage, 12);
        assert_eq!(config.min_link_support, 20);
        assert_eq!(config.skeleton_min_link_support, 10);
        assert_eq!(config.skeleton_min_link_ratio, 0.20);
        assert_eq!(config.skeleton_rescue_link_support, 20);
    }

    #[test]
    fn mito_low_depth_detection_is_narrow() {
        let mito_low = Config::from_args(
            [
                "simple_draft_asm",
                "--organelle",
                "mito",
                "--data-mode",
                "low",
                "-i",
                "reads.fasta.gz",
            ]
            .iter()
            .map(|value| value.to_string())
            .collect(),
        )
        .unwrap();
        let mito_standard = Config::from_args(
            [
                "simple_draft_asm",
                "--organelle",
                "mito",
                "-i",
                "reads.fasta.gz",
            ]
            .iter()
            .map(|value| value.to_string())
            .collect(),
        )
        .unwrap();
        let plastid_low = Config::from_args(
            [
                "simple_draft_asm",
                "--organelle",
                "plastid",
                "--data-mode",
                "low",
                "-i",
                "reads.fasta.gz",
            ]
            .iter()
            .map(|value| value.to_string())
            .collect(),
        )
        .unwrap();

        assert!(is_mito_compact(&mito_low));
        assert!(!is_mito_compact(&mito_standard));
        assert!(!is_mito_compact(&plastid_low));
    }

    #[test]
    fn skeleton_graph_stats_identifies_open_roles() {
        let config = Config::from_args(
            [
                "simple_draft_asm",
                "--organelle",
                "mito",
                "--data-mode",
                "low",
                "-i",
                "reads.fasta.gz",
            ]
            .iter()
            .map(|value| value.to_string())
            .collect(),
        )
        .unwrap();
        let segments = vec![
            SkeletonSegment {
                name: "utg0".to_string(),
                sequence: "AAAA".to_string(),
            },
            SkeletonSegment {
                name: "utg1".to_string(),
                sequence: "CCCC".to_string(),
            },
        ];
        let mut links = HashMap::new();
        links.insert(
            SkeletonLinkKey {
                from: "utg0".to_string(),
                from_orient: '+',
                to: "utg1".to_string(),
                to_orient: '+',
            },
            SkeletonLinkSupport {
                skeleton_support: config.min_link_support,
                ..SkeletonLinkSupport::default()
            },
        );

        let stats = skeleton_graph_stats(&config, &segments, &links);
        assert_eq!(stats.kept_links, 1);
        assert_eq!(
            stats.components.iter().map(Vec::len).collect::<Vec<_>>(),
            vec![2]
        );
        assert!(stats.open_out.contains(&("utg1".to_string(), '+')));
        assert!(stats.open_in.contains(&("utg0".to_string(), '+')));
    }

    #[test]
    fn mito_stable_prunes_closed_redundant_chord() {
        let config = Config::from_args(
            [
                "simple_draft_asm",
                "--organelle",
                "mito",
                "--mito-stable",
                "-i",
                "reads.fastq.gz",
            ]
            .iter()
            .map(|value| value.to_string())
            .collect(),
        )
        .unwrap();
        let segments = vec![
            SkeletonSegment {
                name: "a".to_string(),
                sequence: "AAAA".to_string(),
            },
            SkeletonSegment {
                name: "b".to_string(),
                sequence: "CCCC".to_string(),
            },
            SkeletonSegment {
                name: "c".to_string(),
                sequence: "GGGG".to_string(),
            },
        ];
        let mut links = HashMap::new();
        let support = SkeletonLinkSupport {
            skeleton_support: config.min_link_support,
            ..SkeletonLinkSupport::default()
        };
        for key in [
            SkeletonLinkKey {
                from: "a".to_string(),
                from_orient: '+',
                to: "b".to_string(),
                to_orient: '+',
            },
            SkeletonLinkKey {
                from: "b".to_string(),
                from_orient: '+',
                to: "c".to_string(),
                to_orient: '+',
            },
            SkeletonLinkKey {
                from: "c".to_string(),
                from_orient: '+',
                to: "a".to_string(),
                to_orient: '+',
            },
            SkeletonLinkKey {
                from: "a".to_string(),
                from_orient: '+',
                to: "c".to_string(),
                to_orient: '+',
            },
        ] {
            links.insert(key, support.clone());
        }

        let before = mito_stable_topology_stats(&config, &segments, &links);
        assert_eq!(before.open_ends, 0);
        assert_eq!(before.branch_ends, 2);
        assert_eq!(before.cycle_rank, 2);

        let pruned = prune_mito_stable_redundant_links(&config, &segments, &mut links);
        let after = mito_stable_topology_stats(&config, &segments, &links);
        assert_eq!(pruned.len(), 1);
        assert_eq!(pruned[0].from, "a");
        assert_eq!(pruned[0].to, "c");
        assert_eq!(after.open_ends, 0);
        assert_eq!(after.branch_ends, 0);
        assert_eq!(after.cycle_rank, 1);
    }

    #[test]
    fn mito_stable_bridge_selection_respects_explicit_forbid_three_way() {
        let config = Config::from_args(
            [
                "simple_draft_asm",
                "--organelle",
                "mito",
                "--mito-stable",
                "-i",
                "reads.fastq.gz",
            ]
            .iter()
            .map(|value| value.to_string())
            .collect(),
        )
        .unwrap();
        let segments = ["a", "b", "c"]
            .iter()
            .map(|name| SkeletonSegment {
                name: (*name).to_string(),
                sequence: "AAAA".to_string(),
            })
            .collect::<Vec<_>>();
        let support = SkeletonLinkSupport {
            skeleton_support: config.min_link_support,
            ..SkeletonLinkSupport::default()
        };
        let mut base_links = HashMap::new();
        for (from, to) in [("a", "b"), ("b", "a")] {
            insert_mito_stable_supported_link_pair(
                &config,
                &mut base_links,
                SkeletonLinkKey {
                    from: from.to_string(),
                    from_orient: '+',
                    to: to.to_string(),
                    to_orient: '+',
                },
                support.clone(),
            );
        }
        let candidate = MitoCompactBridgeCandidate {
            key: SkeletonLinkKey {
                from: "c".to_string(),
                from_orient: '+',
                to: "b".to_string(),
                to_orient: '+',
            },
            support: config.skeleton_min_link_support + 1,
            source: MitoCompactBridgeSource::LocalPaf,
        };

        let selected = select_mito_stable_bridges(
            &config,
            &segments,
            &base_links,
            &[candidate.clone()],
            false,
        );
        assert_eq!(selected.len(), 1);

        let mut forbid_config = config.clone();
        forbid_config
            .mito_stable_forbid_three_way
            .insert("b".to_string());
        let selected_forbidden =
            select_mito_stable_bridges(&forbid_config, &segments, &base_links, &[candidate], true);
        assert!(selected_forbidden.is_empty());
    }

    #[test]
    fn mito_stable_scan_detects_three_way_merge_candidates() {
        let config = Config::from_args(
            [
                "simple_draft_asm",
                "--organelle",
                "mito",
                "--mito-stable",
                "-i",
                "reads.fastq.gz",
            ]
            .iter()
            .map(|value| value.to_string())
            .collect(),
        )
        .unwrap();
        let segments = ["t1", "t2", "bridge", "x", "y", "z", "w", "p", "q"]
            .iter()
            .map(|name| SkeletonSegment {
                name: (*name).to_string(),
                sequence: "AAAA".to_string(),
            })
            .collect::<Vec<_>>();
        let mut links = HashMap::new();
        let support = SkeletonLinkSupport {
            skeleton_support: config.min_link_support,
            ..SkeletonLinkSupport::default()
        };
        let add = |links: &mut HashMap<SkeletonLinkKey, SkeletonLinkSupport>,
                   from: &str,
                   from_orient: char,
                   to: &str,
                   to_orient: char| {
            insert_mito_stable_supported_link_pair(
                &config,
                links,
                SkeletonLinkKey {
                    from: from.to_string(),
                    from_orient,
                    to: to.to_string(),
                    to_orient,
                },
                support.clone(),
            );
        };
        for (from, from_orient, to, to_orient) in [
            ("t1", '+', "bridge", '+'),
            ("x", '+', "t1", '+'),
            ("y", '+', "t1", '+'),
            ("t2", '+', "bridge", '-'),
            ("z", '+', "t2", '+'),
            ("w", '+', "t2", '+'),
            ("p", '+', "bridge", '+'),
            ("q", '+', "bridge", '-'),
        ] {
            add(&mut links, from, from_orient, to, to_orient);
        }

        let degrees = mito_stable_physical_side_degrees(&config, &segments, &links);
        assert_eq!(
            mito_stable_node_degree_class(degrees["bridge"].0, degrees["bridge"].1),
            "2-2"
        );
        assert_eq!(
            mito_stable_node_degree_class(degrees["t1"].0, degrees["t1"].1),
            "three_way"
        );
        assert_eq!(
            mito_stable_node_degree_class(degrees["t2"].0, degrees["t2"].1),
            "three_way"
        );

        let bridge_evidence = mito_stable_small_repeat_bridge_evidence(&config, &segments, &links);
        assert_eq!(
            bridge_evidence.get("bridge").cloned().unwrap_or_default(),
            vec!["t1".to_string(), "t2".to_string()]
        );
        let merge_evidence = mito_stable_three_way_merge_evidence(&config, &segments, &links);
        assert!(merge_evidence
            .get("t1")
            .is_some_and(|evidence| evidence.iter().any(|entry| entry.contains("bridge=bridge"))));
        assert!(merge_evidence
            .get("t2")
            .is_some_and(|evidence| evidence.iter().any(|entry| entry.contains("bridge=bridge"))));
    }

    #[test]
    fn mito_stable_expands_low_overlap_repeat_shortcut_when_node_class_is_stable() {
        let config = Config::from_args(
            [
                "simple_draft_asm",
                "--organelle",
                "mito",
                "--mito-stable",
                "-i",
                "reads.fastq.gz",
            ]
            .iter()
            .map(|value| value.to_string())
            .collect(),
        )
        .unwrap();
        let mut segments = vec![
            SkeletonSegment {
                name: "a".to_string(),
                sequence: format!("{}{}", "G".repeat(20), "A".repeat(120)),
            },
            SkeletonSegment {
                name: "b".to_string(),
                sequence: format!("{}{}", "A".repeat(120), "C".repeat(120)),
            },
            SkeletonSegment {
                name: "c".to_string(),
                sequence: format!("{}{}", "C".repeat(120), "T".repeat(20)),
            },
            SkeletonSegment {
                name: "d".to_string(),
                sequence: format!("{}{}", "C".repeat(120), "G".repeat(20)),
            },
        ];
        let key_ab = SkeletonLinkKey {
            from: "a".to_string(),
            from_orient: '+',
            to: "b".to_string(),
            to_orient: '+',
        };
        let key_ac = SkeletonLinkKey {
            from: "a".to_string(),
            from_orient: '+',
            to: "c".to_string(),
            to_orient: '+',
        };
        let key_bc = SkeletonLinkKey {
            from: "b".to_string(),
            from_orient: '+',
            to: "c".to_string(),
            to_orient: '+',
        };
        let key_db = SkeletonLinkKey {
            from: "d".to_string(),
            from_orient: '+',
            to: "b".to_string(),
            to_orient: '+',
        };
        let key_da = SkeletonLinkKey {
            from: "d".to_string(),
            from_orient: '-',
            to: "a".to_string(),
            to_orient: '-',
        };
        let mut links = HashMap::new();
        insert_mito_stable_supported_link_pair(
            &config,
            &mut links,
            key_ab.clone(),
            SkeletonLinkSupport {
                skeleton_support: config.min_link_support,
                ..SkeletonLinkSupport::default()
            },
        );
        insert_mito_stable_supported_link_pair(
            &config,
            &mut links,
            key_db.clone(),
            SkeletonLinkSupport {
                skeleton_support: config.min_link_support,
                ..SkeletonLinkSupport::default()
            },
        );
        insert_mito_stable_supported_link_pair(
            &config,
            &mut links,
            key_da,
            SkeletonLinkSupport {
                skeleton_support: config.min_link_support,
                ..SkeletonLinkSupport::default()
            },
        );
        insert_mito_stable_supported_link_pair(
            &config,
            &mut links,
            key_ac.clone(),
            SkeletonLinkSupport {
                paf_support: config.skeleton_min_link_support + 1,
                ..SkeletonLinkSupport::default()
            },
        );
        refresh_skeleton_link_ratios(&config, &mut links);
        let candidates = vec![MitoCompactBridgeCandidate {
            key: key_bc.clone(),
            support: config.skeleton_min_link_support + 5,
            source: MitoCompactBridgeSource::LocalPaf,
        }];
        let mut depth = HashMap::new();

        let expansions = expand_mito_stable_repeat_shortcuts(
            &config,
            &mut segments,
            &mut depth,
            &mut links,
            &candidates,
        );

        assert_eq!(expansions.len(), 1);
        assert_eq!(expansions[0].template, "b");
        assert_eq!(expansions[0].shortcut_overlap, 0);
        assert_eq!(expansions[0].left_overlap, 120);
        assert_eq!(expansions[0].right_overlap, 120);
        assert!(segments.iter().any(|segment| segment.name == "b_copy0"));
        assert!(!links.contains_key(&key_ac));
        assert!(!links.contains_key(&key_ab));
        assert!(skeleton_link_is_kept(
            &config,
            &links,
            &SkeletonLinkKey {
                from: "a".to_string(),
                from_orient: '+',
                to: "b_copy0".to_string(),
                to_orient: '+',
            },
        ));
        assert!(skeleton_link_is_kept(
            &config,
            &links,
            &SkeletonLinkKey {
                from: "b_copy0".to_string(),
                from_orient: '+',
                to: "c".to_string(),
                to_orient: '+',
            },
        ));
        assert!(skeleton_link_is_kept(&config, &links, &key_db));
    }

    #[test]
    fn mito_stable_detects_new_three_way_nodes() {
        let config = Config::from_args(
            [
                "simple_draft_asm",
                "--organelle",
                "mito",
                "--mito-stable",
                "-i",
                "reads.fastq.gz",
            ]
            .iter()
            .map(|value| value.to_string())
            .collect(),
        )
        .unwrap();
        let before_segments = ["a", "b", "c", "d", "e", "x"]
            .iter()
            .map(|name| SkeletonSegment {
                name: (*name).to_string(),
                sequence: "AAAA".to_string(),
            })
            .collect::<Vec<_>>();
        let mut after_segments = before_segments.clone();
        after_segments.push(SkeletonSegment {
            name: "b_copy0".to_string(),
            sequence: "AAAA".to_string(),
        });

        let mut before_links = HashMap::new();
        let mut after_links = HashMap::new();
        let support = SkeletonLinkSupport {
            skeleton_support: config.min_link_support,
            ..SkeletonLinkSupport::default()
        };
        let add = |links: &mut HashMap<SkeletonLinkKey, SkeletonLinkSupport>,
                   from: &str,
                   from_orient: char,
                   to: &str,
                   to_orient: char| {
            insert_mito_stable_supported_link_pair(
                &config,
                links,
                SkeletonLinkKey {
                    from: from.to_string(),
                    from_orient,
                    to: to.to_string(),
                    to_orient,
                },
                support.clone(),
            );
        };

        for (from, from_orient, to, to_orient) in [
            ("a", '+', "b", '+'),
            ("x", '+', "b", '+'),
            ("b", '+', "d", '+'),
            ("b", '+', "e", '+'),
            ("a", '+', "c", '+'),
        ] {
            add(&mut before_links, from, from_orient, to, to_orient);
        }
        for (from, from_orient, to, to_orient) in [
            ("x", '+', "b", '+'),
            ("b", '+', "d", '+'),
            ("b", '+', "e", '+'),
            ("a", '+', "b_copy0", '+'),
            ("b_copy0", '+', "c", '+'),
        ] {
            add(&mut after_links, from, from_orient, to, to_orient);
        }

        assert!(
            !mito_stable_three_way_nodes(&config, &before_segments, &before_links).contains("b")
        );
        let after_three_way = mito_stable_three_way_nodes(&config, &after_segments, &after_links);
        assert_eq!(after_three_way, HashSet::from(["b".to_string()]));
        assert!(mito_stable_introduces_new_three_way_nodes(
            &config,
            &before_segments,
            &before_links,
            &after_segments,
            &after_links
        ));
        assert!(!mito_stable_node_constraints_ok_for_transition(
            &config,
            &before_segments,
            &before_links,
            &after_segments,
            &after_links
        ));

        let mut allow_config = config.clone();
        allow_config.mito_stable_allow_three_way = after_three_way;
        assert!(mito_stable_node_constraints_ok_for_transition(
            &allow_config,
            &before_segments,
            &before_links,
            &after_segments,
            &after_links
        ));

        let mut forbid_config = config.clone();
        forbid_config
            .mito_stable_forbid_three_way
            .insert("b".to_string());
        assert!(!mito_stable_node_constraints_ok(
            &forbid_config,
            &after_segments,
            &after_links
        ));
        assert!(!mito_stable_node_constraints_ok_for_transition(
            &forbid_config,
            &before_segments,
            &before_links,
            &after_segments,
            &after_links
        ));
    }

    #[test]
    fn read_subsets_are_nested() {
        let subset_25 = ReadSubset { basis_points: 2500 };
        let subset_50 = ReadSubset { basis_points: 5000 };
        let subset_75 = ReadSubset { basis_points: 7500 };
        for i in 0..1000 {
            let name = format!("read-{i}");
            let bucket = read_subset_bucket(&name, i);
            if subset_25.includes_bucket(bucket) {
                assert!(subset_50.includes_bucket(bucket));
            }
            if subset_50.includes_bucket(bucket) {
                assert!(subset_75.includes_bucket(bucket));
            }
        }
    }
}
