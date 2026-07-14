use std::cmp::Ordering;
use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet};
use std::fmt::Write as FmtWrite;
use std::fs::{self, File};
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};

use crate::commands::polish::{evaluate_sv_candidate, SvCandidateEvaluation};
use crate::error::OrgraftError;
use crate::sv_graph::{
    localize_sv_to_unitig_graph, write_unavailable_localization, SvGraphLocalization,
    SvGraphLocalizationRequest,
};

const MAX_GENERATED_CANDIDATES: usize = 100_000;
const DEFAULT_CANDIDATE_EVALUATIONS: usize = 24;
const MIN_TARGET_TYPE1_FRACTION: f64 = 0.80;
const MAX_GLOBAL_SUPPORT_DROP: f64 = 0.005;
const IMPROVEMENT_EPSILON: f64 = 1e-6;
const REPEAT_PAIRING_MIN_OVERLAP_BP: f64 = 500.0;
const REPEAT_PAIRING_MIN_PATH_READS: usize = 10;
const REPEAT_PAIRING_MIN_DOMINANCE_RATIO: f64 = 3.0;
const REPEAT_PAIRING_FLANK_BP: usize = 1_000;

#[derive(Debug, Clone, Copy)]
pub(crate) struct SvMetrics {
    pub reference_support_read_fraction: f64,
    pub reference_support_depth_area_fraction: f64,
    pub low_green_window_fraction: f64,
}

#[derive(Debug)]
pub(crate) struct SvRepairRequest<'a> {
    pub reference: &'a Path,
    pub reads: &'a Path,
    pub soft_paths: &'a Path,
    pub summary: &'a Path,
    pub high_subgroups: &'a Path,
    pub segments: &'a Path,
    pub read_index: &'a Path,
    pub unitig_graph: Option<&'a Path>,
    pub output_dir: &'a Path,
    pub manual_subgroup: Option<&'a str>,
    pub threads: usize,
    pub prior_fastas: &'a [PathBuf],
    pub max_candidate_evaluations: Option<usize>,
}

#[derive(Debug, Clone)]
pub(crate) struct SvRepairResult {
    pub subgroup_spec: String,
    pub corrected_fasta: PathBuf,
    pub score_table: PathBuf,
    pub repair_report: PathBuf,
    pub candidate_count: usize,
    pub evaluated_candidates: usize,
    pub target_reads: usize,
    pub target_type1_reads: usize,
    pub graph_localization: SvGraphLocalization,
    pub evaluation: SvCandidateEvaluation,
}

#[derive(Debug, Clone)]
pub(crate) struct RepeatPairingMismatch {
    pub repeat_node: String,
    pub repeat_length: usize,
    pub representative_group: String,
    pub representative_old_index: usize,
    pub subgroup_specs: Vec<String>,
    pub alternate_path_reads: [usize; 2],
    pub reference_spanning_reads: [usize; 2],
    pub repeat_intervals: [(usize, usize); 2],
    pub dominance_ratio: f64,
    pub target_read_ids: HashSet<String>,
}

#[derive(Debug, Clone)]
struct RepeatPairingSubgroupRow {
    old_index: usize,
    start: usize,
    end: usize,
    reads: usize,
    overlap: f64,
}

#[derive(Debug, Clone)]
struct PhysicalRepeatPath {
    endpoints: (usize, usize),
    overlap: f64,
    rows: Vec<RepeatPairingSubgroupRow>,
}

#[derive(Debug, Clone)]
struct HighSubgroup {
    group_name: String,
    old_index: usize,
    boundary_key: String,
    subgroup_reads: usize,
    is_reference_support: bool,
    auto_highlight: bool,
    min_reference_fraction: f64,
    judgement: String,
}

impl HighSubgroup {
    fn spec(&self) -> String {
        format!("{}:{}", self.group_name, self.old_index)
    }

    fn is_automatic_repair_candidate(&self) -> bool {
        !self.is_reference_support
            && self.auto_highlight
            && self.judgement == "possible_reference_sv_error"
    }
}

#[derive(Debug, Clone)]
struct Segment {
    read_id: String,
    read_class: String,
    group_name: String,
    subgroup_old_index: Option<usize>,
    segment_index: usize,
    segment_count: usize,
    query_start: usize,
    query_end: usize,
    subject_start: usize,
    subject_end: usize,
    strand: char,
}

#[derive(Debug, Clone)]
struct Block {
    id: usize,
    start: usize,
    end: usize,
}

impl Block {
    fn len(&self) -> usize {
        self.end - self.start + 1
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
struct BlockToken {
    id: usize,
    orient: char,
}

#[derive(Debug, Clone)]
struct CandidateDescriptor {
    tokens: Vec<BlockToken>,
    predicted_fl_reads: usize,
    predicted_target_reads: usize,
    original_adjacencies: usize,
    reversed_bases: usize,
}

#[derive(Debug)]
struct CandidateScore {
    rank: usize,
    descriptor: CandidateDescriptor,
    target_type1_reads: usize,
    target_reads: usize,
    evaluation: SvCandidateEvaluation,
    candidate_path: PathBuf,
    evaluation_dir: PathBuf,
}

impl CandidateScore {
    fn target_type1_fraction(&self) -> f64 {
        if self.target_reads == 0 {
            0.0
        } else {
            self.target_type1_reads as f64 / self.target_reads as f64
        }
    }
}

pub(crate) fn read_sv_metrics(path: &Path) -> Result<SvMetrics, OrgraftError> {
    Ok(SvMetrics {
        reference_support_read_fraction: read_metric_f64(path, "reference_support_read_fraction")?,
        reference_support_depth_area_fraction: read_metric_f64(
            path,
            "reference_support_depth_area_fraction",
        )?,
        low_green_window_fraction: read_metric_f64(path, "low_green_window_fraction")?,
    })
}

pub(crate) fn select_sv_subgroup_spec(
    high_subgroups: &Path,
    manual_subgroup: Option<&str>,
) -> Result<Option<String>, OrgraftError> {
    let rows = read_high_subgroups(high_subgroups)?;
    Ok(select_subgroup(&rows, manual_subgroup)?.map(|row| row.spec()))
}

pub(crate) fn detect_repeat_pairing_mismatch(
    subgroup_stats: &Path,
    read_index: &Path,
    read_evidence: &Path,
    unitig_graph: Option<&Path>,
) -> Result<Option<RepeatPairingMismatch>, OrgraftError> {
    if !subgroup_stats.is_file() || !read_index.is_file() || !read_evidence.is_file() {
        return Ok(None);
    }

    let rows = read_repeat_pairing_subgroups(subgroup_stats)?;
    let physical_paths = collapse_reciprocal_repeat_paths(&rows);
    let read_rows = read_index_rows(read_index)?;
    let evidence = read_primary_evidence(read_evidence)?;
    let baseline_type1 = read_rows
        .iter()
        .filter(|row| row.read_class == "FL" && row.group_name == "type_1_subtype_NA")
        .map(|row| row.read_id.clone())
        .collect::<HashSet<_>>();
    let max_reference_position = evidence.iter().map(|row| row.target_end).max().unwrap_or(0);
    let mut mismatches = Vec::new();

    for left_index in 0..physical_paths.len() {
        for right_index in left_index + 1..physical_paths.len() {
            let left = &physical_paths[left_index];
            let right = &physical_paths[right_index];
            let repeat_length = ((left.overlap + right.overlap) / 2.0).round() as usize;
            let tolerance = ((repeat_length as f64 * 0.15).round() as usize).max(50);
            let first_delta = left.endpoints.0.abs_diff(right.endpoints.0);
            let second_delta = left.endpoints.1.abs_diff(right.endpoints.1);
            if first_delta.abs_diff(repeat_length) > tolerance
                || second_delta.abs_diff(repeat_length) > tolerance
            {
                continue;
            }

            let intervals = [
                ordered_pair(left.endpoints.0, right.endpoints.0),
                ordered_pair(left.endpoints.1, right.endpoints.1),
            ];
            if intervals.iter().any(|(start, end)| {
                *start <= REPEAT_PAIRING_FLANK_BP
                    || end.saturating_add(REPEAT_PAIRING_FLANK_BP) > max_reference_position
            }) {
                continue;
            }

            let left_specs = path_specs(left);
            let right_specs = path_specs(right);
            let left_ids = ids_for_subgroup_specs(&read_rows, &left_specs);
            let right_ids = ids_for_subgroup_specs(&read_rows, &right_specs);
            let alternate_path_reads = [left_ids.len(), right_ids.len()];
            if alternate_path_reads
                .iter()
                .any(|count| *count < REPEAT_PAIRING_MIN_PATH_READS)
            {
                continue;
            }

            let reference_spanning_reads = intervals
                .map(|interval| spanning_type1_reads(&evidence, &baseline_type1, interval).len());
            let alternate_total = alternate_path_reads.iter().sum::<usize>();
            let reference_total = reference_spanning_reads.iter().sum::<usize>();
            let dominance_ratio = alternate_total as f64 / reference_total.max(1) as f64;
            if dominance_ratio < REPEAT_PAIRING_MIN_DOMINANCE_RATIO {
                continue;
            }

            let mut all_rows = left.rows.iter().chain(&right.rows).collect::<Vec<_>>();
            all_rows.sort_by(|a, b| {
                b.reads
                    .cmp(&a.reads)
                    .then_with(|| a.old_index.cmp(&b.old_index))
            });
            let representative = all_rows[0];
            let mut subgroup_specs = left_specs;
            subgroup_specs.extend(right_specs);
            subgroup_specs.sort();
            subgroup_specs.dedup();
            let mut target_read_ids = left_ids;
            target_read_ids.extend(right_ids);
            let Some(unitig_graph) = unitig_graph.filter(|path| path.is_file()) else {
                continue;
            };
            let Some(repeat_node) = unique_repeat_node_for_length(unitig_graph, repeat_length)?
            else {
                continue;
            };

            mismatches.push(RepeatPairingMismatch {
                repeat_node,
                repeat_length,
                representative_group: "type_2_subtype_rep_NA".to_string(),
                representative_old_index: representative.old_index,
                subgroup_specs,
                alternate_path_reads,
                reference_spanning_reads,
                repeat_intervals: intervals,
                dominance_ratio,
                target_read_ids,
            });
        }
    }

    mismatches.sort_by(|left, right| {
        right
            .dominance_ratio
            .partial_cmp(&left.dominance_ratio)
            .unwrap_or(Ordering::Equal)
            .then_with(|| right.target_read_ids.len().cmp(&left.target_read_ids.len()))
            .then_with(|| left.repeat_node.cmp(&right.repeat_node))
    });
    Ok(mismatches.into_iter().next())
}

pub(crate) fn write_repeat_pairing_mismatch_report(
    path: &Path,
    mismatch: &RepeatPairingMismatch,
) -> Result<(), OrgraftError> {
    let mut output = String::from("metric\tvalue\n");
    metric(&mut output, "status", "repeat_pairing_mismatch");
    metric(&mut output, "repeat_node", &mismatch.repeat_node);
    metric(
        &mut output,
        "repeat_length",
        &mismatch.repeat_length.to_string(),
    );
    metric(
        &mut output,
        "repeat_intervals",
        &mismatch
            .repeat_intervals
            .iter()
            .map(|(start, end)| format!("{start}-{end}"))
            .collect::<Vec<_>>()
            .join(","),
    );
    metric(
        &mut output,
        "alternate_path_reads",
        &mismatch
            .alternate_path_reads
            .iter()
            .map(ToString::to_string)
            .collect::<Vec<_>>()
            .join(","),
    );
    metric(
        &mut output,
        "reference_spanning_reads",
        &mismatch
            .reference_spanning_reads
            .iter()
            .map(ToString::to_string)
            .collect::<Vec<_>>()
            .join(","),
    );
    metric(
        &mut output,
        "dominance_ratio",
        &format!("{:.6}", mismatch.dominance_ratio),
    );
    metric(&mut output, "subgroups", &mismatch.subgroup_specs.join(","));
    metric(
        &mut output,
        "representative_subgroup",
        &format!(
            "{}:{}",
            mismatch.representative_group, mismatch.representative_old_index
        ),
    );
    metric(
        &mut output,
        "target_reads",
        &mismatch.target_read_ids.len().to_string(),
    );
    fs::write(path, output)?;
    Ok(())
}

pub(crate) fn repair_sv_subgroup(
    request: &SvRepairRequest<'_>,
) -> Result<Option<SvRepairResult>, OrgraftError> {
    let baseline = read_sv_metrics(request.summary)?;
    let high_rows = read_high_subgroups(request.high_subgroups)?;
    let Some(selected) = select_subgroup(&high_rows, request.manual_subgroup)? else {
        return Ok(None);
    };

    if request.output_dir.exists() {
        fs::remove_dir_all(request.output_dir)?;
    }
    fs::create_dir_all(request.output_dir)?;
    let graph_localization = localize_selected_subgroup(request, &selected)?;

    let segments = read_segments(request.segments)?;
    let target_ids = read_subgroup_ids(request.read_index, &selected)?;
    if target_ids.is_empty() {
        return Err(OrgraftError::InvalidArgument(format!(
            "{} contains no reads for SV subgroup `{}`",
            request.read_index.display(),
            selected.spec()
        )));
    }

    let (header, reference) = read_single_fasta(request.reference)?;
    let generated = generate_candidates(&reference, &segments, &selected, &target_ids)?;
    if generated.is_empty() {
        return Err(OrgraftError::InvalidArgument(format!(
            "SV subgroup `{}` did not produce a structurally valid candidate",
            selected.spec()
        )));
    }

    let max_evaluations = request
        .max_candidate_evaluations
        .unwrap_or(DEFAULT_CANDIDATE_EVALUATIONS)
        .clamp(1, generated.len());
    let candidates_dir = request.output_dir.join("candidates");
    fs::create_dir_all(&candidates_dir)?;
    let score_table = request.output_dir.join("sv_candidate_scores.tsv");
    let mut score_text = String::from(
        "candidate\tpredicted_fl_reads\tpredicted_target_reads\toriginal_adjacencies\treversed_bases\ttarget_type1_reads\ttarget_reads\ttarget_type1_fraction\tsv_status\tfl_reads\treference_support_reads\treference_support_read_fraction\treference_support_depth_area_fraction\tlow_green_window_fraction\taccepted\n",
    );
    let mut best: Option<CandidateScore> = None;
    let mut evaluated_candidates = 0usize;
    let prior_sequences = read_prior_sequences(request.prior_fastas)?;

    for (index, descriptor) in generated.iter().take(max_evaluations).enumerate() {
        let rank = index + 1;
        let sequence = materialize_candidate(&reference, descriptor, &selected, &segments)?;
        if is_historical_sequence(&sequence, &prior_sequences) {
            append_skipped_score(&mut score_text, rank, descriptor, "historical_cycle");
            continue;
        }

        let candidate_path = candidates_dir.join(format!("candidate_{rank:03}.fasta"));
        write_fasta(&candidate_path, &header, &sequence)?;
        let evaluation_dir = candidates_dir.join(format!("candidate_{rank:03}.sv_eval"));
        let evaluation = evaluate_sv_candidate(
            &candidate_path,
            request.reads,
            request.soft_paths,
            &evaluation_dir,
            request.threads,
        )?;
        evaluated_candidates += 1;
        let target_type1_reads = count_type1_reads(&evaluation.read_index_path, &target_ids)?;
        let score = CandidateScore {
            rank,
            descriptor: descriptor.clone(),
            target_type1_reads,
            target_reads: target_ids.len(),
            evaluation,
            candidate_path,
            evaluation_dir,
        };
        let accepted = candidate_improves(&score, baseline);
        append_candidate_score(&mut score_text, &score, accepted);

        if accepted
            && best
                .as_ref()
                .is_none_or(|current| candidate_better(&score, current))
        {
            if let Some(previous) = best.take() {
                remove_candidate_artifacts(&previous)?;
            }
            best = Some(score);
        } else {
            remove_candidate_artifacts(&score)?;
        }
    }
    fs::write(&score_table, score_text)?;

    let Some(best) = best else {
        let report = request.output_dir.join("sv_repair.tsv");
        write_no_candidate_report(
            &report,
            &selected,
            baseline,
            generated.len(),
            evaluated_candidates,
            &score_table,
            &graph_localization,
        )?;
        return Ok(None);
    };

    let corrected_fasta = request.output_dir.join("sv_corrected.fasta");
    fs::rename(&best.candidate_path, &corrected_fasta)?;
    let selected_eval_dir = request.output_dir.join("selected_sv_eval");
    fs::rename(&best.evaluation_dir, &selected_eval_dir)?;
    let _ = fs::remove_dir(&candidates_dir);

    let repair_report = request.output_dir.join("sv_repair.tsv");
    write_repair_report(
        &repair_report,
        &selected,
        baseline,
        &best,
        generated.len(),
        evaluated_candidates,
        &corrected_fasta,
        &score_table,
        &graph_localization,
    )?;

    let mut evaluation = best.evaluation;
    evaluation.read_index_path = selected_eval_dir
        .join("sv_candidate/subgraph_001/round_2/03.validate/01.data/sv_read_index.tsv");
    evaluation.high_subgroups_path = selected_eval_dir
        .join("sv_candidate/subgraph_001/round_2/03.validate/03.reports/sv_high_subgroups.tsv");
    evaluation.summary_path = selected_eval_dir
        .join("sv_candidate/subgraph_001/round_2/03.validate/03.reports/sv_snv_indel_summary.tsv");

    Ok(Some(SvRepairResult {
        subgroup_spec: selected.spec(),
        corrected_fasta,
        score_table,
        repair_report,
        candidate_count: generated.len(),
        evaluated_candidates,
        target_reads: best.target_reads,
        target_type1_reads: best.target_type1_reads,
        graph_localization,
        evaluation,
    }))
}

fn localize_selected_subgroup(
    request: &SvRepairRequest<'_>,
    selected: &HighSubgroup,
) -> Result<SvGraphLocalization, OrgraftError> {
    let Some(unitig_graph) = request.unitig_graph else {
        return write_unavailable_localization(
            request.output_dir,
            &selected.spec(),
            None,
            "02.unitig_graph/graph.gfa was not found",
        );
    };
    if !unitig_graph.exists() {
        return write_unavailable_localization(
            request.output_dir,
            &selected.spec(),
            Some(unitig_graph),
            "configured unitig graph does not exist",
        );
    }
    match localize_sv_to_unitig_graph(&SvGraphLocalizationRequest {
        subgroup: &selected.spec(),
        boundary_key: &selected.boundary_key,
        reference: request.reference,
        unitig_gfa: unitig_graph,
        soft_paths: request.soft_paths,
        output_dir: request.output_dir,
        threads: request.threads,
    }) {
        Ok(localization) => Ok(localization),
        Err(error) => write_unavailable_localization(
            request.output_dir,
            &selected.spec(),
            Some(unitig_graph),
            &error.to_string(),
        ),
    }
}

fn read_high_subgroups(path: &Path) -> Result<Vec<HighSubgroup>, OrgraftError> {
    let file = File::open(path).map_err(|error| {
        OrgraftError::InvalidArgument(format!(
            "could not read SV subgroup report {}: {error}",
            path.display()
        ))
    })?;
    let mut lines = BufReader::new(file).lines();
    let header = lines
        .next()
        .transpose()?
        .ok_or_else(|| OrgraftError::InvalidArgument(format!("{} is empty", path.display())))?;
    let columns = header_columns(&header);
    let mut rows = Vec::new();
    for line_result in lines {
        let line = line_result?;
        if line.trim().is_empty() {
            continue;
        }
        let fields = line.split('\t').collect::<Vec<_>>();
        rows.push(HighSubgroup {
            group_name: field(&fields, &columns, "group_name", path)?.to_string(),
            old_index: parse_usize_field(&fields, &columns, "old_index", path)?,
            boundary_key: field(&fields, &columns, "boundary_key", path)?.to_string(),
            subgroup_reads: parse_usize_field(&fields, &columns, "subgroup_reads", path)?,
            is_reference_support: parse_bool_field(
                &fields,
                &columns,
                "is_reference_support_subgroup",
                path,
            )?,
            auto_highlight: parse_bool_field(&fields, &columns, "auto_highlight_default", path)?,
            min_reference_fraction: parse_f64_field(
                &fields,
                &columns,
                "min_window_reference_support_fraction",
                path,
            )?,
            judgement: field(&fields, &columns, "judgement", path)?.to_string(),
        });
    }
    Ok(rows)
}

fn select_subgroup(
    rows: &[HighSubgroup],
    manual: Option<&str>,
) -> Result<Option<HighSubgroup>, OrgraftError> {
    if let Some(spec) = manual {
        let (requested_group, requested_index) = parse_subgroup_spec(spec)?;
        let matching = rows.iter().find(|row| {
            row.old_index == requested_index
                && (row.group_name == requested_group
                    || normalize_group_name(&row.group_name)
                        == normalize_group_name(&requested_group))
        });
        return Ok(matching.cloned());
    }

    let mut candidates = rows
        .iter()
        .filter(|row| row.is_automatic_repair_candidate())
        .cloned()
        .collect::<Vec<_>>();
    candidates.sort_by(|left, right| {
        left.min_reference_fraction
            .partial_cmp(&right.min_reference_fraction)
            .unwrap_or(Ordering::Equal)
            .then_with(|| right.subgroup_reads.cmp(&left.subgroup_reads))
            .then_with(|| left.group_name.cmp(&right.group_name))
            .then_with(|| left.old_index.cmp(&right.old_index))
    });
    Ok(candidates.into_iter().next())
}

fn parse_subgroup_spec(spec: &str) -> Result<(String, usize), OrgraftError> {
    let (group, index) = spec.rsplit_once(':').ok_or_else(|| {
        OrgraftError::InvalidArgument(format!(
            "SV subgroup `{spec}` must use group_name:old_index"
        ))
    })?;
    let index = index.parse::<usize>().map_err(|error| {
        OrgraftError::InvalidArgument(format!(
            "SV subgroup `{spec}` has invalid old_index: {error}"
        ))
    })?;
    Ok((group.to_string(), index))
}

fn normalize_group_name(value: &str) -> String {
    value
        .replace("_subtype_", "_")
        .trim_end_matches("_NA")
        .to_string()
}

#[derive(Debug, Clone)]
struct RepeatReadIndexRow {
    read_id: String,
    read_class: String,
    group_name: String,
    old_index: Option<usize>,
}

#[derive(Debug, Clone)]
struct PrimaryEvidenceRow {
    read_id: String,
    target_start: usize,
    target_end: usize,
}

fn read_repeat_pairing_subgroups(
    path: &Path,
) -> Result<Vec<RepeatPairingSubgroupRow>, OrgraftError> {
    let file = File::open(path)?;
    let mut lines = BufReader::new(file).lines();
    let header = lines
        .next()
        .transpose()?
        .ok_or_else(|| OrgraftError::InvalidArgument(format!("{} is empty", path.display())))?;
    let columns = header_columns(&header);
    let mut rows = Vec::new();
    for line_result in lines {
        let line = line_result?;
        if line.trim().is_empty() {
            continue;
        }
        let fields = line.split('\t').collect::<Vec<_>>();
        if field(&fields, &columns, "group_name", path)? != "type_2_subtype_rep_NA" {
            continue;
        }
        let reads = parse_usize_field(&fields, &columns, "subgroup_reads", path)?;
        let overlap = parse_f64_field(&fields, &columns, "mid_olps", path)?;
        if reads < REPEAT_PAIRING_MIN_PATH_READS || overlap < REPEAT_PAIRING_MIN_OVERLAP_BP {
            continue;
        }
        let boundary_key = field(&fields, &columns, "boundary_key", path)?.to_string();
        let (start, end) = parse_first_boundary_pair(&boundary_key, path)?;
        rows.push(RepeatPairingSubgroupRow {
            old_index: parse_usize_field(&fields, &columns, "old_index", path)?,
            start,
            end,
            reads,
            overlap,
        });
    }
    Ok(rows)
}

fn parse_first_boundary_pair(value: &str, path: &Path) -> Result<(usize, usize), OrgraftError> {
    let first = value.split(';').next().unwrap_or(value);
    let mut se = None;
    let mut ss = None;
    for item in first.split(',') {
        let Some((key, value)) = item.split_once('=') else {
            continue;
        };
        let parsed = value.parse::<usize>().map_err(|error| {
            OrgraftError::InvalidArgument(format!(
                "{} has invalid boundary `{item}`: {error}",
                path.display()
            ))
        })?;
        if key.starts_with("se") {
            se = Some(parsed);
        } else if key.starts_with("ss") {
            ss = Some(parsed);
        }
    }
    match (se, ss) {
        (Some(se), Some(ss)) => Ok((se, ss)),
        _ => Err(OrgraftError::InvalidArgument(format!(
            "{} has invalid repeat boundary key `{value}`",
            path.display()
        ))),
    }
}

fn collapse_reciprocal_repeat_paths(rows: &[RepeatPairingSubgroupRow]) -> Vec<PhysicalRepeatPath> {
    let mut grouped: BTreeMap<(usize, usize), Vec<RepeatPairingSubgroupRow>> = BTreeMap::new();
    for row in rows {
        grouped
            .entry(ordered_pair(row.start, row.end))
            .or_default()
            .push(row.clone());
    }
    grouped
        .into_iter()
        .map(|(endpoints, rows)| PhysicalRepeatPath {
            endpoints,
            overlap: rows.iter().map(|row| row.overlap).sum::<f64>() / rows.len() as f64,
            rows,
        })
        .collect()
}

fn path_specs(path: &PhysicalRepeatPath) -> Vec<String> {
    path.rows
        .iter()
        .map(|row| format!("type_2_subtype_rep_NA:{}", row.old_index))
        .collect()
}

fn read_index_rows(path: &Path) -> Result<Vec<RepeatReadIndexRow>, OrgraftError> {
    let file = File::open(path)?;
    let mut lines = BufReader::new(file).lines();
    let header = lines
        .next()
        .transpose()?
        .ok_or_else(|| OrgraftError::InvalidArgument(format!("{} is empty", path.display())))?;
    let columns = header_columns(&header);
    let mut rows = Vec::new();
    for line_result in lines {
        let line = line_result?;
        if line.trim().is_empty() {
            continue;
        }
        let fields = line.split('\t').collect::<Vec<_>>();
        let old_index = field(&fields, &columns, "subgroup_old_index", path)?;
        rows.push(RepeatReadIndexRow {
            read_id: field(&fields, &columns, "read_id", path)?.to_string(),
            read_class: field(&fields, &columns, "read_class", path)?.to_string(),
            group_name: field(&fields, &columns, "group_name", path)?.to_string(),
            old_index: if old_index == "." {
                None
            } else {
                Some(old_index.parse::<usize>().map_err(|error| {
                    OrgraftError::InvalidArgument(format!(
                        "{} has invalid subgroup_old_index `{old_index}`: {error}",
                        path.display()
                    ))
                })?)
            },
        });
    }
    Ok(rows)
}

fn ids_for_subgroup_specs(rows: &[RepeatReadIndexRow], specs: &[String]) -> HashSet<String> {
    let specs = specs
        .iter()
        .filter_map(|spec| spec.rsplit_once(':'))
        .filter_map(|(group, index)| Some((group, index.parse::<usize>().ok()?)))
        .collect::<HashSet<_>>();
    rows.iter()
        .filter(|row| {
            row.read_class == "FL"
                && row
                    .old_index
                    .is_some_and(|index| specs.contains(&(row.group_name.as_str(), index)))
        })
        .map(|row| row.read_id.clone())
        .collect()
}

fn read_primary_evidence(path: &Path) -> Result<Vec<PrimaryEvidenceRow>, OrgraftError> {
    let file = File::open(path)?;
    let mut lines = BufReader::new(file).lines();
    let header = lines
        .next()
        .transpose()?
        .ok_or_else(|| OrgraftError::InvalidArgument(format!("{} is empty", path.display())))?;
    let columns = header_columns(&header);
    let mut rows = Vec::new();
    for line_result in lines {
        let line = line_result?;
        if line.trim().is_empty() {
            continue;
        }
        let fields = line.split('\t').collect::<Vec<_>>();
        if field(&fields, &columns, "alignment_role", path)? == "secondary"
            || field(&fields, &columns, "target_start", path)? == "."
        {
            continue;
        }
        rows.push(PrimaryEvidenceRow {
            read_id: field(&fields, &columns, "read_id", path)?.to_string(),
            target_start: parse_usize_field(&fields, &columns, "target_start", path)?,
            target_end: parse_usize_field(&fields, &columns, "target_end", path)?,
        });
    }
    Ok(rows)
}

fn spanning_type1_reads(
    evidence: &[PrimaryEvidenceRow],
    type1_ids: &HashSet<String>,
    interval: (usize, usize),
) -> HashSet<String> {
    let required_start = interval.0.saturating_sub(REPEAT_PAIRING_FLANK_BP);
    let required_end = interval.1.saturating_add(REPEAT_PAIRING_FLANK_BP);
    evidence
        .iter()
        .filter(|row| {
            type1_ids.contains(&row.read_id)
                && row.target_start <= required_start
                && row.target_end >= required_end
        })
        .map(|row| row.read_id.clone())
        .collect()
}

fn unique_repeat_node_for_length(
    path: &Path,
    repeat_length: usize,
) -> Result<Option<String>, OrgraftError> {
    let file = File::open(path)?;
    let mut lengths = BTreeMap::new();
    let mut degree = HashMap::<String, usize>::new();
    let mut links = HashSet::new();
    for line_result in BufReader::new(file).lines() {
        let line = line_result?;
        let fields = line.split('\t').collect::<Vec<_>>();
        match fields.first().copied() {
            Some("S") if fields.len() >= 3 => {
                let length = if fields[2] != "*" {
                    fields[2].len()
                } else {
                    fields[3..]
                        .iter()
                        .find_map(|field| field.strip_prefix("LN:i:"))
                        .and_then(|value| value.parse::<usize>().ok())
                        .unwrap_or(0)
                };
                lengths.insert(fields[1].to_string(), length);
            }
            Some("L") if fields.len() >= 5 => {
                let left = format!("{}{}", fields[1], fields[2]);
                let right = format!("{}{}", fields[3], fields[4]);
                let key = if left <= right {
                    (left, right)
                } else {
                    (right, left)
                };
                if links.insert(key) {
                    *degree.entry(fields[1].to_string()).or_default() += 1;
                    *degree.entry(fields[3].to_string()).or_default() += 1;
                }
            }
            _ => {}
        }
    }
    let tolerance = ((repeat_length as f64 * 0.15).round() as usize).max(50);
    let mut candidates = lengths
        .into_iter()
        .filter(|(name, length)| {
            degree.get(name).copied().unwrap_or(0) >= 4
                && length.abs_diff(repeat_length) <= tolerance
        })
        .map(|(name, _)| name)
        .collect::<Vec<_>>();
    candidates.sort();
    Ok((candidates.len() == 1).then(|| candidates.remove(0)))
}

fn ordered_pair(left: usize, right: usize) -> (usize, usize) {
    if left <= right {
        (left, right)
    } else {
        (right, left)
    }
}

fn read_segments(path: &Path) -> Result<Vec<Segment>, OrgraftError> {
    let file = File::open(path).map_err(|error| {
        OrgraftError::InvalidArgument(format!(
            "could not read SV segment table {}: {error}",
            path.display()
        ))
    })?;
    let mut lines = BufReader::new(file).lines();
    let header = lines
        .next()
        .transpose()?
        .ok_or_else(|| OrgraftError::InvalidArgument(format!("{} is empty", path.display())))?;
    let columns = header_columns(&header);
    let mut segments = Vec::new();
    for line_result in lines {
        let line = line_result?;
        if line.trim().is_empty() {
            continue;
        }
        let fields = line.split('\t').collect::<Vec<_>>();
        let subgroup_value = field(&fields, &columns, "subgroup_old_index", path)?;
        let strand = field(&fields, &columns, "strand", path)?
            .chars()
            .next()
            .unwrap_or('+');
        if !matches!(strand, '+' | '-') {
            return Err(OrgraftError::InvalidArgument(format!(
                "{} contains invalid strand `{strand}`",
                path.display()
            )));
        }
        segments.push(Segment {
            read_id: field(&fields, &columns, "read_id", path)?.to_string(),
            read_class: field(&fields, &columns, "read_class", path)?.to_string(),
            group_name: field(&fields, &columns, "group_name", path)?.to_string(),
            subgroup_old_index: if subgroup_value == "." {
                None
            } else {
                Some(subgroup_value.parse::<usize>().map_err(|error| {
                    OrgraftError::InvalidArgument(format!(
                        "{} has invalid subgroup_old_index `{subgroup_value}`: {error}",
                        path.display()
                    ))
                })?)
            },
            segment_index: parse_usize_field(&fields, &columns, "segment_index", path)?,
            segment_count: parse_usize_field(&fields, &columns, "segment_count", path)?,
            query_start: parse_usize_field(&fields, &columns, "query_start", path)?,
            query_end: parse_usize_field(&fields, &columns, "query_end", path)?,
            subject_start: parse_usize_field(&fields, &columns, "subject_start", path)?,
            subject_end: parse_usize_field(&fields, &columns, "subject_end", path)?,
            strand,
        });
    }
    Ok(segments)
}

fn read_subgroup_ids(
    path: &Path,
    selected: &HighSubgroup,
) -> Result<HashSet<String>, OrgraftError> {
    let file = File::open(path)?;
    let mut lines = BufReader::new(file).lines();
    let header = lines
        .next()
        .transpose()?
        .ok_or_else(|| OrgraftError::InvalidArgument(format!("{} is empty", path.display())))?;
    let columns = header_columns(&header);
    let mut ids = HashSet::new();
    for line_result in lines {
        let line = line_result?;
        let fields = line.split('\t').collect::<Vec<_>>();
        let group = field(&fields, &columns, "group_name", path)?;
        let old_index = field(&fields, &columns, "subgroup_old_index", path)?;
        if group == selected.group_name && old_index == selected.old_index.to_string() {
            ids.insert(field(&fields, &columns, "read_id", path)?.to_string());
        }
    }
    Ok(ids)
}

fn generate_candidates(
    reference: &str,
    segments: &[Segment],
    selected: &HighSubgroup,
    target_ids: &HashSet<String>,
) -> Result<Vec<CandidateDescriptor>, OrgraftError> {
    let selected_by_read = group_segments(segments.iter().filter(|segment| {
        segment.group_name == selected.group_name
            && segment.subgroup_old_index == Some(selected.old_index)
    }));
    let representative = choose_representative(&selected_by_read)?;
    let normalized = normalize_representative(representative)?;
    let forced_coordinates = forced_coordinates(&normalized, reference.len())?;
    let outer_start = normalized
        .first()
        .map(|segment| segment.subject_start)
        .ok_or_else(|| {
            OrgraftError::InvalidArgument(
                "selected SV subgroup has no alignment segments".to_string(),
            )
        })?;
    let outer_end = normalized
        .last()
        .map(|segment| segment.subject_end)
        .unwrap_or(outer_start);

    let successor = |position: usize| {
        if position == reference.len() {
            1
        } else {
            position + 1
        }
    };
    let mut starts = vec![1, outer_start, successor(outer_end)];
    for &(source, target) in &forced_coordinates {
        starts.push(successor(source));
        starts.push(target);
    }
    starts.sort_unstable();
    starts.dedup();
    let blocks = starts
        .iter()
        .enumerate()
        .map(|(index, start)| Block {
            id: index,
            start: *start,
            end: if index + 1 == starts.len() {
                reference.len()
            } else {
                starts[index + 1] - 1
            },
        })
        .collect::<Vec<_>>();
    let block_by_start = blocks
        .iter()
        .map(|block| (block.start, block.id))
        .collect::<HashMap<_, _>>();
    let block_by_end = blocks
        .iter()
        .map(|block| (block.end, block.id))
        .collect::<HashMap<_, _>>();
    let mut forced = HashMap::new();
    let mut forced_in = HashMap::new();
    for (source, target) in forced_coordinates {
        let from = *block_by_end.get(&source).ok_or_else(|| {
            OrgraftError::InvalidArgument(format!(
                "SV crossover source {source} is not a block end"
            ))
        })?;
        let to = *block_by_start.get(&target).ok_or_else(|| {
            OrgraftError::InvalidArgument(format!(
                "SV crossover target {target} is not a block start"
            ))
        })?;
        if forced
            .insert(from, to)
            .is_some_and(|existing| existing != to)
            || forced_in
                .insert(to, from)
                .is_some_and(|existing| existing != from)
        {
            return Err(OrgraftError::InvalidArgument(
                "selected SV subgroup implies conflicting crossover edges".to_string(),
            ));
        }
    }
    let components = forced_components(&blocks, &forced, &forced_in)?;
    let read_paths = build_read_paths(segments, &blocks);
    let target_paths = target_ids
        .iter()
        .filter_map(|id| read_paths.get(id))
        .cloned()
        .collect::<Vec<_>>();

    let anchor_component = components
        .iter()
        .position(|component| component.contains(&0))
        .unwrap_or(0);
    let mut anchored_components = components;
    let anchor = anchored_components.remove(anchor_component);
    anchored_components.insert(0, anchor);
    let mut component_orders = Vec::new();
    let mut current = Vec::new();
    let mut used = vec![false; anchored_components.len().saturating_sub(1)];
    permute_component_indices(
        anchored_components.len().saturating_sub(1),
        &mut current,
        &mut used,
        &mut component_orders,
    );

    let current_key = canonical_token_key(
        &(0..blocks.len())
            .map(|id| BlockToken { id, orient: '+' })
            .collect::<Vec<_>>(),
    );
    let mut seen = BTreeSet::new();
    let mut candidates = Vec::new();
    for order in component_orders {
        let mut ordered = vec![anchored_components[0].clone()];
        for index in order {
            ordered.push(anchored_components[index + 1].clone());
        }
        let reversible = ordered
            .iter()
            .enumerate()
            .filter_map(|(index, component)| (component.len() == 1).then_some(index))
            .collect::<Vec<_>>();
        if reversible.len() >= usize::BITS as usize {
            return Err(OrgraftError::InvalidArgument(
                "SV candidate orientation search is too large".to_string(),
            ));
        }
        for mask in 0..(1usize << reversible.len()) {
            let mut tokens = Vec::with_capacity(blocks.len());
            for (component_index, component) in ordered.iter().enumerate() {
                let reverse = reversible
                    .iter()
                    .position(|index| *index == component_index)
                    .is_some_and(|bit| mask & (1usize << bit) != 0);
                if reverse {
                    for id in component.iter().rev() {
                        tokens.push(BlockToken {
                            id: *id,
                            orient: '-',
                        });
                    }
                } else {
                    for id in component {
                        tokens.push(BlockToken {
                            id: *id,
                            orient: '+',
                        });
                    }
                }
            }
            let key = canonical_token_key(&tokens);
            if key == current_key || !seen.insert(key) {
                continue;
            }
            let predicted_fl_reads = read_paths
                .values()
                .filter(|path| path_is_contiguous(path, &tokens))
                .count();
            let predicted_target_reads = target_paths
                .iter()
                .filter(|path| path_is_contiguous(path, &tokens))
                .count();
            let original_adjacencies = count_original_adjacencies(&tokens, blocks.len());
            let reversed_bases = tokens
                .iter()
                .filter(|token| token.orient == '-')
                .map(|token| blocks[token.id].len())
                .sum();
            candidates.push(CandidateDescriptor {
                tokens,
                predicted_fl_reads,
                predicted_target_reads,
                original_adjacencies,
                reversed_bases,
            });
            if candidates.len() > MAX_GENERATED_CANDIDATES {
                return Err(OrgraftError::InvalidArgument(format!(
                    "SV subgroup `{}` generated more than {MAX_GENERATED_CANDIDATES} candidates; specify a simpler subgroup or edit manually",
                    selected.spec()
                )));
            }
        }
    }
    candidates.sort_by(|left, right| {
        right
            .predicted_target_reads
            .cmp(&left.predicted_target_reads)
            .then_with(|| right.predicted_fl_reads.cmp(&left.predicted_fl_reads))
            .then_with(|| right.original_adjacencies.cmp(&left.original_adjacencies))
            .then_with(|| left.reversed_bases.cmp(&right.reversed_bases))
            .then_with(|| left.tokens.cmp(&right.tokens))
    });
    diversify_candidates(candidates)
}

fn diversify_candidates(
    candidates: Vec<CandidateDescriptor>,
) -> Result<Vec<CandidateDescriptor>, OrgraftError> {
    if candidates.len() <= DEFAULT_CANDIDATE_EVALUATIONS {
        return Ok(candidates);
    }
    let mut buckets: BTreeMap<(usize, usize), Vec<CandidateDescriptor>> = BTreeMap::new();
    for candidate in candidates {
        buckets
            .entry((
                candidate.predicted_target_reads,
                candidate.original_adjacencies,
            ))
            .or_default()
            .push(candidate);
    }
    let mut keys = buckets.keys().copied().collect::<Vec<_>>();
    keys.sort_by(|left, right| right.cmp(left));
    let mut output = Vec::new();
    loop {
        let mut added = false;
        for key in &keys {
            if let Some(bucket) = buckets.get_mut(key) {
                if !bucket.is_empty() {
                    output.push(bucket.remove(0));
                    added = true;
                }
            }
        }
        if !added {
            break;
        }
    }
    Ok(output)
}

fn group_segments<'a>(
    segments: impl Iterator<Item = &'a Segment>,
) -> BTreeMap<String, Vec<Segment>> {
    let mut grouped: BTreeMap<String, Vec<Segment>> = BTreeMap::new();
    for segment in segments {
        grouped
            .entry(segment.read_id.clone())
            .or_default()
            .push(segment.clone());
    }
    for rows in grouped.values_mut() {
        rows.sort_by_key(|segment| segment.segment_index);
    }
    grouped
}

fn choose_representative(
    grouped: &BTreeMap<String, Vec<Segment>>,
) -> Result<&[Segment], OrgraftError> {
    let mut complete = grouped
        .values()
        .filter(|rows| {
            rows.first()
                .is_some_and(|first| rows.len() == first.segment_count)
        })
        .collect::<Vec<_>>();
    if complete.is_empty() {
        return Err(OrgraftError::InvalidArgument(
            "selected SV subgroup has no complete read alignment path".to_string(),
        ));
    }
    complete.sort_by_key(|rows| {
        rows.iter()
            .map(|segment| segment.query_end)
            .max()
            .unwrap_or(0)
    });
    Ok(complete[complete.len() / 2])
}

fn normalize_representative(rows: &[Segment]) -> Result<Vec<Segment>, OrgraftError> {
    let strands = rows
        .iter()
        .map(|segment| segment.strand)
        .collect::<BTreeSet<_>>();
    if strands.len() != 1 {
        return Err(OrgraftError::InvalidArgument(
            "automatic SV correction currently requires a subgroup with one common strand; select a simpler subgroup or edit manually"
                .to_string(),
        ));
    }
    if strands.contains(&'+') {
        return Ok(rows.to_vec());
    }
    let read_length = rows
        .iter()
        .map(|segment| segment.query_end)
        .max()
        .unwrap_or(0);
    Ok(rows
        .iter()
        .rev()
        .map(|segment| Segment {
            read_id: segment.read_id.clone(),
            read_class: segment.read_class.clone(),
            group_name: segment.group_name.clone(),
            subgroup_old_index: segment.subgroup_old_index,
            segment_index: segment.segment_count - segment.segment_index + 1,
            segment_count: segment.segment_count,
            query_start: read_length - segment.query_end + 1,
            query_end: read_length - segment.query_start + 1,
            subject_start: segment.subject_end,
            subject_end: segment.subject_start,
            strand: '+',
        })
        .collect())
}

fn forced_coordinates(
    rows: &[Segment],
    reference_len: usize,
) -> Result<Vec<(usize, usize)>, OrgraftError> {
    let mut output = Vec::new();
    for pair in rows.windows(2) {
        let left = &pair[0];
        let right = &pair[1];
        if right.query_start > left.query_end + 1 {
            return Err(OrgraftError::InvalidArgument(format!(
                "selected SV subgroup has an unaligned query gap between segments {} and {}",
                left.segment_index, right.segment_index
            )));
        }
        let (left_query, right_query) = if right.query_start <= left.query_end {
            let midpoint = (right.query_start + left.query_end) / 2;
            (midpoint, midpoint + 1)
        } else {
            (left.query_end, right.query_start)
        };
        let source = map_query_coordinate(left, left_query).clamp(1, reference_len);
        let target = map_query_coordinate(right, right_query).clamp(1, reference_len);
        output.push((source, target));
    }
    Ok(output)
}

fn map_query_coordinate(segment: &Segment, query: usize) -> usize {
    let query_span = segment.query_end.saturating_sub(segment.query_start);
    if query_span == 0 {
        return segment.subject_start;
    }
    let subject_span = segment.subject_end as isize - segment.subject_start as isize;
    let offset = query.saturating_sub(segment.query_start) as f64;
    let mapped = segment.subject_start as f64 + offset * subject_span as f64 / query_span as f64;
    mapped.round().max(1.0) as usize
}

fn forced_components(
    blocks: &[Block],
    forced: &HashMap<usize, usize>,
    forced_in: &HashMap<usize, usize>,
) -> Result<Vec<Vec<usize>>, OrgraftError> {
    let starts = blocks
        .iter()
        .map(|block| block.id)
        .filter(|id| !forced_in.contains_key(id))
        .collect::<Vec<_>>();
    let mut covered = HashSet::new();
    let mut components = Vec::new();
    for start in starts {
        let mut component = Vec::new();
        let mut current = start;
        loop {
            if !covered.insert(current) {
                return Err(OrgraftError::InvalidArgument(
                    "selected SV subgroup implies a forced edge cycle".to_string(),
                ));
            }
            component.push(current);
            let Some(next) = forced.get(&current) else {
                break;
            };
            current = *next;
        }
        components.push(component);
    }
    if covered.len() != blocks.len() {
        return Err(OrgraftError::InvalidArgument(
            "selected SV subgroup does not define valid reference blocks".to_string(),
        ));
    }
    Ok(components)
}

fn permute_component_indices(
    count: usize,
    current: &mut Vec<usize>,
    used: &mut [bool],
    output: &mut Vec<Vec<usize>>,
) {
    if current.len() == count {
        output.push(current.clone());
        return;
    }
    for index in 0..count {
        if used[index] {
            continue;
        }
        used[index] = true;
        current.push(index);
        permute_component_indices(count, current, used, output);
        current.pop();
        used[index] = false;
    }
}

fn build_read_paths(segments: &[Segment], blocks: &[Block]) -> HashMap<String, Vec<BlockToken>> {
    let grouped = group_segments(segments.iter().filter(|segment| segment.read_class == "FL"));
    grouped
        .into_iter()
        .map(|(read_id, rows)| {
            let cuts = rows
                .windows(2)
                .map(|pair| {
                    if pair[1].query_start <= pair[0].query_end {
                        (pair[1].query_start + pair[0].query_end) / 2
                    } else {
                        pair[0].query_end
                    }
                })
                .collect::<Vec<_>>();
            let mut path = Vec::new();
            for (index, segment) in rows.iter().enumerate() {
                let effective_start =
                    if index > 0 && segment.query_start <= rows[index - 1].query_end {
                        cuts[index - 1] + 1
                    } else {
                        segment.query_start
                    };
                let effective_end =
                    if index < cuts.len() && rows[index + 1].query_start <= segment.query_end {
                        cuts[index]
                    } else {
                        segment.query_end
                    };
                let start = map_query_coordinate(segment, effective_start);
                let end = map_query_coordinate(segment, effective_end);
                let low = start.min(end);
                let high = start.max(end);
                let mut ids = blocks
                    .iter()
                    .filter(|block| block.end >= low && block.start <= high)
                    .map(|block| block.id)
                    .collect::<Vec<_>>();
                if segment.strand == '-' {
                    ids.reverse();
                }
                for id in ids {
                    let token = BlockToken {
                        id,
                        orient: segment.strand,
                    };
                    if path.last() != Some(&token) {
                        path.push(token);
                    }
                }
            }
            (read_id, path)
        })
        .collect()
}

fn path_is_contiguous(path: &[BlockToken], candidate: &[BlockToken]) -> bool {
    if path.is_empty() || candidate.is_empty() {
        return false;
    }
    let reverse = path
        .iter()
        .rev()
        .map(|token| BlockToken {
            id: token.id,
            orient: flip_orient(token.orient),
        })
        .collect::<Vec<_>>();
    [path, reverse.as_slice()].iter().any(|oriented| {
        candidate
            .iter()
            .enumerate()
            .filter(|(_, token)| **token == oriented[0])
            .any(|(start, _)| {
                oriented
                    .iter()
                    .enumerate()
                    .all(|(offset, token)| candidate[(start + offset) % candidate.len()] == *token)
            })
    })
}

fn canonical_token_key(tokens: &[BlockToken]) -> Vec<i64> {
    let encode = |token: BlockToken| {
        if token.orient == '+' {
            token.id as i64 + 1
        } else {
            -(token.id as i64 + 1)
        }
    };
    let forward = tokens.iter().copied().map(encode).collect::<Vec<_>>();
    let reverse = tokens
        .iter()
        .rev()
        .map(|token| {
            encode(BlockToken {
                id: token.id,
                orient: flip_orient(token.orient),
            })
        })
        .collect::<Vec<_>>();
    let mut best: Option<Vec<i64>> = None;
    for values in [&forward, &reverse] {
        for offset in 0..values.len() {
            let rotated = values[offset..]
                .iter()
                .chain(values[..offset].iter())
                .copied()
                .collect::<Vec<_>>();
            if best.as_ref().is_none_or(|current| rotated < *current) {
                best = Some(rotated);
            }
        }
    }
    best.unwrap_or_default()
}

fn count_original_adjacencies(tokens: &[BlockToken], block_count: usize) -> usize {
    tokens
        .iter()
        .enumerate()
        .filter(|(index, left)| {
            let right = tokens[(*index + 1) % tokens.len()];
            (left.orient == '+' && right.orient == '+' && right.id == (left.id + 1) % block_count)
                || (left.orient == '-'
                    && right.orient == '-'
                    && right.id == (left.id + block_count - 1) % block_count)
        })
        .count()
}

fn materialize_candidate(
    reference: &str,
    descriptor: &CandidateDescriptor,
    selected: &HighSubgroup,
    segments: &[Segment],
) -> Result<String, OrgraftError> {
    let selected_by_read = group_segments(segments.iter().filter(|segment| {
        segment.group_name == selected.group_name
            && segment.subgroup_old_index == Some(selected.old_index)
    }));
    let representative = choose_representative(&selected_by_read)?;
    let normalized = normalize_representative(representative)?;
    let forced_coordinates = forced_coordinates(&normalized, reference.len())?;
    let outer_start = normalized[0].subject_start;
    let outer_end = normalized.last().unwrap().subject_end;
    let successor = |position: usize| {
        if position == reference.len() {
            1
        } else {
            position + 1
        }
    };
    let mut starts = vec![1, outer_start, successor(outer_end)];
    for &(source, target) in &forced_coordinates {
        starts.push(successor(source));
        starts.push(target);
    }
    starts.sort_unstable();
    starts.dedup();
    let blocks = starts
        .iter()
        .enumerate()
        .map(|(index, start)| Block {
            id: index,
            start: *start,
            end: if index + 1 == starts.len() {
                reference.len()
            } else {
                starts[index + 1] - 1
            },
        })
        .collect::<Vec<_>>();
    let block_by_start = blocks
        .iter()
        .map(|block| (block.start, block.id))
        .collect::<HashMap<_, _>>();
    let block_by_end = blocks
        .iter()
        .map(|block| (block.end, block.id))
        .collect::<HashMap<_, _>>();
    let mut forced_blocks = HashSet::new();
    for (source, target) in forced_coordinates {
        forced_blocks.insert(*block_by_end.get(&source).unwrap());
        forced_blocks.insert(*block_by_start.get(&target).unwrap());
    }

    let token_sequences = descriptor
        .tokens
        .iter()
        .map(|token| {
            let block = &blocks[token.id];
            let sequence = &reference[block.start - 1..block.end];
            if token.orient == '+' {
                sequence.to_string()
            } else {
                reverse_complement(sequence)
            }
        })
        .collect::<Vec<_>>();
    let sequence = token_sequences.join("");
    let origin_token = descriptor
        .tokens
        .iter()
        .enumerate()
        .filter(|(_, token)| !forced_blocks.contains(&token.id))
        .max_by_key(|(index, _)| token_sequences[*index].len())
        .map(|(index, _)| index)
        .or_else(|| {
            token_sequences
                .iter()
                .enumerate()
                .max_by_key(|(_, sequence)| sequence.len())
                .map(|(index, _)| index)
        })
        .unwrap_or(0);
    let origin_offset = token_sequences[..origin_token]
        .iter()
        .map(String::len)
        .sum::<usize>()
        + token_sequences[origin_token].len() / 2;
    Ok(format!(
        "{}{}",
        &sequence[origin_offset..],
        &sequence[..origin_offset]
    ))
}

fn reverse_complement(sequence: &str) -> String {
    sequence
        .bytes()
        .rev()
        .map(|base| match base.to_ascii_uppercase() {
            b'A' => 'T',
            b'C' => 'G',
            b'G' => 'C',
            b'T' | b'U' => 'A',
            _ => 'N',
        })
        .collect()
}

fn read_prior_sequences(paths: &[PathBuf]) -> Result<Vec<String>, OrgraftError> {
    paths
        .iter()
        .filter(|path| path.is_file())
        .map(|path| read_single_fasta(path).map(|(_, sequence)| sequence))
        .collect()
}

fn is_historical_sequence(sequence: &str, history: &[String]) -> bool {
    history
        .iter()
        .any(|previous| circular_equivalent(sequence, previous))
}

fn circular_equivalent(left: &str, right: &str) -> bool {
    if left.len() != right.len() {
        return false;
    }
    let doubled = format!("{right}{right}");
    if doubled.contains(left) {
        return true;
    }
    let reverse = reverse_complement(right);
    format!("{reverse}{reverse}").contains(left)
}

fn count_type1_reads(path: &Path, target_ids: &HashSet<String>) -> Result<usize, OrgraftError> {
    Ok(type1_read_ids(path)?.intersection(target_ids).count())
}

fn type1_read_ids(path: &Path) -> Result<HashSet<String>, OrgraftError> {
    let file = File::open(path)?;
    let mut lines = BufReader::new(file).lines();
    let header = lines
        .next()
        .transpose()?
        .ok_or_else(|| OrgraftError::InvalidArgument(format!("{} is empty", path.display())))?;
    let columns = header_columns(&header);
    let mut ids = HashSet::new();
    for line_result in lines {
        let line = line_result?;
        let fields = line.split('\t').collect::<Vec<_>>();
        let read_id = field(&fields, &columns, "read_id", path)?;
        let group = field(&fields, &columns, "group_name", path)?;
        if group == "type_1_subtype_NA" {
            ids.insert(read_id.to_string());
        }
    }
    Ok(ids)
}

fn candidate_improves(score: &CandidateScore, baseline: SvMetrics) -> bool {
    score.target_type1_fraction() >= MIN_TARGET_TYPE1_FRACTION
        && score.evaluation.low_green_window_fraction + IMPROVEMENT_EPSILON
            < baseline.low_green_window_fraction
        && score.evaluation.reference_support_read_fraction + MAX_GLOBAL_SUPPORT_DROP
            >= baseline.reference_support_read_fraction
        && score.evaluation.reference_support_depth_area_fraction + MAX_GLOBAL_SUPPORT_DROP
            >= baseline.reference_support_depth_area_fraction
}

fn candidate_better(left: &CandidateScore, right: &CandidateScore) -> bool {
    left.evaluation
        .low_green_window_fraction
        .partial_cmp(&right.evaluation.low_green_window_fraction)
        .unwrap_or(Ordering::Equal)
        .then_with(|| {
            right
                .evaluation
                .reference_support_depth_area_fraction
                .partial_cmp(&left.evaluation.reference_support_depth_area_fraction)
                .unwrap_or(Ordering::Equal)
        })
        .then_with(|| {
            right
                .target_type1_fraction()
                .partial_cmp(&left.target_type1_fraction())
                .unwrap_or(Ordering::Equal)
        })
        .then_with(|| {
            right
                .descriptor
                .original_adjacencies
                .cmp(&left.descriptor.original_adjacencies)
        })
        == Ordering::Less
}

fn remove_candidate_artifacts(score: &CandidateScore) -> Result<(), OrgraftError> {
    match fs::remove_file(&score.candidate_path) {
        Ok(()) => {}
        Err(error) if error.kind() == std::io::ErrorKind::NotFound => {}
        Err(error) => return Err(error.into()),
    }
    match fs::remove_dir_all(&score.evaluation_dir) {
        Ok(()) => {}
        Err(error) if error.kind() == std::io::ErrorKind::NotFound => {}
        Err(error) => return Err(error.into()),
    }
    Ok(())
}

fn append_candidate_score(output: &mut String, score: &CandidateScore, accepted: bool) {
    writeln!(
        output,
        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.6}\t{}\t{}\t{}\t{:.6}\t{:.6}\t{:.6}\t{}",
        score.rank,
        score.descriptor.predicted_fl_reads,
        score.descriptor.predicted_target_reads,
        score.descriptor.original_adjacencies,
        score.descriptor.reversed_bases,
        score.target_type1_reads,
        score.target_reads,
        score.target_type1_fraction(),
        score.evaluation.status,
        score.evaluation.fl_reads,
        score.evaluation.reference_support_reads,
        score.evaluation.reference_support_read_fraction,
        score.evaluation.reference_support_depth_area_fraction,
        score.evaluation.low_green_window_fraction,
        accepted,
    )
    .unwrap();
}

fn append_skipped_score(
    output: &mut String,
    rank: usize,
    descriptor: &CandidateDescriptor,
    reason: &str,
) {
    writeln!(
        output,
        "{rank}\t{}\t{}\t{}\t{}\t.\t.\t.\t.\t.\t.\t.\t.\t.\t{reason}",
        descriptor.predicted_fl_reads,
        descriptor.predicted_target_reads,
        descriptor.original_adjacencies,
        descriptor.reversed_bases,
    )
    .unwrap();
}

fn write_repair_report(
    path: &Path,
    selected: &HighSubgroup,
    baseline: SvMetrics,
    score: &CandidateScore,
    candidate_count: usize,
    evaluated_candidates: usize,
    corrected_fasta: &Path,
    score_table: &Path,
    graph_localization: &SvGraphLocalization,
) -> Result<(), OrgraftError> {
    let mut output = String::from("metric\tvalue\n");
    metric(&mut output, "status", "corrected");
    metric(&mut output, "subgroup", &selected.spec());
    metric(&mut output, "candidate_count", &candidate_count.to_string());
    metric(
        &mut output,
        "evaluated_candidates",
        &evaluated_candidates.to_string(),
    );
    metric(&mut output, "selected_candidate", &score.rank.to_string());
    metric(&mut output, "target_reads", &score.target_reads.to_string());
    metric(
        &mut output,
        "target_type1_reads",
        &score.target_type1_reads.to_string(),
    );
    metric(
        &mut output,
        "target_type1_fraction",
        &format!("{:.6}", score.target_type1_fraction()),
    );
    metric(
        &mut output,
        "before_low_green_window_fraction",
        &format!("{:.6}", baseline.low_green_window_fraction),
    );
    metric(
        &mut output,
        "after_low_green_window_fraction",
        &format!("{:.6}", score.evaluation.low_green_window_fraction),
    );
    metric(
        &mut output,
        "before_reference_support_read_fraction",
        &format!("{:.6}", baseline.reference_support_read_fraction),
    );
    metric(
        &mut output,
        "after_reference_support_read_fraction",
        &format!("{:.6}", score.evaluation.reference_support_read_fraction),
    );
    metric(
        &mut output,
        "before_reference_support_depth_area_fraction",
        &format!("{:.6}", baseline.reference_support_depth_area_fraction),
    );
    metric(
        &mut output,
        "after_reference_support_depth_area_fraction",
        &format!(
            "{:.6}",
            score.evaluation.reference_support_depth_area_fraction
        ),
    );
    metric(
        &mut output,
        "corrected_fasta",
        &corrected_fasta.display().to_string(),
    );
    metric(
        &mut output,
        "candidate_scores",
        &score_table.display().to_string(),
    );
    append_graph_localization_metrics(&mut output, graph_localization);
    fs::write(path, output)?;
    Ok(())
}

fn write_no_candidate_report(
    path: &Path,
    selected: &HighSubgroup,
    baseline: SvMetrics,
    candidate_count: usize,
    evaluated_candidates: usize,
    score_table: &Path,
    graph_localization: &SvGraphLocalization,
) -> Result<(), OrgraftError> {
    let mut output = String::from("metric\tvalue\n");
    metric(&mut output, "status", "manual_required");
    metric(&mut output, "subgroup", &selected.spec());
    metric(&mut output, "candidate_count", &candidate_count.to_string());
    metric(
        &mut output,
        "evaluated_candidates",
        &evaluated_candidates.to_string(),
    );
    metric(
        &mut output,
        "before_low_green_window_fraction",
        &format!("{:.6}", baseline.low_green_window_fraction),
    );
    metric(
        &mut output,
        "candidate_scores",
        &score_table.display().to_string(),
    );
    append_graph_localization_metrics(&mut output, graph_localization);
    fs::write(path, output)?;
    Ok(())
}

fn append_graph_localization_metrics(
    output: &mut String,
    graph_localization: &SvGraphLocalization,
) {
    metric(
        output,
        "graph_localization_report",
        &graph_localization.report.display().to_string(),
    );
    metric(
        output,
        "graph_problem_scope",
        &graph_localization.problem_scope,
    );
    metric(
        output,
        "graph_suspect_segments",
        &if graph_localization.suspect_segments.is_empty() {
            ".".to_string()
        } else {
            graph_localization.suspect_segments.join(",")
        },
    );
    metric(output, "graph_guidance", &graph_localization.guidance);
}

fn metric(output: &mut String, key: &str, value: &str) {
    writeln!(output, "{key}\t{}", value.replace(['\t', '\n', '\r'], " ")).unwrap();
}

fn read_single_fasta(path: &Path) -> Result<(String, String), OrgraftError> {
    let mut header = None;
    let mut sequence = String::new();
    for line_result in BufReader::new(File::open(path)?).lines() {
        let line = line_result?;
        if let Some(value) = line.strip_prefix('>') {
            if header.is_some() {
                return Err(OrgraftError::InvalidArgument(format!(
                    "{} must contain exactly one FASTA record for SV correction",
                    path.display()
                )));
            }
            header = Some(value.to_string());
        } else {
            sequence.push_str(line.trim());
        }
    }
    let header = header.ok_or_else(|| {
        OrgraftError::InvalidArgument(format!("{} contains no FASTA record", path.display()))
    })?;
    if sequence.is_empty() {
        return Err(OrgraftError::InvalidArgument(format!(
            "{} contains an empty FASTA sequence",
            path.display()
        )));
    }
    Ok((header, sequence.to_ascii_uppercase()))
}

fn write_fasta(path: &Path, header: &str, sequence: &str) -> Result<(), OrgraftError> {
    let mut output = String::new();
    writeln!(output, ">{header}").unwrap();
    for chunk in sequence.as_bytes().chunks(80) {
        output.push_str(std::str::from_utf8(chunk).map_err(|error| {
            OrgraftError::InvalidArgument(format!("candidate FASTA is not UTF-8: {error}"))
        })?);
        output.push('\n');
    }
    fs::write(path, output)?;
    Ok(())
}

fn read_metric_f64(path: &Path, metric_name: &str) -> Result<f64, OrgraftError> {
    for line_result in BufReader::new(File::open(path)?).lines() {
        let line = line_result?;
        let mut fields = line.split('\t');
        if fields.next() != Some(metric_name) {
            continue;
        }
        return fields.next().unwrap_or("").parse::<f64>().map_err(|error| {
            OrgraftError::InvalidArgument(format!(
                "{} metric `{metric_name}` is invalid: {error}",
                path.display()
            ))
        });
    }
    Err(OrgraftError::InvalidArgument(format!(
        "{} does not contain metric `{metric_name}`",
        path.display()
    )))
}

fn header_columns(header: &str) -> HashMap<&str, usize> {
    header
        .split('\t')
        .enumerate()
        .map(|(index, value)| (value, index))
        .collect()
}

fn field<'a>(
    fields: &'a [&str],
    columns: &HashMap<&str, usize>,
    name: &str,
    path: &Path,
) -> Result<&'a str, OrgraftError> {
    let index = columns.get(name).copied().ok_or_else(|| {
        OrgraftError::InvalidArgument(format!("{} is missing column `{name}`", path.display()))
    })?;
    fields.get(index).copied().ok_or_else(|| {
        OrgraftError::InvalidArgument(format!(
            "{} has a row missing column `{name}`",
            path.display()
        ))
    })
}

fn parse_usize_field(
    fields: &[&str],
    columns: &HashMap<&str, usize>,
    name: &str,
    path: &Path,
) -> Result<usize, OrgraftError> {
    let value = field(fields, columns, name, path)?;
    value.parse::<usize>().map_err(|error| {
        OrgraftError::InvalidArgument(format!(
            "{} column `{name}` has invalid value `{value}`: {error}",
            path.display()
        ))
    })
}

fn parse_f64_field(
    fields: &[&str],
    columns: &HashMap<&str, usize>,
    name: &str,
    path: &Path,
) -> Result<f64, OrgraftError> {
    let value = field(fields, columns, name, path)?;
    value.parse::<f64>().map_err(|error| {
        OrgraftError::InvalidArgument(format!(
            "{} column `{name}` has invalid value `{value}`: {error}",
            path.display()
        ))
    })
}

fn parse_bool_field(
    fields: &[&str],
    columns: &HashMap<&str, usize>,
    name: &str,
    path: &Path,
) -> Result<bool, OrgraftError> {
    match field(fields, columns, name, path)? {
        "true" => Ok(true),
        "false" => Ok(false),
        value => Err(OrgraftError::InvalidArgument(format!(
            "{} column `{name}` has invalid boolean `{value}`",
            path.display()
        ))),
    }
}

fn flip_orient(orient: char) -> char {
    if orient == '+' {
        '-'
    } else {
        '+'
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::time::{SystemTime, UNIX_EPOCH};

    #[test]
    fn shortened_group_name_matches_full_subtype_name() {
        assert_eq!(
            normalize_group_name("type_3_subtype_rep_rep_NA"),
            "type_3_rep_rep"
        );
    }

    #[test]
    fn circular_equivalence_accepts_rotation_and_reverse_complement() {
        assert!(circular_equivalent("CCCGGGAAA", "AAACCCGGG"));
        assert!(circular_equivalent("TTTCCCGGG", "AAACCCGGG"));
        assert!(!circular_equivalent("AAACCCGGA", "AAACCCGGG"));
    }

    #[test]
    fn path_contiguity_accepts_both_read_strands() {
        let candidate = vec![
            BlockToken { id: 0, orient: '+' },
            BlockToken { id: 2, orient: '-' },
            BlockToken { id: 1, orient: '+' },
        ];
        assert!(path_is_contiguous(&candidate[1..], &candidate));
        let reverse = vec![
            BlockToken { id: 1, orient: '-' },
            BlockToken { id: 2, orient: '+' },
        ];
        assert!(path_is_contiguous(&reverse, &candidate));
    }

    #[test]
    fn automatic_selection_requires_low_local_reference_support() {
        let rows = vec![
            HighSubgroup {
                group_name: "type_3_subtype_rep_rep_NA".to_string(),
                old_index: 4,
                boundary_key: "se1=1,ss2=299704;se2=299355,ss3=340100".to_string(),
                subgroup_reads: 69,
                is_reference_support: false,
                auto_highlight: true,
                min_reference_fraction: 0.019,
                judgement: "possible_reference_sv_error".to_string(),
            },
            HighSubgroup {
                group_name: "type_3_subtype_rep_rep_NA".to_string(),
                old_index: 9,
                boundary_key: "se1=10,ss2=20".to_string(),
                subgroup_reads: 90,
                is_reference_support: false,
                auto_highlight: true,
                min_reference_fraction: 0.70,
                judgement: "minor_recombination_or_alternative_configuration".to_string(),
            },
        ];
        let selected = select_subgroup(&rows, None).unwrap().unwrap();
        assert_eq!(selected.old_index, 4);
    }

    #[test]
    fn automatic_selection_skips_high_subgroup_with_green_support() {
        let rows = vec![HighSubgroup {
            group_name: "type_3_subtype_ref_rep_NA".to_string(),
            old_index: 9,
            boundary_key: "se1=10,ss2=20".to_string(),
            subgroup_reads: 90,
            is_reference_support: false,
            auto_highlight: true,
            min_reference_fraction: 0.70,
            judgement: "minor_recombination_or_alternative_configuration".to_string(),
        }];
        assert!(select_subgroup(&rows, None).unwrap().is_none());
    }

    #[test]
    fn manual_selection_accepts_shortened_group_name_and_disappears_after_repair() {
        let rows = vec![HighSubgroup {
            group_name: "type_3_subtype_rep_rep_NA".to_string(),
            old_index: 4,
            boundary_key: "se1=1,ss2=299704;se2=299355,ss3=340100".to_string(),
            subgroup_reads: 69,
            is_reference_support: false,
            auto_highlight: true,
            min_reference_fraction: 0.019,
            judgement: "possible_reference_sv_error".to_string(),
        }];
        assert_eq!(
            select_subgroup(&rows, Some("type_3_rep_rep:4"))
                .unwrap()
                .unwrap()
                .old_index,
            4
        );
        assert!(select_subgroup(&rows, Some("type_3_rep_rep:5"))
            .unwrap()
            .is_none());
    }

    #[test]
    fn repeat_pairing_detector_selects_the_dominant_minor_pairing_repeat() {
        let dir = test_dir("repeat_pairing_detector");
        fs::create_dir_all(&dir).unwrap();
        let subgroup_stats = dir.join("sv_subgroup_stats.tsv");
        let read_index = dir.join("sv_read_index.tsv");
        let read_evidence = dir.join("sv_read_evidence.tsv");
        let graph = dir.join("graph.gfa");

        fs::write(
            &subgroup_stats,
            concat!(
                "group_name\told_index\tboundary_key\tsubgroup_reads\tmid_olps\n",
                "type_2_subtype_rep_NA\t1\tse1=10000,ss2=50000\t12\t1000\n",
                "type_2_subtype_rep_NA\t2\tse1=50000,ss2=10000\t12\t1000\n",
                "type_2_subtype_rep_NA\t3\tse1=11000,ss2=51000\t12\t1000\n",
                "type_2_subtype_rep_NA\t4\tse1=51000,ss2=11000\t12\t1000\n",
            ),
        )
        .unwrap();
        let mut index = String::from("read_id\tread_class\tgroup_name\tsubgroup_old_index\n");
        for old_index in 1..=4 {
            for read in 1..=12 {
                writeln!(
                    index,
                    "alt_{old_index}_{read}\tFL\ttype_2_subtype_rep_NA\t{old_index}"
                )
                .unwrap();
            }
        }
        for read in 1..=4 {
            writeln!(index, "ref_{read}\tFL\ttype_1_subtype_NA\t.").unwrap();
        }
        fs::write(&read_index, index).unwrap();
        fs::write(
            &read_evidence,
            concat!(
                "read_id\talignment_role\ttarget_start\ttarget_end\n",
                "ref_1\tprimary\t8000\t13000\n",
                "ref_2\tprimary\t8000\t13000\n",
                "ref_3\tprimary\t48000\t53000\n",
                "ref_4\tprimary\t48000\t53000\n",
            ),
        )
        .unwrap();
        let repeat = "A".repeat(1_000);
        let flank = "C".repeat(2_000);
        fs::write(
            &graph,
            format!(
                "S\tutg5\t{repeat}\nS\tleft_a\t{flank}\nS\tleft_b\t{flank}\nS\tright_a\t{flank}\nS\tright_b\t{flank}\nL\tleft_a\t+\tutg5\t+\t0M\nL\tleft_b\t+\tutg5\t+\t0M\nL\tutg5\t+\tright_a\t+\t0M\nL\tutg5\t+\tright_b\t+\t0M\n"
            ),
        )
        .unwrap();

        let mismatch = detect_repeat_pairing_mismatch(
            &subgroup_stats,
            &read_index,
            &read_evidence,
            Some(&graph),
        )
        .unwrap()
        .unwrap();
        assert_eq!(mismatch.repeat_node, "utg5");
        assert_eq!(mismatch.repeat_length, 1_000);
        assert_eq!(mismatch.alternate_path_reads, [24, 24]);
        assert_eq!(mismatch.reference_spanning_reads, [2, 2]);
        assert!(mismatch.dominance_ratio >= 3.0);
        let _ = fs::remove_dir_all(dir);
    }

    fn test_dir(name: &str) -> PathBuf {
        let nanos = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .map(|value| value.as_nanos())
            .unwrap_or(0);
        std::env::temp_dir().join(format!("orgraft_{name}_{nanos}"))
    }
}
