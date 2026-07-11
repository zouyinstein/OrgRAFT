use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet};
use std::fs::{self, File};
use std::io::{BufRead, BufReader, Write};
use std::path::{Path, PathBuf};
use std::process::Command;

use crate::error::OrgraftError;

const ENDPOINT_SLOP: usize = 1_000;
const MAX_PROJECTION_DISTANCE: usize = 1_000;

type PhysicalEnd = (String, char);
type PhysicalLink = (PhysicalEnd, PhysicalEnd);

#[derive(Debug, Clone)]
pub(crate) struct SvGraphLocalization {
    pub report: PathBuf,
    pub problem_scope: String,
    pub suspect_segments: Vec<String>,
    pub guidance: String,
}

pub(crate) struct SvGraphLocalizationRequest<'a> {
    pub subgroup: &'a str,
    pub boundary_key: &'a str,
    pub reference: &'a Path,
    pub unitig_gfa: &'a Path,
    pub soft_paths: &'a Path,
    pub output_dir: &'a Path,
    pub threads: usize,
}

#[derive(Debug, Clone)]
struct GfaSegment {
    name: String,
    sequence: String,
    coverage: Option<f64>,
}

#[derive(Debug)]
struct UnitigGraph {
    segments: Vec<GfaSegment>,
    segment_by_name: HashMap<String, GfaSegment>,
    links: HashSet<PhysicalLink>,
    incident: HashMap<PhysicalEnd, BTreeSet<String>>,
}

#[derive(Debug, Clone)]
struct BoundaryPosition {
    pair_index: usize,
    role: String,
    position: usize,
}

#[derive(Debug, Clone)]
struct PafAlignment {
    query: String,
    query_len: usize,
    query_start: usize,
    query_end: usize,
    strand: char,
    target_len: usize,
    target_start: usize,
    target_end: usize,
    matches: usize,
    mapq: u8,
}

#[derive(Debug, Clone)]
struct Projection {
    boundary: BoundaryPosition,
    segment: String,
    segment_offset: usize,
    segment_len: usize,
    segment_side: char,
    strand: char,
    target_distance: usize,
    coverage: Option<f64>,
    incident_links: String,
    edge_status: String,
}

pub(crate) fn localize_sv_to_unitig_graph(
    request: &SvGraphLocalizationRequest<'_>,
) -> Result<SvGraphLocalization, OrgraftError> {
    fs::create_dir_all(request.output_dir)?;
    let graph = read_unitig_graph(request.unitig_gfa)?;
    let boundaries = parse_boundary_key(request.boundary_key)?;
    let segments_fasta = request.output_dir.join("unitig_segments.fasta");
    write_segment_fasta(&segments_fasta, &graph.segments)?;
    let paf_path = request.output_dir.join("unitig_to_polished.paf");
    let stderr_path = request.output_dir.join("unitig_to_polished.minimap2.log");
    run_unitig_projection(
        request.soft_paths,
        request.reference,
        &segments_fasta,
        &paf_path,
        &stderr_path,
        request.threads,
    )?;
    let alignments = read_paf(&paf_path)?;
    let mut projections = boundaries
        .iter()
        .filter_map(|boundary| project_boundary(boundary, &alignments, &graph))
        .collect::<Vec<_>>();
    annotate_edge_status(&mut projections, &graph);
    let (problem_scope, suspect_segments, guidance) = classify_problem(&projections, &boundaries);
    let report = request.output_dir.join("sv_graph_localization.tsv");
    write_localization_report(
        &report,
        request,
        &boundaries,
        &projections,
        &problem_scope,
        &suspect_segments,
        &guidance,
    )?;
    Ok(SvGraphLocalization {
        report,
        problem_scope,
        suspect_segments,
        guidance,
    })
}

pub(crate) fn write_unavailable_localization(
    output_dir: &Path,
    subgroup: &str,
    unitig_gfa: Option<&Path>,
    reason: &str,
) -> Result<SvGraphLocalization, OrgraftError> {
    fs::create_dir_all(output_dir)?;
    let report = output_dir.join("sv_graph_localization.tsv");
    let guidance = format!("unitig graph localization unavailable: {reason}");
    let mut out = File::create(&report)?;
    writeln!(out, "# metric\tvalue")?;
    writeln!(out, "# status\tunavailable")?;
    writeln!(out, "# subgroup\t{subgroup}")?;
    writeln!(
        out,
        "# unitig_graph\t{}",
        unitig_gfa
            .map(|path| path.display().to_string())
            .unwrap_or_else(|| ".".to_string())
    )?;
    writeln!(out, "# problem_scope\tunavailable")?;
    writeln!(out, "# guidance\t{guidance}")?;
    writeln!(
        out,
        "pair_index\trole\treference_position\tsegment\tapprox_segment_offset\tsegment_length\tsegment_side\tstrand\tprojection_distance\tcoverage\tincident_links\tedge_status"
    )?;
    Ok(SvGraphLocalization {
        report,
        problem_scope: "unavailable".to_string(),
        suspect_segments: Vec::new(),
        guidance,
    })
}

fn read_unitig_graph(path: &Path) -> Result<UnitigGraph, OrgraftError> {
    let reader = BufReader::new(File::open(path).map_err(|error| {
        OrgraftError::InvalidArgument(format!(
            "could not read unitig graph {}: {error}",
            path.display()
        ))
    })?);
    let mut segments = Vec::new();
    let mut raw_links = Vec::new();
    for line_result in reader.lines() {
        let line = line_result?;
        let fields = line.split('\t').collect::<Vec<_>>();
        match fields.first().copied() {
            Some("S") if fields.len() >= 3 && fields[2] != "*" => {
                let coverage = fields[3..].iter().find_map(|field| {
                    field
                        .strip_prefix("DP:f:")
                        .or_else(|| field.strip_prefix("RC:f:"))
                        .and_then(|value| value.parse::<f64>().ok())
                });
                segments.push(GfaSegment {
                    name: fields[1].to_string(),
                    sequence: fields[2].to_string(),
                    coverage,
                });
            }
            Some("L") if fields.len() >= 5 => raw_links.push((
                fields[1].to_string(),
                fields[2].chars().next().unwrap_or('+'),
                fields[3].to_string(),
                fields[4].chars().next().unwrap_or('+'),
            )),
            _ => {}
        }
    }
    if segments.is_empty() {
        return Err(OrgraftError::InvalidArgument(format!(
            "{} has no sequence-bearing S records",
            path.display()
        )));
    }
    let segment_by_name = segments
        .iter()
        .cloned()
        .map(|segment| (segment.name.clone(), segment))
        .collect::<HashMap<_, _>>();
    let mut links = HashSet::new();
    let mut incident: HashMap<PhysicalEnd, BTreeSet<String>> = HashMap::new();
    for (from, from_orient, to, to_orient) in raw_links {
        let from_end = (from, if from_orient == '+' { 'R' } else { 'L' });
        let to_end = (to, if to_orient == '+' { 'L' } else { 'R' });
        let physical = canonical_physical_link(from_end.clone(), to_end.clone());
        if !links.insert(physical) {
            continue;
        }
        incident
            .entry(from_end.clone())
            .or_default()
            .insert(format_physical_end(&to_end));
        incident
            .entry(to_end)
            .or_default()
            .insert(format_physical_end(&from_end));
    }
    Ok(UnitigGraph {
        segments,
        segment_by_name,
        links,
        incident,
    })
}

fn parse_boundary_key(value: &str) -> Result<Vec<BoundaryPosition>, OrgraftError> {
    let mut positions = Vec::new();
    for (pair_index, pair) in value.split(';').enumerate() {
        let mut count = 0usize;
        for item in pair.split(',') {
            let (role, position) = item.split_once('=').ok_or_else(|| {
                OrgraftError::InvalidArgument(format!(
                    "invalid SV boundary item `{item}` in `{value}`"
                ))
            })?;
            let position = position.parse::<usize>().map_err(|error| {
                OrgraftError::InvalidArgument(format!(
                    "invalid SV boundary position `{position}` in `{value}`: {error}"
                ))
            })?;
            positions.push(BoundaryPosition {
                pair_index: pair_index + 1,
                role: role.to_string(),
                position,
            });
            count += 1;
        }
        if count != 2 {
            return Err(OrgraftError::InvalidArgument(format!(
                "SV boundary pair `{pair}` must contain two positions"
            )));
        }
    }
    if positions.is_empty() {
        return Err(OrgraftError::InvalidArgument(
            "SV boundary key contains no positions".to_string(),
        ));
    }
    Ok(positions)
}

fn write_segment_fasta(path: &Path, segments: &[GfaSegment]) -> Result<(), OrgraftError> {
    let mut out = File::create(path)?;
    for segment in segments {
        writeln!(out, ">{}", segment.name)?;
        for chunk in segment.sequence.as_bytes().chunks(80) {
            out.write_all(chunk)?;
            out.write_all(b"\n")?;
        }
    }
    Ok(())
}

fn run_unitig_projection(
    soft_paths: &Path,
    reference: &Path,
    segments_fasta: &Path,
    paf_path: &Path,
    stderr_path: &Path,
    threads: usize,
) -> Result<(), OrgraftError> {
    let minimap2 = read_tool_path(soft_paths, "minimap2")?;
    let output = Command::new(&minimap2)
        .arg("-x")
        .arg("asm5")
        .arg("--secondary=no")
        .arg("-t")
        .arg(threads.max(1).to_string())
        .arg(reference)
        .arg(segments_fasta)
        .output()
        .map_err(|error| {
            OrgraftError::InvalidArgument(format!(
                "failed to run {} for SV graph localization: {error}",
                minimap2.display()
            ))
        })?;
    fs::write(paf_path, &output.stdout)?;
    fs::write(stderr_path, &output.stderr)?;
    if !output.status.success() {
        return Err(OrgraftError::InvalidArgument(format!(
            "minimap2 failed during SV graph localization; see {}",
            stderr_path.display()
        )));
    }
    Ok(())
}

fn read_tool_path(soft_paths: &Path, name: &str) -> Result<PathBuf, OrgraftError> {
    let reader = BufReader::new(File::open(soft_paths).map_err(|error| {
        OrgraftError::InvalidArgument(format!(
            "could not read {} for SV graph localization: {error}",
            soft_paths.display()
        ))
    })?);
    for line_result in reader.lines() {
        let line = line_result?;
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        let fields = trimmed.split_whitespace().collect::<Vec<_>>();
        if fields.len() >= 2 && fields[0] == name {
            return Ok(PathBuf::from(fields[1]));
        }
    }
    Err(OrgraftError::InvalidArgument(format!(
        "{} is missing required tool `{name}`",
        soft_paths.display()
    )))
}

fn read_paf(path: &Path) -> Result<Vec<PafAlignment>, OrgraftError> {
    let reader = BufReader::new(File::open(path)?);
    let mut alignments = Vec::new();
    for (line_number, line_result) in reader.lines().enumerate() {
        let line = line_result?;
        if line.trim().is_empty() {
            continue;
        }
        let fields = line.split('\t').collect::<Vec<_>>();
        if fields.len() < 12 {
            return Err(OrgraftError::InvalidArgument(format!(
                "{}:{} has fewer than 12 PAF fields",
                path.display(),
                line_number + 1
            )));
        }
        let parse = |index: usize, name: &str| -> Result<usize, OrgraftError> {
            fields[index].parse::<usize>().map_err(|error| {
                OrgraftError::InvalidArgument(format!(
                    "{}:{} has invalid {name}: {error}",
                    path.display(),
                    line_number + 1
                ))
            })
        };
        alignments.push(PafAlignment {
            query: fields[0].to_string(),
            query_len: parse(1, "query length")?,
            query_start: parse(2, "query start")?,
            query_end: parse(3, "query end")?,
            strand: fields[4].chars().next().unwrap_or('+'),
            target_len: parse(6, "target length")?,
            target_start: parse(7, "target start")?,
            target_end: parse(8, "target end")?,
            matches: parse(9, "matches")?,
            mapq: fields[11].parse::<u8>().unwrap_or(0),
        });
    }
    Ok(alignments)
}

fn project_boundary(
    boundary: &BoundaryPosition,
    alignments: &[PafAlignment],
    graph: &UnitigGraph,
) -> Option<Projection> {
    let mut candidates = alignments
        .iter()
        .filter_map(|alignment| {
            let target_position = boundary
                .position
                .saturating_sub(1)
                .min(alignment.target_len.saturating_sub(1));
            let distance = interval_distance(
                target_position,
                alignment.target_start,
                alignment.target_end,
            );
            (distance <= MAX_PROJECTION_DISTANCE).then_some((alignment, distance, target_position))
        })
        .collect::<Vec<_>>();
    candidates.sort_by(|left, right| {
        left.1
            .cmp(&right.1)
            .then_with(|| right.0.mapq.cmp(&left.0.mapq))
            .then_with(|| right.0.matches.cmp(&left.0.matches))
    });
    let (alignment, target_distance, target_position) = candidates.first().copied()?;
    let segment = graph.segment_by_name.get(&alignment.query)?;
    let segment_offset = projected_query_offset(alignment, target_position);
    let segment_side = classify_segment_side(segment_offset, alignment.query_len);
    Some(Projection {
        boundary: boundary.clone(),
        segment: alignment.query.clone(),
        segment_offset,
        segment_len: alignment.query_len,
        segment_side,
        strand: alignment.strand,
        target_distance,
        coverage: segment.coverage,
        incident_links: incident_links_for_segment(graph, &alignment.query),
        edge_status: ".".to_string(),
    })
}

fn interval_distance(position: usize, start: usize, end: usize) -> usize {
    if position < start {
        start - position
    } else if position >= end {
        position - end.saturating_sub(1)
    } else {
        0
    }
}

fn projected_query_offset(alignment: &PafAlignment, target_position: usize) -> usize {
    let clamped = target_position.clamp(
        alignment.target_start,
        alignment.target_end.saturating_sub(1),
    );
    let target_span = alignment
        .target_end
        .saturating_sub(alignment.target_start)
        .max(1);
    let query_span = alignment
        .query_end
        .saturating_sub(alignment.query_start)
        .max(1);
    let target_offset = clamped.saturating_sub(alignment.target_start);
    let query_delta = target_offset.saturating_mul(query_span) / target_span;
    let query_zero = if alignment.strand == '+' {
        alignment.query_start.saturating_add(query_delta)
    } else {
        alignment
            .query_end
            .saturating_sub(1)
            .saturating_sub(query_delta)
    };
    query_zero.min(alignment.query_len.saturating_sub(1)) + 1
}

fn classify_segment_side(offset: usize, length: usize) -> char {
    let slop = ENDPOINT_SLOP.min((length / 3).max(1));
    if offset <= slop {
        'L'
    } else if length.saturating_sub(offset) < slop {
        'R'
    } else {
        'I'
    }
}

fn annotate_edge_status(projections: &mut [Projection], graph: &UnitigGraph) {
    let mut by_pair: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
    for (index, projection) in projections.iter().enumerate() {
        by_pair
            .entry(projection.boundary.pair_index)
            .or_default()
            .push(index);
    }
    for indices in by_pair.values() {
        if indices.len() != 2 {
            continue;
        }
        let left = &projections[indices[0]];
        let right = &projections[indices[1]];
        let status = if left.segment_side == 'I' || right.segment_side == 'I' {
            "node_split_required"
        } else {
            let physical = canonical_physical_link(
                (left.segment.clone(), left.segment_side),
                (right.segment.clone(), right.segment_side),
            );
            if graph.links.contains(&physical) {
                "link_present"
            } else {
                "link_absent"
            }
        };
        for index in indices {
            projections[*index].edge_status = status.to_string();
        }
    }
}

fn classify_problem(
    projections: &[Projection],
    boundaries: &[BoundaryPosition],
) -> (String, Vec<String>, String) {
    let suspect_segments = projections
        .iter()
        .map(|projection| projection.segment.clone())
        .collect::<BTreeSet<_>>()
        .into_iter()
        .collect::<Vec<_>>();
    if projections.len() < boundaries.len() {
        return (
            "partially_unmapped".to_string(),
            suspect_segments,
            "not every SV breakpoint could be projected; inspect the PAF and GFA manually"
                .to_string(),
        );
    }
    let internal = projections
        .iter()
        .filter(|projection| projection.segment_side == 'I')
        .collect::<Vec<_>>();
    if !internal.is_empty() {
        let mut offsets: BTreeMap<String, BTreeSet<usize>> = BTreeMap::new();
        for projection in &internal {
            offsets
                .entry(projection.segment.clone())
                .or_default()
                .insert(projection.segment_offset);
        }
        let detail = offsets
            .iter()
            .map(|(segment, values)| format!("{}:{}", segment, format_offset_clusters(values)))
            .collect::<Vec<_>>()
            .join(";");
        let scope = if suspect_segments.len() == 1 {
            "intra_unitig"
        } else {
            "mixed_internal_and_link"
        };
        return (
            scope.to_string(),
            suspect_segments,
            format!(
                "breakpoints fall inside unitig sequence ({detail}); split the listed S record(s) near these approximate offsets before editing L links"
            ),
        );
    }
    let absent_pairs = projections
        .iter()
        .filter(|projection| projection.edge_status == "link_absent")
        .map(|projection| projection.boundary.pair_index)
        .collect::<BTreeSet<_>>();
    if !absent_pairs.is_empty() {
        return (
            "inter_unitig_missing_link".to_string(),
            suspect_segments,
            format!(
                "inspect or add the read-supported physical connection for boundary pair(s) {}; remove a competing path only after read-support review",
                absent_pairs
                    .iter()
                    .map(ToString::to_string)
                    .collect::<Vec<_>>()
                    .join(",")
            ),
        );
    }
    (
        "inter_unitig_existing_link".to_string(),
        suspect_segments,
        "the read-supported boundary links already exist in graph.gfa; inspect competing links and copy usage around the listed segment ends"
            .to_string(),
    )
}

fn incident_links_for_segment(graph: &UnitigGraph, segment: &str) -> String {
    ['L', 'R']
        .into_iter()
        .map(|side| {
            let values = graph
                .incident
                .get(&(segment.to_string(), side))
                .map(|neighbors| neighbors.iter().cloned().collect::<Vec<_>>().join(","))
                .filter(|value| !value.is_empty())
                .unwrap_or_else(|| ".".to_string());
            format!("{side}=[{values}]")
        })
        .collect::<Vec<_>>()
        .join(";")
}

fn write_localization_report(
    path: &Path,
    request: &SvGraphLocalizationRequest<'_>,
    boundaries: &[BoundaryPosition],
    projections: &[Projection],
    problem_scope: &str,
    suspect_segments: &[String],
    guidance: &str,
) -> Result<(), OrgraftError> {
    let mut out = File::create(path)?;
    writeln!(out, "# metric\tvalue")?;
    writeln!(out, "# status\tlocalized")?;
    writeln!(out, "# subgroup\t{}", request.subgroup)?;
    writeln!(out, "# boundary_key\t{}", request.boundary_key)?;
    writeln!(out, "# unitig_graph\t{}", request.unitig_gfa.display())?;
    writeln!(out, "# reference\t{}", request.reference.display())?;
    writeln!(
        out,
        "# projection_method\tminimap2 asm5; segment offsets are approximate near alignment indels"
    )?;
    writeln!(out, "# problem_scope\t{problem_scope}")?;
    writeln!(
        out,
        "# suspect_segments\t{}",
        if suspect_segments.is_empty() {
            ".".to_string()
        } else {
            suspect_segments.join(",")
        }
    )?;
    writeln!(out, "# guidance\t{guidance}")?;
    writeln!(
        out,
        "pair_index\trole\treference_position\tsegment\tapprox_segment_offset\tsegment_length\tsegment_side\tstrand\tprojection_distance\tcoverage\tincident_links\tedge_status"
    )?;
    let projection_by_role = projections
        .iter()
        .map(|projection| {
            (
                (
                    projection.boundary.pair_index,
                    projection.boundary.role.as_str(),
                ),
                projection,
            )
        })
        .collect::<HashMap<_, _>>();
    for boundary in boundaries {
        if let Some(projection) =
            projection_by_role.get(&(boundary.pair_index, boundary.role.as_str()))
        {
            writeln!(
                out,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                boundary.pair_index,
                boundary.role,
                boundary.position,
                projection.segment,
                projection.segment_offset,
                projection.segment_len,
                projection.segment_side,
                projection.strand,
                projection.target_distance,
                projection
                    .coverage
                    .map(|value| format!("{value:.3}"))
                    .unwrap_or_else(|| ".".to_string()),
                projection.incident_links,
                projection.edge_status
            )?;
        } else {
            writeln!(
                out,
                "{}\t{}\t{}\t.\t.\t.\t.\t.\t.\t.\t.\tunmapped",
                boundary.pair_index, boundary.role, boundary.position
            )?;
        }
    }
    Ok(())
}

fn canonical_physical_link(left: PhysicalEnd, right: PhysicalEnd) -> PhysicalLink {
    if left <= right {
        (left, right)
    } else {
        (right, left)
    }
}

fn format_offset_clusters(values: &BTreeSet<usize>) -> String {
    let mut clusters: Vec<(usize, usize)> = Vec::new();
    for value in values {
        if let Some((_, end)) = clusters.last_mut() {
            if value.saturating_sub(*end) <= 50 {
                *end = *value;
                continue;
            }
        }
        clusters.push((*value, *value));
    }
    clusters
        .into_iter()
        .map(|(start, end)| {
            if start == end {
                start.to_string()
            } else {
                format!("{start}-{end}")
            }
        })
        .collect::<Vec<_>>()
        .join(",")
}

fn format_physical_end(endpoint: &PhysicalEnd) -> String {
    format!("{}:{}", endpoint.0, endpoint.1)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_multi_breakpoint_key() {
        let positions = parse_boundary_key("se1=1,ss2=299704;se2=299355,ss3=340100").unwrap();
        assert_eq!(positions.len(), 4);
        assert_eq!(positions[0].pair_index, 1);
        assert_eq!(positions[0].role, "se1");
        assert_eq!(positions[3].position, 340100);
    }

    #[test]
    fn projects_plus_and_minus_alignments() {
        let plus = PafAlignment {
            query: "utg1".to_string(),
            query_len: 1_000,
            query_start: 0,
            query_end: 1_000,
            strand: '+',
            target_len: 2_000,
            target_start: 100,
            target_end: 1_100,
            matches: 1_000,
            mapq: 60,
        };
        let mut minus = plus.clone();
        minus.strand = '-';
        assert_eq!(projected_query_offset(&plus, 100), 1);
        assert_eq!(projected_query_offset(&plus, 1_099), 1_000);
        assert_eq!(projected_query_offset(&minus, 100), 1_000);
        assert_eq!(projected_query_offset(&minus, 1_099), 1);
    }

    #[test]
    fn internal_breakpoints_require_unitig_split() {
        let boundaries = parse_boundary_key("se1=10,ss2=20").unwrap();
        let projections = boundaries
            .iter()
            .map(|boundary| Projection {
                boundary: boundary.clone(),
                segment: "utg7".to_string(),
                segment_offset: if boundary.role == "se1" { 8_000 } else { 8_350 },
                segment_len: 58_337,
                segment_side: 'I',
                strand: '+',
                target_distance: 0,
                coverage: Some(58.88),
                incident_links: "L=[utg0:R];R=[utg5:R]".to_string(),
                edge_status: "node_split_required".to_string(),
            })
            .collect::<Vec<_>>();
        let (scope, suspects, guidance) = classify_problem(&projections, &boundaries);
        assert_eq!(scope, "intra_unitig");
        assert_eq!(suspects, vec!["utg7"]);
        assert!(guidance.contains("utg7:8000,8350"));
    }

    #[test]
    fn nearby_offsets_are_reported_as_one_split_interval() {
        let offsets = BTreeSet::from([8_010, 8_360, 48_808, 48_813]);
        assert_eq!(format_offset_clusters(&offsets), "8010,8360,48808-48813");
    }
}
