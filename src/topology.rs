use std::collections::BTreeMap;
use std::io::BufRead;

use crate::error::OrgraftError;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Endpoint {
    Left,
    Right,
}

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct EndpointDegrees {
    pub left: usize,
    pub right: usize,
    pub self_links: usize,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct NodeTaxon {
    pub code: &'static str,
    pub name: &'static str,
    pub interpretation: &'static str,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct NodeTopology {
    pub node_id: String,
    pub degrees: EndpointDegrees,
    pub taxon: NodeTaxon,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TopologyReport {
    pub node_count: usize,
    pub link_count: usize,
    pub nodes: Vec<NodeTopology>,
}

pub const TAXONOMY: &[NodeTaxon] = &[
    NodeTaxon {
        code: "0-0",
        name: "linear or isolated node",
        interpretation: "Chlamydomonas-like linear node or unlinked contig",
    },
    NodeTaxon {
        code: "0-1/1-0",
        name: "open node",
        interpretation: "open endpoint evidence for a non-closed conformation",
    },
    NodeTaxon {
        code: "1-1",
        name: "mergeable connection or self loop",
        interpretation: "basic operation class and a handle for uncycling",
    },
    NodeTaxon {
        code: "1-2/2-1",
        name: "branch node",
        interpretation: "can create non-circular topology; some cases are evidence-resolvable",
    },
    NodeTaxon {
        code: "1-2 self/2-1 self",
        name: "self-associated branch",
        interpretation: "special non-circular type handled by topology-aware review",
    },
    NodeTaxon {
        code: "2-2",
        name: "two-in two-out bridge node",
        interpretation: "foundation for automated graph resolving and repeat bridge handling",
    },
    NodeTaxon {
        code: "other",
        name: "higher-order complex node",
        interpretation: "can be locally decomposed or reduced to the basic classes",
    },
];

pub fn classify_node(degrees: EndpointDegrees) -> NodeTaxon {
    match (degrees.left, degrees.right, degrees.self_links > 0) {
        (0, 0, _) => TAXONOMY[0],
        (0, 1, _) | (1, 0, _) => TAXONOMY[1],
        (1, 1, _) => TAXONOMY[2],
        (1, 2, true) | (2, 1, true) => TAXONOMY[4],
        (1, 2, false) | (2, 1, false) => TAXONOMY[3],
        (2, 2, _) => TAXONOMY[5],
        _ => TAXONOMY[6],
    }
}

pub fn analyze_gfa<R>(reader: R) -> Result<TopologyReport, OrgraftError>
where
    R: BufRead,
{
    let mut nodes: BTreeMap<String, EndpointDegrees> = BTreeMap::new();
    let mut link_count = 0usize;

    for (line_index, line_result) in reader.lines().enumerate() {
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
                nodes.entry((*segment_id).to_string()).or_default();
            }
            Some("L") => {
                let from = fields.get(1).ok_or_else(|| {
                    OrgraftError::InvalidArgument(format!(
                        "GFA line {line_number}: link record is missing from segment"
                    ))
                })?;
                let from_orientation = fields.get(2).ok_or_else(|| {
                    OrgraftError::InvalidArgument(format!(
                        "GFA line {line_number}: link record is missing from orientation"
                    ))
                })?;
                let to = fields.get(3).ok_or_else(|| {
                    OrgraftError::InvalidArgument(format!(
                        "GFA line {line_number}: link record is missing to segment"
                    ))
                })?;
                let to_orientation = fields.get(4).ok_or_else(|| {
                    OrgraftError::InvalidArgument(format!(
                        "GFA line {line_number}: link record is missing to orientation"
                    ))
                })?;

                let from_endpoint = oriented_exit_endpoint(from_orientation, line_number)?;
                let to_endpoint = oriented_entry_endpoint(to_orientation, line_number)?;

                increment_endpoint(nodes.entry((*from).to_string()).or_default(), from_endpoint);
                increment_endpoint(nodes.entry((*to).to_string()).or_default(), to_endpoint);

                if from == to {
                    nodes.entry((*from).to_string()).or_default().self_links += 1;
                }

                link_count += 1;
            }
            _ => {}
        }
    }

    let nodes = nodes
        .into_iter()
        .map(|(node_id, degrees)| NodeTopology {
            node_id,
            degrees,
            taxon: classify_node(degrees),
        })
        .collect::<Vec<_>>();

    Ok(TopologyReport {
        node_count: nodes.len(),
        link_count,
        nodes,
    })
}

pub fn nodes_tsv(report: &TopologyReport) -> String {
    let mut output = String::from(
        "node_id\tleft_degree\tright_degree\tself_links\tclass_code\tclass_name\tinterpretation\n",
    );

    for node in &report.nodes {
        output.push_str(&format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
            node.node_id,
            node.degrees.left,
            node.degrees.right,
            node.degrees.self_links,
            node.taxon.code,
            node.taxon.name,
            node.taxon.interpretation
        ));
    }

    output
}

pub fn summary_tsv(report: &TopologyReport) -> String {
    let mut counts: BTreeMap<&str, usize> = BTreeMap::new();
    for node in &report.nodes {
        *counts.entry(node.taxon.code).or_default() += 1;
    }

    let mut output = String::from("metric\tvalue\tnotes\n");
    output.push_str(&format!(
        "node_count\t{}\tnumber of segment records plus linked implicit nodes\n",
        report.node_count
    ));
    output.push_str(&format!(
        "link_count\t{}\tnumber of GFA L records\n",
        report.link_count
    ));
    for (class_code, count) in counts {
        output.push_str(&format!(
            "class:{class_code}\t{count}\tnode endpoint-degree class count\n"
        ));
    }
    output
}

fn oriented_exit_endpoint(orientation: &str, line_number: usize) -> Result<Endpoint, OrgraftError> {
    match orientation {
        "+" => Ok(Endpoint::Right),
        "-" => Ok(Endpoint::Left),
        _ => Err(OrgraftError::InvalidArgument(format!(
            "GFA line {line_number}: invalid orientation `{orientation}`"
        ))),
    }
}

fn oriented_entry_endpoint(
    orientation: &str,
    line_number: usize,
) -> Result<Endpoint, OrgraftError> {
    match orientation {
        "+" => Ok(Endpoint::Left),
        "-" => Ok(Endpoint::Right),
        _ => Err(OrgraftError::InvalidArgument(format!(
            "GFA line {line_number}: invalid orientation `{orientation}`"
        ))),
    }
}

fn increment_endpoint(degrees: &mut EndpointDegrees, endpoint: Endpoint) {
    match endpoint {
        Endpoint::Left => degrees.left += 1,
        Endpoint::Right => degrees.right += 1,
    }
}

#[cfg(test)]
mod tests {
    use std::io::Cursor;

    use super::*;

    #[test]
    fn classifies_basic_degree_patterns() {
        assert_eq!(
            classify_node(EndpointDegrees {
                left: 0,
                right: 0,
                self_links: 0
            })
            .code,
            "0-0"
        );
        assert_eq!(
            classify_node(EndpointDegrees {
                left: 1,
                right: 2,
                self_links: 0
            })
            .code,
            "1-2/2-1"
        );
        assert_eq!(
            classify_node(EndpointDegrees {
                left: 1,
                right: 2,
                self_links: 1
            })
            .code,
            "1-2 self/2-1 self"
        );
        assert_eq!(
            classify_node(EndpointDegrees {
                left: 2,
                right: 2,
                self_links: 0
            })
            .code,
            "2-2"
        );
    }

    #[test]
    fn analyzes_gfa_endpoint_degrees() {
        let gfa = "\
S\tA\tACGT
S\tB\tACGT
S\tC\tACGT
L\tA\t+\tB\t+\t0M
L\tB\t+\tC\t+\t0M
L\tB\t-\tC\t+\t0M
";
        let report = analyze_gfa(Cursor::new(gfa)).unwrap();
        assert_eq!(report.node_count, 3);
        assert_eq!(report.link_count, 3);

        let b = report
            .nodes
            .iter()
            .find(|node| node.node_id == "B")
            .unwrap();
        assert_eq!(b.degrees.left, 2);
        assert_eq!(b.degrees.right, 1);
        assert_eq!(b.taxon.code, "1-2/2-1");
    }
}
