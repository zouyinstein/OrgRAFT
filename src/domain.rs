#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ToolRole {
    pub source: &'static str,
    pub orgraft_role: &'static str,
    pub boundary: &'static str,
}

pub const TOOL_ROLES: &[ToolRole] = &[
    ToolRole {
        source: "simple_draft_asm",
        orgraft_role: "raw draft graph construction",
        boundary: "Rust algorithm core",
    },
    ToolRole {
        source: "OrgRAFT",
        orgraft_role: "read recruitment, polish evaluation, and compact evidence summaries",
        boundary: "Rust orchestration with external aligner and variant-caller calls",
    },
    ToolRole {
        source: "GFA_Editor",
        orgraft_role: "GUI review, edit history, auto graph operations, and export",
        boundary: "GUI decision layer plus CLI execution layer",
    },
];

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct EvidenceFile {
    pub filename: &'static str,
    pub producer: &'static str,
    pub purpose: &'static str,
}

pub const EVIDENCE_FILES: &[EvidenceFile] = &[
    EvidenceFile {
        filename: "organelle_reads.fastq.gz",
        producer: "recruit",
        purpose: "enriched read set for assembly and validation",
    },
    EvidenceFile {
        filename: "read_classification.tsv",
        producer: "recruit",
        purpose: "per-read organelle/nuclear/ambiguous assignment",
    },
    EvidenceFile {
        filename: "recruitment_summary.tsv",
        producer: "recruit",
        purpose: "read recruitment counts and coverage summary",
    },
    EvidenceFile {
        filename: "graph.gfa",
        producer: "asm",
        purpose: "raw draft assembly graph",
    },
    EvidenceFile {
        filename: "checked_draft.gfa",
        producer: "GFA_Editor GUI",
        purpose: "manual checkpoint 1 graph after low-support edit decisions",
    },
    EvidenceFile {
        filename: "topology_summary.tsv",
        producer: "resolve",
        purpose: "concise topology summary for checked draft GFA",
    },
    EvidenceFile {
        filename: "reference_rotation.tsv",
        producer: "resolve",
        purpose: "reversible rotation record for each linear reference",
    },
    EvidenceFile {
        filename: "auto_repeat_check.tsv",
        producer: "resolve",
        purpose: "GFA_Editor CLI auto-repeat and auto-merge summary",
    },
    EvidenceFile {
        filename: "linearized_subgraph.fasta",
        producer: "resolve",
        purpose: "draft subgraph FASTA aligned and oriented to the rotated reference",
    },
    EvidenceFile {
        filename: "subgraph_reads.fastq.gz",
        producer: "resolve",
        purpose: "reads binned to one linearized subgraph",
    },
    EvidenceFile {
        filename: "polished_aligned.fasta",
        producer: "polish",
        purpose: "polished FASTA aligned to the rotated linear reference",
    },
    EvidenceFile {
        filename: "all_sorted_blastn_alignments.txt",
        producer: "polish",
        purpose: "SV alignment evidence for polished FASTA evaluation",
    },
    EvidenceFile {
        filename: "all_bcftools_calls.txt",
        producer: "polish",
        purpose: "SNP/InDel evidence before compact annotation",
    },
    EvidenceFile {
        filename: "variants_summary.tsv",
        producer: "polish",
        purpose: "compact SNP/InDel annotation and correction candidates",
    },
    EvidenceFile {
        filename: "pos_ref_alt.txt",
        producer: "polish",
        purpose: "accepted or candidate reference correction records",
    },
    EvidenceFile {
        filename: "verified.gfa",
        producer: "rebuild",
        purpose: "final verified graph model",
    },
    EvidenceFile {
        filename: "verified.fasta",
        producer: "rebuild",
        purpose: "final verified sequence model",
    },
    EvidenceFile {
        filename: "graph_coordinate_map.tsv",
        producer: "rebuild",
        purpose: "final linear-to-graph coordinate bridge",
    },
    EvidenceFile {
        filename: "node_coverage.tsv",
        producer: "rebuild",
        purpose: "node-level coverage in the verified graph",
    },
    EvidenceFile {
        filename: "repeat_path_support.tsv",
        producer: "rebuild",
        purpose: "repeat path evidence in the verified graph",
    },
    EvidenceFile {
        filename: "final_report.html",
        producer: "rebuild",
        purpose: "compact handoff report",
    },
];

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn records_three_source_tool_roles() {
        assert_eq!(TOOL_ROLES.len(), 3);
    }

    #[test]
    fn includes_core_validation_tables() {
        let names: Vec<&str> = EVIDENCE_FILES.iter().map(|file| file.filename).collect();
        assert!(names.contains(&"topology_summary.tsv"));
        assert!(names.contains(&"all_sorted_blastn_alignments.txt"));
        assert!(names.contains(&"all_bcftools_calls.txt"));
        assert!(names.contains(&"repeat_path_support.tsv"));
    }
}
