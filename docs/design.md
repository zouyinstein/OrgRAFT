# OrgRAFT Design Notes

This is the detailed design and file-contract document. Keep quick usage in the
README and operational details in `docs/usage.md`.

## Product Thesis

The reliable product of plant organelle genome assembly is not a FASTA sequence
alone. It is a curated graph/sequence model accompanied by read support, graph
evidence, variant evidence, and reproducible edit history.

OrgRAFT should be assembler-aware but not assembler-dependent. It can generate
a conservative draft graph itself, while keeping the interfaces simple enough
to audit graphs and FASTA files produced by other plant organelle assembly
tools.

## Component Boundaries

| Component | Role | Boundary |
| --- | --- | --- |
| OrgRAFT Rust CLI | Recruit reads, assemble draft graphs, resolve checked GFAs, polish, rebuild verified GFAs, and write evidence tables | Rust algorithm core plus delegated aligner calls |
| Workflow layer | Parse TOML, generate scripts, manage checkpoints, and bound correction rounds | Orchestration only; core algorithms stay in subcommands |
| GFA_Editor CLI | Optional reference-coloured GFA PDF/SVG export | Visualization/export layer consuming stable GFA/FASTA files |
| External tools | `minimap2`, `blastn`, `pigz`, Python plotting helpers | Runtime dependencies recorded in `soft_paths.txt` and `requirements.txt` |

## Pipeline Layers

```text
OrgRAFT Rust CLI
  setup       external tool and Python package checks
  workflow    config template, command generation, checkpoints, summaries
  recruit     read recruitment
  asm         conservative graph construction
  resolve     checked graph orientation, conservative merge, subgraph export
  polish      polished linear sequence and variant/read-support validation
  rebuild     final verified graph/FASTA reconstruction and compact reports

GFA_Editor
  topology view
  evidence overlays
  variant markers
  unsupported-link flags
  edit history
  curated GFA/FASTA export
```

## File Contract

Paths are stage-oriented and should change deliberately, because downstream
inspection, reporting, and optional GFA_Editor export depend on them.

| Stage | Core outputs | Notes |
| --- | --- | --- |
| `recruit` | `mito.fastq.gz`, `plastid.fastq.gz`, `logs/bait.fasta`, `logs/recruitment_summary.tsv`, `logs/read_stats.tsv` | `pigz` may be used for fast gzip IO |
| `asm` | `ORGANELLE/03.finalize_graph/graph.gfa`, optional `graph.edited.gfa`, `ORGANELLE/logs/*.tsv` | The finalize graph removes reverse-complement duplicate `L` records by default while preserving intermediate evidence |
| `workflow checkpoint1` | `checkpoint_1.status.tsv`, `checked_draft.gfa`, optional `manual_edit_required.gfa` | Topology and GFA consistency gate before resolve |
| `resolve` | `ORGANELLE/graph/merged_unresolved.gfa`, `merged_unresolved_subgraph_001.gfa`, `ORGANELLE/fasta/rotated_reference.fasta`, `resolved_subgraphs.fasta`, `logs/resolve_details.tsv` | Uses checked draft graph plus reference FASTA |
| `polish` | `ORGANELLE/SUBGRAPH/02.polish/polished_aln.fasta`, `03.validate/round_1/01.data/*.tsv`, `03.validate/round_1/02.plots/*.png`, `logs/report.tsv` | Plot scripts require `matplotlib` |
| `workflow checkpoint2` | `checkpoint_2/round_N/checkpoint_2.status.tsv`, optional `pos_ref_alt.txt`, optional `polish_aln_v{N+1}.fasta` | SV failures stop for manual review; SNV/InDel correction can continue up to `max_rounds` |
| `rebuild` | `OUT/SUBGRAPH/rebuild_SUBGRAPH.gfa`, `rebuild_SUBGRAPH.fasta`, `rebuild_SUBGRAPH_nodes.fasta`, `OUT/logs/*.tsv` | PDF/SVG export is attempted only when an image reference is supplied |

## Rebuild Outputs

`rebuild_SUBGRAPH.gfa` is the compact verified graph. Its S records carry
rebuilt node sequences and restored depth/path support tags when evidence is
available.

`rebuild_SUBGRAPH.fasta` is the complete verified/polished linear sequence for
the subgraph.

`rebuild_SUBGRAPH_nodes.fasta` contains one FASTA record per S node in the
rebuilt GFA. It is intended for node-level comparison, depth checks, and
sequence consistency checks.

Optional graph images:

```text
OUT/SUBGRAPH/rebuild_SUBGRAPH.pdf
OUT/SUBGRAPH/rebuild_SUBGRAPH.svg
```

These are exported through `gfa_editor_cli` when `--image-reference-fasta FILE`
is supplied. Image export failures are recorded in the rebuild run report
without failing core GFA/FASTA output.

Important rebuild tables:

| Table | Shape | Purpose |
| --- | --- | --- |
| `rebuild_SUBGRAPH_extract.tsv` | one row per verified node copy projected onto polished coordinates | Node extraction and sequence projection evidence |
| `rebuild_SUBGRAPH_run_report.tsv` | `section`, `key`, `value` | Run metadata, tool status, input/output paths, image-export status |
| `rebuild_SUBGRAPH_result_stats.tsv` | `section`, `subgraph`, `item`, `metric`, `value`, `extra` | Summary, graph, consistency, depth, and repeat-path metrics |

## Topology Assumption

OrgRAFT treats the primary assembly object as an evidence-constrained graph
model rather than a pre-imposed linear FASTA. After orientation normalization,
topology simplification, and evidence-based filtering of weak links, complex
structures should be reducible to a finite set of node classes.

Endpoint degree is interpreted as left-in / right-out connectivity after
normalization.

| Class | Interpretation | Assembly meaning |
| --- | --- | --- |
| `0-0` | Linear or isolated node | Unlinked contig or linear node |
| `0-1` / `1-0` | Open node | Evidence for an open end and non-closed conformation |
| `1-1` | Ordinary connection or self loop | Basic operation class and writing handle for uncycling |
| `1-2` / `2-1` | Branch node | Can create non-circular topology; may be resolved by evidence |
| `1-2 self` / `2-1 self` | Self-associated branch | Special non-circular branch type |
| `2-2` | Two-in two-out bridge node | Foundation for automated resolving and repeat bridge handling |
| `other` | Higher-order complex node | Requires decomposition, reduction, or manual review |

Checkpoint 1 currently treats higher-order complex nodes and self-associated
branches as manual-review cases unless the workflow config explicitly changes
the allowed simple classes.

## Development Direction

1. Keep command names, evidence files, and report schemas stable.
2. Make the full workflow reproducible from TOML plus generated scripts.
3. Keep durable algorithmic work in Rust.
4. Delegate mature primitives such as alignment, compression, and optional
   graph visualization until replacing them creates clear value.
5. Validate with simulated/toy graphs, public plant HiFi datasets, and complex
   biological cases.
