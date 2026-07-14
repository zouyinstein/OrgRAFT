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
Workflow-generated roots are numbered for the major products:
`01.recruit`, `02.draft_asm`, `03.resolve_gfa`, `04.polish`, and
`05.rebuild`; checkpoint scripts and status files live under `workflow/`.

| Stage | Core outputs | Notes |
| --- | --- | --- |
| `recruit` | `mito.fastq.gz`, `plastid.fastq.gz`, `logs/bait.fasta`, `logs/recruitment_summary.tsv`, `logs/read_stats.tsv` | `pigz` may be used for fast gzip IO |
| `asm` | `ORGANELLE/03.finalize_graph/graph.gfa`, optional `graph.edited.gfa`, `ORGANELLE/logs/*.tsv` | The finalize graph removes reverse-complement duplicate `L` records by default while preserving intermediate evidence |
| `workflow checkpoint1` | `checkpoint_1.status.tsv`, `checked_draft.gfa`, optional `manual_edit_required.gfa` | Topology and GFA consistency gate before resolve |
| `resolve` | `ORGANELLE/graph/merged_unresolved.gfa`, `merged_unresolved_subgraph_001.gfa`, `ORGANELLE/fasta/rotated_reference.fasta`, `resolved_subgraphs.fasta`, `logs/resolve_details.tsv` | Uses checked draft graph plus reference FASTA |
| `polish` | `ORGANELLE/SUBGRAPH/round_N/{01.inputs,02.polish,03.validate,logs}`; round 2+ leaves `02.polish` empty and validates `01.inputs/linear_subgraph.round_N.fasta` | Plot scripts require `matplotlib` |
| `workflow checkpoint2` | `checkpoint_2/round_N/checkpoint_2.status.tsv`, optional `sv_repair/{sv_repair.tsv,sv_candidate_scores.tsv,sv_graph_localization.tsv}`, optional `repeat_pairing_repair/repeat_pairing_mismatch.tsv`, optional `pos_ref_alt.txt`, optional `polish_aln_v{N+1}.fasta` | Prepares one correction at a time; a repeat-pairing mismatch emits `rebuild_ready` and preserves resolve output |
| `rebuild` | `OUT/SUBGRAPH/rebuild_SUBGRAPH.gfa`, `rebuild_SUBGRAPH.fasta`, `rebuild_SUBGRAPH_nodes.fasta`, `OUT/logs/*.tsv` | Workflow passes reads so repeat `P` records carry complete flank-repeat-flank read counts; optional repeat constraint filters graph-valid circular candidates before writing FASTA |

## Repeat-Pairing Constraint

Checkpoint 2 only records the repeat implicated by the validation mismatch.
Rebuild retains the unresolved graph and writes all four real spanning-read path
counts as `P` records. It then uses the uniquely dominant perfect pairing as a
hard filter over complete graph-valid circular candidates. One remaining
candidate is selected directly; multiple candidates retain the existing global
k-mer-chain scoring against the current polished FASTA. No match or a top-score
tie stops for manual review. The constrained FASTA is passed to one additional
validate-only round, while the rebuild GFA remains unresolved and auditable.
The first constrained rebuild's candidate ID is retained as `source_candidate`;
`local_candidate` records any renumbering after later sequence correction.

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
| `rebuild_SUBGRAPH_repeat_path_support.tsv` | repeat path, endpoints, support count, ratio, read IDs | Audit trail for `PM:Z:spanning_reads` and integer `RC:i` path tags |
| `rebuild_SUBGRAPH_repeat_resolution.tsv` | dominant pairing, filtered candidates, source/local candidate IDs, orientation/rotation | Audit trail for optional constrained FASTA generation |

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
