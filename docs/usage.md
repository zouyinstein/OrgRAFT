# OrgRAFT Usage Guide

This is the detailed operating guide. Keep the README short; put workflow,
configuration, checkpoint, and output details here.

## Setup First

Check external tools before generating or running workflow scripts:

```bash
orgraft setup --soft-paths soft_paths.txt --requirements requirements.txt
```

`soft_paths.txt` should provide executable names or paths for `python`,
`minimap2`, `blastn`, `pigz`, and optionally `gfa_editor_cli`.

## Command Map

```text
orgraft setup
orgraft workflow
orgraft recruit
orgraft asm
orgraft resolve
orgraft polish
orgraft rebuild
```

Use `orgraft <command> --help` for the live interface contract.

## Workflow Commands

`orgraft workflow` coordinates project folders, manual checkpoints, and
validation rounds. It does not replace the core algorithms in
`recruit/asm/resolve/polish/rebuild`; it strings them together from one config.
For multi-case configs, the generated master script runs the shared recruit
stage once, then calls case scripts for draft assembly through rebuild.

```bash
orgraft workflow template
orgraft workflow init --out results_workflow/orgraft.workflow.toml
orgraft workflow plan --config results_workflow/orgraft.workflow.toml
orgraft workflow run-script --config results_workflow/orgraft.workflow.toml
orgraft workflow run --config results_workflow/orgraft.workflow.toml
orgraft workflow runtime-summary --config results_workflow/orgraft.workflow.toml --force
```

Checkpoint and correction commands:

```bash
orgraft workflow checkpoint1 --config results_workflow/orgraft.workflow.toml --case mito_subgraph_001
orgraft workflow checkpoint2 --config results_workflow/orgraft.workflow.toml --case mito_subgraph_001 --round 1
orgraft workflow checkpoint2 --config results_workflow/orgraft.workflow.toml --case mito_subgraph_001 --round 1 --sv-subgroup type_3_subtype_rep_rep_NA:4
orgraft workflow correct --input-fasta INPUT.fasta --pos-ref-alt pos_ref_alt.txt --out OUTPUT.fasta
orgraft workflow test-correction --input-fasta INPUT.fasta --pos-ref-alt pos_ref_alt.txt --out-dir results_workflow/correction_test
orgraft workflow test-fake-validate --input-fasta ERROR_POLISH_ALN.fasta --pos-ref-alt pos_ref_alt.txt --out-dir results_workflow/fake_validate --force
```

## Workflow Config

The workflow file is TOML. Generate a starter file instead of hand-writing it:

```bash
orgraft workflow init --out results_workflow/orgraft.workflow.toml
```

Minimal shape:

```toml
[project]
sample = "sample_001"
results_dir = "results_workflow"
# Optional when reusing existing outputs as read-only inputs:
# source_results_dir = "results"

[software]
soft_paths = "soft_paths.txt"

[workflow]
mode = "stepwise"
max_rounds = 3
auto_sv_correction = true
threads = 64
force = false
auto_snv_indel_correction = true
topology_simple_allowed_classes = "0-0,0-1/1-0,1-1,1-2/2-1,2-2"

[commands.recruit]
enabled = true
raw_reads = "/path/to/raw_hifi.fastq.gz"
baits = "mito=/path/to/mito.fasta,plastid=/path/to/plastid.fasta"
out_dir = "${results_dir}/01.recruit"
threads = 16
# platform = "HiFi"
# bait_format = "auto"
# gzip_output = true
# max_reads = "all,20000"

[commands.asm]
enabled = true
out_dir = "${results_dir}/02.draft_asm"
threads = 8
# profile = "standard"

[commands.resolve]
out_dir = "${results_dir}/03.resolve_gfa"
# gfa_editor_mode = "rust"
# max_states = 5000
# max_candidates = 100

[commands.polish]
out_dir = "${results_dir}/04.polish"
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

[commands.rebuild]
enabled = true
out_dir = "${results_dir}/05.rebuild"
threads = 16

[workflow.case.mito_subgraph_001]
enabled = true
name = "mito_subgraph_001"
organelle = "mito"
subgraph = "subgraph_001"
workflow_dir = "${results_dir}/workflow/${organelle}/${subgraph}"
draft_graph = "${results_dir}/02.draft_asm/${organelle}/03.finalize_graph/graph.gfa"
checked_draft_gfa = "${results_dir}/workflow/${organelle}/${subgraph}/checkpoint_1/checked_draft.gfa"
# resolve_out_dir defaults to commands.resolve.out_dir
reads = "${results_dir}/01.recruit/${organelle}.fastq.gz"
# linearized_fasta defaults to commands.resolve.out_dir/${organelle}/fasta/resolved_subgraphs.fasta
# polish_reference defaults to commands.resolve.out_dir/${organelle}/fasta/rotated_reference.fasta
# polish_out_dir defaults to commands.polish.out_dir
# sv_correction_subgroup = "type_3_subtype_rep_rep_NA:4"
rebuild_out_dir = "${results_dir}/05.rebuild/${organelle}"
```

When `reference` is omitted, workflow uses the matching FASTA from
`commands.recruit.baits`. Multiple `[workflow.case.NAME]` sections support
multi-sample, multi-organelle, or multi-subgraph runs. A case with
`enabled = false` stays in the config but is omitted from the generated master
script unless selected explicitly with `--case NAME`.

## Reusing Existing Results

For a clean workflow run that reads old results but writes new outputs, set
`results_dir = "results_workflow"` and point selected inputs at `results`.

Example:

```toml
[project]
results_dir = "results_workflow"
source_results_dir = "results"

[commands.recruit]
enabled = false

[commands.asm]
enabled = false

[commands.resolve]
out_dir = "${results_dir}/03.resolve_gfa"

[commands.polish]
out_dir = "${results_dir}/04.polish"

[workflow.case.mito_subgraph_001]
draft_graph = "${source_results_dir}/draft_asm/${organelle}/03.finalize_graph/graph.gfa"
reads = "${source_results_dir}/recruit/${organelle}.fastq.gz"
rebuild_out_dir = "${results_dir}/05.rebuild/${organelle}"
```

This keeps existing `results/` as a reference input and writes new checkpoint,
resolve, polish, and rebuild products under `results_workflow/`.

## Manual Checkpoints

Checkpoint 1 reads the draft GFA. If the graph topology is simple and all links
reference declared segment records, it writes `checked_draft.gfa` and marks the
checkpoint as checked. If the graph is complex or inconsistent, it writes a
manual-edit copy and stops.

Checkpoint 2 reads polish validation evidence. If SV support fails, it stops
for manual inspection unless `sv_high_subgroups.tsv` contains an automatically
repairable `possible_reference_sv_error`. It repairs one subgroup at a time,
requires most selected reads to become `type_1`, and accepts a candidate only
when `low_green_window_fraction` improves without reducing global reference
support. Set case-level `sv_correction_subgroup = "group_name:old_index"` or
pass `--sv-subgroup` to choose the subgroup manually.

Repeat-pairing mismatch is a separate checkpoint 2 path. When two alternative
flank pairings dominate the complete reads spanning the pairing currently used
by the linear FASTA, workflow records the matching repeat node and returns
`rebuild_ready`. Rebuild measures complete flank-repeat-flank support for all
four paths and uses the dominant perfect pairing to filter complete graph-valid
single-circle candidates. A unique survivor is selected directly; if more than
one survives, the existing global k-mer-chain score against the current polished
FASTA breaks the tie. No candidate or a tied top score requires manual review.
The unresolved rebuild GFA keeps all four `P` records and real read counts, while
the constrained rebuild FASTA is validated with all reads in the next round.

Before changing sequence, SV check projects subgroup breakpoints back to
`02.anchor_graph_core/02.unitig_graph/graph.gfa` and writes
`sv_repair/sv_graph_localization.tsv`. It distinguishes breakpoints inside one
`S` record (split the unitig before editing links) from missing or competing
connections between unitig ends. The graph path is derived automatically;
set case-level `unitig_graph = "/path/to/02.unitig_graph/graph.gfa"` only for a
nonstandard layout. Projection is advisory and never modifies the GFA.

An accepted subgroup-level SV repair or repeat-pairing rebuild adds one full-read
validation round without consuming the ordinary `max_rounds` budget. The
repeat-pairing round is validate-only: its `02.polish` directory remains empty
and `01.inputs/linear_subgraph.round_N.fasta` is the constrained rebuild FASTA.
The hard total across the workflow remains 10. Once no repairable SV subgroup
remains, checkpoint 2
applies SNV/InDel correction under the ordinary budget. If both validations
pass, the case is complete.

## Generated Layout

Typical workflow outputs:

```text
results_workflow/
  workflow.commands.sh
  runtime_summary.md
  01.recruit/
  02.draft_asm/
  03.resolve_gfa/
  04.polish/
    mito/subgraph_001/
      round_1/
        01.inputs/
        02.polish/
          polished_aln.fasta
        03.validate/
          01.data/
          02.plots/
            plot_sv_bubble.py
            bubble_type_2_rep_raw.png
            bubble_type_2_ins_raw.png
          03.reports/
      round_2/
        01.inputs/
          linear_subgraph.round_2.fasta
        02.polish/
        03.validate/
  05.rebuild/
  workflow/
    mito/subgraph_001/
      workflow.commands.sh
      checkpoint_1/
        checkpoint_1.status.tsv
        checked_draft.gfa
        manual_edit_required.gfa
      checkpoint_2/
        round_1/
          checkpoint_2.status.tsv
          repeat_pairing_repair/
            repeat_pairing_mismatch.tsv
          pos_ref_alt.txt
          polish_aln_v2.fasta
```

The generated shell scripts are reproducible runners. Edit the TOML config,
then regenerate scripts with `orgraft workflow plan`; do not hand-edit generated
scripts as the source of truth.

## Individual Stage Synopsis

```bash
orgraft recruit --reads reads.fastq.gz --mito mito.fa --plastid plastid.fa
orgraft asm --reads results_workflow/01.recruit/mito.fastq.gz --organelle mito
orgraft resolve --checked-draft-gfa checked_draft.gfa --reference mito.fa
orgraft polish --organelle mito --subgraph subgraph_001 --draft resolved_subgraphs.fasta --reference rotated_reference.fasta --reads mito.fastq.gz
orgraft rebuild --organelle mito --subgraph subgraph_001 --edited-gfa checked_draft.gfa --polished-fasta polish_aln_v2.fasta
```

### Standalone read-backed GFA repair

Use `--skeleton-gfa FILE --stable` when a candidate or manually edited graph
should enter read-backed open-end link repair without rebuilding anchor walks
and unitigs from reads:

```bash
orgraft asm \
  --reads results_workflow/01.recruit/plastid.fastq.gz \
  --organelle plastid \
  --skeleton-gfa graph.edited.gfa \
  --stable \
  --soft-paths soft_paths.txt \
  --out-dir results/gfa_repair_r01
```

`--out-dir` is a new assembly root; the public result is written to
`results/gfa_repair_r01/plastid/03.finalize_graph/graph.gfa`. Do not point it at
an existing `03.finalize_graph` directory. Repair mode requires the standard
profile and rejects `--subsets`, `--min-graph-coverage`, and `--branch-ratio`.

The supplied GFA replaces de novo Steps 01-02. Reads remain required as remapping
evidence: Steps 01, 03, and 04 are reported as skipped, Step 02 records the
provided GFA, Step 05 recalculates skeleton depth/link evidence, Step 06 evaluates
read-backed links between open endpoints, and Steps 07-08 publish the handoff and
provenance. Automatic triplet/copy-choice, internal-split, repeat-expansion, and
redundant-link-pruning rewrites are disabled in the production Stable route.

Treat depth and link support cautiously. They depend on the supplied node
boundaries, minimap2 filters, and the exact recruitment/splitting/subsetting
history of the evidence reads; reads reused from graph construction are not an
independent validation set. Input links without `RC:i` support must earn enough
read-remapping support to remain. The repair reader accepts GFA1 `S`/`L` topology
with unique segment IDs, inline sequences, valid `+/-` orientations, and declared
link endpoints. The repaired graph is reconstructed and does not preserve `P`,
`W`, or arbitrary custom tags. Run checkpoint 1 again before `resolve`.

Standalone command defaults remain command-local, such as `results/recruit`,
`results/draft_asm`, `resolve_gfa`, `results/polish`, and `results/rebuild`.
The workflow template overrides them with numbered roots under `results_dir`.

## Correction Smoke Test

`test-fake-validate` simulates a swapped `polish_aln.fasta` that contains known
variants, triggers checkpoint2 correction, and then checks that the next round
can finish cleanly.

```bash
orgraft workflow test-fake-validate --input-fasta ERROR_POLISH_ALN.fasta --pos-ref-alt pos_ref_alt.txt --out-dir results_workflow/fake_validate --force
```

The fake validate inputs are explicit by design. Keep local or project-specific
test data paths outside the Rust source tree and pass them through CLI options.
