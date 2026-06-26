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
threads = 8
force = false
auto_snv_indel_correction = true
topology_simple_allowed_classes = "0-0,0-1/1-0,1-1,1-2/2-1,2-2"

[commands.recruit]
enabled = true
raw_reads = "/path/to/raw_hifi.fastq.gz"
baits = "mito=/path/to/mito.fasta,plastid=/path/to/plastid.fasta"
out_dir = "${results_dir}/recruit"
platform = "HiFi"
bait_format = "auto"
gzip_output = true

[commands.asm]
enabled = true
out_dir = "${results_dir}/draft_asm"
profile = "standard"

[commands.rebuild]
enabled = true
out_dir = "${results_dir}/rebuild"
threads = 4

[workflow.case.mito_subgraph_001]
enabled = true
name = "mito_subgraph_001"
organelle = "mito"
subgraph = "subgraph_001"
draft_graph = "${results_dir}/draft_asm/${organelle}/03.finalize_graph/graph.gfa"
checked_draft_gfa = "${results_dir}/workflow/${organelle}/${subgraph}/checkpoint_1/checked_draft.gfa"
resolve_out_dir = "${results_dir}/resolve_gfa"
reads = "${results_dir}/recruit/${organelle}.fastq.gz"
polish_out_dir = "${results_dir}/polish"
rebuild_out_dir = "${results_dir}/rebuild/${organelle}"
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

[workflow.case.mito_subgraph_001]
draft_graph = "${source_results_dir}/draft_asm/${organelle}/03.finalize_graph/graph.gfa"
reads = "${source_results_dir}/recruit/${organelle}.fastq.gz"
resolve_out_dir = "${results_dir}/resolve_gfa"
polish_out_dir = "${results_dir}/polish"
```

This keeps existing `results/` as a reference input and writes new checkpoint,
resolve, polish, and rebuild products under `results_workflow/`.

## Manual Checkpoints

Checkpoint 1 reads the draft GFA. If the graph topology is simple and all links
reference declared segment records, it writes `checked_draft.gfa` and marks the
checkpoint as checked. If the graph is complex or inconsistent, it writes a
manual-edit copy and stops.

Checkpoint 2 reads polish validation evidence. If SV support fails, it stops
for manual inspection. If SV passes and SNV/InDel rows require correction, it
writes the next-round polish input until `max_rounds` is reached. If both SV
and SNV/InDel validation pass, the case is complete.

## Generated Layout

Typical workflow outputs:

```text
results_workflow/
  workflow.commands.sh
  runtime_summary.md
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
          pos_ref_alt.txt
          polish_aln_v2.fasta
  recruit/
  draft_asm/
  resolve_gfa/
  polish/
  rebuild/
```

The generated shell scripts are reproducible runners. Edit the TOML config,
then regenerate scripts with `orgraft workflow plan`; do not hand-edit generated
scripts as the source of truth.

## Individual Stage Synopsis

```bash
orgraft recruit --reads reads.fastq.gz --mito mito.fa --plastid plastid.fa
orgraft asm --reads results_workflow/recruit/mito.fastq.gz --organelle mito
orgraft resolve --checked-draft-gfa checked_draft.gfa --reference mito.fa
orgraft polish --organelle mito --subgraph subgraph_001 --draft resolved_subgraphs.fasta --reference rotated_reference.fasta --reads mito.fastq.gz
orgraft rebuild --organelle mito --subgraph subgraph_001 --edited-gfa checked_draft.gfa --polished-fasta polish_aln_v2.fasta
```

Default output roots are stage-oriented: `results/recruit`,
`results/draft_asm`, `resolve_gfa`, `results/polish`, and `results/rebuild`,
unless overridden by workflow config or command options.

## Correction Smoke Test

`test-fake-validate` simulates a swapped `polish_aln.fasta` that contains known
variants, triggers checkpoint2 correction, and then checks that the next round
can finish cleanly.

```bash
orgraft workflow test-fake-validate --input-fasta ERROR_POLISH_ALN.fasta --pos-ref-alt pos_ref_alt.txt --out-dir results_workflow/fake_validate --force
```

The fake validate inputs are explicit by design. Keep local or project-specific
test data paths outside the Rust source tree and pass them through CLI options.
