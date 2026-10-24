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
orgraft postrepair
orgraft rebuild
orgraft refine
```

Use `orgraft <command> --help` for the live interface contract.

## Standalone postrepair

Use the same ordinary validation round that would enter the workflow repair stage:

```bash
orgraft postrepair --source-dir 04.polish/mito/subgraph_001/round_2 --out-dir 05.postrepair --soft-paths soft_paths.txt --threads 8
```

The source report supplies the component and original input paths. Original resolve
provenance supplies topology; if it is outside the standard sibling directory, pass
`--resolve-details` and `--checked-gfa`. The source stays read-only. Unresolved
candidates remain review requirements. Workflow verifies completed output before
reuse. Direct standalone execution refuses an existing repair stage; use a fresh
output directory after inspecting an interrupted run.

## Workflow Commands

`orgraft workflow` coordinates project folders, manual checkpoints, and
validation rounds. It does not replace the core algorithms in
`recruit/asm/resolve/polish/postrepair/rebuild/refine`; it connects them from one
config. For multi-case configs, the generated run script contains the commands
directly. Recruitment runs once. Components with identical draft and reference
inputs share assembly and resolve after their checked graphs are compared.
Each component then passes its own validation and postrepair checks before
rebuild. Refine runs once for all components of an organelle.

```bash
orgraft workflow template
orgraft workflow init --sample S1 --out results_workflow/orgraft.workflow.toml
orgraft workflow plan --config results_workflow/orgraft.workflow.toml
orgraft workflow run-script --config results_workflow/orgraft.workflow.toml
orgraft workflow run --config results_workflow/orgraft.workflow.toml --case mito_subgraph_001
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
orgraft workflow init --sample S1 --out results_workflow/orgraft.workflow.toml
```

Minimal shape:

```toml
[project]
sample = "S1"
results_dir = "."
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

[commands.postrepair]
enabled = true
out_dir = "${results_dir}/05.postrepair"
frequency = 0.2
threads = 16

[commands.rebuild]
enabled = true
out_dir = "${results_dir}/06.rebuild"
threads = 16

[commands.refine]
enabled = true
out_dir = "${results_dir}/07.refine"
threads = 16
# Required when cases use different split read pools:
# organelle_reads = "${results_dir}/01.recruit/${organelle}.fastq.gz"

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
rebuild_out_dir = "${results_dir}/06.rebuild/${organelle}"
# For the recorded-path rebuild, use the corresponding unresolved component:
rebuild_edited_gfa = "${results_dir}/03.resolve_gfa/${organelle}/graph/merged_unresolved_${subgraph}.gfa"
```

All relative paths, including `results_dir`, reads, references, and
`soft_paths`, are relative to the directory containing the TOML. Use
`results_dir = "."` when the config is in the sample's output folder. Absolute
paths remain absolute. `workflow init --out S1/orgraft.workflow.toml` infers `S1`
when `--sample` is omitted; explicitly set the sample for a route folder such as
`samples/S1/software/curated`.

`workflow plan` writes `run_S1.sh` and `run_S1_manual.sh`. The main script directly
contains each standard command and its checkpoint checks, without calling a
chain of generated shell scripts. The manual script shows one validation round
per case; stop at its checkpoint when another correction round or review is
needed. To run a command separately, retain its displayed input assignments or
replace shell variables with their actual paths. The script header selects the
OrgRAFT binary; neither entrypoint needs an AI service or development script.

New configs enable `commands.postrepair` and `commands.refine`. Existing configs
that omit them keep their previous behavior, including explicit `05.rebuild`
paths and `commands.polish.homorepair`. New postrepair uses only the completed
ordinary round as input and writes to its own root. It preserves unresolved
status and must pass its sequence/evidence checks before rebuild.

Refine runs once per organelle after all enabled components in that organelle
have completed rebuild. `workflow refine-bundle --config FILE --organelle mito`
prepares the bundle from the actual rebuild outputs and verifies their source,
validation, and read bindings. Missing or unreviewed components stop this step.
For components with different split read pools, set
`commands.refine.organelle_reads` to the complete organelle read pool. The path
may contain `${organelle}`. Refine refuses to use the first component's split
pool for the whole organelle.

For completed standalone native rebuilds, use:

```bash
orgraft refine --rebuild-dir 06.rebuild/mito --reads 01.recruit/mito.fastq.gz --sample S1 --organelle mito --out-dir 07.refine/mito --auto
```

This verifies all component reports, graphs, FASTAs, and their original validation
evidence before adding `refine.inputs.tsv` and `refine_source_<component>.tsv`
to `06.rebuild/mito/logs`. Existing rebuild products stay unchanged. Supply the
complete organelle read pool with `--reads`; component read pools do not replace
it. Keep `--out-dir` separate from the rebuild and source directories.


Repeat downgrade preserves recorded segment boundaries by default. To use verified
nonbranching chain compaction, explicitly pass
`--downgrade-representation linear-compaction`. In a workflow, set
`refine_downgrade_representation = "linear-compaction"` in every enabled
`[workflow.case.NAME]` for the chosen organelle; its other organelle can retain
`"recorded"`. Mixed settings within one organelle are rejected. The setting is
bound to the run inputs, so changing it cannot reuse an older completion record.
The candidate, read-support and native validation rules are unchanged.

An explicit `commands.refine.sample_bundle` can supply a reviewed bundle;
`${organelle}` in that path is expanded for each organelle. `refine-check`
requires `confirmed` or `no_change`; it does not turn `review_required` or
`partial_requires_review` into success.

When `reference` is omitted, workflow uses the matching FASTA from
`commands.recruit.baits`. Multiple `[workflow.case.NAME]` sections support
multi-organelle or multi-subgraph runs. With refine enabled, use one sample per
config so components from different samples cannot enter the same bundle. Old
configs without refine can still describe multiple samples. A case with
`enabled = false` stays in the config but is omitted from the generated master
script unless selected explicitly with `--case NAME`. When every case is disabled,
both generated entrypoints explain that state and exit with a nonzero status
without running analysis.

## Reusing Existing Results

To rerun resolve and later stages from an existing draft and recruited reads,
keep the old inputs under `source_results_dir` and use a separate `results_dir`.
Apply these changes to the generated TOML; retain its references and other
required settings:

```toml
[project]
results_dir = "rerun"
source_results_dir = "results"

[workflow]
force = false

[commands.recruit]
enabled = false

[commands.asm]
enabled = false

[commands.resolve]
out_dir = "${results_dir}/03.resolve_gfa"

[commands.polish]
out_dir = "${results_dir}/04.polish"

[commands.refine]
organelle_reads = "${source_results_dir}/01.recruit/${organelle}.fastq.gz"

[workflow.case.mito_subgraph_001]
workflow_dir = "${results_dir}/workflow/${organelle}/${subgraph}"
draft_graph = "${source_results_dir}/02.draft_asm/${organelle}/03.finalize_graph/graph.gfa"
checked_draft_gfa = "${results_dir}/workflow/${organelle}/${subgraph}/checkpoint_1/checked_draft.gfa"
reads = "${source_results_dir}/01.recruit/${organelle}.fastq.gz"
rebuild_out_dir = "${results_dir}/06.rebuild/${organelle}"
```

Use the actual paths in the old run when its layout differs.
`source_results_dir` is a path substitution, not a read-only switch: keep every
output setting and case override under the new `results_dir`, including
postrepair and refine. Checkpoint files belong in the new `workflow_dir`.

To start at 05.postrepair from a completed native validation round, add
`postrepair_source_dir` to the case. `workflow plan` generates main and manual
scripts that run postrepair, rebuild, and refine directly; `workflow run --case`
uses the same entrypoint. Recruit and asm must be disabled. No checkpoint 1 or 2,
resolve, or polish command is run for that case, and no old checkpoint status is
copied or marked complete.

```toml
[workflow.case.mito_subgraph_001]
organelle = "mito"
subgraph = "subgraph_001"
postrepair_source_dir = "/old-run/04.polish/mito/subgraph_001/round_2"
reads = "/fixed-inputs/mito.fastq.gz"
checked_draft_gfa = "/old-run/workflow/mito/subgraph_001/checkpoint_1/checked_draft.gfa"
resolve_out_dir = "/old-run/03.resolve_gfa"
rebuild_edited_gfa = "/old-run/03.resolve_gfa/mito/graph/merged_unresolved_subgraph_001.gfa"
# Optional overrides for original provenance stored elsewhere:
# resolve_details = "/old-run/03.resolve_gfa/mito/logs/resolve_details.tsv"
# first_polish_report = "/old-run/04.polish/mito/subgraph_001/round_1/logs/report.tsv"
workflow_dir = "${results_dir}/workflow/${organelle}/${subgraph}"
rebuild_out_dir = "${results_dir}/06.rebuild/${organelle}"
```

Enable `commands.postrepair`, `commands.rebuild`, and `commands.refine`, and keep
their outputs under the new `results_dir`. The source must retain its native
report, per-read calls, variant table, SV-pass evidence, reference, reads, and
validated FASTA. The case identity and source round are checked before scripts
are written; unresolved ordinary corrections are rejected. The source round's
reference defines downstream coordinates. A configured `polish_reference` must
match it exactly. Postrepair verifies that the configured reads are the same
cohort as the source screening.

Keep the original checked graph, resolve record, and first polish report as
read-only evidence. Circularity and repair decisions still require their normal
native evidence; missing provenance is not replaced with a guessed topology.
For split components, `commands.refine.organelle_reads` must name the full
organelle pool. A review-required postrepair result stops the workflow; starting
from an old round does not approve it. `postrepair-check` defaults to the round
named by `postrepair_source_dir`, and an explicit different `--round` is rejected.

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
  orgraft.workflow.toml
  run_S1.sh
  run_S1_manual.sh
  runtime_summary.md             # written by workflow runtime-summary
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
  05.postrepair/
    mito/subgraph_001/
      postrepair.md
      postrepair.tsv
      round_hr/
      round_N+1/                 # only when the sequence changes
  06.rebuild/
    mito/
      subgraph_001/
        rebuild_subgraph_001.gfa
        rebuild_subgraph_001.fasta
        rebuild_subgraph_001_nodes.fasta
      logs/
  07.refine/
    mito/
      01_assess/
      auto_complete.tsv           # --auto completion and retained bundle
  workflow/
    mito/
      refine.inputs.tsv
    mito/subgraph_001/
      refine_source.tsv
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

For standalone `postrepair --source-dir`, open `postrepair.md` in the component
folder for the final status, selected sequence, and evidence links. Its ordinary
review section records that stage separately when a later repair changes the
sequence. Native probe tables and final candidate evidence remain available;
the redundant per-probe `homorepair_review.tsv` snapshot is omitted. The legacy
`--source-round` layout retains its existing `summary.md` and probe snapshots.

The generated shell scripts are reproducible runners. Edit the TOML config,
then regenerate scripts with `orgraft workflow plan --config FILE --force`. For
`plan`, this replaces only the generated main and manual scripts; it does not run
analysis or enable forced replacement in them. Without `--force`, an identical
plan is reused and a changed script is refused. Do not hand-edit generated scripts
as the source of truth. `--out path/run_S1.sh` puts the main script at that path
and `run_S1_manual.sh` in the same parent directory.

## Individual Stage Synopsis

These are separate invocations, not an unattended script. Use the reviewed
checked graph and completed validation round. Replace `SELECTED_FASTA` with the
selected sequence path reported by `postrepair.md`, keeping its native validation
workspace available. For a special repair, also retain the evidence arguments
shown by the generated manual entrypoint. Standalone relative paths use the
current shell directory.

```bash
orgraft recruit --reads reads.fastq.gz --mito mito.fa --plastid plastid.fa --out-dir 01.recruit
orgraft asm --reads 01.recruit/mito.fastq.gz --organelle mito --out-dir 02.draft_asm
orgraft resolve --checked-draft-gfa checked_draft.gfa --reference mito.fa --organelle mito --out-dir 03.resolve_gfa
orgraft polish --organelle mito --subgraph subgraph_001 --draft 03.resolve_gfa/mito/fasta/resolved_subgraphs.fasta --reference 03.resolve_gfa/mito/fasta/rotated_reference.fasta --reads 01.recruit/mito.fastq.gz --out-dir 04.polish
orgraft postrepair --source-dir 04.polish/mito/subgraph_001/round_2 --out-dir 05.postrepair
orgraft rebuild --organelle mito --subgraph subgraph_001 --edited-gfa 03.resolve_gfa/mito/graph/merged_unresolved_subgraph_001.gfa --checked-gfa checked_draft.gfa --resolve-details 03.resolve_gfa/mito/logs/resolve_details.tsv --polished-fasta SELECTED_FASTA --reads 01.recruit/mito.fastq.gz --out-dir 06.rebuild/mito
orgraft refine --rebuild-dir 06.rebuild/mito --reads 01.recruit/mito.fastq.gz --sample S1 --organelle mito --out-dir 07.refine/mito --auto
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

A human may retain the unchanged source at one unresolved ordinary InDel with
`orgraft postrepair review --result-dir 05.postrepair/mito/subgraph_001 --keep-source --pos 194368 --ref CA --alt CAAACAAAGAGA --reason "Reviewed this event; retain source"`.
This writes only `manual_review.tsv`, bound to the current source and evidence;
the automatic review status stays unchanged. `workflow postrepair-check` and the
generated scripts verify that record before continuing. `--keep-source` requires
an unchanged source and refuses adopted edits or unrelated pending reviews.

To retain the current selected sequence after accepted edits, use
`postrepair review --result-dir DIR --keep-current --reason TEXT` and explicitly
list all remaining events with repeated `--retain-indel STEP:POS:REF:ALT`,
`--retain-circle START:END:MOTIF:COPIES`, or
`--retain-junction LEFT:RIGHT:BASE:CURRENT` options. InDel coordinates belong to
the recorded native step; circle coordinates are 1-based and inclusive on the
pre-circle sequence. Junction endpoints and CURRENT use the original graph
orientation. A junction decision cannot be combined with other event types.
This records a current human choice; it does not turn unresolved evidence into
an automatic pass. All accepted edits and their native evidence must still
replay, and the full selected sequence must match its validation. See
`orgraft postrepair review --help` for the complete argument contract.

To replay an already approved postrepair decision from a clean run, declare a
human-readable TSV per organelle. This is an input decision, not a copied receipt:

```toml
[commands.postrepair.review]
mito = "inputs/mito.postrepair.tsv"
```

The TSV has exactly five columns: `component_id`, `decision`, `event_type`,
`event`, and `reason`. Each row names one event. Rows for a component must share
one decision and one nonempty reason. Only enabled, unambiguous components may
be listed; unlisted components retain the automatic workflow.

| decision | event_type | event syntax |
| --- | --- | --- |
| `keep_source` | `indel` | `POS:REF:ALT` (one event) |
| `keep_source` | `repeat_extra` | `START:END:MOTIF:COPIES` (one or more) |
| `keep_current` | `indel` | `STEP:POS:REF:ALT` |
| `keep_current` | `circle` | `START:END:MOTIF:COPIES` |
| `keep_current` | `junction` | `LEFT:RIGHT:BASE:CURRENT` (one event, no other types) |
| `correct_snv` | `correct_snv` | `POS:REF:ALT` (exactly one) |
| `correct_snv` | `retain_indel` | `POS:REF:ALT` (exactly one) |
| `correct_snv` | `evidence_file` | optional supplementary file path |

The two `keep_source` event types cannot be mixed. Coordinate conventions and
all restrictions are identical to `postrepair review`; duplicate events and
conflicting declarations are rejected. TSV and attachment paths resolve from
the TOML directory. Manual correction uses the configured postrepair threads
(or polish threads) and software paths. It runs the existing native validation;
supplementary attachments are not recomputed or treated as automatic approval.

For a listed component, the generated main/manual scripts use
`orgraft workflow postrepair-run --config orgraft.workflow.toml --case CASE --round N`.
This shares the execution path used by `workflow run`: ordinary postrepair runs
first, and only its explicit terminal review can proceed to the declared native
`postrepair review` command. Other execution errors remain errors. An existing
review result must pass the same complete native replay before a receipt is
created. An existing receipt must match the declared mode, every event, reason,
and correction inputs; its original evidence must still verify. Automatic
status and evidence stay unchanged. This interface needs neither a batch JSON
nor an external orchestration script.

An explicit refine decision can be supplied per organelle in the TOML:

```toml
[commands.refine.adopt]
mito = "inputs/mito.refine.tsv"
```

The decision TSV uses `component_id`, `canonical_repeat_id`, `decision`, and
`decision_note` columns. For an approved intra-repeat upgrade, the decision is
`adopt_upgrade_intra`. Only listed targets receive that explicit choice;
unspecified targets retain automatic rules. Paths resolve from the configuration
file. The file must exist before `workflow plan` regenerates the main and manual
scripts. Both scripts and `workflow run` pass it to native `refine --auto --adopt`.
If a previous refine run used different inputs, preserve that run before starting
the new one; an existing output is not silently replaced.


Intra-repeat measurement defaults to `exact-core`. To explicitly measure p2/p3
at the original complete-copy boundaries for unequal copies in one organelle,
use the following setting (native equivalent: `--intra-measurement full-copy`):

```toml
[commands.refine.intra_measurement]
mito = "full-copy"
```

The main script, manual script and `workflow run` use the same setting. It is
recorded in the automatic input signature and checked by `workflow refine-check`.
It does not change copy identity, graph surgery, SV tolerances or adoption rules.
An omitted setting preserves the existing exact-core measurement. Historical
completion records without the method field are interpreted as exact-core;
explicit full-copy requires its matching recorded method and completion digest.

The duplicate-template full-span method is off by default. Enable it explicitly
for an approved organelle (native equivalent: `--duplicate-template-full-span`):

```toml
[commands.refine.duplicate_template_full_span]
mito = true
```

Only when enabled, if the four default graph-path templates contain identical
complete DNA, refine may derive a complete family span from the current graph and its bound
source walk. This path requires both exact family copies, all four physical
ports and four distinct complete extended template sequences. The query library
is then mapped again; existing counts are not transferred to the new templates.
The existing length-band path and its proof format remain unchanged. The explicit
choice is carried through child assessments, source-proof replay and the automatic
input signature. A historical run without this flag may retain its default-method
status only if its bound template tables did not use the duplicate-template branch.
Preserve an affected old run before starting its replacement; no result is overwritten.

For an explicitly approved graph with a short segmented flank, refine can follow
source-bound single-copy, single-in/single-out `0M` chains to obtain its existing
2000 bp templates and full 1000 bp required outer anchors:

```toml
[commands.refine.linear_chain_flanks]
plastid = true
```

The native option is `--linear-chain-flanks`; it is off by default. The option
applies to initial, step and final assessments and is bound to completion. A
branch, cycle, another repeat, unknown source or insufficient chain stops the
assessment. It does not compact the graph, reuse counts or change support gates.
Combining it with a complete-family span extension is rejected. Shared rebuild
commands retain their default adjacent-node template behavior.

To measure graph repeats over the complete repeat family, enable the option
for the selected organelle (native equivalent: `--graph-full-family`):

```toml
[commands.refine.graph_full_family]
plastid = true
```

This option is off by default. It derives the complete repeat-family boundaries
from the current graph, its bound source sequence and the native family evidence.
Both complete copies, their placement and all four ports must be established;
the four complete templates must be distinct and have the existing 2000 bp outer
flanks and 1000 bp required anchors. Missing or ambiguous evidence stops the new
measurement instead of falling back to core-only templates. The native read
mapping and counting rules, SV checks and adoption thresholds remain unchanged.

The main script, manual script and `workflow run` share this setting. Child and
final assessments carry it, and `workflow refine-check` verifies the recorded
method and completion binding. A historical run without the flag remains a
default-mode result and cannot satisfy an explicit full-family request. Preserve
an old output before rerunning 07 with a different method. Do not combine this
option with `duplicate_template_full_span` or `linear_chain_flanks` for the same
organelle. It does not change `intra_measurement` or standalone rebuild behavior.


To measure the maximal exact graph span while retaining the complete repeat
family as the policy identity, select the explicit measurement mode:

```toml
[commands.refine.graph_measurement]
mito = "maximal-exact-span"
```

The native option is `--graph-measurement maximal-exact-span`; the default remains
`default`. The mode derives the span from the current graph and its fixed native
source, then checks it against the unique exact block in the retained complete
HSP. Both full copies must have equal length, no alignment gaps and at least 99%
identity. Their complete family length, copy coordinates and actual DNA remain
bound separately from the shorter exact measurement span. Four distinct complete
templates and the existing 2000 bp flanks and 1000 bp required outer anchors are
mandatory. Missing, ambiguous or inconsistent evidence stops the assessment.

The complete family continues to determine representation policy and thresholds.
Mapping, counting, SV validation, post-operation target checks and non-target
checks are unchanged. A downgrade still requires its complete-family
post-operation measurement; it cannot reuse pre-operation graph counts. This
option does not change graph DNA or standalone rebuild behavior and does not
predict the resulting counts or automatic decision.

To measure both complete observed family alleles, use the value
`complete-family-alleles` in the same TOML field, or pass
`--graph-measurement complete-family-alleles`. This mode requires two complete,
equal-length, ungapped HSP copies with at least 99% identity, bound to the current
native graph and fixed source walk. The exact shared core and span remain a
separate proof. Both complete copy sequences, their coordinates, all four ports,
full template DNA and required windows are recorded and replayed.

Every template's entire family sequence must equal one of the two observed
copies. A combination of left and right variants that creates an unobserved
allele is rejected. Four distinct templates, the full 2000 bp flanks and 1000 bp
outer anchors are required. Complete-family policy, the fixed read pool and all
post-operation target and non-target checks remain unchanged. The option is off
by default and does not predict the resulting counts or automatic decision.

Main, manual, direct workflow, child and final assessments carry the same mode.
Completion and `workflow refine-check` bind and replay the mode and its actual
source inputs. Preserve an old 07 output before changing the measurement mode.
Do not combine this option with `graph_full_family`, `duplicate_template_full_span`
or `linear_chain_flanks` for the same organelle. `intra_measurement` is independent.

To retain every physical continuation when a repeat has short, branching flanks,
use `--graph-measurement multi-context`, or set `mito = "multi-context"` under
`[commands.refine.graph_measurement]`. This explicit mode requires a native
source-bound complete graph walk, all physical continuations, and full 2000 bp
template flanks with the existing 1000 bp required outer anchors. It maps the
whole component together, counts a read once per physical port pair, and keeps
the existing cross-port ambiguity rule. If any read reaches the mapper’s
secondary-alignment limit, the whole context measurement requires review; no
read is silently dropped to obtain a smaller count. Complete-family identity and policy
thresholds remain unchanged; an independent exact graph unit does not upgrade
a lower-identity HSP into a qualified family. For a native same-component
downgrade, each assessment maps the complete original component query library
again to the same bound read pool. Resolved targets remain in that competing
library. Source contracts must replay the parent audit, both target copies,
compaction, and every non-target core and physical port. A removed continuation
must be explained by the target's resolved port pairing and is recorded as
absent from the candidate; its query hits are not candidate-path support.
Unexplained changes, new paths, merges, changed family identity, or changed
query windows require review. Target and non-target acceptance checks remain
in force. This option cannot be combined with the other graph-template flags.
The default mode is unchanged.

For an explicitly selected `multi-context` or `multi-context-maximal-exact-span`
assessment, a verified current source
may contain two exact complete-family paths represented by shared fragments.
Before building the query library, refine can normalize one uniquely eligible
region per component: it must be acyclic, have exactly two equal-DNA internal
paths, two entry and two exit ports, and no other internal branch. Every source
occurrence, physical edge, complete family copy and P record remains bound to
the original source. The derived graph and fixed walk are saved under
`01_assess/00_source`; consumers replay the same proof without writing files.
The original bundle, FASTA, reads and validation inputs remain the execution
inputs. No SV or old graph counts are copied to the derived graph, and all
measurement and adoption gates still run. Derived or pending child sources are
not normalized again. An existing manual declaration whose original identity
is absent is rejected; it is never silently removed or reassigned. Other graph
measurement modes retain their previous behavior.

To combine fixed original component queries with maximal exact spans for original
single four-path graph repeats, use `--graph-measurement multi-context-maximal-exact-span`
or set `mito = "multi-context-maximal-exact-span"` (or `plastid`) under
`[commands.refine.graph_measurement]`. Multi-path repeats retain their original
context spans; descendants inherit that choice. A selected maximal span must
pass the existing complete-family HSP, exact-copy, four-template, 2000-bp flank
and 1000-bp anchor checks. Failure requires review rather than falling back to
a different span. The entire original component query library is measured
fresh; target and non-target adoption gates remain unchanged.

## Version and Reproducibility

The current package version is 0.5.0. Record the complete source commit, package
version, compiled executable checksum, external-tool versions, TOML, and explicit
review inputs when starting a new analysis. `Cargo.lock` pins the Rust package
resolution; `soft_paths.txt` records host-specific tool locations.

```bash
git rev-parse HEAD
cargo build --release --locked
./target/release/orgraft --version
```

Version 0.5.0 was prepared after the recorded rerun audit. Its scientific Rust
sources and required test fixtures are identical to the audited 0.4.2 source
snapshot, recorded in commit `e40358ecdc01af41eed0fd5080f364e0a07bfe25`.
The 0.5.0 update changes package metadata and user documentation. Historical
outputs retain the executable and source identities recorded by their own
commands, logs, and receipts; a later package version does not relabel those runs.

A rerun from a fixed read pool must preserve that input and its provenance.
This is distinct from repeating the original recruitment from whole-genome
reads. Report the original recruitment version, parameters, and seed only when
they were recorded. A current default is not evidence of a historical setting.
The reference and explicit review decisions remain declared inputs; reference-
assisted agreement is not an independent accuracy benchmark.
