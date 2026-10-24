# OrgRAFT

OrgRAFT means **Organelle Read-backed Assembly with FASTA–Graph Traceability**.

This repository contains a Rust CLI for plant organelle assembly and
post-assembly validation. The intended product is not only a final FASTA: it is
a curated graph/sequence model with read support, graph evidence, variant
evidence, and reproducible edit history.

## Concise Guide

Use this README as the short entry point. The detailed material is consolidated
into two files:

- [docs/usage.md](docs/usage.md): command usage, workflow config, checkpoints,
  output layout, and smoke tests.
- [docs/design.md](docs/design.md): architecture, file contracts, rebuild
  outputs, topology assumptions, and external-tool boundaries.

## Commands

```text
-- Project setup
orgraft setup       # check external software paths and Python packages
orgraft workflow    # generate and run config-driven workflow checkpoints

-- Raw graph generation
orgraft recruit     # organelle HiFi read recruitment
orgraft asm         # conservative draft assembly or read-backed GFA repair

-- High-quality graph generation
orgraft resolve     # resolve checked draft GFA into reference-oriented products
orgraft polish      # polish linearized graph FASTA and evaluate variants
orgraft postrepair  # review residual repeat, InDel, and circular-origin errors
orgraft rebuild     # rebuild final verified graph and compact reports
orgraft refine      # refine repeat representation and component connectivity
```

Run any command without arguments, or with `--help`, to print its current
interface.

`postrepair` starts from a completed validation round and infers its organelle,
component, FASTA, reference, and reads from the native report:

```bash
orgraft postrepair --source-dir 04.polish/mito/subgraph_001/round_2 --out-dir 05.postrepair --soft-paths soft_paths.txt --threads 8
```

The command uses the recorded graph and validation evidence to select applicable
repeat, InDel, and circular-origin checks. It writes a separate repair stage;
changes require fresh read validation. Unresolved evidence remains `review_required`
and stops the handoff to rebuild. It does not need an AI service or a separate
development script. Workflow checks completed output before reusing it. A direct
standalone rerun needs a fresh output directory; it does not replace an existing
repair stage.

The `compressed`, `joint-ir`, and `junction` subcommands expose specific evidence
checks for explicit replay. Compressed-reference decisions retain their separate
label. They do not replace the original strict read-support rules.

`refine` provides repeat audit, automatic proposal, and explicit adoption through
one command. It preserves the accepted sequence while changing repeat
representation or component connectivity. After rebuilding all components of an
organelle, use the complete organelle read pool:

```bash
orgraft refine --rebuild-dir 06.rebuild/mito --reads 01.recruit/mito.fastq.gz --sample S1 --organelle mito --out-dir 07.refine/mito --auto
```

The old `repeat-reclassify` name and
`polish --homorepair-*` options remain hidden compatibility interfaces.

Standalone repair of a candidate or manually edited GFA is available without
rerunning de novo graph construction:

```bash
orgraft asm --reads recruited.fastq.gz --organelle plastid \
  --skeleton-gfa graph.edited.gfa --stable --out-dir results/gfa_repair_r01
```

This mode still remaps reads to the supplied graph. Depth and link support are
conditional on that graph and read set, and the repaired output must pass the
ordinary checkpoint-1 review before `resolve`. Workflow TOML does not currently
generate this standalone repair command.

## Quick Start

```bash
cargo test
cargo build --release
export ORGRAFT_BIN="$PWD/target/release/orgraft"

"$ORGRAFT_BIN" --help
"$ORGRAFT_BIN" setup --soft-paths soft_paths.txt --requirements requirements.txt
"$ORGRAFT_BIN" workflow init --sample S1 --soft-paths "$PWD/soft_paths.txt" --out S1/orgraft.workflow.toml
# Edit the reads and references in S1/orgraft.workflow.toml.
"$ORGRAFT_BIN" workflow plan --config S1/orgraft.workflow.toml
bash S1/run_S1.sh
```

For direct source-tree execution:

```bash
cargo run -- workflow template
cargo run -- polish --help
```

## Workflow Summary

`orgraft workflow` is the orchestration layer. It reads
`orgraft.workflow.toml`, generates runnable scripts, and coordinates:

```text
recruit -> asm -> checkpoint1 -> resolve -> polish/checkpoint2 -> postrepair -> rebuild -> refine
```

The workflow layer stays thin: core algorithms remain in the individual
commands, while workflow owns project layout, command generation, manual
checkpoint status, one-at-a-time read-backed SV repair, and bounded SNV/InDel
correction rounds. SV check also projects abnormal breakpoints onto the
read-only `02.unitig_graph/graph.gfa` so manual review can distinguish an
internal unitig split from a link-level problem. Accepted SV repairs add a
validation round without consuming the ordinary round budget; total validation
rounds are capped at 10. A dominant repeat-pairing mismatch is recorded at
checkpoint 2 without changing resolve output. Rebuild then counts reads that
completely span each flank-repeat-flank path, uses the dominant pairing to
filter graph-valid circular candidates, and writes the constrained rebuild
FASTA. The next round is validate-only, so `02.polish` remains empty.
The workflow passes its completed validation round to `postrepair` as read-only
source evidence. New repair and validation files go to `05.postrepair`.
Evidence or review requirements remain explicit; an unresolved result stops the
handoff to rebuild. Refine receives all rebuilt components of one organelle
together, including the validation and source bindings for each component.

Each sample/assembly route, such as `samples/S1/software/curated`, has
`orgraft.workflow.toml`, `run_<sample>.sh`, and `run_<sample>_manual.sh`. The run
script contains the standard CLI commands directly. The manual script lists
them with the checkpoint checks and input assignments needed to run each stage.
Every subcommand is also available on its own. Relative paths in the TOML are
resolved from its directory. The seven output
roots are `01.recruit`, `02.draft_asm`, `03.resolve_gfa`, `04.polish`,
`05.postrepair`, `06.rebuild`, and `07.refine`. No AI service is needed at runtime.
Existing configs with explicit old output paths and the legacy homorepair option
remain usable; new configs declare postrepair and refine separately.

## Runtime Dependencies

`soft_paths.txt` records executable paths checked by `orgraft setup`:

- `python`
- `minimap2`
- `blastn`
- `pigz`
- optional `gfa_editor_cli` for reference-coloured PDF/SVG graph export

`requirements.txt` currently contains Python packages needed by OrgRAFT helper
scripts, especially polish validation plotting.

## Status

OrgRAFT is an active internal Rust CLI. Heavy alignment and visualization are
delegated to mature external tools where that remains the most reliable
boundary.
