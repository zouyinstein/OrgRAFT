# OrgRAFT

OrgRAFT means **Organelle Graph Read-backed Assembly and FASTA Traceability**.

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
orgraft asm         # conservative draft graph assembly

-- High-quality graph generation
orgraft resolve     # resolve checked draft GFA into reference-oriented products
orgraft polish      # polish linearized graph FASTA and evaluate variants
orgraft rebuild     # rebuild final verified graph and compact reports
```

Run any command without arguments, or with `--help`, to print its current
interface.

## Quick Start

```bash
cargo test
cargo build

./orgraft --help
./orgraft setup --soft-paths soft_paths.txt --requirements requirements.txt
./orgraft workflow init --out results_workflow/orgraft.workflow.toml
./orgraft workflow plan --config results_workflow/orgraft.workflow.toml
bash results_workflow/workflow.commands.sh
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
recruit -> asm -> checkpoint1 -> resolve -> polish/checkpoint2 -> rebuild
```

The workflow layer stays thin: core algorithms remain in the individual
commands, while workflow owns project layout, command generation, manual
checkpoint status, and bounded SNV/InDel correction rounds.

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
