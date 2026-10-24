use crate::error::OrgraftError;
use crate::{commands, workflow};

const HELP: &str = concat!(
    "Program: orgraft (Organelle Read-backed Assembly with FASTA–Graph Traceability)\n",
    "Version: ",
    env!("CARGO_PKG_VERSION"),
    "\n\n",
    "Usage:   orgraft [--help] [--version] <command> <argument>\n\n",
    "Commands:\n\n",
    " -- Project setup\n",
    "    setup       check external software paths and Python packages\n",
    "    workflow    generate and run end-to-end workflow checkpoints\n\n",
    " -- Raw graph generation\n",
    "    recruit     organelle HiFi read recruitment\n",
    "    asm         conservative draft graph assembly\n",
    "\n",
    " -- High-quality graph generation\n",
    "    resolve     resolve checked draft GFA into reference-oriented graph products\n",
    "    polish      polish linearized graph FASTA and evaluate variants\n",
    "    postrepair  review residual repeat, InDel, and circular-origin errors\n",
    "    rebuild     rebuild final verified graph and compact reports\n",
    "    refine      Refine repeat representation and component connectivity.\n\n",
    " Workflow runs are configured through orgraft.workflow.toml.\n"
);

pub fn run<I>(args: I) -> Result<(), OrgraftError>
where
    I: IntoIterator<Item = String>,
{
    let mut args = args.into_iter();
    let _binary = args.next();

    let Some(command) = args.next() else {
        print_help();
        return Ok(());
    };

    let rest: Vec<String> = args.collect();

    match command.as_str() {
        "-h" | "--help" | "help" => {
            print_help();
            Ok(())
        }
        "-V" | "--version" | "version" => {
            print_version();
            Ok(())
        }
        "workflow" => workflow::run(&rest),
        "recruit" => commands::recruit::run(&rest),
        "asm" => commands::asm::run(&rest),
        "resolve" => commands::resolve::run(&rest),
        "polish" => commands::polish::run(&rest),
        "postrepair" => commands::postrepair::run(&rest),
        "rebuild" => commands::rebuild::run(&rest),
        "refine" | "repeat-reclassify" => commands::refine::run(&rest),
        "setup" => commands::setup::run(&rest),
        _ => Err(OrgraftError::UnknownSubcommand(command)),
    }
}

fn print_help() {
    println!("{HELP}");
}

fn print_version() {
    println!("orgraft {}", env!("CARGO_PKG_VERSION"));
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn help_command_succeeds() {
        let args = vec!["orgraft".to_string(), "--help".to_string()];
        assert!(run(args).is_ok());
    }

    #[test]
    fn version_command_succeeds() {
        for flag in ["--version", "-V", "version"] {
            let args = vec!["orgraft".to_string(), flag.to_string()];
            assert!(run(args).is_ok());
        }
    }

    #[test]
    fn unknown_command_fails() {
        let args = vec!["orgraft".to_string(), "unknown".to_string()];
        assert!(matches!(run(args), Err(OrgraftError::UnknownSubcommand(_))));
    }

    #[test]
    fn high_quality_graph_commands_succeed() {
        for command in ["resolve", "polish", "postrepair", "rebuild", "refine", "workflow"] {
            let args = vec!["orgraft".to_string(), command.to_string()];
            assert!(run(args).is_ok());
        }
    }

    #[test]
    fn removed_scaffold_commands_fail() {
        for command in [
            "import", "graph", "graph-qc", "mx", "validate", "project", "curate", "report",
            "build", "config",
        ] {
            let args = vec!["orgraft".to_string(), command.to_string()];
            assert!(matches!(run(args), Err(OrgraftError::UnknownSubcommand(_))));
        }
    }
}
