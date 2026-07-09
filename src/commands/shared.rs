use std::collections::HashMap;
use std::path::{Path, PathBuf};
use std::process::Command;

pub struct CommandContract {
    pub command: &'static str,
    pub origin: &'static str,
    pub purpose: &'static str,
    pub inputs: &'static [&'static str],
    pub outputs: &'static [&'static str],
    pub notes: &'static [&'static str],
}

pub fn print_contract(contract: &CommandContract) {
    println!("orgraft {}", contract.command);
    println!();
    println!("Origin: {}", contract.origin);
    println!("Purpose: {}", contract.purpose);
    println!();
    print_list("Inputs", contract.inputs);
    print_list("Outputs", contract.outputs);
    print_list("Notes", contract.notes);
}

fn print_list(label: &str, values: &[&str]) {
    if values.is_empty() {
        return;
    }

    println!("{label}:");
    for value in values {
        println!("  - {value}");
    }
    println!();
}

#[derive(Debug, Clone)]
pub struct GfaImageExport {
    pub format: String,
    pub output: PathBuf,
    pub command: String,
    pub status: String,
    pub stdout: String,
    pub stderr: String,
}

pub fn resolve_gfa_editor_cli(soft_paths: &HashMap<String, PathBuf>) -> Result<PathBuf, String> {
    soft_paths
        .get("gfa_editor_cli")
        .cloned()
        .ok_or_else(|| "missing gfa_editor_cli in soft paths for optional image export".to_string())
}

pub fn run_gfa_editor_image(
    gfa_editor_cli: &Path,
    soft_paths: &HashMap<String, PathBuf>,
    input_gfa: &Path,
    output_path: &Path,
    image_reference_fasta: &Path,
    format: &str,
) -> GfaImageExport {
    let mut command = command_for_python_script(gfa_editor_cli, soft_paths);
    command.extend([
        "image".to_string(),
        input_gfa.display().to_string(),
        output_path.display().to_string(),
        "--colour".to_string(),
        "blastsolid".to_string(),
        "--query".to_string(),
        image_reference_fasta.display().to_string(),
        "--alignment-tool".to_string(),
        "minimap2".to_string(),
        "--layout".to_string(),
        "bandage".to_string(),
        "--target-role".to_string(),
        "subject".to_string(),
        "--alignment-args".to_string(),
        "-x asm5 -c --secondary=yes".to_string(),
    ]);
    let command_text = command.join(" ");
    let Some((program, args)) = command.split_first() else {
        return GfaImageExport {
            format: format.to_string(),
            output: output_path.to_path_buf(),
            command: command_text,
            status: "failed".to_string(),
            stdout: String::new(),
            stderr: "empty GFA_Editor command".to_string(),
        };
    };
    let output = match Command::new(program).args(args).output() {
        Ok(output) => output,
        Err(error) => {
            return GfaImageExport {
                format: format.to_string(),
                output: output_path.to_path_buf(),
                command: command_text,
                status: "failed_to_run".to_string(),
                stdout: String::new(),
                stderr: format!("failed to run GFA_Editor image export: {error}"),
            };
        }
    };
    let stdout = String::from_utf8_lossy(&output.stdout).to_string();
    let stderr = String::from_utf8_lossy(&output.stderr).to_string();
    let status = if output.status.success() {
        "written"
    } else {
        "failed"
    };
    GfaImageExport {
        format: format.to_string(),
        output: output_path.to_path_buf(),
        command: command_text,
        status: status.to_string(),
        stdout,
        stderr,
    }
}

fn command_for_python_script(path: &Path, soft_paths: &HashMap<String, PathBuf>) -> Vec<String> {
    let path_text = path.display().to_string();
    if path.extension().and_then(|value| value.to_str()) == Some("py") {
        if let Some(parent) = path.parent().and_then(Path::parent) {
            let local_python = parent.join(".venv/bin/python");
            if local_python.exists() {
                return vec![local_python.display().to_string(), path_text];
            }
        }
        if let Some(python) = soft_paths.get("python") {
            return vec![python.display().to_string(), path_text];
        }
        return vec!["python3".to_string(), path_text];
    }
    vec![path_text]
}
