use std::fs;
use std::io::{self, Write};
use std::path::{Path, PathBuf};
use std::process::Command;

use crate::error::OrgraftError;

const HELP: &str = r#"orgraft setup

Check external software paths and Python package requirements.

Usage:
  orgraft setup [--soft-paths soft_paths.txt] [--requirements requirements.txt]
"#;

const VERSION_ARGS: &[&[&str]] = &[&["--version"], &["-version"], &["-V"], &["version"]];
const PY_VERSION_CODE: &str = r#"import importlib.metadata as m, sys
try:
    print(m.version(sys.argv[1]))
except m.PackageNotFoundError:
    sys.exit(1)
"#;

struct Tool {
    name: String,
    path: PathBuf,
}

pub fn run(args: &[String]) -> Result<(), OrgraftError> {
    if args.iter().any(|arg| arg == "-h" || arg == "--help") {
        println!("{HELP}");
        return Ok(());
    }

    let soft_paths = option_value(args, "--soft-paths")?.unwrap_or("soft_paths.txt");
    let requirements = option_value(args, "--requirements")?.unwrap_or("requirements.txt");
    let tools = read_soft_paths(Path::new(soft_paths))?;

    println!("Software:");
    let mut passed = true;
    for tool in &tools {
        match check_tool(tool) {
            ToolStatus::Ok(version) => println!("  OK   {:18} {}", tool.name, version),
            ToolStatus::NoVersion => {
                println!("  OK   {:18} installed; version not detected", tool.name)
            }
            ToolStatus::Problem(message) => {
                passed = false;
                println!("  MISS {:18} {}", tool.name, message);
            }
        }
    }

    println!();
    passed &= check_python_packages(&tools, Path::new(requirements))?;

    println!();
    if passed {
        println!("Setup check passed.");
        Ok(())
    } else {
        println!("Setup check found problems.");
        io::stdout().flush()?;
        Err(OrgraftError::InvalidArgument(
            "one or more setup checks failed".to_string(),
        ))
    }
}

enum ToolStatus {
    Ok(String),
    NoVersion,
    Problem(String),
}

fn check_tool(tool: &Tool) -> ToolStatus {
    if !tool.path.is_absolute() {
        return ToolStatus::Problem(format!("path is not absolute: {}", tool.path.display()));
    }

    let metadata = match fs::metadata(&tool.path) {
        Ok(metadata) => metadata,
        Err(_) => return ToolStatus::Problem(format!("not found: {}", tool.path.display())),
    };

    if !metadata.is_file() {
        return ToolStatus::Problem(format!("not a file: {}", tool.path.display()));
    }

    if !is_executable(&metadata) {
        return ToolStatus::Problem(format!("not executable: {}", tool.path.display()));
    }

    detect_version(&tool.path)
        .map(ToolStatus::Ok)
        .unwrap_or(ToolStatus::NoVersion)
}

fn check_python_packages(tools: &[Tool], requirements: &Path) -> Result<bool, OrgraftError> {
    println!("Python packages:");

    let packages = match read_requirements(requirements)? {
        Some(packages) => packages,
        None => {
            println!(
                "  MISS requirements file not found: {}",
                requirements.display()
            );
            return Ok(false);
        }
    };

    let Some(python) = find_python(tools) else {
        println!("  MISS no python executable found in software paths");
        return Ok(false);
    };

    if !matches!(
        check_tool(python),
        ToolStatus::Ok(_) | ToolStatus::NoVersion
    ) {
        println!(
            "  MISS python executable is not usable: {}",
            python.path.display()
        );
        return Ok(false);
    }

    let mut passed = true;
    for package in packages {
        match package_version(&python.path, &package) {
            Some(version) => println!("  OK   {:18} {}", package, version),
            None => {
                passed = false;
                println!(
                    "  MISS {:18} not installed in {}",
                    package,
                    python.path.display()
                );
            }
        }
    }

    Ok(passed)
}

fn detect_version(path: &Path) -> Option<String> {
    for args in VERSION_ARGS {
        let output = Command::new(path).args(*args).output().ok()?;
        if output.status.success() {
            if let Some(line) = first_line(&output.stdout).or_else(|| first_line(&output.stderr)) {
                return Some(line);
            }
        }
    }
    None
}

fn package_version(python: &Path, package: &str) -> Option<String> {
    let output = Command::new(python)
        .args(["-c", PY_VERSION_CODE, package])
        .output()
        .ok()?;
    output
        .status
        .success()
        .then(|| first_line(&output.stdout))
        .flatten()
}

fn read_soft_paths(path: &Path) -> Result<Vec<Tool>, OrgraftError> {
    let text = fs::read_to_string(path).map_err(|error| {
        OrgraftError::InvalidArgument(format!("cannot read {}: {error}", path.display()))
    })?;

    let mut tools = Vec::new();
    for (index, line) in text.lines().enumerate() {
        let line = strip_comment(line).trim();
        if line.is_empty() {
            continue;
        }

        let (name, path) = split_tool_line(line).ok_or_else(|| {
            OrgraftError::InvalidArgument(format!(
                "{}:{} expected software_name<TAB>absolute_path_to_executable",
                path.display(),
                index + 1
            ))
        })?;

        tools.push(Tool {
            name: name.to_string(),
            path: PathBuf::from(path),
        });
    }

    if tools.is_empty() {
        return Err(OrgraftError::InvalidArgument(format!(
            "no software paths found in {}",
            path.display()
        )));
    }

    Ok(tools)
}

fn read_requirements(path: &Path) -> Result<Option<Vec<String>>, OrgraftError> {
    if !path.exists() {
        return Ok(None);
    }

    let text = fs::read_to_string(path)?;
    Ok(Some(
        text.lines()
            .filter_map(parse_requirement)
            .collect::<Vec<_>>(),
    ))
}

fn split_tool_line(line: &str) -> Option<(&str, &str)> {
    line.split_once('\t')
        .or_else(|| line.split_once(char::is_whitespace))
        .map(|(name, path)| (name.trim(), path.trim()))
        .filter(|(name, path)| !name.is_empty() && !path.is_empty())
}

fn parse_requirement(line: &str) -> Option<String> {
    let requirement = strip_comment(line).trim();
    if requirement.is_empty() {
        return None;
    }

    let name = requirement
        .split(['=', '<', '>', '~', '!', '[', ';'])
        .next()
        .unwrap_or("")
        .trim();

    (!name.is_empty()).then(|| name.to_string())
}

fn strip_comment(line: &str) -> &str {
    line.split_once('#').map(|(value, _)| value).unwrap_or(line)
}

fn find_python(tools: &[Tool]) -> Option<&Tool> {
    tools
        .iter()
        .find(|tool| tool.name.eq_ignore_ascii_case("python"))
        .or_else(|| {
            tools.iter().find(|tool| {
                tool.path
                    .file_name()
                    .and_then(|name| name.to_str())
                    .is_some_and(|name| name.starts_with("python"))
            })
        })
}

fn first_line(output: &[u8]) -> Option<String> {
    String::from_utf8_lossy(output)
        .lines()
        .map(str::trim)
        .find(|line| !line.is_empty())
        .map(str::to_string)
}

#[cfg(unix)]
fn is_executable(metadata: &fs::Metadata) -> bool {
    use std::os::unix::fs::PermissionsExt;
    metadata.permissions().mode() & 0o111 != 0
}

#[cfg(not(unix))]
fn is_executable(_: &fs::Metadata) -> bool {
    true
}

fn option_value<'a>(args: &'a [String], name: &str) -> Result<Option<&'a str>, OrgraftError> {
    for (index, arg) in args.iter().enumerate() {
        if arg == name {
            return args
                .get(index + 1)
                .map(|value| Some(value.as_str()))
                .ok_or_else(|| OrgraftError::InvalidArgument(format!("missing value for {name}")));
        }
    }
    Ok(None)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_tab_separated_tool_line() {
        let (name, path) = split_tool_line("python\t/usr/bin/python3").unwrap();
        assert_eq!(name, "python");
        assert_eq!(path, "/usr/bin/python3");
    }

    #[test]
    fn parses_requirement_name_before_specifier() {
        assert_eq!(
            parse_requirement("numpy>=2.0 # numerical arrays"),
            Some("numpy".to_string())
        );
    }
}
