use std::process::ExitCode;

fn main() -> ExitCode {
    match orgraft::run(std::env::args()) {
        Ok(()) => ExitCode::SUCCESS,
        Err(error) => {
            eprintln!("error: {error}");
            ExitCode::FAILURE
        }
    }
}
