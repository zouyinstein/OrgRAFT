use std::fmt;

#[derive(Debug)]
pub enum OrgraftError {
    InvalidArgument(String),
    UnknownSubcommand(String),
    Io(std::io::Error),
}

impl fmt::Display for OrgraftError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidArgument(message) => write!(formatter, "{message}"),
            Self::UnknownSubcommand(command) => write!(formatter, "unknown subcommand `{command}`"),
            Self::Io(error) => write!(formatter, "{error}"),
        }
    }
}

impl std::error::Error for OrgraftError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Io(error) => Some(error),
            _ => None,
        }
    }
}

impl From<std::io::Error> for OrgraftError {
    fn from(error: std::io::Error) -> Self {
        Self::Io(error)
    }
}
