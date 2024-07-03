//! # Errors
//! Errors emitted by NEOSpy_core

/// Define all errors which may be raise by this crate, as well as optionally provide
/// conversion to pyo3 error types which allow for the errors to be raised in Python.
use std::{error::Error, fmt, io};

/// Possible Errors which may be raised by this crate.
#[derive(Debug, Clone)]
pub enum NEOSpyError {
    /// Numerical method did not converge within the algorithms limits.
    Convergence(String),

    /// Input or variable exceeded expected or allowed bounds.
    ValueError(String),

    /// Querying an SPK file failed due to it missing the requisite data.
    DAFLimits(String),

    /// Attempting to load or convert to/from an Frame of reference which is not known.
    UnknownFrame(usize),

    /// Error related to IO.
    IOError(String),

    /// Propagator detected an impact.
    Impact(i64, f64),
}

impl Error for NEOSpyError {}

impl fmt::Display for NEOSpyError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            NEOSpyError::Convergence(s) => {
                write!(f, "{}", s)
            }
            NEOSpyError::ValueError(s) => {
                write!(f, "{}", s)
            }
            NEOSpyError::DAFLimits(s) => {
                write!(f, "{}", s)
            }
            NEOSpyError::UnknownFrame(_) => {
                write!(f, "This reference frame is not supported.")
            }
            NEOSpyError::IOError(s) => {
                write!(f, "{}", s)
            }
            NEOSpyError::Impact(s, t) => {
                write!(f, "Propagation detected an impact with {} at time {}", s, t)
            }
        }
    }
}

#[cfg(feature = "pyo3")]
use pyo3::{exceptions, PyErr};

#[cfg(feature = "pyo3")]
impl From<NEOSpyError> for PyErr {
    fn from(err: NEOSpyError) -> PyErr {
        match err {
            NEOSpyError::Convergence(s) => PyErr::new::<exceptions::PyValueError, _>(s),

            NEOSpyError::ValueError(s) => PyErr::new::<exceptions::PyValueError, _>(s),

            NEOSpyError::DAFLimits(s) => PyErr::new::<exceptions::PyValueError, _>(s),

            NEOSpyError::UnknownFrame(_) => {
                PyErr::new::<exceptions::PyValueError, _>("This reference frame is not supported.")
            }

            NEOSpyError::IOError(s) => PyErr::new::<exceptions::PyValueError, _>(s),

            NEOSpyError::Impact(s, t) => PyErr::new::<exceptions::PyValueError, _>(format!(
                "Propagation detected an impact with {} at time {}",
                s, t
            )),
        }
    }
}

impl From<io::Error> for NEOSpyError {
    fn from(error: io::Error) -> Self {
        NEOSpyError::IOError(error.to_string())
    }
}

impl From<std::num::ParseIntError> for NEOSpyError {
    fn from(value: std::num::ParseIntError) -> Self {
        NEOSpyError::IOError(value.to_string())
    }
}
impl From<std::num::ParseFloatError> for NEOSpyError {
    fn from(value: std::num::ParseFloatError) -> Self {
        NEOSpyError::IOError(value.to_string())
    }
}
