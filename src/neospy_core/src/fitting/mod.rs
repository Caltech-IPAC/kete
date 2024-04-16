//! # Fitting
//! Fitting tools, including root finding.
mod newton;
mod reduced_chi2;

pub use newton::newton_raphson;
pub use reduced_chi2::{fit_reduced_chi2, reduced_chi2, reduced_chi2_der};
