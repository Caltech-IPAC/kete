mod daf;
mod pck;
mod spk;
mod tle;

pub use daf::*;
pub use pck::*;
pub use spk::*;
pub use tle::*;

use pyo3::pyfunction;

/// Return a list of MPC observatory codes, along with the latitude, longitude (deg),
/// altitude (m above the WGS84 surface), and name.
#[pyfunction]
#[pyo3(name = "observatory_codes")]
pub fn obs_codes() -> Vec<(f64, f64, f64, String, String)> {
    let mut codes = Vec::new();
    for row in kete_core::io::obs_codes::OBS_CODES.iter() {
        codes.push((
            row.lat,
            row.lon,
            row.altitude,
            row.name.clone(),
            row.code.clone(),
        ))
    }
    codes
}
