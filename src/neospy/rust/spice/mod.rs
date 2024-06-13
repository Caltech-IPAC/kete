mod daf;
mod spk;
mod pck;

pub use daf::*;
pub use spk::*;
pub use pck::*;

use pyo3::pyfunction;


/// Return a list of MPC observatory codes, along with the latitude, longitude (deg),
/// altitude (m above the WGS84 surface), and name.
#[pyfunction]
#[pyo3(name = "observatory_codes")]
pub fn obs_codes() -> Vec<(f64, f64, f64, String, String)> {
    let mut codes = Vec::new();
    for row in neospy_core::io::obs_codes::OBS_CODES.iter() {
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
