use kete_core::frames::ecef_to_geodetic_lat_lon;
use kete_core::spice::{get_pck_singleton, get_spk_singleton};
use kete_core::{constants, prelude::*};
use pyo3::{pyfunction, PyResult};

use crate::state::PyState;

/// Load all specified files into the PCK shared memory singleton.
#[pyfunction]
#[pyo3(name = "pck_load")]
pub fn pck_load_py(filenames: Vec<String>) -> PyResult<()> {
    let mut singleton = get_pck_singleton().write().unwrap();
    for filename in filenames.iter() {
        let load = (*singleton).load_file(filename);
        if let Err(err) = load {
            eprintln!("{} failed to load. {}", filename, err);
        }
    }
    Ok(())
}

/// Convert a position vector which is geocentered in the earth frame to a position
/// vector in the sun centered ecliptic frame.
///
/// This requires the `earth_000101_*.pck` file to be loaded which contains the
/// instantaneous earth frame information. The one provided by kete has dates from
/// around 2000 to early 2024. New files may be downloaded from the NAIF website for
/// additional precision or years of epoch. The current file is accurate to ~5 cm
/// precision.
///
/// parameters
/// ----------
/// pos : list[float]
///     Position of the object in the geocentric reference frame in AU.
/// jd : float
///     Time in JD using TDB scaling.
/// new_center: int
///     NAIF ID of the new center point.
/// name : String
///     Optional name of the state.
#[pyfunction]
#[pyo3(name = "pck_earth_frame_to_ecliptic", signature = (pos, jd, new_center, name=None))]
pub fn pck_earth_frame_py(
    pos: [f64; 3],
    jd: f64,
    new_center: i64,
    name: Option<String>,
) -> PyResult<PyState> {
    let desig = {
        match name {
            Some(d) => Desig::Name(d),
            None => Desig::Empty,
        }
    };
    let pcks = get_pck_singleton().try_read().unwrap();
    let frame = pcks.try_get_orientation(3000, jd)?;

    let mut state = State::new(desig, jd, pos.into(), [0.0, 0.0, 0.0].into(), frame, 399);
    state.try_change_frame_mut(Frame::Ecliptic)?;

    let spks = get_spk_singleton().try_read().unwrap();
    spks.try_change_center(&mut state, new_center)?;
    Ok(PyState(state))
}

/// Convert a [`State`] to the Earth's surface lat/lon/height on the WGS84 reference.
///
/// This requires the `earth_000101_*.pck` file to be loaded which contains the
/// instantaneous earth frame information. The one provided by kete has dates from
/// around 2000 to early 2024, along with predicts into the future. New files may be
/// downloaded from the NAIF website for additional precision or years of epoch.
///
/// Parameters
/// ----------
/// state: State
///     Convert the given state to latitude, longitude, and height in km on the WGS84
///     reference.
#[pyfunction]
#[pyo3(name = "state_to_earth_pos")]
pub fn pck_state_to_earth(state: PyState) -> PyResult<(f64, f64, f64)> {
    let pcks = get_pck_singleton().try_read().unwrap();
    let state = state.change_center(399)?.as_ecliptic()?;
    let frame = pcks.try_get_orientation(3000, state.jd())?;
    let mut state = state.0;

    state.try_change_frame_mut(frame)?;
    let [x, y, z] = state.pos;
    let (lat, lon, height) = ecef_to_geodetic_lat_lon(
        x * constants::AU_KM,
        y * constants::AU_KM,
        z * constants::AU_KM,
    );

    Ok((lat.to_degrees(), lon.to_degrees(), height))
}

/// Reset the contents of the PCK shared memory to the default set of PCK kernels.
#[pyfunction]
#[pyo3(name = "pck_reset")]
pub fn pck_reset_py() {
    get_pck_singleton().write().unwrap().reset()
}
