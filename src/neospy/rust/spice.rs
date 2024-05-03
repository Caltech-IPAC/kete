use neospy_core::prelude::*;
use neospy_core::spice::{get_pck_singleton, get_spk_singleton, try_name_from_id, Daf};
use pyo3::{pyfunction, PyResult};

use crate::frame::PyFrames;
use crate::state::PyState;
use crate::vector::Vector;

#[pyfunction]
#[pyo3(name = "spk_load")]
pub fn spk_load_py(filenames: Vec<String>) -> PyResult<()> {
    let mut singleton = get_spk_singleton().write().unwrap();
    if filenames.len() > 100 {
        eprintln!("Loading {} spk files...", filenames.len());
    }
    for filename in filenames.iter() {
        let load = (*singleton).load_file(filename);
        if let Err(err) = load {
            eprintln!("{} failed to load. {}", filename, err);
        }
    }
    (*singleton).build_cache();
    Ok(())
}

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
/// instantaneous earth frame information. The one provided by neospy has dates from
/// around 2000 to early 2024. New files may be downloaded from the NAIF website for
/// additional precision or years of epoch. The current file is accurate to ~5 cm
/// precision.
#[pyfunction]
#[pyo3(name = "pck_earth_frame_to_ecliptic")]
pub fn pck_earth_frame_py(
    pos: [f64; 3],
    jd: f64,
    new_center: isize,
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
    state.try_change_frame(Frame::Ecliptic)?;

    let spks = get_spk_singleton().try_read().unwrap();
    spks.try_change_center(&mut state, new_center)?;
    Ok(PyState(state))
}

/// Convert a [`State`] to the Earth's surface reference frame.
///
/// This requires the `earth_000101_*.pck` file to be loaded which contains the
/// instantaneous earth frame information. The one provided by neospy has dates from
/// around 2000 to early 2024. New files may be downloaded from the NAIF website for
/// additional precision or years of epoch. The current file is accurate to ~5 cm
/// precision.
#[pyfunction]
#[pyo3(name = "pck_state_to_frame")]
pub fn pck_state_to_earth(state: PyState) -> PyResult<Vector> {
    let pcks = get_pck_singleton().try_read().unwrap();
    let state = state.change_center(399)?.as_ecliptic()?;
    let frame = pcks.try_get_orientation(3000, state.jd())?;
    let mut state = state.0;

    state.try_change_frame(frame)?;

    Ok(Vector::new(state.pos, PyFrames::Undefined))
}

#[pyfunction]
#[pyo3(name = "spk_available_info")]
pub fn spk_available_info_py(naif_id: isize) -> Vec<(f64, f64, isize, String, usize)> {
    let singleton = get_spk_singleton().try_read().unwrap();
    let info = singleton.available_info(naif_id);

    info.iter()
        .map(|(t0, t1, center, frame, segment_type)| {
            (*t0, *t1, *center, frame.to_string(), *segment_type)
        })
        .collect()
}

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

#[pyfunction]
#[pyo3(name = "spk_loaded")]
pub fn spk_loaded_objects_py() -> Vec<isize> {
    let spk = get_spk_singleton().try_read().unwrap();
    let loaded = spk.loaded_objects(false);
    let mut loaded: Vec<isize> = loaded.into_iter().collect();
    loaded.sort();
    loaded
}

#[pyfunction]
#[pyo3(name = "spk_get_name_from_id")]
pub fn spk_get_name_from_id_py(id: isize) -> String {
    try_name_from_id(id).unwrap_or(id.to_string())
}

#[pyfunction]
#[pyo3(name = "spk_reset")]
pub fn spk_reset_py() {
    get_spk_singleton().write().unwrap().reset()
}

#[pyfunction]
#[pyo3(name = "pck_reset")]
pub fn pck_reset_py() {
    get_pck_singleton().write().unwrap().reset()
}

#[pyfunction]
#[pyo3(name = "spk_state")]
pub fn spk_state_py(id: isize, jd: f64, center: isize, frame: PyFrames) -> PyResult<PyState> {
    let spk = get_spk_singleton().try_read().unwrap();
    let mut state = spk.try_get_state(id, jd, center, frame.into())?;
    state.try_naif_id_to_name();
    Ok(PyState(state))
}

#[pyfunction]
#[pyo3(name = "spk_raw_state")]
pub fn spk_raw_state_py(id: isize, jd: f64) -> PyResult<PyState> {
    let spk = get_spk_singleton().try_read().unwrap();
    Ok(PyState(spk.try_get_raw_state(id, jd)?))
}

#[pyfunction]
#[pyo3(name = "daf_header_info")]
pub fn daf_header_info_py(filename: &str) -> PyResult<String> {
    let mut file = std::fs::File::open(filename)?;
    let daf = Daf::try_load_header(&mut file)?;
    Ok(daf.comments.join(""))
}
