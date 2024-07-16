use neospy_core::frames::*;
use pyo3::prelude::*;

/// Defined inertial frames supported by the python side of NEOSpy
#[pyclass(frozen, eq, eq_int, name = "Frames", module = "neospy")]
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum PyFrames {
    Ecliptic,
    Equatorial,
    Galactic,
    FK4,
    Undefined,
}

/// Provide a mapping from the python mapping above to the backend frames
impl From<PyFrames> for neospy_core::frames::Frame {
    fn from(value: PyFrames) -> Self {
        match value {
            PyFrames::Ecliptic => neospy_core::frames::Frame::Ecliptic,
            PyFrames::Equatorial => neospy_core::frames::Frame::Equatorial,
            PyFrames::Galactic => neospy_core::frames::Frame::Galactic,
            PyFrames::FK4 => neospy_core::frames::Frame::FK4,
            PyFrames::Undefined => neospy_core::frames::Frame::Unknown(0),
        }
    }
}

/// Provide a mapping from the python mapping above to the backend frames
impl From<neospy_core::frames::Frame> for PyFrames {
    fn from(value: neospy_core::frames::Frame) -> Self {
        match value {
            neospy_core::frames::Frame::Ecliptic => PyFrames::Ecliptic,
            neospy_core::frames::Frame::Equatorial => PyFrames::Equatorial,
            neospy_core::frames::Frame::FK4 => PyFrames::FK4,
            neospy_core::frames::Frame::Galactic => PyFrames::Galactic,
            _ => PyFrames::Undefined,
        }
    }
}

/// Convert a vector in the input frame to the same vector in the output frame.
#[allow(clippy::too_many_arguments)]
#[pyfunction]
#[pyo3(name = "frame_change")]
pub fn frame_change_py(
    states: Vec<[f64; 3]>,
    input_frame: PyFrames,
    output_frame: PyFrames,
) -> Vec<[f64; 3]> {
    states
        .into_iter()
        .map(|vec| {
            std::convert::Into::<Frame>::into(input_frame)
                .try_vec_frame_change(vec.into(), output_frame.into())
                .unwrap()
                .into()
        })
        .collect()
}

/// Compute a ECEF position from WCS84 Geodetic latitude/longitude/height.
///
/// This returns the X/Y/Z coordinates in km from the geocenter of the Earth.
///
/// Parameters
/// ----------
/// lat :
///     Latitude in degrees.
/// lon :
///     Longitude in degrees.
/// h :
///     Height above the surface of the Earth in km.
#[pyfunction]
#[pyo3(name = "wgs_lat_lon_to_ecef")]
pub fn wgs_lat_lon_to_ecef(lat: f64, lon: f64, h: f64) -> (f64, f64, f64) {
    geodetic_lat_lon_to_ecef(lat.to_radians(), lon.to_radians(), h)
}
