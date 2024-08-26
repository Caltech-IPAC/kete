use kete_core::frames::*;
use pyo3::prelude::*;

/// Defined inertial frames supported by the python side of kete
#[pyclass(frozen, eq, eq_int, name = "Frames", module = "kete")]
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum PyFrames {
    Ecliptic,
    Equatorial,
    Galactic,
    FK4,
    Undefined,
}

/// Provide a mapping from the python mapping above to the backend frames
impl From<PyFrames> for Frame {
    fn from(value: PyFrames) -> Self {
        match value {
            PyFrames::Ecliptic => Frame::Ecliptic,
            PyFrames::Equatorial => Frame::Equatorial,
            PyFrames::Galactic => Frame::Galactic,
            PyFrames::FK4 => Frame::FK4,
            PyFrames::Undefined => Frame::Unknown(0),
        }
    }
}

/// Provide a mapping from the python mapping above to the backend frames
impl From<Frame> for PyFrames {
    fn from(value: Frame) -> Self {
        match value {
            Frame::Ecliptic => PyFrames::Ecliptic,
            Frame::Equatorial => PyFrames::Equatorial,
            Frame::FK4 => PyFrames::FK4,
            Frame::Galactic => PyFrames::Galactic,
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
            Into::<Frame>::into(input_frame)
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
