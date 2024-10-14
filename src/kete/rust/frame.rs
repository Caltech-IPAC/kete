//! Python Frame of reference support
use kete_core::frames::*;
use pyo3::prelude::*;

/// Defined inertial frames supported by the python side of kete
#[pyclass(frozen, eq, eq_int, name = "Frames", module = "kete")]
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum PyFrames {
    /// Ecliptic Frame
    Ecliptic,
    /// Equatorial Frame
    Equatorial,
    /// Galactic Frame
    Galactic,
    /// FK4 Frame
    FK4,
    /// Undefined Frame
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

/// Compute WCS84 Geodetic latitude/longitude/height from a ECEF position.
///
/// This returns the lat, lon, and height from the WGS84 oblate Earth.
///
/// Parameters
/// ----------
/// x :
///     ECEF x position in km.
/// y :
///     ECEF y position in km.
/// z :
///     ECEF z position in km.
#[pyfunction]
#[pyo3(name = "ecef_to_wgs_lat_lon")]
pub fn ecef_to_wgs_lat_lon(x: f64, y: f64, z: f64) -> (f64, f64, f64) {
    ecef_to_geodetic_lat_lon(x, y, z)
}

/// Calculate the obliquity angle of the Earth at the specified time.
///
/// This is only valid for several centuries near J2000.
///
/// The equation is from the 2010 Astronomical Almanac.
///
/// Parameters
/// ----------
/// time:
///     Calculate the obliquity angle of the Earth at the specified time.
#[pyfunction]
#[pyo3(name = "compute_obliquity")]
pub fn calc_obliquity_py(time: f64) -> f64 {
    calc_obliquity(time).to_degrees()
}
