//! Python vector support with frame information.
use kete_core::frames::rotate_around;
use pyo3::basic::CompareOp;
use pyo3::exceptions;
use pyo3::exceptions::PyNotImplementedError;
use std::f64::consts::FRAC_PI_2;

use crate::frame::*;
use kete_core::prelude::*;
use nalgebra::Vector3;
use pyo3::prelude::*;
use pyo3::PyResult;

/// Vector class which is a vector along with a reference frame.
///
/// Parameters
/// ----------
/// raw : list
///     3 floats which define the direction of the vector.
/// frame :
///     The frame of reference defining the coordinate frame of the vector, defaults
///     to ecliptic.
#[pyclass(sequence, frozen, module = "kete")]
#[derive(Clone, Debug)]
pub struct Vector {
    /// X/Y/Z numbers of the vector
    pub raw: [f64; 3],

    frame: PyFrames,
}

impl Vector {
    /// Construct a new vector
    pub fn new(raw: [f64; 3], frame: PyFrames) -> Self {
        Self { raw, frame }
    }
}

/// Polymorphic support
#[derive(Debug, FromPyObject, IntoPyObject)]
pub enum VectorLike {
    /// Vector directly
    Vec(Vector),

    /// Vector from x/y/z
    Arr([f64; 3]),
}

impl VectorLike {
    /// Cast VectorLike into a Vector3
    pub fn into_vec(self, target_frame: PyFrames) -> Vector3<f64> {
        match self {
            VectorLike::Arr(arr) => Vector3::from(arr),
            VectorLike::Vec(mut vec) => {
                if vec.frame() != target_frame {
                    vec = vec.change_frame(target_frame);
                };
                Vector3::from(vec.raw)
            }
        }
    }

    /// Cast VectorLike into a python Vector
    pub fn into_vector(self, target_frame: PyFrames) -> Vector {
        let vec = self.into_vec(target_frame);
        Vector {
            raw: vec.into(),
            frame: target_frame,
        }
    }
}

#[pymethods]
impl Vector {
    /// create new vector
    #[new]
    #[pyo3(signature = (raw, frame=None))]
    pub fn py_new(raw: VectorLike, frame: Option<PyFrames>) -> PyResult<Self> {
        match raw {
            VectorLike::Arr(raw) => {
                let frame = frame.unwrap_or(PyFrames::Ecliptic);
                Ok(Self { raw, frame })
            }
            VectorLike::Vec(vec) => {
                if frame.is_some() {
                    return Err(Error::ValueError(
                        "If a vector is provided, then the frame cannot be specified.".into(),
                    )
                    .into());
                }
                Ok(vec)
            }
        }
    }

    /// Create a new Vector from the elevation and azimuthal angle in degrees.
    ///
    /// Parameters
    /// ----------
    /// el : float
    ///     Elevation above the X-Y plane of the frame. (Degrees)
    /// az : float
    ///     Azimuthal angle on the X-Y plane for the frame. (Degrees)
    /// frame : Frames
    ///     Frame of reference which define the coordinate axis.
    /// r :
    ///     Optional length of the vector, defaults to 1.
    #[staticmethod]
    pub fn from_el_az(el: f64, az: f64, r: f64, frame: PyFrames) -> Self {
        let (el_sin, el_cos) = (FRAC_PI_2 - el.to_radians()).sin_cos();
        let (az_sin, az_cos) = az.to_radians().sin_cos();
        let x = r * el_sin * az_cos;
        let y = r * el_sin * az_sin;
        let z = r * el_cos;
        Self::new([x, y, z], frame)
    }

    /// Create a new Ecliptic Vector with the specified latitude/longitude.
    ///
    /// Parameters
    /// ----------
    /// lat : float
    ///     Latitude in the ecliptic frame. (Degrees)
    /// lon : float
    ///     Longitude in the ecliptic frame. (Degrees)
    /// r :
    ///     Optional length of the vector, defaults to 1.
    #[staticmethod]
    #[pyo3(signature = (lat, lon, r=None))]
    pub fn from_lat_lon(lat: f64, lon: f64, r: Option<f64>) -> Self {
        Self::from_el_az(lat, lon, r.unwrap_or(1.0), PyFrames::Ecliptic)
    }

    /// Create a new Equatorial Vector with the specified RA/DEC.
    ///
    /// Parameters
    /// ----------
    /// ra : float
    ///     Right Ascension in the equatorial frame. (Degrees)
    /// dec : float
    ///     Declination in the equatorial frame. (Degrees)
    /// r :
    ///     Optional length of the vector, defaults to 1.
    #[staticmethod]
    #[pyo3(signature = (ra, dec, r=None))]
    pub fn from_ra_dec(ra: f64, dec: f64, r: Option<f64>) -> Self {
        Self::from_el_az(dec, ra, r.unwrap_or(1.0), PyFrames::Equatorial)
    }

    /// The raw vector without the Frame.
    #[getter]
    fn raw(&self) -> [f64; 3] {
        self.raw
    }

    /// The Frame of reference.
    #[getter]
    pub fn frame(&self) -> PyFrames {
        self.frame
    }

    /// Length of the Vector
    #[getter]
    pub fn r(&self) -> f64 {
        let data: &Vector3<f64> = &self.raw.into();
        data.norm()
    }

    /// Azimuth in degrees from the X axis in the X-Y plane of the coordinate frame.
    #[getter]
    pub fn az(&self) -> f64 {
        let data: &Vector3<f64> = &self.raw.into();
        let r = data.norm();
        if r < 1e-8 {
            return 0.0;
        }
        f64::atan2(data.y, data.x).to_degrees().rem_euclid(360.0)
    }

    /// Elevation in degrees from the X-Y plane of the coordinate frame.
    /// Values will be between -180 and 180
    #[getter]
    pub fn el(&self) -> f64 {
        let data: &Vector3<f64> = &self.raw.into();
        let r = data.norm();
        if r < 1e-8 {
            return 0.0;
        }
        ((FRAC_PI_2 - (data.z / r).clamp(-1.0, 1.0).acos()).to_degrees() + 180.0).rem_euclid(360.0)
            - 180.0
    }

    /// Right Ascension in degrees if the frame is Equatorial.
    #[getter]
    pub fn ra(&self) -> PyResult<f64> {
        if self.frame != PyFrames::Equatorial {
            return Err(Error::ValueError(
                "Cannot compute RA as the frame is not equatorial. Change frame to equatorial before calling ra/dec."
                    .into(),
            )
            .into());
        }
        Ok(self.az())
    }

    /// Declination in degrees if the frame is Equatorial.
    #[getter]
    pub fn dec(&self) -> PyResult<f64> {
        if self.frame != PyFrames::Equatorial {
            return Err(Error::ValueError(
                "Cannot compute Dec as the frame is not equatorial. Change frame to equatorial before calling ra/dec."
                    .into(),
            )
            .into());
        }
        Ok(self.el())
    }

    /// Latitude in degrees if the frame is Ecliptic.
    #[getter]
    pub fn lat(&self) -> PyResult<f64> {
        if self.frame != PyFrames::Ecliptic {
            return Err(Error::ValueError(
                "Cannot compute Latitude as the frame is not ecliptic. Change frame to ecliptic."
                    .into(),
            )
            .into());
        }
        Ok(self.el())
    }

    /// Longitude in degrees if the frame is Ecliptic.
    #[getter]
    pub fn lon(&self) -> PyResult<f64> {
        if self.frame != PyFrames::Ecliptic {
            return Err(Error::ValueError(
                "Cannot compute Longitude as the frame is not ecliptic. Change frame to ecliptic."
                    .into(),
            )
            .into());
        }
        Ok(self.az())
    }

    /// Compute the angle in degrees between two vectors in degrees.
    /// This will automatically make a frame change if necessary.
    pub fn angle_between(&self, other: VectorLike) -> f64 {
        let self_vec = Vector3::from(self.raw);
        let other_vec = other.into_vec(self.frame());
        self_vec.angle(&other_vec).to_degrees()
    }

    /// Return the vector in the ecliptic frame, regardless of starting frame.
    #[getter]
    pub fn as_ecliptic(&self) -> Self {
        self.change_frame(PyFrames::Ecliptic)
    }

    /// Return the vector in the equatorial frame, regardless of starting frame.
    #[getter]
    pub fn as_equatorial(&self) -> Self {
        self.change_frame(PyFrames::Equatorial)
    }

    /// Return the vector in the galactic frame, regardless of starting frame.
    #[getter]
    pub fn as_galactic(&self) -> Self {
        self.change_frame(PyFrames::Galactic)
    }

    /// Return the vector in the fk4 frame, regardless of starting frame.
    #[getter]
    pub fn as_fk4(&self) -> Self {
        self.change_frame(PyFrames::FK4)
    }

    /// Return the vector in the target frame, regardless of starting frame.
    pub fn change_frame(&self, target_frame: PyFrames) -> Self {
        let new_dat = Into::<Frame>::into(self.frame)
            .try_vec_frame_change(self.raw.into(), target_frame.into())
            .unwrap();
        Self::new(new_dat.into(), target_frame)
    }

    /// X coordinate in au.
    #[getter]
    pub fn x(&self) -> f64 {
        self.raw[0]
    }

    /// Y coordinate in au.
    #[getter]
    pub fn y(&self) -> f64 {
        self.raw[1]
    }

    /// Z coordinate in au.
    #[getter]
    pub fn z(&self) -> f64 {
        self.raw[2]
    }

    /// Rotate this vector around another vector by the provided angle.
    ///
    ///
    /// Parameters
    /// ----------
    /// other : Vector
    ///     The other vector to rotate around.
    /// angle :
    ///     The angle in degrees of the rotation.
    pub fn rotate_around(&self, other: VectorLike, angle: f64) -> Self {
        let self_vec = Vector3::from(self.raw);
        let other_vec = other.into_vec(self.frame());
        let rotated = rotate_around(&self_vec, other_vec, angle.to_radians());
        Self::new(rotated.into(), self.frame)
    }

    #[allow(missing_docs)]
    pub fn __repr__(&self) -> String {
        // 1e-12 AU is about 15cm, this seems like a reasonable printing resolution
        let x = (self.raw[0] * 1e12).round() / 1e12 + 0.0;
        let y = (self.raw[1] * 1e12).round() / 1e12 + 0.0;
        let z = (self.raw[2] * 1e12).round() / 1e12 + 0.0;
        format!("Vector([{:?}, {:?}, {:?}], {:?})", x, y, z, self.frame)
    }

    #[allow(missing_docs)]
    pub fn __sub__(&self, other: VectorLike) -> Self {
        let self_vec = Vector3::from(self.raw);
        let other_vec = other.into_vec(self.frame());
        let diff = self_vec - other_vec;
        Self::new(diff.into(), self.frame)
    }

    #[allow(missing_docs)]
    pub fn __add__(&self, other: VectorLike) -> Self {
        let self_vec = Vector3::from(self.raw);
        let other_vec = other.into_vec(self.frame());
        let diff = self_vec + other_vec;
        Self::new(diff.into(), self.frame)
    }

    #[allow(missing_docs)]
    pub fn __mul__(&self, other: f64) -> Self {
        let self_vec = Vector3::from(self.raw);
        Self::new((self_vec * other).into(), self.frame)
    }

    #[allow(missing_docs)]
    pub fn __truediv__(&self, other: f64) -> Self {
        let self_vec = Vector3::from(self.raw);
        Self::new((self_vec / other).into(), self.frame)
    }

    #[allow(missing_docs)]
    pub fn __neg__(&self) -> Self {
        Self::new([-self.x(), -self.y(), -self.z()], self.frame)
    }

    #[allow(missing_docs)]
    pub fn __len__(&self) -> usize {
        3
    }

    #[allow(missing_docs)]
    pub fn __getitem__(&self, idx: usize) -> PyResult<f64> {
        if idx >= 3 {
            return Err(PyErr::new::<exceptions::PyIndexError, _>(""));
        }
        Ok(self.raw[idx])
    }

    fn __richcmp__(&self, other: VectorLike, op: CompareOp, _py: Python<'_>) -> PyResult<bool> {
        let self_vec = Vector3::from(self.raw);
        let other_vec = other.into_vec(self.frame());
        match op {
            CompareOp::Eq => Ok((self_vec - other_vec).norm() < 1e-12),
            CompareOp::Ne => Ok((self_vec - other_vec).norm() >= 1e-12),
            _ => Err(PyNotImplementedError::new_err(
                "Vectors can only be checked for equality.",
            )),
        }
    }
}
