use neospy_core::frames::{Ecliptic, Equatorial, Galactic, Vector, FK4};
use pyo3::exceptions;
use std::f64::consts::FRAC_PI_2;

use crate::frame::*;
use nalgebra::Vector3;
use neospy_core::prelude::*;
use pyo3::class::basic::CompareOp;
use pyo3::prelude::*;
use pyo3::{PyResult, Python};

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum FrameVector {
    Ecliptic(Vector<Ecliptic>),
    Equatorial(Vector<Equatorial>),
    Galactic(Vector<Galactic>),
    FK4(Vector<FK4>),
}

impl FrameVector {
    pub fn new(raw: [f64; 3], frame: PyFrames) -> Self {
        match frame {
            PyFrames::Ecliptic => FrameVector::Ecliptic(raw.into()),
            PyFrames::Equatorial => FrameVector::Equatorial(raw.into()),
            PyFrames::FK4 => FrameVector::FK4(raw.into()),
            PyFrames::Galactic => FrameVector::Galactic(raw.into()),
        }
    }

    pub fn frame(&self) -> PyFrames {
        match self {
            FrameVector::Ecliptic(_) => PyFrames::Ecliptic,
            FrameVector::Equatorial(_) => PyFrames::Equatorial,
            FrameVector::FK4(_) => PyFrames::FK4,
            FrameVector::Galactic(_) => PyFrames::Galactic,
        }
    }

    pub fn raw(self) -> [f64; 3] {
        match self {
            FrameVector::Ecliptic(v) => v.into(),
            FrameVector::Equatorial(v) => v.into(),
            FrameVector::FK4(v) => v.into(),
            FrameVector::Galactic(v) => v.into(),
        }
    }

    pub fn r(self) -> f64 {
        match self {
            FrameVector::Ecliptic(v) => v.norm(),
            FrameVector::Equatorial(v) => v.norm(),
            FrameVector::FK4(v) => v.norm(),
            FrameVector::Galactic(v) => v.norm(),
        }
    }

    pub fn into_frame(self, frame: PyFrames) -> Self {
        match frame {
            PyFrames::Ecliptic => self.to_ecliptic(),
            PyFrames::Equatorial => self.to_equatorial(),
            PyFrames::FK4 => self.to_fk4(),
            PyFrames::Galactic => self.to_galactic(),
        }
    }

    pub fn to_equatorial(self) -> Self {
        FrameVector::Equatorial(self.to_equatorial_raw())
    }

    pub fn to_ecliptic(self) -> Self {
        FrameVector::Ecliptic(self.to_ecliptic_raw())
    }

    pub fn to_fk4(self) -> Self {
        FrameVector::FK4(self.to_fk4_raw())
    }

    pub fn to_galactic(self) -> Self {
        FrameVector::Galactic(self.to_galactic_raw())
    }

    pub fn to_equatorial_raw(self) -> Vector<Equatorial> {
        match self {
            FrameVector::Ecliptic(v) => v.into_frame(),
            FrameVector::Equatorial(v) => v,
            FrameVector::FK4(v) => v.into_frame(),
            FrameVector::Galactic(v) => v.into_frame(),
        }
    }

    pub fn to_ecliptic_raw(self) -> Vector<Ecliptic> {
        match self {
            FrameVector::Ecliptic(v) => v,
            FrameVector::Equatorial(v) => v.into_frame(),
            FrameVector::FK4(v) => v.into_frame(),
            FrameVector::Galactic(v) => v.into_frame(),
        }
    }

    pub fn to_fk4_raw(self) -> Vector<FK4> {
        match self {
            FrameVector::Ecliptic(v) => v.into_frame(),
            FrameVector::Equatorial(v) => v.into_frame(),
            FrameVector::FK4(v) => v,
            FrameVector::Galactic(v) => v.into_frame(),
        }
    }

    pub fn to_galactic_raw(self) -> Vector<Galactic> {
        match self {
            FrameVector::Ecliptic(v) => v.into_frame(),
            FrameVector::Equatorial(v) => v.into_frame(),
            FrameVector::FK4(v) => v.into_frame(),
            FrameVector::Galactic(v) => v,
        }
    }
}

impl From<Vector<Equatorial>> for FrameVector {
    fn from(value: Vector<Equatorial>) -> Self {
        FrameVector::new(value.into(), PyFrames::Equatorial)
    }
}

impl From<Vector<Ecliptic>> for FrameVector {
    fn from(value: Vector<Ecliptic>) -> Self {
        FrameVector::new(value.into(), PyFrames::Ecliptic)
    }
}

impl From<Vector<Galactic>> for FrameVector {
    fn from(value: Vector<Galactic>) -> Self {
        FrameVector::new(value.into(), PyFrames::Galactic)
    }
}

impl From<Vector<FK4>> for FrameVector {
    fn from(value: Vector<FK4>) -> Self {
        FrameVector::new(value.into(), PyFrames::FK4)
    }
}

/// Vector class which is a vector along with a reference frame.
#[pyclass(sequence, frozen, module = "neospy")]
#[derive(Clone, Debug, Copy)]
pub struct PyVector(FrameVector);

impl PyVector {
    pub fn new(raw: [f64; 3], frame: PyFrames) -> Self {
        PyVector(FrameVector::new(raw, frame))
    }
}

/// Polymorphic support
#[derive(FromPyObject)]
pub enum VectorLike {
    PyVec(PyVector),
    Arr([f64; 3]),
}

impl VectorLike {
    pub fn into_pyvector(self, target_frame: PyFrames) -> PyVector {
        match self {
            VectorLike::Arr(raw) => PyVector(FrameVector::new(raw, target_frame)),
            VectorLike::PyVec(pyvec) => pyvec.change_frame(target_frame),
        }
    }
    pub fn into_ecliptic_raw(self) -> Vector<Ecliptic> {
        match self {
            VectorLike::Arr(raw) => Vector::<Ecliptic>::new(raw),
            VectorLike::PyVec(pyvec) => pyvec.0.to_ecliptic_raw(),
        }
    }
}

#[pymethods]
impl PyVector {
    #[new]
    pub fn py_new(raw: VectorLike, frame: Option<PyFrames>) -> Self {
        let frame = frame.unwrap_or(PyFrames::Ecliptic);
        raw.into_pyvector(frame)
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
    pub fn from_ra_dec(ra: f64, dec: f64, r: Option<f64>) -> Self {
        Self::from_el_az(dec, ra, r.unwrap_or(1.0), PyFrames::Equatorial)
    }

    /// The raw vector without the Frame.
    #[getter]
    fn raw(&self) -> [f64; 3] {
        self.0.raw()
    }

    /// The Frame of reference.
    #[getter]
    pub fn frame(&self) -> PyFrames {
        self.0.frame()
    }

    /// Length of the Vector
    #[getter]
    pub fn r(&self) -> f64 {
        self.0.r()
    }

    /// Azimuth in degrees from the X axis in the X-Y plane of the coordinate frame.
    #[getter]
    pub fn az(&self) -> f64 {
        let data: &Vector3<f64> = &self.raw().into();
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
        let data: &Vector3<f64> = &self.raw().into();
        let r = data.norm();
        if r < 1e-8 {
            return 0.0;
        }
        ((FRAC_PI_2 - ((data.z / r).clamp(-1.0, 1.0)).acos()).to_degrees() + 180.0)
            .rem_euclid(360.0)
            - 180.0
    }

    /// Right Ascension in degrees if the frame is Equatorial.
    #[getter]
    pub fn ra(&self) -> PyResult<f64> {
        if self.frame() != PyFrames::Equatorial {
            return Err(NEOSpyError::ValueError(
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
        if self.frame() != PyFrames::Equatorial {
            return Err(NEOSpyError::ValueError(
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
        if self.frame() != PyFrames::Ecliptic {
            return Err(NEOSpyError::ValueError(
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
        if self.frame() != PyFrames::Ecliptic {
            return Err(NEOSpyError::ValueError(
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
        let self_vec = self.0.to_ecliptic_raw();
        let other_vec = other.into_ecliptic_raw();
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
        Self(self.0.into_frame(target_frame))
    }

    /// X coordinate in au.
    #[getter]
    pub fn x(&self) -> f64 {
        self.raw()[0]
    }

    /// Y coordinate in au.
    #[getter]
    pub fn y(&self) -> f64 {
        self.raw()[1]
    }

    /// Z coordinate in au.
    #[getter]
    pub fn z(&self) -> f64 {
        self.raw()[2]
    }

    pub fn rotate_around(&self, other: VectorLike, angle: f64) -> Self {
        let self_vec = self.0.to_ecliptic_raw();
        let other_vec = other.into_ecliptic_raw();
        let rotated = self_vec.rotate_around(other_vec, angle.to_radians());
        Self::new(rotated.into(), self.frame())
    }

    pub fn __repr__(&self) -> String {
        // 1e-12 AU is about 15cm, this seems like a reasonable printing resolution
        let x = (self.raw()[0] * 1e12).round() / 1e12 + 0.0;
        let y = (self.raw()[1] * 1e12).round() / 1e12 + 0.0;
        let z = (self.raw()[2] * 1e12).round() / 1e12 + 0.0;
        format!("Vector([{:?}, {:?}, {:?}], {:?})", x, y, z, self.frame())
    }

    pub fn __sub__(&self, other: VectorLike) -> Self {
        let diff = self.0.to_ecliptic_raw() - &other.into_ecliptic_raw();
        PyVector(diff.into())
    }

    pub fn __add__(&self, other: VectorLike) -> Self {
        let diff = self.0.to_ecliptic_raw() + &other.into_ecliptic_raw();
        PyVector(diff.into())
    }

    pub fn __mul__(&self, other: f64) -> Self {
        let self_vec = Vector3::from(self.raw());
        Self::new((self_vec * other).into(), self.frame())
    }

    pub fn __truediv__(&self, other: f64) -> Self {
        let self_vec = Vector3::from(self.raw());
        Self::new((self_vec / other).into(), self.frame())
    }

    pub fn __neg__(&self) -> Self {
        Self::new([-self.x(), -self.y(), -self.z()], self.frame())
    }

    pub fn __len__(&self) -> usize {
        3
    }

    pub fn __getitem__(&self, idx: usize) -> PyResult<f64> {
        if idx >= 3 {
            return Err(PyErr::new::<exceptions::PyIndexError, _>(""));
        }
        Ok(self.0.raw()[idx])
    }

    fn __richcmp__(&self, other: VectorLike, op: CompareOp, py: Python<'_>) -> PyObject {
        let diff = self.__sub__(other).r();
        match op {
            CompareOp::Eq => (diff < 1e-12).into_py(py),
            CompareOp::Ne => (diff >= 1e-12).into_py(py),
            _ => py.NotImplemented(),
        }
    }
}
