use nalgebra::Vector3;
use neospy_core::fov;
use neospy_core::fov::{FovLike, SkyPatch};
use pyo3::{exceptions, prelude::*};

use crate::vector::VectorLike;
use crate::{state::PyState, vector::Vector};

/// Field of view of a WISE CMOS chip.
/// Since all WISE CMOS see the same patch of sky, there is no differentiation
/// of the individual wavelength bands.
#[pyclass(module = "neospy", frozen, name = "WiseCmos")]
#[derive(Clone, Debug)]
#[allow(clippy::upper_case_acronyms)]
pub struct PyWiseCmos(pub fov::WiseCmos);

/// Field of view of a NEOS CMOS chip.
#[pyclass(module = "neospy", frozen, name = "NeosCmos")]
#[derive(Clone, Debug)]
#[allow(clippy::upper_case_acronyms)]
pub struct PyNeosCmos(pub fov::NeosCmos);

/// Field of view of a Single ZTF chips/quad combination.
#[pyclass(module = "neospy", frozen, name = "ZtfCcdQuad")]
#[derive(Clone, Debug)]
#[allow(clippy::upper_case_acronyms)]
pub struct PyZtfCcdQuad(pub fov::ZtfCcdQuad);

/// Field of view of all 64 ZTF chips/quad combinations.
/// This is a meta collection of individual ZTF CCD Quad FOVs.
#[pyclass(module = "neospy", frozen, name = "ZtfField", sequence)]
#[derive(Clone, Debug)]
#[allow(clippy::upper_case_acronyms)]
pub struct PyZtfField(pub fov::ZtfField);

/// Generic Rectangular Field of view.
#[pyclass(module = "neospy", frozen, name = "RectangleFOV")]
#[derive(Clone, Debug)]
pub struct PyGenericRectangle(pub fov::GenericRectangle);

/// Field of views supported by the python interface
#[derive(Clone, FromPyObject)]
#[allow(clippy::upper_case_acronyms)]
pub enum AllowedFOV {
    WISE(PyWiseCmos),
    NEOS(PyNeosCmos),
    Rectangle(PyGenericRectangle),
    ZTF(PyZtfCcdQuad),
    ZTFField(PyZtfField),
}

impl AllowedFOV {
    pub fn jd(&self) -> f64 {
        match self {
            AllowedFOV::NEOS(fov) => fov.0.observer().jd,
            AllowedFOV::WISE(fov) => fov.0.observer().jd,
            AllowedFOV::Rectangle(fov) => fov.0.observer().jd,
            AllowedFOV::ZTF(fov) => fov.0.observer().jd,
            AllowedFOV::ZTFField(fov) => fov.0.observer().jd,
        }
    }

    pub fn get_fov(self, idx: Option<usize>) -> fov::FOV {
        let idx = idx.unwrap_or_default();
        match self {
            AllowedFOV::WISE(fov) => fov.0.get_fov(idx),
            AllowedFOV::Rectangle(fov) => fov.0.get_fov(idx),
            AllowedFOV::NEOS(fov) => fov.0.get_fov(idx),
            AllowedFOV::ZTF(fov) => fov.0.get_fov(idx),
            AllowedFOV::ZTFField(fov) => fov.0.get_fov(idx),
        }
    }

    pub fn unwrap(self) -> fov::FOV {
        match self {
            AllowedFOV::WISE(fov) => fov::FOV::Wise(fov.0),
            AllowedFOV::Rectangle(fov) => fov::FOV::GenericRectangle(fov.0),
            AllowedFOV::NEOS(fov) => fov::FOV::NeosCmos(fov.0),
            AllowedFOV::ZTF(fov) => fov::FOV::ZtfCcdQuad(fov.0),
            AllowedFOV::ZTFField(fov) => fov::FOV::ZtfField(fov.0),
        }
    }

    pub fn __repr__(self) -> String {
        match self {
            AllowedFOV::WISE(fov) => fov.__repr__(),
            AllowedFOV::Rectangle(fov) => fov.__repr__(),
            AllowedFOV::NEOS(fov) => fov.__repr__(),
            AllowedFOV::ZTF(fov) => fov.__repr__(),
            AllowedFOV::ZTFField(fov) => fov.__repr__(),
        }
    }
}

impl IntoPy<PyObject> for AllowedFOV {
    fn into_py(self, py: Python) -> PyObject {
        match self {
            Self::WISE(fov) => fov.into_py(py),
            Self::NEOS(fov) => fov.into_py(py),
            Self::Rectangle(fov) => fov.into_py(py),
            Self::ZTF(fov) => fov.into_py(py),
            Self::ZTFField(fov) => fov.into_py(py),
        }
    }
}

impl From<fov::FOV> for AllowedFOV {
    fn from(value: fov::FOV) -> Self {
        match value {
            fov::FOV::Wise(fov) => AllowedFOV::WISE(PyWiseCmos(fov)),
            fov::FOV::ZtfCcdQuad(fov) => AllowedFOV::ZTF(PyZtfCcdQuad(fov)),
            fov::FOV::NeosCmos(fov) => AllowedFOV::NEOS(PyNeosCmos(fov)),
            fov::FOV::GenericRectangle(fov) => AllowedFOV::Rectangle(PyGenericRectangle(fov)),
            fov::FOV::ZtfField(fov) => AllowedFOV::ZTFField(PyZtfField(fov)),
            _ => {
                unimplemented!("Python interface doesn't support this FOV.")
            }
        }
    }
}

#[pymethods]
impl PyWiseCmos {
    #[new]
    pub fn new(
        pointing: VectorLike,
        rotation: f64,
        observer: PyState,
        frame_num: usize,
        scan_id: String,
    ) -> Self {
        let pointing = pointing.into_vector(observer.frame());
        let pointing = pointing.raw.into();
        let scan_id = scan_id.into();
        PyWiseCmos(fov::WiseCmos::new(
            pointing,
            rotation.to_radians(),
            observer.0,
            frame_num,
            scan_id,
        ))
    }

    /// Position of the observer in this FOV.
    #[getter]
    pub fn observer(&self) -> PyState {
        self.0.observer().clone().into()
    }

    /// Direction that the observer is looking.
    #[getter]
    pub fn pointing(&self) -> Vector {
        let pointing = self.0.patch.pointing().into_inner().into();
        Vector::new(pointing, self.0.patch.frame.into())
    }

    /// WISE Frame number.
    #[getter]
    pub fn frame_num(&self) -> usize {
        self.0.frame_num
    }

    /// WISE Scan ID.
    #[getter]
    pub fn scan_id(&self) -> String {
        self.0.scan_id.to_string()
    }

    /// Rotation angle of the FOV in degrees.
    #[getter]
    pub fn rotation(&self) -> f64 {
        self.0.rotation.to_degrees()
    }

    fn __repr__(&self) -> String {
        format!(
            "WiseCmos(pointing={}, rotation={}, observer={}, frame_num={}, scan_id={:?})",
            self.pointing().__repr__(),
            self.rotation(),
            self.observer().__repr__(),
            self.frame_num(),
            self.scan_id()
        )
    }
}

#[pymethods]
impl PyGenericRectangle {
    #[new]
    pub fn new(
        pointing: VectorLike,
        rotation: f64,
        observer: PyState,
        lon_width: f64,
        lat_width: f64,
    ) -> Self {
        let pointing = pointing.into_vector(observer.frame());
        PyGenericRectangle(fov::GenericRectangle::new(
            pointing.raw.into(),
            rotation.to_radians(),
            lon_width.to_radians(),
            lat_width.to_radians(),
            observer.0,
        ))
    }

    /// Construct a new Rectangle FOV from the corners.
    /// The corners must be provided in clockwise order.
    ///
    /// Parameters
    /// ----------
    /// corners :
    ///     4 Vectors which represent the corners of the FOV, these must be provided in clockwise order.
    /// observer :
    ///     Position of the observer as a State.
    #[staticmethod]
    pub fn from_corners(corners: [VectorLike; 4], observer: PyState) -> Self {
        let corners: [Vector3<f64>; 4] = corners.map(|x| x.into_vec(observer.frame()));
        PyGenericRectangle(fov::GenericRectangle::from_corners(corners, observer.0))
    }

    /// The observer State.
    #[getter]
    pub fn observer(&self) -> PyState {
        self.0.observer().clone().into()
    }

    /// Direction that the observer is looking.
    #[getter]
    pub fn pointing(&self) -> Vector {
        Vector::new(
            self.0.patch.pointing().into_inner().into(),
            self.0.observer().frame.into(),
        )
    }

    /// The longitudinal width of the FOV.
    #[getter]
    pub fn lon_width(&self) -> f64 {
        self.0.lon_width().to_degrees()
    }

    /// The Latitudinal width of the FOV.
    #[getter]
    pub fn lat_width(&self) -> f64 {
        self.0.lat_width().to_degrees()
    }

    fn __repr__(&self) -> String {
        format!(
            "GenericRectangle(pointing={}, observer={}, lon_width={}, lat_width={})",
            self.pointing().__repr__(),
            self.observer().__repr__(),
            self.lon_width(),
            self.lat_width(),
        )
    }
}

#[pymethods]
#[allow(clippy::too_many_arguments)]
impl PyNeosCmos {
    /// Construct a new NEOS FOV.
    ///
    /// Parameters
    /// ----------
    /// pointing :
    ///     Vector defining the center of the FOV.
    /// rotation :
    ///     Rotation of the FOV in degrees.
    /// observer :
    ///     State of the observer.
    /// side_id :
    ///     Side ID indicating where we are in the survey.
    /// stack_id :
    ///     Stack ID indicating where we are in the survey.
    /// quad_id :
    ///     Quad ID indicating where we are in the survey.
    /// loop_id :
    ///     Loop ID indicating where we are in the survey.
    /// subloop_id :
    ///     Subloop ID indicating where we are in the survey.
    /// exposure_id :
    ///     Exposure number indicating where we are in the survey.
    /// cmos_id :
    ///     Which chip of the target band this represents.
    /// band :
    ///     Band, can be either 1 or 2 to represent NC1/NC2.
    #[new]
    pub fn new(
        pointing: VectorLike,
        rotation: f64,
        observer: PyState,
        side_id: u16,
        stack_id: u8,
        quad_id: u8,
        loop_id: u8,
        subloop_id: u8,
        exposure_id: u8,
        cmos_id: u8,
        band: u8,
    ) -> Self {
        let pointing = pointing.into_vector(observer.frame());
        let pointing = pointing.raw.into();
        PyNeosCmos(fov::NeosCmos::new(
            pointing,
            rotation.to_radians(),
            observer.0,
            side_id,
            stack_id,
            quad_id,
            loop_id,
            subloop_id,
            exposure_id,
            cmos_id,
            band,
        ))
    }

    /// The observer State.
    #[getter]
    pub fn observer(&self) -> PyState {
        self.0.observer().clone().into()
    }

    /// Direction that the observer is looking.
    #[getter]
    pub fn pointing(&self) -> Vector {
        Vector::new(
            self.0.patch.pointing().into_inner().into(),
            self.0.observer().frame.into(),
        )
    }

    /// Metadata about where this FOV is in the Survey.
    #[getter]
    pub fn side_id(&self) -> u16 {
        self.0.side_id
    }

    /// Metadata about where this FOV is in the Survey.
    #[getter]
    pub fn stack_id(&self) -> u8 {
        self.0.stack_id
    }

    /// Metadata about where this FOV is in the Survey.
    #[getter]
    pub fn quad_id(&self) -> u8 {
        self.0.quad_id
    }

    /// Metadata about where this FOV is in the Survey.
    #[getter]
    pub fn loop_id(&self) -> u8 {
        self.0.loop_id
    }

    /// Metadata about where this FOV is in the Survey.
    #[getter]
    pub fn subloop_id(&self) -> u8 {
        self.0.subloop_id
    }

    /// Metadata about where this FOV is in the Survey.
    #[getter]
    pub fn exposure_id(&self) -> u8 {
        self.0.exposure_id
    }

    /// Rotation angle of the FOV in degrees.
    #[getter]
    pub fn rotation(&self) -> f64 {
        self.0.rotation.to_degrees()
    }

    fn __repr__(&self) -> String {
        format!(
            "NeosCmos(pointing={}, rotation={}, observer={}, side_id={}, stack_id={}, quad_id={}, loop_id={}, subloop_id={}, exposure_id={})",
            self.pointing().__repr__(),
            self.rotation(),
            self.observer().__repr__(),
            self.side_id(),
            self.stack_id(),
            self.quad_id(),
            self.loop_id(),
            self.subloop_id(),
            self.exposure_id()
        )
    }
}

#[pymethods]
#[allow(clippy::too_many_arguments)]
impl PyZtfCcdQuad {
    /// Construct a new ZTF CCD FOV from the corners.
    /// The corners must be provided in clockwise order.
    ///
    /// Parameters
    /// ----------
    /// corners :
    ///     4 Vectors which represent the corners of the FOV, these must be provided in clockwise order.
    /// observer :
    ///     Position of the observer as a State.
    /// field :
    ///     Field number of the survey.
    /// filefracday :
    ///     Which fraction of a day was this FOV captured.
    /// ccdid :
    ///     CCD ID describes which of the 16 CCDs this represents.
    /// filtercode :
    ///     Which filter was used for this exposure.
    /// imgtypecode :
    ///     Type code describing the data product of this field.
    /// qid :
    ///     Which quadrant of the CCD does this FOV represent.
    /// maglimit :
    ///     Effective magnitude limit of this exposure.
    /// fid :
    ///     The FID of this exposure.
    #[new]
    pub fn new(
        corners: [VectorLike; 4],
        observer: PyState,
        field: u32,
        filefracday: u64,
        ccdid: u8,
        filtercode: String,
        imgtypecode: String,
        qid: u8,
        maglimit: f64,
        fid: usize,
    ) -> Self {
        let corners = corners
            .into_iter()
            .map(|v| v.into_vector(observer.frame()).raw.into())
            .collect::<Vec<_>>()
            .try_into()
            .unwrap();
        PyZtfCcdQuad(fov::ZtfCcdQuad::new(
            corners,
            observer.0,
            field,
            filefracday,
            ccdid,
            filtercode.into(),
            imgtypecode.into(),
            qid,
            maglimit,
            fid,
        ))
    }

    #[getter]
    pub fn observer(&self) -> PyState {
        self.0.observer().clone().into()
    }

    /// Metadata about where this FOV is in the Survey.
    #[getter]
    pub fn field(&self) -> u32 {
        self.0.field
    }

    /// Direction that the observer is looking.
    #[getter]
    pub fn pointing(&self) -> Vector {
        Vector::new(
            self.0.patch.pointing().into_inner().into(),
            self.0.observer().frame.into(),
        )
    }

    /// Metadata about where this FOV is in the Survey.
    #[getter]
    pub fn filefracday(&self) -> u64 {
        self.0.filefracday
    }

    /// Metadata about where this FOV is in the Survey.
    #[getter]
    pub fn ccdid(&self) -> u8 {
        self.0.ccdid
    }

    /// Metadata about where this FOV is in the Survey.
    #[getter]
    pub fn filtercode(&self) -> String {
        self.0.filtercode.to_string()
    }

    /// Metadata about where this FOV is in the Survey.
    #[getter]
    pub fn imgtypecode(&self) -> String {
        self.0.imgtypecode.to_string()
    }

    /// Metadata about where this FOV is in the Survey.
    #[getter]
    pub fn qid(&self) -> u8 {
        self.0.qid
    }

    /// Magnitude limit of this exposure.
    #[getter]
    pub fn maglimit(&self) -> f64 {
        self.0.maglimit
    }

    /// Metadata about where this FOV is in the Survey.
    #[getter]
    pub fn fid(&self) -> usize {
        self.0.fid
    }

    fn __repr__(&self) -> String {
        format!(
            "ZtfCcdQuad(observer={}, filefracday={}, ccdid={}, filtercode={:?}, imgtypecode={:?}, qid={}, maglimit={}, fid={})",
            self.observer().__repr__(),
            self.filefracday(),
            self.ccdid(),
            self.filtercode(),
            self.imgtypecode(),
            self.qid(),
            self.maglimit(),
            self.fid()
        )
    }
}

#[pymethods]
#[allow(clippy::too_many_arguments)]
impl PyZtfField {
    /// Representation of an entire ZTF Field, made up of up to 64 ZTF CCD FOVs.
    ///
    /// Parameters
    /// ----------
    /// ztf_ccd_fields :
    ///     List containing all of the CCD FOVs.
    ///     These must have matching metadata.
    #[new]
    pub fn new(ztf_ccd_fields: Vec<PyZtfCcdQuad>) -> Self {
        PyZtfField(fov::ZtfField::new(ztf_ccd_fields.into_iter().map(|x| x.0).collect()).unwrap())
    }

    /// State of the observer for this FOV.
    #[getter]
    pub fn observer(&self) -> PyState {
        self.0.observer().clone().into()
    }

    /// Metadata about where this FOV is in the Survey.
    #[getter]
    pub fn field(&self) -> u32 {
        self.0.field
    }

    /// Metadata about where this FOV is in the Survey.
    #[getter]
    pub fn filtercode(&self) -> String {
        self.0.filtercode.to_string()
    }

    /// Metadata about where this FOV is in the Survey.
    #[getter]
    pub fn imgtypecode(&self) -> String {
        self.0.imgtypecode.to_string()
    }

    /// Metadata about where this FOV is in the Survey.
    #[getter]
    pub fn fid(&self) -> usize {
        self.0.fid
    }

    /// Return all of the individual CCD quads present in this field.
    #[getter]
    pub fn ccd_quads(&self) -> Vec<PyZtfCcdQuad> {
        (0..self.0.n_patches()).map(|idx| PyZtfCcdQuad(match self.0.get_fov(idx) {
            fov::FOV::ZtfCcdQuad(fov) => fov,
            _ => unreachable!(),
        })).collect()
    }

    pub fn __len__(&self) -> usize {
        self.0.n_patches()
    }

    /// Retrieve a specific CCD FOV from the contained list of FOVs.
    pub fn __getitem__(&self, idx: usize) -> PyResult<PyZtfCcdQuad> {
        if idx >= self.__len__() {
            return Err(PyErr::new::<exceptions::PyIndexError, _>(""));
        }

        Ok(PyZtfCcdQuad(match self.0.get_fov(idx) {
            fov::FOV::ZtfCcdQuad(fov) => fov,
            _ => unreachable!(),
        }))
    }

    fn __repr__(&self) -> String {
        format!(
            "ZtfField(ccd_quads=<{} frames>, observer={}, filtercode={:?}, imgtypecode={:?}, fid={})",
            self.0.n_patches(),
            self.observer().__repr__(),
            self.filtercode(),
            self.imgtypecode(),
            self.fid()
        )
    }
}
