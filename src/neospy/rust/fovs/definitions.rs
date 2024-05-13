use neospy_core::fov;
use neospy_core::fov::{FovLike, SkyPatch};
use pyo3::prelude::*;

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

/// Field of view of a NEOS Visit.
#[pyclass(module = "neospy", frozen, name = "NeosVisit")]
#[derive(Clone, Debug)]
#[allow(clippy::upper_case_acronyms)]
pub struct PyNeosVisit(pub fov::NeosVisit);

/// Field of view of a Single ZTF chips/quad combination.
#[pyclass(module = "neospy", frozen, name = "ZtfCcdQuad")]
#[derive(Clone, Debug)]
#[allow(clippy::upper_case_acronyms)]
pub struct PyZtfCcdQuad(pub fov::ZtfCcdQuad);

/// Field of view of all 64 ZTF chips/quad combinations.
/// This is a meta collection of individual ZTF CCD Quad FOVs.
#[pyclass(module = "neospy", frozen, name = "ZtfField")]
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
    NEOSVisit(PyNeosVisit),
}

impl AllowedFOV {
    pub fn jd(&self) -> f64 {
        match self {
            AllowedFOV::NEOS(fov) => fov.0.observer().jd,
            AllowedFOV::WISE(fov) => fov.0.observer().jd,
            AllowedFOV::Rectangle(fov) => fov.0.observer().jd,
            AllowedFOV::ZTF(fov) => fov.0.observer().jd,
            AllowedFOV::ZTFField(fov) => fov.0.observer().jd,
            AllowedFOV::NEOSVisit(fov) => fov.0.observer().jd,
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
            AllowedFOV::NEOSVisit(fov) => fov.0.get_fov(idx),
        }
    }

    pub fn unwrap(self) -> fov::FOV {
        match self {
            AllowedFOV::WISE(fov) => fov::FOV::Wise(fov.0),
            AllowedFOV::Rectangle(fov) => fov::FOV::GenericRectangle(fov.0),
            AllowedFOV::NEOS(fov) => fov::FOV::NeosCmos(fov.0),
            AllowedFOV::ZTF(fov) => fov::FOV::ZtfCcdQuad(fov.0),
            AllowedFOV::ZTFField(fov) => fov::FOV::ZtfField(fov.0),
            AllowedFOV::NEOSVisit(fov) => fov::FOV::NeosVisit(fov.0),
        }
    }

    pub fn __repr__(self) -> String {
        match self {
            AllowedFOV::WISE(fov) => fov.__repr__(),
            AllowedFOV::Rectangle(fov) => fov.__repr__(),
            AllowedFOV::NEOS(fov) => fov.__repr__(),
            AllowedFOV::ZTF(fov) => fov.__repr__(),
            AllowedFOV::ZTFField(fov) => fov.__repr__(),
            AllowedFOV::NEOSVisit(fov) => fov.__repr__(),
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
            Self::NEOSVisit(fov) => fov.into_py(py),
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
            fov::FOV::NeosVisit(fov) => AllowedFOV::NEOSVisit(PyNeosVisit(fov)),
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
        let pointing = pointing.into_vector(crate::frame::PyFrames::Ecliptic);
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

    #[getter]
    pub fn observer(&self) -> PyState {
        self.0.observer().clone().into()
    }

    #[getter]
    pub fn pointing(&self) -> Vector {
        let pointing = self.0.patch.pointing().into_inner().into();
        Vector::new(pointing, self.0.observer().frame.into())
    }

    #[getter]
    pub fn frame_num(&self) -> usize {
        self.0.frame_num
    }

    #[getter]
    pub fn scan_id(&self) -> String {
        self.0.scan_id.to_string()
    }

    #[getter]
    pub fn rotation(&self) -> f64 {
        self.0.rotation
    }

    fn __repr__(&self) -> String {
        format!(
            "WISEFOV(pointing={}, rotation={}, observer={}, frame_num={}, scan_id={:?})",
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
        let pointing = pointing.into_vector(crate::frame::PyFrames::Ecliptic);
        PyGenericRectangle(fov::GenericRectangle::new(
            pointing.raw.into(),
            rotation.to_radians(),
            lon_width.to_radians(),
            lat_width.to_radians(),
            observer.0,
        ))
    }

    #[getter]
    pub fn observer(&self) -> PyState {
        self.0.observer().clone().into()
    }

    #[getter]
    pub fn pointing(&self) -> Vector {
        Vector::new(
            self.0.patch.pointing().into_inner().into(),
            self.0.observer().frame.into(),
        )
    }

    #[getter]
    pub fn rotation(&self) -> f64 {
        self.0.rotation.to_degrees()
    }

    #[getter]
    pub fn lon_width(&self) -> f64 {
        self.0.lon_width().to_degrees()
    }

    #[getter]
    pub fn lat_width(&self) -> f64 {
        self.0.lat_width().to_degrees()
    }

    fn __repr__(&self) -> String {
        format!(
            "RectangleFOV(pointing={}, rotation={}, observer={}, lon_width={}, lat_width={})",
            self.pointing().__repr__(),
            self.rotation(),
            self.observer().__repr__(),
            self.lon_width(),
            self.lat_width(),
        )
    }
}

#[pymethods]
#[allow(clippy::too_many_arguments)]
impl PyNeosCmos {
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
        let pointing = pointing.into_vector(crate::frame::PyFrames::Ecliptic);
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

    #[getter]
    pub fn observer(&self) -> PyState {
        self.0.observer().clone().into()
    }

    #[getter]
    pub fn pointing(&self) -> Vector {
        Vector::new(
            self.0.patch.pointing().into_inner().into(),
            self.0.observer().frame.into(),
        )
    }

    #[getter]
    pub fn side_id(&self) -> u16 {
        self.0.side_id
    }

    #[getter]
    pub fn stack_id(&self) -> u8 {
        self.0.stack_id
    }

    #[getter]
    pub fn quad_id(&self) -> u8 {
        self.0.quad_id
    }

    #[getter]
    pub fn loop_id(&self) -> u8 {
        self.0.loop_id
    }

    #[getter]
    pub fn subloop_id(&self) -> u8 {
        self.0.subloop_id
    }

    #[getter]
    pub fn exposure_id(&self) -> u8 {
        self.0.exposure_id
    }

    #[getter]
    pub fn rotation(&self) -> f64 {
        self.0.rotation
    }

    fn __repr__(&self) -> String {
        format!(
            "NEOSFOV(pointing={}, rotation={}, observer={}, side_id={}, stack_id={}, quad_id={}, loop_id={}, subloop_id={}, exposure_id={})",
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
impl PyNeosVisit {
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
        let pointing = pointing.into_vector(crate::frame::PyFrames::Ecliptic);
        let pointing = pointing.raw.into();
        PyNeosVisit(fov::NeosVisit::new(
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

    #[getter]
    pub fn observer(&self) -> PyState {
        self.0.observer().clone().into()
    }

    #[getter]
    pub fn pointing(&self) -> Vector {
        Vector::new(
            self.0.patch.pointing().into_inner().into(),
            self.0.observer().frame.into(),
        )
    }

    #[getter]
    pub fn side_id(&self) -> u16 {
        self.0.side_id
    }

    #[getter]
    pub fn stack_id(&self) -> u8 {
        self.0.stack_id
    }

    #[getter]
    pub fn quad_id(&self) -> u8 {
        self.0.quad_id
    }

    #[getter]
    pub fn loop_id(&self) -> u8 {
        self.0.loop_id
    }

    #[getter]
    pub fn subloop_id(&self) -> u8 {
        self.0.subloop_id
    }

    #[getter]
    pub fn exposure_id(&self) -> u8 {
        self.0.exposure_id
    }

    #[getter]
    pub fn rotation(&self) -> f64 {
        self.0.rotation
    }

    fn __repr__(&self) -> String {
        format!(
            "NEOSVisit(pointing={}, rotation={}, observer={}, side_id={}, stack_id={}, quad_id={}, loop_id={}, subloop_id={}, exposure_id={})",
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
            .map(|v| v.into_vector(crate::frame::PyFrames::Ecliptic).raw.into())
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

    #[getter]
    pub fn field(&self) -> u32 {
        self.0.field
    }

    #[getter]
    pub fn ra(&self) -> f64 {
        let pointing = self.pointing();
        pointing.as_equatorial().ra().unwrap()
    }

    #[getter]
    pub fn pointing(&self) -> Vector {
        Vector::new(
            self.0.patch.pointing().into_inner().into(),
            self.0.observer().frame.into(),
        )
    }

    #[getter]
    pub fn dec(&self) -> f64 {
        let pointing = self.pointing();
        pointing.as_equatorial().dec().unwrap()
    }

    #[getter]
    pub fn filefracday(&self) -> u64 {
        self.0.filefracday
    }

    #[getter]
    pub fn ccdid(&self) -> u8 {
        self.0.ccdid
    }

    #[getter]
    pub fn filtercode(&self) -> String {
        self.0.filtercode.to_string()
    }

    #[getter]
    pub fn imgtypecode(&self) -> String {
        self.0.imgtypecode.to_string()
    }

    #[getter]
    pub fn qid(&self) -> u8 {
        self.0.qid
    }

    #[getter]
    pub fn maglimit(&self) -> f64 {
        self.0.maglimit
    }

    #[getter]
    pub fn fid(&self) -> usize {
        self.0.fid
    }

    fn __repr__(&self) -> String {
        format!(
            "ZTFFOV(ra={}, dec={}, observer={}, filefracday={}, ccdid={}, filtercode={}, imgtypecode={}, qid={}, maglimit={}, fid={})",
            self.ra(),
            self.dec(),
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
    #[new]
    pub fn new(ztf_ccd_fields: Vec<PyZtfCcdQuad>) -> Self {
        PyZtfField(fov::ZtfField::new(ztf_ccd_fields.into_iter().map(|x| x.0).collect()).unwrap())
    }

    #[getter]
    pub fn observer(&self) -> PyState {
        self.0.observer().clone().into()
    }

    #[getter]
    pub fn field(&self) -> u32 {
        self.0.field
    }

    #[getter]
    pub fn filtercode(&self) -> String {
        self.0.filtercode.to_string()
    }

    #[getter]
    pub fn imgtypecode(&self) -> String {
        self.0.imgtypecode.to_string()
    }

    #[getter]
    pub fn fid(&self) -> usize {
        self.0.fid
    }

    pub fn get_ccd(&self, idx: usize) -> PyZtfCcdQuad {
        PyZtfCcdQuad(match self.0.get_fov(idx) {
            fov::FOV::ZtfCcdQuad(fov) => fov,
            _ => panic!(),
        })
    }

    fn __repr__(&self) -> String {
        format!(
            "ZTFFOV(ccd_quads=<{} frames>, observer={},filtercode={}, imgtypecode={}, fid={})",
            self.0.n_patches(),
            self.observer().__repr__(),
            self.filtercode(),
            self.imgtypecode(),
            self.fid()
        )
    }
}
