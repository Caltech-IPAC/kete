use crate::elements::PyCometElements;
use crate::frame::*;
use crate::vector::*;
use neospy_core::frames::Equatorial;
use neospy_core::prelude;
use pyo3::prelude::*;

/// Representation of the state of an object at a specific moment in time.
///
#[pyclass(frozen, module = "neospy", name = "State")]
#[derive(Clone, Debug)]
pub struct PyState(pub prelude::State<Equatorial>);

impl From<prelude::State<Equatorial>> for PyState {
    fn from(value: prelude::State<Equatorial>) -> Self {
        Self(value)
    }
}

#[pymethods]
impl PyState {
    #[new]
    #[pyo3(signature = (desig, jd, pos, vel, frame=None, center_id=10))]
    pub fn new(
        desig: Option<String>,
        jd: f64,
        pos: VectorLike,
        vel: VectorLike,
        frame: Option<PyFrames>,
        center_id: i64,
    ) -> Self {
        let desig = desig.map(|x| prelude::Desig::Name(x));

        // if no frame is provided, but pos or vel have a frame, use that one.
        let frame = frame.unwrap_or({
            if let VectorLike::PyVec(v) = &pos {
                v.frame()
            } else if let VectorLike::PyVec(v) = &vel {
                v.frame()
            } else {
                PyFrames::Ecliptic
            }
        });

        // if frames dont match, force them all into the target frame.
        let pos = pos.into_vec(frame);
        let vel = vel.into_vec(frame);

        let state = prelude::State::new(desig, jd, pos.into(), vel.into(), center_id);
        Self(state)
    }

    /// Change the coordinate frame of the object to the target frame.
    ///
    /// Parameters
    /// ----------
    /// frame : Frames
    ///     New frame to change to.
    pub fn change_frame(&self, frame: PyFrames) -> PyResult<Self> {
        let mut state = self.0.clone();
        state.into_frame()?;
        Ok(Self(state))
    }

    /// Change the center ID of the state from the current state to the target state.
    ///
    /// If the desired state is not a known NAIF id this will raise an exception.
    pub fn change_center(&self, naif_id: i64) -> PyResult<Self> {
        let mut state = self.0.clone();
        let spk = neospy_core::prelude::get_spk_singleton()
            .try_read()
            .unwrap();
        spk.try_change_center(&mut state, naif_id)?;
        Ok(Self(state))
    }

    /// Convert state to the Ecliptic Frame.
    #[getter]
    pub fn as_ecliptic(&self) -> PyResult<Self> {
        self.change_frame(PyFrames::Ecliptic)
    }

    /// Convert state to the Equatorial Frame.
    #[getter]
    pub fn as_equatorial(&self) -> PyResult<Self> {
        self.change_frame(PyFrames::Equatorial)
    }

    /// Convert state to the Galactic Frame.
    #[getter]
    pub fn as_galactic(&self) -> PyResult<Self> {
        self.change_frame(PyFrames::Galactic)
    }

    /// Convert state to the FK4 Frame.
    #[getter]
    pub fn as_fk4(&self) -> PyResult<Self> {
        self.change_frame(PyFrames::FK4)
    }

    /// JD of the object's state in TDB scaled time.
    #[getter]
    pub fn jd(&self) -> f64 {
        self.0.jd
    }

    /// Position of the object in AU with respect to the central object.
    #[getter]
    pub fn pos(&self) -> PyVector {
        let v = self.0.pos;
        PyVector::new(v, self.frame())
    }

    /// Velocity of the object in AU/Day.
    #[getter]
    pub fn vel(&self) -> PyVector {
        let v = self.0.vel;
        PyVector::new(v, self.frame())
    }

    /// Frame of reference used to define the coordinate system.
    #[getter]
    pub fn frame(&self) -> PyFrames {
        self.0.frame.into()
    }

    /// Central ID of the object used as reference for the coordinate frame.
    #[getter]
    pub fn center_id(&self) -> i64 {
        self.0.center_id
    }

    /// Cometary orbital elements of the state.
    #[getter]
    pub fn elements(&self) -> PyCometElements {
        PyCometElements::from_state(self)
    }

    /// Designation of the object if defined.
    #[getter]
    pub fn desig(&self) -> String {
        match &self.0.desig {
            prelude::Desig::Name(s) => s.clone(),
            prelude::Desig::Naif(s) => {
                neospy_core::spice::try_name_from_id(*s).unwrap_or(s.to_string())
            }
            prelude::Desig::Perm(s) => format!("{:?}", s),
            prelude::Desig::Prov(s) => s.clone(),
            prelude::Desig::Empty => "None".into(),
        }
    }

    pub fn __repr__(&self) -> String {
        format!(
            "State(desig={:?}, jd={:?}, pos={:?}, vel={:?}, frame={:?}, center_id={:?})",
            self.desig(),
            self.jd(),
            self.pos().raw,
            self.vel().raw,
            self.frame(),
            self.center_id()
        )
    }
}
