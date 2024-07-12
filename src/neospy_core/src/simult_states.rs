//! Collections of [`State`] objects are grouped together into a [`SimultaneousStates`].
//! These primarily exist to allow for easy read/write to binary formats.

use crate::fov::FOV;
use crate::frames::{Ecliptic, Equatorial, Galactic, InertialFrame, FK4};
use crate::io::FileIO;
use crate::prelude::{NEOSpyError, State};
use serde::{Deserialize, Serialize};

/// Collection of [`State`] at the same time.
#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct SimultaneousStates<T: InertialFrame> {
    /// Collection of states
    pub states: Vec<State<T>>,

    /// Common JD time of all states
    pub jd: f64,

    /// Center ID of all states.
    pub center_id: i64,

    /// An optional field of view.
    pub fov: Option<FOV>,
}

impl FileIO for SimultaneousStates<Ecliptic> {}
impl FileIO for SimultaneousStates<Equatorial> {}
impl FileIO for SimultaneousStates<FK4> {}
impl FileIO for SimultaneousStates<Galactic> {}

impl<T: InertialFrame> SimultaneousStates<T> {
    /// Create a new Exact SimultaneousStates
    /// Simultaneous States occur at the same JD, which is defined by either the time
    /// in the optional fov, or the time of the first state.
    pub fn new_exact(states: Vec<State<T>>, fov: Option<FOV>) -> Result<Self, NEOSpyError> {
        if states.is_empty() {
            return Err(NEOSpyError::ValueError(
                "SimultaneousStates must contain at least one state.".into(),
            ));
        }
        let (mut jd, center_id) = {
            let first = states.first().unwrap();
            (first.jd, first.center_id)
        };

        if let Some(f) = &fov {
            jd = f.observer().jd
        }

        if !states
            .iter()
            .all(|state: &State<T>| state.center_id == center_id)
        {
            return Err(NEOSpyError::ValueError(
                "Center IDs do not match expected".into(),
            ));
        };
        if fov.is_none() && !states.iter().all(|state: &State<T>| state.jd == jd) {
            return Err(NEOSpyError::ValueError(
                "Epoch JDs do not match expected".into(),
            ));
        };
        Ok(SimultaneousStates {
            states,
            jd,
            center_id,
            fov,
        })
    }
}
