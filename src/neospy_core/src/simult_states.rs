//! Collections of [`State`] objects are grouped together into a [`SimultaneousStates`].
//! These primarily exist to allow for easy read/write to binary formats.

use crate::fov::FOV;
use crate::io::FileIO;
use crate::prelude::{Frame, NEOSpyError, State};
use serde::{Deserialize, Serialize};

/// Collection of [`State`] at the same time.
#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct SimultaneousStates {
    /// Collection of states
    pub states: Vec<State>,

    /// Common JD time of all states
    pub jd: f64,

    /// Center ID of all states.
    pub center_id: i64,

    /// Frame of reference for all states.
    pub frame: Frame,

    /// An optional field of view.
    pub fov: Option<FOV>,
}

impl FileIO for SimultaneousStates {}

impl SimultaneousStates {
    /// Create a new Exact SimultaneousStates
    /// Simultaneous States occur at the same JD, which is defined by either the time
    /// in the optional fov, or the time of the first state.
    pub fn new_exact(mut states: Vec<State>, fov: Option<FOV>) -> Result<Self, NEOSpyError> {
        if states.is_empty() {
            return Err(NEOSpyError::ValueError(
                "SimultaneousStates must contain at least one state.".into(),
            ));
        }
        let (mut jd, frame, center_id) = {
            let first = states.first().unwrap();
            (first.jd, first.frame, first.center_id)
        };

        if let Some(f) = &fov {
            jd = f.observer().jd
        }

        if !states
            .iter_mut()
            .map(|state| state.try_change_frame_mut(frame))
            .all(|x| x.is_ok())
        {
            return Err(NEOSpyError::ValueError("Failed to change frames".into()));
        };
        if !states
            .iter_mut()
            .all(|state: &mut State| state.center_id == center_id)
        {
            return Err(NEOSpyError::ValueError(
                "Center IDs do not match expected".into(),
            ));
        };
        if fov.is_none() && !states.iter_mut().all(|state: &mut State| state.jd == jd) {
            return Err(NEOSpyError::ValueError(
                "Epoch JDs do not match expected".into(),
            ));
        };
        Ok(SimultaneousStates {
            states,
            jd,
            center_id,
            frame,
            fov,
        })
    }
}
