//! Collections of [`State`] objects are grouped together into a [`SimultaneousStates`].
//! These primarily exist to allow for easy read/write to binary formats.

use crate::fov::FOV;
use crate::frames::to_lat_lon;
use crate::io::FileIO;
use crate::prelude::{Error, Frame, KeteResult, State};
use nalgebra::Vector3;
use rayon::iter::{IntoParallelRefIterator, ParallelIterator};
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
    pub fn new_exact(mut states: Vec<State>, fov: Option<FOV>) -> KeteResult<Self> {
        if states.is_empty() {
            return Err(Error::ValueError(
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
            return Err(Error::ValueError("Failed to change frames".into()));
        };
        if !states
            .iter_mut()
            .all(|state: &mut State| state.center_id == center_id)
        {
            return Err(Error::ValueError("Center IDs do not match expected".into()));
        };
        if fov.is_none() && states.iter_mut().any(|state: &mut State| state.jd != jd) {
            return Err(Error::ValueError(
                "Epoch JDs do not match expected, this is only allowed if there is an associated FOV."
                    .into(),
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

    /// Compute RA/Dec along with linearized rates.
    ///
    /// Returns a vector containing [ra, dec, ra' * cos(dec), dec'], all vectors
    /// are automatically cast into the equatorial frame.
    /// The returned RA rate is scaled by cos(dec) so that it is equivalent to a
    /// linear projection onto the observing plane.
    pub fn ra_dec_with_rates(&self) -> KeteResult<Vec<[f64; 4]>> {
        if self.fov.is_none() {
            return Err(Error::ValueError(
                "Field of view must be specified for the ra/dec to be computed.".into(),
            ));
        }
        let fov = self.fov.as_ref().unwrap();

        let mut obs = fov.observer().clone();
        obs.try_change_frame_mut(Frame::Equatorial)?;

        let obs_pos = Vector3::from(obs.pos);
        let obs_vel = Vector3::from(obs.vel);

        Ok(self
            .states
            .par_iter()
            .map(|state| {
                let mut state = state.clone();
                state.try_change_frame_mut(Frame::Equatorial).unwrap();
                let pos_diff = Vector3::from(state.pos) - obs_pos;
                let vel_diff = Vector3::from(state.vel) - obs_vel;

                let d_ra = (pos_diff.x * vel_diff.y - pos_diff.y * vel_diff.x)
                    / (pos_diff.x.powi(2) + pos_diff.y.powi(2));
                let r2 = pos_diff.norm_squared();

                let d_dec = (vel_diff.z - pos_diff.z * pos_diff.dot(&vel_diff) / r2)
                    / (r2 - pos_diff.z.powi(2)).sqrt();
                let (dec, ra) = to_lat_lon(pos_diff[0], pos_diff[1], pos_diff[2]);

                [ra, dec, d_ra * dec.cos(), d_dec]
            })
            .collect())
    }
}
