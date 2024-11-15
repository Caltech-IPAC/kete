//! Collections of [`State`] objects are grouped together into a [`SimultaneousStates`].
//! These primarily exist to allow for easy read/write to binary formats.

use crate::fov::FOV;
use crate::frames::{Ecliptic, Equatorial, Galactic, InertialFrame, FK4};
use crate::io::FileIO;
use crate::prelude::{Error, KeteResult, State};
use rayon::iter::{IntoParallelRefIterator, ParallelIterator};
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
    pub fn new_exact(mut states: Vec<State<T>>, fov: Option<FOV>) -> KeteResult<Self> {
        if states.is_empty() {
            return Err(Error::ValueError(
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
            .iter_mut()
            .all(|state: &mut State<T>| state.center_id == center_id)
        {
            return Err(Error::ValueError("Center IDs do not match expected".into()));
        };
        if fov.is_none() && states.iter_mut().any(|state: &mut State<T>| state.jd != jd) {
            return Err(Error::ValueError(
                "Epoch JDs do not match expected, this is only allowed if there is an associated FOV."
                    .into(),
            ));
        };
        Ok(SimultaneousStates {
            states,
            jd,
            center_id,
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

        let obs = fov.observer().clone().into_frame::<Equatorial>();

        Ok(self
            .states
            .par_iter()
            .map(|state| {
                let state = state.clone().into_frame::<Equatorial>();
                let pos_diff = state.pos - obs.pos;
                let vel_diff = state.vel - obs.vel;

                let d_ra = (pos_diff.x() * vel_diff.y() - pos_diff.y() * vel_diff.x())
                    / (pos_diff.x().powi(2) + pos_diff.y().powi(2));
                let r2 = pos_diff.norm_squared();

                let d_dec = (vel_diff.z() - pos_diff.z() * pos_diff.dot(&vel_diff) / r2)
                    / (r2 - pos_diff.z().powi(2)).sqrt();

                let (ra, dec) = pos_diff.to_ra_dec();

                [ra, dec, d_ra * dec.cos(), d_dec]
            })
            .collect())
    }
}
