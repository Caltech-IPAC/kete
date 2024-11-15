//! State vector representations.
//!
//! Keeping track of the location and velocity of an object requires more information
//! than just a position and velocity vector. Because there is no universal coordinate
//! system, positions have to be provided with respect to a reference frame.
//! There are two pieces to this, the basis of the reference frame, and the origin.
//!
//! Bringing this all together, the minimum information to know the state of an object
//! is:
//! - Frame of reference
//! - Origin
//! - Position
//! - Velocity
//! - Time
//! - ID - Some unique identifier for the object so that other objects may reference it.
//!
//! Below is the [`State`] which defines this minimum information.
//!
use serde::{Deserialize, Serialize};
use smol_str::SmolStr;
use std::fmt::{Debug, Display};

use crate::errors::{Error, KeteResult};
use crate::frames::{InertialFrame, Vector};
use crate::spice;

/// Designation for an object.
#[derive(Debug, Clone, Deserialize, Serialize, PartialEq, Hash, Eq)]
pub enum Desig {
    /// Permanent ID, an integer.
    Perm(u64),

    /// Provisional Designation
    Prov(SmolStr),

    /// Text name
    Name(SmolStr),

    /// NAIF id for the object.
    /// These are used by SPICE kernels for identification.
    Naif(i64),

    /// No id assigned.
    Empty,
}

impl From<Desig> for String {
    fn from(value: Desig) -> Self {
        match value {
            Desig::Empty => "".to_string(),
            Desig::Prov(s) => s.to_string(),
            Desig::Name(s) => s.to_string(),
            Desig::Perm(i) => i.to_string(),
            Desig::Naif(i) => i.to_string(),
        }
    }
}

impl Display for Desig {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&
            match &self {
                Desig::Empty => "".to_string(),
                Desig::Prov(s) => s.to_string(),
                Desig::Name(s) => s.to_string(),
                Desig::Perm(i) => i.to_string(),
                Desig::Naif(i) => i.to_string(),
            })
    }
}

/// Exact State of an object.
///
/// This represents the id, position, and velocity of an object with respect to a
/// coordinate frame and a center point.
///
/// This state object assumes no uncertainty in its values.
#[derive(Debug, Deserialize, Serialize, Clone)]
pub struct State<T: InertialFrame> {
    /// Designation number which corresponds to the object.
    pub desig: Desig,

    /// JD of the object's state in TDB scaled time.
    pub jd: f64,

    /// Position of the object with respect to the center_id object, units of AU.
    pub pos: Vector<T>,

    /// Velocity of the object with respect to the center_id object, units of AU/Day.
    pub vel: Vector<T>,

    /// Position and velocity are given with respect to the specified center_id.
    /// The only privileged center ID is the Solar System Barycenter 0.
    pub center_id: i64,
}

impl<T: InertialFrame> State<T> {
    /// Construct a new State object.
    #[inline(always)]
    pub fn new(desig: Desig, jd: f64, pos: Vector<T>, vel: Vector<T>, center_id: i64) -> Self {
        State {
            desig,
            jd,
            pos,
            vel,
            center_id,
        }
    }

    /// Construct a new state made of NAN pos and vel vectors but containing the
    /// remaining data. This is primarily useful as a place holder when propagation
    /// has failed and the object needs to be recorded still.
    #[inline(always)]
    pub fn new_nan(desig: Desig, jd: f64, center_id: i64) -> Self {
        Self::new(desig, jd, Vector::new_nan(), Vector::new_nan(), center_id)
    }

    /// Mutate the current state from its current [`Frame`] to the new target [`Frame`].
    ///
    /// # Arguments
    ///
    /// * `target_frame` - Target frame from the [`Frame`] enum.
    #[inline(always)]
    pub fn into_frame<B: InertialFrame>(self) -> State<B> {
        let pos = self.pos.into_frame::<B>();
        let vel = self.vel.into_frame::<B>();
        State::new(self.desig, self.jd, pos, vel, self.center_id)
    }

    /// Trade the center ID and ID values, and flip the direction of the position and
    /// velocity vectors.
    #[inline(always)]
    pub fn try_flip_center_id(&mut self) -> KeteResult<()> {
        if let Desig::Naif(mut id) = self.desig {
            (id, self.center_id) = (self.center_id, id);
            for i in 0..3 {
                self.pos[i] = -self.pos[i];
                self.vel[i] = -self.vel[i];
            }
            self.desig = Desig::Naif(id);
            return Ok(());
        }
        Err(Error::ValueError(
            "Flip center ID is only valid for NAIF ids.".into(),
        ))
    }

    /// Mutate the current state and change its center to the center defined in the
    /// provided state.
    /// For example if the current states center id is 2 (Venus), and a state is
    /// provided which represents 2 (Venus) with its center defined as 10 (Sun), then
    /// this changes the current states center to be 10 (the Sun).
    ///
    /// This will flip the center id and ID of the provided state if necessary.
    ///
    /// # Arguments
    ///
    /// * `state` - [`State`] object which defines the new center point.
    #[inline(always)]
    pub fn try_change_center(&mut self, mut state: Self) -> KeteResult<()> {
        if self.jd != state.jd {
            return Err(Error::ValueError("States don't have matching jds.".into()));
        }

        let state_id = match state.desig {
            Desig::Naif(id) => id,
            _ => {
                return Err(Error::ValueError(
                    "Changing centers only works on states with NAIF Ids.".into(),
                ))
            }
        };

        // target state does not match at all, error
        if self.center_id != state.center_id && self.center_id != state_id {
            return Err(Error::ValueError(
                "States do not reference one another at all, cannot change center.".into(),
            ));
        }

        // Flip center ID if necessary for the state
        if self.center_id == state.center_id {
            state.try_flip_center_id()?;
        }

        // Now they match, and we know the center id will change, update target frame
        // if necessary.
        let state = state.into_frame::<T>();

        // Now the state is where it is supposed to be, update as required.
        self.center_id = state.center_id;
        for i in 0..3 {
            self.pos[i] += state.pos[i];
            self.vel[i] += state.vel[i];
        }
        Ok(())
    }

    /// Attempt to update the designation from a naif id to a name.
    pub fn try_naif_id_to_name(&mut self) -> Option<()> {
        if let Desig::Naif(id) = self.desig {
            self.desig = Desig::Name(spice::try_name_from_id(id)?);
            Some(())
        } else {
            None
        }
    }
}

#[cfg(test)]
mod tests {

    use crate::frames::{Ecliptic, Equatorial};

    use super::*;

    #[test]
    fn flip_center() {
        let mut a = State::<Ecliptic>::new(
            Desig::Naif(1i64),
            0.0,
            [1.0, 0.0, 0.0].into(),
            [0.0, 1.0, 0.0].into(),
            0,
        );
        a.try_flip_center_id().unwrap();

        assert!(a.center_id == 1);
        assert!(a.pos == [-1.0, 0.0, 0.0].into());
        assert!(a.vel == [0.0, -1.0, 0.0].into());
    }

    #[test]
    fn desig_strings() {
        assert!(Desig::Empty.to_string() == "");
        assert!(Desig::Naif(100).to_string() == "100");
        assert!(Desig::Name("Foo".into()).to_string() == "Foo");
        assert!(Desig::Perm(123).to_string() == "123");
        assert!(Desig::Prov("Prov".into()).to_string() == "Prov");
    }

    #[test]
    fn naif_name_resolution() {
        let mut a = State::<Ecliptic>::new(
            Desig::Naif(1),
            0.0,
            [1.0, 0.0, 0.0].into(),
            [0.0, 1.0, 0.0].into(),
            0,
        );
        assert!(a.try_naif_id_to_name().is_some());
        assert!(a.desig == Desig::Name("mercury barycenter".into()));
        assert!(a.desig.to_string() == "mercury barycenter");
    }

    #[test]
    fn change_center() {
        let mut a = State::<Equatorial>::new(
            Desig::Naif(1),
            0.0,
            [1.0, 0.0, 0.0].into(),
            [1.0, 0.0, 0.0].into(),
            0,
        );
        let b = State::<Equatorial>::new(
            Desig::Naif(3),
            0.0,
            [0.0, 1.0, 0.0].into(),
            [0.0, 1.0, 0.0].into(),
            0,
        );
        a.try_change_center(b).unwrap();

        assert!(a.center_id == 3);
        assert!(a.pos[0] == 1.0);
        assert!(a.pos[1] != 0.0);
        assert!(a.pos[2] != 0.0);
        assert!(a.vel[0] == 1.0);

        // try cases which cause errors
        let diff_jd = State::<Equatorial>::new(
            Desig::Naif(3),
            1.0,
            [0.0, 1.0, 0.0].into(),
            [0.0, 1.0, 0.0].into(),
            0,
        );
        assert!(a.try_change_center(diff_jd).is_err());

        let not_naif_id = State::<Equatorial>::new(
            Desig::Empty,
            0.0,
            [0.0, 1.0, 0.0].into(),
            [0.0, 1.0, 0.0].into(),
            0,
        );
        assert!(a.try_change_center(not_naif_id).is_err());

        let no_matching_id = State::<Equatorial>::new(
            Desig::Naif(2),
            0.0,
            [0.0, 1.0, 0.0].into(),
            [0.0, 1.0, 0.0].into(),
            1000000000,
        );
        assert!(a.try_change_center(no_matching_id).is_err());
    }
}
