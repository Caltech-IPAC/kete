//! Spice specific coordinate frames
//!

/// Inertial frames of reference saved in SPICE
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum SpiceFrames {
    /// Equatorial J2000, effectively the same as ICRF
    J2000,

    /// Ecliptic J2000, same as equatorial above with a fixed X rotation.
    ECLIPJ2000,
}

impl From<i32> for SpiceFrames {
    fn from(value: i32) -> Self {
        match value {
            1 => SpiceFrames::J2000,       // Equatorial
            17 => SpiceFrames::ECLIPJ2000, // Ecliptic
            i => panic!("{} is an unsupported frame", i),
        }
    }
}
