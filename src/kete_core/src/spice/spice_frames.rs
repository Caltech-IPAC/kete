//! Spice specific coordinate frames
//!

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum SpiceFrames {
    J2000,
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
