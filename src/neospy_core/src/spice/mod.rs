//! Spice replacement methods
//! Primarily includes the ability to read the contents of SPK files.
mod binary;
mod daf;
mod interpolation;
mod naif_ids;
mod pck;
mod pck_segments;
mod records;
mod spk;
mod spk_segments;

// expose the public methods in spk to the outside world.
pub use daf::*;
pub use naif_ids::try_name_from_id;
pub use pck::*;
pub use spk::*;

/// Convert seconds from J2000 into JD.
fn spice_jds_to_jd(jd_sec: f64) -> f64 {
    // 86400.0 = 60 * 60 * 24
    jd_sec / 86400.0 + 2451545.0
}

/// Convert JD to seconds from J2000.
fn jd_to_spice_jd(jd: f64) -> f64 {
    // 86400.0 = 60 * 60 * 24
    (jd -2451545.0) * 86400.0
}
