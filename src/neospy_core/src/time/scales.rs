//! Available Time Scales
//!
//! TT Varies from TDB by up to about 1.7 ms per year in a period manner.
//! This correction is a complicated relationship due to relativistic motion of
//! the observer on Earth vs the rest of the solar system. For performance reasons,
//! the conversion from TDB to TT is unnecessary, as they are never more than about
//! 2 ms apart over a century.
use super::leap_second::tai_to_utc_offset;

/// Offset from TT to TAI.
/// Definitional offset from TT to TAI, however we treat TT==TDB.
///
/// TT = TAI + TT_TO_TAI
/// TAI = TT - TT_TO_TAI
pub(crate) const TT_TO_TAI: f64 = 32.184 / 86400.0;

/// Offset from JD to MJD
/// MJD = JD + JD_TO_MJD;
pub const JD_TO_MJD: f64 = -2400000.5;

/// Time Scaling support, all time scales must implement this.
pub trait TimeScale {
    /// Convert julian date to TDB scaled julian date.
    fn to_tdb(jd: f64) -> f64;

    /// Convert from TDB scaled julian date.
    fn from_tdb(jd: f64) -> f64;
}

/// TDB Scaled JD time.
///
/// This is in good agreement with "Ephemeris Time" - Which is commonly referred
/// to as the time units used by the JPL Ephemeris, and is used in NEOSpy as the
/// base time scaling.
///
/// This is essentially the rate of time from an observer not on the surface of
/// the Earth (and as a result doesn't feel the relativistic dilation effects of
/// Earth).
///
/// This is in agreement with TT up to about 1.7ms per century.
#[derive(Debug, Clone, Copy)]
pub struct TDB;

impl TimeScale for TDB {
    fn to_tdb(t: f64) -> f64 {
        t
    }
    fn from_tdb(t: f64) -> f64 {
        t
    }
}

/// UTC Scaled JD time.
#[derive(Debug, Clone, Copy)]
pub struct UTC;

impl TimeScale for UTC {
    fn to_tdb(mut jd: f64) -> f64 {
        // move to tai first
        // Guess the tai time that this UTC time came from
        let offset = tai_to_utc_offset(&(jd + JD_TO_MJD));

        // use that guess to update the TAI guess to fix the leap second offset
        let offset = tai_to_utc_offset(&(jd + JD_TO_MJD + offset));
        jd += offset;

        // then tai to tdb
        jd += TT_TO_TAI;
        jd
    }
    fn from_tdb(mut jd: f64) -> f64 {
        // convert from TDB to TAI
        jd -= TT_TO_TAI;

        // Time is now TAI
        // calculate leap seconds for that time to convert from TAI to UTC
        let offset = tai_to_utc_offset(&(jd + JD_TO_MJD));
        jd -= offset;
        jd
    }
}

/// TAI Time
/// This is the international standard for the measurement of time.
/// Atomic clocks around the world keep track of this time, which then gets
/// converted to UTC time which is commonly used. UTC Time is an approximation
/// of the current time with corrections for the rotation rate of the Earth.
///
/// As the Earth actually varies in its rotation speed,
#[derive(Debug, Clone, Copy)]
pub struct TAI;

impl TimeScale for TAI {
    fn from_tdb(jd: f64) -> f64 {
        jd - TT_TO_TAI
    }
    fn to_tdb(jd: f64) -> f64 {
        jd + TT_TO_TAI
    }
}
