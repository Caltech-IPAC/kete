//! Time representation and conversions
pub mod leap_second;

/// Known Time Scales
///
///          TAI             <-  physically realized
///           :
///         offset           <-  observed (nominally +32.184s)
///           :
///        TT / TDB             <-  terrestrial time / TDB

#[derive(Debug, PartialEq)]
pub enum Scale {
    /// TDB Scaled JD time.
    /// This is in agreement with TT up to about 1ns per century.
    TDB,

    /// TAI Time
    /// This is the international standard for the measurement of time.
    TAI,

    /// UT1 Scaled JD time.
    UT1,
}

/// Offset from TDB to TAI.
/// Technically the definitional offset from TT to TAI, however we treat TT==TDB.
const TDB_TO_TAI: f64 = 32.184 / 86400.0;

/// Representation of Time.
///
/// The intended precision is at minimum on the millisecond accuracy scale.
/// Due to this not requiring more precision than this, TDB can be assumed to be TT.
///
/// Underlying time representation is two numbers, one is an Integer Day,
/// the second is a float and is a fraction of a day. Splitting the time into
/// two components like this significantly increases numerical accuracy when the
/// integer portion is very large. This is because at a fundamental level computers
/// represent floats as discrete steps, and when the values are large the step
/// between neighboring numbers is larger.
#[derive(Debug)]
pub struct Time {
    /// Time scaling.
    scale: Scale,

    /// Large component of the date.
    _days: u64,

    /// Sub day fraction of the day.
    day: f64,
}

impl Time {
    /// Rescale time to the desired time scaling.
    pub fn rescale(&mut self, scale: &Scale) {
        if &self.scale == scale {
            return;
        }
        match &self.scale {
            Scale::TDB => {
                self.day -= TDB_TO_TAI;
            }
            Scale::TAI => (),
            Scale::UT1 => todo!(),
        };

        match scale {
            Scale::TDB => self.day += TDB_TO_TAI,
            Scale::TAI => (),
            Scale::UT1 => todo!(),
        };
    }
}
