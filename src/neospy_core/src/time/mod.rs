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
    days: i64,

    /// Fraction of the day.
    frac_day: f64,
}

impl Time {
    /// Rescale time to the desired time scaling.
    pub fn rescale(&mut self, scale: &Scale) {
        if &self.scale == scale {
            return;
        }
        match &self.scale {
            Scale::TDB => {
                self.frac_day -= TDB_TO_TAI;
            }
            Scale::TAI => (),
            Scale::UT1 => todo!(),
        };

        match scale {
            Scale::TDB => self.frac_day += TDB_TO_TAI,
            Scale::TAI => (),
            Scale::UT1 => todo!(),
        };
        self.normalize()
    }

    /// Return the Gregorian year, month, day, and fraction of a day.
    pub fn year_month_day(&self) -> (i64, u8, u8, f64) {
        let mut l = self.days + 68569;

        let mut frac_day = 0.5 + self.frac_day;
        l += frac_day.div_euclid(1.0) as i64;
        frac_day = frac_day.rem_euclid(1.0);

        let n = (4 * l) / 146097;
        l -= (146097 * n + 3) / 4;
        let i = (4000 * (l + 1)) / 1461001;
        l -= (1461 * i) / 4 - 31;
        let k = (80 * l) / 2447;
        let day = l - (2447 * k) / 80;
        l = k / 11;

        let month = k + 2 - 12 * l;
        let year = 100 * (n - 49) + i + l;
        (year, month as u8, day as u8, frac_day)
    }

    /// Shift value of the day frac so that it is between 0-1, any value outside of
    /// that range is added or subtracted to the days count as appropriate.
    #[inline(always)]
    fn normalize(&mut self) {
        self.days += self.frac_day.div_euclid(1.0) as i64;
        self.frac_day = self.frac_day.rem_euclid(1.0);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_leap_second() {
        let t = Time {
            scale: Scale::TAI,
            days: 2451545,
            frac_day: 0.0,
        };
        assert!(t.year_month_day() == (2000, 1, 1, 0.5));

        let t = Time {
            scale: Scale::TAI,
            days: 2000000,
            frac_day: 0.0,
        };
        assert!(t.year_month_day() == (763, 9, 18, 0.5));
    }
}
