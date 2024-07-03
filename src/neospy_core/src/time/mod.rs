//! Time representation and conversions
pub mod leap_second;

/// Known Time Scales
#[derive(Debug, PartialEq, Default)]
pub enum TimeScale {
    /// TDB Scaled JD time.
    /// This is in agreement with TT up to about 1ns per century.
    #[default]
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

/// Offset from JD to MJD
/// MJD = JD + JD_TO_MJD;
const JD_TO_MJD: f64 = -2400000.5;

/// Representation of Time.
///
/// The intended accuracy is at minimum on the 50ms time scale.
/// Due to this not requiring more accuracy than this, TDB can be assumed to be TT.
///
/// Machine precision between float 64s with numbers near J2000 (IE: 2451545.0) is
/// around 23 microseconds (2.7e-10 days). So times near J2000 by necessity can only
/// be represented with about 23 microsecond accuracy. This may be dealt with by
/// splitting time into two components, integer days and a float to represent the
/// fraction of a day. However as long as accuracy below ~30ms is not required, then
/// a single f64 is sufficient.
///
/// Any conversions to a single float will by necessity result in some small accuracy
/// loss due to the nature of the representation of numbers on computers.
///
#[derive(Debug)]
pub struct Time {
    /// Float Date
    julian_date: f64,

    /// Time scaling.
    scale: TimeScale,
}

impl Time {
    /// Create Time from an Modified Julian Date (MJD).
    pub fn from_mjd(&self, mjd: f64, scale: TimeScale) -> Self {
        Self {
            julian_date: mjd - JD_TO_MJD,
            scale,
        }
    }
    /// Convert to an MJD float.
    /// This may result in some precision loss.
    pub fn mjd(&self) -> f64 {
        self.julian_date + JD_TO_MJD
    }

    /// Rescale time to the desired time scaling.
    pub fn rescale(&mut self, scale: &TimeScale) {
        if &self.scale == scale {
            return;
        }
        match &self.scale {
            TimeScale::TDB => (),
            TimeScale::TAI => {
                self.julian_date += TDB_TO_TAI;
            }
            TimeScale::UT1 => todo!(),
        };

        match scale {
            TimeScale::TDB => (),
            TimeScale::TAI => self.julian_date -= TDB_TO_TAI,
            TimeScale::UT1 => todo!(),
        };
    }

    /// Return the Gregorian year, month, day, and fraction of a day.
    ///
    /// Algorithm from:
    /// "A Machine Algorithm for Processing Calendar Dates"
    /// https://doi.org/10.1145/364096.364097
    ///
    pub fn year_month_day(&self) -> (i64, u8, u8, f64) {
        let offset = self.julian_date + 0.5;
        let frac_day = offset.rem_euclid(1.0);

        let mut l = offset.div_euclid(1.0) as i64 + 68569;

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

    /// Create Time from the date in the Gregorian calendar.
    ///
    /// Algorithm from:
    /// "A Machine Algorithm for Processing Calendar Dates"
    /// https://doi.org/10.1145/364096.364097
    ///
    pub fn from_year_month_day(
        year: i64,
        month: u8,
        day: u8,
        frac_day: f64,
        scale: TimeScale,
    ) -> Self {
        let day = day as i64 + frac_day.div_euclid(1.0) as i64;
        let frac_day = frac_day.rem_euclid(1.0);
        let month = month as i64;

        let tmp = (month - 14) / 12;
        let days = day - 32075 + 1461 * (year + 4800 + tmp) / 4 + 367 * (month - 2 - tmp * 12) / 12
            - 3 * ((year + 4900 + tmp) / 100) / 4;

        Self {
            scale,
            julian_date: days as f64 + frac_day,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_leap_second() {
        let t = Time {
            scale: TimeScale::TAI,
            julian_date: 2451545.,
        };
        assert!(t.year_month_day() == (2000, 1, 1, 0.5));

        let t2 = Time::from_year_month_day(2000, 1, 1, 0.5, TimeScale::TAI);
        assert!(t2.julian_date == 2451545.5);

        let t2 = Time::from_year_month_day(2000, 1, 2, -0.5, TimeScale::TAI);
        assert!(t2.julian_date == 2451545.5);

        let t = Time {
            scale: TimeScale::TAI,
            julian_date: 2000000.,
        };
        assert!(t.year_month_day() == (763, 9, 18, 0.5));

        let t2 = Time::from_year_month_day(763, 9, 18, 0., TimeScale::TAI);
        assert!(t2.julian_date == 2000000.);
    }
}
