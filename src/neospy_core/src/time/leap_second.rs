//! Leap Second information
use std::str::FromStr;

use lazy_static::lazy_static;
use serde::Deserialize;

use crate::prelude::NEOSpyError;

/// Leap Second Information
#[derive(Debug, Deserialize)]
pub struct LeapSecond {
    ///  MJD
    pub mjd: f64,

    /// Date in d/m/year
    pub date: (u8, u8, u32),

    /// offset
    pub tai_m_utc: i16,
}

impl FromStr for LeapSecond {
    type Err = NEOSpyError;

    /// Load an LeapSecond from a single string.
    fn from_str(row: &str) -> Result<Self, Self::Err> {
        let row: Vec<_> = row.split_whitespace().collect();
        if row.len() != 5 {
            return Err(NEOSpyError::IOError(
                "Leap Second file incorrectly formatted.".into(),
            ));
        }
        Ok(LeapSecond {
            mjd: row[0].parse()?,
            date: (row[1].parse()?, row[2].parse()?, row[3].parse()?),
            tai_m_utc: row[4].parse()?,
        })
    }
}

const PRELOAD_LEAPSECONDS: &[u8] = include_bytes!("../../data/leap_second.dat");

lazy_static! {
    /// Leap second definitions
    pub static ref LEAP_SECONDS: Vec<LeapSecond> = {
        let mut codes = Vec::new();
        let text = std::str::from_utf8(PRELOAD_LEAPSECONDS).unwrap().split('\n');
        for row in text.filter(|x| !x.starts_with('#') & (!x.trim().is_empty())) {
            let code = LeapSecond::from_str(row).unwrap();
            codes.push(code);
        }
        codes
    };
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_leap_second() {
        let t = &LEAP_SECONDS[0];
        assert!(t.date == (1, 1, 1972));
        assert!(t.tai_m_utc == 10);
        assert!(t.mjd == 41317.0);

        let t = &LEAP_SECONDS[27];
        assert!(t.date == (1, 1, 2017));
        assert!(t.tai_m_utc == 37);
        assert!(t.mjd == 57754.0);
    }
}
