//! Observatory codes used by the MPC
use lazy_static::lazy_static;
use serde::Deserialize;

use crate::prelude::{Error, KeteResult};
use std::str;
use std::str::FromStr;

/// Observatory information
#[derive(Debug, Deserialize)]
pub struct ObsCode {
    /// observatory code
    pub code: String,

    /// longitude in degrees
    pub lon: f64,

    /// latitude in degrees
    pub lat: f64,

    /// altitude in meters
    pub altitude: f64,

    /// name of the observatory
    pub name: String,
}

impl FromStr for ObsCode {
    type Err = Error;

    /// Load an ObsCode from a single string.
    fn from_str(row: &str) -> KeteResult<Self> {
        let code = row[3..6].to_string();
        let lon = f64::from_str(row[8..19].trim()).unwrap();
        let lat = f64::from_str(row[18..29].trim()).unwrap();
        let altitude = f64::from_str(row[29..39].trim()).unwrap() / 1000.0;
        let name = row[79..].trim().to_string();
        Ok(ObsCode {
            code,
            lon,
            lat,
            altitude,
            name,
        })
    }
}

const PRELOAD_OBS: &[u8] = include_bytes!("../../data/mpc_obs.tsv");

lazy_static! {
    /// Observatory Codes
    pub static ref OBS_CODES: Vec<ObsCode> = {
        let mut codes = Vec::new();
        let text = str::from_utf8(PRELOAD_OBS).unwrap().split('\n');
        for row in text.skip(1) {
            let code: ObsCode = ObsCode::from_str(row).unwrap();
            codes.push(code);
        }
        codes
    };
}

#[cfg(test)]
mod tests {

    use super::*;

    #[test]
    fn obs_codes() {
        let codes = &OBS_CODES;
        assert!(!codes.is_empty());
    }
}
