//! List of NAIF ID values.
//! This list is not comprehensive, but is more complete than the C-SPICE
//! implementation.
use lazy_static::lazy_static;
use serde::Deserialize;
use smol_str::SmolStr;

use crate::prelude::{Error, KeteResult};
use std::str;
use std::str::FromStr;

/// NAIF ID information
#[derive(Debug, Deserialize)]
pub struct NaifId {
    /// NAIF id
    pub id: i32,

    /// name of the object
    pub name: SmolStr,
}

impl FromStr for NaifId {
    type Err = Error;

    /// Load an NaifId from a single string.
    fn from_str(row: &str) -> KeteResult<Self> {
        let id = i32::from_str(row[0..10].trim()).unwrap();
        let name = row[11..].trim().to_string().into();
        Ok(NaifId { id, name })
    }
}

const PRELOAD_IDS: &[u8] = include_bytes!("../../data/naif_ids.csv");

lazy_static! {
    /// NAIF Ids
    pub static ref NAIF_IDS: Box<[NaifId]> = {
        let mut ids = Vec::new();
        let text = str::from_utf8(PRELOAD_IDS).unwrap().split('\n');
        for row in text.skip(1) {
            ids.push(NaifId::from_str(row).unwrap());
        }
        ids.into()
    };
}

/// Return the string name of the desired ID if possible.
pub fn try_name_from_id(id: i64) -> Option<SmolStr> {
    let id = id as i32;
    for naif_id in NAIF_IDS.iter() {
        if naif_id.id == id {
            return Some(naif_id.name.clone());
        }
    }
    None
}
