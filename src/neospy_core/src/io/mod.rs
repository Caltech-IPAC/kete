//! File IO related tools

pub mod obs_codes;
pub mod serde_const_arr;

use crate::prelude::NEOSpyError;
use bincode::serde::{decode_from_std_read, encode_into_std_write};
use serde::{Deserialize, Serialize};
use std::fs::File;
use std::io::{BufReader, BufWriter};

/// Support for automatic derivation of Save/Load
pub trait FileIO: Serialize
where
    for<'de> Self: Deserialize<'de>,
{
    /// Save into a file.
    fn save(&self, filename: String) -> Result<usize, NEOSpyError> {
        let mut f = BufWriter::new(File::create(filename)?);
        encode_into_std_write(self, &mut f, bincode::config::legacy())
            .map_err(|_| NEOSpyError::IOError("Failed to write to file".into()))
    }

    /// Load from a file.
    fn load(filename: String) -> Result<Self, NEOSpyError> {
        let mut f = BufReader::new(File::open(filename)?);
        decode_from_std_read(&mut f, bincode::config::legacy())
            .map_err(|_| NEOSpyError::IOError("Failed to read from file".into()))
    }

    /// Save a vector of this object.
    fn save_vec(vec: &[Self], filename: String) -> Result<(), NEOSpyError> {
        let mut f = BufWriter::new(File::create(filename)?);

        let _ = encode_into_std_write(vec, &mut f, bincode::config::legacy())
            .map_err(|_| NEOSpyError::IOError("Failed to write to file".into()))?;
        Ok(())
    }

    /// load a vector of this object.
    fn load_vec(filename: String) -> Result<Vec<Self>, NEOSpyError> {
        let mut f = BufReader::new(File::open(filename)?);

        let res: Vec<Self> = decode_from_std_read(&mut f, bincode::config::legacy())
            .map_err(|_| NEOSpyError::IOError("Failed to load from file".into()))?;
        Ok(res)
    }
}
