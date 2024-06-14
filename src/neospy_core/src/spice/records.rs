use std::ops::Index;
use std::io::{Read, Seek};
use std::slice::SliceIndex;

use crate::errors::NEOSpyError;

use super::binary::read_f64_vec;


/// DAF Arrays are commonly formatted with a set of some metadata, then repeated sets of records of the same format.
/// All arrays are made up of f64s.
#[derive(Debug)]
pub struct DafArray (pub Box<[f64]>);

impl DafArray {
    pub fn try_load_array<T: Read + Seek>(file: &mut T, array_start:u64, array_end: u64, little_endian: bool)-> Result<Self, NEOSpyError> {
        let _ = file.seek(std::io::SeekFrom::Start(8 * (array_start - 1)))?;
        
        let n_floats = (array_end - array_start + 1) as usize;

        let data = read_f64_vec(file, n_floats, little_endian)?;

        Ok(Self(data))
    }


    pub fn len(&self) -> usize {
        self.0.len()
    }
}


impl< Idx> Index<Idx> for DafArray where
Idx: SliceIndex<[f64], Output = f64>
{
    type Output = f64;

    fn index(&self, idx: Idx) -> &Self::Output {
        self.0.index(idx)
    }
}
