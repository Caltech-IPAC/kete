use std::ops::Index;
// use std::io::{Read, Seek};

use crate::errors::NEOSpyError;

// use super::binary::read_f64_vec;

#[derive(Debug)]
pub struct DafRecords {
    data: Vec<f64>,
    n_records: usize,
    record_len: usize,
}

impl DafRecords {
    pub fn with_capacity(n_records: usize, record_len: usize) -> DafRecords {
        DafRecords {
            data: Vec::with_capacity(n_records * record_len),
            n_records: 0,
            record_len,
        }
    }

    pub fn try_push(&mut self, record: &[f64]) -> Result<(), NEOSpyError> {
        if record.len() != self.record_len {
            return Err(NEOSpyError::DAFLimits(
                "Record length does not match expected.".into(),
            ));
        }
        self.n_records += 1;
        self.data.extend(record);
        Ok(())
    }

    // pub fn try_load_record<T: Read + Seek>(&mut self, file: &mut T, idx:u64, little_endian: bool)-> Result<(), NEOSpyError> {
    //     let _ = file.seek(std::io::SeekFrom::Start(1024 * (idx - 1)))?;
    //     let record = read_f64_vec(file, self.record_len, little_endian)?;
    //     self.n_records += 1;
    //     self.data.extend(record.iter());
    //     Ok(())
    // }

    pub unsafe fn get_unchecked(&self, idx: usize) -> &[f64] {
        self.data
            .get_unchecked(idx * self.record_len..(idx + 1) * self.record_len)
    }

    pub fn len(&self) -> usize {
        self.n_records
    }
}

impl Index<usize> for DafRecords {
    type Output = [f64];

    fn index(&self, idx: usize) -> &Self::Output {
        &self.data[idx * self.record_len..(idx + 1) * self.record_len]
    }
}
