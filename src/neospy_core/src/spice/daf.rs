/// Support for arbitrary DAF files
/// DAF is a superset which includes SPK files.
///
/// DAF files are laid out in 1024 Byte "Records"
/// - The first record is header information about the contents of the file.
/// - The following N records are text comments.
/// - Immediately following the comments there is a Summary Record.
///
/// These summary records contain the location information for all the contents
/// of the DAF file.
///
use super::{binary::{
    bytes_to_f64, bytes_to_f64_vec, bytes_to_i32, bytes_to_i32_vec, bytes_to_str, read_bytes_exact,
    read_str,
}, records::DafRecords};
use crate::errors::NEOSpyError;
use std::io::{Read, Seek};

/// DAF Files can contain multiple different types of data.
/// This list contains the supported formats.
#[derive(Debug, PartialEq)]
pub enum DAFType {
    /// SPK files are planetary and satellite ephemeris data.
    Spk,

    /// PCK Files are planetary and satellite orientation data.
    Pck,
}

/// DAF files header information.
/// This contains
#[derive(Debug)]
pub struct DafHeader {
    /// Magic number within the DAF file corresponds to this DAF type.
    pub daf_type: DAFType,

    /// Number of f64 in each array.
    pub n_doubles: i32,

    /// Number of i32s in each array.
    pub n_ints: i32,

    /// Number of chars in the descriptor string of each array.
    pub n_chars: i32,

    /// Is the file little endian.
    pub little_endian: bool,

    /// Internal Descriptor.
    pub internal_desc: String,

    /// Index of initial summary record.
    /// Note that this is 1 indexed and corresponds to record index
    /// not file byte index.
    pub init_summary_record_index: i32,

    /// Index of final summary record.
    /// Note that this is 1 indexed and corresponds to record index
    /// not file byte index.
    pub final_summary_record_index: i32,

    /// First free address of the file.
    /// Index of initial summary record
    /// Note that this is 1 indexed.
    pub first_free: i32,

    /// FTP Validation string
    pub ftp_validation_str: String,

    /// The comment records.
    /// Each record is trimmed to 1000 chars, as that is what SPKs use internally.
    pub comments: String,
}

type DafSummary = (Box<[f64]>, Box<[i32]>);

impl DafHeader {
    /// Try to load a single record from the DAF.
    pub fn try_load_record<T: Read + Seek>(
        file: &mut T,
        idx: u64,
    ) -> Result<Box<[u8]>, NEOSpyError> {
        let _ = file.seek(std::io::SeekFrom::Start(1024 * (idx - 1)))?;
        read_bytes_exact(file, 1024)
    }

    /// Load the header information from the DAF file.
    pub fn try_load_header<T: Read + Seek>(mut file: T) -> Result<Self, NEOSpyError> {
        let bytes = Self::try_load_record(&mut file, 1)?;
        let daf_type = Self::parse_magic(&bytes_to_str(&bytes[0..8]))?;

        let little_endian = match bytes_to_str(&bytes[88..96]).to_lowercase().as_str() {
            "ltl-ieee" => true,
            "big-ieee" => false,
            _ => {
                return Err(NEOSpyError::IOError(
                    "Expected little or big endian in DAF file, found neither".into(),
                ))
            }
        };

        let n_doubles = bytes_to_i32(&bytes[8..12], little_endian)?;
        let n_ints = bytes_to_i32(&bytes[12..16], little_endian)?;
        let n_chars= 8 * (n_doubles + (n_ints + 1) / 2);

        // record index of the first summary record in the file
        // records are 1024 long, and 1 indexed because fortran.
        let init_summary_record_index = bytes_to_i32(&bytes[76..80], little_endian)?;

        // the following values are not used, so are not stored.
        let internal_desc = bytes_to_str(&bytes[16..76]);
        let final_summary_record_index = bytes_to_i32(&bytes[80..84], little_endian)?;
        let first_free = bytes_to_i32(&bytes[84..88], little_endian)?;

        let ftp_validation_str = bytes_to_str(&bytes[966..966 + 28]);

        // after the header, there are comments until the first record index.
        // so read the next (init_summary_record_index-2) records:
        // -1 for fortran indexing
        // -1 for having already read a single record
        let mut comments: Vec<String> = Vec::with_capacity(init_summary_record_index as usize - 2);
        for _ in 0..(init_summary_record_index - 2) {
            // TODO: Check if the 1000 character limit is what other formats use.
            // 1k is used by SPK for sure.
            comments.push(read_str(&mut file, 1024)?.split_at(1000).0.into());
        }

        Ok(DafHeader {
            daf_type,
            n_doubles,
            n_ints,
            n_chars,
            little_endian,
            internal_desc,
            init_summary_record_index,
            final_summary_record_index,
            first_free,
            ftp_validation_str,
            comments: comments.join(""),
        })
    }

    /// Load all summary records from the DAF file.
    /// These are tuples containing a series of f64s and i32s.
    /// The meaning of these values depends on the particular implementation of the DAF,
    /// IE: SPK files have 2 floats and 6 ints.
    pub fn try_load_summaries<T: Read + Seek>(
        &self,
        file: &mut T,
    ) -> Result<Vec<DafSummary>, NEOSpyError> {
        let summary_size = self.n_doubles + (self.n_ints + 1) / 2;

        let mut next_idx = self.init_summary_record_index;
        let mut summaries: Vec<DafSummary> = Vec::new();
        loop {
            if next_idx == 0 {
                break;
            }
            let bytes = DafHeader::try_load_record(file, next_idx as u64)?;

            next_idx = bytes_to_f64(&bytes[0..8], self.little_endian)? as i32;
            // let prev_idx = bytes_to_f64(&bytes[8..16], daf.little_endian)? as i32;
            let n_summaries = bytes_to_f64(&bytes[16..24], self.little_endian)? as i32;

            for idy in 0..n_summaries {
                let sum_start = (3 * 8 + idy * summary_size * 8) as usize;
                let floats = bytes_to_f64_vec(
                    &bytes[sum_start..sum_start + 8 * self.n_doubles as usize],
                    self.little_endian,
                )?;
                let ints = bytes_to_i32_vec(
                    &bytes[sum_start + 8 * self.n_doubles as usize
                        ..sum_start + (8 * self.n_doubles + 4 * self.n_ints) as usize],
                    self.little_endian,
                )?;
                summaries.push((floats, ints));
            }
        }
        Ok(summaries)
    }

    /// Try to parse the magic string of the DAF file into the correct DAFType.
    pub fn parse_magic(magic: &str) -> Result<DAFType, NEOSpyError> {
        match &magic.to_uppercase()[4..7] {
            "SPK" => Ok(DAFType::Spk),
            "PCK" => Ok(DAFType::Pck),
            _ => Err(NEOSpyError::IOError(format!(
                "File type currently unsupported by NEOSpy. Magic Number {:?}",
                magic
            ))),
        }
    }
}


/// Contents of an DAF file.
#[derive(Debug)]
#[allow(dead_code)]
pub struct DafFile{
    header: Box<DafHeader>,

    summaries: Vec<DafSummary>,

    records: Box<DafRecords>,
}

impl  DafFile {
    
    /// Load all summary records from the DAF file.
    /// These are tuples containing a series of f64s and i32s.
    /// The meaning of these values depends on the particular implementation of the DAF,
    /// IE: SPK files have 2 floats and 6 ints.
    pub fn try_load_summaries<T: Read + Seek>(
        &self,
        file: &mut T,
    ) -> Result<Vec<DafSummary>, NEOSpyError> {
        let summary_size = self.header.n_doubles + (self.header.n_ints + 1) / 2;

        let mut next_idx = self.header.init_summary_record_index;
        let mut summaries: Vec<DafSummary> = Vec::new();
        loop {
            if next_idx == 0 {
                break;
            }
            let bytes = DafHeader::try_load_record(file, next_idx as u64)?;

            next_idx = bytes_to_f64(&bytes[0..8], self.header.little_endian)? as i32;
            // let prev_idx = bytes_to_f64(&bytes[8..16], daf.little_endian)? as i32;
            let n_summaries = bytes_to_f64(&bytes[16..24], self.header.little_endian)? as i32;

            for idy in 0..n_summaries {
                let sum_start = (3 * 8 + idy * summary_size * 8) as usize;
                let floats = bytes_to_f64_vec(
                    &bytes[sum_start..sum_start + 8 * self.header.n_doubles as usize],
                    self.header.little_endian,
                )?;
                let ints = bytes_to_i32_vec(
                    &bytes[sum_start + 8 * self.header.n_doubles as usize
                        ..sum_start + (8 * self.header.n_doubles + 4 * self.header.n_ints) as usize],
                        self.header.little_endian,
                )?;
                summaries.push((floats, ints));
            }
        }
        Ok(summaries)
    }
}