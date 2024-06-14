use super::binary::read_f64_vec;
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
use super::binary::{
    bytes_to_f64, bytes_to_f64_vec, bytes_to_i32, bytes_to_i32_vec, bytes_to_string,
    read_bytes_exact, read_str,
};
use super::pck_segments::PckSegment;
use super::spk_segments::SpkSegment;

use crate::errors::NEOSpyError;
use std::io::Cursor;
use std::{
    io::{Read, Seek},
    ops::Index,
    slice::SliceIndex,
};

/// DAF Files can contain multiple different types of data.
/// This list contains the supported formats.
#[derive(Debug)]
pub enum DAFType {
    /// SPK files are planetary and satellite ephemeris data.
    Spk,

    /// PCK Files are planetary and satellite orientation data.
    Pck,

    /// An unrecognized DAF type.
    Unrecognized(String),
}

impl From<&str> for DAFType {
    fn from(magic: &str) -> Self {
        match &magic.to_uppercase()[4..7] {
            "SPK" => DAFType::Spk,
            "PCK" => DAFType::Pck,
            other => DAFType::Unrecognized(other.into()),
        }
    }
}
/// DAF Files can contain multiple different types of data.
/// This list contains the supported formats.
#[derive(Debug)]
pub enum DafSegments {
    /// SPK files are planetary and satellite ephemeris data.
    Spk(SpkSegment),

    /// PCK Files are planetary and satellite orientation data.
    Pck(PckSegment),
}

impl DafSegments {
    /// Return the contained SPK segment.
    /// Panic if not SPK.
    pub fn spk(self) -> SpkSegment {
        if let Self::Spk(seg) = self {
            seg
        } else {
            panic!("Not an SPK segment.")
        }
    }

    /// Return the contained PCK segment.
    /// Panic if not PCK.
    pub fn pck(self) -> PckSegment {
        if let Self::Pck(seg) = self {
            seg
        } else {
            panic!("Not an PCK segment.")
        }
    }
}

/// DAF files header information.
/// This contains
#[derive(Debug)]
pub struct DafFile {
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

    /// DAF segments
    pub segments: Vec<DafSegments>,
}

type DafSummary = (Box<[f64]>, Box<[i32]>);

impl DafFile {
    /// Try to load a single record from the DAF.
    pub fn try_load_record<T: Read + Seek>(
        file: &mut T,
        idx: u64,
    ) -> Result<Box<[u8]>, NEOSpyError> {
        let _ = file.seek(std::io::SeekFrom::Start(1024 * (idx - 1)))?;
        read_bytes_exact(file, 1024)
    }

    /// Load the contents of a DAF file.
    pub fn from_buffer<T: Read + Seek>(mut buffer: T) -> Result<Self, NEOSpyError> {
        let bytes = Self::try_load_record(&mut buffer, 1)?;
        let daf_type: DAFType = bytes_to_string(&bytes[0..8]).as_str().into();

        let little_endian = match bytes_to_string(&bytes[88..96]).to_lowercase().as_str() {
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
        let n_chars = 8 * (n_doubles + (n_ints + 1) / 2);

        // record index of the first summary record in the file
        // records are 1024 long, and 1 indexed because fortran.
        let init_summary_record_index = bytes_to_i32(&bytes[76..80], little_endian)?;

        // the following values are not used, so are not stored.
        let internal_desc = bytes_to_string(&bytes[16..76]);
        let final_summary_record_index = bytes_to_i32(&bytes[80..84], little_endian)?;
        let first_free = bytes_to_i32(&bytes[84..88], little_endian)?;

        let ftp_validation_str = bytes_to_string(&bytes[966..966 + 28]);

        // after the header, there are comments until the first record index.
        // so read the next (init_summary_record_index-2) records:
        // -1 for fortran indexing
        // -1 for having already read a single record
        let mut comments: Vec<String> = Vec::with_capacity(init_summary_record_index as usize - 2);
        for _ in 0..(init_summary_record_index - 2) {
            // TODO: Check if the 1000 character limit is what other formats use.
            // 1k is used by SPK for sure.
            comments.push(read_str(&mut buffer, 1024)?.chars().take(1000).collect());
        }

        let mut daf = DafFile {
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
            segments: Vec::new(),
        };

        daf.try_load_segments(&mut buffer)?;
        Ok(daf)
    }

    /// Load DAF file from the specified filename.
    pub fn from_file(filename: &str) -> Result<Self, NEOSpyError> {
        let mut file = std::fs::File::open(filename)?;
        let mut buffer = Vec::new();
        let _ = file.read_to_end(&mut buffer)?;
        let mut buffer = Cursor::new(&buffer);
        Self::from_buffer(&mut buffer)
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
            let bytes = DafFile::try_load_record(file, next_idx as u64)?;

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

    /// Attempt to load all arrays
    pub fn try_load_segments<T: Read + Seek>(&mut self, file: &mut T) -> Result<(), NEOSpyError> {
        let summaries = self.try_load_summaries(file)?;

        for (floats, ints) in summaries {
            match self.daf_type {
                DAFType::Spk => {
                    let seg = DafSegments::Spk(SpkSegment::from_summary(file, &floats, &ints)?);
                    self.segments.push(seg);
                }
                DAFType::Pck => {
                    let seg = DafSegments::Pck(PckSegment::from_summary(file, &floats, &ints)?);
                    self.segments.push(seg);
                }
                _ => panic!(),
            }
        }

        Ok(())
    }
}

/// DAF Arrays are f64 arrays of structured data.
///
/// Contents of the structure depends on specific file formats, however they are all
/// made up of floats.
#[derive(Debug)]
pub struct DafArray(pub Box<[f64]>);

impl DafArray {
    /// Load array from file
    pub fn try_load_array<T: Read + Seek>(
        file: &mut T,
        array_start: u64,
        array_end: u64,
        little_endian: bool,
    ) -> Result<Self, NEOSpyError> {
        let _ = file.seek(std::io::SeekFrom::Start(8 * (array_start - 1)))?;

        let n_floats = (array_end - array_start + 1) as usize;

        let data = read_f64_vec(file, n_floats, little_endian)?;

        Ok(Self(data))
    }

    /// Total length of the array.
    pub fn len(&self) -> usize {
        self.0.len()
    }

    /// Test if array is empty.
    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }
}

impl<Idx> Index<Idx> for DafArray
where
    Idx: SliceIndex<[f64], Output = f64>,
{
    type Output = f64;

    fn index(&self, idx: Idx) -> &Self::Output {
        self.0.index(idx)
    }
}
