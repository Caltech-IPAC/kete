//! Loading and reading of states from JPL PCK kernel files.
//!
//! PCKs are intended to be loaded into a singleton which is accessible via the
//! [`get_pck_singleton`] function defined below. This singleton is wrapped in a RwLock,
//! meaning before its use it must by unwrapped. A vast majority of intended use cases
//! will only be the read case.
//!
use super::daf::{DAFType, DafFile};
use super::pck_segments::PckSegment;
use crate::errors::{Error, KeteResult};
use crate::frames::EclipticNonInertial;
use crossbeam::sync::ShardedLock;
use lazy_static::lazy_static;

use std::io::Cursor;

const PRELOAD_PCK: &[&[u8]] = &[
    include_bytes!("../../data/earth_000101_240215_231123.bpc"),
    include_bytes!("../../data/earth_200101_990825_predict.bpc"),
];

/// A collection of segments.
#[derive(Debug)]
pub struct PckCollection {
    /// Collection of PCK file information
    pub segments: Vec<PckSegment>,
}

/// Define the PCK singleton structure.
pub type PckSingleton = ShardedLock<PckCollection>;

impl PckCollection {
    /// Given an PCK filename, load all the segments present inside of it.
    /// These segments are added to the PCK singleton in memory.
    pub fn load_file(&mut self, filename: &str) -> KeteResult<()> {
        let file = DafFile::from_file(filename)?;
        if !matches!(file.daf_type, DAFType::Pck) {
            Err(Error::IOError(format!(
                "File {:?} is not a PCK formatted file.",
                filename
            )))?;
        }
        self.segments
            .extend(file.segments.into_iter().map(|x| x.pck()));
        Ok(())
    }

    /// Get the raw orientation from the loaded PCK files.
    /// This orientation will have the frame of what was originally present in the file.
    pub fn try_get_orientation(&self, id: isize, jd: f64) -> KeteResult<EclipticNonInertial> {
        for segment in self.segments.iter() {
            if id == segment.center_id && segment.contains(jd) {
                return segment.try_get_orientation(jd);
            }
        }

        Err(Error::DAFLimits(format!(
            "Object ({}) does not have an PCK record for the target JD.",
            id
        )))
    }

    /// Delete all segments in the PCK singleton, equivalent to unloading all files.
    pub fn reset(&mut self) {
        let files = PckCollection {
            segments: Vec::new(),
        };

        *self = files;

        for preload in PRELOAD_PCK {
            let mut buffer = Cursor::new(preload);
            let file = DafFile::from_buffer(&mut buffer).unwrap();
            self.segments
                .extend(file.segments.into_iter().map(|x| x.pck()))
        }
    }
}

lazy_static! {
    /// PCK singleton.
    /// This is a RwLock protected PCKCollection, and must be `.try_read().unwrapped()` for any
    /// read-only cases.
    pub static ref LOADED_PCK: PckSingleton = {
        let mut files = PckCollection {
            segments: Vec::new(),
        };
        files.reset();
        ShardedLock::new(files)
    };
}
