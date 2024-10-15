//! Loading and reading of states from JPL PCK kernel files.
//!
//! PCKs are intended to be loaded into a singleton which is accessible via the
//! [`get_pck_singleton`] function defined below. This singleton is wrapped in a RwLock,
//! meaning before its use it must by unwrapped. A vast majority of intended use cases
//! will only be the read case.
//!
use super::daf::{DAFType, DafFile};
use super::pck_segments::PckSegment;
use crate::errors::{Error, NeosResult};
use crate::frames::Frame;
use crossbeam::sync::ShardedLock;

use std::io::Cursor;
use std::mem::MaybeUninit;
use std::sync::Once;

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
    pub fn load_file(&mut self, filename: &str) -> NeosResult<()> {
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
    pub fn try_get_orientation(&self, id: isize, jd: f64) -> NeosResult<Frame> {
        for segment in self.segments.iter() {
            if id == segment.center_id && segment.contains(jd) {
                return segment.try_get_orientation(jd);
            }
        }

        Err(Error::DAFLimits(format!(
            "Object ({}) does not have an PCK record for the target JD.",
            id
        )))?
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

/// Get the PCK singleton.
/// This is a RwLock protected PCKCollection, and must be `.try_read().unwrapped()` for any
/// read-only cases.
pub fn get_pck_singleton() -> &'static PckSingleton {
    // Create an uninitialized static
    static mut SINGLETON: MaybeUninit<PckSingleton> = MaybeUninit::uninit();
    static ONCE: Once = Once::new();

    unsafe {
        ONCE.call_once(|| {
            let mut files = PckCollection {
                segments: Vec::new(),
            };
            files.reset();
            let singleton: PckSingleton = ShardedLock::new(files);
            // Store it to the static var, i.e. initialize it
            let _ = SINGLETON.write(singleton);
        });

        // Now we give out a shared reference to the data, which is safe to use concurrently.
        SINGLETON.assume_init_ref()
    }
}
