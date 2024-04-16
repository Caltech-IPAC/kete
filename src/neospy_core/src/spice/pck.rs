/// Loading and reading of states from JPL PCK kernel files.
///
/// PCKs are intended to be loaded into a singleton which is accessible via the
/// [`get_pck_singleton`] function defined below. This singleton is wrapped in a RwLock,
/// meaning before its use it must by unwrapped. A vast majority of intended use cases
/// will only be the read case.
/// ```
///
///
use super::daf::{DAFType, Daf};
use super::pck_segments::*;
use crate::errors::NEOSpyError;
use crate::frames::Frame;

use std::io::{Cursor, Read, Seek};
use std::mem::MaybeUninit;
use std::sync::{Once, RwLock};

const PRELOAD_PCK: &[&[u8]] = &[
    include_bytes!("../../data/earth_000101_240215_231123.bpc"),
    include_bytes!("../../data/earth_200101_990825_predict.bpc"),
];

/// A collection of segments.
#[derive(Debug)]
pub struct PckSegmentCollection {
    segments: Vec<PckSegment>,
}

/// Define the PCK singleton structure.
pub type PckSingleton = RwLock<PckSegmentCollection>;

impl PckSegmentCollection {
    /// Given an PCK filename, load all the segments present inside of it.
    /// These segments are added to the PCK singleton in memory.
    pub fn load_file(&mut self, filename: &str) -> Result<(), NEOSpyError> {
        let mut file = std::fs::File::open(filename)?;
        let mut buffer = Vec::new();
        let _ = file.read_to_end(&mut buffer)?;
        let mut buffer = Cursor::new(&buffer);
        self.load_segments(&mut buffer)
    }

    /// Given an PCK buffer, load all the segments present inside of it.
    /// These segments are added to the PCK singleton in memory.
    pub fn load_segments<T: Read + Seek>(&mut self, mut buffer: T) -> Result<(), NEOSpyError> {
        let daf = Daf::try_load_header(&mut buffer)?;
        if daf.daf_type != DAFType::Pck {
            return Err(NEOSpyError::IOError(
                "Attempted to load a DAF file which is not an PCK as an PCK.".into(),
            ));
        }

        let summaries = daf.try_load_summaries(&mut buffer)?;

        for summary in summaries {
            self.segments
                .push(PckSegment::from_summary(&mut buffer, summary)?)
        }

        Ok(())
    }

    /// Get the raw orientation from the loaded PCK files.
    /// This orientation will have the frame of what was originally present in the file.
    pub fn try_get_orientation(&self, id: isize, jd: f64) -> Result<Frame, NEOSpyError> {
        for segment in self.segments.iter() {
            if id == segment.center_id && segment.contains(jd) {
                return segment.try_get_orientation(jd);
            }
        }
        Err(NEOSpyError::DAFLimits(format!(
            "Object ({}) does not have an PCK record for the target JD.",
            id
        )))
    }

    /// Delete all segments in the PCK singleton, equivalent to unloading all files.
    pub fn reset(&mut self) {
        let segments: PckSegmentCollection = PckSegmentCollection {
            segments: Vec::new(),
        };

        *self = segments;

        for preload in PRELOAD_PCK {
            let mut precise = Cursor::new(preload);
            self.load_segments(&mut precise).unwrap();
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
            let mut segments: PckSegmentCollection = PckSegmentCollection {
                segments: Vec::new(),
            };
            segments.reset();
            let singleton: PckSingleton = RwLock::new(segments);
            // Store it to the static var, i.e. initialize it
            let _ = SINGLETON.write(singleton);
        });

        // Now we give out a shared reference to the data, which is safe to use concurrently.
        SINGLETON.assume_init_ref()
    }
}
