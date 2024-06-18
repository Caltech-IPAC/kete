/// Loading and reading of states from JPL SPK kernel files.
///
/// SPKs are intended to be loaded into a singleton which is accessible via the
/// [`get_spk_singleton`] function defined below. This singleton is wrapped in a RwLock,
/// meaning before its use it must by unwrapped. A vast majority of intended use cases
/// will only be the read case.
///
/// Here is a small worked example:
/// ```
///     use neospy_core::spice::get_spk_singleton;
///     use neospy_core::frames::Frame;
///
///     // get a read-only reference to the [`SegmentCollection`]
///     let singleton = get_spk_singleton().try_read().unwrap();
///
///     // get the state of 399 (Earth) with respect to the Sun (10)
///     let state = singleton.try_get_state(399, 2451545.0, 10, Frame::Ecliptic);
/// ```
///
///
use super::daf::DafFile;
use super::{spk_segments::*, DAFType};
use crate::errors::NEOSpyError;
use crate::frames::Frame;
use crate::state::State;
use pathfinding::prelude::dijkstra;
use std::collections::{HashMap, HashSet};

use crossbeam::sync::ShardedLock;
use std::io::Cursor;
use std::mem::MaybeUninit;
use std::sync::Once;

const PRELOAD_SPKS: &[&[u8]] = &[
    include_bytes!("../../data/de440s.bsp"),
    include_bytes!("../../data/20000001.bsp"),
    include_bytes!("../../data/20000002.bsp"),
    include_bytes!("../../data/20000004.bsp"),
    include_bytes!("../../data/20000010.bsp"),
    include_bytes!("../../data/20000704.bsp"),
    include_bytes!("../../data/wise.bsp"),
];

/// A collection of segments.
/// This does not contain the full SPK file information, it is just a vector containing
/// all of the segments within the SPK files.
#[derive(Debug)]
pub struct SpkCollection {
    /// Collection of SPK Segment information
    segments: Vec<SpkSegment>,

    map_cache: HashMap<(i32, i32), Vec<i32>>,

    /// Map from object id to all connected pairs.
    nodes: HashMap<i32, HashSet<(i32, i32)>>,
}

/// Define the SPK singleton structure.
pub type SpkSingleton = ShardedLock<SpkCollection>;

impl SpkCollection {
    /// Get the raw state from the loaded SPK files.
    /// This state will have the center and frame of whatever was originally loaded
    /// into the file.
    pub fn try_get_raw_state(&self, id: i32, jd: f64) -> Result<State, NEOSpyError> {
        for segment in self.segments.iter() {
            if id == segment.obj_id && segment.contains(jd) {
                return segment.try_get_state(jd);
            }
        }
        Err(NEOSpyError::DAFLimits(
            format!(
                "Object ({}) does not have an SPK record for the target JD.",
                id
            )
            .to_string(),
        ))
    }

    /// Load a state from the file, then attempt to change the center to the center id
    /// specified.
    pub fn try_get_state(
        &self,
        id: i32,
        jd: f64,
        center: isize,
        frame: Frame,
    ) -> Result<State, NEOSpyError> {
        let mut state = self.try_get_raw_state(id, jd)?;
        self.try_change_center(&mut state, center)?;
        state.try_change_frame_mut(frame)?;
        Ok(state)
    }

    /// Use the data loaded in the SPKs to change the center ID of the provided state.
    pub fn try_change_center(&self, state: &mut State, new_center: isize) -> Result<(), NEOSpyError> {
        if state.center_id == new_center {
            return Ok(());
        }

        let path = self.find_path(state.center_id as i32, new_center as i32)?;

        for intermediate in path {
            let next = self.try_get_raw_state(intermediate, state.jd)?;
            state.try_change_center(next)?;
        }
        Ok(())
    }

    /// For a given NAIF ID, return all increments of time which are currently loaded.
    pub fn available_info(&self, id: i32) -> Vec<(f64, f64, i32, Frame, i32)> {
        let mut segment_info = Vec::<(f64, f64, i32, Frame, i32)>::new();
        for segment in self.segments.iter() {
            if id == segment.obj_id {
                let jd_range = segment.jd_range();
                segment_info.push((
                    jd_range.0,
                    jd_range.1,
                    segment.center_id,
                    segment.ref_frame,
                    segment.segment_type,
                ))
            }
        }
        if segment_info.is_empty() {
            return segment_info;
        }

        segment_info.sort_by(|a, b| (a.0).total_cmp(&b.0));

        let mut avail_times = Vec::<(f64, f64, i32, Frame, i32)>::new();

        let mut cur_segment = segment_info[0];
        for segment in segment_info.iter().skip(1) {
            // if the segments are overlapped or nearly overlapped, join them together
            // 1e-8 is approximately a millisecond
            if cur_segment.1 <= (segment.0 - 1e-8) {
                avail_times.push(cur_segment);
                cur_segment = *segment;
            } else {
                cur_segment.1 = segment.1.max(cur_segment.1)
            }
        }
        avail_times.push(cur_segment);

        avail_times
    }

    /// Return a hash set of all unique identifies loaded in the SPKs.
    /// If include centers is true, then this additionally includes the IDs for the
    /// center IDs. For example, if include_centers is false, then `0` will never be
    /// included in the loaded objects set, as 0 is a privileged position at the
    /// barycenter of the solar system. It is not typically defined in relation to
    /// anything else.
    pub fn loaded_objects(&self, include_centers: bool) -> HashSet<i32> {
        let mut found = HashSet::new();

        for segment in self.segments.iter() {
            let _ = found.insert(segment.obj_id);
            if include_centers {
                let _ = found.insert(segment.center_id);
            }
        }
        found
    }

    /// Given a NAIF ID, and a target NAIF ID, find the intermediate SPICE Segments
    /// which need to be loaded to find a path from one object to the other.
    /// Use Dijkstra plus the known segments to calculate a path.
    fn find_path(&self, start: i32, goal: i32) -> Result<Vec<i32>, NEOSpyError> {
        // first we check to see if the cache contains the lookup we need.
        if let Some(path) = self.map_cache.get(&(start, goal)) {
            return Ok(path.clone());
        }

        // not in the cache, manually compute
        let nodes = &self.nodes;
        let result = dijkstra(
            &(start, i32::MIN),
            |&current| match nodes.get(&current.0) {
                Some(set) => set.iter().map(|p| (*p, 1_i32)).collect(),
                None => Vec::<((i32, i32), i32)>::new(),
            },
            |&p| p.0 == goal,
        );
        if let Some((v, _)) = result {
            Ok(v.iter().skip(1).map(|x| x.1).collect())
        } else {
            Err(NEOSpyError::DAFLimits(format!(
                "SPK files are missing information to be able to map from obj {} to obj {}",
                start, goal
            )))
        }
    }

    /// Return all mappings from one object to another.
    ///
    /// These mappings are used to be able to change the center ID from whatever is saved in
    /// the spks to any possible combination.
    pub fn build_cache(&mut self) {
        static PRECACHE: &[i32] = &[0, 10, 399];

        let mut nodes: HashMap<i32, HashSet<(i32, i32)>> = HashMap::new();

        fn update_nodes(segment: &SpkSegment, nodes: &mut HashMap<i32, HashSet<(i32, i32)>>) {
            if let std::collections::hash_map::Entry::Vacant(e) = nodes.entry(segment.obj_id) {
                let mut set = HashSet::new();
                let _ = set.insert((segment.center_id, segment.obj_id));
                let _ = e.insert(set);
            } else {
                let _ = nodes
                    .get_mut(&segment.obj_id)
                    .unwrap()
                    .insert((segment.center_id, segment.obj_id));
            }
            if let std::collections::hash_map::Entry::Vacant(e) = nodes.entry(segment.center_id) {
                let mut set = HashSet::new();
                let _ = set.insert((segment.obj_id, segment.obj_id));
                let _ = e.insert(set);
            } else {
                let _ = nodes
                    .get_mut(&segment.center_id)
                    .unwrap()
                    .insert((segment.obj_id, segment.obj_id));
            }
        }

        for segment in self.segments.iter() {
            update_nodes(segment, &mut nodes);
        }

        let loaded = self.loaded_objects(true);

        for &start in loaded.iter() {
            for &goal in PRECACHE {
                let key = (start, goal);

                if self.map_cache.contains_key(&key) {
                    continue;
                }

                let result = dijkstra(
                    &(start, -100_i32),
                    |&current| match nodes.get(&current.0) {
                        Some(set) => set.iter().map(|p| (*p, 1_i32)).collect(),
                        None => Vec::<((i32, i32), i32)>::new(),
                    },
                    |&p| p.0 == goal,
                );

                if let Some((v, _)) = result {
                    let v: Vec<i32> = v.iter().skip(1).map(|x| x.1).collect();
                    let _ = self.map_cache.insert(key, v);
                }
            }
        }

        self.nodes = nodes;
    }

    /// Given an SPK filename, load all the segments present inside of it.
    /// These segments are added to the SPK singleton in memory.
    pub fn load_file(&mut self, filename: &str) -> Result<(), NEOSpyError> {
        let file = DafFile::from_file(filename)?;

        if !matches!(file.daf_type, DAFType::Spk) {
            return Err(NEOSpyError::IOError(format!(
                "File {:?} is not a PCK formatted file.",
                filename
            )));
        }
        self.segments
            .extend(file.segments.into_iter().map(|x| x.spk()));
        Ok(())
    }

    /// Delete all segments in the SPK singleton, equivalent to unloading all files.
    pub fn reset(&mut self) {
        let spk_files: SpkCollection = SpkCollection {
            map_cache: HashMap::new(),
            nodes: HashMap::new(),
            segments: Vec::new(),
        };

        *self = spk_files;

        for buffer in PRELOAD_SPKS {
            let mut curse = Cursor::new(buffer);
            let file = DafFile::from_buffer(&mut curse).unwrap();
            self.segments
                .extend(file.segments.into_iter().map(|x| x.spk()));
        }
        self.build_cache();
    }
}

/// Get the SPK singleton.
/// This is a RwLock protected SPKCollection, and must be `.try_read().unwrapped()` for any
/// read-only cases.
///
/// This singleton starts initialized with preloaded SPK files for the planets.
pub fn get_spk_singleton() -> &'static SpkSingleton {
    // Create an uninitialized static
    static mut SINGLETON: MaybeUninit<SpkSingleton> = MaybeUninit::uninit();
    static ONCE: Once = Once::new();

    unsafe {
        ONCE.call_once(|| {
            let mut segments: SpkCollection = SpkCollection {
                map_cache: HashMap::new(),
                nodes: HashMap::new(),
                segments: Vec::new(),
            };
            segments.reset();
            let singleton: SpkSingleton = ShardedLock::new(segments);
            // Store it to the static var, i.e. initialize it
            let _ = SINGLETON.write(singleton);
        });

        // Now we give out a shared reference to the data, which is safe to use
        // concurrently.
        SINGLETON.assume_init_ref()
    }
}
