/// Most users should interface with `spk.rs`, not this module.
///
/// SPK Files are collections of `Segments`, which are ranges of times where the state
/// of an object is recorded. These segments are typically made up of many individual
/// `Records`, with an associated maximum and minimum time where they are valid for.
///
/// There are unique structs for each possible segment type, not all are currently
/// supported. Each segment type must implement the SPKSegment trait, which allows for
/// the loading and querying of states contained within.
///
/// <https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/spk.html#Supported%20Data%20Types>
///
/// There is a lot of repetition in this file, as many of the segment types have very
/// similar internal structures.
use super::binary::read_f64_vec;
use super::interpolation::*;
use super::records::DafRecords;
use super::{spice_jds_to_jd, jd_to_spice_jd};
use crate::constants::AU_KM;
use crate::errors::NEOSpyError;
use crate::frames::Frame;
use crate::prelude::Desig;
use crate::state::State;
use std::fmt::Debug;
use std::io::{Read, Seek};

#[derive(Debug)]
enum SPKSegmentType {
    Type1(SpkSegmentType1),
    Type2(SpkSegmentType2),
    Type3(SpkSegmentType3),
    Type13(SpkSegmentType13),
    Type21(SpkSegmentType21),
}

impl SPKSegmentType{
    pub fn try_load<T: Read + Seek>(segment_type:usize, file: &mut T, array_start: usize, array_end:usize )-> Result<Self, NEOSpyError>{
        match segment_type {
        1 => SpkSegmentType1::try_load(file, array_start, array_end),
        2 => SpkSegmentType2::try_load(file, array_start, array_end),
        3 => SpkSegmentType3::try_load(file, array_start, array_end),
        13 => SpkSegmentType13::try_load(file, array_start, array_end),
        21 => SpkSegmentType21::try_load(file, array_start, array_end),
        v => Err(NEOSpyError::IOError(format!(
                "SPK Segment type {:?} not supported.",
                v
            )))
        
    }}
}

#[derive(Debug)]
pub struct SpkSegment {
    /// The NAIF ID of the object recorded in this Segment.
    pub obj_id: isize,

    /// Start time of the segment.
    pub jd_start: f64,

    /// End time of the segment.
    pub jd_end: f64,

    /// The reference center NAIF ID for the position/velocity in this Segment.
    pub center_id: isize,

    /// [`Frame`] of reference for this Segment.
    pub ref_frame: Frame,

    /// Number which defines the segment type as defined by the DAF standard.
    pub segment_type: usize,

    /// Internal data representation.
    segment: SPKSegmentType,
}

impl SpkSegment {
    pub fn jd_range(&self) -> (f64, f64) {
        (self.jd_start, self.jd_end)
    }

    pub fn contains(&self, jd: f64) -> bool {
        (jd >= self.jd_start) && (jd <= self.jd_end)
    }

    /// Load this segment from the provided file and summary.
    pub fn from_summary<T: Read + Seek>(
        file: &mut T,
        summary: (Box<[f64]>, Box<[i32]>),
    ) -> Result<SpkSegment, NEOSpyError> {
        let (floats, ints) = summary;
        let jd_start = spice_jds_to_jd(floats[0]);
        let jd_end = spice_jds_to_jd(floats[1]);
        let obj_id = ints[0] as isize;
        let center_id = ints[1] as isize;
        let frame_num = ints[2];
        let segment_type = ints[3] as usize;
        let array_start = ints[4] as usize;
        let array_end = ints[5] as usize;

        let ref_frame = match frame_num {
            1 => Frame::Equatorial, // J2000
            17 => Frame::Ecliptic,  // ECLIPJ2000
            _ => Frame::Unknown(frame_num as usize),
        };

        let segment = SPKSegmentType::try_load(segment_type, file, array_start, array_end)?;

        Ok(SpkSegment {
            obj_id,
            jd_start,
            jd_end,
            center_id,
            ref_frame,
            segment_type,
            segment,
        })
    }

    /// Return the [`State`] object at the specified JD. If the requested time is
    /// not within the available range, this will fail.
    pub fn try_get_state(&self, jd: f64) -> Result<State, NEOSpyError> {
        if jd < self.jd_start || jd > self.jd_end {
            return Err(NEOSpyError::DAFLimits(
                "JD is not present in this record.".to_string(),
            ));
        }

        let (pos, vel) = match &self.segment {
            SPKSegmentType::Type1(v) => v.try_get_pos_vel(self, jd)?,
            SPKSegmentType::Type2(v) => v.try_get_pos_vel(self, jd)?,
            SPKSegmentType::Type3(v) => v.try_get_pos_vel(self, jd)?,
            SPKSegmentType::Type13(v) => v.try_get_pos_vel(self, jd)?,
            SPKSegmentType::Type21(v) => v.try_get_pos_vel(self, jd)?,
        };

        Ok(State::new(
            Desig::Naif(self.obj_id),
            jd,
            pos.into(),
            vel.into(),
            self.ref_frame,
            self.center_id,
        ))
    }
}

/// Modified Difference Arrays
///
/// <https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/spk.html#Type%201:%20Modified%20Difference%20Arrays>
///
#[derive(Debug)]
struct SpkSegmentType1 {
    records: DafRecords,
    times: Box<[f64]>,
}

impl SpkSegmentType1 {
    fn try_load<T: Read + Seek>(
        file: &mut T,
        array_start: usize,
        array_end: usize,
    ) -> Result<SPKSegmentType, NEOSpyError> {
        let _ = file.seek(std::io::SeekFrom::Start(8 * (array_start as u64 - 1)))?;
        let array_len = array_end - array_start + 1;
        let bin_segment = read_f64_vec(file, array_len, true)?;

        let n_records = bin_segment[array_len - 1] as usize;
        let record_len = 71;

        let mut records = DafRecords::with_capacity(n_records, record_len);
        for idx in 0..n_records {
            records.try_push(&bin_segment[idx * record_len..(idx + 1) * record_len])?;
        }

        let times = {
            let mut tmp = Vec::<f64>::with_capacity(n_records);
            for idy in (record_len * n_records)..(record_len * n_records + n_records) {
                tmp.push(bin_segment[idy]);
            }
            tmp
        };
        Ok(SPKSegmentType::Type1(SpkSegmentType1 {
            records,
            times: times.into(),
        }))
    }

    fn try_get_pos_vel(
        &self,
        _: &SpkSegment,
        jd: f64,
    ) -> Result<([f64; 3], [f64; 3]), NEOSpyError> {
        // Records are laid out as so:
        //
        // Size      Description
        // ----------------------
        // 1          Reference Epoch for the difference line
        // n_coef     Step size function vector
        // 6          Reference state - x, vx, y, vy, z, vz  (interleaved order)
        // 3*n_coef   Modified divided difference arrays
        // 1          Maximum integration order plus 1
        // 3          Integration order array
        // total: 11 + 4*n_coef
        // we need to find the first record which has a time greater than or equal
        // to the target jd.

        let jd = jd_to_spice_jd(jd);
        let start_idx = self
            .times
            .binary_search_by(|probe| probe.total_cmp(&jd))
            .unwrap_or_else(|c| c);

        let record = unsafe { self.records.get_unchecked(start_idx) };
        debug_assert!(record.len() == 71);

        let ref_time = record[0];

        let func_vec = &record[1..16];
        let ref_state = &record[16..22];

        let divided_diff_array = &record[22..67];

        let kq_max1 = record[67] as usize;
        let kq = &record[68..71];

        // in the spice code ref_time is in seconds from j2000
        let dt = jd - ref_time;

        let mut fc = [0.0; 15];
        let mut wc = [0.0; 15];

        let mut tp = dt;
        for idx in 0..(kq_max1 - 2) {
            let f = func_vec[idx];
            if f == 0.0 {
                // don't divide by 0 below, file was built incorrectly.
                return Err(NEOSpyError::IOError(
                    "SPK File contains segments of type 1 has invalid contents.".into(),
                ));
            }

            fc[idx] = tp / f;
            wc[idx] = dt / f;
            tp = dt + f;
        }

        let mut w: Box<[f64]> = { (0..kq_max1).map(|x| (x as f64 + 1.0).recip()).collect() };

        let mut ks = kq_max1 - 1;
        let mut jx = 0;
        let mut ks1 = ks - 1;

        while ks >= 2 {
            jx += 1;
            for j in 0..jx {
                w[j + ks] = fc[j] * w[j + ks1] - wc[j] * w[j + ks];
            }
            ks = ks1;
            ks1 -= 1;
        }

        // position interpolation
        let pos = std::array::from_fn(|idx| {
            let sum: f64 = (1..(kq[idx] as usize + 1))
                .rev()
                .map(|j| divided_diff_array[15 * idx + j - 1] * w[j + ks - 1])
                .sum();
            (ref_state[2 * idx] + dt * (sum * dt + ref_state[2 * idx + 1])) / AU_KM
        });

        // Recompute W for velocities
        for j in 0..jx {
            w[j + ks] = fc[j] * w[j + ks1] - wc[j] * w[j + ks];
        }
        ks -= 1;

        // velocity interpolation
        let vel = std::array::from_fn(|idx| {
            let sum: f64 = (1..(kq[idx] as usize + 1))
                .rev()
                .map(|j| divided_diff_array[15 * idx + j - 1] * w[j + ks - 1])
                .sum();
            (ref_state[2 * idx + 1] + dt * sum) / AU_KM * 86400.0
        });
        Ok((pos, vel))
    }
}

/// Chebyshev Polynomials (Position Only)
///
/// <https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/spk.html#Type%202:%20Chebyshev%20position%20only>
///
#[derive(Debug)]
pub struct SpkSegmentType2 {
    records: DafRecords,
    jd_step: f64,
    n_coef: usize,
}

impl SpkSegmentType2 {
    fn try_load<T: Read + Seek>(
        file: &mut T,
        array_start: usize,
        array_end: usize,
    ) -> Result<SPKSegmentType, NEOSpyError> {
        let _ = file.seek(std::io::SeekFrom::Start(8 * (array_start as u64 - 1)))?;
        let array_len = array_end - array_start + 1;
        let bin_segment = read_f64_vec(file, array_len, true)?;

        let n_records = bin_segment[array_len - 1] as usize;
        let record_len = bin_segment[array_len - 2] as usize;
        let jd_step = bin_segment[array_len - 3];

        let n_coef = (record_len - 2) / 3;

        if 3 * n_coef + 2 != record_len {
            return Err(NEOSpyError::DAFLimits("File incorrectly formatted, found number of Chebyshev coefficients doesn't match expected".into()));
        }

        let mut records = DafRecords::with_capacity(n_records, record_len);
        for idy in 0..n_records {
            records.try_push(&bin_segment[idy * record_len .. (1 + idy) * record_len])?;
        }

        Ok(SPKSegmentType::Type2(SpkSegmentType2 {
            jd_step,
            n_coef,
            records,
        }))
    }

    fn try_get_pos_vel(
        &self,
        segment: &SpkSegment,
        jd: f64,
    ) -> Result<([f64; 3], [f64; 3]), NEOSpyError> {
        let jd = jd_to_spice_jd(jd);
        let jd_start = jd_to_spice_jd(segment.jd_start);
        let record_index = ((jd - jd_start) / self.jd_step).floor() as usize;
        let record = unsafe { self.records.get_unchecked(record_index) };

        let x_coef = unsafe { record.get_unchecked(2..(self.n_coef + 2)) };
        let y_coef = unsafe { record.get_unchecked((self.n_coef + 2)..(2 * self.n_coef + 2)) };
        let z_coef = unsafe { record.get_unchecked((2 * self.n_coef + 2)..(3 * self.n_coef + 2)) };

        let t_mid = unsafe{record.get_unchecked(0)};
        let t_step = unsafe{record.get_unchecked(1)};
        let t = (jd - t_mid) / t_step;

        let t_step_scaled = 86400.0 / t_step / AU_KM;

        let (x, vx) = chebyshev_evaluate_both(t, x_coef)?;
        let (y, vy) = chebyshev_evaluate_both(t, y_coef)?;
        let (z, vz) = chebyshev_evaluate_both(t, z_coef)?;
        Ok(([x / AU_KM, y / AU_KM, z / AU_KM], [vx * t_step_scaled, vy * t_step_scaled, vz * t_step_scaled]))
    }
}

/// Chebyshev Polynomials (Position and Velocity)
///
/// <https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/spk.html#Type%203:%20Chebyshev%20position%20and%20velocity>
///
#[derive(Debug)]
pub struct SpkSegmentType3 {
    records: DafRecords,
    jd_step: f64,
    n_coef: usize,
}

impl SpkSegmentType3 {
    fn try_load<T: Read + Seek>(
        file: &mut T,
        array_start: usize,
        array_end: usize,
    ) -> Result<SPKSegmentType, NEOSpyError> {
        let _ = file.seek(std::io::SeekFrom::Start(8 * (array_start as u64 - 1)))?;
        let array_len = array_end - array_start + 1;
        let bin_segment = read_f64_vec(file, array_len, true)?;

        let n_records = bin_segment[array_len - 1] as usize;
        let record_len = bin_segment[array_len - 2] as usize;
        let jd_step = bin_segment[array_len - 3] / 86400.0;

        let n_coef = (record_len - 2) / 6;

        let mut records = DafRecords::with_capacity(n_records, record_len);
        for idy in 0..n_records {
            records.try_push(&bin_segment[idy * record_len .. (1 + idy) * record_len])?;
        }

        Ok(SPKSegmentType::Type3(SpkSegmentType3 {
            jd_step,
            n_coef,
            records,
        }))
    }

    fn try_get_pos_vel(
        &self,
        segment: &SpkSegment,
        jd: f64,
    ) -> Result<([f64; 3], [f64; 3]), NEOSpyError> {
        let jd = jd_to_spice_jd(jd);
        let jd_start = jd_to_spice_jd(segment.jd_start);
        let record_index = ((jd - jd_start) / self.jd_step).floor() as usize;
        let record = unsafe { self.records.get_unchecked(record_index) };
        let t_mid = record[0];
        let t_step = record[1];
        let t = (jd - t_mid) / t_step;

        let t_step_scaled = 86400.0 / t_step / AU_KM;

        let x_coef = &record[2..(self.n_coef + 2)];
        let y_coef = &record[(self.n_coef + 2)..(2 * self.n_coef + 2)];
        let z_coef = &record[(2 * self.n_coef + 2)..(3 * self.n_coef + 2)];

        let vx_coef = &record[(3 * self.n_coef + 2)..(4 * self.n_coef + 2)];
        let vy_coef = &record[(4 * self.n_coef + 2)..(5 * self.n_coef + 2)];
        let vz_coef = &record[(5 * self.n_coef + 2)..(6 * self.n_coef + 2)];

        let x = chebyshev_evaluate_type1(t, x_coef)?;
        let y = chebyshev_evaluate_type1(t, y_coef)?;
        let z = chebyshev_evaluate_type1(t, z_coef)?;

        let vx = chebyshev_evaluate_type1(t, vx_coef)?;
        let vy = chebyshev_evaluate_type1(t, vy_coef)?;
        let vz = chebyshev_evaluate_type1(t, vz_coef)?;

        Ok(([x / AU_KM, y / AU_KM, z / AU_KM], [vx * t_step_scaled, vy * t_step_scaled, vz * t_step_scaled]))
    }
}

// TODO: SPK Segment type 12 should be a minor variation on type 13. This was not
// implemented here due to missing a valid SPK file to test against.

/// Hermite Interpolation (Uneven Time Steps)
///
/// This uses a collection of individual positions/velocities and interpolates between
/// them using hermite interpolation.
/// <https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/spk.html#Type%2013:%20Hermite%20Interpolation%20---%20Unequal%20Time%20Steps>
#[derive(Debug)]
pub struct SpkSegmentType13 {
    records: DafRecords,
    times: Box<[f64]>,
    window_size: usize,
}
impl SpkSegmentType13 {
    fn try_load<T: Read + Seek>(
        file: &mut T,
        array_start: usize,
        array_end: usize,
    ) -> Result<SPKSegmentType, NEOSpyError> {
        let _ = file.seek(std::io::SeekFrom::Start(8 * (array_start as u64 - 1)))?;
        let array_len = array_end - array_start + 1;
        let bin_segment = read_f64_vec(file, array_len, true)?;

        let n_records = bin_segment[array_len - 1] as usize;
        let window_size = bin_segment[array_len - 2] as usize;
        let record_len = 6;

        let mut states = DafRecords::with_capacity(n_records, record_len);
        for idx in 0..n_records {
            states.try_push(
                & bin_segment[idx * 6..(idx + 1) * record_len])?;
        }

        let mut times = Vec::with_capacity(n_records);
        for idx in (n_records * record_len)..(n_records * record_len + record_len) {
            times.push(bin_segment[idx]);
        }

        Ok(SPKSegmentType::Type13(SpkSegmentType13 {
            records: states,
            times: times.into(),
            window_size,
        }))
    }

    fn try_get_pos_vel(
        &self,
        _: &SpkSegment,
        jd: f64,
    ) -> Result<([f64; 3], [f64; 3]), NEOSpyError> {
        let jd = jd_to_spice_jd(jd);
        let start_idx: isize = match self.times.binary_search_by(|probe| probe.total_cmp(&jd)) {
            Ok(c) => c as isize - (self.window_size as isize) / 2,
            Err(c) => {
                if (jd - self.times[c - 1]).abs() < (jd - self.times[c]).abs() {
                    c as isize - 1 - self.window_size as isize / 2
                } else {
                    c as isize - self.window_size as isize / 2
                }
            }
        };
        let start_idx =
            start_idx.clamp(0, self.records.len() as isize - self.window_size as isize) as usize;

        let mut pos = [0.0; 3];
        let mut vel = [0.0; 3];
        let times: Box<[&f64]> = self
            .times
            .iter()
            .skip(start_idx)
            .take(self.window_size)
            .collect();
        for idx in 0..3 {
            let p: Box<[f64]> = (0..self.window_size)
                .map(|i| self.records[i + start_idx][idx])
                .collect();
            let dp: Box<[f64]> = (0..self.window_size)
                .map(|i| self.records[i + start_idx][idx + 3])
                .collect();
            let (p, v) = hermite_interpolation(&times, &p, &dp, jd);
            pos[idx] = p/ AU_KM;
            vel[idx] = v/ AU_KM * 86400.;
        }

        Ok((pos, vel))
    }
}

/// Extended Modified Difference Arrays
///
/// <https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/spk.html#Type%2021:%20Extended%20Modified%20Difference%20Arrays>
///
#[derive(Debug)]
pub struct SpkSegmentType21 {
    records: DafRecords,
    times: Box<[f64]>,
    n_coef: usize,
}
impl SpkSegmentType21 {
    fn try_load<T: Read + Seek>(
        file: &mut T,
        array_start: usize,
        array_end: usize,
    ) -> Result<SPKSegmentType, NEOSpyError> {
        let _ = file.seek(std::io::SeekFrom::Start(8 * (array_start as u64 - 1)))?;
        let array_len = array_end - array_start + 1;
        let bin_segment = read_f64_vec(file, array_len, true)?;

        let n_records = bin_segment[array_len - 1] as usize;
        let n_coef = bin_segment[array_len - 2] as usize;

        let record_len = 4 * n_coef + 11;

        let mut records = DafRecords::with_capacity(n_records, record_len);
        for idx in 0..n_records {
            records.try_push(
                & bin_segment[idx * record_len..(idx + 1) * record_len])?;
        }

        let times = {
            let mut tmp = Vec::<f64>::with_capacity(n_records);
            for idy in (record_len * n_records)..(record_len * n_records + n_records) {
                tmp.push(bin_segment[idy]);
            }
            tmp
        };

        Ok(SPKSegmentType::Type21(SpkSegmentType21 {
            records,
            times: times.into(),
            n_coef,
        }))
    }

    fn try_get_pos_vel(
        &self,
        _: &SpkSegment,
        jd: f64,
    ) -> Result<([f64; 3], [f64; 3]), NEOSpyError> {
        // Records are laid out as so:
        //
        // Size      Description
        // ----------------------
        // 1          Reference Epoch for the difference line
        // n_coef     Step size function vector
        // 6          Reference state - x, vx, y, vy, z, vz  (interleaved order)
        // 3*n_coef   Modified divided difference arrays
        // 1          Maximum integration order plus 1
        // 3          Integration order array
        // total: 11 + 4*n_coef

        // we need to find the first record which has a time greater than or equal
        // to the target jd.

        let jd = jd_to_spice_jd(jd);
        let start_idx = self
            .times
            .binary_search_by(|probe| probe.total_cmp(&jd))
            .unwrap_or_else(|c| c);

        let record = unsafe { self.records.get_unchecked(start_idx) };

        let ref_time = record[0];

        let func_vec = &record[1..self.n_coef + 1];
        let ref_state = &record[self.n_coef + 1..self.n_coef + 7];

        let divided_diff_array = &record[self.n_coef + 7..4 * self.n_coef + 7];

        let kq_max1 = record[4 * self.n_coef + 7] as usize;
        let kq = &record[4 * self.n_coef + 8..4 * self.n_coef + 11];

        // in the spice code ref_time is in seconds from j2000
        let dt = jd - ref_time;

        let mut fc = Vec::<f64>::with_capacity(self.n_coef);
        let mut wc = Vec::<f64>::with_capacity(self.n_coef);

        let mut tp = dt;
        for f in func_vec.iter().take(kq_max1 - 2) {
            if *f == 0.0 {
                // don't divide by 0 below, file was built incorrectly.
                return Err(NEOSpyError::IOError(
                    "SPK File contains segments of type 21 has invalid contents.".into(),
                ));
            }

            fc.push(tp / f);
            wc.push(dt / f);
            tp = dt + f;
        }

        let mut w: Box<[f64]> = (0..kq_max1).map(|x| (x as f64 + 1.0).recip()).collect();

        let mut ks = kq_max1 - 1;
        let mut jx = 0;
        let mut ks1 = ks - 1;

        while ks >= 2 {
            jx += 1;
            for j in 0..jx {
                w[j + ks] = fc[j] * w[j + ks1] - wc[j] * w[j + ks];
            }
            ks = ks1;
            ks1 -= 1;
        }

        // position interpolation
        let pos = std::array::from_fn(|idx| {
            let sum: f64 = (1..(kq[idx] as usize + 1))
                .rev()
                .map(|j| divided_diff_array[idx * self.n_coef + j - 1] * w[j + ks - 1])
                .sum();
            (ref_state[2 * idx] + dt * (sum * dt + ref_state[2 * idx + 1])) / AU_KM
        });

        // Recompute W for velocities
        for j in 0..jx {
            w[j + ks] = fc[j] * w[j + ks1] - wc[j] * w[j + ks];
        }
        ks -= 1;

        // velocity interpolation
        let vel = std::array::from_fn(|idx| {
            let sum: f64 = (1..(kq[idx] as usize + 1))
                .rev()
                .map(|j| divided_diff_array[idx * self.n_coef + j - 1] * w[j + ks - 1])
                .sum();
            (ref_state[2 * idx + 1] + dt * sum) / AU_KM * 86400.0
        });

        Ok((pos, vel))
    }
}
