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
use super::interpolation::*;
use super::{jd_to_spice_jd, spice_jds_to_jd, DafArray};
use crate::constants::AU_KM;
use crate::errors::NEOSpyError;
use crate::frames::Frame;
use crate::prelude::Desig;
use crate::state::State;
use std::fmt::Debug;

#[derive(Debug)]
pub enum SpkSegmentType {
    Type1(SpkSegmentType1),
    Type2(SpkSegmentType2),
    Type13(SpkSegmentType13),
    Type21(SpkSegmentType21),
}

impl SpkSegmentType {
    /// Create a Segment from a DafArray, the segment type must be specified.
    pub fn from_array(segment_type: i32, array: DafArray) -> Result<Self, NEOSpyError> {
        match segment_type {
            1 => Ok(SpkSegmentType::Type1(array.into())),
            2 => Ok(SpkSegmentType::Type2(array.into())),
            13 => Ok(SpkSegmentType::Type13(array.into())),
            21 => Ok(SpkSegmentType::Type21(array.into())),
            v => Err(NEOSpyError::IOError(format!(
                "SPK Segment type {:?} not supported.",
                v
            ))),
        }
    }
}

impl From<SpkSegmentType> for DafArray {
    fn from(value: SpkSegmentType) -> Self {
        match value {
            SpkSegmentType::Type1(seg) => seg.array,
            SpkSegmentType::Type2(seg) => seg.array,
            SpkSegmentType::Type13(seg) => seg.array,
            SpkSegmentType::Type21(seg) => seg.array,
        }
    }
}

impl From<i32> for Frame {
    fn from(value: i32) -> Self {
        match value {
            1 => Frame::Equatorial, // J2000
            17 => Frame::Ecliptic,  // ECLIPJ2000
            _ => Frame::Unknown(value as usize),
        }
    }
}

#[derive(Debug)]
pub struct SpkSegment {
    /// The NAIF ID of the object recorded in this Segment.
    pub obj_id: i64,

    /// Start time of the segment.
    pub jd_start: f64,

    /// End time of the segment.
    pub jd_end: f64,

    /// The reference center NAIF ID for the position/velocity in this Segment.
    pub center_id: i64,

    /// [`Frame`] of reference for this Segment.
    pub ref_frame: Frame,

    /// Number which defines the segment type as defined by the DAF standard.
    pub segment_type: i32,

    /// Internal data representation.
    segment: SpkSegmentType,
}

impl From<SpkSegment> for DafArray {
    fn from(value: SpkSegment) -> Self {
        value.segment.into()
    }
}

impl TryFrom<DafArray> for SpkSegment {
    type Error = NEOSpyError;

    fn try_from(value: DafArray) -> Result<SpkSegment, Self::Error> {
        let summary_floats = &value.summary_floats;
        let summary_ints = &value.summary_ints;
        let jd_start = spice_jds_to_jd(summary_floats[0]);
        let jd_end = spice_jds_to_jd(summary_floats[1]);
        let obj_id = summary_ints[0] as i64;
        let center_id = summary_ints[1] as i64;
        let frame_num = summary_ints[2];
        let segment_type = summary_ints[3];

        let ref_frame = frame_num.into();

        let segment = SpkSegmentType::from_array(segment_type, value)?;

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
}

impl SpkSegment {
    pub fn jd_range(&self) -> (f64, f64) {
        (self.jd_start, self.jd_end)
    }

    pub fn contains(&self, jd: f64) -> bool {
        (jd >= self.jd_start) && (jd <= self.jd_end)
    }

    /// Return the [`State`] object at the specified JD. If the requested time is
    /// not within the available range, this will fail.
    pub fn try_get_state(&self, jd: f64) -> Result<State, NEOSpyError> {
        // this is faster than calling contains, probably because the || instead of &&
        if jd < self.jd_start || jd > self.jd_end {
            return Err(NEOSpyError::DAFLimits(
                "JD is not present in this record.".to_string(),
            ));
        }

        let (pos, vel) = match &self.segment {
            SpkSegmentType::Type1(v) => v.try_get_pos_vel(self, jd)?,
            SpkSegmentType::Type2(v) => v.try_get_pos_vel(self, jd)?,
            SpkSegmentType::Type13(v) => v.try_get_pos_vel(self, jd)?,
            SpkSegmentType::Type21(v) => v.try_get_pos_vel(self, jd)?,
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
pub struct SpkSegmentType1 {
    array: DafArray,

    n_records: usize,
}

impl SpkSegmentType1 {
    fn get_record(&self, idx: usize) -> &[f64] {
        unsafe { self.array.data.get_unchecked(idx * 71..(idx + 1) * 71) }
    }

    fn get_times(&self) -> &[f64] {
        unsafe {
            self.array
                .data
                .get_unchecked(self.n_records * 71..(self.n_records * 72))
        }
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
            .get_times()
            .binary_search_by(|probe| probe.total_cmp(&jd))
            .unwrap_or_else(|c| c);

        let record = self.get_record(start_idx);

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

impl From<DafArray> for SpkSegmentType1 {
    fn from(array: DafArray) -> Self {
        let n_records = array[array.len() - 1] as usize;

        SpkSegmentType1 { array, n_records }
    }
}

/// Chebyshev Polynomials (Position Only)
///
/// <https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/spk.html#Type%202:%20Chebyshev%20position%20only>
///
#[derive(Debug)]
pub struct SpkSegmentType2 {
    array: DafArray,
    jd_step: f64,
    n_coef: usize,
    record_len: usize,
}

/// Type 2 Record View
/// A view into a record of type 2, provided mainly for clarity to the underlying
/// data structure.
struct Type2RecordView<'a> {
    t_mid: &'a f64,
    t_step: &'a f64,

    x_coef: &'a [f64],
    y_coef: &'a [f64],
    z_coef: &'a [f64],
}

impl SpkSegmentType2 {
    fn get_record(&self, idx: usize) -> Type2RecordView {
        unsafe {
            let vals = self
                .array
                .data
                .get_unchecked(idx * self.record_len..(idx + 1) * self.record_len);

            Type2RecordView {
                t_mid: vals.get_unchecked(0),
                t_step: vals.get_unchecked(1),
                x_coef: vals.get_unchecked(2..(self.n_coef + 2)),
                y_coef: vals.get_unchecked((self.n_coef + 2)..(2 * self.n_coef + 2)),
                z_coef: vals.get_unchecked((2 * self.n_coef + 2)..(3 * self.n_coef + 2)),
            }
        }
    }

    fn try_get_pos_vel(
        &self,
        segment: &SpkSegment,
        jd: f64,
    ) -> Result<([f64; 3], [f64; 3]), NEOSpyError> {
        let jd = jd_to_spice_jd(jd);
        let jd_start = jd_to_spice_jd(segment.jd_start);
        let record_index = ((jd - jd_start) / self.jd_step).floor() as usize;
        let record = self.get_record(record_index);

        let t_step = record.t_step;

        let t = (jd - record.t_mid) / t_step;

        let t_step_scaled = 86400.0 / t_step / AU_KM;

        let (p, v) = chebyshev3_evaluate_both(t, record.x_coef, record.y_coef, record.z_coef)?;
        Ok((
            [p[0] / AU_KM, p[1] / AU_KM, p[2] / AU_KM],
            [
                v[0] * t_step_scaled,
                v[1] * t_step_scaled,
                v[2] * t_step_scaled,
            ],
        ))
    }
}

impl From<DafArray> for SpkSegmentType2 {
    fn from(array: DafArray) -> Self {
        let record_len = array[array.len() - 2] as usize;
        let jd_step = array[array.len() - 3];

        let n_coef = (record_len - 2) / 3;

        if 3 * n_coef + 2 != record_len {
            panic!("File incorrectly formatted, found number of Chebyshev coefficients doesn't match expected");
        }

        SpkSegmentType2 {
            array,
            n_coef,
            record_len,
            jd_step,
        }
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
    array: DafArray,
    window_size: usize,
    n_records: usize,
}

impl From<DafArray> for SpkSegmentType13 {
    fn from(array: DafArray) -> Self {
        let n_records = array[array.len() - 1] as usize;
        let window_size = array[array.len() - 2] as usize;

        Self {
            array,
            window_size,
            n_records,
        }
    }
}

/// Type 13 Record View
/// A view into a record of type 13, provided mainly for clarity to the underlying
/// data structure.
struct Type13RecordView<'a> {
    pos: &'a [f64; 3],
    vel: &'a [f64; 3],
}

impl SpkSegmentType13 {
    fn get_record(&self, idx: usize) -> Type13RecordView {
        unsafe {
            let rec = self.array.data.get_unchecked(idx * 6..(idx + 1) * 6);
            Type13RecordView {
                pos: rec[0..3].try_into().unwrap(),
                vel: rec[3..6].try_into().unwrap(),
            }
        }
    }

    fn get_times(&self) -> &[f64] {
        unsafe {
            self.array
                .data
                .get_unchecked(self.n_records * 6..self.n_records * 7)
        }
    }

    fn try_get_pos_vel(
        &self,
        _: &SpkSegment,
        jd: f64,
    ) -> Result<([f64; 3], [f64; 3]), NEOSpyError> {
        let jd = jd_to_spice_jd(jd);
        let times = self.get_times();
        let start_idx: isize = match times.binary_search_by(|probe| probe.total_cmp(&jd)) {
            Ok(c) => c as isize - (self.window_size as isize) / 2,
            Err(c) => {
                if (jd - times[c - 1]).abs() < (jd - times[c]).abs() {
                    c as isize - 1 - self.window_size as isize / 2
                } else {
                    c as isize - self.window_size as isize / 2
                }
            }
        };
        let start_idx =
            start_idx.clamp(0, self.n_records as isize - self.window_size as isize) as usize;

        let mut pos = [0.0; 3];
        let mut vel = [0.0; 3];
        for idx in 0..3 {
            let p: Box<[f64]> = (0..self.window_size)
                .map(|i| self.get_record(i + start_idx).pos[idx])
                .collect();
            let dp: Box<[f64]> = (0..self.window_size)
                .map(|i| self.get_record(i + start_idx).vel[idx])
                .collect();
            let (p, v) =
                hermite_interpolation(&times[start_idx..start_idx + self.window_size], &p, &dp, jd);
            pos[idx] = p / AU_KM;
            vel[idx] = v / AU_KM * 86400.;
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
    array: DafArray,
    n_coef: usize,
    n_records: usize,
    record_len: usize,
}

impl From<DafArray> for SpkSegmentType21 {
    fn from(array: DafArray) -> Self {
        let n_records = array[array.len() - 1] as usize;
        let n_coef = array[array.len() - 2] as usize;

        let record_len = 4 * n_coef + 11;

        Self {
            array,
            n_coef,
            n_records,
            record_len,
        }
    }
}

impl SpkSegmentType21 {
    fn get_record(&self, idx: usize) -> &[f64] {
        unsafe {
            self.array
                .data
                .get_unchecked(idx * self.record_len..(idx + 1) * self.record_len)
        }
    }

    fn get_times(&self) -> &[f64] {
        unsafe {
            self.array.data.get_unchecked(
                self.n_records * self.record_len..self.n_records * (self.record_len + 1),
            )
        }
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
            .get_times()
            .binary_search_by(|probe| probe.total_cmp(&jd))
            .unwrap_or_else(|c| c);

        let record = self.get_record(start_idx);

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
