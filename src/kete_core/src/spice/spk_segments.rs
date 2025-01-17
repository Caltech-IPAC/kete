//! Most users should interface with `spk.rs`, not this module.
//!
//! SPK Files are collections of `Segments`, which are ranges of times where the state
//! of an object is recorded. These segments are typically made up of many individual
//! `Records`, with an associated maximum and minimum time where they are valid for.
//!
//! There are unique structs for each possible segment type, not all are currently
//! supported. Each segment type must implement the SPKSegment trait, which allows for
//! the loading and querying of states contained within.
//!
//! <https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/spk.html#Supported%20Data%20Types>
//!
//! There is a lot of repetition in this file, as many of the segment types have very
//! similar internal structures.
use super::interpolation::*;
use super::{jd_to_spice_jd, spice_jds_to_jd, DafArray};
use crate::constants::AU_KM;
use crate::errors::Error;
use crate::frames::InertialFrame;
use crate::prelude::{Desig, KeteResult};
use crate::state::State;
use crate::time::scales::TDB;
use crate::time::Time;
use itertools::Itertools;
use sgp4::{julian_years_since_j2000, Constants, Geopotential, MinutesSinceEpoch, Orbit};
use std::fmt::Debug;
use std::marker::PhantomData;

#[derive(Debug)]
pub enum SpkSegmentType {
    Type1(SpkSegmentType1),
    Type2(SpkSegmentType2),
    Type3(SpkSegmentType3),
    Type9(SpkSegmentType9),
    Type10(SpkSegmentType10),
    Type13(SpkSegmentType13),
    Type18(SpkSegmentType18),
    Type21(SpkSegmentType21),
}

impl SpkSegmentType {
    /// Create a Segment from a DafArray, the segment type must be specified.
    pub fn from_array(segment_type: i32, array: DafArray) -> KeteResult<Self> {
        match segment_type {
            1 => Ok(SpkSegmentType::Type1(array.into())),
            2 => Ok(SpkSegmentType::Type2(array.into())),
            3 => Ok(SpkSegmentType::Type3(array.into())),
            9 => Ok(SpkSegmentType::Type9(array.into())),
            13 => Ok(SpkSegmentType::Type13(array.into())),
            10 => Ok(SpkSegmentType::Type10(array.into())),
            18 => Ok(SpkSegmentType::Type18(array.into())),
            21 => Ok(SpkSegmentType::Type21(array.into())),
            v => Err(Error::IOError(format!(
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
            SpkSegmentType::Type3(seg) => seg.array,
            SpkSegmentType::Type9(seg) => seg.array,
            SpkSegmentType::Type10(seg) => seg.array.array,
            SpkSegmentType::Type13(seg) => seg.array,
            SpkSegmentType::Type18(seg) => seg.array,
            SpkSegmentType::Type21(seg) => seg.array,
        }
    }
}

#[derive(Debug)]
pub struct SpkSegment<T: InertialFrame> {
    /// The NAIF ID of the object recorded in this Segment.
    pub obj_id: i64,

    /// Start time of the segment.
    pub jd_start: f64,

    /// End time of the segment.
    pub jd_end: f64,

    /// The reference center NAIF ID for the position/velocity in this Segment.
    pub center_id: i64,

    /// Number which defines the segment type as defined by the DAF standard.
    pub segment_type: i32,

    /// SPICE frame number.
    pub frame_num: i32,

    /// Internal data representation.
    segment: SpkSegmentType,

    frame: PhantomData<T>,
}

impl<T: InertialFrame> From<SpkSegment<T>> for DafArray {
    fn from(value: SpkSegment<T>) -> Self {
        value.segment.into()
    }
}

impl<T: InertialFrame> TryFrom<DafArray> for SpkSegment<T> {
    type Error = Error;

    fn try_from(value: DafArray) -> KeteResult<SpkSegment<T>> {
        let (jd_start, jd_end, obj_id, center_id, segment_type, frame_num) =
            value.spk_read_summary();

        let segment = SpkSegmentType::from_array(segment_type, value)?;

        Ok(SpkSegment::<T> {
            obj_id: obj_id as i64,
            jd_start,
            jd_end,
            center_id: center_id as i64,
            segment_type,
            frame_num,
            segment,
            frame: PhantomData,
        })
    }
}

impl<T: InertialFrame> SpkSegment<T> {
    #[inline(always)]
    pub fn jd_range(&self) -> (f64, f64) {
        (self.jd_start, self.jd_end)
    }

    #[inline(always)]
    pub fn contains(&self, jd: f64) -> bool {
        (jd >= self.jd_start) && (jd <= self.jd_end)
    }

    /// Return the [`State`] object at the specified JD. If the requested time is
    /// not within the available range, this will fail.
    #[inline(always)]
    pub fn try_get_state(&self, jd: f64) -> KeteResult<State<T>> {
        // this is faster than calling contains, probably because the || instead of &&
        if jd < self.jd_start || jd > self.jd_end {
            return Err(Error::DAFLimits(
                "JD is not present in this record.".to_string(),
            ));
        }

        let (pos, vel) = match &self.segment {
            SpkSegmentType::Type1(v) => v.try_get_pos_vel(self, jd)?,
            SpkSegmentType::Type2(v) => v.try_get_pos_vel(self, jd)?,
            SpkSegmentType::Type3(v) => v.try_get_pos_vel(self, jd)?,
            SpkSegmentType::Type9(v) => v.try_get_pos_vel(self, jd)?,
            SpkSegmentType::Type10(v) => v.try_get_pos_vel(self, jd)?,
            SpkSegmentType::Type13(v) => v.try_get_pos_vel(self, jd)?,
            SpkSegmentType::Type18(v) => v.try_get_pos_vel(self, jd)?,
            SpkSegmentType::Type21(v) => v.try_get_pos_vel(self, jd)?,
        };
        let des = Desig::Naif(self.obj_id);
        Ok(State::<T>::new(des, jd, pos.into(), vel.into(), self.center_id).into_frame())
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
    #[inline(always)]
    fn get_record(&self, idx: usize) -> &[f64] {
        unsafe { self.array.data.get_unchecked(idx * 71..(idx + 1) * 71) }
    }

    #[inline(always)]
    fn get_times(&self) -> &[f64] {
        unsafe {
            self.array
                .data
                .get_unchecked(self.n_records * 71..(self.n_records * 72))
        }
    }

    #[inline(always)]
    fn try_get_pos_vel<T: InertialFrame>(
        &self,
        _: &SpkSegment<T>,
        jd: f64,
    ) -> KeteResult<([f64; 3], [f64; 3])> {
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
                Err(Error::IOError(
                    "SPK File contains segments of type 1 has invalid contents.".into(),
                ))?;
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
    #[inline(always)]
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

    #[inline(always)]
    fn try_get_pos_vel<T: InertialFrame>(
        &self,
        segment: &SpkSegment<T>,
        jd: f64,
    ) -> KeteResult<([f64; 3], [f64; 3])> {
        let jd = jd_to_spice_jd(jd);
        let jd_start = jd_to_spice_jd(segment.jd_start);
        let record_index = ((jd - jd_start) / self.jd_step).floor() as usize;
        let record = self.get_record(record_index);

        let t_step = record.t_step;

        let t = (jd - record.t_mid) / t_step;

        let t_step_scaled = 86400.0 / t_step / AU_KM;

        let (p, v) = chebyshev_evaluate_both(t, record.x_coef, record.y_coef, record.z_coef)?;
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

/// Type 3 - Chebyshev Polynomials (Position & Velocity)
///
/// <https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/FORTRAN/req/spk.html#Type%203:%20Chebyshev%20position%20and%20velocity>
///
#[derive(Debug)]
pub struct SpkSegmentType3 {
    array: DafArray,
    jd_step: f64,
    n_coef: usize,
    record_len: usize,
}

/// Type 3 Record View
/// A view into a record of type 3, provided mainly for clarity to the underlying
/// data structure.
struct Type3RecordView<'a> {
    t_mid: &'a f64,
    t_step: &'a f64,

    x_coef: &'a [f64],
    y_coef: &'a [f64],
    z_coef: &'a [f64],

    vx_coef: &'a [f64],
    vy_coef: &'a [f64],
    vz_coef: &'a [f64],
}

impl SpkSegmentType3 {
    #[inline(always)]
    fn get_record(&self, idx: usize) -> Type3RecordView {
        unsafe {
            let vals = self
                .array
                .data
                .get_unchecked(idx * self.record_len..(idx + 1) * self.record_len);

            Type3RecordView {
                t_mid: vals.get_unchecked(0),
                t_step: vals.get_unchecked(1),
                x_coef: vals.get_unchecked(2..(self.n_coef + 2)),
                y_coef: vals.get_unchecked((self.n_coef + 2)..(2 * self.n_coef + 2)),
                z_coef: vals.get_unchecked((2 * self.n_coef + 2)..(3 * self.n_coef + 2)),
                vx_coef: vals.get_unchecked((3 * self.n_coef + 2)..(4 * self.n_coef + 2)),
                vy_coef: vals.get_unchecked((4 * self.n_coef + 2)..(5 * self.n_coef + 2)),
                vz_coef: vals.get_unchecked((5 * self.n_coef + 2)..(6 * self.n_coef + 2)),
            }
        }
    }

    #[inline(always)]
    fn try_get_pos_vel<T: InertialFrame>(&self, segment: &SpkSegment<T>, jd: f64) -> KeteResult<([f64; 3], [f64; 3])> {
        let jd = jd_to_spice_jd(jd);
        let jd_start = jd_to_spice_jd(segment.jd_start);
        let record_index = ((jd - jd_start) / self.jd_step).floor() as usize;
        let record = self.get_record(record_index);

        let t_step = record.t_step;

        let t = (jd - record.t_mid) / t_step;

        let t_scaled = 86400.0 / AU_KM;

        let p = chebyshev_evaluate(t, record.x_coef, record.y_coef, record.z_coef)?;
        let v = chebyshev_evaluate(t, record.vx_coef, record.vy_coef, record.vz_coef)?;
        Ok((
            [p[0] / AU_KM, p[1] / AU_KM, p[2] / AU_KM],
            [v[0] * t_scaled, v[1] * t_scaled, v[2] * t_scaled],
        ))
    }
}

impl From<DafArray> for SpkSegmentType3 {
    fn from(array: DafArray) -> Self {
        let record_len = array[array.len() - 2] as usize;
        let jd_step = array[array.len() - 3];

        let n_coef = (record_len - 2) / 6;

        if 6 * n_coef + 2 != record_len {
            panic!("File incorrectly formatted, found number of Chebyshev coefficients doesn't match expected");
        }

        SpkSegmentType3 {
            array,
            n_coef,
            record_len,
            jd_step,
        }
    }
}

// TODO: SPK Segment type 8 should be a minor variation on type 9. This was not
// implemented here due to missing a valid SPK file to test against.

/// Lagrange Interpolation (Uneven Time Steps)
///
/// This uses a collection of individual positions/velocities and interpolates between
/// them using Lagrange interpolation.
/// <https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/FORTRAN/req/spk.html#Type%209:%20Lagrange%20Interpolation%20---%20Unequal%20Time%20Steps>
#[derive(Debug)]
pub struct SpkSegmentType9 {
    array: DafArray,
    poly_degree: usize,
    n_records: usize,
}

impl From<DafArray> for SpkSegmentType9 {
    fn from(array: DafArray) -> Self {
        let n_records = array[array.len() - 1] as usize;
        let poly_degree = array[array.len() - 2] as usize;

        Self {
            array,
            poly_degree,
            n_records,
        }
    }
}

/// Type 9 Record View
/// A view into a record of type 9, provided mainly for clarity to the underlying
/// data structure.
struct Type9RecordView<'a> {
    pos: &'a [f64; 3],
    vel: &'a [f64; 3],
}

impl SpkSegmentType9 {
    #[inline(always)]
    fn get_record(&self, idx: usize) -> Type9RecordView {
        unsafe {
            let rec = self.array.data.get_unchecked(idx * 6..(idx + 1) * 6);
            Type9RecordView {
                pos: rec[0..3].try_into().unwrap(),
                vel: rec[3..6].try_into().unwrap(),
            }
        }
    }

    #[inline(always)]
    fn get_times(&self) -> &[f64] {
        unsafe {
            self.array
                .data
                .get_unchecked(self.n_records * 6..self.n_records * 7)
        }
    }

    #[inline(always)]
    fn try_get_pos_vel<T: InertialFrame>(&self, _: &SpkSegment<T>, jd: f64) -> KeteResult<([f64; 3], [f64; 3])> {
        let jd = jd_to_spice_jd(jd);
        let times = self.get_times();
        let window_size = self.poly_degree + 1;
        let start_idx: isize = match times.binary_search_by(|probe| probe.total_cmp(&jd)) {
            Ok(c) => c as isize - (window_size as isize) / 2,
            Err(c) => {
                if (jd - times[c - 1]).abs() < (jd - times[c]).abs() {
                    c as isize - 1 - window_size as isize / 2
                } else {
                    c as isize - window_size as isize / 2
                }
            }
        };
        let start_idx = start_idx.clamp(0, (self.n_records - window_size) as isize) as usize;

        let mut pos = [0.0; 3];
        let mut vel = [0.0; 3];
        for idx in 0..3 {
            let mut p: Box<[f64]> = (0..window_size)
                .map(|i| self.get_record(i + start_idx).pos[idx])
                .collect();
            let mut dp: Box<[f64]> = (0..window_size)
                .map(|i| self.get_record(i + start_idx).vel[idx])
                .collect();
            let p = lagrange_interpolation(&times[start_idx..start_idx + window_size], &mut p, jd);
            let v = lagrange_interpolation(&times[start_idx..start_idx + window_size], &mut dp, jd);
            pos[idx] = p / AU_KM;
            vel[idx] = v / AU_KM * 86400.;
        }

        Ok((pos, vel))
    }
}

/// Space Command two-line elements
///
/// <https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/spk.html#Type%2010:%20Space%20Command%20Two-Line%20Elements>
///
#[derive(Debug)]
pub struct SpkSegmentType10 {
    /// Generic Segments are a collection of a few different directories:
    /// `Packets` are where Type 10 stores the TLE values.
    /// `Packet Directory` is unused.
    /// `Reference Directory` is where the 100 step JDs are stored.
    /// `Reference Items` is a list of all JDs
    array: GenericSegment,

    /// spg4 uses a geopotential model which is loaded from the spice kernel.
    /// Unfortunately SGP4 doesn't support custom altitude bounds, but this
    /// probably shouldn't be altered from the defaults.
    geopotential: Geopotential,
}

impl SpkSegmentType10 {
    #[inline(always)]
    fn get_times(&self) -> &[f64] {
        self.array.get_reference_items()
    }

    /// Return the SGP4 record stored within the spice kernel.
    #[inline(always)]
    fn get_record(&self, idx: usize) -> Constants {
        let rec = self.array.get_packet::<15>(idx);
        let [_, _, _, b_star, inclination, right_ascension, eccentricity, argument_of_perigee, mean_anomaly, kozai_mean_motion, epoch, _, _, _, _] =
            *rec;

        let epoch = julian_years_since_j2000(
            &Time::<TDB>::new(spice_jds_to_jd(epoch))
                .utc()
                .to_datetime()
                .unwrap()
                .naive_utc(),
        );

        // use the provided goepotential even if it is not correct.
        let orbit_0 = Orbit::from_kozai_elements(
            &self.geopotential,
            inclination,
            right_ascension,
            eccentricity,
            argument_of_perigee,
            mean_anomaly,
            kozai_mean_motion,
        )
        .expect("Failed to load orbit values");
        Constants::new(
            self.geopotential,
            sgp4::iau_epoch_to_sidereal_time,
            epoch,
            b_star,
            orbit_0,
        )
        .expect("Failed to load orbit values")
    }

    #[inline(always)]
    fn try_get_pos_vel<T: InertialFrame>(
        &self,
        _: &SpkSegment<T>,
        jd: f64,
    ) -> KeteResult<([f64; 3], [f64; 3])> {
        // TODO: this does not yet implement the interpolation between two neighboring states
        // which is present in the cSPICE implementation.
        // This currently matches the cspice implementation to within about 20km, where the error
        // is less near the year 2000.

        // There is also an outstanding small time conversion issue.
        // I am somewhat certain that this conversion is incorrect in cSPICE itself.
        // Much of this error may be fixed by applying a small linear offset to time which
        // causes about a 3 second offset in 2024 vs a 0 second offset in 2000.
        // See #66 for more details.
        let jd = jd_to_spice_jd(jd);
        let times = self.get_times();
        let idx: usize = match times.binary_search_by(|probe| probe.total_cmp(&jd)) {
            Ok(c) => c,
            Err(c) => {
                if c == 0 {
                    c
                } else if c == times.len() || (jd - times[c - 1]).abs() < (jd - times[c]).abs() {
                    c - 1
                } else {
                    c
                }
            }
        };
        let epoch = times[idx];
        let record = self.get_record(idx);
        let prediction = record
            .propagate(MinutesSinceEpoch((jd - epoch) / 60.0))
            .unwrap();

        let [x, y, z] = prediction.position;
        let [vx, vy, vz] = prediction.velocity;
        let v_scale = 86400.0 / AU_KM;
        Ok((
            [x / AU_KM, y / AU_KM, z / AU_KM],
            [vx * v_scale, vy * v_scale, vz * v_scale],
        ))
    }
}

impl From<DafArray> for SpkSegmentType10 {
    fn from(array: DafArray) -> Self {
        let array: GenericSegment = array.try_into().unwrap();
        let constants = array.constants();
        let geopotential = Geopotential {
            j2: constants[0],
            j3: constants[1],
            j4: constants[2],
            ke: constants[3],
            ae: constants[6],
        };

        SpkSegmentType10 {
            array,
            geopotential,
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
    #[inline(always)]
    fn get_record(&self, idx: usize) -> Type13RecordView {
        unsafe {
            let rec = self.array.data.get_unchecked(idx * 6..(idx + 1) * 6);
            Type13RecordView {
                pos: rec[0..3].try_into().unwrap(),
                vel: rec[3..6].try_into().unwrap(),
            }
        }
    }

    #[inline(always)]
    fn get_times(&self) -> &[f64] {
        unsafe {
            self.array
                .data
                .get_unchecked(self.n_records * 6..self.n_records * 7)
        }
    }

    #[inline(always)]
    fn try_get_pos_vel<T: InertialFrame>(
        &self,
        _: &SpkSegment<T>,
        jd: f64,
    ) -> KeteResult<([f64; 3], [f64; 3])> {
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

/// Type 18 Record
///
/// This is actually 2 types in 1, under the stated goal of reducing the number
/// of unique SPICE kernel types.
///
/// Subtype 0 is a Hermite Interpolation of both position and velocity, a record
/// contains 12 numbers, 3 position, 3 derivative of position, 3 velocity, and 3
/// derivative of velocity. Note that it explicitly allows that the 3 velocity
/// values do not have to match the derivative of the position vectors.
/// Subtype 1 is a Lagrange Interpolation, a record of which contains 6 values,
/// 3 position, and 3 velocity.
///
/// <https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/FORTRAN/req/spk.html#Type%2018:%20ESOC/DDID%20Hermite/Lagrange%20Interpolation>
#[derive(Debug)]
pub struct SpkSegmentType18 {
    array: DafArray,
    subtype: usize,
    window_size: usize,
    n_records: usize,
    record_size: usize,
}

impl From<DafArray> for SpkSegmentType18 {
    fn from(array: DafArray) -> Self {
        let n_records = array[array.len() - 1] as usize;
        let window_size = array[array.len() - 2] as usize;
        let subtype = array[array.len() - 3] as usize;
        let record_size = {
            if subtype == 0 {
                12
            } else if subtype == 1 {
                6
            } else {
                panic!("Spk Segment Type 12 only supports subtype of 0 or 1")
            }
        };

        SpkSegmentType18 {
            array,
            subtype,
            window_size,
            n_records,
            record_size,
        }
    }
}

/// Type 18 Record View
/// A view into a record of type 18, provided mainly for clarity to the underlying
/// data structure.
struct Type18RecordView<'a> {
    pos: &'a [f64],
    vel: &'a [f64],
}

impl SpkSegmentType18 {
    #[inline(always)]
    fn get_record(&self, idx: usize) -> Type18RecordView {
        unsafe {
            let rec = self
                .array
                .data
                .get_unchecked(idx * self.record_size..(idx + 1) * self.record_size);
            Type18RecordView {
                pos: &rec[0..self.record_size / 2],
                vel: &rec[self.record_size / 2..self.record_size],
            }
        }
    }

    #[inline(always)]
    fn get_times(&self) -> &[f64] {
        unsafe {
            self.array.data.get_unchecked(
                self.n_records * self.record_size..self.n_records * (self.record_size + 1),
            )
        }
    }

    #[inline(always)]
    fn try_get_pos_vel<T: InertialFrame>(&self, _: &SpkSegment<T>, jd: f64) -> KeteResult<([f64; 3], [f64; 3])> {
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
        match self.subtype {
            0 => {
                for idx in 0..3 {
                    let p: Box<[f64]> = (0..self.window_size)
                        .map(|i| self.get_record(i + start_idx).pos[idx])
                        .collect();
                    let dp: Box<[f64]> = (0..self.window_size)
                        .map(|i| self.get_record(i + start_idx).pos[idx + 3])
                        .collect();
                    let (p, _) = hermite_interpolation(
                        &times[start_idx..start_idx + self.window_size],
                        &p,
                        &dp,
                        jd,
                    );
                    pos[idx] = p / AU_KM;

                    let p: Box<[f64]> = (0..self.window_size)
                        .map(|i| self.get_record(i + start_idx).vel[idx])
                        .collect();
                    let dp: Box<[f64]> = (0..self.window_size)
                        .map(|i| self.get_record(i + start_idx).vel[idx + 3])
                        .collect();
                    let (v, _) = hermite_interpolation(
                        &times[start_idx..start_idx + self.window_size],
                        &p,
                        &dp,
                        jd,
                    );
                    vel[idx] = v / AU_KM * 86400.;
                }
            }
            1 => {
                for idx in 0..3 {
                    let mut p: Box<[f64]> = (0..self.window_size)
                        .map(|i| self.get_record(i + start_idx).pos[idx])
                        .collect();
                    let mut dp: Box<[f64]> = (0..self.window_size)
                        .map(|i| self.get_record(i + start_idx).vel[idx])
                        .collect();
                    let p = lagrange_interpolation(
                        &times[start_idx..start_idx + self.window_size],
                        &mut p,
                        jd,
                    );
                    let v = lagrange_interpolation(
                        &times[start_idx..start_idx + self.window_size],
                        &mut dp,
                        jd,
                    );
                    pos[idx] = p / AU_KM;
                    vel[idx] = v / AU_KM * 86400.;
                }
            }
            _ => {
                unreachable!()
            }
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
    #[inline(always)]
    fn get_record(&self, idx: usize) -> &[f64] {
        unsafe {
            self.array
                .data
                .get_unchecked(idx * self.record_len..(idx + 1) * self.record_len)
        }
    }

    #[inline(always)]
    fn get_times(&self) -> &[f64] {
        unsafe {
            self.array.data.get_unchecked(
                self.n_records * self.record_len..self.n_records * (self.record_len + 1),
            )
        }
    }

    #[inline(always)]
    fn try_get_pos_vel<T: InertialFrame>(
        &self,
        _: &SpkSegment<T>,
        jd: f64,
    ) -> KeteResult<([f64; 3], [f64; 3])> {
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
                return Err(Error::IOError(
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

/// Segments of type 10 and 14 use a "generic segment" definition.
/// This segment type has poor documentation on the NAIF website, a significant amount
/// of reverse engineering was required to understand this.
/// The DAF Array is big flat vector of floats.

#[derive(Debug)]
#[allow(dead_code)]
struct GenericSegment {
    array: DafArray,

    /// Number of metadata value stored in this segment.
    n_meta: usize,

    // Below meta data is guaranteed to exist.
    /// address of the constant values
    const_addr: usize,

    /// Number of constants
    n_consts: usize,

    /// Address of reference directory
    ref_dir_addr: usize,

    /// Number of reference directory items
    n_item_ref_dir: usize,

    /// Type of reference directory
    ref_dir_type: usize,

    /// Address of reference items
    ref_items_addr: usize,

    /// Number of reference items
    n_ref_items: usize,

    /// Address of the data packets
    packet_dir_addr: usize,

    /// Number of data packets
    n_dir_packets: usize,

    /// Packet directory type
    packet_dir_dype: usize,

    /// Packet address
    packet_addr: usize,

    /// Number of data packets
    n_packets: usize,

    /// Address of reserved area
    res_addr: usize,

    /// number of entries in reserved area.
    n_reserved: usize,
}

#[allow(dead_code)]
impl GenericSegment {
    fn constants(&self) -> &[f64] {
        unsafe {
            self.array
                .data
                .get_unchecked(self.const_addr..self.const_addr + self.n_consts)
        }
    }

    fn peak_packet(&self) -> &[f64] {
        unsafe {
            self.array
                .data
                .get_unchecked(self.packet_addr - 3..self.packet_addr + 7)
        }
    }

    /// Slice into the entire reference items array.
    fn get_reference_items(&self) -> &[f64] {
        unsafe {
            self.array
                .data
                .get_unchecked(self.ref_items_addr..self.ref_items_addr + self.n_ref_items)
        }
    }

    fn get_packet<const T: usize>(&self, idx: usize) -> &[f64; T] {
        unsafe {
            self.array
                .data
                .get_unchecked(self.packet_addr + T * idx..self.packet_addr + T * (idx + 1))
                .try_into()
                .unwrap()
        }
    }
}

impl TryFrom<DafArray> for GenericSegment {
    type Error = Error;

    fn try_from(array: DafArray) -> KeteResult<Self> {
        // The very last value of this array is an int (cast to f64) which indicates the number
        // of meta-data values.

        let n_meta = array[array.len() - 1] as usize;

        if n_meta < 15 {
            Err(Error::IOError(
                "SPK File not correctly formatted. There are fewer values found than expected."
                    .into(),
            ))?;
        }
        // there are guaranteed to be 15 meta data values.
        let (
            const_addr,
            n_consts,
            ref_dir_addr,
            n_item_ref_dir,
            ref_dir_type,
            ref_items_addr,
            n_ref_items,
            packet_dir_addr,
            n_dir_packets,
            packet_dir_dype,
            packet_addr,
            n_packets,
        ) = array
            .data
            .get(array.len() - n_meta..array.len() - 1)
            .unwrap()
            .iter()
            .map(|x| *x as usize)
            .next_tuple()
            .unwrap();

        let (res_addr, n_reserved) = array
            .data
            .get(array.len() - n_meta..array.len() - 1)
            .unwrap()
            .iter()
            .map(|x| *x as usize)
            .next_tuple()
            .unwrap();

        Ok(GenericSegment {
            array,
            n_meta,
            const_addr,
            n_consts,
            ref_dir_addr,
            n_item_ref_dir,
            ref_dir_type,
            ref_items_addr,
            n_ref_items,
            packet_dir_addr,
            n_dir_packets,
            packet_dir_dype,
            packet_addr,
            n_packets,
            res_addr,
            n_reserved,
        })
    }
}
