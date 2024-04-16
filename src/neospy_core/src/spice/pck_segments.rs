/// Most users should interface with `pck.rs`, not this module.
///
/// PCK Files are collections of `Segments`, which are ranges of times where the state
/// of an object is recorded. These segments are typically made up of many individual
/// `Records`, with an associated maximum and minimum time where they are valid for.
///
/// There are unique structs for each possible segment type, not all are currently
/// supported. Each segment type must implement the PCKSegment trait, which allows for
/// the loading and querying of states contained within.
///
/// <https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/pck.html>
///
/// There is a lot of repetition in this file, as many of the segment types have very
/// similar internal structures.
///
use super::binary::read_f64_vec;
use super::interpolation::*;
use super::records::DafRecords;
use super::spice_jds_to_jd;
use crate::errors::NEOSpyError;
use crate::frames::Frame;
use std::fmt::Debug;
use std::io::{Read, Seek};

#[derive(Debug)]
enum PCKSegmentType {
    Type2(PckSegmentType2),
}

#[derive(Debug)]
pub struct PckSegment {
    /// The reference center NAIF ID for the central body in this Segment.
    pub center_id: isize,

    /// [`Frame`] of reference for this Segment.
    pub ref_frame: Frame,

    /// Start time of the segment.
    pub jd_start: f64,

    /// End time of the segment.
    pub jd_end: f64,

    /// Type of the segment
    pub segment_type: usize,

    /// Internal data representation.
    segment: PCKSegmentType,
}

impl PckSegment {
    pub fn contains(&self, jd: f64) -> bool {
        (jd >= self.jd_start) && (jd <= self.jd_end)
    }

    /// Load this segment from the provided file and summary.
    pub fn from_summary<T: Read + Seek>(
        file: &mut T,
        summary: (Box<[f64]>, Box<[i32]>),
    ) -> Result<PckSegment, NEOSpyError> {
        let (floats, ints) = summary;
        let jd_start = spice_jds_to_jd(floats[0]);
        let jd_end = spice_jds_to_jd(floats[1]);

        let center_id = ints[0] as isize;
        let frame_id = ints[1];
        let segment_type = ints[2] as usize;
        let array_start = ints[3] as usize;
        let array_end = ints[4] as usize;

        let ref_frame = match frame_id {
            1 => Frame::Equatorial, // J2000
            17 => Frame::Ecliptic,  // ECLIPJ2000
            _ => Frame::Unknown(frame_id as usize),
        };

        let segment = match segment_type {
            2 => PckSegmentType2::try_load(file, array_start, array_end)?,
            v => {
                return Err(NEOSpyError::IOError(format!(
                    "SPK Segment type {:?} not supported.",
                    v
                )));
            }
        };

        Ok(PckSegment {
            jd_start,
            jd_end,
            center_id,
            ref_frame,
            segment_type,
            segment,
        })
    }

    /// Return the [`Frame`] at the specified JD. If the requested time is not within
    /// the available range, this will fail.
    pub fn try_get_orientation(&self, jd: f64) -> Result<Frame, NEOSpyError> {
        if jd < self.jd_start || jd > self.jd_end {
            return Err(NEOSpyError::DAFLimits(
                "JD is not present in this record.".to_string(),
            ));
        }

        match &self.segment {
            PCKSegmentType::Type2(v) => v.try_get_orientation(self, jd),
        }
    }
}

/// Chebyshev polynomials (Euler angles only)
///
/// https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/spk.html#Type%201:%20Modified%20Difference%20Arrays
///
#[derive(Debug)]
pub struct PckSegmentType2 {
    pub jd_step: f64,
    pub n_coef: usize,
    pub records: DafRecords,
}

impl PckSegmentType2 {
    fn try_load<T: Read + Seek>(
        file: &mut T,
        array_start: usize,
        array_end: usize,
    ) -> Result<PCKSegmentType, NEOSpyError> {
        let _ = file.seek(std::io::SeekFrom::Start(8 * (array_start as u64 - 1)))?;
        let array_len = array_end - array_start + 1;
        let bin_segment = read_f64_vec(file, array_len, true)?;

        let n_rec = bin_segment[array_len - 1] as usize;
        let rec_len = bin_segment[array_len - 2] as usize;
        let jd_step = bin_segment[array_len - 3] / 86400.0;
        let n_coef = (rec_len - 2) / 3;

        let mut idy = 0;
        let mut records = DafRecords::with_capacity(n_rec, rec_len);
        for _ in 0..n_rec {
            let mut record = Vec::<f64>::with_capacity(rec_len);
            for idx in 0..rec_len {
                record.push({
                    if idx == 0 {
                        spice_jds_to_jd(bin_segment[idy])
                    } else if idx == 1 {
                        bin_segment[idy] / 86400.0
                    } else {
                        bin_segment[idy]
                    }
                });
                idy += 1;
            }
            records.try_push(record)?;
        }

        Ok(PCKSegmentType::Type2(PckSegmentType2 {
            jd_step,
            n_coef,
            records,
        }))
    }

    /// Return the stored orientation, along with the rate of change of the orientation.
    fn try_get_orientation(&self, segment: &PckSegment, jd: f64) -> Result<Frame, NEOSpyError> {
        // Records in the segment contain information about the central position of the
        // north pole, as well as the position of the prime meridian. These values for
        // type 2 segments are stored as chebyshev polynomials of the first kind, in
        // essentially the exact same format as the Type 2 SPK segments.
        // Records for this type are structured as so:
        // - time at midpoint of record.
        // - (length of time record is valid for) / 2.0
        // - N Chebyshev polynomial coefficients for ra
        // - N Chebyshev polynomial coefficients for dec
        // - N Chebyshev polynomial coefficients for w
        //
        // Rate of change for each of these values can be calculated by using the
        // derivative of chebyshev of the first kind, which is done below.
        let record_index = ((jd - segment.jd_start) / self.jd_step).floor() as usize;
        let record = &self.records[record_index];
        let t_mid = record[0];
        let t_step = record[1];
        let t = (jd - t_mid) / t_step;

        let ra_coef = &record[2..(self.n_coef + 2)];
        let dec_coef = &record[(self.n_coef + 2)..(2 * self.n_coef + 2)];
        let w_coef = &record[(2 * self.n_coef + 2)..(3 * self.n_coef + 2)];

        let (ra, ra_der) = chebyshev_evaluate_both(t, ra_coef)?;
        let (dec, dec_der) = chebyshev_evaluate_both(t, dec_coef)?;
        let (w, w_der) = chebyshev_evaluate_both(t, w_coef)?;

        // rem_euclid is equivalent to the modulo operator, so this maps w to [0, 2pi]
        let w = w.rem_euclid(std::f64::consts::TAU);

        Ok(Frame::EclipticNonInertial(
            segment.center_id,
            [
                ra,
                dec,
                w,
                ra_der / t_step,
                dec_der / t_step,
                w_der / t_step,
            ],
        ))
    }
}
