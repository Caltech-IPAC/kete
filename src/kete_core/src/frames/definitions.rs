//! Define coordinate frames and object states.
//!
//! Distances measured in AU, time is in units of days with TDB scaling.
//!
use lazy_static::lazy_static;
use nalgebra::{Matrix3, Rotation3, UnitVector3, Vector3};
use serde::{Deserialize, Serialize};
use std::f64::consts::{PI, TAU};
use std::fmt::{self, Debug, Display};

use crate::prelude::{Error, NeosResult};
use crate::time::Time;

/// Coordinate frames.
///
/// Frames define the orientation of the coordinate system.
/// These distance units are always AU, and time is measured in days with TDB scaling.
///
#[derive(Debug, PartialEq, Clone, Copy, Deserialize, Serialize)]
pub enum Frame {
    /// Equivalent to Ecliptic J2000
    Ecliptic,

    /// Equivalent to Equatorial J2000
    Equatorial,

    /// FK4 Frame
    FK4,

    /// Galactic Frame as defined by SPICE
    Galactic,

    /// Unknown is to allow SPK files to be loaded with other frames.
    Unknown(usize),

    /// Non-inertial frame as defined by rotations from the ecliptic frame.
    ///
    /// - isize value represents the frame identifier.
    /// - array of 6 floats represent the euler angles and their derivatives to move to
    ///   this frame from the ecliptic frame.
    /// Rotation is done with a ZXZ set of chained rotations.
    EclipticNonInertial(isize, [f64; 6]),
    // Other non inertial frames will require multi-step conversions
}

impl Frame {
    /// Change a vector from the current frame into the target frame.
    pub fn try_vec_frame_change(
        &self,
        mut vec: Vector3<f64>,
        target: Frame,
    ) -> NeosResult<Vector3<f64>> {
        match self {
            Frame::Equatorial => {
                vec = equatorial_to_ecliptic(&vec);
            }
            Frame::Ecliptic => {}
            Frame::EclipticNonInertial(_, _) => {
                 Err(Error::ValueError(
                    "Cannot convert bare vector between non-inertial frames.  This may only be done with a State".into(),
                ))?
            }
            Frame::FK4 => {
                vec = fk4_to_ecliptic(&vec);
            }
            Frame::Galactic => {
                vec = galactic_to_ecliptic(&vec);
            }
            Frame::Unknown(id) =>  Err(Error::UnknownFrame(*id))?,
        }

        // new_vec is now in ecliptic.
        match target {
            Frame::Equatorial => {
                vec = ecliptic_to_equatorial(&vec);
            }
            Frame::Ecliptic => {}
            Frame::EclipticNonInertial(_, _) => {
                 Err(Error::ValueError(
                    "Cannot convert bare vector between non-inertial frames. This may only be done with a State.".into(),
                ))?
            }
            Frame::FK4 => {
                vec = ecliptic_to_fk4(&vec);
            }
            Frame::Galactic => {
                vec = ecliptic_to_galactic(&vec);
            }
            Frame::Unknown(id) =>  Err(Error::UnknownFrame(id))?,
        }
        Ok(vec)
    }
}

impl Display for Frame {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        Debug::fmt(self, f)
    }
}

/// Ecliptic obliquity angle in radians at the J2000 epoch. This is using the definition
/// from the 1984 JPL DE Series. These constants allow the conversion between Ecliptic
/// and Equatorial frames. Note that there are more modern definitions for these values,
/// however these are used for compatibility with JPL Horizons and Spice.
///
/// See:
///     - https://en.wikipedia.org/wiki/Axial_tilt#Short_term
///     - https://ssd.jpl.nasa.gov/horizons/manual.html#defs
const OBLIQUITY: f64 = 0.40909280422232897;

/// Compute the angle of obliquity of Earth.
///
/// This is only valid for several centuries near J2000.
///
/// The equation here is from the 2010 Astronomical Almanac.
///
pub fn calc_obliquity(jd: f64) -> f64 {
    // centuries from j2000
    let c = (jd - Time::j2000().jd) / 365.25 / 100.0;
    (23.439279444444444
        + c * (-0.013010213611111
            + c * (-5.08611111111111e-08
                + c * (5.565e-07 - c * (1.6e-10 + -1.1777777777777779e-11 * c)))))
        .to_radians()
}

lazy_static! {
    static ref ECLIPTIC_EQUATORIAL_ROT: Rotation3<f64> = {
        let x = nalgebra::Unit::new_unchecked(Vector3::new(1.0, 0.0, 0.0));
        Rotation3::from_axis_angle(&x, OBLIQUITY)
    };
    static ref FK4_ECLIPTIC_ROT: Rotation3<f64> = {
        let y = nalgebra::Unit::new_unchecked(Vector3::new(0.0, 1.0, 0.0));
        let z = nalgebra::Unit::new_unchecked(Vector3::new(0.0, 0.0, 1.0));
        let r1 = Rotation3::from_axis_angle(&z, (1152.84248596724 + 0.525) / 3600.0 * PI / 180.0);
        let r2 = Rotation3::from_axis_angle(&y, -1002.26108439117 / 3600.0 * PI / 180.0);
        let r3 = Rotation3::from_axis_angle(&z, 1153.04066200330 / 3600.0 * PI / 180.0);
        (*ECLIPTIC_EQUATORIAL_ROT).inverse() * r3 * r2 * r1
    };
    static ref GALACTIC_ECLIPTIC_ROT: Rotation3<f64> = {
        let x = nalgebra::Unit::new_unchecked(Vector3::new(1.0, 0.0, 0.0));
        let z = nalgebra::Unit::new_unchecked(Vector3::new(0.0, 0.0, 1.0));
        let r1 = Rotation3::from_axis_angle(&z, 1177200.0 / 3600.0 * PI / 180.0);
        let r2 = Rotation3::from_axis_angle(&x, 225360.0 / 3600.0 * PI / 180.0);
        let r3 = Rotation3::from_axis_angle(&z, 1016100.0 / 3600.0 * PI / 180.0);
        (*FK4_ECLIPTIC_ROT) * r3 * r2 * r1
    };
}

/// Convert a vector in the fk4 frame to ecliptic frame.
///
/// # Arguments
///
/// * `vec` - Vector, arbitrary units.
pub fn fk4_to_ecliptic(vec: &Vector3<f64>) -> Vector3<f64> {
    FK4_ECLIPTIC_ROT.transform_vector(vec)
}

/// Convert a vector in the ecliptic to the fk4 frame.
///
/// # Arguments
///
/// * `vec` - Vector, arbitrary units.
pub fn ecliptic_to_fk4(vec: &Vector3<f64>) -> Vector3<f64> {
    FK4_ECLIPTIC_ROT.inverse().transform_vector(vec)
}

/// Convert a vector in the galactic frame to ecliptic frame.
///
/// # Arguments
///
/// * `vec` - Vector, arbitrary units.
pub fn galactic_to_ecliptic(vec: &Vector3<f64>) -> Vector3<f64> {
    GALACTIC_ECLIPTIC_ROT.transform_vector(vec)
}

/// Convert a vector in the ecliptic to the galactic frame.
///
/// # Arguments
///
/// * `vec` - Vector, arbitrary units.
pub fn ecliptic_to_galactic(vec: &Vector3<f64>) -> Vector3<f64> {
    GALACTIC_ECLIPTIC_ROT.inverse().transform_vector(vec)
}

/// Convert a vector from ecliptic to equatorial J2000 frames assuming the 1984 JPL DE
/// definition.
///
/// # Arguments
///
/// * `vec` - Vector, arbitrary units.
pub fn ecliptic_to_equatorial(vec: &Vector3<f64>) -> Vector3<f64> {
    ECLIPTIC_EQUATORIAL_ROT.transform_vector(vec)
}

/// Convert a vector from equatorial to ecliptic J2000 frames assuming the 1984 JPL DE
/// definition.
///
/// # Arguments
///
/// * `vec` - Vector, arbitrary units.
pub fn equatorial_to_ecliptic(vec: &Vector3<f64>) -> Vector3<f64> {
    ECLIPTIC_EQUATORIAL_ROT.inverse().transform_vector(vec)
}

/// Derivative of the z rotation matrix with respect to the rotation angle.
fn rot_z_der(angle: f64) -> Matrix3<f64> {
    let (sin_a, cos_a) = angle.sin_cos();
    Matrix3::<f64>::from([[-sin_a, cos_a, 0.0], [-cos_a, -sin_a, 0.0], [0.0, 0.0, 0.0]])
}

/// Derivative of the x rotation matrix with respect to the rotation angle.
fn rot_x_der(angle: f64) -> Matrix3<f64> {
    let (sin_a, cos_a) = angle.sin_cos();
    Matrix3::<f64>::from([[0.0, 0.0, 0.0], [0.0, -sin_a, cos_a], [0.0, -cos_a, -sin_a]])
}

/// Compute two rotation matrices from a target inertial frame to the frame defined by
/// the provided angles. The first 3 angles here define the rotation ZXZ, the second
/// three values define the derivative of the 3 angles.
///
/// This then provides two rotation matrices, one is the 3x3 rotation matrix, and the
/// second is the derivative of the 3x3 matrix with respect to time. These two matrices
/// may be used to compute the new position and velocities when moving from one frame
/// to another.
pub fn noninertial_rotation(frame_angles: &[f64; 6]) -> (Matrix3<f64>, Matrix3<f64>) {
    let r_z1 = Rotation3::from_axis_angle(&Vector3::z_axis(), frame_angles[0]);
    let r_x = Rotation3::from_axis_angle(&Vector3::x_axis(), frame_angles[1]);
    let r_z0 = Rotation3::from_axis_angle(&Vector3::z_axis(), frame_angles[2]);
    let dr_z1 = rot_z_der(frame_angles[0]);
    let dr_x = rot_x_der(frame_angles[1]);
    let dr_z0 = rot_z_der(frame_angles[2]);

    // math for computing the derivative:
    // r = rot_z(z1) * rot_x(x) * rot_z(z0)
    // dr / dt =
    //  (d rot_z(z1) / d z1 * d z1 / dt) * rot_x(x) * rot_z(z0) +
    //  rot_z(z1) * (d rot_x(x) / d x * d x / dt) * rot_z(z0) +
    //  rot_z(z0) * rot_x(x) * (d rot_z(z1) / d z1 * d z1 / dt)
    let mut dr_dt = dr_z1 * r_x * r_z0 * frame_angles[3];
    dr_dt += r_z1 * dr_x * r_z0 * frame_angles[4];
    dr_dt += r_z1 * r_x * dr_z0 * frame_angles[5];

    ((r_z1 * r_x * r_z0).into(), dr_dt)
}

/// Convert positions and velocities from a non-inertial frame to an inertial frame
/// defined by the provided angles. The first 3 angles here define the rotation ZXZ, the
/// second three values define the derivative of the 3 angles. These angles define the
/// rotation from the inertial to the non-inertial frame.
pub fn noninertial_to_inertial(
    frame_angles: &[f64; 6],
    pos: &Vector3<f64>,
    vel: &Vector3<f64>,
) -> (Vector3<f64>, Vector3<f64>) {
    let (rot_p, rot_dp) = noninertial_rotation(frame_angles);

    let new_pos = rot_p * pos;
    let new_vel = rot_dp * pos + rot_p * vel;

    (new_pos, new_vel)
}

/// Convert positions and velocities from a non-inertial frame to an inertial frame
/// defined by the provided angles. The first 3 angles here define the rotation ZXZ, the
/// second three values define the derivative of the 3 angles. These angles define the
/// rotation from the inertial to the non-inertial frame.
pub fn inertial_to_noninertial(
    frame_angles: &[f64; 6],
    pos: &Vector3<f64>,
    vel: &Vector3<f64>,
) -> (Vector3<f64>, Vector3<f64>) {
    let (rot_p, rot_dp) = noninertial_rotation(frame_angles);

    let new_pos = rot_p.transpose() * pos;
    let new_vel = rot_dp.transpose() * pos + rot_p.transpose() * vel;

    (new_pos, new_vel)
}

/// Create a unit vector from polar spherical theta and phi angles in radians.
///
/// <https://en.wikipedia.org/wiki/Spherical_coordinate_system#Cartesian_coordinates>
pub fn from_polar_spherical(theta: f64, phi: f64) -> UnitVector3<f64> {
    let (theta_sin, theta_cos) = theta.sin_cos();
    let (phi_sin, phi_cos) = phi.sin_cos();
    UnitVector3::<f64>::new_unchecked([theta_sin * phi_cos, theta_sin * phi_sin, theta_cos].into())
}

/// Create a unit vector from latitude and longitude in units of radians.
pub fn from_lat_lon(lat: f64, lon: f64) -> UnitVector3<f64> {
    from_polar_spherical(PI / 2.0 - lat, lon)
}

/// Create a unit vector from ra and dec in units of radians.
pub fn from_ra_dec(ra: f64, dec: f64) -> UnitVector3<f64> {
    from_polar_spherical(PI / 2.0 - dec, ra)
}

/// Convert a unit vector to polar spherical coordinates.
///
/// <https://en.wikipedia.org/wiki/Spherical_coordinate_system#Cartesian_coordinates>
pub fn to_polar_spherical(vec: UnitVector3<f64>) -> (f64, f64) {
    let theta = vec.z.acos();
    let phi = vec.y.atan2(vec.x) % TAU;
    (theta, phi)
}

/// Convert a unit vector to latitude and longitude.
///
/// Input vector needs to be in the [`Frame::Ecliptic`] frame.
pub fn to_lat_lon(vec: UnitVector3<f64>) -> (f64, f64) {
    let (mut lat, mut lon) = to_polar_spherical(vec);
    if lat > PI {
        lat = TAU - lat;
        lon += PI
    }
    (PI / 2.0 - lat, lon)
}

/// Convert a unit vector to ra and dec.
///
/// Input vector needs to be in the [`Frame::Equatorial`] frame.
pub fn to_ra_dec(vec: UnitVector3<f64>) -> (f64, f64) {
    let (mut dec, mut ra) = to_polar_spherical(vec);
    if dec > PI {
        dec = TAU - dec;
        ra += PI
    }
    (ra, PI / 2.0 - dec)
}

/// Rotate a collection of vectors around the specified rotation vector.
///
/// # Arguments
///
/// * `vectors` - A Matrix containing N vectors of length 3.
/// * `rotation_vec` - The single vector around which to rotate the vectors.
/// * `angle` - The angle in radians to rotate the vectors.
///
pub fn rotate_around(
    vector: &Vector3<f64>,
    rotation_vec: Vector3<f64>,
    angle: f64,
) -> Vector3<f64> {
    let rot = Rotation3::from_axis_angle(&UnitVector3::new_normalize(rotation_vec), angle);
    rot.transform_vector(vector)
}

#[cfg(test)]
mod tests {
    use crate::frames::*;

    #[test]
    fn test_ecliptic_rot_roundtrip() {
        let vec = ecliptic_to_equatorial(&[1.0, 2.0, 3.0].into());
        let vec_return = equatorial_to_ecliptic(&vec);
        assert!((1.0 - vec_return[0]).abs() <= 10.0 * f64::EPSILON);
        assert!((2.0 - vec_return[1]).abs() <= 10.0 * f64::EPSILON);
        assert!((3.0 - vec_return[2]).abs() <= 10.0 * f64::EPSILON);
    }
    #[test]
    fn test_fk4_roundtrip() {
        let vec = ecliptic_to_fk4(&[1.0, 2.0, 3.0].into());
        let vec_return = fk4_to_ecliptic(&vec);
        assert!((1.0 - vec_return[0]).abs() <= 10.0 * f64::EPSILON);
        assert!((2.0 - vec_return[1]).abs() <= 10.0 * f64::EPSILON);
        assert!((3.0 - vec_return[2]).abs() <= 10.0 * f64::EPSILON);
    }
    #[test]
    fn test_galactic_rot_roundtrip() {
        let vec = ecliptic_to_galactic(&[1.0, 2.0, 3.0].into());
        let vec_return = galactic_to_ecliptic(&vec);
        assert!((1.0 - vec_return[0]).abs() <= 10.0 * f64::EPSILON);
        assert!((2.0 - vec_return[1]).abs() <= 10.0 * f64::EPSILON);
        assert!((3.0 - vec_return[2]).abs() <= 10.0 * f64::EPSILON);
    }

    #[test]
    fn test_noninertial_rot_roundtrip() {
        let angles = [0.11, 0.21, 0.31, 0.41, 0.51, 0.61];
        let pos = [1.0, 2.0, 3.0].into();
        let vel = [0.1, 0.2, 0.3].into();
        let (r_pos, r_vel) = inertial_to_noninertial(&angles, &pos, &vel);
        let (pos_return, vel_return) = noninertial_to_inertial(&angles, &r_pos, &r_vel);

        assert!((1.0 - pos_return[0]).abs() <= 10.0 * f64::EPSILON);
        assert!((2.0 - pos_return[1]).abs() <= 10.0 * f64::EPSILON);
        assert!((3.0 - pos_return[2]).abs() <= 10.0 * f64::EPSILON);
        assert!((0.1 - vel_return[0]).abs() <= 10.0 * f64::EPSILON);
        assert!((0.2 - vel_return[1]).abs() <= 10.0 * f64::EPSILON);
        assert!((0.3 - vel_return[2]).abs() <= 10.0 * f64::EPSILON);
    }
}
