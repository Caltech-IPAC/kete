//! Basic Geometric shapes on a sphere (typically the celestial sphere)

use crate::frames::{rotate_around, Frame};
use crate::io::serde_const_arr;
use nalgebra::{UnitVector3, Vector3};
use serde::{Deserialize, Serialize};
use std::f64::consts::FRAC_PI_2;

use super::NeosResult;

/// Bounded areas can either contains a vector or not.
/// This enum specifies if the vector is within the area, or
/// the minimum distance the vector must move to be within the area.
#[derive(Debug, Clone)]
pub enum Contains {
    /// Vector is contained within the area.
    Inside,

    /// Vector is outside of the area
    /// The f64 defines the minimum distance required to move into the area.
    Outside(f64),
}

impl Contains {
    /// Returns true if is_inside.
    pub fn is_inside(&self) -> bool {
        matches!(self, Contains::Inside)
    }
}

/// Given an iterable of [`Contains`], find the closest one to being Inside.
pub(super) fn closest_inside(contains: &[Contains]) -> (usize, Contains) {
    let mut best = (usize::MAX, f64::INFINITY);
    for (idx, con) in contains.iter().enumerate() {
        match con {
            Contains::Inside => return (idx, Contains::Inside),
            Contains::Outside(d) => {
                if d < &best.1 {
                    best = (idx, *d)
                }
            }
        }
    }
    (best.0, Contains::Outside(best.1))
}

/// Trait which defines an area on the surface area of a sphere.
pub trait SkyPatch: Sized {
    /// Checks to see if a unit vector is within the bounded area.
    fn contains(&self, obs_to_obj: &Vector3<f64>) -> Contains;

    /// Frame of reference for the bounded area
    fn frame(&self) -> Frame;

    /// Change the frame of the bounded area to the new target frame.
    fn try_frame_change(&self, target_frame: Frame) -> NeosResult<Self>;

    /// Center of the field of view
    fn pointing(&self) -> UnitVector3<f64>;
}

/// A Spherical Polygon as represented by a series of planes through the central axis.
///
/// Conceptually the spherical polygon as defined here can be thought of as a unit sphere
/// with a polygon drawn on it. This polygon has sides defined by great circle lines.
/// If these great circle lines were continued all the way around the sphere, they would
/// define planes which cut the circle. Defining surface normals to each of these planes
/// enables us to define an "inside" and "outside" of the polygon, where inside is
/// defined by having all of the unit vectors pointing toward the inside of the polygon.
///
/// To test if a point on the unit circle is inside of the polygon, all that is
/// required is it multiply each surface normal against the vector which defines the
/// point's position. If all of the dot products are positive, then the vector must lie
/// within the polygon. If any are negative, then the point must be outside of the
/// polygon.
///
/// Note that the surface normals are defined by planes through the center of the sphere
/// meaning that constructed polygons can only have great circle edges.
///
#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct SphericalPolygon<const N_SIDES: usize> {
    /// Normal vectors which define the boundary of a polygon.
    #[serde(with = "serde_const_arr")]
    edge_normals: [[f64; 3]; N_SIDES],

    /// Coordinate frame where the boundary is defined.
    pub frame: Frame,
}

/// A rectangular patch of sky.
pub type OnSkyRectangle = SphericalPolygon<4>;

impl OnSkyRectangle {
    /// Construct a rectangular spherical polygon.
    ///
    /// # Arguments
    ///
    /// * `edge_normals` - Normal vectors which define the boundary of a polygon.
    /// * `frame` - Coordinate frame of the rectangle.
    pub fn from_normals(edge_normals: [[f64; 3]; 4], frame: Frame) -> Self {
        // construct the 4 normal vectors
        Self {
            edge_normals,
            frame,
        }
    }

    /// Construct a rectangular spherical polygon.
    ///
    /// This constructs a new SphericalPolygon made up of a rectangular shape on the unit
    /// sphere. Where the edges of the rectangle are great circle arcs.
    ///
    /// # Arguments
    ///
    /// * `pointing` - A vector pointing to the center of the rectangle.
    /// * `rotation` - Rotation of the center of the rectangle in radians.
    /// * `lon_width` - If the rotation is 0, this defines the width of the rectangle
    ///                 longitudinally in radians.
    /// * `lat_width` - If the rotation is 0, this defines the width of the rectangle
    ///                 latitudinally in radians.
    /// * `frame` - Coordinate frame of the rectangle.
    pub fn new(
        pointing: Vector3<f64>,
        rotation: f64,
        lon_width: f64,
        lat_width: f64,
        frame: Frame,
    ) -> Self {
        // Rotate the Z axis to match the defined rotation angle, this vector is not
        // orthogonal to the pointing vector, but is in the correct plane of the final
        // up vector.
        let up_vec = rotate_around(&Vector3::new(0.0, 0.0, 1.0), pointing, -rotation);

        // construct the vector orthogonal to the pointing and rotate z axis vectors.
        // left = cross(up, pointing)
        let left_vec = pointing.cross(&up_vec);

        // Given the new left vector, and the existing orthogonal pointing vector,
        // construct a new up vector which is in the same plane as it was before, but now
        // orthogonal to the two existing vectors.
        // up = cross(pointing, left)
        let up_vec = pointing.cross(&left_vec);

        // These have to be enumerated in clockwise order for the pointing calculation to be correct.
        let n1: Vector3<f64> = rotate_around(&left_vec, up_vec, -lon_width / 2.0);
        let n2: Vector3<f64> = rotate_around(&up_vec, left_vec, lat_width / 2.0);
        let n3: Vector3<f64> = rotate_around(&(-left_vec), up_vec, lon_width / 2.0);
        let n4: Vector3<f64> = rotate_around(&(-up_vec), left_vec, -lat_width / 2.0);

        // construct the 4 normal vectors
        Self {
            edge_normals: [n1.into(), n2.into(), n3.into(), n4.into()],
            frame,
        }
    }

    /// Construct the patch from the 4 corners of the field of view.
    /// Corners have to be provided in the clockwise order around the center of the FOV.
    pub fn from_corners(corners: [Vector3<f64>; 4], frame: Frame) -> Self {
        let n1 = corners[0].cross(&corners[1]).normalize();
        let n2 = corners[1].cross(&corners[2]).normalize();
        let n3 = corners[2].cross(&corners[3]).normalize();
        let n4 = corners[3].cross(&corners[0]).normalize();

        Self {
            edge_normals: [n1.into(), n2.into(), n3.into(), n4.into()],
            frame,
        }
    }

    /// Latitudinal width of the patch, the assumes the patch is rectangular.
    pub fn lat_width(&self) -> f64 {
        let pointing = self.pointing();
        2.0 * (FRAC_PI_2 - pointing.angle(&Vector3::from(self.edge_normals[1])))
    }

    /// Longitudinal width of the patch, the assumes the patch is rectangular.
    pub fn lon_width(&self) -> f64 {
        let pointing = self.pointing();
        2.0 * (FRAC_PI_2 - pointing.angle(&Vector3::from(self.edge_normals[0])))
    }
}

impl<const D: usize> SkyPatch for SphericalPolygon<D> {
    /// Is the obs_to_obj vector inside of the polygon.
    fn contains(&self, obs_to_obj: &Vector3<f64>) -> Contains {
        let mut closest_edge = f64::NEG_INFINITY;
        for normal in self.edge_normals.iter() {
            let normal = Vector3::from(*normal);
            let d = obs_to_obj.dot(&normal);
            if d.is_nan() {
                return Contains::Outside(d);
            }
            // of all the edges where d is outside, we need to find the maximum
            // distance one. This will be the *minimum* distance which the
            // object must move to be inside of the patch.
            if d.is_sign_negative() && d.abs() > closest_edge {
                closest_edge = d.abs();
            }
        }

        match closest_edge {
            x if x.is_finite() => Contains::Outside(x.min(obs_to_obj.norm())),
            _ => Contains::Inside,
        }
    }

    fn frame(&self) -> Frame {
        self.frame
    }

    fn try_frame_change(&self, target_frame: Frame) -> NeosResult<Self> {
        let new_edges = self
            .edge_normals
            .iter()
            .map(|vec| {
                self.frame
                    .try_vec_frame_change(Vector3::from(*vec), target_frame)
            })
            .collect::<NeosResult<Vec<_>>>()?;
        let new_edges: Vec<[f64; 3]> = new_edges.into_iter().map(|e| e.into()).collect();
        let new_edges: [[f64; 3]; D] = new_edges.try_into().unwrap();

        Ok(Self {
            edge_normals: new_edges,
            frame: target_frame,
        })
    }

    fn pointing(&self) -> UnitVector3<f64> {
        let mut x = 0.0;
        let mut y = 0.0;
        let mut z = 0.0;

        for (idx, idy) in (0..D).zip(1..D) {
            let v = Vector3::from(self.edge_normals[idx])
                .cross(&Vector3::from(self.edge_normals[idy]))
                .normalize();

            x += v.x;
            y += v.y;
            z += v.z;
        }

        let v = Vector3::from(self.edge_normals[D - 1])
            .cross(&Vector3::from(self.edge_normals[0]))
            .normalize();
        x += v.x;
        y += v.y;
        z += v.z;

        UnitVector3::new_normalize([x, y, z].into())
    }
}

/// Represent a cone on a sphere.
#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct SphericalCone {
    /// Unit vector which defines the direction of the cone.
    pointing: [f64; 3],

    /// Size of the cone in degrees.
    pub angle: f64,

    frame: Frame,
}

impl SphericalCone {
    /// Construct a new `SphericalCone` given the central vector and the angle of the
    /// cone.
    pub fn new(pointing: &Vector3<f64>, angle: f64, frame: Frame) -> Self {
        let pointing = pointing.normalize();
        Self {
            pointing: pointing.into(),
            angle,
            frame,
        }
    }
}

impl SkyPatch for SphericalCone {
    /// Is the obs_to_obj vector inside of the cone.
    fn contains(&self, obs_to_obj: &Vector3<f64>) -> Contains {
        let dist = Vector3::from(self.pointing)
            .dot(&obs_to_obj.normalize())
            .acos()
            .abs();
        match dist {
            // if d is less than the angle, it is inside cone
            d if d <= self.angle => Contains::Inside,

            // outside of the cone, but how badly?
            d => {
                let r = obs_to_obj.norm();
                let min_angle = (d - self.angle).abs();

                // The minimum distance required for the object to be within the cone
                // of shame.
                match min_angle {
                    // if the min angle is more than 90 degrees, then the object
                    // can be visible by moving on top of the observer, so going r dist
                    theta if theta > FRAC_PI_2 => Contains::Outside(r),

                    // if the min angle is less than 90 degrees, then the object can
                    // move directly toward the edge of the cone, which is a right
                    // angle triangle, where the hypotenuse is r
                    theta => Contains::Outside(theta.sin() * r),
                }
            }
        }
    }

    fn frame(&self) -> Frame {
        self.frame
    }

    fn try_frame_change(&self, target_frame: Frame) -> NeosResult<Self> {
        let pointing = self
            .frame
            .try_vec_frame_change(Vector3::from(self.pointing), target_frame)?;

        Ok(Self {
            pointing: pointing.into(),
            angle: self.angle,
            frame: target_frame,
        })
    }

    fn pointing(&self) -> UnitVector3<f64> {
        UnitVector3::new_normalize(self.pointing.into())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_rectangular_patch() {
        let rot = (45f64).to_radians();
        let inside = [1.0, 0.01, 0.01].into();
        let outside = [1.0, 0.1, 0.0].into();
        let just_inside = [1.0, (0.05f64).sin() * 0.99, (0.05f64).sin() * 0.99].into();
        let just_outside = [1.0, (0.05f64).sin() * 1.01, (0.05f64).sin() * 1.01].into();
        let fov = OnSkyRectangle::new([1.0, 0.0, 0.0].into(), 0.0, 0.1, 0.1, Frame::Ecliptic);
        let fov_rot = OnSkyRectangle::new([1.0, 0.0, 0.0].into(), rot, 0.1, 0.1, Frame::Ecliptic);

        assert!(fov.contains(&inside).is_inside());
        assert!(fov.contains(&just_inside).is_inside());
        assert!(!fov.contains(&outside).is_inside());
        assert!(!fov.contains(&just_outside).is_inside());

        assert!(fov_rot.contains(&inside).is_inside());
        assert!(!fov_rot.contains(&just_inside).is_inside());
        assert!((fov_rot.pointing().into_inner() - Vector3::from([1.0, 0.0, 0.0])).norm() < 1e-10);
    }

    #[test]
    fn test_rectangular_patch_latlong() {
        let rot = (45f64).to_radians();
        let fov = OnSkyRectangle::new([1.0, 0.0, 0.0].into(), 0.0, 0.1, 0.2, Frame::Ecliptic);
        let fov_rot = OnSkyRectangle::new([1.0, 0.0, 0.0].into(), rot, 0.1, 0.2, Frame::Ecliptic);

        assert!((fov.lat_width() - 0.2).abs() < 1e-10);
        assert!((fov.lon_width() - 0.1).abs() < 1e-10);
        assert!((fov_rot.lat_width() - 0.2).abs() < 1e-10);
        assert!((fov_rot.lon_width() - 0.1).abs() < 1e-10);
    }
}
