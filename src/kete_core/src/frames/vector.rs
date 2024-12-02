use super::{Ecliptic, Equatorial, InertialFrame};
use nalgebra::{Rotation3, UnitVector3, Vector3};
use serde::{Deserialize, Serialize};
use std::f64::consts::{PI, TAU};
use std::fmt::Debug;
use std::marker::PhantomData;
use std::ops::{Add, Div, DivAssign, Index, IndexMut, Mul, MulAssign, Neg, Sub};

/// Vector with frame information.
#[derive(Debug, Clone, Copy, Deserialize, Serialize, PartialEq)]
pub struct Vector<T: InertialFrame> {
    /// Raw Vector
    pub vec: [f64; 3],

    /// PhantomData is used here as the frame is only a record keeping convenience.
    frame: PhantomData<T>,
}

impl<T: InertialFrame> Vector<T> {
    /// New Vector
    #[inline(always)]
    pub fn new(vec: [f64; 3]) -> Self {
        Vector::<T> {
            vec,
            frame: PhantomData,
        }
    }

    /// New Vector of NANs
    #[inline(always)]
    pub fn new_nan() -> Self {
        Vector::<T> {
            vec: [f64::NAN, f64::NAN, f64::NAN],
            frame: PhantomData,
        }
    }

    /// Convert Vector from one frame to another.
    #[inline(always)]
    pub fn into_frame<Target: InertialFrame>(self) -> Vector<Target> {
        let vec = T::convert::<Target>(self.into());
        Vector::<Target>::new(vec.into())
    }

    /// Rotate a vector around the specified rotation vector.
    ///
    /// # Arguments
    ///
    /// * `rotation_vec` - The single vector around which to rotate the vectors.
    /// * `angle` - The angle in radians to rotate the vectors.
    ///
    #[inline(always)]
    pub fn rotate_around(self, rotation_vec: Vector<T>, angle: f64) -> Self {
        let rot =
            Rotation3::from_axis_angle(&UnitVector3::new_normalize(rotation_vec.into()), angle);
        rot.transform_vector(&self.into()).into()
    }

    /// Dot product between two vectors
    #[inline(always)]
    pub fn dot(&self, other: &Vector<T>) -> f64 {
        self.vec
            .iter()
            .zip(other.vec.iter())
            .map(|(a, b)| a * b)
            .sum()
    }

    /// Cross product between two vectors
    #[inline(always)]
    pub fn cross(&self, other: &Vector<T>) -> Vector<T> {
        Vector3::from(self.vec)
            .cross(&Vector3::from(other.vec))
            .into()
    }

    /// Squared euclidean length.
    #[inline(always)]
    pub fn norm_squared(&self) -> f64 {
        self.vec.iter().map(|a| a.powi(2)).sum()
    }

    /// The euclidean length of the vector.
    #[inline(always)]
    pub fn norm(&self) -> f64 {
        self.vec.iter().map(|a| a.powi(2)).sum::<f64>().sqrt()
    }

    /// The angle between two vectors in radians.
    #[inline(always)]
    pub fn angle(&self, other: &Self) -> f64 {
        Vector3::from(self.vec).angle(&Vector3::from(other.vec))
    }

    /// Create a new vector of unit length in the same direction as this vector.
    #[inline(always)]
    pub fn normalize(&self) -> Self {
        self / self.norm()
    }
    /// Normalize the vector to unit length in place.
    #[inline(always)]
    pub fn normalize_mut(&mut self) {
        *self /= self.norm();
    }

    /// Create a unit vector from polar spherical theta and phi angles in radians.
    ///
    /// <https://en.wikipedia.org/wiki/Spherical_coordinate_system#Cartesian_coordinates>
    #[inline(always)]
    pub fn from_polar_spherical(theta: f64, phi: f64) -> Self {
        let (theta_sin, theta_cos) = theta.sin_cos();
        let (phi_sin, phi_cos) = phi.sin_cos();
        [theta_sin * phi_cos, theta_sin * phi_sin, theta_cos].into()
    }

    /// Convert a unit vector to polar spherical coordinates.
    ///
    /// <https://en.wikipedia.org/wiki/Spherical_coordinate_system#Cartesian_coordinates>
    #[inline(always)]
    pub fn to_polar_spherical(&self) -> (f64, f64) {
        let theta = self.vec[2].acos();
        let phi = self.vec[1].atan2(self.vec[0]) % TAU;
        (theta, phi)
    }

    /// X-coordinate of the vector
    #[inline(always)]
    pub fn x(&self) -> f64 {
        unsafe { *self.vec.get_unchecked(0) }
    }

    /// Y-coordinate of the vector
    #[inline(always)]
    pub fn y(&self) -> f64 {
        unsafe { *self.vec.get_unchecked(1) }
    }

    /// Z-coordinate of the vector
    #[inline(always)]
    pub fn z(&self) -> f64 {
        unsafe { *self.vec.get_unchecked(2) }
    }

    /// Are all elements of the vector finite valued.
    #[inline(always)]
    pub fn is_finite(&self) -> bool {
        self.vec.iter().all(|x| x.is_finite())
    }
}

/// Unit length vector.
/// Thin wrapper over [`Vector`] to indicate it is of unit length.
#[repr(transparent)]
#[derive(Debug, Clone, Copy)]
pub struct UnitVector<T: InertialFrame>(Vector<T>);

impl<T: InertialFrame> UnitVector<T> {
    /// Make a new unit vector from a vector assumed to be unit length.
    #[inline(always)]
    pub fn new_unchecked(vec: Vector<T>) -> Self {
        Self(vec)
    }

    /// Make a new unit vector, normalizing it as necessary.
    #[inline(always)]
    pub fn new_checked(mut vec: Vector<T>) -> Self {
        vec.normalize_mut();
        Self(vec)
    }

    /// Unwrap the underlying vector.
    #[inline(always)]
    pub fn into_inner(self) -> Vector<T> {
        self.0
    }
}

impl Vector<Ecliptic> {
    /// Create a unit vector from latitude and longitude in units of radians.
    #[inline(always)]
    pub fn from_lat_lon(lat: f64, lon: f64) -> Self {
        Self::from_polar_spherical(PI / 2.0 - lat, lon)
    }

    /// Convert a unit vector to latitude and longitude.
    ///
    /// Input vector needs to be in the [`Frame::Ecliptic`] frame.
    #[inline(always)]
    pub fn to_lat_lon(self) -> (f64, f64) {
        let (mut lat, mut lon) = self.to_polar_spherical();
        if lat > PI {
            lat = TAU - lat;
            lon += PI
        }
        (PI / 2.0 - lat, lon)
    }
}

impl Vector<Equatorial> {
    /// Create a unit vector from ra and dec in units of radians.
    #[inline(always)]
    pub fn from_ra_dec(ra: f64, dec: f64) -> Self {
        Self::from_polar_spherical(PI / 2.0 - dec, ra)
    }

    /// Convert a unit vector to ra and dec.
    ///
    /// Input vector needs to be in the [`Frame::Equatorial`] frame.
    #[inline(always)]
    pub fn to_ra_dec(self) -> (f64, f64) {
        let (mut dec, mut ra) = self.to_polar_spherical();
        if dec > PI {
            dec = TAU - dec;
            ra += PI
        }
        (ra, PI / 2.0 - dec)
    }
}

impl<T: InertialFrame> Index<usize> for Vector<T> {
    type Output = f64;
    #[inline(always)]
    fn index(&self, index: usize) -> &Self::Output {
        &self.vec[index]
    }
}

impl<T: InertialFrame> IndexMut<usize> for Vector<T> {
    #[inline(always)]
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        &mut self.vec[index]
    }
}

impl<T: InertialFrame> IntoIterator for Vector<T> {
    type Item = f64;
    type IntoIter = std::array::IntoIter<Self::Item, 3>;
    #[inline(always)]
    fn into_iter(self) -> Self::IntoIter {
        self.vec.into_iter()
    }
}

impl<T: InertialFrame> From<[f64; 3]> for Vector<T> {
    #[inline(always)]
    fn from(value: [f64; 3]) -> Self {
        Vector::new(value)
    }
}

impl<T: InertialFrame> From<Vector3<f64>> for Vector<T> {
    #[inline(always)]
    fn from(value: Vector3<f64>) -> Self {
        Vector::new(value.into())
    }
}

impl<T: InertialFrame> From<Vector<T>> for Vector3<f64> {
    #[inline(always)]
    fn from(value: Vector<T>) -> Self {
        value.vec.into()
    }
}

impl<T: InertialFrame> From<Vector<T>> for [f64; 3] {
    #[inline(always)]
    fn from(value: Vector<T>) -> Self {
        value.vec
    }
}

impl<T: InertialFrame> From<Vector<T>> for Vec<f64> {
    #[inline(always)]
    fn from(value: Vector<T>) -> Self {
        value.vec.into()
    }
}

impl<T: InertialFrame> Sub<&Vector<T>> for &Vector<T> {
    type Output = Vector<T>;
    #[inline(always)]
    fn sub(self, rhs: &Vector<T>) -> Self::Output {
        Vector::<T>::new([
            self.vec[0] - rhs.vec[0],
            self.vec[1] - rhs.vec[1],
            self.vec[2] - rhs.vec[2],
        ])
    }
}

impl<T: InertialFrame> Sub<&Vector<T>> for Vector<T> {
    type Output = Vector<T>;
    #[inline(always)]
    fn sub(mut self, rhs: &Vector<T>) -> Self::Output {
        (0..3).for_each(|i| self.vec[i] -= rhs.vec[i]);
        self
    }
}

impl<T: InertialFrame> Sub<Vector<T>> for Vector<T> {
    type Output = Vector<T>;
    #[inline(always)]
    fn sub(mut self, rhs: Vector<T>) -> Self::Output {
        (0..3).for_each(|i| self.vec[i] -= rhs.vec[i]);
        self
    }
}

impl<T: InertialFrame> Add<&Vector<T>> for &Vector<T> {
    type Output = Vector<T>;
    #[inline(always)]
    fn add(self, rhs: &Vector<T>) -> Self::Output {
        Vector::<T>::new([
            self.vec[0] + rhs.vec[0],
            self.vec[1] + rhs.vec[1],
            self.vec[2] + rhs.vec[2],
        ])
    }
}

impl<T: InertialFrame> Add<&Vector<T>> for Vector<T> {
    type Output = Vector<T>;
    #[inline(always)]
    fn add(mut self, rhs: &Vector<T>) -> Self::Output {
        (0..3).for_each(|i| self.vec[i] += rhs.vec[i]);
        self
    }
}

impl<T: InertialFrame> Add<Vector<T>> for Vector<T> {
    type Output = Vector<T>;
    #[inline(always)]
    fn add(mut self, rhs: Vector<T>) -> Self::Output {
        (0..3).for_each(|i| self.vec[i] += rhs.vec[i]);
        self
    }
}

impl<T: InertialFrame> Div<f64> for Vector<T> {
    type Output = Vector<T>;
    #[inline(always)]
    fn div(mut self, rhs: f64) -> Self::Output {
        (0..3).for_each(|i| self.vec[i] /= rhs);
        self
    }
}

impl<T: InertialFrame> Div<f64> for &Vector<T> {
    type Output = Vector<T>;
    #[inline(always)]
    fn div(self, rhs: f64) -> Self::Output {
        let mut vec = self.vec;
        (0..3).for_each(|i| vec[i] /= rhs);
        vec.into()
    }
}

impl<T: InertialFrame> Mul<f64> for Vector<T> {
    type Output = Vector<T>;
    #[inline(always)]
    fn mul(mut self, rhs: f64) -> Self::Output {
        (0..3).for_each(|i| unsafe { *self.vec.get_unchecked_mut(i) *= rhs });
        self
    }
}

impl<T: InertialFrame> Mul<f64> for &Vector<T> {
    type Output = Vector<T>;
    #[inline(always)]
    fn mul(self, rhs: f64) -> Self::Output {
        let mut vec = self.vec;
        (0..3).for_each(|i| unsafe { *vec.get_unchecked_mut(i) *= rhs });
        vec.into()
    }
}

impl<T: InertialFrame> Mul<Vector<T>> for f64 {
    type Output = Vector<T>;
    #[inline(always)]
    fn mul(self, mut rhs: Vector<T>) -> Self::Output {
        (0..3).for_each(|i| unsafe { *rhs.vec.get_unchecked_mut(i) *= self });
        rhs
    }
}

impl<T: InertialFrame> Mul<&Vector<T>> for f64 {
    type Output = Vector<T>;
    #[inline(always)]
    fn mul(self, rhs: &Vector<T>) -> Self::Output {
        let mut vec = rhs.vec;
        (0..3).for_each(|i| unsafe { *vec.get_unchecked_mut(i) *= self });
        vec.into()
    }
}

impl<T: InertialFrame> MulAssign<f64> for &mut Vector<T> {
    #[inline(always)]
    fn mul_assign(&mut self, rhs: f64) {
        self.vec.iter_mut().for_each(|v| *v *= rhs);
    }
}

impl<T: InertialFrame> DivAssign<f64> for &mut Vector<T> {
    #[inline(always)]
    fn div_assign(&mut self, rhs: f64) {
        self.vec.iter_mut().for_each(|v| *v /= rhs);
    }
}

impl<T: InertialFrame> MulAssign<f64> for Vector<T> {
    #[inline(always)]
    fn mul_assign(&mut self, rhs: f64) {
        self.vec.iter_mut().for_each(|v| *v *= rhs);
    }
}

impl<T: InertialFrame> DivAssign<f64> for Vector<T> {
    #[inline(always)]
    fn div_assign(&mut self, rhs: f64) {
        self.vec.iter_mut().for_each(|v| *v /= rhs);
    }
}

impl<T: InertialFrame> Neg for &Vector<T> {
    type Output = Vector<T>;
    #[inline(always)]
    fn neg(self) -> Self::Output {
        let mut vec = self.vec;
        (0..3).for_each(|i| vec[i] = -vec[i]);
        vec.into()
    }
}

impl<T: InertialFrame> Neg for Vector<T> {
    type Output = Vector<T>;
    #[inline(always)]
    fn neg(self) -> Self::Output {
        let mut vec = self.vec;
        (0..3).for_each(|i| vec[i] = -vec[i]);
        vec.into()
    }
}
