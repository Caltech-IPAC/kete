//! Basic Shape Models
use crate::{
    constants::GOLDEN_RATIO,
    frames::{Equatorial, InertialFrame, UnitVector, Vector},
};
use lazy_static::lazy_static;
use std::f64::consts::TAU;

lazy_static! {
    /// Pre-compute a default shape.
    pub static ref DEFAULT_SHAPE: ConvexShape<Equatorial> = ConvexShape::new_fibonacci_lattice(2048);
}

/// Facet of a shape.
#[derive(Debug, Clone)]
pub struct Facet<T: InertialFrame> {
    /// Normal unit vector defining the facets face
    pub normal: UnitVector<T>,

    /// Surface area of the facet
    pub area: f64,
}

/// Convex shape made up of individual facets.
#[derive(Debug, Clone)]
pub struct ConvexShape<T: InertialFrame> {
    /// The facets which make up this shape.
    pub facets: Box<[Facet<T>]>,
}

impl<T: InertialFrame> ConvexShape<T> {
    /// Construct a new ConvexShape using fibonacci lattice spacing.
    ///
    /// Evenly place points on a sphere using the Fibonacci Lattice algorithm.
    ///
    /// This uses a slightly modified method where an epsilon term is added to shift the
    /// points slightly, causing the average spacing between the points to be minimized.
    ///
    /// See:
    /// <http://extremelearning.com.au/how-to-evenly-distribute-points-on-a-sphere-more-effectively-than-the-canonical-fibonacci-lattice/>
    ///
    ///
    /// Total surface area is set to 1.
    pub fn new_fibonacci_lattice(n_facets: usize) -> Self {
        let mut facets: Vec<Facet<T>> = Vec::with_capacity(n_facets);

        const EPSILON: f64 = 0.36;

        let n_normals = n_facets as f64;
        let area = n_normals.recip();

        for idx in 0..n_facets {
            let theta: f64 = TAU * (idx as f64) / GOLDEN_RATIO;
            let phi: f64 =
                (1.0 - 2.0 * ((idx as f64) + EPSILON) / (n_normals - 1.0 + 2.0 * EPSILON)).acos();
            let normal = UnitVector::new_unchecked(Vector::new([
                theta.cos() * phi.sin(),
                theta.sin() * phi.sin(),
                phi.cos(),
            ]));

            facets.push(Facet { normal, area });
        }

        Self {
            facets: facets.into(),
        }
    }

    /// Rescale the total areas to sum to 1.
    pub fn normalize_areas(&mut self) {
        let total_area_inv = self.facets.iter().map(|f| f.area).sum::<f64>().recip();
        self.facets
            .iter_mut()
            .for_each(|x| x.area *= total_area_inv);
    }
}

#[cfg(test)]
mod tests {

    use crate::frames::Equatorial;

    use super::*;

    #[test]
    fn test_convex_shape() {
        let n1024: ConvexShape<Equatorial> = ConvexShape::new_fibonacci_lattice(1024);

        assert!(n1024.facets.len() == 1024);
        assert!(n1024.facets.iter().all(|x| x.area == (1024f64).recip()))
    }
}
