/// Gauss-Radau Spacing Numerical Integrator
/// This solves a second-order initial value problem.
use crate::errors::NEOSpyError;
use itertools::izip;
use lazy_static::lazy_static;
use nalgebra::allocator::Allocator;
use nalgebra::*;
use nalgebra::{DefaultAllocator, Dim, OMatrix, OVector, RowSVector, SMatrix, SVector, U1, U7};

/// Function will be of the form y'' = F(t, y, y', metadata)
/// This is the second-order general solver (class II in the Everhart paper).
pub type RadauFunc<'a, MType, D> = &'a dyn Fn(
    f64,
    &OVector<f64, D>,
    &OVector<f64, D>,
    &mut MType,
    bool,
) -> Result<OVector<f64, D>, NEOSpyError>;

/// Integrator will return a result of this type.
pub type RadauResult<MType, D> = Result<(OVector<f64, D>, OVector<f64, D>, MType), NEOSpyError>;

const GAUSS_RADAU_SPACINGS: [f64; 8] = [
    0.0,
    0.05626256053692215,
    0.18024069173689236,
    0.3526247171131696,
    0.5471536263305554,
    0.7342101772154105,
    0.8853209468390958,
    0.9775206135612875,
];

lazy_static! {
    // this initializes w, u, c, d, r
    static ref W_VEC: RowSVector<f64, 7> = {
        let mut w = RowSVector::<f64, 7>::zeros();
        for (idx, e) in w.iter_mut().enumerate() {
            *e = (((idx + 2) * (idx + 3)) as f64).recip();
        }
        w
    };

    static ref U_VEC: RowSVector<f64, 7> = {
        let mut u = RowSVector::<f64, 7>::zeros();
        for (idx, e) in u.iter_mut().enumerate() {
            *e = ((idx + 2) as f64).recip();
        }
        u
    };


    static ref C_MAT: SMatrix<f64, 7, 7> = {
        let mut c = SMatrix::<f64, 7, 7>::identity();
        for idx in 0..7 {
            if idx > 0 {
                c[(idx, 0)] = -GAUSS_RADAU_SPACINGS[idx] * c[(idx - 1, 0)];
            }
            for idy in 1..idx {
                c[(idx, idy)] = c[(idx - 1, idy - 1)]
                    - GAUSS_RADAU_SPACINGS[idx] * c[(idx - 1, idy)];
            }
        }
        c
    };

}

const MIN_RATIO: f64 = 0.25;
const EPSILON: f64 = 1e-6;
const MIN_STEP: f64 = 0.01;

/// Gauss-Radau Spacing Numerical Integrator
/// This solves a second-order initial value problem.
///
/// References:
/// E. Everhart (1985), 'An efficient integrator that uses Gauss-Radau spacings',
/// A. Carusi and G. B. Valsecchi (eds.),
/// Dynamics of Comets: Their Origin and Evolution (proceedings),
/// Astrophysics and Space Science Library, vol. 115, D. Reidel Publishing Company
///
/// E. Everhart (1974), 'Implicit single-sequence methods for integrating orbits',
/// Celestial Mechanics, vol. 10, no. 1, pp. 35-55
///
/// This uses the 15th-order integrator as seen in the original RADAU code, however
/// many changes and improvements have been made. Some variable names have been chosen
/// to match the original Fortran implementation. After some experimentation it was
/// found that the correction and prediction steps turned out to in general help such to
/// a small degree that they were not worth the added complexity.
#[allow(missing_debug_implementations)]
pub struct RadauIntegrator<'a, MType, D: Dim>
where
    DefaultAllocator: Allocator<f64, D, U1> + Allocator<f64, D, U7>,
{
    func: RadauFunc<'a, MType, D>,
    metadata: MType,

    final_time: f64,

    cur_time: f64,
    cur_state: OVector<f64, D>,
    cur_state_der: OVector<f64, D>,
    cur_state_der_der: OVector<f64, D>,

    cur_b: OMatrix<f64, D, U7>,
    g_scratch: OMatrix<f64, D, U7>,

    state_scratch: OVector<f64, D>,
    state_der_scratch: OVector<f64, D>,
    b_scratch: OVector<f64, D>,
    eval_scratch: OVector<f64, D>,
}

impl<'a, MType, D: Dim> RadauIntegrator<'a, MType, D>
where
    DefaultAllocator: Allocator<f64, D, U1> + Allocator<f64, D, U7>,
{
    fn new(
        func: RadauFunc<'a, MType, D>,
        state_init: OVector<f64, D>,
        state_der_init: OVector<f64, D>,
        time_init: f64,
        final_time: f64,
        metadata: MType,
    ) -> Result<Self, NEOSpyError> {
        let (dim, _) = state_init.shape_generic();
        if state_init.len() != state_der_init.len() {
            return Err(NEOSpyError::ValueError(
                "Input vectors must be the same length".into(),
            ));
        }
        let mut res = Self {
            func,
            metadata,
            final_time,
            cur_time: time_init,
            cur_state: state_init,
            cur_state_der: state_der_init,
            cur_state_der_der: Matrix::zeros_generic(dim, U1),
            cur_b: Matrix::zeros_generic(dim, U7),
            g_scratch: Matrix::zeros_generic(dim, U7),
            b_scratch: Matrix::zeros_generic(dim, U1),
            state_scratch: Matrix::zeros_generic(dim, U1),
            state_der_scratch: Matrix::zeros_generic(dim, U1),
            eval_scratch: Matrix::zeros_generic(dim, U1),
        };

        res.cur_state_der_der = (res.func)(
            time_init,
            &res.cur_state,
            &res.cur_state_der,
            &mut res.metadata,
            true,
        )?;
        Ok(res)
    }

    /// Integrate the functions from the initial time to the final time.
    pub fn integrate(
        func: RadauFunc<'a, MType, D>,
        state_init: OVector<f64, D>,
        state_der_init: OVector<f64, D>,
        time_init: f64,
        final_time: f64,
        metadata: MType,
    ) -> RadauResult<MType, D> {
        let mut integrator = Self::new(
            func,
            state_init,
            state_der_init,
            time_init,
            final_time,
            metadata,
        )?;
        if (final_time - time_init).abs() < 1e-10 {
            return Ok((
                integrator.cur_state,
                integrator.cur_state_der,
                integrator.metadata,
            ));
        }
        let mut next_step_size: f64 = 1.0_f64.copysign(integrator.final_time - integrator.cur_time);

        let mut step_failures = 0;
        loop {
            if (integrator.cur_time - integrator.final_time).abs() <= next_step_size.abs() {
                next_step_size = integrator.final_time - integrator.cur_time;
            }
            match integrator.step(next_step_size) {
                Ok(s) => {
                    next_step_size = s;
                    if integrator.cur_time == integrator.final_time {
                        return Ok((
                            integrator.cur_state,
                            integrator.cur_state_der,
                            integrator.metadata,
                        ));
                    }
                    step_failures = 0;
                }
                Err(error) => match error {
                    NEOSpyError::Impact(_, _) => return Err(error),
                    NEOSpyError::DAFLimits(_) => return Err(error),
                    _ => {
                        step_failures += 1;
                        next_step_size *= 0.7;
                        if step_failures > 10 {
                            return Err(NEOSpyError::Convergence(
                                "Radau failed to converge.".into(),
                            ));
                        }
                    }
                },
            }
            if next_step_size.abs() < MIN_STEP {
                next_step_size = MIN_STEP.copysign(next_step_size);
            }
        }
    }

    /// If this function fails, then the step guess is almost certainly too large.
    ///
    /// This function will update the current b matrices to be correct for the
    /// step guess provided.
    ///
    fn step(&mut self, step_size: f64) -> Result<f64, NEOSpyError> {

        self.g_scratch.fill(0.0);
        self.state_scratch.fill(0.0);
        self.state_der_scratch.fill(0.0);
        self.eval_scratch.set_column(0, &self.cur_state_der_der);

        for _ in 0..10 {
            self.b_scratch.set_column(0, &self.cur_b.column(6));
            // Calculate b and g
            for (idj, gauss_radau_frac) in GAUSS_RADAU_SPACINGS.iter().enumerate().skip(1) {
                // the sample point at the Guass-Radau spacings.
                // Update each parameter using the current B as a guess to estimate the
                // state of the integrator at the current time + the Gauss-Radau spacing.

                let w_pow = SVector::from_iterator(
                    W_VEC
                        .iter()
                        .enumerate()
                        .map(|(idx, w)| gauss_radau_frac.powi(1 + idx as i32) * w),
                );
                let u_pow = SVector::from_iterator(
                    U_VEC
                        .iter()
                        .enumerate()
                        .map(|(idx, u)| gauss_radau_frac.powi(1 + idx as i32) * u),
                );

                izip!(
                    self.state_scratch.iter_mut(),
                    self.cur_state.iter(),
                    self.cur_state_der.iter(),
                    self.cur_state_der_der.iter(),
                    self.cur_b.row_iter(),
                )
                .for_each(|(out, state, der, derder, b)| {
                    *out = state
                        + gauss_radau_frac * der * step_size
                        + gauss_radau_frac.powi(2)
                            * step_size.powi(2)
                            * (derder / 2.0 + unsafe { (b * w_pow).get_unchecked(0) })
                });

                izip!(
                    self.state_der_scratch.iter_mut(),
                    self.cur_state_der.iter(),
                    self.cur_state_der_der.iter(),
                    self.cur_b.row_iter(),
                )
                .for_each(|(out, der, derder, b)| {
                    *out = der
                        + gauss_radau_frac
                            * step_size
                            * (derder + unsafe { (b * u_pow).get_unchecked(0) })
                });

                // Evaluate the function at this new intermediate state.
                self.eval_scratch.set_column(0, &(self.func)(
                    self.cur_time + gauss_radau_frac * step_size,
                    &self.state_scratch,
                    &self.state_der_scratch,
                    &mut self.metadata,
                    false,
                )?);

                let diff = &self.eval_scratch - &self.cur_state_der_der;

                // Use the result of that evaluation to update the current G and B
                // matrices for the next gauss spacing.

                // This is equivalent to equation (4) in everhart's paper.
                // The lookup tables and switch statements he uses were performing
                // ~100x slower than this implementation.
                self.g_scratch.set_column(idj - 1, &{
                    let mut gk = diff / *gauss_radau_frac;

                    for (idz, gr_step) in GAUSS_RADAU_SPACINGS.iter().enumerate().take(idj).skip(1)
                    {
                        gk = (gk - self.g_scratch.column(idz - 1)) / (gauss_radau_frac - gr_step);
                    }
                    gk
                });
            }

            // Compute the largest b value and compare it to the last b value.
            // convergence is decided if B stops changing.
            self.g_scratch.mul_to(&C_MAT, &mut self.cur_b);
            let b_diff = (self.cur_b.column(6) - &self.b_scratch).abs();
            let func_eval_max = self.eval_scratch.abs().add_scalar(1e-6);

            // This is using the convergence criterion as defined in
            // https://arxiv.org/pdf/1409.4779.pdf  equation (8)
            if b_diff.component_div(&func_eval_max).max() < 1e-14 {
                for idx in 0..self.cur_state.len() {
                    unsafe {
                        self.cur_state[idx] += self.cur_state_der.get_unchecked(idx) * step_size
                            + step_size.powi(2)
                                * (self.cur_state_der_der.get_unchecked(idx) * 0.5
                                    + self.cur_b.row(idx).dot(&W_VEC));
                        self.cur_state_der[idx] += step_size
                            * (self.cur_state_der_der.get_unchecked(idx)
                                + self.cur_b.row(idx).dot(&U_VEC));
                    }
                }
                self.cur_time += step_size;
                self.cur_state_der_der = (self.func)(
                    self.cur_time,
                    &self.cur_state,
                    &self.cur_state_der,
                    &mut self.metadata,
                    true,
                )?;
                return Ok(step_size
                    * (EPSILON * func_eval_max.max() / self.cur_b.column(6).abs().max())
                        .powf(1.0 / 7.0)
                        .clamp(MIN_RATIO, MIN_RATIO.recip()));
            }
        }
        Err(NEOSpyError::Convergence(
            "Radau step failed to converge".into(),
        ))
    }
}

#[cfg(test)]
mod tests {
    use nalgebra::Vector3;

    use super::*;
    use crate::propagation::{central_accel, CentralAccelMeta};

    #[test]
    fn basic_two_body() {
        let (pos, vel, _meta) = RadauIntegrator::integrate(
            &central_accel,
            Vector3::new(0.46937657, -0.8829981, 0.),
            Vector3::new(0.01518942, 0.00807426, 0.),
            0.0,
            1000.,
            CentralAccelMeta::default(),
        )
        .unwrap();
        assert!((pos[0] + 0.916350120888658).abs() < 1e-8);
        assert!((pos[1] + 0.4003771936559588).abs() < 1e-8);
        assert_eq!(pos[2], 0.0);

        assert!((vel[0] - 0.006887328686018099).abs() < 1e-8);
        assert!((vel[1] + 0.01576315407302832).abs() < 1e-8);
        assert_eq!(vel[2], 0.0);
    }
}
