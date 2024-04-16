use crate::errors::NEOSpyError;
use nalgebra::SVector;

/// Function will be of the form y' = F(t, y, metadata, exact_eval).
/// Metadata is passed for every evaluation. The `exact_eval` bool indicates to the
/// function that the input parameters are known to be solutions for the IVP.
pub type FirstOrderFunc<'a, MType, const D: usize> =
    &'a dyn Fn(f64, &SVector<f64, D>, &mut MType, bool) -> Result<SVector<f64, D>, NEOSpyError>;

/// Integrator will return a result of this type.
pub type FirstOrderResult<MType, const D: usize> = Result<(SVector<f64, D>, MType), NEOSpyError>;

/// Runge-Kutta-Fehlberg Integrator - RK4(5)
/// <https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta%E2%80%93Fehlberg_method>
#[allow(missing_debug_implementations)]
pub struct RK45Integrator<'a, MType, const D: usize> {
    func: FirstOrderFunc<'a, MType, D>,
    metadata: MType,

    final_time: f64,

    cur_time: f64,
    cur_state: SVector<f64, D>,
    cur_der: SVector<f64, D>,

    tol: f64,
}

impl<'a, MType, const D: usize> RK45Integrator<'a, MType, D> {
    /// Attempt to integrate by the provided step size.
    /// This may take a smaller step that specified if it failed to converge with the
    /// specified tolerance.
    fn step(&mut self, step_size: f64) -> Result<f64, NEOSpyError> {
        if step_size < 1e-30 {
            return Err(NEOSpyError::Convergence(
                "Runge-Kutta step size too small".into(),
            ));
        }

        let k2 = step_size
            * (self.func)(
                self.cur_time + step_size / 4.0,
                &(self.cur_state + self.cur_der * (step_size / 4.0)),
                &mut self.metadata,
                false,
            )?;

        let k3 = step_size
            * (self.func)(
                self.cur_time + step_size * (3.0 / 8.0),
                &(self.cur_state + self.cur_der * (3.0 * step_size / 32.0) + k2 * (9.0 / 32.0)),
                &mut self.metadata,
                false,
            )?;

        let k4 = step_size
            * (self.func)(
                self.cur_time + step_size * (12.0 / 13.0),
                &(self.cur_state + self.cur_der * (step_size * 1932.0 / 2197.0)
                    - k2 * (7200.0 / 2197.0)
                    + k3 * (7296.0 / 2197.0)),
                &mut self.metadata,
                false,
            )?;

        let k5 = step_size
            * (self.func)(
                self.cur_time + step_size,
                &(self.cur_state + self.cur_der * (step_size * 439.0 / 216.0) - k2 * 8.0
                    + k3 * (3680.0 / 513.0)
                    - k4 * (845.0 / 4104.0)),
                &mut self.metadata,
                false,
            )?;

        let k6 = step_size
            * (self.func)(
                self.cur_time + step_size * 0.5,
                &(self.cur_state - self.cur_der * (step_size * 8.0 / 27.0) + k2 * 2.0
                    - k3 * (3544.0 / 2565.0)
                    + k4 * (1859.0 / 4104.0)
                    - k5 * (11.0 / 40.0)),
                &mut self.metadata,
                false,
            )?;

        let trunc_err =
            (self.cur_der * (step_size / -360.) + k3 * (128.0 / 4275.0) + k4 * (2197.0 / 75240.0)
                - k5 / 50.0
                - k6 * (2.0 / 55.0))
                .norm();

        let new_step: f64;
        if trunc_err == 0.0 {
            new_step = 2.0 * step_size;
        } else {
            new_step = 0.9 * step_size * (self.tol / trunc_err).powf(1.0 / 5.0);
            if trunc_err > self.tol {
                return self.step(new_step);
            }
        }

        self.cur_state += self.cur_der * (step_size * 16.0 / 135.0)
            + k3 * (6656.0 / 12825.0)
            + k4 * (28561.0 / 56430.0)
            - k5 * (9.0 / 50.0)
            + k6 * (2.0 / 55.0);
        self.cur_time += step_size;
        self.cur_der = (self.func)(self.cur_time, &self.cur_state, &mut self.metadata, true)?;

        Ok(new_step)
    }

    /// Adaptive step size Runge-Kutta - RK4(5)
    /// This performs integration on first order Initial Value Problems (IVP).
    ///
    /// This is a relatively fast, but not very precise numerical integration method.
    pub fn integrate(
        func: FirstOrderFunc<'a, MType, D>,
        state_init: SVector<f64, D>,
        time_init: f64,
        final_time: f64,
        mut metadata: MType,
        tol: f64,
    ) -> FirstOrderResult<MType, D> {
        if (time_init - final_time).abs() < 1e-30 {
            return Ok((state_init, metadata));
        }
        let cur_der = func(time_init, &state_init, &mut metadata, true)?;
        let mut integrator = Self {
            func,
            metadata,
            final_time,
            cur_time: time_init,
            cur_state: state_init,
            cur_der,
            tol,
        };
        let mut next_step = final_time - integrator.cur_time;
        next_step = integrator.step(next_step)?;
        while (integrator.cur_time - integrator.final_time).abs() > tol {
            if (integrator.cur_time - integrator.final_time).abs() < next_step.abs() {
                next_step = integrator.final_time - integrator.cur_time;
            }
            next_step = integrator.step(next_step)?;
        }

        Ok((integrator.cur_state, integrator.metadata))
    }
}

#[cfg(test)]
mod tests {
    use nalgebra::Vector1;

    use crate::errors::NEOSpyError;
    use crate::propagation::RK45Integrator;

    #[test]
    fn basic_exp() {
        fn f(
            _t: f64,
            state: &Vector1<f64>,
            _meta: &mut (),
            _eval: bool,
        ) -> Result<Vector1<f64>, NEOSpyError> {
            Ok(-state)
        }

        for step in 1..4 {
            let res = RK45Integrator::integrate(
                &f,
                Vector1::<f64>::repeat(1f64),
                0.0,
                step as f64 * 10.0,
                (),
                1e-14,
            );
            assert!(res.is_ok());
            assert!((res.unwrap().0[0] - (-step as f64 * 10.0).exp()).abs() < 1e-14)
        }
    }

    #[test]
    fn basic_line() {
        fn f(
            _t: f64,
            _state: &Vector1<f64>,
            _meta: &mut (),
            _eval: bool,
        ) -> Result<Vector1<f64>, NEOSpyError> {
            Ok(Vector1::<f64>::repeat(2.0))
        }

        let res = RK45Integrator::integrate(&f, Vector1::<f64>::repeat(1f64), 0.0, 10.0, (), 0.001);

        assert!(res.is_ok());
        assert!((res.unwrap().0[0] - 21.0).abs() < 1e-12);
    }
}
