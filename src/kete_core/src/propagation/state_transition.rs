use crate::constants::GMS_SQRT;
use crate::prelude::{KeteResult, State};
use crate::propagation::{central_accel, central_accel_grad, CentralAccelMeta, RK45Integrator};
use nalgebra::{Const, Matrix6, SVector, Vector3, U1, U6};

fn stm_ivp_eqn(
    jd: f64,
    state: &SVector<f64, 42>,
    meta: &mut CentralAccelMeta,
    exact_eval: bool,
) -> KeteResult<SVector<f64, 42>> {
    let mut res = SVector::<f64, 42>::zeros();

    // first 6 values of the state are pos and vel respectively.
    let pos = Vector3::new(state[0], state[1], state[2]);
    let vel = Vector3::new(state[3], state[4], state[5]);
    let accel = central_accel(jd, &pos, &vel, meta, exact_eval)?;

    // the derivative of pos is the velocity, and the derivative of vel is the acceleration
    // set those as appropriate for the output state
    res.fixed_rows_mut::<3>(0).set_column(0, &vel);
    res.fixed_rows_mut::<3>(3).set_column(0, &accel);

    // the remainder of res is the state transition matrix calculation.
    let mut stm = Matrix6::<f64>::zeros();
    stm.fixed_view_mut::<3, 3>(0, 3)
        .set_diagonal(&Vector3::repeat(1.0));

    let grad = central_accel_grad(0.0, &pos, &vel, meta);
    let mut view = stm.fixed_view_mut::<3, 3>(3, 0);
    view.set_row(0, &grad.row(0));
    view.set_row(1, &grad.row(1));
    view.set_row(2, &grad.row(2));

    let vec_reshape = state
        .fixed_rows::<36>(6)
        .into_owned()
        .reshape_generic(U6, U6);
    res.rows_mut(6, 36).set_column(
        0,
        &(stm * vec_reshape)
            .into_owned()
            .reshape_generic(Const::<36>, U1),
    );

    Ok(res)
}

/// Compute a state transition matrix assuming only 2-body mechanics.
///
/// This uses a Runge-Kutta 4/5 algorithm.
pub fn compute_state_transition(
    state: &mut State,
    jd: f64,
    central_mass: f64,
) -> ([[f64; 3]; 2], Matrix6<f64>) {
    let meta = CentralAccelMeta {
        mass_scaling: central_mass,
        ..Default::default()
    };

    let mut initial_state = SVector::<f64, 42>::zeros();

    initial_state.rows_mut(6, 36).set_column(
        0,
        &Matrix6::<f64>::identity().reshape_generic(Const::<36>, U1),
    );

    initial_state
        .fixed_rows_mut::<3>(0)
        .set_column(0, &state.pos.into());
    initial_state
        .fixed_rows_mut::<3>(3)
        .set_column(0, &(Vector3::from(state.vel) / GMS_SQRT));
    let rad = RK45Integrator::integrate(
        &stm_ivp_eqn,
        initial_state,
        state.jd * GMS_SQRT,
        jd * GMS_SQRT,
        meta,
        1e-12,
    )
    .unwrap();

    let vec_reshape = rad
        .0
        .fixed_rows::<36>(6)
        .into_owned()
        .reshape_generic(U6, U6)
        .transpose();

    let scaling_a = Matrix6::<f64>::from_diagonal(
        &[
            1.0,
            1.0,
            1.0,
            1.0 / GMS_SQRT,
            1.0 / GMS_SQRT,
            1.0 / GMS_SQRT,
        ]
        .into(),
    );
    let scaling_b =
        Matrix6::<f64>::from_diagonal(&[1.0, 1.0, 1.0, GMS_SQRT, GMS_SQRT, GMS_SQRT].into());
    (
        [
            rad.0.fixed_rows::<3>(0).into(),
            (rad.0.fixed_rows::<3>(3) * GMS_SQRT).into(),
        ],
        scaling_a * vec_reshape * scaling_b,
    )
}
