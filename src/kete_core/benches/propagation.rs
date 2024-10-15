extern crate criterion;
use std::time::Duration;

use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion};
use kete_core::prelude::*;
use kete_core::*;
use lazy_static::lazy_static;
use pprof::criterion::{Output, PProfProfiler};
use rayon::iter::{IntoParallelIterator, ParallelIterator};

lazy_static! {
    static ref CIRCULAR: State = {
        State::new(
            Desig::Name("Circular".into()),
            2451545.0,
            [0.0, 1., 0.0].into(),
            [-constants::GMS_SQRT, 0.0, 0.0].into(),
            Frame::Ecliptic,
            0,
        )
    };
    static ref ELLIPTICAL: State = {
        State::new(
            Desig::Name("Elliptical".into()),
            2451545.0,
            [0.0, 1.5, 0.0].into(),
            [-constants::GMS_SQRT, 0.0, 0.0].into(),
            Frame::Ecliptic,
            0,
        )
    };
    static ref PARABOLIC: State = {
        State::new(
            Desig::Name("Parabolic".into()),
            2451545.0,
            [0.0, 2., 0.0].into(),
            [-constants::GMS_SQRT, 0.0, 0.0].into(),
            Frame::Ecliptic,
            0,
        )
    };
    static ref HYPERBOLIC: State = {
        State::new(
            Desig::Name("Hyperbolic".into()),
            2451545.0,
            [0.0, 3., 0.0].into(),
            [-constants::GMS_SQRT, 0.0, 0.0].into(),
            Frame::Ecliptic,
            0,
        )
    };
}

fn prop_n_body_radau(state: State, dt: f64) {
    let jd = state.jd + dt;
    propagation::propagate_n_body_spk(state, jd, false, None).unwrap();
}

fn prop_n_body_vec_radau(mut state: State, dt: f64) {
    let spk = get_spk_singleton().read().unwrap();
    spk.try_change_center(&mut state, 10).unwrap();
    state.try_change_frame_mut(Frame::Ecliptic).unwrap();
    let states = vec![state.clone(); 100];
    let non_gravs = vec![None; 100];
    let jd = state.jd + dt;
    propagation::propagate_n_body_vec(states, jd, None, non_gravs).unwrap();
}

fn prop_n_body_radau_par(state: State, dt: f64) {
    let states: Vec<State> = (0..100).map(|_| state.clone()).collect();
    let _tmp: Vec<State> = states
        .into_par_iter()
        .map(|s| {
            let jd = s.jd + dt;
            propagation::propagate_n_body_spk(s, jd, false, None).unwrap()
        })
        .collect();
}

fn prop_2_body_radau(state: State, dt: f64) {
    propagation::propagation_central(&state, state.jd + dt).unwrap();
}

fn prop_2_body_kepler(state: State, dt: f64) {
    propagation::propagate_two_body(&state, state.jd + dt).unwrap();
}

pub fn two_body_numeric(c: &mut Criterion) {
    let mut twobody_num_group = c.benchmark_group("2-Body-Numeric");

    for state in [
        CIRCULAR.clone(),
        ELLIPTICAL.clone(),
        PARABOLIC.clone(),
        HYPERBOLIC.clone(),
    ] {
        let name = match &state.desig {
            Desig::Name(n) => n,
            _ => panic!(),
        };
        twobody_num_group.bench_with_input(BenchmarkId::new("Single", name), &state, |b, s| {
            b.iter(|| prop_2_body_radau(black_box(s.clone()), black_box(1000.0)))
        });
    }
}

pub fn n_body_prop(c: &mut Criterion) {
    let mut nbody_group = c.benchmark_group("N-Body");

    for state in [
        CIRCULAR.clone(),
        ELLIPTICAL.clone(),
        PARABOLIC.clone(),
        HYPERBOLIC.clone(),
    ] {
        let name = match &state.desig {
            Desig::Name(n) => n,
            _ => panic!(),
        };
        nbody_group.bench_with_input(BenchmarkId::new("Single", name), &state, |b, s| {
            b.iter(|| prop_n_body_radau(black_box(s.clone()), black_box(1000.0)))
        });

        nbody_group.bench_with_input(BenchmarkId::new("Parallel", name), &state, |b, s| {
            b.iter(|| prop_n_body_radau_par(black_box(s.clone()), black_box(1000.0)))
        });
    }
}
pub fn n_body_prop_vec(c: &mut Criterion) {
    let mut nbody_group = c.benchmark_group("N-Body-Vec");

    for state in [
        CIRCULAR.clone(),
        ELLIPTICAL.clone(),
        PARABOLIC.clone(),
        HYPERBOLIC.clone(),
    ] {
        let name = match &state.desig {
            Desig::Name(n) => n,
            _ => panic!(),
        };
        nbody_group.bench_with_input(BenchmarkId::new("Single", name), &state, |b, s| {
            b.iter(|| prop_n_body_vec_radau(black_box(s.clone()), black_box(1000.0)))
        });
    }
}

pub fn two_body_analytic(c: &mut Criterion) {
    let mut twobody_group = c.benchmark_group("2-Body-Analytic");

    for state in [
        CIRCULAR.clone(),
        ELLIPTICAL.clone(),
        PARABOLIC.clone(),
        HYPERBOLIC.clone(),
    ] {
        let name = match &state.desig {
            Desig::Name(n) => n,
            _ => panic!(),
        };
        twobody_group.bench_with_input(BenchmarkId::new("Single", name), &state, |b, s| {
            b.iter(|| prop_2_body_kepler(s.clone(), black_box(1000.0)))
        });
    }
}

criterion_group!(name=benches;
                 config = Criterion::default().sample_size(30).measurement_time(Duration::from_secs(15)).with_profiler(PProfProfiler::new(100, Output::Flamegraph(None)));
                 targets=n_body_prop_vec, two_body_analytic, n_body_prop, two_body_numeric);
criterion_main!(benches);
