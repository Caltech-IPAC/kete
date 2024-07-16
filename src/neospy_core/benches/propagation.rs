extern crate criterion;
use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion};
use frames::Ecliptic;
use lazy_static::lazy_static;
use neospy_core::prelude::*;
use neospy_core::*;
use pprof::criterion::{Output, PProfProfiler};
use rayon::iter::{IntoParallelIterator, ParallelIterator};

lazy_static! {
    static ref CIRCULAR: State<Ecliptic> = {
        State::new(
            Some(Desig::Name("Circular".into())),
            2451545.0,
            [0.0, 1., 0.0].into(),
            [-constants::GMS_SQRT, 0.0, 0.0].into(),
            0,
        )
    };
    static ref ELLIPTICAL: State<Ecliptic> = {
        State::<Ecliptic>::new(
            Some(Desig::Name("Elliptical".into())),
            2451545.0,
            [0.0, 1.5, 0.0].into(),
            [-constants::GMS_SQRT, 0.0, 0.0].into(),
            0,
        )
    };
    static ref PARABOLIC: State<Ecliptic> = {
        State::new(
            Some(Desig::Name("Parabolic".into())),
            2451545.0,
            [0.0, 2., 0.0].into(),
            [-constants::GMS_SQRT, 0.0, 0.0].into(),
            0,
        )
    };
    static ref HYPERBOLIC: State<Ecliptic> = {
        State::new(
            Some(Desig::Name("Hyperbolic".into())),
            2451545.0,
            [0.0, 3., 0.0].into(),
            [-constants::GMS_SQRT, 0.0, 0.0].into(),
            0,
        )
    };
}

fn prop_n_body_radau(state: State<Ecliptic>, dt: f64) {
    let jd = state.jd + dt;
    propagation::propagate_n_body_spk(state.into_frame(), jd, false, None).unwrap();
}

fn prop_n_body_radau_par(state: State<Ecliptic>, dt: f64) {
    let states: Vec<State<_>> = (0..100).map(|_| state.clone()).collect();
    let _tmp: Vec<State<_>> = states
        .into_par_iter()
        .map(|s| {
            let jd = s.jd + dt;
            propagation::propagate_n_body_spk(s.into_frame(), jd, false, None).unwrap()
        })
        .collect();
}

fn prop_2_body_radau(state: State<Ecliptic>, dt: f64) {
    let jd = state.jd + dt;
    propagation::propagation_central(&state.into_frame(), jd).unwrap();
}

fn prop_2_body_kepler(state: State<Ecliptic>, dt: f64) {
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
            Some(Desig::Name(n)) => n,
            _ => panic!(),
        };
        twobody_num_group.bench_with_input(BenchmarkId::new("Single", name), &state, |b, s| {
            b.iter(|| prop_2_body_radau(s.clone(), black_box(1000.0)))
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
            Some(Desig::Name(n)) => n,
            _ => panic!(),
        };
        nbody_group.bench_with_input(BenchmarkId::new("Single", name), &state, |b, s| {
            b.iter(|| prop_n_body_radau(s.clone(), black_box(1000.0)))
        });

        nbody_group.bench_with_input(BenchmarkId::new("Parallel", name), &state, |b, s| {
            b.iter(|| prop_n_body_radau_par(black_box(s.clone()), black_box(1000.0)))
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
            Some(Desig::Name(n)) => n,
            _ => panic!(),
        };
        twobody_group.bench_with_input(BenchmarkId::new("Single", name), &state, |b, s| {
            b.iter(|| prop_2_body_kepler(s.clone(), black_box(1000.0)))
        });
    }
}

criterion_group!(name=benches;
                 config = Criterion::default().with_profiler(PProfProfiler::new(100, Output::Flamegraph(None)));
                 targets=two_body_analytic, n_body_prop, two_body_numeric);
criterion_main!(benches);
