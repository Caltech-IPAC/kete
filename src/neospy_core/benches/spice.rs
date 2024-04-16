extern crate criterion;
use criterion::{black_box, criterion_group, criterion_main, Criterion};
use neospy_core::spice::get_spk_singleton;
use pprof::criterion::{Output, PProfProfiler};

fn spice_get_raw_state(jd: f64) {
    let spice = get_spk_singleton().try_read().unwrap();
    for _ in 0..1000 {
        let _ = spice.try_get_raw_state(5, jd).unwrap();
    }
}

fn spice_get_state(jd: f64) {
    let spice = get_spk_singleton().try_read().unwrap();
    for _ in 0..1000 {
        let _ = spice
            .try_get_state(5, jd, 10, neospy_core::frames::Frame::Ecliptic)
            .unwrap();
    }
}

pub fn spice_benchmark(c: &mut Criterion) {
    c.bench_function("spice_get_raw_state", |b| {
        b.iter(|| spice_get_raw_state(black_box(2451545.0)))
    });
    c.bench_function("spice_get_state", |b| {
        b.iter(|| spice_get_state(black_box(2451545.0)))
    });
}

criterion_group!(name=spice;
                config = Criterion::default().with_profiler(PProfProfiler::new(100, Output::Flamegraph(None)));
                targets=spice_benchmark);
criterion_main!(spice);
