use criterion::{Criterion, criterion_group, criterion_main};
use std::hint::black_box;

pub fn cubic_roots(c: &mut Criterion) {
    let poly = poly_cool::Cubic {
        c3: 1.0,
        c2: -6.0,
        c1: 11.0,
        c0: -6.0,
    };

    c.bench_function("roots_between_multiple_searches", |b| {
        b.iter(|| black_box(poly).roots_between_multiple_searches(-1.0, 4.0, 1e-12))
    });
    c.bench_function("roots_between", |b| {
        b.iter(|| black_box(poly).roots_between(-1.0, 4.0, 1e-12))
    });
    c.bench_function("roots_between with one root", |b| {
        b.iter(|| black_box(poly).roots_between(-1.0, 1.5, 1e-12))
    });
    c.bench_function("roots_blinn", |b| b.iter(|| black_box(poly).roots_blinn()));
    c.bench_function("roots_blinn_and_deflate", |b| {
        b.iter(|| black_box(poly).roots_blinn_and_deflate())
    });
    c.bench_function("kurbo", |b| {
        b.iter(|| {
            kurbo::common::solve_cubic(
                black_box(poly.c0),
                black_box(poly.c1),
                black_box(poly.c2),
                black_box(poly.c3),
            )
        })
    });
}

criterion_group!(benches, cubic_roots);
criterion_main!(benches);
