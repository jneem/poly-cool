use criterion::{Criterion, criterion_group, criterion_main};
use std::hint::black_box;

pub fn quadratic_roots(c: &mut Criterion) {
    let poly = poly_cool::Quadratic {
        c2: -6.0,
        c1: 11.0,
        c0: -6.0,
    };

    c.bench_function("full", |b| b.iter(|| black_box(poly).roots()));
    c.bench_function("positive discriminant", |b| {
        b.iter(|| black_box(poly).positive_discriminant_roots())
    });
    c.bench_function("positive discriminant no overflow check", |b| {
        b.iter(|| black_box(poly).positive_discriminant_roots_no_overflow_check())
    });
    c.bench_function("positive discriminant no overflow check", |b| {
        b.iter(|| black_box(poly).positive_discriminant_roots_no_overflow_check_half_b())
    });
    c.bench_function("kurbo", |b| {
        b.iter(|| {
            kurbo::common::solve_quadratic(
                black_box(poly.c0),
                black_box(poly.c1),
                black_box(poly.c2),
            )
        })
    });
}

criterion_group!(benches, quadratic_roots);
criterion_main!(benches);
