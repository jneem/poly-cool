use criterion::{criterion_group, criterion_main, Criterion};
use std::hint::black_box;

pub fn quadratic_roots(c: &mut Criterion) {
    let poly = poly_cool::Quadratic::new([-6.0, 11.0, -6.0]);

    c.bench_function("full", |b| b.iter(|| black_box(poly).roots()));
    c.bench_function("positive discriminant", |b| {
        b.iter(|| black_box(poly).positive_discriminant_roots())
    });
    c.bench_function("kurbo", |b| {
        b.iter(|| {
            let &[c0, c1, c2] = poly.coeffs();
            kurbo::common::solve_quadratic(black_box(c0), black_box(c1), black_box(c2))
        })
    });
}

criterion_group!(benches, quadratic_roots);
criterion_main!(benches);
