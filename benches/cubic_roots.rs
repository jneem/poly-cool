use criterion::{criterion_group, criterion_main, Criterion};
use std::hint::black_box;

pub fn cubic_roots(c: &mut Criterion) {
    let poly = poly_cool::Cubic::new([-6.0, 11.0, -6.0, 1.0]);

    let mut group = c.benchmark_group("simple roots");

    for accuracy in [1e-6, 1e-8, 1e-12] {
        group.bench_with_input(
            format!("roots_between {accuracy:?}"),
            &accuracy,
            |b, accuracy| b.iter(|| black_box(poly).roots_between(-1.0, 4.0, *accuracy)),
        );
    }
    group.bench_function("roots_between with one root", |b| {
        b.iter(|| black_box(poly).roots_between(-1.0, 1.5, 1e-12))
    });
    group.bench_function("roots_blinn", |b| b.iter(|| black_box(poly).roots_blinn()));
    group.bench_function("roots_blinn with preconditioning", |b| {
        b.iter(|| black_box(poly).precondition().roots_blinn())
    });
    group.bench_function("roots_blinn_and_deflate", |b| {
        b.iter(|| black_box(poly).roots_blinn_and_deflate())
    });
    group.bench_function("kurbo", |b| {
        let &[c0, c1, c2, c3] = poly.coeffs();
        b.iter(|| {
            kurbo::common::solve_cubic(black_box(c0), black_box(c1), black_box(c2), black_box(c3))
        })
    });

    let poly = poly_cool::Cubic::new([1.0, -6.0, 11.0, -100.0]);
    group.bench_function("roots_blinn with one root", |b| {
        b.iter(|| black_box(poly).roots_blinn())
    });
}

criterion_group!(benches, cubic_roots);
criterion_main!(benches);
