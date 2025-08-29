# poly-cool

[![GitHub Actions CI status.](https://img.shields.io/github/actions/workflow/status/jneem/linesweeper/ci.yml?logo=github&label=CI)](https://github.com/jneem/poly-cool/actions)
[![Latest published version.](https://img.shields.io/crates/v/poly-cool.svg)](https://crates.io/crates/poly-cool)
[![Documentation build status.](https://img.shields.io/docsrs/poly-cool.svg)](https://docs.rs/poly-cool)

A rust crate for numerically finding roots of low-degree polynomials.

```rust
use poly_cool::Poly;

// The polynomial x^3 - 6x^2 + 11x - 6
let p = Poly::new([-6.0, 11.0, -6.0, 1.0]);

dbg!(p.roots_between(-10.0, 10.0, 1e-6));
// [0.9999999999999996, 2.0000000000000018, 2.9999999999999982]
```

Currently, we implement [Yuksel's] iterative solver for finding roots within a
given interval to a specified target accuracy.


[Yuksel's]: https://www.cemyuksel.com/research/polynomials/
