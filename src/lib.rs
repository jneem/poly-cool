//! This is a crate for numerical polynomial root-finding.
//!
//! Currently, we implement a single solver: Yuksel's iterative algorithm
//! for finding roots in a bounded interval. We aspire to have
//! more, with thorough tests and benchmarks.

mod cubic;
mod poly;
mod poly_dyn;
mod quadratic;
mod yuksel;

#[cfg(any(test, feature = "arbitrary"))]
pub mod arbitrary;

#[cfg(any(test, feature = "dashu-float"))]
pub mod accurate;

pub use poly::{Cubic, Poly, Quadratic, Quartic, Quintic};
pub use poly_dyn::PolyDyn;

fn different_signs(x: f64, y: f64) -> bool {
    (x < 0.0) != (y < 0.0)
}
