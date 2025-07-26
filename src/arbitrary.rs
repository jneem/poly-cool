//! Utilities for fuzz and/or property testing using `arbitrary`.

use arbitrary::Unstructured;

use crate::{Cubic, Quadratic};

fn check_finite(f: f64) -> Result<f64, arbitrary::Error> {
    if f.is_finite() {
        Ok(f)
    } else {
        Err(arbitrary::Error::IncorrectFormat)
    }
}

pub fn finite_float(u: &mut Unstructured<'_>) -> Result<f64, arbitrary::Error> {
    check_finite(u.arbitrary()?)
}

/// Generate a float, but give it a chance to be close to another float.
fn another_finite_float(orig: f64, u: &mut Unstructured<'_>) -> Result<f64, arbitrary::Error> {
    let close: bool = u.arbitrary()?;
    if close {
        let ulps: i32 = u.int_in_range(-32..=32)?;
        let scale = 1.0f64 + ulps as f64 * f64::EPSILON;
        check_finite(orig * scale)
    } else {
        finite_float(u)
    }
}

/// Generate an arbitrary quadratic polynomial.
pub fn quadratic(u: &mut Unstructured<'_>) -> Result<Quadratic, arbitrary::Error> {
    let use_coeffs: bool = u.arbitrary()?;
    if use_coeffs {
        let c2 = finite_float(u)?;
        let c1 = another_finite_float(c2, u)?;
        let c0 = another_finite_float(c1, u)?;

        Ok(Quadratic { c2, c1, c0 })
    } else {
        let r1 = finite_float(u)?;
        let r2 = another_finite_float(r1, u)?;
        let scale = finite_float(u)?;

        Ok(Quadratic {
            c2: scale,
            c1: check_finite(-scale * (r1 + r2))?,
            c0: check_finite(scale * r1 * r2)?,
        })
    }
}

/// Generate an arbitrary cubic polynomial.
pub fn cubic(u: &mut Unstructured<'_>) -> Result<Cubic, arbitrary::Error> {
    let use_coeffs: bool = u.arbitrary()?;
    if use_coeffs {
        let c3 = finite_float(u)?;
        let c2 = another_finite_float(c3, u)?;
        let c1 = another_finite_float(c2, u)?;
        let c0 = another_finite_float(c1, u)?;

        Ok(Cubic { c3, c2, c1, c0 })
    } else {
        // Generate the roots, with a bias towards roots being almost-repeated.
        let r1 = finite_float(u)?;
        let r2 = another_finite_float(r1, u)?;
        let r3 = another_finite_float(r2, u)?;
        let scale = finite_float(u)?;

        Ok(Cubic {
            c3: scale,
            c2: check_finite(-scale * (r1 + r2 + r3))?,
            c1: check_finite(scale * (r1 * r2 + r1 * r3 + r2 * r3))?,
            c0: check_finite(-scale * r1 * r2 * r3)?,
        })
    }
}
