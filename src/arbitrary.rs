//! Utilities for fuzz and/or property testing using `arbitrary`.

use arbitrary::Unstructured;

use crate::{Cubic, Poly, Quadratic};

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

        Ok(Quadratic::new([c0, c1, c2]))
    } else {
        let r1 = finite_float(u)?;
        let r2 = another_finite_float(r1, u)?;
        let scale = finite_float(u)?;

        Ok(Quadratic::new([
            check_finite(scale * r1 * r2)?,
            check_finite(-scale * (r1 + r2))?,
            scale,
        ]))
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

        Ok(Cubic::new([c0, c1, c2, c3]))
    } else {
        // Generate the roots, with a bias towards roots being almost-repeated.
        let r1 = finite_float(u)?;
        let r2 = another_finite_float(r1, u)?;
        let r3 = another_finite_float(r2, u)?;
        let scale = finite_float(u)?;

        Ok(Cubic::new([
            check_finite(-scale * r1 * r2 * r3)?,
            check_finite(scale * (r1 * r2 + r1 * r3 + r2 * r3))?,
            check_finite(-scale * (r1 + r2 + r3))?,
            scale,
        ]))
    }
}

pub fn poly<const N: usize>(u: &mut Unstructured<'_>) -> Result<Poly<N>, arbitrary::Error> {
    assert!(N >= 2);

    let use_coeffs: bool = u.arbitrary()?;
    if use_coeffs {
        let mut coeffs = [0.0; N];
        coeffs[0] = finite_float(u)?;
        for i in 1..coeffs.len() {
            coeffs[i] = another_finite_float(coeffs[i - 1], u)?;
        }

        Ok(Poly::new(coeffs))
    } else {
        // Generate the roots, with a bias towards roots being almost-repeated.
        let mut r = finite_float(u)?;
        let mut coeffs = [0.0; N];
        coeffs[1] = 1.0;
        coeffs[0] = -r;

        for _ in 0..(N - 2) {
            r = another_finite_float(r, u)?;
            mul(&mut coeffs, r);
        }

        let scale = finite_float(u)?;
        for c in &mut coeffs {
            *c *= scale;
            check_finite(*c)?;
        }

        Ok(Poly::new(coeffs))
    }
}

// Takes the polynomial in `coeffs` and multiplies it by (x - root).
// (Only correct if the last coefficient is zero.)
fn mul<const N: usize>(coeffs: &mut [f64; N], root: f64) {
    for i in (1..coeffs.len()).rev() {
        coeffs[i] = coeffs[i - 1] - root * coeffs[i];
    }
}
