use arrayvec::ArrayVec;

use crate::Quadratic;

impl Quadratic {
    pub fn roots(&self) -> ArrayVec<f64, 2> {
        let &[c, b, a] = self.coeffs();
        let disc = b * b - 4.0 * a * c;
        if disc.is_finite() {
            let mut ret = ArrayVec::new();
            let mut push = |r: f64| {
                if r.is_finite() {
                    ret.push(r)
                }
            };
            if disc > 0.0 {
                let q = -0.5 * (b + disc.sqrt().copysign(b));
                let r0 = q / a;
                let r1 = c / q;
                push(r0.min(r1));
                push(r0.max(r1));
            } else if disc == 0.0 {
                let root = -0.5 * b / a;
                if root.is_finite() {
                    push(root);
                } else if c == 0.0 {
                    // This is kurbo's behavior: the intention is that if the
                    // whole thing is zero, return zero as a single root. I'm
                    // not sure I love it.
                    //
                    // Bear in mind that this branch is not *only* for the
                    // identically zero case: if a == c == 0.0 and b * b
                    // underflows then we will end up here. In that case,
                    // zero is the only root.
                    push(0.0);
                }
            } else {
                // No roots.
            }
            ret
        } else {
            // At least one of the coefficients was too large and triggered
            // overflow.
            //
            // The exponent of f64 maxes out at 1023, so scaling down by
            // 2^{-512} is enough to ensure that squaring doesn't overflow. We
            // do an extra factor of 2^{-3} for some wiggle room. This can't
            // completely destroy all the coefficients: because of the overflow,
            // we know that at least one of them was big.
            let scale = 2.0f64.powi(-515);
            // TODO: this can stack overflow if we're infinite. How should
            // we handle that?
            (*self * scale).roots()
        }
    }

    pub fn positive_discriminant_roots(&self) -> Option<(f64, f64)> {
        let &[c, b, a] = self.coeffs();
        let disc = b * b - 4.0 * a * c;
        if disc.is_finite() {
            if disc > 0.0 {
                let q = -0.5 * (b + disc.sqrt().copysign(b));
                let r0 = q / a;
                let r1 = c / q;
                Some((r0.min(r1), r0.max(r1)))
            } else {
                None
            }
        } else {
            self.positive_discriminant_roots_scaled()
        }
    }

    #[cold]
    fn positive_discriminant_roots_scaled(&self) -> Option<(f64, f64)> {
        if self.is_finite() {
            let scale = 2.0f64.powi(-515);
            (*self * scale).positive_discriminant_roots()
        } else {
            None
        }
    }

    pub fn positive_discriminant_roots_no_overflow_check(&self) -> Option<(f64, f64)> {
        let &[c, b, a] = self.coeffs();
        let disc = b * b - 4.0 * a * c;
        if disc > 0.0 {
            let q = -0.5 * (b + disc.sqrt().copysign(b));
            let r0 = q / a;
            let r1 = c / q;
            Some((r0.min(r1), r0.max(r1)))
        } else {
            None
        }
    }

    pub fn positive_discriminant_roots_no_overflow_check_half_b(&self) -> Option<(f64, f64)> {
        let &[c, b, a] = self.coeffs();
        let disc = b * b - a * c;
        if disc > 0.0 {
            let q = -(b + disc.sqrt().copysign(b));
            let r0 = q / a;
            let r1 = c / q;
            Some((r0.min(r1), r0.max(r1)))
        } else {
            None
        }
    }
}

#[cfg(test)]
mod tests {
    #[test]
    fn root_evaluation() {
        arbtest::arbtest(|u| {
            let q = crate::arbitrary::quadratic(u)?;
            // Arbitrary quadratics can have coefficients with wild magnitudes,
            // so we need to adjust our error expectations accordingly.
            let magnitude = q.magnitude().max(1.0);

            for r in q.roots() {
                let y = q.eval(r);
                // To evaluate the polynomial, we need to square r, so our error
                // should be relative to the magnitude of r squared.
                let r_magnitude = r.abs().max(1.0);
                let threshold = r_magnitude * 1e-14 * r_magnitude * magnitude;
                if y.is_finite() && threshold.is_finite() {
                    assert!(y.abs() <= threshold);
                }
            }
            Ok(())
        })
        .budget_ms(5_000);
    }
}
