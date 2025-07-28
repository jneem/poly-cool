use arrayvec::ArrayVec;

use crate::Quadratic;

#[derive(Debug, Copy, Clone)]
pub struct Cubic {
    pub c0: f64,
    pub c1: f64,
    pub c2: f64,
    pub c3: f64,
}

impl std::ops::Div<f64> for Cubic {
    type Output = Cubic;

    fn div(self, rhs: f64) -> Cubic {
        Cubic {
            c0: self.c0 / rhs,
            c1: self.c1 / rhs,
            c2: self.c2 / rhs,
            c3: self.c3 / rhs,
        }
    }
}

struct TerminationData {
    step: f64,
    value: f64,
}

fn different_signs(x: f64, y: f64) -> bool {
    (x < 0.0) != (y < 0.0)
}

impl Cubic {
    pub fn eval(&self, x: f64) -> f64 {
        let xx = x * x;
        let xxx = xx * x;
        self.c0 + self.c1 * x + self.c2 * xx + self.c3 * xxx
    }

    pub fn deriv(&self) -> Quadratic {
        Quadratic {
            c0: self.c1,
            c1: 2.0 * self.c2,
            c2: 3.0 * self.c3,
        }
    }

    pub fn max_coeff(&self) -> f64 {
        self.c0
            .abs()
            .max(self.c1.abs())
            .max(self.c2.abs())
            .max(self.c3.abs())
    }

    fn too_large_for_roots(&self) -> bool {
        self.max_coeff() >= 2.0f64.powi(500)
    }

    fn shrinkage_factor() -> f64 {
        2.0f64.powi(524)
    }

    /// Computes the critical points of this cubic, as long
    /// as the discriminant of the derivative is positive.
    /// The return values are in increasing order.
    ///
    /// Some corner cases worth noting:
    ///   - If the discriminant is zero, returns nothing. That is,
    ///     we don't find double-roots of the derivative.
    ///   - If the derivative is linear or close to it, we might
    ///     return +/- infinity as one of the roots.
    ///   - Unless some input is NaN, we don't return NaN.
    fn critical_points(&self) -> Option<(f64, f64)> {
        let mut a = 3.0 * self.c3;
        let mut b_2 = self.c2;
        let mut c = self.c1;
        let mut disc_4 = b_2 * b_2 - a * c;

        if !disc_4.is_finite() {
            let scale = 2.0f64.powi(-515);
            a = self.c3 * scale * 3.0;
            b_2 *= scale;
            c *= scale;
            disc_4 = b_2 * b_2 - a * c;
            debug_assert!(disc_4.is_finite());
        }

        if disc_4 > 0.0 {
            let q = -(b_2 + disc_4.sqrt().copysign(b_2));
            let r0 = q / a;
            let r1 = c / q;
            Some((r0.min(r1), r0.max(r1)))
        } else {
            None
        }
    }

    fn one_root<Term: Fn(TerminationData) -> bool>(
        &self,
        mut lower: f64,
        mut upper: f64,
        term: Term,
    ) -> f64 {
        let val_lower = self.eval(lower);
        let val_upper = self.eval(upper);
        if !val_lower.is_finite() || !val_upper.is_finite() || !self.deriv().is_finite() {
            return f64::NAN;
        }
        debug_assert!(different_signs(val_lower, val_upper));

        let mut x = lower + (upper - lower) / 2.0;
        let mut val_x = self.eval(x);
        let mut step = (upper - lower) / 2.0;

        while x.is_finite() && !term(TerminationData { step, value: val_x }) {
            let root_in_first_half = different_signs(val_lower, val_x);
            if root_in_first_half {
                upper = x;
            } else {
                lower = x;
            }

            let deriv_x = self.deriv().eval(x);
            debug_assert!(deriv_x.is_finite());
            debug_assert!(val_x.is_finite());

            step = -val_x / deriv_x;
            let mut new_x = x + step;

            if new_x <= lower || new_x >= upper {
                new_x = lower + (upper - lower) / 2.0;

                if new_x == upper || new_x == lower {
                    // This should be rare, but it happens if they ask for more
                    // accuracy than is reasonable. For example, suppse (because
                    // of large coefficients) the output value jumps from -1.0
                    // to 1.0 between adjacent floats and they ask for an output
                    // error of smaller than 0.5. Then we'll eventually shrink
                    // the search interval to a pair of adjacent floats and hit
                    // this case.
                    return new_x;
                }
            }
            step = new_x - x;
            x = new_x;
            val_x = self.eval(x);
        }
        x
    }

    pub fn root_between_with_output_error(self, lower: f64, upper: f64, y_error: f64) -> f64 {
        self.one_root(lower, upper, |TerminationData { value, .. }| {
            value.abs() <= y_error
        })
    }

    pub fn root_between(self, lower: f64, upper: f64, x_error: f64) -> f64 {
        self.one_root(lower, upper, |TerminationData { step, .. }| {
            step.abs() <= x_error
        })
    }

    /// Computes all roots between `lower` and `upper`, to the desired accuracy.
    ///
    /// "Accuracy" is measured with respect to the cubic's value: if this cubic
    /// is called `f` and we find some `x` with `|f(x)| < accuracy` (and `x` is
    /// contained between two endpoints where `f` has opposite signs) then we'll
    /// call `x` a root.
    ///
    /// We make no guarantees about multiplicity. In fact, if there's a
    /// double-root that isn't a triple-root (and therefore has no sign change
    /// nearby) then there's a good chance we miss it altogether. This is
    /// fine if you're using this root-finding to optimize a quartic, because
    /// double-roots of the derivative aren't local extrema.
    pub fn roots_between_with_output_error(
        mut self,
        lower: f64,
        upper: f64,
        mut y_error: f64,
    ) -> ArrayVec<f64, 3> {
        if self.too_large_for_roots() {
            self = self / Cubic::shrinkage_factor();
            y_error /= Cubic::shrinkage_factor();
        }

        let mut possible_endpoints = ArrayVec::<f64, 3>::new();
        if let Some((x0, x1)) = self.critical_points() {
            possible_endpoints.push(x0);
            possible_endpoints.push(x1);
        }
        possible_endpoints.push(upper);

        let mut last = lower;
        let mut last_val = self.eval(last);
        let mut ret = ArrayVec::new();

        for x in possible_endpoints {
            if x > last && x <= upper {
                let val = self.eval(x);
                if different_signs(last_val, val) {
                    ret.push(self.root_between_with_output_error(last, x, y_error));
                }

                last = x;
                last_val = val;
            }
        }
        ret
    }
}

#[cfg(test)]
mod tests {
    #[test]
    fn root_evaluation() {
        arbtest::arbtest(|u| {
            let c = crate::arbitrary::cubic(u)?;
            //dbg!(c);
            // Arbitrary cubics can have coefficients with wild magnitudes,
            // so we need to adjust our error expectations accordingly.
            let magnitude = c.max_coeff().max(1.0);
            let accuracy = magnitude * 1e-12;

            // We could have a wider range of roots, but then we might need
            // to lower the accuracy depending on what the actual root is: the
            // intermediate computations scale like the cube of the root.
            for r in c.roots_between_with_output_error(-10.0, 10.0, accuracy) {
                let y = c.eval(r);
                if y.is_finite() {
                    assert!(y.abs() <= accuracy);
                }
            }
            Ok(())
        })
        .budget_ms(5_000);
    }

    #[test]
    #[ignore]
    fn root_evaluation_kurbo() {
        arbtest::arbtest(|u| {
            let c = crate::arbitrary::cubic(u)?;
            // Arbitrary cubics can have coefficients with wild magnitudes,
            // so we need to adjust our error expectations accordingly.
            let magnitude = c.max_coeff().max(1.0);
            let accuracy = magnitude * 1e-12;

            // We could have a wider range of roots, but then we might need
            // to lower the accuracy depending on what the actual root is: the
            // intermediate computations scale like the cube of the root.
            for r in kurbo::common::solve_cubic(c.c0, c.c1, c.c2, c.c3) {
                let y = c.eval(r);
                if y.is_finite() {
                    assert!(y.abs() <= accuracy);
                }
            }
            Ok(())
        })
        .budget_ms(5_000);
    }
}
