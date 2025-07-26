use arrayvec::ArrayVec;

use crate::Quadratic;

#[derive(Debug, Copy, Clone)]
pub struct Cubic {
    pub c0: f64,
    pub c1: f64,
    pub c2: f64,
    pub c3: f64,
}

struct TerminationData {
    step: f64,
    value: f64,
}

fn different_signs(x: f64, y: f64) -> bool {
    (x < 0.0) == (y < 0.0)
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

    fn critical_points(&self) -> Option<(f64, f64)> {
        let mut a = 3.0 * self.c3;
        let mut b_2 = self.c2;
        let mut c = self.c1;
        let mut disc_4 = b_2 * b_2 - a * c;

        if !disc_4.is_finite() {
            let scale = 2.0f64.powi(-515);
            a *= scale;
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

        // We do one "binary search" step before the truncated Newton
        // iterations. If the range `(lower, upper)` doesn't include any
        // critical points, this guarantees that the Newton method converges
        // (as pointed out by Yuksel): Newton can only oscillate if it crosses
        // an inflection point, and chopping the inter-critical-point range in
        // half guarantees that it doesn't contain an inflection point. (This
        // is specific to cubics: the derivative is a quadratic, and so its
        // extremum is exactly halfway between its roots.)
        let mut x = (upper + lower) / 2.0;
        let mut val_x = self.eval(x);
        let mut step;
        if different_signs(val_x, val_lower) {
            step = x - lower;
            lower = x;
        } else {
            step = upper - x;
            upper = x;
        }

        // TODO: if they ask for unrealistic accuracy, this might loop
        // forever...
        while !term(TerminationData { step, value: val_x }) {
            let deriv_x = self.deriv().eval(x);
            debug_assert!(deriv_x.is_finite());
            debug_assert!(val_x.is_finite());

            step = -val_x / deriv_x;
            let new_x = (x + step).clamp(lower, upper);
            step = new_x - x;
            x = new_x;
            val_x = self.eval(x);
            //dbg!(step, x, val_x);
        }
        x
    }

    pub fn root_between_with_output_error(&self, lower: f64, upper: f64, y_error: f64) -> f64 {
        self.one_root(lower, upper, |TerminationData { value, .. }| {
            value.abs() <= y_error
        })
    }

    pub fn root_between(&self, lower: f64, upper: f64, x_error: f64) -> f64 {
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
    pub fn roots_between_value_accuracy(
        &self,
        lower: f64,
        upper: f64,
        accuracy: f64,
    ) -> ArrayVec<f64, 3> {
        let q = self.deriv();
        if !q.is_finite() {
            // TODO: what's the right response here?
            return ArrayVec::new();
        }

        let mut possible_endpoints = ArrayVec::<f64, 3>::new();
        possible_endpoints.extend(q.roots());
        possible_endpoints.push(upper);

        let mut last = lower;
        let mut last_sign = self.eval(last).signum();
        let mut ret = ArrayVec::new();

        for x in possible_endpoints {
            if x > last && x <= upper {
                let sign = self.eval(x).signum();
                if sign != last_sign {
                    ret.push(self.root_between_with_output_error(last, x, accuracy));
                }

                last = x;
                last_sign = sign;
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
            // Arbitrary cubics can have coefficients with wild magnitudes,
            // so we need to adjust our error expectations accordingly.
            let magnitude =
                c.c0.abs()
                    .max(c.c1.abs())
                    .max(c.c2.abs())
                    .max(c.c3.abs())
                    .max(1.0);
            let accuracy = magnitude * 1e-12;

            // We could have a wider range of roots, but then we might need
            // to lower the accuracy depending on what the actual root is: the
            // intermediate computations scale like the cube of the root.
            for r in c.roots_between_value_accuracy(-10.0, 10.0, accuracy) {
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
