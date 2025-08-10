use dashu_float::{
    FBig,
    ops::{Abs, SquareRoot},
    round::{Rounded, mode::Up},
};

#[derive(Clone, Debug)]
pub struct AccuPoly {
    coeffs: Vec<FBig>,
}

impl<'a> std::ops::Mul<&'a AccuPoly> for &'a AccuPoly {
    type Output = AccuPoly;

    fn mul(self, rhs: &'a AccuPoly) -> Self::Output {
        let mut coeffs = vec![exact(0.0); (self.coeffs.len() + rhs.coeffs.len()).saturating_sub(1)];

        for (i, c) in self.coeffs.iter().enumerate() {
            for (j, d) in rhs.coeffs.iter().enumerate() {
                coeffs[i + j] += c * d;
            }
        }
        AccuPoly { coeffs }
    }
}

impl<'a> std::ops::Mul<&'a AccuPoly> for AccuPoly {
    type Output = AccuPoly;

    fn mul(self, rhs: &'a AccuPoly) -> Self::Output {
        (&self) * rhs
    }
}

fn exact(x: f64) -> FBig {
    let Rounded::Exact(x) = FBig::try_from(x).unwrap().with_precision(0) else {
        unreachable!()
    };
    x
}

fn finite_f64(x: &FBig) -> bool {
    let max = exact(f64::MAX);
    let min = exact(f64::MIN);
    &min <= x && x <= &max
}

fn quadratic_roots(a: &FBig, b: &FBig, c: &FBig) -> Option<(FBig, FBig)> {
    // We need to take the square root of the discriminant, so it has to have
    // finite precision. If we wanted to be real fancy here, we'd check the
    // rounding error and increase the precision until we're sure that our
    // answer is correct to the nearest `f64`. But these are only intermediate
    // values for the algorithm, so it shouldn't really matter.
    let disc = (b.sqr() - a * c * exact(4.0)).with_precision(128).value();

    let zero = exact(0.0);
    if disc > zero {
        let mut sqrt_disc = disc.sqrt();
        if b < &zero {
            sqrt_disc = -sqrt_disc;
        }
        let q = -(b + sqrt_disc) / exact(2.0);
        let r0 = &q / a;
        let r1 = c / q;
        if r0 > r1 {
            Some((r1, r0))
        } else {
            Some((r0, r1))
        }
    } else {
        None
    }
}

impl AccuPoly {
    pub fn new(coeffs: impl IntoIterator<Item = f64>) -> Self {
        let mut ret = Self {
            coeffs: coeffs.into_iter().map(exact).collect(),
        };
        ret.remove_zeros();
        ret
    }

    fn remove_zeros(&mut self) {
        let zero = exact(0.0);
        while let Some(last) = self.coeffs.last() {
            if last == &zero {
                self.coeffs.pop();
            } else {
                return;
            }
        }
    }

    pub fn eval(&self, x: f64) -> f64 {
        self.eval_exact(x).to_f64().value()
    }

    pub fn eval_exact(&self, x: f64) -> FBig {
        let x = exact(x);
        let mut power = exact(1.0);
        let mut ret = exact(0.0);

        for c in &self.coeffs {
            ret += &power * c;
            power *= &x;
        }
        ret
    }

    pub fn deriv(&self) -> Self {
        let mut coeffs = Vec::new();
        for (i, c) in self.coeffs.iter().enumerate().skip(1) {
            coeffs.push(c.clone() * exact(i as f64));
        }
        Self { coeffs }
    }

    pub fn roots(&self) -> Vec<f64> {
        if self.coeffs.len() < 3 {
            unimplemented!()
        } else if self.coeffs.len() == 3 {
            if let Some((r0, r1)) =
                quadratic_roots(&self.coeffs[2], &self.coeffs[1], &self.coeffs[0])
            {
                let mut ret = Vec::new();
                if finite_f64(&r0) {
                    ret.push(r0.to_f64().value());
                }
                if finite_f64(&r1) {
                    ret.push(r1.to_f64().value());
                }
                ret
            } else {
                Vec::new()
            }
        } else {
            let bound = self.effective_domain_bounds();
            let (lower, upper) = if finite_f64(&bound) {
                (-bound.to_f64().value(), bound.to_f64().value())
            } else {
                (f64::MIN, f64::MAX)
            };
            let deriv = self.deriv();
            let mut crits: Vec<f64> = deriv.roots();
            crits.push(upper);
            let mut x0 = lower;
            let mut x0_val = self.eval_exact(x0);
            let mut out = Vec::new();
            for x1 in crits {
                let x1_val = self.eval_exact(x1);
                if x0_val.sign() != x1_val.sign() {
                    out.push(self.one_root(&deriv, x0, x1));
                }
                x0 = x1;
                x0_val = x1_val;
            }
            out
        }
    }

    // Returns a number x such that all roots of this polynomial
    // must lie in the interval [-x, x].
    fn effective_domain_bounds(&self) -> FBig {
        let mut x = exact(0.0);
        let factor = exact(self.coeffs.len().saturating_sub(1) as f64).with_rounding();
        let Some(highest_order_coeff) = self.coeffs.last().cloned() else {
            return x;
        };
        let highest_order_coeff = highest_order_coeff.abs().with_rounding();

        for coeff in self.coeffs.iter().rev().skip(1) {
            let coeff: FBig<Up> = coeff
                .clone()
                .abs()
                .with_rounding()
                .with_precision(128)
                .value();
            let bound = (coeff / &highest_order_coeff) * &factor;
            x = x.max(bound.with_rounding());
        }
        x
    }

    fn one_root(&self, deriv: &AccuPoly, mut lower: f64, mut upper: f64) -> f64 {
        let val_lower = self.eval_exact(lower);
        let val_upper = self.eval_exact(upper);

        debug_assert!(val_lower.repr().is_finite());
        debug_assert!(val_upper.repr().is_finite());
        debug_assert!(val_lower.sign() != val_upper.sign());

        let mut x = lower + (upper - lower) / 2.0;
        let mut val_x = self.eval_exact(x);
        let zero = exact(0.0);

        while val_x != zero {
            let root_in_first_half = val_lower.sign() != val_x.sign();
            if root_in_first_half {
                upper = x;
            } else {
                lower = x;
            }

            let deriv_x = deriv.eval_exact(x);
            debug_assert!(deriv_x.repr().is_finite());
            debug_assert!(val_x.repr().is_finite());

            // TODO: maybe rational is better than float, since we need to divide
            // in a few places? We'd have to figure out what to do with the sqrt
            // in the quadratic formula...
            // TODO: unsure how to choose the precision here. If we don't have enough, there
            // are assertion errors in dashu-float...
            let step = (-val_x.with_precision(128).value() / deriv_x.with_precision(256).value())
                .with_precision(0)
                .value();
            let mut new_x = (exact(x) + &step).to_f64().value();

            if new_x <= lower || new_x >= upper {
                new_x = lower + (upper - lower) / 2.0;
            }
            if new_x == x {
                return new_x;
            }
            x = new_x;
            val_x = self.eval_exact(x);
        }
        x
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn quadratic_roots() {
        let x_minus_one = AccuPoly::new([-1.0, 1.0]);
        let x_minus_two = AccuPoly::new([-2.0, 1.0]);
        let p = &x_minus_one * &x_minus_two;
        assert_eq!(p.roots(), vec![1.0, 2.0]);

        let p = &x_minus_two * &x_minus_two;
        assert_eq!(p.roots(), vec![]);
    }

    #[test]
    fn cubic_roots() {
        let x_minus_one = AccuPoly::new([-1.0, 1.0]);
        let x_minus_two = AccuPoly::new([-2.0, 1.0]);
        let x_minus_three = AccuPoly::new([-3.0, 1.0]);
        let p = &(&x_minus_one * &x_minus_two) * &x_minus_three;
        assert_eq!(p.roots(), vec![1.0, 2.0, 3.0]);
    }

    #[test]
    fn quadric_roots() {
        let x_minus_one = AccuPoly::new([-1.0, 1.0]);
        let x_minus_two = AccuPoly::new([-2.0, 1.0]);
        let x_minus_three = AccuPoly::new([-3.0, 1.0]);
        let x_minus_four = AccuPoly::new([-4.0, 1.0]);
        let p = &x_minus_one * &x_minus_two * &x_minus_three * &x_minus_four;
        assert_eq!(p.roots(), vec![1.0, 2.0, 3.0, 4.0]);
    }
}
