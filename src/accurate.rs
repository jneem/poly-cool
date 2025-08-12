use dashu_float::{
    FBig,
    ops::{Abs, SquareRoot},
    round::{Round, Rounded, mode::Up},
};

use crate::{Cubic, Poly};

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

// Overflow-proof midpoint between two numbers.
fn midpoint(x: f64, y: f64) -> f64 {
    if (x > 0.0) == (y > 0.0) {
        x + (y - x) / 2.0
    } else {
        (x + y) / 2.0
    }
}

fn round<RoundingMode: Round>(x: FBig<RoundingMode>, precision: usize) -> FBig<RoundingMode> {
    // Work around the bug with arbitrary precision (i.e. precision == 0)
    // numbers don't get rounded by `with_precision`.
    x.with_precision(precision + 1)
        .value()
        .with_precision(precision)
        .value()
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

    /// Evaluate this polynomial at the given point, as a correctly-rounded
    /// `f64`.
    pub fn eval(&self, x: f64) -> f64 {
        self.eval_exact(x).to_f64().value()
    }

    /// Evaluate this polynomial at the given point, exactly.
    pub fn eval_exact(&self, x: f64) -> FBig {
        let x = exact(x);
        let mut power = exact(1.0);
        let mut ret = exact(0.0);

        // We're doing this in arbitrary precision, which is hilariously
        // inefficient if the coefficients have very different magnitudes.
        // But it's obviously correct, and if we were constructed from `f64`s
        // (which have only 11 exponent bits) then we're wasting a few kilobytes
        // at most.
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
        match self.coeffs.len() {
            0..2 => Vec::new(),
            2 => {
                let root = -round(self.coeffs[0].clone(), 128) / round(self.coeffs[1].clone(), 256);
                if finite_f64(&root) {
                    vec![root.to_f64().value()]
                } else {
                    Vec::new()
                }
            }
            3 => {
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
            }
            _ => {
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
                    if x1 <= x0 {
                        continue;
                    }

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
    }

    // Returns a number x such that all roots of this polynomial
    // must lie in the interval [-x, x].
    //
    // TODO: also, make sure the result is small enough so that the polynomial
    // is finite (in f64) on the range.
    fn effective_domain_bounds(&self) -> FBig {
        let mut x = exact(0.0);
        let factor = exact(self.coeffs.len().saturating_sub(1) as f64).with_rounding();
        let Some(highest_order_coeff) = self.coeffs.last().cloned() else {
            return x;
        };
        let highest_order_coeff = highest_order_coeff.abs().with_rounding();

        for coeff in self.coeffs.iter().rev().skip(1) {
            let coeff: FBig<Up> = round(coeff.clone().abs().with_rounding(), 128);
            let bound = (coeff / &highest_order_coeff) * &factor;
            x = x.max(bound.with_rounding());
        }
        x
    }

    // Given a bracketing interval, returns a root.
    //
    // The return value is guaranteed to either be a root (in the sense that
    // correctly-rounded exact evaluation of the polynomial gives zero) or be
    // just below a root (in the sense that evaluation of the returned value
    // has a different sign as evaluation of the next float larger tha it).
    fn one_root(&self, deriv: &AccuPoly, mut lower: f64, mut upper: f64) -> f64 {
        debug_assert!(lower < upper);
        let val_lower = self.eval_exact(lower);
        let val_upper = self.eval_exact(upper);

        debug_assert!(val_lower.repr().is_finite());
        debug_assert!(val_upper.repr().is_finite());
        debug_assert!(val_lower.sign() != val_upper.sign());

        let mut x = midpoint(lower, upper);
        let mut val_x = self.eval_exact(x);
        let zero = exact(0.0);

        assert!(x >= lower);
        assert!(x <= upper);
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
            let step = (-round(val_x.clone(), 128) / round(deriv_x.clone(), 256))
                .with_precision(0)
                .value();
            let mut new_x = (exact(x) + &step).to_f64().value();

            if new_x <= lower || new_x >= upper {
                new_x = midpoint(lower, upper);
            }
            if new_x == x {
                break;
            }
            x = new_x;
            val_x = self.eval_exact(x);
        }

        let mut val_x = val_x.to_f64().value();
        // dbg!(lower, x, upper);
        // dbg!(val_x, &val_lower, &val_upper);
        if val_x != 0.0 {
            // Do a linear search among `f64`s to find the exact root.
            // FIXME: sometimes this takes too long. I'm not sure why the
            // Newton stops...
            let root_in_first_half = (val_lower > zero) != (val_x > 0.0);
            //dbg!(root_in_first_half);
            if root_in_first_half {
                let mut prev_x = x.next_down();
                let mut val_prev_x = self.eval(prev_x);
                //dbg!(x, val_x, prev_x, val_prev_x);
                let mut count = 0;
                while val_prev_x != 0.0 && val_prev_x.signum() == val_x.signum() {
                    count += 1;
                    if count > 32 {
                        panic!();
                    }
                    dbg!(x, val_x, prev_x, val_prev_x);
                    x = prev_x;
                    val_x = val_prev_x;
                    prev_x = x.next_down();
                    val_prev_x = self.eval(prev_x);
                }
                return prev_x;
            } else {
                let mut next_x = x.next_up();
                let mut val_next_x = self.eval(next_x);
                let mut count = 0;
                while val_x != 0.0 && val_next_x.signum() == val_x.signum() {
                    //dbg!(x, val_x, next_x, val_next_x);
                    count += 1;
                    if count > 32 {
                        panic!();
                    }
                    x = next_x;
                    val_x = val_next_x;
                    next_x = x.next_up();
                    val_next_x = self.eval(next_x)
                }
                if val_x != 0.0 && val_next_x == 0.0 {
                    return next_x;
                }
                return x;
            }
        }
        dbg!(x)
    }
}

impl From<Poly> for AccuPoly {
    fn from(p: Poly) -> AccuPoly {
        AccuPoly::new(p.coeffs().iter().copied())
    }
}

impl From<Cubic> for AccuPoly {
    fn from(p: Cubic) -> AccuPoly {
        AccuPoly::new([p.c0, p.c1, p.c2, p.c3])
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn linear_roots_once() {
        let x_minus_one = AccuPoly::new([-1.0, 1.0]);
        assert_eq!(x_minus_one.roots(), vec![1.0]);
    }

    #[test]
    fn quadratic_roots_once() {
        let x_minus_one = AccuPoly::new([-1.0, 1.0]);
        let x_minus_two = AccuPoly::new([-2.0, 1.0]);
        let p = &x_minus_one * &x_minus_two;
        assert_eq!(p.roots(), vec![1.0, 2.0]);

        let p = &x_minus_two * &x_minus_two;
        assert_eq!(p.roots(), vec![]);
    }

    #[test]
    fn cubic_roots_once() {
        let x_minus_one = AccuPoly::new([-1.0, 1.0]);
        let x_minus_two = AccuPoly::new([-2.0, 1.0]);
        let x_minus_three = AccuPoly::new([-3.0, 1.0]);
        let p = &(&x_minus_one * &x_minus_two) * &x_minus_three;
        assert_eq!(p.roots(), vec![1.0, 2.0, 3.0]);
    }

    #[test]
    fn quadric_roots_once() {
        let x_minus_one = AccuPoly::new([-1.0, 1.0]);
        let x_minus_two = AccuPoly::new([-2.0, 1.0]);
        let x_minus_three = AccuPoly::new([-3.0, 1.0]);
        let x_minus_four = AccuPoly::new([-4.0, 1.0]);
        let p = &x_minus_one * &x_minus_two * &x_minus_three * &x_minus_four;
        assert_eq!(p.roots(), vec![1.0, 2.0, 3.0, 4.0]);
    }

    #[test]
    fn cubic_roots() {
        arbtest::arbtest(|u| {
            let c: AccuPoly = crate::arbitrary::cubic(u)?.into();
            for r in c.roots() {
                let val = c.eval(r);
                let next_val = c.eval(r.next_up());
                //dbg!(&c, r, val, next_val);
                assert!(val == 0.0 || (next_val != 0.0 && val.signum() != next_val.signum()));
            }
            Ok(())
        })
        //.seed(0xfe5fea2200000028);
        //.seed(0x76b66cbb00000020);
        //.seed(0x3034152d00003215);
        //.seed(0xdb4f193500000060);
        //.seed(0xe2a4a0d300010000);
        //.seed(0xd440507e00010000);
        .budget_ms(10_000);
    }
}
