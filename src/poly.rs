use arrayvec::ArrayVec;

/// Polynomial multiplication is not yet implemented, because doing it "nicely"
/// would require const generic expressions: ideally we'd do something like
///
/// ```ignore
/// impl<N, M> Mul<Poly<M>> for Poly<N> {
///     type Output = Poly<{M + N - 1}>;
/// }
/// ```
///
/// It's possible to work around this with macros, but there are lots of
/// possibilities and I didn't feel like it was worth the trouble (and the hit
/// to compilation time).
#[derive(Clone, Copy, Debug, PartialEq, PartialOrd)]
pub struct Poly<const N: usize> {
    pub(crate) coeffs: [f64; N],
}

pub type Quadratic = Poly<3>;
pub type Cubic = Poly<4>;
pub type Quartic = Poly<5>;
pub type Quintic = Poly<6>;

impl<const N: usize> Poly<N> {
    /// Creates a new polynomial with the provided coefficients.
    ///
    /// The constant coefficient comes first, then the linear coefficient, and
    /// so on. So if you pass `[c, b, a]` you'll get the polynomial
    /// `a x^2 + b x + c`.
    pub const fn new(coeffs: [f64; N]) -> Poly<N> {
        Poly { coeffs }
    }

    /// The coefficients of this polynomial.
    ///
    /// In the returned array, the coefficient of `x^i` is at index `i`.
    pub fn coeffs(&self) -> &[f64; N] {
        &self.coeffs
    }

    /// Evaluates this polynomial at a point.
    pub fn eval(&self, x: f64) -> f64 {
        let mut acc = 0.0;
        for c in self.coeffs.iter().rev() {
            // It would be nice to use `f64::mul_add` here, but it's slow on
            // architectures that don't have a dedicated instruction.
            acc = acc * x + c;
        }
        acc
    }

    /// Returns the largest absolute value of any coefficient.
    ///
    /// Always returns a non-negative number, or NaN if some coefficient is NaN.
    pub fn magnitude(&self) -> f64 {
        let mut max = 0.0f64;
        for c in &self.coeffs {
            max = max.max(c.abs());
        }
        max
    }

    /// Are all the coefficients finite?
    pub fn is_finite(&self) -> bool {
        self.coeffs.iter().all(|c| c.is_finite())
    }
}

macro_rules! impl_deriv_and_deflate {
    ($N:literal, $N_MINUS_ONE:literal) => {
        impl Poly<$N> {
            /// Compute the derivative of this polynomial, as a polynomial with
            /// one less coefficient.
            pub fn deriv(&self) -> Poly<$N_MINUS_ONE> {
                let mut coeffs = [0.0; $N_MINUS_ONE];
                for (i, (d, c)) in coeffs.iter_mut().zip(&self.coeffs[1..]).enumerate() {
                    *d = (i + 1) as f64 * c;
                }
                Poly::new(coeffs)
            }

            /// Divide this polynomial by the polynomial `x - root`, returning the
            /// quotient (as a polynomial with one less coefficient) and ignoring
            /// the remainder.
            ///
            /// If `root` is actually a root of `self` (as the name suggests
            /// it should be, but this is not actually required), the
            /// remainder will be zero. In general, the remainder will be
            /// `self.eval(root)`.
            pub fn deflate(&self, root: f64) -> Poly<$N_MINUS_ONE> {
                let mut acc = 0.0;
                let mut coeffs = [0.0; $N_MINUS_ONE];
                for (d, c) in coeffs.iter_mut().zip(&self.coeffs[1..]).rev() {
                    acc = acc * root + c;
                    *d = acc;
                }
                Poly::new(coeffs)
            }
        }
    };
}

macro_rules! impl_roots_between_recursive {
    ($N:literal, $N_MINUS_ONE:literal) => {
        impl Poly<$N> {
            pub fn roots_between(
                self,
                lower: f64,
                upper: f64,
                x_error: f64,
            ) -> ArrayVec<f64, $N_MINUS_ONE> {
                let mut ret = ArrayVec::new();
                let mut scratch = ArrayVec::new();
                self.roots_between_with_buffer(lower, upper, x_error, &mut scratch, &mut ret);
                ret
            }

            // This would ideally have a `where M >= N - 1` bound on it,
            // but it's private so it shouldn't matter too much.
            // We assume that `scratch` and `out` are both empty.
            fn roots_between_with_buffer<const M: usize>(
                self,
                lower: f64,
                upper: f64,
                x_error: f64,
                scratch: &mut ArrayVec<f64, M>,
                out: &mut ArrayVec<f64, M>,
            ) {
                let deriv = self.deriv();
                if !deriv.is_finite() {
                    return;
                }
                deriv.roots_between_with_buffer(lower, upper, x_error, out, scratch);
                scratch.push(upper);
                out.clear();
                let mut last = lower;
                let mut last_val = self.eval(last);

                // `endpoint` now contains all the critical points (in increasing order)
                // and the upper endpoint of the interval. These are the endpoints
                // of the potential bracketing intervals of our polynomial.
                for &mut x in scratch {
                    let val = self.eval(x);
                    if $crate::different_signs(last_val, val) {
                        out.push($crate::yuksel::find_root(
                            |x| self.eval(x),
                            |x| deriv.eval(x),
                            last,
                            x,
                            last_val,
                            val,
                            x_error,
                        ));
                    }

                    last = x;
                    last_val = val;
                }
            }
        }
    };
}

impl_deriv_and_deflate!(3, 2);
impl_deriv_and_deflate!(4, 3);
impl_deriv_and_deflate!(5, 4);
impl_deriv_and_deflate!(6, 5);
impl_deriv_and_deflate!(7, 6);

impl_roots_between_recursive!(5, 4);
impl_roots_between_recursive!(6, 5);
impl_roots_between_recursive!(7, 6);

impl<const N: usize> std::ops::Mul<f64> for Poly<N> {
    type Output = Poly<N>;

    fn mul(mut self, scale: f64) -> Poly<N> {
        self *= scale;
        self
    }
}

impl<const N: usize> std::ops::MulAssign<f64> for Poly<N> {
    fn mul_assign(&mut self, scale: f64) {
        for c in &mut self.coeffs {
            *c *= scale;
        }
    }
}

impl<const N: usize> std::ops::Mul<f64> for &Poly<N> {
    type Output = Poly<N>;

    fn mul(self, scale: f64) -> Poly<N> {
        (*self) * scale
    }
}

impl<const N: usize> std::ops::Div<f64> for Poly<N> {
    type Output = Poly<N>;

    fn div(mut self, scale: f64) -> Poly<N> {
        self /= scale;
        self
    }
}

impl<const N: usize> std::ops::DivAssign<f64> for Poly<N> {
    fn div_assign(&mut self, scale: f64) {
        for c in &mut self.coeffs {
            *c /= scale;
        }
    }
}

impl<const N: usize> std::ops::Div<f64> for &Poly<N> {
    type Output = Poly<N>;

    fn div(self, scale: f64) -> Poly<N> {
        (*self) / scale
    }
}

impl<const N: usize> std::ops::AddAssign<&Poly<N>> for Poly<N> {
    fn add_assign(&mut self, rhs: &Poly<N>) {
        for (c, d) in self.coeffs.iter_mut().zip(rhs.coeffs) {
            *c += d;
        }
    }
}

impl<const N: usize> std::ops::AddAssign<Poly<N>> for Poly<N> {
    fn add_assign(&mut self, rhs: Poly<N>) {
        *self += &rhs;
    }
}

impl<const N: usize> std::ops::Add<Poly<N>> for Poly<N> {
    type Output = Poly<N>;

    fn add(mut self, rhs: Poly<N>) -> Poly<N> {
        self += rhs;
        self
    }
}

impl<const N: usize> std::ops::Add<&Poly<N>> for Poly<N> {
    type Output = Poly<N>;

    fn add(mut self, rhs: &Poly<N>) -> Poly<N> {
        self += rhs;
        self
    }
}

impl<const N: usize> std::ops::Add<Poly<N>> for &Poly<N> {
    type Output = Poly<N>;

    fn add(self, mut rhs: Poly<N>) -> Poly<N> {
        rhs += self;
        rhs
    }
}

impl<const N: usize> std::ops::SubAssign<&Poly<N>> for Poly<N> {
    fn sub_assign(&mut self, rhs: &Poly<N>) {
        for (c, d) in self.coeffs.iter_mut().zip(rhs.coeffs) {
            *c -= d;
        }
    }
}

impl<const N: usize> std::ops::SubAssign<Poly<N>> for Poly<N> {
    fn sub_assign(&mut self, rhs: Poly<N>) {
        *self -= &rhs;
    }
}

impl<const N: usize> std::ops::Sub<Poly<N>> for Poly<N> {
    type Output = Poly<N>;

    fn sub(mut self, rhs: Poly<N>) -> Poly<N> {
        self -= rhs;
        self
    }
}

impl<const N: usize> std::ops::Sub<&Poly<N>> for Poly<N> {
    type Output = Poly<N>;

    fn sub(mut self, rhs: &Poly<N>) -> Poly<N> {
        self -= rhs;
        self
    }
}

impl<const N: usize> std::ops::Sub<Poly<N>> for &Poly<N> {
    type Output = Poly<N>;

    fn sub(self, mut rhs: Poly<N>) -> Poly<N> {
        rhs -= self;
        rhs
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn smoke() {
        let p = Poly::new([-6.0, 11.0, -6.0, 1.0]);

        let roots = p.roots_between(0.0, 5.0, 1e-6);
        assert_eq!(roots.len(), 3);
        assert!((roots[0] - 1.0).abs() <= 1e-6);
        assert!((roots[1] - 2.0).abs() <= 1e-6);
        assert!((roots[2] - 3.0).abs() <= 1e-6);

        let p = Poly::new([24.0, -50.0, 35.0, -10.0, 1.0]);

        let roots = p.roots_between(0.0, 5.0, 1e-6);
        assert_eq!(roots.len(), 4);
        assert!((roots[0] - 1.0).abs() <= 1e-6);
        assert!((roots[1] - 2.0).abs() <= 1e-6);
        assert!((roots[2] - 3.0).abs() <= 1e-6);
        assert!((roots[3] - 4.0).abs() <= 1e-6);
    }
}
