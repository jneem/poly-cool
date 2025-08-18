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

    pub fn coeffs(&self) -> &[f64; N] {
        &self.coeffs
    }

    pub fn eval(&self, x: f64) -> f64 {
        let mut acc = 0.0;
        for c in self.coeffs.iter().rev() {
            // It would be nice to use `f64::mul_add` here, but it's slow on
            // architectures that don't have a dedicated instruction.
            acc = acc * x + c;
        }
        acc
    }

    pub fn magnitude(&self) -> f64 {
        let mut max = 0.0f64;
        for c in &self.coeffs {
            max = max.max(c.abs());
        }
        max
    }

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

impl_deriv_and_deflate!(3, 2);
impl_deriv_and_deflate!(4, 3);
impl_deriv_and_deflate!(5, 4);
impl_deriv_and_deflate!(6, 5);

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
