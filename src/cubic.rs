use arrayvec::ArrayVec;

use crate::{Cubic, different_signs};

impl Cubic {
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
        let a = 3.0 * self.coeffs[3];
        let b_2 = self.coeffs[2];
        let c = self.coeffs[1];
        let disc_4 = b_2 * b_2 - a * c;

        if !disc_4.is_finite() {
            return self.rescaled_critical_points();
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

    #[cold]
    fn rescaled_critical_points(&self) -> Option<(f64, f64)> {
        let scale = 2.0f64.powi(-515);
        (*self * scale).critical_points()
    }

    fn one_root(
        &self,
        lower: f64,
        upper: f64,
        lower_val: f64,
        upper_val: f64,
        x_error: f64,
    ) -> f64 {
        let deriv = self.deriv();
        if !deriv.is_finite() {
            return f64::NAN;
        }
        crate::yuksel::find_root(
            |x| self.eval(x),
            |x| deriv.eval(x),
            lower,
            upper,
            lower_val,
            upper_val,
            x_error,
        )
    }

    // This has to be pub because we're benchmarking it right now.
    #[doc(hidden)]
    pub fn root_between(self, lower: f64, upper: f64, x_error: f64) -> f64 {
        self.one_root(lower, upper, self.eval(lower), self.eval(upper), x_error)
    }

    fn first_root(self, lower: f64, upper: f64, x_error: f64) -> Option<f64> {
        if let Some((x0, x1)) = self.critical_points() {
            let possible_endpoints: [f64; 3] = [x0, x1, upper];
            let mut last = lower;
            let mut last_val = self.eval(last);
            for x in possible_endpoints {
                if x > last && x <= upper {
                    let val = self.eval(x);
                    if different_signs(last_val, val) {
                        return Some(self.one_root(last, x, last_val, val, x_error));
                    }

                    last = x;
                    last_val = val;
                }
            }
            None
        } else {
            let lower_val = self.eval(lower);
            let upper_val = self.eval(upper);
            if different_signs(lower_val, upper_val) {
                Some(self.one_root(lower, upper, lower_val, upper_val, x_error))
            } else {
                None
            }
        }
    }

    /// Computes all roots between `lower` and `upper`, to the desired accuracy.
    ///
    /// We make no guarantees about multiplicity. In fact, if there's a
    /// double-root that isn't a triple-root (and therefore has no sign change
    /// nearby) then there's a good chance we miss it altogether. This is
    /// fine if you're using this root-finding to optimize a quartic, because
    /// double-roots of the derivative aren't local extrema.
    pub fn roots_between(self, lower: f64, upper: f64, x_error: f64) -> ArrayVec<f64, 3> {
        let mut ret = ArrayVec::new();
        let mut scratch = ArrayVec::new();
        self.roots_between_with_buffer(lower, upper, x_error, &mut scratch, &mut ret);
        ret
    }

    pub(crate) fn roots_between_with_buffer<const M: usize>(
        self,
        lower: f64,
        upper: f64,
        x_error: f64,
        _scratch: &mut ArrayVec<f64, M>,
        out: &mut ArrayVec<f64, M>,
    ) {
        if let Some(r) = self.first_root(lower, upper, x_error) {
            out.push(r);
            let quad = self.deflate(r);
            if let Some((x0, x1)) = quad.positive_discriminant_roots() {
                if lower <= x0 && x0 <= upper {
                    out.push(x0);
                }
                if lower <= x1 && x1 <= upper {
                    out.push(x1);
                }
            }
        }
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
    // This has to be pub because we're benchmarking it right now.
    //
    // We don't really want this public because it appears to be slower than
    // the other method. This one does a Newton search for each bracketing
    // interval, while the other one just does a single Newton search
    // and then deflates to find the other two roots.
    #[doc(hidden)]
    pub fn roots_between_multiple_searches(
        self,
        lower: f64,
        upper: f64,
        x_error: f64,
    ) -> ArrayVec<f64, 3> {
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
                    ret.push(self.one_root(last, x, last_val, val, x_error));
                }

                last = x;
                last_val = val;
            }
        }
        ret
    }

    #[doc(hidden)]
    pub fn precondition(&self) -> Cubic {
        // Truncate coefficients too close to zero, to ensure that there's
        // no underflow when calculating the discriminant.
        let min_coeff = 2.0f64.powi(-256);
        let truncate = |x: &mut f64| {
            if x.abs() <= min_coeff {
                *x = 0.0
            }
        };

        // We can't just truncate, because if some other coefficient is just above
        // min_coeff then it will introduce a big relative error. So we renormalize
        // if things are too close.
        let small_coeff = 2.0f64.powi(-64);

        let large_coeff = 2.0f64.powi(64);

        let mut c = *self;
        if (self.magnitude() != 0.0 && self.magnitude() <= small_coeff)
            || self.magnitude() >= large_coeff
        {
            c /= self.magnitude();
        }

        truncate(&mut c.coeffs[0]);
        truncate(&mut c.coeffs[1]);
        truncate(&mut c.coeffs[2]);
        truncate(&mut c.coeffs[3]);
        c
    }

    pub fn roots_blinn(&self) -> ArrayVec<f64, 3> {
        let mut ret = ArrayVec::new();
        let a = self.coeffs[3];
        let b = self.coeffs[2] * (1.0 / 3.0);
        let c = self.coeffs[1] * (1.0 / 3.0);
        let d = self.coeffs[0];

        let delta_1 = a * c - b * b;
        let delta_2 = a * d - b * c;
        let delta_3 = b * d - c * c;
        let disc = 4.0 * delta_1 * delta_3 - delta_2 * delta_2;

        if !disc.is_finite() {
            return ret;
        }
        // What about disc = 0?
        // For now, we put it in the one-root case, although in principle it could also
        // be a single root and a double root. The issue with the other branch is
        // that we might end up with `atan(0, 0)`, which gives NaN. Blinn says the NaN
        // doesn't matter because you end up multiplying it by \bar C = 0, but (1)
        // floats don't work that way without a little effort, and (2) it's possible to
        // have disc = \bar D = 0.0 (numerically) and \bar C \ne 0.
        let mut push = |x: f64| {
            if x.is_finite() {
                ret.push(x);
            }
        };
        if disc <= 0.0 {
            //dbg!(disc);
            let (tilde_a, tilde_c, tilde_d) = if b * b * b * d >= a * c * c * c {
                (a, delta_1, -2.0 * b * delta_1 + a * delta_2)
            } else {
                (d, delta_3, -d * delta_2 + 2.0 * c * delta_3)
            };
            //dbg!(tilde_a, tilde_c, tilde_d);
            let t_0 = -tilde_a.copysign(tilde_d) * (-disc).sqrt();
            let t_1 = -tilde_d + t_0;
            let p = (t_1 / 2.0).cbrt();

            let q = if t_0 == t_1 { -p } else { -tilde_c / p };
            //dbg!(p, q);
            let tilde_x = if tilde_c <= 0.0 {
                p + q
            } else {
                -tilde_d / (p * p + q * q + tilde_c)
            };

            let (x, w) = if b * b * b * d >= a * c * c * c {
                (tilde_x - b, a)
            } else {
                (-d, tilde_x + c)
            };

            push(x / w);
        } else {
            //dbg!(disc);
            fn one_root(a_or_d: f64, disc: f64, bar_c: f64, bar_d: f64) -> (f64, f64) {
                let sqrt_c = (-bar_c).sqrt();
                let theta = (1.0 / 3.0) * (a_or_d * disc.sqrt()).atan2(-bar_d).abs();
                let (sin_theta, cos_theta) = theta.sin_cos();
                //dbg!(theta, cos_theta);
                let tilde_x_1 = 2.0 * sqrt_c * cos_theta;
                let tilde_x_3 = sqrt_c * (-cos_theta - 3.0f64.sqrt() * sin_theta);
                (tilde_x_1, tilde_x_3)
            }

            let bar_c_a = delta_1;
            let bar_d_a = -2.0 * b * delta_1 + a * delta_2;
            let (tilde_x_1_a, tilde_x_3_a) = one_root(a, disc, bar_c_a, bar_d_a);

            let bar_c_d = delta_3;
            let bar_d_d = -d * delta_2 + 2.0 * c * delta_3;
            let (tilde_x_1_d, tilde_x_3_d) = one_root(d, disc, bar_c_d, bar_d_d);

            let tilde_x_l = if tilde_x_1_a + tilde_x_3_a > 2.0 * b {
                tilde_x_1_a
            } else {
                tilde_x_3_a
            };
            let tilde_x_s = if tilde_x_1_d + tilde_x_3_d < 2.0 * c {
                tilde_x_1_d
            } else {
                tilde_x_3_d
            };

            let (x_l, w_l) = (tilde_x_l - b, a);
            let (x_s, w_s) = (-d, tilde_x_s + c);

            let e = w_l * w_s;
            let f = -x_l * w_s - w_l * x_s;
            let g = x_l * x_s;

            let (x_m, w_m) = (c * f - b * g, c * e - b * f);

            push(x_l / w_l);
            push(x_s / w_s);
            push(x_m / w_m);
        }
        ret
    }

    pub fn roots_blinn_and_deflate(&self) -> ArrayVec<f64, 3> {
        let mut ret = ArrayVec::new();
        let a = self.coeffs[3];
        let b = self.coeffs[2] * (1.0 / 3.0);
        let c = self.coeffs[1] * (1.0 / 3.0);
        let d = self.coeffs[0];

        let delta_1 = a * c - b * b;
        let delta_2 = a * d - b * c;
        let delta_3 = b * d - c * c;
        let disc = 4.0 * delta_1 * delta_3 - delta_2 * delta_2;

        if !disc.is_finite() {
            return ret;
        }
        //dbg!(delta_1, delta_2, delta_3, disc);

        // TODO: what about disc = 0?
        if disc <= 0.0 {
            let (tilde_a, tilde_c, tilde_d) = if b * b * b * d >= a * c * c * c {
                (a, delta_1, -2.0 * b * delta_1 + a * delta_2)
            } else {
                (d, delta_3, -d * delta_2 + 2.0 * c * delta_3)
            };
            let t_0 = -tilde_a.copysign(tilde_d) * (-disc).sqrt();
            let t_1 = -tilde_d + t_0;
            let p = (t_1 / 2.0).cbrt();

            let q = if t_0 == t_1 { -p } else { -tilde_c / p };
            let tilde_x = if tilde_c <= 0.0 {
                p + q
            } else {
                -tilde_d / (p * p + q * q + tilde_c)
            };

            let (x, w) = if b * b * b * d >= a * c * c * c {
                (tilde_x - b, a)
            } else {
                (-d, tilde_x + c)
            };

            if x.is_finite() && w.is_finite() {
                ret.push(x / w);
            }
        } else {
            fn one_root(a_or_d: f64, disc: f64, bar_c: f64, bar_d: f64) -> (f64, f64) {
                let sqrt_c = (-bar_c).sqrt();
                let theta = (1.0 / 3.0) * (a_or_d * disc.sqrt()).atan2(-bar_d).abs();
                let (sin_theta, cos_theta) = theta.sin_cos();
                let tilde_x_1 = 2.0 * sqrt_c * cos_theta;
                let tilde_x_3 = sqrt_c * (-theta.cos() - 3.0f64.sqrt() * sin_theta);
                (tilde_x_1, tilde_x_3)
            }

            // FIXME: I'm confused about which choice is supposed to give me the
            // small-magnitude root...
            let bar_c_a = delta_1;
            let bar_d_a = -2.0 * b * delta_1 + a * delta_2;
            let (tilde_x_1_a, tilde_x_3_a) = one_root(a, disc, bar_c_a, bar_d_a);

            let bar_c_d = delta_3;
            let bar_d_d = -d * delta_2 + 2.0 * c * delta_3;
            let (tilde_x_1_d, tilde_x_3_d) = one_root(d, disc, bar_c_d, bar_d_d);

            let tilde_x_l = if tilde_x_1_a + tilde_x_3_a > 2.0 * b {
                tilde_x_1_a
            } else {
                tilde_x_3_a
            };
            //dbg!(tilde_x_1_a - b, tilde_x_3_a - b, tilde_x_l - b);
            let tilde_x_s = if tilde_x_1_d + tilde_x_3_d < 2.0 * c {
                tilde_x_1_d
            } else {
                tilde_x_3_d
            };

            let (x_l, w_l) = (tilde_x_l - b, a);
            let (x_s, w_s) = (-d, tilde_x_s + c);

            let x = if (x_l * w_s).abs() <= (x_s * w_l).abs() {
                x_l / w_l
            } else {
                x_s / w_s
            };
            if x.is_finite() {
                ret.push(x);
            }
            let q = self.deflate(x);
            //dbg!(&q);
            if q.is_finite() {
                ret.extend(q.roots());
                ret.sort_by(|x, y| x.partial_cmp(y).unwrap());
            }
        }
        ret
    }
}

#[cfg(test)]
mod tests {
    use crate::Cubic;

    const TRICKY_CUBICS: [Cubic; 5] = [
        // This one has infinite discriminant.
        Cubic::new([
            1.6149620090145706e-94,
            1.6149620090145634e-94,
            1.6149620090145663e-94,
            9.66803867245343e272,
        ]),
        // This one has a very large second root (-7e202), which causes roots_blinn to NaN on the last one.
        //
        // When deflating, it gives a quadratic that's basically zero and so it underflows the discriminant
        // and ends up reporting just a single root.
        Cubic::new([
            -6.323283382275869e98,
            3.0957754283429482e-307,
            3.095775428342964e-307,
            3.095775428342951e-307,
        ]),
        // Here's one with sane coefficients, but a similar issue as
        // the last one.
        Cubic::new([
            -8.522348907129e-161,
            4.471145208374078e-67,
            -0.052026185927646074,
            -2.9441090045938734e-57,
        ]),
        // Here's one where the discriminant is numerically zero, causing some stability issues
        // for Blinn.
        Cubic::new([
            -2.5162489269306657e-175,
            -2.516248926930655e-175,
            -2.5162489269306522e-175,
            -0.39205037382350466,
        ]),
        Cubic::new([
            -6.428720163649757e103,
            -6.428720163649766e103,
            -3.3646756114322413e-74,
            -3.3646756114322547e-74,
        ]),
    ];

    #[test]
    fn smoke() {
        // Here's an example where Blinn's method has a large error. The
        // small-magnitude root is about -1.0 and the large magnitude root is
        // apparently of order 1e225, which then causes big errors when trying
        // to compute the third root. I guess in general we have expect that
        // the magnitude of the big root affects the error in the middle root.
        //
        // Correction: the middle root is actually ok here. It has a pretty
        // large magnitude (1e17ish), and so it's allowed to not evaluate
        // super close to zero.
        // let poly = super::Cubic {
        //     c0: -3.565233507454652e74,
        //     c1: -3.5652335074546437e74,
        //     c2: -1.2298855640101194e-17,
        //     c3: 9.133009604987547e-243,
        // };

        // let poly = Cubic {
        //     c0: 5.5174454041519107,
        //     c1: -1.6144740273415798e-245,
        //     c2: -3.892738574215212e-288,
        //     c3: 3.0860491510941517e-292,
        // };
        let poly = Cubic::new([
            -6.428720163649757e103,
            -6.428720163649766e103,
            -3.3646756114322413e-74,
            -3.3646756114322547e-74,
        ]);

        let roots = poly.precondition().roots_blinn();
        //let roots = poly.roots_between_multiple_searches(-10.0, 10.0, 1e-12);
        dbg!(&roots);
        for r in roots {
            dbg!(poly.eval(r));
        }
    }

    #[test]
    fn bad_for_blinn() {
        for c in TRICKY_CUBICS {
            dbg!(c.roots_blinn());
            dbg!(c.roots_blinn_and_deflate());
        }
    }

    // Asserts that the supplied "roots" are close to being roots of the
    // cubic, in the sense that the cubic evaluates to approximately zero
    // on each of the roots.
    fn check_root_values(c: &Cubic, roots: &[f64]) {
        // Arbitrary cubics can have coefficients with wild magnitudes,
        // so we need to adjust our error expectations accordingly.
        let magnitude = c.magnitude().max(1.0);
        let accuracy = magnitude * 1e-12;

        for r in roots {
            // We can't expect great accuracy for very large roots,
            // because the polynomial evaluation will involve very
            // large terms.
            let accuracy = accuracy * r.abs().powi(3).max(1.0);
            let y = c.eval(*r);
            if y.is_finite() {
                assert!(
                    y.abs() <= accuracy,
                    "cubic {c:?} had root {r} evaluate to {y:?}, but expected {accuracy:?}"
                );
            }
        }
    }

    #[test]
    fn root_evaluation() {
        arbtest::arbtest(|u| {
            let c = crate::arbitrary::cubic(u)?;

            // We could have a wider range of roots, but then we might need
            // to lower the accuracy depending on what the actual root is: the
            // intermediate computations scale like the cube of the root.
            let roots = c.roots_between_multiple_searches(-10.0, 10.0, 1e-13);
            check_root_values(&c, &roots);
            let roots = c.roots_between(-10.0, 10.0, 1e-13);
            check_root_values(&c, &roots);

            let preconditioned = c.precondition();
            check_root_values(&c, &preconditioned.roots_blinn());
            check_root_values(&c, &preconditioned.roots_blinn_and_deflate());

            // Even preconditioning is not enough for kurbo's current solver.
            // let Cubic { c0, c1, c2, c3 } = preconditioned;
            // dbg!(preconditioned);
            // check_root_values(&c, &kurbo::common::solve_cubic(c0, c1, c2, c3));

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
            let magnitude = c.magnitude().max(1.0);
            let accuracy = magnitude * 1e-12;

            // We could have a wider range of roots, but then we might need
            // to lower the accuracy depending on what the actual root is: the
            // intermediate computations scale like the cube of the root.
            let &[c0, c1, c2, c3] = c.coeffs();
            for r in kurbo::common::solve_cubic(c0, c1, c2, c3) {
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
