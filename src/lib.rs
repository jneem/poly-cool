mod cubic;
mod poly;
mod quadratic;

#[cfg(any(test, feature = "arbitrary"))]
pub mod arbitrary;

pub use cubic::Cubic;
pub use poly::Poly;
pub use quadratic::Quadratic;

trait TerminationCondition {
    fn stop(&self, last_step: f64, value: f64) -> bool;
}

struct InputError(f64);

impl TerminationCondition for InputError {
    fn stop(&self, last_step: f64, _value: f64) -> bool {
        last_step.abs() <= self.0
    }
}

struct ValueError(f64);

impl TerminationCondition for ValueError {
    fn stop(&self, _last_step: f64, value: f64) -> bool {
        value.abs() <= self.0
    }
}

fn different_signs(x: f64, y: f64) -> bool {
    (x < 0.0) != (y < 0.0)
}
