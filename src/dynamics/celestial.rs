extern crate nalgebra as na;
use self::na::{Matrix6, MatrixMN, U1, U3, U36, U42, U6, Vector6, VectorN};
use super::Dynamics;
use celestia::{CelestialBody, CoordinateFrame, State};
use od::Linearization;
use std::f64;

/// `TwoBody` exposes the equations of motion for a simple two body propagation.
#[derive(Copy, Clone)]
pub struct TwoBody {
    pub mu: f64,
    time: f64,
    pos_vel: Vector6<f64>,
}

impl TwoBody {
    /// Initialize TwoBody dynamics given a provided gravitional parameter (as `mu`)
    pub fn from_state_vec_with_gm(state: &Vector6<f64>, mu: f64) -> TwoBody {
        TwoBody {
            time: 0.0,
            pos_vel: *state,
            mu: mu,
        }
    }

    /// Initialize TwoBody dynamics around a provided `CelestialBody` from the provided state vector (cf. nyx::celestia).
    pub fn from_state_vec<B: CelestialBody>(state: &Vector6<f64>) -> TwoBody {
        TwoBody {
            time: 0.0,
            pos_vel: *state,
            mu: B::gm(),
        }
    }
}

impl Dynamics for TwoBody {
    type StateSize = U6;
    fn time(&self) -> f64 {
        self.time
    }

    fn state(&self) -> VectorN<f64, Self::StateSize> {
        self.pos_vel
    }

    fn set_state(&mut self, new_t: f64, new_state: &VectorN<f64, Self::StateSize>) {
        self.time = new_t;
        self.pos_vel = *new_state;
    }

    fn eom(&self, _t: f64, state: &VectorN<f64, Self::StateSize>) -> VectorN<f64, Self::StateSize> {
        let radius = state.fixed_rows::<U3>(0).into_owned();
        let velocity = state.fixed_rows::<U3>(3).into_owned();
        let body_acceleration = (-self.mu / radius.norm().powi(3)) * radius;
        Vector6::from_iterator(velocity.iter().chain(body_acceleration.iter()).cloned())
    }
}

/// `TwoBody` exposes the equations of motion for a simple two body propagation. It inherently supports
/// the State Transition Matrix for orbital determination.
#[derive(Copy, Clone)]
pub struct TwoBodyWithStm {
    pub stm: Matrix6<f64>,
    pub two_body_dyn: TwoBody,
}

impl TwoBodyWithStm {
    /// Initialize TwoBody dynamics around a provided `CelestialBody` from the provided position and velocity state (cf. nyx::celestia).
    pub fn from_state<B: CelestialBody, F: CoordinateFrame>(state: State<F>) -> TwoBodyWithStm {
        TwoBodyWithStm {
            stm: Matrix6::identity(),
            two_body_dyn: TwoBody::from_state_vec::<B>(&state.to_cartesian_vec()),
        }
    }
}

impl Dynamics for TwoBodyWithStm {
    type StateSize = U42;
    fn time(&self) -> f64 {
        self.two_body_dyn.time()
    }

    fn state(&self) -> VectorN<f64, Self::StateSize> {
        let stm_as_vec = self.stm.fixed_resize::<U1, U36>(0.0);
        VectorN::<f64, Self::StateSize>::from_iterator(self.two_body_dyn.state().iter().chain(stm_as_vec.iter()).cloned())
    }

    fn set_state(&mut self, new_t: f64, new_state: &VectorN<f64, Self::StateSize>) {
        self.two_body_dyn
            .set_state(new_t, &new_state.fixed_rows::<U6>(0).into_owned());
        self.stm = new_state.fixed_rows::<U36>(6).fixed_resize::<U6, U6>(0.0);
    }

    fn eom(&self, t: f64, state: &VectorN<f64, Self::StateSize>) -> VectorN<f64, Self::StateSize> {
        let pos_vel = state.fixed_rows::<U6>(0).into_owned();
        let two_body_dt = self.two_body_dyn.eom(t, &pos_vel);
        let stm_dt = self.stm * self.gradient(t, &pos_vel);
        VectorN::<f64, Self::StateSize>::from_iterator(
            two_body_dt.iter().chain(stm_dt.fixed_resize::<U1, U36>(0.0).iter()).cloned(),
        )
    }
}

impl Linearization for TwoBodyWithStm {
    type StateSize = U6;

    fn gradient(&self, _t: f64, state: &VectorN<f64, Self::StateSize>) -> MatrixMN<f64, Self::StateSize, Self::StateSize> {
        let mut grad = MatrixMN::<f64, U6, U6>::zeros();
        // Top right is Identity 3x3
        grad[(3, 0)] = 1.0;
        grad[(4, 1)] = 1.0;
        grad[(2, 5)] = 1.0;
        // Bottom left is where the magic happens.
        let x = state[(0, 0)];
        let y = state[(1, 0)];
        let z = state[(2, 0)];
        let x2 = x.powi(2);
        let y2 = y.powi(2);
        let z2 = z.powi(2);
        let r232 = (x2 + y2 + z2).powf(3.0 / 2.0);
        let r252 = (x2 + y2 + z2).powf(5.0 / 2.0);

        // Add the body perturbations
        let dax_dx = 3.0 * self.two_body_dyn.mu * x2 / r252 - self.two_body_dyn.mu / r232;
        let dax_dy = 3.0 * self.two_body_dyn.mu * x * y / r252;
        let dax_dz = 3.0 * self.two_body_dyn.mu * x * z / r252;
        let day_dx = 3.0 * self.two_body_dyn.mu * x * y / r252;
        let day_dy = 3.0 * self.two_body_dyn.mu * y2 / r252 - self.two_body_dyn.mu / r232;
        let day_dz = 3.0 * self.two_body_dyn.mu * y * z / r252;
        let daz_dx = 3.0 * self.two_body_dyn.mu * x * z / r252;
        let daz_dy = 3.0 * self.two_body_dyn.mu * y * z / r252;
        let daz_dz = 3.0 * self.two_body_dyn.mu * z2 / r252 - self.two_body_dyn.mu / r232;

        // Let the gradient
        grad[(0, 3)] = dax_dx;
        grad[(0, 4)] = day_dx;
        grad[(0, 5)] = daz_dx;
        grad[(1, 3)] = dax_dy;
        grad[(1, 4)] = day_dy;
        grad[(1, 5)] = daz_dy;
        grad[(2, 3)] = dax_dz;
        grad[(2, 4)] = day_dz;
        grad[(2, 5)] = daz_dz;

        grad
    }
}
