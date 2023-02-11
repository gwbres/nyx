/*
    Nyx, blazing fast astrodynamics
    Copyright (C) 2022 Christopher Rabotin <christopher.rabotin@gmail.com>

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Affero General Public License as published
    by the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Affero General Public License for more details.

    You should have received a copy of the GNU Affero General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

use super::{InterpState, TrajError};
use crate::linalg::allocator::Allocator;
use crate::linalg::DefaultAllocator;
use crate::polyfit::Polynomial;
use crate::time::{Duration, Epoch};
use crate::utils::normalize;
use crate::NyxError;

pub(crate) const SPLINE_DEGREE: usize = 8;
pub(crate) const INTERPOLATION_SAMPLES: usize = 8;

/// Stores a segment of an interpolation, a spline. Each spline is a polynomial
#[derive(Clone)]
pub struct Spline<S: InterpState>
where
    DefaultAllocator:
        Allocator<f64, S::VecLength> + Allocator<f64, S::Size> + Allocator<f64, S::Size, S::Size>,
{
    pub(crate) start_epoch: Epoch,
    pub(crate) duration: Duration,
    // TODO: When rustc is cool with more const generics, switch this to a [Poly<{S::DEGREE}>; S::CURVES]
    pub(crate) polynomials: Vec<Polynomial<SPLINE_DEGREE>>,
    pub(crate) end_state: S,
}

impl<S: InterpState> Spline<S>
where
    DefaultAllocator:
        Allocator<f64, S::VecLength> + Allocator<f64, S::Size> + Allocator<f64, S::Size, S::Size>,
{
    /// Evaluate a specific segment at the provided Epoch, requires an initial state as a "template"
    pub fn evaluate(&self, from: S, epoch: Epoch) -> Result<S, NyxError> {
        // Compute the normalized time
        let dur_into_window = epoch - self.start_epoch;
        if dur_into_window > self.duration {
            return Err(NyxError::Trajectory(TrajError::OutOfSpline {
                req_epoch: epoch,
                req_dur: dur_into_window,
                spline_dur: self.duration,
            }));
        }

        let t_prime = normalize(
            dur_into_window.to_seconds(),
            0.0,
            self.duration.to_seconds(),
        );

        let mut state = from;

        // Rebuild the polynominals
        for (i, param) in S::params().iter().enumerate() {
            let (value, value_dt) = self.polynomials[i].eval_n_deriv(t_prime);
            state.set_value_and_deriv(param, value, value_dt)?;
        }
        state.set_epoch(epoch);

        Ok(state)
    }
}
