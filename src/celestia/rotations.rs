extern crate meval;
use self::meval::{Context, Expr};
use crate::log::error;
use crate::na::Matrix3;
use crate::time::{Epoch, J2000_OFFSET, MJD_OFFSET};
use crate::utils::{r1, r2, r3};
pub use celestia::xb::Identifier as XbId;
use std::cmp::PartialEq;
use std::collections::HashMap;
use std::fmt;

pub trait ParentRotation: fmt::Debug {
    fn dcm_to_parent(&self, datetime: Epoch) -> Option<Matrix3<f64>>;
}

#[derive(Debug)]
pub struct NoRotation;
impl ParentRotation for NoRotation {
    fn dcm_to_parent(&self, _: Epoch) -> Option<Matrix3<f64>> {
        Some(Matrix3::identity())
    }
}

/// Defines an Euler rotation, angle must be in radians
#[derive(Clone, Copy, Debug)]
pub enum EulerRotation {
    R1(f64),
    R2(f64),
    R3(f64),
}

impl EulerRotation {
    pub fn r1_from_degrees(angle_deg: f64) -> Self {
        Self::R1(angle_deg.to_radians())
    }
    pub fn r2_from_degrees(angle_deg: f64) -> Self {
        Self::R2(angle_deg.to_radians())
    }
    pub fn r3_from_degrees(angle_deg: f64) -> Self {
        Self::R3(angle_deg.to_radians())
    }
    /// Get the DCM from this Euler rotation
    pub fn dcm(&self) -> Matrix3<f64> {
        match *self {
            Self::R1(angle) => r1(angle),
            Self::R2(angle) => r2(angle),
            Self::R3(angle) => r3(angle),
        }
    }
}

impl ParentRotation for EulerRotation {
    fn dcm_to_parent(&self, _: Epoch) -> Option<Matrix3<f64>> {
        Some(self.dcm())
    }
}

/// A fixed three-axis Euler rotation
#[derive(Debug)]
pub struct Euler3Axis {
    /// The first rotation (e.g. R3)
    pub first: EulerRotation,
    /// The second rotation (e.g. R1)
    pub second: EulerRotation,
    /// The third and final rotation (e.g. R3, to complete a 3-1-1 rotation)
    pub third: EulerRotation,
}

impl ParentRotation for Euler3Axis {
    fn dcm_to_parent(&self, _: Epoch) -> Option<Matrix3<f64>> {
        Some(self.first.dcm() * self.second.dcm() * self.third.dcm())
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum AngleUnit {
    Degrees,
    Radians,
}

/// A time varying three-axis Euler rotation
pub struct Euler3AxisDt<'a> {
    pub base_context: Context<'a>,
    pub rot_order: [(EulerRotation, Expr); 3],
    pub unit: AngleUnit,
    pub is_ra_dec_w: bool,
}

impl<'a> Euler3AxisDt<'a> {
    /// Initialize a new time varying transformation.
    /// Reserved keywords in the context are "T" for centuries past 2000 Jan 1 12h TBD
    /// epoch (JDE 2451545.0), and "d" for days since that epoch.
    fn new(
        first_rot: (EulerRotation, Expr),
        second_rot: (EulerRotation, Expr),
        third_rot: (EulerRotation, Expr),
        context: &HashMap<String, f64>,
        unit: AngleUnit,
        is_ra_dec_w: bool,
    ) -> Self {
        let mut ctx = Context::default();
        for (var, value) in context {
            ctx.var(var, *value);
        }
        let rot_order = [first_rot, second_rot, third_rot];
        Self {
            base_context: ctx,
            rot_order,
            unit,
            is_ra_dec_w,
        }
    }

    /// Specify how to compute this frame from the provided Euler angles and their time varying expressions.
    /// Note that these angles define how to go from THIS frame TO the PARENT frame (e.g. Sun fixed to ICRF).
    pub fn from_euler_angles(
        first_rot: (EulerRotation, Expr),
        second_rot: (EulerRotation, Expr),
        third_rot: (EulerRotation, Expr),
        context: &HashMap<String, f64>,
        unit: AngleUnit,
    ) -> Self {
        Self::new(first_rot, second_rot, third_rot, context, unit, false)
    }

    /// A time varying Right ascension, Declination, and W frame
    /// Conversion TO parent frame (e.g. Sun body to ICRF) defined as:
    /// R3(-(alpha-90 deg)) * R1(delta - 90 deg) * R3(-W)
    /// Where alpha is the right ascension and delta the declination
    pub fn from_ra_dec_w(
        right_asc: Expr,
        declin: Expr,
        w: Expr,
        context: &HashMap<String, f64>,
        unit: AngleUnit,
    ) -> Self {
        Self::new(
            (EulerRotation::R3(0.0), right_asc),
            (EulerRotation::R1(0.0), declin),
            (EulerRotation::R3(0.0), w),
            context,
            unit,
            true,
        )
    }
}

impl<'a> fmt::Debug for Euler3AxisDt<'a> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "{:?}-{:?}-{:?}",
            self.rot_order[0].0, self.rot_order[1].0, self.rot_order[2].0
        )
    }
}

impl<'a> ParentRotation for Euler3AxisDt<'a> {
    fn dcm_to_parent(&self, datetime: Epoch) -> Option<Matrix3<f64>> {
        let days_d = datetime.as_jde_et_days() - MJD_OFFSET - J2000_OFFSET;
        let centuries_t = days_d / 36_525.0;
        // Now let's clone the context, and add the time variables.
        let mut ctx = self.base_context.clone();
        ctx.var("d", days_d);
        ctx.var("T", centuries_t);
        let mut dcm = Matrix3::identity();
        for (rot, expr) in &self.rot_order {
            // Compute the correct angle
            match expr.eval_with_context(&ctx) {
                Ok(eval_angle) => {
                    // Convert the angle to radians if needed
                    let angle = if self.unit == AngleUnit::Degrees {
                        eval_angle.to_radians()
                    } else {
                        eval_angle
                    };
                    let rot_with_angl = match rot {
                        EulerRotation::R1(_) => EulerRotation::R1(angle),
                        EulerRotation::R2(_) => EulerRotation::R2(angle),
                        EulerRotation::R3(_) => EulerRotation::R3(angle),
                    };
                    dcm *= rot_with_angl.dcm();
                }
                Err(e) => {
                    error!("{}", e);
                    // Stop here if something when wrong
                    return None;
                }
            }
        }
        Some(dcm)
    }
}
