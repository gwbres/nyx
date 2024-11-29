//! Solid body tidal effects

/// Estimates Earth crust displacement at `rx` ECEF coordinates,
/// due to solid tide effects induced by Sun and Moon.
/// ## Input
/// - r_sun: Unit vector from geocenter to Sun (at epoch of estimation)
/// - r_moon: Unit vector from geocenter to Sun (at epoch of estimation)
/// - r_rx: Unit vector from geocenter to calculation coordinates (at epoch of estimation)
/// ## Output
/// - Vector3 displacement vector in mm ECEF
pub fn earth_tidal_crust_displacement_ecef_mm(
    r_rx: Vector3<f64>,
    r_sun: Vector3<f64>,
    r_moon: Vector3<f64>,
) -> Vector3<f64> {
    const H2_LOVE_NUMBER: f64 = 0.6078;
    const L2_SHIDA_NUMBER: f64 = 0.0847;
    const G_M_EARTH: f64 = 3.986E14;
    const G_M_SUN: f64 = 1.327E20;
    const G_M_MOON: f64 = 4.904E12;

    let r_sun_mag = r_sun.mag();
    let r_moon_mag = r_moon.mag();

    let mut sum = Vector3::<f64>::default();

    // to make this generic: 3=2+moons.len()
    // will not work for binary star systems
    for j in 2..3 {
        let is_star = j == 3;

        // Constant numerator
        let num_g = if is_star {
            G_M_SUN * EARTH_RADIUS_M.powi(4)
        } else {
            G_M_MOON * EARTH_RADIUS_M.powi(4)
        };

        // Constant denominator
        let denom_g = if is_star {
            G_M_EARTH * r_sun_mag.powi(3)
        } else {
            G_M_EARTH * r_moon_mag.powi(3)
        };

        // Calculate the Body / Rx dot product
        let body_rx_dot = if is_star {
            r_sun.dot(r_rx)
        } else {
            r_sun.dot(r_rx)
        };

        // Calculate the h2 love term
        let h2_love = H2_LOVE_NUMBER * (3.0 / 2.0 * body_rx_dot.powi(2) - 0.5);

        // Calculate the right hand side
        let rhs = if is_star {
            r_sun - (body_rx_dot) * r_rx
        } else {
            r_moon - (body_rx_dot) * r_rx
        };

        // Calculate the l2 shida term
        let l2_shida = 3.0 * L2_SHIDA_NUMBER * body_rx_dot * rhs;

        // Calculate total term
        sum += num_g / denom_g * (h2_love + l2_shida);
    }
    sum
}

#[cfg(test)]
mod test {
    use super::earth_tidal_crust_displacement_ecef_mm;

    #[test]
    fn earth_tidal_crust_displacement() {
        let r_rx = Vector3::<f64>::default();
        let r_sun = Vector3::<f64>::default();
        let r_moon = Vector3::<f64>::default();
        let r_disp = earth_tidal_crust_displacement_ecef_mm();
    }
}
