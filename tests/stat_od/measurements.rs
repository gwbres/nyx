extern crate hifitime;
extern crate nalgebra as na;
extern crate nyx_space as nyx;

#[test]
fn nil_measurement() {
    use self::hifitime::julian::*;
    use self::nyx::celestia::{Cosm, Geoid, State};
    use self::nyx::od::ranging::GroundStation;
    use self::nyx::od::Measurement;
    use std::f64::EPSILON;
    // Let's create a station and make it estimate the range and range rate of something which is strictly in the same spot.

    let lat = -7.906_635_7;
    let long = 345.5975;
    let height = 0.0;
    let dt = ModifiedJulian::j2000();
    let cosm = Cosm::from_xb("./de438s");
    let earth_geoid = cosm.geoid_from_id(3).unwrap();

    let station = GroundStation::from_noise_values("nil", 0.0, lat, long, height, 0.0, 0.0);

    let at_station = State::<Geoid>::from_geodesic(lat, long, height, dt, earth_geoid);

    let meas = station.measure(at_station, dt.into_instant());

    let h_tilde = *meas.sensitivity();
    println!("{}", h_tilde);
    assert!(h_tilde[(0, 0)].is_nan(), "expected NaN");
    assert!(h_tilde[(0, 1)].is_nan(), "expected NaN");
    assert!(h_tilde[(0, 2)].is_nan(), "expected NaN");
    assert!(h_tilde[(1, 0)].is_nan(), "expected NaN");
    assert!(h_tilde[(1, 1)].is_nan(), "expected NaN");
    assert!(h_tilde[(1, 2)].is_nan(), "expected NaN");
    assert!(h_tilde[(1, 3)].is_nan(), "expected NaN");
    assert!(h_tilde[(1, 4)].is_nan(), "expected NaN");
    assert!(h_tilde[(1, 5)].is_nan(), "expected NaN");

    assert!(meas.observation()[(0, 0)] - 0.0 < EPSILON, "observation is not range=0");
    assert!(meas.observation()[(1, 0)].is_nan(), "observation is not range=0");
}
