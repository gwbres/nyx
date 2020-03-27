extern crate hifitime;
extern crate nalgebra as na;

extern crate nyx_space as nyx;

use hifitime::{Epoch, SECONDS_PER_DAY};
use nyx::celestia::eclipse::{EclipseLocator, EclipseState};
use nyx::celestia::{bodies, Cosm, LTCorr, OrbitState};
use nyx::dynamics::celestial::CelestialDynamics;
use nyx::propagators::{PropOpts, Propagator};
use std::sync::mpsc;
use std::sync::mpsc::{Receiver, Sender};
use std::thread;

#[test]
fn leo_sun_earth_eclipses() {
    let prop_time = 2.0 * SECONDS_PER_DAY;

    let cosm = Cosm::from_xb("./de438s");
    let earth = cosm.frame_by_id(bodies::EARTH);
    let sun = cosm.frame_by_id(bodies::SUN);

    let start_time = Epoch::from_gregorian_tai_at_midnight(2020, 1, 1);

    let leo = OrbitState::keplerian(6778.0, 0.1, 60.0, 0.0, 0.0, 0.0, start_time, earth);

    let (truth_tx, truth_rx): (Sender<OrbitState>, Receiver<OrbitState>) = mpsc::channel();

    let bodies = vec![bodies::SUN, bodies::JUPITER_BARYCENTER];

    thread::spawn(move || {
        let cosm = Cosm::from_xb("./de438s");
        let mut dynamics = CelestialDynamics::new(leo, bodies, &cosm);
        let mut prop = Propagator::default(&mut dynamics, &PropOpts::with_fixed_step(60.0));
        prop.tx_chan = Some(&truth_tx);
        prop.until_time_elapsed(prop_time);
    });

    // Initialize the EclipseLocator
    let e_loc = EclipseLocator {
        light_source_id: bodies::SUN,
        shadow_body_ids: vec![bodies::EARTH],
        cosm: &cosm,
        correction: LTCorr::None,
    };

    // Receive the states on the main thread.
    let mut prev_eclipse_state = EclipseState::Umbra;
    let mut cnt_changes = 0;
    while let Ok(rx_state) = truth_rx.recv() {
        let new_eclipse_state = e_loc.compute(&rx_state);
        if new_eclipse_state != prev_eclipse_state {
            println!(
                "{:.6} now in {:?}",
                rx_state.dt.as_jde_tai_days(),
                new_eclipse_state
            );
            prev_eclipse_state = new_eclipse_state;
            cnt_changes += 1;
        }
    }

    assert_eq!(cnt_changes, 79, "wrong number of eclipse state changes");
}

#[test]
fn geo_sun_earth_eclipses() {
    let prop_time = 2.0 * SECONDS_PER_DAY;

    let cosm = Cosm::from_xb("./de438s");
    let earth = cosm.frame_by_id(bodies::EARTH);
    let sun = cosm.frame_by_id(bodies::SUN);

    // GEO are in shadow or near shadow during the equinoxes.
    let start_time = Epoch::from_gregorian_tai_at_midnight(2020, 3, 19);

    let leo = OrbitState::keplerian(42000.0, 0.1, 0.1, 0.0, 0.0, 0.0, start_time, earth);

    let (truth_tx, truth_rx): (Sender<OrbitState>, Receiver<OrbitState>) = mpsc::channel();

    let bodies = vec![bodies::SUN, bodies::JUPITER_BARYCENTER];

    thread::spawn(move || {
        let cosm = Cosm::from_xb("./de438s");
        let mut dynamics = CelestialDynamics::new(leo, bodies, &cosm);
        let mut prop = Propagator::default(&mut dynamics, &PropOpts::with_fixed_step(60.0));
        prop.tx_chan = Some(&truth_tx);
        prop.until_time_elapsed(prop_time);
    });

    // Initialize the EclipseLocator
    let e_loc = EclipseLocator {
        light_source_id: bodies::SUN,
        shadow_body_ids: vec![bodies::EARTH],
        cosm: &cosm,
        correction: LTCorr::None,
    };

    // Receive the states on the main thread.
    let mut prev_eclipse_state = EclipseState::Umbra;
    let mut cnt_changes = 0;
    while let Ok(rx_state) = truth_rx.recv() {
        let new_eclipse_state = e_loc.compute(&rx_state);
        if new_eclipse_state != prev_eclipse_state {
            println!(
                "{:.6} now in {:?}",
                rx_state.dt.as_jde_tai_days(),
                new_eclipse_state
            );
            prev_eclipse_state = new_eclipse_state;
            cnt_changes += 1;
        }
    }

    assert_eq!(cnt_changes, 62, "wrong number of eclipse state changes");
}
