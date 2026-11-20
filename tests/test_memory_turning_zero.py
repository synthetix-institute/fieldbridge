"""A stable zero that turns much faster than it decays, near the end of an implicit run of a relaxation that has not
settled (correction 35 of the discovery engine).

Correction 34 settled such a run at the zero that Newton reaches from its end, when the zero has no growing direction,
attracts, is not a center, and the run ended closer to it than it began. Two things went wrong on fast-turning foci:

- BDF damps a fast turn numerically. The Bautin point r' = r (mu + r^2 - r^4), mu = -0.2, turning at 2e4, has a stable
  origin, an unstable cycle at r = 0.526 and a stable cycle at r = 0.851. Runs that start beyond the unstable cycle move
  outwards, slowly, but BDF carried them inwards, and they were settled at the origin: the card read every start settled
  and hid the stable cycle. On a zero that turns faster than it decays the endpoints of a run cannot tell an approach
  from that damping. The run is now settled only by the speed criterion, and the zero is still a stored state.
- Newton's step was capped at half the base length of the carrier. A pendulum that turns at w = 2e4 on a carrier of scale
  1 ends its runs at velocities up to 1e4, and 80 capped steps reach no zero more than 40 away. Where the capped search
  does not converge, Newton runs without the cap, its steps halved until the drift decreases.

A zero with real rates, as that of a stiff relaxation, is settled as under correction 34."""
from __future__ import annotations

import numpy as np

from fieldbridge.memory import analysis as an
from fieldbridge.memory import spec


def _spec(name, variables, parameters, drift, control, control_range, scale=1.0):
    return {"schema": "fieldbridge-memory/1", "kind": "equations", "name": name, "question": "what does it store?",
            "assumptions": ["test model"], "provenance": {"source": "test"},
            "carrier": {"kind": "euclid", "name": "state", "variables": variables, "scale": scale},
            "parameters": parameters, "drift": drift, "control": {"name": control, "range": control_range},
            "noise": 0.01, "dt": 0.01}


# the fast Bautin point of set 30: stable origin, unstable cycle at r = 0.526, stable cycle at r = 0.851
BAUTIN = _spec("fast Bautin point", ["x", "y"], {"mu": -0.2, "w": 2e4},
               {"x": "mu*x - w*y + x*((x**2 + y**2) - (x**2 + y**2)**2)",
                "y": "w*x + mu*y + y*((x**2 + y**2) - (x**2 + y**2)**2)"}, "mu", [-0.24, -0.15], scale=0.5)
# the fast pendulum of set 34
PENDULUM = _spec("fast pendulum", ["x", "y"], {"w": 2e4, "g": 0.4}, {"x": "y", "y": "-w**2*sin(x) - g*y"}, "g",
                 [0.2, 0.8])
# the stiff pair of set 17: rates 1e5 and 1e-3, whose runs end unsettled at the end of their time
STIFF = _spec("stiff pair", ["x", "y"], {"k": 1e5, "e": 1e-3}, {"x": "-k*x", "y": "-e*y"}, "k", [5e4, 2e5])


def test_a_run_beyond_the_unstable_cycle_is_not_settled_at_the_origin():
    real = spec.load(BAUTIN)
    assert 300.0 / an.safe_dt(real) > an.RELAX_BUDGET  # its relaxation is implicit
    q, speed = an.relax(real, np.array([[0.6, 0.0]]))
    assert not np.linalg.norm(q[0]) < 1e-6  # it has not been put at the origin (its run ended at the budget: NaN)
    zeros = []
    an.relax(real, np.array([[0.6, 0.0]]), zeros=zeros)
    assert len(zeros) == 1 and np.linalg.norm(zeros[0][1]) < 1e-9  # the origin is the stable zero near it


def test_the_fast_bautin_point_stores_the_origin_and_its_starts_do_not_settle():
    states, counts, spectra, unconverged, continuum = an.stored_states(spec.load(BAUTIN), np.random.default_rng(0), 8)
    assert len(states) == 1 and np.linalg.norm(states[0]) < 1e-9
    assert unconverged == 8


def test_the_fast_pendulum_stores_its_state_from_every_search():
    for seed in (0, 1):
        states, counts, spectra, unconverged, continuum = an.stored_states(spec.load(PENDULUM),
                                                                           np.random.default_rng(seed), 8)
        assert len(states) == 1 and np.linalg.norm(states[0]) < 1e-9


def test_newton_without_the_cap_reaches_a_distant_zero():
    real = spec.load(PENDULUM)
    end = np.array([-0.0156, -351.2382])  # the end of a run of the probe of correction 34
    z, ok = an.newton(real, end)
    assert not ok  # 80 steps of 0.5 reach no zero 351 away
    z, ok = an.newton(real, end, descent=True)
    assert ok and np.linalg.norm(z) < 1e-9


def test_a_zero_with_real_rates_is_settled_as_before():
    real = spec.load(STIFF)
    assert 300.0 / an.safe_dt(real) > an.RELAX_BUDGET
    start = np.array([[0.8, 0.9], [-0.5, -0.7]])
    assert np.all(np.abs(1e-3 * start[:, 1] * np.exp(-0.3)) > 1e-4)  # unsettled by the speed criterion at t = 300
    q, speed = an.relax(real, start)
    assert np.all(np.linalg.norm(q, axis=1) < 1e-9)  # settled at the origin by the rule, as under correction 34
