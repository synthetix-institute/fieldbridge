"""The step of a fourth-order Runge-Kutta integration (correction 35 of the discovery engine).

Correction 34 bounded the explicit step by |a| / |l|^2 over the stable eigenvalues l = a + i b, because explicit Euler
multiplies l by |1 + l dt| and a focus that turns much faster than it decays grew at the step 0.3 / |l|. The attractor
step of the engine integrates by Runge-Kutta with a fixed number of steps, and took the same step: the time it covered
fell with the bound, by a factor of 7e-5 in one corpus model, whose motion then read a drift instead of a torus. Runge-Kutta
multiplies l by R(l dt) = 1 + z + z^2 / 2 + z^3 / 6 + z^4 / 24, and |R(i y)|^2 = 1 - y^6 / 72 + y^8 / 576 stays below 1 for
y < 2 sqrt(2), so at |l| dt = 0.3 it keeps the focus stable. rate_dt is the step bounded by the fastest rate alone."""
from __future__ import annotations

import numpy as np
import pytest

from fieldbridge.memory import analysis as an
from fieldbridge.memory import spec


def _spec(name, variables, parameters, drift, control, control_range, scale=1.0):
    return {"schema": "fieldbridge-memory/1", "kind": "equations", "name": name, "question": "what does it store?",
            "assumptions": ["test model"], "provenance": {"source": "test"},
            "carrier": {"kind": "euclid", "name": "state", "variables": variables, "scale": scale},
            "parameters": parameters, "drift": drift, "control": {"name": control, "range": control_range},
            "noise": 0.01, "dt": 0.01}


# a linear stable focus that decays at 0.3 and turns at w, with w up to 4e4 over its control range
FOCUS = _spec("fast-turning stable focus", ["x", "y"], {"w": 3e4}, {"x": "-0.3*x - w*y", "y": "w*x - 0.3*y"}, "w",
              [2e4, 4e4])
# two real rates, 1 and 50: no spiral
NODE = _spec("stable node", ["x", "y"], {"k": 50.0}, {"x": "-x", "y": "-k*y"}, "k", [40.0, 60.0])


def _rk4_factor(z):
    return 1 + z + z ** 2 / 2 + z ** 3 / 6 + z ** 4 / 24


def test_the_runge_kutta_step_is_set_by_the_fastest_rate_alone():
    real = spec.load(FOCUS)
    fastest = abs(complex(-0.3, 4e4))
    assert an.rate_dt(real) == pytest.approx(0.3 / fastest, rel=1e-6)
    assert an.safe_dt(real) == pytest.approx(0.3 / fastest ** 2, rel=1e-6)  # |a| / |l|^2 at the end of the range
    assert an.safe_dt(real) < 1e-4 * an.rate_dt(real)


def test_runge_kutta_keeps_the_focus_decaying_where_euler_would_not():
    real = spec.load(FOCUS)
    dt = an.rate_dt(real)
    for w in (2e4, 3e4, 4e4):
        z = complex(-0.3, w) * dt
        assert abs(1 + z) > 1  # Euler grows at this step
        assert abs(_rk4_factor(z)) < 1  # Runge-Kutta does not
    # and the integration itself: 20,000 steps at w = 3e4 from a unit start, the amplitude never grows
    J = np.array([[-0.3, -3e4], [3e4, -0.3]])
    q, amplitude = np.array([1.0, 0.0]), [1.0]
    f = lambda v: J @ v
    for _ in range(20_000):
        k1 = f(q)
        k2 = f(q + 0.5 * dt * k1)
        k3 = f(q + 0.5 * dt * k2)
        k4 = f(q + dt * k3)
        q = q + dt * (k1 + 2 * k2 + 2 * k3 + k4) / 6
        amplitude.append(float(np.linalg.norm(q)))
    assert np.all(np.diff(amplitude) <= 0) and amplitude[-1] < 1


def test_the_two_steps_agree_where_no_spiral_binds():
    real = spec.load(NODE)
    assert an.rate_dt(real) == an.safe_dt(real) == pytest.approx(0.3 / 60.0)


@pytest.mark.parametrize("first", ["rate", "safe"])
def test_a_finer_resolution_divides_both_steps(first):
    real = spec.load(FOCUS)
    rate, euler = an.rate_dt(spec.load(FOCUS)), an.safe_dt(spec.load(FOCUS))
    (an.rate_dt if first == "rate" else an.safe_dt)(real)
    an.refine_steps(real, 2)
    assert an.rate_dt(real) == pytest.approx(rate / 2) and an.safe_dt(real) == pytest.approx(euler / 2)
