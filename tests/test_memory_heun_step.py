"""The step of Heun's rule (analysis.heun_dt), which integrates the moments of the completion's linear stage.

Heun's rule multiplies an eigenvalue l by R(l dt) = 1 + l dt + (l dt)^2 / 2, with |R(i y)|^2 = 1 + y^4 / 4 on the
imaginary axis: at a focus that turns much faster than it decays it needs a bound, as explicit Euler does, but a far
milder one. The criterion is that of safe_dt, half the decay kept, taken at the fastest rate of the covariance, 2 l."""
from __future__ import annotations

from pathlib import Path

import numpy as np

from fieldbridge.memory import analysis as an
from fieldbridge.memory import spec

QUANTUM = Path(__file__).resolve().parents[1] / "examples" / "quantum" / "semiclassical"


def _focus(a: float, b: float, dt: float = 0.01):
    return spec.load({"schema": "fieldbridge-memory/1", "kind": "equations", "name": "linear focus", "question": "test",
                      "assumptions": ["test model"], "provenance": {"source": "test"},
                      "carrier": {"kind": "euclid", "name": "plane", "variables": ["x", "y"], "scale": 1.0},
                      "parameters": {"a": a, "b": b}, "drift": {"x": "a*x - b*y", "y": "b*x + a*y"},
                      "noise": 0.01, "dt": dt})


def test_the_heun_step_follows_the_small_step_law_at_a_weakly_damped_focus():
    """For small steps the criterion at 2 l reads 4 |l|^4 dt^4 <= 2 |a| dt (1 + 4 |l|^2 dt^2), dt^3 ~ |a| / (2 |l|^4):
    for a = -0.02, b = 40 it gives 2.5e-3, between Euler's |a| / |l|^2 = 1.25e-5 and the rate step 0.3 / |l| = 7.5e-3."""
    a, b = -0.02, 40.0
    real = _focus(a, b)
    lam = complex(a, b)
    h = an.heun_dt(real)
    assert abs(h / (abs(a) / (2 * abs(lam) ** 4)) ** (1 / 3) - 1) < 0.05
    assert an.safe_dt(real) < h < an.rate_dt(real)
    z = 2 * lam * h
    assert abs(1 + z + z * z / 2) ** 2 <= 1 + z.real + 1e-12        # half the decay kept at the step
    z = 2 * lam * 1.05 * h
    assert abs(1 + z + z * z / 2) ** 2 > 1 + z.real                 # and not beyond it


def test_a_strongly_damped_or_real_spectrum_leaves_the_rate_step():
    assert an.heun_dt(_focus(-5.0, 1.0)) == an.rate_dt(_focus(-5.0, 1.0))


def test_the_linear_stage_of_the_magnon_kerr_element_keeps_its_step():
    """The focus bound of safe_dt shortened the step of the completion's linear stage on the magnon-Kerr element by 48
    (one test went from 9 minutes to 5 hours); Heun's bound lies just above the step it takes, half the rate step."""
    real = spec.load(QUANTUM / "magnon_kerr.json")
    assert an.safe_dt(real) < 0.05 * an.rate_dt(real)
    assert an.heun_dt(real) >= 0.5 * an.rate_dt(real)


def test_refining_divides_every_step():
    real = _focus(-0.02, 40.0)
    rate, euler, heun = an.rate_dt(real), an.safe_dt(real), an.heun_dt(real)
    an.refine_steps(real, 2.0)
    assert np.allclose([an.rate_dt(real), an.safe_dt(real), an.heun_dt(real)], [rate / 2, euler / 2, heun / 2])
