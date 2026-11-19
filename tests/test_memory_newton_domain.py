"""Newton and the least-squares steps of the write-point analysis outside the set where the drift is finite.

The damped Newton of the analysis solves J s = F by least squares. Where the drift or its central-difference Jacobian is
not finite at an iterate (a start of a Euclidean carrier at a negative value under a square root, a stencil that
reaches below zero, a control value that divides the drift), least squares raised "SVD did not converge in Linear Least
Squares", and the card was not derived. Such an iterate now ends the search as not converged, as a search that
diverges does; a zero whose stencil leaves the set is not a converged zero, since every caller takes the eigenvalues of
its Jacobian; the continuation halves its step where its predictor is not finite; and the slow-manifold reduction of the
normal form stops there, so that its window is halved. Since correction 23 a stencil with one side outside the set takes
the one-sided difference from the other side, so a zero on the edge of the set is a converged zero."""
from __future__ import annotations

import numpy as np
import pytest

from fieldbridge.memory import analysis as an
from fieldbridge.memory import spec


def _spec(name, variables, parameters, drift, control, control_range):
    return {"schema": "fieldbridge-memory/1", "kind": "equations", "name": name, "question": "what does it store?",
            "assumptions": ["test model"], "provenance": {"source": "test"},
            "carrier": {"kind": "euclid", "name": "state", "variables": variables, "scale": 1.0},
            "parameters": parameters, "drift": drift, "control": {"name": control, "range": control_range},
            "noise": 0.01, "dt": 0.01}


ROOT = _spec("relaxation under a square root", ["x"], {"g": 1.0}, {"x": "g - sqrt(x)"}, "g", [0.5, 1.5])
EDGE = _spec("zero on the edge of the domain", ["x"], {"g": 1.0}, {"x": "-g*x + x**1.5"}, "g", [0.5, 1.5])
DIVIDED = _spec("relaxation whose control divides the drift", ["x"], {"c": 1.0}, {"x": "c*(1 - x/c)"}, "c",
                [-1.0, 1.0])
ARC = _spec("fold under a square root with a slaved variable", ["x", "y"], {"m": 0.09},
            {"x": "sqrt(0.01 - x**2) - m + 0.1*y", "y": "-y"}, "m", [0.05, 0.12])


def test_a_start_outside_the_domain_is_not_converged():
    real = spec.load(ROOT)
    with np.errstate(all="ignore"):
        q, ok = an.newton(real, np.array([-1.0]))
    assert not ok and np.all(np.isfinite(q))
    q, ok = an.newton(real, np.array([0.5]))
    assert ok and q[0] == pytest.approx(1.0)


def test_a_stencil_outside_the_domain_is_one_sided():
    # at x = 0 the drift is finite (g), but the central difference evaluates sqrt(-1e-6); since correction 23 that entry
    # is the one-sided difference from inside the domain, a large finite slope, and Newton from x = 0 reaches x = g^2
    real = spec.load(ROOT)
    with np.errstate(all="ignore"):
        q, ok = an.newton(real, np.array([0.0]))
    assert ok and q[0] == pytest.approx(1.0)


def test_a_zero_on_the_edge_of_the_domain_has_its_one_sided_rate():
    # x = 0 is a zero of -g x + x^1.5, and its Jacobian takes (-1e-6)^1.5; since correction 23 the slope there is the
    # one-sided difference, -g, so x = 0 is a converged, stable zero beside the unstable zero x = g^2
    real = spec.load(EDGE)
    with np.errstate(all="ignore"):
        q, ok = an.newton(real, np.array([0.0]))
        found = an.fixed_points(real, np.random.default_rng(0), 40)
    assert ok and q[0] == 0.0
    assert an.kappa_spectrum(real, q)[0].real == pytest.approx(1.0, abs=2e-3)
    assert 1.0 in [round(float(z[0]), 6) for z, _ in found]


def test_the_continuation_stops_where_the_control_divides_the_drift():
    # c (1 - x / c) is c - x except at c = 0 itself, where it is 0 * inf: no zero is found there, and the continuation
    # that ends at c = 0 halves its step against it until the step is spent
    real = spec.load(DIVIDED)
    with np.errstate(all="ignore"):
        assert an.fixed_points(real, np.random.default_rng(0), 10, {"c": 0.0}) == []
        branch, ev = an.track(real, "c", np.array([1.0]), 1.0, 0.0)
    assert ev is not None and ev["type"] == "vanishes" and 0.0 < ev["v"] < 1e-3
    assert all(np.isfinite(b["q"]).all() for b in branch)


def test_the_normal_form_of_a_fold_near_the_edge_of_the_domain():
    # x' = sqrt(0.01 - x^2) - m has a fold at m = 0.1, x = 0; the first window of the reduction, +-0.25, reaches
    # |x| > 0.1, where the slow-manifold correction met a Jacobian that is not finite
    real = spec.load(ARC)
    with np.errstate(all="ignore"):
        nf = an.normal_form(real, np.zeros(2), "m", 0.1)
    assert nf["kind"].startswith("saddle-node") and np.isfinite(nf["a2"]) and nf["a2"] < 0
