"""The Jacobian at a state on the edge of the set where the drift is finite (correction 23 of the discovery engine).

The Jacobian is a central difference with a stencil of 1e-6. At a state on the edge of an orthant, under a non-integer
power p^(5/2), the stencil evaluates the drift at p = -1e-6, where it is not defined, and the eigenvalues of the
relaxed state p = 0 could not be taken: the card was not derived. An entry whose central difference is not finite, while
the drift at the state and on one side is, is now the one-sided difference from that side. Where every value is finite
the Jacobian is unchanged."""
from __future__ import annotations

import numpy as np
import pytest

from fieldbridge.memory import analysis as an
from fieldbridge.memory import spec


def _spec(name, variables, parameters, drift, control, control_range, carrier="euclid"):
    return {"schema": "fieldbridge-memory/1", "kind": "equations", "name": name, "question": "what does it store?",
            "assumptions": ["test model"], "provenance": {"source": "test"},
            "carrier": {"kind": carrier, "name": "state", "variables": variables, "scale": 1.0},
            "parameters": parameters, "drift": drift, "control": {"name": control, "range": control_range},
            "noise": 0.01, "dt": 0.01}


# an autoactivating gene with Hill exponent 3/2: p = 0 is stable at every V (rate 1); a fold creates the upper state
GENE = _spec("gene with Hill exponent 3/2", ["p"], {"V": 1.0}, {"p": "-p + V*p**1.5/(0.5**1.5 + p**1.5)"}, "V",
             [0.5, 1.5], carrier="orthant")
SMOOTH = _spec("two coupled cubic relaxations", ["x", "y"], {"a": 0.3},
               {"x": "a*x - x**3 + 0.2*y", "y": "-y + x**2"}, "a", [-0.5, 0.5])
ROOT = _spec("relaxation to zero under a square root", ["x"], {"g": 1.0}, {"x": "-g*sqrt(x)"}, "g", [0.5, 1.5],
             carrier="orthant")


def test_the_state_on_the_edge_has_its_one_sided_rate():
    real = spec.load(GENE)
    # the one-sided difference of p^(3/2) over eps carries V eps^(1/2) / K^(3/2) = 2.8e-3 of its own
    J = an.jacobian(real, np.array([0.0]))
    assert np.all(np.isfinite(J)) and J[0, 0] == pytest.approx(-1.0 + 1e-3 / 0.5**1.5, rel=1e-6)
    assert an.kappa_spectrum(real, np.array([0.0]))[0].real == pytest.approx(1.0, abs=5e-3)


def test_the_gene_on_the_orthant_stores_its_states():
    real = spec.load(GENE)
    states, counts, spectra, unconverged, continuum = an.stored_states(real, np.random.default_rng(0), 48)
    assert sorted(round(float(s[0]), 2) for s in states)[0] == pytest.approx(0.0, abs=1e-2)
    assert len(states) == 2 and not continuum
    assert all(np.all(np.isfinite(k)) for k in spectra)


def test_the_jacobian_is_unchanged_where_every_value_is_finite():
    real = spec.load(SMOOTH)
    eps = 1e-6
    for q in np.random.default_rng(1).normal(size=(20, 2)):
        E = eps * np.eye(2)
        plain = ((real.F(q[None] + E) - real.F(q[None] - E)) / (2 * eps)).T
        assert np.array_equal(an.jacobian(real, q), plain)


def test_an_infinite_one_sided_slope_is_a_large_finite_rate():
    # -g sqrt(x) reaches x = 0 in a finite time; its one-sided slope there is infinite, and the difference gives
    # -g / sqrt(eps), a large restoring rate
    real = spec.load(ROOT)
    J = an.jacobian(real, np.array([0.0]))
    assert J[0, 0] == pytest.approx(-1.0 / np.sqrt(1e-6), rel=1e-9)


def test_a_drift_that_is_not_finite_at_the_state_stays_not_finite():
    # outside the set itself there is no side to take the difference from
    real = spec.load(_spec("square root on a line", ["x"], {"g": 1.0}, {"x": "g - sqrt(x)"}, "g", [0.5, 1.5]))
    with np.errstate(all="ignore"):
        assert not np.any(np.isfinite(an.jacobian(real, np.array([-1.0]))))
