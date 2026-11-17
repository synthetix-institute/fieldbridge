"""A write point whose critical eigenvalue has no dual direction along the mode.

The normal form divides the left eigenvector of the critical eigenvalue by its overlap with the right one. Where the
eigenvalue is defective (a Jordan block) the one left and the one right eigenvector are orthogonal, the overlap is zero,
and every coefficient was NaN: the kind fell through every comparison to "subcritical pitchfork" with critical
dimension 0. Such a point is now a degenerate critical mode whose dimension is the number of repeated eigenvalues; there
are no coefficients along a line, and they stay NaN."""
from __future__ import annotations

import numpy as np

from fieldbridge.memory import analysis as an
from fieldbridge.memory import spec


def _spec(name, variables, parameters, drift, control, control_range):
    return {"schema": "fieldbridge-memory/1", "kind": "equations", "name": name, "question": "what does it store?",
            "assumptions": ["test model"], "provenance": {"source": "test"},
            "carrier": {"kind": "euclid", "name": "state", "variables": variables, "scale": 1.0},
            "parameters": parameters, "drift": drift, "control": {"name": control, "range": control_range},
            "noise": 0.01, "dt": 0.01}


CHAIN = _spec("third-order linear chain", ["x", "y", "z"], {"k": 0.0}, {"x": "y", "y": "z", "z": "k*x"}, "k",
              [-1.0, 1.0])
PITCHFORK = _spec("pitchfork with a relaxing variable", ["x", "y"], {"k": 0.0}, {"x": "k*x - x**3", "y": "-y"}, "k",
                  [-1.0, 1.0])


def test_a_jordan_block_is_a_degenerate_critical_mode_of_its_full_dimension():
    real = spec.load(CHAIN)
    with np.errstate(all="ignore"):
        nf = an.normal_form(real, np.zeros(3), "k", 0.0)
    assert nf["kind"].startswith("degenerate critical mode")
    assert nf["critical_dimension"] == 3 and nf["defective"]
    assert np.isnan(nf["a3"]) and np.isnan(nf["kappa_slope"])


def test_a_simple_zero_keeps_its_coefficients():
    real = spec.load(PITCHFORK)
    nf = an.normal_form(real, np.zeros(2), "k", 0.0)
    assert nf["kind"].startswith("supercritical pitchfork") and nf["critical_dimension"] == 1
    assert "defective" not in nf and np.isfinite(nf["a3"]) and nf["a3"] < 0
