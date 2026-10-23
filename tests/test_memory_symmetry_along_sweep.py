"""A symmetry that holds at the write point and not along the sweep (codiscovery.symmetry_along_sweep): the law check
then tests both signs of the bias, so that the offset the passage writes is separated from the response. The case is
the nonsymmetric pitchfork of arXiv:1804.03264 read in the records of the corpus cells: at epsilon = 0 the drift has
the exchange symmetry, and the term epsilon^2 breaks it as soon as the control moves."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from fieldbridge.memory import codiscovery as cd
from fieldbridge.memory import spec as mspec
from fieldbridge.memory.construct import sweep_ends

ROOT = Path(__file__).resolve().parents[1]
NONSYMMETRIC = {
    "schema": "fieldbridge-memory/1", "kind": "equations", "name": "nonsymmetric pitchfork (test)", "field": "mathematics",
    "question": "test", "assumptions": [], "provenance": {"source": "arXiv:1804.03264", "contributor": "test"},
    "carrier": {"kind": "euclid", "name": "plane", "variables": ["x", "y"], "scale": 1},
    "parameters": {"epsilon": 0.0},
    "drift": {"x": "y**2 - (epsilon + 1)*y - x", "y": "x**2 - (epsilon + 1)*x - y + epsilon**2"},
    "control": {"name": "epsilon", "range": [-0.1, 0.1]}, "noise": 0.01, "dt": 0.01,
}


def test_the_sweep_of_the_write_check_ends_beyond_the_write_point_on_the_unstable_side():
    nf = {"value": 0.0, "kappa_slope": -1.0}   # the critical rate grows with the control
    lo, hi = sweep_ends(nf, rate=0.02)
    assert abs(lo + 4 * np.sqrt(0.02)) < 1e-12 and abs(hi - 4 * np.sqrt(0.02)) < 1e-12
    lo2, hi2 = sweep_ends({"value": 1.0, "kappa_slope": 2.0}, rate=0.02)   # the rate falls with the control: swept down
    assert hi2 < 1.0 < lo2 and abs((lo2 - 1.0) - 2 * np.sqrt(0.02)) < 1e-12


def test_a_symmetry_broken_along_the_sweep_makes_the_law_check_resolve_the_offset():
    real = mspec.load(NONSYMMETRIC)
    res = cd.derive_symmetric_write(real, np.random.default_rng(0), check_law=True, n_traj=400)
    assert res["status"] == "reached" and res["word"].startswith("S")
    S = next(s for s in res["steps"] if s["letter"] == "S")
    assert S["along_sweep"] is False
    assert any("not along the sweep" in n for n in res.get("notes", []))
    if res["law_constant"]:   # a side can saturate at this bias; then the law is reported as not measured
        L = next(s for s in res["steps"] if s["letter"] == "L")
        assert L["signs"] == [1, -1]
        lc = res["law_constant"]
        assert lc["offset_resolved"] and abs(lc["offset"]) > 0.5
        assert abs(lc["constant"] - np.pi ** 0.25) < 4 * lc["stderr"] + 0.02
    else:
        assert any("every trajectory ended on one side" in n for n in res["notes"])


def test_the_symmetric_pitchfork_keeps_its_symmetry_along_the_sweep_and_one_sign():
    pf = mspec.load(json.loads((ROOT / "examples" / "memory" / "pitchfork.json").read_text(encoding="utf-8")))
    res = cd.derive_symmetric_write(pf, np.random.default_rng(0), check_law=True, n_traj=100)
    assert next(s for s in res["steps"] if s["letter"] == "S")["along_sweep"] is True
    assert next(s for s in res["steps"] if s["letter"] == "L")["signs"] == [1]
