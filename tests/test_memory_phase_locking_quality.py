"""The phase-locking derivation chooses its operating point deterministically and reports the Adler constant only when
the fit is a measurement (phase_locking.operating_point, adler_quality). Found on the stirred Belousov-Zhabotinsky model
of the corpus cells, whose constant ranged from 0.95 to 405 across seeds: the random starts of the scan decided which
end of the oscillating interval was its middle, and the fit at a slow relaxation oscillation had no gate."""
from __future__ import annotations

import numpy as np

from fieldbridge.memory import phase_locking as pl
from fieldbridge.memory import spec as mspec

VDP = {
    "schema": "fieldbridge-memory/1", "kind": "equations", "name": "van der Pol below onset (test)", "field": "test",
    "question": "test", "assumptions": [], "provenance": {"source": "test", "contributor": "test"},
    "carrier": {"kind": "euclid", "name": "plane", "variables": ["x", "v"], "scale": 1},
    "parameters": {"mu": -0.5},
    "drift": {"x": "v", "v": "mu*(1 - x**2)*v - x"},
    "control": {"name": "mu", "range": [-1.0, 1.0]}, "noise": 0.0, "dt": 0.01,
}


def _fit(rows_inside_residual=1e-9, slip=1.0, centre=0.0, half_width=1.0, eps=0.1):
    rows = [{"side": "inside", "nu": nu, "rate": 0.5, "residual": rows_inside_residual} for nu in pl.NU_INSIDE]
    rows += [{"side": "outside", "nu": nu, "frequency": slip} for nu in pl.NU_OUTSIDE]
    return {"eps": eps, "rows": rows, "centre": centre, "half_width": half_width}


def test_the_adler_fit_is_a_measurement_only_when_every_part_of_it_succeeded():
    clean = [_fit(eps=0.2), _fit(eps=0.1)]
    q = pl.adler_quality(clean, const=1.0)
    assert q["measured"] and q["quality_notes"] == []
    bad = pl.adler_quality([_fit(eps=0.2), _fit(eps=0.1, rows_inside_residual=1e-2)], const=1.0)
    assert not bad["measured"] and any("no locked orbit" in n for n in bad["quality_notes"])
    bad = pl.adler_quality([_fit(eps=0.2), _fit(eps=0.1, slip=float("nan"))], const=1.0)
    assert not bad["measured"] and any("no slip frequency" in n for n in bad["quality_notes"])
    bad = pl.adler_quality([_fit(eps=0.2), _fit(eps=0.1, centre=-1.2)], const=1.0)
    assert not bad["measured"] and any("centre" in n for n in bad["quality_notes"])
    bad = pl.adler_quality([_fit(eps=0.2, half_width=5.0), _fit(eps=0.1, half_width=1.0)], const=-3.0)
    assert not bad["measured"] and any("extrapolation" in n for n in bad["quality_notes"])


def test_the_operating_point_is_the_scanned_value_nearest_the_middle_of_the_oscillating_run_and_reuses_its_cycle():
    real = mspec.load(VDP)
    p1, cyc, continued, scan = pl.operating_point(real, np.random.default_rng(0))
    assert continued and cyc is not None
    osc = [row["value"] for row in scan if row["oscillates"]]
    assert min(osc) >= 0.0 and len(osc) >= 3         # no cycle for mu < 0; at mu = 0 the linear motion is periodic
    mid = 0.5 * (min(osc) + max(osc))
    assert abs(p1 - min(osc, key=lambda v: abs(v - mid))) < 1e-9
    assert cyc["period"] > 0 and p1 in [row["value"] for row in scan]   # the cycle found by the scan at that value
    p2, _, _, _ = pl.operating_point(real, np.random.default_rng(7))
    assert abs(p2 - p1) < 1e-9                         # and the same with other random starts
