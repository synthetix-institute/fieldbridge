"""The slow spectrum at a write point and the dynamical critical dimension (analysis.slow_spectrum,
analysis.dynamical_critical_dimension, completion.law_at). A direction whose rate at the write point is below the square
root of the sweep rate relaxes less than once during the passage, so the write takes place in more directions than the
continuation counts as critical; found on the rf-SQUID triangle, where a weak bond asymmetry splits the two critical
rates by 0.018 against a passage scale of 0.14 and the one-dimensional law then fails."""
from __future__ import annotations

import json
import math
from math import sqrt
from pathlib import Path

import numpy as np

from fieldbridge.memory import analysis as an
from fieldbridge.memory import completion as cm
from fieldbridge.memory import spec as mspec

ROOT = Path(__file__).resolve().parents[1]
PAIR = {   # two pitchforks weakly coupled: the symmetric and antisymmetric modes cross at a - e and a + e
    "schema": "fieldbridge-memory/1", "kind": "equations", "name": "two weakly coupled pitchforks (test)", "field": "test",
    "question": "test", "assumptions": [], "provenance": {"source": "test", "contributor": "test"},
    "carrier": {"kind": "euclid", "name": "plane", "variables": ["x1", "x2"], "scale": 1},
    "parameters": {"a": 0.0, "e": 0.03},
    "drift": {"x1": "a*x1 - x1**3 + e*(x2 - x1)", "x2": "a*x2 - x2**3 + e*(x1 - x2)"},
    "control": {"name": "a", "range": [-1.0, 1.0]}, "material": ["e"], "noise": 0.002, "dt": 0.01,
}


def test_the_dynamical_dimension_counts_the_directions_too_slow_to_relax_during_the_passage():
    nf = {"critical_dimension": 1, "spectrum": [[0.0, 0.0], [-0.05, 0.0], [-2.0, 0.0]]}
    assert an.dynamical_critical_dimension(nf, rate=0.02) == 2          # 0.05 < sqrt(0.02) = 0.141
    assert an.dynamical_critical_dimension(nf, rate=0.0001) == 1        # 0.05 > sqrt(0.0001) = 0.01
    assert an.dynamical_critical_dimension({"critical_dimension": 2}, rate=0.02) == 2   # no spectrum recorded: the count stands
    assert an.dynamical_critical_dimension(nf, rate=0.02, factor=0.1) == 1


def test_a_write_point_records_its_slow_spectrum_and_second_rate():
    real = mspec.load(PAIR)
    ev = an.locate_writes(real, np.random.default_rng(0), n_starts=12)
    first = min(ev, key=lambda x: x["v"])          # the symmetric mode crosses first, at a = 0 (e cancels on it)
    nf = an.normal_form(real, first["q"], real.control, first["v"])
    assert nf["critical_dimension"] == 1 and nf["kind"].startswith("supercritical pitchfork")
    assert abs(nf["second_rate"] - 2 * 0.03) < 5e-3                   # the antisymmetric mode relaxes at 2e
    assert len(nf["spectrum"]) == 2 and abs(nf["spectrum"][0][0]) < 1e-6
    assert an.dynamical_critical_dimension(nf, rate=0.02) == 2          # 0.06 < 0.141: both directions take part
    assert an.dynamical_critical_dimension(nf, rate=0.0004) == 1        # 0.06 > 0.02: the second relaxes in time
    pf = mspec.load(json.loads((ROOT / "examples" / "memory" / "pitchfork.json").read_text(encoding="utf-8")))
    ev = an.locate_writes(pf, np.random.default_rng(0))
    nf1 = an.normal_form(pf, ev[0]["q"], pf.control, ev[0]["v"])
    assert nf1["second_rate"] is None and an.dynamical_critical_dimension(nf1, 0.02) == 1


def test_the_completion_law_uses_the_dynamical_dimension_and_reports_both():
    real = mspec.load(PAIR)
    rng = np.random.default_rng(0)
    ev = an.locate_writes(real, rng, n_starts=12)
    first = min(ev, key=lambda x: x["v"])
    nf = an.normal_form(real, first["q"], real.control, first["v"])
    law = cm.law_at(real, nf, rng, n_samples=400, n_traj=200, certify=False)
    assert law["status"] == "derived"
    assert law["critical_dimension"] == 1 and law["dynamical_critical_dimension"] == 2
    assert abs(law["second_rate"] - 0.06) < 5e-3 and "relax less than once" in law["slow_modes"]
    assert len(law["rows"]) == 8 and all(len(r["direction"]) == 2 for r in law["rows"])   # eight directions in the plane
    # the sector law of the critical plane: a distribution over the written states for every direction
    for r in law["rows"]:
        assert r["sector_law"] is not None and abs(sum(r["sector_law"]) - 1.0) < 1e-9 and len(r["sector_law"]) == len(law["stable_states"])
    assert all(r.get("cov_critical") is not None for r in law["rows"])
    # with a stricter factor the second direction counts as relaxed, and the law is the one-dimensional one again
    import fieldbridge.memory.completion as cmod
    saved = cmod.SLOW_FACTOR
    cmod.SLOW_FACTOR = 0.1
    try:
        one = cm.law_at(real, nf, np.random.default_rng(0), n_samples=200, n_traj=100, certify=False)
    finally:
        cmod.SLOW_FACTOR = saved
    assert one["dynamical_critical_dimension"] == 1 and "slow_modes" not in one and "one_direction" in one


DRIFT = {"schema": "fieldbridge-memory/1", "kind": "equations", "name": "t", "question": "q", "assumptions": ["a"],
         "closure": "c", "observable": "o", "provenance": {"origin": "test"}, "noise": 0.01}


def _model(variables, drift, parameters, control, scale=1.0):
    return mspec.load({**DRIFT, "carrier": {"kind": "euclid", "name": "x", "variables": variables, "scale": scale},
                      "drift": drift, "parameters": parameters, "control": control,
                      "material": [k for k in parameters if k != control["name"]]})


def test_a_fold_has_one_critical_direction_and_is_not_read_as_transcritical():
    """A washboard of a phase qubit folds where cos(ph) = -ell. The eigenvalue goes as a square root of the distance
    from the fold, and the other, of rate gam, is not critical."""
    ell, gam = 0.1, 0.2
    real = _model(["ph", "v"], {"ph": "v", "v": "-gam*v - ell*(ph - phe) - sin(ph)"},
                  {"ell": ell, "gam": gam, "phe": 0.0}, {"name": "phe", "range": [0.0, 8.0]}, scale=15)
    ph = 4 * math.pi - math.acos(-ell)
    phe = ph + math.sin(ph) / ell
    nf = an.normal_form(real, np.array([ph - 1e-3, 0.0]), "phe", phe)
    assert nf["critical_dimension"] == 1
    assert nf["kind"].startswith("saddle-node")


def test_a_plane_of_equilibria_of_a_dissipationless_top_gives_no_write_point():
    real = _model(["x", "y", "z"], {"x": "2*y*z*(C-B)", "y": "-2*x*z*(C-A)", "z": "2*x*y*(B-A)"},
                  {"A": 9.3, "B": 14.5, "C": 27.9}, {"name": "B", "range": [9.3, 27.9]})
    events = [{"v": 9.3, "q": np.array([-0.576, 0.374, 0.0])}, {"v": 9.3, "q": np.array([0.207, 0.446, 0.0])}]
    kept, aside = an.distinct_write_points(real, events)
    assert kept == [] and len(aside) == 2


def test_samples_of_a_line_of_states_count_once_and_mirror_states_twice():
    line = _model(["ph", "r"], {"ph": "f*sin(ph)", "r": "-r"}, {"f": 0.0}, {"name": "f", "range": [0.0, 2.0]})
    kept, aside = an.distinct_write_points(line, [{"v": 0.0, "q": np.array([p, 0.0])} for p in (0.5, 1.0, 2.0)])
    assert len(kept) == 1 and kept[0]["samples"] == 3 and aside == []
    mirror = _model(["x", "y"], {"x": "a*x - x**3", "y": "-y"}, {"a": 1.0}, {"name": "a", "range": [0.0, 2.0]})
    kept, _ = an.distinct_write_points(mirror, [{"v": 1.0, "q": np.array([s, 0.0])} for s in (1.0, -1.0)])
    assert len(kept) == 2


def test_a_direction_neutral_on_both_sides_is_not_counted_with_the_eigenvalue_that_crosses():
    """A coordinate z that never moves keeps a zero eigenvalue along the branch; the pitchfork of x crosses at a = 0."""
    real = _model(["x", "z"], {"x": "a*x - x**3", "z": "0*z"}, {"a": 1.0}, {"name": "a", "range": [-1.0, 1.0]})
    nf = an.normal_form(real, np.array([0.0, 0.3]), "a", 0.0)
    assert nf["critical_dimension"] == 1 and nf["kind"].startswith("supercritical pitchfork")
    pair = _model(["x", "y"], {"x": "a*x - y**3", "y": "a*y + x**3"}, {"a": 1.0}, {"name": "a", "range": [-1.0, 1.0]})
    assert an.normal_form(pair, np.zeros(2), "a", 0.0)["critical_dimension"] == 2
