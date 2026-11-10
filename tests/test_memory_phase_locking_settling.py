"""The search for a cycle stops a start that has settled on a stable state even where the error of the integration keeps
the drift above the settling threshold (phase_locking.find_cycle, _at_stable_state). Found on the spin ensemble in a
cavity and the bright-state resonator of the hardware survey: at a strong drive their stable state is a fast spiral, the
drift stayed at 1e-5 to 1e-3 of its median over the carrier while the state was within 1e-7 of the scale of the fixed
point, and every start ran to t_max, so that the phase-locking derivation took 4 and 20 hours to find no oscillation."""
from __future__ import annotations

from pathlib import Path

import numpy as np

from fieldbridge.memory import analysis as an
from fieldbridge.memory import phase_locking as pl
from fieldbridge.memory import spec as mspec

SPIN = Path(__file__).resolve().parents[1] / "examples/quantum/semiclassical/spin_ensemble_cavity.json"

SPIRAL = {
    "schema": "fieldbridge-memory/1", "kind": "equations", "name": "stable and unstable focus (test)", "field": "test",
    "question": "test", "assumptions": [], "provenance": {"source": "test", "contributor": "test"},
    "carrier": {"kind": "euclid", "name": "plane", "variables": ["x", "y"], "scale": 1},
    "parameters": {"a": -0.5, "w": 3.0},
    "drift": {"x": "a*(x - 2) - w*y", "y": "w*(x - 2) + a*y"},
    "control": {"name": "a", "range": [-1.0, 1.0]}, "noise": 0.0, "dt": 0.01,
}


def test_a_start_at_the_error_floor_of_a_fast_spiral_stops_after_the_first_chunk(monkeypatch):
    real = mspec.load(SPIN)
    over = {real.control: 18.0}
    calls = []
    solve = pl._solve

    def counted(*a, **k):
        sol = solve(*a, **k)
        calls.append(sol.y[:, -1])
        return sol

    monkeypatch.setattr(pl, "_solve", counted)
    rng = np.random.default_rng(1)
    assert pl.find_cycle(real, over, rng, n_starts=1) is None
    probe = real.carrier.sample(np.random.default_rng(1), 16)
    f_ref = float(np.nanmedian(np.linalg.norm(real.F(probe, **over), axis=-1)))
    drift = float(np.linalg.norm(real.F(calls[-1], **over)))
    assert drift > pl.SETTLED_DRIFT * f_ref  # the drift alone would not have stopped this start
    assert len(calls) <= 3, len(calls)       # it ran 30 chunks (to t_max) before


def test_only_a_stable_state_next_to_the_end_of_a_chunk_counts_as_settled():
    real = mspec.load(SPIRAL)
    scale = an.base_length(real)
    near = np.array([2.0 + 1e-8 * scale, 0.0])
    assert pl._at_stable_state(real, {"a": -0.5}, near)
    assert not pl._at_stable_state(real, {"a": 0.5}, near)                      # the same point, an unstable focus
    assert not pl._at_stable_state(real, {"a": -0.5}, np.array([2.0 + 1e-3 * scale, 0.0]))  # not yet reached
