"""Formula (9) of the completed write law (completion.consistent_selection): version 14's expansion with the exact
remainders of the push and of the outside fluctuation. It is off by default (CONSISTENT_EXPANSION); these tests turn it
on. The samples per seed are lowered in conftest.py; the six seeds are kept."""
from __future__ import annotations

from math import sqrt

import numpy as np

from fieldbridge.memory import analysis as an
from fieldbridge.memory import completion as cm
from fieldbridge.memory import spec as mspec



def _linear_path(N=40, dt=0.1, D=1e-3, bias=0.1):
    """A swept linear passage in two variables: x's rate crosses zero, y relaxes. Biases +bias and -bias along x."""
    n = 2
    ts = dt * np.arange(N + 1)
    Js = np.array([[[0.2 * t - 0.4, 0.0], [0.0, -1.0]] for t in ts])
    dts = np.diff(ts)
    Id = np.eye(n)
    A = np.array([Id + 0.5 * dts[i] * (Js[i] + Js[i + 1]) + 0.5 * dts[i] ** 2 * Js[i + 1] @ Js[i] for i in range(N)])
    Cd = np.zeros((N + 1, n, n))
    for i in range(N):
        Cd[i + 1] = A[i] @ Cd[i] @ A[i].T + 2 * D * dts[i] * Id
    lam = np.zeros((N + 1, n))
    lam[N] = [1.0, 0.0]
    for i in range(N - 1, -1, -1):
        lam[i] = A[i].T @ lam[i + 1]
    var = float(lam[N] @ Cd[N] @ lam[N])
    w = np.zeros(N + 1)
    w[:-1] += 0.5 * dts
    w[1:] += 0.5 * dts
    c = np.einsum("tij,tj->ti", Cd, lam)
    m = np.zeros((N + 1, 2, n))
    m[:, 0, 0] = bias * ts / ts[-1]
    m[:, 1, 0] = -bias * ts / ts[-1]
    mus = list(m[N] @ lam[N])
    nus = [m[:, j, :] - c * (mus[j] / var) for j in range(2)]
    return dict(A=A, dts=dts, D=D, w=w, lam=lam, c=c, var=var, m=m, nus=nus, mus=mus)


def test_the_switch_is_off_by_default():
    assert cm.CONSISTENT_EXPANSION is False


def test_without_curvature_formula_9_is_the_expansion_and_the_overlap_is_carried_with_it():
    p = _linear_path()
    sigma = sqrt(p["var"])
    rows = [{"z0": mu / sigma, "Delta": 0.0, "delta2": 0.0, "third_order": 0.0} for mu in p["mus"]]
    N1, n = p["lam"].shape
    flat = np.zeros((N1, n, n))
    cs = cm.consistent_selection(p["A"], p["dts"], p["D"], p["w"], flat, p["lam"], p["c"], p["var"], p["m"], p["nus"],
                                 p["mus"], 0.0, rows, seeds=(1, 2, 3), n_mc=300)
    assert cs["applied"]
    # no push, no outside rate: the remainder vanishes and both estimators give the first-order statistic
    assert np.allclose(cs["z"], [r["z0"] for r in rows], atol=1e-9)
    for rec in cs["pairs"]:
        for q in ("half_difference", "mean"):
            assert np.allclose(rec[q]["values"]["cv"], rec[q]["values"]["direct"], atol=1e-9)
    # with curvature, the overlap of the push with the third-order mean (O2) is carried through the auxiliary push of
    # the critical part; a port that dropped that push gave O2 = 0 and counted the overlap twice
    curved = np.tile(np.array([[0.1, 0.03], [0.03, 0.05]]), (N1, 1, 1))
    nu = np.stack(p["nus"], axis=1)
    S = cm._outside_fluctuation_sums(p["A"], curved, p["lam"], p["c"], p["var"], p["m"], nu, p["w"], p["D"], p["dts"], 200, 1)
    assert np.max(np.abs(S["O2"])) > 0


DAMPING = {   # a mass whose damping a very slow independent variable modulates (held-out set 4 of the completed write law,
              # KnowledgeParser law_v23): y relaxes at 0.1 < sqrt(0.02) and is decoupled from the choice, but its
              # fluctuation changes the rate at which the critical coordinate grows
    "schema": "fieldbridge-memory/1", "kind": "equations", "name": "mass with a slowly modulated damping (test)",
    "field": "test", "question": "test", "assumptions": [], "provenance": {"source": "test", "contributor": "test"},
    "carrier": {"kind": "euclid", "name": "x, v, y", "variables": ["x", "v", "y"], "scale": 3.0},
    "parameters": {"e": -1.0, "g": 1.5, "c": 2.5, "k": 0.1},
    "drift": {"x": "v", "v": "-(g + c*y)*v + e*x - x**3", "y": "-k*y"},
    "noise": 0.001, "control": {"name": "e", "range": [-1.0, 1.0]}, "material": [],
}


def test_formula_9_gives_the_response_of_the_full_dynamics_on_a_slowly_modulated_damping(monkeypatch):
    monkeypatch.setattr(cm, "CONSISTENT_EXPANSION", True)
    real = mspec.load(DAMPING)
    rng = np.random.default_rng(0)
    ev = an.locate_writes(real, rng)
    nf = an.normal_form(real, ev[0]["q"], real.control, ev[0]["v"])
    law = cm.law_at(real, nf, rng, n_samples=2000, n_traj=1000, certify=False, floors=(1.0,))
    assert law.get("consistent"), law.get("status")
    (p,) = law["consistent"]["per_floor"]
    assert p["in_regime"] and law.get("next_order_not_small") is None
    # direct simulation (400000 trajectories, KnowledgeParser law_v23 set 4) at one noise floor of bias:
    # response constant 1.3051 +- 0.0020, offset -0.0019 +- 0.0015; the gate of the held-out sets is 4 se + 0.02
    gate = lambda e, ea: 4 * sqrt(e ** 2 + ea ** 2) + 0.02          # noqa: E731
    assert abs(p["response_constant"] - 1.3051) <= gate(p["response_constant_error"], 0.0020)
    assert abs(p["offset"] + 0.0019) <= gate(p["offset_error"], 0.0015)
    # and the value of the six-seed estimate at 16000 samples per seed, 1.3047; a sum that loses the overlap of the push
    # with the third-order mean gives 1.288 here
    assert abs(p["response_constant"] - 1.3047) < 0.005
