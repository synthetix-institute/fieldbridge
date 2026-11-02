"""Why the symmetric write of the rf-SQUID phase qubit lies below pi^(1/4) in the law record (Module 15, section 4).

The swept-write check starts an element from rest, with the bias and the noise on, a time 4/sqrt(r) before the crossing. A
weakly damped phase particle rings for a time 2/gam, and at gam = 0.2 that is not short against 4/sqrt(r) = 28. Three
calculations separate this from the write itself:

  measured   the constant of `memory codiscover` at the record's 25,600 trajectories, on independent seeds and at several
             dampings of the specification;
  linear     the linear stage alone, integrated exactly (mean and covariance) from the start of the check and from an
             earlier start, with no simulation;
  full       the nonlinear dynamics at the cusp, read at the well after a hold, with the start of the check and with a
             longer run-up, and at weaker noise.

    python3 -B examples/quantum/semiclassical/survey/phase_qubit_start.py

The results are written beside this file (phase_qubit_start.json). `measured` takes about half an hour on one core and `full`
about ten minutes; --only linear takes a few seconds.
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path
from statistics import NormalDist

import numpy as np
from scipy.integrate import solve_ivp

ROOT = Path(__file__).resolve().parents[4]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from fieldbridge.memory import codiscovery, record, spec  # noqa: E402

SPEC = ROOT / "examples/quantum/semiclassical/phase_qubit_rfsquid.json"
RATE = 0.02           # the rate of swept_write_check
TARGET_ACCURACY = 0.8  # the accuracy its bias aims at
ND = NormalDist()
MEASURED = [(0.2, s) for s in (1, 2, 3, 4, 5, 6)] + [(g, s) for g in (0.1, 0.5, 1.0, 2.0) for s in (1, 2, 3)]
LINEAR_DAMPINGS = (0.05, 0.1, 0.2, 0.3, 0.5, 1.0, 2.0)
FULL = [(0.1, 4.0, 0.05), (0.1, 8.0, 0.05), (0.1, 16.0, 0.05), (0.2, 4.0, 0.05), (0.2, 8.0, 0.05), (0.2, 16.0, 0.05),
        (0.1, 4.0, 0.003), (0.1, 4.0, 0.0008)]


def left_mode(gam: float) -> np.ndarray:
    """The left critical vector at the cusp, one on the right vector (1, 0): the slow coordinate is x + v / gam."""
    return np.array([1.0, 1.0 / gam])


def linear_constant(gam: float, before: float = 4.0, after: float = 8.0, noise_on_phase: bool = True) -> float:
    """The constant of the linear stage alone, x' = v + h + xi1, v' = gam r t x - gam v + xi2, from rest a time before / sqrt(r)
    ahead of the crossing, from the mean and the covariance of the slow coordinate at after / sqrt(r). The slower
    eigenvalue crosses zero with slope r, whatever the damping."""
    t0, t1 = -before / math.sqrt(RATE), after / math.sqrt(RATE)
    w, force = left_mode(gam), np.array([1.0, 0.0])
    noise = np.diag([1.0 if noise_on_phase else 0.0, 1.0])

    def rhs(t, y):
        A = np.array([[0.0, 1.0], [gam * RATE * t, -gam]])
        m, S = y[:2], y[2:].reshape(2, 2)
        return np.concatenate([A @ m + force, (A @ S + S @ A.T + 2.0 * noise).ravel()])

    sol = solve_ivp(rhs, (t0, t1), np.zeros(6), rtol=1e-9, atol=1e-11, method="LSODA")
    m, S = sol.y[:2, -1], sol.y[2:, -1].reshape(2, 2)
    return float(w @ m / math.sqrt(w @ S @ w) * math.sqrt(w @ noise @ w) * RATE ** 0.25 / (w @ force))


def full_constant(gam: float, before: float, gamma: float, n: int = 40000, dt: float = 0.01, seed: int = 1) -> dict:
    """The nonlinear dynamics at the cusp, v' = -gam v - ell x + sin x, ell = 1 - gam r t, from rest a time before / sqrt(r)
    ahead of the crossing, at the noise of the check (noise ratio gamma) and with its bias. The side is read at the
    well after the control has been held at its end value for 6 / gam, when the velocity has decayed."""
    rng = np.random.default_rng(seed)
    w = left_mode(gam)
    ww, a3 = float(w @ w), 1.0 / (6.0 * gam)
    D = gamma * RATE / (a3 * ww)
    D_s = D * ww
    h_s = ND.inv_cdf(TARGET_ACCURACY) * math.sqrt(D_s) * RATE ** 0.25 / math.pi ** 0.25
    x, v = np.zeros(n), np.zeros(n)
    step = math.sqrt(2.0 * D * dt)
    t0, t_end = -before / math.sqrt(RATE), 4.0 / math.sqrt(RATE)
    sweep, hold = int(round((t_end - t0) / dt)), int(round(6.0 / gam / dt))
    for i in range(sweep + hold):
        ell = 1.0 - gam * RATE * (t0 + min(i, sweep) * dt)
        dx = (v + h_s) * dt + step * rng.standard_normal(n)
        dv = (-gam * v - ell * x + np.sin(x)) * dt + step * rng.standard_normal(n)
        x, v = x + dx, v + dv
    p = min(max(float(np.mean(x > 0)), 1e-4), 1.0 - 1e-4)
    scale = math.sqrt(D_s) * RATE ** 0.25 / h_s
    return {"damping": gam, "start_before_crossing_sqrt_r": before, "noise_ratio": gamma, "trajectories": n,
            "constant": ND.inv_cdf(p) * scale, "stderr": math.sqrt(p * (1.0 - p) / n) / ND.pdf(ND.inv_cdf(p)) * scale}


def measured_constant(gam: float, seed: int, n_traj: int = 25600, replicates: int = 8) -> dict:
    data = json.loads(SPEC.read_text(encoding="utf-8"))
    data["parameters"]["gam"] = gam
    row = codiscovery.codiscover([spec.load(data)], seed, target="symmetric-write", check_law=True, n_traj=n_traj,
                                 replicates=replicates)["rows"][0]
    const = row.get("law_constant") or {}
    return {"damping": gam, "seed": seed, "word": row.get("word"), "constant": const.get("constant"),
            "stderr": const.get("stderr"), "offset": const.get("offset")}


def weighted_mean(rows: list) -> dict:
    weights = np.array([1.0 / r["stderr"] ** 2 for r in rows])
    values = np.array([r["constant"] for r in rows])
    mean = float((weights * values).sum() / weights.sum())
    return {"runs": len(rows), "constant": mean, "stderr": float(weights.sum() ** -0.5),
            "chi_squared": float((weights * (values - mean) ** 2).sum())}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--only", choices=("linear", "full", "measured"), default=None)
    ap.add_argument("--out", type=Path, default=Path(__file__).with_suffix(".json"))
    args = ap.parse_args()
    out = {"implementation_sha256": record.implementation_hash(), "versions": record.versions(),
           "law": math.pi ** 0.25, "rate": RATE, "start_of_the_check": 4.0 / math.sqrt(RATE)}
    if args.only in (None, "linear"):
        out["linear"] = [{"damping": g, "from_rest_4": linear_constant(g), "from_rest_8": linear_constant(g, before=8.0),
                          "from_rest_16": linear_constant(g, before=16.0),
                          "no_noise_on_the_phase_8": linear_constant(g, before=8.0, noise_on_phase=False)}
                         for g in LINEAR_DAMPINGS]
    if args.only in (None, "full"):
        out["full"] = [full_constant(*p) for p in FULL]
    if args.only in (None, "measured"):
        rows = [measured_constant(g, s) for g, s in MEASURED]
        out["measured"] = rows
        out["mean_by_damping"] = {str(g): weighted_mean([r for r in rows if r["damping"] == g])
                                  for g in sorted({r["damping"] for r in rows})}
    args.out.write_text(json.dumps(out, indent=1) + "\n", encoding="utf-8")
    print(json.dumps({k: v for k, v in out.items() if k in ("linear", "mean_by_damping")}, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
