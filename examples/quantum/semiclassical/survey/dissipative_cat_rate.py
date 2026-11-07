"""Why the symmetric write of the dissipative cat qubit lies below pi^(1/4) in the law record (Module 15, section 4).

The quadrature that the write does not choose relaxes at kappa2 A + kappa1 / 2 = 0.05 at the write point, which is less than
sqrt(r) = 0.14 at the rate of the check, and it loses stability itself at A = -0.025, only 0.05 below the write point. At the
rate of the check the sweep starts at A = -0.54, beyond that second pitchfork. The constant is measured here at slower sweeps,
at which the start moves toward the stable side and the second rate exceeds sqrt(r).

    python3 -B examples/quantum/semiclassical/survey/dissipative_cat_rate.py

The results are written beside this file (dissipative_cat_rate.json), in about half an hour on one core.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[4]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from fieldbridge.memory import codiscovery, construct, record, spec  # noqa: E402

RATES = (0.02, 0.005, 0.00125, 0.0003125, 0.00015)
SEEDS = (1, 2)


def write_point():
    """The realization and the normal form that the symmetric-write derivation hands to its law check."""
    passed = []
    real_check = codiscovery.swept_write_check
    codiscovery.swept_write_check = lambda real, nf, rng, **kw: passed.append((real, nf)) or {}
    try:
        codiscovery.derive_symmetric_write(spec.load(ROOT / "examples/quantum/semiclassical/dissipative_cat.json"),
                                           np.random.default_rng(0), check_law=True, n_traj=2)
    finally:
        codiscovery.swept_write_check = real_check
    return passed[-1]


def main() -> int:
    real, nf = write_point()
    rows = []
    for rate in RATES:
        for seed in SEEDS:
            law = construct.swept_write_check(real, nf, np.random.default_rng(seed), rate=rate, n_traj=12800, replicates=8)
            const = codiscovery.law_constant(law)
            rows.append({"rate": rate, "seed": seed, "sqrt_rate": rate ** 0.5, "start": construct.sweep_ends(nf, rate)[0],
                         "constant": const.get("constant"), "stderr": const.get("stderr")})
    means = {}
    for rate in RATES:
        part = [r for r in rows if r["rate"] == rate]
        weights = np.array([1.0 / r["stderr"] ** 2 for r in part])
        means[str(rate)] = {"constant": float(sum(w * r["constant"] for w, r in zip(weights, part)) / weights.sum()),
                            "stderr": float(weights.sum() ** -0.5)}
    out = {"implementation_sha256": record.implementation_hash(), "versions": record.versions(), "write_point": nf["value"],
           "second_rate": nf.get("second_rate"), "other_pitchfork": -nf["value"], "rows": rows, "mean_by_rate": means}
    Path(__file__).with_suffix(".json").write_text(json.dumps(out, indent=1) + "\n", encoding="utf-8")
    print(json.dumps(means, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
