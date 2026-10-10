"""The decision card of a population: the predictions from its equations beside, on request, stochastic runs.

- synchronization: K_c, Omega, a1, b1, kappa, mu at the protocol's couplings; with `simulate`, full units at the first
  coupling factor (their growth rate of synchrony against mu).
- collective-write: the law at every size, rate and bias of the protocol (closed form, and along the actual passage
  for collectives); with `simulate`, replicas of the stochastic collective (chemical Langevin noise) or of the units,
  with the gate |P - P_law| < 4 stderr + 0.02.
- seeded-passage: d, the rate at the end of the step and the window law; with `simulate`, passage times at every size,
  the window in the growth exponent against ln q_d and the slope of the median against ln N (law 1/2).
"""
from __future__ import annotations

import time
from typing import Dict

import numpy as np

from . import passage
from .predict import predict


def card(pop, simulate: bool = False, replicas: int = 2000, seed: int = 1) -> Dict:
    pr = predict(pop)
    out = {"name": pop.name, "field": pop.field, "source": pop.source, "path": pop.path,
           "spec_sha256": pop.spec_sha256, "kind": pop.kind, "target": pop.target, "class": pr["class"],
           "prediction": {k: v for k, v in pr.items() if not k.startswith("_")}}
    if pop.expect is not None:
        out["expect"] = pop.expect
    if not simulate:
        return out
    t0 = time.time()
    b, red = pop.body, pr["_red"]
    if pop.target == "synchronization":
        out["simulation"] = (b.simulate(red, pop.body.factors[0]) if red["K_c"] is not None else
                             {"skipped": "no root of the dispersion relation inside the band: there is no onset, and "
                                         "no coupling above one to simulate"})
    elif pop.target == "collective-write":
        rows = []
        for k, row in enumerate(pr["rows"]):
            if pop.kind == "units":
                sim = b.sweep(red, row["h"], int(row["N"]), row["r"], pop.sweep_from, pop.sweep_to, replicas,
                              sample=pop.sample, seed=seed + k)
            else:
                sim = b.sweep(red, row["h"], row["N"], row["r"], pop.sweep_from, pop.sweep_to, replicas,
                              seed=seed + k, dt=pop.dt)
            agrees = abs(sim["P"] - row["P"]) < 4 * sim["stderr"] + 0.02
            rows.append({**row, "measured": sim["P"], "stderr": sim["stderr"], "agrees": bool(agrees)})
        out["simulation"] = {"rows": rows, "within_gate": int(sum(r["agrees"] for r in rows))}
    else:
        st = pr["_step"]
        rows = []
        for k, N in enumerate(pop.sizes):
            times = b.step_times(red, N, pop.step_from, pop.step_to, pop.threshold, pop.duration, replicas, st,
                                 seed=seed + k, dt=pop.dt)
            try:
                rows.append({"N": N, **passage.check(times, st["exponent"], st["d"])})
            except ValueError as error:
                rows.append({"N": N, "error": str(error)})
        good = [r for r in rows if "error" not in r]
        sim = {"rows": rows, "within_gate": int(sum(r["agrees"] for r in good))}
        if len(good) > 1:
            sim["median_slope"] = passage.median_slope([r["N"] for r in good], [r["Lambda_median"] for r in good],
                                                       [r["Lambda_median_se"] for r in good])
        out["simulation"] = sim
    out["simulation"]["seconds"] = time.time() - t0
    if pop.target != "synchronization":
        out["simulation"]["replicas"] = replicas
    return out
