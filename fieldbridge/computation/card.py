"""The computation card of a driven body: predictions from its equations beside the capacities they predict.

The card holds:
- the predictions (``predict.predict``): the class, the slowest rate of the modes that the input reaches and the
  observables see, the rank n_lin of the linear response with its singular values, the degree-1 profile at small
  amplitude, and the symmetry (odd, linear);
- for a body with one or two state variables, the exact capacities of every degree and delay (``exact``);
- on request, capacities estimated from a simulation of the full body (``ipc``), at each amplitude and at a stated
  measurement noise.

Classes (literature names): linear memory (Jaeger 2002), odd degrees only (Dambre et al. 2012), nonlinear capacity
(Dambre et al. 2012), and no fading memory (integrating; the echo state property of Jaeger 2001 fails).
"""
from __future__ import annotations

import time
from typing import Dict, Iterable, Optional

import numpy as np

from . import ipc
from .predict import predict
from .spec import Body

DELAYS = {1: 40, 2: 8, 3: 5}


def simulate(body: Body, delays: Dict[int, int], streams: int = 16, length: Optional[int] = None,
             noise: float = 0.0, seed: int = 1) -> Dict:
    """Capacities estimated from a simulation; with ``noise``, independent white measurement noise of that size
    relative to the standard deviation of each signal is added first."""
    rng = np.random.default_rng(seed)
    U, Y = body.drive(rng, streams, length)
    if noise > 0:
        sd = Y[:, body.washout:].reshape(-1, Y.shape[-1]).std(axis=0)
        Y = Y + noise * sd * rng.standard_normal(Y.shape)
    res = ipc.capacities([Y[s] for s in range(streams)], [U[s] for s in range(streams)], body.amplitude,
                         delays=delays, law=body.law, washout=body.washout)
    res.update(streams=streams, noise=noise, seed=seed)
    return res


def exact(body: Body, delays: Dict[int, int], noise: float = 0.0) -> Optional[Dict]:
    """Exact capacities for one or two state variables (None otherwise)."""
    if body.n > 2 or body.law == "gaussian":
        return None
    from .predict import linear
    if linear(body)["class"] != "fading":          # no stationary distribution without fading memory
        return None
    from .exact import ExactCapacities
    ex = ExactCapacities(body, box=body.exact_box, noise=noise)
    res = ex.capacities(delays)
    res["profile_degree_1"] = [ex.capacity(((k, 1),)) for k in range(delays.get(1, 0) + 1)] if 1 in delays else []
    return res


def card(body: Body, delays: Optional[Dict[int, int]] = None, amplitudes: Optional[Iterable[float]] = None,
         simulate_runs: bool = False, noise: float = 0.0, streams: int = 16, length: Optional[int] = None,
         seed: int = 1) -> Dict:
    delays = dict(delays or DELAYS)
    pr = predict(body)
    out = {"name": body.name, "field": body.field, "source": body.source, "path": body.path,
           "spec_sha256": body.spec_sha256, "input": body.input, "law": body.law, "amplitude": body.amplitude,
           "offset": body.offset, "hold": body.hold, "observables": body.obs_text, "virtual_nodes": body.V,
           "variables": body.n, "class": pr["structure"], "prediction": pr, "delays": {str(k): v for k, v in delays.items()},
           "noise": noise, "runs": []}
    if body.expect is not None:
        out["expect"] = body.expect
    for A in (list(amplitudes) if amplitudes else [body.amplitude]):
        b = body if A == body.amplitude else _with_amplitude(body, A)
        run: Dict = {"amplitude": A}
        ex = exact(b, delays, noise)
        if ex is not None:
            run["exact"] = {"by_degree": ex["by_degree"], "total": ex["total"], "outside": ex["outside"],
                            "profile_degree_1": ex["profile_degree_1"], "grid": ex["grid"]}
        if simulate_runs:
            t0 = time.time()
            res = simulate(b, delays, streams, length, noise, seed)
            run["simulated"] = {"by_degree": res["by_degree"], "total": res["total"], "rank": res["rank"],
                                "T": res["T"], "threshold": res["threshold"],
                                "profile_degree_1": res["profile_degree_1"],
                                "largest": [[[list(p) for p in t], c] for t, c in res["functions"][:20]],
                                "seconds": time.time() - t0}
        out["runs"].append(run)
    return out


def _with_amplitude(body: Body, A: float) -> Body:
    import copy
    b = copy.copy(body)
    b.amplitude = float(A)
    return b
