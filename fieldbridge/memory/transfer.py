"""Construction steps between carriers: a map alpha and the intertwining defect it leaves.

A map alpha: Xi_s -> Xi_t transfers a memory exactly when it carries the source flow onto the target
flow, F_t(alpha(q)) = D alpha(q) F_s(q) for every q (the intertwining condition of Hyperion, Omega_t alpha =
beta Omega_s with beta = D alpha). The defect Delta(q) = F_t(alpha(q)) - D alpha(q) F_s(q) is measured on
sampled source states; a nonzero defect is the obstruction, and -Delta is the term a target material
would have to supply (the Obstruction Principle). A uniform time rescaling c (different mobilities) is
fitted separately, since it changes no state and no barrier.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Dict

import numpy as np

from .identity import Realization


@dataclass
class Transfer:
    name: str
    source: Realization
    target: Realization
    alpha: Callable[[np.ndarray], np.ndarray]
    note: str = ""

    def defect(self, rng, n: int = 256, eps: float = 1e-6) -> Dict[str, float]:
        qs = self.source.carrier.sample(rng, n)
        Fs = self.source.F(qs)
        Ft = self.target.F(self.target.carrier.wrap(self.alpha(qs)))
        DaFs = self.target.carrier.diff(self.alpha(qs - eps * Fs), self.alpha(qs + eps * Fs)) / (2 * eps)
        c = float(np.sum(Ft * DaFs) / max(np.sum(DaFs * DaFs), 1e-300))
        norm = lambda x: float(np.linalg.norm(x))
        raw = norm(Ft - DaFs) / max(norm(Ft) + norm(DaFs), 1e-300)
        res = norm(Ft - c * DaFs) / max(norm(Ft) + norm(c * DaFs), 1e-300)
        return {"relative_defect": raw, "time_rescaling": c, "relative_defect_rescaled": res,
                "exact": bool(res < 1e-6)}


def double_angle(source: Realization, target: Realization, note: str = "") -> Transfer:
    """psi = 2 theta: nematic orientations (period pi) to vectors (period 2 pi)."""
    return Transfer("psi = 2 theta", source, target, lambda q: 2.0 * q, note)


def embed_line(source: Realization, target: Realization, origin: np.ndarray, direction: np.ndarray,
               note: str = "") -> Transfer:
    """A one-dimensional memory (e.g. the pitchfork) placed along a line of the target carrier."""
    v = np.asarray(direction, float) / np.linalg.norm(direction)
    o = np.asarray(origin, float)
    return Transfer(f"line through {np.round(o, 3).tolist()} along {np.round(v, 3).tolist()}", source, target,
                    lambda q: o[None] + np.asarray(q)[:, :1] * v[None], note)
