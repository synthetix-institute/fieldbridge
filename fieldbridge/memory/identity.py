"""The identity of a memory: I_real = ((Omega, Xi); C, R, P; A).

A Realization contains the operation Omega as a drift F(q; A) on a carrier Xi, with an optional potential
(for gradient systems, F = -mobility * grad V), the parameters A through which a material implements the model,
the closure C (what is specified or discarded to close the equations: boundaries, imposed conservation, bath) as text,
the observable R through which the state is read, and the control parameter along which a protocol P can change
the landscape (the parameter a write-hot-fix-cold sweep acts on). A write is always a bounded drive toward a
target state:

  torus           energy -h cos(m (q - target)) per coordinate, m = 2 pi / period (a field aligning angles)
  euclid/orthant  a force of magnitude h pointing from q toward the target (a directed bounded drive)

Nothing here asserts how a realization writes or retains a state; the analysis measures it.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable, Dict, Optional, Tuple

import numpy as np

from .carriers import Carrier

Drift = Callable[[np.ndarray, Dict[str, float]], np.ndarray]


@dataclass
class Realization:
    name: str
    carrier: Carrier
    drift: Drift
    params: Dict[str, float]
    potential: Optional[Callable[[np.ndarray, Dict[str, float]], np.ndarray]] = None
    control: Optional[str] = None
    control_range: Optional[Tuple[float, float]] = None
    noise: float = 0.05
    dt: float = 0.01
    closure: str = ""
    observable: str = "identity of the nearest stored state"
    observable_fn: Optional[Callable[[np.ndarray], np.ndarray]] = None
    provenance: str = ""
    transport: Optional[str] = None
    tags: Tuple[str, ...] = field(default_factory=tuple)
    material: Tuple[str, ...] = field(default_factory=tuple)   # further tunable parameters (to cancel an obstruction)
    null_modes: Optional[np.ndarray] = None                     # directions of a continuous symmetry (always flat)
    canonical: Optional[Callable[[np.ndarray], np.ndarray]] = None  # removes that symmetry when states are compared
    geometry: Optional[Dict] = None                             # positions and bonds, for drawings

    # -- dynamics ------------------------------------------------------------------------------
    def F(self, q: np.ndarray, /, **over: float) -> np.ndarray:
        # q is positional-only, so that a parameter of the specification may also be named q
        p = {**self.params, **over}
        q2 = np.atleast_2d(np.asarray(q, float))
        out = self.drift(q2, p)
        return out if np.ndim(q) > 1 else out[0]

    def V(self, q: np.ndarray, /, **over: float) -> np.ndarray:
        if self.potential is None:
            raise ValueError(f"{self.name} has no potential (non-gradient dynamics)")
        p = {**self.params, **over}
        q2 = np.atleast_2d(np.asarray(q, float))
        out = self.potential(q2, p)
        return out if np.ndim(q) > 1 else out[0]

    @property
    def m(self) -> float:
        return 2.0 * np.pi / self.carrier.period

    def write_force(self, q: np.ndarray, target: np.ndarray, h: float) -> np.ndarray:
        if self.carrier.kind == "torus":
            return -h * self.m * np.sin(self.m * self.carrier.diff(target, q))
        d = self.carrier.diff(q, target)
        norm = np.linalg.norm(d, axis=-1, keepdims=True)
        return h * d / np.maximum(norm, 1e-12)

    def overlap(self, q: np.ndarray, target: np.ndarray) -> np.ndarray:
        """Pattern overlap: mean cos(m dq) on a torus, 1 - |q - t| / |t| elsewhere."""
        if self.carrier.kind == "torus":
            return np.mean(np.cos(self.m * self.carrier.diff(target, q)), axis=-1)
        t = np.asarray(target, float)
        return 1.0 - np.linalg.norm(self.carrier.diff(t, q), axis=-1) / max(np.linalg.norm(t), 1e-12)

    # -- description ---------------------------------------------------------------------------
    def slots(self) -> Dict[str, str]:
        ctrl = f"; control {self.control} in {self.control_range}" if self.control else ""
        return {
            "Omega": self.name + (f" (bond transport: {self.transport})" if self.transport else ""),
            "Xi": self.carrier.describe(),
            "C": self.closure or "thermal bath",
            "R": self.observable,
            "P": "bounded write toward a target, release, measurement" + (f"; sweep of {self.control}" if self.control else ""),
            "A": ", ".join(f"{k}={v:g}" for k, v in self.params.items()) + ctrl + f"; noise D={self.noise:g}",
        }
