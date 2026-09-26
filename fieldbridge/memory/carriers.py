"""Carriers Xi: the state spaces on which a memory lives.

A carrier supplies what the dynamics and the analysis need and nothing physical: how a state is wrapped
back into the space, the difference between two states (the tangent vector), distances, and random
sampling. Three kinds cover the realizations in the library:

  euclid   R^n (displacements, conductances, order parameters)
  torus    angles modulo a period (2 pi / m for an m-fold orientation: pi for nematic rods, 2 pi for dipoles)
  orthant  non-negative concentrations R_+^n (reflected at zero)
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class Carrier:
    name: str
    dim: int
    kind: str = "euclid"
    period: float = 2.0 * np.pi
    scale: float = 1.0

    def wrap(self, q: np.ndarray) -> np.ndarray:
        if self.kind == "torus":
            return np.mod(q, self.period)
        if self.kind == "orthant":
            return np.abs(q)
        return q

    def diff(self, a: np.ndarray, b: np.ndarray) -> np.ndarray:
        """b - a as a tangent vector (shortest angular difference on a torus)."""
        d = np.asarray(b, float) - np.asarray(a, float)
        if self.kind == "torus":
            d = (d + 0.5 * self.period) % self.period - 0.5 * self.period
        return d

    def distance(self, a: np.ndarray, b: np.ndarray) -> np.ndarray:
        return np.linalg.norm(self.diff(a, b), axis=-1)

    def sample(self, rng: np.random.Generator, n: int) -> np.ndarray:
        if self.kind == "torus":
            return rng.uniform(0.0, self.period, (n, self.dim))
        if self.kind == "orthant":
            # half uniform, half log-uniform over three decades: low concentrations are reached too
            u = rng.uniform(0.0, self.scale, (n, self.dim))
            logu = self.scale * 10.0 ** rng.uniform(-3.0, 0.0, (n, self.dim))
            return np.where(rng.random((n, 1)) < 0.5, u, logu)
        return rng.uniform(-self.scale, self.scale, (n, self.dim))

    def describe(self) -> str:
        if self.kind == "torus":
            return f"{self.name}: angles on (S^1)^{self.dim}, period {self.period:.4g}"
        if self.kind == "orthant":
            return f"{self.name}: concentrations in R_+^{self.dim}"
        return f"{self.name}: R^{self.dim}"


def euclid(name: str, dim: int, scale: float = 1.0) -> Carrier:
    return Carrier(name, dim, "euclid", scale=scale)


def torus(name: str, dim: int, period: float) -> Carrier:
    return Carrier(name, dim, "torus", period=period)


def orthant(name: str, dim: int, scale: float = 1.0) -> Carrier:
    return Carrier(name, dim, "orthant", scale=scale)
