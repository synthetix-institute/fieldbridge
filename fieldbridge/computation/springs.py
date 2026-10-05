"""A carrier for networks of nonlinear springs and point masses (Hauser et al., Biol. Cybern. 105, 355 (2011)).

A network of many springs exceeds what written-out expressions can hold, so the body is given as data:

    {"schema": "fieldbridge-springs/1", "name", "field", "question", "provenance", "assumptions",
     "network": {"positions": [[x, y], ...], "springs": [[i, j], ...], "k1": [...], "k3": [...], "d1": [...],
                 "d3": [...], "fixed": [i, ...], "inputs": {"nodes": [...], "weights": [...]}, "mass": 1.0},
     "computation": {"input": "u", "law", "amplitude", "offset", "hold", "observables": "lengths" or [spring indices],
                     "virtual_nodes", "washout", "length", "substeps"}}

Each spring s joining nodes i and j has the length l_s, the resting length l0_s (its length in the given positions,
where the network is at rest), x1 = l_s - l0_s and x2 = dl_s/dt (the exact rate of the length), and pulls its nodes
together with the tension (the source's Eq. 4)

    f_s = k3 x1^3 + k1 x1 + d3 x2^3 + d1 x2.

A free node obeys m d2p/dt2 = sum over its springs of f_s e_s + w_in u e_x (Eqs. 5, 6), with e_s the unit vector
towards the other node. The state is the positions and velocities of the free nodes; the measured signals are the
spring lengths. SpringNetwork is a spec.Body, so the simulation, the card and the estimator apply unchanged; the
linear predictions use finite differences instead of symbolic derivatives.
"""
from __future__ import annotations

from typing import Dict

import numpy as np
from scipy.optimize import root

from .spec import Body, SpecError


class SpringNetwork(Body):
    def __init__(self, spec: Dict, path: str = "<dict>"):
        self.spec, self.path = spec, path
        net = spec.get("network")
        if not isinstance(net, dict):
            raise SpecError("Missing block 'network'")
        self.p0 = np.asarray(net["positions"], float)
        self.N = len(self.p0)
        self.edges = np.asarray(net["springs"], int)
        self.L = len(self.edges)
        self.k1, self.k3, self.d1, self.d3 = (np.asarray(net[k], float) for k in ("k1", "k3", "d1", "d3"))
        self.mass = float(net.get("mass", 1.0))
        fixed = set(int(i) for i in net.get("fixed", []))
        self.free = np.array([i for i in range(self.N) if i not in fixed], int)
        self.F = len(self.free)
        self.inputs = np.asarray(net["inputs"]["nodes"], int)
        self.w_in = np.asarray(net["inputs"]["weights"], float)
        i, j = self.edges[:, 0], self.edges[:, 1]
        self.l0 = np.linalg.norm(self.p0[j] - self.p0[i], axis=1)
        self.incidence = np.zeros((self.N, self.L))
        self.incidence[i, np.arange(self.L)] = 1.0          # a spring pulls node i towards j ...
        self.incidence[j, np.arange(self.L)] = -1.0         # ... and node j towards i
        for name, arr in (("k1", self.k1), ("k3", self.k3), ("d1", self.d1), ("d3", self.d3)):
            if arr.shape != (self.L,) or not np.all(np.isfinite(arr)) or np.any(arr < 0):
                raise SpecError(f"network.{name} must give one non-negative value per spring")
        if self.edges.ndim != 2 or self.edges.shape[1] != 2 or self.edges.min() < 0 or self.edges.max() >= self.N:
            raise SpecError("network.springs must be pairs of node indices")
        if np.any(self.l0 <= 0):
            raise SpecError("two nodes of a spring coincide")
        if self.inputs.shape != self.w_in.shape:
            raise SpecError("network.inputs needs one weight per input node")
        comp = spec.get("computation")
        self._protocol(comp)
        self.input = comp.get("input", "u")
        obs = comp.get("observables", "lengths")
        self.obs_index = np.arange(self.L) if obs == "lengths" else np.asarray(obs, int)
        self.n, self.m = 4 * self.F, len(self.obs_index)
        self.obs_text = ["lengths of the springs"] if obs == "lengths" else [f"length of spring {k}" for k in self.obs_index]
        self.variables = [f"{c}{i}" for i in self.free for c in ("x", "y")] + [f"v{c}{i}" for i in self.free for c in ("x", "y")]
        self.initial = None
        if self.form != "flow":
            raise SpecError("a spring network is a flow")

    # -------------------------------------------------------------------------------------------- state and drift
    def _unpack(self, X: np.ndarray):
        S = X.shape[0]
        P = np.broadcast_to(self.p0, (S, self.N, 2)).copy()
        Vel = np.zeros((S, self.N, 2))
        P[:, self.free] = X[:, : 2 * self.F].reshape(S, self.F, 2)
        Vel[:, self.free] = X[:, 2 * self.F:].reshape(S, self.F, 2)
        return P, Vel

    def _springs(self, P: np.ndarray, Vel: np.ndarray):
        i, j = self.edges[:, 0], self.edges[:, 1]
        d = P[:, j] - P[:, i]
        l = np.linalg.norm(d, axis=-1)
        e = d / l[..., None]
        x1 = l - self.l0
        x2 = np.sum((Vel[:, j] - Vel[:, i]) * e, axis=-1)
        tension = self.k3 * x1 ** 3 + self.k1 * x1 + self.d3 * x2 ** 3 + self.d1 * x2
        return l, e, tension

    def forces(self, P: np.ndarray, Vel: np.ndarray, u) -> np.ndarray:
        _, e, tension = self._springs(P, Vel)
        Fn = np.einsum("nl,sld->snd", self.incidence, tension[..., None] * e)
        Fn[:, self.inputs, 0] += self.w_in[None, :] * np.asarray(u, float).reshape(-1, 1)
        return Fn

    def f(self, X: np.ndarray, u) -> np.ndarray:
        X = np.atleast_2d(X)
        u = np.broadcast_to(np.asarray(u, float), (X.shape[0],))
        P, Vel = self._unpack(X)
        A = self.forces(P, Vel, u) / self.mass
        S = X.shape[0]
        return np.concatenate([Vel[:, self.free].reshape(S, -1), A[:, self.free].reshape(S, -1)], axis=1)

    def observe(self, X: np.ndarray, u) -> np.ndarray:
        X = np.atleast_2d(X)
        P, Vel = self._unpack(X)
        l, _, _ = self._springs(P, Vel)
        return l[:, self.obs_index]

    def steady(self, u0: float | None = None) -> np.ndarray:
        """Rest at zero input; otherwise the static equilibrium of the free nodes under the held input u0."""
        u0 = self.offset if u0 is None else u0
        x0 = np.concatenate([self.p0[self.free].ravel(), np.zeros(2 * self.F)])
        if u0 == 0.0:
            return x0

        def residual(pos):
            X = np.concatenate([pos, np.zeros(2 * self.F)])[None, :]
            P, Vel = self._unpack(X)
            return self.forces(P, Vel, np.array([u0]))[0, self.free].ravel()

        sol = root(residual, x0[: 2 * self.F], tol=1e-12)
        if not sol.success:
            raise SpecError("no static equilibrium found under the held input")
        return np.concatenate([sol.x, np.zeros(2 * self.F)])
