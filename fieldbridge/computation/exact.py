"""Exact capacities of a body with one or two state variables, computed from its equations without simulation.

Driven by independent inputs, the state after each interval is a Markov chain, x_t = M(x_{t-1}, u_t), with M the
interval map (the flow over one interval with the input held, or the map itself). For a target
z = prod_j Pn_{d_j}(u_{t-k_j} / A) with largest delay K, the covariance with a measured signal h(x_t, u_t) follows by a
backward recursion of conditional expectations over the delays:

    phi_1(x) = E_u[ h(M(x, u), u) w_0(u) ],   phi_{j+1}(x) = E_u[ phi_j(M(x, u)) w_j(u) ],   c = E_pi[ phi_{K+1} ],

with w_j = Pn_d at a delay j of the target and 1 elsewhere, pi the stationary distribution of the chain, and the
expectations over u by Gauss-Legendre quadrature (uniform inputs). The state lives on a grid (linear interpolation in
one dimension, bilinear in two); pi is the fixed point of the transfer operator of the same discretization. The
capacity is c^T Sigma^{-1} c with Sigma the covariance of the signals under pi. No expansion in the amplitude and no
sampling noise: the error is that of the grid and the quadrature.

    ExactCapacities(body, N=2001 (one variable) or 161 (two), Q=24, box=None, noise=0).capacity(target)

With measurement noise of relative size `noise` the covariance gains noise^2 times its diagonal; weak directions of
the signals, which time-multiplexed measurements produce in large numbers, then stop counting.
"""
from __future__ import annotations

from typing import Dict, Iterable, Tuple

import numpy as np
from scipy.sparse import csr_matrix

from . import ipc
from .spec import Body


class ExactCapacities:
    def __init__(self, body: Body, N: int | None = None, Q: int = 24, box=None, margin: float = 0.05,
                 rng_seed: int = 0, noise: float = 0.0):
        if body.law not in ("uniform", "binary"):
            raise NotImplementedError("uniform or binary inputs only")
        if body.n not in (1, 2):
            raise NotImplementedError("one or two state variables")
        self.body, self.d = body, body.n
        self.N = N or (2001 if self.d == 1 else 161)
        if body.law == "binary":
            xq, wq, Q = np.array([-1.0, 1.0]), np.array([1.0, 1.0]), 2
        else:
            xq, wq = np.polynomial.legendre.leggauss(Q)
        self.xq, self.wq = xq, wq / 2
        self.uq = body.amplitude * xq + body.offset
        self.Q = Q
        self.box = np.array(box, float) if box is not None else self._box(margin, rng_seed)
        axes = [np.linspace(lo, hi, self.N) for lo, hi in self.box]
        mesh = np.meshgrid(*axes, indexing="ij")
        self.points = np.stack([m.ravel() for m in mesh], axis=1)           # (P, d)
        P = len(self.points)
        X = np.repeat(self.points, Q, axis=0)
        U = np.tile(self.uq, P)
        Mx, nodes = self._interval_map(X, U, record=True)                   # (P Q, d), V x (P Q, d)
        out_mask = np.any((Mx < self.box[:, 0]) | (Mx > self.box[:, 1]), axis=1).reshape(P, Q)
        Mx = np.clip(Mx, self.box[:, 0], self.box[:, 1])
        # h at the measurement times within the interval (the last at its end): (P, Q, m V)
        self.obs_mapped = np.concatenate([body.observe(Xn, U) for Xn in nodes], axis=1).reshape(P, Q, body.m * body.V)
        # interpolation weights of each mapped point onto the grid
        pos = (Mx - self.box[:, 0]) / (self.box[:, 1] - self.box[:, 0]) * (self.N - 1)
        i0 = np.clip(np.floor(pos).astype(int), 0, self.N - 2)
        fr = pos - i0
        rows, cols, vals = [], [], []
        corners = [(0,), (1,)] if self.d == 1 else [(0, 0), (0, 1), (1, 0), (1, 1)]
        for c in corners:
            idx = i0 + np.array(c)
            weight = np.prod(np.where(np.array(c) == 1, fr, 1 - fr), axis=1)
            flat = idx[:, 0] if self.d == 1 else idx[:, 0] * self.N + idx[:, 1]
            rows.append(np.arange(P * Q))
            cols.append(flat)
            vals.append(weight)
        self.W = csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(P * Q, P))
        # stationary distribution: pi_{new}(grid) = sum over sources and nodes of pi(source) wq * weights
        T = (self.W.multiply(np.tile(self.wq, P)[:, None])).tocsr()
        S = csr_matrix((np.ones(P * Q), (np.repeat(np.arange(P), Q), np.arange(P * Q))), shape=(P, P * Q))
        self.transfer = (S @ T).T.tocsr()                                   # pi_new = transfer @ pi
        pi = np.full(P, 1.0 / P)
        for self.iterations in range(1, 100001):
            new = self.transfer @ pi
            new /= new.sum()
            if np.abs(new - pi).sum() < 1e-14:
                break
            pi = new
        self.pi = new
        # the probability, under the stationary distribution, that one interval leaves the grid (clipped back)
        self.outside = float(self.pi @ (out_mask @ self.wq))
        # the measured signals are h(x_t, u_t) with x_t = M(x_{t-1}, u_t): their mean and covariance follow from the
        # mapped observations, averaged over x_{t-1} ~ pi and u_t (correct also when h depends on the input directly)
        weights = self.pi[:, None] * self.wq[None, :]                       # (P, Q)
        self.mean = np.einsum("pq,pqm->m", weights, self.obs_mapped)
        dev = self.obs_mapped - self.mean
        self.Sigma = np.einsum("pq,pqa,pqb->ab", weights, dev, dev)
        self.set_noise(noise)

    def set_noise(self, noise: float) -> None:
        """Measurement noise: independent white noise on every signal with a standard deviation `noise` times that
        of the signal, C = c^T (Sigma + noise^2 diag Sigma)^-1 c. Without noise the capacity counts every independent
        direction of the signals however weak (cut at 1e-12 of the largest eigenvalue)."""
        self.noise = float(noise)
        if noise > 0:
            self.Sigma_inv = np.linalg.inv(self.Sigma + noise ** 2 * np.diag(np.diag(self.Sigma)))
        else:
            self.Sigma_inv = np.linalg.pinv(self.Sigma, rcond=1e-12)

    def _interval_map(self, X: np.ndarray, U: np.ndarray, record: bool = False):
        """The state after one interval; with record, also the states at the V measurement times."""
        return self.body.interval(X, U, record=record)

    def _box(self, margin: float, seed: int) -> np.ndarray:
        """The range of the driven state from one long run and the steady states at the extreme inputs, widened by
        the margin."""
        b = self.body
        states = []
        for u in (-b.amplitude, b.amplitude):
            try:
                states.append(b.steady(b.offset + u))
            except (RuntimeError, ValueError):   # no steady state at the extreme input: the long run bounds the state
                pass
        # follow one long stream to bound the attractor
        rng = np.random.default_rng(seed)
        X = b.steady()[None, :]
        traj = []
        for u in b.draw(rng, 20000) + b.offset:
            X = self._interval_map(X, np.array([u]))
            traj.append(X[0])
        traj = np.array(traj + states)
        lo, hi = traj.min(axis=0), traj.max(axis=0)
        pad = margin * (hi - lo) + 1e-9
        return np.stack([lo - pad, hi + pad], axis=1)

    def covariance(self, target: Tuple[Tuple[int, int], ...]) -> np.ndarray:
        """c = E[(h - mean) z] for one target, by the backward recursion."""
        degs = dict(target)
        K = max(degs)
        w = lambda j: (ipc.polynomial(self.xq, degs[j], self.body.law) if j in degs else np.ones(self.Q)) * self.wq
        phi = np.einsum("pqm,q->pm", self.obs_mapped - self.mean, w(0))
        P = len(self.points)
        for j in range(1, K + 1):
            phi = np.einsum("pqm,q->pm", (self.W @ phi).reshape(P, self.Q, -1), w(j))
        return self.pi @ phi

    def capacity(self, target) -> float:
        c = self.covariance(target)
        return float(c @ self.Sigma_inv @ c)

    def capacities(self, delays: Dict[int, int]) -> Dict:
        by_degree, functions = {}, []
        for D, max_delay in sorted(delays.items()):
            total = 0.0
            for target in ipc.targets(D, max_delay, self.body.law):
                C = self.capacity(target)
                total += C
                functions.append((target, C))
            by_degree[D] = total
        return {"by_degree": by_degree, "total": float(sum(by_degree.values())), "functions": functions,
                "grid": self.N, "Q": self.Q, "outside": self.outside, "iterations": self.iterations}
