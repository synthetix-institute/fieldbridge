"""The capacity of measured signals for functions of the input history (Dambre et al., Sci. Rep. 2, 514 (2012)).

For independent inputs u_t, uniform on [-A, A], the products of normalized Legendre polynomials of delayed inputs,

    z(t) = prod_j sqrt(2 d_j + 1) P_{d_j}(u_{t - k_j} / A),

are orthonormal (Hermite polynomials He_d / sqrt(d!) for Gaussian inputs of deviation A; for binary inputs +-A with
equal probability only degree 1 per delay exists, u^2 = A^2, and the targets are products of distinct delayed inputs).
A target is written as a
tuple of (delay k_j, degree d_j) with distinct delays; its degree is the sum of the d_j. The capacity of the signals
Y (T rows after the washout, N columns) for a target z is the fraction of the variance of z that a least-squares
combination of the signals and a constant reproduces:

    C[z] = |U^T z_c|^2 / |z_c|^2,

with z_c the centred target and U an orthonormal basis of the centred signals (from the singular value decomposition,
so that dependent signals count once; N_eff is the rank). For signals and targets without serial dependence a target
unrelated to the signals has T C distributed as chi^2 with N_eff degrees of freedom; in general the null distribution
is wider (measured here: mean T C 2.2-2.4 against 2 on a nonlinear oscillator). Following Dambre et al. (SI 3.2),
capacities below `factor` times the (1 - p) quantile, divided by T, are set to zero, with factor = 2.

    capacities(Y, u, amplitude, delays={1: 60, 2: 20, 3: 10}, law="uniform", washout=0, p=1e-4, factor=2)
"""
from __future__ import annotations

from itertools import combinations
from math import factorial, sqrt
from typing import Dict, Iterator, List, Sequence, Tuple

import numpy as np
from numpy.polynomial import hermite_e, legendre
from scipy.stats import chi2

Target = Tuple[Tuple[int, int], ...]


def polynomial(x: np.ndarray, d: int, law: str) -> np.ndarray:
    """The normalized polynomial of degree d of a standardized input (uniform on [-1, 1], or unit Gaussian)."""
    c = np.zeros(d + 1)
    c[d] = 1.0
    if law == "uniform":
        return sqrt(2 * d + 1) * legendre.legval(x, c)
    if law == "gaussian":
        return hermite_e.hermeval(x, c) / sqrt(factorial(d))
    if law == "binary":
        if d != 1:
            raise ValueError("a binary input has degree 1 per delay only")
        return x
    raise ValueError(f"unknown input law {law!r}")


def partitions(D: int, largest: int | None = None) -> Iterator[Tuple[int, ...]]:
    """The partitions of D into positive parts, in non-increasing order."""
    largest = D if largest is None else largest
    if D == 0:
        yield ()
        return
    for first in range(min(D, largest), 0, -1):
        for rest in partitions(D - first, first):
            yield (first,) + rest


def targets(D: int, max_delay: int, law: str = "uniform") -> Iterator[Target]:
    """Every target of degree D with delays 0..max_delay: a partition of D placed on distinct delays (for a binary
    input, only the partition into ones)."""
    for parts in partitions(D):
        if law == "binary" and max(parts) > 1:
            continue
        n = len(parts)
        for ks in combinations(range(max_delay + 1), n):
            seen = set()
            for perm in _distinct_permutations(parts):
                key = tuple(sorted(zip(ks, perm)))
                if key not in seen:
                    seen.add(key)
                    yield key


def _distinct_permutations(parts: Tuple[int, ...]) -> Iterator[Tuple[int, ...]]:
    if len(parts) <= 1:
        yield parts
        return
    used = set()
    for i, p in enumerate(parts):
        if p in used:
            continue
        used.add(p)
        for rest in _distinct_permutations(parts[:i] + parts[i + 1:]):
            yield (p,) + rest


def capacities(Y, u, amplitude: float, delays: Dict[int, int] | None = None,
               law: str = "uniform", washout: int = 0, p: float = 1e-4, factor: float = 2.0, chunk: int = 200,
               keep_functions: bool = True) -> Dict:
    """Capacities of the signals Y (rows aligned with the inputs u) for the targets of each degree in `delays`
    (degree -> largest delay). Y and u may be lists of independent segments (streams); every segment loses its
    washout and the largest delay, and targets never reach across segments. Returns the rank, the threshold, the
    totals by degree, the profile of degree 1 by delay, and the targets with a nonzero capacity."""
    delays = {1: 60, 2: 20, 3: 10} if delays is None else delays
    segments = list(zip(Y, u)) if isinstance(Y, (list, tuple)) else [(Y, u)]
    start = washout + max(delays.values())
    Ys, xs = [], []
    for Yk, uk in segments:
        Yk, uk = np.asarray(Yk, float), np.asarray(uk, float)
        Yk = Yk[:, None] if Yk.ndim == 1 else Yk
        if len(Yk) != len(uk):
            raise ValueError("signals and inputs must have the same number of rows")
        if len(Yk) <= start:
            raise ValueError("a segment is shorter than the washout and the largest delay")
        Ys.append(Yk[start:])
        xs.append(uk / amplitude)
    Ycat = np.concatenate(Ys)
    Yc = Ycat - Ycat.mean(axis=0)
    T = len(Yc)
    Uy, s, _ = np.linalg.svd(Yc, full_matrices=False)
    rank = int(np.sum(s > s[0] * 1e-10)) if s.size and s[0] > 0 else 0
    Uy = Uy[:, :rank]
    threshold = float(factor * chi2.ppf(1 - p, max(rank, 1)) / T)
    cache: Dict[Tuple[int, int], np.ndarray] = {}

    def column(k: int, d: int) -> np.ndarray:
        if (k, d) not in cache:
            cache[k, d] = np.concatenate([polynomial(x[start - k: len(x) - k], d, law) for x in xs])
        return cache[k, d]

    by_degree: Dict[int, float] = {}
    profile: List[float] = []
    kept: List[Tuple[Target, float]] = []
    for D, max_delay in sorted(delays.items()):
        total, batch, names = 0.0, [], []

        def flush():
            nonlocal total
            if not batch:
                return
            Z = np.stack(batch, axis=1)
            Z -= Z.mean(axis=0)
            C = np.sum((Uy.T @ Z) ** 2, axis=0) / np.sum(Z * Z, axis=0)
            C = np.where(C > threshold, C, 0.0)
            for name, c in zip(names, C):
                if D == 1:
                    profile.append(float(c))
                if c > 0 and keep_functions:
                    kept.append((name, float(c)))
            total += float(C.sum())
            batch.clear()
            names.clear()

        for target in targets(D, max_delay, law):
            z = np.ones(T)
            for k, d in target:
                z = z * column(k, d)
            batch.append(z)
            names.append(target)
            if len(batch) >= chunk:
                flush()
        flush()
        by_degree[D] = total
    return {"T": T, "rank": rank, "threshold": threshold, "p": p, "factor": factor, "by_degree": by_degree,
            "total": float(sum(by_degree.values())), "profile_degree_1": profile,
            "functions": sorted(kept, key=lambda kc: -kc[1])}
