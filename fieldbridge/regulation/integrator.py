"""Search for an integrator, and for conservation laws, by linear algebra on the sampled drift.

An integrator is a function phi of the state with

    dphi/dt = g(q) (h(q, u) - y0)    at every state q and input u,

so that a steady state can only have y = y0 when g keeps one sign. Two stages are tried.

Stage 1, constant gain g = k_I: phi = w.q + v.ln q (logarithms on the orthant only). The identity
w.F + v.(F/q) - k_I h + k_I y0 = 0 is linear in (w, v, k_I, k_I y0), so these are a null vector of the matrix with
rows [F, F/q, h, 1] at sampled states and inputs. A null vector without the columns h and 1 is a conservation law;
one with only the column 1 is a constant drift, which rules out a steady state.

Stage 2, gain linear in the state (orthant only, tried when stage 1 finds no integrator): y0 is the output at the
steady state and the columns are (h - y0) and q_j (h - y0); the gain must keep one sign on the sampled states. This
finds, for example, the incoherent feedforward loop dx/dt = k3 u - k4 x, dy/dt = k1 u - k2 x y, for which
phi = (k1/k3) x - y has dphi/dt = k2 x (y - y0).

Stage 3, a rate that depends on the state only through the output (tried when stages 1 and 2 find none): a linear
combination w.q whose rate w.F is a function psi(h) of the output alone, with psi(y0) = 0 and psi of one sign on
each side. The condition grad(w.F) parallel to grad h, in the coordinates (q, u), is linear in w:
    sum_k w_k (d_i F_k d_j h - d_j F_k d_i h) = 0   for every pair of coordinates i < j and every sample.
This finds, for example, the methylation level of bacterial chemotaxis, dm/dt = F(a) with F decreasing through a0.

Integrators in other coordinates are not found. A robust return without a detected integrator is reported as such.
"""
from __future__ import annotations

from typing import Callable, Dict, List, Optional, Tuple

import numpy as np

from .spec import Regulated

NULL_TOL = 1e-9        # singular value / largest singular value below which a direction is a null direction
CHECK_TOL = 1e-7       # residual of the identity on fresh samples, relative to the size of its terms


# ---------------------------------------------------------------------------------------------------- sampling
def sample_states(m: Regulated, centers: np.ndarray, count: int, rng, spread: float = 1.0) -> np.ndarray:
    """States around the given centres: within a factor 10**spread on the orthant, Gaussian on R^n."""
    centers = np.atleast_2d(np.asarray(centers, float))
    pick = centers[rng.integers(len(centers), size=count)]
    if m.orthant:
        top = np.max(np.abs(centers))
        floor = 1e-3 * top if top > 0 else 1e-3
        base = np.maximum(np.abs(pick), floor)
        return base * 10.0 ** (spread * rng.uniform(-1, 1, size=pick.shape))
    scale = np.maximum(np.max(np.abs(centers), axis=0), 1.0)
    return pick + spread * scale * rng.standard_normal(pick.shape)


def sample_inputs(m: Regulated, count: int, rng, spread: float = 1.0) -> np.ndarray:
    lo, hi = min(m.steps), max(m.steps)
    if lo > 0:
        return np.exp(rng.uniform(np.log(lo) - spread * np.log(2), np.log(hi) + spread * np.log(2), size=count))
    width = max(hi - lo, abs(hi), 1.0)
    return rng.uniform(lo - spread * width, hi + spread * width, size=count)


def _rows(m: Regulated, p: np.ndarray, Q: np.ndarray, U: np.ndarray):
    P = np.repeat(p[None, :], len(Q), axis=0)
    P[:, m.u_index] = U
    F = m.f_batch(Q, P)
    Y = m.y_batch(Q, P)
    good = np.all(np.isfinite(F), axis=1) & np.isfinite(Y)
    return F[good], Y[good], Q[good]


def _drift_columns(m: Regulated, F: np.ndarray, Q: np.ndarray, logs: bool) -> Tuple[np.ndarray, List[str]]:
    cols, names = [F], [f"w:{v}" for v in m.variables]
    if logs and m.orthant:
        cols.append(F / Q)
        names += [f"v:{v}" for v in m.variables]
    return np.hstack(cols), names


def _null(AF: np.ndarray, AG: np.ndarray, scale_F: Optional[np.ndarray] = None):
    """Directions x_F with AF x_F in the range of AG, with the matching x_G (AF x_F + AG x_G = 0).

    Columns are scaled to unit norm; AF is projected off the range of AG, so that dependences among the columns of
    AG alone are never returned. Returns (X_F, X_G, X_G in scaled columns, singular values / the largest). In scaled
    columns every null direction has unit norm, so an entry of X_G near zero means the error term carries no weight.
    scale_F, when given, sets the size of each column of AF (for a column that is itself a difference of terms)."""
    sF = np.linalg.norm(AF if scale_F is None else scale_F, axis=0)
    sF[sF == 0] = 1.0
    AFs = AF / sF
    if AG.shape[1]:
        sG = np.linalg.norm(AG, axis=0)
        sG[sG == 0] = 1.0
        AGs = AG / sG
        U, s, _ = np.linalg.svd(AGs, full_matrices=False)
        QG = U[:, : int((s > s[0] * 1e-12).sum())] if s.size and s[0] > 0 else U[:, :0]
        B = AFs - QG @ (QG.T @ AFs)
    else:
        B = AFs
    _, sv, Vt = np.linalg.svd(B, full_matrices=False)
    smax = (sv[0] if sv.size and sv[0] > 0 else 1.0) if scale_F is None else 1.0
    idx = [i for i, s in enumerate(sv) if s <= NULL_TOL * smax]
    XF = Vt[idx].T if idx else np.zeros((AF.shape[1], 0))
    if AG.shape[1]:
        XGs = -np.linalg.pinv(AGs) @ (AFs @ XF)
        XG = XGs / sG[:, None]
    else:
        XGs = XG = np.zeros((0, XF.shape[1]))
    return XF / sF[:, None], XG, XGs, sv / smax


def _split(Es: np.ndarray, tol: float = 1e-6):
    """Split null directions (columns, unit norm in scaled columns) into a basis that carries an error term and a
    basis of conservation laws. Es is X_G in scaled columns, so the tolerance is absolute."""
    k = Es.shape[1]
    if k == 0:
        return np.zeros((0, 0)), np.zeros((0, 0))
    if Es.shape[0] == 0:
        return np.zeros((k, 0)), np.eye(k)
    _, s, Vt = np.linalg.svd(Es, full_matrices=True)
    rank = int((s > tol).sum())
    return Vt[:rank].T, Vt[rank:].T


def _sparsest(x: np.ndarray, conserved: np.ndarray) -> np.ndarray:
    """x plus a combination of conserved directions (columns), chosen greedily to make x as sparse as possible:
    an integrator is defined only up to the conserved quantities, and the sparsest form is the one to read."""
    x = np.array(x, float)
    for j in range(conserved.shape[1] if conserved.size else 0):
        c = conserved[:, j]
        cands = [x] + [x - (x[i] / c[i]) * c for i in range(len(c)) if abs(c[i]) > 1e-9 * np.max(np.abs(c))]
        x = min(cands, key=lambda y: (int((np.abs(y) > 1e-9 * np.max(np.abs(y))).sum()), float(np.abs(y).sum())))
    return x


def _coefficients(x: np.ndarray, names: List[str]) -> Dict[str, Dict[str, float]]:
    out: Dict[str, Dict[str, float]] = {"w": {}, "v": {}}
    top = float(np.max(np.abs(x)))
    for val, name in zip(x, names):
        kind, var = name.split(":", 1)
        if abs(val) > 1e-9 * top:
            out[kind][var] = float(val)
    return out


def phi(coeff: Dict[str, Dict[str, float]], variables: List[str], q: np.ndarray) -> float:
    """The integrator coordinate w.q + v.ln q at state q."""
    q = np.asarray(q, float)
    total = sum(c * q[variables.index(v)] for v, c in coeff["w"].items())
    total += sum(c * np.log(q[variables.index(v)]) for v, c in coeff["v"].items())
    return float(total)


def gain_value(info: Dict[str, object], variables: List[str], q: np.ndarray) -> float:
    """The gain g(q) of an integrator of stage 1 or 2."""
    if info.get("stage") == 1:
        return float(info["k_I"])
    g = info["gain_coefficients"]
    return float(g.get("1", 0.0) + sum(c * q[variables.index(v)] for v, c in g.items() if v != "1"))


def rate(info: Dict[str, object], m: Regulated, q: np.ndarray, p: np.ndarray, y: float) -> float:
    """d phi/dt of an integrator found by any stage: g(q) (y - y0) for stages 1 and 2, w.F for stage 3."""
    if info.get("stage") == 3:
        w = np.array([info["coefficients"]["w"].get(v, 0.0) for v in m.variables])
        return float(w @ m.f(q, p))
    return gain_value(info, m.variables, q) * (y - float(info["set_point"]))


def _check(m: Regulated, p, x, rate: Callable, centers, rng, logs) -> Tuple[float, np.ndarray, np.ndarray]:
    """Largest relative residual of the identity AF x = rate(Q, Y) on fresh, wider samples."""
    count = max(200, 20 * len(x))
    Q = sample_states(m, centers, count, rng, spread=1.5)
    F, Y, Q = _rows(m, p, Q, sample_inputs(m, count, rng, spread=1.5))
    AF, _ = _drift_columns(m, F, Q, logs)
    lhs, rhs = AF @ x, rate(Q, Y)
    scale = np.abs(AF * x).sum(axis=1) + np.abs(rhs)
    if np.linalg.norm(rhs) < 1e-6 * np.linalg.norm(scale):
        return float("inf"), Q, Y           # the error term carries no weight: a conservation law, not an integrator
    return float(np.max(np.abs(lhs - rhs) / np.maximum(scale, 1e-300))), Q, Y


# ---------------------------------------------------------------------------------------------------- invariants
def conservation_laws(m: Regulated, p: np.ndarray, centers: np.ndarray, rng, count: Optional[int] = None) -> np.ndarray:
    """Linear conservation laws: orthonormal rows L with L F(q, u) = 0 at every sampled state and input."""
    count = count or max(200, 40 * m.n)
    Q = sample_states(m, centers, count, rng)
    F, _, _ = _rows(m, p, Q, sample_inputs(m, count, rng))
    XF, _, _, _ = _null(F, np.zeros((len(F), 0)))
    if XF.shape[1] == 0:
        return np.zeros((0, m.n))
    L, _ = np.linalg.qr(XF)
    return L.T


def find_integrator(m: Regulated, p: np.ndarray, centers: np.ndarray, rng, y0_steady: Optional[float] = None,
                    count: Optional[int] = None, logs: bool = True, stage1: bool = True, stage2: bool = True,
                    stage3: bool = True) -> Dict[str, object]:
    """Stage 1, then stage 2. Returns found, stage, gain, k_I or gain_coefficients, set_point, coefficients of phi
    (w for q, v for ln q, normalized to a largest coefficient of 1 and a positive gain), the residual of the identity
    on fresh samples, the number of conservation laws and whether a constant drift was found."""
    ncols = (2 if logs and m.orthant else 1) * m.n + m.n + 2
    count = count or max(300, 40 * ncols)
    Q = sample_states(m, centers, count, rng)
    F, Y, Q = _rows(m, p, Q, sample_inputs(m, count, rng))
    AF, names = _drift_columns(m, F, Q, logs)
    out: Dict[str, object] = {"found": False, "stage": None, "conservation_laws": 0, "constant_drift": False,
                              "coordinates": "linear and logarithmic" if logs and m.orthant else "linear"}

    # ---- stage 1: constant gain
    XF, XG, XGs, sv = _null(AF, np.column_stack([Y, np.ones_like(Y)]))
    carrying, conserved = _split(XGs)
    out["conservation_laws"] = int(conserved.shape[1])
    cons = XF @ conserved if conserved.size else np.zeros((AF.shape[1], 0))
    out["smallest_singular_values"] = [float(s) for s in np.sort(sv)[:3]]
    if carrying.shape[1] == 2:
        out["constant_drift"] = True
    if stage1 and carrying.shape[1] >= 1:
        j = max(range(carrying.shape[1]), key=lambda j: abs((XGs @ carrying[:, j])[0]))
        x, (a, b) = XF @ carrying[:, j], XG @ carrying[:, j]
        if abs((XGs @ carrying[:, j])[0]) > 1e-6:
            kI, y0 = -a, -b / a
            resid, _, _ = _check(m, p, x, lambda Q_, Y_: kI * (Y_ - y0), centers, rng, logs)
            out["residual"] = resid
            if resid < CHECK_TOL:
                x = _sparsest(x, cons)
                norm, sign = float(np.max(np.abs(x))), float(np.sign(kI))
                out.update(found=True, stage=1, gain="constant", k_I=abs(kI) / norm, set_point=float(y0),
                           coefficients=_coefficients(sign * x / norm, names))
                return out
        else:
            out["constant_drift"] = True
    if stage2 and m.orthant and y0_steady is not None:
        res = _stage2(m, p, centers, rng, logs, AF, names, Y, Q, y0_steady, out, cons)
        if res is not None:
            return res
    if stage3:
        res = _stage3(m, p, centers, rng, y0_steady, out)
        if res is not None:
            return res
    return out


def _stage2(m, p, centers, rng, logs, AF, names, Y, Q, y0_steady, out, cons):
    """Gain linear in the state, set point from the steady state."""
    err = Y - y0_steady
    XF, XG, XGs, sv = _null(AF, np.column_stack([err] + [Q[:, j] * err for j in range(m.n)]))
    carrying, _ = _split(XGs)
    for j in range(carrying.shape[1]):
        x, g = XF @ carrying[:, j], -(XG @ carrying[:, j])
        resid, Qc, _ = _check(m, p, x, lambda Q_, Y_, g=g: (g[0] + Q_ @ g[1:]) * (Y_ - y0_steady), centers, rng, logs)
        if resid >= CHECK_TOL:
            continue
        gv = g[0] + Qc @ g[1:]
        if not (np.all(gv > 0) or np.all(gv < 0)):
            out["stage2_note"] = "an identity was found, but its gain changes sign on the sampled states"
            continue
        xs = _sparsest(x, cons)                    # the conserved part adds nothing to the rate
        g = g * (np.max(np.abs(x)) / np.max(np.abs(xs))) if np.max(np.abs(xs)) > 0 else g
        x = xs * (np.max(np.abs(x)) / np.max(np.abs(xs))) if np.max(np.abs(xs)) > 0 else xs
        norm, sign = float(np.max(np.abs(x))), float(np.sign(gv[0]))
        gn = sign * g / norm
        out.update(found=True, stage=2, gain="linear in the state", set_point=float(y0_steady), residual=resid,
                   coefficients=_coefficients(sign * x / norm, names),
                   gain_coefficients={"1": float(gn[0]), **{v: float(gn[1 + i]) for i, v in enumerate(m.variables)
                                                           if abs(gn[1 + i]) > 1e-9 * np.max(np.abs(gn))}})
        return out
    return None


def _gradient_rows(m: Regulated, p: np.ndarray, Q: np.ndarray, U: np.ndarray):
    """Rows d_i F . d_j h - d_j F . d_i h (one column per variable k) over pairs of coordinates (q, u) and samples,
    with the matching rows |d_i F . d_j h| + |d_j F . d_i h| that give the size of each entry."""
    P = np.repeat(p[None, :], len(Q), axis=0)
    P[:, m.u_index] = U
    J, Fu, gh, hu = m.derivatives_batch(Q, P)
    D = np.concatenate([J, Fu[:, :, None]], axis=2)            # samples x k x coordinate
    dh = np.concatenate([gh, hu[:, None]], axis=1)             # samples x coordinate
    good = np.all(np.isfinite(D.reshape(len(Q), -1)), axis=1) & np.all(np.isfinite(dh), axis=1)
    D, dh = D[good], dh[good]
    rows, sizes = [], []
    for i in range(m.n + 1):
        for j in range(i + 1, m.n + 1):
            a, b = D[:, :, i] * dh[:, j:j + 1], D[:, :, j] * dh[:, i:i + 1]
            rows.append(a - b)
            sizes.append(np.abs(a) + np.abs(b))
    if not rows:
        return np.zeros((0, m.n)), np.zeros((0, m.n))
    return np.vstack(rows), np.vstack(sizes)


def _stage3(m: Regulated, p, centers, rng, y0_steady, out, count: int = 200):
    """A combination w.q whose rate is a function of the output alone, of one sign on each side of y0."""
    Q = sample_states(m, centers, count, rng)
    U = sample_inputs(m, count, rng)
    M, sizes = _gradient_rows(m, p, Q, U)
    if not len(M):
        return None
    XF, _, _, sv = _null(M, np.zeros((len(M), 0)), scale_F=sizes)
    if XF.shape[1] == 0:
        return None
    # remove conservation laws (w.F identically zero), which satisfy the condition trivially
    F, Y, Qf = _rows(m, p, Q, U)
    for j in range(XF.shape[1]):
        w = XF[:, j] / np.max(np.abs(XF[:, j]))
        Z = F @ w
        if np.linalg.norm(Z) < 1e-8 * np.linalg.norm(np.abs(F) @ np.abs(w)):
            continue                                   # a conservation law
        # fresh, wider samples: Z must be a function of Y with one zero and one sign on each side
        Qc = sample_states(m, centers, count, rng, spread=1.5)
        Fc, Yc, _ = _rows(m, p, Qc, sample_inputs(m, count, rng, spread=1.5))
        Zc = Fc @ w
        order = np.argsort(Yc)
        Ys, Zs = Yc[order], Zc[order]
        scale = np.max(np.abs(Zs))
        crossings = np.nonzero(np.sign(Zs[:-1]) * np.sign(Zs[1:]) < 0)[0]
        if len(crossings) != 1:
            continue
        i = crossings[0]
        y0 = float(Ys[i] - Zs[i] * (Ys[i + 1] - Ys[i]) / (Zs[i + 1] - Zs[i]))
        if y0_steady is not None:
            y0 = float(y0_steady)
        side = np.sign(Zs[np.abs(Ys - y0) > 1e-6 * max(abs(y0), 1e-12)]) * np.sign(Ys[np.abs(Ys - y0) > 1e-6 * max(abs(y0), 1e-12)] - y0)
        if not (np.all(side > 0) or np.all(side < 0)):
            continue
        sign = float(side[0])
        names = [f"w:{v}" for v in m.variables]
        out.update(found=True, stage=3, gain="a function of the output", set_point=y0,
                   coefficients=_coefficients(sign * w, names), residual=float(np.min(sv)),
                   note="d phi/dt = psi(y), psi(y0) = 0, psi of the sign of y - y0")
        return out
    return None
