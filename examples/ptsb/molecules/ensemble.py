"""Thermal ensembles of DFT trajectories from TS1 for entries 31 and 32: which product each start reaches, and what at TS1 decides it.

    python ensemble.py init                                      # candidates, without any DFT; writes ens/init_E_S.npz (8 processes)
    python ensemble.py energy --gpu 0                            # one worker: DFT energy of each classical candidate
    python ensemble.py select                                    # chooses the classical starts from those energies; writes ens/order_E_S.npy
    python ensemble.py work --gpu 0 --set 31cl,32cl,31qc,32qc    # one worker: takes unclaimed trajectories until none is left

Two samplings of the same dividing surface (the imaginary mode of TS1 at zero, in curvilinear modes):
  cl  classical. Positions are drawn from the harmonic energy of the curvilinear modes at 298 K and then resampled by their DFT energy, so that they follow the
      Boltzmann factor on the surface. Momenta are Maxwell-Boltzmann in mass-weighted Cartesian coordinates with the flux along the normal of the surface.
      The weights have a tail index near 1, from soft rotor-like modes that a harmonic proposal covers badly, so the effective size of the raw weights does not grow with the
      number of candidates. They are truncated at sqrt(N) times their mean (Ionides, 2008). The truncated marginal of the torsion is stable under the choice of cap.
  qc  quasiclassical. Each mode is given a level of its quantum Boltzmann distribution with its zero-point energy and a random phase, the usual normal-mode sampling.

Each start runs forward until the product is decided and backward until the molecule is back on one side, so that only reactant-to-product trajectories count.
The starting conditions come only from a generator seeded by entry, sampling and index, so a stopped run keeps its order and its statistics are not selected.
"""
from __future__ import annotations

import argparse
import os
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from internal import Curvilinear  # noqa: E402
from qm import BOHR, BOND, OUT, Calc, project  # noqa: E402
from yamaki_geometries import atoms, coords  # noqa: E402

ENS = HERE / "ens"
KT = 9.44e-4                                      # 298 K in hartree
WAVENUMBER = 5140.4871                            # cm-1 per sqrt(hartree/(amu bohr^2))
HW = 4.556335e-6                                  # hartree per cm-1
TIME_UNIT_FS = 1.03275                            # sqrt(amu) bohr / sqrt(hartree) in fs
DT_FS = 0.5
PHI = {31: (6, 0, 2, 3), 32: (8, 0, 2, 3)}        # torsion through the carbon on the turning end: 180 in Z (P1), 0 in E (P2)
HOLD_FS = 40.0                                    # the torsion must stay in one well this long for the product to count as decided
T_FORWARD_FS = 300.0
T_BACKWARD_FS = 200.0
# a (75,302) grid and an SCF tolerance of 1e-8 give the forces to 5e-6 hartree/bohr (rms) of the (99,590) and 1e-10 values at half the cost
GRID, TOL = (75, 302), 1e-8


class Frame:
    """The curvilinear modes of TS1 and their conversions at any geometry."""

    def __init__(self, entry):
        self.entry = entry
        ts, irc = np.load(OUT / f"ts_{entry}.npz"), np.load(OUT / f"irc_{entry}.npz")
        self.m, self.L, self.q0 = ts["m"], ts["L"], irc["q"][0]
        self.w = np.abs(ts["w"]) / WAVENUMBER     # sqrt(hartree/(amu bohr^2))
        self.wcm = ts["w"]
        self.E0 = float(ts["E"])
        self.sm = np.sqrt(self.m)
        self.C = Curvilinear([a[0] for a in atoms(entry, "TS1")], [coords(entry, n) for n in ("TS1", "PM", "P1", "P2")], BOND[entry], self.q0, self.m, self.L)
        self.iphi = self.C.nb + self.C.na + self.C.g[2].index(PHI[entry])
        a, b = BOND[entry]
        self.ab = (a, b)
        self.R_ts = self.bond(self.q0)
        # the sign of the imaginary mode that shortens the forming bond
        self.sign = -1.0 if self.bond(self.q0 + 0.05 * self.L[:, 0]) > self.R_ts else 1.0

    def bond(self, q):
        x = (q / self.sm).reshape(-1, 3)
        return float(np.linalg.norm(x[self.ab[0]] - x[self.ab[1]]) * BOHR)

    def phi(self, q):
        return float(self.C.xi(q)[self.iphi])

    def bmatrix(self, q):
        B = np.zeros((len(self.C.xi0), len(q)))
        for i in range(len(q)):
            h = np.zeros(len(q))
            h[i] = 1e-3
            B[:, i] = self.C.wrap(self.C.xi(q + h) - self.C.xi(q - h)) / 2e-3
        return B

    def solve(self, target):
        """The geometry with curvilinear mode coordinates target, reached in stages with Newton steps, or None when there is none."""
        q = self.q0.copy()
        for scale in np.linspace(0.2, 1.0, 5):
            for _ in range(12):
                res = scale * target - self.C(q)
                if np.linalg.norm(res) < 1e-7:
                    break
                step = np.linalg.pinv(self.C.pinv @ self.bmatrix(q) @ project(q, self.m), rcond=1e-4) @ res
                n = np.linalg.norm(step)
                q = q + (step if n < 1.0 else step / n)
        return q if np.linalg.norm(target - self.C(q)) < 1e-5 else None

    def jacobian(self, q):
        """A = dQ/dq at q on the vibrational subspace; the volume factor of the map, |grad Q0| and the unit normal of the dividing surface towards products."""
        A = self.C.pinv @ self.bmatrix(q)
        P = project(q, self.m)
        A = A @ P
        s = np.linalg.svd(A, compute_uv=False)
        logj = -np.log(s[:len(self.w)]).sum()
        n = A[0] * self.sign
        return logj, np.linalg.norm(n), n / np.linalg.norm(n)


def sample_start(f, sampling, idx):
    code = {"cl": 0, "qc": 1}[sampling]
    rng = np.random.default_rng([f.entry, code, idx])
    nm = len(f.w)
    w = f.w
    for attempt in range(40):
        Q = np.zeros(nm)
        P = np.zeros(nm)
        if sampling == "cl":
            Q[1:] = rng.standard_normal(nm - 1) * np.sqrt(KT) / w[1:]
        else:
            beta = HW * np.abs(f.wcm[1:]) / KT
            levels = rng.geometric(1 - np.exp(-beta)) - 1
            energy = (levels + 0.5) * HW * np.abs(f.wcm[1:])
            angle = rng.uniform(0, 2 * np.pi, nm - 1)
            Q[1:] = np.sqrt(2 * energy) / w[1:] * np.cos(angle)
            P[1:] = np.sqrt(2 * energy) * np.sin(angle)
        q = f.solve(Q)
        if q is not None:
            break
    else:
        raise RuntimeError("no starting geometry")
    logj, gradn, n = f.jacobian(q)
    Pr = project(q, f.m)
    vharm = 0.5 * float((w[1:] ** 2 * Q[1:] ** 2).sum())
    if sampling == "cl":
        v = Pr @ (rng.standard_normal(len(q)) * np.sqrt(KT))
    else:
        v = Pr @ (np.linalg.pinv(f.bmatrix(q), rcond=1e-6) @ (f.C.BL @ P))
    n = Pr @ n
    n = n / np.linalg.norm(n)
    v = v - (v @ n) * n
    vn = np.sqrt(2 * KT * rng.exponential())
    v = v + vn * n
    return dict(q=q, v=v, Q=Q, P=P, vn=vn, logj=logj, gradn=gradn, attempts=attempt, phi=f.phi(q), n=n, vharm=vharm)


def _one(args):
    e, s, i = args
    return sample_start(_frames.setdefault(e, Frame(e)), s, i)


_frames = {}


def part_init(counts):
    """Candidates by index. A file that already holds some is extended, and the generator of each index makes the earlier ones unchanged."""
    from multiprocessing import Pool
    ENS.mkdir(exist_ok=True)
    with Pool(12) as pool:
        for e, s, count in counts:
            t0 = time.time()
            path = ENS / f"init_{e}_{s}.npz"
            old = dict(np.load(path)) if path.exists() else None
            start = len(old["q"]) if old else 0
            if start >= count:
                continue
            rows = pool.map(_one, [(e, s, i) for i in range(start, count)], chunksize=8)
            new = {k: np.array([r[k] for r in rows]) for k in rows[0]}
            if old:
                new = {k: np.concatenate([old[k], new[k]]) for k in new}
            np.savez(path, **new)
            phi = np.degrees(new["phi"])
            print(f"{e} {s}: {count} candidates; phi mean {phi.mean():.1f} sd {phi.std():.1f}, resampled {new['attempts'].sum()}, {time.time() - t0:.0f} s", flush=True)


def part_energy(gpu_label):
    from qm import Calc as _Calc
    (ENS / "energy").mkdir(parents=True, exist_ok=True)
    (ENS / "claims").mkdir(exist_ok=True)
    data = {e: np.load(ENS / f"init_{e}_cl.npz") for e in (31, 32)}
    count = max(len(data[e]["q"]) for e in data)
    for e in (31, 32):
        p = ENS / "energy" / f"TS{e}.txt"
        if claim(f"TS{e}"):
            ts = np.load(OUT / f"ts_{e}.npz")
            p.write_text(repr(_Calc(e, GRID, TOL).energy(ts["x"])))
    for i in range(count):
        for e in (31, 32):
            name = f"E{e}_{i:05d}"
            if i >= len(data[e]["q"]) or not claim(name):
                continue
            sm = np.sqrt(np.load(OUT / f"ts_{e}.npz")["m"])
            try:
                t0 = time.time()
                en = _Calc(e, GRID, TOL).energy(data[e]["q"][i] / sm)
                (ENS / "energy" / f"{name}.txt").write_text(repr(en))
                print(f"{name} {en:.6f} {time.time() - t0:.1f} s", flush=True)
            except Exception as err:
                print(f"{name} FAILED {err!r}", flush=True)


def part_select(draws):
    for e in (31, 32):
        f = Frame(e)
        d = np.load(ENS / f"init_{e}_cl.npz")
        n = len(d["q"])
        en = np.full(n, np.nan)
        for i in range(n):
            p = ENS / "energy" / f"E{e}_{i:05d}.txt"
            if p.exists():
                en[i] = float(p.read_text())
        ok = np.isfinite(en)
        e0 = float((ENS / "energy" / f"TS{e}.txt").read_text())
        lw = np.where(ok, -((en - e0) - d["vharm"]) / KT + d["logj"] + np.log(d["gradn"]), -np.inf)
        raw = np.exp(lw - lw[ok].max())
        w = np.minimum(raw, np.sqrt(ok.sum()) * raw[ok].mean())
        rng = np.random.default_rng([e, 7])
        order = rng.choice(n, size=draws, p=w / w.sum())
        np.savez(ENS / f"order_{e}_cl.npz", order=order, w=w, raw=raw, en=en)
        phi = np.degrees(d["phi"])
        ess = lambda u: u.sum() ** 2 / (u ** 2).sum()  # noqa: E731
        print(f"{e}: {ok.sum()} of {n} energies; effective size raw {ess(raw):.0f}, truncated {ess(w):.0f}; excess over harmonic mean {np.nanmean((en - e0 - d['vharm']) * 627.509):.2f} kcal/mol; "
              f"phi harmonic {phi[ok].mean():.1f} +- {phi[ok].std():.1f}; truncated {np.average(phi, weights=w):.1f} +- {np.sqrt(np.cov(phi, aweights=w)):.1f}; distinct starts {len(set(order.tolist()))} of {draws}; "
              f"share below -87: harmonic {(phi[ok] < -87).mean():.3f}, raw {raw[phi < -87].sum() / raw.sum():.3f}, truncated {w[phi < -87].sum() / w.sum():.3f}", flush=True)


def claim(name):
    try:
        os.mkdir(ENS / "claims" / name)
        return True
    except FileExistsError:
        return False


def integrate(calc, f, q, v, stop, t_end_fs, path, meta):
    """Velocity Verlet until stop(t, q) gives a verdict; written every 20 steps. A partial file is continued with a fresh density."""
    dt = DT_FS / TIME_UNIT_FS
    Q, V, E, verdict = [], [], [], None
    if path.exists():
        o = np.load(path, allow_pickle=False)
        Q, V, E = list(o["q"]), list(o["v"]), list(o["E"])
        q, v = Q[-1].copy(), V[-1].copy()
        for k, qq in enumerate(Q):
            verdict = stop(k * DT_FS, qq)
    if verdict is not None:
        return verdict, len(Q)
    e, gc = calc(q / f.sm)
    g = gc / f.sm
    if not Q:
        Q, V, E = [q.copy()], [v.copy()], [e]
    n = int(round(t_end_fs / DT_FS))
    for k in range(len(Q) - 1, n):
        v = v - 0.5 * dt * g
        q = q + dt * v
        e, gc = calc(q / f.sm)
        g = gc / f.sm
        v = v - 0.5 * dt * g
        Q.append(q.copy())
        V.append(v.copy())
        E.append(e)
        verdict = stop((k + 1) * DT_FS, q)
        if k % 20 == 19 or verdict is not None or k == n - 1:
            tmp = path.with_suffix(".tmp.npz")
            np.savez(tmp, q=np.array(Q), v=np.array(V), E=np.array(E), verdict="" if verdict is None else verdict, **meta)
            os.replace(tmp, path)
        if verdict is not None:
            return verdict, len(Q)
    return "undecided", len(Q)


def forward_stop(f):
    """Z or E when the torsion, after the bond has formed, has come within 45 degrees of planar and stayed within 73 of it for HOLD_FS. Reactant when the bond opens again."""
    state = dict(side=None, since=0.0)

    def stop(t, q):
        R = f.bond(q)
        if R > f.R_ts + 0.5:
            return "reactant"
        c = np.cos(f.phi(q))
        side = state["side"]
        if side == "Z" and c > -0.3 or side == "E" and c < 0.3 or R > 1.8:
            side = None
        if side is None and R < 1.8:
            side = "Z" if c < -0.7 else "E" if c > 0.7 else None
        if side != state["side"]:
            state.update(side=side, since=t)
        if side is not None and t - state["since"] >= HOLD_FS:
            return side
        return None

    return stop


def backward_stop(f):
    """Reactant when the forming bond has opened by 0.5 angstrom, product when it has closed below 1.75."""

    def stop(t, q):
        R = f.bond(q)
        return "reactant" if R > f.R_ts + 0.5 else "product" if R < 1.75 else None

    return stop


def tasks(sets, count):
    """All forward segments in rounds, so that every set has the same number done at any time, and then the backward segments."""
    out = []
    for direction in ("f", "b"):
        for i in range(count):
            for s in sets:
                out.append((direction, int(s[:2]), s[2:], i))
    return out


def part_work(gpu_label, sets, count):
    (ENS / "claims").mkdir(parents=True, exist_ok=True)
    (ENS / "traj").mkdir(exist_ok=True)
    frames = {e: Frame(e) for e in {int(s[:2]) for s in sets}}
    inits = {s: np.load(ENS / f"init_{s[:2]}_{s[2:]}.npz") for s in sets}
    orders = {s: (np.load(ENS / f"order_{s[:2]}_cl.npz")["order"] if s.endswith("cl") else np.arange(count)) for s in sets}
    for d, e, s, draw in tasks(sets, count):
        i = int(orders[f"{e}{s}"][draw])
        name = f"{e}{s}_{i:04d}_{d}"
        path = ENS / "traj" / f"{name}.npz"
        if (ENS / "traj" / f"{name}.done").exists() or not claim(name):
            continue
        f = frames[e]
        ini = inits[f"{e}{s}"]
        q, v = ini["q"][i].copy(), ini["v"][i].copy()
        if d == "b":
            v = -v
        t0 = time.time()
        meta = dict(entry=e, sampling=s, index=i, direction=d, dt_fs=DT_FS)
        try:
            calc = Calc(e, GRID, TOL)
            stop, t_end = (forward_stop(f), T_FORWARD_FS) if d == "f" else (backward_stop(f), T_BACKWARD_FS)
            verdict, frames_done = integrate(calc, f, q, v, stop, t_end, path, meta)
            o = np.load(path)
            np.savez(path, **{k: o[k] for k in o.files if k != "verdict"}, verdict=verdict)
            (ENS / "traj" / f"{name}.done").write_text(verdict)
            print(f"{name} {verdict} {frames_done * DT_FS:.0f} fs, {time.time() - t0:.0f} s, worker {gpu_label}", flush=True)
        except Exception as err:  # a failed SCF leaves one trajectory out, which is reported by the analysis
            (ENS / "traj" / f"{name}.failed").write_text(repr(err))
            print(f"{name} FAILED {err!r}", flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("part", choices=["init", "energy", "select", "work"])
    ap.add_argument("count", type=int, nargs="?", default=200)
    ap.add_argument("--gpu", default="0")
    ap.add_argument("--set", default="31cl,32cl,31qc,32qc")
    args = ap.parse_args()
    if args.part == "init":
        part_init([(31, "qc", 400), (32, "qc", 400), (32, "cl", 3000), (31, "cl", 10000)])
    elif args.part == "energy":
        part_energy(args.gpu)
    elif args.part == "select":
        part_select(args.count)
    else:
        part_work(args.gpu, args.set.split(","), args.count)


if __name__ == "__main__":
    main()
