"""Memory in fields: how a conservation law and the dimension of space set the law of loss.

For a field the rates of loss form a spectrum kappa(k) over spatial modes, and the information of a write is a sum
over them: in a Gaussian (linear) field the signal-to-noise ratio of the best measurement of the field is

    SNR(t) = (1/N) sum_k |dphi_k|^2 exp(-2 kappa(k) t) / sigma_k^2 ,   I(t) = (1/2) ln(1 + SNR(t)),

so the late-time law is set by the density of slow modes near kappa = 0, not by one rate.

  conserved density (Model B: dphi/dt = M lap(dF/dphi) + conserved noise): kappa(k) = M c k^2 -> 0 at long
      wavelength. A write that changes the conserved total is kept for ever by the total, the Noether charge; its
      spatial profile spreads, and the information it carries decays as t^{-d/2}. A write that only moves density
      (a dipole) carries no charge; its information decays as t^{-d/2-1}.
  non-conserved density (Model A): kappa(k) = Gamma (kappa0 + D k^2) >= Gamma kappa0 > 0: exponential loss.

These exponents follow from the symmetry (conservation), the dimension and the shape of the write, without
solving the dynamics. The spectral sum gives them exactly for the linear field; simulate() tests them on a
nonlinear field (F = sum phi^2/2 + g phi^4/4), where no closed form is available. For the conserved nonlinear field
the late-time law is the linear one with two numbers renormalized by the static site variance <phi^2>: the noise, and
the clock (collective diffusion M T / <phi^2>, gradient dynamics with an on-site free energy). With r = T / <phi^2>,
    SNR(t) = r SNR_linear(r t),
a prediction without a fitted parameter (renormalized_snr). For the non-conserved field the quartic term only
stiffens the potential: the linear rate 2 Gamma kappa0 is a lower bound on the rate of loss (Bakry-Emery); the
self-consistent (Hartree) mass m^2 = kappa0 + 3 g <phi^2> estimates the rate itself (hartree_mass).
"""
from __future__ import annotations

import itertools
from dataclasses import dataclass
from typing import Dict, Sequence

import numpy as np


@dataclass
class FieldModel:
    d: int = 1
    L: int = 128
    conserved: bool = True
    write: str = "charge"        # "charge" (changes the conserved total) or "dipole" (moves density)
    amplitude: float = 4.0
    M: float = 1.0               # mobility (conserved) or Gamma (non-conserved)
    T: float = 0.5               # noise temperature
    g: float = 0.0               # quartic coupling (nonlinearity) for the simulation
    kappa0: float = 0.2          # mass of the non-conserved field
    Dgrad: float = 1.0           # gradient stiffness of the non-conserved field


def predicted_law(model: FieldModel) -> Dict[str, object]:
    """The law of loss from structure alone: conservation, dimension, shape of the write."""
    if not model.conserved:
        return {"law": "exponential", "rate": 2 * model.M * model.kappa0,
                "rate_is": "exact" if model.g == 0 else "lower bound",
                "reason": "no conservation: every mode relaxes at least at Gamma kappa0; a quartic term only stiffens "
                          "the potential"}
    n = 0 if model.write == "charge" else 1
    return {"law": "power", "exponent": -(model.d / 2 + n),
            "global": "the total keeps a charge-changing write for ever" if model.write == "charge" else "no charge",
            "reason": f"conserved density in d = {model.d}: slow modes kappa ~ k^2 near k = 0; the write's spectrum "
                      f"starts at order k^{2 * n} in |dphi_k|^2"}


def _lambda(d: int, L: int) -> np.ndarray:
    ks = 2 * np.pi * np.arange(L) / L
    grids = np.meshgrid(*([ks] * d), indexing="ij")
    return sum(2 * (1 - np.cos(k)) for k in grids), grids


def spectral_snr(model: FieldModel, times: Sequence[float], modes: str = "profile") -> Dict[str, object]:
    """Exact SNR of the write for the linear field on the ring of L^d sites.

    modes="profile": all modes except the conserved total, as measured in a finite closed system (the total, whose
    share A^2 / (N T) never decays, is reported separately). modes="all": every mode; for t << L^2 / M this is the
    infinite-lattice result, in which the total has no weight."""
    lam, grids = _lambda(model.d, model.L)
    N = model.L ** model.d
    if model.write == "charge":
        amp2 = np.full(lam.shape, model.amplitude ** 2)
    else:
        amp2 = model.amplitude ** 2 * 2 * (1 - np.cos(grids[0]))
    if model.conserved:
        rate = model.M * lam
        var = np.full(lam.shape, model.T)
    else:
        rate = model.M * (model.kappa0 + model.Dgrad * lam)
        var = model.T / (model.kappa0 + model.Dgrad * lam)
    mask = lam > 1e-12 if (model.conserved and modes == "profile") else np.ones(lam.shape, bool)
    out = []
    for t in times:
        out.append(float(np.sum(amp2[mask] * np.exp(-2 * rate[mask] * t) / var[mask]) / N))
    snr = np.array(out)
    total = model.amplitude ** 2 / (N * model.T) if (model.conserved and model.write == "charge") else 0.0
    return {"times": list(map(float, times)), "snr_profile": snr.tolist(),
            "information_profile": (0.5 * np.log1p(snr)).tolist(), "total_charge_snr": total}


def _lap(f: np.ndarray, d: int) -> np.ndarray:
    axes = range(1, d + 1)
    return sum(np.roll(f, 1, axis=a) + np.roll(f, -1, axis=a) - 2 * f for a in axes)


def stable_dt(model: FieldModel, phi_max: float = None) -> float:
    """Explicit step below the stability limit at the largest field value phi_max (by default the write plus three
    standard deviations of the noise)."""
    if phi_max is None:
        phi_max = model.amplitude + 3.0 * np.sqrt(model.T)
    stiff = 1.0 + 3.0 * model.g * phi_max ** 2
    if model.conserved:
        rate = model.M * 4 * model.d * stiff
    else:
        rate = model.M * (model.kappa0 + 4 * model.d * model.Dgrad + 3.0 * model.g * phi_max ** 2)
    return float(min(0.05, 0.5 / rate))


def simulate(model: FieldModel, times: Sequence[float], rng, pairs: int = 200, dt: float = None,
             t_equil: float = 20.0) -> Dict[str, object]:
    """Pairs of fields, one written and one not, evolved with the same noise. The mean difference over pairs is the
    trace of the write; its SNR (per-site noise from the stationary field) is summed over sites after removing the
    uniform (total) part. Without a fixed dt the step follows the stability limit at the current largest field value,
    so it lengthens once the write has spread."""
    d, L = model.d, model.L
    fixed = dt
    shape = (pairs,) + (L,) * d
    phi = np.sqrt(model.T) * rng.standard_normal(shape)
    if model.conserved:
        phi -= phi.mean(axis=tuple(range(1, d + 1)), keepdims=True)

    def drift(f):
        if model.conserved:
            return model.M * _lap(f + model.g * f ** 3, d)
        return -model.M * (model.kappa0 * f + model.g * f ** 3 - model.Dgrad * _lap(f, d))

    def step(f, noise):
        """Stochastic Heun step (predictor-corrector with the same noise): the stationary statistics are correct to
        second order in the step, so the step can follow the stability limit."""
        if model.conserved:
            kick = 0.0
            for a in range(1, d + 1):
                J = np.sqrt(2 * model.M * model.T * dt) * noise[a - 1]
                kick = kick + (np.roll(J, 1, axis=a) - J)
        else:
            kick = np.sqrt(2 * model.M * model.T * dt) * noise[0]
        a0 = drift(f)
        guess = f + dt * a0 + kick
        return f + 0.5 * dt * (a0 + drift(guess)) + kick

    def noise():
        return [rng.standard_normal(shape) for _ in range(d if model.conserved else 1)]

    def step_size(*fs):
        return fixed if fixed is not None else stable_dt(model, max(float(np.abs(f).max()) for f in fs))

    t = 0.0
    while t < t_equil - 1e-9:
        dt = min(step_size(phi), t_equil - t)
        phi = step(phi, noise())
        t += dt
    var_site = float(np.var(phi))
    written = phi.copy()
    origin = (slice(None),) + (0,) * d
    written[origin] += model.amplitude
    if model.write == "dipole":
        second = (slice(None), 1) + (0,) * (d - 1)
        written[second] -= model.amplitude
    t, out, steps = 0.0, [], 0
    for target in times:
        while t < target - 1e-9:
            dt = min(step_size(phi, written), target - t)
            nz = noise()
            phi, written = step(phi, nz), step(written, nz)
            t += dt
            steps += 1
        diff = (written - phi).mean(axis=0)
        if model.conserved:
            diff -= diff.mean()  # the conserved total is reported separately: it is never lost
        out.append(float(np.sum(diff ** 2) / var_site))
    return {"times": list(map(float, times)), "snr_profile": out, "site_variance": var_site,
            "steps": steps, "mean_dt": float(times[-1]) / max(steps, 1)}


def fit_exponent(times: Sequence[float], snr: Sequence[float], t_range=None) -> float:
    """Slope of ln SNR against ln t, over t_range (default: all positive points)."""
    t, s = np.asarray(times, float), np.asarray(snr, float)
    ok = s > 0
    if t_range is not None:
        ok &= (t >= t_range[0]) & (t <= t_range[1])
    return float(np.polyfit(np.log(t[ok]), np.log(s[ok]), 1)[0])


def exponent_check(model: FieldModel, L_large: int = None) -> Dict[str, object]:
    """The predicted exponent against the exact linear field in the infinite-lattice limit: every mode on a lattice
    large enough that the periodic images do not contribute, t from 20 to L^2 / 200 (diffusion length well below L)."""
    L = L_large or (4096 if model.d == 1 else 512)
    big = FieldModel(**{**model.__dict__, "L": L})
    times = np.geomspace(20.0, L ** 2 / 200.0, 12)
    exact = spectral_snr(big, times, modes="all")
    return {"L": L, "times": times.tolist(), "snr": exact["snr_profile"],
            "exponent": fit_exponent(times, exact["snr_profile"])}


def site_variance(model: FieldModel) -> float:
    """Equilibrium variance per site of the conserved field with the on-site free energy phi^2/2 + g phi^4/4 (the
    equilibrium is a product measure; the constraint on the total changes it by 1/N)."""
    x = np.linspace(-1.0, 1.0, 40001) * (12.0 * np.sqrt(model.T) + 1.0)
    w = np.exp(-(x ** 2 / 2 + model.g * x ** 4 / 4) / model.T)
    return float(np.sum(x ** 2 * w) / np.sum(w))


def renormalized_snr(model: FieldModel, times: Sequence[float]) -> Dict[str, object]:
    """Prediction for the conserved nonlinear field without a fitted parameter: the linear curve with the noise and the
    clock renormalized by the static site variance, SNR(t) = r SNR_linear(r t), r = T / <phi^2>."""
    r = model.T / site_variance(model)
    t = np.asarray(times, float)
    lin = np.asarray(spectral_snr(model, r * t)["snr_profile"])
    return {"r": float(r), "site_variance": model.T / r, "snr_profile": (r * lin).tolist()}


def hartree_mass(model: FieldModel, iterations: int = 200) -> float:
    """Self-consistent mass of the non-conserved field, m^2 = kappa0 + 3 g <phi^2>, with <phi^2> the variance of the
    Gaussian field of mass m^2 on the same lattice (an approximation; exact for g = 0)."""
    lam, _ = _lambda(model.d, model.L)
    m2 = model.kappa0
    for _ in range(iterations):
        var = float(np.mean(model.T / (m2 + model.Dgrad * lam)))
        new = model.kappa0 + 3 * model.g * var
        if abs(new - m2) < 1e-12:
            break
        m2 = 0.5 * (m2 + new)
    return float(m2)


def compare(times, sim_snr, predicted, t_min: float = 0.0) -> Dict[str, float]:
    """Root-mean-square difference of ln SNR between simulation and prediction for t >= t_min, and the ratio at the
    last time."""
    t, s, p = (np.asarray(v, float) for v in (times, sim_snr, predicted))
    ok = (t >= t_min) & (s > 0) & (p > 0)
    d = np.log(s[ok]) - np.log(p[ok])
    return {"t_min": float(t_min), "points": int(ok.sum()), "rms_log_residual": float(np.sqrt(np.mean(d ** 2))),
            "mean_log_residual": float(np.mean(d)), "last_ratio": float(s[ok][-1] / p[ok][-1])}


def _tail(times, snr, floor: float) -> np.ndarray:
    """Indices of the last five points above floor * max (empty if fewer than five)."""
    s = np.asarray(snr, float)
    idx = np.flatnonzero(s > floor * s.max())
    return idx[-5:] if idx.size >= 5 else idx[:0]


def _rate(times, snr, idx) -> float:
    """-d ln SNR / dt over the points idx (None if there are none)."""
    if len(idx) == 0:
        return None
    t, s = np.asarray(times, float), np.asarray(snr, float)
    return float(-np.polyfit(t[idx], np.log(s[idx]), 1)[0])


def run(model: FieldModel, rng, times=None, pairs: int = 120, t_equil: float = 20.0) -> Dict[str, object]:
    if times is None:
        if model.conserved:
            times = np.geomspace(2.0, 0.02 * model.L ** 2, 14)
        else:
            times = np.linspace(0.5, 10.0 / (2 * model.M * model.kappa0), 16)
    times = np.asarray(times, float)
    pred = predicted_law(model)
    exact = spectral_snr(model, times)
    sim = simulate(model, times, rng, pairs=pairs, t_equil=t_equil)
    out = {"model": model.__dict__, "prediction": pred, "exact_same_size": exact, "simulation": sim}
    if pred["law"] == "power":
        out["exponent_check"] = exponent_check(model)
        ren = renormalized_snr(model, times)
        # linear response once the written site has spread over a few sites
        out["renormalized"] = {**ren, **compare(times, sim["snr_profile"], ren["snr_profile"], t_min=10.0)}
        out["total_charge_information"] = float(0.5 * np.log1p(model.amplitude ** 2 / (model.L ** model.d * ren["site_variance"]))) \
            if model.write == "charge" else 0.0
    else:
        # rates over the window where the simulated signal is still above its noise floor, the same for every curve
        idx = _tail(times, sim["snr_profile"], 1e-6)
        m2 = hartree_mass(model)
        hartree = spectral_snr(FieldModel(**{**model.__dict__, "kappa0": m2, "g": 0.0}), times)
        out["window"] = [float(times[idx[0]]), float(times[idx[-1]])] if len(idx) else None
        out["rate_linear"] = _rate(times, exact["snr_profile"], idx)
        out["rate_simulated"] = _rate(times, sim["snr_profile"], idx)
        out["hartree"] = {"mass2": m2, "rate": _rate(times, hartree["snr_profile"], idx), "asymptotic_rate": 2 * model.M * m2}
    return out


def from_spec(spec: Dict) -> FieldModel:
    """A FieldModel from a specification of kind 'field'; entries are checked by name, type and range."""
    from .spec import SCHEMA, SpecError
    if spec.get("schema") != SCHEMA or spec.get("kind") != "field":
        raise SpecError(f"a field specification needs schema {SCHEMA!r} and kind 'field'")
    for key in ("question", "assumptions", "provenance"):
        if key not in spec:
            raise SpecError(f"Missing required field {key!r}")
    f = spec.get("field", {})
    allowed = FieldModel.__dataclass_fields__
    unknown = set(f) - set(allowed)
    if unknown:
        raise SpecError(f"unknown field entries {sorted(unknown)}; allowed: {sorted(allowed)}")
    for key, value in f.items():
        want = allowed[key].type
        if want in ("int", int) and not (isinstance(value, int) and not isinstance(value, bool)):
            raise SpecError(f"field.{key} must be an integer")
        if want in ("bool", bool) and not isinstance(value, bool):
            raise SpecError(f"field.{key} must be true or false")
        if want in ("float", float) and not (isinstance(value, (int, float)) and not isinstance(value, bool)):
            raise SpecError(f"field.{key} must be a number")
    model = FieldModel(**f)
    if model.d not in (1, 2, 3):
        raise SpecError("field.d must be 1, 2 or 3")
    if model.write not in ("charge", "dipole"):
        raise SpecError("field.write must be 'charge' (changes the conserved total) or 'dipole' (moves density)")
    if not (8 <= model.L <= 4096) or (model.d > 1 and model.L ** model.d > 512 ** 2):
        raise SpecError("field.L must lie between 8 and 4096 sites per side, at most 512^2 sites in all")
    for key in ("M", "T"):
        if getattr(model, key) <= 0:
            raise SpecError(f"field.{key} must be positive")
    if model.g < 0 or (not model.conserved and model.kappa0 <= 0):
        raise SpecError("field.g must be non-negative, and a non-conserved field needs kappa0 > 0")
    return model
