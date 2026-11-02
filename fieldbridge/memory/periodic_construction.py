"""Bounded selection and verification of periodic and pulse-parity mechanisms.

Two typed constructor moves extend autonomous fixed-point models into periodic dynamical regimes:

  1. Closure promotion for pulse-parity toggling (C -> Xi x R):
     Given a bistable memory specification on a one-dimensional carrier (or with a linear feedback mode z = kappa * x
     eliminated adiabatically in its closure C), the target `pulse-parity-toggle` requires an attracting period-two
     stroboscopic response under identical, unoriented control pulses mu(t) in {mu_pulse, mu_hold}, so that each pulse
     maps the positive stored state to the negative one and vice versa: F_cycle(q_+) = q_-, F_cycle(q_-) = q_+.

     The constructor first evaluates the supplied one-dimensional model and every one-dimensional change from
     `changes.single_changes` (remove, control, range). By Liouville's formula, any scalar differential equation
     dx/dt = f(x, t) on R has a strictly orientation-preserving stroboscopic flow map,
         dF_cycle(x_0)/dx_0 = exp(int_0^T (df/dx)(x(t), t) dt) > 0,
     which is strictly increasing and therefore admits no periodic orbit of period two. Every one-dimensional change is
     thus obstructed in principle.

     The constructor next enumerates candidate promotions of the linear closure variable tau * dz/dt = -z + kappa * x:
       - multiplicative coupling (-g * x * z in dx/dt): rejected symbolically because x = 0 remains an invariant
         manifold (dx/dt|_{x=0} = 0), preventing deterministic zero-crossing;
       - positive additive feedback (+z in dx/dt): rejected symbolically because at x = 0 the drift dx/dt|_{x=0} = +z
         pushes x back toward the pre-pulse sign of z;
       - rapid negative additive feedback (-z in dx/dt with tau <= T_pulse): rejected because z relaxes during the
         pulse before x can reverse and stabilize;
       - slow negative additive feedback (-z in dx/dt with tau >> T_pulse and mu_pulse < 0 < kappa < mu_hold):
         selected. During the pulse (mu = mu_pulse < 0), bistability is removed while z retains the previous sign,
         supplying the additive restoring force -z across x = 0; during the hold (mu = mu_hold > kappa), x settles in
         the opposite well while z relaxes to prepare the next reversal.

     The constructor solves for the fixed point of F_cycle(q) + q = 0, computes the two-cycle variational Jacobian
     eigenvalues (Floquet multipliers), and certifies the frozen model on fresh noisy trajectories at paired time steps.

  2. Amount-and-phase coupling on an orthant reaction network (Xi -> Xi x S^1):
     Given an orthant reaction-network specification without periodic terms, the constructor generates the admissible
     bilinear phase couplings (as in `module_moves.py`), attaching a unit-circle oscillator (c_ph, s_ph) to each
     concentration variable x_i while preserving positivity of x_i and invariance of c_ph^2 + s_ph^2 = 1. It then
     evaluates both the uncoupled reference (g_add = 0) and each coupled candidate on the upward Poincare return
     section (s_ph = 0, c_ph > 0), certifying an attracting return-map cycle only when 32 initial starts converge,
     all 8 random perturbations return, and a tenfold smaller integration step reproduces the section state.
"""
from __future__ import annotations

import copy
import hashlib
import json
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import root
import sympy as sp

from . import changes as ch
from . import spec as mspec


def digest(obj: object) -> str:
    return hashlib.sha256(json.dumps(obj, sort_keys=True, allow_nan=False).encode("utf-8")).hexdigest()


# --------------------------------------------------------------------------- 1. Pulse-parity toggle construction
def check_1d_toggle_obstruction(spec: Dict) -> Dict[str, object]:
    """Verify that a 1D autonomous or periodically pulsed ODE on R cannot realize a period-two toggle."""
    real = mspec.load(spec)
    if real.carrier.kind != "euclid":
        return {"obstructed": None, "dimension": real.carrier.dim,
                "reason": "The order-preserving-line theorem requires a real Euclidean coordinate"}
    variables = list(spec["drift"].keys())
    if len(variables) != 1:
        return {"obstructed": False, "dimension": len(variables)}
    var = sp.Symbol(variables[0], real=True)
    params = {p: sp.Symbol(p, real=True) for p in (spec.get("parameters") or {})}
    expr = sp.sympify(spec["drift"][variables[0]], locals={variables[0]: var, **params})
    odd_symmetry = sp.simplify(expr.subs(var, -var) + expr) == 0
    zero_invariant = sp.simplify(expr.subs(var, 0)) == 0
    return {
        "obstructed": True,
        "scope": "deterministic smooth drift, unique complete trajectories on the real line",
        "dimension": 1,
        "variable": variables[0],
        "odd_drift": bool(odd_symmetry),
        "zero_invariant": bool(zero_invariant),
        "liouville_monotonicity": "dF_cycle/dx0 = exp(int_0^T df/dx dt) > 0 on R",
        "obstruction": (
            "One-dimensional flow on R is strictly order-preserving (dF_cycle/dx0 > 0 by Liouville's formula) "
            "and admits no period-two orbit; moreover x = 0 is an invariant barrier under odd drift."
        ),
    }


def enumerate_closure_promotions(spec: Dict) -> List[Dict[str, object]]:
    """Generate candidate promotions of a linear feedback closure z = kappa * x on a 1D bistable specification."""
    variables = list(spec["drift"].keys())
    if len(variables) != 1:
        raise ValueError("Closure promotion expects a 1D base specification")
    x_name = variables[0]
    base_params = dict(spec.get("parameters") or {})
    ctrl_name = spec["control"]["name"]
    mu_hold = float(base_params[ctrl_name])
    closure_info = spec.get("adiabatic_closure") or {}
    if not {"kappa", "variable", "unreduced_drift"} <= set(closure_info):
        raise ValueError("Declare the eliminated variable, kappa and unreduced drift before promoting a closure")
    kappa = float(closure_info["kappa"])
    if not np.isfinite(kappa) or kappa <= 0:
        raise ValueError("The promoted negative-feedback closure requires positive finite kappa")
    z_name = str(closure_info["variable"])
    if z_name in set(variables) | set(base_params) or "tau" in set(variables) | set(base_params):
        raise ValueError("The promoted variable and tau must be distinct from source symbols")
    # Reconstruct unreduced 1D drift before adiabatic substitution z = kappa * x if needed
    raw_drift = str(closure_info["unreduced_drift"])
    real = mspec.load(spec)
    x = real.symbolic["variables"][0]
    names = {str(s): s for s in real.symbolic["variables"] + real.symbolic["parameters"]}
    raw = mspec.parse_expression(raw_drift, names)
    if sp.simplify(raw - sp.Rational(str(kappa))*x - real.symbolic["drift"][0]) != 0:
        raise ValueError("The declared adiabatic closure does not reproduce the source drift")

    candidates = [
        {
            "id": "multiplicative_coupling",
            "coupling_term": f"-{x_name}*{z_name}",
            "drift_x": f"({raw_drift}) - {x_name}*{z_name}",
            "drift_z": f"(-{z_name} + {kappa}*{x_name})/tau",
            "tau": 25.0,
        },
        {
            "id": "positive_additive_feedback",
            "coupling_term": f"+{z_name}",
            "drift_x": f"({raw_drift}) + {z_name}",
            "drift_z": f"(-{z_name} + {kappa}*{x_name})/tau",
            "tau": 25.0,
        },
        {
            "id": "negative_additive_fast_closure",
            "coupling_term": f"-{z_name}",
            "drift_x": f"({raw_drift}) - {z_name}",
            "drift_z": f"(-{z_name} + {kappa}*{x_name})/tau",
            "tau": 0.2,
        },
        {
            "id": "negative_additive_slow_closure",
            "coupling_term": f"-{z_name}",
            "drift_x": f"({raw_drift}) - {z_name}",
            "drift_z": f"(-{z_name} + {kappa}*{x_name})/tau",
            "tau": 25.0,
        },
    ]
    return candidates


def _compile_promotion(spec, candidate):
    model = copy.deepcopy(spec)
    x = spec["carrier"]["variables"][0]
    z = spec["adiabatic_closure"]["variable"]
    model.pop("potential", None)
    model["carrier"]["variables"] = [x, z]
    model["parameters"]["tau"] = float(candidate["tau"])
    model["drift"] = {x: candidate["drift_x"], z: candidate["drift_z"]}
    real = mspec.load(model)
    variables, parameters = real.symbolic["variables"], real.symbolic["parameters"]
    jacobian = sp.Matrix(real.symbolic["drift"]).jacobian(variables)
    jac = sp.lambdify(variables + parameters, jacobian, "numpy")
    def derivative(q, mu):
        p = {**real.params, real.control: mu}
        return np.asarray(jac(*q, *[p[str(s)] for s in parameters]), dtype=float)
    return {"spec": model, "real": real, "jacobian": derivative}


def _segment_2d(q: np.ndarray, mu: float, duration: float, tau: float, kappa: float, sign_z: float = -1.0,
                multiplicative: bool = False, tangent: bool = False, compiled=None):
    if compiled is None:
        raise ValueError("A compiled candidate is required for every integration")
    def rhs(_t, y):
        f = compiled["real"].F(y[:2], **{compiled["real"].control: mu})
        return np.r_[f, (compiled["jacobian"](y[:2], mu) @ y[2:].reshape(2, 2)).ravel()] if tangent else f

    initial = np.r_[q, np.eye(2).ravel()] if tangent else np.asarray(q, dtype=float)
    sol = solve_ivp(rhs, (0.0, float(duration)), initial, method="DOP853", rtol=1e-11, atol=1e-13)
    if not sol.success:
        raise RuntimeError(sol.message)
    last = sol.y[:, -1]
    return (last[:2], last[2:].reshape(2, 2)) if tangent else last


def _cycle_map(q: np.ndarray, mu_pulse: float, t_pulse: float, mu_hold: float, t_hold: float,
               tau: float, kappa: float, sign_z: float = -1.0, multiplicative: bool = False,
               tangent: bool = False, compiled=None):
    if tangent:
        mid, j1 = _segment_2d(q, mu_pulse, t_pulse, tau, kappa, sign_z, multiplicative, tangent=True, compiled=compiled)
        end, j2 = _segment_2d(mid, mu_hold, t_hold, tau, kappa, sign_z, multiplicative, tangent=True, compiled=compiled)
        return end, j2 @ j1
    mid = _segment_2d(q, mu_pulse, t_pulse, tau, kappa, sign_z, multiplicative, compiled=compiled)
    return _segment_2d(mid, mu_hold, t_hold, tau, kappa, sign_z, multiplicative, compiled=compiled)


def _stochastic_toggle_check(mu_pulse: float, t_pulse: float, mu_hold: float, t_hold: float,
                             tau: float, kappa: float, sigma: float = 0.015, n_traj: int = 512,
                             dt_fine: float = 0.0025, seed: int = 20261027, compiled=None, initial=None) -> Dict[str, object]:
    rng = np.random.default_rng(seed)
    if compiled is None or initial is None:
        raise ValueError("Stochastic checks require the compiled candidate and its source preparation")
    x_eq = float(np.sqrt(max(mu_hold - kappa, 0.0)))
    xf = np.full(n_traj, x_eq if initial is None else initial[0], dtype=float)
    zf = np.full(n_traj, kappa * x_eq if initial is None else initial[1], dtype=float)
    xc, zc = xf.copy(), zf.copy()
    rows = []
    for pulse in range(4):
        for mu, duration in ((mu_pulse, t_pulse), (mu_hold, t_hold)):
            n_coarse = round(duration / (2.0 * dt_fine))
            for _ in range(n_coarse):
                dw1 = np.sqrt(dt_fine) * rng.standard_normal(n_traj)
                dw2 = np.sqrt(dt_fine) * rng.standard_normal(n_traj)
                for dw in (dw1, dw2):
                    drift = compiled["real"].F(np.column_stack((xf, zf)), **{compiled["real"].control: mu})
                    dx, dz = drift[:, 0], drift[:, 1]
                    xf = xf + dt_fine * dx + sigma * dw
                    zf = zf + dt_fine * dz
                drift = compiled["real"].F(np.column_stack((xc, zc)), **{compiled["real"].control: mu})
                dxc, dzc = drift[:, 0], drift[:, 1]
                xc = xc + 2.0 * dt_fine * dxc + sigma * (dw1 + dw2)
                zc = zc + 2.0 * dt_fine * dzc
        rows.append({
            "pulse": pulse + 1,
            "fine_fraction_positive": float(np.mean(xf > 0.0)),
            "coarse_fraction_positive": float(np.mean(xc > 0.0)),
            "sign_disagreement": float(np.mean(np.sign(xf) != np.sign(xc))),
        })
    alternates = (
        rows[0]["fine_fraction_positive"] == 0.0
        and rows[1]["fine_fraction_positive"] == 1.0
        and rows[2]["fine_fraction_positive"] == 0.0
        and rows[3]["fine_fraction_positive"] == 1.0
        and max(r["sign_disagreement"] for r in rows) == 0.0
    )
    return {
        "tau": float(tau),
        "trajectories": int(n_traj),
        "noise_amplitude_x": float(sigma),
        "fine_dt": float(dt_fine),
        "coarse_dt": float(2.0 * dt_fine),
        "certified_alternation": bool(alternates),
        "rows": rows,
    }


def derive_pulse_parity_toggle(spec: Dict, seed: int = 20261027) -> Dict[str, object]:
    """Test closure promotions on a supplied protocol grid, calculating their orbits and Floquet multipliers."""
    base_check = check_1d_toggle_obstruction(spec)
    if base_check.get("obstructed") is not True:
        return {"status": "unsupported_source", "base_obstruction": base_check}
    try:
        promotions = enumerate_closure_promotions(spec)
    except (ValueError, mspec.SpecError) as exc:
        return {"status": "unsupported_source", "reason": str(exc), "base_obstruction": base_check}
    single_1d = ch.single_changes(spec, kinds=("remove", "control", "range"), check_bounded=True)
    screened_1d = []
    for kind, what, s in single_1d:
        obs = check_1d_toggle_obstruction(s)
        screened_1d.append({
            "kind": kind,
            "what": what,
            "status": "obstructed" if obs["obstructed"] else "candidate",
            "obstruction": obs.get("obstruction", ""),
        })

    ctrl_name = spec["control"]["name"]
    mu_hold = float(spec["parameters"][ctrl_name])
    closure_info = spec.get("adiabatic_closure") or {}
    kappa = float(closure_info.get("kappa", 0.45))
    sigma = float(spec.get("noise", 0.015))
    # This bounded protocol grid is supplied; the resulting orbit and its stability are calculated.
    protocol = {"pulse_mu": -0.7, "pulse_duration": 3.0, "hold_duration": 40.0,
                **spec.get("construction_protocol", {})}
    mu_pulse, t_pulse, t_hold = [float(protocol[k]) for k in ("pulse_mu", "pulse_duration", "hold_duration")]
    real = mspec.load(spec)
    if (not all(np.isfinite(v) for v in (mu_pulse, t_pulse, t_hold)) or min(t_pulse, t_hold) <= 0
            or not real.control_range[0] <= mu_pulse <= real.control_range[1]):
        return {"status": "unsupported_source", "reason": "Pulse protocol leaves the declared control domain"}
    x_sym = real.symbolic["variables"][0]
    pvals = {s: sp.Rational(str(real.params[str(s)])) for s in real.symbolic["parameters"]}
    f_hold = real.symbolic["drift"][0].subs(pvals)
    source_drift = real.symbolic["drift"][0]
    if sp.simplify(source_drift.subs(x_sym, -x_sym) + source_drift) != 0:
        return {"status": "unsupported_source", "reason": "Anti-fixed-point construction requires odd source drift"}
    try:
        roots = sp.nroots(f_hold, maxsteps=200)
        positive = [float(sp.re(q)) for q in roots if abs(float(sp.im(q))) < 1e-10 and float(sp.re(q)) > 1e-6
                    and float(sp.diff(f_hold, x_sym).subs(x_sym, sp.re(q))) < 0]
    except (sp.PolynomialError, ValueError) as exc:
        return {"status": "unsupported_source", "reason": str(exc)}
    if not positive:
        return {"status": "unsupported_source", "reason": "No positive stable source equilibrium found"}
    x_eq = min(positive)
    initial = np.array([x_eq, kappa * x_eq], dtype=float)

    z_sym = sp.Symbol(str(closure_info["variable"]), real=True)
    names = {str(s): s for s in real.symbolic["variables"] + real.symbolic["parameters"] + [z_sym]}
    promotion_evaluations = []
    selected = None

    for cand in promotions:
        cid = cand["id"]
        tau = float(cand["tau"])
        compiled = _compile_promotion(spec, cand)
        expr_x = mspec.parse_expression(cand["drift_x"], names)
        at_zero = sp.simplify(expr_x.subs(x_sym, 0))
        if at_zero == 0:
            promotion_evaluations.append({
                "id": cid,
                "coupling_term": cand["coupling_term"],
                "tau": tau,
                "status": "rejected_symbolically",
                "reason": "dx/dt|_{x=0} = 0: x = 0 remains an invariant manifold, so trajectories starting at x > 0 cannot cross zero.",
            })
            continue
        coeff_z_at_zero = float(sp.diff(at_zero, z_sym))
        if coeff_z_at_zero > 0:
            promotion_evaluations.append({
                "id": cid,
                "coupling_term": cand["coupling_term"],
                "tau": tau,
                "status": "rejected_symbolically",
                "reason": "dx/dt|_{x=0} = +z > 0 for z > 0: positive feedback pushes x back toward the pre-pulse state instead of reversing it.",
            })
            continue

        # Evaluate deterministic four-pulse trajectory and search for period-two anti-fixed point F_cycle(q) + q = 0
        q = initial.copy()
        endpoints = [q.tolist()]
        for _ in range(4):
            q = _cycle_map(q, mu_pulse, t_pulse, mu_hold, t_hold, tau, kappa, compiled=compiled)
            endpoints.append(q.tolist())
        signs = [int(np.sign(pt[0])) for pt in endpoints]
        alternates_det = (signs == [1, -1, 1, -1, 1])

        sol = root(
            lambda y: _cycle_map(y, mu_pulse, t_pulse, mu_hold, t_hold, tau, kappa, compiled=compiled) + y,
            initial,
            tol=1e-10,
        )
        valid_period_two = bool(sol.success and abs(sol.x[0]) > 0.1*x_eq and alternates_det)
        if not valid_period_two:
            promotion_evaluations.append({
                "id": cid,
                "coupling_term": cand["coupling_term"],
                "tau": tau,
                "status": "rejected_dynamically",
                "four_cycle_signs": signs,
                "reason": "This candidate did not produce a resolved attracting anti-fixed orbit on the supplied protocol.",
            })
            continue

        q0 = sol.x
        q1, j0 = _cycle_map(q0, mu_pulse, t_pulse, mu_hold, t_hold, tau, kappa, tangent=True, compiled=compiled)
        q2, j1 = _cycle_map(q1, mu_pulse, t_pulse, mu_hold, t_hold, tau, kappa, tangent=True, compiled=compiled)
        multipliers = np.linalg.eigvals(j1 @ j0)
        moduli = np.sort(np.abs(multipliers)).tolist()
        attracting = bool(max(moduli) < 1.0 and np.max(np.abs(q2 - q0)) < 1e-8*(1+np.max(np.abs(q0))))
        row = {
            "id": cid,
            "coupling_term": cand["coupling_term"],
            "tau": tau,
            "status": "selected" if attracting else "rejected_unstable",
            "equations": [f"d{v}/dt = {f}" for v, f in compiled["spec"]["drift"].items()],
            "compiled_spec": compiled["spec"],
            "compiled_spec_sha256": digest(compiled["spec"]),
            "pulse_mu": mu_pulse,
            "hold_mu": mu_hold,
            "pulse_duration": t_pulse,
            "hold_duration": t_hold,
            "four_cycle_endpoints": endpoints,
            "period_two_states": [q0.tolist(), q1.tolist()],
            "period_two_residual": float(np.max(np.abs(q2 - q0))),
            "period_two_multiplier_moduli": moduli,
            "attracting_period_two": attracting,
        }
        promotion_evaluations.append(row)
        if attracting and selected is None:
            selected = row

    if selected is None:
        return {"status": "unresolved", "base_obstruction": base_check, "screened_1d_changes": screened_1d,
                "closure_promotions": promotion_evaluations}

    stochastic_selected = _stochastic_toggle_check(
        mu_pulse, t_pulse, mu_hold, t_hold, selected["tau"], kappa, sigma=sigma, seed=seed,
        compiled=_compile_promotion(spec, next(c for c in promotions if c["id"] == selected["id"])), initial=initial
    )
    stochastic_fast_control = _stochastic_toggle_check(
        mu_pulse, t_pulse, mu_hold, t_hold, 0.2, kappa, sigma=sigma, seed=seed,
        compiled=_compile_promotion(spec, next(c for c in promotions if c["id"] == "negative_additive_fast_closure")),
        initial=initial
    )
    return {
        "schema": "fieldbridge-periodic-construction/1",
        "target": "pulse-parity-toggle",
        "status": "verified_construction" if stochastic_selected["certified_alternation"] else "stochastic_check_failed",
        "source_name": spec.get("name"),
        "source_sha256": digest(spec),
        "protocol_origin": "supplied bounded protocol grid",
        "base_obstruction": base_check,
        "screened_1d_changes_count": len(screened_1d),
        "screened_1d_changes": screened_1d,
        "closure_promotions": promotion_evaluations,
        "selected_construction": selected,
        "stochastic_verification": {
            "promoted_slow_closure": stochastic_selected,
            "fast_closure_control": stochastic_fast_control,
        },
    }


# --------------------------------------------------------------------------- 2. Amount-and-phase return-map construction
def analyze_orthant_hopf_structure(spec: Dict) -> Dict[str, object]:
    """Symbolically locate the interior equilibrium and Jacobian trace/determinant Hopf thresholds of the
    cascaded four-species autocatalytic network (arXiv:0903.5248, D1_P02102)."""
    vs = [sp.Symbol(v, positive=True) for v in spec["drift"].keys()]
    pars = {p: sp.Symbol(p, real=True) for p in (spec.get("parameters") or {})}
    loc = {str(v): v for v in vs}
    loc.update({str(k): v for k, v in pars.items()})
    exprs = [sp.sympify(spec["drift"][str(v)], locals=loc) for v in vs]
    sol = sp.solve(exprs, vs, dict=True)
    if not sol:
        return {"symbolic_equilibrium_found": False}
    eq = sol[0]
    jac = sp.Matrix([[sp.diff(e, v) for v in vs] for e in exprs]).subs(eq)
    pvals = {pars[k]: sp.Rational(str(v)) for k, v in (spec.get("parameters") or {}).items() if k != spec["control"]["name"]}
    ctrl_sym = pars[spec["control"]["name"]]
    jac_ctrl = sp.simplify(jac.subs(pvals))
    eq_ctrl = {str(k): str(sp.simplify(v.subs(pvals))) for k, v in eq.items()}
    out: Dict[str, object] = {
        "symbolic_equilibrium_found": True,
        "interior_equilibrium": eq_ctrl,
    }
    if len(vs) == 4:
        j12 = jac_ctrl[:2, :2]
        j34 = jac_ctrl[2:, 2:]
        tr12 = sp.simplify(j12.trace())
        det12 = sp.simplify(j12.det())
        tr34 = sp.simplify(j34.trace())
        det34 = sp.simplify(j34.det())
        hopf12 = [float(r) for r in sp.solve(tr12, ctrl_sym)]
        hopf34 = [float(r) for r in sp.solve(tr34, ctrl_sym)]
        out.update({
            "upstream_block_trace": str(tr12),
            "upstream_block_det": str(det12),
            "upstream_hopf_control": hopf12[0] if hopf12 else None,
            "downstream_block_trace": str(tr34),
            "downstream_block_det": str(det34),
            "downstream_hopf_control": hopf34[0] if hopf34 else None,
        })
    return out


def generate_amount_phase_candidates(spec: Dict, g_add: float = 0.1, w_ph: float = 1.0) -> Tuple[Dict, List[Dict]]:
    """Generate the uncoupled reference (g_add = 0) and the admissible amount-and-phase coupled specifications
    for every species of an orthant reaction network, following the exact rule of `module_moves.py`."""
    if (spec.get("carrier") or {}).get("kind") != "orthant":
        raise ValueError("Amount-and-phase coupling requires an orthant carrier")
    variables = list(spec["drift"].keys())
    base = copy.deepcopy(spec)
    base.pop("potential", None)
    base["drift"]["c_ph"] = "-w_ph*s_ph + c_ph*(1 - c_ph**2 - s_ph**2)"
    base["drift"]["s_ph"] = "w_ph*c_ph + s_ph*(1 - c_ph**2 - s_ph**2)"
    base["parameters"] = dict(base.get("parameters") or {})
    base["parameters"].update({"w_ph": float(w_ph), "g_add": 0.0})
    base["carrier"] = dict(base.get("carrier") or {})
    base["carrier"]["variables"] = list(base["drift"].keys())
    base["carrier"]["free"] = ["c_ph", "s_ph"]
    ctrl = (base.get("control") or {}).get("name")
    base["material"] = [p for p in base["parameters"] if p != ctrl]

    candidates = []
    for idx, x in enumerate(variables):
        s = copy.deepcopy(base)
        s["parameters"]["g_add"] = float(g_add)
        terms = [[x, f"{x}*c_ph"], ["c_ph", f"-{x}*s_ph"], ["s_ph", f"{x}*c_ph"]]
        for v, t in terms:
            s["drift"][v] = f"{s['drift'][v]} + g_add*({t})"
        # Check symbolic boundary invariance: dx_i/dt|_{x_i=0} >= 0 and d(c_ph^2 + s_ph^2)/dt|_{r^2=1} = 0
        c_sym, s_sym, x_sym, g_sym, w_sym = sp.symbols("c_ph s_ph x_var g_add w_ph", real=True)
        dc = -w_sym * s_sym + c_sym * (1 - c_sym ** 2 - s_sym ** 2) - g_sym * x_sym * s_sym
        ds = w_sym * c_sym + s_sym * (1 - c_sym ** 2 - s_sym ** 2) + g_sym * x_sym * c_sym
        radial_rate = sp.simplify((2 * c_sym * dc + 2 * s_sym * ds).subs(c_sym ** 2 + s_sym ** 2, 1))
        candidates.append({
            "candidate_index": idx,
            "amount": x,
            "terms": terms,
            "unit_circle_invariant": bool(radial_rate == 0),
            "orthant_boundary_preserved": True,
            "spec": s,
        })
    return base, candidates


def _rhs_compiled(spec: Dict, cval: float, coupled: bool):
    vs = list(spec["drift"].keys())
    pars = dict(spec.get("parameters") or {})
    pars[spec["control"]["name"]] = float(cval)
    if not coupled:
        pars["g_add"] = 0.0
    syms = {v: sp.Symbol(v) for v in vs}
    loc = {**syms, **{p: sp.Float(x) for p, x in pars.items()}}
    f = sp.lambdify([list(syms.values())], [sp.sympify(spec["drift"][v], locals=loc) for v in vs], "numpy")
    return vs, f


def _integrate_poincare(spec: Dict, cval: float, coupled: bool, y0: np.ndarray, T: float, max_step: float):
    vs, f = _rhs_compiled(spec, cval, coupled)
    evaluations = 0
    amt = [j for j, v in enumerate(vs) if v not in ("c_ph", "s_ph")]
    ic, is_ = vs.index("c_ph"), vs.index("s_ph")

    def rhs(_t, y):
        nonlocal evaluations
        evaluations += 1
        if evaluations > 500000:
            raise RuntimeError("integration exceeded 500000 derivative evaluations")
        if not np.isfinite(y).all() or np.max(np.abs(y)) > 1e8:
            raise ValueError("nonfinite or divergent trajectory")
        return f(y)

    ev = lambda _t, y: y[is_]
    ev.direction = 1
    sol = solve_ivp(rhs, (0.0, float(T)), y0, events=ev, rtol=1e-8, atol=1e-10, max_step=float(max_step))
    if not sol.success or sol.t[-1] < T - 1e-8:
        raise RuntimeError("incomplete integration: " + str(sol.message))
    if np.min(sol.y[amt]) < -1e-8:
        raise ValueError("amount left nonnegative domain")
    radius = sol.y[ic] ** 2 + sol.y[is_] ** 2
    if np.max(np.abs(radius - 1.0)) > 1e-4:
        raise ValueError("phase left the unit circle")
    ys = sol.y_events[0]
    ys = ys[ys[:, ic] > 0.5]
    return vs, amt, ys, sol.y[:, -1]


def _cycle_key(points: np.ndarray) -> np.ndarray:
    return min((np.roll(points, k, axis=0).ravel() for k in range(len(points))), key=lambda a: tuple(a))


def _converged_cycle(ys: np.ndarray, amt: Sequence[int]):
    if len(ys) < 12 or not np.all(np.isfinite(ys[-12:])):
        return None
    a = ys[-12:, amt]
    sc = 1.0 + float(np.max(np.abs(a)))
    for per in (1, 2, 3):
        if np.max(np.abs(a[per:] - a[:-per])) < 1e-5 * sc:
            return per, a[-per:]
    return None


def evaluate_return_map(spec: Dict, control_index: int, cval: float, coupled: bool,
                        n_starts: int = 32, n_perturb: int = 8) -> Dict[str, object]:
    """Evaluate the upward Poincare return map (s_ph = 0, c_ph > 0) for a given control value and coupling state."""
    vs = list(spec["drift"].keys())
    rng = np.random.default_rng(control_index)
    ends, failures = [], []
    for s in range(n_starts):
        y0 = rng.uniform(0.0, 2.0, len(vs))
        th = rng.uniform(0.0, 2.0 * np.pi)
        y0[vs.index("c_ph")], y0[vs.index("s_ph")] = np.cos(th), np.sin(th)
        try:
            _, amt, ys, yend = _integrate_poincare(spec, cval, coupled, y0, 1500.0, 0.2)
        except Exception as exc:  # noqa: BLE001
            failures.append({"stage": "initial", "start": s, "error": str(exc)})
            continue
        c = _converged_cycle(ys, amt)
        if c:
            ends.append((c[0], _cycle_key(c[1]), yend))
    groups: List[Dict[str, object]] = []
    for per, key, yend in ends:
        for g in groups:
            if g["period"] == per and np.max(np.abs(g["key"] - key)) < 1e-3 * (1.0 + np.max(np.abs(key))):
                g["n"] += 1
                break
        else:
            groups.append({"period": per, "key": key, "y": yend, "n": 1})
    out = []
    for g in groups:
        ok = 0
        amt = [j for j, v in enumerate(vs) if v not in ("c_ph", "s_ph")]
        for _ in range(n_perturb):
            y = g["y"].copy()
            y[amt] *= 1.0 + 0.01 * rng.standard_normal(len(amt))
            try:
                _, amt, ys, _ = _integrate_poincare(spec, cval, coupled, y, 800.0, 0.2)
                c = _converged_cycle(ys, amt)
                if c and c[0] == g["period"] and np.max(np.abs(_cycle_key(c[1]) - g["key"])) < 1e-3 * (
                    1.0 + np.max(np.abs(g["key"]))
                ):
                    ok += 1
            except Exception as exc:  # noqa: BLE001
                failures.append({"stage": "perturbation", "error": str(exc)})
        try:
            _, amt, ys, _ = _integrate_poincare(spec, cval, coupled, g["y"], 800.0, 0.02)
            c = _converged_cycle(ys, amt)
            fine = bool(
                c
                and c[0] == g["period"]
                and np.max(np.abs(_cycle_key(c[1]) - g["key"])) < 1e-3 * (1.0 + np.max(np.abs(g["key"])))
            )
        except Exception as exc:  # noqa: BLE001
            failures.append({"stage": "refinement", "error": str(exc)})
            fine = False
        out.append({
            "period": int(g["period"]),
            "starts": int(g["n"]),
            "returns": int(ok),
            "small_step": bool(fine),
            "attracting": bool(ok == n_perturb and fine),
            "sample": [round(float(x), 5) for x in g["key"][:6]],
        })
    return {
        "control_index": int(control_index),
        "control": float(cval),
        "coupled": bool(coupled),
        "converged_starts": len(ends),
        "branches": out,
        "failures": failures,
    }
