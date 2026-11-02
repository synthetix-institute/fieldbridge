"""Residual-driven extension of a constructor's equation moves.

The first backend treats autonomous, deterministic polynomial vector fields on
Euclidean state spaces. It starts with supplied observables and proposed dynamics,
then solves Dh F = G(h) for missing polynomial terms. When those coordinates do
not close, it tries observables present in the residual, then source coordinates. No canonical target or response
constant is supplied. An accepted recipe remembers monomial support and coordinate
attachments, not fitted coefficients: every application derives and checks them
again. This is bounded symbolic program synthesis, not neural-network training.

Selection uses development interventions only. One selected candidate is frozen
before withheld interventions are evaluated; a failure there ends the campaign.
The SQLite library stores successful recipes and unresolved attempts separately.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
from itertools import combinations, combinations_with_replacement
import hashlib
import json
from pathlib import Path
import platform
import sqlite3
import time
import uuid

import numpy as np
import scipy
from scipy.integrate import solve_ivp
import sympy as sp

from ..verification import ConstructionError, expression
from . import spec as mspec

SCHEMA = "fieldbridge-learning/1"
ENGINE = "polynomial-intertwining/1"


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, allow_nan=False).encode()).hexdigest()


@dataclass(frozen=True)
class Policy:
    max_degree: int = 6
    max_added_coordinates: int = 2
    max_terms: int = 84
    max_candidates: int = 40
    prediction_tolerance: float = 1e-6
    convergence_tolerance: float = 1e-7
    minimum_relative_improvement: float = 0.9

    def __post_init__(self):
        if not (1 <= self.max_degree <= 8 and 0 <= self.max_added_coordinates <= 3
                and 1 <= self.max_terms <= 120 and 1 <= self.max_candidates <= 100):
            raise ValueError("Search limits exceed the bounded polynomial backend")
        if not (0 < self.convergence_tolerance <= self.prediction_tolerance <= 1e-4
                and 0.9 <= self.minimum_relative_improvement < 1):
            raise ValueError("Invalid verification tolerances")


class MoveLibrary:
    """Transactional storage of data-only recipes and their calculation records."""

    def __init__(self, path):
        self.path = Path(path)
        self.path.parent.mkdir(parents=True, exist_ok=True)
        with self.connect() as db:
            db.execute("CREATE TABLE IF NOT EXISTS attempts (id TEXT PRIMARY KEY, report TEXT NOT NULL)")
            db.execute("CREATE TABLE IF NOT EXISTS moves (id TEXT PRIMARY KEY, recipe TEXT NOT NULL, "
                       "evidence_attempt TEXT NOT NULL REFERENCES attempts(id))")
            db.execute("CREATE TABLE IF NOT EXISTS departures (id TEXT PRIMARY KEY, record TEXT NOT NULL)")
            db.execute("CREATE TABLE IF NOT EXISTS withheld_uses (model TEXT, intervention TEXT, "
                       "PRIMARY KEY(model, intervention))")
            db.execute("CREATE TABLE IF NOT EXISTS repair_campaigns (departure TEXT PRIMARY KEY, "
                       "report TEXT NOT NULL)")

    def connect(self):
        db = sqlite3.connect(self.path, timeout=30)
        db.execute("PRAGMA foreign_keys=ON")
        return db

    def recipes(self):
        with self.connect() as db:
            rows = db.execute("SELECT moves.id, recipe, report FROM moves JOIN attempts "
                              "ON evidence_attempt=attempts.id ORDER BY moves.id").fetchall()
        out = []
        for key, raw, evidence in rows:
            recipe, report = json.loads(raw), json.loads(evidence)
            if (digest(recipe) == key and recipe.get("engine") == ENGINE
                    and report.get("status") == "verified_extension"
                    and report.get("recipe_sha256") == key):
                out.append((key, recipe))
        return out

    def store(self, report, recipe=None):
        report = {**report, "attempt_id": uuid.uuid4().hex}
        if recipe is not None:
            if report.get("status") != "verified_extension":
                raise ValueError("Only verified extensions enter the move library")
            report["recipe_sha256"] = digest(recipe)
        with self.connect() as db:
            db.execute("INSERT INTO attempts VALUES (?,?)", (report["attempt_id"], json.dumps(report, allow_nan=False)))
            if recipe is not None:
                db.execute("INSERT OR IGNORE INTO moves VALUES (?,?,?)",
                           (digest(recipe), json.dumps(recipe), report["attempt_id"]))
        return report

    def retain(self, row):
        key = digest(row)
        with self.connect() as db:
            db.execute("INSERT OR IGNORE INTO departures VALUES (?,?)", (key, json.dumps(row, allow_nan=False)))
        return key

    def reserve_withheld(self, model, cases):
        """Consume tests before seeing their responses, including unsuccessful tests."""
        keys = [digest({k: r[k] for k in ("initial", "parameters")}) for r in cases]
        with self.connect() as db:
            db.execute("BEGIN IMMEDIATE")
            if any(db.execute("SELECT 1 FROM withheld_uses WHERE model=? AND intervention=?",
                              (model, k)).fetchone() for k in keys):
                return False
            db.executemany("INSERT INTO withheld_uses VALUES (?,?)", [(model, k) for k in keys])
        return True

    def departures(self):
        with self.connect() as db:
            return [(key, json.loads(raw)) for key, raw in db.execute(
                "SELECT id, record FROM departures ORDER BY id")]

    def campaign(self, key):
        with self.connect() as db:
            row = db.execute("SELECT report FROM repair_campaigns WHERE departure=?", (key,)).fetchone()
        return json.loads(row[0]) if row else None

    def store_campaign(self, key, report):
        with self.connect() as db:
            db.execute("UPDATE repair_campaigns SET report=? WHERE departure=?",
                       (json.dumps(report, allow_nan=False), key))

    def claim_campaign(self, key):
        with self.connect() as db:
            result = db.execute("INSERT OR IGNORE INTO repair_campaigns VALUES (?,?)",
                                (key, json.dumps({"status": "running"})))
            return result.rowcount == 1


def _parse_task(task):
    if task.get("schema") != SCHEMA:
        raise ConstructionError(f"Expected {SCHEMA}")
    model = task["model"]
    if (model.get("kind") != "equations" or model.get("carrier", {}).get("kind") != "euclid"
            or model.get("noise") != 0):
        raise ConstructionError("This backend requires deterministic equations on a Euclidean carrier; "
                                "stochastic, periodic and boundary-constrained models remain unsupported")
    real = mspec.load(model)
    x = real.symbolic["variables"]
    p = real.symbolic["parameters"]
    if not 1 <= len(x) <= 4:
        raise ConstructionError("Supply one to four source coordinates")
    names = {str(s): s for s in x + p}
    F = [expression(model["drift"][str(s)], names) for s in x]
    coordinates = task["coordinates"]
    if not isinstance(coordinates, dict) or not 1 <= len(coordinates) <= 4:
        raise ConstructionError("Declare one to four target coordinates")
    if any(not k.isidentifier() or k in names or k in mspec.RESERVED or k in mspec.FUNCTIONS
           or k.startswith("retained_") for k in coordinates):
        raise ConstructionError("Target names must be new, non-reserved identifiers")
    z = [sp.Symbol(k, real=True) for k in coordinates]
    h = [expression(v, names) for v in coordinates.values()]
    if set(task["candidate_drift"]) != set(coordinates):
        raise ConstructionError("Supply a proposed drift for every target coordinate")
    target_names = {str(s): s for s in z + p}
    G = [expression(task["candidate_drift"][str(s)], target_names) for s in z]
    try:
        polys = [sp.Poly(v, *x) for v in F + h] + [sp.Poly(v, *z) for v in G]
        if any(v.total_degree() > 8 or len(v.terms()) > 120 for v in polys):
            raise ConstructionError("Polynomial degree or number of terms exceeds the search budget")
    except sp.PolynomialError as error:
        raise ConstructionError("Non-polynomial dynamics require another derivation backend") from error
    if any(not item.free_symbols.intersection(x) for item in h):
        raise ConstructionError("An observable must depend on source coordinates")
    # Parameter-only poles need a domain proof not provided by this backend.
    if any(sp.denom(sp.cancel(c)).free_symbols for poly in polys for c in poly.coeffs()):
        raise ConstructionError("Parameter denominators require explicit nonzero-domain conditions")
    return real, x, p, F, z, h, G


def _validate_cases(task, real):
    groups = []
    for group in ("development", "withheld"):
        rows = task.get(group, [])
        if not 2 <= len(rows) <= 24:
            raise ConstructionError(f"{group} needs two to 24 independent interventions")
        normalized = []
        for row in rows:
            if set(row.get("parameters", {})) - set(real.params):
                raise ConstructionError("Unknown intervention parameter")
            params = {**real.params, **row.get("parameters", {})}
            if any(isinstance(v, bool) or not isinstance(v, (int, float)) for v in params.values()):
                raise ConstructionError("Intervention parameters must be finite real numbers")
            params = {k: float(v) for k, v in params.items()}
            initial = np.asarray(row["initial"], dtype=float)
            duration = float(row["duration"])
            if (initial.shape != (real.carrier.dim,) or not np.isfinite(initial).all()
                    or not all(np.isfinite(v) for v in params.values()) or not 0 < duration <= 100):
                raise ConstructionError("Invalid preparation, parameters or duration")
            for name, (lo, hi) in real.bounds.items():
                if (lo is not None and params[name] < lo) or (hi is not None and params[name] > hi):
                    raise ConstructionError("Intervention crosses a declared parameter bound")
            normalized.append({"id": str(row["id"]), "initial": initial.tolist(), "parameters": params,
                               "duration": duration})
        groups.append(normalized)
    flat = groups[0] + groups[1]
    keys = [digest({k: row[k] for k in ("initial", "parameters")}) for row in flat]
    if len(set(keys)) != len(keys) or len({r["id"] for r in flat}) != len(flat):
        raise ConstructionError("Development and withheld interventions must be distinct, without duplicate preparations")
    if real.params and not any(r["parameters"] not in [a["parameters"] for a in groups[0]] for r in groups[1]):
        raise ConstructionError("Withheld tests must include an unseen parameter setting")
    return groups


def _residual(x, F, z, h, G):
    mapped = dict(zip(z, h))
    return [sp.expand(sum(sp.diff(obs, var) * f for var, f in zip(x, F))
                      - g.subs(mapped, simultaneous=True)) for obs, g in zip(h, G)]


def _powers(n, degree):
    yield (0,) * n
    for order in range(1, degree + 1):
        for indices in combinations_with_replacement(range(n), order):
            yield tuple(indices.count(i) for i in range(n))


def _monomial(variables, powers):
    return sp.prod(v ** k for v, k in zip(variables, powers))


def _solve_support(x, p, F, z, h, baseline, powers):
    """Find corrections by coefficient equality, without using trajectory labels."""
    pulled = [sp.Poly(sp.expand(_monomial(h, k)), *x) for k in powers]
    residuals = [sp.Poly(v, *x) for v in _residual(x, F, z, h, baseline)]
    monomials = sorted(set().union(*(set(poly.monoms()) for poly in pulled + residuals)))
    A = sp.Matrix([[poly.coeff_monomial(m) for poly in pulled] for m in monomials])
    B = sp.Matrix([[poly.coeff_monomial(m) for poly in residuals] for m in monomials])
    try:
        coeffs, free = A.gauss_jordan_solve(B)
    except ValueError:
        return None
    coeffs = coeffs.subs(dict.fromkeys(free, 0))
    if any(sp.denom(sp.cancel(c)).free_symbols for c in coeffs):
        return None
    used = [i for i in range(len(powers)) if any(sp.simplify(c) != 0 for c in coeffs.row(i))]
    corrections = [sp.expand(sum(coeffs[j, i] * _monomial(z, powers[j]) for j in used))
                   for i in range(len(z))]
    result = [sp.expand(g + c) for g, c in zip(baseline, corrections)]
    return result, [powers[i] for i in used]


def _candidate(real, x, p, F, z, h, G, extra, powers):
    extra_z = [sp.Symbol(f"retained_{i}", real=True) for i in range(len(extra))]
    zz, hh = z + extra_z, h + [_monomial(x, k) for k in extra]
    baseline = G + [sp.S.Zero] * len(extra)
    solved = _solve_support(x, p, F, zz, hh, baseline, powers)
    if solved is None:
        return None
    law, support = solved
    return {"coordinates": {str(v): str(obs) for v, obs in zip(zz, hh)},
            "drift": {str(v): str(g) for v, g in zip(zz, law)},
            "recipe": {"engine": ENGINE, "source_dimension": len(x), "original_dimension": len(z),
                       "extra_powers": [list(k) for k in extra], "powers": [list(k) for k in support]}}


def verify_identity(task, candidate):
    """Reparse and recompute every component, independently of the search solution."""
    real, x, p, F, z, h, G = _parse_task(task)
    c = candidate["coordinates"]
    if list(c)[:len(z)] != [str(v) for v in z] or set(candidate["drift"]) != set(c):
        raise ConstructionError("A correction must retain every original observable and equation")
    names = {str(v): v for v in x + p}
    zz = [sp.Symbol(v, real=True) for v in c]
    if any(not v.isidentifier() or v in names for v in c):
        raise ConstructionError("Invalid target coordinate name")
    hh = [expression(value, names) for value in c.values()]
    if any(sp.simplify(a - b) != 0 for a, b in zip(h, hh)):
        raise ConstructionError("A correction cannot change the observable being tested")
    gg = [expression(candidate["drift"][str(v)], {str(s): s for s in zz + p}) for v in zz]
    residual = _residual(x, F, zz, hh, gg)
    return {"passed": all(sp.simplify(r) == 0 for r in residual), "residuals": list(map(str, residual))}


def _trajectory(rhs, initial, duration, *, fine=False, reduced=False):
    times = np.linspace(0, duration, 81)
    def bounded_rhs(t, y):
        if not np.isfinite(y).all() or np.max(np.abs(y)) > 1e8:
            raise ArithmeticError("Trajectory left the numerical range")
        value = np.asarray(rhs(t, y), dtype=float)
        if not np.isfinite(value).all():
            raise ArithmeticError("Non-finite derivative")
        return value
    sol = solve_ivp(bounded_rhs, (0, duration), initial, t_eval=times,
                    method="RK45" if reduced else "DOP853", rtol=1e-11 if fine else 1e-9,
                    atol=1e-12 if fine else 1e-11, max_step=duration / (160 if fine else 80))
    if not sol.success or sol.y.shape[1] != len(times) or not np.isfinite(sol.y).all():
        raise ArithmeticError(f"Integration did not complete: {sol.message}")
    return sol.y.T


def _calculate_cases(task, candidate, cases, policy):
    real, x, p, F, z, h, G = _parse_task(task)
    zz = [sp.Symbol(k, real=True) for k in candidate["coordinates"]]
    hh = [expression(v, {str(s): s for s in x + p}) for v in candidate["coordinates"].values()]
    gg = [expression(candidate["drift"][str(v)], {str(s): s for s in zz + p}) for v in zz]
    target_fn = sp.lambdify(zz + p, gg, "numpy")
    map_fn = sp.lambdify(x + p, hh, "numpy")
    base_fn = sp.lambdify(z + p, G, "numpy")
    rows = []
    for case in cases:
        params = case["parameters"]
        pv = [params[str(v)] for v in p]
        mapped = lambda q: np.asarray(map_fn(*q, *pv), dtype=float)
        full_rhs = lambda t, q: real.F(q, **params)
        target_rhs = lambda t, q: np.asarray(target_fn(*q, *pv), dtype=float)
        base_rhs = lambda t, q: np.asarray(base_fn(*q, *pv), dtype=float)
        T, initial = case["duration"], case["initial"]
        full = _trajectory(full_rhs, initial, T)
        finer = _trajectory(full_rhs, initial, T, fine=True)
        truth, fine_truth = np.array([mapped(q) for q in full]), np.array([mapped(q) for q in finer])
        target = _trajectory(target_rhs, mapped(initial), T, reduced=True)
        refined = _trajectory(target_rhs, mapped(initial), T, fine=True, reduced=True)
        base = _trajectory(base_rhs, mapped(initial)[:len(z)], T, fine=True, reduced=True)
        scale = 1 + np.max(np.abs(fine_truth), axis=0)
        error = float(np.max(np.abs(refined - fine_truth) / scale))
        convergence = float(max(np.max(np.abs(truth - fine_truth) / scale),
                                np.max(np.abs(target - refined) / scale)))
        omission = float(np.max(np.abs(base - fine_truth[:, :len(z)]) / scale[:len(z)]))
        rows.append({"id": case["id"], "intervention_sha256": digest(case), "error": error,
                     "solver_refinement_error": convergence, "omission_error": omission,
                     "passed": error <= policy.prediction_tolerance
                               and convergence <= policy.convergence_tolerance})
    return {"rows": rows, "passed": all(r["passed"] for r in rows),
            "max_error": max(r["error"] for r in rows),
            "max_omission_error": max(r["omission_error"] for r in rows)}


def _recipes(x, z, h, residual, policy, library):
    if library is not None:
        for key, r in library.recipes():
            extra, powers = r.get("extra_powers", []), r.get("powers", [])
            if r.get("source_dimension") != len(x) or r.get("original_dimension") != len(z):
                continue
            if (len(extra) > policy.max_added_coordinates or len({tuple(k) for k in extra}) != len(extra)
                    or any(len(k) != len(x) or any(type(i) is not int or i < 0 for i in k)
                           or not 1 <= sum(k) <= 8 for k in extra)):
                continue
            n = len(z) + len(extra)
            if (not powers or len(powers) > policy.max_terms
                    or any(len(k) != n or any(type(v) is not int or v < 0 for v in k)
                           or sum(k) > policy.max_degree for k in powers)):
                continue
            yield "learned", key, extra, powers
    eligible = []
    monomials = [k for r in residual for k in sp.Poly(r, *x).monoms()]
    monomials += [tuple(int(i == j) for j in range(len(x))) for i in range(len(x))]
    for k in monomials:
        if sum(k) and k not in eligible and all(sp.simplify(_monomial(x, k) - obs) != 0 for obs in h):
            eligible.append(k)
    eligible = eligible[:12]
    for count in range(min(policy.max_added_coordinates, len(eligible)) + 1):
        for extra in combinations(eligible, count):
            for degree in range(1, policy.max_degree + 1):
                powers = list(_powers(len(z) + count, degree))
                if len(powers) <= policy.max_terms:
                    yield "synthesized", None, extra, powers


def _valid_proposal(extra, powers, source_dim, retained_dim, policy):
    n = retained_dim + len(extra)
    return (len(extra) <= policy.max_added_coordinates
            and len({tuple(m) for m in extra}) == len(extra)
            and all(len(m) == source_dim and all(type(i) is int and i >= 0 for i in m)
                    and 1 <= sum(m) <= 8 for m in extra)
            and 0 < len(powers) <= policy.max_terms
            and len({tuple(m) for m in powers}) == len(powers)
            and all(len(m) == n and all(type(i) is int and i >= 0 for i in m)
                    and sum(m) <= policy.max_degree for m in powers))


def learn(task, library=None, policy=None, *, proposer=None):
    """Derive at most one candidate, freeze it, test it once, then optionally store it."""
    policy = policy or Policy()
    # The input hash can also identify an invalid legacy JSON containing NaN;
    # validation still rejects such numbers before any calculation.
    input_hash = hashlib.sha256(json.dumps(task, sort_keys=True).encode()).hexdigest()
    report = {"engine": ENGINE, "input_sha256": input_hash, "policy": asdict(policy),
              "novelty_established": False, "source_alignment_verified": False,
              "validation_scope": "supplied deterministic polynomial equations and declared interventions",
              "trials": [], "status": "unresolved"}
    recipe = None
    started = time.monotonic()
    try:
        real, x, p, F, z, h, G = _parse_task(task)
        development, withheld = _validate_cases(task, real)
        residual = _residual(x, F, z, h, G)
        report["initial_residual"] = list(map(str, residual))
        if all(sp.simplify(r) == 0 for r in residual):
            report["status"] = "already_closed"
        else:
            # The proposer receives equations, policy and old recipes only, never withheld outcomes.
            proposals = (proposer(task, policy, report) if proposer is not None
                         else _recipes(x, z, h, residual, policy, library))
            for origin, key, extra, powers in proposals:
                if not _valid_proposal(extra, powers, len(x), len(z), policy):
                    report.setdefault("rejected_proposals", []).append({
                        "origin": origin, "reason": "proposal_exceeds_policy_or_has_invalid_support"})
                    continue
                if len(report["trials"]) >= policy.max_candidates:
                    report["status"] = "search_budget_exhausted"
                    break
                trial = {"origin": origin, "library_move": key, "added_observables": [str(_monomial(x, k)) for k in extra],
                         "basis_size": len(powers)}
                report["trials"].append(trial)
                candidate = _candidate(real, x, p, F, z, h, G, extra, powers)
                if candidate is None:
                    trial["outcome"] = "no_polynomial_closure"
                    continue
                identity = verify_identity(task, candidate)
                if not identity["passed"]:
                    trial["outcome"] = "identity_failed"
                    continue
                checks = _calculate_cases(task, candidate, development, policy)
                trial["development"] = checks
                if not checks["passed"]:
                    report["status"] = "numerical_validation_failed"
                    trial["outcome"] = report["status"]
                    break
                if checks["max_omission_error"] <= 10 * policy.prediction_tolerance:
                    report["status"] = "insufficient_excitation"
                    trial["outcome"] = report["status"]
                    break
                trial["outcome"] = "selected"
                report.update(candidate=candidate, candidate_sha256=digest(candidate), identity=identity,
                              selection_origin=origin, reused_move=key,
                              development=checks, withheld_commitment=digest(withheld))
                # No further candidate is considered after reading withheld responses.
                model_key = digest({"variables": list(map(str, x)), "parameters": list(map(str, p)),
                                    "drift": list(map(str, F)), "carrier": "euclid"})
                if library is not None and not library.reserve_withheld(model_key, withheld):
                    report["status"] = "withheld_already_used"
                    report["reason"] = "Use fresh interventions; earlier evaluation outcomes cannot select another move"
                    break
                validation = _calculate_cases(task, candidate, withheld, policy)
                report["withheld"] = validation
                improvement = 1 - validation["max_error"] / max(validation["max_omission_error"], 1e-30)
                report["relative_improvement"] = improvement
                ok = (validation["passed"] and validation["max_omission_error"] > 10 * policy.prediction_tolerance
                      and improvement >= policy.minimum_relative_improvement)
                report["status"] = "verified_extension" if ok else "withheld_failed"
                if ok:
                    recipe = candidate["recipe"]
                break
    except (ConstructionError, mspec.SpecError, KeyError, TypeError, sp.PolynomialError) as error:
        report.update(status="unsupported_or_invalid", reason=f"{type(error).__name__}: {error}")
    except (ArithmeticError, np.linalg.LinAlgError, ValueError) as error:
        report.update(status="numerical_failure", reason=f"{type(error).__name__}: {error}")
    report["seconds"] = round(time.monotonic() - started, 4)
    report["candidate_evaluations"] = len(report["trials"])
    report["implementation_sha256"] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    files = [Path(__file__), Path(mspec.__file__), Path(mspec.__file__).with_name("identity.py"),
             Path(__file__).parents[1] / "verification.py"]
    report["implementation_files"] = {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in files}
    report["versions"] = {"python": platform.python_version(), "numpy": np.__version__,
                          "scipy": scipy.__version__, "sympy": sp.__version__}
    if library is not None:
        report = library.store(report, recipe)
    return report


def retain_record(library, record, *, source):
    """Import old target results as questions, never as successful learned moves."""
    nonfinite = []
    def preserve(value, path=""):
        if isinstance(value, float) and not np.isfinite(value):
            nonfinite.append(path)
            return {"nonfinite_value": str(value)}
        if isinstance(value, dict):
            return {k: preserve(v, f"{path}/{k}") for k, v in value.items()}
        if isinstance(value, list):
            return [preserve(v, f"{path}/{i}") for i, v in enumerate(value)]
        return value
    record = preserve(record)
    outcomes = []
    items = list(record.get("targets", {}).items())
    for section in ("completion", "write_points", "attractor"):
        value = record.get(section)
        if isinstance(value, dict):
            items.append((section, value))
            for i, point in enumerate(value.get("write_points", [])):
                if isinstance(point.get("law"), dict):
                    items.append((f"{section}/{i}/law", point["law"]))
                for j, unfolding in enumerate(point.get("unfoldings", [])):
                    if isinstance(unfolding.get("law"), dict):
                        items.append((f"{section}/{i}/unfolding/{j}", unfolding["law"]))
    for target, row in items:
        status = str(row.get("status", ""))
        detail = str(row.get("note", row.get("obstruction", "")))
        law = row.get("law_constant") or {}
        if row.get("certified") is False:
            kind = "unexplained_law_departure_requires_convergence_checks"
        elif status.startswith("reached"):
            kind = "target_reached_requires_independent_check"
        elif status.startswith("not derived") or row.get("error"):
            kind = "numerical_or_execution_failure"
        elif status.startswith("obstructed"):
            kind = "target_mismatch_requires_equation_analysis"
        else:
            kind = "unclassified"
        entry = {"source": str(source), "record_sha256": digest(record), "target": target,
                 "paper_id": record.get("paper_id"), "model": record.get("spec"), "kind": kind,
                 "original_status": status, "detail": detail, "law": law,
                 "original_result": row, "nonfinite_fields_in_record": nonfinite,
                 "learning_task": row.get("learning_task", record.get("learning_task")),
                 "eligible_for_promotion": False}
        outcomes.append({"id": library.retain(entry), "kind": kind})
    return outcomes


def cmd_learn(args):
    task = json.loads(Path(args.task).read_text(encoding="utf-8"))
    report = learn(task, MoveLibrary(args.library), Policy(max_degree=args.max_degree,
                                                          max_added_coordinates=args.max_added_coordinates))
    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    (out / "learning.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({"status": report["status"], "candidate_evaluations": report["candidate_evaluations"],
                      "selection_origin": report.get("selection_origin"), "output": str(out / "learning.json")}))
    return 0 if report["status"] in {"verified_extension", "already_closed"} else 2


def cmd_departures(args):
    from collections import Counter
    library = MoveLibrary(args.library)
    counts = Counter()
    for path in args.records:
        with Path(path).open() as stream:
            for line in stream:
                if line.strip():
                    counts.update(r["kind"] for r in retain_record(library, json.loads(line), source=path))
    print(json.dumps({"retained_outcomes": dict(counts), "library": args.library,
                      "promoted_from_records": 0}, indent=2))
    return 0


