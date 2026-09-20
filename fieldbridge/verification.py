"""Exact calculations for supplied models, separate from analogue retrieval.

The Ito calculation proves an identity of local differential expressions.
The finite-dimensional quantum calculation proves closure under a supplied
Hamiltonian. Neither result establishes corpus provenance or physical novelty.
Expressions use a restricted arithmetic parser, never eval or sympify on input.
"""

from __future__ import annotations

import ast
import hashlib
import json
from pathlib import Path

import sympy as sp


class ConstructionError(ValueError):
    """The supplied model is outside the supported calculation contract."""


def expression(value, symbols):
    text = str(value)
    if len(text) > 4000:
        raise ConstructionError("Expression is too long")
    try:
        tree = ast.parse(text, mode="eval")
    except SyntaxError as error:
        raise ConstructionError("Expected an arithmetic expression, not an equation") from error
    if sum(1 for _ in ast.walk(tree)) > 300:
        raise ConstructionError("Expression is too complex")
    functions = {"sqrt": sp.sqrt, "sin": sp.sin, "cos": sp.cos, "exp": sp.exp, "log": sp.log}

    def visit(node):
        if isinstance(node, ast.Constant) and type(node.value) in (int, float):
            return sp.Rational(str(node.value))
        if isinstance(node, ast.Name) and node.id in symbols:
            return symbols[node.id]
        if isinstance(node, ast.UnaryOp) and isinstance(node.op, (ast.USub, ast.UAdd)):
            item = visit(node.operand)
            return -item if isinstance(node.op, ast.USub) else item
        if isinstance(node, ast.BinOp):
            left, right = visit(node.left), visit(node.right)
            if isinstance(node.op, ast.Add):
                return left + right
            if isinstance(node.op, ast.Sub):
                return left - right
            if isinstance(node.op, ast.Mult):
                return left * right
            if isinstance(node.op, ast.Div):
                return left / right
            if isinstance(node.op, ast.Pow) and right.is_Rational and abs(right) <= 16:
                return left ** right
        if (isinstance(node, ast.Call) and isinstance(node.func, ast.Name)
                and node.func.id in functions and len(node.args) == 1 and not node.keywords):
            return functions[node.func.id](visit(node.args[0]))
        raise ConstructionError("Unsupported expression syntax or undeclared symbol")

    result = visit(tree.body)
    if result.has(sp.zoo, sp.oo, -sp.oo, sp.nan):
        raise ConstructionError("Expression is not finite")
    return result


def _symbols(spec):
    parameters = spec.get("parameters", {})
    if not isinstance(parameters, dict):
        raise ConstructionError("parameters must map names to real or positive")
    names = {}
    for name, assumption in parameters.items():
        if not name.isidentifier() or name in {"x", "y", "t", "I"}:
            raise ConstructionError("Invalid or reserved parameter name")
        if assumption not in {"real", "positive"}:
            raise ConstructionError("Parameter assumptions must be real or positive")
        names[name] = sp.Symbol(name, **{assumption: True})
    return names


def _zero(value):
    return sp.simplify(value) == 0


def ito_transfer(spec):
    """Derive an Ito generator from an explicitly interpreted scalar SDE."""
    names = _symbols(spec)
    convention = spec.get("convention")
    if not isinstance(convention, str) or convention not in {"ito", "stratonovich"}:
        raise ConstructionError("Declare convention explicitly as ito or stratonovich; no default is assumed")
    if not isinstance(spec.get("domain"), str) or spec["domain"] not in {"positive", "real"}:
        raise ConstructionError("Declare domain as positive or real")
    if not isinstance(spec.get("target_domain"), str) or spec["target_domain"] not in {"positive", "real"}:
        raise ConstructionError("Declare target_domain as positive or real")
    x = sp.Symbol("x", **{spec["domain"]: True})
    y = sp.Symbol("y", **{spec["target_domain"]: True})
    names.update(x=x, y=y)
    parse = lambda value: expression(value, names)
    drift, noise = parse(spec["drift"]), parse(spec["noise"])
    h, inverse = parse(spec["state_map"]), parse(spec["inverse_map"])
    if any(y in item.free_symbols for item in (drift, noise, h)) or x in inverse.free_symbols:
        raise ConstructionError("Source expressions use x; inverse_map uses y")
    if not (_zero(h.subs(x, inverse) - y) and _zero(inverse.subs(y, h) - x)):
        raise ConstructionError("The supplied inverse is not established on the declared domain")
    if (getattr(h, "is_" + spec["target_domain"]) is not True
            or getattr(inverse, "is_" + spec["domain"]) is not True):
        raise ConstructionError("The maps are not established within the declared source and target domains")
    convention_correction = (sp.simplify(noise * sp.diff(noise, x) / 2)
                             if convention == "stratonovich" else sp.S.Zero)
    ito_drift = sp.simplify(drift + convention_correction)
    correction = sp.simplify(noise ** 2 * sp.diff(h, x, 2) / 2)
    transformed_drift = sp.simplify(sp.diff(h, x) * ito_drift + correction)
    transformed_variance = sp.simplify((sp.diff(h, x) * noise) ** 2)
    target_drift = sp.simplify(transformed_drift.subs(x, inverse))
    target_variance = sp.simplify(transformed_variance.subs(x, inverse))
    # Both coefficients are needed: agreement of the mean alone is insufficient.
    residuals = [sp.simplify(transformed_drift - target_drift.subs(y, h)),
                 sp.simplify((transformed_variance - target_variance.subs(y, h)) / 2)]
    if not all(_zero(item) for item in residuals):
        raise ConstructionError("The target differential expression did not close")
    candidate = spec.get("candidate")
    candidate_result = None
    if candidate is not None:
        b, v = parse(candidate["drift"]), parse(candidate["variance"])
        if x in b.free_symbols | v.free_symbols:
            raise ConstructionError("Candidate coefficients must use target coordinate y")
        terms = [sp.simplify(transformed_drift - b.subs(y, h)),
                 sp.simplify((transformed_variance - v.subs(y, h)) / 2)]
        candidate_result = {"passes_local_identity": all(_zero(z) for z in terms),
                            "residual_d_phi": str(terms[0]), "residual_d2_phi": str(terms[1])}
    moments = {}
    for order in (1, 2):
        f = y ** order
        moments[str(order)] = str(sp.simplify(target_drift * sp.diff(f, y)
                                  + target_variance * sp.diff(f, y, 2) / 2))
    mean_law = None
    if _zero(sp.diff(target_drift, y, 2)):
        slope = sp.simplify(sp.diff(target_drift, y))
        offset = sp.simplify(target_drift - slope * y)
        t, m0 = sp.symbols("t m0", real=True)
        mean = m0 + offset * t if _zero(slope) else (m0 + offset / slope) * sp.exp(slope * t) - offset / slope
        mean_law = {"expression": str(mean), "initial_mean": "m0",
                    "condition": "finite moments and a boundary realization with no additional moment contribution"}
    return {"calculation": "scalar_ito_transfer", "status": "verified_local_generator_identity",
            "source_convention": convention, "target_convention": "ito",
            "source_ito_drift": str(ito_drift),
            "source_convention_correction": str(convention_correction),
            "domains": {"source": spec["domain"], "target": spec["target_domain"]},
            "target": {"drift": str(target_drift), "variance": str(target_variance)},
            "residual_coefficients": [str(item) for item in residuals],
            "candidate": candidate_result,
            "omission_control": {"removed": "quadratic_variation_drift",
                                 "residual_d_phi": str(correction),
                                 "detects_omission": not _zero(correction)},
            "observable_generators": moments, "mean_prediction": mean_law,
            "scope": "Interior C2 test functions; boundary conditions, units and global stochastic well-posedness require separate checks."}


def quantum_closure(spec):
    """Minimal invariant observable span for finite Hamiltonian dynamics."""
    names = _symbols(spec)
    names["I"] = sp.I

    def matrix(value):
        if not isinstance(value, list) or not value or len(value) > 8:
            raise ConstructionError("Matrices must have between 1 and 8 rows")
        if any(not isinstance(row, list) or len(row) != len(value) for row in value):
            raise ConstructionError("Expected a square matrix")
        return sp.Matrix([[expression(item, names) for item in row] for row in value])

    H = matrix(spec["hamiltonian"])
    O = matrix(spec["observable"])
    if H.shape != O.shape or not all(_zero(z) for z in H - H.conjugate().T):
        raise ConstructionError("Hamiltonian must be Hermitian and match the observable size")
    if not all(_zero(z) for z in O - O.conjugate().T) or all(_zero(z) for z in O):
        raise ConstructionError("Observable must be nonzero and Hermitian")
    n = H.rows
    action = lambda A: (sp.I * (H * A - A * H)).applyfunc(sp.simplify)
    row = lambda A: sp.Matrix(1, n * n, list(A))
    basis = [O]
    index = 0
    while index < len(basis):
        proposed = action(basis[index])
        rows = sp.Matrix.vstack(*(row(A) for A in basis))
        if rows.col_join(row(proposed)).rank() > rows.rank():
            basis.append(proposed)
        index += 1
    rows = sp.Matrix.vstack(*(row(A) for A in basis))
    coefficients = []
    for A in basis:
        solution, free = rows.T.gauss_jordan_solve(row(action(A)).T)
        if free.rows:
            raise ConstructionError("Observable basis unexpectedly dependent")
        coefficients.append(list(solution))
    G = sp.Matrix(coefficients).applyfunc(sp.simplify)
    identities = []
    for i, A in enumerate(basis):
        residual = action(A) - sum((G[i, j] * B for j, B in enumerate(basis)), sp.zeros(n))
        identities.append(all(_zero(z) for z in residual))
    encode = lambda A: [[str(A[i, j]) for j in range(A.cols)] for i in range(A.rows)]
    return {"calculation": "finite_quantum_observable_closure",
            "status": "verified_finite_dimensional_closure", "units": "hbar = 1",
            "hilbert_dimension": n, "observable_dimension": len(basis),
            "basis": [encode(A) for A in basis], "evolution_matrix": encode(G),
            "closed_identities": identities,
            "omission_control": {"single_observable_closed": len(basis) == 1,
                                 "initial_derivative_operator": encode(action(O))},
            "prediction": "For m_i=Tr(rho O_i), m(t)=exp(G t)m(0).",
            "scope": "Exact for the supplied finite Hermitian Hamiltonian and all density matrices; generic parameter rank may drop at exceptional values."}


def verify_construction(spec):
    if spec.get("schema") != "fieldbridge-construction/1":
        raise ConstructionError("Unsupported construction schema")
    for key in ("question", "assumptions", "provenance"):
        if not spec.get(key):
            raise ConstructionError(f"Missing {key}")
    handlers = {"ito_transfer": ito_transfer, "quantum_closure": quantum_closure}
    if spec.get("kind") not in handlers:
        raise ConstructionError("Unsupported calculation kind")
    result = handlers[spec["kind"]](spec)
    return {"schema": "fieldbridge-calculation/1", "question": spec["question"],
            "input_sha256": hashlib.sha256(json.dumps(spec, sort_keys=True).encode()).hexdigest(),
            "assumptions": spec["assumptions"], "provenance": spec["provenance"],
            "source_alignment_verified": False, "novelty_established": False,
            **result}


def render_calculation(report):
    return "# Calculated construction\n\n" + report["question"] + "\n\n```json\n" + json.dumps(report, indent=2) + "\n```\n"


def run_file(input_path, out_dir):
    spec = json.loads(Path(input_path).read_text(encoding="utf-8"))
    report = verify_construction(spec)
    destination = Path(out_dir)
    destination.mkdir(parents=True, exist_ok=True)
    (destination / "calculation.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    (destination / "calculation.md").write_text(render_calculation(report), encoding="utf-8")
    return report
