"""Export the tutorial's equations as a restricted, executable arithmetic tree."""
from __future__ import annotations

import ast
import json
from pathlib import Path

import numpy as np

from .memory.spec import load


def arithmetic(text: str):
    def visit(node):
        if isinstance(node, ast.Constant) and type(node.value) in (int, float):
            return ["number", node.value]
        if isinstance(node, ast.Name):
            return ["name", node.id]
        if isinstance(node, ast.UnaryOp) and isinstance(node.op, (ast.USub, ast.UAdd)):
            return ["neg", visit(node.operand)] if isinstance(node.op, ast.USub) else visit(node.operand)
        if isinstance(node, ast.BinOp):
            ops = {ast.Add: "add", ast.Sub: "sub", ast.Mult: "mul", ast.Div: "div", ast.Pow: "pow"}
            if type(node.op) in ops:
                return [ops[type(node.op)], visit(node.left), visit(node.right)]
        if isinstance(node, ast.Call) and isinstance(node.func, ast.Name):
            if node.func.id in ("sqrt", "exp", "log", "sin", "cos", "tanh", "abs") and len(node.args) == 1 and not node.keywords:
                return [node.func.id, visit(node.args[0])]
        raise ValueError(f"Unsupported browser expression: {text}")
    return visit(ast.parse(text, mode="eval").body)


def additive_terms(node, sign=1):
    if isinstance(node, ast.BinOp) and isinstance(node.op, (ast.Add, ast.Sub)):
        yield from additive_terms(node.left, sign)
        yield from additive_terms(node.right, sign * (-1 if isinstance(node.op, ast.Sub) else 1))
    else:
        text = ast.unparse(node)
        yield text if sign == 1 else f"-({text})"


def export_models(root: Path, gallery: list[dict]) -> list[dict]:
    cards = {c["specification"]: c for c in gallery}
    paths = sorted((root / "examples/memory").glob("*.json"))
    paths += sorted((root / "examples/memory/oscillators").glob("*.json"))
    out = []
    for path in paths:
        spec = json.loads(path.read_text(encoding="utf-8"))
        relative = str(path.relative_to(root))
        out.append(model_record(spec, path.stem, relative, cards.get(relative, {})))
    return out


def model_record(spec: dict, stem: str, relative: str, card: dict | None = None) -> dict:
    """The browser form of one memory specification: arithmetic trees of its drift terms, or its network."""
    card = card or {}
    real = load(spec)
    kind = spec.get("network", {}).get("unit", "equations")
    variables = getattr(real, "variables", [f"q{i+1}" for i in range(real.carrier.dim)])
    if kind == "rotor":
        variables = [f"θ{i+1}" for i in range(real.carrier.dim)]
    terms, links = [], []
    graph = None
    if kind == "equations":
        for i, variable in enumerate(variables):
            for j, text in enumerate(additive_terms(ast.parse(spec["drift"][variable], mode="eval").body)):
                term_id = f"{variable}-{j}"
                names = {n.id for n in ast.walk(ast.parse(text, mode="eval")) if isinstance(n, ast.Name)}
                label = "Decay" if text.replace(" ", "") in (f"-({variable})", f"-{variable}") else "Production" if "/" in text and kind == "equations" else text
                labels = {"pitchfork": ["Linear feedback", "Nonlinear saturation"],
                          "schlogl": ["Reverse cubic reaction", "Autocatalytic production", "Linear removal", "Feed"],
                          "toggle": ["Repressed production", "Dilution"],
                          "toggle_unequal": ["Repressed production", "Dilution"]}
                if stem in labels and j < len(labels[stem]):
                    label = labels[stem][j]
                terms.append({"id": term_id, "variable": i, "label": f"{variable}: {label.lower()}", "text": text, "tree": arithmetic(text)})
                for source, name in enumerate(variables):
                    if name in names:
                        links.append({"source": source, "target": i, "term": term_id})
    elif kind == "gene":
        graph = real.network
        links = [{"source": i, "target": j, "term": f"edge-{k}", "sign": sign}
                 for k, (i, j, sign) in enumerate(graph["edges"])]
        terms = [{"id": f"edge-{k}", "label": f"{variables[i]} {'activates' if sign > 0 else 'represses'} {variables[j]}", "text": "Hill regulation"}
                 for k, (i, j, sign) in enumerate(graph["edges"])]
        terms.append({"id": "decay", "label": "Dilution", "text": "−qᵢ"})
    elif kind == "rotor":
        g = real.graph
        graph = {k: np.asarray(g[k]).tolist() for k in ("src", "tgt", "phi", "wa", "wc")}
        graph["m"] = g["m"]
        terms = [{"id": "align", "label": "Alignment", "text": "−J cos[m(θᵢ − θⱼ)]"},
                 {"id": "reflect", "label": "Bond-directional coupling", "text": "−G cos[m(θᵢ + θⱼ − 2φᵢⱼ)]"}]
        links = [{"source": i, "target": j, "term": f"edge-{k}"} for k, (i, j) in enumerate(zip(graph["src"], graph["tgt"]))]
    else:
        raise ValueError(f"No browser model for {kind}")
    ranges = {}
    for name, value in real.params.items():
        if name == real.control:
            lo, hi = real.control_range
        elif name == "n":
            lo, hi = 1, 6
        elif value > 0:
            lo, hi = max(0.01, value * 0.25), value * 2
        else:
            lo, hi = -1, 1
        ranges[name] = [lo, hi, 1 if name == "n" else (hi - lo) / 100]
    n, scale = real.carrier.dim, real.carrier.scale
    initial = [scale * (0.8 if i % 2 == 0 else 0.2) for i in range(n)]
    alternate = list(reversed(initial)) if n > 1 else [scale * 0.2]
    if real.carrier.kind == "euclid":
        initial = [0.4] * n
        alternate = [-0.4] + [0.4] * (n - 1)
    if real.carrier.kind == "torus":
        initial = [real.carrier.period * (i + 1) / (n + 2) for i in range(n)]
        alternate = [v + real.carrier.period / 3 for v in initial]
    if stem == "compartments":
        initial, alternate = [0.2, 0.1], [0.2, -0.1]
    return {"id": stem, "name": spec["name"], "field": spec.get("field", "Model systems"),
            "question": spec["question"], "kind": kind, "variables": variables,
            "carrier": {"kind": real.carrier.kind, "name": real.carrier.name, "period": real.carrier.period, "scale": scale},
            "params": real.params, "ranges": ranges, "terms": terms, "links": links, "graph": graph,
            "positions": (real.geometry or {}).get("pos"), "initial": initial, "alternate": alternate,
            "dt": min(real.dt, 0.002 if kind == "rotor" else 0.01),
            "duration": 30 if "oscillators" in relative or kind == "gene" else 15,
            "noise": real.noise, "closure": spec.get("closure", real.closure), "observable": real.observable,
            "assumptions": spec["assumptions"], "source": spec["provenance"]["source"], "source_url": card.get("source_url"),
            "specification": relative,
            "tutorial": "docs/tutorial/20_memory_phase.md" if "oscillators" in relative else "docs/tutorial/21_memory_new_material.md#3-the-gallery-of-realizations",
            "saved": card or None}
