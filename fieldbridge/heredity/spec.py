"""Specifications of bodies that grow and divide (schema ``fieldbridge-heredity/1``).

A specification has the header of every FieldBridge example (``name``, ``field``, ``question``, ``assumptions``,
``provenance``), a key ``body`` (``equations``, ``field`` or ``chain``, see ``bodies``), the payload of that body, and a
block ``lineage``:

    size            a declared parameter: the extensive size that growth raises and division halves
    threshold       [L_lo, L_hi]: the interval in which the threshold of the symmetric state is sought
    symmetric_state equations: a starting state for the symmetric state; field: a value per variable (one per cell);
                    chain: not used (the straight filament)
    growth          {"law": "exponential" | "linear", "rate": a number or a declared parameter, "mode": "uniform" |
                    "faces" (field only: deposition on both faces)}
    division        {"at": L_div} or {"at_relative": L_div / L_c}; for equations also "variables" (keep or halve per
                    variable), "partition" ("exact" or "binomial") and "volume" (the parameter that counts molecules)
    noise           equations: a number (additive) or one expression per variable (the variance rate); field: the
                    intensity per unit length; chain: k_B T
    dt              the time step (defaults: equations 1e-2, field 4e-3, chain 1e-4)
    start_spread    the spread of the order with which lineages start (default 0.05)
    expect          the class that the source, a closed form or the construction gives
"""
from __future__ import annotations

import copy
import hashlib
import json
from pathlib import Path
from typing import Dict, Tuple, Union

import numpy as np

from ..memory import spec as memory_spec
from .bodies import ChainBody, EquationsBody, FieldBody, Growth, SpecError, constant

SCHEMA = "fieldbridge-heredity/1"
BODIES = ("equations", "field", "chain")
DT = {"equations": 1e-2, "field": 4e-3, "chain": 1e-4}


def _read(source: Union[str, Path, Dict]) -> Tuple[Dict, str]:
    if isinstance(source, dict):
        return copy.deepcopy(source), "<dict>"
    path = Path(source)
    try:
        return json.loads(path.read_text(encoding="utf-8")), str(path)
    except (OSError, json.JSONDecodeError) as error:
        raise SpecError(f"Cannot read {path}: {error}") from error


class Lineage:
    """A body with its protocol of growth and division."""

    def __init__(self, spec: Dict, path: str = "<dict>"):
        if spec.get("schema") != SCHEMA:
            raise SpecError(f"schema must be {SCHEMA!r}")
        for key in ("question", "assumptions", "provenance"):
            if key not in spec:
                raise SpecError(f"Missing field {key!r}")
        kind = spec.get("body")
        if kind not in BODIES:
            raise SpecError(f"body must be one of {', '.join(BODIES)}")
        lin = memory_spec._require(spec, "lineage", dict)
        self.spec, self.path, self.kind = spec, path, kind
        self.size = memory_spec._require(lin, "size", str)
        params = memory_spec._require(spec, "parameters", dict)
        g = memory_spec._require(lin, "growth", dict)
        rate = g.get("rate")
        if isinstance(rate, str):
            if rate not in params:
                raise SpecError("lineage.growth.rate names an undeclared parameter")
            rate = params[rate]
        growth = Growth(g.get("law", "exponential"), memory_spec._number(rate, "lineage.growth.rate"))
        dt = memory_spec._number(lin.get("dt", DT[kind]), "lineage.dt")
        noise = memory_spec._require(lin, "noise")
        division = memory_spec._require(lin, "division", dict)
        if kind == "equations":
            self.body = EquationsBody(spec, self.size, growth, noise, division, dt)
        elif kind == "field":
            self.body = FieldBody(spec, self.size, growth, noise, dt)
        else:
            self.body = ChainBody(spec, self.size, growth, noise, dt)
        self.growth = growth
        lo, hi = (memory_spec._number(v, "lineage.threshold") for v in memory_spec._require(lin, "threshold", list))
        if not 0 < lo < hi:
            raise SpecError("lineage.threshold must be an increasing pair of positive sizes")
        self.search = (lo, hi)
        if ("at" in division) == ("at_relative" in division):
            raise SpecError("lineage.division needs exactly one of 'at' and 'at_relative'")
        self.division_at = memory_spec._number(division["at"], "lineage.division.at") if "at" in division else None
        self.division_relative = (memory_spec._number(division["at_relative"], "lineage.division.at_relative")
                                  if "at_relative" in division else None)
        self.start_spread = memory_spec._number(lin.get("start_spread", 0.05), "lineage.start_spread")
        self.expect = lin.get("expect")
        self.q_sym0 = self._symmetric_guess(lin.get("symmetric_state"))
        self.name = str(spec.get("name", spec["question"]))
        self.field = str(spec.get("field", ""))
        prov = spec["provenance"]
        self.source = str(prov.get("source", "")) if isinstance(prov, dict) else str(prov)
        self.spec_sha256 = hashlib.sha256(json.dumps(spec, sort_keys=True).encode()).hexdigest()
        self.noise_text = json.dumps(noise) if isinstance(noise, dict) else str(noise)

    def _symmetric_guess(self, given) -> np.ndarray:
        b = self.body
        if self.kind == "chain":
            return np.zeros(b.n)
        if self.kind == "field":
            vals = given or {}
            if not isinstance(vals, dict):
                raise SpecError("lineage.symmetric_state of a field gives one value per variable")
            return np.concatenate([np.full(b.N, constant(vals.get(v, b.u0[i]), f"symmetric_state.{v}"))
                                   for i, v in enumerate(b.variables)])
        if given is None:
            return np.zeros(b.n)
        if not isinstance(given, list) or len(given) != b.n:
            raise SpecError("lineage.symmetric_state of an equations body gives one value per variable")
        return np.array([memory_spec._number(v, "lineage.symmetric_state") for v in given])

    def L_div(self, L_c: float) -> float:
        return self.division_at if self.division_at is not None else self.division_relative * L_c


def load(source: Union[str, Path, Dict]) -> Lineage:
    spec, path = _read(source)
    return Lineage(spec, path)
