"""Specifications of populations that decide (schema ``fieldbridge-decision/1``).

A specification has the header of every FieldBridge example (``name``, ``field``, ``question``, ``assumptions``,
``provenance``), a ``kind``, the equations payload of the memory module (``carrier``, ``parameters``, ``drift``) where
the kind needs one, a block ``population`` and a block ``protocol``:

kind ``oscillators`` (target synchronization): the payload is one limit-cycle unit;
    population.heterogeneity  {parameter, density: gaussian, centre, spread}: spread is the frequency spread
                              sigma_omega / omega that the parameter spread produces
    population.coupling       {through: c(x), acting_on: variable, factors: [K/K_c, ...]}

kind ``collective`` (targets collective-write and seeded-passage): the payload is the mean-field equations of the
population (the drift may be left out when transitions are given: F = sum nu a);
    population.size           a declared parameter: N (units, molecules per source volume, or a stability factor)
    population.control        {name, range}: the symmetric state loses stability inside the range
    population.bias           {name}
    population.noise          {transitions: [{change: {variable: +-1}, rate: expr}]} (D = (1/2N) sum nu nu^T a)
                              or {variance: {variable: expr}} (D = diag(expr) / 2N)
    population.symmetric_state  a starting point for the symmetric state

kind ``units`` (target collective-write): N units with a diverse parameter coupled through their mean;
    population.unit           {variable, mean, drift}: the unit drift in its variable, the diverse parameter and the mean
    population.heterogeneity  {parameter, density: gaussian, centre, width}
    population.control, population.bias as above; population.noise a number D (dx = f dt + sqrt(2 D) dW)

protocol (target collective-write): {"sweep": {"from", "to", "rates": [...]}, "sizes": [...], "z": [...],
    "sample": "quantiles" | "random" (units)}; z are biases in units of the thermal spread of the law.
protocol (target seeded-passage): {"step": {"from", "to"}, "threshold", "duration", "sizes": [...], "dt"}.
expect: the class that the source, a closed form or the construction gives.
"""
from __future__ import annotations

import copy
import hashlib
import json
from pathlib import Path
from typing import Dict, Tuple, Union

from ..memory import spec as memory_spec
from .collective import Collective
from .oscillators import Oscillators
from .units import Units

SCHEMA = "fieldbridge-decision/1"
SpecError = memory_spec.SpecError
KINDS = {"oscillators": Oscillators, "collective": Collective, "units": Units}
TARGETS = {"oscillators": ("synchronization",), "collective": ("collective-write", "seeded-passage"),
           "units": ("collective-write",)}
# spec key -> slot of the identity (used by the page's edit check)
SLOT_KEYS = {"drift": "Omega", "carrier": "Xi", "population.coupling": "C", "population.noise": "C",
             "population.heterogeneity": "Xi", "population.size": "Xi", "population.control": "P", "protocol": "P",
             "population.bias": "A"}


def _read(source: Union[str, Path, Dict]) -> Tuple[Dict, str]:
    if isinstance(source, dict):
        return copy.deepcopy(source), "<dict>"
    path = Path(source)
    try:
        return json.loads(path.read_text(encoding="utf-8")), str(path)
    except (OSError, json.JSONDecodeError) as error:
        raise SpecError(f"Cannot read {path}: {error}") from error


class Population:
    def __init__(self, spec: Dict, path: str = "<dict>"):
        if spec.get("schema") != SCHEMA:
            raise SpecError(f"schema must be {SCHEMA!r}")
        for key in ("question", "assumptions", "provenance"):
            if key not in spec:
                raise SpecError(f"Missing field {key!r}")
        kind = spec.get("kind")
        if kind not in KINDS:
            raise SpecError(f"kind must be one of {', '.join(KINDS)}")
        self.target = spec.get("target", TARGETS[kind][0])
        if self.target not in TARGETS[kind]:
            raise SpecError(f"a {kind} population serves the target(s) {', '.join(TARGETS[kind])}")
        self.kind, self.spec, self.path = kind, spec, path
        self.body = KINDS[kind](spec)
        self.protocol = spec.get("protocol", {})
        if self.target == "collective-write" and getattr(self.body, "bias", None) is None:
            raise SpecError("the collective write needs population.bias")
        if self.target == "collective-write":
            sw = memory_spec._require(self.protocol, "sweep", dict)
            self.sweep_from = memory_spec._number(memory_spec._require(sw, "from"), "protocol.sweep.from")
            self.sweep_to = memory_spec._number(memory_spec._require(sw, "to"), "protocol.sweep.to")
            self.rates = [memory_spec._number(r, "protocol.sweep.rates") for r in memory_spec._require(sw, "rates",
                                                                                                     list)]
            if not self.sweep_from < self.sweep_to or min(self.rates) <= 0:
                raise SpecError("protocol.sweep: from < to and positive rates")
            self.sample = self.protocol.get("sample", "quantiles")
            if self.sample not in ("quantiles", "random") or (self.sample == "random" and kind != "units"):
                raise SpecError("protocol.sample: quantiles, or random for units")
        elif self.target == "seeded-passage":
            st = memory_spec._require(self.protocol, "step", dict)
            self.step_from = memory_spec._number(memory_spec._require(st, "from"), "protocol.step.from")
            self.step_to = memory_spec._number(memory_spec._require(st, "to"), "protocol.step.to")
            self.threshold = memory_spec._number(memory_spec._require(self.protocol, "threshold"),
                                                 "protocol.threshold")
            self.duration = memory_spec._number(memory_spec._require(self.protocol, "duration"), "protocol.duration")
            if self.threshold <= 0 or self.duration <= 0:
                raise SpecError("protocol.threshold and protocol.duration must be positive")
        if self.target != "synchronization":
            self.sizes = [memory_spec._number(n, "protocol.sizes") for n in memory_spec._require(self.protocol,
                                                                                                   "sizes", list)]
            if min(self.sizes) <= 0:
                raise SpecError("protocol.sizes must be positive")
        self.z = [memory_spec._number(z, "protocol.z") for z in self.protocol.get("z", [0.8])]
        self.dt = memory_spec._number(self.protocol.get("dt", 1e-2), "protocol.dt")
        self.expect = spec.get("expect")
        self.name = str(spec.get("name", spec["question"]))
        self.field = str(spec.get("field", ""))
        prov = spec["provenance"]
        self.source = str(prov.get("source", "")) if isinstance(prov, dict) else str(prov)
        self.spec_sha256 = hashlib.sha256(json.dumps(spec, sort_keys=True).encode()).hexdigest()


def load(source: Union[str, Path, Dict]) -> Population:
    spec, path = _read(source)
    return Population(spec, path)
