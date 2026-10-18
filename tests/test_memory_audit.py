"""Audit of a derived writing law (memory audit): each specific possible error and its test."""
import json
from pathlib import Path

import numpy as np

from fieldbridge.memory import audit, spec as mspec

ROOT = Path(__file__).resolve().parents[1]
PITCHFORK = json.loads((ROOT / "examples/memory/pitchfork.json").read_text(encoding="utf-8"))


def negative_rate_pitchfork():
    s = dict(PITCHFORK)
    s.pop("potential", None)
    s.update(name="pitchfork written only at a negative rate", drift={"x": "-(k + 1)*x - x**3"},
             parameters={"k": 1.0}, control={"name": "k", "range": [-3.0, 2.0]}, material=[])
    return mspec.load(s)


def test_a_write_at_a_rate_of_the_wrong_sign_is_confirmed_and_one_of_the_right_sign_refuted():
    real = negative_rate_pitchfork()
    bad = audit.check_published_sign(real, real, {"value": -1.0}, None, {})
    assert bad["status"] == audit.CONFIRMED
    good = audit.check_published_sign(real, real, {"value": 0.5}, None, {})
    assert good["status"] == audit.REFUTED


def test_the_few_quanta_check_applies_only_to_semiclassical_descriptions_of_quantum_elements():
    real = mspec.load(PITCHFORK)
    item = audit.check_few_quanta(real, real, {"state": [0.0], "value": 0.0, "kappa_slope": -1.0}, None, {"law": {}})
    assert item["status"] == audit.NOT_APPLICABLE
    kerr = mspec.load(json.loads((ROOT / "examples/quantum/semiclassical/kerr_cat_kpo.json").read_text()))
    assert audit._quantum_description(kerr)


def test_the_pitchfork_is_ready_for_review():
    """A valid law must pass: no check may be confirmed on the canonical pitchfork."""
    res = audit.audit_symmetric_write(mspec.load(PITCHFORK), np.random.default_rng(5), n_traj=800, replicates=2)
    assert res["status"].startswith("reached")
    assert res["verdict"] == "ready for review", res["ledger"]
    statuses = {x["check"]: x["status"] for x in res["ledger"]}
    assert statuses["window"] == audit.REFUTED and statuses["one sign"] == audit.REFUTED
    assert statuses["few quanta"] == audit.NOT_APPLICABLE
