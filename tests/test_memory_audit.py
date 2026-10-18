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


def test_a_flux_phase_is_compared_modulo_two_pi_and_a_rate_is_not():
    """The external flux phase of a fluxonium is defined modulo 2 pi together with the phase of the loop: -3 pi at the
    write is the published pi. The rate of the pitchfork has no such partner, so -1 stays a wrong sign."""
    flux = mspec.load(json.loads((ROOT / "examples/quantum/semiclassical/fluxonium.json").read_text()))
    assert audit.periodic_partner(flux, "phe") == 0
    item = audit.check_published_sign(flux, flux, {"value": -3 * np.pi, "_control_at_write": -3 * np.pi}, None, {})
    assert flux.control == "phe" and item["status"] == audit.REFUTED
    (row,) = [v for k, v in item["evidence"].items() if k.startswith("phe")]
    assert np.isclose(row[2], np.pi)
    real = negative_rate_pitchfork()
    assert audit.periodic_partner(real, "k") is None


def test_a_write_whose_law_was_not_measured_is_not_ready_for_review():
    refuted = [{"check": "window", "status": audit.REFUTED, "error": "e"}]
    assert audit.verdict(refuted, {}).startswith("not testable")
    assert audit.verdict(refuted, {"constant": 1.33}) == "ready for review"
    confirmed = refuted + [{"check": "one sign", "status": audit.CONFIRMED, "error": "offset"}]
    assert audit.verdict(confirmed, {}) == "errors found: offset"


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
