"""The four review checks (fieldbridge/review): each check passes a sound claim and fails the fault it is meant to
find, and the worked example finds the faults of arXiv:2404.07492 recorded in its review file."""
import json
from pathlib import Path

from fieldbridge.review import checks, theorems
from fieldbridge.review.cli import add_parser, cmd_check

ROOT = Path(__file__).resolve().parents[1]
SYSTEMS = {
    "fixed": {"carrier": {"variables": ["positions"], "fixed_parameters": ["radii"]}, "closure": ["periodic"],
              "realization": {"radii": "given"}},
    "free": {"carrier": {"variables": ["positions", "radii"], "free_parameters": ["radii"]}, "closure": ["periodic"],
             "realization": {"radii": "tuned"}},
    "tuned": {"carrier": {"variables": ["positions"], "fixed_parameters": ["radii"]}, "closure": ["periodic"],
              "realization": {"radii": "tuned"}, "counts": {"equations": "3*N", "unknowns": "2*N + 1"}},
}


def test_same_system_passes_on_one_system_and_fails_across_systems():
    assert checks.same_system({"computed_on": "fixed", "claimed_for": "fixed"}, SYSTEMS)["status"] == "pass"
    r = checks.same_system({"computed_on": "free", "claimed_for": "fixed"}, SYSTEMS)
    assert r["status"] == "fail" and "radii" in r["explanation"]


def test_overdetermined_state_fails_a_transfer_and_warns_on_its_own_system():
    r = checks.same_system({"computed_on": "free", "claimed_for": "tuned"}, SYSTEMS)
    assert r["status"] == "fail" and r["counts"][-1]["overdetermined"]
    assert checks.same_system({"computed_on": "tuned", "claimed_for": "tuned"}, SYSTEMS)["status"] == "warn"


def test_protocol_independence_detects_a_trend_with_run_length():
    flat = {"protocol": {"quantity": "T_m", "setting": "run length", "estimates": [
        {"setting": s, "value": 0.002, "error": 1e-4} for s in (1e4, 1e5, 1e6)]}}
    assert checks.protocol_independence(flat)["status"] == "pass"
    drift = {"protocol": {"quantity": "T_m", "setting": "run length", "estimates": [
        {"setting": s, "value": v, "error": 5e-5} for s, v in ((1e4, 0.004), (1e5, 0.003), (1e6, 0.002))]}}
    assert checks.protocol_independence(drift)["status"] == "fail"


def test_a_divergence_inside_the_observation_window_fails():
    c = {"protocol": {"window": {"estimate": 0.0029, "warmest_without_relaxation": 0.0025,
                                 "coldest_with_relaxation": 0.003}}}
    assert checks.protocol_independence(c)["status"] == "fail"
    c = {"protocol": {"window": {"estimate": 0.001, "warmest_without_relaxation": 0.0005,
                                 "coldest_with_relaxation": 0.003}}}
    assert checks.protocol_independence(c)["status"] == "warn"
    c["protocol"]["window"]["error"] = 1e-4
    assert checks.protocol_independence(c)["status"] == "fail"
    c = {"protocol": {"window": {"estimate": 0.001, "warmest_without_relaxation": 0.0009,
                                 "coldest_with_relaxation": 0.0011}}}
    assert checks.protocol_independence(c)["status"] == "warn"


def test_thermodynamic_claim_from_a_kinetic_observable_needs_equilibrium_evidence():
    c = {"kind": "thermodynamic", "observable_type": "relaxation time"}
    assert checks.protocol_independence(c)["status"] == "fail"
    c["equilibrium_evidence"] = "free energies of liquid and solid"
    assert checks.protocol_independence(c)["status"] == checks.NA


def test_theorem_scope_hypotheses_and_inferences():
    ok = theorems.check_use("mermin_wagner_hohenberg", {"dimension": 2, "short_range_interactions": True,
                                                        "temperature_positive": True, "continuous_symmetry": True},
                            "no_long_range_order_of_continuous_symmetry")
    assert ok["status"] == "pass"
    bad = theorems.check_use("mermin_wagner_hohenberg", {"dimension": 2, "short_range_interactions": True,
                                                         "temperature_positive": True, "continuous_symmetry": True},
                             "zero_melting_temperature")
    assert bad["status"] == "fail"
    outside = theorems.check_use("bender_gao_wormald_local_limit",
                                 {"graph_class": "labelled 2-connected planar", "edges_per_vertex": 3.0})
    assert outside["status"] == "fail"
    inside = theorems.check_use("bender_gao_wormald_local_limit",
                                {"graph_class": "labelled 2-connected planar", "edges_per_vertex": 2.0},
                                "gaussian_count_inside_interval")
    assert inside["status"] == "pass"
    assert theorems.check_use("no_such_theorem", {})["status"] == "unknown"


def test_prior_art_requires_a_recorded_search():
    assert checks.prior_art({"novelty_claimed": True})["status"] == "fail"
    good = {"novelty_claimed": True, "prior_art": [{"statement": "s", "queries": ["q1", "q2"], "verdict": "not found"}]}
    assert checks.prior_art(good)["status"] == "pass"
    weak = {"novelty_claimed": True, "prior_art": [{"statement": "s", "queries": ["q1"], "verdict": "not found"}]}
    assert checks.prior_art(weak)["status"] == "warn"
    taken = {"novelty_claimed": True, "prior_art": [{"statement": "s", "queries": ["q"], "verdict": "published",
                                                     "sources": ["X 2020"]}]}
    assert checks.prior_art(taken)["status"] == "fail"
    contra = {"prior_art": [{"statement": "s", "queries": ["q"], "verdict": "published", "sources": ["Y"],
                             "contradicts": True}]}
    assert checks.prior_art(contra)["status"] == "fail"


def test_worked_example_finds_the_recorded_faults(tmp_path):
    doc = json.loads((ROOT / "examples/review/ideal_glass_2404.07492.json").read_text())
    rep = checks.review(doc)
    status = {c["id"]: {r["check"]: r["status"] for r in c["results"]} for c in rep["claims"]}
    assert status["C1"]["same_system"] == "fail" and status["C1"]["theorem_scope"] == "fail"
    assert status["C2"]["prior_art"] == "fail"
    assert status["C3"]["protocol_independence"] == "fail"
    assert status["C4"]["theorem_scope"] == "fail"
    assert status["C5"]["protocol_independence"] == "fail"

    class A:
        review = str(ROOT / "examples/review/ideal_glass_2404.07492.json")
        out_dir = str(tmp_path)
        strict = True
    assert cmd_check(A()) == 1
    assert (tmp_path / "review.json").exists() and (tmp_path / "review.md").exists()


def test_review_command_is_mounted():
    import argparse
    sub = argparse.ArgumentParser().add_subparsers()
    add_parser(sub)
