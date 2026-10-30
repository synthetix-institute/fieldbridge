"""The loop from a retained law departure to a tested description of smaller dimension."""
import json
from pathlib import Path

import numpy as np
import pytest

ln = pytest.importorskip("fieldbridge.memory.learning", reason="the learning modules are not committed yet")
rf = pytest.importorskip("fieldbridge.memory.repair_families", reason="the learning modules are not committed yet")
from fieldbridge.memory import repair_loop as rl  # noqa: E402
from fieldbridge.memory import spec as ms  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
RECORDS = ROOT / "docs/records/quantum_semiclassical.jsonl"
FAST = rl.OrderPolicy(max_order=5, radii=(0.25, 0.125, 0.0625), per_radius=2)


def model(drift, parameters, variables, kind="euclid", **carrier):
    return {"schema": "fieldbridge-memory/1", "kind": "equations", "name": "test", "question": "test",
            "assumptions": ["test"], "provenance": {"source": "test"}, "noise": 0, "dt": 0.01,
            "carrier": {"kind": kind, "variables": variables, **carrier}, "parameters": parameters, "drift": drift,
            "control": {"name": list(parameters)[0], "range": [-1, 1]}}


FAST_SLOW = model({"x": "e*x-x*z", "z": "-z+x**2"}, {"e": 0.0}, ["x", "z"])


def closure(spec=FAST_SLOW, state=(0.0, 0.0), order=5, swept=("e",), at=None):
    real = ms.load(spec)
    return rl._solve_closure(real, list(state), at or {"e": 0.0}, list(swept), order, rl.OrderPolicy())


def records():
    return [json.loads(line) for line in RECORDS.read_text().splitlines() if line.strip()]


def test_the_exact_identity_fails_where_the_closure_modulo_order_holds():
    # the exact search cannot close x alone: at best it adds source coordinates and restates the model
    task = {"schema": ln.SCHEMA, "model": FAST_SLOW, "coordinates": {"u": "x"}, "candidate_drift": {"u": "e*x"},
            "development": [{"id": f"d{i}", "initial": [0.2 + 0.1 * i, 0.01 * i], "parameters": {"e": 0.1 + 0.05 * i},
                             "duration": 2.0} for i in range(3)],
            "withheld": [{"id": f"w{i}", "initial": [-0.2 - 0.1 * i, 0.02 * i], "parameters": {"e": -0.1 - 0.05 * i},
                          "duration": 2.0} for i in range(3)]}
    report = ln.learn(task, None, ln.Policy(max_added_coordinates=1, max_candidates=8))
    assert report["status"] != "verified_extension" or rl.assess_exact_report(task, report)["verdict"] != "smaller"
    c = closure()
    assert rl.verify_modulo_order({"model": FAST_SLOW}, c, rl.OrderPolicy())["passed"]


def test_the_reduced_drift_matches_the_slow_manifold_on_the_fibre_of_z_zero():
    # x' = e x - x z, z' = -z + x^2: with h(x, 0, e) = x the cubic term is -x^3 and the critical rate is e
    g = closure()["g"]
    assert g[(1, 1)] == pytest.approx(1.0, abs=1e-12)
    assert g[(3, 0)] == pytest.approx(-1.0, abs=1e-12)
    assert g.get((2, 0), 0.0) == pytest.approx(0.0, abs=1e-12)


@pytest.mark.parametrize("order", [1, 2, 3, 5])
def test_the_verifier_reports_the_first_omitted_order(order):
    check = rl.verify_modulo_order({"model": FAST_SLOW}, rl.truncate(closure(order=6), order), rl.OrderPolicy())
    assert check["passed"]
    assert check["first_omitted_order"] == order + 1


def test_the_verifier_rejects_an_altered_coefficient():
    c = closure()
    altered = {**c, "g": {**c["g"], (3, 0): -1.1}}
    assert not rl.verify_modulo_order({"model": FAST_SLOW}, altered, rl.OrderPolicy())["passed"]
    c2 = dict(c["h"])
    key = next(k for k in c2 if sum(k) == 2)
    assert not rl.verify_modulo_order({"model": FAST_SLOW}, {**c, "h": {**c2, key: c2[key] + 0.2}},
                                      rl.OrderPolicy())["passed"]


def test_a_preparation_away_from_the_centre_manifold_is_predicted():
    # the fibres of h carry the transient: z(0) is far from x^2
    c = closure(order=7)
    ev = rl.Evaluator(c)
    case = {"id": "off", "initial": [0.15, 0.1], "parameters": {"e": 0.02}, "duration": 30.0, "radius": 1.0}
    real = ms.load(FAST_SLOW)
    runs = rl._source_trajectories(real, [case])
    row = rl._rows(c, [case], runs, rl.OrderPolicy())[0]
    plain = rl._rows(rl._plain_projection(c), [case], runs, rl.OrderPolicy())[0]
    assert row["error"] < 1e-3
    assert plain["error"] > 10 * row["error"]
    assert ev.h(np.array(case["initial"]), case["parameters"]) != pytest.approx(0.15)


def test_a_nonpolynomial_drift_outside_the_polynomial_backend_is_reduced():
    spec = model({"x": "e*sin(x)-z*x", "z": "-z+1-cos(x)"}, {"e": 0.0}, ["x", "z"])
    with pytest.raises(ln.ConstructionError):
        ln._parse_task({"schema": ln.SCHEMA, "model": spec, "coordinates": {"u": "x"}, "candidate_drift": {"u": "x"}})
    c = closure(spec)
    assert rl.verify_modulo_order({"model": spec}, c, rl.OrderPolicy())["passed"]


def test_more_than_four_source_coordinates():
    # the exact backend stops at four source coordinates
    names = ["x", "a", "b", "c", "d", "f"]
    drift = {"x": "e*x-x*a", "a": "-a+x**2", "b": "-2*b+a*x", "c": "-3*c+b", "d": "-d+c*x", "f": "-4*f+d"}
    spec = model(drift, {"e": 0.0}, names)
    with pytest.raises(ln.ConstructionError):
        ln._parse_task({"schema": ln.SCHEMA, "model": spec, "coordinates": {"u": "x"}, "candidate_drift": {"u": "x"}})
    c = closure(spec, state=[0.0] * 6, order=4)
    assert rl.verify_modulo_order({"model": spec}, c, rl.OrderPolicy())["passed"]
    assert c["g"][(3, 0)] == pytest.approx(-1.0, abs=1e-10)


def test_a_torus_carrier_in_its_covering_space():
    spec = model({"t1": "e*sin(t1)+sin(t2)*sin(t1)", "t2": "-sin(t2)+1-cos(t1)"}, {"e": 0.0}, ["t1", "t2"],
                 kind="torus", period=6.283185307179586)
    c = closure(spec)
    assert rl.verify_modulo_order({"model": spec}, c, rl.OrderPolicy())["passed"]


def test_an_orthant_point_on_the_boundary_is_reported_and_an_interior_one_reduced():
    spec = model({"x": "e*x+x*z-x**2", "z": "-z+x"}, {"e": 0.0}, ["x", "z"], kind="orthant")
    real = ms.load(spec)
    assert rl._carrier_limit(real, np.array([0.0, 0.0]))
    assert rl._carrier_limit(real, np.array([0.5, 0.5])) is None


def test_a_second_critical_direction_and_a_complex_pair_are_refused():
    spec = model({"x": "e*x-x**3", "y": "e*y-y**3"}, {"e": 0.0}, ["x", "y"])
    with pytest.raises(ValueError, match="second eigenvalue"):
        closure(spec)
    spec = model({"x": "e*x-y", "y": "x+e*y"}, {"e": 0.0}, ["x", "y"])
    with pytest.raises(ValueError, match="complex"):
        closure(spec)


def test_a_resonance_of_the_non_critical_modes_is_reported():
    # eigenvalues 0, +1, -1: the monomial s1 s2 has a vanishing denominator
    spec = model({"x": "x*a*b+e*x", "a": "a", "b": "-b"}, {"e": 0.0}, ["x", "a", "b"])
    with pytest.raises(ValueError, match="resonance"):
        closure(spec, state=[0.0, 0.0, 0.0], order=3)


# ------------------------------------------------------------------------------------------- the size rule
@pytest.mark.parametrize("source,retained,rank,verdict", [
    (2, 1, 1, "smaller"), (2, 2, 2, "restates_model"), (2, 3, 2, "restates_model"), (3, 3, 2, "not_smaller"),
    (3, 2, 2, "smaller")])
def test_the_size_rule(source, retained, rank, verdict):
    assert rl.description_size(source, retained, rank)["verdict"] == verdict


def test_a_verified_extension_that_reattaches_every_source_coordinate_restates_the_model():
    spec = model({"x": "-a*x+y", "y": "-y"}, {"a": 1.0}, ["x", "y"])
    task = {"schema": ln.SCHEMA, "model": spec, "coordinates": {"u": "x"}, "candidate_drift": {"u": "-a*u"},
            "development": [{"id": f"d{i}", "initial": [1.0 + i, 1.0 - 0.2 * i], "parameters": {"a": 1.0 + 0.1 * i},
                             "duration": 2.0} for i in range(3)],
            "withheld": [{"id": f"w{i}", "initial": [-1.0 - i, 2.0 + 0.2 * i], "parameters": {"a": 1.5 + 0.1 * i},
                          "duration": 2.0} for i in range(3)]}
    report = ln.learn(task, None, ln.Policy(max_degree=1))
    assert report["status"] == "verified_extension"
    assert rl.assess_exact_report(task, report)["verdict"] == "restates_model"


def test_a_verified_extension_that_stays_below_the_source_dimension_is_smaller():
    # r = x^2 + y^2 closes on itself in a rotation with damping
    spec = model({"x": "-a*x-y", "y": "x-a*y"}, {"a": 1.0}, ["x", "y"])
    task = {"schema": ln.SCHEMA, "model": spec, "coordinates": {"r": "x**2+y**2"}, "candidate_drift": {"r": "0*r"},
            "development": [{"id": f"d{i}", "initial": [1.0 + i, 1.0 - 0.2 * i], "parameters": {"a": 0.1 + 0.1 * i},
                             "duration": 2.0} for i in range(3)],
            "withheld": [{"id": f"w{i}", "initial": [-1.0 - i, 2.0 + 0.2 * i], "parameters": {"a": 0.5 + 0.1 * i},
                          "duration": 2.0} for i in range(3)]}
    report = ln.learn(task, None, ln.Policy(max_degree=1))
    assert report["status"] == "verified_extension"
    assert rl.assess_exact_report(task, report)["verdict"] == "smaller"


# ------------------------------------------------------------------------------------------- generator and protocol
def certified_case():
    for record in records():
        if record["source"].endswith("cim_dopo.json"):
            return record


def test_the_generator_writes_a_valid_task_from_a_record():
    record = certified_case()
    task = rl.learning_task(record, "completion/0/law", root=ROOT, policy=FAST)
    assert task["schema"] == ln.SCHEMA and task["model"]["noise"] == 0
    assert task["noise_set_to_zero"] == json.loads((ROOT / record["source"]).read_text())["noise"]
    real, x, p, F, z, h, G = ln._parse_task(task)
    # the exact path accepts the same task; the proposed drift carries the linear-stage terms only
    assert len(z) == 1 and len(G) == 1
    development = {tuple(r["initial"]) for r in task["development"]}
    assert not development & {tuple(r["initial"]) for r in task["withheld"]}
    assert task["failed_law"]["target"] == "completion/0/law"
    assert rl.learning_task(record, "completion/0/law", root=ROOT, policy=FAST)["withheld"] == task["withheld"]


def test_an_unknown_law_of_a_record_is_an_error():
    with pytest.raises(KeyError):
        rl.learning_task(certified_case(), "completion/9/law", root=ROOT, policy=FAST)


def test_a_withheld_set_is_consumed_once_and_a_failure_is_kept(tmp_path):
    library = ln.MoveLibrary(tmp_path / "library.sqlite")
    task = rl.learning_task(certified_case(), "completion/0/law", root=ROOT, policy=FAST)
    first = rl.close_modulo_order(task, library, FAST)
    assert first["novelty_established"] is False
    assert first["status"] != "unsupported_or_resonant"
    second = rl.close_modulo_order(task, library, FAST)
    assert second["status"] == "withheld_already_used"
    with library.connect() as db:
        attempts = db.execute("SELECT COUNT(*) FROM attempts").fetchone()[0]
    assert attempts == 2


def test_a_closure_does_not_enter_the_recipes_of_the_exact_search(tmp_path):
    library = ln.MoveLibrary(tmp_path / "library.sqlite")
    recipe = {"engine": rl.ORDER_ENGINE, "source_dimension": 2, "retained_dimension": 1, "order": 3, "reach": 1.0,
              "kind": "nonlinear stage"}
    library.store({"status": "verified_extension"}, recipe)
    assert library.recipes() == []
    assert rl.order_recipes(library) == [recipe]


def test_every_retained_outcome_of_the_records_gets_a_row(tmp_path):
    library = ln.MoveLibrary(tmp_path / "library.sqlite")
    for record in records():
        ln.retain_record(library, record, source=str(RECORDS))
    rows = rl.close_departures(library, root=ROOT, policy=FAST, select=set())
    assert len(rows) == len(library.departures()) >= 227
    assert not [r for r in rows if r["status"] == "needs_equation_task"]
    failed = [r for r in rows if r["status"] == "outside_this_campaign"]
    assert sum(r["write_class"] == "nonlinear stage" for r in failed) == 9
    assert sum(r["write_class"] == "equilibrium write" for r in failed) == 4


def test_the_exact_baseline_and_the_loop_do_not_share_a_withheld_key():
    record = certified_case()
    task = rl.learning_task(record, "completion/0/law", root=ROOT, policy=FAST)
    keys = {ln.digest({k: r[k] for k in ("initial", "parameters")}) for r in task["withheld"]}
    assert len(keys) == len(task["withheld"])


def test_the_power_control_finds_the_term_that_shapes_the_trajectory():
    c = closure(order=5)
    key = rl._strongest_nonlinear(c, 0.5, 1.0, 0.1)
    assert key == (3, 0)
