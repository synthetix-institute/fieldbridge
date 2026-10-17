"""The balanced counterpart of a hysteron network and the states measured against it.

The counterpart reverses the smallest set of signs, couplings or actions of the drive, that makes the network with
the drive balanced; by the structural rule its return to every turning point is then exact. The controls: the
enumeration of states reproduces the Preisach counts and stays closed under the dynamics of Run, the count does not
depend on a relabeling, the counterpart keeps magnitudes, fields and thresholds, local search never undercuts the exact
index, and the measured identities C = C_bal + gain and C_bal = C_0 + loss hold.
"""
import json
import subprocess
import sys
from pathlib import Path

import pytest

np = pytest.importorskip("numpy")
pytest.importorskip("scipy")

from fieldbridge.memory import counterpart as cp
from fieldbridge.memory import hysterons as hy

ROOT = Path(__file__).resolve().parents[1]
EX = ROOT / "examples" / "counterpart"


def make(J, h, b, eta, name="test"):
    J = np.asarray(J, float)
    n = J.shape[0]
    bonds = [(i, j, float(J[i, j])) for i in range(n) for j in range(n) if i != j and J[i, j] != 0]
    return hy.Hysterons(name=name, n=n, bonds=bonds, h=np.asarray(h, float), b=np.asarray(b, float),
                        eta=np.asarray(eta, float), reciprocal=bool(np.allclose(J, J.T)))


def random_network(rng, n=8, scale=0.6):
    A = rng.normal(0.0, scale, (n, n))
    J = np.triu(A, 1)
    return make(J + J.T, rng.normal(0.0, 1.0, n), rng.uniform(0.1, 1.0, n), np.ones(n))


# ------------------------------------------------------------------------------------------------ enumeration
def test_the_enumeration_reproduces_the_preisach_counts():
    # loops one inside the other reach all four states; shifted loops cannot hold the second up and the first down
    inside = make(np.zeros((2, 2)), [0.0, -1.0], [1.0, 2.5], [1.0, 1.0])
    shifted = make(np.zeros((2, 2)), [0.0, -1.0], [1.0, 0.5], [1.0, 1.0])
    assert cp.reachable(inside)["count"] == 4 and cp.reachable(shifted)["count"] == 3
    rng = np.random.default_rng(2)
    for n in (3, 5, 7):
        free = make(np.zeros((n, n)), rng.normal(0, 1, n), rng.uniform(0.1, 1.5, n), np.ones(n))
        tg = cp.reachable(free)
        assert n + 1 <= tg["count"] <= 2 ** n
        assert cp.exact_return(tg)["failed"] == 0


def test_the_dynamics_of_run_stays_inside_the_enumerated_states():
    rng = np.random.default_rng(5)
    model = random_network(rng)
    keys = {cp._key(s) for s in cp.reachable(model)["states"]}
    sat = abs(model.saturation())
    for _ in range(10):
        run = hy.Run(model)
        for H in rng.uniform(-sat, sat, 10):
            run.ramp(float(H))
            assert cp._key(run.s) in keys


def test_a_relabeling_leaves_the_count_unchanged():
    rng = np.random.default_rng(9)
    model = random_network(rng, n=7)
    g = rng.choice([-1.0, 1.0], model.n)
    assert cp.reachable(model)["count"] == cp.reachable(hy.relabel(model, g))["count"]


def test_an_enumeration_that_stops_says_so():
    model = random_network(np.random.default_rng(4), n=10)
    assert cp.measure(model, max_states=3) == {"censored": True, "max_states": 3}


# ------------------------------------------------------------------------------------------------ counterpart
def test_a_balanced_network_is_its_own_counterpart():
    ferro = make(np.ones((4, 4)) - np.eye(4), [0.1, -0.2, 0.3, 0.0], [0.5] * 4, [1.0] * 4)
    fr = cp.frustration(ferro)
    assert fr["F"] == 0 and fr["exact"]
    assert cp.balanced_counterpart(ferro, fr["relabeling"]).bonds == ferro.bonds
    assert "its own counterpart" in cp.construct(ferro, np.random.default_rng(0))["statement"]


def test_an_antiferromagnetic_pair_in_a_uniform_drive_has_one_unsatisfiable_edge():
    pair = make([[0.0, -1.0], [-1.0, 0.0]], [0.0, 0.0], [0.5, 0.5], [1.0, 1.0])
    assert cp.frustration(pair)["F"] == 1


def test_the_counterpart_keeps_magnitudes_fields_and_thresholds_and_is_balanced():
    rng = np.random.default_rng(13)
    for _ in range(4):
        model = random_network(rng)
        fr = cp.frustration(model)
        twin = cp.balanced_counterpart(model, fr["relabeling"])
        assert hy.predict(twin)["guaranteed"] and cp.frustration(twin)["F"] == 0
        assert [(i, j, abs(J)) for i, j, J in twin.bonds] == [(i, j, abs(J)) for i, j, J in model.bonds]
        assert np.array_equal(twin.h, model.h) and np.array_equal(twin.b, model.b)
        assert np.array_equal(np.abs(twin.eta), np.abs(model.eta))
        reversed_edges = sum(J != K for (_, _, J), (_, _, K) in zip(model.undirected(), twin.undirected()))
        reversed_edges += int(np.sum(twin.eta != model.eta))
        assert reversed_edges == fr["F"]


def test_local_search_never_undercuts_the_exact_index():
    rng = np.random.default_rng(17)
    for _ in range(4):
        model = random_network(rng, n=9)
        exact = cp.frustration(model)
        bound = cp.frustration(model, np.random.default_rng(1), exact_up_to=0)
        assert exact["exact"] and not bound["exact"] and bound["F"] >= exact["F"]


def test_the_counterpart_returns_to_every_turning_point_and_the_parts_add_up():
    rng = np.random.default_rng(21)
    for _ in range(3):
        r = cp.construct(random_network(rng), np.random.default_rng(0))
        assert r["balanced"]["return"]["failed"] == 0
        assert r["network"]["C"] == pytest.approx(r["balanced"]["C"] + r["gain"])
        assert r["balanced"]["C"] == pytest.approx(r["uncoupled"]["C"] + r["loss"])


# ------------------------------------------------------------------------------------------------ example, command
def test_the_example_gains_states_over_a_counterpart_that_restores_the_return():
    spec = json.loads((EX / "frustrated_network.json").read_text(encoding="utf-8"))
    assert spec["question"] and spec["assumptions"] and spec["provenance"]["source"]
    r = cp.construct(hy.from_spec(spec), np.random.default_rng(0))
    assert r["frustration"]["F"] > 0 and r["frustration"]["exact"]
    assert r["network"]["return"]["failed"] > 0 and r["balanced"]["return"]["failed"] == 0
    assert r["gain"] > 0 and r["loss"] < 0


def test_command_line_writes_a_report_with_provenance(tmp_path):
    out = tmp_path / "counterpart"
    cmd = [sys.executable, "-B", "-m", "fieldbridge", "memory", "counterpart",
           str(EX / "frustrated_network.json"), "--out-dir", str(out)]
    subprocess.run(cmd, check=True, cwd=ROOT, capture_output=True)
    report = json.loads((out / "counterpart.json").read_text(encoding="utf-8"))
    assert report["command"] == "counterpart" and report["novelty_established"] is False
    assert report["balanced"]["return"]["failed"] == 0 and "gain" in report
    assert "Balanced counterpart" in (out / "counterpart.md").read_text(encoding="utf-8")


# ------------------------------------------------------------------------------------------------ the reference
def test_the_ferromagnetic_reference_is_unique_balanced_and_returns():
    rng = np.random.default_rng(5)
    m = random_network(rng, n=8)
    ref = cp.reference(m)
    assert all(J > 0 for _, _, J in ref.bonds) and (ref.eta > 0).all()
    assert cp.frustration(ref)["F"] == 0
    r = cp.measure(ref, return_test=True)
    assert r["return"]["failed"] == 0


def test_the_construction_reports_both_references_and_the_identities_hold():
    rng = np.random.default_rng(11)
    m = random_network(rng, n=8)
    r = cp.construct(m, rng)
    assert r["reference"]["unique"] and not r["counterpart"]["unique"]
    assert abs(r["network"]["C"] - (r["ferromagnetic_reference"]["C"] + r["gain"])) < 1e-12
    assert abs(r["ferromagnetic_reference"]["C"] - (r["uncoupled"]["C"] + r["loss"])) < 1e-12
    assert r["ferromagnetic_reference"]["return"]["failed"] == 0
    assert "typical_in_ensemble" in r and r["typical_in_ensemble"]["gain"] >= 0.0


def test_a_balanced_network_has_zero_gain_over_its_reference_up_to_a_relabeling():
    # a ferromagnet with a staggered drive relabeled on one sublattice is balanced; its states equal the reference's
    rng = np.random.default_rng(3)
    n = 8
    J = np.zeros((n, n))
    for i in range(n - 1):
        J[i, i + 1] = J[i + 1, i] = -0.5          # antiferromagnetic chain ...
    eta = np.array([1.0 if i % 2 == 0 else -1.0 for i in range(n)])   # ... with a staggered drive: balanced
    m = make(J, rng.normal(0, 0.7, n), rng.uniform(0.1, 0.8, n), eta)
    assert cp.frustration(m)["F"] == 0
    r = cp.construct(m, rng)
    assert r["network"]["count"] == r["balanced"]["count"]


def test_the_rules_registry_names_the_implemented_rule():
    assert cp.RULES["R8"]["implemented"] and cp.RULES["holonomy"]["implemented"]
    assert cp.RULES["write point"]["implemented"] == "construction" and cp.RULES["closure"]["implemented"]
    assert cp.RULES["closure"]["command"] == "quantum closure"


def test_predictions_command_writes_a_report(tmp_path):
    out = tmp_path / "pred"
    cmd = [sys.executable, "-B", "-m", "fieldbridge", "memory", "predictions", str(EX), "--out-dir", str(out),
           "--enumerate-up-to", "12"]
    res = subprocess.run(cmd, capture_output=True, text=True, cwd=ROOT)
    assert res.returncode == 0, res.stderr
    report = json.loads((out / "predictions.json").read_text(encoding="utf-8"))
    assert report["rule"] == "R8" and report["rows"]
    assert any("gain" in r for r in report["rows"])


def test_the_write_point_counterpart_holds_the_control_short_of_the_point():
    spec = json.loads((ROOT / "examples" / "memory" / "pitchfork.json").read_text(encoding="utf-8"))
    r = cp.write_point_check(spec, np.random.default_rng(0))
    assert r["write_points"], r
    assert r["counterpart"] is not None and r["verified"], r
    lo, hi = r["counterpart"]["control_range"]
    assert not (lo <= r["write_points"][0]["value"] <= hi)
