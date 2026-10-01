"""Return to a turning point: interacting hysterons under a slow drive.

The structural prediction (no frustrated loop in the network that counts the drive as an element) is compared with
subloops of the quasi-static dynamics. The controls are changes that must change the answer: reversing the couplings
in a uniform field, staggering the field, and reversing the response of one element to the drive with and without
its couplings.
"""
import json
import subprocess
import sys
from dataclasses import replace
from pathlib import Path

import pytest

np = pytest.importorskip("numpy")
pytest.importorskip("scipy")

from fieldbridge.memory import hysterons as hy
from fieldbridge.memory.spec import SpecError

ROOT = Path(__file__).resolve().parents[1]
EX = ROOT / "examples" / "memory" / "hysterons"
SCHEMA = "fieldbridge-memory/1"


def lattice(shape="square", L=8, p=0.0, drive="uniform", width=2.0, b=0.0, seed=1):
    return hy.from_spec({"schema": SCHEMA, "kind": "hysterons", "name": "test",
                         "hysterons": {"lattice": {"shape": shape, "L": L},
                                       "couplings": {"value": 1.0, "negative_fraction": p},
                                       "fields": {"distribution": "gaussian", "width": width}, "half_widths": b,
                                       "drive": drive, "seed": seed}})


def gauged_graph(rng, n=30, extra=30):
    """A random graph that becomes cooperative, with a uniform drive, after the relabeling g."""
    g = rng.choice([-1.0, 1.0], n)
    edges = {(int(rng.integers(i)), i) for i in range(1, n)}            # a spanning tree
    while len(edges) < n - 1 + extra:
        i, j = sorted(int(x) for x in rng.choice(n, 2, replace=False))
        edges.add((i, j))
    couplings = [[i, j, float(g[i] * g[j] * rng.uniform(0.2, 1.5))] for i, j in sorted(edges)]
    spec = {"schema": SCHEMA, "kind": "hysterons", "name": "gauged",
            "hysterons": {"units": n, "couplings": couplings, "fields": {"distribution": "gaussian", "width": 1.5},
                          "half_widths": {"distribution": "uniform", "low": 0.0, "high": 0.5},
                          "drive": [float(x * rng.uniform(0.5, 1.5)) for x in g], "seed": int(rng.integers(1000))}}
    return hy.from_spec(spec), g


# ------------------------------------------------------------------------------------------------ structure
def test_a_network_cooperative_after_relabeling_returns_to_every_turning_point():
    rng = np.random.default_rng(7)
    for _ in range(3):
        model, g = gauged_graph(rng)
        pred = hy.predict(model)
        assert pred["guaranteed"] and pred["relabeling"] == [int(x) for x in g]
        assert pred["canonical"]["min_coupling"] > 0 and pred["canonical"]["drive_positive"]
        result = hy.check(model, rng, count=6)
        assert result["consistent"] and result["failed"] == 0 and result["did_not_end"] == 0


def test_the_relabeling_maps_the_dynamics_onto_the_cooperative_form():
    model, g = gauged_graph(np.random.default_rng(11))
    canonical = hy.relabel(model, g)
    assert all(J >= 0 for _, _, J in canonical.bonds) and np.all(canonical.eta > 0)
    a, b = hy.Run(model), hy.Run(canonical)
    for H in (-0.5, 1.2, -0.8, 0.4, -0.2, 1.2, 3.0, -3.0):
        a.ramp(H)
        b.ramp(H)
        assert np.array_equal(b.s, g * a.s)


def test_an_antiferromagnet_frustrates_the_loops_through_the_drive_not_its_plaquettes():
    anti = hy.predict(lattice(p=1.0))
    assert anti["couplings_balanced"] and anti["frustrated_plaquettes"] == 0.0
    assert not anti["guaranteed"] and anti["frustrated_with_drive"] == 1.0
    # the staggered field is the same network relabeled; the ferromagnet in a uniform field is its canonical form
    assert hy.predict(lattice(p=1.0, drive="staggered"))["guaranteed"]
    assert hy.predict(lattice())["guaranteed"]
    # positive couplings under a drive of random sign: the couplings alone are not frustrated, the drive is
    mixed = hy.predict(lattice(drive="random"))
    assert mixed["couplings_balanced"] and not mixed["guaranteed"]


def test_a_frustrated_loop_through_the_drive_can_break_the_return_and_need_not():
    rng = np.random.default_rng(3)
    assert hy.check(lattice(L=16, p=1.0), rng, count=8)["failed"] > 0
    assert hy.check(lattice(L=16, p=1.0, drive="staggered"), rng, count=8)["failed"] == 0
    # random antiferromagnetic chains started from a large field return exactly (Deutsch, Dhar and Narayan, 2004)
    chain = lattice(shape="chain", L=128, p=1.0)
    result = hy.check(chain, rng, count=8)
    assert not result["prediction"]["guaranteed"] and result["failed"] == 0
    assert "allows a failure without forcing one" in result["verdict"]


def test_reversing_one_element_with_all_its_couplings_is_a_relabeling():
    model = lattice()
    flipped = model.eta.copy()
    flipped[0] = -1.0
    only_drive = replace(model, eta=flipped)
    pred = hy.predict(only_drive)
    assert not pred["guaranteed"] and pred["bonds_frustrated_with_drive"] == 4
    bonds = [(i, j, -J if 0 in (i, j) else J) for i, j, J in model.bonds]
    assert hy.predict(replace(model, eta=flipped, bonds=bonds))["guaranteed"]


def test_independent_hysterons_are_the_preisach_model():
    model = hy.from_spec({"schema": SCHEMA, "kind": "hysterons", "name": "pores",
                          "hysterons": {"units": 60, "couplings": [], "fields": {"distribution": "gaussian", "width": 1},
                                        "half_widths": {"distribution": "uniform", "low": 0.1, "high": 0.8},
                                        "seed": 5}})
    pred = hy.predict(model)
    assert pred["guaranteed"] and pred["canonical"]["form"] == "independent hysterons (Preisach)"
    assert hy.check(model, np.random.default_rng(1), count=4)["failed"] == 0


def test_an_avalanche_that_does_not_end_is_reported():
    # element 0 drives 1 up, element 1 drives 0 down: at fixed field the pair switches in a cycle
    model = hy.from_spec({"schema": SCHEMA, "kind": "hysterons", "name": "cycle",
                          "hysterons": {"units": 2, "directed": True,
                                        "couplings": [{"from": 0, "to": 1, "J": 2.0}, {"from": 1, "to": 0, "J": -2.0}]}})
    assert not model.reciprocal and not hy.predict(model)["guaranteed"]
    with pytest.raises(hy.AvalancheError):
        hy.Run(model).ramp(0.0)
    result = hy.check(model, np.random.default_rng(0), count=2)
    assert result["did_not_end"] and result["consistent"] and "did not end" in result["verdict"]


def test_directed_cooperative_couplings_still_return():
    rng = np.random.default_rng(2)
    n = 20
    couplings = [{"from": int(i), "to": int(j), "J": float(rng.uniform(0.2, 1.0))}
                 for i, j in {tuple(sorted(rng.choice(n, 2, replace=False))) for _ in range(40)}]
    model = hy.from_spec({"schema": SCHEMA, "kind": "hysterons", "name": "directed",
                          "hysterons": {"units": n, "directed": True, "couplings": couplings,
                                        "fields": {"distribution": "gaussian", "width": 1.0}, "seed": 3}})
    assert hy.predict(model)["guaranteed"] and hy.check(model, rng, count=4)["failed"] == 0


@pytest.mark.parametrize("bad, message", [
    ({"kind": "network"}, "kind must be"),
    ({"hysterons": None}, "needs a 'hysterons' object"),
    ({"hysterons": {"lattice": {"shape": "hexagonal", "L": 8}}}, "lattice.shape"),
    ({"hysterons": {"lattice": {"shape": "square", "L": 2}}}, "at least 3"),
    ({"hysterons": {"lattice": {"shape": "square", "L": 4}, "couplings": {"negative_fraction": 2}}}, "between 0 and 1"),
    ({"hysterons": {"units": 3, "drive": "staggered"}}, "staggered drive needs a lattice"),
    ({"hysterons": {"units": 3, "half_widths": -1}}, "must not be negative"),
    ({"hysterons": {"units": 3, "couplings": [[1, 1, 0.5]]}}, "does not act on itself"),
    ({"hysterons": {"units": 3, "couplings": [[0, 3, 0.5]]}}, "numbered 0 to 2"),
    ({"hysterons": {"units": 3, "fields": [1, 2]}}, "expected 3 values"),
    ({"hysterons": {"units": 3, "fields": {"distribution": "cauchy"}}}, "distribution must be"),
    ({"hysterons": {}}, "Give either 'lattice'"),
])
def test_specifications_outside_the_contract_are_refused(bad, message):
    spec = {"schema": SCHEMA, "kind": "hysterons", "name": "bad", "hysterons": {"units": 3}}
    spec.update(bad)
    with pytest.raises(SpecError, match=message):
        hy.from_spec(spec)


def test_a_drive_that_acts_on_no_element_is_refused():
    model = hy.from_spec({"schema": SCHEMA, "kind": "hysterons", "name": "undriven",
                          "hysterons": {"units": 2, "couplings": [[0, 1, 1.0]], "drive": [0, 0]}})
    with pytest.raises(SpecError, match="No element is driven"):
        hy.Run(model)


# ------------------------------------------------------------------------------------------------ examples, command
EXPECTED = {"rfim_ferromagnet": True, "rfim_antiferromagnet": False, "rfim_antiferromagnet_staggered": True,
            "antiferromagnetic_chain": False, "adsorption_pores": True, "soft_spots": False}


def test_the_examples_state_their_question_and_are_predicted_from_structure():
    found = {p.stem: json.loads(p.read_text(encoding="utf-8")) for p in sorted(EX.glob("*.json"))}
    assert set(found) == set(EXPECTED)
    fields = set()
    for name, spec in found.items():
        assert spec["question"] and spec["assumptions"] and spec["provenance"]["source"] and spec["id"] == name
        fields.add(spec["field"])
        assert hy.predict(hy.from_spec(spec))["guaranteed"] is EXPECTED[name], name
    assert len(fields) == 3


def test_command_line_writes_a_report_with_provenance(tmp_path):
    out = tmp_path / "staggered"
    cmd = [sys.executable, "-B", "-m", "fieldbridge", "memory", "hysterons",
           str(EX / "rfim_antiferromagnet_staggered.json"), "--out-dir", str(out), "--subloops", "2"]
    subprocess.run(cmd, check=True, cwd=ROOT, capture_output=True)
    report = json.loads((out / "hysterons.json").read_text(encoding="utf-8"))
    assert report["command"] == "hysterons" and report["novelty_established"] is False
    assert report["prediction"]["guaranteed"] and report["failed"] == 0 and report["total"] == 8
    assert "Return to a turning point" in (out / "hysterons.md").read_text(encoding="utf-8")
