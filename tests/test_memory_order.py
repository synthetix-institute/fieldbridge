"""The order of writes measured against the uncoupled reference (fieldbridge.memory.order).

Controls: the ring of eight of the memory manuscript reproduces the composition census's control case state for
state; the uncoupled reference records the last write that moved each element and has no coupling; the commutation
index is zero when the writes act on separate components and counts the edits otherwise; the commands write reports.
"""
import json
import subprocess
import sys
from pathlib import Path

import pytest

np = pytest.importorskip("numpy")
pytest.importorskip("scipy")

from fieldbridge.memory import order as od
from fieldbridge.memory import counterpart as cp

ROOT = Path(__file__).resolve().parents[1]
EX = ROOT / "examples" / "counterpart" / "order"

CENSUS_RING = {"ABC": "--------", "ACB": "--------", "BAC": "-++++---", "BCA": "-+++++--", "CAB": "---+++--",
               "CBA": "--------"}


def spec(units, graph, couplings="attractive", coupling=0.13, writes=None):
    writes = writes or {"A": {"k": units - 1}, "B": {"k": units - 1}}
    return {"schema": "fieldbridge-memory/1", "kind": "order", "name": "t",
            "order": {"units": units, "graph": graph, "couplings": couplings, "coupling": coupling, "writes": writes}}


def test_the_ring_of_eight_reproduces_the_census_control():
    net = od.from_spec(json.loads((EX / "ring8.json").read_text(encoding="utf-8")))
    m = od.measure(net)
    assert m["states"] == CENSUS_RING
    assert m["distinct"] == 4 and m["split"] == 3


def test_the_uncoupled_reference_has_no_coupling_and_the_identity_holds():
    net = od.from_spec(spec(4, "ring"))
    ref = od.uncoupled(net)
    assert ref.coupling == 0.0 and ref.edges == net.edges
    r = od.construct(net)
    assert abs(r["gain"] - (r["network"]["capacity"] - r["uncoupled_reference"]["capacity"])) < 1e-12


def test_the_commutation_index_is_zero_for_writes_on_separate_components_and_counts_edits():
    # two nodes, no edge, one write each: the writes commute
    s = spec(2, [], writes={"A": {"support": [0]}, "B": {"support": [1]}})
    assert od.commutation_index(od.from_spec(s))["H"] == 0
    # the same with one coupling: one cut, or one write removed
    s = spec(2, [[0, 1]], writes={"A": {"support": [0]}, "B": {"support": [1]}})
    assert od.commutation_index(od.from_spec(s))["H"] == 1
    # both writes on both nodes of an edge: each node loses one write, and the edge is cut or both nodes take one write
    s = spec(2, [[0, 1]], writes={"A": {"support": [0, 1]}, "B": {"support": [0, 1]}})
    assert od.commutation_index(od.from_spec(s))["H"] == 2


def test_a_specification_error_names_the_field():
    bad = spec(3, "ring"); bad["order"]["couplings"] = "random"
    with pytest.raises(Exception) as e:
        od.from_spec(bad)
    assert "couplings" in str(e.value)


def test_the_registry_marks_the_order_rule_implemented():
    assert cp.RULES["holonomy"]["implemented"] and cp.RULES["holonomy"]["module"] == "fieldbridge.memory.order"


def test_order_and_predictions_commands_write_reports(tmp_path):
    out = tmp_path / "order"
    res = subprocess.run([sys.executable, "-B", "-m", "fieldbridge", "memory", "order", str(EX / "ring8.json"),
                          "--out-dir", str(out)], capture_output=True, text=True, cwd=ROOT)
    assert res.returncode == 0, res.stderr
    report = json.loads((out / "order.json").read_text(encoding="utf-8"))
    assert report["network"]["distinct"] == 4 and report["network"]["split"] == 3
    out2 = tmp_path / "pred"
    res = subprocess.run([sys.executable, "-B", "-m", "fieldbridge", "memory", "predictions", str(EX), "--out-dir",
                          str(out2), "--enumerate-up-to", "8"], capture_output=True, text=True, cwd=ROOT)
    assert res.returncode == 0, res.stderr
    rows = json.loads((out2 / "predictions.json").read_text(encoding="utf-8"))["rows"]
    assert all(r["rule"] == "holonomy" for r in rows) and any("gain" in r for r in rows)
