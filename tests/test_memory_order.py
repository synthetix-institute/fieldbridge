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


# ------------------------------------------------------------------------------------------------ equation networks
from fieldbridge.memory import order_equations as oe


def amplitude_network(n, edges, J, noise=0.002):
    """Parametric amplitudes x_i' = a_i x_i - x_i^3 + h_i + J sum_j x_j (the coherent Ising machine's normal form)."""
    drift = {f"x{i}": f"a{i}*x{i} - x{i}**3 + h{i}" + "".join(f" + J*x{j}" for a, b in edges for j in
                                                               ((b,) if a == i else (a,) if b == i else ())) for i in range(n)}
    params = {"J": J, **{f"a{i}": -1.0 for i in range(n)}, **{f"h{i}": 0.0 for i in range(n)}}
    return {"schema": "fieldbridge-memory/1", "kind": "equations", "name": "amplitudes", "question": "q",
            "assumptions": [], "provenance": {"source": "test"},
            "carrier": {"kind": "euclid", "variables": [f"x{i}" for i in range(n)], "scale": 1},
            "parameters": params, "drift": drift, "noise": noise, "dt": 0.01}


def equation_spec(network, n, writes, mode="ramp", couplings="attractive", **write):
    w = {"mode": mode, "rest": -1.0, "on": 1.0, "ramp": 8.0, "bias": 0.05, "hold": 2.0, "settle": 2.0, "release": 20.0}
    w.update(write)
    return {"schema": "fieldbridge-memory/1", "kind": "order", "name": "t",
            "order": {"material": "equations", "network": network,
                      "nodes": {"count": n, "variables": ["x{i}"], "readout": "x{i}", "control": "a{i}", "bias": "h{i}"},
                      "coupling": ["J"], "couplings": couplings, "write": w, "writes": writes}}


def test_sequential_ramps_record_the_order_only_where_a_coupling_is_against_the_write_signs():
    # attractive pair written with opposite signs: the second node follows the first (two states, reference one)
    net = od.from_spec(equation_spec(amplitude_network(2, [(0, 1)], 0.2), 2,
                                     {"A": {"support": [0], "sign": 1}, "B": {"support": [1], "sign": -1}}))
    assert net.edges == [(0, 1)]
    r = od.construct(net)
    assert r["network"]["states"] == {"AB": "++", "BA": "--"} and r["uncoupled_reference"]["distinct"] == 1
    assert r["index"]["H"] == 1 and abs(r["gain"] - 1.0) < 1e-12
    # the same signs on an antiferromagnetic pair satisfy the coupling: the writes commute although H = 1
    net = od.from_spec(equation_spec(amplitude_network(2, [(0, 1)], -0.2), 2,
                                     {"A": {"support": [0], "sign": 1}, "B": {"support": [1], "sign": -1}},
                                     couplings="antiferromagnetic"))
    r = od.construct(net)
    assert r["network"]["distinct"] == 1 and r["index"]["H"] == 1 and r["gain"] == 0.0
    # a coupling weaker than the write bias is not recorded
    net = od.from_spec(equation_spec(amplitude_network(2, [(0, 1)], 0.01), 2,
                                     {"A": {"support": [0], "sign": 1}, "B": {"support": [1], "sign": -1}}))
    assert od.measure(net)["distinct"] == 1


def test_the_graph_of_an_equation_network_is_read_from_its_drift_and_unwritten_nodes_read_zero():
    net = od.from_spec(equation_spec(amplitude_network(3, [(0, 1), (1, 2)], 0.2), 3,
                                     {"A": {"support": [0], "sign": 1}, "B": {"support": [2], "sign": -1}}))
    assert net.edges == [(0, 1), (1, 2)]
    m = od.measure(net)
    assert all(s[1] == "0" for s in m["states"].values())     # node 1 is never written and holds no state
    assert od.commutation_index(net)["H"] == 1


def test_pulse_writes_on_bistable_nodes_and_a_specification_error_names_the_field():
    # pulses far above the coercive field on separate nodes of a weakly coupled pair: each node keeps its last pulse
    s = equation_spec(amplitude_network(2, [(0, 1)], 0.02), 2,
                      {"A": {"support": [0, 1], "sign": 1}, "B": {"support": [1], "sign": -1}}, mode="pulse",
                      bias=1.5, pulse=4.0)
    m = od.measure(od.from_spec(s))
    assert m["states"] == {"AB": "+-", "BA": "++"} and m["split"] == 0
    bad = equation_spec(amplitude_network(2, [(0, 1)], 0.2), 2, {"A": {"support": [0]}, "B": {"support": [1]}})
    bad["order"]["nodes"]["control"] = "pump{i}"
    with pytest.raises(Exception) as e:
        od.from_spec(bad)
    assert "control" in str(e.value)


def test_the_noisy_record_is_the_capacity_of_the_order_to_pattern_channel():
    assert abs(oe.capacity(np.eye(3)) - np.log2(3)) < 1e-9
    assert oe.capacity(np.full((4, 2), 0.5)) < 1e-9
    assert abs(oe.capacity(np.array([[0.9, 0.1], [0.1, 0.9]])) - (1 + 0.9 * np.log2(0.9) + 0.1 * np.log2(0.1))) < 1e-6
    s = equation_spec(amplitude_network(2, [(0, 1)], 0.2), 2,
                      {"A": {"support": [0], "sign": 1}, "B": {"support": [1], "sign": -1}})
    s["order"]["noise"] = {"trajectories": 60, "seed": 3}
    r = od.construct(od.from_spec(s))
    a = r["noisy"]["network"]
    assert a["capacity"] > 0.8 and a["p_value"] < 0.05 and r["noisy"]["uncoupled_reference"]["capacity"] < 0.2
