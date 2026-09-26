"""The language of mechanisms on quantum carriers: detachment, attachment and co-discovery of the Bloch rotation.

Expected values are analytic: rotation rates are the lengths of the field vectors, representations follow from the
carrier (N bosons: j = N/2; a chain of N spins with one flipped: j = (N-1)/2), and the couplings that carry a spin
rotation on a chain are proportional to sqrt((j + 1)(N - 1 - j)). The controls are terms that must break the rotation.
"""
import json
import math
import subprocess
import sys
from pathlib import Path

import pytest

np = pytest.importorskip("numpy")
pytest.importorskip("sympy")

from fieldbridge.quantum import carriers, language as ql

ROOT = Path(__file__).resolve().parents[1]
EX = ROOT / "examples" / "quantum"


def load(name):
    return ql.load(EX / f"{name}.json")


def test_carrier_operators_obey_their_algebras():
    J = carriers.spin_matrices(1.5)
    assert np.allclose(J["Jx"] @ J["Jy"] - J["Jy"] @ J["Jx"], 1j * J["Jz"])
    f = carriers.make({"kind": "fermions", "modes": 2}).tokens
    assert np.allclose(f["c0"] @ f["cd0"] + f["cd0"] @ f["c0"], np.eye(4))
    assert np.allclose(f["c0"] @ f["c1"] + f["c1"] @ f["c0"], 0)
    b = carriers.make({"kind": "bosons", "modes": ["a", "b"], "max_quanta": 3}).tokens
    comm = b["a"] @ b["ad"] - b["ad"] @ b["a"]
    low = np.array([i // 4 < 3 for i in range(16)])  # states below the truncation of mode a
    assert np.allclose(np.diag(comm)[low], 1.0)


def test_the_module_11_spins_carry_a_spin_made_of_correlations():
    row = ql.derive_bloch_rotation(load("two_spins"))
    assert row["word"] == "AKL" and row["status"] == "reached"
    a = next(s for s in row["steps"] if s["letter"] == "A")
    assert set(a["closing_operators"]) == {"X0", "Y0 Z1", "Z0 Z1"}
    assert row["representation"] == [{"j": 0.5, "copies": 2}]
    g, h = 1.0, 0.5
    assert row["signature"]["rate"] == pytest.approx(2 * math.hypot(g, h))  # H = 2g (Z0Z1/2) + 2h (X0/2)
    assert row["signature"]["theta_deg"] == pytest.approx(math.degrees(math.atan2(g, h)))
    assert row["law"]["residual"] < 1e-10
    plain = dict(load("two_spins").spec, parameters={"g": 1.0, "h": 0.0})  # Module 11 without the field
    row0 = ql.derive_bloch_rotation(ql.load(plain))
    assert row0["status"] == "reached" and row0["signature"]["theta_deg"] == pytest.approx(90.0)
    assert row0["law"]["value_at_inversion"] == pytest.approx(-1.0)  # x(t) = x(0) cos(2 g t): full inversion


@pytest.mark.parametrize("name, word, reps", [
    ("nmr_spin", "AKL", [{"j": 0.5, "copies": 1}]),
    ("bose_josephson", "SAKL", [{"j": 2.0, "copies": 1}]),
    ("state_transfer_chain", "SAKL", [{"j": 1.5, "copies": 1}]),
    ("cooper_pair", "AKL", [{"j": 0.5, "copies": 1}, {"j": 0.0, "copies": 2}]),
    ("spin1_atom", "AKL", [{"j": 1.0, "copies": 1}]),
])
def test_one_rotation_on_carriers_from_five_fields(name, word, reps):
    row = ql.derive_bloch_rotation(load(name))
    assert row["word"] == word and row["representation"] == reps
    assert row["law"]["residual"] < 1e-10


def test_obstructions_name_the_term_that_breaks_the_rotation():
    rows = {n: ql.derive_bloch_rotation(load(n)) for n in
            ("interacting_bosons", "nv_centre", "module14_chain", "heteronuclear_spins")}
    assert all(r["status"] == "obstructed" for r in rows.values())
    assert "U/2 [na na + nb nb]" in rows["interacting_bosons"]["obstruction"]  # one-axis twisting
    assert "D [Jz Jz]" in rows["nv_centre"]["obstruction"]  # zero-field splitting
    assert rows["module14_chain"]["algebra_dimension"] == 8  # su(3): unequal bonds 3 and 4
    assert rows["heteronuclear_spins"]["algebra_dimension"] == 6  # two rotations at different frequencies


def test_attaching_the_rotation_designs_the_transfer_chain():
    src = load("nmr_spin")
    a = ql.attach(src, "chain", 4)
    assert a["attached"]
    assert a["native"]["coupling_ratios"] == pytest.approx([1.0, 2 / math.sqrt(3), 1.0])  # sqrt(3) : 2 : sqrt(3)
    for key in ("rate", "weight", "theta_deg"):
        assert a["preserved"][key]["target"] == pytest.approx(a["preserved"][key]["source"])
    assert a["changed"]["representation"]["target"] == [{"j": 1.5, "copies": 1}]
    assert a["law"]["target_residual"] < 1e-9
    three = ql.attach(load("two_spins"), "chain", 3)
    assert three["native"]["coupling_ratios"] == pytest.approx([1.0, 1.0])  # Module 14's bonds would have to be equal


@pytest.mark.parametrize("to, size, rep", [("bosons", 4, [{"j": 2.0, "copies": 1}]),
                                           ("fermion-pair", None, [{"j": 0.5, "copies": 1}, {"j": 0.0, "copies": 2}]),
                                           ("correlated-pair", None, [{"j": 0.5, "copies": 2}]),
                                           ("spin", 3, [{"j": 1.5, "copies": 1}])])
def test_the_detached_rotation_attaches_to_every_carrier(to, size, rep):
    a = ql.attach(load("nmr_spin"), to, size)
    assert a["attached"] and a["changed"]["representation"]["target"] == rep
    assert a["law"]["target_residual"] < 1e-9


def test_an_obstructed_realization_has_nothing_to_detach():
    d = ql.detach(load("nv_centre"))
    assert d["detached"] is False and "su(3)" in d["reason"]


def test_quantum_commands_write_reports(tmp_path):
    run = lambda *a: subprocess.run([sys.executable, "-B", "-m", "fieldbridge", "quantum", *a],  # noqa: E731
                                    check=True, cwd=ROOT, capture_output=True)
    run("detach", str(EX / "two_spins.json"), "--out-dir", str(tmp_path / "d"))
    d = json.loads((tmp_path / "d" / "detach.json").read_text(encoding="utf-8"))
    assert d["detached"] and d["novelty_established"] is False and d["implementation_sha256"]
    run("attach", "--from", str(EX / "nmr_spin.json"), "--to", "bosons", "--size", "3", "--no-figure",
        "--out-dir", str(tmp_path / "a"))
    attached = ql.load(tmp_path / "a" / "attached.json")
    assert ql.derive_bloch_rotation(attached)["status"] == "reached"
    run("codiscover", str(EX / "nmr_spin.json"), str(EX / "nv_centre.json"), "--no-figure",
        "--out-dir", str(tmp_path / "c"))
    c = json.loads((tmp_path / "c" / "codiscover.json").read_text(encoding="utf-8"))
    assert [r["status"] for r in c["rows"]] == ["reached", "obstructed"]


@pytest.mark.parametrize("bad, message", [
    ({"schema": "other"}, "schema"),
    ({"carrier": {"kind": "qubits", "n": 9}}, "between 1 and 6"),
    ({"hamiltonian": [{"coefficient": 1, "operator": "Q0"}]}, "unknown operator"),
    ({"hamiltonian": [{"coefficient": "k", "operator": "X0"}]}, "undeclared"),
    ({"observable": [{"coefficient": 1, "operator": "X0 Z0"}]}, "not Hermitian"),
])
def test_specifications_outside_the_contract_are_refused(bad, message):
    spec = {**json.loads((EX / "nmr_spin.json").read_text(encoding="utf-8")), **bad}
    with pytest.raises(ql.SpecError, match=message):
        ql.load(spec)
