"""The language of mechanisms on quantum carriers: detachment, attachment and co-discovery of the Bloch rotation.

Expected values are analytic: rotation rates are the lengths of the field vectors, representations follow from the
carrier (N bosons: j = N/2; a chain of N spins with one flipped: j = (N-1)/2), and the couplings that carry a spin
rotation on a chain are proportional to sqrt((j + 1)(N - 1 - j)). The controls are terms that must break the rotation.

The rotation is reached through the algebra (the Hamiltonian and the observable generate su(2)) or through the
closure of the observable (the algebra is larger, the commutators of the Hamiltonian with the observable close on
three operators or fewer, and the prepared state lies along the observable). The controls of the second route are a
closure of five operators, and closures of three operators in which the prepared state does not lie along the
observable: the signal then has one frequency and another amplitude.
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


def spins(hamiltonian, observable=((1, "X0"),), sector=None, base="two_spins", name=None):
    """A realization on the carrier of an example, with the Hamiltonian and the observable given as pairs
    (coefficient, operator); an operator ending in ' + h.c.' is completed by its Hermitian conjugate. The parameters
    are those of the example (two spins: g = 1, h = 0.5)."""
    def terms(pairs):
        return [{"coefficient": c, "operator": op.removesuffix(" + h.c."), **({"hc": True} if op.endswith(" + h.c.") else {})}
                for c, op in pairs]

    spec = {**json.loads((EX / f"{base}.json").read_text(encoding="utf-8")), "hamiltonian": terms(hamiltonian),
            "observable": terms(observable)}
    if name:
        spec["name"] = name
    if sector:
        spec["sector"] = sector
    return ql.load(spec)


MOVED = (("g", "Z0 Z1"), ("h", "X1"))                     # the field of the two spins on the second spin
BOTH = (("g", "Z0 Z1"), ("h", "X0"), ("h", "X1"))          # a field on both spins
COUNTER = dict(hamiltonian=((1, "Z0"), (1, "Z1")), observable=((1, "X0"), (0.5, "X0 Z1"), (0.7, "Z0")))


def exact_amplitude(real, psi0=None):
    """b of the exact signal f(t) = 1 - b (1 - cos |Omega| t) from a top eigenstate of the observable (the one eigh
    returns, or psi0), and the largest deviation of the signal from this form, by direct evolution."""
    O = real.O - np.trace(real.O) / real.O.shape[0] * np.eye(real.O.shape[0])
    w, U = np.linalg.eigh(O)
    E, V = np.linalg.eigh(real.H)
    c = V.conj().T @ (U[:, -1] if psi0 is None else np.asarray(psi0, complex) / np.linalg.norm(psi0))
    rate = ql.observable_frequencies(real.H, O)[0]
    t = np.linspace(0.0, 4 * math.pi / rate, 241)
    f = np.array([np.real((V @ (np.exp(-1j * E * tk) * c)).conj() @ O @ (V @ (np.exp(-1j * E * tk) * c))) for tk in t]) / w[-1]
    b = (1.0 - f.min()) / 2
    return b, float(np.max(np.abs(f - (1 - b * (1 - np.cos(rate * t))))))


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


def test_the_rotation_is_reached_through_the_closure_of_the_observable():
    """H = g Z0 Z1 + h X1 with the observable X0: the algebra is so(4), of dimension 6, but i[H, .] carries X0,
    Y0 Z1 and Y0 Y1 into one another. The rotation has the rate, the weight and the angle of the two spins."""
    two = ql.derive_bloch_rotation(load("two_spins"))
    row = ql.derive_bloch_rotation(spins(MOVED))
    assert row["word"] == "AOKL" and row["status"] == ql.REACHED_CLOSURE == "reached through the closure of the observable"
    assert ql.reached(row) and ql.reached(two) and two["status"] == ql.REACHED
    assert (row["algebra_dimension"], row["closure_dimension"]) == (6, 3)
    assert row["frequencies"] == pytest.approx([2 * math.hypot(1.0, 0.5)])
    step = {s["letter"]: s for s in row["steps"]}
    assert set(step["O"]["closure_operators"]) == {"X0", "Y0 Z1", "Y0 Y1"}
    assert set(step["A"]["closing_operators"]) == {"X0", "X1", "Y0 Y1", "Y0 Z1", "Z0 Y1", "Z0 Z1"}
    for key in ("rate", "weight", "theta_deg"):
        assert row["signature"][key] == pytest.approx(two["signature"][key], rel=1e-12)
    assert "algebra" not in row["signature"] and row["signature"]["relations"] == "i[H, J_a] = eps_abc Omega_b J_c"
    assert step["K"]["closes"] is False and row["representation"] is None   # the three operators are not an algebra
    assert step["K"]["rotation_residual"] < 1e-12 and abs(step["K"]["amplitude_offset"]) < 1e-12
    assert row["law"]["residual"] < 1e-10
    assert "without the term h [X1] the algebra is su(2)" in row["algebra"]
    # a conserved term beside the su(2): dimension 4, and the closure is the su(2) of the two spins itself
    kept = ql.derive_bloch_rotation(spins((("g", "Z0 Z1"), ("h", "X0"), ("h", "Z1"))))
    assert kept["word"] == "AOKL" and kept["algebra_dimension"] == 4
    assert kept["representation"] == [{"j": 0.5, "copies": 2}] and kept["law"]["residual"] < 1e-10
    assert kept["signature"]["rate"] == pytest.approx(two["signature"]["rate"])


def test_the_controls_of_the_closure_give_the_results_of_the_tutorial():
    """The three specifications of examples/quantum/controls, with the table of Section 6 of the chapter."""
    rows = {p.stem: ql.derive_bloch_rotation(ql.load(p)) for p in sorted((EX / "controls").glob("*.json"))}
    assert set(rows) == {"field_on_second_spin", "field_on_both_spins", "weighted_observable"}
    moved, both, weighted = rows["field_on_second_spin"], rows["field_on_both_spins"], rows["weighted_observable"]
    assert (moved["word"], moved["algebra_dimension"], moved["closure_dimension"]) == ("AOKL", 6, 3)
    assert moved["signature"]["rate"] == pytest.approx(math.sqrt(5)) and moved["signature"]["weight"] == pytest.approx(2.0)
    assert moved["signature"]["theta_deg"] == pytest.approx(math.degrees(math.atan2(1.0, 0.5)))
    assert moved["law"]["residual"] < 1e-12
    assert (both["word"], both["algebra_dimension"], both["closure_dimension"]) == ("AO", 6, 5)
    assert both["frequencies"] == pytest.approx([2.0, 2 * math.sqrt(2)])
    assert (weighted["word"], weighted["algebra_dimension"], weighted["closure_dimension"]) == ("AO", 7, 3)
    assert weighted["frequencies"] == pytest.approx([2.0])
    assert "the amplitude of the signal is 0.8212, not sin^2 theta = 0.7184" in weighted["obstruction"]
    # they are the models built in these tests from the two spins
    for name, (ham, obs) in {"field_on_second_spin": (MOVED, ((1, "X0"),)), "field_on_both_spins": (BOTH, ((1, "X0"),)),
                             "weighted_observable": (COUNTER["hamiltonian"], COUNTER["observable"])}.items():
        file, built = ql.load(EX / "controls" / f"{name}.json"), spins(ham, obs)
        assert np.allclose(file.H, built.H) and np.allclose(file.O, built.O), name


def test_a_closure_of_five_operators_has_two_frequencies_and_stops_the_derivation():
    row = ql.derive_bloch_rotation(spins(BOTH))
    assert row["word"] == "AO" and row["status"] == "obstructed" and not ql.reached(row)
    assert (row["algebra_dimension"], row["closure_dimension"]) == (6, 5)
    assert row["frequencies"] == pytest.approx([2.0, 2 * math.sqrt(2)])     # 2g and 2 sqrt(g^2 + 4 h^2)
    assert row["obstruction_short"] == "closure of 5 operators, 2 frequencies"
    assert "5 operators, which move with 2 frequencies" in row["obstruction"]
    assert "without the term h [X1] the algebra is su(2)" in row["obstruction"]


def test_one_frequency_does_not_fix_the_amplitude_of_the_signal():
    """H = Z0 + Z1 with R = X0 (1 + 0.5 Z1) + 0.7 Z0. The closure of R is R, Y0 (1 + 0.5 Z1) and Z0: three operators
    and the frequency 2, with sin^2 theta = 5 / 6.96 from their norms. The top eigenstate of R lies in the sector
    Z1 = +1, where R = 1.5 X0 + 0.7 Z0 and the amplitude is 2.25 / 2.74. The rotation is not reached; with the
    sector declared it is, through the algebra."""
    law, signal = 5 / 6.96, 2.25 / 2.74
    real = spins(**COUNTER)
    row = ql.derive_bloch_rotation(real)
    assert row["word"] == "AO" and row["status"] == "obstructed"
    assert row["closure_dimension"] == 3 and row["frequencies"] == pytest.approx([2.0])
    assert row["obstruction_short"] == "one frequency, amplitude not sin^2 theta"
    assert f"the amplitude of the signal is {signal:.4g}, not sin^2 theta = {law:.4g}" in row["obstruction"]
    O0 = ql._traceless(real.O)
    canon = ql.canonical_closure(ql.closure_basis(real.H, O0), real.H, O0)
    assert math.sin(canon["theta"]) ** 2 == pytest.approx(law) and canon["amplitude_offset"] == pytest.approx(signal - law)
    b, deviation = exact_amplitude(real)
    assert b == pytest.approx(signal) and deviation < 1e-12          # one frequency, with the other amplitude
    sector = ql.derive_bloch_rotation(spins(**COUNTER, sector={"operator": [{"coefficient": 1, "operator": "Z1"}], "value": 1}))
    assert sector["word"] == "SAKL" and sector["status"] == "reached"
    assert math.sin(math.radians(sector["signature"]["theta_deg"])) ** 2 == pytest.approx(signal)
    # a spin 1 with a quadrupole part in the observable: again three operators, one frequency, another amplitude
    spin1 = spins(((1, "Jz"),), ((0.7, "Jz"), (1, "Jx"), (0.4, "Jx Jz + h.c.")), base="spin1_atom")
    row1 = ql.derive_bloch_rotation(spin1)
    assert row1["word"] == "AO" and row1["closure_dimension"] == 3 and row1["frequencies"] == pytest.approx([1.0])
    O1 = ql._traceless(spin1.O)
    canon1 = ql.canonical_closure(ql.closure_basis(spin1.H, O1), spin1.H, O1)
    b1, deviation1 = exact_amplitude(spin1)
    assert b1 == pytest.approx(math.sin(canon1["theta"]) ** 2 + canon1["amplitude_offset"]) and deviation1 < 1e-12
    assert abs(canon1["amplitude_offset"]) > 5e-3


def test_a_degenerate_top_eigenspace_is_tested_as_a_whole():
    """H = Z0 + Z1 with R = X0 X1: the closure has three operators and the frequency 4. The product state |++> gives
    the amplitude 1/2 = sin^2 theta, the state (|00> + |11>)/sqrt 2 of the same eigenvalue gives 1: the preparation
    'top eigenstate' does not fix the signal, and the rotation is not reached."""
    real = spins(((1, "Z0"), (1, "Z1")), ((1, "X0 X1"),))
    row = ql.derive_bloch_rotation(real)
    assert row["word"] == "AO" and row["closure_dimension"] == 3 and row["frequencies"] == pytest.approx([4.0])
    assert row["obstruction_short"] == "one frequency, amplitude depends on the prepared state"
    assert "the amplitude of the signal is between 0 and 1, depending on the state prepared, and sin^2 theta = 0.5" \
        in row["obstruction"]
    canon = ql.canonical_closure(ql.closure_basis(real.H, real.O), real.H, real.O)
    assert canon["top_dimension"] == 2 and canon["amplitude_range"] == pytest.approx([0.0, 1.0], abs=1e-12)
    assert math.sin(canon["theta"]) ** 2 == pytest.approx(0.5) and abs(canon["amplitude_offset"]) == pytest.approx(0.5)
    # exact evolution from three states with X0 X1 = +1, in the basis |00>, |01>, |10>, |11>
    for psi0, b in (([1, 1, 1, 1], 0.5), ([1, 0, 0, 1], 1.0), ([0, 1, 1, 0], 0.0)):
        got, deviation = exact_amplitude(real, psi0)
        assert got == pytest.approx(b, abs=1e-12) and deviation < 1e-12
    # the two spins with the field on the second spin also have two top eigenstates, and both lie along the observable
    moved = spins(MOVED)
    ok = ql.canonical_closure(ql.closure_basis(moved.H, moved.O), moved.H, moved.O)
    assert ok["top_dimension"] == 2 and ok["amplitude_range"] == pytest.approx([0.8, 0.8], abs=1e-12)


def test_a_closure_of_two_operators_is_a_rotation_at_ninety_degrees():
    """H = g Z0 Z1 + h Z1: X0 and Y0 Z1 rotate into one another at 2g, Z1 is conserved, and the algebra has dimension 4."""
    row = ql.derive_bloch_rotation(spins((("g", "Z0 Z1"), ("h", "Z1"))))
    assert row["word"] == "AOKL" and (row["algebra_dimension"], row["closure_dimension"]) == (4, 2)
    assert row["signature"]["rate"] == pytest.approx(2.0) and row["signature"]["theta_deg"] == pytest.approx(90.0)
    assert row["law"]["residual"] < 1e-10 and row["law"]["value_at_inversion"] == pytest.approx(-1.0)


def test_a_conserved_observable_is_named_whatever_the_algebra():
    row = ql.derive_bloch_rotation(spins((("g", "Z0 Z1"), ("h", "X0")), ((1, "Z1"),)))
    assert row["word"] == "AO" and row["status"] == "obstructed" and row["closure_dimension"] == 1
    assert row["algebra_dimension"] == 2 and row["obstruction_short"] == "conserved observable"
    assert row["obstruction"].startswith("the observable commutes with the Hamiltonian")


@pytest.mark.parametrize("name", ["nmr_spin", "two_spins", "bose_josephson", "state_transfer_chain", "cooper_pair",
                                  "spin1_atom"])
def test_the_frame_of_the_closure_is_the_frame_of_the_algebra_where_both_exist(name):
    """In the realizations that reach the rotation through su(2), the closure of the observable gives the same
    rate, weight and angle, and the prepared state lies along the observable."""
    real = load(name)
    V, _ = ql.restrict(real)
    Hs, Os = V.conj().T @ real.H @ V, V.conj().T @ real.O @ V
    H0, O0 = ql._traceless(Hs), ql._traceless(Os)
    algebra = ql.canonical_su2(ql.lie_closure([Hs, Os])[0], H0, O0)
    closure = ql.canonical_closure(ql.closure_basis(Hs, O0), Hs, O0)
    assert closure["rate"] == pytest.approx(algebra["rate"], rel=1e-10)
    assert closure["weight"] == pytest.approx(algebra["weight"], rel=1e-10)
    assert math.sin(closure["theta"]) ** 2 == pytest.approx(math.sin(algebra["theta"]) ** 2, abs=1e-10)
    assert abs(closure["amplitude_offset"]) < 1e-12 and closure["rotation_residual"] < 1e-10


def test_a_rotation_reached_through_the_closure_detaches_and_attaches():
    moved = spins(MOVED, name="two coupled spins, field on the second spin")
    d = ql.detach(moved)
    assert d["detached"] and d["derivation"]["word"] == "AOKL"
    left = d["left_behind"]
    assert set(left["closure_operators"]) == {"X0", "Y0 Z1", "Y0 Y1"} and left["algebra_dimension"] == 6
    assert left["representation"] is None and "closing_operators" not in left
    for to, size in (("qubit", None), ("chain", 4), ("bosons", 4), ("correlated-pair", None)):
        a = ql.attach(moved, to, size)
        assert a["attached"] and a["target_derivation"]["status"] == "reached", to   # on the new carrier: through su(2)
        for key in ("rate", "weight", "theta_deg"):
            assert a["preserved"][key]["target"] == pytest.approx(a["preserved"][key]["source"]), (to, key)
        assert a["law"]["target_residual"] < 1e-9
    # written on the correlations of two spins, the rotation is the model of Chapter 11 again: g = 1, h = 0.5
    back = ql.attach(moved, "correlated-pair")["native"]
    assert back["ising_coupling_g"] == pytest.approx(1.0) and back["transverse_field_h"] == pytest.approx(0.5)
    assert ql.detach(spins(BOTH))["detached"] is False


def test_codiscovery_counts_the_two_routes():
    rep = ql.codiscover([load("two_spins"), spins(MOVED, name="field on the second spin"),
                         spins(BOTH, name="field on both spins"), spins(**COUNTER, name="another amplitude")])
    s = rep["summary"]
    assert (s["reached"], s["reached_through_closure"], s["obstructed"]) == (2, 1, 2)
    assert s["derivation_classes"] == {"A K L": ["two coupled spins (Module 11)"], "A O K L": ["field on the second spin"]}
    assert s["structure_residual_max"] < 1e-9 and s["rotation_residual_max"] < 1e-9 and s["law_residual_max"] < 1e-9
    text = ql.markdown(rep)
    assert "in 1 of them through the closure of the observable" in text
    assert "Reached in 2 realizations from 1 field (quantum information)" in text
    assert "| field on the second spin | quantum information | 2 spins-1/2 | AOKL | 6 | - |" in text
    assert "amplitude of the signal is 0.8212, not sin^2 theta = 0.7184" in text


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
    assert [r["word"] for r in c["rows"]] == ["AKL", "AO"]
    # a rotation reached through the closure: the report of the detachment and of the attachment say so
    moved = tmp_path / "moved.json"
    moved.write_text(json.dumps(spins(MOVED, name="field on the second spin").spec), encoding="utf-8")
    run("detach", str(moved), "--out-dir", str(tmp_path / "dm"))
    text = (tmp_path / "dm" / "detach.md").read_text(encoding="utf-8")
    assert "Derivation: AOKL; the rotation is reached through the closure of the observable" in text
    assert "closure of the observable (O): X0, Y0 Y1, Y0 Z1" in text or "closure of the observable (O): X0, Y0 Z1, Y0 Y1" in text
    run("attach", "--from", str(moved), "--to", "qubit", "--no-figure", "--out-dir", str(tmp_path / "am"))
    report = (tmp_path / "am" / "attach.md").read_text(encoding="utf-8")
    assert "In the source it is reached through the closure of the observable (AOKL)" in report
    assert ql.derive_bloch_rotation(ql.load(tmp_path / "am" / "attached.json"))["word"] == "AKL"


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
