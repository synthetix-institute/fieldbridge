"""The closure rule as a construction: the index (Hamiltonian terms removed until the observable's closure rotates),
the closed counterpart, the excess of the signal over it, and the check that no line of the signal lies outside the
closure spectrum.

Expected values are analytic where possible: a realization that rotates has index 0, one line and no spread; the
nitrogen-vacancy centre loses its rotation to the single anisotropy term D Jz^2 (index 1, closure of seven operators
with three frequencies); two nuclear spins of different species need both Zeeman terms removed (index 2); a control
with a weighted observable rotates with one frequency but an amplitude offset, the 'not necessary' side of the rule.
The property that no Bohr line of the signal lies outside the closure spectrum is checked on random Hamiltonians."""
import json
import subprocess
import sys
from pathlib import Path

import pytest

np = pytest.importorskip("numpy")
pytest.importorskip("sympy")

from fieldbridge.quantum import closure as qc, language as ql

ROOT = Path(__file__).resolve().parents[1]
EX = ROOT / "examples" / "quantum"


def construct(name):
    return qc.construct(ql.load(EX / f"{name}.json"))


def test_rotating_realization_has_index_zero_one_line_and_no_spread():
    r = construct("nmr_spin")
    assert r["index"] == 0 and r["closure"]["dimension"] in (2, 3) and len(r["closure"]["frequencies"]) == 1
    e = r["excess"]
    assert e["lines"] == 1 and e["spread"] < 1e-9 and abs(e["mean_shift"]) < 1e-7
    assert abs(e["amplitude_offset"]) < ql.ALONG_TOL and r["rule_holds"]
    assert "scaling" not in r


def test_nv_centre_loses_the_rotation_to_the_anisotropy():
    r = construct("nv_centre")
    assert r["index"] == 1 and r["counterparts"][0]["removed"] == ["D [Jz Jz]"]
    assert r["closure"]["dimension"] == 7 and len(r["closure"]["frequencies"]) == 3
    e = r["excess"]
    assert e["lines"] == 3 and e["spread"] > 1.0 and e["weight_outside_counterpart_line"] > 0.3
    assert r["rule_holds"] and r["signal"]["lines_outside_closure"] == 0
    pts = r["scaling"]["points"]
    spreads = [p["spread"] for p in pts]
    assert spreads == sorted(spreads, reverse=True) and spreads[-1] < 0.1 * spreads[0]   # the excess vanishes with the term
    assert 0.7 < r["scaling"]["exponent"] < 1.3                                        # a first-order splitting


def test_two_species_need_both_zeeman_terms_removed():
    r = construct("heteronuclear_spins")
    assert r["index"] == 2 and sorted(r["counterparts"][0]["removed"]) == ["w1/2 [Z0]", "w2/2 [Z1]"]
    assert r["rule_holds"] and r["excess"]["lines"] == 2


def test_weighted_observable_rotates_with_an_amplitude_offset():
    r = qc.construct(ql.load(EX / "controls" / "weighted_observable.json"))
    assert r["index"] == 0 and r["excess"]["lines"] == 1 and r["excess"]["spread"] < 1e-9
    assert r["excess"]["amplitude_offset"] > 1e-3          # one frequency, the prepared state not along the observable
    assert r["rule_holds"]


def test_conserved_observable_has_no_index():
    spec = json.loads((EX / "two_spins.json").read_text(encoding="utf-8"))
    spec["hamiltonian"] = [{"coefficient": 1, "operator": "Z0"}]
    spec["observable"] = [{"coefficient": 1, "operator": "Z0"}]
    r = qc.construct(ql.load(spec))
    assert r["index"] is None and r["closure"]["dimension"] == 1 and "conserved" in r["obstruction"]
    assert r["excess"] is None and r["signal"]["frequencies"] == []


@pytest.mark.parametrize("seed", range(5))
def test_no_line_of_the_signal_lies_outside_the_closure_spectrum(seed):
    rng = np.random.default_rng(seed)
    d = 5
    A = rng.normal(size=(d, d)) + 1j * rng.normal(size=(d, d)); H = A + A.conj().T
    B = rng.normal(size=(d, d)) + 1j * rng.normal(size=(d, d)); O = B + B.conj().T
    O0 = O - np.trace(O) / d * np.eye(d)
    _, n, freqs = qc.closure_spectrum(H, O0)
    lines, _ = qc.signal_spectrum(H, O, qc.prepared_state(O))
    scale = max(freqs)
    for w, a in lines:
        if 2 * abs(a) > 1e-9:
            assert min(abs(w - f) for f in freqs) < 1e-6 * scale


def test_closure_command_writes_the_report(tmp_path):
    out = tmp_path / "closure"
    cmd = [sys.executable, "-B", "-m", "fieldbridge", "quantum", "closure", str(EX / "nv_centre.json"), str(EX / "nmr_spin.json"),
           "--out-dir", str(out)]
    res = subprocess.run(cmd, cwd=ROOT, capture_output=True, text=True, check=True)
    summary = json.loads(res.stdout)
    assert summary["realizations"] == 2 and summary["rotating"] == 1 and summary["index_one"] == 1 and summary["rule_holds"] == 2
    report = json.loads((out / "closure.json").read_text(encoding="utf-8"))
    assert report["rows"][0]["index"] == 1 and "implementation_sha256" in report
    assert "| nitrogen-vacancy centre ground-state spin |" in (out / "closure.md").read_text(encoding="utf-8")


def test_predictions_read_quantum_specifications(tmp_path):
    out = tmp_path / "pred"
    cmd = [sys.executable, "-B", "-m", "fieldbridge", "memory", "predictions", str(EX), "--out-dir", str(out)]
    res = subprocess.run(cmd, cwd=ROOT, capture_output=True, text=True, check=True)
    summary = json.loads(res.stdout)
    rows = json.loads((out / "predictions.json").read_text(encoding="utf-8"))["rows"]
    closure_rows = [r for r in rows if r["rule"] == "closure"]
    assert len(closure_rows) >= 10 and summary["rotating"] == sum(r["rotating"] for r in closure_rows)
    assert all(r["rule_holds"] for r in closure_rows)
