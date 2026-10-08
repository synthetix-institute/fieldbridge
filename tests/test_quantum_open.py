"""The quantum write on an open carrier: analytic limits of the law, exactness checks, refusals, the pinned examples
and their classes, and the command."""
import json
import math
import subprocess
import sys
from pathlib import Path

import pytest

np = pytest.importorskip("numpy")
pytest.importorskip("sympy")
pytest.importorskip("scipy")
from scipy.special import ndtr

from fieldbridge.quantum import open as qo

ROOT = Path(__file__).resolve().parents[1]
EX = ROOT / "examples" / "quantum" / "open"


def spec(**over):
    base = json.loads((EX / "kerr_parametric_oscillator.json").read_text(encoding="utf-8"))
    base = json.loads(json.dumps(base))
    for key, value in over.items():
        if isinstance(value, dict) and isinstance(base.get(key), dict):
            base[key].update(value)
        else:
            base[key] = value
    return base


def test_examples_are_pinned_and_load():
    names = sorted(p.name for p in EX.glob("*.json"))
    assert names == ["few_photon_oscillator.json", "kerr_parametric_oscillator.json", "thermal_control.json"]
    for name in names:
        real = qo.load(EX / name)
        assert real.carrier.dim == real.spec["carrier"]["max_quanta"] + 1


def test_the_law_has_the_classical_and_the_closed_limits():
    # kappa >> sqrt(r): Phi(pi^(1/4) h / (D^(1/2) r^(1/4))) with D = kappa (2 nbar + 1)/8
    h, r, kappa = 0.02, 1e-3, 1.0
    canon = {"r": r, "h": h, "two_d": kappa / 4, "sigma0_sq": 0.25, "kappa": kappa, "gain_from": 0.0, "gain_to": 1.2}
    pr = {"to": 1.2, "from": 0.0, "rate": r, "hold": 0.0}
    expected = ndtr(math.pi ** 0.25 * h / (math.sqrt(kappa / 8) * r ** 0.25))
    assert qo.protocol_probability(canon, pr) == pytest.approx(expected, abs=1e-6)
    # kappa -> 0 with the vacuum as the seed: Phi(h sqrt(2 pi / r)) once the ramp has amplified the seed
    r = 0.05
    canon = {"r": r, "h": 0.05, "two_d": 0.0, "sigma0_sq": 0.25, "kappa": 0.0, "gain_from": 0.0, "gain_to": 3.0}
    pr = {"to": 3.0, "from": 0.0, "rate": r, "hold": 0.0}
    assert qo.protocol_probability(canon, pr) == pytest.approx(ndtr(0.05 * math.sqrt(2 * math.pi / r)), abs=1e-6)


def test_without_a_bias_the_probability_is_one_half_by_parity():
    s = spec(carrier={"max_quanta": 23}, parameters={"K": 0.3}, protocol={"sweep": {"parameter": "eps2", "from": 0.0, "to": 1.5, "rate": 0.5},
                                                                           "hold": 2.0, "bias": {"parameter": "h", "value": 0.0}})
    row = qo.derive_quantum_write(qo.load(s))
    assert row["status"] == "reached" and row["word"] == "GKL"
    assert row["law"]["P_exact"] == pytest.approx(0.5, abs=1e-9)
    assert row["law"]["P_law"] == 0.5


def test_a_quadratic_generator_follows_the_law_exactly():
    # K = 0: the evolution is Gaussian and the law holds at every rate
    s = spec(carrier={"max_quanta": 59}, parameters={"K": 0.0},
             protocol={"sweep": {"parameter": "eps2", "from": 0.0, "to": 1.2, "rate": 0.5}, "hold": 0.0,
                       "bias": {"parameter": "h", "value": 0.05}})
    row = qo.derive_quantum_write(qo.load(s))
    assert row["status"] == "reached"
    assert abs(row["law"]["difference"]) < 1e-6
    assert row["law"]["top_population"] < 1e-7


def test_the_examples_fall_in_their_classes():
    reached = qo.derive_quantum_write(qo.load(EX / "kerr_parametric_oscillator.json"))
    assert reached["status"] == "reached" and reached["class"] == "linear-stage write"
    assert abs(reached["law"]["difference"]) < 2e-4
    assert reached["canonical"]["threshold"] == pytest.approx(0.5, abs=1e-9)
    few = qo.derive_quantum_write(qo.load(EX / "few_photon_oscillator.json"))
    assert few["status"] == "obstructed at L" and few["class"] == "equilibrium write"
    assert few["law"]["P_exact"] < few["law"]["P_law"] - 0.05
    assert abs(few["law"]["P_exact"] - few["equilibrium"]["P_eq"]) <= 0.02


def test_the_thermal_control_lowers_the_probability_and_keeps_the_law():
    cold = qo.derive_quantum_write(qo.load(EX / "kerr_parametric_oscillator.json"), law=False)
    warm = qo.derive_quantum_write(qo.load(EX / "thermal_control.json"))
    assert warm["status"] == "reached"
    assert warm["law"]["P_law"] < cold["law"]["P_law"]
    assert abs(warm["law"]["difference"]) < 1e-3
    assert warm["canonical"]["two_d"] == pytest.approx(0.25 * (2 * 0.6 + 1))


def test_refusals():
    with pytest.raises(qo.SpecError):
        qo.load(spec(schema="fieldbridge-quantum/1"))
    with pytest.raises(qo.SpecError):
        qo.load(spec(dissipators=[]))
    with pytest.raises(qo.SpecError):
        qo.load(spec(hamiltonian=[{"coefficient": "eps2/2", "operator": "ad ad", "imaginary": True}]))
    with pytest.raises(qo.SpecError):
        qo.load(spec(hamiltonian=[{"coefficient": "-K", "operator": "b b"}]))
    with pytest.raises(qo.SpecError):
        qo.load(spec(hamiltonian=[{"coefficient": "__import__('os').system('true')", "operator": "na"}]))
    with pytest.raises(qo.SpecError):
        qo.load(spec(parameters={"K": 0.02, "kappa": 1.0, "nbar": 0.0, "eps2": 1.0}))


def test_obstructions_are_named():
    short = qo.derive_quantum_write(qo.load(spec(protocol={"sweep": {"parameter": "eps2", "from": 0.0, "to": 0.3, "rate": 0.1}})), law=False)
    assert short["status"] == "obstructed at G" and short["class"] == "no threshold"
    # a bias on the squeezed quadrature (h (a + a†) instead of i h (a† - a)) does not select
    off = spec(hamiltonian=[{"coefficient": "-K", "operator": "ad ad a a"},
                            {"coefficient": "eps2/2", "operator": "ad ad", "imaginary": True, "hc": True},
                            {"coefficient": "h", "operator": "ad", "hc": True}])
    row = qo.derive_quantum_write(qo.load(off), law=False)
    assert row["status"] == "obstructed at G" and row["class"] == "bias off the quadrature"


def test_the_command_writes_a_report(tmp_path):
    cmd = [sys.executable, "-B", "-m", "fieldbridge", "quantum", "write", str(EX / "few_photon_oscillator.json"),
           "--out-dir", str(tmp_path), "--no-law"]
    proc = subprocess.run(cmd, cwd=ROOT, capture_output=True, text=True, check=True)
    summary = json.loads(proc.stdout)
    assert summary["models"] == 1
    report = json.loads((tmp_path / "write.json").read_text(encoding="utf-8"))
    assert report["novelty_established"] is False
    assert report["rows"][0]["word"] == "GK"
    assert "The quantum write on open carriers" in (tmp_path / "write.md").read_text(encoding="utf-8")


def test_the_device_specifications_load_and_reach_the_canonical_form():
    # published parameters (a Josephson oscillator with a t^5 ramp; a Kerr-cat resonator with a tanh ramp): the
    # letters G and K hold without evolving the density operator, which takes minutes for these devices
    devices = sorted((EX / "devices").glob("*.json"))
    assert [p.name for p in devices] == ["grimm2020_kerr_cat.json", "yamaji2025_jpo.json"]
    for path in devices:
        real = qo.load(path)
        assert real.spec["provenance"]["origin"].startswith("published parameters")
        row = qo.derive_quantum_write(real, law=False)
        assert row["word"] == "GK", (path.name, row.get("obstruction"))
        assert row["canonical"]["threshold"] > 0 and row["law"]["P_law"] > 0.5
