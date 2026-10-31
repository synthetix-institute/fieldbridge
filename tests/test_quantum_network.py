"""Open quantum networks: pulses through a line into a detector, catch and release by a resonator, a coherent drive
against the mean-field equation, the one-mode swept write against the open-carrier solver, the decay of a stored
photon, refusals, the example and the command."""
import json
import math
import subprocess
import sys
from pathlib import Path

import pytest

np = pytest.importorskip("numpy")
pytest.importorskip("sympy")
pytest.importorskip("scipy")
from scipy.integrate import solve_ivp

from fieldbridge.quantum import network as qn
from fieldbridge.quantum import open as qo

ROOT = Path(__file__).resolve().parents[1]
EX = ROOT / "examples" / "quantum" / "network"
BASE = {"schema": qn.SCHEMA, "name": "test", "question": "q", "assumptions": ["a"], "provenance": {"origin": "test"}}


def network(**fields):
    return {**BASE, "hamiltonian": [], **fields}


def test_a_photon_sent_straight_into_the_detector_is_found_there():
    shape = "exp(-(t-3)**2/(2*0.5**2))"
    spec = network(subsystems=[{"name": "c", "kind": "mode", "max_quanta": 1}], duration=6.0,
                   ports=[{"name": "line", "operator": "c", "rate": 0,
                           "input": {"kind": "fock", "photons": 1, "shape": shape, "from": 0, "to": 6},
                           "output": {"shape": shape, "from": 0, "to": 6}}])
    rep = qn.measure(qn.load(spec), samples=0)
    port = rep["ports"][0]
    assert port["efficiency"] > 1 - 1e-6
    assert abs(rep["budget"]["difference"]) < 1e-9


def test_a_resonator_catches_a_rising_exponential_and_releases_it_into_the_output_mode():
    ta = 8.0
    spec = network(subsystems=[{"name": "r", "kind": "mode", "max_quanta": 1}], duration=2 * ta,
                   parameters={"kappa": 1.0}, observables={"nr": "nr"},
                   windows=[{"port": "line", "from": 0, "to": ta}],
                   ports=[{"name": "line", "operator": "r", "rate": "kappa",
                           "input": {"kind": "fock", "photons": 1, "shape": f"exp(kappa*(t-{ta})/2)",
                                     "from": 0, "to": ta},
                           "output": {"shape": f"exp(-kappa*(t-{ta})/2)", "from": ta, "to": 2 * ta}}])
    rep = qn.measure(qn.load(spec), samples=16)
    at_ta = [row for row in rep["trace"] if abs(row["t"] - ta) < 1e-9][0]
    reflected = math.exp(-ta)        # the part of the truncated exponential that the empty resonator reflects
    assert at_ta["nr"] == pytest.approx(1 - reflected, abs=1e-6)
    assert rep["windows"][0]["photons"] == pytest.approx(reflected, rel=1e-3)
    assert rep["ports"][0]["efficiency"] > 1 - 2 * reflected - 1e-5
    assert abs(rep["budget"]["difference"]) < 1e-9


def test_a_coherent_pulse_drives_the_mean_field_of_an_empty_resonator():
    kappa = 2.0
    beta = lambda t: 0.6 * np.exp(-(t - 3) ** 2 / 2)
    spec = network(subsystems=[{"name": "r", "kind": "mode", "max_quanta": 8}], duration=7.0,
                   parameters={"kappa": kappa}, observables={"x": "r + rd", "n": "nr"},
                   ports=[{"name": "line", "operator": "r", "rate": "kappa",
                           "input": {"kind": "coherent", "amplitude": "0.6*exp(-(t-3)**2/2)", "from": 0, "to": 7}}])
    rep = qn.measure(qn.load(spec), samples=14)
    times = [row["t"] for row in rep["trace"]]
    rhs = lambda t, a: -kappa / 2 * a - math.sqrt(kappa) * beta(t)    # Gardiner-Collett: b_in enters as -sqrt(kappa)
    sol = solve_ivp(rhs, (0, 7), [0.0], t_eval=times, rtol=1e-11, atol=1e-13)
    assert np.abs(np.array([row["x"] for row in rep["trace"]]) - 2 * sol.y[0]).max() < 1e-6
    assert np.abs(np.array([row["n"] for row in rep["trace"]]) - sol.y[0] ** 2).max() < 1e-6
    fine = np.linspace(0, 7, 7001)
    a = solve_ivp(rhs, (0, 7), [0.0], t_eval=fine, rtol=1e-11, atol=1e-13).y[0]
    out = (beta(fine) + math.sqrt(kappa) * a) ** 2
    emitted = float(np.sum(0.5 * (out[1:] + out[:-1]) * np.diff(fine)))
    assert rep["ports"][0]["emitted_into_line"] == pytest.approx(emitted, rel=1e-5)


def test_the_swept_write_of_the_open_carrier_is_a_network_with_one_mode():
    spec = json.loads((ROOT / "examples/quantum/open/few_photon_oscillator.json").read_text(encoding="utf-8"))
    spec["carrier"]["max_quanta"] = 10
    spec["protocol"]["sweep"].update({"to": 1.2, "rate": 0.4})
    spec["protocol"]["hold"] = 2.0
    spec["protocol"]["bias"]["off_at"] = 2.0
    ev_open = qo.evolve(qo.load(spec))
    ev_net = qn.evolve(qn.load(qn.from_open(spec)))
    assert np.abs(ev_open["rho_end"] - ev_net["rho"]).max() < 1e-6
    assert qo.prob_positive(ev_net["rho"]) == pytest.approx(qo.prob_positive(ev_open["rho_end"]), abs=1e-6)


def storage(wait: float) -> dict:
    """The photon starts in a lossless transmon, is swapped into a lossy storage mode, waits, and is swapped back to
    the transmon and on to the readout resonator, which emits it into the measured mode. Only the storage mode holds
    the excitation during the wait, so the photons found in the output mode scale exactly as exp(-ks wait)."""
    g = 20.0
    tau = math.pi / (2 * g)
    t_ret = tau + wait
    te = t_ret + 2 * tau
    return network(
        subsystems=[{"name": "r", "kind": "mode", "max_quanta": 1}, {"name": "q", "kind": "levels", "levels": 2},
                    {"name": "s", "kind": "mode", "max_quanta": 1}],
        parameters={"kappa": 1.0, "g": g, "ks": 0.05},
        controls={"gqs": [{"from": 0, "to": tau, "value": "g"}, {"from": t_ret, "to": t_ret + tau, "value": "g"}],
                  "grq": [{"from": t_ret + tau, "to": te, "value": "g"}]},
        hamiltonian=[{"coefficient": "grq", "operator": "rd q", "hc": True},
                     {"coefficient": "gqs", "operator": "qd s", "hc": True}],
        dissipators=[{"name": "storage loss", "rate": "ks", "operator": "s"}],
        ports=[{"name": "line", "operator": "r", "rate": "kappa",
                "output": {"shape": f"exp(-kappa*(t-{te})/2)", "from": te, "to": te + 8}}],
        initial={"q": 1}, duration=te + 8)


def test_a_stored_photon_decays_at_the_storage_rate_and_the_budget_closes():
    short, long_ = qn.measure(qn.load(storage(2.0)), samples=0), qn.measure(qn.load(storage(6.0)), samples=0)
    ratio = long_["ports"][0]["photons_in_output_mode"] / short["ports"][0]["photons_in_output_mode"]
    assert ratio == pytest.approx(math.exp(-0.05 * 4.0), rel=1e-6)
    lost = long_["removed_by_losses"]["storage loss"] - short["removed_by_losses"]["storage loss"]
    assert lost == pytest.approx(math.exp(-0.05 * 2.0) - math.exp(-0.05 * 6.0), rel=5e-3)   # up to the swap
    for rep in (short, long_):
        assert abs(rep["budget"]["difference"]) < 1e-9
        assert rep["checks"]["smallest_eigenvalue"] > -1e-9


def test_refusals():
    good = storage(1.0)
    with pytest.raises(qn.SpecError, match="Hermitian"):
        qn.load({**good, "hamiltonian": [{"coefficient": "g", "operator": "rd q"}]})
    with pytest.raises(qn.SpecError, match="unknown operator"):
        qn.load({**good, "hamiltonian": [{"coefficient": "g", "operator": "rd x", "hc": True}]})
    with pytest.raises(qn.SpecError, match="dimension"):
        qn.load({**good, "subsystems": [{"name": f"m{k}", "kind": "mode", "max_quanta": 3} for k in range(5)]})
    with pytest.raises(qn.SpecError, match="overlap"):
        qn.load({**good, "controls": {"gqs": [{"from": 0, "to": 2, "value": "g"}, {"from": 1, "to": 3, "value": "g"}]}})
    with pytest.raises(qn.SpecError, match="reserved"):
        qn.load({**good, "parameters": {**good["parameters"], "t": 1.0}})
    with pytest.raises(qn.SpecError, match="window"):
        qn.load({**good, "windows": [{"port": "line", "from": 0, "to": 1e6}]})


def test_the_example_loads_with_its_virtual_modes():
    names = sorted(p.name for p in EX.glob("*.json"))
    assert names == ["transmon_storage.json"]
    net = qn.load(EX / "transmon_storage.json")
    assert net.carrier.dim == 2 * 3 * 2 * 2 * 2 * 2       # r, q, s1, s2 and the virtual source and detector modes
    assert [s["name"] for s in net.subsystems][-2:] == ["lineIn", "lineOut"]


def test_command_writes_the_report(tmp_path):
    spec = tmp_path / "storage.json"
    spec.write_text(json.dumps(storage(1.0)), encoding="utf-8")
    out = tmp_path / "out"
    res = subprocess.run([sys.executable, "-B", "-m", "fieldbridge", "quantum", "network", str(spec),
                          "--out-dir", str(out), "--samples", "10"], cwd=ROOT, capture_output=True, text=True)
    assert res.returncode == 0, res.stderr
    rep = json.loads((out / "network.json").read_text(encoding="utf-8"))
    assert rep["rows"][0]["ports"][0]["photons_in_output_mode"] > 0.5
    assert "Budget" in (out / "network.md").read_text(encoding="utf-8")
