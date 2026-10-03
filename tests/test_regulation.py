"""REGULATE: the return of an output to its set point after a step of an input.

Every expected value is derived in closed form for the benchmark in question or stated by the source of a published
model; no test pins a solver residual. The controls are changes that must change the class: a leak, a removed
integrator, a cancellation away from its parameter surface, a turnover of the methylation level.
"""
import json
import subprocess
import sys
from pathlib import Path

import pytest

np = pytest.importorskip("numpy")
pytest.importorskip("scipy")
pytest.importorskip("sympy")

from fieldbridge.regulation import card as cardmod
from fieldbridge.regulation import gains, integrator, steady, step
from fieldbridge.regulation import spec as model

ROOT = Path(__file__).resolve().parents[1]
EX = ROOT / "examples" / "regulation"
BENCH = EX / "benchmarks"
SPECS = EX


def load(name: str) -> model.Regulated:
    return model.load(BENCH / f"{name}.json")


def settled(m, **over):
    p = m.pvec(m.u0, **over)
    rng = np.random.default_rng(1)
    L = integrator.conservation_laws(m, p, m.initial[None, :], rng)
    st = steady.steady_state(m, p, m.initial, L)
    assert st["converged"] and st["stable"]
    return p, st, L


# ------------------------------------------------------------------------------------------------ static gains
def test_static_gain_equals_bordered_determinant_and_finite_difference():
    # leaky integrator: dy/du = 1 / (lam + kI/delta)
    m = load("leaky_integrator")
    p, st, L = settled(m)
    g = gains.static_gain(m, p, st["q"], L)
    assert g["G"] == pytest.approx(1.0 / (1.0 + 0.5 / 0.05), rel=1e-9)
    assert g["bordered"] == pytest.approx(g["G"], rel=1e-9)
    assert cardmod._fd_gain(m, p, st["q"], L) == pytest.approx(g["G"], rel=1e-5)


def test_clamping_the_leaky_integrator_gives_the_remaining_fraction():
    # closed 1/(lam + kI/delta), open (z clamped) 1/lam: remaining fraction delta lam / (delta lam + kI)
    m = load("leaky_integrator")
    p, st, L = settled(m)
    z = next(r for r in gains.attenuation(m, p, st["q"], L) if r["variable"] == "z")
    assert z["role"] == "feedback"
    assert z["remaining_fraction"] == pytest.approx(0.05 / (0.05 + 0.5), rel=1e-9)


def test_output_dependent_on_the_input_enters_the_gain(tmp_path):
    # y = x + u with dx/dt = -x: dy/du = 1 through h_u alone
    spec = json.loads((BENCH / "no_path.json").read_text(encoding="utf-8"))
    spec["drift"] = {"x": "-x", "y": "-y"}
    spec["regulation"]["output"] = "x + u"
    path = tmp_path / "hu.json"
    path.write_text(json.dumps(spec), encoding="utf-8")
    m = model.load(path)
    p, st, L = settled(m)
    assert gains.static_gain(m, p, st["q"], L)["G"] == pytest.approx(1.0, rel=1e-12)


# ------------------------------------------------------------------------------------------------ integrators
def test_antithetic_pair_is_found_as_a_difference():
    # d(z2 - z1)/dt = theta x2 - mu = theta (x2 - mu/theta)
    m = load("antithetic")
    p, st, _ = settled(m)
    info = integrator.find_integrator(m, p, st["q"][None, :], np.random.default_rng(2), y0_steady=st["y"])
    assert info["found"] and info["stage"] == 1
    w = info["coefficients"]["w"]
    assert set(w) == {"z1", "z2"} and w["z2"] == pytest.approx(1.0) and w["z1"] == pytest.approx(-1.0)
    assert info["k_I"] == pytest.approx(m.params["theta"], rel=1e-8)
    assert info["set_point"] == pytest.approx(m.params["mu"] / m.params["theta"], rel=1e-8)


def test_logarithmic_integrator_needs_the_logarithmic_columns():
    # d ln z/dt = k (y - y0)
    m = load("autocatalytic_integrator")
    p, st, _ = settled(m)
    centers = st["q"][None, :]
    lin = integrator.find_integrator(m, p, centers, np.random.default_rng(3), y0_steady=st["y"], logs=False,
                                     stage2=False)
    assert not lin["found"]
    log = integrator.find_integrator(m, p, centers, np.random.default_rng(3), y0_steady=st["y"])
    assert log["found"] and log["stage"] == 1 and log["coefficients"] == {"w": {}, "v": {"z": pytest.approx(1.0)}}
    assert log["k_I"] == pytest.approx(m.params["k"], rel=1e-8)


def test_concentration_robustness_carries_a_logarithmic_integrator_and_a_conservation_law():
    # A + B -> 2B, B -> A, A <-> C: d ln B/dt = k1 (A - k2/k1); A + B + C is conserved
    m = load("acr_network")
    p, st, L = settled(m)
    assert len(L) == 1 and np.allclose(np.abs(L[0]) / np.abs(L[0]).max(), 1.0)
    assert st["y"] == pytest.approx(m.params["k2"] / m.params["k1"], rel=1e-9)
    info = integrator.find_integrator(m, p, st["q"][None, :], np.random.default_rng(4), y0_steady=st["y"])
    assert info["found"] and info["coefficients"]["v"] == {"B": pytest.approx(1.0)}
    assert info["conservation_laws"] == 1


def test_feedforward_product_needs_a_gain_linear_in_the_state():
    # phi = (k1/k3) x - y has dphi/dt = k2 x (y - y0), y0 = k1 k4 / (k2 k3)
    m = load("feedforward_product")
    p, st, _ = settled(m)
    centers = st["q"][None, :]
    one = integrator.find_integrator(m, p, centers, np.random.default_rng(5), y0_steady=st["y"], stage2=False)
    assert not one["found"]
    two = integrator.find_integrator(m, p, centers, np.random.default_rng(5), y0_steady=st["y"])
    assert two["found"] and two["stage"] == 2
    assert two["coefficients"]["w"] == {"x": pytest.approx(1.0), "y": pytest.approx(-1.0)}
    g = two["gain_coefficients"]
    assert g["x"] == pytest.approx(m.params["k2"]) and abs(g["1"]) < 1e-8


def test_subtractive_feedforward_has_an_integrator_only_where_it_cancels():
    # phi = a x - y integrates y - 1 only if a k3 = k1 and a k4 = k2, i.e. k1 k4 = k2 k3
    m = load("feedforward_subtractive")
    p, st, _ = settled(m)
    assert integrator.find_integrator(m, p, st["q"][None, :], np.random.default_rng(6), y0_steady=st["y"])["found"]
    p2, st2, L2 = settled(m, k1=1.2)
    assert not integrator.find_integrator(m, p2, st2["q"][None, :], np.random.default_rng(6),
                                          y0_steady=st2["y"])["found"]
    assert abs(gains.static_gain(m, p2, st2["q"], L2)["S"]) > 1e-3


def test_conservation_law_is_not_an_integrator():
    m = load("closed_exchange")
    p, st, _ = settled(m)
    info = integrator.find_integrator(m, p, st["q"][None, :], np.random.default_rng(7), y0_steady=st["y"])
    assert not info["found"] and info["conservation_laws"] == 1


# ------------------------------------------------------------------------------------------------ dynamics
def test_integral_gain_threshold_of_the_second_order_plant():
    # s^3 + 2 lam s^2 + lam^2 s + lam kI = 0 is stable iff kI < 2 lam^2
    m = load("third_order_loop")
    lam = m.params["lam"]
    for kI, stable in ((2 * lam ** 2 * 0.99, True), (2 * lam ** 2 * 1.01, False)):
        p = m.pvec(m.u0, kI=kI)
        q = np.array([m.params["y0"], m.params["y0"], m.u0 - lam * m.params["y0"]])   # the equilibrium
        assert np.allclose(m.f(q, p), 0.0)
        assert bool(steady.rates(m, p, q).real.max() < 0) == stable


def test_response_integral_equals_the_change_of_the_integrator():
    # PI loop: kI int (y - y0) dt = delta z = delta u (an identity: a calibration of the simulation)
    m = load("pi_loop")
    p, st, L = settled(m)
    info = integrator.find_integrator(m, p, st["q"][None, :], np.random.default_rng(8), y0_steady=st["y"])
    r = step.step_response(m, p, st["q"], 3.0, info, L)
    assert r["settled"] and r["calibration_ratio"] == pytest.approx(1.0, abs=1e-8)
    assert m.params["kI"] * r["integral_of_deviation"] == pytest.approx(3.0 - m.u0, rel=1e-8)
    assert r["final"] == pytest.approx(0.0, abs=1e-8) and abs(r["peak"]) > 0.1


@pytest.mark.parametrize("path", sorted(BENCH.glob("*.json")), ids=lambda p: p.stem)
def test_benchmark_classes(path):
    m = model.load(path)
    res = cardmod.card(m, samples=6, seed=0)
    assert res["class"] == m.spec["regulation"]["expect"]


# ------------------------------------------------------------------------------------------------ refusals
def _write(tmp_path, change):
    spec = json.loads((BENCH / "pi_loop.json").read_text(encoding="utf-8"))
    change(spec)
    path = tmp_path / "spec.json"
    path.write_text(json.dumps(spec), encoding="utf-8")
    return path


@pytest.mark.parametrize("change", [
    lambda s: s.pop("regulation"),
    lambda s: s["regulation"].update(input="w"),
    lambda s: s["regulation"].update(output="y + undeclared"),
    lambda s: s["regulation"].update(output="__import__('os').getcwd()"),
    lambda s: s["regulation"].update(steps=[2.0]),
    lambda s: s["regulation"].update(initial={"nope": 1.0}),
])
def test_refusals(tmp_path, change):
    with pytest.raises(model.SpecError):
        model.load(_write(tmp_path, change))


def test_zero_set_point_is_judged_by_the_open_loop_gain():
    # y0 = 0: S = (u/y) G is undefined; the gain is zero against the open-loop gain 1/lam of the clamp of z
    m = load("pi_loop_zero_set_point")
    p, st, L = settled(m)
    assert abs(st["y"]) < 1e-12
    assert gains.reference_gain(m, p, st["q"], L) == pytest.approx(1.0 / m.params["lam"], rel=1e-9)
    res = cardmod.card(m, samples=6)
    assert res["S"] is None and res["class"] == "perfect-adaptation"


def test_dilution_of_a_logarithmic_integrator_moves_its_set_point():
    # dz/dt = k z (y - y0) - delta z: d ln z/dt = k (y - y0 - delta/k), so the return is exact to y0 + delta/k,
    # whereas a leak of a linear integrator leaves a residual (test_clamping_the_leaky_integrator...)
    m = load("autocatalytic_diluted")
    p, st, _ = settled(m)
    target = m.params["y0"] + m.params["delta"] / m.params["k"]
    assert st["y"] == pytest.approx(target, rel=1e-9)
    info = integrator.find_integrator(m, p, st["q"][None, :], np.random.default_rng(9), y0_steady=st["y"])
    assert info["found"] and info["coefficients"]["v"] == {"z": pytest.approx(1.0)}
    assert info["set_point"] == pytest.approx(target, rel=1e-8)
    assert cardmod.card(m, samples=6)["class"] == "perfect-adaptation"


def test_methylation_level_is_an_integrator_of_a_function_of_the_activity():
    # Tu, Shimizu and Berg 2008: dm/dt = F(a), F decreasing through a0 = 1/3 (piecewise linear, written with abs);
    # no constant-gain integrator exists, stage 3 finds phi = -m with set point a0
    m = model.load(SPECS / "chemotaxis_tu2008.json")
    p, st, _ = settled(m)
    assert st["y"] == pytest.approx(1 / 3, rel=1e-9)
    info = integrator.find_integrator(m, p, st["q"][None, :], np.random.default_rng(10), y0_steady=st["y"])
    assert info["found"] and info["stage"] == 3
    assert info["coefficients"]["w"] == {"m": pytest.approx(-1.0)}
    no3 = integrator.find_integrator(m, p, st["q"][None, :], np.random.default_rng(10), y0_steady=st["y"],
                                     stage3=False)
    assert not no3["found"]


def test_methylation_turnover_removes_the_integrator():
    m = model.load(SPECS / "controls" / "chemotaxis_tu2008_control1.json")
    p, st, _ = settled(m)
    assert not integrator.find_integrator(m, p, st["q"][None, :], np.random.default_rng(11),
                                          y0_steady=st["y"])["found"]


def test_piecewise_methylation_rate_matches_the_source():
    # F(0) = 1, F(a0) = 0 and F(1) = -2 (Fig. 2 legend of the source); F is written through the activity, so it is
    # evaluated at the methylation levels where the activity is 0+, a0 and 1-
    from scipy.optimize import brentq
    m = model.load(SPECS / "chemotaxis_tu2008.json")
    p = m.pvec(m.u0)
    for target, expected in ((1e-9, 1.0), (1 / 3, 0.0), (1 - 1e-9, -2.0)):
        mval = brentq(lambda x: m.y(np.array([x]), p) - target, -50.0, 50.0, xtol=1e-14)
        assert m.f(np.array([mval]), p)[0] == pytest.approx(expected, abs=1e-6)


def test_quasi_integral_error_matches_the_source():
    # Qian and Del Vecchio 2018, Eq. (2.8): at steady state u - x = eps gamma (z1 - z2) / k, and the leak
    # leaves an error that shrinks with eps; without dilution (gamma = 0) the error vanishes
    for name, eps in (("qian2018_quasi.json", 0.02), ("controls/qian2018_quasi_control1.json", 1.0)):
        m = model.load(SPECS / name)
        p, st, _ = settled(m, d=5.0)
        x, z1, z2 = st["q"]
        pr = m.params
        assert pr["u"] - x == pytest.approx(eps * pr["gamma"] * (z1 - z2) / pr["k"], rel=1e-8)
    m = model.load(SPECS / "controls" / "qian2018_quasi_control2.json")
    p, st, _ = settled(m, d=5.0)
    assert st["q"][0] == pytest.approx(m.params["u"], rel=1e-9)


# ------------------------------------------------------------------------------------------------ published models
@pytest.mark.parametrize("path", sorted(EX.glob("*.json")) + sorted((EX / "controls").glob("*.json")),
                         ids=lambda p: p.stem)
def test_published_models_and_controls_get_their_class(path):
    m = model.load(path)
    res = cardmod.card(m, samples=4, seed=0)
    expect = m.spec["regulation"]["expect"]
    assert res["class"] in (expect if isinstance(expect, list) else [expect])


def test_every_example_names_its_source_and_its_assumptions():
    for path in EX.rglob("*.json"):
        s = json.loads(path.read_text(encoding="utf-8"))
        assert s["schema"] == "fieldbridge-regulation/1" and s["provenance"]["source"] and s["assumptions"]
        assert s["regulation"].get("expect")


# ------------------------------------------------------------------------------------------------ command line
def _run(*args):
    return subprocess.run([sys.executable, "-B", "-m", "fieldbridge", "regulation", *args], cwd=ROOT,
                          capture_output=True, text=True)


def test_command_line_writes_a_card_whose_certificate_passes(tmp_path):
    out = tmp_path / "antithetic"
    assert _run("card", str(EX / "antithetic_briat2016.json"), "--out-dir", str(out), "--samples", "4").returncode == 0
    report = json.loads((out / "regulation.json").read_text(encoding="utf-8"))
    assert report["command"] == "regulation card" and report["novelty_established"] is False
    assert report["certificate"]["class"] == "perfect-adaptation"
    assert report["certificate"]["integrator"]["coefficients"]["w"] == {"z1": pytest.approx(-1.0),
                                                                       "z2": pytest.approx(1.0)}
    assert "Return to a set point" in (out / "regulation.md").read_text(encoding="utf-8")
    checked = _run("check", str(out / "regulation.json"))
    assert checked.returncode == 0 and "simplifies to zero" in checked.stdout


def test_a_certificate_fails_against_a_changed_specification(tmp_path):
    out = tmp_path / "pi"
    assert _run("card", str(BENCH / "pi_loop.json"), "--out-dir", str(out), "--samples", "4").returncode == 0
    spec = json.loads((BENCH / "pi_loop.json").read_text(encoding="utf-8"))
    spec["drift"]["z"] = "kI*(y - y0) - 0.05*z"              # a leak: the recorded integrator no longer holds
    changed = tmp_path / "pi_leaky.json"
    changed.write_text(json.dumps(spec), encoding="utf-8")
    checked = _run("check", str(out / "regulation.json"), "--spec", str(changed))
    assert checked.returncode == 1 and "FAILED" in checked.stdout


def test_command_line_survey_compares_with_the_expected_classes(tmp_path):
    res = _run("survey", str(BENCH / "pi_loop.json"), str(BENCH / "leaky_integrator.json"), "--out-dir",
               str(tmp_path), "--samples", "4")
    assert res.returncode == 0 and "0 not as expected" in res.stdout
    survey = json.loads((tmp_path / "survey.json").read_text(encoding="utf-8"))
    assert [r["card"]["class"] for r in survey["results"]] == ["perfect-adaptation", "partial-adaptation"]
