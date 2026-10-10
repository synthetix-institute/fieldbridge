"""Tests of the decision family: closed forms of the onset and of passage statistics, the reductions of published
populations at their thresholds, refusals, controls that change the answer, and the command line."""
from __future__ import annotations

import json
import subprocess
import sys
from math import log, pi, sqrt
from pathlib import Path

import numpy as np
import pytest
from scipy.stats import norm

pytest.importorskip("sympy")

from fieldbridge.decision import load, onset, passage  # noqa: E402
from fieldbridge.decision.predict import predict  # noqa: E402
from fieldbridge.decision.spec import SpecError  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
EX = ROOT / "examples" / "decision"


def spec_of(rel):
    return json.loads((EX / rel).read_text())


# ------------------------------------------------------------------------------------------------ closed forms
def test_lorentzian_onset_and_growth_rate():
    # 1 = (K/2)(b1 - i a1)/(lambda + gamma): K_c = 2 gamma/b1, mu = K b1/2 - gamma, Omega = K a1/2
    g = onset.lorentzian(0.1)
    assert onset.critical_coupling(g, 0.0, 1.0, scale=0.1)["K_c"] == pytest.approx(0.2, rel=1e-6)
    out = onset.growth_rate(g, 0.3, 0.8, 0.5, scale=0.1)
    assert out["mu"] == pytest.approx(0.1, rel=1e-4)
    assert out["Omega"] == pytest.approx(0.075, rel=1e-4)


def test_gaussian_onset_and_no_onset_without_b1():
    s = 0.2
    assert onset.critical_coupling(onset.gaussian(s), 0.0, 1.0, scale=s)["K_c"] == pytest.approx(sqrt(8 / pi) * s,
                                                                                                   rel=1e-6)
    assert onset.critical_coupling(onset.uniform(s), 0.4, 0.0, scale=s) is None       # a cosine coupling shifts only


def test_window_ratios():
    assert passage.q_ratio(1) == pytest.approx(norm.ppf(0.95) / norm.ppf(0.55), rel=1e-10)
    assert passage.q_ratio(1) == pytest.approx(13.0896, abs=1e-4)
    assert passage.q_ratio(2) == pytest.approx(sqrt(log(10) / log(10 / 9)), rel=1e-10)
    zs = [0.0, 0.5, 1.0, 2.0]
    ratios = [passage.folded_ratio(z) for z in zs]
    assert ratios[0] == pytest.approx(passage.q_ratio(1))
    assert all(a > b for a, b in zip(ratios, ratios[1:]))                              # a bias narrows the window


@pytest.mark.parametrize("d", [1, 2])
def test_window_of_exponential_growth_from_a_gaussian_seed(d):
    rng = np.random.default_rng(d)
    lam = 0.7
    xi = np.linalg.norm(rng.standard_normal((200000, d)), axis=1) * 1e-6
    out = passage.check(np.log(1.0 / xi) / lam, passage.Exponent(np.array([0.0, 1e3]), np.array([lam, lam])), d)
    assert out["ratio"] == pytest.approx(1.0, abs=0.01)


def test_leading_dimension():
    assert passage.leading_dimension(np.array([[0.1, -1.0], [1.0, 0.1]])) == 2               # a rotation that grows
    assert passage.leading_dimension(np.array([[0.3, 0.0], [0.0, -1.0]])) == 1


def test_fold_normal_form():
    # x' = -h - x^2 has its fold at h = 0: a_h = -1, b = -1, ghost passage pi/(h)^(1/2)
    nf = passage.fold_normal_form(lambda x, h: np.array([-h - x[0] ** 2, -x[1]]),
                                  lambda x, h: np.array([[-2 * x[0], 0.0], [0.0, -1.0]]), np.array([0.0, 0.0]), 0.0)
    assert nf["coefficient"] == pytest.approx(pi, rel=1e-5)


# ------------------------------------------------------------------------------------------------ reductions
def test_published_thresholds():
    ising = load(EX / "write/ising_glauber1963.json").body.reduce()
    assert ising["c_star"] == pytest.approx(1.0, abs=1e-7)
    assert ising["a"] == pytest.approx(1.0, rel=1e-3)                                     # finite differences
    v = 5.0
    bees = load(EX / "write/honeybees_pais2013.json").body.reduce()
    assert bees["c_star"] == pytest.approx(4 * v ** 3 / (v ** 2 - 1) ** 2, rel=1e-8)    # Pais et al. 2013, eq. (4)
    chiral = load(EX / "write/chiral_autocatalysis_saito2007.json").body.reduce(s_max=0.5)
    assert chiral["c_star"] == pytest.approx(2 * (1 + 1 / 4) / sqrt(4e-3), rel=1e-8)    # Saito et al. 2007, eq. (31)
    assert chiral["b"] > 0
    ww = load(EX / "write/decision_network_wong_wang2006.json").body.reduce(s_max=0.01)
    assert ww["c_star"] == pytest.approx(10.68, abs=0.01)
    assert ww["b"] < 0                                                                     # subcritical
    dopo = load(EX / "write/coherent_ising_machine_wang2013.json").body.reduce()
    assert dopo["c_star"] == pytest.approx(0.4, abs=1e-8)


def test_chemical_langevin_noise_of_the_transitions():
    b = load(EX / "write/chiral_autocatalysis_saito2007.json").body
    red = b.reduce(s_max=0.5)
    r, s = red["y_sym"]
    a = red["c_star"] - r - s
    D = np.diag([(1 + 4e-3 * r * r) * a + r, (1 + 4e-3 * s * s) * a + s]) / 2.0
    assert b.D(red["y_sym"], red["c_star"], 1.0) == pytest.approx(D, rel=1e-10)


def test_history_law_tends_to_the_closed_form_for_slow_sweeps():
    b = load(EX / "write/ising_glauber1963.json").body
    red = b.reduce()
    # lambda = K - 1, h_s and D_s are constant at the symmetric state: only the start could differ, and it is
    # forgotten as exp(-a (c* - c0)^2 / 2r)
    lin = b.law(red, 1e-3, 800, 1e-3)
    hist = b.law_history(red, 1e-3, 800, 1e-3, 0.5, 1.6, points=2000)
    assert hist["probit"] == pytest.approx(lin["probit"], rel=0.01)


def test_diverse_units_threshold_and_frozen_bias():
    b = load(EX / "write/diverse_units_tessone2006.json").body
    red = b.reduce()
    Cs = red["c_star"]
    a = np.linspace(-9, 9, 40001)
    x = b.profile(a, 0.0, Cs)
    trap = getattr(np, "trapezoid", None) or np.trapz
    w = norm.pdf(a)
    assert Cs * trap(w / (Cs - 1 + 3 * x ** 2), a) == pytest.approx(1.0, abs=2e-4)
    # df/dX = C and df/dh = 1 for every unit: s_q N^(1/2) = C* <x*^2>^(1/2)
    assert red["frozen_per_sqrtN"] == pytest.approx(Cs * sqrt(trap(w * x ** 2, a)), rel=1e-3)


def test_the_mean_state_lags_a_fast_ramp():
    b = load(EX / "write/decision_network_wong_wang2006.json").body
    red = b.reduce(s_max=0.01)
    N, r = 16000, 5.0
    h = 0.8 * sqrt(b.D_s(red, N)) * (red["a"] * r) ** 0.25 / (pi ** 0.25 * red["h_s"])
    assert b.law_history(red, h, N, r, 0.0, 25.0)["P"] - b.law(red, h, N, r)["P"] == pytest.approx(0.023, abs=0.005)


# ------------------------------------------------------------------------------------------------ classes and controls
def test_controls_change_the_class():
    assert predict(load(EX / "controls/diverse_units_random_sample.json"))["class"] == "set-by-the-sample"
    assert predict(load(EX / "write/diverse_units_tessone2006.json"))["class"] == "follows-the-bias"
    assert predict(load(EX / "passage/macrospin_stoner_wohlfarth.json"))["class"] == "rotation-seed"
    pr = predict(load(EX / "passage/coherent_ising_machine_step.json"))
    assert pr["class"] == "reflection-seed" and pr["rate_end"] == pytest.approx(0.4, rel=1e-4)
    josephson = predict(load(EX / "controls/josephson_no_onset.json"))
    assert josephson["class"] == "no-onset"
    assert abs(josephson["reduction"]["b1"]) < 1e-6


def test_van_der_pol_coupling_function():
    cf = load(EX / "synchronization/van_der_pol.json").body.coupling_function()
    assert cf["a1"] == pytest.approx(-0.093, abs=2e-3)
    assert cf["b1"] == pytest.approx(0.499, abs=2e-3)


def test_examples_are_pinned():
    names = sorted(str(p.relative_to(EX)) for p in EX.rglob("*.json"))
    assert names == ["controls/diverse_units_random_sample.json", "controls/josephson_no_onset.json",
                     "passage/chiral_autocatalysis_step.json", "passage/coherent_ising_machine_step.json",
                     "passage/decision_network_step.json", "passage/macrospin_stoner_wohlfarth.json",
                     "synchronization/brusselator.json", "synchronization/fitzhugh_nagumo.json",
                     "synchronization/goodwin.json", "synchronization/predator_prey.json",
                     "synchronization/van_der_pol.json", "write/chiral_autocatalysis_saito2007.json",
                     "write/coherent_ising_machine_wang2013.json", "write/decision_network_wong_wang2006.json",
                     "write/diverse_units_tessone2006.json", "write/honeybees_pais2013.json",
                     "write/ising_glauber1963.json"]
    for p in EX.rglob("*.json"):
        s = json.loads(p.read_text())
        assert s["provenance"]["source"] and s["assumptions"] and s["expect"]


# ------------------------------------------------------------------------------------------------ refusals
@pytest.mark.parametrize("edit, message", [
    (lambda s: s.update(schema="fieldbridge-memory/1"), "schema"),
    (lambda s: s.pop("question"), "question"),
    (lambda s: s.update(kind="lattice"), "kind"),
    (lambda s: s.update(target="synchronization"), "serves"),
    (lambda s: s["population"].pop("bias"), "bias"),
    (lambda s: s["population"]["control"].update(range=[1.2, 1.5]), "does not lose stability"),
    (lambda s: s["population"].update(size="Nunits"), "declared parameter"),
    (lambda s: s["population"]["noise"]["transitions"][0]["change"].update(z=1), "undeclared variable"),
    (lambda s: s["protocol"]["sweep"].update(rates=[-0.01]), "positive rates"),
])
def test_refusals(edit, message):
    s = spec_of("write/ising_glauber1963.json")
    edit(s)
    with pytest.raises(SpecError, match=message):
        pop = load(s)
        pop.body.reduce()


def test_a_coefficient_cannot_import():
    s = spec_of("write/ising_glauber1963.json")
    s["population"]["noise"]["transitions"][0]["rate"] = "__import__('os').getcwd()"
    with pytest.raises(SpecError):
        load(s)


# ------------------------------------------------------------------------------------------------ command line
def run(*args):
    return subprocess.run([sys.executable, "-m", "fieldbridge", "decision", *args], cwd=ROOT, capture_output=True,
                          text=True, timeout=600)


def test_cli_predict():
    out = run("predict", "examples/decision/passage/coherent_ising_machine_step.json")
    assert out.returncode == 0, out.stderr
    assert "one-component seed" in out.stdout


def test_cli_card(tmp_path):
    out = run("card", "examples/decision/controls/diverse_units_random_sample.json", "--out-dir", str(tmp_path))
    assert out.returncode == 0, out.stderr
    rep = json.loads((tmp_path / "diverse_units_random_sample.json").read_text())
    assert rep["class"] == "set-by-the-sample"
    assert rep["provenance"]["implementation_sha256"]
    assert "Boundary" in (tmp_path / "diverse_units_random_sample.md").read_text()


def test_card_without_an_onset_does_not_simulate():
    from fieldbridge.decision.card import card
    res = card(load(EX / "controls/josephson_no_onset.json"), simulate=True)
    assert res["class"] == "no-onset" and "skipped" in res["simulation"]


def test_cli_survey(tmp_path):
    folder = tmp_path / "specs"
    folder.mkdir()
    for rel in ("passage/coherent_ising_machine_step.json", "controls/diverse_units_random_sample.json"):
        (folder / Path(rel).name).write_text((EX / rel).read_text())
    out = run("survey", str(folder), "--out-dir", str(tmp_path / "out"))
    assert out.returncode == 0, out.stdout + out.stderr
    rows = json.loads((tmp_path / "out" / "survey.json").read_text())["rows"]
    assert {r["class"] for r in rows} == {"reflection-seed", "set-by-the-sample"}
