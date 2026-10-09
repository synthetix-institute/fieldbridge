"""The heredity module: the law of inheritance through a threshold against its normal form, the reduction at the
threshold, the classes of the examples and their controls, refusals, and the command line.

Statistical comparisons use the gate of the module, |P - P_law| < 4 stderr + 0.02; no solver residual is pinned.
"""
from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

pytest.importorskip("sympy")

from fieldbridge.heredity import SpecError, law, load  # noqa: E402
from fieldbridge.heredity.lineage import agrees, condition  # noqa: E402
from fieldbridge.heredity.predict import CLASSES, predict, reduction  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
EX = ROOT / "examples" / "heredity"
# every example with the class it must have; the slow bodies (fields with many cells, the polarity) are checked in
# the survey test only
FAST = ["benchmarks/normal_form.json", "controls/normal_form_no_dip.json", "controls/normal_form_lost.json",
        "chiral_autocatalysis_saito2007.json", "controls/chiral_autocatalysis_no_dip.json",
        "filament_baczynski2007.json", "controls/filament_short_division.json"]


def _within(p_sim, p_law, n):
    se = np.sqrt(max(p_sim * (1 - p_sim), 1e-4) / n)
    return abs(p_sim - p_law) < 4 * se + 0.02


@pytest.mark.parametrize("phi0, mu0, r, D", [(0.3, 0.6, 0.2, 1e-3), (1.0, 1.5, 0.5, 4e-3), (0.3, 0.3, 1.0, 4e-3)])
def test_the_law_with_the_cubic_decay_matches_the_normal_form(phi0, mu0, r, D):
    n = 4000
    p = law.simulate_normal_form(phi0, mu0, r, D, n, np.random.default_rng(1))
    assert _within(p, law.probability(phi0, mu0, r, D), n)


def test_the_linear_decay_misses_a_large_order():
    n = 4000
    p = law.simulate_normal_form(1.0, 1.5, 0.5, 1e-3, n, np.random.default_rng(2))
    assert law.probability(1.0, 1.5, 0.5, 1e-3, cubic=False) - p > 0.04
    assert _within(p, law.probability(1.0, 1.5, 0.5, 1e-3), n)


def test_without_noise_the_sign_at_the_crossing_decides():
    assert list(law.probabilities([1, 1, -1], [0.2, -0.1, -0.3], 0.0)) == [1.0, 0.0, 1.0]


def test_the_reduction_of_the_normal_form_is_the_identity():
    red = reduction(load(EX / "benchmarks" / "normal_form.json"))
    assert red["L_c"] == pytest.approx(1.0) and red["a"] == pytest.approx(1.0, rel=1e-4)
    assert red["b"] == pytest.approx(1.0, rel=1e-3)


def test_the_filament_threshold_is_the_euler_length_and_supercritical():
    red = reduction(load(EX / "filament_baczynski2007.json"))
    assert red["L_c"] == pytest.approx(1.0, abs=0.01)      # pi (kappa/F)^(1/2) = 1, ten segments
    assert red["b"] > 0


@pytest.mark.parametrize("name", FAST)
def test_examples_have_their_expected_class(name):
    lin = load(EX / name)
    assert predict(lin)["class"] == lin.expect


def test_lineages_of_the_normal_form_follow_the_law_and_keep_the_sign_without_noise():
    lin = load(EX / "benchmarks" / "normal_form.json")
    red = reduction(lin)
    row = condition(lin, red, 1.5, n=1500, generations=3, burn=1, seed=4)
    assert row["lnG_generation"] > 0.5 and row["Lambda"] > 10
    assert agrees(row, "law_body") and agrees(row, "law_history")
    clean = condition(lin, red, 1.5, noise=0.0, n=200, generations=3, burn=1, seed=5)
    assert clean["measured"] == 1.0


def test_a_binomial_partition_keeps_the_mean_and_adds_the_molecular_variance():
    lin = load(EX / "chiral_autocatalysis_saito2007.json")
    Q = np.full((20000, 2), 30.0)
    D = lin.body.divide(Q, np.random.default_rng(0))
    assert D.mean(0) == pytest.approx([15.0, 15.0], abs=0.02)
    assert D.var(0) == pytest.approx([30.0 / (4 * 16.0)] * 2, rel=0.05)     # x / (4 omega)


def test_specification_errors_are_refused():
    spec = json.loads((EX / "benchmarks" / "normal_form.json").read_text())
    for change in ({"size": "W"}, {"threshold": [1.5, 0.5]}, {"growth": {"law": "logistic", "rate": 0.1}},
                   {"division": {"at": 1.5, "at_relative": 1.5}}):
        bad = json.loads(json.dumps(spec))
        bad["lineage"].update(change)
        with pytest.raises(SpecError):
            load(bad)
    bad = json.loads(json.dumps(spec))
    bad["lineage"]["threshold"] = [1.2, 1.5]                  # no threshold in the interval
    with pytest.raises(SpecError):
        reduction(load(bad))


def test_the_commands_predict_and_survey(tmp_path):
    out = subprocess.run([sys.executable, "-B", "-m", "fieldbridge", "heredity", "predict",
                          str(EX / "controls" / "normal_form_no_dip.json")], capture_output=True, text=True, cwd=ROOT)
    assert out.returncode == 0 and "kept above the threshold" in out.stdout
    folder = tmp_path / "ex"
    folder.mkdir()
    for name in ("benchmarks/normal_form.json", "controls/normal_form_lost.json"):
        (folder / Path(name).name).write_text((EX / name).read_text())
    out = subprocess.run([sys.executable, "-B", "-m", "fieldbridge", "heredity", "survey", str(folder),
                          "--out-dir", str(tmp_path / "s")], capture_output=True, text=True, cwd=ROOT)
    assert out.returncode == 0 and "lost in the dip" in out.stdout
    assert {r["class"] for r in json.loads((tmp_path / "s" / "survey.json").read_text())["rows"]} <= set(CLASSES)


def test_the_card_command_writes_a_report_with_provenance(tmp_path):
    out = subprocess.run([sys.executable, "-B", "-m", "fieldbridge", "heredity", "card",
                          str(EX / "controls" / "normal_form_lost.json"), "--out-dir", str(tmp_path)],
                         capture_output=True, text=True, cwd=ROOT)
    assert out.returncode == 0
    rep = json.loads((tmp_path / "normal_form_lost.json").read_text())
    assert rep["class"] == "lost-in-the-dip" and rep["provenance"]["novelty_established"] is False
