"""The computation module against bodies whose capacities are known in closed form or by symmetry.

Statistical comparisons use four times the measured spread of the estimator (0.017 for one signal at 40000 samples,
scaled with the square root of the number of signals and of the sample size); no solver residual is pinned.
"""
from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

pytest.importorskip("sympy")

from fieldbridge.computation import ipc, load  # noqa: E402
from fieldbridge.computation.card import card, simulate  # noqa: E402
from fieldbridge.computation.exact import ExactCapacities  # noqa: E402
from fieldbridge.computation.predict import linear, predict  # noqa: E402
from fieldbridge.computation.spec import SpecError  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
EX = ROOT / "examples" / "computation"
SD1 = 0.017                     # spread of the total of one signal at T = 40000


def _tol(signals: int, T: int) -> float:
    return 4 * SD1 * signals ** 0.5 * (40000 / T) ** 0.5


def test_targets_are_partitions_on_distinct_delays():
    assert len(list(ipc.targets(1, 4))) == 5
    assert len(list(ipc.targets(2, 2))) == 6
    assert len(list(ipc.targets(3, 2))) == 10
    assert len(list(ipc.targets(2, 3, "binary"))) == 6        # binary inputs: products of distinct delays only


def test_basis_is_orthonormal():
    rng = np.random.default_rng(0)
    for law in ("uniform", "gaussian"):
        x = rng.uniform(-1, 1, 200000) if law == "uniform" else rng.standard_normal(200000)
        P = np.array([ipc.polynomial(x, d, law) for d in range(1, 4)])
        assert np.allclose(P @ P.T / len(x), np.eye(3), atol=3e-2)


def test_one_mode_has_the_closed_form_profile_exactly_and_in_simulation():
    body = load(EX / "benchmarks" / "one_mode_map.json")
    pr = predict(body)
    assert pr["structure"] == "linear-memory" and pr["n_lin"] == 1
    assert np.allclose(pr["profile_degree_1"][:10], 0.36 * 0.64 ** np.arange(10), atol=1e-9)
    ex = ExactCapacities(body)
    assert abs(ex.capacity(((3, 1),)) - 0.36 * 0.64 ** 3) < 2e-5
    assert ex.capacity(((0, 2),)) < 1e-12
    res = simulate(body, {1: 30, 2: 4}, streams=8, length=2500)
    assert res["by_degree"][1] == pytest.approx(1.0, abs=_tol(1, res["T"])) and res["by_degree"][2] == 0.0


def test_benchmarks_and_controls_have_their_expected_classes():
    for path in sorted((EX / "benchmarks").glob("*.json")) + sorted((EX / "controls").glob("*.json")):
        body = load(path)
        assert predict(body)["structure"] == body.expect, path.name


def test_a_linear_chain_holds_three_and_one_measurement_time_of_its_last_stage_holds_one():
    spec = json.loads((EX / "benchmarks" / "linear_chain.json").read_text())
    assert linear(load(spec))["n_lin"] == 3
    spec["computation"]["observables"] = ["x3"]
    assert linear(load(spec))["n_lin"] == 1
    spec["computation"]["virtual_nodes"] = 3                 # three times per interval recover the hidden stages
    assert linear(load(spec))["n_lin"] == 3


def test_an_odd_body_has_no_even_degree_until_a_bias_breaks_it():
    # driven at amplitude 2, where degree 3 holds about 0.13 (at amplitude 1 it is spread over small targets)
    specs = {}
    for name in ("benchmarks/odd_oscillator.json", "controls/odd_oscillator_bias.json"):
        specs[name] = json.loads((EX / name).read_text())
        specs[name]["computation"]["amplitude"] = 2.0
    res = simulate(load(specs["benchmarks/odd_oscillator.json"]), {1: 20, 2: 6, 3: 4}, streams=8, length=1500)
    assert res["by_degree"][2] == 0.0 and res["by_degree"][3] > 0.03
    biased = simulate(load(specs["controls/odd_oscillator_bias.json"]), {1: 20, 2: 6}, streams=8, length=1500)
    assert biased["by_degree"][2] > 0.02


def test_a_square_cascade_ends_at_degree_two():
    body = load(EX / "benchmarks" / "square_cascade.json")
    ex = ExactCapacities(body, N=101).capacities({1: 30, 2: 8, 3: 4})
    assert ex["by_degree"][1] == pytest.approx(1.0, abs=5e-3)
    assert ex["by_degree"][2] == pytest.approx(1.0, abs=2e-2) and ex["by_degree"][3] < 1e-6


def test_measurement_noise_lowers_the_exact_capacity_and_matches_a_simulation():
    body = load(EX / "chemotaxis_methylation_tu2008.json")
    clean, noisy = ExactCapacities(body), ExactCapacities(body, noise=0.3)
    assert noisy.capacity(((1, 1),)) < clean.capacity(((1, 1),))
    sim = simulate(body, {1: 6}, streams=8, length=2500, noise=0.3)
    assert sim["profile_degree_1"][1] == pytest.approx(noisy.capacity(((1, 1),)), abs=0.03)


def test_one_spring_on_a_fixed_node_is_the_damped_oscillator_of_hauser_eq_1():
    spec = {"schema": "fieldbridge-springs/1", "name": "one spring", "question": "q", "assumptions": ["benchmark"],
            "provenance": {"source": "Hauser et al. 2011, Eq. (1)"},
            "network": {"positions": [[0.0, 0.0], [1.0, 0.0]], "springs": [[0, 1]], "k1": [4.0], "k3": [0.0],
                        "d1": [1.0], "d3": [0.0], "fixed": [0], "inputs": {"nodes": [1], "weights": [1.0]}},
            "computation": {"amplitude": 0.5, "hold": 0.5, "substeps": 10}}
    body = load(spec)
    x = body.steady()
    X = x.copy()
    X[0] += 0.1
    assert body.f(X[None, :], np.array([0.0]))[0][2] == pytest.approx(-0.4)
    pr = linear(body)
    assert pr["modes_reached_and_seen"] == 2 and pr["n_lin"] == 1
    assert pr["slowest_rate"] == pytest.approx(0.5, rel=1e-4)            # d/2m: the transverse motion is not seen


def test_specification_errors_are_refused():
    spec = json.loads((EX / "benchmarks" / "linear_chain.json").read_text())
    for change in ({"input": "w"}, {"virtual_nodes": 5}, {"law": "cauchy"}):
        bad = json.loads(json.dumps(spec))
        bad["computation"].update(change)
        with pytest.raises(SpecError):
            load(bad)
    bad = json.loads(json.dumps(spec))
    del bad["computation"]
    with pytest.raises(SpecError):
        load(bad)


def test_the_card_of_a_one_variable_body_has_exact_capacities():
    res = card(load(EX / "benchmarks" / "one_mode_map.json"), delays={1: 20, 2: 4})
    assert res["class"] == "linear-memory"
    assert res["runs"][0]["exact"]["by_degree"][1] == pytest.approx(1.0, abs=1e-3)


def test_the_command_predicts_a_class():
    out = subprocess.run([sys.executable, "-m", "fieldbridge", "computation", "predict",
                          str(EX / "benchmarks" / "running_sum.json")], capture_output=True, text=True, cwd=ROOT)
    assert out.returncode == 0 and "no fading memory" in out.stdout
