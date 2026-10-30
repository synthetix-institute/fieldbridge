"""The agent comparison of Module 17: its exact chains and its saved results against the equations."""
import importlib.util
import json
from pathlib import Path

import pytest

np = pytest.importorskip("numpy")
pytest.importorskip("scipy")

ROOT = Path(__file__).resolve().parents[1]
FOLDER = ROOT / "examples" / "population_agents"
spec = importlib.util.spec_from_file_location("population_agents_compare", FOLDER / "compare.py")
compare = importlib.util.module_from_spec(spec)
spec.loader.exec_module(compare)
RESULTS = json.loads((FOLDER / "results.json").read_text(encoding="utf-8"))


def test_exit_time_of_a_symmetric_walk():
    # at unit rates with reflection at 1, the mean time from n to n + 1 is n
    tau = compare.exit_time(lambda n: 1.0, lambda n: 1.0 if n > 1 else 0.0, 1, 5)
    assert tau == pytest.approx(1.0 + 2.0 + 3.0 + 4.0)


def test_extinction_time_of_a_pure_death_chain():
    # with no births the time from n to 0 is the sum of the waiting times
    log_tau = compare.extinction_time(lambda n: 0.0, lambda n: float(n), 4, 4)
    assert np.exp(log_tau) == pytest.approx(1 + 1 / 2 + 1 / 3 + 1 / 4)


def test_constant_noise_misses_the_exponent_and_a_state_dependent_noise_recovers_it():
    actions = compare.part_actions(None, True)
    for escape in actions.values():
        assert escape["diffusion"] == pytest.approx(escape["rates"], rel=0.02)
        assert abs(escape["constant"] / escape["rates"] - 1) > 0.3


def test_saved_agents_agree_with_the_equations():
    k = RESULTS["kuramoto"]
    assert np.allclose(k["r"]["10000"][4:], k["theory_r"][4:], atol=0.005)
    cycle = RESULTS["cycle"]
    ratio = [r * r * 100000 / t for r, t in zip(cycle["rms"]["100000"], cycle["linear_noise_variance_times_N"]) if t]
    assert all(0.8 < x < 1.2 for x in ratio)
    equation = RESULTS["closed"]["budworm_sweep_equation"]
    large = [x for x in RESULTS["fold_passage"]["jumps"] if x["N"] == 100000]
    for agents, eq in zip(large, equation):
        assert agents["up"] == pytest.approx(eq["jump_up"], abs=0.05)
        assert agents["down"] == pytest.approx(eq["jump_down"], abs=0.02)
    for row in RESULTS["closed"]["budworm_activated_jumps"]:
        assert row["up_median"] == pytest.approx(row["up_predicted"], abs=3 * row["up_sd"] / np.sqrt(row["reps"]))
