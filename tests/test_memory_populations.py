"""Population models whose transitions are known exactly: the write points must lie where the equations put them."""
import json
from pathlib import Path

import pytest

np = pytest.importorskip("numpy")
pytest.importorskip("scipy")
pytest.importorskip("sympy")

from fieldbridge.memory import analysis, construct, discovery, spec  # noqa: E402

POP = Path(__file__).resolve().parents[1] / "examples" / "memory" / "populations"


def control_events(name):
    real = spec.load(POP / f"{name}.json")
    events = construct.construct(real, np.random.default_rng(1), check_write=False)["events"]
    return real, [e for e in events if e["source"] == "control"]


def test_synchronization_sets_in_at_twice_the_width_of_the_frequencies():
    """Ott-Antonsen: dz/dt = (K/2 - Delta) z - (K/2) |z|^2 z, so that K_c = 2 Delta in two directions at once."""
    real, events = control_events("kuramoto_population")
    assert len(events) == 1 and events[0]["kind"].startswith("degenerate critical mode")
    assert events[0]["value"] == pytest.approx(2 * real.params["Delta"], abs=1e-2)


def test_a_circle_of_states_is_a_continuum_whatever_the_seed():
    """Starts on a circle merge: with the seed of the catalog 23 states remain of 48 starts, which was read as discrete
    states with a zero threshold field."""
    real = spec.load(POP / "kuramoto_population.json")
    for seed in (discovery.card_rng(20260923, real.name), np.random.default_rng(1)):
        states, _, _, _, continuum = analysis.stored_states(real, seed, 48)
        assert len(states) >= 6 and continuum


def test_cyclic_dominance_starts_to_cycle_below_the_critical_mutation_rate():
    """At the symmetric state the eigenvalues are (b - a)/6 - mu +- i sqrt(3) (a + b)/6."""
    real, events = control_events("cyclic_dominance_mutation")
    a, b = real.params["a"], real.params["b"]
    assert len(events) == 1 and events[0]["kind"].startswith("Hopf")
    assert events[0]["value"] == pytest.approx((b - a) / 6, abs=1e-3)
    assert events[0]["omega"] == pytest.approx(np.sqrt(3) * (a + b) / 6, abs=1e-3)


def test_the_fittest_genotype_is_lost_at_the_error_threshold():
    real, events = control_events("error_threshold")
    assert len(events) == 1 and events[0]["kind"].startswith("transcritical")
    assert events[0]["value"] == pytest.approx(1 / real.params["f0"], abs=5e-3)


def test_budworm_outbreak_and_collapse_are_two_folds_of_one_hysteresis_loop():
    """The folds lie on r = 2 u^3 / (1 + u^2)^2, k = 2 u^3 / (u^2 - 1); for r = 0.55 at k = 6.596 and 11.582."""
    _, events = control_events("spruce_budworm")
    assert all(e["kind"].startswith("saddle-node") for e in events)
    assert sorted(e["value"] for e in events) == pytest.approx([6.596, 11.582], abs=0.03)


def test_population_specifications_state_a_source():
    for path in POP.glob("*.json"):
        assert json.loads(path.read_text(encoding="utf-8"))["provenance"]["source"]


def test_the_length_of_the_chaos_test_can_be_set_for_a_slow_system(tmp_path):
    """The food chain has lambda_1 of about 0.007: the default run of 1500 holds too few e-foldings to call it chaotic."""
    import subprocess
    import sys
    root = POP.parents[2]
    assert spec.load(root / "examples/dynamics/food_chain.json").control == "b1"
    run = subprocess.run([sys.executable, "-B", "-m", "fieldbridge", "memory", "dynamics", "examples/dynamics/lorenz.json",
                          "--time", "200", "--transient", "50", "--starts", "2", "--out-dir", str(tmp_path)],
                         cwd=root, capture_output=True, text=True)
    assert run.returncode == 0, run.stderr
    assert json.loads((tmp_path / "dynamics.json").read_text(encoding="utf-8"))["chaos"]["lyapunov"]["median"] > 0
