"""Birth and death rates of a population: the exponent of the time to leave a state, against the exact chain."""
import json
from pathlib import Path

import pytest

np = pytest.importorskip("numpy")
pytest.importorskip("scipy")
pytest.importorskip("sympy")

from fieldbridge.memory import population, spec  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
FOLDER = ROOT / "examples" / "memory" / "populations"


def path(report, start, stop):
    return next(p for p in report["paths"] if p["from"] == pytest.approx(start, abs=1e-2)
                and p["to"] == pytest.approx(stop, abs=1e-2))


def epidemic(R0):
    return {"schema": spec.SCHEMA, "question": "q", "assumptions": ["none"], "provenance": {"source": "s"}, "kind": "equations", "name": "epidemic",
            "carrier": {"kind": "euclid", "variables": ["x"]}, "parameters": {"R0": R0},
            "drift": {"x": "R0*x*(1 - x) - x"}, "noise": 0.01,
            "rates": {"birth": "R0*x*(1 - x)", "death": "x", "upper": "1", "scale": "hosts"}}


def test_budworm_exponents_and_exact_times_agree_with_the_agent_comparison():
    report = population.escape(spec.load(FOLDER / "spruce_budworm.json"), N=100)
    up, down = path(report, 0.8386, 1.769), path(report, 5.392, 1.769)
    assert up["rates"] == pytest.approx(0.02586, abs=1e-4) and down["rates"] == pytest.approx(0.20325, abs=1e-4)
    assert up["ln_time_exact"] == pytest.approx(6.3159, abs=1e-3) and down["ln_time_exact"] == pytest.approx(22.905, abs=0.01)
    assert report["equal_occupation"]["rates"] == pytest.approx(7.130, abs=0.002)
    assert report["equal_occupation"]["constant_noise"] == pytest.approx(6.947, abs=0.002)


@pytest.mark.parametrize("R0", [1.5, 10.0, 30.0])
def test_the_rates_give_the_slope_of_the_exact_chain_where_the_state_noise_does_not(R0):
    real = spec.load(epidemic(R0))
    x = 1 - 1 / R0
    b, d, upper = population._functions(real, None)
    Ns = (400, 800) if R0 < 2 else (100, 200)
    times = [population.log_exit_time(b, d, upper, N, x, 0.0) for N in Ns]
    slope = (times[1] - times[0]) / (Ns[1] - Ns[0])
    row = population.escape(real)["paths"][0]
    assert row["rates"] == pytest.approx(slope, rel=0.02)
    if R0 > 5:
        assert abs(row["state_noise"] / row["rates"] - 1) > 0.15


def test_the_exponent_from_the_rates_is_independent_of_the_unit_of_time_and_the_constant_noise_is_not():
    plain = epidemic(3.0)
    scaled = json.loads(json.dumps(plain))   # time rescaled by a function of the state
    scaled["drift"]["x"] = "(R0*x*(1 - x) - x)*(1 + x)"
    scaled["rates"].update(birth="R0*x*(1 - x)*(1 + x)", death="x*(1 + x)")
    a, b = (population.escape(spec.load(s))["paths"][0] for s in (plain, scaled))
    assert a["rates"] == pytest.approx(b["rates"], rel=1e-6) and a["state_noise"] == pytest.approx(b["state_noise"], rel=1e-6)
    assert abs(a["constant_noise"] / b["constant_noise"] - 1) > 0.1


def test_rates_must_add_up_to_the_drift_and_need_one_variable():
    wrong = epidemic(2.0)
    wrong["rates"]["death"] = "2*x"
    with pytest.raises(spec.SpecError, match="must equal the drift"):
        spec.load(wrong)
    two = epidemic(2.0)
    two["carrier"]["variables"] = ["x", "y"]
    two["drift"]["y"] = "-y"
    with pytest.raises(spec.SpecError, match="one variable"):
        spec.load(two)


def test_a_model_without_rates_is_refused():
    with pytest.raises(ValueError, match="no birth and death rates"):
        population.escape(spec.load(FOLDER / "kuramoto_population.json"))
