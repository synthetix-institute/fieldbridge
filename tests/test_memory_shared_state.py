"""Two copies exchanging one variable: the limit of their difference and the part held by the shared state."""
from pathlib import Path

from fieldbridge.memory import shared_state as ss
from fieldbridge.memory.spec import load

EX = Path(__file__).resolve().parents[1] / "examples" / "memory" / "shared_state"


def test_limit_through_the_shared_state_only():
    real = load(EX / "second_variable_limit.json")
    b = ss.balance(real, "y", 3.0, 1.0)
    assert abs(b["g_own"]) < 1e-6 * abs(b["g"]) and b["g_shared"] < 0
    free = ss.integrate(real, "y", 3.0, 1.0, b["equal_state"])
    assert free["endpoint"] == "rest" and abs(free["size"] / b["free_difference"] - 1) < 0.05
    assert ss.integrate(real, "y", 3.0, 1.0, b["equal_state"], hold=True)["endpoint"] == "runs away"


def test_shared_state_against_the_own_limit():
    real = load(EX / "fitzhugh_nagumo_units.json")
    b = ss.balance(real, "v", 3.0, 0.6125)
    assert b["g_own"] < 0 < b["g_shared"] and b["free_valid"]
    held = ss.integrate(real, "v", 3.0, 0.6125, b["equal_state"], hold=True)
    assert held["endpoint"] == "rest" and abs(held["size"] / b["held_difference"] - 1) < 0.1
    free = ss.integrate(real, "v", 3.0, 0.6125, b["equal_state"])
    assert free["endpoint"] == "rest" and free["size"] > held["size"]
