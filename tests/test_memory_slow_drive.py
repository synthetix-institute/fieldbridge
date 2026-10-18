"""Return to a turning point of a material given by equations (fieldbridge.memory.slow_drive).

Two bistable coordinates, each in its own quartic potential, are coupled linearly and driven by the same slow field.
Their thresholds differ by a small field of opposite sign on each. The transition graph of the equations, the exact
return test on it and the reduced hysteron model are compared with what the signs of the coupling require: a
cooperative coupling reduces to a network balanced with the drive, an antiferromagnetic one to a network with a
frustrated loop through the drive. Identical elements with an antiferromagnetic coupling lose stability by a choice
that the drive does not make, and the construction must say so instead of returning a transition graph.
"""
import pytest

np = pytest.importorskip("numpy")
pytest.importorskip("sympy")

from fieldbridge.memory import slow_drive as sd
from fieldbridge.memory import spec as spec_mod
from fieldbridge.memory.spec import SpecError

FAST = {"steps": 200, "n_H": 8, "n_starts": 16}


def pair(J, a=0.05, H_range=(-1.5, 1.5)):
    return {"schema": "fieldbridge-memory/1", "kind": "equations", "name": f"pair J={J}", "question": "test",
            "assumptions": [], "provenance": "test", "noise": 0.01,
            "carrier": {"kind": "euclid", "variables": ["x1", "x2"], "scale": 1.0},
            "parameters": {"J": J, "H": 0.0, "a": a},
            "drift": {"x1": "x1 - x1**3 + J*x2 + H + a", "x2": "x2 - x2**3 + J*x1 + H - a"},
            "slow_drive": {"parameter": "H", "range": list(H_range), "elements": ["x1", "x2"]}}


def run(spec):
    return sd.check(spec, spec_mod.load(spec), **FAST)


def test_cooperative_pair_reduces_to_a_balanced_network():
    r = run(pair(0.1))
    J = np.array(r["reduced"]["couplings"])
    assert J[0, 1] > 0 and J[1, 0] > 0
    assert r["prediction"]["guaranteed"]
    assert r["return"]["failed"] == 0
    # one avalanche switches both elements: two states, and the reduced model has the same two
    assert r["count"] == 2 and r["reduced_measure"]["count"] == 2
    assert r["agreement"]["count"] and r["agreement"]["return"]


def test_antiferromagnetic_pair_has_a_frustrated_loop_through_the_drive():
    r = run(pair(-0.1))
    J = np.array(r["reduced"]["couplings"])
    assert J[0, 1] < 0 and J[1, 0] < 0
    assert not r["prediction"]["guaranteed"]
    # the elements switch one at a time: --, +-, ++; the state -+ exists but no history reaches it
    assert [s["signs"] for s in r["states"]] == ["--", "+-", "++"]
    assert r["stable_states_not_reached"] == ["-+"]
    assert r["return"]["failed"] == 0
    assert r["agreement"]["count"] and r["agreement"]["return"]
    # the fields of opposite sign on the two elements are recovered, in units of the drive
    h = np.array(r["reduced"]["fields"])
    assert h[0] > 0 > h[1] and abs(h[0] + h[1]) < 1e-6


def test_identical_antiferromagnetic_elements_are_a_choice_the_drive_does_not_make():
    with pytest.raises(sd.ReductionError, match="does not decide"):
        run(pair(-0.1, a=0.0))


def test_the_range_must_saturate_the_material():
    with pytest.raises(SpecError, match="does not saturate"):
        run(pair(0.1, H_range=(-0.1, 0.1)))


def test_the_drive_must_be_a_declared_parameter():
    spec = pair(0.1)
    spec["slow_drive"]["parameter"] = "B"
    with pytest.raises(SpecError, match="declared parameter"):
        run(spec)
