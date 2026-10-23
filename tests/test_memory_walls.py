"""Chains of bistable elements (fieldbridge.memory.walls): the carriers, their stable patterns and the kinetic chain.

The Kerr-oscillator chain is Hamiltonian plus a uniform contraction, so the eigenvalues at every equilibrium pair
about -1/2; an isolated site of a ferromagnetically coupled chain stops being a stable state at a finite coupling,
while a wall stays pinned up to nearly the coercive field of one element. The kinetic chain with Glauber's rates must
reproduce Glauber's exact overlap, and the first error of a word is the sum of the rates of its sites.
"""
import pytest

np = pytest.importorskip("numpy")
pytest.importorskip("scipy")

from fieldbridge.memory import walls as w
from fieldbridge.memory.spec import SpecError


def _numerical_jacobian(ch, z):
    x = ch.to_real(z)
    F = np.zeros((x.size, x.size))
    for k in range(x.size):
        e = np.zeros(x.size)
        e[k] = 1e-6
        F[:, k] = (ch.to_real(ch.drift(ch.from_real(x + e))) - ch.to_real(ch.drift(ch.from_real(x - e)))) / 2e-6
    return F


@pytest.mark.parametrize("ch", [w.KPOChain(N=5, J=0.3, E=0.2, Delta=0.1), w.KPOChain(N=2, J=0.3, ring=False),
                                w.KPOChain(N=4, J=0.2, ring=False, ordered=True),
                                w.RFSquidChain(N=6, J=0.1, theta_x=0.05), w.RFSquidChain(N=2, J=0.2, ring=False)])
def test_the_jacobians_are_the_derivatives_of_the_drifts(ch):
    rng = np.random.default_rng(0)
    z = ch.start(rng.choice([-1, 1], ch.N))
    z = z + 0.1 * (rng.standard_normal(ch.N) + (1j * rng.standard_normal(ch.N) if ch.complex_state else 0))
    assert np.abs(ch.jacobian(z) - _numerical_jacobian(ch, z)).max() < 1e-7


def test_the_kerr_chain_is_hamiltonian_plus_a_uniform_contraction():
    ch = w.KPOChain(N=8, J=0.4, E=0.1)
    fp = w.fixed_point(ch, [1, 1, 1, -1, -1, 1, -1, -1])
    ev = fp["eigenvalues"]
    # eigenvalues come in pairs lambda, -1 - conj(lambda) (real parts symmetric about -kappa/2 = -1/2)
    mirrored = -1 - np.conj(ev)
    assert np.abs(mirrored[:, None] - ev[None, :]).min(1).max() < 1e-8
    assert abs(np.mean(ev.real) + 0.5) < 1e-10


def test_an_isolated_site_is_lost_at_a_finite_coupling_and_a_wall_is_not():
    ch = w.KPOChain(N=12, p=4.0, J=0.0)
    iso = w.loss_threshold(ch, [1] * 11 + [-1], "J", 0.0, 1.0)
    assert iso["reached"] and abs(iso["value"] - 0.728) < 0.005
    wall = w.loss_threshold(ch, [1] * 6 + [-1] * 6, "J", 0.0, 2.0)
    assert not wall["reached"]


def test_a_wall_is_pinned_up_to_nearly_the_coercive_field_of_one_element():
    single = w.loss_threshold(w.KPOChain(N=1, p=4.0, J=0.0), [-1], "E", 0.0, 4.0)
    wall = w.loss_threshold(w.KPOChain(N=12, p=4.0, J=0.25), [1] * 6 + [-1] * 6, "E", 0.0, 4.0)
    assert single["reached"] and wall["reached"]
    assert 0.9 * single["value"] < wall["value"] < single["value"]
    assert w.wall_velocity(w.KPOChain(N=24, p=4.0, J=0.25), 0.5 * wall["value"], t_max=60.0)["pinned"]


def test_classes_are_canonical_under_reflection_and_reversal():
    assert w.class_name([1, 1, -1, 1, -1], True) == w.class_name([-1, 1, -1, 1, 1], True)
    assert w.class_name([1, 1, -1, 1, -1], True) == w.class_name([-1, -1, 1, -1, 1], True)
    assert w.class_name([1, 1, -1, 1, -1], False) != w.class_name([-1, -1, 1, -1, 1], False)
    assert w.nn_name([1, 1, -1, 1, -1], True) == "m2"
    assert w.nn_name([0, 0, 1, -1, 1], False) == "m1e+"


def _rates(g0, g1, g2):
    return {"five_site": {}, "nearest": {c: {"rate": g, "events": 100} for c, g in (("m0", g0), ("m1", g1), ("m2", g2))}}


def test_the_kinetic_chain_with_glauber_rates_has_glauber_overlap():
    g0, g1 = 0.01, 0.1
    kc = w.KineticChain(400, True, _rates(g0, g1, 2 * g1 - g0), True, use_five=False)
    rng = np.random.default_rng(1)
    times = np.array([0.0, 2.0, 5.0, 10.0, 20.0])
    s = rng.choice([-1, 1], 400)
    out = kc.evolve(s, times, 40, rng)
    exact = w.glauber_overlap(times, g1, g0)
    assert np.allclose(out["overlap"], exact, atol=0.02)


def test_the_first_error_of_a_word_is_at_the_sum_of_the_rates_of_its_sites():
    kc = w.KineticChain(8, True, _rates(0.001, 0.01, 0.1), True, use_five=False)
    s = [1, 1, 1, -1, -1, 1, -1, -1]          # site 5 is isolated (m = 2), four sites sit at walls (m = 1)
    expected = sum({0: 0.001, 1: 0.01, 2: 0.1}[int(s[j - 1] != s[j]) + int(s[(j + 1) % 8] != s[j])] for j in range(8))
    assert w.word_rate(kc, s) == pytest.approx(expected)
    rng = np.random.default_rng(2)
    out = kc.evolve(s, [0.0, 5.0], 20000, rng)
    assert out["no_switch"][1] == pytest.approx(np.exp(-expected * 5.0), abs=0.01)


def test_the_gradient_chain_has_symmetric_langer_rates_at_a_wall():
    ch = w.RFSquidChain(N=8, J=0.05, T=0.1)
    s = [1, 1, 1, 1, -1, -1, -1, -1]
    left, right = w.langer_rate(ch, s, 3), w.langer_rate(ch, s, 4)
    single = w.langer_rate(w.RFSquidChain(N=8, J=0.0, T=0.1), s, 3)
    assert left["rate"] == pytest.approx(right["rate"], rel=1e-8)
    # a wall site feels no net field from its neighbours: its barrier is that of a free loop to first order in J
    assert abs(left["barrier"] - single["barrier"]) < 0.02 * single["barrier"]


def test_a_specification_builds_a_chain_and_refuses_unknown_entries():
    spec = {"schema": "fieldbridge-memory/1", "kind": "bistable-chain", "carrier": "kpo", "question": "q",
            "assumptions": [], "provenance": "test", "chain": {"N": 8, "J": 0.2, "K": 0.5}}
    assert isinstance(w.from_spec(spec), w.KPOChain)
    with pytest.raises(SpecError, match="unknown chain entries"):
        w.from_spec({**spec, "chain": {"N": 8, "g": 1.0}})
    with pytest.raises(SpecError, match="carrier"):
        w.from_spec({**spec, "carrier": "spin"})
