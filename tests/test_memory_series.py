"""The order-by-order reduction at a write point (analysis.series_reduction, memory/series.py) and its use by
analysis.normal_form. The window fit biases a3 at a pitchfork by up to 29 % (ring of 4: -0.385 against -0.541) and a2 at a
fold by up to 24 %, and it does not resolve the fifth-order term of a parametric oscillator whose rate is small against
its loss, which it recorded as a line of states."""
from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import pytest
import sympy as sp

from fieldbridge.memory import analysis as an
from fieldbridge.memory import codiscovery, construct, library
from fieldbridge.memory import spec as mspec
from fieldbridge.memory.series import NotAnalytic, Series

ROOT = Path(__file__).resolve().parents[1]
QUANTUM = "examples/quantum/semiclassical/"


def load(path):
    return mspec.load(json.loads((ROOT / path).read_text()))


def write_points(real, rng=0):
    return an.locate_writes(real, np.random.default_rng(rng), n_starts=24 if real.carrier.dim <= 4 else 12)


def reduced(real, value, rng=0):
    ev = min(write_points(real, rng), key=lambda e: abs(e["v"] - value))
    return an.normal_form(real, ev["q"], real.control, ev["v"])


def equations(drift, parameters, control, scale=1.0):
    variables = list(drift)
    return mspec.load({"schema": "fieldbridge-memory/1", "kind": "equations", "name": "test", "question": "test",
                       "assumptions": [], "provenance": {"source": "test"},
                       "carrier": {"kind": "euclid", "variables": variables, "scale": scale},
                       "parameters": parameters, "drift": drift, "noise": 0.01,
                       "control": {"name": control, "range": [-1.0, 1.0]}})


@pytest.mark.parametrize("name, value", [("sqrt", 2.0), ("exp", 0.3), ("log", 1.7), ("sin", 0.6), ("cos", 0.6),
                                         ("tanh", 0.4), ("abs", -1.3)])
def test_the_series_of_every_function_of_a_specification_is_its_taylor_series(name, value):
    assert name in mspec.FUNCTIONS
    s = sp.Symbol("s", real=True)
    fn = sp.Abs if name == "abs" else getattr(sp, name)
    expected = sp.series(fn(value + 2.0 * s + 0.5 * s ** 2), s, 0, 8).removeO()
    x = Series([value, 2.0, 0.5] + [0.0] * 5)
    got = abs(x) if name == "abs" else getattr(x, name)()
    assert got.c == pytest.approx([float(expected.coeff(s, k)) for k in range(8)], rel=1e-10, abs=1e-12)


def test_the_operators_and_powers_of_a_series():
    s = sp.Symbol("s", real=True)
    y = 1.5 + s - 0.25 * s ** 3
    x = Series([1.5, 1.0, 0.0, -0.25, 0.0, 0.0])
    for got, form in [(x * x - 3 * x + 2.0, y * y - 3 * y + 2), (1.0 / x, 1 / y), (2.0 - x / 3.0, 2 - y / 3),
                      (x ** 3, y ** 3), (x ** 2.5, y ** 2.5), (x ** -2.0, y ** -2), (2.0 ** x, sp.exp(y * sp.log(2))), (x ** x, sp.exp(y * sp.log(y)))]:
        expected = [float(sp.diff(form, s, k).subs(s, 0)) / math.factorial(k) for k in range(6)]
        assert got.c == pytest.approx(expected, rel=1e-10, abs=1e-12)


@pytest.mark.parametrize("op", [lambda x: abs(x), lambda x: x.sqrt(), lambda x: x.log(), lambda x: x ** 0.5, lambda x: 1.0 / x])
def test_a_function_that_is_not_analytic_at_zero_is_refused(op):
    with pytest.raises(NotAnalytic):
        op(Series([0.0, 1.0, 0.0]))


def test_a_kink_is_reported_as_not_covered_and_the_fit_is_kept():
    real = equations({"x": "a*x - abs(x)*x"}, {"a": 0.0}, "a")
    nf = an.normal_form(real, np.array([0.0]), "a", 0.0)
    assert not nf["series"]["covered"] and "abs at a kink" in nf["series"]["reason"]
    assert nf["a2"] == nf["fit"]["a2"] and "a4" not in nf


@pytest.mark.parametrize("drift, expected", [
    ("a*x - x**3", {2: 0.0, 3: -1.0, 4: 0.0, 5: 0.0}),
    ("a + x**2 - x**3 + 2*x**5", {2: 1.0, 3: -1.0, 4: 0.0, 5: 2.0}),
    ("a*x + x**2*sqrt(1 + x)", {2: 1.0, 3: 0.5, 4: -0.125, 5: 0.0625}),
    ("a*x + exp(x) - 1 - x", {2: 0.5, 3: 1 / 6, 4: 1 / 24, 5: 1 / 120}),
    ("a*x - x*log(1 + x**2)", {2: 0.0, 3: -1.0, 4: 0.0, 5: 0.5}),
    ("a*x - tanh(x)**3", {2: 0.0, 3: -1.0, 4: 0.0, 5: 1.0}),
    ("a*x + sin(x) - x", {2: 0.0, 3: -1 / 6, 4: 0.0, 5: 1 / 120}),
    ("a*x + cos(x) - 1", {2: -0.5, 3: 0.0, 4: 1 / 24, 5: 0.0}),
])
def test_the_reduction_of_one_coordinate_is_the_taylor_series_of_its_drift(drift, expected):
    real = equations({"x": drift}, {"a": 0.0}, "a")
    a = an.series_reduction(real, np.array([0.0]), np.array([1.0]), np.array([1.0]), {"a": 0.0})["coefficients"]
    for n, c in expected.items():
        assert a[n] == pytest.approx(c, rel=1e-9, abs=1e-12)


def test_the_centre_manifold_is_closed_order_by_order():
    """x' = a x - x z, z' = -L (z - x^2): on z = x^2 + 2 x^4 / L + ... (a = 0), x' = -x^3 - (2 / L) x^5 + ..., checked
    down to a transverse rate 1e7 times the critical one, where a rank decision on the Jacobian would drop the mode."""
    for rate in (1.0, 1e7):
        real = equations({"x": "a*x - x*z", "z": "-L*(z - x**2)"}, {"a": 0.0, "L": rate}, "a")
        a = an.series_reduction(real, np.zeros(2), np.array([1.0, 0.0]), np.array([1.0, 0.0]), {"a": 0.0})["coefficients"]
        assert a[3] == pytest.approx(-1.0, rel=1e-12) and a[5] == pytest.approx(-2.0 / rate, rel=1e-9)
        assert a[2] == 0.0 and a[4] == 0.0


@pytest.mark.parametrize("path, value, a3", [
    ("examples/memory/laser.json", 1.0, -0.5),
    ("examples/memory/toggle.json", 2.0, -0.125),
    ("examples/memory/stoner_wohlfarth.json", 1.0, -0.5),
    ("examples/memory/repressor_ring4.json", 1.013, -0.5413),
    ("examples/memory/convection_roll.json", 1.0, -0.1641),
])
def test_the_cubic_coefficient_of_a_pitchfork_has_no_window_bias(path, value, a3):
    nf = reduced(load(path), value)
    assert nf["series"]["covered"] and nf["kind"].startswith("supercritical pitchfork")
    assert nf["a3"] == pytest.approx(a3, rel=2e-4 if a3 == -0.5413 else 2e-3)
    assert abs(nf["a2"]) < 1e-7 and nf["asymmetry"] < 1e-6


@pytest.mark.parametrize("path, value", [(QUANTUM + "magnon_kerr.json", 1.112), (QUANTUM + "transmon_driven.json", 0.6102),
                                         ("examples/memory/schlogl.json", 1.49)])
def test_the_quadratic_coefficient_of_a_fold_is_half_the_second_derivative_of_the_drift(path, value):
    """a2 is the same in every coordinate with this tangent, so that it can be checked against the drift itself; a3 is
    not, and it is not compared with anything at a fold."""
    real = load(path)
    nf = reduced(real, value)
    assert nf["kind"].startswith("saddle-node")
    q, vec, w = (np.asarray(nf[k]) for k in ("state", "mode", "left_mode"))
    over = {real.control: nf["value"]}

    def curvature(h):
        return 0.5 * float(w @ (real.F(q + h * vec, **over) + real.F(q - h * vec, **over) - 2 * real.F(q, **over))) / h ** 2

    assert nf["a2"] == pytest.approx((4 * curvature(5e-4) - curvature(1e-3)) / 3, rel=1e-4)


def test_stored_values_of_the_folds_named_in_the_records():
    assert reduced(load(QUANTUM + "magnon_kerr.json"), 1.112)["a2"] == pytest.approx(0.36823, rel=1e-4)
    assert reduced(load(QUANTUM + "transmon_driven.json"), 0.6102)["a2"] == pytest.approx(-1.18684, rel=1e-4)


def test_a_fold_with_large_higher_terms_is_read_as_a_fold():
    """At a fold of the fluxonium, a3 = -19.6 and a7 = -5e7 in the coordinate of the series: its radius is 0.03, far
    below the window of the carrier (1.5), and over that window a3 s^3 is not a measure of the drift. The comparison of
    the quadratic with the cubic term, and of the bias with the growth of the rate, is made at the reach of the series:
    a state appears at the point, where the window read two states exchanging stability."""
    nf = reduced(load(QUANTUM + "fluxonium.json"), 5.701)
    assert nf["a2"] == pytest.approx(1.6908, rel=1e-3) and nf["half"] > 0.5
    assert an.series_reach(nf["series"], nf["half"], nf["second_rate"]) < 0.1
    assert nf["asymmetry"] > construct.TOL and nf["bias"] > construct.TOL and nf["kind"].startswith("saddle-node")


def test_the_reduced_law_is_withdrawn_where_the_eliminated_voltage_is_not_fast():
    """Fluxonium is a phase with inertia, ph' = v, v' = -0.3 v - U'(ph). Eliminating the voltage gives the overdamped drift
    -U'/0.3, which at the written states varies faster than the damping it assumes (ratio about 1.8, above
    REDUCED_SLOW). The one-coordinate law, which there reads 1.09 against 1.32 of the full dynamics, is withdrawn, and the
    linear stage with its outcome map is reported and certified."""
    from fieldbridge.memory import completion as cm
    real = load(QUANTUM + "fluxonium.json")
    entry = cm.complete(real, np.random.default_rng(1), n_samples=2000, n_traj=2000)["write_points"][0]
    law = next(u["law"] for u in entry["unfoldings"] if "reduction_validity" in u.get("law", {}))
    assert not law["reduction_validity"]["slow"] and law["reduction_validity"]["ratio"] > 1.0
    assert law["law_kind"] == "linear stage and outcome map" and "reduction_not_slow" in law
    assert law["certified"] and law["write_class"] == "linear-stage write"


@pytest.mark.parametrize("path, value, coefficient", [(QUANTUM + "kerr_cat_kpo.json", 0.5, -25.0),
                                                      (QUANTUM + "kpo_pair.json", 0.5, -6.25),
                                                      (QUANTUM + "pplo_oscillator.json", 1.0, -5.76e-4)])
def test_the_leading_term_of_a_vanishing_cubic(path, value, coefficient):
    nf = reduced(load(path), value)
    assert nf["kind"].startswith("degenerate cubic")
    assert nf["leading_term"]["order"] == 5
    assert nf["leading_term"]["coefficient"] == pytest.approx(coefficient, rel=1e-5)


def test_the_phase_locked_oscillator_is_a_vanishing_cubic_and_not_a_line_of_states():
    """The Kerr rate K = 12 gnl = 0.024 rotates the amplitude, as in a KPO at zero detuning, so that the cubic term
    vanishes and the drift saturates at -K^2 / kappa s^5, too small over the default window for the amplitudes of
    analysis.leading_order to resolve."""
    real = load(QUANTUM + "pplo_oscillator.json")
    nf = reduced(real, 1.0)
    K, kappa = 12 * real.params["gnl"], real.params["kappa"]
    assert nf["leading_term"]["coefficient"] == pytest.approx(-K ** 2 / kappa, rel=1e-5)
    ob = construct.obstruction(nf)
    assert ob["degenerate_cubic"] and "vanishing cubic (leading term of order 5)" in ob["terms"]
    row = codiscovery.derive_symmetric_write(real, np.random.default_rng(0), check_law=False)
    assert row["obstruction_short"] == "vanishing cubic: leading term of order 5"


def test_a_drift_that_vanishes_along_the_critical_direction_is_a_line_of_states():
    real = equations({"x": "a*x - a*x", "y": "-y + x*y"}, {"a": 0.0}, "a")
    a = an.series_reduction(real, np.zeros(2), np.array([1.0, 0.0]), np.array([1.0, 0.0]), {"a": 0.0})
    assert a["covered"] and max(abs(c) for c in a["coefficients"]) < 1e-15
    assert an.leading_series_term(a, 0.25, 1.0) is None


def test_several_critical_directions_are_not_covered():
    nf = reduced(load(QUANTUM + "kpo_triangle.json"), 0.5)
    assert nf["critical_dimension"] == 2 and not nf["series"]["covered"]
    assert "2 critical directions" in nf["series"]["reason"] and nf["a3"] == nf["fit"]["a3"]


def test_the_series_runs_through_the_realizations_that_wrap_a_drift():
    """The vectorized drifts of the library, the write field, and a combined sweep take an object array of Series."""
    ring = library.bistable_ring(4)
    assert an.normal_form(ring, np.zeros(4), "eps", 0.0)["series"]["covered"]
    capillary = library.rotor_patch("capillary")
    a = an.series_reduction(capillary, np.ones(12), np.eye(12)[0], np.eye(12)[0], {})
    assert a["covered"], a.get("reason")
    toggle = library.toggle()
    q = np.array([1.0, 1.0])
    forced = an.with_write_field(toggle, np.array([3.0, 0.5]), (0.0, 1.0))
    a = an.series_reduction(forced, q, np.array([1.0, 0.0]), np.array([1.0, 0.0]), {an.WRITE_FIELD: 0.5})
    assert a["covered"], a.get("reason")
    path = construct.along_path(toggle, {"alpha": 2.0, "n": 2.0}, {"alpha": 1.0, "n": 0.0})
    assert an.normal_form(path, q, "t", 0.0)["series"]["covered"]
    torus = an.with_write_field(capillary, np.zeros(12), (0.0, 1.0))
    a = an.series_reduction(torus, np.ones(12), np.eye(12)[0], np.eye(12)[0], {an.WRITE_FIELD: 0.5})
    assert a["covered"], a.get("reason")


def test_the_cusp_of_the_bifurcation_readout_needs_no_narrowed_window():
    """At the cusp of the JBA the series has a2 and a4 about 1e-11 in the coordinate s = w.(x - q), and the sweep reads as
    a pitchfork at the default window. The cubic fitted in the coordinate along the mode, kept beside it, reads the even
    term of the drift in that coordinate as a quadratic coefficient of order 1e-3."""
    real = load(QUANTUM + "jba_resonator.json")
    ev = write_points(real)[0]
    nf = an.normal_form(real, ev["q"], real.control, ev["v"])
    c = construct.cancel_asymmetry(real, nf, "chi")
    assert c["success"] and c["normal_form_along_sweep"]["kind"].startswith("supercritical pitchfork")
    path = construct.along_path(real, {real.control: c["control"], "chi": c["value"]}, c["direction"])
    sweep = an.normal_form(path, np.asarray(c["state"]), "t", 0.0)
    assert sweep["kind"].startswith("supercritical pitchfork")
    assert sweep["a3"] == pytest.approx(-0.057735, rel=1e-4) and abs(sweep["a2"]) < 1e-9 and abs(sweep["a4"]) < 1e-9
    assert abs(sweep["fit"]["a2"]) > 1e-4


def test_the_reduction_is_fast_enough_to_run_at_every_write_point():
    """Not a timing bound on a loaded machine: the number of evaluations of the drift. One per coordinate at first order,
    and one per order, for the critical direction."""
    real = library.bistable_ring(8)
    calls = []
    base = real.drift

    def counted(q, p):
        calls.append(q.dtype)
        return base(q, p)

    real.drift = counted
    an.series_reduction(real, np.zeros(8), np.ones(8) / math.sqrt(8), np.ones(8) / math.sqrt(8), {"eps": 0.0})
    assert len(calls) == 8 + an.SERIES_ORDER
