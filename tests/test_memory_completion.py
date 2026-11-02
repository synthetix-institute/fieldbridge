"""The law at an obstruction (memory.completion): the completion that each obstruction names, the law derived from the
sweep, and its certificate by the full stochastic dynamics.

Statistical thresholds only (the probabilities are Monte Carlo estimates on a few thousand samples); no solver
residual is pinned."""
from __future__ import annotations

import json
from math import pi
from pathlib import Path

import numpy as np
import pytest

from fieldbridge.memory import analysis as an, completion, construct, spec

EX = Path(__file__).resolve().parents[1] / "examples" / "memory"


def _law(name: str, **kw):
    real = spec.load(EX / f"{name}.json")
    rep = completion.complete(real, np.random.default_rng(3), floors=(1.0,), **kw)
    return rep["write_points"]


def test_completion_kind_follows_the_obstruction():
    gen = {"kind": "supercritical pitchfork", "critical_dimension": 1}
    assert completion.completion_kind(gen, {"realizes_generator": True}) == "none"
    assert completion.completion_kind({"critical_dimension": 2}, {"degenerate": True}) == "critical space"
    assert completion.completion_kind({"critical_dimension": 1}, {"degenerate_cubic": True}) == "leading order"
    assert completion.completion_kind({"critical_dimension": 1}, {"asymmetric": True}) == "second parameter"
    assert completion.completion_kind({"critical_dimension": 1}, {"subcritical": True}) == "distant branch"
    assert completion.completion_kind({"hopf": True}, {"hopf": True}) == "phase"


def test_diverged_trajectories_are_not_assigned_to_a_state():
    real = spec.load(EX / "pitchfork.json")
    stable = [np.array([1.0]), np.array([-1.0])]
    finals = np.array([[0.9], [-1.1], [np.nan], [np.inf]])
    lab = completion.assign(real, finals, stable)
    assert lab.tolist() == [0, 1, -1, -1]
    assert completion._probs(lab, 2) == [0.5, 0.5]


def test_the_rung_derives_the_swept_pitchfork_law_where_nothing_obstructs():
    """Control: the law derived from the sweep alone (no constant put in) is the transferred one, pi^(1/4), with no
    offset, and the full dynamics certifies it."""
    pts = _law("pitchfork", n_samples=3000, n_traj=3000)
    assert len(pts) == 1 and pts[0]["completion"] == "none"
    L = pts[0]["law"]
    assert L["status"] == "derived" and L["certified"]
    for src in ("law", "certificate"):
        x = L["one_direction"][src][0]
        assert abs(x["response_constant"] - pi ** 0.25) < 0.12
        assert abs(x["offset"]) < 0.1


def test_a_subcritical_write_point_gets_a_certified_law():
    """Two equal tubes: a positive cubic used to stop the derivation; the decision is still made in the linear
    stage, and the law of the completed model (with the distant branch) is certified."""
    pts = _law("tubes", n_samples=3000, n_traj=3000)
    sub = [p for p in pts if p["completion"] == "distant branch"]
    assert sub, [p["kind"] for p in pts]
    L = sub[0]["law"]
    assert L["status"] == "derived" and L["certified"]
    assert abs(L["one_direction"]["certificate"][0]["response_constant"] - pi ** 0.25) < 0.15


def test_two_critical_directions_give_a_law_over_all_written_states():
    """A square-symmetric two-component model has two critical directions; the law over its four written states,
    for biases in eight directions of the critical plane, agrees with the full dynamics."""
    pts = _law("two_component_square", n_samples=2000, n_traj=2000)
    deg = [p for p in pts if p["completion"] == "critical space"]
    assert deg, [p["kind"] for p in pts]
    L = deg[0]["law"]
    assert L["critical_dimension"] == 2 and len(L["stable_states"]) == 4 and len(L["rows"]) == 8
    for row in L["rows"]:
        assert abs(sum(row["law"]) - 1) < 1e-9
        assert row["deviation"] < 0.06
    # a bias favours the states on its side of the plane: the most probable state changes with the direction
    assert len({int(np.argmax(r["law"])) for r in L["rows"]}) == 4


def test_a_degenerate_critical_pair_keeps_a_left_vector_dual_to_its_mode():
    """The write point that the default seed of memory complete finds on the square model: the Jacobian there is eps
    times the identity to rounding, and the eigensolver returned the left eigenvectors of the pair in the other order,
    so that the left vector picked was orthogonal to the mode. The coefficients became NaN and the point read as a
    subcritical pitchfork with no critical direction; the left vector is now the dual one."""
    real = spec.load(EX / "two_component_square.json")
    q = np.array([4.588679728723345e-19, -4.58914073268652e-19])
    nf = an.normal_form(real, q, "eps", 1.000006995815892e-07)
    assert nf["critical_dimension"] == 2 and nf["kind"].startswith("degenerate critical mode")
    assert np.all(np.isfinite(nf["left_mode"])) and np.isfinite(nf["a3"])
    assert float(np.dot(nf["left_mode"], nf["mode"])) == pytest.approx(1.0, abs=1e-9)
    assert completion.completion_kind(nf, construct.obstruction(nf)) == "critical space"


def test_a_choice_the_noise_undoes_during_the_hold_is_an_equilibrium_write():
    """The pitchfork swept at the standard rate with a noise far above the depth of the wells at the end of the sweep:
    the linear-stage law predicts the usual response, the full dynamics loses the choice to hopping, and the rung names
    the equilibrium write by the trajectories that change state during a hold at the end value of the control."""
    real = spec.load(EX / "pitchfork.json")
    rep = completion.complete(real, np.random.default_rng(3), floors=(1.0,), noise=0.5, n_samples=400, n_traj=300)
    laws = [e["law"] for e in rep["write_points"] if "law" in e and e["law"].get("status") == "derived"]
    assert laws, rep["write_points"]
    L = laws[0]
    assert L["certified"] is False and L["write_class"] == "equilibrium write"
    assert L["hold"]["switching_fraction"] > completion.SWITCHING_EQUILIBRIUM
    assert "equilibrium write" in L["obstruction_in_nonlinear_stage"]
    # the standard noise rule keeps the pitchfork's write: the class is then the linear-stage write
    rep = completion.complete(real, np.random.default_rng(3), floors=(1.0,), n_samples=400, n_traj=300)
    L = next(e["law"] for e in rep["write_points"] if "law" in e and e["law"].get("status") == "derived")
    assert L["write_class"] == "linear-stage write" and "hold" not in L


def test_the_sweep_through_the_kerr_cusp_reads_as_a_pitchfork_on_the_narrowed_window():
    """At the cusp of a driven Kerr resonator the quartic term of the drift appears as a quadratic one on the default
    window, which reads as transcritical and would refuse the law; the completion reads the sweep on the narrowed
    window of construct.sweep_normal_form, where it is a pitchfork."""
    import json
    from pathlib import Path
    from fieldbridge.memory import analysis as an, construct, spec as mspec
    root = Path(__file__).resolve().parents[1]
    real = mspec.load(json.loads((root / "examples/quantum/semiclassical/jba_resonator.json").read_text()))
    ev = an.locate_writes(real, np.random.default_rng(0), n_starts=24)[0]
    nf = an.normal_form(real, ev["q"], real.control, ev["v"])
    c = construct.cancel_asymmetry(real, nf, "chi")
    assert c.get("success")
    path = construct.along_path(real, {real.control: c["control"], "chi": c["value"]}, c["direction"])
    assert "transcritical" in an.normal_form(path, np.asarray(c["state"]), "t", 0.0)["kind"]
    narrowed = construct.sweep_normal_form(path, np.asarray(c["state"]))
    assert "pitchfork" in narrowed["kind"] and narrowed["window_halvings"] >= 1


def test_the_curved_sweep_moves_the_parameters_along_a_parabola():
    """along_curved_path: at sweep parameter t the two parameters sit at origin + t d + c2 t^2 n, with n the direction
    d turned by a right angle in their plane."""
    import json
    from pathlib import Path
    from fieldbridge.memory import construct, spec as mspec
    root = Path(__file__).resolve().parents[1]
    real = mspec.load(json.loads((root / "examples/quantum/semiclassical/jba_resonator.json").read_text()))
    seen = {}
    base = real.drift

    def spy(q, p):
        seen.update(F=p["F"], chi=p["chi"])
        return base(q, p)

    from dataclasses import replace
    real = replace(real, drift=spy)
    origin, d = {"F": 2.0, "chi": 0.3}, {"F": 0.6, "chi": 0.8}
    path = construct.along_curved_path(real, origin, d, 0.5)
    t = 0.4
    path.drift(np.zeros((1, 2)), {**path.params, "t": t})
    assert abs(seen["F"] - (2.0 + t * 0.6 + 0.5 * t * t * (-0.8))) < 1e-12
    assert abs(seen["chi"] - (0.3 + t * 0.8 + 0.5 * t * t * 0.6)) < 1e-12
    straight = construct.along_curved_path(real, origin, d, 0.0)
    straight.drift(np.zeros((1, 2)), {**straight.params, "t": t})
    assert abs(seen["F"] - (2.0 + t * 0.6)) < 1e-12


def _pair():
    """Two pitchforks coupled through their mean."""
    base = json.loads((EX / "pitchfork.json").read_text())
    base.pop("potential", None)
    base["carrier"] = {"kind": "euclid", "name": "2 copies", "variables": ["x_1", "x_2"], "scale": 2}
    base["parameters"] = {"eps": 1.0, "g_pop": 0.05}
    base["drift"] = {"x_1": "eps*x_1 - x_1**3 + g_pop*((x_1 + x_2)/2 - x_1)",
                     "x_2": "eps*x_2 - x_2**3 + g_pop*((x_1 + x_2)/2 - x_2)"}
    return spec.load(base)


def test_the_saturation_time_grows_at_weak_noise_and_is_absent_at_strong_noise():
    assert completion.saturation_time(0.05) == pytest.approx(1.46, abs=0.02)
    assert completion.saturation_time(1e-9) == pytest.approx(4.58, abs=0.02)
    assert completion.saturation_time(0.2) is None


def test_final_states_still_at_the_unwritten_state_are_unsettled():
    real = spec.load(EX / "pitchfork.json")
    stable = [np.array([1.0]), np.array([-1.0])]
    reference = np.array([0.0])
    assert completion.unsettled_fraction(real, np.array([[0.1], [-0.2], [0.9], [-1.1]]), stable, reference) == 0.5
    assert completion.unsettled_fraction(real, np.array([[0.9], [-1.1], [np.nan]]), stable, reference) == 0.0


def test_a_sweep_too_short_for_very_weak_noise_is_extended_and_gets_the_share_right():
    """Two directions at Gamma = 1e-9: the fluctuation saturates at tau = 4.6, after the window of tau = 4. Over the old
    window every trajectory is still at the unwritten state, and the mixed states, which take about 0.12 of the
    trajectories, would take about a half of them by the labels."""
    real = _pair()
    nf = an.normal_form(real, np.zeros(2), real.control, 0.0)
    a3 = abs(nf["a3"])
    _, L = completion.critical_space(real, nf, 2)
    D = 1e-9 * 0.02 / (a3 * float(np.mean(np.diag(L.T @ L))))
    sweep = completion.Sweep(real, nf, 0.02)
    rng = np.random.default_rng(1)
    q0 = completion._relaxed(real, np.zeros(2), sweep.over(0.0))
    R, _ = completion.critical_space(real, nf, 2)
    stable = completion.stable_states(real, rng, sweep.over(sweep.T), seeds=completion.swept_seeds(real, sweep, q0, R, sweep.T))
    old = completion.certificate(real, sweep, q0, [q0 + np.array([1.0, 0.0]) * 1e3], [0.0], D, 200, rng)[0]
    assert completion.unsettled_fraction(real, old, stable, np.zeros(2)) > 0.9

    r = completion.law_at(real, nf, np.random.default_rng(5), rate=0.02, noise=D, n_samples=200, certify=False)
    assert r["sweep_extended_for_saturation"]["saturation_time"] == pytest.approx(4.58, abs=0.02)
    assert "unsettled_at_end" not in r
    st = np.array(r["stable_states"])
    mixed = [i for i, s in enumerate(st) if len({x > 0 for x in s}) > 1]
    assert np.mean([sum(row["law"][i] for i in mixed) for row in r["rows"]]) < 0.25


def test_the_standard_noise_rule_stays_inside_the_window_of_the_sweep():
    # Gamma = 0.05 saturates at tau = 1.46 and the window ends at tau = 4, so reference records keep their sweep
    assert completion.saturation_time(0.05) + completion.SETTLE_MARGIN < 4.0


def test_the_crossover_reads_the_flow_until_the_spread_reaches_the_written_amplitude():
    """At Gamma = 0.002 the pulled-back spread of the pair at x = 0.35 is about two thirds of the written amplitude at
    the mapping time. The sector law there puts the mixed states at about 0.37, the full dynamics at about 0.26; the
    flow from the mapping time follows the full dynamics. At the reference noise the spread exceeds the amplitude and
    the sector law is kept."""
    real = _pair()
    nf = an.normal_form(real, np.zeros(2), real.control, 0.0)
    a3 = abs(nf["a3"])
    _, L = completion.critical_space(real, nf, 2)
    dfac = float(np.mean(np.diag(L.T @ L)))
    weak = completion.law_at(real, nf, np.random.default_rng(3), rate=0.02, noise=0.002 * 0.02 / (a3 * dfac),
                             n_samples=1000, certify=False)
    assert weak["map_regime"] == "linear at the mapping time"
    st = np.array(weak["stable_states"])
    mixed = [i for i, s in enumerate(st) if len({x > 0 for x in s}) > 1]
    share = np.mean([sum(row["law"][i] for i in mixed) for row in weak["rows"]])
    sector = np.mean([sum(row["sector_law"][i] for i in mixed) for row in weak["rows"]])
    assert share < sector - 0.05
    ref = completion.law_at(real, nf, np.random.default_rng(3), rate=0.02, n_samples=400, certify=False)
    assert ref["map_regime"] == "saturated at the mapping time"
