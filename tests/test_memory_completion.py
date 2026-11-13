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


def test_the_critical_direction_does_not_depend_on_the_sign_the_eigensolver_gives(monkeypatch):
    """The sign of an eigenvector is arbitrary, and two LAPACK builds returned opposite signs at the cusp of the driven
    trapped ion. The swept-write check gives the bias toward +mode its own random streams, so the two machines
    measured independent samples (law constants 1.47 and 0.90) and offsets of opposite sign. The critical direction
    is oriented so that its largest component (the first of near-equal ones) is positive."""
    real = spec.load(Path(__file__).resolve().parents[1] / "examples" / "quantum" / "semiclassical" / "ion_duffing.json")
    q, v = np.array([2.752267826816943, -1.806257226422566]), 1.3037475667148828   # its first fold in the record
    nf = an.normal_form(real, q, "Delta", v)
    eig = np.linalg.eig
    monkeypatch.setattr(np.linalg, "eig", lambda M: (lambda e, V: (e, -V))(*eig(M)))
    flipped = an.normal_form(real, q, "Delta", v)
    assert nf["mode"][an.orientation_index(nf["mode"])] > 0
    for key in ("mode", "left_mode", "a2", "bias_along_mode"):
        assert np.allclose(flipped[key], nf[key], rtol=0, atol=1e-12), key


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
    # the standard noise rule keeps the pitchfork's write: the hold, which is always read, retains the choice, and the class is
    # then the linear-stage write
    rep = completion.complete(real, np.random.default_rng(3), floors=(1.0,), n_samples=400, n_traj=300)
    L = next(e["law"] for e in rep["write_points"] if "law" in e and e["law"].get("status") == "derived")
    assert L["write_class"] == "linear-stage write"
    assert L["hold"]["switching_fraction"] <= completion.SWITCHING_EQUILIBRIUM


def test_the_sector_law_in_one_direction_is_the_sign_of_the_sample():
    z = np.array([[0.3], [-2.0], [5.0]])
    assert completion._sector_labels(z, np.array([[1.0], [-1.0]])).tolist() == [0, 1, 0]
    # a state near the unwritten one is no sector
    assert completion._sector_labels(z, np.array([[1.0], [0.1], [-1.0]])).tolist() == [0, 2, 0]


def test_written_states_in_one_direction_have_no_sector_law():
    """The direction of a sample does not tell apart two written states on the same side of the unwritten one."""
    assert completion._sector_labels(np.array([[1.0]]), np.array([[1.0], [1.6]])) is None
    plane = np.array([[10.75, -1.07], [16.38, -1.62], [-10.0, 0.0]])
    assert completion._sector_labels(np.array([[1.0, 0.0]]), plane) is None
    assert completion._sector_labels(np.array([[1.0, 0.0]]), np.array([[1.0], [-1.0]]) * np.array([[1.0, 0.0]])) is not None


def test_a_saturated_deviation_in_one_direction_keeps_the_law_of_the_sign():
    """The phase-locked oscillator: at the mapping time the pulled-back spread is several times the written amplitude,
    and a flow from samples of that size crosses to the other state along the Kerr rotation (response 0.81 at a bias of
    one floor, where the full dynamics gives 1.33). The sign of the sample decides, and the law is the pitchfork's."""
    real = spec.load(EX.parent / "quantum" / "semiclassical" / "pplo_oscillator.json")
    nf = an.normal_form(real, np.zeros(2), real.control, 1.0)
    L = completion.law_at(real, nf, np.random.default_rng(0), n_samples=4000, certify=False)
    assert L["map_regime"] == "saturated at the mapping time"
    assert abs(L["one_direction"]["law"][0]["response_constant"] - pi ** 0.25) < 0.1


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


QUANTUM = EX.parents[0] / "quantum" / "semiclassical"


def test_the_one_direction_law_keeps_its_constant_at_every_bias():
    """The phase-locked parametric oscillator: the full dynamics keep the constant of the pitchfork at 0.5, 1 and 2 noise
    floors. At the mapping time the linear stage gives the samples up to four written amplitudes, where no trajectory is,
    and a noise-free flow from that size follows the amplitude-dependent rotation of the Kerr term to either state: the
    law read from it fell with the bias (1.03, 0.79, 0.21). Only the sign of a sample counts in one direction."""
    real = spec.load(QUANTUM / "pplo_oscillator.json")
    rng = np.random.default_rng(0)
    events, _ = an.distinct_write_points(real, an.locate_writes(real, rng, n_starts=24))
    nf = an.normal_form(real, events[0]["q"], real.control, events[0]["v"])
    L = completion.law_at(real, nf, rng, floors=(0.5, 1.0, 2.0), n_samples=3000, certify=False)
    assert L["dynamical_critical_dimension"] == 1
    for row in L["one_direction"]["law"]:
        assert abs(row["response_constant"] - pi ** 0.25) < 0.15, row


def test_a_flow_in_one_direction_starts_at_a_size_that_keeps_the_sign():
    """At a Kerr rate a quarter of that of the phase-locked oscillator no written amplitude is defined at the mapping
    time, the deviation is not counted as saturated and the flow is the map. Started at the size the linear stage gives
    the samples it read 1.18 at two noise floors, where the full dynamics give 1.34; started small it keeps the sign."""
    source = json.loads((QUANTUM / "pplo_oscillator.json").read_text(encoding="utf-8"))
    source["parameters"]["gnl"] = 0.0005
    real = spec.load(source)
    rng = np.random.default_rng(1)
    events, _ = an.distinct_write_points(real, an.locate_writes(real, rng, n_starts=24))
    nf = an.normal_form(real, events[0]["q"], real.control, events[0]["v"])
    L = completion.law_at(real, nf, rng, floors=(2.0,), n_samples=8000, certify=False)
    assert L["rows"][0]["regime"] == "linear at the mapping time"
    assert abs(L["one_direction"]["law"][0]["response_constant"] - pi ** 0.25) < 0.08


def test_the_kerr_cat_qubit_is_not_certified_by_a_response_near_zero():
    """At the Kerr strength of the Kerr-cat qubit a flow from the size the linear stage gives the samples scrambles them,
    and the law fell to the 0.05 of the full dynamics: every probability within 0.02 of one half, certified, while half of
    the trajectories change state during a hold. The law of the linear stage keeps its constant, and the rung names the
    difference an equilibrium write."""
    real = spec.load(QUANTUM / "kerr_cat_kpo.json")
    rng = np.random.default_rng(0)
    events, _ = an.distinct_write_points(real, an.locate_writes(real, rng, n_starts=24))
    nf = an.normal_form(real, events[0]["q"], real.control, events[0]["v"])
    L = completion.law_at(real, nf, rng, n_samples=1000, n_traj=500)
    assert abs(L["one_direction"]["law"][0]["response_constant"] - pi ** 0.25) < 0.25
    assert L["certified"] is False and L["write_class"] == "equilibrium write"


def test_past_a_fold_there_is_no_choice_to_saturate():
    """A fold with a slow second mode has two directions that take part, but at the mapping time the unwritten branch is
    gone, and the written states do not lie in distinct directions from it. The sector law gave them probabilities they
    do not have (deviation 1.0), while the deterministic motion after the fold ends in one state."""
    source = json.loads((QUANTUM / "phase_qubit_rfsquid.json").read_text(encoding="utf-8"))
    source["parameters"]["gam"] = 0.1
    real = spec.load(source)
    rng = np.random.default_rng(0)
    events, _ = an.distinct_write_points(real, an.locate_writes(real, rng, n_starts=24))
    nf = an.normal_form(real, events[0]["q"], real.control, events[0]["v"])
    assert "fold" in nf["kind"] and an.dynamical_critical_dimension(nf, 0.02, completion.SLOW_FACTOR) == 2
    L = completion.law_at(real, nf, rng, n_samples=500, n_traj=200, floors=(1.0,))
    if L["status"].startswith("the sweep does not reach the write point"):
        # the sweep starts in another well and never comes to the fold: nothing is chosen, and there is no law to saturate
        return
    assert "saturated at the mapping time" not in {r["regime"] for r in L["rows"]}
    assert max(r["deviation"] for r in L["rows"]) < 0.1


def test_the_swept_seeds_reach_every_frustrated_state_of_a_triangle():
    """The axes of the critical plane are 30 degrees from two of the six frustrated states of the rf-SQUID triangle and
    lead to the four others; the two left were found only when a random start fell into their basin. Where it did not,
    their trajectories were counted to the state of equal phases, and the law differed from the full dynamics by 0.355."""
    real = spec.load(QUANTUM / "rfsquid_triangle.json")
    rng = np.random.default_rng(0)
    events, _ = an.distinct_write_points(real, an.locate_writes(real, rng, n_starts=12))
    nf = an.normal_form(real, events[0]["q"], real.control, events[0]["v"])
    R, _ = completion.critical_space(real, nf, 2)
    sweep = completion.Sweep(real, nf, 0.02)
    q0 = completion._relaxed(real, np.asarray(nf["state"], float), sweep.over(0.0))
    seeds = completion.swept_seeds(real, sweep, q0, R, sweep.T)
    found = an.zeros_from(real, seeds, sweep.over(sweep.T))
    stable = [q for q, kap in found if an.n_unstable(kap) == 0]
    frustrated = {tuple(np.sign(q[::2]).astype(int)) for q in stable if len(set(np.sign(q[::2]))) == 2}
    assert len(frustrated) == 6, frustrated


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


def test_written_states_closer_than_the_merge_tolerance_are_found_at_the_mapping_time():
    # a pitchfork shortly after the crossing on a carrier of large scale: the written states at +-0.1 lie closer to the
    # unwritten branch than the tolerance within which stored states are one, and the branch is found first
    real = spec.load({"schema": "fieldbridge-memory/1", "kind": "equations", "name": "test", "question": "test",
                      "assumptions": [], "provenance": {"source": "test"},
                      "carrier": {"kind": "euclid", "variables": ["x"], "scale": 4.0}, "parameters": {"a": 0.01},
                      "drift": {"x": "a*x - x**3"}, "noise": 0.001, "control": {"name": "a", "range": [-1.0, 1.0]}})
    seeds = np.array([[0.001], [0.12], [-0.12]])
    rng = np.random.default_rng(0)
    assert completion.stable_states(real, rng, {}, n_starts=0, seeds=seeds) == []
    found = sorted(float(q[0]) for q in completion.stable_zeros(real, rng, {}, n_starts=0, seeds=seeds))
    assert found == pytest.approx([-0.1, 0.1], abs=1e-6)


def test_written_states_closer_than_the_merge_tolerance_are_kept_at_the_end_of_the_sweep():
    # a pitchfork whose written states at +-0.05 lie closer together than the tolerance of stored states (0.2 on a carrier
    # of scale 4): stable_states keeps one of them, end_states keeps both and merges only zeros that are the same root
    real = spec.load({"schema": "fieldbridge-memory/1", "kind": "equations", "name": "test", "question": "test",
                      "assumptions": [], "provenance": {"source": "test"},
                      "carrier": {"kind": "euclid", "variables": ["x"], "scale": 4.0}, "parameters": {"a": 0.0025},
                      "drift": {"x": "a*x - x**3"}, "noise": 0.001, "control": {"name": "a", "range": [-1.0, 1.0]}})
    seeds = np.array([[0.06], [0.07], [-0.06]])
    rng = np.random.default_rng(0)
    assert len(completion.stable_states(real, rng, {}, n_starts=0, seeds=seeds)) == 1
    found = sorted(float(q[0]) for q in completion.end_states(real, rng, {}, n_starts=0, seeds=seeds))
    assert found == pytest.approx([-0.05, 0.05], abs=1e-6)


def test_the_sign_in_one_direction_is_read_from_the_unwritten_branch_at_the_mapping_time():
    """The magnon-Kerr element: the sweep through the unfolded cusp moves the unwritten branch by more than the distance
    to the written states, and from the state at the write point both states lie on one side (critical coordinates
    0.35 and 2.85; -1.39 and 1.10 from the branch at the mapping time). Read from the write point there was no sector
    law, and the flow from small samples sent every sample to one state (a response of 0, where the full dynamics give
    1.26); read from the branch, the sign of the sample gives the response of the pitchfork."""
    real = spec.load(QUANTUM / "magnon_kerr.json")
    rep = completion.complete(real, np.random.default_rng(0), n_samples=1000, certify=False)
    # an unfolding along a parameter that finds no cusp has no law (on Linux one of the two does not)
    laws = [u["law"] for e in rep["write_points"] for u in e.get("unfoldings", []) if (u.get("law") or {}).get("status") == "derived"]
    assert laws
    for L in laws:
        assert L["map_regime"] == "saturated at the mapping time"
        assert abs(L["one_direction"]["law"][0]["response_constant"] - pi ** 0.25) < 0.2


def _slow_fast(c=0.7, k=2.0, m=1.0):
    return spec.load({"schema": "fieldbridge-memory/1", "name": "slow and fast", "kind": "equations", "question": "test",
                      "provenance": {"source": "test", "kind": "authored example"}, "closure": "test", "observable": "x", "assumptions": ["test"],
                      "field": "test", "carrier": {"kind": "euclid", "name": "plane", "variables": ["x", "y"], "scale": 2},
                      "parameters": {"eps": 0.0, "c": c, "k": k, "m": m},
                      "drift": {"x": "eps*x - x**3 + c*y**2", "y": "-k*y + m*x"},
                      "control": {"name": "eps", "range": [-0.5, 0.5]}, "noise": 0.05, "dt": 0.01})


def test_the_noise_of_the_fast_direction_gives_the_slow_coordinate_the_drift_of_its_curvature():
    """x' = eps x - x^3 + c y^2, y' = -k y + m x at eps = 0: the slow direction is v = (1, m/k), its left vector
    w = (1, 0). With y = (m/k) x + z, z' = -k z + xi_y - (m/k) xi_x: <z^2> = D (1 + (m/k)^2)/k and <x z> = -2 D m/k^2,
    and c y^2 gives x the drift c <z^2> + 2 c (m/k) <x z> = c D (1 - 3 (m/k)^2)/k beside its own (m/k)^2 c x^2.
    Leaving it out made the offset of the write at a Kerr cusp 1.00 where the full dynamics give 0.78; the fast part
    alone, c D (1 + (m/k)^2)/k, broke the law of the phonon laser, where the slow and the fast rate come close."""
    c, k, m, D = 0.7, 2.0, 1.0, 0.05
    real = _slow_fast(c, k, m)
    q = np.zeros(2)
    v = np.array([1.0, m / k]); v /= np.linalg.norm(v)
    w = np.array([1.0, 0.0]); w /= w @ v
    b = completion.fast_noise_drift(real, q, {"eps": 0.0}, v, w, D)
    assert abs(float(b[0]) - c * D * (1 - 3 * (m / k) ** 2) / k) < 1e-6 and abs(float(b[1])) < 1e-9
    assert np.allclose(completion.fast_noise_drift(real, q, {"eps": 0.0}, v, w, 0.0), 0.0)
    one = spec.load(EX / "pitchfork.json")
    assert np.allclose(completion.fast_noise_drift(one, np.zeros(1), {"eps": 0.0}, np.ones(1), np.ones(1), D), 0.0)


def test_the_noise_across_the_critical_space_pushes_the_mean_through_the_curvature():
    # x' = e x - x^3 + c y^2, y' = -k y: the fluctuations of y give x the mean force c D / k; a weakly damped mass,
    # whose relaxing mode lies almost along the critical one, gets none from its velocity, which the force does not bend
    def load(drift, variables, parameters):
        return spec.load({"schema": "fieldbridge-memory/1", "kind": "equations", "name": "test", "question": "test",
                          "assumptions": [], "provenance": {"source": "test"},
                          "carrier": {"kind": "euclid", "variables": variables, "scale": 3.0},
                          "parameters": parameters, "drift": drift, "noise": 0.001,
                          "control": {"name": "e", "range": [-1.0, 1.0]}})
    D = 0.01
    pumped = load({"x": "e*x - x**3 + c*y**2", "y": "-k*y"}, ["x", "y"], {"e": 0.0, "c": 3.0, "k": 2.0})
    g = completion.noise_curvature_drift(pumped, np.zeros(2), {}, np.array([[1.0], [0.0]]), D)
    assert g == pytest.approx([3.0 * D / 2.0, 0.0], rel=1e-3, abs=1e-9)
    mass = load({"x": "v", "v": "-g*v + e*x - x**3 + 2*x**2"}, ["x", "v"], {"e": 0.0, "g": 0.3})
    g = completion.noise_curvature_drift(mass, np.zeros(2), {}, np.array([[1.0], [0.0]]), D)
    assert g == pytest.approx([0.0, 0.0], abs=1e-9)


def test_the_second_order_shift_is_the_mean_of_the_curvature_among_the_undecided_samples():
    # Two cases with a closed form. The critical coordinate diffuses (rate zero) for a time T, so sigma^2 = 2 D T.
    # (1) x' = a x^2: the curvature acts on the critical coordinate itself. Among all samples <x^2> = 2 D t and the
    # push is a D T^2; among the samples that end at x = 0 the variance is that of a bridge, 2 D t (1 - t / T), and the
    # push is a third of it. (2) x' = c y^2, y' = -k y: y is not correlated with x, the condition leaves it alone, and
    # the push is c times the integral of <y^2> = (D / k)(1 - exp(-2 k t)), with or without the condition.
    class Still:
        def over(self, t):
            return {}

    def load(drift, variables, parameters):
        return spec.load({"schema": "fieldbridge-memory/1", "kind": "equations", "name": "test", "question": "test",
                          "assumptions": [], "provenance": {"source": "test"},
                          "carrier": {"kind": "euclid", "variables": variables, "scale": 3.0},
                          "parameters": parameters, "drift": drift, "noise": 0.001,
                          "control": {"name": "e", "range": [-1.0, 1.0]}})
    D, T, N = 0.01, 5.0, 500
    ts = np.linspace(0.0, T, N + 1)
    sigma = np.sqrt(2 * D * T)
    a = 0.4
    one = load({"x": "a*x**2 + 0*e"}, ["x"], {"e": 0.0, "a": a})
    so = completion.second_order_selection(one, Still(), ts, np.zeros((N + 1, 1)), np.zeros((N + 1, 1, 1)),
                                           (2 * D * ts)[:, None, None], np.zeros((N + 1, 1, 1)), np.ones(1), 0.0, D, n_mc=4000)
    assert so["Delta_unconditioned"] == pytest.approx(a * D * T ** 2 / sigma, rel=1e-3)
    assert so["Delta_cov"] == pytest.approx(a * D * T ** 2 / (3 * sigma), rel=1e-3)
    assert so["Delta"][0] == pytest.approx(so["Delta_cov"], abs=1e-12)          # no bias, no forcing
    assert so["sampled"]["rows"][0]["conditional_mean_over_sigma"] == pytest.approx(so["Delta_cov"], rel=0.1)
    c, k = 3.0, 2.0
    two = load({"x": "c*y**2 + 0*e", "y": "-k*y"}, ["x", "y"], {"e": 0.0, "c": c, "k": k})
    J = np.array([[0.0, 0.0], [0.0, -k]])
    C = np.zeros((N + 1, 2, 2))
    C[:, 0, 0] = 2 * D * ts
    C[:, 1, 1] = (D / k) * (1 - np.exp(-2 * k * ts))
    so = completion.second_order_selection(two, Still(), ts, np.zeros((N + 1, 2)), np.repeat(J[None], N + 1, axis=0), C,
                                           np.zeros((N + 1, 1, 2)), np.array([1.0, 0.0]), 0.0, D, n_mc=4000)
    exact = c * (D / k) * (T - (1 - np.exp(-2 * k * T)) / (2 * k)) / sigma
    assert so["Delta_cov"] == pytest.approx(exact, rel=1e-3)
    assert so["Delta_unconditioned"] == pytest.approx(exact, rel=1e-3)
    assert 0.0 < so["eta"] < completion.SECOND_ORDER_SPREAD


def test_a_write_in_one_direction_carries_the_second_order_shift_and_its_spread():
    # a pitchfork whose slow partner pumps it: the law keeps the quantities of the shift, and a spread at or above
    # the limit withdraws the claim
    law = completion.complete(spec.load(EX / "pitchfork.json"), np.random.default_rng(0), n_samples=500,
                              certify=False)["write_points"][0]["law"]
    so = law["second_order"]
    assert so["eta"] < completion.SECOND_ORDER_SPREAD and "second_order_not_small" not in law
    assert len(so["Delta"]) == len(law["rows"])


def test_a_symmetric_critical_direction_is_oriented_by_its_first_component_not_by_rounding():
    from fieldbridge.memory import analysis as an
    assert an.orientation_index(np.array([0.7071067811865475, -0.7071067811865476])) == 0
    assert an.orientation_index(np.array([-0.7071067811865476, 0.7071067811865475])) == 0
    assert an.orientation_index(np.array([0.3, -0.9, 0.2])) == 1
