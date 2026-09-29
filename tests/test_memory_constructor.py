"""The memory constructor: specifications, structural predictions, calculations and their agreement.

Expected values are analytic where one exists (write points, cusps, kernels, loop frustration), and the controls
are cases whose predictions must change: a broken symmetry, a flipped loop sign, a mixed bond.
"""
import json
import math
import subprocess
import sys
from pathlib import Path

import pytest

np = pytest.importorskip("numpy")
pytest.importorskip("scipy")
pytest.importorskip("sympy")

from fieldbridge.memory import analysis as an
from fieldbridge.memory import compose, construct, networks, predict, spec

ROOT = Path(__file__).resolve().parents[1]
EX = ROOT / "examples" / "memory"


def load(name):
    return spec.load(EX / f"{name}.json")


def events(real, seed=0):
    return an.locate_writes(real, np.random.default_rng(seed), n_starts=16)


# ------------------------------------------------------------------------------------------------ specifications
def test_every_example_loads_with_its_provenance():
    for path in sorted(EX.glob("*.json")):
        real = spec.load(path)
        assert real.spec_sha256 and real.spec["question"] and real.spec["assumptions"]


@pytest.mark.parametrize("bad, message", [
    ({"drift": {"u": "__import__('os')", "v": "v"}}, "Unsupported syntax"),
    ({"drift": {"u": "alpha/(1 + w**n) - u", "v": "alpha/(1 + u**n) - v"}}, "Undeclared symbol"),
    ({"drift": {"u": "u.real", "v": "v"}}, "Unsupported syntax"),
    ({"control": {"name": "beta", "range": [0, 1]}}, "must be a declared parameter"),
    ({"schema": "other"}, "schema must be"),
])
def test_specifications_outside_the_contract_are_refused(bad, message):
    base = json.loads((EX / "toggle.json").read_text(encoding="utf-8"))
    base.update(bad)
    with pytest.raises(spec.SpecError, match=message):
        spec.load(base)


def test_a_potential_must_generate_the_drift():
    base = json.loads((EX / "schlogl.json").read_text(encoding="utf-8"))
    base["potential"] = "x**4/4"
    with pytest.raises(spec.SpecError, match="drift = -grad"):
        spec.load(base)


# ------------------------------------------------------------------------------------------------ predictions
SYMMETRIC_WRITE = "a symmetric state that loses stability"


def structure(name):
    return predict.predict(load(name), np.random.default_rng(0))


def test_symmetry_restricts_the_kind_of_write_before_any_calculation():
    assert structure("toggle")["predictions"]["write"]["prediction"].startswith(SYMMETRIC_WRITE)
    assert structure("toggle_unequal")["predictions"]["write"]["prediction"].startswith("one-sided")
    assert structure("repressor_ring4")["predictions"]["write"]["prediction"].startswith(SYMMETRIC_WRITE)
    assert structure("repressilator")["predictions"]["write"]["prediction"].startswith("no symmetric pitchfork")
    assert structure("schlogl")["predictions"]["write"]["prediction"].startswith("one-sided")


def test_loop_signs_decide_storage_and_oscillation():
    ring3, ring4 = structure("repressilator")["predictions"], structure("repressor_ring4")["predictions"]
    assert ring3["multistability"]["prediction"].startswith("impossible")
    assert ring3["oscillation"]["prediction"].startswith("possible")
    assert ring4["multistability"]["prediction"].startswith("possible")
    assert ring4["oscillation"]["prediction"].startswith("not expected")


def test_control_roles_and_reciprocity():
    colloid = structure("colloid_patch")
    assert colloid["control_role"] == "scale" and colloid["reciprocal"]
    assert colloid["predictions"]["oscillation"]["prediction"] == "impossible"
    assert any(g.get("shift") for g in colloid["symmetries"])  # 90-degree partner states of nematic rods
    assert structure("schlogl")["control_role"] == "bias"
    assert structure("compartments").get("linear")


@pytest.mark.parametrize("name", ["pitchfork", "toggle", "toggle_unequal", "schlogl", "tubes", "tubes_unequal",
                                  "repressilator", "repressor_ring4", "compartments"])
def test_structural_predictions_agree_with_the_calculation(name):
    from fieldbridge.memory import discovery
    real = load(name)
    pred = predict.predict(real, np.random.default_rng(0))
    card = discovery.evaluate(real, np.random.default_rng(1), quick=True)
    comparison = predict.compare(pred, card)
    assert comparison["all_consistent"]
    if name in ("pitchfork", "toggle", "tubes", "repressor_ring4"):
        # the prediction is conditional; these models have a write point to which the condition applies
        assert comparison["write"]["tested_write_points"] >= 1


# ------------------------------------------------------------------------------------------------ calculations
def test_pitchfork_write_point_and_barrier():
    real = load("pitchfork")
    (ev,) = events(real)
    nf = an.normal_form(real, ev["q"], "eps", ev["v"])
    assert ev["v"] == pytest.approx(0.0, abs=1e-6) and nf["a3"] == pytest.approx(-1.0, rel=1e-6)
    states, *_ = an.stored_states(real, np.random.default_rng(1), 24)
    assert an.neb_barrier(real, states[0], states[1], iters=2000)["barrier"] == pytest.approx(0.25, abs=2e-3)


def test_toggle_pitchfork_and_transferred_write_law():
    real = load("toggle")
    (ev,) = events(real)
    nf = an.normal_form(real, ev["q"], "alpha", ev["v"])
    assert nf["value"] == pytest.approx(2.0, abs=1e-6) and nf["kappa_slope"] == pytest.approx(-0.25, abs=1e-4)
    law = construct.swept_write_check(real, nf, np.random.default_rng(7), target_p=0.8, n_traj=600)
    assert abs(law["measured"] - 0.8) < 4 * law["stderr"] + 0.02


def test_schlogl_cusp_removes_the_one_sided_write():
    real = load("schlogl")
    ev = events(real)
    assert [round(e["v"], 3) for e in ev] == [1.490, 2.260]
    nf = an.normal_form(real, ev[0]["q"], "b", ev[0]["v"])
    c = construct.cancel_asymmetry(real, nf, "a")
    a_c = math.sqrt(3 * 5.75)
    assert c["success"] and c["value"] == pytest.approx(a_c, abs=1e-6)
    assert c["direction"]["b"] / c["direction"]["a"] == pytest.approx(-(a_c / 3) ** 2, rel=1e-4)
    assert c["normal_form_along_sweep"]["kind"].startswith("supercritical pitchfork")


def test_unequal_tubes_are_cancelled_only_by_equal_lengths():
    from scipy.optimize import brentq
    mu_c = brentq(lambda m: m - 1 - 2 ** -m, 1.0, 2.0)
    assert any(abs(e["v"] - mu_c) < 1e-6 for e in events(load("tubes")))
    real = load("tubes_unequal")
    nf = an.normal_form(real, events(real)[0]["q"], "mu", events(real)[0]["v"])
    c = construct.cancel_asymmetry(real, nf, "L2")
    assert c["success"] and c["value"] == pytest.approx(1.0, abs=1e-6) and c["control"] == pytest.approx(mu_c, abs=1e-6)


def test_repressor_rings_hopf_and_pitchfork():
    ev3 = events(load("repressilator"))[0]
    nf3 = an.normal_form(load("repressilator"), ev3["q"], "alpha", ev3["v"])
    assert nf3["hopf"] and nf3["value"] == pytest.approx(2.0, abs=1e-6) and nf3["omega"] == pytest.approx(math.sqrt(3), rel=1e-5)
    ev4 = events(load("repressor_ring4"))[0]
    u = 3 ** -0.25
    assert ev4["v"] == pytest.approx(u * (1 + u ** 4), abs=1e-6)


def test_compartment_kernel():
    from fieldbridge.memory.library import compartments
    t = np.linspace(0, 1, 5)
    assert np.allclose(an.memory_kernel(compartments(), t), 1.5 * 4.5 * np.exp(-3.9 * t))


# ------------------------------------------------------------------------------------------------ networks and loops
def test_rotor_constraints_exact_cases():
    tri = compose._ring_positions(3)
    sq = compose._ring_positions(4)
    phis = lambda pos, n: [math.atan2(*(pos[(k + 1) % n] - pos[k])[::-1]) for k in range(n)]
    ring = lambda n, kind: [(k, (k + 1) % n, kind) for k in range(n)]
    t = networks.rotor_constraints(3, ring(3, "reflect"), phis(tri, 3), 2)
    assert t["satisfiable"] and not t["continuous_symmetry"]
    s = networks.rotor_constraints(4, ring(4, "reflect"), phis(sq, 4), 2)
    assert s["satisfiable"] and s["continuous_symmetry"]  # a square: zero flux, free signed rotation
    anti3 = networks.rotor_constraints(3, ring(3, "anti"), phis(tri, 3), 1)
    assert not anti3["satisfiable"]  # the antiferromagnetic triangle


def test_holonomy_predicts_loops_in_three_fields():
    rng = np.random.default_rng(5)
    assert compose.evaluate_loop(compose.gene_loop(3, [True] * 3), rng)["measured"] == "oscillation"
    assert compose.evaluate_loop(compose.gene_loop(2, [True, True]), rng)["stable_states"] == 2
    s3 = compose.evaluate_loop(compose.spin_loop(3, [-1, 1, 1]), rng)
    assert s3["unsatisfied_bonds"] == 1
    for _ in range(3):
        r = compose.evaluate_loop(compose.rotor_loop(4, ["reflect"] * 4, compose._ring_positions(4, rng, 0.3)), rng)
        assert r["measured"] == pytest.approx(1 - math.cos(r["holonomy"]["flux"] / 4), abs=1e-6)


def test_random_networks_follow_their_structure():
    rows = networks.benchmark(np.random.default_rng(3), {"gene": 4, "spin": 4, "rotor": 4}, sizes=(4, 5))
    assert all(r["consistent"] for r in rows)


# ------------------------------------------------------------------------------------------------ command line
def test_command_line_writes_reports_with_provenance(tmp_path):
    out = tmp_path / "toggle"
    cmd = [sys.executable, "-B", "-m", "fieldbridge", "memory", "predict", str(EX / "toggle.json"), "--out-dir", str(out)]
    subprocess.run(cmd, check=True, cwd=ROOT, capture_output=True)
    report = json.loads((out / "predict.json").read_text(encoding="utf-8"))
    assert report["novelty_established"] is False and report["input_sha256"] and report["implementation_sha256"]
    assert report["structure"]["predictions"]["write"]["prediction"].startswith(SYMMETRIC_WRITE)


# ------------------------------------------------------------------------------------------------ phase memory
def test_a_parameter_may_be_named_like_the_state_argument():
    """A copying fidelity or a charge is often called q, the name of the state argument of Realization.F."""
    s = {"schema": "fieldbridge-memory/1", "kind": "equations", "name": "master sequence",
         "question": "At which copying fidelity is the master sequence lost?", "assumptions": ["single peak"],
         "provenance": {"source": "Eigen, Naturwissenschaften 58, 465 (1971)"},
         "carrier": {"kind": "euclid", "variables": ["x"], "scale": 1.0},
         "parameters": {"f0": 2.0, "q": 0.8}, "drift": {"x": "x*(f0*q - 1 - (f0 - 1)*x)"},
         "control": {"name": "q", "range": [0.2, 1.0]}, "noise": 0.002}
    real = spec.load(s)
    assert real.F(np.array([0.5]), q=0.5)[0] == pytest.approx(-0.25)
    pred = predict.predict(real, np.random.default_rng(0))  # failed with "multiple values for argument 'q'"
    assert pred["control_role"] in ("shape", "bias")


def test_a_fold_at_a_scan_value_is_located():
    """The overdamped Josephson junction loses its stable phase at the critical current 1, which is one of the
    values of the scan along the control; the fold must be found there, not missed."""
    from fieldbridge.memory import construct
    real = spec.load(EX / "oscillators" / "josephson.json")
    events = construct.construct(real, np.random.default_rng(1), check_write=False)["events"]
    folds = [e for e in events if e["source"] == "control" and e["kind"].startswith("saddle-node")]
    assert folds and abs(folds[0]["value"] - 1.0) < 1e-3


def test_a_card_separates_a_limit_cycle_from_a_family_of_neutral_cycles():
    """Lotka-Volterra orbits are neutral (a conserved quantity); the van der Pol cycle is isolated."""
    from fieldbridge.memory import discovery
    osc = EX / "oscillators"
    lv = discovery.evaluate(spec.load(osc / "lotka_volterra.json"), np.random.default_rng(1), quick=True)
    vdp = discovery.evaluate(spec.load(osc / "van_der_pol.json"), np.random.default_rng(1), quick=True)
    assert lv["states"]["neutral_cycles"] and lv["states"]["transverse_floquet_multiplier"] > 0.99
    assert discovery.mechanism_class(lv) == "neutral cycles: no isolated phase"
    assert "limit cycle" in vdp["loss_law"] and vdp["states"]["transverse_floquet_multiplier"] < 0.1
    assert discovery.mechanism_class(vdp) == "limit cycle: phase memory"


def test_an_oscillator_remembers_a_pulse_in_its_phase():
    """No stable state, yet memory: the phase is a flat direction (time-translation symmetry)."""
    from fieldbridge.memory import phase
    real = load("repressilator")
    rng = np.random.default_rng(1)
    cyc = phase.limit_cycle(real, rng)
    assert cyc["period"] == pytest.approx(5.93, abs=0.05)
    h = phase.hold(real, cyc, rng, shift=np.pi / 2, D=0.005, n_traj=120, n_periods=20)
    sep, var, t = np.array(h["separation"]), np.array(h["variance"]), np.array(h["time"])
    assert abs(sep[-1] - sep[4]) < 0.1 and abs(sep[-1] - np.pi / 2) < 0.15  # the written shift persists
    assert np.std(var[4:] / t[4:]) / np.mean(var[4:] / t[4:]) < 0.2  # variance linear in time: Law 2
    L = phase.lock(real, cyc, rng, strengths=[0.01, 0.03, 0.1], target=np.pi / 2, D=0.005, n_traj=80)
    assert 0.85 < L["log_log_slope"] < 1.1  # flat-direction lock: ratio linear in the write strength


def test_equal_promoters_remove_the_toggles_one_sided_write():
    real = load("toggle_unequal")
    ev = events(real)[0]
    nf = an.normal_form(real, ev["q"], "alpha", ev["v"])
    assert nf["kind"].startswith("saddle-node")
    c = construct.cancel_asymmetry(real, nf, "gamma")
    assert c["success"] and c["value"] == pytest.approx(1.0, abs=1e-6) and c["control"] == pytest.approx(2.0, abs=1e-6)
    assert c["normal_form_along_sweep"]["kind"].startswith("supercritical pitchfork")


def test_attaching_the_colloid_memory_to_dipoles_keeps_its_class(tmp_path):
    out = tmp_path / "attach"
    cmd = [sys.executable, "-B", "-m", "fieldbridge", "memory", "attach", "--from", str(EX / "colloid_patch.json"),
           "--to", str(EX / "dipole_patch.json"), "--out-dir", str(out)]
    subprocess.run(cmd, check=True, cwd=ROOT, capture_output=True)
    report = json.loads((out / "attach.json").read_text(encoding="utf-8"))
    assert not report["changed"]  # the memory signature survives the step
    assert "carrier" in report["carrier_and_material"]  # period pi -> 2 pi


# ------------------------------------------------------------------------------------------------ time and fields
def test_one_formula_covers_the_three_regimes_of_kappa():
    """I(t) = (1/2) ln[1 + (kappa s0^2/D)/(exp(2 kappa t) - 1)]: exponential loss, 1/t loss, and a plateau."""
    from fieldbridge.memory import regimes
    s0, D, k = 0.3, 0.02, 1.0
    t = np.array([20.0, 40.0])
    W = k * s0 ** 2 / D
    assert regimes.gaussian_information(-k, s0, D, t) == pytest.approx(0.5 * np.log1p(W), rel=1e-12)
    flat = regimes.gaussian_information(0.0, s0, D, np.array([1e3, 2e3]))
    assert flat[0] / flat[1] == pytest.approx(2.0, rel=1e-3)  # Law 2: information ~ 1/t
    fast = regimes.gaussian_information(k, s0, D, t)
    assert np.log(fast[0] / fast[1]) / 20.0 == pytest.approx(2 * k, rel=1e-6)  # Law 1: rate 2 kappa
    small = regimes.gaussian_information(1e-6, s0, D, np.array([5.0]))  # the flat case is the limit of both
    assert small[0] == pytest.approx(regimes.gaussian_information(0.0, s0, D, np.array([5.0]))[0], rel=1e-4)


def test_an_expanding_direction_freezes_the_write():
    from fieldbridge.memory import regimes
    times = [0.5, 3.0, 8.0]
    rng = np.random.default_rng(4)
    up = regimes.simulate(-1.0, 0.3, 0.02, times, rng, n=4000)
    down = regimes.simulate(1.0, 0.3, 0.02, times, rng, n=4000)
    assert abs(up["accuracy"][-1] - up["accuracy"][1]) < 0.02 and up["accuracy"][-1] > 0.95  # kept after the ridge
    assert down["accuracy"][-1] < 0.6  # lost in the contracting case


@pytest.mark.parametrize("d, write, exponent", [(1, "charge", -0.5), (1, "dipole", -1.5), (2, "charge", -1.0),
                                                (2, "dipole", -2.0)])
def test_conservation_and_dimension_set_the_power_law(d, write, exponent):
    from fieldbridge.memory import fields
    model = fields.FieldModel(d=d, write=write, amplitude=3.0, T=0.5)
    assert fields.predicted_law(model)["exponent"] == exponent
    assert fields.exponent_check(model)["exponent"] == pytest.approx(exponent, abs=0.01)


def test_the_finite_total_keeps_a_charge_write():
    """On a closed ring the profile is lost but the total is not: the uniform mode never decays."""
    from fieldbridge.memory import fields
    model = fields.FieldModel(d=1, L=32, write="charge", amplitude=3.0, T=0.5)
    profile = fields.spectral_snr(model, [1e4])
    everything = fields.spectral_snr(model, [1e4], modes="all")
    assert profile["snr_profile"][0] < 1e-12
    assert everything["snr_profile"][0] == pytest.approx(profile["total_charge_snr"], rel=1e-9)
    assert fields.spectral_snr(fields.FieldModel(d=1, L=32, write="dipole"), [1e4], modes="all")["snr_profile"][0] < 1e-12


def test_a_nonlinear_conserved_field_follows_the_renormalized_linear_law():
    """The quartic field keeps the law; its noise and clock are set by the site variance, with no fitted number."""
    from fieldbridge.memory import fields
    model = fields.FieldModel(d=1, L=48, write="charge", amplitude=3.0, T=0.5, g=0.3)
    times = np.geomspace(4.0, 40.0, 6)
    r = fields.run(model, np.random.default_rng(7), times=times, pairs=40, t_equil=20.0)
    ren = r["renormalized"]
    assert ren["r"] == pytest.approx(1.3068, abs=1e-3)
    assert ren["rms_log_residual"] < 0.06
    bare = fields.compare(times, r["simulation"]["snr_profile"], r["exact_same_size"]["snr_profile"], t_min=10.0)
    assert bare["rms_log_residual"] > ren["rms_log_residual"]


def test_without_conservation_the_loss_is_exponential_and_the_quartic_term_only_speeds_it():
    from fieldbridge.memory import fields
    model = fields.FieldModel(d=1, L=64, conserved=False, amplitude=3.0, T=0.5, g=0.3, kappa0=0.2)
    law = fields.predicted_law(model)
    assert law["law"] == "exponential" and law["rate"] == pytest.approx(0.4) and law["rate_is"] == "lower bound"
    assert fields.hartree_mass(model) > model.kappa0
    r = fields.run(model, np.random.default_rng(3), times=np.linspace(0.5, 8.0, 10), pairs=30, t_equil=10.0)
    assert r["rate_simulated"] > law["rate"]


@pytest.mark.parametrize("bad, message", [
    ({"kind": "equations"}, "kind 'field'"),
    ({"field": {"d": 4}}, "field.d"),
    ({"field": {"write": "quadrupole"}}, "field.write"),
    ({"field": {"L": 100000}}, "field.L"),
    ({"field": {"conserved": False, "kappa0": 0.0}}, "kappa0"),
    ({"field": {"mass": 1.0}}, "unknown field entries"),
    ({"field": {"d": 1.5}}, "must be an integer"),
])
def test_field_specifications_outside_the_contract_are_refused(bad, message):
    from fieldbridge.memory import fields
    base = json.loads((EX / "fields" / "conserved_1d_charge.json").read_text(encoding="utf-8"))
    base.update(bad)
    with pytest.raises(spec.SpecError, match=message):
        fields.from_spec(base)


def test_every_field_example_loads():
    from fieldbridge.memory import fields
    paths = sorted((EX / "fields").glob("*.json"))
    assert len(paths) == 5
    for path in paths:
        model = fields.from_spec(json.loads(path.read_text(encoding="utf-8")))
        assert fields.predicted_law(model)["law"] in ("power", "exponential")


def test_regimes_command_writes_a_report(tmp_path):
    out = tmp_path / "regimes"
    cmd = [sys.executable, "-B", "-m", "fieldbridge", "memory", "regimes", "--n", "400", "--out-dir", str(out)]
    subprocess.run(cmd, check=True, cwd=ROOT, capture_output=True)
    report = json.loads((out / "regimes.json").read_text(encoding="utf-8"))
    assert report["novelty_established"] is False and report["implementation_sha256"]
    assert report["plateau_gaussian"] == pytest.approx(0.5 * math.log1p(4.5))
    p = 0.5 * (1 + math.erf(math.sqrt(4.5 / 2)))  # the sign is decided in the linear stage: Phi(W^1/2)
    assert report["binary_plateau"]["accuracy"] == pytest.approx(p)
    assert report["regimes"]["expanding"]["simulation"]["accuracy"][-1] == pytest.approx(p, abs=0.03)


# ------------------------------------------------------------------------------------------------ co-discovery
def test_models_from_different_fields_reach_the_symmetric_write_by_different_derivations():
    from fieldbridge.memory import codiscovery as cd
    rng = np.random.default_rng(1)
    rows = {n: cd.derive_symmetric_write(load(n), rng, check_law=False)
            for n in ("laser", "toggle", "repressor_ring4", "schlogl", "repressilator", "tubes")}
    laser = rows["laser"]
    assert laser["status"] == "reached" and laser["word"] == "SCRK"
    assert laser["write_point"]["value"] == pytest.approx(1.0, abs=1e-3)  # threshold P = gamma kappa / g
    assert laser["class"] == "S(reflection) C R(pitchfork) K"
    assert rows["toggle"]["class"] == "S(exchange) C R(pitchfork) K"
    assert rows["repressor_ring4"]["class"] == "S(cyclic, order 4) C R(pitchfork) K"
    assert rows["schlogl"]["status"] == "reached after unfolding"
    assert rows["schlogl"]["class"] == "C R(fold) U R(pitchfork) K"
    for n in ("laser", "toggle", "repressor_ring4", "schlogl"):
        canon = rows[n]["canonical"]
        assert canon["cubic_sign"] == -1 and canon["even_part"] < 1e-6
    assert rows["repressilator"]["status"] == "obstructed" and "Hopf" in rows["repressilator"]["obstruction"]
    assert rows["tubes"]["class"] == "S(exchange) C R(subcritical pitchfork)"
    assert "subcritical" in rows["tubes"]["obstruction"]


def test_the_write_law_has_the_same_constant_in_a_laser_and_a_toggle():
    from fieldbridge.memory import codiscovery as cd
    rep = cd.codiscover([load("laser"), load("toggle")], np.random.default_rng(5), n_traj=1000)
    s = rep["summary"]
    assert s["reached"] == 2 and s["fields_reached"] == ["laser physics", "synthetic biology"]
    assert len(s["derivation_classes"]) == 2
    assert abs(s["law_constant_mean"] - cd.LAW_CONSTANT) < 3.5 * s["law_constant_stderr"]


def test_the_delay_constant_is_the_first_zero_of_the_airy_derivative():
    from scipy.special import ai_zeros
    from fieldbridge.memory import codiscovery as cd
    assert cd.DELAY_CONSTANT == pytest.approx(-ai_zeros(1)[1][0], abs=1e-12)


def test_a_field_writes_the_landau_model_at_the_coercive_field_with_the_airy_delay():
    from fieldbridge.memory import codiscovery as cd
    row = cd.derive_threshold_write(load("pitchfork"), np.random.default_rng(1))
    assert row["status"] == "reached by a write field"
    assert row["class"] == "S C R(pitchfork) W R(fold) K L"  # the symmetry forbids a threshold along the control
    assert row["write_point"]["value"] == pytest.approx(2 / (3 * math.sqrt(3)), abs=1e-6)  # h_c of x - x^3 - h
    assert row["canonical"]["cubic"] == pytest.approx(-1 / 3, abs=1e-3)
    assert row["law_constant"]["constant"] == pytest.approx(cd.DELAY_CONSTANT, abs=1e-4)


def test_threshold_write_by_the_control_by_a_field_and_its_obstruction():
    from fieldbridge.memory import codiscovery as cd
    rng = np.random.default_rng(1)
    rows = {n: cd.derive_threshold_write(load(n), rng, check_law=False)
            for n in ("schlogl", "toggle_unequal", "colloid_patch", "compartments")}
    # with a = 4.5 and k3 = 5.75 the Schloegl drift is -y^3 + y + (b - 15/8), y = x - 3/2: a Landau model in a field
    assert rows["schlogl"]["class"] == "C R(fold) K"
    assert rows["schlogl"]["write_point"]["value"] == pytest.approx(15 / 8 - 2 / (3 * math.sqrt(3)), abs=1e-6)
    assert rows["toggle_unequal"]["class"] == "C R(fold) K"
    assert rows["colloid_patch"]["class"] == "W R(fold) K"  # the control only rescales the drift
    for n in ("schlogl", "toggle_unequal", "colloid_patch"):
        canon = rows[n]["canonical"]
        assert canon["linear_part"] < 1e-6 and canon["c2"] == pytest.approx(1.0, abs=1e-3)
    assert rows["compartments"]["status"] == "obstructed"
    assert "single stable state" in rows["compartments"]["obstruction"]


OSC = EX / "oscillators"


def test_every_oscillator_example_loads_with_its_field():
    paths = sorted(OSC.glob("*.json"))
    assert len(paths) == 8
    for path in paths:
        real = spec.load(path)
        assert real.spec["field"] and real.control and real.spec["assumptions"]


def test_the_symmetry_of_the_cycle_sets_the_locking_ratio():
    from fieldbridge.memory import phase_locking as pl
    rng = np.random.default_rng(1)
    rows = {n: pl.derive(spec.load(path), rng, check_law=False)
            for n, path in (("van der Pol", OSC / "van_der_pol.json"), ("pumped", OSC / "parametron.json"),
                            ("ring of 3", EX / "repressilator.json"))}
    assert rows["van der Pol"]["class"] == "R(phase) K(1:1)"
    assert rows["pumped"]["class"] == "S(reflection) R(phase) K(2:1)"
    assert rows["ring of 3"]["class"] == "S(cyclic, order 3) R(phase) K(3:1)"
    # the reflection maps the van der Pol cycle onto itself half a period later; the bias reverses under it, so the
    # phase response along the bias has odd harmonics only, while the pumped stiffness keeps it and has even ones only
    h = rows["van der Pol"]["canonical"]["harmonics"]
    assert max(h[1], h[3]) < 1e-6 * h[0] and h[2] > 0.05 * h[0]
    h = rows["pumped"]["canonical"]["harmonics"]
    assert max(h[0], h[2]) < 1e-5 * h[1]
    h = rows["ring of 3"]["canonical"]["harmonics"]
    assert max(h[0], h[1], h[3], h[4]) < 1e-5 * h[2]
    for r in rows.values():
        assert r["canonical"]["deviation"] < 0.05  # averaged drift of the driven realization against K cos(n psi - phi_n)


def test_the_junction_has_a_one_harmonic_phase_response_and_the_adler_width():
    """dphi/dt = I - sin phi: Z(theta) = (I^2 + cos theta + w0 sin theta) / (w0 I), w0 = (I^2 - 1)^(1/2), up to the
    origin of theta; hence |Z_1| = 1 / w0, the mean of Z is I / w0, and no higher harmonic."""
    from fieldbridge.memory import phase_locking as pl
    row = pl.derive(spec.load(OSC / "josephson.json"), np.random.default_rng(1), check_law=True)
    assert row["class"] == "C R(phase) K(1:1) L"
    cur = row["write_point"]["value"]
    assert cur == pytest.approx(1.6)  # middle of the oscillating values 1.2 ... 2.0 on the grid over [0, 2]
    w0 = math.sqrt(cur ** 2 - 1)
    c = row["canonical"]
    assert c["omega"] == pytest.approx(w0, rel=1e-9)
    assert c["coupling_per_unit_drive"] == pytest.approx(1 / (2 * w0), rel=1e-7)
    assert max(c["harmonics"][1:]) < 1e-7 * c["harmonics"][0]
    assert np.mean(c["Zp"]) == pytest.approx(cur / w0, rel=1e-7)
    law = row["law_constant"]
    assert abs(law["constant"] - 1.0) < max(3 * law["stderr"], 2e-3)


def test_phase_locking_obstructions_name_their_reason():
    from fieldbridge.memory import phase_locking as pl
    rng = np.random.default_rng(1)
    lv = pl.derive(spec.load(OSC / "lotka_volterra.json"), rng, check_law=False)
    assert lv["status"] == "obstructed" and "neutral cycles" in lv["obstruction"] and lv["word"] == "R"
    assert "rate of time" in pl.derive(load("colloid_patch"), rng, check_law=False)["obstruction"]
    assert "no oscillation" in pl.derive(load("toggle"), rng, check_law=False)["obstruction"]


def test_codiscover_command_writes_report_and_figure(tmp_path):
    out = tmp_path / "codiscover"
    cmd = [sys.executable, "-B", "-m", "fieldbridge", "memory", "codiscover", str(EX / "laser.json"),
           str(EX / "tubes.json"), "--no-law", "--out-dir", str(out)]
    subprocess.run(cmd, check=True, cwd=ROOT, capture_output=True)
    report = json.loads((out / "codiscover.json").read_text(encoding="utf-8"))
    assert report["novelty_established"] is False and report["implementation_sha256"]
    assert [r["status"] for r in report["rows"]] == ["reached", "obstructed"]
    assert "Classes of derivation" in (out / "codiscover.md").read_text(encoding="utf-8")
    assert (out / "codiscover.png").stat().st_size > 0
    out2 = tmp_path / "threshold"
    cmd = [sys.executable, "-B", "-m", "fieldbridge", "memory", "codiscover", str(EX / "laser.json"),
           str(EX / "compartments.json"), "--target", "threshold-write", "--no-law", "--out-dir", str(out2)]
    subprocess.run(cmd, check=True, cwd=ROOT, capture_output=True)
    report = json.loads((out2 / "codiscover.json").read_text(encoding="utf-8"))
    assert report["summary"]["target_key"] == "threshold-write"
    assert [r["status"] for r in report["rows"]] == ["reached by a write field", "obstructed"]
    assert (out2 / "codiscover.png").stat().st_size > 0
    out3 = tmp_path / "phase"
    cmd = [sys.executable, "-B", "-m", "fieldbridge", "memory", "codiscover", str(OSC / "van_der_pol.json"),
           str(OSC / "lotka_volterra.json"), "--target", "phase-locking", "--no-law", "--out-dir", str(out3)]
    subprocess.run(cmd, check=True, cwd=ROOT, capture_output=True)
    report = json.loads((out3 / "codiscover.json").read_text(encoding="utf-8"))
    assert report["summary"]["ratios"] == {"1:1": ["van der Pol oscillator"]}
    assert [r["status"] for r in report["rows"]] == ["reached", "obstructed"]
    assert (out3 / "codiscover.png").stat().st_size > 0


def test_the_magnet_is_a_single_domain_particle_with_the_astroid():
    """The Stoner-Wohlfarth realization writes by a pitchfork, a fold at the astroid field, or a subcritical
    pitchfork, as the field turns from the hard axis to the easy axis."""
    base = json.loads((EX / "stoner_wohlfarth.json").read_text(encoding="utf-8"))
    for deg, kind in ((90, "supercritical pitchfork"), (60, "saddle-node"), (20, "saddle-node"), (0, "subcritical")):
        s = json.loads(json.dumps(base))
        s["parameters"]["psi"] = float(np.deg2rad(deg))
        real = spec.load(s)
        events = an.locate_writes(real, np.random.default_rng(1), values=np.linspace(0.0, 1.5, 7))
        psi = np.deg2rad(deg)
        astroid = (np.cos(psi) ** (2 / 3) + np.sin(psi) ** (2 / 3)) ** -1.5
        assert abs(events[0]["v"] - astroid) < 1e-4, deg
        nf = an.normal_form(real, np.asarray(events[0]["q"]), "h", events[0]["v"])
        assert nf["kind"].startswith(kind), (deg, nf["kind"])


def test_a_parameter_named_h_is_not_replaced_by_the_write_field():
    """The write field of the analysis has its own reserved name, so a material parameter h is kept."""
    real = spec.load(EX / "stoner_wohlfarth.json")
    rw = an.with_write_field(real, np.array([1.0]), (0.0, 1.0))
    assert rw.params["h"] == real.params["h"] and rw.params[an.WRITE_FIELD] == 0.0
    changed = json.loads((EX / "stoner_wohlfarth.json").read_text(encoding="utf-8"))
    changed["parameters"][an.WRITE_FIELD] = 0.1
    with pytest.raises(spec.SpecError):
        spec.load(changed)
