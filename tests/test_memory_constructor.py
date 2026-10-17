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
    assert h["phase_observable"] == "cycle_counted_crossing_time" and h["phase_topology"] == "unwrapped"
    sep, var, t = np.array(h["separation"]), np.array(h["variance"]), np.array(h["time"])
    assert abs(sep[-1] - sep[4]) < 0.1 and abs(sep[-1] - np.pi / 2) < 0.15  # the written shift persists
    assert np.std(var[4:] / t[4:]) / np.mean(var[4:] / t[4:]) < 0.2  # variance linear in time: Law 2
    L = phase.lock(real, cyc, rng, strengths=[0.01, 0.03, 0.1], target=np.pi / 2, D=0.005, n_traj=80)
    assert L["retention_criterion"] == "unwrapped_phase_variance_equals_squared_shift"
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


def square_amplitude(split):
    """The amplitude equation of a mode that the symmetry of a square (D4) makes degenerate,
    dz/dt = eps z - |z|^2 z + g zbar^3 with z = x + i y; the eigenvalue along y is lowered by an imposed split."""
    return spec.load({
        "schema": "fieldbridge-memory/1", "kind": "equations", "name": f"amplitude with square symmetry, split {split:g}",
        "question": "Which states does a degenerate critical mode store?", "assumptions": ["cubic normal form"],
        "provenance": {"source": "M. Golubitsky, I. Stewart and D. G. Schaeffer, Singularities and Groups in "
                                 "Bifurcation Theory, Vol. II (Springer, 1988)"},
        "carrier": {"kind": "euclid", "variables": ["x", "y"], "scale": 1.0},
        "parameters": {"eps": 0.5, "g": 0.3, "d": split},
        "drift": {"x": "eps*x - (x**2 + y**2)*x + g*(x**3 - 3*x*y**2)",
                  "y": "(eps - d)*y - (x**2 + y**2)*y + g*(y**3 - 3*x**2*y)"},
        "control": {"name": "eps", "range": [-0.6, 1.0]}, "noise": 0.01, "dt": 0.05})


@pytest.mark.parametrize("split, dimension", [(0.0, 2), (1e-4, 2), (0.05, 1)])
def test_eigenvalues_that_cross_zero_together_are_counted(split, dimension):
    """The continuation does not separate write points closer than 2e-3 of the range of the control (3.2e-3 here,
    where the critical eigenvalue grows at rate 1). Two eigenvalues split by less cross together, and the critical
    eigenspace is a plane; a larger split leaves one critical direction and the pitchfork along x."""
    real = square_amplitude(split)
    ev = next(e for e in events(real) if abs(e["v"]) < 1e-3)
    nf = an.normal_form(real, ev["q"], "eps", ev["v"])
    ob = construct.obstruction(nf)
    assert nf["critical_dimension"] == dimension
    if dimension == 2:
        assert nf["kind"].startswith("degenerate critical mode") and ob["degenerate"]
        assert not ob["realizes_generator"] and not (ob["asymmetric"] or ob["biased"] or ob["subcritical"])
    else:
        assert nf["kind"].startswith("supercritical pitchfork") and ob["realizes_generator"]


def test_a_degenerate_pair_in_rotated_coordinates_is_not_a_hopf_point():
    """In coordinates rotated by 30 degrees the degenerate pair of the square symmetry is not diagonal, and rounding
    gives it imaginary parts of order 1e-13 while the whole spectrum is critical; the write point is still a degenerate
    critical mode of two directions, as without the rotation (found by the held-out set of the discovery engine)."""
    c, s = math.cos(math.pi / 6), math.sin(math.pi / 6)
    x, y = f"({c}*u - {s}*w)", f"({s}*u + {c}*w)"
    fx = f"eps*{x} - ({x}**2 + {y}**2)*{x} + g*({x}**3 - 3*{x}*{y}**2)"
    fy = f"eps*{y} - ({x}**2 + {y}**2)*{y} + g*({y}**3 - 3*{x}**2*{y})"
    real = spec.load({
        "schema": "fieldbridge-memory/1", "kind": "equations", "name": "amplitude with square symmetry, rotated",
        "question": "Which states does a degenerate critical mode store?", "assumptions": ["cubic normal form"],
        "provenance": {"source": "M. Golubitsky, I. Stewart and D. G. Schaeffer, Singularities and Groups in "
                                 "Bifurcation Theory, Vol. II (Springer, 1988)"},
        "carrier": {"kind": "euclid", "variables": ["u", "w"], "scale": 1.0},
        "parameters": {"eps": 0.5, "g": 0.3},
        "drift": {"u": f"{c}*({fx}) + {s}*({fy})", "w": f"-{s}*({fx}) + {c}*({fy})"},
        "control": {"name": "eps", "range": [-0.6, 1.0]}, "noise": 0.01, "dt": 0.05})
    ev = next(e for e in events(real) if abs(e["v"]) < 1e-3)
    nf = an.normal_form(real, ev["q"], "eps", ev["v"])
    assert not nf["hopf"] and nf["kind"].startswith("degenerate critical mode") and nf["critical_dimension"] == 2


def test_a_degenerate_critical_mode_stops_the_derivation_at_the_reduction():
    from fieldbridge.memory import codiscovery as cd
    rng = np.random.default_rng(1)
    square = cd.derive_symmetric_write(square_amplitude(0.0), rng, check_law=False)
    split = cd.derive_symmetric_write(square_amplitude(0.05), rng, check_law=False)
    assert square["status"] == "obstructed" and square["word"] == "CR"
    assert square["obstruction_short"] == "degenerate critical mode: 2 directions"
    reduction = square["steps"][-1]
    assert reduction["letter"] == "R" and reduction["critical_dimension"] == 2
    assert reduction["eliminated_directions"] == 0
    assert square["route"].endswith("R[degenerate critical mode: 2 directions]")
    assert split["status"] == "reached" and split["class"] == "S(reflection) C R(pitchfork) K"
    assert "critical_dimension" not in split["steps"][2] and split["steps"][2]["eliminated_directions"] == 1


def test_the_card_names_a_degenerate_critical_mode_and_its_four_states():
    """The square symmetry leaves four stable states on the axes (0 < g < 1), not the two of a pitchfork; the
    structural prediction of a pitchfork concerns a single reversed mode and is not tested at this write point."""
    from fieldbridge.memory import discovery
    real = square_amplitude(0.0)
    card = discovery.evaluate(real, np.random.default_rng(1), quick=True)
    assert card["states"]["count"] == 4
    assert discovery.mechanism_class(card).startswith("write into a degenerate critical mode")
    assert card["verdict"]["constructs"].startswith("degenerate critical mode at eps = 0: 2 directions")
    comparison = predict.compare(predict.predict(real, np.random.default_rng(0)), card)
    assert comparison["all_consistent"] and comparison["write"]["tested_write_points"] == 0


def test_the_card_figure_says_whether_the_field_removes_the_state_or_destabilizes_it():
    """At a fold the stored state disappears. At a transcritical point it remains an equilibrium and loses stability,
    as for the capillary rotors, each at 90 degrees to its target angle, where the field exerts no torque."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from fieldbridge.memory import visual
    branch = [{"h": h, "kappa_min": 1.0 - h, "overlap": -1.0} for h in (0.0, 0.5, 0.99)]
    labels = {}
    for kind in ("saddle-node (fold): a state appears or vanishes", "transcritical: two states exchange stability"):
        card = {"construct": {"events": [{"source": "write field", "kind": kind, "value": 1.0}],
                              "write_field": {"branch": branch, "h_c": 1.0}}}
        fig, ax = plt.subplots()
        visual.panel_bifurcation(ax, card)
        labels[kind.split(":")[0]] = [t.get_text() for t in ax.texts]
        plt.close(fig)
    assert labels == {"saddle-node (fold)": ["state 1 ceases to exist\nat h = 1"],
                      "transcritical": ["state 1 loses stability\nat h = 1"]}


def ring_of_cells(n, c=0.2):
    """A ring of n convection cells (three modes each) sharing heat-conducting walls. Wall i, between cells i and
    i + 1, adds -c (y_i + y_(i+1)) to both, so a mode y_j ~ cos(q j) is damped by c (2 + 2 cos q): an odd ring has no
    alternating mode, and its critical mode is the degenerate pair q = pi (n -+ 1)/n, at R = 1 + 2c (1 - cos(pi/n))."""
    drift = {}
    for i in range(n):
        left, right = (i - 1) % n, (i + 1) % n
        drift[f"x{i}"] = f"p*(y{i} - x{i})"
        drift[f"y{i}"] = f"R*x{i} - x{i}*z{i} - y{i} - c*(2*y{i} + y{left} + y{right})"
        drift[f"z{i}"] = f"x{i}*y{i} - z{i} - eta*(2*z{i} - z{left} - z{right})"
    return spec.load({
        "schema": "fieldbridge-memory/1", "kind": "equations", "name": f"ring of {n} convection cells",
        "question": "Which states does a ring of convection cells store?", "assumptions": ["three modes per cell"],
        "provenance": {"source": "B. Saltzman, J. Atmos. Sci. 19, 329 (1962); E. N. Lorenz, J. Atmos. Sci. 20, "
                                 "130 (1963); walls conducting heat between neighbouring cells"},
        "carrier": {"kind": "euclid", "variables": [f"{a}{i}" for a in "xyz" for i in range(n)], "scale": 1.0},
        "parameters": {"R": 1.3, "p": 10.0, "c": c, "eta": 0.2}, "drift": drift,
        "control": {"name": "R", "range": [0.5, 2.0]}, "noise": 0.001, "dt": 0.002})


def test_a_ring_of_five_cells_is_written_into_a_degenerate_pair_found_by_the_card_and_the_derivation():
    """Fifteen variables: the card searches the control as the continuation C does, above twelve as below."""
    from fieldbridge.memory import codiscovery as cd
    real = ring_of_cells(5)
    row = cd.derive_symmetric_write(real, np.random.default_rng(1), check_law=False)
    card = construct.construct(real, np.random.default_rng(1), check_write=False)
    (ev,) = [e for e in card["events"] if e["source"] == "control"]
    assert ev["value"] == pytest.approx(1 + 0.4 * (1 - math.cos(math.pi / 5)), abs=1e-4)
    assert ev["critical_dimension"] == 2 and ev["obstruction"]["degenerate"]
    assert [w["value"] for w in row["steps"][0]["write_points"]] == pytest.approx([ev["value"]], abs=1e-9)
    assert row["status"] == "obstructed" and row["obstruction_short"] == "degenerate critical mode: 2 directions"
    assert row["steps"][1]["eliminated_directions"] == 13


def test_above_six_variables_a_ring_is_turned_and_reflected_cell_by_cell():
    """The symmetries of a ring of four cells act on the three variables of every cell together: the group of the
    square, with or without the reversal of every circulation (x, y -> -x, -y), 16 elements with the identity. Shifts
    of the list of variables, the candidates tried before, are not among them. The reversal and the turn by one cell
    both reverse the alternating mode that becomes critical at R = 1, so the derivation starts with S."""
    from fieldbridge.memory import codiscovery as cd
    real = ring_of_cells(4)
    group = predict.symmetries(real, np.random.default_rng(0))
    elements = {(tuple(g["permutation"]), tuple(g["signs"])) for g in group}
    assert len(group) == 15 and len(elements) == 15
    turn = tuple(4 * b + (i + 1) % 4 for b in range(3) for i in range(4))
    assert (turn, (1,) * 12) in elements
    assert (tuple(range(12)), (-1,) * 8 + (1,) * 4) in elements
    assert (tuple(int(i) for i in np.roll(np.arange(12), 1)), (1,) * 12) not in elements
    J = an.jacobian(real, np.zeros(12), {"R": 1.0})
    ev, vecs = np.linalg.eig(J)
    mode = np.real(vecs[:, np.argmin(np.abs(ev))])
    sym = cd.reversing_symmetry(real, np.random.default_rng(0), mode)
    assert sym is not None and cd.describe_symmetry(real, sym).startswith("reflection x0 -> -x0")
    row = cd.derive_symmetric_write(real, np.random.default_rng(1), check_law=False)
    assert row["status"] == "reached" and row["class"] == "S(reflection) C R(pitchfork) K"
    assert row["write_point"]["value"] == pytest.approx(1.0, abs=1e-4)


def test_the_write_law_has_the_same_constant_in_a_laser_and_a_toggle():
    from fieldbridge.memory import codiscovery as cd
    rep = cd.codiscover([load("laser"), load("toggle")], 5, n_traj=1000)
    s = rep["summary"]
    assert s["reached"] == 2 and s["fields_reached"] == ["laser physics", "synthetic biology"]
    assert len(s["derivation_classes"]) == 2
    assert abs(s["law_constant_mean"] - cd.LAW_CONSTANT) < 3.5 * s["law_constant_stderr"]


def test_realizations_that_simulate_the_same_equations_count_as_one_model():
    """The unequal toggle, unfolded to equal promoters, is the toggle: its law constant is a replicate of the
    toggle's, not a second model. The laser reduces to the same canonical form by other equations."""
    from fieldbridge.memory import codiscovery as cd
    rep = cd.codiscover([load(n) for n in ("toggle", "toggle_unequal", "laser")], 3, n_traj=400, replicates=4)
    s = rep["summary"]
    assert s["law_constant_realizations"] == 3 and s["law_constant_dof"] == 2
    assert s["same_equations"] == [[load("toggle").name, load("toggle_unequal").name]]


def test_the_error_of_a_swept_write_comes_from_independent_runs():
    from fieldbridge.memory import analysis as an, construct
    real = load("toggle")
    ev = an.locate_writes(real, np.random.default_rng(7), n_starts=24)[0]
    nf = an.normal_form(real, ev["q"], real.control, ev["v"])
    law = construct.swept_write_check(real, nf, np.random.default_rng(1), n_traj=800, replicates=4)
    assert law["replicates"] == 4 and len(law["runs"]) == 4 and law["trajectories_per_run"] == 200
    assert law["measured"] == pytest.approx(np.mean(law["runs"]))
    assert law["stderr"] == max(law["stderr_binomial"], law["stderr_runs"])
    # the runs use streams spawned from the generator, not the generator itself: a few numbers drawn from it first
    # give other runs, not the same noise assigned to the trajectories in another order
    shifted = np.random.default_rng(1)
    shifted.standard_normal(6)
    other = construct.swept_write_check(real, nf, shifted, n_traj=800, replicates=4)
    assert other["runs"] != law["runs"]


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


def test_the_delay_fit_recovers_the_constant_when_the_fold_is_mislocated():
    """An error delta of the located fold adds delta r^(-2/3) to the measured constant. A quadratic extrapolation
    over the five slowest rates, without this term, is then biased; the fit of delay_constant is not."""
    from fieldbridge.memory import codiscovery as cd
    r = np.array(cd.DELAY_RATES)
    e = r ** (1 / 3)
    clean = cd.DELAY_CONSTANT - 0.0935 * e + 0.0411 * e ** 2 - 0.0269 * e ** 3 + 0.0179 * e ** 4   # the Landau model
    exact = cd.delay_constant(r, clean)
    assert abs(exact["constant"] - cd.DELAY_CONSTANT) < 1e-9 and abs(exact["fold_offset"]) < 1e-12
    assert exact["fit_terms"]["powers"] == [0, 1, 2, 3, 4, -2] and exact["stderr"] == cd.DELAY_FLOOR
    delta = 9e-9
    shifted = cd.delay_constant(r, clean + delta * e ** -2.0)
    assert abs(shifted["constant"] - cd.DELAY_CONSTANT) < 1e-9
    assert shifted["fold_offset"] == pytest.approx(delta, rel=1e-6)
    slow = np.argsort(e)[:5]
    biased = np.polyfit(e[slow], (clean + delta * e ** -2.0)[slow], 2)[2]
    assert biased - cd.DELAY_CONSTANT > 3e-4          # the extrapolation without the term
    # with five rates the polynomial has first order, and the uncertainty covers the truncation
    few = cd.delay_constant(r[-5:], clean[-5:])
    assert few["fit_terms"]["powers"] == [0, 1, -2] and abs(few["constant"] - cd.DELAY_CONSTANT) < few["stderr"]


def test_a_fold_shifted_by_hand_appears_as_the_offset_of_the_delay_law():
    """The Schloegl reactor swept through its fold, with the fold deliberately placed at a wrong value of the
    control: the crossings move by the same mu, the fit reports the shift, and the constant is unchanged."""
    from fieldbridge.memory import analysis as an
    from fieldbridge.memory import codiscovery as cd
    real = load("schlogl")
    ev = next(e for e in an.locate_writes(real, np.random.default_rng(1), n_starts=24))
    q_f, p_f = cd.refine_fold(real, real.control, ev["q"], ev["v"])
    canon = cd.canonical_fold(real, real.control, q_f, p_f, an.normal_form(real, q_f, real.control, p_f))
    right = cd.fold_delay_law(real, real.control, q_f, p_f, canon)
    assert abs(right["constant"] - cd.DELAY_CONSTANT) < max(right["stderr"], 1e-6)
    delta = 9e-9                                         # in canonical units: mu = a2 b (p - p_f)
    wrong = cd.fold_delay_law(real, real.control, q_f, p_f - delta / canon["mu_per_unit_parameter"], canon)
    # compared with the sweeps through the located fold, so that the accuracy of its location does not enter
    assert wrong["fold_offset"] - right["fold_offset"] == pytest.approx(delta, rel=0.05)
    assert abs(wrong["constant"] - right["constant"]) < 1e-6
    slowest = min(wrong["rows"], key=lambda row: row["rate"])
    assert slowest["constant"] - min(right["rows"], key=lambda row: row["rate"])["constant"] == pytest.approx(
        delta * slowest["rate"] ** (-2 / 3), rel=0.05)


def _quintic_fold(A=0.14, c=1.2):
    """dx/dt = h - g(x), g' = A (x^2 - 1)(x^2 - c^2). The occupied state, -1 < x < 1, exists for |h| < g(1) and
    disappears at the fold x = 1, h = g(1); a second stable state, x < -c, exists for h < g(-c). With A = 0.14 the
    coefficient a2 b of the canonical form is 0.062, small as for a pair of states of the caged dipoles."""
    from fieldbridge.memory import codiscovery as cd
    base = json.loads((EX / "pitchfork.json").read_text(encoding="utf-8"))
    base.pop("potential")
    base.update(name="quintic fold", parameters={"h": 0.0, "A": A, "c": c},
                drift={"x": "h - A*(x**5/5 - (1 + c**2)*x**3/3 + c**2*x)"}, control={"name": "h", "range": [-0.3, 0.3]})
    real = spec.load(base)
    g1 = A * (1 / 5 - (1 + c ** 2) / 3 + c ** 2)
    q_f, p_f = cd.refine_fold(real, "h", np.array([0.999]), g1 - 1e-6)
    canon = cd.canonical_fold(real, "h", q_f, p_f, an.normal_form(real, q_f, "h", p_f))
    return real, q_f, p_f, canon, g1


def test_a_sweep_that_would_start_beyond_the_branch_of_the_occupied_state_is_left_out(monkeypatch):
    """The sweeps start at mu = -5 r^(2/3), far from the fold in the parameter where a2 b is small. For the two
    fastest rates the start lies below h = -g(1), where the branch of the occupied state has ended at its other fold,
    and the root search converges on the second stable state, on the side of the occupied state. Without the test of
    the branch these sweeps are used; with it they are left out, and the other sweeps and the fit are those of a run
    without the two rates."""
    from fieldbridge.memory import codiscovery as cd
    real, q_f, p_f, canon, g1 = _quintic_fold()
    assert p_f == pytest.approx(g1, abs=1e-9) and canon["mu_per_unit_parameter"] == pytest.approx(0.44 * 0.14)
    beyond = [r for r in cd.DELAY_RATES if p_f - 5 * r ** (2 / 3) / canon["mu_per_unit_parameter"] < -g1]
    assert beyond == list(cd.DELAY_RATES[:2])
    law = cd.fold_delay_law(real, "h", q_f, p_f, canon)
    assert [x["rate"] for x in law["left_out"]] == sorted(beyond)
    assert [x["reason"] for x in law["left_out"]] == ["the branch of the occupied state ends before the start",
                                                      "a slower sweep was left out"]
    alone = cd.fold_delay_law(real, "h", q_f, p_f, canon, rates=cd.DELAY_RATES[2:])
    assert alone["left_out"] == [] and law["rows"] == alone["rows"] and law["constant"] == alone["constant"]
    monkeypatch.setattr(cd, "off_branch", lambda *args, **kwargs: None)
    unchecked = cd.fold_delay_law(real, "h", q_f, p_f, canon)
    assert [row["rate"] for row in unchecked["rows"]] == list(cd.DELAY_RATES)  # the root search found a stable state
    assert unchecked["rows"][2:] == law["rows"]


def test_a_sweep_left_out_ends_the_series_so_that_the_fit_rests_on_the_slowest_rates(monkeypatch):
    """The start of the sweep at r = 3.2e-5 is refused by hand. The faster sweep is left out with it rather than
    leaving a gap, and the constant is fitted over the five slower sweeps, which are unchanged."""
    from fieldbridge.memory import codiscovery as cd
    real, q_f, p_f, canon, _ = _quintic_fold()
    rates = cd.DELAY_RATES[2:]
    full = cd.fold_delay_law(real, "h", q_f, p_f, canon, rates=rates)
    refused = p_f - 5 * rates[1] ** (2 / 3) / canon["mu_per_unit_parameter"]
    original = cd.off_branch

    def refuse(real_, param, q0, p0, *args, **kwargs):
        if abs(p0 - refused) < 1e-9:
            return "refused by the test"
        return original(real_, param, q0, p0, *args, **kwargs)

    monkeypatch.setattr(cd, "off_branch", refuse)
    law = cd.fold_delay_law(real, "h", q_f, p_f, canon, rates=rates)
    assert law["left_out"] == [{"rate": rates[1], "reason": "refused by the test"},
                               {"rate": rates[0], "reason": "a slower sweep was left out"}]
    assert law["rows"] == full["rows"][2:]
    fit = cd.delay_constant([row["rate"] for row in law["rows"]], [row["constant"] for row in law["rows"]])
    assert law["constant"] == fit["constant"] and law["stderr"] == fit["stderr"]


def test_a_root_search_for_the_fold_that_stops_early_is_restarted(monkeypatch):
    """A root search can stop before it converges, and its result can still pass the test of the residual: on one
    machine this displaced the fold of the rotor patch by 9e-9 in canonical units. Here the first search for the
    fold of the Schloegl reactor is cut off after a few evaluations, and the search is repeated from there."""
    import scipy.optimize
    from fieldbridge.memory import analysis as an
    from fieldbridge.memory import codiscovery as cd
    real = load("schlogl")
    ev = next(e for e in an.locate_writes(real, np.random.default_rng(1), n_starts=24))
    exact = 15 / 8 - 2 / (3 * math.sqrt(3))
    original, searches = scipy.optimize.root, []

    def cut_off(fun, x0, **kwargs):
        if not searches:
            kwargs["options"] = {**kwargs["options"], "maxfev": 4}
        searches.append(original(fun, x0, **kwargs))
        return searches[-1]

    monkeypatch.setattr(scipy.optimize, "root", cut_off)
    q_f, p_f = cd.refine_fold(real, real.control, ev["q"], ev["v"])
    assert not searches[0].success and abs(searches[0].x[-1] - exact) > 1e-10   # 1.5e-6 where it was cut off
    assert len(searches) > 1 and abs(p_f - exact) < 1e-12


def test_the_delay_constant_of_the_rotor_patch_lies_within_its_uncertainty():
    from fieldbridge import site_data as sd
    from fieldbridge.memory import codiscovery as cd
    law = cd.derive_threshold_write(load("colloid_patch"), sd._rng("colloid_patch", "threshold-write"))["law"]
    assert len(law["rows"]) == 9
    assert abs(law["constant"] - cd.DELAY_CONSTANT) < max(law["stderr"], 1e-6)


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


def test_adding_a_specification_leaves_the_random_streams_of_the_other_derivations_unchanged(tmp_path, monkeypatch):
    """Each realization draws the stream of its own material. With one stream shared in the order of the files, a
    specification that sorted first re-drew every later realization, and the law constants of the default run of
    Module 9 moved although no equation had changed."""
    import argparse
    import shutil
    from fieldbridge.memory import cli, codiscovery as cd
    # the derivation itself: the runs of the law, spawned from the stream of the material, are the same
    alone = cd.codiscover([load("pitchfork")], 7, n_traj=64, replicates=2)["rows"][0]
    behind = cd.codiscover([load("laser"), load("pitchfork")], 7, n_traj=64, replicates=2)["rows"][1]
    assert alone["law"]["runs"] == behind["law"]["runs"] and alone["law_constant"] == behind["law_constant"]
    drawn = {}

    def derive(real, rng, check_law=True, n_traj=0, replicates=1):
        drawn[real.name] = rng.random(4).tolist()
        return {"name": real.name, "field": "", "status": "obstructed", "obstruction": "no control parameter",
                "obstruction_short": "no control parameter", "word": "", "route": "", "class": "", "steps": []}

    monkeypatch.setitem(cd.TARGETS, "symmetric-write", derive)
    folder = tmp_path / "memory"
    folder.mkdir()
    for name in ("schlogl.json", "toggle.json"):
        shutil.copy(EX / name, folder / name)

    def run():
        drawn.clear()
        cli.cmd_codiscover(argparse.Namespace(specs=[], examples=str(folder), target="symmetric-write", seed=20260923,
                                              no_law=True, trajectories=400, replicates=8,
                                              out_dir=str(tmp_path / "out"), no_figure=True))
        return dict(drawn)

    before = run()
    shutil.copy(EX / "laser.json", folder / "a_laser.json")  # sorts before the other two
    after = run()
    assert set(after) == set(before) | {"single-mode laser"}
    assert all(after[name] == before[name] for name in before)
    assert len({tuple(d) for d in after.values()}) == 3  # different materials draw different streams


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


def test_codiscovery_markdown_without_reached_realization(tmp_path):
    """With no realization reached, the report and its figure name the obstruction instead of failing on absent
    invariants."""
    pytest.importorskip("matplotlib")
    from fieldbridge.memory import codiscovery
    from fieldbridge.memory.visual import codiscovery_figure
    real = spec.load(EX / "oscillators" / "van_der_pol.json")
    for target in ("symmetric-write", "threshold-write"):
        rep = codiscovery.codiscover([real], 0, target=target, check_law=False)
        assert rep["summary"]["reached"] == 0
        md = codiscovery.markdown(rep)
        assert "No realization reaches the target" in md
        assert rep["rows"][0]["obstruction"] in md
        codiscovery_figure(rep, tmp_path / f"{target}.png")
        assert (tmp_path / f"{target}.png").exists()


# the Kerr parametric oscillator of a Kerr-cat qubit at zero detuning (Puri et al., npj Quantum Inf. 3, 18 (2017)):
# mean field and vacuum noise of the loss
KPO = {
    "schema": "fieldbridge-memory/1", "question": "Does a Kerr parametric oscillator at zero detuning write as a pitchfork?",
    "provenance": {"source": "Puri et al., npj Quantum Inf. 3, 18 (2017)", "kind": "authored test"},
    "assumptions": ["single mode, rotating frame of the pump"], "field": "superconducting quantum circuits",
    "kind": "equations", "name": "Kerr parametric oscillator at zero detuning",
    "carrier": {"kind": "euclid", "variables": ["x", "y"], "scale": 1},
    "parameters": {"K": 5.0, "Delta": 0.0, "kappa": 1.0, "p": 6.0},
    "drift": {"x": "(K*(x**2 + y**2) + p + Delta)*y - 0.5*kappa*x",
              "y": "-(K*(x**2 + y**2) - p + Delta)*x - 0.5*kappa*y"},
    "control": {"name": "p", "range": [0.0, 12.0]}, "noise": 0.125, "dt": 0.005, "material": ["Delta"]}


def _write_point(real, value):
    ev = min(an.locate_writes(real, np.random.default_rng(0), n_starts=24), key=lambda e: abs(e["v"] - value))
    return an.normal_form(real, ev["q"], real.control, ev["v"])


def test_vanishing_cubic_is_an_obstruction():
    """At zero detuning the Kerr term only rotates the amplitude: the reduced drift at the threshold p = kappa/2 has
    no cubic term and saturates at fifth order, -K^2/kappa s^5 (the amplitude grows as (p - p_c)^(1/4)). The
    symmetric write is obstructed and no law constant is reported; a detuning restores a negative cubic term."""
    from fieldbridge.memory import codiscovery
    real = spec.load(KPO)
    nf = _write_point(real, 0.5)
    assert nf["value"] == pytest.approx(0.5, abs=1e-4)
    assert nf["kind"].startswith("degenerate cubic")
    assert nf["leading_term"]["order"] == 5
    assert nf["leading_term"]["coefficient"] == pytest.approx(-25.0, rel=1e-3)
    ob = construct.obstruction(nf)
    assert ob["degenerate_cubic"] and not ob["realizes_generator"]
    assert construct.swept_write_check(real, nf, np.random.default_rng(0)) == {}
    row = codiscovery.derive_symmetric_write(real, np.random.default_rng(0), n_traj=40)
    assert row["status"] == "obstructed" and not row.get("law_constant")
    assert row["class"] == "S(reflection) C R(degenerate cubic)"
    assert row["obstruction_short"] == "vanishing cubic: leading term of order 5"
    detuned = _write_point(spec.load(dict(KPO, parameters=dict(KPO["parameters"], Delta=1.0))), 1.118)
    assert detuned["value"] == pytest.approx(math.sqrt(1.25), abs=1e-4)
    assert detuned["kind"].startswith("supercritical pitchfork") and detuned["a3"] < 0


def _second_order_bias(b):
    return {"schema": "fieldbridge-memory/1", "question": "Does a sweep that biases the write at second order offset it?",
            "provenance": {"source": "pitchfork normal form with a bias b eps^2", "kind": "authored test"},
            "assumptions": ["one variable"], "field": "test", "kind": "equations",
            "name": f"pitchfork with a bias of second order in the sweep (b = {b})",
            "carrier": {"kind": "euclid", "variables": ["x"], "scale": 1}, "parameters": {"eps": -0.5, "b": b},
            "drift": {"x": "eps*x - x**3 + b*eps**2"}, "control": {"name": "eps", "range": [-1.0, 1.0]},
            "noise": 0.01, "dt": 0.01}


def test_swept_write_law_separates_the_offset_from_the_slope():
    """A bias b eps^2 vanishes with its derivative at the write point eps = 0, which therefore realizes the generator;
    along the sweep it favours the sign of b whatever the bias of the test, as a straight sweep through a cusp does.
    Without a reversing symmetry the law is tested with the bias toward both states: the constants of the two signs
    are slope + offset and slope - offset, the slope is pi^(1/4) and the offset has the sign of b. With b = 0 the
    reflection x -> -x forbids an offset and one sign is tested."""
    from fieldbridge.memory import codiscovery
    sym = codiscovery.derive_symmetric_write(spec.load(_second_order_bias(0.0)), np.random.default_rng(3), n_traj=800,
                                             replicates=4)
    assert sym["word"] == "SCRKL" and "offset" not in sym["law_constant"] and "by_sign" not in sym["law"]
    for b in (0.3, -0.3):
        row = codiscovery.derive_symmetric_write(spec.load(_second_order_bias(b)), np.random.default_rng(3),
                                                 n_traj=800, replicates=4)
        c = row["law_constant"]
        assert row["word"] == "CRKL" and [s["sign"] for s in row["law"]["by_sign"]] == [1, -1]
        assert c["offset_resolved"] and np.sign(c["offset"]) == np.sign(b) and abs(c["offset"]) > 0.5
        assert abs(c["constant"] - codiscovery.LAW_CONSTANT) < 4 * c["stderr"]
        k_plus, k_minus = (s["constant"] for s in c["by_sign"])
        assert c["constant"] == pytest.approx((k_plus + k_minus) / 2) and c["offset"] == pytest.approx((k_plus - k_minus) / 2)
        assert any("offset" in note for note in row["notes"])


def test_law_constant_of_both_signs_and_of_a_saturated_test():
    """The slope and the offset from the accuracies of the two signs, and no constant where every trajectory of a
    test ended on one side (the accuracy then only bounds the constant)."""
    from statistics import NormalDist
    from fieldbridge.memory import codiscovery
    nd = NormalDist()
    law = {"rate": 0.02, "D_along_mode": 1e-3, "h_along_mode": 0.01, "measured": 0.75, "stderr": 0.01,
           "by_sign": [{"sign": 1, "measured": 0.9, "stderr": 0.01}, {"sign": -1, "measured": 0.6, "stderr": 0.01}]}
    scale = math.sqrt(1e-3) * 0.02 ** 0.25 / 0.01
    k = [nd.inv_cdf(p) * scale for p in (0.9, 0.6)]
    e = [0.01 / nd.pdf(nd.inv_cdf(p)) * scale for p in (0.9, 0.6)]
    c = codiscovery.law_constant(law)
    assert c["constant"] == pytest.approx((k[0] + k[1]) / 2) and c["offset"] == pytest.approx((k[0] - k[1]) / 2)
    assert c["stderr"] == pytest.approx(math.hypot(*e) / 2) and c["offset_resolved"]
    assert codiscovery.law_constant(dict(law, by_sign=[law["by_sign"][0], dict(law["by_sign"][1], measured=1.0)])) == {}
    single = {k_: v for k_, v in law.items() if k_ != "by_sign"}
    assert codiscovery.law_constant(dict(single, measured=0.0, stderr=0.0)) == {}
    assert codiscovery.law_constant(single)["constant"] == pytest.approx(nd.inv_cdf(0.75) * scale)
