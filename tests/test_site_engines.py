"""The calculations of the web page against FieldBridge: the browser engines in fieldbridge/web reproduce the Python
results for the realizations the page shows."""
import copy
import json
import shutil
import subprocess
from pathlib import Path

import pytest

np = pytest.importorskip("numpy")
pytest.importorskip("scipy")
pytest.importorskip("sympy")

from fieldbridge import site_data as sd  # noqa: E402
from fieldbridge import site_registry as reg  # noqa: E402
from fieldbridge.web_models import model_record  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
NODE = shutil.which("node")
needs_node = pytest.mark.skipif(NODE is None, reason="node is not installed")
QUANTUM = [n["id"] for n in reg.NODES
           if n.get("family") == "unitary" or "attach" in n or n.get("base") in
           ("two_spins", "two_spins_h0", "spin1_transverse", "nmr_collective", "nmr_chain", "module14_chain",
            "heteronuclear_spins")]


def run(job, tmp_path):
    f = tmp_path / "job.json"
    f.write_text(json.dumps(job, default=float), encoding="utf-8")
    done = subprocess.run([NODE, str(ROOT / "tests" / "site_engine_harness.cjs"), str(f)], capture_output=True,
                          text=True, check=True)
    return json.loads(done.stdout)


@pytest.fixture(scope="module")
def quantum():
    nodes = sd.resolve(only=QUANTUM)
    records = {}
    for nid in sorted(nodes, key=lambda n: nodes[n]["def"].get("parent") is not None):
        pid = nodes[nid]["def"].get("parent")
        parent = records.get(pid) or (sd.unitary_record(sd.resolve(only=[pid])[pid]) if pid else None)
        if pid and pid not in records:
            records[pid] = parent
        records[nid] = sd.unitary_record(nodes[nid], parent)
    return nodes, records


# ------------------------------------------------------------------------------------------------ quantum
def test_observable_closure_is_the_chapter_11_calculation():
    """The closure the page reports is the symbolic calculation of Chapter 11 (verification.quantum_closure)."""
    from fieldbridge.quantum import language as ql
    from fieldbridge.verification import quantum_closure
    chapter = json.loads((ROOT / "examples/construction/quantum_correlations.json").read_text(encoding="utf-8"))
    z_first = {**chapter, "observable": [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, -1, 0], [0, 0, 0, -1]]}
    assert quantum_closure(chapter)["observable_dimension"] == 2
    assert quantum_closure(z_first)["observable_dimension"] == 1
    two = json.loads((ROOT / "examples/quantum/two_spins.json").read_text(encoding="utf-8"))
    two["parameters"]["h"] = 0.0
    real = ql.load(two)
    assert sd.observable_closure(real.H, real.O) == 2
    real_z = ql.load({**two, "observable": [{"coefficient": 1, "operator": "Z0"}]})
    assert sd.observable_closure(real_z.H, real_z.O) == 1


@needs_node
def test_browser_eigendecomposition_matches_numpy(tmp_path):
    rng = np.random.default_rng(3)
    for n in (2, 3, 5, 8, 16):
        A = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
        H = A + A.conj().T
        got = run({"op": "eigh", "matrix": sd._mat(H)}, tmp_path)
        assert np.allclose(got, np.linalg.eigvalsh(H), atol=1e-9 * np.abs(H).max())


@needs_node
def test_browser_closures_match_the_language(quantum, tmp_path):
    nodes, records = quantum
    from fieldbridge.quantum import language as ql
    for nid in QUANTUM:
        rec = records[nid]
        got = run({"op": "unitary_analyze", "engine": rec["engine"]}, tmp_path)
        row = ql.derive_bloch_rotation(ql.load(nodes[nid]["spec"]))
        assert got["closure"] == rec["facts"]["closure"], nid
        assert len(got["frequencies"]) == rec["facts"]["frequencies"], nid
        assert np.allclose(got["frequencies"][:6], rec["facts"]["frequency_list"], rtol=1e-5), nid
        assert got["dim"] == row["algebra_dimension"], nid
        assert got["reached"] == (row["status"] == "reached"), nid
        if got["reached"]:
            assert abs(got["rate"] - row["signature"]["rate"]) < 1e-9, nid
            assert abs(got["theta"] - row["signature"]["theta_deg"]) < 1e-7, nid
        if rec["facts"].get("cause"):
            term = rec["engine"]["terms"][got["cause"]]
            assert sd.operator_html(term["operator"]) in rec["facts"]["cause"], nid


@needs_node
def test_browser_evolution_is_exact_and_follows_the_rabi_law(quantum, tmp_path):
    nodes, records = quantum
    from scipy.linalg import expm
    from fieldbridge.quantum import language as ql
    times = [0.0, 0.3, 1.1, 2.7, 5.0]
    for nid in ("two_spins", "two_spins_x1", "nmr_chain", "interacting_bosons", "spin1_easy_axis"):
        rec, real = records[nid], ql.load(nodes[nid]["spec"])
        V, _ = ql.restrict(real)
        H, O = V.conj().T @ real.H @ V, V.conj().T @ real.O @ V
        w, U = np.linalg.eigh(O)
        psi0, o0 = U[:, -1], np.trace(O).real / O.shape[0]
        exact = [((psi0.conj() @ expm(1j * H * t) @ O @ expm(-1j * H * t) @ psi0).real - o0) / (w[-1] - o0)
                 for t in times]
        got = run({"op": "unitary_signal", "engine": rec["engine"], "times": times}, tmp_path)
        assert np.allclose([r["f"] for r in got], exact, atol=1e-10), nid
        if rec["class"] == "rotation":
            th, W = np.deg2rad(rec["facts"]["theta"]), rec["facts"]["rate"]
            assert np.allclose(exact, np.cos(th) ** 2 + np.sin(th) ** 2 * np.cos(W * np.array(times)), atol=1e-9)
            assert np.allclose([np.linalg.norm(r["m"]) for r in got], 1.0, atol=1e-9), nid  # on the sphere
    # the anisotropy term takes the expectation of J off the sphere
    got = run({"op": "unitary_signal", "engine": records["spin1_easy_axis"]["engine"], "times": [0.0, 1.5]}, tmp_path)
    assert abs(np.linalg.norm(got[0]["m"]) - 1) < 1e-9 and np.linalg.norm(got[1]["m"]) < 0.99


@needs_node
def test_a_larger_algebra_changes_the_signal_only_when_the_closure_of_the_observable_grows(quantum, tmp_path):
    """With the field moved to the second spin, H and X0 generate an algebra of dimension 6 and the derivation of a
    rotation stops, but the closure of X0 keeps three operators and the signal is that of the two spins of Chapter
    11. With the field on both spins the closure has five operators and two frequencies: the law is left."""
    nodes, records = quantum
    times = list(np.linspace(0.0, 12.0, 61))
    f = {nid: np.array([r["f"] for r in run({"op": "unitary_signal", "engine": records[nid]["engine"], "times": times},
                                             tmp_path)]) for nid in ("two_spins", "two_spins_x1", "two_spins_both")}
    facts = {nid: records[nid]["facts"] for nid in f}
    assert [records[n]["class"] for n in f] == ["rotation", "obstructed", "obstructed"]
    assert (facts["two_spins_x1"]["dim"], facts["two_spins_x1"]["closure"], facts["two_spins_x1"]["frequencies"]) == (6, 3, 1)
    assert facts["two_spins_x1"]["frequency_list"] == pytest.approx([5 ** 0.5])
    assert np.max(np.abs(f["two_spins_x1"] - f["two_spins"])) < 1e-10
    assert (facts["two_spins_both"]["closure"], facts["two_spins_both"]["frequencies"]) == (5, 2)
    assert np.max(np.abs(f["two_spins_both"] - f["two_spins"])) > 0.5
    # every other realization of the page with a larger algebra has more than one frequency
    for nid, rec in records.items():
        if rec["class"] == "obstructed" and nid != "two_spins_x1":
            assert rec["facts"]["frequencies"] > 1, nid
        if rec["class"] == "rotation":
            assert rec["facts"]["frequencies"] == 1, nid


@needs_node
def test_the_two_preparations_of_chapter_11_give_opposite_signals(quantum, tmp_path):
    rec = quantum[1]["two_spins_h0"]
    plus, minus = (p for p in rec["engine"]["preparations"] if p["id"] in ("plus_y", "minus_y"))
    t = [0.0, 0.2, 0.5]
    fp = [r["f"] for r in run({"op": "unitary_signal", "engine": rec["engine"], "times": t, "psi0": plus["state"]}, tmp_path)]
    fm = [r["f"] for r in run({"op": "unitary_signal", "engine": rec["engine"], "times": t, "psi0": minus["state"]}, tmp_path)]
    assert abs(fp[0]) < 1e-12 and abs(fm[0]) < 1e-12
    assert np.allclose(fp, -np.array(fm), atol=1e-12) and abs(fp[2]) > 0.5
    assert np.allclose(fp, -np.sin(2 * np.array(t)), atol=1e-12)  # x(t) = -c(0) sin 2gt with c(0) = +1


def test_the_classical_limit_of_the_anisotropic_spin_is_the_stoner_wohlfarth_energy():
    """<H>/(2|D| j^2) = e(phi) - 1/2 - cos^2(phi)/(4j) for every spin; the correction vanishes as j grows."""
    nodes = sd.resolve(only=["spin1_easy_axis", "stoner_wohlfarth"])
    q, d = nodes["spin1_easy_axis"]["spec"], nodes["stoner_wohlfarth"]["spec"]
    assert sd.classical_limit(q, d) < 1e-12
    for twoj in range(2, 8):
        j = twoj / 2
        spec = copy.deepcopy(q)
        spec["carrier"] = {"kind": "spin", "j": f"{twoj}/2"}
        # the same reduced field h = |b| / (2|D| j) at every j
        spec["hamiltonian"][0]["coefficient"] = 2.0 * j
        assert sd.classical_limit(spec, d, j) < 1e-12
    with pytest.raises(sd.SiteError):
        wrong = copy.deepcopy(d)
        wrong["parameters"]["psi"] = 0.3
        sd.classical_limit(q, wrong)


# ------------------------------------------------------------------------------------------------ memory
def model(nid):
    node = sd.resolve(only=[nid])[nid]
    rec = model_record(node["spec"], nid, node["path"])
    return node, rec


@needs_node
def test_browser_eigenvalues_of_real_matrices_match_numpy(tmp_path):
    rng = np.random.default_rng(4)
    mats = [rng.normal(size=(n, n)) for n in (3, 4, 5, 8)]
    Q = np.linalg.qr(rng.normal(size=(4, 4)))[0]
    mats.append(Q @ np.diag([-1.0, -1.0 + 1e-6, -1.0 - 1e-6, -1.0 + 4e-6]) @ np.linalg.inv(Q))  # a cluster
    mats.append(np.array([[0.0, -1, 0], [1, 0, 0], [0, 0, -2]]))  # a rotation
    for A in mats:
        got = run({"op": "eigenvalues", "matrix": A.tolist()}, tmp_path)
        want = np.sort_complex(np.linalg.eigvals(A))
        z = np.sort_complex(np.array([v["re"] + 1j * v["im"] for v in got]))
        assert np.allclose(z, want, atol=1e-8 * max(1.0, np.abs(A).max()))


@needs_node
@pytest.mark.parametrize("nid", ["toggle", "toggle_unequal", "stoner_wohlfarth", "sw_oblique", "laser",
                                 "repressor_ring4", "repressilator_activation"])
def test_browser_equilibria_match_fixed_points(nid, tmp_path):
    from fieldbridge.memory import analysis as an
    from fieldbridge.memory import spec as mspec
    node, rec = model(nid)
    real = mspec.load(node["spec"])
    py = an.fixed_points(real, np.random.default_rng(1), 40)
    py_stable = [q for q, k in py if an.n_unstable(k) == 0 and an.classify(k)["flat"] == 0]
    got = run({"op": "equilibria", "model": rec, "n": 40, "seed": 5}, tmp_path)
    js_stable = [e["q"] for e in got if e["stable"]]
    assert len(js_stable) == len(py_stable), nid
    for q in py_stable:
        assert min(np.linalg.norm(real.carrier.diff(np.asarray(p), q)) for p in js_stable) < 1e-6, nid
    for e in got:  # eigenvalues against the Python Jacobian
        k = an.kappa_spectrum(real, np.asarray(e["q"]))
        assert np.allclose(sorted(-k.real), sorted(v["re"] for v in e["eig"]), atol=1e-5), nid


@needs_node
@pytest.mark.parametrize("nid,param,below,above,change", [
    ("toggle", "alpha", 1.9, 2.1, (1, 2)), ("stoner_wohlfarth", "h", 0.95, 1.05, (2, 1)),
    ("sw_oblique", "h", 0.55, 0.60, (2, 1)), ("toggle_unequal", "alpha", 3.2, 3.35, (1, 2))])
def test_browser_scan_changes_the_number_of_states_at_the_write_point(nid, param, below, above, change, tmp_path):
    _, rec = model(nid)
    got = run({"op": "scan", "model": rec, "name": param, "values": [below, above], "n": 24}, tmp_path)
    assert tuple(sum(p["stable"] for p in row["points"]) for row in got) == change


@needs_node
def test_a_drive_locks_the_oscillator_inside_the_adler_range_and_slips_outside(tmp_path):
    data = sd.build(only=["van_der_pol"], strict=False, log=lambda *a: None)
    rec = data["nodes"]["van_der_pol"]
    drive = rec["engine"]["drive"]
    n, w0 = drive["ratio"], drive["omega"]
    K, eps = 10 * drive["K"], 10 * drive["eps"]  # ten times the derivation's weak drive, to lock within ~20 periods
    q0 = [2.0, 0.0]
    phases = {}
    for nu in (0.4, 3.0):
        wf = n * (w0 - nu * K)
        got = run({"op": "driven_phase", "model": rec["engine"], "drive": {"param": drive["param"], "p1": drive["p1"],
                                                                          "eps": eps, "omega": wf},
                   "q0": q0, "T": 60 * 2 * np.pi / w0, "ratio": n}, tmp_path)
        phases[nu] = np.unwrap(2 * np.pi * np.array([p[1] for p in got]))
    assert np.ptp(phases[0.4][-15:]) < 0.02          # locked: the phase against the drive is constant
    assert abs(phases[3.0][-1] - phases[3.0][-15]) > 0.5   # outside the range the phase slips


@needs_node
@pytest.mark.parametrize("nid", ["colloid_patch", "toggle", "stoner_wohlfarth"])
def test_every_kind_of_dynamics_runs_in_the_browser_as_exported_for_the_page(nid, tmp_path):
    """The record of the page keeps the kind of dynamics (equations, gene, rotor) that the browser integrator reads."""
    node = sd.resolve(only=[nid])[nid]
    engine = sd.dissipative_record(node, {})["engine"]
    got = run({"op": "attractors", "model": engine, "opts": {"n": 6, "T": 30, "seed": 3}}, tmp_path)
    assert got["unbounded"] == 0 and len(got["points"]) >= 1, nid


# ------------------------------------------------------------------------------------------------ fields
@needs_node
def test_browser_field_loss_is_the_exact_spectral_sum(tmp_path):
    from fieldbridge.memory import fields as mf
    times = [0.5, 3.0, 20.0, 150.0]
    for name in ("nonconserved_1d", "conserved_1d_charge", "conserved_1d_dipole", "conserved_2d_charge"):
        spec = json.loads((ROOT / "examples/memory/fields" / f"{name}.json").read_text(encoding="utf-8"))
        model_ = mf.from_spec(spec)
        engine = {**spec["field"], "kappa0": model_.kappa0, "Dgrad": model_.Dgrad}
        got = run({"op": "fields_snr", "engine": engine, "times": times}, tmp_path)
        assert np.allclose(got, mf.spectral_snr(model_, times)["snr_profile"], rtol=1e-9), name
        prof = run({"op": "fields_profile", "engine": engine, "t": 10.0}, tmp_path)
        if spec["field"]["conserved"] and spec["field"]["d"] == 1:
            assert abs(sum(prof) - (spec["field"]["amplitude"] if spec["field"]["write"] == "charge" else 0)) < 1e-9


@needs_node
def test_equations_are_written_as_mathml(tmp_path):
    from fieldbridge.web_models import arithmetic
    tree = run({"op": "mathml_tree", "tree": arithmetic("alpha/(1 + v**n) - u")}, tmp_path)
    assert "<mfrac>" in tree and "<msup><mi>v</mi><mi>n</mi></msup>" in tree and "α" in tree
    term = run({"op": "mathml_term", "tree": arithmetic("U/2"), "operator": "na na + nb nb", "hc": False}, tmp_path)
    assert "<msup><msub><mi>n</mi><mi>a</mi></msub><mn>2</mn></msup>" in term


@needs_node
def test_a_constant_is_printed_to_the_decimal_place_of_its_uncertainty(tmp_path):
    """The page printed the delay constant as 1.0188 +- 5.0e-7: a value rounded 7e-6 away from the one measured."""
    pairs = [[1.0187931203, 5e-7], [1.0187733792, 2.87e-4], [1.3271, 0.014], [1.0000066, 0.0025286], [-0.5, 0.02],
             [2.0, 0.0]]
    shown = run({"op": "value_with_uncertainty", "pairs": pairs}, tmp_path)
    assert shown == ["1.01879312 ± 0.00000050", "1.01877 ± 0.00029", "1.327 ± 0.014", "1.0000 ± 0.0025",
                     "−0.500 ± 0.020", "2"]
