"""Compare the browser's calculations with independent Python consequences."""
import json
import math
from pathlib import Path
import shutil
import subprocess

import pytest

ROOT = Path(__file__).resolve().parents[1]
ENGINE = ROOT / "fieldbridge/web/physics.js"
NODE = shutil.which("node")

def browser(expression):
    if not NODE:
        pytest.skip("Node is needed for browser-equation comparisons")
    source = f"const P=require({json.dumps(str(ENGINE))});process.stdout.write(JSON.stringify({expression}));"
    return json.loads(subprocess.check_output([NODE, "-e", source], text=True))

def memory(**changes):
    return {"feedback": True, "saturation": True, "field": True, "thermal": False,
            "eps": 1, "gamma": 1, "h": 0, "D": .04, **changes}

def test_offline_build_embeds_verified_reports(tmp_path):
    pytest.importorskip("sympy")
    from fieldbridge.web_demo import build_studio
    page = build_studio(tmp_path).read_text()
    reports = json.loads((tmp_path / "verified_examples.json").read_text())
    assert len(reports) == 3
    assert all(r["calculation"]["residual_coefficients"] == ["0", "0"] for r in reports.values())
    assert reports["ito_square"]["calculation"]["candidate"]["passes_local_identity"] is False
    assert "__CONFIG__" not in page
    assert "https://github.com/synthetix-institute/fieldbridge/blob/main/docs/tutorial/" in page
    assert "https://doi.org/10.3792/pia/1195572786" in page
    assert (tmp_path / "studio/physics.js").exists()
    assert (tmp_path / "studio/lucide.min.js").exists()
    assert "Equations Derived & Verified!" not in page
    assert 'studio/collections.css?v=' in page
    assert 'studio/collections.js?v=' in page
    assert len(list((tmp_path / "gallery").glob("card*.json"))) == 12
    assert len(list((tmp_path / "gallery").glob("card*.png"))) == 12
    assert "Memory in model materials" in page
    assert "two coupled spins" in page
    quantum = json.loads((tmp_path / "quantum_examples.json").read_text())
    assert len(quantum["attachments"]) == 7
    assert all(r["attached"] and r["law"]["target_residual"] < 1e-8
               for r in quantum["attachments"].values())
    assert "index.html#gallery" in (tmp_path / "gallery.html").read_text()

def test_memory_states_barrier_and_fold():
    s=memory()
    result=browser(f"P.memory({json.dumps(s)})")
    assert [e["x"] for e in result["stable"]] == pytest.approx([-1,1])
    assert result["barrier"] == pytest.approx(.25)
    assert result["threshold"] == pytest.approx(2/(3*math.sqrt(3)))
    for e in result["equilibria"]:
        assert e["x"]-e["x"]**3 == pytest.approx(0, abs=1e-12)

@pytest.mark.parametrize("part", ["feedback", "saturation"])
def test_omission_breaks_bistability(part):
    s=memory();s[part]=False
    assert browser(f"P.memory({json.dumps(s)})")["bistable"] is False

def test_nonlinear_attraction_at_critical_point():
    s=memory();s["eps"]=0
    result=browser(f"P.memory({json.dumps(s)})")
    assert len(result["stable"]) == 1
    assert result["stable"][0]["x"] == pytest.approx(0)

def test_writing_pulse_against_independent_ode_integrator():
    scipy=pytest.importorskip("scipy.integrate")
    # Integrate three intervals separately so no solver interpolates across a discontinuity.
    q=-1.0
    for start,end,field in [(0,1,0),(1,3,.8),(3,8,0)]:
        result=scipy.solve_ivp(lambda t,y:y-y**3+field,(start,end),[q],rtol=1e-10,atol=1e-12)
        q=float(result.y[0,-1])
    output=browser(f"P.trajectory(P.memory({json.dumps(memory())}),-1,.8)")
    assert not output["diverged"]
    assert output["data"][-1][1] > .99
    assert output["data"][-1][1] == pytest.approx(q, abs=2e-5)
    unforced=browser(f"P.trajectory(P.memory({json.dumps(memory())}),-1,0)")
    assert unforced["data"][-1][1] == pytest.approx(-1)

@pytest.mark.parametrize("state_map", ["square", "log"])
@pytest.mark.parametrize("noise", [True, False])
@pytest.mark.parametrize("correction", [True, False])
@pytest.mark.parametrize("convention", ["ito", "stratonovich"])
def test_exported_model_and_browser_agree_with_symbolic_verifier(state_map,noise,correction,convention):
    sp=pytest.importorskip("sympy")
    from fieldbridge.verification import verify_construction
    s=dict(map=state_map,noise=noise,correction=correction,convention=convention,
           theta=1.2,alpha=.7,sigma=.8,mu=.2)
    js=browser(f"({{spec:P.constructionSpec({json.dumps(s)}),result:P.stochastic({json.dumps(s)})}})")
    report=verify_construction(js["spec"])
    values=js["spec"]["numeric_parameters"]
    substitute=lambda text:float(sp.sympify(text).subs({**values,"y":0}))
    assert js["result"]["exact"] == pytest.approx(substitute(report["target"]["drift"]))
    assert js["result"]["residual"] == pytest.approx(substitute(report["candidate"]["residual_d_phi"]))
    assert report["candidate"]["passes_local_identity"] == (abs(js["result"]["residual"]) < 1e-10)

def test_convention_changes_mean_but_not_variance():
    s=dict(map="log",noise=True,correction=True,theta=1,alpha=1,sigma=.8,mu=.2,convention="ito")
    ito=browser(f"P.stochastic({json.dumps(s)})")
    s["convention"]="stratonovich"
    strat=browser(f"P.stochastic({json.dumps(s)})")
    assert strat["exact"]-ito["exact"] == pytest.approx(.32)
    assert strat["variance"] == ito["variance"] == pytest.approx(.64)

@pytest.mark.parametrize("carrier", ["correlated-pair", "qubit", "spin", "bosons", "chain", "fermion-pair", "collective"])
@pytest.mark.parametrize("g,h", [(1,.5),(.7,.2),(1.5,0),(0,.5),(0,0)])
def test_browser_spin_construction_against_exact_hamiltonian(carrier,g,h):
    np=pytest.importorskip("numpy")
    from fieldbridge.quantum.language import load, restrict
    s=dict(coupling=True,field=True,g=g,h=h)
    data=browser(f"({{spec:P.spinSpec({json.dumps(s)},{json.dumps(carrier)}),values:[0,.2,.7,1.3,2].map(t=>P.spin({json.dumps(s)}).signal(t))}})")
    real=load(data["spec"])
    V,_=restrict(real)
    H=V.conj().T@real.H@V;O=V.conj().T@real.O@V
    baseline=np.trace(O).real/O.shape[0]
    _,states=np.linalg.eigh(O);psi0=states[:,-1]
    normalizer=(psi0.conj()@O@psi0).real-baseline
    E,U=np.linalg.eigh(H);coefficients=U.conj().T@psi0
    exact=[]
    for t in [0,.2,.7,1.3,2]:
        psi=U@(np.exp(-1j*E*t)*coefficients)
        exact.append(float(((psi.conj()@O@psi).real-baseline)/normalizer))
    assert data["values"] == pytest.approx(exact, abs=1e-10)

def test_removing_spin_coupling_conserves_observable():
    s=dict(coupling=False,field=True,g=1,h=.5)
    assert browser(f"[0,1,2,3].map(t=>P.spin({json.dumps(s)}).signal(t))") == [1,1,1,1]

def studio(actions, initial_hash=""):
    if not NODE:
        pytest.skip("Node is needed for studio event-handler tests")
    return json.loads(subprocess.check_output([
        NODE, str(ROOT / "tests/web_studio_harness.cjs"), json.dumps(actions), initial_hash
    ], text=True))

def test_memory_slider_recalculates_equation_without_build_action():
    initial, changed = studio([{"id":"param-eps", "value":"-0.8", "type":"input"}])
    assert "-0.8x" in changed["equation"]
    assert changed["metrics"] != initial["metrics"]
    assert "pending" not in changed["status"]

def test_memory_detach_recalculates_equation_and_prediction():
    initial, changed = studio([{"id":"part-feedback", "checked":False, "type":"change"}])
    assert "1x " not in changed["equation"]
    assert changed["metrics"] != initial["metrics"]
    assert "One attracting state" in changed["consequence"]

def test_graph_part_buttons_change_the_actual_model():
    initial, changed = studio([{"selector":"[data-graph-part=feedback]", "type":"click"}])
    assert changed["equation"] != initial["equation"]
    assert "One attracting state" in changed["consequence"]

def test_pulse_parameter_controls_switching():
    strong, weak = studio([{"id":"param-pulse", "value":"0.1", "type":"input"}])
    assert "switches the negative" in strong["consequence"]
    assert "remain in opposite" in weak["consequence"]
    assert "0.1" in weak["drive"]
    assert weak["metrics"] != strong["metrics"]

def test_coordinate_change_updates_equation_and_final_states():
    initial, changed = studio([{"id":"realization", "value":"displacement", "type":"change"}])
    assert "dq/dt" in changed["equation"]
    assert "0.25q" in changed["equation"]
    assert "Final state from −2" in changed["metrics"]

def test_plot_modes_do_not_require_a_build_after_parameter_change():
    initial, changed, potential = studio([
        {"id":"param-eps", "value":"-0.8", "type":"input"},
        {"selector":"[data-plot=potential]", "type":"click"},
    ])
    assert "Barrier at h = 0" in potential["metrics"]
    assert potential["equation"] == changed["equation"]

def test_spin_slider_updates_hamiltonian_without_calculate_action():
    initial, changed = studio([{"id":"spin-g", "value":"1.5", "type":"input"}])
    assert "1.5 Z0 Z1" in changed["spinEquation"]
    assert changed["spinMetrics"] != initial["spinMetrics"]

def test_legacy_application_route_opens_the_constructor():
    result=studio([], "#applications")[0]
    assert result["hash"] == "#memory"
    assert "dx/dt" in result["equation"]

def test_attached_spin_carrier_keeps_live_parameter_updates():
    snapshots=studio([
        {"id":"spin-detach", "type":"click"},
        {"id":"spin-carrier", "value":"qubit", "type":"change"},
        {"id":"spin-attach", "type":"click"},
        {"id":"spin-g", "value":"1.5", "type":"input"},
    ])
    assert "1.5 X0" in snapshots[-1]["spinEquation"]
    assert "One spin-½" in snapshots[-1]["spinMetrics"]

@pytest.mark.parametrize("pulse_duration", [.2,.73,2,4.1])
def test_editable_pulse_matches_independent_piecewise_integration(pulse_duration):
    scipy=pytest.importorskip("scipy.integrate")
    x=-1.0
    for start,end,field in [(0,1,0),(1,1+pulse_duration,.8),(1+pulse_duration,8,0)]:
        solution=scipy.solve_ivp(lambda t,y:y-y**3+field,(start,end),[x],rtol=1e-10,atol=1e-12)
        x=float(solution.y[0,-1])
    js=browser(f"P.trajectory(P.memory({json.dumps(memory())}),-1,.8,8,800,{{start:1,end:{1+pulse_duration}}})")
    assert js["data"][-1][1] == pytest.approx(x, abs=2e-7)
