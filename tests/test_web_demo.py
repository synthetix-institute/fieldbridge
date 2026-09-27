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
