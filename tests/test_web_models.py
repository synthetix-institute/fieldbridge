"""Browser dynamics checked against the tutorial's Python realizations."""
import json
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pytest
from scipy.integrate import solve_ivp

from fieldbridge.memory.spec import load
from fieldbridge.web_models import arithmetic, export_models

ROOT = Path(__file__).resolve().parents[1]
NODE = shutil.which("node")
MODELS = export_models(ROOT, [])


def browser(model, action, **kwargs):
    program = "const M=require(process.argv[1]);const x=JSON.parse(process.argv[2]);process.stdout.write(JSON.stringify(M[x.action](x.model,...x.args)));"
    args = [kwargs["q"], kwargs.get("params", model["params"]), kwargs.get("removed", [])] if action == "drift" else [kwargs]
    return json.loads(subprocess.check_output([NODE, "-e", program, str(ROOT / "fieldbridge/web/model-physics.js"), json.dumps({"model": model, "action": action, "args": args})], text=True))


@pytest.mark.parametrize("model", MODELS, ids=lambda m: m["id"])
def test_each_exported_model_matches_its_original_python_drift(model):
    real = load(ROOT / model["specification"])
    rng = np.random.default_rng(137)
    for multiplier in (0.8, 1.0, 1.2):
        q = rng.uniform(.15, 1.2, real.carrier.dim)
        params = {k: v * multiplier for k, v in model["params"].items()}
        expected = real.F(q, **params)
        np.testing.assert_allclose(browser(model, "drift", q=q.tolist(), params=params), expected, rtol=1e-11, atol=1e-11)


@pytest.mark.parametrize("name", ["toggle", "schlogl", "tubes", "repressilator", "colloid_patch", "van_der_pol"])
def test_browser_trajectory_matches_independent_scipy_solver(name):
    model = next(m for m in MODELS if m["id"] == name)
    real = load(ROOT / model["specification"])
    expected = solve_ivp(lambda _, q: real.F(q), [0, 3], model["initial"], rtol=1e-10, atol=1e-11).y[:, -1]
    result = browser(model, "integrate", duration=3)
    assert result["status"] == "complete"
    np.testing.assert_allclose(result["final"], expected, atol=2e-6, rtol=2e-6)


def test_detaching_repression_changes_the_actual_equations():
    model = next(m for m in MODELS if m["id"] == "toggle")
    intact = browser(model, "drift", q=[1, 2])
    detached = browser(model, "drift", q=[1, 2], removed=["u-0"])
    assert intact != detached
    assert detached[0] == -1
    assert detached[1] == intact[1]


def test_source_expression_export_refuses_executable_syntax():
    with pytest.raises(ValueError):
        arithmetic("__import__('os').system('touch unsafe')")


def test_scalar_dependency_graph_retains_self_couplings():
    model = next(m for m in MODELS if m["id"] == "schlogl")
    assert {link["term"] for link in model["links"]} == {"x-0", "x-1", "x-2"}
    assert all(link["source"] == link["target"] == 0 for link in model["links"])


@pytest.mark.parametrize("name", ["schlogl", "toggle_unequal", "tubes_unequal"])
def test_design_is_the_tutorials_symmetry_restoration(name):
    model = next(m for m in MODELS if m["id"] == name)
    result = browser(model, "design", **model["params"])
    if name == "schlogl":
        x = result["a"] / 3
        assert abs(-x**3 + result["a"]*x**2 - result["k3"]*x + result["b"]) < 1e-12
        assert abs(-3*x**2 + 2*result["a"]*x - result["k3"]) < 1e-12
        assert abs(-6*x + 2*result["a"]) < 1e-12
    elif name == "toggle_unequal":
        assert result["gamma"] == 1
    else:
        assert result["L2"] == result["L1"]


def handlers(tmp_path, actions):
    path = tmp_path / "models.json"
    path.write_text(json.dumps(MODELS))
    return json.loads(subprocess.check_output([NODE, str(ROOT / "tests/web_studio_harness.cjs"), json.dumps(actions), "#gallery", str(path)], text=True))


def test_material_choice_loads_native_equations_and_calculates(tmp_path):
    rows = handlers(tmp_path, [{"id": "model-select", "value": "schlogl", "type": "change"}])
    assert "genetic toggle" in rows[0]["modelName"]
    assert "Schl" in rows[1]["modelName"]
    assert rows[0]["modelEquations"] != rows[1]["modelEquations"]
    assert rows[0]["modelFinal"] != rows[1]["modelFinal"]


def test_real_model_term_detaches_and_changes_prediction(tmp_path):
    rows = handlers(tmp_path, [{"selector": "[data-model-term=u-0]", "checked": False, "type": "change"}])
    assert 'removed' in rows[1]["modelEquations"]
    assert rows[0]["modelFinal"] != rows[1]["modelFinal"]


def test_parameter_recalculates_real_model_not_a_saved_card(tmp_path):
    rows = handlers(tmp_path, [{"id": "model-coefficient-alpha", "value": 1, "type": "input"}])
    assert rows[0]["modelFinal"] != rows[1]["modelFinal"]
    assert rows[0]["modelMetrics"] != rows[1]["modelMetrics"]


def test_both_preparations_are_used_by_the_calculation(tmp_path):
    rows = handlers(tmp_path, [{"id": "model-alternate-0", "value": 3, "type": "input"}])
    assert rows[0]["modelFinal"] == rows[1]["modelFinal"]
    assert rows[0]["modelMetrics"] != rows[1]["modelMetrics"]


def test_detaching_repressor_production_updates_equation_and_dynamics(tmp_path):
    rows = handlers(tmp_path, [{"id": "model-select", "value": "repressilator", "type": "change"},
                               {"selector": "[data-model-term=u0-0]", "checked": False, "type": "change"}])
    assert rows[-2]["modelEquations"] != rows[-1]["modelEquations"]
    assert rows[-2]["modelFinal"] != rows[-1]["modelFinal"]


def test_worker_runs_the_same_three_trajectories(tmp_path):
    model = next(m for m in MODELS if m["id"] == "toggle")
    program = """
const fs=require('fs'),vm=require('vm'),path=require('path');
const state={params:input.model.params,initial:input.model.initial,
 alternate:input.model.alternate,duration:3,removed:['u-0'],compare:true};
const context={URL,Math,Set,console};context.globalThis=context;
context.self={location:{href:'https://example.org/studio/model-worker.js'},
 postMessage:result=>process.stdout.write(JSON.stringify(result))};
context.importScripts=()=>{vm.runInContext(fs.readFileSync(path.join(input.assets,'model-physics.js'),'utf8'),context);context.self.FieldBridgeModels=context.FieldBridgeModels;};
vm.createContext(context);
vm.runInContext(fs.readFileSync(path.join(input.assets,'model-worker.js'),'utf8'),context);
context.self.onmessage({data:{id:17,model:input.model,state}});
"""
    output = subprocess.check_output([NODE, "-e", "const input=" + json.dumps({"model": model, "assets": str(ROOT / "fieldbridge/web")}) + ";" + program], text=True)
    result = json.loads(output)
    assert result["id"] == 17
    assert len(result["results"]) == 3
    assert all(row["status"] == "complete" for row in result["results"])
    expected = browser(model, "integrate", duration=3, removed=["u-0"])
    np.testing.assert_allclose(result["results"][0]["final"], expected["final"])


def test_symmetry_design_changes_parameters_and_calculation(tmp_path):
    rows = handlers(tmp_path, [{"id": "model-select", "value": "schlogl", "type": "change"}, {"id": "model-solve", "type": "click"}])
    assert "Triple root" in rows[-1]["modelDesign"]
    assert "− -" not in rows[-1]["modelEquations"]
    assert rows[-1]["modelFinal"] != rows[-2]["modelFinal"]
