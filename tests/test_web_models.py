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
