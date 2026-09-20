import json
import math
from pathlib import Path
import random
import shutil
import statistics

import pytest

sp = pytest.importorskip("sympy")
from fieldbridge.calculation_adapter import canonical_source_equation
from fieldbridge.constructor import construct_transfer
from fieldbridge.verification import ConstructionError, expression, verify_construction

ROOT = Path(__file__).resolve().parents[1]


def spec(convention="ito"):
    return json.loads((ROOT / f"examples/construction/log_signal_{convention}.json").read_text())


def symbolic(value):
    return expression(value, {"mu": sp.Symbol("mu", real=True),
                              "sigma": sp.Symbol("sigma", positive=True)})


def test_same_display_different_convention_changes_only_the_log_drift():
    ito = verify_construction(spec())
    strat = verify_construction(spec("stratonovich"))
    assert sp.simplify(symbolic(strat["target"]["drift"]) - symbolic(ito["target"]["drift"])
                       - symbolic("sigma**2/2")) == 0
    assert ito["target"]["variance"] == strat["target"]["variance"] == "sigma**2"
    assert ito["target_convention"] == strat["target_convention"] == "ito"
    assert ito["domains"] == {"source": "positive", "target": "real"}


@pytest.mark.parametrize("convention", [None, "", "unspecified", "Ito", [], {}])
def test_unknown_or_missing_convention_never_defaults(convention):
    source = spec()
    if convention is None:
        source.pop("convention")
    else:
        source["convention"] = convention
    with pytest.raises(ConstructionError, match="convention"):
        verify_construction(source)


def test_convention_conversion_preserves_the_same_process():
    source = spec("stratonovich")
    source["drift"] = "(mu-sigma**2/2)*x"
    source["assumptions"] = ["The source is the Stratonovich representation of the Ito benchmark."]
    assert verify_construction(source)["target"] == verify_construction(spec())["target"]


def test_consistent_coefficient_renaming_is_an_invariance_control():
    source = spec()
    source.update(parameters={"growth": "real", "amplitude": "positive"},
                  drift="growth*x", noise="amplitude*x")
    renamed = verify_construction(source)["target"]
    names = {"growth": symbolic("mu"), "amplitude": symbolic("sigma")}
    target = verify_construction(spec())["target"]
    for coefficient in ("drift", "variance"):
        assert sp.simplify(expression(renamed[coefficient], names) - symbolic(target[coefficient])) == 0


def test_logarithm_requires_distinct_source_and_target_domains():
    for key, value in (("domain", "real"), ("target_domain", "positive")):
        source = spec()
        source[key] = value
        with pytest.raises(ConstructionError, match="domain"):
            verify_construction(source)
    source = spec()
    source.pop("target_domain")
    with pytest.raises(ConstructionError, match="target_domain"):
        verify_construction(source)


def _retrieved_log(tmp_path, convention):
    data = tmp_path / str(convention)
    shutil.copytree(ROOT / "examples/calculated_transfer/data", data)
    path = data / "index/core_examples.json"
    records = json.loads(path.read_text())
    source = records[0]["calculation_source"]
    source.update(family="scalar_sde", parameters={"mu": "real", "sigma": "positive"},
                  domain="positive", drift="mu*x", noise="sigma*x",
                  assumptions=[f"The source convention is {convention}.", "X is positive."])
    if convention is None:
        source.pop("convention", None)
        source["assumptions"] = ["X is positive; convention unstated."]
    else:
        source["convention"] = convention
    records[0]["equations"][0] = canonical_source_equation(source)
    records[0]["title"] = "Multiplicative Brownian fluctuations"
    records[0]["summary"] = "A positive stochastic coordinate driven by multiplicative Brownian noise."
    records[0]["keywords"] = ["stochastic", "Brownian", "multiplicative", "diffusion"]
    path.write_text(json.dumps(records))
    request = json.loads((ROOT / "examples/calculated_transfer/squared_signal.json").read_text())
    request.update(state_map="log(x)", inverse_map="exp(y)", target_domain="real",
                   question="What mean growth rate and variance follow for log X?",
                   assumptions=["Y=log(X) on the positive source domain."])
    request.pop("candidate")
    result = construct_transfer("Brownian stochastic diffusion", "stochastic_dynamics",
                                data_dir=data, include_hyperion=False, calculation_request=request)
    return result, records[0]["equations"][0]


def test_convention_is_used_along_the_retrieval_to_prediction_path(tmp_path):
    ito, equation = _retrieved_log(tmp_path, "ito")
    strat, other_equation = _retrieved_log(tmp_path, "stratonovich")
    assert equation == other_equation
    assert ito.calculation["status"] == strat.calculation["status"] == "calculated"
    assert ito.calculation["binding"]["source_model_sha256"] != strat.calculation["binding"]["source_model_sha256"]
    for result, convention in ((ito, "ito"), (strat, "stratonovich")):
        assert result.calculation["report"]["target"] == verify_construction(spec(convention))["target"]
    assert ito.predictions != strat.predictions
    missing, _ = _retrieved_log(tmp_path, None)
    assert missing.calculation["status"] == "refused"
    assert not missing.predictions


def _simulate_source(mu, sigma, steps=512, paths=4000):
    """Integrate X directly: Euler-Maruyama and stochastic Heun share increments."""
    rng = random.Random(20260920)
    dt = 1 / steps
    sqrt_dt = math.sqrt(dt)
    ito, strat = [], []
    for _ in range(paths):
        x, z = 1.0, 1.0
        for _ in range(steps):
            dw = rng.gauss(0, sqrt_dt)
            x += mu*x*dt + sigma*x*dw
            predictor = z + mu*z*dt + sigma*z*dw
            z += mu*(z+predictor)*dt/2 + sigma*(z+predictor)*dw/2
            assert x > 0 and z > 0
        ito.append(math.log(x))
        strat.append(math.log(z))
    return ito, strat


def test_derived_log_moments_against_independent_source_integrators():
    mu, sigma = 0.2, 0.8
    ito, strat = _simulate_source(mu, sigma)
    targets = []
    for convention, values in (("ito", ito), ("stratonovich", strat)):
        target = verify_construction(spec(convention))["target"]
        numeric = lambda value: float(symbolic(value).subs({symbolic("mu"): mu, symbolic("sigma"): sigma}))
        mean, variance = numeric(target["drift"]), numeric(target["variance"])
        targets.append(mean)
        assert abs(statistics.mean(values)-mean) < 4*math.sqrt(variance/len(values))
        assert abs(statistics.variance(values)-variance) < 4*variance*math.sqrt(2/(len(values)-1))
    assert abs(statistics.mean(b-a for a, b in zip(ito, strat)) - (targets[1]-targets[0])) < 0.005
    assert abs(statistics.variance(strat)-statistics.variance(ito)) < 0.005
