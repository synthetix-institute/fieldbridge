import copy
from dataclasses import asdict
import hashlib
import importlib.util
import json
from pathlib import Path

import numpy as np
import pytest
import sympy as sp

from fieldbridge.memory import stochastic_passage as st

loader=importlib.util.spec_from_file_location("passage_design",Path(__file__).resolve().parents[1]/"scripts/run_stochastic_passage_loop.py")
design=importlib.util.module_from_spec(loader)
loader.loader.exec_module(design)

loader=importlib.util.spec_from_file_location("passage_refinement",Path(__file__).resolve().parents[1]/"scripts/refine_stochastic_protocol.py")
refinement=importlib.util.module_from_spec(loader)
loader.loader.exec_module(refinement)

loader=importlib.util.spec_from_file_location("passage_reuse",Path(__file__).resolve().parents[1]/"scripts/learn_verified_protocol.py")
reuse=importlib.util.module_from_spec(loader)
loader.loader.exec_module(reuse)


def task():
    return copy.deepcopy(design.design()[0]["task"])


def test_generator_repair_checks_quadratic_as_well_as_linear_observables():
    c=st.compile_task(task())
    w=st.generator_witness(c)
    assert w["exact_generator_identity"]
    assert len(w["probes"])==9
    assert any(r["linear_defect"]!="0" for r in w["probes"])
    assert all(r["omit_diffusion_defect"]=="0" for r in w["probes"][:3])
    assert any(r["omit_diffusion_defect"]=="2*D" for r in w["probes"][3:])


def test_scaled_generator_is_derived_exactly_on_both_model_families():
    for family in ("coupled_amplitudes","cubic_anisotropy"):
        t=task();t["model"]=design.model(family,t["diffusion"])
        witness=st.scaling_witness(st.compile_task(t))
        assert witness["homogeneous_cubic_scaling_verified"]
        assert witness["residuals"]==["0","0","0"]


@pytest.mark.parametrize("edit",["convention","noise","carrier","moving_branch"])
def test_unsupported_or_ambiguous_physics_is_rejected(edit):
    t=task()
    if edit=="convention":
        t.pop("convention")
    elif edit=="noise":
        t["model"]["noise"]*=2
    elif edit=="carrier":
        t["model"]["carrier"]["kind"]="torus"
    else:
        t["model"]["drift"]["x"]+="+1"
    with pytest.raises(ValueError):
        st.compile_task(t)


def test_final_states_are_actual_stable_equilibria_not_direction_bins():
    for family,count in (("coupled_amplitudes",8),("cubic_anisotropy",4)):
        t=task();t["model"]=design.model(family,t["diffusion"])
        c=st.compile_task(t)
        states=st.final_states(c)
        assert len(states["states"])==count
        assert np.max(np.abs(c["drift"](np.asarray(states["states"]),c["end"])))<1e-8
        assert all(max(e)<0 for e in states["eigenvalues"])
        assert st.labels_of(np.array([[100.,100.,100.]]),states)[0]==-1


def test_rates_do_not_change_the_final_landscape_or_holding_conditions():
    jobs=design.design()
    for family in {j["task"]["model"]["name"] for j in jobs}:
        selected=[j["task"] for j in jobs if j["task"]["model"]["name"]==family]
        ends=[st.compile_task(t)["F"].subs(st.compile_task(t)["control_symbol"],.5) for t in selected]
        assert all(F==ends[0] for F in ends)
        assert len({(t["hold_diffusion"],t["hold_time"],t["relax_time"]) for t in selected})==1


def test_gaussian_prefix_matching_is_computed_from_the_source():
    c=st.compile_task(task())
    states=st.final_states(c)
    p=st.linear_prediction(c,states,st.Policy())
    assert 0<p["matching_time"]<p["t_sector"]
    assert p["matching_ratio"]>=.05
    assert sum(p["sector_probabilities"])==pytest.approx(1)
    assert np.linalg.eigvalsh(p["matching_covariance"]).min()>0


def test_torch_field_matches_source_and_changed_coefficients():
    pytest.importorskip("torch")
    t=task();t["model"]["drift"]["x"]="mu*x-7*x**3+kap*(y+z-2*x)"
    c=st.compile_task(t)
    f=st.torch_drift(c,"cpu")
    import torch
    q=np.array([[.2,.3,-.1]])
    out=f(torch.tensor(q),torch.tensor(.5,dtype=torch.float64)).numpy()
    assert np.allclose(out,c["drift"](q,.5),atol=1e-12)


def test_a_prediction_is_committed_before_reference_calculation(tmp_path,monkeypatch):
    calls=[]
    def simulate(c,states,policy,seed,**kw):
        calls.append("prediction" if kw.get("prediction") is not None else "reference")
        n=len(states["states"])
        return {"probabilities":[1/n]*n,"valid":True,"n":4096}
    monkeypatch.setattr(st,"simulate",simulate)
    out=tmp_path/"case"
    st.prepare_prediction(task(),st.Policy(),3,out)
    assert calls==["prediction","prediction"]
    assert (out/"prediction.json").exists() and not (out/"comparison.json").exists()
    report=st.verify_prediction(out,5)
    assert calls[-1]=="reference" and report["repair_verified"]
    with pytest.raises(ValueError,match="already"):
        st.verify_prediction(out,7)


def test_an_invalid_reference_cannot_accept_a_repair():
    result=st.compare([.5,.5],{"probabilities":[.5,.5],"n":4096,"valid":False},4096,st.Policy())
    assert not result["passed"]


def test_refinement_preserves_the_derived_scaled_noise_and_signal():
    old=task()
    new=refinement.rescale(old,old["protocol"]["rate"]/2)
    r0,r1=old["protocol"]["rate"],new["protocol"]["rate"]
    assert new["diffusion"]/r1==pytest.approx(old["diffusion"]/r0)
    assert np.allclose(np.asarray(new["bias"])/r1**.75,np.asarray(old["bias"])/r0**.75)
    assert new["model"]["drift"]==old["model"]["drift"]
    assert new["model"]["parameters"]==old["model"]["parameters"]
    assert new["model"]["noise"]==new["diffusion"]
    assert all(new[k]==old[k] for k in ("hold_diffusion","hold_time","relax_time"))


def test_rate_only_recipe_distinguishes_the_factorial_controls():
    old=task()
    new=copy.deepcopy(old)
    new["protocol"]["rate"]/=2
    assert reuse.changed_settings(old,new)==["rate"]
    scaled=refinement.rescale(old,old["protocol"]["rate"]/2)
    assert reuse.changed_settings(old,scaled)==["rate","diffusion","bias"]


@pytest.mark.parametrize("field",["drift","hold","end"])
def test_a_protocol_recipe_rejects_a_changed_final_system(field):
    old=task();new=copy.deepcopy(old)
    if field=="drift":
        new["model"]["drift"]["x"]+="-x**3"
    elif field=="hold":
        new["hold_diffusion"]*=2
    else:
        new["protocol"]["end"]+=.1
    with pytest.raises(ValueError):
        reuse.changed_settings(old,new)


def test_the_goal_uses_a_lower_bound_not_the_point_estimate():
    assert refinement.lower_bound(.601,4096)<.60
    assert refinement.lower_bound(.667,4096)>.60


def test_an_incompatible_frozen_policy_cannot_be_used(tmp_path):
    rule={"schema":st.RULE,"policy":asdict(st.Policy(matching_fraction=.025)),"repair":"wrong"}
    p=tmp_path/"rule.json";p.write_text(json.dumps(rule))
    with pytest.raises(ValueError,match="authorize"):
        st.prepare_prediction(task(),st.Policy(),1,tmp_path/"case",frozen_rule=p)


def test_public_command_keeps_prediction_and_reference_separate(tmp_path,monkeypatch):
    from fieldbridge.cli import main
    calls=[]
    def fake_prepare(t,policy,seed,out,**kwargs):
        calls.append("predict")
        Path(out).mkdir()
        (Path(out)/"prediction.json").write_text(json.dumps({"repaired_prediction":{"seed":seed}}))
        return {"prediction_valid":True,"prediction_sha256":"test"}
    def fake_verify(out,seed,**kwargs):
        calls.append(("verify",seed))
        return {"repair_verified":True}
    monkeypatch.setattr(st,"prepare_prediction",fake_prepare)
    monkeypatch.setattr(st,"verify_prediction",fake_verify)
    p=tmp_path/"task.json";p.write_text(json.dumps(task()))
    out=tmp_path/"out"
    assert main(["memory","passage","--task",str(p),"--out-dir",str(out),"--seed","123"])==0
    assert calls==["predict"]
    assert main(["memory","passage","--stage","verify","--out-dir",str(out)])==0
    assert calls==["predict",("verify",125)]
