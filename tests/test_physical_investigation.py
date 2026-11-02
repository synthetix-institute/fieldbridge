import copy
import hashlib
import importlib.util
import json
from pathlib import Path

import numpy as np
import pytest

from fieldbridge.memory import investigation as inv
from fieldbridge.memory import pulse_inquiry as pulse
from fieldbridge.memory import stochastic_passage as st
from fieldbridge.memory import paired_occupations as paired

loader=importlib.util.spec_from_file_location("passage_study",Path(__file__).resolve().parents[1]/"scripts/run_stochastic_passage_loop.py")
study=importlib.util.module_from_spec(loader)
loader.loader.exec_module(study)


def record(tmp_path,name,observed=.51,seed=1,valid=True,numerical=0.):
    out=tmp_path/name
    out.mkdir()
    pred={"task":{},"states":{"states":[[-1],[1]]},"implementation_sha256":"fixed",
          "policy":{"n":1000000},"repaired_prediction":{"probabilities":[.5,.5],
          "coarse_probabilities":[.5,.5],"n":1000000,"seed":seed,"valid":True}}
    raw=json.dumps(pred)
    (out/"prediction.json").write_text(raw)
    comparison={"prediction_sha256":hashlib.sha256(raw.encode()).hexdigest(),"task_sha256":"same_task",
                "repair_verified":True,"reference":{"probabilities":[observed,1-observed],
                "coarse_probabilities":[observed+numerical,1-observed-numerical],
                "n":1000000,"seed":seed+100,"valid":valid}}
    (out/"comparison.json").write_text(json.dumps(comparison))
    return out/"comparison.json"


def test_small_departures_are_investigated_even_when_original_prediction_passes(tmp_path):
    paths=[record(tmp_path,"a",seed=1),record(tmp_path,"b",observed=.511,seed=2)]
    result=inv.investigate(paths)
    assert result["repeatable_departures"]==2
    assert all(q["accepted_original_cases"]==2 for q in result["questions"])
    assert result["questions"][0]["adjusted_p_value"]<.05


def test_copied_records_and_reused_seed_pairs_are_not_independent_repeats(tmp_path):
    a=record(tmp_path,"a")
    b=record(tmp_path,"b")
    result=inv.investigate([a,a,b])
    assert len(result["cases"])==1
    assert result["repeatable_departures"]==0


def test_reversing_residuals_do_not_establish_a_repeatable_effect(tmp_path):
    a=record(tmp_path,"a",observed=.51,seed=1)
    b=record(tmp_path,"b",observed=.49,seed=2)
    result=inv.investigate([a,b])
    assert result["repeatable_departures"]==0
    assert all(not q["same_sign"] for q in result["questions"])


@pytest.mark.parametrize("valid,numerical",[(False,0.),(True,.02)])
def test_unresolved_numerics_prevent_a_physical_flag(tmp_path,valid,numerical):
    result=inv.investigate([record(tmp_path,"a",valid=valid,numerical=numerical)])
    assert all(q["status"]=="check_numerical_resolution" for q in result["questions"])


def test_modified_predictions_are_rejected(tmp_path):
    path=record(tmp_path,"a")
    path.with_name("prediction.json").write_text("{}")
    with pytest.raises(ValueError,match="another prediction"):
        inv.investigate([path])


def test_multiple_comparison_adjustment_and_bounded_sampling():
    assert inv.adjusted_probabilities([.01,.04,.03])==pytest.approx([.03,.04,.04])
    assert inv.probability_variance(0.,1000)>0


def test_no_discordance_has_a_finite_uncertainty_interval():
    c=paired.contrast(0,0,4096)
    assert c["effect"]==0 and c["p_value"]==1
    assert c["lower_95"]<0<c["upper_95"]
    assert c["upper_95"]<.02


def test_five_rare_changes_do_not_become_significant_by_gaussian_approximation():
    c=paired.contrast(5,0,4096)
    assert c["p_value"]==pytest.approx(.0625)
    assert c["lower_95"]<0
    recovered=paired.from_summary(c["effect"],c["standard_error"],4096)
    assert recovered==c


def test_public_investigation_command_keeps_accepted_predictions(tmp_path):
    from fieldbridge.cli import main
    a=record(tmp_path,"a")
    b=record(tmp_path,"b",seed=2)
    out=tmp_path/"followup"
    assert main(["memory","investigate",str(a),str(b),"--out-dir",str(out)])==0
    report=json.loads((out/"investigation.json").read_text())
    assert report["repeatable_departures"]==2


def task():
    return copy.deepcopy(study.design()[3]["task"])


def test_source_equations_prove_mean_invariance_and_variance_contraction():
    c,M,w=pulse.mixing_witness(task())
    assert w["direct_mean_force"]=="0"
    assert w["complete_graph_variance_identity"]
    assert w["noncommuting_mean_dynamics"]
    assert np.allclose(np.linalg.eigvalsh(M),[-3,-3,0])


def test_a_pulse_with_direct_mean_force_is_rejected():
    t=task()
    t["model"]["drift"]["x"]+="+kap*x"
    with pytest.raises(ValueError,match="direct force"):
        pulse.mixing_witness(t)


def test_mixed_sign_branch_has_a_nondegenerate_fold_below_the_applied_pulse():
    c,M,w=pulse.mixing_witness(task())
    folds=pulse.mixed_state_folds(c,M,.16)
    assert len(folds)>=1
    assert any(.01<f["coupling"]<.03 for f in folds)
    assert all(f["equation_residual"]<1e-8 and f["equilibria_on_lower_coupling_side"] for f in folds)


def test_shared_noise_makes_identical_schedules_identical():
    pytest.importorskip("torch")
    t=task()
    c=st.compile_task(t)
    q=st.final_states(c)["states"][-1]
    t["protocol"]={"start":.49,"end":.5,"rate":1.}
    t["preparation"]={"mean":q,"covariance":(np.eye(3)*1e-8).tolist()}
    t["relax_time"]=0.;t["hold_time"]=0.
    schedules=[{"id":name,"lower":.491,"upper":.499,"amplitude":0.} for name in ["one","two"]]
    result=pulse.run(t,schedules,n=256,seed=12)
    assert result["results"][0]["probabilities"]==result["results"][1]["probabilities"]
    assert result["contrasts"][0]["effect"]==0
    assert result["contrasts"][0]["standard_error"]==0
    assert len(result["all_state_contrasts"])==len(result["states"]["states"])
    assert result["all_valid"]
