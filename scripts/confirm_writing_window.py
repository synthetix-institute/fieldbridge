#!/usr/bin/env python3
"""Test whether the selected timing benefit is lost with a late equal-dose pulse."""
import argparse
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from fieldbridge.memory import pulse_inquiry as pulse
from fieldbridge.memory import investigation as inquiry
from fieldbridge.memory import stochastic_passage as st
from fieldbridge.memory import paired_occupations

loader=importlib.util.spec_from_file_location("writing_inquiry",Path(__file__).with_name("investigate_writing_path.py"))
study=importlib.util.module_from_spec(loader)
loader.loader.exec_module(study)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--parent",type=Path,required=True)
    parser.add_argument("--out-dir",type=Path,required=True)
    parser.add_argument("--n",type=int,default=4096)
    args=parser.parse_args()
    parent=json.loads((args.parent/"report.json").read_text(encoding="utf-8"))
    if not parent["timing_effect_confirmed"]:
        raise SystemExit("The first timing contrast must be resolved before this follow-up")
    plan=json.loads((args.parent/"plan.json").read_text(encoding="utf-8"))
    design=json.loads((args.parent/"committed_design.json").read_text(encoding="utf-8"))
    high=next(s for s in design["selected"] if s["id"]==design["high"])
    late=max(plan["schedules"],key=lambda s:s["lower"])
    baseline=plan["schedules"][0]
    if late["lower"]<=high["lower"]:
        raise SystemExit("The exploratory result did not select an interior window")
    schedules=[baseline,high,late]
    args.out_dir.mkdir(parents=True,exist_ok=False)
    hashes={str(p.relative_to(ROOT)):inquiry.file_hash(p) for p in [Path(__file__),Path(pulse.__file__),Path(inquiry.__file__),Path(st.__file__),Path(study.__file__),Path(paired_occupations.__file__)]}
    frozen={"question":"Does an equal-dose late pulse lose the selected occupation benefit?",
            "parent_report_sha256":inquiry.file_hash(args.parent/"report.json"),"task":plan["task"],
            "schedules":schedules,"n":args.n,"seeds":[404000,404002],"source_hashes":hashes,
            "late_equivalence_margin":.02,
            "all_states":"Report all paired occupation changes, not only the original all-positive target."}
    study.save(args.out_dir/"plan.json",frozen)
    prediction=pulse.run(plan["task"],schedules,n=args.n,seed=404000,scheme="euler")
    study.save(args.out_dir/"prediction.json",prediction)
    commitment=inquiry.file_hash(args.out_dir/"prediction.json")
    study.save(args.out_dir/"prediction_commitment.json",{"prediction_sha256":commitment})
    idx=prediction["positive_state_index"]
    print(json.dumps({"stage":"late_test_prediction_committed","probabilities":{r["id"]:r["probabilities"][idx] for r in prediction["results"]}}),flush=True)
    reference=pulse.run(plan["task"],schedules,n=args.n,seed=404002,scheme="heun")
    study.save(args.out_dir/"reference.json",reference)
    paths=study.export_comparisons(prediction,reference,args.out_dir/"comparisons",hashes["fieldbridge/memory/pulse_inquiry.py"])
    inquiry.write_report(paths,args.out_dir/"residual_investigation")
    def contrast(left,right):
        row=next(r for r in reference["contrasts"] if {r["left"],r["right"]}=={left,right})
        sign=1 if row["left"]==left else -1
        effect=sign*row["effect"]
        se=row["standard_error"]
        return {"effect":effect,"standard_error":se,
                "lower_95":row["lower_95"] if sign==1 else -row["upper_95"],
                "upper_95":row["upper_95"] if sign==1 else -row["lower_95"],"p_value":row["p_value"]}
    middle_late=contrast(high["id"],late["id"])
    late_base=contrast(late["id"],baseline["id"])
    all_states=reference["all_state_contrasts"]
    pvalues=[c["p_value"] for c in all_states]
    for row,q in zip(all_states,inquiry.adjusted_probabilities(pvalues)):
        row["adjusted_p_value"]=q
        row["resolved_occupation_change"]=q<.05 and abs(row["effect"])>2*max(r["step_refinement_difference"] for r in reference["results"])
    c,M,w=pulse.mixing_witness(plan["task"])
    folds={str(mu):pulse.mixed_state_folds(c,M,mu) for mu in (.12,.16,.20)}
    success=(prediction["all_valid"] and reference["all_valid"] and
             middle_late["lower_95"]>0 and late_base["lower_95"]>-.02 and late_base["upper_95"]<.02)
    report={"complete":True,"window_effect_confirmed":success,
            "probabilities":{r["id"]:r["probabilities"][idx] for r in reference["results"]},
            "middle_minus_late":middle_late,"late_minus_no_pulse":late_base,
            "all_state_changes":all_states,"mixing_witness":w,"mixed_state_folds":folds,
            "source_hashes":hashes,"novelty_established":False}
    if commitment!=inquiry.file_hash(args.out_dir/"prediction.json") or any(inquiry.file_hash(ROOT/p)!=h for p,h in hashes.items()):
        raise RuntimeError("A frozen input changed during the follow-up")
    study.save(args.out_dir/"report.json",report)
    print(json.dumps({k:v for k,v in report.items() if k in {"complete","window_effect_confirmed","probabilities","middle_minus_late","late_minus_no_pulse"}},indent=2),flush=True)
    return 0 if reference["all_valid"] and prediction["all_valid"] else 2


if __name__=="__main__":
    raise SystemExit(main())
