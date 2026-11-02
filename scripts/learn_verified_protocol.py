#!/usr/bin/env python3
"""Retain the simplest goal-reaching protocol from factorial controls, then test its reuse on a fresh task."""
import argparse
import copy
import hashlib
import importlib.util
import json
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from fieldbridge.memory import stochastic_passage as st

loader=importlib.util.spec_from_file_location("protocol_refinement",Path(__file__).with_name("refine_stochastic_protocol.py"))
refinement=importlib.util.module_from_spec(loader)
loader.loader.exec_module(refinement)


def changed_settings(before,after):
    fixed=("drift","parameters","carrier")
    if any(before["model"][k]!=after["model"][k] for k in fixed):
        raise ValueError("A protocol-only recipe must preserve the supplied generator drift")
    if any(before[k]!=after[k] for k in ("hold_diffusion","hold_time","relax_time")):
        raise ValueError("Holding conditions changed")
    if any(before["protocol"][k]!=after["protocol"][k] for k in ("start","end")):
        raise ValueError("The final landscape or start control changed")
    return [name for name,a,b in (("rate",before["protocol"]["rate"],after["protocol"]["rate"]),
                                 ("diffusion",before["diffusion"],after["diffusion"]),
                                 ("bias",before["bias"],after["bias"])) if a!=b]


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--parent",type=Path,required=True)
    parser.add_argument("--design",type=Path,required=True)
    parser.add_argument("--out-dir",type=Path,required=True)
    parser.add_argument("--device",default="cpu")
    args=parser.parse_args()
    if args.out_dir.exists():
        raise SystemExit("Choose a new output directory")
    args.out_dir.mkdir(parents=True)
    report=json.loads((args.design/"report.json").read_text())
    committed=json.loads((args.design/"committed_design.json").read_text())
    parent=json.loads((args.parent/"report.json").read_text())
    source=json.loads((args.parent/parent["designed_protocol"]["chosen"]["id"]/"prediction.json").read_text())["task"]
    target=report["target_probability"]
    if not report["target_verified"]:
        raise SystemExit("The construction did not reach its stated target")
    options=[]
    for name,result in report["controls"].items():
        path=Path(committed["selected"]["out"]) if name=="selected" else args.design/name
        prediction=json.loads((path/"prediction.json").read_text())
        idx=next(i for i,q in enumerate(prediction["states"]["states"]) if min(q)>0)
        changes=changed_settings(source,prediction["task"])
        p=result["reference"]["probabilities"][idx]
        lower=refinement.lower_bound(p,result["reference"]["n"])
        pp=prediction["repaired_prediction"]["probabilities"][idx]
        eligible=(result["repair_verified"] and lower>=target and
                  refinement.lower_bound(pp,prediction["repaired_prediction"]["n"])>=target)
        options.append({"name":name,"changes":changes,"eligible":eligible,"observed_probability":p,
                        "lower_95":lower,"task":prediction["task"],
                        "comparison_sha256":hashlib.sha256((path/"comparison.json").read_bytes()).hexdigest()})
    selected=min((r for r in options if r["eligible"]),key=lambda r:len(r["changes"]))
    if selected["changes"]!=["rate"]:
        raise SystemExit("The controls did not support a rate-only recipe")
    factor=selected["task"]["protocol"]["rate"]/source["protocol"]["rate"]
    recipe={"schema":"fieldbridge-verified-protocol/1","operation":"multiply writing rate",
            "factor":factor,"preserve":["drift","diffusion","bias","holding_conditions"],
            "source_drift_sha256":st.digest(source["model"]["drift"]),"target_probability":target,
            "scope":"verified on the recorded source conditions; a new application is a proposal until tested",
            "evidence":selected,"comparison":options}
    raw=json.dumps(recipe,sort_keys=True,indent=2)+"\n"
    (args.out_dir/"learned_recipe.json").write_text(raw)
    recipe_hash=hashlib.sha256(raw.encode()).hexdigest()
    # A new base rate and preparation were not part of either previous comparison.
    task=copy.deepcopy(source)
    task["protocol"]["rate"]=.002
    task["preparation"]={"mean":[-.015,.0125,-.005],
                         "covariance":[[2e-4,0,0],[0,1e-4,0],[0,0,2e-4]]}
    task["protocol"]["rate"]*=recipe["factor"]
    (args.out_dir/"new_task.json").write_text(json.dumps(task,indent=2)+"\n")
    frozen=args.parent/"frozen_rule.json"
    policy=st.Policy(**json.loads(frozen.read_text())["policy"])
    out=args.out_dir/"transfer"
    st.prepare_prediction(task,policy,2026103100,out,device=args.device,frozen_rule=frozen)
    prediction=json.loads((out/"prediction.json").read_text())
    idx=next(i for i,q in enumerate(prediction["states"]["states"]) if min(q)>0)
    pp=prediction["repaired_prediction"]["probabilities"][idx]
    print(json.dumps({"stage":"committed_reuse_prediction","rate":task["protocol"]["rate"],
                      "predicted_probability":pp,"learned_recipe_sha256":recipe_hash}),flush=True)
    result=st.verify_prediction(out,2026103102,device=args.device)
    p=result["reference"]["probabilities"][idx]
    lower=refinement.lower_bound(p,result["reference"]["n"])
    verified=(result["repair_verified"] and lower>=target and
              refinement.lower_bound(pp,policy.n_prediction)>=target)
    if hashlib.sha256((args.out_dir/"learned_recipe.json").read_bytes()).hexdigest()!=recipe_hash:
        raise RuntimeError("Recipe changed during its test")
    summary={"learned_recipe_sha256":recipe_hash,"changes_retained":recipe["preserve"],
             "rate_multiplier":factor,"new_rate":task["protocol"]["rate"],"predicted_probability":pp,
             "observed_probability":p,"lower_95":lower,"verified_reuse":verified,
             "comparison":result,"novelty_established":False}
    (args.out_dir/"report.json").write_text(json.dumps(summary,indent=2,allow_nan=False)+"\n")
    print(json.dumps({k:v for k,v in summary.items() if k!="comparison"},indent=2))
    return 0 if verified else 2


if __name__=="__main__":
    raise SystemExit(main())
