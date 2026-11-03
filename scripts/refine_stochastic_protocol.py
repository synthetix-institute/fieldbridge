#!/usr/bin/env python3
"""Prediction-led rate refinement after a frozen protocol search misses its unchanged physical target.

At most two lower-rate proposals are allowed. Noise and signal follow the
derived scaling; a selected proposal and two factorial controls are committed
before any new reference trajectories are evaluated. A failed reference stops
the campaign rather than selecting another candidate from its outcomes.
"""
import argparse
import copy
from dataclasses import asdict
import hashlib
import json
from multiprocessing import get_context
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from fieldbridge.memory import stochastic_passage as st


def lower_bound(p,n):
    z=1.96
    return (p+z*z/(2*n)-z*((p*(1-p)+z*z/(4*n))/n)**.5)/(1+z*z/n)


def rescale(task,rate):
    new=copy.deepcopy(task)
    factor=rate/task["protocol"]["rate"]
    new["protocol"]["rate"]=rate
    new["diffusion"]*=factor
    new["model"]["noise"]=new["diffusion"]
    new["bias"]=[h*factor**.75 for h in task["bias"]]
    return new


def prepare_job(args):
    t,policy,seed,out,device,rule=args
    return st.prepare_prediction(t,policy,seed,out,device=device,frozen_rule=rule)


def verify_job(args):
    out,seed,device=args
    return st.verify_prediction(out,seed,device=device)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--parent",type=Path,required=True)
    parser.add_argument("--out-dir",type=Path,required=True)
    parser.add_argument("--device",default="cpu")
    parser.add_argument("--workers",type=int,default=2)
    args=parser.parse_args()
    parent=json.loads((args.parent/"report.json").read_text(encoding="utf-8"))
    if not parent["closed_loop_verified"]:
        raise SystemExit("The frozen repair must pass its independent transfer tests first")
    if args.out_dir.exists():
        raise SystemExit("Choose a new output directory")
    args.out_dir.mkdir(parents=True)
    original=parent["designed_protocol"]
    chosen=original["chosen"]
    prior=json.loads((args.parent/chosen["id"]/"prediction.json").read_text(encoding="utf-8"))
    policy=st.Policy(**prior["policy"])
    rule=args.parent/"frozen_rule.json"
    frozen=hashlib.sha256(rule.read_bytes()).hexdigest()
    target=original["target_probability"]
    plan={"parent_report_sha256":hashlib.sha256((args.parent/"report.json").read_bytes()).hexdigest(),
          "frozen_rule_sha256":frozen,"unchanged_target":target,
          "initial_task":prior["task"],"policy":asdict(policy),"maximum_new_candidates":2,
          "minimum_rate":.000625,"proposal":"halve the best predicted rate; retain scaled D/r and h/r^(3/4)",
          "controls":["new rate with old writing noise and signal","old rate with new writing noise and signal"],
          "implementation_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    (args.out_dir/"plan.json").write_text(json.dumps(plan,indent=2)+"\n",encoding="utf-8")
    attempts=[];selected=None
    for k in range(plan["maximum_new_candidates"]):
        rate=chosen["rate"]/(2**(k+1))
        if rate<plan["minimum_rate"]:
            break
        task=rescale(prior["task"],rate)
        seed=2026102900+17*k
        out=args.out_dir/f"proposal_{k}"
        st.prepare_prediction(task,policy,seed,out,device=args.device,frozen_rule=rule)
        prediction=json.loads((out/"prediction.json").read_text(encoding="utf-8"))
        idx=next(i for i,q in enumerate(prediction["states"]["states"]) if min(q)>0)
        p=prediction["repaired_prediction"]["probabilities"][idx]
        row={"rate":rate,"probability":p,"prediction_lower_95":lower_bound(p,policy.n_prediction),
             "task_sha256":st.digest(task),"prediction_sha256":hashlib.sha256((out/"prediction.json").read_bytes()).hexdigest()}
        attempts.append(row)
        print(json.dumps(row),flush=True)
        if row["prediction_lower_95"]>=target:
            selected={"task":task,"seed":seed,"out":str(out),"index":idx,**row}
            break
    if selected is None:
        (args.out_dir/"report.json").write_text(json.dumps({"status":"target_not_reached_by_predictions", "attempts":attempts,
            "target":target,"complete":True,"novelty_established":False},indent=2)+"\n")
        return 2
    old=prior["task"]
    rate_only=copy.deepcopy(old);rate_only["protocol"]["rate"]=selected["rate"]
    forcing_only=copy.deepcopy(selected["task"]);forcing_only["protocol"]["rate"]=old["protocol"]["rate"]
    controls=[("rate_only",rate_only,2026103001),("noise_signal_only",forcing_only,2026103002)]
    with get_context("spawn").Pool(args.workers) as pool:
        pool.map(prepare_job,[(t,policy,seed,str(args.out_dir/name),args.device,rule) for name,t,seed in controls])
    design={"selected":selected,"attempts":attempts,"target":target,"controls":[
        {"name":name,"task":t,"seed":seed,"prediction_sha256":hashlib.sha256(
            (args.out_dir/name/"prediction.json").read_bytes()).hexdigest()} for name,t,seed in controls]}
    (args.out_dir/"committed_design.json").write_text(json.dumps(design,indent=2)+"\n",encoding="utf-8")
    jobs=[("selected",selected["out"],selected["seed"]+2)]+[(name,str(args.out_dir/name),seed+2) for name,_,seed in controls]
    with get_context("spawn").Pool(args.workers) as pool:
        rows=pool.map(verify_job,[(out,seed,args.device) for _,out,seed in jobs])
    results={name:r for (name,_,_),r in zip(jobs,rows)}
    result=results["selected"]
    p=result["reference"]["probabilities"][selected["index"]]
    target_verified=result["repair_verified"] and lower_bound(p,policy.n_reference)>=target
    report={"status":"verified_design" if target_verified else "design_not_verified", "complete":True,
            "target_probability":target,"selected_rate":selected["rate"],"predicted_probability":selected["probability"],
            "observed_probability":p,"observed_lower_95":lower_bound(p,policy.n_reference),"target_verified":target_verified,
            "controls":results,"parent_observed_probability":original["observed_probability"],"novelty_established":False}
    if hashlib.sha256(rule.read_bytes()).hexdigest()!=frozen:
        raise RuntimeError("Rule changed during design verification")
    (args.out_dir/"report.json").write_text(json.dumps(report,indent=2,allow_nan=False)+"\n",encoding="utf-8")
    print(json.dumps({k:v for k,v in report.items() if k!="controls"},indent=2))
    return 0 if target_verified else 2


if __name__=="__main__":
    raise SystemExit(main())
