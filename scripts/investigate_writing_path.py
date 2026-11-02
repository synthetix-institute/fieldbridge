#!/usr/bin/env python3
"""Investigate equal-dose coupling pulses without changing the final model.

The question and pulse family are supplied. A bounded exploratory calculation
selects a timing contrast; its probability predictions are frozen before fresh
reference trajectories at the original and a new rate/preparation are used.
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
from fieldbridge.memory import pulse_inquiry as pulse
from fieldbridge.memory import stochastic_passage as st
from fieldbridge.memory import investigation as inquiry
from fieldbridge.memory import paired_occupations


def save(path,value):
    Path(path).write_text(json.dumps(value,indent=2,allow_nan=False)+"\n",encoding="utf-8")


def worker(args):
    task,schedules,n,seed,scheme=args
    return pulse.run(task,schedules,n=n,seed=seed,scheme=scheme)


def export_comparisons(prediction,reference,out,implementation):
    paths=[]
    for p,r in zip(prediction["results"],reference["results"]):
        task=copy.deepcopy(prediction["task"])
        task["mixing_pulse"]=p["pulse"]
        record={"task":task,"task_sha256":st.digest(task),"states":prediction["states"],
                "probability_prediction":{**p,"seed":prediction["seed"],"n":prediction["n"]},
                "implementation_sha256":implementation,"policy":{"dt":prediction["dt"],"n":prediction["n"]}}
        directory=Path(out)/p["id"]
        directory.mkdir(parents=True,exist_ok=False)
        save(directory/"prediction.json",record)
        ref={**r,"seed":reference["seed"],"n":reference["n"]}
        comparison=st.compare(p["probabilities"],ref,prediction["n"],st.Policy())
        save(directory/"comparison.json",{"task_sha256":record["task_sha256"],
             "prediction_sha256":inquiry.file_hash(directory/"prediction.json"),"reference":ref,
             "prediction_verified":comparison["passed"] and p["valid"],"probability_comparison":comparison})
        paths.append(directory/"comparison.json")
    return paths


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--task",type=Path,required=True)
    parser.add_argument("--out-dir",type=Path,required=True)
    parser.add_argument("--pilot-n",type=int,default=1024)
    parser.add_argument("--n",type=int,default=4096)
    parser.add_argument("--workers",type=int,default=2)
    args=parser.parse_args()
    args.out_dir.mkdir(parents=True,exist_ok=False)
    task=json.loads(args.task.read_text(encoding="utf-8"))
    # A single physical question, not a search over equations or final states.
    centres=[-.08,.04,.16,.28,.40]
    schedules=[{"id":"no_pulse","lower":-.12,"upper":-.04,"amplitude":0.}]
    schedules += [{"id":f"pulse_{i}","lower":round(c-.04,8),"upper":round(c+.04,8),"amplitude":.02}
                  for i,c in enumerate(centres)]
    files=[Path(__file__),Path(pulse.__file__),Path(st.__file__),Path(inquiry.__file__),Path(paired_occupations.__file__)]
    hashes={str(p.relative_to(ROOT)):inquiry.file_hash(p) for p in files}
    plan={"question":"Can the timing of a zero-direct-mean-force mixing pulse change final pattern selection?",
          "task":task,"schedules":schedules,"pilot_n":args.pilot_n,"validation_n":args.n,
          "selection":"largest and smallest predicted all-positive occupation among equal-dose nonzero pulses",
          "fresh_tests":["same task, new trajectories","half the rate and a new preparation"],
          "seeds":[401000,402000,403000,402002,403002],"source_hashes":hashes,
          "claim":"numerical timing effect under identical final drift; novelty assessed separately"}
    save(args.out_dir/"plan.json",plan)
    pilot=pulse.run(task,schedules,n=args.pilot_n,seed=401000)
    save(args.out_dir/"pilot.json",pilot)
    idx=pilot["positive_state_index"]
    eligible=[r for r in pilot["results"] if r["id"]!="no_pulse" and r["valid"]]
    if len(eligible)<2:
        save(args.out_dir/"report.json",{"status":"pilot_requires_numerical_resolution","pilot":pilot["all_valid"]})
        return 2
    low=min(eligible,key=lambda r:r["probabilities"][idx])
    high=max(eligible,key=lambda r:r["probabilities"][idx])
    if low["id"]==high["id"]:
        save(args.out_dir/"report.json",{"status":"no_timing_contrast_in_pilot"})
        return 0
    selected=[schedules[0],low["pulse"],high["pulse"]]
    transfer=copy.deepcopy(task)
    transfer["protocol"]["rate"]/=2
    transfer["preparation"]={"mean":[-.015,.0125,-.005],"covariance":[[2e-4,0,0],[0,1e-4,0],[0,0,2e-4]]}
    design={"selected":selected,"low":low["id"],"high":high["id"],"tasks":[task,transfer],
            "source_hashes":hashes,"pilot":inquiry.file_hash(args.out_dir/"pilot.json")}
    save(args.out_dir/"committed_design.json",design)
    with get_context("spawn").Pool(args.workers) as pool:
        predictions=pool.map(worker,[(t,selected,args.n,seed,"euler")
                                    for t,seed in zip(design["tasks"],[402000,403000])])
    # These full records, including all predicted probabilities, predate fresh references.
    for i,p in enumerate(predictions):
        save(args.out_dir/f"prediction_{i}.json",p)
    save(args.out_dir/"prediction_commitment.json",{f"prediction_{i}.json":inquiry.file_hash(args.out_dir/f"prediction_{i}.json")
                                                  for i in range(len(predictions))})
    print(json.dumps({"stage":"predictions_committed","low":low["id"],"high":high["id"],
          "probabilities":[{r["id"]:r["probabilities"][idx] for r in p["results"]} for p in predictions]}),flush=True)
    with get_context("spawn").Pool(args.workers) as pool:
        references=pool.map(worker,[(t,selected,args.n,seed,"heun")
                                   for t,seed in zip(design["tasks"],[402002,403002])])
    comparisons=[];result=[]
    for i,(p,r) in enumerate(zip(predictions,references)):
        save(args.out_dir/f"reference_{i}.json",r)
        comparisons+=export_comparisons(p,r,args.out_dir/f"comparisons_{i}",hashes["fieldbridge/memory/pulse_inquiry.py"])
        contrast=next(c for c in r["contrasts"] if {c["left"],c["right"]}=={high["id"],low["id"]})
        sign=1 if contrast["left"]==high["id"] else -1
        effect=sign*contrast["effect"]
        z=abs(effect)/contrast["standard_error"] if contrast["standard_error"] else 0.
        valid=p["all_valid"] and r["all_valid"]
        result.append({"rate":r["task"]["protocol"]["rate"],"timing_gain":effect,
                       "paired_standard_error":contrast["standard_error"],"z":z,"numerically_valid":valid,
                       "timing_gain_resolved":bool(valid and effect>0 and contrast["p_value"]<.005),
                       "probabilities":{s["id"]:s["probabilities"][idx] for s in r["results"]},
                       "equal_nonzero_pulse_dose":high["pulse"]["amplitude"]*(high["pulse"]["upper"]-high["pulse"]["lower"])/r["task"]["protocol"]["rate"]})
    if any(inquiry.file_hash(ROOT/p)!=h for p,h in hashes.items()):
        raise RuntimeError("Implementation changed during the calculation")
    inquiry.write_report(comparisons,args.out_dir/"residual_investigation")
    c,M,witness=pulse.mixing_witness(task)
    instant={s["id"]:{"without_pulse":pulse.instantaneous_states(c,M,(s["lower"]+s["upper"])/2,0.),
                       "during_pulse":pulse.instantaneous_states(c,M,(s["lower"]+s["upper"])/2,s["amplitude"])}
             for s in selected[1:]}
    report={"complete":True,"results":result,"mixing_witness":witness,"instantaneous_states":instant,
            "timing_effect_confirmed":all(r["timing_gain_resolved"] for r in result),
            "source_hashes":hashes,"novelty_established":False,
            "autonomy_scope":"The question and pulse family are supplied; timing selection, fresh predictions and checks are automatic."}
    save(args.out_dir/"report.json",report)
    print(json.dumps({k:v for k,v in report.items() if k not in {"mixing_witness","instantaneous_states","source_hashes"}},indent=2),flush=True)
    return 0 if all(r["numerically_valid"] for r in result) else 2


if __name__=="__main__":
    raise SystemExit(main())
