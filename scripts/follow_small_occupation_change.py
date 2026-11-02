#!/usr/bin/env python3
"""Follow a rare non-target occupation change with a fixed fresh-sample budget."""
import argparse
import json
import math
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from fieldbridge.memory import pulse_inquiry as pulse
from fieldbridge.memory import paired_occupations as paired
from fieldbridge.memory import investigation as inquiry
from fieldbridge.memory import stochastic_passage as st


def save(path,value):
    Path(path).write_text(json.dumps(value,indent=2,allow_nan=False)+"\n",encoding="utf-8")


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--parent",type=Path,required=True)
    parser.add_argument("--out-dir",type=Path,required=True)
    parser.add_argument("--maximum-n",type=int,default=16384)
    args=parser.parse_args()
    reference=json.loads((args.parent/"reference.json").read_text(encoding="utf-8"))
    plan=json.loads((args.parent/"plan.json").read_text(encoding="utf-8"))
    baseline=next(s for s in plan["schedules"] if s["amplitude"]==0)
    late=max(plan["schedules"],key=lambda s:s["lower"])
    candidates=[]
    for row in reference["all_state_contrasts"]:
        if row["left"]!=late["id"] or row["right"]!=baseline["id"] or row["state_index"]==reference["positive_state_index"]:
            continue
        counts=paired.from_summary(row["effect"],row["standard_error"],reference["n"])
        if counts["effect"] and counts["standard_error"]:
            candidates.append({"state_index":row["state_index"],**counts})
    if not candidates:
        raise SystemExit("No non-target occupation change to follow")
    selected=sorted(candidates,key=lambda r:(-abs(r["effect"])/r["standard_error"],r["state_index"]))[0]
    required=math.ceil(16*reference["n"]*selected["standard_error"]**2/selected["effect"]**2)
    n=min(args.maximum_n, max(256,2**math.ceil(math.log2(required))))
    args.out_dir.mkdir(parents=True,exist_ok=False)
    sources=[Path(__file__),Path(pulse.__file__),Path(paired.__file__),Path(st.__file__)]
    hashes={str(p.relative_to(ROOT)):inquiry.file_hash(p) for p in sources}
    frozen={"question":"Does the rare late-pulse redistribution persist outside the original target state?",
            "parent_reference_sha256":inquiry.file_hash(args.parent/"reference.json"),
            "exploratory_hint":selected,"task":plan["task"],"schedules":[baseline,late],
            "reference_n":n,"requested_precision_n":required,"maximum_n":args.maximum_n,
            "seeds":[406000,406002],"source_hashes":hashes,
            "criterion":"same direction and fresh exact paired test adjusted across all final states"}
    save(args.out_dir/"plan.json",frozen)
    prediction=pulse.run(plan["task"],frozen["schedules"],n=4096,seed=406000,scheme="euler")
    save(args.out_dir/"prediction.json",prediction)
    commitment=inquiry.file_hash(args.out_dir/"prediction.json")
    save(args.out_dir/"prediction_commitment.json",{"prediction_sha256":commitment})
    print(json.dumps({"stage":"rare_change_prediction_committed","reference_n":n,"selected_state":selected["state_index"]}),flush=True)
    fresh=pulse.run(plan["task"],frozen["schedules"],n=n,seed=406002,scheme="heun")
    save(args.out_dir/"reference.json",fresh)
    changes=fresh["all_state_contrasts"]
    for row,qvalue in zip(changes,inquiry.adjusted_probabilities([r["p_value"] for r in changes])):
        row["adjusted_p_value"]=qvalue
    found=next(r for r in changes if r["state_index"]==selected["state_index"])
    probability_checks=[st.compare(p["probabilities"],{**r,"n":n},4096,st.Policy())["passed"]
                        for p,r in zip(prediction["results"],fresh["results"])]
    resolved=(fresh["all_valid"] and prediction["all_valid"] and
              found["effect"]*selected["effect"]>0 and found["adjusted_p_value"]<.05)
    report={"complete":True,"selected_state":fresh["states"]["states"][selected["state_index"]],
            "exploratory_hint":selected,"fresh_change":found,"all_state_changes":changes,
            "small_occupation_change_replicated":bool(resolved),
            "all_original_probability_checks_pass":all(probability_checks),
            "reference_n":n,"source_hashes":hashes,"novelty_established":False}
    if commitment!=inquiry.file_hash(args.out_dir/"prediction.json") or any(inquiry.file_hash(ROOT/p)!=h for p,h in hashes.items()):
        raise RuntimeError("A frozen input changed")
    save(args.out_dir/"report.json",report)
    print(json.dumps({k:v for k,v in report.items() if k not in {"all_state_changes","source_hashes"}},indent=2),flush=True)
    return 0 if fresh["all_valid"] else 2


if __name__=="__main__":
    raise SystemExit(main())
