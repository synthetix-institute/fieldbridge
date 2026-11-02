#!/usr/bin/env python3
"""Freeze and test one finite-rate writing candidate, learn a conditional repair, then test new problems.

The two drift families are the coupled pitchfork amplitudes and cubic-anisotropy
vector model already used in the 3D-passage screen. This script changes the
writing rate, preparation and bias, never the final landscape within a family.
No reference outcome enters a predicted equation or its coefficients.
"""
import argparse
from dataclasses import asdict
import hashlib
import json
from multiprocessing import get_context
from pathlib import Path
import platform
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from fieldbridge.memory import stochastic_passage as st


def model(family,D):
    common = {"schema":"fieldbridge-memory/1","kind":"equations","name":family,
              "question":"Predict the final basin probabilities after a finite-rate write and a force-free hold",
              "assumptions":["overdamped polynomial amplitudes","additive Ito Brownian noise"],
              "provenance":{"kind":"constructed drift used in the 20261018 autonomous-search follow-up",
                            "source":"run_nb_and_fb3d_fast.py, eval_3d_regime_job; corrected anisotropy interpretation"},
              "carrier":{"kind":"euclid","variables":["x","y","z"],"scale":2.},
              "parameters":{"mu":.5,"g1":.02,"g2":.25},
              "control":{"name":"mu","range":[-.3,.5]},"noise":D}
    if family=="coupled_amplitudes":
        common["parameters"]["kap"] = .01
        common["drift"]={"x":"mu*x-x**3+kap*(y+z-2*x)",
                         "y":"(mu-g1)*y-y**3+kap*(x+z-2*y)",
                         "z":"(mu-g2)*z-z**3+kap*(x+y-2*z)"}
    elif family=="cubic_anisotropy":
        common["parameters"]["gamma"] = .5
        common["drift"]={"x":"mu*x-(x**2+y**2+z**2)*x-gamma*x**3",
                         "y":"(mu-g1)*y-(x**2+y**2+z**2)*y-gamma*y**3",
                         "z":"(mu-g2)*z-(x**2+y**2+z**2)*z-gamma*z**3"}
    else:
        raise ValueError(family)
    return common


def task(family,rate,preparation,direction):
    import numpy as np
    from statistics import NormalDist
    D = .05*rate
    h = NormalDist().inv_cdf(.8)*np.sqrt(D)*rate**.25/np.pi**.25
    direction = np.asarray(direction,float)
    bias = h*direction/np.linalg.norm(direction)
    return {"schema":st.SCHEMA,"model":model(family,D),"convention":"ito",
            "protocol":{"start":-.3,"end":.5,"rate":rate},"diffusion":D,
            "bias":bias.tolist(),"preparation":preparation,"hold_diffusion":1e-5,
            "relax_time":80.,"hold_time":80.}


def design():
    zero = {"mean":[0.,0.,0.],"covariance":[[0.,0.,0.],[0.,0.,0.],[0.,0.,0.]]}
    prepared = {"mean":[.02,-.015,.01],"covariance":[[1e-4,0.,0.],[0.,2e-4,0.],[0.,0.,1e-4]]}
    jobs = [
        ("dev_fast","development",task("coupled_amplitudes",.02,zero,[1,0,0])),
        ("dev_slow","development",task("coupled_amplitudes",.005,zero,[1,0,0])),
        ("transfer_same_fast","withheld",task("coupled_amplitudes",.0125,prepared,[1,0,0])),
        ("transfer_same_slow","withheld",task("coupled_amplitudes",.0025,prepared,[1,0,0])),
        ("transfer_other_fast","withheld",task("cubic_anisotropy",.0125,prepared,[1,0,0])),
        ("transfer_other_slow","withheld",task("cubic_anisotropy",.0025,prepared,[1,0,0])),
    ]
    return [{"id":name,"split":split,"task":t,"seed":2026102800+17*i} for i,(name,split,t) in enumerate(jobs)]


def worker(args):
    job,policy,out,device,rule,stage = args
    started=time.monotonic()
    try:
        if stage=="prepare":
            result=st.prepare_prediction(job["task"],policy,job["seed"],out,device=device,frozen_rule=rule)
        elif stage=="verify":
            result=st.verify_prediction(out,job["seed"]+2,device=device)
        else:
            result=st.evaluate(job["task"],policy,job["seed"],out,device=device,frozen_rule=rule)
        return {"id":job["id"],"split":job["split"],"seconds":time.monotonic()-started,**result}
    except Exception as exc:
        return {"id":job["id"],"split":job["split"],"error":f"{type(exc).__name__}: {exc}",
                "repair_verified":False,"seconds":time.monotonic()-started}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-dir",type=Path,required=True)
    parser.add_argument("--device",choices=("cpu","cuda"),default="cpu")
    parser.add_argument("--workers",type=int,default=1)
    parser.add_argument("--n",type=int,default=4096)
    parser.add_argument("--development-only",action="store_true")
    args=parser.parse_args()
    if args.out_dir.exists():
        raise SystemExit("Choose a new output directory")
    args.out_dir.mkdir(parents=True)
    policy=st.Policy(n_prediction=args.n,n_reference=args.n)
    jobs=design()
    files=[Path(st.__file__),Path(__file__),ROOT/"fieldbridge/memory/spec.py"]
    source_hashes={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in files}
    plan={"schema":"fieldbridge-stochastic-loop/1","policy":asdict(policy),"jobs":jobs,
          "source_hashes":source_hashes,"python":platform.python_version(),
          "selection":"retain repair only after a resolved linear-prediction failure and verified development repairs",
          "construction_target":{"family":"coupled_amplitudes","state":"all amplitudes positive",
                                 "probability":.60,"confidence":"Wilson 95% lower bound"},
          "novelty_established":False,"complete":False}
    (args.out_dir/"plan.json").write_text(json.dumps(plan,indent=2,allow_nan=False)+"\n", encoding="utf-8")
    plan_hash=st.digest(plan)
    results=[]
    design_record=None
    for split in ("development","withheld"):
        if split=="withheld" and args.development_only:
            break
        if split=="withheld":
            dev=[r for r in results if r["split"]=="development"]
            if not all(r["repair_verified"] for r in dev) or not any(not r.get("linear_test",{}).get("passed",True) for r in dev):
                break
            rule={"schema":st.RULE,"policy":asdict(policy),"source_hashes":source_hashes,
                  "plan_sha256":plan_hash,"development_ids":[r["id"] for r in dev],
                  "applicability":["polynomial real diffusion with a fixed zero reference branch",
                                   "constant isotropic Ito diffusion","matching rms remainder below 0.05",
                                   "checked multistable final drift and force-free retention"],
                  "repair":"retain the calculated nonlinear drift and diffusion after Gaussian-prefix matching"}
            raw=json.dumps(rule,sort_keys=True,indent=2)+"\n"
            (args.out_dir/"frozen_rule.json").write_text(raw, encoding="utf-8")
            frozen_hash=hashlib.sha256(raw.encode()).hexdigest()
        batch=[j for j in jobs if j["split"]==split]
        parameters=[(job,policy,str(args.out_dir/job["id"]),
                     f"cuda:{i%args.workers}" if args.device=="cuda" else "cpu",
                     str(args.out_dir/"frozen_rule.json") if split=="withheld" else None,
                     "prepare" if split=="withheld" else "evaluate") for i,job in enumerate(batch)]
        if split=="withheld":
            with get_context("spawn").Pool(args.workers) as pool:
                prepared=list(pool.map(worker,parameters))
            if any("error" in row for row in prepared):
                results+=prepared
                break
            # Select a protocol solely from committed predictions, before any withheld reference trajectory exists.
            options=[]
            for job in batch:
                if job["task"]["model"]["name"]!="coupled_amplitudes":
                    continue
                record=json.loads((args.out_dir/job["id"]/"prediction.json").read_text(encoding="utf-8"))
                index=next(i for i,q in enumerate(record["states"]["states"]) if min(q)>0)
                options.append({"id":job["id"],"state_index":index,"rate":job["task"]["protocol"]["rate"],
                                "predicted_probability":record["repaired_prediction"]["probabilities"][index],
                                "prediction_sha256":prepared[[r["id"] for r in prepared].index(job["id"])]["prediction_sha256"]})
            chosen=max(options,key=lambda row:row["predicted_probability"])
            design_record={"typed_move":"change the writing rate, keeping the final drift and holding conditions fixed",
                           "target":"all three amplitudes positive with probability at least 0.60",
                           "target_probability":.60,"options":options,"chosen":chosen,
                           "prediction_met_target":chosen["predicted_probability"]>=.60}
            (args.out_dir/"protocol_design.json").write_text(json.dumps(design_record,indent=2)+"\n", encoding="utf-8")
            parameters=[(*p[:-1],"verify") for p in parameters]
        with get_context("spawn").Pool(args.workers) as pool:
            for result in pool.imap_unordered(worker,parameters):
                results.append(result)
                (args.out_dir/"progress.json").write_text(json.dumps(results,indent=2,allow_nan=False)+"\n", encoding="utf-8")
                print(json.dumps({k:result[k] for k in ("id","seconds","repair_verified")},default=str),flush=True)
        if split=="withheld" and hashlib.sha256((args.out_dir/"frozen_rule.json").read_bytes()).hexdigest()!=frozen_hash:
            raise RuntimeError("Rule changed during withheld evaluation")
    withheld=[r for r in results if r["split"]=="withheld"]
    report={"plan_sha256":plan_hash,"results":results,
            "complete":len(results)==len(jobs) and all("reference" in r for r in results),
            "closed_loop_verified":len(withheld)==4 and all(r["repair_verified"] for r in results),
            "novelty_established":False,"source_hashes":source_hashes}
    if design_record is not None and withheld:
        chosen=design_record["chosen"]
        result=next(r for r in withheld if r["id"]==chosen["id"])
        if "reference" in result:
            n=result["reference"]["n"]
            p=result["reference"]["probabilities"][chosen["state_index"]]
            z=1.96
            lower=(p+z*z/(2*n)-z*((p*(1-p)+z*z/(4*n))/n)**.5)/(1+z*z/n)
            report["designed_protocol"]={**design_record,"observed_probability":p,"wilson_95_lower":lower,
                                         "target_verified":design_record["prediction_met_target"] and lower>=.60
                                                           and result["reference"]["valid"]}
    if args.development_only:
        report["scope"]="development only; no withheld outcomes evaluated"
    report["new_physics_claim"]="not established by this validation"
    (args.out_dir/"report.json").write_text(json.dumps(report,indent=2,allow_nan=False)+"\n", encoding="utf-8")
    print(json.dumps({k:v for k,v in report.items() if k!="results"},indent=2))
    return 0 if report["closed_loop_verified"] or args.development_only else 2


if __name__=="__main__":
    raise SystemExit(main())
