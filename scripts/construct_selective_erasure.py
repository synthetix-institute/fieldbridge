#!/usr/bin/env python3
"""Construct a low-noise selective-erasure pulse from a split pair of folds.

The source equations and desired distinction between mixed-sign states are
given. Fold thresholds set the pulse amplitude, local passage time and noise
scale before independent trajectory tests.
"""
import argparse
import copy
import json
import math
from pathlib import Path
import sys

import numpy as np
from scipy.integrate import solve_ivp

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from fieldbridge.memory import pulse_inquiry as pulse
from fieldbridge.memory import stochastic_passage as st
from fieldbridge.memory import investigation as inquiry


def save(path,value):
    Path(path).write_text(json.dumps(value,indent=2,allow_nan=False)+"\n",encoding="utf-8")


def deterministic(c,M,initial,peak,duration,mu,method):
    bias=c["bias"]
    def evolve(q,time,field):
        sol=solve_ivp(field,(0,time),q,method=method,rtol=1e-10,atol=1e-12)
        if not sol.success:
            raise ArithmeticError(sol.message)
        return sol.y[:,-1]
    q=evolve(initial,duration,lambda t,q:c["drift"](q,mu)[0]+peak*(M@q)+bias)
    ramp=(c["end"]-mu)/c["rate"]
    q=evolve(q,ramp,lambda t,q:c["drift"](q,mu+c["rate"]*t)[0]+bias)
    return evolve(q,c["relax_T"],lambda t,q:c["drift"](q,c["end"])[0])


def simulate(c,M,initial,schedules,mu,duration,short,D,n,seed,scheme):
    import torch
    torch.set_num_threads(1)
    field=st.torch_drift(c,"cpu")
    states=st.final_states(c)
    k,s=len(schedules),len(initial)
    q=torch.tensor(np.broadcast_to(np.asarray(initial)[None,:,None,:],(k,s,n,3)).copy(),dtype=torch.float64)
    fine,coarse=q.clone(),q.clone()
    matrix=torch.tensor(M,dtype=torch.float64)
    bias=torch.tensor(c["bias"],dtype=torch.float64)
    generator=torch.Generator().manual_seed(seed)
    def rhs(q,control,amp,writing):
        base=field(q.reshape(-1,3),torch.tensor(control,dtype=torch.float64)).reshape(q.shape)
        return base+amp*(q@matrix.T)+(bias if writing else 0.)
    def step(q,control,dt,dw,amp,writing,moving):
        f=rhs(q,control,amp,writing)
        if scheme=="euler":
            return q+dt*f+dw
        proposal=q+dt*f+dw
        return q+.5*dt*(f+rhs(proposal,control+c["rate"]*dt if moving else control,amp,writing))+dw
    phases=[("pause",short,0.),("pause",duration-short,short),
            ("ramp",(c["end"]-mu)/c["rate"],0.),("relax",c["relax_T"],0.),("hold",c["hold_T"],0.)]
    pre=None
    with torch.no_grad():
        for phase,length,offset in phases:
            amps=[r["amplitude"] if phase=="pause" and (not r["short"] or offset==0.) else 0. for r in schedules]
            amp=torch.tensor(amps,dtype=torch.float64)[:,None,None,None]
            writing=phase in {"pause","ramp"}
            noise=D if phase!="relax" else 0.
            count=max(1,math.ceil(length/.02))
            dt=length/count
            for j in range(count):
                w1=torch.randn((1,s,n,3),generator=generator,dtype=torch.float64)*np.sqrt(noise*dt)
                w2=torch.randn((1,s,n,3),generator=generator,dtype=torch.float64)*np.sqrt(noise*dt)
                control=mu+c["rate"]*j*dt if phase=="ramp" else (mu if phase=="pause" else c["end"])
                fine=step(fine,control,dt/2,w1,amp,writing,phase=="ramp")
                fine=step(fine,control+c["rate"]*dt/2 if phase=="ramp" else control,dt/2,w2,amp,writing,phase=="ramp")
                coarse=step(coarse,control,dt,w1+w2,amp,writing,phase=="ramp")
            if phase=="relax":
                pre=fine.numpy().copy()
    labels=np.array([[st.labels_of(q,states) for q in group] for group in fine.numpy()])
    old=np.array([[st.labels_of(q,states) for q in group] for group in pre])
    rough=np.array([[st.labels_of(q,states) for q in group] for group in coarse.numpy()])
    rows=[]
    for i,protocol in enumerate(schedules):
        for j,source in enumerate(initial):
            pf=np.asarray(st.probabilities(labels[i,j],len(states["states"])))
            pc=np.asarray(st.probabilities(rough[i,j],len(states["states"])))
            pair=[int(np.sign(source[0])),int(np.sign(source[1]))]
            match=np.array([[int(np.sign(q[0])),int(np.sign(q[1]))]==pair for q in states["states"]])
            positive=np.array([q[0]>0 and q[1]>0 for q in states["states"]])
            unassigned=float(max(np.mean(labels[i,j]<0),np.mean(old[i,j]<0)))
            switching=float(np.mean(labels[i,j]!=old[i,j]))
            difference=float(np.max(abs(pf-pc)))
            rows.append({"protocol":protocol["id"],"initial_signs":pair,"probabilities":pf.tolist(),
                         "original_pair_retained":float(pf[match].sum()),"positive_pair":float(pf[positive].sum()),
                         "unassigned":unassigned,"hold_switching":switching,"refinement_difference":difference,
                         "valid":bool(unassigned<=.001 and switching<=.01 and difference<=.015)})
    return {"seed":seed,"n":n,"scheme":scheme,"states":states,"results":rows,
            "all_valid":all(r["valid"] for r in rows)}


def lower(p,n):
    z=1.96
    return (p+z*z/(2*n)-z*np.sqrt((p*(1-p)+z*z/(4*n))/n))/(1+z*z/n)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--task",type=Path,required=True)
    parser.add_argument("--parent",type=Path,help="Completed timing inquiry; use its selected pulse centre")
    parser.add_argument("--out-dir",type=Path,required=True)
    parser.add_argument("--n",type=int,default=256)
    args=parser.parse_args()
    args.out_dir.mkdir(parents=True,exist_ok=False)
    source=json.loads(args.task.read_text(encoding="utf-8"))
    c,M,witness=pulse.mixing_witness(source)
    mu=.16
    if args.parent:
        parent=json.loads((args.parent/"report.json").read_text(encoding="utf-8"))
        design=json.loads((args.parent/"committed_design.json").read_text(encoding="utf-8"))
        if not parent["timing_effect_confirmed"] or st.digest(source)!=st.digest(design["tasks"][0]):
            raise SystemExit("A completed inquiry on the same source is required")
        selected=next(s for s in design["selected"] if s["id"]==design["high"])
        mu=(selected["lower"]+selected["upper"])/2
    folds=sorted(pulse.mixed_state_folds(c,M,mu),key=lambda f:f["coupling"])
    if len(folds)!=2:
        raise SystemExit("Two distinct mixed-state folds are required")
    lo,hi=folds
    peak=(lo["coupling"]+hi["coupling"])/2
    gap=hi["coupling"]-lo["coupling"]
    if gap<=1e-6:
        raise SystemExit("No resolved selective interval")
    delay=np.pi/np.sqrt(abs(lo["parameter_coefficient"]*lo["quadratic_coefficient"]*(peak-lo["coupling"])))
    duration=2*delay
    barrier=(4/3)*abs(hi["parameter_coefficient"])**1.5*(hi["coupling"]-peak)**1.5/np.sqrt(abs(hi["quadratic_coefficient"]))
    D=barrier/100
    base=c["real"].params["kap"]
    schedules=[{"id":"no_pulse","amplitude":0.,"short":False},
               {"id":"below_both","amplitude":lo["coupling"]-gap/2-base,"short":False},
               {"id":"selective","amplitude":peak-base,"short":False},
               {"id":"above_both","amplitude":hi["coupling"]+gap/2-base,"short":False},
               {"id":"short_selective","amplitude":peak-base,"short":True}]
    initial=pulse.instantaneous_states(c,M,mu,0.)["stable_states"]
    short=.25*delay
    files=[Path(__file__),Path(pulse.__file__),Path(st.__file__)]
    hashes={str(p.relative_to(ROOT)):inquiry.file_hash(p) for p in files}
    plan={"source":source,"mu_pause":mu,"folds":folds,"selected_coupling":peak,
          "fold_gap":gap,"normal_form_passage_time":delay,"pause_duration":duration,
          "short_pulse_duration":short,"protected_normal_form_barrier_estimate":barrier,
          "writing_and_holding_diffusion":D,"initial_states":initial,"schedules":schedules,
          "source_hashes":hashes,"seeds":[505000,505002],"n":args.n,
          "question":"Can the split folds select one mixed-sign state for erasure while protecting the other?",
          "noise_scope":"A new low-noise regime selected from the local barrier estimate, not the earlier writing diffusion.",
          "falsifiers":["protected state erased","target state retained","loss of assignment or convergence"]}
    save(args.out_dir/"plan.json",plan)
    exact=[]
    for protocol in schedules[:4]:
        for q in initial:
            a=deterministic(c,M,q,protocol["amplitude"],duration,mu,"DOP853")
            b=deterministic(c,M,q,protocol["amplitude"],duration,mu,"Radau")
            exact.append({"protocol":protocol["id"],"initial":q,"final":a.tolist(),
                          "solver_difference":float(np.max(abs(a-b))),"label":int(st.labels_of(a[None,:],st.final_states(c))[0])})
    save(args.out_dir/"deterministic.json",exact)
    prediction=simulate(c,M,initial,schedules,mu,duration,short,D,args.n,505000,"euler")
    save(args.out_dir/"prediction.json",prediction)
    commitment=inquiry.file_hash(args.out_dir/"prediction.json")
    save(args.out_dir/"prediction_commitment.json",{"sha256":commitment})
    print(json.dumps({"stage":"selective_erasure_prediction_committed","peak":peak,"diffusion":D,
                      "results":[{k:r[k] for k in ("protocol","initial_signs","original_pair_retained","positive_pair")} for r in prediction["results"]]}),flush=True)
    reference=simulate(c,M,initial,schedules,mu,duration,short,D,args.n,505002,"heun")
    save(args.out_dir/"reference.json",reference)
    target=next(r for r in reference["results"] if r["protocol"]=="selective" and r["initial_signs"]==[1,-1])
    protected=next(r for r in reference["results"] if r["protocol"]=="selective" and r["initial_signs"]==[-1,1])
    selective=(reference["all_valid"] and prediction["all_valid"] and
               lower(target["positive_pair"],args.n)>.95 and lower(protected["original_pair_retained"],args.n)>.95)
    verified_probabilities=[]
    for p,r in zip(prediction["results"],reference["results"]):
        check=st.compare(p["probabilities"],{"probabilities":r["probabilities"],"n":args.n,"valid":r["valid"]},args.n,st.Policy())
        verified_probabilities.append(check["passed"] and p["valid"])
    report={"complete":True,"selective_erasure_confirmed":bool(selective and all(verified_probabilities)),
            "selected_coupling":peak,"fold_gap":gap,"pulse_duration":duration,"diffusion":D,
            "target_to_positive_pair":target["positive_pair"],"protected_pair_retained":protected["original_pair_retained"],
            "all_probability_predictions_verified":all(verified_probabilities),
            "deterministic_solver_max_difference":max(r["solver_difference"] for r in exact),
            "mixing_witness":witness,"source_hashes":hashes,"novelty_established":False,
            "scope":"A constructed low-noise operating regime of the supplied equations; finite-noise behavior tested at the stated diffusion."}
    def retention(protocol,pair):
        return next(r["original_pair_retained"] for r in reference["results"] if r["protocol"]==protocol and r["initial_signs"]==pair)
    controls=(retention("below_both",[1,-1])>.95 and retention("below_both",[-1,1])>.95 and
              retention("short_selective",[1,-1])>.95 and retention("above_both",[-1,1])<.05)
    report["control_predictions_verified"]=controls
    if report["selective_erasure_confirmed"] and controls:
        save(args.out_dir/"construction_recipe.json",{
            "operation":"selectively erase a mixed-sign state using split coupling folds",
            "source_task_sha256":st.digest(source),"stage_control":mu,
            "amplitude_rule":"midpoint between two verified nondegenerate mixed-state folds",
            "duration_rule":"twice the lower-fold normal-form passage time",
            "noise_rule":"protected local barrier estimate divided by 100; verify trajectories afresh",
            "selected_coupling":peak,"pulse_duration":duration,"verified_diffusion":D,
            "preserves":["final drift","final state catalogue","holding conditions across protocols"],
            "scope":"conditional recipe verified on this source; another source requires new folds and verification",
            "controls_verified":controls,"source_hashes":hashes})
    if commitment!=inquiry.file_hash(args.out_dir/"prediction.json") or any(inquiry.file_hash(ROOT/p)!=h for p,h in hashes.items()):
        raise RuntimeError("A frozen input changed")
    save(args.out_dir/"report.json",report)
    print(json.dumps({k:v for k,v in report.items() if k not in {"mixing_witness","source_hashes"}},indent=2),flush=True)
    return 0 if reference["all_valid"] else 2


if __name__=="__main__":
    raise SystemExit(main())
