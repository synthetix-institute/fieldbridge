"""Calculate equal-dose mixing pulses with the same initial and final equations."""
from __future__ import annotations

import copy
from itertools import product
import math

import numpy as np
from scipy.optimize import root
import sympy as sp

from . import stochastic_passage as st
from . import paired_occupations


def mixing_witness(task, parameter="kap"):
    c = st.compile_task(task)
    real = c["real"]
    variables = c["variables"]
    original = sp.Matrix(real.symbolic["drift"])
    parameters = real.symbolic["parameters"]
    symbol = next((p for p in parameters if str(p)==parameter), None)
    if symbol is None:
        raise ValueError("The source equations do not contain the requested coupling")
    B = original.diff(symbol)
    if B.diff(symbol)!=sp.zeros(len(variables),1):
        raise ValueError("The pulsed parameter must enter affinely")
    values = {p:sp.Rational(str(real.params[str(p)])) for p in parameters if str(p)!=real.control}
    B = sp.simplify(B.subs(values))
    M = B.jacobian(variables)
    if any(sp.simplify(e)!=0 for e in B-M*sp.Matrix(variables)) or any(e.free_symbols for e in M):
        raise ValueError("The mixing contribution must be a fixed linear map")
    n = len(variables)
    if sp.ones(1,n)*M != sp.zeros(1,n):
        raise ValueError("The proposed pulse has a direct force on the mean")
    if M!=M.T or max(np.linalg.eigvalsh(np.asarray(M,float)))>1e-12:
        raise ValueError("The proposed mixing must be symmetric and dissipative")
    q = sp.Matrix(variables)
    mean = sum(variables)/n
    V = sum((v-mean)**2 for v in variables)/2
    variance_action = sp.expand(sum(sp.diff(V,x)*b for x,b in zip(variables,B)))
    F = c["F"]
    bracket = sp.simplify(B.jacobian(variables)*F-F.jacobian(variables)*B)
    projected = sp.simplify(sum(bracket)/n)
    return c, np.asarray(M,float), {
        "parameter":parameter, "mixing_field":list(map(str,B)), "matrix":np.asarray(M,float).tolist(),
        "direct_mean_force":str(sp.simplify(sum(B)/n)),
        "internal_variance":str(sp.expand(V)), "variance_action_per_unit_coupling":str(variance_action),
        "complete_graph_variance_identity":sp.simplify(variance_action+2*n*V)==0,
        "bracket_convention":"DB.F - DF.B", "bracket_on_mean":str(projected),
        "noncommuting_mean_dynamics":projected!=0,
        "final_generator_unchanged":"The pulse is zero before the sweep ends and throughout relaxation and holding.",
    }


def instantaneous_states(c, M, mu, amplitude):
    x = c["variables"]
    F = c["F"].subs(c["control_symbol"],sp.Rational(str(mu)))
    F += sp.Rational(str(amplitude))*sp.Matrix(M)*sp.Matrix(x)+sp.Matrix(c["bias"])
    jacobian = sp.lambdify([x],F.jacobian(x),"numpy")
    drift = sp.lambdify([x],F,"numpy")
    states=[]
    for start in product((-1.,0.,1.), repeat=c["dimension"]):
        solution=root(lambda q:np.asarray(drift(q),float).ravel(),start,jac=jacobian,tol=1e-11)
        if np.linalg.norm(drift(solution.x))>1e-8:
            continue
        if np.linalg.eigvalsh(jacobian(solution.x)).max()>=-1e-7:
            continue
        if not any(np.linalg.norm(solution.x-q)<1e-5 for q in states):
            states.append(solution.x)
    return {"stable_states":[q.tolist() for q in states],"count":len(states),
            "scope":"stable equilibria found from the declared grid; not an exhaustive root proof"}


def mixed_state_folds(c, M, mu):
    """Locate nondegenerate folds of mixed-sign states as coupling increases."""
    x=c["variables"]
    a=sp.Symbol("excess_coupling",real=True)
    field=c["F"].subs(c["control_symbol"],sp.Rational(str(mu)))+a*sp.Matrix(M)*sp.Matrix(x)+sp.Matrix(c["bias"])
    J=field.jacobian(x)
    equations=sp.Matrix([*field,J.det()])
    arguments=[*x,a]
    fn=sp.lambdify([arguments],equations,"numpy")
    jf=sp.lambdify([arguments],equations.jacobian(arguments),"numpy")
    jstate=sp.lambdify([arguments],J,"numpy")
    seeds=instantaneous_states(c,M,mu,0.)["stable_states"]
    folds=[]
    for state in seeds:
        if min(state)>=0 or max(state)<=0:
            continue
        for amplitude in (.005,.015,.025,.05):
            solution=root(lambda q:np.asarray(fn(q),float).ravel(),[*state,amplitude],jac=jf,tol=1e-11)
            if np.linalg.norm(fn(solution.x))>1e-8 or not 0<solution.x[-1]<.1:
                continue
            q,amp=solution.x[:-1],solution.x[-1]
            eigenvalues,eigenvectors=np.linalg.eigh(jstate(solution.x))
            if abs(eigenvalues[-1])>1e-7 or eigenvalues[-2]>=-1e-6:
                continue
            if any(abs(amp-r["excess_coupling"])<1e-6 and np.linalg.norm(q-r["state"])<1e-5 for r in folds):
                continue
            v=eigenvectors[:,-1]
            substitutions=dict(zip(arguments,solution.x))
            b=float(v @ (M@q))
            second=np.array([float((sp.Matrix(v).T*sp.hessian(f,x)*sp.Matrix(v))[0].subs(substitutions)) for f in field])
            quadratic=float(.5*v@second)
            if abs(b)<1e-8 or abs(quadratic)<1e-8:
                continue
            folds.append({"mu":mu,"state":q.tolist(),"excess_coupling":float(amp),
                          "coupling":float(c["real"].params["kap"]+amp),
                          "equation_residual":float(np.linalg.norm(fn(solution.x))),
                          "jacobian_eigenvalues":eigenvalues.tolist(),"parameter_coefficient":b,
                          "quadratic_coefficient":quadratic,
                          "equilibria_on_lower_coupling_side":bool(b/quadratic>0)})
    return folds


def run(task, schedules, *, n=1024, seed=100, dt=.02, device="cpu", scheme="heun"):
    """Shared Brownian trajectories for paired pulse comparisons and step refinement."""
    import torch
    if scheme not in ("heun","euler") or n<256 or not 0<dt<=.02:
        raise ValueError("Invalid pulse calculation policy")
    c,M,witness=mixing_witness(task)
    states=st.final_states(c)
    k=len(schedules)
    if not k or len({s["id"] for s in schedules})!=k:
        raise ValueError("Supply distinct pulse identifiers")
    start,end,rate=c["start"],c["end"],c["rate"]
    for s in schedules:
        if not start<s["lower"]<s["upper"]<end or not np.isfinite(s["amplitude"]) or s["amplitude"]<0:
            raise ValueError("A pulse must end strictly before the final control and have nonnegative finite strength")
    boundaries=sorted({start,end,*[v for s in schedules for v in (s["lower"],s["upper"])]})
    field=st.torch_drift(c,device)
    matrix=torch.tensor(M,device=device,dtype=torch.float64)
    initial=st.gaussian_points(c["mean"],c["covariance"],n,seed)
    q=torch.tensor(np.broadcast_to(initial,(k,*initial.shape)).copy(),device=device,dtype=torch.float64)
    fine,coarse=q.clone(),q.clone()
    generator=torch.Generator(device=device).manual_seed(seed+1000003)
    bias=torch.tensor(c["bias"],device=device,dtype=torch.float64)
    def rhs(q,mu,amplitudes,writing):
        value=field(q.reshape(-1,c["dimension"]),torch.tensor(mu,device=device,dtype=torch.float64)).reshape(q.shape)
        return value+amplitudes*(q@matrix.T)+(bias if writing else 0.)
    def step(q,mu,delta,noise,amplitudes,writing):
        first=rhs(q,mu,amplitudes,writing)
        if scheme=="euler":
            return q+delta*first+noise
        proposed=q+delta*first+noise
        last=rhs(proposed,mu+rate*delta if writing else mu,amplitudes,writing)
        return q+.5*delta*(first+last)+noise
    moments=[]
    segments=[("write",(b-a)/rate,a,b) for a,b in zip(boundaries[:-1],boundaries[1:])]
    segments += [("relax",c["relax_T"],end,end),("hold",c["hold_T"],end,end)]
    before_hold=None
    with torch.no_grad():
        for phase,duration,lower,upper in segments:
            writing=phase=="write"
            middle=(lower+upper)/2
            amps=[s["amplitude"] if writing and s["lower"]<=middle<s["upper"] else 0. for s in schedules]
            amplitude=torch.tensor(amps,device=device,dtype=torch.float64)[:,None,None]
            D=c["D"] if writing else (c["hold_D"] if phase=="hold" else 0.)
            count=max(1,int(math.ceil(duration/dt)))
            delta=duration/count
            for j in range(count):
                noise1=torch.randn((1,n,c["dimension"]),generator=generator,device=device,dtype=q.dtype)*np.sqrt(D*delta)
                noise2=torch.randn((1,n,c["dimension"]),generator=generator,device=device,dtype=q.dtype)*np.sqrt(D*delta)
                mu=lower+rate*j*delta if writing else end
                fine=step(fine,mu,delta/2,noise1,amplitude,writing)
                fine=step(fine,mu+rate*delta/2 if writing else mu,delta/2,noise2,amplitude,writing)
                coarse=step(coarse,mu,delta,noise1+noise2,amplitude,writing)
            if writing:
                mean=fine.mean(-1)
                V=((fine-mean[...,None])**2).sum(-1)/2
                moments.append({"mu":upper,"mean":mean.mean(-1).cpu().tolist(),
                                "internal_variance":V.mean(-1).cpu().tolist()})
            if phase=="relax":
                before_hold=fine.cpu().numpy().copy()
    final=fine.cpu().numpy()
    coarse=coarse.cpu().numpy()
    labels=np.array([st.labels_of(q,states) for q in final])
    coarse_labels=np.array([st.labels_of(q,states) for q in coarse])
    old_labels=np.array([st.labels_of(q,states) for q in before_hold])
    rows=[]
    for i,s in enumerate(schedules):
        p=np.asarray(st.probabilities(labels[i],len(states["states"])))
        pc=np.asarray(st.probabilities(coarse_labels[i],len(states["states"])))
        unassigned=float(max(np.mean(labels[i]<0),np.mean(old_labels[i]<0)))
        switching=float(np.mean(old_labels[i]!=labels[i]))
        refinement=float(np.max(np.abs(p-pc)))
        finite=bool(np.isfinite(final[i]).all() and np.isfinite(coarse[i]).all())
        rows.append({"id":s["id"],"pulse":s,"probabilities":p.tolist(),"coarse_probabilities":pc.tolist(),
                     "unassigned":unassigned,"hold_switching":switching,"step_refinement_difference":refinement,
                     "valid":bool(finite and unassigned<=.001 and switching<=.01 and refinement<=.015),
                     "integrated_excess_coupling":s["amplitude"]*(s["upper"]-s["lower"])/rate})
    index=next(i for i,q in enumerate(states["states"]) if min(q)>0)
    contrasts=[];all_state_contrasts=[]
    for i in range(k):
        for j in range(i):
            for state_index in range(len(states["states"])):
                differences=(labels[i]==state_index).astype(float)-(labels[j]==state_index)
                effect=float(differences.mean())
                se=float(differences.std(ddof=1)/np.sqrt(n))
                row={"left":schedules[i]["id"],"right":schedules[j]["id"],"state_index":state_index,
                     **paired_occupations.contrast(int(np.sum(differences>0)),int(np.sum(differences<0)),n)}
                all_state_contrasts.append(row)
                if state_index==index:
                    contrasts.append(row)
    return {"schema":"fieldbridge-pulse-inquiry/1","task":copy.deepcopy(task),"seed":seed,"n":n,"dt":dt,
            "scheme":scheme,"states":states,"positive_state_index":index,"mixing_witness":witness,
            "results":rows,"contrasts":contrasts,"all_state_contrasts":all_state_contrasts,"moments":moments,
            "all_valid":all(r["valid"] for r in rows),"novelty_established":False}
