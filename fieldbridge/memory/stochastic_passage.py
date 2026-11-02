"""Test and repair a Gaussian-sector prediction for a driven polynomial diffusion.

The supplied prediction uses the Gaussian of the linearized drift to assign
angular sectors. The alternative keeps a Gaussian prefix only while the rms
nonlinear drift is small, then propagates the complete Ito generator. This is
a supplied approximation and repair backend, not a learned physical law.
Every probability refers to fixed points of one fixed final drift, after the
writing force is removed. Reuse freezes the matching and acceptance policy,
not the coefficients of a different physical system.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
from itertools import product
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import root
import sympy as sp

from . import spec as mspec

SCHEMA = "fieldbridge-stochastic-passage/1"
RULE = "gaussian-prefix-full-generator/1"


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, allow_nan=False).encode()).hexdigest()


@dataclass(frozen=True)
class Policy:
    matching_fraction: float = .05
    moment_step: float = .1
    dt: float = .02
    n_prediction: int = 4096
    n_reference: int = 4096
    probability_floor: float = .02
    sigma: float = 4.
    refinement_floor: float = .015
    max_unassigned: float = .001
    max_hold_switching: float = .01

    def __post_init__(self):
        if not (0 < self.matching_fraction <= .1 and 0 < self.moment_step <= .2 and 0 < self.dt <= .02
                and min(self.n_prediction, self.n_reference) >= 256 and 3 <= self.sigma <= 6):
            raise ValueError("Invalid fixed passage policy")


def compile_task(task):
    if task.get("schema") != SCHEMA or task.get("convention") != "ito":
        raise ValueError("Supply the passage schema and explicit Ito convention")
    real = mspec.load(task["model"])
    if task["model"]["kind"] != "equations" or real.carrier.kind != "euclid" or not 1 <= real.carrier.dim <= 4:
        raise ValueError("This backend supports one to four real polynomial coordinates")
    x, p, F = [real.symbolic[k] for k in ("variables", "parameters", "drift")]
    if not real.control:
        raise ValueError("Declare the swept control")
    degree = max(sp.Poly(f,*x).total_degree() for f in F)
    if degree > 5:
        raise ValueError("Polynomial degree exceeds the passage backend")
    params = {s: sp.Rational(str(real.params[str(s)])) for s in p if str(s) != real.control}
    F = sp.Matrix([sp.expand(f.subs(params)) for f in F])
    mu = next(s for s in p if str(s) == real.control)
    at_zero = dict.fromkeys(x, 0)
    if any(sp.simplify(f.subs(at_zero)) != 0 for f in F):
        raise ValueError("The current backend requires a fixed zero reference branch")
    J = F.jacobian(x).subs(at_zero)
    if sp.simplify(J - J.T) != sp.zeros(len(x)):
        raise ValueError("The supplied angular prediction currently requires a symmetric linearization")
    if any(sp.simplify(f) != 0 for f in J.diff(mu) - sp.eye(len(x))):
        raise ValueError("The control must sweep the linear growth rate with unit slope")
    nonlinear = sp.simplify(F - J*sp.Matrix(x))
    protocol = task["protocol"]
    start, end, rate = [float(protocol[k]) for k in ("start", "end", "rate")]
    lo, hi = real.control_range
    if not (lo <= start < end <= hi and rate > 0):
        raise ValueError("Sweep leaves the declared control range")
    D = float(task["diffusion"])
    if not np.isclose(real.noise,D,rtol=1e-12,atol=1e-15):
        raise ValueError("The writing diffusion must match the source specification")
    hold_D, hold_T, relax_T = [float(task[k]) for k in ("hold_diffusion", "hold_time", "relax_time")]
    bias = np.asarray(task["bias"], dtype=float)
    mean = np.asarray(task["preparation"]["mean"], dtype=float)
    cov = np.asarray(task["preparation"]["covariance"], dtype=float)
    n = len(x)
    if (bias.shape != (n,) or mean.shape != (n,) or cov.shape != (n,n)
            or not np.allclose(cov, cov.T) or np.linalg.eigvalsh(cov).min() < -1e-12
            or min(D, hold_D, hold_T, relax_T) < 0
            or not all(np.isfinite(a).all() for a in (bias, mean, cov))
            or not np.isfinite([D, hold_D, hold_T, relax_T, start, end, rate]).all()):
        raise ValueError("Invalid preparation, noise or protocol")
    fn = sp.lambdify([*x, mu], list(F), "numpy")
    jfn = sp.lambdify(mu, J, "numpy")
    def drift(q, control):
        q = np.atleast_2d(q)
        vals = fn(*q.T, control)
        return np.stack([np.broadcast_to(v, (len(q),)) for v in vals], axis=-1)
    return {"real": real, "variables": x, "control_symbol": mu, "F": F, "J": J,
            "nonlinear": nonlinear, "drift": drift, "jacobian": jfn, "start": start, "end": end,
            "rate": rate, "D": D, "bias": bias, "mean": mean, "covariance": cov,
            "hold_D": hold_D, "hold_T": hold_T, "relax_T": relax_T, "dimension": n,"degree":degree}


def generator_witness(compiled):
    """Coordinate and quadratic probes determine drift and diffusion of a second-order Ito generator."""
    x, F, J = compiled["variables"], compiled["F"], compiled["J"]
    D = sp.Symbol("D", nonnegative=True)
    linear = J*sp.Matrix(x)
    def apply(f, field, diffusion):
        return sp.expand(sum(sp.diff(f, xi)*fi for xi,fi in zip(x,field))
                         + diffusion*sum(sp.diff(f, xi, 2) for xi in x))
    probes = [*x, *[x[i]*x[j] for i in range(len(x)) for j in range(i,len(x))]]
    rows = []
    for probe in probes:
        full = apply(probe, F, D)
        rows.append({"observable": str(probe), "linear_defect": str(sp.expand(full-apply(probe, linear, D))),
                     "omit_diffusion_defect": str(sp.expand(full-apply(probe, F, 0))),
                     "repaired_defect": str(sp.expand(full-apply(probe, linear+(F-linear), D)))})
    return {"operator": "F.grad + D*Laplacian", "added_drift": list(map(str, F-linear)),
            "restored_diffusion": "D*Laplacian", "probes": rows,
            "exact_generator_identity": all(r["repaired_defect"] == "0" for r in rows)}


def scaling_witness(compiled):
    """Exact finite-rate scaling for a homogeneous cubic drift, without fitting a response curve."""
    x=compiled["variables"]
    mu=compiled["control_symbol"]
    rate=sp.Symbol("r",positive=True)
    tau=sp.Symbol("tau",real=True)
    v=sp.symbols(f"v0:{len(x)}",real=True)
    scaled=compiled["F"].subs({**dict(zip(x,[rate**sp.Rational(1,4)*s for s in v])),
                              mu:sp.sqrt(rate)*tau},simultaneous=True)/rate**sp.Rational(3,4)
    expected=(tau*sp.eye(len(x))+compiled["J"].subs(mu,0)/sp.sqrt(rate))*sp.Matrix(v)
    expected+=compiled["nonlinear"].subs(dict(zip(x,v)),simultaneous=True)
    residual=list(map(sp.simplify,scaled-expected))
    return {"state_scaling":"x=r^(1/4)*v","time_scaling":"tau=sqrt(r)*t",
            "scaled_drift":list(map(str,expected)),"scaled_bias":"H=h/r^(3/4)",
            "scaled_diffusion":"Gamma=D/r","residuals":list(map(str,residual)),
            "homogeneous_cubic_scaling_verified":all(r==0 for r in residual)}


def final_states(compiled):
    """Find and check final equilibria; unmatched trajectories remain unassigned rather than entering forced bins."""
    F = compiled["F"].subs(compiled["control_symbol"], compiled["end"])
    x = compiled["variables"]
    jac = F.jacobian(x)
    jf = sp.lambdify([x], jac, "numpy")
    f = lambda q: compiled["drift"](q, compiled["end"])[0]
    scale = max(1., np.sqrt(abs(compiled["end"])))
    states, spectra = [], []
    for q in product((-scale, 0., scale), repeat=len(x)):
        sol = root(f, q, jac=jf, tol=1e-11)
        if np.linalg.norm(f(sol.x)) > 1e-8:
            continue
        ev = np.linalg.eigvalsh(np.asarray(jf(sol.x),float))
        if ev.max() >= -1e-7 or any(np.linalg.norm(sol.x-a)<1e-5 for a in states):
            continue
        states.append(sol.x); spectra.append(ev.tolist())
    order = sorted(range(len(states)), key=lambda i:tuple(np.round(states[i],8)))
    states = np.array([states[i] for i in order])
    if len(states) < 2:
        raise ValueError("No checked multistable final state set")
    distances = np.linalg.norm(states[:,None]-states[None,:],axis=-1)
    distances[distances == 0] = np.inf
    radius = min(.15, .24*distances.min())
    return {"states": states.tolist(), "eigenvalues": [spectra[i] for i in order], "assignment_radius": radius,
            "scope": "stable equilibria found from the declared grid; unmatched samples fail the check"}


def gaussian_points(mean, covariance, count, seed):
    ev, vec = np.linalg.eigh((covariance+covariance.T)/2)
    return mean + np.random.default_rng(seed).standard_normal((count,len(mean))) @ (vec*np.sqrt(np.maximum(ev,0))).T


def linear_prediction(compiled, states, policy):
    n, r, D = compiled["dimension"], compiled["rate"], compiled["D"]
    start, end = compiled["start"], compiled["end"]
    crossing = -float(np.linalg.eigvalsh(compiled["jacobian"](0)).max())
    t_cross = (crossing-start)/r
    if not 0 < t_cross < (end-start)/r:
        raise ValueError("The declared sweep does not cross the zero-branch instability")
    t_sector = min((end-start)/r, t_cross+np.sqrt(2/r))
    initial = np.r_[compiled["mean"],compiled["covariance"].ravel()]
    def rhs(t,z):
        J = np.asarray(compiled["jacobian"](start+r*t),float)
        m,C = z[:n],z[n:].reshape(n,n)
        return np.r_[J@m+compiled["bias"],(J@C+C@J.T+2*D*np.eye(n)).ravel()]
    sol = solve_ivp(rhs,(0,t_sector),initial,rtol=1e-10,atol=1e-12,dense_output=True,max_step=.2)
    if not sol.success:
        raise ArithmeticError(sol.message)
    # Tensor Gaussian quadrature integrates the squared polynomial remainder exactly.
    order = compiled["degree"]+1
    nodes,weights = np.polynomial.hermite.hermgauss(order)
    grid = np.array(list(product(range(order),repeat=n)))
    standard = np.sqrt(2)*nodes[grid]
    weights = np.prod(weights[grid]/np.sqrt(np.pi),axis=1)
    def ratio(t):
        z = sol.sol(t); m,C = z[:n],z[n:].reshape(n,n)
        ev,V = np.linalg.eigh((C+C.T)/2)
        q = m + standard @ (V*np.sqrt(np.maximum(ev,0))).T
        J = np.asarray(compiled["jacobian"](start+r*t),float)
        N = compiled["drift"](q,start+r*t)-q@J.T
        numerator = np.sqrt(np.sum(weights*np.sum(N*N,axis=1)))
        denominator = np.sqrt(r)*np.sqrt(np.sum(weights*np.sum(q*q,axis=1)))
        return numerator/max(denominator,1e-30)
    times = np.linspace(0,t_sector,int(np.ceil(t_sector/policy.moment_step))+1)
    match = next((float(t) for t in times if ratio(t)>=policy.matching_fraction),float(t_sector))
    if match == 0:
        raise ValueError("Preparation is already outside the Gaussian-prefix applicability condition")
    zmatch,zsector = sol.sol(match),sol.sol(t_sector)
    samples = gaussian_points(zsector[:n],zsector[n:].reshape(n,n),65536,1729)
    S = np.asarray(states["states"])
    units = S/np.linalg.norm(S,axis=1)[:,None]
    labels = np.argmax(samples@units.T,axis=1)
    return {"crossing": crossing, "t_cross": t_cross, "t_sector": t_sector,
            "sector_probabilities": probabilities(labels,len(S)), "sector_samples": len(samples),
            "matching_time": match, "matching_control": start+r*match, "matching_ratio": ratio(match),
            "matching_mean": zmatch[:n].tolist(), "matching_covariance": zmatch[n:].reshape(n,n).tolist(),
            "matching_policy": asdict(policy), "witness": generator_witness(compiled)}


def probabilities(labels,n):
    return (np.bincount(labels[labels>=0],minlength=n)/len(labels)).tolist()


def labels_of(q,states):
    S = np.asarray(states["states"])
    distance = np.linalg.norm(q[:,None]-S[None,:],axis=-1)
    labels = np.argmin(distance,axis=1)
    labels[np.min(distance,axis=1)>states["assignment_radius"]] = -1
    labels[~np.isfinite(q).all(axis=1)] = -1
    return labels


def torch_drift(compiled,device):
    import torch
    x,mu = compiled["variables"],compiled["control_symbol"]
    function = sp.lambdify([*x,mu],list(compiled["F"]),modules="math")
    def field(q,control):
        vals = function(*[q[:,i] for i in range(q.shape[1])],control)
        return torch.stack([v+torch.zeros_like(q[:,0]) for v in vals],dim=1)
    example = torch.zeros((16,len(x)),device=device,dtype=torch.float64)
    traced = torch.jit.trace(field,(example,torch.tensor(0.,device=device,dtype=torch.float64)))
    q = np.random.default_rng(109).normal(scale=.3,size=(16,len(x)))
    for muval in (compiled["start"],compiled["end"]):
        actual = traced(torch.tensor(q,device=device),torch.tensor(muval,device=device,dtype=torch.float64)).cpu().numpy()
        if not np.allclose(actual,compiled["drift"](q,muval),atol=1e-12,rtol=1e-12):
            raise ArithmeticError("Torch equations disagree with independently compiled source equations")
    return traced


def simulate(compiled,states,policy,seed,*,device="cpu",prediction=None,omit_diffusion=False):
    """Paired Brownian refinement: reference uses stochastic Heun, repaired prediction Euler-Maruyama."""
    import torch
    field = torch_drift(compiled,device)
    n = policy.n_reference if prediction is None else policy.n_prediction
    t0 = 0. if prediction is None else prediction["matching_time"]
    mean = compiled["mean"] if prediction is None else np.asarray(prediction["matching_mean"])
    cov = compiled["covariance"] if prediction is None else np.asarray(prediction["matching_covariance"])
    q = torch.tensor(gaussian_points(mean,cov,n,seed),device=device,dtype=torch.float64)
    fine,coarse = q.clone(),q.clone()
    gen = torch.Generator(device=device).manual_seed(seed+1000003)
    total = (compiled["end"]-compiled["start"])/compiled["rate"]
    bias = torch.tensor(compiled["bias"],device=device)
    control = lambda value:torch.tensor(value,device=device,dtype=torch.float64)
    def step(q,t,dt,dw,phase):
        mu = compiled["start"]+compiled["rate"]*t if phase=="write" else compiled["end"]
        f = field(q,control(mu))+(bias if phase=="write" else 0.)
        if prediction is not None:
            return q+dt*f+dw
        proposed = q+dt*f+dw
        mu1 = compiled["start"]+compiled["rate"]*(t+dt) if phase=="write" else compiled["end"]
        return q+.5*dt*(f+field(proposed,control(mu1))+(bias if phase=="write" else 0.))+dw
    samples = {}
    with torch.no_grad():
        for phase,duration,D,offset in (("write",total-t0,0. if omit_diffusion else compiled["D"],t0),
                                       ("relax",compiled["relax_T"],0.,0.),
                                       ("hold",compiled["hold_T"],compiled["hold_D"],0.)):
            steps = max(1,int(np.ceil(duration/policy.dt)))
            dt = duration/steps
            for k in range(steps):
                dw1 = torch.randn(q.shape,generator=gen,device=device,dtype=q.dtype)*np.sqrt(D*dt)
                dw2 = torch.randn(q.shape,generator=gen,device=device,dtype=q.dtype)*np.sqrt(D*dt)
                t = offset+k*dt
                fine = step(fine,t,dt/2,dw1,phase)
                fine = step(fine,t+dt/2,dt/2,dw2,phase)
                coarse = step(coarse,t,dt,dw1+dw2,phase)
            samples[phase] = (fine.cpu().numpy(),coarse.cpu().numpy())
    labels = {phase:(labels_of(a,states),labels_of(b,states)) for phase,(a,b) in samples.items()}
    pf = np.asarray(probabilities(labels["hold"][0],len(states["states"])))
    pc = np.asarray(probabilities(labels["hold"][1],len(states["states"])))
    pre,post = labels["relax"][0],labels["hold"][0]
    assigned = (pre>=0)&(post>=0)
    switching = float(np.mean(pre[assigned]!=post[assigned])) if assigned.any() else 1.
    unassigned = float(max(np.mean(pre<0),np.mean(post<0)))
    finite = bool(all(np.isfinite(a).all() and np.isfinite(b).all() for a,b in samples.values()))
    convergence = float(np.max(np.abs(pf-pc)))
    return {"seed": seed, "n": n, "scheme": "Heun" if prediction is None else "Euler-Maruyama",
            "probabilities": pf.tolist(), "coarse_probabilities": pc.tolist(),
            "before_hold": probabilities(pre,len(states["states"])), "hold_switching": switching,
            "unassigned": unassigned, "finite": finite, "refinement_error": convergence,
            "refinement_label_disagreement": float(np.mean(labels["hold"][0]!=labels["hold"][1])),
            "valid": finite and unassigned<=policy.max_unassigned and convergence<=policy.refinement_floor
                     and switching<=policy.max_hold_switching}


def compare(predicted,reference,n_prediction,policy):
    p,q = np.asarray(predicted),np.asarray(reference["probabilities"])
    error = np.abs(p-q)
    se = np.sqrt(p*(1-p)/n_prediction+q*(1-q)/reference["n"])
    bound = np.maximum(policy.probability_floor,policy.sigma*se)
    return {"passed": bool(reference["valid"] and np.all(error<=bound)),
            "max_error": float(error.max()), "total_variation": float(error.sum()/2),
            "bounds": bound.tolist(), "reference_valid": reference["valid"]}


def prepare_prediction(task,policy,seed,output,device="cpu",frozen_rule=None):
    output = Path(output);output.mkdir(parents=True,exist_ok=False)
    rule_hash = None
    if frozen_rule is not None:
        raw = Path(frozen_rule).read_bytes()
        rule = json.loads(raw)
        if (rule.get("schema") != RULE or rule.get("policy") != asdict(policy)
                or rule.get("repair") != "retain the calculated nonlinear drift and diffusion after Gaussian-prefix matching"):
            raise ValueError("The frozen rule does not authorize this calculation policy")
        expected = rule["source_hashes"].get("fieldbridge/memory/stochastic_passage.py")
        if hashlib.sha256(Path(__file__).read_bytes()).hexdigest() != expected:
            raise ValueError("The frozen rule belongs to another implementation")
        rule_hash = hashlib.sha256(raw).hexdigest()
    compiled = compile_task(task)
    states = final_states(compiled)
    prediction = linear_prediction(compiled,states,policy)
    repaired = simulate(compiled,states,policy,seed,device=device,prediction=prediction)
    omitted = simulate(compiled,states,policy,seed+1,device=device,prediction=prediction,omit_diffusion=True)
    commitment = {"task":task,"task_sha256":digest(task),"states":states,"linear":prediction,
                  "scaling":scaling_witness(compiled),
                  "repaired_prediction":repaired,"omit_diffusion":omitted,"policy":asdict(policy),
                  "frozen_rule_sha256":rule_hash,"frozen_rule_path":str(frozen_rule) if frozen_rule else None,
                  "implementation_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    # No reference outcomes have been evaluated when the prediction is committed.
    path = output/"prediction.json"
    path.write_text(json.dumps(commitment,indent=2,allow_nan=False)+"\n", encoding="utf-8")
    return {"prediction_sha256":hashlib.sha256(path.read_bytes()).hexdigest(),
            "prediction_valid":repaired["valid"]}


def verify_prediction(output,seed,device="cpu"):
    output=Path(output)
    if (output/"comparison.json").exists():
        raise ValueError("Reference evaluation has already been used")
    path=output/"prediction.json"
    prediction_hash=hashlib.sha256(path.read_bytes()).hexdigest()
    commitment=json.loads(path.read_text(encoding="utf-8"))
    if hashlib.sha256(Path(__file__).read_bytes()).hexdigest()!=commitment["implementation_sha256"]:
        raise ValueError("Implementation changed after prediction")
    task,states,prediction=[commitment[k] for k in ("task","states","linear")]
    repaired,omitted=[commitment[k] for k in ("repaired_prediction","omit_diffusion")]
    policy=Policy(**commitment["policy"])
    compiled=compile_task(task)
    reference=simulate(compiled,states,policy,seed,device=device)
    result = {"task_sha256":digest(task),"prediction_sha256":prediction_hash,"reference":reference,
              "frozen_rule_sha256":commitment["frozen_rule_sha256"],
              "linear_test":compare(prediction["sector_probabilities"],reference,prediction["sector_samples"],policy),
              "repair_test":compare(repaired["probabilities"],reference,repaired["n"],policy),
              "diffusion_omission_test":compare(omitted["probabilities"],reference,omitted["n"],policy),
              "prediction_valid":repaired["valid"],"generator_identity":prediction["witness"]["exact_generator_identity"]}
    result["repair_verified"] = result["repair_test"]["passed"] and result["prediction_valid"] and result["generator_identity"]
    if (commitment["frozen_rule_path"] is not None and
            hashlib.sha256(Path(commitment["frozen_rule_path"]).read_bytes()).hexdigest()!=commitment["frozen_rule_sha256"]):
        raise RuntimeError("The rule changed during its application")
    if hashlib.sha256(path.read_bytes()).hexdigest()!=prediction_hash:
        raise RuntimeError("The prediction changed during reference evaluation")
    (output/"comparison.json").write_text(json.dumps(result,indent=2,allow_nan=False)+"\n", encoding="utf-8")
    return result


def evaluate(task,policy,seed,output,device="cpu",frozen_rule=None):
    prepare_prediction(task,policy,seed,output,device,frozen_rule)
    return verify_prediction(output,seed+2,device)
