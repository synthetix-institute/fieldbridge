/* Small, deterministic consequences of the displayed models. No retrieval or inference. */
(function (root, factory) {
  const api = factory();
  if (typeof module === 'object' && module.exports) module.exports = api;
  else root.FieldBridgePhysics = api;
})(typeof globalThis !== 'undefined' ? globalThis : this, function () {
  'use strict';
  function roots(a, b, h) {
    if (b === 0) return a !== 0 ? [-h / a] : [];
    const p = -a / b, q = -h / b;
    const d = q*q/4 + p*p*p/27;
    if (d > 1e-12) return [Math.cbrt(-q/2 + Math.sqrt(d)) + Math.cbrt(-q/2 - Math.sqrt(d))];
    if (Math.abs(d) <= 1e-12) {
      const u = Math.cbrt(-q/2);
      return [...new Set([2*u, -u])].sort((x,y) => x-y);
    }
    const radius = 2*Math.sqrt(-p/3);
    const angle = Math.acos(Math.max(-1, Math.min(1, (3*q/(2*p))*Math.sqrt(-3/p))))/3;
    return [0,1,2].map(k => radius*Math.cos(angle-2*Math.PI*k/3)).sort((x,y) => x-y);
  }
  function memory(s) {
    const a = s.feedback ? s.eps : 0;
    const b = s.saturation ? s.gamma : 0;
    const h = s.field ? s.h : 0;
    const equilibria = roots(a,b,h).map(x => ({x, slope:a-3*b*x*x,
      stable:a-3*b*x*x < -1e-9 || (Math.abs(x)<1e-9 && a===0 && b>0)}));
    const stable = equilibria.filter(e => e.stable);
    const bistable = stable.length === 2;
    const barrier = bistable && h === 0 ? a*a/(4*b) : null;
    const threshold = a > 0 && b > 0 ? 2*Math.pow(a,1.5)/(3*Math.sqrt(3*b)) : null;
    const tau = s.thermal && barrier !== null && s.D > 0 && s.D/barrier <= .2
      ? 2*Math.PI/(Math.sqrt(2)*a)*Math.exp(barrier/s.D) : null;
    return {a,b,h,equilibria,stable,bistable,barrier,threshold,tau,
      bounded: b > 0 || a < 0, neutral: a === 0 && b === 0 && h === 0};
  }
  function trajectory(m, x0, pulse, duration=8, steps=800, window={start:1,end:3}) {
    const dt=duration/steps, data=[[0,x0]];
    let x=x0;
    for(let i=0;i<steps;i++) {
      const start=i*dt,end=(i+1)*dt;
      // Split exactly at pulse boundaries; RK4 must not average across a field jump.
      const cuts=[start,...[window.start,window.end].filter(t=>t>start&&t<end),end];
      for(let j=0;j<cuts.length-1;j++){
        const left=cuts[j],right=cuts[j+1],step=right-left,mid=(left+right)/2;
        const field=mid>=window.start&&mid<window.end?pulse:m.h;
        const f=q=>m.a*q-m.b*q*q*q+field;
        const k1=f(x),k2=f(x+step*k1/2),k3=f(x+step*k2/2),k4=f(x+step*k3);
        x+=step*(k1+2*k2+2*k3+k4)/6;
      }
      if(!Number.isFinite(x) || Math.abs(x)>20) return {data,diverged:true};
      data.push([(i+1)*dt,x]);
    }
    return {data,diverged:false};
  }
  function stochastic(s) {
    const noise=s.noise ? s.sigma : 0;
    if(s.map==='square') {
      // Unit-noise radial model; noise detach changes the model, not the map.
      const correction=noise*noise;
      const exact=2*s.theta+1+correction;
      const candidate=2*s.theta+1+(s.correction ? correction : 0);
      const rate=2*s.alpha;
      const mean=(t,correct=true)=>{
        const c=(correct ? exact : candidate)/rate;
        return c+(1-c)*Math.exp(-rate*t);
      };
      return {exact,candidate,rate,correction,residual:exact-candidate,
        variance:4*noise*noise,mean,limit:exact/rate,candidateLimit:candidate/rate};
    }
    const correction=s.convention==='ito' ? -noise*noise/2 : 0;
    const exact=s.mu+correction;
    const candidate=s.mu+(s.correction ? correction : 0);
    return {exact,candidate,correction,residual:exact-candidate,variance:noise*noise,
      mean:(t,correct=true)=>(correct ? exact : candidate)*t,limit:null,candidateLimit:null};
  }
  function constructionSpec(s) {
    if(s.map==='square') return {
      schema:'fieldbridge-construction/1',kind:'ito_transfer',question:'Which drift preserves Y=X^2?',
      parameters:{theta:'positive',alpha:'positive',sigma:'positive'},domain:'positive',target_domain:'positive',
      convention:'ito',drift:'(theta + 1/2)/x - alpha*x',noise:s.noise?'sigma':'0',state_map:'x**2',inverse_map:'sqrt(y)',
      candidate:{drift:s.correction&&s.noise?'2*theta + 1 + sigma**2 - 2*alpha*y':'2*theta + 1 - 2*alpha*y',variance:s.noise?'4*sigma**2*y':'0'},
      assumptions:['Local identity on C2 observables in x>0; boundary conditions are separate.','Mean evolution requires finite moments and no boundary contribution.'],
      provenance:{origin:'interactive worked example',novelty:'known stochastic transformation'},
      numeric_parameters:{theta:s.theta,alpha:s.alpha,sigma:s.sigma}
    };
    return {schema:'fieldbridge-construction/1',kind:'ito_transfer',question:'What is the mean growth of log X?',
      parameters:{mu:'real',sigma:'positive'},domain:'positive',target_domain:'real',convention:s.convention,
      drift:'mu*x',noise:s.noise?'sigma*x':'0',state_map:'log(x)',inverse_map:'exp(y)',
      candidate:{drift:s.correction&&s.noise&&s.convention==='ito'?'mu-sigma**2/2':'mu',variance:s.noise?'sigma**2':'0'},
      assumptions:['X(0)>0; constant coefficients; standard Brownian noise.'],
      provenance:{origin:'interactive worked example',novelty:'known stochastic transformation'},
      numeric_parameters:{mu:s.mu,sigma:s.sigma}};
  }
  function spin(s) {
    const g=s.coupling?s.g:0,h=s.field?s.h:0;
    const rate=2*Math.hypot(g,h),offset=rate?4*h*h/(rate*rate):1;
    return {g,h,rate,offset,amplitude:1-offset,
      theta:rate?Math.atan2(Math.abs(g),Math.abs(h))*180/Math.PI:0,
      signal:t=>offset+(1-offset)*Math.cos(rate*t)};
  }
  function spinSpec(s,carrier='correlated-pair') {
    const r=spin(s),drive=2*Math.abs(r.g),detuning=2*Math.abs(r.h),size=4;
    const base={schema:'fieldbridge-quantum/1',name:'Browser rotation on '+carrier,
      question:'Which Hamiltonian on this carrier preserves the detached rotation?',
      field:'worked quantum construction',assumptions:['Closed unitary evolution; hbar=1.','Prepare a top eigenstate of the observable.'],
      provenance:{origin:'interactive mechanism attachment',novelty:'known representation of su(2)'}};
    const term=(coefficient,operator)=>({coefficient,operator});
    if(carrier==='correlated-pair')return {...base,carrier:{kind:'qubits',n:2},hamiltonian:[term(r.g,'Z0 Z1'),term(r.h,'X0')],observable:[term(1,'X0')]};
    if(carrier==='qubit')return {...base,carrier:{kind:'qubits',n:1},hamiltonian:[term(drive/2,'X0'),term(detuning/2,'Z0')],observable:[term(1,'Z0')]};
    if(carrier==='spin')return {...base,carrier:{kind:'spin',j:'3/2'},hamiltonian:[term(drive,'Jx'),term(detuning,'Jz')],observable:[term(2,'Jz')]};
    if(carrier==='bosons')return {...base,carrier:{kind:'bosons',modes:['a','b'],max_quanta:4},sector:{operator:[term(1,'na + nb')],value:4},hamiltonian:[{...term(drive/2,'ad b'),hc:true},term(detuning/2,'na'),term(-detuning/2,'nb')],observable:[term(1,'na'),term(-1,'nb')]};
    if(carrier==='fermion-pair')return {...base,carrier:{kind:'fermions',modes:2},hamiltonian:[term(detuning/2,'n0'),term(detuning/2,'n1'),{...term(drive/2,'cd0 cd1'),hc:true}],observable:[term(1,'n0'),term(1,'n1')]};
    if(carrier==='collective')return {...base,carrier:{kind:'qubits',n:3},hamiltonian:[term(drive/2,'X0 + X1 + X2'),term(detuning/2,'Z0 + Z1 + Z2')],observable:[term(1,'Z0 + Z1 + Z2')]};
    if(carrier!=='chain')throw new Error('Unknown spin carrier');
    const ham=[],obs=[];
    for(let j=0;j<size-1;j++)ham.push(term(drive*Math.sqrt((j+1)*(size-1-j))/4,`X${j} X${j+1} + Y${j} Y${j+1}`));
    for(let j=0;j<size;j++){ham.push(term(-detuning*(j-1.5)/2,`Z${j}`));obs.push(term(-(j-1.5),`Z${j}`));}
    return {...base,carrier:{kind:'qubits',n:4},sector:{operator:Array.from({length:4},(_,j)=>term(1,`Z${j}`)),value:2},hamiltonian:ham,observable:obs};
  }
  return {roots,memory,trajectory,stochastic,constructionSpec,spin,spinSpec};
});
