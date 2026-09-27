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
  function trajectory(m, x0, pulse, duration=8, steps=800) {
    const dt=duration/steps, data=[[0,x0]];
    let x=x0;
    const f=(t,q) => m.a*q-m.b*q*q*q + (t >= 1 && t < 3 ? pulse : m.h);
    for(let i=0;i<steps;i++) {
      const t=i*dt, k1=f(t,x), k2=f(t+dt/2,x+dt*k1/2), k3=f(t+dt/2,x+dt*k2/2), k4=f(t+dt,x+dt*k3);
      x += dt*(k1+2*k2+2*k3+k4)/6;
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
  return {roots,memory,trajectory,stochastic,constructionSpec};
});
