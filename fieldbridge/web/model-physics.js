(function(root){
  'use strict';
  const functions={sqrt:Math.sqrt,exp:Math.exp,log:Math.log,sin:Math.sin,cos:Math.cos,tanh:Math.tanh,abs:Math.abs};
  function value(tree,env){
    const [op,a,b]=tree;
    if(op==='number')return a;
    if(op==='name')return a==='pi'?Math.PI:env[a];
    if(op==='neg')return -value(a,env);
    if(functions[op])return functions[op](value(a,env));
    const x=value(a,env),y=value(b,env);
    if(op==='add')return x+y;if(op==='sub')return x-y;if(op==='mul')return x*y;
    if(op==='div')return x/y;if(op==='pow')return x**y;
    throw new Error('Unsupported arithmetic operation: '+op);
  }
  function drift(model,q,params=model.params,removed=[]){
    const off=new Set(removed),n=q.length,out=Array(n).fill(0);
    if(model.kind==='equations'){
      const env={...params,...Object.fromEntries(model.variables.map((v,i)=>[v,q[i]]))};
      for(const term of model.terms)if(!off.has(term.id))out[term.variable]+=value(term.tree,env);
    }else if(model.kind==='gene'){
      out.fill(params.alpha);
      model.graph.edges.forEach(([i,j,s],k)=>{if(off.has('edge-'+k))return;const x=q[i]**params.n;out[j]*=s>0?x/(1+x):1/(1+x);});
      if(!off.has('decay'))for(let i=0;i<n;i++)out[i]-=q[i];
    }else if(model.kind==='rotor'){
      const g=model.graph,m=g.m;
      g.src.forEach((i,k)=>{
        if(off.has('edge-'+k))return;
        const j=g.tgt[k],a=off.has('align')?0:params.J0*params.lam*g.wa[k]*m*Math.sin(m*(q[i]-q[j]));
        const c=off.has('reflect')?0:params.G0*params.lam*g.wc[k]*m*Math.sin(m*(q[i]+q[j]-2*g.phi[k]));
        out[i]-=params.mobility*(a+c);out[j]+=params.mobility*(a-c);
      });
    }
    return out;
  }
  function evaluator(model,params,removed){
    const off=new Set(removed),n=model.variables.length;
    if(model.kind==='rotor'){
      const g=model.graph,m=g.m,edges=g.src.map((i,k)=>[i,g.tgt[k],g.phi[k],off.has('align')?0:params.mobility*params.J0*params.lam*g.wa[k]*m,off.has('reflect')?0:params.mobility*params.G0*params.lam*g.wc[k]*m]).filter((_,k)=>!off.has('edge-'+k));
      return q=>{const out=Array(n).fill(0);for(const [i,j,phi,wa,wc] of edges){const a=wa*Math.sin(m*(q[i]-q[j])),c=wc*Math.sin(m*(q[i]+q[j]-2*phi));out[i]-=a+c;out[j]+=a-c;}return out;};
    }
    if(model.kind==='gene'){
      const edges=model.graph.edges.filter((_,k)=>!off.has('edge-'+k)),decay=!off.has('decay');
      return q=>{const out=Array(n).fill(params.alpha),powers=q.map(x=>x**params.n);for(const [i,j,s] of edges)out[j]*=s>0?powers[i]/(1+powers[i]):1/(1+powers[i]);if(decay)for(let i=0;i<n;i++)out[i]-=q[i];return out;};
    }
    const terms=model.terms.filter(t=>!off.has(t.id)),env={...params};
    return q=>{model.variables.forEach((v,i)=>{env[v]=q[i];});const out=Array(n).fill(0);for(const term of terms)out[term.variable]+=value(term.tree,env);return out;};
  }
  function distance(model,a,b){
    const P=model.carrier.period;
    return Math.hypot(...a.map((v,i)=>model.carrier.kind==='torus'?((v-b[i]+P/2)%P+P)%P-P/2:v-b[i]));
  }
  function integrate(model,options={}){
    const params=options.params||model.params,removed=options.removed||[],T=options.duration||model.duration;
    const pulse=options.pulse||0,start=options.pulseStart??1,end=start+(options.pulseDuration??2),axis=options.axis||0;
    let q=[...(options.initial||model.initial)],t=0,steps=0;
    const native=evaluator(model,params,removed);
    const data=[[0,...q]],interval=T/350;
    let next=interval;
    const project=a=>model.carrier.kind==='orthant'?a.map(x=>Math.max(0,x)):a;
    const f=(a,drive)=>{const out=native(project(a));out[axis]+=drive;return out;};
    while(t<T-1e-10&&steps<60000){
      let dt=Math.min(options.dt||model.dt,T-t);
      if(t<start&&t+dt>start)dt=start-t;
      if(t<end&&t+dt>end)dt=end-t;
      const drive=t+dt/2>=start&&t+dt/2<end?pulse:0;
      const a=f(q,drive),b=f(q.map((x,i)=>x+dt*a[i]/2),drive),c=f(q.map((x,i)=>x+dt*b[i]/2),drive),d=f(q.map((x,i)=>x+dt*c[i]),drive);
      q=project(q.map((x,i)=>x+dt*(a[i]+2*b[i]+2*c[i]+d[i])/6));t+=dt;steps++;
      if(q.some(x=>!Number.isFinite(x)||Math.abs(x)>1e5))return {data,final:q,status:'unbounded',time:t};
      if(t>=next||t>=T-1e-10){data.push([t,...q]);next+=interval;}
    }
    return {data,final:q,status:steps>=60000?'step_limit':'complete',time:t};
  }
  function design(model,params){
    if(model.design==='schlogl'){
      const a=Math.sqrt(3*params.k3);return {...params,a,b:a**3/27};
    }
    if(model.design==='promoters')return {...params,gamma:1};
    if(model.design==='tubes')return {...params,L2:params.L1};
    return {...params};
  }
  const api={value,drift,integrate,distance,design};
  if(typeof module!=='undefined')module.exports=api;else root.FieldBridgeModels=api;
})(typeof window==='undefined'?globalThis:window);
