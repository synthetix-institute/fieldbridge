(function(){
  'use strict';
  const raw=document.getElementById('studio-config').textContent;
  if(raw.trim().startsWith('__'))return;
  const config=JSON.parse(raw),models=config.models;
  if(!models)return;
  const M=window.FieldBridgeModels,$=id=>document.getElementById(id);
  const fmt=(v,n=3)=>Number.isFinite(v)?Number(v.toFixed(n)).toString():'—';
  const esc=s=>String(s).replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
  const repo=p=>config.repo+'/blob/main/'+p;
  const pretty=s=>esc(s).replace(/\*\*\s*([\w.]+|\([^)]*\))/g,'<sup>$1</sup>').replace(/\*\*/g,'^').replace(/\*/g,'·').replace(/\b(alpha|gamma|kappa|eps|mu|phi)\b/g,k=>({alpha:'α',gamma:'γ',kappa:'κ',eps:'ε',mu:'μ',phi:'φ'}[k]));
  const colors=['#287c57','#356eb6','#b84e40','#9160ae','#bd8625','#378d91'];
  let model=models.find(m=>m.id==='toggle')||models[0],state,results,cursor=0,playing=false,frame,last=0;
  let worker=null,revision=0,failedWorker=false;
  function makeWorker(){
    if(typeof Worker==='undefined'||failedWorker)return;
    try{worker=new Worker(config.model_worker||'studio/model-worker.js');}catch(error){failedWorker=true;worker=null;return;}
    worker.onmessage=event=>{if(event.data.id!==revision)return;if(event.data.error){$('model-result').classList.remove('calculating');$('model-status').textContent='Calculation failed: '+event.data.error;return;}results=event.data.results;finish();};
    worker.onerror=()=>{failedWorker=true;worker.terminate();worker=null;calculate();};
  }
  function reset(){
    state={params:{...model.params},removed:[],initial:[...model.initial],alternate:[...model.alternate],duration:model.duration,pulse:0,pulseDuration:2,axis:0,compare:true};cursor=0;playing=false;
    controls();calculate();
  }
  const fields=[...new Set(models.map(m=>m.field))];
  $('model-select').innerHTML=fields.map(field=>`<optgroup label="${esc(field)}">${models.filter(m=>m.field===field).map(m=>`<option value="${m.id}">${esc(m.name)}</option>`).join('')}</optgroup>`).join('');
  function parameter(key,label,value,min,max,step,group){return `<div class="parameter"><div class="parameter-top"><label for="model-${group}-${key}">${esc(label)}</label><output id="model-value-${group}-${key}">${fmt(value)}</output></div><input id="model-${group}-${key}" data-model-param="${key}" data-group="${group}" aria-label="${esc(label)}" type="range" min="${min}" max="${max}" step="${step===1?1:'any'}" value="${value}"></div>`;}
  function controls(){
    $('model-compare').checked=state.compare;
    $('model-select').value=model.id;$('model-name').textContent=model.name;$('model-question').textContent=model.question;
    $('model-terms').innerHTML=model.terms.map(term=>`<label class="model-term"><input type="checkbox" data-model-term="${term.id}" ${state.removed.includes(term.id)?'':'checked'} aria-label="Attach ${esc(term.label)}"><span>${esc(term.label)}<small>${pretty(term.text)}</small></span></label>`).join('');
    $('model-parameters').innerHTML=Object.entries(state.params).map(([k,v])=>{const [lo,hi,step]=model.ranges[k];return parameter(k,k,v,Math.min(lo,v),Math.max(hi,v),step,'coefficient');}).join('');
    $('model-preparation').innerHTML=['initial','alternate'].map((group,j)=>`<h3>Preparation ${j+1}</h3>`+model.variables.map((v,i)=>parameter(i,v+'(0)',state[group][i],model.carrier.kind==='orthant'?0:-model.carrier.scale*2,model.carrier.kind==='torus'?model.carrier.period:model.carrier.scale*2,.02,group)).join('')).join('');
    $('model-pulse-controls').innerHTML=`<label for="model-pulse-axis">Forced coordinate</label><select id="model-pulse-axis">${model.variables.map((v,i)=>`<option value="${i}">${esc(v)}</option>`).join('')}</select>`+parameter('pulse','Force amplitude',state.pulse,-2,2,.05,'protocol')+parameter('pulseDuration','Force duration',state.pulseDuration,.2,6,.2,'protocol')+parameter('duration','Observation time',state.duration,5,60,1,'protocol');
    $('model-pulse-axis').value=state.axis;
    $('model-design').hidden=!model.design;
    $('model-design').innerHTML=model.design?`<button id="model-solve"><i data-lucide="wand-sparkles"></i>${model.design==='schlogl'?'Construct symmetric write':model.design==='promoters'?'Equalize promoters':'Equalize tube lengths'}</button><a href="${repo('docs/tutorial/19_memory_transfer_and_design.md#4-removing-the-obstruction')}" target="_blank" rel="noopener">Derivation ↗</a>`:'';
    if($('model-solve'))$('model-solve').addEventListener('click',()=>{state.params=M.design(model,state.params);controls();calculate();$('model-design-result').textContent=model.design==='schlogl'?`Triple root: x = ${fmt(state.params.a/3)}, a = ${fmt(state.params.a)}, b = ${fmt(state.params.b)}. F(x) = F′(x) = F″(x) = 0.`:'Exchange symmetry restored. The trajectories are recalculated with equal physical parameters.';});
    document.querySelectorAll('[data-model-param]').forEach(el=>el.addEventListener('input',()=>{
      const key=el.dataset.modelParam,value=Number(el.value),group=el.dataset.group;
      if(group==='coefficient')state.params[key]=value;else if(group==='initial'||group==='alternate')state[group][Number(key)]=value;else state[key]=value;
      $('model-value-'+group+'-'+key).textContent=fmt(value);$('model-design-result').textContent='';calculate();
    }));
    document.querySelectorAll('[data-model-term]').forEach(el=>el.addEventListener('change',()=>{remove(el.dataset.modelTerm,!el.checked);calculate();}));
    $('model-pulse-axis').addEventListener('change',e=>{state.axis=Number(e.target.value);calculate();});
    $('model-observable').innerHTML=model.variables.map((v,i)=>`<option value="${i}">${esc(v)}</option>`).join('');
    $('model-assumptions').innerHTML=`<p>${esc(model.closure)}</p><p>${esc(model.assumptions.join('; '))}. Deterministic drift; the source noise strength is D = ${fmt(model.noise)}.</p>`;
    $('model-sources').innerHTML=`<a href="${repo(model.specification)}" target="_blank" rel="noopener">Source equations</a> · ${model.source_url?`<a href="${model.source_url}" target="_blank" rel="noopener">${esc(model.source)}</a>`:esc(model.source)} · <a href="${repo(model.tutorial)}" target="_blank" rel="noopener">Tutorial &amp; Python calculation ↗</a>`;
    $('model-saved').hidden=!model.saved;
    $('model-saved').parentElement.hidden=!model.saved;
    $('model-saved').innerHTML=model.saved?`<img src="${model.saved.image}" alt="Full Python calculation for the unmodified ${esc(model.name)}"><p><a href="${model.saved.record}">Full calculation record</a> · Original parameter values; stable-state search and retention calculation.</p>`:'';
    window.lucide.createIcons();
  }
  function remove(id,off){state.removed=state.removed.filter(t=>t!==id);if(off)state.removed.push(id);}
  function equations(){
    if(model.kind==='equations')return model.variables.map((v,i)=>`<div><span>d${esc(v)}/dt = </span>${model.terms.filter(t=>t.variable===i).map((t,j)=>{const negative=t.tree[0]==='neg',text=negative?(t.text.startsWith('-(')?t.text.slice(2,-1):t.text.slice(1)):t.text;return `<button class="equation-term ${state.removed.includes(t.id)?'removed':''}" data-equation-term="${t.id}" aria-pressed="${!state.removed.includes(t.id)}" title="${state.removed.includes(t.id)?'Attach':'Detach'} ${esc(t.label)}">${negative?'− ':j?'+ ':''}${pretty(text)}</button>`;}).join(' ')}${state.pulse&&state.axis===i?' + f(t)':''}</div>`).join('');
    if(model.kind==='gene')return model.variables.map((v,i)=>'<div>d'+esc(v)+'/dt = α'+model.graph.edges.map(([j,target,sign],k)=>target===i&&!state.removed.includes('edge-'+k)?' H'+(sign>0?'₊':'₋')+'('+esc(model.variables[j])+')':'').join('')+(state.removed.includes('decay')?'':' − '+esc(v))+(state.pulse&&state.axis===i?' + f(t)':'')+'</div>').join('')+'<small>H₋(q) = 1/(1 + qⁿ) · H₊(q) = qⁿ/(1 + qⁿ)</small>';
    const active=(!state.removed.includes('align')||!state.removed.includes('reflect'))&&model.links.some(link=>!state.removed.includes(link.term));
    return '<div>E = '+(active?'−Σ<sub>attached bonds</sub> ['+(state.removed.includes('align')?'': 'Jᵢⱼ cos m(θᵢ − θⱼ)')+(!state.removed.includes('align')&&!state.removed.includes('reflect')?' + ':'')+(state.removed.includes('reflect')?'':'Gᵢⱼ cos m(θᵢ + θⱼ − 2φᵢⱼ)')+']':'0')+'</div><small>dθᵢ/dt = −μ ∂E/∂θᵢ'+(state.pulse?' + fᵢ(t)':'')+' · m = '+model.graph.m+'</small>';
  }
  function calculate(){
    playing=false;cancelAnimationFrame(frame);cursor=state.duration;
    revision++;if(worker){worker.terminate();makeWorker();}
    $('model-equations').innerHTML=equations();
    document.querySelectorAll('[data-equation-term]').forEach(el=>el.addEventListener('click',()=>{remove(el.dataset.equationTerm,!state.removed.includes(el.dataset.equationTerm));controls();calculate();}));
    $('model-drive').textContent=state.pulse?`f(t) = ${fmt(state.pulse)} on ${model.variables[state.axis]} for 1 ≤ t < ${fmt(1+state.pulseDuration)}; zero otherwise.`:'No external force.';
    $('model-time').max=state.duration;$('model-time').value=cursor;
    if(worker){results=null;$('model-play').disabled=true;$('model-time').disabled=true;$('model-export').disabled=true;$('model-final').textContent='';$('model-metrics').innerHTML='';$('model-scene').innerHTML='';$('model-plot').width=$('model-plot').width;$('model-status').textContent='Integrating the current equations…';$('model-result').classList.add('calculating');worker.postMessage({id:revision,model,state});return;}
    results=[M.integrate(model,state),M.integrate(model,{...state,initial:state.alternate})];
    if(state.compare)results.push(M.integrate(model,{...state,params:model.params,removed:[]}));
    finish();
  }
  function finish(){
    $('model-result').classList.remove('calculating');$('model-play').disabled=false;$('model-time').disabled=false;$('model-export').disabled=false;
    const first=results[0],second=results[1],problem=results.find(r=>r.status!=='complete'),complete=!problem,velocity=M.drift(model,first.final,state.params,state.removed);
    if(state.duration>=1&&state.duration<1+state.pulseDuration)velocity[state.axis]+=state.pulse;
    const speed=Math.hypot(...velocity);
    $('model-metrics').innerHTML=[['Preparation separation at T',complete?fmt(M.distance(model,first.final,second.final)):'—'],['Speed at T',first.status==='complete'?fmt(speed):'—'],['Change from original model',complete&&state.compare?fmt(M.distance(model,first.final,results[2].final)):'—'],['Detached contributions',String(state.removed.length)]].map(([k,v])=>`<dt>${k}</dt><dd>${v}</dd>`).join('');
    $('model-status').textContent=complete?`t = ${fmt(state.duration)} · drift integration`:`A trajectory stopped at t = ${fmt(problem.time)}: ${problem.status==='unbounded'?'trajectory grows beyond the plotted range':'integration limit reached'}`;
    $('model-final').textContent=model.variables.map((v,i)=>`${v}(${first.status==='complete'?'T':fmt(first.time)}) = ${fmt(first.final[i])}`).join(' · ');
    paint();
  }
  function points(){
    if(model.positions){const xs=model.positions.map(p=>p[0]),ys=model.positions.map(p=>p[1]),xmin=Math.min(...xs),ymin=Math.min(...ys);return model.positions.map(p=>[.12+.76*(p[0]-xmin)/(Math.max(...xs)-xmin||1),.16+.68*(p[1]-ymin)/(Math.max(...ys)-ymin||1)]);}
    const n=model.variables.length;return model.variables.map((_,i)=>n===1?[.5,.5]:[.5+.31*Math.cos(-Math.PI/2+2*Math.PI*i/n),.5+.31*Math.sin(-Math.PI/2+2*Math.PI*i/n)]);
  }
  function bond(link,positions){
    const a=positions[link.source].map((x,i)=>x*(i?270:500)),b=positions[link.target].map((x,i)=>x*(i?270:500));
    const off=state.removed.includes(link.term),term=model.terms.find(t=>t.id===link.term),label=term?term.label:'bond '+model.variables[link.source]+' to '+model.variables[link.target];
    if(link.source===link.target){
      const siblings=model.links.filter(item=>item.source===link.source&&item.target===link.target),index=siblings.indexOf(link),angle=model.variables.length===1?index*120:link.source%2?180:0,radius=55+index*12;
      return `<g class="model-bond ${off?'removed':''}" data-bond="${link.term}" tabindex="0" role="button" aria-label="${off?'Attach':'Detach'} ${esc(label)}"><path transform="rotate(${angle} ${a[0]} ${a[1]})" d="M${a[0]-20},${a[1]-12} C${a[0]-radius*2},${a[1]-radius} ${a[0]-radius*2},${a[1]+radius} ${a[0]-20},${a[1]+12}" marker-end="url(#model-arrow)"/><title>${off?'Attach':'Detach'} ${esc(label)}</title></g>`;
    }
    const directed=model.kind!=='rotor',dx=b[0]-a[0],dy=b[1]-a[1],length=Math.hypot(dx,dy),radius=directed?25:0;
    const x1=a[0]+radius*dx/length,y1=a[1]+radius*dy/length,x2=b[0]-radius*dx/length,y2=b[1]-radius*dy/length;
    const bend=directed?22:0,cx=(x1+x2)/2-bend*dy/length,cy=(y1+y2)/2+bend*dx/length;
    return `<g class="model-bond ${off?'removed':''}" data-bond="${link.term}" tabindex="0" role="button" aria-label="${off?'Attach':'Detach'} ${esc(label)}"><path d="M${x1},${y1} Q${cx},${cy} ${x2},${y2}" ${directed?`marker-end="url(#${link.sign<0?'model-repression':'model-arrow'})"`:''}/><title>${off?'Attach':'Detach'} ${esc(label)}</title></g>`;
  }
  function scene(){
    const data=results[0].data,row=data.find(r=>r[0]>=cursor)||data.at(-1),q=row.slice(1),positions=points();
    $('model-scene').innerHTML=`<svg viewBox="0 0 500 270" role="img" aria-label="${esc(model.carrier.name)} at time ${fmt(cursor)}"><defs><marker id="model-arrow" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="5" markerHeight="5" orient="auto"><path d="M0,0 L10,5 L0,10" fill="#598771"/></marker><marker id="model-repression" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="5" markerHeight="5" orient="auto"><path d="M9,0 L9,10" stroke="#598771" stroke-width="2"/></marker></defs>${model.links.map(link=>bond(link,positions)).join('')}${positions.map(([x,y],i)=>{const cx=x*500,cy=y*270,angle=q[i],rotor=model.kind==='rotor';return `<g><circle cx="${cx}" cy="${cy}" r="${rotor?7:22}" fill="${colors[i%colors.length]}" opacity=".18"/>${rotor?`<line x1="${cx-15*Math.cos(angle)}" y1="${cy-15*Math.sin(angle)}" x2="${cx+15*Math.cos(angle)}" y2="${cy+15*Math.sin(angle)}" stroke="${colors[i%colors.length]}" stroke-width="6" stroke-linecap="round"/>${model.graph.m===1?`<circle cx="${cx+15*Math.cos(angle)}" cy="${cy+15*Math.sin(angle)}" r="4" fill="white" stroke="${colors[i%colors.length]}" stroke-width="2"/>`:''}`:`<text x="${cx}" y="${cy+5}" text-anchor="middle" fill="${colors[i%colors.length]}" font-size="14">${esc(model.variables[i])}</text>`}<text x="${cx}" y="${cy+37}" text-anchor="middle" font-size="12" fill="#35453c">${fmt(q[i],2)}</text></g>`;}).join('')}</svg>`;
    document.querySelectorAll('[data-bond]').forEach(el=>{const action=()=>{remove(el.dataset.bond,!state.removed.includes(el.dataset.bond));controls();calculate();};el.addEventListener('click',action);el.addEventListener('keydown',e=>{if(['Enter',' '].includes(e.key)){e.preventDefault();action();}});});
    $('model-clock').textContent='t = '+fmt(cursor,2);
  }
  function plot(){
    if(!results)return;
    const canvas=$('model-plot'),box=canvas.getBoundingClientRect();if(!box.width)return;
    const dpr=window.devicePixelRatio||1;canvas.width=box.width*dpr;canvas.height=box.height*dpr;
    const ctx=canvas.getContext('2d');ctx.scale(dpr,dpr);const W=box.width,H=box.height,l=42,t=20,w=W-56,h=H-53,axis=Number($('model-observable').value||0);
    const values=results.flatMap(r=>r.data.map(row=>row[axis+1])).filter(Number.isFinite),min=Math.min(...values),max=Math.max(...values),span=Math.max(max-min,.2),lo=min-.1*span,hi=max+.1*span;
    const X=x=>l+w*x/state.duration,Y=y=>t+h*(hi-y)/(hi-lo);ctx.clearRect(0,0,W,H);ctx.font='11px Arial';
    for(let i=0;i<=4;i++){const v=lo+(hi-lo)*i/4;ctx.strokeStyle='#e5ebe7';ctx.beginPath();ctx.moveTo(l,Y(v));ctx.lineTo(l+w,Y(v));ctx.stroke();ctx.fillStyle='#59675f';ctx.fillText(fmt(v,1),2,Y(v)+4);ctx.fillText(fmt(state.duration*i/4,1),X(state.duration*i/4)-7,H-12);}
    if(state.pulse){ctx.fillStyle='#287c5715';ctx.fillRect(X(1),t,X(Math.min(state.duration,1+state.pulseDuration))-X(1),h);}
    results.forEach((r,i)=>{ctx.strokeStyle=['#356eb6','#b84e40','#78847d'][i];ctx.lineWidth=i===2?1.5:2;ctx.setLineDash(i===2?[4,4]:[]);ctx.beginPath();r.data.forEach((row,k)=>k?ctx.lineTo(X(row[0]),Y(row[axis+1])):ctx.moveTo(X(row[0]),Y(row[axis+1])));ctx.stroke();});
    ctx.setLineDash([]);ctx.strokeStyle='#303d35';ctx.beginPath();ctx.moveTo(X(cursor),t);ctx.lineTo(X(cursor),t+h);ctx.stroke();ctx.fillStyle='#35453c';ctx.fillText(model.variables[axis],l,13);
  }
  function paint(){if(!results)return;scene();plot();$('model-play').setAttribute('aria-label',playing?'Pause evolution':'Play evolution');$('model-play').innerHTML=`<i data-lucide="${playing?'pause':'play'}"></i>`;window.lucide.createIcons();}
  function animate(time){if(!playing)return;const delta=last?(time-last)/1000:0;last=time;cursor=Math.min(state.duration,cursor+delta*state.duration/8);$('model-time').value=cursor;paint();if(cursor<state.duration)frame=requestAnimationFrame(animate);else{playing=false;paint();}}
  $('model-select').addEventListener('change',e=>{model=models.find(m=>m.id===e.target.value);$('model-design-result').textContent='';reset();});
  $('model-reset').addEventListener('click',reset);
  $('model-observable').addEventListener('change',plot);
  $('model-compare').addEventListener('change',e=>{state.compare=e.target.checked;calculate();});
  $('model-time').addEventListener('input',e=>{playing=false;cancelAnimationFrame(frame);cursor=Number(e.target.value);paint();});
  $('model-play').addEventListener('click',()=>{playing=!playing;if(cursor>=state.duration)cursor=0;last=0;paint();if(playing)frame=requestAnimationFrame(animate);else cancelAnimationFrame(frame);});
  $('model-export').addEventListener('click',()=>{const record={model:model.specification,params:state.params,removed:state.removed,initial:state.initial,alternate:state.alternate,protocol:{amplitude:state.pulse,start:1,duration:state.pulseDuration,coordinate:model.variables[state.axis]},duration:state.duration,results,approximation:'deterministic drift; original stochastic noise omitted',tutorial:repo(model.tutorial)};const url=URL.createObjectURL(new Blob([JSON.stringify(record,null,2)],{type:'application/json'}));const a=document.createElement('a');a.href=url;a.download='fieldbridge-'+model.id+'-construction.json';a.click();setTimeout(()=>URL.revokeObjectURL(url),1000);});
  window.addEventListener('fieldbridge-module',e=>{if(e.detail==='gallery')paint();else{playing=false;cancelAnimationFrame(frame);}});
  new ResizeObserver(()=>{if(results)plot();}).observe($('model-plot'));makeWorker();reset();
})();
