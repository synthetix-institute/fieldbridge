(function () {
  'use strict';
  if(document.getElementById('studio-config').textContent.trim().startsWith('__'))return;
  const P=window.FieldBridgePhysics, config=JSON.parse(document.getElementById('studio-config').textContent);
  const $=id=>document.getElementById(id), icon=name=>`<i data-lucide="${name}" aria-hidden="true"></i>`;
  const defaults={memory:{feedback:true,saturation:true,field:true,thermal:false,eps:1,gamma:1,h:0,D:.04,pulse:.8,pulseDuration:2,realization:'normalized',scale:2},
    stochastic:{noise:true,correction:true,map:'square',convention:'ito',theta:1,alpha:1,sigma:1,mu:.2}};
  const states={memory:structuredClone(defaults.memory),stochastic:structuredClone(defaults.stochastic)};
  let lab='memory', state=states.memory, result=null, plotMode='trajectory';
  const refreshIcons=()=>window.lucide.createIcons();
  const fmt=(v,n=3)=>v===null?'—':Number.isFinite(v)?Number(v.toFixed(n)).toString():'∞';
  function sourceLink(key){const s=config.sources[key];return `<a href="${s.url}" target="_blank" rel="noopener">${s.authors.split(',')[0]} (${s.year})</a>`;}
  function toggle(key,label,math){return `<div class="part-row"><label for="part-${key}">${label}<small>${math}</small></label><label class="switch"><input id="part-${key}" type="checkbox" data-part="${key}" ${state[key]?'checked':''} aria-label="Attach ${label.toLowerCase()}"><span></span></label></div>`;}
  function slider(key,label,min,max,step,enabled=true){return `<div class="parameter ${enabled?'':'inactive'}"><div class="parameter-top"><label for="param-${key}">${label}</label><output id="value-${key}" for="param-${key}">${fmt(state[key],2)}</output></div><input id="param-${key}" data-parameter="${key}" type="range" min="${min}" max="${max}" step="${step}" value="${state[key]}" ${enabled?'':'disabled'}></div>`;}
  function setPlotMode(mode){plotMode=mode;document.querySelectorAll('[data-plot]').forEach(e=>{e.classList.toggle('selected',e.dataset.plot===mode);e.setAttribute('aria-pressed',String(e.dataset.plot===mode));});}
  function renderControls(){
    if(lab==='memory'){
      $('parts').innerHTML=toggle('feedback','Linear feedback','εx')+toggle('saturation','Nonlinear saturation','−γx³')+toggle('field','Writing field','h')+toggle('thermal','Thermal fluctuations','√(2D) dW');
      $('parameters').innerHTML=slider('eps','Feedback ε',-.8,1.6,.05,state.feedback)+slider('gamma','Saturation γ',.4,1.8,.05,state.saturation)+slider('pulse','Pulse amplitude',-1.5,1.5,.05,state.field)+slider('pulseDuration','Pulse duration',.2,5,.2,state.field)+'<details class="extra-parameters" '+(state.thermal?'open':'')+'><summary>Bias and fluctuations</summary>'+slider('h','Constant field h',-.8,.8,.02,state.field)+(state.thermal?slider('D','Noise strength D',.01,.2,.01):'')+'</details>';
      $('realization-control').innerHTML='<label for="realization">Physical coordinate</label><select id="realization"><option value="normalized">Order parameter x</option><option value="displacement">Displacement q = q₀x</option></select>'+(state.realization==='displacement'?slider('scale','Length scale q₀',.5,4,.1):'');
      $('realization').value=state.realization;
      $('realization').addEventListener('change',e=>{state.realization=e.target.value;renderControls();calculate();});
    }else{
      $('parts').innerHTML=toggle('noise','Brownian fluctuations','dW² = dt')+toggle('correction','Derived drift term',state.map==='square'?'+σ²':'−σ²/2 (Itô)');
      $('parameters').innerHTML=state.map==='square'?slider('theta','Radial parameter θ',.2,2,.1)+slider('alpha','Restoring rate α',.2,2,.1)+slider('sigma','Noise amplitude σ',.2,1,.1):slider('mu','Growth rate μ',-.4,.6,.05)+slider('sigma','Noise amplitude σ',.2,1.5,.1);
      $('realization-control').innerHTML='<label for="state-map">State transformation α</label><select id="state-map"><option value="square">Y = X² · radial → intensity</option><option value="log">Y = log X · growth → log signal</option></select>'+(state.map==='log'?'<label for="convention" style="margin-top:12px">Source convention</label><select id="convention"><option value="ito">Itô</option><option value="stratonovich">Stratonovich</option></select>':'');
      $('state-map').value=state.map;
      $('state-map').addEventListener('change',e=>{state.map=e.target.value;if(state.map==='square')state.sigma=Math.min(1,state.sigma);renderControls();calculate();});
      if($('convention')){$('convention').value=state.convention;$('convention').addEventListener('change',e=>{state.convention=e.target.value;calculate();});}
    }
    document.querySelectorAll('[data-part]').forEach(el=>el.addEventListener('change',()=>{const id=el.id;state[el.dataset.part]=el.checked;renderControls();calculate();$(id).focus();}));
    document.querySelectorAll('[data-parameter]').forEach(el=>el.addEventListener('input',()=>{state[el.dataset.parameter]=Number(el.value);$('value-'+el.dataset.parameter).textContent=fmt(Number(el.value),2);calculate();}));
    refreshIcons();
  }
  function token(text,on,key){return `<button class="term-token ${on?'':'off'}" data-graph-part="${key}" aria-pressed="${on}" title="${on?'Detach':'Attach'} ${key}">${text}${icon(on?'minus':'plus')}</button>`;}
  function graph(){
    const mem=lab==='memory';
    const tokens=mem?token(fmt(state.eps)+'x',state.feedback,'feedback')+token('−'+fmt(state.gamma)+'x³',state.saturation,'saturation')+token('h(t)',state.field,'field')+token('√(2D)dW',state.thermal,'thermal'):token('½σ²∂²',state.noise,'noise')+token(state.map==='square'?'+σ²∂':'δb∂',state.correction,'correction');
    const attached=mem?[state.feedback,state.saturation,state.field,state.thermal].filter(Boolean).length:[state.noise,state.correction].filter(Boolean).length;
    $('part-count').textContent=`${attached} / ${mem?4:2}`;
    $('graph').innerHTML=`<div class="core-group"><div class="group-title">Mechanism core M = (Ω, Ξ)</div><div class="graph-core"><div class="graph-node"><span class="node-symbol">Ω</span><strong>${mem?'Overdamped evolution':'Stochastic generator'}</strong><div class="term-tokens">${tokens}</div></div><div class="graph-connection">${icon('arrow-right')}</div><div class="graph-node"><span class="node-symbol">Ξ</span><strong>${mem?'x ∈ ℝ':state.map==='square'?'X > 0 → Y > 0':'X > 0 → Y ∈ ℝ'}</strong><small>${mem?'One continuous coordinate':'State map: '+(state.map==='square'?'Y = X²':'Y = log X')}</small></div></div><div class="clauses"><span>${mem?'C: V(x), mobility = 1':'C: '+(state.map==='square'?'Itô':state.convention)+', interior domain'}</span><span>${mem?'R: sign(x)':'R: E[Y(t)]'}</span></div></div><div class="attachment-path"></div><div class="clauses"><span>${mem?(state.field?'P: h = '+fmt(state.pulse)+' for '+fmt(state.pulseDuration)+' time units; then '+fmt(state.h):'P: unforced evolution'):'P: X(0) = 1'}</span></div><div class="attachment-path"></div><div class="realization-node">A · ${mem?(state.realization==='displacement'?'Displacement q = '+fmt(state.scale)+'x':'Order parameter x'):(state.map==='square'?'Intensity coordinate':'Logarithmic signal')}</div>`;
    document.querySelectorAll('[data-graph-part]').forEach(b=>b.addEventListener('click',()=>{state[b.dataset.graphPart]=!state[b.dataset.graphPart];renderControls();calculate();}));
    refreshIcons();
  }
  function metrics(rows){$('metrics').innerHTML=rows.map(([k,v])=>`<dt>${k}</dt><dd>${v}</dd>`).join('');}
  function memoryOutcome(){
    const m=result, scale=state.realization==='displacement'?state.scale:1, symbol=state.realization==='displacement'?'q':'x';
    let terms=[];if(state.feedback)terms.push(`${fmt(state.eps)}${symbol}`);if(state.saturation)terms.push(`− ${fmt(state.gamma/(scale*scale))}${symbol}<sup>3</sup>`);if(state.field)terms.push('+ h(t)');
    $('equation').innerHTML=`d${symbol}/dt = ${terms.join(' ')||'0'}${state.thermal?` + ${fmt(Math.sqrt(2*state.D)*scale)} ξ(t)`:''}`;
    $('drive-equation').textContent=state.field?`h(t) = ${fmt(state.pulse*scale)} for 1 ≤ t < ${fmt(1+state.pulseDuration)}; ${fmt(state.h*scale)} otherwise`:'h(t) = 0';
    $('equation-note').textContent=state.realization==='displacement'?'q = q₀x; q₀ = '+fmt(state.scale)+'. Mobility in q is q₀² times mobility in x; the coordinate change preserves the dynamics.':'V(x) = −('+fmt(m.a)+')x²/2 + ('+fmt(m.b)+')x⁴/4 − ('+fmt(m.h)+')x. Noise strength D = '+fmt(state.thermal?state.D:0)+'. All quantities are dimensionless.';
    const text=m.bistable?(state.thermal?'Two attracting states encode opposite signs; thermal escape makes retention finite.':'Two attracting states retain opposite signs after the writing field is removed.'):m.neutral?'Every coordinate is stationary. There is no restoring force that corrects a perturbation.':!m.bounded?'Removing saturation leaves outward growth without finite attracting states.':m.stable.length===1?'One attracting state erases the distinction between the two initial preparations.':'At this parameter value the restoring slope vanishes; the linear relaxation test is marginal.';
    $('consequence').className='consequence'+(m.bistable?'':' warning');$('consequence').querySelector('p').textContent=text;
    const pulse=state.field?state.pulse:0;
    const trajectories=[-1,1].map(x=>P.trajectory(m,x,pulse,8,800,{start:1,end:1+state.pulseDuration}));
    const ends=trajectories.map(t=>t.diverged?'diverges':fmt(t.data[t.data.length-1][1]*scale));
    const fold=m.threshold===null?null:m.threshold*scale;
    metrics(plotMode==='trajectory'?[['Final state from −'+fmt(scale),ends[0]],['Final state from +'+fmt(scale),ends[1]],state.thermal?['Escape time τ',m.tau===null?'weak-noise limit unmet':fmt(m.tau,1)]:['Writing threshold |h꜀|',fmt(fold)],['Stable states',m.stable.length?m.stable.map(e=>fmt(e.x*scale)).join(', '):'none']]:[['Stable states',m.stable.length?m.stable.map(e=>fmt(e.x*scale)).join(', '):'none'],['Barrier at h = 0',m.barrier===null?'—':fmt(m.barrier)],['Fold field |h꜀|',fmt(fold)],['Escape time τ',!state.thermal?'deterministic':m.tau===null?'weak-noise limit unmet':fmt(m.tau,1)]]);
    if(plotMode==='trajectory'&&m.bistable&&m.h===0&&!trajectories.some(t=>t.diverged)){
      const final=trajectories.map(t=>t.data[t.data.length-1][1]);
      $('consequence').querySelector('p').textContent=final.every(x=>x>0)?'The pulse switches the negative preparation into the positive basin. The written state persists after the pulse.':final.every(x=>x<0)?'The pulse switches the positive preparation into the negative basin. The written state persists after the pulse.':'The two preparations remain in opposite basins. This pulse has not switched the stored state.';
      if(state.thermal)$('consequence').querySelector('p').textContent+=' Retention is limited by thermal escape.';
    }
    $('detail-title').textContent=plotMode==='potential'?'Calculated equilibrium structure':'Predicted writing response';
    $('detail-text').innerHTML=(state.field?`h(t) = ${fmt(state.pulse*scale)} from t = 1 to ${fmt(1+state.pulseDuration)}, then ${fmt(state.h*scale)}. Blue and red curves start at ${symbol} = −${fmt(scale)} and +${fmt(scale)}. `:'h(t) = 0. Blue and red curves start in opposite preparations. ')+(state.thermal?'Curves show deterministic drift; τ is a weak-noise escape estimate.':'')+` <a href="${config.repo}/blob/main/tests/test_web_demo.py" target="_blank" rel="noopener">Calculation tests ↗</a>`;
    $('source-links').innerHTML=sourceLink('kramers')+' · Model: scalar normal form';
    $('lesson-link').href=config.repo+'/blob/main/docs/tutorial/18_memory_writing_and_retention.md#scalar-constructor-demonstration';
  }
  function stochasticOutcome(){
    const r=result,square=state.map==='square', sigma=state.noise?state.sigma:0;
    $('equation').innerHTML=square?`dY = (${fmt(r.candidate)} − ${fmt(2*state.alpha)}Y)dt<br>+ ${fmt(2*sigma)}√Y dW`:`dY = (${fmt(r.candidate)})dt + ${fmt(sigma)} dW`;
    $('drive-equation').textContent='';
    $('equation-note').textContent=square?'Source: dX = [(θ + ½)/X − αX]dt + σdW; Y = X². Positive interior domain; finite first moment and no boundary contribution.':'Source: dX = μXdt + σXdW ('+(state.convention==='ito'?'Itô':'Stratonovich')+'); Y = log X. X(0) = 1; the target is written in Itô form.';
    const closed=Math.abs(r.residual)<1e-10;
    $('consequence').className='consequence'+(closed?'':' warning');$('consequence').querySelector('p').textContent=closed?(square?'Quadratic variation supplies the extra drift. The target generator preserves the transformed source dynamics.':'The convention fixes the mean logarithmic growth; both conventions give the same variance rate.'):'Removing the derived drift leaves a nonzero generator remainder and shifts the predicted mean. The source process has not changed.';
    metrics([['Required drift',fmt(r.exact)],['Generator remainder Δ',fmt(r.residual)],['Target mean at t = 2',fmt(r.mean(2,false))],['Source − target mean',fmt(r.mean(2,true)-r.mean(2,false))]]);
    $('detail-title').textContent=closed?'Verified construction':'Omission control detected';
    $('detail-text').innerHTML=`${square?'The complete target has Δ = 0. Removing quadratic variation gives Δ = σ², changing the mean without changing the supplied diffusion coefficient.':'Independent source integrators check the derived log mean and variance; swapping the convention changes the mean by σ²/2.'} <a class="test-link" href="${config.repo}/blob/main/tests/${square?'test_verification.py':'test_stochastic_conventions.py'}" target="_blank" rel="noopener">Python tests ↗</a>`;
    $('source-links').innerHTML=sourceLink('ito')+' · Model: worked stochastic transformation';
    $('lesson-link').href=config.repo+'/blob/main/docs/tutorial/10_stochastic_construction.md';
  }
  function draw(){
    if(!['memory','stochastic'].includes(lab)||!result)return;
    const canvas=$('plot'), box=canvas.getBoundingClientRect(), dpr=window.devicePixelRatio||1;
    if(!box.width)return;
    canvas.width=Math.round(box.width*dpr);canvas.height=Math.round(box.height*dpr);
    const ctx=canvas.getContext('2d');ctx.scale(dpr,dpr);
    const W=box.width,H=box.height,pad={l:43,r:14,t:27,b:37},w=W-pad.l-pad.r,h=H-pad.t-pad.b;
    let lines=[],points=[],xlim,ylim,xlabel,ylabel;
    if(lab==='memory'&&plotMode==='potential'){
      const scale=state.realization==='displacement'?state.scale:1;
      xlim=[-2*scale,2*scale];const m=result;lines=[{color:'#287652',data:Array.from({length:201},(_,i)=>{const x=-2+i*.02;return [x*scale,-m.a*x*x/2+m.b*x**4/4-m.h*x];})}];
      points=m.equilibria.map(e=>({x:e.x*scale,y:-m.a*e.x*e.x/2+m.b*e.x**4/4-m.h*e.x,stable:e.stable}));
      const ys=lines[0].data.map(e=>e[1]);ylim=[Math.min(...ys)-.3,Math.max(.5,Math.max(...ys))+.25];xlabel=scale===1?'x':'q';ylabel='V';
      $('plot-heading').textContent='Energy landscape';$('plot-legend').textContent='● attracting   ○ unstable';canvas.setAttribute('aria-label','Calculated potential and equilibrium states');
    }else if(lab==='memory'){
      const scale=state.realization==='displacement'?state.scale:1;
      xlim=[0,8];const pulse=state.field?state.pulse:0;lines=[{color:'#336eb5',data:P.trajectory(result,-1,pulse,8,800,{start:1,end:1+state.pulseDuration}).data.map(([t,x])=>[t,x*scale])},{color:'#b5483b',data:P.trajectory(result,1,pulse,8,800,{start:1,end:1+state.pulseDuration}).data.map(([t,x])=>[t,x*scale]),dash:true}];
      const ys=lines.flatMap(l=>l.data.map(e=>e[1]));ylim=[Math.max(-3*scale,Math.min(-1.2*scale,...ys)-.1),Math.min(3*scale,Math.max(1.2*scale,...ys)+.1)];xlabel='t';ylabel=state.realization==='displacement'?'q(t)':'x(t)';
      $('plot-heading').textContent=state.thermal?'Drift response':'Writing response';$('plot-legend').textContent=`blue −${fmt(scale)} · red +${fmt(scale)} preparation${state.thermal?' · drift only':''}`;canvas.setAttribute('aria-label','Calculated drift trajectories under a writing pulse');
    }else{
      xlim=[0,3];lines=[{color:'#336eb5',data:Array.from({length:101},(_,i)=>[i*.03,result.mean(i*.03,true)])},{color:'#b5483b',dash:true,data:Array.from({length:101},(_,i)=>[i*.03,result.mean(i*.03,false)])}];
      const ys=lines.flatMap(l=>l.data.map(e=>e[1]));ylim=[Math.min(0,...ys)-.1,Math.max(.3,...ys)+.2];xlabel='t';ylabel='E[Y(t)]';
      $('plot-heading').textContent='Predicted mean';$('plot-legend').textContent='blue source   red constructed';canvas.setAttribute('aria-label','Source mean compared with constructed target mean');
    }
    const X=x=>pad.l+(x-xlim[0])/(xlim[1]-xlim[0])*w,Y=y=>pad.t+h-(y-ylim[0])/(ylim[1]-ylim[0])*h;
    ctx.clearRect(0,0,W,H);ctx.font='10px Arial';ctx.textBaseline='middle';
    for(let i=0;i<=4;i++){
      const y=ylim[0]+i*(ylim[1]-ylim[0])/4;ctx.strokeStyle='#e7ede9';ctx.beginPath();ctx.moveTo(pad.l,Y(y));ctx.lineTo(W-pad.r,Y(y));ctx.stroke();ctx.fillStyle='#78827c';ctx.textAlign='right';ctx.fillText(fmt(y,1),pad.l-8,Y(y));
      const x=xlim[0]+i*(xlim[1]-xlim[0])/4;ctx.textAlign='center';ctx.fillText(fmt(x,1),X(x),H-pad.b+14);
    }
    if(lab==='memory'&&plotMode==='trajectory'&&state.field){ctx.fillStyle='#28765212';ctx.fillRect(X(1),pad.t,X(1+state.pulseDuration)-X(1),h);ctx.fillStyle='#287652';ctx.textAlign='center';ctx.fillText('pulse',X(1+state.pulseDuration/2),pad.t+10);}
    ctx.save();ctx.beginPath();ctx.rect(pad.l,pad.t,w,h);ctx.clip();
    lines.forEach(l=>{ctx.strokeStyle=l.color;ctx.lineWidth=2;ctx.setLineDash(l.dash?[5,4]:[]);ctx.beginPath();l.data.forEach(([x,y],i)=>{if(i===0)ctx.moveTo(X(x),Y(y));else ctx.lineTo(X(x),Y(y));});ctx.stroke();});ctx.setLineDash([]);
    points.forEach(p=>{ctx.beginPath();ctx.arc(X(p.x),Y(p.y),4,0,2*Math.PI);ctx.fillStyle=p.stable?'#287652':'white';ctx.fill();ctx.strokeStyle='#287652';ctx.lineWidth=1.5;ctx.stroke();});ctx.restore();
    ctx.strokeStyle='#b6c3bb';ctx.lineWidth=1;ctx.beginPath();ctx.moveTo(pad.l,pad.t);ctx.lineTo(pad.l,pad.t+h);ctx.lineTo(pad.l+w,pad.t+h);ctx.stroke();ctx.fillStyle='#4a5b50';ctx.font='italic 13px Georgia';ctx.textAlign='center';ctx.fillText(xlabel,pad.l+w/2,H-9);ctx.fillText(ylabel,29,12);
  }
  function calculate(){result=lab==='memory'?P.memory(state):P.stochastic(state);$('build-state').textContent='Calculated with the current parts';$('build-state').className='build-state';graph();if(lab==='memory')memoryOutcome();else stochasticOutcome();draw();}
  function switchLab(next){
    if(next===lab)return;
    lab=next;document.querySelectorAll('[data-lab]').forEach(b=>{b.setAttribute('aria-selected',String(b.dataset.lab===lab));b.tabIndex=b.dataset.lab===lab?0:-1;});
    $('workbench').hidden=!['memory','stochastic'].includes(lab);$('applications').hidden=lab!=='applications';
    $('gallery').hidden=lab!=='gallery';$('spins').hidden=lab!=='spins';
    history.replaceState(null,'','#'+lab);
    if(lab==='applications'){renderApplication();return;}
    if(lab==='gallery'||lab==='spins'){window.dispatchEvent(new CustomEvent('fieldbridge-module',{detail:lab}));return;}
    $('workbench').setAttribute('aria-labelledby','tab-'+lab);state=states[lab];plotMode=lab==='memory'?'trajectory':'potential';
    document.querySelectorAll('[data-plot]').forEach(e=>{e.classList.toggle('selected',e.dataset.plot===plotMode);e.setAttribute('aria-pressed',String(e.dataset.plot===plotMode));});
    $('lab-category').textContent=lab==='memory'?'Scalar normal form':'Generator correspondence';$('lab-title').textContent=lab==='memory'?'Memory constructor':'Stochastic transformation';$('plot-modes').hidden=lab!=='memory';renderControls();calculate();
  }
  const targets={
    mechanical:{title:'Mechanical event latch',proposal:'A transient force writes a persistent displacement. Two stable configurations encode whether an event occurred.',requirements:['Two attracting configurations','Write with a finite pulse','Retain after the drive is removed'],source:'snap',known:'Snap-through already provides a strong mechanical response to a small fluidic input.',open:'A useful event latch additionally needs calibrated switching and reset thresholds, a retention interval and robustness to vibration. Those properties are not supplied by the cited actuator study.',tests:[['Calibrate the force curve','Measure F(q), locate stable equilibria and fit the local coefficients ε and γ.'],['Predict a withheld pulse','Use the fitted fold field to select a pulse above threshold, then test whether its state persists after removal.'],['Break the mechanism','Remove the bistable element. The same pulse should produce relaxation rather than retained displacement.']]},
    biological:{title:'Cellular exposure memory',proposal:'A short exposure changes which repressor remains high. Mutual repression stores the event after the exposure ends.',requirements:['Two stable concentration states','Write with a transient exposure','Retain through cell division'],source:'toggle',known:'A genetic toggle with two stable expression states was constructed experimentally in 2000.',open:'A new exposure recorder would require a specified input receptor, a measured switching dose, inheritance and leakage data. Reusing bistability alone supplies no originality.',tests:[['Fit the two-variable model','Measure production, repression and dilution rates; retain the positive concentration domain.'],['Predict a withheld dose','Derive a switching boundary for the selected input, then test a pulse whose outcome was withheld.'],['Break the mechanism','Disable one repression link and compare persistence of expression under the same exposure.']]},
    stochastic:{title:'Stochastic growth measurement',proposal:'Measure the logarithmic signal to distinguish growth drift from the bias introduced by multiplicative fluctuations.',requirements:['Positive source signal','Declared stochastic convention','Measure mean and variance separately'],source:'ito',known:'For an Itô process dX = μXdt + σXdW, the log drift is μ − σ²/2. The Stratonovich reading has log drift μ.',open:'The transformed law is established mathematics. The application question is whether a measured physical signal and its preparation satisfy this model.',tests:[['Identify the source process','Estimate drift and noise on the original positive signal, rather than fitting the transformed answer.'],['Predict a withheld interval','Calculate the log mean and variance and compare them with withheld trajectories.'],['Break the convention','Swap the convention while keeping the displayed coefficients. The mean prediction must shift; the variance rate remains σ².']]}
  };
  function renderApplication(){
    const key=$('application-target').value,t=targets[key];
    $('requirements').innerHTML=t.requirements.map((r,i)=>`<label class="requirement"><input type="checkbox" checked data-requirement="${i}">${r}</label>`).join('');
    const update=()=>{
      const chosen=[...document.querySelectorAll('[data-requirement]:checked')].map(e=>t.requirements[Number(e.dataset.requirement)]);
      if(key==='stochastic'){
        const s={...states.stochastic,map:'log',correction:true},ito=P.stochastic({...s,convention:'ito'}),strat=P.stochastic({...s,convention:'stratonovich'});
        $('application-model').innerHTML=`<strong>Log-coordinate consequence</strong><br>At μ = ${fmt(s.mu)}, σ = ${fmt(s.noise?s.sigma:0)}:<br>drift = ${fmt(ito.exact)} (Itô)<br>drift = ${fmt(strat.exact)} (Stratonovich)`;
      }else{
        const m=P.memory(states.memory);
        $('application-model').innerHTML=`<strong>Current memory construction</strong><br>${m.stable.length} attracting state(s)<br>|h꜀| = ${fmt(m.threshold)}<br>ΔV = ${fmt(m.barrier)}<br><br>${m.bistable&&states.memory.field?'Bistability and a writing input are attached.':'A bistable, writable starting model is still required.'}`;
      }
      $('application-result').innerHTML=`<p class="proposed">${t.proposal}</p><div class="prior-item"><h3>Already demonstrated</h3><p>${t.known}</p><a href="${config.sources[t.source].url}" target="_blank" rel="noopener">${config.sources[t.source].title}${icon('arrow-up-right')}</a></div><div class="prior-item"><h3>Application question</h3><p>${t.open}</p></div><div class="prior-item"><h3>Selected requirements</h3><p>${chosen.length?chosen.join(' · '):'No requirements selected.'}</p><p style="margin-top:9px">Comparison covers the cited examples. Establishing originality requires a wider literature search for this specific combination.</p></div>`;
      $('application-experiment').innerHTML=t.tests.map(([a,b],i)=>`<div class="test-step"><b>${i+1}. ${a}</b><p>${b}</p></div>`).join('');refreshIcons();
    };
    document.querySelectorAll('[data-requirement]').forEach(e=>e.addEventListener('change',update));update();
  }
  function download(payload,name){const blob=new Blob([JSON.stringify(payload,null,2)+'\n'],{type:'application/json'}),url=URL.createObjectURL(blob),a=document.createElement('a');a.href=url;a.download=name;a.click();setTimeout(()=>URL.revokeObjectURL(url),1000);}
  $('reset').addEventListener('click',()=>{state=structuredClone(defaults[lab]);states[lab]=state;renderControls();calculate();});
  document.querySelectorAll('[data-mobile]').forEach(b=>b.addEventListener('click',()=>{document.querySelector('.work-grid').dataset.mobileView=b.dataset.mobile;document.querySelectorAll('[data-mobile]').forEach(e=>e.setAttribute('aria-pressed',String(e===b)));draw();}));
  $('export').addEventListener('click',()=>{download(lab==='stochastic'?P.constructionSpec(state):{schema:config.schema,model:'overdamped-quartic',parameters:state,calculated:result,source:config.sources.kramers.url},'fieldbridge-'+lab+'.json');});
  $('application-export').addEventListener('click',()=>{const t=targets[$('application-target').value];download({schema:'fieldbridge-application-question/1',target:t.title,construction:t.proposal,requirements:[...document.querySelectorAll('[data-requirement]:checked')].map(e=>t.requirements[Number(e.dataset.requirement)]),prior_work:config.sources[t.source],open_question:t.open,tests:t.tests,novelty:'not established'},'fieldbridge-application.json');});
  document.querySelectorAll('[data-lab]').forEach(b=>b.addEventListener('click',()=>switchLab(b.dataset.lab)));
  $('application-target').addEventListener('change',renderApplication);
  document.querySelectorAll('[data-plot]').forEach(b=>b.addEventListener('click',()=>{setPlotMode(b.dataset.plot);memoryOutcome();draw();}));
  $('sources-content').innerHTML=Object.values(config.sources).map(s=>`<article class="source-record"><a href="${s.url}" target="_blank" rel="noopener">${s.title}</a><small>${s.authors} · ${s.journal} (${s.year})</small><p>${s.supports}</p></article>`).join('')+'<article class="source-record"><a href="verified_examples.json" target="_blank" rel="noopener">Symbolic calculation records</a><p>Source specifications, generator identities and omission controls generated by the Python verifier.</p></article>';
  $('sources-open').addEventListener('click',()=>$('sources-dialog').showModal());$('sources-close').addEventListener('click',()=>$('sources-dialog').close());
  document.querySelector('[role=tablist]').addEventListener('keydown',e=>{if(!['ArrowLeft','ArrowRight','Home','End'].includes(e.key))return;e.preventDefault();const tabs=[...document.querySelectorAll('[data-lab]')];let i=tabs.findIndex(t=>t.dataset.lab===lab);i=e.key==='Home'?0:e.key==='End'?tabs.length-1:(i+(e.key==='ArrowRight'?1:-1)+tabs.length)%tabs.length;switchLab(tabs[i].dataset.lab);tabs[i].focus();});
  window.addEventListener('hashchange',()=>{const key=location.hash.slice(1);if(['memory','stochastic','applications','gallery','spins'].includes(key))switchLab(key);});
  new ResizeObserver(()=>draw()).observe($('plot'));
  setPlotMode('trajectory');renderControls();calculate();refreshIcons();
  const initial=location.hash.slice(1);if(['stochastic','applications','gallery','spins'].includes(initial))switchLab(initial);
})();
