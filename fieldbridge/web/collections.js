(function () {
  'use strict';
  const raw=document.getElementById('studio-config').textContent;
  if(raw.trim().startsWith('__'))return;
  const config=JSON.parse(raw),P=window.FieldBridgePhysics;
  const $=id=>document.getElementById(id),fmt=(x,n=3)=>Number(x.toFixed(n)).toString();
  const esc=s=>String(s).replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
  const link=(path,label)=>`<a href="${config.repo}/blob/main/${path}" target="_blank" rel="noopener">${label} ↗</a>`;
  if(!config.models&&$('material-select')){
  let selected=config.gallery[0].id;
  $('material-select').innerHTML=config.gallery.map(c=>`<option value="${c.id}">${esc(c.name)}</option>`).join('');
  $('material-select').addEventListener('change',e=>{selected=e.target.value;materialList();materialDetail();});
  function materialList(){
    const query=$('material-search').value.toLowerCase();
    const matches=config.gallery.filter(c=>(c.name+' '+c.tags.join(' ')).toLowerCase().includes(query));
    $('material-list').innerHTML=matches.map(c=>`<button data-material="${c.id}" aria-pressed="${c.id===selected}"><strong>${esc(c.name)}</strong><small>${esc(c.verdict.stores)}</small></button>`).join('')||'<p class="empty-state">No matching realization.</p>';
    document.querySelectorAll('[data-material]').forEach(b=>b.addEventListener('click',()=>{selected=b.dataset.material;materialList();materialDetail();}));
  }
  function materialDetail(){
    const c=config.gallery.find(c=>c.id===selected);
    $('material-select').value=selected;
    const names={Omega:'Ω · operation',Xi:'Ξ · carrier',C:'C · closure',R:'R · observable',P:'P · preparation',A:'A · realization'};
    const labels={stores:'Stored states',writes:'Writing',constructs:'Construction',holds:'Retention',lock:'Writing trials'};
    $('material-detail').innerHTML=`<div class="surface-heading"><h2>${esc(c.name)}</h2><span>Saved calculation</span></div><div class="material-actions">${selected==='card04'?'<button class="primary" id="open-scalar"><i data-lucide="sliders-horizontal"></i>Open memory constructor</button>':''}${link('docs/tutorial/21_memory_new_material.md#3-the-gallery-of-realizations','Run this model in Python')}</div><div class="material-layout"><figure><img src="${c.image}" alt="Calculated states, writing and retention for ${esc(c.name)}"><figcaption>${Object.entries(c.params).map(([k,v])=>esc(k)+' = '+esc(v)).join(' · ')}</figcaption></figure><div class="material-findings">${['stores','writes','holds'].map(k=>`<div><h3>${labels[k]}</h3><p>${esc(c.verdict[k])}</p></div>`).join('')}</div></div><details class="material-roles"><summary>Mechanism, preparation and physical realization</summary><dl>${Object.entries(c.slots).map(([k,v])=>`<dt>${names[k]||esc(k)}</dt><dd>${esc(v)}</dd>`).join('')}</dl>${['constructs','lock'].map(k=>`<h3>${labels[k]}</h3><p>${esc(c.verdict[k])}</p>`).join('')}</details><footer class="evidence-strip"><div><span>${c.source_url?`<a href="${c.source_url}" target="_blank" rel="noopener">${esc(c.provenance)}</a>`:esc(c.provenance)}<br>${link(c.specification,'Source model')} · <a href="${c.record}" download>Calculation record ↓</a></span></div>${link('docs/tutorial/23_memory_codiscovery.md','Shared mechanism derivation')}</footer>`;
    if($('open-scalar'))$('open-scalar').addEventListener('click',()=>{$('tab-memory').click();$('reset').click();});
    window.lucide.createIcons();
  }
  $('material-search').addEventListener('input',materialList);
  materialList();materialDetail();
  }

  const labels={chain:'Four-spin exchange chain',qubit:'One spin-½',spin:'One spin-3/2',bosons:'Four bosons in two wells','fermion-pair':'Cooper-pair level',collective:'Three collective spins','correlated-pair':'Two correlated spins'};
  const dimensions={chain:'4 states in the one-flip sector',qubit:'2 states',spin:'4 states',bosons:'5 states at fixed particle number','fermion-pair':'Pair sector and two unchanged single-particle states',collective:'8 states', 'correlated-pair':'4 states; two spin-½ representations'};
  let state={coupling:true,field:true,g:1,h:.5},computed={...state},rotation=P.spin(computed),phase='source',target='correlated-pair';
  function spinControls(){
    $('spin-parts').innerHTML=[['coupling','Ising coupling','g Z₀Z₁'],['field','Transverse field','h X₀']].map(([k,n,e])=>`<div class="part-row"><label for="spin-${k}">${n}<small>${e}</small></label><label class="switch"><input type="checkbox" id="spin-${k}" ${state[k]?'checked':''} aria-label="Attach ${n.toLowerCase()}"><span></span></label></div>`).join('');
    $('spin-parameters').innerHTML=[['g','Coupling g',.1,2,.1],['h','Field h',0,1.5,.05]].map(([k,n,min,max,step])=>`<div class="parameter"><div class="parameter-top"><label for="spin-${k}">${n}</label><output id="spin-value-${k}">${fmt(state[k])}</output></div><input id="spin-${k}" type="range" min="${min}" max="${max}" step="${step}" value="${state[k]}"></div>`).join('');
    ['coupling','field'].forEach(k=>$('spin-'+k).addEventListener('change',e=>{state[k]=e.target.checked;updateSpin();}));
    ['g','h'].forEach(k=>$('spin-'+k).addEventListener('input',e=>{state[k]=Number(e.target.value);$('spin-value-'+k).textContent=fmt(state[k]);updateSpin();}));
  }
  function updateSpin(){computed={...state};rotation=P.spin(computed);spinRender();}
  function carrierDiagram(carrier){
    let count=carrier==='chain'?4:carrier==='collective'?3:carrier==='correlated-pair'||carrier==='bosons'||carrier==='fermion-pair'?2:1;
    const xs=Array.from({length:count},(_,i)=>count===1?160:50+220*i/(count-1));
    const lines=xs.slice(1).map((x,i)=>`<line x1="${xs[i]+15}" x2="${x-15}" y1="40" y2="40" stroke="#336eb5" stroke-width="${carrier==='correlated-pair'&&!computed.coupling?0:3}"/>${carrier==='chain'?`<text x="${(xs[i]+x)/2}" y="22">${fmt(2*Math.abs(rotation.g)*Math.sqrt((i+1)*(3-i))/4)}</text>`:''}`).join('');
    const symbols=carrier==='bosons'?['a','b']:carrier==='fermion-pair'?['k','−k']:xs.map((_,i)=>'↑');
    return `<svg class="carrier-diagram" viewBox="0 0 320 76" role="img" aria-label="${labels[carrier]}">${lines}${xs.map((x,i)=>`<circle cx="${x}" cy="40" r="15" fill="#eef6f1" stroke="#287652"/><text x="${x}" y="45">${symbols[i]}</text>`).join('')}</svg>`;
  }
  function spinRender(){
    const r=rotation,attached=phase==='attached',detached=phase==='detached';
    $('spin-status').textContent=detached?'Rotation detached · choose a carrier':attached?'Target Hamiltonian constructed':'Source evolution calculated';$('spin-status').className='build-state';
    $('spin-detach').disabled=!r.g||detached;$('spin-attach').disabled=!detached;
    const carrier=attached?target:'correlated-pair';
    $('spin-graph').innerHTML=`<div class="core-group"><div class="group-title">Detached mechanism · rotation algebra su(2)</div><div class="spin-algebra"><span>X₀</span><span>Y₀Z₁</span><span>Z₀Z₁</span></div><p class="algebra-relation">J = (X₀, Y₀Z₁, Z₀Z₁)/2<br>[Jₐ, Jᵦ] = i εₐᵦ𝒸 J𝒸</p><div class="clauses"><span>ω = ${fmt(r.rate)}</span><span>θ = ${fmt(r.theta,1)}°</span></div></div><div class="attachment-path"></div><div class="realization-node ${detached?'detached-carrier':''}">${detached?'A · carrier detached':`A · ${labels[carrier]}`}<small>${detached?'Choose the realization to attach':dimensions[carrier]}</small>${detached?'':carrierDiagram(carrier)}</div>`;
    $('spin-graph').querySelector('.group-title').textContent=detached?'Detached mechanism · rotation algebra su(2)':'Rotation algebra su(2)';
    if(attached){$('spin-graph').querySelector('.spin-algebra').innerHTML='<span>J₁</span><span>J₂</span><span>J₃</span>';$('spin-graph').querySelector('.algebra-relation').innerHTML='[Jₐ, Jᵦ] = i εₐᵦ𝒸 J𝒸';}
    if(!r.g){$('spin-graph').querySelector('.group-title').textContent='Conserved X₀ · no rotation to detach';$('spin-graph').querySelector('.spin-algebra').innerHTML='<span>X₀</span>';$('spin-graph').querySelector('.algebra-relation').textContent='[H, X₀] = 0';}
    const spec=P.spinSpec(computed,carrier);
    const terms=spec.hamiltonian.filter(t=>t.coefficient!==0).map((t,i)=>`${t.coefficient<0?'− ':i?'+ ':''}${fmt(Math.abs(t.coefficient))} ${t.operator.includes('+')?'('+esc(t.operator)+')':esc(t.operator)}${t.hc?' + h.c.':''}`);
    $('spin-equation').innerHTML=detached?'H = ω(sin θ J₁ + cos θ J₃)':`H = ${terms.join(' ')||'0'}`;
    $('spin-native').textContent=detached?'The rotation rate, angle and normalized prediction are retained. The target carrier supplies its own operators.':carrier==='chain'?'Exchange bonds are proportional to √3 : 2 : √3. The one-flipped-spin sector realizes j = 3/2.':carrier==='bosons'?'Fixed total particle number N = 4; tunnelling and population imbalance realize j = 2.':carrier==='fermion-pair'?'The empty and paired states carry the pseudospin; the singly occupied states remain unchanged.':'Closed unitary dynamics. H and the measured operator determine the rotating subspace.';
    $('spin-consequence').textContent=!r.g?'Removing the interaction makes X₀ conserved: the predicted signal stays at its initial value.':detached?'The correlation Y₀Z₁ carries the change in X₀. Detaching keeps that rotation while releasing the physical carrier.':attached?`The ${labels[carrier].toLowerCase()} carries the same normalized rotation, with its own Hamiltonian and preparation.`:'The magnetization and two spin correlations close under the Hamiltonian. Their rotation predicts the transverse signal.';
    $('spin-metrics').innerHTML=[['Rotation rate ω',fmt(r.rate)],['Signal at t = π/ω',r.rate?fmt(r.offset-r.amplitude):'1'],['Signal amplitude',fmt(r.amplitude)],['Target carrier',detached?'not attached':labels[carrier]]].map(([k,v])=>`<dt>${k}</dt><dd>${v}</dd>`).join('');
    $('spin-legend').textContent=attached?'blue source · red target':'blue source · red coupling omitted';
    spinPlot();window.lucide.createIcons();
  }
  function spinPlot(){
    const canvas=$('spin-plot'),box=canvas.getBoundingClientRect();if(!box.width)return;
    const dpr=window.devicePixelRatio||1;canvas.width=box.width*dpr;canvas.height=box.height*dpr;
    const ctx=canvas.getContext('2d');ctx.scale(dpr,dpr);
    const W=box.width,H=box.height,l=38,t=30,w=W-52,h=H-66,period=rotation.rate?2*Math.PI/rotation.rate:4;
    const X=x=>l+x/(period*1.5)*w,Y=y=>t+(1.1-y)/2.2*h;
    ctx.clearRect(0,0,W,H);ctx.font='10px Arial';
    for(let i=0;i<5;i++){const y=-1+i*.5;ctx.strokeStyle='#e7ede9';ctx.beginPath();ctx.moveTo(l,Y(y));ctx.lineTo(l+w,Y(y));ctx.stroke();ctx.fillStyle='#68716e';ctx.fillText(fmt(y),8,Y(y)+3);const x=period*1.5*i/4;ctx.fillText(fmt(x,1),X(x)-7,H-18);}
    ['#336eb5','#b5483b'].forEach((color,index)=>{ctx.strokeStyle=color;ctx.lineWidth=2;ctx.setLineDash(index?[6,4]:[]);ctx.beginPath();for(let i=0;i<=240;i++){const x=period*1.5*i/240,y=index&&phase!=='attached'?1:rotation.signal(x);if(i===0)ctx.moveTo(X(x),Y(y));else ctx.lineTo(X(x),Y(y));}ctx.stroke();});
    ctx.setLineDash([]);ctx.fillStyle='#4a5b50';ctx.font='italic 12px Georgia';ctx.fillText('f(t)',8,16);ctx.fillText('t',W/2,H-2);
  }
  $('spin-build').addEventListener('click',()=>{state={coupling:true,field:true,g:1,h:.5};computed={...state};rotation=P.spin(computed);phase='source';target='correlated-pair';spinControls();spinRender();});
  $('spin-detach').addEventListener('click',()=>{phase='detached';spinRender();});
  $('spin-attach').addEventListener('click',()=>{target=$('spin-carrier').value;phase='attached';spinRender();});
  $('spin-carrier').addEventListener('change',()=>{if(phase==='attached'){phase='detached';spinRender();}});
  document.querySelectorAll('[data-spin-view]').forEach(b=>b.addEventListener('click',()=>{document.querySelector('.spin-grid').dataset.view=b.dataset.spinView;document.querySelectorAll('[data-spin-view]').forEach(e=>e.setAttribute('aria-pressed',String(e===b)));spinPlot();}));
  $('spin-export').addEventListener('click',()=>{const payload={specification:P.spinSpec(computed,phase==='attached'?target:'correlated-pair'),prediction:{rate:rotation.rate,angle_degrees:rotation.theta,offset:rotation.offset,amplitude:rotation.amplitude},tutorial:config.quantum.tutorial};const url=URL.createObjectURL(new Blob([JSON.stringify(payload,null,2)],{type:'application/json'}));const a=document.createElement('a');a.href=url;a.download='fieldbridge-spin-construction.json';a.click();setTimeout(()=>URL.revokeObjectURL(url),1000);});
  $('quantum-list').innerHTML=config.quantum.examples.rows.map(r=>`<article><strong>${esc(r.name)}</strong><span>${r.status==='reached'?`Rotation closes · ω = ${fmt(r.signature.rate)}`:esc(r.obstruction)}</span></article>`).join('')+`<p>${link('docs/tutorial/24_spin_language.md#6-where-a-derivation-stops','Derivations and stopping conditions')} · ${link('docs/tutorial/14_inverse_construction.md','Design a cancellation time')}</p>`;
  window.addEventListener('fieldbridge-module',e=>{if(e.detail==='spins')spinRender();});
  new ResizeObserver(spinPlot).observe($('spin-plot'));
  spinControls();spinRender();
})();
