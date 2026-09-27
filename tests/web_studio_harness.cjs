// Exercises the real studio event handlers. Layout and canvas are tested in-browser.
const fs = require('node:fs');
const path = require('node:path');
const vm = require('node:vm');
const root = path.resolve(__dirname, '..');
const elements = new Map();
let serial = 0;
function matches(el, selector) {
  if (selector.startsWith('.')) return el.classes.has(selector.slice(1));
  const match = /^\[([^=\]]+)(?:=([^\]]+))?\]$/.exec(selector);
  if (match) return match[2] ? el.attributes[match[1]] === match[2] : match[1] in el.attributes;
  return el.tag === selector;
}
class Element {
  constructor(tag, attributes = {}) {
    this.tag = tag; this.attributes = attributes; this.id = attributes.id || '';
    this.dataset = Object.fromEntries(Object.entries(attributes).filter(([k]) => k.startsWith('data-')).map(([k,v]) => [k.slice(5).replace(/-([a-z])/g, (_,c) => c.toUpperCase()),v]));
    this.classes = new Set((attributes.class || '').split(' ')); this.listeners = {};
    this.value = attributes.value || ''; this.checked = 'checked' in attributes; this.children = [];
    this.classList = {toggle: (name,on) => on ? this.classes.add(name) : this.classes.delete(name),add:name=>this.classes.add(name),remove:name=>this.classes.delete(name)};
    this.parentElement = {hidden:false};
    elements.set(this.id || 'node-' + serial++, this);
  }
  set innerHTML(html) {
    const remove = el => {el.children.forEach(remove); for (const [key,value] of elements) if (value === el) elements.delete(key);};
    this.children.forEach(remove); this.children = parse(html); this.html = html;
  }
  get innerHTML() {return this.html || '';}
  addEventListener(type, callback) {(this.listeners[type] ||= []).push(callback);}
  setAttribute(name,value) {this.attributes[name] = value;}
  querySelector(selector) {return this.children.find(el => matches(el,selector));}
  getBoundingClientRect() {return {width:0,height:0};}
  focus() {} scrollIntoView() {}
  fire(type) {(this.listeners[type] || []).forEach(f => f({target:this,preventDefault(){}}));}
}
function parse(html) {
  return [...html.matchAll(/<([a-z][\w-]*)\b([^>]*)>/gi)].map(([,tag,raw]) => {
    const attributes = Object.fromEntries([...raw.matchAll(/([\w-]+)(?:="([^"]*)")?/g)].map(([,key,value]) => [key,value || '']));
    return new Element(tag,attributes);
  });
}
parse(fs.readFileSync(path.join(root,'fieldbridge/web/index.html'),'utf8'));
elements.get('consequence').children.push(new Element('p'));
const config = {repo:'https://github.com/synthetix-institute/fieldbridge',schema:'test',sources:{kramers:{authors:'Kramers',url:'https://example.org',year:1940},ito:{authors:'Ito',url:'https://example.org',year:1944}},gallery:[{id:'card04',name:'pitchfork',tags:[],verdict:{stores:'two states',writes:'pulse',holds:'retention',constructs:'normal form',lock:'test'},slots:{},params:{eps:1},image:'test.png',specification:'test.json'}],quantum:{tutorial:'https://example.org',examples:{rows:[]}}};
if(process.argv[4])config.models=JSON.parse(fs.readFileSync(process.argv[4],'utf8'));
elements.get('studio-config').textContent = JSON.stringify(config);
const document = {
  getElementById: id => elements.get(id),
  querySelectorAll: selector => [...elements.values()].filter(e => matches(e,selector)),
  querySelector: selector => [...elements.values()].find(e => matches(e,selector)),
};
const window = {FieldBridgePhysics:require(path.join(root,'fieldbridge/web/physics.js')),lucide:{createIcons(){}},addEventListener(){},dispatchEvent(){}};
const location = {hash:process.argv[3] || ''};
vm.runInNewContext(fs.readFileSync(path.join(root,'fieldbridge/web/studio.js'),'utf8'),{
  document,window,history:{replaceState(_,__,hash){location.hash=hash;}},location,structuredClone,
  ResizeObserver:class {observe(){}},CustomEvent:class {},setTimeout,
});
vm.runInNewContext(fs.readFileSync(path.join(root,'fieldbridge/web/collections.js'),'utf8'),{
  document,window,ResizeObserver:class {observe(){}},setTimeout,
});
if(config.models){
  window.FieldBridgeModels=require(path.join(root,'fieldbridge/web/model-physics.js'));
  vm.runInNewContext(fs.readFileSync(path.join(root,'fieldbridge/web/model-studio.js'),'utf8'),{
    document,window,ResizeObserver:class {observe(){}},setTimeout,cancelAnimationFrame(){},requestAnimationFrame(){},
  });
}
const snapshots = [];
function snapshot() {
  return {equation:elements.get('equation').innerHTML,drive:elements.get('drive-equation').textContent,
    metrics:elements.get('metrics').innerHTML,status:elements.get('build-state').textContent,
    graph:elements.get('graph').innerHTML,consequence:elements.get('consequence').querySelector('p').textContent,
    spinEquation:elements.get('spin-equation').innerHTML,spinMetrics:elements.get('spin-metrics').innerHTML,
    spinConsequence:elements.get('spin-consequence').textContent,hash:location.hash,
    modelName:elements.get('model-name').textContent,modelEquations:elements.get('model-equations').innerHTML,
    modelFinal:elements.get('model-final').textContent,modelMetrics:elements.get('model-metrics').innerHTML,
    modelDesign:elements.get('model-design-result').textContent};
}
snapshots.push(snapshot());
for (const action of JSON.parse(process.argv[2] || '[]')) {
  const element = action.id ? elements.get(action.id) : document.querySelector(action.selector);
  if (!element) throw new Error('Missing control: ' + JSON.stringify(action));
  if ('value' in action) element.value = action.value;
  if ('checked' in action) element.checked = action.checked;
  element.fire(action.type);
  snapshots.push(snapshot());
}
process.stdout.write(JSON.stringify(snapshots));
