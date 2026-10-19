// Runs the page scripts of fieldbridge/web against a stubbed document, for tests/test_site_ui.py: the navigation of
// the instrument (edges, sequences, the path to a realization, undo) and the consequences it reports. Drawing is
// replaced by no-op canvases; layout is checked in a browser.
const fs = require('node:fs');
const path = require('node:path');
const vm = require('node:vm');
const [dataFile, scenarioFile] = process.argv.slice(2);
const web = path.resolve(__dirname, '..', 'fieldbridge', 'web');

const noop = () => {};
const ctx2d = new Proxy({}, {get: (t, k) => (k === 'measureText' ? () => ({width: 10}) : k in t ? t[k] : noop), set: (t, k, v) => { t[k] = v; return true; }});
const elements = new Map();
function element(id) {
  if (id && elements.has(id)) return elements.get(id);
  const listeners = {};
  const el = {
    id, innerHTML: '', textContent: '', hidden: false, disabled: false, value: '', dataset: {}, style: {}, width: 400, height: 400,
    clientWidth: 900, attributes: {}, listeners,
    classList: {add: noop, remove: noop, toggle: noop, contains: () => false},
    setAttribute(k, v) { this.attributes[k] = v; }, getAttribute(k) { return this.attributes[k]; },
    addEventListener(type, fn) { (listeners[type] ||= []).push(fn); }, removeEventListener: noop,
    querySelector: () => element(), querySelectorAll: () => [], scrollIntoView: noop, setPointerCapture: noop, focus: noop,
    getBoundingClientRect: () => ({width: 400, height: 400, left: 0, top: 0, right: 400, bottom: 400}),
    getContext: () => ctx2d, appendChild: noop,
    dispatch(type, event = {}) { (listeners[type] || []).forEach(f => f(event)); },
    click() { this.dispatch('click'); },
  };
  el.parentElement = el;
  if (id) elements.set(id, el);
  return el;
}
const docListeners = {};
const document = {
  getElementById: element, querySelector: () => element(), querySelectorAll: () => [], body: element('body'),
  documentElement: element('html'), hidden: false,
  addEventListener: (type, fn) => { (docListeners[type] ||= []).push(fn); },
};
const timers = [];
const context = {
  document, console, Math, JSON, Number, String, Array, Object, Set, Map, Float64Array, Promise, Proxy, URLSearchParams, Error,
  location: {search: '', hash: ''}, history: {replaceState: noop}, addEventListener: noop,
  matchMedia: () => ({matches: false, addEventListener: noop}),
  requestAnimationFrame: () => 0, cancelAnimationFrame: noop, setTimeout: fn => { timers.push(fn); return timers.length; },
  clearTimeout: noop, getComputedStyle: () => ({getPropertyValue: () => '#000', fontFamily: 'sans-serif'}),
  ResizeObserver: class { observe() {} }, devicePixelRatio: 1, scrollY: 0, scrollTo: noop,
};
context.window = context;
context.globalThis = context;
vm.createContext(context);
const load = file => vm.runInContext(fs.readFileSync(file, 'utf8'), context, {filename: file});
load(dataFile);
for (const name of ['model-physics.js', 'unitary.js', 'dissipative.js', 'fields.js', 'mathml.js', 'views.js', 'mechanisms.js', 'expression.js', 'site.js'])
  load(path.join(web, name));
(docListeners.DOMContentLoaded || []).forEach(f => f());
const flush = () => { while (timers.length) timers.shift()(); };
flush();

const I = context.FieldBridgeInstrument, text = id => element(id).innerHTML.replace(/<[^>]+>/g, '');
const snapshot = () => ({node: I.state.node, path: I.state.path.map(p => p.node), lead: text('consequence-lead'),
                         facts: text('facts'), sequence: text('sequence-text'), count: element('path-count').textContent,
                         mechanism: element('m-name').textContent, changes: text('slots'),
                         reference: I.reference(), legend: text('view-legend'),
                         text: text('consequence-text'), derivation: text('derivation'), axis: I.axis(),
                         memory: text('memory-summary') + ' | ' + text('memory-groups'), laws: text('memory-laws'),
                         memoryHtml: element('memory-summary').innerHTML + element('memory-groups').innerHTML,
                         absent: (element('mechanisms').innerHTML.match(/class="mech absent"/g) || []).length,
                         cards: (element('mechanisms').innerHTML.match(/class="mech( absent)?"/g) || []).length,
                         module: context.FieldBridgeMechanisms.module(),
                         planned: text('planned-grid'),
                         chooser: {hidden: element('material-chooser').hidden, expanded: element('material-choose').attributes['aria-expanded'],
                           options: element('material-options').innerHTML, count: element('material-count').textContent,
                           mechanism: element('material-mechanism').value},
                         selected: element('now-name').textContent + ' | ' + text('now-where'),
                         params: {...I.state.params},
                         sequences: {expanded: element('sequence-start').attributes['aria-expanded'],
                           menuHidden: element('sequence-menu').hidden, picker: element('sequence-picker').innerHTML,
                           title: element('sequence-title').textContent, playerHidden: element('sequence-player').hidden},
                         guide: {hidden: element('guide').hidden, expanded: element('guide-toggle').attributes['aria-expanded'],
                           collection: element('guide-collection').textContent, map: element('map-collection').textContent},
                         back: {hidden: element('back-to').hidden, text: element('back-to').textContent}});
const out = [];
for (const step of JSON.parse(fs.readFileSync(scenarioFile, 'utf8'))) {
  if (step.do === 'start') I.start(step.node);
  if (step.do === 'edge') I.applyEdge(I.byId[step.edge], !!step.reverse);
  if (step.do === 'walk') I.walkTo(step.node);
  if (step.do === 'undo') I.undo();
  if (step.do === 'toolbar-undo') element('constructor-undo').click();
  if (step.do === 'reset') element('constructor-reset').click();
  if (step.do === 'choose-material') element('material-choose').click();
  if (step.do === 'material-filter') element('material-' + step.key).dispatch(step.key === 'search' ? 'input' : 'change', {target: {value: step.value}});
  if (step.do === 'material-select') element('material-options').dispatch('click', {target: {closest: () => ({dataset: {materialId: step.node}})}});
  if (step.do === 'material-escape') element('material-chooser').dispatch('keydown', {key: 'Escape'});
  if (step.do === 'material-close') element('material-close').click();
  if (step.do === 'module') context.FieldBridgeMechanisms.module(step.key);
  if (step.do === 'sequence') context.FieldBridgeSite.begin(step.id);
  if (step.do === 'next') element('sequence-next').click();
  if (step.do === 'prev') element('sequence-prev').click();
  if (step.do === 'open-sequences') element('sequence-start').click();
  if (step.do === 'leave') element('sequence-close').click();
  if (step.do === 'open') context.FieldBridgeSite.open(step.node, element('row-' + step.node));
  if (step.do === 'back') element('back-to').click();
  if (step.do === 'guide-close') element('guide-close').click();
  if (step.do === 'guide-toggle') element('guide-toggle').click();
  if (step.do === 'param') { I.state.params[step.name] = step.value; I.refresh(); }
  if (step.do === 'term') { I.state.active[step.k] = step.on; I.refresh(); }
  flush();
  out.push(snapshot());
}
process.stdout.write(JSON.stringify(out));
