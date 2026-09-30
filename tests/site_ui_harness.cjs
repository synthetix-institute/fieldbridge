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
    querySelector: () => element(), querySelectorAll: () => [], scrollIntoView: noop, setPointerCapture: noop,
    getBoundingClientRect: () => ({width: 400, height: 400, left: 0, top: 0, right: 400, bottom: 400}),
    getContext: () => ctx2d, appendChild: noop, click() { (listeners.click || []).forEach(f => f({})); },
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
  location: {search: ''}, history: {replaceState: noop},
  matchMedia: () => ({matches: false, addEventListener: noop}),
  requestAnimationFrame: () => 0, cancelAnimationFrame: noop, setTimeout: fn => { timers.push(fn); return timers.length; },
  clearTimeout: noop, getComputedStyle: () => ({getPropertyValue: () => '#000', fontFamily: 'sans-serif'}),
  ResizeObserver: class { observe() {} }, devicePixelRatio: 1,
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
                         selected: element('now-name').textContent + ' | ' + text('now-where')});
const out = [];
for (const step of JSON.parse(fs.readFileSync(scenarioFile, 'utf8'))) {
  if (step.do === 'start') I.start(step.node);
  if (step.do === 'edge') I.applyEdge(I.byId[step.edge], !!step.reverse);
  if (step.do === 'walk') I.walkTo(step.node);
  if (step.do === 'undo') I.undo();
  if (step.do === 'sequence') context.FieldBridgeSite.begin(step.id);
  if (step.do === 'next') element('sequence-next').click();
  if (step.do === 'prev') element('sequence-prev').click();
  if (step.do === 'param') { I.state.params[step.name] = step.value; I.refresh(); }
  flush();
  out.push(snapshot());
}
process.stdout.write(JSON.stringify(out));
