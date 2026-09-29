/* Equations of the page as MathML: the arithmetic trees of the specifications and the operators of the carriers. */
(function (root, factory) {
  const api = factory();
  if (typeof module === 'object' && module.exports) module.exports = api;
  else root.FieldBridgeMath = api;
})(typeof globalThis !== 'undefined' ? globalThis : this, function () {
  'use strict';
  const GREEK = {alpha: 'α', beta: 'β', gamma: 'γ', delta: 'δ', Delta: 'Δ', eps: 'ε', epsilon: 'ε', kappa: 'κ',
                 lam: 'λ', lambda: 'λ', mu: 'μ', nu: 'ν', phi: 'φ', psi: 'ψ', sigma: 'σ', tau: 'τ', theta: 'θ',
                 omega: 'ω', Omega: 'Ω', pi: 'π', rho: 'ρ', xi: 'ξ', chi: 'χ', eta: 'η', zeta: 'ζ'};
  const SPECIAL = {rabi: ['Ω', 'R'], detuning: ['δ', ''], rf: ['Ω', 'rf'], w1: ['ω', '1'], w2: ['ω', '2'], hill: ['n', '']};
  const esc = s => String(s).replace(/[&<>"]/g, c => ({'&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;'}[c]));
  const mi = (s, normal) => `<mi${normal ? ' mathvariant="normal"' : ''}>${esc(s)}</mi>`;
  const sub = (base, index) => `<msub>${base}${/^\d+$/.test(index) ? `<mn>${index}</mn>` : mi(index, index.length > 1)}</msub>`;

  function name(n) {
    if (SPECIAL[n]) return SPECIAL[n][1] ? sub(mi(SPECIAL[n][0]), SPECIAL[n][1]) : mi(SPECIAL[n][0]);
    if (GREEK[n]) return mi(GREEK[n]);
    let m = /^([A-Za-z]+?)_?(\d+)$/.exec(n);
    if (m) return sub(GREEK[m[1]] ? mi(GREEK[m[1]]) : mi(m[1], m[1].length > 1), m[2]);
    m = /^([A-Za-z])_([A-Za-z]+)$/.exec(n);
    if (m) return sub(mi(m[1]), m[2]);
    return mi(n, n.length > 1);
  }
  function number(v) {
    const s = Number.isInteger(v) ? String(v) : String(+v.toPrecision(6));
    return `<mn>${s.replace('-', '−')}</mn>`;
  }
  const PREC = {add: 1, sub: 1, neg: 2, mul: 3, div: 3, pow: 4};
  const prec = t => PREC[t[0]] ?? (t[0] === 'number' && t[1] < 0 ? 2 : 5);
  const paren = s => `<mrow><mo>(</mo>${s}<mo>)</mo></mrow>`;
  /** An arithmetic tree ["add", a, b] ... as MathML (no <math> wrapper). */
  function tree(t) {
    const [op, a, b] = t;
    const wrap = (x, min) => prec(x) < min ? paren(tree(x)) : tree(x);
    switch (op) {
      case 'number': return number(a);
      case 'name': return name(a);
      case 'neg': return `<mrow><mo>−</mo>${wrap(a, 3)}</mrow>`;
      case 'add': return `<mrow>${tree(a)}<mo>+</mo>${wrap(b, 1)}</mrow>`;
      case 'sub': {
        if (b[0] === 'neg') return `<mrow>${tree(a)}<mo>+</mo>${wrap(b[1], 2)}</mrow>`;
        return `<mrow>${tree(a)}<mo>−</mo>${wrap(b, 2)}</mrow>`;
      }
      case 'mul': {
        const left = wrap(a, 3), right = wrap(b, 3);
        const dot = b[0] === 'number' || (a[0] === 'number' && b[0] === 'number') ? '<mo>·</mo>' : '<mo>&#x2062;</mo>';
        return `<mrow>${left}${dot}${right}</mrow>`;
      }
      case 'div': return `<mfrac>${tree(a)}${tree(b)}</mfrac>`;
      case 'pow': {
        if (b[0] === 'number' && b[1] === 0.5) return `<msqrt>${tree(a)}</msqrt>`;
        return `<msup>${prec(a) < 5 ? paren(tree(a)) : tree(a)}${tree(b)}</msup>`;
      }
      case 'sqrt': return `<msqrt>${tree(a)}</msqrt>`;
      case 'abs': return `<mrow><mo>|</mo>${tree(a)}<mo>|</mo></mrow>`;
      default: return `<mrow>${mi(op, true)}<mo>&#x2061;</mo>${paren(tree(a))}</mrow>`;
    }
  }
  /** A product of carrier operators: 'X0 X1', 'Jz Jz', 'ad b', 'cd0 cd1', 'na na'. */
  function token(t) {
    let m = /^([XYZ])(\d)$/.exec(t); if (m) return sub(mi(m[1]), m[2]);
    m = /^J([xyzpm])$/.exec(t); if (m) return sub(mi('J'), {p: '+', m: '−'}[m[1]] || m[1]);
    m = /^cd(\d)$/.exec(t); if (m) return `<msubsup>${mi('c')}<mn>${m[1]}</mn><mo>†</mo></msubsup>`;
    m = /^([cn])(\d)$/.exec(t); if (m) return sub(mi(m[1]), m[2]);
    m = /^n([a-z])$/.exec(t); if (m) return sub(mi('n'), m[1]);
    m = /^([a-z])d$/.exec(t); if (m) return `<msup>${mi(m[1])}<mo>†</mo></msup>`;
    return mi(t);
  }
  function product(p) {
    const toks = p.trim().split(/\s+/), out = [];
    for (let k = 0; k < toks.length;) {
      let n = 1; while (k + n < toks.length && toks[k + n] === toks[k]) n++;
      out.push(n > 1 ? `<msup>${token(toks[k])}<mn>${n}</mn></msup>` : token(toks[k])); k += n;
    }
    return `<mrow>${out.join('')}</mrow>`;
  }
  function operator(text) {
    const parts = String(text).split('+');
    return parts.length > 1 ? `<mrow>${parts.map(product).join('<mo>+</mo>')}</mrow>` : product(parts[0]);
  }
  /** One term of a Hamiltonian: coefficient, operator, and + h.c. */
  function term(coefTree, op, hc) {
    const sumOp = String(op).includes('+');
    let c = tree(coefTree);
    if (coefTree[0] === 'number' && coefTree[1] === 1) c = '';
    else if (coefTree[0] === 'number' && coefTree[1] === -1) c = '<mo>−</mo>';
    else if (prec(coefTree) < 3) c = paren(c);
    const o = sumOp ? paren(operator(op)) : operator(op);
    return `<mrow>${c}${o}${hc ? '<mo>+</mo><mi mathvariant="normal">h.c.</mi>' : ''}</mrow>`;
  }
  /** A coefficient written in a specification: a number, a simple fraction or a parameter name, as a tree. */
  function coefficient(c) {
    const s = String(c).trim();
    if (/^-?\d+(\.\d+)?$/.test(s)) return ['number', Number(s)];
    const m = /^(-?\d+)\/(\d+)$/.exec(s);
    if (m) return ['div', ['number', Number(m[1])], ['number', Number(m[2])]];
    return ['name', s];
  }
  const math = (inner, block) => `<math${block ? ' display="block"' : ''}>${inner}</math>`;
  return {tree, name, operator, term, coefficient, math, GREEK};
});
