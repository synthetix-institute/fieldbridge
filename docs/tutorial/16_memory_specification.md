# Module 2. Describing a material as a realization

**Learning objectives.** After this module you can

1. name the six components of a realization and fill them for a material from your field;
2. write a specification as equations or as a network of coupled units;
3. read and correct the messages with which an invalid specification is refused.

**Prerequisites.** Module 1. **Time.** About 30 minutes.

## 1. Concepts

Memory is a property of a complete model of a material, not of a substance. FieldBridge writes such a model as a
realization,

$$
I_{\mathrm{real}} = \bigl((\Omega,\Xi);\,C,\,R,\,P;\,A\bigr),
$$

and a memory specification is this description in one JSON file. The same fields describe a colloidal monolayer,
a gene circuit, a chemical reactor or a flow network.

| Component | Meaning | Field in the specification |
| --- | --- | --- |
| Carrier $\Xi$ | the space of states that can hold information | `carrier`: `kind` (`euclid` for $\mathbb{R}^n$, `torus` for angles, `orthant` for concentrations), `variables`, `scale`, `period` for angles |
| Operation $\Omega$ | the deterministic dynamics | `drift` (one expression per variable) and optionally `potential`; or a `network` of units and bonds |
| Closure $C$ | what is specified or discarded to obtain closed equations: the bath, boundaries, imposed conservation, admissible states | `closure` (text) and `noise` (the intensity $D$ of the bath) |
| Observable $R$ | what is measured | `observable` (text) |
| Protocol $P$ | how the material is driven | `control`: the parameter an experiment varies, with its range |
| Parameters $A$ | the material that implements the model, entered as parameter values | `parameters`; `material` lists those an experimenter could also tune |

Every file also declares `schema: "fieldbridge-memory/1"`, a physical `question`, its `assumptions` and its
`provenance`. These are copied into every report together with a hash of the input, so that each result can be
traced to the exact specification that produced it.

## 2. Two forms of the operation

**Equations.** One drift expression per variable, written in the declared variables and parameters. The parser
accepts numbers, declared names, `+ - * / **`, the functions `sqrt exp log sin cos tanh abs` and the constant `pi`.
Text is never evaluated as Python. If a potential is given, the program verifies symbolically that it generates
the drift, $F = -\nabla V$.

**Networks.** For many units coupled in pairs, give the unit type and the bonds. Each bond carries a transport:
the map from the state of one unit to the state that the bond prefers for the other.

| Unit | Bond | Transport |
| --- | --- | --- |
| `rotor` (angles with $m$-fold symmetry) | `align`, `anti`, `reflect`, or a mixture `{"align": a, "reflect": g}` | keep the angle; turn it by $\pi/m$; reflect it across the bond axis |
| `spin` (bistable coordinates) | sign $+1$ or $-1$ | keep or invert |
| `gene` (concentrations) | regulation $+1$ (activation) or $-1$ (repression) | keep or invert, in one direction only |

[colloid_patch.json](../../examples/memory/colloid_patch.json) is a network of twelve caged anisotropic colloids
at a fluid interface. Each bond combines an alignment of strength 4 with a capillary reflection of strength
$5\,(r_0/r)^4$.

## 3. Worked example: a bistable reaction

The Schlögl reaction scheme has two stable concentrations of an autocatalytic species. Its specification,
[schlogl.json](../../examples/memory/schlogl.json), has the drift

```json
"drift": {"x": "-x**3 + a*x**2 - k3*x + b"},
"potential": "x**4/4 - a*x**3/3 + k3*x**2/2 - b*x",
"control": {"name": "b", "range": [1.075, 2.675]}
```

with the autocatalytic feed `a`, the rate constant `k3` and the feed `b` as control. The program verifies that the
potential generates the drift. Run the structural analysis of the next module on it:

```bash
python3 -B -m fieldbridge memory predict examples/memory/schlogl.json --out-dir build/schlogl
```

The report `predict.md` lists what the structure implies: a gradient drift cannot oscillate
(`oscillation: impossible`), and without a symmetry the write is a fold (`one-sided writes at folds; a symmetric
write needs one tuned parameter (a cusp)`). The feed `b` tilts the landscape, so writing by the control alone is
one-sided.

## 4. What is refused, and why

Each refusal names what to change. All of them are covered by tests.

| Input | Message | Reason |
| --- | --- | --- |
| `"u": "__import__('os')"` | Unsupported syntax | Only arithmetic is evaluated |
| an undeclared symbol | Undeclared symbol 'w' | Every name must be a variable or a parameter |
| `"u": "u.real"` | Unsupported syntax | No attribute access |
| a control that is not a parameter | must be a declared parameter | The control is varied by the program |
| a potential that does not generate the drift | drift = -grad(potential) | A gradient claim must hold |
| `kind: "field"` given to `memory card` | read by 'fieldbridge memory field' | Fields have their own command (Module 8) |

## 5. Exercises

1. Copy `schlogl.json`, replace the rate law by one from your own field, and keep `question`, `assumptions` and
   `provenance` accurate. Run `memory predict` on it.
2. Add a potential that does not generate your drift and read the refusal.
3. For a network of three rotors with bonds `reflect`, `reflect` and `align`, which transport composition decides
   whether the loop can satisfy all bonds?

<details><summary>Answers</summary>

1. The report lists the structural predictions for your drift; each names its rule and the result that would
   falsify it.
2. The refusal reads `drift = -grad(potential)`; correct the potential or remove it.
3. The composition of two reflections and one identity around the loop: two reflections compose to a rotation,
   so the loop is satisfiable only if that rotation is a multiple of $2\pi/m$ (Module 3).

</details>

## Summary

- A specification is a realization: carrier, operation, closure, observable, protocol and parameters.
- The operation is given as equations or as a network of units with bond transports.
- Uncertain modelling choices belong in `assumptions`; the program records them but does not verify them.

## Reference

Code: [`spec.py`](../../fieldbridge/memory/spec.py). Tests: `test_specifications_outside_the_contract_are_refused`,
`test_a_potential_must_generate_the_drift`, `test_every_example_loads_with_its_provenance`.

[Previous: Module 1](15_memory_first_card.md) · [Next: Module 3, predictions from structure](17_memory_predictions.md) · [Tutorial index](index.md) · [Glossary](memory_glossary.md)
