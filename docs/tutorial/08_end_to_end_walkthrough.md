# Chapter 8: End-to-End Walkthrough

This walkthrough follows one mechanism from source equation to target-field
construction.

## 1. Inspect the source operation

```bash
fieldbridge fingerprint examples/brownian_probability_flow.tex
fieldbridge extract examples/brownian_probability_flow.tex
```

The fingerprint identifies transport and diffusion evidence. The mechanism
sheet adds the carrier, input, boundary, output, equations, measurements, and
controls needed to interpret those scores.

## 2. Find existing target receptors

```bash
fieldbridge search examples/brownian_probability_flow.tex \
  --target-field stochastic_optimization
```

`find_analogs` compares route and fiber vectors with records loaded by
`database.load_all`. Keyword overlap is reported as supporting evidence; it is
not the primary representation.

## 3. Build the target mechanism

```bash
fieldbridge construct examples/brownian_probability_flow.tex \
  --to stochastic_optimization \
  --no-hyperion
```

Read the result as a sequence of typed edits:

```text
retain Omega      gradient drift plus diffusion
replace Xi        particle position -> parameter space
attach C          normalized domain and zero-current closure
attach R          stationary density and free-energy readout
attach P          stochastic-gradient and noise protocol
test I_op         reverse drift, remove noise, shuffle gradients
```

The target equations are only one part of the result. Closure tells us when
they define an admissible model, readout connects them to observations, and the
protocol states how the mechanism is executed.

## 4. Add a new target field

```bash
fieldbridge build-field-adapter /path/to/papers \
  --field-id active_matter \
  --label "Active Matter" \
  --out-dir build/active_matter
```

Then point the constructor at the generated data tree:

```bash
fieldbridge construct examples/brownian_probability_flow.tex \
  --to active_matter \
  --data-dir build/active_matter
```

The adapter does not decide that the transfer is valid. It supplies target-side
carriers, operations, closures, readouts, protocols, falsifiers, and source
passages from which a candidate can be derived.

## 5. Evaluate before promotion

Use `validate-zero-shot` to test whether the representation retrieves matching
mechanisms from complete unseen papers. Use `validate-continuation` when the
claim concerns the next mathematical move or destination state. A proposed
construction still requires dimensional analysis, closure and residual checks,
negative controls, and independent target-system evidence.

You have now traversed the complete public code path. Return to the
[tutorial index](index.md) or inspect `fieldbridge/constructor.py` to extend the
constructor operations.
