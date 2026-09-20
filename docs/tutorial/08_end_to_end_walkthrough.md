# Your first calculated construction

We will predict a signal proportional to the square of a positive noisy
coordinate. The supplied model gives the coordinate's drift and Brownian
noise. The question is what equation governs the squared signal.

This walkthrough assumes you installed the `construction` extra as described
in the [tutorial index](index.md). All commands run from the repository root.

## 1. Read the input before running it

Open [ito_square.json](../../examples/construction/ito_square.json). Its
mathematical fields are:

```json
{
  "convention": "ito",
  "domain": "positive",
  "target_domain": "positive",
  "drift": "(theta + 1/2)/x - alpha*x",
  "noise": "1",
  "state_map": "x**2",
  "inverse_map": "sqrt(y)",
  "candidate": {
    "drift": "2*theta + 1 - 2*alpha*y",
    "variance": "4*y"
  }
}
```

The full file also declares positive coordinates and parameters, dimensionless
units, assumptions and provenance. The candidate is intentionally incomplete.
It is an equation to test, not an answer used to derive the target.

## 2. Derive and inspect

```bash
python3 -B -m fieldbridge verify-construction \
  examples/construction/ito_square.json --out-dir build/ito_square
python3 -m json.tool build/ito_square/calculation.json
```

The status should be `verified_local_generator_identity`. In `target`,
the drift is equivalent to `2*theta + 2 - 2*alpha*y` and the quadratic
variation rate (the output key `variance`) is `4*y`. This is a local noise
coefficient, not the variance of Y at a finite time. The two corrected residual coefficients are zero. The candidate
fails, leaving `residual_d_phi: "1"`.

Both statements belong together: the derived equation verifies, while the
deliberately incomplete equation is rejected. Brownian quadratic variation
supplies the missing unit of drift. [The next derivation](10_stochastic_construction.md)
shows why its coefficient is fixed.

## 3. Read the observable consequence

Under the stated moment and boundary assumptions, the mean approaches
$(\theta+1)/\alpha$. The incomplete candidate instead approaches
$(\theta+1/2)/\alpha$. Its error persists at long times; it is not a
transient discrepancy caused by the initial condition.

The report's `mean_prediction` gives the time-dependent expression. This
turns the algebraic remainder into a specific difference between predicted
signals.

## 4. Compare a map that requires no extra drift

```bash
python3 -B -m fieldbridge verify-construction \
  examples/construction/affine_transfer.json --out-dir build/affine_transfer
```

The affine map has zero second derivative. Its omission control should report
`detects_omission: false`: removing a term that is already zero changes nothing.
That is the correct control outcome, not a failed test.

## 5. Connect the calculation to retrieval

```bash
python3 -B -m fieldbridge --data-dir examples/calculated_transfer/data \
  construct 'radial Brownian Ito diffusion' --to stochastic_dynamics \
  --no-hyperion --calculate \
  --correspondence examples/calculated_transfer/squared_signal.json \
  --out-dir build/calculated_transfer
```

Open `build/calculated_transfer/transfer.md`. The source drift now comes from
the selected retrieved record; the correspondence file supplies the map.
The same verifier derives the same target. The
[adapter chapter](13_retrieval_to_calculation.md) explains the saved source
binding and how to alter the map without altering the source law.

## Troubleshooting

| Symptom | Check |
| --- | --- |
| `No module named fieldbridge` | Activate the environment and run `python3 -m pip install -e '.[construction]'` from this repository |
| SymPy import fails | Install the `construction` extra in the same Python environment |
| Input file not found | Run from the repository root, or use absolute paths |
| `unrecognized arguments: --data-dir` | Put the global option before the subcommand |
| Calculated source is refused | Read `run_status.json` and `calculation.json`; check the selected identifier, canonical annotation and map |
| Convention or domain declaration is missing | Specify `convention`, source `domain` and `target_domain`; the calculation does not assume them |
| No atlas snapshot | The examples use `--no-hyperion`; no atlas is needed |

A repeated run uses the same output directory. Choose a new directory when
comparing variants so the two inputs and results remain available.

[Next: equations and assumptions](09_equations_and_assumptions.md) · [Tutorial](index.md)
