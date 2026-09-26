# Predict the next written step

Retrieval asks where a related equation already exists. Continuation asks which
move follows the current one in an unseen paper. FieldBridge includes a
count-based evaluation for this second question; it uses supplied mechanism
labels rather than generating new equations.

A transition record contains a current operation, carrier and completion state,
an observed move, and the withheld next move and destination:

```text
(current_omega, current_xi, current_completion, first_move)
    -> (next_move, destination_omega, destination_xi, destination_completion)
```

## Run the tutorial evaluation

```bash
python3 examples/build_tutorial_validation_data.py
python3 -B -m fieldbridge validate-continuation \
  build/tutorial_data/transitions.jsonl \
  --min-evaluation-transitions 12 \
  --out-json build/future_state_validation.json \
  --out-md build/future_state_validation.md
```

The generated records are synthetic. They check the split, counting model and
report, and do not reproduce historical large-corpus results. Refer to a
versioned corpus report for those results rather than transferring percentages
from a different codebook or split into this run.

## What information does the current state add?

[validate_future_state](../../fieldbridge/continuation.py) assigns complete
papers to fitting or evaluation by a stable hash. The full model conditions
on the current state and first move. The baseline conditions only on the first
move. Their comparison asks whether knowing the current mechanism labels adds
information beyond the usual successor of a move.

The code-length gain for a true future label $z$ compares its probabilities:

$$
\log_2 p_{\mathrm{full}}(z\mid\text{current state, first move})
-\log_2 p_{\mathrm{base}}(z\mid\text{first move}).
$$

Averaging this quantity measures how much better the full distribution predicts
the withheld labels. It can improve while top-1 accuracy declines: the model
may assign better probabilities without choosing the most frequent label
more often. Inspect both metrics and the baseline, not accuracy alone.

## What to hold fixed in a scientific comparison

Use one label system, one paper split and the same target for both models.
When codebooks are refitted, an absolute gain in bits may change because the
target distribution changed. A claim about the value of the role division
also needs alternatives of comparable capacity, including controls on the
feature assignment.

**Exercise.** Identify the destination targets separately from the next move in
the JSON report. Explain a positive code-length gain in one sentence without
claiming that the program generated a correct physical theory.

[Next: index a paper collection](07_pdf_field_adapter.md) · [Tutorial](index.md)
