# Evaluate whether retrieval finds the right mechanism

A retrieval system can look successful when a fragment of the query paper
remains in its gallery. The complete-paper evaluation removes that paper and
asks whether the remaining candidates contain an independently labelled
counterpart in another field.

First generate the small tutorial dataset:

```bash
python3 examples/build_tutorial_validation_data.py
python3 -B -m fieldbridge validate-zero-shot build/tutorial_data/papers.json \
  --top-k 10 --min-eligible-queries 4 \
  --out-json build/full_paper_zero_shot.json \
  --out-md build/full_paper_zero_shot.md
```

Open the Markdown report. The tutorial dataset exists to exercise the software;
its small, authored collection cannot support a claim about scientific
retrieval at corpus scale.

## How relevance is defined

A manifest row supplies `paper_id`, `path`, `mechanism_id` and `field_id`.
A relevant candidate has the same mechanism label and a different field label.
The labels must therefore be assigned independently of the fingerprint being
evaluated. Defining the “correct mechanism” using the same matching rules
would make the comparison circular.

[validate_full_paper_zero_shot](../../fieldbridge/zero_shot.py) compares
operational retrieval with a lexical TF-IDF baseline. It computes ranking
metrics and a paired bootstrap interval for their difference. Read the number
of eligible queries alongside the scores: a query without a cross-field
counterpart cannot measure the intended task.

## Read the result without changing the question

Precision asks how many returned records are relevant; recall asks how many
relevant records were recovered; reciprocal rank rewards finding one early.
None of these evaluates the mathematical validity of a transfer. That requires
the equation-based tests in the construction chapters.

For a research evaluation, remove duplicate papers and near-duplicate
derivations across partitions, use independently assigned labels, and fix the
sample and thresholds before inspecting the result. The command's default
minimum is 100 eligible queries. The lower tutorial setting tests execution,
not the adequacy of a small research sample.

**Exercise.** Inspect one poorly ranked query and its top lexical and
operational matches. Is the disagreement due to a missing recognizer,
a disputed mechanism label, or genuinely different equations? Those cases
require different repairs.

[Next: evaluate continuation](06_future_state_prediction.md) · [Tutorial](index.md)
