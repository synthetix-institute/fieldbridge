# What a fingerprint can recognize

A probability-density equation and a noisy update law may both contain drift
and diffusion. Before comparing them mathematically, it helps to find those
terms in the equations and nearby prose. FieldBridge uses a small collection
of matching rules for this first selection.

```bash
python3 -B -m fieldbridge fingerprint examples/brownian_probability_flow.tex
```

The report displays the scores and matched expressions. Six routes distinguish
transport, spectral operators, constraints, boundary or weak formulations,
ordered procedures, and incompatibility. Five evidence fibers describe
additional kinds of mathematical or physical context. Their code definitions
are in [routes.py](../../fieldbridge/routes.py).

These scores answer “which cues occur here?” rather than “which equation has
been proved?”. A mention of normalization can activate a closure cue even when
the paper has not established normalization. This is why the matching text is
useful to inspect alongside the score.

## Remove a condition and compare

Run both inputs:

```bash
python3 -B -m fieldbridge fingerprint \
  'partial_t q + nabla dot J = 0; boundary flux; observable'
python3 -B -m fieldbridge fingerprint \
  'partial_t q + nabla dot J = 0'
```

The second input removes explicit boundary and measurement words. Compare the
matched expressions, not only the largest score: rules may recognize several
cues in the same notation, and normalization of scores can change more than
one coordinate.

Physically, the continuity equation states local balance. Conservation of the
total amount also requires the boundary flux. Removing the boundary sentence
from a document removes evidence about that condition; it does not itself
change the physical boundary. The fingerprint detects the textual change,
while a later mathematical specification must supply the actual condition.

## From cues to retrieval

[Fingerprint.vector](../../fieldbridge/models.py) fixes the order of numerical
coordinates. [score_record](../../fieldbridge/search.py) compares this vector
with stored records, then adds route, fiber and keyword contributions. A record
with different notation can therefore be retrieved through common cues, but
a change of notation can also defeat a matching rule.

**Exercise.** Replace a familiar word by a mathematically equivalent expression.
Does its score remain the same? A change reveals a recognizer limitation,
not a difference between physical mechanisms. Keep the two inputs as a
regression test when improving the rules.

[Next: inspect the mechanism sheet](02_mechanism_sheets.md) · [Tutorial](index.md)
