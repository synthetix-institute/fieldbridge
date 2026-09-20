# Recover the equation and its physical conditions

A retrieved diffusion equation becomes a usable source only when its
variables and conditions are known. Is its density defined in physical space
or configuration space? Is the boundary reflecting or absorbing? Which
quantity will be measured? The mechanism sheet puts possible answers next to
the extracted equation so they can be checked against the document.

```bash
python3 -B -m fieldbridge extract examples/brownian_probability_flow.tex \
  --title 'Brownian probability flow'
```

The source is a small readable TeX example. Open
[brownian_probability_flow.tex](../../examples/brownian_probability_flow.tex)
and compare it with the output: which equation is present verbatim, which
state or boundary has been inferred, and which controls were suggested?

## What extraction actually does

[extract_mechanism](../../fieldbridge/extract.py) first collects
equation-bearing fragments and computes the fingerprint.
`guess_roles` then chooses state, input, boundary and output descriptions
through rules. Measurements and controls include generic suggestions.

There is an important fallback: when `extract_equations` finds no equation,
`default_equations` supplies a route-dependent example. That formula is not
a recovered source equation. Inspect the input before treating a sheet as
evidence, especially when it was built from sparse prose or a difficult PDF.

## Give the equation a mathematical domain

The tutorials use the notation

```math
M=(\Omega,\Xi),\qquad I_{\mathrm{op}}=(M;C,R,P).
```

For diffusion, $\Omega$ describes a diffusion operation and $\Xi$ the
space of admissible fields. Conditions $C$ specify the constitutive law and
boundary domain. The observable $R$ might be total mass or concentration
at a detector; $P$ gives preparation and the sequence of applied changes.

A reflecting and an absorbing boundary can share the bulk expression
$D\nabla^2$ while defining different generators. Their total-mass predictions
differ because the boundary flux differs. Retaining a stored bulk formula is
therefore weaker than preserving the physical evolution.

## Compare two sheets

```bash
python3 -B -m fieldbridge compare \
  examples/brownian_probability_flow.tex examples/material_memory.tex
```

The comparison locates shared and changed cues. To turn one of these suggestions
into a calculated relation, recover the equations and assumptions and then use
the [mathematical specification](09_equations_and_assumptions.md).
The exact verifier never treats a guessed role or a default equation as proof
that the original paper contained it.

[Next: retrieve related records](03_cross_field_retrieval.md) · [Tutorial](index.md)
