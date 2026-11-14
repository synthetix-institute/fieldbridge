# Module 18: Two copies that exchange one variable

Two identical copies of a material that exchange one variable can become different from each other: two reactors
that share a reagent, two oscillators coupled through one coordinate. Turing's condition says when their equal state
is lost. This module calculates what holds the two apart afterwards.

When the copies differ by `d = (u_a - u_b)/2`, their shared state `m = (u_a + u_b)/2` moves by an amount proportional
to `d^2`, the same whichever copy is high, and the moved shared state acts back on the difference. The amplitude of
the difference obeys `dA/dt = lam A + g A^3` with two parts in `g`:

- `g_own`, the copy's own cubic term;
- `g_shared`, the back-action of the shared state, `-w . H(J^{-1} H(v, v), v)/2`.

With `g < 0` the copies settle apart at `|u_a - u_b| = 2 (lam/|g|)^(1/2)`. Holding the shared state fixed removes
`g_shared`: the held difference is `2 (lam/|g_own|)^(1/2)`, or there is none near the equal state when `g_own > 0`.

```bash
python3 -B -m fieldbridge memory pair examples/memory/shared_state/second_variable_limit.json \
    --exchange y --rate 3 --control 1 --out-dir /tmp/pair_limit
python3 -B -m fieldbridge memory pair examples/memory/shared_state/fitzhugh_nagumo_units.json \
    --exchange v --rate 3 --control 0.6125 --out-dir /tmp/pair_fhn
```

The first unit has no cubic term: its only nonlinearity is a quadratic coupling to a relaxing variable `s`. The
whole limit comes from the shared state (`g_own = 0`, `g_shared = -1.00`); the free pair settles apart at 1.458
(predicted 1.455), and with the shared state held the copies run apart.

In two FitzHugh-Nagumo units exchanging their inhibitor the units' own cubic term holds the difference
(`g_own = -0.90`) and the shared state pushes it outward (`g_shared = +0.38`). Holding the shared state lowers the
difference from 1.63 to 1.56 (predicted 1.55). The free prediction, 2.03, is less accurate here, because the
difference grows (0.54) almost as fast as the shared state relaxes (0.63).

The calculation is local. It does not see an unequal state far from the equal one, an oscillating difference, or how
long an unequal state survives noise; the integration of the pair, which the command also runs, does.
