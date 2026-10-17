# Module 13. The write law across the noise ratio

**Learning objectives.** After this module you can

1. state which number decides whether a swept switch follows the Kondepudi-Nelson law, and why;
2. read the constant of the write law, $c(\Gamma)$, from its small-noise value $\pi^{1/4}$ to its activated limit;
3. explain why a first-order switch, which jumps instead of growing continuously, obeys the same law, and where it
   stops doing so.

**Prerequisites.** [Module 4](18_memory_writing_and_retention.md) and [Module 9](23_memory_codiscovery.md).
**Time.** About 30 minutes.

## 1. The question

A switch swept slowly through a pitchfork chooses the state favoured by a weak bias $h$ with probability

$$
P = \Phi\!\left(\frac{c\,h}{D^{1/2} r^{1/4}}\right), \qquad (1)
$$

where $D$ is the noise along the critical mode and $r$ the rate at which its growth rate rises. Kondepudi and Nelson
derived $c = \pi^{1/4}$. Module 4 states the condition: the barrier behind the choice must grow faster than the noise
can act, so that $\Gamma = D_s|a_3|/r$ is small. This module asks what $c$ is when that condition fails, and what happens to a switch
whose transition is first order, which the symmetric-write derivation of Module 9 stops at.

## 2. One number decides the law

In the scaled variables $\tau = r^{1/2} t$, $\xi = x|a_3|^{1/2}/r^{1/4}$ and $\eta = h|a_3|^{1/2}/r^{3/4}$ the
normal form reads

$$
\mathrm d\xi = (\tau\xi \mp \xi^3 + \eta)\,\mathrm d\tau + (2\Gamma)^{1/2}\,\mathrm dW,
\qquad \Gamma = \frac{D|a_3|}{r}, \qquad (2)
$$

with the upper sign for a continuous (supercritical) and the lower for a first-order (subcritical) switch. The
outcome depends on $\Gamma$ and on the bias alone, so $c$ is a function $c(\Gamma)$. The barrier between the two
states grows as $\varepsilon^2/4|a_3|$ while the growth rate $\varepsilon$ rises at the rate $r$; $\Gamma$ compares
the noise with this barrier at the moment the choice is made.

- For $\Gamma \ll 1$ the barrier outgrows the noise before the state can cross it. The choice is made in the linear
  stage at the bifurcation, and $c = \pi^{1/4}$.
- For $\Gamma \gg 1$ the state crosses back and forth between the wells after the bifurcation. It follows the
  thermal equilibrium of the tilted double well until the barrier freezes it, and the outcome is a Boltzmann ratio.
  Kramers rates at that moment give the activated limit
  $c \simeq \pi^{1/2}[\ln(1.6\,\Gamma)]^{1/4}\,\Gamma^{-1/4}$, and a slower sweep improves the accuracy only as
  $[\ln(1/r)]^{1/4}$.

```bash
python3 -m fieldbridge memory noise-ratio --law --out-dir results/noise_ratio/law
```

| $\Gamma$ | 0.1 | 1 | 3 | 10 | 30 | 100 |
| --- | --- | --- | --- | --- | --- | --- |
| $c$, continuous | 1.331 | 1.331 | 1.268 | 1.144 | 0.949 | 0.784 |
| $c$, first order | 1.331 | 1.331 | 1.250 | 1.110 | 0.962 | 0.791 |
| regime | linear stage | linear stage | crossover | crossover | activated | activated |

The values were simulated on Eq. (2), with the bias toward each side in turn: 40 000 sweeps per side up to
$\Gamma = 3$ and 20 000 above, with standard errors of 0.012 to 0.016. The activated limit lies about 10% above them
at $\Gamma = 30$ to 100.

## 3. A first-order switch decides before its bifurcation

The derivation of Module 9 stops at two equal tubes of an adaptive transport network because their cubic term is
positive: the symmetric state loses stability in a first-order transition, and the state jumps to a distant branch.
Before the bifurcation the symmetric state is metastable behind a barrier $\varepsilon^2/4a_3$ that shrinks to zero.
At small $\Gamma$ the state stays until the bifurcation and the bias acts in the linear stage, as in a continuous
switch. At larger $\Gamma$ it escapes over the shrinking barrier before the bifurcation: 3% of the sweeps at
$\Gamma = 0.1$, 83% at $\Gamma = 1$ and all of them from $\Gamma = 3$. The bias then acts through the Arrhenius ratio
of the two escape routes, and $c(\Gamma)$ is the same as for a continuous switch within 4%.

A first-order switch has a second scale, its bistable range: the range of the control in which the far states exist,
which ends at a fold. The law holds while the choice is made inside it. When stronger noise carries the choice toward
the fold, the far states decide: shallow far states let the state hop and lose the bias, deep ones hold it.

## 4. Running it on a material

```bash
python3 -m fieldbridge memory noise-ratio examples/memory/pitchfork.json --gammas 0.05,1,10 \
    --out-dir results/noise_ratio/pitchfork
python3 -m fieldbridge memory noise-ratio examples/noise_ratio/first_order_landau.json --gammas 1 \
    --out-dir results/noise_ratio/first_order_landau
```

The command finds the symmetric critical point along the control, continuous or first order, and sweeps the full
model through it at each noise ratio, with the bias toward each side in turn. The product $D_s|a_3|$ is fixed (option
`--q`, default $5\times10^{-4}$), so that the choice is made near the bifurcation, and the rate is
$r = D_s|a_3|/\Gamma$. The report gives the accuracy, the constant $c$ with its standard error, the expected
$c(\Gamma)$, the offset of the sweep (half the difference of the two one-sided accuracies) and, for a first-order
switch, the share of sweeps decided before the bifurcation. It also states $\Gamma$ for the material's own noise at a
sweep rate of 0.02 and the rate below which the switch leaves the linear stage. The cost grows as
$\Gamma/q^{1/2}$: $\Gamma = 10$ takes $10^5$ to $10^6$ steps.

Example results (1 000 sweeps per side): the Landau model at $\Gamma = 0.05$ gives $c = 1.33$; the first-order
Landau model at $\Gamma = 1$ gives $c = 1.27 \pm 0.08$ with 84% of the sweeps decided before the bifurcation
(normal form: 83%).

## 5. What is and is not new

Kondepudi stated the weak-noise condition of the law, the parameter that separates weak from strong noise in a swept
pitchfork is known (Berglund and Gentz 2002), and the freezing of populations at a Boltzmann ratio is the blocking of
fine magnetic grains (Néel 1949). The function $c(\Gamma)$, its activated limit and the equivalence of first-order
switches were not found in the literature checked. Tests: `tests/test_memory_noise_ratio.py`.
