# Module 16. Chaos and blow-up

**Learning objectives.** After this module you can

1. tell apart the three things a material without a stable state can do: settle on a limit cycle, settle on a
   chaotic attractor, or reach infinity in finite time;
2. measure the largest Lyapunov exponent of a chaotic material and test the predictability law;
3. find the exponents of a blow-up from the balance of the dominant terms of the drift, and see when the degree of the
   drift gives them and when it does not.

**Prerequisites.** [Module 4](18_memory_writing_and_retention.md) for the three laws by which a stored state is lost,
and [Module 6](20_memory_phase.md) for a material whose only memory is the phase of a limit cycle.
**Time.** About 30 minutes of reading; each example runs in one to three minutes.

## 1. Two outcomes without a stable state

Module 4 classifies a material by its stable states: where one is written, and how it is lost by relaxation, by
diffusion along a neutral direction or by activation over a barrier. A material without a stable state keeps nothing
in this sense. Module 6 treats the case in which its preparations settle on a limit cycle, whose phase is a neutral
direction that a pulse can shift and a periodic drive can fix. Two other outcomes are common in published models and
have laws of their own.

On a **chaotic attractor** the preparations settle on a bounded set on which neighbouring states separate
exponentially, as $\delta\,e^{\lambda_1 t}$, with $\lambda_1 > 0$ the largest Lyapunov exponent. A preparation
known to the resolution $\delta$ can be told from its neighbours only until the separation reaches the size $\Delta$
of the attractor, after the predictability time

$$T(\delta) = \frac{1}{\lambda_1}\ln\frac{\Delta}{\delta} \qquad \text{(Lorenz, 1969)}.$$

The initial state is lost at $\lambda_1/\ln 2$ bits per unit time, without noise. This is not one of the three laws
of loss: those describe a stored state under noise, whereas here there is no stored state and the loss comes from the
stretching of the deterministic flow.

In a **blow-up** every start reaches infinity at a finite time $t^*$. Near $t^*$ each growing variable follows a power
of the remaining time, $x_j \approx c_j (t^*-t)^{-p_j}$. What the material keeps of its start is the time $t^*$, and
the exponents $p_j$ follow from the drift alone.

`memory dynamics` decides which of the two applies and measures its law:

```bash
python3 -m fieldbridge memory dynamics examples/dynamics/lorenz.json --out-dir results/lorenz
python3 -m fieldbridge memory dynamics examples/dynamics/cubic_fall.json --out-dir results/cubic_fall
```

The results of the four examples of this module are saved in
[examples/dynamics/results](../../examples/dynamics/results). On the web page, these examples, and materials
without a stable state whose oscillation does not lock to a periodic drive, are classed by the same tests.

## 2. Chaotic loss of the initial state

**The test.** Starts drawn from the carrier are integrated through a transient. The largest Lyapunov exponent of
each is measured from a neighbour at distance $10^{-8}$, renormalized every ten steps (Benettin et al., 1980). The
motion is called chaotic when $\lambda_1$ is positive at every start, its spread between starts is below half its
value, its value over the first half of the run agrees with the full run to within a quarter, and the run spans at
least ten e-foldings. The last two conditions exclude neutral cycles: neighbouring closed orbits of different periods
separate linearly in time, so that their finite-time exponent, $\ln t/t$, is positive but halves when the run
doubles. The Lotka–Volterra cycles of Module 11 give $\lambda_1 \approx 0.02$ over a run of 300 time units and are
not called chaotic.

The predictability law is then measured directly. From 64 points spread along the attractor, pairs separated by
$\delta = 10^{-9}, 10^{-7}, 10^{-5}$ and $10^{-3}$ of its size are integrated until they are a tenth of its size
apart. The slope of the mean separation time against $\ln(1/\delta)$, multiplied by $\lambda_1$, must be 1.

**The Lorenz model of convection** (`examples/dynamics/lorenz.json`). Lorenz (1963) truncated the equations of a
fluid layer heated from below to three modes: the amplitude $x$ of a convection roll and two modes $y$, $z$ of the
temperature. At $\sigma = 10$, $b = 8/3$ and Rayleigh number $r = 28$ the model has no stable state, and

| | $\lambda_1$ | bits per unit time | slope $\times\,\lambda_1$ |
| --- | --- | --- | --- |
| Lorenz model, $r = 28$ | 0.906 | 1.31 | 0.98 |
| laser above its second threshold | 0.481 | 0.69 | 0.99 |

The value $\lambda_1 = 0.906$ found by Sprott (2003) for these parameters is recovered. On the web page, lowering $r$
to 0.5, below the onset of convection, leaves one state, the layer at rest, which every preparation reaches: the same
equations then keep nothing of a preparation, for a different reason.

**A single-mode laser** (`examples/dynamics/laser_lorenz_haken.json`). Haken (1975) showed that a homogeneously
broadened single-mode laser obeys the same equations: the field $E$, the polarization $P$ and the inversion $D$ of the
Maxwell–Bloch equations,

$$\dot E = k(P - E), \qquad \dot P = ED - P, \qquad \dot D = g(r - D) - EP,$$

with time in units of the polarization decay time, become the Lorenz variables with $x = E$, $y = P$, $z = r - D$,
$\sigma = k$ and $b = g$. The laser is chaotic only in a bad cavity, $k > g + 1$, and above a second threshold of the
pump, $r_c = k(k+g+3)/(k-g-1)$; with $k = 3$ and $g = 1$ this is 21 times the first lasing threshold, and the example
pumps at 25. The page joins the two realizations as one mechanism in two fields; the table gives the laser's
exponent in units of the polarization decay rate.

## 3. Blow-up in finite time

**The balance of the dominant terms.** Insert $x_j \approx c_j \tau^{-p_j}$, with $\tau = t^* - t$, into a
polynomial drift $\dot x_i = \sum_k a_{ik} \prod_j x_j^{m_{kj}}$. Near $t^*$ the largest powers must balance in every
equation (Kowalevski, 1889; Goriely and Hyde, 2000):

$$\max_k\, m_k\cdot p = p_i + 1, \qquad p_i\, c_i = \sum_{\text{dominant } k} a_{ik} \prod_j c_j^{m_{kj}},$$

with a real solution $c$. The state then grows with the exponent $p = \max_j p_j$, and along the escape the drift
grows as a power $n_{\mathrm{eff}} = 1 + 1/p$ of the state. A variable that does not grow as a power must have a
drift that grows at most as $\tau^{-1}$: below that it stays bounded, and at $\tau^{-1}$ it grows as $\ln(1/\tau)$,
which changes the other variables only by logarithmic factors ($\dot x = x^3$, $\dot y = x^2 - y$ gives
$x \sim \tau^{-1/2}$ and $y \sim \tfrac12 \ln(1/\tau)$).

If the escape runs along the terms of highest degree $n$, the balance gives $p = 1/(n-1)$, and the degree predicts
the exponent. It need not run along them.

**A particle falling in a cubic potential** (`examples/dynamics/cubic_fall.json`). A unit mass in $V = -x^3/3$ obeys
$\ddot x = x^2$, written as $\dot x = y$, $\dot y = x^2$. The drift has degree 2, which suggests $p = 1$. The balance
instead gives $x \approx 6\,\tau^{-2}$ and $y \approx 12\,\tau^{-3}$, because the first equation is linear and the
escape runs through it: the velocity grows with exponent 3. The measurement gives $n_{\mathrm{eff}} =$ 1.334
along the escape, so $p =$ 2.997.

**Hyperbolic growth of a population** (`examples/dynamics/hyperbolic_growth.json`). von Foerster, Mora and Amiot
(1960) fitted the world population of the preceding two millennia with a per-capita growth rate proportional to the
population, $\dot n = n^2/C$, whose solution $n = C/(t^*-t)$ reaches infinity at a finite date; their fit put $t^*$
in November 2026, which gave the paper its title. Here the degree gives the exponent: the balance is $p = 1$ with
$c = C$, and the measurement gives $p =$ 1.000. After the 1960s the world population grew more slowly than the
model, which has no limit of resources.

**Scaled starts.** The dominant part of the drift does not change under $x_j \to s^{p_j} x_j$, $t \to t/s$, so starts
scaled with the weights $p_j$ of the balance reach infinity after a time proportional to $1/s$. The slope of the
logarithm of the escape time against $\ln s$ is −1.000 for the population and −1.064 for the falling
particle, whose subdominant term still shifts the slope by a few per cent at the amplitudes of the test.

## 4. Applying it to your model

`memory dynamics SPEC` reports, for any specification of equations:

- the chaos test, with $\lambda_1$ at each start, the predictability times and the slope times $\lambda_1$;
- the blow-up test when a start escapes: the measured exponent, the balances of the drift, whether one of them
  predicts the exponent, the degree of the drift, and the escape times of scaled starts.

A growing state that does not blow up, such as $\dot x = x$, is reported as growth without blow-up; a material in
which no start escapes and the motion is not chaotic is left to Modules 4 and 6.

**Limits.** The balances are found for polynomial drifts only; a drift with other functions is measured but not
predicted. The tests use the deterministic part of the drift and leave out the noise of the specification. The chaos
test is numerical: an attractor with a small exponent needs a run long enough for ten e-foldings, and a long
transient can hide a chaotic attractor behind a slow approach. A blow-up exponent needs the escape to reach
amplitudes at which the dominant terms rule.
