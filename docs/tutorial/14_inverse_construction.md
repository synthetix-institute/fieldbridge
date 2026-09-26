# How to use FieldBridge for discovery: solve for the interactions

An end spin in an exchange-coupled chain transfers its polarization to its
neighbours and can later recover it. Changing the exchange strengths usually
changes the time dependence. Suppose the intended device instead needs a
polarization cancellation at a time set by a separate interaction, despite
uncertain exchange strengths. Which interactions can produce that separation
of time scales?

This example starts with that requirement and solves for coupling coefficients.
It then calculates the polarization, checks it by full Hamiltonian evolution,
and changes an assumption to determine why the result holds. It also imposes
a nearest-neighbour restriction and reports when the chosen construction has
no nonzero solution. The calculation uses a familiar collective-spin mechanism;
its purpose is to make a discovery search executable, not to label that known
mechanism as a new discovery.

## Run the complete example

From the FieldBridge repository root, using Python 3.10 or newer:

```bash
python3 -m pip install -e '.[construction]'
python3 -B -m fieldbridge design-spin-cancellation \
  examples/construction/spin_cancellation_design.json \
  --out-dir build/spin_cancellation
```

Open `build/spin_cancellation/design.md` for the explanation and
`build/spin_cancellation/design.json` for the coupling basis, numerical
comparisons and 129 time samples. `input.json` preserves the input. The report
records the input and implementation hashes and the numerical-library versions.
No language model, GPU or paper download is used.

The input specifies exchange energies 3 and 4, collective energy 5/4 and a
single-pair error of 1/10. It does **not** specify the relative Ising couplings
that solve the commutation equations. Those are calculated.

## From the physical requirement to equations

Let $X_j,Y_j,Z_j$ be Pauli operators on spin $j$. On an open chain of $N$
spins, the exchange Hamiltonian and the unknown Ising interaction are

$$
H_0=\sum_{j=0}^{N-2}J_j A_j,\qquad
A_j=X_jX_{j+1}+Y_jY_{j+1},\qquad
Q=\sum_{i<j}q_{ij}Z_iZ_j.
$$

The search chooses interactions for which exchange and Ising evolution can
be applied separately, exactly: $[H_0,Q]=0$. Since the exchange energies $J_j$
may change independently, this equality must hold separately for each bond,
$[A_j,Q]=0$. Cancellation between two commutators at one special choice of
$J_j$ would not establish the requested robustness.

For three spins there are three unknown coefficients. The first exchange
bond requires $q_{02}=q_{12}$; the second requires $q_{01}=q_{02}$. Their
intersection is

$$
q_{01}=q_{02}=q_{12}=\lambda.
$$

The program obtains this result from exact Pauli matrices. Each matrix entry
of each commutator gives a homogeneous linear equation in the unknown
coefficients. Stacking those equations gives a matrix whose nullspace is the
allowed interaction family. For the three-spin input the nullspace has one
basis vector, $(1,1,1)$, so its free scale is the energy $\lambda$.

This is implemented in [`solve_couplings`](../../fieldbridge/inverse_spin.py).
It solves the same problem for two to five spins; the complete command uses
three to five. The tests check exact constraint rank and nullspace at each
of those sizes. They do not infer an arbitrary-size theorem from that finite
enumeration.

## Why the resulting interaction separates the dynamics

Exchange conserves the total magnetization $M=\sum_j Z_j$. The uniform
interaction found by the solver can be written

$$
Q=\frac{\lambda}{2}(M^2-NI).
$$

It therefore commutes with exchange for any values of $J_j$. This identity
explains the numerical nullspace and proves that the uniform family works
for an arbitrary chain length. The propagator separates as
$e^{-i(H_0+Q)t/\hbar}=e^{-iH_0t/\hbar}e^{-iQt/\hbar}$.

To obtain a measured response, prepare the end spin along $+x$ and the
remaining spins maximally mixed:

$$
\rho(0)=\frac{I+X_0}{2^N}.
$$

For the end spin, the Jordan-Wigner representation has no preceding string
of spin operators. Exchange becomes single-particle hopping with matrix
$h_{j,j+1}=h_{j+1,j}=2J_j$ and zero diagonal. The collective interaction
adds a phase set by the magnetization of the other spins. Each maximally
mixed spin contributes the average of two opposite phases,
$\cos(2\lambda t/\hbar)$. Taking the trace leaves the hopping return
amplitude multiplied by these $N-1$ phase averages:

$$
\langle X_0(t)\rangle
=\operatorname{Re}[e^{-iht/\hbar}]_{00}
  \cos^{N-1}(2\lambda t/\hbar).
$$

This response formula is an explicit analytic construction encoded in the
example, not a formula invented by the nullspace solver. The separate
numerical calculation builds the full $2^N$-dimensional Hamiltonian and
evaluates the trace directly; it does not use the factorized expression.

Collective phase cancellation occurs at

$$
t_m=\frac{(2m+1)\pi\hbar}{4|\lambda|},\qquad m=0,1,\ldots.
$$

Exchange changes the response between these times, but not these zeros.
Exchange may introduce additional, earlier zeros, so the output calls $t_0$
the **first collective zero**, not necessarily the first zero of the signal.
The quantum state remains present; a zero of this observable is not erasure
or protected storage of quantum information.

## Read the calculation and the comparison

For the supplied three-spin input, the three-spin expression is

$$
\langle X_0(t)\rangle
=\frac{16+9\cos(10t)}{25}\cos^2(5t/2),\qquad \hbar=1.
$$

| Output | Expected result | What it establishes |
| --- | --- | --- |
| `construction.basis` | `[["1", "1", "1"]]` | Uniform pair coupling follows from the constraints |
| `prediction.first_collective_zero_time` | $\pi/5\simeq0.628319$ | The specified collective cancellation time |
| `checks.max_full_hamiltonian_error` | Below $10^{-10}$ | Direct evolution agrees with the response formula |
| `checks.exchange_disorder_zero_values` | Both near zero | Two different exchange realizations preserve that time |
| `checks.perturbed_signal_at_original_zero` | Approximately $0.001689$ | Changing one Ising coupling breaks this cancellation in the example |
| `nearest_neighbor_comparison.solution_dimension` | `0` | The restricted commuting construction has no nonzero solution |

The last row answers a separate design question. For three spins, restricting
Ising interactions to neighbours sets $q_{02}=0$. Combined with the solved
equalities, this forces every $q_{ij}$ to zero. Thus the nonzero collective
construction needs the end-to-end interaction. The exact solver also finds
no nonzero nearest-neighbour solution at four and five spins. This conclusion
concerns the specified pairwise, commuting construction; it does not rule out
other Hamiltonians, pulses or isolated cancellation times.

The coupling-error comparison changes only $q_{01}$ by 1/10, leaving the
preparation, observable and exchange unchanged. Increasing that error to 1/2
gives a signal near $0.0418$ at the original collective zero. For other inputs,
the control may leave a zero intact; the report records whether it actually
changes the signal rather than asserting that every perturbation must do so.

## Change the physical question

Create another input without replacing the example:

```bash
python3 - <<'PY'
import json
from pathlib import Path

spec = json.loads(Path('examples/construction/spin_cancellation_design.json').read_text())
spec['bonds'] = [3, 4, 2]
spec['question'] = 'Does the four-spin construction retain its collective zero with unequal exchange?'
out = Path('build/four_spin_design.json')
out.parent.mkdir(parents=True, exist_ok=True)
out.write_text(json.dumps(spec, indent=2))
PY
python3 -B -m fieldbridge design-spin-cancellation \
  build/four_spin_design.json --out-dir build/four_spin_design
```

The interaction space now has six pair coefficients, whose solution has one
free scale. The collective zero stays at $\pi/5$, while the response between
zeros changes and the collective factor becomes $\cos^3(5t/2)$.

Tests: [`test_inverse_spin.py`](../../tests/test_inverse_spin.py). The tests
include unequal and negative exchange, a disconnected numerical control,
inadmissible inputs and failure of a rerun after an earlier successful result.

## From this example to a discovery candidate

A useful next question is whether a physically accessible interaction or
control sequence produces a comparable response with less demanding coupling
requirements, or what quantitative error remains when those requirements
cannot be met. The constructor would solve for the permitted coefficients,
calculate the resulting signal and compare it with the collective reference.
Noise, pulse design and automatic selection of a new interaction family are
not implemented by this command.

The proposal becomes scientifically interesting through a result: an unknown
admissible family, an observable consequence beyond an existing construction,
a useful robustness bound, or a demonstrated obstruction within a specified
class. A smaller observable closure can arise from degeneracy, selection of
observables or coincident frequencies; it is a clue to examine, not itself
evidence of a new symmetry or a novelty decision.

Collective spin dynamics and inverse symmetry construction have substantial
precedents: [Kitagawa and Ueda](https://doi.org/10.1103/PhysRevA.47.5138) and
[Moudgalya and Motrunich](https://arxiv.org/abs/2302.03028). Comparison with
those results and the relevant target-system literature is a separate part
of a discovery study. A successful experiment would test the physical
realization; it would not substitute for that comparison.

The collective interaction $Q$ is quadratic in the total magnetization and is
constant within a sector of fixed magnetization. [Chapter 24](24_spin_language.md)
uses it as an example of a term that enlarges the algebra of a spin rotation,
and derives the exchange couplings, proportional to $\sqrt{(j+1)(N-1-j)}$, with
which a chain carries such a rotation.

[Tutorial](index.md) | [Retrieved source calculations](13_retrieval_to_calculation.md)
