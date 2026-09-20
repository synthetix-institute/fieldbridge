# Which observables determine a quantum response?

Two preparations can have the same magnetization of one spin and different
correlations with another spin. If the interaction converts that correlation
into magnetization, the first measurement alone does not determine the later
signal. The missing information can be found from the Hamiltonian rather than
chosen from a catalogue of additional variables.

Take two spin-1/2 degrees of freedom. Let X, Y and Z denote Pauli matrices and
let I denote the identity. In units with hbar=1, choose

```math
H=g\,Z\otimes Z,\qquad O_0=X\otimes I,\qquad g>0.
```

The Heisenberg equation gives the observable derivative i[H,O]. Repeatedly
applying this operation to O0 produces

```math
O_1=i[H,O_0]=-2g\,Y\otimes Z,
\qquad i[H,O_1]=-4g^2 O_0.
```

The second derivative introduces no new independent observable. Thus two
expectation values, m_i=Tr(rho O_i), obey a closed equation for every initial
density matrix rho:

```math
\frac{d}{dt}\begin{pmatrix}m_0\\m_1\end{pmatrix}
=\begin{pmatrix}0&1\\-4g^2&0\end{pmatrix}
\begin{pmatrix}m_0\\m_1\end{pmatrix}.
```

Writing x(t)=E[X tensor I] and c(t)=E[Y tensor Z] gives
x(t)=x(0) cos(2gt)-c(0) sin(2gt). The correlation supplies an independent
initial condition. In the product states with the first spin along +y or -y
and the second along +z, x(0)=0 in both cases, but c(0)=+1 or -1. The two
signals therefore depart in opposite directions. Entanglement is not required
for this failure of a one-observable description.

## How the program finds the additional observable

[`quantum_closure`](../../fieldbridge/verification.py) starts with the supplied
observable, forms its commutator with H and appends the result only if it is
linearly independent of the existing span. Every appended observable is
processed in turn. The algorithm terminates because the operator space is
finite-dimensional. It then solves for the coefficients of the commutator
within the final basis and checks each matrix identity exactly.

The resulting span is the smallest invariant linear subspace containing the
specified observable for this generator, at generic declared parameter values.
It is not a reconstruction of the full density matrix and does not make a
claim about the smallest nonlinear realization of all possible experiments.
Special parameter values can reduce its dimension; setting g to zero in the
input restores one-observable closure.

```bash
python3 -B -m fieldbridge verify-construction \
  examples/construction/quantum_correlations.json \
  --out-dir build/quantum_correlations
```

The input is the full four-dimensional Hamiltonian matrix and the measured
operator, not the two-dimensional answer. The tests compare the derived
evolution with full unitary evolution and check the uncoupled control. A
non-Hermitian Hamiltonian is rejected by this closed-system handler.

This calculation complements the quantum book: its state-factorization,
generator and observable chapters become parts of one predictive problem.
The correlation matters because the interaction carries it into an observed
quantity. This is a known quantum example used to test the construction
procedure, not a claim of a new quantum effect.

For the underlying Heisenberg equation, see
[Cambridge's quantum mechanics notes](https://www.damtp.cam.ac.uk/user/tong/qm/qmhtml/S3.html).
Continue with [reproduction and discovery](12_reproduction_and_discovery.md).

## Exercise: let a field compete with the interaction

Add a transverse field on the first spin, so that
$H=g Z\otimes Z+h X\otimes I$, with both parameters positive. The field
commutes with the initially measured operator, so it does not change its first
derivative. It does change the evolution of the correlation: $Y\otimes Z$
now couples to $Z\otimes Z$. The required span therefore grows from two to
three independent operators for generic nonzero parameters.

Create and calculate that variant:

```bash
python3 - <<'PY'
import json
from pathlib import Path

spec = json.loads(Path('examples/construction/quantum_correlations.json').read_text(encoding='utf-8'))
spec['parameters']['h'] = 'positive'
spec['hamiltonian'] = [['g',0,'h',0], [0,'-g',0,'h'],
                       ['h',0,'-g',0], [0,'h',0,'g']]
spec['question'] = 'Which expectations close when a transverse field competes with the spin interaction?'
spec['assumptions'][0] = 'Two spins with H=g sigma_z tensor sigma_z + h sigma_x tensor identity.'
spec['assumptions'].append('Both g and h are positive.')
spec['provenance']['origin'] = 'tutorial perturbation of an authored benchmark'
spec['provenance']['equation_ids'] = ['tutorial_transverse_field:H', 'tutorial_transverse_field:observable']
out = Path('build/tutorial_transverse_field.json')
out.parent.mkdir(parents=True, exist_ok=True)
out.write_text(json.dumps(spec, indent=2), encoding='utf-8')
PY
python3 -B -m fieldbridge verify-construction \
  build/tutorial_transverse_field.json --out-dir build/tutorial_transverse_field
```

Expect `observable_dimension: 3` and three true `closed_identities`.
The program's basis can contain scaled linear combinations rather than the
three Pauli products used in this explanation; compare their spans, not their
printed names. The physical reason for the enlarged description is the new
coupling in the observable equations.

For a second control, set the supplied Hamiltonian matrix to zero. The measured
observable is then constant and the span has dimension one. Changing the
generator, rather than manually changing the expected answer, is what makes
these useful tests.

[Tutorial](index.md)
