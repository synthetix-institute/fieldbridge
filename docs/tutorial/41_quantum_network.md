# Storing a photon behind a transmon: open networks with pulses in and out

**Objectives.** Write a network of resonators and a transmon as an open realization, with exchange couplings that are
switched on and off in time, the losses of every element, and a transmission line through which a single photon
arrives and leaves; compute what fraction of the photon comes back in the expected temporal mode; and account for
every lost part of it, loss by loss.

**Prerequisites.** [The quantum write](30_quantum_write.md) for open realizations on one mode and the Lindblad
equation; the input-output relation of a resonator coupled to a line (Gardiner and Collett, Phys. Rev. A 31, 3761
(1985)). Python 3.10 or newer with numpy and scipy.

**Command.** `fieldbridge quantum network [SPEC ...] --out-dir OUT` evolves every specification of
`examples/quantum/network` and writes `network.json` and `network.md` with the hashes of the inputs and the
implementation.

## 1. The question

A quantum memory in a superconducting circuit catches a photon that travels along a transmission line, holds it in a
long-lived resonator, and releases it later into the same line. The photon reaches the storage resonator only through
a nonlinear element, a transmon, which exchanges the excitation with the readout resonator on one side and with one
of several storage resonators on the other. Every element loses photons at its own rate, and the readout resonator
loses them into the line, which is how the photon arrives and leaves in the first place. What fraction of the photon
comes back in the temporal mode in which it is expected, and where does the rest go?

## 2. The mechanism

The realization ((Ω, Ξ); C, R, P; A) is written in the schema `fieldbridge-quantum-network/1`.

| Component | In the specification |
| --- | --- |
| Ξ, the carrier | `subsystems`: modes truncated at `max_quanta`, and multilevel systems with `levels` states; the product space of all of them |
| Ω, the operation | `hamiltonian`: terms with coefficients in the parameters, time and the `controls`; `dissipators`: Lindblad losses rate × D[operator] |
| C, the closures | the truncation of every mode and of the transmon; the rotating frame stated in the assumptions |
| R, the observable | the photons found in the output mode of each port, the photons emitted into each line in each window, the excitations each loss removes, and the `observables` recorded in time |
| P, the protocol | `controls` (piecewise in time), the `input` pulse of each port, `initial` and `duration` |
| A, the realization | `parameters` |

A mode named `r` gives the operators `r`, `rd` (creation) and `nr` (number). A multilevel system named `q` gives `q`,
`qd` and `nq` for its ladder, the operators of a transmon, and the projectors `q_j_k` = |j⟩⟨k| for any other
multilevel system.

**Lines.** A port couples an operator c of the network to a transmission line at a rate κ. The field that leaves the
network is b_out = b_in + √κ c. A port with nothing attached is the loss κ D[c]. A coherent input b_in = β(t) adds the
displacement of the line. A single photon in the temporal mode u(t) is emitted into the line by a virtual source mode
whose coupling releases exactly that mode, and the photon in a chosen output mode v(t) is absorbed by a virtual
detector mode (Kiilerich and Mølmer, Phys. Rev. Lett. 123, 123604 (2019)); the occupation of the detector at the end
is the number of photons the line carried in v. Source, network and detector form a cascade, in which each element
drives the next and none drives the one before it.

**The budget.** Photons that enter must leave or stay. The report adds the photons found in the output mode, those
emitted into the line outside it, those each loss removed and those left in the network, and compares the sum with
the photons supplied. Pure dephasing removes none. The sum equals the supply whenever the Hamiltonian conserves the
number of excitations, as the exchange couplings here do.

## 3. The example

`examples/quantum/network/transmon_storage.json` has a readout resonator r, a transmon q with three levels, two
storage resonators s1 and s2, and a line coupled to r at the rate κ, which sets the unit of time.

1. A single photon arrives with an exponentially rising amplitude, the temporal mode that an empty resonator absorbs
   completely, and is caught by r.
2. Two exchange pulses, each g t = π/2 with g = 50κ, move it from r to q and from q to s1.
3. After ten units of time, two more pulses move it from s1 through q to s2.
4. After ten more, two pulses bring it back through q to r, and r emits it as a decaying exponential, the measured
   output mode.

```bash
python3 -B -m fieldbridge quantum network --out-dir build/network
```

The run takes about ten seconds. Of one photon sent in, 0.930 comes back in the output mode. The rest is accounted
for to 10<sup>-15</sup>:

| Where the photon went | Photons |
| --- | --- |
| output mode | 0.930 |
| into the line while the swaps moved it out of and into r | 0.030 |
| internal loss of r | 0.020 |
| loss of s1 and of s2, ten units of time each | 0.010 and 0.010 |
| decay of the transmon | 0.0005 |
| left in the network | 0.0001 |

Each entry has a simple origin. While an exchange pulse moves the photon out of r, r keeps leaking into the line at
the rate κ; a pulse of length π/(2g) loses about κπ/(4g) on average, 0.016 per swap at g = 50κ, twice. The photon
spends about 1/κ in r while it is caught and 1/κ while it is released, so the internal loss of r takes about 2 κ_int
/ κ = 0.02. A storage resonator holding the photon for a time W keeps exp(−κ_s W). The transmon holds the photon
only during the swaps. Its third level is never reached with one photon, and the report records its population (zero)
as a check.

## 4. What the calculation checks

`tests/test_quantum_network.py` holds the limits in which the answer is known exactly.

- A photon sent from the source straight into a detector with the same mode is found there, to 10<sup>-6</sup>.
- An empty resonator catches a rising exponential truncated at a time T up to the part exp(−κT) that it reflects,
  and releases it into the decaying exponential.
- A coherent pulse drives the mean amplitude of an empty resonator as α' = −(κ/2) α − √κ β(t), and the photons it
  sends into the line are ∫|β + √κ α|<sup>2</sup> dt.
- The swept write of [the quantum write](30_quantum_write.md) is a network with one mode: `network.from_open`
  rewrites a `fieldbridge-quantum-open/1` specification with its sweep and bias as controls, and the density operator
  agrees with that of the one-mode solver to 10<sup>-9</sup> on the examples.
- When only a storage resonator holds the photon during a wait W, the photons in the output mode scale as
  exp(−κ_s W) exactly. A swap performed while a loss acts on an element leaves a small part of the photon behind;
  that part decays at the rate of the element where it was left, and it is why the scaling is exact only in this case.

## 5. Exercises

1. Lower g from 50κ to 10κ in the example and predict the photons lost into the line during the swaps before you
   run it.
2. Store the photon in s2 directly, without the move from s1, and compare the efficiency with that of the example.
3. Replace the single photon by a coherent pulse with 0.3 photons on average (`"kind": "coherent"`, with the same
   rising exponential as amplitude, the resonators truncated at two quanta and the output mode at
   `"max_photons": 2`; the product space then has 243 states, and a coherent state with 0.3 photons has three or
   more with probability 0.004). The transmon passes one excitation per swap: see how the component with two photons
   is lost or left behind, and read the population of the transmon's third level.

## 6. Limits

The solver evolves the full density operator, so the product of all dimensions, including one virtual mode for each
pulse in and out, is limited to 400. Rates of the losses and of the ports are constant; the couplings and drives
change in time through the controls. The rotating frame, the truncations and the values of the parameters are
assumptions of each specification, and the report says nothing about whether a device realizes them.
