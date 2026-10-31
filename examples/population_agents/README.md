# Populations of agents against their mean-field equations

`compare.py` simulates the four populations of `examples/memory/populations` as individual events and compares them with
the equations and, for the two birth-death models, with the exact mean times of the chain. `figure.py` draws
`results.json`. [Module 17](../../docs/tutorial/42_memory_populations.md), section 7, explains the results.

The full run takes about four hours on four cores, nearly all of it the budworm sweeps at N = 100000. Each part is saved as
it finishes, and `--parts` reruns one and keeps the others in the output file. `--quick` leaves out the largest populations.
Numba makes the event loops about a hundred times faster, and the code runs without it.
