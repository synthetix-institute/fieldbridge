"""Decision: when a population of units acts as one.

Three targets: the onset and growth of synchrony of limit-cycle units, reduced to phase oscillators (g(omega), H); the
swept collective write, P = Phi(pi^(1/4) h_s h / (D_s^(1/2) (a r)^(1/4))) with D_s proportional to 1/N, and the frozen
bias of a diverse sample; the passage from a seed after a step, whose 10-90% window is ln q_d in the growth exponent
(q_1 = 13.09, q_2 = 4.675).
"""
from .spec import SCHEMA, SpecError, load

__all__ = ["SCHEMA", "SpecError", "load"]
