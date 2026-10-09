"""Heredity: which written states pass to daughters through growth and division.

Module 1, inheritance through a threshold: a body whose order exists only above a critical size divides below it; the
daughter's order decays in the dip, and growth carries it back through a supercritical pitchfork, where it keeps its
parent's sign with probability P = Phi(phi_c / sigma_c).
"""
from .spec import SCHEMA, SpecError, load

__all__ = ["SCHEMA", "SpecError", "load"]
