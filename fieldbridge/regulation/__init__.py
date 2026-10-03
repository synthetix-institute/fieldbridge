"""REGULATE: the return of an output to its set point after a step of an input.

A realization regulates when its steady output does not depend on a constant input. The module decides this from the
equations: the static gain of the output with and without a clamped variable, an integrator found by linear algebra
on the sampled drift (a function phi of the state with dphi/dt = g (y - y0)), and the response to steps. An
integrator together with a stable steady state makes the return exact at every parameter point where that steady
state exists (the internal model principle); a leaky integrator leaves a fraction of the step.
"""
