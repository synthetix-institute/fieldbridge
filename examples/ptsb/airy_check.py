"""Check of the conservative local law (`ptsb.airy_law`) against the linear transverse equation.

y'' = R t y - f, entered at t0 = -T0 R^{-1/3} in thermal equilibrium of the local well (mean -f/omega0^2, variances
kT/omega0^2 and kT, action J = kT/omega0), integrated to t = 6 R^{-1/3}. The law P = Phi(-pi^{1/2} f/(J R)^{1/2})
assumes an adiabatic entry; T0 = 6 is close to it and T0 = 20 nearly exact.

    python airy_check.py
"""
import sys
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

sys.path.insert(0, str(Path(__file__).resolve().parent))
from ptsb import Phi, airy_law  # noqa: E402

worst = {6.0: 0.0, 20.0: 0.0}
for R in (0.1, 1.0, 10.0):
    for f, kT in ((0.01, 0.01), (0.03, 0.01), (0.01, 0.1)):
        for T0 in (6.0, 20.0):
            t0 = -T0 / R ** (1 / 3)
            om0 = np.sqrt(-R * t0)

            def rhs(t, z):
                y, p, m00, m01, m10, m11 = z
                w2 = -R * t
                return [p, -w2 * y - f, m10, m11, -w2 * m00, -w2 * m01]

            s = solve_ivp(rhs, (t0, 6 / R ** (1 / 3)), [-f / om0 ** 2, 0, 1, 0, 0, 1], rtol=1e-10, atol=1e-13,
                          method="DOP853")
            y, m00, m01 = s.y[0, -1], s.y[2, -1], s.y[3, -1]
            P_lin = float(Phi(y / np.sqrt(m00 ** 2 * kT / om0 ** 2 + m01 ** 2 * kT)))
            P_airy = airy_law(f, kT / om0, R)
            worst[T0] = max(worst[T0], abs(P_lin - P_airy))
            print(f"R {R:5.1f}  f {f}  kT {kT}  T0 {T0:4.0f}:  linear equation {P_lin:.4f}   local law {P_airy:.4f}")
print("largest difference:", {k: round(v, 4) for k, v in worst.items()})
