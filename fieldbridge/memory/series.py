"""Truncated power series in one variable, the arithmetic of the order-by-order reduction at a write point.

A Series holds the coefficients c_0 .. c_N of c_0 + c_1 s + ... + c_N s^N and supports what a drift is written in: the
operators + - * / **, abs, the comparisons, % and the functions of spec.FUNCTIONS. The compiled expressions of
spec.py and the vectorized drifts of the library run on it unchanged: an object array holding Series is passed where
the drift takes a float array, and NumPy calls the method of the same name on each element (np.sin calls .sin(),
np.absolute calls __abs__). Nothing is parsed or evaluated here.

Compositions use the standard recurrences for the coefficients of a power series (Knuth, The Art of Computer
Programming, vol. 2, section 4.7; J.C.P. Miller's recurrence for a power). A function that is not analytic at the
constant term (abs at a kink, sqrt, log or a non-integer power at zero, division by zero) raises NotAnalytic. There is
no __float__, so that a drift that casts its state to float raises TypeError rather than losing the higher orders.
"""
from __future__ import annotations

import numpy as np


class NotAnalytic(ValueError):
    """The drift has no power series at this point; the message names the function."""


KINK = 1e-12  # abs is at a kink where its argument at s = 0 is below this share of its coefficients


class Series:
    __slots__ = ("c",)

    def __init__(self, c):
        self.c = np.asarray(c, float)

    def _lift(self, other):
        if isinstance(other, Series):
            return other.c
        out = np.zeros_like(self.c)
        out[0] = float(other)
        return out

    def __add__(self, other):
        return Series(self.c + self._lift(other))

    __radd__ = __add__

    def __neg__(self):
        return Series(-self.c)

    def __pos__(self):
        return self

    def __sub__(self, other):
        return Series(self.c - self._lift(other))

    def __rsub__(self, other):
        return Series(self._lift(other) - self.c)

    def __mul__(self, other):
        if isinstance(other, Series):
            return Series(np.convolve(self.c, other.c)[:self.c.size])
        return Series(self.c * float(other))

    __rmul__ = __mul__

    def inverse(self):
        a = self.c
        if a[0] == 0.0:
            raise NotAnalytic("division by a quantity that vanishes at the write point")
        b = np.zeros_like(a)
        b[0] = 1.0 / a[0]
        for k in range(1, a.size):
            b[k] = -np.dot(a[1:k + 1], b[k - 1::-1]) / a[0]
        return Series(b)

    def __truediv__(self, other):
        if isinstance(other, Series):
            return self * other.inverse()
        return Series(self.c / float(other))

    def __rtruediv__(self, other):
        return self.inverse() * other

    def __pow__(self, n):
        if isinstance(n, Series):
            return (n * self.log()).exp()
        n = float(n)
        if n == int(n) and n >= 0:
            out = Series(self._lift(1.0))
            base, m = self, int(n)
            while m:
                if m & 1:
                    out = out * base
                base = base * base
                m >>= 1
            return out
        a = self.c
        if a[0] <= 0.0 and n != int(n) or a[0] == 0.0:
            raise NotAnalytic(f"the power {n:g} of a quantity {a[0]:g} at the write point")
        b = np.zeros_like(a)
        b[0] = a[0] ** n
        for k in range(1, a.size):
            b[k] = sum((n * j - (k - j)) * a[j] * b[k - j] for j in range(1, k + 1)) / (k * a[0])
        return Series(b)

    def __rpow__(self, base):
        return (self * float(np.log(float(base)))).exp()

    def __abs__(self):
        if abs(self.c[0]) <= KINK * float(np.max(np.abs(self.c[1:]), initial=0.0)):
            raise NotAnalytic("abs at a kink: its argument vanishes at the write point")
        return Series(np.sign(self.c[0]) * self.c)

    def __mod__(self, period):
        # the wrapping of an angle: the same function of s, with the constant term in the fundamental period
        out = self.c.copy()
        out[0] = out[0] % float(period)
        return Series(out)

    # the comparisons read the value at s = 0, which selects the branch of a maximum, a minimum or a clip
    def _value(self, other):
        return other.c[0] if isinstance(other, Series) else float(other)

    def __lt__(self, other):
        return self.c[0] < self._value(other)

    def __le__(self, other):
        return self.c[0] <= self._value(other)

    def __gt__(self, other):
        return self.c[0] > self._value(other)

    def __ge__(self, other):
        return self.c[0] >= self._value(other)

    def conjugate(self):
        return self

    def exp(self):
        a = self.c
        e = np.zeros_like(a)
        e[0] = np.exp(a[0])
        for k in range(1, a.size):
            e[k] = sum(j * a[j] * e[k - j] for j in range(1, k + 1)) / k
        return Series(e)

    def log(self):
        a = self.c
        if a[0] <= 0.0:
            raise NotAnalytic(f"log of a quantity {a[0]:g} at the write point")
        b = np.zeros_like(a)
        b[0] = np.log(a[0])
        for k in range(1, a.size):
            b[k] = (a[k] - sum(j * b[j] * a[k - j] for j in range(1, k)) / k) / a[0]
        return Series(b)

    def sqrt(self):
        if self.c[0] <= 0.0:
            raise NotAnalytic(f"sqrt of a quantity {self.c[0]:g} at the write point")
        return self ** 0.5

    def _sin_cos(self):
        a = self.c
        s, c = np.zeros_like(a), np.zeros_like(a)
        s[0], c[0] = np.sin(a[0]), np.cos(a[0])
        for k in range(1, a.size):
            s[k] = sum(j * a[j] * c[k - j] for j in range(1, k + 1)) / k
            c[k] = -sum(j * a[j] * s[k - j] for j in range(1, k + 1)) / k
        return Series(s), Series(c)

    def sin(self):
        return self._sin_cos()[0]

    def cos(self):
        return self._sin_cos()[1]

    def tanh(self):
        # T' = (1 - T^2) a'
        a = self.c
        t, t2 = np.zeros_like(a), np.zeros_like(a)
        t[0] = np.tanh(a[0])
        t2[0] = t[0] * t[0]
        for k in range(1, a.size):
            t[k] = a[k] - sum(j * a[j] * t2[k - j] for j in range(1, k + 1)) / k
            t2[k] = np.dot(t[:k + 1], t[k::-1])
        return Series(t)


def coefficients(value, size: int) -> np.ndarray:
    """The coefficients of an entry of the result of a drift, which is a constant where the drift does not depend on
    the state."""
    if isinstance(value, Series):
        return value.c
    out = np.zeros(size)
    out[0] = float(value)
    return out
