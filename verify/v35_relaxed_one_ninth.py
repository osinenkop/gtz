#!/usr/bin/env python3
"""Exact arithmetic checks for the determinant-sum relaxed GTZ(6,3) proof.

The proof itself is analytic.  This script records the constants used:

* t_* = (1 - sqrt(3/5))/2 > 1/9;
* t_* < 2/17, which is the trace-separation ceiling;
* shifted determinant sum:
      sum_T det(P_TT - t I) = 1 - 12t + 30t^2 - 20t^3,
  whose first positive root is t_*.
"""

from dataclasses import dataclass
from fractions import Fraction


@dataclass(frozen=True)
class Qr:
    """Element a + b*r with r^2 = 3/5."""

    a: Fraction
    b: Fraction = Fraction(0)

    def __add__(self, other):
        other = as_qr(other)
        return Qr(self.a + other.a, self.b + other.b)

    def __sub__(self, other):
        other = as_qr(other)
        return Qr(self.a - other.a, self.b - other.b)

    def __mul__(self, other):
        other = as_qr(other)
        return Qr(
            self.a * other.a + Fraction(3, 5) * self.b * other.b,
            self.a * other.b + self.b * other.a,
        )

    def __rmul__(self, other):
        return self * other

    def __neg__(self):
        return Qr(-self.a, -self.b)

    def __pow__(self, n: int):
        out = Qr(Fraction(1))
        for _ in range(n):
            out *= self
        return out

    def __str__(self):
        return f"{self.a} + ({self.b})*sqrt(3/5)"


def as_qr(value):
    if isinstance(value, Qr):
        return value
    return Qr(Fraction(value))


def shifted_sum(t):
    t = as_qr(t)
    return Qr(Fraction(1)) - 12 * t + 30 * (t ** 2) - 20 * (t ** 3)


def main() -> int:
    t_one_ninth = Fraction(1, 9)
    t_one_eighth = Fraction(1, 8)
    t_star = Qr(Fraction(1, 2), Fraction(-1, 2))

    print("RELAXED GTZ(6,3) DETERMINANT-SUM CONSTANT CHECK")
    print(f"t_* = {t_star}")
    print(f"p(t_*) = {shifted_sum(t_star)}")
    print(f"p(1/9) = {shifted_sum(t_one_ninth)}")
    print(f"p(1/8) = {shifted_sum(t_one_eighth)}")

    assert shifted_sum(t_star) == Qr(Fraction(0))
    assert shifted_sum(t_one_ninth) == Qr(Fraction(7, 729))
    assert shifted_sum(t_one_eighth) == Qr(Fraction(-9, 128))

    # t_* > 1/9 iff sqrt(3/5) < 7/9, checked by squaring positive sides.
    assert Fraction(3, 5) < Fraction(7, 9) ** 2
    # t_* < 2/17 iff sqrt(3/5) > 13/17.
    assert Fraction(3, 5) > Fraction(13, 17) ** 2
    print("all exact checks passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
