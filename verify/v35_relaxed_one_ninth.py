#!/usr/bin/env python3
"""Exact arithmetic checks for the relaxed F(P) >= 1/9 proof.

The proof itself is analytic.  This script records the rational constants used:

* low leverage cutoff: (1 - ell) / 5 >= 1/9 iff ell <= 4/9;
* high-leverage trace separation: 1 + 2/9 < 4/3;
* shifted determinant sum:
      sum_T det(P_TT - t I) = 1 - 12t + 30t^2 - 20t^3,
  which is positive at t = 1/9.
"""

from fractions import Fraction


def main() -> int:
    t = Fraction(1, 9)
    leverage_cutoff = 1 - 5 * t
    two_small_trace_cap = 1 + 2 * t
    core_trace_floor = Fraction(4, 3)
    shifted_sum = 1 - 12 * t + 30 * t * t - 20 * t * t * t

    print("RELAXED GTZ(6,3) CONSTANT CHECK")
    print(f"t = {t}")
    print(f"low leverage cutoff 1 - 5t = {leverage_cutoff}")
    print(f"two-small trace cap 1 + 2t = {two_small_trace_cap}")
    print(f"core trace floor = {core_trace_floor}")
    print(f"shifted determinant sum = {shifted_sum}")

    assert leverage_cutoff == Fraction(4, 9)
    assert two_small_trace_cap < core_trace_floor
    assert shifted_sum == Fraction(7, 729)
    assert shifted_sum > 0
    print("all exact checks passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
