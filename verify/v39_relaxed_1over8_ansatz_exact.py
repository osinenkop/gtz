#!/usr/bin/env python3
"""Exact symbolic check for the lower relaxed-1/8 core-boundary ansatz.

The SLSQP polish in ``v38_relaxed_1over8_slsqp_polish.py`` finds a structured
boundary point with

    F = q = (7 - sqrt(17))/16.

After row permutation its projector has the five-variable form

    [ 5/8  -x  -3/8  -y  -z  -z ]
    [ -x    a   -x    w   y   y ]
    [ -3/8 -x   5/8  -y  -z  -z ]
    [ -y    w   -y   1-a -x  -x ]
    [ -z    y   -z   -x  3/8 3/8]
    [ -z    y   -z   -x  3/8 3/8].

This script verifies, over exact rationals, that ``P^2=P`` reduces to four
quadratic equations and that adding active determinant equations for a block
with eigenvalues ``q`` and ``1-q`` forces
``(2*q - 1)*(8*q^2 - 7*q + 1) = 0`` in the ansatz ideal.  This does not prove
the global relaxed theorem; it only certifies the algebra of the numerical
boundary target.
"""

from __future__ import annotations

import sympy as sp


def ansatz_matrix():
    a, x, y, z, w, q = sp.symbols("a x y z w q")
    p = sp.Matrix([
        [sp.Rational(5, 8), -x, -sp.Rational(3, 8), -y, -z, -z],
        [-x, a, -x, w, y, y],
        [-sp.Rational(3, 8), -x, sp.Rational(5, 8), -y, -z, -z],
        [-y, w, -y, 1 - a, -x, -x],
        [-z, y, -z, -x, sp.Rational(3, 8), sp.Rational(3, 8)],
        [-z, y, -z, -x, sp.Rational(3, 8), sp.Rational(3, 8)],
    ])
    return p, (a, x, y, z, w, q)


def idempotence_equations(p, variables):
    a, x, y, z, w, _ = variables
    equations = [
        32 * x**2 + 32 * y**2 + 64 * z**2 - 3,
        4 * a * x + 4 * w * y - 3 * x + 8 * y * z,
        4 * a * y - 4 * w * x + 8 * x * z - y,
        a**2 - a + w**2 + 2 * x**2 + 2 * y**2,
    ]
    groebner = sp.groebner(equations, *variables[:-1], order="grevlex")
    residual = p * p - p
    for i in range(6):
        for j in range(i, 6):
            _, reduced = groebner.reduce(sp.factor(32 * residual[i, j]))
            assert sp.factor(reduced) == 0
    return equations


def det_eq(p, triple, q):
    block = p.extract(triple, triple)
    return sp.factor((block - q * sp.eye(3)).det())


def main() -> int:
    p, variables = ansatz_matrix()
    a, x, y, z, w, q = variables
    idem = idempotence_equations(p, variables)

    active_equations = [
        det_eq(p, (0, 1, 3), q),
        det_eq(p, (0, 1, 4), q),
        det_eq(p, (0, 1, 4), 1 - q),
    ]
    groebner = sp.groebner(idem + active_equations, a, x, y, z, w, q, order="grevlex")
    basis = [sp.factor(poly.as_expr()) for poly in groebner.polys]

    q_poly = sp.factor((2 * q - 1) * (8 * q**2 - 7 * q + 1))
    assert any(sp.factor(expr * 16) == q_poly for expr in basis) or q_poly in basis
    assert sp.factor(8 * det_eq(p, (0, 1, 3), q)).subs(q, sp.Rational(1, 8)) != 0

    q_low = (sp.Integer(7) - sp.sqrt(17)) / 16
    slack = sp.simplify(q_low - sp.Rational(1, 8))
    assert sp.simplify(8 * q_low**2 - 7 * q_low + 1) == 0
    assert sp.simplify(q_poly.subs(q, q_low)) == 0
    assert slack == (sp.Integer(5) - sp.sqrt(17)) / 16
    assert q_low > sp.Rational(1, 8)

    print("=" * 78)
    print("EXACT RELAXED-1/8 ANSATZ CHECK")
    print("idempotence equations:")
    for eq in idem:
        print(f"  {sp.factor(eq)} = 0")
    print("\nGroebner basis contains the q-only factor:")
    for expr in basis:
        if sp.factor(expr * 16) == q_poly or expr == q_poly:
            print(f"  {expr}")
    print(f"\nq = (7 - sqrt(17))/16 = {sp.N(q_low, 18)}")
    print(f"q - 1/8 = (5 - sqrt(17))/16 = {sp.N(slack, 18)}")
    print("all exact ansatz checks passed")
    print("=" * 78)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
