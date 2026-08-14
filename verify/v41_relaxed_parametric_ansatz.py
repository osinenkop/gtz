#!/usr/bin/env python3
"""Parametric exact check for the relaxed-bound boundary ansatz.

The numerical sweep in ``v40_relaxed_target_sweep.py`` suggests that the same
structured boundary family controls the reduced relaxed obstruction for targets
slightly above ``1/8``.  Put

    l = 1 - 5t,   u = 5t = 1 - l.

The family is

    [  u  -x  -l  -y  -z  -z ]
    [ -x   a  -x   w   y   y ]
    [ -l  -x   u  -y  -z  -z ]
    [ -y   w  -y 1-a -x  -x ]
    [ -z   y  -z  -x   l   l ]
    [ -z   y  -z  -x   l   l ].

This script verifies over QQ(l) that idempotence reduces to four quadratic
equations and that adding the active determinant equations

    det(P_{013} - qI) = 0,
    det(P_{014} - qI) = 0,
    det(P_{014} - (1-q)I) = 0

forces the q-only relation

    (2q - 1) * (3q^2 - (3-l)q + l) = 0.

Equivalently, with l = 1 - 5t, the lower ansatz branch is

    q(t) = (2 + 5t - sqrt(25t^2 + 80t - 8)) / 6.

This is an exact check of the boundary family seen numerically; it is not a
global certificate for the relaxed theorem.
"""

from __future__ import annotations

import sympy as sp


def ansatz_matrix():
    l, a, x, y, z, w, q = sp.symbols("l a x y z w q")
    u = 1 - l
    p = sp.Matrix([
        [u, -x, -l, -y, -z, -z],
        [-x, a, -x, w, y, y],
        [-l, -x, u, -y, -z, -z],
        [-y, w, -y, 1 - a, -x, -x],
        [-z, y, -z, -x, l, l],
        [-z, y, -z, -x, l, l],
    ])
    return p, (l, a, x, y, z, w, q)


def idempotence_equations(p, variables):
    l, a, x, y, z, w, _ = variables
    equations = [
        2 * l**2 - l + x**2 + y**2 + 2 * z**2,
        -a * x + 2 * l * x - w * y - 2 * y * z,
        a * y + 2 * l * y - w * x + 2 * x * z - y,
        a**2 - a + w**2 + 2 * x**2 + 2 * y**2,
    ]

    groebner = sp.groebner(
        equations,
        a,
        x,
        y,
        z,
        w,
        order="grevlex",
        domain=sp.QQ.frac_field(l),
    )
    residual = p * p - p
    for i in range(6):
        for j in range(i, 6):
            _, reduced = groebner.reduce(sp.factor(residual[i, j]))
            assert sp.factor(reduced) == 0
    return equations


def det_eq(p, triple, value):
    block = p.extract(triple, triple)
    return sp.factor((block - value * sp.eye(3)).det())


def exact_relation(p, variables, idem):
    l, a, x, y, z, w, q = variables
    active = [
        det_eq(p, (0, 1, 3), q),
        det_eq(p, (0, 1, 4), q),
        det_eq(p, (0, 1, 4), 1 - q),
    ]
    groebner = sp.groebner(
        idem + active,
        a,
        x,
        y,
        z,
        w,
        q,
        order="grevlex",
        domain=sp.QQ.frac_field(l),
    )
    relation = sp.factor((2 * q - 1) * (l * q + l + 3 * q**2 - 3 * q))
    _, reduced = groebner.reduce(relation)
    assert sp.factor(reduced) == 0

    q_only = []
    for poly in groebner.polys:
        expr = sp.factor(poly.as_expr())
        if not any(var in expr.free_symbols for var in (a, x, y, z, w)):
            q_only.append(expr)
    assert any(sp.factor(6 * expr - relation) == 0 for expr in q_only)
    return relation, q_only


def lower_branch_formula():
    t = sp.symbols("t")
    q_lower = sp.factor((2 + 5 * t - sp.sqrt(25 * t**2 + 80 * t - 8)) / 6)
    return t, q_lower


def check_specializations(relation, variables):
    l, _, _, _, _, _, q = variables
    t, q_lower = lower_branch_formula()
    cases = [
        (sp.Rational(1, 8), (sp.Integer(7) - sp.sqrt(17)) / 16),
        (sp.Rational(16, 125), (sp.Integer(11) - sp.sqrt(46)) / 25),
        (sp.Rational(9, 70), sp.Rational(1, 6)),
    ]
    for target, expected in cases:
        l_value = 1 - 5 * target
        quad = sp.factor((relation / (2 * q - 1)).subs(l, l_value))
        assert sp.simplify(q_lower.subs(t, target) - expected) == 0
        assert sp.factor(quad.subs(q, expected)) == 0

    q16 = cases[1][1]
    assert q16 > sp.Rational(16, 125)
    assert cases[2][1] == sp.Rational(1, 6)


def main() -> int:
    p, variables = ansatz_matrix()
    l, _, _, _, _, _, q = variables
    idem = idempotence_equations(p, variables)
    relation, q_only = exact_relation(p, variables, idem)
    check_specializations(relation, variables)
    t, q_lower = lower_branch_formula()

    print("=" * 78)
    print("PARAMETRIC RELAXED ANSATZ CHECK")
    print("idempotence equations:")
    for eq in idem:
        print(f"  {sp.factor(eq)} = 0")
    print("\nGroebner q-only relation over QQ(l):")
    for expr in q_only:
        print(f"  {expr}")
    print("\nEquivalently:")
    print(f"  {sp.factor(relation)} = 0")
    print(f"  lower q(t) = {q_lower}")
    print("\nSpecializations:")
    for target in [sp.Rational(1, 8), sp.Rational(16, 125), sp.Rational(9, 70)]:
        q_value = sp.factor(q_lower.subs(t, target))
        slack = sp.factor(q_value - target)
        quad = sp.factor((l * q + l + 3 * q**2 - 3 * q).subs(l, 1 - 5 * target))
        print(f"  t={target}: q={q_value} = {sp.N(q_value, 18)}, slack={sp.N(slack, 18)}")
        print(f"       quadratic factor: {quad}")
    print("all exact parametric ansatz checks passed")
    print("=" * 78)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
