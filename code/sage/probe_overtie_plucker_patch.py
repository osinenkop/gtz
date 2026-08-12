#!/usr/bin/env python3
"""Patch-parametrized Plucker ansatz probes.

This is a more aggressive local version of ``probe_overtie_plucker_ansatz.py``.
It solves the observed determinant-one relation on a nonzero coordinate patch
that contains the refined numerical root, substitutes into the active/tie
system over a fraction field, and clears denominators.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from sage.all import GF, QQ, Matrix, PolynomialRing

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from gtz63_semialgebraic import TRIPLES  # noqa: E402


CASES = {
    "s7_78612_a_patch": {
        "active": (0, 1, 5, 8, 15, 16, 17),
        "ties": (2, 7, 10, 19),
        "variables": ("a", "b", "c", "x", "y", "e", "q"),
        "denominator": "a",
        "description": (
            "z0=a,z1=b,z6=c,z7=(b*c-1)/a,z3=x,z4=y,z5=e,"
            " z2=-a*y+b*x,z8=x*z7-y*c"
        ),
    },
    "s8_79656_c_patch": {
        "active": (0, 1, 2, 6, 14, 16, 18, 19),
        "ties": (10, 11, 15),
        "variables": ("a", "b", "c", "d", "x", "e", "q"),
        "denominator": "c",
        "description": (
            "z0=a,z1=b,z2=e,z3=c,z4=d,z6=x,z7=(1+d*x)/c,"
            " z5=-a*d+b*c,z8=-a*z7+b*x"
        ),
    },
}


def make_ring(case_name, characteristic, order):
    base = QQ if characteristic == 0 else GF(characteristic)
    return PolynomialRing(base, CASES[case_name]["variables"], order=order)


def substituted_z(case_name, ring):
    frac = ring.fraction_field()
    gens = [frac(g) for g in ring.gens()]
    if case_name == "s7_78612_a_patch":
        a, b, c, x, y, e, q = gens
        d = (b * c - 1) / a
        z = (
            a,
            b,
            -a * y + b * x,
            x,
            y,
            e,
            c,
            d,
            x * d - y * c,
        )
    elif case_name == "s8_79656_c_patch":
        a, b, c, d, x, e, q = gens
        y = (1 + d * x) / c
        z = (
            a,
            b,
            e,
            c,
            d,
            -a * d + b * c,
            x,
            y,
            -a * y + b * x,
        )
    else:  # pragma: no cover - argparse restricts choices.
        raise ValueError(case_name)
    return z, q


def chart_from_z(frac, z):
    y = Matrix(
        frac,
        [
            [1, 0, 0],
            [0, 1, 0],
            [0, 0, 1],
            [z[0], z[1], z[2]],
            [z[3], z[4], z[5]],
            [z[6], z[7], z[8]],
        ],
    )
    gram = y.transpose() * y
    d = gram.det()
    n = y * gram.adjugate() * y.transpose()
    return y, d, n


def block(n, d, frac, triple_index, theta):
    triple = TRIPLES[triple_index]
    return Matrix(
        frac,
        3,
        3,
        lambda i, j: 6 * n[triple[i], triple[j]]
        - (theta * d if i == j else frac.zero()),
    )


def cleared_saturated_det(mat, d, ring):
    expr = mat.det() / (d**2)
    numerator = expr.numerator()
    try:
        return ring(numerator)
    except TypeError:
        return ring(str(numerator))


def build_system(case_name, characteristic, order):
    case = CASES[case_name]
    ring = make_ring(case_name, characteristic, order)
    frac = ring.fraction_field()
    z, q = substituted_z(case_name, ring)
    _, d, n = chart_from_z(frac, z)
    polys = []
    labels = []
    for idx in case["active"]:
        polys.append(cleared_saturated_det(block(n, d, frac, idx, frac(1)), d, ring))
        labels.append(f"active:{idx}:{''.join(str(x) for x in TRIPLES[idx])}")
    for idx in case["ties"]:
        polys.append(cleared_saturated_det(block(n, d, frac, idx, q), d, ring))
        labels.append(f"tie:{idx}:{''.join(str(x) for x in TRIPLES[idx])}")
    return ring, z, d, polys, labels


def factor_summary(poly):
    try:
        return str(poly.factor())
    except Exception:
        return str(poly)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--case", choices=sorted(CASES), required=True)
    parser.add_argument("--characteristic", type=int, default=32003)
    parser.add_argument("--order", default="degrevlex")
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    ring, z, d, polys, labels = build_system(args.case, args.characteristic, args.order)
    ideal = ring.ideal(polys)
    print("=" * 78)
    print("OVERTIE PLUCKER PATCH PROBE")
    print(f"case: {args.case}")
    print(f"field: {'QQ' if args.characteristic == 0 else 'F_'+str(args.characteristic)}")
    print(f"variables: {ring.variable_names()}")
    print(f"order: {args.order}")
    print(f"ansatz: {CASES[args.case]['description']}")
    print(f"denominator patch: {CASES[args.case]['denominator']} != 0")
    print(f"equations: {len(polys)}")
    print(f"degrees: {[int(poly.total_degree()) for poly in polys]}")
    print("computing Groebner basis...", flush=True)
    gb = ideal.groebner_basis()
    gb_ideal = ring.ideal(gb)
    dim = gb_ideal.dimension()
    degree = None
    if dim == 0:
        try:
            degree = gb_ideal.vector_space_dimension()
        except Exception:
            degree = None
    print(f"basis size: {len(gb)}")
    print(f"dimension:  {dim}")
    print(f"degree:     {degree}")

    univariate = []
    variables = tuple(str(name) for name in ring.variable_names())
    for i, poly in enumerate(gb):
        used = [name for name in variables if ring(name) in poly.variables()]
        if len(used) == 1:
            univariate.append({"index": i, "variable": used[0], "polynomial": str(poly)})
            print(f"univariate {used[0]}: {poly}")

    q_polys = [row for row in univariate if row["variable"] == "q"]
    payload = {
        "case": args.case,
        "characteristic": args.characteristic,
        "order": args.order,
        "active": list(CASES[args.case]["active"]),
        "active_triples": [list(TRIPLES[idx]) for idx in CASES[args.case]["active"]],
        "ties": list(CASES[args.case]["ties"]),
        "tie_triples": [list(TRIPLES[idx]) for idx in CASES[args.case]["ties"]],
        "ansatz": CASES[args.case]["description"],
        "denominator": CASES[args.case]["denominator"],
        "z_substitution": [str(entry) for entry in z],
        "labels": labels,
        "variables": list(ring.variable_names()),
        "n_variables": len(ring.gens()),
        "n_equations": len(polys),
        "degrees": [int(poly.total_degree()) for poly in polys],
        "det_gram": str(d),
        "basis_size": len(gb),
        "dimension": int(dim),
        "degree": int(degree) if degree is not None else None,
        "univariate": univariate,
        "q_factorizations": [factor_summary(ring(row["polynomial"])) for row in q_polys],
        "basis": [str(poly) for poly in gb],
    }
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(payload, indent=1) + "\n")
    print(f"wrote {out}")
    print("=" * 78)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
