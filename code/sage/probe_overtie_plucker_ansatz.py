#!/usr/bin/env python3
"""Reduced exact probes for Plucker-pattern over-tied roots.

The numerical relation screen found that two hard low-active roots satisfy
height-one Plucker equalities.  This script substitutes the corresponding
lower-dimensional chart directly before building the active/tie equations.
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
    "s7_78612": {
        "active": (0, 1, 5, 8, 15, 16, 17),
        "ties": (2, 7, 10, 19),
        "variables": ("a", "b", "c", "d", "x", "y", "e", "q"),
        "description": (
            "z0=a,z1=b,z6=c,z7=d,z3=x,z4=y,z5=e,"
            " z2=-a*y+b*x,z8=x*d-y*c, a*d-b*c=-1"
        ),
    },
    "s8_79656": {
        "active": (0, 1, 2, 6, 14, 16, 18, 19),
        "ties": (10, 11, 15),
        "variables": ("a", "b", "c", "d", "x", "y", "e", "q"),
        "description": (
            "z0=a,z1=b,z2=e,z3=c,z4=d,z6=x,z7=y,"
            " z5=-a*d+b*c,z8=-a*y+b*x, c*y-d*x=1"
        ),
    },
}


def make_ring(case_name, characteristic, order):
    base = QQ if characteristic == 0 else GF(characteristic)
    return PolynomialRing(base, CASES[case_name]["variables"], order=order)


def substituted_z(case_name, ring):
    a, b, c, d, x, y, e, q = ring.gens()
    if case_name == "s7_78612":
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
        relation = a * d - b * c + 1
    elif case_name == "s8_79656":
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
        relation = c * y - d * x - 1
    else:  # pragma: no cover - argparse restricts choices.
        raise ValueError(case_name)
    return z, q, relation


def chart_from_z(ring, z):
    y = Matrix(
        ring,
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


def saturated_det(mat, d, ring):
    raw = mat.det()
    quotient, remainder = raw.quo_rem(d**2)
    if remainder != ring.zero():
        raise ArithmeticError("determinant is not divisible by d^2")
    return quotient


def block(n, d, ring, triple_index, theta):
    triple = TRIPLES[triple_index]
    return Matrix(
        ring,
        3,
        3,
        lambda i, j: 6 * n[triple[i], triple[j]]
        - (theta * d if i == j else ring.zero()),
    )


def build_system(case_name, characteristic, order):
    case = CASES[case_name]
    ring = make_ring(case_name, characteristic, order)
    z, q, relation = substituted_z(case_name, ring)
    _, d, n = chart_from_z(ring, z)
    polys = []
    labels = []
    for idx in case["active"]:
        polys.append(saturated_det(block(n, d, ring, idx, ring.one()), d, ring))
        labels.append(f"active:{idx}:{''.join(str(x) for x in TRIPLES[idx])}")
    for idx in case["ties"]:
        polys.append(saturated_det(block(n, d, ring, idx, q), d, ring))
        labels.append(f"tie:{idx}:{''.join(str(x) for x in TRIPLES[idx])}")
    polys.append(relation)
    labels.append(f"ansatz:{relation}")
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
    print("OVERTIE PLUCKER ANSATZ PROBE")
    print(f"case: {args.case}")
    print(f"field: {'QQ' if args.characteristic == 0 else 'F_'+str(args.characteristic)}")
    print(f"variables: {ring.variable_names()}")
    print(f"order: {args.order}")
    print(f"ansatz: {CASES[args.case]['description']}")
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
