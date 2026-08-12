#!/usr/bin/env python3
"""Exact probe for the structured size-6 low-active margin ansatz.

The refined over-tied size-6 root in
``refine_overtie_s6_78593_t5_7_12_14_18_p600.json`` has the chart pattern

    Z = [ -a   a  -b
           c   d  -a
          -d  -c  -a ].

This script substitutes that four-parameter ansatz into the active determinant
equations and the five outside tie determinant equations.  The goal is to test
whether the hard singular over-tied local root is explained by a small exact
subsystem rather than by the full ten-variable chart.
"""
from __future__ import annotations

import argparse
import itertools
import json
import sys
from pathlib import Path

from sage.all import GF, QQ, Matrix, PolynomialRing

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from gtz63_semialgebraic import TRIPLES  # noqa: E402


ACTIVE = (0, 1, 9, 15, 16, 17)
TIES = (5, 7, 12, 14, 18)


def block(chart_n, d, ring, triple_index: int, theta):
    triple = TRIPLES[triple_index]
    return Matrix(
        ring,
        3,
        3,
        lambda i, j: 6 * chart_n[triple[i], triple[j]]
        - (theta * d if i == j else ring.zero()),
    )


def saturated_det(mat, d, ring):
    raw = mat.det()
    quotient, remainder = raw.quo_rem(d**2)
    if remainder != ring.zero():
        raise ArithmeticError("determinant is not divisible by d^2")
    return quotient


def build_system(characteristic: int, include_q_stationarity: bool, branch: str, order: str):
    base = QQ if characteristic == 0 else GF(characteristic)
    ring = PolynomialRing(base, ("a", "b", "c", "d", "q"), order=order)
    a, b, c, e, q = ring.gens()

    y = Matrix(
        ring,
        [
            [1, 0, 0],
            [0, 1, 0],
            [0, 0, 1],
            [-a, a, -b],
            [c, e, -a],
            [-e, -c, -a],
        ],
    )
    gram = y.transpose() * y
    det_gram = gram.det()
    n = y * gram.adjugate() * y.transpose()

    polys = []
    labels = []
    for idx in ACTIVE:
        polys.append(saturated_det(block(n, det_gram, ring, idx, ring.one()), det_gram, ring))
        labels.append(f"active:{idx}:{''.join(str(x) for x in TRIPLES[idx])}")
    for idx in TIES:
        polys.append(saturated_det(block(n, det_gram, ring, idx, q), det_gram, ring))
        labels.append(f"tie:{idx}:{''.join(str(x) for x in TRIPLES[idx])}")

    if branch == "cplusd_sq5":
        polys.append((c + e) ** 2 - 5)
        labels.append("branch:(c+d)^2-5")
    elif branch == "cdiff_sq5":
        polys.append((c - e) ** 2 - 5)
        labels.append("branch:(c-d)^2-5")
    elif branch == "a_zero":
        polys.append(a)
        labels.append("branch:a")
    elif branch != "none":
        raise ValueError(f"unknown branch {branch!r}")

    if include_q_stationarity:
        # On the ansatz curve, q is minimized when the active/tie equations have
        # a nonzero tangent annihilating dq.  A cheap necessary algebraic proxy is
        # that the Jacobian of all equations with respect to (a,b,c,d) drops rank.
        jac = Matrix(ring, len(polys), 4, [[poly.derivative(var) for var in (a, b, c, e)] for poly in polys])
        minors = list(jac.minors(4))
        polys.extend(minors)
        labels.extend(f"rank4minor:{i}" for i in range(len(minors)))

    return ring, det_gram, polys, labels


def factor_summary(poly):
    try:
        return str(poly.factor())
    except Exception:
        return str(poly)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--characteristic", type=int, default=0)
    parser.add_argument("--include-q-stationarity", action="store_true")
    parser.add_argument(
        "--branch",
        choices=("none", "cplusd_sq5", "cdiff_sq5", "a_zero"),
        default="none",
        help="Optional exact branch relation suggested by the grevlex tail.",
    )
    parser.add_argument("--order", default="degrevlex")
    parser.add_argument("--timeout-note", default="")
    parser.add_argument("--out", default="code/sage/out/s6_ansatz_QQ.json")
    args = parser.parse_args()

    ring, det_gram, polys, labels = build_system(
        args.characteristic,
        args.include_q_stationarity,
        args.branch,
        args.order,
    )
    ideal = ring.ideal(polys)
    print("=" * 78)
    print("SIZE-6 STRUCTURED ANSATZ")
    print(f"field: {'QQ' if args.characteristic == 0 else 'F_'+str(args.characteristic)}")
    print(f"variables: {ring.variable_names()}")
    print(f"order: {args.order}")
    print(f"branch: {args.branch}")
    print(f"equations: {len(polys)}")
    print(f"degrees: {[int(poly.total_degree()) for poly in polys[:len(ACTIVE)+len(TIES)]]}")
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
        "characteristic": args.characteristic,
        "active": list(ACTIVE),
        "active_triples": [list(TRIPLES[idx]) for idx in ACTIVE],
        "ties": list(TIES),
        "tie_triples": [list(TRIPLES[idx]) for idx in TIES],
        "ansatz": "Z=[[-a,a,-b],[c,d,-a],[-d,-c,-a]]",
        "include_q_stationarity": args.include_q_stationarity,
        "branch": args.branch,
        "order": args.order,
        "n_variables": len(ring.gens()),
        "n_equations": len(polys),
        "degrees": [int(poly.total_degree()) for poly in polys],
        "labels": labels,
        "det_gram": str(det_gram),
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
