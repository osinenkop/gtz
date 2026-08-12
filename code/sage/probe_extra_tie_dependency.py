#!/usr/bin/env python3
"""Exact low-degree dependency checks for over-tied margin roots.

Given a square active-plus-tie subsystem and an extra tied outside determinant,
this script tests the cheapest possible exact explanation:

    extra determinant in the QQ-linear span of the selected determinant equations?

Failure here does not rule out local algebraic dependence.  It only says that the
extra equality is not explained by a global linear identity among the determinant
polynomials.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from types import SimpleNamespace

from sage.all import Matrix, QQ

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from gtz63_semialgebraic import TRIPLES, parse_indices  # noqa: E402
from probe_margin_tie_system import resolve_active, shifted_determinant  # noqa: E402
from refine_margin_tie_root import build_system  # noqa: E402


def parse_case(text: str):
    parts = text.split(":")
    if len(parts) != 4:
        raise ValueError("case must have form size:canon:tie0,tie1:extra")
    size = int(parts[0])
    canon = int(parts[1])
    ties = parse_indices(parts[2])
    extra = int(parts[3])
    return size, canon, ties, extra


def tlabel(idx: int) -> str:
    return "".join(str(x) for x in TRIPLES[idx])


def check_case(size: int, canon: int, ties: tuple[int, ...], extra: int, order: str):
    active = resolve_active(SimpleNamespace(size=size, canon=canon, active_indices=None))
    ring, chart, _variables, polynomials, _jacobian = build_system(active, ties, order)
    q = dict(zip(ring.variable_names(), ring.gens()))["q"]
    extra_poly = shifted_determinant(chart, extra, q)
    all_polys = list(polynomials) + [extra_poly]
    monomials = sorted(set().union(*(poly.dict().keys() for poly in all_polys)))
    rows = []
    for poly in all_polys:
        coeffs = poly.dict()
        rows.append([coeffs.get(monomial, QQ(0)) for monomial in monomials])
    matrix = Matrix(QQ, rows)
    base_rank = matrix[:-1, :].rank()
    extended_rank = matrix.rank()
    linear_member = bool(base_rank == extended_rank)
    coeffs = None
    if linear_member:
        solution = Matrix(QQ, rows[:-1]).transpose().solve_right(
            Matrix(QQ, rows[-1]).transpose()
        )
        coeffs = [str(c) for c in solution.list()]
    return {
        "size": size,
        "canon": canon,
        "active": list(active),
        "active_triples": [list(TRIPLES[i]) for i in active],
        "ties": list(ties),
        "tie_labels": [tlabel(i) for i in ties],
        "extra": extra,
        "extra_label": tlabel(extra),
        "order": order,
        "monomial_count": len(monomials),
        "base_rank": int(base_rank),
        "extended_rank": int(extended_rank),
        "linear_member": linear_member,
        "coefficients": coeffs,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--case",
        action="append",
        default=[],
        help="case as size:canon:tie0,tie1,...:extra; may be repeated",
    )
    parser.add_argument("--size", type=int)
    parser.add_argument("--canon", type=int)
    parser.add_argument("--tie-indices", default="")
    parser.add_argument("--extra-index", type=int)
    parser.add_argument("--order", default="degrevlex")
    parser.add_argument("--out", default="")
    args = parser.parse_args()

    cases = [parse_case(item) for item in args.case]
    if args.size is not None or args.canon is not None or args.tie_indices or args.extra_index is not None:
        if args.size is None or args.canon is None or not args.tie_indices or args.extra_index is None:
            raise SystemExit("single-case mode needs --size, --canon, --tie-indices, --extra-index")
        cases.append((args.size, args.canon, parse_indices(args.tie_indices), args.extra_index))
    if not cases:
        raise SystemExit("provide at least one --case or a complete single case")

    results = [check_case(*case, order=args.order) for case in cases]
    payload = {"results": results}
    if args.out:
        Path(args.out).parent.mkdir(parents=True, exist_ok=True)
        Path(args.out).write_text(json.dumps(payload, indent=1) + "\n")
    print("=" * 78)
    print("EXTRA TIE DEPENDENCY PROBE")
    for row in results:
        print(
            f"({row['size']},{row['canon']}) ties={row['ties']} extra={row['extra']} "
            f"rank {row['base_rank']}->{row['extended_rank']} "
            f"linear_member={row['linear_member']}"
        )
    if args.out:
        print(f"wrote {args.out}")
    print("=" * 78)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
