#!/usr/bin/env python3
"""Screen square row bases inside an over-tied determinant system.

For an overdetermined active-plus-tie root with m > 10 equations, enumerate
10-row subsystems and rank them by the smallest singular value of the Jacobian
at the refined root.  This identifies good square bases for later interval
certification/local-membership work.
"""
from __future__ import annotations

import argparse
import itertools
import json
import sys
from pathlib import Path

from sage.all import Matrix, RDF, RealField, vector

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from gtz63_semialgebraic import TRIPLES  # noqa: E402
from refine_margin_tie_root import build_system, eval_matrix, eval_vector, max_abs  # noqa: E402


def label(prefix: str, idx: int) -> str:
    return f"{prefix}:{idx}:{''.join(str(a) for a in TRIPLES[idx])}"


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--precision", type=int, default=500)
    parser.add_argument("--order", default="degrevlex")
    parser.add_argument("--top", type=int, default=20)
    parser.add_argument("--require-all-ties", action="store_true")
    parser.add_argument("--out", default="")
    args = parser.parse_args()

    data = json.loads(Path(args.input).read_text())
    active = tuple(data["active"])
    ties = tuple(data["ties"])
    ring, chart, variables, polynomials, jacobian = build_system(active, ties, args.order)
    RF = RealField(args.precision)
    x = vector(RF, [RF(s) for s in list(data["z"]) + [data["q"]]])
    J = eval_matrix(jacobian, x, RF)
    F = eval_vector(polynomials, x, RF)
    labels = [label("active", i) for i in active] + [label("tie", i) for i in ties]
    tie_rows = set(range(len(active), len(active) + len(ties)))

    rows = []
    for subset in itertools.combinations(range(len(polynomials)), 10):
        subset_set = set(subset)
        if args.require_all_ties and not tie_rows <= subset_set:
            continue
        Js = Matrix(RDF, 10, 10, [RDF(J[i, j]) for i in subset for j in range(10)])
        sv = [float(v) for v in Js.singular_values()]
        min_sv = min(sv)
        max_sv = max(sv)
        omitted = [i for i in range(len(polynomials)) if i not in subset_set]
        rows.append({
            "rows": list(subset),
            "labels": [labels[i] for i in subset],
            "omitted_rows": omitted,
            "omitted_labels": [labels[i] for i in omitted],
            "min_singular_value": min_sv,
            "max_singular_value": max_sv,
            "condition": float(max_sv / min_sv) if min_sv else None,
            "residual_inf_on_subset": str(max(abs(F[i]) for i in subset)),
            "omits_active_count": sum(1 for i in omitted if i < len(active)),
            "omits_tie_count": sum(1 for i in omitted if i >= len(active)),
        })
    rows.sort(key=lambda row: row["min_singular_value"], reverse=True)
    payload = {
        "input": args.input,
        "active": list(active),
        "active_triples": [list(TRIPLES[i]) for i in active],
        "ties": list(ties),
        "tie_triples": [list(TRIPLES[i]) for i in ties],
        "equation_count": len(polynomials),
        "require_all_ties": bool(args.require_all_ties),
        "full_residual_inf": str(max_abs(F)),
        "top_bases": rows[: args.top],
        "basis_count": len(rows),
    }
    if args.out:
        Path(args.out).parent.mkdir(parents=True, exist_ok=True)
        Path(args.out).write_text(json.dumps(payload, indent=1) + "\n")

    print("=" * 78)
    print("OVERTIE SQUARE BASIS SCREEN")
    print(f"input: {args.input}")
    print(f"equations: {len(polynomials)}, bases screened: {len(rows)}")
    for row in rows[: min(args.top, 8)]:
        print(
            f"min_sv={row['min_singular_value']:.3e} cond={row['condition']:.3e} "
            f"omitted={row['omitted_labels']}"
        )
    if args.out:
        print(f"wrote {args.out}")
    print("=" * 78)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
