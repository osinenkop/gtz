#!/usr/bin/env python3
"""Screen minimal-polynomial degrees of linear forms in a quotient algebra."""
from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from compute_q_power_relation import (  # noqa: E402
    make_ring,
    monomial_key,
    relation_for_columns,
    relation_to_string,
    vector_from_poly,
)


def parse_expressions(items: list[str]) -> list[str]:
    expressions: list[str] = []
    for item in items:
        for part in item.split(","):
            part = part.strip()
            if part:
                expressions.append(part)
    return expressions


def expression_poly(ring, text: str):
    return ring(text.replace("^", "**"))


def relation_for_expression(
    *,
    ring,
    gb,
    normal_basis,
    basis_index,
    expression: str,
    max_degree: int,
    check_every: int,
    store_coefficients: bool,
    max_polynomial_degree: int,
) -> dict:
    field = ring.base_ring()
    elt = expression_poly(ring, expression)
    vectors = []
    powers = []
    current = ring(1)
    found = None
    found_at = None
    last_checked = 0

    for power in range(max_degree + 1):
        powers.append(current)
        vectors.append(vector_from_poly(current, basis_index, field))
        should_check = (
            power == max_degree
            or power + 1 > len(normal_basis)
            or (power > 0 and power % check_every == 0)
        )
        if should_check:
            print(f"{expression}: checking powers <= {power}", flush=True)
            relation = relation_for_columns(field, vectors)
            if relation is not None:
                start = max(last_checked + 1, 1)
                for end_power in range(start, power + 1):
                    refined = relation_for_columns(field, vectors[: end_power + 1])
                    if refined is not None:
                        found = refined
                        found_at = end_power
                        break
                if found is None:
                    found = relation
                    found_at = len(relation) - 1
                break
            last_checked = power
        if power < max_degree:
            current = (current * elt).reduce(gb)

    row = {
        "expression": expression,
        "relation_found": found is not None,
        "relation_degree": int(found_at) if found_at is not None else None,
        "searched_degree": max_degree,
    }
    if found is None:
        return row

    residual = sum(coeff * powers[i] for i, coeff in enumerate(found)).reduce(gb)
    nonzero = [i for i, coeff in enumerate(found) if coeff != 0]
    row.update(
        {
            "residual_zero": residual == 0,
            "residual_terms": len(residual.dict()),
            "nonzero_terms": len(nonzero),
            "density": len(nonzero) / len(found),
            "degree_log10_height_proxy": math.log10(
                max(1, max(abs(int(coeff)) for coeff in found))
            ),
        }
    )
    if store_coefficients:
        row["coefficients_low_to_high"] = [int(coeff) for coeff in found]
    if found_at is not None and found_at <= max_polynomial_degree:
        row["polynomial"] = relation_to_string(found, expression)
    return row


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--basis-json", required=True)
    parser.add_argument(
        "--expression",
        action="append",
        default=[],
        help="Expression to screen. Comma-separated lists are accepted.",
    )
    parser.add_argument("--max-degree", type=int, default=0)
    parser.add_argument("--check-every", type=int, default=100)
    parser.add_argument("--store-coefficients", action="store_true")
    parser.add_argument("--max-polynomial-degree", type=int, default=20)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    data = json.loads(Path(args.basis_json).read_text())
    ring = make_ring(data)
    gb = [ring(poly) for poly in data["basis"]]
    ideal = ring.ideal(gb)
    quotient_degree = int(data.get("degree") or ideal.vector_space_dimension())
    max_degree = args.max_degree or quotient_degree
    expressions = parse_expressions(args.expression)
    if not expressions:
        expressions = list(ring.variable_names())

    print("=" * 78)
    print("LINEAR FORM RELATION SCREEN")
    print(f"basis json: {args.basis_json}")
    print(f"field: {ring.base_ring()}")
    print(f"variables: {ring.variable_names()}")
    print(f"quotient degree: {quotient_degree}")
    print(f"searched max degree: {max_degree}")
    print("computing normal basis...", flush=True)
    normal_basis = list(ideal.normal_basis())
    basis_index = {monomial_key(monomial): i for i, monomial in enumerate(normal_basis)}
    print(f"normal basis size: {len(normal_basis)}")

    rows = []
    for expression in expressions:
        print("-" * 78)
        print(f"expression: {expression}", flush=True)
        row = relation_for_expression(
            ring=ring,
            gb=gb,
            normal_basis=normal_basis,
            basis_index=basis_index,
            expression=expression,
            max_degree=max_degree,
            check_every=args.check_every,
            store_coefficients=args.store_coefficients,
            max_polynomial_degree=args.max_polynomial_degree,
        )
        if row["relation_found"]:
            print(
                f"{expression}: relation degree {row['relation_degree']} "
                f"residual_zero={row['residual_zero']}",
                flush=True,
            )
        else:
            print(f"{expression}: no relation through degree {max_degree}", flush=True)
        rows.append(row)

    payload = {
        "basis_json": args.basis_json,
        "characteristic": int(data["characteristic"]),
        "variables": list(ring.variable_names()),
        "quotient_degree": quotient_degree,
        "normal_basis_size": len(normal_basis),
        "max_degree": max_degree,
        "check_every": args.check_every,
        "rows": rows,
    }
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(payload, indent=1) + "\n")
    print("=" * 78)
    print(f"wrote {out}")
    print("=" * 78)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
