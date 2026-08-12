#!/usr/bin/env python3
"""Screen active 2x2 principal minors on a lex determinant slice."""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from sage.all import AA, PolynomialRing, QQ

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from certify_lex_slice_minor import active_minor_polynomial, to_univariate  # noqa: E402
from classify_lex_slice_roots import find_univariate, read_polynomials, solve_shape_relations  # noqa: E402
from gtz63_semialgebraic import TRIPLES  # noqa: E402


def sign_at_real_roots(poly, witness) -> tuple[int, int, int]:
    neg = zero = pos = 0
    for root, _multiplicity in poly.roots(AA, multiplicities=True):
        value = witness(root)
        if value < 0:
            neg += 1
        elif value > 0:
            pos += 1
        else:
            zero += 1
    return neg, zero, pos


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--lex-json", required=True)
    parser.add_argument("--parameter", default="z0")
    parser.add_argument(
        "--factor-mode",
        choices=["residual", "whole"],
        default="residual",
        help="residual removes z^2-5/9; whole screens every root of the parameter polynomial",
    )
    parser.add_argument("--out")
    args = parser.parse_args()

    lex_json_path = Path(args.lex_json)
    data = json.loads(lex_json_path.read_text())
    basis_path = Path(data["lex_basis"]["path"])
    variables = tuple(data["summary"]["variables"])
    active = tuple(data["summary"]["active_indices"])

    multi_ring, polynomials = read_polynomials(basis_path, variables)
    _parameter_index, parameter_poly_multi = find_univariate(polynomials, args.parameter)
    relations = solve_shape_relations(polynomials, multi_ring, args.parameter)

    ring = PolynomialRing(QQ, args.parameter)
    parameter = ring.gen()
    parameter_poly = ring(str(parameter_poly_multi).replace("^", "**"))
    known_quadratic = parameter**2 - QQ(5) / QQ(9)
    if args.factor_mode == "residual":
        if parameter_poly % known_quadratic != 0:
            raise SystemExit("known quadratic does not divide the parameter polynomial")
        screened_factor = parameter_poly // known_quadratic
    else:
        screened_factor = parameter_poly

    z_exprs = []
    for i in range(9):
        name = f"z{i}"
        if name == args.parameter:
            z_exprs.append(parameter)
        else:
            z_exprs.append(to_univariate(relations[name], ring))

    records = []
    for active_index in active:
        for minor_rows in ((0, 1), (0, 2), (1, 2)):
            witness, _d_poly = active_minor_polynomial(
                z_exprs,
                triple_index=active_index,
                minor_rows=minor_rows,
                ring=ring,
            )
            remainder = witness % screened_factor
            neg, zero, pos = sign_at_real_roots(screened_factor, remainder)
            records.append({
                "active_index": active_index,
                "active_triple": list(TRIPLES[active_index]),
                "minor_rows": list(minor_rows),
                "remainder_degree": int(remainder.degree()),
                "remainder": str(remainder) if remainder.degree() <= 4 else None,
                "negative_roots": neg,
                "zero_roots": zero,
                "positive_roots": pos,
                "uniform_negative": bool(neg > 0 and zero == 0 and pos == 0),
                "uniform_positive": bool(pos > 0 and zero == 0 and neg == 0),
            })

    result = {
        "source": str(lex_json_path),
        "parameter": args.parameter,
        "factor_mode": args.factor_mode,
        "screened_factor_degree": int(screened_factor.degree()),
        "screened_factor_real_roots": len(screened_factor.roots(AA, multiplicities=False)),
        "uniform_negative_count": sum(record["uniform_negative"] for record in records),
        "uniform_positive_count": sum(record["uniform_positive"] for record in records),
        "records": records,
    }
    out_path = Path(args.out) if args.out else lex_json_path.with_name(f"{lex_json_path.stem}_minor_screen.json")
    out_path.write_text(json.dumps(result, indent=1) + "\n")
    print(
        "screened",
        len(records),
        "minors;",
        result["uniform_negative_count"],
        "uniform negative",
    )
    print(f"wrote {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
