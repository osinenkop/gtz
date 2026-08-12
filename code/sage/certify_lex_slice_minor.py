#!/usr/bin/env python3
"""Exact algebraic-real sign check for a lex slice minor witness."""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from sage.all import AA, Matrix, PolynomialRing, QQ

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from classify_lex_slice_roots import (  # noqa: E402
    find_univariate,
    read_polynomials,
    solve_shape_relations,
)
from gtz63_semialgebraic import TRIPLES  # noqa: E402


def parse_minor_rows(text: str) -> tuple[int, int]:
    rows = tuple(int(x.strip()) for x in text.split(",") if x.strip())
    if len(rows) != 2 or len(set(rows)) != 2 or any(row < 0 or row >= 3 for row in rows):
        raise SystemExit("--minor-rows must be two distinct local row indices, e.g. 0,2")
    return rows


def to_univariate(expr, ring):
    return ring(str(expr).replace("^", "**"))


def active_minor_polynomial(
    z_exprs: list,
    triple_index: int,
    minor_rows: tuple[int, int],
    ring,
):
    y = Matrix(ring, 6, 3, 0)
    for i in range(3):
        y[i, i] = ring.one()
    for r in range(3):
        for c in range(3):
            y[3 + r, c] = z_exprs[3 * r + c]
    gram = y.transpose() * y
    d = gram.det()
    n = y * gram.adjugate() * y.transpose()
    triple = TRIPLES[triple_index]
    block = Matrix(
        ring,
        2,
        2,
        lambda a, b: 6 * n[triple[minor_rows[a]], triple[minor_rows[b]]]
        - (d if minor_rows[a] == minor_rows[b] else ring.zero()),
    )
    return block.det(), d


def root_signs(poly, witness, precision_digits: int) -> list[dict]:
    roots = poly.roots(AA, multiplicities=True)
    out = []
    for index, (root, multiplicity) in enumerate(roots, start=1):
        value = witness(root)
        sign = -1 if value < 0 else 1 if value > 0 else 0
        out.append({
            "root_index": index,
            "root_approx": root.n(digits=precision_digits).str(no_sci=False),
            "multiplicity": int(multiplicity),
            "witness_sign": sign,
            "witness_approx": value.n(digits=precision_digits).str(no_sci=False),
        })
    return out


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--lex-json", required=True)
    parser.add_argument("--parameter", default="z0")
    parser.add_argument("--active-index", type=int, default=1)
    parser.add_argument("--minor-rows", default="0,2")
    parser.add_argument("--out")
    parser.add_argument("--witness-out")
    parser.add_argument("--precision-digits", type=int, default=20)
    args = parser.parse_args()

    lex_json_path = Path(args.lex_json)
    data = json.loads(lex_json_path.read_text())
    basis_path = Path(data["lex_basis"]["path"])
    variables = tuple(data["summary"]["variables"])

    multi_ring, polynomials = read_polynomials(basis_path, variables)
    _parameter_index, parameter_poly_multi = find_univariate(polynomials, args.parameter)
    relations = solve_shape_relations(polynomials, multi_ring, args.parameter)

    ring = PolynomialRing(QQ, args.parameter)
    parameter = ring.gen()
    parameter_poly = ring(str(parameter_poly_multi).replace("^", "**"))
    known_quadratic = parameter**2 - QQ(5) / QQ(9)
    if parameter_poly % known_quadratic != 0:
        raise SystemExit("known quadratic does not divide the parameter polynomial")
    residual_factor = parameter_poly // known_quadratic

    z_exprs = []
    for i in range(9):
        name = f"z{i}"
        if name == args.parameter:
            z_exprs.append(parameter)
        else:
            z_exprs.append(to_univariate(relations[name], ring))
    known_z_remainders = {
        f"z{i}_minus_{args.parameter}": str((z_exprs[i] - parameter) % known_quadratic)
        for i in range(9)
        if f"z{i}" != args.parameter
    }
    known_u0_remainder = (
        str(to_univariate(relations["u0"], ring) % known_quadratic)
        if "u0" in relations
        else None
    )

    minor_rows = parse_minor_rows(args.minor_rows)
    witness, d_poly = active_minor_polynomial(
        z_exprs,
        triple_index=args.active_index,
        minor_rows=minor_rows,
        ring=ring,
    )
    witness_remainder = witness % residual_factor
    d_remainder = d_poly % residual_factor
    residual_signs = root_signs(residual_factor, witness_remainder, args.precision_digits)
    known_signs = root_signs(known_quadratic, witness % known_quadratic, args.precision_digits)
    d_signs = root_signs(residual_factor, d_remainder, args.precision_digits)

    out_path = Path(args.out) if args.out else lex_json_path.with_name(f"{lex_json_path.stem}_minor_cert.json")
    witness_path = (
        Path(args.witness_out)
        if args.witness_out
        else lex_json_path.with_name(f"{lex_json_path.stem}_minor_{args.active_index}_{''.join(map(str, minor_rows))}.txt")
    )
    witness_path.write_text(str(witness_remainder) + "\n")

    result = {
        "source": str(lex_json_path),
        "basis_path": str(basis_path),
        "parameter": args.parameter,
        "active_index": args.active_index,
        "active_triple": list(TRIPLES[args.active_index]),
        "minor_rows": list(minor_rows),
        "parameter_degree": int(parameter_poly.degree()),
        "residual_factor_degree": int(residual_factor.degree()),
        "residual_factor_real_roots": len(residual_signs),
        "witness_polynomial_degree": int(witness.degree()),
        "witness_remainder_degree": int(witness_remainder.degree()),
        "witness_remainder_path": str(witness_path),
        "all_residual_roots_have_negative_witness": all(item["witness_sign"] < 0 for item in residual_signs),
        "residual_root_signs": residual_signs,
        "known_quadratic_root_signs": known_signs,
        "known_quadratic_z_remainders": known_z_remainders,
        "known_quadratic_u0_remainder": known_u0_remainder,
        "known_quadratic_is_all_equal_z": all(value == "0" for value in known_z_remainders.values()),
        "d_positive_on_residual_roots": all(item["witness_sign"] > 0 for item in d_signs),
        "d_signs_on_residual_roots": d_signs,
    }
    out_path.write_text(json.dumps(result, indent=1) + "\n")
    print(
        "residual roots",
        len(residual_signs),
        "negative witness:",
        result["all_residual_roots_have_negative_witness"],
    )
    print(f"wrote {out_path}")
    print(f"wrote {witness_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
