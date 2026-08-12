#!/usr/bin/env python3
"""Search bounded-degree local separators for over-tied determinant equations.

For a selected square basis f_i and an omitted equation h, solve the homogeneous
linear system

    s h^r = sum_i a_i f_i,

where deg(a_i) <= D and deg(s) <= S.  A solution with s nonzero at the refined
root is a local-component separator candidate.  Over QQ it is an exact
certificate; over a finite field it is modular evidence that should be lifted.
"""
from __future__ import annotations

import argparse
import json
import sys
from math import comb
from pathlib import Path

from sage.all import GF, Matrix, QQ, RealField

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from gtz63_semialgebraic import TRIPLES  # noqa: E402
from probe_bounded_ideal_membership import (  # noqa: E402
    coefficient_field,
    exponent_tuples,
    monomial_from_exp,
    parse_rows,
)
from refine_margin_tie_root import build_system  # noqa: E402


def label(prefix: str, idx: int) -> str:
    return f"{prefix}:{idx}:{''.join(str(a) for a in TRIPLES[idx])}"


def lift_coeff(coeff, characteristic: int):
    if characteristic == 0:
        return QQ(coeff)
    value = int(coeff)
    if value > characteristic // 2:
        value -= characteristic
    return QQ(value)


def support_from_candidate(path: str, candidate_index: int):
    data = json.loads(Path(path).read_text())
    candidates = data["passing_candidates"]
    if candidate_index < 0 or candidate_index >= len(candidates):
        raise SystemExit(f"--support-candidate-index out of range for {len(candidates)} candidates")
    return [term["monomial"] for term in candidates[candidate_index]["separator_terms_prefix"]]


def build_separator_matrix(
    generators,
    target,
    multiplier_degree: int,
    separator_degree: int,
    field,
    separator_support: list[str] | None = None,
):
    ring = target.parent()
    gens = ring.gens()
    multiplier_monomials = [
        monomial_from_exp(ring, gens, exp)
        for exp in exponent_tuples(len(gens), multiplier_degree)
    ]
    if separator_support is None:
        separator_monomials = [
            monomial_from_exp(ring, gens, exp)
            for exp in exponent_tuples(len(gens), separator_degree)
        ]
    else:
        separator_monomials = [ring(monomial) for monomial in separator_support]

    columns = []
    monomial_set = set()
    for i, gen in enumerate(generators):
        for j, monomial in enumerate(multiplier_monomials):
            poly = monomial * gen
            columns.append(("generator", i, j, poly))
            monomial_set.update(poly.dict().keys())
    for j, monomial in enumerate(separator_monomials):
        poly = -monomial * target
        columns.append(("separator", 0, j, poly))
        monomial_set.update(poly.dict().keys())

    monomials = sorted(monomial_set)
    row_index = {monomial: row for row, monomial in enumerate(monomials)}
    entries = {}
    for col, (_kind, _i, _j, poly) in enumerate(columns):
        for monomial, coeff in poly.dict().items():
            if coeff:
                entries[(row_index[monomial], col)] = field(coeff)
    matrix = Matrix(field, len(monomials), len(columns), entries, sparse=True)
    return matrix, columns, multiplier_monomials, separator_monomials


def eval_separator(separator_coeffs, separator_monomials, root, characteristic: int, precision: int):
    RF = RealField(precision)
    value = RF(0)
    for coeff, monomial in zip(separator_coeffs, separator_monomials):
        if not coeff:
            continue
        c = RF(str(lift_coeff(coeff, characteristic)))
        term = RF(1)
        for gen_power, x in zip(monomial.exponents()[0], root):
            if gen_power:
                term *= x ** gen_power
        value += c * term
    return value


def separator_terms(separator_coeffs, separator_monomials, characteristic: int, limit: int):
    terms = []
    for coeff, monomial in zip(separator_coeffs, separator_monomials):
        if coeff:
            terms.append({
                "coefficient": str(lift_coeff(coeff, characteristic)),
                "monomial": str(monomial),
            })
            if len(terms) >= limit:
                break
    return terms


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, help="over-tie refinement JSON")
    parser.add_argument("--rows", required=True, help="selected generator rows")
    parser.add_argument("--target-row", type=int, required=True)
    parser.add_argument("--target-power", type=int, default=1)
    parser.add_argument("--degree", type=int, default=4, help="multiplier degree bound D")
    parser.add_argument("--separator-degree", type=int, default=1, help="separator degree bound S")
    parser.add_argument(
        "--characteristic",
        type=int,
        default=32003,
        help="0 for QQ, otherwise run over GF(p)",
    )
    parser.add_argument("--order", default="degrevlex")
    parser.add_argument("--precision", type=int, default=200)
    parser.add_argument("--eval-threshold-exp", type=int, default=40)
    parser.add_argument("--term-limit", type=int, default=40)
    parser.add_argument(
        "--candidate-limit",
        type=int,
        default=10,
        help="maximum number of candidates to save in the JSON output; use 0 for all",
    )
    parser.add_argument(
        "--support-from",
        default="",
        help="JSON output from this script; restrict separator monomials to one saved candidate support",
    )
    parser.add_argument("--support-candidate-index", type=int, default=0)
    parser.add_argument(
        "--known-generator-rank",
        type=int,
        default=-1,
        help="rank of the generator-multiple block, if already known",
    )
    parser.add_argument(
        "--rank-only",
        action="store_true",
        help="only report whether separator columns cause a rank defect",
    )
    parser.add_argument("--out", default="")
    args = parser.parse_args()

    data = json.loads(Path(args.input).read_text())
    active = tuple(data["active"])
    ties = tuple(data["ties"])
    rows = parse_rows(args.rows)
    ring, _chart, _variables, polynomials, _jacobian = build_system(active, ties, args.order)
    if args.target_row in rows:
        raise SystemExit("--target-row is already among selected rows")
    if min(rows + (args.target_row,)) < 0 or max(rows + (args.target_row,)) >= len(polynomials):
        raise SystemExit(f"row index out of range for {len(polynomials)} equations")
    if args.target_power < 1:
        raise SystemExit("--target-power must be positive")

    labels = [label("active", i) for i in active] + [label("tie", i) for i in ties]
    generators = [polynomials[i] for i in rows]
    target = polynomials[args.target_row] ** args.target_power
    field = coefficient_field(args.characteristic)
    separator_support = None
    if args.support_from:
        separator_support = support_from_candidate(args.support_from, args.support_candidate_index)
    multiplier_count_estimate = comb(len(ring.gens()) + args.degree, args.degree)
    if separator_support is None:
        separator_count_estimate = comb(len(ring.gens()) + args.separator_degree, args.separator_degree)
    else:
        separator_count_estimate = len(separator_support)
    print("=" * 78, flush=True)
    print("LOCAL SEPARATOR MEMBERSHIP", flush=True)
    print(f"input: {args.input}", flush=True)
    print(f"target: row {args.target_row} {labels[args.target_row]}", flush=True)
    print(f"target power: {args.target_power}", flush=True)
    print(f"multiplier degree: {args.degree}", flush=True)
    print(f"separator degree: {args.separator_degree}", flush=True)
    if args.support_from:
        print(
            f"support: candidate {args.support_candidate_index} from {args.support_from}",
            flush=True,
        )
    print(f"field: {'QQ' if args.characteristic == 0 else 'GF(%d)' % args.characteristic}", flush=True)
    print(
        "building sparse matrix: "
        f"{len(generators)} generators x {multiplier_count_estimate} multiplier monomials "
        f"+ {separator_count_estimate} separator monomials "
        f"= {len(generators) * multiplier_count_estimate + separator_count_estimate} columns",
        flush=True,
    )
    matrix, columns, multiplier_monomials, separator_monomials = build_separator_matrix(
        generators, target, args.degree, args.separator_degree, field, separator_support
    )
    separator_offset = len(generators) * len(multiplier_monomials)
    print(f"matrix: {matrix.nrows()} x {matrix.ncols()}", flush=True)
    print("computing ranks...", flush=True)
    matrix_rank = int(matrix.rank())
    if args.known_generator_rank >= 0:
        generator_rank = args.known_generator_rank
    else:
        generator_rank = int(matrix.matrix_from_columns(range(separator_offset)).rank())
    separator_rank_defect = generator_rank + len(separator_monomials) - matrix_rank

    RF = RealField(args.precision)
    root = [RF(str(v)) for v in data["z"]] + [RF(str(data["q"]))]
    threshold = RF(10) ** (-args.eval_threshold_exp)
    candidates = []
    if separator_rank_defect > 0 and not args.rank_only:
        kernel = matrix.right_kernel()
        for basis_index, vec in enumerate(kernel.basis()):
            coeffs = list(vec[separator_offset:])
            nonzero_terms = sum(1 for coeff in coeffs if coeff)
            if not nonzero_terms:
                continue
            value = eval_separator(coeffs, separator_monomials, root, args.characteristic, args.precision)
            candidates.append({
                "basis_index": basis_index,
                "separator_nonzero_terms": nonzero_terms,
                "separator_value": str(value),
                "separator_abs_value": str(abs(value)),
                "passes_eval_threshold": bool(abs(value) > threshold),
                "separator_terms_prefix": separator_terms(
                    coeffs,
                    separator_monomials,
                    args.characteristic,
                    args.term_limit,
                ),
            })

    passing = [candidate for candidate in candidates if candidate["passes_eval_threshold"]]
    if args.candidate_limit < 0:
        raise SystemExit("--candidate-limit must be nonnegative")
    saved_candidate_limit = None if args.candidate_limit == 0 else args.candidate_limit
    top_candidates = sorted(
        candidates,
        key=lambda candidate: RF(candidate["separator_abs_value"]),
        reverse=True,
    )
    best_candidate = top_candidates[0] if top_candidates else None

    payload = {
        "input": args.input,
        "active": list(active),
        "ties": list(ties),
        "selected_rows": list(rows),
        "selected_labels": [labels[i] for i in rows],
        "target_row": args.target_row,
        "target_label": labels[args.target_row],
        "target_power": args.target_power,
        "multiplier_degree": args.degree,
        "separator_degree": args.separator_degree,
        "support_from": args.support_from,
        "support_candidate_index": args.support_candidate_index if args.support_from else None,
        "characteristic": args.characteristic,
        "generator_count": len(generators),
        "multiplier_monomial_count": len(multiplier_monomials),
        "separator_monomial_count": len(separator_monomials),
        "unknown_count": matrix.ncols(),
        "monomial_equation_count": matrix.nrows(),
        "matrix_rank": int(matrix_rank),
        "generator_rank": int(generator_rank),
        "separator_rank_defect": int(separator_rank_defect),
        "kernel_dimension": int(matrix.ncols() - matrix_rank),
        "rank_only": bool(args.rank_only),
        "candidate_count": len(candidates),
        "passing_candidate_count": len(passing),
        "saved_candidate_limit": args.candidate_limit,
        "max_separator_abs_value": best_candidate["separator_abs_value"] if best_candidate else None,
        "max_separator_abs_basis_index": best_candidate["basis_index"] if best_candidate else None,
        "passing_candidates": passing[:saved_candidate_limit],
        "top_candidates_by_abs": top_candidates[:saved_candidate_limit],
        "all_candidate_prefix": candidates[:saved_candidate_limit],
    }
    if args.out:
        Path(args.out).parent.mkdir(parents=True, exist_ok=True)
        Path(args.out).write_text(json.dumps(payload, indent=1) + "\n")

    print(
        f"rank: {payload['matrix_rank']}, generator rank: {payload['generator_rank']}, "
        f"separator rank defect: {payload['separator_rank_defect']}"
    )
    print(f"kernel dimension: {payload['kernel_dimension']}")
    print(f"separator candidates: {len(candidates)}, passing at root: {len(passing)}")
    if passing:
        best = passing[0]
        print(f"first passing basis index: {best['basis_index']}")
        print(f"|s(root)|: {best['separator_abs_value']}")
    elif best_candidate:
        print(f"largest |s(root)| basis index: {best_candidate['basis_index']}")
        print(f"largest |s(root)|: {best_candidate['separator_abs_value']}")
    if args.out:
        print(f"wrote {args.out}")
    print("=" * 78)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
