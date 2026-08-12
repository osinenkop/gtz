#!/usr/bin/env python3
"""Recover short integer separators from a stable modular separator row space.

The usual reconstruction script first computes modular RREF bases and then
tries rational reconstruction coordinate by coordinate.  That can fail even for
a stable rational subspace, because RREF coordinates are often much higher
height than a good basis of the lattice itself.

This script instead builds the integer lattice

    {v in Z^n : v mod p lies in the saved separator row space for every p}

using the common modular RREF pivot pattern, then applies LLL and exports the
shortest rows as separator candidates.
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

from sage.all import CRT_list, Integer, Matrix, PolynomialRing, QQ, ZZ, vector

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from reconstruct_separator_subspace import (  # noqa: E402
    VARIABLE_NAMES,
    eval_terms,
    format_monomial,
    load_probe,
    monomial_key,
    row_space_basis,
)


def center_residue(value: Integer, modulus: Integer) -> Integer:
    value = Integer(value) % modulus
    if value > modulus // 2:
        value -= modulus
    return value


def crt_center(residues: list[int], primes: list[int], modulus: Integer) -> Integer:
    value = Integer(CRT_list([Integer(v) for v in residues], [Integer(p) for p in primes]))
    return center_residue(value, modulus)


def vector_terms(row, monomials: list[str], ring) -> list[dict]:
    terms = []
    for coeff, monomial in zip(row, monomials):
        coeff = ZZ(coeff)
        if coeff:
            terms.append(
                {
                    "coefficient": str(coeff),
                    "monomial": format_monomial(ring, monomial),
                }
            )
    return terms


def gcd_of_row(row) -> Integer:
    out = ZZ(0)
    for coeff in row:
        out = math.gcd(int(out), int(coeff))
    return ZZ(abs(out))


def primitive_row(row) -> list[Integer]:
    gcd = gcd_of_row(row)
    if gcd <= 1:
        return [ZZ(coeff) for coeff in row]
    return [ZZ(coeff) // gcd for coeff in row]


def canonical_sign(row: list[Integer]) -> list[Integer]:
    for coeff in row:
        if coeff:
            return [-x for x in row] if coeff < 0 else row
    return row


def row_norm2(row) -> Integer:
    return ZZ(sum(ZZ(coeff) * ZZ(coeff) for coeff in row))


def row_max_abs(row) -> Integer:
    return max((abs(ZZ(coeff)) for coeff in row), default=ZZ(0))


def parse_rational(text: str):
    if "." not in text:
        return QQ(text)
    sign = -1 if text.startswith("-") else 1
    unsigned = text[1:] if text.startswith("-") else text
    whole, frac = unsigned.split(".", 1)
    numerator = ZZ(whole + frac)
    denominator = ZZ(10) ** len(frac)
    return QQ(sign * numerator) / QQ(denominator)


def verify_with_rref(row, modular_bases, pivots, primes) -> bool:
    for basis, prime in zip(modular_bases, primes):
        field = basis[0][0].parent()
        target = [field(coeff) for coeff in row]
        combination = [field(0)] * len(row)
        for basis_row_index, pivot_col in enumerate(pivots):
            scalar = target[pivot_col]
            if not scalar:
                continue
            basis_row = basis[basis_row_index]
            for col, value in enumerate(basis_row):
                combination[col] += scalar * value
        if combination != target:
            return False
    return True


def build_lattice_basis(modular_bases, pivots, primes, modulus: Integer):
    rank = len(pivots)
    ncols = len(modular_bases[0][0])
    pivot_set = set(pivots)
    nonpivots = [col for col in range(ncols) if col not in pivot_set]
    rows = []

    for row_index, pivot_col in enumerate(pivots):
        row = [ZZ(0)] * ncols
        row[pivot_col] = ZZ(1)
        for col in nonpivots:
            residues = [int(basis[row_index][col]) for basis in modular_bases]
            row[col] = ZZ(crt_center(residues, primes, modulus))
        rows.append(row)

    for col in nonpivots:
        row = [ZZ(0)] * ncols
        row[col] = ZZ(modulus)
        rows.append(row)

    if len(rows) != ncols:
        raise SystemExit(f"internal error: expected square lattice basis, got {len(rows)} x {ncols}")
    if rank != len(pivots):
        raise SystemExit("internal error: rank mismatch")
    return rows, nonpivots


def comparable_meta_value(key: str, value):
    if key == "support_from" and not value:
        return ""
    if key in {"input", "support_from"} and value:
        return Path(str(value)).name
    return value


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--inputs", nargs="+", required=True)
    parser.add_argument("--root", default="")
    parser.add_argument("--precision", type=int, default=500)
    parser.add_argument("--out-prefix", required=True)
    parser.add_argument("--max-output", type=int, default=30)
    parser.add_argument("--lll-delta", default="0.99")
    args = parser.parse_args()

    paths = [Path(path) for path in args.inputs]
    probes = [load_probe(path) for path in paths]
    primes = [int(data["characteristic"]) for data in probes]
    if len(set(primes)) != len(primes):
        raise SystemExit("input characteristics must be distinct")
    modulus = Integer(math.prod(primes))

    meta_keys = [
        "input",
        "selected_rows",
        "target_row",
        "target_power",
        "multiplier_degree",
        "separator_degree",
        "support_from",
        "support_candidate_index",
    ]
    baseline = probes[0]
    for data in probes[1:]:
        for key in meta_keys:
            if comparable_meta_value(key, data.get(key)) != comparable_meta_value(key, baseline.get(key)):
                raise SystemExit(f"metadata mismatch for key {key}")

    qq_ring = PolynomialRing(QQ, VARIABLE_NAMES, order="degrevlex")
    monomials = sorted(
        {
            term["monomial"]
            for data in probes
            for candidate in data["passing_candidates"]
            for term in candidate["separator_terms_prefix"]
        },
        key=lambda monomial: monomial_key(qq_ring, monomial),
    )

    modular_bases = []
    pivot_patterns = []
    for data in probes:
        rows, pivots = row_space_basis(data, monomials, qq_ring)
        modular_bases.append(rows)
        pivot_patterns.append(pivots)
        print(
            f"p={data['characteristic']}: saved candidates={len(data['passing_candidates'])}, "
            f"separator row rank={len(rows)}, pivots={[monomials[i] for i in pivots]}",
            flush=True,
        )

    if len(set(pivot_patterns)) != 1:
        raise SystemExit("pivot patterns differ across primes; cannot align row spaces")
    pivots = pivot_patterns[0]
    rank = len(pivots)

    lattice_rows, nonpivots = build_lattice_basis(modular_bases, pivots, primes, modulus)
    lattice = Matrix(ZZ, lattice_rows)
    print("=" * 78, flush=True)
    print("SHORT SEPARATOR LATTICE", flush=True)
    print(f"rank over finite fields: {rank}", flush=True)
    print(f"ambient monomials: {len(monomials)}", flush=True)
    print(f"nonpivots: {len(nonpivots)}", flush=True)
    print(f"CRT modulus: {modulus}", flush=True)
    print(f"lattice shape: {lattice.nrows()} x {lattice.ncols()}", flush=True)
    print("running LLL...", flush=True)
    reduced = lattice.LLL(delta=parse_rational(args.lll_delta))
    print("LLL done", flush=True)

    root_path = Path(args.root) if args.root else None
    seen = set()
    candidates = []
    for i in range(reduced.nrows()):
        row = [ZZ(reduced[i, j]) for j in range(reduced.ncols())]
        if not any(row):
            continue
        row = canonical_sign(primitive_row(row))
        key = tuple(row)
        if key in seen:
            continue
        seen.add(key)
        if not verify_with_rref(row, modular_bases, pivots, primes):
            raise SystemExit(f"LLL row {i} failed modular row-space verification")
        terms = vector_terms(row, monomials, qq_ring)
        value = eval_terms(terms, root_path, args.precision)
        entry = {
            "lll_row_index": i,
            "term_count": len(terms),
            "norm2": str(row_norm2(row)),
            "max_abs_coefficient": str(row_max_abs(row)),
            "terms": terms,
        }
        if value is not None:
            entry["separator_value"] = str(value)
            entry["separator_abs_value"] = str(abs(value))
        candidates.append((row_norm2(row), row_max_abs(row), entry))

    candidates.sort(key=lambda item: (item[0], item[1], item[2]["term_count"]))
    payload_candidates = [entry for _norm, _height, entry in candidates[: args.max_output]]
    payload = {
        "sources": [str(path) for path in paths],
        "characteristics": primes,
        "modulus": str(modulus),
        "input": baseline["input"],
        "selected_rows": baseline["selected_rows"],
        "target_row": baseline["target_row"],
        "target_power": baseline["target_power"],
        "multiplier_degree": baseline["multiplier_degree"],
        "separator_degree": baseline["separator_degree"],
        "support_from": baseline.get("support_from"),
        "support_candidate_index": baseline.get("support_candidate_index"),
        "monomial_count": len(monomials),
        "rank": rank,
        "pivot_monomials": [monomials[i] for i in pivots],
        "lll_delta": args.lll_delta,
        "candidate_count_scanned": len(candidates),
        "short_separators": payload_candidates,
    }

    print("-" * 78, flush=True)
    for index, entry in enumerate(payload_candidates):
        extra = ""
        if "separator_abs_value" in entry:
            extra = f", |s(root)|={entry['separator_abs_value']}"
        print(
            f"candidate {index}: terms={entry['term_count']}, "
            f"max_abs={entry['max_abs_coefficient']}, norm2={entry['norm2']}{extra}",
            flush=True,
        )

    prefix = Path(args.out_prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    summary_path = prefix.with_suffix(".json")
    summary_path.write_text(json.dumps(payload, indent=1) + "\n")
    for index, entry in enumerate(payload_candidates):
        sep_path = prefix.parent / f"{prefix.name}_short{index}.json"
        sep_payload = {
            "source_summary": str(summary_path),
            "short_index": index,
            "lll_row_index": entry["lll_row_index"],
            "terms": entry["terms"],
        }
        sep_path.write_text(json.dumps(sep_payload, indent=1) + "\n")
    print(f"wrote {summary_path}", flush=True)
    print(f"wrote {len(payload_candidates)} short separator files with prefix {prefix}", flush=True)
    print("=" * 78, flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
