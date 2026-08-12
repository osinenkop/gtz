#!/usr/bin/env python3
"""Reconstruct a rational separator subspace from modular separator probes.

The input files are JSON outputs from probe_local_separator_membership.py.  The
script uses only the saved separator parts of passing modular kernel vectors:
it row-reduces those separator vectors over each finite field, checks that the
pivot pattern is stable, then CRT/rational-reconstructs the reduced basis.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from sage.all import CRT_list, GF, Integer, Matrix, PolynomialRing, QQ, RealField


VARIABLE_NAMES = [f"z{i}" for i in range(9)] + ["q"]


def load_probe(path: Path):
    data = json.loads(path.read_text())
    characteristic = int(data["characteristic"])
    if characteristic <= 0:
        raise SystemExit(f"{path}: expected positive characteristic")
    candidates = data.get("passing_candidates", [])
    if not candidates:
        raise SystemExit(f"{path}: no passing candidates saved")
    for candidate in candidates:
        saved = len(candidate.get("separator_terms_prefix", []))
        expected = int(candidate.get("separator_nonzero_terms", saved))
        if saved < expected:
            raise SystemExit(
                f"{path}: candidate {candidate.get('basis_index')} is truncated "
                f"({saved} saved, {expected} nonzero)"
            )
    return data


def monomial_key(ring, monomial: str):
    return ring(monomial).exponents()[0]


def row_space_basis(data, monomials: list[str], ring):
    field = GF(int(data["characteristic"]))
    index = {monomial: i for i, monomial in enumerate(monomials)}
    rows = []
    for candidate in data["passing_candidates"]:
        row = [field(0)] * len(monomials)
        for term in candidate["separator_terms_prefix"]:
            row[index[term["monomial"]]] = field(term["coefficient"])
        rows.append(row)
    matrix = Matrix(field, rows)
    rref = matrix.rref()
    nonzero_rows = []
    for i in range(rref.nrows()):
        row = [rref[i, j] for j in range(rref.ncols())]
        if any(row):
            nonzero_rows.append(row)
    pivots = []
    for row in nonzero_rows:
        pivots.append(next(j for j, coeff in enumerate(row) if coeff))
    return nonzero_rows, tuple(pivots)


def reconstruct_entry(residues, primes, row_index: int, col_index: int, monomial: str):
    modulus = Integer(math.prod(primes))
    value = Integer(CRT_list([Integer(v) for v in residues], [Integer(p) for p in primes]))
    try:
        return QQ(value.rational_reconstruction(modulus))
    except ArithmeticError as exc:
        raise SystemExit(
            "rational reconstruction failed for "
            f"basis row {row_index}, column {col_index}, monomial {monomial}, "
            f"residues={residues}, modulus={modulus}"
        ) from exc


def eval_terms(terms, root_path: Path | None, precision: int):
    if root_path is None:
        return None
    data = json.loads(root_path.read_text())
    root = [RealField(precision)(str(v)) for v in data["z"]]
    root.append(RealField(precision)(str(data["q"])))
    ring = PolynomialRing(QQ, VARIABLE_NAMES, order="degrevlex")
    RF = RealField(precision)
    value = RF(0)
    for term in terms:
        coeff = QQ(term["coefficient"])
        monomial = ring(term["monomial"])
        product = RF(1)
        for exponent, coordinate in zip(monomial.exponents()[0], root):
            if exponent:
                product *= coordinate ** exponent
        value += (RF(coeff.numerator()) / RF(coeff.denominator())) * product
    return value


def format_monomial(ring, monomial: str) -> str:
    return str(PolynomialRing(QQ, VARIABLE_NAMES, order="degrevlex")(monomial))


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
    parser.add_argument("--out-prefix", default="")
    args = parser.parse_args()

    paths = [Path(path) for path in args.inputs]
    probes = [load_probe(path) for path in paths]
    primes = [int(data["characteristic"]) for data in probes]
    if len(set(primes)) != len(primes):
        raise SystemExit("input characteristics must be distinct")

    meta_keys = [
        "input",
        "selected_rows",
        "target_row",
        "target_power",
        "multiplier_degree",
        "separator_degree",
        "support_from",
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
        raise SystemExit("pivot patterns differ across primes; cannot align RREF bases")
    rank = len(modular_bases[0])
    pivots = pivot_patterns[0]

    reconstructed = []
    for row_index in range(rank):
        terms = []
        for col_index, monomial in enumerate(monomials):
            residues = [int(basis[row_index][col_index]) for basis in modular_bases]
            coeff = reconstruct_entry(residues, primes, row_index, col_index, monomial)
            if coeff:
                terms.append({"coefficient": str(coeff), "monomial": format_monomial(qq_ring, monomial)})
        reconstructed.append(terms)

    root_path = Path(args.root) if args.root else None
    root_values = [eval_terms(terms, root_path, args.precision) for terms in reconstructed]
    payload = {
        "sources": [str(path) for path in paths],
        "characteristics": primes,
        "modulus": str(math.prod(primes)),
        "input": baseline["input"],
        "selected_rows": baseline["selected_rows"],
        "target_row": baseline["target_row"],
        "target_power": baseline["target_power"],
        "multiplier_degree": baseline["multiplier_degree"],
        "separator_degree": baseline["separator_degree"],
        "monomial_count": len(monomials),
        "rank": rank,
        "pivot_monomials": [monomials[i] for i in pivots],
        "basis": [],
    }
    print("=" * 78)
    print("RECONSTRUCTED SEPARATOR SUBSPACE")
    print(f"rank: {rank}")
    print(f"monomials: {len(monomials)}")
    for index, terms in enumerate(reconstructed):
        coeffs = [abs(QQ(term["coefficient"])) for term in terms]
        max_height = max(coeffs) if coeffs else QQ(0)
        entry = {
            "basis_index": index,
            "term_count": len(terms),
            "max_abs_coefficient": str(max_height),
            "terms": terms,
        }
        if root_values[index] is not None:
            value = root_values[index]
            entry["separator_value"] = str(value)
            entry["separator_abs_value"] = str(abs(value))
        payload["basis"].append(entry)
        print(
            f"basis {index}: terms={len(terms)}, max_abs_coeff={max_height}"
            + (
                f", |s(root)|={abs(root_values[index])}"
                if root_values[index] is not None
                else ""
            )
        )

    if args.out_prefix:
        prefix = Path(args.out_prefix)
        prefix.parent.mkdir(parents=True, exist_ok=True)
        summary_path = prefix.with_suffix(".json")
        summary_path.write_text(json.dumps(payload, indent=1) + "\n")
        for entry in payload["basis"]:
            sep_path = prefix.parent / f"{prefix.name}_basis{entry['basis_index']}.json"
            sep_payload = {
                "source_summary": str(summary_path),
                "basis_index": entry["basis_index"],
                "terms": entry["terms"],
            }
            sep_path.write_text(json.dumps(sep_payload, indent=1) + "\n")
        print(f"wrote {summary_path}")
        print(f"wrote {rank} separator basis files with prefix {prefix}")
    print("=" * 78)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
