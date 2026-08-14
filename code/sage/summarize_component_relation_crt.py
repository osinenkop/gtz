#!/usr/bin/env python3
"""CRT-lift coordinate relation polynomials from component probe JSON files."""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from sage.all import QQ, ZZ, PolynomialRing, crt


def centered_crt(residues: list[int], primes: list[int]) -> int:
    modulus = ZZ.prod(ZZ(p) for p in primes)
    value = ZZ(crt([ZZ(r) for r in residues], [ZZ(p) for p in primes]))
    if value > modulus // 2:
        value -= modulus
    return int(value)


def rational_reconstruct(residues: list[int], primes: list[int]):
    modulus = ZZ.prod(ZZ(p) for p in primes)
    value = ZZ(crt([ZZ(r) for r in residues], [ZZ(p) for p in primes]))
    try:
        return value.rational_reconstruction(modulus)
    except ArithmeticError:
        return None


def relation_map(component: dict) -> dict[str, dict]:
    return {
        str(row["expression"]): row
        for row in component.get("relations", [])
        if row.get("relation_found")
    }


def summarize_expression(expression: str, records: list[dict], factor: bool) -> dict:
    primes = [int(row["prime"]) for row in records]
    relations = [row["relations"][expression] for row in records]
    degrees = [int(row["relation_degree"]) for row in relations]
    if len(set(degrees)) != 1:
        return {"expression": expression, "ok": False, "reason": "degree_mismatch", "degrees": degrees}
    degree = degrees[0]
    coeff_lists = [[int(c) for c in row["coefficients_low_to_high"]] for row in relations]
    integer_coeffs = []
    rational_coeffs = []
    failures = []
    for exponent in range(degree + 1):
        residues = [coeffs[exponent] for coeffs in coeff_lists]
        integer_coeffs.append(centered_crt(residues, primes))
        reconstructed = rational_reconstruct(residues, primes)
        rational_coeffs.append(reconstructed)
        if reconstructed is None:
            failures.append(exponent)

    ring = PolynomialRing(QQ, "u")
    u = ring.gen()
    integer_poly = sum(QQ(integer_coeffs[i]) * u**i for i in range(degree + 1))
    rational_poly = None
    rational_factorization = []
    if not failures:
        rational_poly = sum(QQ(rational_coeffs[i]) * u**i for i in range(degree + 1))
        if factor:
            rational_factorization = [
                {
                    "degree": int(poly.degree()),
                    "multiplicity": int(multiplicity),
                    "factor": str(poly),
                }
                for poly, multiplicity in rational_poly.factor()
            ]
    successful = [QQ(value) for value in rational_coeffs if value is not None]
    stats = {
        "success_count": len(successful),
        "failure_count": len(failures),
        "failed_exponents": failures,
    }
    if successful:
        stats["max_abs_numerator"] = int(max(abs(value.numerator()) for value in successful))
        stats["max_abs_denominator"] = int(max(abs(value.denominator()) for value in successful))
    return {
        "expression": expression,
        "ok": True,
        "degree": degree,
        "integer_coefficients_low_to_high": integer_coeffs,
        "integer_max_abs_coeff": int(max(abs(c) for c in integer_coeffs)),
        "integer_max_log10_abs_coeff": math.log10(max(1, max(abs(c) for c in integer_coeffs))),
        "integer_lift_polynomial": str(integer_poly),
        "rational_reconstruction": stats,
        "rational_polynomial": None if rational_poly is None else str(rational_poly),
        "rational_factorization": rational_factorization,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", required=True)
    parser.add_argument("--factor-rational", action="store_true")
    parser.add_argument("component_jsons", nargs="+")
    args = parser.parse_args()

    records = []
    for text in args.component_jsons:
        path = Path(text)
        data = json.loads(path.read_text())
        records.append(
            {
                "path": str(path),
                "prime": int(data["characteristic"]),
                "degree": data.get("degree"),
                "augmented_basis_size": data.get("augmented_basis_size"),
                "relations": relation_map(data),
            }
        )
    records.sort(key=lambda row: row["prime"])
    primes = [row["prime"] for row in records]
    if len(set(primes)) != len(primes):
        raise SystemExit(f"duplicate primes: {primes}")

    expression_sets = [set(row["relations"]) for row in records]
    common = sorted(set.intersection(*expression_sets))
    summaries = [summarize_expression(expression, records, args.factor_rational) for expression in common]
    modulus = ZZ.prod(ZZ(p) for p in primes)
    payload = {
        "input_count": len(records),
        "primes": primes,
        "modulus": int(modulus),
        "modulus_log10": math.log10(int(modulus)),
        "inputs": [{k: v for k, v in row.items() if k != "relations"} for row in records],
        "common_expressions": common,
        "summaries": summaries,
    }
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(payload, indent=1) + "\n")

    print("=" * 78)
    print("COMPONENT RELATION CRT")
    print(f"primes: {primes}")
    print(f"common expressions: {common}")
    print(f"modulus log10: {payload['modulus_log10']:.3f}")
    for row in summaries:
        if not row.get("ok"):
            print(f"{row['expression']}: {row['reason']} {row.get('degrees')}")
            continue
        rr = row["rational_reconstruction"]
        total = rr["success_count"] + rr["failure_count"]
        print(
            f"{row['expression']}: degree={row['degree']} "
            f"rr={rr['success_count']}/{total} "
            f"int_log10={row['integer_max_log10_abs_coeff']:.3f} "
            f"{row['rational_polynomial'] or '<no rational lift>'}"
        )
    print(f"wrote {out}")
    print("=" * 78)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
