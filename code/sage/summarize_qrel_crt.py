#!/usr/bin/env python3
"""Summarize and CRT-lift modular q-relation JSON files."""
from __future__ import annotations

import argparse
import json
import math
import re
from pathlib import Path

from sage.all import QQ, ZZ, PolynomialRing, crt


def infer_prime(path: Path, data: dict) -> int:
    characteristic = int(data.get("characteristic") or 0)
    if characteristic:
        return characteristic
    match = re.search(r"_p(\d+)", path.name)
    if not match:
        raise ValueError(f"could not infer prime from {path}")
    return int(match.group(1))


def factor_path(path: Path) -> Path:
    return path.with_name(path.stem + "_factors" + path.suffix)


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


def read_record(path: Path) -> dict:
    data = json.loads(path.read_text())
    coeffs = [int(c) for c in data["coefficients_low_to_high"]]
    prime = infer_prime(path, data)
    factors = None
    fpath = factor_path(path)
    if fpath.exists():
        fdata = json.loads(fpath.read_text())
        factors = {
            "path": str(fpath),
            "factor_count": int(fdata["factor_count"]),
            "factor_degrees": [int(x) for x in fdata["factor_degrees"]],
        }
    return {
        "path": str(path),
        "prime": prime,
        "degree": int(data["relation_degree"]),
        "quotient_degree": int(data["quotient_degree"]),
        "residual_zero": bool(data["residual_zero"]),
        "coefficients": coeffs,
        "factors": factors,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--variable", default="q")
    parser.add_argument("--out")
    parser.add_argument("--factor-rational", action="store_true")
    parser.add_argument("json_paths", nargs="+")
    args = parser.parse_args()

    records = [read_record(Path(path)) for path in args.json_paths]
    records.sort(key=lambda row: row["prime"])
    primes = [row["prime"] for row in records]
    if len(set(primes)) != len(primes):
        raise SystemExit(f"duplicate primes: {primes}")
    degrees = [row["degree"] for row in records]
    quotient_degrees = [row["quotient_degree"] for row in records]
    if len(set(degrees)) != 1:
        raise SystemExit(f"inconsistent q-relation degrees: {degrees}")
    if len(set(quotient_degrees)) != 1:
        raise SystemExit(f"inconsistent quotient degrees: {quotient_degrees}")

    degree = degrees[0]
    coeff_lists = [row["coefficients"] for row in records]
    if any(len(coeffs) != degree + 1 for coeffs in coeff_lists):
        raise SystemExit("coefficient length does not match relation degree")

    modulus = ZZ.prod(ZZ(p) for p in primes)
    integer_coeffs = []
    rational_coeffs = []
    rational_failures = []
    for exponent in range(degree + 1):
        residues = [coeffs[exponent] for coeffs in coeff_lists]
        integer_coeffs.append(centered_crt(residues, primes))
        reconstructed = rational_reconstruct(residues, primes)
        rational_coeffs.append(reconstructed)
        if reconstructed is None:
            rational_failures.append(exponent)

    rational_stats = {
        "success_count": len(rational_coeffs) - len(rational_failures),
        "failure_count": len(rational_failures),
        "failed_exponents_sample": rational_failures[:30],
    }
    successful = [x for x in rational_coeffs if x is not None]
    if successful:
        rational_stats["max_abs_numerator"] = int(max(abs(QQ(x).numerator()) for x in successful))
        rational_stats["max_abs_denominator"] = int(max(abs(QQ(x).denominator()) for x in successful))
        rational_stats["max_log10_abs_numerator"] = math.log10(
            max(1, rational_stats["max_abs_numerator"])
        )
        rational_stats["max_log10_abs_denominator"] = math.log10(
            max(1, rational_stats["max_abs_denominator"])
        )

    factor_summaries = [
        {
            "prime": row["prime"],
            "factor_count": None if row["factors"] is None else row["factors"]["factor_count"],
            "factor_degrees": None if row["factors"] is None else row["factors"]["factor_degrees"],
        }
        for row in records
    ]

    rational_poly = None
    rational_factorization = []
    if not rational_failures:
        ring = PolynomialRing(QQ, args.variable)
        q = ring.gen()
        rational_poly = sum(QQ(rational_coeffs[i]) * q**i for i in range(degree + 1))
        if args.factor_rational:
            rational_factorization = [
                {
                    "degree": int(factor.degree()),
                    "multiplicity": int(multiplicity),
                    "factor": str(factor),
                }
                for factor, multiplicity in rational_poly.factor()
            ]

    payload = {
        "variable": args.variable,
        "input_count": len(records),
        "primes": primes,
        "modulus": int(modulus),
        "modulus_log10": math.log10(int(modulus)),
        "quotient_degree": quotient_degrees[0],
        "relation_degree": degree,
        "all_residuals_zero": all(row["residual_zero"] for row in records),
        "integer_coefficients_low_to_high": integer_coeffs,
        "integer_max_abs_coeff": max(abs(c) for c in integer_coeffs),
        "integer_max_log10_abs_coeff": math.log10(max(1, max(abs(c) for c in integer_coeffs))),
        "rational_reconstruction": rational_stats,
        "rational_polynomial": None if rational_poly is None else str(rational_poly),
        "rational_factorization": rational_factorization,
        "factor_summaries": factor_summaries,
        "inputs": [{k: v for k, v in row.items() if k != "coefficients"} for row in records],
    }

    print("=" * 78)
    print("QREL CRT SUMMARY")
    print(f"inputs: {len(records)}")
    print(f"primes: {primes}")
    print(f"quotient degree: {quotient_degrees[0]}")
    print(f"relation degree: {degree}")
    print(f"modulus log10: {payload['modulus_log10']:.3f}")
    print(
        "rational reconstruction: "
        f"{rational_stats['success_count']} success, {rational_stats['failure_count']} fail"
    )
    print(f"integer max log10 coeff: {payload['integer_max_log10_abs_coeff']:.3f}")
    if args.out:
        out = Path(args.out)
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(json.dumps(payload, indent=1) + "\n")
        print(f"wrote {out}")
    print("=" * 78)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
