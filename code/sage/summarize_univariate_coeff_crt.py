#!/usr/bin/env python3
"""CRT-lift univariate coefficient-list JSON files."""
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


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--variable", default=None)
    parser.add_argument("--factor-rational", action="store_true")
    parser.add_argument("--out", required=True)
    parser.add_argument("json_paths", nargs="+")
    args = parser.parse_args()

    records = []
    for text in args.json_paths:
        path = Path(text)
        data = json.loads(path.read_text())
        coeffs = [int(c) for c in data["coefficients_low_to_high"]]
        variable = args.variable or data.get("variable", "t")
        records.append(
            {
                "path": str(path),
                "prime": infer_prime(path, data),
                "degree": len(coeffs) - 1,
                "variable": variable,
                "coefficients": coeffs,
            }
        )
    records.sort(key=lambda row: row["prime"])

    primes = [row["prime"] for row in records]
    degrees = [row["degree"] for row in records]
    variables = [row["variable"] for row in records]
    if len(set(primes)) != len(primes):
        raise SystemExit(f"duplicate primes: {primes}")
    if len(set(degrees)) != 1:
        raise SystemExit(f"inconsistent degrees: {degrees}")
    if len(set(variables)) != 1:
        raise SystemExit(f"inconsistent variables: {variables}")

    degree = degrees[0]
    variable = variables[0]
    coeff_lists = [row["coefficients"] for row in records]
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

    modulus = ZZ.prod(ZZ(p) for p in primes)
    ring = PolynomialRing(QQ, variable)
    t = ring.gen()
    integer_poly = sum(QQ(integer_coeffs[i]) * t**i for i in range(degree + 1))
    rational_poly = None
    rational_factorization = []
    if not rational_failures:
        rational_poly = sum(QQ(rational_coeffs[i]) * t**i for i in range(degree + 1))
        if args.factor_rational:
            rational_factorization = [
                {
                    "degree": int(factor.degree()),
                    "multiplicity": int(multiplicity),
                    "factor": str(factor),
                }
                for factor, multiplicity in rational_poly.factor()
            ]

    successful = [QQ(x) for x in rational_coeffs if x is not None]
    rational_stats = {
        "success_count": len(successful),
        "failure_count": len(rational_failures),
        "failed_exponents": rational_failures,
    }
    if successful:
        rational_stats["max_abs_numerator"] = int(max(abs(x.numerator()) for x in successful))
        rational_stats["max_abs_denominator"] = int(max(abs(x.denominator()) for x in successful))
        rational_stats["max_log10_abs_numerator"] = math.log10(
            max(1, rational_stats["max_abs_numerator"])
        )
        rational_stats["max_log10_abs_denominator"] = math.log10(
            max(1, rational_stats["max_abs_denominator"])
        )

    payload = {
        "variable": variable,
        "input_count": len(records),
        "primes": primes,
        "modulus": int(modulus),
        "modulus_log10": math.log10(int(modulus)),
        "degree": degree,
        "integer_coefficients_low_to_high": integer_coeffs,
        "integer_max_abs_coeff": max(abs(c) for c in integer_coeffs),
        "integer_max_log10_abs_coeff": math.log10(max(1, max(abs(c) for c in integer_coeffs))),
        "integer_lift_polynomial": str(integer_poly),
        "rational_reconstruction": rational_stats,
        "rational_polynomial": None if rational_poly is None else str(rational_poly),
        "rational_factorization": rational_factorization,
        "inputs": [{k: v for k, v in row.items() if k != "coefficients"} for row in records],
    }

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(payload, indent=1) + "\n")

    print("=" * 78)
    print("UNIVARIATE COEFFICIENT CRT")
    print(f"inputs: {len(records)}")
    print(f"primes: {primes}")
    print(f"degree: {degree}")
    print(f"modulus log10: {payload['modulus_log10']:.3f}")
    print(
        "rational reconstruction: "
        f"{rational_stats['success_count']} success, {rational_stats['failure_count']} fail"
    )
    print(f"integer max log10 coeff: {payload['integer_max_log10_abs_coeff']:.3f}")
    print(f"wrote {out}")
    print("=" * 78)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
