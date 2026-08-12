#!/usr/bin/env python3
"""Reconstruct a small univariate eliminant from modular lex JSON files."""
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

from sage.all import GF, PolynomialRing, QQ, ZZ, crt


def infer_prime(path: Path, data: dict) -> int:
    characteristic = int(data.get("characteristic") or 0)
    if characteristic:
        return characteristic
    match = re.search(r"_p(\d+)_", path.name)
    if not match:
        raise ValueError(f"could not infer prime from {path}")
    return int(match.group(1))


def univariate_text(data: dict, variable: str) -> str:
    for record in data.get("univariate", []):
        if record.get("variable") == variable:
            return str(record["polynomial"])
    for record in data.get("lex_basis", {}).get("univariate_polynomials", []):
        if record.get("variable") == variable:
            return str(record["polynomial"])
    parsed = data.get("run", {}).get("parsed", {})
    poly = parsed.get("lex_first_poly")
    if isinstance(poly, str) and re.search(rf"(?<![A-Za-z0-9_]){re.escape(variable)}(?![A-Za-z0-9_])", poly):
        return poly
    raise ValueError(f"no univariate polynomial in {variable}")


def centered_mod_coeffs(poly_text: str, variable: str, prime: int) -> list[int]:
    ring = PolynomialRing(ZZ, variable)
    poly = ring(poly_text.replace("^", "**"))
    degree = int(poly.degree())
    out = []
    gen = ring.gen()
    for exponent in range(degree + 1):
        coeff = int(poly.monomial_coefficient(gen**exponent))
        out.append(coeff % prime)
    return out


def modular_coeffs(poly_text: str, variable: str, prime: int, divide_factor: str | None) -> list[int]:
    if divide_factor is None:
        return centered_mod_coeffs(poly_text, variable, prime)

    ring = PolynomialRing(GF(prime), variable)
    poly = ring(poly_text.replace("^", "**"))
    factor = ring(divide_factor.replace("^", "**"))
    quotient, remainder = poly.quo_rem(factor)
    if remainder != 0:
        raise ValueError(f"divisor does not divide polynomial modulo {prime}")
    lead = quotient.leading_coefficient()
    if lead != 1:
        quotient = quotient / lead
    gen = ring.gen()
    return [int(quotient.monomial_coefficient(gen**exponent)) for exponent in range(quotient.degree() + 1)]


def reconstruct_coeff(residues: list[int], primes: list[int]) -> int:
    modulus = ZZ.prod(primes)
    value = ZZ(crt([ZZ(r) for r in residues], [ZZ(p) for p in primes]))
    if value > modulus // 2:
        value -= modulus
    return int(value)


def reconstruct_rational_coeff(residues: list[int], primes: list[int]):
    modulus = ZZ.prod(primes)
    value = ZZ(crt([ZZ(r) for r in residues], [ZZ(p) for p in primes]))
    try:
        return value.rational_reconstruction(modulus)
    except ArithmeticError:
        return None


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--variable", default="z0")
    parser.add_argument("--divide-factor")
    parser.add_argument("--out")
    parser.add_argument("json_paths", nargs="+")
    args = parser.parse_args()

    paths = [Path(path) for path in args.json_paths]
    records = []
    primes = []
    degrees = []
    coeff_lists = []
    for path in paths:
        data = json.loads(path.read_text())
        prime = infer_prime(path, data)
        poly_text = univariate_text(data, args.variable)
        coeffs = modular_coeffs(poly_text, args.variable, prime, args.divide_factor)
        records.append({"path": str(path), "prime": prime, "polynomial": poly_text})
        primes.append(prime)
        degrees.append(len(coeffs) - 1)
        coeff_lists.append(coeffs)

    if len(set(primes)) != len(primes):
        raise SystemExit(f"duplicate primes: {primes}")
    if len(set(degrees)) != 1:
        raise SystemExit(f"inconsistent degrees: {degrees}")

    degree = degrees[0]
    coeffs = []
    rational_coeffs = []
    rational_failures = []
    for exponent in range(degree + 1):
        residues = [coeff_list[exponent] for coeff_list in coeff_lists]
        coeffs.append(reconstruct_coeff(residues, primes))
        reconstructed = reconstruct_rational_coeff(residues, primes)
        rational_coeffs.append(reconstructed)
        if reconstructed is None:
            rational_failures.append(exponent)

    ring = PolynomialRing(QQ, args.variable)
    gen = ring.gen()
    integer_lift = sum(QQ(coeffs[i]) * gen**i for i in range(degree + 1))
    rational_lift = None
    rational_factorization = []
    if not rational_failures:
        rational_lift = sum(QQ(rational_coeffs[i]) * gen**i for i in range(degree + 1))
        rational_factorization = [
            {"factor": str(factor), "degree": int(factor.degree()), "multiplicity": int(multiplicity)}
            for factor, multiplicity in rational_lift.factor()
        ]
    integer_factorization = [
        {"factor": str(factor), "degree": int(factor.degree()), "multiplicity": int(multiplicity)}
        for factor, multiplicity in integer_lift.factor()
    ]
    result = {
        "variable": args.variable,
        "input_count": len(paths),
        "primes": primes,
        "modulus": int(ZZ.prod(primes)),
        "degree": degree,
        "integer_coefficients_low_to_high": coeffs,
        "integer_lift_polynomial": str(integer_lift),
        "integer_lift_factorization": integer_factorization,
        "rational_reconstruction_success": not rational_failures,
        "rational_reconstruction_failed_exponents": rational_failures,
        "rational_coefficients_low_to_high": [
            None if coeff is None else str(coeff) for coeff in rational_coeffs
        ],
        "rational_lift_polynomial": None if rational_lift is None else str(rational_lift),
        "rational_lift_factorization": rational_factorization,
        "inputs": records,
    }

    if args.out:
        out_path = Path(args.out)
        out_path.write_text(json.dumps(result, indent=1) + "\n")
        print(f"wrote {out_path}")
    if rational_lift is not None:
        print(result["rational_lift_polynomial"])
        for item in rational_factorization:
            print(f"{item['factor']}  multiplicity={item['multiplicity']}")
    else:
        print("rational reconstruction failed for exponents", rational_failures)
        print(result["integer_lift_polynomial"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
