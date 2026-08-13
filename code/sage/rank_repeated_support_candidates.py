#!/usr/bin/env python3
"""Rank screened linear forms by repeated-support CRT reconstruction quality."""
from __future__ import annotations

import argparse
import collections
import json
import math
from pathlib import Path

from sage.all import GF, QQ, ZZ, PolynomialRing, crt


def coeffs_low_to_high(poly, variable):
    return [int(poly.monomial_coefficient(variable**i)) for i in range(poly.degree() + 1)]


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


def relation_polynomial(row: dict, prime: int, variable: str):
    coeffs = row.get("coefficients_low_to_high")
    if coeffs is None:
        raise ValueError("row has no coefficients_low_to_high")
    ring = PolynomialRing(GF(prime), variable)
    t = ring.gen()
    poly = sum(GF(prime)(coeff) * t**i for i, coeff in enumerate(coeffs))
    return ring, t, poly


def repeated_support(poly, variable, min_multiplicity: int) -> dict:
    ring = poly.parent()
    product = ring.one()
    selected = []
    profile = collections.Counter()
    for factor, multiplicity in poly.factor():
        degree = int(factor.degree())
        multiplicity = int(multiplicity)
        profile[(degree, multiplicity)] += 1
        if multiplicity < min_multiplicity:
            continue
        factor = factor / factor.leading_coefficient()
        product *= factor
        selected.append({"degree": degree, "multiplicity": multiplicity})
    product = product / product.leading_coefficient()
    return {
        "degree": int(product.degree()),
        "coefficients_low_to_high": coeffs_low_to_high(product, variable),
        "selected_count": len(selected),
        "selected": selected,
        "factor_profile": [
            {"degree": degree, "multiplicity": multiplicity, "count": count}
            for (degree, multiplicity), count in sorted(profile.items())
        ],
    }


def crt_stats(products: list[dict], primes: list[int]) -> dict:
    degrees = [int(product["degree"]) for product in products]
    if len(set(degrees)) != 1:
        return {"ok": False, "degrees": degrees}
    degree = degrees[0]
    coeff_lists = [product["coefficients_low_to_high"] for product in products]
    integer_coeffs = []
    rational_coeffs = []
    failed = []
    for exponent in range(degree + 1):
        residues = [coeffs[exponent] for coeffs in coeff_lists]
        integer_coeffs.append(centered_crt(residues, primes))
        reconstructed = rational_reconstruct(residues, primes)
        rational_coeffs.append(reconstructed)
        if reconstructed is None:
            failed.append(exponent)

    successful = [QQ(value) for value in rational_coeffs if value is not None]
    stats = {
        "ok": True,
        "degree": degree,
        "integer_max_abs_coeff": int(max(abs(c) for c in integer_coeffs)),
        "integer_max_log10_abs_coeff": math.log10(max(1, max(abs(c) for c in integer_coeffs))),
        "rational_reconstruction": {
            "success_count": len(successful),
            "failure_count": len(failed),
            "failed_exponents": failed,
        },
    }
    if successful:
        stats["rational_reconstruction"].update(
            {
                "max_abs_numerator": int(max(abs(x.numerator()) for x in successful)),
                "max_abs_denominator": int(max(abs(x.denominator()) for x in successful)),
            }
        )
    return stats


def load_labels(path_text: str | None) -> dict[int, str]:
    if not path_text:
        return {}
    payload = json.loads(Path(path_text).read_text())
    candidates = payload.get("candidates", payload) if isinstance(payload, dict) else payload
    labels = {}
    for index, candidate in enumerate(candidates):
        if isinstance(candidate, dict) and candidate.get("label"):
            labels[index] = str(candidate["label"])
    return labels


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--target-degree", type=int, default=2800)
    parser.add_argument("--min-multiplicity", type=int, default=2)
    parser.add_argument("--variable", default="t")
    parser.add_argument("--labels-json")
    parser.add_argument("--out", required=True)
    parser.add_argument("screen_jsons", nargs="+")
    args = parser.parse_args()

    screens = [json.loads(Path(path).read_text()) for path in args.screen_jsons]
    primes = [int(screen["characteristic"]) for screen in screens]
    if len(set(primes)) != len(primes):
        raise SystemExit(f"duplicate primes: {primes}")
    row_counts = {len(screen["rows"]) for screen in screens}
    if len(row_counts) != 1:
        raise SystemExit(f"inconsistent row counts: {sorted(row_counts)}")

    labels = load_labels(args.labels_json)
    candidates = []
    for trial in range(next(iter(row_counts))):
        rows = [screen["rows"][trial] for screen in screens]
        degrees = [row.get("relation_degree") for row in rows]
        expression = rows[0].get("expression")
        row = {
            "trial": trial,
            "label": labels.get(trial),
            "expression": expression,
            "degrees": degrees,
            "eligible": all(degree == args.target_degree for degree in degrees),
        }
        if not row["eligible"]:
            candidates.append(row)
            continue

        products = []
        for screen, screen_row in zip(screens, rows):
            prime = int(screen["characteristic"])
            _, variable, poly = relation_polynomial(screen_row, prime, args.variable)
            product = repeated_support(poly, variable, args.min_multiplicity)
            product["prime"] = prime
            products.append(product)
        row["products"] = products
        row["crt"] = crt_stats(products, primes)
        candidates.append(row)

    ranked = [
        row
        for row in candidates
        if row.get("eligible") and row.get("crt", {}).get("ok")
    ]
    ranked.sort(
        key=lambda row: (
            -row["crt"]["rational_reconstruction"]["success_count"],
            row["crt"]["integer_max_log10_abs_coeff"],
            row["trial"],
        )
    )

    payload = {
        "screen_jsons": args.screen_jsons,
        "primes": primes,
        "target_degree": args.target_degree,
        "min_multiplicity": args.min_multiplicity,
        "variable": args.variable,
        "ranked_trials": [row["trial"] for row in ranked],
        "candidates": candidates,
    }
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(payload, indent=1) + "\n")

    print("=" * 78)
    print("REPEATED-SUPPORT CRT RANKING")
    print(f"primes: {primes}")
    print(f"eligible: {len(ranked)} / {len(candidates)}")
    for row in ranked[:20]:
        crt = row["crt"]
        rr = crt["rational_reconstruction"]
        label = row["label"] or "-"
        print(
            f"trial {row['trial']:2d} {label:18s} "
            f"support_deg={crt['degree']:2d} "
            f"rr={rr['success_count']:2d}/{rr['success_count'] + rr['failure_count']:2d} "
            f"int_log10={crt['integer_max_log10_abs_coeff']:.3f} "
            f"{row['expression']}"
        )
    print(f"wrote {out}")
    print("=" * 78)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
