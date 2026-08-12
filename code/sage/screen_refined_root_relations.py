#!/usr/bin/env python3
"""Screen a refined root for small integer relations.

This is a diagnostic tool for turning a high-precision numerical root into a
small exact ansatz.  It looks for PSLQ/PARI ``lindep`` relations among
``q,z0,...,z8`` and, optionally, among quadratic monomials in those variables.
"""
from __future__ import annotations

import argparse
import itertools
import json
import math
import sys
from pathlib import Path

from sage.all import RealField, ZZ, pari


def parse_real_field(data, precision):
    rf = RealField(precision)
    values = {"q": rf(str(data["q"]))}
    if "margin" in data:
        values["margin"] = rf(str(data["margin"]))
    for i, text in enumerate(data.get("z", [])):
        values[f"z{i}"] = rf(str(text))
    return rf, values


def relation_from_pari(numbers):
    try:
        relation = pari(numbers).lindep()
    except Exception:
        return None
    try:
        coeffs = [ZZ(c) for c in relation]
    except Exception:
        return None
    if not coeffs or all(c == 0 for c in coeffs):
        return None
    gcd = ZZ(0)
    for coeff in coeffs:
        gcd = gcd.gcd(ZZ(abs(coeff)))
    if gcd > 1:
        coeffs = [coeff // gcd for coeff in coeffs]
    first = next((coeff for coeff in coeffs if coeff != 0), ZZ(1))
    if first < 0:
        coeffs = [-coeff for coeff in coeffs]
    return coeffs


def relation_stats(coeffs, numbers):
    residual = abs(sum(coeff * number for coeff, number in zip(coeffs, numbers)))
    height = max(abs(int(coeff)) for coeff in coeffs)
    if residual == 0:
        digits = math.inf
    else:
        digits = float(-residual.log10())
    return height, residual, digits


def format_relation(labels, coeffs):
    parts = []
    for label, coeff in zip(labels, coeffs):
        if coeff == 0:
            continue
        sign = "+" if coeff > 0 else "-"
        mag = abs(int(coeff))
        term = label if mag == 1 else f"{mag}*{label}"
        if not parts:
            parts.append(term if coeff > 0 else f"-{term}")
        else:
            parts.append(f"{sign} {term}")
    return " ".join(parts) + " = 0"


def screen_relations(named_values, max_terms, height_cap, min_digits, include_constant):
    rows = []
    seen = set()
    labels_all = [name for name, _ in named_values]
    values_all = [value for _, value in named_values]

    for r in range(1, max_terms + 1):
        for indices in itertools.combinations(range(len(values_all)), r):
            labels = [labels_all[i] for i in indices]
            numbers = [values_all[i] for i in indices]
            if include_constant:
                labels = ["1"] + labels
                numbers = [numbers[0].parent()(1)] + numbers
            coeffs = relation_from_pari(numbers)
            if coeffs is None:
                continue
            supported = [
                (label, coeff, number)
                for label, coeff, number in zip(labels, coeffs, numbers)
                if coeff != 0
            ]
            if len(supported) < 2:
                continue
            labels = [label for label, _, _ in supported]
            coeffs = [coeff for _, coeff, _ in supported]
            numbers = [number for _, _, number in supported]
            key = (tuple(labels), tuple(int(coeff) for coeff in coeffs))
            if key in seen:
                continue
            seen.add(key)
            height, residual, digits = relation_stats(coeffs, numbers)
            if height > height_cap or digits < min_digits:
                continue
            rows.append(
                {
                    "labels": labels,
                    "coefficients": [int(c) for c in coeffs],
                    "height": int(height),
                    "residual": str(residual),
                    "residual_digits": "inf" if digits == math.inf else digits,
                    "relation": format_relation(labels, coeffs),
                }
            )
    rows.sort(
        key=lambda row: (
            len(row["labels"]),
            row["height"],
            -math.inf if row["residual_digits"] == "inf" else -float(row["residual_digits"]),
            row["relation"],
        )
    )
    return rows


def build_monomials(values, variable_names, degree):
    rf = next(iter(values.values())).parent()
    monomials = [("1", rf(1))]
    if degree >= 1:
        monomials.extend((name, values[name]) for name in variable_names)
    if degree >= 2:
        for i, name_i in enumerate(variable_names):
            for name_j in variable_names[i:]:
                monomials.append((f"{name_i}*{name_j}", values[name_i] * values[name_j]))
    return monomials


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("inputs", nargs="+")
    parser.add_argument("--precision", type=int, default=1200)
    parser.add_argument("--max-linear-terms", type=int, default=4)
    parser.add_argument("--max-monomial-terms", type=int, default=0)
    parser.add_argument("--monomial-degree", type=int, default=2)
    parser.add_argument("--height-cap", type=int, default=1000000)
    parser.add_argument("--min-digits", type=float, default=80.0)
    parser.add_argument("--limit", type=int, default=80)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    payload = {
        "precision_bits": args.precision,
        "max_linear_terms": args.max_linear_terms,
        "max_monomial_terms": args.max_monomial_terms,
        "monomial_degree": args.monomial_degree,
        "height_cap": args.height_cap,
        "min_digits": args.min_digits,
        "inputs": [],
    }

    print("=" * 78)
    print("REFINED ROOT INTEGER RELATION SCREEN")
    print(f"precision: {args.precision} bits")
    print(f"height cap: {args.height_cap}")
    print(f"min residual digits: {args.min_digits:g}")

    for input_name in args.inputs:
        path = Path(input_name)
        data = json.loads(path.read_text())
        _, values = parse_real_field(data, args.precision)
        base_names = ["q"] + [f"z{i}" for i in range(len(data.get("z", [])))]
        named_base = [(name, values[name]) for name in base_names]

        homogeneous = screen_relations(
            named_base,
            args.max_linear_terms,
            args.height_cap,
            args.min_digits,
            include_constant=False,
        )
        affine = screen_relations(
            named_base,
            args.max_linear_terms,
            args.height_cap,
            args.min_digits,
            include_constant=True,
        )

        monomial_rows = []
        if args.max_monomial_terms:
            monomials = build_monomials(values, base_names, args.monomial_degree)
            # The monomial list already contains the constant.
            monomial_rows = screen_relations(
                monomials,
                args.max_monomial_terms,
                args.height_cap,
                args.min_digits,
                include_constant=False,
            )

        item = {
            "input": str(path),
            "q_decimal_30": data.get("q_decimal_30"),
            "margin_decimal_30": data.get("margin_decimal_30"),
            "active": data.get("active"),
            "ties": data.get("ties"),
            "homogeneous_linear": homogeneous[: args.limit],
            "affine_linear": affine[: args.limit],
            "monomial_relations": monomial_rows[: args.limit],
            "counts": {
                "homogeneous_linear": len(homogeneous),
                "affine_linear": len(affine),
                "monomial_relations": len(monomial_rows),
            },
        }
        payload["inputs"].append(item)

        print("-" * 78)
        print(path)
        print(f"q={data.get('q_decimal_30')} margin={data.get('margin_decimal_30')}")
        print(f"homogeneous linear: {len(homogeneous)}")
        for row in homogeneous[: min(10, args.limit)]:
            print(f"  H={row['height']:<8} digits={row['residual_digits']}  {row['relation']}")
        print(f"affine linear: {len(affine)}")
        for row in affine[: min(10, args.limit)]:
            print(f"  H={row['height']:<8} digits={row['residual_digits']}  {row['relation']}")
        if args.max_monomial_terms:
            print(f"monomial relations: {len(monomial_rows)}")
            for row in monomial_rows[: min(10, args.limit)]:
                print(f"  H={row['height']:<8} digits={row['residual_digits']}  {row['relation']}")

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(payload, indent=1) + "\n")
    print("-" * 78)
    print(f"wrote {out}")
    print("=" * 78)
    return 0


if __name__ == "__main__":
    sys.exit(main())
