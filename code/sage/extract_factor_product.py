#!/usr/bin/env python3
"""Extract selected factor products from a finite-field factorization JSON."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from sage.all import GF, PolynomialRing


def coeffs_low_to_high(poly, variable):
    return [int(poly.monomial_coefficient(variable**i)) for i in range(poly.degree() + 1)]


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--min-multiplicity", type=int, default=2)
    parser.add_argument(
        "--with-multiplicity",
        action="store_true",
        help="include selected factors with their full multiplicity instead of once",
    )
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    data = json.loads(Path(args.input).read_text())
    prime = int(data["characteristic"])
    variable = data.get("variable", "t")
    ring = PolynomialRing(GF(prime), variable)
    t = ring.gen()

    selected = []
    product = ring.one()
    for index, item in enumerate(data["factors"]):
        multiplicity = int(item["multiplicity"])
        if multiplicity < args.min_multiplicity:
            continue
        factor = ring(str(item["polynomial"]).replace("^", "**"))
        factor = factor / factor.leading_coefficient()
        repeat = multiplicity if args.with_multiplicity else 1
        for _ in range(repeat):
            product *= factor
        selected.append(
            {
                "index": index,
                "degree": int(item["degree"]),
                "multiplicity": multiplicity,
                "repeat_in_product": repeat,
            }
        )

    product = product / product.leading_coefficient()
    payload = {
        "source": args.input,
        "characteristic": prime,
        "variable": variable,
        "selection": {
            "min_multiplicity": args.min_multiplicity,
            "with_multiplicity": args.with_multiplicity,
            "selected_count": len(selected),
            "selected": selected,
        },
        "degree": int(product.degree()),
        "relation_degree": int(product.degree()),
        "coefficients_low_to_high": coeffs_low_to_high(product, t),
        "polynomial": str(product),
    }

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(payload, indent=1) + "\n")
    print("=" * 78)
    print("FACTOR PRODUCT EXTRACTION")
    print(f"input: {args.input}")
    print(f"selected factors: {len(selected)}")
    print(f"degree: {product.degree()}")
    print(f"wrote {out}")
    print("=" * 78)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
