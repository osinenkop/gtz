#!/usr/bin/env python3
"""Factor a finite-field q-relation JSON produced by compute_q_power_relation."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from sage.all import GF, PolynomialRing


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    data = json.loads(Path(args.input).read_text())
    p = int(data["characteristic"])
    variable = data.get("variable", "q")
    ring = PolynomialRing(GF(p), variable)
    q = ring.gen()
    coeffs = data["coefficients_low_to_high"]
    poly = sum(GF(p)(coeff) * q**i for i, coeff in enumerate(coeffs))
    factorization = poly.factor()
    factors = [
        {
            "degree": int(factor.degree()),
            "multiplicity": int(multiplicity),
            "polynomial": str(factor),
        }
        for factor, multiplicity in factorization
    ]
    payload = {
        "input": args.input,
        "characteristic": p,
        "variable": variable,
        "degree": int(poly.degree()),
        "factor_count": len(factors),
        "factor_degrees": [item["degree"] for item in factors],
        "factor_degree_multiplicities": [
            [item["degree"], item["multiplicity"]] for item in factors
        ],
        "factors": factors,
    }
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(payload, indent=1) + "\n")
    print("=" * 78)
    print("Q-RELATION FACTORIZATION")
    print(f"input: {args.input}")
    print(f"degree: {poly.degree()}")
    print(f"factors: {len(factors)}")
    print([(item["degree"], item["multiplicity"]) for item in factors])
    print(f"wrote {out}")
    print("=" * 78)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
