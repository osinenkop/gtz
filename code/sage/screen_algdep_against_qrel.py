#!/usr/bin/env python3
"""Test numerical algdep q-candidates against exact modular q-relations."""
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

from sage.all import GF, PolynomialRing


def qrel_poly(path: Path, variable: str):
    data = json.loads(path.read_text())
    prime = int(data["characteristic"])
    ring = PolynomialRing(GF(prime), variable)
    q = ring.gen()
    coeffs = [int(c) for c in data["coefficients_low_to_high"]]
    poly = sum(ring(c) * q**i for i, c in enumerate(coeffs))
    return prime, ring, poly


def candidate_poly(text: str, ring, variable: str):
    text = re.sub(r"\bx\b", variable, text)
    return ring(text.replace("^", "**"))


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--algdep-json", required=True)
    parser.add_argument("--coordinate", default="q")
    parser.add_argument("--best", type=int, default=12)
    parser.add_argument("--variable", default="q")
    parser.add_argument("--out", required=True)
    parser.add_argument("qrel_jsons", nargs="+")
    args = parser.parse_args()

    algdep = json.loads(Path(args.algdep_json).read_text())
    candidates = algdep["results"][args.coordinate]["best"][: args.best]
    qrels = [qrel_poly(Path(path), args.variable) for path in args.qrel_jsons]

    rows = []
    for index, candidate in enumerate(candidates):
        gcd_rows = []
        for prime, ring, qpoly in qrels:
            cpoly = candidate_poly(candidate["polynomial"], ring, args.variable)
            gcd = qpoly.gcd(cpoly)
            gcd_rows.append(
                {
                    "prime": prime,
                    "gcd_degree": int(gcd.degree()) if gcd else -1,
                    "gcd": str(gcd),
                }
            )
        rows.append(
            {
                "candidate_index": index,
                "degree_bound": candidate["degree_bound"],
                "polynomial_degree": candidate["polynomial_degree"],
                "height": candidate["height"],
                "residual_digits": candidate["residual_digits"],
                "gcds": gcd_rows,
                "stable_nontrivial": all(row["gcd_degree"] > 0 for row in gcd_rows),
                "polynomial": candidate["polynomial"],
            }
        )

    payload = {
        "algdep_json": args.algdep_json,
        "coordinate": args.coordinate,
        "qrel_jsons": args.qrel_jsons,
        "candidate_count": len(rows),
        "rows": rows,
    }
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(payload, indent=1) + "\n")

    print("=" * 78)
    print("ALGDEP VS QREL SCREEN")
    print(f"algdep: {args.algdep_json}")
    print(f"qrels: {len(args.qrel_jsons)}")
    for row in rows:
        degrees = [entry["gcd_degree"] for entry in row["gcds"]]
        print(
            f"candidate {row['candidate_index']}: "
            f"polydeg={row['polynomial_degree']} gcd_degrees={degrees}"
        )
    print(f"wrote {out}")
    print("=" * 78)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
