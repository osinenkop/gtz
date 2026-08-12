#!/usr/bin/env python3
"""Run stronger algdep screens on a saved high-precision refined root."""
from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

from sage.all import RealField


def coeff_height(poly):
    coeffs = poly.coefficients(sparse=False)
    try:
        return max(abs(int(c)) for c in coeffs)
    except TypeError:
        return None


def log10_height(height):
    if height is None or height <= 0:
        return None
    return math.log10(height)


def residual_digits(value):
    if value == 0:
        return "inf"
    return float(-value.log10())


def algdep_rows(x, max_degree):
    rows = []
    seen = set()
    for degree in range(1, max_degree + 1):
        try:
            poly = x.algdep(degree)
        except Exception as exc:  # pragma: no cover - diagnostic path
            rows.append({"degree": degree, "error": repr(exc)})
            continue
        poly_text = str(poly)
        height = coeff_height(poly)
        value = abs(poly(x))
        duplicate = poly_text in seen
        seen.add(poly_text)
        rows.append(
            {
                "degree_bound": degree,
                "polynomial_degree": int(poly.degree()),
                "height": height,
                "log10_height": log10_height(height),
                "residual": str(value),
                "residual_digits": residual_digits(value),
                "duplicate_previous": duplicate,
                "polynomial": poly_text,
            }
        )
    return rows


def best_rows(rows, max_rows):
    usable = [row for row in rows if "error" not in row and not row["duplicate_previous"]]
    usable.sort(
        key=lambda row: (
            -float("inf") if row["residual_digits"] == "inf" else -row["residual_digits"],
            row["polynomial_degree"],
            row["log10_height"] if row["log10_height"] is not None else float("inf"),
        )
    )
    return usable[:max_rows]


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--precision", type=int, default=1400)
    parser.add_argument("--max-degree", type=int, default=64)
    parser.add_argument("--coordinates", default="q,margin")
    parser.add_argument("--best", type=int, default=12)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    data = json.loads(Path(args.input).read_text())
    RF = RealField(args.precision)
    values = {
        "q": RF(str(data["q"])),
        "margin": (RF(str(data["q"])) - 1) / 6,
    }
    for i, text in enumerate(data.get("z", [])):
        values[f"z{i}"] = RF(str(text))

    selected = [part.strip() for part in args.coordinates.split(",") if part.strip()]
    missing = [name for name in selected if name not in values]
    if missing:
        raise SystemExit(f"unknown coordinate(s): {', '.join(missing)}")

    results = {}
    for name in selected:
        rows = algdep_rows(values[name], args.max_degree)
        results[name] = {
            "rows": rows,
            "best": best_rows(rows, args.best),
        }

    payload = {
        "input": args.input,
        "precision_bits": args.precision,
        "max_degree": args.max_degree,
        "coordinates": selected,
        "results": results,
    }

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(payload, indent=1) + "\n")

    print("=" * 78)
    print("REFINED ROOT ALGDEP SCREEN")
    print(f"input: {args.input}")
    print(f"precision: {args.precision}")
    print(f"max degree: {args.max_degree}")
    for name in selected:
        best = results[name]["best"][:3]
        print(f"{name}:")
        for row in best:
            print(
                "  deg<=",
                row["degree_bound"],
                "polydeg=",
                row["polynomial_degree"],
                "log10H=",
                row["log10_height"],
                "digits=",
                row["residual_digits"],
            )
    print(f"wrote {out}")
    print("=" * 78)
    return 0


if __name__ == "__main__":
    sys.exit(main())
