#!/usr/bin/env python3
"""Generate reproducible primitive low-height linear-form coefficient banks."""
from __future__ import annotations

import argparse
import json
import math
import random
from pathlib import Path


DEFAULT_VARIABLES = ("a", "b", "c", "d", "x", "y", "e", "q")


def canonicalize(coeffs: list[int]) -> tuple[int, ...] | None:
    if not coeffs or all(coeff == 0 for coeff in coeffs):
        return None
    divisor = 0
    for coeff in coeffs:
        divisor = math.gcd(divisor, abs(coeff))
    if divisor > 1:
        coeffs = [coeff // divisor for coeff in coeffs]
    for coeff in coeffs:
        if coeff < 0:
            coeffs = [-value for value in coeffs]
            break
        if coeff > 0:
            break
    return tuple(coeffs)


def load_candidates(path: Path) -> list[dict]:
    payload = json.loads(path.read_text())
    candidates = payload.get("candidates", payload) if isinstance(payload, dict) else payload
    rows = []
    for index, candidate in enumerate(candidates):
        if isinstance(candidate, dict):
            coeffs = candidate.get("coefficients")
            label = candidate.get("label", f"{path.stem}_{index}")
        else:
            coeffs = candidate
            label = f"{path.stem}_{index}"
        if coeffs is None:
            raise ValueError(f"candidate {index} in {path} has no coefficients")
        rows.append({"label": str(label), "coefficients": [int(c) for c in coeffs]})
    return rows


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", required=True)
    parser.add_argument("--count", type=int, default=96)
    parser.add_argument("--height", type=int, default=6)
    parser.add_argument("--seed", type=int, default=20260814)
    parser.add_argument("--include-json", action="append", default=[])
    parser.add_argument("--variable-order", default=",".join(DEFAULT_VARIABLES))
    parser.add_argument(
        "--allow-zero-coefficients",
        action="store_true",
        help="allow random candidates with zero coordinates",
    )
    args = parser.parse_args()

    variable_order = [part.strip() for part in args.variable_order.split(",") if part.strip()]
    if not variable_order:
        raise SystemExit("empty variable order")
    if args.height < 1:
        raise SystemExit("--height must be positive")
    if args.count < 1:
        raise SystemExit("--count must be positive")

    rng = random.Random(args.seed)
    seen: set[tuple[int, ...]] = set()
    rows: list[dict] = []

    def add(label: str, coeffs: list[int]) -> None:
        if len(coeffs) != len(variable_order):
            raise SystemExit(
                f"{label}: expected {len(variable_order)} coefficients, got {len(coeffs)}"
            )
        canonical = canonicalize(coeffs)
        if canonical is None or canonical in seen:
            return
        seen.add(canonical)
        rows.append(
            {
                "label": label,
                "coefficients": list(canonical),
                "height": max(abs(value) for value in canonical),
            }
        )

    for text in args.include_json:
        for candidate in load_candidates(Path(text)):
            add(str(candidate["label"]), candidate["coefficients"])

    attempts = 0
    max_attempts = max(10000, args.count * 1000)
    while len(rows) < args.count and attempts < max_attempts:
        attempts += 1
        coeffs = [rng.randint(-args.height, args.height) for _ in variable_order]
        if not args.allow_zero_coefficients and any(coeff == 0 for coeff in coeffs):
            continue
        canonical = canonicalize(coeffs)
        if canonical is None or canonical in seen:
            continue
        label = f"rand_h{args.height}_s{args.seed}_{len(rows):04d}"
        add(label, list(canonical))

    if len(rows) < args.count:
        raise SystemExit(f"generated only {len(rows)} candidates after {attempts} attempts")

    payload = {
        "description": (
            "Primitive normalized low-height linear-form candidates. "
            "Variable order is " + ",".join(variable_order) + "."
        ),
        "variable_order": variable_order,
        "seed": args.seed,
        "height": args.height,
        "count": len(rows),
        "include_json": args.include_json,
        "allow_zero_coefficients": args.allow_zero_coefficients,
        "candidates": rows,
    }

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(payload, indent=1) + "\n")
    print("=" * 78)
    print("LOW-HEIGHT LINEAR-FORM BANK")
    print(f"out: {out}")
    print(f"count: {len(rows)}")
    print(f"height: {args.height}")
    print(f"seed: {args.seed}")
    print(f"attempts: {attempts}")
    print("=" * 78)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
