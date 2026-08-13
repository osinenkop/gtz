#!/usr/bin/env python3
"""Count modular factor-product subsets by target degree."""
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path


def parse_csv_ints(text: str) -> list[int]:
    values = []
    for part in text.split(","):
        part = part.strip()
        if part:
            values.append(int(part))
    if not values:
        raise argparse.ArgumentTypeError("expected at least one integer")
    return values


def infer_prime(path: Path, data: dict) -> int:
    characteristic = int(data.get("characteristic") or 0)
    if characteristic:
        return characteristic
    match = re.search(r"_p(\d+)", path.name)
    if not match:
        raise ValueError(f"cannot infer prime from {path}")
    return int(match.group(1))


def expanded_degrees(data: dict) -> list[int]:
    degrees = []
    for item in data["factors"]:
        degrees.extend([int(item["degree"])] * int(item["multiplicity"]))
    return degrees


def product_counts(degrees: list[int], targets: list[int], cap: int) -> dict[int, int]:
    max_target = max(targets)
    dp = {0: 1}
    for degree in degrees:
        next_dp = dict(dp)
        for total, count in dp.items():
            next_total = total + degree
            if next_total <= max_target:
                next_dp[next_total] = min(cap, next_dp.get(next_total, 0) + count)
        dp = next_dp
    return {target: dp.get(target, 0) for target in targets}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--targets", type=parse_csv_ints, required=True)
    parser.add_argument("--cap", type=int, default=20000)
    parser.add_argument("--out", required=True)
    parser.add_argument("factor_jsons", nargs="+")
    args = parser.parse_args()

    rows = []
    for path_text in args.factor_jsons:
        path = Path(path_text)
        data = json.loads(path.read_text())
        degrees = expanded_degrees(data)
        counts = product_counts(degrees, args.targets, args.cap)
        rows.append(
            {
                "path": str(path),
                "prime": infer_prime(path, data),
                "expanded_factor_count": len(degrees),
                "target_counts": {str(target): counts[target] for target in args.targets},
            }
        )
    rows.sort(key=lambda row: row["prime"])

    summary = {}
    for target in args.targets:
        values = [(row["prime"], row["target_counts"][str(target)]) for row in rows]
        summary[str(target)] = {
            "min_count": min(count for _, count in values),
            "max_count": max(count for _, count in values),
            "capped_count": sum(count >= args.cap for _, count in values),
            "best_primes": [
                {"prime": prime, "count": count}
                for prime, count in sorted(values, key=lambda item: item[1])[:10]
            ],
        }

    payload = {
        "targets": args.targets,
        "cap": args.cap,
        "input_count": len(rows),
        "summary": summary,
        "rows": rows,
    }
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(payload, indent=1) + "\n")

    print("=" * 78)
    print("FACTOR PRODUCT COUNTS")
    print(f"inputs: {len(rows)}")
    print(f"targets: {args.targets}")
    for target in args.targets:
        item = summary[str(target)]
        best = ", ".join(
            f"{row['prime']}:{row['count']}" for row in item["best_primes"][:5]
        )
        print(
            f"degree {target}: min={item['min_count']} max={item['max_count']} "
            f"capped={item['capped_count']}/{len(rows)} best={best}"
        )
    print(f"wrote {out}")
    print("=" * 78)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
