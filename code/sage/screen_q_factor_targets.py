#!/usr/bin/env python3
"""Screen q-factor target degrees after removing known rational fibers.

The q-eliminant factors very differently modulo different primes.  Before
trying an expensive factor-product CRT recovery for a fixed degree, this helper
counts how many modular factor products of each total degree remain after
discarding already-certified rational factors such as q=1 and q=5.
"""
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path


def parse_prime(path: Path, data: dict) -> int:
    characteristic = int(data.get("characteristic") or 0)
    if characteristic:
        return characteristic
    match = re.search(r"_p(\d+)", path.name)
    if not match:
        raise ValueError(f"cannot infer prime from {path}")
    return int(match.group(1))


def parse_int_list(text: str) -> tuple[int, ...]:
    values = []
    for part in text.split(","):
        part = part.strip()
        if part:
            values.append(int(part))
    return tuple(values)


def monic_linear_root_mod_prime(poly_text: str, prime: int) -> int | None:
    poly = poly_text.replace(" ", "")
    if poly == "q":
        return 0
    match = re.fullmatch(r"q\+([0-9]+)", poly)
    if match:
        return (-int(match.group(1))) % prime
    match = re.fullmatch(r"q-([0-9]+)", poly)
    if match:
        return int(match.group(1)) % prime
    return None


def filtered_degrees(data: dict, prime: int, excluded_roots: tuple[int, ...]) -> tuple[list[int], list[dict]]:
    excluded = {root % prime for root in excluded_roots}
    degrees: list[int] = []
    removed: list[dict] = []
    for index, item in enumerate(data["factors"]):
        degree = int(item["degree"])
        multiplicity = int(item["multiplicity"])
        root = monic_linear_root_mod_prime(str(item["polynomial"]), prime) if degree == 1 else None
        if root is not None and root in excluded:
            removed.append(
                {
                    "index": index,
                    "root": int(root),
                    "degree": degree,
                    "multiplicity": multiplicity,
                    "polynomial": item["polynomial"],
                }
            )
            continue
        degrees.extend([degree] * multiplicity)
    return degrees, removed


def counts_up_to(degrees: list[int], max_degree: int, cap: int) -> dict[int, int]:
    dp = {0: 1}
    for degree in degrees:
        next_dp = dict(dp)
        for total, count in dp.items():
            next_total = total + degree
            if next_total <= max_degree:
                next_dp[next_total] = min(cap, next_dp.get(next_total, 0) + count)
        dp = next_dp
    return {degree: dp.get(degree, 0) for degree in range(1, max_degree + 1)}


def score_summary(row_counts: list[tuple[int, int]], cap: int) -> tuple[int, int, int, int]:
    positive = [count for _prime, count in row_counts if count > 0]
    if not positive:
        return (10**18, 10**18, 10**18, len(row_counts))
    capped = sum(count >= cap for _prime, count in row_counts)
    return (
        capped,
        max(positive),
        sum(positive),
        -len(positive),
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-degree", type=int, default=120)
    parser.add_argument("--cap", type=int, default=20000)
    parser.add_argument("--exclude-rational-roots", default="1,5")
    parser.add_argument("--top", type=int, default=30)
    parser.add_argument("--out", required=True)
    parser.add_argument("factor_jsons", nargs="+")
    args = parser.parse_args()

    excluded_roots = parse_int_list(args.exclude_rational_roots)
    rows = []
    for path_text in args.factor_jsons:
        path = Path(path_text)
        data = json.loads(path.read_text())
        prime = parse_prime(path, data)
        degrees, removed = filtered_degrees(data, prime, excluded_roots)
        counts = counts_up_to(degrees, args.max_degree, args.cap)
        rows.append(
            {
                "path": str(path),
                "prime": prime,
                "remaining_factor_count": len(degrees),
                "remaining_total_degree": sum(degrees),
                "removed": removed,
                "counts": {str(degree): counts[degree] for degree in range(1, args.max_degree + 1)},
            }
        )
    rows.sort(key=lambda row: row["prime"])

    summaries = []
    for degree in range(1, args.max_degree + 1):
        values = [(int(row["prime"]), int(row["counts"][str(degree)])) for row in rows]
        positive = [(prime, count) for prime, count in values if count > 0]
        summaries.append(
            {
                "degree": degree,
                "positive_primes": len(positive),
                "zero_primes": len(values) - len(positive),
                "min_positive_count": min((count for _prime, count in positive), default=0),
                "max_count": max((count for _prime, count in values), default=0),
                "capped_primes": sum(count >= args.cap for _prime, count in values),
                "sum_count": sum(count for _prime, count in values),
                "best_primes": [
                    {"prime": prime, "count": count}
                    for prime, count in sorted(positive, key=lambda item: item[1])[:10]
                ],
                "worst_primes": [
                    {"prime": prime, "count": count}
                    for prime, count in sorted(positive, key=lambda item: item[1], reverse=True)[:10]
                ],
                "score": list(score_summary(values, args.cap)),
            }
        )
    ranked = sorted(summaries, key=lambda item: tuple(item["score"]))
    payload = {
        "input_count": len(rows),
        "max_degree": args.max_degree,
        "cap": args.cap,
        "excluded_rational_roots": list(excluded_roots),
        "ranked_degrees": ranked,
        "rows": rows,
    }
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(payload, indent=1) + "\n")

    print("=" * 78)
    print("Q-FACTOR TARGET SCREEN")
    print(f"inputs: {len(rows)}")
    print(f"excluded rational roots: {list(excluded_roots)}")
    print(f"max degree: {args.max_degree}")
    print(f"cap: {args.cap}")
    print("best target degrees:")
    for item in ranked[: args.top]:
        best = ", ".join(
            f"{entry['prime']}:{entry['count']}" for entry in item["best_primes"][:5]
        )
        print(
            f"  D={item['degree']:3d} pos={item['positive_primes']:2d}/{len(rows)} "
            f"zero={item['zero_primes']:2d} max={item['max_count']:5d} "
            f"capped={item['capped_primes']:2d} sum={item['sum_count']:7d} "
            f"best={best}"
        )
    print(f"wrote {out}")
    print("=" * 78)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
