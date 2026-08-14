#!/usr/bin/env python3
"""Canonicalize active triple-constraint patterns from relaxed target sweeps.

The relaxed minimax sweeps in ``v40_relaxed_target_sweep.py`` store the best
and top feasible SLSQP candidates.  This helper quotients their active
low/high triple constraints by row permutations.  The output is diagnostic: it
identifies which KKT active-pattern orbits are worth attacking exactly.
"""

from __future__ import annotations

import argparse
import itertools
import json
import os
from collections import Counter, defaultdict
from pathlib import Path


TRIPLES = list(itertools.combinations(range(6), 3))
TRIPLE_INDEX = {triple: i for i, triple in enumerate(TRIPLES)}
PERMUTATIONS = list(itertools.permutations(range(6)))


def constraint_code(kind: str, triple_index: int) -> int:
    offset = 0 if kind == "low" else len(TRIPLES)
    return offset + triple_index


def decode_constraint(code: int) -> tuple[str, tuple[int, int, int]]:
    if code < len(TRIPLES):
        return "low", TRIPLES[code]
    return "high", TRIPLES[code - len(TRIPLES)]


def permute_code(code: int, permutation: tuple[int, ...]) -> int:
    kind, triple = decode_constraint(code)
    image = tuple(sorted(permutation[i] for i in triple))
    return constraint_code(kind, TRIPLE_INDEX[image])


def canonical_pattern(codes: set[int]) -> tuple[int, ...]:
    return min(tuple(sorted(permute_code(code, perm) for code in codes)) for perm in PERMUTATIONS)


def pattern_labels(pattern: tuple[int, ...]) -> list[str]:
    labels = []
    for code in pattern:
        kind, triple = decode_constraint(code)
        prefix = "L" if kind == "low" else "H"
        labels.append(prefix + "".join(str(i) for i in triple))
    return labels


def row_pattern(row: dict) -> tuple[int, ...]:
    codes = {
        constraint_code(item["kind"], int(item["index"]))
        for item in row.get("active_triple_constraints", [])
    }
    return canonical_pattern(codes)


def cluster_key(value: float, ndigits: int = 10) -> str:
    return f"{value:.{ndigits}f}"


def summarize_file(path: Path) -> dict:
    data = json.load(open(path))
    targets = []
    for target in data.get("targets", []):
        rows = list(target.get("top", []))
        best = target.get("best")
        if best is not None:
            rows.append(best)

        by_pattern: dict[tuple[int, ...], list[dict]] = defaultdict(list)
        by_s = Counter()
        for row in rows:
            pattern = row_pattern(row)
            by_pattern[pattern].append(row)
            by_s[cluster_key(float(row["s"]))] += 1

        pattern_rows = []
        for pattern, members in sorted(
            by_pattern.items(),
            key=lambda item: (min(float(row["s"]) for row in item[1]), -len(item[1])),
        ):
            best_member = min(members, key=lambda row: float(row["s"]))
            pattern_rows.append({
                "count": len(members),
                "active_count": len(pattern),
                "best_s": float(best_member["s"]),
                "best_bound_value": float(best_member["bound_value"]),
                "best_label": best_member["label"],
                "labels": pattern_labels(pattern),
            })

        targets.append({
            "target": target.get("target"),
            "n_rows_seen": len(rows),
            "s_clusters": dict(sorted(by_s.items())),
            "patterns": pattern_rows,
        })

    return {
        "input": str(path),
        "targets": targets,
    }


def write_json(path: str, payload: dict) -> None:
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    tmp = path + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(payload, fh, indent=1)
    os.replace(tmp, path)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("inputs", nargs="+")
    parser.add_argument("--out", default="verify/out/v42_relaxed_active_patterns.json")
    args = parser.parse_args()

    summaries = [summarize_file(Path(path)) for path in args.inputs]
    payload = {"summaries": summaries}
    write_json(args.out, payload)

    print("=" * 78)
    print("RELAXED ACTIVE-PATTERN SUMMARY")
    for summary in summaries:
        print(f"\n{summary['input']}")
        for target in summary["targets"]:
            print(f"  target={target['target']} rows_seen={target['n_rows_seen']}")
            print(f"  s clusters: {target['s_clusters']}")
            for idx, pattern in enumerate(target["patterns"][:8], 1):
                print(
                    f"    #{idx}: count={pattern['count']} active={pattern['active_count']} "
                    f"best_s={pattern['best_s']:.12f} label={pattern['best_label']}"
                )
    print(f"\nwrote {args.out}")
    print("=" * 78)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
