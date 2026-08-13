#!/usr/bin/env python3
"""Checkpointed parallel repeated-support CRT ranking for screened forms."""
from __future__ import annotations

import argparse
import concurrent.futures
import json
import os
from pathlib import Path

from rank_repeated_support_candidates import (
    crt_stats,
    load_labels,
    relation_polynomial,
    repeated_support,
)


SCREENS = None
PRIMES = None
LABELS = None
TARGET_DEGREE = None
MIN_MULTIPLICITY = None
VARIABLE = None


def init_worker(screens, primes, labels, target_degree, min_multiplicity, variable):
    global SCREENS, PRIMES, LABELS, TARGET_DEGREE, MIN_MULTIPLICITY, VARIABLE
    SCREENS = screens
    PRIMES = primes
    LABELS = labels
    TARGET_DEGREE = target_degree
    MIN_MULTIPLICITY = min_multiplicity
    VARIABLE = variable


def rank_trial(trial: int) -> dict:
    rows = [screen["rows"][trial] for screen in SCREENS]
    degrees = [row.get("relation_degree") for row in rows]
    expression = rows[0].get("expression")
    row = {
        "trial": trial,
        "label": LABELS.get(trial),
        "expression": expression,
        "degrees": degrees,
        "eligible": all(degree == TARGET_DEGREE for degree in degrees),
    }
    if not row["eligible"]:
        return row

    products = []
    for screen, screen_row in zip(SCREENS, rows):
        prime = int(screen["characteristic"])
        _, variable, poly = relation_polynomial(screen_row, prime, VARIABLE)
        product = repeated_support(poly, variable, MIN_MULTIPLICITY)
        product["prime"] = prime
        products.append(product)
    row["products"] = products
    row["crt"] = crt_stats(products, PRIMES)
    return row


def ranked_trials(candidates: list[dict | None]) -> list[dict]:
    ranked = [
        row
        for row in candidates
        if row and row.get("eligible") and row.get("crt", {}).get("ok")
    ]
    ranked.sort(
        key=lambda row: (
            -row["crt"]["rational_reconstruction"]["success_count"],
            row["crt"]["integer_max_log10_abs_coeff"],
            row["trial"],
        )
    )
    return ranked


def write_payload(
    path: Path,
    args,
    screen_jsons: list[str],
    primes: list[int],
    candidates: list[dict | None],
) -> None:
    complete_candidates = [row for row in candidates if row is not None]
    ranked = ranked_trials(candidates)
    payload = {
        "screen_jsons": screen_jsons,
        "primes": primes,
        "target_degree": args.target_degree,
        "min_multiplicity": args.min_multiplicity,
        "variable": args.variable,
        "workers": args.workers,
        "completed_trials": sorted(row["trial"] for row in complete_candidates),
        "ranked_trials": [row["trial"] for row in ranked],
        "candidates": candidates,
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(payload, indent=1) + "\n")
    tmp.replace(path)


def load_checkpoint(path: Path, row_count: int) -> list[dict | None]:
    candidates: list[dict | None] = [None] * row_count
    if not path.exists():
        return candidates
    payload = json.loads(path.read_text())
    for row in payload.get("candidates", []):
        if not row:
            continue
        trial = int(row["trial"])
        if 0 <= trial < row_count:
            candidates[trial] = row
    return candidates


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--target-degree", type=int, default=2800)
    parser.add_argument("--min-multiplicity", type=int, default=2)
    parser.add_argument("--variable", default="t")
    parser.add_argument("--labels-json")
    parser.add_argument("--out", required=True)
    parser.add_argument("--checkpoint-out")
    parser.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 2) // 2))
    parser.add_argument("--checkpoint-every", type=int, default=1)
    parser.add_argument("screen_jsons", nargs="+")
    args = parser.parse_args()

    screens = [json.loads(Path(path).read_text()) for path in args.screen_jsons]
    primes = [int(screen["characteristic"]) for screen in screens]
    if len(set(primes)) != len(primes):
        raise SystemExit(f"duplicate primes: {primes}")
    row_counts = {len(screen["rows"]) for screen in screens}
    if len(row_counts) != 1:
        raise SystemExit(f"inconsistent row counts: {sorted(row_counts)}")
    row_count = next(iter(row_counts))

    labels = load_labels(args.labels_json)
    out = Path(args.out)
    checkpoint = Path(args.checkpoint_out) if args.checkpoint_out else out
    candidates = load_checkpoint(checkpoint, row_count)
    pending = [trial for trial, row in enumerate(candidates) if row is None]

    print("=" * 78, flush=True)
    print("PARALLEL REPEATED-SUPPORT CRT RANKING", flush=True)
    print(f"primes: {primes}", flush=True)
    print(f"workers: {args.workers}", flush=True)
    print(f"pending: {len(pending)} / {row_count}", flush=True)

    completed_since_checkpoint = 0
    with concurrent.futures.ProcessPoolExecutor(
        max_workers=args.workers,
        initializer=init_worker,
        initargs=(
            screens,
            primes,
            labels,
            args.target_degree,
            args.min_multiplicity,
            args.variable,
        ),
    ) as executor:
        futures = {executor.submit(rank_trial, trial): trial for trial in pending}
        for future in concurrent.futures.as_completed(futures):
            row = future.result()
            trial = int(row["trial"])
            candidates[trial] = row
            completed_since_checkpoint += 1

            if row.get("eligible") and row.get("crt", {}).get("ok"):
                rr = row["crt"]["rational_reconstruction"]
                total = rr["success_count"] + rr["failure_count"]
                print(
                    f"trial {trial:2d}: eligible rr={rr['success_count']:2d}/{total:2d} "
                    f"support_deg={row['crt']['degree']:2d}",
                    flush=True,
                )
            elif row.get("eligible"):
                print(f"trial {trial:2d}: eligible crt_mismatch", flush=True)
            else:
                print(f"trial {trial:2d}: ineligible degrees={row['degrees']}", flush=True)

            if completed_since_checkpoint >= args.checkpoint_every:
                write_payload(checkpoint, args, args.screen_jsons, primes, candidates)
                completed_since_checkpoint = 0

    write_payload(out, args, args.screen_jsons, primes, candidates)
    ranked = ranked_trials(candidates)
    print("-" * 78, flush=True)
    eligible = sum(1 for row in candidates if row and row.get("eligible"))
    print(f"eligible: {eligible} / {len(candidates)}", flush=True)
    for row in ranked[:20]:
        crt = row["crt"]
        rr = crt["rational_reconstruction"]
        label = row["label"] or "-"
        print(
            f"trial {row['trial']:2d} {label:18s} "
            f"support_deg={crt['degree']:2d} "
            f"rr={rr['success_count']:2d}/{rr['success_count'] + rr['failure_count']:2d} "
            f"int_log10={crt['integer_max_log10_abs_coeff']:.3f} "
            f"{row['expression']}",
            flush=True,
        )
    print(f"wrote {out}", flush=True)
    print("=" * 78, flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
