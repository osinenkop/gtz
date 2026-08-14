#!/usr/bin/env python3
"""Core-constrained numerical probe for the relaxed ``F(P) >= 1/8`` target.

``v36_relaxed_1over8_obstruction.py`` minimizes the maximum violation over all
closed necessary conditions for a 1/8 counterexample.  Its best point is the
known seventh extremal, which misses the deletion core by only 1/56 but misses
the triple 1/8 conditions by 1/24.  This companion probe asks a sharper
question:

    Inside the deletion core and pair-straddling region, how small can the
    triple-only obstruction become?

The core and pair constraints are enforced by a large quadratic penalty, while
the objective reports the true max violation separately.  Numerical only.
"""

from __future__ import annotations

import argparse
import json
import os
from multiprocessing import Pool

import numpy as np
from scipy.optimize import minimize

from v36_relaxed_1over8_obstruction import (
    conference_start,
    objective_parts,
    projector,
    summarize,
    violations_from_p,
)


CORE_KINDS = {"leverage_low", "leverage_high", "pair_low", "pair_high"}
TRIPLE_KINDS = {"triple_low", "triple_high"}


def split_violations(v: np.ndarray) -> tuple[float, float, float, float]:
    rows = violations_from_p(projector(v))
    core = np.array([row["violation"] for row in rows if row["kind"] in CORE_KINDS])
    triple = np.array([row["violation"] for row in rows if row["kind"] in TRIPLE_KINDS])
    return (
        float(core.max(initial=0.0)),
        float(core @ core),
        float(triple.max(initial=0.0)),
        float(triple @ triple),
    )


def objective(v: np.ndarray, core_weight: float, triple_sum_weight: float) -> float:
    _, core_sumsq, triple_max, triple_sumsq = split_violations(v)
    return triple_max + triple_sum_weight * triple_sumsq + core_weight * core_sumsq


def run_one(args: tuple[int, bool, int, float, float, int]) -> dict:
    seed, use_conference, maxiter, core_weight, triple_sum_weight, rounds = args
    rng = np.random.default_rng(seed)
    v = conference_start() if use_conference else rng.standard_normal(18)
    if use_conference:
        v = v + 1e-3 * rng.standard_normal(18)

    for _ in range(rounds):
        res = minimize(
            lambda z: objective(z, core_weight, triple_sum_weight),
            v,
            method="Nelder-Mead",
            options={
                "maxiter": maxiter,
                "maxfev": maxiter,
                "xatol": 1e-13,
                "fatol": 1e-15,
            },
        )
        v = res.x
    row = summarize(v, seed, "conference" if use_conference else "random")
    core_max, core_sumsq, triple_max, triple_sumsq = split_violations(v)
    row.update({
        "core_max_violation": core_max,
        "core_sum_square_violation": core_sumsq,
        "triple_max_violation": triple_max,
        "triple_sum_square_violation": triple_sumsq,
        "penalty_objective": objective(v, core_weight, triple_sum_weight),
    })
    return row


def write_json(path: str, payload: dict) -> None:
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    tmp = path + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(payload, fh, indent=1)
    os.replace(tmp, path)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--starts", type=int, default=80)
    parser.add_argument("--conference-starts", type=int, default=20)
    parser.add_argument("--seed", type=int, default=20260815)
    parser.add_argument("--jobs", type=int, default=max(1, (os.cpu_count() or 4) - 2))
    parser.add_argument("--maxiter", type=int, default=6000)
    parser.add_argument("--rounds", type=int, default=3)
    parser.add_argument("--core-weight", type=float, default=1e4)
    parser.add_argument("--triple-sum-weight", type=float, default=0.02)
    parser.add_argument("--out", default="verify/out/v37_relaxed_1over8_core_probe.json")
    args = parser.parse_args()

    random_seeds = [
        int(x) for x in np.random.SeedSequence(args.seed).generate_state(args.starts, dtype=np.uint32)
    ]
    conference_seeds = [
        int(x) for x in np.random.SeedSequence(args.seed + 1).generate_state(
            args.conference_starts, dtype=np.uint32,
        )
    ]
    jobs = [
        (seed, False, args.maxiter, args.core_weight, args.triple_sum_weight, args.rounds)
        for seed in random_seeds
    ]
    jobs.extend(
        (seed, True, args.maxiter, args.core_weight, args.triple_sum_weight, args.rounds)
        for seed in conference_seeds
    )

    print("=" * 78)
    print("CORE-CONSTRAINED RELAXED 1/8 PROBE")
    print(f"starts={args.starts} conference_starts={args.conference_starts} jobs={args.jobs}")
    print(f"core_weight={args.core_weight:g} rounds={args.rounds} maxiter={args.maxiter}")
    print("=" * 78, flush=True)

    rows: list[dict] = []
    with Pool(args.jobs) as pool:
        for done, row in enumerate(pool.imap_unordered(run_one, jobs, chunksize=1), 1):
            rows.append(row)
            rows.sort(key=lambda item: (item["core_max_violation"], item["triple_max_violation"]))
            if done % 10 == 0 or done == len(jobs):
                feasibleish = [r for r in rows if r["core_max_violation"] < 1e-5]
                best_core = rows[0]
                best_triple = min(feasibleish or rows, key=lambda item: item["triple_max_violation"])
                print(
                    f"[{done:4d}/{len(jobs)}] best core={best_core['core_max_violation']:.3e} "
                    f"triple={best_core['triple_max_violation']:.3e}; "
                    f"best feasible-ish triple={best_triple['triple_max_violation']:.3e} "
                    f"core={best_triple['core_max_violation']:.3e}",
                    flush=True,
                )
                write_json(args.out, {
                    "complete": False,
                    "parameters": vars(args),
                    "top_by_core": rows[:20],
                    "top_feasibleish_by_triple": sorted(
                        feasibleish, key=lambda item: item["triple_max_violation"],
                    )[:20],
                })

    rows_by_core = sorted(rows, key=lambda item: (item["core_max_violation"], item["triple_max_violation"]))
    feasibleish = [r for r in rows if r["core_max_violation"] < 1e-5]
    rows_by_feasible_triple = sorted(
        feasibleish or rows,
        key=lambda item: (item["triple_max_violation"], item["core_max_violation"]),
    )
    payload = {
        "complete": True,
        "parameters": vars(args),
        "n_runs": len(rows),
        "n_core_violation_lt_1e-5": len(feasibleish),
        "best_by_core": rows_by_core[0],
        "best_feasibleish_by_triple": rows_by_feasible_triple[0],
        "top_by_core": rows_by_core[:40],
        "top_feasibleish_by_triple": rows_by_feasible_triple[:40],
    }
    write_json(args.out, payload)

    best = rows_by_feasible_triple[0]
    print("\nsummary:")
    print(f"  feasible-ish runs:       {len(feasibleish)} / {len(rows)}")
    print(f"  best core violation:     {rows_by_core[0]['core_max_violation']:.12e}")
    print(f"  best triple violation:   {best['triple_max_violation']:.12e}")
    print(f"  core violation there:    {best['core_max_violation']:.12e}")
    print(f"  F there:                 {best['F']:.12f}")
    print("  worst constraints at best feasible-ish point:")
    for row in best["worst_violations"][:8]:
        print(
            f"    {row['kind']:<14} {row['index']}  "
            f"violation={row['violation']:.9e} value={row['value']:.12f}"
        )
    print(f"\nwrote {args.out}")
    print("=" * 78)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
