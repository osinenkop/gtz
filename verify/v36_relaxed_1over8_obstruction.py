#!/usr/bin/env python3
"""Numerical feasibility probe for the reduced relaxed GTZ(6,3) obstruction.

The analytic one-row/two-row deletion reductions show that a counterexample to
``F(P) >= 1/8`` must satisfy the following closed necessary conditions:

    3/8 <= ell_i <= 5/8                         for all rows i,
    lambda_min(P_ij) <= 1/2 <= lambda_max(P_ij)  for all pairs {i,j},
    lambda_min(P_T) <= 1/8 and lambda_max(P_T) >= 7/8
                                                    for all triples T.

This script tries to solve those inequalities on Gr(3,6), with
``P = A A^T`` and ``A in St(6,3)`` imposed by polar retraction.  It is only a
numerical diagnostic.  A robust positive floor would identify a much smaller
semialgebraic certificate target than the original sharp GTZ problem; a near-zero
point would identify the active pattern to certify or refute exactly.
"""

from __future__ import annotations

import argparse
import itertools
import json
import os
from multiprocessing import Pool
from typing import Iterable

import numpy as np
from scipy.optimize import minimize


TRIPLES = list(itertools.combinations(range(6), 3))
PAIRS = list(itertools.combinations(range(6), 2))
T = 1.0 / 8.0
LEV_LOW = 3.0 / 8.0
LEV_HIGH = 5.0 / 8.0
PAIR_CUT = 1.0 / 2.0
TRI_HIGH = 1.0 - T


def retract(x: np.ndarray) -> np.ndarray:
    """Nearest 6x3 Stiefel point by polar retraction."""
    u, _, vt = np.linalg.svd(x, full_matrices=False)
    return u @ vt


def projector(v: np.ndarray) -> np.ndarray:
    a = retract(v.reshape(6, 3))
    return a @ a.T


def eigvals_block(p: np.ndarray, idx: Iterable[int]) -> np.ndarray:
    block = p[np.ix_(idx, idx)]
    return np.linalg.eigvalsh(block)


def violations_from_p(p: np.ndarray) -> list[dict]:
    out: list[dict] = []
    ell = np.diag(p)
    for i, value in enumerate(ell):
        out.append({
            "kind": "leverage_low",
            "index": [i],
            "violation": float(max(LEV_LOW - value, 0.0)),
            "value": float(value),
        })
        out.append({
            "kind": "leverage_high",
            "index": [i],
            "violation": float(max(value - LEV_HIGH, 0.0)),
            "value": float(value),
        })

    for pair in PAIRS:
        ev = eigvals_block(p, pair)
        out.append({
            "kind": "pair_low",
            "index": list(pair),
            "violation": float(max(ev[0] - PAIR_CUT, 0.0)),
            "value": float(ev[0]),
        })
        out.append({
            "kind": "pair_high",
            "index": list(pair),
            "violation": float(max(PAIR_CUT - ev[-1], 0.0)),
            "value": float(ev[-1]),
        })

    for triple in TRIPLES:
        ev = eigvals_block(p, triple)
        out.append({
            "kind": "triple_low",
            "index": list(triple),
            "violation": float(max(ev[0] - T, 0.0)),
            "value": float(ev[0]),
        })
        out.append({
            "kind": "triple_high",
            "index": list(triple),
            "violation": float(max(TRI_HIGH - ev[-1], 0.0)),
            "value": float(ev[-1]),
        })
    return out


def objective_parts(v: np.ndarray) -> tuple[float, float, np.ndarray]:
    p = projector(v)
    vals = np.array([row["violation"] for row in violations_from_p(p)], dtype=float)
    return float(vals.max(initial=0.0)), float(vals @ vals), vals


def soft_objective(v: np.ndarray, beta: float, sum_weight: float) -> float:
    max_v, sumsq, vals = objective_parts(v)
    if beta <= 0:
        return max_v + sum_weight * sumsq
    shifted = vals - max_v
    smooth_max = max_v + float(np.log(np.exp(beta * shifted).sum()) / beta)
    return smooth_max + sum_weight * sumsq


def hard_objective(v: np.ndarray, sum_weight: float) -> float:
    max_v, sumsq, _ = objective_parts(v)
    return max_v + sum_weight * sumsq


def summarize(v: np.ndarray, seed: int, source: str) -> dict:
    p = projector(v)
    violations = violations_from_p(p)
    violations.sort(key=lambda row: row["violation"], reverse=True)
    tri_eigs = np.array([eigvals_block(p, triple) for triple in TRIPLES])
    pair_eigs = np.array([eigvals_block(p, pair) for pair in PAIRS])
    ell = np.diag(p)
    max_v, sumsq, _ = objective_parts(v)

    h = 2.0 * p - np.eye(6)
    h_tri_eigs = np.array([eigvals_block(h, triple) for triple in TRIPLES])
    complement_pairs = []
    seen = set()
    for triple in TRIPLES:
        comp = tuple(i for i in range(6) if i not in triple)
        key = tuple(sorted((triple, comp)))
        if key in seen:
            continue
        seen.add(key)
        i = TRIPLES.index(triple)
        j = TRIPLES.index(comp)
        complement_pairs.append({
            "triple": list(triple),
            "complement": list(comp),
            "lambda_min_pair": [float(tri_eigs[i, 0]), float(tri_eigs[j, 0])],
            "lambda_max_pair": [float(tri_eigs[i, -1]), float(tri_eigs[j, -1])],
            "both_low_violation": float(max(tri_eigs[i, 0] - T, tri_eigs[j, 0] - T, 0.0)),
            "both_high_violation": float(max(TRI_HIGH - tri_eigs[i, -1],
                                             TRI_HIGH - tri_eigs[j, -1], 0.0)),
            "h_spectra": [
                [float(x) for x in h_tri_eigs[i]],
                [float(x) for x in h_tri_eigs[j]],
            ],
        })
    complement_pairs.sort(
        key=lambda row: max(row["both_low_violation"], row["both_high_violation"]),
        reverse=True,
    )

    return {
        "seed": int(seed),
        "source": source,
        "max_violation": max_v,
        "sum_square_violation": sumsq,
        "F": float(tri_eigs[:, 0].max()),
        "min_triple_lambda_min": float(tri_eigs[:, 0].min()),
        "max_triple_lambda_min": float(tri_eigs[:, 0].max()),
        "min_triple_lambda_max": float(tri_eigs[:, -1].min()),
        "max_triple_lambda_max": float(tri_eigs[:, -1].max()),
        "min_pair_lambda_min": float(pair_eigs[:, 0].min()),
        "max_pair_lambda_min": float(pair_eigs[:, 0].max()),
        "min_pair_lambda_max": float(pair_eigs[:, -1].min()),
        "max_pair_lambda_max": float(pair_eigs[:, -1].max()),
        "min_leverage": float(ell.min()),
        "max_leverage": float(ell.max()),
        "leverages": [float(x) for x in ell],
        "worst_violations": violations[:16],
        "worst_complement_pairs": complement_pairs[:10],
        "v": [float(x) for x in v],
    }


def conference_start() -> np.ndarray:
    c = np.array([
        [0, 1, 1, 1, 1, 1],
        [1, 0, 1, -1, -1, 1],
        [1, 1, 0, 1, -1, -1],
        [1, -1, 1, 0, 1, -1],
        [1, -1, -1, 1, 0, 1],
        [1, 1, -1, -1, 1, 0],
    ], dtype=float)
    h = c / np.sqrt(5.0)
    p = (np.eye(6) + h) / 2.0
    w, u = np.linalg.eigh(p)
    a = u[:, np.argsort(w)[-3:]]
    return a.reshape(-1)


def worker(args: tuple[int, int, float, int, bool]) -> dict:
    seed, maxiter, sum_weight, rounds, use_conference = args
    rng = np.random.default_rng(seed)
    v = conference_start() if use_conference else rng.standard_normal(18)
    if use_conference:
        v = v + 1e-3 * rng.standard_normal(18)

    for beta in (10.0, 40.0, 160.0, 0.0):
        for _ in range(rounds):
            if beta:
                fun = lambda z: soft_objective(z, beta, sum_weight)
            else:
                fun = lambda z: hard_objective(z, sum_weight)
            res = minimize(
                fun,
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
    return summarize(v, seed, "conference" if use_conference else "random")


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
    parser.add_argument("--maxiter", type=int, default=5000)
    parser.add_argument("--rounds", type=int, default=1)
    parser.add_argument("--sum-weight", type=float, default=0.05)
    parser.add_argument("--jobs", type=int, default=max(1, (os.cpu_count() or 4) - 2))
    parser.add_argument("--out", default="verify/out/v36_relaxed_1over8_obstruction.json")
    args = parser.parse_args()

    random_seeds = [
        int(x) for x in np.random.SeedSequence(args.seed).generate_state(args.starts, dtype=np.uint32)
    ]
    conference_seeds = [
        int(x) for x in np.random.SeedSequence(args.seed + 1).generate_state(
            args.conference_starts, dtype=np.uint32,
        )
    ]
    jobs = [(seed, args.maxiter, args.sum_weight, args.rounds, False) for seed in random_seeds]
    jobs.extend((seed, args.maxiter, args.sum_weight, args.rounds, True) for seed in conference_seeds)

    print("=" * 78)
    print("RELAXED 1/8 REDUCED-OBSTRUCTION FEASIBILITY PROBE")
    print(f"starts={args.starts}  conference_starts={args.conference_starts}  jobs={args.jobs}")
    print(f"target: max violation of closed necessary conditions for F(P) <= 1/8")
    print("=" * 78, flush=True)

    rows = []
    with Pool(args.jobs) as pool:
        for done, row in enumerate(pool.imap_unordered(worker, jobs, chunksize=1), 1):
            rows.append(row)
            rows.sort(key=lambda item: item["max_violation"])
            if done % 10 == 0 or done == len(jobs):
                best = rows[0]
                print(
                    f"[{done:4d}/{len(jobs)}] best max_violation="
                    f"{best['max_violation']:.9e}  F={best['F']:.9f}  "
                    f"lev=[{best['min_leverage']:.4f},{best['max_leverage']:.4f}]",
                    flush=True,
                )
                write_json(args.out, {
                    "complete": False,
                    "parameters": vars(args),
                    "best": best,
                    "top": rows[:20],
                })

    rows.sort(key=lambda item: item["max_violation"])
    best = rows[0]
    max_values = np.array([row["max_violation"] for row in rows], dtype=float)
    payload = {
        "complete": True,
        "parameters": vars(args),
        "n_runs": len(rows),
        "max_violation_min": float(max_values.min()),
        "max_violation_p05": float(np.percentile(max_values, 5)),
        "max_violation_median": float(np.median(max_values)),
        "best": best,
        "top": rows[:40],
    }
    write_json(args.out, payload)

    print("\nsummary:")
    print(f"  min max-violation: {payload['max_violation_min']:.12e}")
    print(f"  p05 max-violation: {payload['max_violation_p05']:.12e}")
    print(f"  median violation:  {payload['max_violation_median']:.12e}")
    print(f"  best F:            {best['F']:.12f}")
    print(f"  leverages:         min={best['min_leverage']:.12f} max={best['max_leverage']:.12f}")
    print("  worst constraints:")
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
