#!/usr/bin/env python3
"""SLSQP polish for the core-constrained relaxed 1/8 minimax problem.

This takes a candidate from ``v37_relaxed_1over8_core_probe.py`` and solves the
explicit constrained minimax problem

    minimise s

subject to

    A^T A = I,
    3/8 <= diag(A A^T) <= 5/8,
    lambda_min(P_ij) <= 1/2 <= lambda_max(P_ij) for all pairs,
    lambda_min(P_T) - 1/8 <= s,
    7/8 - lambda_max(P_T) <= s                  for all triples.

The minimum value ``s`` is the remaining triple obstruction after imposing the
closed deletion core.  Numerical only; the goal is to identify a structured
boundary point for exact follow-up.
"""

from __future__ import annotations

import argparse
import itertools
import json
import os

import numpy as np
from scipy.optimize import minimize

from v36_relaxed_1over8_obstruction import eigvals_block, retract, summarize
from v37_relaxed_1over8_core_probe import split_violations


TRIPLES = list(itertools.combinations(range(6), 3))
PAIRS = list(itertools.combinations(range(6), 2))
T = 1.0 / 8.0
LEV_LOW = 3.0 / 8.0
LEV_HIGH = 5.0 / 8.0
PAIR_CUT = 1.0 / 2.0
TRI_HIGH = 1.0 - T


def p_of(x: np.ndarray) -> np.ndarray:
    a = x[:18].reshape(6, 3)
    return a @ a.T


def stiefel_equalities(x: np.ndarray) -> np.ndarray:
    a = x[:18].reshape(6, 3)
    g = a.T @ a - np.eye(3)
    return np.array([g[i, j] for i in range(3) for j in range(i, 3)])


def inequalities(x: np.ndarray) -> np.ndarray:
    p = p_of(x)
    s = x[-1]
    out: list[float] = []

    ell = np.diag(p)
    out.extend(ell - LEV_LOW)
    out.extend(LEV_HIGH - ell)

    for pair in PAIRS:
        ev = eigvals_block(p, pair)
        out.append(PAIR_CUT - ev[0])
        out.append(ev[-1] - PAIR_CUT)

    for triple in TRIPLES:
        ev = eigvals_block(p, triple)
        out.append(s - (ev[0] - T))
        out.append(s - (TRI_HIGH - ev[-1]))

    return np.array(out, dtype=float)


def load_candidate(path: str, key: str, index: int) -> dict:
    data = json.load(open(path))
    item = data
    for part in key.split("."):
        if isinstance(item, list):
            item = item[int(part)]
        else:
            item = item[part]
    if isinstance(item, list):
        item = item[index]
    return item


def active_triple_constraints(p: np.ndarray, s: float, tol: float) -> list[dict]:
    out: list[dict] = []
    for idx, triple in enumerate(TRIPLES):
        ev = eigvals_block(p, triple)
        low_value = ev[0] - T
        high_value = TRI_HIGH - ev[-1]
        if abs(low_value - s) <= tol:
            out.append({
                "kind": "triple_low",
                "index": idx,
                "triple": list(triple),
                "lambda": float(ev[0]),
                "slack_value": float(low_value),
            })
        if abs(high_value - s) <= tol:
            out.append({
                "kind": "triple_high",
                "index": idx,
                "triple": list(triple),
                "lambda": float(ev[-1]),
                "slack_value": float(high_value),
            })
    return out


def write_json(path: str, payload: dict) -> None:
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    tmp = path + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(payload, fh, indent=1)
    os.replace(tmp, path)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", default="verify/out/v37_relaxed_1over8_core_probe_100.json")
    parser.add_argument("--key", default="best_feasibleish_by_triple")
    parser.add_argument("--index", type=int, default=0)
    parser.add_argument("--maxiter", type=int, default=2000)
    parser.add_argument("--ftol", type=float, default=1e-12)
    parser.add_argument("--active-tol", type=float, default=1e-7)
    parser.add_argument("--out", default="verify/out/v38_relaxed_1over8_slsqp_polish.json")
    args = parser.parse_args()

    candidate = load_candidate(args.input, args.key, args.index)
    a0 = retract(np.array(candidate["v"], dtype=float).reshape(6, 3))
    s0 = float(candidate.get("triple_max_violation", candidate.get("max_violation", 0.1)))
    x0 = np.r_[a0.reshape(-1), s0]

    print("=" * 78)
    print("SLSQP POLISH: CORE-CONSTRAINED RELAXED 1/8 MINIMAX")
    print(f"input={args.input} key={args.key}")
    print(f"initial s={s0:.12e}")
    print(
        f"initial eqnorm={np.linalg.norm(stiefel_equalities(x0)):.3e} "
        f"minineq={inequalities(x0).min():+.3e}",
        flush=True,
    )

    result = minimize(
        lambda x: x[-1],
        x0,
        method="SLSQP",
        constraints=[
            {"type": "eq", "fun": stiefel_equalities},
            {"type": "ineq", "fun": inequalities},
        ],
        bounds=[(None, None)] * 18 + [(-0.2, 0.3)],
        options={"maxiter": args.maxiter, "ftol": args.ftol, "disp": True},
    )

    p = p_of(result.x)
    a = result.x[:18].reshape(6, 3)
    row = summarize(a.reshape(-1), int(candidate.get("seed", -1)), "slsqp-polish")
    core_max, core_sumsq, triple_max, triple_sumsq = split_violations(a.reshape(-1))
    row.update({
        "s": float(result.x[-1]),
        "core_max_violation": core_max,
        "core_sum_square_violation": core_sumsq,
        "triple_max_violation": triple_max,
        "triple_sum_square_violation": triple_sumsq,
        "stiefel_eq_norm": float(np.linalg.norm(stiefel_equalities(result.x))),
        "min_inequality": float(inequalities(result.x).min()),
        "active_triple_constraints": active_triple_constraints(
            p, float(result.x[-1]), args.active_tol,
        ),
    })

    try:
        import sympy as sp

        basis = [sp.sqrt(10), sp.sqrt(17)]
        row["nsimplify"] = {
            "F": str(sp.nsimplify(row["F"], basis, tolerance=1e-12)),
            "s": str(sp.nsimplify(row["s"], basis, tolerance=1e-12)),
            "leverages_sorted": [
                str(sp.nsimplify(x, basis, tolerance=1e-12)) for x in sorted(row["leverages"])
            ],
        }
    except Exception as exc:  # pragma: no cover - diagnostic only
        row["nsimplify_error"] = repr(exc)

    payload = {
        "input": args.input,
        "key": args.key,
        "index": args.index,
        "success": bool(result.success),
        "message": str(result.message),
        "nit": int(result.nit),
        "objective": float(result.fun),
        "candidate": row,
    }
    write_json(args.out, payload)

    print("\nsummary:")
    print(f"  success:          {result.success} ({result.message})")
    print(f"  s:                {row['s']:.15f}")
    print(f"  F:                {row['F']:.15f}")
    print(f"  core max viol.:   {row['core_max_violation']:.3e}")
    print(f"  min inequality:   {row['min_inequality']:+.3e}")
    print(f"  eq norm:          {row['stiefel_eq_norm']:.3e}")
    print(f"  leverages sorted: {np.round(sorted(row['leverages']), 12)}")
    if "nsimplify" in row:
        print(f"  nsimplify F:      {row['nsimplify']['F']}")
        print(f"  nsimplify s:      {row['nsimplify']['s']}")
    print(f"  active triples:   {len(row['active_triple_constraints'])}")
    print(f"\nwrote {args.out}")
    print("=" * 78)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
