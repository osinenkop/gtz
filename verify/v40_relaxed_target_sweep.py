#!/usr/bin/env python3
"""SLSQP sweep for reduced relaxed GTZ(6,3) targets.

For a proposed bound ``F(P) >= t`` with ``0 < t < 1/6``, the one-row and
two-row deletion reductions imply that a counterexample must lie in the closed
reduced obstruction

    1 - 5t <= ell_i <= 5t,
    lambda_max(P_ij) >= 1 - 4t,
    lambda_min(P_ij) <= 4t,
    lambda_min(P_T) <= t and lambda_max(P_T) >= 1 - t.

This script solves the explicit constrained minimax relaxation

    minimise s

subject to the Stiefel equations and the first two lines above, with ``s``
bounding the two triple violations:

    lambda_min(P_T) - t <= s,
    1 - t - lambda_max(P_T) <= s.

If the global minimum has ``s > 0``, then the reduced obstruction is empty for
that target.  This is numerical only, but it identifies which targets and active
patterns are plausible for exact certification.
"""

from __future__ import annotations

import argparse
import itertools
import json
import os
from fractions import Fraction
from pathlib import Path

import numpy as np
from scipy.optimize import minimize

from v17_active_orbits import TTSP_CASES, ttsp_projector
from v36_relaxed_1over8_obstruction import conference_start, eigvals_block, retract


TRIPLES = list(itertools.combinations(range(6), 3))
PAIRS = list(itertools.combinations(range(6), 2))


def parse_target(text: str) -> float:
    return float(Fraction(text)) if "/" in text else float(text)


def p_of(x: np.ndarray) -> np.ndarray:
    a = x[:18].reshape(6, 3)
    return a @ a.T


def stiefel_equalities(x: np.ndarray) -> np.ndarray:
    a = x[:18].reshape(6, 3)
    g = a.T @ a - np.eye(3)
    return np.array([g[i, j] for i in range(3) for j in range(i, 3)])


def inequalities(x: np.ndarray, t: float) -> np.ndarray:
    p = p_of(x)
    s = x[-1]
    lev_low = 1.0 - 5.0 * t
    lev_high = 5.0 * t
    pair_low_cut = 4.0 * t
    pair_high_cut = 1.0 - 4.0 * t
    out: list[float] = []

    ell = np.diag(p)
    out.extend(ell - lev_low)
    out.extend(lev_high - ell)

    for pair in PAIRS:
        ev = eigvals_block(p, pair)
        out.append(pair_low_cut - ev[0])
        out.append(ev[-1] - pair_high_cut)

    for triple in TRIPLES:
        ev = eigvals_block(p, triple)
        out.append(s - (ev[0] - t))
        out.append(s - ((1.0 - t) - ev[-1]))

    return np.array(out, dtype=float)


def triple_violations(p: np.ndarray, t: float) -> np.ndarray:
    vals = []
    for triple in TRIPLES:
        ev = eigvals_block(p, triple)
        vals.append(ev[0] - t)
        vals.append((1.0 - t) - ev[-1])
    return np.array(vals, dtype=float)


def summarize(x: np.ndarray, t: float, label: str) -> dict:
    p = p_of(x)
    tri = np.array([eigvals_block(p, triple) for triple in TRIPLES])
    pair = np.array([eigvals_block(p, pair) for pair in PAIRS])
    ell = np.diag(p)
    ineq = inequalities(x, t)
    triple_vals = triple_violations(p, t)
    active = []
    for idx, triple in enumerate(TRIPLES):
        ev = tri[idx]
        low = ev[0] - t
        high = (1.0 - t) - ev[-1]
        if abs(low - x[-1]) <= 1e-6:
            active.append({"kind": "low", "index": idx, "triple": list(triple), "value": float(ev[0])})
        if abs(high - x[-1]) <= 1e-6:
            active.append({"kind": "high", "index": idx, "triple": list(triple), "value": float(ev[-1])})
    return {
        "label": label,
        "target": t,
        "s": float(x[-1]),
        "bound_value": float(t + x[-1]),
        "success_gap": float(x[-1]),
        "F": float(tri[:, 0].max()),
        "eq_norm": float(np.linalg.norm(stiefel_equalities(x))),
        "min_inequality": float(ineq.min()),
        "max_triple_violation": float(triple_vals.max(initial=-np.inf)),
        "min_leverage": float(ell.min()),
        "max_leverage": float(ell.max()),
        "leverages": [float(v) for v in ell],
        "min_pair_lambda_min": float(pair[:, 0].min()),
        "max_pair_lambda_min": float(pair[:, 0].max()),
        "min_pair_lambda_max": float(pair[:, -1].min()),
        "max_pair_lambda_max": float(pair[:, -1].max()),
        "min_triple_lambda_min": float(tri[:, 0].min()),
        "max_triple_lambda_min": float(tri[:, 0].max()),
        "min_triple_lambda_max": float(tri[:, -1].min()),
        "max_triple_lambda_max": float(tri[:, -1].max()),
        "active_triple_constraints": active,
        "x": [float(v) for v in x],
    }


def stiefel_vector_from_projector(p: np.ndarray) -> np.ndarray:
    w, u = np.linalg.eigh(p)
    a = u[:, np.argsort(w)[-3:]]
    return a.reshape(-1)


def initial_vectors(seed: int, random_starts: int, include_outputs: bool) -> list[tuple[str, np.ndarray]]:
    rng = np.random.default_rng(seed)
    starts: list[tuple[str, np.ndarray]] = []
    for i in range(random_starts):
        starts.append((f"random:{i}", retract(rng.standard_normal((6, 3))).reshape(-1)))
    starts.append(("conference", conference_start()))
    for label, tree, weights in TTSP_CASES:
        starts.append((f"ttsp:{label}", stiefel_vector_from_projector(ttsp_projector(tree, weights))))
    seventh = Path("verify/data/P514_seventh.npy")
    if seventh.exists():
        starts.append(("seventh", stiefel_vector_from_projector(np.load(seventh))))

    if include_outputs:
        for path in [
            "verify/out/v38_relaxed_1over8_slsqp_polish_from_v37.json",
            "verify/out/v38_relaxed_1over8_slsqp_polish_from_bestcore.json",
            "verify/out/v38_relaxed_1over8_slsqp_polish_top2.json",
            "verify/out/v38_relaxed_1over8_slsqp_polish_top3.json",
        ]:
            p = Path(path)
            if not p.exists():
                continue
            data = json.load(open(p))
            cand = data["candidate"]
            starts.append((p.stem, retract(np.array(cand["v"]).reshape(6, 3)).reshape(-1)))
    return starts


def polish_start(v0: np.ndarray, t: float, label: str, maxiter: int, ftol: float) -> dict:
    p0 = v0.reshape(6, 3) @ v0.reshape(6, 3).T
    s0 = max(0.0, float(triple_violations(p0, t).max(initial=0.0)))
    x0 = np.r_[v0, s0]
    result = minimize(
        lambda x: x[-1],
        x0,
        method="SLSQP",
        constraints=[
            {"type": "eq", "fun": stiefel_equalities},
            {"type": "ineq", "fun": lambda x: inequalities(x, t)},
        ],
        bounds=[(None, None)] * 18 + [(-0.1, 0.2)],
        options={"maxiter": maxiter, "ftol": ftol, "disp": False},
    )
    row = summarize(result.x, t, label)
    row.update({
        "success": bool(result.success),
        "message": str(result.message),
        "nit": int(result.nit),
        "objective": float(result.fun),
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
    parser.add_argument("--targets", default="1/8,1/7,3/20,0.16")
    parser.add_argument("--random-starts", type=int, default=24)
    parser.add_argument("--seed", type=int, default=20260816)
    parser.add_argument("--maxiter", type=int, default=2000)
    parser.add_argument("--ftol", type=float, default=1e-12)
    parser.add_argument("--no-output-starts", action="store_true")
    parser.add_argument("--out", default="verify/out/v40_relaxed_target_sweep.json")
    args = parser.parse_args()

    targets = [parse_target(item) for item in args.targets.split(",") if item]
    starts = initial_vectors(args.seed, args.random_starts, not args.no_output_starts)

    print("=" * 78)
    print("REDUCED RELAXED-TARGET SLSQP SWEEP")
    print(f"targets={targets}")
    print(f"starts per target={len(starts)}")
    print("=" * 78, flush=True)

    all_rows = []
    for t in targets:
        rows = []
        print(f"\ntarget t={t:.12f}", flush=True)
        for idx, (label, v0) in enumerate(starts, 1):
            row = polish_start(v0, t, label, args.maxiter, args.ftol)
            rows.append(row)
            best = min(rows, key=lambda item: (item["s"], item["min_inequality"] < -1e-7))
            print(
                f"  [{idx:3d}/{len(starts)}] {label:<35} "
                f"s={row['s']:+.9e} minineq={row['min_inequality']:+.1e} "
                f"best={best['s']:+.9e} ({best['label']})",
                flush=True,
            )
        feasible = [row for row in rows if row["success"] and row["min_inequality"] >= -1e-7]
        ranked = sorted(feasible or rows, key=lambda item: item["s"])
        best = ranked[0]
        print(
            f"  BEST t={t:.12f}: s={best['s']:+.12e}, "
            f"F={best['F']:.12f}, active={len(best['active_triple_constraints'])}, "
            f"label={best['label']}",
            flush=True,
        )
        all_rows.append({
            "target": t,
            "n_starts": len(starts),
            "n_feasible": len(feasible),
            "best": best,
            "top": ranked[:20],
        })
        write_json(args.out, {
            "complete": False,
            "parameters": vars(args),
            "targets": all_rows,
        })

    write_json(args.out, {
        "complete": True,
        "parameters": vars(args),
        "targets": all_rows,
    })
    print(f"\nwrote {args.out}")
    print("=" * 78)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
