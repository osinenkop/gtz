#!/usr/bin/env python3
"""
v18_kkt_screen.py -- symmetry-reduced numerical KKT screen for GTZ(6,3).

This is a candidate-finding and triage script, not a proof.  It searches active
sets modulo S_6 for points P in Gr(3,6) where:

  * the selected triples share the same value t = lambda_min(P_TT);
  * all unselected triples have lambda_min <= t, so the selected set can be the
    active set of F(P);
  * a convex combination of selected active gradients has zero tangent projection
    on Gr(3,6), the Clarke/KKT stationarity condition for a local minimizer of F.

Two modes are supported:

  counterexample  Search for stationary candidates with t < 1/6.
  equality        Search for stationary candidates on the equality level t = 1/6.

Any successful counterexample candidate must still be reconstructed and checked
independently in exact/certified arithmetic before it means anything.  The
counterexample tripwire is armed: if a polished point with F < 1/6 - 1e-6 is
found, the script writes it separately and exits with code 2.
"""
from __future__ import annotations

import argparse
import itertools
import json
import os
import sys
from multiprocessing import Pool

import numpy as np
from scipy.optimize import least_squares

from v17_active_orbits import S6_ACTIONS, TRIPLES, apply_action, ttsp_projector

TH = 1 / 6
CX = 1e-6
E = ("e",)

TTSP_CASES = [
    ("P(S(e,e,e),e,e,e)", ("P", ("S", E, E, E), E, E, E), [1, 1, 1, 5 / 9, 5 / 9, 5 / 9]),
    ("P(S(e,e),S(e,e),e,e)", ("P", ("S", E, E), ("S", E, E), E, E), [1, 1, 1, 1, 5 / 8, 5 / 8]),
    ("P(S(P(e,e,e),e),S(e,e))", ("P", ("S", ("P", E, E, E), E), ("S", E, E)),
     [1, 1, 1, 9 / 5, 9 / 5, 9 / 5]),
    ("P(S(P(S(e,e),e,e),e),e)", ("P", ("S", ("P", ("S", E, E), E, E), E), E),
     [1, 1, 5 / 8, 5 / 8, 1, 1]),
    ("P(S(P(e,e),P(e,e)),S(e,e))", ("P", ("S", ("P", E, E), ("P", E, E)), ("S", E, E)),
     [1, 1, 1, 1, 8 / 5, 8 / 5]),
    ("P(S(P(e,e),e),S(P(e,e),e))", ("P", ("S", ("P", E, E), E), ("S", ("P", E, E), E)),
     [1, 1, 8 / 5, 1, 1, 8 / 5]),
]


def retract(x: np.ndarray) -> np.ndarray:
    u, _, vt = np.linalg.svd(x.reshape(6, 3), full_matrices=False)
    return u @ vt


def projector_from_y(y: np.ndarray, m: int) -> tuple[np.ndarray, np.ndarray]:
    a = retract(y[:18])
    logits = y[18:18 + m]
    logits = logits - np.max(logits)
    alpha = np.exp(logits)
    alpha /= np.sum(alpha)
    return a @ a.T, alpha


def triple_data(p: np.ndarray) -> tuple[np.ndarray, list[np.ndarray]]:
    lams = np.zeros(len(TRIPLES))
    vecs = []
    for i, triple in enumerate(TRIPLES):
        w, u = np.linalg.eigh(p[np.ix_(triple, triple)])
        lams[i] = w[0]
        vecs.append(u[:, 0])
    return lams, vecs


def embed_gradient(active_index: int, vec: np.ndarray) -> np.ndarray:
    h = np.zeros((6, 6))
    triple = TRIPLES[active_index]
    block = np.outer(vec, vec)
    for a, ia in enumerate(triple):
        for b, ib in enumerate(triple):
            h[ia, ib] = block[a, b]
    return h


def tangent_residual(p: np.ndarray, h: np.ndarray) -> np.ndarray:
    q = np.eye(6) - p
    # For symmetric H, QHP = 0 is equivalent to zero Grassmann tangent gradient.
    return (q @ h @ p).reshape(-1)


def residual(
    y: np.ndarray,
    active: tuple[int, ...],
    mode: str,
    margin: float,
    weights: dict[str, float],
) -> np.ndarray:
    p, alpha = projector_from_y(y, len(active))
    lams, vecs = triple_data(p)
    active_lams = lams[list(active)]
    t = float(np.mean(active_lams))
    inactive = [i for i in range(len(TRIPLES)) if i not in active]

    parts = []
    parts.append(weights["spread"] * (active_lams - t))
    if mode == "equality":
        parts.append(np.array([weights["level"] * (t - TH)]))
    else:
        # This is a one-sided penalty: a true sub-1/6 KKT point has zero here.
        parts.append(np.array([weights["below"] * max(t - (TH - margin), 0.0)]))
    if inactive:
        parts.append(weights["inactive"] * np.maximum(lams[inactive] - t, 0.0))

    h = np.zeros((6, 6))
    for coeff, idx in zip(alpha, active):
        h += coeff * embed_gradient(idx, vecs[idx])
    parts.append(weights["stationarity"] * tangent_residual(p, h))

    # Any genuine counterexample can be assumed core by the existing Case-A
    # reduction.  In equality mode this is only a useful stabilizer.
    parts.append(weights["core"] * np.maximum(TH - np.diag(p), 0.0))
    return np.concatenate(parts)


def metrics(y: np.ndarray, active: tuple[int, ...], mode: str, margin: float) -> dict:
    p, alpha = projector_from_y(y, len(active))
    lams, vecs = triple_data(p)
    active_lams = lams[list(active)]
    t = float(np.mean(active_lams))
    inactive = [i for i in range(len(TRIPLES)) if i not in active]
    h = np.zeros((6, 6))
    for coeff, idx in zip(alpha, active):
        h += coeff * embed_gradient(idx, vecs[idx])
    stationarity = float(np.linalg.norm(tangent_residual(p, h)))
    spread = float(np.max(np.abs(active_lams - t))) if len(active_lams) else 0.0
    inactive_excess = float(np.max(lams[inactive] - t)) if inactive else -float("inf")
    core_deficit = float(np.max(TH - np.diag(p)))
    fval = float(np.max(lams))
    actual_active = [i for i, v in enumerate(lams) if abs(v - fval) < 1e-7]
    selected = set(active)
    actual = set(actual_active)
    loss = float(np.linalg.norm(residual(y, active, mode, margin, DEFAULT_WEIGHTS)))
    return {
        "loss": loss,
        "t": t,
        "F": fval,
        "spread": spread,
        "inactive_excess": inactive_excess,
        "stationarity": stationarity,
        "core_deficit": core_deficit,
        "min_leverage": float(np.min(np.diag(p))),
        "max_leverage": float(np.max(np.diag(p))),
        "alpha_min": float(np.min(alpha)),
        "alpha_max": float(np.max(alpha)),
        "actual_active_size": len(actual_active),
        "actual_active": actual_active,
        "exact_active": bool(actual == selected),
        "selected_inside_actual": bool(selected.issubset(actual)),
        "tripwire": bool(fval < TH - CX),
        "mode_ok": bool((mode == "equality" and abs(t - TH) < 1e-7)
                        or (mode == "counterexample" and t < TH - margin / 2)),
    }


DEFAULT_WEIGHTS = {
    "spread": 50.0,
    "level": 50.0,
    "below": 50.0,
    "inactive": 25.0,
    "stationarity": 10.0,
    "core": 10.0,
}


def mask_from_indices(indices: tuple[int, ...]) -> int:
    mask = 0
    for i in indices:
        mask |= 1 << i
    return mask


def indices_from_mask(mask: int) -> tuple[int, ...]:
    return tuple(i for i in range(len(TRIPLES)) if mask & (1 << i))


def orbit_representatives(size_min: int, size_max: int, max_orbits: int | None = None) -> list[int]:
    reps = []
    for size in range(size_min, size_max + 1):
        seen = set()
        for combo in itertools.combinations(range(len(TRIPLES)), size):
            mask = mask_from_indices(combo)
            if mask in seen:
                continue
            orbit = {apply_action(mask, action) for action in S6_ACTIONS}
            seen.update(orbit)
            reps.append(min(orbit))
            if max_orbits is not None and len(reps) >= max_orbits:
                return reps
    return reps


def initial_y(rng: np.random.Generator, m: int, p0: np.ndarray | None = None) -> np.ndarray:
    if p0 is None:
        a0 = rng.standard_normal((6, 3))
    else:
        w, u = np.linalg.eigh(p0)
        a0 = u[:, np.argsort(-w)[:3]]
        a0 += 0.01 * rng.standard_normal((6, 3))
    logits = 0.1 * rng.standard_normal(m)
    return np.concatenate([a0.reshape(-1), logits])


def fit_one(
    active: tuple[int, ...],
    mode: str,
    margin: float,
    starts: int,
    seed: int,
    max_nfev: int,
    p0: np.ndarray | None = None,
) -> tuple[dict, list[float]]:
    rng = np.random.default_rng(seed)
    best_y = None
    best_loss = float("inf")
    losses = []
    for _ in range(starts):
        y0 = initial_y(rng, len(active), p0)
        res = least_squares(
            residual,
            y0,
            args=(active, mode, margin, DEFAULT_WEIGHTS),
            max_nfev=max_nfev,
            xtol=1e-10,
            ftol=1e-10,
            gtol=1e-10,
        )
        loss = float(np.linalg.norm(res.fun))
        losses.append(loss)
        if loss < best_loss:
            best_loss = loss
            best_y = res.x
    assert best_y is not None
    rec = metrics(best_y, active, mode, margin)
    rec["best_y"] = best_y.tolist()
    rec["losses"] = losses
    return rec, losses


def worker(args: tuple[int, str, float, int, int, int]) -> dict:
    mask, mode, margin, starts, seed, max_nfev = args
    active = indices_from_mask(mask)
    try:
        rec, _ = fit_one(active, mode, margin, starts, seed, max_nfev)
        rec["error"] = None
    except Exception as exc:
        rec = {
            "loss": float("inf"),
            "error": repr(exc),
            "tripwire": False,
            "mode_ok": False,
            "mask": mask,
            "active_size": len(active),
            "active": list(active),
            "active_triples": [list(TRIPLES[i]) for i in active],
        }
        return rec
    rec["mask"] = mask
    rec["active_size"] = len(active)
    rec["active"] = list(active)
    rec["active_triples"] = [list(TRIPLES[i]) for i in active]
    return rec


def known_smoke(mode: str, margin: float, starts: int, max_nfev: int) -> list[dict]:
    out = []
    for label, tree, weights in TTSP_CASES:
        p0 = ttsp_projector(tree, weights)
        lams, _ = triple_data(p0)
        active = tuple(i for i, v in enumerate(lams) if abs(v - TH) < 1e-8)
        rec, _ = fit_one(active, mode, margin, starts, 12345 + len(out), max_nfev, p0=p0)
        rec.pop("best_y", None)
        rec["label"] = label
        rec["active_size"] = len(active)
        out.append(rec)

    seventh = "verify/data/P514_seventh.npy"
    if os.path.exists(seventh):
        p0 = np.load(seventh)
        lams, _ = triple_data(p0)
        active = tuple(i for i, v in enumerate(lams) if abs(v - TH) < 1e-8)
        rec, _ = fit_one(active, mode, margin, starts, 54321, max_nfev, p0=p0)
        rec.pop("best_y", None)
        rec["label"] = "OUT-OF-FAMILY (5/14,9/14)"
        rec["active_size"] = len(active)
        out.append(rec)
    return out


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", choices=["counterexample", "equality"], default="counterexample")
    parser.add_argument("--size-min", type=int, default=1)
    parser.add_argument("--size-max", type=int, default=4)
    parser.add_argument("--starts", type=int, default=3)
    parser.add_argument("--max-orbits", type=int, default=None)
    parser.add_argument("--max-nfev", type=int, default=2000)
    parser.add_argument("--seed", type=int, default=20260803)
    parser.add_argument("--jobs", type=int, default=max(1, (os.cpu_count() or 4) - 2))
    parser.add_argument("--margin", type=float, default=1e-4)
    parser.add_argument("--known-smoke", action="store_true")
    args = parser.parse_args()

    os.makedirs("verify/out", exist_ok=True)
    print("=" * 78)
    print("NUMERICAL ACTIVE-SET KKT SCREEN")
    print(f"mode={args.mode} sizes={args.size_min}..{args.size_max} starts={args.starts}")
    print(f"max_nfev={args.max_nfev} jobs={args.jobs} seed={args.seed}")
    print("=" * 78, flush=True)

    smoke = []
    if args.known_smoke:
        print("\nknown-extremal smoke test:")
        smoke = known_smoke("equality", args.margin, max(1, args.starts), args.max_nfev)
        for rec in smoke:
            print(f"  {rec['label']:<42} |A|={rec['active_size']:2d} "
                  f"loss={rec['loss']:.3e} stationarity={rec['stationarity']:.3e} "
                  f"t-1/6={rec['t'] - TH:+.3e}")

    reps = orbit_representatives(args.size_min, args.size_max, args.max_orbits)
    print(f"\norbit representatives to screen: {len(reps)}", flush=True)
    tasks = [(mask, args.mode, args.margin, args.starts, args.seed + 1009 * i, args.max_nfev)
             for i, mask in enumerate(reps)]

    records = []
    if tasks:
        with Pool(args.jobs) as pool:
            for i, rec in enumerate(pool.imap_unordered(worker, tasks, chunksize=1), 1):
                records.append(rec)
                if rec["tripwire"]:
                    path = "verify/out/v18_KKT_TRIPWIRE.json"
                    with open(path, "w") as fh:
                        json.dump(rec, fh, indent=1)
                    print(f"\n*** TRIPWIRE: F={rec['F']:.16f} at mask {rec['mask']} ***")
                    print(f"wrote {path}")
                    pool.terminate()
                    return 2
                if i % 25 == 0 or i == len(tasks):
                    best = min(records, key=lambda r: r["loss"])
                    print(f"  [{i}/{len(tasks)}] best loss={best['loss']:.3e} "
                          f"t={best.get('t', float('nan')):.9f} "
                          f"F={best.get('F', float('nan')):.9f} "
                          f"|A|={best['active_size']} "
                          f"actual={best.get('actual_active_size')}",
                          flush=True)
                    partial_path = f"verify/out/v18_kkt_{args.mode}_{args.size_min}_{args.size_max}.partial.json"
                    with open(partial_path, "w") as fh:
                        json.dump({"done": i, "total": len(tasks), "records": records}, fh, indent=1)

    records.sort(key=lambda r: r["loss"])
    candidates = [
        r for r in records
        if r["loss"] < 1e-5 and r["inactive_excess"] < 1e-7 and r["stationarity"] < 1e-6
    ]
    if args.mode == "counterexample":
        candidates = [r for r in candidates if r["t"] < TH - args.margin / 2]

    print("\nsummary:")
    print(f"  screened:   {len(records)} orbit reps")
    errors = [r for r in records if r.get("error")]
    print(f"  worker errors: {len(errors)}")
    print(f"  candidates: {len(candidates)}")
    if records:
        print("  best records:")
        for rec in records[:10]:
            print(f"    loss={rec['loss']:.3e} |A|={rec['active_size']:2d} "
                  f"t={rec['t']:.9f} F={rec['F']:.9f} "
                  f"spread={rec['spread']:.1e} inactive={rec['inactive_excess']:.1e} "
                  f"stat={rec['stationarity']:.1e} actual={rec['actual_active_size']} "
                  f"exact={rec['exact_active']}")

    out = {
        "mode": args.mode,
        "size_min": args.size_min,
        "size_max": args.size_max,
        "starts": args.starts,
        "max_nfev": args.max_nfev,
        "seed": args.seed,
        "margin": args.margin,
        "n_orbits": len(reps),
        "n_candidates": len(candidates),
        "known_smoke": smoke,
        "records": records,
        "candidates": candidates,
    }
    path = f"verify/out/v18_kkt_{args.mode}_{args.size_min}_{args.size_max}.json"
    with open(path, "w") as fh:
        json.dump(out, fh, indent=1)
    print(f"\nwrote {path}")
    print("=" * 78)
    return 0


if __name__ == "__main__":
    sys.exit(main())
