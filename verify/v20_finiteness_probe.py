#!/usr/bin/env python3
"""
v20_finiteness_probe.py -- classify equality-KKT hits for the finiteness route.

The obstruction to finiteness would be a non-sharp equality extremal.  The KKT
screen v18 searches selected active-set strata, but a selected set can be a
proper subset of the actual active set at the polished point.  This script
therefore reloads v18 outputs, reconstructs the projector from best_y, computes
the actual active set, canonicalizes it under S_6, compares it with the known
certified extremal orbits, and runs a numerical sharpness test on the actual
active gradients.

This is a diagnostic, not a proof.  Any genuinely new or non-sharp candidate
reported here must be reconstructed and certified independently.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from collections import Counter, defaultdict

import numpy as np
from scipy.optimize import linprog

sys.path.append(os.path.dirname(__file__))
from v17_active_orbits import TRIPLES, canonical_mask, known_extremal_orbits
from v18_kkt_screen import indices_from_mask, mask_from_indices, projector_from_y, triple_data

TH = 1.0 / 6.0


DEFAULT_INPUTS = [
    "verify/out/v18_kkt_equality_1_4.json",
    "verify/out/v18_kkt_equality_5_6.json",
    "verify/out/v18_kkt_equality_7_8.json",
    "verify/out/v18_kkt_equality_9_9.json",
    "verify/out/v18_kkt_equality_10_13.json",
    "verify/out/v18_kkt_equality_14_20.json",
]


def tangent_basis(projector: np.ndarray) -> list[np.ndarray]:
    eigvals, eigvecs = np.linalg.eigh(projector)
    order = np.argsort(-eigvals)
    range_basis = eigvecs[:, order[:3]]
    kernel_basis = eigvecs[:, order[3:]]
    basis = []
    for a in range(3):
        for b in range(3):
            x = np.outer(range_basis[:, a], kernel_basis[:, b])
            x = x + x.T
            norm = np.linalg.norm(x, "fro")
            basis.append(x / norm)
    return basis


def gradient_matrix(projector: np.ndarray, active: list[int]) -> tuple[np.ndarray, list[float]]:
    basis = tangent_basis(projector)
    rows = []
    gaps = []
    for idx in active:
        triple = TRIPLES[idx]
        values, vecs = np.linalg.eigh(projector[np.ix_(triple, triple)])
        gaps.append(float(values[1] - values[0]) if len(values) > 1 else 0.0)
        vec = vecs[:, 0]
        rows.append([
            float(vec @ tangent[np.ix_(triple, triple)] @ vec)
            for tangent in basis
        ])
    return np.array(rows), gaps


def positive_multiplier_margin(lmat: np.ndarray) -> dict:
    """Maximize min_i lambda_i subject to L^T lambda=0, sum lambda=1."""
    m = lmat.shape[0]
    if m == 0:
        return {"success": False, "margin": None, "message": "empty active set"}
    # variables: lambda_0..lambda_{m-1}, s
    c = np.zeros(m + 1)
    c[-1] = -1.0
    a_eq = np.zeros((lmat.shape[1] + 1, m + 1))
    a_eq[:lmat.shape[1], :m] = lmat.T
    a_eq[-1, :m] = 1.0
    b_eq = np.zeros(lmat.shape[1] + 1)
    b_eq[-1] = 1.0
    a_ub = np.zeros((m, m + 1))
    for i in range(m):
        a_ub[i, i] = -1.0
        a_ub[i, -1] = 1.0
    b_ub = np.zeros(m)
    bounds = [(0.0, 1.0)] * m + [(0.0, 1.0)]
    result = linprog(c, A_ub=a_ub, b_ub=b_ub, A_eq=a_eq, b_eq=b_eq, bounds=bounds, method="highs")
    if not result.success:
        return {"success": False, "margin": None, "message": result.message}
    return {
        "success": True,
        "margin": float(result.x[-1]),
        "message": result.message,
        "residual": float(np.linalg.norm(lmat.T @ result.x[:m])),
    }


def analyze_projector(projector: np.ndarray, active_tol: float) -> dict:
    lams, _vecs = triple_data(projector)
    fval = float(np.max(lams))
    actual = [i for i, value in enumerate(lams) if abs(float(value) - fval) <= active_tol]
    actual_mask = mask_from_indices(tuple(actual))
    actual_canonical = canonical_mask(actual_mask)
    lmat, gaps = gradient_matrix(projector, actual)
    sv = np.linalg.svd(lmat, compute_uv=False) if lmat.size else np.array([])
    rank = int(np.sum(sv > (sv[0] * 1e-8 if len(sv) else 0.0)))
    mult = positive_multiplier_margin(lmat)
    return {
        "F": fval,
        "actual_active_size": len(actual),
        "actual_active": actual,
        "actual_active_triples": [list(TRIPLES[i]) for i in actual],
        "actual_canonical_mask": int(actual_canonical),
        "rank_L": rank,
        "singular_values": [float(x) for x in sv],
        "gap_min": float(min(gaps)) if gaps else None,
        "positive_multiplier": mult,
        "sharp_numeric": bool(rank == 9 and mult["success"] and (mult["margin"] or 0.0) > 1e-8),
    }


def load_records(paths: list[str]) -> list[dict]:
    rows = []
    for path in paths:
        if not os.path.exists(path):
            continue
        with open(path) as fh:
            doc = json.load(fh)
        for rec in doc.get("records", []):
            rec = dict(rec)
            rec["_source"] = path
            rows.append(rec)
    return rows


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("inputs", nargs="*", default=DEFAULT_INPUTS)
    parser.add_argument("--loss-tol", type=float, default=1e-5)
    parser.add_argument("--level-tol", type=float, default=1e-6)
    parser.add_argument("--active-tol", type=float, default=1e-7)
    parser.add_argument("--output", default="verify/out/v20_finiteness_probe.json")
    args = parser.parse_args()

    known_by_mask = defaultdict(list)
    for rec in known_extremal_orbits():
        known_by_mask[int(rec["canonical_mask"])].append(rec["label"])

    records = load_records(args.inputs)
    hits = []
    for rec in records:
        if "best_y" not in rec:
            continue
        if rec.get("loss", float("inf")) > args.loss_tol:
            continue
        if abs(rec.get("F", 999.0) - TH) > args.level_tol:
            continue
        active = tuple(rec.get("active", indices_from_mask(rec["mask"])))
        projector, _alpha = projector_from_y(np.array(rec["best_y"], dtype=float), len(active))
        analysis = analyze_projector(projector, args.active_tol)
        selected_mask = mask_from_indices(active)
        hit = {
            "source": rec["_source"],
            "selected_active_size": len(active),
            "selected_canonical_mask": int(canonical_mask(selected_mask)),
            "selected_inside_actual": rec.get("selected_inside_actual"),
            "loss": rec.get("loss"),
            "t": rec.get("t"),
            "F_record": rec.get("F"),
            **analysis,
        }
        hit["known_orbit_labels"] = known_by_mask.get(hit["actual_canonical_mask"], [])
        hit["new_active_orbit"] = not bool(hit["known_orbit_labels"])
        hit["nonsharp_candidate"] = bool(
            abs(hit["F"] - TH) <= args.level_tol
            and (not hit["sharp_numeric"])
        )
        hits.append(hit)

    counts = Counter((h["actual_active_size"], h["actual_canonical_mask"]) for h in hits)
    new_hits = [h for h in hits if h["new_active_orbit"]]
    nonsharp = [h for h in hits if h["nonsharp_candidate"]]
    small_actual = [h for h in hits if h["actual_active_size"] <= 9]

    print("=" * 78)
    print("FINITENESS PROBE: equality KKT hits classified by actual active set")
    print("=" * 78)
    print(f"input records:       {len(records)}")
    print(f"low-res equality hits: {len(hits)}")
    print(f"distinct actual active orbits: {len(counts)}")
    print(f"new active orbits:   {len(new_hits)}")
    print(f"nonsharp candidates: {len(nonsharp)}")
    print(f"actual |A| <= 9:     {len(small_actual)}")
    print()
    for i, hit in enumerate(sorted(hits, key=lambda h: (h["new_active_orbit"], h["actual_active_size"], h["loss"]))):
        labels = "; ".join(hit["known_orbit_labels"]) if hit["known_orbit_labels"] else "NEW"
        margin = hit["positive_multiplier"]["margin"]
        print(
            f"[{i:02d}] actual |A|={hit['actual_active_size']:2d} "
            f"rank={hit['rank_L']} sharp={hit['sharp_numeric']} "
            f"margin={margin if margin is not None else float('nan'):.3e} "
            f"new={hit['new_active_orbit']} labels={labels}"
        )
        print(f"     source={hit['source']} selected |A|={hit['selected_active_size']} loss={hit['loss']:.3e}")

    out = {
        "inputs": args.inputs,
        "loss_tol": args.loss_tol,
        "level_tol": args.level_tol,
        "active_tol": args.active_tol,
        "n_records": len(records),
        "n_hits": len(hits),
        "distinct_actual_active_orbits": len(counts),
        "n_new_active_orbits": len(new_hits),
        "n_nonsharp_candidates": len(nonsharp),
        "n_actual_active_le_9": len(small_actual),
        "hits": hits,
    }
    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    with open(args.output, "w") as fh:
        json.dump(out, fh, indent=1)
    print(f"\nwrote {args.output}")
    print("=" * 78)
    return 1 if new_hits or nonsharp or small_actual else 0


if __name__ == "__main__":
    raise SystemExit(main())
