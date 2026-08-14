#!/usr/bin/env python3
"""Numerical KKT diagnostic for relaxed reduced-obstruction minimizers.

For a target ``t``, ``v40_relaxed_target_sweep.py`` minimizes ``s`` subject to
the closed one-row/two-row deletion core

    1 - 5t <= P_ii <= 5t,
    lambda_min(P_ij) <= 4t,
    lambda_max(P_ij) >= 1 - 4t,

and the triple constraints

    lambda_min(P_T) - t <= s,
    1 - t - lambda_max(P_T) <= s.

This script takes a stored SLSQP candidate and checks the chart-level KKT
balance

    sum beta_j grad(phi_j) - sum mu_k grad(core_k) = 0,
    sum beta_j = 1,
    beta_j, mu_k >= 0,

where the active objective functions are
``phi_j = lambda_min(P_T)-t`` and ``phi_j = 1-t-lambda_max(P_T)``.  Gradients
are central finite differences in a well-conditioned Grassmann chart.  The
result is diagnostic, not a proof; its purpose is to identify which active
orbit is worth exact algebraic certification.
"""

from __future__ import annotations

import argparse
import itertools
import json
import os
from pathlib import Path

import numpy as np
from scipy.optimize import lsq_linear


TRIPLES = list(itertools.combinations(range(6), 3))
PAIRS = list(itertools.combinations(range(6), 2))


def projector_from_stiefel_vector(x: list[float]) -> np.ndarray:
    a = np.array(x[:18], dtype=float).reshape(6, 3)
    return a @ a.T


def stiefel_from_projector(p: np.ndarray) -> np.ndarray:
    values, vectors = np.linalg.eigh(p)
    return vectors[:, np.argsort(values)[-3:]]


def choose_chart(a: np.ndarray) -> tuple[tuple[int, int, int], float]:
    best_rows: tuple[int, int, int] | None = None
    best_det = -1.0
    for rows in TRIPLES:
        det = abs(float(np.linalg.det(a[list(rows), :])))
        if det > best_det:
            best_rows = rows
            best_det = det
    assert best_rows is not None
    return best_rows, best_det


def z_from_projector(p: np.ndarray) -> tuple[np.ndarray, tuple[int, int, int], float]:
    a = stiefel_from_projector(p)
    pivot, det = choose_chart(a)
    y = a @ np.linalg.inv(a[list(pivot), :])
    free = [i for i in range(6) if i not in pivot]
    return y[free, :].reshape(9), pivot, det


def projector_from_z(z: np.ndarray, pivot: tuple[int, int, int]) -> np.ndarray:
    free = [i for i in range(6) if i not in pivot]
    y = np.zeros((6, 3))
    for col, row in enumerate(pivot):
        y[row, col] = 1.0
    y[free, :] = z.reshape(3, 3)
    gram = y.T @ y
    return y @ np.linalg.inv(gram) @ y.T


def eigvals_block(p: np.ndarray, indices: tuple[int, ...]) -> np.ndarray:
    return np.linalg.eigvalsh(p[np.ix_(indices, indices)])


def objective_value(p: np.ndarray, target: float, kind: str, triple_index: int) -> float:
    ev = eigvals_block(p, TRIPLES[triple_index])
    if kind == "low":
        return float(ev[0] - target)
    if kind == "high":
        return float(1.0 - target - ev[-1])
    raise ValueError(f"unknown objective kind {kind!r}")


def core_constraints(p: np.ndarray, target: float) -> list[dict]:
    lev_low = 1.0 - 5.0 * target
    lev_high = 5.0 * target
    pair_low_cut = 4.0 * target
    pair_high_cut = 1.0 - 4.0 * target
    rows: list[dict] = []

    for i, value in enumerate(np.diag(p)):
        rows.append({"kind": "lev_low", "index": [i], "value": float(value - lev_low)})
        rows.append({"kind": "lev_high", "index": [i], "value": float(lev_high - value)})

    for idx, pair in enumerate(PAIRS):
        ev = eigvals_block(p, pair)
        rows.append({"kind": "pair_low", "pair_index": idx, "index": list(pair), "value": float(pair_low_cut - ev[0])})
        rows.append({"kind": "pair_high", "pair_index": idx, "index": list(pair), "value": float(ev[-1] - pair_high_cut)})
    return rows


def core_value(p: np.ndarray, target: float, row: dict) -> float:
    kind = row["kind"]
    if kind == "lev_low":
        i = row["index"][0]
        return float(p[i, i] - (1.0 - 5.0 * target))
    if kind == "lev_high":
        i = row["index"][0]
        return float(5.0 * target - p[i, i])
    pair = tuple(row["index"])
    ev = eigvals_block(p, pair)
    if kind == "pair_low":
        return float(4.0 * target - ev[0])
    if kind == "pair_high":
        return float(ev[-1] - (1.0 - 4.0 * target))
    raise ValueError(f"unknown core kind {kind!r}")


def central_gradient(fun, z: np.ndarray, h: float) -> np.ndarray:
    grad = np.zeros_like(z)
    for j in range(len(z)):
        dz = np.zeros_like(z)
        dz[j] = h
        grad[j] = (fun(z + dz) - fun(z - dz)) / (2.0 * h)
    return grad


def parse_row(data: dict, target_index: int, source: str) -> tuple[dict, dict]:
    target = data["targets"][target_index]
    if source == "best":
        return target, target["best"]
    if source.startswith("top:"):
        idx = int(source.split(":", 1)[1])
        return target, target["top"][idx]
    raise ValueError("source must be 'best' or 'top:<index>'")


def kkt_diagnostic(row: dict, target: float, h: float, active_tol: float, core_tol: float) -> dict:
    p = projector_from_stiefel_vector(row["x"])
    z, pivot, chart_det = z_from_projector(p)
    pz = projector_from_z(z, pivot)
    chart_error = float(np.linalg.norm(p - pz))
    s = float(row["s"])

    objective_rows = []
    for triple_index, triple in enumerate(TRIPLES):
        for kind in ("low", "high"):
            value = objective_value(pz, target, kind, triple_index)
            gap = s - value
            if abs(gap) <= active_tol:
                ev = eigvals_block(pz, triple)
                objective_rows.append({
                    "kind": kind,
                    "triple_index": triple_index,
                    "triple": list(triple),
                    "value": value,
                    "gap_to_s": float(gap),
                    "spectrum": [float(x) for x in ev],
                    "edge_gap": float(ev[1] - ev[0] if kind == "low" else ev[-1] - ev[-2]),
                })

    core_rows = core_constraints(pz, target)
    active_core = [item for item in core_rows if item["value"] <= core_tol]
    min_core = min(item["value"] for item in core_rows)

    obj_grads = []
    for item in objective_rows:
        obj_grads.append(central_gradient(
            lambda zz, it=item: objective_value(
                projector_from_z(zz, pivot),
                target,
                it["kind"],
                it["triple_index"],
            ),
            z,
            h,
        ))

    core_grads = []
    for item in active_core:
        core_grads.append(central_gradient(
            lambda zz, it=item: core_value(projector_from_z(zz, pivot), target, it),
            z,
            h,
        ))

    if not obj_grads:
        raise RuntimeError("no active objective rows found; increase --active-tol")

    obj_mat = np.column_stack(obj_grads)
    if core_grads:
        core_mat = np.column_stack(core_grads)
        top = np.hstack([obj_mat, -core_mat])
        bottom = np.hstack([np.ones((1, obj_mat.shape[1])), np.zeros((1, core_mat.shape[1]))])
    else:
        top = obj_mat
        bottom = np.ones((1, obj_mat.shape[1]))
    mat = np.vstack([top, bottom])
    rhs = np.zeros(10)
    rhs[-1] = 1.0

    sol = lsq_linear(mat, rhs, bounds=(0.0, np.inf), tol=1e-12, lsmr_tol="auto", max_iter=10000)
    coeffs = sol.x
    beta = coeffs[:len(objective_rows)]
    mu = coeffs[len(objective_rows):]
    residual = mat @ coeffs - rhs

    return {
        "target": target,
        "source_label": row.get("label"),
        "s": s,
        "bound_value": float(target + s),
        "F": float(row.get("F", np.nan)),
        "pivot_rows": list(pivot),
        "chart_det_abs": chart_det,
        "chart_reconstruction_error": chart_error,
        "finite_difference_h": h,
        "active_tol": active_tol,
        "core_tol": core_tol,
        "min_core_slack": float(min_core),
        "active_objective_count": len(objective_rows),
        "active_core_count": len(active_core),
        "active_objectives": objective_rows,
        "active_core": active_core,
        "kkt_residual_norm": float(np.linalg.norm(residual)),
        "kkt_residual_inf": float(np.max(np.abs(residual))),
        "lsq_status": int(sol.status),
        "lsq_success": bool(sol.success),
        "lsq_message": str(sol.message),
        "beta_sum": float(beta.sum()),
        "beta_min": float(beta.min(initial=np.inf)),
        "beta_max": float(beta.max(initial=-np.inf)),
        "mu_sum": float(mu.sum()) if len(mu) else 0.0,
        "mu_min": float(mu.min(initial=np.inf)) if len(mu) else None,
        "mu_max": float(mu.max(initial=-np.inf)) if len(mu) else None,
        "beta": [float(x) for x in beta],
        "mu": [float(x) for x in mu],
        "singular_values": [float(x) for x in np.linalg.svd(mat, compute_uv=False)],
    }


def write_json(path: str, payload: dict) -> None:
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    tmp = path + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(payload, fh, indent=1)
    os.replace(tmp, path)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", default="verify/out/v40_relaxed_target_sweep_16over125.json")
    parser.add_argument("--target-index", type=int, default=0)
    parser.add_argument("--source", default="best", help="'best' or 'top:<index>'")
    parser.add_argument("--h", type=float, default=1e-6)
    parser.add_argument("--active-tol", type=float, default=2e-6)
    parser.add_argument("--core-tol", type=float, default=2e-6)
    parser.add_argument("--out", default="verify/out/v43_relaxed_kkt_diagnostic.json")
    args = parser.parse_args()

    data = json.load(open(args.input))
    target_row, row = parse_row(data, args.target_index, args.source)
    target = float(target_row["target"])
    payload = {
        "input": args.input,
        "target_index": args.target_index,
        "source": args.source,
        "diagnostic": kkt_diagnostic(row, target, args.h, args.active_tol, args.core_tol),
    }
    write_json(args.out, payload)

    diag = payload["diagnostic"]
    print("=" * 78)
    print("RELAXED MINIMAX KKT DIAGNOSTIC")
    print(f"input: {args.input}")
    print(f"target={diag['target']:.15f} s={diag['s']:.15f} target+s={diag['bound_value']:.15f}")
    print(f"pivot={diag['pivot_rows']} |det|={diag['chart_det_abs']:.6e}")
    print(f"active objectives={diag['active_objective_count']} active core={diag['active_core_count']}")
    print(f"min core slack={diag['min_core_slack']:+.3e}")
    print(f"KKT residual norm={diag['kkt_residual_norm']:.3e} inf={diag['kkt_residual_inf']:.3e}")
    print(f"beta sum={diag['beta_sum']:.15f}, beta min={diag['beta_min']:+.3e}")
    print(f"mu sum={diag['mu_sum']:.15f}, mu min={diag['mu_min']}")
    print("active core constraints:")
    for item, value in zip(diag["active_core"], diag["mu"]):
        print(f"  {item['kind']:9s} {item['index']} slack={item['value']:+.3e} mu={value:+.3e}")
    print(f"wrote {args.out}")
    print("=" * 78)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
