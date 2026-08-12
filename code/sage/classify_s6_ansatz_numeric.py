#!/usr/bin/env python3
"""Numerically classify real roots of the structured size-6 ansatz branch.

This is a fast companion to the exact modular eliminant.  It reads the lifted
degree-11 residual q-polynomial, adds the obvious real nuisance q-values, solves
the branch equations with q fixed, and classifies roots by the actual GTZ
eigenvalue inequalities.
"""
from __future__ import annotations

import argparse
import itertools
import json
import math
from fractions import Fraction
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares


TRIPLES = list(itertools.combinations(range(6), 3))
ACTIVE = (0, 1, 9, 15, 16, 17)
TIES = (5, 7, 12, 14, 18)
TH = 1.0 / 6.0


def parse_fraction(text: str) -> Fraction:
    return Fraction(text)


def residual_q_roots(path: Path, imag_tol: float) -> list[float]:
    data = json.loads(path.read_text())
    raw = data["rational_coefficients_low_to_high"]
    coeffs = [float(parse_fraction(text)) for text in raw]
    roots = np.roots(list(reversed(coeffs)))
    real = sorted(float(root.real) for root in roots if abs(root.imag) <= imag_tol)
    return real


def q_values(path: Path, imag_tol: float) -> list[float]:
    values = [3.0, 5.0]
    values.extend(residual_q_roots(path, imag_tol))
    out = []
    for value in sorted(values):
        if not out or abs(value - out[-1]) > 1e-8:
            out.append(value)
    return out


def projector(a: float, b: float, c: float, d: float) -> np.ndarray:
    y = np.array(
        [
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [-a, a, -b],
            [c, d, -a],
            [-d, -c, -a],
        ],
        dtype=float,
    )
    gram = y.T @ y
    return y @ np.linalg.inv(gram) @ y.T


def determinant_equations(x: np.ndarray, q: float) -> np.ndarray:
    a, b, c, d = x
    p = projector(a, b, c, d)
    values = []
    for idx in ACTIVE:
        block = p[np.ix_(TRIPLES[idx], TRIPLES[idx])] - TH * np.eye(3)
        values.append(np.linalg.det(block))
    for idx in TIES:
        block = p[np.ix_(TRIPLES[idx], TRIPLES[idx])] - (q * TH) * np.eye(3)
        values.append(np.linalg.det(block))
    values.append((c + d) ** 2 - 5.0)
    return np.array(values, dtype=float)


def lambdas(p: np.ndarray) -> np.ndarray:
    return np.array(
        [np.linalg.eigvalsh(p[np.ix_(triple, triple)])[0] for triple in TRIPLES],
        dtype=float,
    )


def z_from_x(x: np.ndarray) -> list[float]:
    a, b, c, d = x
    return [-a, a, -b, c, d, -a, -d, -c, -a]


def root_seen(roots: list[np.ndarray], root: np.ndarray, tol: float) -> bool:
    return any(np.linalg.norm(root - old) <= tol * max(1.0, np.linalg.norm(old)) for old in roots)


def classify_root(x: np.ndarray, q: float) -> dict:
    p = projector(*x)
    lam = lambdas(p)
    active_set = set(ACTIVE)
    inactive = [idx for idx in range(len(TRIPLES)) if idx not in active_set]
    active_shifts = {idx: float(lam[idx] - TH) for idx in ACTIVE}
    inactive_shifts = {idx: float(lam[idx] - TH) for idx in inactive}
    active_min_index = min(active_shifts, key=active_shifts.get)
    inactive_max_index = max(inactive_shifts, key=inactive_shifts.get)
    active_min = active_shifts[active_min_index]
    inactive_max = inactive_shifts[inactive_max_index]
    tie_max_abs = max(abs(float(lam[idx] - q * TH)) for idx in TIES)
    actual_active = [idx for idx, value in enumerate(lam) if abs(float(value - TH)) <= 1e-7]
    return {
        "q": float(q),
        "x": [float(v) for v in x],
        "z": z_from_x(x),
        "residual_inf": float(np.max(np.abs(determinant_equations(x, q)))),
        "F_minus_1_6": float(np.max(lam) - TH),
        "active_min_lambda_minus_1_6": active_min,
        "active_min_index": active_min_index,
        "active_min_triple": list(TRIPLES[active_min_index]),
        "inactive_max_lambda_minus_1_6": inactive_max,
        "inactive_max_index": inactive_max_index,
        "inactive_max_triple": list(TRIPLES[inactive_max_index]),
        "tie_max_abs_lambda_minus_q_over_6": tie_max_abs,
        "lambda_minus_1_6": [float(value - TH) for value in lam],
        "passes_active_psd": bool(active_min >= -1e-7),
        "passes_inactive": bool(inactive_max <= 1e-7),
        "passes_equality_inequalities": bool(active_min >= -1e-7 and inactive_max <= 1e-7),
        "actual_active": actual_active,
        "actual_active_size": len(actual_active),
    }


def starts_for_q(q: float, count: int, rng: np.random.Generator, target: np.ndarray | None) -> list[np.ndarray]:
    starts: list[np.ndarray] = []
    if target is not None and abs(q - target[4]) <= 1e-5:
        starts.append(target[:4].copy())
        starts.append(np.array([-target[0], target[1], target[2], target[3]], dtype=float))
        starts.append(np.array([target[0], -target[1], target[2], target[3]], dtype=float))
    root5 = math.sqrt(5.0)
    for sign in (-1.0, 1.0):
        for t in np.linspace(-3.0, 3.0, 9):
            c = (sign * root5 - t) / 2.0
            d = (sign * root5 + t) / 2.0
            starts.append(np.array([0.6, 1.3, c, d], dtype=float))
            starts.append(np.array([-0.6, -1.3, c, d], dtype=float))
    while len(starts) < count:
        sign = -1.0 if rng.random() < 0.5 else 1.0
        t = rng.uniform(-4.0, 4.0)
        c = (sign * root5 - t) / 2.0
        d = (sign * root5 + t) / 2.0
        starts.append(np.array([rng.uniform(-2.5, 2.5), rng.uniform(-3.0, 3.0), c, d], dtype=float))
    return starts[:count]


def solve_for_q(q: float, starts: list[np.ndarray], max_nfev: int, accept_residual: float, root_tol: float) -> list[dict]:
    roots: list[np.ndarray] = []
    records: list[dict] = []
    for start in starts:
        sol = least_squares(
            lambda x: determinant_equations(x, q),
            start,
            max_nfev=max_nfev,
            xtol=1e-12,
            ftol=1e-12,
            gtol=1e-12,
        )
        residual = float(np.max(np.abs(sol.fun)))
        if residual > accept_residual:
            continue
        if root_seen(roots, sol.x, root_tol):
            continue
        roots.append(sol.x.copy())
        records.append(classify_root(sol.x, q))
    records.sort(key=lambda row: (row["inactive_max_lambda_minus_1_6"], row["active_min_lambda_minus_1_6"]))
    return records


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--q-residual",
        default="code/sage/out/s6_ansatz_cplusd_q_residual_crt_6primes.json",
    )
    parser.add_argument("--starts-per-q", type=int, default=400)
    parser.add_argument("--max-nfev", type=int, default=5000)
    parser.add_argument("--accept-residual", type=float, default=1e-10)
    parser.add_argument("--root-tol", type=float, default=1e-7)
    parser.add_argument("--imag-tol", type=float, default=1e-7)
    parser.add_argument("--seed", type=int, default=9001)
    parser.add_argument("--out", default="code/sage/out/classify_s6_ansatz_cplusd_numeric.json")
    args = parser.parse_args()

    q_roots = q_values(Path(args.q_residual), args.imag_tol)
    target = np.array(
        [
            0.6270504207739223,
            1.3086725496015859,
            0.7524391608884586,
            1.4836288166113012,
            1.1806958367665793,
        ],
        dtype=float,
    )
    rng = np.random.default_rng(args.seed)
    all_records = []
    per_q = []
    for q in q_roots:
        starts = starts_for_q(q, args.starts_per_q, rng, target)
        records = solve_for_q(q, starts, args.max_nfev, args.accept_residual, args.root_tol)
        all_records.extend(records)
        per_q.append(
            {
                "q": float(q),
                "roots": len(records),
                "passes_active_psd": int(sum(row["passes_active_psd"] for row in records)),
                "passes_inactive": int(sum(row["passes_inactive"] for row in records)),
                "passes_both": int(sum(row["passes_equality_inequalities"] for row in records)),
                "best_inactive_max_lambda_minus_1_6": None
                if not records
                else min(row["inactive_max_lambda_minus_1_6"] for row in records),
                "best_active_min_lambda_minus_1_6": None
                if not records
                else max(row["active_min_lambda_minus_1_6"] for row in records),
            }
        )

    payload = {
        "q_residual": args.q_residual,
        "q_roots": q_roots,
        "starts_per_q": args.starts_per_q,
        "accepted_roots": len(all_records),
        "passes_active_psd": int(sum(row["passes_active_psd"] for row in all_records)),
        "passes_inactive": int(sum(row["passes_inactive"] for row in all_records)),
        "passes_both": int(sum(row["passes_equality_inequalities"] for row in all_records)),
        "per_q": per_q,
        "roots": all_records,
    }
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(payload, indent=1) + "\n")

    print("=" * 78)
    print("S6 ANSATZ C+ D NUMERIC CLASSIFICATION")
    print(f"q roots: {len(q_roots)}")
    print(f"accepted roots: {len(all_records)}")
    print(f"passes active PSD: {payload['passes_active_psd']}")
    print(f"passes inactive:   {payload['passes_inactive']}")
    print(f"passes both:       {payload['passes_both']}")
    for row in per_q:
        print(
            f"q={row['q']:.12g} roots={row['roots']} "
            f"active={row['passes_active_psd']} inactive={row['passes_inactive']} both={row['passes_both']}"
        )
    print(f"wrote {out}")
    print("=" * 78)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
