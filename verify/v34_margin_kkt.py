#!/usr/bin/env python3
"""Numerical KKT diagnostic for v30 low-active margin minimizers.

For a best point of

    minimise max_{T outside A} (lambda_min(P_TT)-1/6)
    subject to lambda_min(P_TT)=1/6 for T in A,

the nonsmooth KKT condition says that the gradients of the outside triples tied
at the maximum should have a convex combination lying in the span of the active
constraint gradients.  This script checks that condition in the standard chart
Y=[I;Z] by finite differences.
"""
import argparse
import itertools
import json
import os

import numpy as np

TRIPLES = list(itertools.combinations(range(6), 3))
TH = 1.0 / 6.0


def retract(X):
    U, _, Vt = np.linalg.svd(X, full_matrices=False)
    return U @ Vt


def projector_from_z(z):
    Y = np.vstack([np.eye(3), z.reshape(3, 3)])
    G = Y.T @ Y
    return Y @ np.linalg.inv(G) @ Y.T


def z_from_point(v):
    A = retract(v.reshape(6, 3))
    top = A[:3, :]
    det = float(np.linalg.det(top))
    if abs(det) < 1e-8:
        raise ValueError(f"top chart is ill-conditioned, det={det}")
    Y = A @ np.linalg.inv(top)
    return Y[3:, :].reshape(9), det


def lam_all(P):
    return np.array([np.linalg.eigvalsh(P[np.ix_(T, T)])[0] for T in TRIPLES])


def lam_index_from_z(z, idx):
    P = projector_from_z(z)
    T = TRIPLES[idx]
    return float(np.linalg.eigvalsh(P[np.ix_(T, T)])[0])


def grad_lam(z, idx, h):
    g = np.zeros_like(z)
    for j in range(len(z)):
        dz = np.zeros_like(z)
        dz[j] = h
        g[j] = (lam_index_from_z(z + dz, idx) - lam_index_from_z(z - dz, idx)) / (2 * h)
    return g


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", default="verify/out/v30_margin_7_14113.json")
    ap.add_argument("--h", type=float, default=1e-6)
    ap.add_argument("--tie-tol", type=float, default=1e-6)
    ap.add_argument("--out", default="")
    args = ap.parse_args()

    data = json.load(open(args.input))
    active = list(data["active"])
    z, chart_det = z_from_point(np.array(data["best_point"], dtype=float))
    P = projector_from_z(z)
    lam = lam_all(P)
    outside = [i for i in range(20) if i not in set(active)]
    F = float(lam.max())
    ties = [i for i in outside if abs(float(lam[i]) - F) <= args.tie_tol]
    if not ties:
        raise SystemExit("no outside max ties found; increase --tie-tol")

    grad_active = np.column_stack([grad_lam(z, i, args.h) for i in active])
    grad_ties = np.column_stack([grad_lam(z, i, args.h) for i in ties])
    top = np.hstack([grad_active, grad_ties])
    bottom = np.hstack([np.zeros((1, len(active))), np.ones((1, len(ties)))])
    M = np.vstack([top, bottom])
    rhs = np.zeros(10)
    rhs[-1] = 1.0
    sol, _, _, sing = np.linalg.lstsq(M, rhs, rcond=None)
    alpha = sol[:len(active)]
    beta = sol[len(active):]
    resid = M @ sol - rhs
    payload = {
        "input": args.input,
        "h": args.h,
        "tie_tol": args.tie_tol,
        "chart_det_top": chart_det,
        "F": F,
        "margin": F - TH,
        "active": active,
        "active_triples": [list(TRIPLES[i]) for i in active],
        "ties": ties,
        "tie_triples": [list(TRIPLES[i]) for i in ties],
        "kkt_residual_norm": float(np.linalg.norm(resid)),
        "kkt_singular_values": [float(x) for x in sing],
        "alpha": [float(x) for x in alpha],
        "beta": [float(x) for x in beta],
        "beta_sum": float(beta.sum()),
        "beta_min": float(beta.min()),
    }

    print("=" * 78)
    print(f"MARGIN KKT DIAGNOSTIC {args.input}")
    print(f"chart det(top 3 rows) = {chart_det:+.6e}")
    print(f"F = {F:.15f}   F-1/6 = {F - TH:.15e}")
    print(f"active: {active}")
    print(f"  {[TRIPLES[i] for i in active]}")
    print(f"ties:   {ties}")
    print(f"  {[TRIPLES[i] for i in ties]}")
    print(f"finite-difference h = {args.h:g}")
    print(f"KKT residual norm = {np.linalg.norm(resid):.3e}")
    print(f"smallest singular value of KKT matrix = {sing[-1]:.3e}")
    print("\nactive multipliers alpha:")
    for idx, val in zip(active, alpha):
        print(f"  {idx:2d} {TRIPLES[idx]}  {val:+.12e}")
    print("\noutside tie convex weights beta:")
    for idx, val in zip(ties, beta):
        print(f"  {idx:2d} {TRIPLES[idx]}  {val:+.12e}")
    print(f"  sum beta = {float(beta.sum()):.15f}")
    print(f"  min beta = {float(beta.min()):+.3e}")
    if args.out:
        os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
        with open(args.out, "w") as f:
            json.dump(payload, f, indent=1)
        print(f"\nwrote {args.out}")
    print("=" * 78)


if __name__ == "__main__":
    main()
