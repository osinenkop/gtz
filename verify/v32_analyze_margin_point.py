#!/usr/bin/env python3
"""Inspect numerical minimizers of the low-active margin problem.

The v30 run stores the best chart vector for the unique size-6 orbit.  This script
reconstructs the projector and prints the spectral/combinatorial pattern: which
outside triples attain the max, which triples are nearest to 1/6, row leverage
symmetries, and low-degree rational-looking values.
"""
import argparse
import itertools
import json
from fractions import Fraction

import numpy as np

TRIPLES = list(itertools.combinations(range(6), 3))
TH = 1.0 / 6.0


def retract(X):
    U, _, Vt = np.linalg.svd(X, full_matrices=False)
    return U @ Vt


def lam_all(P):
    return np.array([np.linalg.eigvalsh(P[np.ix_(T, T)])[0] for T in TRIPLES])


def approx(x, max_den=500):
    return str(Fraction(float(x)).limit_denominator(max_den))


def shifted_invariants(P, triple):
    M = P[np.ix_(triple, triple)] - TH * np.eye(3)
    minors1 = [float(M[i, i]) for i in range(3)]
    minors2 = []
    for ij in itertools.combinations(range(3), 2):
        minors2.append(float(np.linalg.det(M[np.ix_(ij, ij)])))
    return {
        "trace": float(np.trace(M)),
        "min_1x1": min(minors1),
        "min_2x2": min(minors2),
        "det": float(np.linalg.det(M)),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", default="verify/out/v30_margin_6_78593.json")
    ap.add_argument("--point", choices=["best_point"], default="best_point")
    ap.add_argument("--tol", type=float, default=1e-7)
    args = ap.parse_args()

    data = json.load(open(args.input))
    active = set(data["active"])
    v = np.array(data[args.point], dtype=float)
    A = retract(v.reshape(6, 3))
    P = A @ A.T
    lam = lam_all(P)

    print("=" * 78)
    print(f"ANALYSE MARGIN POINT from {args.input}")
    print(f"active indices: {sorted(active)}")
    print("=" * 78)
    print(f"projector residual ||P^2-P||_F = {np.linalg.norm(P @ P - P):.3e}")
    print(f"trace(P) = {np.trace(P):.15f}")
    print(f"F = {lam.max():.15f}   F-1/6 = {lam.max() - TH:.15e}")
    print(f"active equality residual = {sum((lam[i] - TH) ** 2 for i in active):.3e}")

    lev = np.diag(P)
    print("\nrow leverages:")
    for i, x in enumerate(lev):
        print(f"  {i}: {x:.15f}  ~ {approx(x)}")

    print("\nleast eigenvalues by triple:")
    order = np.argsort(lam)
    for i in order:
        tag = "A" if int(i) in active else "out"
        print(f"  {i:2d} {TRIPLES[int(i)]} {tag:3s}  "
              f"lambda_min={lam[i]:.15f}  delta={lam[i]-TH:+.6e}  ~ {approx(lam[i])}")

    maxv = float(lam.max())
    max_triples = [int(i) for i, x in enumerate(lam) if abs(float(x) - maxv) <= args.tol]
    near_th = [int(i) for i, x in enumerate(lam) if abs(float(x) - TH) <= args.tol]
    print(f"\ntriples attaining F within {args.tol:g}: {max_triples}")
    print(f"  {[TRIPLES[i] for i in max_triples]}")
    print(f"triples at 1/6 within {args.tol:g}: {near_th}")
    print(f"  {[TRIPLES[i] for i in near_th]}")

    print("\nshifted invariants for outside triples attaining F:")
    for i in max_triples:
        inv = shifted_invariants(P, TRIPLES[i])
        print(f"  {i:2d} {TRIPLES[i]}  "
              f"tr={inv['trace']:+.6e} min1={inv['min_1x1']:+.6e} "
              f"min2={inv['min_2x2']:+.6e} det={inv['det']:+.6e}")

    print("\nP matrix:")
    for row in P:
        print("  " + " ".join(f"{x:+.12f}" for x in row))
    print("=" * 78)


if __name__ == "__main__":
    main()
