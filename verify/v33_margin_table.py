#!/usr/bin/env python3
"""Summarise v30 margin-minimisation outputs.

For each saved v30 JSON, reconstruct the best projector and report the numerical
margin, the outside triples attaining the max, and the triples numerically pinned
at 1/6.  This is a bookkeeping tool for the low-active semialgebraic route.
"""
import argparse
import glob
import itertools
import json
import os

import numpy as np

TRIPLES = list(itertools.combinations(range(6), 3))
TH = 1.0 / 6.0


def retract(X):
    U, _, Vt = np.linalg.svd(X, full_matrices=False)
    return U @ Vt


def lam_all(P):
    return np.array([np.linalg.eigvalsh(P[np.ix_(T, T)])[0] for T in TRIPLES])


def tlabel(i):
    return "".join(str(x) for x in TRIPLES[int(i)])


def summarise(path, tol):
    data = json.load(open(path))
    if "best_point" not in data:
        return None
    active = set(data["active"])
    v = np.array(data["best_point"], dtype=float)
    A = retract(v.reshape(6, 3))
    P = A @ A.T
    lam = lam_all(P)
    maxv = float(lam.max())
    max_tie = [i for i, x in enumerate(lam) if abs(float(x) - maxv) <= tol]
    near_th = [i for i, x in enumerate(lam) if abs(float(x) - TH) <= tol]
    return dict(
        file=path,
        size=data["size"],
        canon=data["canon"],
        starts=data.get("completed_starts") or data.get("starts"),
        complete=data.get("complete"),
        n_on_locus=data["n_on_locus"],
        verdict=data["verdict"],
        g_min=data["g_min"],
        F=maxv,
        active=sorted(active),
        max_tie=max_tie,
        max_tie_out=[i for i in max_tie if i not in active],
        near_th=near_th,
        active_resid=float(sum((lam[i] - TH) ** 2 for i in active)),
    )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("inputs", nargs="*", default=[])
    ap.add_argument("--glob", default="verify/out/v30_margin_*.json")
    ap.add_argument("--tol", type=float, default=1e-6)
    ap.add_argument("--json", action="store_true")
    args = ap.parse_args()

    paths = args.inputs or sorted(glob.glob(args.glob))
    rows = [summarise(path, args.tol) for path in paths]
    rows = [row for row in rows if row is not None]
    rows.sort(key=lambda r: (r["size"], r["canon"]))

    if args.json:
        json.dump(rows, fp=os.sys.stdout, indent=1)
        print()
        return

    print("| orbit | starts | locus hits | g_min | F | max-tie outside | at 1/6 | verdict |")
    print("|---|---:|---:|---:|---:|---|---|---|")
    for r in rows:
        outside = ",".join(tlabel(i) for i in r["max_tie_out"]) or "-"
        near = ",".join(tlabel(i) for i in r["near_th"]) or "-"
        starts = str(r["starts"])
        if r["complete"] is False:
            starts += "*"
        print(f"| ({r['size']},{r['canon']}) | {starts} | {r['n_on_locus']} | "
              f"{r['g_min']:.6e} | {r['F']:.12f} | {outside} | {near} | "
              f"{r['verdict']} |")


if __name__ == "__main__":
    main()
