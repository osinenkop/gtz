#!/usr/bin/env python3
"""
v27_obstruction_split.py -- decompose WHY the low-active strata are empty, by
separating the two possible obstructions.

THE QUESTION THIS ANSWERS.  v24 showed no low-active orbit admits a feasible point
(penalty floors 9e-4 .. 2.8e-2, never 1e-14).  v26 then showed something sharper and
initially surprising for the |A|=6 orbit: ON the active-equality locus (reached to
residual 9e-11), the active PSD conditions are NOT violated -- every active 2x2
principal minor is positive by at least 1.2e-2.  So the obstruction is NOT active
PSD.  By elimination it must be the INACTIVE conditions.

That is a much better statement than "the stratum is empty", because it names a
mechanism:

    on the active-equality locus of A, some triple OUTSIDE A also reaches 1/6,
    i.e. THE ACTIVE SET CANNOT BE EXACTLY A -- IT ALWAYS GROWS.

If that holds for every A with |A| <= 9, then no extremal has |A| <= 9, which is
exactly the finiteness statement we want (combined with sharp => |A| >= 10).

WHAT IS MEASURED, per orbit A:
  step 1 -- can we land on the active-equality locus at all?
            minimise   R(P) = sum_{T in A} (lambda_min(P_TT) - 1/6)^2
            (equalities ONLY; no inequality terms, so the optimizer is not fighting
            two objectives at once -- this was the flaw that made v26 report "no
            point reached the locus" for the size-7 orbits).
  step 2 -- on the points reached, report both candidate obstructions:
            (a) active PSD slack:  min over T in A of the smallest 2x2 principal
                minor of (P_TT - I/6).  Negative => active PSD fails.
            (b) inactive excess:   max over T not in A of (lambda_min(P_TT) - 1/6).
                Positive => an outside triple also hits the threshold => the actual
                active set is strictly larger than A.
  step 3 -- classify the orbit by which obstruction is operative, and report the
            actual active set size that the optimizer lands on.

A uniform verdict "(b) is positive, with margin" across all 124 orbits is the
certificate target for the exact layer, and it is a statement about ONE polynomial
inequality per orbit -- no root isolation, no degree-1880 elimination.

Numerical.  NUMERICALLY SUPPORTED, never PROVED.
"""
import argparse, itertools, json, os, sys
import numpy as np
from scipy.optimize import minimize

TRIPLES = list(itertools.combinations(range(6), 3))
IDX = {t: i for i, t in enumerate(TRIPLES)}
TH = 1.0 / 6.0


def retract(X):
    U, _, Vt = np.linalg.svd(X, full_matrices=False)
    return U @ Vt


def lam_all(P):
    return np.array([np.linalg.eigvalsh(P[np.ix_(T, T)])[0] for T in TRIPLES])


def eq_resid(v, active):
    A = retract(v.reshape(6, 3))
    P = A @ A.T
    lam = lam_all(P)
    return float(sum((lam[i] - TH) ** 2 for i in active))


def diagnose(v, active):
    A = retract(v.reshape(6, 3))
    P = A @ A.T
    lam = lam_all(P)
    # (a) active PSD slack: smallest 2x2 principal minor of (P_TT - I/6)
    worst_minor = np.inf
    for i in active:
        M = P[np.ix_(TRIPLES[i], TRIPLES[i])] - TH * np.eye(3)
        for idx in itertools.combinations(range(3), 2):
            worst_minor = min(worst_minor, float(np.linalg.det(M[np.ix_(idx, idx)])))
    # (b) inactive excess
    inact = [i for i in range(20) if i not in active]
    excess = max(float(lam[i] - TH) for i in inact) if inact else -np.inf
    # actual active set by the gap rule
    dev = np.abs(lam - TH)
    order = np.argsort(dev)
    sd = dev[order]
    bk, br = 1, 0.0
    for k in range(1, 20):
        r = sd[k] / max(sd[k - 1], 1e-16)
        if r > br:
            bk, br = k, r
    return dict(worst_active_minor=worst_minor, inactive_excess=excess,
                actual_active=sorted(int(i) for i in order[:bk]),
                actual_size=int(bk), gap_ratio=float(br),
                F=float(np.max(lam)))


def probe(indices, n_starts, seed):
    rng = np.random.default_rng(seed)
    act = set(indices)
    onloc = []
    best_r = np.inf
    for _ in range(n_starts):
        v = rng.standard_normal(18)
        for _ in range(4):
            v = minimize(eq_resid, v, args=(act,), method="Nelder-Mead",
                         options=dict(maxiter=8000, maxfev=8000,
                                      xatol=1e-14, fatol=1e-17)).x
        r = eq_resid(v, act)
        best_r = min(best_r, r)
        if r < 1e-12:
            onloc.append((r, diagnose(v, act)))
    return best_r, onloc


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--sizes", default="6,7")
    ap.add_argument("--starts", type=int, default=40)
    ap.add_argument("--seed", type=int, default=20260804)
    ap.add_argument("--limit", type=int, default=0)
    a = ap.parse_args()
    sizes = [int(x) for x in a.sizes.split(",")]

    d = json.load(open("verify/out/v22_low_active.json"))
    orbits = [x for x in d["full_pair_cover"] if x["size"] in sizes]
    orbits.sort(key=lambda z: (z["size"], z["canon"]))
    if a.limit:
        orbits = orbits[:a.limit]

    print("=" * 78)
    print("OBSTRUCTION SPLIT: why are the low-active strata empty?")
    print(f"  sizes={sizes}  orbits={len(orbits)}  starts={a.starts}")
    print("  (a) active PSD slack  < 0  => active PSD is the obstruction")
    print("  (b) inactive excess   > 0  => an OUTSIDE triple also reaches 1/6,")
    print("                                so the active set always GROWS")
    print("=" * 78, flush=True)

    rows, kind = [], {"psd": 0, "inactive": 0, "both": 0, "neither": 0, "off": 0}
    for n, orb in enumerate(orbits, 1):
        idx = sorted(IDX[tuple(t)] for t in orb["triples"])
        br, onloc = probe(idx, a.starts, a.seed + 99991 * n)
        if not onloc:
            kind["off"] += 1
            rows.append(dict(size=orb["size"], canon=orb["canon"], on_locus=0,
                             best_eq_residual=float(br)))
            print(f"[{n}/{len(orbits)}] size {orb['size']} canon {orb['canon']}: "
                  f"never reached the locus (best eq residual {br:.2e})", flush=True)
            continue
        # take the point with the LEAST obstruction: most favourable to feasibility
        best = min(onloc, key=lambda z: max(-z[1]["worst_active_minor"],
                                            z[1]["inactive_excess"]))
        dg = best[1]
        psd_bad = dg["worst_active_minor"] < -1e-9
        inact_bad = dg["inactive_excess"] > 1e-9
        k = ("both" if psd_bad and inact_bad else
             "psd" if psd_bad else "inactive" if inact_bad else "neither")
        kind[k] += 1
        rows.append(dict(size=orb["size"], canon=orb["canon"], on_locus=len(onloc),
                         best_eq_residual=float(br), kind=k, **dg))
        print(f"[{n}/{len(orbits)}] size {orb['size']} canon {orb['canon']}: "
              f"on-locus {len(onloc)}/{a.starts}  "
              f"PSDslack {dg['worst_active_minor']:+.3e}  "
              f"inact.excess {dg['inactive_excess']:+.3e}  "
              f"actual|A| {dg['actual_size']}  => {k}", flush=True)

    print("\n" + "=" * 78)
    print("classification of the operative obstruction:")
    for k, v in kind.items():
        print(f"   {k:<9}: {v}")
    got = [r for r in rows if r.get("on_locus")]
    if got:
        ex = np.array([r["inactive_excess"] for r in got])
        ps = np.array([r["worst_active_minor"] for r in got])
        print(f"\n inactive excess over orbits reaching the locus:")
        print(f"   min={ex.min():+.4e}  median={np.median(ex):+.4e}  max={ex.max():+.4e}")
        print(f" active PSD slack:")
        print(f"   min={ps.min():+.4e}  median={np.median(ps):+.4e}  max={ps.max():+.4e}")
        sizes_hit = sorted({r["actual_size"] for r in got})
        print(f"\n actual active sizes landed on: {sizes_hit}")
        if all(r["inactive_excess"] > 1e-9 for r in got):
            print("\n *** UNIFORM MECHANISM: on every reachable active-equality locus,")
            print("     some OUTSIDE triple also reaches 1/6.  The active set cannot")
            print("     stay at A -- it always grows.  This, not active PSD, is the")
            print("     obstruction, and it is the right certificate target. ***")
    os.makedirs("verify/out", exist_ok=True)
    json.dump(dict(sizes=sizes, starts=a.starts, seed=a.seed,
                   classification=kind, rows=rows),
              open("verify/out/v27_obstruction_split.json", "w"), indent=1)
    print("\nwrote verify/out/v27_obstruction_split.json")
    print("=" * 78)
