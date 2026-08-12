#!/usr/bin/env python3
"""
v29_witness_cover.py -- shrink the §3i.5 disjunction by SET COVER.

THE PROBLEM v28 LEFT.  On the active-equality locus of a low-active orbit A, some
outside triple always overshoots 1/6, but WHICH one varies (12+ witnesses on the
|A|=6 orbit; at the worst points only one overshoots).  Naively the exact certificate
is a disjunction over all 20-|A| outside triples, i.e. ~124*14 ~ 1700 regions.

THE QUESTION HERE.  Do we actually need all 14?  If a SMALL subset S of outside
triples has the property

    at EVERY point of the locus, some T in S overshoots,

then the certificate is a disjunction over |S| branches, not 14.  That is a classic
set-cover problem over the sampled locus points, and it is cheap to compute.

WHAT IS COMPUTED.
  1. sample the active-equality locus densely (independent seeds from v28);
  2. build the incidence matrix  point x outside-triple  of "overshoots by >= tau";
  3. exact minimum cover by ILP-free branch and bound (the instance is tiny), plus
     a greedy upper bound for sanity;
  4. report the cover, its margin (the minimum overshoot achieved by the covering
     triple, over all points), and how the cover behaves as tau rises -- a cover that
     survives a LARGER tau is a stronger certificate.

  A cover of size 1 would resurrect the single-inequality target of §3i.3 (v28 says
  it will not).  A cover of size 2-4 makes the exact layer perhaps 10x smaller than
  the 1700-region estimate and is the realistic good outcome.

HONESTY.  The cover is computed against SAMPLED points, so it is a lower bound on
what the true cover must be: a cover valid on the samples may fail elsewhere on the
locus.  It is therefore a *hypothesis generator* for the exact layer, not a proof.
The exact layer must then verify, for the proposed S, that
    {active equalities} AND {no T in S overshoots}  is empty.
"""
import argparse, itertools, json, os, sys
from collections import Counter
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


def sample_locus(active, n_starts, seed, maxiter=8000):
    rng = np.random.default_rng(seed)
    act = set(active)
    pts = []
    for _ in range(n_starts):
        v = rng.standard_normal(18)
        for _ in range(4):
            v = minimize(eq_resid, v, args=(act,), method="Nelder-Mead",
                         options=dict(maxiter=maxiter, maxfev=maxiter,
                                      xatol=1e-14, fatol=1e-17)).x
        if eq_resid(v, act) < 1e-12:
            A = retract(v.reshape(6, 3))
            pts.append(lam_all(A @ A.T))
    return pts


def min_cover(rows, universe, npts):
    """Exact minimum set cover by deepening; instance is tiny (<=14 sets, a few
    hundred points), so plain deepening is fine and exact.

    BUG FIXED: `full` must be the mask of all POINTS, not of the sets.  Using
    len(rows) (= number of outside triples) made `full` far too small, so the
    search compared against the wrong target and returned None -- which printed as
    'EXACT minimum cover size -' in the first run."""
    full = (1 << npts) - 1
    # rows[i] = bitmask of points covered by outside-triple i
    best = None
    for k in range(1, len(universe) + 1):
        found = None
        for combo in itertools.combinations(range(len(universe)), k):
            m = 0
            for c in combo:
                m |= rows[universe[c]]
            if m == full:
                found = combo
                break
        if found:
            best = [universe[c] for c in found]
            break
    return best


def greedy_cover(rows, universe, npts):
    full = (1 << npts) - 1
    have, pick = 0, []
    while have != full:
        b, bi = -1, None
        for j in universe:
            g = bin(rows[j] & ~have).count("1")
            if g > b:
                b, bi = g, j
        if bi is None or b == 0:
            return None
        pick.append(bi)
        have |= rows[bi]
    return pick


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--canon", type=int, default=78593)
    ap.add_argument("--size", type=int, default=6)
    ap.add_argument("--starts", type=int, default=200)
    ap.add_argument("--seed", type=int, default=20260806)
    a = ap.parse_args()

    d = json.load(open("verify/out/v22_low_active.json"))
    orb = next(x for x in d["full_pair_cover"]
               if x["canon"] == a.canon and x["size"] == a.size)
    active = sorted(IDX[tuple(t)] for t in orb["triples"])
    inact = [i for i in range(20) if i not in active]

    print("=" * 78)
    print(f"WITNESS SET COVER   orbit size {a.size} canon {a.canon}")
    print(f"  active indices: {active}")
    print(f"  outside triples: {len(inact)}")
    print(f"  sampling locus with {a.starts} starts", flush=True)
    print("=" * 78)

    pts = sample_locus(active, a.starts, a.seed)
    print(f"locus points obtained: {len(pts)}/{a.starts}", flush=True)
    if len(pts) < 10:
        print("too few points to cover meaningfully")
        sys.exit(1)

    out = dict(canon=a.canon, size=a.size, active=active,
               n_points=len(pts), starts=a.starts, seed=a.seed, covers={})

    for tau in (1e-9, 1e-3, 1e-2, 2e-2, 3e-2, 4e-2):
        rows = {j: 0 for j in inact}
        uncovered = 0
        for pi, lam in enumerate(pts):
            any_j = False
            for j in inact:
                if lam[j] - TH >= tau:
                    rows[j] |= 1 << pi
                    any_j = True
            if not any_j:
                uncovered += 1
        if uncovered:
            print(f"\ntau = {tau:g}: {uncovered} point(s) covered by NO outside "
                  f"triple -- no cover exists at this margin")
            out["covers"][str(tau)] = dict(exists=False, uncovered=uncovered)
            continue
        g = greedy_cover(rows, inact, len(pts))
        mc = min_cover(rows, inact, len(pts)) if len(inact) <= 16 else None
        cov = mc if mc else g
        # margin of the cover: min over points of the max overshoot within the cover
        margins = []
        for lam in pts:
            margins.append(max(lam[j] - TH for j in cov))
        print(f"\ntau = {tau:g}:")
        print(f"   greedy cover size {len(g) if g else '-'}"
              f"   EXACT minimum cover size {len(mc) if mc else '-'}")
        if cov:
            print(f"   cover = {[list(TRIPLES[j]) for j in cov]}")
            print(f"   cover margin: min={min(margins):+.4e} "
                  f"median={np.median(margins):+.4e}")
        out["covers"][str(tau)] = dict(exists=True,
                                       greedy=len(g) if g else None,
                                       exact=len(mc) if mc else None,
                                       cover=[list(TRIPLES[j]) for j in cov],
                                       cover_indices=[int(j) for j in cov],
                                       margin_min=float(min(margins)),
                                       margin_median=float(np.median(margins)))

    # how many triples overshoot per point, and the singleton points
    singles = Counter()
    for lam in pts:
        over = [j for j in inact if lam[j] - TH > 1e-9]
        if len(over) == 1:
            singles[over[0]] += 1
    print(f"\npoints where EXACTLY ONE outside triple overshoots: {sum(singles.values())}"
          f" of {len(pts)}")
    if singles:
        print("   those forced witnesses (must be in every cover):")
        for j, c in singles.most_common():
            print(f"     {TRIPLES[j]} (idx {j}): {c} point(s)")
        print(f"   => any valid cover contains at least these {len(singles)} triples")
    out["forced_witnesses"] = {str(TRIPLES[j]): c for j, c in singles.items()}

    os.makedirs("verify/out", exist_ok=True)
    json.dump(out, open(f"verify/out/v29_cover_{a.size}_{a.canon}.json", "w"), indent=1)
    print(f"\nwrote verify/out/v29_cover_{a.size}_{a.canon}.json")
    print("=" * 78)
