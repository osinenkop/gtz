#!/usr/bin/env python3
"""
v24_low_active_feasibility.py -- decide, numerically but sharply, whether any of
the 124 surviving low-active orbits from v22 can host a REAL FEASIBLE point.

WHY THIS IS THE RIGHT NEXT STEP.
  For the unique |A|=6 orbit, the saturated determinant ideal has top dimension 3
  over F_32003 (dim 0 with 3 zero-sum sections, empty with 4).  So the EQUALITY part
  of that stratum is a 3-fold: nonempty.  Emptiness of the stratum therefore has to
  come from the INEQUALITY layer, exactly as code/sage/RESULTS.md predicted for the
  known-base system.  The inequalities are:

      active   T in A : lambda_min(P_TT) = 1/6   AND  P_TT - I/6 PSD
      inactive T not in A : lambda_min(P_TT) < 1/6

  This script attacks that directly: minimize a penalty that is zero exactly on the
  feasible set of the stratum, from many starts, over the 9-dim chart.  A converged
  zero would be a genuine low-active extremal -- which (being |A| <= 9 < 10) is
  automatically NON-sharp and would break the simple-active finiteness route.  It is
  the single most consequential object this project could find, so the test is run
  with the tripwire armed and with the result reported either way.

  Failure to find a zero is NOT a proof of emptiness.  It is evidence, and it tells
  the exact-algebra layer which orbits deserve the expensive CAD/lex treatment.

CRITICAL DESIGN POINT -- do not let the optimizer cheat.
  A naive penalty "sum over active T of (lambda_min - 1/6)^2" is minimized by
  configurations where the ACTIVE SET IS NOT the requested one (e.g. the optimizer
  drifts to a known 10-active extremal that happens to contain a 6-subset).  Every
  hit is therefore re-validated: we recompute the ACTUAL active set from scratch and
  keep the point only if it equals the requested orbit (up to S_6 canonicalization).
  This is the same trap that produced the false "8-active" hit in the v18 screens.

Deterministic.  Output is NUMERICALLY SUPPORTED, never PROVED.
"""
import argparse, itertools, json, os, sys
import numpy as np
from scipy.optimize import minimize

TRIPLES = list(itertools.combinations(range(6), 3))
IDX = {t: i for i, t in enumerate(TRIPLES)}
TH = 1.0 / 6.0
PERMS = list(itertools.permutations(range(6)))


def perm_act(p):
    return tuple(IDX[tuple(sorted((p[a], p[b], p[c])))] for (a, b, c) in TRIPLES)


PERM_ACTION = [perm_act(p) for p in PERMS]


def canon_mask(mask):
    best = None
    for act in PERM_ACTION:
        m = 0
        for i in range(20):
            if mask >> i & 1:
                m |= 1 << act[i]
        if best is None or m < best:
            best = m
    return best


def retract(X):
    U, _, Vt = np.linalg.svd(X, full_matrices=False)
    return U @ Vt


def blocks(P):
    return [P[np.ix_(T, T)] for T in TRIPLES]


def lam_min(P):
    return np.array([np.linalg.eigvalsh(B)[0] for B in blocks(P)])


def penalty(v, active, w_in=1.0):
    """Zero exactly on the feasible set of the stratum with active set `active`.

    active   T : lambda_min(P_TT) = 1/6           -> squared residual
    inactive T : lambda_min(P_TT) <= 1/6 - margin -> hinge
    """
    A = retract(v.reshape(6, 3))
    P = A @ A.T
    lam = lam_min(P)
    s = 0.0
    for i in range(20):
        if i in active:
            s += (lam[i] - TH) ** 2
        else:
            # want strictly below; penalize only the violation
            s += w_in * max(0.0, lam[i] - TH) ** 2
    return s


def actual_active(P, tol_mode="gap"):
    """Recompute the actual active set, cutting at the largest multiplicative gap
    in the deviations (never a fixed threshold -- that was the v12 artifact)."""
    lam = lam_min(P)
    dev = np.abs(lam - TH)
    order = np.argsort(dev)
    sd = dev[order]
    best_k, best_r = 1, 0.0
    for k in range(1, 20):
        lo = max(sd[k - 1], 1e-16)
        r = sd[k] / lo
        if r > best_r:
            best_k, best_r = k, r
    return sorted(int(i) for i in order[:best_k]), float(best_r), float(sd[best_k - 1])


def try_orbit(indices, n_starts, seed):
    rng = np.random.default_rng(seed)
    act = set(indices)
    best = (np.inf, None)
    hits = []
    for _ in range(n_starts):
        v = rng.standard_normal(18)
        for _ in range(4):
            v = minimize(penalty, v, args=(act,), method="Nelder-Mead",
                         options=dict(maxiter=8000, maxfev=8000,
                                      xatol=1e-13, fatol=1e-16)).x
        p = penalty(v, act)
        if p < best[0]:
            best = (p, v.copy())
        if p < 1e-14:
            A = retract(v.reshape(6, 3))
            P = A @ A.T
            aa, gap, cut = actual_active(P)
            hits.append(dict(penalty=float(p), actual_active=aa,
                             actual_size=len(aa), gap_ratio=gap,
                             matches_request=bool(sorted(act) == aa),
                             F=float(np.max(lam_min(P))),
                             lev=sorted(np.diag(P).tolist())))
    return best, hits


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--sizes", default="6,7")
    ap.add_argument("--starts", type=int, default=60)
    ap.add_argument("--seed", type=int, default=20260803)
    ap.add_argument("--limit", type=int, default=0)
    a = ap.parse_args()
    sizes = [int(x) for x in a.sizes.split(",")]

    d = json.load(open("verify/out/v22_low_active.json"))
    orbits = [x for x in d["full_pair_cover"] if x["size"] in sizes]
    orbits.sort(key=lambda z: (z["size"], z["canon"]))
    if a.limit:
        orbits = orbits[:a.limit]

    print("=" * 78)
    print(f"LOW-ACTIVE FEASIBILITY PROBE   sizes={sizes}  orbits={len(orbits)}")
    print(f"  starts per orbit: {a.starts}   seed: {a.seed}")
    print("  A converged zero with the REQUESTED active set would be a genuine")
    print("  low-active extremal => automatically non-sharp => breaks the")
    print("  simple-active finiteness route.  Reported either way.")
    print("=" * 78, flush=True)

    results, genuine = [], []
    for n, orb in enumerate(orbits, 1):
        idx = sorted(IDX[tuple(t)] for t in orb["triples"])
        (bp, bv), hits = try_orbit(idx, a.starts, a.seed + 7919 * n)
        matched = [h for h in hits if h["matches_request"]]
        drift = [h for h in hits if not h["matches_request"]]
        rec = dict(size=orb["size"], canon=orb["canon"], indices=idx,
                   best_penalty=float(bp), n_zero=len(hits),
                   n_matching=len(matched), n_drifted=len(drift),
                   drifted_sizes=sorted({h["actual_size"] for h in drift}))
        results.append(rec)
        flag = ""
        if matched:
            flag = "   *** GENUINE LOW-ACTIVE POINT ***"
            genuine.append(dict(rec, sample=matched[0]))
        print(f"[{n}/{len(orbits)}] size {orb['size']} canon {orb['canon']}: "
              f"best penalty {bp:.3e}  zeros {len(hits)} "
              f"(matching {len(matched)}, drifted {len(drift)}"
              f"{' -> sizes ' + str(rec['drifted_sizes']) if drift else ''})"
              f"{flag}", flush=True)

    print("\n" + "=" * 78)
    print(f"orbits probed:                        {len(results)}")
    print(f"orbits with a genuine matching point: {len(genuine)}")
    if genuine:
        print("\n*** THE SIMPLE-ACTIVE FINITENESS ROUTE IS IN DOUBT. ***")
        print("*** Certify these exactly before drawing any conclusion.  ***")
        for g in genuine[:5]:
            print(f"    size {g['size']} canon {g['canon']}  "
                  f"actual={g['sample']['actual_active']}  F={g['sample']['F']!r}")
    else:
        print("\nNo orbit admitted a feasible point with its REQUESTED active set.")
        print("Every numerical zero drifted to a larger actual active set, i.e. the")
        print("optimizer could only satisfy the equalities by landing on a")
        print("higher-active (known) extremal.  This is evidence -- NOT proof --")
        print("that the low-active strata are empty.")
    os.makedirs("verify/out", exist_ok=True)
    json.dump(dict(sizes=sizes, starts=a.starts, seed=a.seed,
                   results=results, genuine=genuine),
              open("verify/out/v24_low_active_feasibility.json", "w"), indent=1)
    print("\nwrote verify/out/v24_low_active_feasibility.json")
    print("=" * 78)
    sys.exit(2 if genuine else 0)
