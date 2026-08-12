#!/usr/bin/env python3
"""
v31_trace_uniformity.py -- the decisive test for dissolving the §3i.5 disjunction.

THE REFORMULATION.  For a triple T put M_T = P_TT - (1/6)I, eigenvalues mu_i.
  det M_T = mu1 mu2 mu3, so
     exactly one mu <= 0  =>  det <= 0
     exactly two  mu <= 0  =>  det >= 0     <-- the ONLY case breaking the encoding
     all three   mu <= 0  =>  det <= 0
  Two eigenvalues of P_TT at or below 1/6 force the third to satisfy
     tr(P_TT) <= 1/6 + 1/6 + 1 = 4/3      (using P <= I, so every eigenvalue <= 1).
  Hence:

      tr(P_TT) > 4/3   ==>   [ lambda_min(P_TT) <= 1/6  <==>  det(P_TT - I/6) <= 0 ]

  ONE polynomial inequality per inactive triple, NO 3-way Descartes disjunction.
  (This is the off-slice form of Lemma 2 of slice-framework.md, which the corpus
  states but which we had not exploited for the inactive constraints.)

WHY THIS COULD DISSOLVE ~1700 REGIONS INTO 124 SYSTEMS.  The disjunction of §3i.5
came from encoding "NOT PSD" via Descartes sign patterns, which branches. If instead
every triple on the active-equality locus is HEAVY (tr(P_TT) > 4/3), the encoding
above is exact everywhere and each orbit becomes a single system:

      T in A     : det(6 N_TT - d I) = 0          (degree 6)
      T not in A : det(6 N_TT - d I) <= 0         (degree 6)
      d > 0,  plus tr conditions

with no case split at all.  If instead LIGHT triples can occur on the locus, the
branching returns -- over subsets of light triples, which would be worse (2^14), so
this test decides between a big simplification and a real obstruction.

WHAT IS MEASURED.  Minimise  min_T tr(P_TT)  over the active-equality locus of an
orbit, from many starts.  The threshold is 4/3.
  * min stays > 4/3  =>  heavy everywhere; the encoding is exact; DISJUNCTION GONE.
  * min dips <= 4/3  =>  light triples occur; branching returns.

Context: sum_T tr(P_TT) = 10 tr(P) = 30 identically, so the MEAN trace is 1.5 and
4/3 = 1.333 is only 11% below the mean -- the test is genuinely tight, not a formality.
Random projectors DO have light triples (observed min 0.13), so this is a real
question about the locus, not about Gr(3,6).
"""
import argparse, itertools, json, os, sys
import numpy as np
from scipy.optimize import minimize

TRIPLES = list(itertools.combinations(range(6), 3))
IDX = {t: i for i, t in enumerate(TRIPLES)}
TH = 1.0 / 6.0
FOUR3 = 4.0 / 3.0


def retract(X):
    U, _, Vt = np.linalg.svd(X, full_matrices=False)
    return U @ Vt


def lam_all(P):
    return np.array([np.linalg.eigvalsh(P[np.ix_(T, T)])[0] for T in TRIPLES])


def traces(P):
    return np.array([float(np.trace(P[np.ix_(T, T)])) for T in TRIPLES])


def eq_resid_P(P, active):
    lam = lam_all(P)
    return float(sum((lam[i] - TH) ** 2 for i in active))


def obj(v, active, w):
    A = retract(v.reshape(6, 3))
    P = A @ A.T
    return float(np.min(traces(P))) + w * eq_resid_P(P, active)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--canon", type=int, default=78593)
    ap.add_argument("--size", type=int, default=6)
    ap.add_argument("--starts", type=int, default=120)
    ap.add_argument("--seed", type=int, default=20260808)
    a = ap.parse_args()

    d = json.load(open("verify/out/v22_low_active.json"))
    orb = next(x for x in d["full_pair_cover"]
               if x["canon"] == a.canon and x["size"] == a.size)
    active = sorted(IDX[tuple(t)] for t in orb["triples"])

    print("=" * 78)
    print(f"TRACE UNIFORMITY TEST  orbit size {a.size} canon {a.canon}")
    print(f"  active {active}")
    print("  minimising  min_T tr(P_TT)  ON the active-equality locus")
    print(f"  threshold 4/3 = {FOUR3:.6f};  identity: sum_T tr(P_TT) = 30, mean 1.5")
    print("  min > 4/3  =>  encoding det(M_T)<=0 is EXACT  =>  disjunction dissolves")
    print("=" * 78, flush=True)

    rng = np.random.default_rng(a.seed)
    rows = []
    for k in range(a.starts):
        v = rng.standard_normal(18)
        for w in (1e3, 1e5, 1e6):
            for _ in range(3):
                v = minimize(obj, v, args=(active, w), method="Nelder-Mead",
                             options=dict(maxiter=8000, maxfev=8000,
                                          xatol=1e-15, fatol=1e-18)).x
        A = retract(v.reshape(6, 3))
        P = A @ A.T
        r = eq_resid_P(P, active)
        tr = traces(P)
        if r < 1e-12:
            j = int(np.argmin(tr))
            rows.append(dict(min_trace=float(tr.min()), resid=r,
                             argmin_triple=list(TRIPLES[j]),
                             heavy=bool(tr.min() > FOUR3),
                             F=float(np.max(lam_all(P)))))
        if (k + 1) % 30 == 0:
            ms = [x["min_trace"] for x in rows]
            print(f"  [{k+1}/{a.starts}] on-locus {len(rows)}  "
                  f"min min_trace so far "
                  f"{min(ms) if ms else float('nan'):.6f}", flush=True)

    if not rows:
        print("no on-locus points")
        sys.exit(1)
    rows.sort(key=lambda z: z["min_trace"])
    mt = np.array([x["min_trace"] for x in rows])
    print(f"\non-locus points: {len(rows)}")
    print(f"  min_T tr(P_TT):  min={mt.min():.6f}  p5={np.percentile(mt,5):.6f}  "
          f"median={np.median(mt):.6f}")
    print(f"  threshold 4/3 =  {FOUR3:.6f}")
    light = [x for x in rows if not x["heavy"]]
    print(f"\npoints with a LIGHT triple (tr <= 4/3): {len(light)} of {len(rows)}")
    for x in rows[:6]:
        print(f"   min_trace={x['min_trace']:.6f} at {x['argmin_triple']}  "
              f"heavy={x['heavy']}  resid={x['resid']:.1e}  F={x['F']:.8f}")

    margin = float(mt.min() - FOUR3)
    verdict = "HEAVY-UNIFORM" if margin > 0 else "LIGHT-OCCURS"
    print(f"\nmargin above 4/3: {margin:+.6f}")
    print(f"VERDICT: {verdict}")
    if verdict == "HEAVY-UNIFORM":
        print("  Every locus point has all 20 traces > 4/3, so on this locus")
        print("    lambda_min(P_TT) <= 1/6  <==>  det(P_TT - I/6) <= 0")
        print("  is an EXACT equivalence.  The inactive constraints become ONE")
        print("  degree-6 polynomial inequality each -- no Descartes disjunction,")
        print("  no 1700 regions.  Per orbit: |A| equalities + (20-|A|)")
        print("  inequalities + d>0, which is exactly what Singular/CAD wants.")
        print(f"  CAVEAT: the margin is only {margin:.4f} -- thin, so the exact")
        print("  layer must certify tr > 4/3 on the locus, not assume it.")
    else:
        print("  Light triples occur on the locus, so the det-encoding is not exact")
        print("  there and branching returns (over subsets of light triples).")

    os.makedirs("verify/out", exist_ok=True)
    json.dump(dict(canon=a.canon, size=a.size, active=active, starts=a.starts,
                   seed=a.seed, n_on_locus=len(rows), verdict=verdict,
                   min_trace=float(mt.min()), margin_above_4_3=margin,
                   threshold=FOUR3, rows=rows[:40]),
              open(f"verify/out/v31_trace_{a.size}_{a.canon}.json", "w"), indent=1)
    print(f"\nwrote verify/out/v31_trace_{a.size}_{a.canon}.json")
    print("=" * 78)
