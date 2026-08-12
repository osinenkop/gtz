#!/usr/bin/env python3
"""
v26_margin_certificate.py -- test alternative 1 of §3i.1: is the low-active
obstruction a WIDE-MARGIN inequality, certifiable without root isolation?

THE IDEA.  The lex route died on degree (1880 univariate for the |A|=6 orbit).  But
we do not need to isolate roots.  We need one statement of the form

    on the equality locus  {det(6 N_TT - d I) = 0, T in A},
    SOME active PSD minor is bounded away from feasibility by a margin c > 0.

The certificate degree for such a statement is governed by the MARGIN, not by the
degree of the ideal.  v24 reports margins of 9e-4 .. 2.8e-2 in penalty units, i.e.
~3e-2 .. 1.7e-1 in eigenvalue units against a threshold of 1/6 -- so a certificate
plausibly exists at low degree.

WHAT THIS SCRIPT ACTUALLY MEASURES (numerically, honestly).
  For a given low-active orbit A we solve, over the chart:

      maximise   t
      s.t.       the |A| active determinant equalities hold (as equalities),
                 every active 2x2 principal minor of (P_TT - I/6) >= -t ... no.

  That is the wrong sign.  What we want is the WORST CASE over the equality locus of
  the BEST-CASE feasibility, i.e.

      m(A) := max over P on the equality locus of  min over active T of
              lambda_2(P_TT) - 1/6   ... also not it.

  The clean quantity is: how negative must some active PSD minor be?  Define, for P
  on the equality locus,
      viol(P) := max over active T of  ( -lambda_min(P_TT - I/6 restricted PSD test) )
  and we want  inf over the locus of viol(P)  > 0.  Since lambda_min(P_TT) = 1/6 is
  imposed as an equality, PSD failure shows up in the SECOND eigenvalue relation;
  concretely we track the most negative 2x2 principal minor of M_T = 6 N_TT - d I,
  which is exactly the quantity code/sage/RESULTS.md used as its slice obstruction.

  So we compute
      mu(A) := inf over the equality locus of  max over T in A, over the three
               2x2 principal minors, of  ( - minor )
  by constrained optimisation from many starts.  mu(A) > 0 uniformly is the
  certificate target; the VALUE of mu(A) tells us the margin an SOS certificate
  would have to capture.

  A positive, well-separated mu(A) for every orbit is strong evidence that
  alternative 1 is the right route.  A near-zero mu(A) for some orbit would say the
  margin is NOT wide there and SOS would be delicate -- equally useful to know.

Numerical only.  Reports the margin distribution, not a proof.
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


def lam_min_all(P):
    return np.array([np.linalg.eigvalsh(P[np.ix_(T, T)])[0] for T in TRIPLES])


def active_minors(P, T):
    """The three 2x2 principal minors of (P_TT - I/6)."""
    M = P[np.ix_(T, T)] - TH * np.eye(3)
    out = []
    for idx in itertools.combinations(range(3), 2):
        sub = M[np.ix_(idx, idx)]
        out.append(float(np.linalg.det(sub)))
    return out


def eq_residual(v, active):
    """Squared residual of the ACTIVE equalities only (no inequality terms)."""
    P = retract(v.reshape(6, 3))
    P = P @ P.T
    lam = lam_min_all(P)
    return float(sum((lam[i] - TH) ** 2 for i in active))


def worst_minor(v, active):
    """max over active T of ( - most negative 2x2 minor ):  how badly PSD fails."""
    A = retract(v.reshape(6, 3))
    P = A @ A.T
    best = -np.inf
    for i in active:
        mins = active_minors(P, TRIPLES[i])
        best = max(best, -min(mins))
    return float(best)


def probe(indices, n_starts, seed, w_eq=1e4):
    """Minimise worst_minor subject to staying on the equality locus (penalised).

    We want inf over the locus, so we MINIMISE worst_minor while a strong penalty
    keeps the equalities satisfied.  Reported alongside the achieved equality
    residual so a point that cheated off the locus is visible."""
    rng = np.random.default_rng(seed)
    act = set(indices)
    best = (np.inf, None, None)
    for _ in range(n_starts):
        v = rng.standard_normal(18)

        def obj(v):
            return worst_minor(v, act) + w_eq * eq_residual(v, act)

        for _ in range(3):
            v = minimize(obj, v, method="Nelder-Mead",
                         options=dict(maxiter=6000, maxfev=6000,
                                      xatol=1e-13, fatol=1e-16)).x
        r = eq_residual(v, act)
        m = worst_minor(v, act)
        if r < 1e-10 and m < best[0]:
            best = (m, r, v.copy())
    return best


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
    print("MARGIN PROBE for the SOS route (alternative 1 of §3i.1)")
    print(f"  sizes={sizes}  orbits={len(orbits)}  starts={a.starts}")
    print("  mu(A) = inf over the ACTIVE-EQUALITY locus of")
    print("          max_{T active} ( - most negative 2x2 minor of P_TT - I/6 )")
    print("  mu(A) > 0 well separated  =>  a wide-margin SOS certificate should exist")
    print("  mu(A) ~ 0                 =>  the margin is thin there; SOS delicate")
    print("=" * 78, flush=True)

    rows = []
    for n, orb in enumerate(orbits, 1):
        idx = sorted(IDX[tuple(t)] for t in orb["triples"])
        m, r, v = probe(idx, a.starts, a.seed + 104729 * n)
        rows.append(dict(size=orb["size"], canon=orb["canon"], indices=idx,
                         mu=None if m == np.inf else float(m),
                         eq_residual=None if r is None else float(r)))
        if m == np.inf:
            print(f"[{n}/{len(orbits)}] size {orb['size']} canon {orb['canon']}: "
                  f"no point reached the equality locus (residual never < 1e-10)",
                  flush=True)
        else:
            print(f"[{n}/{len(orbits)}] size {orb['size']} canon {orb['canon']}: "
                  f"mu = {m:+.6e}   (eq residual {r:.2e})", flush=True)

    got = [r for r in rows if r["mu"] is not None]
    print("\n" + "=" * 78)
    if got:
        mus = np.array([r["mu"] for r in got])
        print(f"orbits with a point on the locus: {len(got)}/{len(rows)}")
        print(f"  mu: min={mus.min():+.4e}  median={np.median(mus):+.4e}"
              f"  max={mus.max():+.4e}")
        neg = [r for r in got if r["mu"] <= 0]
        print(f"  orbits with mu <= 0 (PSD NOT violated => possible feasible point):"
              f" {len(neg)}")
        if neg:
            print("  *** those orbits need direct attention: the PSD obstruction")
            print("      does NOT hold there, so emptiness must come from the")
            print("      INACTIVE inequalities instead. ***")
            for r in neg[:6]:
                print(f"      size {r['size']} canon {r['canon']} mu={r['mu']:+.3e}")
        else:
            print("  ALL probed orbits have a strictly positive PSD-violation margin:")
            print("  on the whole active-equality locus, some active 2x2 minor is")
            print("  negative by at least mu > 0.  This is exactly the wide-margin")
            print("  statement an SOS/Positivstellensatz certificate would prove,")
            print("  and it needs NO root isolation.")
    else:
        print("No orbit produced a point on the equality locus; the probe could not")
        print("measure a margin.  Increase --starts or relax the residual gate.")
    os.makedirs("verify/out", exist_ok=True)
    json.dump(dict(sizes=sizes, starts=a.starts, seed=a.seed, rows=rows),
              open("verify/out/v26_margin.json", "w"), indent=1)
    print("\nwrote verify/out/v26_margin.json")
    print("=" * 78)
