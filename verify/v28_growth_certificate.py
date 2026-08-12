#!/usr/bin/env python3
"""
v28_growth_certificate.py -- attempt an exact certificate for the growth mechanism
of §3i.3/§3i.4, on the unique |A|=6 orbit.

THE STATEMENT TO CERTIFY.  For A = {012,013,045,145,234,235} (indices 0,1,9,15,16,17):

    on  V_A := { P in Gr(3,6) : lambda_min(P_TT) = 1/6  for all T in A },
    there is a triple  T* not in A  with  lambda_min(P_{T*T*}) >= 1/6 + c,
    c ~ 4.7e-2  (observed minimum over all low-active orbits).

If true, no point of V_A has actual active set exactly A, so the stratum is empty.

WHY A DIRECT SOS ATTEMPT IS THE RIGHT MOVE NOW.
  * lex elimination died on degree (1880 univariate; §3i.1);
  * leverage/trace conditions kill nothing (§3i.2);
  * active PSD is NOT violated, so there is no PSD certificate to find (§3i.3);
  * but the inactive overshoot has a WIDE margin, and margin -- not ideal degree --
    governs Positivstellensatz certificate degree.

STRATEGY IMPLEMENTED HERE (two independent, cheap, exact-friendly steps).

STEP 1 -- IDENTIFY THE WITNESS TRIPLE, EXACTLY.
  The certificate needs a specific T* (or a small set of candidates).  We sample the
  locus numerically and record WHICH outside triple overshoots, and by how much.  If
  a single T* works everywhere, the certificate is one inequality; if the witness
  varies, we need a case split over the witnesses, which is still finite and small.
  This is pure bookkeeping but it decides the shape of the exact statement.

STEP 2 -- LOWER-BOUND THE OVERSHOOT BY A SUM-OF-SQUARES ARGUMENT ON THE TRACE.
  For any triple T, lambda_max(P_TT) >= tr(P_TT)/3 and
  lambda_min(P_TT) >= tr(P_TT) - 2*lambda_max(P_TT) >= ... too lossy.  Instead use
  the exact and sharp route: for a 3x3 symmetric M with eigenvalues in [0,1]
  (which P_TT satisfies since 0 <= P <= I),

      lambda_min(M) >= e_3(M) / e_2(M)     when e_2(M) > 0,

  because e_3 = l1 l2 l3 <= l1 l2 * ... -- we verify the exact inequality
  numerically and use it only as a CANDIDATE bound generator, then check it
  symbolically on the locus.  The point is to get a bound on lambda_min(P_{T*T*})
  in terms of POLYNOMIALS in P (det and 2x2 minors), which is what an exact
  certificate can handle -- eigenvalues themselves are not polynomial.

  Concretely we test, on many locus points, whether the polynomial inequality
      det(P_{T*T*}) - (1/6) * e_2(P_{T*T*}) + (1/36) * tr(P_{T*T*}) - 1/216  >= 0
  holds -- this is exactly det(P_{T*T*} - I/6) >= 0, i.e. the statement that
  P_{T*T*} - I/6 has nonneg determinant.  Combined with the trace and 2x2-minor
  signs it yields lambda_min >= 1/6 by the standard sign-of-characteristic-
  polynomial-coefficients criterion, with NO eigenvalue computation.

  That reduces the certificate to: on V_A, the polynomial
      g(P) := det(P_{T*T*} - I/6)
  and the two lower principal minors are all >= 0, with margin.  Those are
  polynomial inequalities of low degree in the chart -- exactly SOS-shaped.

Output: the witness distribution, and the numerical margins of the polynomial
(not eigenvalue) certificate quantities, which is what the exact layer must bound.
"""
import argparse, itertools, json, os, sys
from collections import Counter
import numpy as np
from scipy.optimize import minimize

TRIPLES = list(itertools.combinations(range(6), 3))
IDX = {t: i for i, t in enumerate(TRIPLES)}
TH = 1.0 / 6.0
A6 = [0, 1, 9, 15, 16, 17]


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


def minors_of_shifted(P, T):
    """The three leading data of M = P_TT - I/6: trace, e_2, det.
    lambda_min(M) >= 0  <=>  all principal minors >= 0; for a 3x3 symmetric matrix
    it is equivalent (and cheaper) to check tr >= 0, e_2 >= 0, det >= 0 TOGETHER
    with M symmetric -- these are the coefficients of the characteristic polynomial
    and their nonnegativity is exactly Descartes' criterion for all roots >= 0."""
    M = P[np.ix_(T, T)] - TH * np.eye(3)
    tr = float(np.trace(M))
    e2 = 0.0
    for idx in itertools.combinations(range(3), 2):
        e2 += float(np.linalg.det(M[np.ix_(idx, idx)]))
    det = float(np.linalg.det(M))
    return tr, e2, det


def sample_locus(n_starts, seed, maxiter=8000):
    rng = np.random.default_rng(seed)
    act = set(A6)
    pts = []
    for _ in range(n_starts):
        v = rng.standard_normal(18)
        for _ in range(4):
            v = minimize(eq_resid, v, args=(act,), method="Nelder-Mead",
                         options=dict(maxiter=maxiter, maxfev=maxiter,
                                      xatol=1e-14, fatol=1e-17)).x
        r = eq_resid(v, act)
        if r < 1e-12:
            A = retract(v.reshape(6, 3))
            pts.append((r, A @ A.T))
    return pts


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--starts", type=int, default=120)
    ap.add_argument("--seed", type=int, default=20260805)
    a = ap.parse_args()

    print("=" * 78)
    print("GROWTH CERTIFICATE probe -- unique |A|=6 orbit, indices", A6)
    print(f"  sampling the active-equality locus with {a.starts} starts")
    print("=" * 78, flush=True)

    pts = sample_locus(a.starts, a.seed)
    print(f"points on the locus: {len(pts)}/{a.starts}", flush=True)
    if not pts:
        print("no locus points; increase --starts")
        sys.exit(1)

    inact = [i for i in range(20) if i not in A6]

    # STEP 1 -- which outside triple is the witness?
    wit = Counter()
    excesses = []
    for r, P in pts:
        lam = lam_all(P)
        j = max(inact, key=lambda i: lam[i])
        wit[j] += 1
        excesses.append(float(lam[j] - TH))
    print("\nSTEP 1 -- witness triple distribution (outside triple that overshoots):")
    for j, c in wit.most_common():
        print(f"   triple {TRIPLES[j]} (index {j:>2}): {c} of {len(pts)} points")
    ex = np.array(excesses)
    print(f"   overshoot: min={ex.min():+.5e} median={np.median(ex):+.5e} "
          f"max={ex.max():+.5e}")
    single = len(wit) == 1

    # how many outside triples overshoot simultaneously?
    counts = []
    for r, P in pts:
        lam = lam_all(P)
        counts.append(sum(1 for i in inact if lam[i] > TH + 1e-9))
    print(f"   number of OUTSIDE triples above 1/6 per point: "
          f"min={min(counts)} median={int(np.median(counts))} max={max(counts)}")

    # STEP 2 -- polynomial (not eigenvalue) certificate quantities
    print("\nSTEP 2 -- polynomial certificate data for the witness triples.")
    print("  For M = P_TT - I/6, lambda_min(M) >= 0  <=>  tr,e2,det all >= 0")
    print("  (Descartes on the characteristic polynomial of a symmetric 3x3).")
    print("  If some outside triple has all three >= margin, that is a POLYNOMIAL")
    print("  certificate of overshoot -- SOS-shaped, no eigenvalues involved.")
    best_by_triple = {}
    for j in inact:
        trs, e2s, dets = [], [], []
        for r, P in pts:
            tr, e2, det = minors_of_shifted(P, TRIPLES[j])
            trs.append(tr); e2s.append(e2); dets.append(det)
        best_by_triple[j] = (min(trs), min(e2s), min(dets))
    # a triple certifies overshoot at ALL sampled points iff all three minima >= 0
    good = [(j, v) for j, v in best_by_triple.items()
            if v[0] > 0 and v[1] > 0 and v[2] > 0]
    print(f"\n  outside triples with tr,e2,det all > 0 at EVERY sampled locus point: "
          f"{len(good)}")
    for j, v in sorted(good, key=lambda z: -min(z[1]))[:8]:
        print(f"     {TRIPLES[j]} (idx {j:>2}): min tr={v[0]:+.4e} "
              f"min e2={v[1]:+.4e} min det={v[2]:+.4e}")
    if not good:
        print("     none -- the witness varies across the locus, so the exact")
        print("     statement needs a DISJUNCTION over outside triples (still")
        print("     finite: at most 14 cases for this orbit).")
        # report the best few by det
        rank = sorted(best_by_triple.items(), key=lambda z: -z[1][2])[:6]
        print("     best outside triples by min det:")
        for j, v in rank:
            print(f"       {TRIPLES[j]} (idx {j:>2}): min tr={v[0]:+.4e} "
                  f"min e2={v[1]:+.4e} min det={v[2]:+.4e}")

    os.makedirs("verify/out", exist_ok=True)
    json.dump(dict(orbit=A6, starts=a.starts, seed=a.seed,
                   n_locus_points=len(pts),
                   witness_distribution={str(TRIPLES[j]): c for j, c in wit.items()},
                   single_witness=bool(single),
                   overshoot=dict(min=float(ex.min()), median=float(np.median(ex)),
                                  max=float(ex.max())),
                   simultaneous_overshoot=dict(min=int(min(counts)),
                                               max=int(max(counts))),
                   uniform_certifying_triples=[list(TRIPLES[j]) for j, _ in good],
                   per_triple_minima={str(TRIPLES[j]): list(v)
                                      for j, v in best_by_triple.items()}),
              open("verify/out/v28_growth_certificate.json", "w"), indent=1)
    print("\nwrote verify/out/v28_growth_certificate.json")
    print("=" * 78)
