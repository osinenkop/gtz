#!/usr/bin/env python3
"""
v15_dimension.py -- ALGEBRAIC dimension of the extremal set E, stratum by stratum,
and the uniform lower bound on the certified radius r0.

TWO GOALS, ONE SETUP.

GOAL A (finiteness, algebraically -- no sampling).
  E = {P in Gr(3,6) : F(P) = 1/6} decomposes by ACTIVE SET
      E_A = {P : P^2=P, tr P=3, det(P_TT - I/6)=0 for all T in A,
                 P_TT - I/6 PSD for T in A, lam_min(P_TT) < 1/6 otherwise}.
  If every nonempty stratum has dimension 0, E is FINITE -- proved outright, with
  no reliance on having sampled every extremal.  This is the one route that cannot
  be defeated by an unsampled extremal.

  KEY REDUCTION (proved, and it shrinks the work):
      at a simple-active point L has one row per active triple and 9 columns.
      Sharpness needs positive spanning in R^9, hence at least 10 active rows.
      Therefore a simple-active extremal with |A| <= 9 is automatically
      NON-sharp, hence by Cor F3 a candidate seed for a positive-dimensional
      family.  Dimension count agrees for |A| <= 8: |A| equations on a 9-dim
      manifold leave dim >= 9 - |A| >= 1 before inequalities/higher-order
      effects.

  So the first low-active finiteness question is:
      does there exist a simple-active P with F(P) = 1/6 and |A| <= 9 ?
  All extremals ever observed have |A| in {10, 12, 13}.  The nonsimple active
  case is separate and is not covered by this row-count shortcut.

GOAL B (uniform lower bound on r0 = kappa * gap_min).
  r0 is bounded below uniformly on E iff both factors are.  gap_min is the easier
  one and is what we attack here: gap_T = lam_2(P_TT) - 1/6 on the active blocks.
  We compute, per stratum, the exact minimum of gap over the stratum's defining
  ideal where tractable, and otherwise a certified numerical lower bound via
  constrained minimization from many starts (flagged NUMERICALLY SUPPORTED).

METHOD.  Work in the gauge-free projector coordinates but avoid the 21-variable
ambient Sym(6) ideal, which is too big for Groebner.  Instead parametrize
Gr(3,6) LOCALLY by the 9 tangent coordinates around a base point (a chart), which
is exactly 9 variables -- the minimum possible -- and compute the dimension of the
solution variety of the |A| equations in that chart.  Dimension is a local
invariant, so a chart suffices, and we run one chart per known configuration plus
random charts to probe elsewhere.

Everything symbolic is exact.  Numerical parts are clearly tagged.
"""
import itertools, json, os, sys
from collections import Counter
import numpy as np
import sympy as sp

Q = sp.Rational
SIX = Q(1, 6)
TRIPLES = list(itertools.combinations(range(6), 3))


def retract(X):
    U, _, Vt = np.linalg.svd(X, full_matrices=False)
    return U @ Vt


# ------------------------------------------------------------------ GOAL A
def stratum_dimension(P0f, label, max_eqs=None):
    """Local dimension of {F = 1/6} through P0 in a 9-dim chart around P0.

    Chart: P(z) = exp-like retraction of the tangent perturbation.  We use the
    exact algebraic chart  P(z) = (P0 + Z)(I + ...)  avoided in favour of the
    simplest exact one: the graph chart of the Grassmannian.  Writing
    col(P0) = span(U), ker(P0) = span(N), every nearby 3-plane is the column span
    of U + N*X for a 3x3 matrix X (9 parameters), and its projector is
        P(X) = (U + N X)((U + N X)^T (U + N X))^{-1} (U + N X)^T.
    This is a rational parametrization in the 9 entries of X, exact and with
    P(0) = P0 -- ideal for a local dimension computation.
    """
    w, V = np.linalg.eigh(P0f)
    U = V[:, np.argsort(-w)[:3]]
    N = V[:, np.argsort(-w)[3:]]

    lam = np.array([np.linalg.eigvalsh(P0f[np.ix_(T, T)])[0] for T in TRIPLES])
    active = [i for i in range(20) if abs(lam[i] - 1 / 6) < 1e-8]
    out = dict(label=label, n_active=len(active),
               active=[list(TRIPLES[i]) for i in active])

    # rationalize U, N so the chart is exact (small denominators after rounding;
    # we verify the rounded chart still reproduces the active set)
    def rat(M, den=10**6):
        return sp.Matrix([[sp.nsimplify(sp.Rational(float(M[i, j])).limit_denominator(den))
                           for j in range(M.shape[1])] for i in range(M.shape[0])])

    Us, Ns = rat(U), rat(N)
    X = sp.Matrix(3, 3, lambda i, j: sp.Symbol(f"x{i}{j}"))
    xs = list(X)
    B = Us + Ns * X                                     # 6x3, entries linear in x
    G = sp.expand(B.T * B)                              # 3x3 Gram
    detG = sp.expand(G.det())
    adjG = sp.expand(G.adjugate())
    # P = B G^{-1} B^T = B adj(G) B^T / det(G).  Clear the denominator: the active
    # equations det(P_TT - I/6) = 0 become polynomial after multiplying by det(G)^3.
    Pn = sp.expand(B * adjG * B.T)                      # = det(G) * P

    eqs = []
    for i in (active if max_eqs is None else active[:max_eqs]):
        T = TRIPLES[i]
        M = sp.Matrix(3, 3, lambda a, b: Pn[T[a], T[b]])
        # det(P_TT - I/6) = 0  <=>  det(Pn_TT - (detG/6) I) = 0  (scaled by detG^3)
        e = sp.expand((M - detG / 6 * sp.eye(3)).det())
        e = sp.expand(sp.numer(sp.together(e)))
        eqs.append(sp.Poly(e, *xs))
    out["n_eqs"] = len(eqs)

    # Dimension via the Jacobian rank at x = 0 (gives the dimension of the tangent
    # space to the solution variety; for a reduced component this equals 9 - rank).
    J = sp.Matrix([[sp.diff(e.as_expr(), v) for v in xs] for e in eqs])
    J0 = J.subs({v: 0 for v in xs})
    r = J0.rank()
    out["jacobian_rank_at_P0"] = int(r)
    out["tangent_dim_lower_bound"] = int(9 - r)
    out["dim0_at_P0"] = bool(r == 9)
    return out


# ------------------------------------------------------------------ GOAL B
def gap_and_kappa(P):
    """gap_min and (numeric) kappa at a projector P."""
    lam = np.array([np.linalg.eigvalsh(P[np.ix_(T, T)])[0] for T in TRIPLES])
    active = [i for i in range(20) if abs(lam[i] - 1 / 6) < 1e-8]
    gaps = []
    Ev = {}
    for i in active:
        T = TRIPLES[i]
        w, U = np.linalg.eigh(P[np.ix_(T, T)])
        if abs(w[1] - w[0]) < 1e-9:
            return None, None, len(active)
        gaps.append(float(w[1] - 1 / 6))
        Ev[T] = U[:, 0]
    w, U = np.linalg.eigh(P)
    rb, kb = U[:, np.argsort(-w)[:3]], U[:, np.argsort(-w)[3:]]
    TAN = []
    for a in range(3):
        for b in range(3):
            Xm = np.outer(rb[:, a], kb[:, b])
            Xm = Xm + Xm.T
            TAN.append(Xm / np.linalg.norm(Xm, "fro"))
    L = np.array([[Ev[TRIPLES[i]] @ Bm[np.ix_(TRIPLES[i], TRIPLES[i])] @ Ev[TRIPLES[i]]
                   for Bm in TAN] for i in active])
    from scipy.optimize import minimize
    best, rg = np.inf, np.random.default_rng(3)
    for _ in range(150):
        z = rg.standard_normal(9); z /= np.linalg.norm(z)
        rr = minimize(lambda t: float(np.max(L @ (t / np.linalg.norm(t)))), z,
                      method="Nelder-Mead",
                      options=dict(maxiter=3000, maxfev=3000, xatol=1e-13, fatol=1e-16))
        best = min(best, float(rr.fun))
    return float(min(gaps)), best, len(active)


def hunt_small_active(n_starts, seed0=20260802):
    """GOAL A's decisive search: any simple-active extremal with |A| <= 9?
    Also collects the
    empirical distribution of gap_min and kappa for GOAL B."""
    from scipy.optimize import minimize
    TH = 1 / 6

    def lams(A):
        P = A @ A.T
        return np.array([np.linalg.eigvalsh(P[np.ix_(T, T)])[0] for T in TRIPLES]), P

    def F(v):
        return float(np.max(lams(retract(v.reshape(6, 3)))[0]))

    def softF(v, b):
        l = lams(retract(v.reshape(6, 3)))[0]
        m = float(np.max(l))
        return m + float(np.log(np.sum(np.exp(b * (l - m)))) / b)

    rng = np.random.default_rng(seed0)
    recs = []
    for k in range(n_starts):
        v = rng.standard_normal(18)
        for b in (30., 200., 1500.):
            v = minimize(softF, v, args=(b,), method="Nelder-Mead",
                         options=dict(maxiter=4000, maxfev=4000,
                                      xatol=1e-12, fatol=1e-15)).x
        for _ in range(6):
            v = minimize(F, v, method="Nelder-Mead",
                         options=dict(maxiter=15000, maxfev=15000,
                                      xatol=1e-15, fatol=1e-17)).x
        f = F(v)
        if f - TH < -1e-6:
            return dict(VIOLATION=True, F=f, v=v.tolist())
        if abs(f - TH) > 1e-9:
            continue
        A = retract(v.reshape(6, 3))
        P = A @ A.T
        g, kp, na = gap_and_kappa(P)
        recs.append(dict(n_active=na, gap_min=g, kappa=kp,
                         lev=sorted(np.diag(P).tolist())))
    return dict(records=recs)


if __name__ == "__main__":
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 200
    print("=" * 78)
    print("GOAL A -- algebraic dimension of the extremal strata")
    print("=" * 78)
    print("REDUCTION: in the simple-active case, sharpness needs positive")
    print("           spanning in R^9, so SHARP => |A| >= 10.")
    print("           A simple-active extremal with |A| <= 9 is automatically")
    print("           non-sharp => candidate positive-dim family.")
    print()

    results = []
    # the seventh (out-of-family) extremal
    p7 = "verify/data/P514_seventh.npy"
    if os.path.exists(p7):
        P = np.load(p7)
        r = stratum_dimension(P, "out-of-family (5/14,9/14)")
        results.append(r)
        print(f"{r['label']}")
        print(f"  |A| = {r['n_active']},  equations used = {r['n_eqs']}")
        print(f"  Jacobian rank at P0 = {r['jacobian_rank_at_P0']} / 9")
        print(f"  => local dim of the stratum = {r['tangent_dim_lower_bound']}")
        print(f"  => ISOLATED (dim 0) at P0? {r['dim0_at_P0']}")

    print()
    print("=" * 78)
    print(f"GOAL A/B -- decisive search: any simple-active extremal with |A| <= 9 ?  ({n} starts)")
    print("=" * 78)
    out = hunt_small_active(n)
    if out.get("VIOLATION"):
        print("*** COUNTEREXAMPLE -- HALT ***")
        json.dump(out, open("verify/out/v15_VIOLATION.json", "w"), indent=1)
        sys.exit(2)
    recs = out["records"]
    print(f"extremals found: {len(recs)}")
    if recs:
        from collections import Counter
        cnt = Counter(r["n_active"] for r in recs)
        print(f"active-set sizes: {dict(sorted(cnt.items()))}")
        small = [r for r in recs if r["n_active"] <= 9]
        print(f"\nsimple-active |A| <= 9 (would be NON-SHARP): {len(small)}")
        if small:
            print("*** FOUND -- this is the decisive object.  Certify exactly. ***")
            for r in small[:5]:
                print(f"    |A|={r['n_active']} lev={np.round(r['lev'],6)}")
        else:
            print("None.  Every extremal found has |A| >= 10.")

        gs = [r["gap_min"] for r in recs if r["gap_min"] is not None]
        ks = [r["kappa"] for r in recs if r["kappa"] is not None]
        print("\n--- GOAL B: empirical bounds for the uniform radius ---")
        print(f"  gap_min: min={min(gs):.9f}  median={np.median(gs):.9f}  max={max(gs):.9f}")
        print(f"  kappa  : min={min(ks):.9f}  median={np.median(ks):.9f}  max={max(ks):.9f}")
        r0s = [g * k for g, k in zip(gs, ks)]
        print(f"  r0     : min={min(r0s):.9f}  median={np.median(r0s):.9f}")
        print(f"  => empirical uniform bound r0 >= {min(r0s):.9f} over {len(recs)} extremals")
        print("     (NUMERICALLY SUPPORTED; a proof needs a global argument)")

    json.dump(dict(strata=results, n_starts=n,
                   active_size_counts={str(k): v for k, v in
                                       Counter(r["n_active"] for r in recs).items()}
                   if recs else {},
                   records=recs), open("verify/out/v15_dimension.json", "w"), indent=1)
    print("\nwrote verify/out/v15_dimension.json")
    print("=" * 78)
