#!/usr/bin/env python3
"""
v3_extremal.py -- identify and algebraically characterize the true minimizers of
F(A) = max_{|T|=3} lam_min(P_TT) over the FULL St(6,3), not just the slice.

Motivation (a finding of this session): the slice minimum is G_PMC = (1-sin36)/2
= 0.20610..., but unconstrained search over St(6,3) reaches F ~ 0.1849, which is
much closer to the 1/6 = 0.16667 threshold.  So the slice is NOT where GTZ(6,3)
is tight, and Track B's target -- while a fine standalone lemma -- is far from
the real extremal structure.  The binding configurations live off the slice.

This script:
  1. runs many starts, clusters the local minima by value,
  2. for the lowest cluster, extracts the leverage profile (diag P), the active
     triple set, and the Gram spectrum,
  3. tries to rationalize / algebraically identify the minimizing value with PSLQ
     (mpmath.pslq / findpoly) to guess its minimal polynomial,
  4. reports the true margin over 1/6.

NOTHING here is PROVED; all output is NUMERICALLY SUPPORTED.  The point is to
localize the real extremal structure so a certificate can target it.

Deterministic: master seed fixed.
"""
import itertools, sys, os, json
import numpy as np
from scipy.optimize import minimize
from multiprocessing import Pool
from mpmath import mp, mpf, findpoly, nstr, identify

MASTER = 20260731
TRIPLES = list(itertools.combinations(range(6), 3))
THRESH = 1.0 / 6.0
mp.dps = 40


def retract(X):
    U, _, Vt = np.linalg.svd(X, full_matrices=False)
    return U @ Vt


def lams(A):
    P = A @ A.T
    return np.array([np.linalg.eigvalsh(P[np.ix_(T, T)])[0] for T in TRIPLES]), P


def F(v):
    return float(np.max(lams(retract(v.reshape(6, 3)))[0]))


def softF(v, beta):
    """Numerically stable log-sum-exp (shift by the max; naive exp overflows for
    beta >~ 700/max|lam| and silently poisons the optimizer with inf/nan)."""
    l = lams(retract(v.reshape(6, 3)))[0]
    m = float(np.max(l))
    return m + float(np.log(np.sum(np.exp(beta * (l - m)))) / beta)


def run(seed):
    rng = np.random.default_rng(seed)
    v = rng.standard_normal(18)
    for beta in (30.0, 150.0, 800.0, 4000.0):
        v = minimize(softF, v, args=(beta,), method="Nelder-Mead",
                     options=dict(maxiter=6000, maxfev=6000,
                                  xatol=1e-12, fatol=1e-15)).x
    for _ in range(4):
        v = minimize(F, v, method="Nelder-Mead",
                     options=dict(maxiter=8000, maxfev=8000,
                                  xatol=1e-13, fatol=1e-16)).x
    return float(F(v)), v.tolist(), int(seed)


if __name__ == "__main__":
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 2000
    ncpu = max(1, min(20, os.cpu_count() - 2))
    seeds = [int(x) for x in np.random.SeedSequence(MASTER).generate_state(n, dtype=np.uint32)]
    print(f"master={MASTER} starts={n} cores={ncpu}", flush=True)

    with Pool(ncpu) as pool:
        out = pool.map(run, seeds, chunksize=4)
    out.sort(key=lambda z: z[0])

    Fmin = out[0][0]
    print(f"\nglobal min F  = {Fmin!r}")
    print(f"1/6           = {THRESH!r}")
    print(f"margin        = {Fmin - THRESH:+.8e}")

    viol = [z for z in out if z[0] < THRESH - 1e-6]
    print(f"counterexamples (F < 1/6 - 1e-6): {len(viol)}")

    # cluster local minima
    vals = sorted(z[0] for z in out)
    clusters, cur = [], [vals[0]]
    for x in vals[1:]:
        if x - cur[-1] < 1e-7:
            cur.append(x)
        else:
            clusters.append(cur); cur = [x]
    clusters.append(cur)
    print(f"\ndistinct local minima: {len(clusters)}")
    for c in clusters[:10]:
        print(f"   value {c[0]!r}   multiplicity {len(c)}")

    # ---- structure of the global minimizer
    A = retract(np.array(out[0][1]).reshape(6, 3))
    l, P = lams(A)
    lev = np.diag(P)
    order = np.argsort(-l)
    active = [TRIPLES[i] for i in order if abs(l[i] - Fmin) < 1e-9]
    print(f"\n--- structure of the global minimizer (seed {out[0][2]}) ---")
    print(f"leverages diag(P)   = {np.round(np.sort(lev), 10)}")
    print(f"  sum = {lev.sum():.12f} (must be 3)   all > 1/6? {bool(np.all(lev > 1/6))}")
    print(f"  on the slice (all=1/2)? {bool(np.allclose(lev, 0.5, atol=1e-7))}")
    print(f"active triples ({len(active)}): {active}")
    print(f"top 6 lam_min values: {np.round(np.sort(l)[-6:], 12)}")
    Pev = np.linalg.eigvalsh(P)
    print(f"spec(P) = {np.round(Pev, 10)}  (should be 1,1,1,0,0,0)")
    off = sorted(abs(P[i, j]) for i, j in itertools.combinations(range(6), 2))
    print(f"|c_ij| sorted = {np.round(off, 8)}")

    # ---- try to identify Fmin algebraically
    print("\n--- algebraic identification of F_min (PSLQ / findpoly) ---")
    xm = mpf(repr(Fmin))
    for deg in (1, 2, 3, 4, 5, 6, 8):
        pc = findpoly(xm, deg, maxcoeff=10**6, tol=mpf(10) ** -13)
        if pc:
            print(f"  degree {deg}: candidate minpoly coeffs {pc}")
            break
    else:
        print("  no low-degree integer polynomial found at 1e-13 tolerance")
    idn = identify(xm)
    print(f"  identify() -> {idn}")

    os.makedirs("verify/out", exist_ok=True)
    with open("verify/out/v3_extremal.json", "w") as fh:
        json.dump(dict(master=MASTER, n=n, Fmin=Fmin, margin=Fmin - THRESH,
                       n_clusters=len(clusters),
                       cluster_values=[c[0] for c in clusters[:40]],
                       cluster_mult=[len(c) for c in clusters[:40]],
                       leverages=sorted(lev.tolist()),
                       active_triples=active, spec_P=Pev.tolist(),
                       offdiag_abs=off, A=A.tolist(), seed=out[0][2],
                       n_violations=len(viol)), fh, indent=1)
    print("\nwrote verify/out/v3_extremal.json")

    if viol:
        print("\n*** COUNTEREXAMPLE PROTOCOL TRIGGERED -- HALT AND ESCALATE ***")
        sys.exit(2)
