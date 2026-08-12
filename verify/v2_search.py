#!/usr/bin/env python3
"""
v2_search.py -- large-scale minimization of F(A)=max_{|T|=3} lam_min(P_TT) over
(a) the equal-leverage slice and (b) all of St(6,3), to firm up Track B's target
and to look for counterexamples to GTZ(6,3).

Per `GTZ 6 3 Closure Attempt.md` sec.0 rule 3 (COUNTEREXAMPLE PROTOCOL): any sample
with F < 1/6 - 1e-6 halts everything and is reported, never averaged away.

Everything is deterministic: seeds are derived from a fixed master seed.
Results are NUMERICALLY SUPPORTED by construction -- this script certifies nothing.

Usage:  .venv/bin/python -u verify/v2_search.py [n_starts_slice] [n_starts_full]
"""
import itertools, sys, os, json
import numpy as np
from scipy.optimize import minimize
from multiprocessing import Pool

MASTER_SEED = 20260730
TRIPLES = list(itertools.combinations(range(6), 3))
THRESH = 1.0 / 6.0
CX_TOL = 1e-6                       # counterexample tolerance from the brief
G_PMC = 0.20610737385376343541564702268   # (1 - sin36)/2, verified exactly in v1


# ------------------------------------------------------------------ objective
def F_from_A(A):
    """F(A) = max_T lam_min(P_TT), P = A A^T."""
    P = A @ A.T
    best = -np.inf
    for T in TRIPLES:
        blk = P[np.ix_(T, T)]
        lm = np.linalg.eigvalsh(blk)[0]
        if lm > best:
            best = lm
    return best


def retract(X):
    """Nearest point of St(6,3) (polar retraction)."""
    U, _, Vt = np.linalg.svd(X, full_matrices=False)
    return U @ Vt


# ---- full St(6,3): minimize F over A, parametrized by an unconstrained 6x3
def negF_full(v):
    A = retract(v.reshape(6, 3))
    return F_from_A(A)


# ---- slice: A in St(6,3) with diag(A A^T) == 1/2.
# Parametrize by an unconstrained 6x3 V, orthonormalize columns EXACTLY via polar
# retraction (so A^T A = I_3 holds to machine precision, no penalty needed), then
# penalize ONLY the 6 diagonal deviations diag(P) - 1/2.  This is far better
# conditioned than penalizing a frame condition after row-renormalization.
#
# The max objective is nonsmooth; we smooth it with log-sum-exp at temperature
# beta and anneal beta upward, polishing with Nelder-Mead at the end.
def _lams(A):
    P = A @ A.T
    return np.array([np.linalg.eigvalsh(P[np.ix_(T, T)])[0] for T in TRIPLES]), P


def project_slice(v, iters=250):
    """Alternating projection onto  St(6,3) INTERSECT {diag(A A^T) = 1/2}.
    Row-scaling to norm^2=1/2 then polar retraction; converges fast in practice."""
    A = retract(v.reshape(6, 3))
    for _ in range(iters):
        nrm = np.linalg.norm(A, axis=1, keepdims=True)
        nrm[nrm < 1e-300] = 1e-300
        A = A * (np.sqrt(0.5) / nrm)          # force row norms to 1/sqrt2
        A = retract(A)                        # back onto the Stiefel manifold
    return A


def slice_obj(v, pen=0.0, beta=None):
    """Objective ON the slice: project first, so feasibility is never traded away.
    log-sum-exp is shifted by the max -- the naive form overflows for large beta."""
    A = project_slice(v)
    lam, _ = _lams(A)
    if beta is None:
        return float(np.max(lam))
    m = float(np.max(lam))
    return m + float(np.log(np.sum(np.exp(beta * (lam - m)))) / beta)


def slice_G(v):
    """True slice G and the residual AFTER projection."""
    A = project_slice(v)
    lam, P = _lams(A)
    return float(np.max(lam)), float(np.linalg.norm(np.diag(P) - 0.5))


# ------------------------------------------------------------------- workers
def run_slice(seed):
    rng = np.random.default_rng(seed)
    v = rng.standard_normal(18)
    # annealed smoothing of the nonsmooth max: soft -> hard
    for beta in (40.0, 200.0, 1000.0, None):
        r = minimize(slice_obj, v, args=(0.0, beta), method="Nelder-Mead",
                     options=dict(maxiter=4000, maxfev=4000, xatol=1e-11, fatol=1e-14))
        v = r.x
    G, resid = slice_G(v)
    return G, resid, v.tolist(), int(seed)


def run_full(seed):
    rng = np.random.default_rng(seed)
    v0 = rng.standard_normal(18)
    r = minimize(negF_full, v0, method="Nelder-Mead",
                 options=dict(maxiter=20000, maxfev=20000, xatol=1e-12, fatol=1e-15))
    r = minimize(negF_full, r.x, method="L-BFGS-B",
                 options=dict(maxiter=4000, ftol=1e-16, gtol=1e-14))
    A = retract(r.x.reshape(6, 3))
    return float(F_from_A(A)), A.tolist(), int(seed)


# ---------------------------------------------------------------------- main
if __name__ == "__main__":
    n_slice = int(sys.argv[1]) if len(sys.argv) > 1 else 4000
    n_full = int(sys.argv[2]) if len(sys.argv) > 2 else 4000
    ncpu = max(1, min(20, os.cpu_count() - 2))
    ss = np.random.SeedSequence(MASTER_SEED)
    seeds_s = [int(x) for x in ss.spawn(1)[0].generate_state(n_slice, dtype=np.uint32)]
    seeds_f = [int(x) for x in ss.spawn(1)[0].generate_state(n_full, dtype=np.uint32)]

    print(f"master seed = {MASTER_SEED};  cores = {ncpu}")
    print(f"G_PMC (exact, from v1) = {G_PMC!r}")
    print(f"threshold 1/6          = {THRESH!r}")
    print("=" * 74)

    # ---------------- (A) slice search
    print(f"[A] slice minimization, {n_slice} random starts ...", flush=True)
    with Pool(ncpu) as pool:
        out_s = pool.map(run_slice, seeds_s, chunksize=8)
    valid = [(G, r, v, sd) for G, r, v, sd in out_s if r < 1e-7]
    print(f"    {len(valid)}/{n_slice} starts converged onto the slice (resid<1e-7)")
    if valid:
        valid.sort(key=lambda z: z[0])
        Gmin, rmin, vmin, sdmin = valid[0]
        print(f"    slice min G      = {Gmin!r}   (frame resid {rmin:.2e}, seed {sdmin})")
        print(f"    G_PMC            = {G_PMC!r}")
        print(f"    G_min - G_PMC    = {Gmin - G_PMC:+.3e}")
        below = [z for z in valid if z[0] < G_PMC - 1e-9]
        print(f"    starts landing strictly BELOW G_PMC (by >1e-9): {len(below)}")
        near = [z for z in valid if abs(z[0] - G_PMC) < 1e-7]
        print(f"    starts landing AT G_PMC (within 1e-7):          {len(near)}")
        # histogram of distinct local minima
        lv = sorted({round(z[0], 9) for z in valid})
        print(f"    distinct local minima (9dp): {len(lv)}; smallest 8: {lv[:8]}")

    # ---------------- (B) full St(6,3) search + counterexample watch
    print(f"\n[B] full St(6,3) minimization, {n_full} random starts ...", flush=True)
    with Pool(ncpu) as pool:
        out_f = pool.map(run_full, seeds_f, chunksize=8)
    out_f.sort(key=lambda z: z[0])
    Fmin, Amin, sdf = out_f[0]
    print(f"    global min F over all starts = {Fmin!r}  (seed {sdf})")
    print(f"    1/6                          = {THRESH!r}")
    print(f"    F_min - 1/6                  = {Fmin - THRESH:+.6e}")

    viol = [z for z in out_f if z[0] < THRESH - CX_TOL]
    print(f"    samples with F < 1/6 - 1e-6  = {len(viol)}   <-- counterexample count")

    lvf = sorted({round(z[0], 9) for z in out_f})
    print(f"    distinct local minima (9dp): {len(lvf)}; smallest 8: {lvf[:8]}")

    # save artifacts
    os.makedirs("verify/out", exist_ok=True)
    with open("verify/out/v2_results.json", "w") as fh:
        json.dump(dict(master_seed=MASTER_SEED, n_slice=n_slice, n_full=n_full,
                       G_PMC=G_PMC, threshold=THRESH,
                       slice_min=(valid[0][0] if valid else None),
                       slice_min_seed=(valid[0][3] if valid else None),
                       slice_converged=len(valid),
                       full_min=Fmin, full_min_seed=sdf, full_min_A=Amin,
                       n_violations=len(viol),
                       slice_local_minima=lv[:40] if valid else [],
                       full_local_minima=lvf[:40]), fh, indent=1)
    print("\n    wrote verify/out/v2_results.json")

    print("=" * 74)
    if viol:
        print("*** COUNTEREXAMPLE PROTOCOL TRIGGERED (brief sec.0 rule 3) ***")
        print(f"*** {len(viol)} sample(s) with F < 1/6 - 1e-6. HALT AND ESCALATE. ***")
        print(f"*** worst F = {viol[0][0]!r} at seed {viol[0][2]} ***")
        sys.exit(2)
    else:
        print("No counterexample found: every sample has F >= 1/6 - 1e-6.")
        print("GTZ(6,3) survives this search.  Status: NUMERICALLY SUPPORTED (not proved).")
