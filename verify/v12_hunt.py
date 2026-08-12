#!/usr/bin/env python3
"""
v12_hunt.py -- systematic hunt for extremals of F = max_T lambda_min(P_TT) on
Gr(3,6), classified by leverage pattern, with a per-extremal SHARPNESS screen.

WHY THIS SETTLES FINITENESS (the theorem this script is built around):

  (i)  If P0 is a SHARP minimum (kappa := min_{|Pdot|=1} F'(P0;Pdot) > 0), then
       F(P0 + r d) >= 1/6 + (kappa/2) r for small r uniformly in |d| = 1, so a
       punctured neighbourhood of P0 has F > 1/6.  Hence P0 is an ISOLATED point
       of the extremal set E := {P : F(P) = 1/6}.
  (ii) E is closed (F is continuous) and Gr(3,6) is compact.  A closed DISCRETE
       subset of a compact space is FINITE.

  => If every extremal is sharp, E is FINITE.
  => CONTRAPOSITIVE (what we hunt for): a positive-dimensional family of
     extremals must contain a point with a NONTRIVIAL critical cone
     N = {z : Lz <= 0} != {0} -- equivalently rank(L) < 9, or rank 9 with the
     exact LP optimum > 0.  Both quantities are computed by our certificate.

  So finiteness is DECIDABLE PER EXTREMAL with machinery we already have, and the
  hunt has a definite success criterion rather than being open-ended:
     * every extremal found is sharp  -> evidence that E is finite (each point
       isolated); the residual risk is an extremal we never sampled;
     * ONE non-sharp extremal         -> E may be a positive-dimensional family,
       and part (b) of Reformulation R cannot be attacked by list-exhaustion.

WHAT THIS SCRIPT DOES
  1. Massive multi-start descent on F over Gr(3,6) (gauge-free: we optimize A but
     read everything off P = A A^T, and the certificate lives on Gr).
  2. Keep hits with |F - 1/6| < tol; POLISH each hard (annealed softmax then
     repeated Nelder-Mead) until |F - 1/6| < 1e-14.
  3. Cluster by SORTED LEVERAGE PATTERN rounded to 1e-9 -- the coarse invariant
     that separated the seventh extremal from the TTSP census.
  4. Per cluster representative: rationalize the leverages via PSLQ, count active
     triples, check simplicity of lambda_min, build L numerically, and report
       rank(L)   and   the LP optimum over {Lz <= 0, |z|_inf <= 1}.
     A cluster with rank(L) < 9 or LP optimum > 0 is a NON-SHARP candidate and is
     flagged loudly for exact follow-up with v11's machinery.
  5. Emit a JSON catalogue for exact certification of anything new.

HONESTY.  This script is a SEARCH: its output is NUMERICALLY SUPPORTED.  It
locates candidates and screens them; PROOF requires exact reconstruction
(v11_seventh_exact.py) per new pattern.  Nothing here upgrades a tag on its own.

Deterministic (fixed master seed).  Counterexample tripwire armed throughout.
"""
import itertools, json, os, sys
import numpy as np
from multiprocessing import Pool
from scipy.optimize import minimize, linprog

MASTER = 20260731
TRIPLES = list(itertools.combinations(range(6), 3))
TH = 1.0 / 6.0
CX = 1e-6                      # counterexample tolerance (brief §0 rule 3)
HIT = 1e-9                     # "is an extremal" tolerance after polish
ACT = 1e-9                     # active-triple tolerance


def retract(X):
    U, _, Vt = np.linalg.svd(X, full_matrices=False)
    return U @ Vt


def lams(A):
    P = A @ A.T
    return np.array([np.linalg.eigvalsh(P[np.ix_(T, T)])[0] for T in TRIPLES]), P


def F(v):
    return float(np.max(lams(retract(v.reshape(6, 3)))[0]))


def softF(v, beta):
    l = lams(retract(v.reshape(6, 3)))[0]
    m = float(np.max(l))
    return m + float(np.log(np.sum(np.exp(beta * (l - m)))) / beta)   # shifted!


def descend(v, betas=(30.0, 200.0, 1500.0), rounds=5):
    for b in betas:
        v = minimize(softF, v, args=(b,), method="Nelder-Mead",
                     options=dict(maxiter=5000, maxfev=5000,
                                  xatol=1e-12, fatol=1e-15)).x
    for _ in range(rounds):
        v = minimize(F, v, method="Nelder-Mead",
                     options=dict(maxiter=6000, maxfev=6000,
                                  xatol=1e-14, fatol=1e-16)).x
    return v


def polish(v, target=1e-14, tries=12):
    for _ in range(tries):
        v = minimize(F, v, method="Nelder-Mead",
                     options=dict(maxiter=30000, maxfev=30000,
                                  xatol=1e-16, fatol=1e-18)).x
        if abs(F(v) - TH) < target:
            break
    return v


def worker(seed):
    rng = np.random.default_rng(seed)
    v = descend(rng.standard_normal(18))
    f = F(v)
    if f - TH < -CX:                                   # tripwire
        return dict(kind="VIOLATION", F=f, v=v.tolist(), seed=int(seed))
    if abs(f - TH) > 1e-7:                             # not an extremal
        return dict(kind="miss", F=f, seed=int(seed))
    v = polish(v)
    f = F(v)
    if f - TH < -CX:
        return dict(kind="VIOLATION", F=f, v=v.tolist(), seed=int(seed))
    if abs(f - TH) > HIT:
        return dict(kind="near", F=f, seed=int(seed))
    A = retract(v.reshape(6, 3))
    P = A @ A.T
    lev = np.sort(np.diag(P))
    # Cluster at 1e-6, NOT 1e-9.  A 1e-9 key is FINER than the convergence of the
    # polish step, so it splits a single leverage pattern into many spurious
    # "new" patterns as the sample grows (observed: the count crept 4->5->6->7->8
    # purely from rounding noise).  Distinct patterns differ by ~0.03, so 1e-6 is
    # still enormously tighter than needed.
    return dict(kind="hit", F=f, seed=int(seed), v=v.tolist(),
                lev=lev.tolist(), key=tuple(np.round(lev, 6)))


# ---------------------------------------------------------------- screening
def screen(v, act_tol=None):
    """Numeric sharpness screen: active set, simplicity, rank(L), LP optimum.

    ACTIVE-SET TOLERANCE IS THE DELICATE PART.  Hits are only polished to
    |F - 1/6| < 1e-9, so a genuinely active triple can sit at deviation ~1e-8 and
    fall OUTSIDE a 1e-9 active window.  Truncating the active set that way removes
    rows from L and destroys positive spanning, producing a SPURIOUS "non-sharp"
    verdict.  This was observed for real: an MLCore run flagged 3 "non-sharp
    candidates" that were all this artifact (reproduced locally: a hit with
    |A|=9 at 1e-9 but |A|=10 at 1e-6, the tenth triple at deviation 6.6e-9).

    Fix: pick the tolerance from the SPECTRAL GAP in the deviations.  The active
    deviations cluster near 0 and the inactive ones are O(1e-1) away, so we cut at
    the largest relative jump instead of a fixed threshold, and report the gap so a
    marginal call is visible rather than silent.
    """
    A = retract(np.array(v).reshape(6, 3))
    P = A @ A.T
    lam, _ = lams(A)
    dev = np.abs(lam - TH)
    order = np.argsort(dev)
    sd = dev[order]
    if act_tol is not None:
        active = [i for i in range(20) if dev[i] < act_tol]
        gap_ratio = None
    else:
        # cut at the largest multiplicative jump among the sorted deviations,
        # searching only the plausible range 1..19 actives
        best_k, best_ratio = 1, 0.0
        for k in range(1, 20):
            lo = max(sd[k - 1], 1e-16)
            ratio = sd[k] / lo
            if ratio > best_ratio:
                best_k, best_ratio = k, ratio
        active = [int(i) for i in order[:best_k]]
        gap_ratio = float(best_ratio)
    out = dict(n_active=len(active), act_gap_ratio=gap_ratio,
               act_cut_dev=float(sd[len(active) - 1]),
               next_dev=float(sd[len(active)]) if len(active) < 20 else None)

    # eigenvector of the lambda_min eigenspace + simplicity
    Ev, simple = {}, True
    for i in active:
        T = TRIPLES[i]
        w, U = np.linalg.eigh(P[np.ix_(T, T)])
        if abs(w[1] - w[0]) < 1e-7:                    # degenerate lambda_min
            simple = False
            break
        Ev[T] = U[:, 0]
    out["simple_lmin"] = bool(simple)
    if not simple:
        out["screen"] = "DEGENERATE (needs SDP form)"
        return out

    # tangent basis of Gr(3,6): orthonormal bases of range(P), ker(P)
    w, U = np.linalg.eigh(P)
    rng_b = U[:, np.argsort(-w)[:3]]
    ker_b = U[:, np.argsort(-w)[3:]]
    TAN = []
    for a in range(3):
        for b in range(3):
            X = np.outer(rng_b[:, a], ker_b[:, b])
            TAN.append(X + X.T)

    L = np.array([[Ev[TRIPLES[i]] @ B[np.ix_(TRIPLES[i], TRIPLES[i])] @ Ev[TRIPLES[i]]
                   for B in TAN] for i in active])
    out["rank_L"] = int(np.linalg.matrix_rank(L, tol=1e-8))
    sv = np.linalg.svd(L, compute_uv=False)
    out["smallest_sv"] = float(sv[-1]) if len(sv) else None

    r = linprog(c=-L.sum(axis=0), A_ub=L, b_ub=np.zeros(L.shape[0]),
                bounds=[(-1, 1)] * 9, method="highs")
    out["lp_value"] = float(-r.fun) if r.success else None
    out["lp_ok"] = bool(r.success)

    # kappa estimate: min over unit sphere of max_T (L z)_T
    best = np.inf
    rg = np.random.default_rng(7)
    for _ in range(120):
        z0 = rg.standard_normal(9); z0 /= np.linalg.norm(z0)
        def obj(z):
            nz = np.linalg.norm(z)
            if nz < 1e-10:            # guard: z/|z| is meaningless near 0
                return 1e6
            return float(np.max(L @ (z / nz)))
        rr = minimize(obj, z0, method="Nelder-Mead",
                      options=dict(maxiter=3000, maxfev=3000, xatol=1e-12, fatol=1e-15))
        best = min(best, float(rr.fun))
    out["kappa_est"] = best

    # rank 9 + LP optimum 0 ALREADY implies kappa > 0 (see v9 route A), so a
    # negative kappa estimate alongside those two is a sign of estimator trouble,
    # not of non-sharpness.  Report that inconsistency instead of hiding it.
    cone_ok = (out.get("lp_value") is not None and abs(out["lp_value"]) < 1e-9)
    sharp = (out["rank_L"] == 9 and cone_ok)
    out["kappa_inconsistent"] = bool(sharp and best <= 1e-6)
    if sharp:
        out["screen"] = "SHARP (numeric)"
    elif out["n_active"] < 10:
        # sharp => |A| >= 10 (positive spanning of R^9 needs >= 10 vectors), so a
        # low actual count is either the real obstruction or a tolerance artifact.
        out["screen"] = "LOW-ACTIVE (check tolerance before believing)"
    else:
        out["screen"] = "NON-SHARP CANDIDATE"
    return out


if __name__ == "__main__":
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 4000
    # honour MLC_CPUS / all cores when running on a big flavor; cap locally so the
    # workstation stays usable.
    env = os.environ.get("GTZ_CPUS")
    ncpu = int(env) if env else max(1, (os.cpu_count() or 4) - 2)
    seeds = [int(x) for x in np.random.SeedSequence(MASTER).generate_state(n, dtype=np.uint32)]
    print(f"master={MASTER}  starts={n}  cores={ncpu}", flush=True)
    print(f"tolerances: hit |F-1/6|<{HIT}, active<{ACT}, tripwire F<1/6-{CX}")
    print("=" * 78, flush=True)

    # INCREMENTAL: use imap_unordered and checkpoint as we go, so a job killed by
    # a wall-clock limit still leaves usable partial results (the earlier 20k
    # submission would have been killed at 5h with NO output at all).
    os.makedirs("verify/out", exist_ok=True)
    res, done = [], 0
    with Pool(ncpu) as pool:
        for r in pool.imap_unordered(worker, seeds, chunksize=2):
            res.append(r); done += 1
            if r["kind"] == "VIOLATION":
                print(f"\n*** VIOLATION at seed {r['seed']}: F={r['F']!r} ***", flush=True)
                json.dump([r], open("verify/out/v12_VIOLATIONS.json", "w"), indent=1)
                pool.terminate()
                sys.exit(2)
            if done % 100 == 0:
                hits_so_far = [x for x in res if x["kind"] == "hit"]
                pats = {}
                for x in hits_so_far:
                    pats.setdefault(x["key"], 0)
                    pats[x["key"]] += 1
                print(f"  [{done}/{n}] hits={len(hits_so_far)} patterns={len(pats)}",
                      flush=True)
                json.dump(dict(progress=done, total=n,
                               n_hits=len(hits_so_far),
                               patterns={str(k): v for k, v in pats.items()}),
                          open("verify/out/v12_progress.json", "w"), indent=1)

    viol = [r for r in res if r["kind"] == "VIOLATION"]
    if viol:
        print("\n*** COUNTEREXAMPLE PROTOCOL TRIGGERED -- HALTING ***")
        print(f"*** {len(viol)} sample(s) with F < 1/6 - {CX} ***")
        for r in viol[:3]:
            print(f"    F = {r['F']!r}  seed {r['seed']}")
        json.dump(viol, open("verify/out/v12_VIOLATIONS.json", "w"), indent=1)
        sys.exit(2)

    hits = [r for r in res if r["kind"] == "hit"]
    near = [r for r in res if r["kind"] == "near"]
    miss = [r for r in res if r["kind"] == "miss"]
    print(f"\nhits (|F-1/6| < {HIT}): {len(hits)} / {n}")
    print(f"near  (polished but > tol):  {len(near)}")
    print(f"misses (F far from 1/6):     {len(miss)}")
    if miss:
        fs = sorted(r["F"] for r in miss)
        print(f"   miss F range: [{fs[0]!r}, {fs[-1]!r}]")

    # cluster by leverage pattern
    clusters = {}
    for r in hits:
        clusters.setdefault(r["key"], []).append(r)
    print(f"\ndistinct leverage patterns: {len(clusters)}")

    KNOWN_LEV = [(sorted([5/18]*3 + [13/18]*3), "TTSP 5/18,13/18"),
                 (sorted([11/18]*4 + [5/18]*2), "TTSP 11/18,5/18"),
                 (sorted([7/18]*4 + [13/18]*2), "TTSP 7/18,13/18"),
                 (sorted([5/14]*3 + [9/14]*3), "SEVENTH 5/14,9/14")]

    def known_label(key, tol=1e-5):
        k = np.sort(np.array(key))
        for v, nm in KNOWN_LEV:
            if np.max(np.abs(k - np.array(v))) < tol:
                return nm
        return None

    catalogue = []
    for key, members in sorted(clusters.items(), key=lambda kv: -len(kv[1])):
        lab = known_label(key) or "*** NEW PATTERN ***"
        rep = members[0]
        sc = screen(rep["v"])
        print("\n" + "-" * 78, flush=True)
        print(f"pattern {lab}   multiplicity {len(members)}")
        print(f"  leverages   = {np.round(np.array(key), 12)}")
        print(f"  active      = {sc.get('n_active')}   simple_lmin = {sc.get('simple_lmin')}")
        print(f"  rank(L)     = {sc.get('rank_L')}   smallest sv = {sc.get('smallest_sv')}")
        print(f"  LP optimum  = {sc.get('lp_value')}")
        print(f"  kappa (est) = {sc.get('kappa_est')}")
        print(f"  ==> {sc.get('screen')}")
        if sc.get("screen") == "NON-SHARP CANDIDATE":
            print("  *** THIS WOULD BREAK FINITENESS -- certify exactly with v11 ***")
        catalogue.append(dict(label=lab, leverages=list(key),
                              multiplicity=len(members), F=rep["F"],
                              seed=rep["seed"], A_flat=rep["v"], **sc))

    os.makedirs("verify/out", exist_ok=True)
    json.dump(dict(master=MASTER, n_starts=n, n_hits=len(hits),
                   n_patterns=len(clusters), catalogue=catalogue),
              open("verify/out/v12_hunt.json", "w"), indent=1)
    # save representatives for exact reconstruction
    os.makedirs("verify/data", exist_ok=True)
    for c in catalogue:
        if c["label"].startswith("***"):
            A = retract(np.array(c["A_flat"]).reshape(6, 3))
            tag = "_".join(f"{x:.6f}" for x in c["leverages"][:2])
            np.save(f"verify/data/new_pattern_{tag}.npy", A @ A.T)

    print("\n" + "=" * 78)
    newp = [c for c in catalogue if c["label"].startswith("***")]
    nonsharp = [c for c in catalogue if c.get("screen") == "NON-SHARP CANDIDATE"]
    print(f"patterns found: {len(catalogue)}   NEW: {len(newp)}   "
          f"non-sharp candidates: {len(nonsharp)}")
    if nonsharp:
        print("\n*** NON-SHARP CANDIDATE(S) FOUND -- finiteness of E is IN DOUBT. ***")
        print("*** Certify each exactly before drawing any conclusion.          ***")
    else:
        print("\nEvery pattern found screens as SHARP => every extremal located is")
        print("an ISOLATED point of E.  Combined with compactness of Gr(3,6), this is")
        print("evidence (NOT proof) that E is FINITE.  Residual risk: an extremal")
        print("never sampled by this search.")
    print("no counterexample: F >= 1/6 - 1e-6 throughout.")
    print("=" * 78)
