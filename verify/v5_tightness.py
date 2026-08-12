#!/usr/bin/env python3
"""
v5_tightness.py -- test the "no-uniform-gap" picture directly, now that v4 has
established F(A_0) = 1/6 EXACTLY at the Nesterenko extremal.

Three experiments:

 (E1) LOCAL: start AT the exact Nesterenko extremal A_0 (F = 1/6 exactly) and
      run a descent.  If GTZ(6,3) is true, no descent can go below 1/6.  Any
      reproducible descent below 1/6 - 1e-6 is a counterexample (brief §0 rule 3).
      This is the single sharpest test available: it starts on the boundary.

 (E2) PERTURBATION SWEEP: random perturbations of A_0 at geometrically
      decreasing radius r, measuring min F over each shell.  If A_0 is a strict
      local min (Reformulation R part (a), currently NUMERICALLY SUPPORTED),
      then min over each shell should be >= 1/6, approaching 1/6 as r -> 0, and
      the deficit should scale like r^2 (a genuine minimum), not like r
      (a saddle / first-order descent direction).

 (E3) GLOBAL FLOOR: many random starts, recording how close the unconstrained
      minimum gets to 1/6 from above.  v2 reached only ~0.1849; the question is
      whether generic descent finds the 1/6 floor at all, or whether the
      extremals are hard to reach by descent (which matters a lot for how a
      branch-and-bound / Track C search must be seeded).

Deterministic.  Everything here is NUMERICALLY SUPPORTED; nothing is proved.
"""
import itertools, sys, os, json
import numpy as np
import sympy as sp
from scipy.optimize import minimize
from multiprocessing import Pool

MASTER = 20260801
TRIPLES = list(itertools.combinations(range(6), 3))
THRESH = 1.0 / 6.0
CX = 1e-6


# ---------- exact A_0 from v4, converted to float64
def nesterenko_A():
    Bred = sp.Matrix([[-1, -1, -1, 0, 0, -1],
                      [0, 0, 0, -1, 1, 0],
                      [0, 0, 0, 0, -1, 1]])
    w = [sp.Integer(1)] * 3 + [sp.Rational(9, 5)] * 3
    Mt = sp.diag(*[sp.sqrt(x) for x in w]) * Bred.T
    P0 = sp.simplify(Mt * (Mt.T * Mt).inv() * Mt.T)
    Pf = np.array([[float(sp.N(P0[i, j], 40)) for j in range(6)] for i in range(6)])
    ev, V = np.linalg.eigh(Pf)
    A = V[:, np.argsort(-ev)[:3]]                  # top-3 eigenvectors: A A^T = P0
    return A, Pf


def retract(X):
    U, _, Vt = np.linalg.svd(X, full_matrices=False)
    return U @ Vt


def lams(A):
    P = A @ A.T
    return np.array([np.linalg.eigvalsh(P[np.ix_(T, T)])[0] for T in TRIPLES])


def F(v):
    return float(np.max(lams(retract(v.reshape(6, 3)))))


def softF(v, beta):
    l = lams(retract(v.reshape(6, 3)))
    m = float(np.max(l))
    return m + float(np.log(np.sum(np.exp(beta * (l - m)))) / beta)


def descend(v, rounds=6, betas=(60.0, 400.0, 2500.0)):
    for b in betas:
        v = minimize(softF, v, args=(b,), method="Nelder-Mead",
                     options=dict(maxiter=4000, maxfev=4000, xatol=1e-13, fatol=1e-16)).x
    for _ in range(rounds):
        v = minimize(F, v, method="Nelder-Mead",
                     options=dict(maxiter=5000, maxfev=5000, xatol=1e-14, fatol=1e-16)).x
    return v, float(F(v))


def _e2(args):
    seed, r, v0 = args
    rng = np.random.default_rng(seed)
    d = rng.standard_normal(18)
    d /= np.linalg.norm(d)
    return float(F(np.asarray(v0) + r * d))


def _e3(seed):
    rng = np.random.default_rng(seed)
    v, f = descend(rng.standard_normal(18), rounds=3, betas=(40.0, 300.0, 2000.0))
    return f, v.tolist(), int(seed)


if __name__ == "__main__":
    ncpu = max(1, min(20, os.cpu_count() - 2))
    A0, P0f = nesterenko_A()
    v0 = A0.reshape(-1).copy()
    F0 = F(v0)
    print("=" * 74)
    print(f"A_0 = Nesterenko extremal;  F(A_0) = {F0!r}")
    print(f"1/6 = {THRESH!r};   F(A_0) - 1/6 = {F0 - THRESH:+.3e}")
    print(f"(v4 proved F(A_0) = 1/6 EXACTLY; float error above is roundoff)")
    print(f"leverages = {np.round(np.sort(np.diag(A0 @ A0.T)), 12)}")
    print(f"           exact = 5/18={5/18:.12f}, 13/18={13/18:.12f}")
    print("=" * 74)

    # ---------------- E1: descend FROM the extremal
    print("\n[E1] descent starting AT the exact extremal (F = 1/6)")
    v1, f1 = descend(v0.copy(), rounds=8)
    print(f"     after descent: F = {f1!r}")
    print(f"     F - 1/6 = {f1 - THRESH:+.6e}")
    e1_viol = f1 < THRESH - CX
    print(f"     descended below 1/6 - 1e-6?  {e1_viol}   <-- must be False")

    # also try many randomized descents from the extremal (kick then descend)
    print("\n[E1b] 200 kicked descents from the extremal (kick 1e-3, then descend)")
    rng = np.random.default_rng(MASTER)
    worst, nviol = np.inf, 0
    for _ in range(200):
        d = rng.standard_normal(18); d /= np.linalg.norm(d)
        _, f = descend(v0 + 1e-3 * d, rounds=2, betas=(200.0, 2000.0))
        worst = min(worst, f)
        if f < THRESH - CX:
            nviol += 1
    print(f"     worst F over 200 kicked descents = {worst!r}")
    print(f"     F - 1/6 = {worst - THRESH:+.6e};  violations = {nviol}   <-- must be 0")

    # ---------------- E2: perturbation shells
    print("\n[E2] perturbation shells around A_0 (is it a strict local min?)")
    print(f"     {'radius':>10}  {'min F on shell':>22}  {'minF - 1/6':>14}  {'ratio':>9}")
    shells, prev = [], None
    ss = np.random.SeedSequence(MASTER + 1)
    for k, r in enumerate([10.0 ** -e for e in range(1, 8)]):
        seeds = [int(x) for x in ss.spawn(1)[0].generate_state(600, dtype=np.uint32)]
        with Pool(ncpu) as pool:
            vals = pool.map(_e2, [(s, r, v0.tolist()) for s in seeds], chunksize=20)
        mn = float(np.min(vals))
        dev = mn - THRESH
        ratio = (dev / prev) if prev not in (None, 0) else float("nan")
        print(f"     {r:>10.1e}  {mn!r:>22}  {dev:>+14.3e}  {ratio:>9.3f}")
        shells.append(dict(r=r, minF=mn, dev=dev))
        prev = dev
    print("     interpretation: dev ~ r^2 (ratio ~ 0.01 per decade) => strict min;")
    print("                     dev ~ r   (ratio ~ 0.1)             => first-order descent exists;")
    print("                     dev < 0                            => COUNTEREXAMPLE.")
    neg = [s for s in shells if s["dev"] < -CX]
    print(f"     shells with min F < 1/6 - 1e-6: {len(neg)}   <-- must be 0")

    # ---------------- E3: how close does generic descent get to 1/6?
    n3 = int(sys.argv[1]) if len(sys.argv) > 1 else 600
    print(f"\n[E3] {n3} generic random-start descents: does descent find the 1/6 floor?")
    seeds = [int(x) for x in np.random.SeedSequence(MASTER + 2).generate_state(n3, dtype=np.uint32)]
    with Pool(ncpu) as pool:
        out = pool.map(_e3, seeds, chunksize=4)
    out.sort(key=lambda z: z[0])
    fs = np.array([z[0] for z in out])
    print(f"     min      F = {fs[0]!r}   (F - 1/6 = {fs[0] - THRESH:+.4e})")
    print(f"     5th pct  F = {np.percentile(fs, 5)!r}")
    print(f"     median   F = {np.median(fs)!r}")
    print(f"     # within 1e-4 of 1/6 : {int(np.sum(fs < THRESH + 1e-4))}")
    print(f"     # within 1e-2 of 1/6 : {int(np.sum(fs < THRESH + 1e-2))}")
    viol3 = int(np.sum(fs < THRESH - CX))
    print(f"     # below 1/6 - 1e-6   : {viol3}   <-- must be 0")

    os.makedirs("verify/out", exist_ok=True)
    with open("verify/out/v5_tightness.json", "w") as fh:
        json.dump(dict(master=MASTER, F_A0=F0, E1_after_descent=f1,
                       E1b_worst=worst, E1b_violations=nviol,
                       E2_shells=shells, E3_n=n3, E3_min=float(fs[0]),
                       E3_median=float(np.median(fs)),
                       E3_within_1e4=int(np.sum(fs < THRESH + 1e-4)),
                       E3_violations=viol3,
                       total_violations=int(e1_viol) + nviol + len(neg) + viol3), fh, indent=1)
    print("\nwrote verify/out/v5_tightness.json")

    total = int(e1_viol) + nviol + len(neg) + viol3
    print("=" * 74)
    if total:
        print(f"*** COUNTEREXAMPLE PROTOCOL TRIGGERED: {total} violation(s). HALT. ***")
        sys.exit(2)
    print("No violation anywhere: F >= 1/6 - 1e-6 in every experiment,")
    print("including descents launched exactly AT the tight extremal.")
    print("Status: NUMERICALLY SUPPORTED (not proved).")
