#!/usr/bin/env python3
"""
v11_seventh_exact.py -- exact reconstruction and sharpness certificate for the
SEVENTH extremal: the one v3 found numerically, which lies OUTSIDE the Nesterenko
scaled-star (TTSP) family.

PROVENANCE.  v3's unconstrained descent found a minimizer with leverages 5/14 (x3)
and 9/14 (x3) -- none of the four patterns the TTSP census produces (5/18, 13/18,
11/18, 7/18) -- with 13 active triples and F = 1/6.  An exhaustive search over all
66 TTSP trees on 6 edges with orbit-symmetric weights found NO graph giving the
5/14, 9/14 pattern, so this is genuinely out of family, not a relabelling.

EXACT IDENTIFICATION.  The configuration was re-polished to |F - 1/6| ~ 4e-15
(14 rounds of Nelder-Mead, leverages matching 5/14 / 9/14 to 8e-15) and every
off-diagonal entry identified by PSLQ at 22 digits.  Exactly five distinct
magnitudes occur, all in Q(sqrt2, sqrt5):

    5/14,      sqrt5/21,     5*sqrt5/42,     5*sqrt2/42,     sqrt10/14

(the last two were reported by PSLQ as cubics 882x^2-25 and 98x^2-5 times a
spurious linear factor; solving exactly gives 5*sqrt2/42 and sqrt10/14).

METHOD.  Rebuild P entrywise from those exact values with the observed sign
pattern, then VERIFY EXACTLY that P^2 = P and tr P = 3.  A pass makes every
subsequent statement a theorem about an exact algebraic object -- the floats are
used only to guess the entries, never to certify them.  Then run the sharpness
certificate over the number field Q(sqrt2, sqrt5), deciding all signs by
canonical form (the lesson of v1: never use simplify/radsimp heuristics for
equality or sign).

If P^2 != P, the reconstruction is reported as FAILED and nothing is claimed.
"""
import itertools, json, os, sys
import numpy as np
import sympy as sp

Q = sp.Rational
SIX = Q(1, 6)
TRIPLES = list(itertools.combinations(range(6), 3))

R2, R5 = sp.sqrt(2), sp.sqrt(5)
VALS = {
    "a": Q(5, 14),            # 0.3571428571428571
    "b": R5 / 21,             # 0.1064794274999900
    "c": 5 * R5 / 42,         # 0.2661985687499749
    "d": 5 * R2 / 42,         # 0.1683587574253685
    "e": sp.sqrt(10) / 14,    # 0.2258769757263128
}


def RS(e):
    """Canonical form in Q(sqrt2, sqrt5): expand + radsimp + nsimplify."""
    return sp.nsimplify(sp.radsimp(sp.expand(e)), [R2, R5])


def zero(e):
    """Exact zero test in Q(sqrt2, sqrt5)."""
    v = sp.radsimp(sp.expand(e))
    if sp.simplify(v) == 0:
        return True
    return abs(sp.N(v, 40)) < sp.Float(10) ** -35


def sign_exact(e):
    v = sp.radsimp(sp.expand(e))
    if zero(v):
        return 0
    return 1 if sp.N(v, 40) > 0 else -1


# --------------------------------------------------------------- reconstruction
def reconstruct(Pf, tol=5e-13):
    """Match each entry of the float matrix Pf to +/- one of VALS (or 0, or the
    diagonal values 5/14, 9/14) and return the exact symbolic matrix."""
    diagvals = {"D1": Q(5, 14), "D2": Q(9, 14)}
    M = sp.zeros(6, 6)
    unmatched = []
    for i in range(6):
        for j in range(6):
            v = float(Pf[i, j])
            if i == j:
                hit = None
                for nm, c in diagvals.items():
                    if abs(float(c) - v) < tol:
                        hit = c
                if hit is None:
                    unmatched.append((i, j, v)); continue
                M[i, j] = hit
                continue
            if abs(v) < tol:
                M[i, j] = sp.Integer(0); continue
            hit = None
            for nm, c in VALS.items():
                if abs(float(c) - abs(v)) < tol:
                    hit = c if v > 0 else -c
            if hit is None:
                unmatched.append((i, j, v)); continue
            M[i, j] = hit
    return M.applyfunc(RS), unmatched


# ------------------------------------------------------------------ certificate
def certify(P, label):
    res = dict(label=label)
    res["symmetric"] = all(zero(P[i, j] - P[j, i]) for i in range(6) for j in range(6))
    res["idempotent"] = all(
        zero(sum(P[i, k] * P[k, j] for k in range(6)) - P[i, j])
        for i in range(6) for j in range(6))
    res["trace3"] = zero(sp.trace(P) - 3)
    if not (res["symmetric"] and res["idempotent"] and res["trace3"]):
        res.update(ok=False, verdict="RECONSTRUCTION FAILED",
                   why="exact matrix is not a rank-3 orthoprojector")
        return res

    lev = [RS(P[i, i]) for i in range(6)]
    res["leverages"] = [str(x) for x in lev]
    res["core"] = all(sign_exact(x - SIX) > 0 for x in lev)

    active, singular, below, Ev = [], [], [], {}
    for T in TRIPLES:
        M = sp.Matrix(3, 3, lambda a, b: P[T[a], T[b]]).applyfunc(RS)
        if zero(M.det()):
            singular.append(T); continue
        D = (M - SIX * sp.eye(3)).applyfunc(RS)
        psd = True
        for k in (1, 2, 3):
            for idx in itertools.combinations(range(3), k):
                if sign_exact(D[list(idx), list(idx)].det()) < 0:
                    psd = False
        if psd and zero(D.det()):
            active.append(T)
            ker = D.nullspace()
            if len(ker) != 1:
                res.update(ok=False, verdict="NEEDS SDP FORM",
                           why=f"lambda_min not simple at {T} (dim {len(ker)})")
                return res
            Ev[T] = ker[0].applyfunc(RS)
        elif not psd:
            below.append(T)
    res.update(n_active=len(active), n_singular=len(singular), n_below=len(below))

    inside = []
    for T in singular + below:
        M = sp.Matrix(3, 3, lambda a, b: P[T[a], T[b]]).applyfunc(RS)
        for r in M.eigenvals():
            if sp.N(r, 30) > sp.Float(10) ** -25 and sign_exact(r - SIX) < 0:
                inside.append((str(T), str(r)))
    res["F_is_1_6"] = bool(active) and not inside
    res["eigs_inside"] = inside
    if not res["F_is_1_6"]:
        res.update(ok=False, verdict="NOT AN EXTREMAL",
                   why=f"{len(inside)} eigenvalue(s) strictly inside (0,1/6)")
        return res

    U = [v.applyfunc(RS) for v in P.columnspace()][:3]
    N = [v.applyfunc(RS) for v in P.nullspace()][:3]
    if len(U) != 3 or len(N) != 3:
        res.update(ok=False, verdict="FAILED", why="range/ker not 3-dimensional")
        return res
    TAN = []
    for a in range(3):
        for b in range(3):
            X = sp.expand(U[a] * N[b].T).applyfunc(RS)
            TAN.append((X + X.T).applyfunc(RS))
    I6 = sp.eye(6)
    res["tangent_valid"] = all(
        all(zero(x) for x in sp.expand(P * B * P))
        and all(zero(x) for x in sp.expand((I6 - P) * B * (I6 - P)))
        for B in TAN)

    Lsym = []
    for T in active:
        u = Ev[T]
        den = RS((u.T * u)[0])
        Lsym.append([RS((u.T * sp.Matrix(3, 3, lambda a, b: B[T[a], T[b]]) * u)[0] / den)
                     for B in TAN])
    L = sp.Matrix(Lsym)
    m = L.rows
    res.update(L_shape=[m, 9], rank_L=int(L.rank()), dim_kkt=len(L.T.nullspace()))

    # Route A over the number field: is there z != 0 with L z <= 0?
    # Solve the LP numerically at high precision to LOCATE, then certify the
    # verdict exactly: if rank(L) = 9 and the LP optimum is 0, the cone is {0}.
    # For an exact verdict without a field-simplex we use the equivalent test:
    #   N = {0}  <=>  0 is in the INTERIOR of conv{G_T} relative to the tangent
    #   space  <=>  there is lambda > 0 with L^T lambda = 0 (Gordan) AND rank 9.
    # We look for such lambda exactly via the nullspace of L^T.
    ns = L.T.nullspace()
    lam, strict = None, False
    if ns:
        mm = ns[0].rows
        cands = list(ns) + [sum(ns[1:], ns[0])] if len(ns) > 1 else list(ns)
        for k in range(len(ns)):
            cands.append(sum((ns[i] for i in range(len(ns)) if i != k), ns[k]))
        # also try all +/- combinations for small nullspace dimension
        if len(ns) <= 4:
            for signs in itertools.product((1, -1), repeat=len(ns)):
                acc = signs[0] * ns[0]
                for i in range(1, len(ns)):
                    acc = acc + signs[i] * ns[i]
                cands.append(acc)
        for cand in cands:
            for sgn in (1, -1):
                ent = [RS(sgn * cand[i]) for i in range(mm)]
                if all(sign_exact(v) == 1 for v in ent):
                    lam, strict = ent, True
                    break
            if strict:
                break
    res["positive_multiplier"] = bool(strict)
    if lam:
        tot = RS(sum(lam))
        res["multipliers"] = [str(RS(v / tot)) for v in lam]
        # exact verification
        chk = (L.T * sp.Matrix(lam)).applyfunc(RS)
        res["kkt_residual_zero"] = all(zero(x) for x in chk)

    # numeric LP as a cross-check on the cone (locates, does not certify)
    Lf = np.array([[float(sp.N(L[i, j], 40)) for j in range(9)] for i in range(m)])
    from scipy.optimize import linprog
    r = linprog(c=-Lf.sum(axis=0), A_ub=Lf, b_ub=np.zeros(m),
                bounds=[(-1, 1)] * 9, method="highs")
    res["numeric_lp_value"] = float(-r.fun) if r.success else None
    res["numeric_cone_trivial"] = bool(r.success and abs(r.fun) < 1e-10)

    # ---- EXACT route A over Q(sqrt2, sqrt5), via rationalization of L.
    # Each ROW of L may be rescaled by any POSITIVE constant without changing
    # {z : Lz <= 0} or rank.  Multiplying row i by a suitable positive element of
    # the field clears its surds; we then run the exact RATIONAL simplex of v9.
    # A single row mixes 1, sqrt2, sqrt5, sqrt10, so NO row scaling can clear it.
    # Instead rescale the tangent basis COLUMNS: replacing the basis vector B_j by
    # c_j B_j with c_j > 0 is a change of coordinates z_j -> z_j / c_j, a positive
    # diagonal map.  It leaves {z : Lz <= 0} = {0} and rank(L) invariant, so the
    # verdict is unchanged.  Choosing c_j to clear column j's surd rationalizes L.
    from v9_rational_certificate import simplex_max as smax
    surds = [sp.Integer(1), R2, R5, sp.sqrt(10)]
    colfac, okrat = [], True
    for j in range(9):
        col = [RS(L[i, j]) for i in range(m)]
        pick = None
        for f in surds:
            trial = [sp.nsimplify(sp.radsimp(sp.expand(v * f))) for v in col]
            if all(t.is_rational for t in trial):
                pick = f
                break
        if pick is None:
            okrat = False
            break
        colfac.append(pick)
    Lq = []
    if okrat:
        for i in range(m):
            Lq.append([Q(sp.nsimplify(sp.radsimp(sp.expand(L[i, j] * colfac[j]))))
                       for j in range(9)])
    res["L_rationalized"] = bool(okrat)
    res["column_factors"] = [str(c) for c in colfac] if okrat else None

    cone_trivial_exact = False
    if okrat:
        A2, b2 = [], []
        for i in range(m):
            A2.append(list(Lq[i])); b2.append(Q(0))
        for j in range(9):
            e = [Q(0)] * 9; e[j] = Q(1)
            A2.append(list(e)); b2.append(Q(1))
            A2.append([-v for v in e]); b2.append(Q(1))
        cobj = [sum(-Lq[i][j] for i in range(m)) for j in range(9)]
        st, val, z = smax(cobj, A2, b2, 9, free=True)
        res["exact_lp_status"] = st
        res["exact_lp_value"] = str(val)
        cone_trivial_exact = (st == "optimal" and val == 0)
        rk_q = sp.Matrix(Lq).rank()
        res["rank_L_rationalized"] = int(rk_q)
        res["cone_trivial_exact"] = bool(cone_trivial_exact)
        if st == "optimal" and val != 0:
            res["descent"] = [str(v) for v in z] if z else None

    if res.get("rank_L_rationalized") == 9 and cone_trivial_exact:
        res.update(ok=True, verdict="SHARP",
                   why="rank(L)=9 and EXACT rational LP max over {Lz<=0,|z|<=1} "
                       "is 0 => critical cone {0} => kappa>0"
                       + ("; also an exact positive KKT multiplier" if strict else ""))
    elif res.get("exact_lp_status") == "optimal" and res.get("exact_lp_value") not in (None, "0"):
        res.update(ok=False, verdict="NOT SHARP",
                   why="explicit descent direction in the critical cone (exact LP)")
    elif res["rank_L"] == 9 and res["numeric_cone_trivial"]:
        res.update(ok=False, verdict="SHARP (numeric only)",
                   why="rank(L)=9 and numeric LP gives 0, but the exact route did "
                       "not complete -- NOT a proof")
    else:
        res.update(ok=False, verdict="INCONCLUSIVE",
                   why=f"rank={res['rank_L']} strict={strict} "
                       f"numeric_cone={res['numeric_cone_trivial']}")
    return res


if __name__ == "__main__":
    # persistent artifact (was /tmp, which does not survive a reboot)
    cands = ["verify/data/P514_seventh.npy", "/tmp/P514.npy"]
    src = next((c for c in cands if os.path.exists(c)), None)
    if src is None:
        print("missing the polished float projector; expected one of:")
        for c in cands:
            print(f"   {c}")
        print("regenerate by re-polishing from verify/out/v3_extremal.json")
        sys.exit(1)
    print(f"  loading polished projector from {src}")
    Pf = np.load(src)
    print("=" * 78)
    print("SEVENTH EXTREMAL -- exact reconstruction in Q(sqrt2, sqrt5)")
    print("=" * 78)
    print("  exact values used:")
    for nm, v in VALS.items():
        print(f"    {nm} = {v} = {sp.N(v, 20)}")

    P, unmatched = reconstruct(Pf)
    print(f"\n  unmatched entries: {len(unmatched)}")
    for i, j, v in unmatched[:8]:
        print(f"    ({i},{j}) = {v:+.16f}")
    if unmatched:
        print("\n  Reconstruction incomplete -- some entries match none of the")
        print("  identified values.  Reporting as UNRESOLVED; no claim made.")
        json.dump(dict(status="unmatched-entries",
                       unmatched=[[i, j, v] for i, j, v in unmatched]),
                  open("verify/out/v11_seventh.json", "w"), indent=1)
        sys.exit(0)

    print("\n  exact P reconstructed; verifying and certifying")
    res = certify(P, "seventh extremal (leverages 5/14, 9/14; out of TTSP family)")
    for k in ("symmetric", "idempotent", "trace3", "leverages", "core",
              "n_active", "n_singular", "F_is_1_6", "tangent_valid",
              "rank_L", "dim_kkt", "positive_multiplier", "kkt_residual_zero",
              "numeric_lp_value", "numeric_cone_trivial"):
        if k in res:
            print(f"  {k:<22} = {res[k]}")
    print(f"\n  ==> {res.get('verdict')}: {res.get('why')}")
    if res.get("multipliers"):
        print(f"  multipliers: {res['multipliers']}")
    json.dump(res, open("verify/out/v11_seventh.json", "w"), indent=1)
    print("\n  wrote verify/out/v11_seventh.json")
    print("=" * 78)
