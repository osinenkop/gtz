#!/usr/bin/env python3
"""
v10_seventh.py -- exact identification and sharpness certificate for the extremal
that v3 found NUMERICALLY and that lies OUTSIDE the Nesterenko scaled-star family.

WHAT WAS FOUND.  v3's unconstrained descent (1500 starts) returned a minimizer with
  * leverages exactly 5/14 (x3) and 9/14 (x3)  -- NOT 5/18, 13/18, 11/18, 7/18,
    the four patterns produced by the TTSP scaled-star census (v7/v9),
  * F = 1/6 exactly (v3's PSLQ gave minpoly 36x^2-12x+1 = (6x-1)^2),
  * THIRTEEN active triples (vs 10 or 12 in the census),
  * lambda_min SIMPLE on every active block (so the nullspace form of the
    certificate applies; no SDP needed).

A search over all 66 TTSP trees on 6 edges (with orbit-symmetric weights) found NO
graph producing the 5/14, 9/14 leverage pattern.  So this is a genuinely different
configuration, not a relabelling of a census member.  It matters because the
"every extremal" scope of proofs/sharp-cone-at-extremal.md §4b is stated relative
to the census; a configuration outside the census widens what must be certified.

STRATEGY.  Rationalize the numerical projector exactly, then certify:
  1. take v3's float P, round every entry to a rational with small denominator
     (the leverages are 5/14 and 9/14, so denominators divide 14 or 28);
  2. VERIFY EXACTLY that the rounded matrix is a genuine rank-3 orthoprojector
     (P^2 = P, tr P = 3) -- if rounding broke it, the candidate is rejected, so a
     PASS here is a proof about an exact object, not about a float;
  3. run the v9 certificate (active set, simplicity, tangent space, exact rational
     LP / Gordan) on that exact P.

If step 2 fails, we report that honestly and fall back to lattice reduction
(PSLQ per entry) rather than fudging the rounding.

Deterministic; no floating point is load-bearing in the verdict.
"""
import itertools, json, os, sys
import numpy as np
import sympy as sp

Q = sp.Rational
SIX = Q(1, 6)
TRIPLES = list(itertools.combinations(range(6), 3))


def RQ(e):
    return sp.nsimplify(sp.cancel(sp.together(sp.expand(e))))


# reuse the exact simplex + LP helpers from v9
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from v9_rational_certificate import simplex_max, lp_feasible_eq_ge1  # noqa: E402


def rationalize(Pf, max_den=2520):
    """Round each entry to a rational with denominator dividing max_den."""
    M = sp.zeros(6, 6)
    for i in range(6):
        for j in range(6):
            M[i, j] = sp.nsimplify(sp.Rational(Pf[i, j]).limit_denominator(max_den))
    # symmetrize exactly
    return ((M + M.T) / 2).applyfunc(RQ)


def certify(P, label):
    res = dict(label=label)
    res["is_projector"] = (sp.expand(P * P - P).applyfunc(RQ) == sp.zeros(6, 6))
    res["trace3"] = (RQ(sp.trace(P)) == 3)
    res["symmetric"] = (sp.simplify(P - P.T) == sp.zeros(6, 6))
    if not (res["is_projector"] and res["trace3"] and res["symmetric"]):
        res.update(ok=False, verdict="REJECTED",
                   why="rationalized matrix is not an exact rank-3 orthoprojector")
        return res

    lev = [RQ(P[i, i]) for i in range(6)]
    res["leverages"] = [str(x) for x in lev]
    res["core"] = all(x > SIX for x in lev)

    active, singular, below, Ev = [], [], [], {}
    for T in TRIPLES:
        M = sp.Matrix(3, 3, lambda a, b: P[T[a], T[b]]).applyfunc(RQ)
        if RQ(M.det()) == 0:
            singular.append(T); continue
        D = (M - SIX * sp.eye(3)).applyfunc(RQ)
        psd = all(RQ(D[list(idx), list(idx)].det()) >= 0
                  for k in (1, 2, 3) for idx in itertools.combinations(range(3), k))
        if psd and RQ(D.det()) == 0:
            active.append(T)
            ker = D.nullspace()
            if len(ker) != 1:
                res.update(ok=False, verdict="NEEDS SDP FORM",
                           why=f"lambda_min not simple at {T}")
                return res
            Ev[T] = ker[0].applyfunc(RQ)
        elif not psd:
            below.append(T)
    res.update(n_active=len(active), n_singular=len(singular), n_below=len(below))

    inside = []
    for T in singular + below:
        M = sp.Matrix(3, 3, lambda a, b: P[T[a], T[b]]).applyfunc(RQ)
        for r in M.eigenvals():
            rr = sp.nsimplify(r)
            if sp.N(rr, 30) > sp.Float(10) ** -25 and sp.simplify(rr - SIX) < 0:
                inside.append((str(T), str(rr)))
    res["F_is_1_6"] = bool(active) and not inside
    res["eigs_inside_0_16"] = inside
    if not res["F_is_1_6"]:
        res.update(ok=False, verdict="NOT AN EXTREMAL",
                   why=f"F != 1/6; {len(inside)} eigenvalue(s) strictly inside (0,1/6)")
        return res

    # tangent space of Gr(3,6): rational bases, no Gram-Schmidt (P is rational here)
    U = [v.applyfunc(RQ) for v in P.columnspace()][:3]
    N = [v.applyfunc(RQ) for v in P.nullspace()][:3]
    if len(U) != 3 or len(N) != 3:
        res.update(ok=False, verdict="REJECTED", why="range/ker not 3-dimensional")
        return res
    TAN = []
    for a in range(3):
        for b in range(3):
            X = sp.expand(U[a] * N[b].T).applyfunc(RQ)
            TAN.append((X + X.T).applyfunc(RQ))
    I6 = sp.eye(6)
    res["tangent_valid"] = all(
        sp.expand(P * B * P).applyfunc(RQ) == sp.zeros(6, 6)
        and sp.expand((I6 - P) * B * (I6 - P)).applyfunc(RQ) == sp.zeros(6, 6)
        for B in TAN)

    Lq = []
    for T in active:
        u = Ev[T]
        den = RQ((u.T * u)[0])
        Lq.append([RQ((u.T * sp.Matrix(3, 3, lambda a, b: B[T[a], T[b]]) * u)[0] / den)
                   for B in TAN])
    res["L_rational"] = all(sp.nsimplify(v).is_rational for row in Lq for v in row)
    if not res["L_rational"]:
        res.update(ok=False, verdict="REJECTED", why="L not rational")
        return res

    L = sp.Matrix(Lq)
    m = L.rows
    res.update(L_shape=[m, 9], rank_L=int(L.rank()), dim_kkt=len(L.T.nullspace()))

    # exact primal LP (Gordan route A) -- self-sufficient given rank(L) = 9
    A2, b2 = [], []
    for i in range(m):
        A2.append([Lq[i][j] for j in range(9)]); b2.append(Q(0))
    for j in range(9):
        e = [Q(0)] * 9; e[j] = Q(1)
        A2.append(list(e)); b2.append(Q(1))
        A2.append([-v for v in e]); b2.append(Q(1))
    cobj = [sum(-Lq[i][j] for i in range(m)) for j in range(9)]
    st, val, z = simplex_max(cobj, A2, b2, 9, free=True)
    res.update(lp_status=st, lp_value=str(val))
    cone_trivial = (st == "optimal" and val == 0)
    res["cone_trivial"] = bool(cone_trivial)

    posB, lam = lp_feasible_eq_ge1(Lq, m, ncol=9)
    res["positive_multiplier"] = bool(posB)
    if lam:
        tot = sum(lam)
        res["multipliers"] = [str(RQ(v / tot)) for v in lam]

    if res["rank_L"] == 9 and cone_trivial:
        res.update(ok=True, verdict="SHARP",
                   why="rank(L)=9 and exact LP max is 0 => critical cone {0} => kappa>0")
    elif st == "optimal" and val != 0:
        res.update(ok=False, verdict="NOT SHARP",
                   why="explicit descent direction in the critical cone",
                   descent=[str(v) for v in z] if z else None)
    else:
        res.update(ok=False, verdict="INCONCLUSIVE",
                   why=f"rank={res['rank_L']} cone_trivial={cone_trivial} status={st}")
    return res


if __name__ == "__main__":
    d = json.load(open("verify/out/v3_extremal.json"))
    A = np.array(d["A"])
    Pf = A @ A.T
    print("=" * 78)
    print("v3's out-of-family extremal: leverages 5/14 (x3), 9/14 (x3)")
    print(f"  float F_min = {d['Fmin']!r}   F - 1/6 = {d['Fmin'] - 1/6:+.3e}")
    print(f"  float leverages = {np.round(np.sort(np.diag(Pf)), 12)}")
    print("=" * 78)

    best = None
    for den in (14, 28, 42, 56, 84, 126, 252, 504, 2520):
        P = rationalize(Pf, den)
        ok = (sp.expand(P * P - P).applyfunc(RQ) == sp.zeros(6, 6)
              and RQ(sp.trace(P)) == 3)
        print(f"  limit_denominator({den:>4}): exact projector? {ok}")
        if ok:
            best = P
            break

    if best is None:
        print("\n  Rounding did not yield an exact projector at any tried denominator.")
        print("  The configuration is real (v3 found it to 1e-10) but its exact")
        print("  algebraic form is not a small-denominator rational matrix.")
        print("  NEXT STEP: PSLQ each entry against Q(sqrt d) for small d, or")
        print("  reconstruct it from its 13-triple active-set equations directly.")
        print("  Reporting as UNRESOLVED rather than guessing.")
        json.dump(dict(status="unresolved-rationalization"),
                  open("verify/out/v10_seventh.json", "w"), indent=1)
        sys.exit(0)

    print("\n  exact P found; running the certificate")
    r = certify(best, "v3 out-of-family extremal (leverages 5/14, 9/14)")
    for k in ("leverages", "core", "n_active", "n_singular", "F_is_1_6",
              "tangent_valid", "L_rational", "rank_L", "dim_kkt",
              "lp_status", "lp_value", "cone_trivial", "positive_multiplier"):
        if k in r:
            print(f"  {k:<20} = {r[k]}")
    print(f"\n  ==> {r.get('verdict')}: {r.get('why')}")
    if r.get("multipliers"):
        print(f"  multipliers: {r['multipliers']}")
    if r.get("descent"):
        print(f"  DESCENT: {r['descent']}")
    json.dump(r, open("verify/out/v10_seventh.json", "w"), indent=1)
    print("\n  wrote verify/out/v10_seventh.json")
    print("=" * 78)
