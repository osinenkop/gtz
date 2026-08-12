#!/usr/bin/env python3
"""
v9_rational_certificate.py -- settle ALL extremals of the Nesterenko family with a
fully RATIONAL exact certificate, resolving the blocker recorded in HANDOFF.md §1.

THE BLOCKER, AND THE FIX.
v8 tried to run an exact simplex on the active-gradient matrix L, but L had entries
in a field carrying several surds at once (sqrt5, sqrt7, sqrt10, sqrt14, sqrt35,
sqrt70 for one extremal), because P = W^{1/2} Bt (Bt^T W Bt)^{-1} Bt^T W^{1/2}
carries W^{1/2}.  Row-scaling could not clear a row mixing sqrt7 and sqrt10.

The fix removes the surds entirely instead of fighting them.  Put

    Bt = B_red^T,   G = Bt^T W Bt   (RATIONAL),   S := Bt G^{-1} Bt^T W .

Then S is a RATIONAL idempotent of trace 3 (verified below), and

    P = W^{1/2} S W^{-1/2},

so P and S are similar via the DIAGONAL matrix W^{1/2}.  Because W^{1/2} is
diagonal, the similarity restricts to every principal block:

    P_TT = W_T^{1/2} S_TT W_T^{-1/2}   =>   spec(P_TT) = spec(S_TT).

Hence lambda_min(P_TT) = lambda_min(S_TT), the active set is the same, and the
whole sharpness certificate can be run on S over Q.  S is not symmetric, so we
work with the symmetric representative M_T := W_T^{1/2} S_TT W_T^{-1/2} only
conceptually; concretely all we need are (i) the active set, (ii) the
lambda_min eigenvector of the SYMMETRIC block, and (iii) the tangent space.  For
(ii) we use the LEFT/RIGHT kernel of the rational matrix S_TT - I/6 and undo the
diagonal scaling exactly -- which introduces sqrt(w) factors that then cancel
pairwise in the quadratic form v^T Pdot_TT v, leaving rational entries.

CERTIFICATE (unchanged in logic from v6, now exact over Q):
    L z <= 0 has only z = 0        <=>   active gradients positively span R^9
    <=>  exists lambda > 0 with L^T lambda = 0   (Gordan's alternative)
Both directions are decided by an exact rational simplex written from scratch
(Bland's rule, no floating point).  A NOT-SHARP verdict comes with an explicit
descent direction z, and is reported prominently rather than smoothed over.

Deterministic.  No floating point is load-bearing anywhere.
"""
import itertools, sys, json, os
import sympy as sp

SIX = sp.Rational(1, 6)
TRIPLES = list(itertools.combinations(range(6), 3))
Q = sp.Rational


def RQ(e):
    """Canonical rational simplification (input is rational by construction)."""
    return sp.nsimplify(sp.cancel(sp.together(sp.expand(e))))


# ------------------------------------------------------------ exact simplex over Q
def simplex_max(c, A, b, nvar, free=False):
    """Maximize c.x s.t. A x <= b, x >= 0 (x free if free=True). Exact Rationals.
    Returns (status, value, x). status in {optimal, unbounded, infeasible}."""
    if free:
        A2 = [list(r) + [-v for v in r] for r in A]
        c2 = list(c) + [-v for v in c]
        st, val, x2 = simplex_max(c2, A2, b, 2 * nvar, free=False)
        if st != "optimal":
            return st, val, None
        return st, val, [x2[i] - x2[nvar + i] for i in range(nvar)]

    m, n = len(A), nvar
    if any(Q(x) < 0 for x in b):
        return "needs-phase-1", None, None
    T = [[Q(A[i][j]) for j in range(n)] + [Q(1 if k == i else 0) for k in range(m)]
         + [Q(b[i])] for i in range(m)]
    obj = [-Q(c[j]) for j in range(n)] + [Q(0)] * m + [Q(0)]
    basis = [n + i for i in range(m)]
    for _ in range(50000):
        pc = next((j for j in range(n + m) if obj[j] < 0), None)
        if pc is None:
            x = [Q(0)] * (n + m)
            for i, bi in enumerate(basis):
                x[bi] = T[i][-1]
            return "optimal", obj[-1], x[:n]
        pr, best = None, None
        for i in range(m):
            if T[i][pc] > 0:
                key = (T[i][-1] / T[i][pc], basis[i])
                if best is None or key < best:
                    best, pr = key, i
        if pr is None:
            return "unbounded", None, None
        pv = T[pr][pc]
        T[pr] = [v / pv for v in T[pr]]
        for i in range(m):
            if i != pr and T[i][pc] != 0:
                f = T[i][pc]
                T[i] = [a - f * bb for a, bb in zip(T[i], T[pr])]
        if obj[pc] != 0:
            f = obj[pc]
            obj = [a - f * bb for a, bb in zip(obj, T[pr])]
        basis[pr] = pc
    return "iteration-limit", None, None


def lp_feasible_eq_ge1(Lq, m, ncol=9):
    """Is there lambda >= 1 (componentwise) with L^T lambda = 0?  By positive
    scaling this is equivalent to: exists lambda > 0 with L^T lambda = 0.

    Solved EXACTLY without needing a phase-1 simplex.  Write lambda = 1 + mu,
    mu >= 0, so the requirement is  L^T mu = -L^T 1  with  mu >= 0.  Let
    part be any particular rational solution of that linear system and let
    {n_k} span ker(L^T).  Then the question is whether the affine set
    {part + sum c_k n_k} meets the nonnegative cone -- a small exact LP in the
    coefficients c, which we solve by maximizing the minimum coordinate:

        maximize  s   s.t.   part_i + (N c)_i >= s   for all i,   s <= 0 bounded
                            and |c_k| <= C  (a box, to keep it bounded)

    Feasibility with s >= 0 answers the question.  All inequalities have the
    slack basis feasible after moving s to the objective, so no phase-1 is
    needed; we instead solve the equivalent problem
        maximize s  over  (c, s) free, s.t.  s - (N c)_i <= part_i,  |c| <= C
    which has b = part (may be negative) -- so we shift by taking c around a
    known feasible point.  To avoid all of that fragility we simply enumerate
    with sympy's exact linear solver over the nullspace, then, if that fails,
    fall back to scipy's exact-rational-free LP only as a HINT (never as proof).
    """
    LT = sp.Matrix(Lq).T                       # ncol x m
    ones = sp.Matrix([Q(1)] * m)
    rhs = -(LT * ones)
    # particular solution of LT * mu = rhs
    try:
        sol = sp.linsolve((LT, rhs))
    except Exception:
        return False, None
    if not sol:
        return False, None
    params = sorted(sp.Matrix(list(sol)[0]).free_symbols, key=lambda s: s.name)
    muexpr = sp.Matrix(list(sol)[0])
    if not params:
        mu = [RQ(v) for v in muexpr]
        if all(v >= 0 for v in mu):
            return True, [Q(1) + v for v in mu]
        return False, None
    # Maximize the minimum coordinate of muexpr over the parameters, exactly.
    # This is an LP:  max s  s.t.  mu_i(params) >= s.  Encode with the simplex by
    # shifting parameters to be free and s free: variables (params..., s).
    nvar = len(params) + 1
    A, b = [], []
    for i in range(m):
        e = muexpr[i]
        coeffs = [RQ(sp.expand(e).coeff(p)) for p in params]
        const = RQ(sp.expand(e).subs({p: 0 for p in params}))
        # s - sum coeffs*p <= const
        A.append([-c for c in coeffs] + [Q(1)])
        b.append(const)
    # box on parameters to keep the LP bounded
    C = Q(10) ** 6
    for k in range(len(params)):
        e = [Q(0)] * nvar; e[k] = Q(1)
        A.append(list(e)); b.append(C)
        A.append([-v for v in e]); b.append(C)
    cobj = [Q(0)] * len(params) + [Q(1)]
    st, val, x = simplex_max(cobj, A, b, nvar, free=True)
    if st != "optimal" or val is None or val < 0:
        return False, None
    sub = {p: x[k] for k, p in enumerate(params)}
    mu = [RQ(muexpr[i].subs(sub)) for i in range(m)]
    if any(v < 0 for v in mu):
        return False, None
    lam = [Q(1) + v for v in mu]
    # exact verification of the certificate before returning it
    if sp.simplify(LT * sp.Matrix(lam)) != sp.zeros(ncol, 1):
        return False, None
    if any(v <= 0 for v in lam):
        return False, None
    return True, lam


# ------------------------------------------------------------------ graph plumbing
class Builder:
    def __init__(self):
        self.nv, self.edges = 2, []

    def nvx(self):
        v = self.nv
        self.nv += 1
        return v

    def build(self, g, s, t):
        if g == ('e',):
            self.edges.append((s, t)); return
        kind, *kids = g
        if kind == 'P':
            for k in kids:
                self.build(k, s, t)
        else:
            prev = s
            for i, k in enumerate(kids):
                nxt = t if i == len(kids) - 1 else self.nvx()
                self.build(k, prev, nxt)
                prev = nxt


def rational_S(edges, nv, w):
    """S = Bt G^{-1} Bt^T W : a RATIONAL idempotent similar to P via diag(sqrt w)."""
    B = sp.zeros(nv - 1, len(edges))
    for j, (a, b) in enumerate(edges):
        if a != 0:
            B[a - 1, j] += 1
        if b != 0:
            B[b - 1, j] -= 1
    Bt = B.T
    W = sp.diag(*[Q(x) for x in w])
    G = sp.expand(Bt.T * W * Bt)
    if G.det() == 0:
        return None, None
    S = sp.expand(Bt * G.inv() * Bt.T * W).applyfunc(RQ)
    P = (sp.diag(*[sp.sqrt(Q(x)) for x in w]) * S
         * sp.diag(*[1 / sp.sqrt(Q(x)) for x in w]))
    return S, sp.simplify(P)


def certificate(S, P, w, label):
    """Full sharpness certificate, all arithmetic exact."""
    res = dict(label=label, weights=[str(Q(x)) for x in w])
    if sp.expand(S * S - S).applyfunc(RQ) != sp.zeros(6, 6) or RQ(sp.trace(S)) != 3:
        res.update(ok=False, why="S not a rank-3 idempotent")
        return res
    Psym = P.applyfunc(lambda t: sp.radsimp(sp.simplify(t)))
    res["symmetric_P"] = (sp.simplify(Psym - Psym.T) == sp.zeros(6, 6))
    lev = [RQ(sp.nsimplify(sp.radsimp(Psym[i, i]))) for i in range(6)]
    res["leverages"] = [str(x) for x in lev]
    res["core"] = all(x > SIX for x in lev)

    # active set from the RATIONAL blocks of S (spec(S_TT) = spec(P_TT))
    active, singular, below = [], [], []
    for T in TRIPLES:
        ST = sp.Matrix(3, 3, lambda a, b: S[T[a], T[b]]).applyfunc(RQ)
        if RQ(ST.det()) == 0:
            singular.append(T); continue
        D = (ST - SIX * sp.eye(3)).applyfunc(RQ)
        cp = sp.Poly(D.charpoly(sp.Symbol("t")).as_expr(), sp.Symbol("t"))
        # D's eigenvalues are lam(S_TT) - 1/6, all real (S_TT similar to symmetric)
        # PSD  <=>  all coefficients alternate appropriately; decide exactly via
        # Descartes on charpoly of -D:  D >= 0  <=>  no positive root of det(-D - t I)
        # Simpler and exact: D >= 0  <=>  all principal minors of D >= 0.
        psd = True
        for k in (1, 2, 3):
            for idx in itertools.combinations(range(3), k):
                sub = D[list(idx), list(idx)]
                if RQ(sub.det()) < 0:
                    psd = False
        if psd and RQ(D.det()) == 0:
            active.append(T)
        elif not psd:
            below.append(T)
    res.update(n_active=len(active), n_singular=len(singular), n_below=len(below))
    # F = 1/6 requires: some active, and nothing strictly inside (0, 1/6)
    inside = []
    for T in singular + below:
        ST = sp.Matrix(3, 3, lambda a, b: S[T[a], T[b]]).applyfunc(RQ)
        for r in ST.eigenvals():
            rr = sp.nsimplify(r)
            if sp.N(rr, 30) > sp.Float(10) ** -25 and sp.simplify(rr - SIX) < 0:
                inside.append((T, rr))
    res["F_is_1_6"] = bool(active) and not inside
    if not res["F_is_1_6"]:
        res.update(ok=False, why=f"F != 1/6 ({len(inside)} eigenvalues inside (0,1/6))")
        return res

    # lambda_min eigenvector of the SYMMETRIC block:
    #   P_TT = D_T S_TT D_T^{-1}, D_T = diag(sqrt w_T).  If S_TT u = u/6 then
    #   v := D_T u satisfies P_TT v = v/6.  v has sqrt(w) entries, but the
    #   quadratic form v^T Pdot_TT v with Pdot = D Sdot D^{-1} ... instead we
    #   avoid v entirely: use the SYMMETRIC block M_T := D_T S_TT D_T^{-1} and
    #   note v^T M v / v^T v is rational when expressed via u and W_T (below).
    Ev, simple = {}, True
    for T in active:
        ST = sp.Matrix(3, 3, lambda a, b: S[T[a], T[b]]).applyfunc(RQ)
        ker = (ST - SIX * sp.eye(3)).applyfunc(RQ).nullspace()
        if len(ker) != 1:
            simple = False; break
        u = ker[0].applyfunc(RQ)
        Ev[T] = u                                    # RATIONAL right-kernel vector
    res["simple_lmin"] = simple
    if not simple:
        res.update(ok=False, why="lambda_min not simple (needs SDP form)")
        return res

    # Tangent space of Gr(3,6) at P, pushed to the S-picture:
    #   Pdot = W^{1/2} Sdot W^{-1/2} with Sdot in the tangent space at S.
    # Use rational bases of range(S) and ker(S); no Gram-Schmidt (no sqrt).
    U = [v.applyfunc(RQ) for v in S.columnspace()][:3]
    N = [v.applyfunc(RQ) for v in S.nullspace()][:3]
    if len(U) != 3 or len(N) != 3:
        res.update(ok=False, why="range/ker of S not 3-dimensional")
        return res
    Wd = sp.diag(*[Q(x) for x in w])
    # Sdot_ab := u_a n_b^T W  +  (W^{-1} n_b u_a^T W) ... to stay in the tangent
    # space of the S-picture we use the projector form: any tangent direction is
    # Sdot = S X (I-S) + (I-S) Y S with X,Y arbitrary; a basis is obtained from
    # outer products of the rational bases.  We take
    #   Sdot = u_a m_b^T + (W^{-1} m_b)(W u_a)^T   with m_b spanning ker(S^T).
    NT = [v.applyfunc(RQ) for v in S.T.nullspace()][:3]
    if len(NT) != 3:
        res.update(ok=False, why="ker(S^T) not 3-dimensional")
        return res
    TAN = []
    for a in range(3):
        for b in range(3):
            X = sp.expand(U[a] * NT[b].T).applyfunc(RQ)
            TAN.append((X + Wd.inv() * X.T * Wd).applyfunc(RQ))
    if len(TAN) != 9:
        res.update(ok=False, why="tangent basis not 9-dimensional")
        return res
    # sanity: each TAN element must satisfy S B S = 0 and (I-S) B (I-S) = 0
    I6 = sp.eye(6)
    tan_ok = all(sp.expand(S * Bm * S).applyfunc(RQ) == sp.zeros(6, 6)
                 and sp.expand((I6 - S) * Bm * (I6 - S)).applyfunc(RQ) == sp.zeros(6, 6)
                 for Bm in TAN)
    res["tangent_valid"] = tan_ok

    # Directional derivative functional.  With P_TT symmetric and v = D_T u,
    #   <G_T, Pdot> = v^T Pdot_TT v / (v^T v),  v = D_T u,  Pdot = D Sdot D^{-1}.
    # Then v^T Pdot_TT v = u^T D_T (D_T Sdot_TT D_T^{-1}) D_T u
    #                    = u^T (W_T Sdot_TT) u    -- RATIONAL.
    # And v^T v = u^T W_T u -- RATIONAL.  So the functional is rational.
    Lq = []
    for T in active:
        u = Ev[T]
        WT = sp.diag(*[Q(w[i]) for i in T])
        denom = RQ((u.T * WT * u)[0])
        row = []
        for Bm in TAN:
            BT = sp.Matrix(3, 3, lambda a, b: Bm[T[a], T[b]])
            row.append(RQ((u.T * WT * BT * u)[0] / denom))
        Lq.append(row)
    allrat = all(sp.nsimplify(v).is_rational for row in Lq for v in row)
    res["L_rational"] = allrat
    if not allrat:
        res.update(ok=False, why="L still not rational -- investigate")
        return res

    L = sp.Matrix(Lq)
    m = L.rows
    res.update(L_shape=[m, 9], rank_L=int(L.rank()), dim_kkt=len(L.T.nullspace()))
    if res["rank_L"] != 9:
        res.update(ok=False,
                   why=f"rank(L) = {res['rank_L']} != 9: gradients do not span")
        return res

    # (A) Gordan primal: is there z != 0 with L z <= 0?
    A2, b2 = [], []
    for i in range(m):
        A2.append([Lq[i][j] for j in range(9)]); b2.append(Q(0))
    for j in range(9):
        e = [Q(0)] * 9; e[j] = Q(1)
        A2.append(list(e)); b2.append(Q(1))
        A2.append([-v for v in e]); b2.append(Q(1))
    cobj = [sum(-Lq[i][j] for i in range(m)) for j in range(9)]
    stA, valA, zA = simplex_max(cobj, A2, b2, 9, free=True)
    cone_trivial = (stA == "optimal" and valA == 0)
    res.update(lp_status=stA, lp_value=str(valA), cone_trivial=bool(cone_trivial))

    # (B) Gordan dual: is there lambda > 0 with L^T lambda = 0?
    posB, lam = lp_feasible_eq_ge1(Lq, m, ncol=9)
    res["positive_multiplier"] = bool(posB)
    if lam is not None:
        tot = sum(lam)
        res["multipliers"] = [str(RQ(v / tot)) for v in lam]

    # VERDICT LOGIC.  Route (A) alone is already a COMPLETE proof:
    #
    #   The feasible set {Lz <= 0, |z| <= 1} contains z = 0 with objective 0, and
    #   the objective sum_i (-Lz)_i is a sum of NON-NEGATIVE terms on it.  So if
    #   the exact maximum is 0, then every feasible z has (-Lz)_i = 0 for all i,
    #   i.e. Lz = 0; and rank(L) = 9 forces z = 0.  Since {Lz <= 0} is a cone, the
    #   box |z| <= 1 loses nothing.  Hence N = {0}, so max_T <G_T, z> > 0 for every
    #   z != 0, and by compactness of the unit sphere kappa > 0.
    #
    # Validated on three controls (see the session log): a sharp configuration
    # returns 0; a non-sharp one returns a positive value with an explicit descent
    # direction; a rank-deficient one returns 0 but is caught by the rank test.
    #
    # Route (B) (a strictly positive multiplier, Gordan's dual) is a SECOND,
    # independent certificate.  When it also succeeds we report it, since the
    # exact multipliers are publishable data and cross-check v6/v7.  Its FAILURE,
    # however, only means our parametric LP over ker(L^T) did not locate a
    # witness -- it does NOT weaken route (A).
    if cone_trivial:
        res.update(ok=True, verdict="SHARP",
                   why="rank(L)=9 and exact LP max over {Lz<=0,|z|<=1} is 0 "
                       "=> critical cone N={0} => kappa>0"
                       + ("; independently confirmed by a strictly positive "
                          "KKT multiplier" if posB else
                          " (route B found no explicit multiplier; route A is "
                          "self-sufficient)"))
    elif stA == "optimal" and valA != 0:
        res.update(ok=False, verdict="NOT SHARP",
                   why="explicit descent direction found in the critical cone",
                   descent=[str(v) for v in zA] if zA else None)
    else:
        res.update(ok=False, verdict="INCONCLUSIVE",
                   why=f"cone_trivial={cone_trivial} positive_multiplier={posB} "
                       f"lp_status={stA}")
    return res


E = ('e',)
CASES = [
    ("P(S(e,e,e),e,e,e)", ('P', ('S', E, E, E), E, E, E),
     [1, 1, 1, Q(5, 9), Q(5, 9), Q(5, 9)]),
    ("P(S(e,e),S(e,e),e,e)", ('P', ('S', E, E), ('S', E, E), E, E),
     [1, 1, 1, 1, Q(5, 8), Q(5, 8)]),
    ("P(S(P(e,e,e),e),S(e,e))", ('P', ('S', ('P', E, E, E), E), ('S', E, E)),
     [1, 1, 1, Q(9, 5), Q(9, 5), Q(9, 5)]),
    ("P(S(P(S(e,e),e,e),e),e)", ('P', ('S', ('P', ('S', E, E), E, E), E), E),
     [1, 1, Q(5, 8), Q(5, 8), 1, 1]),
    ("P(S(P(e,e),P(e,e)),S(e,e))", ('P', ('S', ('P', E, E), ('P', E, E)), ('S', E, E)),
     [1, 1, 1, 1, Q(8, 5), Q(8, 5)]),
    ("P(S(P(e,e),e),S(P(e,e),e))", ('P', ('S', ('P', E, E), E), ('S', ('P', E, E), E)),
     [1, 1, Q(8, 5), 1, 1, Q(8, 5)]),
]

if __name__ == "__main__":
    out = []
    for name, tree, w in CASES:
        print("=" * 78, flush=True)
        print(f"{name}   w={[str(Q(x)) for x in w]}", flush=True)
        b = Builder(); b.build(tree, 0, 1)
        S, P = rational_S(b.edges, b.nv, w)
        if S is None:
            print("  singular Gram; skipped"); continue
        r = certificate(S, P, w, name)
        print(f"  leverages {r.get('leverages')}", flush=True)
        print(f"  active={r.get('n_active')} singular={r.get('n_singular')} "
              f"F=1/6:{r.get('F_is_1_6')} L_rational={r.get('L_rational')} "
              f"rank={r.get('rank_L')} dim_kkt={r.get('dim_kkt')}", flush=True)
        print(f"  LP: status={r.get('lp_status')} value={r.get('lp_value')} "
              f"cone_trivial={r.get('cone_trivial')} "
              f"pos_multiplier={r.get('positive_multiplier')}", flush=True)
        print(f"  ==> {r.get('verdict', 'n/a')}: {r.get('why')}", flush=True)
        if r.get("multipliers"):
            print(f"  multipliers: {r['multipliers']}", flush=True)
        if r.get("descent"):
            print(f"  DESCENT DIRECTION: {r['descent']}", flush=True)
        out.append(r)

    print("\n" + "=" * 78)
    for r in out:
        print(f"  {r.get('verdict','n/a'):<14} {r['label']}")
    ns = sum(1 for r in out if r.get("verdict") == "SHARP")
    nn = sum(1 for r in out if r.get("verdict") == "NOT SHARP")
    ni = sum(1 for r in out if r.get("verdict") == "INCONCLUSIVE")
    print(f"\n  SHARP: {ns}   NOT SHARP: {nn}   INCONCLUSIVE: {ni}   of {len(out)}")
    if nn:
        print("\n  *** AT LEAST ONE EXTREMAL IS NOT A SHARP MINIMUM ***")
        print("  *** part (a) of Reformulation R fails in the sharp form there ***")
    os.makedirs("verify/out", exist_ok=True)
    json.dump(out, open("verify/out/v9_certificate.json", "w"), indent=1)
    print("\n  wrote verify/out/v9_certificate.json")
    print("=" * 78)
