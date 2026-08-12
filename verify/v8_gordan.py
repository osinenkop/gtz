#!/usr/bin/env python3
"""
v8_gordan.py -- settle the three extremals that v7 left UNDECIDED, using an exact
LP / Gordan alternative instead of the witness-guessing heuristic v7 used.

WHY v7 WAS INCONCLUSIVE, NOT NEGATIVE.  For each extremal v7 computes the 10-or-12
active gradients as rows of L (shape m x 9) and looks for lambda >= 0 with
L^T lambda = 0 and every lambda_T > 0.  When dim ker L^T = 1 the nullspace is a
single ray and reading off signs settles it.  Three extremals have
dim ker L^T = 3, and there v7 only tried the basis vectors and a few sums; failing
to find a strictly positive witness among those proves NOTHING.

THE CORRECT TEST.  Sharpness at P_0 is exactly
    N := {z : L z <= 0} = {0},
i.e. the active gradients POSITIVELY span R^9.  Two exact routes, both used here
so they cross-check:

  (A) Gordan's theorem.  Exactly one of the following holds:
        (i)  exists z != 0 with L z <= 0        [a non-trivial critical cone]
        (ii) exists lambda > 0 with L^T lambda = 0   [positive spanning]
      We decide (i) directly as an exact rational LP feasibility problem:
      minimize 0 s.t. L z <= 0, and check whether z = 0 is the only solution by
      maximizing sum_T (-L z)_T subject to L z <= 0 and a box -1 <= z <= 1.
      If the optimum is 0 the cone is trivial; if positive we have an explicit
      descent direction and the extremal is NOT sharp.

  (B) The relative interior of the nullspace.  Search ker L^T for a strictly
      positive vector by exact LP (Fourier-Motzkin free): solve
      L^T lambda = 0, lambda >= 1 (componentwise), which is feasible iff a
      strictly positive multiplier vector exists (scaling is free).

Both are solved in EXACT RATIONAL ARITHMETIC with a from-scratch simplex over
sympy Rationals -- no floating point anywhere, so the verdicts are PROVED either
way.  If an extremal turns out NOT sharp, that is reported prominently: it would
mean part (a) of Reformulation R is FALSE as stated for that extremal, which is a
substantive finding and must not be smoothed over.
"""
import itertools, sys, json, os
import sympy as sp

SIX = sp.Rational(1, 6)
TRIPLES = list(itertools.combinations(range(6), 3))


def R(e):
    return sp.radsimp(sp.simplify(sp.expand(e)))


# ----------------------------------------------------------------- exact simplex
def simplex_max(c, A, b, nvar, free=False):
    """Maximize c.x s.t. A x <= b, x >= 0 (or x free if free=True), exact Rationals.
    Returns (status, value, x).  status in {'optimal','unbounded','infeasible'}.
    Implemented via the standard tableau method with Bland's rule (anti-cycling).
    If free=True, substitute x = x+ - x- and solve the nonnegative problem."""
    if free:
        A2 = [row + [-v for v in row] for row in A]
        c2 = list(c) + [-v for v in c]
        st, val, x2 = simplex_max(c2, A2, b, 2 * nvar, free=False)
        if st != "optimal":
            return st, val, None
        x = [x2[i] - x2[nvar + i] for i in range(nvar)]
        return st, val, x

    m, n = len(A), nvar
    # tableau: [A | I | b] with objective row
    T = [[sp.Rational(A[i][j]) for j in range(n)] +
         [sp.Rational(1 if k == i else 0) for k in range(m)] +
         [sp.Rational(b[i])] for i in range(m)]
    obj = [-sp.Rational(c[j]) for j in range(n)] + [sp.Rational(0)] * m + [sp.Rational(0)]
    basis = [n + i for i in range(m)]

    # all b >= 0 in our uses, so the slack basis is feasible
    if any(sp.Rational(x) < 0 for x in b):
        return "infeasible", None, None

    for _ in range(20000):
        # Bland: smallest index with negative objective coefficient
        piv_c = None
        for j in range(n + m):
            if obj[j] < 0:
                piv_c = j
                break
        if piv_c is None:
            x = [sp.Rational(0)] * (n + m)
            for i, bi in enumerate(basis):
                x[bi] = T[i][-1]
            return "optimal", obj[-1], x[:n]
        # ratio test, Bland tie-break on smallest basis index
        piv_r, best = None, None
        for i in range(m):
            if T[i][piv_c] > 0:
                ratio = T[i][-1] / T[i][piv_c]
                key = (ratio, basis[i])
                if best is None or key < best:
                    best, piv_r = key, i
        if piv_r is None:
            return "unbounded", None, None
        # pivot
        pv = T[piv_r][piv_c]
        T[piv_r] = [v / pv for v in T[piv_r]]
        for i in range(m):
            if i != piv_r and T[i][piv_c] != 0:
                f = T[i][piv_c]
                T[i] = [a - f * bb for a, bb in zip(T[i], T[piv_r])]
        if obj[piv_c] != 0:
            f = obj[piv_c]
            obj = [a - f * bb for a, bb in zip(obj, T[piv_r])]
        basis[piv_r] = piv_c
    return "iteration-limit", None, None


# ------------------------------------------------------- rebuild an extremal
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
                self.build(k, prev, nxt); prev = nxt


def projector(edges, nv, w):
    B = sp.zeros(nv - 1, len(edges))
    for j, (a, b) in enumerate(edges):
        if a != 0:
            B[a - 1, j] += 1
        if b != 0:
            B[b - 1, j] -= 1
    Mt = sp.diag(*[sp.sqrt(x) for x in w]) * B.T
    return sp.simplify(Mt * (Mt.T * Mt).inv() * Mt.T).applyfunc(R)


def build_L(P):
    """Return (active, L) with L the m x 9 matrix of active gradients."""
    active, Ev = [], {}
    for T in TRIPLES:
        M = sp.Matrix(3, 3, lambda a, b: P[T[a], T[b]]).applyfunc(R)
        if R(M.det()) == 0:
            continue
        D = (M - SIX * sp.eye(3)).applyfunc(R)
        if R(D.det()) != 0:
            continue
        ker = D.nullspace()
        if len(ker) != 1:
            return None, None
        v = ker[0].applyfunc(R)
        Ev[T] = (v / sp.sqrt(R((v.T * v)[0]))).applyfunc(R)
        active.append(T)

    def gs(vs):
        o = []
        for v in vs:
            u = v
            for q in o:
                u = (u - R((q.T * v)[0]) * q).applyfunc(R)
            o.append((u / sp.sqrt(R((u.T * u)[0]))).applyfunc(R))
        return o

    U = gs([v.applyfunc(R) for v in P.columnspace()][:3])
    N = gs([v.applyfunc(R) for v in P.nullspace()][:3])
    TAN = [(U[a] * N[b].T + N[b] * U[a].T).applyfunc(R) for a in range(3) for b in range(3)]
    L = sp.Matrix([[R((Ev[T].T * sp.Matrix(3, 3, lambda a, b: B[T[a], T[b]]) * Ev[T])[0])
                    for B in TAN] for T in active])
    return active, L


E = ('e',)
CASES = [
    # (name, expression tree, weights)  -- the three v7 left undecided, plus one
    # already-SHARP case as a positive control.
    ("P(S(e,e,e),e,e,e)  [control, v7 says SHARP]",
     ('P', ('S', E, E, E), E, E, E),
     [1, 1, 1, sp.Rational(5, 9), sp.Rational(5, 9), sp.Rational(5, 9)]),
    ("P(S(P(S(e,e),e,e),e),e)",
     ('P', ('S', ('P', ('S', E, E), E, E), E), E),
     [1, 1, sp.Rational(5, 8), sp.Rational(5, 8), 1, 1]),
    ("P(S(P(e,e),P(e,e)),S(e,e))",
     ('P', ('S', ('P', E, E), ('P', E, E)), ('S', E, E)),
     [1, 1, 1, 1, sp.Rational(8, 5), sp.Rational(8, 5)]),
    ("P(S(P(e,e),e),S(P(e,e),e))",
     ('P', ('S', ('P', E, E), E), ('S', ('P', E, E), E)),
     [1, 1, sp.Rational(8, 5), 1, 1, sp.Rational(8, 5)]),
]

if __name__ == "__main__":
    out = []
    for name, tree, w in CASES:
        print("=" * 78)
        print(name)
        print(f"  weights {[str(x) for x in w]}")
        b = Builder()
        b.build(tree, 0, 1)
        P = projector(b.edges, b.nv, w)
        active, L = build_L(P)
        if L is None:
            print("  lambda_min not simple -> needs SDP form; skipped")
            continue
        m = L.rows
        print(f"  active = {m},  rank(L) = {L.rank()},  dim ker L^T = {len(L.T.nullspace())}")

        Lq = [[sp.nsimplify(L[i, j]) for j in range(9)] for i in range(m)]
        allrat = all(v.is_rational for row in Lq for v in row)
        print(f"  all entries rational? {allrat}")

        if not allrat:
            # L lives in a real quadratic field Q(sqrt d) (the weights bring in one
            # surd).  Two exact escapes, tried in order:
            #  (1) each ROW may be rescaled by any POSITIVE constant without
            #      changing {z : Lz<=0} or the existence of a strictly positive
            #      multiplier.  Divide each row by the positive square root of the
            #      surd it carries when that lands the whole row in Q.
            #  (2) otherwise, embed Q(sqrt d) -> Q by picking a rational
            #      approximation is NOT exact, so instead reduce over Q(sqrt d) by
            #      writing each entry as a + b sqrt(d) and using the exact order on
            #      that field via sp.sign (which is exact for algebraic numbers).
            surds = set()
            for row in Lq:
                for v in row:
                    for at in sp.preorder_traversal(sp.radsimp(v)):
                        if at.is_Pow and at.exp == sp.Rational(1, 2) and at.base.is_Rational:
                            surds.add(sp.nsimplify(at.base))
            print(f"  surds present: {sorted(surds, key=str)}")
            fixed, okrows = [], True
            for row in Lq:
                cand = row
                for s_ in list(surds) + [1]:
                    trial = [sp.nsimplify(sp.radsimp(sp.expand(v / sp.sqrt(s_))))
                             for v in row]
                    if all(t.is_rational for t in trial):
                        cand = trial
                        break
                if not all(sp.nsimplify(x).is_rational for x in cand):
                    okrows = False
                    break
                fixed.append(cand)
            if okrows:
                Lq, allrat = fixed, True
                print("  -> each row divided by a POSITIVE surd; L now rational "
                      "(both LP questions are invariant under positive row scaling)")

        if not allrat:
            # The active gradients are only defined up to a POSITIVE scale (each
            # v_T is a unit eigenvector; rescaling row i by c_i > 0 changes neither
            # {z : Lz <= 0} nor the existence of a strictly positive multiplier --
            # it just reweights lambda_i by 1/c_i).  So we may clear radicals
            # ROW BY ROW with a positive factor and land in Q, exactly.
            newL, ok = [], True
            for i in range(m):
                row = Lq[i]
                # positive common factor: gcd-free approach -- multiply by the
                # product of the distinct surds appearing, then verify rationality
                fac = sp.Integer(1)
                for v in row:
                    if v == 0:
                        continue
                    d = sp.denom(sp.radsimp(v))
                    n_ = sp.numer(sp.radsimp(v))
                    for expr in (d, n_):
                        for at in sp.preorder_traversal(expr):
                            if at.is_Pow and at.exp == sp.Rational(1, 2):
                                cand = sp.sqrt(at.base)
                                if sp.simplify(fac % 1) == 0:
                                    pass
                                if not (fac / cand).is_rational:
                                    fac = sp.simplify(fac * cand)
                scaled = [sp.nsimplify(sp.radsimp(sp.expand(v * fac))) for v in row]
                if not all(s.is_rational for s in scaled):
                    # fallback: scale by 1/max|entry| times a big rational approx is
                    # NOT exact, so instead scale by the row's own norm-squared
                    nrm2 = sp.nsimplify(sp.radsimp(sum(v * v for v in row)))
                    trial = [sp.nsimplify(sp.radsimp(sp.expand(v * sp.sqrt(nrm2))))
                             for v in row]
                    if all(t.is_rational for t in trial):
                        scaled = trial
                    else:
                        ok = False
                        break
                newL.append(scaled)
            if ok:
                Lq = newL
                allrat = True
                print("  -> rationalized row-by-row with POSITIVE factors "
                      "(sign structure and both LP answers are unchanged)")
            else:
                print("  -> could not rationalize exactly; skipping this case")
                out.append(dict(name=name, weights=[str(x) for x in w], m=m,
                                verdict="SKIPPED (irrational L)"))
                continue

        # ---- (B) strictly positive multiplier: L^T lam = 0, lam >= 1
        # equality as two inequalities; variables lam_i = 1 + mu_i, mu >= 0
        A, bb = [], []
        for j in range(9):
            row = [Lq[i][j] for i in range(m)]
            rhs = -sum(Lq[i][j] for i in range(m))
            A.append(row); bb.append(rhs)
            A.append([-v for v in row]); bb.append(-rhs)
        stB, valB, xB = simplex_max([0] * m, A, bb, m)
        posB = (stB == "optimal")
        print(f"  (B) exists lambda > 0 with L^T lambda = 0 : {posB}   [{stB}]")

        # ---- (A) Gordan: is there z != 0 with L z <= 0 ?
        # maximize sum_i (-Lz)_i  s.t.  Lz <= 0,  -1 <= z <= 1
        A2, b2 = [], []
        for i in range(m):
            A2.append([Lq[i][j] for j in range(9)]); b2.append(sp.Rational(0))
        for j in range(9):
            e = [sp.Rational(0)] * 9; e[j] = sp.Rational(1)
            A2.append(list(e)); b2.append(sp.Rational(1))
            A2.append([-v for v in e]); b2.append(sp.Rational(1))
        cobj = [sum(-Lq[i][j] for i in range(m)) for j in range(9)]
        stA, valA, zA = simplex_max(cobj, A2, b2, 9, free=True)
        print(f"  (A) max sum(-Lz) over {{Lz<=0, |z|<=1}} = {valA}   [{stA}]")
        coneA = (stA == "optimal" and valA == 0)
        print(f"      critical cone trivial (N = {{0}})? {coneA}")

        verdict = "SHARP" if (posB and coneA) else (
            "NOT SHARP" if stA == "optimal" and valA != 0 else "INCONCLUSIVE")
        print(f"  ==> VERDICT: {verdict}")
        if verdict == "NOT SHARP":
            print(f"      explicit descent direction z = {zA}")
        out.append(dict(name=name, weights=[str(x) for x in w], m=m,
                        rank=int(L.rank()), dim_kkt=len(L.T.nullspace()),
                        gordan_positive_multiplier=bool(posB),
                        cone_trivial=bool(coneA), verdict=verdict,
                        lp_value=str(valA)))

    print("\n" + "=" * 78)
    for o in out:
        print(f"  {o['verdict']:<12} {o['name']}")
    sharp = sum(1 for o in out if o["verdict"] == "SHARP")
    notsharp = sum(1 for o in out if o["verdict"] == "NOT SHARP")
    inc = sum(1 for o in out if o["verdict"] == "INCONCLUSIVE")
    print(f"\n  SHARP: {sharp}   NOT SHARP: {notsharp}   INCONCLUSIVE: {inc}")
    if notsharp:
        print("\n  *** AT LEAST ONE EXTREMAL IS NOT A SHARP MINIMUM. ***")
        print("  *** Part (a) of Reformulation R does NOT hold in the sharp form ***")
        print("  *** at every extremal.  Report this to Pavel prominently.        ***")
    os.makedirs("verify/out", exist_ok=True)
    json.dump(out, open("verify/out/v8_gordan.json", "w"), indent=1)
    print("\n  wrote verify/out/v8_gordan.json")
    print("=" * 78)
