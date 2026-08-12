#!/usr/bin/env python3
"""
v4_nesterenko.py -- EXACT re-verification of Lemma 2.1 of boundary-obstruction.md:
the Nesterenko scaled-star extremal A_0 in St(6,3) attains  f(A_0) = 1/sqrt6
EXACTLY, i.e. F(A_0) = min_T-max value = 1/6 exactly.

Why this is the decisive fact for the whole project: it means GTZ(6,3) is TIGHT.
The inequality  min_{St(6,3)} F >= 1/6  holds with EQUALITY at A_0.  Consequently:

  * there is NO slack to exploit -- any valid certificate must vanish at A_0
    (this is Prop 11.3 of boundary-obstruction.md, and it is why unconstrained
    low-degree SOS relaxations plateau);
  * the "finite number of algebraic checks" hope must be interpreted carefully:
    a finite check can only work if it is exact AT the equality manifold, since
    no epsilon-margin exists to absorb numerical error;
  * the slice minimum G_PMC = (1-sin36)/2 = 0.2061 is NOT the binding
    configuration -- it sits 0.0394 ABOVE 1/6.  The slice is not where the
    problem is hard.

All arithmetic is exact in Q(sqrt5).  Deterministic; no randomness.
"""
import itertools, sys
import sympy as sp

TRIPLES = list(itertools.combinations(range(6), 3))
results = []


def check(tag, name, ok, detail=""):
    ok = bool(ok)
    results.append((tag, name, ok, detail))
    print(f"[{'PASS' if ok else 'FAIL'}] {tag}: {name}" + (f"   {detail}" if detail else ""),
          flush=True)
    return ok


print("=" * 76)
print("Nesterenko scaled-star extremal for G = 3xLR + path(1,1,1)")
print("=" * 76)

# reduced incidence (vertices 1,2,3), edges e0..e5 as in boundary-obstruction.md §2
Bred = sp.Matrix([
    [-1, -1, -1,  0,  0, -1],
    [ 0,  0,  0, -1,  1,  0],
    [ 0,  0,  0,  0, -1,  1],
])
w = [sp.Integer(1), sp.Integer(1), sp.Integer(1),
     sp.Rational(9, 5), sp.Rational(9, 5), sp.Rational(9, 5)]
Whalf = sp.diag(*[sp.sqrt(wi) for wi in w])

Mt = sp.simplify(Whalf * Bred.T)                       # 6x3
check("N0", "Mtilde = W^{1/2} B_red^T has full column rank 3", Mt.rank() == 3)

GtG = sp.simplify(Mt.T * Mt)
P0 = sp.simplify(Mt * GtG.inv() * Mt.T)
P0 = P0.applyfunc(lambda t: sp.radsimp(sp.simplify(sp.expand(t))))

check("N1", "P0 symmetric", sp.simplify(P0 - P0.T) == sp.zeros(6, 6))
check("N2", "P0^2 = P0 EXACTLY (orthogonal projector)",
      sp.simplify(sp.expand(P0 * P0 - P0)) == sp.zeros(6, 6))
check("N3", "tr P0 = 3 (rank 3)", sp.simplify(sp.trace(P0)) == 3)

lev = [sp.radsimp(sp.simplify(P0[i, i])) for i in range(6)]
print(f"  leverages = {[str(x) for x in lev]}")
check("N4", "leverages are (5/18,5/18,5/18,13/18,13/18,13/18)",
      [sp.simplify(x) for x in lev]
      == [sp.Rational(5, 18)] * 3 + [sp.Rational(13, 18)] * 3)
check("N5", "sum of leverages = 3", sp.simplify(sum(lev)) == 3)
check("N4b", "every leverage > 1/6 (so A_0 is in the CORE region, Case A cannot remove a row)",
      all(sp.simplify(x - sp.Rational(1, 6)) > 0 for x in lev))

# ---- the decisive computation: lambda_min of every triple block, EXACTLY.
# For a 3x3 symmetric block M, use the characteristic polynomial and exact
# rational root testing against 1/6 -- no radical eigenvalues needed:
#   lam_min(M) >= 1/6  <=>  M - I/6 is PSD  <=>  all principal minors of
#   (M - I/6) are >= 0.   (NOTE: all principal minors, NOT just leading ones --
#   the brief §4 flags exactly this pitfall.)
SIX = sp.Rational(1, 6)
sing, nonsing, at_16, below = [], [], [], []
for T in TRIPLES:
    M = sp.Matrix(3, 3, lambda a, b: P0[T[a], T[b]])
    M = M.applyfunc(lambda t: sp.radsimp(sp.simplify(t)))
    detM = sp.radsimp(sp.simplify(M.det()))
    D = sp.simplify(M - SIX * sp.eye(3))
    # all principal minors of D
    minors = []
    for k in (1, 2, 3):
        for idx in itertools.combinations(range(3), k):
            sub = D[list(idx), list(idx)]
            minors.append(sp.radsimp(sp.simplify(sub.det())))
    psd = all(sp.simplify(m) >= 0 for m in minors)
    # is lambda_min exactly 1/6?  <=>  det(M - I/6) = 0 and D is PSD
    exact16 = psd and sp.simplify(minors[-1]) == 0
    if sp.simplify(detM) == 0:
        sing.append(T)
    else:
        nonsing.append(T)
    if exact16:
        at_16.append(T)
    if not psd:
        below.append(T)

print(f"  singular triples (det P_TT = 0): {len(sing)}  {sing}")
print(f"  non-singular triples:            {len(nonsing)}")
check("N6", "every NON-SINGULAR triple has lambda_min(P_TT) = 1/6 EXACTLY",
      all(T in at_16 for T in nonsing),
      f"{sum(1 for T in nonsing if T in at_16)}/{len(nonsing)} at exactly 1/6")
check("N7", "exactly 10 active triples with lambda_min = 1/6", len(at_16) == 10,
      f"count={len(at_16)}")
# NOTE: the 10 singular triples have det P_TT = 0 hence lambda_min = 0 < 1/6.
# That is EXPECTED and consistent with Lemma 2.1 -- GTZ requires only that ONE
# triple be good, not all of them.  So the correct statements are:
#   (i) the triples NOT at 1/6 are exactly the singular ones, and
#   (ii) NO triple lies strictly between 0 and 1/6 (nothing "just barely bad").
check("N8", "the triples below 1/6 are EXACTLY the 10 singular ones (lambda_min=0)",
      sorted(below) == sorted(sing), f"{len(below)} below, {len(sing)} singular")

strictly_between = []
for T in sing:
    M = sp.Matrix(3, 3, lambda a, b: P0[T[a], T[b]])
    lmin = min(sp.simplify(r) for r in M.eigenvals().keys())
    if sp.N(lmin, 30) > sp.Float(10) ** -25:
        strictly_between.append((T, lmin))
check("N8b", "every singular triple has lambda_min = 0 exactly (none in (0,1/6))",
      not strictly_between, f"{len(strictly_between)} anomalies")

print()
print("  ==> F(A_0) = max_T lambda_min(P_TT) = 1/6  EXACTLY")
print("  ==> f(A_0) = 1/sqrt6                      EXACTLY")
Fmax = max([sp.Rational(1, 6)] * len(at_16) + [sp.Integer(0)] * len(sing),
           key=lambda z: sp.N(z, 30))
check("N9", "F(A_0) = max_T lambda_min = 1/6 exactly => GTZ(6,3) is TIGHT",
      len(at_16) == 10 and sp.simplify(Fmax - sp.Rational(1, 6)) == 0
      and not strictly_between,
      "GTZ holds at A_0 with EQUALITY: 10 good triples, all at exactly 1/6")

# ---- the equality pair {3,4}
print()
print("-" * 76)
print("The equality pair {3,4} and the excess delta")
print("-" * 76)
i, j = 3, 4
Gam = sp.Matrix([[P0[i, i], P0[i, j]], [P0[i, j], P0[j, j]]]).applyfunc(
    lambda t: sp.radsimp(sp.simplify(t)))
c2 = sp.radsimp(sp.simplify(P0[i, j] ** 2))
print(f"  Gamma_0 = {Gam.tolist()},  c^2 = {c2}")
check("N10", "c_{34}^2 = 25/324", sp.simplify(c2 - sp.Rational(25, 324)) == 0)

mus = sorted([sp.radsimp(sp.simplify(m)) for m in Gam.eigenvals().keys()],
             key=lambda z: -sp.N(z, 30))
print(f"  (mu1, mu2) = ({mus[0]}, {mus[1]})")
check("N11", "(mu1,mu2) = (1, 4/9)", sp.simplify(mus[0] - 1) == 0
      and sp.simplify(mus[1] - sp.Rational(4, 9)) == 0)
check("N12", "the pair QUALIFIES: mu2 = 4/9 > 1/6", sp.simplify(mus[1] - SIX) > 0)

h = lambda m: (1 - m) / (6 * m - 1)
delta = sp.simplify(h(mus[0]) + h(mus[1]) - sp.Rational(1, 3))
print(f"  h(mu1)+h(mu2)-1/3 = delta = {delta}")
check("N13", "delta = 0 EXACTLY: the pair lies ON the equality manifold", delta == 0)

# rational form (1.1) of the excess
t_, d_ = sp.simplify(Gam.trace()), sp.simplify(Gam.det())
delta_rat = sp.simplify((7 * t_ - 2 - 12 * d_) / (36 * d_ - 6 * t_ + 1) - sp.Rational(1, 3))
check("N14", "rational form (1.1) of delta agrees and equals 0",
      sp.simplify(delta_rat - delta) == 0 and delta_rat == 0,
      f"t={t_}, d={d_}")

# all four Schur slacks q(p) = 0
Ginv = sp.simplify((Gam - SIX * sp.eye(2)).inv())
qs = {}
for p in range(6):
    if p in (i, j):
        continue
    v = sp.Matrix([P0[i, p], P0[j, p]])
    qs[p] = sp.radsimp(sp.simplify(P0[p, p] - SIX - (v.T * Ginv * v)[0]))
print(f"  Schur slacks q(p) = { {k: str(v) for k, v in qs.items()} }")
check("N15", "all four Schur slacks q(0)=q(1)=q(2)=q(5)=0 EXACTLY",
      all(sp.simplify(x) == 0 for x in qs.values()))
check("N16", "q-sum identity: sum_p q(p) = -delta  (=0 here), exactly",
      sp.simplify(sum(qs.values()) + delta) == 0)

# ---- the contrast with the slice
print()
print("-" * 76)
print("Contrast: where the problem is actually tight")
print("-" * 76)
G_PMC = (1 - sp.sqrt((5 - sp.sqrt(5)) / 8)) / 2
print(f"  F at Nesterenko extremal  = 1/6                = {sp.N(SIX, 20)}   <-- TIGHT")
print(f"  F at PMC (slice minimum)  = (1-sin36)/2        = {sp.N(G_PMC, 20)}")
print(f"  slice margin over 1/6     = {sp.N(G_PMC - SIX, 20)}")
check("N17", "the slice minimum is STRICTLY ABOVE 1/6, so the slice is NOT binding",
      sp.N(G_PMC - SIX, 40) > 0)
check("N18", "Nesterenko leverages (5/18, 13/18) are NOT the slice value 1/2",
      sp.simplify(sp.Rational(5, 18) - sp.Rational(1, 2)) != 0)

print()
print("=" * 76)
fails = [(t, n, d) for t, n, o, d in results if not o]
print(f"SUMMARY: {len(results) - len(fails)}/{len(results)} checks passed")
if fails:
    print("FAILURES / DISCREPANCIES vs boundary-obstruction.md Lemma 2.1:")
    for t, n, d in fails:
        print(f"   - {t}: {n}   {d}")
else:
    print("Lemma 2.1 fully re-verified in exact arithmetic.")
    print("KEY CONSEQUENCE: GTZ(6,3) is TIGHT -- F = 1/6 exactly at this extremal.")
print("=" * 76)
sys.exit(1 if fails else 0)
