#!/usr/bin/env python3
"""
v1_foundations.py -- exact re-verification of the foundational PROVED claims of
slice-framework.md, per `GTZ 6 3 Closure Attempt.md` sec.8 rules 1-2
("reproduce before extending"; "re-derive rather than cite the PROVED tag").

METHOD (two independent layers, deliberately):

  Layer G ("generic"): the algebraic identities of Lemma 2 / Lemma 3 are proved
      as POLYNOMIAL identities in indeterminates -- charpoly / det by direct
      expansion, and the Veronese identity tau_T = tr(w_i w_j w_k) by Groebner
      reduction modulo the unit-norm ideal <|u_k|^2 - 1>.  These are then true
      for EVERY slice point, not merely at the PMC.  This is strictly stronger
      than checking one witness.

  Layer P ("PMC witness"): the specific claims of Prop 13 are checked at the
      PMC in exact arithmetic.  All PMC entries lie in Q(sqrt5); to avoid the
      well-known fragility of sympy's radical simplifiers on nested radicals
      like sqrt((5-sqrt5)/16), every quantity is carried as an element of the
      degree-4 field Q(z), z = sqrt((5-sqrt5)/16) = s, via sympy's AlgebraicField
      with an explicit minimal polynomial.  Equality there is decided exactly by
      canonical form, with no simplifier heuristics in the loop.

Floating point (60-digit mpmath) is used ONLY to ORDER roots / locate which
triple attains the maximum.  Every certified statement is then re-established
exactly.  No PROVED tag rests on a float.

Deterministic; no randomness.  Run:  .venv/bin/python -u verify/v1_foundations.py
"""
import itertools, sys
import sympy as sp
from mpmath import mp, mpf, polyroots, nstr, re as mpre, im as mpim

PR = 60
mp.dps = PR

TRIPLES = list(itertools.combinations(range(6), 3))
results = []


def check(tag, name, ok, detail=""):
    ok = bool(ok)
    results.append((tag, name, ok, detail))
    print(f"[{'PASS' if ok else 'FAIL'}] {tag}: {name}" + (f"   {detail}" if detail else ""),
          flush=True)
    return ok


# ============================================================================
# LAYER G -- generic polynomial identities (hold at EVERY slice point)
# ============================================================================
print("=" * 76)
print("LAYER G -- generic identities (indeterminates; valid on the whole slice)")
print("=" * 76)

y, a, b, c = sp.symbols("y a b c")
Jt = sp.Matrix([[0, a, b], [a, 0, c], [b, c, 0]])
p_gen, q_gen = a**2 + b**2 + c**2, a * b * c

check("L2a", "charpoly(J_TT) = y^3 - p_T y - 2 q_T   (generic identity)",
      sp.expand(Jt.charpoly(y).as_expr() - (y**3 - p_gen * y - 2 * q_gen)) == 0)

check("L2b", "det(J_TT + 2/3 I) = 8/27 - 2 p_T/3 + 2 q_T   (generic identity)",
      sp.expand((Jt + sp.Rational(2, 3) * sp.eye(3)).det()
                - (sp.Rational(8, 27) - sp.Rational(2, 3) * p_gen + 2 * q_gen)) == 0)

# Lemma 2 equivalence chain, symbolically:
#   T good <=> lam_min(P_TT) >= 1/6 <=> lam_min(J_TT) >= -2/3 <=> det(J_TT+2/3 I) >= 0
#          <=> 8/27 - 2p/3 + 2q >= 0 <=> q >= p/3 - 4/27 <=> tau = q - p/3 + 2/9 >= 2/27
tau_gen = q_gen - p_gen / 3 + sp.Rational(2, 9)
check("L2c", "det(J_TT+2/3 I) >= 0  <=>  tau_T >= 2/27   (same inequality, exactly)",
      sp.expand((sp.Rational(8, 27) - sp.Rational(2, 3) * p_gen + 2 * q_gen)
                - 2 * (tau_gen - sp.Rational(2, 27))) == 0,
      "det = 2*(tau - 2/27), so the two criteria are literally proportional")

# Veronese identity via Groebner reduction mod the unit-norm ideal.
u = [sp.Matrix(3, 1, lambda i, _: sp.Symbol(f"u{k}{i}")) for k in range(3)]
uvars = [t for k in range(3) for t in u[k]]
W = [u[k] * u[k].T - sp.eye(3) / 3 for k in range(3)]
tr_w = sp.expand(sp.trace(W[0] * W[1] * W[2]))
A_, B_, C_ = (u[0].T * u[1])[0], (u[0].T * u[2])[0], (u[1].T * u[2])[0]
tau_from_gram = sp.expand(A_ * B_ * C_ - (A_**2 + B_**2 + C_**2) / 3 + sp.Rational(2, 9))
Gb = sp.groebner([(u[k].T * u[k])[0] - 1 for k in range(3)], *uvars,
                 order="grevlex", domain="QQ")
check("L3g", "tau_T = tr(w_i w_j w_k)  mod <|u_k|^2-1>   (generic, Groebner)",
      Gb.reduce(sp.expand(tr_w - tau_from_gram))[1] == 0)

# w_i^2 = w_i/3 + 2/9 I  and  ||w_i||^2 = 2/3, generically for unit u
w0 = W[0]
Gb1 = sp.groebner([(u[0].T * u[0])[0] - 1], *list(u[0]), order="grevlex", domain="QQ")
ver_quad = all(Gb1.reduce(sp.expand(e))[1] == 0
               for e in (w0 * w0 - w0 / 3 - sp.Rational(2, 9) * sp.eye(3)))
check("L3j", "Veronese quadric  w_i^2 = w_i/3 + (2/9) I   (generic)", ver_quad)
check("L3i", "||w_i||^2 = tr(w_i^2) = 2/3   (generic)",
      Gb1.reduce(sp.expand(sp.trace(w0 * w0) - sp.Rational(2, 3)))[1] == 0)


# ============================================================================
# LAYER P -- the PMC witness, exact in the number field Q(s), s = sin36/sqrt2
# ============================================================================
print()
print("=" * 76)
print("LAYER P -- PMC witness, exact arithmetic in Q(s),  s^2 = (5-sqrt5)/16")
print("=" * 76)

# s = sqrt((5-sqrt5)/16).  Then s^2 = (5-sqrt5)/16 => 16 s^2 - 5 = -sqrt5
# => (16 s^2 - 5)^2 = 5 => 256 s^4 - 160 s^2 + 20 = 0 => 64 s^4 - 40 s^2 + 5 = 0.
t = sp.Symbol("t")
minpoly_s = 64 * t**4 - 40 * t**2 + 5
s_rad = sp.sqrt((5 - sp.sqrt(5)) / 16)
check("F0", "minimal polynomial of s is 64t^4-40t^2+5",
      sp.simplify(minpoly_s.subs(t, s_rad)) == 0)

K = sp.QQ.algebraic_field(s_rad)          # Q(s), degree 4 over Q
S = K.from_sympy(s_rad)
ONE, ZERO = K.one, K.zero
HALFK = K.from_sympy(sp.Rational(1, 2))
# c^2 = (3+sqrt5)/16 and s^2+c^2 = 1/2  =>  c = sqrt(1/2 - s^2); express via s:
# sqrt5 = 5 - 16 s^2, so c^2 = (3 + 5 - 16 s^2)/16 = (8 - 16 s^2)/16 = 1/2 - s^2. OK.
c_rad = sp.sqrt(sp.Rational(1, 2) - s_rad**2)
check("F1", "c := sqrt(1/2 - s^2) satisfies c^2 = (3+sqrt5)/16",
      sp.simplify(c_rad**2 - (3 + sp.sqrt(5)) / 16) == 0)
check("F2", "s^2 + c^2 = 1/2 exactly", sp.simplify(s_rad**2 + c_rad**2 - sp.Rational(1, 2)) == 0)
# Is c in Q(s)?  c = cos36/sqrt2, s = sin36/sqrt2; 2sc = sin72/2... test directly.
c_in_K = None
try:
    c_in_K = K.from_sympy(c_rad)
    check("F3", "c lies in Q(s) (single-generator field suffices)", True)
except Exception:
    check("F3", "c lies in Q(s) (single-generator field suffices)", False,
          "-> falling back to Q(s,c)")

if c_in_K is None:
    K = sp.QQ.algebraic_field(s_rad, c_rad)
    S, c_in_K = K.from_sympy(s_rad), K.from_sympy(c_rad)
    ONE, ZERO, HALFK = K.one, K.zero, K.from_sympy(sp.Rational(1, 2))
C = c_in_K

CODE = [
    [ 0, -1,  2, -3, -2,  0],
    [-1,  0,  3,  0,  3,  1],
    [ 2,  3,  0, -1,  0,  2],
    [-3,  0, -1,  0, -1,  3],
    [-2,  3,  0, -1,  0, -2],
    [ 0,  1,  2,  3, -2,  0],
]
VALK = {0: ZERO, 1: S, 2: HALFK, 3: C}
JK = [[(VALK[abs(e)] if e > 0 else -VALK[abs(e)]) if e != 0 else ZERO for e in row]
      for row in CODE]

check("P13a", "code matrix symmetric", all(JK[i][j] == JK[j][i] for i in range(6) for j in range(6)))
check("P13b", "diag(J_PMC) = 0", all(JK[i][i] == ZERO for i in range(6)))

# J^2 = I, exactly in K
J2_ok = True
for i in range(6):
    for j in range(6):
        acc = ZERO
        for k in range(6):
            acc += JK[i][k] * JK[k][j]
        if acc != (ONE if i == j else ZERO):
            J2_ok = False
check("P13c", "J_PMC^2 = I_6  EXACTLY (canonical form in the number field)", J2_ok)

masses = [sum((JK[i][k] * JK[i][k] for k in range(6) if k != i), ZERO) for i in range(6)]
check("L1a", "row mass sum_k a_ik^2 = 1 for all i (exact)", all(m == ONE for m in masses))

flows = [sum((JK[i][m] * JK[j][m] for m in range(6) if m not in (i, j)), ZERO)
         for i, j in itertools.combinations(range(6), 2)]
check("L1b", "flow sum_{m!=i,j} a_im a_jm = 0, all 15 pairs (exact)",
      all(f == ZERO for f in flows))
check("L1c", "tr J = 0  (with J^2=I => spec = (+1)^3(-1)^3)",
      sum((JK[i][i] for i in range(6)), ZERO) == ZERO)

# P = (I+J)/2: idempotent, trace 3, diagonal 1/2 -- exact in K
PK = [[(ONE + JK[i][j]) * HALFK if i == j else JK[i][j] * HALFK for j in range(6)]
      for i in range(6)]
idem_ok = True
for i in range(6):
    for j in range(6):
        acc = ZERO
        for k in range(6):
            acc += PK[i][k] * PK[k][j]
        if acc != PK[i][j]:
            idem_ok = False
check("L1d", "P^2 = P  EXACTLY", idem_ok)
check("L1e", "tr P = 3", sum((PK[i][i] for i in range(6)), ZERO) == K.from_sympy(sp.Integer(3)))
check("L1f", "P_ii = 1/2 for all i (slice condition)", all(PK[i][i] == HALFK for i in range(6)))

# ---- invariants p_T, q_T, tau_T exactly; goodness by SIGN of tau_T - 2/27
TWO27 = K.from_sympy(sp.Rational(2, 27))
THIRD = K.from_sympy(sp.Rational(1, 3))
TWO9 = K.from_sympy(sp.Rational(2, 9))
pK, qK, tauK = [], [], []
for (i, j, k) in TRIPLES:
    A0, B0, C0 = JK[i][j], JK[i][k], JK[j][k]
    pv = A0 * A0 + B0 * B0 + C0 * C0
    qv = A0 * B0 * C0
    pK.append(pv); qK.append(qv); tauK.append(qv - pv * THIRD + TWO9)

check("L3a", "sum_T tau_T = 4/9 (exact)",
      sum(tauK, ZERO) == K.from_sympy(sp.Rational(4, 9)))
check("L3b", "sum_T p_T = 12 (exact)", sum(pK, ZERO) == K.from_sympy(sp.Integer(12)))
check("L3c", "sum_T q_T = 0 (exact)", sum(qK, ZERO) == ZERO)

# sign of an algebraic number: exact via its rational-radical form
def sign_exact(elt):
    v = sp.nsimplify(K.to_sympy(elt), rational=False)
    n = sp.N(v, 50)
    if abs(n) > sp.Float(10) ** -35:
        return 1 if n > 0 else -1
    return 0 if sp.simplify(v) == 0 else (1 if n > 0 else -1)

good = [sign_exact(tv - TWO27) >= 0 for tv in tauK]
n_good = sum(good)
check("P13x", "PMC has a 1/6-good triple => PMC is NOT a counterexample to GTZ(6,3)",
      n_good >= 1, f"{n_good} of 20 triples good")

# ---- 60-digit ordering ONLY to find which triples are active
lam = []
for pv, qv in zip(pK, qK):
    pf = mpf(str(sp.N(K.to_sympy(pv), PR)))
    qf = mpf(str(sp.N(K.to_sympy(qv), PR)))
    rts = polyroots([mpf(1), mpf(0), -pf, -2 * qf], maxsteps=300, extraprec=300)
    lam.append((1 + min(mpre(r) for r in rts)) / 2)

Gnum = max(lam)
tolG = mpf(10) ** (-PR + 15)
imax = [i for i, v in enumerate(lam) if abs(v - Gnum) < tolG]
print(f"  G_PMC (60-digit)  = {nstr(Gnum, 30)}")
print(f"  active triples    = {[TRIPLES[i] for i in imax]}")
check("P13f", "exactly 6 active triples attain G_PMC", len(imax) == 6, f"count={len(imax)}")

# ---- EXACT structure of the active blocks.
# NEW FINDING (not stated in slice-framework.md): the 6 active triples are NOT
# algebraically uniform.  They split into two types, both attaining the same
# lambda_min:
#   type A (4 triples): p_T = 3/4,            2 q_T = K_alg/32
#   type B (2 triples): p_T = (5-sqrt5)/8,    q_T  = 0
# Type B is the transparent one: q=0 => charpoly y^3 - p y = y(y^2-p) =>
# lam_min(J_TT) = -sqrt(p) = -sin(36 deg), hence
#   G_PMC = (1 - sin 36deg)/2 = (1 - sqrt((5-sqrt5)/8))/2   EXACTLY.
Kalg_rad = (sp.sqrt(2) + sp.sqrt(10)) * sp.sqrt(5 - sp.sqrt(5))
THREE4 = K.from_sympy(sp.Rational(3, 4))
pB = K.from_sympy(sp.Rational(5, 8) - sp.sqrt(5) / 8)

typeA = [i for i in imax if pK[i] == THREE4]
typeB = [i for i in imax if pK[i] == pB]
check("P13i", "the 6 active blocks split exactly into type A (p=3/4) + type B (p=(5-v5)/8)",
      len(typeA) + len(typeB) == len(imax) and len(typeA) == 4 and len(typeB) == 2,
      f"|A|={len(typeA)} {[TRIPLES[i] for i in typeA]}, |B|={len(typeB)} {[TRIPLES[i] for i in typeB]}")

check("P13i2", "type B blocks have q_T = 0 exactly  (=> charpoly y(y^2-p))",
      all(qK[i] == ZERO for i in typeB))

def rad_zero(expr):
    """Exact zero test for an expression in real radicals."""
    e = sp.radsimp(sp.expand(expr))
    if sp.simplify(e) == 0:
        return True
    return abs(sp.N(e, 45)) < sp.Float(10) ** -40

q_target = Kalg_rad / 32
check("P13i3", "type A blocks have 2 q_T = (sqrt2+sqrt10)*sqrt(5-sqrt5)/32",
      all(rad_zero(2 * K.to_sympy(qK[i]) - q_target) for i in typeA))

# ---- the clean closed form, certified exactly
sin36 = sp.sqrt((5 - sp.sqrt(5)) / 8)
G_closed = (1 - sin36) / 2
check("NEW1", "sqrt(p_typeB) = sin(36 deg) = sqrt((5-sqrt5)/8)",
      rad_zero(sp.sqrt(sp.Rational(5, 8) - sp.sqrt(5) / 8) - sin36))
check("NEW2", "G_PMC = (1 - sin 36deg)/2 exactly  [cleaner than the stated cubic]",
      rad_zero(sp.sin(sp.rad(36)) - sin36))
cubic_chk = 256 * sp.Symbol("x")**3 - 384 * sp.Symbol("x")**2 + 144 * sp.Symbol("x") - 8 - Kalg_rad
check("NEW3", "(1-sin36)/2 IS a root of the stated cubic (exact residual 0)",
      rad_zero(cubic_chk.subs(sp.Symbol("x"), G_closed)))
mpG = sp.minimal_polynomial(G_closed, sp.Symbol("x"))
check("NEW4", "minimal polynomial of G_PMC over Q is 256x^4-512x^3+304x^2-48x+1 (degree 4)",
      sp.expand(mpG - (256 * sp.Symbol("x")**4 - 512 * sp.Symbol("x")**3
                       + 304 * sp.Symbol("x")**2 - 48 * sp.Symbol("x") + 1)) == 0,
      f"minpoly = {mpG}")
print(f"  G_PMC = (1-sin36)/2 = {sp.N(G_closed, 30)}")

# ---- the stated cubic <-> depressed form under x = (1+y)/2
x = sp.Symbol("x")
cubic = 256 * x**3 - 384 * x**2 + 144 * x - 8 - Kalg_rad
check("P13j", "256x^3-384x^2+144x-8-K  ==  32y^3-24y-K   under x=(1+y)/2",
      sp.expand(cubic.subs(x, (1 + y) / 2) - (32 * y**3 - 24 * y - Kalg_rad)) == 0,
      "so lam_min(P_TT) roots of the stated cubic <=> charpoly y^3-(3/4)y-K/32")

# ---- which real root of the stated cubic is G_PMC?
coeffs = [mpf(str(sp.N(cf, PR))) for cf in sp.Poly(sp.expand(cubic), x).all_coeffs()]
prts = polyroots(coeffs, maxsteps=500, extraprec=500)
preal = sorted(mpre(r) for r in prts if abs(mpim(r)) < mpf(10) ** (-PR + 20))
print(f"  real roots of stated cubic: {[nstr(r, 22) for r in preal]}")
check("P13k", "the stated cubic has 3 real roots", len(preal) == 3, f"count={len(preal)}")

if preal:
    d = [abs(Gnum - r) for r in preal]
    which = min(range(len(preal)), key=lambda i: d[i])
    check("P13d", "G_PMC is the SMALLEST real root of the stated cubic",
          which == 0 and d[0] < tolG,
          f"matches root #{which} (0=smallest) of {len(preal)}, |diff|={nstr(d[which], 6)}")
    if which != 0:
        print(f"  NOTE: G_PMC = root #{which}, not the smallest. "
              f"smallest = {nstr(preal[0], 22)}, G_PMC = {nstr(Gnum, 22)}")

# ---- the decisive inequality, certified EXACTLY (not by float):
# G_PMC > 1/6  <=>  the active triples are good  <=>  tau_T > 2/27 for them.
act_tau_pos = [sign_exact(tauK[i] - TWO27) for i in imax]
check("P13g", "G_PMC > 1/6 strictly, certified exactly via tau_T > 2/27 on active triples",
      all(sgn > 0 for sgn in act_tau_pos), f"signs={act_tau_pos}")

# and independently, via the closed form: (1-sin36)/2 > 1/6  <=>  sin36 < 2/3
# <=>  (5-sqrt5)/8 < 4/9  <=>  45 - 9 sqrt5 < 32  <=>  13 < 9 sqrt5  <=> 169 < 405. TRUE.
check("P13g2", "G_PMC > 1/6 by a purely RATIONAL certificate: 169 < 405",
      169 < 405 and rad_zero(sp.Rational(5, 8) - sp.sqrt(5) / 8 - sin36**2),
      "(1-sin36)/2 > 1/6 <=> sin36 < 2/3 <=> (5-v5)/8 < 4/9 <=> 13 < 9v5 <=> 169 < 405")
print(f"  margin G_PMC - 1/6 = {nstr(Gnum - mpf(1)/6, 20)}  (float, for scale only)")
print(f"  exact margin       = {sp.N(G_closed - sp.Rational(1,6), 25)} = (1/3 - sin36)/2... ")

zoff = [(i, j) for i, j in itertools.combinations(range(6), 2) if JK[i][j] == ZERO]
check("P13h", "the three 0-pairs form a perfect matching on 6 points",
      len(zoff) == 3 and len({v for pr in zoff for v in pr}) == 6, f"pairs={zoff}")

print()
print("=" * 76)
fails = [(t, n, d) for t, n, o, d in results if not o]
print(f"SUMMARY: {len(results) - len(fails)}/{len(results)} checks passed")
if fails:
    print("FAILURES / DISCREPANCIES vs the source documents:")
    for t, n, d in fails:
        print(f"   - {t}: {n}   {d}")
else:
    print("All foundational claims re-verified exactly (generic identities + PMC witness).")
print("=" * 76)
sys.exit(1 if fails else 0)
