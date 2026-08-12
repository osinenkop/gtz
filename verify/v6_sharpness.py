#!/usr/bin/env python3
"""
v6_sharpness.py -- EXACT certificate for the sharp-minimum claim at the
Nesterenko extremal A_0, per sec.6 of proofs/sharp-cone-at-extremal.md.

SETUP.  Work on the Grassmannian Gr(3,6) = {P sym : P^2=P, tr P=3} (9-dimensional,
no O(3) gauge redundancy -- this matters, see the WARNING below).  At P_0 the
active set is the 10 triples with lambda_min(P_0,TT) = 1/6 exactly (v4).  The
other 10 triples are singular with lambda_min = 0, far below 1/6, so they do not
participate in the max near P_0.

For a tangent direction Pdot, Rellich/Danskin gives the one-sided derivative
    F'(P_0; Pdot) = max_{T active}  min_{v in E_T, |v|=1}  v^T Pdot_TT v,
where E_T is the lambda_min-eigenspace of P_0,TT.  When each E_T is 1-dimensional
(checked exactly below) this is
    F'(P_0; Pdot) = max_{T active} <G_T, Pdot>_F
for explicit symmetric matrices G_T -- a max of 10 LINEAR functionals on a
9-dimensional space.

THE CERTIFICATE.  Define the critical cone N = {Pdot : <G_T,Pdot> <= 0 for all
active T}.  Then:
    * N = {0}   <=>  max_T <G_T,Pdot> > 0 for every Pdot != 0
                <=>  kappa := min_{|Pdot|=1} max_T <G_T,Pdot> > 0    (compactness)
                <=>  A_0 is a SHARP (first-order) local minimum.
    * A finite set {G_T} positively spans a 9-dim space  <=>  it spans linearly
      (rank 9) AND 0 is a STRICTLY positive combination sum lambda_T G_T = 0,
      lambda_T > 0.  With 10 vectors in 9 dimensions the null space is generically
      1-dimensional, so this reduces to: rank = 9, nullity = 1, and the null
      vector has all 10 entries of the same sign.
This is an exact rational/algebraic linear-algebra check -- NOT a
Positivstellensatz problem.  Note sum lambda_T G_T = 0 is precisely the KKT
condition (5.1) of boundary-obstruction.md; STRICT positivity of every lambda_T
is the extra ingredient that upgrades "KKT point" to "sharp minimum".

WARNING THIS SCRIPT EXISTS TO SETTLE.  The kappa ~ 6.59e-3 reported in
sharp-cone-at-extremal.md sec.4 was minimized over ambient directions in the
18-dimensional A-space.  That space contains the 3-dimensional O(3) GAUGE
directions A -> A*Omega, along which P (hence F) is exactly CONSTANT, so the true
ambient minimum is 0, attained on a 3-dim subspace, and any ambient slope
minimization should drive to 0, not to a positive floor.  The honest question is
whether kappa > 0 holds TRANSVERSALLY, i.e. on the 9-dim Grassmannian tangent
space with the gauge quotiented out.  That is what is computed here, exactly.
Section 4's numeric kappa must be treated as unreliable until this settles it.

All arithmetic exact (sympy, Q(sqrt5)).  Deterministic; no randomness.
"""
import itertools, sys
import sympy as sp

TRIPLES = list(itertools.combinations(range(6), 3))
SIX = sp.Rational(1, 6)
results = []


def check(tag, name, ok, detail=""):
    ok = bool(ok)
    results.append((tag, name, ok, detail))
    print(f"[{'PASS' if ok else 'FAIL'}] {tag}: {name}" + (f"   {detail}" if detail else ""),
          flush=True)
    return ok


def R(e):
    return sp.radsimp(sp.simplify(sp.expand(e)))


print("=" * 78)
print("EXACT SHARPNESS CERTIFICATE AT THE NESTERENKO EXTREMAL")
print("=" * 78)

# ---------------------------------------------------------------- exact P_0
Bred = sp.Matrix([[-1, -1, -1, 0, 0, -1],
                  [0, 0, 0, -1, 1, 0],
                  [0, 0, 0, 0, -1, 1]])
w = [sp.Integer(1)] * 3 + [sp.Rational(9, 5)] * 3
Mt = sp.simplify(sp.diag(*[sp.sqrt(x) for x in w]) * Bred.T)
P0 = sp.simplify(Mt * (Mt.T * Mt).inv() * Mt.T).applyfunc(R)

check("S0", "P0^2 = P0, tr P0 = 3 (exact)",
      sp.simplify(sp.expand(P0 * P0 - P0)) == sp.zeros(6, 6)
      and sp.simplify(sp.trace(P0)) == 3)

# ------------------------------------------------- active set and eigenspaces
active, singular = [], []
for T in TRIPLES:
    M = sp.Matrix(3, 3, lambda a, b: P0[T[a], T[b]]).applyfunc(R)
    if R(M.det()) == 0:
        singular.append(T)
    else:
        active.append(T)

print(f"\n  active (non-singular, lambda_min = 1/6) triples: {len(active)}")
print(f"  singular (lambda_min = 0) triples:               {len(singular)}")
check("S1", "exactly 10 active triples", len(active) == 10)

# E_T = kernel of (P_0,TT - I/6);  verify lambda_min = 1/6 and its multiplicity
Evecs, mults = {}, {}
for T in active:
    M = sp.Matrix(3, 3, lambda a, b: P0[T[a], T[b]]).applyfunc(R)
    D = (M - SIX * sp.eye(3)).applyfunc(R)
    ker = D.nullspace()
    mults[T] = len(ker)
    if len(ker) == 1:
        v = ker[0].applyfunc(R)
        nrm2 = R((v.T * v)[0])
        Evecs[T] = (v / sp.sqrt(nrm2)).applyfunc(R)      # unit eigenvector
    else:
        Evecs[T] = [k.applyfunc(R) for k in ker]

print(f"  multiplicities of the 1/6 eigenvalue: {sorted(set(mults.values()))}")
check("S2", "1/6 is a SIMPLE eigenvalue of every active block (dim E_T = 1)",
      all(m == 1 for m in mults.values()),
      f"multiplicities = { {str(k): v for k, v in mults.items()} }"
      if not all(m == 1 for m in mults.values()) else "")

if not all(m == 1 for m in mults.values()):
    print("\n  dim E_T > 1 for some T: the derivative is a min-eigenvalue of a")
    print("  compression, not a single linear functional.  Handled below only if")
    print("  all dims are 1; otherwise this script stops (needs the SDP form).")
    sys.exit(1)

# also confirm 1/6 is the SMALLEST eigenvalue (not just an eigenvalue)
lmin_ok = True
for T in active:
    M = sp.Matrix(3, 3, lambda a, b: P0[T[a], T[b]]).applyfunc(R)
    D = (M - SIX * sp.eye(3)).applyfunc(R)
    # D PSD via ALL principal minors (the brief's sec.4 pitfall)
    for k in (1, 2, 3):
        for idx in itertools.combinations(range(3), k):
            if not (sp.simplify(D[list(idx), list(idx)].det()) >= 0):
                lmin_ok = False
check("S3", "P_0,TT - I/6 is PSD for every active T (all principal minors >= 0)",
      lmin_ok, "so 1/6 really is lambda_min, not an interior eigenvalue")

# --------------------------------------- tangent space to Gr(3,6) at P_0 (9-dim)
# Pdot = sum_{a,b} z_ab (u_a n_b^T + n_b u_a^T),  u_a in range(P_0), n_b in ker(P_0).
rng_basis = [v.applyfunc(R) for v in P0.columnspace()][:3]
ker_basis = [v.applyfunc(R) for v in P0.nullspace()][:3]
check("S4", "range(P0) is 3-dim and ker(P0) is 3-dim",
      len(rng_basis) == 3 and len(ker_basis) == 3)


def gram_schmidt(vs):
    out = []
    for v in vs:
        u = v
        for o in out:
            u = (u - (R((o.T * v)[0])) * o).applyfunc(R)
        n2 = R((u.T * u)[0])
        out.append((u / sp.sqrt(n2)).applyfunc(R))
    return out


U = gram_schmidt(rng_basis)
Nn = gram_schmidt(ker_basis)
orth_ok = all(R((U[i].T * U[j])[0]) == (1 if i == j else 0) for i in range(3) for j in range(3)) \
    and all(R((Nn[i].T * Nn[j])[0]) == (1 if i == j else 0) for i in range(3) for j in range(3)) \
    and all(R((U[i].T * Nn[j])[0]) == 0 for i in range(3) for j in range(3))
check("S5", "orthonormal bases of range(P0) and ker(P0), mutually orthogonal", orth_ok)

TAN = []
for a in range(3):
    for b in range(3):
        B = (U[a] * Nn[b].T + Nn[b] * U[a].T).applyfunc(R)
        TAN.append(B)
check("S6", "tangent space has dimension 9 = 3*(6-3)", len(TAN) == 9)

# tangent vectors are symmetric, trace-free, and satisfy the Grassmann condition
# Pdot = P Pdot (I-P) + (I-P) Pdot P  (equivalently P Pdot P = 0 and (I-P)Pdot(I-P)=0)
I6 = sp.eye(6)
tan_ok = True
for B in TAN:
    if sp.simplify(B - B.T) != sp.zeros(6, 6):
        tan_ok = False
    if R(sp.trace(B)) != 0:
        tan_ok = False
    if sp.simplify(sp.expand(P0 * B * P0)) != sp.zeros(6, 6):
        tan_ok = False
    if sp.simplify(sp.expand((I6 - P0) * B * (I6 - P0))) != sp.zeros(6, 6):
        tan_ok = False
check("S7", "every basis tangent vector: symmetric, traceless, P B P = 0, (I-P)B(I-P) = 0",
      tan_ok)

# Frobenius-orthonormality of the tangent basis (so |Pdot| is the coefficient norm)
fro = sp.Matrix(9, 9, lambda i, j: R(sp.trace(TAN[i] * TAN[j])))
check("S8", "tangent basis is Frobenius-orthogonal with norm^2 = 2",
      fro == 2 * sp.eye(9), f"Gram diag = {[fro[i, i] for i in range(9)]}")

# ------------------------------------- the 10 linear functionals L_T(Pdot)
# L_T(Pdot) = v_T^T Pdot_TT v_T ; express in tangent coordinates z (9 vector).
Lrows = []
for T in active:
    v = Evecs[T]
    row = []
    for B in TAN:
        Bt = sp.Matrix(3, 3, lambda a, b: B[T[a], T[b]])
        row.append(R((v.T * Bt * v)[0]))
    Lrows.append(row)

Lmat = sp.Matrix(Lrows)                      # 10 x 9 ; L(z) = Lmat * z
print(f"\n  functional matrix L is {Lmat.shape[0]} x {Lmat.shape[1]}")

rk = Lmat.rank()
check("S9", "rank(L) = 9: the active gradients span the tangent space", rk == 9,
      f"rank = {rk}")

# ---------------- KKT: is there lambda >= 0, lambda != 0, with L^T lambda = 0 ?
ns = Lmat.T.nullspace()                       # lambda with sum_T lambda_T G_T = 0
print(f"  dim of left-null space (KKT multipliers): {len(ns)}")
check("S10", "the KKT multiplier space is exactly 1-dimensional", len(ns) == 1,
      f"dim = {len(ns)}")

if len(ns) == 1:
    lam = ns[0].applyfunc(R)
    ent = [R(lam[i]) for i in range(lam.rows)]
    signs = [sp.sign(sp.nsimplify(e)) for e in ent]
    print(f"  multiplier signs: {signs}")
    allpos = all(s == 1 for s in signs)
    allneg = all(s == -1 for s in signs)
    if allneg:
        lam = (-lam).applyfunc(R)
        ent = [R(lam[i]) for i in range(lam.rows)]
    strict = allpos or allneg
    check("S11", "KKT: sum_T lambda_T G_T = 0 with ALL 10 multipliers STRICTLY nonzero "
                 "and of one sign", strict,
          "=> the active gradients POSITIVELY span => critical cone N = {0}")
    # normalize to a probability vector and show it
    tot = R(sum(ent))
    if strict and tot != 0:
        prob = [R(e / tot) for e in ent]
        print(f"  normalized multipliers (sum = 1):")
        for T, p in zip(active, prob):
            print(f"     {T}: {p}")
        # verify the KKT identity exactly
        acc = sp.zeros(6, 6)
        for T, p in zip(active, prob):
            v = Evecs[T]
            Gfull = sp.zeros(6, 6)
            for a in range(3):
                for b in range(3):
                    Gfull[T[a], T[b]] = R(v[a] * v[b])
            acc = (acc + p * Gfull).applyfunc(R)
        # projection of acc onto the tangent space must vanish
        proj = [R(sp.trace(acc * B)) for B in TAN]
        check("S12", "exact verification: sum_T lambda_T G_T projects to 0 on the "
                     "tangent space", all(p == 0 for p in proj),
              f"residuals = {proj}")

    # ------------------------------------------------ the sharpness constant
    if strict:
        print("\n" + "-" * 78)
        print("SHARPNESS: kappa = min_{|z|=1} max_T (L z)_T  > 0")
        print("-" * 78)
        print("  Positive spanning (S9+S11) already PROVES kappa > 0:")
        print("  if max_T (Lz)_T <= 0 for some z, then 0 = <sum lam_T G_T, z>")
        print("  = sum lam_T (Lz)_T < 0 unless every (Lz)_T = 0, and rank L = 9")
        print("  forces z = 0.  Hence max_T (Lz)_T > 0 for all z != 0, and the min")
        print("  over the compact unit sphere is attained and strictly positive.")
        check("S13", "PROVED: A_0 is a SHARP first-order local minimum of F on Gr(3,6)",
              rk == 9 and strict)

        # numeric value of kappa (LP over the sphere -> solve numerically, exact
        # sign already certified above)
        import numpy as np
        Lf = np.array([[float(sp.N(Lmat[i, j], 40)) for j in range(9)]
                       for i in range(Lmat.rows)])
        # kappa = min_{|z|=1} max_T (L z)_T ; solve by SLSQP from many starts
        from scipy.optimize import minimize
        best = None
        rng = np.random.default_rng(20260802)
        for _ in range(400):
            z0 = rng.standard_normal(9)
            z0 /= np.linalg.norm(z0)
            f = lambda z: float(np.max(Lf @ (z / np.linalg.norm(z))))
            r = minimize(f, z0, method="Nelder-Mead",
                         options=dict(maxiter=4000, maxfev=4000,
                                      xatol=1e-13, fatol=1e-16))
            if best is None or r.fun < best:
                best = float(r.fun)
        print(f"\n  numeric kappa (min over unit sphere, 400 starts) = {best!r}")
        print(f"  (the SIGN is proved exactly above; this value is numeric only)")
        # scale note: |Pdot|_F^2 = 2|z|^2 by S8
        print(f"  in Frobenius normalization |Pdot|_F = 1: kappa_F = {best/sp.sqrt(2)!s}"
              f" ~ {best/2**0.5!r}")

print("\n" + "=" * 78)
fails = [(t, n, d) for t, n, o, d in results if not o]
print(f"SUMMARY: {len(results) - len(fails)}/{len(results)} checks passed")
if fails:
    print("FAILURES:")
    for t, n, d in fails:
        print(f"   - {t}: {n}   {d}")
else:
    print("SHARP MINIMUM CERTIFIED EXACTLY at the Nesterenko extremal.")
    print("This is a PROVED upgrade of Reformulation R part (a) at A_0.")
print("=" * 78)
sys.exit(1 if fails else 0)
