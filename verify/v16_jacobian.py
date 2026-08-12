#!/usr/bin/env python3
"""
v16_jacobian.py -- LOCAL DIMENSION of the extremal set E at each known extremal,
via the Jacobian rank of the active-set equations in a 9-variable Grassmann chart.

WHY THIS REPLACES v15's SYMBOLIC ROUTE.  v15 built the equations
det(P_TT - I/6) = 0 symbolically in the chart; after clearing denominators these
are degree-9 polynomials in 9 variables, and sympy's exact Jacobian rank on 13 of
them does not terminate in reasonable time.  But dimension only needs the rank of
the Jacobian AT THE BASE POINT -- a 9-column matrix whose entries we can get to
central-difference accuracy directly from the eigenvalue functions, with no
polynomial expansion at all.

CHART (the standard Grassmann graph chart).  With U an orthonormal basis of
col(P0) and N one of ker(P0), every nearby 3-plane is col(U + N X) for a 3x3 real
matrix X (9 parameters), and
    P(X) = (U + N X) ((U + N X)^T (U + N X))^{-1} (U + N X)^T,   P(0) = P0.
The active equations are g_T(X) := lambda_min(P(X)_TT) - 1/6 = 0 for T in A.
Because lambda_min is a SIMPLE eigenvalue at the base point (verified exactly in
v9/v11), each g_T is analytic near X = 0, so the Jacobian is well defined and its
rows are exactly the functionals of the matrix L.

WHAT THE RANK MEANS.
    rank(J) = 9      => the differentials cut the tangent space to {0}
                     => the stratum is 0-dimensional => P0 is ISOLATED in E.
    rank(J) = r < 9  => tangent space to the stratum has dim >= 9 - r > 0, i.e. a
                        CANDIDATE positive-dimensional family (exact follow-up).

This is an INDEPENDENT route to isolation: v9/v11 obtain it from the critical cone
(an inequality / LP argument), this obtains it from a rank condition (an equality /
linear-algebra argument).  Agreement between the two is a genuine cross-check, not
a restatement.

Numerical (central differences), but a rank verdict is reported only alongside the
spectral gap ratio, so an ambiguous rank is visible rather than hidden.  Tagged
NUMERICALLY SUPPORTED.
"""
import itertools, json, os, sys
import numpy as np

TRIPLES = list(itertools.combinations(range(6), 3))
TH = 1 / 6


def PX(U, N, x):
    B = U + N @ x.reshape(3, 3)
    G = B.T @ B
    return B @ np.linalg.solve(G, B.T)


def gvec(U, N, x, active):
    P = PX(U, N, x)
    return np.array([np.linalg.eigvalsh(P[np.ix_(TRIPLES[i], TRIPLES[i])])[0] - TH
                     for i in active])


def jacobian(U, N, active, h=1e-6):
    J = np.zeros((len(active), 9))
    for j in range(9):
        e = np.zeros(9); e[j] = h
        J[:, j] = (gvec(U, N, e, active) - gvec(U, N, -e, active)) / (2 * h)
    return J


def analyze(P0, label):
    w, V = np.linalg.eigh(P0)
    U = V[:, np.argsort(-w)[:3]]
    N = V[:, np.argsort(-w)[3:]]
    lam = np.array([np.linalg.eigvalsh(P0[np.ix_(T, T)])[0] for T in TRIPLES])
    active = [i for i in range(20) if abs(lam[i] - TH) < 1e-8]

    simple = True
    for i in active:
        ev = np.linalg.eigvalsh(P0[np.ix_(TRIPLES[i], TRIPLES[i])])
        if abs(ev[1] - ev[0]) < 1e-9:
            simple = False

    J = jacobian(U, N, active)
    sv = np.linalg.svd(J, compute_uv=False)
    rank = int(np.sum(sv > sv[0] * 1e-7))
    gap = float(sv[rank - 1] / sv[rank]) if rank < len(sv) else float("inf")
    return dict(label=label, n_active=len(active), simple_lmin=bool(simple),
                singular_values=[float(s) for s in sv], rank_J=rank,
                rank_gap=gap, dim_stratum=int(9 - rank), isolated=bool(rank == 9))


def ttsp_P(tree, w):
    class B:
        def __init__(s):
            s.nv, s.edges = 2, []

        def nvx(s):
            v = s.nv; s.nv += 1; return v

        def build(s, g, a, b):
            if g == ('e',):
                s.edges.append((a, b)); return
            k, *ks = g
            if k == 'P':
                for x in ks:
                    s.build(x, a, b)
            else:
                p = a
                for i, x in enumerate(ks):
                    n = b if i == len(ks) - 1 else s.nvx()
                    s.build(x, p, n); p = n

    b = B(); b.build(tree, 0, 1)
    Bm = np.zeros((b.nv - 1, 6))
    for j, (x, y) in enumerate(b.edges):
        if x != 0:
            Bm[x - 1, j] += 1
        if y != 0:
            Bm[y - 1, j] -= 1
    Mt = np.diag([float(t) ** 0.5 for t in w]) @ Bm.T
    return Mt @ np.linalg.inv(Mt.T @ Mt) @ Mt.T


E = ('e',)
CASES = [
    ("P(S(e,e,e),e,e,e)", ('P', ('S', E, E, E), E, E, E), [1, 1, 1, 5/9, 5/9, 5/9]),
    ("P(S(e,e),S(e,e),e,e)", ('P', ('S', E, E), ('S', E, E), E, E), [1, 1, 1, 1, 5/8, 5/8]),
    ("P(S(P(e,e,e),e),S(e,e))", ('P', ('S', ('P', E, E, E), E), ('S', E, E)),
     [1, 1, 1, 9/5, 9/5, 9/5]),
    ("P(S(P(S(e,e),e,e),e),e)", ('P', ('S', ('P', ('S', E, E), E, E), E), E),
     [1, 1, 5/8, 5/8, 1, 1]),
    ("P(S(P(e,e),P(e,e)),S(e,e))", ('P', ('S', ('P', E, E), ('P', E, E)), ('S', E, E)),
     [1, 1, 1, 1, 8/5, 8/5]),
    ("P(S(P(e,e),e),S(P(e,e),e))", ('P', ('S', ('P', E, E), E), ('S', ('P', E, E), E)),
     [1, 1, 8/5, 1, 1, 8/5]),
]

if __name__ == "__main__":
    out = []
    print("=" * 78)
    print("LOCAL DIMENSION of E via Jacobian rank in a 9-variable Grassmann chart")
    print("  rank(J) = 9  =>  stratum 0-dimensional  =>  extremal ISOLATED")
    print("=" * 78)
    for nm, tree, w in CASES:
        r = analyze(ttsp_P(tree, w), nm)
        out.append(r)
        print(f"\n{nm}")
        print(f"  |A|={r['n_active']}  simple_lmin={r['simple_lmin']}")
        print(f"  sing.values: {np.array2string(np.array(r['singular_values']), precision=4)}")
        print(f"  rank(J) = {r['rank_J']}/9   (gap ratio {r['rank_gap']:.3e})")
        print(f"  => dim stratum = {r['dim_stratum']}    ISOLATED: {r['isolated']}")

    p7 = "verify/data/P514_seventh.npy"
    if os.path.exists(p7):
        r = analyze(np.load(p7), "OUT-OF-FAMILY (5/14,9/14)")
        out.append(r)
        print(f"\nOUT-OF-FAMILY (5/14,9/14)")
        print(f"  |A|={r['n_active']}  simple_lmin={r['simple_lmin']}")
        print(f"  sing.values: {np.array2string(np.array(r['singular_values']), precision=4)}")
        print(f"  rank(J) = {r['rank_J']}/9   (gap ratio {r['rank_gap']:.3e})")
        print(f"  => dim stratum = {r['dim_stratum']}    ISOLATED: {r['isolated']}")

    print("\n" + "=" * 78)
    iso = [r for r in out if r["isolated"]]
    print(f"isolated (rank 9): {len(iso)} / {len(out)}")
    if len(iso) == len(out):
        print("Every known extremal has a 0-DIMENSIONAL stratum: no curve of")
        print("extremals passes through any of them.  This AGREES with the")
        print("critical-cone certificates of v9/v11 by a completely independent")
        print("argument (rank/equality vs cone/LP) -- a genuine cross-check.")
    else:
        print("*** SOME STRATUM HAS POSITIVE DIMENSION -- certify exactly ***")
    json.dump(out, open("verify/out/v16_jacobian.json", "w"), indent=1)
    print("wrote verify/out/v16_jacobian.json")
    print("=" * 78)
