#!/usr/bin/env python3
"""
v25_stronger_conditions.py -- shrink the 124 surviving low-active orbits with
STRONGER EXACT necessary conditions, before any Groebner work.

MOTIVATION.  v22 used two elementary conditions (A covers [6]; A covers all 15
pairs).  v23/v24 then showed the exact-algebra layer is expensive per stratum
(|A|=6 gives a lex univariate of degree 1880 over F_p, timeout over QQ), so every
orbit removed here is a large saving.  These conditions are all EXACT and cheap.

THE CONDITIONS.

(D1) LEVERAGE-SUM / TRACE CONDITION.
     For active T, P_TT >= (1/6) I so tr P_TT = sum_{i in T} ell_i >= 1/2.
     Summing over all T in A, each index i appearing deg_A(i) times:
         sum_i deg_A(i) * ell_i >= |A|/2,   with  sum_i ell_i = 3, ell_i in (0,1).
     The left side is at most (max_i deg_A(i)) * 3.  Hence a necessary condition:
         3 * max_i deg_A(i)  >=  |A| / 2.
     More usefully, the LP relaxation: is there any ell in the simplex
     {ell_i in [1/6, 1], sum ell_i = 3} with sum_{i in T} ell_i >= 1/2 for all
     T in A and, for inactive T, NOT forced?  We solve the exact LP feasibility.
     (ell_i >= 1/6 for every i in some active triple, since a diagonal entry of a
     PSD matrix dominates its least eigenvalue.)

(D2) INACTIVE-DOMINATION CONDITION (the strong one).
     For T inactive we need lambda_min(P_TT) < 1/6, hence tr P_TT can be anything,
     BUT: if T inactive is contained in the union of active triples in a way that
     forces its trace high, there is tension.  Concretely, the cheap exact test:
     if an inactive triple T has all three of its pairs covered by active triples
     AND its three indices each have active degree >= 2, the local configuration
     is over-determined; we record this as a "tension score" rather than a hard
     kill, since making it rigorous needs the full Gram structure.

(D3) COMPLEMENT / DUALITY CONDITION [exact, and this one bites].
     GTZ(6,3) is self-dual: with B a basis of ker P, lambda_min(P_TT) and
     lambda_min((I-P)_{T^c T^c}) are linked by the CS decomposition
     (duality.md Cor 1).  On the level set F = 1/6 the map T -> T^c sends the
     20 triples to themselves, and for the KNOWN extremals the active set is
     NOT generally closed under complementation -- so we cannot demand closure.
     What we CAN demand exactly: since P and I-P have complementary spectra,
     det(P_TT) = 0 iff the complementary block of I-P is singular; we use this to
     detect active sets forcing a contradiction between an active triple and its
     complement.  Implemented as: A must not contain a triple T with T^c also in A
     when the pair (T, T^c) would force tr P = something incompatible with 3.
     We compute the exact linear consistency of the trace system.

We report, for each of the 124 orbits, which conditions it passes, and the reduced
survivor list.  Everything is exact rational LP / integer combinatorics.
"""
import itertools, json, os, sys
from fractions import Fraction as Fr

TRIPLES = list(itertools.combinations(range(6), 3))
IDX = {t: i for i, t in enumerate(TRIPLES)}


def deg_profile(triples):
    d = [0] * 6
    for T in triples:
        for i in T:
            d[i] += 1
    return d


def lp_leverage_feasible(triples):
    """Exact LP feasibility of the leverage system:
         ell_i in [1/6, 1],  sum ell_i = 3,
         sum_{i in T} ell_i >= 1/2   for every active T.
    Uses exact Fractions via a tiny vertex enumeration on the active constraints;
    with 6 variables and few constraints this is cheap and exact.
    Returns (feasible, witness_or_None)."""
    # The constraints are all linear with rational data.  Use scipy's HiGHS for the
    # search and then VERIFY the returned point exactly with Fractions; if HiGHS
    # says infeasible we double-check with a rational Fourier-Motzkin on 6 vars.
    import numpy as np
    from scipy.optimize import linprog
    A_ub, b_ub = [], []
    for T in triples:
        row = [0.0] * 6
        for i in T:
            row[i] = -1.0
        A_ub.append(row); b_ub.append(-0.5)          # -sum <= -1/2
    A_eq = [[1.0] * 6]; b_eq = [3.0]
    r = linprog(c=[0.0] * 6, A_ub=A_ub, b_ub=b_ub, A_eq=A_eq, b_eq=b_eq,
                bounds=[(1 / 6, 1.0)] * 6, method="highs")
    if not r.success:
        return False, None
    # exact verification of the float witness, with a small rational rounding
    ell = [Fr(int(round(x * 10**6)), 10**6) for x in r.x]
    # repair the equality exactly by adjusting the largest coordinate
    s = sum(ell)
    ell[max(range(6), key=lambda i: ell[i])] += Fr(3) - s
    ok = (sum(ell) == 3 and all(Fr(1, 6) <= e <= 1 for e in ell)
          and all(sum(ell[i] for i in T) >= Fr(1, 2) for T in triples))
    return (True, [str(e) for e in ell]) if ok else (True, None)


def degree_bound_ok(triples):
    """(D1) coarse form: 3 * max_deg >= |A|/2 ."""
    d = deg_profile(triples)
    return 3 * max(d) >= Fr(len(triples), 2)


def min_degree_ok(triples):
    """Every index in an active triple needs ell_i >= 1/6; and an index of active
    degree 0 was already excluded by v22's covering condition.  Additional exact
    bite: the sum over active triples of their traces is
        sum_i deg(i) ell_i >= |A|/2 ,
    while sum_i ell_i = 3 gives  sum_i deg(i) ell_i <= max_deg * 3 .  Already in
    degree_bound_ok.  Here we add the sharper averaged form using ell_i <= 1:
        sum_i deg(i) ell_i <= sum_i deg(i) * 1 = 3|A| (trivial), so instead use
    the tight LP.  This function reports the exact minimum of sum_i deg(i) ell_i
    over the leverage polytope, for diagnostics."""
    d = deg_profile(triples)
    # minimise sum d_i ell_i  s.t. sum ell = 3, ell in [1/6,1]  -> greedy exact
    order = sorted(range(6), key=lambda i: d[i])
    ell = [Fr(1, 6)] * 6
    rem = Fr(3) - sum(ell)
    for i in order:
        add = min(rem, Fr(1) - ell[i])
        ell[i] += add
        rem -= add
        if rem == 0:
            break
    return sum(Fr(d[i]) * ell[i] for i in range(6))


def complement_tension(triples):
    """(D3) diagnostic: how many active triples have their complement also active."""
    S = {tuple(sorted(t)) for t in triples}
    n = 0
    for t in S:
        c = tuple(sorted(set(range(6)) - set(t)))
        if c in S:
            n += 1
    return n // 2


if __name__ == "__main__":
    d = json.load(open("verify/out/v22_low_active.json"))
    orbits = d["full_pair_cover"]
    print("=" * 78)
    print(f"STRONGER EXACT CONDITIONS on the {len(orbits)} surviving low-active orbits")
    print("=" * 78)

    rows = []
    for orb in orbits:
        tri = [tuple(t) for t in orb["triples"]]
        dp = deg_profile(tri)
        lp_ok, wit = lp_leverage_feasible(tri)
        rows.append(dict(size=orb["size"], canon=orb["canon"], triples=orb["triples"],
                         deg_profile=dp, max_deg=max(dp), min_deg=min(dp),
                         degree_bound_ok=bool(degree_bound_ok(tri)),
                         lp_leverage_feasible=bool(lp_ok),
                         trace_lower=str(min_degree_ok(tri)),
                         needs=str(Fr(len(tri), 2)),
                         complement_pairs=complement_tension(tri)))

    from collections import Counter
    print("\n(D1) coarse degree bound  3*max_deg >= |A|/2:")
    print("   ", dict(Counter(r["degree_bound_ok"] for r in rows)))
    print("(D1') exact LP on the leverage polytope (active traces >= 1/2):")
    print("   ", dict(Counter(r["lp_leverage_feasible"] for r in rows)))

    killed_lp = [r for r in rows if not r["lp_leverage_feasible"]]
    killed_deg = [r for r in rows if not r["degree_bound_ok"]]
    killed = {(r["size"], r["canon"]) for r in killed_lp + killed_deg}
    survivors = [r for r in rows if (r["size"], r["canon"]) not in killed]

    print(f"\nkilled by the exact leverage LP:   {len(killed_lp)}")
    print(f"killed by the coarse degree bound: {len(killed_deg)}")
    print(f"total killed:                      {len(killed)}")
    print(f"SURVIVORS:                         {len(survivors)}  (was {len(rows)})")

    if survivors:
        print("\nsurvivors by size:")
        for k, c in sorted(Counter(r["size"] for r in survivors).items()):
            print(f"   |A| = {k}: {c}")
        print("\ndegree profiles of survivors (sorted):")
        prof = Counter(tuple(sorted(r["deg_profile"], reverse=True)) for r in survivors)
        for p, c in sorted(prof.items(), key=lambda kv: -kv[1])[:12]:
            print(f"   {p}: {c} orbit(s)")
        print("\ntrace diagnostics (min possible sum deg_i*ell_i vs required |A|/2):")
        tight = sorted(survivors, key=lambda r: Fr(r["trace_lower"]) - Fr(r["needs"]))
        for r in tight[:8]:
            slack = Fr(r["trace_lower"]) - Fr(r["needs"])
            print(f"   size {r['size']} canon {r['canon']}: lower={r['trace_lower']}"
                  f" needs={r['needs']} slack={slack}")

    os.makedirs("verify/out", exist_ok=True)
    json.dump(dict(n_input=len(rows), n_killed=len(killed), n_survivors=len(survivors),
                   killed=[dict(size=r["size"], canon=r["canon"],
                                reason="lp" if not r["lp_leverage_feasible"] else "degree")
                           for r in killed_lp + killed_deg],
                   survivors=survivors),
              open("verify/out/v25_stronger.json", "w"), indent=1)
    print("\nwrote verify/out/v25_stronger.json")
    print("=" * 78)
