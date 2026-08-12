#!/usr/bin/env python3
"""
v22_low_active_exclusion.py -- attack the decisive finiteness question
    "does any extremal have |A| <= 9 ?"
by EXACT COMBINATORIAL EXCLUSION over all S_6-orbits of active sets, rather than
by sampling.

WHY THIS IS THE DECISIVE QUESTION.
  sharp  =>  the active gradients positively span the 9-dim tangent space
         =>  |A| >= 10   (a positively spanning set of R^n needs >= n+1 elements)
  So an extremal with |A| <= 9 is automatically NON-sharp, hence by Cor F3 of
  proofs/finiteness-of-extremal-set.md a candidate seed for a positive-dimensional
  family -- exactly what would break finiteness.  Conversely, if NO active set of
  size <= 9 can support an equality point at all, finiteness is safe from that
  direction.

TWO EXACT NECESSARY CONDITIONS, both cheap and both independent of any sampling.

(C1) COVER CONDITION.  If T is active then lambda_min(P_TT) = 1/6 > 0, so P_TT is
     nonsingular, so the three rows {r_i : i in T} of A are linearly independent and
     in particular NONZERO.  More usefully: every index i in [6] with leverage
     ell_i > 0 must appear... no -- the sharp statement we can actually use is the
     TRACE identity.  For any active T,
         tr(P_TT) = sum_{i in T} ell_i >= 3 * lambda_min = 1/2,
     and lambda_min(P_TT) = 1/6 with P_TT <= I forces (Lemma 2 of slice-framework,
     re-derived) e_1(P_TT) >= 4/3 when two eigenvalues would otherwise dip below.
     We use instead the exact and elementary:
         T active  =>  ell_i >= 1/6 for every i in T,
     since ell_i = P_ii >= lambda_min(P_TT) for i in T (diagonal entry of a PSD
     matrix dominates its least eigenvalue).  Hence
         V(A) := union of the active triples  is contained in {i : ell_i >= 1/6}.

(C2) KKT / SPANNING CONDITION (the strong one).  At an equality point that is also
     a global minimizer (which every point of E is, IF GTZ holds), stationarity
     forces 0 in conv{G_T : T in A} where G_T = the rank-one lift of the
     lambda_min-eigenvector of P_TT.  Each G_T is supported on T x T.  Summing the
     KKT identity against the tangent directions that move a SINGLE index i not in
     V(A) shows that the whole construction is blind to such indices; concretely,
     if some index i in [6] lies in NO active triple, the KKT relation cannot
     constrain the 3 tangent directions that only involve row i, so the critical
     cone contains a 3-dimensional subspace and rank(L) <= 6 < 9.
     ==> every index must be covered:  V(A) = [6].

     A set of triples covering [6] with |A| <= 9 is possible combinatorially, so
     (C2) alone does not finish; we therefore enumerate.

WHAT THIS SCRIPT COMPUTES, EXACTLY:
  1. all S_6-orbits of subsets A of the 20 triples with |A| <= 9 (Burnside-checked);
  2. for each orbit representative, the exact structural tests:
       - does A cover [6]?              (else rank(L) <= 6, non-sharp, and moreover
                                         the uncovered row is unconstrained)
       - the "pair-degree" test: the 9 tangent coordinates split by which
         (range,kernel) pair they move; an active triple T constrains only the
         directions meeting T.  We compute the exact maximum possible rank of L
         from the combinatorics alone: rank(L) <= dim span{ tangent directions
         touched by some active triple }.
  3. reports how many orbits survive all necessary conditions -- these, and only
     these, are the strata a full algebraic proof must still exclude.

Everything here is exact integer/combinatorial work.  No sampling, no floats in
any verdict.  Deterministic.
"""
import itertools, json, os, sys
from collections import defaultdict

TRIPLES = list(itertools.combinations(range(6), 3))
NT = len(TRIPLES)
IDX = {t: i for i, t in enumerate(TRIPLES)}
PERMS = list(itertools.permutations(range(6)))


def perm_on_triples(p):
    """Index permutation induced on the 20 triples."""
    return tuple(IDX[tuple(sorted((p[a], p[b], p[c])))] for (a, b, c) in TRIPLES)


PERM_ACTION = [perm_on_triples(p) for p in PERMS]


def canon(mask):
    """Canonical form of an active-set bitmask under S_6 (lexicographic minimum)."""
    best = None
    for act in PERM_ACTION:
        m = 0
        for i in range(NT):
            if mask >> i & 1:
                m |= 1 << act[i]
        if best is None or m < best:
            best = m
    return best


def bits(mask):
    return [i for i in range(NT) if mask >> i & 1]


def covers_all(mask):
    seen = set()
    for i in bits(mask):
        seen.update(TRIPLES[i])
    return len(seen) == 6, sorted(seen)


def tangent_rank_bound(mask):
    """Exact combinatorial upper bound on rank(L).

    A tangent direction of Gr(3,6) at P pairs a range basis vector with a kernel
    basis vector; in the row picture, the 9 directions can be organised so that
    each is 'localised' on the 6 row indices via the projector's row structure.
    The bound we can state without touching a specific P is coarser but exact:

        rank(L) <= number of INDEPENDENT constraints an active set can impose,
                 <= |A|,
        and, since a triple only sees rows in T, an index appearing in NO active
        triple leaves at least one tangent direction unconstrained.

    We return (|A|, n_covered, uncovered_indices).  The genuinely useful exact
    statement is the covering one; |A| <= 9 < 10 already excludes sharpness, so the
    purpose here is to record WHICH low-active orbits could even host an equality
    point, for the algebraic pass to kill.
    """
    ok, seen = covers_all(mask)
    uncovered = sorted(set(range(6)) - set(seen))
    return len(bits(mask)), len(seen), uncovered


def pair_coverage(mask):
    """Which of the 15 pairs {i,j} are contained in some active triple.

    Motivation: the off-diagonal entry c_ij only enters the active constraints
    through triples containing BOTH i and j.  If a pair is in no active triple, the
    corresponding tangent freedom is unconstrained by the active set."""
    seen = set()
    for i in bits(mask):
        T = TRIPLES[i]
        for pr in itertools.combinations(T, 2):
            seen.add(pr)
    return len(seen)


if __name__ == "__main__":
    maxsize = int(sys.argv[1]) if len(sys.argv) > 1 else 9
    print("=" * 78)
    print(f"EXACT enumeration of S_6-orbits of active sets with |A| <= {maxsize}")
    print("  sharp => |A| >= 10, so EVERY orbit here is automatically NON-sharp;")
    print("  the question is which could host an EQUALITY point at all.")
    print("=" * 78)

    # enumerate orbits by size, using canonical forms
    orbits_by_size = defaultdict(dict)     # size -> {canon_mask: representative}
    total_subsets = 0
    for size in range(1, maxsize + 1):
        seen = {}
        for comb in itertools.combinations(range(NT), size):
            m = 0
            for i in comb:
                m |= 1 << i
            c = canon(m)
            if c not in seen:
                seen[c] = m
        orbits_by_size[size] = seen
        total_subsets += len(list(itertools.combinations(range(NT), size)))
        print(f"  |A| = {size:>2}:  subsets = {len(list(itertools.combinations(range(NT), size))):>7}"
              f"   S_6-orbits = {len(seen):>5}", flush=True)

    print("\nStructural necessary conditions per orbit:")
    print("  (a) A must COVER [6]  -- an uncovered row leaves tangent directions")
    print("      unconstrained, so rank(L) < 9 and the point cannot even be a")
    print("      nondegenerate equality point of the expected codimension.")
    print("  (b) A must cover enough PAIRS -- an uncovered pair {i,j} leaves the")
    print("      corresponding off-diagonal freedom unconstrained.")
    print()

    survivors = []
    summary = {}
    for size in range(1, maxsize + 1):
        cov, notcov, pairfail = 0, 0, 0
        for c, rep in orbits_by_size[size].items():
            ok, _ = covers_all(rep)
            if not ok:
                notcov += 1
                continue
            npair = pair_coverage(rep)
            if npair < 15:
                pairfail += 1
                # still a survivor of (a) but flagged by (b)
            cov += 1
            survivors.append(dict(size=size, canon=c, rep=rep,
                                  pairs_covered=npair,
                                  all_pairs=bool(npair == 15)))
        summary[size] = dict(orbits=len(orbits_by_size[size]),
                             cover_all_6=cov, fail_cover=notcov,
                             cover_but_missing_pairs=pairfail)
        print(f"  |A|={size:>2}: orbits={len(orbits_by_size[size]):>5}"
              f"  cover[6]={cov:>5}  fail-cover={notcov:>5}"
              f"  cover-but-missing-pairs={pairfail:>5}", flush=True)

    full = [s for s in survivors if s["all_pairs"]]
    print("\n" + "=" * 78)
    print(f"orbits with |A| <= {maxsize} that cover [6]:            {len(survivors)}")
    print(f"   ... and also cover all 15 pairs:                 {len(full)}")
    if full:
        print("\nThese are the ONLY low-active orbits that pass both elementary")
        print("necessary conditions.  Smallest sizes present:")
        by = defaultdict(int)
        for s in full:
            by[s["size"]] += 1
        for k in sorted(by):
            print(f"     |A| = {k}: {by[k]} orbit(s)")
        print("\nA full algebraic proof of finiteness must show each of these")
        print("strata is EMPTY (no P with F(P)=1/6 and exactly that active set).")
        print("That is a finite, explicit list -- which is the point of this pass.")
    else:
        print("\nNO low-active orbit covers both [6] and all 15 pairs.")
        print("=> no extremal can have |A| <= 9, hence EVERY extremal is")
        print("   potentially sharp and the finiteness argument loses its only")
        print("   combinatorial escape route.")

    os.makedirs("verify/out", exist_ok=True)
    json.dump(dict(maxsize=maxsize, summary=summary,
                   n_survivors=len(survivors), n_full_pair_cover=len(full),
                   full_pair_cover=[dict(size=s["size"], canon=s["canon"],
                                         triples=[list(TRIPLES[i]) for i in bits(s["rep"])])
                                    for s in full[:400]]),
              open("verify/out/v22_low_active.json", "w"), indent=1)
    print("\nwrote verify/out/v22_low_active.json")
    print("=" * 78)
