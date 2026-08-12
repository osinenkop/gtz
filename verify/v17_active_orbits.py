#!/usr/bin/env python3
"""
v17_active_orbits.py -- symmetry-reduced active-set bookkeeping for GTZ(6,3).

This is a combinatorial sizing pass for active-set/KKT routes.  It does NOT prove
GTZ(6,3), finiteness of the extremal set, or infeasibility of any stratum.

The 20 triples of [6] are acted on by S_6.  Any exact enumeration by active set
should work modulo this action, not over all 2^20 subsets.  Burnside's lemma gives
the exact number of active-subset orbits by cardinality.  We also record the
active-set orbit of each of the seven certified extremals so future scripts can
avoid re-discovering this bookkeeping.
"""
import itertools
import json
import os
from collections import defaultdict

import numpy as np

TRIPLES = list(itertools.combinations(range(6), 3))
TRIPLE_INDEX = {t: i for i, t in enumerate(TRIPLES)}
TH = 1 / 6
E = ("e",)

TTSP_CASES = [
    ("P(S(e,e,e),e,e,e)", ("P", ("S", E, E, E), E, E, E), [1, 1, 1, 5 / 9, 5 / 9, 5 / 9]),
    ("P(S(e,e),S(e,e),e,e)", ("P", ("S", E, E), ("S", E, E), E, E), [1, 1, 1, 1, 5 / 8, 5 / 8]),
    ("P(S(P(e,e,e),e),S(e,e))", ("P", ("S", ("P", E, E, E), E), ("S", E, E)),
     [1, 1, 1, 9 / 5, 9 / 5, 9 / 5]),
    ("P(S(P(S(e,e),e,e),e),e)", ("P", ("S", ("P", ("S", E, E), E, E), E), E),
     [1, 1, 5 / 8, 5 / 8, 1, 1]),
    ("P(S(P(e,e),P(e,e)),S(e,e))", ("P", ("S", ("P", E, E), ("P", E, E)), ("S", E, E)),
     [1, 1, 1, 1, 8 / 5, 8 / 5]),
    ("P(S(P(e,e),e),S(P(e,e),e))", ("P", ("S", ("P", E, E), E), ("S", ("P", E, E), E)),
     [1, 1, 8 / 5, 1, 1, 8 / 5]),
]


def s6_actions():
    actions = []
    for p in itertools.permutations(range(6)):
        perm = []
        for triple in TRIPLES:
            image = tuple(sorted(p[i] for i in triple))
            perm.append(TRIPLE_INDEX[image])
        actions.append(tuple(perm))
    return actions


S6_ACTIONS = s6_actions()


def apply_action(mask, action):
    out = 0
    for i, j in enumerate(action):
        if mask & (1 << i):
            out |= 1 << j
    return out


def canonical_mask(mask):
    return min(apply_action(mask, action) for action in S6_ACTIONS)


def orbit_size(mask):
    return len({apply_action(mask, action) for action in S6_ACTIONS})


def cycle_lengths(action):
    seen = [False] * len(TRIPLES)
    lengths = []
    for i in range(len(TRIPLES)):
        if seen[i]:
            continue
        j, n = i, 0
        while not seen[j]:
            seen[j] = True
            n += 1
            j = action[j]
        lengths.append(n)
    return lengths


def burnside_counts():
    counts = [0] * (len(TRIPLES) + 1)
    for action in S6_ACTIONS:
        poly = [1] + [0] * len(TRIPLES)
        for length in cycle_lengths(action):
            for degree in range(len(TRIPLES) - length, -1, -1):
                if poly[degree]:
                    poly[degree + length] += poly[degree]
        counts = [a + b for a, b in zip(counts, poly)]
    return [c // len(S6_ACTIONS) for c in counts]


def ttsp_projector(tree, weights):
    class Builder:
        def __init__(self):
            self.nv = 2
            self.edges = []

        def new_vertex(self):
            vertex = self.nv
            self.nv += 1
            return vertex

        def build(self, graph, source, sink):
            if graph == E:
                self.edges.append((source, sink))
                return
            kind, *kids = graph
            if kind == "P":
                for kid in kids:
                    self.build(kid, source, sink)
            else:
                prev = source
                for idx, kid in enumerate(kids):
                    nxt = sink if idx == len(kids) - 1 else self.new_vertex()
                    self.build(kid, prev, nxt)
                    prev = nxt

    builder = Builder()
    builder.build(tree, 0, 1)
    incidence = np.zeros((builder.nv - 1, 6))
    for col, (x, y) in enumerate(builder.edges):
        if x != 0:
            incidence[x - 1, col] += 1
        if y != 0:
            incidence[y - 1, col] -= 1
    matrix = np.diag([w ** 0.5 for w in weights]) @ incidence.T
    return matrix @ np.linalg.inv(matrix.T @ matrix) @ matrix.T


def active_mask(projector):
    mask = 0
    values = []
    for i, triple in enumerate(TRIPLES):
        lam = float(np.linalg.eigvalsh(projector[np.ix_(triple, triple)])[0])
        values.append(lam)
        if abs(lam - TH) < 1e-8:
            mask |= 1 << i
    return mask, values


def mask_triples(mask):
    return [TRIPLES[i] for i in range(len(TRIPLES)) if mask & (1 << i)]


def known_extremal_orbits():
    records = []
    for label, tree, weights in TTSP_CASES:
        mask, values = active_mask(ttsp_projector(tree, weights))
        records.append({
            "label": label,
            "active_size": int(mask.bit_count()),
            "orbit_size": orbit_size(mask),
            "canonical_mask": canonical_mask(mask),
            "active_triples": [list(t) for t in mask_triples(mask)],
            "F": max(values),
        })

    seventh = "verify/data/P514_seventh.npy"
    if os.path.exists(seventh):
        mask, values = active_mask(np.load(seventh))
        records.append({
            "label": "OUT-OF-FAMILY (5/14,9/14)",
            "active_size": int(mask.bit_count()),
            "orbit_size": orbit_size(mask),
            "canonical_mask": canonical_mask(mask),
            "active_triples": [list(t) for t in mask_triples(mask)],
            "F": max(values),
        })
    return records


if __name__ == "__main__":
    counts = burnside_counts()
    print("=" * 78)
    print("S_6 ORBITS OF ACTIVE SUBSETS OF THE 20 TRIPLES")
    print("=" * 78)
    for size, count in enumerate(counts):
        if count:
            print(f"  |A|={size:2d}: {count:4d}")
    print("-" * 78)
    print(f"total orbits:        {sum(counts)}")
    print(f"nonempty |A| <= 8:   {sum(counts[1:9])}")
    print(f"|A| >= 9:            {sum(counts[9:])}")

    known = known_extremal_orbits()
    by_orbit = defaultdict(list)
    for record in known:
        by_orbit[record["canonical_mask"]].append(record["label"])

    print("\nknown certified extremals:")
    for record in known:
        print(f"  |A|={record['active_size']:2d}  orbit={record['orbit_size']:3d}  {record['label']}")
    print(f"\ndistinct active-set orbits among known extremals: {len(by_orbit)}")

    os.makedirs("verify/out", exist_ok=True)
    out = {
        "orbit_counts_by_active_size": {str(i): c for i, c in enumerate(counts) if c},
        "total_orbits": sum(counts),
        "nonempty_orbits_active_size_le_8": sum(counts[1:9]),
        "orbits_active_size_ge_9": sum(counts[9:]),
        "known_extremals": known,
        "distinct_known_active_orbits": len(by_orbit),
    }
    with open("verify/out/v17_active_orbits.json", "w") as fh:
        json.dump(out, fh, indent=1)
    print("\nwrote verify/out/v17_active_orbits.json")
    print("=" * 78)
