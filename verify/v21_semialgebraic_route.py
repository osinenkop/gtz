#!/usr/bin/env python3
"""
v21_semialgebraic_route.py -- polynomial formulation for the finiteness obstruction.

The finiteness route needs to exclude non-sharp equality points.  In a Grassmann
chart Y=[I;Z], write

    P(Z) = Y (Y^T Y)^(-1) Y^T = N(Z) / d(Z),
    d = det(Y^T Y) > 0.

For an active triple T, lambda_min(P_TT)=1/6 and simplicity/rank two of
P_TT-I/6 means

    M_T(Z) := 6 N_TT(Z) - d(Z) I_3

is PSD of rank two.  The raw determinant det(M_T) has the universal factor d^2,
so on the full-rank chart d>0 the active determinant equation is the saturated
degree-6 quotient det(M_T)/d^2=0.  Algebraic ideal computations should still
encode d != 0, e.g. by adding an inverse variable u and u*d-1=0.  On a cofactor
patch, a nonzero kernel vector w_T(Z) is the cross product of two rows of M_T.
Thus the active gradient
inequality in a tangent witness h is polynomial after clearing the positive
denominator d^2:

    w_T^T [ sum_j h_j (d * dN/dz_j - N * dd/dz_j)_TT ] w_T <= 0.

The non-sharp obstruction is therefore semialgebraic in just 18 variables
(9 chart variables z and 9 tangent-witness variables h), plus cofactor-patch
inequalities.  This script does not solve the system.  It generates the system
shape, audits degrees, and records the finite active-set/cofactor-patch problem
that an external CAD/SOS/Groebner pipeline would need to discharge.
"""
from __future__ import annotations

import argparse
import itertools
import json
import os
import sys
from dataclasses import dataclass

import sympy as sp

sys.path.append(os.path.dirname(__file__))
from v17_active_orbits import TRIPLES, burnside_counts


ROW_PAIRS = [(0, 1), (0, 2), (1, 2)]


@dataclass
class ChartPolys:
    z: list[sp.Symbol]
    h: list[sp.Symbol]
    y: sp.Matrix
    d: sp.Expr
    n: sp.Matrix
    ddir: sp.Expr
    ndir: sp.Matrix


def standard_chart(with_directional: bool = False) -> ChartPolys:
    z = list(sp.symbols("z0:9"))
    h = list(sp.symbols("h0:9"))
    y = sp.zeros(6, 3)
    for i in range(3):
        y[i, i] = 1
    for r in range(3):
        for c in range(3):
            y[3 + r, c] = z[3 * r + c]
    gram = y.T * y
    d = sp.factor(gram.det())
    n = y * gram.adjugate() * y.T
    ddir = sp.Integer(0)
    ndir = sp.zeros(6, 6)
    if with_directional:
        ddir = sum(h[j] * sp.diff(d, z[j]) for j in range(9))
        for j in range(9):
            nd = n.diff(z[j])
            ndir += h[j] * (d * nd - n * sp.diff(d, z[j]))
    return ChartPolys(z=z, h=h, y=y, d=d, n=n, ddir=ddir, ndir=ndir)


def principal_minors_1_2(mat: sp.Matrix) -> list[sp.Expr]:
    minors = []
    for i in range(3):
        minors.append(mat[i, i])
    for i, j in itertools.combinations(range(3), 2):
        minors.append(mat.extract([i, j], [i, j]).det())
    return minors


def cross_row_kernel(mat: sp.Matrix, pair: tuple[int, int]) -> sp.Matrix:
    r1 = [mat[pair[0], j] for j in range(3)]
    r2 = [mat[pair[1], j] for j in range(3)]
    return sp.Matrix([
        r1[1] * r2[2] - r1[2] * r2[1],
        r1[2] * r2[0] - r1[0] * r2[2],
        r1[0] * r2[1] - r1[1] * r2[0],
    ])


def total_degree(expr: sp.Expr, vars_: list[sp.Symbol]) -> int:
    return sp.Poly(sp.expand(expr), *vars_).total_degree()


def active_block(chart: ChartPolys, triple: tuple[int, int, int]) -> sp.Matrix:
    block = sp.Matrix(3, 3, lambda a, b: chart.n[triple[a], triple[b]])
    return 6 * block - chart.d * sp.eye(3)


def derivative_block(chart: ChartPolys, triple: tuple[int, int, int]) -> sp.Matrix:
    return sp.Matrix(3, 3, lambda a, b: chart.ndir[triple[a], triple[b]])


def audit_one_triple(triple_index: int, row_pair: tuple[int, int], exact_expansion: bool = False) -> dict:
    chart = standard_chart(with_directional=exact_expansion)
    triple = TRIPLES[triple_index]
    m = active_block(chart, triple)
    w = cross_row_kernel(m, row_pair)
    vars_z = chart.z
    result = {
        "triple_index": triple_index,
        "triple": list(triple),
        "row_pair": list(row_pair),
        "degree_d": total_degree(chart.d, vars_z),
        "degree_N_entry_max": max(total_degree(chart.n[i, j], vars_z) for i in range(6) for j in range(6)),
        "structural_degree_upper_bounds": {
            "active_det_raw": 18,
            "active_det_saturated": 6,
            "psd_2x2_minor": 12,
            "cofactor_patch_norm": 24,
            "cleared_directional_ineq": 35,
        },
        "exact_expansion": exact_expansion,
    }
    if exact_expansion:
        dblock = derivative_block(chart, triple)
        patch_norm = sum(x * x for x in w)
        det_eq = m.det()
        sat_det_poly, sat_det_rem = sp.div(
            sp.Poly(sp.expand(det_eq), *vars_z),
            sp.Poly(sp.expand(chart.d ** 2), *vars_z),
        )
        if sat_det_rem.as_expr() != 0:
            raise ArithmeticError("active determinant is not divisible by d^2")
        sat_det = sat_det_poly.as_expr()
        deriv_ineq = (w.T * dblock * w)[0]
        minors = principal_minors_1_2(m)
        vars_zh = chart.z + chart.h
        result.update({
            "degree_active_det_raw": total_degree(det_eq, vars_z),
            "degree_active_det_saturated": total_degree(sat_det, vars_z),
            "degree_psd_minor_max": max(total_degree(x, vars_z) for x in minors),
            "degree_patch_norm": total_degree(patch_norm, vars_z),
            "degree_derivative_ineq": total_degree(deriv_ineq, vars_zh),
            "n_terms_active_det_raw": len(sp.Poly(sp.expand(det_eq), *vars_z).terms()),
            "n_terms_active_det_saturated": len(sp.Poly(sp.expand(sat_det), *vars_z).terms()),
            "n_terms_derivative_ineq": len(sp.Poly(sp.expand(deriv_ineq), *vars_zh).terms()),
        })
    return result


def system_size(active_size: int, include_inactive_flags: bool = False) -> dict:
    """Counts for one cofactor-patch branch of a simple-active obstruction."""
    return {
        "active_size": active_size,
        "variables": {
            "chart_z": 9,
            "tangent_h": 9,
            "total": 18,
        },
        "equalities": {
            "active_det_M_T_saturated": active_size,
            "tangent_normalization_sum_h2_eq_1": 1,
            "total": active_size + 1,
        },
        "inequalities": {
            "chart_denominator_d_gt_0": 1,
            "cofactor_patch_norm_gt_0": active_size,
            "active_psd_1x1_2x2_minors_ge_0": 6 * active_size,
            "nonsharp_directional_derivatives_le_0": active_size,
            "inactive_F_le_1_6_disjunctive": "not expanded" if include_inactive_flags else "relaxed away",
            "total_without_inactive_disjunctions": 1 + 8 * active_size,
        },
        "cofactor_patch_branches_per_active_set": 3 ** active_size,
        "note": (
            "Simple-active branch. Multiple least-eigenvalue active blocks need "
            "a separate SDP/subgradient formulation."
        ),
    }


def orbit_work_table() -> list[dict]:
    counts = burnside_counts()
    rows = []
    for size, n_orbits in enumerate(counts):
        if size == 0 or n_orbits == 0:
            continue
        rows.append({
            "active_size": size,
            "active_set_orbits": n_orbits,
            "simple_sharp_possible": size >= 10,
            "simple_nonsharp_automatic": size <= 9,
            "cofactor_patch_branches_per_orbit": 3 ** size,
        })
    return rows


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--triple-index", type=int, default=0)
    parser.add_argument("--row-pair", choices=["01", "02", "12"], default="01")
    parser.add_argument("--sizes", default="8,9,10,12,13")
    parser.add_argument("--exact-expansion", action="store_true")
    parser.add_argument("--output", default="verify/out/v21_semialgebraic_route.json")
    args = parser.parse_args()

    pair = {"01": (0, 1), "02": (0, 2), "12": (1, 2)}[args.row_pair]
    audit = audit_one_triple(args.triple_index, pair, exact_expansion=args.exact_expansion)
    sizes = [int(x) for x in args.sizes.split(",") if x.strip()]
    size_rows = [system_size(size) for size in sizes]
    work = orbit_work_table()

    print("=" * 78)
    print("SEMIALGEBRAIC FINITENESS OBSTRUCTION: system shape")
    print("=" * 78)
    print("Cofactor-patch model variables: 9 chart variables z + 9 tangent variables h")
    print("Simple-active sharpness obstruction: det equations + PSD + nonsharp h")
    print()
    print("degree audit for one active triple:")
    for key, value in audit.items():
        print(f"  {key}: {value}")
    if not args.exact_expansion:
        print("  exact term counts skipped; pass --exact-expansion to compute them")
    print()
    print("system sizes:")
    for row in size_rows:
        print(
            f"  |A|={row['active_size']:2d}: vars={row['variables']['total']} "
            f"eqs={row['equalities']['total']} "
            f"ineqs={row['inequalities']['total_without_inactive_disjunctions']} "
            f"patches/orbit=3^{row['active_size']}"
        )
    print()
    print("simple-active automatic nonsharp range: |A| <= 9")
    print("simple-active sharpness can only start at |A| >= 10")

    out = {
        "model": "simple-active cofactor-patch nonsharp obstruction",
        "degree_audit": audit,
        "system_sizes": size_rows,
        "active_orbit_work_table": work,
        "proof_target": (
            "Show no real solution of the equality + PSD + nonsharp witness "
            "system outside the known sharp active-set orbits, with a separate "
            "case split for multiple least-eigenvalue active blocks."
        ),
    }
    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    with open(args.output, "w") as fh:
        json.dump(out, fh, indent=1)
    print(f"\nwrote {args.output}")
    print("=" * 78)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
