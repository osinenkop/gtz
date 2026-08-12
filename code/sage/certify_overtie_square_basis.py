#!/usr/bin/env python3
"""Krawczyk certificate for a square basis of an over-tied system.

The full over-tied determinant system has more equations than variables.  This
script selects 10 equations, certifies a unique root of that square subsystem by
Krawczyk, and evaluates every omitted equation on the certified root box.

This is deliberately not reported as a full overdetermined proof unless all
omitted residual intervals collapse to exact zero.  In practice it gives the
right local certificate plus a quantified local-membership gap.
"""
from __future__ import annotations

import argparse
import itertools
import json
import sys
from pathlib import Path

from sage.all import Matrix, RealField, RealIntervalField, vector

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from certify_margin_tie_krawczyk import (  # noqa: E402
    interval_from_center_radius,
    interval_poly,
    interval_subset_strict,
    lower_as_str,
    max_upper,
    min_lower,
    principal_minor_polys,
    shifted_block,
    upper_as_str,
)
from gtz63_semialgebraic import TRIPLES, parse_indices  # noqa: E402
from probe_margin_tie_system import shifted_determinant  # noqa: E402
from refine_margin_tie_root import build_system, eval_matrix, eval_vector, max_abs  # noqa: E402


def parse_rows(text: str) -> tuple[int, ...]:
    rows = tuple(int(part.strip()) for part in text.split(",") if part.strip())
    if len(rows) != 10:
        raise ValueError("expected exactly 10 row indices")
    if len(set(rows)) != 10:
        raise ValueError("row indices must be distinct")
    return rows


def label(prefix: str, idx: int) -> str:
    return f"{prefix}:{idx}:{''.join(str(a) for a in TRIPLES[idx])}"


def spectral_branch_polys(block):
    trace = sum(block[i, i] for i in range(3))
    e2 = sum(
        block[i, i] * block[j, j] - block[i, j] * block[j, i]
        for i, j in itertools.combinations(range(3), 2)
    )
    return trace, e2


def krawczyk_attempt(data, polynomials, jacobian, precision, radius_exp):
    RF = RealField(precision)
    RIF = RealIntervalField(precision)
    center_strings = list(data["z"]) + [data["q"]]
    center_rf = vector(RF, [RF(s) for s in center_strings])
    center = vector(RIF, [RIF(s) for s in center_strings])
    radius = RIF(10) ** (-radius_exp)
    X = vector(RIF, [interval_from_center_radius(RIF, c, radius) for c in center_strings])
    delta = vector(RIF, [RIF(-radius, radius) for _ in center_strings])

    F_center = eval_vector(polynomials, center, RIF)
    J_center = eval_matrix(jacobian, center_rf, RF)
    C = J_center.inverse()
    C_interval = Matrix(RIF, C.nrows(), C.ncols(), [RIF(x) for x in C.list()])
    J_box = eval_matrix(jacobian, X, RIF)

    eye = Matrix(RIF, len(center_strings), len(center_strings), 1)
    newton_center = center - C_interval * F_center
    contraction = eye - C_interval * J_box
    K = newton_center + contraction * delta

    included = True
    max_k_radius_ratio = RIF(0)
    coordinate_rows = []
    for i, (ki, xi, ci) in enumerate(zip(K, X, center)):
        ok = interval_subset_strict(ki, xi)
        included = included and ok
        k_radius = max(abs(ki.lower() - ci.center()), abs(ki.upper() - ci.center()))
        ratio = RIF(k_radius) / radius
        if ratio > max_k_radius_ratio:
            max_k_radius_ratio = ratio
        coordinate_rows.append({
            "index": i,
            "K_lower": lower_as_str(ki),
            "K_upper": upper_as_str(ki),
            "X_lower": lower_as_str(xi),
            "X_upper": upper_as_str(xi),
            "strict_subset": bool(ok),
            "radius_ratio_upper": str(ratio),
        })

    return {
        "precision": precision,
        "radius_exp": radius_exp,
        "radius": str(radius),
        "included": bool(included),
        "max_k_radius_ratio": str(max_k_radius_ratio),
        "F_center_inf": str(max_abs(F_center)),
        "coordinates": coordinate_rows,
        "root_box": K,
        "RIF": RIF,
    }


def branch_and_residual_checks(data, chart, RIF, X, all_polynomials, selected_rows, labels):
    active = tuple(data["active"])
    ties = tuple(data["ties"])
    selected = set(selected_rows)
    q = dict(zip(chart.ring.variable_names(), chart.ring.gens()))["q"]

    d_int = interval_poly(chart.d, X, RIF)
    q_int = X[-1]
    margin_int = (q_int - 1) / 6

    active_minor_intervals = []
    tie_minor_intervals = []
    active_trace_e2_intervals = []
    tie_trace_e2_intervals = []
    for idx in active:
        block = shifted_block(chart, idx, chart.ring.one())
        for poly in principal_minor_polys(block):
            active_minor_intervals.append(interval_poly(poly, X, RIF))
        for poly in spectral_branch_polys(block):
            active_trace_e2_intervals.append(interval_poly(poly, X, RIF))
    for idx in ties:
        block = shifted_block(chart, idx, q)
        for poly in principal_minor_polys(block):
            tie_minor_intervals.append(interval_poly(poly, X, RIF))
        for poly in spectral_branch_polys(block):
            tie_trace_e2_intervals.append(interval_poly(poly, X, RIF))

    strict_outside = []
    for idx in range(len(TRIPLES)):
        if idx in active or idx in ties:
            continue
        strict_outside.append((idx, interval_poly(shifted_determinant(chart, idx, q), X, RIF)))

    omitted = []
    for i, poly in enumerate(all_polynomials):
        if i in selected:
            continue
        val = interval_poly(poly, X, RIF)
        omitted.append({
            "row": i,
            "label": labels[i],
            "lower": lower_as_str(val),
            "upper": upper_as_str(val),
            "contains_zero": bool(val.lower() <= 0 <= val.upper()),
            "abs_upper": str(max(abs(val.lower()), abs(val.upper()))),
        })

    min_active = min_lower(active_minor_intervals)
    min_tie = min_lower(tie_minor_intervals)
    min_active_trace_e2 = min_lower(active_trace_e2_intervals)
    min_tie_trace_e2 = min_lower(tie_trace_e2_intervals)
    max_strict_outside = max_upper(strict_outside)
    max_omitted_abs = max((RIF(row["abs_upper"]).upper() for row in omitted), default=RIF(0))
    omitted_all_contain_zero = all(row["contains_zero"] for row in omitted)
    return {
        "chart_d_lower": lower_as_str(d_int),
        "chart_d_upper": upper_as_str(d_int),
        "q_lower": lower_as_str(q_int),
        "q_upper": upper_as_str(q_int),
        "margin_lower": lower_as_str(margin_int),
        "margin_upper": upper_as_str(margin_int),
        "active_minor_min_lower": str(min_active),
        "tie_minor_min_lower": str(min_tie),
        "active_trace_e2_min_lower": str(min_active_trace_e2),
        "tie_trace_e2_min_lower": str(min_tie_trace_e2),
        "strict_outside_shifted_det_max_upper": str(max_strict_outside),
        "chart_d_positive": bool(d_int.lower() > 0),
        "q_gt_one": bool(q_int.lower() > 1),
        "active_minors_positive": bool(min_active is not None and min_active > 0),
        "tie_minors_positive": bool(min_tie is not None and min_tie > 0),
        "active_trace_e2_positive": bool(min_active_trace_e2 is not None and min_active_trace_e2 > 0),
        "tie_trace_e2_positive": bool(min_tie_trace_e2 is not None and min_tie_trace_e2 > 0),
        "strict_outside_shifted_dets_negative": bool(max_strict_outside is None or max_strict_outside < 0),
        "omitted_residuals": omitted,
        "omitted_residuals_all_contain_zero": bool(omitted_all_contain_zero),
        "omitted_residual_abs_max_upper": str(max_omitted_abs),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--rows", required=True, help="10 comma-separated equation row indices")
    parser.add_argument("--precision", type=int, default=800)
    parser.add_argument("--radius-exp-start", type=int, default=30)
    parser.add_argument("--radius-exp-stop", type=int, default=140)
    parser.add_argument("--radius-exp-step", type=int, default=10)
    parser.add_argument("--order", default="degrevlex")
    parser.add_argument("--out", default="")
    args = parser.parse_args()

    data = json.loads(Path(args.input).read_text())
    active = tuple(data["active"])
    ties = tuple(data["ties"])
    rows = parse_rows(args.rows)
    ring, chart, variables, all_polynomials, all_jacobian = build_system(active, ties, args.order)
    if max(rows) >= len(all_polynomials) or min(rows) < 0:
        raise SystemExit(f"row indices out of range for {len(all_polynomials)} equations")
    selected_polynomials = [all_polynomials[i] for i in rows]
    selected_jacobian = [all_jacobian[i] for i in rows]
    labels = [label("active", i) for i in active] + [label("tie", i) for i in ties]

    attempts = []
    winner = None
    for radius_exp in range(args.radius_exp_start, args.radius_exp_stop + 1, args.radius_exp_step):
        attempt = krawczyk_attempt(data, selected_polynomials, selected_jacobian, args.precision, radius_exp)
        root_box = attempt.pop("root_box")
        RIF = attempt.pop("RIF")
        attempts.append(attempt)
        if attempt["included"]:
            branch = branch_and_residual_checks(
                data, chart, RIF, root_box, all_polynomials, rows, labels
            )
            attempt["branch_and_residual_checks"] = branch
            winner = attempt
            break

    basis_certified = bool(winner is not None)
    full_branch_checks = False
    if winner is not None:
        b = winner["branch_and_residual_checks"]
        full_branch_checks = bool(
            b["chart_d_positive"]
            and b["q_gt_one"]
            and b["active_trace_e2_positive"]
            and b["tie_trace_e2_positive"]
            and b["strict_outside_shifted_dets_negative"]
        )

    payload = {
        "input": args.input,
        "active": list(active),
        "active_triples": [list(TRIPLES[i]) for i in active],
        "ties": list(ties),
        "tie_triples": [list(TRIPLES[i]) for i in ties],
        "selected_rows": list(rows),
        "selected_labels": [labels[i] for i in rows],
        "omitted_rows": [i for i in range(len(all_polynomials)) if i not in set(rows)],
        "omitted_labels": [labels[i] for i in range(len(all_polynomials)) if i not in set(rows)],
        "precision_bits": args.precision,
        "attempts": attempts,
        "basis_certified": basis_certified,
        "full_branch_checks": full_branch_checks,
        "full_system_certified": False,
        "full_system_certified_note": "omitted determinant equalities are only interval-enclosed, not proven exact",
        "winner": winner,
    }
    out = args.out or Path(args.input).with_name(Path(args.input).stem + "_basis_krawczyk.json")
    Path(out).parent.mkdir(parents=True, exist_ok=True)
    Path(out).write_text(json.dumps(payload, indent=1) + "\n")

    print("=" * 78)
    print("OVERTIE SQUARE BASIS KRAWCZYK")
    print(f"input: {args.input}")
    print(f"rows:  {list(rows)}")
    print(f"labels: {payload['selected_labels']}")
    if winner is None:
        print("BASIS CERTIFIED: no")
    else:
        b = winner["branch_and_residual_checks"]
        print("BASIS CERTIFIED: yes")
        print(f"radius: 1e-{winner['radius_exp']}")
        print(f"K radius ratio: {winner['max_k_radius_ratio']}")
        print(f"margin interval: [{b['margin_lower']}, {b['margin_upper']}]")
        print(f"omitted residual abs max upper: {b['omitted_residual_abs_max_upper']}")
        print(f"omitted residuals contain zero: {b['omitted_residuals_all_contain_zero']}")
        print("spectral branch checks:",
              b["chart_d_positive"],
              b["q_gt_one"],
              b["active_trace_e2_positive"],
              b["tie_trace_e2_positive"],
              b["strict_outside_shifted_dets_negative"])
    print(f"wrote {out}")
    print("=" * 78)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
