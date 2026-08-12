#!/usr/bin/env python3
"""Interval Krawczyk certification for refined active-plus-tie roots.

Input is a JSON file produced by ``refine_margin_tie_root.py``.  The script
certifies, in interval arithmetic, that the square active-plus-tie determinant
system has a unique root in a small box around the refined point.  It also checks
the spectral branch inequalities needed by the low-active margin argument:

* the chart determinant d is positive;
* q > 1, so the outside tie level is strictly above 1/6;
* active blocks 6P_TT-I and equation-tie blocks 6P_TT-qI have positive 1x1
  and 2x2 principal minors at the certified root branch;
* every other outside triple has det(6P_TT-qI) < 0, hence lambda_min(P_TT)<q/6.

For over-tied numerical points one may pass ``--nonstrict-outside-indices``.
Those outside triples are not part of the square system, but the checker will
try to prove the weaker valid inequality det(6P_TT-qI) <= 0 for them.  If this
fails, no branch certificate is claimed.

The branch checks use interval enclosures on the whole Krawczyk box.  They are
stronger than needed at the root, but convenient and auditable.
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

from gtz63_semialgebraic import TRIPLES, parse_indices  # noqa: E402
from probe_margin_tie_system import shifted_determinant  # noqa: E402
from refine_margin_tie_root import build_system, eval_matrix, eval_vector, max_abs  # noqa: E402


def interval_from_center_radius(RIF, center, radius):
    return RIF(center) + RIF(-radius, radius)


def interval_subset_strict(inner, outer):
    return inner.lower() > outer.lower() and inner.upper() < outer.upper()


def interval_poly(poly, x, RIF):
    return RIF(poly(*x))


def shifted_block(chart, triple_index, theta):
    triple = TRIPLES[triple_index]
    ring = chart.ring
    return Matrix(
        ring,
        3,
        3,
        lambda a, b: 6 * chart.n[triple[a], triple[b]]
        - (theta * chart.d if a == b else ring.zero()),
    )


def principal_minor_polys(block):
    out = []
    for i in range(3):
        out.append(block[i, i])
    for i, j in itertools.combinations(range(3), 2):
        out.append(block[i, i] * block[j, j] - block[i, j] * block[j, i])
    return out


def lower_as_str(x):
    return str(x.lower())


def upper_as_str(x):
    return str(x.upper())


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

    coordinate_rows = []
    included = True
    max_k_radius_ratio = RIF(0)
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
        "box": X,
        "root_box": K,
        "RIF": RIF,
    }


def min_lower(intervals):
    if not intervals:
        return None
    return min(x.lower() for x in intervals)


def max_upper(indexed_intervals):
    if not indexed_intervals:
        return None
    return max(x.upper() for _, x in indexed_intervals)


def branch_checks(data, chart, RIF, X, nonstrict_outside=()):
    active = set(data["active"])
    equation_ties = set(data["ties"])
    nonstrict_outside = set(nonstrict_outside)
    bad_nonstrict = sorted(nonstrict_outside & (active | equation_ties))
    if bad_nonstrict:
        raise ValueError(f"nonstrict outside indices overlap active/equation ties: {bad_nonstrict}")
    q = dict(zip(chart.ring.variable_names(), chart.ring.gens()))["q"]

    d_int = interval_poly(chart.d, X, RIF)
    q_int = X[-1]
    margin_int = (q_int - 1) / 6

    active_minor_intervals = []
    equation_tie_minor_intervals = []
    for idx in sorted(active):
        block = shifted_block(chart, idx, chart.ring.one())
        for poly in principal_minor_polys(block):
            active_minor_intervals.append(interval_poly(poly, X, RIF))
    for idx in sorted(equation_ties):
        block = shifted_block(chart, idx, q)
        for poly in principal_minor_polys(block):
            equation_tie_minor_intervals.append(interval_poly(poly, X, RIF))

    outside_det_intervals = []
    nonstrict_det_intervals = []
    for idx in range(len(TRIPLES)):
        if idx in active or idx in equation_ties:
            continue
        val = interval_poly(shifted_determinant(chart, idx, q), X, RIF)
        if idx in nonstrict_outside:
            nonstrict_det_intervals.append((idx, val))
        else:
            outside_det_intervals.append((idx, val))

    min_active_minor = min_lower(active_minor_intervals)
    min_equation_tie_minor = min_lower(equation_tie_minor_intervals)
    max_outside_det = max_upper(outside_det_intervals)
    max_nonstrict_det = max_upper(nonstrict_det_intervals)
    outside_negative = bool(max_outside_det is None or max_outside_det < 0)
    nonstrict_nonpositive = bool(max_nonstrict_det is None or max_nonstrict_det <= 0)
    active_minors_positive = bool(min_active_minor is not None and min_active_minor > 0)
    equation_tie_minors_positive = bool(
        min_equation_tie_minor is not None and min_equation_tie_minor > 0
    )

    result = {
        "chart_d_lower": lower_as_str(d_int),
        "chart_d_upper": upper_as_str(d_int),
        "q_lower": lower_as_str(q_int),
        "q_upper": upper_as_str(q_int),
        "margin_lower": lower_as_str(margin_int),
        "margin_upper": upper_as_str(margin_int),
        "active_minor_min_lower": str(min_active_minor),
        "equation_tie_minor_min_lower": str(min_equation_tie_minor),
        "tie_minor_min_lower": str(min_equation_tie_minor),
        "outside_shifted_det_max_upper": str(max_outside_det),
        "nonstrict_outside_shifted_det_max_upper": str(max_nonstrict_det),
        "chart_d_positive": bool(d_int.lower() > 0),
        "q_gt_one": bool(q_int.lower() > 1),
        "active_minors_positive": active_minors_positive,
        "equation_tie_minors_positive": equation_tie_minors_positive,
        "tie_minors_positive": equation_tie_minors_positive,
        "outside_shifted_dets_negative": outside_negative,
        "nonstrict_outside_shifted_dets_nonpositive": nonstrict_nonpositive,
        "strict_outside_indices": [idx for idx, _ in outside_det_intervals],
        "nonstrict_outside_indices": sorted(nonstrict_outside),
        "outside_shifted_dets": [
            {
                "index": idx,
                "triple": list(TRIPLES[idx]),
                "lower": lower_as_str(val),
                "upper": upper_as_str(val),
            }
            for idx, val in outside_det_intervals
        ],
        "nonstrict_outside_shifted_dets": [
            {
                "index": idx,
                "triple": list(TRIPLES[idx]),
                "lower": lower_as_str(val),
                "upper": upper_as_str(val),
            }
            for idx, val in nonstrict_det_intervals
        ],
    }
    result["all_branch_checks"] = bool(
        result["chart_d_positive"]
        and result["q_gt_one"]
        and result["active_minors_positive"]
        and result["equation_tie_minors_positive"]
        and result["outside_shifted_dets_negative"]
        and result["nonstrict_outside_shifted_dets_nonpositive"]
    )
    return result


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--precision", type=int, default=700)
    parser.add_argument("--radius-exp-start", type=int, default=30)
    parser.add_argument("--radius-exp-stop", type=int, default=120)
    parser.add_argument("--radius-exp-step", type=int, default=10)
    parser.add_argument("--order", default="degrevlex")
    parser.add_argument(
        "--nonstrict-outside-indices",
        default="",
        help="outside triples allowed to satisfy det(6P_TT-qI)<=0 instead of <0",
    )
    parser.add_argument("--out", default="")
    args = parser.parse_args()

    data = json.loads(Path(args.input).read_text())
    active = tuple(data["active"])
    ties = tuple(data["ties"])
    nonstrict_outside = parse_indices(args.nonstrict_outside_indices)
    if len(active) + len(ties) != 10:
        raise SystemExit("Krawczyk certificate expects a square 10-equation system")

    ring, chart, variables, polynomials, jacobian = build_system(active, ties, args.order)
    attempts = []
    winner = None
    for radius_exp in range(args.radius_exp_start, args.radius_exp_stop + 1, args.radius_exp_step):
        attempt = krawczyk_attempt(data, polynomials, jacobian, args.precision, radius_exp)
        box = attempt.pop("box")
        root_box = attempt.pop("root_box")
        RIF = attempt.pop("RIF")
        attempts.append(attempt)
        if attempt["included"]:
            branch = branch_checks(data, chart, RIF, root_box, nonstrict_outside=nonstrict_outside)
            attempt["branch_checks"] = branch
            attempt["branch_certified"] = bool(branch["all_branch_checks"])
            winner = attempt
            break

    branch_certified = bool(winner is not None and winner.get("branch_certified"))
    payload = {
        "input": args.input,
        "active": list(active),
        "active_triples": [list(TRIPLES[i]) for i in active],
        "ties": list(ties),
        "tie_triples": [list(TRIPLES[i]) for i in ties],
        "equation_ties": list(ties),
        "equation_tie_triples": [list(TRIPLES[i]) for i in ties],
        "nonstrict_outside_indices": list(nonstrict_outside),
        "nonstrict_outside_triples": [list(TRIPLES[i]) for i in nonstrict_outside],
        "precision_bits": args.precision,
        "attempts": attempts,
        "krawczyk_certified": bool(winner is not None),
        "branch_certified": branch_certified,
        "certified": branch_certified,
        "winner": winner,
    }
    out = args.out or Path(args.input).with_name(Path(args.input).stem + "_krawczyk.json")
    Path(out).parent.mkdir(parents=True, exist_ok=True)
    Path(out).write_text(json.dumps(payload, indent=1) + "\n")

    print("=" * 78)
    print("MARGIN TIE KRAWCZYK CERTIFICATE")
    print(f"input:  {args.input}")
    print(f"active: {list(active)}")
    print(f"ties:   {list(ties)}")
    if nonstrict_outside:
        print(f"nonstrict outside: {list(nonstrict_outside)}")
    if winner is None:
        print("CERTIFIED: no")
        print(f"attempts: {len(attempts)}")
    else:
        b = winner["branch_checks"]
        print("KRAWCZYK: yes")
        print(f"BRANCH CERTIFIED: {'yes' if winner['branch_certified'] else 'no'}")
        print(f"radius: 1e-{winner['radius_exp']}")
        print(f"K radius ratio: {winner['max_k_radius_ratio']}")
        print(f"F(center)_inf: {winner['F_center_inf']}")
        print(f"margin interval: [{b['margin_lower']}, {b['margin_upper']}]")
        print(f"active minor lower: {b['active_minor_min_lower']}")
        print(f"equation-tie minor lower: {b['equation_tie_minor_min_lower']}")
        print(f"strict outside det upper: {b['outside_shifted_det_max_upper']}")
        print(f"nonstrict outside det upper: {b['nonstrict_outside_shifted_det_max_upper']}")
        print("branch checks:",
              b["chart_d_positive"],
              b["q_gt_one"],
              b["active_minors_positive"],
              b["equation_tie_minors_positive"],
              b["outside_shifted_dets_negative"],
              b["nonstrict_outside_shifted_dets_nonpositive"])
    print(f"wrote {out}")
    print("=" * 78)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
