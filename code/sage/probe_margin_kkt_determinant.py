#!/usr/bin/env python3
"""Probe determinant-gradient KKT systems for low-active margin minimizers.

For a low-active set A and an outside max-tie set B, this builds the square
polynomial system for critical points of

    minimise q
    subject to f_T(z)=0             for T in A,
               h_S(z,q)=0          for S in B,
               det(Y^T Y) != 0.

Here Y=[I;Z] is the standard chart, f_T=det(6*N_TT-d*I)/d^2, and
h_S=det(6*N_SS-q*d*I)/d^2.  Thus q is six times the tied outside eigenvalue,
and the margin is (q-1)/6.

The KKT equations are

    sum a_T grad_z f_T + sum b_S grad_z h_S = 0,
    1 + sum b_S partial_q h_S = 0.

Inequality checks (PSD, ordering, multiplier signs after converting back from
determinant multipliers) are deliberately external.
"""
from __future__ import annotations

import argparse
import json
import shutil
import sys
from pathlib import Path

from sage.all import PolynomialRing, QQ

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from gtz63_semialgebraic import (  # noqa: E402
    TRIPLES,
    SingularTextSystem,
    active_determinant,
    standard_chart,
    text_system_summary,
)
from probe_determinant_ideals import (  # noqa: E402
    random_linear_forms,
    run_singular,
    write_probe_script,
    zero_sum_z_linear_forms,
)
from probe_margin_tie_system import (  # noqa: E402
    resolve_active,
    shifted_determinant,
)


def margin_kkt_text_system(
    active_indices: tuple[int, ...],
    tie_indices: tuple[int, ...],
    order: str,
    multiplier_order: str,
) -> SingularTextSystem:
    multiplier_names = (
        [f"a{i}" for i in range(len(active_indices))]
        + [f"b{i}" for i in range(len(tie_indices))]
    )
    chart_names = [f"z{i}" for i in range(9)] + ["q", "u0"]
    names = (
        multiplier_names + chart_names
        if multiplier_order == "first"
        else chart_names + multiplier_names
    )
    ring = PolynomialRing(QQ, names, order=order)
    gens = dict(zip(ring.variable_names(), ring.gens()))
    chart = standard_chart(ring, include_h=False, inverse_count=1)
    z = chart.z
    q = gens["q"]
    alpha = [gens[f"a{i}"] for i in range(len(active_indices))]
    beta = [gens[f"b{i}"] for i in range(len(tie_indices))]

    active_polys = [active_determinant(chart, idx) for idx in active_indices]
    tie_polys = [shifted_determinant(chart, idx, q) for idx in tie_indices]

    equalities = []
    degrees = []
    for poly in active_polys + tie_polys:
        equalities.append(poly)
        degrees.append(int(poly.total_degree()))

    inv_eq = chart.inv[0] * chart.d - 1
    equalities.append(inv_eq)
    degrees.append(int(inv_eq.total_degree()))

    for zj in z:
        stationarity = ring.zero()
        for coeff, poly in zip(alpha, active_polys):
            stationarity += coeff * poly.derivative(zj)
        for coeff, poly in zip(beta, tie_polys):
            stationarity += coeff * poly.derivative(zj)
        equalities.append(stationarity)
        degrees.append(int(stationarity.total_degree()))

    q_stationarity = ring.one()
    for coeff, poly in zip(beta, tie_polys):
        q_stationarity += coeff * poly.derivative(q)
    equalities.append(q_stationarity)
    degrees.append(int(q_stationarity.total_degree()))

    return SingularTextSystem(
        mode="margin_kkt_det_inverted",
        variables=tuple(str(name) for name in ring.variable_names()),
        active_indices=active_indices,
        row_pairs=(),
        equalities=tuple(str(poly) for poly in equalities),
        nonzero=(),
        equality_degrees=tuple(degrees),
        nonzero_degrees=(),
        metadata={
            "chart": "Y=[I;Z]",
            "description": (
                "determinant-gradient KKT system for minimising outside tied "
                "q subject to active and tie determinant equations"
            ),
            "active_triples": [list(TRIPLES[i]) for i in active_indices],
            "tie_indices": list(tie_indices),
            "tie_triples": [list(TRIPLES[i]) for i in tie_indices],
            "q_meaning": "q = 6 * tied outside eigenvalue; margin = (q-1)/6",
            "alpha_variables": [str(x) for x in alpha],
            "beta_variables": [str(x) for x in beta],
            "ring_order": order,
            "multiplier_order": multiplier_order,
            "invert_d": True,
        },
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--active-indices")
    parser.add_argument("--size", type=int)
    parser.add_argument("--canon", type=int)
    parser.add_argument("--tie-indices", required=True)
    parser.add_argument("--order", default="degrevlex")
    parser.add_argument("--multiplier-order", choices=["last", "first"], default="last")
    parser.add_argument("--methods", default="slimgb")
    parser.add_argument("--timeout", type=int, default=300)
    parser.add_argument("--characteristic", type=int, default=32003)
    parser.add_argument("--singular-bin", default=shutil.which("Singular") or "Singular")
    parser.add_argument("--out-prefix", default="code/sage/out/margin_kkt")
    parser.add_argument("--linear-sections", type=int, default=0)
    parser.add_argument(
        "--linear-section-mode",
        choices=["random", "zero-sum-z"],
        default="random",
    )
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--coefficient-bound", type=int, default=17)
    parser.add_argument("--export-only", action="store_true")
    args = parser.parse_args()

    active = resolve_active(args)
    from gtz63_semialgebraic import parse_indices  # local to reuse resolve_active parser shape
    ties = parse_indices(args.tie_indices)
    if set(active) & set(ties):
        raise SystemExit("tie indices must be outside the active set")

    system = margin_kkt_text_system(active, ties, args.order, args.multiplier_order)
    summary = text_system_summary(system)
    if args.linear_section_mode == "random":
        linear_forms = random_linear_forms(
            system.variables,
            count=args.linear_sections,
            seed=args.seed,
            characteristic=args.characteristic,
            coefficient_bound=args.coefficient_bound,
        )
    else:
        linear_forms = zero_sum_z_linear_forms(
            system.variables,
            count=args.linear_sections,
            seed=args.seed,
            coefficient_bound=args.coefficient_bound,
        )

    out_prefix = Path(args.out_prefix)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    results = {
        "summary": summary,
        "timeout_seconds": args.timeout,
        "characteristic": args.characteristic,
        "linear_sections": args.linear_sections,
        "linear_section_mode": args.linear_section_mode,
        "linear_section_seed": args.seed,
        "linear_forms": list(linear_forms),
        "singular_bin": args.singular_bin,
        "methods": {},
    }

    methods = tuple(x.strip() for x in args.methods.split(",") if x.strip())
    for method in methods:
        script_path = out_prefix.with_name(f"{out_prefix.name}_{method}.sing")
        write_probe_script(
            system,
            script_path,
            method=method,
            order=args.order,
            characteristic=args.characteristic,
            linear_forms=linear_forms,
        )
        if args.export_only:
            results["methods"][method] = {"script": str(script_path), "not_run": True}
            print(f"{method}: exported {script_path}")
            continue
        run = run_singular(script_path, args.singular_bin, args.timeout)
        results["methods"][method] = {"script": str(script_path), **run}
        parsed = run.get("parsed", {})
        dim = parsed.get("dimension", "?") if isinstance(parsed, dict) else "?"
        degree = parsed.get("degree", "?") if isinstance(parsed, dict) else "?"
        status = "timeout" if run["timed_out"] else f"exit {run['exit_code']}"
        print(f"{method}: {status}, elapsed={run['elapsed_seconds']:.1f}s, "
              f"dimension={dim}, degree={degree}")

    result_path = out_prefix.with_suffix(".json")
    result_path.write_text(json.dumps(results, indent=1) + "\n")
    print(f"wrote {result_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
