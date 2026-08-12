#!/usr/bin/env python3
"""Multiplier-free rank probes for low-active margin KKT strata.

The determinant-gradient KKT equations are linear in the multiplier variables.
For active equations f_T(z)=0 and outside tie equations h_S(z,q)=0, write

    M(z,q) c = (0,...,0,-1)^T

where M has one column for each active/tie constraint, the first nine rows are
z-gradients, and the final row is the q-derivative.  This script eliminates the
multipliers by adding rank conditions on the augmented matrix [M | rhs].

For m=|A|+|B|<10 and on a branch where M has rank m, existence of multipliers is
equivalent to all (m+1)x(m+1) minors of [M | rhs] vanishing.  These minors are
high degree, so this is an experiment, not the final certificate format.
"""
from __future__ import annotations

import argparse
import itertools
import json
import shutil
import sys
from pathlib import Path

from sage.all import Matrix, PolynomialRing, QQ

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from gtz63_semialgebraic import (  # noqa: E402
    TRIPLES,
    SingularTextSystem,
    active_determinant,
    parse_indices,
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


def parse_rows(text: str, expected: int | None = None) -> tuple[int, ...]:
    rows = tuple(int(x.strip()) for x in text.split(",") if x.strip())
    bad = [r for r in rows if r < 0 or r >= 10]
    if bad:
        raise ValueError(f"row indices must be in 0..9, got {bad}")
    if expected is not None and len(rows) != expected:
        raise ValueError(f"expected {expected} rows, got {len(rows)}")
    if len(set(rows)) != len(rows):
        raise ValueError("duplicate row indices")
    return rows


def margin_kkt_rank_text_system(
    active_indices: tuple[int, ...],
    tie_indices: tuple[int, ...],
    order: str,
    rank_rows: tuple[int, ...] | None,
    minor_limit: int,
) -> SingularTextSystem:
    m = len(active_indices) + len(tie_indices)
    if m >= 10:
        raise ValueError("rank-minor elimination is useful only for |A|+|B| < 10")

    names = [f"z{i}" for i in range(9)] + ["q", "u0"]
    if rank_rows is not None:
        names.append("r0")
    ring = PolynomialRing(QQ, names, order=order)
    gens = dict(zip(ring.variable_names(), ring.gens()))
    chart = standard_chart(ring, include_h=False, inverse_count=1)
    z = chart.z
    q = gens["q"]

    active_polys = [active_determinant(chart, idx) for idx in active_indices]
    tie_polys = [shifted_determinant(chart, idx, q) for idx in tie_indices]
    all_polys = active_polys + tie_polys

    equalities = []
    degrees = []
    for poly in all_polys:
        equalities.append(poly)
        degrees.append(int(poly.total_degree()))

    inv_eq = chart.inv[0] * chart.d - 1
    equalities.append(inv_eq)
    degrees.append(int(inv_eq.total_degree()))

    rows = []
    for zj in z:
        rows.append([poly.derivative(zj) for poly in all_polys] + [ring.zero()])
    rows.append(
        [ring.zero() for _ in active_polys]
        + [poly.derivative(q) for poly in tie_polys]
        + [ring.one()]
    )
    aug = Matrix(ring, rows)
    mat = aug[:, :m]

    minor_rows = list(itertools.combinations(range(10), m + 1))
    if minor_limit > 0:
        minor_rows = minor_rows[:minor_limit]
    for row_set in minor_rows:
        poly = aug[list(row_set), :].det()
        equalities.append(poly)
        degrees.append(int(poly.total_degree()))

    if rank_rows is not None:
        if len(rank_rows) != m:
            raise ValueError(f"rank branch needs {m} rows, got {len(rank_rows)}")
        rank_minor = mat[list(rank_rows), :].det()
        branch_eq = gens["r0"] * rank_minor - 1
        equalities.append(branch_eq)
        degrees.append(int(branch_eq.total_degree()))

    return SingularTextSystem(
        mode="margin_kkt_rank_inverted",
        variables=tuple(str(name) for name in ring.variable_names()),
        active_indices=active_indices,
        row_pairs=(),
        equalities=tuple(str(poly) for poly in equalities),
        nonzero=(),
        equality_degrees=tuple(degrees),
        nonzero_degrees=(),
        metadata={
            "chart": "Y=[I;Z]",
            "description": "multiplier-free rank conditions for determinant-gradient KKT",
            "active_triples": [list(TRIPLES[i]) for i in active_indices],
            "tie_indices": list(tie_indices),
            "tie_triples": [list(TRIPLES[i]) for i in tie_indices],
            "rank_rows": list(rank_rows) if rank_rows is not None else None,
            "n_augmented_minors": len(minor_rows),
            "augmented_minor_size": m + 1,
            "q_meaning": "q = 6 * tied outside eigenvalue; margin = (q-1)/6",
            "ring_order": order,
            "invert_d": True,
        },
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--active-indices")
    parser.add_argument("--size", type=int)
    parser.add_argument("--canon", type=int)
    parser.add_argument("--tie-indices", required=True)
    parser.add_argument("--rank-rows", default="")
    parser.add_argument("--minor-limit", type=int, default=0)
    parser.add_argument("--order", default="degrevlex")
    parser.add_argument("--methods", default="slimgb")
    parser.add_argument("--timeout", type=int, default=300)
    parser.add_argument("--characteristic", type=int, default=32003)
    parser.add_argument("--singular-bin", default=shutil.which("Singular") or "Singular")
    parser.add_argument("--out-prefix", default="code/sage/out/margin_kkt_rank")
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
    ties = parse_indices(args.tie_indices)
    if set(active) & set(ties):
        raise SystemExit("tie indices must be outside the active set")
    m = len(active) + len(ties)
    rank_rows = parse_rows(args.rank_rows, expected=m) if args.rank_rows else None

    system = margin_kkt_rank_text_system(
        active,
        ties,
        args.order,
        rank_rows=rank_rows,
        minor_limit=args.minor_limit,
    )
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
