#!/usr/bin/env python3
"""Probe active-plus-outside-tie determinant systems for margin minimizers.

The v30 minimisation suggests that the best point on a low-active locus often has
only a few outside triples tied at the maximum.  For a candidate active set A and
outside tie set B, this script builds the algebraic system

    det(6 P_TT - I) = 0          for T in A,
    det(6 P_TT - q I) = 0        for T in B,

in the standard chart Y=[I;Z], with d=det(Y^T Y) inverted.  The variable q is
six times the tied outside eigenvalue, so the margin is (q-1)/6.

This is only an equality probe.  PSD, ordering, and KKT/criticality inequalities
are intentionally external at this stage.
"""
from __future__ import annotations

import argparse
import json
import shutil
import sys
from pathlib import Path

from sage.all import Matrix, PolynomialRing, QQ

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from gtz63_semialgebraic import (  # noqa: E402
    TRIPLE_INDEX,
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


def active_from_size_canon(size: int, canon: int, path: str = "verify/out/v22_low_active.json"):
    data = json.loads(Path(path).read_text())
    matches = [
        row for row in data["full_pair_cover"]
        if row["size"] == size and row["canon"] == canon
    ]
    if len(matches) != 1:
        raise ValueError(f"expected one active set for size={size}, canon={canon}; got {len(matches)}")
    return tuple(TRIPLE_INDEX[tuple(triple)] for triple in matches[0]["triples"])


def resolve_active(args) -> tuple[int, ...]:
    if args.active_indices:
        if args.size is not None or args.canon is not None:
            raise SystemExit("use either --active-indices or --size/--canon, not both")
        return parse_indices(args.active_indices)
    if args.size is None or args.canon is None:
        raise SystemExit("provide --active-indices or both --size and --canon")
    return active_from_size_canon(args.size, args.canon)


def shifted_block(chart, triple_index: int, q):
    triple = TRIPLES[triple_index]
    return Matrix(
        chart.ring,
        3,
        3,
        lambda a, b: 6 * chart.n[triple[a], triple[b]]
        - (q * chart.d if a == b else chart.ring.zero()),
    )


def shifted_determinant(chart, triple_index: int, q):
    """Saturated determinant det(6*N_TT - q*d*I)/d^2."""
    raw = shifted_block(chart, triple_index, q).det()
    quotient, remainder = raw.quo_rem(chart.d ** 2)
    if remainder != chart.ring.zero():
        raise ArithmeticError("shifted determinant is not divisible by d^2")
    return quotient


def margin_tie_text_system(
    active_indices: tuple[int, ...],
    tie_indices: tuple[int, ...],
    order: str,
) -> SingularTextSystem:
    names = [f"z{i}" for i in range(9)] + ["q", "u0"]
    ring = PolynomialRing(QQ, names, order=order)
    chart = standard_chart(ring, include_h=False, inverse_count=1)
    q = dict(zip(ring.variable_names(), ring.gens()))["q"]

    equalities = []
    degrees = []
    for idx in active_indices:
        poly = active_determinant(chart, idx)
        equalities.append(str(poly))
        degrees.append(int(poly.total_degree()))
    for idx in tie_indices:
        poly = shifted_determinant(chart, idx, q)
        equalities.append(str(poly))
        degrees.append(int(poly.total_degree()))
    inv_eq = chart.inv[0] * chart.d - 1
    equalities.append(str(inv_eq))
    degrees.append(int(inv_eq.total_degree()))

    return SingularTextSystem(
        mode="margin_tie_det_inverted",
        variables=tuple(str(name) for name in ring.variable_names()),
        active_indices=active_indices,
        row_pairs=(),
        equalities=tuple(equalities),
        nonzero=(),
        equality_degrees=tuple(degrees),
        nonzero_degrees=(),
        metadata={
            "chart": "Y=[I;Z]",
            "description": "saturated active det(6P-I)=0 plus outside tie det(6P-qI)=0, with d inverted",
            "active_triples": [list(TRIPLES[i]) for i in active_indices],
            "tie_indices": list(tie_indices),
            "tie_triples": [list(TRIPLES[i]) for i in tie_indices],
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
    parser.add_argument("--order", default="degrevlex")
    parser.add_argument("--methods", default="slimgb")
    parser.add_argument("--timeout", type=int, default=300)
    parser.add_argument("--characteristic", type=int, default=32003)
    parser.add_argument("--singular-bin", default=shutil.which("Singular") or "Singular")
    parser.add_argument("--out-prefix", default="code/sage/out/margin_tie")
    parser.add_argument("--linear-sections", type=int, default=0)
    parser.add_argument(
        "--linear-section-mode",
        choices=["random", "zero-sum-z"],
        default="random",
    )
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--coefficient-bound", type=int, default=17)
    args = parser.parse_args()

    active = resolve_active(args)
    ties = parse_indices(args.tie_indices)
    if set(active) & set(ties):
        raise SystemExit("tie indices must be outside the active set")

    system = margin_tie_text_system(active, ties, args.order)
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
