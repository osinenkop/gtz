#!/usr/bin/env python3
"""Probe determinant ideals after adding a candidate minor relation."""
from __future__ import annotations

import argparse
import json
import shutil
import sys
from dataclasses import replace
from pathlib import Path

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from sage.all import Matrix  # noqa: E402

from gtz63_semialgebraic import (  # noqa: E402
    determinant_text_system,
    make_ring,
    standard_chart,
    active_block,
)
from probe_determinant_ideals import (  # noqa: E402
    random_linear_forms,
    resolve_active,
    run_singular,
    write_probe_script,
    zero_sum_z_linear_forms,
)


def parse_minor_rows(text: str) -> tuple[int, int]:
    rows = tuple(int(x.strip()) for x in text.split(",") if x.strip())
    if len(rows) != 2 or len(set(rows)) != 2 or any(row < 0 or row >= 3 for row in rows):
        raise SystemExit("--minor-rows must be two distinct local row indices, e.g. 0,2")
    return rows


def minor_relation_text(
    active_index: int,
    minor_rows: tuple[int, int],
    coefficient: int,
    order: str,
    invert_d: bool,
    variable: str,
) -> tuple[str, int]:
    ring = make_ring(include_h=False, inverse_count=1 if invert_d else 0, order=order)
    chart = standard_chart(ring, include_h=False, inverse_count=1 if invert_d else 0)
    block = active_block(chart, active_index)
    minor = Matrix(
        ring,
        2,
        2,
        lambda a, b: block[minor_rows[a], minor_rows[b]],
    ).det()
    gens = dict(zip(ring.variable_names(), ring.gens()))
    relation = minor + coefficient * gens[variable] ** 2
    return str(relation), int(relation.total_degree())


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--active-indices")
    parser.add_argument("--mask")
    parser.add_argument("--known-label")
    parser.add_argument("--order", default="degrevlex")
    parser.add_argument("--invert-d", action="store_true")
    parser.add_argument("--method", default="slimgb")
    parser.add_argument("--timeout", type=int, default=180)
    parser.add_argument("--characteristic", type=int, default=32003)
    parser.add_argument("--singular-bin", default=shutil.which("Singular") or "Singular")
    parser.add_argument("--out-prefix", default="code/sage/out/probe_minor_relation")
    parser.add_argument("--linear-sections", type=int, default=0)
    parser.add_argument(
        "--linear-section-mode",
        choices=["random", "zero-sum-z"],
        default="zero-sum-z",
    )
    parser.add_argument("--seed", type=int, default=401)
    parser.add_argument("--coefficient-bound", type=int, default=17)
    parser.add_argument("--minor-active-index", type=int, default=1)
    parser.add_argument("--minor-rows", default="0,2")
    parser.add_argument("--coefficient", type=int, default=46656)
    parser.add_argument("--variable", default="z0")
    args = parser.parse_args()

    active = resolve_active(args)
    system = determinant_text_system(active, order=args.order, invert_d=args.invert_d)
    relation, relation_degree = minor_relation_text(
        active_index=args.minor_active_index,
        minor_rows=parse_minor_rows(args.minor_rows),
        coefficient=args.coefficient,
        order=args.order,
        invert_d=args.invert_d,
        variable=args.variable,
    )
    system = replace(
        system,
        mode=f"{system.mode}_minor_relation",
        equalities=system.equalities + (relation,),
        equality_degrees=system.equality_degrees + (relation_degree,),
        metadata={
            **system.metadata,
            "minor_relation": {
                "active_index": args.minor_active_index,
                "minor_rows": list(parse_minor_rows(args.minor_rows)),
                "coefficient": args.coefficient,
                "variable": args.variable,
                "text": relation,
            },
        },
    )
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
    script_path = out_prefix.with_name(f"{out_prefix.name}_{args.method}.sing")
    write_probe_script(
        system,
        script_path,
        method=args.method,
        order=args.order,
        characteristic=args.characteristic,
        linear_forms=linear_forms,
    )
    run = run_singular(script_path, args.singular_bin, args.timeout)
    results = {
        "active_indices": list(active),
        "characteristic": args.characteristic,
        "linear_sections": args.linear_sections,
        "linear_section_mode": args.linear_section_mode,
        "linear_section_seed": args.seed,
        "linear_forms": list(linear_forms),
        "minor_relation": system.metadata["minor_relation"],
        "method": args.method,
        "script": str(script_path),
        "run": run,
    }
    result_path = out_prefix.with_suffix(".json")
    result_path.write_text(json.dumps(results, indent=1) + "\n")
    status = "timeout" if run["timed_out"] else f"exit {run['exit_code']}"
    parsed = run.get("parsed", {})
    dim = parsed.get("dimension", "?") if isinstance(parsed, dict) else "?"
    degree = parsed.get("degree_value", "?") if isinstance(parsed, dict) else "?"
    print(f"{args.method}: {status}, elapsed={run['elapsed_seconds']:.1f}s, dimension={dim}, degree={degree}")
    print(f"wrote {result_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
