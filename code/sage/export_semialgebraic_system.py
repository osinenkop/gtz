#!/usr/bin/env python3
"""Export GTZ(6,3) semialgebraic active-set systems for Sage/Singular."""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from gtz63_semialgebraic import (  # noqa: E402
    determinant_system,
    determinant_text_system,
    indices_from_mask,
    kernel_nonsharp_system,
    known_active_indices,
    nonsharp_text_system,
    nonsharp_system,
    parse_indices,
    parse_row_pairs,
    system_summary,
    text_system_summary,
    write_singular_script,
    write_singular_text_script,
)


def resolve_active(args) -> tuple[int, ...]:
    sources = [
        args.active_indices is not None,
        args.mask is not None,
        args.known_label is not None,
    ]
    if sum(sources) != 1:
        raise SystemExit("provide exactly one of --active-indices, --mask, --known-label")
    if args.active_indices is not None:
        return parse_indices(args.active_indices)
    if args.mask is not None:
        return indices_from_mask(int(args.mask, 0))
    return known_active_indices(args.known_label)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", choices=["det", "nonsharp", "kernel"], default="det")
    parser.add_argument("--active-indices")
    parser.add_argument("--mask")
    parser.add_argument("--known-label")
    parser.add_argument("--row-pairs", default="01")
    parser.add_argument("--order", default="degrevlex")
    parser.add_argument(
        "--expand-det",
        action="store_true",
        help="expand determinant polynomials in Sage; needed for Sage dimension/Groebner",
    )
    parser.add_argument("--patch-inverses", action="store_true")
    parser.add_argument(
        "--expand-nonsharp",
        action="store_true",
        help="expand nonsharp directional polynomials in Sage; slow, needed for Sage dimension/Groebner",
    )
    parser.add_argument("--dimension", action="store_true")
    parser.add_argument("--groebner", action="store_true")
    parser.add_argument("--singular-compute", action="store_true")
    parser.add_argument(
        "--invert-d",
        action="store_true",
        help="for determinant mode, add u0*d-1=0 to remove the chart boundary d=0",
    )
    parser.add_argument("--out-prefix", default="code/sage/out/system")
    args = parser.parse_args()

    active = resolve_active(args)
    if args.mode == "det":
        if (args.dimension or args.groebner) and not args.expand_det:
            args.expand_det = True
        if args.expand_det:
            system = determinant_system(active, order=args.order, invert_d=args.invert_d)
            summary = system_summary(system)
            text_system = False
        else:
            system = determinant_text_system(active, order=args.order, invert_d=args.invert_d)
            summary = text_system_summary(system)
            text_system = True
    elif args.mode == "nonsharp":
        if (args.dimension or args.groebner) and not args.expand_nonsharp:
            raise SystemExit("--dimension/--groebner for nonsharp requires --expand-nonsharp")
        row_pairs = parse_row_pairs(args.row_pairs, len(active))
        if args.expand_nonsharp:
            system = nonsharp_system(
                active,
                row_pairs,
                order=args.order,
                patch_inverses=args.patch_inverses,
            )
            summary = system_summary(system)
            text_system = False
        else:
            system = nonsharp_text_system(
                active,
                row_pairs,
                order=args.order,
                patch_inverses=args.patch_inverses,
            )
            summary = text_system_summary(system)
            text_system = True
    else:
        system = kernel_nonsharp_system(active, order=args.order)
        summary = system_summary(system)
        text_system = False

    ideal = None
    if args.dimension or args.groebner:
        ideal = system.ring.ideal(system.equalities)
    if args.dimension:
        summary["dimension"] = int(ideal.dimension())
    if args.groebner:
        gb = ideal.groebner_basis()
        summary["groebner_basis_size"] = len(gb)
        summary["groebner_basis_degrees"] = [int(poly.total_degree()) for poly in gb]

    out_prefix = Path(args.out_prefix)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    summary_path = out_prefix.with_suffix(".json")
    singular_path = out_prefix.with_suffix(".sing")
    summary_path.write_text(json.dumps(summary, indent=1) + "\n")
    if text_system:
        write_singular_text_script(system, singular_path, order=args.order, compute=args.singular_compute)
    else:
        write_singular_script(system, singular_path, compute=args.singular_compute)

    print("=" * 78)
    print("GTZ(6,3) semialgebraic system export")
    print("=" * 78)
    print(f"mode:             {summary['mode']}")
    print(f"active size:      {summary['active_size']}")
    print(f"variables:        {summary['n_variables']}")
    print(f"equalities:       {summary['n_equalities']}")
    print(f"max eq degree:    {summary['max_equality_degree']}")
    if "dimension" in summary:
        print(f"dimension:        {summary['dimension']}")
    if "groebner_basis_size" in summary:
        print(f"groebner size:    {summary['groebner_basis_size']}")
    print(f"wrote:            {summary_path}")
    print(f"wrote:            {singular_path}")
    print("=" * 78)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
