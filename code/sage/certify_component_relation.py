#!/usr/bin/env python3
"""Test whether a candidate relation vanishes on top-dimensional components.

The two checks are:

1. ``I + <separator>`` has smaller dimension than ``I``.  This means no
   top-dimensional component is contained in the separator hyperplane.
2. ``I + <a*separator-1, b*relation-1>`` is empty.  This means relation is zero
   on the separator-open part of the determinant locus.

Together, these are a computational certificate that the relation vanishes on
the top-dimensional components, subject to the coefficient field/backend used.
"""
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

from gtz63_semialgebraic import determinant_text_system, text_system_summary  # noqa: E402
from probe_determinant_ideals import resolve_active, run_singular, write_probe_script  # noqa: E402
from probe_minor_relation import minor_relation_text, parse_minor_rows  # noqa: E402


def separator_degree(separator: str) -> int:
    if "^" in separator or "*" in separator:
        return 2
    return 1


def system_with_extra_equalities(system, mode_suffix: str, variables, equalities, degrees, metadata):
    return replace(
        system,
        mode=f"{system.mode}_{mode_suffix}",
        variables=tuple(variables),
        equalities=system.equalities + tuple(equalities),
        equality_degrees=system.equality_degrees + tuple(degrees),
        metadata={**system.metadata, **metadata},
    )


def parse_methods(text: str) -> tuple[str, ...]:
    return tuple(item.strip() for item in text.split(",") if item.strip())


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--active-indices")
    parser.add_argument("--mask")
    parser.add_argument("--known-label")
    parser.add_argument("--order", default="degrevlex")
    parser.add_argument("--invert-d", action="store_true")
    parser.add_argument("--methods", default="slimgb")
    parser.add_argument("--timeout", type=int, default=180)
    parser.add_argument("--characteristic", type=int, default=32003)
    parser.add_argument("--singular-bin", default=shutil.which("Singular") or "Singular")
    parser.add_argument("--out-prefix", default="code/sage/out/cert_component_relation")
    parser.add_argument("--separator", default="z1-z0")
    parser.add_argument("--minor-active-index", type=int, default=1)
    parser.add_argument("--minor-rows", default="0,2")
    parser.add_argument("--coefficient", type=int, default=46656)
    parser.add_argument("--variable", default="z0")
    args = parser.parse_args()

    active = resolve_active(args)
    base = determinant_text_system(active, order=args.order, invert_d=args.invert_d)
    relation, relation_degree = minor_relation_text(
        active_index=args.minor_active_index,
        minor_rows=parse_minor_rows(args.minor_rows),
        coefficient=args.coefficient,
        order=args.order,
        invert_d=args.invert_d,
        variable=args.variable,
    )
    relation_metadata = {
        "component_relation": {
            "separator": args.separator,
            "minor_relation": {
                "active_index": args.minor_active_index,
                "minor_rows": list(parse_minor_rows(args.minor_rows)),
                "coefficient": args.coefficient,
                "variable": args.variable,
                "text": relation,
            },
        }
    }
    sep_degree = separator_degree(args.separator)
    separator_system = system_with_extra_equalities(
        base,
        "separator",
        base.variables,
        [args.separator],
        [sep_degree],
        relation_metadata,
    )
    open_relation_system = system_with_extra_equalities(
        base,
        "separator_open_relation_nonzero",
        base.variables + ("a_sep", "a_rel"),
        [f"a_sep*({args.separator}) - 1", f"a_rel*({relation}) - 1"],
        [sep_degree + 1, relation_degree + 1],
        relation_metadata,
    )

    out_prefix = Path(args.out_prefix)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    results = {
        "summary": text_system_summary(base),
        "characteristic": args.characteristic,
        "timeout_seconds": args.timeout,
        "singular_bin": args.singular_bin,
        "separator": args.separator,
        "relation": relation_metadata["component_relation"]["minor_relation"],
        "tests": {},
    }
    for test_name, system in (
        ("separator_section", separator_system),
        ("separator_open_relation_nonzero", open_relation_system),
    ):
        results["tests"][test_name] = {}
        for method in parse_methods(args.methods):
            script_path = out_prefix.with_name(f"{out_prefix.name}_{test_name}_{method}.sing")
            write_probe_script(
                system,
                script_path,
                method=method,
                order=args.order,
                characteristic=args.characteristic,
                linear_forms=(),
            )
            run = run_singular(script_path, args.singular_bin, args.timeout)
            results["tests"][test_name][method] = {"script": str(script_path), **run}
            parsed = run.get("parsed", {})
            dim = parsed.get("dimension", "?") if isinstance(parsed, dict) else "?"
            degree = parsed.get("degree_value", "?") if isinstance(parsed, dict) else "?"
            status = "timeout" if run["timed_out"] else f"exit {run['exit_code']}"
            print(f"{test_name}/{method}: {status}, elapsed={run['elapsed_seconds']:.1f}s, dim={dim}, degree={degree}")

    result_path = out_prefix.with_suffix(".json")
    result_path.write_text(json.dumps(results, indent=1) + "\n")
    print(f"wrote {result_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
