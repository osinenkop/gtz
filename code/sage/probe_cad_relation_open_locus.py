#!/usr/bin/env python3
"""Check whether a CAD obstruction relation is nonzero only off the top stratum.

For a premise ideal B and target relation R, this computes dimensions of

    B
    B + <a*R - 1>

If the second dimension is strictly smaller than the first, then R vanishes on
all top-dimensional components of B.  This is weaker than ideal membership but
is exactly the component-level certificate needed for the PSD-forcing cascade.
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
from probe_cad_obstruction_strata import equality_degree, parse_relation_specs  # noqa: E402
from probe_cad_relation_normal_forms import parse_one_relation, relation_text_from_spec  # noqa: E402
from probe_determinant_ideals import resolve_active, run_singular, write_probe_script  # noqa: E402


def parse_equalities(text: str) -> tuple[str, ...]:
    return tuple(item.strip() for item in text.split(",") if item.strip())


def parse_methods(text: str) -> tuple[str, ...]:
    return tuple(item.strip() for item in text.split(",") if item.strip())


def with_premises(system, relation_specs, equalities, order: str, invert_d: bool):
    relation_records = []
    relations = []
    relation_degrees = []
    for spec in relation_specs:
        relation, degree, record = relation_text_from_spec(spec, order, invert_d)
        relations.append(relation)
        relation_degrees.append(degree)
        relation_records.append(record)
    return replace(
        system,
        mode=f"{system.mode}_open_premise",
        equalities=system.equalities + tuple(relations) + tuple(equalities),
        equality_degrees=system.equality_degrees
        + tuple(relation_degrees)
        + tuple(equality_degree(eq) for eq in equalities),
        metadata={
            **system.metadata,
            "premise_relations": relation_records,
            "premise_equalities": list(equalities),
        },
    )


def with_target_inverse(system, target: str, target_degree: int):
    return replace(
        system,
        mode=f"{system.mode}_target_nonzero",
        variables=system.variables + ("a_rel",),
        equalities=system.equalities + (f"a_rel*({target}) - 1",),
        equality_degrees=system.equality_degrees + (target_degree + 1,),
        metadata={**system.metadata, "target_nonzero_encoded_by": "a_rel*target-1"},
    )


def run_system(name: str, system, args, out_prefix: Path) -> dict:
    results = {}
    for method in parse_methods(args.methods):
        script_path = out_prefix.with_name(f"{out_prefix.name}_{name}_{method}.sing")
        write_probe_script(
            system,
            script_path,
            method=method,
            order=args.order,
            characteristic=args.characteristic,
            linear_forms=(),
        )
        run = run_singular(script_path, args.singular_bin, args.timeout)
        results[method] = {"script": str(script_path), **run}
        parsed = run.get("parsed", {})
        dim = parsed.get("dimension", "?") if isinstance(parsed, dict) else "?"
        degree = parsed.get("degree_value", "?") if isinstance(parsed, dict) else "?"
        status = "timeout" if run["timed_out"] else f"exit {run['exit_code']}"
        print(f"{name}/{method}: {status}, elapsed={run['elapsed_seconds']:.1f}s, dim={dim}, degree={degree}")
    return results


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--active-indices")
    parser.add_argument("--mask")
    parser.add_argument("--known-label")
    parser.add_argument("--order", default="degrevlex")
    parser.add_argument("--invert-d", action="store_true")
    parser.add_argument("--premise-relations", default="")
    parser.add_argument("--premise-equalities", default="")
    parser.add_argument("--target-relation", required=True)
    parser.add_argument("--methods", default="slimgb")
    parser.add_argument("--timeout", type=int, default=240)
    parser.add_argument("--characteristic", type=int, default=32003)
    parser.add_argument("--singular-bin", default=shutil.which("Singular") or "Singular")
    parser.add_argument("--out-prefix", default="code/sage/out/cad_relation_open")
    parser.add_argument("--skip-base", action="store_true")
    args = parser.parse_args()

    active = resolve_active(args)
    base = determinant_text_system(active, order=args.order, invert_d=args.invert_d)
    premise = with_premises(
        base,
        parse_relation_specs(args.premise_relations),
        parse_equalities(args.premise_equalities),
        args.order,
        args.invert_d,
    )
    target, target_degree, target_record = relation_text_from_spec(
        parse_one_relation(args.target_relation),
        args.order,
        args.invert_d,
    )
    target_open = with_target_inverse(premise, target, target_degree)

    out_prefix = Path(args.out_prefix)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    results = {
        "summary": text_system_summary(premise),
        "characteristic": args.characteristic,
        "timeout_seconds": args.timeout,
        "premise_relations": premise.metadata.get("premise_relations", []),
        "premise_equalities": premise.metadata.get("premise_equalities", []),
        "target_relation": target_record,
        "tests": {},
    }
    if not args.skip_base:
        results["tests"]["base"] = run_system("base", premise, args, out_prefix)
    results["tests"]["target_nonzero"] = run_system("target_nonzero", target_open, args, out_prefix)

    result_path = out_prefix.with_suffix(".json")
    result_path.write_text(json.dumps(results, indent=1) + "\n")
    print(f"wrote {result_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
