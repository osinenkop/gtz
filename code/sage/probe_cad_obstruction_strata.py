#!/usr/bin/env python3
"""Probe CAD-sized obstruction strata inside the determinant locus.

The first target is the residual component seen in exact lex sections of the
known-base active set.  On that component, the candidate relation

    principal_minor((0,1,3), rows 0,2) + 46656*z0^2 = 0

turns an active-PSD inequality into the equality ``z0=0``.  This script measures
the algebraic cost of imposing that equality before attempting a real CAD/QE
certificate on the resulting smaller stratum.
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


def parse_methods(text: str) -> tuple[str, ...]:
    return tuple(item.strip() for item in text.split(",") if item.strip())


def parse_equalities(text: str) -> tuple[str, ...]:
    return tuple(item.strip() for item in text.split(",") if item.strip())


def parse_compact_rows(text: str) -> tuple[int, int]:
    if "," in text:
        return parse_minor_rows(text)
    rows = tuple(int(char) for char in text.strip())
    if len(rows) != 2 or len(set(rows)) != 2 or any(row < 0 or row >= 3 for row in rows):
        raise SystemExit("minor rows must be two distinct local row indices, e.g. 02 or 0,2")
    return rows


def parse_relation_specs(text: str) -> tuple[tuple[int, tuple[int, int], int, str], ...]:
    """Parse active:rows:coefficient:variable entries separated by semicolons."""
    specs = []
    for item in text.split(";"):
        item = item.strip()
        if not item:
            continue
        parts = item.split(":")
        if len(parts) != 4:
            raise SystemExit(
                "--extra-minor-relations entries must have active:rows:coefficient:variable, "
                "for example 1:12:46656:z1"
            )
        active_index, rows_text, coefficient, variable = parts
        specs.append((int(active_index), parse_compact_rows(rows_text), int(coefficient), variable))
    return tuple(specs)


def parse_cases(text: str) -> tuple[str, ...]:
    valid = {"determinant", "relation", "relation_forced"}
    cases = tuple(item.strip() for item in text.split(",") if item.strip())
    bad = [case for case in cases if case not in valid]
    if bad:
        raise SystemExit(f"unknown --cases entries {bad}; valid entries are {sorted(valid)}")
    return cases


def equality_degree(text: str) -> int:
    """Small metadata-only degree heuristic for generated Singular summaries."""
    if "^" not in text and "*" not in text:
        return 1
    degree = 1
    for chunk in text.replace("-", "+").split("+"):
        factors = [factor.strip() for factor in chunk.split("*") if factor.strip()]
        term_degree = 0
        for factor in factors:
            if "^" in factor:
                try:
                    term_degree += int(factor.rsplit("^", 1)[1])
                except ValueError:
                    term_degree += 1
            elif factor:
                term_degree += 1
        degree = max(degree, term_degree)
    return degree


def with_extra_equalities(system, mode_suffix: str, equalities: tuple[str, ...], metadata: dict):
    return replace(
        system,
        mode=f"{system.mode}_{mode_suffix}",
        equalities=system.equalities + equalities,
        equality_degrees=system.equality_degrees + tuple(equality_degree(eq) for eq in equalities),
        metadata={**system.metadata, **metadata},
    )


def run_case(case_name: str, system, args, out_prefix: Path) -> dict:
    results = {}
    for method in parse_methods(args.methods):
        script_path = out_prefix.with_name(f"{out_prefix.name}_{case_name}_{method}.sing")
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
        basis = parsed.get("basis_size", "?") if isinstance(parsed, dict) else "?"
        status = "timeout" if run["timed_out"] else f"exit {run['exit_code']}"
        print(
            f"{case_name}/{method}: {status}, elapsed={run['elapsed_seconds']:.1f}s, "
            f"dim={dim}, degree={degree}, basis={basis}"
        )
    return results


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
    parser.add_argument("--out-prefix", default="code/sage/out/cad_obstruction_strata")
    parser.add_argument(
        "--cases",
        default="determinant,relation,relation_forced",
        help="comma-separated subset of determinant,relation,relation_forced",
    )
    parser.add_argument("--minor-active-index", type=int, default=1)
    parser.add_argument("--minor-rows", default="0,2")
    parser.add_argument("--coefficient", type=int, default=46656)
    parser.add_argument("--variable", default="z0")
    parser.add_argument(
        "--extra-minor-relations",
        default="",
        help="semicolon-separated active:rows:coefficient:variable entries, e.g. 1:12:46656:z1",
    )
    parser.add_argument(
        "--forced-equalities",
        default="z0",
        help="comma-separated equalities forced by PSD on the relation component",
    )
    args = parser.parse_args()

    active = resolve_active(args)
    base = determinant_text_system(active, order=args.order, invert_d=args.invert_d)
    minor_rows = parse_minor_rows(args.minor_rows)
    relation_specs = (
        (args.minor_active_index, minor_rows, args.coefficient, args.variable),
        *parse_relation_specs(args.extra_minor_relations),
    )
    relations = []
    relation_degrees = []
    relation_records = []
    for active_index, rows, coefficient, variable in relation_specs:
        relation, relation_degree = minor_relation_text(
            active_index=active_index,
            minor_rows=rows,
            coefficient=coefficient,
            order=args.order,
            invert_d=args.invert_d,
            variable=variable,
        )
        relations.append(relation)
        relation_degrees.append(relation_degree)
        relation_records.append({
            "active_index": active_index,
            "minor_rows": list(rows),
            "coefficient": coefficient,
            "variable": variable,
            "text": relation,
            "degree": relation_degree,
        })
    relation_metadata = {
        "cad_obstruction": {
            "minor_relations": relation_records,
            "forced_equalities": list(parse_equalities(args.forced_equalities)),
            "interpretation": (
                "If these relations hold on the relevant residual component and "
                "the active blocks are PSD, then each represented principal minor "
                "is nonnegative and the corresponding square variable is forced "
                "to zero over the reals."
            ),
        }
    }
    relation_system = replace(
        base,
        mode=f"{base.mode}_minor_relation",
        equalities=base.equalities + tuple(relations),
        equality_degrees=base.equality_degrees + tuple(relation_degrees),
        metadata={**base.metadata, **relation_metadata},
    )
    forced = parse_equalities(args.forced_equalities)
    forced_system = with_extra_equalities(
        relation_system,
        "forced_psd_equalities",
        forced,
        relation_metadata,
    )

    out_prefix = Path(args.out_prefix)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    results = {
        "summary": text_system_summary(base),
        "characteristic": args.characteristic,
        "timeout_seconds": args.timeout,
        "singular_bin": args.singular_bin,
        "methods": list(parse_methods(args.methods)),
        "requested_cases": list(parse_cases(args.cases)),
        "relation": relation_metadata["cad_obstruction"]["minor_relations"][0],
        "relations": relation_metadata["cad_obstruction"]["minor_relations"],
        "forced_equalities": list(forced),
        "cases": {},
    }
    systems = {
        "determinant": base,
        "relation": relation_system,
        "relation_forced": forced_system,
    }
    for case_name in parse_cases(args.cases):
        system = systems[case_name]
        results["cases"][case_name] = {
            "mode": system.mode,
            "n_equalities": len(system.equalities),
            "max_equality_degree": max(system.equality_degrees) if system.equality_degrees else 0,
            "runs": run_case(case_name, system, args, out_prefix),
        }

    result_path = out_prefix.with_suffix(".json")
    result_path.write_text(json.dumps(results, indent=1) + "\n")
    print(f"wrote {result_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
