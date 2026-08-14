#!/usr/bin/env python3
"""Batch-probe components cut out by repeated factors of a linear-form relation."""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from compute_q_power_relation import make_ring  # noqa: E402
from probe_linear_factor_component import (  # noqa: E402
    parse_expression,
    relation_for_expression,
    substituted_factor,
)


def write_payload(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(payload, indent=1) + "\n")
    tmp.replace(path)


def load_existing(path: Path) -> dict | None:
    if not path.exists():
        return None
    return json.loads(path.read_text())


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--basis-json", required=True)
    parser.add_argument("--factor-json", required=True)
    parser.add_argument("--min-multiplicity", type=int, default=2)
    parser.add_argument("--degree", type=int, default=0, help="optional selected factor degree")
    parser.add_argument(
        "--relation-expression",
        action="append",
        default=[],
        help="expression to eliminate on each component; may be repeated",
    )
    parser.add_argument("--max-relation-degree", type=int, default=0)
    parser.add_argument("--store-basis", action="store_true")
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    basis_data = json.loads(Path(args.basis_json).read_text())
    factor_data = json.loads(Path(args.factor_json).read_text())
    ring = make_ring(basis_data)
    prime = int(basis_data["characteristic"])
    if int(factor_data["characteristic"]) != prime:
        raise SystemExit("basis and factor JSON characteristics differ")

    relation_expressions = args.relation_expression or list(ring.variable_names())
    linear_form_text = factor_data.get("expression") or factor_data.get("variable", "q")
    linear_expression = parse_expression(ring, linear_form_text)
    gb = [ring(poly) for poly in basis_data["basis"]]
    out = Path(args.out)

    selected_indices = []
    for index, item in enumerate(factor_data["factors"]):
        if int(item["multiplicity"]) < args.min_multiplicity:
            continue
        if args.degree and int(item["degree"]) != args.degree:
            continue
        selected_indices.append(index)

    existing = load_existing(out)
    components = [] if existing is None else existing.get("components", [])
    done = {int(row["factor_index"]) for row in components}
    payload = {
        "basis_json": args.basis_json,
        "factor_json": args.factor_json,
        "characteristic": prime,
        "variables": list(ring.variable_names()),
        "linear_form": linear_form_text,
        "min_multiplicity": args.min_multiplicity,
        "degree_filter": args.degree or None,
        "relation_expressions": relation_expressions,
        "selected_factor_indices": selected_indices,
        "components": components,
    }

    print("=" * 78)
    print("REPEATED FACTOR COMPONENT BATCH")
    print(f"basis json: {args.basis_json}")
    print(f"factor json: {args.factor_json}")
    print(f"selected factors: {selected_indices}")
    print(f"already done: {sorted(done)}")
    print(f"relations: {relation_expressions}", flush=True)

    for factor_index in selected_indices:
        if factor_index in done:
            continue
        factor_record = factor_data["factors"][factor_index]
        factor, component_equation = substituted_factor(
            ring,
            factor_record["polynomial"],
            factor_data.get("variable", "t"),
            linear_expression,
        )
        print("-" * 78)
        print(
            f"factor {factor_index}: degree={factor_record['degree']} "
            f"multiplicity={factor_record['multiplicity']} {factor}",
            flush=True,
        )
        ideal = ring.ideal(gb + [component_equation])
        augmented_basis = list(ideal.groebner_basis())
        quotient = ring.ideal(augmented_basis)
        dimension = int(quotient.dimension())
        degree = None
        normal_basis = []
        if dimension == 0:
            degree = int(quotient.vector_space_dimension())
            normal_basis = list(quotient.normal_basis())
        print(f"  basis_size={len(augmented_basis)} dimension={dimension} degree={degree}", flush=True)

        relations = []
        if dimension == 0:
            max_degree = args.max_relation_degree or degree
            for expression_text in relation_expressions:
                relation = relation_for_expression(
                    quotient_basis=augmented_basis,
                    normal_basis=normal_basis,
                    expression=parse_expression(ring, expression_text),
                    max_degree=max_degree,
                    display_variable="u",
                )
                relation["expression"] = expression_text
                relations.append(relation)
                if relation.get("relation_found"):
                    print(
                        f"    {expression_text}: degree={relation['relation_degree']} "
                        f"{relation['polynomial']}",
                        flush=True,
                    )
                else:
                    print(f"    {expression_text}: no relation", flush=True)

        component = {
            "factor_index": factor_index,
            "factor": {
                **factor_record,
                "substituted_equation_total_degree": int(component_equation.total_degree()),
                "substituted_equation_terms": len(component_equation.dict()),
            },
            "augmented_basis_size": len(augmented_basis),
            "dimension": dimension,
            "degree": degree,
            "normal_basis_sample": [str(monomial) for monomial in normal_basis[:20]],
            "relations": relations,
        }
        if args.store_basis:
            component["augmented_basis"] = [str(poly) for poly in augmented_basis]
        components.append(component)
        payload["components"] = components
        write_payload(out, payload)
        print(f"  checkpointed {out}", flush=True)

    write_payload(out, payload)
    print("=" * 78)
    print(f"wrote {out}")
    print("=" * 78)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
