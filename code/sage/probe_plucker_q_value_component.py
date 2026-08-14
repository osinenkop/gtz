#!/usr/bin/env python3
"""Probe a Plucker ansatz component after fixing q to a rational value."""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from sage.all import GF, Matrix, QQ

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from compute_q_power_relation import monomial_key, vector_from_poly  # noqa: E402
from probe_overtie_plucker_ansatz import build_system  # noqa: E402


def relation_for_columns(field, vectors):
    if not vectors:
        return None
    nrows = len(vectors[0])
    ncols = len(vectors)
    matrix = Matrix(field, nrows, ncols, lambda i, j: vectors[j][i])
    if matrix.rank() == ncols:
        return None
    kernel = matrix.right_kernel().basis()
    if not kernel:
        return None
    relation = list(kernel[0])
    last = max((i for i, coeff in enumerate(relation) if coeff != 0), default=-1)
    if last < 0:
        return None
    lead = relation[last]
    return [coeff / lead for coeff in relation[: last + 1]]


def relation_to_string(coeffs, variable_name: str) -> str:
    parts = []
    for power in range(len(coeffs) - 1, -1, -1):
        coeff = coeffs[power]
        if coeff == 0:
            continue
        if power == 0:
            monomial = "1"
        elif power == 1:
            monomial = variable_name
        else:
            monomial = f"{variable_name}^{power}"
        if coeff == 1 and power != 0:
            term = monomial
        else:
            term = f"{coeff}*{monomial}"
        parts.append(term)
    return " + ".join(parts) if parts else "0"


def relation_for_expression(quotient_basis, normal_basis, expression, max_degree: int) -> dict:
    ring = expression.parent()
    field = ring.base_ring()
    basis_index = {monomial_key(monomial): i for i, monomial in enumerate(normal_basis)}
    vectors = []
    powers = []
    current = ring(1)
    found = None
    for power in range(max_degree + 1):
        powers.append(current)
        vectors.append(vector_from_poly(current, basis_index, field))
        relation = relation_for_columns(field, vectors)
        if relation is not None:
            found = relation
            break
        if power < max_degree:
            current = (current * expression).reduce(quotient_basis)
    if found is None:
        return {"relation_found": False, "max_degree": max_degree}
    residual = sum(coeff * powers[i] for i, coeff in enumerate(found)).reduce(quotient_basis)
    return {
        "relation_found": True,
        "relation_degree": len(found) - 1,
        "coefficients_low_to_high": [str(coeff) for coeff in found],
        "polynomial": relation_to_string(found, "u"),
        "residual_zero": bool(residual == 0),
        "residual_terms": len(residual.dict()),
    }


def parse_q_value(ring, text: str):
    base = ring.base_ring()
    if base == QQ:
        return QQ(text)
    return base(QQ(text))


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--case", default="s8_79656")
    parser.add_argument("--characteristic", type=int, default=0)
    parser.add_argument("--q-value", default="5")
    parser.add_argument(
        "--relation-expression",
        action="append",
        default=[],
        help="expression to eliminate; defaults to all variables",
    )
    parser.add_argument("--max-relation-degree", type=int, default=0)
    parser.add_argument("--store-basis", action="store_true")
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    ring, z, det_gram, polys, labels = build_system(args.case, args.characteristic, "degrevlex")
    q = ring("q")
    q_value = parse_q_value(ring, args.q_value)
    polys = list(polys) + [q - q_value]
    labels = list(labels) + [f"q:{q - q_value}"]
    relation_expressions = args.relation_expression or list(ring.variable_names())

    print("=" * 78)
    print("PLUCKER Q-VALUE COMPONENT PROBE")
    print(f"case: {args.case}")
    print(f"field: {'QQ' if args.characteristic == 0 else 'F_'+str(args.characteristic)}")
    print(f"q value: {q_value}")
    print(f"equations: {len(polys)}")
    print("computing Groebner basis...", flush=True)
    ideal = ring.ideal(polys)
    gb = list(ideal.groebner_basis())
    quotient = ring.ideal(gb)
    dimension = int(quotient.dimension())
    degree = None
    normal_basis = []
    if dimension == 0:
        degree = int(quotient.vector_space_dimension())
        print("computing normal basis...", flush=True)
        normal_basis = list(quotient.normal_basis())
    print(f"basis size: {len(gb)}")
    print(f"dimension: {dimension}")
    print(f"degree: {degree}")

    relations = []
    if dimension == 0:
        max_degree = args.max_relation_degree or degree
        for expression_text in relation_expressions:
            print(f"computing relation for {expression_text}", flush=True)
            relation = relation_for_expression(
                gb,
                normal_basis,
                ring(expression_text.replace("^", "**")),
                max_degree,
            )
            relation["expression"] = expression_text
            relations.append(relation)
            if relation.get("relation_found"):
                print(
                    f"  degree={relation['relation_degree']} "
                    f"residual_zero={relation['residual_zero']} {relation['polynomial']}",
                    flush=True,
                )
            else:
                print("  no relation found", flush=True)

    payload = {
        "case": args.case,
        "characteristic": args.characteristic,
        "field": "QQ" if args.characteristic == 0 else f"F_{args.characteristic}",
        "q_value": str(q_value),
        "variables": list(ring.variable_names()),
        "labels": labels,
        "dimension": dimension,
        "degree": degree,
        "basis_size": len(gb),
        "normal_basis_sample": [str(monomial) for monomial in normal_basis[:20]],
        "relations": relations,
    }
    if args.store_basis:
        payload["basis"] = [str(poly) for poly in gb]
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(payload, indent=1) + "\n")
    print(f"wrote {out}")
    print("=" * 78)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
