#!/usr/bin/env python3
"""Probe quotient components cut out by a factor of a linear-form relation."""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from sage.all import GF, Matrix, PolynomialRing

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from compute_q_power_relation import make_ring, monomial_key, vector_from_poly  # noqa: E402


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
            term = f"{int(coeff)}*{monomial}"
        parts.append(term)
    return " + ".join(parts) if parts else "0"


def relation_for_expression(
    *,
    quotient_basis,
    normal_basis,
    expression,
    max_degree: int,
    display_variable: str,
) -> dict:
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
        "coefficients_low_to_high": [int(coeff) for coeff in found],
        "polynomial": relation_to_string(found, display_variable),
        "residual_zero": bool(residual == 0),
        "residual_terms": len(residual.dict()),
    }


def substituted_factor(ring, factor_text: str, variable: str, expression):
    univariate = PolynomialRing(ring.base_ring(), variable)
    t = univariate.gen()
    factor = univariate(factor_text.replace("^", "**"))
    value = ring.zero()
    for exponent in range(factor.degree() + 1):
        coeff = ring.base_ring()(factor.monomial_coefficient(t**exponent))
        if coeff:
            value += coeff * expression**exponent
    return factor, value


def parse_expression(ring, expression_text: str):
    return ring(expression_text.replace("^", "**"))


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--basis-json", required=True)
    parser.add_argument("--factor-json", required=True)
    parser.add_argument("--factor-index", type=int, required=True)
    parser.add_argument(
        "--relation-expression",
        action="append",
        default=[],
        help="second expression to eliminate on the component; may be repeated",
    )
    parser.add_argument("--max-relation-degree", type=int, default=0)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    basis_data = json.loads(Path(args.basis_json).read_text())
    factor_data = json.loads(Path(args.factor_json).read_text())
    ring = make_ring(basis_data)
    prime = int(basis_data["characteristic"])
    if int(factor_data["characteristic"]) != prime:
        raise SystemExit("basis and factor JSON characteristics differ")
    factors = factor_data["factors"]
    if args.factor_index < 0 or args.factor_index >= len(factors):
        raise SystemExit(f"factor index {args.factor_index} out of range for {len(factors)}")

    linear_form_text = factor_data.get("expression") or factor_data.get("variable", "q")
    linear_expression = parse_expression(ring, linear_form_text)
    factor_record = factors[args.factor_index]
    factor, component_equation = substituted_factor(
        ring,
        factor_record["polynomial"],
        factor_data.get("variable", "t"),
        linear_expression,
    )

    gb = [ring(poly) for poly in basis_data["basis"]]
    augmented_generators = gb + [component_equation]
    ideal = ring.ideal(augmented_generators)

    print("=" * 78)
    print("LINEAR FACTOR COMPONENT PROBE")
    print(f"basis json: {args.basis_json}")
    print(f"factor json: {args.factor_json}")
    print(f"factor index: {args.factor_index}")
    print(f"factor degree: {factor_record['degree']}")
    print(f"factor multiplicity: {factor_record['multiplicity']}")
    print(f"linear form: {linear_form_text}")
    print(f"factor: {factor}")
    print("computing augmented Groebner basis...", flush=True)
    augmented_basis = list(ideal.groebner_basis())
    quotient = ring.ideal(augmented_basis)
    dimension = int(quotient.dimension())
    degree = None
    normal_basis = []
    if dimension == 0:
        degree = int(quotient.vector_space_dimension())
        print("computing normal basis...", flush=True)
        normal_basis = list(quotient.normal_basis())
    print(f"augmented basis size: {len(augmented_basis)}")
    print(f"dimension: {dimension}")
    print(f"degree: {degree}")

    relations = []
    if dimension == 0:
        max_degree = args.max_relation_degree or degree
        for expression_text in args.relation_expression:
            expression = parse_expression(ring, expression_text)
            print(f"computing relation for {expression_text} up to degree {max_degree}", flush=True)
            relation = relation_for_expression(
                quotient_basis=augmented_basis,
                normal_basis=normal_basis,
                expression=expression,
                max_degree=max_degree,
                display_variable="u",
            )
            relation["expression"] = expression_text
            relations.append(relation)
            if relation.get("relation_found"):
                print(
                    f"  relation degree {relation['relation_degree']} "
                    f"residual_zero={relation['residual_zero']}",
                    flush=True,
                )
            else:
                print("  no relation found", flush=True)

    payload = {
        "basis_json": args.basis_json,
        "factor_json": args.factor_json,
        "factor_index": args.factor_index,
        "characteristic": prime,
        "variables": list(ring.variable_names()),
        "linear_form": linear_form_text,
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
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(payload, indent=1) + "\n")
    print(f"wrote {out}")
    print("=" * 78)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
