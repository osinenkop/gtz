#!/usr/bin/env python3
"""Compute a modular q-eliminant from powers of q in a quotient algebra."""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from sage.all import GF, Matrix, PolynomialRing


def make_ring(data):
    characteristic = int(data["characteristic"])
    if characteristic <= 0:
        raise SystemExit("this helper is intended for finite-field basis JSON")
    variables = tuple(data.get("variables") or data.get("lex_variable_order") or ())
    if not variables:
        # Fall back to the known variable orders used by the Plucker probes.
        n_variables = int(data["n_variables"])
        if n_variables == 8 and "ansatz" in data:
            variables = ("a", "b", "c", "d", "x", "y", "e", "q")
        elif n_variables == 10:
            variables = tuple(f"z{i}" for i in range(9)) + ("q",)
        elif n_variables == 7 and data.get("case") == "s7_78612_a_patch":
            variables = ("a", "b", "c", "x", "y", "e", "q")
        elif n_variables == 7 and data.get("case") == "s8_79656_c_patch":
            variables = ("a", "b", "c", "d", "x", "e", "q")
        else:
            raise SystemExit("basis JSON does not contain a recoverable variable order")
    return PolynomialRing(GF(characteristic), variables, order=data.get("order", "degrevlex"))


def monomial_key(poly):
    exps = poly.exponents()
    if len(exps) != 1:
        raise ValueError(f"normal basis entry is not a monomial: {poly}")
    return tuple(exps[0])


def vector_from_poly(poly, basis_index, field):
    coeffs = [field(0)] * len(basis_index)
    for exp, coeff in poly.dict().items():
        try:
            coeffs[basis_index[tuple(exp)]] = field(coeff)
        except KeyError as exc:
            raise KeyError(f"normal form contains monomial outside basis: {exp}") from exc
    return coeffs


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
    relation = [coeff / lead for coeff in relation[: last + 1]]
    return relation


def relation_to_string(coeffs, variable):
    parts = []
    for power in range(len(coeffs) - 1, -1, -1):
        coeff = coeffs[power]
        if coeff == 0:
            continue
        if power == 0:
            monomial = "1"
        elif power == 1:
            monomial = variable
        else:
            monomial = f"{variable}^{power}"
        if coeff == 1 and power != 0:
            term = monomial
        else:
            term = f"{int(coeff)}*{monomial}"
        parts.append(term)
    return " + ".join(parts) if parts else "0"


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--basis-json", required=True)
    parser.add_argument("--variable", default="q")
    parser.add_argument("--max-degree", type=int, default=0)
    parser.add_argument("--check-every", type=int, default=100)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    data = json.loads(Path(args.basis_json).read_text())
    ring = make_ring(data)
    field = ring.base_ring()
    if args.variable not in ring.variable_names():
        raise SystemExit(f"unknown variable {args.variable!r}; variables={ring.variable_names()}")
    q = ring(args.variable)
    gb = [ring(poly) for poly in data["basis"]]
    ideal = ring.ideal(gb)
    degree = int(data.get("degree") or ideal.vector_space_dimension())
    max_degree = args.max_degree or degree

    print("=" * 78)
    print("QUOTIENT POWER RELATION")
    print(f"basis json: {args.basis_json}")
    print(f"field: {field}")
    print(f"variables: {ring.variable_names()}")
    print(f"quotient degree: {degree}")
    print(f"variable: {args.variable}")
    print(f"max degree: {max_degree}")
    print("computing normal basis...", flush=True)
    normal_basis = list(ideal.normal_basis())
    basis_index = {monomial_key(monomial): i for i, monomial in enumerate(normal_basis)}
    print(f"normal basis size: {len(normal_basis)}")

    vectors = []
    powers = []
    current = ring(1)
    found = None
    found_at = None
    last_checked = 0
    for power in range(max_degree + 1):
        powers.append(current)
        vectors.append(vector_from_poly(current, basis_index, field))
        should_check = (
            power == max_degree
            or power + 1 > len(normal_basis)
            or (power > 0 and power % args.check_every == 0)
        )
        if should_check:
            print(f"checking powers <= {power}", flush=True)
            relation = relation_for_columns(field, vectors)
            if relation is not None:
                start = max(last_checked + 1, 1)
                for end_power in range(start, power + 1):
                    refined = relation_for_columns(field, vectors[: end_power + 1])
                    if refined is not None:
                        found = refined
                        found_at = end_power
                        break
                if found is None:
                    found = relation
                    found_at = len(relation) - 1
                break
            last_checked = power
        if power < max_degree:
            current = (current * q).reduce(gb)

    residual_zero = None
    residual_terms = None
    if found is not None:
        residual = sum(coeff * powers[i] for i, coeff in enumerate(found)).reduce(gb)
        residual_zero = residual == 0
        residual_terms = len(residual.dict())
        print(f"found relation degree {found_at}")
        print(f"residual zero: {residual_zero}")
        print(relation_to_string(found, args.variable))
    else:
        print("no relation found")

    payload = {
        "basis_json": args.basis_json,
        "characteristic": int(data["characteristic"]),
        "variables": list(ring.variable_names()),
        "variable": args.variable,
        "quotient_degree": degree,
        "normal_basis_size": len(normal_basis),
        "max_degree": max_degree,
        "relation_found": found is not None,
        "relation_degree": int(found_at) if found_at is not None else None,
        "coefficients_low_to_high": [int(c) for c in found] if found is not None else [],
        "polynomial": relation_to_string(found, args.variable) if found is not None else "",
        "residual_zero": residual_zero,
        "residual_terms": residual_terms,
        "normal_basis_sample": [str(m) for m in normal_basis[:20]],
    }
    if found is not None:
        nonzero = [i for i, coeff in enumerate(found) if coeff != 0]
        payload["nonzero_terms"] = len(nonzero)
        payload["density"] = len(nonzero) / len(found)
        payload["degree_log10_height_proxy"] = math.log10(max(1, max(abs(int(c)) for c in found)))

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(payload, indent=1) + "\n")
    print(f"wrote {out}")
    print("=" * 78)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
