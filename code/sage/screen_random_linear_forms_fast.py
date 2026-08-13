#!/usr/bin/env python3
"""Search random linear forms using quotient multiplication matrices."""
from __future__ import annotations

import argparse
import json
import math
import random
import sys
from pathlib import Path

from sage.all import Matrix, vector

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from compute_q_power_relation import make_ring, monomial_key, relation_for_columns  # noqa: E402
from screen_random_linear_forms import expression_from_coeffs  # noqa: E402


def multiplication_matrix(*, field, ring, gb, basis_index, normal_basis, elt):
    n = len(normal_basis)
    entries = {}
    for col, monomial in enumerate(normal_basis):
        image = (monomial * elt).reduce(gb)
        for exp, coeff in image.dict().items():
            row = basis_index[tuple(exp)]
            value = field(coeff)
            if value != 0:
                entries[(row, col)] = value
    return Matrix(field, n, n, entries, sparse=True)


def relation_for_matrix_powers(
    *,
    field,
    matrix,
    one_index: int,
    max_degree: int,
    check_every: int,
) -> dict:
    n = matrix.nrows()
    current = vector(field, n)
    current[one_index] = field(1)
    vectors = []
    found = None
    found_at = None
    last_checked = 0

    for power in range(max_degree + 1):
        vectors.append(current)
        should_check = (
            power == max_degree
            or power + 1 > n
            or (power > 0 and power % check_every == 0)
        )
        if should_check:
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
            current = matrix * current

    row = {
        "relation_found": found is not None,
        "relation_degree": int(found_at) if found_at is not None else None,
        "searched_degree": max_degree,
    }
    if found is None:
        row["span_rank_lower_bound"] = len(vectors)
        return row

    residual = vector(field, n)
    for index, coeff in enumerate(found):
        if coeff != 0:
            residual += coeff * vectors[index]
    nonzero = [coeff for coeff in found if coeff != 0]
    row.update(
        {
            "residual_zero": all(coeff == 0 for coeff in residual),
            "nonzero_terms": len(nonzero),
            "density": len(nonzero) / len(found),
            "degree_log10_height_proxy": math.log10(
                max(1, max(abs(int(coeff)) for coeff in found))
            ),
        }
    )
    return row


def berlekamp_massey(sequence, field):
    """Return the monic recurrence C for s_n + C_1 s_{n-1}+...+C_L s_{n-L}=0."""
    one = field(1)
    zero = field(0)
    C = [one]
    B = [one]
    L = 0
    m = 1
    b = one

    for n, value in enumerate(sequence):
        discrepancy = value
        for i in range(1, L + 1):
            discrepancy += C[i] * sequence[n - i]
        if discrepancy == 0:
            m += 1
            continue

        previous = list(C)
        scale = discrepancy / b
        if len(C) < len(B) + m:
            C.extend([zero] * (len(B) + m - len(C)))
        for j, coeff in enumerate(B):
            C[j + m] -= scale * coeff

        if 2 * L <= n:
            L = n + 1 - L
            B = previous
            b = discrepancy
            m = 1
        else:
            m += 1

    return C[: L + 1]


def relation_for_wiedemann(
    *,
    field,
    matrix,
    one_index: int,
    max_degree: int,
    check_every: int,
    rng: random.Random,
) -> dict:
    n = matrix.nrows()
    terms = 2 * max_degree + 1
    probe = vector(field, [field(rng.randrange(int(field.characteristic()))) for _ in range(n)])
    if all(entry == 0 for entry in probe):
        probe[0] = field(1)

    current = vector(field, n)
    current[one_index] = field(1)
    sequence = []
    print(f"  computing {terms} Wiedemann terms", flush=True)
    for index in range(terms):
        sequence.append(probe.dot_product(current))
        if index > 0 and index % check_every == 0:
            print(f"  Wiedemann terms <= {index}", flush=True)
        if index + 1 < terms:
            current = matrix * current

    connection = berlekamp_massey(sequence, field)
    relation = list(reversed(connection))
    degree = len(relation) - 1
    print(f"  scalar recurrence degree {degree}", flush=True)

    residual = vector(field, n)
    current = vector(field, n)
    current[one_index] = field(1)
    for index, coeff in enumerate(relation):
        if coeff != 0:
            residual += coeff * current
        if index + 1 < len(relation):
            current = matrix * current

    residual_zero = all(coeff == 0 for coeff in residual)
    nonzero = [coeff for coeff in relation if coeff != 0]
    row = {
        "relation_found": residual_zero and degree <= max_degree,
        "relation_degree": int(degree) if residual_zero and degree <= max_degree else None,
        "scalar_relation_degree": int(degree),
        "searched_degree": max_degree,
        "residual_zero": residual_zero,
        "nonzero_terms": len(nonzero),
        "density": len(nonzero) / len(relation),
        "degree_log10_height_proxy": math.log10(
            max(1, max(abs(int(coeff)) for coeff in relation))
        ),
    }
    return row


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--basis-json", required=True)
    parser.add_argument("--method", choices=("wiedemann", "rank"), default="wiedemann")
    parser.add_argument("--trials", type=int, default=16)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--coeff-min", type=int, default=1)
    parser.add_argument("--coeff-max", type=int, default=29)
    parser.add_argument("--max-degree", type=int, default=0)
    parser.add_argument("--check-every", type=int, default=400)
    parser.add_argument("--target-degree", type=int, default=0)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    data = json.loads(Path(args.basis_json).read_text())
    ring = make_ring(data)
    field = ring.base_ring()
    gb = [ring(poly) for poly in data["basis"]]
    ideal = ring.ideal(gb)
    quotient_degree = int(data.get("degree") or ideal.vector_space_dimension())
    max_degree = args.max_degree or quotient_degree
    target_degree = args.target_degree or quotient_degree
    variable_names = tuple(ring.variable_names())

    print("=" * 78)
    print("FAST RANDOM LINEAR FORM SCREEN")
    print(f"basis json: {args.basis_json}")
    print(f"field: {field}")
    print(f"variables: {variable_names}")
    print(f"quotient degree: {quotient_degree}")
    print(f"max degree: {max_degree}")
    print(f"target degree: {target_degree}")
    print("computing normal basis...", flush=True)
    normal_basis = list(ideal.normal_basis())
    basis_index = {monomial_key(monomial): i for i, monomial in enumerate(normal_basis)}
    one_index = basis_index[(0,) * len(variable_names)]
    print(f"normal basis size: {len(normal_basis)}")

    print("precomputing variable multiplication matrices...", flush=True)
    variable_matrices = {}
    for name in variable_names:
        print(f"  {name}", flush=True)
        variable_matrices[name] = multiplication_matrix(
            field=field,
            ring=ring,
            gb=gb,
            basis_index=basis_index,
            normal_basis=normal_basis,
            elt=ring(name),
        )

    rng = random.Random(args.seed)
    rows = []
    for trial in range(args.trials):
        coeffs = [rng.randint(args.coeff_min, args.coeff_max) for _ in variable_names]
        expression = expression_from_coeffs(variable_names, coeffs)
        print("-" * 78)
        print(f"trial {trial}: {expression}", flush=True)
        matrix = Matrix(field, quotient_degree, quotient_degree, sparse=True)
        for coeff, name in zip(coeffs, variable_names):
            if coeff:
                matrix += field(coeff) * variable_matrices[name]
        if args.method == "rank":
            row = relation_for_matrix_powers(
                field=field,
                matrix=matrix,
                one_index=one_index,
                max_degree=max_degree,
                check_every=args.check_every,
            )
        else:
            row = relation_for_wiedemann(
                field=field,
                matrix=matrix,
                one_index=one_index,
                max_degree=max_degree,
                check_every=args.check_every,
                rng=rng,
            )
        row["trial"] = trial
        row["expression"] = expression
        row["coefficients"] = dict(zip(variable_names, coeffs))
        rows.append(row)
        print(
            f"trial {trial}: found={row['relation_found']} "
            f"degree={row['relation_degree']}",
            flush=True,
        )
        if row["relation_found"] and int(row["relation_degree"]) >= target_degree:
            print(f"target reached at trial {trial}", flush=True)
            break

    payload = {
        "basis_json": args.basis_json,
        "characteristic": int(data["characteristic"]),
        "variables": list(variable_names),
        "quotient_degree": quotient_degree,
        "normal_basis_size": len(normal_basis),
        "max_degree": max_degree,
        "check_every": args.check_every,
        "seed": args.seed,
        "trials_requested": args.trials,
        "method": f"sparse_variable_multiplication_matrices_{args.method}",
        "rows": rows,
    }
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(payload, indent=1) + "\n")
    print("=" * 78)
    print(f"wrote {out}")
    print("=" * 78)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
