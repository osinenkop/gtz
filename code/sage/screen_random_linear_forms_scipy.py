#!/usr/bin/env python3
"""Search random linear forms with SciPy CSR quotient actions over a prime field."""
from __future__ import annotations

import argparse
import json
import random
import sys
from pathlib import Path

import numpy as np
from scipy import sparse

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from compute_q_power_relation import make_ring, monomial_key  # noqa: E402
from screen_random_linear_forms import expression_from_coeffs  # noqa: E402


def bm_mod(sequence: list[int], prime: int) -> list[int]:
    """Berlekamp-Massey over GF(prime), returning connection [1,c1,...,cL]."""
    c = [1]
    b = [1]
    length = 0
    shift = 1
    last_discrepancy = 1

    for n, value in enumerate(sequence):
        discrepancy = value
        for i in range(1, length + 1):
            discrepancy = (discrepancy + c[i] * sequence[n - i]) % prime
        if discrepancy == 0:
            shift += 1
            continue

        previous = list(c)
        scale = discrepancy * pow(last_discrepancy, -1, prime) % prime
        if len(c) < len(b) + shift:
            c.extend([0] * (len(b) + shift - len(c)))
        for j, coeff in enumerate(b):
            c[j + shift] = (c[j + shift] - scale * coeff) % prime

        if 2 * length <= n:
            length = n + 1 - length
            b = previous
            last_discrepancy = discrepancy
            shift = 1
        else:
            shift += 1

    return c[: length + 1]


def cache_paths(prefix: Path, variable_names: tuple[str, ...]) -> tuple[Path, dict[str, Path]]:
    return (
        prefix.with_suffix(".json"),
        {name: prefix.with_name(f"{prefix.name}_{name}.npz") for name in variable_names},
    )


def load_actions(prefix: Path, variable_names: tuple[str, ...], prime: int, n: int):
    meta_path, matrix_paths = cache_paths(prefix, variable_names)
    if not meta_path.exists() or not all(path.exists() for path in matrix_paths.values()):
        return None
    meta = json.loads(meta_path.read_text())
    if (
        int(meta.get("prime", -1)) != prime
        or int(meta.get("dimension", -1)) != n
        or tuple(meta.get("variables", ())) != variable_names
    ):
        return None
    matrices = {name: sparse.load_npz(matrix_paths[name]).astype(np.int64) for name in variable_names}
    return int(meta["one_index"]), matrices


def save_actions(
    *,
    prefix: Path,
    variable_names: tuple[str, ...],
    prime: int,
    n: int,
    one_index: int,
    matrices: dict[str, sparse.csr_matrix],
) -> None:
    prefix.parent.mkdir(parents=True, exist_ok=True)
    meta_path, matrix_paths = cache_paths(prefix, variable_names)
    for name, matrix in matrices.items():
        sparse.save_npz(matrix_paths[name], matrix)
    meta_path.write_text(
        json.dumps(
            {
                "prime": prime,
                "dimension": n,
                "one_index": one_index,
                "variables": list(variable_names),
                "matrix_files": {name: str(matrix_paths[name]) for name in variable_names},
            },
            indent=1,
        )
        + "\n"
    )


def build_actions(data: dict, *, cache_prefix: Path | None):
    ring = make_ring(data)
    prime = int(data["characteristic"])
    gb = [ring(poly) for poly in data["basis"]]
    ideal = ring.ideal(gb)
    quotient_degree = int(data.get("degree") or ideal.vector_space_dimension())
    variable_names = tuple(ring.variable_names())

    if cache_prefix is not None:
        cached = load_actions(cache_prefix, variable_names, prime, quotient_degree)
        if cached is not None:
            one_index, matrices = cached
            print(f"loaded cached actions from {cache_prefix}", flush=True)
            return prime, quotient_degree, variable_names, one_index, matrices

    print("computing normal basis...", flush=True)
    normal_basis = list(ideal.normal_basis())
    basis_index = {monomial_key(monomial): i for i, monomial in enumerate(normal_basis)}
    one_index = basis_index[(0,) * len(variable_names)]
    print(f"normal basis size: {len(normal_basis)}", flush=True)

    matrices: dict[str, sparse.csr_matrix] = {}
    for name in variable_names:
        print(f"building action {name}", flush=True)
        rows: list[int] = []
        cols: list[int] = []
        values: list[int] = []
        elt = ring(name)
        for col, monomial in enumerate(normal_basis):
            image = (monomial * elt).reduce(gb)
            for exp, coeff in image.dict().items():
                value = int(coeff) % prime
                if value:
                    rows.append(basis_index[tuple(exp)])
                    cols.append(col)
                    values.append(value)
        matrix = sparse.csr_matrix(
            (np.array(values, dtype=np.int64), (rows, cols)),
            shape=(quotient_degree, quotient_degree),
            dtype=np.int64,
        )
        matrix.sum_duplicates()
        matrix.data %= prime
        matrix.eliminate_zeros()
        matrices[name] = matrix
        print(f"  nnz={matrix.nnz}", flush=True)

    if cache_prefix is not None:
        save_actions(
            prefix=cache_prefix,
            variable_names=variable_names,
            prime=prime,
            n=quotient_degree,
            one_index=one_index,
            matrices=matrices,
        )
        print(f"saved cached actions to {cache_prefix}", flush=True)

    return prime, quotient_degree, variable_names, one_index, matrices


def combine_action(
    matrices: dict[str, sparse.csr_matrix],
    variable_names: tuple[str, ...],
    coeffs: list[int],
    prime: int,
) -> sparse.csr_matrix:
    result = None
    for coeff, name in zip(coeffs, variable_names):
        coeff %= prime
        if coeff == 0:
            continue
        term = matrices[name] * coeff
        result = term if result is None else result + term
    if result is None:
        n = next(iter(matrices.values())).shape[0]
        result = sparse.csr_matrix((n, n), dtype=np.int64)
    result = result.tocsr()
    result.sum_duplicates()
    result.data %= prime
    result.eliminate_zeros()
    return result


def matvec_mod(matrix: sparse.csr_matrix, vector: np.ndarray, prime: int) -> np.ndarray:
    out = matrix.dot(vector)
    out %= prime
    return np.asarray(out, dtype=np.int64)


def dot_mod(left: np.ndarray, right: np.ndarray, prime: int) -> int:
    return int(np.dot(left, right) % prime)


def wiedemann_relation(
    *,
    matrix: sparse.csr_matrix,
    prime: int,
    one_index: int,
    max_degree: int,
    check_every: int,
    rng: random.Random,
    store_coefficients: bool,
) -> dict:
    n = matrix.shape[0]
    terms = 2 * max_degree + 1
    probe = np.array([rng.randrange(prime) for _ in range(n)], dtype=np.int64)
    if not np.any(probe):
        probe[0] = 1
    current = np.zeros(n, dtype=np.int64)
    current[one_index] = 1
    sequence: list[int] = []

    print(f"  computing {terms} Wiedemann terms", flush=True)
    for index in range(terms):
        sequence.append(dot_mod(probe, current, prime))
        if index > 0 and index % check_every == 0:
            print(f"  Wiedemann terms <= {index}", flush=True)
        if index + 1 < terms:
            current = matvec_mod(matrix, current, prime)

    connection = bm_mod(sequence, prime)
    relation = list(reversed(connection))
    degree = len(relation) - 1
    print(f"  scalar recurrence degree {degree}", flush=True)

    residual = np.zeros(n, dtype=np.int64)
    current = np.zeros(n, dtype=np.int64)
    current[one_index] = 1
    for index, coeff in enumerate(relation):
        if coeff:
            residual += coeff * current
            residual %= prime
        if index + 1 < len(relation):
            current = matvec_mod(matrix, current, prime)
    residual_zero = not np.any(residual % prime)

    nonzero_terms = sum(1 for coeff in relation if coeff)
    row = {
        "relation_found": bool(residual_zero and degree <= max_degree),
        "relation_degree": int(degree) if residual_zero and degree <= max_degree else None,
        "scalar_relation_degree": int(degree),
        "searched_degree": int(max_degree),
        "residual_zero": bool(residual_zero),
        "nonzero_terms": int(nonzero_terms),
        "density": nonzero_terms / len(relation),
    }
    if store_coefficients and residual_zero and degree <= max_degree:
        row["coefficients_low_to_high"] = [int(coeff) for coeff in relation]
    return row


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--basis-json", required=True)
    parser.add_argument("--trials", type=int, default=16)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--coeff-min", type=int, default=1)
    parser.add_argument("--coeff-max", type=int, default=29)
    parser.add_argument("--max-degree", type=int, default=0)
    parser.add_argument("--check-every", type=int, default=400)
    parser.add_argument("--target-degree", type=int, default=0)
    parser.add_argument("--cache-prefix")
    parser.add_argument("--store-coefficients", action="store_true")
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    data = json.loads(Path(args.basis_json).read_text())
    cache_prefix = Path(args.cache_prefix) if args.cache_prefix else None
    prime, quotient_degree, variable_names, one_index, matrices = build_actions(
        data,
        cache_prefix=cache_prefix,
    )
    max_degree = args.max_degree or quotient_degree
    target_degree = args.target_degree or quotient_degree

    print("=" * 78)
    print("SCIPY WIEDEMANN RANDOM LINEAR FORM SCREEN")
    print(f"basis json: {args.basis_json}")
    print(f"prime: {prime}")
    print(f"variables: {variable_names}")
    print(f"quotient degree: {quotient_degree}")
    print(f"max degree: {max_degree}")
    print(f"target degree: {target_degree}")
    print(f"action nnz: { {name: matrices[name].nnz for name in variable_names} }", flush=True)

    rng = random.Random(args.seed)
    rows = []
    for trial in range(args.trials):
        coeffs = [rng.randint(args.coeff_min, args.coeff_max) for _ in variable_names]
        expression = expression_from_coeffs(variable_names, coeffs)
        print("-" * 78)
        print(f"trial {trial}: {expression}", flush=True)
        matrix = combine_action(matrices, variable_names, coeffs, prime)
        print(f"  combined nnz={matrix.nnz}", flush=True)
        row = wiedemann_relation(
            matrix=matrix,
            prime=prime,
            one_index=one_index,
            max_degree=max_degree,
            check_every=args.check_every,
            rng=rng,
            store_coefficients=args.store_coefficients,
        )
        row["trial"] = trial
        row["expression"] = expression
        row["coefficients"] = dict(zip(variable_names, coeffs))
        rows.append(row)
        print(
            f"trial {trial}: found={row['relation_found']} "
            f"degree={row['relation_degree']} scalar_degree={row['scalar_relation_degree']}",
            flush=True,
        )
        if row["relation_found"] and int(row["relation_degree"]) >= target_degree:
            print(f"target reached at trial {trial}", flush=True)
            break

    payload = {
        "basis_json": args.basis_json,
        "characteristic": int(prime),
        "variables": list(variable_names),
        "quotient_degree": int(quotient_degree),
        "normal_basis_size": int(quotient_degree),
        "max_degree": int(max_degree),
        "check_every": int(args.check_every),
        "seed": int(args.seed),
        "trials_requested": int(args.trials),
        "method": "scipy_csr_wiedemann",
        "variable": "t",
        "action_nnz": {name: int(matrices[name].nnz) for name in variable_names},
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
