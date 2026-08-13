#!/usr/bin/env python3
"""Search random linear forms for high-degree quotient power relations."""
from __future__ import annotations

import argparse
import json
import random
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from compute_q_power_relation import make_ring, monomial_key  # noqa: E402
from screen_linear_form_relations import relation_for_expression  # noqa: E402


def expression_from_coeffs(variable_names: tuple[str, ...], coeffs: list[int]) -> str:
    parts = []
    for coeff, name in zip(coeffs, variable_names):
        if coeff == 0:
            continue
        if coeff == 1:
            parts.append(name)
        elif coeff == -1:
            parts.append(f"-{name}")
        else:
            parts.append(f"{coeff}*{name}")
    return "+".join(parts).replace("+-", "-") if parts else "0"


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--basis-json", required=True)
    parser.add_argument("--trials", type=int, default=8)
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
    gb = [ring(poly) for poly in data["basis"]]
    ideal = ring.ideal(gb)
    quotient_degree = int(data.get("degree") or ideal.vector_space_dimension())
    max_degree = args.max_degree or quotient_degree
    target_degree = args.target_degree or quotient_degree
    variable_names = tuple(ring.variable_names())

    print("=" * 78)
    print("RANDOM LINEAR FORM SCREEN")
    print(f"basis json: {args.basis_json}")
    print(f"field: {ring.base_ring()}")
    print(f"variables: {variable_names}")
    print(f"quotient degree: {quotient_degree}")
    print(f"max degree: {max_degree}")
    print(f"target degree: {target_degree}")
    print("computing normal basis...", flush=True)
    normal_basis = list(ideal.normal_basis())
    basis_index = {monomial_key(monomial): i for i, monomial in enumerate(normal_basis)}
    print(f"normal basis size: {len(normal_basis)}")

    rng = random.Random(args.seed)
    rows = []
    for trial in range(args.trials):
        coeffs = [rng.randint(args.coeff_min, args.coeff_max) for _ in variable_names]
        expression = expression_from_coeffs(variable_names, coeffs)
        print("-" * 78)
        print(f"trial {trial}: {expression}", flush=True)
        row = relation_for_expression(
            ring=ring,
            gb=gb,
            normal_basis=normal_basis,
            basis_index=basis_index,
            expression=expression,
            max_degree=max_degree,
            check_every=args.check_every,
            store_coefficients=False,
            max_polynomial_degree=0,
        )
        row["trial"] = trial
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
