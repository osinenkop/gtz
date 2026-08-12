#!/usr/bin/env python3
"""Factor a saved separator candidate from probe_local_separator_membership.py."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from sage.all import GF, PolynomialRing, QQ


def base_field(characteristic: int):
    return QQ if characteristic == 0 else GF(characteristic)


def lift_coeff(coeff, characteristic: int):
    if characteristic == 0:
        return QQ(coeff)
    value = int(coeff)
    if value > characteristic // 2:
        value -= characteristic
    return QQ(value)


def monomial_string(names, exponents):
    pieces = []
    for name, exponent in zip(names, exponents):
        if exponent == 1:
            pieces.append(name)
        elif exponent:
            pieces.append(f"{name}^{exponent}")
    return "*".join(pieces) if pieces else "1"


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--candidate-index", type=int, default=0)
    parser.add_argument("--all", action="store_true", help="factor every passing candidate")
    parser.add_argument("--characteristic", type=int, default=-1)
    parser.add_argument("--factor-index", type=int, default=0)
    parser.add_argument(
        "--export-scale",
        default="1",
        help="multiply the selected factor by this scalar before exporting",
    )
    parser.add_argument("--out", default="", help="write the selected nonconstant factor as separator JSON")
    args = parser.parse_args()

    data = json.loads(Path(args.input).read_text())
    characteristic = data["characteristic"] if args.characteristic < 0 else args.characteristic
    field = base_field(characteristic)
    ring = PolynomialRing(field, [f"z{i}" for i in range(9)] + ["q"], order="degrevlex")
    print("=" * 78)
    print("SEPARATOR CANDIDATE FACTOR")
    print(f"input: {args.input}")
    print(f"field: {'QQ' if characteristic == 0 else 'GF(%d)' % characteristic}")
    indices = range(len(data["passing_candidates"])) if args.all else [args.candidate_index]
    selected_factor_payload = None
    for index in indices:
        candidate = data["passing_candidates"][index]
        poly = ring.zero()
        for term in candidate["separator_terms_prefix"]:
            poly += field(term["coefficient"]) * ring(term["monomial"])
        factorization = poly.factor()
        print("-" * 78)
        print(f"candidate index: {index}")
        print(f"basis index: {candidate['basis_index']}")
        print(f"terms: {len(poly.dict())}")
        print(f"total degree: {poly.total_degree()}")
        print(f"degree in q: {poly.degree(ring.gen(9))}")
        print(f"|s(root)|: {candidate['separator_abs_value']}")
        print("factorization:")
        print(factorization)
        if args.out and index == args.candidate_index:
            nonconstant = [(factor, exp) for factor, exp in factorization if factor.total_degree() > 0]
            if not nonconstant:
                raise SystemExit("selected candidate has no nonconstant factor")
            if args.factor_index < 0 or args.factor_index >= len(nonconstant):
                raise SystemExit(f"--factor-index out of range for {len(nonconstant)} nonconstant factors")
            factor, exp = nonconstant[args.factor_index]
            if exp != 1:
                factor = factor ** exp
            factor *= field(args.export_scale)
            selected_factor_payload = {
                "source": args.input,
                "candidate_index": index,
                "basis_index": candidate["basis_index"],
                "field_characteristic": characteristic,
                "factor_index": args.factor_index,
                "terms": [
                    {
                        "coefficient": str(lift_coeff(coeff, characteristic)),
                        "monomial": monomial_string(ring.variable_names(), monomial),
                    }
                    for monomial, coeff in sorted(factor.dict().items())
                    if coeff
                ],
            }
    print("=" * 78)
    if args.out:
        if selected_factor_payload is None:
            raise SystemExit("--out requested, but no selected factor was exported")
        Path(args.out).parent.mkdir(parents=True, exist_ok=True)
        Path(args.out).write_text(json.dumps(selected_factor_payload, indent=1) + "\n")
        print(f"wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
