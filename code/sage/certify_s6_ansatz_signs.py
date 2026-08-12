#!/usr/bin/env python3
"""Exact real-sign classification for structured size-6 ansatz lex branches."""
from __future__ import annotations

import argparse
import itertools
import json
import re
from pathlib import Path

from sage.all import AA, Matrix, PolynomialRing, QQ

from gtz63_semialgebraic import TRIPLES
from lexify_s6_ansatz_basis import split_singular_ideal


VARIABLES = ("a", "b", "c", "d", "q")
ACTIVE = (0, 1, 9, 15, 16, 17)


def parse_compact_singular_poly(ring, text: str):
    """Parse Singular's compact output, e.g. ``3d2q-5`` as ``3*d^2*q-5``."""
    text = text.strip().replace(" ", "")
    if not text:
        return ring.zero()
    if text[0] not in "+-":
        text = "+" + text
    gens = {str(gen): gen for gen in ring.gens()}
    out = ring.zero()
    for match in re.finditer(r"([+-])([^+-]+)", text):
        sign, body = match.groups()
        coeff_match = re.match(r"(\d+(?:/\d+)?)?(.*)$", body)
        if coeff_match is None:
            raise ValueError(f"cannot parse term {body!r}")
        coeff_text, monomial = coeff_match.groups()
        coeff = QQ(coeff_text) if coeff_text else QQ(1)
        if sign == "-":
            coeff = -coeff
        term = ring(coeff)
        i = 0
        while i < len(monomial):
            name = monomial[i]
            if name not in gens:
                raise ValueError(f"unknown variable {name!r} in {text!r}")
            i += 1
            j = i
            while j < len(monomial) and monomial[j].isdigit():
                j += 1
            exponent = int(monomial[i:j] or "1")
            term *= gens[name] ** exponent
            i = j
        out += term
    return out


def read_lex_polynomials(path: Path):
    ring = PolynomialRing(QQ, VARIABLES, order="lex")
    polynomials = [
        parse_compact_singular_poly(ring, item)
        for item in split_singular_ideal(path.read_text().strip())
    ]
    return ring, polynomials


def projector_numerator(a, b, c, d):
    y = Matrix(
        AA,
        [
            [1, 0, 0],
            [0, 1, 0],
            [0, 0, 1],
            [-a, a, -b],
            [c, d, -a],
            [-d, -c, -a],
        ],
    )
    gram = y.transpose() * y
    det_gram = gram.det()
    return y * gram.adjugate() * y.transpose(), det_gram


def shifted_block(numer, det_gram, triple_index: int):
    triple = TRIPLES[triple_index]
    return Matrix(
        AA,
        3,
        3,
        lambda i, j: 6 * numer[triple[i], triple[j]] - (det_gram if i == j else AA(0)),
    )


def principal_minor(mat, rows: tuple[int, ...]):
    return mat.matrix_from_rows_and_columns(rows, rows).det()


def decimal(value, digits: int = 30) -> str:
    return value.n(digits=digits).str()


def first_active_failure(numer, det_gram):
    for active_index in ACTIVE:
        mat = shifted_block(numer, det_gram, active_index)
        for size in (1, 2, 3):
            for rows in itertools.combinations(range(3), size):
                value = principal_minor(mat, rows)
                if value < 0:
                    return {
                        "active_index": active_index,
                        "active_triple": list(TRIPLES[active_index]),
                        "rows": list(rows),
                        "minor_value": str(value),
                        "minor_decimal": decimal(value),
                    }
    return None


def first_inactive_pd_witness(numer, det_gram):
    active_set = set(ACTIVE)
    for index, triple in enumerate(TRIPLES):
        if index in active_set:
            continue
        mat = shifted_block(numer, det_gram, index)
        leading = [
            principal_minor(mat, (0,)),
            principal_minor(mat, (0, 1)),
            principal_minor(mat, (0, 1, 2)),
        ]
        if all(value > 0 for value in leading):
            return {
                "inactive_index": index,
                "inactive_triple": list(triple),
                "leading_minors": [str(value) for value in leading],
                "leading_minors_decimal": [decimal(value) for value in leading],
            }
    return None


def certify_branch(path: Path) -> dict:
    ring, polynomials = read_lex_polynomials(path)
    ideal = ring.ideal(polynomials)
    points = ideal.variety(ring=AA)
    records = []
    counts = {
        "real_points": len(points),
        "active_failure": 0,
        "inactive_pd_failure": 0,
        "unclassified": 0,
    }
    for point in points:
        values = {name: AA(point[ring(name)]) for name in VARIABLES}
        numer, det_gram = projector_numerator(values["a"], values["b"], values["c"], values["d"])
        if det_gram <= 0:
            raise ArithmeticError(f"nonpositive chart determinant at {path}: {det_gram}")
        active_failure = first_active_failure(numer, det_gram)
        inactive_witness = None
        if active_failure is None:
            inactive_witness = first_inactive_pd_witness(numer, det_gram)
        if active_failure is not None:
            counts["active_failure"] += 1
        elif inactive_witness is not None:
            counts["inactive_pd_failure"] += 1
        else:
            counts["unclassified"] += 1
        records.append(
            {
                "q": str(values["q"]),
                "q_decimal": decimal(values["q"]),
                "a_decimal": decimal(values["a"], digits=20),
                "b_decimal": decimal(values["b"], digits=20),
                "c_decimal": decimal(values["c"], digits=20),
                "d_decimal": decimal(values["d"], digits=20),
                "det_gram_decimal": decimal(det_gram, digits=20),
                "active_failure": active_failure,
                "inactive_pd_failure": inactive_witness,
                "certified_infeasible": active_failure is not None or inactive_witness is not None,
            }
        )
    counts["certified_infeasible"] = counts["active_failure"] + counts["inactive_pd_failure"]
    return {
        "lex_basis_path": str(path),
        "polynomial_count": len(polynomials),
        "ideal_dimension": int(ideal.dimension()),
        "ideal_vdim": int(ideal.vector_space_dimension()),
        "counts": counts,
        "roots": records,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--lex-basis",
        action="append",
        required=True,
        help="Lex basis text file generated by lexify_s6_ansatz_basis.py.",
    )
    parser.add_argument("--out", default="code/sage/out/s6_ansatz_cplusd_sign_cert_QQ.json")
    args = parser.parse_args()

    branches = [certify_branch(Path(path)) for path in args.lex_basis]
    total = {
        "real_points": sum(branch["counts"]["real_points"] for branch in branches),
        "active_failure": sum(branch["counts"]["active_failure"] for branch in branches),
        "inactive_pd_failure": sum(branch["counts"]["inactive_pd_failure"] for branch in branches),
        "unclassified": sum(branch["counts"]["unclassified"] for branch in branches),
        "certified_infeasible": sum(
            branch["counts"]["certified_infeasible"] for branch in branches
        ),
    }
    payload = {
        "active": list(ACTIVE),
        "active_triples": [list(TRIPLES[index]) for index in ACTIVE],
        "branches": branches,
        "total": total,
    }
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(payload, indent=1) + "\n")

    print("=" * 78)
    print("S6 ANSATZ EXACT SIGN CERTIFICATION")
    for branch in branches:
        counts = branch["counts"]
        print(
            f"{branch['lex_basis_path']}: real={counts['real_points']} "
            f"active_fail={counts['active_failure']} "
            f"inactive_pd={counts['inactive_pd_failure']} "
            f"unclassified={counts['unclassified']}"
        )
    print(f"total: {total}")
    print(f"wrote {out}")
    print("=" * 78)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
