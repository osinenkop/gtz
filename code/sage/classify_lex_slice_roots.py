#!/usr/bin/env python3
"""Classify real roots reconstructed from a lex determinant slice."""
from __future__ import annotations

import argparse
import itertools
import json
import sys
from pathlib import Path

import numpy as np
from sage.all import ComplexField, PolynomialRing, QQ, RealField

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from gtz63_semialgebraic import TRIPLES  # noqa: E402
from lexify_determinant_slice import split_singular_ideal  # noqa: E402

TH = 1.0 / 6.0


def parse_rows(text: str) -> tuple[int, int, int]:
    rows = tuple(int(x.strip()) for x in text.split(",") if x.strip())
    if len(rows) != 3 or len(set(rows)) != 3 or any(row < 0 or row >= 6 for row in rows):
        raise SystemExit("--chart-rows must be three distinct row indices, e.g. 0,1,2")
    return rows


def projector_from_z(z: np.ndarray, chart_rows: tuple[int, int, int]) -> np.ndarray:
    zmat = z.reshape(3, 3)
    y = np.zeros((6, 3))
    y[list(chart_rows), :] = np.eye(3)
    complement = tuple(row for row in range(6) if row not in chart_rows)
    y[list(complement), :] = zmat
    gram = y.T @ y
    return y @ np.linalg.inv(gram) @ y.T


def lambdas(projector: np.ndarray) -> np.ndarray:
    return np.array([
        np.linalg.eigvalsh(projector[np.ix_(triple, triple)])[0]
        for triple in TRIPLES
    ])


def determinant_values(projector: np.ndarray, active: tuple[int, ...]) -> np.ndarray:
    return np.array([
        np.linalg.det(projector[np.ix_(TRIPLES[idx], TRIPLES[idx])] - TH * np.eye(3))
        for idx in active
    ])


def classify_z(z: np.ndarray, active: tuple[int, ...], chart_rows: tuple[int, int, int]) -> dict:
    projector = projector_from_z(z, chart_rows)
    lam = lambdas(projector)
    active_set = set(active)
    inactive = [idx for idx in range(len(TRIPLES)) if idx not in active_set]
    active_min_shift = min(float(lam[idx] - TH) for idx in active)
    inactive_max_shift = max(float(lam[idx] - TH) for idx in inactive) if inactive else -float("inf")
    actual_active = [idx for idx, value in enumerate(lam) if abs(float(value - TH)) <= 1e-7]
    return {
        "z": [float(x) for x in z],
        "det_residual_max": float(np.max(np.abs(determinant_values(projector, active)))) if active else 0.0,
        "F_minus_1_6": float(np.max(lam) - TH),
        "active_min_lambda_minus_1_6": active_min_shift,
        "inactive_max_lambda_minus_1_6": inactive_max_shift,
        "passes_active_psd": bool(active_min_shift >= -1e-7),
        "passes_inactive": bool(inactive_max_shift <= 1e-7),
        "passes_equality_inequalities": bool(active_min_shift >= -1e-7 and inactive_max_shift <= 1e-7),
        "actual_active": actual_active,
        "actual_active_triples": [list(TRIPLES[idx]) for idx in actual_active],
        "actual_active_size": len(actual_active),
        "selected_equals_actual": bool(set(actual_active) == active_set),
    }


def read_polynomials(path: Path, variables: tuple[str, ...]):
    ring = PolynomialRing(QQ, variables)
    return ring, [ring(text.replace("^", "**")) for text in split_singular_ideal(path.read_text())]


def find_univariate(polynomials, variable_name: str):
    matches = []
    for index, poly in enumerate(polynomials, start=1):
        names = {str(var) for var in poly.variables()}
        if names and names <= {variable_name}:
            matches.append((index, poly))
    if len(matches) != 1:
        raise ValueError(f"expected one univariate polynomial in {variable_name}, got {len(matches)}")
    return matches[0]


def variable_names(expr) -> set[str]:
    variables = getattr(expr, "variables", None)
    if variables is None:
        return set()
    return {str(var) for var in variables()}


def solve_shape_relations(polynomials, ring, parameter_name: str) -> dict[str, object]:
    gens = {str(gen): gen for gen in ring.gens()}
    parameter = gens[parameter_name]
    relations = {}

    def unique_candidates(candidates):
        unique = []
        seen = set()
        for candidate in candidates:
            key = str(candidate)
            if key not in seen:
                seen.add(key)
                unique.append(candidate)
        return unique

    for name, gen in gens.items():
        if name == parameter_name:
            continue
        for poly in polynomials:
            names = {str(var) for var in poly.variables()}
            if names == {name} and poly.degree(gen) == 1:
                coeff = poly.monomial_coefficient(gen)
                rest = poly - coeff * gen
                relations[name] = -rest / coeff
                break
            if names == {name} and poly.degree(gen) == 2 and poly.monomial_coefficient(gen**2) != 0:
                rest = poly - poly.monomial_coefficient(gen**2) * gen**2
                if rest == 0:
                    relations[name] = ring.zero()
                    break

    remaining = [name for name in gens if name != parameter_name and name not in relations]
    while remaining:
        progress = False
        known = set(relations)
        for name in list(remaining):
            gen = gens[name]
            candidates = []
            for poly in polynomials:
                names = {str(var) for var in poly.variables()}
                if name not in names or poly.degree(gen) != 1:
                    continue
                if not names <= ({name, parameter_name} | known):
                    continue
                expr = poly
                for rel_name, rel_expr in relations.items():
                    expr = expr.subs({gens[rel_name]: rel_expr})
                coeff = expr.monomial_coefficient(gen)
                rest = expr - coeff * gen
                coeff_names = variable_names(coeff)
                rest_names = variable_names(rest)
                if coeff_names <= {parameter_name} and rest_names <= {parameter_name}:
                    candidates.append(-rest / coeff)
            candidates = unique_candidates(candidates)
            if len(candidates) == 1:
                relations[name] = candidates[0]
                remaining.remove(name)
                progress = True
        if not progress:
            unresolved = {name: len([
                poly for poly in polynomials
                if name in {str(var) for var in poly.variables()}
            ]) for name in remaining}
            raise ValueError(f"could not solve triangular shape relations: {unresolved}")
    return relations


def factor_summary(poly, known_quadratic) -> list[dict]:
    out = []
    for factor, multiplicity in poly.factor():
        record = {
            "degree": int(factor.degree()),
            "multiplicity": int(multiplicity),
            "is_known_quadratic": bool(factor.monic() == known_quadratic.monic()),
        }
        if factor.degree() <= 2:
            record["factor"] = str(factor)
        out.append(record)
    return out


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--lex-json", required=True)
    parser.add_argument("--out")
    parser.add_argument("--parameter", default="z0")
    parser.add_argument("--chart-rows", default="0,1,2")
    parser.add_argument("--real-prec", type=int, default=100)
    args = parser.parse_args()

    lex_json_path = Path(args.lex_json)
    data = json.loads(lex_json_path.read_text())
    basis_path = Path(data["lex_basis"]["path"])
    variables = tuple(data["summary"]["variables"])
    active = tuple(data["summary"]["active_indices"])
    chart_rows = parse_rows(args.chart_rows)

    ring, polynomials = read_polynomials(basis_path, variables)
    parameter_index, parameter_poly = find_univariate(polynomials, args.parameter)
    relations = solve_shape_relations(polynomials, ring, args.parameter)
    gens = {str(gen): gen for gen in ring.gens()}
    parameter = gens[args.parameter]
    univariate_ring = PolynomialRing(QQ, args.parameter)
    univariate_parameter = univariate_ring.gen()
    parameter_univariate = univariate_ring(str(parameter_poly).replace("^", "**"))
    known_quadratic = univariate_parameter**2 - QQ(5) / QQ(9)

    rr = RealField(args.real_prec)
    cc = ComplexField(args.real_prec)
    real_roots = parameter_univariate.roots(rr, multiplicities=True)
    complex_roots = parameter_univariate.roots(cc, multiplicities=True)

    roots = []
    z_names = [f"z{i}" for i in range(9)]
    for root_index, (root, multiplicity) in enumerate(real_roots, start=1):
        values = {args.parameter: root}
        for name, expr in relations.items():
            values[name] = expr.subs({parameter: root})
        z = np.array([float(values[name]) for name in z_names], dtype=float)
        record = {
            "root_index": root_index,
            "parameter": float(root),
            "multiplicity": int(multiplicity),
            **classify_z(z, active, chart_rows),
        }
        roots.append(record)

    result = {
        "source": str(lex_json_path),
        "basis_path": str(basis_path),
        "parameter": args.parameter,
        "parameter_polynomial_index": parameter_index,
        "parameter_polynomial_degree": int(parameter_univariate.degree()),
        "known_quadratic_divides": bool(parameter_univariate % known_quadratic == 0),
        "factor_summary": factor_summary(parameter_univariate, known_quadratic),
        "complex_root_count_with_multiplicity": int(sum(m for _root, m in complex_roots)),
        "real_root_count_with_multiplicity": int(sum(m for _root, m in real_roots)),
        "real_roots": roots,
        "passes_equality_inequalities": int(sum(root["passes_equality_inequalities"] for root in roots)),
        "passes_active_psd": int(sum(root["passes_active_psd"] for root in roots)),
        "passes_inactive": int(sum(root["passes_inactive"] for root in roots)),
    }

    out_path = Path(args.out) if args.out else lex_json_path.with_name(f"{lex_json_path.stem}_classified.json")
    out_path.write_text(json.dumps(result, indent=1) + "\n")
    print(
        "classified",
        result["real_root_count_with_multiplicity"],
        "real roots;",
        result["passes_equality_inequalities"],
        "pass inequalities",
    )
    print(f"wrote {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
