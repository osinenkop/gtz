#!/usr/bin/env python3
"""Recover a rational q-factor by matching modular factor products of fixed degree."""
from __future__ import annotations

import argparse
import json
import math
import re
from dataclasses import dataclass
from pathlib import Path

from sage.all import GF, PolynomialRing


@dataclass(frozen=True)
class ProductCandidate:
    degree: int
    indices: tuple[int, ...]
    coeffs: tuple[int, ...]


@dataclass(frozen=True)
class BeamState:
    primes: tuple[int, ...]
    choices: tuple[tuple[int, ...], ...]
    coeffs: tuple[int, ...]
    modulus: int
    max_abs_coeff: int
    sum_abs_coeff: int


def parse_prime(path: Path, data: dict) -> int:
    characteristic = int(data.get("characteristic") or 0)
    if characteristic:
        return characteristic
    match = re.search(r"_p(\d+)", path.name)
    if not match:
        raise ValueError(f"cannot infer prime from {path}")
    return int(match.group(1))


def parse_int_list(text: str) -> tuple[int, ...]:
    values = []
    for part in text.split(","):
        part = part.strip()
        if part:
            values.append(int(part))
    return tuple(values)


def monic_linear_root_mod_prime(poly_text: str, prime: int) -> int | None:
    poly = poly_text.replace(" ", "")
    if poly == "q":
        return 0
    match = re.fullmatch(r"q\+([0-9]+)", poly)
    if match:
        return (-int(match.group(1))) % prime
    match = re.fullmatch(r"q-([0-9]+)", poly)
    if match:
        return int(match.group(1)) % prime
    return None


def center_coeff(value: int, modulus: int) -> int:
    value %= modulus
    if value > modulus // 2:
        value -= modulus
    return value


def center_coeffs(coeffs: tuple[int, ...], modulus: int) -> tuple[int, ...]:
    return tuple(center_coeff(int(c), modulus) for c in coeffs)


def combine_centered_coeffs(
    coeffs: tuple[int, ...], modulus: int, residues: tuple[int, ...], prime: int
) -> tuple[int, ...]:
    inv_modulus = pow(modulus % prime, -1, prime)
    next_modulus = modulus * prime
    half = next_modulus // 2
    out = []
    for coeff, residue in zip(coeffs, residues):
        step = ((int(residue) - (coeff % prime)) * inv_modulus) % prime
        value = coeff + modulus * step
        if value > half:
            value -= next_modulus
        out.append(value)
    return tuple(out)


def score_centered_coeffs(coeffs: tuple[int, ...]) -> tuple[int, int]:
    return max(abs(c) for c in coeffs), sum(abs(c) for c in coeffs)


def read_expanded_factors(path: Path, variable: str, excluded_roots: tuple[int, ...]):
    data = json.loads(path.read_text())
    prime = parse_prime(path, data)
    ring = PolynomialRing(GF(prime), variable)
    excluded = {root % prime for root in excluded_roots}
    entries = []
    removed = []
    for index, item in enumerate(data["factors"]):
        poly = ring(str(item["polynomial"]).replace("^", "**"))
        lead = poly.leading_coefficient()
        if lead != 1:
            poly = poly / lead
        root = (
            monic_linear_root_mod_prime(str(item["polynomial"]), prime)
            if int(item["degree"]) == 1
            else None
        )
        if root is not None and root in excluded:
            removed.append(
                {
                    "index": index,
                    "root": int(root),
                    "degree": int(item["degree"]),
                    "multiplicity": int(item["multiplicity"]),
                    "polynomial": item["polynomial"],
                }
            )
            continue
        for _ in range(int(item["multiplicity"])):
            entries.append((int(item["degree"]), poly))
    return prime, ring, entries, removed


def coeff_tuple(poly, q, degree: int) -> tuple[int, ...]:
    lead = poly.leading_coefficient()
    if lead != 1:
        poly = poly / lead
    return tuple(int(poly.monomial_coefficient(q**i)) for i in range(degree + 1))


def enumerate_products(
    path: Path,
    target_degree: int,
    variable: str,
    max_candidates: int,
    excluded_roots: tuple[int, ...],
):
    prime, ring, entries, removed = read_expanded_factors(path, variable, excluded_roots)
    q = ring.gen()
    entries = tuple(entries)
    suffix = [0] * (len(entries) + 1)
    for i in range(len(entries) - 1, -1, -1):
        suffix[i] = suffix[i + 1] + entries[i][0]

    out: list[ProductCandidate] = []

    def rec(pos: int, remaining: int, chosen: list[int], product):
        if len(out) >= max_candidates:
            return
        if remaining == 0:
            out.append(
                ProductCandidate(
                    target_degree,
                    tuple(chosen),
                    coeff_tuple(product, q, target_degree),
                )
            )
            return
        if pos >= len(entries) or remaining < 0 or suffix[pos] < remaining:
            return

        degree, factor = entries[pos]
        if degree <= remaining:
            chosen.append(pos)
            rec(pos + 1, remaining - degree, chosen, product * factor)
            chosen.pop()
        rec(pos + 1, remaining, chosen, product)

    rec(0, target_degree, [], ring.one())
    return prime, out, len(entries), removed


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--degree", type=int, required=True)
    parser.add_argument("--variable", default="q")
    parser.add_argument("--beam", type=int, default=200)
    parser.add_argument("--max-candidates-per-prime", type=int, default=20000)
    parser.add_argument("--exclude-rational-roots", default="")
    parser.add_argument("--sort-by-candidate-count", action="store_true")
    parser.add_argument("--out", required=True)
    parser.add_argument("factor_jsons", nargs="+")
    args = parser.parse_args()

    per_prime = []
    excluded_roots = parse_int_list(args.exclude_rational_roots)
    for factor_json in args.factor_jsons:
        path = Path(factor_json)
        prime, candidates, expanded_count, removed = enumerate_products(
            path,
            args.degree,
            args.variable,
            args.max_candidates_per_prime,
            excluded_roots,
        )
        if not candidates:
            raise SystemExit(f"no degree-{args.degree} products in {path}")
        per_prime.append(
            {
                "path": str(path),
                "prime": prime,
                "expanded_factor_count": expanded_count,
                "removed": removed,
                "candidate_count": len(candidates),
                "candidates": candidates,
            }
        )
        print(
            f"{path}: prime={prime} expanded_factors={expanded_count} "
            f"degree-{args.degree} candidates={len(candidates)}",
            flush=True,
        )
    if args.sort_by_candidate_count:
        per_prime.sort(key=lambda record: (record["candidate_count"], record["prime"]))

    first = per_prime[0]
    states = []
    for candidate in first["candidates"]:
        coeffs = center_coeffs(candidate.coeffs, first["prime"])
        max_abs, sum_abs = score_centered_coeffs(coeffs)
        states.append(
            BeamState(
                primes=(first["prime"],),
                choices=(candidate.indices,),
                coeffs=coeffs,
                modulus=first["prime"],
                max_abs_coeff=max_abs,
                sum_abs_coeff=sum_abs,
            )
        )
    states.sort(key=lambda item: (item.max_abs_coeff, item.sum_abs_coeff))
    states = states[: args.beam]

    for record in per_prime[1:]:
        next_states: list[BeamState] = []
        for state in states:
            for candidate in record["candidates"]:
                primes = (*state.primes, record["prime"])
                coeffs = combine_centered_coeffs(
                    state.coeffs, state.modulus, candidate.coeffs, record["prime"]
                )
                max_abs, sum_abs = score_centered_coeffs(coeffs)
                next_states.append(
                    BeamState(
                        primes=primes,
                        choices=(*state.choices, candidate.indices),
                        coeffs=coeffs,
                        modulus=state.modulus * record["prime"],
                        max_abs_coeff=max_abs,
                        sum_abs_coeff=sum_abs,
                    )
                )
        next_states.sort(key=lambda item: (item.max_abs_coeff, item.sum_abs_coeff))
        states = next_states[: args.beam]
        print(
            f"after prime {record['prime']}: best max_abs={states[0].max_abs_coeff} "
            f"log10={math.log10(max(1, states[0].max_abs_coeff)):.3f}",
            flush=True,
        )

    best = states[0]
    coeffs = list(best.coeffs)
    payload = {
        "variable": args.variable,
        "degree": args.degree,
        "beam": args.beam,
        "max_candidates_per_prime": args.max_candidates_per_prime,
        "excluded_rational_roots": list(excluded_roots),
        "sort_by_candidate_count": bool(args.sort_by_candidate_count),
        "inputs": [
            {
                "path": record["path"],
                "prime": record["prime"],
                "expanded_factor_count": record["expanded_factor_count"],
                "removed": record["removed"],
                "candidate_count": record["candidate_count"],
            }
            for record in per_prime
        ],
        "best": [
            {
                "rank": rank,
                "primes": list(state.primes),
                "choices": [list(choice) for choice in state.choices],
                "max_abs_coeff": state.max_abs_coeff,
                "sum_abs_coeff": state.sum_abs_coeff,
                "log10_max_abs_coeff": math.log10(max(1, state.max_abs_coeff)),
            }
            for rank, state in enumerate(states[: min(20, len(states))])
        ],
        "best_integer_coefficients_low_to_high": coeffs,
    }
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(payload, indent=1) + "\n")
    print("=" * 78)
    print("RECOVER Q FACTOR BY DEGREE")
    print(f"degree: {args.degree}")
    print(f"primes: {list(best.primes)}")
    print(f"best max abs coeff: {best.max_abs_coeff}")
    print(f"best log10 max abs coeff: {math.log10(max(1, best.max_abs_coeff)):.3f}")
    print(f"wrote {out}")
    print("=" * 78)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
