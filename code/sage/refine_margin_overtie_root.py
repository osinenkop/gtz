#!/usr/bin/env python3
"""High-precision Gauss-Newton refinement for over-tied margin roots.

This complements ``refine_margin_tie_root.py``.  The square solver handles
``|active|+|ties|=10``.  Over-tied numerical minimizers have more tie equations
than variables, so this script solves the full active-plus-all-ties system in a
least-squares/Gauss-Newton sense and reports whether the residual actually drops
to high precision with full column-rank Jacobian.

The output is diagnostic evidence, not an interval certificate.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from sage.all import Matrix, RDF, RealField, vector

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from gtz63_semialgebraic import TRIPLES, parse_indices  # noqa: E402
from probe_margin_tie_system import resolve_active  # noqa: E402
from refine_margin_tie_root import build_system, eval_matrix, eval_vector, max_abs, z_from_v30  # noqa: E402


def gauss_newton_step(J, F):
    normal = J.transpose() * J
    rhs = -(J.transpose() * F)
    return normal.solve_right(rhs)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--active-indices")
    parser.add_argument("--size", type=int)
    parser.add_argument("--canon", type=int)
    parser.add_argument("--tie-indices", required=True)
    parser.add_argument("--precision", type=int, default=500)
    parser.add_argument("--steps", type=int, default=30)
    parser.add_argument("--tol-exp", type=int, default=220)
    parser.add_argument("--order", default="degrevlex")
    parser.add_argument("--out", default="")
    args = parser.parse_args()

    active = resolve_active(args)
    ties = parse_indices(args.tie_indices)
    if len(active) + len(ties) < 10:
        raise SystemExit("this over-tie refiner expects at least 10 equations")

    z0, source = z_from_v30(args.input)
    q0 = 6.0 * float(source["best"][0]["F"])
    _ring, _chart, _variables, polynomials, jacobian = build_system(active, ties, args.order)

    RF = RealField(args.precision)
    x = vector(RF, [RF(str(v)) for v in z0] + [RF(str(q0))])
    target = RF(10) ** (-args.tol_exp)
    history = []
    for step in range(args.steps):
        F = eval_vector(polynomials, x, RF)
        J = eval_matrix(jacobian, x, RF)
        fnorm = max_abs(F)
        delta = gauss_newton_step(J, F)
        dnorm = max_abs(delta)
        history.append({"step": step, "residual_inf": str(fnorm), "delta_inf": str(dnorm)})
        x = x + delta
        if fnorm < target and dnorm < target:
            break

    F = eval_vector(polynomials, x, RF)
    J = eval_matrix(jacobian, x, RF)
    q = x[-1]
    margin = (q - 1) / 6
    singular_values = [float(v) for v in Matrix(RDF, J).singular_values()]
    min_sv = min(singular_values) if singular_values else None
    max_sv = max(singular_values) if singular_values else None
    residuals = [str(v) for v in F]
    labels = [f"active:{i}:{''.join(str(a) for a in TRIPLES[i])}" for i in active]
    labels.extend(f"tie:{i}:{''.join(str(a) for a in TRIPLES[i])}" for i in ties)

    payload = {
        "input": args.input,
        "active": list(active),
        "active_triples": [list(TRIPLES[i]) for i in active],
        "ties": list(ties),
        "tie_triples": [list(TRIPLES[i]) for i in ties],
        "equation_count": len(polynomials),
        "variable_count": 10,
        "precision_bits": args.precision,
        "steps_requested": args.steps,
        "steps_used": len(history),
        "history": history,
        "residual_inf": str(max_abs(F)),
        "residuals": [{"label": label, "value": value} for label, value in zip(labels, residuals)],
        "q": str(q),
        "margin": str(margin),
        "q_decimal_30": f"{float(q):.30g}",
        "margin_decimal_30": f"{float(margin):.30g}",
        "jacobian_singular_values_float": singular_values,
        "jacobian_min_singular_value_float": min_sv,
        "jacobian_condition_float": float(max_sv / min_sv) if min_sv else None,
        "full_column_rank_numeric": bool(min_sv is not None and min_sv > 1e-10),
        "z": [str(v) for v in x[:-1]],
    }

    out = args.out or f"code/sage/out/refine_overtie_{len(active)}_{source.get('canon','x')}.json"
    Path(out).parent.mkdir(parents=True, exist_ok=True)
    Path(out).write_text(json.dumps(payload, indent=1) + "\n")

    print("=" * 78)
    print("MARGIN OVERTIE ROOT REFINEMENT")
    print(f"active: {list(active)}")
    print(f"ties:   {list(ties)}")
    print(f"equations/variables: {len(polynomials)}/10")
    print(f"steps:  {len(history)}")
    print(f"residual_inf: {payload['residual_inf']}")
    print(f"q:      {payload['q_decimal_30']}")
    print(f"margin: {payload['margin_decimal_30']}")
    print(f"min sv / cond(J) float: {min_sv:.3e} / {payload['jacobian_condition_float']:.3e}")
    print(f"wrote {out}")
    print("=" * 78)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
