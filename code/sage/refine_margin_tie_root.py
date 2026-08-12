#!/usr/bin/env python3
"""High-precision Newton refinement for active-plus-tie determinant roots.

This is the "specialized solver" fallback for the low-active margin route.  It
does not try to Groebner the whole ideal.  Instead it takes a numerical v30
minimizer, converts it to the standard chart, and refines the square system

    f_T(z) = 0          for T in A,
    h_S(z,q) = 0        for outside tied triples S,

where f_T=det(6*N_TT-d*I)/d^2 and h_S=det(6*N_SS-q*d*I)/d^2.
The variable q is six times the tied outside eigenvalue.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
from sage.all import Matrix, PolynomialRing, QQ, RDF, RealField, vector

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from gtz63_semialgebraic import (  # noqa: E402
    TRIPLES,
    active_determinant,
    parse_indices,
    standard_chart,
)
from probe_margin_tie_system import (  # noqa: E402
    resolve_active,
    shifted_determinant,
)


def retract(X):
    U, _, Vt = np.linalg.svd(X, full_matrices=False)
    return U @ Vt


def z_from_v30(path: str):
    data = json.loads(Path(path).read_text())
    A = retract(np.array(data["best_point"], dtype=float).reshape(6, 3))
    top = A[:3, :]
    det = float(np.linalg.det(top))
    if abs(det) < 1e-10:
        raise ValueError(f"top chart is ill-conditioned, det={det}")
    Y = A @ np.linalg.inv(top)
    return Y[3:, :].reshape(9), data


def build_system(active, ties, order):
    names = [f"z{i}" for i in range(9)] + ["q"]
    ring = PolynomialRing(QQ, names, order=order)
    gens = dict(zip(ring.variable_names(), ring.gens()))
    chart = standard_chart(ring, include_h=False)
    q = gens["q"]
    polynomials = [active_determinant(chart, idx) for idx in active]
    polynomials.extend(shifted_determinant(chart, idx, q) for idx in ties)
    variables = list(chart.z) + [q]
    jacobian = [[poly.derivative(var) for var in variables] for poly in polynomials]
    return ring, chart, variables, polynomials, jacobian


def eval_vector(polynomials, x, RF):
    return vector(RF, [RF(poly(*x)) for poly in polynomials])


def eval_matrix(polynomials, x, RF):
    return Matrix(RF, len(polynomials), len(x), [RF(poly(*x)) for row in polynomials for poly in row])


def max_abs(vec):
    return max(abs(v) for v in vec) if len(vec) else 0


def coeff_height(poly):
    coeffs = poly.coefficients(sparse=False)
    try:
        return max(abs(int(c)) for c in coeffs)
    except TypeError:
        return None


def algdep_screen(x, max_degree, max_height):
    rows = []
    for deg in range(1, max_degree + 1):
        try:
            poly = x.algdep(deg)
        except Exception as exc:  # pragma: no cover - diagnostic path
            rows.append({"degree": deg, "error": repr(exc)})
            continue
        height = coeff_height(poly)
        value = abs(poly(x))
        rows.append({
            "degree": deg,
            "polynomial": str(poly),
            "height": height,
            "residual": str(value),
            "accepted": bool(height is not None and height <= max_height),
        })
    return rows


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", default="verify/out/v30_margin_7_14113.json")
    parser.add_argument("--active-indices")
    parser.add_argument("--size", type=int)
    parser.add_argument("--canon", type=int)
    parser.add_argument("--tie-indices", required=True)
    parser.add_argument("--precision", type=int, default=180)
    parser.add_argument("--steps", type=int, default=20)
    parser.add_argument("--tol-exp", type=int, default=120)
    parser.add_argument("--algdep-degree", type=int, default=12)
    parser.add_argument("--algdep-height", type=int, default=10**9)
    parser.add_argument("--order", default="degrevlex")
    parser.add_argument("--out", default="")
    args = parser.parse_args()

    active = resolve_active(args)
    ties = parse_indices(args.tie_indices)
    if len(active) + len(ties) != 10:
        raise SystemExit("this Newton solver expects a square 10-equation system")

    z0, source = z_from_v30(args.input)
    q0 = 6.0 * float(source["best"][0]["F"])
    ring, chart, variables, polynomials, jacobian = build_system(active, ties, args.order)

    RF = RealField(args.precision)
    x = vector(RF, [RF(str(v)) for v in z0] + [RF(str(q0))])
    target = RF(10) ** (-args.tol_exp)
    history = []
    for step in range(args.steps):
        F = eval_vector(polynomials, x, RF)
        J = eval_matrix(jacobian, x, RF)
        fnorm = max_abs(F)
        try:
            delta = J.solve_right(-F)
        except Exception as exc:
            raise RuntimeError(f"Newton solve failed at step {step}: {exc}") from exc
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
    algdep_q = algdep_screen(q, args.algdep_degree, args.algdep_height)
    algdep_margin = algdep_screen(margin, args.algdep_degree, args.algdep_height)

    payload = {
        "input": args.input,
        "active": list(active),
        "active_triples": [list(TRIPLES[i]) for i in active],
        "ties": list(ties),
        "tie_triples": [list(TRIPLES[i]) for i in ties],
        "precision_bits": args.precision,
        "steps_requested": args.steps,
        "steps_used": len(history),
        "history": history,
        "residual_inf": str(max_abs(F)),
        "q": str(q),
        "margin": str(margin),
        "q_decimal_30": f"{float(q):.30g}",
        "margin_decimal_30": f"{float(margin):.30g}",
        "jacobian_singular_values_float": singular_values,
        "jacobian_condition_float": float(max(singular_values) / min(singular_values)),
        "algdep_q": algdep_q,
        "algdep_margin": algdep_margin,
        "z": [str(v) for v in x[:-1]],
    }

    out = args.out or f"code/sage/out/refine_tie_{len(active)}_{source.get('canon','x')}.json"
    Path(out).parent.mkdir(parents=True, exist_ok=True)
    Path(out).write_text(json.dumps(payload, indent=1) + "\n")

    print("=" * 78)
    print("MARGIN TIE ROOT REFINEMENT")
    print(f"active: {list(active)}")
    print(f"ties:   {list(ties)}")
    print(f"steps:  {len(history)}")
    print(f"residual_inf: {payload['residual_inf']}")
    print(f"q:      {payload['q_decimal_30']}")
    print(f"margin: {payload['margin_decimal_30']}")
    print(f"cond(J) float: {payload['jacobian_condition_float']:.3e}")
    accepted = [row for row in algdep_q if row.get("accepted")]
    if accepted:
        print("accepted q algdep candidates:")
        for row in accepted:
            print(f"  deg {row['degree']} height {row['height']}: {row['polynomial']}")
    else:
        print("no low-height q algdep candidate accepted")
    print(f"wrote {out}")
    print("=" * 78)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
