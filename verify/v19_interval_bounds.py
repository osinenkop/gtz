#!/usr/bin/env python3
"""
v19_interval_bounds.py -- first interval-arithmetic primitives for a certified
branch-and-bound attack on GTZ(6,3).

The target inequality is

    F(P) := max_{|T|=3} lambda_min(P_TT) >= 1/6

for every rank-3 orthogonal projector P in R^6.  This script provides two
certified building blocks:

  1.  An Arb enclosure for P(Z) = Y (Y^T Y)^{-1} Y^T on an affine Grassmann
      chart Y = [I; Z].
  2.  Certified lower bounds for lambda_min of every 3x3 interval block, using
      Arb eigenvalue enclosures, centered Weyl/Frobenius bounds, and Gershgorin
      as a transparent fallback.
  3.  An optional chart-Lipschitz enclosure for P(Z) using the exact derivative
      of the graph-chart projector map.

The branch-and-bound driver here is intentionally a smoke-test scaffold.  It can
certify boxes where interval overestimation is mild, and it records unresolved
boxes rather than treating numerical evidence as proof.  Boxes containing sharp
equality points should not be expected to certify by entrywise interval spectral
bounds alone; they need to be paired with the exact local sharpness certificates
from v6/v9/v11 and a quantitative radius argument.
"""
from __future__ import annotations

import argparse
import heapq
import itertools
import json
import math
import os
import sys
from dataclasses import dataclass, field
from typing import Iterable

import numpy as np
from flint import arb, arb_mat, ctx

sys.path.append(os.path.dirname(__file__))
from v17_active_orbits import TRIPLES, TTSP_CASES, ttsp_projector

TH_FLOAT = 1.0 / 6.0
SIX = arb(1) / 6
SQRT19 = math.sqrt(19.0)


def aconst(x: float | int | str) -> arb:
    """Create an Arb ball from a decimal string, avoiding repr surprises."""
    return arb(str(x))


def aball(center: float, radius: float) -> arb:
    if radius == 0:
        return aconst(center)
    return aconst(center) + arb(f"0 +/- {radius:.18e}")


def radius_ball(radius: arb) -> arb:
    return arb("+/- " + radius.upper().str(30, radius=False))


def eye_mat(n: int) -> arb_mat:
    out = arb_mat(n, n)
    for i in range(n):
        out[i, i] = arb(1)
    return out


def arb_lower_number(x: arb) -> arb:
    """A point Arb ball at a certified lower endpoint of x."""
    return x.lower()


def arb_upper_number(x: arb) -> arb:
    """A point Arb ball at a certified upper endpoint of x."""
    return x.upper()


def arb_min(values: Iterable[arb]) -> arb:
    it = iter(values)
    try:
        out = next(it)
    except StopIteration as exc:
        raise ValueError("empty min") from exc
    for value in it:
        out = out.min(value)
    return out


def arb_max(values: Iterable[arb]) -> arb:
    it = iter(values)
    try:
        out = next(it)
    except StopIteration as exc:
        raise ValueError("empty max") from exc
    for value in it:
        out = out.max(value)
    return out


def arb_to_json(x: arb, digits: int = 24) -> dict[str, str]:
    return {
        "ball": x.str(digits),
        "lower": x.lower().str(digits, radius=False),
        "upper": x.upper().str(digits, radius=False),
    }


def parse_arb_number(text: str) -> arb:
    text = text.strip()
    if "/" in text:
        num, den = text.split("/", 1)
        return arb(int(num.strip())) / int(den.strip())
    return arb(text)


def finite_arb(x: arb) -> bool:
    return bool(x.is_finite()) and not bool(x.is_nan())


def finite_acb_real(z) -> bool:
    return finite_arb(z.real)


def matrix_get(mat: arb_mat, i: int, j: int) -> arb:
    return mat[i, j]


def symmetrize(mat: arb_mat) -> arb_mat:
    n = mat.nrows()
    out = arb_mat(n, n)
    half = arb(1) / 2
    for i in range(n):
        for j in range(n):
            out[i, j] = half * (mat[i, j] + mat[j, i])
    return out


def frobenius_abs_upper(mat: arb_mat) -> arb:
    total = arb(0)
    for i in range(mat.nrows()):
        for j in range(mat.ncols()):
            a = mat[i, j].abs_upper()
            total += a * a
    return total.sqrt().upper()


def point_symmetric_eig_bounds(mat: arb_mat) -> list[tuple[arb, arb]]:
    eigs = mat.eig(multiple=True, algorithm="rump")
    bounds = []
    for eig in eigs:
        if not finite_acb_real(eig):
            raise ValueError("non-finite Arb eigenvalue enclosure")
        bounds.append((eig.real.lower(), eig.real.upper()))
    return bounds


def gershgorin_lmin_lower(block: arb_mat) -> arb:
    """Certified lower bound on lambda_min for a symmetric interval matrix."""
    rows = block.nrows()
    candidates = []
    for i in range(rows):
        radius_sum = arb(0)
        for j in range(rows):
            if i != j:
                radius_sum += block[i, j].abs_upper()
        candidates.append(arb_lower_number(block[i, i]) - radius_sum)
    return arb_lower_number(arb_min(candidates))


def arb_eig_lmin_lower(block: arb_mat, algorithm: str = "rump") -> arb | None:
    """Certified lower bound from Arb eigenvalue enclosures, or None on failure."""
    try:
        eigs = block.eig(multiple=True, algorithm=algorithm)
    except Exception:
        return None
    lows = []
    for eig in eigs:
        if not finite_acb_real(eig):
            return None
        lows.append(arb_lower_number(eig.real))
    if len(lows) != block.nrows():
        return None
    return arb_lower_number(arb_min(lows))


def point_lmin_lower(block: arb_mat) -> tuple[arb, str]:
    eig = arb_eig_lmin_lower(block)
    if eig is not None:
        return eig, "point_arb_eig"
    return gershgorin_lmin_lower(block), "point_gershgorin"


def triple_lmin_lower(block: arb_mat) -> tuple[arb, str, dict[str, dict[str, str] | None]]:
    """Best certified lower bound available for a 3x3 interval block."""
    gersh = gershgorin_lmin_lower(block)
    eig = arb_eig_lmin_lower(block)
    if eig is None:
        return gersh, "gershgorin", {"gershgorin": arb_to_json(gersh), "arb_eig": None}
    best = arb_lower_number(eig.max(gersh))
    source = "arb_eig+gershgorin" if best.overlaps(eig) else "gershgorin"
    return best, source, {"gershgorin": arb_to_json(gersh), "arb_eig": arb_to_json(eig)}


def centered_triple_lmin_lower(block: arb_mat, center_block: arb_mat) -> tuple[arb, str, dict]:
    """Certified Weyl/Frobenius lower bound around a point center block."""
    center_lb, center_source = point_lmin_lower(center_block)
    error = frobenius_abs_upper(block - center_block)
    lower = arb_lower_number(center_lb - error)
    return lower, "centered_weyl_frobenius", {
        "center_lmin": arb_to_json(center_lb),
        "center_source": center_source,
        "error_frobenius_upper": error.str(24),
    }


def point_projector_to_interval(projector: np.ndarray, radius: float) -> arb_mat:
    out = arb_mat(6, 6)
    for i in range(6):
        for j in range(6):
            out[i, j] = aball(float(projector[i, j]), radius)
    return symmetrize(out)


def chart_y_interval(
    pivot: tuple[int, int, int],
    center: tuple[float, ...],
    radii: tuple[float, ...],
) -> tuple[arb_mat, arb_mat]:
    """Interval and midpoint matrices Y for the affine chart with pivot rows = I."""
    if len(center) != 9 or len(radii) != 9:
        raise ValueError("chart boxes use nine Z coordinates")
    pivot_pos = {row: col for col, row in enumerate(pivot)}
    complement = [row for row in range(6) if row not in pivot_pos]

    y = arb_mat(6, 3)
    y_mid = arb_mat(6, 3)
    z = [[aball(center[3 * r + c], radii[3 * r + c]) for c in range(3)] for r in range(3)]
    z_mid = [[aball(center[3 * r + c], 0.0) for c in range(3)] for r in range(3)]
    for row in range(6):
        if row in pivot_pos:
            y[row, pivot_pos[row]] = arb(1)
            y_mid[row, pivot_pos[row]] = arb(1)
        else:
            rr = complement.index(row)
            for col in range(3):
                y[row, col] = z[rr][col]
                y_mid[row, col] = z_mid[rr][col]
    return y, y_mid


def chart_gram_interval(
    pivot: tuple[int, int, int],
    center: tuple[float, ...],
    radii: tuple[float, ...],
) -> tuple[arb_mat, arb_mat, arb_mat, arb_mat]:
    y, y_mid = chart_y_interval(pivot, center, radii)
    return y, y_mid, y.transpose() * y, y_mid.transpose() * y_mid


def inverse_interval_neumann(gram: arb_mat, gram_mid: arb_mat) -> tuple[arb_mat, dict]:
    """Enclose gram^{-1} using a midpoint inverse and a Frobenius Neumann bound."""
    mid_inv = gram_mid.inv()
    delta = gram - gram_mid
    delta_norm = frobenius_abs_upper(delta)
    try:
        eig_bounds = point_symmetric_eig_bounds(gram_mid)
        lambda_min_lower = arb_min([lo for lo, _hi in eig_bounds])
        inv_norm = (arb(1) / lambda_min_lower).upper()
        inv_norm_source = "midpoint_lambda_min"
    except Exception:
        inv_norm = frobenius_abs_upper(mid_inv)
        inv_norm_source = "midpoint_inverse_frobenius"
    q = (inv_norm * delta_norm).upper()
    meta = {
        "method": "neumann",
        "delta_frobenius_upper": delta_norm.str(20),
        "mid_inv_norm_upper": inv_norm.str(20),
        "mid_inv_norm_source": inv_norm_source,
        "q_upper": q.str(20),
    }
    if not bool(q < arb(1)):
        raise ZeroDivisionError(f"Neumann inverse enclosure failed: q={q}")
    inv_radius = (inv_norm * q / (arb(1) - q)).upper()
    out = arb_mat(3, 3)
    for i in range(3):
        for j in range(3):
            out[i, j] = mid_inv[i, j] + radius_ball(inv_radius)
    meta["entry_radius_upper"] = inv_radius.str(20)
    return out, meta


def inverse_interval(gram: arb_mat, gram_mid: arb_mat) -> tuple[arb_mat, dict]:
    try:
        return gram.inv(), {"method": "arb_direct"}
    except Exception as direct_exc:
        inv, meta = inverse_interval_neumann(gram, gram_mid)
        meta["direct_error"] = repr(direct_exc)
        return inv, meta


def chart_projector_interval(
    pivot: tuple[int, int, int],
    center: tuple[float, ...],
    radii: tuple[float, ...],
) -> tuple[arb_mat, dict]:
    """Interval enclosure of P(Z) for the affine chart with pivot rows = I."""
    y, _y_mid, gram, gram_mid = chart_gram_interval(pivot, center, radii)
    inv, meta = inverse_interval(gram, gram_mid)
    return symmetrize(y * inv * y.transpose()), meta


def chart_projector_center(
    pivot: tuple[int, int, int],
    center: tuple[float, ...],
) -> tuple[arb_mat, dict]:
    return chart_projector_interval(pivot, center, tuple([0.0] * 9))


def chart_projector_lipschitz_interval(
    pivot: tuple[int, int, int],
    center: tuple[float, ...],
    radii: tuple[float, ...],
) -> tuple[arb_mat, arb_mat, dict]:
    """Mean-value enclosure of P(Z) from interval derivatives over the box."""
    y, _y_mid, gram, gram_mid = chart_gram_interval(pivot, center, radii)
    inv, meta = inverse_interval(gram, gram_mid)
    p_interval = symmetrize(y * inv * y.transpose())
    p_center, center_meta = chart_projector_center(pivot, center)

    pivot_pos = {row: col for col, row in enumerate(pivot)}
    complement = [row for row in range(6) if row not in pivot_pos]
    q_interval = eye_mat(6) - p_interval
    yt = y.transpose()

    entry_radii = [[arb(0) for _ in range(6)] for _ in range(6)]
    max_deriv_abs = arb(0)
    for rr, row in enumerate(complement):
        for col in range(3):
            coord = 3 * rr + col
            dy = arb_mat(6, 3)
            dy[row, col] = arb(1)
            deriv = symmetrize(q_interval * dy * inv * yt + y * inv * dy.transpose() * q_interval)
            radius = aconst(radii[coord])
            for i in range(6):
                for j in range(6):
                    deriv_abs = deriv[i, j].abs_upper()
                    max_deriv_abs = max_deriv_abs.max(deriv_abs)
                    entry_radii[i][j] += deriv_abs * radius

    out = arb_mat(6, 6)
    max_entry_radius = arb(0)
    for i in range(6):
        for j in range(6):
            radius = entry_radii[i][j].upper()
            max_entry_radius = max_entry_radius.max(radius)
            out[i, j] = p_center[i, j] + radius_ball(radius)

    meta.update({
        "method": "chart_lipschitz",
        "center_method": center_meta,
        "max_derivative_abs_upper": max_deriv_abs.str(20),
        "max_entry_radius_upper": max_entry_radius.str(20),
    })
    return symmetrize(out), p_center, meta


def chart_projector_affine_model(
    pivot: tuple[int, int, int],
    center: tuple[float, ...],
    radii: tuple[float, ...],
) -> tuple[dict, dict]:
    """First-order model P0 + sum h_k D_k plus a certified interval remainder."""
    y, _y_mid, gram, gram_mid = chart_gram_interval(pivot, center, radii)
    inv, meta = inverse_interval(gram, gram_mid)
    p_interval = symmetrize(y * inv * y.transpose())

    y0, _y0_mid, gram0, gram0_mid = chart_gram_interval(pivot, center, tuple([0.0] * 9))
    inv0, center_meta = inverse_interval(gram0, gram0_mid)
    p0 = symmetrize(y0 * inv0 * y0.transpose())

    pivot_pos = {row: col for col, row in enumerate(pivot)}
    complement = [row for row in range(6) if row not in pivot_pos]
    q_interval = eye_mat(6) - p_interval
    q0 = eye_mat(6) - p0
    yt = y.transpose()
    y0t = y0.transpose()

    derivs = []
    remainder_radii = [[arb(0) for _ in range(6)] for _ in range(6)]
    max_deriv_variation = arb(0)
    for rr, row in enumerate(complement):
        for col in range(3):
            coord = 3 * rr + col
            dy = arb_mat(6, 3)
            dy[row, col] = arb(1)
            deriv_interval = symmetrize(q_interval * dy * inv * yt + y * inv * dy.transpose() * q_interval)
            deriv0 = symmetrize(q0 * dy * inv0 * y0t + y0 * inv0 * dy.transpose() * q0)
            derivs.append(deriv0)
            radius = aconst(radii[coord])
            for i in range(6):
                for j in range(6):
                    variation = (deriv_interval[i, j] - deriv0[i, j]).abs_upper()
                    max_deriv_variation = max_deriv_variation.max(variation)
                    remainder_radii[i][j] += variation * radius

    remainder = arb_mat(6, 6)
    max_remainder = arb(0)
    for i in range(6):
        for j in range(6):
            radius = remainder_radii[i][j].upper()
            max_remainder = max_remainder.max(radius)
            remainder[i, j] = radius_ball(radius)

    model = {
        "center": p0,
        "derivatives": derivs,
        "radii": tuple(aconst(r) for r in radii),
        "remainder": symmetrize(remainder),
    }
    meta.update({
        "method": "chart_affine_taylor",
        "center_method": center_meta,
        "max_derivative_variation_abs_upper": max_deriv_variation.str(20),
        "max_remainder_entry_radius_upper": max_remainder.str(20),
    })
    return model, meta


def chart_domain_status(
    pivot: tuple[int, int, int],
    center: tuple[float, ...],
    radii: tuple[float, ...],
) -> dict:
    """Classify det(P_II) >= 1/20 via det(I+Z^T Z) <= 20."""
    _y, _y_mid, gram, gram_mid = chart_gram_interval(pivot, center, radii)
    delta_norm = frobenius_abs_upper(gram - gram_mid)
    try:
        eig_bounds = point_symmetric_eig_bounds(gram_mid)
        lower_factors = []
        upper_factors = []
        for lo, hi in eig_bounds:
            low = lo - delta_norm
            # For this chart Gram = I + Z^T Z, so every true eigenvalue is >= 1.
            if not bool(low > arb(1)):
                low = arb(1)
            lower_factors.append(low.lower())
            upper_factors.append((hi + delta_norm).upper())
        det_lower = arb(1)
        det_upper = arb(1)
        for value in lower_factors:
            det_lower *= value
        for value in upper_factors:
            det_upper *= value
        method = "weyl_midpoint"
    except Exception:
        det = gram.det()
        det_lower = det.lower()
        det_upper = det.upper()
        method = "interval_det_fallback"

    if bool(det_lower > arb(20)):
        status = "outside"
    elif bool(det_upper <= arb(20)):
        status = "inside"
    else:
        status = "unknown"
    return {
        "status": status,
        "method": method,
        "det_gram_lower": arb_to_json(det_lower),
        "det_gram_upper": arb_to_json(det_upper),
        "delta_frobenius_upper": delta_norm.str(20),
    }


def interval_F_lower(
    projector_interval: arb_mat,
    projector_center: arb_mat | None = None,
    bound_mode: str = "direct",
    threshold: arb = SIX,
) -> dict:
    """Certified lower bound for F over an interval projector enclosure."""
    if bound_mode not in {"direct", "centered", "hybrid"}:
        raise ValueError(f"unknown bound mode {bound_mode}")
    if bound_mode in {"centered", "hybrid"} and projector_center is None:
        if bound_mode == "centered":
            raise ValueError("centered bound mode requires projector_center")
        bound_mode = "direct"

    records = []
    lowers = []
    for idx, triple in enumerate(TRIPLES):
        block = arb_mat([[projector_interval[i, j] for j in triple] for i in triple])
        direct_lb = direct_source = direct_detail = None
        centered_lb = centered_source = centered_detail = None
        if bound_mode in {"direct", "hybrid"}:
            direct_lb, direct_source, direct_detail = triple_lmin_lower(block)
        if bound_mode in {"centered", "hybrid"}:
            center_block = arb_mat([[projector_center[i, j] for j in triple] for i in triple])
            centered_lb, centered_source, centered_detail = centered_triple_lmin_lower(block, center_block)

        if bound_mode == "direct":
            lb, source, detail = direct_lb, direct_source, direct_detail
        elif bound_mode == "centered":
            lb, source, detail = centered_lb, centered_source, centered_detail
        else:
            if bool(centered_lb > direct_lb):
                lb, source, detail = centered_lb, centered_source, centered_detail
            else:
                lb, source, detail = direct_lb, direct_source, direct_detail
            detail = {
                "selected": detail,
                "direct": direct_detail,
                "centered": centered_detail,
                "direct_lower": arb_to_json(direct_lb),
                "centered_lower": arb_to_json(centered_lb),
            }
        lowers.append(lb)
        records.append({
            "index": idx,
            "triple": list(triple),
            "lower": arb_to_json(lb),
            "source": source,
            "detail": detail,
        })
    f_lower = arb_lower_number(arb_max(lowers))
    argmax = max(range(len(lowers)), key=lambda i: float(lowers[i].lower()))
    return {
        "bound_mode": bound_mode,
        "F_lower": f_lower,
        "F_lower_json": arb_to_json(f_lower),
        "certifies_GTZ": bool(f_lower.lower() >= SIX),
        "certifies_threshold": bool(f_lower.lower() >= threshold),
        "argmax_lower_index": argmax,
        "argmax_lower_triple": list(TRIPLES[argmax]),
        "triples": records,
    }


def affine_triple_lmin_lower(
    model: dict,
    triple: tuple[int, int, int],
    early_threshold: arb | None = SIX,
) -> tuple[arb, dict]:
    center = model["center"]
    derivatives = model["derivatives"]
    radii = model["radii"]
    remainder = model["remainder"]
    block0 = arb_mat([[center[i, j] for j in triple] for i in triple])
    deriv_blocks = [
        arb_mat([[deriv[i, j] for j in triple] for i in triple])
        for deriv in derivatives
    ]
    rem_block = arb_mat([[remainder[i, j] for j in triple] for i in triple])
    rem_norm = frobenius_abs_upper(rem_block)

    center_lb, center_source = point_lmin_lower(block0)
    linear_radius = arb(0)
    for radius, deriv in zip(radii, deriv_blocks):
        linear_radius += radius * frobenius_abs_upper(deriv)
    cheap_lower = arb_lower_number(center_lb - linear_radius - rem_norm)
    if early_threshold is not None and bool(cheap_lower.lower() >= early_threshold):
        return cheap_lower, {
            "center_lmin": arb_to_json(center_lb),
            "center_source": center_source,
            "linear_frobenius_radius_upper": linear_radius.str(24),
            "remainder_frobenius_upper": rem_norm.str(24),
            "cheap_lower": arb_to_json(cheap_lower),
            "vertices_checked": 0,
            "vertex_min": None,
        }

    vertex_min = None
    vertex_count = 0
    for signs in itertools.product((-1, 1), repeat=len(radii)):
        block = arb_mat([[block0[i, j] for j in range(3)] for i in range(3)])
        for sign, radius, deriv in zip(signs, radii, deriv_blocks):
            step = radius if sign > 0 else -radius
            for i in range(3):
                for j in range(3):
                    block[i, j] += step * deriv[i, j]
        lb, source = point_lmin_lower(block)
        vertex_min = lb if vertex_min is None else vertex_min.min(lb)
        vertex_count += 1

    vertex_lower = arb_lower_number(vertex_min - rem_norm)
    lower = arb_lower_number(vertex_lower.max(cheap_lower))
    return lower, {
        "center_lmin": arb_to_json(center_lb),
        "center_source": center_source,
        "linear_frobenius_radius_upper": linear_radius.str(24),
        "cheap_lower": arb_to_json(cheap_lower),
        "vertex_min": arb_to_json(vertex_min),
        "vertex_lower": arb_to_json(vertex_lower),
        "remainder_frobenius_upper": rem_norm.str(24),
        "vertices_checked": vertex_count,
    }


def affine_model_F_lower(model: dict, early_threshold: arb | None = SIX) -> dict:
    """Certified lower bound using affine triple pencils and interval remainders."""
    center_np = np.array([[float(model["center"][i, j].mid()) for j in range(6)] for i in range(6)])
    center_lams = [
        float(np.linalg.eigvalsh(center_np[np.ix_(triple, triple)])[0])
        for triple in TRIPLES
    ]
    order = sorted(range(len(TRIPLES)), key=lambda i: -center_lams[i])

    records = []
    lowers = []
    certified_early = False
    for idx in order:
        triple = TRIPLES[idx]
        lb, detail = affine_triple_lmin_lower(model, triple, early_threshold)
        lowers.append(lb)
        records.append({
            "index": idx,
            "triple": list(triple),
            "center_lmin_approx": center_lams[idx],
            "lower": arb_to_json(lb),
            "source": "affine_vertex_min_plus_remainder",
            "detail": detail,
        })
        if early_threshold is not None and bool(lb.lower() >= early_threshold):
            certified_early = True
            break

    f_lower = arb_lower_number(arb_max(lowers))
    argmax_record = max(records, key=lambda r: float(arb(r["lower"]["lower"])))
    return {
        "bound_mode": "affine",
        "F_lower": f_lower,
        "F_lower_json": arb_to_json(f_lower),
        "certifies_GTZ": bool(f_lower.lower() >= SIX),
        "certifies_threshold": bool(
            f_lower.lower() >= (early_threshold if early_threshold is not None else SIX)
        ),
        "certified_early": certified_early,
        "triples_evaluated": len(records),
        "argmax_lower_index": argmax_record["index"],
        "argmax_lower_triple": argmax_record["triple"],
        "triples": records,
    }


def point_F(projector: np.ndarray) -> float:
    return float(max(np.linalg.eigvalsh(projector[np.ix_(triple, triple)])[0] for triple in TRIPLES))


def random_projector(rng: np.random.Generator) -> np.ndarray:
    q, _ = np.linalg.qr(rng.standard_normal((6, 3)))
    return q @ q.T


def known_extremals() -> list[tuple[str, np.ndarray]]:
    out = []
    for label, tree, weights in TTSP_CASES:
        out.append((f"TTSP {label}", ttsp_projector(tree, weights)))
    seventh = "verify/data/P514_seventh.npy"
    if os.path.exists(seventh):
        out.append(("OUT-OF-FAMILY (5/14,9/14)", np.load(seventh)))
    return out


def run_known_smoke(radii: list[float], bound_mode: str, threshold: arb) -> list[dict]:
    records = []
    print(f"\nknown-extremal entrywise interval smoke ({bound_mode} bounds):")
    for label, projector in known_extremals():
        f_center = point_F(projector)
        center = point_projector_to_interval(projector, 0.0)
        per_radius = []
        for radius in radii:
            rec = interval_F_lower(point_projector_to_interval(projector, radius), center, bound_mode, threshold)
            per_radius.append({
                "radius": radius,
                "F_center": f_center,
                "F_lower": rec["F_lower_json"],
                "certifies_GTZ": rec["certifies_GTZ"],
                "certifies_threshold": rec["certifies_threshold"],
                "argmax_lower_triple": rec["argmax_lower_triple"],
            })
            print(
                f"  {label:<42} rad={radius:.1e} "
                f"F_center={f_center:.16f} "
                f"F_lb={float(rec['F_lower'].lower()):.16f} "
                f"cert_threshold={rec['certifies_threshold']} "
                f"cert_gtz={rec['certifies_GTZ']}"
            )
        records.append({"label": label, "records": per_radius})
    return records


def run_random_smoke(count: int, radius: float, seed: int, bound_mode: str, threshold: arb) -> list[dict]:
    rng = np.random.default_rng(seed)
    records = []
    print(f"\nrandom-projector entrywise interval smoke ({bound_mode} bounds):")
    for i in range(count):
        projector = random_projector(rng)
        f_center = point_F(projector)
        center = point_projector_to_interval(projector, 0.0)
        rec = interval_F_lower(point_projector_to_interval(projector, radius), center, bound_mode, threshold)
        row = {
            "sample": i,
            "radius": radius,
            "F_center": f_center,
            "F_lower": rec["F_lower_json"],
            "certifies_GTZ": rec["certifies_GTZ"],
            "certifies_threshold": rec["certifies_threshold"],
            "argmax_lower_triple": rec["argmax_lower_triple"],
        }
        records.append(row)
        print(
            f"  sample={i:02d} rad={radius:.1e} "
            f"F_center={f_center:.16f} "
            f"F_lb={float(rec['F_lower'].lower()):.16f} "
            f"cert_threshold={rec['certifies_threshold']} "
            f"cert_gtz={rec['certifies_GTZ']}"
        )
    return records


@dataclass(order=True)
class QueueBox:
    priority: float
    serial: int
    pivot: tuple[int, int, int] = field(compare=False)
    center: tuple[float, ...] = field(compare=False)
    radii: tuple[float, ...] = field(compare=False)
    depth: int = field(compare=False)


def box_to_json(box: QueueBox) -> dict:
    return {
        "priority": box.priority,
        "serial": box.serial,
        "pivot": list(box.pivot),
        "center": list(box.center),
        "radii": list(box.radii),
        "depth": box.depth,
    }


def box_from_json(record: dict, serial: int, priority_mode: str) -> QueueBox:
    pivot = tuple(int(x) for x in record["pivot"])
    center = tuple(float(x) for x in record["center"])
    radii = tuple(float(x) for x in record["radii"])
    depth = int(record["depth"])
    return QueueBox(box_priority(pivot, center, radii, priority_mode, depth), serial, pivot, center, radii, depth)


def box_priority(
    pivot: tuple[int, int, int],
    center: tuple[float, ...],
    radii: tuple[float, ...],
    mode: str,
    depth: int = 0,
) -> float:
    if mode == "deep":
        return -depth
    if mode == "depth":
        return -max(radii)
    center_f = point_F(center_projector_in_chart(pivot, center))
    if mode == "center-high":
        return -center_f
    if mode == "center-low":
        return center_f
    raise ValueError(f"unknown priority mode {mode}")


def center_projector_in_chart(pivot: tuple[int, int, int], center: tuple[float, ...]) -> np.ndarray:
    pivot_pos = {row: col for col, row in enumerate(pivot)}
    complement = [row for row in range(6) if row not in pivot_pos]
    y = np.zeros((6, 3))
    for row in range(6):
        if row in pivot_pos:
            y[row, pivot_pos[row]] = 1.0
        else:
            rr = complement.index(row)
            y[row, :] = center[3 * rr:3 * rr + 3]
    return y @ np.linalg.inv(y.T @ y) @ y.T


def split_box(box: QueueBox, serial_start: int, priority_mode: str) -> tuple[QueueBox, QueueBox]:
    dim = max(range(9), key=lambda i: box.radii[i])
    half = box.radii[dim] / 2.0
    centers = []
    for sign in (-1.0, 1.0):
        c = list(box.center)
        c[dim] += sign * half
        centers.append(tuple(c))
    radii = list(box.radii)
    radii[dim] = half
    children = []
    for k, center in enumerate(centers):
        depth = box.depth + 1
        priority = box_priority(box.pivot, center, tuple(radii), priority_mode, depth)
        children.append(QueueBox(priority, serial_start + k, box.pivot, center, tuple(radii), depth))
    return children[0], children[1]


def evaluate_chart_box(
    box: QueueBox,
    domain_filter: bool,
    bound_mode: str,
    projector_mode: str,
    threshold: arb,
) -> dict:
    max_radius = max(box.radii)
    inverse_meta = None
    domain = {"status": "unchecked", "det_gram": None}
    if domain_filter:
        try:
            domain = chart_domain_status(box.pivot, box.center, box.radii)
        except Exception as exc:
            domain = {"status": "unknown", "det_gram": None, "error": repr(exc)}
        if domain["status"] == "outside":
            p_center = center_projector_in_chart(box.pivot, box.center)
            center_f = point_F(p_center)
            return {
                "pivot": list(box.pivot),
                "depth": box.depth,
                "max_radius": max_radius,
                "center_F": center_f,
                "domain": domain,
                "inverse_ok": False,
                "inverse_meta": None,
                "error": None,
                "F_lower": None,
                "F_lower_json": None,
                "certifies_GTZ": False,
                "certifies_threshold": False,
                "argmax_lower_triple": None,
            }
    try:
        if projector_mode == "cascade":
            pint, interval_meta = chart_projector_interval(box.pivot, box.center, box.radii)
            pcenter = None
            if bound_mode in {"centered", "hybrid"}:
                pcenter, center_meta = chart_projector_center(box.pivot, box.center)
                interval_meta["center_method"] = center_meta
            interval_bound = interval_F_lower(pint, pcenter, bound_mode, threshold)
            if interval_bound["certifies_threshold"]:
                bound = interval_bound
                inverse_meta = {
                    "method": "cascade",
                    "selected": "interval",
                    "interval": interval_meta,
                    "affine_evaluated": False,
                }
            else:
                model, affine_meta = chart_projector_affine_model(box.pivot, box.center, box.radii)
                affine_bound = affine_model_F_lower(model, threshold)
                if bool(affine_bound["F_lower"] > interval_bound["F_lower"]):
                    bound = affine_bound
                    selected = "affine"
                else:
                    bound = interval_bound
                    selected = "interval"
                inverse_meta = {
                    "method": "cascade",
                    "selected": selected,
                    "interval": interval_meta,
                    "affine": affine_meta,
                    "affine_evaluated": True,
                    "interval_F_lower": interval_bound["F_lower_json"],
                    "affine_F_lower": affine_bound["F_lower_json"],
                }
        elif projector_mode == "affine":
            model, inverse_meta = chart_projector_affine_model(box.pivot, box.center, box.radii)
            bound = affine_model_F_lower(model, threshold)
        elif projector_mode == "both":
            pint, interval_meta = chart_projector_interval(box.pivot, box.center, box.radii)
            pcenter = None
            if bound_mode in {"centered", "hybrid"}:
                pcenter, center_meta = chart_projector_center(box.pivot, box.center)
                interval_meta["center_method"] = center_meta
            interval_bound = interval_F_lower(pint, pcenter, bound_mode, threshold)

            lpint, lpcenter, lipschitz_meta = chart_projector_lipschitz_interval(box.pivot, box.center, box.radii)
            lipschitz_bound = interval_F_lower(lpint, lpcenter, bound_mode, threshold)

            if bool(lipschitz_bound["F_lower"] > interval_bound["F_lower"]):
                bound = lipschitz_bound
                selected = "lipschitz"
            else:
                bound = interval_bound
                selected = "interval"
            inverse_meta = {
                "method": "both",
                "selected": selected,
                "interval": interval_meta,
                "lipschitz": lipschitz_meta,
                "interval_F_lower": interval_bound["F_lower_json"],
                "lipschitz_F_lower": lipschitz_bound["F_lower_json"],
            }
        elif projector_mode == "interval":
            pint, inverse_meta = chart_projector_interval(box.pivot, box.center, box.radii)
            pcenter = None
            if bound_mode in {"centered", "hybrid"}:
                pcenter, center_meta = chart_projector_center(box.pivot, box.center)
                inverse_meta["center_method"] = center_meta
        elif projector_mode == "lipschitz":
            pint, pcenter, inverse_meta = chart_projector_lipschitz_interval(box.pivot, box.center, box.radii)
        else:
            raise ValueError(f"unknown projector mode {projector_mode}")

        if projector_mode not in {"both", "affine", "cascade"} and bound_mode in {"centered", "hybrid"} and pcenter is None:
            pcenter, center_meta = chart_projector_center(box.pivot, box.center)
            inverse_meta["center_method"] = center_meta
        if projector_mode not in {"both", "affine", "cascade"}:
            bound = interval_F_lower(pint, pcenter, bound_mode, threshold)
        inverse_ok = True
        error = None
    except Exception as exc:
        bound = None
        inverse_ok = False
        error = repr(exc)
    p_center = center_projector_in_chart(box.pivot, box.center)
    center_f = point_F(p_center)
    out = {
        "pivot": list(box.pivot),
        "depth": box.depth,
        "max_radius": max_radius,
        "center_F": center_f,
        "domain": domain,
        "inverse_ok": inverse_ok,
        "inverse_meta": inverse_meta,
        "error": error,
    }
    if bound is not None:
        affine_vertex_checks = 0
        if bound.get("bound_mode") == "affine":
            affine_vertex_checks = sum(
                rec.get("detail", {}).get("vertices_checked", 0)
                for rec in bound.get("triples", [])
            )
        out.update({
            "F_lower": bound["F_lower"],
            "F_lower_json": bound["F_lower_json"],
            "certifies_GTZ": bound["certifies_GTZ"],
            "certifies_threshold": bound["certifies_threshold"],
            "argmax_lower_triple": bound["argmax_lower_triple"],
            "bound_summary": {
                "bound_mode": bound.get("bound_mode"),
                "certified_early": bound.get("certified_early"),
                "triples_evaluated": bound.get("triples_evaluated"),
                "affine_vertex_checks": affine_vertex_checks,
            },
        })
    else:
        out.update({
            "F_lower": None,
            "F_lower_json": None,
            "certifies_GTZ": False,
            "certifies_threshold": False,
            "argmax_lower_triple": None,
            "bound_summary": None,
        })
    return out


def run_chart_bnb(args: argparse.Namespace) -> dict:
    queue = []
    serial = 0
    if args.frontier_input:
        with open(args.frontier_input) as fh:
            frontier_doc = json.load(fh)
        frontier = frontier_doc.get("chart_bnb", frontier_doc).get("frontier")
        if frontier is None:
            raise ValueError(f"no chart_bnb.frontier in {args.frontier_input}")
        for record in frontier:
            heapq.heappush(queue, box_from_json(record, serial, args.priority))
            serial += 1
        pivots = sorted({tuple(box.pivot) for box in queue})
    else:
        pivots = list(itertools.combinations(range(6), 3))
        if args.charts is not None:
            pivots = pivots[:args.charts]
        zero = tuple([0.0] * 9)
        root_radii = tuple([SQRT19] * 9)
        for pivot in pivots:
            priority = box_priority(pivot, zero, root_radii, args.priority, 0)
            heapq.heappush(queue, QueueBox(priority, serial, pivot, zero, root_radii, 0))
            serial += 1

    stats = {
        "charts": len(pivots),
        "frontier_input": args.frontier_input,
        "max_boxes": args.max_boxes,
        "min_radius": args.min_radius,
        "priority": args.priority,
        "domain_filter": args.domain_filter,
        "bound_mode": args.bound_mode,
        "projector_mode": args.projector_mode,
        "threshold": args.threshold,
        "processed": 0,
        "certified": 0,
        "certified_gtz": 0,
        "split": 0,
        "unresolved": 0,
        "inverse_failures": 0,
        "domain_outside": 0,
        "domain_unknown": 0,
        "affine_evaluations": 0,
        "affine_vertex_checks": 0,
        "best_lower": None,
        "best_center_F": None,
        "smallest_unresolved_radius": None,
        "unresolved_examples": [],
    }

    print("\nchart branch-and-bound smoke:")
    print(
        f"  charts={len(pivots)} max_boxes={args.max_boxes} "
        f"min_radius={args.min_radius:.2e} priority={args.priority} "
        f"domain_filter={args.domain_filter} bound_mode={args.bound_mode} "
        f"projector_mode={args.projector_mode} threshold={args.threshold}"
    )
    if args.frontier_input:
        print(f"  loaded frontier boxes={len(queue)} from {args.frontier_input}")

    while queue and stats["processed"] < args.max_boxes:
        box = heapq.heappop(queue)
        rec = evaluate_chart_box(
            box,
            args.domain_filter,
            args.bound_mode,
            args.projector_mode,
            args.threshold_arb,
        )
        stats["processed"] += 1
        if rec["domain"]["status"] == "outside":
            stats["domain_outside"] += 1
            continue
        if rec["domain"]["status"] == "unknown":
            stats["domain_unknown"] += 1
        if not rec["inverse_ok"]:
            stats["inverse_failures"] += 1
        meta = rec.get("inverse_meta") or {}
        if meta.get("method") == "chart_affine_taylor":
            stats["affine_evaluations"] += 1
        elif meta.get("method") == "cascade" and meta.get("affine_evaluated"):
            stats["affine_evaluations"] += 1
        if bound := rec.get("bound_summary"):
            stats["affine_vertex_checks"] += bound.get("affine_vertex_checks", 0)
        if stats["best_center_F"] is None or rec["center_F"] > stats["best_center_F"]:
            stats["best_center_F"] = rec["center_F"]
        if rec["F_lower"] is not None:
            lb = float(rec["F_lower"].lower())
            if stats["best_lower"] is None or lb > stats["best_lower"]:
                stats["best_lower"] = lb

        if rec["certifies_threshold"]:
            stats["certified"] += 1
            if rec["certifies_GTZ"]:
                stats["certified_gtz"] += 1
        elif rec["max_radius"] <= args.min_radius:
            stats["unresolved"] += 1
            r = rec["max_radius"]
            if stats["smallest_unresolved_radius"] is None or r < stats["smallest_unresolved_radius"]:
                stats["smallest_unresolved_radius"] = r
            if len(stats["unresolved_examples"]) < 10:
                lite = dict(rec)
                lite.pop("F_lower", None)
                stats["unresolved_examples"].append(lite)
        else:
            left, right = split_box(box, serial, args.priority)
            serial += 2
            heapq.heappush(queue, left)
            heapq.heappush(queue, right)
            stats["split"] += 1

        if stats["processed"] % args.progress_every == 0 or stats["processed"] == 1:
            print(
                f"  [{stats['processed']:5d}] certified={stats['certified']:5d} "
                f"gtz={stats['certified_gtz']:5d} "
                f"split={stats['split']:5d} unresolved={stats['unresolved']:5d} "
                f"outside={stats['domain_outside']:5d} queue={len(queue):5d} "
                f"best_lb={stats['best_lower']}"
            )

    stats["queue_remaining"] = len(queue)
    if args.save_frontier:
        stats["frontier"] = [box_to_json(box) for box in sorted(queue)]
    print(
        f"  done: processed={stats['processed']} certified={stats['certified']} "
        f"gtz={stats['certified_gtz']} "
        f"split={stats['split']} unresolved={stats['unresolved']} "
        f"outside={stats['domain_outside']} inverse_failures={stats['inverse_failures']} "
        f"affine_evals={stats['affine_evaluations']} queue={len(queue)}"
    )
    print(f"  best certified lower bound seen: {stats['best_lower']}")
    print(f"  best center F seen:              {stats['best_center_F']}")
    return stats


def parse_radii(text: str) -> list[float]:
    return [float(part) for part in text.split(",") if part.strip()]


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--prec", type=int, default=160, help="Arb precision in bits")
    parser.add_argument("--known-smoke", action="store_true")
    parser.add_argument("--radii", default="0,1e-16,1e-14,1e-12,1e-10,1e-8")
    parser.add_argument("--random-smoke", type=int, default=0)
    parser.add_argument("--random-radius", type=float, default=1e-8)
    parser.add_argument("--bound-mode", choices=["direct", "centered", "hybrid"], default="hybrid")
    parser.add_argument("--threshold", default="1/6", help="certification target, e.g. 1/6 or 0.15")
    parser.add_argument("--seed", type=int, default=20260804)
    parser.add_argument("--chart-bnb", action="store_true")
    parser.add_argument("--charts", type=int, default=None, help="prefix of the 20 affine charts to test")
    parser.add_argument("--frontier-input", default=None, help="JSON output containing chart_bnb.frontier")
    parser.add_argument("--save-frontier", action="store_true")
    parser.add_argument("--max-boxes", type=int, default=200)
    parser.add_argument("--min-radius", type=float, default=0.5)
    parser.add_argument(
        "--projector-mode",
        choices=["interval", "lipschitz", "both", "affine", "cascade"],
        default="interval",
    )
    parser.add_argument("--priority", choices=["center-high", "center-low", "depth", "deep"], default="center-high")
    parser.add_argument("--domain-filter", action="store_true")
    parser.add_argument("--progress-every", type=int, default=25)
    parser.add_argument("--output", default="verify/out/v19_interval_bounds.json")
    args = parser.parse_args()

    ctx.prec = args.prec
    args.threshold_arb = parse_arb_number(args.threshold)
    os.makedirs("verify/out", exist_ok=True)

    print("=" * 78)
    print("INTERVAL LOWER BOUNDS FOR GTZ(6,3)")
    print(f"Arb precision: {ctx.prec} bits")
    print(f"Certification threshold: {args.threshold}")
    print("=" * 78)

    out = {
        "precision_bits": args.prec,
        "threshold": args.threshold,
        "known_smoke": None,
        "random_smoke": None,
        "chart_bnb": None,
    }

    if args.known_smoke:
        out["known_smoke"] = run_known_smoke(
            parse_radii(args.radii),
            args.bound_mode,
            args.threshold_arb,
        )
    if args.random_smoke:
        out["random_smoke"] = run_random_smoke(
            args.random_smoke,
            args.random_radius,
            args.seed,
            args.bound_mode,
            args.threshold_arb,
        )
    if args.chart_bnb:
        out["chart_bnb"] = run_chart_bnb(args)

    path = args.output
    with open(path, "w") as fh:
        json.dump(out, fh, indent=1)
    print(f"\nwrote {path}")
    print("=" * 78)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
