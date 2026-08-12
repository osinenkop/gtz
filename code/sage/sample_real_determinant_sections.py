#!/usr/bin/env python3
"""Numerically sample real sections of a saturated determinant active locus.

This is a probe, not a proof.  It solves

    det(P(Z)_TT - I/6) = 0  for T in A

plus deterministic random affine linear sections through a known point, then
classifies the real roots by the actual GTZ inequalities:

    active blocks:   lambda_min(P_TT) >= 1/6
    inactive blocks: lambda_min(P_TT) <= 1/6.

The intended use is to see whether the semialgebraic inequalities collapse the
positive-dimensional determinant locus before attempting a certified proof.
"""
from __future__ import annotations

import argparse
import itertools
import json
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

ROOT = Path(__file__).resolve().parents[2]
VERIFY = ROOT / "verify"
if str(VERIFY) not in sys.path:
    sys.path.insert(0, str(VERIFY))
from v17_active_orbits import TRIPLES, TTSP_CASES, active_mask, mask_triples, ttsp_projector  # noqa: E402
from v20_finiteness_probe import gradient_matrix, positive_multiplier_margin  # noqa: E402

TH = 1.0 / 6.0


def known_projector(label_substring: str) -> tuple[str, np.ndarray]:
    matches = [
        (label, tree, weights)
        for label, tree, weights in TTSP_CASES
        if label_substring.lower() in label.lower()
    ]
    if len(matches) != 1:
        labels = [label for label, *_ in TTSP_CASES]
        raise SystemExit(f"expected one TTSP label match for {label_substring!r}; labels={labels}")
    label, tree, weights = matches[0]
    return label, ttsp_projector(tree, weights)


def retract(matrix: np.ndarray) -> np.ndarray:
    u, _, vt = np.linalg.svd(matrix, full_matrices=False)
    return u @ vt


def indices_from_mask(mask: int) -> tuple[int, ...]:
    return tuple(idx for idx in range(len(TRIPLES)) if mask & (1 << idx))


def parse_active_indices(text: str) -> tuple[int, ...]:
    return tuple(int(part.strip()) for part in text.split(",") if part.strip())


def load_v30_projector(path: str) -> tuple[str, np.ndarray, tuple[int, ...]]:
    data = json.loads(Path(path).read_text())
    label = f"v30 size {data.get('size')} canon {data.get('canon')}"
    projector = retract(np.array(data["best_point"], dtype=float).reshape(6, 3))
    return label, projector @ projector.T, tuple(int(idx) for idx in data["active"])


def load_base_projector(args) -> tuple[str, np.ndarray]:
    if args.v30_input:
        label, projector, _active = load_v30_projector(args.v30_input)
        return label, projector
    if args.projector_npy:
        path = Path(args.projector_npy)
        return path.stem, np.load(path)
    return known_projector(args.known_label)


def resolve_active_indices(args, projector: np.ndarray) -> tuple[int, ...]:
    sources = [
        bool(args.active_indices),
        bool(args.active_mask),
        bool(args.v30_input),
    ]
    if sum(sources) > 1:
        raise SystemExit("use at most one of --active-indices, --active-mask, --v30-input")
    if args.active_indices:
        return parse_active_indices(args.active_indices)
    if args.active_mask:
        return indices_from_mask(int(args.active_mask, 0))
    if args.v30_input:
        _label, _projector, active = load_v30_projector(args.v30_input)
        return active
    return active_indices_from_projector(projector)


def active_indices_from_projector(projector: np.ndarray) -> tuple[int, ...]:
    mask, _values = active_mask(projector)
    triples = mask_triples(mask)
    index = {triple: idx for idx, triple in enumerate(TRIPLES)}
    return tuple(index[triple] for triple in triples)


def parse_chart_rows(text: str, projector: np.ndarray) -> tuple[int, int, int]:
    if text == "auto":
        return auto_chart_rows(projector)
    rows = tuple(int(x.strip()) for x in text.split(",") if x.strip())
    if len(rows) != 3 or len(set(rows)) != 3 or any(row < 0 or row >= 6 for row in rows):
        raise SystemExit("--chart-rows must be 'auto' or three distinct row indices, e.g. 0,1,2")
    return rows


def auto_chart_rows(projector: np.ndarray) -> tuple[int, int, int]:
    eigvals, eigvecs = np.linalg.eigh(projector)
    basis = eigvecs[:, np.argsort(-eigvals)[:3]]
    best_rows = None
    best_det = -1.0
    for rows in itertools.combinations(range(6), 3):
        det = abs(np.linalg.det(basis[list(rows), :]))
        if det > best_det:
            best_det = det
            best_rows = rows
    if best_rows is None or best_det < 1e-10:
        raise ValueError(f"no nonsingular 3-row chart found; best determinant={best_det}")
    return tuple(best_rows)


def chart_z_from_projector(projector: np.ndarray, chart_rows: tuple[int, int, int]) -> np.ndarray:
    eigvals, eigvecs = np.linalg.eigh(projector)
    basis = eigvecs[:, np.argsort(-eigvals)[:3]]
    top = basis[list(chart_rows), :]
    det_top = np.linalg.det(top)
    if abs(det_top) < 1e-10:
        raise ValueError(f"chart rows {chart_rows} nearly singular: det={det_top}")
    chart_basis = basis @ np.linalg.inv(top)
    complement = tuple(row for row in range(6) if row not in chart_rows)
    return chart_basis[list(complement), :].reshape(9)


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


def make_sections(z0: np.ndarray, count: int, seed: int) -> tuple[np.ndarray, np.ndarray]:
    rng = np.random.default_rng(seed)
    a = rng.standard_normal((count, z0.size))
    a /= np.linalg.norm(a, axis=1, keepdims=True)
    b = a @ z0
    return a, b


def make_coordinate_sections(z0: np.ndarray, indices: tuple[int, ...]) -> tuple[np.ndarray, np.ndarray]:
    a = np.zeros((len(indices), z0.size))
    for row, idx in enumerate(indices):
        if idx < 0 or idx >= z0.size:
            raise SystemExit(f"coordinate section index out of range: {idx}")
        a[row, idx] = 1.0
    b = z0[list(indices)]
    return a, b


def make_equal_sum_sections(z0: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Three simple sections through the all-equal base chart point.

    These are chosen to be interpretable rather than generic:

      z00 = z01,  z00 = z10,  sum(z) = sum(z0).
    """
    a = np.zeros((3, z0.size))
    a[0, 0] = 1.0
    a[0, 1] = -1.0
    a[1, 0] = 1.0
    a[1, 3] = -1.0
    a[2, :] = 1.0
    b = a @ z0
    return a, b


def make_zero_sum_integer_sections(z0: np.ndarray, count: int, seed: int, bound: int = 5) -> tuple[np.ndarray, np.ndarray]:
    """Random integer sections with coefficient sum zero and exact RHS zero at z0.

    This is tailored to the base chart where z0 is a scalar multiple of the
    all-ones vector.  The resulting equations have rational coefficients:

        a_0 z_0 + ... + a_8 z_8 = 0,    sum_i a_i = 0.
    """
    rng = np.random.default_rng(seed)
    rows = []
    seen = set()
    while len(rows) < count:
        coeff = rng.integers(-bound, bound + 1, size=z0.size)
        coeff[-1] = -int(np.sum(coeff[:-1]))
        if np.all(coeff == 0):
            continue
        gcd = int(np.gcd.reduce(np.abs(coeff)))
        if gcd > 1:
            coeff = coeff // gcd
        key = tuple(int(x) for x in coeff)
        if key in seen:
            continue
        seen.add(key)
        rows.append(coeff.astype(float))
    a = np.vstack(rows)
    a /= np.linalg.norm(a, axis=1, keepdims=True)
    b = np.zeros(count)
    return a, b


def residual(
    z: np.ndarray,
    active: tuple[int, ...],
    a: np.ndarray,
    b: np.ndarray,
    section_weight: float,
    chart_rows: tuple[int, int, int],
) -> np.ndarray:
    p = projector_from_z(z, chart_rows)
    dets = determinant_values(p, active)
    sections = a @ z - b
    return np.concatenate([dets, section_weight * sections])


def classify_root(
    z: np.ndarray,
    active: tuple[int, ...],
    a: np.ndarray,
    b: np.ndarray,
    chart_rows: tuple[int, int, int],
) -> dict:
    p = projector_from_z(z, chart_rows)
    lam = lambdas(p)
    active_set = set(active)
    inactive = [idx for idx in range(len(TRIPLES)) if idx not in active_set]
    active_min_shift = min(float(lam[idx] - TH) for idx in active)
    inactive_max_shift = max(float(lam[idx] - TH) for idx in inactive) if inactive else -float("inf")
    f_minus = float(np.max(lam) - TH)
    actual_active = [idx for idx, value in enumerate(lam) if abs(float(value - TH)) <= 1e-7]
    lmat, gaps = gradient_matrix(p, actual_active)
    mult = positive_multiplier_margin(lmat)
    return {
        "z": [float(x) for x in z],
        "section_residual_max": float(np.max(np.abs(a @ z - b))) if len(a) else 0.0,
        "det_residual_max": float(np.max(np.abs(determinant_values(p, active)))) if active else 0.0,
        "F_minus_1_6": f_minus,
        "active_min_lambda_minus_1_6": active_min_shift,
        "inactive_max_lambda_minus_1_6": inactive_max_shift,
        "passes_active_psd": bool(active_min_shift >= -1e-7),
        "passes_inactive": bool(inactive_max_shift <= 1e-7),
        "passes_equality_inequalities": bool(active_min_shift >= -1e-7 and inactive_max_shift <= 1e-7),
        "actual_active": actual_active,
        "actual_active_triples": [list(TRIPLES[idx]) for idx in actual_active],
        "actual_active_size": len(actual_active),
        "selected_equals_actual": bool(set(actual_active) == active_set),
        "gap_min": float(min(gaps)) if gaps else None,
        "rank_L": int(np.linalg.matrix_rank(lmat, tol=1e-8)) if lmat.size else 0,
        "positive_multiplier": mult,
    }


def root_seen(roots: list[np.ndarray], z: np.ndarray, tol: float) -> bool:
    return any(np.linalg.norm(z - old) <= tol * max(1.0, np.linalg.norm(old)) for old in roots)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--known-label", default="P(S(e,e,e),e,e,e)")
    parser.add_argument("--projector-npy", help="load a 6x6 projector instead of a TTSP known label")
    parser.add_argument("--v30-input", help="load projector and active set from verify/v30_margin_*.json")
    parser.add_argument("--active-indices", default="", help="comma-separated active triple indices")
    parser.add_argument("--active-mask", default="", help="active-set bitmask")
    parser.add_argument("--sections", type=int, default=3)
    parser.add_argument(
        "--section-mode",
        choices=["random", "coordinate", "equal-sum", "zero-sum-integer"],
        default="random",
    )
    parser.add_argument(
        "--coordinate-indices",
        default="",
        help="comma-separated z-coordinate indices for --section-mode coordinate",
    )
    parser.add_argument("--chart-rows", default="auto")
    parser.add_argument("--seed", type=int, default=101)
    parser.add_argument("--starts", type=int, default=600)
    parser.add_argument("--section-weight", type=float, default=1.0)
    parser.add_argument("--root-tol", type=float, default=1e-7)
    parser.add_argument("--accept-residual", type=float, default=1e-8)
    parser.add_argument("--max-nfev", type=int, default=6000)
    parser.add_argument("--out", default="code/sage/out/real_section_known_base.json")
    args = parser.parse_args()

    label, p0 = load_base_projector(args)
    active = resolve_active_indices(args, p0)
    chart_rows = parse_chart_rows(args.chart_rows, p0)
    chart_complement = tuple(row for row in range(6) if row not in chart_rows)
    z0 = chart_z_from_projector(p0, chart_rows)
    if args.section_mode == "random":
        a, b = make_sections(z0, args.sections, args.seed)
    elif args.section_mode == "coordinate":
        indices = tuple(int(x.strip()) for x in args.coordinate_indices.split(",") if x.strip())
        if len(indices) != args.sections:
            raise SystemExit("--coordinate-indices length must match --sections")
        a, b = make_coordinate_sections(z0, indices)
    else:
        if args.section_mode == "zero-sum-integer":
            a, b = make_zero_sum_integer_sections(z0, args.sections, args.seed)
        else:
            if args.sections != 3:
                raise SystemExit("--section-mode equal-sum requires --sections 3")
            a, b = make_equal_sum_sections(z0)
    if args.section_mode == "equal-sum":
        if args.sections != 3:
            raise SystemExit("--section-mode equal-sum requires --sections 3")

    rng = np.random.default_rng(args.seed + 1)
    starts = [z0]
    scales = [0.02, 0.08, 0.25, 0.8, 2.0, 5.0]
    while len(starts) < args.starts:
        scale = scales[(len(starts) - 1) % len(scales)]
        starts.append(z0 + scale * rng.standard_normal(9))

    roots: list[np.ndarray] = []
    attempts = []
    for start_idx, start in enumerate(starts):
        sol = least_squares(
            residual,
            start,
            args=(active, a, b, args.section_weight, chart_rows),
            method="trf",
            x_scale="jac",
            ftol=1e-12,
            xtol=1e-12,
            gtol=1e-12,
            max_nfev=args.max_nfev,
        )
        res = residual(sol.x, active, a, b, args.section_weight, chart_rows)
        res_norm = float(np.linalg.norm(res))
        res_max = float(np.max(np.abs(res)))
        accepted = bool(res_max <= args.accept_residual)
        duplicate = bool(accepted and root_seen(roots, sol.x, args.root_tol))
        if accepted and not duplicate:
            roots.append(sol.x.copy())
        attempts.append({
            "start": start_idx,
            "success": bool(sol.success),
            "cost": float(sol.cost),
            "residual_norm": res_norm,
            "residual_max": res_max,
            "accepted": accepted,
            "duplicate": duplicate,
            "nfev": int(sol.nfev),
        })

    records = [classify_root(root, active, a, b, chart_rows) for root in roots]
    records.sort(key=lambda row: (
        not row["passes_equality_inequalities"],
        row["F_minus_1_6"],
        row["inactive_max_lambda_minus_1_6"],
    ))

    out = {
        "known_label": label,
        "active_indices": list(active),
        "active_triples": [list(TRIPLES[idx]) for idx in active],
        "z0": [float(x) for x in z0],
        "chart_rows": list(chart_rows),
        "chart_complement": list(chart_complement),
        "sections": args.sections,
        "section_mode": args.section_mode,
        "section_seed": args.seed,
        "section_matrix": a.tolist(),
        "section_rhs": b.tolist(),
        "starts": args.starts,
        "accepted_roots": len(records),
        "passes_equality_inequalities": sum(1 for row in records if row["passes_equality_inequalities"]),
        "attempt_summary": {
            "accepted_attempts": sum(1 for row in attempts if row["accepted"]),
            "duplicate_attempts": sum(1 for row in attempts if row["duplicate"]),
            "best_residual_max": min(row["residual_max"] for row in attempts) if attempts else None,
        },
        "roots": records,
        "attempts": attempts,
    }
    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(out, indent=1) + "\n")

    print("=" * 78)
    print("REAL DETERMINANT SECTION SAMPLE")
    print("=" * 78)
    print(f"known label:       {label}")
    print(f"active size:       {len(active)}")
    print(f"chart rows:        {chart_rows}")
    print(f"sections:          {args.sections}")
    print(f"section mode:      {args.section_mode}")
    print(f"starts:            {args.starts}")
    print(f"accepted roots:    {len(records)}")
    print(f"passes inequalities: {out['passes_equality_inequalities']}")
    print(f"best residual max: {out['attempt_summary']['best_residual_max']:.3e}")
    for idx, row in enumerate(records[:20]):
        print(
            f"[{idx:02d}] pass={row['passes_equality_inequalities']} "
            f"F-1/6={row['F_minus_1_6']:+.3e} "
            f"act_min={row['active_min_lambda_minus_1_6']:+.3e} "
            f"inact_max={row['inactive_max_lambda_minus_1_6']:+.3e} "
            f"|Aactual|={row['actual_active_size']} rankL={row['rank_L']}"
        )
    print(f"wrote {out_path}")
    print("=" * 78)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
