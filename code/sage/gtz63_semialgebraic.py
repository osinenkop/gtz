#!/usr/bin/env python3
"""Sage helpers for GTZ(6,3) semialgebraic active-set systems.

This module is meant to be run with the Sage conda environment:

    ~/miniforge3/bin/mamba run -n sage python code/sage/export_semialgebraic_system.py ...

It deliberately uses ordinary Python plus ``sage.all`` rather than ``.sage``
preprocessor syntax.
"""
from __future__ import annotations

import itertools
import json
from dataclasses import dataclass
from pathlib import Path

from sage.all import Matrix, PolynomialRing, QQ


TRIPLES = list(itertools.combinations(range(6), 3))
TRIPLE_INDEX = {triple: i for i, triple in enumerate(TRIPLES)}
ROW_PAIR_MAP = {"01": (0, 1), "02": (0, 2), "12": (1, 2)}


@dataclass(frozen=True)
class Chart:
    ring: object
    z: tuple
    h: tuple
    kernels: tuple
    slack: tuple
    inv: tuple
    y: object
    d: object
    n: object


@dataclass(frozen=True)
class PolynomialSystem:
    mode: str
    ring: object
    active_indices: tuple[int, ...]
    row_pairs: tuple[tuple[int, int], ...]
    equalities: tuple
    nonzero: tuple
    inequalities: tuple
    metadata: dict

    @property
    def variables(self) -> tuple[str, ...]:
        return tuple(str(name) for name in self.ring.variable_names())


@dataclass(frozen=True)
class SingularTextSystem:
    mode: str
    variables: tuple[str, ...]
    active_indices: tuple[int, ...]
    row_pairs: tuple[tuple[int, int], ...]
    equalities: tuple[str, ...]
    nonzero: tuple[str, ...]
    equality_degrees: tuple[int, ...]
    nonzero_degrees: tuple[int, ...]
    metadata: dict


def parse_indices(text: str) -> tuple[int, ...]:
    if not text.strip():
        return ()
    out = tuple(int(x.strip()) for x in text.split(",") if x.strip())
    bad = [i for i in out if i < 0 or i >= len(TRIPLES)]
    if bad:
        raise ValueError(f"triple indices out of range: {bad}")
    return out


def indices_from_mask(mask: int) -> tuple[int, ...]:
    return tuple(i for i in range(len(TRIPLES)) if mask & (1 << i))


def mask_from_indices(indices: tuple[int, ...]) -> int:
    mask = 0
    for idx in indices:
        mask |= 1 << idx
    return mask


def parse_row_pairs(text: str, n_active: int) -> tuple[tuple[int, int], ...]:
    raw = [x.strip() for x in text.split(",") if x.strip()]
    if not raw:
        raw = ["01"]
    for item in raw:
        if item not in ROW_PAIR_MAP:
            raise ValueError(f"unknown row pair {item!r}; use one of 01, 02, 12")
    if len(raw) == 1:
        raw *= n_active
    if len(raw) != n_active:
        raise ValueError(f"got {len(raw)} row pairs for {n_active} active triples")
    return tuple(ROW_PAIR_MAP[item] for item in raw)


def known_active_indices(label_substring: str, path: str | Path = "verify/out/v17_active_orbits.json") -> tuple[int, ...]:
    data = json.loads(Path(path).read_text())
    matches = [
        record for record in data["known_extremals"]
        if label_substring.lower() in record["label"].lower()
    ]
    if not matches:
        labels = [record["label"] for record in data["known_extremals"]]
        raise ValueError(f"no known extremal matching {label_substring!r}; labels={labels}")
    if len(matches) > 1:
        labels = [record["label"] for record in matches]
        raise ValueError(f"ambiguous known extremal {label_substring!r}; matches={labels}")
    return tuple(TRIPLE_INDEX[tuple(triple)] for triple in matches[0]["active_triples"])


def make_ring(
    include_h: bool,
    kernel_count: int = 0,
    slack_count: int = 0,
    inverse_count: int = 0,
    order: str = "degrevlex",
):
    names = [f"z{i}" for i in range(9)]
    if include_h:
        names.extend(f"h{i}" for i in range(9))
    for i in range(kernel_count):
        names.extend(f"w{i}{j}" for j in range(3))
    names.extend(f"s{i}" for i in range(slack_count))
    names.extend(f"u{i}" for i in range(inverse_count))
    return PolynomialRing(QQ, names, order=order)


def standard_chart(
    ring,
    include_h: bool,
    kernel_count: int = 0,
    slack_count: int = 0,
    inverse_count: int = 0,
) -> Chart:
    gens = dict(zip(ring.variable_names(), ring.gens()))
    z = tuple(gens[f"z{i}"] for i in range(9))
    h = tuple(gens[f"h{i}"] for i in range(9)) if include_h else ()
    kernels = tuple(
        tuple(gens[f"w{i}{j}"] for j in range(3))
        for i in range(kernel_count)
    )
    slack = tuple(gens[f"s{i}"] for i in range(slack_count))
    inv = tuple(gens[f"u{i}"] for i in range(inverse_count))

    y = Matrix(ring, 6, 3, 0)
    for i in range(3):
        y[i, i] = ring.one()
    for r in range(3):
        for c in range(3):
            y[3 + r, c] = z[3 * r + c]

    gram = y.transpose() * y
    d = gram.det()
    n = y * gram.adjugate() * y.transpose()
    return Chart(ring=ring, z=z, h=h, kernels=kernels, slack=slack, inv=inv, y=y, d=d, n=n)


def active_block(chart: Chart, triple_index: int):
    triple = TRIPLES[triple_index]
    return Matrix(
        chart.ring,
        3,
        3,
        lambda a, b: 6 * chart.n[triple[a], triple[b]]
        - (chart.d if a == b else chart.ring.zero()),
    )


def active_determinant(chart: Chart, triple_index: int):
    """Saturated active determinant equation on the full-rank chart."""
    raw = active_block(chart, triple_index).det()
    quotient, remainder = raw.quo_rem(chart.d ** 2)
    if remainder != chart.ring.zero():
        raise ArithmeticError("active determinant is not divisible by d^2")
    return quotient


def cross_row_kernel(mat, pair: tuple[int, int]):
    r1 = [mat[pair[0], j] for j in range(3)]
    r2 = [mat[pair[1], j] for j in range(3)]
    return (
        r1[1] * r2[2] - r1[2] * r2[1],
        r1[2] * r2[0] - r1[0] * r2[2],
        r1[0] * r2[1] - r1[1] * r2[0],
    )


def patch_norm(kernel: tuple) -> object:
    return sum(entry * entry for entry in kernel)


def directional_block(chart: Chart, triple_index: int):
    if not chart.h:
        raise ValueError("directional block requires h variables")
    triple = TRIPLES[triple_index]
    block = Matrix(chart.ring, 3, 3, 0)
    for j, zj in enumerate(chart.z):
        dd = chart.d.derivative(zj)
        for a in range(3):
            for b in range(3):
                nij = chart.n[triple[a], triple[b]]
                block[a, b] += chart.h[j] * (chart.d * nij.derivative(zj) - nij * dd)
    return block


def quadratic_form(vec: tuple, mat) -> object:
    total = mat.base_ring().zero()
    for i in range(3):
        for j in range(3):
            total += vec[i] * mat[i, j] * vec[j]
    return total


def _poly_text(poly) -> str:
    return f"({poly})"


def _sum_text(parts: list[str]) -> str:
    if not parts:
        return "0"
    return " + ".join(f"({part})" for part in parts)


def patch_norm_text(kernel: tuple) -> str:
    return _sum_text([f"{_poly_text(entry)}*{_poly_text(entry)}" for entry in kernel])


def det3_text(mat) -> str:
    return (
        f"{_poly_text(mat[0, 0])}*("
        f"{_poly_text(mat[1, 1])}*{_poly_text(mat[2, 2])}"
        f" - {_poly_text(mat[1, 2])}*{_poly_text(mat[2, 1])})"
        f" - {_poly_text(mat[0, 1])}*("
        f"{_poly_text(mat[1, 0])}*{_poly_text(mat[2, 2])}"
        f" - {_poly_text(mat[1, 2])}*{_poly_text(mat[2, 0])})"
        f" + {_poly_text(mat[0, 2])}*("
        f"{_poly_text(mat[1, 0])}*{_poly_text(mat[2, 1])}"
        f" - {_poly_text(mat[1, 1])}*{_poly_text(mat[2, 0])})"
    )


def directional_q_text(chart: Chart, triple_index: int, kernel: tuple) -> str:
    """Factored text for w^T D_T(h) w, avoiding polynomial expansion in Sage."""
    if not chart.h:
        raise ValueError("directional text requires h variables")
    triple = TRIPLES[triple_index]
    h_terms = []
    for j, zj in enumerate(chart.z):
        dd = chart.d.derivative(zj)
        block_terms = []
        for a in range(3):
            for b in range(3):
                nij = chart.n[triple[a], triple[b]]
                dn = nij.derivative(zj)
                middle = f"{_poly_text(chart.d)}*{_poly_text(dn)} - {_poly_text(nij)}*{_poly_text(dd)}"
                block_terms.append(f"{_poly_text(kernel[a])}*({middle})*{_poly_text(kernel[b])}")
        h_terms.append(f"{chart.h[j]}*({_sum_text(block_terms)})")
    return _sum_text(h_terms)


def determinant_system(
    active_indices: tuple[int, ...],
    order: str = "degrevlex",
    invert_d: bool = False,
) -> PolynomialSystem:
    ring = make_ring(include_h=False, inverse_count=1 if invert_d else 0, order=order)
    chart = standard_chart(ring, include_h=False, inverse_count=1 if invert_d else 0)
    equalities = list(active_determinant(chart, idx) for idx in active_indices)
    if invert_d:
        equalities.append(chart.inv[0] * chart.d - 1)
    return PolynomialSystem(
        mode="det_saturated_inverted" if invert_d else "det_saturated",
        ring=ring,
        active_indices=active_indices,
        row_pairs=(),
        equalities=tuple(equalities),
        nonzero=(),
        inequalities=(),
        metadata={
            "chart": "Y=[I;Z]",
            "description": "saturated active determinant equations det(6*N_TT-d*I)/d^2=0",
            "mask": mask_from_indices(active_indices),
            "saturated_by": "d^2, where d=det(Y^T Y)",
            "invert_d": invert_d,
        },
    )


def determinant_text_system(
    active_indices: tuple[int, ...],
    order: str = "degrevlex",
    invert_d: bool = False,
) -> SingularTextSystem:
    ring = make_ring(include_h=False, inverse_count=1 if invert_d else 0, order=order)
    chart = standard_chart(ring, include_h=False, inverse_count=1 if invert_d else 0)
    equalities = list(active_determinant(chart, idx) for idx in active_indices)
    if invert_d:
        equalities.append(chart.inv[0] * chart.d - 1)
    return SingularTextSystem(
        mode="det_saturated_inverted" if invert_d else "det_saturated",
        variables=tuple(str(name) for name in ring.variable_names()),
        active_indices=active_indices,
        row_pairs=(),
        equalities=tuple(str(poly) for poly in equalities),
        nonzero=(),
        equality_degrees=tuple(int(poly.total_degree()) for poly in equalities),
        nonzero_degrees=(),
        metadata={
            "chart": "Y=[I;Z]",
            "description": "saturated active determinant equations det(6*N_TT-d*I)/d^2=0",
            "mask": mask_from_indices(active_indices),
            "ring_order": order,
            "saturated_by": "d^2, where d=det(Y^T Y)",
            "invert_d": invert_d,
        },
    )


def nonsharp_text_system(
    active_indices: tuple[int, ...],
    row_pairs: tuple[tuple[int, int], ...],
    order: str = "degrevlex",
    patch_inverses: bool = False,
) -> SingularTextSystem:
    """Factored nonsharp system for Singular export.

    This avoids expanding the degree-35 directional polynomials in Sage.  It is
    the default export path for nonsharp systems; use ``nonsharp_system`` only
    when an actual Sage polynomial ideal is needed.
    """
    inv_count = len(active_indices) if patch_inverses else 0
    ring = make_ring(
        include_h=True,
        slack_count=len(active_indices),
        inverse_count=inv_count,
        order=order,
    )
    chart = standard_chart(
        ring,
        include_h=True,
        slack_count=len(active_indices),
        inverse_count=inv_count,
    )
    equalities = [_sum_text([f"{h}^2" for h in chart.h]) + " - 1"]
    nonzero = []
    equality_degrees = [2]
    nonzero_degrees = []

    for local_idx, (triple_idx, pair) in enumerate(zip(active_indices, row_pairs)):
        m = active_block(chart, triple_idx)
        w = cross_row_kernel(m, pair)
        q_text = directional_q_text(chart, triple_idx, w)
        wnorm_text = patch_norm_text(w)

        det_poly = active_determinant(chart, triple_idx)
        equalities.append(str(det_poly))
        equality_degrees.append(int(det_poly.total_degree()))
        equalities.append(f"({q_text}) + {chart.slack[local_idx]}^2")
        equality_degrees.append(35)
        if patch_inverses:
            equalities.append(f"{chart.inv[local_idx]}*({wnorm_text}) - 1")
            equality_degrees.append(25)
        else:
            nonzero.append(wnorm_text)
            nonzero_degrees.append(24)

    return SingularTextSystem(
        mode="nonsharp_factored",
        variables=tuple(str(name) for name in ring.variable_names()),
        active_indices=active_indices,
        row_pairs=row_pairs,
        equalities=tuple(equalities),
        nonzero=tuple(nonzero),
        equality_degrees=tuple(equality_degrees),
        nonzero_degrees=tuple(nonzero_degrees),
        metadata={
            "chart": "Y=[I;Z]",
            "description": "factored simple-active cofactor-patch nonsharp witness system",
            "mask": mask_from_indices(active_indices),
            "patch_inverses": patch_inverses,
            "ring_order": order,
            "note": "directional inequalities are intentionally not expanded by Sage",
        },
    )


def nonsharp_system(
    active_indices: tuple[int, ...],
    row_pairs: tuple[tuple[int, int], ...],
    order: str = "degrevlex",
    patch_inverses: bool = False,
) -> PolynomialSystem:
    inv_count = len(active_indices) if patch_inverses else 0
    ring = make_ring(
        include_h=True,
        slack_count=len(active_indices),
        inverse_count=inv_count,
        order=order,
    )
    chart = standard_chart(
        ring,
        include_h=True,
        slack_count=len(active_indices),
        inverse_count=inv_count,
    )
    equalities = []
    nonzero = []
    inequalities = []

    equalities.append(sum(x * x for x in chart.h) - 1)

    for local_idx, (triple_idx, pair) in enumerate(zip(active_indices, row_pairs)):
        m = active_block(chart, triple_idx)
        w = cross_row_kernel(m, pair)
        wnorm = patch_norm(w)
        q = quadratic_form(w, directional_block(chart, triple_idx))
        equalities.append(active_determinant(chart, triple_idx))
        equalities.append(q + chart.slack[local_idx] ** 2)
        if patch_inverses:
            equalities.append(chart.inv[local_idx] * wnorm - 1)
        else:
            nonzero.append(wnorm)
        inequalities.append(("patch_norm_gt_0", wnorm))
        inequalities.append(("directional_le_0_encoded_by_slack", q))

    return PolynomialSystem(
        mode="nonsharp",
        ring=ring,
        active_indices=active_indices,
        row_pairs=row_pairs,
        equalities=tuple(equalities),
        nonzero=tuple(nonzero),
        inequalities=tuple(inequalities),
        metadata={
            "chart": "Y=[I;Z]",
            "description": "simple-active cofactor-patch nonsharp witness system",
            "mask": mask_from_indices(active_indices),
            "patch_inverses": patch_inverses,
        },
    )


def kernel_nonsharp_system(
    active_indices: tuple[int, ...],
    order: str = "degrevlex",
) -> PolynomialSystem:
    """Lower-degree simple-active nonsharp system with explicit kernel vectors.

    For each active triple T this uses variables w_T and equations
    M_T w_T = 0, ||w_T||^2 = 1.  This avoids the degree blow-up and branch count
    of cofactor patches.  It is an algebraic relaxation until PSD/simplicity and
    inactive inequalities are handled separately.
    """
    ring = make_ring(
        include_h=True,
        kernel_count=len(active_indices),
        slack_count=len(active_indices),
        order=order,
    )
    chart = standard_chart(
        ring,
        include_h=True,
        kernel_count=len(active_indices),
        slack_count=len(active_indices),
    )
    equalities = [sum(x * x for x in chart.h) - 1]
    for local_idx, triple_idx in enumerate(active_indices):
        m = active_block(chart, triple_idx)
        w = chart.kernels[local_idx]
        for row in range(3):
            equalities.append(sum(m[row, col] * w[col] for col in range(3)))
        equalities.append(sum(entry * entry for entry in w) - 1)
        q = quadratic_form(w, directional_block(chart, triple_idx))
        equalities.append(q + chart.slack[local_idx] ** 2)

    return PolynomialSystem(
        mode="kernel_nonsharp",
        ring=ring,
        active_indices=active_indices,
        row_pairs=(),
        equalities=tuple(equalities),
        nonzero=(),
        inequalities=(),
        metadata={
            "chart": "Y=[I;Z]",
            "description": "explicit-kernel simple-active nonsharp witness system",
            "mask": mask_from_indices(active_indices),
            "patch_inverses": False,
            "note": (
                "lower-degree algebraic relaxation; PSD, simplicity, and inactive "
                "inequalities are external semialgebraic constraints"
            ),
        },
    )


def system_summary(system: PolynomialSystem) -> dict:
    degrees = [int(poly.total_degree()) for poly in system.equalities]
    nonzero_degrees = [int(poly.total_degree()) for poly in system.nonzero]
    return {
        "mode": system.mode,
        "variables": list(system.variables),
        "n_variables": len(system.variables),
        "active_indices": list(system.active_indices),
        "active_triples": [list(TRIPLES[i]) for i in system.active_indices],
        "active_size": len(system.active_indices),
        "row_pairs": [list(pair) for pair in system.row_pairs],
        "n_equalities": len(system.equalities),
        "equality_degrees": degrees,
        "max_equality_degree": max(degrees) if degrees else None,
        "n_nonzero_constraints": len(system.nonzero),
        "nonzero_degrees": nonzero_degrees,
        "metadata": system.metadata,
    }


def text_system_summary(system: SingularTextSystem) -> dict:
    return {
        "mode": system.mode,
        "variables": list(system.variables),
        "n_variables": len(system.variables),
        "active_indices": list(system.active_indices),
        "active_triples": [list(TRIPLES[i]) for i in system.active_indices],
        "active_size": len(system.active_indices),
        "row_pairs": [list(pair) for pair in system.row_pairs],
        "n_equalities": len(system.equalities),
        "equality_degrees": list(system.equality_degrees),
        "max_equality_degree": max(system.equality_degrees) if system.equality_degrees else None,
        "n_nonzero_constraints": len(system.nonzero),
        "nonzero_degrees": list(system.nonzero_degrees),
        "metadata": system.metadata,
    }


def singular_order(order: str) -> str:
    if order in {"degrevlex", "degrevlex(1)"}:
        return "dp"
    if order in {"lex", "lexicographic"}:
        return "lp"
    return "dp"


def write_singular_text_script(
    system: SingularTextSystem,
    path: str | Path,
    order: str = "degrevlex",
    compute: bool = False,
) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    varlist = ",".join(system.variables)
    lines = [
        "// Generated by code/sage/export_semialgebraic_system.py",
        f"// mode: {system.mode}",
        f"// active_indices: {list(system.active_indices)}",
        "// Directional expressions are factored intentionally.",
        f"ring r = 0,({varlist}),({singular_order(order)});",
        "",
    ]
    for i, expr in enumerate(system.equalities):
        lines.append(f"poly f{i} = {expr};")
    if system.equalities:
        lines.append(f"ideal I = {','.join(f'f{i}' for i in range(len(system.equalities)))};")
    else:
        lines.append("ideal I = 0;")
    lines.extend([
        "",
        '"nvars";',
        f"{len(system.variables)};",
        '"neqs";',
        f"{len(system.equalities)};",
    ])
    if system.nonzero:
        for i, expr in enumerate(system.nonzero):
            lines.append(f"poly nz{i} = {expr};")
        lines.append(f"ideal NZ = {','.join(f'nz{i}' for i in range(len(system.nonzero)))};")
    if compute:
        lines.extend([
            "option(redSB);",
            "ideal G = groebner(I);",
            '"groebner_size";',
            "size(G);",
            '"dimension";',
            "dim(G);",
        ])
    path.write_text("\n".join(lines) + "\n")


def write_singular_script(system: PolynomialSystem, path: str | Path, compute: bool = False) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    varlist = ",".join(system.variables)
    lines = [
        "// Generated by code/sage/export_semialgebraic_system.py",
        f"// mode: {system.mode}",
        f"// active_indices: {list(system.active_indices)}",
        f"ring r = 0,({varlist}),({singular_order(system.ring.term_order().name())});",
        "",
    ]
    for i, poly in enumerate(system.equalities):
        lines.append(f"poly f{i} = {poly};")
    if system.equalities:
        lines.append(f"ideal I = {','.join(f'f{i}' for i in range(len(system.equalities)))};")
    else:
        lines.append("ideal I = 0;")
    lines.extend([
        "",
        '"nvars";',
        f"{len(system.variables)};",
        '"neqs";',
        f"{len(system.equalities)};",
    ])
    if system.nonzero:
        for i, poly in enumerate(system.nonzero):
            lines.append(f"poly nz{i} = {poly};")
        lines.append(f"ideal NZ = {','.join(f'nz{i}' for i in range(len(system.nonzero)))};")
    if compute:
        lines.extend([
            'option(redSB);',
            'ideal G = groebner(I);',
            '"groebner_size";',
            'size(G);',
            '"dimension";',
            'dim(G);',
        ])
    path.write_text("\n".join(lines) + "\n")
