# v21 semialgebraic-route summary

Date: 2026-07-31.

Status: **formulation and sizing only**.  This is not a proof of GTZ(6,3) or of
finiteness.

## Goal

The finiteness obstruction is a non-sharp equality point of
`F(P)=max_T lambda_min(P_TT)`.  A positive-dimensional equality family would
force such a point.

`verify/v21_semialgebraic_route.py` writes the simple-active obstruction as a
semialgebraic system in a Grassmann chart:

- 9 chart variables `z`;
- 9 tangent-witness variables `h`;
- no per-active-triple eigenvector variables, using cofactor patches instead.

## Polynomial model

In the chart `Y=[I;Z]`,

```text
P(Z) = Y (Y^T Y)^(-1) Y^T = N(Z) / d(Z),
d = det(Y^T Y) > 0.
```

For an active triple `T`, set

```text
M_T(Z) = 6 N_TT(Z) - d(Z) I_3.
```

The simple-active equality condition is `M_T >= 0`, `det M_T = 0`, and
`rank M_T = 2`.  On a cofactor patch, a kernel vector is obtained as a cross
product of two rows of `M_T`.  The non-sharp tangent witness is encoded by
cleared polynomial inequalities

```text
w_T^T [sum_j h_j (d dN/dz_j - N dd/dz_j)_TT] w_T <= 0.
```

## Fast audit result

Command:

```bash
.venv/bin/python -u verify/v21_semialgebraic_route.py
```

Output file:

```text
verify/out/v21_semialgebraic_route.json
```

Degree audit:

| quantity | degree / upper bound |
|---|---:|
| `d` | 6 |
| entries of `N` | 6 |
| raw active determinant | <= 18 |
| saturated active determinant | 6 |
| active PSD 2x2 minors | <= 12 |
| cofactor-patch norm | <= 24 |
| cleared directional inequality | <= 35 |

One cofactor-patch branch has:

| active size | variables | equalities | inequalities, inactive relaxed |
|---:|---:|---:|---:|
| 8 | 18 | 9 | 65 |
| 9 | 18 | 10 | 73 |
| 10 | 18 | 11 | 81 |
| 12 | 18 | 13 | 97 |
| 13 | 18 | 14 | 105 |

## Interpretation

The key point is that the simple-active obstruction has fixed ambient dimension
18, independent of active-set size.  The hard part is branch count:
`3^|A|` cofactor patches per active set before symmetry/pruning.

For simple-active points, sharpness requires positive spanning in a
9-dimensional tangent space, so at least 10 active gradients are needed.
Therefore a simple-active equality point with actual `|A| <= 9` is automatically
non-sharp.

This count does not cover nonsimple active blocks.  Those require a separate
subgradient/SDP-style semialgebraic branch.

## Next proof target

Show that the semialgebraic non-sharp obstruction is empty outside the known
sharp active-set orbits, with nonsimple active eigenvalues handled separately.

For an infeasibility proof, inactive constraints may be relaxed away: if the
relaxed obstruction is empty, the true obstruction is empty.  If the relaxed
system has solutions, inactive inequalities must be restored to classify them.
