# v19 interval-bound route summary

Date: 2026-07-31.

Status: **partial infrastructure only**.  This is not a proof of GTZ(6,3).

## What was added

- `verify/v19_interval_bounds.py`
- Arb interval lower bounds for `lambda_min(P_TT)` on all 20 triples.
- Entrywise interval smoke tests around known extremals and random projectors.
- Affine-chart interval enclosure for
  `P(Z)=Y(Y^T Y)^{-1}Y^T`, with `Y=[I;Z]`.
- Centered Weyl/Frobenius spectral lower bounds:
  `lambda_min(M) >= lambda_min(M0) - ||M-M0||_F`.
- A chart-Lipschitz projector enclosure using the exact differential
  `dP=(I-P)dY(Y^TY)^{-1}Y^T + Y(Y^TY)^{-1}dY^T(I-P)`.
- A first-order affine/Taylor projector model
  `P(Z0+h)=P0+sum h_i D_i+R(h)`, with `R` bounded by interval derivative
  variation over the box.  Triple lower bounds minimize the affine matrix pencil
  over all `2^9=512` vertices and subtract the certified nonlinear remainder.
- Frontier save/resume support via `--save-frontier` and `--frontier-input`,
  so cheap global passes can be followed by expensive affine refinement only on
  surviving boxes.
- Branch-and-bound smoke driver over the 20 affine charts.
- Optional chart-domain filter using the Cauchy-Binet cover condition
  `det(P_II) >= 1/20`, equivalently `det(I+Z^T Z) <= 20`.

## Smoke observations

Known equality extremals:

- At radius `0`, some known extremals certify numerically and some miss because
  this smoke test uses floating projector centers rather than the exact
  `v9`/`v11` algebraic projectors.
- At any positive entrywise radius tested (`1e-16` and above), the entrywise
  interval lower bound drops below `1/6`.
- Interpretation: boxes containing equality extremals cannot be cleared by
  naive entrywise spectral bounds alone.  They must be covered by exact local
  sharpness/radius certificates.

Random projectors:

- Five random projector boxes with entrywise radius `1e-8` all certified
  `F(P) > 1/6`, with lower bounds between about `0.298` and `0.485`.
- Interpretation: the Arb triple-bound primitive works well away from the sharp
  boundary.

Unfiltered chart branch-and-bound:

- Command:
  `.venv/bin/python -u verify/v19_interval_bounds.py --chart-bnb --charts 20 --max-boxes 60000 --min-radius 0.125 --priority center-high`
- Processed boxes: `60000`
- Certified boxes: `29060`
- Split boxes: `30940`
- Inverse failures: `167`
- Queue remaining: `1900`
- Best certified lower bound seen: `0.28516037323345245`

This confirms that the interval machinery certifies many small high-margin chart
boxes, but does not approach a complete cover at this subdivision scale.

Matched `10000`-box chart calibration, same priority and radius:

| projector enclosure | spectral bound | certified | split | queue | inverse failures |
|---|---|---:|---:|---:|---:|
| interval | direct | 4692 | 5308 | 636 | 167 |
| interval | hybrid centered/direct | 4738 | 5262 | 544 | 167 |
| Lipschitz derivative | hybrid centered/direct | 4658 | 5342 | 704 | 167 |
| best of interval and Lipschitz | hybrid centered/direct | 4658 | 5342 | 704 | 167 |

The centered spectral bound gives a small improvement over direct interval
eigenvalue bounds.  The derivative/Lipschitz projector enclosure gives stronger
lower bounds on some high-margin boxes, but does not improve the overall
certification count in this calibration.

Affine/Taylor calibration:

| run | certified | split | queue | inverse failures | notes |
|---|---:|---:|---:|---:|---|
| one chart, 300 boxes, direct | 112 | 188 | 77 | 34 | baseline |
| one chart, 300 boxes, affine | 117 | 183 | 67 | 34 | stronger bounds, slower |
| one chart, 300 boxes, cascade with affine precheck | 117 | 183 | 67 | 34 | `266` affine evaluations, `522752` vertex eigen checks |
| all charts, 500 boxes, cascade | 152 | 348 | 216 | 167 | slight improvement over the earlier direct 500-box count `142` |

This is the first version that genuinely keeps first-order chart correlations.
It improves the frontier modestly and produces much stronger local lower bounds,
but the improvement is not large enough to make exhaustive B&B look imminent.

Two-stage frontier run:

| stage | input | processed | certified | split | queue | affine evals | vertex checks |
|---|---|---:|---:|---:|---:|---:|---:|
| 1 | roots, cheap interval/hybrid | 5000 | 2283 | 2717 | 454 | 0 | 0 |
| 2 | stage-1 frontier, cascade | 300 | 271 | 29 | 212 | 186 | 92160 |
| 3 | stage-2 frontier, cascade | 500 | 252 | 248 | 208 | 500 | 985600 |
| deep diagnostic | stage-3 frontier, cascade/deep | 1000 | 505 | 495 | 198 | 1000 | 2002944 |

Stage 2 is the clearest success: applying affine checks to the cheap frontier
certified most of the first 300 boxes.  Later refinement stabilized near a
roughly 200-box frontier where each pass certifies about half the boxes but
splits the other half, so the queue size stops decreasing quickly.

The surviving frontier is not simply an equality-neighborhood frontier.  It
still contains broad chart boxes with center `F` around `0.31`, well above
`1/6`, and large coordinate radius.  These look more like chart/inverse
conditioning artifacts than likely counterexamples.

Chart-domain-filtered runs:

- The determinant-domain filter discards poor-pivot boxes as intended.
- In low-`F` priority mode on one chart, the filter discarded `1458` boxes in a
  `5000`-box run.
- Remaining hard boxes tend to have determinant lower bounds very close to `20`
  and large upper bounds, so they sit near the overlap boundary between charts.

## Interpretation

The route is viable as infrastructure, but not yet as a proof.  The current
blockers are:

1. The affine chart interval enclosure is too conservative near low-`F` regions.
2. The chart-domain boundary `det(I+Z^T Z)=20` needs a sharper covering
   strategy, otherwise many boxes remain undecided.
3. Equality neighborhoods need to be removed from the global search by exact
   local certificates; entrywise interval spectral bounds lose the sharp
   equality immediately.
4. First-order Lipschitz enclosure alone is not sharp enough.
5. The affine/Taylor model helps, but its vertex enumeration is expensive and
   only modestly improves certification counts at the tested scale.
6. The current chart cover leaves broad low-priority boxes that are difficult
   for the interval inverse, even when their center value is far above `1/6`.

## Relaxed-threshold calibration

`verify/v19_interval_bounds.py` now accepts `--threshold`, so the same machinery
can target fixed relaxed bounds `F >= c` with `c < 1/6`.

A calibration at `c = 0.15 = 1/6 - 1/60`:

```bash
.venv/bin/python -u verify/v19_interval_bounds.py --chart-bnb --charts 20 \
  --max-boxes 1000 --min-radius 0.25 --priority center-high \
  --projector-mode cascade --bound-mode hybrid --threshold 0.15 \
  --output verify/out/v19_relaxed_015_cascade_1000.json
```

processed 1000 boxes, certified 396 at the relaxed threshold, certified 368 of
those at the sharp `1/6` threshold, and left a queue of 228 boxes.  Thus the
relaxed target helps only modestly with the current cover.  The bottleneck is
still interval overestimation/chart conditioning, not only the sharp equality
level.

## Next technical move

The most realistic next step is a two-stage global run:

1. Use cheap interval/hybrid bounds to clear high-margin boxes.
2. Apply the affine/Taylor vertex bound only to the remaining frontier.
3. Remove neighborhoods of the seven known equality points using exact
   local sharpness/radius certificates rather than global interval boxes.

If this still leaves a large frontier, the next upgrade should be an optimized
affine bound: vectorized approximate vertex screening followed by rigorous Arb
verification only where the approximate margin is small.

The other high-value upgrade is a better chart cover or chart-domain split.  The
current affine chart cube `[ -sqrt(19), sqrt(19) ]^9` is too blunt; many hard
boxes appear to be coordinate artifacts of this cover rather than genuinely
near the sharp set.
