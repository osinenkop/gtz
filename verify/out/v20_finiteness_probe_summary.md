# v20 finiteness-probe summary

Date: 2026-07-31.

Status: **numerical diagnostic only**.  This does not prove finiteness.

## Goal

The finiteness obstruction is a non-sharp equality extremal.  A
positive-dimensional equality family would force some point on the level set
`F=1/6` to have a nontrivial critical cone.

`verify/v20_finiteness_probe.py` classifies equality-mode KKT hits from `v18` by:

- reconstructing the projector from `best_y`;
- recomputing the actual active set, not just the selected active set;
- canonicalizing the actual active set under `S_6`;
- comparing with the seven known certified extremal orbits;
- computing numerical rank and positive-multiplier margin for the actual active
  gradients.

## Inputs

- `verify/out/v18_kkt_equality_1_4.json`
- `verify/out/v18_kkt_equality_5_6.json`
- `verify/out/v18_kkt_equality_7_8.json`
- `verify/out/v18_kkt_equality_9_9.json`
- `verify/out/v18_kkt_equality_10_13.json`
- `verify/out/v18_kkt_equality_14_20.json`

Together these cover all nonempty active-set orbit sizes `1..20` in equality
mode.

## Result

- input records: `2135`
- low-residual equality hits: `3`
- distinct actual active-set orbits: `2`
- new actual active-set orbits: `0`
- numerical non-sharp candidates: `0`
- actual active size `<=9`: `0`

The three equality hits are all selected-subset manifestations of known sharp
orbits:

| actual active size | rank | positive margin | known orbit |
|---:|---:|---:|---|
| 12 | 9 | `5.000e-02` | TTSP 12-active orbit |
| 12 | 9 | `5.000e-02` | TTSP 12-active orbit |
| 13 | 9 | `2.222e-02` | out-of-family `(5/14,9/14)` orbit |

## Interpretation

This is stronger evidence for finiteness than another random descent batch,
because it searches directly for the kind of equality point that would obstruct
finiteness.  The result is negative: no new actual active orbit, no low-active
equality point, and no numerical non-sharp signature.

The proof target suggested by this is:

> Every equality KKT point has actual active set in one of the known sharp
> active-set orbits, or at least has rank-9 positive spanning active gradients.

Proving that would imply the extremal set is finite by the compact-discrete
argument already recorded in `proofs/finiteness-of-extremal-set.md`.
