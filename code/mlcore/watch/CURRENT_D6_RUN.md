# Current D6/S6 MLCore Run

Updated: `2026-08-11T14:20Z`.

## Local processes

- `1596911`: `d6_fullsave` MLCore watcher, still alive with `--hours 9.25`.
  Its watched jobs have all finished; it is only polling until the deadline.
- `1604860`: `d6_support` MLCore watcher, still alive with `--hours 8`.
  Its watched jobs have all finished; it is only polling until the deadline.
- `1683035`: `d6_followup` local reconstruction loop, still alive with
  `--hours 9.25`.
- `1692177`: `d6_preempt_s8` MLCore watcher, still alive with `--hours 9`.
  Its five preemptible size-8 support-6 `_allcand` clone jobs have all
  succeeded; it is only polling until the deadline.

The scoped cleanup helper is:

```bash
bash code/mlcore/watch/kill_gtz_watchers.sh --dry-run
bash code/mlcore/watch/kill_gtz_watchers.sh
bash code/mlcore/watch/kill_gtz_watchers.sh --dry-run --dir d6_fullsave
bash code/mlcore/watch/kill_gtz_watchers.sh --dry-run --dir d6_followup
bash code/mlcore/watch/kill_gtz_watchers.sh --dry-run --dir d6_preempt_s8
```

Use `--followup-only` to stop only the local reconstruction loop. Use `--dir`
with a run directory name to target one watcher/follow-up directory.

## Watched outputs

- Full-candidate watcher:
  `code/mlcore/watch/d6_fullsave/latest_summary.md`
- Support watcher:
  `code/mlcore/watch/d6_support/latest_summary.md`
- Preemptible size-8 watcher:
  `code/mlcore/watch/d6_preempt_s8/latest_summary.md`
- Reconstruction follow-up:
  `code/mlcore/watch/d6_followup/latest_summary.md`

## Current artifact status

- Size-7 support-19 `_allcand` has all six primes available:
  `32003,32009,32027,32029,32051,32057`.
- Alternate size-7 support-0 `_allcand` has three primes available:
  `32003,32009,32027`.
- Size-7 support-8 limited has three primes available:
  `32003,32009,32027`.
- Size-8 support-6 `_allcand` has now completed over all six primes:
  `32003,32009,32027,32029,32051,32057`.
- For size-8 support-6 `_allcand`, every prime has the same modular profile:
  matrix rank `79327`, generator rank `79071`, support defect `63`, kernel
  dimension `1072`, and `262/262` passing candidates.

## Current conclusions

- Six-prime size-7 support-19 `_allcand` reconstructs a rank-`42` rational row
  space on `272` monomials. Both RREF and LLL recoveries show only vanishing
  rows at the refined root; the next rows are modular artifacts.
- Size-7 support-0 `_allcand` reconstructs `35` short low-height rows, again
  all vanishing at the refined root.
- Size-7 support-8 limited reconstructs a rank-`2` row space, also vanishing at
  the refined root.
- Combining the size-7 support-19, support-0, and support-8 vanishing rows as
  extra generators is still negative for bounded membership over `F_32003`
  through multiplier degree `4`; the degree-4 matrix is `189785 x 89089`,
  rank `39443`, augmented rank `39444`.
- Six-prime size-8 support-6 `_allcand` RREF reconstruction now succeeds:

  ```text
  code/sage/out/reconstruct_s8_79656_D6_S6_support6_allcand_32003_32009_32027_32029_32051_32057_p1400.json
  ```

- The reconstructed size-8 row space has rank `63` on `319` monomials. Its
  RREF basis evaluates at the refined root only at residual scale
  (`~1e-148` to `~1e-150` in the 500-bit evaluation).
- The automatic LLL recovery also succeeds:

  ```text
  code/sage/out/short_s8_79656_D6_S6_support6_allcand_p32003_32009_32027_32029_32051_32057_p1400.json
  ```

- The first `63` size-8 LLL rows are low-height and vanish at the 1400-bit
  refined root; the largest absolute value among them is about `3.9e-419`.
  The next row jumps to modular-artifact scale, with coefficient bound about
  `2.6e25` and a huge nonzero root value.

The practical interpretation is that the support-restricted separator route is
now much weaker: size-7 support `19`, support `0`, support `8`, and size-8
support `6` all expose stable rational row spaces, but those spaces vanish at
the refined root rather than giving a nonzero rational open separator.

The current proof target is the semialgebraic/CAD layer. The residual-component
minor relation plus active PSD forces `z0=0`; this cuts the determinant residual
component from dimension `3` to dimension `2`, degree `16`, over both
`F_32003` and `QQ`. Exact rational two-sections of this branch (`501`, `502`,
`503`) have zero GTZ-feasible real roots and expose a second square obstruction,
forcing `z1=0`. This cuts to a curve of degree `24`. An exact one-section of the
curve exposes a third square obstruction, forcing `z7=0`, and the resulting
finite stratum has dimension `0`, degree `56`, over both `F_32003` and `QQ`.
Its lex basis is
`z2^2-5,216*u0-1,z8^2,z7,z6^2-5,z5^2,z4^2-5,z3*z5*z8,z3^2,z1,z0`;
the eight real radical points all fail active PSD and inactive inequalities.
This is the strongest CAD-route signal so far.

The component-relation obstacle was substantially beaten by open-locus checks
over `QQ`.  For each target relation `R`, compute the branch plus
`a*R-1`.  Results: `R0 != 0` on the full determinant ideal has dimension `0`,
degree `608`, while the determinant branch has dimension `3`; `R1 != 0` on
`I+R0+z0` is empty; `R2 != 0` on `I+R0+z0+R1+z1` is empty.  Normal-form checks
also show `R2` is ordinary ideal membership.  The remaining task is no longer
search; it is to package these Groebner/open-locus computations into a
hand-checkable or reproducible proof certificate for the partial paper.
