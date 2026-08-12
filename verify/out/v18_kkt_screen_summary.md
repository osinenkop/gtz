# v18 KKT active-set screen summary

Date: 2026-07-31.

Status: **NUMERICALLY SUPPORTED** only.  These runs are heuristic
least-squares/KKT screens over symmetry-reduced active-set orbits; they do not
prove infeasibility of any stratum.

## Local validation

Smoke test:

- `verify/v18_kkt_screen.py --mode equality --known-smoke ...`
- Seeded at the seven certified extremals, the residuals dropped to about
  `1e-13` to `1e-14`, confirming that the stationarity formulation is consistent
  with the exact `v9`/`v11` certificates.

Low-active sweeps:

| mode | active sizes | orbit reps | candidates | notes |
|---|---:|---:|---:|---|
| counterexample | 1..4 | 32 | 0 | no tripwire |
| equality | 1..4 | 32 | 0 | no exact low-active candidate |
| counterexample | 5..6 | 137 | 0 | no tripwire |
| equality | 5..6 | 137 | 0 | no exact low-active candidate |
| counterexample | 7..8 | 410 | 0 | no tripwire |
| equality | 7..8 | 410 | 1 selected-subset hit | not genuinely 8-active |
| counterexample | 9 | 312 | 0 | one near-threshold non-candidate |
| equality | 9 | 312 | 0 | no equality stratum; best actual active size 1 |
| counterexample | 14..20 | 170 | 0 | no tripwire; all best records collapsed to actual active size 1 |
| equality | 14..20 | 170 | 0 | no equality stratum; all best records collapsed to actual active size 1 |

The equality `|A|=8` hit had residual `7.09e-13`, but its actual active set had
size 12.  Its canonical actual active mask was `32243`, matching the known TTSP
12-active orbit (`P(S(e,e),S(e,e),e,e)` / `P(S(P(S(e,e),e,e),e),e)`).  Therefore
it is a selected subset of a known equality extremal, not a genuine low-active
stratum.

The size-9 counterexample near-miss polished to:

- `loss = 0.004956328630263608`
- `t = 0.16666492665787122`
- `F = 0.16666994562575838`
- actual active size `1`

It stayed above `1/6` and did not trigger the counterexample protocol.

## MLCore sweep

MLCore job:

- job name: `gtz63-kkt-x79fea`
- project: `test-paos`
- region: `ix-m5-sm12`
- preset: `mlcore_job/kkt.mlcore-preset.yaml`
- command inside container:
  `python -u /work/v18_kkt_screen.py --mode counterexample --size-min 10 --size-max 13 --starts 4 --max-nfev 3500 --jobs 19`

Result:

- screened `1074` orbit representatives
- worker errors: `0`
- candidates: `0`
- no tripwire

Best MLCore record:

- `loss = 4.926e-03`
- selected `|A| = 12`
- `t = 0.166663747`
- `F = 0.166677230`
- spread `1.3e-05`
- inactive excess `-2.9e-06`
- stationarity `1.8e-12`
- actual active size `1`
- exact active set: `False`

Artifact reported by MLCore:

`https://devplatform.tcsbank.ru/mlcore/job-artifacts/test-paos/gtz63-kkt-x79fea/artifacts/v18_out.zip`

## High-active local sweep

The remaining counterexample active-set sizes were small enough to run locally:

- command:
  `.venv/bin/python -u verify/v18_kkt_screen.py --mode counterexample --size-min 14 --size-max 20 --starts 4 --max-nfev 3500 --jobs 19`
- screened `170` orbit representatives
- worker errors: `0`
- candidates: `0`
- no tripwire

Best record:

- selected `|A| = 15`
- `loss = 5.596e+00`
- `t = 0.152512215`
- `F = 0.217821904`
- stationarity `5.0e-07`
- actual active size `1`
- exact active set: `False`

Together with the earlier local and MLCore sweeps, this gives numerical
counterexample-screen coverage of all nonempty active sizes `1..20`.  It remains
heuristic evidence, not an exact infeasibility certificate.

## Equality sweep and finiteness probe

Additional equality-mode runs:

- `9..9`: screened `312` orbit representatives, candidates `0`.
- `10..13`: screened `1074` orbit representatives, candidates `2`.
- `14..20`: screened `170` orbit representatives, candidates `0`.

The two `10..13` candidates and the earlier `7..8` selected-subset hit were
classified by `verify/v20_finiteness_probe.py`, which reconstructs the projector
from `best_y`, computes the **actual** active set, canonicalizes it under `S_6`,
and runs a numerical sharpness test on the actual active gradients.

Result over all equality outputs:

- input records: `2135`
- low-residual equality hits: `3`
- distinct actual active-set orbits: `2`
- new active-set orbits: `0`
- non-sharp candidates: `0`
- actual active size `<= 9`: `0`

The three hits are:

- actual `|A|=12`, rank `9`, sharp, matching the known TTSP 12-active orbit
  (`P(S(e,e),S(e,e),e,e)` / `P(S(P(S(e,e),e,e),e),e)`);
- the same actual `|A|=12` orbit found from a selected `|A|=10` subset;
- actual `|A|=13`, rank `9`, sharp, matching the out-of-family
  `(5/14,9/14)` orbit.

This deepens the search in the direction relevant to finiteness: no equality
KKT hit with a new actual active orbit, low active size, or numerical non-sharp
signature was found.

## Interpretation

These runs did not find any sub-`1/6` KKT candidate.  They also show why selected
active-set screens must track the **actual** active set: a subset of a larger
known active set can satisfy stationarity and equality to high precision.

The current result is evidence only.  A proof would require an exact infeasibility
certificate, interval-certified exclusion, or a symbolic/semialgebraic argument
per active-set stratum.
