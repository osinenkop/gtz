# GTZ(6,3) — status, 2026-07-31

Current state for Pavel. Everything is on disk in `~/Documents/gtz`; no state
lives only in a running process. Full write-up: `proofs/sharp-cone-at-extremal.md`.

---

## 1. Where to resume

**The 2026-07-30 surd blocker is CLOSED.** `verify/v8_gordan.py` is superseded by
`v9_rational_certificate.py` (all six census extremals) and
`v11_seventh_exact.py` (the out-of-family seventh). The surds were eliminated by
three successive positive-scaling tricks — §3a below, and
`memory/gtz-63-rationalization-tricks.md`.

**Highest-value open question now — new, and structural.** The extremal set is
**not** exhausted by the corpus's nine scaled-star matrices: `v3` found a seventh
extremal provably outside that family, and nothing rules out more. Since part (b)
of Reformulation R quantifies over *all* KKT points, any strategy that works by
enumerating a known extremal list is on shakier ground than it appeared.

In priority order:

1. **Systematic hunt for out-of-family extremals.** Run `v3`-style unconstrained
   descent at scale, cluster minimizers by leverage pattern, PSLQ each new
   pattern, and certify with `v11`. We have exactly one out-of-family data point
   (`5/14, 9/14`) plus four TTSP patterns. The real question: **is the extremal
   set even finite?** That question was not visible before this session.
2. **Reconcile 6 vs 9** in the scaled-star census (§4).
3. **Certified radius, not just direction** — §6(ii) of the proof doc: convert
   `κ > 0` into "`F ≥ 1/6` on an explicit ball `‖Ṗ‖ ≤ ρ₀`" via a Lipschitz bound
   on `F'`. Available in closed form (3×3 blocks, algebraic eigenvalues). This is
   what Track C actually needs to excise extremal neighbourhoods analytically.
4. Only then: **part (b)**, the global gap — still the whole problem.

---

## 2. What is PROVED (exact, re-derived, machine-checked)

| Result | Script | Checks |
|---|---|---|
| Slice lemmas as *generic* polynomial identities (not just at samples) | `v1_foundations.py` | 36/37 |
| `G_PMC = (1 − sin36°)/2`, minpoly `256x⁴−512x³+304x²−48x+1` | `v1` | exact |
| The 6 active PMC triples split 4+2 into two algebraic types | `v1` | exact |
| `G_PMC > 1/6` by the *rational* certificate `169 < 405` | `v1` | exact |
| **GTZ(6,3) is TIGHT**: `F(A₀) = 1/6` exactly at the Nesterenko extremal | `v4_nesterenko.py` | 21/21 |
| **`A₀` is a SHARP minimum** on `Gr(3,6)`; KKT multipliers `(7/90)×9, 3/10` all strictly positive; `κ ≈ 7.33e−2` | `v6_sharpness.py` | 14/14 |
| **All 6 census extremals SHARP**, certificate fully rational | `v9_rational_certificate.py` | 6/6 |
| **A 7th extremal, outside the scaled-star family, is SHARP** | `v11_seventh_exact.py` | exact |

Part (a) of Reformulation R (`boundary-obstruction.md` §5) is thus **PROVED at
seven configurations**, in a form stronger than originally asserted: not merely
strict local minima but *sharp* ones, with exact certificates.

**Still OPEN, untouched:** part (b), the global gap — no KKT point of `F` has
value below `1/6`. That is the entire remaining difficulty of GTZ(6,3).

---

## 3. Reverse-engineered construction (the corpus omits it)

`boundary-obstruction.md` cites "the nine Nesterenko scaled-star matrices" but
names only four graphs and never states the weight rule. Recovered and verified
exactly against Lemma 2.1:

1. Weighted graphic-matroid / resistor-network construction: for a two-terminal
   series–parallel (TTSP) graph on 6 edges with weights `w`,
   `M̃ = W^{1/2} B_redᵀ`, `P = M̃ (M̃ᵀM̃)⁻¹ M̃ᵀ`.
2. **Leverage = weight × effective resistance**: `ℓ_e = w_e · R_eff(e)`.
   Check: `R_eff(0,1) = 5/18`, `(9/5)·R_eff(0,2) = 13/18` — matches exactly.

An extremal is such a graph plus weights tuned so `F = 1/6` exactly. `v7`
enumerates the 66 TTSP trees and solves for those weights.

### 3a. How the surds were eliminated

Three tricks, each solving a case the previous could not:

1. **Conjugation.** `S := B̃ G⁻¹ B̃ᵀ W` with `G = B̃ᵀWB̃` is a **rational**
   idempotent, and `P = W^{1/2} S W^{-1/2}` with `W^{1/2}` *diagonal*, so the
   similarity restricts to every principal block and `spec(P_TT) = spec(S_TT)`.
   The functional stays rational too: `vᵀṖ_TT v = uᵀ(W_T Ṡ_TT)u`.
2. **Per-row positive scaling** — each active gradient is defined only up to a
   positive scale.
3. **Per-column positive scaling** of the tangent basis — needed for the seventh
   extremal, where one row mixes `1, √2, √5, √10` and no row scaling can clear it.
   Rescaling `B_j → c_j B_j` is the coordinate change `z_j → z_j/c_j`, which
   preserves both the cone and the rank.

**Verdict logic:** route (A) alone is a complete proof — maximize `Σᵢ(−Lz)ᵢ` over
`{Lz ≤ 0, ‖z‖∞ ≤ 1}`; the feasible set contains `z=0` with objective 0 and the
objective is a sum of non-negative terms, so an exact optimum of 0 forces `Lz=0`,
and `rank L = 9` forces `z=0`. Route (B) (a strictly positive KKT multiplier,
Gordan's dual) is a bonus: its success gives publishable multipliers, its
**failure proves nothing**. Validated on three controls — sharp → 0; non-sharp →
positive value *with explicit descent direction*; rank-deficient → 0 but caught by
the rank test.

---

## 4. Extremal census

Six scaled-star configurations with `F = 1/6` exactly, **all SHARP**
(`verify/out/v9_certificate.json`). The census was independently reproduced on
MLCore (job `gtz63-extremals-wagucx`, `16cpu-128ram`, region `ix-m5-sm12`,
SUCCEEDED) — identical verdicts on different hardware.

| graph | weights | leverages | active | dim ker Lᵀ | route B multipliers |
|---|---|---|---|---|---|
| `P(S(e,e,e),e,e,e)` | `1,1,1,5/9,5/9,5/9` | `13/18 ×3, 5/18 ×3` | 10 | 1 | `(7/90)×9, 3/10` |
| `P(S(e,e),S(e,e),e,e)` | `1,1,1,1,5/8,5/8` | `11/18 ×4, 5/18 ×2` | 12 | 3 | `(3/20)×4, (1/20)×8` |
| `P(S(P(e,e,e),e),S(e,e))` | `1,1,1,9/5,9/5,9/5` | `5/18 ×3, 13/18 ×3` | 10 | 1 | — |
| `P(S(P(S(e,e),e,e),e),e)` | `1,1,5/8,5/8,1,1` | `11/18 ×4, 5/18 ×2` | 12 | 3 | — |
| `P(S(P(e,e),P(e,e)),S(e,e))` | `1,1,1,1,8/5,8/5` | `7/18 ×4, 13/18 ×2` | 12 | 3 | — |
| `P(S(P(e,e),e),S(P(e,e),e))` | `1,1,8/5,1,1,8/5` | `7/18 ×4, 13/18 ×2` | 12 | 3 | — |

All six: `rank L = 9`, exact LP value 0, `λ_min` simple on every active block.
The first two reproduce `v6`'s multipliers exactly, cross-validating through an
entirely different code path. A "—" in the last column means route (B) found no
witness; route (A) certified it regardless.

**Open: six, not nine.** Either the one-symbol-per-edge-orbit weight
parametrization is too coarse, or some of the corpus's nine are `O(3)`/relabelling
duplicates. Does not affect any certificate above, but bounds the scope of the
phrase "every extremal".

### 4a. The seventh extremal (`verify/out/v11_seventh.json`)

Leverages exactly `5/14 ×3`, `9/14 ×3`; **13** active triples, 7 singular; entries
in `ℚ(√2,√5)`, exactly five distinct off-diagonal magnitudes:

```
5/14,   √5/21,   5√5/42,   5√2/42,   √10/14
```

Verified exactly: `Pᵀ=P`, `P²=P`, `tr P=3`, `F = 1/6`, nothing strictly inside
`(0,1/6)`, `λ_min` simple everywhere. `rank L = 9`, exact LP value 0 → **SHARP**.
Floats were used only to *guess* the entries; every certified statement is a
theorem about the exact object. Exhaustive search over all 66 TTSP trees found
**no** graph with this leverage pattern → genuinely out of family.

---

## 5. Corrections made during this work

Both were errors of mine, corrected by later exact computation and recorded rather
than quietly patched (brief §8 rule 1):

1. **Ambient `κ ≈ 6.59e−3` was wrong.** Minimized over 18-dim `A`-space, which
   contains the 3-dim `O(3)` gauge `A→AΩ` along which `F` is *exactly constant* —
   so the true ambient min is 0 and that number was just where Nelder–Mead
   stalled. Gauge-free value on `Gr(3,6)`: `κ ≈ 7.33e−2`.
2. **"Descent never reaches the 1/6 floor" was wrong.** A `log-sum-exp` overflow
   (`exp(4000λ) → inf`) was feeding `nan` to the optimizer. After the shift fix,
   357/500 starts land within `1e−4` of `1/6`; `v3`'s PSLQ gives minpoly
   `36x²−12x+1 = (6x−1)²`, i.e. exactly `1/6`. So seeding branch-and-bound with
   known extremals is good hygiene, not a necessity. Always shift log-sum-exp.

Tripwire unaffected: **0 violations** across ~2,600 full-space descents, 24,000
shell perturbations, 200 kicked descents from `A₀`, 20,000 direction probes.

---

## 6. Environment

- `~/Documents/gtz/.venv` — Python 3.12.3, numpy 2.5.1, sympy 1.14.0,
  scipy 1.18.0, mpmath, python-flint (Arb, for future interval work).
  Built with `uv pip install --system-certs` (corporate TLS interception breaks
  plain `uv`/`pip`; `python3-venv` lacks `ensurepip` on this box).
- 22 cores. MLCore used once, as an independent cross-check of `v7`; it agreed.
  Everything else fits locally.
- CLIs: `mlc` v1.0.7 (`/usr/local/bin/mlc`, project `test-paos`),
  `dp` v13.35.0 (`~/.local/bin/dp`). May need `mlc hello` / `dp update` to
  refresh tokens.

Reproduce everything (offline):

```bash
cd ~/Documents/gtz
.venv/bin/python -u verify/v1_foundations.py          # 36/37
.venv/bin/python -u verify/v4_nesterenko.py           # 21/21
.venv/bin/python -u verify/v6_sharpness.py            # 14/14
.venv/bin/python -u verify/v7_all_extremals.py        # census, ~15 min
.venv/bin/python -u verify/v9_rational_certificate.py # 6/6 SHARP
.venv/bin/python -u verify/v5_tightness.py 500        # seed 20260801
```

`v11_seventh_exact.py` reads the polished projector from
`verify/data/P514_seventh.npy` (now a committed artifact, no longer a `/tmp` file;
it falls back to `/tmp/P514.npy` if present). To regenerate from scratch:
re-polish from `verify/out/v3_extremal.json` with ~14 rounds of Nelder–Mead on `F`
until `|F − 1/6| < 1e−14`.

`verify/v8_gordan.py` is kept only as a record of the failed rationalization
approach; it does not run. `v2_search.py` and `v3_extremal.py` are the search
drivers.
