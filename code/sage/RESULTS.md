# Sage/Singular probe results

These are local smoke-test results for the GTZ(6,3) semialgebraic route.

## Environment

- SageMath: `10.9` in `~/miniforge3/envs/sage`
- Singular in that environment: `4.4.1`

Run commands as:

```bash
~/miniforge3/bin/mamba run -n sage python code/sage/export_semialgebraic_system.py ...
~/miniforge3/bin/mamba run -n sage Singular -q path/to/script.sing
```

## Completed probes

| System | File prefix | Variables | Equations | Max degree | Result |
| --- | --- | ---: | ---: | ---: | --- |
| One active saturated determinant, active index `0` | `code/sage/out/smoke_det_A0` | 9 | 1 | 6 | Sage dimension `8` |
| Two active saturated determinants, active indices `0,1` | `code/sage/out/prefix_det_A0_1` | 9 | 2 | 6 | Sage dimension `7` |
| Two active saturated determinants, active indices `0,12` | `code/sage/out/pair_det_A0_12` | 9 | 2 | 6 | Sage dimension `7` |
| Three active saturated determinants, active indices `0,1,2` | `code/sage/out/prefix_det_A0_2` | 9 | 3 | 6 | Sage dimension `6` |
| One active cofactor nonsharp, active index `0`, row patch `01` | `code/sage/out/smoke_nonsharp_A0` | 20 | 4 | 35 | Export succeeds, but file is about 3.5 MB and Singular parsing/computation is already heavy |
| One active explicit-kernel nonsharp, active index `0` | `code/sage/out/smoke_kernel_A0` | 22 | 6 | 12 | Export succeeds; Singular parses the generated file |
| Known base saturated determinant active set `P(S(e,e,e),e,e,e)` | `code/sage/out/known_base_det` | 9 | 10 | 6 | Export succeeds; Singular parses the generated file |
| Known base saturated determinant active set with `d` inverted | `code/sage/out/known_base_det_invd` | 10 | 11 | 7 | Export succeeds; Singular parses the generated file |
| Known base explicit-kernel nonsharp active set `P(S(e,e,e),e,e,e)` | `code/sage/out/known_base_kernel` | 58 | 51 | 13 | Export succeeds; Singular parses the generated file |

## Saturation finding

The raw determinant numerator has a universal factor:

```text
det(6*N_TT-d*I) = d^2 * F_T
```

where `d=det(Y^T Y)` and `deg(F_T)=6`.  Since the standard chart assumes
`d!=0`, the active determinant equation should be `F_T=0`.  For exact algebraic
probes of the chart itself, add the inverse equation `u0*d-1=0`; otherwise the
polynomial ideal may still contain components on the boundary `d=0`.

Before this saturation, two raw determinant equations still reported dimension
`8`; after saturation, sampled two-active systems report dimension `7`.

## Modular section probes

The helper `probe_determinant_ideals.py` can run Singular over a prime field and
add deterministic random affine linear sections.  These are evidence/probes,
not a characteristic-zero proof unless later certified.

For the known base active set `P(S(e,e,e),e,e,e)`, over `F_32003` with `d`
inverted:

| Probe | Result |
| --- | --- |
| 3 generic sections | dimension `0`, section degree `16` |
| 4 generic sections | dimension `-1` |
| 3 zero-sum integer sections, seed `401` | dimension `0`, section degree `18` |

This indicates that the determinant-only active locus in the actual chart has
top dimension `3`; determinant equalities alone are not a finiteness proof.

For four-subsets of the same 10-active set over `F_32003` without `d` inversion,
using five generic sections:

| Inferred top dimension | Count |
| ---: | ---: |
| 5 | 204 |
| 6 | 6 |

The six exceptional four-subsets are exactly the evident dependent patterns:
the four triples through one fixed pair, or the four triples on one four-element
vertex set inside the known-base pattern.  Spot checks with `d` inverted show
the same behavior for one exceptional and one generic four-subset.

For five-subsets of the same active set over `F_32003` without `d` inversion,
using four generic sections:

| Inferred top dimension | Count |
| ---: | ---: |
| 4 | 216 |
| 5 | 36 |

## Exact lex slice certificate

For the known base active set with `d` inverted, the zero-sum integer section
seed `401` was also computed over `QQ` with Singular's modular Groebner backend:

```bash
~/miniforge3/bin/mamba run -n sage python code/sage/lexify_determinant_slice.py \
  --known-label 'P(S(e,e,e),e,e,e)' --invert-d --method modslimgb \
  --characteristic 0 --linear-sections 3 --linear-section-mode zero-sum-z \
  --seed 401 --timeout 600 \
  --out-prefix code/sage/out/lex_known_base_invd_QQ_zerosum401_sec3_z0
```

The lex order was `z1,z2,z3,z4,z5,z6,z7,z8,u0,z0`.  The computation finished in
about 60 seconds, with quotient dimension `18` and a lex basis of size `10`.
The first lex polynomial is univariate in `z0`, degree `18`, and factors over
`QQ` as

```text
(z0^2 - 5/9) * q16(z0)
```

up to a nonzero rational scalar, where `q16` is irreducible over `QQ`.  The
degree-16 factor has 12 real roots and 4 non-real roots.

The reconstructed real roots of the triangular basis give 14 real section
points.  Exactly two pass the active/inactive GTZ inequalities.  They are the
two roots of `z0^2 - 5/9`; modulo this quadratic factor the lex relations imply

```text
z1=z2=...=z8=z0,   u0=1/6.
```

For the other real roots, all lying on `q16`, Sage's algebraic real field
certifies a uniform active-PSD failure: for active triple index `1`, i.e.
`(0,1,3)`, the principal minor on local rows `(0,2)` of
`6*N_TT-d*I` reduces modulo `q16` to

```text
-46656*z0^2.
```

Since `q16(0) != 0`, this is strictly negative at every real root of `q16`.
The same certificate file also checks `d>0` at those roots.  Thus this rational
three-section has no semialgebraically feasible points beyond the two sign
copies of the base extremal.  This is a slice certificate, not yet a global
certificate for the whole three-dimensional determinant locus.

The same exact computation was repeated for two further zero-sum integer
sections.  All three runs show the same factor and sign pattern:

| Seed | Lex quotient dimension | Factor degrees | Real roots | Feasible roots | Residual real roots | Minor witness modulo residual factor |
| ---: | ---: | --- | ---: | ---: | ---: | --- |
| `401` | `18` | `2+16` | `14` | `2` | `12` | `-46656*z0^2` |
| `402` | `18` | `2+16` | `14` | `2` | `12` | `-46656*z0^2` |
| `403` | `18` | `2+16` | `14` | `2` | `12` | `-46656*z0^2` |
| `404` | `18` | `2+16` | `10` | `2` | `8` | `-46656*z0^2` |

Screening all active `2x2` principal minors on these four lex slices gives an
identical obstruction pattern.  For every active triple containing one lower
chart row, both `2x2` minors involving that lower row are strictly negative at
all real roots of the residual degree-16 factor.  This accounts for 18 of the
30 active `2x2` principal minors.  The active triple `(0,1,2)` contributes no
such obstruction because it is the fixed identity-chart triple.

Symbolically, the representative obstructing minor for triple `(0,1,3)` and
local rows `(0,2)` factors as `-d*C`, where `d=det(Y^T Y)` and `C` is a
degree-6 polynomial not equal to any one active determinant equation.  Thus the
slice obstruction is genuinely from the active PSD inequalities, not merely
from reusing one of the determinant equalities.

A modular probe over `F_32003` also supports the component interpretation.  If
the candidate relation

```text
minor((0,1,3); rows 0,2) + 46656*z0^2 = 0
```

is added to the known-base determinant ideal with `d` inverted, then the same
four zero-sum three-sections have degree `16`, not `18`.  This is exactly the
degree of the non-known factor in the lex computations, so the relation appears
to select the residual component sampled by these sections.

The no-section modular computation makes this interpretation sharper.  Over
`F_32003`, the determinant ideal with `d` inverted has top dimension `3` and
top degree `16`; adding the minor relation still gives top dimension `3` and
top degree `16`:

| System over `F_32003` | Linear sections | Dimension | Top degree | Basis size |
| --- | ---: | ---: | ---: | ---: |
| determinant ideal | `0` | `3` | `16` | `1472` |
| determinant ideal plus minor relation | `0` | `3` | `16` | `498` |

Thus the exact zero-sum section degree `18` should be read as `16` points from
the positive-dimensional top component plus the two isolated known sign points
forced into the section.  The candidate minor relation appears to vanish on the
top component while excluding those isolated sign points.

The corresponding no-section comparison was also completed over `QQ` with
Singular's modular Groebner backend:

| System over `QQ` | Linear sections | Dimension | Top degree | Basis size | Elapsed |
| --- | ---: | ---: | ---: | ---: | ---: |
| determinant ideal | `0` | `3` | `16` | `1472` | `370s` |
| determinant ideal plus minor relation | `0` | `3` | `16` | `498` | `367s` |

This gives the same component-degree signal in characteristic zero: the
relation-augmented ideal has the same top-dimensional degree as the determinant
ideal itself, while being a closed sublocus of it.

One direct open-set test was also attempted.  With separator `z1-z0`, the
separator section over `QQ` has dimension `2` and degree `16`, so no
top-dimensional component is contained in `z1=z0`.  The harder rational test
with both `(z1-z0)^{-1}` and `relation^{-1}` timed out after about 913 seconds.
Over `F_32003`, the same open-relation-nonzero test finishes with dimension
`0`, so the relation is nonzero only on lower-dimensional leftovers in that
modular model.

## CAD obstruction stratum

The helper `probe_cad_obstruction_strata.py` now probes the smaller
semialgebraic branch suggested by the active-PSD obstruction.  If the component
relation

```text
minor((0,1,3); rows 0,2) + 46656*z0^2 = 0
```

holds on the residual top component, then active PSD forces `z0=0` over the
reals.  Adding this forced equality gives a genuine algebraic dimension drop:

| System | Field/backend | Dimension | Top degree | Basis size | Elapsed |
| --- | --- | ---: | ---: | ---: | ---: |
| determinant ideal | `F_32003`/`slimgb` | `3` | `16` | `1472` | `35s` |
| determinant ideal plus relation | `F_32003`/`slimgb` | `3` | `16` | `498` | `37s` |
| determinant ideal plus relation plus `z0=0` | `F_32003`/`slimgb` | `2` | `16` | `102` | `0.4s` |
| determinant ideal plus relation plus `z0=0` | `QQ`/`modslimgb` | `2` | `16` | `102` | `4s` |

The exact outputs are:

```text
code/sage/out/cad_obstruction_known_base_invd_p32003_z0.json
code/sage/out/cad_obstruction_known_base_invd_QQ_z0_forced.json
```

`lexify_determinant_slice.py` was also extended with `--add-minor-relation` and
`--extra-equalities`, so exact slices of this forced stratum can be computed
with the same FGLM pipeline.  Three independent rational zero-sum two-sections
of the dimension-2 forced stratum were computed over `QQ`.  Each has lex
quotient dimension `16`; none has a real point passing the GTZ active/inactive
inequalities:

| Seed | Lex quotient dimension | Real roots | Pass active PSD | Pass inactive | Pass both |
| ---: | ---: | ---: | ---: | ---: | ---: |
| `501` | `16` | `12` | `0` | `0` | `0` |
| `502` | `16` | `8` | `0` | `0` | `0` |
| `503` | `16` | `8` | `0` | `0` | `0` |

The best active-PSD margin among these real slice roots is still
`lambda_min-1/6 = -1/6`, and the inactive maximum is at least `2/3` above the
threshold.  This is still a section probe, not a global certificate, but it
turns the CAD target from the original three-dimensional determinant locus into
a two-dimensional boundary stratum with stable exact slices.  The next screen
therefore searched these slices for another PSD-forced relation.

The exact minor screen on the forced two-section slices found a stable second
square obstruction.  On all three seeds, the same 16 active `2x2` principal
minors are uniformly negative on all real slice roots; in particular

```text
minor((0,1,3); rows 1,2) = -46656*z1^2
```

on the sampled forced branch.  Adding the relation
`minor((0,1,3); rows 1,2) + 46656*z1^2 = 0` preserves the `z0=0` branch
dimension/degree, while active PSD then forces `z1=0` and cuts to a curve:

| System | Field/backend | Dimension | Top degree | Basis size |
| --- | --- | ---: | ---: | ---: |
| relation `z0` branch plus second relation | `F_32003`/`slimgb` | `2` | `16` | `66` |
| relation `z0` branch plus second relation | `QQ`/`modslimgb` | `2` | `16` | `66` |
| plus forced `z1=0` | `F_32003`/`slimgb` | `1` | `24` | `19` |
| plus forced `z1=0` | `QQ`/`modslimgb` | `1` | `24` | `19` |

An exact one-section of this curve over `QQ` has lex quotient dimension `24`;
after changing the lex parameter to `z7`, the real-root classifier finds 16
real roots and none pass the GTZ inequalities.  Screening active minors on this
curve section produces a third square obstruction:

```text
minor((0,1,4); rows 0,2) = -46656*z7^2.
```

Adding
`minor((0,1,4); rows 0,2) + 46656*z7^2 = 0` preserves the curve, while active
PSD forces `z7=0` and cuts to a finite stratum:

| System | Field/backend | Dimension | Top degree | Basis size |
| --- | --- | ---: | ---: | ---: |
| curve branch plus third relation | `F_32003`/`slimgb` | `1` | `24` | `19` |
| plus forced `z7=0` | `F_32003`/`slimgb` | `0` | `56` | `11` |
| plus forced `z7=0` | `QQ`/`modslimgb` | `0` | `56` | `11` |

The exact lex basis of the final forced finite stratum is especially simple:

```text
z2^2-5, 216*u0-1, z8^2, z7, z6^2-5, z5^2,
z4^2-5, z3*z5*z8, z3^2, z1, z0
```

Thus its real radical has

```text
z0=z1=z3=z5=z7=z8=0,    z2^2=z4^2=z6^2=5,    u0=1/216.
```

The eight real sign choices for `(z2,z4,z6)` all fail active PSD
(`lambda_min-1/6 = -1/6`) and also fail the inactive inequalities
(`inactive_max_lambda-1/6 = 2/3`).  The CAD route has therefore reduced the
sampled residual branch to a fully checked finite endpoint.  What remains for a
formal proof is to certify the three component relations on the corresponding
successive branches, ideally by saturation/open-set computations over `QQ` or
by replacing those computations with hand-verifiable Groebner reductions.

### Component-relation certificates

The open-locus formulation beats the main obstacle.  For a branch ideal `B` and
target relation `R`, the helper `probe_cad_relation_open_locus.py` computes
`B + <a*R-1>`.  If this has dimension smaller than `dim(B)`, then `R` vanishes
on every top-dimensional component of `B`.  This avoids the older two-inverse
separator test that timed out.

Over `QQ` with Singular's modular backend:

| Target relation | Branch `B` | `dim(B)` | `dim(B + <aR-1>)` | Degree of nonzero locus | Meaning |
| --- | --- | ---: | ---: | ---: | --- |
| `R0 = minor((0,1,3); rows 0,2)+46656*z0^2` | determinant ideal | `3` | `0` | `608` | `R0` holds on all top-dimensional determinant components |
| `R1 = minor((0,1,3); rows 1,2)+46656*z1^2` | `I+R0+z0` | `2` | `-1` | `0` | `R1` holds on the whole branch |
| `R2 = minor((0,1,4); rows 0,2)+46656*z7^2` | `I+R0+z0+R1+z1` | `1` | `-1` | `0` | `R2` holds on the whole branch |

The corresponding artifacts are:

```text
code/sage/out/cad_open_R0_on_I_QQ.json
code/sage/out/cad_open_R1_on_z0_QQ.json
code/sage/out/cad_open_R2_on_z0z1_QQ.json
```

Normal-form probes provide a useful cross-check:

| Target | Branch | Normal form over `QQ` |
| --- | --- | --- |
| `R1` | `I+R0+z0` | nonzero, 5 terms: `88992*z2^2+88992*z5^2+88992*z8^2+821192256/7*u0-6916536/7` |
| `R2` | `I+R0+z0+R1+z1` | `0` |

Thus `R2` is ordinary ideal membership, while `R1` is a genuine
component/open-locus statement.  `R0` also needs the component statement:
its nonzero locus is finite, not empty.  This is enough for the PSD-forcing
cascade, since the proof only needs the relations on top-dimensional real
branches; the finite `R0 != 0` leftovers do not contribute a residual
top-dimensional semialgebraic component.

### Finite `R0 != 0` residue

The finite residue from the first open-locus step has now been split.  The
helper `lexify_determinant_slice.py` was extended with
`--target-nonzero-relation`, and `reconstruct_univariate_crt.py` reconstructs
small univariate eliminants from modular lex JSON files.

The open residue

```text
I0 + <a*(minor((0,1,3); rows 0,2) + 46656*z0^2)-1>
```

was computed directly over `QQ` with Singular's modular backend:

```text
code/sage/out/lex_cad_R0_nonzero_QQmod_z0.json
```

It has lex quotient dimension `608`, and its first lex polynomial is

```text
12252303*z0^10 - 22643145*z0^8 - 17737650*z0^6
  + 61988625*z0^4 - 41649375*z0^2 + 8556250
```

or, monically,

```text
z0^10 - 815/441*z0^8 - 93850/64827*z0^6
  + 2295875/453789*z0^4 - 13883125/4084101*z0^2
  + 8556250/12252303.
```

The same eliminant was reconstructed independently from six modular lex bases
over `32003,32009,32027,32029,32051,32057` using

```text
code/sage/out/lex_cad_R0_nonzero_z0_eliminant_crt_6primes.json
```

Its factorization is

```text
(z0^2 - 5/9)
(z0^2 - 10/21)
(z0^2 + 5/3)
(z0^4 - 365/147*z0^2 + 34225/21609).
```

The last two factors have no real zeros: the quartic is a quadratic in
`z0^2` with discriminant `-25/147`.

The real factor `z0^2=5/9` was computed directly over `QQ` with Singular's
modular backend:

```text
code/sage/out/lex_cad_R0_nonzero_QQmod_z0sq_5_9_z1.json
```

It has lex quotient dimension `32` and lex basis

```text
9*z1^2-5, 25980*a_rel-1, 6*u0-1,
9*z8^2-5, 9*z7^2-5, 9*z6^2-5, 9*z5^2-5,
5*z4-9*z5*z7*z8, 5*z3-9*z5*z6*z8,
5*z2-9*z7*z8*z1, 5*z0-9*z6*z7*z1.
```

Thus this real fiber is exactly a 32-point sign family.  Numerical substitution
using the exact formulas finds `32` distinct projectors; all have actual active
set

```text
[0,1,2,3,4,5,6,10,11,12]
```

and all pass the active and inactive GTZ inequalities with `F-1/6 = 0`.

The other real `z0` factor `z0^2=10/21` was also computed directly over `QQ`:

```text
code/sage/out/lex_cad_R0_nonzero_QQmod_z0sq_10_21_z1.json
```

It has quotient dimension `256`; its `z1` eliminant is

```text
1361367*z1^8 - 1759590*z1^6 - 2948400*z1^4
  + 5249625*z1^2 - 1711250
```

and factors as

```text
(z1^2 - 10/21)
(z1^2 + 5/3)
(z1^4 - 365/147*z1^2 + 34225/21609).
```

Again the latter two factors have no real zeros.  On the only real subfactor
`z1^2=10/21`, the exact `QQ` lex computation

```text
code/sage/out/lex_cad_R0_nonzero_QQmod_z0z1sq_10_21_z2.json
```

has quotient dimension `64` and first lex polynomial

```text
3*z2^2+5.
```

Therefore the `z0^2=10/21` branch has no real points.  Combining these checks,
the finite `R0 != 0` residue contributes only the known 32-point sign family in
the standard chart.  Together with the CAD cascade above, the base active set is
now locally accounted for in this chart; the remaining work is propagation to
the other active-set orbits and charts.

## Margin tie-system probe

The helper `probe_margin_tie_system.py` was added for low-active margin
minimizers.  It builds the determinant equality system

```text
det(6 P_TT - I) = 0      for active triples T,
det(6 P_TT - q I) = 0    for outside triples tied at the margin maximum,
```

in the standard chart with `d=det(Y^T Y)` inverted.  Here `q` is six times the
tied outside eigenvalue, so the margin is `(q-1)/6`.

The symbolic shifted determinant still has the universal `d^2` factor:

| Polynomial | Raw degree | Saturated degree |
| --- | ---: | ---: |
| `det(6*N_TT - q*d*I)` | 21 | 9 |

This makes the exact target substantially smaller than the raw determinant
system.  However, the bare active-plus-tie ideal is still too broad.  For the
tightest numerical margin point found so far, size-7 canon `14113` with active
indices `[0,1,2,6,14,16,18]` and outside ties `[10,11,12]`, all of the following
modular probes over `F_32003` timed out:

| Linear sections | Timeout |
| ---: | ---: |
| 0 | 300s |
| 1 | 180s |
| 2 | 180s |
| 3 | 180s |

The conclusion is not that the route is dead, but that tie determinants alone
are too broad for plain Groebner.

The determinant-gradient KKT system was then implemented in
`probe_margin_kkt_determinant.py`.  It adds multiplier variables and the
stationarity equations for minimising `q` subject to the active and tie
determinant equations.  The numerical KKT diagnostic in
`verify/v34_margin_kkt.py` shows strict convex outside weights at the minimizer,
but the direct algebraic KKT systems still time out over `F_32003`:

| System | Variables | Equations | Result |
| --- | ---: | ---: | --- |
| size-7 `14113`, tie `[10]` | 19 | 19 | timeout, 240s |
| size-7 `14113`, ties `[10,11]` | 20 | 20 | timeout, 240s |
| size-7 `14113`, ties `[10,11,12]` | 21 | 21 | timeout, 300s |
| size-8 `14117`, ties `[5,6]` | 21 | 21 | timeout, 240s |

Reordering multiplier variables first did not change the size-7 two-tie result.
For strata with `|A|+|B|=10`, the active-plus-tie equations are already square in
`(z,q)`, so the KKT multiplier block is largely automatic at regular isolated
points.  The next exact implementation should exploit the linear multiplier
block explicitly, or solve the active-plus-tie zero-dimensional systems with a
more specialized method.

The explicit rank-minor version of multiplier elimination was started as
`probe_margin_kkt_rank.py`.  For the size-7 `14113` two-tie case, even exporting
the augmented determinant did not finish in several minutes, so this is not the
right representation either: the minors are too high degree.

The practical specialized solver is `refine_margin_tie_root.py`.  It keeps the
square active-plus-tie determinant system in `(z,q)`, evaluates exact symbolic
Jacobians, and performs high-precision Newton refinement from the `v30` point.

| System | Precision | Residual | Jacobian condition | Margin |
| --- | ---: | ---: | ---: | ---: |
| size-7 `14113`, ties `[10,11,12]` | 500 bits | `1.37e-148` | `5.53` | `0.0219593612977106471453390091112` |
| size-8 `14117`, ties `[5,6]` | 400 bits | `8.98e-119` | `6.02` | `0.0329381178422993768053283604786` |

The refined roots pass the intended branch check: active triples have least
eigenvalue exactly `1/6`, the listed outside triples tie at `q/6`, and all other
outside triples are below.  Low-height `algdep` did not find a simple polynomial
for `q` in these tests (degree up to 16/14, height at most `1e6`).  This points to
interval/Krawczyk or rational-univariate certification of isolated numerical
roots rather than closed-form recognition.

The interval certification layer is `certify_margin_tie_krawczyk.py`.  It uses
Sage interval arithmetic to verify a Krawczyk inclusion for the same square
systems and then checks the spectral branch inequalities on the certified root
box `K(X)`.  The batch driver `batch_margin_tie_certify.py` enumerates the
square tie-subsystems for the saved `v30` points.

The six non-over-tied tested size-7/8 minimizers now have interval certificates:

| System | Equation ties | Box radius | Branch checks | Certified margin interval |
| --- | --- | ---: | --- | --- |
| size-7 `13105` | `[9,10,18]` | `1e-30` | all true | `[0.04464692543123093116213202070744, 0.04464692543123093116213202070777]` |
| size-7 `14113` | `[10,11,12]` | `1e-30` | all true | `[0.02195936129771064714533900911105, 0.02195936129771064714533900911138]` |
| size-7 `78595` | `[4,13,19]` | `1e-30` | all true | `[0.04521367233867756678261877566381, 0.04521367233867756678261877566414]` |
| size-7 `78601` | `[4,7,14]` | `1e-30` | all true | `[0.05150617318963458614354853172380, 0.05150617318963458614354853172414]` |
| size-8 `13107` | `[11,13]` | `1e-30` | all true | `[0.09820775013166748461073489955964, 0.09820775013166748461073489955997]` |
| size-8 `14117` | `[5,6]` | `1e-30` | all true | `[0.03293811784229937680532836047841, 0.03293811784229937680532836047875]` |

The previous two certificates have tiny Krawczyk images relative to the box
(`2.7e-28` and `3.7e-28`); the new cases show the same qualitative behavior and
certify that each root lies on the intended positive-margin spectral branch.

The remaining observed low-active margin points are over-tied:

| Case | Numerical outside ties | Square-subset result |
| --- | --- | --- |
| size-6 `78593` | `[5,7,12,14,18]` | no Krawczyk inclusion for any 4-tie subset; Jacobians are very ill-conditioned |
| size-7 `78612` | `[2,7,10,19]` | two 3-tie subsets have Krawczyk roots, but the fourth tie is only enclosed as a tiny interval around zero (`upper ~ 1e-54` after checking on `K(X)`) |
| size-8 `78613` | `[14,15,17,18]` | no 2-tie subset gives a usable Krawczyk inclusion |
| size-8 `79656` | `[10,11,15]` | two 2-tie subsets have Krawczyk roots, but the third tie is only enclosed as a tiny interval around zero (`upper ~ 1e-54` after checking on `K(X)`) |

This is a useful failure mode: it does not refute the positive-margin mechanism;
it says that over-tied cases need a genuine overdetermined/local-membership
certificate, rather than certification of an arbitrary square subset alone.
A direct linear-span test over `QQ` also failed for the regular over-tied cases,
using `probe_extra_tie_dependency.py`:

| Case | Selected ties | Extra tie | Rank change |
| --- | --- | --- | --- |
| size-7 `78612` | `[2,7,10]` | `19` | `10 -> 11` |
| size-7 `78612` | `[7,10,19]` | `2` | `10 -> 11` |
| size-8 `79656` | `[10,11]` | `15` | `10 -> 11` |
| size-8 `79656` | `[11,15]` | `10` | `10 -> 11` |

Thus the extra determinant is not a global linear combination of the selected
determinant equations.  The output is
`code/sage/out/extra_tie_dependency_regular_overtied.json`.

The overdetermined numerical layer was then made explicit in
`refine_margin_overtie_root.py`.  It runs high-precision Gauss-Newton on all
active equations plus all numerical outside ties.  This separates the over-tied
cases into full-rank common roots and a genuinely singular case:

| Case | Equations / variables | Residual | Min singular value | Condition | Margin |
| --- | ---: | ---: | ---: | ---: | ---: |
| size-7 `78612`, ties `[2,7,10,19]` | `11 / 10` | `2.29e-178` | `4.22e1` | `1.23e1` | `0.0426033908353660797918571745413` |
| size-8 `79656`, ties `[10,11,15]` | `11 / 10` | `1.85e-178` | `9.79e0` | `5.43e1` | `0.0525669822295852240401359267707` |
| size-8 `78613`, ties `[14,15,17,18]` | `12 / 10` | `3.51e-178` | `7.23e1` | `2.27e1` | `0.0529979913139675476441148305184` |
| size-6 `78593`, ties `[5,7,12,14,18]` | `11 / 10` | `2.01e-25` | `1.85e-13` | `7.81e15` | `0.0301159727944298842572390384476` |

The global all-tie Singular probes still time out at 120 seconds for size-7
`78612`, size-8 `79656`, size-8 `78613`, and size-6 `78593`, so this did not turn
into a global Groebner/RUR computation.

For the three full-rank over-tied size-7/8 roots, `screen_overtie_square_bases.py`
finds well-conditioned 10-row square bases that include all outside tie
equations.  `certify_overtie_square_basis.py` then proves a Krawczyk inclusion
for those bases and evaluates the omitted equations on the certified root box:

| Case | Certified square basis | Omitted equations | Omitted residual max | Spectral branch |
| --- | --- | --- | ---: | --- |
| size-7 `78612` | all four ties + all active except `013` | `active:013` | `1.31e-54` | true |
| size-8 `79656` | all three ties + all active except `012` | `active:012` | `1.44e-54` | true |
| size-8 `78613` | all four ties + active `013,024,035,045,125,234` | `active:012`, `active:134` | `8.41e-54` | active branch still degenerate |

For `78613`, the failed branch inequality is not macroscopic; the lowest active
principal minor and the trace/e2 spectral coefficient are both at the `1e-52`
scale.  This points to an additional PSD-boundary equality, not to an outside
overshoot.

The two clean over-tied cases were stress-tested at higher precision.  Refining
the full over-tied systems at 1400 bits reduced the full residual to about
`1e-419`; rerunning the square-basis Krawczyk certificates at 1600 bits with
initial radius `1e-120` reduced the omitted determinant intervals to about
`1e-234`, while still containing zero:

| Case | High-precision residual | Omitted residual interval scale |
| --- | ---: | ---: |
| size-7 `78612` | `5.55e-419` | `1.31e-234` |
| size-8 `79656` | `3.24e-419` | `1.44e-234` |

This is strong numerical-algebraic evidence that the omitted determinants vanish
exactly at these local roots, but it is not by itself an exact proof.

The bounded ideal-membership probe `probe_bounded_ideal_membership.py` tests
whether the omitted determinant, or a separated multiple of it, is a polynomial
combination of the selected ten equations with multiplier degree at most `D`.
The clean cases are negative for the simple global/separator targets below:

| Target | Field | Max multiplier degree | Result |
| --- | --- | ---: | --- |
| omitted determinant `h` | `QQ` | 2 | no membership |
| omitted determinant `h` | `F_32003` | 5 | no membership |
| `h^2` | `F_32003` | 4 | no membership |
| chart-denominator multiple `d*h` | `F_32003` | 6 | no membership |
| chart-denominator multiple `d^2*h` | `F_32003` | 5 | no membership |

So the remaining proof is probably not a low-degree global ideal-membership
identity with an obvious separator.  However, the new
`probe_local_separator_membership.py` search does find genuine modular local
separators.  It solves

```text
s*h = sum_i a_i f_i,   deg(a_i) <= D,   deg(s) <= S,
```

with `s` unknown.  For both clean over-tied cases, `D=4` and `S<=3` are
negative, as are `D=5` with `S<=3`.  At `D=5, S=4`, over `F_32003`, the
separator rank defect becomes positive:

| Case | Rank defect | Candidate separators | Nonzero at refined root |
| --- | ---: | ---: | ---: |
| size-7 `78612`, omitted `active:013` | 5 | 12 | 12 |
| size-8 `79656`, omitted `active:012` | 6 | 21 | 21 |

For size-7, one modular quartic candidate has two useful lifts to compare:

| Lift | Exact membership over `QQ` | Separator on certified root box |
| --- | --- | --- |
| scaled rational factor (`candidate9_factor_scaled5`) | yes: rank `29920`, augmented rank `29920` | inconclusive/vanishing; center value about `5.5e-420`, interval contains zero |
| raw residue integer lift (`candidate9_rawlift`) | no: rank `29920`, augmented rank `29921` | rigorously nonzero; interval is near `-7606.296` |

So the first characteristic-zero lift attempt did not yet close the proof.  It
does show that the modular signal is not noise, but the usable nonzero separator
has not been lifted to a `QQ` membership identity.

The modular separator spaces were then recomputed over `F_32009` and `F_32027`
using ML Core project `aida`.  For size-7 the rank defect remains `5` over all
three primes; for size-8 the rank defect remains `6`.  Row-reducing the saved
separator parts gives identical pivot patterns across `32003,32009,32027`, so
CRT/rational reconstruction is stable:

| Case | Reconstructed separator-space rank | Pivot monomials | Root evaluation at 1400 bits |
| --- | ---: | --- | --- |
| size-7 `78612`, `D=5,S=4` | 5 | `1,z8^2,z7^2,z6^2,z5^2` | all basis vectors are at `~1e-420` residual scale |
| size-8 `79656`, `D=5,S=4` | 6 | `1,z8^2,z7^2,z6^2,z5^2,z5^2*z7^2` | all basis vectors are at `~1e-420` residual scale |

Thus the characteristic-zero separator subspaces visible at `D=5,S=4` appear
to vanish on the refined local roots.  This explains why low-height
reconstructions such as `candidate4_crt115`, `candidate5_factor_scaled20`,
`candidate6_factor_scaled10`, `candidate7_factor_scaled5`,
`candidate8_factor_scaled5`, and `candidate9_factor_scaled5` do not close the
localization proof.  The nonzero modular root values came from raw residue lifts,
not from a rational separator that is nonzero at the root.  The unrestricted
exact `QQ` separator job on ML Core (`gtz-sage-local-sep-vnvipz`) also OOM-killed
at about `316.8 GiB`, so the next search should stay modular until a stronger
candidate is isolated.

The first one-step larger separator probe is substantially more promising for
the size-8 clean over-tied case.  Over `F_32003`, `D=5,S=5` gives:

| Case | Matrix | Generator rank | Total rank | Separator rank defect | Kernel candidates | Passing at root |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| size-8 `79656`, omitted `active:012`, `D=5,S=5` | `304099 x 33033` | `29865` | `32802` | `66` | `231` | `231` |

The saved separator projections are not truncated; their union contains `407`
monomials and has row rank `66`, exactly matching the rank defect.  The next
target is therefore multi-prime reconstruction of this rank-66 separator
subspace, not exact `QQ` rank computation.

The corresponding `F_32009` and `F_32027` runs completed with the same rank
data (`separator rank defect = 66`, `231` candidates, all passing at the
refined root).  CRT/RREF reconstruction over `32003,32009,32027` gives a stable
rank-66 rational separator subspace with `407` monomials and the same pivot
pattern over all three primes.  However, evaluating the reconstructed rational
basis at the 1400-bit refined root again gives only residual-scale values: the
largest basis value is about `2.3e-420`, and four basis vectors evaluate to zero
at the printed precision.  Thus `D=5,S=5` for size-8 is also a vanishing
separator space, not the needed open separator.

The analogous size-7 `D=5,S=5` multi-prime computation also stabilized, but
with the same negative interpretation:

| Case | Matrix | Generator rank | Total rank | Separator rank defect | Kernel candidates | Passing at root |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| size-7 `78612`, omitted `active:013`, `D=5,S=5`, over each of `F_32003,F_32009,F_32027` | `306649 x 33033` | `29920` | `32868` | `55` | `132` | `132` |

CRT/RREF reconstruction over these three primes gives a stable rank-55 rational
separator subspace with `407` monomials.  Evaluating the reconstructed basis at
the 1400-bit refined root again gives only residual-scale values: the largest
printed value is about `2.9e-420`, and four basis vectors evaluate to zero at
the printed precision.  Thus the size-7 `D=5,S=5` separator space also appears
to vanish on the local root.

One further rank-only screen for size-8 at `D=5,S=6` over `F_32003` did not
increase the separator defect:

| Case | Matrix | Generator rank | Total rank | Separator rank defect | Kernel dimension |
| --- | ---: | ---: | ---: | ---: | ---: |
| size-8 `79656`, omitted `active:012`, `D=5,S=6`, rank-only | `378585 x 38038` | `29865` | `37807` | `66` | `231` |

This makes a full `D=5,S=6` candidate job low priority: the rank-only screen
shows no larger separator space than `D=5,S=5`.

As a first component/saturation follow-up, the reconstructed vanishing separator
bases were added as extra ideal generators and the omitted determinant was
tested for bounded-degree membership over `F_32003`.  For both clean over-tied
cases, the omitted determinant remains outside the ideal generated by the
selected ten equations plus the full reconstructed `D=5,S=5` separator basis
for multiplier degree `0,1,2,3,4`:

| Case | Extra separator generators | Max multiplier degree tested | Largest matrix | Result |
| --- | ---: | ---: | ---: | --- |
| size-7 `78612`, omitted `active:013` | 55 | 4 | `147220 x 65065` | no membership; rank `19910`, augmented rank `19911` |
| size-8 `79656`, omitted `active:012` | 66 | 4 | `147220 x 76076` | no membership; rank `21857`, augmented rank `21858` |

This does not rule out the component route, but it rules out the easiest
low-degree certificate from the reconstructed separator equations through
multiplier degree 4.  Degree-5 modular screens were launched on ML Core as
`gtz-sage-aug-member-nty05z` (size-7) and
`gtz-sage-aug-member-5rdc2o` (size-8).

A parallel radical-membership check with target `h^2` is also negative through
multiplier degree 3 after adjoining the full reconstructed separator basis:

| Case | Target power | Extra separator generators | Max multiplier degree tested | Largest matrix | Result |
| --- | ---: | ---: | ---: | ---: | --- |
| size-7 `78612`, omitted `active:013` | 2 | 55 | 3 | `56027 x 18590` | no membership; rank `6425`, augmented rank `6426` |
| size-8 `79656`, omitted `active:012` | 2 | 66 | 3 | `56027 x 21736` | no membership; rank `7135`, augmented rank `7136` |

Degree-4 modular screens for this radical-membership check were launched on
ML Core as `gtz-sage-aug-member-kbivq0` (size-7) and
`gtz-sage-aug-member-nukyc7` (size-8).

Thus the remaining proof is more likely a local component separator/saturation
with rational reconstruction across primes, a direct `QQ` local-separator solve,
or an exact RUR/lex computation for the selected local component.

The first local-component diagnostic was run in
`probe_local_component_quotient.py`.  Numerically, at the 1400-bit refined
roots, the selected square subsystem, the selected subsystem plus the full
reconstructed separator basis, and that augmented system plus the omitted
determinant all have Jacobian rank 10:

| Case | Selected rows | Extra separator generators | `rank(J_selected)` | `rank(J_augmented)` | `rank(J_augmented+h)` | Target residual at 500-bit evaluation |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| size-7 `78612`, omitted `active:013` | `0,2,3,4,5,6,7,8,9,10` | 55 | 10 | 10 | 10 | `1.96e-149` |
| size-8 `79656`, omitted `active:012` | `1,2,3,4,5,6,7,8,9,10` | 66 | 10 | 10 | 10 | `0` at 500-bit print precision |

Thus the extra reconstructed equations and the omitted determinant are
first-order compatible with the same smooth isolated local root.  The direct
global quotient probe for the augmented ideal over `F_32003` using Singular
`slimgb` timed out after 180 seconds for both size-7 and size-8, before a
dimension or normal form was produced.  The global Groebner route remains
unattractive.

As a more local denominator search, rank-only `D=6,S=6` modular screens for
`s h = sum_i a_i f_i` were launched on ML Core over `F_32003` as
`gtz-sage-local-sep-cfappr` (size-7) and
`gtz-sage-local-sep-uspefh` (size-8).  The key comparison is whether the
separator rank defect exceeds the previously reconstructed vanishing ranks
55 and 66.

Both rank-only screens completed with real jumps in separator space:

| Case | `D,S` | Matrix | Generator rank | Total rank | Separator rank defect | Kernel dimension |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| size-7 `78612`, omitted `active:013` | `6,6` | `659931 x 88088` | `79399` | `87073` | `334` | `1015` |
| size-8 `79656`, omitted `active:012` | `6,6` | `655873 x 88088` | `79071` | `86680` | `399` | `1408` |

Since `334` is far above the previous rank-55 vanishing space, a full candidate
extraction was launched as `gtz-sage-local-sep-o9fr7w`, saving complete
separator supports for up to 20 passing candidates.  Since `399` is also far
above the previous size-8 rank-66 vanishing space, the analogous size-8 full
candidate extraction was launched as `gtz-sage-local-sep-b321g6`.  Four
multi-prime rank-only stability screens were also launched:

| Case | Prime | ML Core job |
| --- | ---: | --- |
| size-7 `78612` | `32009` | `gtz-sage-local-sep-uskwjl` |
| size-7 `78612` | `32027` | `gtz-sage-local-sep-6nm5vj` |
| size-8 `79656` | `32009` | `gtz-sage-local-sep-0l3ssu` |
| size-8 `79656` | `32027` | `gtz-sage-local-sep-lpvaxl` |

All four stability screens have completed and exactly match the `F_32003`
rank-only results:

| Case | Prime | Matrix | Generator rank | Total rank | Separator rank defect | Kernel dimension |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| size-7 `78612`, omitted `active:013` | `32009` | `659931 x 88088` | `79399` | `87073` | `334` | `1015` |
| size-7 `78612`, omitted `active:013` | `32027` | `659931 x 88088` | `79399` | `87073` | `334` | `1015` |
| size-8 `79656`, omitted `active:012` | `32009` | `655873 x 88088` | `79071` | `86680` | `399` | `1408` |
| size-8 `79656`, omitted `active:012` | `32027` | `655873 x 88088` | `79071` | `86680` | `399` | `1408` |

Both full `F_32003` candidate extractions then completed and found separators
nonzero at the refined root:

| Case | ML Core job | Candidate count | Passing at root | First passing basis index | First `|s(root)|` |
| --- | --- | ---: | ---: | ---: | ---: |
| size-7 `78612`, omitted `active:013`, `D=6,S=6` | `gtz-sage-local-sep-o9fr7w` | `834` | `832` | `0` | `56893.732694294523888969643558830655805705414653933790915097` |
| size-8 `79656`, omitted `active:012`, `D=6,S=6` | `gtz-sage-local-sep-b321g6` | `1408` | `1408` | `0` | `17127.824611337000254977154717255327572103414811313214978664` |

The extracted JSON files are:

- `code/sage/out/local_separator_s7_78612_omit_active1_D6_S6_32003.json`
- `code/sage/out/local_separator_s8_79656_omit_active0_D6_S6_32003.json`

This is modular evidence for actual local-component denominators, not yet a
characteristic-zero certificate.  The full extractions currently save complete
supports for selected `F_32003` passing candidates only; the other primes are
rank-only stability checks.  The next step is support-restricted extraction for
one or more of those passing supports over `F_32009` and `F_32027`, then
CRT/rational reconstruction, exact QQ verification of `s h = sum_i a_i f_i`,
and a numerical/algebraic check that `s` is nonzero on the isolated local
component.

The full-extraction helper now also records the largest evaluated
`|s(root)|` among all tested kernel-basis separators, so a no-passing-candidate
result is easier to diagnose on reruns.

An 8-hour support-restricted follow-up was then run from
`code/mlcore/watch/d6_support`.  The support restriction used two saved
`F_32003` candidates per case: size-7 support indices `19` and `8`, and size-8
support indices `6` and `8`.  All jobs except size-7 support `8` at `F_32003`
completed within the watcher window; the missing job was still pending for
quota.  The completed support families give consistent rank data across
available primes:

| Case / support | Primes completed | Support rank defect | Passing candidates | Interpretation |
| --- | --- | ---: | ---: | --- |
| size-7 `78612`, support `19` | `32003,32009,32027` | `42` | `153` | stable high-value family; 3-prime CRT modulus too small for rational reconstruction |
| size-7 `78612`, support `8` | `32009,32027` | `2` | `9` | missing `F_32003` support job |
| size-8 `79656`, support `6` | `32003,32009,32027` | `63` | `262` | stable high-value family; 3-prime CRT modulus too small for rational reconstruction |
| size-8 `79656`, support `8` | `32003,32009,32027` | `5` | `20` | reconstructs over `QQ`, but reconstructed rational separators vanish at residual scale |

The size-8 support-8 reconstruction is
`code/sage/out/reconstruct_s8_79656_D6_S6_support8_p32003_32009_32027_p1400.json`;
its reconstructed basis evaluates at about `1e-150` or exactly zero at the
refined root, so it is another vanishing rational separator space.  For the
two high-value supports (`s7` support `19`, `s8` support `6`), the modular pivot
patterns agree but rational reconstruction fails with modulus
`32003*32009*32027`, so the next step is simply more primes.  A second watcher
batch was launched in `code/mlcore/watch/d6_more_primes` over
`32029,32051,32057` for these two supports.

As of `2026-08-10T09:20Z`, the second watcher has extracted
`code/sage/out/local_separator_s7_78612_omit_active1_D6_S6_32029_support19.json`.
This fourth prime preserves the same 10-row saved-candidate pivot pattern, but
the four-prime reconstruction still fails:

```text
rational reconstruction failed for basis row 0, column 8, monomial z5^2*q^2,
residues=[15929, 10023, 19529, 16951],
modulus=1050805741917077141
```

The failure is still compatible with a stable rational subspace of larger
height.  However, all current support-restricted jobs saved only the first 10
passing candidates despite support rank defects `42` and `63`.  A follow-up
batch was therefore queued in `code/mlcore/watch/d6_fullsave`, using
`CANDIDATE_LIMIT=0` and distinct `_allcand` output filenames for the same two
support families over primes `32003,32009,32027,32029,32051,32057`.  A local
reconstruction loop in `code/mlcore/watch/d6_followup` runs for 8.5 hours and
automatically attempts limited and all-candidate reconstructions as enough
artifacts appear.

The follow-up loop was corrected to reconstruct from any available prime set,
not only an initial prefix of the requested prime list.  This matters because
`32057` started before `32051`.  At `2026-08-10T09:51Z`, the size-8 support-6
job over `F_32057` completed with the same rank profile as the earlier primes:
matrix rank `79327`, generator rank `79071`, support rank defect `63`, and
`262/262` passing candidates.  The four-prime limited reconstruction over
`32003,32009,32027,32057` still fails with stable pivots:

```text
rational reconstruction failed for basis row 0, column 13, monomial z5^2*z8^2,
residues=[10358, 16769, 26082, 9345],
modulus=1051724364439593553
```

Thus both high-value limited families now have stable four-prime row spaces but
no rational reconstruction at roughly `1e18` CRT modulus.  This is not negative
evidence against the route; it mainly says the saved 10-row slice is not yet a
small-height rational certificate.  The `_allcand` batch remains the more
diagnostic next test.

At `2026-08-10T10:12Z`, the size-7 support-19 job over `F_32057` also
completed with the same rank profile as the earlier support-19 primes: matrix
rank `79629`, generator rank `79399`, support rank defect `42`, and `153/153`
passing candidates.  The five-prime limited reconstruction over
`32003,32009,32027,32029,32057` still fails with the same stable pivot pattern:

```text
rational reconstruction failed for basis row 0, column 7, monomial z5^2*q,
residues=[9350, 24626, 13060, 8937, 6777],
modulus=33685679668635741909037
```

So adding one more prime raised the CRT modulus by four orders of magnitude
without producing a small-height rational basis for the saved 10-dimensional
slice.  The current priority is therefore the queued all-candidate support
outputs, not more limited-candidate repetitions.

At `2026-08-10T11:04Z`, the size-8 support-6 job over `F_32051` completed, also
with matrix rank `79327`, generator rank `79071`, support rank defect `63`, and
`262/262` passing candidates.  The five-prime limited reconstruction over
`32003,32009,32027,32051,32057` still has stable pivots but fails:

```text
rational reconstruction failed for basis row 0, column 21,
monomial z5^2*z6^2*z7^2,
residues=[24832, 28028, 23262, 16062, 7894],
modulus=33708817604653412967203
```

Thus both high-value limited supports have now failed rational reconstruction
at the five-prime level.  The next useful computation is still the full
all-candidate support row space, because the modular rank defects are `42` and
`63` while the limited artifacts save only a 10-row slice.

At `2026-08-10T11:45Z`, the last queued size-8 support-6 limited prime,
`F_32029`, completed.  It again had matrix rank `79327`, generator rank
`79071`, support rank defect `63`, and `262/262` passing candidates.  The
six-prime limited reconstruction over
`32003,32009,32027,32029,32051,32057` still fails:

```text
rational reconstruction failed for basis row 0, column 13, monomial z5^2*z8^2,
residues=[10358, 16769, 26082, 4967, 23057, 9345],
modulus=1079659719059444163926544887
```

This is strong evidence that the saved 10-row size-8 slice is not a small-height
rational certificate.  It does not contradict the modular rank evidence,
because the support-restricted separator quotient has defect `63` and the
limited artifact intentionally saves only 10 passing rows.

At `2026-08-10T12:37Z`, the last queued size-7 support-19 limited prime,
`F_32051`, completed.  It again had matrix rank `79629`, generator rank
`79399`, support rank defect `42`, and `153/153` passing candidates.  The
six-prime limited reconstruction over
`32003,32009,32027,32029,32051,32057` still has stable pivots but fails:

```text
rational reconstruction failed for basis row 0, column 16, monomial z5^4,
residues=[16145, 14063, 17492, 2500, 27480, 12810],
modulus=1079659719059444163926544887
```

Thus both high-value limited supports have now failed rational reconstruction
over all six extra primes.  The limited rank stability remains useful, but the
decisive test has shifted to the full `_allcand` support row spaces.

At `2026-08-10T13:12Z`, the first full `_allcand` support extraction completed:
size-8 `79656`, omitted `active:012`, support `6`, over `F_32003`.  It matches
the limited rank profile but saves the full passing set:

| Case | Prime | Matrix rank | Generator rank | Support defect | Kernel | Passing candidates | Max `|s(root)|` |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| size-8 support `6` `_allcand` | `32003` | `79327` | `79071` | `63` | `1072` | `262/262` | `335532.38032476428718769369276921539404366250403911816963332` |

This confirms that the full-save path is producing useful complete separator
row-space data.  It is not yet reconstructible over `QQ` because only one
all-candidate prime is available for this support.

At `2026-08-10T21:17Z`, the size-7 support-19 `_allcand` family had completed
over all six primes `32003,32009,32027,32029,32051,32057`.  Each prime has the
same rank profile: matrix rank `79629`, generator rank `79399`, support defect
`42`, kernel dimension `723`, and `153/153` passing candidates.  The automatic
six-prime RREF reconstruction succeeded:

```text
code/sage/out/reconstruct_s7_78612_D6_S6_support19_allcand_32003_32009_32027_32029_32051_32057_p1400.json
```

The reconstructed rational row space has rank `42` on `272` monomials.  Its
basis evaluates at the refined local root only at numerical residual scale
(`~1e-149` in the 500-bit follow-up run), so the reconstructed rational
separators appear to vanish on the local component.

Because RREF coordinates can have unnecessarily large height, the same
six-prime row space was also processed by the lattice/LLL recovery helper
`recover_short_separator_lattice.py`.  The first `42` LLL rows give low-height
integer separators; for example the first row has `66` terms, coefficient bound
`10`, and squared norm `768`.  At 1400-bit evaluation all of these short rows
again vanish at residual scale (`~1e-419`).  The next LLL row jumps to modular
artifact size, with coefficient bound about `1.1e26`, and has a huge nonzero
root value.  This cleanly separates the genuine rational row space from the
ambient modular lattice artifacts.

Thus the earlier modular "nonzero at root" signal for this support does not
lift to a low-height rational open denominator.  It comes from residue-basis
choices outside the short rational separator space.  This does not refute the
local-component route, but size-7 support `19` is no longer the promising
certificate support.

As a component/saturation check, the `42` short vanishing separators were
adjoined as extra generators and the omitted determinant was tested over
`F_32003`.  Bounded membership remains negative through multiplier degree `4`;
the degree-4 matrix has `183780` equations and `52052` unknowns, with rank
`37452` and augmented rank `37453`.  A local exact `QQ` degree-6 membership
attempt for the first short separator was launched under a `900s` timeout, but
the process ended without producing an output JSON, so there is no exact
characteristic-zero certificate from that check.

An alternate size-7 support-0 `_allcand` three-prime test was then run over
`32003,32009,32027`.  It is stable but has the same negative interpretation:
each prime has matrix rank `79580`, generator rank `79399`, support defect
`35`, kernel dimension `716`, and `142/142` passing candidates.  Coordinate
RREF reconstruction fails at the three-prime modulus, but LLL recovery cleanly
finds the short rational-looking part:

```text
code/sage/out/short_s7_78612_D6_S6_support0_allcand_p32003_32009_32027_p1400.json
```

The first `35` LLL rows are low-height and vanish at 1400-bit residual scale.
The first row again has `66` terms, coefficient bound `10`, and squared norm
`768`; the last short row has coefficient bound `590`.  Row `35` jumps to
coefficient bound about `1.9e12` and has a huge nonzero root value.  Thus
support `0`, like support `19`, does not expose a plausible nonzero rational
denominator.

The two vanishing spaces were also combined as extra generators: the first `42`
support-19 short rows plus the first `35` support-0 short rows.  Over
`F_32003`, the omitted determinant is still not in the bounded ideal through
multiplier degree `4`:

| Extra generators | Degree | Matrix | Rank | Augmented rank | Result |
| ---: | ---: | ---: | ---: | ---: | --- |
| `77` | `0` | `393 x 87` | `52` | `53` | no membership |
| `77` | `1` | `4149 x 957` | `563` | `564` | no membership |
| `77` | `2` | `21426 x 5742` | `3163` | `3164` | no membership |
| `77` | `3` | `72625 x 24882` | `12451` | `12452` | no membership |
| `77` | `4` | `187867 x 87087` | `38964` | `38965` | no membership |

The delayed size-7 support-8 follow-up also completed over
`32003,32009,32027`.  The support-restricted row space has rank `2` on `34`
monomials, with pivot monomials `z8^3` and `z5^2*z8`.  The reconstructed
rational basis again vanishes at the 1400-bit refined root, with basis values
around `1e-420`.

Adding these two support-8 generators to the previous `42+35` support-19 and
support-0 generators gives `79` extra generators.  The omitted determinant is
still not in the bounded ideal over `F_32003` through multiplier degree `4`:

| Extra generators | Degree | Matrix | Rank | Augmented rank | Result |
| ---: | ---: | ---: | ---: | ---: | --- |
| `79` | `0` | `427 x 89` | `54` | `55` | no membership |
| `79` | `1` | `4430 x 979` | `581` | `582` | no membership |
| `79` | `2` | `22395 x 5874` | `3241` | `3242` | no membership |
| `79` | `3` | `74349 x 25454` | `12673` | `12674` | no membership |
| `79` | `4` | `189785 x 89089` | `39443` | `39444` | no membership |

This makes the low-degree saturation/component certificate unlikely for the
currently tested size-7 support families.  The search should either wait for
the pending size-8 support-6 `_allcand` primes, try a qualitatively different
support family, or move to a more intrinsically local exact representation.

As a direct exact-representation probe, `probe_refined_root_algdep.py` was added
to rerun stronger `algdep` screens on high-precision refined roots.  At 1400
bits, size-7 `z0` appeared to have a degree-9 relation, and size-8 `z2` appeared
to have a degree-22 relation.  Both were unstable under higher-precision
validation: after refining to 2200 bits, the size-7 `z0` best candidates moved
to degrees `14--16`, and after refining size-8 to 2200 and 3000 bits, the
size-8 `z2` best candidates moved from degree `22` to degree `19` and then to
degree `30`.  Thus simple algebraic recognition of a coordinate is currently
overfitting the available precision rather than exposing a usable low-degree
number field.

The full `_allcand` watcher is still useful for the size-8 support-6 family.
As of the same poll, only the `F_32003` size-8 `_allcand` artifact is local; the
remaining size-8 `_allcand` jobs are still queued in ML Core with
`WAITING_RESOURCES_IN_QUOTA`.  The next informative computation is the same
LLL/RREF recovery on size-8 support `6` once at least several more full
all-candidate primes arrive.

To reduce waiting on the non-preemptible queue, preemptible ML Core clones were
submitted for all five missing size-8 support-6 `_allcand` primes
`32009,32027,32029,32051,32057`.  The documented ML Core route is
`/jobs/extra_parameters/preemptible`: the jobs use CLI
`--preemption=allowed`, request the same `20cpu-320ram` flavor, and started
immediately with `qos: Preemptible`.  They are monitored in
`code/mlcore/watch/d6_preempt_s8/latest_summary.md`.

The preemptible and original watcher paths completed the size-8 support-6
`_allcand` family over all six primes `32003,32009,32027,32029,32051,32057`.
Every prime has the same rank profile: matrix rank `79327`, generator rank
`79071`, support defect `63`, kernel dimension `1072`, and `262/262` passing
candidates.  Six-prime RREF reconstruction now succeeds:

```text
code/sage/out/reconstruct_s8_79656_D6_S6_support6_allcand_32003_32009_32027_32029_32051_32057_p1400.json
```

The reconstructed rational row space has rank `63` on `319` monomials.  Its
RREF basis evaluates at the refined root only at residual scale (`~1e-148` to
`~1e-150` in the 500-bit evaluation).

The automatic LLL recovery also succeeded:

```text
code/sage/out/short_s8_79656_D6_S6_support6_allcand_p32003_32009_32027_32029_32051_32057_p1400.json
```

The first `63` short rows are low-height and vanish at the 1400-bit root
residual scale; the largest absolute value among those rows is about
`3.9e-419`.  The next row jumps to modular-artifact scale, with coefficient
bound about `2.6e25` and a huge nonzero root value.  Thus size-8 support `6`,
like the tested size-7 supports, does not expose a plausible nonzero rational
open denominator.

Even after switching Singular to the well-conditioned all-ties square bases,
`slimgb` over `F_32003` still timed out after 240 seconds for the clean size-7
`78612` and size-8 `79656` cases.

## Real section probes

The helper `sample_real_determinant_sections.py` numerically solves real
determinant equations plus three random affine sections through a known
extremal, then tests the actual GTZ inequalities:

- active triples: `lambda_min(P_TT) >= 1/6`;
- inactive triples: `lambda_min(P_TT) <= 1/6`.

This is not a proof, but it directly probes the next semialgebraic layer.

For `P(S(e,e,e),e,e,e)`, four independent section seeds found 32 real section
roots in total.  Exactly four roots passed the inequalities, one per section,
and each was the known base point.  Every other sampled determinant root already
failed active PSD; their inactive maximum was about `1/6 + 2/3 = 5/6`.

The same section-through-known-point probe was run on all seven known
configurations: the six TTSP cases and the stored out-of-family projector
`verify/data/P514_seventh.npy`.  Across the 11 recorded section runs, 51 real
determinant roots were accepted by the nonlinear solver; exactly 11 passed the
GTZ inequalities, again one per section run and always the known point.

Representative output files:

| File | Active size | Accepted roots | Roots passing inequalities |
| --- | ---: | ---: | ---: |
| `code/sage/out/real_section_known_base_seed101.json` | 10 | 11 | 1 |
| `code/sage/out/real_section_known_10b_seed203.json` | 10 | 7 | 1 |
| `code/sage/out/real_section_known_12a_seed201.json` | 12 | 1 | 1 |
| `code/sage/out/real_section_known_12b_seed202.json` | 12 | 1 | 1 |
| `code/sage/out/real_section_known_12c_seed204.json` | 12 | 1 | 1 |
| `code/sage/out/real_section_known_12d_seed205.json` | 12 | 1 | 1 |
| `code/sage/out/real_section_known_seventh_seed206.json` | 13 | 1 | 1 |

## Negative probe

The four-active prefix `0,1,2,3` did not finish a Sage dimension computation
within a 180-second local timeout.

The four-active prefix remains dimension `6` in modular section probes even
after adding `u0*d-1=0`; this is an interior dependency, not just a boundary
artifact.

The known-base saturated determinant system with Singular Groebner/dimension
enabled did not finish within a 120-second local timeout:

```bash
timeout 120s ~/miniforge3/bin/mamba run -n sage Singular -q \
  code/sage/out/known_base_det_compute.sing
```

This does not imply the route is blocked; it means naive exact Groebner
computation on the full known-base determinant system is not a cheap local
smoke test.

## Follow-up known-orbit chart probes

After closing the base active-set chart, the same saturated determinant probe
was run on the remaining canonical active-set masks from the known equality
configurations.  Over `F_32003`:

| Canonical mask | Active size | Representative | Dimension | Degree |
| ---: | ---: | --- | ---: | ---: |
| `32243` | 12 | `P(S(e,e),S(e,e),e,e)` orbit | `1` | `16` |
| `31215` | 12 | `P(S(P(e,e),P(e,e)),S(e,e))` orbit | `1` | `16` |
| `32253` | 13 | out-of-family `(5/14,9/14)` point | `0` | `32` |

Artifacts:

```text
code/sage/out/probe_det_mask32243_p32003.json
code/sage/out/probe_det_mask31215_p32003.json
code/sage/out/probe_det_mask32253_p32003.json
```

The size-13 mask is already finite in the standard chart.  The exact
characteristic-zero lex run

```text
code/sage/out/lex_det_mask32253_QQmod_z0.json
```

has quotient dimension `32` and basis

```text
8*z0^2-5, 63*u0-8, 8*z8^2-5, 64*z7^2-45,
64*z6^2-45, 8*z5^2-5, 5*z4-8*z5*z7*z8,
5*z3-8*z5*z6*z8, z2, 45*z1-64*z6*z7*z0.
```

Thus it is a 32-point real sign family.  Numerical substitution from the exact
formulas gives 32 distinct projectors, all passing the GTZ inequalities with
`F-1/6` at roundoff scale and all having actual active set

```text
[0,2,3,4,5,6,7,8,10,11,12,13,14].
```

The next symbolic targets are therefore the two size-12 curve ideals.  They
have the same modular degree (`16`) as the intermediate branches in the base
CAD cascade, so they are plausible candidates for the same relation/forced
coordinate strategy.

A first arbitrary one-section of mask `32243` used the zero-sum linear form

```text
-5*z0 + 17*z1 - 15*z2 - 11*z3 + 12*z4 + 7*z5 + 14*z6 - 17*z7 - 2*z8.
```

Over each of `32003,32009,32027,32029,32051,32057` it has lex quotient
dimension `16`; over `F_32003` the basis is shape-like with `z0=0`, `z1=0`, and
a degree-16 eliminant in `z8`.  However, six-prime rational reconstruction of
that eliminant still fails for two coefficients, and the direct exact
`QQ/modslimgb` section timed out after `300s`:

```text
code/sage/out/lex_det_mask32243_p32003_sec601_z0.json
code/sage/out/lex_det_mask32243_sec601_z8_eliminant_crt_6primes.json
code/sage/out/lex_det_mask32243_QQmod_sec601_z0.json
```

This suggests that arbitrary sections are not the cleanest route for the
size-12 curves.  The next attempt should use sections or normal forms adapted
to the known point, or screen modular slices for low-height minor relations
before trying exact reconstruction.

### Exact real branch split for the size-12 curves

The arbitrary-section obstruction was bypassed by using the inverse Gram
determinant `u0`.  For mask `32243`, an exact `QQ/slimgb` Groebner basis of the
full determinant ideal was computed:

```text
code/sage/out/gb_det_mask32243_QQ.json
code/sage/out/gb_det_mask32243_QQ_basis.txt
```

The basis contains

```text
729*u0^3 + 81*u0^2 - 9*u0 - 1
  = 729*(u0 - 1/9)*(u0 + 1/9)^2.
```

For mask `31215`, exact `QQ/modslimgb` was needed; the sandboxed run hits
Singular's local socket restriction, while the approved unsandboxed run
completed:

```text
code/sage/out/gb_det_mask31215_QQmod.json
code/sage/out/gb_det_mask31215_QQmod_basis.txt
```

The basis contains the monic equivalent

```text
u0^3 + 1/9*u0^2 - 1/81*u0 - 1/729
  = (u0 - 1/9)*(u0 + 1/9)^2.
```

Since `u0*d(Z)-1=0` and `d(Z)=det(Y^T Y)>0` on the real full-rank chart, all
real points must lie on the branch `u0=1/9`; the `u0=-1/9` curve branch has no
real chart points.

Adding `9*u0-1` gives exact finite lex bases.

For mask `32243`:

```text
code/sage/out/lex_det_mask32243_QQmod_u_pos_z0.json
```

with basis

```text
z0^2-1, 9*u0-1, z8, 8*z7^2-5, 8*z6^2-5, z5,
8*z4^2-5, 5*z3-8*z4*z6*z7, z2^2-1, 5*z1-8*z6*z7*z0.
```

For mask `31215`:

```text
code/sage/out/lex_det_mask31215_QQmod_u_pos_z2.json
```

with basis

```text
z2^2-1, 9*u0-1, z8^2-1, 8*z7^2-5, 8*z6^2-5, z5^2-1,
z4-z5*z7*z8, z3-z5*z6*z8, z1, z0.
```

Both branches have lex quotient dimension `32`.  Explicit enumeration of the
five independent signs gives 32 distinct projectors in each case.  All pass the
active and inactive GTZ inequalities with `F-1/6` at roundoff scale, and in each
case the actual active set is exactly the selected size-12 mask:

```text
code/sage/out/classify_det_mask32243_u_pos_sign_family.json
code/sage/out/classify_det_mask31215_u_pos_sign_family.json
```

Auxiliary modular decomposition checks also support the interpretation of the
discarded branch.  For `32243`, `facstd` over `F_32003` and `F_32009` splits the
positive-dimensional locus into eight degree-2 components, all sitting over
`u0=-1/9`:

```text
code/sage/out/fac_dump_mask32243_p32003.json
code/sage/out/fac_dump_mask32243_p32009.json
```

## Practical next target

The base ten-active chart is now accounted for by the CAD cascade plus finite
residue split; the two size-12 known-orbit charts are finite on their only real
branch; and the thirteen-active out-of-family chart is already finite.  In the
standard chart, this accounts for all canonical active-set orbits arising from
the seven known equality projectors.

The next certifiable targets are no longer the known equality active sets, but
the remaining active-set orbits that could support a sub-threshold Clarke
critical point.  The practical route is to automate the same branch-polynomial
screen used above:

1. For each symmetry-canonical active set of plausible size, compute modular
   determinant ideals with `d` inverted.
2. Search exact or CRT-reconstructed univariate relations in `u0=1/d`.
3. Discard real-impossible branches (`u0 <= 0`) and lexify only the positive
   branches.
4. On any remaining positive-dimensional real branch, return to the KKT or
   active-minor PSD constraints.

The relaxed-bound idea remains useful only as a companion gap certificate:
prove a computable bound `F >= 1/6 - eps` away from the equality neighborhoods,
then use the exact active-set/CAD machinery in the remaining near-threshold
region.  A relaxed bound by itself does not imply the sharp GTZ statement.

## Low-active follow-up after the known charts

The direct determinant-ideal strategy becomes much harder on low-active
survivors.  The first no-section modular probes over `F_32003`, with `d`
inverted, timed out before a dimension was available:

| Active-set mask | Active size | Linear sections | Timeout | Result |
| ---: | ---: | ---: | ---: | --- |
| `78593` | 6 | 0 | 180s | timeout |
| `4040` | 7 | 0 | 180s | timeout |
| `8065` | 7 | 0 | 180s | timeout |
| `13105` | 7 | 0 | 180s | timeout |

Adding enough zero-sum affine sections to make the expected dimension zero helps
only for the unique size-six survivor.  The same three-section system has stable
degree `2040` over two primes:

| Active-set mask | Prime | Sections | Dimension | Degree | Basis size | Elapsed |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `78593` | `32003` | 3 | 0 | 2040 | 1343 | 150s |
| `78593` | `32009` | 3 | 0 | 2040 | 1343 | 141s |
| `4040` | `32003` | 2 | ? | ? | ? | timeout, 180s |
| `8065` | `32003` | 2 | ? | ? | ? | timeout, 180s |
| `13105` | `32003` | 2 | ? | ? | ? | timeout, 180s |

The size-six section degree explains why a naive lex/RUR follow-up is likely too
large: it is over two orders of magnitude bigger than the known-base degree-18
section.

The real-section sampler was extended so it can take a saved `v30` margin point
as the section center.  Two independent random three-sections through the
size-six margin minimizer find many real determinant roots, but none satisfying
both active PSD and inactive `<=1/6` inequalities:

| File | Accepted roots | Passing roots | Minimum inactive overshoot |
| --- | ---: | ---: | ---: |
| `code/sage/out/real_section_low_s6_78593_v30_seed801.json` | 79 | 0 | `0.030152235779059622` |
| `code/sage/out/real_section_low_s6_78593_v30_seed802.json` | 70 | 0 | `0.030226783505172466` |

This supports the same picture as the margin minimisation: the size-six
low-active equality locus is real and large, but sampled real points force
inactive growth by about `3e-2` or more.  It is evidence for a positive-margin
certificate, not a proof.

The clean size-seven active-plus-tie square determinant systems also timed out
over `F_32003` before producing a dimension:

| Active-set orbit | Outside ties | Timeout | Result |
| --- | --- | ---: | --- |
| size-7 canon `13105` | `[9,10,18]` | 180s | timeout |
| size-7 canon `78595` | `[4,13,19]` | 180s | timeout |
| size-7 canon `78601` | `[4,7,14]` | 180s | timeout |

Together with the earlier failures for `14113`, `78612`, and the over-tied
size-eight cases, this makes plain Groebner/RUR on low-active margin systems a
poor next target.  The more realistic exact layer is local: Krawczyk-certified
roots plus local component/separator certificates, or a different parametrised
representation of the relevant equality component.

The interval driver `verify/v19_interval_bounds.py` now accepts `--threshold`.
A calibration at the relaxed target `0.15 = 1/6 - 1/60`:

```text
verify/out/v19_relaxed_015_cascade_1000.json
```

processed 1000 boxes over all 20 charts with the cascade affine bound.  It
certified 396 boxes at the relaxed threshold, of which 368 also cleared the
sharp `1/6` threshold, and left a queue of 228 boxes.  Lowering the target thus
does not by itself remove the interval bottleneck; the current chart cover and
inverse enclosure remain the limiting factor.

## Structured size-six ansatz

The singular over-tied size-six root for mask `78593`, ties
`[5,7,12,14,18]`, has the numerical chart pattern

```text
Z = [ -a   a  -b
       c   d  -a
      -d  -c  -a ].
```

Substituting this four-parameter ansatz into the six active determinant
equations and five outside-tie determinant equations gives a much smaller exact
system.  Over both `F_32003` and `QQ`, the unrestricted ansatz ideal is
zero-dimensional of degree `280`:

```text
code/sage/out/s6_ansatz_p32003.json
code/sage/out/s6_ansatz_QQ.json
```

The exact `QQ` grevlex basis has the tail relations

```text
a^2((c+d)^2-5) = 0,
((c+d)^2-5)((c-d)^2-5) = 0.
```

Thus every real point lies on either `a=0` or `(c+d)^2=5`.  The `a=0` branch
also finishes exactly:

```text
code/sage/out/s6_ansatz_azero_QQ.json
code/sage/out/s6_ansatz_azero_fglm_QQ.json
```

It has quotient degree `96`; FGLM gives a six-polynomial lex basis, and exact
enumeration over Sage's algebraic real field finds no real points:

```text
code/sage/out/s6_ansatz_azero_sign_cert_QQ.json
```

As a redundant check, the `(c-d)^2=5` branch also has quotient degree `96`, its
grevlex basis contains `a^2`, and its lex real enumeration is empty:

```text
code/sage/out/s6_ansatz_cdiff_QQ.json
code/sage/out/s6_ansatz_cdiff_fglm_QQ.json
code/sage/out/s6_ansatz_cdiff_sign_cert_QQ.json
```

It remains to inspect the `(c+d)^2=5` branch, which contains the refined
over-tied numerical root.  On this branch, the exact `QQ` grevlex computation
also finishes:

```text
code/sage/out/s6_ansatz_cplusd_QQ.json
```

It has quotient degree `184` and basis size `51`, matching the modular
`F_32003` grevlex run.  Lex computations over six primes all have the same
shape: quotient degree `184`, basis size `10`, and a degree-22 eliminant in
`q`:

```text
code/sage/out/s6_ansatz_cplusd_lex_p32003.json
code/sage/out/s6_ansatz_cplusd_lex_p32009.json
code/sage/out/s6_ansatz_cplusd_lex_p32027.json
code/sage/out/s6_ansatz_cplusd_lex_p32029.json
code/sage/out/s6_ansatz_cplusd_lex_p32051.json
code/sage/out/s6_ansatz_cplusd_lex_p32057.json
```

After dividing the visible nuisance factor

```text
(q-3)(q-5)^2(q^2-3q+6)^2(q^2-9q+24)^2
```

the six-prime CRT reconstruction gives a rational lift:

```text
code/sage/out/s6_ansatz_cplusd_q_residual_crt_6primes.json
```

The residual degree-11 polynomial is

```text
q^11 - 575/16*q^10 + 562849/1024*q^9 - 4768527/1024*q^8
+ 97733649/4096*q^7 - 311529979/4096*q^6
+ 154321449/1024*q^5 - 374696325/2048*q^4
+ 543960225/4096*q^3 - 225237375/4096*q^2
+ 24046875/2048*q - 253125/256
```

and factors as

```text
(q^3 - 9*q^2 + 81/4*q - 45/4)
*
(q^8 - 431/16*q^7 + 293857/1024*q^6 - 776859/512*q^5
 + 2094513/512*q^4 - 1353145/256*q^3
 + 3301125/1024*q^2 - 453375/512*q + 5625/64).
```

The target over-tied root has
`q = 1.180695836766579...`, the largest real root of the octic, giving inactive
overshoot `(q-1)/6 = 0.0301159727944...`.

The CRT lift has now been checked over characteristic zero.  The helper
`verify_s6_ansatz_eliminant.py` reduces the full degree-22 product

```text
(q-3)(q-5)^2(q^2-3q+6)^2(q^2-9q+24)^2 * residual(q)
```

against the exact `QQ` grevlex basis for the `(c+d)^2=5` branch.  The normal
form is zero:

```text
code/sage/out/s6_ansatz_cplusd_eliminant_verify_QQ.json
```

Thus this degree-22 polynomial is an exact element of the branch elimination
ideal over `QQ`; the remaining issue on this branch is no longer the modular
reconstruction, but the real-root sign certification.

The real-root sign certification is now also exact on this branch.  Starting
from the exact `QQ` grevlex branch basis, the helper
`lexify_s6_ansatz_basis.py` converts to lex order by FGLM without recomputing
the determinant ideal.  The full branch conversion has lex quotient dimension
`184`:

```text
code/sage/out/s6_ansatz_cplusd_fglm_QQ.json
code/sage/out/s6_ansatz_cplusd_fglm_QQ_lex_basis.txt
```

Splitting by residual factors gives small exact lex systems:

| Branch | Lex dimension | Real interpretation |
| --- | ---: | --- |
| cubic `q^3-9q^2+81/4 q-45/4` | `24` | no real points, since the lex basis contains `b^2+1` and `a^2` |
| octic residual factor | `32` | `16` real points |
| nuisance `q=3` | `32` | `8` real points |
| nuisance `q=5` | `32` | `16` real points |

The exact lex bases are:

```text
code/sage/out/s6_ansatz_cplusd_cubic_fglm_QQ.json
code/sage/out/s6_ansatz_cplusd_octic_fglm_QQ.json
code/sage/out/s6_ansatz_cplusd_q3_fglm_QQ.json
code/sage/out/s6_ansatz_cplusd_q5_fglm_QQ.json
```

The helper `certify_s6_ansatz_signs.py` enumerates the real algebraic points
from these lex bases in Sage's algebraic real field and checks principal-minor
certificates.  Its output is:

```text
code/sage/out/s6_ansatz_cplusd_sign_cert_QQ.json
code/sage/out/s6_ansatz_structured_sign_cert_QQ.json
```

It certifies all `40` real points in the octic, `q=3`, and `q=5` branches as
infeasible for the GTZ inequalities:

| Branch | Real points | Active-PSD failure | Inactive positive-definite failure | Unclassified |
| --- | ---: | ---: | ---: | ---: |
| octic residual factor | `16` | `12` | `4` | `0` |
| `q=3` | `8` | `8` | `0` | `0` |
| `q=5` | `16` | `0` | `16` | `0` |
| total | `40` | `20` | `20` | `0` |

For the octic branch, the real `q` values are approximately
`0.249409271817`, `0.342580053484`, `0.605330032909`, and
`1.180695836767`.  The first three are excluded by negative active principal
minors; the last is excluded by an inactive positive-definite block.  The
`q=3` branch is excluded by active minors, and the `q=5` branch by inactive
positive-definite blocks.

A deterministic numerical classifier using the lifted q-values and the ansatz
equations found the same real-root profile with 250 and 800 starts:

| Run | Accepted roots | Active PSD | Inactive feasible | Both |
| --- | ---: | ---: | ---: | ---: |
| `classify_s6_ansatz_cplusd_numeric_s250.json` | 40 | 20 | 0 | 0 |
| `classify_s6_ansatz_cplusd_numeric_s800.json` | 40 | 20 | 0 | 0 |

Per q-value, the active-PSD real roots occur only at `q=1.180695836766579...`
and `q=5`; all fail inactive inequalities.  The smaller real q-values fail
active PSD.  The exact sign certificate above confirms this numerical
classification on the certified `(c+d)^2=5` branch.

Combining the empty `a=0` real branch with the certified `(c+d)^2=5` branch
excludes the full four-parameter structured ansatz for the selected six-active
and five-outside-tie threshold pattern.

## Plucker-pattern probes for the hard size-7 and size-8 over-tied roots

The high-precision over-tied roots

```text
code/sage/out/refine_overtie_s7_78612_t2_7_10_19_p2200.json
code/sage/out/refine_overtie_s8_79656_t10_11_15_p3000.json
```

were screened for small integer relations by
`screen_refined_root_relations.py`.  No affine linear relations among
`q,z0,...,z8` were found up to four terms and coefficient height `10^6`.
However, both roots satisfy five height-one quadratic monomial relations at the
precision floor:

```text
code/sage/out/relations_overtie_s7_s8_quad3_p1400.json
```

For `s7_78612` these are

```text
z8 - z3*z7 + z4*z6 = 0
1 + z0*z7 - z1*z6 = 0
z2 + z0*z4 - z1*z3 = 0
z3 + z0*z8 - z2*z6 = 0
z4 + z1*z8 - z2*z7 = 0
```

For `s8_79656` they are

```text
z8 + z0*z7 - z1*z6 = 0
z0 + z3*z8 - z5*z6 = 0
1 - z3*z7 + z4*z6 = 0
z1 + z4*z8 - z5*z7 = 0
z5 + z0*z4 - z1*z3 = 0
```

These are naturally interpreted as equalities among raw Plucker coordinates of
the standard `Y=[I;Z]` chart.

Adding these Plucker relations to the active/tie determinant systems gives
finite modular loci.  The first full basis computations over `F_32003` were:

| Case | Output | Dimension | Degree | Basis size |
| --- | --- | ---: | ---: | ---: |
| `s7_78612` | `code/sage/out/plucker_locus_s7_78612_p32003.json` | `0` | `4992` | `2728` |
| `s8_79656` | `code/sage/out/plucker_locus_s8_79656_p32003.json` | `0` | `2832` | `1952` |
| `s8_79656` reduced ansatz | `code/sage/out/plucker_ansatz_s8_79656_p32003.json` | `0` | `2832` | `1901` |

The denominator-cleared `s8` patch probe produced a positive-dimensional
artifact (`dimension=5`) and its large raw output was discarded; this patch
model should not be used as a certificate without saturation/inversion cleanup.

A second-prime check over `F_32009` reproduced the two relevant quotient
degrees without adding the large raw Groebner basis dumps to the repository:

| Case | Formulation | Dimension | Degree |
| --- | --- | ---: | ---: |
| `s7_78612` | full Plucker locus | `0` | `4992` |
| `s8_79656` | reduced Plucker ansatz | `0` | `2832` |

The grevlex bases above do not expose a univariate `q` polynomial directly.
The helper `compute_q_power_relation.py` instead computes the first exact
relation among `1,q,q^2,...` in the quotient algebra.  This gives compact
modular `q` eliminants:

| Case | Prime | Quotient degree | q-relation degree | Residual | Factor summary |
| --- | ---: | ---: | ---: | --- | --- |
| `s7_78612` | `32003` | `4992` | `510` | zero | `44` factors; largest degrees `62,64` |
| `s7_78612` | `32009` | `4992` | `510` | zero | `38` factors; largest degrees `88,112` |
| `s8_79656` | `32003` | `2832` | `326` | zero | `23` factors; largest degrees `92,101` |
| `s8_79656` | `32009` | `2832` | `326` | zero | `27` factors; largest degrees `70,130` |
| `s8_79656` | `32027` | `2832` | `326` | zero | `28` factors; largest degrees `60,114` |
| `s8_79656` | `32029` | `2832` | `326` | zero | `30` factors; largest degrees `77,132` |

The output files are:

```text
code/sage/out/qrel_plucker_locus_s7_78612_p32003.json
code/sage/out/qrel_plucker_locus_s7_78612_p32003_factors.json
code/sage/out/qrel_plucker_locus_s7_78612_p32009.json
code/sage/out/qrel_plucker_locus_s7_78612_p32009_factors.json
code/sage/out/qrel_plucker_ansatz_s8_79656_p32003.json
code/sage/out/qrel_plucker_ansatz_s8_79656_p32003_factors.json
code/sage/out/qrel_plucker_ansatz_s8_79656_p32009.json
code/sage/out/qrel_plucker_ansatz_s8_79656_p32009_factors.json
code/sage/out/qrel_plucker_ansatz_s8_79656_p32027.json
code/sage/out/qrel_plucker_ansatz_s8_79656_p32027_factors.json
code/sage/out/qrel_plucker_ansatz_s8_79656_p32029.json
code/sage/out/qrel_plucker_ansatz_s8_79656_p32029_factors.json
```

The helper `summarize_qrel_crt.py` does coefficientwise CRT lifts of these
modular q-polynomials.  With primes `32003,32009,32027,32029` for `s8_79656`,
the modulus has decimal size about `10^18`, but full-polynomial rational
reconstruction is still incomplete (`201` coefficients reconstruct and `126`
fail).  This suggests that reconstructing the entire degree-326 q-eliminant is
less realistic than isolating stable lower-degree factors/components.

The previous high-precision `algdep` hints were also tested directly against
the exact modular q-eliminants.  The script `screen_algdep_against_qrel.py`
reduces each candidate polynomial modulo the available primes and computes its
gcd with the exact modular q-relation.  For the top twelve `q` candidates from
both hard roots, all gcd degrees were zero: for `s7_78612` over
`32003,32009`, and for `s8_79656` over `32003,32009,32027`.  Thus the apparent
low-degree `algdep` relations do not define components of the exact Plucker
q-eliminants.

The additional-prime batch pushes this conclusion further.  For the reduced
`s8_79656` ansatz, 37 primes through `32341` all reproduce quotient degree
`2832`, q-relation degree `326`, and zero residual.  The CRT snapshot

```text
code/sage/out/qrel_plucker_s8_79656_crt_37primes_32003_32341.json
```

has modulus size `10^166.776`; full-polynomial rational reconstruction still
fails on `108` of `327` coefficients, and the centered integer lift has maximum
coefficient size `10^166.474`.  For the full `s7_78612` Plucker locus, the
26-prime CRT snapshot

```text
code/sage/out/qrel_plucker_s7_78612_crt_26primes_32003_32341.json
```

has modulus size `10^117.197`; it keeps quotient degree `4992`, q-relation
degree `510`, and still has `197` failed coefficient reconstructions.

The top-twelve numerical `algdep` q-candidates were also rechecked against the
same 37-prime s8 snapshot and the 26-prime s7 snapshot:

```text
code/sage/out/algdep_vs_qrel_s8_79656_q_best12_37primes.json
code/sage/out/algdep_vs_qrel_s7_78612_q_best12_26primes.json
```

Every listed gcd degree is zero.  Thus the high-precision `algdep` polynomials
remain numerical mirages, not exact factors of the Plucker q-eliminants.

The new helper `summarize_factor_product_counts.py` gives a cheap
factor-degree pre-screen for this narrower attack.  For the s8 top-twelve
`algdep` degrees, the artifact

```text
code/sage/out/factor_product_counts_s8_79656_algdep_degrees_37primes.json
```

shows that `p32051` has zero factor-product subsets for ten of the twelve
candidate degrees; only degree `43` (11 products) and degree `44` (1 product)
survive at that anchor.  The corresponding s7 artifact

```text
code/sage/out/factor_product_counts_s7_78612_algdep_degrees_26primes.json
```

is fully capped at `20000` products for every tested degree and prime, so this
simple beam strategy is not selective for the full s7 locus.

The helper `recover_q_factor_by_degree.py` then tries the factor-level CRT:
enumerate modular factor products of a fixed degree, keep a beam ordered by
centered coefficient size, and look for a stable small integer lift.  Both
surviving s8 anchored degrees are negative:

```text
code/sage/out/recover_s8_79656_qfactor_deg43_37primes_beam50.json
code/sage/out/recover_s8_79656_qfactor_deg44_37primes_beam50.json
```

For degree `43`, the 37-prime best candidate has maximum coefficient size
`10^166.390`; for degree `44`, the 37-prime best has size `10^166.406`.
Both track the CRT modulus (`10^166.776`) rather than stabilizing.  These are
capped beam screens (`max_candidates_per_prime=5000` at most primes), not
complete no-factor theorems, but they are strong negative evidence for the only
two s8 degrees from the old numerical `algdep` output that passed the `p32051`
factor-product pre-screen.

The reduced s7 Plucker ansatz was then checked directly:

```text
code/sage/out/qrel_plucker_ansatz_s7_78612_p32003.json
code/sage/out/qrel_plucker_ansatz_s7_78612_p32009.json
```

At both primes the ansatz basis has dimension `0`, quotient degree `4992`, and
q-relation degree `510`.  The q-polynomial coefficients agree exactly with the
corresponding full-locus artifacts
`qrel_plucker_locus_s7_78612_p32003.json` and
`qrel_plucker_locus_s7_78612_p32009.json`.  This validates the 8-variable
ansatz as an equivalent representation of the s7 Plucker component for these
modular checks, although basis construction was not faster in this test
(`basis_size=3282`, about 25--31 minutes per prime).

A different quotient-structure screen is in

```text
code/sage/out/linear_form_screen_s8_79656_p32009_deg900.json
```

It measures first power-relation degrees for coordinates and linear forms in
the s8 ansatz quotient.  The simple coordinates have substantially larger
minimal-polynomial degrees than `q`:

| Expression | Relation degree |
| --- | ---: |
| `q` | `326` |
| `e` | `632` |
| `a` | `647` |
| `b` | `671` |
| `c,d` | `678` |
| `x,y` | `682` |
| `q+a`, `q+b` | `677` |
| `q+c`, `q+d` | `686` |
| `q+x`, `q+y` | `690` |
| `q+e` | `656` |

The random-looking linear form
`q+a+2*b+3*c+5*d+7*x+11*y+13*e` has no relation through degree `900`.
This is the first positive indication that a generic primitive element may be
reachable even though the bare `q`-eliminant is too collapsed for direct
component isolation.

Running the same linear form to the full quotient limit gives:

```text
code/sage/out/linear_form_screen_s8_79656_p32009_random_full.json
```

The exact relation degree is `2800`, with zero residual, compared with quotient
degree `2832`.  Thus this linear form is not fully primitive, but it collapses
only a `32`-dimensional part of the quotient.  This is a much sharper
elimination coordinate than `q` and justifies a short random search for a true
degree-`2832` primitive element.

The first polynomial-power random search was too slow at the full quotient
degree, so a CSR/Wiedemann backend was added:

```text
code/sage/screen_random_linear_forms_scipy.py
```

It caches the eight quotient multiplication actions over `F_32009` under
`/tmp/gtz_plucker/plucker_ansatz_s8_79656_p32009_actions`.  The action matrices
have about `4.1e5`--`5.4e5` nonzeros each, and a generic combined linear action
has about `3.06e6` nonzeros.  On the cached action matrices, a 32-trial random
screen gives:

```text
code/sage/out/random_linear_forms_scipy_s8_79656_p32009_seed20260813_32trials.json
```

| Relation degree | Count |
| ---: | ---: |
| `2800` | `31` |
| `2796` | `1` |

Every listed relation has exact residual zero.  No trial reached the full
quotient degree `2832`, so the evidence now points away from a degree-`2832`
primitive element for this nonradical quotient and toward a generic
degree-`2800` radical/separable support plus extra multiplicity structure.

For the first random form

```text
6*a+9*b+18*c+29*d+23*x+15*y+2*e+5*q
```

the stored degree-`2800` relation factors over `F_32009` as:

```text
code/sage/out/random_linear_forms_scipy_s8_79656_p32009_seed20260813_trial0_coeffs.json
code/sage/out/random_linear_forms_scipy_s8_79656_p32009_seed20260813_trial0_factors.json
```

The factorization has `204` factors.  The squarefree degree is `2704`, while
the degree with multiplicities is `2800`; the multiplicity excess is `96`.
Relative to the quotient degree `2832`, the gaps are:

| Quantity | Value |
| --- | ---: |
| quotient degree minus generic relation degree | `32` |
| quotient degree minus squarefree degree | `128` |
| relation degree minus squarefree degree | `96` |

The repeated factors are concentrated in `16` linear factors and `16` quadratic
factors, all with multiplicity `3`.  This is useful new structure: the
finite-field Plucker ansatz quotient is not behaving as `2832` reduced points
under generic linear projection, but the radical support appears much more
compressed than the original quotient degree.

The same screen was repeated at `p=32003`:

```text
code/sage/out/random_linear_forms_scipy_s8_79656_p32003_seed20260813_8trials.json
code/sage/out/random_linear_forms_scipy_s8_79656_p32003_seed20260813_trial0_coeffs.json
code/sage/out/random_linear_forms_scipy_s8_79656_p32003_seed20260813_trial0_factors.json
```

All eight random forms had exact relation degree `2800`.  The first-form
factorization again has degree `2800`, squarefree degree `2704`, and
multiplicity excess `96`; the quotient-minus-squarefree gap is again `128`.
The factor degrees differ over `F_32003`, but the degree/multiplicity invariants
match `F_32009`, making the degree-`2800` generic-projection behavior unlikely
to be a one-prime artifact.

The repeated-support factor product was then extracted at both primes by taking
the squarefree product of all factors with multiplicity at least `2`:

```text
code/sage/out/repeated_support_s8_79656_p32003_linform_seed20260813_trial0.json
code/sage/out/repeated_support_s8_79656_p32009_linform_seed20260813_trial0.json
```

In both characteristics this product has degree `48`.  A two-prime CRT lift is
stored in

```text
code/sage/out/repeated_support_s8_79656_linform_seed20260813_trial0_crt_32003_32009.json
```

The modulus has decimal size about `10^9.010`; rational reconstruction succeeds
for `36` of the `49` coefficients and fails for `13`.  The centered integer
lift has maximum coefficient size about `10^8.687`, so more primes are needed
before interpreting the lift.  This is now a small, concrete CRT target: a
degree-`48` polynomial associated with the multiple support of the generic
linear projection.

A third-prime check at `p=32027` was then computed from a fresh raw ansatz basis
in `/tmp`:

```text
code/sage/out/random_linear_forms_scipy_s8_79656_p32027_seed20260813_trial0_coeffs.json
code/sage/out/random_linear_forms_scipy_s8_79656_p32027_seed20260813_trial0_factors.json
code/sage/out/repeated_support_s8_79656_p32027_linform_seed20260813_trial0.json
code/sage/out/repeated_support_s8_79656_linform_seed20260813_trial0_crt_32003_32009_32027.json
```

The same trial-0 linear form again has relation degree `2800`, squarefree
degree `2704`, multiplicity excess `96`, and repeated-support degree `48`.
The three-prime CRT modulus has decimal size about `10^13.516`; rational
reconstruction succeeds for `33` of `49` coefficients and fails for `16`, with
centered integer coefficients still tracking the modulus.  Thus the
degree-`48` repeated-support polynomial is stable across three primes, but more
prime data is needed before attempting a rational lift.

The same compact pipeline was run at `p=32029`:

```text
code/sage/out/random_linear_forms_scipy_s8_79656_p32029_seed20260813_trial0_coeffs.json
code/sage/out/random_linear_forms_scipy_s8_79656_p32029_seed20260813_trial0_factors.json
code/sage/out/repeated_support_s8_79656_p32029_linform_seed20260813_trial0.json
code/sage/out/repeated_support_s8_79656_linform_seed20260813_trial0_crt_32003_32009_32027_32029.json
```

The invariants remain unchanged: relation degree `2800`, squarefree degree
`2704`, multiplicity excess `96`, repeated-support degree `48`.  The four-prime
CRT modulus has decimal size about `10^18.022`; rational reconstruction still
fails for `19` of `49` coefficients, and the centered integer lift remains near
the modulus.  This suggests that the chosen random linear form is a high-height
coordinate for the degree-`48` repeated support; adding primes is useful, but a
lower-height linear form should also be tested.

The first rectangular span test,

```text
code/sage/out/generator_span_s8_79656_p32009_q_a_325_20_partial.json
```

is negative for the simplest two-generator tower: the span of
`q^i a^j`, `0 <= i <= 325`, plateaued at rank `680` by `j=6` and remained at
rank `680` through `j=12`, far below the full quotient degree `2832`.  Thus
`q` together with powers of the single coordinate `a` does not provide a small
finite-algebra model of the full s8 Plucker quotient.

This is not yet a proof of infeasibility for these Plucker loci.  It is a
substantial compression: the next realistic step is to treat the degree-`2800`
generic linear-form relation as a radical/RUR target for the s8 ansatz, and to
separate or saturate the lower-dimensional multiplicity structure rather than
continuing to hunt blindly for a degree-`2832` primitive element.

### Low-height repeated-support bank

The CSR/Wiedemann screen was extended to accept explicit coefficient banks via
`--coefficients-json`.  A first low-height bank of `16` linear forms was tested
on the s8 ansatz quotient:

```text
code/sage/linear_form_candidates_s8_79656_lowheight_v1.json
code/sage/out/linear_form_lowheight_bank_v1_s8_79656_p32003_coeffs.json
code/sage/out/linear_form_lowheight_bank_v1_s8_79656_p32009_coeffs.json
code/sage/out/linear_form_lowheight_bank_v1_s8_79656_p32027_coeffs.json
code/sage/out/linear_form_lowheight_bank_v1_s8_79656_p32029_coeffs.json
code/sage/out/linear_form_lowheight_bank_v1_s8_79656_p32051_coeffs.json
```

Across `p=32003,32009,32027,32029`, `12` of the `16` candidates have exact
relation degree `2800`; the same `12` remain degree `2800` after adding
`p=32051`.  The persistent excluded forms are the all-ones control, the
structural degree-`2788` form
`2*a-3*b+c+2*d-2*x+3*y+e+q`, and two forms that drop at `p=32029`.

For every eligible candidate, the squarefree product of repeated factors again
has degree `48`.  The ranking script

```text
code/sage/rank_repeated_support_candidates.py
```

factors the stored relations, extracts those degree-`48` repeated-support
products, and ranks candidates by CRT rational-reconstruction success.  The
five-prime ranking is stored in

```text
code/sage/out/repeated_support_lowheight_bank_v1_s8_79656_rank_32003_32009_32027_32029_32051.json
```

Top five-prime rows:

| Trial | Linear form | Reconstructed coefficients | Integer max log10 |
| ---: | --- | ---: | ---: |
| `12` | `a-2*b+3*c-d+2*x+3*y+e+2*q` | `36/49` | `22.225` |
| `8` | `3*a+b+2*c-d+2*x-3*y+e+2*q` | `35/49` | `22.201` |
| `6` | `a+3*b-2*c+2*d+x-y+3*e+2*q` | `33/49` | `22.200` |
| `10` | `2*a+3*b+c-2*d+x+2*y-e+2*q` | `33/49` | `22.218` |
| `5` | `3*a-2*b+c+d-x+2*y+e+q` | `32/49` | `22.214` |

The best materialized CRT summary is

```text
code/sage/out/repeated_support_lowheight_bank_v1_trial12_s8_79656_crt_32003_32009_32027_32029_32051.json
```

This improves the earlier random four-prime repeated support (`30/49`) and the
small-prime low-height test (`31/49` over four primes), but it is not yet a
usable rational lift: after five primes, `13` of the `49` coefficients still
fail rational reconstruction and the centered integer lift remains near the
full modulus.  The current interpretation is that low-height linear forms help,
but the degree-`48` repeated support is still not exposing a genuinely small
rational polynomial with this bank.  A next computational step would be either a
larger/search-driven coefficient bank or extending the best candidates to more
primes now that the action-cache pipeline is established.

The next pass used a generated height-`6` bank of `96` primitive coefficient
vectors, preserving the original `16` v1 candidates as a prefix:

```text
code/sage/generate_linear_form_bank.py
code/sage/linear_form_candidates_s8_79656_lowheight_v2_96_h6.json
code/sage/out/linear_form_lowheight_bank_v2_96_h6_s8_79656_p32003_coeffs.json
code/sage/out/linear_form_lowheight_bank_v2_96_h6_s8_79656_p32009_coeffs.json
code/sage/out/linear_form_lowheight_bank_v2_96_h6_s8_79656_p32027_coeffs.json
code/sage/out/linear_form_lowheight_bank_v2_96_h6_s8_79656_p32029_coeffs.json
code/sage/out/linear_form_lowheight_bank_v2_96_h6_s8_79656_p32051_coeffs.json
code/sage/out/repeated_support_lowheight_bank_v2_96_h6_s8_79656_rank_32003_32009_32027_32029_32051_parallel.json
```

There were `77` five-prime eligible rows out of `96`.  Trial `70`,
`2*a+3*b+3*c+3*d+4*x-4*y+2*e-q`, is a repeatable structural exclusion with
degree `2776` on all five primes.  The remaining exclusions are prime-specific
degree drops or nonfinds.  For every eligible row, the repeated-support product
again has degree `48`.

Top five-prime v2 rows:

| Trial | Linear form | Reconstructed coefficients | Integer max log10 |
| ---: | --- | ---: | ---: |
| `26` | `3*a-6*b+3*c+3*d+x-6*y+3*e+q` | `37/49` | `22.190` |
| `93` | `5*a+4*b-6*c-5*d-x+3*y+4*e-4*q` | `37/49` | `22.204` |
| `23` | `5*a+3*b-2*c-6*d+4*x-4*y+4*e+q` | `37/49` | `22.212` |
| `72` | `a-4*b-c+d-5*x+3*y-5*e-2*q` | `37/49` | `22.212` |
| `76` | `6*a-b+2*c-3*d-6*x+y-4*e-2*q` | `37/49` | `22.225` |

This is a small improvement over the v1 best `36/49`, but not yet a stable
rational lift.  To test whether the five-prime modulus was simply too small,
the five best candidates were screened at a regenerated sixth prime, `p=32057`:

```text
/tmp/gtz_plucker/plucker_ansatz_s8_79656_p32057.json
/tmp/gtz_plucker/plucker_ansatz_s8_79656_p32057_actions
code/sage/linear_form_candidates_s8_79656_lowheight_v2_top5_37of49.json
code/sage/out/linear_form_lowheight_bank_v2_top5_37of49_s8_79656_p32057_coeffs.json
code/sage/out/repeated_support_lowheight_bank_v2_top5_37of49_s8_79656_rank_32003_32009_32027_32029_32051_32057.json
```

All five top candidates remain exact degree `2800` at `p=32057`, but the
six-prime CRT reconstruction gets worse rather than better:

| Source trial | Linear form | Five primes | Six primes |
| ---: | --- | ---: | ---: |
| `72` | `a-4*b-c+d-5*x+3*y-5*e-2*q` | `37/49` | `34/49` |
| `26` | `3*a-6*b+3*c+3*d+x-6*y+3*e+q` | `37/49` | `32/49` |
| `23` | `5*a+3*b-2*c-6*d+4*x-4*y+4*e+q` | `37/49` | `32/49` |
| `76` | `6*a-b+2*c-3*d-6*x+y-4*e-2*q` | `37/49` | `29/49` |
| `93` | `5*a+4*b-6*c-5*d-x+3*y+4*e-4*q` | `37/49` | `28/49` |

Moreover, the failed exponent sets change substantially after adding `32057`.
This makes the current height-`6` repeated-support lift look like a useful
diagnostic but not a viable certificate path in its present form.  The practical
next step is either to add a different structural filter before CRT ranking, or
to move from one-variable repeated support to a multivariate/RUR reconstruction
that uses more of the quotient action data.

### The `q=5` repeated fiber

The multivariate/RUR follow-up found a stable structural component behind part
of the repeated-support signal.  In the available `s8_79656` q-eliminant
factorizations, the linear factor `q-5` occurs with multiplicity `3` at all `37`
recorded primes from `32003` through `32341`.

Direct component probes were run by adjoining `q-5` to the reduced Plucker
ansatz quotient at six primes:

```text
code/sage/probe_linear_factor_component.py
code/sage/out/component_q_factor6_qeq5_s8_79656_p32003.json
code/sage/out/component_qeq5_s8_79656_p32009.json
code/sage/out/component_qeq5_s8_79656_p32027.json
code/sage/out/component_qeq5_s8_79656_p32029.json
code/sage/out/component_qeq5_s8_79656_p32051.json
code/sage/out/component_qeq5_s8_79656_p32057.json
```

Each prime gives the same finite-field structure: augmented Groebner basis size
`29`, dimension `0`, quotient degree `64`, and coordinate eliminant degrees

```text
a,b,c,d: 6;  x,y: 8;  e: 4;  q: 1.
```

The six-prime CRT lift is stored in

```text
code/sage/out/component_qeq5_s8_79656_crt_32003_32009_32027_32029_32051_32057.json
```

All coordinate eliminants rationally reconstruct over `QQ`:

| Coordinate | Rational eliminant |
| --- | --- |
| `a` | `u^6 + 1/4*u^2` |
| `b` | `u^6 + 1/100*u^2` |
| `c` | `u^6 - 17/2*u^4 + 285/16*u^2 - 25/16` |
| `d` | `u^6 - 5/2*u^4 - 51/16*u^2 + 5/16` |
| `e` | `u^4 + 4/5*u^2 - 1/5` |
| `q` | `u - 5` |
| `x` | `u^8 - 17/10*u^6 + 69/80*u^4 - 17/40*u^2 + 49/320` |
| `y` | `u^8 - 17/10*u^6 + 249/400*u^4 - 17/1000*u^2 + 49/8000` |

Thus the unstable degree-`48` one-variable support search exposed, at least in
part, a concrete nonreduced q-fiber.  For source trial `26`, the eight repeated
quadratic factors over `p=32003` all cut degree-`2` components inside this
fiber:

```text
code/sage/out/linear_form_lowheight_bank_v2_trial26_s8_79656_p32003_factors.json
code/sage/out/components_trial26_repeated_deg2_s8_79656_p32003.json
```

Every such component has `q=5`, `e^2 = 1/5`, `c^2 = 5`, and `d^2 = -1` over the
appropriate finite-field extensions; the `a,x` and `b,y` eliminants split into
sign/factor choices.  This is a substantially better target for exact
decomposition than the full degree-`48` support product.

The exact `QQ` Groebner check is now completed via the Singular-backed helper

```text
code/sage/probe_plucker_q_value_singular.py
```

using Singular's modular rational backend `modGB("slimgb", ...)`.  The full
`q=5` fiber over `QQ` has dimension `0`, degree `64`, and Groebner basis size
`29`.  The run also reduces all eight reconstructed coordinate eliminants above
to zero in the exact quotient:

```text
code/sage/out/plucker_qeq5_singular_s8_79656_QQ_relations.json
```

This turns the `q=5` observation from modular evidence into an exact rational
certificate for this ansatz fiber.

There is also a short real-exclusion certificate.  Over `R`, the exact
eliminants

```text
a^2*(a^4 + 1/4) = 0,   b^2*(b^4 + 1/100) = 0
```

force `a=b=0`.  Adding `a=b=0` to the exact `q=5` fiber gives a
zero-dimensional degree-`16` quotient, and the exact normal-form check reduces
`e^2+1` to zero:

```text
code/sage/out/plucker_qeq5_ab0_singular_s8_79656_QQ.json
```

Hence the `q=5` fiber has no real points in the reduced Plucker ansatz.  It is
a rigorously removable nuisance branch, not a source of real extremizers.

The same exact fixed-fiber probe also eliminates the other stable rational
linear q-factor.  Across all `37` available q-factorizations, the only rational
linear q-values of height at most `30` are `q=1` and `q=5`; `q=1` occurs with
multiplicity `1` at every prime.  The exact `QQ` run

```text
code/sage/out/plucker_qeq1_singular_s8_79656_QQ.json
```

has dimension `0`, degree `16`, and Groebner basis size `23`.  Its stored basis
contains `y^2+1`, so the `q=1` fiber also has no real points in the reduced
Plucker ansatz.  After these two checks, there are no remaining stable rational
q-fibers in the recorded q-eliminant data; the next exact targets would have to
come from non-rational algebraic q-factors or from a different structural
coordinate.

A first non-rational q-factor screen was then run after removing the exact
`q=1` and `q=5` branches from every modular factorization:

```text
code/sage/screen_q_factor_targets.py
code/sage/out/qfactor_target_screen_s8_79656_exclude_q1_q5_D120.json
```

For target degrees `2,3,4,5`, the maximum number of modular factor-products per
prime is still manageable:

| Degree | Min products at a prime | Max products at a prime | Sum over 37 primes |
| ---: | ---: | ---: | ---: |
| `2` | `6` | `155` | `1616` |
| `3` | `10` | `853` | `6095` |
| `4` | `16` | `3386` | `18383` |
| `5` | `22` | `10391` | `46749` |

The factor-product CRT recovery was rerun with `q=1,5` excluded:

```text
code/sage/out/recover_s8_79656_qfactor_deg2_exclude_q1_q5_37primes_beam1000.json
code/sage/out/recover_s8_79656_qfactor_deg3_exclude_q1_q5_37primes_beam1000.json
code/sage/out/recover_s8_79656_qfactor_deg4_exclude_q1_q5_37primes_beam1000_sorted.json
code/sage/out/recover_s8_79656_qfactor_deg5_exclude_q1_q5_37primes_beam500_sorted.json
```

All four are negative in the same sense as the earlier degree-`43`/`44` tests:
the best centered coefficient sizes at 37 primes are `10^165.587`, `10^165.827`,
`10^165.437`, and `10^165.805`, respectively, i.e. they track the CRT modulus
instead of stabilizing.  Thus, after removing the two rational nuisance fibers,
there is no evidence for an exact non-rational q-factor of degree `2` through
`5`.  Degree `6` and above rapidly becomes a high-multiplicity product search,
so the next step should not be blind degree-by-degree q-factor CRT unless a new
component invariant narrows the candidates.

The same Singular-backed helper was also extended to run without fixing `q`.
As a smoke check, the full reduced Plucker ansatz over `F_32003` reproduces the
known dimension-`0`, degree-`2832`, basis-size-`1901` quotient in about one
minute:

```text
code/sage/out/plucker_full_singular_s8_79656_p32003_smoke.json
```

The analogous exact `QQ` monolithic run was tried with the same modular backend
and a `3600` second timeout:

```text
code/sage/out/plucker_full_singular_s8_79656_QQ_probe.json
```

It timed out without producing a Groebner basis.  Thus fixed rational fibers are
now easy to certify exactly, but the full degree-`2832` quotient still needs
pre-decomposition, modular reconstruction, or a different structural coordinate
before attempting characteristic-zero certification.

## Relaxed interval bound calibration

There is already a clean analytic relaxed baseline from the standard
volume-sampling expectation formula.  If `A in St(n,k)` and one samples `k`
rows with probability proportional to `det(A_I^T A_I)`, then applying the
Derezinski--Warmuth inverse moment identity to `X=A^T` gives

```text
E[(A_I^T A_I)^(-1)] <= (n-k+1) I_k.
```

Hence some sampled set satisfies
`tr((A_I^T A_I)^(-1)) <= k(n-k+1)`.  Since
`A_I^T A_I <= I_k`, its nonzero eigenvalues are at most one, so the other
`k-1` inverse eigenvalues contribute at least `k-1` to this trace.  Therefore

```text
sigma_min(A_I)^2 >= 1/(k(n-k)+1).
```

For `(n,k)=(6,3)` this is exactly `F(P) >= 1/10`.  This is not new to the
present project; it is a SOTA baseline consequence of volume sampling.  The
computational relaxed-threshold experiments below should therefore be read as
diagnostics for whether the interval machinery can approach the sharp `1/6`
frontier, not as the natural proof of the `1/10` relaxation.

There is now a stronger analytic relaxed result for the first open case.  Let
`t_* = (1 - sqrt(3/5))/2 = 0.112701665...`.  Then

```text
GTZ(6,3) relaxed with alpha < 3/2:
    F(P) >= t_* > 1/9  for every rank-three projector P in R^6.
```

The proof is short.  For any `t < t_*`, if some leverage satisfies
`ell_i <= 1 - 5t`, delete that row
and use the known exact `(5,3)` case, dual to `k=2`; scaling back gives
`sigma_min^2 >= (1-ell_i)/5 >= t`.  Otherwise all leverages are `>1 - 5t`, so
every triple has trace `>3(1 - 5t)`.  Since `t_* < 2/17`, a bad triple cannot
have two eigenvalues at or below `t`: `1 + 2t < 3(1 - 5t)`.  Hence on this
high-leverage core, for every `t < t_*`,

```text
lambda_min(P_TT) <= t    iff    det(P_TT - t I) <= 0.
```

But

```text
sum_T det(P_TT - t I) = 1 - 12 t + 30 t^2 - 20 t^3,
```

using `sum e3=1`, `sum e2=12`, and `sum e1=30` over the twenty triples.  At
every `t < t_*` the right-hand side is positive, so the determinants cannot all
be nonpositive.  Therefore some triple has `lambda_min > t` in the core case;
letting `t` increase to `t_*` gives `F(P) >= t_*`.
The exact constant check is recorded in `verify/v35_relaxed_one_ninth.py`.

The same determinant-sum method cannot reach `1/8`: at `t=1/8`,
`1 - 12t + 30t^2 - 20t^3 = -9/128`, and the trace-separation step also fails
because `1/8 > 2/17`.

For the next target `t=1/8`, deletion gives a useful reduced obstruction.  A
counterexample must satisfy all of the following:

```text
3/8 < ell_i < 5/8                         for every row i,
lambda_max(P_{ij,ij}) > 1/2                for every row pair {i,j},
lambda_min(P_{ij,ij}) < 1/2                for every row pair {i,j},
lambda_min(P_TT) <= 1/8 and lambda_max(P_TT) >= 7/8
                                              for every triple T.
```

The first line is one-row deletion plus duality.  The pair lines come from
deleting two rows and using the exact `(4,3)` case, again plus duality.  The
triple line is the target and its complementary-dual form.  Pair-spectrum and
complement-pair moment relaxations with these bounds remain feasible, so `1/8`
needs compatibility constraints between the different principal blocks
(projector/Plucker constraints), not just the global `e1,e2,e3` moments.

A new set of numerical probes now makes this target much more concrete.
`verify/v36_relaxed_1over8_obstruction.py` minimizes the maximum violation of
the closed reduced obstruction above over actual rank-three projectors.  On a
run with 60 random starts and 20 conference-matrix starts,

```text
verify/out/v36_relaxed_1over8_obstruction_80.json
```

the best point has maximum violation `0.0416946719` and lands at the known
out-of-family extremal with leverage pattern close to `(5/14)^3,(9/14)^3`.
Evaluating the certified seventh extremal directly gives the exact value
`1/24`; the TTSP extremals violate the leverage core more strongly, by `7/72`.
Thus the nearest obvious obstruction is still sharp-threshold geometry, not a
plausible `1/8` counterexample.

To remove that distraction, `verify/v37_relaxed_1over8_core_probe.py` enforces
the leverage core and pair-straddling constraints by a large quadratic penalty
and minimizes only the triple obstruction.  On

```text
verify/out/v37_relaxed_1over8_core_probe_100.json
```

68 of 100 local runs reached core violation below `1e-5`, and the best such
point still had triple violation about `5.734e-2` (`F = 0.182342191...`).
Polishing the best candidates with the explicit constrained SLSQP minimax
formulation in `verify/v38_relaxed_1over8_slsqp_polish.py` gives a cleaner
boundary:

```text
s = (5 - sqrt(17))/16 = 0.054805898398...
F = 1/8 + s = (7 - sqrt(17))/16 = 0.179805898398...
```

with Stiefel equalities and all closed core/pair inequalities satisfied to
machine precision.  Multiple independent starts polish to the same lower
boundary; several others polish to a higher local boundary
`F = 1/2 - sqrt(10)/10 = 0.183772233983...`.

The lower SLSQP boundary has a small symmetry pattern after row permutation:

```text
P =
[ 5/8  -x  -3/8  -y  -z  -z ]
[ -x    a   -x    w   y   y ]
[ -3/8 -x   5/8  -y  -z  -z ]
[ -y    w   -y   1-a -x  -x ]
[ -z    y   -z   -x  3/8 3/8]
[ -z    y   -z   -x  3/8 3/8]
```

Substituting this ansatz into `P^2=P` leaves only the equations

```text
32*x^2 + 32*y^2 + 64*z^2 - 3 = 0,
4*a*x + 4*w*y - 3*x + 8*y*z = 0,
4*a*y - 4*w*x + 8*x*z - y = 0,
a^2 - a + w^2 + 2*x^2 + 2*y^2 = 0.
```

The exact ansatz algebra is recorded in
`verify/v39_relaxed_1over8_ansatz_exact.py`.  It verifies symbolically that
adding the active determinant equations for triples `(0,1,3)` and `(0,1,4)` to
the ansatz ideal gives a Groebner basis containing

```text
8*q^2 - 7*q + 1.
```

The lower root is exactly `q = (7 - sqrt(17))/16`, with slack
`q - 1/8 = (5 - sqrt(17))/16`.

This suggests a realistic exact relaxed route:

```text
core + pair-straddling  ==>  F(P) >= (7 - sqrt(17))/16 > 1/8.
```

Together with the deletion reduction, such a certificate would prove the relaxed
`1/8` theorem.  This is still numerical evidence, not a proof, but it is a much
smaller semialgebraic target than the original exact `1/6` hypothesis.

Lowering the target threshold did not close the interval branch-and-bound
route.  With cascade projectors, hybrid bounds, all 20 charts, `max_boxes=2000`,
and `min_radius=0.25`, the center-high runs gave:

| Output | Threshold | Certified | GTZ-certified | Split | Queue | Inverse failures |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `verify/out/v19_relaxed_010_cascade_2000.json` | `0.1` | `840` | `628` | `1160` | `340` | `167` |
| `verify/out/v19_relaxed_0125_cascade_2000.json` | `0.125` | `836` | `700` | `1164` | `348` | `167` |

A matched small calibration was run against the same cascade/hybrid setup, all
20 charts, `max_boxes=500`, `min_radius=0.25`, and center-high priority:

| Output | Threshold | Certified | GTZ-certified | Split | Queue | Inverse failures |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `verify/out/v19_exact_cascade_500_min025.json` | `1/6` | `152` | `152` | `348` | `216` | `167` |
| `verify/out/v19_relaxed_1over12_cascade_500.json` | `1/12` | `152` | `110` | `348` | `216` | `167` |
| `verify/out/v19_relaxed_010_cascade_500.json` | `0.1` | `151` | `116` | `349` | `218` | `167` |
| `verify/out/v19_relaxed_1over8_cascade_500.json` | `1/8` | `152` | `127` | `348` | `216` | `167` |

At this granularity even the analytic `1/10` target is essentially no easier
for the interval cover than the sharp target: the same inverse/enclosure
failures dominate.

Adding the determinant-domain filter did not help at this granularity:

| Output | Priority | Threshold | Certified | Outside | Queue | Inverse failures |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| `verify/out/v19_relaxed_010_domain_high_2000.json` | center-high | `0.1` | `0` | `879` | `262` | `1071` |
| `verify/out/v19_relaxed_010_domain_low_2000.json` | center-low | `0.1` | `0` | `472` | `408` | `1528` |

Thus the present interval route is limited by chart/enclosure conditioning near
the chart-domain frontier, not by the exact value of the relaxed threshold.
