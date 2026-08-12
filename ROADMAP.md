# GTZ(6,3) — the whole story as a ladder

Written for an engineer's eye: what the problem is, what rungs are already
nailed down, what rung we are standing on, and exactly what lies between here and
the hypothesis. Statuses are the strict tags of the brief: **PROVED** /
**NUMERICAL** / **OPEN**.

---

## 0. The problem in one box

```
   Given:  A is 6x3 with orthonormal columns  (A^T A = I_3)
   Claim:  you can always pick 3 of the 6 rows so that the 3x3 submatrix
           has smallest singular value >= 1/sqrt(6)
```

Engineering reading: out of 6 sensors, some 3 of them always observe the full
3-dimensional state with a guaranteed conditioning floor. No arrangement of
sensors can make every triple ill-conditioned.

Standard reformulation used throughout (removes the `O(3)` gauge freedom):

```
   P = A A^T                      rank-3 projector, 6x6, lives on Gr(3,6), dim 9
   F(P) = max over the 20 triples T of  lambda_min(P_TT)
   GTZ(6,3)  <=>  min over Gr(3,6) of F  >=  1/6
```

---

## 1. The ladder

```
 ┌────────────────────────────────────────────────────────────────────────────┐
 │ RUNG 0   Problem posed. GTZ(n,k) hypothesis.                               │
 │          k=2 case proved (Sengupta–Pautov); n-2 by duality.                │
 │          => (6,3) is the smallest open case.              [LITERATURE]     │
 └────────────────────────────────────────────────────────────────────────────┘
                                    │
 ┌────────────────────────────────────────────────────────────────────────────┐
 │ RUNG 1   Reductions and dead ends mapped by the earlier campaign.          │
 │          Case A: a light row can be removed -> core region only.           │
 │          Single-pair / second-moment program: PROVABLY DEAD (Thm BND).     │
 │          Slice framework built; PMC found.                 [PROVED]        │
 └────────────────────────────────────────────────────────────────────────────┘
                                    │
 ┌────────────────────────────────────────────────────────────────────────────┐
 │ RUNG 2   Independent exact re-verification of that inheritance.            │
 │          Slice lemmas re-proved as GENERIC identities (stronger than the   │
 │          source's sample checks). Two improvements found:                  │
 │            * G_PMC = (1 - sin36)/2 exactly  (was: root of a messy cubic)   │
 │            * G_PMC > 1/6 by the rational certificate 169 < 405             │
 │                                                            [PROVED]       │
 └────────────────────────────────────────────────────────────────────────────┘
                                    │
 ┌────────────────────────────────────────────────────────────────────────────┐
 │ RUNG 3   THE GOVERNING FACT: the bound is TIGHT.                           │
 │          F = 1/6 EXACTLY at the Nesterenko extremal.                       │
 │          => no slack anywhere; any certificate must vanish there.          │
 │          => the equal-leverage slice is NOT the hard region (sits 0.039     │
 │             above 1/6) -- Track B is a standalone lemma, not a stepping     │
 │             stone.                                          [PROVED]       │
 └────────────────────────────────────────────────────────────────────────────┘
                                    │
 ┌────────────────────────────────────────────────────────────────────────────┐
 │ RUNG 4   Local behaviour at an extremal fully understood.                  │
 │          F has a KINK, not a smooth bowl: it grows LINEARLY away from the   │
 │          extremal (rate kappa > 0), not quadratically. Certified by an      │
 │          exact rational KKT certificate: multipliers (7/90)x9, 3/10, all    │
 │          STRICTLY positive.                                 [PROVED]       │
 │          => upgrades part (a) of Reformulation R from NUMERICAL to PROVED.  │
 └────────────────────────────────────────────────────────────────────────────┘
                                    │
 ┌────────────────────────────────────────────────────────────────────────────┐
 │ RUNG 5   Same certificate at ALL SEVEN known extremals.                    │
 │          Six from the graph ("scaled-star") family + one found by search   │
 │          that is PROVABLY OUTSIDE that family (leverages 5/14, 9/14, 13     │
 │          active triples, entries in Q(sqrt2,sqrt5)).        [PROVED]       │
 │          Structural news: the corpus's catalogue of nine is INCOMPLETE, and │
 │          the out-of-family one is the MOST COMMON basin in search.         │
 └────────────────────────────────────────────────────────────────────────────┘
                                    │
 ┌────────────────────────────────────────────────────────────────────────────┐
 │ RUNG 6   Finiteness of the extremal set, CONDITIONALLY.                    │
 │          sharp => isolated => (with compactness) FINITE.    [PROVED]       │
 │          Explicit radius r0 = kappa * gap_min around each of the seven;     │
 │          gives a packing bound |E| <= (2/r0)^9, finite and explicit.        │
 │          Conditional ONLY on a uniform lower bound for r0 -- NOT on         │
 │          enumerating E.                                     [PROVED at 7]  │
 └────────────────────────────────────────────────────────────────────────────┘
                                    │
                          ==== WE ARE HERE ====
                                    │
 ┌ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ┐
 │ RUNG 7   Uniform lower bound on r0 over ALL of E.              [OPEN]     │
 │          Needs kappa and gap_min bounded below globally. Both are          │
 │          semialgebraic in P, so this is a Positivstellensatz question --    │
 │          but on a MUCH smaller object than Rung 9, since it only concerns   │
 │          local data at the level set.                                      │
 │          PAYOFF: an unconditional finite bound on |E|.                     │
 └ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ┘
                                    │
 ┌ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ┐
 │ RUNG 8   Algebraic dimension of E = 0.                   [PARTLY DONE]    │
 │          E is a real algebraic set. If every stratum (one per active-set   │
 │          pattern) has dimension 0, then E is FINITE outright, with NO      │
 │          sampling -- immune to an unsampled extremal.                     │
 │                                                                           │
 │          DONE: local dim = 0 at all seven known extremals -- Jacobian      │
 │          rank 9/9 in a 9-variable chart, smallest singular value >= 0.5.   │
 │          This is an INDEPENDENT route to isolation (rank/equality) vs      │
 │          Rung 6's (cone/inequality), so the two cross-check.  [PROVED at 7]│
 │                                                                           │
 │          KEY REDUCTION (simple-active): sharpness needs positive spanning  │
 │          in R^9, so SHARP => |A| >= 10. A simple-active extremal with      │
 │          <= 9 active triples is automatically non-sharp.                   │
 │          Observed so far: |A| in {10,12,13}. The nonsimple active case     │
 │          must be split off separately.                                     │
 │                                                                           │
 │          NEW: the low-active route is now a FINITE LIST. Of the 431,909    │
 │          subsets with |A| <= 9, exactly 124 S_6-orbits survive two exact   │
 │          necessary conditions (cover [6]; cover all 15 pairs). None with   │
 │          |A| <= 5 survives; the minimum is a UNIQUE |A|=6 orbit with a     │
 │          perfect-matching structure.                                       │
 │                                                                           │
 │          MECHANISM IDENTIFIED: low-active orbits split into two numerical  │
 │          certificate targets. For reachable active-equality loci, active   │
 │          PSD is FINE (slack > 0) but some triple OUTSIDE A overshoots      │
 │          1/6 by roughly 5e-2 to 1.2e-1, so the active set must grow.        │
 │          Other orbits do not reach the active-equality locus numerically,  │
 │          making equality-locus emptiness the likely certificate target.    │
 │          ALL 124 ORBITS NOW CLASSIFIED (size 9 finished on MLCore):        │
 │            70 reach their locus -> ALL 70 are inactive-growth              │
 │             0 active-PSD obstructions,  0 feasible points                  │
 │            54 never reach their locus numerically (real-emptiness likely,  │
 │               not established; complex loci ARE nonempty by Singular)      │
 │          Margin verified by direct minimisation on reachable loci: unique  │
 │          size-6 orbit, all five reachable size-7 orbits, and first four    │
 │          smallest-excess size-8 orbits all return MARGIN-HOLDS. Current    │
 │          tested minimum is 2.20e-2 at size-7 canon 14113 (13% of the       │
 │          threshold).                                            [OPEN]     │
 │          NEW (v34-v36): the tightest margin minimizers satisfy strict      │
 │          nonsmooth KKT balance numerically (residual ~5e-16; outside       │
 │          weights positive). Direct Groebner KKT still times out, and       │
 │          explicit rank minors are worse, but high-precision Newton on the  │
 │          square active+tie systems gives residuals down to 1e-148. NEW:    │
 │          Krawczyk interval certificates now prove unique local roots +     │
 │          branch inequalities for all six non-over-tied tested size-7/8     │
 │          margin minimizers. For over-tied cases, full Gauss-Newton now     │
 │          separates three full-rank size-7/8 roots from one singular size-6 │
 │          root. All-tie square bases are Krawczyk-certified for the         │
 │          full-rank roots, but omitted determinant equalities remain only   │
 │          interval-enclosed at ~1e-54. New modular local-separator search   │
 │          finds quartic separators at D=5 for both clean cases. For the     │
 │          size-7 clean case, one QQ lift is an exact member but vanishes at │
 │          the root, while the nonzero raw lift is not a QQ member; the      │
 │          missing step is now rational reconstruction/direct QQ separator.  │
 │          The disjunction is IRREDUCIBLE (v29 set cover): needs 10 of 14    │
 │          outside triples at tau=1e-9, all 14 at tau=3e-2; three DIFFERENT  │
 │          forced witnesses. So ~1700 regions stands; it cannot be shrunk.   │
 │          NEW (v28): the growth certificate is a DISJUNCTION, not one       │
 │          inequality. The witness triple VARIES over the locus (12+         │
 │          witnesses on the |A|=6 orbit; at worst points only ONE outside    │
 │          triple overshoots), and NO single outside triple has all three    │
 │          char-poly coefficients positive everywhere. So the target is      │
 │          ~124 x 14 ~ 1700 regions needing a CAD-style covering, not a      │
 │          single SOS identity. The old 4.5e-2 margin was too optimistic;    │
 │          the observed positive-margin mechanism itself still survives.     │
 │          Three routes died first: lex elimination (degree 1880 vs 18 for   │
 │          known-base), leverage/trace conditions (0 of 124 killed), and     │
 │          the single-inequality SOS target (witness not unique).            │
 └ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ┘
                                    │
 ┌ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ┐
 │ RUNG 9   THE GLOBAL GAP -- part (b) of Reformulation R.         [OPEN]    │
 │          "No critical point of F anywhere on Gr(3,6) has value < 1/6."     │
 │          THIS IS THE HYPOTHESIS. Everything else is scaffolding.          │
 └ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ┘
                                    │
                              GTZ(6,3) PROVED
```

---

## 2. What actually lies between us and the hypothesis

Only **one** thing, stated three ways:

| view | statement |
|---|---|
| optimization | no critical point of `F` on `Gr(3,6)` has value `< 1/6` |
| geometry | the sublevel set `{F < 1/6}` is empty |
| algebra | a Positivstellensatz certificate for `F >= 1/6` on a 9-dim variety |

Rungs 7 and 8 are **not** on the critical path to Rung 9 — they are the two
routes that *sharpen* it, and both are strictly easier than it:

```
   Rung 7 (uniform r0)  ──┐
                          ├──>  |E| is finite and explicitly bounded
   Rung 8 (dim E = 0)   ──┘     => the extremals become a FINITE LIST
                                => Rung 9 becomes "check a finite list +
                                   exclude explicit neighbourhoods", which is
                                   the FIRST version of Rung 9 that is even
                                   plausibly a finite computation
```

That is the whole strategic point: **as stated, Rung 9 is an infinite
verification. Rungs 7–8 are what could turn it into a finite one.** Pavel's
original hope ("a finite number of algebraic checks") is exactly right — but it
only becomes available *after* Rung 7 or 8.

---

## 3. Why it is hard, in engineering terms

1. **No stability margin.** The inequality is *tight* (Rung 3): `F = 1/6` exactly
   at seven known points. Any numerical scheme with an `epsilon` tolerance cannot
   distinguish "holds with equality" from "just fails". Everything must be exact
   at the boundary.
2. **The cost function has kinks.** `F` is a max of 20 eigenvalue functions, so it
   is not differentiable at the interesting points. Smooth machinery (gradients,
   Hessians, ordinary SOS) does not directly apply; we needed the nonsmooth
   first-order theory (critical cones, KKT with strict multipliers).
3. **The known extremal catalogue is wrong.** Not slightly — the *most common*
   extremal found in search is absent from it (Rung 5). Any proof that works by
   "check the nine known cases" is unsound as stated.
4. **Sampling cannot close it.** Rung 6's conditional and Rung 9 have the same
   logical shape: `for all P` over a 9-dimensional manifold. No number of random
   starts settles a `forall`. This is why Rung 8 (algebraic, sampling-free) is the
   most promising direction.

---

## 4. Current state of the two live computations

| target | status | what would count as success |
|---|---|---|
| **Rung 7** — uniform `r0` | radii computed at all seven: `r0 >= 0.00996`, `kappa >= 0.046`, `gap_min >= 4/21` (`v13`) | a *proved* global lower bound (currently only at seven points) |
| **Rung 8** — `dim E = 0` | **local dimension = 0 at all seven**, Jacobian rank 9/9 with smallest singular value >= 0.5 (`v16`) | positive spanning at every simple-active stratum, no simple-active extremal with `|A| <= 9`, and a separate nonsimple-active exclusion |

**Rung 8 progress worth highlighting.** Two results moved it:

1. **A proved simple-active reduction.** At a simple-active equality point,
   sharpness means the active gradients positively span the 9-dimensional tangent
   space. Therefore **sharp => |A| >= 10**, i.e. any simple-active extremal with
   `<= 9` active triples is automatically non-sharp. Observed so far: always
   10, 12, or 13. The nonsimple active case is separate.
2. **An independent confirmation of isolation.** Rung 6 got isolation from the
   critical cone (an inequality/LP argument). `v16` gets it again from a **rank**
   condition (Jacobian rank 9 in a 9-variable chart) — a different kind of
   argument, so the agreement is a real cross-check, not a restatement. 7/7 at
   full rank, smallest singular value `>= 0.5`, nowhere near borderline.

3. **Low-active obstruction split now has witnesses for sizes 6 and 7.** The
   new artifact `verify/out/v27_obstruction_split_s6s7.json` records 60-start
   runs for the 8 size-6/7 orbits. Classification: `inactive=6`, `off=2`,
   `psd=0`, `both=0`, `neither=0`. Inactive-growth certificate targets:
   `78593 -> 4=(0,2,3)`, `13105 -> 12=(1,2,5)`,
   `14113 -> 7=(0,3,4)`, `78595 -> 7=(0,3,4)`,
   `78601 -> 14=(1,3,5)`, `78612 -> 11=(1,2,4)`.
   The two off-locus size-7 canons are `4040` and `8065`, with best equality
   residuals `9.62e-5` and `1.81e-3`. The existing size-8 artifact
   `verify/out/v27_obstruction_split.json` has classification `inactive=22`,
   `off=10`, `psd=0`; it predates witness-field output.

4. **Size 9 has started.** The checkpointed artifact
   `verify/out/v27_obstruction_split_s9_000_011.json` covers the first 12
   sorted size-9 canons with `starts=20`, `maxiter=5000`, and `complete=true`.
   All 12 are off-locus in this first pass: `2029`, `4043`, `4046`, `4058`,
   `4060`, `8071`, `8083`, `8085`, `8086`, `8092`, `12121`, `12124`.
   Residual range is `1.86e-4` to `2.05e-2`. Canon `4060` was rerun at
   `starts=80`, `maxiter=8000` in
   `verify/out/v27_obstruction_split_s9_canon4060_deep.json` and remained
   off-locus with the same best residual `1.86e-4`. A direct modular no-section
   determinant-ideal probe for this active set timed out after 600 seconds
   (`code/sage/out/A9_4060_probe.json`), so the exact route should avoid
   brute-force no-section Grobner here.

Supporting searches running: `v12_hunt` on MLCore (3000 starts; at last check
400/3000, **91 hits, still exactly 4 leverage patterns**, no non-sharp candidate),
and `v14_biased` locally (981 biased starts: perturbations of known extremals,
leverage-targeted profiles, icosahedral orbits) to probe where uniform sampling
is thin.

Counterexample tripwire: **never fired**, across every run of this project.

**Compute note (answering the question about MLCore limits).** The 26-CPU wall we
hit is the *default project quota*, not a platform limit — the KB (`/jobs/limits`)
documents defaults of 24-32 CPU per region and says a bigger quota is requested via
the `~ml-core-ask` channel, including **temporary quota** for resource-heavy
projects. So the student is right that this is liftable. Two other documented
escapes, both relevant here: **Preemptible jobs are NOT quota-limited at all**
(30-day duration, max 7 days without preemption) — ideal for an embarrassingly
parallel search that tolerates restarts — and **batch jobs** are the documented
pattern for "many similar experiments differing by parameters", which is exactly
the search shape. Worth pursuing before assuming 20 cores is the ceiling.

---

## 5. One-line summary

We have nailed down *where* the problem is tight and proved that the local
behaviour there is as good as it could be (sharp, certified, isolated). What
remains is a single global statement over a 9-dimensional manifold, and the two
live computations are attempts to make that statement **finite** rather than to
verify it directly.
