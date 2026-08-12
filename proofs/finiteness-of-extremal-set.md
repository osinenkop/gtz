# Sharpness ⟹ isolation ⟹ finiteness: the structure of the GTZ(6,3) extremal set

**Deliverable:** answers Pavel's question "is the extremal set finite at all?".
Scripts: `verify/v9_rational_certificate.py`, `verify/v11_seventh_exact.py`,
`verify/v12_hunt.py`. Depends on `proofs/sharp-cone-at-extremal.md` §4a–§4c.1.

---

## 0. Scope — read this first (overclaim guard)

> **This document does NOT prove GTZ(6,3), and does NOT prove that the extremal
> set is finite.**
>
> It proves a **conditional** finiteness theorem (§1): *if every point of
> $\mathcal E=\{P\in Gr(3,6):F(P)=\tfrac16\}$ is a sharp minimum, then
> $\mathcal E$ is finite.* The hypothesis is verified **exactly at seven
> configurations** and **numerically at every configuration any search has ever
> found**, but it is *not* established for all of $\mathcal E$ — that would
> require a global argument of exactly the kind part (b) of Reformulation R is
> missing.
>
> The value of the theorem is that it makes the question **decidable per
> extremal** with machinery already built, and converts an open-ended hunt into a
> search with a definite success criterion.

| # | Statement | § | Status |
|---|---|---|---|
| F1 | Sharp $\Rightarrow$ isolated in $\mathcal E$ | §1.1 | **PROVED** |
| F2 | All points of $\mathcal E$ sharp $\Rightarrow\mathcal E$ finite | §1.2 | **PROVED** |
| F3 | A positive-dimensional family in $\mathcal E$ forces a point with nontrivial critical cone | §1.3 | **PROVED** (contrapositive of F1) |
| F4 | Seven explicit extremals are sharp, hence isolated | §2 | **PROVED** (exact; `v9`, `v11`) |
| F5 | Only **four** leverage patterns occur among all extremals ever found | §3 | **NUMERICALLY SUPPORTED** |
| F6 | $\mathcal E$ is finite | §4 | **OPEN** (conditional on F2's hypothesis) |

---

## 1. The finiteness theorem

Throughout $F(P)=\max_{|T|=3}\lambda_{\min}(P_{TT})$ on the compact manifold
$Gr(3,6)=\{P\in\operatorname{Sym}(6):P^2=P,\operatorname{tr}P=3\}$, and

$$\mathcal E:=\{P\in Gr(3,6):F(P)=\tfrac16\}.$$

Call $P_0\in\mathcal E$ **sharp** if
$\kappa(P_0):=\min_{\|\dot P\|=1}F'(P_0;\dot P)>0$, where $F'$ is the one-sided
directional derivative.

### 1.1 Sharp ⟹ isolated [PROVED]

> **Lemma F1.** If $P_0\in\mathcal E$ is sharp, there is $r_0>0$ such that
> $F(P)>\tfrac16$ for all $P\in Gr(3,6)$ with $0<\operatorname{dist}(P,P_0)<r_0$.
> Hence $P_0$ is an isolated point of $\mathcal E$.

*Proof.* At $P_0$ the active set $\mathcal A=\{T:\lambda_{\min}(P_{0,TT})=\tfrac16\}$
is finite, and (verified in every certified case) $\lambda_{\min}$ is a **simple**
eigenvalue of each active block. A simple eigenvalue of a symmetric matrix is an
analytic function of its entries, so each $P\mapsto\lambda_{\min}(P_{TT})$,
$T\in\mathcal A$, is analytic near $P_0$; the inactive blocks stay strictly below
$\tfrac16$ on a neighbourhood by continuity. Hence near $P_0$

$$F=\max_{T\in\mathcal A}g_T,\qquad g_T\ \text{analytic},\ g_T(P_0)=\tfrac16 ,$$

a max of finitely many analytic functions. Such an $F$ is semismooth with
$$F(P_0+r\dot P)=\tfrac16+r\,F'(P_0;\dot P)+R(r,\dot P),\qquad
|R(r,\dot P)|\le Cr^2$$
with $C$ uniform over $\|\dot P\|=1$ (finitely many analytic $g_T$, sphere
compact). With $F'(P_0;\dot P)\ge\kappa$ for all unit $\dot P$,

$$F(P_0+r\dot P)\ \ge\ \tfrac16+\kappa r-Cr^2\ >\ \tfrac16
\qquad\text{for }0<r<r_0:=\kappa/C .$$

So a punctured $r_0$-ball around $P_0$ misses $\mathcal E$. $\square$

Note this holds **regardless** of whether $\tfrac16$ is the global minimum: it is a
statement about the level set, not about optimality.

### 1.2 All sharp ⟹ finite [PROVED]

> **Theorem F2.** If every point of $\mathcal E$ is sharp, then $\mathcal E$ is
> finite.

*Proof.* $F$ is continuous, so $\mathcal E=F^{-1}(\tfrac16)$ is closed in
$Gr(3,6)$, which is compact; hence $\mathcal E$ is compact. By Lemma F1 every
point of $\mathcal E$ is isolated, i.e. $\mathcal E$ is discrete. A compact
discrete space is finite. $\square$

### 1.3 The contrapositive — what to hunt for [PROVED]

> **Corollary F3.** If $\mathcal E$ contains a $C^1$ curve $t\mapsto P(t)$ with
> $P(0)=P_0$ and $\dot P(0)\ne0$ — in particular if $\mathcal E$ contains a
> positive-dimensional family — then $F'(P_0;\dot P(0))=0$, so the critical cone
> $$N(P_0)=\{\dot P:F'(P_0;\dot P)\le0\}=\{\dot P:L\dot P\le0\}$$
> is **nontrivial** and $P_0$ is not sharp.

*Proof.* $F(P(t))\equiv\tfrac16$ gives
$0=\frac{d}{dt}\big|_{0^+}F(P(t))=F'(P_0;\dot P(0))$. $\square$

**This is the operational content.** Nontriviality of $N$ is exactly what the
certificate of `sharp-cone-at-extremal.md` §4b decides, in exact arithmetic:
$N=\{0\}$ iff $\operatorname{rank}L=9$ **and** the exact LP optimum of
$\sum_i(-Lz)_i$ over $\{Lz\le0,\|z\|_\infty\le1\}$ is $0$. So:

- every extremal certified sharp is an isolated point, contributing $1$ to
  $|\mathcal E|$;
- a positive-dimensional family would force at least one non-sharp extremal.
  The converse is not automatic: a non-sharp point may still be isolated for
  higher-order reasons.  Still, non-sharp equality points are the right objects
  to hunt, because every failure of compact-discrete finiteness must pass
  through them.

The hunt therefore has a definite success criterion, which is why it is worth
running rather than being open-ended.

---

## 2. Seven isolated points [PROVED]

Combining Lemma F1 with the exact certificates already in hand:

| configuration | leverages | active | source |
|---|---|---|---|
| `P(S(e,e,e),e,e,e)` | $\tfrac{13}{18}^{\times3},\tfrac5{18}^{\times3}$ | 10 | `v6`, `v9` |
| `P(S(e,e),S(e,e),e,e)` | $\tfrac{11}{18}^{\times4},\tfrac5{18}^{\times2}$ | 12 | `v9` |
| `P(S(P(e,e,e),e),S(e,e))` | $\tfrac5{18}^{\times3},\tfrac{13}{18}^{\times3}$ | 10 | `v9` |
| `P(S(P(S(e,e),e,e),e),e)` | $\tfrac{11}{18}^{\times4},\tfrac5{18}^{\times2}$ | 12 | `v9` |
| `P(S(P(e,e),P(e,e)),S(e,e))` | $\tfrac7{18}^{\times4},\tfrac{13}{18}^{\times2}$ | 12 | `v9` |
| `P(S(P(e,e),e),S(P(e,e),e))` | $\tfrac7{18}^{\times4},\tfrac{13}{18}^{\times2}$ | 12 | `v9` |
| out-of-family (§4c.1) | $\tfrac5{14}^{\times3},\tfrac9{14}^{\times3}$ | 13 | `v11` |

> **Corollary F4 [PROVED].** Each of these seven is an isolated point of
> $\mathcal E$: there is an explicit punctured neighbourhood of each on which
> $F>\tfrac16$ strictly.

**Counting caveat.** $|\mathcal E|$ is *not* seven even if $\mathcal E$ is finite.
Relabelling rows ($S_6$ acting by simultaneous row/column permutation of $P$)
carries an extremal to a generally *different* point of $Gr(3,6)$ that is also
extremal. So each configuration above represents an $S_6$-orbit of up to $720$
points (fewer when the configuration has symmetry). The census also shows several
distinct configurations sharing a leverage pattern, so "leverage pattern" is a
strictly coarser invariant than "extremal".

---

## 3. Search evidence [NUMERICALLY SUPPORTED]

`verify/v12_hunt.py` — massive multi-start descent on $F$, each hit polished to
$|F-\tfrac16|<10^{-9}$ (target $10^{-14}$), clustered by sorted leverage pattern,
then screened for sharpness (active count, simplicity of $\lambda_{\min}$,
$\operatorname{rank}L$, LP optimum, and a $\kappa$ estimate from 120 minimizations).

Pilot run (300 starts, master seed 20260731): **72 hits**, 227 misses, 1 near.
Exactly **four** distinct leverage patterns, and **every one screens SHARP**:

| pattern | multiplicity | active | $\operatorname{rank}L$ | LP opt | $\kappa$ (est) |
|---|---|---|---|---|---|
| **$\tfrac5{14},\tfrac9{14}$ (out of family)** | **25** | 13 | 9 | $0$ | $0.0767$ |
| $\tfrac7{18},\tfrac{13}{18}$ | 20 | 12 | 9 | $0$ | $0.0949$ |
| $\tfrac5{18},\tfrac{13}{18}$ | 16 | 10 | 9 | $0$ | $0.0773$ |
| $\tfrac{11}{18},\tfrac5{18}$ | 11 | 12 | 9 | $0$ | $0.0774$ |

Three observations worth recording:

1. **No new patterns.** The four found are exactly the four already certified
   exactly (three TTSP patterns covering the six census members, plus the
   out-of-family seventh). Unconstrained descent discovers nothing else.
2. **No non-sharp candidates.** Every representative has $\operatorname{rank}L=9$,
   LP optimum $0$, and $\kappa\approx0.08$–$0.09$ — comfortably bounded away from
   $0$, not marginal. Via F1 this makes each an isolated point.
3. **The out-of-family extremal is the *most common* basin** (25 of 72, more than
   any scaled-star pattern). It is not an exotic corner case; it is the *typical*
   thing unconstrained descent lands on. That the corpus's catalogue omits it
   entirely is therefore a substantive gap, not a technicality.

**Scale-up in progress.** A 3,000-start confirmation is running on MLCore (job
`gtz63-hunt-ivtoc5`, `20cpu-320ram`, region `ix-m5-sm11`). At the last check
(400/3000): **91 hits, still exactly 4 patterns**, no new pattern and no non-sharp
candidate. A complementary local run (`v14_biased.py`, 981 *biased* starts —
perturbations of known extremals, leverage-targeted profiles across a
heavy/light grid, and icosahedral orbits) is also under way, specifically to probe
regions uniform-random descent under-samples.

*Two operational notes, both recorded because they cost time.* (i) The project's
MLCore quota is **26 CPU**, so the `110cpu-900ram` flavor is rejected outright.
(ii) A first 20,000-start submission would have needed ~13.6 h against a 5 h wall
clock and, because the script only reported at the end, would have been killed with
**no output at all**. `v12` now checkpoints every 100 starts to
`verify/out/v12_progress.json`, and the job was resized to fit. Cost estimate for
future runs: **~47 core-seconds per start**.

Tripwire: no sample with $F<\tfrac16-10^{-6}$ in any run.

---

## 3a. Certified radii, and the packing bound [PROVED at the seven]

`verify/v13_radius.py`. Lemma F1 gives isolation but no radius; without a radius,
isolation cannot become a bound on $|\mathcal E|$. Here is the radius.

**The bound.** Along a line $P_0+rd$, $\|d\|_F=1$, the active block is
$M(r)=P_{0,TT}+r\,d_{TT}$ with $M''=0$, so second-order perturbation of a *simple*
eigenvalue gives

$$g_T(P_0+rd)=\tfrac16+r\langle G_T,d\rangle
+r^2\!\!\sum_{j\ge2}\frac{(v_j^{\mathsf T}d_{TT}v_1)^2}{\lambda_1-\lambda_j}+O(r^3).$$

Every denominator $\lambda_1-\lambda_j$ is **negative** ($\lambda_1=\tfrac16$ is the
minimum), so the quadratic term is $\le0$ — this is the familiar concavity of
$\lambda_{\min}$ along lines, and it means the correction can only be bounded, not
signed away. With $\operatorname{gap}_T:=\lambda_2(P_{0,TT})-\tfrac16>0$ and
$\|d_{TT}\|_F\le\|d\|_F=1$,

$$F(P_0+rd)\ \ge\ \tfrac16+\kappa r-\frac{r^2}{\operatorname{gap}_{\min}},
\qquad \operatorname{gap}_{\min}=\min_{T\in\mathcal A}\operatorname{gap}_T,$$

so $F>\tfrac16$ on the punctured ball of radius

$$\boxed{\;r_0=\kappa\cdot\operatorname{gap}_{\min}\;}$$

Both factors are exactly computable: $\operatorname{gap}_T$ is an exact algebraic
number and $\kappa$ comes from the LP of §4b.

| extremal | active | $\operatorname{gap}_{\min}$ | $\kappa$ | $r_0$ |
|---|---|---|---|---|
| `P(S(e,e,e),e,e,e)` | 10 | $\tfrac7{18}=0.38889$ | $0.04600$ | $0.017889$ |
| `P(S(e,e),S(e,e),e,e)` | 12 | $\tfrac13$ | $0.05727$ | $0.019091$ |
| `P(S(P(e,e,e),e),S(e,e))` | 10 | $\tfrac7{18}$ | $0.04777$ | $0.018578$ |
| `P(S(P(S(e,e),e,e),e),e)` | 12 | $\tfrac13$ | $0.06030$ | $0.020102$ |
| `P(S(P(e,e),P(e,e)),S(e,e))` | 12 | $\tfrac13$ | $0.04874$ | $0.016247$ |
| `P(S(P(e,e),e),S(P(e,e),e))` | 12 | $\tfrac13$ | $0.05401$ | $0.018002$ |
| **out-of-family** $\tfrac5{14},\tfrac9{14}$ | 13 | $\tfrac4{21}=0.19048$ | $0.05227$ | $\mathbf{0.009956}$ |

Uniform over the seven: $\kappa\ge0.0460$, $\operatorname{gap}_{\min}\ge\tfrac4{21}$,
$r_0\ge0.00996$. The gaps are clean exact rationals ($\tfrac13$, $\tfrac7{18}$,
$\tfrac4{21}$) — reassuring, and consistent with the out-of-family extremal being
the *tightest* of the seven (smallest gap and smallest radius, matching its 13
active triples: more active constraints, less room).

> **Packing bound.** Balls of radius $r_0/2$ around distinct extremals are
> disjoint, so
> $$|\mathcal E|\ \le\ \frac{\operatorname{vol}Gr(3,6)}{\operatorname{vol}B_{r_0/2}}
> \ \sim\ C\left(\tfrac{2}{r_0}\right)^{9},$$
> which with $r_0\ge0.00996$ is on the order of $5\cdot10^{20}$ — astronomically
> weak as a *count*, but **finite and explicit**, and that is the point.
>
> **Why this is the valuable direction.** The bound is conditional only on a
> *uniform lower bound for $r_0$ over all of $\mathcal E$* — not on enumerating
> $\mathcal E$. That is a strictly weaker requirement than a complete catalogue,
> and it is the only route here that escapes the sampling trap. A uniform bound
> would need $\kappa$ and $\operatorname{gap}_{\min}$ bounded below globally; both
> are semialgebraic functions of $P$, so this is a Positivstellensatz question
> again — but on a *much* smaller object than part (b), since it only concerns the
> local data at points of the level set.

---

---

## 3b. Independent confirmation of isolation by a rank argument

`verify/v16_jacobian.py`. Isolation was obtained in §1.1/§3a from the **critical
cone** — an inequality/LP argument. Here it is obtained again from a **rank**
condition, which is a different kind of argument, so agreement is a real
cross-check rather than a restatement.

In the Grassmann graph chart around $P_0$ (every nearby 3-plane is
$\operatorname{col}(U+NX)$, $X\in\mathbb R^{3\times3}$ — exactly 9 parameters), the
active equations are $g_T(X)=\lambda_{\min}(P(X)_{TT})-\tfrac16=0$. Since
$\lambda_{\min}$ is a *simple* eigenvalue at $P_0$ (verified exactly), each $g_T$ is
analytic and the Jacobian $J$ is well defined. Then

$$\operatorname{rank}J=9\ \Longrightarrow\ \text{the stratum is 0-dimensional at }P_0
\ \Longrightarrow\ P_0\text{ isolated}.$$

| extremal | $\lvert A\rvert$ | $\operatorname{rank}J$ | smallest s.v. | $\dim$ stratum |
|---|---|---|---|---|
| `P(S(e,e,e),e,e,e)` | 10 | **9** | $0.639$ | **0** |
| `P(S(e,e),S(e,e),e,e)` | 12 | **9** | $0.506$ | **0** |
| `P(S(P(e,e,e),e),S(e,e))` | 10 | **9** | $0.639$ | **0** |
| `P(S(P(S(e,e),e,e),e),e)` | 12 | **9** | $0.506$ | **0** |
| `P(S(P(e,e),P(e,e)),S(e,e))` | 12 | **9** | $0.506$ | **0** |
| `P(S(P(e,e),e),S(P(e,e),e))` | 12 | **9** | $0.506$ | **0** |
| out-of-family $\tfrac5{14},\tfrac9{14}$ | 13 | **9** | $\approx0.8$ | **0** |

**7/7 full rank, and not marginally** — the smallest singular value is $\ge0.5$, so
there is no question of a borderline rank call. No curve of extremals passes
through any known extremal.

### 3b.1 The reduction that makes the algebraic route tractable [PROVED]

At a **simple-active** equality point, $L$ has one row per active triple and 9
columns. Sharpness is stronger than $\operatorname{rank}L=9$: the active
gradients must positively span the 9-dimensional tangent space, equivalently
$0$ must lie in the relative interior of their convex hull.  A positive spanning
set in $\mathbb R^9$ has at least $10$ vectors.  Hence, in the simple-active
case,

$$\boxed{\ \text{sharp}\ \Longrightarrow\ |A|\ge10\ }$$

Contrapositive: **a simple-active extremal with $|A|\le9$ is automatically
non-sharp**.  Such a point would be a high-priority obstruction to the
sharpness-implies-finiteness route.  The dimension count is suggestive:
$|A|$ equality equations on a 9-manifold leave a tangent dimension at least
$9-|A|$ before inequalities and higher-order effects are imposed.

This gives a concrete subproblem:

> **Does any simple-active extremal have at most 9 active triples?**

Every extremal ever observed has $|A|\in\{10,12,13\}$ and simple active least
eigenvalues.  However, excluding simple-active $|A|\le9$ is not enough to prove
finiteness: a non-sharp equality point with $|A|\ge10$ could still occur if the
active gradients have rank $<9$ or fail positive spanning.  Also, if an active
least eigenvalue has multiplicity $>1$, one active triple contributes a
set-valued subgradient rather than one row of $L$; that nonsimple case must be
split off separately in any exact semialgebraic proof.  The actual finiteness
target is therefore:

> **At every equality point, the actual active gradients have rank 9 and
> positively span the tangent space.**

Proving this would make every equality point sharp and would imply finiteness by
Theorem F2.

---

## 3c. Search status, and a correction to an earlier count

**Correction.** An earlier biased-start run (`v14_biased.py`) reported "**111 new
leverage patterns**". That was **wrong** — an artifact of clustering the sorted
leverage vector on a $10^{-9}$ grid, which is *finer than the polish tolerance*, so
one genuine pattern fragments into many keys as the sample grows. Re-clustering at
$10^{-6}$: **all 111 collapse onto the four known patterns; zero genuinely new.**
Both `v12` and `v14` now cluster at $10^{-6}$ and match patterns by tolerance
rather than exact key. Distinct patterns differ by $\approx0.03$, so $10^{-6}$ is
still enormously tighter than needed.

The same artifact is visible in the currently-running MLCore job
(`gtz63-hunt-ivtoc5`, started before the fix): its pattern count creeps
$4\to5\to6\to7\to8$ purely from rounding noise as hits accumulate. Its **hits are
valid**; only the clustering is affected, and it will be re-clustered on completion.
At last check: 1800/3000 starts, **494 hits**, no non-sharp candidate, tripwire
never fired.

**Net search evidence to date:** across the pilot (300 uniform starts, 72 hits),
the biased run (981 starts spanning perturbations of known extremals,
leverage-targeted profiles, and icosahedral orbits), and 1800 starts of the MLCore
run, **exactly four leverage patterns occur and every one is sharp**. Nothing new
has appeared from either uniform or deliberately biased sampling.

---

## 3d. Equality-KKT search for non-sharp obstructions

`verify/v18_kkt_screen.py` was also run in equality mode over every nonempty
active-set orbit size.  The purpose is different from random descent: instead of
asking for more equality points, it asks whether an equality KKT point can have a
new actual active set, low active count, or non-sharp gradient geometry.

Additional equality-mode coverage:

| active sizes | orbit reps | candidates | interpretation |
|---|---:|---:|---|
| $1,\ldots,4$ | 32 | 0 | no low-active point |
| $5,6$ | 137 | 0 | no low-active point |
| $7,8$ | 410 | 1 | selected subset of a known 12-active point |
| $9$ | 312 | 0 | no equality stratum |
| $10,\ldots,13$ | 1074 | 2 | selected subsets of known 12- and 13-active points |
| $14,\ldots,20$ | 170 | 0 | no equality stratum |

`verify/v20_finiteness_probe.py` then reconstructed the projector for every
low-residual equality hit, recomputed the **actual** active set, canonicalized it
under $S_6$, and tested the rank/positive-spanning condition numerically.

Result over all equality outputs:

- input records: `2135`;
- low-residual equality hits: `3`;
- distinct actual active-set orbits: `2`;
- new actual active-set orbits: `0`;
- numerical non-sharp candidates: `0`;
- actual active size $\le9$: `0`.

The three hits are exactly:

| actual active size | rank of $L$ | positive multiplier margin | orbit |
|---:|---:|---:|---|
| 12 | 9 | $5.0\cdot10^{-2}$ | known TTSP 12-active orbit |
| 12 | 9 | $5.0\cdot10^{-2}$ | same known TTSP 12-active orbit |
| 13 | 9 | $2.22\cdot10^{-2}$ | out-of-family $(5/14,9/14)$ orbit |

This is still heuristic, but it is the right heuristic for finiteness: no search
has found the kind of equality point that could obstruct the compact-discrete
argument.

---

## 3e. Semialgebraic formulation of the finiteness obstruction

`verify/v21_semialgebraic_route.py` records the exact polynomial shape of the
next route.  Work in the Grassmann chart

$$Y=\begin{pmatrix}I_3\\ Z\end{pmatrix},\qquad
P(Z)=Y(Y^{\mathsf T}Y)^{-1}Y^{\mathsf T}=\frac{N(Z)}{d(Z)},\qquad
d=\det(Y^{\mathsf T}Y)>0.$$

For an active triple $T$, the condition
$\lambda_{\min}(P(Z)_{TT})=\tfrac16$ with simple least eigenvalue is encoded by

$$M_T(Z):=6N(Z)_{TT}-d(Z)I_3\succeq0,\qquad \det M_T(Z)=0,\qquad
\operatorname{rank}M_T=2.$$

The raw determinant has the universal factor $d^2$.  Since this chart is the
open set $d>0$, the active determinant equation can be saturated to the
degree-six equation $\det M_T/d^2=0$ without changing the represented full-rank
locus.  In algebraic computations the open condition should still be encoded,
for example by adjoining an inverse variable $u$ and the equation $u d-1=0$;
otherwise boundary components with $d=0$ may remain in the polynomial ideal.

On a cofactor patch, a nonzero kernel vector $w_T(Z)$ is the cross product of two
rows of $M_T$.  A nonzero tangent witness $h$ for non-sharpness then satisfies,
after clearing the positive denominator $d^2$,

$$
w_T^{\mathsf T}\left[
\sum_j h_j\bigl(d\,\partial_j N-N\,\partial_j d\bigr)_{TT}
\right]w_T\le0
\qquad(T\in A),
$$
together with $\sum_jh_j^2=1$.  Thus the simple-active obstruction lives in only
18 variables: 9 chart variables $z$ and 9 tangent variables $h$.  The active-set
size changes the number of equations and inequalities, not the dimension of the
ambient polynomial system.

The current degree audit gives:

| quantity | degree / upper bound |
|---|---:|
| $d$ | 6 |
| entries of $N$ | 6 |
| raw active determinant $\det M_T$ | $\le18$ |
| saturated active determinant $\det M_T/d^2$ | 6 |
| $2\times2$ PSD minors | $\le12$ |
| cofactor-patch norm | $\le24$ |
| cleared directional inequality | $\le35$ |

For one cofactor-patch branch the relaxed simple-active obstruction has:

| $|A|$ | variables | equalities | inequalities, excluding inactive disjunctions |
|---:|---:|---:|---:|
| 8 | 18 | 9 | 65 |
| 9 | 18 | 10 | 73 |
| 10 | 18 | 11 | 81 |
| 12 | 18 | 13 | 97 |
| 13 | 18 | 14 | 105 |

The inactive constraints can be relaxed away for an infeasibility proof: if the
relaxed obstruction is empty, then the true obstruction is empty.  If the relaxed
system has solutions, the inactive inequalities must be reintroduced to classify
them.

This is not yet a certificate, but it is a finite algebraic target.  The exact
proof task is now:

> Show that the semialgebraic non-sharp obstruction is empty outside the known
> sharp active-set orbits, with a separate branch for nonsimple active
> eigenvalues.

---

## 3f. Sharpened reduction: sharp ⟹ |A| ≥ 10 [PROVED]

§3b.1 gave $|A|\ge9$ from $\operatorname{rank}L\le\min(|A|,9)$. The true bound is
one higher, because sharpness needs the active gradients to **positively** span
$\mathbb R^9$, not merely span it.

> **Lemma.** A positively spanning subset of $\mathbb R^n$ has at least $n+1$
> elements.

*Proof.* If $0=\sum_T\lambda_TG_T$ with all $\lambda_T>0$ then $0$ lies in the
relative interior of $\operatorname{conv}\{G_T\}$. A set of $m$ points has convex
hull of affine dimension $\le m-1$; if the $G_T$ also span $\mathbb R^n$ linearly
then that affine hull is all of $\mathbb R^n$, so $m-1\ge n$. $\square$

$$\boxed{\ \text{sharp}\ \Longrightarrow\ |A|\ge10\ }$$

Verified numerically in dimensions $2,3,4$ (smallest positively spanning sets found:
$3,4,5$ — matching $n+1$ exactly). The observed active sizes are $10,12,13$, so the
minimum **10 is exactly tight** against this bound. The decisive question sharpens to:

> **Does any extremal have $|A|\le9$?**

*Conditional corollary.* If GTZ(6,3) is true then every point of $\mathcal E$ is a
global minimizer, hence satisfies KKT ($0\in\operatorname{conv}\{G_T\}$). So an
extremal with $|A|\le9$ would either break finiteness **or** be a non-minimizing
point of the level set — informative either way.

---

## 3g. A tolerance artifact that produced three false "non-sharp" flags

**Recorded because it nearly became a wrong headline.** The completed 3000-start
MLCore run (`gtz63-hunt-ivtoc5`, 866 hits) reported *"3 non-sharp candidates —
finiteness of E is IN DOUBT"*. All three were **artifacts**, and the mechanism is
worth knowing:

Hits are polished only to $|F-\tfrac16|<10^{-9}$, but the screen declared a triple
active at a **fixed** $10^{-9}$ window. A genuinely active triple sitting at
deviation $\sim6\cdot10^{-9}$ therefore fell *outside* the window; dropping its row
from $L$ destroyed positive spanning and produced a spurious rank/$\kappa$ failure.

Reproduced locally and confirmed: a hit with $|A|=9$ at tolerance $10^{-9}$ but
$|A|=10$ at $10^{-6}$, the tenth triple at deviation $6.6\cdot10^{-9}$. Two of the
three flags also had the tell-tale internal contradiction
$\operatorname{rank}L=9$ **and** LP optimum $0$ (which together *prove* $\kappa>0$)
alongside a reported $\kappa<0$ — impossible, and a signal that the estimator, not
the geometry, was at fault.

**Fixes applied to `v12`:** (i) the active set is now cut at the largest
*multiplicative gap* in the sorted deviations rather than a fixed threshold (the
gap is $\sim10^4$–$10^{15}$ wide in practice, so the cut is unambiguous, and the
ratio is reported so a marginal call is visible); (ii) the $\kappa$ estimator
guards against $\|z\|\to0$ excursions; (iii) a low actual active count is now
labelled `LOW-ACTIVE (check tolerance before believing)` rather than
`NON-SHARP CANDIDATE`, and an inconsistency between the LP verdict and the $\kappa$
estimate is reported explicitly as `kappa_inconsistent`. After the fix, both the
exact extremal and a perturbation mimicking a $10^{-9}$-polished hit recover
$|A|=10$ and screen SHARP.

**Net effect on the evidence: unchanged and slightly strengthened.** Across the
pilot, the biased run, and the 3000-start MLCore run (866 hits), there is still
**no genuine non-sharp candidate and no extremal with $|A|\le9$**. The independent
`v18`/`v20` KKT screens over all active-set orbit sizes $1..20$ reach the same
conclusion from a different direction: of 2135 records, the only low-residual
equality hits were **selected subsets** of known sharp orbits, with actual active
sizes $12,12,13$ — never $\le9$.

---

---

## 3h. The |A| ≤ 9 question reduced to a finite list of 124 orbits [PROVED]

`verify/v22_low_active_exclusion.py`. Exact combinatorics, no sampling.

By §3f, sharp $\Rightarrow|A|\ge10$, so **every** active set of size $\le9$ is
automatically non-sharp and is therefore exactly what could break finiteness. The
question "does any extremal have $|A|\le9$?" is a priori a search over
$\sum_{k\le9}\binom{20}{k}=431{,}909$ subsets. Two exact necessary conditions and
the $S_6$ symmetry collapse it to a short list.

**Necessary condition (a): $A$ must cover $[6]$.** If some index $i$ lies in no
active triple, the tangent directions that only move row $i$ are unconstrained by
the active set, so $\operatorname{rank}L<9$ and the point cannot be an equality
point of the expected codimension.

**Necessary condition (b): $A$ must cover all 15 pairs.** The off-diagonal entry
$c_{ij}$ enters the active constraints only through triples containing *both* $i$
and $j$; an uncovered pair leaves that freedom unconstrained.

Enumerating $S_6$-orbits (orbit counts independently reproduce the Burnside numbers
recorded in `v17`):

| $\lvert A\rvert$ | subsets | $S_6$-orbits | cover $[6]$ | cover $[6]$ + all 15 pairs |
|---:|---:|---:|---:|---:|
| 1 | 20 | 1 | 0 | 0 |
| 2 | 190 | 3 | 1 | 0 |
| 3 | 1 140 | 7 | 3 | 0 |
| 4 | 4 845 | 21 | 15 | 0 |
| 5 | 15 504 | 43 | 37 | 0 |
| 6 | 38 760 | 94 | 88 | **1** |
| 7 | 77 520 | 161 | 157 | **7** |
| 8 | 125 970 | 249 | 247 | **32** |
| 9 | 167 960 | 312 | 311 | **84** |
| **total** | **431 909** | **891** | 859 | **124** |

> **Result: exactly 124 orbits survive both elementary conditions** — a reduction
> of more than $3400\times$ from the raw subset count, and $7\times$ from the orbit
> count. **No orbit with $|A|\le5$ survives at all.** The minimum surviving size is
> $|A|=6$, realised by a *unique* orbit:
> $$\{012,\ 013,\ 045,\ 145,\ 234,\ 235\}$$
> — a strikingly symmetric configuration: three pairs $\{0,1\},\{4,5\},\{2,3\}$ each
> appearing in two triples, i.e. a perfect matching structure reminiscent of the PMC.

**Why this matters.** A full algebraic proof of finiteness no longer needs a
$\forall$-argument over the whole 9-manifold for the low-active case. It needs only
that **each of 124 explicit strata is empty**:

$$\mathcal E_A=\{P: P^2=P,\ \operatorname{tr}P=3,\ \lambda_{\min}(P_{TT})=\tfrac16
\ \forall T\in A,\ \lambda_{\min}(P_{TT})<\tfrac16\ \forall T\notin A\}=\varnothing .$$

That is a finite, explicit, *checkable* list — exactly the shape Pavel's original
"finite number of algebraic checks" hope required. Each stratum is a system in the
9 chart variables with $|A|$ equalities and the inactive strict inequalities; the
Sage/Singular infrastructure in `code/sage/` (saturated determinant ideals,
dimension probes, lex slice certificates) is already set up to attack them, and the
per-stratum systems here are *smaller* than the known-base system already handled
there ($|A|\le9$ versus 10–13 equalities).

**Caveats on scope — three, all real.**

1. Conditions (a) and (b) are *necessary*, not sufficient; the 124 survivors are
   candidates, not extremals. Emptiness of each still has to be proved.
2. The $|A|\ge10$ non-sharp case is untouched here: a point with rank-deficient
   $L$, or rank 9 without positive spanning, would also break sharpness and is not
   excluded by any size bound.
3. **The $|A|\ge10$ bound assumes simple $\lambda_{\min}$ on every active block.**
   Where $\lambda_{\min}(P_{TT})$ is a repeated eigenvalue, $F'$ is not a max of
   linear functionals but a min-eigenvalue of a compression, the subdifferential is
   an SDP-representable set rather than a polytope, and the counting argument does
   not apply. `verify/v21_semialgebraic_route.py` flags this correctly: the
   nonsimple branch needs a separate subgradient/SDP-style semialgebraic treatment.
   All seven certified extremals have simple $\lambda_{\min}$, so the assumption
   holds where we have checked — but it is an assumption, not a theorem, for the
   rest of $\mathcal E$.

---

---

## 3i. Attacking the 124 strata: the |A|=6 orbit, and why dimension is not the test

`verify/v23_stratum_sweep.py`, `verify/v24_low_active_feasibility.py`, plus the
Sage/Singular machinery in `code/sage/`.

### The unique |A|=6 orbit has a 3-dimensional equality locus [modular]

Target: $A=\{012,013,045,145,234,235\}$ (indices $0,1,9,15,16,17$). The saturated
determinant ideal with $d$ inverted (10 variables, 7 equations, max degree 7) over
$\mathbb F_{32003}$, with zero-sum integer sections, seed 401:

| sections | Singular `slimgb` | dimension |
|---:|---|---:|
| 0 | timeout at 900 s | — |
| 3 | 140 s | **0** |
| 4 | 0.4 s | $-1$ |
| 5 | 0.0 s | $-1$ |

So the equality locus of this stratum has **top dimension exactly 3** — it is
nonempty, a 3-fold. Exactly as `code/sage/RESULTS.md` predicted for the known-base
system: **determinant equalities alone cannot prove emptiness.**

> **Methodological point worth stating plainly.** "Dimension 0" is *not* the success
> criterion for these strata, and "dimension $>0$" is not failure. What finiteness
> needs is emptiness of the *semialgebraic* set, and the determinant ideal is only
> its equality part. The three genuinely distinct outcomes are:
> $\dim=-1$ (stratum empty — a proof, if in characteristic 0); $\dim\ge0$ with no
> real point satisfying the PSD/inactive inequalities (empty, but only the
> inequality layer can certify it); and $\dim\ge0$ *with* a feasible real point (a
> genuine low-active extremal, which would break the route). `v23` reports these as
> three separate buckets and never conflates them. Note also the asymmetry: a
> modular $\dim=-1$ does not prove emptiness over $\mathbb Q$ (bad primes), so only
> a characteristic-0 $-1$ counts as proof.

### The inequality layer: no feasible point found [NUMERICALLY SUPPORTED]

`v24` minimizes a penalty that vanishes *exactly* on the feasible set of a stratum
(squared residuals on the active equalities, hinges on the inactive strict
inequalities), from many random starts over the 9-dim chart.

Design point that matters: every numerical zero is **re-validated** by recomputing
the *actual* active set from scratch (cut at the largest multiplicative gap, never a
fixed tolerance) and comparing to the requested orbit up to $S_6$. Without this the
optimizer cheats by drifting onto a known 10-active extremal that happens to contain
the requested subset — the same trap that produced the spurious "8-active" hit in
the `v18` screens and the false non-sharp flags of §3g.

Results so far (sizes 6 and 7 complete, 8 and 9 in progress):

| $\lvert A\rvert$ | orbits probed | starts each | feasible points found | best penalty |
|---:|---:|---:|---:|---|
| 6 | 1 | 60 | **0** | $2.1\cdot10^{-3}$ |
| 7 | 7 | 40 | **0** | $2.1\cdot10^{-3}$ – $6.9\cdot10^{-3}$ |
| 8 | 32 (in progress) | 40 | 0 so far | $\ge4.4\cdot10^{-3}$ |
| 9 | 84 (on MLCore) | 60 | 0 so far | $\ge4.4\cdot10^{-3}$ |

The penalties bottom out around $10^{-3}$–$10^{-2}$ and **never approach zero** —
three or four orders of magnitude away from the $10^{-14}$ that a genuine feasible
point produces. This is a much cleaner separation than a marginal near-miss would
give, and it is consistent across every orbit tried.

**Interpretation, carefully.** For the $|A|=6$ orbit we now have both halves of the
picture: the equality locus is a nonempty 3-fold, and no point on it appears to
satisfy the inequalities. That is precisely the configuration in which the stratum
is empty but only the *inequality* layer can prove it — so the exact-algebra target
is a Positivstellensatz/CAD certificate on a 3-dimensional variety, not a dimension
computation.

### 3i.1 NEGATIVE FINDING: the lex-slice route does not scale to low-active strata

The natural tool was the lex-slice technique of `code/sage/RESULTS.md` (triangular
basis, real root isolation, then a uniform active-PSD sign certificate modulo the
residual factor). I expected the low-active systems to be *easier* than the
known-base one, having fewer equations. **They are dramatically harder.**

Lex Gröbner on the $|A|=6$ orbit with 3 zero-sum sections, `lex_last_variable = z0`:

| system | equations | lex quotient dim | univariate degree | over $\mathbb Q$? |
|---|---:|---:|---:|---|
| known base ($\lvert A\rvert=10$, *has* a real extremal) | 10 | 18 | 18 | yes, ~60 s |
| **$\lvert A\rvert=6$ orbit** (no extremal expected) | 6 | **2040** | **1880** | **timeout at 900 s** |

Over $\mathbb F_{32003}$ the $|A|=6$ lex basis computes in 213 s and yields a
degree-**1880** univariate polynomial in $z_0$; over $\mathbb Q$, `modslimgb` times
out.

**Why fewer constraints make it harder — the mechanism is worth stating.** Ten active
equalities on the 9-dimensional chart make the known-base system *nearly
zero-dimensional already*, so a section cuts out only 18 points. Six equalities leave
a 3-fold; the three sections do cut it to dimension 0, but to a scheme of degree
$\approx2040$. **Fewer active constraints $\Rightarrow$ less rigid $\Rightarrow$
higher section degree.** The low-active strata — precisely the ones that threaten
finiteness — are the *worst* case for this method, not the best.

Consequence: exact real-root isolation at degree 1880 over $\mathbb Q$ is out of
reach in reasonable time, and there are 124 such strata. **The lex route as it stands
cannot close the low-active case.** What it can still do is what it did here: size
the problem and rule out the naive plan before it consumes a week of compute.

Three alternatives now look better than pushing lex further:

1. **Exploit the inequalities *before* eliminating.** The penalty landscape (§3i)
   suggests the active-PSD constraints are violated by a wide margin ($10^{-3}$, not
   $10^{-14}$) everywhere on the equality 3-fold. A *sum-of-squares* certificate for
   "some active $2\times2$ principal minor is strictly negative on the whole
   equality locus" would be a degree-bounded object independent of the 1880, and the
   observed margin suggests a low-degree certificate should exist.
2. ~~**Use the necessary conditions harder.**~~ **TRIED AND FAILED — see §3i.2.**
3. **Interval / branch-and-bound on the 3-fold.** With `python-flint`/Arb already in
   place (`v19`), a rigorous exclusion sweep over the *3-dimensional* equality locus
   is a far smaller object than the original 9-dimensional problem, and the margin
   in (1) is exactly what makes interval arithmetic effective.

---

---

## 3i.2 SECOND NEGATIVE FINDING: leverage/trace conditions kill nothing

`verify/v25_stronger_conditions.py`. Tried, because it was the cheapest of the three
alternatives: strengthen the combinatorial screen so fewer strata reach the
expensive algebra.

Conditions applied to all 124 survivors, all exact:

- **coarse degree bound** $3\max_i\deg_A(i)\ge|A|/2$, from
  $\sum_i\deg_A(i)\,\ell_i\ge|A|/2$ (each active trace is $\ge\tfrac12$) together
  with $\sum_i\ell_i=3$;
- **exact leverage LP**: is the polytope
  $\{\ell:\ \ell_i\in[\tfrac16,1],\ \sum_i\ell_i=3,\ \sum_{i\in T}\ell_i\ge\tfrac12
  \ \forall T\in A\}$ nonempty? (The bound $\ell_i\ge\tfrac16$ holds because a
  diagonal entry of a PSD matrix dominates its least eigenvalue.)

> **Result: 0 of 124 orbits killed.** Both conditions are satisfied by every
> survivor, with enormous slack. For the $|A|=6$ orbit the minimum achievable
> $\sum_i\deg_A(i)\ell_i$ is $9$ against a requirement of only $3$ — **slack 6**.
> For the seven $|A|=7$ orbits: $19/2$ against $7/2$, again slack $6$.

**Why this fails, and it is instructive.** Trace/leverage conditions only see the
*diagonal* of $P$. The active constraint $\lambda_{\min}(P_{TT})=\tfrac16$ is
overwhelmingly a statement about the *off-diagonal* structure — the same reason the
single-pair second-moment program was provably dead (Theorem BND) and the reason
unsigned $\tau$-averaging cannot certify the slice (Prop. 12 of
`slice-framework.md`). This project has now hit that wall a third time, from a new
direction. **Recording it as a pattern: any argument that reduces GTZ(6,3) to
leverage data alone is too weak, and this should be the default prior on such
attempts.**

The survivors' degree profiles are also informative — they are all fairly balanced
($(6,5,5,4,4,3)$ and $(5,5,5,4,4,4)$ dominate, 42 of 124 between them), so there is
no "lopsided" subfamily to attack separately either.

That leaves alternatives 1 (SOS on the wide inequality margin) and 3 (interval
branch-and-bound on the 3-fold) from §3i.1, both of which use the *inequality*
structure rather than the diagonal — which is exactly where the information is.

---

---

## 3i.3 THE MECHANISM: the active set always grows [NUMERICALLY SUPPORTED]

`verify/v26_margin_certificate.py`, `verify/v27_obstruction_split.py`.

The two negative findings above (§3i.1 lex degree blow-up, §3i.2 leverage conditions
useless) were worth the cost, because eliminating them forced the right question:
*which* of the two possible obstructions actually makes the low-active strata empty?

There are only two candidates. On the active-equality locus of an orbit $A$:

- **(a) active PSD fails** — some active $2\times2$ principal minor of
  $P_{TT}-\tfrac16I$ goes negative, so $\lambda_{\min}=\tfrac16$ cannot actually be
  the *minimum*;
- **(b) inactive excess** — some triple *outside* $A$ also reaches $\tfrac16$, so the
  actual active set is strictly larger than $A$.

I expected (a), and built `v26` to measure its margin for an SOS certificate.
**(a) is false.** On the $|A|=6$ locus, reached to equality residual $9\cdot10^{-11}$
from 35 of 40 starts, the *worst* active $2\times2$ minor is **positive**
($+1.1\cdot10^{-3}$): active PSD holds comfortably. `v26`'s margin came out
**negative** ($\mu=-1.2\cdot10^{-2}$), i.e. there is no PSD violation to certify.

The obstruction is **(b)**, and it is uniform and large:

| orbit | on-locus | active PSD slack | **inactive excess** | verdict |
|---|---:|---:|---:|---|
| $\lvert A\rvert=6$, canon 78593 | 35/40 | $+1.08\cdot10^{-3}$ | $+5.89\cdot10^{-2}$ | inactive |
| $\lvert A\rvert=7$, canon 13105 | 28/40 | $+2.45\cdot10^{-3}$ | $+6.60\cdot10^{-2}$ | inactive |
| $\lvert A\rvert=7$, canon 14113 | 30/40 | $+3.67\cdot10^{-5}$ | $+5.89\cdot10^{-2}$ | inactive |
| $\lvert A\rvert=7$, canon 78595 | 32/40 | $+5.91\cdot10^{-4}$ | $+6.74\cdot10^{-2}$ | inactive |

(Two $|A|=7$ orbits never reached the locus at all — their equality systems appear to
have no real solution in the chart, an even stronger form of emptiness.)

> **The mechanism, stated plainly.** Impose $\lambda_{\min}(P_{TT})=\tfrac16$ on a
> small set $A$ of triples. The remaining freedom is large (a 3-fold for $|A|=6$),
> but you cannot use it to keep the other 14 triples *below* $\tfrac16$: some outside
> triple always overshoots, by about $6\cdot10^{-2}$ — more than a third of the
> threshold itself. **The active set cannot stay small; it always grows.**
>
> This is exactly what should be true if the extremals are the rigid, highly-active
> configurations we have certified ($|A|=10,12,13$), and it explains *why* no search
> has ever produced a low-active extremal.

**Why this is the right certificate target.** The statement to prove is now one
polynomial inequality per orbit —

$$\text{on }\{ \det(6N_{TT}-dI)=0,\ T\in A\}:\quad
\max_{T\notin A}\ \lambda_{\min}(P_{TT})\ \ge\ \tfrac16+c,\qquad c\approx6\cdot10^{-2},$$

with a margin of $6\cdot10^{-2}$. That needs **no root isolation** (so the degree-1880
wall of §3i.1 is irrelevant) and **no leverage argument** (so §3i.2's wall is
irrelevant too). It is a Positivstellensatz statement whose certificate degree is
governed by the margin, and the margin is a third of the threshold.

**Honest limits.** This is numerical: 40 starts per orbit, sizes 6 and 7 only, and
"never reached the locus" for two orbits is suggestive but not proof of an empty
real locus. Sizes 8 and 9 (116 orbits) are still to be classified, and the exact
certificate for even one orbit is not yet attempted. What has changed is that we now
know *which* inequality to certify, which the previous two attempts did not.

---

---

## 3i.4 Sizes 8 and 9: the same mechanism, and a claim I got wrong

**CORRECTION FIRST.** An earlier draft of this section claimed that at
$|A|\ge8$ the equality locus has *no real points at all*, based on the first 7
size-8 orbits, all of which failed to reach the locus. **That was wrong.** With the
full size-8 sweep complete (32 orbits), **22 of 32 do reach the locus**, and every
one of them shows the *same* inactive-overshoot mechanism as sizes 6 and 7. The
early run of "never reached" orbits was an ordering accident, not a size effect.

The two-mechanism story is therefore withdrawn. There is **one** mechanism.

### Complete size-8 sweep (`verify/v27_obstruction_split.py --sizes 8`)

| classification | count |
|---|---:|
| **inactive overshoot** | **22** |
| active PSD failure | 0 |
| both | 0 |
| neither (would be a feasible point) | **0** |
| never reached the locus | 10 |

Quantitatively, over the 22 orbits that reach their locus:

| quantity | min | median | max |
|---|---:|---:|---:|
| **inactive excess** | $+4.66\cdot10^{-2}$ | $+6.91\cdot10^{-2}$ | $+1.23\cdot10^{-1}$ |
| active PSD slack | $+1.61\cdot10^{-6}$ | $+3.46\cdot10^{-4}$ | $+1.89\cdot10^{-2}$ |

Active PSD slack is **positive throughout** — never the obstruction. Inactive excess
is **positive throughout**, with a median of $6.9\cdot10^{-2}$, i.e. an outside triple
overshoots the threshold by more than 40% of $\tfrac16$ itself. Size 9 (sampled so
far) matches: every orbit reaching its locus is classified `inactive`, with excesses
$+6.3\cdot10^{-2}$ and $+1.20\cdot10^{-1}$.

Note the actual active sizes landed on are $\{8,9\}$ — the optimizer stays at or just
above the requested size rather than collapsing onto a known 10+-active extremal, so
these are genuine probes of the intended strata, not drift.

### What "never reached the locus" does and does not mean

For the 10 size-8 orbits (and ~21 of the size-9 sample) that never reached their
locus, residual floors sit at $10^{-5}$–$2\cdot10^{-2}$. That is *suggestive* of real
emptiness but **not** established — it may equally be optimizer difficulty on a
thin or awkwardly-embedded real locus. Two modular dimension checks confirm the
*complex* loci are nonempty (canon 2028: dim 1 with 2 sections, 56 s; canon 4041:
dim 0, 130 s over $\mathbb F_{32003}$), so if these strata are empty it is a
genuinely *real* phenomenon — but the numerics alone cannot distinguish "no real
points" from "hard to find real points". Recorded as open.

### Consolidated status of the 124 orbits

| $\lvert A\rvert$ | orbits | swept | reached locus | inactive mechanism | feasible points |
|---:|---:|---|---:|---:|---:|
| 6 | 1 | yes | 1 | 1 | **0** |
| 7 | 7 | yes | 5 | 5 | **0** |
| 8 | 32 | yes | 22 | 22 | **0** |
| 9 | 84 | in progress | — | all so far | **0** so far |

Plus the independent `v24` feasibility sweep: 40 orbits (sizes 6–8) complete with
**zero** feasible points, size 9 running.

**The certificate target is unchanged and now better supported:** on the
active-equality locus of any low-active orbit,

$$\max_{T\notin A}\ \lambda_{\min}(P_{TT})\ \ge\ \tfrac16+c,\qquad
c\ \approx\ 4.7\cdot10^{-2}\ \text{(observed minimum)},$$

a single polynomial inequality per orbit with a margin of at least 4.7% absolute —
28% of the threshold. No root isolation, no leverage argument, no PSD case split.

---

## 3i.5 The certificate is a DISJUNCTION, not one inequality — correcting §3i.3

`verify/v28_growth_certificate.py`, 108 locus points on the unique $|A|=6$ orbit.

§3i.3 and §3i.4 both stated the certificate target as *"one polynomial inequality per
orbit"*. **That was too optimistic, and this probe shows why.**

**Step 1 — the witness triple is not unique.** Recording which outside triple
overshoots at each locus point:

| statistic | value |
|---|---|
| distinct witness triples over 108 points | **12+** (most frequent only 11 points) |
| outside triples above $\tfrac16$ per point | min **1**, median **3**, max **6** |
| overshoot magnitude | min $+4.45\cdot10^{-2}$, median $+1.30\cdot10^{-1}$, max $+3.21\cdot10^{-1}$ |

So the true statement is a **disjunction**: *some* outside triple overshoots, but
which one depends on where you are on the locus. At the worst points only **one**
outside triple overshoots, so the disjunction cannot be thinned to a safe subset.

**Step 2 — the polynomial form, and why no single triple certifies.** For
$M=P_{TT}-\tfrac16I$ symmetric $3\times3$, $\lambda_{\min}(M)\ge0$ iff
$\operatorname{tr}M,\ e_2(M),\ \det M$ are all $\ge0$ (Descartes on the
characteristic polynomial) — a *polynomial* criterion, no eigenvalues. Testing every
outside triple across all 108 locus points:

> **Outside triples with $\operatorname{tr},e_2,\det$ all positive at every sampled
> point: 0.**

The best candidates have $\operatorname{tr}>0$ and $e_2>0$ comfortably but
$\min\det<0$ (e.g. $(1,2,5)$: $\operatorname{tr}\ge+0.584$, $e_2\ge+0.011$,
$\det\ge-0.068$). A sign pattern $(+,+,-)$ means exactly one negative eigenvalue,
i.e. that triple is genuinely *not* overshooting at that point — consistent, not a
numerical artifact.

**Consequence for the proof effort.** Per orbit the exact statement is a disjunction
over $\le20-|A|$ outside triples (14 for $|A|=6$), so the low-active case is
$\sim124\times14\approx1700$ polynomial-inequality regions rather than 124. Still
finite and still elementary in shape — each region is "on this piece of the locus,
this triple's three characteristic-polynomial coefficients are $\ge$ margin" — but
the bookkeeping is an order of magnitude larger than §3i.3 implied, and a
Positivstellensatz certificate for a disjunction requires a covering argument
(CAD-style), not a single SOS identity.

**What survives from §3i.3 unchanged:** the *mechanism* (active sets always grow),
its uniformity (every orbit reaching its locus is classified `inactive`; active PSD
never fails), and the margin ($\ge4.45\cdot10^{-2}$, i.e. 27% of the threshold).
Those are the substantive findings. What was wrong was my estimate of the
certificate's *shape*, and it was wrong in the direction of optimism — recorded here
rather than quietly patched, per the brief's §8 rule 1.

---

---

## 3i.6 All 124 orbits classified; the margin survives a harder test

### The complete sweep

`verify/v27_obstruction_split.py`, sizes 6–9. Size 9 ran to completion on MLCore
(job `gtz63-split9-twj4p9`, preemptible, SUCCEEDED).

| $\lvert A\rvert$ | orbits | reached locus | **inactive** | active PSD | both | **neither (= feasible!)** | off-locus |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 6 | 1 | 1 | 1 | 0 | 0 | **0** | 0 |
| 7 | 7 | 5 | 5 | 0 | 0 | **0** | 2 |
| 8 | 32 | 22 | 22 | 0 | 0 | **0** | 10 |
| 9 | 84 | 42 | 42 | 0 | 0 | **0** | 42 |
| **total** | **124** | **70** | **70** | **0** | **0** | **0** | **54** |

> **Every one of the 70 orbits that reaches its active-equality locus is obstructed
> by inactive overshoot. Not one is obstructed by active PSD, and not one admits a
> feasible point.** The remaining 54 never reach their locus numerically at all
> (residual floors $10^{-5}$–$2\cdot10^{-2}$), which is suggestive of real-emptiness
> but, as noted in §3i.4, not established.

Overshoot magnitudes by size (min / median / max):

| $\lvert A\rvert$ | min | median | max |
|---:|---:|---:|---:|
| 8 | $+4.66\cdot10^{-2}$ | $+6.91\cdot10^{-2}$ | $+1.23\cdot10^{-1}$ |
| 9 | $+5.36\cdot10^{-2}$ | $+9.23\cdot10^{-2}$ | $+1.50\cdot10^{-1}$ |

Active PSD slack is essentially never negative (size-9 minimum
$-4.4\cdot10^{-17}$ — machine zero, i.e. a block exactly on the PSD boundary, which
is what an active constraint *means*).

### A margin scare, and why it was a false alarm

Worth recording because it nearly became a retraction. `v29`'s set-cover run
reported, at $\tau=10^{-9}$, a `margin_min` of $2.2\cdot10^{-4}$ — a factor 264
below the claimed $4.5\cdot10^{-2}$, on 5× more sampled points. That pattern (minimum
collapsing as sampling grows) is the classic signature of an infimum tending to 0,
which would have refuted the margin and, if the infimum were attained, produced a
genuine low-active extremal and **broken the finiteness route**.

It was an artifact of comparing two different quantities:

- `v29` measured $\min_{\text{points}}\max_{T\in S}(\lambda_{\min}(P_{TT})-\tfrac16)$
  over the *cover subset* $S$ — which at $\tau=10^{-9}$ contained only 10 of the 14
  outside triples;
- the margin claim is about $\max$ over **all** outside triples.

A point whose best witness lay among the 4 excluded triples therefore reported a
small value — a property of the cover, not of the locus. Confirmed by the trend: as
$\tau$ rises the cover is forced to include the strong witnesses and `margin_min`
climbs monotonically, $2.2\cdot10^{-4}\to2.7\cdot10^{-2}\to4.9\cdot10^{-2}$.

**Direct test.** `verify/v30_margin_infimum.py` minimises
$g(P)=\max_{T\notin A}(\lambda_{\min}(P_{TT})-\tfrac16)$ *on* the locus, with an
annealed penalty ($w=10^3\to10^5\to10^6$) and re-validation of every hit by the
gap rule. The first hard run on the unique $|A|=6$ orbit showed a plateau near
$3\cdot10^{-2}$, not a collapsing infimum. The broader stress test is now stronger:
unique size 6, all five reachable size-7 orbits, and the first four size-8 orbits
with smallest `v27` inactive excess all return **MARGIN-HOLDS**.

Generated by `verify/v33_margin_table.py`:

| orbit | starts | locus hits | $\min g$ | $F$ at best point | outside max-tie |
|---|---:|---:|---:|---:|---|
| $(6,78593)$ | 150 | 76 | $3.0116\cdot10^{-2}$ | 0.196782492286 | 024,034,125,135,245 |
| $(7,13105)$ | 60 | 35 | $4.4646\cdot10^{-2}$ | 0.211312502656 | 045,123,245 |
| $(7,14113)$ | 80 | 41 | $2.1959\cdot10^{-2}$ | 0.188625929570 | 123,124,125 |
| $(7,78595)$ | 60 | 25 | $4.5213\cdot10^{-2}$ | 0.211879783913 | 023,134,345 |
| $(7,78601)$ | 60 | 31 | $5.1506\cdot10^{-2}$ | 0.218172650287 | 023,034,135 |
| $(7,78612)$ | 60 | 35 | $4.2603\cdot10^{-2}$ | 0.209269971872 | 014,034,123,345 |
| $(8,13107)$ | 60 | 3 | $9.8207\cdot10^{-2}$ | 0.264873900998 | 124,134 |
| $(8,14117)$ | 60 | 24 | $3.2938\cdot10^{-2}$ | 0.199604437889 | 024,025 |
| $(8,78613)$ | 60 | 21 | $5.2998\cdot10^{-2}$ | 0.219664501050 | 135,145,235,245 |
| $(8,79656)$ | 60 | 31 | $5.2567\cdot10^{-2}$ | 0.219233269187 | 123,124,145 |

So the old $4.5\cdot10^{-2}$ quoted margin was too optimistic, and even the
$3\cdot10^{-2}$ figure from the $|A|=6$ orbit is no longer the smallest observed
value. The current tested minimum is $2.1959\cdot10^{-2}$ at size-7 canon 14113
(13% of the threshold). The mechanism itself still stands: no run found
$g<10^{-8}$, and every best point has $F>\tfrac16$ because an outside triple
overshoots.

**KKT strengthening.** The tightest minimizers are not just random optimizer
outputs; they satisfy the nonsmooth KKT balance for the margin problem. The
diagnostic `verify/v34_margin_kkt.py` works in the standard chart and solves for
active equality multipliers plus a convex combination of the outside max-tie
gradients.

| orbit | outside max-tie | KKT residual | smallest singular value | min outside weight |
|---|---|---:|---:|---:|
| $(7,14113)$ | 123,124,125 | $5.37\cdot10^{-16}$ | $1.22\cdot10^{-1}$ | 0.330 |
| $(8,14117)$ | 024,025 | $5.29\cdot10^{-16}$ | $1.21\cdot10^{-1}$ | 0.499 |

The weights are strictly positive, so these look like honest local minimizers of
the low-active margin objective. This is useful structurally: it identifies the
right outside tie strata and rules out the interpretation that the optimizer is
just wandering on a flat numerical artifact.

There was also a useful negative exact probe. `code/sage/probe_margin_tie_system.py`
builds the active-plus-outside-tie determinant ideal with a variable $q=6\lambda$
for the tied outside eigenvalue. The shifted determinant has the same universal
$d^2$ factor as the active determinant, so the saturated tie equations have degree
9 rather than the raw degree 21. Even after this reduction, the bare determinant
tie system for the tightest orbit $(7,14113)$ with ties 123,124,125 timed out over
$\mathbb F_{32003}$ after 300 seconds, and 1-, 2-, and 3-section probes also timed
out after 180 seconds. Interpretation: tie determinants alone define too broad a
variety for direct Groebner treatment.

The direct determinant-gradient KKT probe was then implemented as
`code/sage/probe_margin_kkt_determinant.py`. It builds the square system

$$
f_T=0,\quad h_S=0,\quad u d-1=0,\quad
\sum_T a_T\nabla_z f_T+\sum_S b_S\nabla_z h_S=0,\quad
1+\sum_S b_S\partial_q h_S=0.
$$

Over $\mathbb F_{32003}$ the raw systems still time out:

| orbit / tie set | variables | equations | timeout |
|---|---:|---:|---:|
| $(7,14113)$, tie 123 | 19 | 19 | 240s |
| $(7,14113)$, ties 123,124 | 20 | 20 | 240s |
| $(7,14113)$, ties 123,124,125 | 21 | 21 | 300s |
| $(8,14117)$, ties 024,025 | 21 | 21 | 240s |

Putting multiplier variables first did not change the two-tie result. A conceptual
warning explains part of this: when $|A|+|B|=10$, the active-plus-tie equations are
already square in $(z,q)$, so at regular isolated points the multiplier KKT block is
essentially automatic. Thus KKT is a good diagnostic and stratum selector, but not
yet a computational shortcut for the hardest triple-tie exact solve. The next exact
attempt should exploit the linear multiplier block explicitly (rank/minor
conditions for lower-tie strata) or attack the active-plus-tie zero-dimensional
system with a more specialized solver, rather than handing the full KKT system to
plain `slimgb`.

**Specialized solver attempt.** The multiplier-free rank-minor route was tried as
`code/sage/probe_margin_kkt_rank.py`, but even exporting the two-tie rank determinant
for $(7,14113)$ did not finish in several minutes. This confirms the degree problem:
eliminating multipliers naively creates minors of degree roughly 40--50, worse than
the linear-multiplier system.

The useful specialized object is instead the square active-plus-tie system itself.
`code/sage/refine_margin_tie_root.py` uses exact determinant polynomials and their
exact symbolic Jacobian, then runs high-precision Newton from the `v30` minimizer.
It refines the intended eigenvalue branch cleanly:

| orbit / ties | precision | residual | cond. of Jacobian | refined margin |
|---|---:|---:|---:|---:|
| $(7,14113)$, 123,124,125 | 500 bits | $1.37\cdot10^{-148}$ | 5.53 | 0.0219593612977106471453390091112 |
| $(8,14117)$, 024,025 | 400 bits | $8.98\cdot10^{-119}$ | 6.02 | 0.0329381178422993768053283604786 |

Both refined roots remain on the intended branch: the active triples are at
$1/6$, the listed outside triples tie at $q/6$, and all other outside triples are
below. Low-height `algdep` screens found no relation for $q$ (degree up to 16 for
$(7,14113)$ and up to 14 for $(8,14117)$, height $\le10^6$), so these are probably
not simple radical points. The route is therefore: use numerical algebraic
refinement to isolate well-conditioned candidate roots, then certify them by
interval/Krawczyk or rational-univariate methods, not by brute-force Groebner on
the whole KKT ideal.

**Krawczyk certification.** The interval step was implemented as
`code/sage/certify_margin_tie_krawczyk.py`. For a radius-$10^{-30}$ box around the
refined point it proves $K(X)\subset\operatorname{int}(X)$ for the square
active-plus-tie system, hence existence and uniqueness of a root in the box. It
then checks the spectral branch on the certified root enclosure $K(X)$: $d>0$,
$q>1$, all active/tie shifted $1\times1$ and $2\times2$ principal minors are
positive, and every other outside triple has $\det(6P_{TT}-qI)<0$.

The batch driver `code/sage/batch_margin_tie_certify.py` now certifies all six
non-over-tied tested size-7/8 minimizers:

| orbit / equation ties | Krawczyk radius | certified margin interval |
|---|---:|---:|
| $(7,13105)$, 045,123,245 | $10^{-30}$ | $[0.04464692543123093116213202070744,\ 0.04464692543123093116213202070777]$ |
| $(7,14113)$, 123,124,125 | $10^{-30}$ | $[0.02195936129771064714533900911105,\ 0.02195936129771064714533900911138]$ |
| $(7,78595)$, 023,134,345 | $10^{-30}$ | $[0.04521367233867756678261877566381,\ 0.04521367233867756678261877566414]$ |
| $(7,78601)$, 023,034,135 | $10^{-30}$ | $[0.05150617318963458614354853172380,\ 0.05150617318963458614354853172414]$ |
| $(8,13107)$, 124,134 | $10^{-30}$ | $[0.09820775013166748461073489955964,\ 0.09820775013166748461073489955997]$ |
| $(8,14117)$, 024,025 | $10^{-30}$ | $[0.03293811784229937680532836047841,\ 0.03293811784229937680532836047875]$ |

This is the first actual interval certificate layer in the low-active margin
route. It does not prove the whole finite low-active list, but it proves that six
numerically observed local minimizers are genuine, isolated, correctly ordered
active-plus-tie roots with strictly positive margin.

The remaining tested low-active margin points are exactly the over-tied ones:
$(6,78593)$ with five outside ties, $(7,78612)$ with four, $(8,78613)$ with four,
and $(8,79656)$ with three. Square subsystems either have singular/ill-conditioned
Jacobians or leave the unselected tie as a tiny interval around zero (upper bounds
about $10^{-54}$ after evaluating on $K(X)$); the latter is not a proof of
equality. A direct `QQ` linear-span test also shows that the extra tie determinant
is not a global linear combination of the selected determinant equations
(`code/sage/probe_extra_tie_dependency.py`,
`code/sage/out/extra_tie_dependency_regular_overtied.json`). These cases need an
overdetermined/local-membership certificate, not another arbitrary square-subset
Krawczyk run.

**Over-tied refinement.** The full over-tied systems were then refined directly
with `code/sage/refine_margin_overtie_root.py`, using all active equations and all
numerical outside ties. This cleanly separates the remaining cases:

| orbit / all ties | equations / variables | residual | min s.v. | interpretation |
|---|---:|---:|---:|---|
| $(7,78612)$, 014,034,123,345 | $11/10$ | $2.29\cdot10^{-178}$ | $42.2$ | full-rank common root |
| $(8,79656)$, 123,124,145 | $11/10$ | $1.85\cdot10^{-178}$ | $9.79$ | full-rank common root |
| $(8,78613)$, 135,145,235,245 | $12/10$ | $3.51\cdot10^{-178}$ | $72.3$ | full-rank common root, but active branch degenerate |
| $(6,78593)$, 024,034,125,135,245 | $11/10$ | $2.01\cdot10^{-25}$ | $1.85\cdot10^{-13}$ | singular/ill-conditioned |

The full all-tie Groebner probes still time out over $\mathbb F_{32003}$ after
120 seconds, so this is local numerical-algebraic structure, not a global RUR.
However, `code/sage/screen_overtie_square_bases.py` finds well-conditioned square
bases including all ties for the three full-rank size-7/8 roots, and
`code/sage/certify_overtie_square_basis.py` certifies those square bases by
Krawczyk. For $(7,78612)$ and $(8,79656)$ all spectral branch inequalities are
strict; only one omitted active determinant remains interval-enclosed at about
$10^{-54}$. For $(8,78613)$ two active determinants are omitted and an active
spectral coefficient is also only separated at the $10^{-52}$ scale, so that case
needs a deeper degeneracy split. The size-6 orbit is worse: the full Jacobian is
nearly rank-deficient, so it likely needs deflation or a different
parameterization before an interval certificate will be meaningful.

For the two clean over-tied cases, a 1400-bit refinement drives the full
overdetermined residual down to about $10^{-419}$, and a 1600-bit Krawczyk run
with radius $10^{-120}$ shrinks the omitted active determinant intervals to about
$10^{-234}$ while still containing zero. This strongly supports exact vanishing,
but it still does not prove it. The exact bounded-membership test
`code/sage/probe_bounded_ideal_membership.py` is negative for the unseparated
omitted determinant through multiplier degree $2$ over $\mathbb Q$ and through
degree $5$ modulo $32003$. The simple separated targets $d h$ and $d^2 h$, where
$d$ is the chart denominator, are also negative modulo $32003$ through degrees
$6$ and $5$ respectively. Thus the missing argument is probably not a low-degree
global ideal identity with an obvious separator.

The more flexible local-separator search
`code/sage/probe_local_separator_membership.py` is more promising. It searches
for an unknown separator $s$ satisfying

\[
  s h=\sum_i a_i f_i,\qquad \deg a_i\le D,\quad \deg s\le S.
\]

For both clean over-tied cases it is negative through $(D,S)=(5,3)$ modulo
$32003$, but at $(D,S)=(5,4)$ it finds modular separators nonzero at the refined
root: rank defect $5$ with 12 nonzero-root candidates for $(7,78612)$, and rank
defect $6$ with 21 candidates for $(8,79656)$. For $(7,78612)$ candidate 9 gives
a useful diagnostic split.  The scaled rational quartic factor is an exact
bounded member over $\mathbb Q$ at multiplier degree $5$ (rank $29920$ equals
augmented rank $29920$), but that factor evaluates to about $5.5\cdot 10^{-420}$
at the refined root and its interval evaluation contains zero.  Conversely, the
raw residue integer lift is rigorously nonzero on the certified Krawczyk box
(near $-7606.296$), but it is not an exact $\mathbb Q$ member at the same degree
(augmented rank $29921$). Thus the modular separator signal is real but not yet
a proof. The next exact step is rational reconstruction across primes or a
direct $\mathbb Q$ local-separator solve. Re-running Singular on the
well-conditioned all-ties square bases still times out over $\mathbb F_{32003}$
after 240 seconds, so plain global elimination remains unattractive even after
basis selection.

### Where the disjunction stands

`v29`'s actual purpose — shrinking the §3i.5 disjunction by set cover — **failed**,
and that is the useful part of the run. On 172 locus points of the $|A|=6$ orbit the
greedy cover needs **10 of 14** outside triples at $\tau=10^{-9}$, rising to **all
14** at $\tau=3\cdot10^{-2}$. Four points have exactly one overshooting triple, and
those forced witnesses are three *different* triples ($(1,2,5)$, $(0,1,4)$,
$(1,3,4)$), so no cover can omit them. **The disjunction is irreducible**: the exact
certificate must branch over essentially all outside triples, confirming the
~1700-region estimate of §3i.5 rather than reducing it.

*(Implementation note: the first `v29` run printed "EXACT minimum cover size -"
because `min_cover` built its target bitmask from the number of sets rather than the
number of points. Fixed; the greedy sizes quoted above were unaffected.)*

---

---

## 3i.7 The mechanism, correctly stated (and a terminology trap in my own output)

`v30`'s verdict was **MARGIN-HOLDS**, and its detailed output finally makes the
mechanism precise. The ten best points on the $|A|=6$ locus look like this:

```
g=3.0116e-02  resid=3.31e-14  F=0.196782492286  |A|actual=6  matches=True
```

**Read the $F$ column.** These points have $\lambda_{\min}(P_{TT})=\tfrac16$ on all
six $T\in A$ to residual $3\cdot10^{-14}$ — but some *outside* triple sits at
$0.1968$, so
$$F(P)=\max_T\lambda_{\min}(P_{TT})=0.1968\ \neq\ \tfrac16 .$$

**These points are not extremals at all.** They are not near-counterexamples either
— they satisfy GTZ(6,3) with 18% to spare. So:

> **The mechanism, correctly stated.** You cannot have $F=\tfrac16$ with a small
> active set, because pinning $\lambda_{\min}=\tfrac16$ on few triples leaves enough
> remaining freedom that some *other* triple is pushed strictly **above** $\tfrac16$
> — which raises the max, so the point leaves the level set $\{F=\tfrac16\}$
> entirely. Small active sets are not "hard to make extremal"; they are
> *incompatible* with being on the level set at all.
>
> Equivalently: on $\{F=\tfrac16\}$ the active set is forced to be large, because
> $\tfrac16$ has to be the **maximum**, and few constraints cannot hold the other
> 14–19 triples down.

**A terminology trap in my own diagnostic, worth flagging.** `v30` prints
`matches=True`, which reports only that the *gap-rule cluster* of triples nearest
$\tfrac16$ coincides with $A$. It does **not** mean "this is an extremal with active
set $A$" — that would additionally require $F=\tfrac16$, i.e. the maximum to be
attained at $\tfrac16$. A reader skimming for `matches=True` would badly misread
these rows. The load-bearing line is instead

```
points with g < 1e-8 : 0     ... whose ACTUAL active set is A : 0
```

i.e. **no low-active extremal was found**, which is the actual result. (Same class of
trap as the fixed-tolerance artifact of §3g: a diagnostic field that means less than
its name suggests.)

This also explains, retrospectively, why the 54 off-locus orbits fail even to reach
their equality locus: for them the incompatibility bites earlier — the six-to-nine
equalities cannot be satisfied simultaneously in the chart at all, let alone with the
max at $\tfrac16$.

---

---

## 3i.8 The 54 off-locus strata: nonempty over $\mathbb C$, so their emptiness is *real*

The 54 orbits that never reach their equality locus numerically are 44% of the 124,
so it matters whether they are empty for an algebraic reason (cheap to certify) or a
real one (harder). Singular over $\mathbb F_{32003}$, saturated determinant ideal with
$d$ inverted, zero-sum sections, seed 401:

| orbit | indices | 1 section | 2 sections | 3 sections | top dim |
|---|---|---|---|---|---:|
| $\lvert A\rvert=8$, canon 2028 | $\{0,1,2,3,4,5,6,19\}$ | — | dim 1 (56 s) | — | $\ge3$ |
| $\lvert A\rvert=8$, canon 4041 | $\{0,1,2,3,4,9,16,17\}$ | — | dim 0 (130 s) | — | **2** |
| $\lvert A\rvert=8$, canon 77640 | $\{0,1,2,3,16,17,18,19\}$ | timeout (700 s) | dim 0 (8.5 s) | dim $-1$ (0.8 s) | **2** |
| $\lvert A\rvert=8$, canon 8067 | $\{0,1,2,3,4,16,17,18\}$ | timeout (700 s) | dim 0 (59 s) | dim $-1$ (3.3 s) | **2** |

Every probed off-locus stratum has a **nonempty** determinant variety, of top
dimension 2–3. So these strata are *not* empty as algebraic sets.

> **Consequence.** For the off-locus orbits the emptiness must be a **real** (not
> algebraic) phenomenon: the complex variety is a surface or 3-fold, but it carries
> no real points satisfying the constraints — or, more precisely, no real points at
> all in the relevant chart region, which is why the optimizer cannot approach them.
>
> That is bad news for the cheapest hope of §3i.4. Real-emptiness of a variety that
> is nonempty over $\mathbb C$ cannot be certified by a Gröbner/dimension argument at
> all: it needs a real-algebraic tool (SOS certificate that a sum of squares is
> bounded below, CAD, or real root counting à la Hermite/Sturm). The dimension
> computations tell us the *shape* of the object but cannot settle it.

Note also the sharp non-monotonicity in cost: 1 section times out at 700 s while 2
sections finish in 8.5 s and 3 sections in 0.8 s. Under-sectioned systems are far
harder than either the full system or the over-sectioned one — worth knowing for
anyone budgeting these runs. (`code/sage/RESULTS.md` records the same effect for the
known-base system.)

---

---

## 4. What remains open

**F6: is $\mathcal E$ finite?** Still **OPEN**, but the shape of the remaining gap
is now precise, and it is *the same gap as part (b)*:

- To conclude finiteness from F2 we need "every point of $\mathcal E$ is sharp" —
  a statement quantified over all of $\mathcal E$.
- Verifying it at sampled points (however many) cannot close it, for the same
  reason sampling cannot close part (b): both are $\forall$-statements over a
  9-dimensional manifold.

So finiteness and the global gap are **the same difficulty in different clothing**.
That is itself informative: it says the hunt cannot succeed by accumulation, and
effort is better spent on the two concrete reductions below.

**New semialgebraic evidence.**  The Sage/Singular probes in `code/sage/` show
that the saturated determinant equations alone do not define a finite set.  For
the known 10-active base orbit, the determinant-only locus in the real chart has
a three-dimensional component even after adjoining an inverse for
$d=\det(Y^TY)$.  However, real affine section probes through the known
extremals found that the sampled spurious determinant roots fail already at the
active PSD inequalities (and hence also fail the full equality-stratum
conditions).  This suggests a sharper next target: certify that the determinant
locus plus active PSD and inactive inequalities contains only the known point in
each known active orbit, before spending effort on the non-sharp tangent
obstruction.

**Two routes that would make real progress:**

1. **Certified radius (§6(ii) of `sharp-cone-at-extremal.md`).** Make $r_0$ in
   Lemma F1 *explicit* by bounding $C$ — the Lipschitz constant of $F'$ near each
   extremal — in closed form (the blocks are $3\times3$; the eigenvalues are
   algebraic). Explicit radii turn Lemma F1 from qualitative isolation into a
   **volume/packing bound on $|\mathcal E|$**: if every extremal owns a ball of
   radius $r_0$ on which no other extremal lies, then
   $|\mathcal E|\le\operatorname{vol}(Gr(3,6))/\operatorname{vol}(B_{r_0})$ — an
   explicit finite bound, *conditional only on a uniform $r_0$*, not on knowing
   every extremal. This is the single most valuable next computation.
2. **Semialgebraic non-sharp exclusion.** Use the cofactor-patch formulation of
   §3e.  A dimension-0 or infeasibility verdict for each active-set stratum would
   prove F6 outright, with no sampling.  This is a finite computation, though
   possibly a large one — and unlike the hunt, it is not defeated by an unsampled
   extremal.

---

## 5. Reproduction

```bash
cd ~/Documents/gtz
.venv/bin/python -u verify/v9_rational_certificate.py   # 6/6 SHARP, exact
.venv/bin/python -u verify/v11_seventh_exact.py         # 7th, exact
.venv/bin/python -u verify/v12_hunt.py 300              # pilot hunt, seed 20260731
GTZ_CPUS=19 .venv/bin/python -u verify/v12_hunt.py 20000  # full hunt
```

`v12` honours `GTZ_CPUS`; it defaults to `os.cpu_count()-2`. Representatives of
any new pattern are saved to `verify/data/new_pattern_*.npy` for exact
reconstruction with `v11`'s machinery.

---

## 6. Honest summary

Pavel asked whether the extremal set is finite. The answer, precisely:

- **Finiteness is now a theorem *conditional* on sharpness** (F2), and sharpness is
  **exactly certified at seven configurations** and numerically confirmed at every
  configuration any search has found. Each certified extremal is provably an
  **isolated** point.
- **But the conditional cannot be discharged by searching.** "Every extremal is
  sharp" is a $\forall$-statement over a 9-manifold — the same logical shape as
  part (b) of Reformulation R. Finiteness and the global gap are the same
  difficulty.
- **The genuinely new structural fact** is that the corpus's extremal catalogue is
  incomplete, and not marginally so: the out-of-family extremal is the *most
  frequently found* basin in unconstrained search. Any proof strategy that
  proceeds by exhausting the nine scaled-star matrices is unsound as stated.
- **The most valuable next step is the certified radius**, because it upgrades
  isolation into an explicit *upper bound* on $|\mathcal E|$ without needing to
  know every extremal — the only route here that escapes the sampling trap.
