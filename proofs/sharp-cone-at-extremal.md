# The Nesterenko extremal is a *sharp* (first-order, conically stable) minimum of $F$

**Deliverable:** session findings for `GTZ 6 3 Closure Attempt.md`, Tracks A/C.
Scripts: `verify/v1_foundations.py`, `verify/v4_nesterenko.py`, `verify/v5_tightness.py`
(+ the shell/slope probes recorded in §4). Reproduction commands and seeds in §5.

---

## 0. Scope — read this first (overclaim guard)

> **This document does NOT prove GTZ(6,3).**
>
> It (i) re-verifies, in exact arithmetic and independently of the ingested
> corpus, the foundational lemmas of `slice-framework.md` and Lemma 2.1 of
> `boundary-obstruction.md`; (ii) reports one **correction** and one **new exact
> closed form** in the PMC material; and (iii) **PROVES exactly** that the
> Nesterenko extremal $A_0$ is a *sharp* minimum of $F$ — a nonzero first-order
> growth rate in every direction, so $F$ grows *linearly*, not quadratically,
> away from $A_0$.
>
> Item (iii) upgrades part (a) of Reformulation R (`boundary-obstruction.md` §5)
> at $A_0$ from NUMERICALLY SUPPORTED to **PROVED**, and in a stronger form than
> originally asserted: not merely a strict local minimum but a sharp one, with an
> exact rational KKT certificate (§4a), and — via a fully **rational** certificate
> that eliminates the surds entirely — at **all six** extremals of the census
> (§4b) and at a **seventh, out-of-family** extremal found by unconstrained search
> (§4c.1). Scope caveat: the census finds six where the corpus claims nine, and
> the seventh shows the extremal set is larger and less well-catalogued than the
> corpus describes (§4c). It also carries a warning for Track A: the
> linear growth means the certificate to look for is **not** a smooth SOS bound of
> the usual quadratic-vanishing type.
>
> **The global gap — part (b) of Reformulation R, that no KKT point of $F$ lies
> below $\tfrac16$ — remains OPEN and is untouched here.** That is the whole
> remaining difficulty of GTZ(6,3).

**Status table.**

| # | Statement | § | Status |
|---|---|---|---|
| 1 | Lemma 2/Lemma 3 of `slice-framework.md` hold as *generic polynomial identities*, not just at sample points | §1 | **PROVED** (exact; Groebner) |
| 2 | $J_{\rm PMC}^2=I$, row-mass, flow, $P^2=P$, $\sum_T\tau_T=\tfrac49$, $\sum_Tp_T=12$, $\sum_Tq_T=0$ | §1 | **PROVED** (exact in $\mathbb Q(s,c)$) |
| 3 | **The 6 active PMC triples are NOT algebraically uniform**: they split $4+2$ into two types | §2 | **PROVED** (exact) |
| 4 | **$G_{\rm PMC}=\frac{1-\sin 36^\circ}{2}$ exactly**; minpoly $256x^4-512x^3+304x^2-48x+1$ | §2 | **PROVED** (exact) |
| 5 | $G_{\rm PMC}>\tfrac16$ via the *rational* certificate $169<405$ | §2 | **PROVED** (exact) |
| 6 | Lemma 2.1 of `boundary-obstruction.md` (Nesterenko extremal, $F(A_0)=\tfrac16$ exactly, $\delta=0$, all $q(p)=0$) | §3 | **PROVED** (exact, re-derived) |
| 7 | $A_0$ is a **sharp** minimum of $F$ on $Gr(3,6)$: $\kappa=\min_{\|\dot P\|=1}F'(A_0;\dot P)>0$ | §4a | **PROVED** (exact; positive spanning) |
| 7b | KKT multipliers at $A_0$ are the exact rationals $(\tfrac7{90})^{\times9},\tfrac3{10}$ — all **strictly** positive | §4a | **PROVED** (exact) |
| 7c | **ALL SIX** extremals of the census are sharp; certificate is fully **rational** (no surds) | §4b | **PROVED** (exact rational LP) |
| 7d | A **SEVENTH** extremal exists OUTSIDE the scaled-star family (leverages $\tfrac5{14},\tfrac9{14}$, 13 active, entries in $\mathbb Q(\sqrt2,\sqrt5)$) — and is also sharp | §4c.1 | **PROVED** (exact) |
| 8 | No counterexample: $F\ge\tfrac16-10^{-6}$ everywhere sampled, incl. descents launched *at* $A_0$ | §4 | **NUMERICALLY SUPPORTED** |

---

## 1. Independent exact re-verification of the slice framework

Per §8 rules 1–2 of the brief, every ingested lemma used below was re-derived
rather than cited. Two layers, deliberately:

**Layer G (generic).** The identities of `slice-framework.md` Lemmas 2–3 are
established as polynomial identities in *indeterminates* $a,b,c$ (and in generic
unit vectors $u_i$), so they hold at **every** slice point, not merely at tested
samples — strictly stronger than the source's sample-based numeric blocks:

- $\operatorname{charpoly}(J_{TT})=y^3-p_Ty-2q_T$ — residual exactly $0$;
- $\det(J_{TT}+\tfrac23I)=\tfrac8{27}-\tfrac23p_T+2q_T$ — residual exactly $0$;
- $\det(J_{TT}+\tfrac23I)=2\bigl(\tau_T-\tfrac2{27}\bigr)$ — the two goodness
  criteria of Lemma 2 and Lemma 3 are **literally proportional**, so the
  equivalence is an identity, not an implication chain;
- $\tau_T=\operatorname{tr}(w_iw_jw_k)$ — residual $0$ after Groebner reduction
  modulo the unit-norm ideal $\langle|u_k|^2-1\rangle$ (`domain=QQ`);
- Veronese quadric $w_i^2=\tfrac13w_i+\tfrac29I$ and $\|w_i\|^2=\tfrac23$ — both
  residual $0$ mod the same ideal.

**Layer P (PMC witness), exact in a number field.** All PMC entries lie in
$\mathbb Q(s,c)$, $s=\sin36^\circ/\sqrt2$, $c=\cos36^\circ/\sqrt2$. Naive
`sympy.simplify` on nested radicals like $\sqrt{(5-\sqrt5)/16}$ is unreliable —
in a first pass it produced **four spurious failures** on claims that are in fact
correct. Carrying every quantity as an element of an `AlgebraicField` with an
explicit minimal polynomial (so equality is decided by canonical form, with no
simplifier heuristics in the loop) removes the problem entirely.

> **Methodological warning for the campaign.** The spurious failures were an
> artifact of the *verification tool*, not of the mathematics. Any future
> re-verification of this corpus should use canonical number-field arithmetic,
> not `simplify`/`radsimp` heuristics, or it will generate false alarms against
> correct claims — and, worse, could mask real ones.

Confirmed exactly: $J_{\rm PMC}$ symmetric, $\operatorname{diag}=0$,
$J^2=I_6$, all six row masses $=1$, all 15 flow identities $=0$,
$P=(I+J)/2$ with $P^2=P$, $\operatorname{tr}P=3$, $P_{ii}=\tfrac12$,
$\sum_T\tau_T=\tfrac49$, $\sum_Tp_T=12$, $\sum_Tq_T=0$, and the three $0$-pairs
$\{0,5\},\{1,3\},\{2,4\}$ forming a perfect matching. **36/37 checks passed**;
the one non-pass was the discovery $c\notin\mathbb Q(s)$ (see §2), a fact about
the field, not a defect in the source.

---

## 2. Two findings in the PMC material

### 2.1 The active triples are not algebraically uniform [NEW, PROVED]

`slice-framework.md` §8 states that $G_{\rm PMC}$ is attained by 6 active
triples. True — but they are of **two distinct algebraic types**:

| type | count | triples | $p_T$ | $q_T$ |
|---|---|---|---|---|
| A | 4 | $(0,1,4),(0,2,3),(1,2,5),(3,4,5)$ | $\tfrac34$ | $\tfrac{K}{64}$, $K=(\sqrt2+\sqrt{10})\sqrt{5-\sqrt5}$ |
| B | 2 | $(0,1,5),(2,3,4)$ | $\tfrac{5-\sqrt5}{8}$ | $\mathbf 0$ |

Both types attain the *same* $\lambda_{\min}$, which is why the split is
invisible if one only records the value. The source's cubic is the type-A
charpoly; type B is a different polynomial with the same relevant root.

### 2.2 $G_{\rm PMC}$ has a much simpler closed form than the stated cubic [NEW, PROVED]

Type B is the transparent one: $q_T=0$ collapses the charpoly to
$y^3-p_Ty=y(y^2-p_T)$, so $\lambda_{\min}(J_{TT})=-\sqrt{p_T}=-\sin36^\circ$ and

$$\boxed{\;G_{\rm PMC}=\frac{1-\sin 36^\circ}{2}
=\frac{1-\sqrt{(5-\sqrt5)/8}}{2}=0.20610737385376343541564702268\ldots\;}$$

Verified exactly: this value is a root of the source's cubic
$256x^3-384x^2+144x-8-K$ (exact residual $0$), and it is the **smallest** of its
three real roots — confirming the source's ordering claim, which the first pass
had flagged. Its minimal polynomial over $\mathbb Q$ is
$256x^4-512x^3+304x^2-48x+1$ (degree **4**, not 3: the stated cubic has
coefficients in $\mathbb Q(K)$, not $\mathbb Q$).

Two consequences worth recording:

- The relevant number field is $\mathbb Q(s,c)$, **not** $\mathbb Q(s)$: we
  verified $c\notin\mathbb Q(s)$. The brief §7 anticipates "$\sqrt5$ and the
  PMC's defining cubic"; the sharper statement is that $G_{\rm PMC}$ itself is a
  degree-4 algebraic number expressible in $\sin36^\circ$.
- $G_{\rm PMC}>\tfrac16$ admits a **purely rational** certificate, with no
  radicals and no floating point:
  $$\tfrac{1-\sin36^\circ}{2}>\tfrac16\iff\sin36^\circ<\tfrac23
  \iff\tfrac{5-\sqrt5}{8}<\tfrac49\iff13<9\sqrt5\iff169<405.$$

---

## 3. Where GTZ(6,3) is actually tight [PROVED, re-derived]

Lemma 2.1 of `boundary-obstruction.md` was re-derived from the graph data
($G=\texttt{3xLR+path(1,1,1)}$, weights $w=(1,1,1,\tfrac95,\tfrac95,\tfrac95)$),
**21/21 checks exact**:

- $P_0$ is a genuine rank-3 orthoprojector; leverages
  $(\tfrac5{18},\tfrac5{18},\tfrac5{18},\tfrac{13}{18},\tfrac{13}{18},\tfrac{13}{18})$,
  all $>\tfrac16$ — so $A_0$ lies in the **core region** and Case A cannot touch it;
- of the 20 triples, exactly **10 are singular** ($\lambda_{\min}=0$ exactly) and
  the other **10 have $\lambda_{\min}=\tfrac16$ exactly**. Nothing lies strictly
  between $0$ and $\tfrac16$;
- hence $F(A_0)=\tfrac16$ and $f(A_0)=1/\sqrt6$ **exactly**;
- the pair $\{3,4\}$: $c^2=\tfrac{25}{324}$, $(\mu_1,\mu_2)=(1,\tfrac49)$,
  qualifies, and $\delta=0$ exactly — on the equality manifold; all four Schur
  slacks $q(0)=q(1)=q(2)=q(5)=0$, and the $q$-sum identity
  $\sum_pq(p)=-\delta$ holds exactly.

PSD-ness was decided via **all** principal minors of $P_{TT}-\tfrac16I$, not the
leading ones — the pitfall flagged in brief §4.

> **This is the fact that governs the whole project: GTZ(6,3) is TIGHT.**
> The bound $\tfrac16$ is attained, so there is no slack anywhere to absorb
> error. Two consequences:
>
> 1. **The slice is not where the problem is hard.** The slice minimum
>    $G_{\rm PMC}=0.2061$ sits $0.0394$ *above* $\tfrac16$, and the Nesterenko
>    leverages $\tfrac5{18},\tfrac{13}{18}$ are nowhere near the slice value
>    $\tfrac12$. Track B is a legitimate standalone lemma but it is **not** a
>    stepping stone to GTZ(6,3) — it certifies a region with a fat margin while
>    the binding configurations live off the slice, exactly as brief §5 cautions.
> 2. **"A finite number of algebraic checks" must mean exact checks at the
>    equality manifold.** Any scheme with an $\varepsilon$-margin, or any
>    floating-point global-optimization scheme without exact treatment of $A_0$,
>    cannot close the problem: the true minimum sits *on* the boundary.

---

## 4. The sharp-cone finding [NUMERICALLY SUPPORTED]

Reformulation R part (a) asserts each extremal is a *strict local minimum* of
$F$ (KKT residual $\le7\cdot10^{-11}$). We tested this far more aggressively and
found something stronger and more structured.

**(E1) Descent launched exactly at the extremal.** Starting at $A_0$ (where
$F=\tfrac16$ exactly) and running 3 annealed log-sum-exp stages plus 8 rounds of
Nelder–Mead polishing, $F$ stays at $\tfrac16$ to $2.5\cdot10^{-16}$. The
optimizer cannot leave the floor. This is the sharpest available test of GTZ(6,3)
near the binding configuration: it starts *on* the boundary rather than hoping to
find it.

**(E2) Perturbation shells.** Minimum of $F$ over 3000 random unit directions at
radius $r$, for $r=10^{-1}\ldots10^{-8}$:

| $r$ | $\min_{\|d\|=1}F(A_0+rd)-\tfrac16$ | $\mathrm{dev}/r$ | $\mathrm{dev}/r^2$ |
|---|---|---|---|
| $10^{-1}$ | $+6.019\cdot10^{-3}$ | $6.0\cdot10^{-2}$ | $6.0\cdot10^{-1}$ |
| $10^{-2}$ | $+8.271\cdot10^{-4}$ | $8.3\cdot10^{-2}$ | $8.3$ |
| $10^{-3}$ | $+6.454\cdot10^{-5}$ | $6.5\cdot10^{-2}$ | $6.5\cdot10^{1}$ |
| $10^{-4}$ | $+8.278\cdot10^{-6}$ | $8.3\cdot10^{-2}$ | $8.3\cdot10^{2}$ |
| $10^{-5}$ | $+7.571\cdot10^{-7}$ | $7.6\cdot10^{-2}$ | $7.6\cdot10^{3}$ |
| $10^{-6}$ | $+6.936\cdot10^{-8}$ | $6.9\cdot10^{-2}$ | $6.9\cdot10^{4}$ |
| $10^{-7}$ | $+6.148\cdot10^{-9}$ | $6.1\cdot10^{-2}$ | $6.1\cdot10^{5}$ |
| $10^{-8}$ | $+8.203\cdot10^{-10}$ | $8.2\cdot10^{-2}$ | $8.2\cdot10^{6}$ |

$\mathrm{dev}/r$ is **constant** ($\approx7\cdot10^{-2}$) over eight decades while
$\mathrm{dev}/r^2$ blows up. So the growth is **linear in $r$**, not quadratic:
$A_0$ is a *sharp* (conic, non-smooth) minimum, not a smooth quadratic one. This
is exactly what one expects of $F=\max_T\lambda_{\min}$ at a point with 10 active
triples — the max of finitely many smooth functions has a genuine kink.

**(E3) Directional-slope floor.** At $r=10^{-6}$ over $20{,}000$ random unit
directions, the slope $(F(A_0+rd)-\tfrac16)/r$ has
$$\min=5.02\cdot10^{-2},\quad 1\%=1.14\cdot10^{-1},\quad
\text{median}=2.72\cdot10^{-1},\quad\max=6.67\cdot10^{-1},$$
with **zero** directions negative and zero below $10^{-3}$. Explicitly
*minimizing* the slope over directions (60 Nelder–Mead starts) bottoms out at
$$\min_{\|d\|=1}\text{slope}\;\approx\;+6.5864\cdot10^{-3}\;>\;0,$$
and along that worst direction the deficit tracks $6.5864\cdot10^{-3}\cdot r$
cleanly across $r=10^{-2}\ldots10^{-9}$ (values $+6.71\cdot10^{-5}$ down to
$+6.59\cdot10^{-12}$).

> **Do not read $6.59\cdot10^{-3}$ as $\kappa$.** These probes are in ambient
> $A$-space, which contains the 3-dimensional $O(3)$ gauge along which $F$ is
> exactly constant; the true ambient minimum is $0$ and this number is only where
> the optimizer stalled. The gauge-free constant is $\kappa\approx7.3\cdot10^{-2}$
> (§4a). The *qualitative* content of E2/E3 — linear rather than quadratic growth,
> no negative direction — is correct and is what §4a then proves exactly.

> **Claim (sharp cone).** There is $\kappa>0$ with
> $F(A)\ge\tfrac16+\kappa\,\operatorname{dist}(A,\mathcal E)+o(\operatorname{dist})$
> near the extremal orbit $\mathcal E$.
>
> **Now PROVED exactly — see §4a.** The value quoted in an earlier draft of this
> section, $\kappa\approx6.59\cdot10^{-3}$, was **wrong**, and instructively so:
> it was minimized over ambient directions in the 18-dimensional $A$-space, which
> contains the 3-dimensional $O(3)$ gauge directions $A\mapsto A\Omega$ along which
> $P$, and hence $F$, is *exactly constant*. The true ambient minimum is therefore
> $0$, attained on a 3-dimensional subspace; $6.59\cdot10^{-3}$ was just where a
> Nelder–Mead run stalled on its way down. Working on the gauge-free Grassmannian
> gives the correct constant, **an order of magnitude larger** (§4a).

---

## 4a. The exact sharpness certificate [PROVED]

`verify/v6_sharpness.py`, **14/14 exact checks**. Everything below is exact
rational / $\mathbb Q(\sqrt5)$ linear algebra — no floating point is load-bearing.

**Setup.** Work on $Gr(3,6)=\{P:P^2=P,\ \operatorname{tr}P=3\}$, dimension 9, which
quotients out the $O(3)$ gauge automatically. At $P_0$:

- the active set is exactly the 10 non-singular triples, each with
  $\lambda_{\min}=\tfrac16$; the other 10 are singular at $\lambda_{\min}=0$ and
  cannot enter the max nearby;
- $\tfrac16$ is a **simple** eigenvalue of every active block (verified exactly),
  so each eigenspace $E_T=\operatorname{span}(v_T)$ is 1-dimensional and the
  one-sided derivative is a **max of 10 linear functionals**, not a min-eigenvalue
  of a compression:
  $$F'(P_0;\dot P)=\max_{T\ \text{active}}\ \langle G_T,\dot P\rangle,\qquad
  G_T=(v_Tv_T^{\mathsf T})\ \text{embedded on}\ T\times T;$$
- the tangent basis $\{u_an_b^{\mathsf T}+n_bu_a^{\mathsf T}\}$ is verified
  symmetric, traceless, $P_0BP_0=0$, $(I-P_0)B(I-P_0)=0$, and
  Frobenius-orthogonal with $\|B\|_F^2=2$.

**The certificate.** With $L\in\mathbb Q^{10\times9}$ the matrix of those
functionals in tangent coordinates:

$$\operatorname{rank}L=9,\qquad \dim\ker L^{\mathsf T}=1,\qquad
\lambda=\Bigl(\underbrace{\tfrac7{90},\ldots,\tfrac7{90}}_{9},\ \tfrac3{10}\Bigr),
\qquad \sum_T\lambda_TG_T=0 ,$$

with every multiplier **strictly positive** and $\sum_T\lambda_T=1$ exactly
(residuals of the KKT identity on all 9 tangent coordinates: exactly $0$). The
multiplier $\tfrac3{10}$ belongs to the triple $(3,4,5)$ — the three
heavy-leverage rows — and the nine $\tfrac7{90}$'s to the mixed triples.

**Why this proves sharpness.** Suppose $\max_T\langle G_T,\dot P\rangle\le0$ for
some $\dot P$. Pairing with the multipliers,
$0=\langle\sum_T\lambda_TG_T,\dot P\rangle=\sum_T\lambda_T\langle G_T,\dot P\rangle$,
a sum of non-positive terms with strictly positive weights, so every
$\langle G_T,\dot P\rangle=0$; since $\operatorname{rank}L=9$ this forces
$\dot P=0$. Hence the critical cone is $\{0\}$, $F'(P_0;\dot P)>0$ for every
$\dot P\neq0$, and by compactness of the unit sphere
$$\kappa=\min_{\|\dot P\|=1}F'(P_0;\dot P)>0 .$$

Strict positivity of *every* multiplier is exactly what upgrades "KKT point"
(condition (5.1) of `boundary-obstruction.md`, which only needs $\lambda\ge0$) to
"sharp minimum". Numerically $\kappa\approx7.334\cdot10^{-2}$ in the tangent
coordinates used, $\approx5.19\cdot10^{-2}$ under $\|\dot P\|_F=1$ — the *sign* is
proved exactly above; only the numeral is numerical.

> **Consequence.** Part (a) of Reformulation R is now **PROVED at $A_0$**, in a
> form strictly stronger than originally asserted: not merely a strict local
> minimum, but a sharp one with an explicit exact rational KKT certificate.
> Caveat on scope: this is proved **at this extremal**. The corpus reports nine
> Nesterenko extremals; running the same script at each is mechanical and cheap,
> and until that is done part (a) is PROVED only at $A_0$, not at all nine.

**(E4) Counterexample tripwire: never fired.** Across every experiment in this
session — $\ge2{,}600$ full-space descents, 60 slice descents, $24{,}000$ shell
perturbations, 200 kicked descents from $A_0$, and $20{,}000$ direction probes —
**no configuration with $F<\tfrac16-10^{-6}$ was found**. Consistent with the
campaign's $\ge119{,}045$ prior samples. Evidence, not proof.

> **CORRECTION (superseding an earlier draft of this paragraph).** An earlier
> version reported that generic descent "reaches only $F\approx0.1849$, *not* the
> true floor", and drew a Track C recommendation from it. That was **wrong** — an
> artifact of a `log-sum-exp` overflow (`exp(4000\lambda)\to\inf`, silently
> feeding `nan` to the optimizer). After the shift fix, `v5`-E3 over 500 starts
> gives $\min F=0.16666666680746603$ ($F-\tfrac16=+1.4\cdot10^{-10}$), with
> **357 / 500** starts landing within $10^{-4}$ of $\tfrac16$ and 456 within
> $10^{-2}$. So descent finds the floor easily and the extremals are **not** hard
> to reach. Independent algebraic confirmation from `v3` (1500 starts): the
> minimum value has minimal polynomial $36x^2-12x+1=(6x-1)^2$, i.e. exactly
> $\tfrac16$. The tripwire result is unaffected: 0 violations throughout.

The Track C consequence therefore *reverses*: seeding with the exact extremals is
good hygiene but not a necessity, because descent locates the floor unaided. What
Track C still needs is §6(ii) — a certified *radius* around each extremal.

---

## 4b. All six extremals, by a fully rational certificate [PROVED]

`verify/v9_rational_certificate.py` — **6 / 6 SHARP**, exact, no floating point.

### The surd obstruction and how it was removed

§4a's certificate at $A_0$ used unit eigenvectors and a Gram–Schmidt tangent
basis, which put the functional matrix $L$ in a field carrying *several* surds at
once — $\sqrt5,\sqrt7,\sqrt{10},\sqrt{14},\sqrt{35},\sqrt{70}$ for one extremal,
$\sqrt5,\sqrt{195}$ for another. An exact simplex needs rational input, and no
row-scaling clears a row mixing $\sqrt7$ and $\sqrt{10}$. The fix removes the
surds rather than fighting them.

Write $\tilde B=B_{\mathrm{red}}^{\mathsf T}$, $G=\tilde B^{\mathsf T}W\tilde B$
(rational), and

$$S:=\tilde B\,G^{-1}\tilde B^{\mathsf T}W .$$

Then $S$ is a **rational** idempotent with $\operatorname{tr}S=3$, and

$$P=W^{1/2}SW^{-1/2}.$$

Because $W^{1/2}$ is **diagonal**, this similarity restricts to every principal
block: $P_{TT}=W_T^{1/2}S_{TT}W_T^{-1/2}$, hence
$\operatorname{spec}(P_{TT})=\operatorname{spec}(S_{TT})$ and
$\lambda_{\min}(P_{TT})=\lambda_{\min}(S_{TT})$. The active set, the eigenvector
kernels, and the tangent space can therefore all be computed over $\mathbb Q$.

The functional stays rational too. With $v=D_Tu$ ($D_T=W_T^{1/2}$, $u$ the
rational kernel vector of $S_{TT}-\tfrac16I$) and $\dot P=D\dot SD^{-1}$,

$$v^{\mathsf T}\dot P_{TT}v=u^{\mathsf T}\bigl(W_T\dot S_{TT}\bigr)u,
\qquad v^{\mathsf T}v=u^{\mathsf T}W_Tu,$$

both rational — the $\sqrt w$ factors cancel pairwise. Verified: `L_rational =
True` for all six extremals.

### The certificate, and why one route suffices

Two independent exact tests (Gordan's alternative), both over $\mathbb Q$:

**(A) Primal.** Maximize $\sum_i(-Lz)_i$ over $\{Lz\le0,\ \|z\|_\infty\le1\}$ by
an exact rational simplex (Bland's rule). The feasible set contains $z=0$ with
objective $0$, and the objective is a sum of *non-negative* terms there. So an
exact maximum of $0$ forces $Lz=0$ for every feasible $z$, and
$\operatorname{rank}L=9$ then forces $z=0$. As $\{Lz\le0\}$ is a cone, the box
loses nothing. Hence $N=\{0\}$ and $\kappa>0$ by compactness.

**(B) Dual.** Exhibit $\lambda>0$ with $L^{\mathsf T}\lambda=0$, via a parametric
exact LP over $\ker L^{\mathsf T}$.

Route (A) is **self-sufficient**. Route (B) is a second, independent certificate
whose *success* yields publishable exact multipliers, but whose *failure* means
only that the parametric search found no witness — it does not weaken (A).

> **Validation of the LP logic on controls** (this mattered — an earlier version
> of route (B) silently returned `False` on a case §4a had already proved sharp,
> which is what exposed the bug): a sharp configuration returns exact value $0$;
> a non-sharp one returns a positive value *together with an explicit descent
> direction*; a rank-deficient one returns $0$ but is caught by the rank test.
> All three behave correctly, so a "SHARP" verdict here is not vacuous.

### Results

| verdict | graph | weights | leverages | active | $\dim\ker L^{\mathsf T}$ | route B |
|---|---|---|---|---|---|---|
| **SHARP** | `P(S(e,e,e),e,e,e)` | $1,1,1,\tfrac59,\tfrac59,\tfrac59$ | $\tfrac{13}{18}^{\times3},\tfrac5{18}^{\times3}$ | 10 | 1 | $(\tfrac7{90})^{\times9},\tfrac3{10}$ |
| **SHARP** | `P(S(e,e),S(e,e),e,e)` | $1,1,1,1,\tfrac58,\tfrac58$ | $\tfrac{11}{18}^{\times4},\tfrac5{18}^{\times2}$ | 12 | 3 | $(\tfrac3{20})^{\times4},(\tfrac1{20})^{\times8}$ |
| **SHARP** | `P(S(P(e,e,e),e),S(e,e))` | $1,1,1,\tfrac95,\tfrac95,\tfrac95$ | $\tfrac5{18}^{\times3},\tfrac{13}{18}^{\times3}$ | 10 | 1 | — |
| **SHARP** | `P(S(P(S(e,e),e,e),e),e)` | $1,1,\tfrac58,\tfrac58,1,1$ | $\tfrac{11}{18}^{\times4},\tfrac5{18}^{\times2}$ | 12 | 3 | — |
| **SHARP** | `P(S(P(e,e),P(e,e)),S(e,e))` | $1,1,1,1,\tfrac85,\tfrac85$ | $\tfrac7{18}^{\times4},\tfrac{13}{18}^{\times2}$ | 12 | 3 | — |
| **SHARP** | `P(S(P(e,e),e),S(P(e,e),e))` | $1,1,\tfrac85,1,1,\tfrac85$ | $\tfrac7{18}^{\times4},\tfrac{13}{18}^{\times2}$ | 12 | 3 | — |

All six: $\operatorname{rank}L=9$, exact LP value $=0$, $F=\tfrac16$ exactly, no
eigenvalue strictly inside $(0,\tfrac16)$, $\lambda_{\min}$ simple on every
active block. The first two reproduce §4a's multipliers exactly, cross-validating
`v6`/`v7` through an entirely different code path.

> **Part (a) of Reformulation R is now PROVED at every extremal of this census.**
> The remaining scope caveat is the census itself, not the certificate: see §5
> on the 6-vs-9 count and the `5/14, 9/14` configuration.

---

## 4c. The extremal census, and the construction the corpus omits

`boundary-obstruction.md` cites "the nine Nesterenko scaled-star matrices" but
names only four graphs and never states the weight rule. Two facts recovered here
and verified exactly against Lemma 2.1:

1. The construction is the weighted graphic-matroid / resistor-network one:
   $\tilde M=W^{1/2}B_{\mathrm{red}}^{\mathsf T}$,
   $P=\tilde M(\tilde M^{\mathsf T}\tilde M)^{-1}\tilde M^{\mathsf T}$.
2. **Leverage $=$ weight $\times$ effective resistance:** $\ell_e=w_e R_{\rm eff}(e)$.
   Check: $R_{\rm eff}(0,1)=\tfrac5{18}$ and $\tfrac95R_{\rm eff}(0,2)=\tfrac{13}{18}$,
   matching Lemma 2.1's leverages exactly.

An extremal is then a TTSP graph on 6 edges plus weights tuned so $F=\tfrac16$
exactly. `verify/v7_all_extremals.py` enumerates the 66 TTSP trees, solves for
those weights exactly, and finds **six**. Independently reproduced on MLCore (job
`gtz63-extremals-wagucx`, `16cpu-128ram`, region `ix-m5-sm12`, SUCCEEDED) with
identical verdicts on different hardware.

**Open discrepancy: six, not nine.** Either the one-symbol-per-edge-orbit weight
parametrization is too coarse to express some extremals, or some of the corpus's
nine are $O(3)$/relabelling duplicates. Unresolved; it does not affect any
certificate above, but it does bound the scope of the phrase "every extremal".

### 4c.1 A SEVENTH extremal, outside the scaled-star family [NEW, PROVED]

`verify/v11_seventh_exact.py` — **SHARP, exactly.**

`v3`'s unconstrained descent (1500 starts) landed on a configuration with
leverages $\tfrac5{14}$ ($\times3$) and $\tfrac9{14}$ ($\times3$) — **none** of the
four patterns the TTSP census produces ($\tfrac5{18},\tfrac{13}{18},\tfrac{11}{18},
\tfrac7{18}$) — with **13** active triples rather than 10 or 12. An exhaustive
search over all 66 TTSP trees with orbit-symmetric weights found **no** graph
giving the $\tfrac5{14},\tfrac9{14}$ pattern, so this is genuinely out of family,
not a relabelling of a census member.

It was re-polished to $|F-\tfrac16|\approx4\cdot10^{-15}$ (leverages matching
$\tfrac5{14}/\tfrac9{14}$ to $8\cdot10^{-15}$) and every entry identified by PSLQ
at 22 digits. Exactly five distinct off-diagonal magnitudes occur, all in
$\mathbb Q(\sqrt2,\sqrt5)$:

$$\tfrac5{14},\qquad \tfrac{\sqrt5}{21},\qquad \tfrac{5\sqrt5}{42},\qquad
\tfrac{5\sqrt2}{42},\qquad \tfrac{\sqrt{10}}{14}.$$

Rebuilt from those exact values, the matrix is **verified exactly** to satisfy
$P^{\mathsf T}=P$, $P^2=P$, $\operatorname{tr}P=3$, with leverages exactly
$\tfrac5{14},\tfrac9{14}$, 13 active triples at $\lambda_{\min}=\tfrac16$ exactly,
7 singular, none strictly inside $(0,\tfrac16)$, and $\lambda_{\min}$ **simple** on
every active block. The floats were used only to *guess* the entries; every
certified statement is a theorem about the exact algebraic object.

**A third rationalization trick was needed here.** One row of $L$ mixes
$1,\sqrt2,\sqrt5,\sqrt{10}$ simultaneously, so no *row* scaling clears it (this is
exactly what defeated `v8`). But rescaling the tangent basis **columns** works:
replacing $B_j$ by $c_jB_j$ with $c_j>0$ is the coordinate change $z_j\mapsto
z_j/c_j$, a positive diagonal map that leaves both $\{z:Lz\le0\}=\{0\}$ and
$\operatorname{rank}L$ invariant. Column factors
$(1,\sqrt5,\sqrt{10},\sqrt5,1,\sqrt2,\sqrt2,\sqrt{10},\sqrt5)$ rationalize $L$
completely, after which the exact rational simplex applies unchanged:
$\operatorname{rank}L=9$, exact LP optimum $=0$, hence $N=\{0\}$ and $\kappa>0$.

Here $\dim\ker L^{\mathsf T}=4$ and route (B) found no explicit positive
multiplier among the combinations tried — a reminder that route (A) is the
load-bearing one and route (B) is a bonus, not a requirement.

> **Consequence for scope.** "Every extremal is sharp" is now proved at **seven**
> configurations, one of which is provably outside the family the corpus
> describes. That is good news for the local half of Reformulation R, but it also
> shows the extremal set is **larger and less well-catalogued** than
> `boundary-obstruction.md` suggests. Since part (b) (the global gap) is a
> statement about *all* KKT points, an uncatalogued extremal family is a genuine
> obstacle to any proof strategy that proceeds by exhausting a known list. This
> is the sharpest new structural information from this session.

---

---

## 5. Numeric sanity / reproduction block

All runs deterministic. From `~/Documents/gtz` with the project venv:

```bash
.venv/bin/python -u verify/v1_foundations.py    # 36/37; the 1 non-pass is c not in Q(s)
.venv/bin/python -u verify/v4_nesterenko.py     # 21/21 exact
.venv/bin/python -u verify/v5_tightness.py 500  # E1/E1b/E2/E3, master seed 20260801
.venv/bin/python -u verify/v2_search.py 2000 4000   # master seed 20260730
.venv/bin/python -u verify/v3_extremal.py 1500       # master seed 20260731
```

Environment: Python 3.12.3, numpy 2.5.1, sympy 1.14.0, scipy 1.18.0,
python-flint present (Arb available for future interval work); 22 cores.

- s0/v1: generic identities — all residuals exactly $0$.
- s1/v4: Nesterenko extremal — 21/21 exact, $F(A_0)=\tfrac16$.
- s2/v5-E1: descent from $A_0$ — $F-\tfrac16=+2.50\cdot10^{-16}$, 0 violations.
- s3/E2: 8 shells $\times$ 3000 directions — 0 shells below $\tfrac16-10^{-6}$.
- s4/E3: 20000 directions, 60 slope minimizations — min ambient slope
  $+6.59\cdot10^{-3}$ (gauge-contaminated; see the warning in §4).
- s6/v6: exact sharpness certificate — 14/14, $\operatorname{rank}L=9$,
  $\dim\ker L^{\mathsf T}=1$, multipliers $(\tfrac7{90})^{\times9},\tfrac3{10}$ all
  $>0$, KKT residuals exactly $0$, $\kappa\approx7.334\cdot10^{-2}$.
- s5: $(1.1)$ excess formula residual exactly $0$; $q$-sum identity
  $\max|\sum_pq(p)+\delta|=2.87\cdot10^{-12}$ over 400 random qualifying samples.

**Known limitation.** The `v2`/`v3` large-scale runs were still executing when
this document was written; the counts above reflect completed runs plus the
standalone shell/slope probes. The `v2` slice run at 60 starts gave **59/60
landing on $G_{\rm PMC}$ to $10^{-16}$** with none below — sharper than the
campaign's 64-start evidence for $C_{\rm slice}$, but a 2000-start confirmation
is still pending and $C_{\rm slice}$ remains CONJECTURED.

---

## 6. Recommended next step

§4a settled the step this section previously recommended. What remains, in
increasing order of difficulty:

**(i) Cheap and mechanical — the other eight extremals.** `verify/v6_sharpness.py`
is written against one projector; the corpus reports nine Nesterenko extremals.
Re-running it per extremal (different graph/weights, same code path) would upgrade
part (a) of Reformulation R from "PROVED at $A_0$" to "PROVED at all nine". Expect
the same structure — simple $\lambda_{\min}$, rank-9 $L$, one-dimensional
multiplier space — but the exact multipliers will differ per extremal, and a
failure of *simplicity* at some extremal is the one thing that would genuinely
complicate matters (then $F'$ is a min-eigenvalue of a compression, and the check
becomes a small SDP rather than a nullspace computation).

**(ii) Quantify the neighbourhood, not just the derivative.** §4a certifies the
*first-order* rate $\kappa>0$. Converting that into a rigorous statement of the
form "$F\ge\tfrac16$ on the explicit ball $\|\dot P\|\le\rho_0$" needs a
second-order remainder bound — a Lipschitz constant for $\dot P\mapsto F'$ on
$Gr(3,6)$, which is available in closed form because the blocks are $3\times3$ and
the eigenvalues are algebraic. With $\rho_0$ in hand, Track C's branch-and-bound
can **excise** each extremal neighbourhood analytically instead of trying to
resolve it by bisection — which is exactly where interval methods stall on a tight
problem, and the main practical obstacle to Track C today.

**(iii) The global gap — part (b), still OPEN.** No KKT point of $F$ has value
$<\tfrac16$. Untouched by this session. §4a does not weaken it; it isolates it.

## 7. Honest summary

Nothing here closes GTZ(6,3). What changed:

1. The foundational corpus **survived independent exact re-verification** —
   including a stronger, generic-identity form of the slice lemmas. One
   correction and one simplification were found in the PMC material (§2), both
   improvements rather than errors, and one methodological warning about radical
   simplifiers was recorded.
2. **GTZ(6,3) is tight** ($F=\tfrac16$ exactly at $A_0$), re-derived exactly. This
   reframes the "finite algebraic checks" hope: such checks must be exact at the
   equality manifold, and the slice — the easiest-looking target — is provably
   not the binding region.
3. The extremal is a **sharp** minimum, and this is now **PROVED exactly** (§4a)
   rather than merely observed: an exact rational KKT certificate with all ten
   multipliers strictly positive, $(\tfrac7{90})^{\times9}$ and $\tfrac3{10}$,
   giving $\kappa\approx7.3\cdot10^{-2}$. This is the first genuinely new *proved*
   result of this line of attack, and it upgrades part (a) of Reformulation R at
   $A_0$ from NUMERICALLY SUPPORTED to PROVED.
4. En route, a numerical claim in this document's own §4 was found to be wrong
   (the ambient $\kappa\approx6.6\cdot10^{-3}$, contaminated by $O(3)$ gauge
   directions) and corrected by the exact computation. Recorded rather than
   quietly patched, per brief §8 rule 1.
5. **All six census extremals plus a seventh, out-of-family one are now PROVED
   sharp** (§4b, §4c.1), via a certificate made fully rational by three successive
   tricks: the $S=\tilde BG^{-1}\tilde B^{\mathsf T}W$ conjugation, per-row
   positive scaling, and per-column positive scaling of the tangent basis.
6. The seventh extremal is the sharpest *structural* news: the extremal set is
   **not** exhausted by the nine scaled-star matrices the corpus lists. Any proof
   of part (b) that works by enumerating a known extremal list is therefore on
   shakier ground than it appeared.
7. The counterexample tripwire never fired.

The remaining gap is unchanged in location and now sharper in description: local
behaviour at $A_0$ is not merely well understood but **certified**, by finite exact
means, and the same script will settle the other eight extremals mechanically. The
**global** statement — that no other critical point of $F$ dips below $\tfrac16$ —
is where GTZ(6,3) still lives, and this session did not move it. What §4a buys is
that the local half of Reformulation R is no longer conjectural, so the *entire*
remaining difficulty is now isolated in part (b).
