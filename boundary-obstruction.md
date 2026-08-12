# The single-pair boundary obstruction and the minimax/KKT reformulation for GTZ(6,3)

**Deliverable:** standalone, referee-grade document for target **VAL-THM-BND-001**,
promoted from attack round 4 (`attacks/round-4-boundary.md`;
`numerics/study-attack-r4-boundary/`).  Exact-arithmetic certificate and float
sanity: `numerics/study-boundary-obstruction/` (deterministic, seed 20260722; §6).

**Status summary.**

| Statement | Status |
|---|---|
| Theorem BND (single-pair boundary impossibility, §3) | **PROVED** (exact arithmetic; §3.2–§3.4) |
| Corollary NG (no single-pair small-excess extension lemma; the h-branch is the maximal reach of the Pair→Extension program, §4) | **PROVED** |
| Reformulation R, part (a): every extremal is a strict local min of $F=\max_T\lambda_{\min}(P_{TT})=f^2$ (§5.2) | **NUMERICALLY SUPPORTED** (KKT residual $\le 7\cdot10^{-11}$ at all nine exact extremals) — NOT proved |
| Reformulation R, part (b): no KKT point of $F$ has value $<1/6$ (the global SOS gap, §5.3) | **OPEN** |
| Logical equivalence GTZ(6,3) $\iff$ (a) $\wedge$ (b) (§5.1) | **PROVED** (equivalence only; (a),(b) not both established) |
| SOS structural consequence: a certificate must vanish on the equality manifold (§5.4) | **PROVED** (consequence of §3 first-order tightness) |

No hidden gaps: every claim carries one of the statuses above; nothing marked
PROVED depends on a NUMERICALLY SUPPORTED or OPEN item.

---

## 0. Scope — read this first (overclaim guard)

> **This document does NOT prove GTZ(6,3).**
>
> It proves a **no-go**: the single-pair second-moment program (`proofs/pair-lemma.md`
> → `proofs/extension-lemma.md`, the *Pair Lemma → Extension Lemma* route)
> **cannot** close GTZ(6,3).  Concretely (Theorem BND, §3): for **every**
> threshold $\delta_0>0$ there is a configuration $A\in\mathrm{St}(6,3)$ carrying
> a qualifying pair (Gram eigenvalue $\mu_2>1/6$) of arbitrarily small excess
> $\delta\in(0,\delta_0]$ **all of whose triples are bad**, while GTZ still holds
> at $A$ through a triple *not* containing the pair.  Hence no lemma of the shape
> "small-excess qualifying pair $\Rightarrow$ $1/6$-good triple through it" can
> exist for any threshold (Corollary NG, §4).  The Extension Lemma's h-branch
> (`extension-lemma.md`, Lemma E) is the **maximal reach** of that program.
>
> The correct replacement (Reformulation R, §5) recasts GTZ(6,3) near the
> equality manifold as a **minimax/KKT** statement for
> $F(A)=\max_T\lambda_{\min}(P_{TT})=f(A)^2$ (GTZ(6,3) $\iff\min F\ge\frac16$):
> *(a)* every extremal is a strict local minimum of $F$
> — **NUMERICALLY SUPPORTED**, KKT stationarity residual $\le 7\cdot10^{-11}$ at
> all nine exact extremals, not proved — **and** *(b)* no KKT point of $F$ has
> value $<\frac16$ — the **global gap**, a Positivstellensatz / real-algebraic-
> geometry problem that remains **OPEN**.
>
> No configuration with $f<1/\sqrt6-10^{-6}$ was found anywhere (GTZ tripwire
> over $\ge 119{,}045$ configurations never fired; `study-attack-r4-boundary`).
> This is evidence, not proof; part (a) of the reformulation is **numerical**,
> the global gap (b) is **open**.

**Dependencies.**  The construction is built at a validated Nesterenko extremal
and uses the extension machinery of **VAL-THM-EXT-001**
(`proofs/extension-lemma.md`: (M0)–(M3), the Schur criterion, the q-sum
identity, the equality case §3.1) and the extremal structure of
**VAL-THM-SLICE-001** (`proofs/slice-framework.md`) / the validated catalogs
**VAL-NUM-002/003**.  Everything used from those is PROVED there; §3 below adds
only exact linear-algebra at one explicit configuration and is self-contained
modulo that machinery.

---

## 1. Setting and notation

Shared project notation (MEMORY.md; `extension-lemma.md` §1–§2).  $A\in
\mathrm{St}(6,3)$ (so $A^TA=I_3$, $\sum_i\ell_i=3$); $P=AA^T$ is the orthogonal
projector onto $\mathrm{col}(A)\subset\mathbb R^6$ ($P^T=P$, $P^2=P$,
$\operatorname{tr}P=3$).  For $I\subseteq[6]$, $\mathrm{Gram}(I)=P_{II}$; for a
triple $|I|=3$, $\lambda_{\min}(P_{II})=\sigma_{\min}(A_I)^2$, and $I$ is
**$1/6$-good** iff $\lambda_{\min}(P_{II})\ge\frac16$.  The GTZ objective is
$$f(A)=\max_{|I|=3}\sigma_{\min}(A_I)=\sqrt{\max_{|I|=3}\lambda_{\min}(P_{II})},
\qquad \text{GTZ}(6,3)\iff \min_{\mathrm{St}(6,3)} f\ge\tfrac1{\sqrt6}.$$

Fix a pair $\{i,j\}$ with Gram $\Gamma=P_{\{i,j\}\{i,j\}}=\begin{pmatrix}
\ell_i&c_{ij}\\ c_{ij}&\ell_j\end{pmatrix}$, eigenvalues $\mu_1\ge\mu_2$,
and $h(\mu)=\frac{1-\mu}{6\mu-1}$.  The pair **qualifies** iff $\mu_2>\frac16$;
its **excess** is
$$\delta \;=\; h(\mu_1)+h(\mu_2)-\tfrac13 .$$
For $p\notin\{i,j\}$ the **Schur slack** (`extension-lemma.md` (M2)) is
$$q(p)=\ell_p-\tfrac16-v_p^T\bigl(\Gamma-\tfrac16 I_2\bigr)^{-1}v_p,\qquad
v_p=(c_{ip},c_{jp})^T,$$
and by (M2) the triple $\{i,j,p\}$ is $1/6$-good $\iff q(p)\ge0$; the pair
**extends** iff $\max_p q(p)\ge0$.  By the q-sum identity (M3),
$\sum_{p\notin\{i,j\}}q(p)=-\delta$.

**A rational form of the excess.**  Writing $t=\operatorname{tr}\Gamma=
\ell_i+\ell_j$ and $d=\det\Gamma=\ell_i\ell_j-c_{ij}^2$, the symmetric function
$h(\mu_1)+h(\mu_2)$ is rational in $(t,d)$:
$$h(\mu_1)+h(\mu_2)=\frac{(1-\mu_1)(6\mu_2-1)+(1-\mu_2)(6\mu_1-1)}
{(6\mu_1-1)(6\mu_2-1)}=\frac{7t-2-12d}{36d-6t+1},
\qquad
\boxed{\ \delta=\frac{7t-2-12d}{36d-6t+1}-\frac13\ }\tag{1.1}$$
(using $\mu_1+\mu_2=t$, $\mu_1\mu_2=d$).  In particular $\delta$ is a **rational
function of the Gram entries**, with no eigenvalue radicals; likewise $q(p)$ is
rational in the Gram entries wherever $\Gamma-\frac16 I_2$ is invertible.  This
is what makes the exact-arithmetic certificate of §3 possible.

---

## 2. The Nesterenko extremal and the equality pair (exact)

Let $A_0\in\mathrm{St}(6,3)$ be the Nesterenko scaled-star matrix
(arXiv:2511.02387; `numerics/study-6-3/nesterenko.py`, validated in VAL-NUM-002)
for the two-terminal series–parallel graph
$$G=\texttt{3xLR+path(1,1,1)}=(\mathrm P,\,e,\,e,\,e,\,(\mathrm S,\,e,\,e,\,e)),$$
i.e. three parallel edges $e_0,e_1,e_2$ between the terminals $0,1$ and a series
path $e_3=(0,2),\,e_4=(2,3),\,e_5=(3,1)$.  Its scaled-star projector is
$$P_0=\tilde M\,(\tilde M^T\tilde M)^{-1}\tilde M^T,\qquad
\tilde M=W^{1/2}B_{\mathrm{red}}^T\in\mathbb R^{6\times3},$$
with graph-induced edge weights $w=(1,1,1,\tfrac95,\tfrac95,\tfrac95)$ and
reduced incidence (vertices $1,2,3$)
$$B_{\mathrm{red}}=\begin{pmatrix}-1&-1&-1&0&0&-1\\ 0&0&0&-1&1&0\\
0&0&0&0&-1&1\end{pmatrix}.$$
Because $w_3=w_4=w_5=\frac95$ and $w_0=w_1=w_2=1$, all entries of $P_0$ lie in
$\mathbb Q(\sqrt5)$ (off-diagonal entries carry a factor $\sqrt{w_ew_f}\in
\{1,\tfrac3{\sqrt5},\tfrac95\}$).

> **Lemma 2.1 (extremal data, exact).**  $P_0$ is a rank-$3$ orthogonal
> projector.  Its leverages are
> $$\ell_0=\ell_1=\ell_2=\tfrac5{18},\qquad \ell_3=\ell_4=\ell_5=\tfrac{13}{18},$$
> every non-singular triple has $\lambda_{\min}(P_{0,II})=\frac16$ and the
> remaining triples are singular, so $f(A_0)=1/\sqrt6$ **exactly**; there are
> exactly $10$ active ($1/6$-good, $\lambda_{\min}=\frac16$) triples.  The pair
> $\{3,4\}$ has
> $$\Gamma_0=\begin{pmatrix}\tfrac{13}{18}&c\\ c&\tfrac{13}{18}\end{pmatrix},
> \quad c^2=\tfrac{25}{324},\qquad (\mu_1,\mu_2)=\Bigl(1,\tfrac49\Bigr),$$
> so it **qualifies** ($\mu_2=\frac49>\frac16$) and lies **exactly on the
> equality manifold**: $\delta=0$ (this is the E-I spectrum $(1,\frac{n+2}{3n})$
> of `extension-lemma.md` §8.2).  Moreover all four Schur slacks vanish,
> $q(0)=q(1)=q(2)=q(5)=0$, so all four triples through $\{3,4\}$ are active with
> $\lambda_{\min}=\frac16$.
>
> **Status: PROVED** (exact; sympy certificate `exact_obstruction.py`, part (A)/(B)).

*Proof.*  $\tilde M$ has full column rank $3$ (the reduced incidence of a
connected graph on $4$ vertices has rank $3$), so $P_0$ is the orthogonal
projector onto $\mathrm{col}(\tilde M)=$ the scaled-star space; $P_0^T=P_0$,
$P_0^2=P_0$, $\operatorname{tr}P_0=3$ are verified exactly.  The leverages and
$\Gamma_0$ are read off $P_0$; $(\mu_1,\mu_2)$ are the eigenvalues of $\Gamma_0$,
here $\frac{13}{18}\pm|c|=\frac{13}{18}\pm\frac5{18}=1,\frac49$.  Then by (1.1)
with $t=\frac{13}9$, $d=\mu_1\mu_2=\frac49$, the numerator is
$7\cdot\frac{13}9-2-12\cdot\frac49=\frac{91}9-\frac{18}9-\frac{48}9=\frac{25}9$
and the denominator is $36\cdot\frac49-6\cdot\frac{13}9+1=16-\frac{78}9+1=
\frac{75}9=\frac{25}3$, so
$$h(\mu_1)+h(\mu_2)=\frac{25/9}{25/3}=\frac13$$
(equivalently $h(1)=0$, $h(\frac49)=\frac{5/9}{6\cdot4/9-1}=\frac{5/9}{5/3}
=\frac13$), so $\delta=h(\mu_1)+h(\mu_2)-\frac13=0$.  The four $q(p)=0$ follow either from the equality
case of Lemma E (`extension-lemma.md` §3.1: $\delta=0$ with no strict extension
forces every $q(p)=0$) or by direct exact evaluation.  The two-valued spectrum
$\{\frac16,0\}$ and the count of $10$ active triples are a finite exact check
over all $\binom63=20$ triples.  $\square$

The ten active triples are the four through $\{3,4\}$ —
$\{3,4,0\},\{3,4,1\},\{3,4,2\},\{3,4,5\}$ — and six others, all containing row
$5$: $\{0,3,5\},\{0,4,5\},\{1,3,5\},\{1,4,5\},\{2,3,5\},\{2,4,5\}$.

---

## 3. Theorem BND: the single-pair boundary lemma is impossible

> **Theorem BND (single-pair boundary impossibility).**
> For every $\delta_0>0$ there exists $A\in\mathrm{St}(6,3)$ carrying a
> qualifying pair $\{i,j\}$ (Gram eigenvalue $\mu_2>\frac16$) of excess
> $\delta\in(0,\delta_0]$ such that **every** triple $\{i,j,p\}$,
> $p\notin\{i,j\}$, is bad ($q(p)<0$), while GTZ(6,3) holds at $A$ through a
> triple **not** containing $\{i,j\}$ (indeed $f(A)>1/\sqrt6$).
>
> **Status: PROVED** (exact arithmetic; the witnessing family and all four
> descent signs, the excess-increase sign, and the external-triple sign are
> certified in $\mathbb Q(\sqrt5)$, §3.2–§3.4).

The witnesses form an **explicit one-parameter family** obtained by moving off
$A_0$ (Lemma 2.1) along a Stiefel tangent common-descent direction at the
equality pair $\{3,4\}$.  This extends the icosahedral obstruction
(`extension-lemma.md` §7, where every pair has h-sum $\frac{13}{11}$, i.e. excess
$\delta=\frac{13}{11}-\frac13=\frac{28}{33}$) all the way down to excess $\to0$.

### 3.1 The construction

Because $f$, $\delta$ and each $q(p)$ depend on $A$ only through $P=AA^T$ (they
are invariant under $A\mapsto AQ$, $Q\in O(3)$), we move on the Grassmannian.
The tangent space at $P_0$ is
$$T_{P_0}=\{\dot P=\dot P^T:\ \dot P=P_0\dot P+\dot P P_0\}
=\{P_0 S(I-P_0)+(I-P_0)S P_0:\ S=S^T\},\tag{3.1}$$
of dimension $k(n-k)=9$ (differentiate $P^2=P$: $\dot P=\dot PP+P\dot P$, i.e.
$P_0\dot P P_0=0$ and $(I-P_0)\dot P(I-P_0)=0$).  Fix the **integer symmetric
generator**
$$S=\begin{pmatrix}
 2& 2& 2&-1& 0&-2\\ 2& 2&-1& 1& 2& 0\\ 2&-1&-2& 2& 2& 1\\
-1& 1& 2&-4&-2& 1\\ 0& 2& 2&-2&-2& 0\\ -2& 0& 1& 1& 0& 2\end{pmatrix},
\qquad \dot P:=P_0S(I-P_0)+(I-P_0)SP_0\in T_{P_0}.\tag{3.2}$$
Let $P(\varepsilon)$ be any smooth curve of rank-$3$ projectors with $P(0)=P_0$,
$P'(0)=\dot P$ (e.g. the Grassmann geodesic, or $P(\varepsilon)=A(\varepsilon)
A(\varepsilon)^T$ with $A(\varepsilon)=\mathrm{qr\_retract}(A_0,\varepsilon\dot P
A_0)$, the family built in §6).  Write $\delta(\varepsilon)$,
$q_p(\varepsilon)$, $\mu_2(\varepsilon)$ for the excess, Schur slacks and
smaller pair-eigenvalue of $\{3,4\}$ along the curve, and
$\lambda^*(\varepsilon)=\lambda_{\min}(P(\varepsilon)_{T^*T^*})$ for the external
triple $T^*=\{0,3,5\}$.  All are smooth near $\varepsilon=0$ (see §3.4).

### 3.2 The five descent signs (exact)

> **Lemma 3.1 (exact first-order data).**  With $A_0$, $\{3,4\}$, $\dot P$ as
> above, all evaluated at $\varepsilon=0$:
> $$\delta'(0)=\frac{59}{9}+\frac{4\sqrt5}{45}>0,\qquad
> \begin{aligned}
> q_0'(0)&=\tfrac29-\tfrac{34\sqrt5}{45}<0, & q_1'(0)&=-\tfrac79-\tfrac{16\sqrt5}{45}<0,\\
> q_2'(0)&=-\tfrac{19}9-\tfrac{4\sqrt5}{45}<0, & q_5'(0)&=-\tfrac{35}9+\tfrac{10\sqrt5}{9}<0.
> \end{aligned}$$
> **Status: PROVED** (exact; `exact_obstruction.py`, part (C)).

*Proof.*  $\dot P$ is a genuine tangent: $\dot P^T=\dot P$ and $\dot P=P_0\dot P
+\dot PP_0$ hold identically in $\mathbb Q(\sqrt5)$ (both verified symbolically).
Along $P(\varepsilon)=P_0+\varepsilon\dot P+O(\varepsilon^2)$ the Gram trace and
determinant of $\{3,4\}$ are polynomials in $\varepsilon$, so by (1.1) the excess
$\delta(\varepsilon)$ is a rational function of $\varepsilon$ whose denominator
$36d-6t+1$ equals $\frac{25}3\ne0$ at $\varepsilon=0$; differentiating (1.1) at
$\varepsilon=0$ gives the stated $\delta'(0)$.  Likewise each
$q_p(\varepsilon)=P_{pp}-\frac16-v_p^T(\Gamma-\frac16I_2)^{-1}v_p$ is rational in
$\varepsilon$ (the $2\times2$ inverse is $\frac1{\det}\mathrm{adj}$, and
$\det(\Gamma-\frac16I_2)=(\mu_1-\frac16)(\mu_2-\frac16)=\frac{25}{108}\ne0$ at
$\varepsilon=0$ since $\mu_2=\frac49>\frac16$); differentiating gives the four
$q_p'(0)$.  Each sign
is an **exact comparison in $\mathbb Q(\sqrt5)$**: e.g. $q_0'(0)<0\iff
\frac{10}{45}<\frac{34\sqrt5}{45}\iff 10<34\sqrt5\iff 100<5780$, and
$q_5'(0)<0\iff10\sqrt5<35\iff 500<1225$; $\delta'(0)>0$ and
$q_1'(0),q_2'(0)<0$ are sign-definite termwise.  sympy decides all five
inequalities exactly.  $\square$

Thus $\dot P$ is a **common strict-descent direction** for the four Schur slacks
through $\{3,4\}$ *and* strictly **increases the excess**.

### 3.3 GTZ is preserved by an external triple (exact)

> **Lemma 3.2 (external triple rises).**  The external triple $T^*=\{0,3,5\}$
> (which does not contain both $3$ and $4$) has $\lambda_{\min}(P_{0,T^*T^*})=
> \frac16$ a **simple** eigenvalue, with exact spectrum
> $\{\frac16,\frac59,1\}$, and
> $$\lambda^{*\prime}(0)=\psi^T\dot P_{T^*T^*}\psi=\frac{115}{63}
> -\frac{13\sqrt5}{63}>0,$$
> where $\psi$ is the unit $\frac16$-eigenvector of $P_{0,T^*T^*}$.
> **Status: PROVED** (exact; `exact_obstruction.py`, part (D)).

*Proof.*  The three eigenvalues $\frac16,\frac59,1$ are distinct, so $\frac16$
is simple and $\lambda^*(\varepsilon)$ is analytic near $\varepsilon=0$ with
$\lambda^{*\prime}(0)=\psi^T\dot P_{T^*T^*}\psi$ (first-order perturbation of a
simple eigenvalue; $\psi$ computed exactly from
$(P_{0,T^*T^*}-\frac16I_3)\psi=0$).  The value $\frac{115}{63}-\frac{13\sqrt5}
{63}>0\iff13\sqrt5<115\iff845<13225$, decided exactly.  $\square$

### 3.4 Proof of Theorem BND

All of $\delta(\varepsilon)$, $q_p(\varepsilon)$ ($p\in\{0,1,2,5\}$),
$\mu_2(\varepsilon)$ and $\lambda^*(\varepsilon)$ are smooth on a neighborhood of
$\varepsilon=0$: the first two are rational with nonvanishing denominators at
$0$ (§3.2), $\mu_2$ is a simple eigenvalue of $\Gamma(\varepsilon)$ at $0$
($\mu_1=1\ne\frac49=\mu_2$) hence smooth, and $\lambda^*$ is smooth by Lemma 3.2.
Their values and exact first derivatives at $0$ are
$$\delta(0)=0,\ \delta'(0)>0;\quad q_p(0)=0,\ q_p'(0)<0;\quad
\mu_2(0)=\tfrac49>\tfrac16;\quad \lambda^*(0)=\tfrac16,\ \lambda^{*\prime}(0)>0.$$
By first-order Taylor expansion there is $\varepsilon_1>0$ such that for all
$\varepsilon\in(0,\varepsilon_1]$:

* $\mu_2(\varepsilon)>\frac16$ — the pair $\{3,4\}$ **qualifies**;
* $\delta(\varepsilon)>0$ and $\delta$ is strictly increasing from $\delta(0)=0$
  — the pair is **strictly diffuse**, with $\delta(\varepsilon)\to0^+$;
* $q_p(\varepsilon)<0$ for every $p\in\{0,1,2,5\}$ — **all four triples through
  $\{3,4\}$ are bad**; by (M2) none of them is $1/6$-good;
* $\lambda^*(\varepsilon)>\frac16$, hence
  $f(A(\varepsilon))\ge\sqrt{\lambda^*(\varepsilon)}>\frac1{\sqrt6}$ — **GTZ
  holds at $A(\varepsilon)$ through $T^*=\{0,3,5\}$**, a triple not containing
  the pair.

Finally, given any $\delta_0>0$, since $\delta$ is continuous, strictly
increasing on $(0,\varepsilon_1]$ with $\delta(0^+)=0$, choose
$\varepsilon\in(0,\varepsilon_1]$ with $\delta(\varepsilon)\le\delta_0$ (possible
as $\delta(\varepsilon)\to0$); then $A=A(\varepsilon)$ is the required witness
with excess $\delta\in(0,\delta_0]$.  $\blacksquare$

**Remark 3.3 (both scalings linear; the deficit is spread evenly).**  Near
$\varepsilon=0$, $\delta(\varepsilon)=\delta'(0)\varepsilon+O(\varepsilon^2)$ and
$f(A(\varepsilon))-\frac1{\sqrt6}$ scale **linearly** in $\varepsilon$; by the
q-sum identity $\sum_p q_p(\varepsilon)=-\delta(\varepsilon)$, the four slacks
average $-\delta/4$, i.e. the deficit is dumped essentially evenly onto the four
triples through the pair.  The boundary band is therefore **first-order tight**
(cf. §5.4).

**Remark 3.4 (four of nine equality pairs are unprotected).**  The same LP for a
common-descent tangent is feasible at an equality pair for **4 of the 9**
built-in extremal graphs (`3xLR+path(1,1,1)` $\{3,4\}$, `LR+path(3,1,1)`
$\{0,4\}$, `LR+path(2,2,1)` $\{0,5\}$, `2path(2,1)+2path(2,1)` $\{2,5\}$;
`study-attack-r4-boundary/s4`); the exact certificate here treats the first.
One unprotected pair suffices for Theorem BND.  At the remaining equality pairs
the four gradients $Dq(p)$ already have $0$ in their convex hull (Gordan), so no
common descent exists — but this does **not** rescue a single-pair lemma, since
the *unprotected* pairs already produce the obstruction.

---

## 4. Corollary NG: no single-pair small-excess extension lemma

> **Corollary NG.**  There is **no** threshold $\delta_0>0$ and no lemma of the
> form
> $$\text{"pair }\{i,j\}\text{ qualifies with excess }\delta\le\delta_0
> \ \Longrightarrow\ \text{some triple }\{i,j,p\}\text{ is }1/6\text{-good"}$$
> valid on all of $\mathrm{St}(6,3)$.  Equivalently, the Extension Lemma's
> h-branch $h(\mu_1)+h(\mu_2)\le\frac13$ (`extension-lemma.md`, Lemma E) is the
> **maximal reach** of the single-pair second-moment program: the boundary
> $\delta=0$ cannot be crossed by any argument that looks only at one pair's
> second moments and asks for a good triple *through that pair*.
>
> **Status: PROVED** (immediate from Theorem BND).

*Proof.*  A lemma as stated would assert, for its $\delta_0$, that every
qualifying pair of excess $\le\delta_0$ extends ($\max_p q(p)\ge0$).  Theorem BND
produces a configuration with a qualifying pair of excess $\delta\in(0,\delta_0]$
and $q(p)<0$ for **all** $p$, i.e. $\max_p q(p)<0$ — the pair does not extend.
Contradiction.  Since $\delta_0>0$ was arbitrary, no such threshold exists.  The
"maximal reach" reading is the contrapositive: Lemma E extends every pair with
$\delta\le0$, and Theorem BND shows the closed half $\delta\le0$ cannot be
enlarged to any $\delta\le\delta_0$ with $\delta_0>0$.  $\blacksquare$

**Consequence for the mission.**  Combined with the icosahedral obstruction
(`extension-lemma.md` §7: the diffuse ETF where every pair has h-sum
$\frac{13}{11}$ — excess $\frac{28}{33}$ — and GTZ holds only through sign
coherence), Corollary NG closes
the entire **Pair Lemma → Extension Lemma** route (rounds 1–3) as a path to
GTZ(6,3).  The single-pair second-moment invariants $(\mu_1,\mu_2)$ of a pair
carry, provably, *insufficient* information to certify a good triple in the
diffuse regime — a genuinely multi-triple / global mechanism is required (§5).
The compensation is not even cleanly two-pair: at a failed small-excess pair the
certifying good triple shares exactly one row with it, yet **disjoint** good
triples also occur ($618/2684$ witnesses; `study-attack-r4-boundary/s3`), so no
"second-pair" localization holds either.

---

## 5. Reformulation R: the minimax / KKT global gap

The correct replacement for the single-pair program is a **minimax** view of the
same objective.  This section states the equivalence (PROVED), the local-minimum
half (**NUMERICALLY SUPPORTED**, not proved) and the global half (**OPEN**), and
records the structural consequence for SOS certificates.

### 5.1 The equivalence

Work with the minimax objective in $\lambda$-units,
$$F(A)=\max_{|I|=3}\lambda_{\min}(P_{II})=f(A)^2,$$
a maximum of $\binom63=20$ functions, each smooth where its $\lambda_{\min}$ is a
simple eigenvalue.  GTZ(6,3) is exactly $\min_{\mathrm{St}(6,3)}F\ge\frac16$
(equivalently $\min f\ge\frac1{\sqrt6}$).  A **KKT point** of $F$ (equivalently a
Clarke-stationary point of the locally-Lipschitz $F$) with value $\frac16$ we
call an **extremal**; the nine Nesterenko scaled-star matrices are extremals.

> **Reformulation R (logical equivalence).**  GTZ(6,3) holds **iff**
> *(a)* every extremal is a local minimum of $F$, **and**
> *(b)* no KKT point of $F$ has value $<\frac16$.
>
> **Status: PROVED** (as a logical equivalence).

*Proof.*  ($\Rightarrow$)  If $\min F\ge\frac16$ then every point with
value $\frac16$ is a global — hence local — minimum, giving (a); and no
point, KKT or otherwise, has value $<\frac16$, giving (b).
($\Leftarrow$)  $F$ is continuous on the compact manifold $\mathrm{St}(6,3)$, so
it attains a global minimum $m$ at some $A_\star$, which is a KKT point.  By (b),
$m\ge\frac16$.  Thus $\min F\ge\frac16$, i.e. GTZ(6,3).  (Part (a)
is not needed for this direction; it is the *local* certificate that the
extremal values $\frac16$ are not saddles feeding a lower branch, and is
the numerically-verified content of §5.2.)  $\blacksquare$

The value of the reformulation is that it **isolates the crux**: (a) is a
finite, local, checkable condition at the nine extremals; (b) is a single global
polynomial-optimization statement.  Theorem BND (§3) is precisely the statement
that (b) cannot be reached by single-pair second moments — a single pair sees
only $4$ of the $10$–$13$ active gradients at an extremal (§5.2), and those $4$
need not positively span.

### 5.2 Part (a): extremals are strict local minima — NUMERICALLY SUPPORTED

The first-order necessary condition for $A_0$ ($F$-value $\frac16$, active set
$\mathcal A$ = the $1/6$-good triples) to be a local min of $F$ is: **for every
tangent $X$, $\max_{T\in\mathcal A} D\lambda_{\min}(P_{TT})[X]\ge0$**, i.e. no
common descent lowering all active triples — equivalently
$$0\in\operatorname{conv}\{\nabla\lambda_{\min}(P_{TT}):T\in\mathcal A\}\subset
T^*_{A_0}\mathrm{St}(6,3).\tag{5.1}$$

> **Numerical finding (`study-attack-r4-boundary/s5`).**  At **all nine** exact
> Nesterenko extremals, (5.1) holds: the minimum-norm convex combination of the
> active-triple gradients has residual $\le 7\cdot10^{-11}$, and the active
> gradients **positively span** the tangent dual (the LP "all active gradients
> $\le-1$" is infeasible).  Hence each extremal is a **strict** local minimum
> with a boundary band of positive radius.
>
> **Status: NUMERICALLY SUPPORTED — NOT PROVED.**  The residual $\le7\cdot
> 10^{-11}$ is consistent with exact $0$-in-hull but is a floating-point KKT
> computation; a proof would require exact positive-spanning of the active
> gradients at each extremal (finite but not done here).

This reconfirms the round-3-variational KKT finding independently.  The
single-pair-vs-set contrast is explicit (`s5`): at every extremal the $2$–$4$
active triples through some equality pair do **not** positively span (a
single-pair descent exists — this is exactly the direction $\dot P$ of §3), yet
the **full** active set does (0-in-hull restored).  Theorem BND is the exact,
proved shadow of this: §3 exhibits the single-pair descent in $\mathbb Q(\sqrt5)$
and Lemma 3.2 shows the set-level compensation ($\lambda^*$ rising).

### 5.3 Part (b): no KKT point below $1/6$ — the global gap, OPEN

> **Part (b) [OPEN].**  No KKT point of $F=\max_T\lambda_{\min}(P_{TT})$ on
> $\mathrm{St}(6,3)$ has value $<\frac16$ (equivalently, no KKT point of $f$ has
> value $<\frac1{\sqrt6}$).

This is a real-algebraic-geometry / Positivstellensatz statement over the
compact variety $\mathrm{St}(6,3)$ — the true remaining crux for GTZ(6,3).  It is
exactly what the single-pair program cannot see (§3–§4) and what the SOS /
symmetry-reduced-Lasserre track targets.  Concrete finite attack (round-4 §next,
not executed here): enumerate the KKT/active-set **combinatorial types** (which
subsets of the $20$ triples can be simultaneously active at a critical point —
heavily constrained by the spanning-tree / matroid structure of
`slice-framework.md` and §2 here) and rule out sub-$1/6$ critical values
type-by-type.  **No configuration with $f<1/\sqrt6-10^{-6}$ was found** across
$\ge119{,}045$ configurations (Haar, perturbed extremals, slice, adversarial
descent; tripwire never fired) — evidence for (b), not a proof.

### 5.4 SOS consequence: a certificate must vanish on the equality manifold

> **Proposition 5.1 (flat certificate necessity).**  Any Positivstellensatz
> certificate of $F=f^2\ge\frac16$ over $\mathrm{St}(6,3)$ must have its
> nonnegativity multipliers **vanish on the equality manifold** (the extremal
> variety): a strictly-positive Putinar certificate cannot exist.
>
> **Status: PROVED** (consequence of the first-order tightness of §3 / Remark 3.3).

*Proof sketch.*  By Theorem BND / Remark 3.3 the GTZ margin $f-\frac1{\sqrt6}$ is
**first-order tight** at the extremal variety: along the family $A(\varepsilon)$
it vanishes linearly, $f(A(\varepsilon))-\frac1{\sqrt6}=\Theta(\varepsilon)$, and
there are directions (the four bad triples) along which individual triple margins
go strictly negative.  A certificate expressing $f^2-\frac16=\sigma_0+\sum_i
\sigma_i g_i$ (SOS $\sigma$'s, constraint polynomials $g_i$) with all
$\sigma\succeq0$ strictly positive at an extremal would force a strictly positive
lower bound on $f^2-\frac16$ there, contradicting $f^2-\frac16=0$ on the
(positive-dimensional) extremal variety.  Hence the multipliers must vanish on
that variety — a degenerate/flat certificate.  $\square$

This explains why generic Lasserre levels stall at the averaging identity
(slice level-2 $=-\frac7{135}$ exactly, reproduced in
`study-attack-r3-diffuse-63/s5_sos`; `study-attack-r4-boundary` §4) and predicts
that a working certificate needs the equality-manifold ideal built in
(symmetry-reduced level $\ge4$ with the extremal variety as equality
constraints).  That reduction ($S_6\ltimes$switching, order $23040$) is heavy
infrastructure and is **not** carried out here.

---

## 6. Numeric sanity block

All checks are deterministic (**seed 20260722**), run from the workspace root
with `.venv/bin/python`.  Directory `numerics/study-boundary-obstruction/`.

**Exact certificate — `exact_obstruction.py` (sympy, field $\mathbb Q(\sqrt5)$;
`ALL EXACT CHECKS PASSED`).**  Every sign below is decided by exact algebraic
comparison, no floating point:

| Claim (this doc) | Exact result |
|---|---|
| Lemma 2.1: $P_0$ projector, $\operatorname{tr}=3$; leverages $\frac5{18},\frac{13}{18}$ | verified in $\mathbb Q(\sqrt5)$ |
| Lemma 2.1: every non-singular triple $\lambda_{\min}=\frac16$; $f(A_0)=1/\sqrt6$; $10$ active | exact (all $20$ triples) |
| Lemma 2.1: pair $\{3,4\}$ $(\mu_1,\mu_2)=(1,\frac49)$, $\mu_2>\frac16$, $\delta=0$, four $q(p)=0$ | exact |
| Lemma 3.1: $\dot P$ symmetric, $\dot P=P_0\dot P+\dot PP_0$ (valid tangent) | exact identity |
| Lemma 3.1: $\delta'(0)=\frac{59}9+\frac{4\sqrt5}{45}>0$ | exact sign |
| Lemma 3.1: $q_0'(0)=\frac29-\frac{34\sqrt5}{45}<0$, $q_1',q_2',q_5'<0$ | exact signs |
| Lemma 3.2: $P_{0,T^*T^*}$ spectrum $\{\frac16,\frac59,1\}$, $\frac16$ simple; $\lambda^{*\prime}(0)=\frac{115}{63}-\frac{13\sqrt5}{63}>0$ | exact |

**Float sanity — `float_family.py` (`SANITY PASSED`).**  Builds the actual
Stiefel family $A(\varepsilon)=\mathrm{qr\_retract}(A_0,\varepsilon\dot PA_0)$:

| $\varepsilon$ | excess $\delta$ | $\max_p q(p)$ | all 4 bad | $f-1/\sqrt6$ |
|---|---|---|---|---|
| $10^{-2}$ | $+7.74\cdot10^{-2}$ | $-1.71\cdot10^{-2}$ | yes | $+1.64\cdot10^{-2}$ |
| $10^{-3}$ | $+6.84\cdot10^{-3}$ | $-1.43\cdot10^{-3}$ | yes | $+1.67\cdot10^{-3}$ |
| $3\cdot10^{-4}$ | $+2.03\cdot10^{-3}$ | $-4.24\cdot10^{-4}$ | yes | $+5.01\cdot10^{-4}$ |
| $10^{-4}$ | $+6.76\cdot10^{-4}$ | $-1.41\cdot10^{-4}$ | yes | $+1.67\cdot10^{-4}$ |

Observed slopes match the exact derivatives of Lemma 3.1/3.2 to $\le7\cdot
10^{-10}$ ($\delta'(0)=6.75431715$, $q_0'(0)=-1.46725136$, …,
$\lambda^{*\prime}(0)=+1.36398597$).  Excess and GTZ margin both scale linearly,
all four triples through $\{3,4\}$ stay bad, $f>1/\sqrt6$ throughout — the
constructed family realizes Theorem BND.

**Broader evidence — `numerics/study-attack-r4-boundary/` (seed 20260722).**
`s1`/`s3`: $\ge119{,}045$ configs, GTZ tripwire ($f<1/\sqrt6-10^{-6}$) never
fired; non-extending diffuse pairs with excess down to $2.4\cdot10^{-4}$.  `s2`:
pure Haar produces **zero** such pairs in $2\cdot10^5$ tries (the phenomenon is
near-extremal only — the boundary band).  `s4`: the LP common-descent family for
$4/9$ graphs (float precursor of §3).  `s5`: KKT residual $\le7\cdot10^{-11}$ at
all nine extremals, active gradients positively span (§5.2).

---

## 7. Provenance and citations

- Promoted from `attacks/round-4-boundary.md` (attack round 4) with the §3
  obstruction family **upgraded to exact arithmetic** here (the round-4
  recommendation "target 3: exact-arithmetic upgrade"); the float precursor is
  `study-attack-r4-boundary/s4`.
- Uses the extension machinery of **`proofs/extension-lemma.md`** (VAL-THM-EXT-001):
  the Schur criterion (M2), the q-sum identity (M3), the equality case §3.1, and
  the h-branch Lemma E; and the extremal structure of
  **`proofs/slice-framework.md`** (VAL-THM-SLICE-001) with the validated catalogs
  **VAL-NUM-002/003**.  The icosahedral obstruction extended to excess $\to0$ is
  `extension-lemma.md` §7.
- The Nesterenko scaled-star construction: arXiv:2511.02387, via
  `numerics/study-6-3/nesterenko.py` (validated VAL-NUM-002).
- Standard facts used by name: differentiation of $P^2=P$ for the Grassmann
  tangent (3.1); first-order perturbation of a **simple** eigenvalue
  ($\lambda'=\psi^T\dot A\psi$) (Lemma 3.2); Gordan's theorem / $0$-in-convex-hull
  for common-descent feasibility (Remarks 3.4, §5.2); Sylvester inertia via the
  Schur criterion (inherited from `extension-lemma.md` (M2)); a continuous
  function on a compact set attains its minimum (§5.1).

---

## 8. Status recap

- **Theorem BND** (single-pair boundary impossibility) — **PROVED** (exact
  arithmetic, $\mathbb Q(\sqrt5)$).
- **Corollary NG** (no single-pair small-excess lemma; h-branch is maximal
  reach) — **PROVED**.
- **Reformulation R** equivalence — **PROVED**; part (a) local-min —
  **NUMERICALLY SUPPORTED** (KKT residual $\le7\cdot10^{-11}$, not proved);
  part (b) no KKT point below $1/6$ — **OPEN** (the global SOS gap).
- **Proposition 5.1** (flat-certificate necessity) — **PROVED**.
- **GTZ(6,3): NOT proved.**  No sub-$1/\sqrt6$ configuration found
  ($\ge119{,}045$ configs).
