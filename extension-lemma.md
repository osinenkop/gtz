# The Extension Lemma and coverage branches for the pair → triple bridge ($k = 3$)

**Deliverable:** standalone, referee-grade proof of the Extension Lemma package,
promoted from attack round 2 (`attacks/round-2-pair-to-triple.md`), target
VAL-THM-EXT-001.
Numerics: `numerics/study-extension-lemma/` (deterministic, seed 20260718; §9).

**Status summary.**

| Statement | Status |
|---|---|
| (M0)–(M4) moment/Schur/convexity machinery (§2) | PROVED |
| Lemma E (Extension Lemma, $k=3$, unconditional) (§3) | PROVED |
| Lemma E equality case (§3.1) | PROVED |
| General-$k$ q-sum identity, Lemma E$_k$, $k=2$ infeasibility (§3.2) | PROVED |
| Corollary E′ (strict interior: $f > 1/\sqrt n$ inside the h-branch) (§8.1) | PROVED |
| Corollary T (trace threshold $(4n+2)/(3n)$, unconditional) (§4) | PROVED |
| Corollary S (spectral special cases) (§4.1) | PROVED |
| Proposition U (near-unit row, $\varepsilon_n = \frac{n-2}{n^2(n-1)}$) (§5) | PROVED CONDITIONAL ON arXiv:2604.05944, $k=2$ theorem (verified in `proofs/sp-verification.md`, VAL-SP-001) |
| Proposition D (dual branch at $(6,3)$, $\ell_i+\ell_j \le 5/9$) (§6) | PROVED CONDITIONAL ON `proofs/duality.md` Theorem 1 / Corollary 1 (VAL-DUAL-001) |
| Icosahedral obstruction: every pair has h-sum $= 13/11 > 1/3$ (§7) | PROVED (exact) |
| At the ico, coherent triples are $1/6$-good, each pair in exactly 2 (§7.2) | PROVED (exact + finite verification) |
| Equality families E-I, E-II lie on the curve $h(\mu_1)+h(\mu_2)=2/n$ (§8.2) | PROVED (symbolic) |
| C_eq (extremals lie on the equality manifold) (§8.3) | CONJECTURED (numerically supported; NOT proved) |
| The bridge C_ext / the diffuse regime (§8.4) | OPEN |

No hidden gaps: every claim below carries one of the statuses above; nothing
marked PROVED depends on a CONJECTURED or OPEN item.

---

## 0. Scope — read this first (overclaim guard)

> **This document does NOT prove GTZ$(n,3)$.**
>
> It proves an *unconditional sufficient condition*: whenever some pair of
> rows is "spectrally heavy" — its $2\times2$ Gram eigenvalues
> $\mu_1 \ge \mu_2 > 1/n$ satisfy $h(\mu_1) + h(\mu_2) \le 2/n$ with
> $h(\mu) = \frac{1-\mu}{n\mu-1}$ (the **h-branch**) — a $1/n$-good triple
> exists (Lemma E), together with three coverage corollaries (T, U, D) and an
> exact obstruction (§7) showing the h-branch is **empty** at the icosahedral
> ETF$(6,3)$.  Consequently the pair → triple bridge is *reduced to*, not
> closed on, the **diffuse regime** (every pair strictly above the h-curve),
> which **remains OPEN** (§8.4).  The icosahedral obstruction proves that no
> argument confined to a single pair's second moments can close the diffuse
> regime — a genuinely different mechanism is required there.
>
> The empirical law that all validated $k=3$ extremal configurations lie
> exactly on the boundary manifold $h(\mu_1)+h(\mu_2) = 2/n$ is
> **CONJECTURED** (C_eq, §8.3): it is supported by all validated extremal
> data but is NOT proved here, and nothing in this document depends on it.

---

## 1. Setting and notation

Shared project notation (MEMORY.md).  For integers $n \ge k \ge 1$,
$$\mathrm{St}(n,k) = \{A \in \mathbb R^{n\times k}:\ A^TA = I_k\}.$$
Rows are $r_1,\dots,r_n \in \mathbb R^k$; leverages $\ell_i = \|r_i\|^2$;
inner products $c_{ij} = \langle r_i, r_j\rangle$; $P = AA^T$.  The Stiefel
constraint is equivalently $\sum_{i=1}^n r_ir_i^T = I_k$; taking traces,
$\sum_i \ell_i = k$.  $P$ is the orthogonal projector onto the column span
of $A$: $P^T = P$ and $P^2 = A(A^TA)A^T = P$, so $0 \preceq P \preceq I_n$.

For $I \subseteq [n]$, $\mathrm{Gram}(I) = P_{II}$ (principal submatrix; row
order inside $I$ is immaterial — a permutation acts on $P_{II}$ by orthogonal
conjugation, preserving eigenvalues).  For $|I| = k$,
$\mathrm{Gram}(I) = A_IA_I^T$, so its eigenvalues are the squared singular
values of $A_I$ and
$$\lambda_{\min}(\mathrm{Gram}(I)) = \sigma_{\min}(A_I)^2 .$$
At $k = 3$, GTZ$(n,3)$ at $A$ asks for a triple $I$ with
$\lambda_{\min}(\mathrm{Gram}(I)) \ge 1/n$; such a triple is called
**$1/n$-good** (equivalently $\sigma_{\min}(A_I) \ge 1/\sqrt n$).
Default below: $k = 3$ (so $\sum_i \ell_i = 3$); §3.2 records the general-$k$
variants explicitly.

Throughout, for $n \ge 2$ define on $\mu \in (\tfrac1n, \infty)$
$$g(\mu) = \frac{(n-1)\mu}{n\mu - 1}, \qquad
  h(\mu) = g(\mu) - 1 = \frac{1-\mu}{n\mu - 1}.$$

---

## 2. Extension machinery

Fix $A \in \mathrm{St}(n,3)$ and a pair $\{i,j\}$, $i \ne j$.  Let $J \in
\mathbb R^{2\times 3}$ have rows $r_i, r_j$, and set
$$\Gamma = JJ^T = \begin{pmatrix}\ell_i & c_{ij}\\ c_{ij} & \ell_j
\end{pmatrix} = P_{\{i,j\}\{i,j\}},$$
with eigenvalues $\mu_1 \ge \mu_2 \ge 0$, and for every row index $p \in [n]$
$$v_p = Jr_p = (c_{ip},\, c_{jp})^T \in \mathbb R^2 .$$

**(M0) $0 \preceq \Gamma \preceq I_2$; in particular $\mu_1 \le 1$.**

*Proof.* $\Gamma$ is a principal $2\times2$ submatrix of the orthoprojector
$P$: for a unit $x \in \mathbb R^2$, extend it by zeros to $\tilde x \in
\mathbb R^n$ supported on $\{i,j\}$; then $x^T\Gamma x = \tilde x^TP\tilde x
\in [0,1]$ because $0 \preceq P \preceq I_n$. $\square$

**Status: PROVED.**

**(M1) Moment identities:**
$$\sum_{p=1}^n v_pv_p^T = \Gamma, \qquad
  v_iv_i^T + v_jv_j^T = \Gamma^2, \qquad
  \text{hence}\quad \sum_{p \notin \{i,j\}} v_pv_p^T = \Gamma - \Gamma^2 .$$

*Proof.* First, $\sum_p v_pv_p^T = J\bigl(\sum_p r_pr_p^T\bigr)J^T =
JI_3J^T = \Gamma$.  Second, $v_i = Jr_i = (\langle r_i,r_i\rangle,
\langle r_j,r_i\rangle)^T = (\ell_i, c_{ij})^T = \Gamma e_1$ and likewise
$v_j = \Gamma e_2$, so $v_iv_i^T + v_jv_j^T = \Gamma(e_1e_1^T +
e_2e_2^T)\Gamma = \Gamma I_2 \Gamma = \Gamma^2$ (using $\Gamma^T = \Gamma$).
Subtract. $\square$

**Status: PROVED.**  (Float check F1: $5.0\cdot10^{-16}$.)

**(M2) Schur criterion.**  *Assume $\mu_2 > 1/n$, so $\Gamma - \frac1nI_2
\succ 0$.  For $p \notin \{i,j\}$ define the **Schur slack**
$$q(p) = \ell_p - \frac1n - v_p^T\Bigl(\Gamma - \frac1n I_2\Bigr)^{-1}v_p .$$
Then*
$$\lambda_{\min}\bigl(\mathrm{Gram}\{i,j,p\}\bigr) \ge \frac1n
\iff q(p) \ge 0,$$
*and if $q(p) > 0$ then $\lambda_{\min}(\mathrm{Gram}\{i,j,p\}) > 1/n$
strictly; if $q(p) = 0$ then $\lambda_{\min}(\mathrm{Gram}\{i,j,p\}) = 1/n$
exactly.*

*Proof.* In row order $(i,j,p)$,
$$M := \mathrm{Gram}\{i,j,p\} - \frac1n I_3 =
\begin{pmatrix} \Gamma - \frac1n I_2 & v_p\\[2pt] v_p^T & \ell_p - \frac1n
\end{pmatrix}.$$
Let $L = \begin{pmatrix} I_2 & 0\\ -v_p^T(\Gamma-\frac1nI_2)^{-1} & 1
\end{pmatrix}$ (invertible).  A direct block computation gives
$LML^T = \mathrm{diag}\bigl(\Gamma - \tfrac1nI_2,\ q(p)\bigr)$.
By Sylvester's law of inertia $M$ and $LML^T$ have the same numbers of
positive/zero/negative eigenvalues.  Since $\Gamma - \frac1nI_2 \succ 0$:
$M \succeq 0 \iff q(p) \ge 0$, i.e. $\lambda_{\min}(\mathrm{Gram}) \ge 1/n
\iff q(p) \ge 0$; $M \succ 0 \iff q(p) > 0$ (all three inertia-eigenvalues
positive), i.e. $\lambda_{\min} > 1/n$; and $q(p) = 0$ makes $M \succeq 0$
singular, i.e. $\lambda_{\min}(\mathrm{Gram}) = 1/n$ exactly. $\square$

**Status: PROVED.**

**(M3) q-sum identity ($k = 3$).**  *Assume $\mu_2 > 1/n$.  Then*
$$\sum_{p \notin \{i,j\}} q(p)
 = \frac{2n+2}{n} - g(\mu_1) - g(\mu_2),$$
*and consequently*
$$\sum_{p \notin \{i,j\}} q(p) \;\ge\; 0
\quad\Longleftrightarrow\quad h(\mu_1) + h(\mu_2) \;\le\; \frac2n ,$$
*with equality on the left iff equality on the right.*

*Proof.*  Summing the definition of $q$ over $p \notin \{i,j\}$ and using
$v_p^TXv_p = \operatorname{Tr}(Xv_pv_p^T)$ with linearity of the trace and
(M1):
$$\sum_{p\notin\{i,j\}} q(p)
 = \Bigl(\sum_{p\notin\{i,j\}}\ell_p\Bigr) - \frac{n-2}{n}
 - \operatorname{Tr}\Bigl[\bigl(\Gamma - \tfrac1nI_2\bigr)^{-1}
   (\Gamma - \Gamma^2)\Bigr].$$
Here $\sum_{p\notin\{i,j\}}\ell_p = 3 - \ell_i - \ell_j = 3 -
\operatorname{Tr}\Gamma = 3 - \mu_1 - \mu_2$ (using $\sum_p \ell_p = 3$).
Both $\Gamma - \Gamma^2$ and $(\Gamma - \frac1nI_2)^{-1}$ are functions of
$\Gamma$, hence share its eigenbasis: writing $\Gamma = \sum_{a=1,2}
\mu_a u_au_a^T$ with orthonormal $u_a$,
$$\operatorname{Tr}\Bigl[\bigl(\Gamma - \tfrac1nI\bigr)^{-1}(\Gamma -
\Gamma^2)\Bigr] = \sum_{a=1,2}\frac{\mu_a(1-\mu_a)}{\mu_a - 1/n}.$$
Therefore, with $3 - \frac{n-2}n = \frac{2n+2}{n}$,
$$\sum_{p\notin\{i,j\}} q(p) = \frac{2n+2}{n}
 - \sum_{a=1,2}\Bigl[\mu_a + \frac{\mu_a(1-\mu_a)}{\mu_a - 1/n}\Bigr].$$
The bracket simplifies:
$$\mu + \frac{\mu(1-\mu)}{\mu - 1/n}
 = \frac{\mu(\mu - \tfrac1n) + \mu(1-\mu)}{\mu - \tfrac1n}
 = \frac{\mu\,(1 - \tfrac1n)}{\mu - \tfrac1n}
 = \frac{(n-1)\mu}{n\mu - 1} = g(\mu).$$
This proves the identity.  For the equivalence: $g = 1 + h$, so
$\frac{2n+2}n - g(\mu_1) - g(\mu_2) = \frac2n - h(\mu_1) - h(\mu_2)$, and
the displayed sum is $\ge 0$ (resp. $= 0$) iff $h(\mu_1)+h(\mu_2) \le \frac2n$
(resp. $=$). $\square$

**Status: PROVED.**  (Float check F1: $1.8\cdot10^{-11}$,
conditioning-limited near $\mu_2 \approx 1/n$; symbolic check S1.)

**(M4) Properties of $h$.**  *On the domain $(\frac1n, \infty)$, $n \ge 2$:*

1. *$h$ is strictly decreasing and strictly convex:*
   $h'(\mu) = \dfrac{1-n}{(n\mu-1)^2} < 0$, $\;h''(\mu) =
   \dfrac{2n(n-1)}{(n\mu-1)^3} > 0$.
2. *$h(1) = 0$ and $h \ge 0$ on $(\frac1n, 1]$ (numerator $1 - \mu \ge 0$,
   denominator $n\mu - 1 > 0$); $h < 0$ on $(1,\infty)$.*
3. *For $\mu \in (\frac1n,\infty)$:\; $h(\mu) \le \frac2n \iff \mu \ge
   \frac{n+2}{3n}$.  Special values:*
   $$h\Bigl(\frac{n+2}{3n}\Bigr) = \frac2n, \qquad
     h\Bigl(\frac{2n+1}{3n}\Bigr) = \frac1{2n}, \qquad
     h\Bigl(\frac{n+1}{2n}\Bigr) = \frac1n .$$

*Proof.*  1., 2.: direct differentiation and inspection of signs
($(n\mu-1)>0$ on the domain).  3.: since $n\mu - 1 > 0$,
$h(\mu) \le \frac2n \iff n(1-\mu) \le 2(n\mu - 1) \iff n + 2 \le 3n\mu$.
The three values follow by substitution, e.g.
$h\bigl(\tfrac{2n+1}{3n}\bigr) = \frac{(n-1)/(3n)}{(2n-2)/3} = \frac1{2n}$.
$\square$

**Status: PROVED.**  (Symbolic checks S1–S2.)

---

## 3. The Extension Lemma

> **Lemma E (pair → triple extension, $k = 3$, unconditional).**
> *Let $n \ge 3$, $A \in \mathrm{St}(n,3)$, and let $\{i,j\}$ be a pair whose
> Gram eigenvalues satisfy $\mu_2 > \frac1n$ and*
> $$h(\mu_1) + h(\mu_2) \;\le\; \frac2n, \qquad
>   h(\mu) = \frac{1-\mu}{n\mu-1}.$$
> *Then there exists $p \notin \{i,j\}$ with
> $\lambda_{\min}\bigl(\mathrm{Gram}\{i,j,p\}\bigr) \ge \frac1n$; i.e. the
> triple $\{i,j,p\}$ witnesses GTZ$(n,3)$ at $A$.*
>
> **Status: PROVED** (unconditional; no assumption on the other rows, no
> Pair-Lemma "goodness", no core/leverage assumption anywhere).

*Proof.*  By (M3) the hypothesis gives $\sum_{p\notin\{i,j\}} q(p) \ge 0$.
The index set $[n]\setminus\{i,j\}$ has $n - 2 \ge 1$ elements, so the sum is
over a nonempty set; if every term were $< 0$ the sum would be $< 0$.  Hence
some $p \notin \{i,j\}$ has $q(p) \ge 0$, and (M2) converts this into
$\lambda_{\min}(\mathrm{Gram}\{i,j,p\}) \ge 1/n$. $\blacksquare$

Numeric verification of the conclusion (§9): 2906 planted boundary-stressed
instances + 294 Haar h-branch instances, 0 violations.

### 3.1 Equality case

**Claim.**  *If $h(\mu_1) + h(\mu_2) = \frac2n$ exactly (and $\mu_2 > 1/n$)
and no third row extends strictly (i.e. $q(p) \le 0$ for all
$p \notin \{i,j\}$), then $q(p) = 0$ for **every** $p \notin \{i,j\}$: all
$n-2$ triples through $\{i,j\}$ satisfy
$\lambda_{\min}(\mathrm{Gram}\{i,j,p\}) = \frac1n$ **exactly**.*

*Proof.*  By (M3), equality in the hypothesis gives
$\sum_{p\notin\{i,j\}} q(p) = 0$; combined with $q(p) \le 0$ for all $p$,
every $q(p) = 0$.  The last part of (M2) then gives
$\lambda_{\min} = 1/n$ exactly for each triple. $\square$

**Status: PROVED.**  (This is precisely the active-triple structure observed
at the sharp $(6,3)$ catalog class 0 — see §8.3; that *observation* is
catalog data, not used anywhere in this document's proofs.)

### 3.2 General $k$, and infeasibility at $k = 2$

The machinery of §2 uses $k = 3$ only twice: $\sum_p \ell_p = k$ in (M3),
and the fact that $\mathrm{Gram}\{i,j,p\}$ is $3\times3$ (true for every
$k$).  Repeating the computation verbatim for $A \in \mathrm{St}(n,k)$,
$k \ge 2$, $\mu_2 > 1/n$:
$$\sum_{p \notin \{i,j\}} q(p) = k - \frac{n-2}{n} - g(\mu_1) - g(\mu_2),$$
whence:

> **Lemma E$_k$.**  *Let $n \ge 3$, $A \in \mathrm{St}(n,k)$, $k \ge 2$, and
> let $\{i,j\}$ be a pair with $\mu_2 > \frac1n$ and*
> $$h(\mu_1) + h(\mu_2) \;\le\; k - 3 + \frac2n .$$
> *Then some triple $\{i,j,p\}$ has
> $\lambda_{\min}(\mathrm{Gram}\{i,j,p\}) \ge \frac1n$.*
>
> **Status: PROVED** (same proof as Lemma E; note for $k \ge 4$ a triple is
> only a *fragment* of a $k$-subset, so this does **not** produce a GTZ
> witness for $k \ge 4$).

*Consistency check at $k = 2$:* the condition reads $h(\mu_1)+h(\mu_2) \le
-1 + \frac2n < 0$ for $n \ge 3$, while (M0) gives $\mu_1 \le 1$, so
$h(\mu_1), h(\mu_2) \ge 0$ by (M4.2) — the hypothesis is **infeasible**.
This is forced: for $k = 2$, any three rows lie in a rank-$\le2$ frame, so
$\mathrm{Gram}\{i,j,p\} = A_{\{i,j,p\}}A_{\{i,j,p\}}^T$ has rank $\le 2$ and
$\lambda_{\min} = 0 < 1/n$ — a correct theory must make the $k=2$ condition
empty, and it does.
**Status: PROVED.**  (Float checks F2: general-$k$ identity to
$4.2\cdot10^{-13}$; on 300 Haar $\mathrm{St}(n,2)$ samples every qualifying
pair had h-sum $\ge 0$ and $\sum q < 0$, max $\sum q = -0.559$.)

---

## 4. The trace-only corollary

> **Corollary T (trace threshold, unconditional).**
> *Let $n \ge 3$, $A \in \mathrm{St}(n,3)$.  If some pair satisfies*
> $$\ell_i + \ell_j \;\ge\; T^* := \frac{4n+2}{3n},$$
> *then a $1/n$-good triple through $\{i,j\}$ exists.*
>
> **Status: PROVED** (unconditional).

*Proof.*  Let $T = \ell_i + \ell_j = \operatorname{Tr}\Gamma = \mu_1 +
\mu_2 \ge T^*$.

*Step 1 ($\mu_2$ qualifies).*  By (M0), $\mu_1 \le 1$, so
$$\mu_2 = T - \mu_1 \;\ge\; T^* - 1 = \frac{4n+2-3n}{3n} = \frac{n+2}{3n}
\;>\; \frac1n \qquad (n + 2 > 3 \text{ for } n \ge 2).$$

*Step 2 (reduce to the threshold $T^*$).*  Both $\mu_1$ and $T^* - \mu_1$
lie in $(\frac1n,\infty)$: indeed $\mu_1 \ge \mu_2 > \frac1n$, and
$T^* - \mu_1 \ge T^* - 1 = \frac{n+2}{3n} > \frac1n$.  Since $h$ is strictly
decreasing on $(\frac1n,\infty)$ (M4.1) and $\mu_2 = T - \mu_1 \ge T^* -
\mu_1$,
$$h(\mu_1) + h(\mu_2) \;\le\; h(\mu_1) + h(T^* - \mu_1) =: f(\mu_1).$$

*Step 3 (convexity, endpoint evaluation).*  The value $\mu_1$ lies in
$[\,T^*/2,\ 1\,]$: $\mu_1 \ge \frac{\mu_1+\mu_2}2 = \frac T2 \ge
\frac{T^*}2$, and $\mu_1 \le 1$ by (M0); it suffices to bound $f$ on this
interval.  On this interval both arguments
$\mu_1$ and $T^*-\mu_1$ stay in $(\frac1n,\infty)$ (Step 2 and $T^*/2 =
\frac{2n+1}{3n} > \frac1n$), and $f$ is a sum of compositions of the convex
function $h$ (M4.1) with affine maps, hence convex.  A convex function on a
closed interval attains its maximum at an endpoint.  By (M4.3):
$$f\Bigl(\frac{T^*}2\Bigr) = 2\,h\Bigl(\frac{2n+1}{3n}\Bigr) = \frac1n,
\qquad
  f(1) = h(1) + h\Bigl(\frac{n+2}{3n}\Bigr) = 0 + \frac2n = \frac2n .$$
Hence $f \le \frac2n$ on $[T^*/2, 1]$, so $h(\mu_1)+h(\mu_2) \le \frac2n$,
and Lemma E (its hypothesis now verified: Step 1 and this bound) yields the
good triple through $\{i,j\}$. $\blacksquare$

**Remark (where the bound is tight).**  On the exact threshold $T = T^*$,
$f$ is convex on $[T^*/2,1]$ with $f(T^*/2) = \frac1n < \frac2n = f(1)$, so
$f(\mu_1) \le \frac1n + \frac{\mu_1 - T^*/2}{1 - T^*/2}\cdot\frac1n <
\frac2n$ strictly for $\mu_1 < 1$: the trace route saturates the h-curve
**only** at $\mu_1 = 1$.  Geometrically, $\mu_1 = 1$ with unit eigenvector
$y$ of $\Gamma$ means $u := y_1r_i + y_2r_j$ is a unit vector ($\|u\|^2 =
y^T\Gamma y = 1$); then $u\cdot r_i = (\Gamma y)_1 = y_1$, $u\cdot r_j =
y_2$, so $\sum_{p\notin\{i,j\}}(u\cdot r_p)^2 = u^T\bigl(\sum_p
r_pr_p^T\bigr)u - y_1^2 - y_2^2 = \|u\|^2 - 1 = 0$: all other rows are
orthogonal to $u$.  This matches the validated extremal heavy pairs, which
have $\mu_1 = 1$ exactly (§8.3).

### 4.1 Handy spectral special cases

> **Corollary S.**  *Let $n \ge 3$, $A \in \mathrm{St}(n,3)$.  Each of the
> following also implies a $1/n$-good triple through $\{i,j\}$:*
> *(a) $\mu_2 \ge \frac{n+1}{2n}$;\quad (b) $\mu_1 = 1$ and $\mu_2 \ge
> \frac{n+2}{3n}$.*
>
> **Status: PROVED.**

*Proof.*  (a) $\frac{n+1}{2n} > \frac1n$ for $n \ge 2$, so $\mu_2$
qualifies; $h$ decreasing gives $h(\mu_1) \le h(\mu_2) \le
h\bigl(\frac{n+1}{2n}\bigr) = \frac1n$ (M4.3), so the h-sum is $\le \frac2n$.
(b) $h(\mu_1) = h(1) = 0$ and $h(\mu_2) \le \frac2n$ by (M4.3).  In both
cases Lemma E applies. $\square$

---

## 5. Near-unit-row branch (Proposition U)

> **Proposition U.**  *Let $n \ge 4$, $A \in \mathrm{St}(n,3)$, and suppose
> (after relabeling) the largest leverage satisfies*
> $$\ell_1 = \ell_{\max} \;\ge\; 1 - \varepsilon_n, \qquad
>   \varepsilon_n := \frac{n-2}{n^2(n-1)} .$$
> *Then a $1/n$-good triple (containing row 1) exists, i.e. GTZ$(n,3)$ holds
> at $A$.*
>
> **Status: PROVED CONDITIONAL ON** the Sengupta–Pautov $k=2$ theorem
> (arXiv:2604.05944, "Hypothesis 1 for $k = 2$", proved in its Section 2),
> independently verified in `proofs/sp-verification.md` (VAL-SP-001) and
> trusted per AGENTS.md rule 4.  No other external input.

*Proof.*  Write $\varepsilon = 1 - \ell_1$, so $0 \le \varepsilon \le
\varepsilon_n$.  Note $\varepsilon_n \le \frac12$ for all $n \ge 3$ (indeed
$\varepsilon_n \le \frac1{24}$ for $n \ge 4$), so $\ell_1 \ge \frac12 > 0$.

*Step 1 (project off the heavy row).*  Let $\hat r = r_1/\sqrt{\ell_1}$
(unit vector) and $\Pi = I_3 - \hat r\hat r^T$, the orthoprojector onto the
plane $\hat r^\perp \subset \mathbb R^3$.  Set $\rho_p = \Pi r_p$ for
$p = 2,\dots,n$.  Then $\Pi r_1 = r_1 - \hat r\,(\hat r\cdot r_1) = r_1 -
\sqrt{\ell_1}\,\hat r = 0$, and applying $\Pi(\cdot)\Pi$ to
$\sum_{p=1}^n r_pr_p^T = I_3$ gives
$$\sum_{p=2}^n \rho_p\rho_p^T = \Pi I_3\Pi = \Pi .$$
Fix an orthonormal basis $u_1, u_2$ of $\hat r^\perp$ and let $A' \in
\mathbb R^{(n-1)\times 2}$ have rows $(\langle\rho_p, u_1\rangle,
\langle\rho_p, u_2\rangle)$, $p = 2,\dots,n$.  Then $(A'^TA')_{ab} =
u_a^T\bigl(\sum_{p\ge2}\rho_p\rho_p^T\bigr)u_b = u_a^T\Pi u_b =
\delta_{ab}$, i.e. $A' \in \mathrm{St}(n-1, 2)$ **exactly**.

*Step 2 (apply $k = 2$).*  Since $n - 1 \ge 3$, the Sengupta–Pautov theorem
applies to $A'$ (its base case $n = 3$ is covered: a self-contained proof is
supplied in `proofs/sp-verification.md` §3.1): there exists a pair
$\{i,j\} \subseteq \{2,\dots,n\}$ with
$\sigma_{\min}(A'_{\{i,j\}}) \ge 1/\sqrt{n-1}$, i.e., writing $G_\rho$ for
the $2\times2$ Gram matrix of $(\rho_i, \rho_j)$,
$$G_\rho = A'_{\{i,j\}}A_{\{i,j\}}'^{\,T} \succeq \frac1{n-1}\,I_2 .$$
(Gram matrices are basis-independent: $\langle\rho_i,\rho_j\rangle$ in
$\mathbb R^3$ equals the coordinate inner product in the basis $u_1,u_2$.)

*Step 3 (lift: the triple $\{1,i,j\}$).*  Let $a_p = \hat r\cdot r_p$, so
$r_p = \rho_p + a_p\hat r$ is the orthogonal split along
$\hat r^\perp \oplus \mathbb R\hat r$.  From $\sum_p r_pr_p^T = I_3$,
$$\sum_{p=1}^n a_p^2 = \hat r^T I_3 \hat r = 1, \qquad a_1^2 = \ell_1,
\qquad\text{so}\quad \sum_{p=2}^n a_p^2 = 1 - \ell_1 = \varepsilon .$$
Put $a = (a_i, a_j)^T$; then $\|a\|^2 \le \varepsilon$.  The Gram data of
the triple $\{1,i,j\}$ (order $(i,j,1)$):
$$\Gamma_{ij} := \begin{pmatrix}\ell_i & c_{ij}\\ c_{ij} & \ell_j
\end{pmatrix} = G_\rho + aa^T \;\succeq\; G_\rho, \qquad
  w := (c_{1i}, c_{1j})^T = \sqrt{\ell_1}\,a,$$
using $c_{ij} = \langle\rho_i,\rho_j\rangle + a_ia_j$, $\ell_i =
\|\rho_i\|^2 + a_i^2$, and $c_{1p} = \langle\sqrt{\ell_1}\hat r,\ \rho_p +
a_p\hat r\rangle = \sqrt{\ell_1}\,a_p$.

*Step 4 (Schur).*  From Steps 2–3,
$$\Gamma_{ij} - \frac1nI_2 \;\succeq\; \Bigl(\frac1{n-1} -
\frac1n\Bigr)I_2 = \frac1{n(n-1)}\,I_2 \;\succ\; 0,$$
so (M2) applies to the pair $\{i,j\}$ with third row $1$; for a symmetric
$X \succeq cI \succ 0$ one has $X^{-1} \preceq c^{-1}I$ (spectral mapping),
hence
$$w^T\Bigl(\Gamma_{ij} - \tfrac1nI_2\Bigr)^{-1}w \;\le\; n(n-1)\,\|w\|^2
 = n(n-1)\,\ell_1\|a\|^2 \;\le\; n(n-1)\,\ell_1\,\varepsilon .$$
Thus $q(1) = \ell_1 - \frac1n - w^T(\Gamma_{ij}-\frac1nI_2)^{-1}w \ge 0$ is
implied by
$$n(n-1)\,\ell_1\,\varepsilon \;\le\; \ell_1 - \frac1n
\quad\Longleftrightarrow\quad
  n(n-1)\,\varepsilon \;\le\; 1 - \frac1{n\ell_1}
  \qquad(\text{divide by } \ell_1 > 0).$$
Since $\ell_1 \ge \frac12$, $\frac1{n\ell_1} \le \frac2n$, so it suffices
that
$$n(n-1)\,\varepsilon \;\le\; 1 - \frac2n = \frac{n-2}{n}
\quad\Longleftrightarrow\quad
  \varepsilon \;\le\; \frac{n-2}{n^2(n-1)} = \varepsilon_n,$$
which is the hypothesis.  Then $q(1) \ge 0$ and (M2) gives
$\lambda_{\min}(\mathrm{Gram}\{1,i,j\}) \ge \frac1n$. $\blacksquare$

*Degenerate check:* $\varepsilon = 0$ forces $a_p = 0$ for $p \ge 2$, so
$w = 0$, $q(1) = 1 - \frac1n > 0$ — covered by the argument.
(For $n = 3$ the statement is vacuous-trivial: $A \in O(3)$, every triple
Gram is $I_3$.)  Numeric check F5: 1500 planted configs, $n = 4..10$,
$0$ violations, worst margin $+9.7\cdot10^{-2}$.

---

## 6. Dual branch at $(6,3)$ (Proposition D)

> **Proposition D.**  *Let $A \in \mathrm{St}(6,3)$ and suppose some pair
> satisfies $\ell_i + \ell_j \le \frac59$ (equivalently: the two smallest
> leverages sum to $\le \frac59$).  Then GTZ$(6,3)$ holds at $A$: some
> triple $I$ has $\sigma_{\min}(A_I) \ge 1/\sqrt6$.*
>
> **Status: PROVED CONDITIONAL ON** `proofs/duality.md` Theorem 1 /
> Corollary 1 (VAL-DUAL-001, validated in-project); the other ingredient is
> Corollary T of this document (unconditional).

*Proof.*  Extend the orthonormal columns of $A$ to an orthonormal basis of
$\mathbb R^6$: $[A \mid B] \in O(6)$ with $B \in \mathrm{St}(6,3)$.  Then
$BB^T = I_6 - AA^T = I - P$, so the leverages of $B$ are
$\ell^B_p = 1 - \ell_p$.  The hypothesis gives
$$\ell^B_i + \ell^B_j = 2 - (\ell_i + \ell_j) \;\ge\; 2 - \frac59 =
\frac{13}9 = \frac{4\cdot6+2}{3\cdot6} = T^*\big|_{n=6}.$$
Corollary T applied to $B \in \mathrm{St}(6,3)$ produces a triple $J \supseteq
\{i,j\}$ with $\lambda_{\min}(\mathrm{Gram}_B(J)) \ge \frac16$, i.e.
$\sigma_{\min}(B_J) \ge 1/\sqrt6$.  By `proofs/duality.md` Corollary 1(a)
(applied to $[A\mid B] \in O(6)$ at $n = 6$, $k = 3 = n - k$, with the
triple $I = J^c$): $\sigma_{\min}(A_{J^c}) = \sigma_{\min}(B_J) \ge
1/\sqrt6$, and $|J^c| = 3$.  The triple $J^c$ witnesses GTZ$(6,3)$ at $A$.
$\blacksquare$

*Remark (scope of the trick).*  Self-duality of the triple structure is
specific to $n = 2k = 6$: for general $(n,3)$ the dual lives in
$\mathrm{St}(n, n-3)$, where the complement of a triple is not a triple.
Numeric check F6: 1500 planted $(6,3)$ configs with $\ell_0 + \ell_1 \le
5/9$, 0 violations, worst GTZ margin $+8.2\cdot10^{-2}$.

---

## 7. Obstruction: the h-branch is empty at the icosahedral ETF$(6,3)$

### 7.1 The configuration, exactly

Let $\varphi = \frac{1+\sqrt5}2$ (so $\varphi^2 = \varphi + 1$).  Take one
vertex from each antipodal pair of a regular icosahedron:
$$v_1 = (0,1,\varphi),\ v_2 = (0,1,-\varphi),\ v_3 = (1,\varphi,0),\
  v_4 = (1,-\varphi,0),\ v_5 = (\varphi,0,1),\ v_6 = (-\varphi,0,1),$$
and set $a_p = v_p/\sqrt{2(1+\varphi^2)}$, $A_{\mathrm{ico}} \in
\mathbb R^{6\times3}$ with rows $a_p$.

**Facts (all exact).**

1. $\sum_p v_pv_p^T = 2(1+\varphi^2)\,I_3$: the diagonal of the sum is
   $(2 + 2\varphi^2)(1,1,1)$ and every off-diagonal entry cancels in pairs
   (direct computation: e.g. the $(1,2)$-entries are $0,0,\varphi,-\varphi,
   0,0$).  Hence $A_{\mathrm{ico}}^TA_{\mathrm{ico}} = I_3$:
   $A_{\mathrm{ico}} \in \mathrm{St}(6,3)$.
2. $\ell_p = \frac{1+\varphi^2}{2(1+\varphi^2)} = \frac12$ for all $p$.
3. Every one of the 15 pairs has $v_p\cdot v_q = \pm\varphi$ (direct
   inspection of the 15 inner products), so, using $1 + \varphi^2 = 2 +
   \varphi$ and $(2+\varphi)^2 = 4 + 4\varphi + \varphi^2 = 5 + 5\varphi =
   5\varphi^2$:
   $$c_{pq}^2 = \frac{\varphi^2}{4(1+\varphi^2)^2}
   = \frac{\varphi^2}{4(2+\varphi)^2} = \frac{\varphi^2}{20\,\varphi^2}
   = \frac1{20} \qquad \text{for all } p \ne q.$$
   (This is the equiangular tight frame ETF$(6,3)$; equality holds in the
   Welch bound, not needed below.)

**Status: PROVED** (exact; symbolic check S7 re-verifies 1.–3. in sympy).

### 7.2 Pair spectra, the obstruction, and the coherent-triple rescue

**Obstruction.**  *Every pair of $A_{\mathrm{ico}}$ has Gram
$\Gamma = \begin{pmatrix}1/2 & c\\ c & 1/2\end{pmatrix}$ with $c^2 =
\frac1{20}$, hence eigenvalues*
$$\mu_{1,2} = \frac12 \pm \frac1{\sqrt{20}} = \frac{5\pm\sqrt5}{10},$$
*both $> \frac16$ (so every pair qualifies for the h-test), and*
$$h(\mu_1) + h(\mu_2) = \frac{13-5\sqrt5}{22} + \frac{13+5\sqrt5}{22}
 = \frac{13}{11} \;>\; \frac13 = \frac2n\Big|_{n=6}.$$
*The h-branch is **empty** at the icosahedral ETF: Lemma E applies to no
pair, and (by (M3)) for every pair even the full Schur-slack **sum** is
negative — no argument using only a single pair's second moments can
produce the good triple here.*

*Proof.*  The eigenvalues of a $2\times2$ matrix with equal diagonal
$\frac12$ and off-diagonal $c$ are $\frac12 \pm |c| = \frac12 \pm
\frac1{\sqrt{20}}$.  Both exceed $\frac16$: $\frac12 - \frac1{\sqrt{20}} >
\frac16 \iff \frac13 > \frac1{\sqrt{20}} \iff 20 > 9$.  With $n = 6$,
$h(\mu) = \frac{1-\mu}{6\mu-1}$:
$$h\Bigl(\frac{5+\sqrt5}{10}\Bigr)
 = \frac{(5-\sqrt5)/10}{(20+6\sqrt5)/10}
 = \frac{5-\sqrt5}{2(10+3\sqrt5)}
 = \frac{(5-\sqrt5)(10-3\sqrt5)}{2(100-45)}
 = \frac{65-25\sqrt5}{110} = \frac{13-5\sqrt5}{22},$$
and symmetrically $h\bigl(\frac{5-\sqrt5}{10}\bigr) =
\frac{13+5\sqrt5}{22}$; the sum is $\frac{26}{22} = \frac{13}{11}$.
Finally $\frac{13}{11} > \frac13$ since $39 > 11$.  The last sentence is
(M3): h-sum $> \frac2n$ $\iff$ $\sum_p q(p) < 0$. $\blacksquare$

**Status: PROVED** (exact; symbolic check S7: h-sum $= 13/11$ for all 15
pairs).

**Rescue (GTZ still holds at the ico — via sign coherence, not via
Lemma E).**  *Call a triple $\{p,q,r\}$ of $A_{\mathrm{ico}}$ **coherent**
if $c_{pq}c_{pr}c_{qr} > 0$.  Then:*

1. *A coherent triple has Gram spectrum
   $\bigl\{\frac12 + \frac2{\sqrt{20}},\ \frac12 - \frac1{\sqrt{20}},\
   \frac12 - \frac1{\sqrt{20}}\bigr\}$, so*
   $$\lambda_{\min} = \frac12 - \frac1{\sqrt{20}} = \frac{5-\sqrt5}{10}
     \;>\; \frac16; $$
   *an incoherent triple has spectrum $\{\frac12 - \frac2{\sqrt{20}},\
   \frac12 + \frac1{\sqrt{20}},\ \frac12 + \frac1{\sqrt{20}}\}$ with
   $\lambda_{\min} = \frac{5-2\sqrt5}{10} < \frac16$.  So a triple of
   $A_{\mathrm{ico}}$ is $\frac16$-good **iff** it is coherent.*
2. *Coherent triples exist — e.g. $\{1,3,5\}$: $v_1\cdot v_3 = v_1\cdot v_5
   = v_3\cdot v_5 = \varphi > 0$.  Hence GTZ$(6,3)$ holds at
   $A_{\mathrm{ico}}$ with margin: $\frac{5-\sqrt5}{10} \approx 0.276 >
   \frac16$.*
3. *(Finite verification.)  Exactly 10 of the 20 triples are coherent, and
   every one of the 15 pairs lies in exactly 2 coherent triples.*

*Proof.*  1.  The triple Gram is $\frac12I_3 + C$ with $C$ symmetric,
zero-diagonal, off-diagonal entries $c_{pq} = \pm\frac1{\sqrt{20}}$.
Conjugating by $D = \mathrm{diag}(\epsilon_p, \epsilon_q, \epsilon_r)$,
$\epsilon \in \{\pm1\}^3$ (an orthogonal matrix; preserves eigenvalues and
the sign product $c_{pq}c_{pr}c_{qr} \mapsto
\epsilon_p^2\epsilon_q^2\epsilon_r^2\,c_{pq}c_{pr}c_{qr}$), one can
normalize the sign pattern: if the product is positive the pattern is
$(+,+,+)$ or has exactly two minus signs, and flipping the row shared by
the two negative entries yields $(+,+,+)$; if negative, the pattern has one
or three minus signs, and flipping the row shared by the two positive
entries (or no flip) yields $(-,-,-)$.  With all off-diagonals equal to
$c = \frac1{\sqrt{20}}$: the matrix is $(\frac12 - c)I_3 + cE$ ($E$ =
all-ones, spectrum $\{3,0,0\}$), giving eigenvalues $\frac12 + 2c,\
\frac12 - c,\ \frac12 - c$; with all equal to $-c$, eigenvalues
$\frac12 - 2c,\ \frac12 + c,\ \frac12 + c$.  The comparisons with
$\frac16$: $\frac12 - \frac1{\sqrt{20}} > \frac16 \iff 20 > 9$ (above), and
$\frac12 - \frac2{\sqrt{20}} < \frac16 \iff \frac13 < \frac2{\sqrt{20}}
\iff 20 < 36$.
2.  Direct computation of the three inner products; then apply 1.
3.  Finite check over all 20 triples, performed in exact arithmetic
(symbolic check S7: all 20 characteristic polynomials match the two exact
spectra above; counts $10$ coherent, per-pair multiplicity exactly $2$).
$\blacksquare$

**Status: PROVED** (exact + finite verification).

**Consequence for the bridge.**  The obstruction shows Lemma E **cannot**
be the whole pair → triple bridge: at $A_{\mathrm{ico}}$ every pair fails
the h-test, yet GTZ holds there through a different mechanism (sign
coherence).  Any complete proof of the bridge must supply a second,
diffuse-regime argument (§8.4).

---

## 8. The equality manifold, the conjecture C_eq, and what remains open

### 8.1 Strict interior of the h-branch

> **Corollary E′.**  *Let $n \ge 3$, $A \in \mathrm{St}(n,3)$.  If some pair
> has $\mu_2 > \frac1n$ and $h(\mu_1) + h(\mu_2) < \frac2n$ (strict), then
> some triple through it has $\lambda_{\min}(\mathrm{Gram}) > \frac1n$
> strictly; in particular $f(A) = \max_{|I|=3}\sigma_{\min}(A_I) >
> 1/\sqrt n$.*
>
> **Status: PROVED.**

*Proof.*  By (M3) the strict inequality gives $\sum_{p\notin\{i,j\}}q(p) >
0$, so some $q(p) > 0$, and the strict part of (M2) gives
$\lambda_{\min} > 1/n$. $\square$

So GTZ has *slack* strictly inside the h-branch; extremal configurations
(where $f = 1/\sqrt n$ exactly) can only sit on its boundary or outside.

### 8.2 The equality manifold and two exact families on it

Define the **equality manifold**
$$\mathcal E_n = \Bigl\{A \in \mathrm{St}(n,3):\ \min_{\{i,j\}:\,\mu_2>1/n}
\bigl[h(\mu_1)+h(\mu_2)\bigr] = \tfrac2n\Bigr\}$$
— the boundary of the h-branch.  Two one-parameter (in $n$) spectra lie
**exactly** on the curve $h(\mu_1)+h(\mu_2) = \frac2n$:
$$\text{E-I} = \Bigl(1,\ \frac{n+2}{3n}\Bigr), \qquad
  \text{E-II} = \Bigl(\frac{n^2-n-4}{2n(n-3)},\ \frac12\Bigr)\ (n \ge 4).$$

*Proof.*  E-I: $h(1) = 0$ and $h\bigl(\frac{n+2}{3n}\bigr) = \frac2n$
(M4.3).  E-II: $h(\frac12) = \frac{1/2}{n/2-1} = \frac1{n-2}$, and solving
$h(\mu) = \frac2n - \frac1{n-2} = \frac{n-4}{n(n-2)}$ via
$\mu = \frac{1+c}{1+cn}$ (the inverse of $h$, valid for $c > -\frac1n$)
gives $\mu = \frac{n^2-n-4}{2n(n-3)}$.  (Symbolic check S4; at $n = 6$
these are $(1, \frac49)$ and $(\frac{13}{18}, \frac12)$.) $\square$

**Status: PROVED** (the curve membership of E-I/E-II — a pure algebraic
identity).

### 8.3 Conjecture C_eq — CONJECTURED, not proved

> **Conjecture C_eq [CONJECTURED — numerically supported, NOT proved].**
> Every GTZ-extremal $(n,3)$ configuration (i.e. $f(A) = 1/\sqrt n$) lies on
> $\mathcal E_n$, with its minimizing pair spectrum at E-I or E-II; the open
> h-branch contains no extremal configuration (this last part **is** a
> theorem — Corollary E′); and $f$ admits a positive gap above $1/\sqrt n$
> uniformly outside any neighborhood of the h-branch.

Supporting data (validated catalogs; *evidence, not proof*): all 996
catalog $(6,3)$ near-extremal matrices (VAL-NUM-002) have min-pair h-sum
$= \frac13$ to within the catalog's own optimizer error ($\le 6.4\cdot
10^{-5}$), never strictly inside; all $k=3$ sweep argmins $n = 6..10$
(VAL-NUM-003) are on the manifold to $\le 5\cdot10^{-7}$ with equality-pair
spectra at E-I ($n = 6, 10$) or E-II ($n = 7,8,9$); catalog classes 0/1
have their heavy pairs at $\ell_1+\ell_2 = \frac{13}9 = T^*|_{n=6}$
exactly, and the equality case of Lemma E (§3.1) reproduces class 0's
observed 10 active triples.  Adversarial forced-diffuse minimization
(round-2 ledger §10(B)) could not push $f$ below $1/\sqrt6 +
1.3\cdot10^{-2}$.  See `attacks/round-2-pair-to-triple.md` §6, §10.

### 8.4 What remains OPEN

- **The bridge C_ext / GTZ$(n,3)$ itself: OPEN.**  This document covers the
  closed h-branch (Lemma E), the trace branch (Corollary T $\Rightarrow$
  h-branch), a near-unit-row sliver (Proposition U), and — at $(6,3)$
  only — a small-leverage-pair branch (Proposition D).  The **diffuse
  regime** (every pair with $\mu_2 > \frac1n$ strictly above the h-curve,
  e.g. the icosahedral ETF) is untouched by all of them, and §7 proves it
  is nonempty and cannot be reached by single-pair second-moment arguments.
- The naive diffuse selection rule "the min-h-sum pair extends" is
  **adversarially refuted** (round-2 ledger §10(A): margin
  $-2.8\cdot10^{-2}$ deep in the diffuse region, harmless to Lemma E whose
  hypothesis fails there) — pair multiplicity or a sign-coherence mechanism
  (§7.2) is needed.
- C_eq (§8.3) is CONJECTURED.  Nothing above uses it.

---

## 9. Numeric sanity block

All checks: `numerics/study-extension-lemma/check_extension.py`
(deterministic, seed 20260718; run from the workspace root with
`.venv/bin/python`; full output in `check_extension.log`, **ALL CHECKS
PASSED, exit 0**).  Claim-by-claim map:

| Claim | Check | Observed |
|---|---|---|
| $g$ closed form, $h$ formula, $h' < 0$, $h''$ (M4) | S1 | exact (sympy) |
| $h$ special values (M4.3) | S2 | exact |
| Corollary T endpoints $f(1) = \frac2n$, $f(T^*/2) = \frac1n$ | S3 | exact |
| E-I, E-II on the curve; $n=6$ instances $(1,\frac49), (\frac{13}{18},\frac12)$ | S4 | exact |
| Prop U bookkeeping: $n(n-1)\varepsilon_n = \frac{n-2}n$, $\varepsilon_n \le \frac12$ | S5 | exact |
| Prop D threshold $2 - \frac59 = \frac{13}9 = T^*(6)$ | S6 | exact |
| Ico: $A^TA = I_3$, $\ell \equiv \frac12$, $c^2 \equiv \frac1{20}$, h-sum $\equiv \frac{13}{11} > \frac13$; 10 coherent triples, each pair in 2; exact spectra of all 20 triples | S7 | exact (sympy) |
| Moment identity $\sum_{p\ne i,j} v_pv_p^T = \Gamma - \Gamma^2$ | F1 | worst err $5.0\cdot10^{-16}$ (205 Haar) |
| q-sum identity ($k=3$) | F1 | worst err $1.8\cdot10^{-11}$ (conditioning-limited) |
| General-$k$ q-sum, $k = 2..6$ | F2 | worst err $4.2\cdot10^{-13}$ (375 samples) |
| $k=2$ infeasibility (h-sum $\ge 0$; $\sum q < 0$ always) | F2 | max $\sum q = -0.559$ (300 Haar) |
| Lemma E conclusion | F3 | 2906 planted (boundary-stressed) + 294 Haar h-branch, 0 violations, worst margin $+2.9\cdot10^{-2}$ |
| Corollary T conclusion | F4 | 3000 planted, 0 violations, worst margin $+6.4\cdot10^{-2}$ |
| Proposition U conclusion | F5 | 1500 planted, $n = 4..10$, 0 violations, worst margin $+9.7\cdot10^{-2}$ |
| Proposition D conclusion | F6 | 1500 planted $(6,3)$, 0 violations, worst margin $+8.2\cdot10^{-2}$ |

Independent earlier evidence (different code path, seed 20260718):
`numerics/study-attack-r2-bridge/` — identities, planted/Haar Lemma E and
Corollary T checks, catalog equality-manifold survey, adversarial probes;
see its README and `attacks/round-2-pair-to-triple.md`.

---

## 10. Provenance and citations

- Promoted from `attacks/round-2-pair-to-triple.md` (attack round 2);
  proofs reconstructed and completed here — this document is
  self-contained modulo the two explicitly conditional inputs.
- External/in-project inputs: **Proposition U only** uses arXiv:2604.05944
  ($k=2$ theorem; verified in `proofs/sp-verification.md`); **Proposition D
  only** uses `proofs/duality.md` (Theorem 1 / Corollary 1).  Lemma E,
  Corollaries T/S/E′, and §7 are unconditional and self-contained.
- Standard facts used by name: Sylvester's law of inertia (M2); convex
  function on an interval attains its maximum at an endpoint (Corollary T);
  spectral mapping $X \succeq cI \succ 0 \Rightarrow X^{-1} \preceq
  c^{-1}I$ (Proposition U).
