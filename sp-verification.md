# Independent verification of arXiv:2604.05944 (Sengupta–Pautov): GTZ for k = 2, arbitrary n

**Source verified:** R. Sengupta, M. Pautov, *On the submatrices with the best-bounded
inverses*, arXiv:2604.05944 (LaTeX source fetched from
`https://arxiv.org/e-print/2604.05944` on 2026-07-17; single file
`richik_and_misha.tex`, dated 2026-04-16).

**Claim verified (the paper's Hypothesis 1, proved there for k = 2):**
for every $n > 2$ and every $A \in \mathbb{R}^{n\times 2}$ with $A^T A = I_2$,
there exists a $2\times 2$ submatrix $\tilde A$ (two rows of $A$) with
$\sigma_2(\tilde A) \ge \tfrac{1}{\sqrt n}$, equivalently
$\|\tilde A^{-1}\|_2 \le \sqrt n$.

**Method:** every step of the paper is reproduced below with all omitted
computations filled in; every numerically checkable claim is tested by
`numerics/study-sp-verification/check_sp.py` (seeds recorded; all checks pass —
see the Numeric sanity block).

**Notation** (shared project notation, MEMORY.md): rows $r_i = (x_i, y_i) \in
\mathbb{R}^2$; leverages $\ell_i = \|r_i\|^2 = x_i^2 + y_i^2$; $\sigma_1 \ge
\sigma_2$ the singular values of a $2\times2$ matrix. The paper writes
$A_{i1}, A_{i2}$ for $x_i, y_i$.

---

## 1. Verdict

### Overall: **CORRECT WITH FIXABLE GAPS** — every gap is filled in this report; none threatens the result.

The theorem **GTZ(n, k=2) for all n > 2 is PROVED** by the paper modulo the
gaps G1–G5 listed below, all of which are routine and are closed in full in
§§3–7. After the fills, the status of the theorem is **PROVED**. Downstream
consumers (VAL-DUAL-002 and anything relying on the k = 2 case) lose nothing.

### Per-component verdict

| Component | Verdict | Gaps (all filled here) |
|---|---|---|
| **Case A** (rotation, row deletion, renormalization $t^2 = \frac{1}{1-b^2}$, constant self-reproduction) | **CORRECT** | G2: $\sigma_2(YZ) \ge \sigma_2(Y)\sigma_2(Z)$ used without proof (standard; proved in §4.4) |
| **Case B, pair ⇒ bound** (pair inequality ⇒ char. poly of Gram ≥ 0 at $\frac1n$ ⇒ $\sigma_2^2 \ge \frac1n$) | **CORRECT** | G5: index typo "$\|r_i\|^2 + \|r_k\|^2$" for "$\|r_i\|^2 + \|r_j\|^2$" (trivial) |
| **Pair existence** (squaring map, eigenvalue count of $G = WW^T - zz^T$, $\operatorname{Tr} G = \frac4n$, Perron–Frobenius contradiction) | **CORRECT** | G3: Perron step glosses over why $\lambda_1(G)$, as an eigenvalue of $M$, cannot equal $\frac2n$, and why an eigenvalue with positive eigenvector is the Perron root (both filled in §6.5); G4: the excluded case $i = j$ is handled correctly (verified, §6.6) |
| **Induction base** | **CORRECT WITH FIXABLE GAP** | G1: base case $n = 3$ is asserted "trivial" with no proof, and $n = 4$ is outsourced to Nesterenko (arXiv:2303.07492). Both filled in §3: a self-contained proof of $n = 3$ is supplied, and it is shown that the induction step itself already covers $n = 4$ from the $n=3$ base, so the external dependency on Nesterenko is removable. |

No FATAL step was found. The logical skeleton (exhaustive Case A / Case B split,
strong induction on $n$, contradiction via Perron–Frobenius) is sound, and all
constants track exactly — in particular the bound self-reproduces without loss:
$\sqrt{\tfrac{n-1}{n}} \cdot \tfrac{1}{\sqrt{n-1}} = \tfrac{1}{\sqrt n}$ in Case
A, and $\tfrac1n$ appears exactly as the evaluation point of the characteristic
polynomial in Case B.

---

## 2. Statement and proof skeleton

**Theorem (paper's Hypothesis 1 for k = 2).** For every $n > 2$ and every
$A \in \mathbb{R}^{n\times2}$ with $A^TA = I_2$ there are indices $i < j$ such
that $\sigma_2\begin{pmatrix} r_i \\ r_j\end{pmatrix} \ge \frac{1}{\sqrt n}$.

**Skeleton of the paper's proof.** Induction on $n$. Bases: $n = 3$
("trivial"), $n = 4$ (citing Nesterenko). Induction step for a fixed $n$
(assuming the statement for $n-1$), by the exhaustive dichotomy:

- **Case A:** some row has $\ell_i \le \frac1n$. Reduce to an $(n-1)\times2$
  orthonormal-column matrix by a rotation, a row deletion and a column
  renormalization; the constant self-reproduces.
- **Case B:** all rows have $\ell_i > \frac1n$. No induction used: find a pair
  $i \ne j$ with $(r_i, r_j)^2 \le (\ell_i - \frac1n)(\ell_j - \frac1n)$
  (existence via the squaring map and a Perron–Frobenius contradiction), then
  read the bound off the characteristic polynomial of the $2\times2$ Gram
  matrix.

The dichotomy is exhaustive: either $\min_i \ell_i \le \frac1n$ (Case A) or
$\min_i \ell_i > \frac1n$ (Case B). Each of §§3–7 below verifies one component.

---

## 3. Induction base — verdict CORRECT WITH FIXABLE GAP (gap filled)

**What the paper says:** "the statement for $n=3$ is trivial and the proof for
$n=4$ is due to Y. Nesterenko." No proof of either is given.

**Gap G1.** "Trivial" is not a proof, and outsourcing $n = 4$ imports an
external dependency (note: MEMORY.md records a typo in system (2.4) of
arXiv:2303.07492, so blind reliance is undesirable). Both are fixable; two
independent fills follow.

### 3.1 Fill 1: direct proof for n = 3

Let $A \in \mathbb{R}^{3\times2}$, $A^TA = I_2$. Complete $A$ to an orthogonal
matrix $Q = [A \mid c] \in O(3)$, where $c \in \mathbb{R}^3$ is a unit vector
orthogonal to both columns of $A$ (exists: the two columns span a plane in
$\mathbb{R}^3$; take $c$ a unit normal). Let $A_{\hat i}$ denote the $2\times2$
submatrix of $A$ obtained by deleting row $i$.

*Claim: $|\det A_{\hat i}| = |c_i|$.* Since $Q^{-1} = Q^T$ and
$Q^{-1} = \frac{1}{\det Q}\operatorname{adj}(Q)$ with $\det Q = \pm1$, the
$(3,i)$ entry of $\operatorname{adj}(Q)$ equals $(\det Q)\, (Q^T)_{3i} = (\det
Q)\, c_i$. By definition of the adjugate,
$\operatorname{adj}(Q)_{3i} = (-1)^{3+i} M_{i3}$, where $M_{i3}$ is the minor
of $Q$ deleting row $i$ and column 3 — i.e. exactly $\det A_{\hat i}$. Taking
absolute values gives the claim.

Since $\sum_i c_i^2 = 1$, some $i$ has $c_i^2 \ge \frac13$. For that $i$:
$\sigma_1(A_{\hat i}) \le \sigma_1(A) = 1$ (a submatrix's largest singular value
is at most the whole matrix's, and $\sigma_1(A) = 1$ because $A^TA = I$), hence

$$\sigma_2(A_{\hat i}) = \frac{|\det A_{\hat i}|}{\sigma_1(A_{\hat i})} \ge |c_i| \ge \frac{1}{\sqrt3}.$$

$\blacksquare$ (Status: PROVED.)

### 3.2 Fill 2: the induction step already covers n = 4 (Nesterenko not needed)

The induction step (§§4–7) is valid for any $n \ge 4$ given the statement for
$n - 1 \ge 3$: Case B (§§5–6) uses no induction at all, and Case A (§4) applies
the hypothesis to an $(n-1)\times2$ orthonormal-column matrix, which for
$n - 1 = 3$ is the base of §3.1. So the citation of Nesterenko for $n=4$ is a
convenience, not a dependency: base $n = 3$ + the step give all $n \ge 3$.

(Remark, not needed for the verdict: the step even works from $n = 3$ down to
$n = 2$ if one grants the trivially true $n = 2$ statement — a $2\times2$
matrix with orthonormal columns is orthogonal and is its own submatrix with
$\sigma_2 = 1 \ge \frac{1}{\sqrt2}$. Under that convention no separate base
beyond $n=2$ is needed at all.)

**Component verdict: CORRECT WITH FIXABLE GAP — filled; base cases $n=3$
(and $n=4$ via the step) now stand proved self-containedly.**

---

## 4. Case A — verdict CORRECT

**Hypothesis of the case:** $\exists\, i:\ \ell_i \le \frac1n$; WLOG $i = 1$
(row permutations act on the left and permute the set of $2\times2$ row-pair
submatrices without changing their singular values, and preserve $A^TA = I$).

### 4.1 Rotation to (b, 0)

Choose a rotation $P \in SO(2)$ with $r_1 P = (b, 0)$, $b^2 = \ell_1 \le
\frac1n$. Explicitly, if $b > 0$, $P = \frac1b\begin{pmatrix} x_1 & -y_1 \\ y_1
& x_1\end{pmatrix}$; check: $r_1 P = \frac1b(x_1^2 + y_1^2,\; -x_1 y_1 + y_1
x_1) = (b, 0)$, and $P^TP = \frac{1}{b^2}(x_1^2 + y_1^2) I = I$. Set $B = AP$.
Then:

- $B^TB = P^T A^T A P = P^T P = I$ — orthonormal columns are preserved;
- for any rows $i, j$, the submatrix of $B$ is $\begin{pmatrix} r_i \\
  r_j\end{pmatrix} P$, and right multiplication by an orthogonal matrix
  preserves singular values ($(\tilde A P)(\tilde A P)^T = \tilde A \tilde
  A^T$). So proving the bound for $B$ proves it for $A$. ✓ (paper's claim,
  verified).

**Sub-case b = 0** (paper handles it separately, correctly): delete row 1 of
$B$; the remaining $(n-1)\times2$ matrix $C$ satisfies $C^TC = B^TB - (0,0)^T
(0,0) = I$. By the induction hypothesis it has a $2\times2$ submatrix with
$\sigma_2 \ge \frac{1}{\sqrt{n-1}} > \frac{1}{\sqrt n}$, and that submatrix is
a submatrix of $B$. ✓ Henceforth $0 < b^2 \le \frac1n$.

### 4.2 Row deletion and the three column identities

Let $C$ be rows $2..n$ of $B$. Writing out $B^TB = I$ and using $B_{11} = b$,
$B_{12} = 0$:

$$\sum_{i=2}^n B_{i1} B_{i2} = 0 - b\cdot 0 = 0, \qquad
\sum_{i=2}^n B_{i2}^2 = 1 - 0^2 = 1, \qquad
\sum_{i=2}^n B_{i1}^2 = 1 - b^2 .$$

All three match the paper's display (their eq. (7)). ✓

### 4.3 Column renormalization: t² = 1/(1 − b²)

Since $b^2 \le \frac1n < 1$, we have $1 - b^2 > 0$; set $t =
\frac{1}{\sqrt{1-b^2}} > 1$ (indeed $t^2 > 1$ because $b^2 > 0$ — paper's
remark, correct and *used later* for identifying $\sigma_2$ of the scaling
matrix). Multiply the first column of $C$ by $t$ to get $\hat C$. Then

$$\hat C^T \hat C = \begin{pmatrix} t^2 \sum_{i\ge2} B_{i1}^2 & t \sum_{i\ge2} B_{i1}B_{i2} \\ t \sum_{i\ge2} B_{i1}B_{i2} & \sum_{i\ge2} B_{i2}^2 \end{pmatrix} = \begin{pmatrix} t^2(1-b^2) & 0 \\ 0 & 1 \end{pmatrix} = I_2 . \checkmark$$

So $\hat C$ is $(n-1)\times2$ with orthonormal columns; the induction
hypothesis yields rows $i, j$ (indices in $2..n$) with

$$\tilde C = \begin{pmatrix} tB_{i1} & B_{i2} \\ tB_{j1} & B_{j2}\end{pmatrix}, \qquad \sigma_2^2(\tilde C) \ge \frac{1}{n-1}.$$

### 4.4 Back-scaling and constant self-reproduction (gap G2 filled)

$\tilde C = \tilde B \operatorname{diag}(t, 1)$ where $\tilde B =
\begin{pmatrix} B_{i1} & B_{i2} \\ B_{j1} & B_{j2}\end{pmatrix}$, i.e. $\tilde
B = \tilde C \operatorname{diag}(t^{-1}, 1)$.

**Gap G2 (fill).** The paper uses $\sigma_2(YZ) \ge \sigma_2(Y)\sigma_2(Z)$
for $2\times2$ matrices without proof. Proof (valid for any square $Y, Z$ of
equal size): for a unit vector $u$ realizing the minimum,
$$\sigma_2(YZ) = \min_{\|u\|=1} \|YZu\| \ge \min_{\|u\|=1} \sigma_2(Y)\, \|Zu\| = \sigma_2(Y)\, \sigma_2(Z),$$
using $\|Yv\| \ge \sigma_2(Y)\|v\|$ for every $v$ (with $v = Zu$). $\square$

Since $t > 1$, the singular values of $\operatorname{diag}(t^{-1}, 1)$ are
$\{t^{-1}, 1\}$ with $\sigma_2 = t^{-1}$. Hence

$$\sigma_2(\tilde B) \ge \sigma_2(\tilde C)\,\sigma_2\!\big(\operatorname{diag}(t^{-1},1)\big) \ge \frac{1}{\sqrt{n-1}}\cdot\frac1t,$$

and squaring, with $t^{-2} = 1 - b^2 \ge 1 - \frac1n$:

$$\sigma_2^2(\tilde B) \ge \frac{1-b^2}{n-1} \ge \frac{1 - \frac1n}{n-1} = \frac{(n-1)/n}{n-1} = \frac1n .$$

This is exactly the claimed **constant self-reproduction**
$\sqrt{\tfrac{n-1}{n}}\cdot\tfrac{1}{\sqrt{n-1}} = \tfrac{1}{\sqrt n}$: the
loss factor $\frac1t = \sqrt{1-b^2} \ge \sqrt{\tfrac{n-1}{n}}$ from
renormalization exactly cancels the gain $\frac{1}{\sqrt{n-1}}$ vs
$\frac{1}{\sqrt n}$ from having one row fewer — this is precisely where the
Case-A threshold $\ell_1 \le \frac1n$ is used, and the chain is tight when
$b^2 = \frac1n$. Finally $\tilde B$ is a submatrix of $B$ (rows $i, j \in
\{2..n\}$), so by §4.1 the corresponding submatrix of $A$ has the same
$\sigma_2 \ge \frac{1}{\sqrt n}$. ✓

**Component verdict: CORRECT.** All computations of the paper reproduce
exactly; the only fill was the standard inequality of G2. Numeric confirmation:
check **N7** (1499 random low-leverage instances: rotation, deletion,
renormalization orthonormality to 1e-9, and the full chain $\sigma_2^2(\tilde
B) \ge \frac{1-b^2}{n-1} \ge \frac1n$; 2000 random checks of G2).

---

## 5. Case B, pair ⇒ bound — verdict CORRECT

**Hypothesis of the case:** $\ell_i > \frac1n$ for all $i$. Assume (proved in
§6) that there exist $i \ne j$ with

$$(r_i, r_j)^2 \le \Big(\ell_i - \frac1n\Big)\Big(\ell_j - \frac1n\Big). \tag{pair}$$

Let $\tilde A = \begin{pmatrix} r_i \\ r_j \end{pmatrix}$ and $\tilde G =
\tilde A \tilde A^T = \begin{pmatrix} \ell_i & (r_i,r_j) \\ (r_j,r_i) & \ell_j
\end{pmatrix}$, with eigenvalues $\tilde\lambda_1 \ge \tilde\lambda_2 \ge 0$
(real: symmetric; nonnegative: Gram). Its characteristic polynomial, in the
paper's convention,

$$P_{\tilde G}(\lambda) = \det(\tilde G - \lambda I) = (\ell_i - \lambda)(\ell_j - \lambda) - (r_i, r_j)^2 = (\tilde\lambda_1 - \lambda)(\tilde\lambda_2 - \lambda),$$

a monic quadratic in $\lambda$ (leading coefficient $(-\lambda)(-\lambda) =
\lambda^2$: sign convention consistent ✓). Evaluate at $\lambda = \frac1n$:

$$P_{\tilde G}\Big(\frac1n\Big) = \Big(\ell_i - \frac1n\Big)\Big(\ell_j - \frac1n\Big) - (r_i, r_j)^2 \ \ge\ 0 \quad\text{by (pair)},$$

hence $(\tilde\lambda_1 - \frac1n)(\tilde\lambda_2 - \frac1n) \ge 0$. Also
$\tilde\lambda_1 + \tilde\lambda_2 = \operatorname{Tr}\tilde G = \ell_i +
\ell_j > \frac2n$ (Case-B hypothesis, twice), so $\tilde\lambda_1 \ge
\frac{\ell_i + \ell_j}{2} > \frac1n$, i.e. $\tilde\lambda_1 - \frac1n > 0$
strictly. Dividing, $\tilde\lambda_2 - \frac1n \ge 0$. The singular values of
$\tilde A$ squared are the eigenvalues of $\tilde A \tilde A^T$, so

$$\sigma_2(\tilde A) = \sqrt{\tilde\lambda_2} \ge \frac{1}{\sqrt n}. \checkmark$$

The strictness bookkeeping is exactly right: the *strict* inequality
$\tilde\lambda_1 > \frac1n$ (from the Case-B hypothesis) is what licenses
dividing the product inequality; without it, $\tilde\lambda_1 =
\tilde\lambda_2 = $ anything with product condition could not be split.

**Gap G5 (typo).** The paper writes
"$Tr(\tilde G) = \|r_i\|^2 + \|r_k\|^2$" — the second index should be $j$
(and the bound quantifier line uses $k$ as a running index in the same
sentence). Purely typographical; the argument is unaffected.

**Component verdict: CORRECT** (one index typo). Numeric confirmation: check
**N6/N6b** (on 515 Case-B instances, the pair maximizing the (pair)-slack
always satisfies $\sigma_2^2 \ge \frac1n$ to 1e-9).

---

## 6. Pair existence — verdict CORRECT

Setting: Case B, so $A^TA = I$ gives $\sum x_i^2 = \sum y_i^2 = 1$, $\sum x_i
y_i = 0$, and $\sum_i \ell_i = \operatorname{Tr}(AA^T) = \operatorname{Tr}(A^TA)
= 2$. ✓ (paper's eq. (11)).

### 6.1 The algebraic identity and the squaring map

**Identity (paper's eq. (13)), verified by expansion:**

$$(x_ix_j + y_iy_j)^2 = \tfrac12(x_i^2+y_i^2)(x_j^2+y_j^2) + \tfrac12\big[(x_i^2-y_i^2)(x_j^2-y_j^2) + 4x_iy_ix_jy_j\big].$$

Expansion of the right side:
$\tfrac12[\,x_i^2x_j^2 + x_i^2y_j^2 + y_i^2x_j^2 + y_i^2y_j^2\,] +
\tfrac12[\,x_i^2x_j^2 - x_i^2y_j^2 - y_i^2x_j^2 + y_i^2y_j^2 +
4x_iy_ix_jy_j\,] = x_i^2x_j^2 + y_i^2y_j^2 + 2x_ix_jy_iy_j = (x_ix_j +
y_iy_j)^2$. ✓ (This is the real form of $\operatorname{Re}(u\bar v)^2 =
\tfrac12|u|^2|v|^2 + \tfrac12\operatorname{Re}(u^2\bar v^2)$ for $u = x_i +
\mathrm{i}y_i$, $v = x_j + \mathrm{i}y_j$ — the squaring map below is $u
\mapsto u^2$; this is the step that genuinely uses $k = 2$ and realness.)

**Squaring map:** $w_i = (x_i^2 - y_i^2,\; 2x_iy_i)$. Two properties, both
verified:

- $\|w_i\|^2 = (x_i^2-y_i^2)^2 + 4x_i^2y_i^2 = (x_i^2+y_i^2)^2 = \ell_i^2$,
  i.e. $\|w_i\| = \ell_i$. ✓
- $\sum_i w_i = (\sum x_i^2 - \sum y_i^2,\; 2\sum x_iy_i) = (1-1,\, 0) =
  (0,0)$. ✓

The identity becomes $(r_i,r_j)^2 = \tfrac12 \ell_i\ell_j + \tfrac12 (w_i,
w_j)$. ✓

### 6.2 Equivalent form of the pair inequality

Substituting and multiplying by 2, (pair) $\iff$

$$\ell_i\ell_j + (w_i,w_j) \le 2\Big(\ell_i - \frac1n\Big)\Big(\ell_j - \frac1n\Big) = 2\ell_i\ell_j - \frac2n(\ell_i + \ell_j) + \frac{2}{n^2}$$
$$\iff (w_i,w_j) \le \ell_i\ell_j - \frac2n(\ell_i+\ell_j) + \frac{2}{n^2} \iff (w_i,w_j) + \frac{2}{n^2} \le \ell_i\ell_j - \frac2n(\ell_i+\ell_j) + \frac{4}{n^2},$$

and the right side is exactly $\big(\ell_i - \frac2n\big)\big(\ell_j -
\frac2n\big)$ ✓ (check: cross terms $-\frac2n(\ell_i+\ell_j)$, constant
$\frac4{n^2}$). With $z_i = \ell_i - \frac2n$ (so $\sum z_i = 2 - n\cdot\frac2n
= 0$ ✓):

$$\text{(pair)} \iff (w_i, w_j) - z_iz_j + \frac{2}{n^2} \le 0. \tag{pair′}$$

All three displayed equivalences of the paper reproduce exactly. Numeric
confirmation: **N1** (2000 random tuples), **N2** (all sweep matrices).

### 6.3 Contradiction setup and the eigenvalue count

Assume (pair′) fails for **all** pairs $(i,j)$, including $i = j$ (this is the
correct negation of "∃ pair"; the case $i=j$ is dispatched in §6.6):

$$(w_i, w_j) - z_i z_j + \frac{2}{n^2} > 0 \quad \forall\, i, j. \tag{¬}$$

Let $G = WW^T - zz^T$ ($W \in \mathbb{R}^{n\times2}$ with rows $w_i$, $z =
(z_i)$), so $G_{ij} = (w_i,w_j) - z_iz_j$, symmetric.

**Eigenvalue count.** $\operatorname{rank}(WW^T) \le \operatorname{rank}(W)
\le 2$ and $WW^T \succeq 0$, so its eigenvalues sorted decreasingly satisfy
$\lambda_3(WW^T) = 0 \le 0$. ✓ The paper then uses: for symmetric $X$ and
symmetric PSD $Y$, $\lambda_p(X - Y) \le \lambda_p(X)$ for all $p$. Proof
(filling the one-line justification): by Courant–Fischer,
$$\lambda_p(X-Y) = \max_{\dim V = p}\ \min_{0 \ne v \in V} \frac{v^T(X-Y)v}{v^Tv} \le \max_{\dim V = p}\ \min_{0 \ne v \in V} \frac{v^TXv}{v^Tv} = \lambda_p(X),$$
using $v^T(X-Y)v \le v^TXv$ pointwise (as $v^TYv \ge 0$). $\square$ With $X =
WW^T$, $Y = zz^T \succeq 0$, $p = 3$: $\lambda_3(G) \le 0$, so $G$ has **at
most two positive eigenvalues**. ✓

### 6.4 Independent recomputation of the trace, and λ₁ ≥ 2/n

**Tr G, recomputed from scratch** (as required): using $\|w_i\|^2 = \ell_i^2$
(§6.1) and $\sum \ell_i = 2$:

$$\operatorname{Tr} G = \sum_{i=1}^n \big(\|w_i\|^2 - z_i^2\big)
= \sum_{i=1}^n \Big(\ell_i^2 - \Big(\ell_i - \frac2n\Big)^2\Big)
= \sum_{i=1}^n \Big(\frac4n \ell_i - \frac{4}{n^2}\Big)
= \frac4n \cdot 2 - n\cdot\frac{4}{n^2} = \frac8n - \frac4n = \boxed{\frac4n}.$$

This **matches the paper's value** $\operatorname{Tr} G = \frac4n$. (The
paper's intermediate line "$\frac4n(2) - \frac4n$" is the same computation;
correct.) Note this uses only $A^TA = I$ — it is unconditional, independent of
(¬) and of the Case-B hypothesis.

**λ₁ ≥ 2/n (fill of the one-line deduction).** Sort the eigenvalues
$\lambda_1 \ge \lambda_2 \ge \lambda_3 \ge \dots \ge \lambda_n$ of $G$. By
§6.3, $\lambda_p \le 0$ for $p \ge 3$, so

$$\lambda_1 + \lambda_2 = \operatorname{Tr} G - \sum_{p\ge3}\lambda_p \ge \operatorname{Tr} G = \frac4n,$$

and since $\lambda_1 \ge \lambda_2$: $2\lambda_1 \ge \lambda_1 + \lambda_2 \ge
\frac4n$, i.e. $\lambda_1 \ge \frac2n > 0$. ✓ (Also unconditional.) Numeric
confirmation: **N3** — on 4400 random matrices, $\operatorname{Tr}G = \frac4n$
to 1e-10, $\lambda_3(G) \le 10^{-9}$, and $\min$ of $\lambda_1\cdot\frac n2 =
1.0350 \ge 1$.

### 6.5 Perron–Frobenius contradiction (gap G3 filled)

Let $E$ be the all-ones $n\times n$ matrix and $M = G + \frac{2}{n^2}E$, i.e.
$M_{ij} = G_{ij} + \frac{2}{n^2}$. Under (¬), $M$ is **entrywise strictly
positive** (all entries, diagonal included). ✓

**$G\mathbf{1} = 0$:** $(G\mathbf1)_i = \sum_j (w_i, w_j) - z_i \sum_j z_j =
(w_i, \sum_j w_j) - z_i\cdot0 = (w_i, 0) = 0$, by §6.1–6.2. ✓ Hence

$$M\mathbf1 = G\mathbf1 + \frac{2}{n^2} E\mathbf1 = 0 + \frac{2}{n^2}\, n\,\mathbf1 = \frac2n \mathbf1,$$

so $\frac2n$ is an eigenvalue of $M$ with strictly positive eigenvector
$\mathbf1$. ✓

**Gap G3 (fill).** The paper invokes Perron–Frobenius to say the eigenspace of
$\frac2n$ is one-dimensional and every *other* eigenvalue $\lambda$ of $M$ has
$|\lambda| < \frac2n$. Two sub-steps deserve justification:

*(G3a) Why $\frac2n$ is the Perron root of $M$.* Perron's theorem for
entrywise-positive matrices states: $\rho = \rho(M) > 0$ is a simple
eigenvalue with an entrywise-positive eigenvector, and every other eigenvalue
$\lambda$ satisfies $|\lambda| < \rho$. It remains to see that our eigenvalue
$\frac2n$ — which comes with a positive eigenvector — *is* $\rho$. Standard
lemma: a positive matrix has exactly one eigenvalue admitting a positive
eigenvector, namely $\rho$. Proof: let $u > 0$ be the Perron eigenvector of
$M^T$ (also entrywise positive, $M^Tu = \rho u$; here $M^T = M$, but the
argument is general). If $Mv = \mu v$ with $v > 0$, then $\rho\, u^Tv = (M^T
u)^T v = u^T M v = \mu\, u^Tv$, and $u^Tv > 0$ (both positive), so $\mu =
\rho$. $\square$ Hence $\rho(M) = \frac2n$, its eigenspace is
$\operatorname{span}(\mathbf1)$ (simplicity), and every eigenvalue $\lambda \ne
\frac2n$ has $|\lambda| < \frac2n$.

*(G3b) Why $\lambda_1$ contradicts this.* Since $\lambda_1 \ge \frac2n > 0$
(§6.4) and $G\mathbf1 = 0$, $\lambda_1 \ne 0$; let $v_1$ be an eigenvector of
$G$ for $\lambda_1$. Eigenvectors of a symmetric matrix for distinct
eigenvalues are orthogonal, and $\mathbf1$ is an eigenvector for $0 \ne
\lambda_1$, so $v_1 \perp \mathbf1$ automatically (the paper says this via the
spectral theorem; correct). Then $Ev_1 = \mathbf1(\mathbf1^Tv_1) = 0$, so $Mv_1
= Gv_1 = \lambda_1 v_1$: $\lambda_1$ is an eigenvalue of $M$ with eigenvector
$v_1 \perp \mathbf1$. Now the case split the paper leaves implicit: if
$\lambda_1 = \frac2n$, then $v_1$ lies in the eigenspace of $\frac2n$, which is
$\operatorname{span}(\mathbf1)$ — impossible for a nonzero vector orthogonal to
$\mathbf1$. Hence $\lambda_1 \ne \frac2n$, so by (G3a) $|\lambda_1| < \frac2n$,
contradicting $\lambda_1 \ge \frac2n$. $\square$

Therefore (¬) is false: **some pair $(i,j)$ satisfies (pair′)**, i.e.
$\min_{i,j} M_{ij} \le 0$. Numeric confirmation: **N4** (unconditionally, on
all 4400 sweep matrices, $\min_{ij} [G_{ij} + \frac{2}{n^2}] \le 0$).

### 6.6 Ruling out i = j (gap G4 — verified, correctly handled)

If the pair found has $i = j$, then (pair′) reads $\|w_i\|^2 - z_i^2 +
\frac{2}{n^2} \le 0$. Substituting $\|w_i\|^2 = \ell_i^2$, $z_i = \ell_i -
\frac2n$:

$$\ell_i^2 - \ell_i^2 + \frac4n\ell_i - \frac{4}{n^2} + \frac{2}{n^2} \le 0 \iff \frac4n \ell_i \le \frac{2}{n^2} \iff \ell_i \le \frac{1}{2n},$$

which contradicts the Case-B hypothesis $\ell_i > \frac1n$. This matches the
paper's "after simplification we obtain $\|r_i\|^2 \le \frac{1}{2n}$" ✓ —
the omitted simplification is exactly the display above. Consequently every
pair satisfying (pair′) has $i \ne j$, and at least one such pair exists.
Combined with §5, Case B is complete. The logical order is sound: existence
first (§6.5), then the $i = j$ possibility excluded for *any* satisfying pair
(§6.6), so an $i \ne j$ pair survives.

**Component verdict: CORRECT.** The two glosses (G3a, G3b) are filled; the
trace value $\frac4n$ is independently confirmed both symbolically and
numerically. Numeric confirmation of the end-to-end claim: **N5/N5b** — on all
515 Case-B instances encountered, an $i \ne j$ pair with $(r_i,r_j)^2 \le
(\ell_i - \frac1n)(\ell_j - \frac1n)$ exists (slack $\ge -10^{-9}$).

---

## 7. Assembly of the induction — verdict CORRECT

Fix $n \ge 4$ and assume the theorem for $n - 1$ ($\ge 3$, covered by the base
§3.1 when $n - 1 = 3$). Cases A and B are exhaustive and each produces a
$2\times2$ submatrix of $A$ with $\sigma_2 \ge \frac{1}{\sqrt n}$ (§4 uses the
hypothesis at $n-1$; §§5–6 use nothing beyond $A^TA = I$ and the case
hypothesis). With the base $n = 3$, strong induction gives the theorem for all
$n \ge 3$. ✓ No circularity: Case B never invokes the induction hypothesis,
and Case A invokes it only at size $n-1$.

**Sharpness (paper §3, checked though not required for the theorem):** the
family with $n-2$ rows $X = (a, 0)$, plus $Y = (b, c)$, $Z = (b, -c)$, $a =
\sqrt{\frac{n-1}{n(n-2)}}$, $b = \frac{1}{\sqrt{2n}}$, $c = \frac{1}{\sqrt2}$
(the MEMORY.md equality family). Orthonormality: $(n-2)a^2 + 2b^2 =
\frac{n-1}{n} + \frac1n = 1$; $2c^2 = 1$; cross terms $(n-2)\cdot a\cdot0 + bc
- bc = 0$. ✓ Numerically (check **N9**, $n = 4..100$): $f(A) =
\max_{i<j}\sigma_2 = \frac{1}{\sqrt n}$ to 1e-12, attained on $(X,Y), (X,Z),
(Y,Z)$ as the paper states, all rows in Case B. So $\frac{1}{\sqrt n}$ cannot
be improved — consistent with MEMORY.md (Nesterenko arXiv:2511.02387:
sharp for all $n > k > 0$).

---

## 8. Gap ledger (complete list of everything the paper omits or gets wrong)

| # | Location in paper | Nature | Severity | Fix |
|---|---|---|---|---|
| G1 | §2 opening: "$n=3$ is trivial", $n=4$ cited to Nesterenko | Missing base-case proof; external dependency (whose write-up has a known typo per MEMORY.md) | Gap (fixable) | §3.1 direct proof of $n=3$; §3.2 shows $n=4$ follows from the step, dependency removable |
| G2 | Case A, eq. (10): $\sigma_2(YZ) \ge \sigma_2(Y)\sigma_2(Z)$ | Standard fact used without proof | Cosmetic | Proved in §4.4 |
| G3 | Case B Perron step | (a) eigenvalue-with-positive-eigenvector ⇒ Perron root not justified; (b) case $\lambda_1 = \frac2n$ vs "any other eigenvalue" not addressed | Gap (fixable) | Both filled in §6.5 |
| G4 | Case B, exclusion of $i=j$: "after simplification" | Omitted computation | Cosmetic | Reproduced in §6.6; paper's conclusion $\ell_i \le \frac{1}{2n}$ is exactly right |
| G5 | Case B: "$Tr(\tilde G) = \|r_i\|^2 + \|r_k\|^2$" | Index typo ($k$ for $j$) | Typo | Noted in §5 |

Nothing else was found: all displayed identities, all constants, all
inequality directions, and the strictness bookkeeping (where $>$ vs $\ge$
matters: $t^2 > 1$ in §4.3, $\tilde\lambda_1 > \frac1n$ in §5, strict
positivity of $M$ in §6.5) check out. Degenerate cases examined: $b = 0$
(handled by the paper explicitly, §4.1), repeated eigenvalues of $G$
(harmless: G3b uses only eigen*space* dimension), zero rows (a zero row has
$\ell = 0 \le \frac1n$ → Case A with $b = 0$), $n = 3$ in Case B (argument
runs verbatim; $G$ is $3\times3$ and "at most 2 positive eigenvalues" is still
a real constraint by Tr).

---

## 9. Numeric sanity block

Script: `numerics/study-sp-verification/check_sp.py`; run
`.venv/bin/python numerics/study-sp-verification/check_sp.py`; master seed
`20260717`, per-trial seeds `default_rng([20260717, <coords>])` as printed in
the script; exit code 0. Full output of the run of 2026-07-17:

```
[N1] identity (r_i,r_j)^2 = l_i l_j/2 + (w_i,w_j)/2 : OK (2000 random tuples)
[N2-N4,N8] sweep: n in [3, 4, 5, ..., 20, 30, 40, 64, 100], 200 trials each (4400 matrices) : OK
[N3] min over trials of lambda_1(G)/(2/n) = 1.034965  (claim: >= 1)
[N5,N6] Case-B instances encountered: 178 : OK
[N8] min over trials of f(A)*sqrt(n) = 1.049482  (claim: >= 1)
[N5b,N6b] perturbed-equality-family Case-B instances: 337 : OK
[N7] Case-A chain verified on 1499 low-leverage instances; submultiplicativity on 2000 random 2x2 pairs : OK
[N9] equality family n=4..100: A^T A=I, f=1/sqrt(n) exactly, attained on (X,Y),(X,Z),(Y,Z), all rows in Case B : OK

RESULT: ALL CHECKS PASSED
```

Mapping to claims: N1 → §6.1 identity; N2 → squaring-map properties §6.1–6.2;
N3 → $\operatorname{Tr}G = \frac4n$, $G\mathbf1 = 0$, $\lambda_3(G) \le 0$,
$\lambda_1(G) \ge \frac2n$ (§§6.3–6.4, unconditional); N4 → pair existence
(any $(i,j)$), §6.5; N5/N5b → **Case-B pair existence with $i \ne j$**
(contract-required test, 515 instances); N6/N6b → pair ⇒ bound §5; N7 → Case-A
chain §4 incl. $t^2 = \frac{1}{1-b^2}$ renormalization and G2; N8 → **the
k=2 GTZ bound itself** (contract-required test, 4400 random instances, min of
$f\sqrt n = 1.0495 \ge 1$); N9 → **equality family reproduces $f = 1/\sqrt n$**
(contract-required test, $n = 4..100$, to 1e-12).

No numeric violation of any claim was observed.

---

## 10. Status lines

- Theorem (GTZ, $k=2$, all $n > 2$): **PROVED** — by arXiv:2604.05944 with the
  gaps G1–G5 filled above (all fills are contained in this report; no step of
  the paper remains unverified).
- §3.1 base case $n = 3$: PROVED (this report).
- §4 Case A: PROVED (paper + fill G2).
- §5 pair ⇒ bound: PROVED (paper; typo G5 noted).
- §6 pair existence: PROVED (paper + fills G3a, G3b; trace $\frac4n$
  independently recomputed and confirmed).
- §7 sharpness of $\frac{1}{\sqrt n}$: NUMERICALLY SUPPORTED here (N9, exact to
  1e-12); a symbolic proof is in the paper's §3 ("direct evaluation") and in
  arXiv:2511.02387 — not required for, and not used by, the theorem.

**Consequence for the project:** per AGENTS.md rule 4, Sengupta–Pautov
2604.05944 may now be trusted to the extent verified here — namely its main
theorem in full. VAL-DUAL-002 and any duality-based transfer $\mathrm{GTZ}(n,2)
\Leftrightarrow \mathrm{GTZ}(n, n-2)$ can rely on the $k = 2$ case without
reservation.
