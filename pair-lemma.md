# The General-k Pair Lemma

**Deliverable:** standalone, referee-grade proof of the General-k Pair Lemma,
promoted from attack round 1 (`attacks/round-1-squaring-k3.md`).
**Overall status: PROVED** (unconditional; every statement below carries its
own status line).
Numerics: `numerics/study-pair-lemma/` (deterministic, seed 20260717; §7).

---

## 0. Scope — read this first (overclaim guard)

> **This lemma produces a good PAIR of rows — a 2×2 Gram submatrix bound —
> and nothing more.  It does NOT produce a good k×k submatrix, and therefore
> it does NOT prove the GTZ hypothesis for any k ≥ 3.**
>
> For k = 3 the lemma yields (under a solidity hypothesis, Corollary 5.2 and
> Remark 5.3) two rows whose 2×2 Gram matrix has λ_min ≥ 6/(5n) > 1/n.  GTZ
> requires a 3×3 (in general k×k) Gram submatrix with λ_min ≥ 1/n.  The
> **pair → k-subset bridge is explicitly OUT of scope of this document and
> remains open** (conjecture C_ext in the round-1 ledger is numerically
> supported but unproved).  Only at k = 2, where "pair" = "k-subset", does
> the lemma close the corresponding case of GTZ — and there it reduces
> exactly to the Sengupta–Pautov argument (§4).

---

## 1. Setting and statement

Throughout, $n \ge k \ge 2$ are integers and
$$\mathrm{St}(n,k) = \{A \in \mathbb R^{n\times k} : A^TA = I_k\}$$
(note $A^TA = I_k$ forces $n \ge k$).  For $A \in \mathrm{St}(n,k)$ write
$r_1,\dots,r_n \in \mathbb R^k$ for its rows,
$$\ell_i = \|r_i\|^2, \qquad c_{ij} = \langle r_i, r_j\rangle .$$
The Stiefel constraint is $\sum_{i=1}^n r_ir_i^T = A^TA = I_k$; taking traces,
$\sum_i \ell_i = k$.  Define the constants
$$s = \frac{k+2}{k^2}, \qquad t = \frac{2k}{(k+2)\,n}, \qquad c = \frac tn
  = \frac{2k}{(k+2)n^2},$$
and, for $i \ne j$, the **pair margin**
$$\mathrm{marg}(i,j) = s\,(\ell_i - t)(\ell_j - t) - c_{ij}^2 .$$
A pair $\{i,j\}$, $i \ne j$, is **good** when $\mathrm{marg}(i,j) \ge 0$.

> **Theorem 1 (General-k Pair Lemma).**
> *For every $n \ge k \ge 2$ and every $A \in \mathrm{St}(n,k)$ there exist
> indices $i \ne j$ with*
> $$c_{ij}^2 \;\le\; \frac{k+2}{k^2}\,(\ell_i - t)(\ell_j - t),
> \qquad t = \frac{2k}{(k+2)\,n}.$$
> *That is, a good pair always exists.*
>
> **Status: PROVED** (§§2–3; unconditional — no assumption on the leverages
> is made anywhere).

At $k = 2$: $s = 1$, $t = 1/n$, and Theorem 1 is *verbatim* the
Sengupta–Pautov pair inequality (§4).  At $k = 3$: $s = 5/9$,
$t = 6/(5n)$, and a good pair with both leverages $\ge t$ has 2×2 Gram
$\lambda_{\min} \ge 6/(5n) > 1/n$ (§5).  The constant $t$ is sharp at
$(n,k) = (6,3)$ (§6).

---

## 2. The traceless-symmetric embedding

Let $\mathrm{Sym}(\mathbb R^k)$ be the space of symmetric $k\times k$
matrices with the Frobenius inner product
$\langle X, Y\rangle = \operatorname{Tr}(XY)$, and
$$\mathrm{Sym}_0(\mathbb R^k) = \{X \in \mathrm{Sym}(\mathbb R^k):
\operatorname{Tr}X = 0\}.$$
Since $\dim \mathrm{Sym}(\mathbb R^k) = k(k+1)/2$ and the trace is a nonzero
linear functional, its kernel has
$$\dim \mathrm{Sym}_0(\mathbb R^k) = \frac{k(k+1)}2 - 1
 = \frac{(k+2)(k-1)}2 =: m .$$
(For $k \ge 2$, $m \ge 2$.)

Define, for $i = 1,\dots,n$,
$$w_i = r_ir_i^T - \frac{\ell_i}{k}\,I_k \;\in\; \mathrm{Sym}_0(\mathbb R^k),
\qquad
z_i = \frac{\sqrt2}{k}\Bigl(\ell_i - \frac kn\Bigr) \in \mathbb R,$$
($w_i$ is the traceless part of $r_ir_i^T$; the trace of $r_i r_i^T$ is
$\ell_i$, so $\operatorname{Tr} w_i = 0$ indeed).  Fix any orthonormal basis
of $\mathrm{Sym}_0(\mathbb R^k)$ and let $W \in \mathbb R^{n\times m}$ be the
matrix whose $i$-th row holds the coordinates of $w_i$; then
$(WW^T)_{ij} = \langle w_i, w_j\rangle$ regardless of the basis chosen.
Let $z = (z_1,\dots,z_n)^T$ and
$$G \;=\; WW^T - zz^T \;\in\; \mathbb R^{n\times n}
\qquad (\text{symmetric}).$$

We now prove the seven framework properties (P1)–(P7).  Throughout,
$\mathbf 1 \in \mathbb R^n$ is the all-ones vector and $Q = \sum_i \ell_i^2$.

**(P1) $\sum_i w_i = 0$ (equivalently $W^T\mathbf 1 = 0$) and
$z^T\mathbf 1 = 0$.**
*Proof.* $\sum_i w_i = \sum_i r_ir_i^T - \frac{\sum_i \ell_i}{k} I_k
= I_k - \frac kk I_k = 0$, using the Stiefel constraint and
$\sum_i \ell_i = k$.  And
$z^T\mathbf 1 = \frac{\sqrt2}{k}\bigl(\sum_i \ell_i - n\cdot\frac kn\bigr)
= \frac{\sqrt2}{k}(k - k) = 0$. $\square$  **Status: PROVED.**

**(P2) $\langle w_i, w_j\rangle = c_{ij}^2 - \dfrac{\ell_i\ell_j}{k}$ for all
$i, j$ (including $i = j$: $\|w_i\|^2 = \frac{k-1}{k}\ell_i^2$).**
*Proof.* Using $\operatorname{Tr}(r_ir_i^T r_jr_j^T)
= r_j^T r_i\, r_i^T r_j = c_{ij}^2$ and
$\operatorname{Tr}(r_ir_i^T) = \ell_i$:
$$\langle w_i, w_j\rangle
= \operatorname{Tr}\Bigl[\Bigl(r_ir_i^T - \tfrac{\ell_i}k I\Bigr)
                         \Bigl(r_jr_j^T - \tfrac{\ell_j}k I\Bigr)\Bigr]
= c_{ij}^2 - \frac{\ell_i\ell_j}{k} - \frac{\ell_i\ell_j}{k}
  + \frac{\ell_i\ell_j}{k^2}\cdot k
= c_{ij}^2 - \frac{\ell_i\ell_j}{k}.$$
For $i = j$, $c_{ii} = \ell_i$ gives $\|w_i\|^2 = \ell_i^2 - \ell_i^2/k$.
$\square$  **Status: PROVED.**

**(P3) For $i \ne j$:
$G_{ij} = c_{ij}^2 - s(\ell_i - t)(\ell_j - t) - c$; consequently**
$$c_{ij}^2 \le s(\ell_i - t)(\ell_j - t) \iff G_{ij} \le -c .$$
*Proof.* By (P2) and the definition of $z$,
$$G_{ij} = c_{ij}^2 - \frac{\ell_i\ell_j}{k}
 - \frac{2}{k^2}\Bigl(\ell_i - \frac kn\Bigr)\Bigl(\ell_j - \frac kn\Bigr).$$
Expand both this and $s(\ell_i-t)(\ell_j-t) + c$ as bilinear polynomials in
$(\ell_i, \ell_j)$ and compare coefficients:
* $\ell_i\ell_j$: $\dfrac1k + \dfrac2{k^2} = \dfrac{k+2}{k^2} = s$;
* $\ell_i$ and $\ell_j$: $\dfrac2{k^2}\cdot\dfrac kn = \dfrac2{kn}$, while
  $s\,t = \dfrac{k+2}{k^2}\cdot\dfrac{2k}{(k+2)n} = \dfrac2{kn}$;
* constant: $\dfrac{2}{k^2}\cdot\dfrac{k^2}{n^2} = \dfrac2{n^2}$, while
  $s\,t^2 + c = st\cdot t + \dfrac tn
  = \dfrac{2}{kn}\cdot\dfrac{2k}{(k+2)n} + \dfrac{2k}{(k+2)n^2}
  = \dfrac{4 + 2k}{(k+2)n^2} = \dfrac2{n^2}$.
Hence $\frac{\ell_i\ell_j}{k} + \frac2{k^2}(\ell_i-\frac kn)(\ell_j-\frac kn)
= s(\ell_i-t)(\ell_j-t) + c$ identically, which is the claim.  The
equivalence follows by rearranging. $\square$  **Status: PROVED.**
(Symbolic verification of all three coefficient identities: §7, check C7.)

**(P4) $G\mathbf 1 = 0$.**
*Proof.* $G\mathbf 1 = W(W^T\mathbf 1) - z(z^T\mathbf 1) = 0$ by (P1).
$\square$  **Status: PROVED.**

**(P5) $\operatorname{Tr} G = \dfrac{(k-2)(k+1)}{k^2}\,Q + \dfrac2n
\;\ge\; \dfrac{k(k-1)}{n} = m\,t > 0$.**
*Proof.*  By (P2), $\operatorname{Tr}(WW^T) = \sum_i \|w_i\|^2
= \frac{k-1}{k} Q$, and
$$\sum_i z_i^2 = \frac{2}{k^2}\sum_i \Bigl(\ell_i - \frac kn\Bigr)^2
= \frac{2}{k^2}\Bigl(Q - \frac{2k}{n}\sum_i\ell_i + n\cdot\frac{k^2}{n^2}\Bigr)
= \frac{2}{k^2}\Bigl(Q - \frac{k^2}{n}\Bigr).$$
Hence
$$\operatorname{Tr}G = \frac{k-1}{k}Q - \frac{2}{k^2}\Bigl(Q - \frac{k^2}n\Bigr)
= \frac{k^2 - k - 2}{k^2}\,Q + \frac2n
= \frac{(k-2)(k+1)}{k^2}\,Q + \frac2n .$$
By Cauchy–Schwarz, $k^2 = (\sum_i\ell_i)^2 \le n\sum_i\ell_i^2 = nQ$, so
$Q \ge k^2/n$; and $(k-2)(k+1) \ge 0$ for $k \ge 2$.  Therefore
$$\operatorname{Tr}G \;\ge\; \frac{(k-2)(k+1)}{k^2}\cdot\frac{k^2}{n}
 + \frac2n = \frac{(k-2)(k+1) + 2}{n} = \frac{k^2-k}{n} = \frac{k(k-1)}n,$$
and $m\,t = \frac{(k+2)(k-1)}{2}\cdot\frac{2k}{(k+2)n} = \frac{k(k-1)}n$.
Finally $k(k-1)/n > 0$ since $k \ge 2$. $\square$  **Status: PROVED.**
(At $k = 2$ the $Q$-coefficient vanishes: $\operatorname{Tr}G = 2/n$
always — cf. §4.)

**(P6) $G$ has at most $m$ positive eigenvalues.**
*Proof.*  If $n \le m$ this is trivial.  Assume $n > m$.  Since
$zz^T \succeq 0$ we have $G = WW^T - zz^T \preceq WW^T$.  The rows of $W$
lie in $\mathbb R^m$, so $\operatorname{rank}(WW^T) = \operatorname{rank}(W)
\le m$; as $WW^T \succeq 0$, its eigenvalues in decreasing order satisfy
$\lambda_{m+1}(WW^T) = 0$.  By the Weyl monotonicity principle (for symmetric
$X \preceq Y$, Courant–Fischer gives $\lambda_j(X) \le \lambda_j(Y)$ for
every $j$),
$$\lambda_{m+1}(G) \le \lambda_{m+1}(WW^T) = 0,$$
so at most the top $m$ eigenvalues of $G$ can be positive. $\square$
**Status: PROVED.**
*(Remark: since $G\mathbf 1 = 0$ by (P4), zero is an eigenvalue and the cap
can be sharpened to $\min(m, n-1)$; this refinement is not used in this
document.)*

**(P7) $\lambda_1(G) \ge \operatorname{Tr}G / m \ge t$.**
*Proof.*  Let $\lambda_1 \ge \dots \ge \lambda_n$ be the eigenvalues of $G$
and $p \in \{0,\dots,m\}$ the number of positive ones ((P6)).  If $p = 0$
then $\operatorname{Tr}G = \sum_i \lambda_i \le 0$, contradicting (P5); so
$p \ge 1$.  Dropping the nonpositive eigenvalues can only increase the sum:
$$\operatorname{Tr}G \le \sum_{i \le p}\lambda_i \le p\,\lambda_1
\le m\,\lambda_1,$$
where the last step uses $\lambda_1 > 0$ and $p \le m$.  Hence
$\lambda_1 \ge \operatorname{Tr}G/m \ge (m t)/m = t$ by (P5). $\square$
**Status: PROVED.**

---

## 3. The Perron–Frobenius step and the proof of Theorem 1

**Lemma 3.1 (Perron package for symmetric nonnegative matrices).**
*Let $M \in \mathbb R^{n\times n}$, $n \ge 2$, be symmetric, entrywise
nonnegative, and irreducible.  Then:*

* *(a) the spectral radius $\rho(M)$ is an algebraically simple eigenvalue
  of $M$ possessing an entrywise positive eigenvector $p$;*
* *(b) $\lambda_{\max}(M) = \rho(M)$;*
* *(c) if $Mx = \mu x$ for some entrywise **positive** $x$, then
  $\mu = \rho(M)$ and the full eigenspace of $\rho(M)$ is
  $\mathbb R x$ — i.e. a nonnegative irreducible symmetric matrix with a
  positive eigenvector has that eigenvalue as its simple spectral radius;*
* *(d) if $v^TMv = \lambda_{\max}(M)$ for some unit vector $v$, then
  $Mv = \lambda_{\max}(M)\,v$.*

*Proof.* (a) is the Perron–Frobenius theorem for irreducible nonnegative
matrices (Horn–Johnson, *Matrix Analysis*, 2nd ed., Theorem 8.4.4).

(b) $M$ is symmetric, so all eigenvalues are real and
$\rho(M) = \max_i|\lambda_i|$; in particular $\lambda_{\max}(M) \le \rho(M)$.
By (a), $\rho(M)$ is itself an eigenvalue, so $\rho(M) \le \lambda_{\max}(M)$.

(c) Let $p > 0$ be the Perron vector from (a).  If $\mu \ne \rho(M)$, then
$x$ and $p$ are eigenvectors of the symmetric $M$ for distinct eigenvalues,
hence orthogonal — impossible, since $\langle x, p\rangle > 0$ for two
entrywise positive vectors.  So $\mu = \rho(M)$; by algebraic (hence
geometric) simplicity the eigenspace is one-dimensional and contains $x$,
so it equals $\mathbb Rx$.

(d) Write the spectral decomposition $M = \sum_q \lambda_q u_qu_q^T$
(orthonormal eigenvectors $u_q$) and $v = \sum_q a_q u_q$,
$\sum_q a_q^2 = 1$.  Then $v^TMv = \sum_q \lambda_q a_q^2
= \lambda_{\max}(M)$ forces $a_q = 0$ for every $q$ with
$\lambda_q < \lambda_{\max}(M)$; hence $v$ lies in the
$\lambda_{\max}$-eigenspace. $\square$  **Status: PROVED** (external input:
Horn–Johnson Thm 8.4.4 only).

**Proof of Theorem 1.**
Suppose, for contradiction, that **no** good pair exists; by (P3) this says
$$G_{ij} > -c \qquad \text{for all } i \ne j. \tag{¬}$$
Let $E = \mathbf 1\mathbf 1^T$,
$$\tau = \max\Bigl(0,\; -\min_i\,(G_{ii} + c)\Bigr) \;\ge\; 0,
\qquad M = G + cE + \tau I .$$
Then:

1. *$M$ is symmetric and entrywise nonnegative:* off the diagonal
   $M_{ij} = G_{ij} + c > 0$ strictly by (¬); on the diagonal
   $M_{ii} = G_{ii} + c + \tau \ge 0$ by the choice of $\tau$.
2. *$M$ is irreducible:* all off-diagonal entries are strictly positive and
   $n \ge k \ge 2$, so the associated digraph is complete, hence strongly
   connected.
3. *$\mathbf 1$ is an eigenvector:* by (P4),
   $$M\mathbf 1 = G\mathbf 1 + c\,E\mathbf 1 + \tau\mathbf 1
   = (cn + \tau)\,\mathbf 1 = (t + \tau)\,\mathbf 1,$$
   since $cn = t$ by the definition $c = t/n$.

$\mathbf 1$ is entrywise positive, so Lemma 3.1(c) applies:
$$\rho(M) = t + \tau \text{ is a simple eigenvalue of } M
\text{ with eigenspace } \mathbb R\mathbf 1,$$
and by Lemma 3.1(b), $\lambda_{\max}(M) = t + \tau$.

Now let $v$ be a **unit** eigenvector of $G$ for its largest eigenvalue
$\lambda_1 = \lambda_1(G)$; by (P7), $\lambda_1 \ge t$.  Compute the Rayleigh
quotient of $M$ at $v$:
$$v^TMv = v^TGv + c\,(\mathbf 1^Tv)^2 + \tau\|v\|^2
       = \lambda_1 + c\,(\mathbf 1^Tv)^2 + \tau
\;\le\; \lambda_{\max}(M) = t + \tau .$$
Subtracting $\tau$:
$$\lambda_1 + c\,(\mathbf 1^Tv)^2 \le t .$$
Since $\lambda_1 \ge t$ and $c > 0$, **both** of the following are forced:
$$\mathbf 1^Tv = 0 \qquad\text{and}\qquad \lambda_1 = t .$$

*Equality case $\lambda_1 = t$, handled explicitly.*  With these values,
$v^TMv = t + \tau = \lambda_{\max}(M)$ exactly, so by Lemma 3.1(d) $v$ is an
eigenvector of $M$ for $\lambda_{\max}(M) = \rho(M)$.  But the
$\rho(M)$-eigenspace is $\mathbb R\mathbf 1$ (simplicity, above), so
$v = \pm\mathbf 1/\sqrt n$ and $\mathbf 1^Tv = \pm\sqrt n \ne 0$ —
contradicting $\mathbf 1^Tv = 0$.

Thus (¬) is impossible, and a good pair exists. $\blacksquare$
**Status: PROVED.**

**Remark 3.2 (what was *not* assumed).**  No lower bound on the leverages is
used anywhere: Theorem 1 holds for every Stiefel matrix, including ones with
zero rows (for which $w_i = 0$ and the framework degrades gracefully).  In
particular the lemma is *unconditional*, unlike the Case-B setting of
arXiv:2604.05944, which assumes $\ell_i > 1/n$ for all $i$.  Also note that
the negation (¬) only concerns $i \ne j$; the diagonal of $M$ is controlled
by the shift $\tau$, so no separate "diagonal pair" case arises (compare
`proofs/sp-verification.md` §6.6, where the S–P normalization requires an
explicit $i = j$ dispatch).

---

## 4. Exact reduction to Sengupta–Pautov at k = 2

**Proposition 4.1.**  *At $k = 2$, Theorem 1 is verbatim the
Sengupta–Pautov pair inequality: $s = 1$, $t = 1/n$, and the statement
reads: every $A \in \mathrm{St}(n,2)$ has a pair $i \ne j$ with*
$$c_{ij}^2 \;\le\; \Bigl(\ell_i - \frac1n\Bigr)\Bigl(\ell_j - \frac1n\Bigr).$$
*Moreover the proof apparatus of §§2–3 coincides with the S–P apparatus
(arXiv:2604.05944, Case B; verified in `proofs/sp-verification.md` §6) up to
the exact global scaling $G = \tfrac12\,G^{\mathrm{SP}}$.*

*Proof.*  Substituting $k = 2$:
$$s = \frac{2+2}{2^2} = 1, \qquad t = \frac{2\cdot2}{(2+2)n} = \frac1n,
\qquad c = \frac{t}{n} = \frac1{n^2}, \qquad m = \frac{(2+2)(2-1)}2 = 2,$$
so the statement of Theorem 1 becomes exactly the display above — the pair
inequality of arXiv:2604.05944 (their Case B; statement (pair) in
`proofs/sp-verification.md` §5).

For the apparatus, write $r_i = (x_i, y_i)$ and recall the S–P squaring map
$w^{\mathrm{SP}}_i = (x_i^2 - y_i^2,\, 2x_iy_i)$,
$z^{\mathrm{SP}}_i = \ell_i - \frac2n$,
$G^{\mathrm{SP}} = W^{\mathrm{SP}}(W^{\mathrm{SP}})^T
- z^{\mathrm{SP}}(z^{\mathrm{SP}})^T$.  Choose the orthonormal basis of
$\mathrm{Sym}_0(\mathbb R^2)$
$$B_1 = \tfrac1{\sqrt2}\begin{pmatrix}1&0\\0&-1\end{pmatrix}, \qquad
B_2 = \tfrac1{\sqrt2}\begin{pmatrix}0&1\\1&0\end{pmatrix}.$$
Our embedding is
$$w_i = r_ir_i^T - \frac{\ell_i}2 I
= \begin{pmatrix}\frac{x_i^2-y_i^2}2 & x_iy_i\\
                 x_iy_i & \frac{y_i^2-x_i^2}2\end{pmatrix},$$
with coordinates
$\langle w_i, B_1\rangle = \tfrac1{\sqrt2}(x_i^2-y_i^2)$ and
$\langle w_i, B_2\rangle = \tfrac1{\sqrt2}\cdot 2x_iy_i$; that is,
$W = \tfrac1{\sqrt2}\,W^{\mathrm{SP}}$.  Likewise
$z_i = \tfrac{\sqrt2}{2}(\ell_i - \tfrac2n) = \tfrac1{\sqrt2}\,
z^{\mathrm{SP}}_i$.  Hence
$$G = WW^T - zz^T = \tfrac12\Bigl(W^{\mathrm{SP}}(W^{\mathrm{SP}})^T
- z^{\mathrm{SP}}(z^{\mathrm{SP}})^T\Bigr) = \tfrac12\,G^{\mathrm{SP}} .$$
Consequently: $\operatorname{Tr}G = \tfrac2n$ corresponds to the S–P trace
$\operatorname{Tr}G^{\mathrm{SP}} = \tfrac4n$ (verified independently in
`proofs/sp-verification.md` §6.4); the good-pair criterion
$G_{ij} \le -c = -\tfrac1{n^2}$ (P3) corresponds to
$G^{\mathrm{SP}}_{ij} \le -\tfrac2{n^2}$, which is exactly the S–P
contradiction hypothesis (their matrix $G^{\mathrm{SP}} + \tfrac2{n^2}E$,
`proofs/sp-verification.md` §6.4–6.5); and the positive-eigenvalue cap is
$m = 2 = k$.  So at $k = 2$ the present proof *is* the S–P proof up to the
factor $\tfrac12$ and the isometric identification
$\mathrm{Sym}_0(\mathbb R^2) \cong \mathbb R^2$. $\square$
**Status: PROVED.**  (Numeric check of $2G = G^{\mathrm{SP}}$ at
$2.2\cdot10^{-16}$: §7, check C5.)

**Remark 4.2.**  S–P state and use the pair inequality only inside their
Case B (all $\ell_i > 1/n$); as recorded in `proofs/sp-verification.md`
§6.4, their argument never uses that hypothesis, and Theorem 1 confirms the
unconditional statement for all $k \ge 2$.

---

## 5. The k = 3 consequence: a spectrally good pair (conditional on solidity)

**Corollary 5.1 (solid/thin dichotomy).**  *If $\{i,j\}$ is a good pair,
then $(\ell_i - t)(\ell_j - t) \ge 0$: either both $\ell_i, \ell_j \ge t$
(a **solid** pair) or both $\ell_i, \ell_j \le t$ (a **thin** pair).*

*Proof.*  $0 \le c_{ij}^2 \le s(\ell_i - t)(\ell_j - t)$ and $s > 0$, so the
product $(\ell_i - t)(\ell_j - t)$ is nonnegative; hence it is not the case
that one factor is strictly positive and the other strictly negative, which
is exactly the stated dichotomy.  (If a factor vanishes, then necessarily
$c_{ij} = 0$; a pair with $\ell_i = t$ is solid or thin according to the
sign of $\ell_j - t$, and both when $\ell_j = t$.) $\square$
**Status: PROVED.**

**Corollary 5.2 (solid good pairs are spectrally good).**  *Let $k \ge 2$
and let $\{i,j\}$ be a good pair with $\ell_i \ge t$ and $\ell_j \ge t$.
Then its $2\times2$ Gram matrix satisfies*
$$\Gamma = \begin{pmatrix}\ell_i & c_{ij}\\ c_{ij} & \ell_j\end{pmatrix}
\;\succeq\; t\,I_2, \qquad\text{in particular}\quad
\lambda_{\min}(\Gamma) \ge t = \frac{2k}{(k+2)n}.$$
*For $k = 3$: $\lambda_{\min}(\Gamma) \ge \dfrac{6}{5n} > \dfrac1n$.
More generally $t > 1/n$ for every $k \ge 3$.*

*Proof.*  First, $s \le 1$ for all $k \ge 2$:
$1 - s = \frac{k^2 - k - 2}{k^2} = \frac{(k-2)(k+1)}{k^2} \ge 0$.
Consider $H = \Gamma - tI_2$.  Its diagonal entries $\ell_i - t,\,
\ell_j - t$ are nonnegative by solidity, so $\operatorname{Tr}H \ge 0$; and
$$\det H = (\ell_i - t)(\ell_j - t) - c_{ij}^2
\;\ge\; (1 - s)(\ell_i - t)(\ell_j - t) \;\ge\; 0,$$
using goodness ($c_{ij}^2 \le s(\ell_i-t)(\ell_j-t)$) and then $s \le 1$
with solidity.  A symmetric $2\times2$ matrix with $\operatorname{Tr} \ge 0$
and $\det \ge 0$ has both eigenvalues $\ge 0$ (they are
$\frac{\operatorname{Tr}H \pm \sqrt{\operatorname{Tr}^2H - 4\det H}}2$;
nonnegative determinant makes them equal-signed or zero, nonnegative trace
makes them nonnegative).  So $H \succeq 0$, i.e. $\Gamma \succeq tI_2$.
Finally $t > 1/n \iff \frac{2k}{k+2} > 1 \iff k > 2$; at $k = 3$,
$t = \frac{6}{5n}$. $\square$  **Status: PROVED.**

**Remark 5.3 (the solidity hypothesis is genuinely needed and is NOT
supplied by Theorem 1 — open).**  Theorem 1 guarantees a good pair, and by
Corollary 5.1 it is solid or thin, but it does **not** guarantee a *solid*
good pair.  A thin good pair can be spectrally worthless: e.g. if $A$ has
two zero rows $i, j$ (allowed on $\mathrm{St}(n,k)$ for $n \ge k+2$), then
$c_{ij} = 0 \le s t^2 = s(\ell_i - t)(\ell_j - t)$ makes $\{i,j\}$ a good
pair with $\lambda_{\min}(\Gamma) = 0$.  Whether some configuration has
**all** of its good pairs thin is an **open question** (not claimed here
either way): adversarial searches in round 1 found no such configuration
(`attacks/round-1-squaring-k3.md`, ledger item "thin pairs"), and in every
Haar sample tested a solid good pair existed (§7, check C4, fraction 1.000
across the grid).  Note the Case-A reduction (`proofs/case-a.md`) does *not*
close this: it removes rows with $\ell \le 1/n$, while thinness is
$\ell \le t = \frac{2k}{(k+2)n} > \frac1n$ for $k \ge 3$, so the band
$\ell \in (1/n, t]$ survives.  Consequently the *usable* $k = 3$ statement
is exactly Corollary 5.2: **any** solid good pair — and Theorem 1 plus
"not all good pairs are thin" would always provide one — has
$\lambda_{\min} \ge 6/(5n) > 1/n$.

> **Scope reminder.**  Even granting a solid good pair, this is a
> $2\times2$ statement.  GTZ$(n,3)$ needs a $3\times3$ submatrix with
> $\lambda_{\min} \ge 1/n$; the pair → triple bridge (conjecture C_ext of
> the round-1 ledger) is open and is not addressed by this document.

---

## 6. Sharpness at (6,3): the icosahedral ETF kills every stronger threshold

### 6.1 The one-parameter family

For $u > 0$ define
$$s_u = \frac{1+u}{k}, \qquad t_u = \frac{uk}{(1+u)\,n}, \qquad
c_u = \frac{t_u}{n},$$
and the statement
$$P_u(n,k):\quad \text{every } A \in \mathrm{St}(n,k) \text{ has a pair }
i \ne j \text{ with } c_{ij}^2 \le s_u(\ell_i - t_u)(\ell_j - t_u).$$
At $u = 2/k$: $s_{2/k} = \frac{k+2}{k^2} = s$ and
$t_{2/k} = \frac{2k}{(k+2)n} = t$, so $P_{2/k}$ **is** Theorem 1.

*Why larger $u$ is the "stronger threshold".*  $t_u = \frac kn\cdot
\frac{u}{1+u}$ is strictly increasing in $u$, and for $u \le k - 1$ one has
$s_u \le 1$, so the proof of Corollary 5.2 applies verbatim at level $u$:
a solid good pair at level $u$ has $2\times2$ Gram
$\lambda_{\min} \ge t_u$.  Thus within the family, larger $u$ yields a
strictly stronger spectral conclusion, and $P_u$ for $u > 2/k$ are the
strictly stronger thresholds.  (For $k = 3$: $u^* = 2/3$, and levels
$u \in [\,1/2,\, 2\,]$ are the ones with $t_u \ge 1/n$ and $s_u \le 1$.)

**Proposition 6.1 (the provable range of the family).**  *For every
$n \ge k \ge 2$ and every $0 < u \le 2/k$, $P_u(n,k)$ holds.*

*Proof.*  Repeat §§2–3 with
$z_{u,i} = \sqrt{u/k}\,\bigl(\ell_i - \frac kn\bigr)$ and
$G_u = WW^T - z_uz_u^T$ ($W$ unchanged).  Each step is replaced as follows.

*(P1$_u$)/(P4$_u$):* $z_u^T\mathbf 1 = \sqrt{u/k}\,(k - k) = 0$, so
$G_u\mathbf 1 = 0$.

*(P3$_u$):* for $i \ne j$,
$(G_u)_{ij} = c_{ij}^2 - \frac{\ell_i\ell_j}{k}
- \frac uk(\ell_i - \frac kn)(\ell_j - \frac kn)
= c_{ij}^2 - s_u(\ell_i - t_u)(\ell_j - t_u) - c_u$, by the coefficient
identities
$$\frac1k + \frac uk = s_u, \qquad
\frac uk\cdot\frac kn = \frac un = s_u t_u, \qquad
\frac uk\cdot\frac{k^2}{n^2} = \frac{uk}{n^2} = s_u t_u^2 + c_u,$$
each verified by direct substitution (the last:
$s_ut_u^2 + c_u = \frac un t_u + \frac{t_u}n = \frac{t_u(1+u)}{n}
= \frac{uk}{n^2}$).  Hence goodness at level $u$
$\iff (G_u)_{ij} \le -c_u$.

*(P5$_u$):* $\operatorname{Tr}G_u = \frac{k-1}{k}Q - \frac uk\bigl(Q -
\frac{k^2}{n}\bigr) = \frac{k-1-u}{k}\,Q + \frac{uk}{n}$.  Since
$u \le 2/k \le 1 \le k-1$, the coefficient $\frac{k-1-u}{k} \ge 0$, so with
$Q \ge k^2/n$:
$$\operatorname{Tr}G_u \ge \frac{k-1-u}{k}\cdot\frac{k^2}n + \frac{uk}n
= \frac{k(k-1)}{n} \;\ge\; m\,t_u,$$
where the last inequality is
$\frac{k(k-1)}{n} \ge \frac{(k+2)(k-1)}{2}\cdot\frac{uk}{(1+u)n}
\iff 2(1+u) \ge (k+2)u \iff ku \le 2$, true by hypothesis.  Also
$\operatorname{Tr}G_u \ge k(k-1)/n > 0$.

*(P6$_u$)/(P7$_u$):* unchanged ($G_u \preceq WW^T$), giving
$\lambda_1(G_u) \ge \operatorname{Tr}G_u/m \ge t_u$.

*Perron step:* identical, with $c_u = t_u/n > 0$ and
$M_u\mathbf 1 = (c_u n + \tau_u)\mathbf 1 = (t_u + \tau_u)\mathbf 1$.
$\square$  **Status: PROVED.**  (Symbolic verification of every displayed
identity and of the boundary equivalence $ku \le 2$: §7, check C7.)

### 6.2 The icosahedral equiangular tight frame is on the Stiefel manifold

Let $\varphi = \frac{1+\sqrt5}2$ (so $\varphi^2 = \varphi + 1$) and take the
six icosahedron-diagonal directions
$$v_1 = (0, 1, \varphi),\; v_2 = (0, 1, -\varphi),\; v_3 = (1, \varphi, 0),\;
v_4 = (1, -\varphi, 0),\; v_5 = (\varphi, 0, 1),\; v_6 = (-\varphi, 0, 1),$$
each of squared norm $\|v_i\|^2 = 1 + \varphi^2 = 2 + \varphi$ (each $v_i$
has one zero coordinate, one $\pm1$ coordinate, and one $\pm\varphi$
coordinate; here $\varphi^2 = \varphi + 1$).  Let $A_{\mathrm{ico}} \in \mathbb R^{6\times3}$ have rows
$$a_i = \frac{v_i}{\sqrt{2}\,\sqrt{1+\varphi^2}}, \qquad i = 1,\dots,6 .$$

**Proposition 6.2.**  *(a) $A_{\mathrm{ico}} \in \mathrm{St}(6,3)$.
(b) $\ell_i = \tfrac12$ for all $i$ and $c_{ij}^2 = \tfrac1{20}$ for all
$i \ne j$.*

*Proof.*  (a) Compute $S = \sum_i v_iv_i^T$ entry by entry.  Diagonal: the
squares of the first coordinates are $0,0,1,1,\varphi^2,\varphi^2$, of the
second $1,1,\varphi^2,\varphi^2,0,0$, of the third
$\varphi^2,\varphi^2,0,0,1,1$ — each set sums to $2 + 2\varphi^2$.
Off-diagonal $(1,2)$: only
$v_3, v_4$ have both coordinates nonzero, contributing
$\varphi + (-\varphi) = 0$; similarly $(1,3)$ gets $\varphi - \varphi = 0$
from $v_5, v_6$ and $(2,3)$ gets $\varphi - \varphi = 0$ from $v_1, v_2$.
Hence $S = 2(1+\varphi^2)\,I_3$ and
$$A_{\mathrm{ico}}^TA_{\mathrm{ico}}
= \frac{1}{2(1+\varphi^2)}\,S = I_3 .$$
(b) $\ell_i = \frac{\|v_i\|^2}{2(1+\varphi^2)} = \frac12$.  For the inner
products, every one of the 15 pairs has
$\langle v_i, v_j\rangle = \pm\varphi$.  Indeed, the six vectors form three
"coplanar" pairs $(v_1,v_2), (v_3,v_4), (v_5,v_6)$, one supported on each
coordinate plane; for each such pair the product is
$1\cdot1 + \varphi\cdot(-\varphi) = 1 - \varphi^2 = -\varphi$ (using
$\varphi^2 = \varphi + 1$).  For the twelve cross pairs, the supports
$\{2,3\}, \{2,3\}, \{1,2\}, \{1,2\}, \{1,3\}, \{1,3\}$ of
$v_1,\dots,v_6$ from *different* coplanar pairs intersect in exactly one
coordinate, on which one vector has entry $\pm1$ and the other $\pm\varphi$
(the $\pm1$ entry of $v_1, v_2$ sits in coordinate 2, of $v_3, v_4$ in
coordinate 1, of $v_5, v_6$ in coordinate 3, never matching the other
vector's $\pm1$ slot), so the product is a single term $\pm\varphi$.
Using $(1+\varphi^2)^2 = (2+\varphi)^2 = 4 + 4\varphi + \varphi^2
= 5 + 5\varphi = 5\varphi^2$:
$$c_{ij}^2 = \frac{\langle v_i, v_j\rangle^2}{4(1+\varphi^2)^2}
= \frac{\varphi^2}{4\cdot 5\varphi^2} = \frac1{20}. \qquad\square$$
**Status: PROVED.**  (Exact sympy verification of (a) and of all 15 margins:
§7, check C6.  This is the standard 6-line equiangular tight frame in
$\mathbb R^3$, $\cos^2\theta = 1/5$.)

### 6.3 All 15 margins vanish at u* = 2/3, and are negative beyond

**Proposition 6.3.**  *At $A_{\mathrm{ico}}$ (so $(n,k) = (6,3)$), the
level-$u$ margin of **every** one of the 15 pairs equals*
$$\mathrm{marg}_u(i,j) = s_u\Bigl(\tfrac12 - t_u\Bigr)^2 - \tfrac1{20}
= \frac{2 - 3u}{60\,(1+u)} \qquad (u > 0).$$
*Hence all 15 margins vanish simultaneously at $u^* = 2/3 = 2/k$
(the level of Theorem 1, $t = \tfrac{6}{5n} = \tfrac15$), and every margin
is strictly negative for every $u > 2/3$.  In particular $P_u(6,3)$ is
**FALSE for every $u > 2/3$**.*

*Proof.*  By Proposition 6.2 all pairs have $\ell_i = \ell_j = \frac12$ and
$c_{ij}^2 = \frac1{20}$, so all 15 margins are equal.  With $k = 3$,
$n = 6$: $t_u = \frac{3u}{6(1+u)} = \frac{u}{2(1+u)}$, hence
$\frac12 - t_u = \frac{1+u-u}{2(1+u)} = \frac{1}{2(1+u)}$ and
$$s_u\Bigl(\tfrac12 - t_u\Bigr)^2 = \frac{1+u}{3}\cdot\frac1{4(1+u)^2}
= \frac{1}{12(1+u)},$$
so
$$\mathrm{marg}_u = \frac1{12(1+u)} - \frac1{20}
= \frac{5 - 3(1+u)}{60(1+u)} = \frac{2-3u}{60(1+u)} .$$
The numerator $2 - 3u$ is zero at $u = 2/3$ and strictly negative for
$u > 2/3$, while $60(1+u) > 0$.  A negative margin at every pair means no
good pair at level $u$ exists at $A_{\mathrm{ico}}$, i.e. $P_u(6,3)$ fails.
$\square$  **Status: PROVED.**  (Exact arithmetic: §7, check C6 — margins at
$u = \frac34, 1, 2, 10$ are $-\frac1{420}, -\frac1{120}, -\frac1{45},
-\frac7{165}$.)

**Corollary 6.4 (sharpness at (6,3)).**  *For $k = 3$ the family threshold
of Theorem 1 is optimal at $n = 6$: $P_u(6,3)$ holds for all
$0 < u \le 2/3$ (Proposition 6.1) and fails for all $u > 2/3$
(Proposition 6.3).  The provable and refuted ranges meet exactly at
$u^* = 2/3$, where the icosahedral ETF attains equality in **all** 15 pairs
simultaneously.  The constant $t = \frac{6}{5n}$ is therefore sharp at
$(6,3)$ — it is a property of the statement, not an artifact of the proof.*
**Status: PROVED** (immediate from 6.1 + 6.3).

**Remark 6.5 (method ceiling, for context).**  The proof mechanism itself
cannot pass $u = 2/k$ at *any* $(n,k)$: at any uniform-leverage
configuration ($\ell_i \equiv k/n$, i.e. $Q = k^2/n$) the trace identity of
(P5$_u$) gives $\operatorname{Tr}G_u = \frac{k(k-1)}{n}$ exactly
(independent of $u$), while $m\,t_u > \frac{k(k-1)}{n}$ for $u > 2/k$; so
the chain $\lambda_1 \ge \operatorname{Tr}G_u/m \ge t_u$ breaks.
(Equal-norm tight frames — uniform-leverage Stiefel points — exist for
every $n \ge k$, e.g. by the frame-potential argument of Benedetto–Fickus,
*Finite normalized tight frames*, Adv. Comput. Math. 18 (2003);
at $(6,3)$ the icosahedron is one.)  This alone would only indict the
method; Proposition 6.3 shows that at $(6,3)$ the *statement* itself fails
beyond $u = 2/3$.  For $n \ge 7$, $k = 3$ the true optimal pair threshold
may exceed $6/(5n)$ (round-1 adversarial data: min margin at $u = 1$ stays
positive at $(8,3)$), but no such improvement is claimed here.
**Status of the displayed computation: PROVED;** the last sentence is an
observation, not a claim.

---

## 7. Numeric sanity block

Worker sanity study: `numerics/study-pair-lemma/check_pair_lemma.py`
(deterministic, seed 20260717; library: `numerics/gtz/` Haar sampler).
Command and observed output (run 2026-07-17, workspace `.venv`,
numpy 2.5.1 / scipy 1.18.0 / sympy 1.14.0):

```
$ .venv/bin/python numerics/study-pair-lemma/check_pair_lemma.py
(C1)-(C4) float checks, 300 Haar samples per (n,k):
     (n,k)   maxIdErr  #posOK  lam1>=t  minBestMargin  solidFrac  minSolidLam
     (2,2)   1.11e-15    True     True   2.500000e-01          -            -
     (4,2)   1.11e-15    True     True   2.484367e-02          -            -
     (8,2)   8.19e-16    True     True   2.297070e-02          -            -
     (6,3)   7.77e-16    True     True   5.611325e-02      1.000     1.217442
     (7,3)   7.42e-16    True     True   4.903009e-02      1.000     1.209291
    (10,3)   7.36e-16    True     True   3.681588e-02      1.000     1.207032
    (30,3)   6.80e-16    True     True   9.283709e-03      1.000     1.202986
     (8,4)   8.88e-16    True     True   5.436995e-02      1.000     1.345245
    (10,5)   1.33e-15    True     True   4.987188e-02      1.000     1.439554
    (14,6)   1.33e-15    True     True   3.975761e-02      1.000     1.517865
    (20,6)   1.11e-15    True     True   2.500332e-02      1.000     1.506780
(C5) k = 2 dictionary  2*G == G_SP:
  s = 1, t = 1/n exactly; max |2G - G_SP| = 2.22e-16 over 400 samples, n in {3,5,9,17}
(C6) icosahedral ETF(6,3), exact arithmetic:
  A^T A = I_3 exact: True;  all 15 margins == 0 at u* = 2/3: True
  margin(u) == (2-3u)/(60(1+u)): True
  margin at u = 3/4, 1, 2, 10: [-1/420, -1/120, -1/45, -7/165] (all < 0)
(C7) symbolic parametric identities (sympy):
  all parametric identities + boundary (2 - ku) + u = 2/k specialization: True

ALL CHECKS PASSED
```

Interpretation against the claims:

* **C1 ↔ (P1)–(P5):** all framework identities hold to $\le 1.4\cdot10^{-15}$
  over 3300 Haar samples, $k$ up to 6, $n$ up to 30.
* **C2 ↔ (P6)–(P7):** $\operatorname{rank}W \le m$, at most $m$ positive
  eigenvalues, $\lambda_1(G) \ge \operatorname{Tr}G/m \ge t$ on every sample.
* **C3 ↔ Theorem 1:** best pair margin $\ge 0$ on every sample (column
  `minBestMargin` is the *minimum over samples of the best* margin — always
  strictly positive; smallest observed $9.3\cdot10^{-3}$ at $(30,3)$).
* **C4 ↔ §5:** every sample had a solid good pair (`solidFrac` = 1.000),
  and every solid good pair had $n\,\lambda_{\min}(\Gamma) \ge n\,t
  = \frac{2k}{k+2}$ (observed minima 1.2174/1.2093/1.2070/1.2030 against
  the bound 1.2 at $k=3$; 1.3452 vs 4/3 at $k=4$; 1.4396 vs 10/7 at $k=5$;
  1.5179/1.5068 vs 3/2 at $k=6$).  This tests the *inequality* of
  Corollary 5.2; the hypothesis coverage (`solidFrac`) is empirical only —
  see Remark 5.3.
* **C5 ↔ §4:** the k = 2 dictionary $2G = G^{\mathrm{SP}}$ holds to machine
  precision; $s = 1$, $t = 1/n$ exactly.
* **C6 ↔ §6.2–6.3:** exact (sympy) — $A_{\mathrm{ico}}^TA_{\mathrm{ico}}
  = I_3$, all 15 margins identically zero at $u^* = 2/3$, the closed form
  $\mathrm{marg}_u = \frac{2-3u}{60(1+u)}$, and strict negativity at sample
  levels $u > 2/3$.
* **C7 ↔ (P3), (P5), Prop 6.1:** the parametric coefficient identities, the
  trace identity, the boundary equivalence
  $\frac{k(k-1)}{n} \ge m t_u \iff ku \le 2$, and the $u = 2/k$
  specialization, verified *symbolically* in $(k, n, u)$.

Additional independent support from round 1 (not rerun here):
`numerics/study-attack-r1-squaring-k3/check_framework.py` (400 Haar samples
per $(n,k)$ on an 11-point grid, identities at $10^{-15}$, pair existence
never violated) and `final_checks.py` (exact icosahedron, C_ext scans).

---

## 8. Status summary

| statement | status |
|---|---|
| (P1)–(P7) framework identities and bounds | **PROVED** |
| Lemma 3.1 (Perron package) | **PROVED** (external input: Horn–Johnson, Matrix Analysis 2nd ed., Thm 8.4.4) |
| **Theorem 1 (General-k Pair Lemma)** | **PROVED**, unconditional |
| Prop 4.1 (exact reduction to Sengupta–Pautov at k = 2; $G = \frac12 G^{\mathrm{SP}}$) | **PROVED** |
| Cor 5.1 (solid/thin dichotomy) | **PROVED** |
| Cor 5.2 (solid good pair ⇒ $\Gamma \succeq tI_2$; k = 3: $\lambda_{\min} \ge \frac6{5n} > \frac1n$) | **PROVED** (conditional on solidity by its very statement) |
| "Some good pair is solid" | **OPEN — not claimed** (Remark 5.3; numerically supported only) |
| Prop 6.1 (family $P_u$ holds for $0 < u \le 2/k$) | **PROVED** |
| Prop 6.2 (icosahedral ETF ∈ St(6,3), $\ell = \frac12$, $c^2 = \frac1{20}$) | **PROVED** (exact arithmetic) |
| Prop 6.3 + Cor 6.4 (all 15 margins vanish at $u^* = \frac23$; $\mathrm{marg}_u = \frac{2-3u}{60(1+u)}$; every $u > \frac23$ FALSE at (6,3); sharpness) | **PROVED** (exact arithmetic) |
| Pair → k-subset bridge (GTZ(n,k), k ≥ 3, from this lemma) | **NOT CLAIMED — OPEN** (see §0; C_ext remains conjectural) |

No GAP markers: every claim labeled PROVED above is proved in full in this
document.  The two OPEN lines are explicitly outside what is claimed.
