# Case A for all $k$: reduction of small-leverage rows

**Deliverable:** generalized Case-A reduction lemma — if $A \in \mathbb{R}^{n\times k}$, $A^TA = I_k$, has a row with leverage $\ell_i \le 1/n$, then $\mathrm{GTZ}(n-1,k) \Rightarrow \mathrm{GTZ}(n,k)$ on such configurations — for **all** $1 \le k < n$, together with the corollary reducing $\mathrm{GTZ}$ to *core configurations*.

**Provenance.** For $k=2$ this is Case A of Sengupta–Pautov, arXiv:2604.05944. The proof below is fully self-contained (it does not rely on any statement of that preprint) and works for every $k$; the $k=2$ case of 2604.05944 is recovered as a special case.

**Notation** (shared, per project MEMORY.md): rows of $A$ are $r_1,\dots,r_n \in \mathbb{R}^k$ (so $A$ has $r_i^T$ as its $i$-th row); leverages $\ell_i = \|r_i\|^2$, with $\sum_i \ell_i = \operatorname{tr}(A A^T) = k$; for $I \subseteq \{1,\dots,n\}$ with $|I| = k$, $A_I$ is the $k\times k$ submatrix of the rows indexed by $I$ (all $k$ columns). Singular values of a matrix $M$ with $k$ columns are written $\sigma_1(M) \ge \dots \ge \sigma_k(M)$, and $\sigma_{\min} = \sigma_k$, $\sigma_{\max} = \sigma_1$.

$\mathrm{GTZ}(n,k)$ denotes the statement: *every $A \in \mathbb{R}^{n\times k}$ with $A^TA = I_k$ has a $k$-subset $I$ with $\sigma_{\min}(A_I) \ge 1/\sqrt{n}$.* For $n = k$ the only choice is $I = \{1,\dots,k\}$ and $A$ itself is orthogonal, so $\sigma_{\min}(A) = 1 \ge 1/\sqrt{k}$: $\mathrm{GTZ}(k,k)$ holds trivially. This trivial case serves as the base of the induction in Corollary 4.

---

## Statements

**Lemma 1 (invariance under row permutation and right rotation).**
Let $A \in \mathbb{R}^{n\times k}$ with $A^TA = I_k$, let $\Pi \in \mathbb{R}^{n\times n}$ be a permutation matrix, realizing a bijection $\pi$ of $\{1,\dots,n\}$ with the convention that row $j$ of $\Pi A$ is row $\pi(j)$ of $A$, and let $Q \in O(k)$. Then $A' = \Pi A Q$ satisfies $A'^TA' = I_k$; the leverages of $A'$ are the leverages of $A$ permuted ($\ell'_j = \ell_{\pi(j)}$); and for every $k$-subset $I$, $A'_I = A_{\pi(I)}\, Q$, so $\sigma_j(A'_I) = \sigma_j(A_{\pi(I)})$ for all $j$. In particular the validity of $\mathrm{GTZ}(n,k)$-type conclusions is invariant under $A \mapsto \Pi A Q$.

**Lemma 2 (multiplicative two-sided bounds on all singular values).**
Let $M \in \mathbb{R}^{m\times k}$ ($m \ge 1$, $k \ge 1$) and let $D \in \mathbb{R}^{k\times k}$ be invertible. Then for every $j = 1,\dots,\min(m,k)$:
$$\sigma_j(M)\,\sigma_{\min}(D) \;\le\; \sigma_j(MD) \;\le\; \sigma_j(M)\,\sigma_{\max}(D).$$

**Lemma 3 (Case-A reduction, all $k$).**
Let $1 \le k < n$ and let $A \in \mathbb{R}^{n\times k}$ satisfy $A^TA = I_k$. Suppose some row $i$ has $\ell_i = \|r_i\|^2 \le 1/n$, and set $b = \sqrt{\ell_i} \ge 0$. Then there exist $Q \in O(k)$ and $\tilde{B} \in \mathbb{R}^{(n-1)\times k}$ with $\tilde{B}^T\tilde{B} = I_k$ such that:

**(i) (construction)** Row $i$ of $A' := AQ$ equals $(b, 0, \dots, 0)$, i.e. $b\,e_1^T$; and, writing $B \in \mathbb{R}^{(n-1)\times k}$ for $A'$ with row $i$ deleted,
$$\tilde{B} = B\,D_t, \qquad D_t = \operatorname{diag}(t, 1, \dots, 1), \qquad t = \frac{1}{\sqrt{1-\ell_i}},$$
which is well defined because $\ell_i \le 1/n < 1$.

**(ii) (pullback with two-sided control of all singular values)** Let $\varphi : \{1,\dots,n-1\} \to \{1,\dots,n\}\setminus\{i\}$ be the order-preserving bijection ($\varphi(j) = j$ for $j < i$, $\varphi(j) = j+1$ for $j \ge i$), matching each row of $\tilde{B}$ (equivalently of $B$) to the row of $A$ it came from. Then every $k$-subset $J \subseteq \{1,\dots,n-1\}$ pulls back to the $k$-subset $\varphi(J) \subseteq \{1,\dots,n\}\setminus\{i\}$ — in particular $i \notin \varphi(J)$ — and for **every** $j = 1,\dots,k$:
$$\sqrt{1-\ell_i}\;\sigma_j(\tilde{B}_J) \;\le\; \sigma_j\!\big(A_{\varphi(J)}\big) \;\le\; \sigma_j(\tilde{B}_J).$$

**(iii) (composition; the constant self-reproduces exactly)** If $\sigma_{\min}(\tilde{B}_J) \ge 1/\sqrt{n-1}$, then
$$\sigma_{\min}\!\big(A_{\varphi(J)}\big) \;\ge\; \sqrt{1-\ell_i}\cdot\frac{1}{\sqrt{n-1}} \;\ge\; \sqrt{\frac{n-1}{n}}\cdot\frac{1}{\sqrt{n-1}} \;=\; \frac{1}{\sqrt{n}}.$$

Consequently: if $\mathrm{GTZ}(n-1,k)$ holds, then every $A$ as above (orthonormal columns, some $\ell_i \le 1/n$) has a $k$-subset $I \subseteq \{1,\dots,n\}\setminus\{i\}$ with $\sigma_{\min}(A_I) \ge 1/\sqrt{n}$.

**Corollary 4 (reduction of GTZ to core configurations).**
Fix $k \ge 1$. Call $A \in \mathbb{R}^{m\times k}$ with $A^TA = I_k$, $m > k$, a **core configuration** if every row satisfies $\ell_i > 1/m$. Assume:

> **(Core hypothesis for $k$.)** For every integer $m > k$ and every core configuration $A \in \mathbb{R}^{m\times k}$, there exists a $k$-subset $I$ with $\sigma_{\min}(A_I) \ge 1/\sqrt{m}$.

Then $\mathrm{GTZ}(n,k)$ holds for **every** $n \ge k$. Moreover, for a fixed target $n$, only the instances of the core hypothesis with $k < m \le n$ are used.

---

## Proofs

### Proof of Lemma 1

$A'^TA' = Q^T A^T \Pi^T \Pi A Q = Q^T (A^TA) Q = Q^TQ = I_k$, since $\Pi^T\Pi = I_n$ and $Q \in O(k)$.

With the stated convention (row $j$ of $\Pi A$ is row $\pi(j)$ of $A$), row $j$ of $A'$ is $r_{\pi(j)}^T Q$, whose squared norm is $\|Q^T r_{\pi(j)}\|^2 = \|r_{\pi(j)}\|^2 = \ell_{\pi(j)}$ ($Q$ orthogonal); so leverages are permuted by $\pi$. For a $k$-subset $I$, selecting rows commutes with right multiplication, so $A'_I = (\Pi A)_I\, Q = A_{\pi(I)}\, Q$ where $\pi(I) = \{\pi(j) : j \in I\}$. Right multiplication by an orthogonal matrix preserves all singular values ($(A'_I)(A'_I)^T = A_{\pi(I)} Q Q^T A_{\pi(I)}^T = A_{\pi(I)}A_{\pi(I)}^T$, and the nonzero spectrum of $MM^T$ determines $\sigma(M)$; here both are $k\times k$ so all singular values agree). Since $I \mapsto \pi(I)$ is a bijection on $k$-subsets, the family $\{\sigma(A'_I)\}_I$ equals the family $\{\sigma(A_I)\}_I$ up to relabeling. $\blacksquare$

### Proof of Lemma 2

Both bounds follow from the Courant–Fischer (min–max) characterization of singular values (standard; e.g. Horn & Johnson, *Matrix Analysis*, §7.3): for $X \in \mathbb{R}^{m\times k}$ and $1 \le j \le \min(m,k)$,
$$\sigma_j(X) \;=\; \max_{\substack{S \subseteq \mathbb{R}^k \\ \dim S = j}}\; \min_{0 \ne v \in S} \frac{\|Xv\|}{\|v\|}.$$
This characterization requires no assumption of distinct singular values.

*Lower bound.* Let $S^\ast \subseteq \mathbb{R}^k$, $\dim S^\ast = j$, attain the maximum for $M$: $\min_{0\ne w \in S^\ast} \|Mw\|/\|w\| = \sigma_j(M)$. Since $D$ is invertible, $S := D^{-1}S^\ast$ has dimension $j$. For any $0 \ne v \in S$, the vector $w := Dv$ lies in $S^\ast$ and is nonzero, so
$$\frac{\|MDv\|}{\|v\|} \;=\; \frac{\|Mw\|}{\|w\|}\cdot\frac{\|Dv\|}{\|v\|} \;\ge\; \sigma_j(M)\,\sigma_{\min}(D),$$
using $\|Dv\| \ge \sigma_{\min}(D)\|v\|$. Taking the min over $0 \ne v \in S$ and then noting that $\sigma_j(MD)$ is the max over all $j$-dimensional subspaces, $\sigma_j(MD) \ge \sigma_j(M)\,\sigma_{\min}(D)$.

*Upper bound.* Apply the lower bound to the pair $(MD,\, D^{-1})$: $\sigma_j(M) = \sigma_j\big((MD)D^{-1}\big) \ge \sigma_j(MD)\,\sigma_{\min}(D^{-1}) = \sigma_j(MD)/\sigma_{\max}(D)$, since $\sigma_{\min}(D^{-1}) = 1/\sigma_{\max}(D)$ for invertible $D$. Rearranging ($\sigma_{\max}(D) > 0$) gives $\sigma_j(MD) \le \sigma_j(M)\,\sigma_{\max}(D)$. $\blacksquare$

### Proof of Lemma 3

Throughout, $n \ge 2$ (from $n > k \ge 1$), hence $\ell_i \le 1/n \le 1/2 < 1$: all divisions by $1-\ell_i$ below are legitimate, and $t \in [1, \sqrt{2}]$ is finite.

**Step 1 (rotation of row $i$ onto the first coordinate axis).**
*Case $r_i \ne 0$:* let $q_1 = r_i/\|r_i\|$ and extend to an orthonormal basis $q_1, q_2, \dots, q_k$ of $\mathbb{R}^k$ (possible for any $k \ge 1$; for $k = 1$ the basis is just $q_1 = \pm 1$). Let $Q = [\,q_1 \mid \cdots \mid q_k\,] \in O(k)$. Row $i$ of $A' := AQ$ is
$$r_i^T Q = \big(r_i^Tq_1,\, r_i^Tq_2,\, \dots,\, r_i^Tq_k\big) = \big(\|r_i\|, 0, \dots, 0\big) = b\,e_1^T,$$
since $r_i^Tq_1 = \|r_i\|$ and $q_s \perp q_1 \parallel r_i$ for $s \ge 2$.
*Case $r_i = 0$ (i.e. $\ell_i = 0$):* take $Q = I_k$; row $i$ of $A' = A$ is $0 = b\,e_1^T$ with $b = 0$. (A zero row already lies "on the axis"; no rotation is needed and none is possible in the sense of normalizing $r_i$, which is why this case is split off.)

In both cases $A'^TA' = I_k$ and, by Lemma 1 (with $\Pi = I$), for every $k$-subset $I$: $A'_I = A_I\,Q$ and $\sigma_j(A'_I) = \sigma_j(A_I)$ for all $j$. So it suffices to prove the singular-value claims for submatrices of $A'$.

**Step 2 (deletion and column renormalization; $\tilde{B}$ has orthonormal columns).**
Let $B \in \mathbb{R}^{(n-1)\times k}$ be $A'$ with row $i$ deleted, and let $c_1, \dots, c_k \in \mathbb{R}^{n-1}$ be the columns of $B$. For each pair of column indices $s, s'$, split the inner product of the corresponding columns of $A'$ into the contribution of row $i$ and the rest:
$$\delta_{ss'} \;=\; (A'^TA')_{ss'} \;=\; c_s^Tc_{s'} \;+\; A'_{is}A'_{is'}.$$
By Step 1, $A'_{i1} = b$ and $A'_{is} = 0$ for $s \ge 2$. Hence:
- $\|c_1\|^2 = 1 - b^2 = 1 - \ell_i > 0$;
- $\|c_s\|^2 = 1$ for $s \ge 2$;
- $c_s^Tc_{s'} = 0$ for all $s \ne s'$ (the row-$i$ term $A'_{is}A'_{is'}$ vanishes whenever $s \ne s'$, because at least one factor has column index $\ge 2$).

Define
$$\tilde{B} := B\,D_t, \qquad D_t = \operatorname{diag}(t, 1, \dots, 1), \qquad t = (1-\ell_i)^{-1/2} \ge 1.$$
The columns of $\tilde{B}$ are $t\,c_1, c_2, \dots, c_k$; by the three bullet points, they are orthonormal: $\|t c_1\|^2 = t^2(1-\ell_i) = 1$, and all other norms and inner products are unchanged. So $\tilde{B}^T\tilde{B} = I_k$, and $\tilde{B}$ is a legitimate $\mathrm{GTZ}(n-1,k)$ instance because $k \le n-1$.

Degenerate subcases, explicitly:
- **$\ell_i = 0$:** then $b = 0$, $t = 1$, $D_t = I_k$, and $\tilde{B} = B$ is simply $A$ with the zero row deleted (deleting a zero row does not change any column inner product). All formulas below hold verbatim with $\sqrt{1-\ell_i} = 1$; the reduction is pure deletion and part (ii) becomes the equality $\sigma_j(A_{\varphi(J)}) = \sigma_j(\tilde{B}_J)$ for all $j$.
- **$k = 1$:** the list $c_2, \dots, c_k$ is empty; $A$ is a unit vector in $\mathbb{R}^n$ with entries $a_1,\dots,a_n$, $\ell_i = a_i^2$, $B \in \mathbb{R}^{(n-1)\times 1}$, $\tilde{B} = t\,B$ has unit norm. Everything degenerates consistently (see also the remark on $\sigma_{\max}(D_t^{-1})$ in Step 4).

**Step 3 (pullback bookkeeping).**
Row $j$ of $B$ is row $\varphi(j)$ of $A'$, by definition of deletion and of the order-preserving bijection $\varphi$ ($\varphi(j) = j$ for $j < i$, $\varphi(j) = j+1$ for $j \ge i$). Hence for every $k$-subset $J \subseteq \{1,\dots,n-1\}$:
$$B_J = A'_{\varphi(J)}, \qquad \varphi(J) \subseteq \{1,\dots,n\}\setminus\{i\}, \quad |{\varphi(J)}| = k,$$
and in particular $i \notin \varphi(J)$: the pulled-back submatrix **never contains the deleted row**, because $J$ indexes rows of the reduced matrix and the range of $\varphi$ omits $i$ by construction. Since right multiplication acts on columns and commutes with row selection, $\tilde{B}_J = (B\,D_t)_J = B_J\,D_t$, so
$$A'_{\varphi(J)} \;=\; B_J \;=\; \tilde{B}_J\,D_t^{-1}, \qquad D_t^{-1} = \operatorname{diag}\big(\sqrt{1-\ell_i},\, 1, \dots, 1\big).$$

*Remark (excluding row $i$ is not merely convenient but essentially forced).* For any $k$-subset $I$ and any $j \in I$, taking the unit coordinate vector picking out row $j$ gives $\sigma_{\min}(A_I) = \sigma_{\min}(A_I^T) \le \|A_I^T e_{(j)}\| = \|r_j\|$, so $\sigma_{\min}(A_I) \le \min_{j\in I}\|r_j\|$. If $\ell_i < 1/n$ strictly, every submatrix containing row $i$ has $\sigma_{\min} \le \sqrt{\ell_i} < 1/\sqrt{n}$ and is useless for GTZ; only in the boundary case $\ell_i = 1/n$ could a submatrix through row $i$ conceivably achieve the bound. The reduction's bookkeeping (search only over $I \not\ni i$) is therefore consistent with the target.

**Step 4 (two-sided control of all singular values).**
$D_t^{-1}$ is invertible with diagonal entries $\sqrt{1-\ell_i} \in (0,1]$ and $1$; hence
$$\sigma_{\min}(D_t^{-1}) = \sqrt{1-\ell_i}, \qquad \sigma_{\max}(D_t^{-1}) \le 1$$
(for $k \ge 2$, $\sigma_{\max}(D_t^{-1}) = 1$ exactly; for $k = 1$, $D_t^{-1} = (\sqrt{1-\ell_i})$ and $\sigma_{\max} = \sigma_{\min} = \sqrt{1-\ell_i} \le 1$ — the stated inequalities hold in both cases). Apply Lemma 2 with $M = \tilde{B}_J \in \mathbb{R}^{k\times k}$ and $D = D_t^{-1}$: for every $j = 1,\dots,k$,
$$\sqrt{1-\ell_i}\;\sigma_j(\tilde{B}_J) \;\le\; \sigma_j\big(\tilde{B}_J D_t^{-1}\big) \;=\; \sigma_j\big(A'_{\varphi(J)}\big) \;\le\; \sigma_{\max}(D_t^{-1})\,\sigma_j(\tilde{B}_J) \;\le\; \sigma_j(\tilde{B}_J).$$
Finally $\sigma_j(A'_{\varphi(J)}) = \sigma_j(A_{\varphi(J)})$ by Step 1. This proves (ii); note it controls **all** $k$ singular values of the pulled-back submatrix, not only the smallest.

**Step 5 (composition of the constants).**
Since $\ell_i \le 1/n$ and $x \mapsto \sqrt{1-x}$ is decreasing on $[0,1]$,
$$\sqrt{1-\ell_i} \;\ge\; \sqrt{1-\tfrac{1}{n}} \;=\; \sqrt{\tfrac{n-1}{n}}.$$
If $\sigma_{\min}(\tilde{B}_J) \ge 1/\sqrt{n-1}$, then by (ii) with $j = k$:
$$\sigma_{\min}\big(A_{\varphi(J)}\big) \;\ge\; \sqrt{1-\ell_i}\cdot\sigma_{\min}(\tilde{B}_J) \;\ge\; \sqrt{\frac{n-1}{n}}\cdot\frac{1}{\sqrt{n-1}} \;=\; \sqrt{\frac{n-1}{n(n-1)}} \;=\; \frac{1}{\sqrt{n}}.$$
The constant $1/\sqrt{\cdot}$ self-reproduces **exactly** — the reduction loses nothing at $\ell_i = 1/n$. This proves (iii).

**Conclusion of Lemma 3.** Assume $\mathrm{GTZ}(n-1,k)$. Since $\tilde{B} \in \mathbb{R}^{(n-1)\times k}$ has $\tilde{B}^T\tilde{B} = I_k$ and $n-1 \ge k$, there is a $k$-subset $J \subseteq \{1,\dots,n-1\}$ with $\sigma_{\min}(\tilde{B}_J) \ge 1/\sqrt{n-1}$. (If $n-1 = k$, this is the trivial statement $\mathrm{GTZ}(k,k)$ with $J = \{1,\dots,k\}$, $\sigma_{\min}(\tilde{B}) = 1$.) By (iii), $I := \varphi(J)$ satisfies $i \notin I$ and $\sigma_{\min}(A_I) \ge 1/\sqrt{n}$. $\blacksquare$

### Proof of Corollary 4

Induction on $n \ge k$, for the fixed $k$.

*Base $n = k$:* $A$ is a $k\times k$ orthogonal matrix, $\sigma_{\min}(A) = 1 \ge 1/\sqrt{k}$; $\mathrm{GTZ}(k,k)$ holds with no hypothesis used.

*Inductive step.* Let $n > k$ and assume $\mathrm{GTZ}(n-1,k)$ (available: $n - 1 \ge k$, and it is either the base case or the previous step). Let $A \in \mathbb{R}^{n\times k}$, $A^TA = I_k$, be arbitrary. Two exhaustive cases:

- **(core)** Every row of $A$ has $\ell_i > 1/n$. Then $A$ is a core configuration with $m = n$, and the core hypothesis (instance $m = n$) directly provides $I$ with $\sigma_{\min}(A_I) \ge 1/\sqrt{n}$.
- **(reducible)** Some row has $\ell_i \le 1/n$. Then Lemma 3, whose hypothesis $\mathrm{GTZ}(n-1,k)$ is the induction hypothesis, provides $I$ with $\sigma_{\min}(A_I) \ge 1/\sqrt{n}$.

Hence $\mathrm{GTZ}(n,k)$ holds, completing the induction.

*Quantifier bookkeeping.* The induction on $n$ threads through **all** intermediate sizes: proving $\mathrm{GTZ}(n,k)$ for a fixed target $n$ consumes the core hypothesis at every $m \in \{k+1, \dots, n\}$ (a reducible configuration at size $m$ defers to size $m-1$, where either case may occur again), and nothing beyond $m = n$. So the correct statement is: $\big[\forall m \in (k, n]:$ core hypothesis at $m\big] \Rightarrow \mathrm{GTZ}(n,k)$; and with the hypothesis for all $m > k$, $\mathrm{GTZ}(n,k)$ follows for all $n \ge k$ simultaneously. Note the corollary does **not** claim any implication in the reverse direction, and the core hypothesis is used only for core matrices — the reduction guarantees no more than that. Note also that a single application of Lemma 3 removes one small-leverage row; the reduced matrix $\tilde{B}$ need not be core (it may again have a row of leverage $\le 1/(n-1)$), which is precisely why the corollary runs a full induction rather than applying Lemma 3 once. $\blacksquare$

---

## Status

| Statement | Status |
|---|---|
| Lemma 1 (invariance) | **PROVED** |
| Lemma 2 (multiplicative singular-value bounds) | **PROVED** |
| Lemma 3 (Case-A reduction, all $1 \le k < n$) | **PROVED** |
| Corollary 4 (reduction to core configurations) | **PROVED** |

No GAP markers. The proof nowhere restricts $k$: the rotation (Step 1), the column bookkeeping (Step 2), the pullback (Step 3), and the singular-value bounds (Step 4) are stated and proved for arbitrary $1 \le k < n$, with the degenerate cases $k = 1$, $\ell_i = 0$, and $n - 1 = k$ treated explicitly. Realness is not essential to this reduction (the same argument works over $\mathbb{C}$ with $Q$ unitary), but the statement is given over $\mathbb{R}$ as that is the GTZ setting.

---

## Numeric sanity block

Every numerically checkable claim above was tested with the workspace `.venv`
(`numpy 2.5.1`, `scipy 1.18.0`). Runnable script (canonical copy, with README):
`numerics/study-case-a-sanity/case_a_sanity.py`. Exact command:

```
.venv/bin/python numerics/study-case-a-sanity/case_a_sanity.py
```

Master seed **20260717** (`numpy.random.default_rng`); deterministic. What is checked:

1. **Lemma 2** directly: 200 random pairs $(M, D)$ — $M \in \mathbb{R}^{m\times k}$ Gaussian, $D \in \mathbb{R}^{k\times k}$ Gaussian invertible (also ill-conditioned diagonal $D$) — verify $\sigma_j(M)\sigma_{\min}(D) - \varepsilon \le \sigma_j(MD) \le \sigma_j(M)\sigma_{\max}(D) + \varepsilon$ for all $j$.
2. **Lemma 3 reduction**, $k \in \{2,3,4\}$ (plus $k = 1$ as a degenerate check), $n \in \{k+1, k+2, k+5, 12, 20\}$, 25 instances per $(n,k)$: build $A$ with a planted small-leverage row via the reverse construction (random orthonormal $\tilde{B}_0$ by QR, random $b \in [0, 1/\sqrt{n}]$ — forced to $b = 0$ in a fifth of instances to hit $\ell_i = 0$ exactly — insert row $b\,e_1^T$ at a random position, undo the renormalization, scramble by a random $Q_r \in O(k)$ and a row permutation). Then **apply the documented reduction from scratch** (find a row with $\ell_i \le 1/n$, build $Q$ from Step 1, delete, renormalize per Step 2) and assert, for all $k$-subsets $J$ when $\binom{n-1}{k} \le 300$ and 200 random $J$ otherwise, **for all** $j = 1,\dots,k$:
$$\sqrt{1-\ell_i}\,\sigma_j(\tilde{B}_J) - \varepsilon \;\le\; \sigma_j(A_{\varphi(J)}) \;\le\; \sigma_j(\tilde{B}_J) + \varepsilon,$$
with $\varepsilon = 10^{-10}$; also $\|A^TA - I\|_\infty, \|\tilde{B}^T\tilde{B} - I\|_\infty < 10^{-12}$, row $i$ of $AQ$ equals $b\,e_1^T$ within $10^{-12}$, and the composition implication ($\sigma_{\min}(\tilde{B}_J) \ge 1/\sqrt{n-1} \Rightarrow \sigma_{\min}(A_{\varphi(J)}) \ge 1/\sqrt{n} - \varepsilon$) on every tested $J$.
3. **Composition identity** $\sqrt{(n-1)/n}\cdot(1/\sqrt{n-1}) = 1/\sqrt{n}$ to machine precision for $n = 2,\dots,50$.
4. **Row-norm remark** of Step 3: $\sigma_{\min}(A_I) \le \min_{j \in I}\|r_j\|$ on random submatrices.

**Observed result (run of 2026-07-17, exit code 0):** all checks pass — verbatim output:

```
Lemma 2: 200/200 pairs OK (max violation beyond tol 0.0e+00)
Lemma 3: k in [1, 2, 3, 4], 500 instances, 91175 (J,j) inequality checks OK (worst lower-slack 3.33e-16, worst upper-slack 4.44e-16)
         ell_i = 0 instances: 100, max |sigma_j(A_phiJ)-sigma_j(Bt_J)| = 0.00e+00 (< 1e-12)
Composition identity: max |sqrt((n-1)/n)/sqrt(n-1) - 1/sqrt(n)| = 1.1e-16 (n<=50)
Row-norm bound sigma_min(A_I) <= min_{j in I} ||r_j||: 100/100 OK
ALL SANITY CHECKS PASSED
```

The residual slacks ($\sim 4\cdot 10^{-16}$) are floating-point noise, far inside the $10^{-10}$ tolerance. The script exits nonzero on any violation; see `numerics/study-case-a-sanity/README.md` for commands and seed.

**Label:** the numeric block is a sanity check of PROVED statements, not evidence for them; per project rules the proofs above stand on the written arguments alone.
