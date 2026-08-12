# Theorem C2: the $p^2$-weighted average inequality at $k = 2$

**Target:** VAL-THM-C2-001. Promoted from attack round 3
(`attacks/round-3-variational.md` §§1–3; original numerics
`numerics/study-attack-r3-variational/`). Fresh numeric sanity for this
document: `numerics/study-theorem-c2/` (§8). Notation follows MEMORY.md.

---

## 0. Scope, status, inputs

**Scope guard (binding).** Every statement in this document concerns $k = 2$
only: matrices $A \in \mathrm{St}(n,2)$, i.e. real $n \times 2$ matrices with
$A^TA = I_2$. The general-$k$ analog C2k
($\sum_{i<j}\|r_i \wedge r_j\|^2[s(\ell_i-t)(\ell_j-t)-c_{ij}^2] \ge 0$)
is **CONJECTURED**, lives exclusively in the attack ledger
(`attacks/round-3-variational.md` §7), and is **not** claimed or used here.
Moreover, the volume-weighted $k$-subset analogs of C2 are provably **FALSE**
for $k \ge 3$ (ledger §6, obstructions O1/O2, stored counterexamples), so no
extension to $k \ge 3$ should be presumed from anything below.

**Status.**

| statement | status |
|---|---|
| Lemmas 1.1, 2.1–2.3, 3.1 (identities, spectra) | PROVED |
| Theorem 4.1 (C2, master identity, equality characterization) | PROVED (unconditional) |
| Proposition 4.3 (equality at the S–P extremal family) | PROVED |
| Corollary 5.1 (Perron-free pair existence, $k=2$) | PROVED (unconditional) |
| Remark 5.3 (Perron-free core-case bound for GTZ$(n,2)$) | PROVED |
| Corollary 6.1 (planar lemma, realizability-free) | PROVED |
| Propositions 7.1–7.2 (uniform-average failure and its anatomy) | PROVED |
| Remark 7.3 (expected-charpoly form of C2) | PROVED |
| Remark 7.4 (apolar-invariant identity; aside, used nowhere) | PROVED |
| C2k (general $k$) | CONJECTURED — ledger only, not part of this document |

**Inputs used (complete list).** Column orthonormality of $A$; elementary
complex-number algebra; the spectral theorem for real symmetric matrices in
the form of the Rayleigh bound $x^TKx \le \lambda_{\max}(K)\,\|x\|^2$ with its
standard equality case; the fact that $WW^T$ and $W^TW$ have the same nonzero
eigenvalues; Cauchy–Schwarz in the form $(\sum_i \ell_i)^2 \le n\sum_i
\ell_i^2$. **No Perron–Frobenius theory, no induction on $n$, and no external
paper is used anywhere in the proofs** (non-circularity is spelled out in
§5.2).

---

## 1. Setup: the squared-row dictionary

Fix $n \ge 2$ and $A = [\,a \;\; b\,] \in \mathrm{St}(n,2)$, so
$$\sum_i a_i^2 = 1, \qquad \sum_i b_i^2 = 1, \qquad \sum_i a_ib_i = 0. \tag{1.1}$$
Rows $r_i = (a_i, b_i) \in \mathbb R^2$, leverages $\ell_i = \|r_i\|^2 = a_i^2
+ b_i^2$, inner products $c_{ij} = \langle r_i, r_j\rangle$, Plücker minors
$p_{ij} = a_ib_j - b_ia_j$ (so $p_{ji} = -p_{ij}$, $p_{ii} = 0$), and
$P = AA^T$ (so $P_{ii} = \ell_i$, $P_{ij} = c_{ij}$). Write
$$Q := \sum_i \ell_i^2 .$$

Complexify: $z_i := a_i + \mathrm i\,b_i \in \mathbb C$ and apply the
squaring map
$$\zeta_i := z_i^2 .$$
Let $w_i := (\operatorname{Re}\zeta_i, \operatorname{Im}\zeta_i) \in \mathbb
R^2$, let $W \in \mathbb R^{n\times 2}$ have rows $w_i$, and set
$$K := WW^T, \qquad L := \ell\ell^T, \qquad
M_2 := \sum_i \zeta_i|\zeta_i|, \qquad M_4 := \sum_i \zeta_i^2 .$$

**Lemma 1.1 (dictionary).** For all $i, j$:

1. $\bar z_i z_j = c_{ij} + \mathrm i\,p_{ij}$; in particular $|z_i|^2 =
   \ell_i$ and $|\zeta_i| = |z_i|^2 = \ell_i$.
2. $\sum_i \zeta_i = 0$ and $\sum_i |\zeta_i| = 2$.
3. $K_{ij} = \operatorname{Re}(\bar\zeta_i\zeta_j) = c_{ij}^2 - p_{ij}^2$,
   and $\ell_i\ell_j = c_{ij}^2 + p_{ij}^2$; in particular $K_{ii} =
   \ell_i^2$.
4. $p_{ij}^2 = \tfrac12(L - K)_{ij}$ for **all** $i,j$ (both sides vanish on
   the diagonal).
5. $W^T\mathbf 1 = 0$, hence $K\mathbf 1 = 0$.

*Proof.* (1) $\bar z_i z_j = (a_i - \mathrm i b_i)(a_j + \mathrm i b_j) =
(a_ia_j + b_ib_j) + \mathrm i(a_ib_j - b_ia_j) = c_{ij} + \mathrm i p_{ij}$.
Taking $j = i$: $|z_i|^2 = \ell_i$ (and $p_{ii} = 0$). Then $|\zeta_i| =
|z_i|^2 = \ell_i$.

(2) $\sum_i \zeta_i = \sum_i (a_i^2 - b_i^2) + 2\mathrm i \sum_i a_ib_i =
(1 - 1) + 2\mathrm i\cdot 0 = 0$ by (1.1); $\sum_i|\zeta_i| = \sum_i \ell_i =
\sum_i a_i^2 + \sum_i b_i^2 = 2$.

(3) Under the standard isometry $\mathbb C \cong \mathbb R^2$
($\langle u, v\rangle_{\mathbb R^2} = \operatorname{Re}\bar uv$), $K_{ij} =
\langle w_i, w_j \rangle = \operatorname{Re}(\bar\zeta_i \zeta_j)$. By (1),
$\bar\zeta_i\zeta_j = (\bar z_iz_j)^2 = (c_{ij} + \mathrm i p_{ij})^2 =
(c_{ij}^2 - p_{ij}^2) + 2\mathrm i\, c_{ij}p_{ij}$, so
$\operatorname{Re}\bar\zeta_i\zeta_j = c_{ij}^2 - p_{ij}^2$. Also
$\ell_i\ell_j = |z_i|^2|z_j|^2 = |\bar z_iz_j|^2 = c_{ij}^2 + p_{ij}^2$.
With $j = i$: $K_{ii} = c_{ii}^2 - 0 = \ell_i^2$.

(4) Subtract the two identities of (3): $(L - K)_{ij} = \ell_i\ell_j -
(c_{ij}^2 - p_{ij}^2) = 2p_{ij}^2$. On the diagonal both sides are $0$ by
(3).

(5) $W^T\mathbf 1 = \sum_i w_i$, which is $\sum_i \zeta_i = 0$ read in
$\mathbb R^2$, by (2). Hence $K\mathbf 1 = W(W^T\mathbf 1) = 0$.
$\blacksquare$

---

## 2. The averaging identities

Define the **Sengupta–Pautov bracket** of a pair $i < j$:
$$b_{ij} := p_{ij}^2 - \frac{\ell_i + \ell_j}{n} + \frac1{n^2}. \tag{2.1}$$

**Lemma 2.1 (row marginal and total mass).** For every $i$:
$\;\sum_{j \ne i} p_{ij}^2 = \ell_i$, and consequently
$\sum_{i<j} p_{ij}^2 = 1$.

*Proof.* Expand and use (1.1) (the $j = i$ term is $0$, so it may be
included):
$$\sum_j p_{ij}^2 = \sum_j (a_ib_j - b_ia_j)^2
= a_i^2\sum_j b_j^2 - 2a_ib_i\sum_j a_jb_j + b_i^2\sum_j a_j^2
= a_i^2 + b_i^2 = \ell_i .$$
Summing over $i$: $\sum_{i \ne j} p_{ij}^2 = \sum_i \ell_i = 2$, and since
$p_{ij}^2 = p_{ji}^2$, $\sum_{i<j} p_{ij}^2 = 1$. (This is also Binet–Cauchy,
$\sum_{i<j}p_{ij}^2 = \det A^TA = 1$; the direct proof keeps the document
self-contained.) $\blacksquare$

**Lemma 2.2 (weighted-average identity; ledger (I1), second half).**
$$\sum_{i<j} p_{ij}^2\, b_{ij} \;=\; \sum_{i<j} p_{ij}^4 \;-\; \frac Qn
\;+\; \frac1{n^2}. \tag{2.2}$$
In particular, since $\sum_{i<j}p_{ij}^2 = 1$ with all weights $\ge 0$, the
left side is a genuine weighted average of the brackets $b_{ij}$.

*Proof.* Expand using (2.1):
$$\sum_{i<j} p_{ij}^2 b_{ij} = \sum_{i<j} p_{ij}^4
- \frac1n \sum_{i<j} p_{ij}^2(\ell_i + \ell_j)
+ \frac1{n^2}\sum_{i<j} p_{ij}^2 .$$
Middle sum, by symmetry of $p_{ij}^2$ and Lemma 2.1:
$$\sum_{i<j} p_{ij}^2(\ell_i + \ell_j) = \tfrac12\sum_{i\ne j}
p_{ij}^2(\ell_i + \ell_j) = \sum_{i\ne j} p_{ij}^2\,\ell_i
= \sum_i \ell_i \sum_{j \ne i} p_{ij}^2 = \sum_i \ell_i^2 = Q .$$
Last sum $= 1$ by Lemma 2.1. $\blacksquare$

**Lemma 2.3 (Frobenius form; ledger (I2)).**
$$\sum_{i<j} p_{ij}^4 = \tfrac18\,\|L - K\|_F^2, \qquad\text{hence}\qquad
\sum_{i<j} p_{ij}^2 b_{ij} = \tfrac18\|L-K\|_F^2 - \frac Qn + \frac1{n^2}.
\tag{2.3}$$

*Proof.* The matrix $(p_{ij}^2)_{ij} = \tfrac12(L-K)$ (Lemma 1.1(4)) is
symmetric with zero diagonal, so
$\sum_{i<j}p_{ij}^4 = \tfrac12 \sum_{i,j} p_{ij}^4 = \tfrac12\,
\bigl\|\tfrac12(L-K)\bigr\|_F^2 = \tfrac18\|L-K\|_F^2$. Combine with (2.2).
$\blacksquare$

**Remark 2.4 (two faces of the bracket).** Using $p_{ij}^2 = \ell_i\ell_j -
c_{ij}^2$ (Lemma 1.1(3)):
$$b_{ij} \;=\; \Bigl(\ell_i - \frac1n\Bigr)\Bigl(\ell_j - \frac1n\Bigr) -
c_{ij}^2 \;=\; \chi_{ij}\!\Bigl(\frac1n\Bigr), \tag{2.4}$$
where $\chi_{ij}(x) = \det(xI_2 - P_{\{i,j\}\{i,j\}}) = x^2 - (\ell_i +
\ell_j)x + p_{ij}^2$ is the characteristic polynomial of the pair Gram block
$P_{\{i,j\}\{i,j\}} = \begin{pmatrix}\ell_i & c_{ij}\\ c_{ij} &
\ell_j\end{pmatrix}$. (First equality: expand $(\ell_i - \tfrac1n)(\ell_j -
\tfrac1n) - c_{ij}^2 = \ell_i\ell_j - c_{ij}^2 - \tfrac{\ell_i+\ell_j}n +
\tfrac1{n^2} = b_{ij}$. Second: evaluate $\chi_{ij}$ at $1/n$.) So
$b_{ij} \ge 0$ is exactly the Sengupta–Pautov "good pair" condition
$c_{ij}^2 \le (\ell_i - 1/n)(\ell_j - 1/n)$.

---

## 3. The spectrum of the squared-row Gram

**Lemma 3.1 (closed form of $W^TW$ and its eigenvalues).**
$$W^TW = \frac12\begin{pmatrix} Q + \operatorname{Re}M_4 &
\operatorname{Im}M_4\\ \operatorname{Im}M_4 & Q - \operatorname{Re}M_4
\end{pmatrix},$$
whose eigenvalues $\mu_1 \ge \mu_2$ satisfy
$$\mu_{1,2} = \tfrac12\bigl(Q \pm |M_4|\bigr), \qquad \mu_1 + \mu_2 = Q,
\qquad \mu_1 - \mu_2 = |M_4|, \qquad \mu_2 \ge 0,$$
and $\mu_1 \ge Q/2 \ge 2/n > 0$. Moreover the eigenvalues of $K = WW^T$ are
$\mu_1, \mu_2$ together with $n - 2$ zeros; in particular
$\lambda_{\max}(K) = \mu_1$.

*Proof.* Entrywise, using $(\operatorname{Re}\zeta)^2 = \tfrac12(|\zeta|^2 +
\operatorname{Re}\zeta^2)$, $(\operatorname{Im}\zeta)^2 = \tfrac12(|\zeta|^2
- \operatorname{Re}\zeta^2)$, $\operatorname{Re}\zeta\operatorname{Im}\zeta =
\tfrac12\operatorname{Im}\zeta^2$ (valid for any $\zeta\in\mathbb C$; expand
$\zeta = x + \mathrm iy$), and $\sum_i|\zeta_i|^2 = \sum_i \ell_i^2 = Q$
(Lemma 1.1(1)):
$$(W^TW)_{11} = \sum_i(\operatorname{Re}\zeta_i)^2 = \tfrac12(Q +
\operatorname{Re}M_4), \quad (W^TW)_{22} = \tfrac12(Q - \operatorname{Re}
M_4), \quad (W^TW)_{12} = \tfrac12\operatorname{Im}M_4 .$$
This matrix is $\tfrac12 Q I_2$ plus $\tfrac12\begin{pmatrix}\alpha & \beta\\
\beta & -\alpha\end{pmatrix}$ with $\alpha = \operatorname{Re}M_4$, $\beta =
\operatorname{Im}M_4$; the traceless part has eigenvalues
$\pm\tfrac12\sqrt{\alpha^2 + \beta^2} = \pm\tfrac12|M_4|$. Hence $\mu_{1,2} =
\tfrac12(Q \pm |M_4|)$ and the sum/difference formulas. $\mu_2 \ge 0$ because
$W^TW \succeq 0$. By Cauchy–Schwarz on $\sum_i \ell_i = 2$ (Lemma 1.1(2)):
$4 = (\sum_i\ell_i)^2 \le nQ$, so $Q \ge 4/n$ and $\mu_1 \ge Q/2 \ge 2/n >
0$. Finally, $WW^T$ and $W^TW$ have the same nonzero eigenvalues (standard),
$K = WW^T \succeq 0$ has rank $\le 2$, and $\mu_1 > 0$, so
$\lambda_{\max}(K) = \mu_1$. $\blacksquare$

---

## 4. Theorem C2

**Theorem 4.1 (C2).** Let $n \ge 2$ and $A \in \mathrm{St}(n,2)$. Set
$v := \ell - \tfrac2n\mathbf 1 \in \mathbb R^n$ and let $\mu_1 \ge \mu_2$ be
the eigenvalues of $W^TW$ (Lemma 3.1). Then the **master identity**
$$\|L - K\|_F^2 - \frac{8Q}{n} + \frac{8}{n^2}
\;=\; 2\Bigl(\mu_2 - \frac2n\Bigr)^2 \;+\; 2s, \qquad
s := \mu_1\|v\|^2 - v^TKv, \tag{$\star$}$$
holds, with $s \ge 0$. Consequently
$$\sum_{i<j} p_{ij}^4 \;\ge\; \frac1n\sum_i \ell_i^2 - \frac1{n^2},
\qquad\text{equivalently}\qquad \sum_{i<j} p_{ij}^2\,b_{ij} \;\ge\; 0,$$
with the exact residual
$$\sum_{i<j} p_{ij}^2 b_{ij} = \sum_{i<j} p_{ij}^4 - \frac Qn + \frac1{n^2}
= \frac14\Bigl(\mu_2 - \frac2n\Bigr)^2 + \frac{s}{4}. \tag{4.1}$$
Equality holds if and only if
$$\mu_2 = \frac2n \quad\text{and}\quad Kv = \mu_1 v \;\;(\text{which includes
} v = 0).$$
In particular: after the single Cauchy–Schwarz step $s \ge 0$, the residual
of the Frobenius-form inequality $\|L-K\|_F^2 \ge 8Q/n - 8/n^2$ is **exactly**
$2(\mu_2 - 2/n)^2$.

*Proof.* By Lemmas 2.2–2.3, the three displayed forms of the inequality are
equivalent, and it suffices to establish $(\star)$ and $s \ge 0$; (4.1) is
$(\star)$ divided by $8$, via (2.3).

**Step 1 (expansion).** Since $\|L\|_F^2 = \|\ell\ell^T\|_F^2 = (\ell^T
\ell)^2 = Q^2$, $\langle L, K\rangle_F = \ell^TK\ell$, and $\|K\|_F^2 =
\operatorname{Tr}(WW^TWW^T) = \operatorname{Tr}\bigl((W^TW)^2\bigr) = \mu_1^2
+ \mu_2^2$ (cyclicity of trace):
$$\|L-K\|_F^2 = Q^2 - 2\,\ell^TK\ell + \mu_1^2 + \mu_2^2 . \tag{4.2}$$

**Step 2 (centering).** $K\mathbf 1 = 0$ (Lemma 1.1(5)) and $K = K^T$ give
$$\ell^TK\ell = (v + \tfrac2n\mathbf 1)^TK(v + \tfrac2n\mathbf 1) = v^TKv,$$
and
$$\|v\|^2 = \ell^T\ell - \tfrac4n\,\mathbf 1^T\ell + \tfrac4{n^2}\,n
= Q - \tfrac8n + \tfrac4n = Q - \tfrac4n \;\ge\; 0,$$
the last inequality being $Q \ge 4/n$ from Lemma 3.1.

**Step 3 (Rayleigh/Cauchy–Schwarz).** $v^TKv = \|W^Tv\|^2 \le
\lambda_{\max}(K)\|v\|^2 = \mu_1\|v\|^2$ by Lemma 3.1, i.e.
$$s = \mu_1\|v\|^2 - v^TKv \ge 0,$$
with equality iff $v$ is an eigenvector of $K$ for the eigenvalue
$\lambda_{\max}(K) = \mu_1$, or $v = 0$ (standard equality case of the
Rayleigh bound for symmetric matrices, applicable since $\mu_1 > 0$).

**Step 4 (algebraic completion).** Substitute $v^TKv = \mu_1\|v\|^2 - s =
\mu_1(Q - \tfrac4n) - s$ into (4.2):
$$\|L-K\|_F^2 - \frac{8Q}n + \frac8{n^2}
= Q^2 + \mu_1^2 + \mu_2^2 - 2\mu_1\Bigl(Q - \frac4n\Bigr) - \frac{8Q}n +
\frac8{n^2} + 2s .$$
Now eliminate $\mu_1 = Q - \mu_2$ (Lemma 3.1) in the non-$s$ part $T$:
$$\begin{aligned}
T &= Q^2 + (Q-\mu_2)^2 + \mu_2^2 - 2(Q - \mu_2)\Bigl(Q - \frac4n\Bigr)
- \frac{8Q}n + \frac8{n^2}\\
&= Q^2 + Q^2 - 2Q\mu_2 + \mu_2^2 + \mu_2^2
- 2Q^2 + \frac{8Q}n + 2Q\mu_2 - \frac{8\mu_2}n - \frac{8Q}n + \frac8{n^2}\\
&= 2\mu_2^2 - \frac{8\mu_2}n + \frac8{n^2}
\;=\; 2\Bigl(\mu_2 - \frac2n\Bigr)^2 .
\end{aligned}$$
This proves $(\star)$; both residual terms are squares/nonnegative, so the
inequality follows, and equality holds iff both vanish: $\mu_2 = 2/n$ and
$s = 0$, the latter being the stated eigenvector condition (Step 3).
$\blacksquare$

**Remark 4.2 (where each hypothesis enters).** Orthonormality is used
exactly twice: through Lemma 1.1(2) ($\sum\zeta_i = 0$, giving $K\mathbf 1 =
0$, i.e. Step 2's centering) and $\sum|\zeta_i| = 2$ (normalizing $\|v\|^2 =
Q - 4/n$ and, via Lemma 2.1, making the weights $p_{ij}^2$ sum to $1$). The
proof never needs the finer structure of $A$ — this is what makes the planar
lemma (§6) possible.

**Proposition 4.3 (equality at the Sengupta–Pautov extremal family).** For
$n \ge 3$ let $A^{(n)} \in \mathrm{St}(n,2)$ have $n-2$ rows $(a, 0)$ and two
rows $(b, c), (b, -c)$, where
$$a^2 = \frac{n-1}{n(n-2)}, \qquad b^2 = \frac1{2n}, \qquad c^2 = \frac12 .$$
Then $A^{(n)} \in \mathrm{St}(n,2)$, and C2 holds **with equality**:
$\mu_2 = 2/n$ exactly, and $v$ is a $\mu_1$-eigenvector of $K$ (for $n = 3$,
$v = 0$).

*Proof.* Orthonormality: $(n-2)a^2 + 2b^2 = \tfrac{n-1}n + \tfrac1n = 1$;
$2c^2 = 1$; $\sum a_ib_i = bc + b(-c) = 0$. Squared rows: for the $n-2$
parallel rows, $z = a$, $\zeta = a^2$ (real); for the conjugate pair, $z = b
\pm \mathrm ic$, $\zeta = (b^2 - c^2) \pm 2\mathrm i\,bc$. Therefore
$$W^TW = \begin{pmatrix} (n-2)a^4 + 2(b^2-c^2)^2 & 0\\ 0 & 8b^2c^2
\end{pmatrix},$$
the off-diagonal entries cancelling between the two conjugate rows. Compute:
$$8b^2c^2 = 8\cdot\frac1{2n}\cdot\frac12 = \frac2n, \qquad
(n-2)a^4 + 2(b^2-c^2)^2 = \frac{(n-1)^2}{n^2(n-2)} +
2\cdot\frac{(n-1)^2}{4n^2} = \frac{(n-1)^2}{2n(n-2)} .$$
Moreover
$$\frac{(n-1)^2}{2n(n-2)} - \frac2n = \frac{(n-1)^2 - 4(n-2)}{2n(n-2)}
= \frac{(n-3)^2}{2n(n-2)} \ge 0,$$
so $\mu_1 = \tfrac{(n-1)^2}{2n(n-2)}$ (the first diagonal entry) and $\mu_2 =
2/n$ exactly.

Eigenvector condition: the $\mu_1$-eigenvector of $W^TW$ is $e_1$, so the
$\mu_1$-eigenspace of $K = WW^T$ contains $We_1 = (\operatorname{Re}
\zeta_i)_i$ (indeed $K(We_1) = W(W^TW)e_1 = \mu_1 We_1$). Entries of $v =
\ell - \tfrac2n\mathbf 1$ vs entries of $We_1$:
$$\text{parallel rows:}\quad v_i = a^2 - \frac2n = \frac{3-n}{n(n-2)}, \qquad
(We_1)_i = a^2 = \frac{n-1}{n(n-2)};$$
$$\text{conjugate rows:}\quad v_i = \frac{n+1}{2n} - \frac2n =
\frac{n-3}{2n}, \qquad (We_1)_i = b^2 - c^2 = \frac{1-n}{2n} .$$
Both ratios equal $\tfrac{3-n}{n-1}$, so $v = \tfrac{3-n}{n-1}\,We_1$ is a
$\mu_1$-eigenvector of $K$ (and $v = 0$ when $n = 3$). By Theorem 4.1,
equality holds in C2. $\blacksquare$

*(Direct cross-check of equality at $n = 4$: $\sum p^4 = 13/64 = Q/4 - 1/16$
with $Q = 17/16$; verified exactly for symbolic $n$ in §8.)*

---

## 5. Corollary A: Perron-free pair existence at $k = 2$

**Corollary 5.1 (unconditional pair existence).** For every $n \ge 2$ and
every $A \in \mathrm{St}(n,2)$ there is a pair $i < j$ with
$$p_{ij} \ne 0 \qquad\text{and}\qquad c_{ij}^2 \;\le\; \Bigl(\ell_i -
\frac1n\Bigr)\Bigl(\ell_j - \frac1n\Bigr),$$
i.e. a pair with $b_{ij} \ge 0$ (by Remark 2.4) and nonsingular $2\times2$
submatrix.

*Proof.* By Theorem 4.1, $\sum_{i<j} p_{ij}^2 b_{ij} \ge 0$. The weights
$p_{ij}^2$ are nonnegative and not all zero, since $\sum_{i<j} p_{ij}^2 = 1$
(Lemma 2.1). Suppose for contradiction that every pair with $p_{ij} \ne 0$
had $b_{ij} < 0$. Then every term $p_{ij}^2 b_{ij}$ is $\le 0$ (terms with
$p_{ij} = 0$ vanish, the others are negative), and at least one term is
strictly negative (some $p_{ij}^2 > 0$), so $\sum_{i<j}p_{ij}^2 b_{ij} < 0$
— contradiction. Hence some pair has $p_{ij} \ne 0$ and $b_{ij} \ge 0$; by
(2.4), $b_{ij} \ge 0 \iff c_{ij}^2 \le (\ell_i - 1/n)(\ell_j - 1/n)$.
$\blacksquare$

**§5.2 Non-circularity.** The claim that Corollary 5.1 is a *new,
Perron-free* proof of the $k=2$ pair-existence step requires checking that
its proof does not (directly or transitively) use that step or comparable
machinery. The complete dependency chain is:

$$\text{(1.1) orthonormality} \;\to\; \text{Lemmas 1.1, 2.1–2.3, 3.1}
\;\to\; \text{Theorem 4.1} \;\to\; \text{Corollary 5.1}.$$

Audit of every input (the list in §0):

- Lemmas 1.1, 2.1–2.3 are direct expansions from (1.1) — no matrix
  positivity theory, no existence statements.
- Lemma 3.1 uses the spectral theorem for one explicit symmetric $2\times2$
  matrix and Cauchy–Schwarz on $(\sum\ell_i)^2 \le nQ$.
- Theorem 4.1 adds one Rayleigh bound $v^TKv \le \mu_1\|v\|^2$ and pure
  algebra.
- **Not used anywhere:** Perron–Frobenius theory (entrywise-positive
  matrices, Perron vectors/roots); induction on $n$; any statement of the
  form "a good pair exists"; the Sengupta–Pautov paper (arXiv:2604.05944) or
  its verification `proofs/sp-verification.md`; the round-1 Pair Lemma
  (`proofs/pair-lemma.md`); GTZ for any $(n,k)$; the Case-A reduction. The
  numerics of §8 are sanity checks, not proof inputs.

The dependency graph is acyclic and grounded in (1.1), so the corollary is a
genuinely independent proof of pair existence. For contrast, the known proof
(S–P Case B heart, as verified in `proofs/sp-verification.md` §6, and its
general-$k$ form `proofs/pair-lemma.md`) proceeds via the matrix $G = WW^T -
zz^T$, a count of its positive eigenvalues, and a Perron–Frobenius
contradiction for the entrywise-positive matrix $M = G + \tfrac2{n^2}E$;
none of those objects appears above. Corollary 5.1 moreover *strengthens*
the conclusion from bare existence to a quantitative statement: the
$p^2$-weighted **average** of the brackets is $\ge 0$, with the exact surplus
(4.1).

**Remark 5.3 (Perron-free core case of GTZ$(n,2)$).** If in addition
$\ell_i + \ell_j \ge 2/n$ holds for the pair of Corollary 5.1 — automatic in
the core case $\min_i \ell_i > 1/n$, which is the only case left after the
Case-A reduction (`proofs/case-a.md`) — then the pair Gram $P_{II}$, $I =
\{i,j\}$, satisfies $\lambda_{\min}(P_{II}) \ge 1/n$, i.e.
$\sigma_{\min}(A_I) \ge 1/\sqrt n$. *Proof:* $\chi_{ij}(1/n) = b_{ij} \ge 0$
means $1/n$ does not lie strictly between the two eigenvalues
$\lambda_{\min} \le \lambda_{\max}$ of $P_{II}$. If $1/n \ge \lambda_{\max}$
then $2/n \ge 2\lambda_{\max} \ge \lambda_{\min} + \lambda_{\max} = \ell_i +
\ell_j \ge 2/n$, forcing $\lambda_{\min} = \lambda_{\max} = 1/n$. In either
case $\lambda_{\min} \ge 1/n$. $\square$
Thus C2 alone settles the core case of GTZ$(n,2)$ Perron-free; the full
GTZ$(n,2)$ still routes through the Case-A induction exactly as in
`proofs/sp-verification.md` — this document does **not** re-prove
GTZ$(n,2)$ from scratch, it replaces the pair-existence engine.

---

## 6. Corollary B: the planar lemma (realizability-free)

**Corollary 6.1 (planar lemma).** Let $n \ge 2$ and let $\zeta_1, \dots,
\zeta_n \in \mathbb C$ satisfy only
$$\sum_i \zeta_i = 0, \qquad \sum_i |\zeta_i| = 2 .$$
Then, with $Q := \sum_i|\zeta_i|^2$, $M_2 := \sum_i \zeta_i|\zeta_i|$,
$M_4 := \sum_i \zeta_i^2$:
$$\Bigl|M_4\Bigr|^2 - 4\Bigl|M_2\Bigr|^2 + 3Q^2 \;\ge\; \frac{16}{n}Q -
\frac{16}{n^2},$$
with equality iff $\mu_2 = 2/n$ and $K v = \mu_1 v$ (or $v = 0$), where
$\ell_i := |\zeta_i|$, $v := \ell - \tfrac2n\mathbf 1$, and $W, K, \mu_1,
\mu_2$ are built from the $\zeta_i$ as in §1/§3. The minimum $0$ of the
slack is attained (image of the family of Proposition 4.3).

*Proof.* Define $\ell_i := |\zeta_i|$, $w_i, W, K, L, Q$ exactly as in §1
(no matrix $A$ needed). Inspect the inputs of Lemma 3.1 and Theorem 4.1:
they use only (i) $\sum_i\zeta_i = 0$ (through $K\mathbf 1 = 0$ and
$W^T\mathbf 1 = 0$), (ii) $\sum_i \ell_i = 2$ (through $\|v\|^2 = Q - 4/n$
and $Q \ge 4/n$), (iii) $|\zeta_i|^2 = \ell_i^2$ (through $K_{ii} =
\ell_i^2$ and $\operatorname{Tr}W^TW = Q$), and linear algebra of $W, K$.
All three hold here by hypothesis/definition. Hence the master identity
$(\star)$ and both slack terms carry over verbatim:
$$\|L - K\|_F^2 \;\ge\; \frac{8Q}n - \frac8{n^2},$$
with the stated equality case. It remains to translate $\|L-K\|_F^2$ into
the moment form:
$$\|L\|_F^2 = Q^2; \qquad
\langle L, K\rangle_F = \ell^TK\ell = \|W^T\ell\|^2 = |M_2|^2,$$
since $W^T\ell = \sum_i \ell_i w_i = (\operatorname{Re}M_2,
\operatorname{Im}M_2)$; and
$$\|K\|_F^2 = \mu_1^2 + \mu_2^2 = \tfrac12\bigl[(\mu_1+\mu_2)^2 +
(\mu_1-\mu_2)^2\bigr] = \tfrac12\bigl(Q^2 + |M_4|^2\bigr)$$
by Lemma 3.1 (whose proof of the closed form of $W^TW$ used none of the
Stiefel structure). Therefore
$$\|L-K\|_F^2 = Q^2 - 2|M_2|^2 + \tfrac12 Q^2 + \tfrac12|M_4|^2
= \tfrac32 Q^2 - 2|M_2|^2 + \tfrac12|M_4|^2,$$
and multiplying the inequality by $2$ gives the statement. Equality
attainment: the squared rows of Proposition 4.3 satisfy the two constraints
and give equality. $\blacksquare$

**Remark 6.2 (relation to Theorem C2; surjectivity of the squaring map).**
For $A \in \mathrm{St}(n,2)$, $16\sum_{i<j}p_{ij}^4 = 2\|L-K\|_F^2 = 3Q^2 -
4|M_2|^2 + |M_4|^2$ (Lemma 2.3 + the translation above), so Corollary 6.1
applied to $\zeta_i = z_i^2$ **is** Theorem C2. Conversely the squaring map
is *onto* the constraint set: given $\zeta_i$ with $\sum\zeta_i = 0$,
$\sum|\zeta_i| = 2$, choose any square roots $z_i$ ($z_i^2 = \zeta_i$) and
set $a = \operatorname{Re}z$, $b = \operatorname{Im}z$; then $\sum z_i^2 = 0$
gives $\|a\|^2 = \|b\|^2$ and $\langle a,b\rangle = 0$, while $\sum|z_i|^2 =
2$ gives $\|a\|^2 + \|b\|^2 = 2$; hence $A = [a\;b] \in \mathrm{St}(n,2)$.
So the two statements are *equivalent* — but the planar formulation shows
that Stiefel realizability plays no role: the $k=2$ average inequality is a
fact about centered planar vector families with total variation $2$.

---

## 7. How the $p^2$-weights repair the failed uniform average

**Proposition 7.1 (the uniform average fails).** For every $A \in
\mathrm{St}(n,2)$, $n \ge 2$:
$$\sum_{i<j} b_{ij} \;=\; -\frac12 + \frac3{2n},$$
which is $< 0$ for every $n \ge 4$ (and $= 0$ at $n = 3$). Hence for $n \ge
4$ the *uniform* average of the brackets can never certify a good pair.

*Proof.* By Lemma 2.1, $\sum_{i<j}p_{ij}^2 = 1$. Also $\sum_{i<j}(\ell_i +
\ell_j) = (n-1)\sum_i\ell_i = 2(n-1)$ (each $\ell_i$ appears in $n-1$
pairs), and there are $\binom n2$ pairs. So
$$\sum_{i<j}b_{ij} = 1 - \frac{2(n-1)}n + \frac{n(n-1)}{2n^2}
= 1 - 2 + \frac2n + \frac12 - \frac1{2n} = -\frac12 + \frac3{2n}.
\qquad\blacksquare$$
This is the failed average recorded in MEMORY.md; it shows a *selection*
(or reweighting) argument is mandatory at $k = 2$.

**Proposition 7.2 (anatomy of the failure at the extremal family).** For
the family $A^{(n)}$ of Proposition 4.3 ($n \ge 3$), the three pair types
have exactly:

| pair type | count | $p_{ij}^2$ | $b_{ij}$ |
|---|---|---|---|
| parallel–parallel $\{(a,0),(a,0)\}$ | $\binom{n-2}2$ | $0$ | $-\dfrac1{n(n-2)}$ |
| mixed $\{(a,0),(b,\pm c)\}$ | $2(n-2)$ | $a^2c^2 = \dfrac{n-1}{2n(n-2)}$ | $0$ |
| conjugate $\{(b,c),(b,-c)\}$ | $1$ | $4b^2c^2 = \dfrac1n$ | $0$ |

Consequently the entire uniform-average deficit is carried by the
zero-Plücker pairs:
$$\sum_{\text{parallel pairs}} b_{ij} = \binom{n-2}2\cdot
\Bigl(-\frac1{n(n-2)}\Bigr) = -\frac{n-3}{2n} = -\frac12 + \frac3{2n}
= \sum_{i<j} b_{ij},$$
while every pair with $p_{ij} \ne 0$ has $b_{ij} = 0$ exactly.

*Proof.* Parallel–parallel: both rows are multiples of $(1,0)$, so $p_{ij}
= a\cdot 0 - 0\cdot a = 0$ and, by (2.1), $b = 0 - \tfrac{2a^2}n +
\tfrac1{n^2}$; with $a^2 = \tfrac{n-1}{n(n-2)}$,
$$b = \frac1{n^2} - \frac{2(n-1)}{n^2(n-2)} = \frac{(n-2) - 2(n-1)}
{n^2(n-2)} = \frac{-n}{n^2(n-2)} = -\frac1{n(n-2)} .$$
Mixed: $p = a(\pm c) - 0\cdot b = \pm ac$, $p^2 = a^2c^2$, and with $\ell_i
= a^2$, $\ell_j = b^2 + c^2 = \tfrac{n+1}{2n}$:
$$b = \frac{n-1}{2n(n-2)} - \frac1n\Bigl(\frac{n-1}{n(n-2)} +
\frac{n+1}{2n}\Bigr) + \frac1{n^2}
= \frac{n(n-1) - (n^2+n-4)}{2n^2(n-2)} + \frac1{n^2}
= \frac{-(2n-4)}{2n^2(n-2)} + \frac1{n^2} = 0 .$$
Conjugate: $p = b(-c) - cb = -2bc$, $p^2 = 4b^2c^2 = \tfrac1n$, $\ell_i =
\ell_j = \tfrac{n+1}{2n}$:
$$b = \frac1n - \frac1n\cdot\frac{n+1}n + \frac1{n^2} = \frac1n -
\frac{n+1}{n^2} + \frac1{n^2} = 0 .$$
The deficit sum: $\binom{n-2}2/\bigl(n(n-2)\bigr) = \tfrac{(n-2)(n-3)}
{2n(n-2)} = \tfrac{n-3}{2n}$, and $-\tfrac{n-3}{2n} = -\tfrac12 +
\tfrac3{2n}$, which equals $\sum_{i<j}b_{ij}$ by Proposition 7.1 —
consistent, since all remaining brackets vanish. $\blacksquare$

**The repair, precisely.** Propositions 7.1–7.2 exhibit the mechanism:

1. At the extremal configuration, every negative bracket sits on a pair with
   $p_{ij} = 0$ (linearly dependent rows), and every pair with $p_{ij} \ne
   0$ has bracket exactly $0$ (these are the *active* pairs: by Remark 2.4,
   $b_{ij} = \chi_{ij}(1/n) = 0$ means $1/n$ is an eigenvalue of the pair
   Gram — the pair sits exactly at the GTZ threshold).
2. The uniform average charges the $\asymp n^2/2$ dependent pairs their
   $-1/(n(n-2)) \asymp -1/n^2$ each, accumulating the fatal deficit
   $-\tfrac12 + \tfrac3{2n}$; the weights $p_{ij}^2$ annihilate exactly
   those pairs and redistribute all mass (total $1$, Lemma 2.1) onto the
   threshold pairs, where the average becomes $0$ — the equality case of
   Theorem 4.1.
3. Theorem C2 asserts this reweighting works not just at the extremal family
   but for **every** $A \in \mathrm{St}(n,2)$: the weighted average is
   globally nonnegative, with residual (4.1) measuring the distance from the
   extremal spectral data ($\mu_2 = 2/n$, $v$ aligned with the top of $K$).

**Remark 7.3 (expected-charpoly form; why the scope guard matters).** By
Remark 2.4 and the sums in Lemma 2.2, the $p^2$-weighted average is the
value at $1/n$ of the volume-weighted expected characteristic polynomial
$$E(x) := \sum_{i<j} p_{ij}^2\,\chi_{ij}(x) = x^2 - Qx + \sum_{i<j}
p_{ij}^4, \qquad E(1/n) = \sum_{i<j}p_{ij}^2 b_{ij} \ge 0 \;(\text{C2}).$$
Since the vertex of the monic parabola $E$ lies at $Q/2 \ge 2/n > 1/n$
(Lemma 3.1), $E(1/n) \ge 0$ means $E$ has **no real root below $1/n$**: if a
root $\lambda_- < 1/n$ existed, $1/n$ would lie strictly between the two
roots (as $1/n <$ vertex $\le \lambda_+$), forcing $E(1/n) < 0$. So C2 is
exactly the $k = 2$ case of "$\mathrm{minroot}(E) \ge 1/n$". The $k \ge 3$
analog of this statement over $k$-subsets is **FALSE** (certified core
counterexamples at $(10,5)$, $(12,4)$: obstruction O2, ledger §6) — which
is why this document claims nothing beyond $k = 2$.

**Remark 7.4 (canonicity of the weights: the apolar invariant; aside, used
nowhere).** The certificate $\sum_{i<j}p_{ij}^4$ is not an ad-hoc choice: it
is the classical quadratic invariant of the empirical binary quartic of the
rows. Let $\Phi(u) = \sum_i \langle r_i, u\rangle^4$, a binary quartic with
coefficients $m_{jk} = \sum_i a_i^jb_i^k$ ($j + k = 4$) in the expansion
$\Phi = m_{40}u_1^4 + 4m_{31}u_1^3u_2 + 6m_{22}u_1^2u_2^2 + 4m_{13}u_1u_2^3
+ m_{04}u_2^4$, and let $I(q) = m_{40}m_{04} - 4m_{31}m_{13} + 3m_{22}^2$
be its apolar (Hankel) invariant. Then
$$I(\Phi) = \sum_{i<j} p_{ij}^4 .$$
*Proof.* $I$ is a quadratic form in the coefficients; let $B$ be its
polarization, $B(q, q') = \tfrac12(m_{40}m'_{04} + m'_{40}m_{04}) -
2(m_{31}m'_{13} + m'_{31}m_{13}) + 3m_{22}m'_{22}$, so that $I(\sum_i q_i)
= \sum_{i,j} B(q_i, q_j)$. For pure fourth powers $q = \langle
(a_1,b_1),u\rangle^4$, $q' = \langle (a_2,b_2),u\rangle^4$ the coefficients
are $m_{jk} = a_1^jb_1^k$, $m'_{jk} = a_2^jb_2^k$; substituting and setting
$x = a_1b_2$, $y = b_1a_2$:
$$B(q,q') = \tfrac12(x^4 + y^4) - 2(x^3y + xy^3) + 3x^2y^2 =
\tfrac12(x-y)^4 = \tfrac12\,p_{12}^4 .$$
In particular $B(q,q) = 0$ (fourth powers are $I$-null), so $I(\Phi) =
\sum_{i \ne j}\tfrac12 p_{ij}^4 = \sum_{i<j}p_{ij}^4$. $\square$
Thus $\sum p^4$ is (up to scale) *the* fundamental $SL_2$-invariant of
degree $2$ in the degree-$4$ row data — the $p^2$-weighting is the canonical
one, which is why it, and not the uniform weighting, is the right average.

---

## 8. Numeric sanity block

Fresh, deterministic study for this document: `numerics/study-theorem-c2/`
(script `check_c2.py`, numpy seed 20260718; the sympy blocks are exact with
symbolic $n$). Command and observed result:

```
$ cd numerics/study-theorem-c2 && ../../.venv/bin/python check_c2.py ; echo "exit=$?"
```

Observed results (run 2026-07-18, exit 0, full log in
`numerics/study-theorem-c2/check_c2.log`; all checks PASS):

- **A. 6000 Haar samples of $\mathrm{St}(n,2)$, $n \in [2, 41)$** — every
  identity of §§1–4 verified at machine precision, and both inequality
  slacks nonnegative:
  - A1 constraints $\sum\zeta = 0$, $\sum|\zeta| = 2$: max error `1.8e-15`;
  - A2 dictionary $p_{ij}^2 = \tfrac12(L-K)_{ij}$ entrywise: `4.4e-16`;
  - A3 row marginal $\sum_{j\ne i}p_{ij}^2 = \ell_i$ and $\sum_{i<j}p^2 =
    1$: `1.4e-15`;
  - A4 bracket two faces (2.1) $=$ (2.4): `8.9e-16`;
  - A5 weighted identity (2.2) and Frobenius form (2.3): `8.9e-16`;
  - A6 spectral identities ($K\mathbf 1 = 0$, $\mu_1 + \mu_2 = Q$, $\mu_1 -
    \mu_2 = |M_4|$): `1.8e-15`;
  - A7 master identity $(\star)$ residual: `7.1e-15`; Cauchy–Schwarz slack
    $s$ never negative;
  - **A8 C2 slack $\sum p^4 - Q/n + 1/n^2 \ge 0$: min `+7.40e-04` over all
    6000 samples, 0 violations**;
  - A9 Corollary 5.1 witness: in every sample, some pair with $p_{ij}^2 >
    10^{-12}$ has $b_{ij} \ge 0$ (min over samples of the best bracket:
    `+6.61e-03`);
  - A10 planar translation $16\sum p^4 = 3Q^2 - 4|M_2|^2 + |M_4|^2$:
    `2.1e-14`;
  - A11 core-case check (Remark 5.3): among the 303 samples with
    $\min_i\ell_i > 1/n$, the best-bracket pair always has
    $\lambda_{\min}(P_{II}) \ge 1/n$ (min `n·λ_min = 1.1071044`).
- **B. Exact (sympy, symbolic $n$), Proposition 4.3 / §7:** all of the
  following simplify to $0$ identically in $n$: $\mu_2 - 2/n$;
  $\mu_1 - \mu_2 - \tfrac{(n-3)^2}{2n(n-2)}$; the C2 equality residual
  $\sum p^4 - Q/n + 1/n^2$ at the family; $(W^TW)_{12}$; the collinearity of
  $v$ and $We_1$; $b_{\text{mixed}}$, $b_{\text{conj}}$, and
  $b_{\text{par}} + \tfrac1{n(n-2)}$; the uniform-average identity
  $\binom{n-2}2 b_{\text{par}} - (-\tfrac12 + \tfrac3{2n})$; the mixed-pair
  "two faces" identity. **PASS (all exact zeros)**.
- **C. Planar lemma (Corollary 6.1):** 4000 random centered normalized
  planar configs, $n \in [2,40)$: min slack `+2.61e-03`, no violation; slack
  at the S–P image, $n \in \{4, 7, 12, 30\}$: max $|{\cdot}|$ `3.9e-16`
  (equality attained); adversarial Nelder–Mead minimization ($n \in \{5, 9,
  16\}$, 8 seeded starts each): minimum found `-1.4e-15` $\approx 0$
  (machine-precision zero, i.e. the adversarial minimum is $0$, attained).

Reproducibility of the promoted claims: the original round-3 script also
still passes in the current workspace,

```
$ cd numerics/study-attack-r3-variational && ../../.venv/bin/python s3_c2_proof_check.py ; echo "exit=$?"   # -> ALL PASS, exit=0
```

(Rerun 2026-07-18: ALL PASS, exit 0 — 5000 Haar samples $n \in [3,30)$ at
$\le 1.8\cdot10^{-15}$, S–P equality family slack $\le 6.9\cdot10^{-17}$ with
$\mu_2 - 2/n \le 2.3\cdot10^{-16}$ at $n \in \{4,6,10,25\}$, planar
adversarial minimum $-1.8\cdot10^{-15} \approx 0$ — matching the round-3 log
`s3.log`.)

---

## 9. Deliverables and status lines

- **Theorem 4.1 (C2):** $\sum_{i<j}p_{ij}^4 \ge \tfrac1n\sum_i\ell_i^2 -
  \tfrac1{n^2}$ for every $A \in \mathrm{St}(n,2)$, $n \ge 2$; equivalently
  the $p^2$-weighted average of the S–P brackets is nonnegative; master
  identity $(\star)$ with post–Cauchy–Schwarz residual exactly $2(\mu_2 -
  2/n)^2$; equality iff $\mu_2 = 2/n$ and $v = \ell - \tfrac2n\mathbf 1$ is
  a top eigenvector of $K$ (or $0$). **Status: PROVED (unconditional), no
  gaps.**
- **Proposition 4.3:** equality holds at the S–P extremal family, with
  $\mu_2 = 2/n$ exact. **Status: PROVED.**
- **Corollary 5.1 + §5.2:** Perron-free, non-circular proof of the
  unconditional $k=2$ pair-existence step; Remark 5.3: Perron-free core-case
  bound $\sigma_{\min}(A_I) \ge 1/\sqrt n$. **Status: PROVED
  (unconditional), no gaps.**
- **Corollary 6.1:** realizability-free planar lemma, equivalent to C2 via
  the (onto) squaring map. **Status: PROVED.**
- **§7:** exact computation of the failed uniform average $-\tfrac12 +
  \tfrac3{2n}$ and its anatomy: at the extremal family the entire deficit is
  carried by zero-Plücker pairs, which the $p^2$-weights annihilate.
  **Status: PROVED.**
- **Not claimed:** any $k \ge 3$ statement. C2k remains CONJECTURED in
  `attacks/round-3-variational.md` §7; volume-weighted $k$-subset analogs
  are FALSE for $k \ge 3$ (obstructions O1/O2, ibid. §6).
