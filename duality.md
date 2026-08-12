# Duality for the GTZ hypothesis: $\mathrm{GTZ}(n,k) \Leftrightarrow \mathrm{GTZ}(n,n-k)$, and the case $k = n-2$

**Deliverables.** Theorem 1 + Corollary 1 (target VAL-DUAL-001, spectral duality and the
equivalence $\mathrm{GTZ}(n,k) \Leftrightarrow \mathrm{GTZ}(n,n-k)$ with the same threshold
$1/\sqrt{n}$) and Corollary 2 (target VAL-DUAL-002, $\mathrm{GTZ}(n,n-2)$ for all $n \ge 4$,
conditional on arXiv:2604.05944).

**Status summary.**

| Statement | Status |
|---|---|
| Lemma 1 (equal nonzero Gram spectra) | PROVED |
| Lemma 2 (submatrix singular values $\le 1$) | PROVED |
| Lemma 3 (three Gram identities) | PROVED |
| Theorem 1 (spectral duality identity, both orientations) | PROVED |
| Corollary 1 ($\sigma_{\min}$ transfer; $f(A)=f(B)$; $\mathrm{GTZ}(n,k)\Leftrightarrow\mathrm{GTZ}(n,n-k)$) | PROVED |
| Corollary 2 ($\mathrm{GTZ}(n,n-2)$ for all $n\ge 4$) | PROVED CONDITIONAL ON arXiv:2604.05944 (Hypothesis 1 for $k=2$, proved in its Section 2) |

No GAP markers appear in this document.

---

## 0. Notation and conventions

We use the shared project notation (MEMORY.md). Throughout, $n \ge 2$ and $1 \le k \le n-1$
are integers, and $[n] = \{1,\dots,n\}$.

- For a matrix $M$ with $n$ rows and a subset $S \subseteq [n]$, $M_S$ denotes the
  $|S| \times (\#\text{cols})$ submatrix formed by the rows of $M$ indexed by $S$, in
  increasing index order. The choice of row order is immaterial for everything below:
  permuting rows replaces $M_S$ by $\Pi M_S$ with $\Pi$ a permutation (hence orthogonal)
  matrix, which changes neither $M_S^T M_S$ nor the singular values.
- $\sigma(M)$ denotes the **multiset** of singular values of a square matrix $M \in
  \mathbb{R}^{m\times m}$: the $m$ nonnegative square roots of the eigenvalues (with
  multiplicity) of the symmetric positive-semidefinite matrix $M^T M$. $\sigma_{\min}$,
  $\sigma_{\max}$ are its smallest and largest elements. No invertibility is assumed
  anywhere; $\sigma_{\min}(M) = 0$ exactly when $M$ is singular.
- $\mathrm{spec}(X)$ is the multiset of eigenvalues of a symmetric matrix $X$.
- Multiset conventions: $\cup$ is multiset union (multiplicities add); $\{1\}^m$ is the
  multiset of $m$ ones ($\{1\}^0 = \emptyset$); a multiset identity "$\mathcal{S} =
  \mathcal{T}$" asserts equality of elements **with multiplicities**. For $c \in \mathbb{R}$
  and a multiset $\mathcal{S}$, $c - \mathcal{S} := \{c - s : s \in \mathcal{S}\}$ (with
  multiplicities). Since $x \mapsto \sqrt{x}$ is a bijection of $[0,\infty)$, two multisets
  of nonnegative reals are equal iff their elementwise square-root multisets are equal.
- $\mathrm{GTZ}(n,k)$, for $1 \le k \le n-1$, denotes the statement: *for every
  $A \in \mathbb{R}^{n\times k}$ with $A^T A = I_k$ there exists $I \subseteq [n]$,
  $|I| = k$, with $\sigma_{\min}(A_I) \ge 1/\sqrt{n}$.* Equivalently (Formulations 1–4 of
  arXiv:2303.07492), writing $f(A) := \max_{|I|=k} \sigma_{\min}(A_I)$:
  $\mathrm{GTZ}(n,k)$ says $f(A) \ge 1/\sqrt{n}$ for all such $A$.
- **Standing setup for §§1–2:** $Q = [A \mid B] \in O(n)$, i.e. $Q^T Q = I_n$, where
  $A \in \mathbb{R}^{n\times k}$ comprises the first $k$ columns and
  $B \in \mathbb{R}^{n\times(n-k)}$ the last $n-k$ columns. Then automatically the columns
  of $A$ are orthonormal, the columns of $B$ are orthonormal, and every column of $B$ is
  orthogonal to every column of $A$; since $\dim(\operatorname{span} A)^\perp = n - k$ and
  the $n-k$ orthonormal columns of $B$ lie in $(\operatorname{span} A)^\perp$, they form an
  orthonormal basis of $(\operatorname{span} A)^\perp$. Conversely, if $A$ has orthonormal
  columns and the columns of $B$ are an orthonormal basis of $(\operatorname{span} A)^\perp$,
  then $[A\mid B] \in O(n)$. So the two descriptions in the assignment coincide.
- We exclude the boundary values $k \in \{0, n\}$ from the theorems: they involve empty
  matrices ($0\times 0$ blocks), for which $\sigma_{\min}$ requires an arbitrary convention.
  For the record, $\mathrm{GTZ}(n,n)$ is trivially true ($A_{[n]} = A \in O(n)$ has
  $\sigma_{\min} = 1 \ge 1/\sqrt n$), and the range $1 \le k \le n-1$ is closed under
  $k \mapsto n-k$, which is all the duality needs.

---

## 1. Three lemmas

### Lemma 1 (equal nonzero spectra of the two Gram matrices)

*Let $X \in \mathbb{R}^{p\times q}$. Then the multisets $\mathrm{spec}(XX^T)$ (size $p$) and
$\mathrm{spec}(X^T X)$ (size $q$) have identical nonzero parts, including multiplicities;
consequently they differ exactly by $|p-q|$ extra zeros on the larger side, and more
precisely the multiplicity of the eigenvalue $0$ is $p - \operatorname{rank} X$ in
$\mathrm{spec}(XX^T)$ and $q - \operatorname{rank} X$ in $\mathrm{spec}(X^T X)$.*

**Proof.** Both $XX^T$ and $X^TX$ are symmetric positive semidefinite, so each is
diagonalizable and every eigenvalue's algebraic multiplicity equals the dimension of its
eigenspace. Fix $\lambda \ne 0$ and let $E_\lambda(X^TX) = \ker(X^TX - \lambda I_q)$,
$E_\lambda(XX^T) = \ker(XX^T - \lambda I_p)$. The linear map $v \mapsto Xv$ sends
$E_\lambda(X^TX)$ into $E_\lambda(XX^T)$: if $X^TXv = \lambda v$ then
$XX^T(Xv) = X(X^TXv) = \lambda\, Xv$. It is injective on $E_\lambda(X^TX)$: if $Xv = 0$ for
$v \in E_\lambda(X^TX)$ then $\lambda v = X^TXv = 0$, so $v = 0$ (as $\lambda \ne 0$).
Hence $\dim E_\lambda(X^TX) \le \dim E_\lambda(XX^T)$. Applying the same argument to $X^T$
in place of $X$ gives the reverse inequality, so the dimensions are equal. Since this holds
for every $\lambda \neq 0$ and multiplicities are eigenspace dimensions, the nonzero parts
of the two spectra coincide as multisets. The zero-eigenvalue multiplicities are
$\dim\ker(XX^T) = p - \operatorname{rank}(XX^T) = p - \operatorname{rank} X$ and likewise
$q - \operatorname{rank} X$, using $\operatorname{rank}(XX^T) = \operatorname{rank}(X^TX)
= \operatorname{rank} X$ (from $\ker(X^TX) = \ker X$: $X^TXv=0 \Rightarrow \|Xv\|^2 =
v^TX^TXv = 0$). $\blacksquare$

### Lemma 2 (row-submatrices of a column-orthonormal matrix are contractions)

*Let $M \in \mathbb{R}^{n\times m}$ satisfy $M^T M = I_m$, and let $S \subseteq [n]$. Then
$\sigma_{\max}(M_S) \le 1$, where for a (possibly rectangular) matrix
$\sigma_{\max}(X) := \max_{\|x\|=1}\|Xx\|$ is the spectral norm — for square $X$ this equals
$\max \sigma(X)$, since $\|Xx\|^2 = x^TX^TXx$ is maximized over unit vectors at
$\lambda_{\max}(X^TX)$. Equivalently $0 \preceq M_S^T M_S \preceq I_m$, i.e. all
eigenvalues of $M_S^T M_S$ lie in $[0,1]$.*

**Proof.** Write the rows of $M$ as $r_1^T,\dots,r_n^T$ ($r_i \in \mathbb{R}^m$). For any
$x \in \mathbb{R}^m$,
$$\|M_S x\|^2 = \sum_{i \in S} (r_i^T x)^2 \le \sum_{i=1}^{n} (r_i^T x)^2 = \|Mx\|^2
 = x^T M^T M x = \|x\|^2 ,$$
since the discarded terms $(r_i^Tx)^2$, $i \notin S$, are nonnegative. Hence
$\sigma_{\max}(M_S)^2 = \max_{\|x\|=1} \|M_S x\|^2 \le 1$. The matrix $M_S^T M_S$ is
positive semidefinite by construction, and $x^T(I_m - M_S^TM_S)x = \|x\|^2 - \|M_Sx\|^2
\ge 0$ gives the upper bound. $\blacksquare$

### Lemma 3 (the three Gram identities)

*In the standing setup $Q = [A\mid B] \in O(n)$, for every $I \subseteq [n]$ with
complement $I^c = [n]\setminus I$:*

$$\text{(C1)}\quad A_I^T A_I + A_{I^c}^T A_{I^c} = I_k, \qquad
  \text{(C2)}\quad B_I^T B_I + B_{I^c}^T B_{I^c} = I_{n-k}, \qquad
  \text{(R)}\quad A_I A_I^T + B_I B_I^T = I_{|I|}.$$

**Proof.** Write the rows of $A$ as $a_1^T,\dots,a_n^T$ ($a_i \in \mathbb{R}^k$). Then
$A^T A = \sum_{i=1}^n a_i a_i^T$, and splitting the sum over $I$ and $I^c$ gives
$A^TA = A_I^TA_I + A_{I^c}^TA_{I^c}$; with $A^TA = I_k$ this is (C1). Identity (C2) is the
same computation for $B$, using $B^TB = I_{n-k}$.

For (R): $Q^TQ = I_n$ means $Q$ is invertible with $Q^{-1} = Q^T$; a two-sided inverse in
$\mathbb{R}^{n\times n}$ gives also $QQ^T = I_n$, i.e. the rows of $Q$ are orthonormal.
The $(i,j)$ entry of $QQ^T$ for $i,j \in I$ reads $a_i^Ta_j + b_i^Tb_j = \delta_{ij}$
(with $b_i^T$ the rows of $B$), and collecting these entries for $i,j \in I$ yields
exactly $Q_I Q_I^T = A_IA_I^T + B_IB_I^T = I_{|I|}$. $\blacksquare$

---

## 2. The duality theorem

### Theorem 1 (spectral duality)

*Let $n \ge 2$, $1 \le k \le n-1$, $Q = [A \mid B] \in O(n)$ as in §0, and let
$I \subseteq [n]$ with $|I| = k$, so that $A_I$ is $k\times k$ and $B_{I^c}$ is
$(n-k)\times(n-k)$. Then, as multisets:*

- *(i) if $n \ge 2k$:* $\;\sigma(B_{I^c}) = \sigma(A_I) \cup \{1\}^{\,n-2k}$;
- *(ii) if $n < 2k$:* $\;\sigma(A_I) = \sigma(B_{I^c}) \cup \{1\}^{\,2k-n}$.

*Equivalently, in both orientations: the singular spectrum of the wider of the two square
matrices $A_I,\ B_{I^c}$ equals that of the narrower one together with $|n-2k|$ additional
singular values equal to $1$; for $n = 2k$ the two spectra coincide exactly.*

**Proof.** All matrices $A_I A_I^T$, $A_I^TA_I$, $B_IB_I^T$, $B_I^TB_I$,
$B_{I^c}^TB_{I^c}$ below are symmetric positive semidefinite, so their spectra are
multisets of nonnegative reals of sizes equal to their dimensions.

**Step 1 (spectrum of the $k\times k$ Gram matrix $B_IB_I^T$).** $A_I$ is square
($k \times k$), so $\mathrm{spec}(A_IA_I^T) = \mathrm{spec}(A_I^TA_I) =
\{\sigma_i(A_I)^2\}_{i=1}^k =: \mathcal{S}$ — indeed by Lemma 1 the two spectra share their
nonzero parts and both have size $k$, hence they have equal zero multiplicities too, and
$\mathrm{spec}(A_I^TA_I)$ is the definition of the squared singular values. By identity
(R) of Lemma 3 (here $|I| = k$),
$$B_IB_I^T = I_k - A_IA_I^T,$$
so
$$\mathcal{T} := \mathrm{spec}(B_IB_I^T) = 1 - \mathcal{S}
 = \{\,1 - \sigma_i(A_I)^2\,\}_{i=1}^k$$
(a size-$k$ multiset; here we use that for a symmetric $X$ with spectral decomposition
$X = \sum_\lambda \lambda\, \Pi_\lambda$, one has $I - X = \sum_\lambda (1-\lambda)\Pi_\lambda$,
so $\mathrm{spec}(I - X) = 1 - \mathrm{spec}(X)$ with multiplicities). By Lemma 2 applied
to $A$, $\mathcal{S} \subseteq [0,1]$, hence also $\mathcal{T} \subseteq [0,1]$.

**Step 2 (transfer to the $(n-k)\times(n-k)$ Gram matrix $B_I^TB_I$).** Apply Lemma 1 to
$X = B_I \in \mathbb{R}^{k\times(n-k)}$. Let $\mathcal{T}^{+}$ be the nonzero part of
$\mathcal{T}$ and $m := |\mathcal{T}^{+}| = \operatorname{rank} B_I$. Note
$m \le \min(k,\,n-k)$, since the rank of a $k\times(n-k)$ matrix is at most either
dimension. Lemma 1 gives, as multisets,
$$\mathcal{T} = \mathrm{spec}(B_IB_I^T) = \mathcal{T}^{+} \cup \{0\}^{\,k-m},
 \qquad
 \mathrm{spec}(B_I^TB_I) = \mathcal{T}^{+} \cup \{0\}^{\,n-k-m}.$$

**Step 3 (spectrum of $B_{I^c}^TB_{I^c}$).** By identity (C2) of Lemma 3,
$B_{I^c}^TB_{I^c} = I_{n-k} - B_I^TB_I$, and $B_{I^c}$ is square of size $n-k$, so — by the
definition of $\sigma$ via the Gram matrix $B_{I^c}^TB_{I^c}$ (§0) and again
$\mathrm{spec}(I - X) = 1 - \mathrm{spec}(X)$ —
$$\sigma(B_{I^c})^2 := \{\sigma_j(B_{I^c})^2\}_{j=1}^{n-k}
 = \mathrm{spec}(B_{I^c}^TB_{I^c})
 = 1 - \mathrm{spec}(B_I^TB_I)
 = (1 - \mathcal{T}^{+}) \cup \{1\}^{\,n-k-m}.$$

**Step 4 (case $n \ge 2k$).** Here $n - k - m \ge k - m \ge 0$, so
$$\mathrm{spec}(B_I^TB_I) = \mathcal{T}^{+} \cup \{0\}^{\,k-m} \cup \{0\}^{\,n-2k}
 = \mathcal{T} \cup \{0\}^{\,n-2k},$$
and therefore, from Step 3 and Step 1,
$$\sigma(B_{I^c})^2 = (1-\mathcal{T}) \cup \{1\}^{\,n-2k}
 = \mathcal{S} \cup \{1\}^{\,n-2k}
 = \{\sigma_i(A_I)^2\}_{i=1}^k \cup \{1\}^{\,n-2k}.$$
Both sides are multisets of nonnegative reals of size $n-k$; taking elementwise square
roots (a multiset bijection on $[0,\infty)$, cf. §0) yields
$\sigma(B_{I^c}) = \sigma(A_I) \cup \{1\}^{\,n-2k}$, which is (i).

**Step 5 (case $n < 2k$).** Here $m \le n-k < k$, so $k - m \ge (n-k-m) + (2k-n)$ with
$2k - n > 0$, and
$$\mathcal{T} = \mathcal{T}^{+} \cup \{0\}^{\,k-m}
 = \left(\mathcal{T}^{+} \cup \{0\}^{\,n-k-m}\right) \cup \{0\}^{\,2k-n}
 = \mathrm{spec}(B_I^TB_I) \cup \{0\}^{\,2k-n}.$$
Hence, using Step 1 and Step 3,
$$\{\sigma_i(A_I)^2\}_{i=1}^k = \mathcal{S} = 1 - \mathcal{T}
 = \left(1 - \mathrm{spec}(B_I^TB_I)\right) \cup \{1\}^{\,2k-n}
 = \sigma(B_{I^c})^2 \cup \{1\}^{\,2k-n},$$
and taking square roots gives $\sigma(A_I) = \sigma(B_{I^c}) \cup \{1\}^{\,2k-n}$, which is
(ii). (In particular, when $n < 2k$ at least $2k - n$ singular values of $A_I$ equal $1$
automatically, for **every** $I$ — the padding lands on whichever side is wider.)

For $n = 2k$ both cases give $\sigma(B_{I^c}) = \sigma(A_I)$ with no padding.
$\blacksquare$

**Status: PROVED.**

### Remarks on degenerate and boundary cases

1. **$k = n/2$ (self-dual size, $n$ even).** $n - 2k = 0$; Theorem 1 reads
   $\sigma(B_{I^c}) = \sigma(A_I)$ exactly, with no padding — both matrices are
   $k \times k$. Nothing in the proof changes; Steps 4 and 5 coincide.
2. **Singular $A_I$ / rank drops.** No step of the proof divides by anything or inverts
   anything: only Gram identities and spectra of positive-semidefinite matrices are used.
   Counting zeros in the two multisets of Theorem 1 (the padded ones are nonzero) gives the
   **corank identity**
   $$k - \operatorname{rank} A_I = (n-k) - \operatorname{rank} B_{I^c},$$
   valid in both orientations. In particular $A_I$ is singular iff $B_{I^c}$ is singular,
   and $\sigma_{\min}(A_I) = 0 \Leftrightarrow \sigma_{\min}(B_{I^c}) = 0$.
3. **Repeated singular values.** All spectral statements are multiset statements; Lemma 1
   was proved with multiplicities (via eigenspace dimensions, using diagonalizability of
   symmetric matrices, so algebraic = geometric multiplicity), and the maps
   $\mathcal{X} \mapsto 1 - \mathcal{X}$ and $\mathcal{X} \mapsto \sqrt{\mathcal{X}}$
   preserve multiplicities. No genericity or distinctness is assumed anywhere.
4. **Extreme admissible $k$.** $k = 1$ and $k = n-1$ are covered by the theorem
   ($1 \le k \le n-1$); e.g. for $k = 1$, $A_I$ is the $1\times1$ matrix $(a_i)$ and
   $\sigma(B_{I^c})$ consists of $|a_i|$ together with $n-2$ ones.
5. **Row order.** As noted in §0, singular values of $A_I$, $B_{I^c}$ do not depend on the
   ordering of the selected rows.

### Corollary 1 ($\sigma_{\min}$ transfer, $f$-invariance, and the GTZ equivalence)

*(a) In the setting of Theorem 1, for every $I$ with $|I| = k$ ($1 \le k \le n-1$):*
$$\sigma_{\min}(A_I) = \sigma_{\min}(B_{I^c}).$$

*(b) Consequently $f(A) = f(B)$, where $f(A) = \max_{|I|=k}\sigma_{\min}(A_I)$ and
$f(B) = \max_{|J|=n-k}\sigma_{\min}(B_J)$.*

*(c) For all $n \ge 2$ and $1 \le k \le n-1$: $\mathrm{GTZ}(n,k)$ holds if and only if
$\mathrm{GTZ}(n,n-k)$ holds — with the* **same** *threshold $1/\sqrt{n}$ on both sides
(no rescaling of the constant occurs).*

**Proof.** (a) Both multisets $\sigma(A_I)$ (size $k \ge 1$) and $\sigma(B_{I^c})$ (size
$n-k \ge 1$) are nonempty, so both minima exist. If $n \ge 2k$, Theorem 1(i) gives
$$\sigma_{\min}(B_{I^c}) = \min\left(\sigma(A_I) \cup \{1\}^{\,n-2k}\right)
 = \min\left(\sigma_{\min}(A_I),\, 1\right) = \sigma_{\min}(A_I),$$
where the last equality holds because $\sigma_{\min}(A_I) \le \sigma_{\max}(A_I) \le 1$ by
Lemma 2 (if $n = 2k$ the padding is empty and the claim is immediate). If $n < 2k$,
Theorem 1(ii) gives symmetrically $\sigma_{\min}(A_I) = \min(\sigma_{\min}(B_{I^c}),1) =
\sigma_{\min}(B_{I^c})$, using $\sigma_{\max}(B_{I^c}) \le 1$ (Lemma 2 applied to $B$).
This also covers singular cases: the common value may be $0$.

(b) The map $I \mapsto I^c$ is a bijection from $\{I \subseteq [n] : |I| = k\}$ onto
$\{J \subseteq [n] : |J| = n-k\}$. By (a), the finite value-multisets
$\{\sigma_{\min}(A_I)\}_{|I|=k}$ and $\{\sigma_{\min}(B_J)\}_{|J|=n-k}$ coincide, so their
maxima coincide: $f(A) = f(B)$.

(c) Assume $\mathrm{GTZ}(n,n-k)$; we prove $\mathrm{GTZ}(n,k)$. Let
$A \in \mathbb{R}^{n\times k}$ with $A^TA = I_k$ be arbitrary. Extend the orthonormal
columns $u_1,\dots,u_k$ of $A$ to an orthonormal basis $u_1,\dots,u_n$ of $\mathbb{R}^n$
(possible for any orthonormal set: complete to a basis and apply Gram–Schmidt, which leaves
the initial orthonormal vectors unchanged), and let $B := [u_{k+1} \cdots u_n] \in
\mathbb{R}^{n\times(n-k)}$. Then $Q = [A\mid B]$ has orthonormal columns, i.e.
$Q \in O(n)$, and $B$ satisfies $B^TB = I_{n-k}$. By $\mathrm{GTZ}(n,n-k)$ applied to $B$,
$f(B) \ge 1/\sqrt{n}$; by (b), $f(A) = f(B) \ge 1/\sqrt{n}$; that is, there exists $I$,
$|I| = k$, with $\sigma_{\min}(A_I) \ge 1/\sqrt{n}$. Since $A$ was arbitrary,
$\mathrm{GTZ}(n,k)$ holds. The converse implication is the same argument with the roles of
$k$ and $n-k$ exchanged (note $[B \mid A]$ is again in $O(n)$, and Theorem 1 applied to it
with the subset $J = I^c$, $|J| = n-k$, returns the pair $(B_J, A_{J^c}) = (B_{I^c}, A_I)$
— the statement is symmetric). The threshold $1/\sqrt{n}$ refers to the same $n$ throughout:
the transfer in (a)–(b) is an exact equality of $\sigma_{\min}$-values, so the constant is
carried over unchanged. $\blacksquare$

**Status: PROVED.**

**Remark (well-definedness on the Grassmannian).** The complement $B$ in (c) is not
unique, but any two orthonormal bases $B, B'$ of $(\operatorname{span} A)^\perp$ satisfy
$B' = BO$ for some $O \in O(n-k)$, whence $B'_J = B_JO$ and $\sigma(B'_J) = \sigma(B_J)$
for every $J$. So $f(B)$ depends only on the subspace $(\operatorname{span} A)^\perp$;
likewise $f(A)$ depends only on $\operatorname{span} A$ (replacing $A$ by $AO'$,
$O' \in O(k)$). Duality is thus the classical involution
$\operatorname{Gr}(k,n) \to \operatorname{Gr}(n-k,n)$, $V \mapsto V^\perp$, and $f$ is a
$\perp$-invariant function on the Grassmannian. This explains structurally why the GTZ
threshold depends only on $n$ and not on $k$ separately.

---

## 3. The case $k = n-2$

### Corollary 2 ($\mathrm{GTZ}(n,n-2)$ for all $n \ge 4$)

*For every $n \ge 4$, $\mathrm{GTZ}(n,n-2)$ holds: every $A \in \mathbb{R}^{n\times(n-2)}$
with $A^TA = I_{n-2}$ has a subset $I \subseteq [n]$, $|I| = n-2$, with
$\sigma_{\min}(A_I) \ge 1/\sqrt{n}$ — the threshold being $1/\sqrt{n}$ of the* **same**
*$n$ (not $1/\sqrt{n-2}$).*

**External input (used exactly once).** The Sengupta–Pautov theorem
[arXiv:2604.05944]: $\mathrm{GTZ}(n,2)$ holds for all $n > 2$. **Exact reference:** the
statement proved is Hypothesis 1 of arXiv:2604.05944 restricted to $k = 2$; the proof
occupies the whole of its Section 2 ("Proof for $k=2$", Cases A and B). We note for the
validator that Hypothesis 1 is the paper's only numbered statement — it contains no
separately numbered theorem for the $k=2$ result, so "Hypothesis 1 for $k=2$, proved in
Section 2" is the sharpest possible citation. (Within that proof, the base case $n = 4$ is
in turn due to Y. Nesterenko, arXiv:2303.07492, cited as [2] there; the case $n = 3$ is
noted as trivial, and $n \ge 5$ is by induction on $n$.)

**Proof of Corollary 2.** Fix $n \ge 4$ and set $k' := 2$, so that $n - k' = n - 2$ and
$1 \le n-2 \le n-1$; the pair $(k', n-k') = (2, n-2)$ is admissible for Corollary 1(c).

*Case $n = 4$.* Here $n - 2 = 2 = n/2$: the statement $\mathrm{GTZ}(4,2)$ **is** the
$k = 2$ theorem at $n = 4$ (self-dual size, cf. Remark 1 after Theorem 1 — duality maps the
problem to itself). It holds by the Sengupta–Pautov theorem (whose $n = 4$ base case is
Nesterenko's result arXiv:2303.07492 — so this particular case was already known
independently of the new content of arXiv:2604.05944). Threshold: $1/\sqrt{4} = 1/2$.

*Case $n \ge 5$.* By the Sengupta–Pautov theorem, $\mathrm{GTZ}(n,2)$ holds — with
threshold $1/\sqrt{n}$ for this same $n$. By Corollary 1(c) with $k = n-2$ (equivalently,
applied to the pair $\{2, n-2\}$), $\mathrm{GTZ}(n, n-2) \Leftrightarrow \mathrm{GTZ}(n,2)$,
and the threshold transfers as an exact equality of $\sigma_{\min}$-values — it remains
$1/\sqrt{n}$, the $n$ being the common number of rows. Hence $\mathrm{GTZ}(n,n-2)$ holds.
$\blacksquare$

**Threshold audit (against the $1/\sqrt{n-2}$ trap).** Both $\mathrm{GTZ}(n,2)$ and
$\mathrm{GTZ}(n,n-2)$ concern matrices with the same number of rows $n$; Corollary 1(a)
equates $\sigma_{\min}(A_I)$ ($A$ being $n\times(n-2)$, $|I| = n-2$) with
$\sigma_{\min}(B_{I^c})$ ($B$ being $n\times 2$, $|I^c| = 2$) — no renormalization, no
change of $n$. At no point is $\mathrm{GTZ}$ invoked at $(n-2)$ rows. The sharp witness in
§4 (Part 3) confirms numerically that $f = 1/\sqrt{n}$ (not $1/\sqrt{n-2}$) is the exact
extremal value on the $(n, n-2)$ side.

**Status: PROVED CONDITIONAL ON arXiv:2604.05944, Hypothesis 1 for $k = 2$ as proved in
its Section 2 (the paper's sole numbered statement; no separately numbered theorem exists
for it). The dependency was independently verified in `proofs/sp-verification.md` (verdict
recorded there). This conditional formulation is permanent: this document requires no
status-upgrade edit; the final report presents the unconditional claim if and only if the
verdict recorded in `proofs/sp-verification.md` is CORRECT or CORRECT WITH FIXABLE GAPS.**

All other ingredients of Corollary 2 (Theorem 1, Corollary 1) are unconditionally PROVED
above; the single conditional dependency is the Sengupta–Pautov $k=2$ theorem.

---

## 4. Numeric sanity block

Everything below was run in the workspace venv
(`.venv`, numpy 2.5.1 / scipy 1.18.0) via

```
.venv/bin/python numerics/study-duality-sanity/check_duality.py
```

(deterministic; single seed `20260717` for `numpy.random.default_rng`, consumed by all
`scipy.stats.ortho_group` draws; full script and README in
`numerics/study-duality-sanity/`). Exit code 0; final line `ALL CHECKS PASSED`. Summary of
the four parts, with observed results from the recorded run of 2026-07-17:

**Part 1 — Theorem 1 identity on random instances.** For each
$(n,k) \in \{(5,1),(6,2),(7,3),(8,4),(6,4),(7,5),(9,6),(5,4),(9,2)\}$ — covering
$n > 2k$, $n = 2k$ (self-dual $(8,4)$), $n < 2k$ ($(6,4),(7,5),(9,6),(5,4)$), $k = 1$ and
$k = n-1$ — 200 trials each: draw $Q \in O(n)$ Haar-random, a uniform random $I$,
$|I| = k$, and compare the sorted multiset $\sigma(B_{I^c})$ against $\sigma(A_I)$ padded
with $|n-2k|$ ones on the side prescribed by Theorem 1. Also checked:
$|\sigma_{\min}(A_I) - \sigma_{\min}(B_{I^c})|$ (Corollary 1(a)) and the corank identity
(Remark 2). Observed: worst multiset discrepancy over all 1800 trials $2.054\times10^{-15}$,
worst $\sigma_{\min}$ gap $7.8\times10^{-16}$, corank identity true in all trials.

**Part 2 — singular and near-singular $A_I$, both orientations.** Constructed completions
$A = \begin{pmatrix} M \\ C \\ 0\end{pmatrix}$ with prescribed top block $M$ ($C^TC = I - M^TM$
via eigendecomposition), $B$ an orthonormal null-space basis of $A^T$: $(7,3)$ with
$\sigma(M) \ni 10^{-12}$ and with $\sigma(M) \ni 0$ exactly; $(7,5)$ and $(6,4)$
($n < 2k$, hence $M$ built with $2k-n$ unit singular values as Theorem 1(ii) forces,
plus a $10^{-12}$ resp. exact $0$); self-dual $(8,4)$ near-singular. Observed: all
orthonormality residuals $\le 1.3\times10^{-15}$, all multiset gaps
$\le 4.5\times10^{-16}$, all $\sigma_{\min}$ gaps $\le 2\times10^{-16}$, corank identity
true in every case (including the exactly singular ones).

**Part 3 — sharp witness for Corollary 2.** For $n = 5,\dots,9$: built the S–P $k=2$
extremal family from the MEMORY.md formula ($n-2$ rows $(a,0)$, rows $(b,c),(b,-c)$,
$a = \sqrt{(n-1)/(n(n-2))}$, $b = 1/\sqrt{2n}$, $c = 1/\sqrt2$; verified
$\|A^TA - I_2\| \le 2.2\times10^{-16}$), took its orthonormal complement
$B \in \mathbb{R}^{n\times(n-2)}$, and computed both extremal functionals by exhaustive
enumeration of subsets. Observed: $|f(A) - 1/\sqrt n| \le 5.6\times10^{-17}$ and
$|f(B) - 1/\sqrt n| \le 2.8\times10^{-16}$ for every $n$ — the dualized family attains the
threshold $1/\sqrt{n}$ (of the same $n$) exactly, confirming sharpness on the $(n,n-2)$
side and the threshold audit above. The known $(4,2)$ extremal matrix from MEMORY.md gives
$f(A) = f(B) = 1/2$ to $1.2\times10^{-16}$ (self-dual case).

**Part 4 — random spot check of Corollary 2's assertion.** For $n = 4,\dots,9$, 100 Haar
draws each of $A \in \mathbb{R}^{n\times(n-2)}$ with orthonormal columns: minimum observed
margin $f(A) - 1/\sqrt n$ per $n$: $+7.23\times10^{-2}$ ($n{=}4$), $+5.93\times10^{-2}$
($5$), $+1.09\times10^{-1}$ ($6$), $+1.21\times10^{-1}$ ($7$), $+1.37\times10^{-1}$ ($8$),
$+1.24\times10^{-1}$ ($9$) — no violation, consistent with (and of course not a proof of)
the corollary; the proof is §3.

Per the project rigour rules: these checks support but never substitute for the proofs;
the status of each statement is determined solely by the arguments of §§1–3.

---

## 5. Statement index for validators

- **Theorem 1** — spectral duality $\sigma(B_{I^c})$ vs $\sigma(A_I)$, both orientations
  $n \ge 2k$ and $n < 2k$, $|n-2k|$ padded ones on the wider side. PROVED.
- **Corollary 1(a)** — $\sigma_{\min}(A_I) = \sigma_{\min}(B_{I^c})$ for all $I$, $|I|=k$,
  $1 \le k \le n-1$. PROVED.
- **Corollary 1(b)** — $f(A) = f(B)$. PROVED.
- **Corollary 1(c)** — $\mathrm{GTZ}(n,k) \Leftrightarrow \mathrm{GTZ}(n,n-k)$, same
  threshold $1/\sqrt n$. PROVED.
- **Corollary 2** — $\mathrm{GTZ}(n,n-2)$ for all $n \ge 4$. PROVED CONDITIONAL ON
  arXiv:2604.05944 (Hypothesis 1 for $k=2$, proved in its Section 2); dependency
  independently verified in `proofs/sp-verification.md` (verdict recorded there).
