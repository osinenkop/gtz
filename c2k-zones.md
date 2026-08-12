# C2k zones: Master Identity, frame bound, and the spiked/flat zones

> **STATUS (updated):** the overarching conjecture **C2k is refuted** in
> `proofs/c2k-middle.md` (VAL-THM-MID-001; explicit Stiefel counterexample
> $(259,5)$ in the open middle zone). The **Master Identity (Theorem MI),
> Lemma Δ, Theorem K, and Propositions S and F below are unconditional and
> unaffected** — they are identities / zone theorems that hold regardless
> (MI holds *at* the counterexample); only the *global* Conjecture C2k and its
> middle‑zone candidates $m_Y$, $m_Y^{(2)}$ are false. The proofs below are
> unchanged.

**Target:** VAL-THM-C2K-001. Promoted from attack round 3
(`attacks/round-3-c2k.md`; numerics `numerics/study-attack-r3-c2k/`, scripts
s4–s7). Fresh numeric sanity for this document:
`numerics/study-c2k-zones/` (§9). Notation follows MEMORY.md. The $k=2$
anchor is Theorem C2 (`proofs/theorem-c2.md`, VAL-THM-C2-001).

---

## 0. Scope, status, what is and is not proved

**This document proves the *unconditional* half of the C2k program.** It does
**not** prove C2k in full, and C2k is **not** a proof of the GTZ hypothesis.
The following three guards are binding and stated up front.

> **SCOPE GUARD 1 — C2k is only partially proved.** C2k (Conjecture below) is
> established here in two explicit regimes — the **spiked zone** $Q \ge
> Q^*(n,k)$ (Proposition S) and the **flat window** $x \le 1/(1-\omega_k)$
> around uniform leverage (Proposition F, $k \ge 3$). The complementary
> **middle zone** ($x > 1/(1-\omega_k)$ **and** $Q < Q^*(n,k)$) is **OPEN**.

> **SCOPE GUARD 2 — the middle-zone candidate is CONJECTURED.** The
> inequality $m_Y \ge 0$ of §8, which would close the middle zone (and hence
> all of C2k), is **CONJECTURED**: it is numerically supported (no violation
> found in extensive adversarial search) but **not proved**. It is labelled
> CONJECTURED wherever it appears. In contrast, the pure-$\kappa$ sufficient
> condition $H + \Delta \ge c_k\kappa$ is **FALSE** in the middle zone
> (proved by counterexample, §8).

> **SCOPE GUARD 3 — C2k is a weighted average, not GTZ.** Even fully proved,
> C2k is the general-$k$ *weighted-average* pair inequality $\Sigma(A) \ge
> 0$ (Conjecture below). Nonnegativity of this average yields the *existence*
> of one good pair (as at $k=2$, Theorem C2 / Cor 5.1 there), **not** the
> existence of a good $k$-subset. It does **not** prove GTZ$(n,k)$. The
> volume-weighted $k$-subset analogs are provably FALSE for $k \ge 3$
> (`attacks/round-3-variational.md` §6, obstructions O1/O2). Do not overclaim.

### Notation and the conjecture

Fix $A \in \mathrm{St}(n,k)$, i.e. a real $n\times k$ matrix with $A^TA =
I_k$, $k \ge 2$. Rows $r_i \in \mathbb R^k$, leverages $\ell_i = \|r_i\|^2$,
inner products $c_{ij} = \langle r_i, r_j\rangle$ (so $c_{ii} = \ell_i$),
unit rows $\hat r_i = r_i/\sqrt{\ell_i}$ (for $\ell_i > 0$), and $Q :=
\sum_i \ell_i^2$. Throughout,
$$s := \frac{k+2}{k^2}, \quad t := \frac{2k}{(k+2)n}, \quad
\theta := \frac{k-1}{k}, \quad
\Sigma(A) := \sum_{i<j}(\ell_i\ell_j - c_{ij}^2)\bigl[s(\ell_i-t)(\ell_j-t)
- c_{ij}^2\bigr].$$

> **Conjecture C2k.** $\Sigma(A) \ge 0$ for every $A \in \mathrm{St}(n,k)$,
> $k \ge 2$. (At $k = 2$ this is Theorem C2, PROVED in
> `proofs/theorem-c2.md`; here it is the weighted-average form of the
> general-$k$ Pair Lemma `proofs/pair-lemma.md`.)

Two standard facts used repeatedly: $\sum_i r_ir_i^T = A^TA = I_k$ (hence
$\sum_i \ell_i = \operatorname{Tr}I_k = k$); and $AA^T$ is an orthogonal
projector, so $\ell_i = (AA^T)_{ii} \in [0,1]$ and, by Cauchy–Schwarz
$(\sum_i\ell_i)^2 \le n\sum_i\ell_i^2$, $Q \ge k^2/n$.

### Status table

| statement | status |
|---|---|
| Lemma 1.1 (isotypic projection constants) | PROVED |
| Eq. (1.2) (spin-4 energy $H = \|\pi_4M_4\|^2 \ge 0$) | PROVED |
| **Theorem MI** (master identity $2\Sigma = H + \Delta - c_kV$) + Cor MI (crux) | PROVED |
| Lemma Δ (factorization + split (3.1)) | PROVED |
| Lemma 4.1 (frame bound $T_\nu \preceq \mathrm{Id}$) | PROVED |
| **Theorem K** ($V \le \kappa = \sum\ell^3 - Q^2/k$) | PROVED |
| **Proposition S** (spiked zone $Q \ge Q^*(n,k)$; closed-form $Q^*_\infty$) | PROVED |
| Lemma B, Lemma T (moment bounds) | PROVED |
| Lemma Y (harmonic Rayleigh grab $H \ge \mathrm{num}^2/\mathrm{den}$) | PROVED |
| **Proposition F** (flat window $x \le 1/(1-\omega_k)$, $k \ge 3$) | PROVED |
| Equality manifold / icosahedral ETF exact | PROVED |
| $H + \Delta \ge c_k\kappa$ globally | **FALSE** (proved by counterexample) |
| $m_Y \ge 0$ globally (would close C2k) | **CONJECTURED** (numerically supported) |
| C2k on the middle zone | **OPEN** |
| C2k in full; GTZ$(n,k)$ | **not proved / not claimed** (Scope Guards 1, 3) |

**Inputs used.** Column orthonormality of $A$; the Frobenius inner product on
symmetric tensors; Cauchy–Schwarz and the Rayleigh characterization of
$\lambda_{\max}$ (spectral theorem for real symmetric matrices); elementary
$O(k)$-representation theory of symmetric 4-tensors (self-contained, Lemma
1.1); standard extremal-trace bounds for a symmetric traceless matrix (Lemma
T). No Perron–Frobenius theory, no induction, no external preprint. Theorem
C2 is cited only as the $k=2$ specialization (§2, Remark 2.2).

---

## 1. The isotypic decomposition of symmetric 4-tensors

Work in $\mathrm{Sym}^4(\mathbb R^k)$, the space of symmetric 4-tensors on
$\mathbb R^k$, with the Frobenius inner product inherited from $(\mathbb
R^k)^{\otimes 4}$; for $a, b \in \mathbb R^k$, $\langle a^{\otimes 4},
b^{\otimes 4}\rangle = \langle a, b\rangle^4$. Under $O(k)$ this space splits
orthogonally into the trivial (spin-0), the spin-2, and the harmonic
(spin-4) isotypic pieces. Concrete carriers:

- **spin-0:** $E := \mathrm{Sym}(\delta\otimes\delta)$, i.e. $E_{abcd} =
  \tfrac13(\delta_{ab}\delta_{cd} + \delta_{ac}\delta_{bd} +
  \delta_{ad}\delta_{bc})$;
- **spin-2:** for symmetric traceless $B$, $F_B := \mathrm{Sym}(\delta\otimes
  B) = \tfrac16[\delta_{ab}B_{cd} + \delta_{cd}B_{ab} + \delta_{ac}B_{bd} +
  \delta_{bd}B_{ac} + \delta_{ad}B_{bc} + \delta_{bc}B_{ad}]$.

**Lemma 1.1 (projection constants).** For $a \in \mathbb R^k$ and symmetric
traceless $B, C$:

1. $\|E\|^2 = \dfrac{k(k+2)}3$;
2. $\langle F_B, F_C\rangle = \dfrac{k+4}6\langle B, C\rangle$ and $\langle
   E, F_B\rangle = 0$;
3. $\langle a^{\otimes 4}, E\rangle = |a|^4$ and $\langle a^{\otimes 4},
   F_B\rangle = |a|^2\,a^TBa$.

*Proof.* Since $E, F_B$ are already symmetric, $\langle X,
\mathrm{Sym}(T)\rangle = \langle X, T\rangle$ for such $X$; so we may pair
against the *unsymmetrized* $\delta\otimes\delta$, $\delta\otimes B$,
$a^{\otimes 4}$.

(3) $\langle a^{\otimes4}, E\rangle = \tfrac13\cdot 3\cdot(a^Ta)(a^Ta) =
|a|^4$; each of the three pairings of $E$ contracts $a^{\otimes4}$ to
$|a|^2\cdot|a|^2$. Likewise $\langle a^{\otimes4}, F_B\rangle = \tfrac16\cdot
6\cdot|a|^2(a^TBa) = |a|^2a^TBa$ (each of the six pairing terms gives
$|a|^2\,a^TBa$).

For (1)–(2) contract one index pair first. Summing $E$ over its first pair,
$\sum_a E_{aacd} = \tfrac13(k\delta_{cd} + \delta_{cd} + \delta_{cd}) =
\tfrac{k+2}3\delta_{cd}$. Hence $\|E\|^2 = \langle E, \delta\otimes\delta
\rangle = \sum_{cd}\bigl(\sum_a E_{aacd}\bigr)\delta_{cd} = \tfrac{k+2}3
\sum_{cd}\delta_{cd}^2 = \tfrac{k+2}3 k$, and $\langle E, F_B\rangle$: summing
$\delta_{ab}B_{cd}$'s first pair against $E$ gives $\tfrac{k+2}3\sum_{cd}
\delta_{cd}B_{cd} = \tfrac{k+2}3\operatorname{Tr}B = 0$. For $\langle F_B,
F_C\rangle = \langle F_B, \delta\otimes C\rangle = \sum_{cd}C_{cd}\sum_a
(F_B)_{aacd}$: setting $b=a$ and summing the six terms of $F_B$, the term
$\delta_{ab}B_{cd}$ gives $kB_{cd}$, the term $\delta_{cd}B_{ab}$ gives
$\delta_{cd}\operatorname{Tr}B = 0$, and each of the four crossed terms gives
$B_{cd}$; total $\tfrac{k+4}6 B_{cd}$, so $\langle F_B, F_C\rangle =
\tfrac{k+4}6\langle B, C\rangle$. $\blacksquare$

*(Exact-arithmetic cross-check for $k=2,3$: `check_c2k_zones.py` block G;
$k=2,3,4$ in `s5_exact_constants.py` block B.)*

**The 4th-moment tensor.** Put $M_4 := \sum_i r_i^{\otimes 4}$, so $\|M_4\|^2
= \sum_{ij}c_{ij}^4 =: \Phi_4$. Let
$$S := \sum_i \ell_i\, r_ir_i^T = A^TD_\ell A, \qquad
S_0 := S - \frac Qk I_k, \qquad V := \|S_0\|_F^2 ,$$
where $D_\ell = \operatorname{diag}(\ell_i)$. Since $\operatorname{Tr}S =
\sum_i\ell_i\|r_i\|^2 = \sum_i\ell_i^2 = Q$, we have $\operatorname{Tr}S_0 =
0$. By Lemma 1.1(3),
$$\langle M_4, E\rangle = \sum_i|r_i|^4 = Q, \qquad
\langle M_4, F_B\rangle = \sum_i \ell_i\, r_i^TBr_i = \operatorname{Tr}(SB) =
\langle S_0, B\rangle \ \ (\operatorname{Tr}B = 0),$$
the last using $\operatorname{Tr}(\tfrac Qk B) = 0$. Orthogonally project
$M_4 = \pi_0M_4 + \pi_2M_4 + \pi_4M_4$ onto spin-0 / spin-2 / spin-4:
$$\|\pi_0M_4\|^2 = \frac{\langle M_4,E\rangle^2}{\|E\|^2} = \frac{3Q^2}{k(k+2)},
\qquad
\pi_2M_4 = F_{B^*},\ B^* = \tfrac6{k+4}S_0,\ \|\pi_2M_4\|^2 = \frac{6V}{k+4},$$
where $B^*$ is forced by $\tfrac{k+4}6\langle B^*,B\rangle = \langle
S_0,B\rangle$ (Lemma 1.1(2)) and $\|\pi_2M_4\|^2 = \tfrac{k+4}6\|B^*\|^2 =
\tfrac6{k+4}V$. By Pythagoras the **spin-4 (harmonic) energy** is
$$H := \|\pi_4M_4\|^2 = \Phi_4 - \frac{3Q^2}{k(k+2)} - \frac{6V}{k+4}
\;\ge\; 0. \tag{1.2}$$

*Status: **PROVED**. (Nonnegativity is Bessel/Pythagoras;
`check_c2k_zones.py` block B verifies $H = \|\pi_4M_4\|^2$ by building $M_4$
and projecting explicitly, rel err $\le 1.1\cdot10^{-16}$, and confirms $H
\ge 0$ on 4000 Haar samples, min $+6.7\cdot10^{-6}$.)*

---

## 2. Theorem MI: the master identity and the crux

**Theorem MI.** For every $A \in \mathrm{St}(n,k)$, with $\Delta(Q) :=
\theta s(Q - tk)^2 - \theta Q^2/(k+2)$ and
$$c_k := (1+s) - \frac6{k+4} = \frac{k^3 - k^2 + 6k + 8}{k^2(k+4)}
= \frac{(k+1)(k^2 - 2k + 8)}{k^2(k+4)},$$
$$\boxed{\,2\,\Sigma(A) = \Phi_4 + \theta s(Q - tk)^2 - \frac{Q^2}k -
(1+s)V \;=\; H + \Delta(Q) - c_kV.\,} \tag{MI}$$

*Proof.* Set $w_i := r_ir_i^T - \frac{\ell_i}k I_k$ (symmetric, traceless)
and $K_{ij} := \langle w_i, w_j\rangle_F$. A direct expansion gives
$$K_{ij} = \operatorname{Tr}(r_ir_i^Tr_jr_j^T) - \tfrac{2\ell_i\ell_j}k
+ \tfrac{\ell_i\ell_j}{k^2}\operatorname{Tr}I = c_{ij}^2 - \frac{\ell_i\ell_j}k ,$$
and $\sum_i w_i = \sum_i r_ir_i^T - \tfrac1k(\sum_i\ell_i)I = I_k - I_k = 0$,
whence $(K\mathbf 1)_i = \langle w_i, \sum_j w_j\rangle = 0$, i.e.
$K\mathbf 1 = 0$. Let $L := \ell\ell^T$ and $\theta L - K$. Off the diagonal
$(\theta L - K)_{ij} = (\theta + \tfrac1k)\ell_i\ell_j - c_{ij}^2 =
\ell_i\ell_j - c_{ij}^2$ (as $\theta + \tfrac1k = 1$); on the diagonal
$(\theta L - K)_{ii} = \theta\ell_i^2 - (\ell_i^2 - \tfrac{\ell_i^2}k) =
\theta\ell_i^2 - \theta\ell_i^2 = 0$. The weight $\ell_i\ell_j - c_{ij}^2$
also vanishes on the diagonal ($c_{ii} = \ell_i$), so the full matrix
$(\ell_i\ell_j - c_{ij}^2)_{i,j}$ (diagonal $:=0$) **equals** $\theta L - K$.
The summand of $\Sigma$ is symmetric in $i,j$ and vanishes at $i=j$, so with
$u := \ell - t\mathbf 1$ and $c_{ij}^2 = K_{ij} + \tfrac1k L_{ij}$,
$$2\Sigma = \sum_{i,j}(\theta L - K)_{ij}\bigl[s\,u_iu_j - K_{ij} -
\tfrac1k L_{ij}\bigr]
= s\,u^T(\theta L - K)u - \langle\theta L - K, K\rangle - \tfrac1k\langle
\theta L - K, L\rangle .$$
Evaluate each piece using $u^T\ell = Q - tk$ (from $\sum_i\ell_i = k$),
$u^TL u = (u^T\ell)^2$, $u^TKu = \ell^TK\ell$ (as $K\mathbf 1 = 0$),
$\langle L, K\rangle = \ell^TK\ell$, $\langle L, L\rangle = (\ell^T\ell)^2 =
Q^2$:
$$2\Sigma = s\theta(Q-tk)^2 - \frac\theta k Q^2 + \|K\|_F^2 -
\tau\,\ell^TK\ell, \qquad \tau := s + \theta - \tfrac1k = \frac{k^2-k+2}{k^2}.$$
Two exact translations close it:

- **(i) $\ell^TK\ell = V$.** Indeed $\ell^TK\ell = \sum_{ij}\ell_i\ell_j
  (c_{ij}^2 - \tfrac{\ell_i\ell_j}k) = \operatorname{Tr}(S^2) - \tfrac{Q^2}k$,
  because $\sum_{ij}\ell_i\ell_jc_{ij}^2 = \operatorname{Tr}\bigl[(\sum_i\ell_i
  r_ir_i^T)(\sum_j\ell_j r_jr_j^T)\bigr] = \operatorname{Tr}(S^2)$; and
  $\operatorname{Tr}(S^2) = \|S\|_F^2 = \|S_0 + \tfrac Qk I\|^2 = V + \tfrac{Q^2}
  k$, so $\ell^TK\ell = V$. *(This is the exact "overlap" quantity — it
  replaces the round-3 relaxation $v^TKv \le \nu_1\|v\|^2$; no bound needed.)*
- **(ii) $\|K\|_F^2 = \Phi_4 - \tfrac{2V}k - \tfrac{Q^2}{k^2}$.** From
  $\sum_{ij}(c_{ij}^2 - \tfrac{\ell_i\ell_j}k)^2 = \Phi_4 - \tfrac2k
  \operatorname{Tr}(S^2) + \tfrac{Q^2}{k^2}$ and $\operatorname{Tr}(S^2) = V +
  Q^2/k$.

Substituting (i)–(ii): the $Q^2$ coefficient (outside the retained
$\theta s(Q-tk)^2$) is $-\tfrac\theta k - \tfrac1{k^2} = -\tfrac1k$, and the
$V$ coefficient is $-\tfrac2k - \tau = -(1+s)$; this is the **first** equality
of (MI). For the **second**, substitute $\Phi_4 = H + \tfrac{3Q^2}{k(k+2)} +
\tfrac{6V}{k+4}$ from (1.2) and write $\theta s(Q-tk)^2 = \Delta(Q) +
\tfrac{\theta}{k+2}Q^2$: the total $Q^2$-coefficient becomes
$\tfrac{\theta}{k+2} + \tfrac3{k(k+2)} - \tfrac1k = 0$ (numerator
$\theta k + 3 - (k+2) = (k-1) + 3 - k - 2 = 0$), and the $V$-coefficient
becomes $\tfrac6{k+4} - (1+s) = -c_k$. $\blacksquare$

*(Constants exact for symbolic $k$ and MI exact on rational Cayley–Stiefel
matrices: `check_c2k_zones.py` block G, `s5_exact_constants.py` blocks A/C;
numeric MI on 4000 Haar $\le 1.4\cdot10^{-15}$ relative, block A.)*

**Corollary MI (crux).** C2k $\iff H + \Delta(Q) \ge c_kV$ on
$\mathrm{St}(n,k)$. In particular any configuration with **uniform leverage**
($\ell_i \equiv k/n$, forcing $S = \tfrac kn I$, $Q = k^2/n$, hence $V = 0$
and $\Delta = 0$) and **spin-4-free** 4th moment ($H = 0$) has $\Sigma = 0$.

*Status: **PROVED** (both equalities of (MI); $c_k > 0$ since its numerator is
positive for $k\ge1$).* 

**Remark 2.2 ($k=2$ anchor).** At $k=2$, $s = 1$, $t = 1/n$, $\theta =
1/2$, and $\ell_i\ell_j - c_{ij}^2 = p_{ij}^2$ (the squared Plücker minor),
so $\Sigma = \sum_{i<j}p_{ij}^2[(\ell_i - \tfrac1n)(\ell_j - \tfrac1n) -
c_{ij}^2] = \sum_{i<j}p_{ij}^2 b_{ij}$ in the notation of
`proofs/theorem-c2.md` (Remark 2.4 there). Theorem C2 (VAL-THM-C2-001) proves
$\Sigma \ge 0$ at $k=2$ **unconditionally**; the present document extends the
identity machinery to all $k$ and closes two zones for $k \ge 3$.

---

## 3. The Δ-factorization

**Lemma Δ.** With $a_k := \dfrac{4(k^2-1)}{k^3(k+2)}$ and $W := Q -
\dfrac{k^2}n \ge 0$,
$$\Delta(Q) = a_k\Bigl(Q - \frac{k^2}n\Bigr)\Bigl(Q - \frac{k^2}{(k+1)n}\Bigr)
\;=\; \frac{4\theta}{k(k+2)}\,WQ + \frac{a_k}{k+1}W^2 . \tag{3.1}$$
Hence $\Delta \ge 0$ on the admissible range $Q \ge k^2/n$, vanishing iff
$Q = k^2/n$ (uniform leverage); and $\Delta \ge \tfrac{4\theta}{k(k+2)}WQ$.

*Proof.* $\Delta$ is a quadratic in $Q$ with leading coefficient
$\theta(s - \tfrac1{k+2}) = \theta\cdot\tfrac{4(k+1)}{k^2(k+2)} = a_k$. Its
roots solve $\sqrt s\,(Q - tk) = \pm Q/\sqrt{k+2}$, i.e. (using $\sqrt s =
\sqrt{k+2}/k$) $Q - tk = \pm\tfrac{k}{k+2}Q$, giving $Q = k^2/n$ (sign $+$)
and $Q = k^2/((k+1)n)$ (sign $-$); this is the first equality. For the split,
$\tfrac{4\theta}{k(k+2)} = a_k\tfrac k{k+1}$, so
$\tfrac{4\theta}{k(k+2)}WQ + \tfrac{a_k}{k+1}W^2 = \tfrac{a_k}{k+1}W(kQ + W) =
\tfrac{a_k}{k+1}W\bigl((k+1)Q - \tfrac{k^2}n\bigr) = a_kW\bigl(Q -
\tfrac{k^2}{(k+1)n}\bigr)$, which is the factored form. Both roots are $\le
k^2/n$ and both $W, Q > 0$ on the admissible range, so both bracket factors
are $\ge 0$. $\blacksquare$

*Status: **PROVED** (polynomial identities, exact for symbolic $(Q,n,k)$:
`check_c2k_zones.py` block G, `s5_exact_constants.py` A2/A3).*

---

## 4. Theorem K: the frame bound and the κ bound

**Lemma 4.1 (frame bound).** For every symmetric $B \in \mathbb R^{k\times k}$,
$$\sum_{i:\,\ell_i > 0} \ell_i\,(\hat r_i^TB\hat r_i)^2 \;\le\;
\operatorname{Tr}(B^2). \tag{4.1}$$
Equivalently the leverage-measure frame operator $T_\nu := \sum_i\ell_i\,
\operatorname{vec}(\hat r_i\hat r_i^T)\operatorname{vec}(\hat r_i\hat
r_i^T)^T$ satisfies $T_\nu \preceq \mathrm{Id}$ — an exact,
$\ell_{\max}$-free bound.

*Proof.* For each nonzero row, Cauchy–Schwarz on the unit vector $\hat r_i$:
$(\hat r_i^TB\hat r_i)^2 \le \|\hat r_i\|^2\|B\hat r_i\|^2 = \hat r_i^TB^2
\hat r_i$. Sum with weights $\ell_i$ and use the exact frame identity
$\sum_i\ell_i\hat r_i\hat r_i^T = \sum_i r_ir_i^T = I_k$:
$\sum_i\ell_i\hat r_i^TB^2\hat r_i = \operatorname{Tr}(B^2\sum_i r_ir_i^T) =
\operatorname{Tr}(B^2)$. The frame-operator form is (4.1) applied to
$\operatorname{vec}(B)$, since $\operatorname{vec}(B)^TT_\nu\operatorname{vec}
(B) = \sum_i\ell_i(\hat r_i^TB\hat r_i)^2$ and $\|\operatorname{vec}(B)\|^2 =
\operatorname{Tr}(B^2)$. $\blacksquare$

**Theorem K.** $V \le \kappa := \sum_i\ell_i^3 - \dfrac{Q^2}k$. Moreover
$\kappa \ge 0$ always, with $\kappa = 0$ iff all nonzero leverages are equal;
and $V = \kappa$ at every orthonormal-block configuration ($j$ rows forming an
orthonormal set at $\ell = 1$, the rest zero).

*Proof.* Let $\beta_i := r_i^TS_0r_i$. Then $\sum_i\beta_i =
\operatorname{Tr}\bigl(S_0\sum_i r_ir_i^T\bigr) = \operatorname{Tr}S_0 = 0$,
and $\sum_i\ell_i\beta_i = \operatorname{Tr}(S_0S) = \operatorname{Tr}(S_0(S_0
+ \tfrac Qk I)) = V$. Write $\beta_i = \ell_i\hat\beta_i$ with $\hat\beta_i =
\hat r_i^TS_0\hat r_i$ (and $\beta_i = 0$ when $\ell_i = 0$). From
$\sum_i\ell_i\hat\beta_i = \sum_i\beta_i = 0$,
$$V = \sum_i\ell_i^2\hat\beta_i = \sum_i\ell_i\Bigl(\ell_i - \tfrac Qk\Bigr)
\hat\beta_i \;\le\; \Bigl(\sum_i\ell_i(\ell_i - \tfrac Qk)^2\Bigr)^{1/2}
\Bigl(\sum_i\ell_i\hat\beta_i^2\Bigr)^{1/2} = \sqrt\kappa\cdot\sqrt V ,$$
using Cauchy–Schwarz with weights $\ell_i$, the identity $\sum_i\ell_i(\ell_i
- \tfrac Qk)^2 = \sum\ell^3 - \tfrac{2Q}k\sum\ell^2 + \tfrac{Q^2}{k^2}\sum\ell
= \sum\ell^3 - \tfrac{Q^2}k = \kappa$ (as $\sum\ell = k$), and Lemma 4.1 with
$B = S_0$ for the second factor. $\kappa \ge 0$ is Cauchy–Schwarz
$(\sum_i\ell_i^2)^2 \le (\sum_i\ell_i)(\sum_i\ell_i^3)$, i.e. $Q^2 \le
k\sum\ell^3$; so $\sqrt\kappa$ is real and $V > 0 \Rightarrow V \le \kappa$,
while $V = 0 \Rightarrow V = 0 \le \kappa$. Equality $\kappa = 0$ is the
Cauchy–Schwarz equality $\ell^{3/2}\propto\ell^{1/2}$ on the support (all
nonzero $\ell_i$ equal). At an orthonormal block of $j$ unit rows, $S =
P_{\mathrm{block}}$ (rank-$j$ projector), $Q = j$, $V = \operatorname{Tr}(S^2)
- \tfrac{Q^2}k = j - \tfrac{j^2}k$, $\sum\ell^3 = j$, $\kappa = j - \tfrac{j^2}
k = V$. $\blacksquare$

$\kappa$ is the variance of the leverage function under the leverage measure
$\nu = \sum_i\ell_i\delta_{\hat r_i}$.

*Status: **PROVED**. (`check_c2k_zones.py` block C: frame bound holds on 3000
Haar samples with random $B$, max excess $\le -7.7\cdot10^{-4}$; $\min(\kappa
- V) = +1.6\cdot10^{-3}$ on Haar; $|\kappa - V| = 0$ exactly at orthonormal
blocks. Exact $\kappa - V \in \mathbb Q_{\ge0}$: block G, `s5` block C.)*

---

## 5. Proposition S: the spiked zone

**Proposition S.** If
$$\Delta(Q) \;\ge\; c_k\,Q\Bigl(1 - \frac Qk\Bigr), \tag{5.1}$$
then $\Sigma(A) \ge 0$. Condition (5.1) holds exactly for $Q \ge Q^*(n,k)$,
the larger root of the upward parabola
$$\bigl(a_k + \tfrac{c_k}k\bigr)Q^2 - \bigl(2\theta stk + c_k\bigr)Q +
\theta st^2k^2 = 0, \tag{5.2}$$
and
$$Q^*(n,k) \;\xrightarrow[n\to\infty]{}\; Q^*_\infty(k) =
\frac{kc_k}{ka_k + c_k} = \frac{k^3 + 4k + 16}{k^2 + 4k + 16}
\qquad (\text{monotone from above}),$$
with $Q^*_\infty(k) < k$ for every $k \ge 2$ (equivalently $4(k+4)(k-1) > 0$).
Numerically $Q^*_\infty/k = 0.496,\,0.500,\,0.528,\,0.561,\,0.594,\,0.625$
for $k = 3,\dots,8$.

*Proof.* Every leverage satisfies $\ell_i \le 1$, so $\sum_i\ell_i^3 \le
\sum_i\ell_i^2 = Q$ and $\kappa \le Q - Q^2/k = Q(1 - Q/k)$. If (5.1) holds
then, since $c_k > 0$, Theorem K gives $c_kV \le c_k\kappa \le c_k Q(1-Q/k)
\le \Delta$, so by Theorem MI $2\Sigma = H + \Delta - c_kV \ge H \ge 0$.
Condition (5.1) is $\Delta(Q) - c_kQ(1-Q/k) \ge 0$; expanding, the left side
is the parabola (5.2) with positive leading coefficient $a_k + c_k/k$ and
positive constant term $\theta st^2k^2$, hence $\ge 0$ exactly outside its two
positive roots, and in particular for $Q \ge Q^*(n,k) :=$ (larger root). As
$n\to\infty$, $t = 2k/((k+2)n)\to 0$, so $\Delta(Q)\to a_kQ^2$ and (5.1)
becomes $(a_k + c_k/k)Q \ge c_k$, i.e. $Q \ge kc_k/(ka_k + c_k)$; the closed
form and $Q^*_\infty < k$ are elementary algebra. $\blacksquare$

This settles, **unconditionally and for every $n$**, the entire spiked regime
on which the original step-(d) relaxation leaked (finite-$n$ thresholds
$Q^*(6,3) = 2.31$, $Q^*(12,3) = 1.92$, $Q^*(100,3) = 1.54$ vs $Q^*_\infty(3)
= 1.486$; orthonormal-block limits included, with strict margin).

*Status: **PROVED**. (`check_c2k_zones.py` block D: on 2493 spiked samples
with $Q \ge Q^*$, $\min(\Delta - c_k\kappa) = +0.124$ and $\min 2\Sigma =
+0.694$, no violation. Closed form $Q^*_\infty$ exact: block G,
`s5` A5.)*

---

## 6. Lemmas for the flat zone

Throughout this section $V > 0$ and $W = Q - k^2/n > 0$ (the boundary cases
$V = 0$ or $W = 0$ fall under Proposition F's trivial branch). Let $B_q :=
\sum_i\beta_i^2$ with $\beta_i = r_i^TS_0r_i$, $T_3 :=
\operatorname{Tr}S_0^3$, $P_4 := \operatorname{Tr}S_0^4$.

**Lemma B (moment lower bound).** $B_q \ge V^2/W$.

*Proof.* $\sum_i\beta_i = 0$ and $\sum_i\ell_i\beta_i = V$ (proof of Theorem
K). For $v := \ell - \tfrac kn\mathbf 1$ (centered leverage, mean $k/n$),
$\sum_i v_i\beta_i = V - \tfrac kn\cdot0 = V$ and $\|v\|^2 = Q - k^2/n = W$.
Cauchy–Schwarz: $V = \sum_i v_i\beta_i \le \|v\|\,\|\beta\| = \sqrt W\sqrt
{B_q}$. $\blacksquare$

**Lemma T (traceless trace bounds; standard).** For the symmetric traceless
$k\times k$ matrix $S_0$ with $\operatorname{Tr}S_0^2 = V$:
$$|T_3| \le \gamma_k V^{3/2}, \quad \gamma_k := \frac{k-2}{\sqrt{k(k-1)}},
\qquad \frac{V^2}k \le P_4 \le V^2 .$$

*Proof.* With eigenvalues $\lambda_1,\dots,\lambda_k$ ($\sum\lambda = 0$,
$\sum\lambda^2 = V$): the extreme of $\sum\lambda^3$ is at the spectrum
$((k-1)a, -a,\dots,-a)$ (one eigenvalue distinct — the Lagrange stationary
configuration of $\sum\lambda^3$ on $\{\sum\lambda = 0,\ \sum\lambda^2 =
V\}$), where $a^2 = V/(k(k-1))$ and $|\sum\lambda^3| = k(k-1)(k-2)a^3 =
\gamma_k V^{3/2}$. The quartic upper bound is $\sum\lambda^4 \le
(\sum\lambda^2)^2 = V^2$; the lower is Cauchy–Schwarz $k\sum\lambda^4 \ge
(\sum\lambda^2)^2$. $\blacksquare$

**Lemma Y (harmonic Rayleigh grab).** Let $Y := \mathrm{Sym}(S_0\otimes S_0)$
(quartic form $u\mapsto (u^TS_0u)^2$). Then $\langle Y, E\rangle = \tfrac23V$,
$\langle Y, F_B\rangle = \tfrac23\operatorname{Tr}(S_0^2B)$
($\operatorname{Tr}B = 0$), $\|Y\|^2 = \tfrac{V^2}3 + \tfrac23P_4$, and with
$$\mathrm{num} := \langle\pi_4M_4, \pi_4Y\rangle = B_q - \frac{2QV}{k(k+2)} -
\frac{4T_3}{k+4}, \qquad
\mathrm{den} := \|\pi_4Y\|^2 = \frac{V^2}3 + \frac23P_4 - \frac{4V^2}{3k(k+2)}
- \frac8{3(k+4)}\Bigl(P_4 - \frac{V^2}k\Bigr) \ge 0,$$
one has $H \ge \mathrm{num}^2/\mathrm{den}$ whenever $\mathrm{den} > 0$, and
$\mathrm{den} = 0 \Rightarrow \mathrm{num} = 0$.

*Proof.* Pairing $Y$ (via its three symmetrization terms) against $E$ gives
$\tfrac13\bigl((\operatorname{Tr}S_0)^2 + 2\|S_0\|^2\bigr) = \tfrac23V$;
against $\delta\otimes B$ gives $\sum_{cd}B_{cd}\sum_a Y_{aacd} = \tfrac23
\operatorname{Tr}(S_0^2B)$ (since $\sum_a Y_{aacd} = \tfrac13(\operatorname
{Tr}S_0\,(S_0)_{cd} + 2(S_0^2)_{cd}) = \tfrac23(S_0^2)_{cd}$); and $\|Y\|^2 =
\tfrac13(V^2 + 2P_4)$ (the three same-pairing terms give $V^2$, the crossed
terms $\operatorname{Tr}S_0^4 = P_4$). Now $\langle M_4, Y\rangle =
\sum_i(r_i^TS_0r_i)^2 = B_q$; subtracting the spin-0 and spin-2 overlaps
computed from Lemma 1.1 (with $\pi_2Y = F_C$, $C = \tfrac4{k+4}(S_0^2 -
\tfrac Vk I)$),
$$\langle\pi_0M_4,\pi_0Y\rangle = \frac{Q\cdot\frac23V}{k(k+2)/3} =
\frac{2QV}{k(k+2)}, \qquad
\langle\pi_2M_4,\pi_2Y\rangle = \tfrac{k+4}6\langle\tfrac6{k+4}S_0,
C\rangle = \frac{4T_3}{k+4},$$
gives the stated $\mathrm{num}$; the same subtraction on $\|Y\|^2$ gives
$\mathrm{den}$. Finally Cauchy–Schwarz in the spin-4 subspace, $\mathrm{num}^2
= \langle\pi_4M_4,\pi_4Y\rangle^2 \le \|\pi_4M_4\|^2\|\pi_4Y\|^2 =
H\cdot\mathrm{den}$; if $\mathrm{den} = 0$ then $\pi_4Y = 0$ so $\mathrm{num}
= 0$. $\blacksquare$

*Status: Lemma B **PROVED**; Lemma T **PROVED** (standard extremal-trace
facts); Lemma Y **PROVED** (projection constants exact, `check_c2k_zones.py`
block G and `s5_exact_constants.py` B4/B5; the cross-multiplied Cauchy–Schwarz
$\mathrm{num}^2 \le H\cdot\mathrm{den}$ verified on Haar, block B, max
$\le -1.9\cdot10^{-5}$).*

---

## 7. Proposition F: the flat window

**Proposition F.** Let $k \ge 3$ and set
$$\rho_k := \frac{2\theta}{c_k} > 1, \qquad \alpha_k := 1 - \frac1{\rho_k} =
1 - \frac{c_k}{2\theta} \in (0,1), \qquad
\omega_k := \min\Bigl\{\frac{\theta\alpha_k^2}{c_k^2\,k(k+2)},\;
\frac{\theta\alpha_k^2(k+4)^2}{16\gamma_k^2\,c_k\,k(k+2)}\Bigr\} > 0.$$
If $W \le \omega_k Q$ — equivalently $x := nQ/k^2 \le 1/(1-\omega_k)$ —
then $\Sigma(A) \ge 0$.

*(That $\rho_k > 1$, $\alpha_k > 0$ for $k\ge3$ is $2\theta - c_k =
\dfrac{(k-2)(k^2+9k+4)}{k^2(k+4)} > 0$, $=0$ at $k=2$. Constants are
deliberately **unoptimized**; only positivity of the window is asserted.
$\omega_3 = \tfrac5{242}$, $\omega_4 = \tfrac{49}{1800}$, $\omega_5 =
0.0231$; windows $x \le 1.0211,\,1.0280,\,1.0237$.)*

*Proof.* **Trivial branch.** If $V \le \Delta/c_k$ (in particular $V = 0$, or
$W = 0$ which forces uniform leverage $S = \tfrac kn I$, $V = \Delta = 0$),
then $c_kV \le \Delta$ and $2\Sigma = H + \Delta - c_kV \ge H \ge 0$.

**Grab branch:** assume $V > \Delta/c_k$ and $W > 0$. By (3.1), $\Delta \ge
\tfrac{4\theta}{k(k+2)}WQ$, so
$$\frac VW > \frac{\Delta}{c_kW} \ge \frac{4\theta}{c_k\,k(k+2)}Q =
\frac{2\rho_k Q}{k(k+2)}, \qquad\text{i.e.}\qquad
\frac{2Q}{k(k+2)} < \frac1{\rho_k}\frac VW. \tag{7.1}$$
By Lemmas B and T, and then (7.1),
$$\mathrm{num} \ge \frac{V^2}W - \frac{2QV}{k(k+2)} - \frac{4\gamma_k}{k+4}
V^{3/2} \;>\; \Bigl(1 - \frac1{\rho_k}\Bigr)\frac{V^2}W -
\frac{4\gamma_k}{k+4}V^{3/2} = \alpha_k\frac{V^2}W - \frac{4\gamma_k}{k+4}
V^{3/2}. \tag{7.2}$$
Next, $V > \Delta/c_k \ge \tfrac{4\theta}{c_kk(k+2)}WQ$ and $Q \ge W/\omega_k$
(from $W \le \omega_kQ$) give $V \ge \tfrac{4\theta}{c_kk(k+2)\omega_k}W^2$;
with $\omega_k \le \tfrac{\theta\alpha_k^2(k+4)^2}{16\gamma_k^2c_kk(k+2)}$ this
yields $V \ge \tfrac{64\gamma_k^2}{\alpha_k^2(k+4)^2}W^2$, equivalently
$\tfrac{4\gamma_k}{k+4}V^{3/2} \le \tfrac{\alpha_k}2\tfrac{V^2}W$. Hence by
(7.2), $\mathrm{num} > \tfrac{\alpha_k}2\tfrac{V^2}W > 0$. Also $\mathrm{den}
\le \tfrac{V^2}3 + \tfrac23P_4 \le V^2$ (the two subtracted terms in
$\mathrm{den}$ are $\ge 0$ since $P_4 \ge V^2/k$; and $\mathrm{den} > 0$ else
$\mathrm{num} = 0$). Lemma Y then gives
$$H \ge \frac{\mathrm{num}^2}{\mathrm{den}} \ge \frac{(\tfrac{\alpha_k}2
\tfrac{V^2}W)^2}{V^2} = \frac{\alpha_k^2}4\frac{V^2}{W^2}.$$
Finally $V > \tfrac{4\theta}{c_kk(k+2)\omega_k}W^2 \ge \tfrac{4c_k}{\alpha_k^2}
W^2$ (using $\omega_k \le \tfrac{\theta\alpha_k^2}{c_k^2k(k+2)}$), i.e.
$\tfrac{\alpha_k^2}4\tfrac V{W^2} \ge c_k$, so $H \ge c_kV$ and $2\Sigma = H +
\Delta - c_kV \ge \Delta \ge 0$. $\blacksquare$

**Remark 7.1 ($k=2$ boundary).** $\alpha_2 = 0$: at $k=2$ the flat window
collapses, consistent with Theorem C2's equality manifold containing
*non*-uniform-leverage configurations (the S–P family). The mechanism making
$k \ge 3$ strictly easier near the manifold is precisely $2\theta > c_k$.

*Status: **PROVED** (validity of the window; constants unoptimized).
(`check_c2k_zones.py` block E: on 1808 near-uniform window samples, every
chain piece holds — $|T_3| - \gamma_kV^{3/2} \le -1.6\cdot10^{-16}$, $B_q -
V^2/W \ge +4.5\cdot10^{-12}$ — and the proved branch margin is $\ge
+9.7\cdot10^{-17}$, no violation.)*

---

## 8. The middle zone: OPEN

What remains is the **middle zone**
$$\frac1{1-\omega_k} \;<\; x = \frac{nQ}{k^2}, \qquad Q \;<\; Q^*(n,k)$$
— leverage profiles that are neither near-uniform nor mass-dominated
(block-cluster territory, $x \approx 1.1$–$1.5$ at moderate $Q$). Its status
is **OPEN**.

- **Obstruction (proved by counterexample) — do not retry $\kappa$-only.**
  The pure-$\kappa$ sufficient condition $H + \Delta \ge c_k\kappa$ is
  **FALSE** in the middle zone: adversarial minimization finds
  $H + \Delta - c_k\kappa = -8.9\cdot10^{-3}$ at $(12,3)$ ($x \approx 1.17$),
  and $-1.8\cdot10^{-3}$ at $(30,5)$ (`s4_kappa_and_map.py`, part C). In the
  middle, $\kappa$ overestimates $V$ by a factor $\sim2.3$ exactly where $H$
  is small.

- **Candidate — CONJECTURED, not proved.** The Y-grab margin
  $$\boxed{\ m_Y := \frac{\mathrm{num}^2}{\mathrm{den}} + \Delta - c_kV \;\ge\;
  0\ } \qquad\textbf{(CONJECTURED)}$$
  implies C2k by Theorem MI + Lemma Y (since $H \ge \mathrm{num}^2/\mathrm{den}
  \ge c_kV - \Delta$). It survived every family and adversarial minimization
  at all tested $(n,k)$ — the full Haar grid, $j$-heavy block families
  (including the s2 failure region at $(30,5)$), single-/two-spike families,
  and Nelder–Mead at $(6,3),(9,3),(10,4),(12,3),(30,5),(12,6)$ — with worst
  value $-4.5\cdot10^{-9}$ at the $(6,3)$ icosahedral equality point (i.e.
  numerical zero, $m_Y\to0$ there), and $\ge +10^{-5}$ away from equality
  manifolds. **Status: NUMERICALLY SUPPORTED (0 true violations), CONJECTURED
  — not part of the proved package.**

- **Why the middle should be provable (program, not a proof).** The two
  Cauchy–Schwarz steps of Theorem K have exact defect identities; $V\approx
  \kappa$ forces every row to be a near-eigenvector of $S_0$ *and* leverage to
  align with the $\hat\beta$-values — i.e. cluster structure — under which
  $B_q \gtrsim V$, making $\mathrm{num}$ large and the Y-grab strong. A
  quantitative stability version of Theorem K along these defects is the
  missing middle-zone ingredient.

*Status of §8: OPEN; obstruction PROVED (FALSE by counterexample); $m_Y$
CONJECTURED. (Candidate re-checked in `check_c2k_zones.py` block F,
cross-multiplied to avoid the $\mathrm{den}\sim V^2$ epsilon trap: $\mathrm
{num}^2 \ge (c_kV - \Delta)\mathrm{den}$ where $c_kV - \Delta > 0$, min
relative margin $+6.4\cdot10^{-2}$, 0 violations. This does **not** upgrade
$m_Y$ to PROVED.)*

---

## 9. Numeric sanity block

Fresh, deterministic, self-contained study for this document:
`numerics/study-c2k-zones/check_c2k_zones.py` (numpy seed `20260718`; sympy
blocks exact). No import of the round-3 study scripts.

```
$ cd numerics/study-c2k-zones && ../../.venv/bin/python check_c2k_zones.py; echo "exit=$?"
```

Observed (run 2026-07-18, **exit 0**, all PASS; full log
`check_c2k_zones.log`):

| block | claim | observed |
|---|---|---|
| A | Master Identity, both forms, 4000 Haar St(n,k), $k=2..7$ | rel err $\le 1.4\cdot10^{-15}$; $H\ge0$ (min $+6.7\cdot10^{-6}$); $\Sigma\ge0$ (min $+2.9\cdot10^{-3}$) |
| B | $H = \|\pi_4M_4\|^2$ via explicit projection; Y-grab $\mathrm{num}^2\le H\,\mathrm{den}$ | rel err $\le 1.1\cdot10^{-16}$; $\mathrm{num}^2 - H\,\mathrm{den} \le -1.9\cdot10^{-5}$ |
| C | Δ-factorization; frame bound (4.1); $V\le\kappa$; block tightness | rel $\le3.6\cdot10^{-14}$; frame excess $\le-7.7\cdot10^{-4}$; $\min(\kappa-V)=+1.6\cdot10^{-3}$; block $|\kappa-V|=0$ |
| D | Zone S: $Q\ge Q^*\Rightarrow \Delta - c_k\kappa\ge0$, $2\Sigma\ge0$ (2493 spiked) | min $+0.124$ / $+0.694$ |
| E | Zone F: window samples $\Rightarrow$ chain pieces + branch margin $\ge0$ (1808) | $|T_3|-\gamma V^{1.5}\le-1.6\cdot10^{-16}$; $B_q-V^2/W\ge+4.5\cdot10^{-12}$; margin $\ge+9.7\cdot10^{-17}$ |
| F | Candidate $m_Y\ge0$ (CONJECTURED), $\mathrm{num}^2\ge(c_kV-\Delta)\mathrm{den}$ | min rel margin $+6.4\cdot10^{-2}$, 0 violations |
| G | Exact (sympy): proj constants, $c_k$, Δ-factor, $Q^*_\infty$ closed form, MI+crux on exact Stiefel, icosahedral ETF | all exact zeros |

**Relative-vs-absolute epsilon gotcha (a6 ledger).** Because $\mathrm{den}
\sim V^2$ and $V$ can be legitimately $\sim10^{-8}$, this study never divides
by $\mathrm{den}$: the Y-grab is checked as $\mathrm{num}^2 \le H\,\mathrm{den}$
(block B) and the candidate as $\mathrm{num}^2 \ge (c_kV-\Delta)\mathrm{den}$
(block F), both with $\mathrm{den}\ge0$, and inequality margins are read
relative to a $V^4$ scale. An absolute guard on $\mathrm{den}$ silently zeroes
the grab and fakes a negative margin (round-3 s7 diagnosis).

The original round-3 scripts also still reproduce in the current workspace
(`s5_exact_constants.py` exit 0, all exact; `s6_zones_and_witnesses.py` exit 0
after the s7 den-guard fix).

---

## 10. Deliverables and status lines

- **Theorem MI + Corollary MI:** $2\Sigma(A) = H + \Delta(Q) - c_kV$ for every
  $A\in\mathrm{St}(n,k)$, with $H\ge0$ the spin-4 energy, $\Delta\ge0$ the
  $Q$-gap factor, $c_k = (k+1)(k^2-2k+8)/(k^2(k+4))$; crux C2k $\iff H+\Delta
  \ge c_kV$; equality manifold = uniform leverage + spin-4-free (icosahedral
  ETF exact). **Status: PROVED (unconditional).**
- **Lemma Δ:** factorization and split (3.1). **Status: PROVED.**
- **Theorem K (+ Lemma 4.1):** $V\le\kappa=\sum\ell^3 - Q^2/k$ via the exact
  frame bound $T_\nu\preceq\mathrm{Id}$; tight at orthonormal blocks. **Status:
  PROVED.**
- **Proposition S:** C2k for all $Q\ge Q^*(n,k)$, closed form $Q^*_\infty =
  (k^3+4k+16)/(k^2+4k+16)$. **Status: PROVED (unconditional, every $n$).**
- **Proposition F:** C2k in the flat window $x\le1/(1-\omega_k)$, $k\ge3$
  ($\omega_3 = 5/242$; constants unoptimized). **Status: PROVED (validity of
  window).**
- **Middle zone:** OPEN. $H+\Delta\ge c_k\kappa$ **FALSE** (counterexample);
  $m_Y\ge0$ **CONJECTURED** (numerically supported, 0 violations). **Status:
  OPEN; $m_Y$ not proved.**
- **Not claimed:** C2k in full (middle zone open); GTZ$(n,k)$ (C2k is the
  weighted-average pair form, not a good-subset statement — Scope Guard 3).
  The $k=2$ anchor Theorem C2 is proved separately in `proofs/theorem-c2.md`
  (VAL-THM-C2-001).
