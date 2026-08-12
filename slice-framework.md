# The $(6,3)$ slice framework: involution reformulation, pinned block spectra, $K_4$ coherence, and two obstructions

**Deliverable:** standalone, referee-grade write-up of the PROVED exact results of
attack round 3 (`attacks/round-3-diffuse-63.md`; numerics
`numerics/study-attack-r3-diffuse-63/`), target **VAL-THM-SLICE-001**. This is the
structural framework for the *diffuse regime* of $\mathrm{GTZ}(6,3)$.

---

## 0. Scope — read this first (overclaim guard)

> **This document does NOT prove $\mathrm{GTZ}(6,3)$.** It develops the exact structure
> of the *equal-leverage slice* $\ell_i \equiv \tfrac12$ (the model of the deep-diffuse
> regime, which contains the icosahedral obstruction), proves pinned spectra of $(6,3)$
> involutions and the $K_4$ coherence/hot-graph/Ramsey lemmas, and records two proved
> obstructions that kill specific proof strategies. It is a *framework*, not a solution.

**What is PROVED here** (exact, adversarially checkable):

| # | Statement | §  | Status |
|---|---|---|---|
| 1 | Slice $\ell_i\equiv\tfrac12 \Leftrightarrow$ zero-diagonal symmetric orthogonal involution $J=2P-I$; row/flow identities | §2 | PROVED |
| 2 | At-most-one-below lemma; goodness $\Leftrightarrow \det(J_{TT}+\tfrac23I)\ge0 \Leftrightarrow \tau_T\ge\tfrac2{27}$ (Veronese); $\sum_T\tau_T=\tfrac49$ | §3 | PROVED |
| 3 | Pinned $5\times5$ spectrum $\{1,1,-d_v,-1,-1\}$ (deleted-column eigenvector) for ANY $(6,3)$ involution | §4.1 | PROVED |
| 4 | Pinned $4\times4$ spectrum $\{1,\beta_2,\beta_3,-1\}$ (explicit symmetric functions); slice $\{1,t,-t,-1\}$, $t=|a_{F^c}|$ | §4.2 | PROVED |
| 5 | $K_4$ coherence balance $\sum_{T\subset F}q_T=0$ for every $4$-set $F$; coherence vector in a $5$-dim space | §4.3 | PROVED |
| 6 | Hot-pair lemma (kills 8 triples), no hot triangle, hot graph $\le6$ edges deg $\le2$; Ramsey lemma | §5 | PROVED |
| 7 | Obstruction A: min-pair-$K_4$ selection rule dead (exact two-ONB witness, all-bad $K_4$ at $\lambda_{\min}=-1/\sqrt2$) | §6 | PROVED |
| 8 | Obstruction B: unsigned $\tau$-averaging (global and $K_4$-local) cannot certify a good triple | §7.1 | PROVED |
| 9 | PMC (pentagonal matching configuration) is an exact slice point; $G_{\mathrm{PMC}}=$ root of the stated cubic | §8 | PROVED |

**What is CONJECTURED / NUMERICALLY SUPPORTED** (must NOT be read as proved):

| Statement | § | Status |
|---|---|---|
| **No-uniform-gap:** $\inf\{f(A):A$ in the open diffuse regime$\}=1/\sqrt6$ (not attained) | §7.2 | **NUMERICALLY SUPPORTED** |
| **$C_{\text{slice}}$:** PMC is the slice minimizer; $\min_{\text{slice}}G=G_{\mathrm{PMC}}$, margin $\approx0.0394$ | §8 | **CONJECTURED / NUMERICALLY SUPPORTED** |
| SOS/SDP certificate for the slice theorem | §7.3 | **NO CERTIFICATE FOUND** at degree $\le6$ (level 2 $=-7/135$, level 3 $\approx-0.047$) |

**Dependencies.** `proofs/duality.md` Theorem 1 / Corollary 1 (VAL-DUAL-001): $(6,3)$ is
self-dual ($n-k=3=k$), and complementation of triples is the slice involution $J\mapsto-J$
(§3, N10). `proofs/extension-lemma.md` (VAL-THM-EXT-001): the h-branch context — the slice
lies *entirely* in the diffuse regime (§1).

---

## 1. Setting, notation, and the diffuse context

Throughout $A\in\mathrm{St}(6,3)$ ($A^TA=I_3$), rows $r_1,\dots,r_6\in\mathbb R^3$,
$\ell_i=\|r_i\|^2$, $c_{ij}=\langle r_i,r_j\rangle$, $P=AA^T$ the rank-$3$ orthoprojector
($P^2=P$, $\operatorname{tr}P=3$, $P_{ii}=\ell_i$, $P_{ij}=c_{ij}$). A $3$-subset (triple)
$T$ is **$1/6$-good** iff $\lambda_{\min}(P_{TT})\ge\tfrac16$ (equivalently
$\sigma_{\min}(A_T)\ge1/\sqrt6$). $\mathrm{GTZ}(6,3)$ asserts every $A$ has a good triple.

We use the notation of `proofs/extension-lemma.md`: for a pair with Gram eigenvalues
$\mu_1\ge\mu_2$, set $h(\mu)=(1-\mu)/(6\mu-1)$; a pair **qualifies** iff $\mu_2>1/6$, and
its **excess** is $h(\mu_1)+h(\mu_2)-\tfrac13\ (\ge0)$. The Extension Lemma
(`extension-lemma.md` §3) closes the *h-branch* excess $\le0$; the **diffuse regime** is
"every qualifying pair has strictly positive excess," where $\mathrm{GTZ}(6,3)$ remains open.

**The slice is entirely diffuse [PROVED].** The polynomial form of the h-branch
(`extension-lemma.md`; verified here in the numeric block, s0 L1) is, for $\mu_2>1/6$,
$$h(\mu_1)+h(\mu_2)\le\tfrac13 \iff R_6(s,g):=54\,s-144\,g-14\le0,\qquad
(h_1+h_2-\tfrac13)\cdot6(6\mu_1-1)(6\mu_2-1)=R_6,$$
with $s=\mu_1+\mu_2=\ell_i+\ell_j$ and $g=\mu_1\mu_2=\ell_i\ell_j-c_{ij}^2$. For an
equal-leverage pair ($\ell_i=\ell_j=\tfrac12$, so $s=1$, $g=\tfrac14-c_{ij}^2$):
$$R_6=54-144(\tfrac14-c_{ij}^2)-14=4+144\,c_{ij}^2>0.$$
Every equal-leverage pair is strictly off the h-branch: **the slice is a model of the deep
diffuse world**, and (§8) it contains the icosahedral ETF. **Status: PROVED.**

---

## 2. The slice is the stratum of zero-diagonal symmetric orthogonal involutions

> **Lemma 1 (slice $\Leftrightarrow$ involution).** Let $J:=2P-I_6$. Then
> $$\ell_i\equiv\tfrac12 \iff \operatorname{diag}(J)=0.$$
> On the slice, $J$ is symmetric, $J^2=I_6$, and $\operatorname{spec}(J)=(+1)^3(-1)^3$.
> Conversely, every zero-diagonal symmetric orthogonal $J\in\mathbb R^{6\times6}$ yields
> $P=(I+J)/2$, a rank-$3$ orthoprojector with $P_{ii}=\tfrac12$, i.e. a slice point (unique
> mod $O(3)$). Writing $a_{ij}:=J_{ij}$ and $u_i:=\sqrt2\,r_i$, the $u_i$ are unit vectors,
> $\sum_i u_iu_i^T=2I_3$, $a_{ij}=\langle u_i,u_j\rangle$, and for all $i$ and all $i\ne j$
> $$\textstyle\sum_{k\ne i}a_{ik}^2=1\ \ \text{(row mass)},\qquad
> \sum_{m\notin\{i,j\}}a_{im}a_{jm}=0\ \ \text{(flow)}.$$

**Proof.** $J_{ii}=2\ell_i-1$, which is $0$ iff $\ell_i=\tfrac12$; and
$\sum_i J_{ii}=2\operatorname{tr}P-6=0$ always. $J$ is symmetric ($P$ is). Since $P^2=P$,
$J^2=(2P-I)^2=4P^2-4P+I=I$; a symmetric involution has eigenvalues $\pm1$, and
$\operatorname{tr}J=0$ forces multiplicities $(+1)^3(-1)^3$ (three $+1$'s from
$\operatorname{rank}P=3$). Conversely, symmetric orthogonal $J$ with $\operatorname{diag}J=0$
has $\operatorname{tr}J=0$, so $P=(I+J)/2$ is symmetric, idempotent
($P^2=\tfrac14(I+J)^2=\tfrac14(2I+2J)=P$), rank $\tfrac12(6+\operatorname{tr}J)=3$, with
$P_{ii}=\tfrac12$; any factorization $P=AA^T$, $A\in\mathrm{St}(6,3)$, is unique mod $O(3)$.
The identities: $|u_i|^2=2\ell_i=1$; since $\sum_i r_ir_i^T=A^TA=I_3$ we get
$\sum_i u_iu_i^T=2\sum_i r_ir_i^T=2I_3$ (a unit tight frame). Row mass:
$\sum_{k\ne i}a_{ik}^2=4\sum_{k\ne i}c_{ik}^2=4\big((P^2)_{ii}-\ell_i^2\big)
=4(\ell_i-\ell_i^2)=4(\tfrac12-\tfrac14)=1$ (using $P^2=P$). Flow: the general identity
$\sum_{m\notin\{i,j\}}c_{im}c_{jm}=(P^2)_{ij}-P_{ii}P_{ij}-P_{ij}P_{jj}+P_{ij}\cdot0
=c_{ij}(1-\ell_i-\ell_j)$ (s0 L5) equals $0$ on the slice; multiply by $4$. **Status: PROVED**
(s2 N1 confirms $J^2=I$, $\operatorname{diag}J=0$ to $\le10^{-12}$ on 200 slice samples). $\square$

Everything below depends only on the six **lines** $\pm u_i$ and the Gram entries $a_{ij}$.

---

## 3. Goodness criterion and the Veronese form

For a triple $T=\{i,j,k\}$ write its three $J$-entries $a=a_{ij}$, $b=a_{ik}$, $c=a_{jk}$ and
$$p_T:=a^2+b^2+c^2,\qquad q_T:=abc.$$
The zero-diagonal block $J_{TT}=\begin{pmatrix}0&a&b\\a&0&c\\b&c&0\end{pmatrix}$ has
$\operatorname{charpoly}(x)=x^3-p_Tx-2q_T$ and $\det(J_{TT}+\tfrac23I)=\tfrac8{27}-\tfrac23p_T+2q_T$
(both exact, s0 L3).

> **Lemma 2 (at-most-one-below).** Any principal block $J_{TT}$ of a slice involution has
> $-I\preceq J_{TT}\preceq I$ and $\operatorname{tr}J_{TT}=0$; hence **at most one** eigenvalue is
> $<-\tfrac23$, and $\lambda_{\mathrm{mid}}=-\tfrac23$ together with $\lambda_{\min}<-\tfrac23$
> is impossible. Consequently
> $$T\text{ good}\iff\lambda_{\min}(J_{TT})\ge-\tfrac23\iff\det(J_{TT}+\tfrac23I)\ge0
> \iff q_T\ge\tfrac{p_T}3-\tfrac4{27}.$$

**Proof.** $J$ is a symmetric involution, so $-I\preceq J\preceq I$, inherited by principal
submatrices; all three eigenvalues lie in $[-1,1]$ and sum to $0$. If two were $\le-\tfrac23$
their sum is $\le-\tfrac43$, forcing the third $\ge\tfrac43>1$ — impossible; the borderline
case $\lambda_{\min}<-\tfrac23,\ \lambda_{\mathrm{mid}}=-\tfrac23$ gives sum $<-\tfrac43$,
same contradiction. So among the three factors $(\lambda_r+\tfrac23)$ at most one is negative;
their product $\det(J_{TT}+\tfrac23I)$ is $\ge0$ iff none is negative iff
$\lambda_{\min}\ge-\tfrac23$. Finally goodness: $P_{TT}=\tfrac12(I+J_{TT})$ so
$\lambda_{\min}(P_{TT})=\tfrac12(1+\lambda_{\min}(J_{TT}))\ge\tfrac16\Leftrightarrow
\lambda_{\min}(J_{TT})\ge-\tfrac23$. The algebraic form follows from
$\det(J_{TT}+\tfrac23I)=\tfrac8{27}-\tfrac23p_T+2q_T\ge0$. **Status: PROVED**
(s0 L7 / s2 N2, 0 violations). $\square$

*(Off the slice the same argument needs $e_1(P_{TT})\ge\tfrac43$: two eigenvalues of $P_{TT}$
below $1/6$ force the third above $1$, impossible since $P_{TT}\preceq I$; then good
$\Leftrightarrow e_3(P_{TT}-\tfrac16I)\ge0$. Verified on $13691$ heavy triples, s0 L4.)*

### 3.1 The Veronese form

Let $w_i:=u_iu_i^T-\tfrac13I_3\in\operatorname{Sym}_0(\mathbb R^3)\cong\mathbb R^5$. Since the
$u_i$ are unit vectors:
$$w_i^2=\tfrac13w_i+\tfrac29I_3\quad(\text{Veronese quadric}),\qquad
\|w_i\|^2=\operatorname{tr}w_i^2=\tfrac23,\qquad
\langle w_i,w_j\rangle=\operatorname{tr}(w_iw_j)=a_{ij}^2-\tfrac13,$$
and tightness $\sum_iw_i=\sum_iu_iu_i^T-2I=0$. (All exact mod the unit-norm ideal; s2 N8.)

> **Lemma 3 (Veronese goodness criterion and totals).** For every triple,
> $$\tau_T:=\operatorname{tr}(w_iw_jw_k)=q_T-\tfrac{p_T}3+\tfrac29
> \quad\Longrightarrow\quad\boxed{\,T\text{ good}\iff\tau_T\ge\tfrac2{27}.\,}$$
> Moreover, over the $20$ triples,
> $$\sum_T\tau_T=\tfrac49,\qquad \sum_T p_T=12,\qquad \sum_T q_T=0.$$

**Proof.** The trace identity $\operatorname{tr}(w_iw_jw_k)=q_T-\tfrac{p_T}3+\tfrac29$ is an
exact polynomial identity in the unit vectors modulo $\{|u|^2-1\}$ (s2 N6, remainder $0$;
note for symmetric matrices all six trace orderings agree, since
$\operatorname{tr}(ABC)=\operatorname{tr}((ABC)^T)=\operatorname{tr}(CBA)$). The criterion:
$\tau_T\ge\tfrac2{27}\Leftrightarrow q_T-\tfrac{p_T}3\ge\tfrac2{27}-\tfrac29=-\tfrac4{27}
\Leftrightarrow q_T\ge\tfrac{p_T}3-\tfrac4{27}$, which is Lemma 2's criterion.

$\sum_T p_T=12$: each of the $15$ pairs $\{i,j\}$ lies in exactly $4$ triples, so
$\sum_Tp_T=4\sum_{i<j}a_{ij}^2=4\cdot3=12$ (total off-diagonal mass
$\sum_{i<j}a_{ij}^2=\tfrac12\sum_i(\text{row mass})=\tfrac12\cdot6=3$).

$\sum_T\tau_T=\tfrac49$: expand $0=\operatorname{tr}\big((\sum_iw_i)^3\big)$ by the multiset
type of the ordered index triple $(i,j,k)$:
$$0=\underbrace{\sum_i\operatorname{tr}(w_i^3)}_{\text{all equal}}
+3\underbrace{\sum_{i\ne j}\operatorname{tr}(w_i^2w_j)}_{\text{two equal}}
+6\underbrace{\sum_{T}\tau_T}_{\text{all distinct}} .$$
Using $w_i^2=\tfrac13w_i+\tfrac29I$: $\operatorname{tr}(w_i^3)=\operatorname{tr}(w_iw_i^2)
=\tfrac13\operatorname{tr}(w_i^2)+\tfrac29\operatorname{tr}(w_i)=\tfrac13\cdot\tfrac23=\tfrac29$,
so the first sum is $6\cdot\tfrac29=\tfrac43$. Next
$\operatorname{tr}(w_i^2w_j)=\tfrac13\operatorname{tr}(w_iw_j)+\tfrac29\operatorname{tr}(w_j)
=\tfrac13\langle w_i,w_j\rangle$, so
$\sum_{i\ne j}\operatorname{tr}(w_i^2w_j)=\tfrac13\big(\operatorname{tr}((\sum w)^2)-\sum_i\|w_i\|^2\big)
=\tfrac13(0-6\cdot\tfrac23)=-\tfrac43.$
Hence $0=\tfrac43+3(-\tfrac43)+6\sum_T\tau_T=-\tfrac83+6\sum_T\tau_T$, giving
$\sum_T\tau_T=\tfrac49$. Then from $\tau_T=q_T-\tfrac{p_T}3+\tfrac29$:
$\tfrac49=\sum_Tq_T-\tfrac13\cdot12+20\cdot\tfrac29=\sum_Tq_T-4+\tfrac{40}9$, so
$\sum_Tq_T=0$. **Status: PROVED** (s2 N6–N8, N11, all exact/$\le10^{-12}$). $\square$

The mean $\bar\tau=\tfrac49/20=\tfrac1{45}<\tfrac2{27}$: unsigned averaging fails by the factor
$\tfrac{2/27}{1/45}=\tfrac{10}3$ (obstruction B, §7.1). Slice $\mathrm{GTZ}(6,3)$ becomes:
*six points on the Veronese surface in $S^4$ with barycenter $0$ always contain a triple with
$\operatorname{tr}(w_iw_jw_k)\ge\tfrac2{27}$.*

**Slice duality [PROVED, cites VAL-DUAL-001].** At $(6,3)$, $n-k=3=k$: the configuration is
self-dual (`duality.md` Cor 1: $f(A)=f(B)$, $\mathrm{GTZ}(6,3)\Leftrightarrow\mathrm{GTZ}(6,3)$).
On the slice, complementation $T\mapsto T^c$ acts as $J\mapsto-J$, giving $p_{T^c}=p_T$ and
$q_{T^c}=-q_T$ (s2 N10) — the shadow of the spectral duality $\sigma_{\min}(A_T)=\sigma_{\min}(B_{T^c})$.

---

## 4. Pinned block spectra of $(6,3)$ involutions

Let $J$ be **any** $(6,3)$ involution ($J=2P-I$, symmetric, $J^2=I$,
$\operatorname{spec}=(+1)^3(-1)^3$), diagonal $d_i:=J_{ii}=2\ell_i-1$, $\sum_id_i=0$. §4.1–4.2
do **not** assume the slice.

### 4.1 The $5\times5$ blocks

> **Lemma 4.** For $F=[6]\setminus\{v\}$, $\operatorname{spec}(J_{FF})=\{1,1,-d_v,-1,-1\}$
> exactly, and the middle eigenvector is the deleted column $a^v:=(J_{iv})_{i\ne v}$.

**Proof.** Order $\operatorname{spec}(J)=(1,1,1,-1,-1,-1)=\lambda_1\ge\dots\ge\lambda_6$ and
$\beta_1\ge\dots\ge\beta_5=\operatorname{spec}(J_{FF})$. Cauchy interlacing
($\lambda_i\ge\beta_i\ge\lambda_{i+1}$) pins $\beta_1\in[\lambda_2,\lambda_1]=\{1\}$,
$\beta_2\in[\lambda_3,\lambda_2]=\{1\}$, $\beta_4\in[\lambda_5,\lambda_4]=\{-1\}$,
$\beta_5\in[\lambda_6,\lambda_5]=\{-1\}$. The trace gives
$\beta_3=\operatorname{tr}J_{FF}-2=\big(\sum_{i\ne v}d_i\big)=-d_v$ (since $\sum_id_i=0$). The
eigenvector: for $i\ne v$,
$(J_{FF}a^v)_i=\sum_{j\ne v}J_{ij}J_{jv}=(J^2)_{iv}-J_{iv}J_{vv}=\delta_{iv}-J_{iv}d_v=-d_v\,(a^v)_i$
(using $J^2=I$ and $i\ne v$). **Status: PROVED** (s2 N4, dev $\le10^{-12}$). $\square$

*(Slice: $d_v=0$, so $\operatorname{spec}=\{1,1,0,-1,-1\}$.)*

### 4.2 The $4\times4$ blocks

> **Lemma 5.** For $F=[6]\setminus\{v,w\}$, $\operatorname{spec}(J_{FF})=\{1,\beta_2,\beta_3,-1\}$
> with $$\beta_2+\beta_3=-(d_v+d_w),\qquad \beta_2^2+\beta_3^2=\textstyle\sum_{i\in F}d_i^2+2m_F-2,$$
> $m_F:=\sum_{i<j,\,i,j\in F}J_{ij}^2$. **On the slice** ($d\equiv0$): $m_F=1+a_{vw}^2$ and
> $$\operatorname{spec}(J_{FF})=\{1,\,t,\,-t,\,-1\},\qquad t=|a_{vw}|=|\text{the complementary pair entry}|.$$

**Proof.** Interlacing with a two-step drop ($\lambda_i\ge\beta_i\ge\lambda_{i+2}$) pins
$\beta_1\in[\lambda_3,\lambda_1]=\{1\}$, $\beta_4\in[\lambda_6,\lambda_4]=\{-1\}$. Trace:
$\beta_2+\beta_3=\operatorname{tr}J_{FF}-(1-1)=\sum_{i\in F}d_i=-(d_v+d_w)$. Frobenius:
$1+\beta_2^2+\beta_3^2+1=\|J_{FF}\|_F^2=\sum_{i\in F}d_i^2+2m_F$, giving the second identity.
On the slice $d\equiv0$: $\beta_2+\beta_3=0$ and $\beta_2^2+\beta_3^2=2m_F-2$, so
$\beta_3=-\beta_2$ and $\beta_2^2=m_F-1$. The mass count: total off-diagonal mass is $3$
(Lemma 3 proof); the mass touching $v$ or $w$ is (row-$v$ mass)$+$(row-$w$ mass)$-a_{vw}^2
=1+1-a_{vw}^2$, so $m_F=3-(2-a_{vw}^2)=1+a_{vw}^2$ and $\beta_2^2=a_{vw}^2$, i.e.
$\beta_2=|a_{vw}|=t$. **Status: PROVED** (s2 N3, dev $\le10^{-12}$; s2 N13 gives the exact
icosahedral case $t=1/\sqrt5$). $\square$

*(Equivalent slice form: $J_{FF}^2=I-bb^T-cc^T$, $b,c$ the two deleted columns restricted to
$F$, $\|b\|^2=\|c\|^2=1-a_{vw}^2$, $\langle b,c\rangle=0$ by the flow identity.)*

### 4.3 $K_4$ coherence balance

> **Lemma 6 (slice).** For every $4$-set $F$, $\sum_{T\subset F}q_T=0$. Off the slice,
> $\sum_{T\subset F}\det(J_{TT})=d_v+d_w$ with $\{v,w\}=F^c$.

**Proof.** For a $4\times4$ matrix, the third elementary symmetric function of the eigenvalues
equals the sum of $3\times3$ principal minors: $e_3(J_{FF})=\sum_{T\subset F}\det(J_{TT})$. On
the slice $\det(J_{TT})=2q_T$ (the zero-diagonal $3\times3$ determinant is $2abc$), while for
$\operatorname{spec}=\{1,t,-t,-1\}$,
$$e_3=1\!\cdot\!t\!\cdot\!(-t)+1\!\cdot\!t\!\cdot\!(-1)+1\!\cdot\!(-t)\!\cdot\!(-1)+t\!\cdot\!(-t)\!\cdot\!(-1)
=-t^2-t+t+t^2=0.$$
Hence $2\sum_{T\subset F}q_T=0$. Off the slice, $e_3(\{1,\beta_2,\beta_3,-1\})=-(\beta_2+\beta_3)
=d_v+d_w$ by Lemma 5. **Status: PROVED** (s2 N5, worst $10^{-12}$). $\square$

> **Corollary 7 (dimension of coherence).** The $15$ balance equations $\sum_{T\subset F}q_T=0$
> ($|F|=4$) have rank $15$, so the admissible coherence vectors $q\in\mathbb R^{20}$ lie in a
> $5$-dimensional subspace.

**Proof.** $T\subset F\iff F^c\subset T^c$ identifies the incidence operator ($4$-sets vs the
$3$-sets they contain) with the inclusion matrix $W_{2,3}(6)$ ($2$-subsets in $3$-subsets),
which has full rank $\binom62=15$ (Gottlieb's theorem; valid since $2+3\le6$). So the solution
space has dimension $20-15=5$. **Status: PROVED** (rank of $W_{2,3}$ cited as standard
combinatorial linear algebra; the balance identities themselves proved above). $\square$

---

## 5. Hot/cold structure and the Ramsey lemma

Call a pair $\{i,j\}$ **hot** if $|a_{ij}|>\tfrac23$. Work on the slice.

> **Lemma 8 (hot pair kills 8).** A hot pair $\{i,j\}$ makes every triple containing it and
> every triple avoiding it bad — all $4+4=8$ such triples have $\lambda_{\min}(J_{TT})<-\tfrac23$.

**Proof.** If $T\supseteq\{i,j\}$, interlacing gives
$\lambda_{\min}(J_{TT})\le\lambda_{\min}(J_{\{i,j\}\{i,j\}})=-|a_{ij}|<-\tfrac23$. If
$T\subseteq F:=\{i,j\}^c$ (a $4$-set), by Lemma 5 the $4\times4$ block $J_{FF}$ has middle
eigenvalue $\beta_3=-|a_{ij}|$ (here $t=|a_{F^c}|=|a_{ij}|$), and interlacing of the $3\times3$
$J_{TT}\subset J_{FF}$ gives $\lambda_{\min}(J_{TT})\le\beta_3=-|a_{ij}|<-\tfrac23$.
Good triples therefore meet a hot pair in exactly one vertex. **Status: PROVED** (s2 N9, 42
hot pairs, 0 escapes). $\square$

> **Lemma 9 (no hot triangle; hot graph is sparse).** No three vertices are pairwise hot. The
> hot graph has maximum degree $\le2$ and at most $6$ edges.

**Proof.** Suppose $\{1,2\},\{1,3\},\{2,3\}$ all hot. The flow identity for $\{2,3\}$ is
$a_{12}a_{13}=-\sum_{m\ge4}a_{2m}a_{3m}$ (moving the $m=1$ term to the left). Row masses give
$\sum_{m\ge4}a_{2m}^2=1-a_{12}^2-a_{23}^2<1-\tfrac49-\tfrac49=\tfrac19$ and likewise
$\sum_{m\ge4}a_{3m}^2<\tfrac19$, so by Cauchy–Schwarz
$|a_{12}a_{13}|\le\big(\sum_{m\ge4}a_{2m}^2\big)^{1/2}\big(\sum_{m\ge4}a_{3m}^2\big)^{1/2}<\tfrac19$.
But $|a_{12}a_{13}|>\tfrac49>\tfrac19$ — contradiction. Degree: three hot edges at a vertex give
row mass $\ge3\cdot\tfrac49>1$, impossible, so degree $\le2$. Edge count: each hot edge
contributes $>\tfrac49$ to $\sum_{i<j}a_{ij}^2=3$, so $\#\text{hot}<3/\tfrac49=\tfrac{27}4=6.75$,
hence $\le6$. **Status: PROVED** (§4 of the ledger; consistent with s2 N9). $\square$

> **Lemma 10 (Ramsey).** If the support graph $G$ (edge $ij$ iff $a_{ij}\ne0$) is
> triangle-free, then an exactly-orthogonal (perfect) triple exists: some $T$ has $J_{TT}=0$,
> $\lambda_{\min}=0\ge-\tfrac23$, $\tau_T=\tfrac29$ — good.

**Proof.** Two-color the edges of $K_6$ by membership in $G$. Since $R(3,3)=6$, every such
coloring has a monochromatic triangle; if $G$ is triangle-free the triangle lies in the
complement, i.e. a triple $T$ with all three $a$-entries $0$, whence $J_{TT}=0$. **Status:
PROVED.** $\square$

A hypothetical all-bad configuration therefore has nonzero coherences of both signs (using
$\sum_Tq_T=0$, Lemma 3) and a support graph with triangles.

---

## 6. Obstruction A: the min-pair-$K_4$ selection-rule family is dead [PROVED]

Lemma 6 ($K_4$ balance) suggested the rule **R1**: *in the $K_4$ complementary to the pair of
smallest $|a|$ — where §4.2 pins the middle spectrum to $\pm|a|_{\min}$ and coherence balances
— pick the best triple.* R1 (and every tie-broken variant) is refuted by an exact witness.

> **Proposition 11 (two-ONB witness).** Let rows $\{0,1,5\}$ and $\{2,3,4\}$ each be an
> orthonormal basis of $\mathbb R^3$, related by the orthogonal cross-Gram
> $R=\begin{pmatrix}\tfrac12&\tfrac12&\tfrac1{\sqrt2}\\[1pt]\tfrac12&\tfrac12&-\tfrac1{\sqrt2}\\[1pt]\tfrac1{\sqrt2}&-\tfrac1{\sqrt2}&0\end{pmatrix}$,
> giving the slice point
> $$J=\begin{pmatrix}
> 0&0&\tfrac12&\tfrac12&\tfrac1{\sqrt2}&0\\
> 0&0&\tfrac12&\tfrac12&-\tfrac1{\sqrt2}&0\\
> \tfrac12&\tfrac12&0&0&0&\tfrac1{\sqrt2}\\
> \tfrac12&\tfrac12&0&0&0&-\tfrac1{\sqrt2}\\
> \tfrac1{\sqrt2}&-\tfrac1{\sqrt2}&0&0&0&0\\
> 0&0&\tfrac1{\sqrt2}&-\tfrac1{\sqrt2}&0&0
> \end{pmatrix}.$$
> Then (i) $J^2=I$, $\operatorname{diag}J=0$ (valid slice point); (ii) the pair $\{4,5\}$ has
> $a_{45}=0$ (tied minimal), and **all four** triangles of its complementary $K_4=\{0,1,2,3\}$
> have $\lambda_{\min}=-1/\sqrt2<-\tfrac23$ (the $K_4$ is entirely bad); (iii) good triples exist
> elsewhere: $\{0,1,5\}$ and $\{2,3,4\}$ are exactly orthogonal ($J_{TT}=0$).

**Proof (hand-checkable).** *(Validity.)* $R\in O(3)$: each row has norm$^2=\tfrac14+\tfrac14+\tfrac12=1$
and the three rows are pairwise orthogonal ($\tfrac14+\tfrac14-\tfrac12=0$;
$\tfrac1{2\sqrt2}-\tfrac1{2\sqrt2}=0$). With $\{u_0,u_1,u_5\}$ and $\{u_2,u_3,u_4\}$ orthonormal
frames and cross-Gram $R$, the six unit vectors form a tight frame $\sum u_iu_i^T=I+I=2I_3$, so
$J_{ij}=\langle u_i,u_j\rangle$ (off-diagonal) with $J^2=I$ (Lemma 1 converse; equivalently
$U^TU\cdot U^TU=U^T(UU^T)U=2\,U^TU$ for $U=[u_0\cdots u_5]$, $UU^T=2I$). The displayed $J$ has
the within-ONB blocks zero and the cross block $=R$; diagonal $0$. *(All-bad $K_4$.)* Inside
$\{0,1,2,3\}$: $a_{01}=0$, $a_{23}=0$, $a_{02}=a_{03}=a_{12}=a_{13}=\tfrac12$. Each triangle has
entry-multiset $\{0,\tfrac12,\tfrac12\}$ up to sign, so $p_T=\tfrac12$, $q_T=0$, charpoly
$x^3-\tfrac12x$, eigenvalues $\{0,\pm\tfrac1{\sqrt2}\}$, $\lambda_{\min}=-\tfrac1{\sqrt2}\approx-0.7071
<-\tfrac23$. *(Good triples.)* $\{0,1,5\}$: $a_{01}=a_{05}=a_{15}=0\Rightarrow J_{TT}=0$
(orthonormal frame), likewise $\{2,3,4\}$. **Status: PROVED** (s7 W1: $\|J^2-I\|<10^{-14}$,
$\operatorname{diag}J=0$ exactly, four $K_4$-triangle $\lambda_{\min}=-0.707107$, the two frames
orthogonal). $\square$

**Consequences.** (i) The entire "select inside the coldest-pair $K_4$" family is closed: any
tie-breaking of R1 can be forced onto $\{4,5\}$, whose $K_4$ has no good triple. (ii) Since a
*legitimate* slice point contains an all-bad $K_4$, **no argument local to a single $K_4$**
(even using its pinned spectrum and balance) can prove the slice theorem — at least two $K_4$'s
must be coupled. Adversarial minimization of the tie-windowed R1 margin reaches $-5.99\cdot10^{-2}$
(s3, seed 20260718) — robustly negative, confirming the obstruction beyond the single witness.

---

## 7. Obstruction B, and the SOS status

### 7.1 Unsigned $\tau$-averaging is dead [PROVED]

> **Proposition 12.** No unsigned average of the goodness margin certifies a good triple:
> the global mean and every $K_4$-local mean of $\tau_T$ lie strictly below the threshold
> $\tfrac2{27}$.

**Proof.** *(Global.)* $\sum_T\tau_T=\tfrac49$ (Lemma 3), so $\bar\tau=\tfrac1{45}<\tfrac2{27}$
(short by the exact factor $\tfrac{10}3$). *($K_4$-local.)* For a $4$-set $F$ with $F^c=\{v,w\}$,
Lemma 6 ($\sum_{T\subset F}q_T=0$) and $\sum_{T\subset F}p_T=2m_F=2(1+a_{vw}^2)$ (each of the $6$
pairs in $F$ lies in $2$ triples of $F$) give
$$\sum_{T\subset F}\tau_T=\sum q_T-\tfrac13\sum p_T+4\cdot\tfrac29
=0-\tfrac23(1+a_{vw}^2)+\tfrac89=\tfrac29-\tfrac23a_{vw}^2,$$
so the $K_4$-mean is $\tfrac1{18}-\tfrac16a_{vw}^2\le\tfrac1{18}<\tfrac2{27}$ (short by the
factor $\ge\tfrac43$, even at $a_{vw}=0$). Averaging — global or over any $K_4$ — can never reach
$\tfrac2{27}$; a proof must **locate** the positive tail of the $\tau$-distribution, not average
it. **Status: PROVED** (s7 W3 confirms $\sum_{T\subset F}\tau_T=\tfrac29-\tfrac23a_{F^c}^2$ to
$10^{-10}$; consistent with the a2 obstruction). $\square$

### 7.2 No uniform diffuse gap [NUMERICALLY SUPPORTED — NOT PROVED]

> **Claim (no-uniform-gap), NUMERICALLY SUPPORTED.** $\inf\{f(A):A$ in the *open* diffuse
> regime at $(6,3)\}=1/\sqrt6$, not attained; the infimum is approached at the h-branch
> boundary along perturbations of the validated extremal (Nesterenko) classes.

**Evidence (not a proof).** Catalog-seeded warm-start continuation down the diffuse-$\delta$
ladder from the VAL-NUM-002 extremal class representatives (s4b/s4c) produces *strictly diffuse*
witnesses with GTZ margin $f-1/\sqrt6$ shrinking to $3.0\cdot10^{-4}$ (class-1 chain, every
qualifying pair at excess $\ge1.1\cdot10^{-3}$) and $7.6\cdot10^{-4}$ (class-2, excess
$\ge3.6\cdot10^{-4}$). Each extremal class has its equality pairs at excess $\approx0$ and **all**
other qualifying pairs at excess $\ge0.5$ (enormous slack), so pushing the few equality pairs off
the h-curve is locally unobstructed and $f\to1/\sqrt6$ by continuity. The earlier
"forced-diffuse floor" $+1.3\cdot10^{-2}$ (round 2) is an **optimizer artifact of generic
starts** (reproduced by s4, dissolved by s4c). **Status: NUMERICALLY SUPPORTED** — this is
*not* proved and must not be cited as such; it says any future diffuse proof must be exact
(first-order tight) on the boundary band, and only the deep-diffuse region admits gap methods.

### 7.3 SOS/SDP track: no certificate found

Lasserre relaxation of $\rho=\min_{\text{slice}}\max_T\varphi_T$ ($\varphi_T=\tau_T-\tfrac2{27}$;
true $\rho\approx+0.033$), epigraph form over the $15$ entry variables with row-mass and flow
equations (s5/s5b, cvxpy+SCS):

| level | moment matrix | extras | bound | wall |
|---|---|---|---|---|
| 2 | $153\times153$ | — | $-0.0518528$ | 12 s |
| 2 | $153\times153$ | $+K_4$ balance | $-0.0518728$ | 1 s |
| 3 | $969\times969$ | $+K_4$ eqs $+$ pairwise products | $-0.0474236$ | 264 s |

The level-2 value is **exactly** the averaging bound $-\tfrac7{135}=-0.051\overline{851}$
($\sum_T\varphi_T=-\tfrac{28}{27}$, mean $-\tfrac7{135}$): degree-4 moments carry no sign
information beyond the mean. Level 3 closes only $8.6\%$ of the gap to $0$. **Report: NO
CERTIFICATE FOUND at degree $\le6$** — this is *not* evidence for or against the slice theorem;
a certificate needs level $\ge4$ with $S_6\ltimes$-switching symmetry reduction, deferred.

---

## 8. The pentagonal matching configuration (PMC)

> **Proposition 13 (PMC is an exact slice point) [PROVED].** With
> $s=\sin36^\circ/\sqrt2$ ($s^2=\tfrac{5-\sqrt5}{16}$), $c=\cos36^\circ/\sqrt2$
> ($c^2=\tfrac{3+\sqrt5}{16}$, note $s^2+c^2=\tfrac12$), the coded matrix
> ($1=s,\,2=\tfrac12,\,3=c$, signs as shown)
> $$\operatorname{code}(J_{\mathrm{PMC}})=\begin{pmatrix}
> 0&-1&2&-3&-2&0\\-1&0&3&0&3&1\\2&3&0&-1&0&2\\-3&0&-1&0&-1&3\\-2&3&0&-1&0&-2\\0&1&2&3&-2&0
> \end{pmatrix}$$
> satisfies $J_{\mathrm{PMC}}^2=I$ exactly. It is six lines in $\mathbb R^3$ forming three
> exactly-orthogonal pairs (the three $0$-pairs are a perfect matching), with cross-angles
> $\arccos(2s),60^\circ,\arccos(2c)$. Its critical value $G_{\mathrm{PMC}}=\max_T\lambda_{\min}(P_{TT})$
> is the smallest real root of
> $$256x^3-384x^2+144x-8-(\sqrt2+\sqrt{10})\sqrt{5-\sqrt5}=0,\qquad
> G_{\mathrm{PMC}}=0.206107373853763\ldots,$$
> attained by $6$ active triples.

**Proof.** Each row of the code has one $0$, two $\pm\tfrac12$, one $\pm s$, one $\pm c$, so row
mass $=2\cdot\tfrac14+s^2+c^2=\tfrac12+\tfrac12=1$; $J_{\mathrm{PMC}}^2=I$ is verified exactly in
sympy (s6) and to $<10^{-14}$ in float (s7 W2), and $P=(I+J)/2$ is the corresponding slice
projector. The eigenvalues of the $20$ blocks $P_{TT}$ are computed exactly; the maximum of their
minima is the stated cubic root (s7 W2: computed $G$ vs cubic root agree to $10^{-13}$), with
$6$ triples active. **Status: PROVED** (exact sympy $J^2=I$; s7 W2 confirms $G=$ cubic root,
$6$ active triples). $\square$

> **Conjecture $C_{\text{slice}}$ [CONJECTURED / NUMERICALLY SUPPORTED].**
> $\min_{\text{slice}}G=G_{\mathrm{PMC}}$; i.e. the slice theorem holds with margin
> $G_{\mathrm{PMC}}-\tfrac16\approx0.0394$, and the icosahedron ($G_{\mathrm{ico}}=0.27639$) is
> **not** the slice minimizer.

**Status: CONJECTURED / NUMERICALLY SUPPORTED** — 64-start slice minimization + a 427-point
survey (s3), polished to $2\cdot10^{-8}$ (s6); nothing below $G_{\mathrm{PMC}}$ observed and the
icosahedral start itself descends to the PMC. This is *not* proved; the fat margin over a
$4$-dimensional slice moduli space makes it a certified-numerics-friendly future target.

---

## 9. Numeric sanity block

All checks run from `numerics/study-attack-r3-diffuse-63/` via `../../.venv/bin/python <script>`
(deterministic; seeds in the study README). Re-run in the current workspace on 2026-07-18 —
**all three scripts exit 0, `ALL CHECKS PASSED`.**

**`s0_exact.py` (exact sympy + float, §§1–3):**
- L1 `R_n identity` PASS; `R_6 = 2(27s-72g-7)` PASS (h-branch polynomial form, §1).
- L2 `slice pair R_6 = 4+144c^2 > 0` PASS (slice strictly diffuse, §1).
- L3 charpoly $x^3-p x-2q$ PASS; $\det(M+\tfrac23I)=\tfrac8{27}-\tfrac23p+2q$ PASS (§3).
- L4 parity guard: $13691$ heavy triples, `0 violations` (off-slice criterion, §3).
- L5 coherence flow identity: worst err $4.0\cdot10^{-16}$ (§2).
- L6 ico exact: Stiefel PASS; coherent Schur slack $=9\sqrt5/55-7/33>0$; ico slice margin
  $=\sqrt5/25-7/135=0.037591>0$ (§3).
- L7 slice goodness $\Leftrightarrow\det$-shift $\ge0$: `0 violations` (Lemma 2).

**`s2_identities.py` (exact + 200 slice samples, §§2–4):**
- N1 $J^2=I$, $\operatorname{diag}J=0$: orth $9.4\cdot10^{-15}$, diag $1.1\cdot10^{-12}$ (Lemma 1).
- N2 at-most-one-below: `0 violations` (Lemma 2).
- N3 $4\times4$ spectrum $\{1,t,-t,-1\}$, $t=|a_{F^c}|$: worst dev $1.1\cdot10^{-12}$ (Lemma 5).
- N4 $5\times5$ spectrum $\{1,1,0,-1,-1\}$: worst dev $1.1\cdot10^{-12}$ (Lemma 4).
- N5 $K_4$ balance $\sum_{T\subset F}q_T=0$: worst $1.1\cdot10^{-12}$ (Lemma 6).
- N6 $\tau=q-p/3+2/9$ (mod unit norms): remainder $0$; good $\Leftrightarrow\tau\ge2/27$: `0 mismatches` (Lemma 3).
- N7 $\sum_T\tau_T=4/9$: worst $9.9\cdot10^{-15}$ (Lemma 3).
- N8 Veronese algebra $w^2=w/3+2I/9$, $\langle w_i,w_j\rangle=a^2-1/3$: exact (Lemma 3).
- N9 hot pair kills its 8 triples: 42 hot pairs, all confirmed (Lemma 8).
- N10 $p_{T^c}=p_T$, $q_{T^c}=-q_T$: $3.8\cdot10^{-15}$ / $8.6\cdot10^{-13}$ (duality, §3).
- N11 $\sum p=12$, $\sum q=0$: $3.0\cdot10^{-14}$ / $5.7\cdot10^{-17}$ (Lemma 3).
- N13 ico $4\times4$ charpoly $=(x^2-1)(x^2-t^2)$, $t=1/\sqrt5$ exact (Lemma 5).

**`s7_witnesses.py` (hand witnesses, §§6, 8):**
- W1 two-ONB: $J^2=I$ ($<10^{-14}$), $\operatorname{diag}J=0$ exact; $K_4\{0,1,2,3\}$ all-bad
  $\lambda_{\min}=-0.707107=-1/\sqrt2$; $\{0,1,5\},\{2,3,4\}$ orthogonal (Proposition 11).
- W2 PMC: $J^2=I$; $G=0.206107373853763=$ smallest cubic root (agree $10^{-13}$); $6$ active
  triples; slice margin $0.039440707$ (Proposition 13).
- W3 $K_4$ mean $\sum_{T\subset F}\tau_T=2/9-\tfrac23a_{F^c}^2$: worst $9.9\cdot10^{-12}$ (Prop. 12).

No configuration with $f<1/\sqrt6-10^{-6}$ was encountered in any script (no tripwire event).

---

## 10. Status summary

| Statement | Status |
|---|---|
| §2 Lemma 1: slice $\Leftrightarrow$ zero-diag symmetric orthogonal $J$; row/flow identities | **PROVED** |
| §3 Lemma 2: at-most-one-below; goodness $\Leftrightarrow\det(J_{TT}+\tfrac23I)\ge0$ | **PROVED** |
| §3 Lemma 3: $\tau_T=q_T-p_T/3+2/9$; good $\Leftrightarrow\tau_T\ge2/27$; $\sum\tau=4/9$, $\sum p=12$, $\sum q=0$ | **PROVED** |
| §3 slice duality $p_{T^c}=p_T$, $q_{T^c}=-q_T$ | **PROVED** (cites VAL-DUAL-001) |
| §4.1 Lemma 4: $5\times5$ spectrum $\{1,1,-d_v,-1,-1\}$, deleted-column eigenvector (any involution) | **PROVED** |
| §4.2 Lemma 5: $4\times4$ spectrum $\{1,\beta_2,\beta_3,-1\}$; slice $\{1,t,-t,-1\}$, $t=|a_{F^c}|$ | **PROVED** |
| §4.3 Lemma 6 / Cor 7: $K_4$ balance $\sum_{T\subset F}q_T=0$; $q$-space $5$-dim | **PROVED** |
| §5 Lemmas 8–10: hot pair kills 8; no hot triangle; hot graph deg $\le2$, $\le6$ edges; Ramsey | **PROVED** |
| §6 Prop 11: two-ONB witness — min-pair-$K_4$ rule family dead | **PROVED** |
| §7.1 Prop 12: unsigned $\tau$-averaging (global + $K_4$-local) cannot certify a good triple | **PROVED** |
| §8 Prop 13: PMC is an exact slice point; $G_{\mathrm{PMC}}=$ cubic root, $6$ active triples | **PROVED** |
| §7.2 no-uniform-gap: $\inf_{\text{diffuse}}f=1/\sqrt6$ at the boundary | **NUMERICALLY SUPPORTED** |
| §8 $C_{\text{slice}}$: PMC is the slice minimizer (margin $\approx0.0394$); ico not minimizer | **CONJECTURED / NUMERICALLY SUPPORTED** |
| §7.3 SOS/SDP slice certificate | **NO CERTIFICATE FOUND** (degree $\le6$) |
| $\mathrm{GTZ}(6,3)$ | **NOT proved here** (explicitly out of scope, §0) |

**Bottom line.** This is the proved structural scaffold of the $(6,3)$ diffuse regime, not a
resolution of it. The min-pair-$K_4$ and unsigned-averaging strategies are permanently closed
(§§6–7.1); the boundary band of the diffuse regime has no uniform gap (§7.2, numerical) and
requires an exact (second-order Extension-Lemma-type) mechanism, while the deep-diffuse slice is
a fat-margin, low-dimensional target for certified numerics.
