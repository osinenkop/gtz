# C2k is FALSE: an explicit high-precision Stiefel counterexample

**Target:** VAL-THM-MID-001. Refutes **Conjecture C2k** of
`proofs/c2k-zones.md` (VAL-THM-C2K-001). Promoted from attack round 4
(`attacks/round-4-c2k-middle.md`; numerics
`numerics/study-attack-r4-c2k-middle/`, witness
`c2k_counterexample.npy`). Notation follows `proofs/c2k-zones.md` and
MEMORY.md.

---

## 0. Scope: what this document does and does NOT touch

This document settles **one** conjecture — the general‑$k$ weighted‑average
pair inequality **C2k** — **negatively**, by an explicit counterexample. It
is an OBSTRUCTION result. The three guards below are binding.

> **SCOPE GUARD A — only C2k dies; no GTZ result is affected.** The
> refutation of C2k does **not** touch, weaken, or call into question **any**
> GTZ result. In particular the following remain valid exactly as proved:
> **GTZ** itself; the general‑$k$ **Pair Lemma** (`proofs/pair-lemma.md`,
> Perron‑based *existence* of a good pair); the **Extension Lemma**
> (`proofs/extension-lemma.md`); **duality** (`proofs/duality.md`); **Case A**
> (`proofs/case-a.md`); and **Theorem C2** — the $k=2$ case
> (`proofs/theorem-c2.md`, VAL-THM-C2-001) — which is a *separately proved
> theorem*, not a corollary of C2k. C2k was only the hoped‑for
> *weighted‑average* generalization of the Pair Lemma to a single clean
> averaged inequality at general $k$; its failure closes **that specific
> route** and nothing else (§4).

> **SCOPE GUARD B — Propositions S and F of `c2k-zones.md` are untouched.**
> C2k is refuted *in the OPEN middle zone*. The two PROVED zone theorems of
> `proofs/c2k-zones.md` — **Proposition S** (spiked zone $Q \ge Q^*(n,k)$)
> and **Proposition F** (flat window $x \le 1/(1-\omega_k)$) — assert
> $\Sigma(A) \ge 0$ **only on their zones**, and the witness lies in
> *neither* (§3). They therefore remain valid unconditional theorems; what is
> refuted is the *global* Conjecture C2k, not its two proved zones. Likewise
> **Theorem MI**, **Lemma Δ**, **Theorem K**, and the frame bound of
> `c2k-zones.md` are identities/inequalities that continue to hold — indeed
> Theorem MI holds *as an exact identity at the counterexample* and is one of
> the three ways $2\Sigma < 0$ is computed below.

> **SCOPE GUARD C — C2k was never GTZ.** Per Scope Guard 3 of
> `c2k-zones.md`, C2k is the weighted‑*average* pair inequality $\Sigma(A) \ge
> 0$; even if true it would give only the *existence of one good pair*, never
> a good $k$‑subset, and hence never GTZ$(n,k)$. Its refutation removes a
> hoped‑for shortcut, not a proved bound.

### The conjecture under refutation

Fix $A \in \mathrm{St}(n,k)$ ($A^TA = I_k$, rows $r_i \in \mathbb R^k$,
leverages $\ell_i = \|r_i\|^2$, $c_{ij} = \langle r_i, r_j\rangle$, $Q =
\sum_i\ell_i^2$), and set $s = (k+2)/k^2$, $t = 2k/((k+2)n)$. Define
$$\Sigma(A) := \sum_{i<j}(\ell_i\ell_j - c_{ij}^2)\bigl[s(\ell_i - t)
(\ell_j - t) - c_{ij}^2\bigr].$$

> **Conjecture C2k** (`proofs/c2k-zones.md`, §0). $\Sigma(A) \ge 0$ for every
> $A \in \mathrm{St}(n,k)$, $k \ge 2$.

By the Master Identity (Theorem MI of `c2k-zones.md`, PROVED), $2\Sigma = H +
\Delta(Q) - c_kV$ with $H \ge 0$ the spin‑4 energy, $\Delta \ge 0$, $V =
\|S_0\|_F^2 \ge 0$, and $c_k = (k+1)(k^2-2k+8)/(k^2(k+4)) > 0$; so C2k $\iff H
+ \Delta \ge c_kV$.

### Status table

| statement | status |
|---|---|
| **Theorem 1** (C2k is FALSE for $k \ge 4$; explicit witness $(259,5)$) | **PROVED (by counterexample)** |
| Second witness $(244,4)$, $2\Sigma = -1.6\cdot10^{-4}$ | PROVED (by counterexample) |
| **Cor 2** (candidates $m_Y \ge 0$ and $m_Y^{(2)} \ge 0$ are FALSE) | **PROVED** |
| **Cor 3** (the C2k middle zone is UNCLOSABLE) | **PROVED** |
| **Prop 4** (witness is in the OPEN middle; Props S, F untouched) | PROVED |
| Remark 5 (sampling / sign‑masking pitfall) | methodological finding |
| Scope guards A–C (no GTZ / Pair‑Lemma / C2 / S,F damage) | PROVED |

---

## 1. Theorem 1: C2k is false

> **Theorem 1.** Conjecture C2k is **FALSE** for $k \ge 4$. Explicitly, the
> Stiefel matrix $A^\star \in \mathrm{St}(259,5)$ saved to
> `numerics/study-attack-r4-c2k-middle/c2k_counterexample.npy` satisfies
> $$\Sigma(A^\star) < 0, \qquad 2\Sigma(A^\star) = -1.234\,325\,513\cdot
> 10^{-3}.$$

**Why a single explicit witness is a complete proof.** C2k is a *universally
quantified* statement ("$\Sigma(A) \ge 0$ for **every** $A \in
\mathrm{St}(n,k)$"). Its negation is *existential*: one genuine Stiefel matrix
with $\Sigma < 0$ disproves it outright. There is no averaging, no asymptotics,
no numerical extrapolation — only the verification that (i) $A^\star$ is a
Stiefel matrix and (ii) $\Sigma(A^\star) < 0$, both to a precision that
excludes floating‑point artifact. We now certify (i)–(ii).

### (i) $A^\star$ is a genuine Stiefel matrix

$A^\star$ is $259 \times 5$ with
$$\|A^{\star T}A^\star - I_5\|_{\max} = 3.3\cdot10^{-16}\ \text{(float64,
entrywise max)}, \qquad \sum_i\ell_i = 5.0 = k \ \text{(exact)}.$$
(The spectral norm of the same defect is $7.3\cdot10^{-16}$; both are at the
float64 unit‑roundoff floor.) To rule out that the sign of $\Sigma$ is an
artifact of a $\sim10^{-16}$ orthonormality defect, $A^\star$ is
**re‑orthonormalized in 80‑digit arithmetic** (high‑precision Gram–Schmidt on
its columns); the rebuilt matrix has
$$\|A^{\star T}A^\star - I_5\|_{\max} = 4.2\cdot10^{-81} \quad (80\ \text{dps}),
\qquad \sum_i\ell_i = 5.0,$$
i.e. it is Stiefel to full 80‑digit precision, and $\Sigma$ is recomputed on
*this* exact‑Stiefel matrix in (ii). The leverage profile is a **single
dominant cluster** ($\ell_{\max} = 0.539$) plus a light tail.

### (ii) $\Sigma(A^\star) < 0$, computed three independent ways

The value $2\Sigma(A^\star) = -1.234\,325\,513\cdot10^{-3}$ is obtained by
three mutually independent computations, all agreeing to the displayed digits:

1. **Direct pair sum** — evaluate $2\Sigma = \sum_{i\ne j}(\ell_i\ell_j -
   c_{ij}^2)[s(\ell_i-t)(\ell_j-t) - c_{ij}^2]$ straight from the definition
   of $\Sigma$, using only the Gram matrix $G = A^\star A^{\star T}$. **This
   path uses none of the isotypic machinery**, so the refutation does not
   depend on Theorem MI being correct.
2. **Master‑Identity form** — evaluate $H + \Delta(Q) - c_kV$ (Theorem MI of
   `c2k-zones.md`) from $H = \Phi_4 - 3Q^2/(k(k+2)) - 6V/(k+4)$, $\Delta$, and
   $V$. Agreement with path 1 to $11$ digits *simultaneously* re‑confirms
   Theorem MI as an identity **and** the sign of $\Sigma$.
3. **80‑dps exact‑Stiefel rebuild** — recompute the direct pair sum on the
   80‑digit re‑orthonormalized matrix of (i). Result:
   $2\Sigma = -0.001\,234\,325\,513\,35$.

The violation is $\bigl|2\Sigma\bigr| = 1.2\cdot10^{-3}$, about $10^{13}\times$
the float64 unit roundoff and $\sim10^{78}\times$ the 80‑dps Stiefel defect.
There is no catastrophic cancellation (the summands are $O(10^{-1})$ and the
result $O(10^{-3})$; the 80‑dps sum lands on the same value). Hence the sign is
real: **$\Sigma(A^\star) < 0$, and C2k is false.** $\blacksquare$

*Reproducible in a fresh process directly from the saved witness:*
```
cd numerics/study-attack-r4-c2k-middle
../../.venv/bin/python -c "import numpy as np; A=np.load('c2k_counterexample.npy'); \
n,k=A.shape; G=A@A.T; l=np.diag(G).copy(); s=(k+2)/k**2; t=2*k/((k+2)*n); \
M=(np.outer(l,l)-G**2)*(s*np.outer(l-t,l-t)-G**2); \
print('2Sigma =', float(np.sum(M)-np.sum(np.diag(M))))"   # -> -1.234e-3 < 0
```

### Second witness and the failure region

A second, independent counterexample sits at $(n,k) = (244,4)$:
$$2\Sigma = -1.5695\cdot10^{-4} < 0 \qquad (x = 10.8,\ \ell_{\max} = 0.573),$$
reproduced by `s11_c2k_false.py` (seed `20260718`, deterministic; the deepest
per‑$k$ hit). Both violations share the geometry that the earlier searches
never entered: **one dominant leverage cluster** ($\ell_{\max} \approx 0.5$–
$0.57$) **+ a light tail**, at **large $n$** ($\approx 240$–$260$) and **large
$x = nQ/k^2$** ($\approx 7$–$11$), just below the spiked threshold. In the same
deterministic search, $k = 3, 6, 7, 8$ stayed nonnegative (smallest
$+1.66\cdot10^{-4}$ at $k = 3$, near the boundary); a violation at $k = 4$ and
$k = 5$ already falsifies C2k, which is claimed for all $k \ge 2$. **$k = 2$ is
unaffected: Theorem C2 (`proofs/theorem-c2.md`) is a separate PROVED theorem
and is not a special case of anything refuted here.**

*Status: **PROVED (by explicit counterexample).***

---

## 2. Consequences: both middle-zone candidates are false; the zone is unclosable

> **Corollary 2 (both candidates are FALSE).** The two middle‑zone candidates
> of the C2k program are false:
> - the a6 single‑grab candidate $m_Y := \mathrm{num}^2/\mathrm{den} + \Delta
>   - c_kV \ge 0$ (Scope Guard 2 / §8 of `c2k-zones.md`), and
> - this round's corrected 2‑plane candidate $m_Y^{(2)} := H_2 + \Delta - c_kV
>   \ge 0$, where $H_2 = [b_1\,b_2]\,\mathrm{Gram}^{-1}[b_1\,b_2]^T \le H$ is
>   the harmonic energy recovered by the closed‑form 2‑dimensional grab over
>   $\mathrm{span}\{\pi_4Y,\ \pi_4 Z_1\}$, $Z_1 = \sum_i\ell_i r_i^{\otimes4}$.

*Proof.* Each candidate is a **sufficient condition for C2k**: by Lemma Y of
`c2k-zones.md`, $H \ge \mathrm{num}^2/\mathrm{den}$, and by construction $H \ge
H_2$ (a Gram/least‑squares energy over a $\le 2$‑dimensional subspace of the
spin‑4 space is $\le$ the full spin‑4 energy $H = \|\pi_4M_4\|^2$). Hence
*pointwise*, by Theorem MI,
$$2\Sigma = H + \Delta - c_kV \;\ge\; \frac{\mathrm{num}^2}{\mathrm{den}} +
\Delta - c_kV = m_Y, \qquad
2\Sigma = H + \Delta - c_kV \;\ge\; H_2 + \Delta - c_kV = m_Y^{(2)}.$$
Therefore $m_Y \ge 0 \Rightarrow 2\Sigma \ge 0$ and $m_Y^{(2)} \ge 0
\Rightarrow 2\Sigma \ge 0$: **each candidate, if true everywhere, would imply
C2k.** Contrapositively, at the witness $A^\star$ where $2\Sigma < 0$ we have
$m_Y \le 2\Sigma < 0$ and $m_Y^{(2)} \le 2\Sigma < 0$. Both candidates are
false. $\blacksquare$

**Independent direct refutation of $m_Y$.** Before C2k itself was seen to be
false, the round refuted $m_Y$ directly and more strongly: **$m_Y \ge 0$ is
FALSE for every $k \ge 4$** — verified at $(212,10)$ (`s6_mY_refute.py`, 60‑dps
exact‑Stiefel rebuild): there $\mathrm{num} = +6.8\cdot10^{-3}$, $\mathrm{den}
= +2.1\cdot10^{-3} \ge 0$, deficit $c_kV - \Delta = +2.8\cdot10^{-2}$, and
$\mathrm{num}^2/\mathrm{den} = +2.2\cdot10^{-2} < $ deficit, so $m_Y =
-5.9\cdot10^{-3} < 0$ while, *at that configuration*, $2\Sigma = +0.038$ still
held. The single harmonic Y‑grab there recovers only $\approx 39\%$ of $H$. The
per‑$k$ floor of $m_Y/\text{deficit}$ is $+0.11$ ($k=3$) then negative for $k =
4,\dots,12$ ($-0.07, -0.09, -0.21, -0.17, -0.27, -0.39, -0.44, -0.51, -0.57$).
So $m_Y$ was already dead as a lemma independently of C2k; Corollary 2 adds
that it (and $m_Y^{(2)}$) must be false wherever C2k is.

> **Corollary 3 (the C2k middle zone is UNCLOSABLE).** There is no proof that
> $\Sigma(A) \ge 0$ on the middle zone, because the statement is **false**
> there: $A^\star$ lies in the open middle zone (Prop 4) and has $\Sigma < 0$.

*Proof.* Immediate from Theorem 1 and Prop 4: a would‑be middle‑zone theorem
would assert $\Sigma \ge 0$ for all middle‑zone $A$, contradicted by $A^\star$.
Consequently the a6 research program of "prove any sufficient condition
($m_Y$, $m_Y^{(2)}$, $\kappa$‑only) to close the middle" cannot succeed — no
true sufficient condition for a false statement exists on that set.
$\blacksquare$

*Status: **PROVED** (Cor 2, Cor 3).*

---

## 3. The witness lies in the OPEN middle zone

> **Proposition 4.** $A^\star$ is in the **open middle zone** of
> `c2k-zones.md`: it is in **neither** the spiked zone (Proposition S) **nor**
> the flat window (Proposition F). Hence it does **not** contradict either
> proved zone theorem.

*Proof.* The zones of `c2k-zones.md` are: spiked zone $= \{Q \ge Q^*(n,k)\}$
(Prop S), flat window $= \{x \le 1/(1-\omega_k)\}$ (Prop F), and the open
middle $= \{x > 1/(1-\omega_k)\ \text{and}\ Q < Q^*(n,k)\}$. For $A^\star$ at
$(n,k) = (259,5)$, computed from the witness:
$$Q = 0.6679, \qquad Q^*(259,5) = 2.6922, \qquad
x = \frac{nQ}{k^2} = 6.9196, \qquad \frac1{1-\omega_5} = 1.0237.$$
Thus $Q = 0.668 < 2.692 = Q^*$ (so $A^\star \notin$ spiked zone; **Prop S does
not apply and is not contradicted** — it only claims $\Sigma \ge 0$ where $Q
\ge Q^*$), and $x = 6.92 > 1.024 = 1/(1-\omega_5)$ (so $A^\star \notin$ flat
window; **Prop F does not apply and is not contradicted** — it only claims
$\Sigma \ge 0$ where $x \le 1/(1-\omega_k)$). Both defining inequalities of the
open middle hold, so $A^\star$ is an open‑middle configuration. What C2k
asserted, and what is refuted, is $\Sigma \ge 0$ on **all** of $\mathrm{St}
(n,k)$ — in particular on the open middle, which S and F never covered.
$\blacksquare$

*Status: **PROVED**. Proposition S and Proposition F of `c2k-zones.md` remain
valid unconditional theorems; only the global Conjecture C2k is settled
negatively.*

---

## 4. What the refutation does and does not close (scope, expanded)

Because a refutation is easy to over‑read, we state precisely what dies and
what survives.

**Dies.** (a) Conjecture C2k (the global weighted‑average pair inequality
$\Sigma \ge 0$). (b) Both middle‑zone candidates $m_Y$, $m_Y^{(2)}$ (Cor 2),
and the hope of *any* sufficient condition closing the middle (Cor 3). (c) The
specific hoped‑for consequence "C2k $\Rightarrow$ Perron‑free pair existence at
$k \ge 3$" (an a6 analogue of Cor 5.1 of `theorem-c2.md`) — that inference
started from C2k, which is false, so the inference is vacated. (d) The idea
that a single clean *averaged* inequality replaces the Perron argument at
general $k$ (the "no averaged form" barrier is now proved, not conjectured).

**Survives (untouched).** Every genuinely proved result of the mission:
- **GTZ** and its status map — the refutation is about an average, not about
  the existence of a good subset (Scope Guard C).
- The general‑$k$ **Pair Lemma** (`proofs/pair-lemma.md`): its good‑pair
  *existence* is proved by the **Perron–Frobenius** argument on $G = WW^T -
  zz^T$, which never uses C2k. Pair existence is *stronger* information than an
  average sign and is unaffected by the average going negative.
- The **Extension Lemma** (`proofs/extension-lemma.md`), **duality**
  (`proofs/duality.md`), **Case A** (`proofs/case-a.md`), and **Theorem C2**
  at $k=2$ (`proofs/theorem-c2.md`, VAL-THM-C2-001).
- Within `c2k-zones.md`: **Theorem MI**, **Lemma Δ**, **Theorem K**, the
  **frame bound**, and **Propositions S and F** (Scope Guard B, §3). These are
  identities/zone theorems; MI in fact holds *at the counterexample* (§1). They
  are re‑scoped from "the unconditional half of a true conjecture" to "the two
  zones on which the (now‑false) average nevertheless stays $\ge 0$."

The correct redirection of the pair‑existence effort is back to the proved Pair
Lemma line (round 1) and the KKT/active‑set program
(`attacks/round-4-boundary.md`), not to any averaged inequality.

---

## 5. Remark: the sampling / sign-masking pitfall

Every prior "C2k: 0 violations" report — a6 §8, the VAL-THM-C2K-001 numerics,
and this round's own early sweeps (`s1`, `s5`, `s8`) — missed this
counterexample. The reasons are methodological and worth recording, because
they invalidate a whole class of "no‑violation" evidence:

1. **Too small $n$.** The Haar and two‑scale middle maps capped at $n \le
   40$–$46$. At $n \le 42$, `s1_map_middle.py` found $\min 2\Sigma =
   +2.7\cdot10^{-3}$ and $2\Sigma$ *growing* with $x$ — an artifact that made
   C2k look like it had uniform positive slack in the middle. The violation
   requires $n \gtrsim 240$; below that the deficit $c_kV - \Delta$ has not yet
   overtaken $H$.
2. **Single‑cluster geometry is required.** The failure locus is one dominant
   leverage cluster ($\ell_{\max} \approx 1/2$) plus a light tail, at $x
   \approx 7$–$11$ and $k \in \{4,5\}$, just below the spiked threshold $Q^*$.
   Generic Haar leverage profiles (spread, no dominant single cluster) never
   land there; clustered warm starts ($V/\kappa \to 1$) are needed.
3. **Sign masking by a `min` tracker.** An adversarial search keyed on
   $\min(2\Sigma,\, m_Y)$ (as in `s5_adversarial_hunt.py`, designed to hunt
   the $m_Y$ refutation) reports the *more negative* of the two. Since $m_Y \le
   2\Sigma$ always (Cor 2), a small $2\Sigma < 0$ is **hidden behind** a
   more‑negative $m_Y$. A C2k check must track $2\Sigma$ *directly and on its
   own*, never bundled under a `min` with a candidate.

> **Gotcha (binding for any future C2k‑type check).** A "C2k: no violation"
> claim is **meaningless** unless the search reaches $n \gtrsim 260$ with
> single‑cluster leverage geometry ($\ell_{\max} \approx 0.5$, light tail, $x
> \approx 7$–$11$, $k \in \{4,5\}$) and tracks $2\Sigma$ directly. Small‑$n$ /
> Haar / bundled‑`min` searches provably pass a false statement.

Note this does **not** impugn the `numerics/study-c2k-zones/` baseline of
`c2k-zones.md`: that study samples *in‑zone* (spiked / flat / near‑uniform) and
small $n$, so it legitimately does not see the open‑middle counterexample; its
PROVED zone claims (Props S, F) are unaffected. It reruns exit 0 (§6).

---

## 6. Numeric sanity block

All quantities below were recomputed **in a fresh process from the saved
witness** `numerics/study-attack-r4-c2k-middle/c2k_counterexample.npy` (no
import of the search scripts), via the workspace `.venv`.

| check | claim | observed |
|---|---|---|
| shape / Stiefel | $A^\star \in \mathrm{St}(259,5)$ | shape $(259,5)$; $\|A^{\star T}A^\star - I\|_{\max} = 3.33\cdot10^{-16}$ (float), spectral $7.31\cdot10^{-16}$; $\sum\ell = 5.0$ |
| 80‑dps rebuild | exact‑Stiefel after 80‑dps Gram–Schmidt | $\|A^{\star T}A^\star - I\|_{\max} = 4.22\cdot10^{-81}$; $\sum\ell = 5.0$ |
| $2\Sigma$ path 1 | direct pair sum (no MI) | $-1.234\,325\,513\,35\cdot10^{-3}$ |
| $2\Sigma$ path 2 | Master‑Identity form $H+\Delta-c_kV$ | $-1.234\,325\,513\,35\cdot10^{-3}$ (agrees) |
| $2\Sigma$ path 3 | 80‑dps exact‑Stiefel rebuild | $-0.001\,234\,325\,513\,35$ (agrees) |
| zone | $Q, Q^*, x, 1/(1-\omega_5)$ | $Q = 0.6679 < Q^* = 2.6922$; $x = 6.9196 > 1.0237$ ⇒ open middle (inS = inF = False) |
| 2nd witness | $(244,4)$, deterministic $s11$ | $2\Sigma = -1.5695\cdot10^{-4} < 0$, $x = 10.8$ |
| baseline intact | `numerics/study-c2k-zones/check_c2k_zones.py` | exit 0, all PASS (samples in‑zone / small‑$n$; does not see the witness) |

Observed run (fresh process, 2026‑07‑18):
```
(n,k)= (259, 5)
Stiefel defect: max-entry float=3.331e-16   spectral float=7.309e-16   80dps=4.217e-81
sum ell = 5.0  (80dps: 5.0)
2Sig (1) direct  =-1.23432551335e-03
2Sig (2) MI form =-1.23432551335e-03
2Sig (3) 80dps   = -0.00123432551335
Q=0.667914  Q*=2.6922  x=6.9196  flat-thr=1.0237
inS(Q>=Q*)= False   inF(x<=thr)= False
```

The deterministic search that located and saved the witness is
`s11_c2k_false.py` (seed `20260718`; its log `s11.log` shows the per‑$k$ hits,
the 80‑dps confirmation, and the zone check). The saved‑witness reproduction
command of §1 and the fresh‑process table above are the validator entry points.

**Relative‑vs‑absolute epsilon note.** The sign of $2\Sigma$ here is robust
(magnitude $\sim10^{-3}$, $\sim10^{13}\times$ machine roundoff), so no epsilon
subtlety arises for Theorem 1. The `num`/`den`/$H_2$ epsilon trap (never divide
by $\mathrm{den} \sim V^2$) matters only for the *candidate* margins of §2 and
is handled there by cross‑multiplication, exactly as in `c2k-zones.md` §9.

---

## 7. Deliverables and status lines

- **Theorem 1 — C2k is FALSE for $k \ge 4$.** Explicit Stiefel witness
  $(259,5)$ with $2\Sigma = -1.234\cdot10^{-3} < 0$, three independent
  computations + 80‑dps exact‑Stiefel rebuild, reproducible from the saved
  `.npy`; second witness $(244,4)$, $2\Sigma = -1.6\cdot10^{-4}$. **Status:
  PROVED (by counterexample).**
- **Corollary 2 — $m_Y$ and $m_Y^{(2)}$ are FALSE.** Each implies C2k
  pointwise ($2\Sigma \ge m_Y$, $2\Sigma \ge m_Y^{(2)}$), so both are $< 0$ at
  the witness; $m_Y$ is independently false for all $k \ge 4$ (verified
  $(212,10)$, 60 dps). **Status: PROVED.**
- **Corollary 3 — the middle zone is UNCLOSABLE.** Its target inequality is
  false there. **Status: PROVED.**
- **Proposition 4 — witness in the OPEN middle zone.** $Q = 0.668 < Q^* =
  2.69$, $x = 6.92 > 1/(1-\omega_5) = 1.024$; Propositions S and F of
  `c2k-zones.md` are not contradicted. **Status: PROVED.**
- **Scope guards A–C — no collateral damage.** GTZ, the general‑$k$ Pair
  Lemma, the Extension Lemma, duality, Case A, Theorem C2 ($k=2$), and
  Theorem MI / Lemma Δ / Theorem K / Props S, F of `c2k-zones.md` are all
  independent of C2k and remain valid. **Status: PROVED.**
- **Not claimed:** nothing about GTZ is asserted or refuted here; C2k was the
  weighted‑average pair form only (Scope Guard C).

Cites: VAL-THM-C2K-001 (`proofs/c2k-zones.md` — defines C2k, the Master
Identity, and Zones S, F); VAL-THM-C2-001 (`proofs/theorem-c2.md` — the
separate $k=2$ theorem).
