# GTZ(6,3) Closure Attempt — Agentic Research Brief

**Target:** Resolve (prove or refute) the first open case of the Goreinov–Tyrtyshnikov–
Zamarashkin (GTZ) row-selection hypothesis: every `A ∈ St(6,3)` (real 6×3, `AᵀA = I₃`)
has a 3×3 row submatrix `A_I` with `σ_min(A_I) ≥ 1/√6`.

**Inputs to ingest first, in this order:** `main.pdf` (self-contained draft template),
then `gtz.pdf` (campaign report — read the abstract and §12 "Open questions" before
anything else), then the `proofs/*.md` files: `duality.md`, `case-a.md`,
`sp-verification.md`, `theorem-c2.md`, `pair-lemma.md`, `extension-lemma.md`,
`c2k-zones.md`, `c2k-middle.md`, `slice-framework.md`, `boundary-obstruction.md`.

---

## 0. Mission and epistemic contract

This is a genuine open problem (verified against the literature: Sengupta–Pautov,
arXiv:2604.05944, proves the k=2 case; the n−2 case follows by duality; (6,3) is
confirmed as the smallest open case). A prior AI-assisted campaign (documented in
the ingested files) produced real, checkable reductions and one genuine no-go
result, and localized the remaining difficulty precisely. Do not re-derive from
scratch — build on this. But do not trust it uncritically either: several of its
claims were spot-checked by hand and held up, but the corpus is large, some of its
numerical claims (e.g. an 80-digit St(259,5) counterexample to a side conjecture)
were not independently re-verified before being handed to you, and AI-generated
math of this style is exactly where subtle, confidently-stated errors hide.

**Non-negotiable rules for the whole engagement:**

1. Every claim you produce gets one of exactly three tags: **PROVED** (complete
   argument, machine-checked in exact or certified arithmetic, no numerical step
   load-bearing), **NUMERICALLY SUPPORTED** (true in every test run but not proved),
   or **OPEN**. Never upgrade a tag without new justification. Never let a PROVED
   claim depend on a NUMERICALLY SUPPORTED one.
2. Before using any lemma from the ingested `.md` files as a building block for new
   work, re-derive or re-verify it yourself (symbolically, exactly) rather than
   citing its PROVED tag on faith. If you find an error in something tagged PROVED
   in the source material, stop, flag it prominently, and do not silently patch
   around it — this is more important than making progress.
3. If at any point a numerical search finds `f(A) < 1/√6 − 10⁻⁶` (i.e. a candidate
   counterexample to GTZ itself), **stop all other work immediately**, do not
   average it away or dismiss it as noise, reproduce it in at least two independent
   exact-arithmetic implementations, and escalate to Pavel before doing anything
   else. Hundreds of thousands of prior samples never found one; if you do, treat
   it with maximal suspicion of a bug first, but do not bury it either way.
4. A full proof of GTZ(6,3) requires a certificate you can hand to a skeptical
   reader for independent verification — exact rational/algebraic arithmetic, or a
   rigorously certified numerical method (verified interval arithmetic with proven
   error bounds), not floating-point SDP output alone. Floating-point SDP results
   are a fine first pass but are NUMERICALLY SUPPORTED, not PROVED, until the
   certificate is extracted and checked exactly.

---

## 1. Self-contained problem statement

For `A ∈ St(6,3)`, write rows `r₁,…,r₆ ∈ ℝ³`, `P = AAᵀ` (rank-3 orthoprojector,
`P² = P`, `tr P = 3`, `0 ⪯ P ⪯ I`), `ℓᵢ = P_ii = ‖rᵢ‖²`, `c_ij = P_ij = ⟨rᵢ,rⱼ⟩`. For a
triple `T ⊂ {1,…,6}`, `|T|=3`, `P_TT` is the eigenvalues-of-`A_TA_Tᵀ` principal
submatrix; `T` is **good** iff `λ_min(P_TT) ≥ 1/6`. Define
`F(A) = max_{|T|=3} λ_min(P_TT) = f(A)²`. GTZ(6,3) is the claim `min_{St(6,3)} F ≥ 1/6`.

**Key structural fact (useful reformulation):** the map `A ↦ P` identifies `St(6,3)`
modulo the right-`O(3)`-action with the real Grassmannian `Gr(3,6)`, realized as
`{P ∈ Sym(6,ℝ) : P² = P, tr P = 3}` — a smooth compact real algebraic variety of
dimension `k(n−k) = 9` inside the 21-dimensional space of symmetric 6×6 matrices.
Working directly with `P` (21 ambient variables, degree-2 defining equations) avoids
the `O(3)` gauge redundancy of working with `A` directly, and is the natural
coordinate system for the SOS/Positivstellensatz work in Track A below.

---

## 2. What you may treat as a trusted starting point (after your own re-check)

- **Duality** (`duality.md`): `GTZ(n,k) ⇔ GTZ(n,n−k)` via CS-decomposition of a
  completed orthogonal matrix. Elementary, low risk.
- **Case A reduction** (`case-a.md`): a row with `ℓᵢ ≤ 1/6` can be removed and the
  problem reduces to `St(5,3)` (itself resolved via duality to the proved `k=2`
  case). Reduces the search to the **core region** `ℓᵢ > 1/6 ∀i`.
- **k=2 case** (Sengupta–Pautov, independently verified in `sp-verification.md`):
  trust this; it is peer-reviewed, published literature, independently confirmed
  in `2604_05944v5.pdf`.
- **Extension Lemma / Schur criterion / trace branch** (`extension-lemma.md`,
  and §3 of `main.pdf`): the moment identity `Σ_{p∉{i,j}} q(p) = 1/3 − h(µ₁) − h(µ₂)`
  and the resulting `h`-branch sufficient condition. This was hand-verified
  independently (symbolically and numerically) before this brief was written and
  is solid.
- **Single-pair boundary obstruction** (`boundary-obstruction.md`, Theorem BND):
  a genuine no-go — no lemma of the form "small-excess qualifying pair ⇒ good
  triple through it" can exist at any threshold. **Do not attempt to patch or
  extend the single-pair/second-moment program — it is provably dead.** Any new
  attack must use either the full active-set structure or a genuinely different
  mechanism (this is why Tracks A–C below do not build on pairs).
- **Slice framework** (`slice-framework.md`): on the equal-leverage locus
  `ℓᵢ ≡ 1/2` (⇔ `J = 2P − I` a zero-diagonal symmetric involution), goodness of
  `T={i,j,k}` is `τ_T := a_{ij}a_{ik}a_{jk} − (a_{ij}²+a_{ik}²+a_{jk}²)/3 + 2/9 ≥ 2/27`.
  Pinned block spectra, `K₄` coherence balance, and two further obstructions
  (min-pair-`K₄` selection dead; unsigned `τ`-averaging cannot certify) are proved
  there. The pentagonal matching configuration (PMC) is an exact slice point with
  `G_PMC ≈ 0.2061074` (root of an explicit cubic), conjectured (not proved) to be
  the slice minimizer, margin `≈ 0.0394` over `1/6`.
- **C2k side-conjecture is false** (`c2k-middle.md`) — this was a candidate route
  toward a general-`k` pair-averaging bound; it's dead for `k ≥ 4` and irrelevant to
  (6,3) specifically (which goes through the `k=3`-specific Extension Lemma, not
  C2k). No need to revisit unless you find independent use for it.

---

## 3. Confirmed dead ends — do not retry

- Any single-pair / second-moment argument for the diffuse regime (Theorem BND).
- Unsigned averaging of `τ_T` over all 20 triples, globally or `K₄`-locally, on the
  slice (Prop. 12 in `slice-framework.md`) — the mean is `1/45 < 2/27`.
- Low-degree (≤6) generic Lasserre/SOS relaxation on the slice without symmetry
  reduction and without the equality-manifold vanishing condition imposed: level 2
  reproduces exactly the (too-weak) averaging bound `−7/135`; level 3 (969×969)
  only closes 8.6% of the gap.

---

## 4. Track A (primary target): Positivstellensatz certificate on `Gr(3,6)`

This is the real prize: a certificate valid on the *entire* manifold, not just the
equal-leverage slice, closing part (b) of the KKT reformulation in
`boundary-obstruction.md` §5 (no KKT point of `F` below `1/6`) — but phrased to
avoid nonsmooth-KKT machinery entirely by working with an infeasibility
certificate instead.

**Formal target.** Show infeasibility, on the compact variety
`V = {P ∈ Sym(6,ℝ) : P² − P = 0, tr P − 3 = 0}`, of the system "for every one of
the 20 triples `T`, `P_TT − (1/6)I₃` is not PSD."

**Correctness pitfall to avoid:** encoding "`M ⪰ 0`" via *leading* principal minors
`≥ 0` (Sylvester) is valid for strict positive-definiteness but is **not**
sufficient for PSD at the boundary (`diag(0,−1)` has both leading minors `≥ 0` but
isn't PSD). Use the correct fact — symmetric `M ⪰ 0 ⟺ all principal minors ≥ 0`
(not just leading ones) — or work directly with the characteristic polynomial of
`P_TT` and a Sturm-sequence / discriminant argument for "smallest root `< 1/6`."
Get this right before building anything on top of it; it is the kind of error that
would invalidate an entire certificate silently.

**Symmetry reduction.** `S₆` acts on `V` by simultaneously permuting rows/columns
of `P`; the objective and constraint set are invariant. Use Gatermann–Parrilo-style
symmetry reduction (isotypic decomposition of `Sym(6,ℝ)` as an `S₆`-module) to
block-diagonalize the moment/SOS matrices before solving — this is the single
biggest lever for making level-4+ relaxations tractable at this dimension (9-dim
variety, 21 ambient variables, `|S₆| = 720`). Tools: TSSOS or SPARSEPOP for
sparsity/symmetry-exploiting moment-SOS; a symbolic pass in GAP or Sympy to build
the symmetry-adapted basis if off-the-shelf tooling doesn't handle this
representation directly.

**Vanishing constraint.** `boundary-obstruction.md` Prop. 11.3 proves any valid
certificate must vanish on the equality manifold (the known extremal configurations
— catalog them from `main.pdf` and the campaign's numerics first). Impose this
explicitly (e.g. via the ideal generated by the known extremal points, or by
subtracting off a term that vanishes there) rather than hoping the solver finds it —
this is likely why unconstrained low-level relaxations plateau at the averaging
bound.

**Staging:** start at level 4 on the symmetry-reduced problem; if infeasible within
budget, characterize exactly where the SDP bound saturates (is it still hitting the
averaging value, or has it moved?) before escalating to level 5+, since that
diagnostic tells you whether more degree or a different multiplier structure is
needed.

---

## 5. Track B (fallback / faster win): close the slice theorem alone

Lower-dimensional (4-dim moduli space after quotienting), self-contained, and
worth pursuing **even as a standalone result independent of whether Track A
succeeds** — the PMC-displaces-icosahedron finding is already a nice, checkable,
apparently novel observation.

- Same symmetry-reduction approach as Track A, but now use the *full* symmetry of
  the slice: `S₆` row-permutations combined with the sign-switching group
  `(ℤ/2)⁶` (`a_ij ↦ ε_iε_j a_ij`, a symmetry of `J² = I` and the zero-diagonal
  condition — note the global flip `ε ≡ −1` acts trivially, so the effective group
  has order `2⁵·720 = 23040`). This connects directly to Seidel-switching /
  two-graph theory; worth checking whether existing two-graph classification tools
  give useful structural shortcuts (e.g. a finite list of switching classes on 6
  points) before setting up the SDP.
- Target: certify `τ_T ≥ 2/27` fails to hold simultaneously for all 20 triples,
  at Lasserre level ≥ 4 with the `K₄`-balance identity (`Σ_{T⊂F} q_T = 0`) imposed
  as an equality constraint, and with vanishing forced at the PMC point (the
  presumed near-tight extremal — verify tightness numerically first with a finer
  search than the campaign's 64-start run before committing compute to a
  degree-6+ certificate).
- If PROVED: this is a clean, publishable, self-contained lemma regardless of
  Track A's outcome, but flag clearly in any writeup that it does **not** by
  itself imply GTZ(6,3) — the reduction from general `A` to the slice is not
  established (Case A only reaches the core region `ℓᵢ > 1/6`, not `ℓᵢ ≡ 1/2`).

---

## 6. Track C (orthogonal attack): certified global optimization, no SOS required

Given large compute is available, run this in parallel with A/B rather than as a
last resort — it's a genuinely different route to the same conclusion and doesn't
depend on finding an elegant algebraic certificate.

- Rigorous branch-and-bound over `St(6,3)` (or the `Gr(3,6)` projector
  representation) using verified interval arithmetic (e.g. Arb/MPFI, or
  interval-Newton with a rigorous C¹ bound on `F`), adaptively refined, to
  directly certify `F(A) ≥ 1/6` everywhere or isolate a genuine counterexample.
- Handle non-smoothness of `F = max_T λ_min(P_TT)` by branching on the active
  triple within each box (finitely many `T`, so this is tractable) and using
  interval eigenvalue bounds for `3×3` symmetric matrices (closed-form via the
  trigonometric solution of the cubic, which interval-arithmetizes cleanly) rather
  than a generic nonsmooth solver.
- Use the campaign's ~9 known exact extremal configurations as required refinement
  centers — the box cover must resolve arbitrarily close to each one and show the
  local bound doesn't dip below `1/6 − ε` for shrinking `ε`, which is a strong
  practical check of second-order optimality at each even before global coverage
  completes.
- This can succeed even where no low-degree SOS certificate exists (SOS
  Positivstellensatz certificates can have arbitrarily high minimal degree even
  for true statements), so treat it as a genuinely independent shot, not a
  consolation prize.

---

## 7. Tooling requirements

- Exact/symbolic: Sympy or equivalent, for re-deriving ingested lemmas and for
  final certificate verification. PSLQ/LLL (e.g. via PARI/GP or Sympy's
  `nsimplify`/`mpmath.pslq`) to rationalize floating SDP output into exact
  algebraic form — expect the relevant number field to involve `√5` (icosahedral
  points) and the PMC's defining cubic, not just `ℚ`.
- SDP: a high-precision solver (SDPA-GMP or similar) for final certificate
  extraction once a floating-point solve (SCS/MOSEK/mainstream solver) has found a
  candidate — floating precision is fine for search, not for the final PROVED tag.
- Symmetry-exploiting SOS: TSSOS, SPARSEPOP, or a custom Gatermann–Parrilo
  reduction; GAP or Sympy for constructing symmetry-adapted bases if needed.
- Certified interval arithmetic: Arb (via python-flint) or MPFI, for Track C.
- Deterministic seeding and full reproducibility for every numerical run, matching
  the house style already used in the ingested `.md` files (script name, seed,
  exact pass/fail count reported).

---

## 8. Verification protocol

1. **Reproduce before extending.** Re-run (or re-derive symbolically) every
   numeric sanity check cited as already-passing in the ingested files before
   using its conclusion. Report any discrepancy immediately, however small.
2. **Adversarial second pass.** For any claim you intend to tag PROVED, run a
   second, independent verification attempt specifically trying to break it —
   different implementation, different simplification path, or (if feasible)
   route the verification through a separate reasoning process from the one that
   produced the claim. Do not let the same derivation grade its own homework.
3. **No claim stronger than its source.** If a result depends on a
   NUMERICALLY SUPPORTED input, it is at most NUMERICALLY SUPPORTED itself —
   propagate tags honestly through the dependency graph.
4. **Counterexample protocol.** See rule 3 in §0 — this is the highest-priority
   interrupt condition in the entire brief.

---

## 9. Deliverables

- An updated status table in the exact style of `gtz.pdf` §12 / `slice-framework.md`
  §10, with every new claim tagged and cross-referenced to its proof file.
- Any new lemmas as standalone `proofs/*.md` documents matching the existing house
  style (overclaim-guard header, dependency list, exact-arithmetic proof sketch,
  numeric sanity block with deterministic seeds).
- A final honest summary addressed to Pavel: what closed, what didn't, and — if
  nothing closed — the sharpest possible restatement of exactly where the
  remaining gap is, in the spirit of `boundary-obstruction.md` §0. A well-localized
  "still open, and here is precisely why" is a legitimate and useful outcome; do
  not let deliverable pressure turn into overclaiming.

---

## 10. Suggested compute staging

1. Cheap first: reproduce existing numerics (few CPU-hours), finer slice search
   near PMC to firm up Track B's target before committing SDP budget to it.
2. Track B (slice, 4-dim, symmetry-reduced) — likely the best compute-to-payoff
   ratio; attempt first, both as a real target and as a testbed for the
   symmetry-reduction tooling before scaling to Track A's larger ambient space.
3. Track C in parallel (embarrassingly parallelizable branch-and-bound) — start
   coarse, refine adaptively; good use of "large computational resources" that
   doesn't bottleneck on Track A/B's algebraic setup work.
4. Track A (full 9-dim Grassmannian, level ≥4 symmetry-reduced Lasserre) — the
   heaviest lift; sequence after B's tooling is validated.
