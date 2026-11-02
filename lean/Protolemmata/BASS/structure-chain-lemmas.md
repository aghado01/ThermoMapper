# BASS — the structure chain

Grounds what the mixing results rest on, and the checks a new reversible-jump kernel owes before it lands (F-X4). Gap G1 of the [lineage-gaps register](../../../issues/regression/bass/lineage-gaps-register.md), and the KHAOS pieces of E4 and E5.

**Status (2026-11-01).** SC-1, SC-2, SC-3 and SC-5 proved in `lean/Enthymemata/BassStructureChain.lean`, with no apologies. SC-4's walks are done, below, and they correct what SC-4 first expected. [`lineage-gap-checks.py`](lineage-gap-checks.py) keeps the numeric checks, and SC-0's measurement is [reported](../../../issues/regression/bass/spectral-gap-brief.md) apart: γ falls exponentially in N on every replicate cell, with all the growth at modes at the count cap, and one mode needing k = 2 nearly stops deepening at two free slots, as SC-5 gives. Weak directions, where a term's gain is comparable to ½ ln N, are unread.

**Statements that changed.**
- `gap_le_cut` drops `∑ π = 1`. The inequality never uses it: normalization is what makes the left side a variance, and `varianceOf_indicator` holds without it.
- SC-2 also holds for a stationary kernel, reversible or not (`gap_le_cut_of_stationary`).
- SC-3's claim that the precision is "neither" `GForm`'s nor `RidgeForm`'s unless g ∝ 1 needs a Gram with off-diagonal entries. On a diagonal Gram it is diagonal for any g (`khaos_precision_of_diagonal`), which is where orthogonal columns land. `khaos_precision_of_const` is the g ∝ 1 case.
- `khaos_precision_eq_add_inv` needs A invertible, and `khaos_precision_posDef` needs A positive definite and g₀ > 0. There, g₀ stands for the paper's g₀².
- `not_isMhPair_of_mul_ne_one` is added: SC-1's contrapositive as a statement.

---

## SC-0 — at fixed τ, the structure chain is a finite reversible chain

**Claim.** With knots at distinct data values and a capped count and order, F is finite for given data. With the coefficients and σ² collapsed and τ held, each structure move is a Metropolis–Hastings kernel on F, reversible with respect to π(S | τ, y). So Aldous & Fill's machinery applies directly, which it does not to BARS's continuous knots. A sampled τ adds one continuous coordinate; a sweep's block scan composes kernels, and a composition need not be reversible even when its parts are.

**Already in the code.** `BassKernelEnumeration` (tests) builds K(S, S′) and π(S | τ, y) at a fixed τ, from the moves' statements, on three fixtures, the continuous one with 2,263 states. It checks detailed balance and πK = π. Nothing computes a gap.

**Triage:** validate empirically. Compute 1 − λ₂ of the symmetrized kernel per move mix, and SC-2's bound over the basins. The fixtures are small and F grows fast with N, so these readings show scaling, not Friedman's barrier itself.

---

## SC-1 — reciprocity, and what KHAOS's printed ratios do

**Claim.** On any set with positive π and positive q both ways, r(s, t) = π(t)q(t, s) / (π(s)q(s, t)) satisfies r(s, t)·r(t, s) = 1, and the kernel q(s, t)·min(1, r(s, t)) satisfies detailed balance. So a printed pair whose product is not 1 is not a Metropolis–Hastings pair for any π and q. Reciprocity is cheap to assert, and belongs in F-X4's contract as a property test: for every enumerated move and its reverse, the product of the two log ratios is 0.

**Reading of arXiv 2510.25036v1 (manuscript only), in the paper's symbols, since this checks its displays.** The acceptance is log α = log-likelihood ratio + log prior ratio + log A_X, so A_X must be reverse over forward. KHAOS's birth draws q₀ with weight q₀^(−s_q)/c_q, coins with product Π_j η_j(q₀)^{χ_j}(1−η_j(q₀))^{1−χ_j} =: H(q₀), a degree with weight d^(−s_d)/c_d, and one of binom(d−1, q−1) partitions. With H̄ = (1/c_q)Σ_{q₀} q₀^(−s_q) H(q₀), the Metropolis–Hastings forms are:

| | Metropolis–Hastings | Printed |
|---|---|---|
| A_Birth(M) | P_D·binom·d^{s_d}c_d / (P_B·(M+1)·H̄) | P_D·binom·H(q₀) / (P_B·(M+1)·q₀^{s_q}c_q·d^{s_d}c_d) |
| A_Death(M) | P_B·M·H̄ / (P_D·binom·d^{s_d}c_d) | P_B·d^{s_d}c_d·H̄ / (P_D·M·binom) |
| A_Mutate1 | d_cand^{s_d}·binom_cand / (d_curr^{s_d}·binom_curr) | d_curr^{s_d}·binom_curr / (d_cand^{s_d}·binom_cand) |

So, as printed:
- the birth's coin product and q₀ weight are inverted, and q₀ is conditioned on where the death marginalizes it;
- the degree weight is inverted in birth, death and Mutate1, and the partition count in Mutate1;
- the death's M is inverted.

The printed pair's product is H(q₀)·H̄ / ((M+1)²·q₀^{s_q}c_q), not 1 (the check matches this form to 16 digits). The (M+1) double count with the prior ratio is already recorded in the port hub's register.

**Also in the manuscript.** P₁(x) = 1 in the recurrence (x is meant). q = Σ_{j=1}^M χ_j (j ≤ p is meant). The partition probability binom(d−1, d−1)⁻¹ (d−1, q−1 is meant). Degrees sampled up to p, where the prior's set stops at d_max. The g₀² conditional has no data term (the extensions brief notes it).

**What this is not.** Not a claim about the `khaos` R package, whose source has not been read, nor about the supplement's Appendix C (the η construction and delayed rejection), which was not read either. η reads the current model's activity, so the reverse birth must read it from the state after the death, which the printed form leaves implicit. Read the code before calling any of this more than a manuscript reading. By the no-parity rule, a KHAOS port derives its ratios and tests reciprocity, whatever the print says.

```lean
import Mathlib.Data.Real.Basic

noncomputable def mhRatio {S : Type*} (π : S → ℝ) (q : S → S → ℝ) (s t : S) : ℝ :=
  π t * q t s / (π s * q s t)

theorem mhRatio_mul_swap {S : Type*} (π : S → ℝ) (q : S → S → ℝ) (s t : S)
    (hs : π s ≠ 0) (ht : π t ≠ 0) (hst : q s t ≠ 0) (hts : q t s ≠ 0) :
    mhRatio π q s t * mhRatio π q t s = 1 := by
  sorry

theorem mh_detailed_balance {S : Type*} (π : S → ℝ) (q : S → S → ℝ) (s t : S)
    (hs : 0 < π s) (ht : 0 < π t) (hst : 0 < q s t) (hts : 0 < q t s) :
    π s * (q s t * min 1 (mhRatio π q s t)) = π t * (q t s * min 1 (mhRatio π q t s)) := by
  sorry
```

---

## SC-2 — a gap is bounded by any set's flow

**Claim.** For a kernel P reversible with respect to π, and any γ with γ·Var f ≤ E(f, f) for all f, every set A has

γ·π(A)(1 − π(A)) ≤ flow(A) = Σ_{s∈A, t∉A} π(s)P(s, t).

**Proof.** Take f = 1_A. Var f = π(A)(1 − π(A)), and reversibility gives E(1_A, 1_A) = ½Σ π(s)P(s, t)(1_A(s) − 1_A(t))² = flow(A).

**Why it matters.** It is the easy half of Cheeger, in the same Dirichlet-form language as [SG-1](sensitivity-generator-lemmas.md). So relaxation time is at least π(A)(1 − π(A)) / flow(A), and any basin gives a lower bound computable on an enumerated kernel. The check has a true gap of 0.295, with flow/variance above it on 20,000 random sets.

```lean
import Mathlib.Algebra.BigOperators.Field
import Mathlib.Data.Fintype.BigOperators
import Mathlib.Data.Real.Basic

open Finset

variable {S : Type*} [Fintype S] [DecidableEq S]

noncomputable def dirichlet (π : S → ℝ) (P : S → S → ℝ) (f : S → ℝ) : ℝ :=
  (1 / 2) * ∑ s, ∑ t, π s * P s t * (f s - f t) ^ 2

noncomputable def varianceOf (π : S → ℝ) (f : S → ℝ) : ℝ :=
  ∑ s, π s * f s ^ 2 - (∑ s, π s * f s) ^ 2

noncomputable def cutFlow (π : S → ℝ) (P : S → S → ℝ) (A : Finset S) : ℝ :=
  ∑ s ∈ A, ∑ t ∈ Aᶜ, π s * P s t

theorem gap_le_cut (π : S → ℝ) (P : S → S → ℝ) (hπ : ∑ s, π s = 1)
    (hrev : ∀ s t, π s * P s t = π t * P t s) (γ : ℝ)
    (hgap : ∀ f : S → ℝ, γ * varianceOf π f ≤ dirichlet π P f) (A : Finset S) :
    γ * ((∑ s ∈ A, π s) * (1 - ∑ s ∈ A, π s)) ≤ cutFlow π P A := by
  sorry
```

**Fishable:** no Cheeger, spectral-gap or conductance result turned up in Mathlib (`rg -il` on those words finds nothing), so this is self-contained. The proof is two finite-sum computations.

---

## SC-3 — KHAOS's posterior precision

**Claim, in the paper's symbols.** With S₀ = g₀²D(g)(ΨᵀΨ)⁻¹D(g) and G_ml = (g₀²g_m g_l + 1)/(g₀²g_m g_l), G∘ΨᵀΨ = ΨᵀΨ + S₀⁻¹ (checked to 1e-13). So the paper's Σ_n = (G∘ΨᵀΨ)⁻¹ is the ordinary conjugate form, and positive definite as a sum of two positive definite matrices.

**Why it matters.** Unless g ∝ 1, this precision is neither `GForm`'s scalar times the Gram nor `RidgeForm`'s Gram plus a diagonal on the Gram's spectrum. It is a third `CoefficientForm`, with its own fast path to find.

```lean
import Mathlib.Data.Real.Basic
import Mathlib.LinearAlgebra.Matrix.Defs

theorem khaos_precision_entry {m : Type*} (A : Matrix m m ℝ) (g : m → ℝ) (g₀ : ℝ) (i j : m)
    (hg₀ : g₀ ≠ 0) (hi : g i ≠ 0) (hj : g j ≠ 0) :
    (g₀ * g i * g j + 1) / (g₀ * g i * g j) * A i j = A i j + A i j / (g₀ * g i * g j) := by
  sorry
```

---

## SC-4 — the proof-walks

Both walks read the papers on 2026-11-01. Each verdict names the lines it rests on.

### Ronen, Saarinen, Tan, Duncan & Yu (arXiv 2210.09352v1), onto BASS at fixed τ

**The result, as printed** (`03_lb.tex`, Theorem 1). For a simplified Bayesian CART, one tree with grow and prune only, discrete covariates, y bounded and σ² a fixed hyperparameter: with probability at least 1 − 1/n, t_mix ≥ exp(n·Var f₀/(2aσ²) − O(√(n log n))).

| Step | Their proof | BASS at fixed τ | Verdict |
|---|---|---|---|
| The chain | Finite tree space, Metropolis–Hastings on the collapsed π(T \| y) | SC-0: finite, reversible, collapsed | **Transports.** BASS's F grows with N, since knots sit at data values; their space is fixed in n, because covariates are discrete. So their constant C(T, T′, α, β, m, d) has no analogue, and BASS's proposal probabilities carry polynomial factors of N |
| Mixing from conductance | t_mix ≥ 1/(4Φ\*), Levin–Peres–Wilmer Theorem 7.3 | Same, or SC-2 in gap form | **Transports** |
| **The bottleneck set** | S(v, τ), the trees whose root splits on (v, τ). Every move adds or removes a single terminal split, so S and Sᶜ meet at **one edge**, between the trivial tree T₀ and the depth-1 tree (`06_proofs.tex`, Lemma 3) | A structure is an unordered collection of terms ([D01](../../../docs/Maths/BayesianRegression.md)). Any term can be removed at any time, and S(t), the structures holding t, can be left from every one of its states | **Fails. This is the step that uses trees.** There is no root, so no single gate |
| The flow | π(T₀)Q(T″ \| T₀) ≤ π(T₀)/(2md): one state's mass | `gap_le_separator` is the form it takes when a separator W replaces the single state | **Posed** as SC-5 |
| The likelihood ratio | Lemma 4: log L(T₀)/L(T) concentrates at −nΔ/(2aσ²), Δ the impurity decrease, for one fixed T | One fixed pair transports. A separator W needs it uniformly over W, and the term measure helps: Σ_W p(s) ≤ 1 pays for W's size at a cost of log(1/p(S_A)) ≈ M·K·log(pN) | **Transports only if W exists** |
| Two good trees | T, T′ with f₀ piecewise constant on both, differing at the root. Lemma 3 takes the smaller of two likelihood ratios because it cannot say which side has π ≤ ½, so both trees need a large Δ | No exact pair in general: f₀ is approximated, and F-M4's modes sit tens of nats apart | **Changes.** Equal posteriors are not needed, only a large explained variance on both sides of the cut |

**The lemma the walk poses.** The first failure is not a technicality. Because a BASS state can hold both a term and its substitute, it has a path a tree does not.

**SC-5 (the superset path).** Adding k columns cannot lower R², so `collapsedLogMarginal_superset_ge`: the collapsed log marginal falls by at most k·½ ln(1 + g), whatever N and whatever variance the columns explain. Removing the old terms afterwards is the same step read from the other side. So S → S ∪ T → T has no likelihood valley of the size this section's first draft assumed. The depth along it is the Occam term, which is constant in N at fixed τ and (k/2)·ln N at g ≈ N. On top sit the count factor and the proposals' entropy, about K·ln(pN) per exact term (fewer where nearby knots serve), since a birth has to draw it.

**Evidence already in the record.** F-M2 measured an admitted birth's Δℓ at −3.1, −2.6 and −3.8 at N = 500, 300 and 1,200, near −½ ln N and about −5 with the prior and proposal parts, while a death's median on Friedman is −234 nats ([measurements](../../../issues/regression/bass/local-move-measurements-brief.md)). So births are how a mode is left cheaply, and deaths are not. F-M4's sentence that single births and deaths each pay tens of nats overstates births against F-M2's own reading; the deaths' tens to hundreds stand.

**What the barrier's growth would have to be.** An exponential bound by this method needs a separator W with ln(π(A)/π(W)) ≥ cN: every structure in W must have lost a fixed per-case explained variance relative to A, uniformly. The superset path rules that out unless something closes it. So:
- *Expected at fixed structures:* polynomial in N, of degree growing with the number k of terms exchanged at once: between N^{k/2} and about N^{k(K + ½)}. Polynomial does not mean fast. At N = 500 the births alone cost about e^{−5} each, so k = 4 costs e^{−20} before any proposal entropy. At fixed τ the N-dependence drops to the proposals. The [path profile](../../../issues/regression/bass/path-profile-brief.md#reports) adds that a term's measure stays in π for as long as the term is held unpaid, so it accumulates along a route: 12.9 nats with one factor and 21.3 with two at N = 500, Occam included. Four unpaid two-factor terms are up to 85 nats by the superset route, and the best routes found on Friedman avoid them by passing a likelihood valley instead.
- *What could still give a superpolynomial bound:* k or M growing with N, since M_n grows for a truth outside the span and the barrier is then M_n·ln N. The count cap M_max and the numerical gate D06 can close a superset. The default cap is min(1000, N − 1), which rarely binds. D06 refuses near-exact collinearity, where a removal usually costs little, with one exception the path profile found. A hinge at a variable's lowest knot and its mirror at the highest are one linear column with the intercept, so D06 refuses a state holding both, and removing the one held is load-bearing: −991 nats on x4 at N = 1000. Production crosses by one redraw, so the sampler is not stuck there, but the superset path is closed.
- *This corrects the first draft,* which expected linear growth from the load-bearing deaths and so the same exponential form. That reasoning looked at deaths only.

**The test.** SC-0's N sweep with SC-2's bound: the flow over a basin should fall polynomially in N. A fall exponential in N would falsify this reading, and would point to a separator worth looking for. The [spectral-gap brief](../../../issues/regression/bass/spectral-gap-brief.md), written before this walk, reads both: linear as the register's first expectation, and logarithmic as a rival, "a birth goes first", with c near ½ per spare term. This walk gives a reason to expect the rival at τ held at unit information. The brief's replicate ladder holds F fixed, so it reads the Occam part alone; its level ladder is where the proposals' entropy enters.

**Kim & Roćková (arXiv 2306.00126v1): what they add.**
- *Upper bound.* With grow and prune only and a signal connected along a tree (their Assumption (a)), τ_ε is polynomial, cubic in n at L = L_max, by canonical paths (`main.tex`, Theorem `eq:theo1`). Its proof needs tree selection consistency first.
- *Lower bound.* For an isolated deep signal, τ_ε ≳ ((n^{c−3/2}/4 − 1)/C)^{L−2} with split prior n^{−c}, c > 5/2 (Theorem `theo:poor`). This is the same shape as SC-5: a per-node Occam factor raised to the number of nodes grown through noise before the likelihood pays. The paper reads it as exponential in n at L ~ log n; as printed it is n^{Θ(log n)}, so superpolynomial. Its bottleneck is also hierarchical: a node is reachable only from trees holding its parent.
- *A fix.* Twiggy proposals attach or detach whole twigs, and mix polynomially without connectivity. The analogue here is a move exchanging several terms at once, recorded as a finding and not designed here.
- *The nearest upper bound for BASS is theirs by comparison only.* They restate Yang, Wainwright & Jordan (2016): rapid mixing of one-site Metropolis–Hastings for a g-prior linear model, by canonical paths after stepwise selection. That is BASS's shape, a set on a dictionary under a g-prior, and it rests on asymptotic unimodality, which adjacent knots' near-duplicate columns and F-M2's modes do not give. That paper is not on disk. On Friedman at N = 500 (one run, one data seed), near-duplicates are not what sets the depth: the [path profile](../../../issues/regression/bass/path-profile-brief.md#reports)'s lattices refuse no state, and swaps between near twins add steps and k without lowering the route's minimum.

### Woodard, Schmidler & Huber (EJP 2009)

**Persistence, exactly** (eq. 4). For A with 0 < π[A] < 1 and a density φ, γ(A, φ) = min{1, φ[A]/π[A]}. The persistence of A at β is γ(A, π_β), with π_β ∝ π^β. A has low persistence when tempering shrinks it: a tall, narrow peak beside lower, wider ones of comparable mass (Section 5).

**The bound** (Corollary 3.1, for any choice of temperatures). Gap(P_pt) ≤ 12·sup over β ∈ [β\*, 1], β′ ∈ [0, β\*) of γ(A, β)·max{Φ_{T_β}(A), δ(A, β, β′), δ(Aᶜ, β, β′)}, and Gap(P_st) ≤ 192·[the same without the Aᶜ overlap]^{1/4}, where Φ_{T_β} is A's conductance at β and δ the overlap of the two tempered laws on A. Corollary 3.2: tempering is torpid if A's conductance on [β\*, 1] and its persistence on [β′, β\*) both fall exponentially in the problem size, and β′ is the lowest temperature or the overlaps from [β\*, 1] to [0, β′) also fall. Theorem 3.1 takes any densities φ_k, tempered or not.

**Why it matters to the record.** The paper's conclusion is that the failure can pass unnoticed: swaps accept, the hot chain mixes, and R̂ works poorly when modes have small basins (Section 6). F-R1's fits fail R̂ outright, so that part is a different symptom. The reference arms' missing a1, whose basins hold at most 1.3%, is the part that resembles it.

**Measuring it on BASS from a flat-histogram run.**
1. BASS's power posterior is on the likelihood path, π_b ∝ prior·L^b, not π^β, so take Theorem 3.1's densities φ_k. γ(A, φ_k) = min{1, π_b[A]/π_1[A]} has the same form.
2. A flat-histogram run gives ln g(ℓ) over the prior. For a set A of structures it does not give ln g_A(ℓ). What it adds is the occupancy: given ℓ, the multicanonical law is the prior restricted to level ℓ, so P̂(A | ℓ), read per level from draws labelled by basin, gives g_A(ℓ) = g(ℓ)·P̂(A | ℓ).
3. Then π_b[A] = ∫g_A(ℓ)e^{bℓ}dℓ / ∫g(ℓ)e^{bℓ}dℓ, for every b from one run. That needs a basin label on structures and a per-level tally, and neither is built. The level's conditional has to be sampled well for P̂(A | ℓ) to be right, so the reading inherits the mixing it measures.
4. Φ_T(A) at b ≥ b\* is a separate reading: flow over mass, which SC-2 prices on an enumerated kernel and a chain's own A-crossings estimate.
5. "Exponentially decreasing" asks for a size sweep. One fixture gives γ's value at one N, which reads as a bound through Corollary 3.1 and not as torpidity.

**The Potts toy's q = 3 question: resolved, the sizes are too short.**
- *What the theorem says, corrected.* Bhatnagar & Randall proved q = 3 at α = α_c, and Woodard et al. note their argument may extend (Section 4.2). Woodard et al. state q ≥ 3 and α ≥ α_c, with the single-site proposal, for any number and choice of temperatures. Their proofs of Propositions 4.2–4.4 are in Woodard's thesis.
- *The settings match.* `MeanFieldTunnellingReading` runs a coupling J = J_c/0.6, so the b = 1 target has α ≥ α_c. It uses the lumped occupation chain, whose Hastings ratio the tests pin to the single-site configuration chain's, and a self-sized ladder, which the theorem's "any choice" covers. It reads replica round trips in ℓ₁, which is colour-blind where the theorem's A = {σ₁ > M/2} is not, and a round-trip rate rather than a gap. Those differences make it a proxy; they do not explain the fit.
- *The sizes.* The barrier per site at α_c, from f(a) along the symmetric path, is 0.00113 nats at q = 3 and 0.0065 at q = 4 (`lineage-gap-checks.py`), which reproduce the test's own header (0.0011, 0.0065). So at q = 3, e^{cN} moves by only about e^{2.2} over the whole sweep from N = 64 to 2,048, where a fitted N³ moves by e^{10.4}. The test's header says as much: the forms may not separate there. An exponential at this rate overtakes N³ near N ≈ 27,000 at q = 3, and near N ≈ 3,800 at q = 4, which is why q = 4 stops crossing within the sweep.
- *So* N^(−3.0) at q = 3 is the pre-asymptotic diffusion, not a counterexample. Reading q = 3 would take N past 10⁵. The line in [Samplers](../../../docs/Maths/Samplers.md) that the reading agrees with the theorem holds at q = 4; at q = 3 the reading cannot speak to it. That line is outside this chip, so it is reported rather than changed.
- *What this implies for BASS.* The same arithmetic applies to any barrier measured at one size: a small per-case rate reads as polynomial until N passes about (ln N)/rate.

---

## SC-5 — the superset path and the separator

Posed by the Ronen walk above, and proved in Lean.

- `collapsedLogMarginal_superset_ge`: for any profile h increasing in R², L(M + k, R²′) ≥ L(M, R²) − k·c, with c the Occam cost per column. `jeffreys_superset_ge` (in `BassContractionCover`, for its imports) discharges the hypothesis on D03's Jeffreys form: `jeffreysProfile_mono` for the profile, and `rSquaredOf_mono` for the nesting, that R² cannot fall when the column space grows. The proper IG(c, c·s²) form has the same shape, with c added to the exponent and c·s² to the sum, so it should follow the same way, though it is not written out. This bounds the log marginal, not π. The [path profile](../../../issues/regression/bass/path-profile-brief.md#reports) read it on 32,577,747 superset edges with no violation, and the least slack is 0, where a term adds nothing.
- `cutFlow_le_separator`, `gap_le_separator`: if every step out of A lands in W or stays in A, then γ·π(A)(1 − π(A)) ≤ e·π(W), with e the largest chance of stepping from W into A. Ronen et al.'s Lemma 3 is W = {T₀}.
- `gap_le_cut_of_stationary`: the bound for a stationary kernel that is not reversible, the case SC-0 raises for a sweep's block scan. Its γ is the constant of the kernel's Dirichlet form, which is that of the additive symmetrization.
