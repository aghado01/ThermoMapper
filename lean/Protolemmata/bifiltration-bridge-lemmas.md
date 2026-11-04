# The (β, θ) thermodynamic bifiltration — physics→topology bridge, refined

The seam between the Swendsen–Wang spin model and the persistent-homology engine, stated
so that the **deep probabilistic content is cited** (Grimmett 2006, *The Random-Cluster
Model*) and **everything else is a provable obligation**. This sharpens the earlier
"axiomatize all of `hmono`" framing from [fable-bifiltration-lemmas.md](fable-bifiltration-lemmas.md):
the temperature-monotonicity of spin agreement is *not* one opaque axiom — it factors into
two cited theorems and a chain of calculus/combinatorics we discharge ourselves.

Same methodology as [spc-pkwang-lemmas.md](../Archeion/PKWang-legacy-2026-08-09/spc-pkwang-lemmas.md)
(Lemma B took the closed-form bond probability as given and proved the algebra from there). Here we
quarantine exactly two inputs and prove the plumbing around them.

**Status (2026-11-04).** A1 is proved in [RandomCluster.lean](../Lemmas/RandomCluster.lean) with no
added axioms, along with P1–P4. The rest compiles in [Bifiltration.lean](../Enthymemata/Bifiltration.lean)
over A2. Those two files are the current form. The block below is the draft they came from, and it
does not compile as written. The compiled files differ in these ways:

- `IsRandomCluster` is defined rather than opaque. The measure is `rcProb`, `randomCluster` shows the
  hypothesis can be met, and `rcProb_one_edge` checks the definition against Grimmett's (3.3).
- `pe` is `noncomputable`, the cut graphs take the `Std.Symm` and `Std.Irrefl` fields, and P3 names
  `SimpleGraph.Reachable.mono`.
- The measure and A2 take `IsEdgeParameter` (symmetric, in `[0, 1)`), which is what Grimmett's
  theorems cover, so the theorems take symmetric `J` and `β₁ ≥ 0`.
- T1 has its threshold half, `cutGraph_mono`, and `cutGraph_reachable_mono`, the nesting of clusters.
- `connCutGraph` is the same family for connection probability, which is what SPC's co-membership
  columns estimate. It needs no axiom.
- T2 and T3 are not declared. T2's type is not a proposition and T3's statement is `True`, so
  both wait for their codomains.

## What is axiom vs. obligation

| # | Statement | Status | Source / strategy |
|---|-----------|--------|-------------------|
| A1 | `comparison_inequality` — coordinatewise `p₁ ≤ p₂` ⇒ `φ_{p₁,q} ≤_st φ_{p₂,q}` (q ≥ 1) | **proved** | Mathlib's `holley`, the lattice inequality per edge, and Grimmett's (3.11) by the modular law for subspaces. Grimmett **Thm 3.21** states it for one `p`, and **§11.4** says it stays valid per edge. `q ≥ 1` load-bearing: positive association fails in general for `q < 1` (see before the proof of Thm 3.8). |
| A2 | `edwards_sokal_twopoint` — `⟨δ_{sᵢsⱼ}⟩ = 1/q + (1−1/q)·φ(i↔j)` | **AXIOM** | Grimmett **Thm 1.16** (Edwards–Sokal coupling); with parameters per edge, **(11.25)** in §11.4. |
| P1 | `pe_mono_in_beta` — `1 − exp(−βJ)` ↑ in β (J ≥ 0) | **proved** | `Real.exp_le_exp`. |
| P2 | `p_vector_mono` — coordinatewise p-monotonicity from P1 | **proved** | pointwise P1. |
| P3 | `connection_event_increasing` — `{i↔j}` upward-closed | **proved** | `SimpleGraph.Reachable.mono`. |
| P4 | `connProb_mono` — A1 applied to the increasing event P3 | **proved** | `comparison_inequality` + P2 + P3. |
| P5 | `spin_agreement_mono` (= `hmono`) — affine∘monotone | **proved** | A2 rewrite + P4 + slope `1−1/q>0`. |
| T1 | `cutGraph_mono` — β↑/θ↓ ⇒ edge-set grows | **proved** | `le_trans` + P5. |
| T2 | `flag_homology_monotone` — edge inclusion ⇒ induced `Hₖ` map | **OBLIGATION** | your flag-complex + homology layer (no Mathlib PH). |
| T3 | `barcode_reparam_invariant` — monotone reparam ⇒ iso module | **OBLIGATION** | certifies the Cᵥ-warp is topologically free. |

One axiom you *cite* (A2). Two you *build* (T2, T3, in your own TDA layer). Everything between the
physics and the edge sets is discharged.

```lean
/-
  bridge.lean — (β, θ) thermodynamic bifiltration, physics→topology.
  AXIOMS: comparison_inequality (Grimmett 3.21), edwards_sokal_twopoint (§11.2).
  All other obligations are calculus / SimpleGraph combinatorics.
-/
import Mathlib.Combinatorics.SimpleGraph.Connectivity
import Mathlib.Analysis.SpecialFunctions.Exp
import Mathlib.Order.Monotone.Basic

universe u
variable {V : Type u} [Fintype V] [DecidableEq V]

/-! ### 0. Configurations, events, the abstract random-cluster measure. -/

/-- An edge configuration *is* its open subgraph; ordered by inclusion (`SimpleGraph`'s
    own lattice `≤`). "More open edges" = larger. -/
abbrev Config (V : Type u) := SimpleGraph V

/-- The connection event: configurations in which `i` and `j` share an open cluster. -/
def Conn (i j : V) : Set (Config V) := {ω | ω.Reachable i j}

/-- An event is *increasing* iff upward-closed under adding open edges. -/
def Increasing (A : Set (Config V)) : Prop :=
  ∀ ⦃ω ω' : Config V⦄, ω ≤ ω' → ω ∈ A → ω' ∈ A

/-- Opaque specs tying the abstract functionals to the genuine objects. These stand in
    for the unformalized measure-theoretic definitions — the axioms are stated *under*
    them, exactly as Lemma B was stated under the closed-form bond probability. -/
opaque IsRandomCluster (φ : (V → V → ℝ) → Set (Config V) → ℝ) (q : ℝ) : Prop
opaque IsPottsTwoPoint
  (spinAgree : (V → V → ℝ) → V → V → ℝ)
  (φ : (V → V → ℝ) → Set (Config V) → ℝ) (q : ℝ) : Prop

/-! ### 1. The two CITED inputs (Grimmett 2006). -/

/-- **AXIOM — Grimmett Thm 3.21 (inhomogeneous comparison, q ≥ 1).**
    Coordinatewise increase of the edge-parameter vector stochastically dominates:
    on every *increasing* event the probability cannot decrease. Bundles Holley (2.1) →
    FKG lattice / positive association (3.8, finite graphs) → comparison (3.21).
    `1 ≤ q` is the hypothesis Grimmett §3.3 shows is necessary. -/
axiom comparison_inequality
    (q : ℝ) (hq : 1 ≤ q)
    (φ : (V → V → ℝ) → Set (Config V) → ℝ) (hφ : IsRandomCluster φ q)
    {p₁ p₂ : V → V → ℝ} (hp : ∀ i j, p₁ i j ≤ p₂ i j)
    {A : Set (Config V)} (hA : Increasing A) :
    φ p₁ A ≤ φ p₂ A

/-- **AXIOM — Grimmett §11.2 (Edwards–Sokal coupling).**
    The Potts two-point function is the baseline `1/q` plus the rescaled random-cluster
    connection probability. The single algebraic fact the bridge needs from the coupling. -/
axiom edwards_sokal_twopoint
    (q : ℝ) (hq : 1 ≤ q)
    (φ : (V → V → ℝ) → Set (Config V) → ℝ)
    (spinAgree : (V → V → ℝ) → V → V → ℝ)
    (hES : IsPottsTwoPoint spinAgree φ q)
    (p : V → V → ℝ) (i j : V) :
    spinAgree p i j = 1 / q + (1 - 1 / q) * φ p (Conn i j)

/-! ### 2. Provable plumbing — the FK edge parameter and the increasing event. -/

/-- FK edge-open probability `p_e(β) = 1 − exp(−β·J_e)`. -/
def pe (J : V → V → ℝ) (β : ℝ) (i j : V) : ℝ := 1 - Real.exp (-(β * J i j))

/-- **P1.** With non-negative couplings, `p_e` is non-decreasing in β. -/
theorem pe_mono_in_beta (J : V → V → ℝ) (hJ : ∀ i j, 0 ≤ J i j) (i j : V) :
    Monotone (fun β => pe J β i j) := by
  intro β₁ β₂ hβ
  have hineq : -(β₂ * J i j) ≤ -(β₁ * J i j) := by nlinarith [hJ i j]
  have := Real.exp_le_exp.mpr hineq
  simp only [pe]; linarith

/-- **P2.** β-monotonicity, transported to the coordinatewise order on the p-vector. -/
theorem p_vector_mono (J : V → V → ℝ) (hJ : ∀ i j, 0 ≤ J i j)
    {β₁ β₂ : ℝ} (hβ : β₁ ≤ β₂) : ∀ i j, pe J β₁ i j ≤ pe J β₂ i j :=
  fun i j => pe_mono_in_beta J hJ i j hβ

/-- **P3.** Connectivity is upward-closed: adding open edges cannot break a path. -/
theorem connection_event_increasing (i j : V) : Increasing (Conn i j) := by
  intro ω ω' hle hmem
  exact hmem.mono hle   -- SimpleGraph.Reachable.mono

/-! ### 3. Composition → temperature-monotonicity of spin agreement (`hmono`). -/

/-- **P4.** Connection probability is non-decreasing in β — A1 fired on the P3 event. -/
theorem connProb_mono
    (J : V → V → ℝ) (hJ : ∀ i j, 0 ≤ J i j)
    (q : ℝ) (hq : 1 ≤ q)
    (φ : (V → V → ℝ) → Set (Config V) → ℝ) (hφ : IsRandomCluster φ q)
    {β₁ β₂ : ℝ} (hβ : β₁ ≤ β₂) (i j : V) :
    φ (pe J β₁) (Conn i j) ≤ φ (pe J β₂) (Conn i j) :=
  comparison_inequality q hq φ hφ (p_vector_mono J hJ hβ) (connection_event_increasing i j)

/-- **P5 = `hmono`.** The estimand the bifiltration needs monotone: ferromagnetic
    (`J ≥ 0`, `q > 1`) spin agreement is non-decreasing as the system cools (β ↑). -/
theorem spin_agreement_mono
    (J : V → V → ℝ) (hJ : ∀ i j, 0 ≤ J i j)
    (q : ℝ) (hq : 1 < q)
    (φ : (V → V → ℝ) → Set (Config V) → ℝ) (hφ : IsRandomCluster φ q)
    (spinAgree : (V → V → ℝ) → V → V → ℝ) (hES : IsPottsTwoPoint spinAgree φ q)
    {β₁ β₂ : ℝ} (hβ : β₁ ≤ β₂) (i j : V) :
    spinAgree (pe J β₁) i j ≤ spinAgree (pe J β₂) i j := by
  rw [edwards_sokal_twopoint q (le_of_lt hq) φ spinAgree hES (pe J β₁) i j,
      edwards_sokal_twopoint q (le_of_lt hq) φ spinAgree hES (pe J β₂) i j]
  have hslope : 0 < 1 - 1 / q := by
    have : 1 / q < 1 := by
      rw [div_lt_one (lt_trans one_pos hq)]; exact hq
    linarith
  have hc := connProb_mono J hJ q (le_of_lt hq) φ hφ hβ i j
  nlinarith [hc, hslope]

/-! ### 4. The filtration bridge: edge-set monotone, then the homology functor. -/

/-- The cut graph at (β, θ): edges whose spin agreement clears the threshold. Requires
    `spinAgree` symmetric (⟨δ_{ij}⟩ = ⟨δ_{ji}⟩). -/
def cutGraph
    (J : V → V → ℝ) (spinAgree : (V → V → ℝ) → V → V → ℝ)
    (hsym : ∀ p i j, spinAgree p i j = spinAgree p j i)
    (β θ : ℝ) : SimpleGraph V where
  Adj i j := i ≠ j ∧ θ ≤ spinAgree (pe J β) i j
  symm := by
    intro i j h; exact ⟨h.1.symm, by rw [hsym]; exact h.2⟩
  loopless := by intro i h; exact h.1 rfl

/-- **T1.** Cooling grows the complex: `β₁ ≤ β₂ ⇒ K(β₁,θ) ⊆ K(β₂,θ)`. The temperature
    half of the bifiltration's monotonicity — rides entirely on `spin_agreement_mono`.
    (The threshold half, `θ₂ ≤ θ₁ ⇒ K(β,θ₁) ⊆ K(β,θ₂)`, is `le_trans` alone — omitted.) -/
theorem cutGraph_mono_beta
    (J : V → V → ℝ) (hJ : ∀ i j, 0 ≤ J i j)
    (q : ℝ) (hq : 1 < q)
    (φ : (V → V → ℝ) → Set (Config V) → ℝ) (hφ : IsRandomCluster φ q)
    (spinAgree : (V → V → ℝ) → V → V → ℝ) (hES : IsPottsTwoPoint spinAgree φ q)
    (hsym : ∀ p i j, spinAgree p i j = spinAgree p j i)
    {β₁ β₂ : ℝ} (hβ : β₁ ≤ β₂) (θ : ℝ) :
    cutGraph J spinAgree hsym β₁ θ ≤ cutGraph J spinAgree hsym β₂ θ := by
  intro i j hadj
  exact ⟨hadj.1, le_trans hadj.2 (spin_agreement_mono J hJ q hq φ hφ spinAgree hES hβ i j)⟩

/-! ### 5. OBLIGATIONS in your TDA layer (no Mathlib PH — build, don't cite). -/

/-- Abstract degree-`k` homology of the flag (clique) complex of a graph, valued in a
    field 𝔽 as a vector space. Stand-in for `Tda.Primitives.FlagComplex` + `PersistentHomology`. -/
opaque FlagHomology (k : ℕ) : SimpleGraph V → Type u

/-- **T2 — functorial monotonicity (clique functor + induced map).** A graph inclusion
    induces a linear map on flag-complex homology. The structural soundness of the
    filtration: every `(β,θ) ⪯ (β',θ')` yields an arrow in the persistence module.
    OBLIGATION — proved once the flag complex + simplicial boundary maps are defined; no
    Mathlib persistence layer to cite. -/
theorem flag_homology_monotone (k : ℕ) {G G' : SimpleGraph V} (h : G ≤ G') :
    FlagHomology k G → FlagHomology k G' := by
  sorry  -- [OBLIGATION] clique-functor + induced boundary map

/-- **T3 — reparameterization invariance (the Cᵥ-warp is free).** A strictly monotone
    reparameterization of the filtration index induces an *isomorphism* of persistence
    modules — so warping the slice speed to spend resolution at `T_c` cannot change the
    barcode, only the axis labels. The formal certificate that Cᵥ-calibration is a
    discretization choice, not a topological one.
    OBLIGATION — once `PersistenceModule`/barcode are defined in your layer. -/
theorem barcode_reparam_invariant
    (h : ℝ → ℝ) (hmono : StrictMono h)
    (K : ℝ → SimpleGraph V) (k : ℕ) :
    True := by   -- placeholder codomain: "PersistenceModule (K ∘ h) ≅ PersistenceModule K"
  sorry  -- [OBLIGATION] persistence module iso under monotone reindex
```

## Notes for the proof pass

- **The honest count is 3**, and they are not equal in kind: **A2 is a citation** (Grimmett Thm 1.16
  and (11.25)). A1 was one too, and is now proved: Mathlib has Holley's inequality, so it did not
  need a project of its own. **T2, T3 are construction obligations** in your TDA layer, not
  Mathlib gaps — they wait on `FlagComplex`/`PersistentHomology` getting Lean definitions.
  The physics→edge-set spine (P1–P5, T1) is **fully discharged** given A2.
- **The inhomogeneous-comparison check** lives inside A1: `comparison_inequality` is stated for
  a *vector* `p : V → V → ℝ`, which is the version SPC needs (heterogeneous curved-metric
  couplings). Checked 2026-11-04: Grimmett 3.8 and 3.21 are written for one `p`, and §11.4 says the
  comparison inequalities stay valid with parameters per edge, there on a finite box of a lattice.
  `comparison_inequality` now proves it on any finite graph with parameters per edge, so the citation
  is no longer needed.
- **`hmono` is the whole point.** Once `spin_agreement_mono` (P5) is green, the noisy-estimate
  monotonization (isotonic / monotone-BARS) is justified *not* as a hack but as recovering the
  property this file proves the estimand has — structure from outside the estimator
  (validation-independence).
- `q > 1` (P5/T1) vs `q ≥ 1` (A1): the comparison inequality needs `q ≥ 1`; the affine slope
  `1 − 1/q > 0` needs `q > 1`. Integer Potts (`q ≥ 2`) satisfies both. `q = 1` is percolation —
  the slope degenerates and "spin agreement" is pure connection probability, which is fine but
  a different object.
