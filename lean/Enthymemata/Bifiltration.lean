import Mathlib.Analysis.SpecialFunctions.Exp
import Mathlib.Combinatorics.SimpleGraph.Connectivity.Connected
import Mathlib.Order.Monotone.Basic
import Mathlib.Tactic.Linarith

/-!
# The (β, θ) thermodynamic bifiltration — physics → topology bridge

Cooling a ferromagnetic Potts model only adds edges to its cut graph, and so does lowering the
threshold. So the exact edge family is monotone on both axes and its clusters nest.

Two results are cited as axioms, over opaque specs of the random-cluster measure and the Potts
two-point function. Everything between them and the edge sets is proved here.

* `comparison_inequality` — Grimmett 2006, Thm 3.21, with per-edge parameters as in his §11.4.
* `edwards_sokal_twopoint` — Grimmett 2006, Thm 1.16, with per-edge parameters as in his (11.25).

All of this is about the exact expectation. A finite-sample estimate need not be monotone, so
nothing here says a sampled stack nests. It says a monotonizer recovers a property the target has.

Protolemmata:

* `lean/Protolemmata/bifiltration-bridge-lemmas.md` — the split between axioms and obligations.
* `lean/Protolemmata/fable-bifiltration-lemmas.md` — the corrections: FK rather than Griffiths, the
  ideal expectation against its estimate, and the C_v warp against the slope.
* `lean/Protolemmata/fable-potts-griffiths-lemmas.md` — the theorem chain through Grimmett.
-/

universe u

namespace Bifiltration

variable {V : Type u}

/-! ### Configurations and events -/

/-- An edge configuration is its open subgraph, ordered by inclusion. -/
abbrev Config (V : Type u) := SimpleGraph V

/-- The event that `i` and `j` share an open cluster. -/
def Conn (i j : V) : Set (Config V) := {ω | ω.Reachable i j}

/-- An event is increasing if adding open edges keeps it. -/
def Increasing (A : Set (Config V)) : Prop :=
  ∀ ⦃ω ω' : Config V⦄, ω ≤ ω' → ω ∈ A → ω' ∈ A

/-- The parameters both cited theorems cover: symmetric, and in `[0, 1)` on every pair. -/
structure IsEdgeParameter (p : V → V → ℝ) : Prop where
  nonneg : ∀ i j, 0 ≤ p i j
  lt_one : ∀ i j, p i j < 1
  symm : ∀ i j, p i j = p j i

/-! ### The two cited inputs -/

/-- Stands in for: `φ p A` is the probability of `A` under the random-cluster measure with
parameters `p` and `q` on the finite graph. -/
opaque IsRandomCluster [Fintype V] (φ : (V → V → ℝ) → Set (Config V) → ℝ) (q : ℝ) : Prop

/-- Stands in for: `spinAgree p i j` is the probability that `i` and `j` agree under the Potts model
coupled to the random-cluster measure `φ p`. The Potts model exists only at integer `q ≥ 2`, so this
can hold only there. -/
opaque IsPottsTwoPoint [Fintype V]
    (spinAgree : (V → V → ℝ) → V → V → ℝ)
    (φ : (V → V → ℝ) → Set (Config V) → ℝ) (q : ℝ) : Prop

/-- **A1, comparison.** Raising edge parameters cannot lower the probability of an increasing event
(Grimmett Thm 3.21, through positive association, Thm 3.8, and Holley, Thm 2.1). He states it for
one `p` on any finite graph. §11.4 says it stays valid with parameters per edge, there on a box of a
lattice, and the proof uses neither homogeneity nor the lattice: `p` enters only through a product
over edges. `1 ≤ q` carries the positive association, which fails in general below 1. -/
axiom comparison_inequality [Fintype V]
    (q : ℝ) (hq : 1 ≤ q)
    (φ : (V → V → ℝ) → Set (Config V) → ℝ) (hφ : IsRandomCluster φ q)
    {p₁ p₂ : V → V → ℝ} (hp₁ : IsEdgeParameter p₁) (hp₂ : IsEdgeParameter p₂)
    (hp : ∀ i j, p₁ i j ≤ p₂ i j)
    {A : Set (Config V)} (hA : Increasing A) :
    φ p₁ A ≤ φ p₂ A

/-- **A2, two-point identity.** Spin agreement is the baseline `1 / q` plus the rescaled connection
probability (Grimmett Thm 1.16; with parameters per edge, (11.25), where `p_e = 1 - exp (-β J_e)`,
which is `pe`). -/
axiom edwards_sokal_twopoint [Fintype V]
    (q : ℝ) (hq : 1 ≤ q)
    (φ : (V → V → ℝ) → Set (Config V) → ℝ)
    (spinAgree : (V → V → ℝ) → V → V → ℝ)
    (hES : IsPottsTwoPoint spinAgree φ q)
    (p : V → V → ℝ) (hp : IsEdgeParameter p) (i j : V) :
    spinAgree p i j = 1 / q + (1 - 1 / q) * φ p (Conn i j)

/-! ### The FK edge parameter and the increasing event -/

/-- FK edge-open probability `p_e(β) = 1 - exp (-β J_e)`. -/
noncomputable def pe (J : V → V → ℝ) (β : ℝ) (i j : V) : ℝ := 1 - Real.exp (-(β * J i j))

/-- **P1.** With non-negative couplings, `p_e` is non-decreasing in β. -/
theorem pe_mono_in_beta (J : V → V → ℝ) (hJ : ∀ i j, 0 ≤ J i j) (i j : V) :
    Monotone (fun β => pe J β i j) := by
  intro β₁ β₂ hβ
  have hineq : -(β₂ * J i j) ≤ -(β₁ * J i j) := by nlinarith [hJ i j]
  have := Real.exp_le_exp.mpr hineq
  simp only [pe]; linarith

/-- **P2.** β-monotonicity in the coordinatewise order on the parameters. -/
theorem p_vector_mono (J : V → V → ℝ) (hJ : ∀ i j, 0 ≤ J i j)
    {β₁ β₂ : ℝ} (hβ : β₁ ≤ β₂) : ∀ i j, pe J β₁ i j ≤ pe J β₂ i j :=
  fun i j => pe_mono_in_beta J hJ i j hβ

/-- Symmetric non-negative couplings at `β ≥ 0` give parameters the cited theorems cover. -/
theorem pe_isEdgeParameter (J : V → V → ℝ) (hJ : ∀ i j, 0 ≤ J i j) (hJs : ∀ i j, J i j = J j i)
    {β : ℝ} (hβ : 0 ≤ β) : IsEdgeParameter (pe J β) where
  nonneg i j := by
    have hle : Real.exp (-(β * J i j)) ≤ 1 :=
      Real.exp_le_one_iff.mpr (neg_nonpos.mpr (mul_nonneg hβ (hJ i j)))
    simp only [pe]; linarith
  lt_one i j := by
    have := Real.exp_pos (-(β * J i j))
    simp only [pe]; linarith
  symm i j := by
    simp only [pe]; rw [hJs i j]

/-- **P3.** Adding open edges cannot break a path. -/
theorem connection_event_increasing (i j : V) : Increasing (Conn i j) :=
  fun _ _ hle hmem => SimpleGraph.Reachable.mono hle hmem

/-! ### Composition: monotonicity in temperature -/

/-- **P4.** Connection probability is non-decreasing in β: A1 on the event of P3. -/
theorem connProb_mono [Fintype V]
    (J : V → V → ℝ) (hJ : ∀ i j, 0 ≤ J i j) (hJs : ∀ i j, J i j = J j i)
    (q : ℝ) (hq : 1 ≤ q)
    (φ : (V → V → ℝ) → Set (Config V) → ℝ) (hφ : IsRandomCluster φ q)
    {β₁ β₂ : ℝ} (hβ₁ : 0 ≤ β₁) (hβ : β₁ ≤ β₂) (i j : V) :
    φ (pe J β₁) (Conn i j) ≤ φ (pe J β₂) (Conn i j) :=
  comparison_inequality q hq φ hφ (pe_isEdgeParameter J hJ hJs hβ₁)
    (pe_isEdgeParameter J hJ hJs (hβ₁.trans hβ)) (p_vector_mono J hJ hβ)
    (connection_event_increasing i j)

/-- **P5.** Ferromagnetic spin agreement is non-decreasing as the system cools. `1 < q` gives the
positive slope `1 - 1 / q`. -/
theorem spin_agreement_mono [Fintype V]
    (J : V → V → ℝ) (hJ : ∀ i j, 0 ≤ J i j) (hJs : ∀ i j, J i j = J j i)
    (q : ℝ) (hq : 1 < q)
    (φ : (V → V → ℝ) → Set (Config V) → ℝ) (hφ : IsRandomCluster φ q)
    (spinAgree : (V → V → ℝ) → V → V → ℝ) (hES : IsPottsTwoPoint spinAgree φ q)
    {β₁ β₂ : ℝ} (hβ₁ : 0 ≤ β₁) (hβ : β₁ ≤ β₂) (i j : V) :
    spinAgree (pe J β₁) i j ≤ spinAgree (pe J β₂) i j := by
  rw [edwards_sokal_twopoint q hq.le φ spinAgree hES (pe J β₁)
        (pe_isEdgeParameter J hJ hJs hβ₁) i j,
      edwards_sokal_twopoint q hq.le φ spinAgree hES (pe J β₂)
        (pe_isEdgeParameter J hJ hJs (hβ₁.trans hβ)) i j]
  have hslope : 0 < 1 - 1 / q := by
    have : 1 / q < 1 := by
      rw [div_lt_one (lt_trans one_pos hq)]; exact hq
    linarith
  have hc := connProb_mono J hJ hJs q hq.le φ hφ hβ₁ hβ i j
  nlinarith [hc, hslope]

/-! ### The cut graph: monotonicity on both axes -/

/-- The cut graph at `(β, θ)`: edges whose spin agreement clears the threshold. -/
def cutGraph
    (J : V → V → ℝ) (spinAgree : (V → V → ℝ) → V → V → ℝ)
    (hsym : ∀ p i j, spinAgree p i j = spinAgree p j i)
    (β θ : ℝ) : SimpleGraph V where
  Adj i j := i ≠ j ∧ θ ≤ spinAgree (pe J β) i j
  symm := ⟨fun i j h => ⟨h.1.symm, by rw [hsym]; exact h.2⟩⟩
  loopless := ⟨fun i h => h.1 rfl⟩

/-- **T1, temperature half.** Cooling grows the cut graph. -/
theorem cutGraph_mono_beta [Fintype V]
    (J : V → V → ℝ) (hJ : ∀ i j, 0 ≤ J i j) (hJs : ∀ i j, J i j = J j i)
    (q : ℝ) (hq : 1 < q)
    (φ : (V → V → ℝ) → Set (Config V) → ℝ) (hφ : IsRandomCluster φ q)
    (spinAgree : (V → V → ℝ) → V → V → ℝ) (hES : IsPottsTwoPoint spinAgree φ q)
    (hsym : ∀ p i j, spinAgree p i j = spinAgree p j i)
    {β₁ β₂ : ℝ} (hβ₁ : 0 ≤ β₁) (hβ : β₁ ≤ β₂) (θ : ℝ) :
    cutGraph J spinAgree hsym β₁ θ ≤ cutGraph J spinAgree hsym β₂ θ :=
  fun i j hadj =>
    ⟨hadj.1, le_trans hadj.2 (spin_agreement_mono J hJ hJs q hq φ hφ spinAgree hES hβ₁ hβ i j)⟩

/-- **T1, threshold half.** Lowering the threshold grows the cut graph. No physics is needed. -/
theorem cutGraph_anti_theta
    (J : V → V → ℝ) (spinAgree : (V → V → ℝ) → V → V → ℝ)
    (hsym : ∀ p i j, spinAgree p i j = spinAgree p j i)
    (β : ℝ) {θ₁ θ₂ : ℝ} (hθ : θ₂ ≤ θ₁) :
    cutGraph J spinAgree hsym β θ₁ ≤ cutGraph J spinAgree hsym β θ₂ :=
  fun _ _ hadj => ⟨hadj.1, le_trans hθ hadj.2⟩

/-- **T1.** The cut graphs grow along the product order on `(β, -θ)`, so they are a bifiltration. -/
theorem cutGraph_mono [Fintype V]
    (J : V → V → ℝ) (hJ : ∀ i j, 0 ≤ J i j) (hJs : ∀ i j, J i j = J j i)
    (q : ℝ) (hq : 1 < q)
    (φ : (V → V → ℝ) → Set (Config V) → ℝ) (hφ : IsRandomCluster φ q)
    (spinAgree : (V → V → ℝ) → V → V → ℝ) (hES : IsPottsTwoPoint spinAgree φ q)
    (hsym : ∀ p i j, spinAgree p i j = spinAgree p j i)
    {β₁ β₂ : ℝ} (hβ₁ : 0 ≤ β₁) (hβ : β₁ ≤ β₂) {θ₁ θ₂ : ℝ} (hθ : θ₂ ≤ θ₁) :
    cutGraph J spinAgree hsym β₁ θ₁ ≤ cutGraph J spinAgree hsym β₂ θ₂ :=
  (cutGraph_mono_beta J hJ hJs q hq φ hφ spinAgree hES hsym hβ₁ hβ θ₁).trans
    (cutGraph_anti_theta J spinAgree hsym β₂ hθ)

/-- The clusters nest: two points together at `(β₁, θ₁)` are together at every colder or
lower-threshold corner, so each partition refines the next. -/
theorem cutGraph_reachable_mono [Fintype V]
    (J : V → V → ℝ) (hJ : ∀ i j, 0 ≤ J i j) (hJs : ∀ i j, J i j = J j i)
    (q : ℝ) (hq : 1 < q)
    (φ : (V → V → ℝ) → Set (Config V) → ℝ) (hφ : IsRandomCluster φ q)
    (spinAgree : (V → V → ℝ) → V → V → ℝ) (hES : IsPottsTwoPoint spinAgree φ q)
    (hsym : ∀ p i j, spinAgree p i j = spinAgree p j i)
    {β₁ β₂ : ℝ} (hβ₁ : 0 ≤ β₁) (hβ : β₁ ≤ β₂) {θ₁ θ₂ : ℝ} (hθ : θ₂ ≤ θ₁) {i j : V}
    (h : (cutGraph J spinAgree hsym β₁ θ₁).Reachable i j) :
    (cutGraph J spinAgree hsym β₂ θ₂).Reachable i j :=
  h.mono (cutGraph_mono J hJ hJs q hq φ hφ spinAgree hES hsym hβ₁ hβ hθ)

/-! ### Obligations, not yet stated

T2 and T3 stay in the protolemma until their codomains are defined: a graph inclusion induces a map
on flag-complex homology (`flag_homology_monotone`), and a strictly monotone reindexing gives an
isomorphic persistence module (`barcode_reparam_invariant`). `Lemmas.IntervalDecomposition` now has
filtered complexes and interval modules to state them against. -/

end Bifiltration
