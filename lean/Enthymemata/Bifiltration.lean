import Lemmas.RandomCluster
import Mathlib.Tactic.Linarith

/-!
# The (β, θ) thermodynamic bifiltration — physics → topology bridge

Cooling a ferromagnetic Potts model only adds edges to its spin-agreement cut graph, and so does
lowering the threshold. So the exact edge family is monotone on both axes and its clusters nest.

`Lemmas.RandomCluster` proves the same for connection probability, comparison inequality included.
Spin agreement needs one more result, cited here as an axiom over an opaque spec of the Potts
two-point function, and everything else is proved.

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

/-! ### The cited input -/

/-- Stands in for: `spinAgree p i j` is the probability that `i` and `j` agree under the Potts model
coupled to the random-cluster measure `φ p`. The Potts model exists only at integer `q ≥ 2`, so this
can hold only there. -/
opaque IsPottsTwoPoint [Fintype V]
    (spinAgree : (V → V → ℝ) → V → V → ℝ)
    (φ : (V → V → ℝ) → Set (Config V) → ℝ) (q : ℝ) : Prop

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

/-! ### Composition: monotonicity in temperature -/

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
