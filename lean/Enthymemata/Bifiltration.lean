import Lemmas.EdwardsSokal
import Mathlib.Tactic.Linarith

/-!
# The (β, θ) thermodynamic bifiltration — physics → topology bridge

Cooling a ferromagnetic Potts model only adds edges to its spin-agreement cut graph, and so does
lowering the threshold. So the exact edge family is monotone on both axes and its clusters nest.

`Lemmas.RandomCluster` proves the same for connection probability, comparison inequality included,
and `Lemmas.EdwardsSokal` gives the identity that carries it to spin agreement. Nothing is left to
cite. `potts_agreement_mono` is the result with every hypothesis discharged.

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

/-! ### Composition: monotonicity in temperature -/

/-- **P5.** Ferromagnetic spin agreement is non-decreasing as the system cools. `1 < q` gives the
positive slope `1 - 1 / q`. -/
theorem spin_agreement_mono [Fintype V]
    (J : V → V → ℝ) (hJ : ∀ i j, 0 ≤ J i j) (hJs : ∀ i j, J i j = J j i)
    (q : ℝ) (hq : 1 < q)
    (φ : (V → V → ℝ) → Set (Config V) → ℝ)
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
  have hc := connProb_mono J hJ hJs q hq.le φ hES.1 hβ₁ hβ i j
  nlinarith [hc, hslope]

/-- P5 with nothing left to discharge: in the Potts model on any finite graph with symmetric
non-negative couplings, `i` and `j` agree at least as often at a lower temperature. `pottsWeight_pe`
says this is the model with Hamiltonian `∑_e J_e (1 - δ_e)`. -/
theorem potts_agreement_mono [Fintype V]
    (J : V → V → ℝ) (hJ : ∀ i j, 0 ≤ J i j) (hJs : ∀ i j, J i j = J j i)
    {n : ℕ} (hn : 2 ≤ n) {β₁ β₂ : ℝ} (hβ₁ : 0 ≤ β₁) (hβ : β₁ ≤ β₂) (i j : V) :
    pottsAgree (pe_isEdgeParameter J hJ hJs hβ₁) n i j ≤
      pottsAgree (pe_isEdgeParameter J hJ hJs (hβ₁.trans hβ)) n i j := by
  have hq : (1 : ℝ) < n := by exact_mod_cast hn
  have h := spin_agreement_mono J hJ hJs n hq (randomCluster n) (pottsTwoPoint n)
    (isPottsTwoPoint_pottsTwoPoint n) hβ₁ hβ i j
  rwa [pottsTwoPoint, pottsTwoPoint, dif_pos (pe_isEdgeParameter J hJ hJs hβ₁),
    dif_pos (pe_isEdgeParameter J hJ hJs (hβ₁.trans hβ))] at h

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
    (φ : (V → V → ℝ) → Set (Config V) → ℝ)
    (spinAgree : (V → V → ℝ) → V → V → ℝ) (hES : IsPottsTwoPoint spinAgree φ q)
    (hsym : ∀ p i j, spinAgree p i j = spinAgree p j i)
    {β₁ β₂ : ℝ} (hβ₁ : 0 ≤ β₁) (hβ : β₁ ≤ β₂) (θ : ℝ) :
    cutGraph J spinAgree hsym β₁ θ ≤ cutGraph J spinAgree hsym β₂ θ :=
  fun i j hadj =>
    ⟨hadj.1, le_trans hadj.2 (spin_agreement_mono J hJ hJs q hq φ spinAgree hES hβ₁ hβ i j)⟩

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
    (φ : (V → V → ℝ) → Set (Config V) → ℝ)
    (spinAgree : (V → V → ℝ) → V → V → ℝ) (hES : IsPottsTwoPoint spinAgree φ q)
    (hsym : ∀ p i j, spinAgree p i j = spinAgree p j i)
    {β₁ β₂ : ℝ} (hβ₁ : 0 ≤ β₁) (hβ : β₁ ≤ β₂) {θ₁ θ₂ : ℝ} (hθ : θ₂ ≤ θ₁) :
    cutGraph J spinAgree hsym β₁ θ₁ ≤ cutGraph J spinAgree hsym β₂ θ₂ :=
  (cutGraph_mono_beta J hJ hJs q hq φ spinAgree hES hsym hβ₁ hβ θ₁).trans
    (cutGraph_anti_theta J spinAgree hsym β₂ hθ)

/-- The clusters nest: two points together at `(β₁, θ₁)` are together at every colder or
lower-threshold corner, so each partition refines the next. -/
theorem cutGraph_reachable_mono [Fintype V]
    (J : V → V → ℝ) (hJ : ∀ i j, 0 ≤ J i j) (hJs : ∀ i j, J i j = J j i)
    (q : ℝ) (hq : 1 < q)
    (φ : (V → V → ℝ) → Set (Config V) → ℝ)
    (spinAgree : (V → V → ℝ) → V → V → ℝ) (hES : IsPottsTwoPoint spinAgree φ q)
    (hsym : ∀ p i j, spinAgree p i j = spinAgree p j i)
    {β₁ β₂ : ℝ} (hβ₁ : 0 ≤ β₁) (hβ : β₁ ≤ β₂) {θ₁ θ₂ : ℝ} (hθ : θ₂ ≤ θ₁) {i j : V}
    (h : (cutGraph J spinAgree hsym β₁ θ₁).Reachable i j) :
    (cutGraph J spinAgree hsym β₂ θ₂).Reachable i j :=
  h.mono (cutGraph_mono J hJ hJs q hq φ spinAgree hES hsym hβ₁ hβ hθ)

/-! ### Obligations, not yet stated

T2 and T3 stay in the protolemma until their codomains are defined: a graph inclusion induces a map
on flag-complex homology (`flag_homology_monotone`), and a strictly monotone reindexing gives an
isomorphic persistence module (`barcode_reparam_invariant`). `Lemmas.IntervalDecomposition` now has
filtered complexes and interval modules to state them against. -/

end Bifiltration
