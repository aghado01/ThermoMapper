import Mathlib.Combinatorics.SimpleGraph.Basic
import Mathlib.Data.Real.Basic

/-!
# Cumulative fields — graph cuts

Threshold a symmetric edge-score field on an explicit base graph. Keeping the
base graph in the definition makes the construction applicable to sparse graph
artifacts rather than silently completing the vertex set.
-/

namespace CumulativeField

variable {V : Type*}

/-- Retain exactly the base-graph edges whose symmetric score exceeds `tau`. -/
def cutGraph (base : SimpleGraph V) (score : V → V → ℝ) (tau : ℝ)
    (hsymm : ∀ u v, score u v = score v u) : SimpleGraph V where
  Adj u v := base.Adj u v ∧ tau < score u v
  symm := by
    constructor
    intro u v huv
    exact ⟨base.adj_symm huv.1, by simpa only [hsymm u v] using huv.2⟩
  loopless := by
    constructor
    intro u huu
    exact base.loopless.irrefl u huu.1

/-- Raising the cut can only remove edges. -/
theorem cutGraph_antitone (base : SimpleGraph V) (score : V → V → ℝ)
    (hsymm : ∀ u v, score u v = score v u) {tau₁ tau₂ : ℝ}
    (hcut : tau₁ ≤ tau₂) :
    cutGraph base score tau₂ hsymm ≤ cutGraph base score tau₁ hsymm := by
  intro u v huv
  exact ⟨huv.1, lt_of_le_of_lt hcut huv.2⟩

/-- Increasing the score field can only add edges at a fixed cut. -/
theorem cutGraph_mono_score (base : SimpleGraph V)
    (lower upper : V → V → ℝ) (tau : ℝ)
    (hlower : ∀ u v, lower u v = lower v u)
    (hupper : ∀ u v, upper u v = upper v u)
    (hscore : ∀ u v, lower u v ≤ upper u v) :
    cutGraph base lower tau hlower ≤ cutGraph base upper tau hupper := by
  intro u v huv
  exact ⟨huv.1, lt_of_lt_of_le huv.2 (hscore u v)⟩

end CumulativeField
