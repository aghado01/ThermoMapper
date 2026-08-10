import Mathlib.Data.Real.Basic
import Mathlib.Tactic.Linarith

/-!
# Cumulative fields — endpoint symmetrization

Three parameter-free reconciliations of two directed endpoint scores. They are
defined in score space, before any nonlinear affinity link.
-/

namespace CumulativeField

/-- Bilateral evidence: the weaker endpoint score controls. -/
noncomputable def mutualAggregate (a b : ℝ) : ℝ := min a b

/-- Balanced evidence: arithmetic mean in cumulative-score (hazard) space. -/
noncomputable def hazardMean (a b : ℝ) : ℝ := (a + b) / 2

/-- Unilateral evidence: the stronger endpoint score controls. -/
noncomputable def inclusiveAggregate (a b : ℝ) : ℝ := max a b

theorem mutualAggregate_le_hazardMean (a b : ℝ) :
    mutualAggregate a b ≤ hazardMean a b := by
  dsimp [mutualAggregate, hazardMean]
  linarith [min_le_left a b, min_le_right a b]

theorem hazardMean_le_inclusiveAggregate (a b : ℝ) :
    hazardMean a b ≤ inclusiveAggregate a b := by
  dsimp [hazardMean, inclusiveAggregate]
  linarith [le_max_left a b, le_max_right a b]

variable {V : Type*}

noncomputable def mutualScore (directed : V → V → ℝ) (u v : V) : ℝ :=
  mutualAggregate (directed u v) (directed v u)

noncomputable def hazardMeanScore (directed : V → V → ℝ) (u v : V) : ℝ :=
  hazardMean (directed u v) (directed v u)

noncomputable def inclusiveScore (directed : V → V → ℝ) (u v : V) : ℝ :=
  inclusiveAggregate (directed u v) (directed v u)

theorem mutualScore_symm (directed : V → V → ℝ) (u v : V) :
    mutualScore directed u v = mutualScore directed v u := by
  simp [mutualScore, mutualAggregate, min_comm]

theorem hazardMeanScore_symm (directed : V → V → ℝ) (u v : V) :
    hazardMeanScore directed u v = hazardMeanScore directed v u := by
  simp [hazardMeanScore, hazardMean, add_comm]

theorem inclusiveScore_symm (directed : V → V → ℝ) (u v : V) :
    inclusiveScore directed u v = inclusiveScore directed v u := by
  simp [inclusiveScore, inclusiveAggregate, max_comm]

theorem mutualScore_le_hazardMeanScore (directed : V → V → ℝ) (u v : V) :
    mutualScore directed u v ≤ hazardMeanScore directed u v :=
  mutualAggregate_le_hazardMean _ _

theorem hazardMeanScore_le_inclusiveScore (directed : V → V → ℝ) (u v : V) :
    hazardMeanScore directed u v ≤ inclusiveScore directed u v :=
  hazardMean_le_inclusiveAggregate _ _

end CumulativeField
