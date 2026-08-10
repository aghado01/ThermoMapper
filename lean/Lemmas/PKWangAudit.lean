import Mathlib.Analysis.SpecialFunctions.Log.Basic
import Mathlib.Data.Fintype.Basic
import Lemmas.CumulativeField.Cut

/-!
# PKWang audit — verified correction fragments

Paper-specific results for the 2020 parallel-SPC proposal. This is the sole
active Lean module that carries the paper name; the recovered mathematics
lives under `CumulativeField` and does not depend on this audit.

The complete critical narrative and remaining obligations are specified in
`Protolemmata/PKWangAudit.md`.
-/

open Real

namespace PKWangAudit

/-- Boundary between adjacent cumulative totals under the paper's literal
"closest cumulative energy" matching rule. -/
noncomputable def closestBoundary (previous current : ℝ) : ℝ :=
  (previous + current) / 2

/-- The exponential surrogate assigns positive survival probability beyond
every finite energy bound. A finite Potts Hamiltonian has no such tail. -/
theorem exponentialSurrogate_has_positive_tail (Hmax T : ℝ) :
    0 < exp (-(Hmax / T)) :=
  exp_pos _

/-- A pairwise mask with `01` and `12` equal but `02` unequal cannot be
realized by any spin labeling: equality is transitive. -/
theorem nontransitiveMask_not_spinRealizable :
    ¬ ∃ spin : Fin 3 → Nat,
      spin 0 = spin 1 ∧ spin 1 = spin 2 ∧ spin 0 ≠ spin 2 := by
  rintro ⟨spin, h01, h12, h02⟩
  exact h02 (h01.trans h12)

/-- The paper's finite-draw affinity is an empirical CDF. Thresholding it is
therefore exactly a comparison with the count below the boundary; no generated
pairwise masks are needed to perform this reduction. -/
noncomputable def empiricalAffinity {M : ℕ} (sample : Fin M → ℝ)
    (boundary : ℝ) : ℝ := by
  classical
  exact ((Finset.univ.filter (fun k => sample k < boundary)).card : ℝ) / M

theorem empiricalAffinity_gt_iff_count {M : ℕ} (hM : 0 < M)
    (sample : Fin M → ℝ) (boundary theta : ℝ) :
    empiricalAffinity sample boundary > theta ↔
      theta * (M : ℝ) <
        ((Finset.univ.filter (fun k => sample k < boundary)).card : ℝ) := by
  classical
  unfold empiricalAffinity
  rw [gt_iff_lt, lt_div_iff₀ (Nat.cast_pos.mpr hM)]

/-- At the population level the paper's fixed `0.5` threshold is the familiar
`T * log 2` cut, applied to whichever boundary its matching rule supplies. -/
theorem surrogate_half_cut (boundary T : ℝ) (hT : 0 < T) :
    CumulativeField.affinity boundary T > 1 / 2 ↔
      boundary > T * log 2 :=
  CumulativeField.affinity_gt_half_iff boundary T hT

end PKWangAudit
