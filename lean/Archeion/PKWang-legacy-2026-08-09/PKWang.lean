import Mathlib.Analysis.SpecialFunctions.Log.Basic
import Mathlib.Combinatorics.SimpleGraph.Basic
import Mathlib.Data.Fintype.Basic
import Mathlib.Tactic.NormNum

/-!
# PKWang — threshold reduction and the deterministic cut

The proved core of Lemma B behind the `PKWang` solver (Wang 2020): in the
`M → ∞` limit the exponential-draw apparatus reduces to the deterministic
threshold `H_cum > T · ln 2`, so the solver is thermal single-linkage — cut
every edge whose cumulative energy fails the threshold, read clusters off
connected components.

Open obligations (Lemma A pooling, single-linkage equivalence) apologize in
`Enthymemata.PKWangB`. Protolemma: `lean/Protolemmata/spc-pkwang-lemmas.md`.
-/

open Real

namespace PKWang

/-! ## Lemma B — analytical reduction and threshold equivalence -/

/-- The closed-form survival probability `1 - exp (-H_cum / T)` exceeds `1/2`
iff the cumulative energy exceeds `T * log 2`. This is the algebraic heart of
Lemma B: the Monte Carlo draw thresholded at `1/2` is a deterministic cut. -/
theorem pk_wang_closed_form_reduction (Hcum T : ℝ) (hT : 0 < T) :
    1 - exp (-(Hcum / T)) > 1 / 2 ↔ Hcum > T * log 2 := by
  rw [gt_iff_lt, gt_iff_lt, lt_sub_comm,
    show (1 : ℝ) - 1 / 2 = 1 / 2 by norm_num,
    ← lt_log_iff_exp_lt (by norm_num : (0 : ℝ) < 1 / 2),
    show Real.log (1 / 2) = -Real.log 2 by rw [one_div, Real.log_inv],
    neg_lt_neg_iff, lt_div_iff₀ hT, mul_comm]

variable {V : Type*} [Fintype V]

/-- The graph formed by the deterministic cut: keep edge `(u, v)` iff its
cumulative energy strictly exceeds `T * log 2`. Symmetry of `H_cum` is taken
as a hypothesis (it holds for the undirected coupling field). -/
def deterministicCutGraph (Hcum : V → V → ℝ) (T : ℝ)
    (hsymm : ∀ u v, Hcum u v = Hcum v u) : SimpleGraph V where
  Adj u v := u ≠ v ∧ T * log 2 < Hcum u v
  symm := by
    constructor
    intro u v ⟨hne, h⟩
    exact ⟨hne.symm, by rw [hsymm v u]; exact h⟩
  loopless := by
    constructor
    intro u h
    exact h.1 rfl

end PKWang
