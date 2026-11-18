import Lemmas.TimeFrequency.PhasePropagation
import Mathlib.Algebra.Ring.Parity

/-! # Recovering two contrasts from products on sampled context pairs -/

namespace TimeFrequency.LinkedProfiles
open SimpleGraph
noncomputable section

variable {V : Type*} {G : SimpleGraph V}

/-- For successive measured edge moments this reads m₀ / (m₁ / (m₂ / ...)). -/
def alternatingRatio (m : V → V → ℝ) {u v : V} : G.Walk u v → ℝ
  | .nil => 1
  | .cons (v := w) _ tail => m u w / alternatingRatio m tail

private theorem ratio_ne_zero {m : V → V → ℝ} (hm : ∀ u v, G.Adj u v → m u v ≠ 0)
    {u v : V} (w : G.Walk u v) : alternatingRatio m w ≠ 0 := by
  induction w with
  | nil => exact one_ne_zero
  | cons h _ ih => exact div_ne_zero (hm _ _ h) ih

theorem alternatingRatio_spec {m : V → V → ℝ} {a : V → ℝ}
    (hm : ∀ u v, G.Adj u v → m u v ≠ 0) (ha : ∀ u v, G.Adj u v → m u v = a u * a v)
    {u v : V} (w : G.Walk u v) :
    (Even w.length → alternatingRatio m w * a v = a u) ∧
    (Odd w.length → alternatingRatio m w = a u * a v) := by
  induction w with
  | nil => simp [alternatingRatio]
  | @cons i j k hij tail ih =>
    have hr := ratio_ne_zero hm tail
    simp only [Walk.length_cons, Nat.even_add_one, Nat.odd_add_one, Nat.not_even_iff_odd,
      Nat.not_odd_iff_even, alternatingRatio, ha i j hij]
    constructor
    · intro h
      have ht := ih.2 h
      field_simp [hr]
      rw [ht]
      ring
    · intro h
      apply (div_eq_iff hr).mpr
      rw [← ih.1 h]
      ring

/-- Once a signed root contrast is chosen, each reached contrast can be read from the walk. -/
theorem contrast_from_walk {m : V → V → ℝ} {a : V → ℝ}
    (hm : ∀ u v, G.Adj u v → m u v ≠ 0) (ha : ∀ u v, G.Adj u v → m u v = a u * a v)
    {u v : V} (hu : a u ≠ 0) (w : G.Walk u v) :
    a v = if Even w.length then a u / alternatingRatio m w else alternatingRatio m w / a u := by
  obtain ⟨he,ho⟩ := alternatingRatio_spec hm ha w
  split_ifs with h
  · exact (eq_div_iff (ratio_ne_zero hm w)).mpr (by rw [mul_comm, he h])
  · exact (eq_div_iff hu).mpr (by rw [ho (Nat.not_even_iff_odd.mp h), mul_comm])

theorem odd_walk_magnitude {m : V → V → ℝ} {a : V → ℝ}
    (hm : ∀ u v, G.Adj u v → m u v ≠ 0) (ha : ∀ u v, G.Adj u v → m u v = a u * a v)
    {root : V} (w : G.Walk root root) (hw : Odd w.length) :
    Real.sqrt (alternatingRatio m w) = |a root| := by
  rw [(alternatingRatio_spec hm ha w).2 hw, ← sq, Real.sqrt_sq_eq_abs]

/-- An independent repeat at one context supplies the squared contrast used here. -/
theorem identify_from_square (hG : G.Preconnected) {a b : V → ℝ} (ha : ∀ v, a v ≠ 0)
    (edges : ∀ u v, G.Adj u v → b u * b v = a u * a v) (root : V)
    (anchor : b root ^ 2 = a root ^ 2) : b = a ∨ b = -a := by
  have hn : |b root| = |a root| := by
    simpa only [Real.sqrt_sq_eq_abs] using congrArg Real.sqrt anchor
  obtain ⟨c, _, hc⟩ := TimeFrequency.exists_phase_of_vertex (G := G)
    (x := fun v ↦ (a v : ℂ)) (y := fun v ↦ (b v : ℂ)) hG
    (fun v ↦ by exact_mod_cast ha v) root (by simpa using hn)
    (fun u v h ↦ by simpa using congrArg (fun r : ℝ ↦ (r : ℂ)) (edges u v h))
  have eqn (v : V) : b v = c.re * a v := by
    simpa [Complex.mul_re] using congrArg Complex.re (congrFun hc v)
  rcases sq_eq_sq_iff_eq_or_eq_neg.mp anchor with h | h
  · have hc1 : c.re = 1 := mul_right_cancel₀ (ha root)
      (by rw [one_mul, ← eqn, h])
    exact Or.inl (funext fun v ↦ by simp [eqn v, hc1])
  · have hc1 : c.re = -1 := mul_right_cancel₀ (ha root)
      (by rw [neg_one_mul, ← eqn, h])
    exact Or.inr (funext fun v ↦ by simp [eqn v, hc1])

theorem identify_from_odd_walk (hG : G.Preconnected) {a b : V → ℝ} (ha : ∀ v, a v ≠ 0)
    (edges : ∀ u v, G.Adj u v → b u * b v = a u * a v)
    {root : V} (w : G.Walk root root) (hw : Odd w.length) : b = a ∨ b = -a := by
  let m := fun u v ↦ a u * a v
  have hm (u v : V) (_ : G.Adj u v) : m u v ≠ 0 := mul_ne_zero (ha u) (ha v)
  have h₁ := (alternatingRatio_spec hm (a := a) (fun _ _ _ ↦ rfl) w).2 hw
  have h₂ := (alternatingRatio_spec hm (a := b) (fun u v e ↦ (edges u v e).symm) w).2 hw
  exact identify_from_square hG ha edges root (by simpa [sq] using h₂.symm.trans h₁)

/-- Even without a magnitude anchor, relative labels align along nonzero-contrast paths. -/
theorem informative_alignment {a b : V → ℝ}
    (edges : ∀ u v, G.Adj u v → b u * b v = a u * a v)
    {u v : {i // a i ≠ 0}} (path : (G.induce {i | a i ≠ 0}).Reachable u v) :
    (0 < a u * b u ↔ 0 < a v * b v) := by
  obtain ⟨w⟩ := path
  induction w with
  | nil => rfl
  | @cons i j k h _ ih =>
    apply Iff.trans _ ih
    have he := edges i j h
    have hpos : 0 < (a i * b i) * (a j * b j) := by
      calc
        0 < (a i * a j)^2 := sq_pos_of_ne_zero (mul_ne_zero i.property j.property)
        _ = (a i * b i) * (a j * b j) := by linear_combination -(a i * a j) * he
    rcases mul_pos_iff.mp hpos with h | h
    · exact iff_of_true h.1 h.2
    · exact iff_of_false (not_lt.mpr h.1.le) (not_lt.mpr h.2.le)

/-- A label exchange on any union of informative components preserves all observed products. -/
theorem component_flip {a : V → ℝ} (flip : V → Bool)
    (hc : ∀ u v, G.Adj u v → a u ≠ 0 → a v ≠ 0 → flip u = flip v)
    {u v : V} (h : G.Adj u v) :
    (if flip u then -a u else a u) * (if flip v then -a v else a v) = a u * a v := by
  by_cases hu : a u = 0
  · simp [hu]
  by_cases hv : a v = 0
  · simp [hv]
  rw [hc u v h hu hv]
  cases flip v <;> simp

def bipartiteRescale (side : V → Bool) (c : ℝ) (a : V → ℝ) (v : V) :=
  if side v then c * a v else a v / c

theorem bipartite_products (side : V → Bool) {c : ℝ} (hc : 0 < c) (a : V → ℝ)
    (color : ∀ u v, G.Adj u v → side u ≠ side v) {u v : V} (h : G.Adj u v) :
    bipartiteRescale side c a u * bipartiteRescale side c a v = a u * a v := by
  have hne := color u v h
  cases hu : side u <;> cases hv : side v <;> simp [hu, hv] at hne
  all_goals
    simp [bipartiteRescale, hu, hv]
    field_simp

#print axioms alternatingRatio_spec
#print axioms contrast_from_walk
#print axioms odd_walk_magnitude
#print axioms identify_from_square
#print axioms identify_from_odd_walk
#print axioms informative_alignment
#print axioms component_flip
#print axioms bipartite_products

end
end TimeFrequency.LinkedProfiles
