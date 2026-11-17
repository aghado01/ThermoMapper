import Enthymemata.TimeFrequency.DepthFourSharpness

/-! # Exact predictive cost of the stationary Fourier ambiguity -/

namespace TimeFrequency.ContextTree.DepthFour

open Finset

noncomputable def excess (t q : Word 4 → ℝ) : ℝ := ∑ h, (1 / 16 : ℝ) * (t h - q h) ^ 2
noncomputable def predict (a b : ℝ) (h : Word 4) := kernel a b h true

private theorem triple_orthogonal (r s : Fin 4) :
    ∑ h : Word 4, triple r h * triple s h = if r = s then 16 else 0 := by decide +revert

theorem prediction_distance (a b c d : ℝ) :
    excess (predict a b) (predict c d) =
      ((a - c) ^ 2 + (b - d) ^ 2 + ((a + b) - (c + d)) ^ 2) / 4 := by
  have orth (r s : Fin 4) :
      (∑ h : Word 4, (triple r h : ℝ) * (triple s h : ℝ)) = if r = s then 16 else 0 := by
    exact_mod_cast triple_orthogonal r s
  have square (h : Word 4) :
      (1 / 16 : ℝ) * (predict a b h - predict c d h) ^ 2 =
      ((a - c) ^ 2 * (triple 3 h : ℝ) * (triple 3 h : ℝ) +
        (b - d) ^ 2 * (triple 2 h : ℝ) * (triple 2 h : ℝ) +
        ((a + b) - (c + d)) ^ 2 * (triple 1 h : ℝ) * (triple 1 h : ℝ) +
        2 * (a - c) * (b - d) * (triple 3 h : ℝ) * (triple 2 h : ℝ) -
        2 * (a - c) * ((a + b) - (c + d)) * (triple 3 h : ℝ) * (triple 1 h : ℝ) -
        2 * (b - d) * ((a + b) - (c + d)) * (triple 2 h : ℝ) * (triple 1 h : ℝ)) / 64 := by
    simp only [predict, prediction_formula]
    ring
  simp only [excess, square, ← sum_div, sum_sub_distrib, sum_add_distrib, mul_assoc, ← mul_sum]
  rw [orth, orth, orth, orth, orth, orth]
  norm_num [Fin.ext_iff]
  ring

theorem midpoint_risk (q0 q1 : Word 4 → ℝ) :
    excess (fun h ↦ (q0 h + q1 h) / 2) q0 = excess q0 q1 /4 ∧
    excess (fun h ↦ (q0 h + q1 h) / 2) q1 = excess q0 q1 /4 := by
  constructor <;> simp only [excess, sum_div] <;> apply sum_congr rfl <;> intros <;> ring

theorem unavoidable_excess (t q0 q1 : Word 4 → ℝ) :
    excess q0 q1 /4 ≤ max (excess t q0) (excess t q1) := by
  have identity :
      excess t q0 + excess t q1 = 2 * excess t (fun h ↦ (q0 h + q1 h) / 2) + excess q0 q1 / 2 := by
    simp only [excess, mul_sum, sum_div, ← sum_add_distrib]
    apply sum_congr rfl; intros; ring
  have hp : 0 ≤ excess t (fun h ↦ (q0 h + q1 h) / 2) := sum_nonneg (by intro h _; positivity)
  have h0 := le_max_left (excess t q0) (excess t q1)
  have h1 := le_max_right (excess t q0) (excess t q1)
  linarith

theorem stationary_minimax (a b c d : ℝ) (t : Word 4 → ℝ) :
    ((a - c) ^ 2 + (b - d) ^ 2 + ((a + b) - (c + d)) ^ 2) / 16 ≤
      max (excess t (predict a b)) (excess t (predict c d)) ∧
    let mid := fun h ↦ (predict a b h + predict c d h) / 2
    excess mid (predict a b) = ((a - c) ^ 2 + (b - d) ^ 2 + ((a + b) - (c + d)) ^ 2) / 16 ∧
    excess mid (predict c d) = ((a - c) ^ 2 + (b - d) ^ 2 + ((a + b) - (c + d)) ^ 2) / 16 := by
  have h := unavoidable_excess t (predict a b) (predict c d)
  have hm := midpoint_risk (predict a b) (predict c d)
  rw [prediction_distance] at h hm
  have he (x : ℝ) : x / 4 / 4 = x / 16 := by ring
  simpa only [he] using And.intro h hm

theorem one_taper_predictive_cost (r : Fin 5) (ε : ℝ) (t : Word 4 → ℝ) :
    ε ^ 2 / 8 ≤ max (excess t (predict 0 0))
      (excess t (predict (ε * invisibleA r) (ε * invisibleB r))) ∧
    let mid := fun h ↦ (predict 0 0 h + predict (ε * invisibleA r) (ε * invisibleB r) h) / 2
    excess mid (predict 0 0) = ε ^ 2 / 8 ∧
    excess mid (predict (ε * invisibleA r) (ε * invisibleB r)) = ε ^ 2 / 8 := by
  have he : ((0-ε * invisibleA r) ^ 2 + (0-ε * invisibleB r) ^ 2+
      ((0 + 0) - (ε * invisibleA r + ε * invisibleB r)) ^ 2) / 16 = ε ^ 2 / 8 := by
    fin_cases r <;> simp [invisibleA, invisibleB] <;> ring
  simpa only [he] using stationary_minimax 0 0 (ε * invisibleA r) (ε * invisibleB r) t

theorem symmetric_pair_cost (r : Fin 5) (ε : ℝ) (t : Word 4 → ℝ) :
    ε ^ 2 / 2 ≤ max (excess t (predict (ε * invisibleA r) (ε * invisibleB r)))
      (excess t (predict (-ε * invisibleA r) (-ε * invisibleB r))) ∧
    excess (predict 0 0) (predict (ε * invisibleA r) (ε * invisibleB r)) = ε ^ 2 / 2 ∧
    excess (predict 0 0) (predict (-ε * invisibleA r) (-ε * invisibleB r)) = ε ^ 2 / 2 := by
  have he : ((ε * invisibleA r - (-ε * invisibleA r)) ^ 2+
      (ε * invisibleB r - (-ε * invisibleB r)) ^ 2+
      ((ε * invisibleA r + ε * invisibleB r) -
        (-ε * invisibleA r + -ε * invisibleB r)) ^ 2) / 16 = ε ^ 2 / 2 := by
    fin_cases r <;> simp [invisibleA, invisibleB] <;> ring
  refine ⟨?_, ?_, ?_⟩
  · simpa only [he] using (stationary_minimax (ε * invisibleA r) (ε * invisibleB r)
      (-ε * invisibleA r) (-ε * invisibleB r) t).1
  · rw [prediction_distance]
    fin_cases r <;> simp [invisibleA, invisibleB] <;> ring
  · rw [prediction_distance]
    fin_cases r <;> simp [invisibleA, invisibleB] <;> ring

theorem sharp_eighth_lower (r : Fin 5) (t : Word 4 → ℝ) {c : ℝ} (hc : c < 1 / 8) :
    ∃ ε : ℝ, 0 < ε ∧ ε < 1 / 2 ∧
      c < max (excess t (predict (ε * invisibleA r) (ε * invisibleB r)))
        (excess t (predict (-ε * invisibleA r) (-ε * invisibleB r))) := by
  let δ := min (1 / 4 : ℝ) (1 / 8 - c)
  have hd : 0 < δ := lt_min (by norm_num) (by linarith)
  have hsmall : δ ≤ 1 / 4 := min_le_left _ _
  have hgap : δ ≤ 1 / 8 - c := min_le_right _ _
  refine ⟨1 / 2-δ, by linarith, by linarith, ?_⟩
  have he := (symmetric_pair_cost r (1 / 2-δ) t).1
  nlinarith [sq_nonneg δ]

/-- Actual indistinguishable positive depth-four laws with predictive cost arbitrarily
close to one eighth. Each source uses its uniform stationary history law. -/
theorem observed_excess_lower (r : Fin 5) (t : Word 4 → ℝ) {c : ℝ} (hc : c < 1 / 8) :
    ∃ a b : ℝ, Admissible a b ∧ Admissible (-a) (-b) ∧
      (∀ v, bin (weight a b) (taper r) v = bin (weight (-a) (-b)) (taper r) v) ∧
      (∀ w o (ho : o + w ≤ 5) (v : Fin w → ℝ),
        bin (weight a b) (fun x ↦ powerVector (slice x o w ho)) v =
          bin (weight (-a) (-b)) (fun x ↦ powerVector (slice x o w ho)) v) ∧
      c < max (excess t (predict a b)) (excess t (predict (-a) (-b))) := by
  obtain ⟨ε, hp, hu, hr⟩ := sharp_eighth_lower r t hc
  obtain ⟨ha, _, hz⟩ := one_taper_direction r hp hu
  refine ⟨ε * invisibleA r, ε * invisibleB r, ha, admissible_neg ha, ?_, ?_, ?_⟩
  · intro v
    rw [taper_law, taper_law, coefficient_neg, hz]
    ring
  · intro w o ho v
    rw [all_fourier_laws _ _ ho, all_fourier_laws _ _ ho]
  · simpa only [neg_mul] using hr

theorem fair_predictor_upper {a b : ℝ} (ha : Admissible a b) :
    excess (predict 0 0) (predict a b) < 1 / 8 := by
  let s := |a|+|b|+|a + b|
  have hs : s < 1 := ha
  have hs0 : 0 ≤ s := by dsimp [s]; positivity
  have ha0 : |a| ≤ |b|+|a + b| := by
    have h := abs_sub (a + b) b
    simpa [add_sub_cancel_right, add_comm] using h
  have hb0 : |b| ≤ |a|+|a + b| := by
    have h := abs_sub (a + b) a
    simpa [add_sub_cancel_left, add_comm] using h
  have hc0 : |a + b| ≤ |a|+|b| := abs_add_le _ _
  have hA : 2*|a| ≤ s := by dsimp [s]; linarith
  have hB : 2*|b| ≤ s := by dsimp [s]; linarith
  have hC : 2*|a + b| ≤ s := by dsimp [s]; linarith
  have h := add_le_add (add_le_add (mul_le_mul_of_nonneg_right hA (abs_nonneg a))
    (mul_le_mul_of_nonneg_right hB (abs_nonneg b)))
    (mul_le_mul_of_nonneg_right hC (abs_nonneg (a + b)))
  have he : s*|a|+s*|b|+s*|a + b| = s ^ 2 := by dsimp [s]; ring
  rw [he] at h
  have hsq (x : ℝ) : |x| * |x| = x ^ 2 := by nlinarith [sq_abs x]
  simp only [mul_assoc, hsq] at h
  rw [prediction_distance]
  nlinarith

theorem suffix_prediction (a b : ℝ) (h : Word 3) :
    (predict a b (Fin.cons false h) + predict a b (Fin.cons true h)) / 2 = 1 / 2 := by
  have ht (r : Fin 4) (hr : r ≠ 0) :
      triple r (Fin.cons false h) + triple r (Fin.cons true h) = 0 := by decide +revert
  have hA : (triple 3 (Fin.cons false h) : ℝ) + (triple 3 (Fin.cons true h) : ℝ) = 0 := by
    exact_mod_cast ht 3 (by decide)
  have hB : (triple 2 (Fin.cons false h) : ℝ) + (triple 2 (Fin.cons true h) : ℝ) = 0 := by
    exact_mod_cast ht 2 (by decide)
  have hC : (triple 1 (Fin.cons false h) : ℝ) + (triple 1 (Fin.cons true h) : ℝ) = 0 := by
    exact_mod_cast ht 1 (by decide)
  simp only [predict, prediction_formula]
  nlinarith [congrArg (fun z : ℝ ↦ a * z) hA, congrArg (fun z : ℝ ↦ b * z) hB,
    congrArg (fun z : ℝ ↦ (a + b) * z) hC]

noncomputable def functional (f q : Word 4 → ℝ) := ∑ h, (1 / 16 : ℝ) * f h * q h
noncomputable def contrast (f : Word 4 → ℝ) (r : Fin 4) :=
  ∑ h, (1 / 16 : ℝ) * f h * ((triple r h : ℝ) - (triple 1 h : ℝ))

theorem functional_formula (f : Word 4 → ℝ) (a b : ℝ) :
    functional f (predict a b) = functional f (predict 0 0) +
      (a * contrast f 3 + b * contrast f 2) / 2 := by
  simp only [functional, contrast, predict, prediction_formula, mul_sum, sum_div, ← sum_add_distrib]
  apply sum_congr rfl; intros; ring

theorem average_prediction (a b : ℝ) : functional (fun _ ↦ 1) (predict a b) = 1 / 2 := by
  have hz (r : Fin 4) : ∑ h : Word 4, triple r h = 0 := by decide +revert
  have hc (r : Fin 4) : contrast (fun _ ↦ 1) r = 0 := by
    have h0 : (∑ h : Word 4, (triple r h : ℝ)) = 0 := by exact_mod_cast hz r
    have h1 : (∑ h : Word 4, (triple 1 h : ℝ)) = 0 := by exact_mod_cast hz 1
    simp only [contrast, mul_one, ← mul_sum, sum_sub_distrib, h0, h1, sub_self, mul_zero]
  rw [functional_formula, hc, hc]
  norm_num [functional, predict, prediction_formula, Fintype.card_fun]

theorem identifiable_linear_functional (f : Word 4 → ℝ) :
    (∀ a b, Admissible a b → functional f (predict a b) = functional f (predict 0 0)) ↔
    contrast f 3 = 0 ∧ contrast f 2 = 0 := by
  constructor
  · intro h
    have hA := h (1 / 4) 0 (by norm_num [Admissible])
    have hB := h 0 (1 / 4) (by norm_num [Admissible])
    rw [functional_formula] at hA hB
    constructor <;> linarith
  · rintro ⟨ha, hb⟩ a b _
    rw [functional_formula, ha, hb]
    ring

theorem single_taper_functional (f : Word 4 → ℝ) (r : Fin 5) :
    (∀ ε : ℝ, 0 < ε → ε < 1 / 2 →
      functional f (predict (ε * invisibleA r) (ε * invisibleB r)) = functional f (predict 0 0)) ↔
    invisibleA r * contrast f 3 + invisibleB r * contrast f 2 = 0 := by
  constructor
  · intro h
    have he := h (1 / 4) (by norm_num) (by norm_num)
    rw [functional_formula] at he
    linarith
  · intro h ε _ _
    rw [functional_formula]
    nlinarith

def TargetVisible (f : Word 4 → ℝ) (r : Fin 5) : Prop :=
  if r = 3 then contrast f 2 = 0 else if r = 2 then contrast f 3 = 0 else
    if r = 1 then contrast f 3 = contrast f 2 else contrast f 3 = 0 ∧ contrast f 2 = 0

theorem one_taper_functionals (f : Word 4 → ℝ) (r : Fin 5) :
    (∀ a b c d, Admissible a b → Admissible c d → coefficient a b r = coefficient c d r →
      functional f (predict a b) = functional f (predict c d)) ↔ TargetVisible f r := by
  constructor
  · intro h
    have ht (a b : ℝ) (ha : Admissible a b) (he : coefficient a b r = 0) :
        functional f (predict a b) = functional f (predict 0 0) :=
      h a b 0 0 ha (by norm_num [Admissible]) (by simpa [coefficient] using he)
    fin_cases r
    · have ha := ht (1 / 4) 0 (by norm_num [Admissible]) (by norm_num [coefficient, Fin.ext_iff])
      have hb := ht 0 (1 / 4) (by norm_num [Admissible]) (by norm_num [coefficient, Fin.ext_iff])
      rw [functional_formula] at ha hb
      simp [TargetVisible, Fin.ext_iff]
      constructor <;> linarith
    · have he := ht (1 / 4) (-1 / 4) (by norm_num [Admissible])
        (by norm_num [coefficient, Fin.ext_iff])
      rw [functional_formula] at he
      simp [TargetVisible, Fin.ext_iff]
      linarith
    · have he := ht (1 / 4) 0 (by norm_num [Admissible]) (by norm_num [coefficient, Fin.ext_iff])
      rw [functional_formula] at he
      simp [TargetVisible, Fin.ext_iff]
      linarith
    · have he := ht 0 (1 / 4) (by norm_num [Admissible]) (by norm_num [coefficient, Fin.ext_iff])
      rw [functional_formula] at he
      simp [TargetVisible]
      linarith
    · have ha := ht (1 / 4) 0 (by norm_num [Admissible]) (by norm_num [coefficient, Fin.ext_iff])
      have hb := ht 0 (1 / 4) (by norm_num [Admissible]) (by norm_num [coefficient, Fin.ext_iff])
      rw [functional_formula] at ha hb
      simp [TargetVisible, Fin.ext_iff]
      constructor <;> linarith
  · intro hv a b c d _ _ he
    rw [functional_formula f a b, functional_formula f c d]
    fin_cases r
    · change contrast f 3 = 0 ∧ contrast f 2 = 0 at hv
      simp [hv.1, hv.2]
    · change contrast f 3 = contrast f 2 at hv
      change -(a + b) = -(c + d) at he
      rw [hv]
      nlinarith [congrArg (fun z : ℝ ↦ z * contrast f 2) he]
    · change contrast f 3 = 0 at hv
      change b = d at he
      simp [hv, he]
    · change contrast f 2 = 0 at hv
      change a = c at he
      simp [hv, he]
    · change contrast f 3 = 0 ∧ contrast f 2 = 0 at hv
      simp [hv.1, hv.2]

theorem augmented_prediction {P Q : Source} (h : SameAugmented P Q) :
    excess (fun s ↦ nextKernel 4 (P.law 5) s true)
      (fun s ↦ nextKernel 4 (Q.law 5) s true) = 0 := by
  rw [(completed_kernel_tree h).1]
  simp [excess]

theorem prediction_range {a b : ℝ} (ha : Admissible a b) (h : Word 4) :
    predict a b h ∈ Set.Ioo 0 1 := by
  obtain ⟨hp, ht⟩ := stochastic ha h
  simp only [Fintype.sum_bool] at ht
  exact ⟨hp true, by dsimp [predict]; linarith [hp false]⟩

theorem midpoint_range {a b c d : ℝ} (ha : Admissible a b) (hc : Admissible c d) (h : Word 4) :
    (predict a b h + predict c d h) / 2 ∈ Set.Ioo 0 1 := by
  have h0 := prediction_range ha h
  have h1 := prediction_range hc h
  constructor <;> linarith [h0.1, h0.2, h1.1, h1.2]

/-- Brier excess for a Bernoulli next-positive indicator. -/
theorem brier_excess (t q : ℝ) :
    q * (t - 1) ^ 2 + (1 - q) * t ^ 2 - (q * (q - 1) ^ 2 + (1 - q) * q ^ 2) = (t - q) ^ 2 := by ring

#print axioms prediction_range
#print axioms midpoint_range
#print axioms brier_excess
#print axioms prediction_distance
#print axioms midpoint_risk
#print axioms unavoidable_excess
#print axioms stationary_minimax
#print axioms one_taper_predictive_cost
#print axioms symmetric_pair_cost
#print axioms sharp_eighth_lower
#print axioms observed_excess_lower
#print axioms fair_predictor_upper
#print axioms suffix_prediction
#print axioms functional_formula
#print axioms average_prediction
#print axioms identifiable_linear_functional
#print axioms single_taper_functional
#print axioms one_taper_functionals
#print axioms augmented_prediction

end TimeFrequency.ContextTree.DepthFour
