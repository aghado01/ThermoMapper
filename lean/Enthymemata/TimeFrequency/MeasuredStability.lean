import Enthymemata.TimeFrequency.WeightedSignCut
import Enthymemata.TimeFrequency.FinitePowerMeasurements

/-! # Stability in the retained finite power and log-power vectors -/

namespace TimeFrequency.MeasuredStability

open SignCut FinitePower Finset

variable {T J : Type*} [Fintype T] [Fintype J]

theorem l2_scale (c : ℝ) (z : J → ℝ) : l2 (c • z) = |c| * l2 z := by
  change ‖c • (WithLp.toLp 2 z : EuclideanSpace ℝ J)‖ = _
  rw [norm_smul, Real.norm_eq_abs]
  rfl

theorem cut_data_sq (a : J → T → ℂ) (u v : T → ℝ) :
    l2 (powers a (u+v)-powers a (u-v))^2 =
      16 * ∑ j, ((amplitude (a j) u) * (starRingEnd ℂ) (amplitude (a j) v)).re^2 := by
  rw [l2_sq, mul_sum]
  apply sum_congr rfl
  intro j _
  rw [Pi.sub_apply, cut_difference]
  ring

theorem cut_distance_sq (u v : T → ℝ) :
    distance (u+v) (u-v)^2 = 4 * min (l2 u ^ 2) (l2 v ^ 2) := by
  have hm : u+v-(u-v) = (2 : ℝ) • v := by ext i; simp; ring
  have hp : u+v+(u-v) = (2 : ℝ) • u := by ext i; simp; ring
  rw [distance, hm, hp, l2_scale, l2_scale]
  have hu : 0 ≤ l2 u := norm_nonneg _
  have hv : 0 ≤ l2 v := norm_nonneg _
  rcases le_total (l2 u) (l2 v) with h | h
  · rw [min_eq_right (mul_le_mul_of_nonneg_left h (abs_nonneg _)),
      min_eq_left (pow_le_pow_left₀ hu h 2)]
    norm_num
    ring
  · rw [min_eq_left (mul_le_mul_of_nonneg_left h (abs_nonneg _)),
      min_eq_right (pow_le_pow_left₀ hv h 2)]
    norm_num
    ring

/-- Ratios use the sign-quotient waveform distance. A positive denominator excludes
equivalent pairs. Nonemptiness is explicit wherever an infimum is used. -/
def ratios (M : (T → ℝ) → J → ℝ) (C : Set (T → ℝ)) : Set ℝ :=
  {b | ∃ x ∈ C, ∃ y ∈ C, 0 < distance x y ∧ b = l2 (M x - M y) / distance x y}

noncomputable def lowerConstant (M : (T → ℝ) → J → ℝ) (C : Set (T → ℝ)) := sInf (ratios M C)

theorem lowerConstant_spec (M : (T → ℝ) → J → ℝ) (C : Set (T → ℝ)) (hn : (ratios M C).Nonempty)
    (k : ℝ) : k ≤ lowerConstant M C ↔
      ∀ x ∈ C, ∀ y ∈ C, k * distance x y ≤ l2 (M x - M y) := by
  have hb : BddBelow (ratios M C) := ⟨0, by
    rintro _ ⟨x,_,y,_,hd,rfl⟩
    exact div_nonneg (norm_nonneg _) hd.le⟩
  constructor
  · intro hk x hx y hy
    by_cases hd : 0 < distance x y
    · exact (le_div_iff₀ hd).mp (hk.trans (csInf_le hb ⟨x,hx,y,hy,hd,rfl⟩))
    · have h0 : distance x y = 0 := le_antisymm (le_of_not_gt hd)
        (le_min (norm_nonneg _) (norm_nonneg _))
      rw [h0,mul_zero]
      exact norm_nonneg _
  · intro h
    apply le_csInf hn
    rintro _ ⟨x,hx,y,hy,hd,rfl⟩
    exact (le_div_iff₀ hd).mpr (h x hx y hy)

theorem cut_upper_bound (a : J → T → ℂ) (C : Set (T → ℝ)) (u v : T → ℝ)
    (hp : u + v ∈ C) (hm : u - v ∈ C) (hd : 0 < distance (u + v) (u - v)) :
    lowerConstant (powers a) C ≤ l2 (powers a (u+v)-powers a (u-v)) / distance (u+v) (u-v) := by
  apply csInf_le
  · refine ⟨0, ?_⟩
    rintro _ ⟨_,_,_,_,hq,rfl⟩
    exact div_nonneg (norm_nonneg _) hq.le
  · exact ⟨_,hp,_,hm,hd,rfl⟩

/-- Transfer from a fixed linear extraction of the included data, in the stated Euclidean norms. -/
theorem extraction_bound {I : Type*} [Fintype I]
    (B : EuclideanSpace ℝ J →L[ℝ] EuclideanSpace ℝ I) (p p' : J → ℝ) (q q' : I → ℝ)
    (hp : B (WithLp.toLp 2 p) = WithLp.toLp 2 q) (hp' : B (WithLp.toLp 2 p') = WithLp.toLp 2 q') :
    l2 (q-q') ≤ ‖B‖ * l2 (p-p') := by
  have h := B.le_opNorm ((WithLp.toLp 2 p : EuclideanSpace ℝ J) - WithLp.toLp 2 p')
  rwa [map_sub,hp,hp'] at h

theorem fixedAmplitude_transfer [LinearOrder T] (G : SimpleGraph T) (x y : T → ℝ)
    (hs : (properCuts x).Nonempty) (hy : ∀ t, |y t| = |x t|) (M : (T → ℝ) → J → ℝ) {b : ℝ}
    (hex : l2 (products G x - products G y) ≤ b * l2 (M x - M y)) :
    cutConstant G x hs * distance x y ^ 2 ≤ b^2 * l2 (M x - M y)^2 := by
  have h := (sharp_constant_iff G x hs _).mpr le_rfl y hy
  have hn : 0 ≤ l2 (products G x - products G y) := norm_nonneg _
  exact h.trans (by simpa only [mul_pow] using pow_le_pow_left₀ hn hex 2)

/-- DC and Nyquist extraction on each edge. The operator sends (p₀,p₁) to (p₀−p₁)/4. -/
noncomputable def pairExtract {I : Type*} (p : I × Fin 2 → ℝ) (i : I) :=
  (p (i,0) - p (i,1))/4

theorem pairExtract_exact {I : Type*} (u v : I → ℝ) :
    pairExtract (fun ik : I × Fin 2 ↦ framePower ![u ik.1, v ik.1] ik.2.val) =
      fun i ↦ u i*v i := by
  funext i
  simpa [pairExtract, twoSample_power_difference] using
    (show 4 * u i * v i / 4 = u i*v i by ring)

theorem pairExtract_bound {I : Type*} [Fintype I] (p q : I × Fin 2 → ℝ) :
    8*l2 (pairExtract p - pairExtract q)^2 ≤ l2 (p-q)^2 := by
  rw [l2_sq, l2_sq, Fintype.sum_prod_type, mul_sum]
  apply sum_le_sum
  intro i _
  simp only [Pi.sub_apply, pairExtract, Fin.sum_univ_two]
  nlinarith [sq_nonneg ((p (i,0)-q (i,0)) + (p (i,1)-q (i,1)))]

/-- A single retained DC power on three unit coefficients cannot extract adjacent products. -/
theorem dc_extraction_obstruction :
    framePower (![1,1,-1] : Fin 3 → ℝ) 0 = framePower ![1,-1,1] 0 ∧
      (1:ℝ)*1 ≠ 1*(-1) := by
  rw [framePower_dc, framePower_dc]
  norm_num [Fin.sum_univ_three]
  change (2 + (-1 : ℝ)) ^ 2 = (1 : ℝ) ^ 2
  norm_num

omit [Fintype J] in
theorem powers_scale (a : J → T → ℂ) (s : ℝ) (x : T → ℝ) :
    powers a (s • x) = s^2 • powers a x := by
  funext j
  have h : amplitude (a j) (s • x) = (s : ℂ) * amplitude (a j) x := by
    simp only [amplitude, Pi.smul_apply, smul_eq_mul, Complex.ofReal_mul, mul_sum]
    exact sum_congr rfl fun i _ ↦ by ring
  simp [powers,h,Complex.normSq_ofReal, sq]

theorem no_positive_bound_on_ray (a : J → T → ℂ) (C : Set (T → ℝ)) (x : T → ℝ)
    (hx : 0 < l2 x) (h0 : 0 ∈ C)
    (hr : ∀ δ : ℝ, 0 < δ → ∃ s : ℝ, 0 < s ∧ s < δ ∧ s • x ∈ C) :
    ¬ ∃ k : ℝ, 0 < k ∧ ∀ u ∈ C, ∀ v ∈ C, k * distance u v ≤ l2 (powers a u-powers a v) := by
  rintro ⟨k,hk,h⟩
  have hp : 0 ≤ l2 (powers a x) := norm_nonneg _
  obtain ⟨s,hs,hl,hc⟩ := hr (k*l2 x/(l2 (powers a x)+1)) (by positivity)
  have hz : powers a 0 = 0 := by ext j; simp [powers, amplitude]
  have he := h _ hc _ h0
  rw [distance, sub_zero, add_zero, min_self, hz, sub_zero, powers_scale, l2_scale,
    l2_scale, abs_of_pos hs, abs_of_nonneg (sq_nonneg _)] at he
  have hh := (lt_div_iff₀ (by linarith : 0 < l2 (powers a x)+1)).mp hl
  nlinarith [mul_lt_mul_of_pos_left hh hs]

/-- A single real sample has a positive varying-amplitude control above a fixed floor. -/
theorem scalar_lower {x y a : ℝ} (hx : a ≤ |x|) (hy : a ≤ |y|) :
    2*a*min |x-y| |x+y| ≤ |x^2-y^2| := by
  have hl : 2*a ≤ max |x-y| |x+y| := by
    rcases le_total 0 x with h | h <;> rcases le_total 0 y with g | g
    · apply le_trans _ (le_max_right _ _)
      rw [abs_of_nonneg h] at hx
      rw [abs_of_nonneg g] at hy
      linarith [le_abs_self (x+y)]
    · apply le_trans _ (le_max_left _ _)
      rw [abs_of_nonneg h] at hx
      rw [abs_of_nonpos g] at hy
      linarith [le_abs_self (x-y)]
    · apply le_trans _ (le_max_left _ _)
      rw [abs_of_nonpos h] at hx
      rw [abs_of_nonneg g] at hy
      linarith [neg_le_abs (x-y)]
    · apply le_trans _ (le_max_right _ _)
      rw [abs_of_nonpos h] at hx
      rw [abs_of_nonpos g] at hy
      linarith [neg_le_abs (x+y)]
  refine (mul_le_mul_of_nonneg_right hl (le_min (abs_nonneg _) (abs_nonneg _))).trans_eq ?_
  rw [show x^2-y^2=(x-y)*(x+y) by ring, abs_mul]
  rcases le_total |x-y| |x+y| with h | h <;> simp [h,mul_comm]

noncomputable def edgePowers [LinearOrder T] (G : SimpleGraph T) (x : T → ℝ)
    (e : (T × T) × Fin 2) : ℝ := open Classical in
  if e.1.1 < e.1.2 ∧ G.Adj e.1.1 e.1.2 then framePower ![x e.1.1,x e.1.2] e.2.val else 0

omit [Fintype T] in
theorem edge_extraction [LinearOrder T] (G : SimpleGraph T) (x : T → ℝ) :
    pairExtract (edgePowers G x) = products G x := by
  funext e
  by_cases he : e.1 < e.2 ∧ G.Adj e.1 e.2
  · simp only [pairExtract, edgePowers, products, if_pos he]
    change (framePower _ 0 - framePower _ 1)/4 = _
    rw [twoSample_power_difference]
    ring
  · simp [pairExtract,edgePowers,products,he]

theorem fixedAmplitude_edgePowers [LinearOrder T] (G : SimpleGraph T) (x y : T → ℝ)
    (hs : (properCuts x).Nonempty) (hy : ∀ t, |y t| = |x t|) :
    8 * cutConstant G x hs * distance x y ^ 2 ≤ l2 (edgePowers G x-edgePowers G y)^2 := by
  have h := pairExtract_bound (edgePowers G x) (edgePowers G y)
  rw [edge_extraction,edge_extraction] at h
  have hg := (sharp_constant_iff G x hs _).mpr le_rfl y hy
  linarith

private theorem log_ordered {a b : ℝ} (ha : 0 < a) (hab : a ≤ b) :
    (b-a)/b ≤ Real.log b-Real.log a ∧ Real.log b-Real.log a ≤ (b-a)/a := by
  have hb := lt_of_lt_of_le ha hab
  rw [← Real.log_div hb.ne' ha.ne']
  have hr := div_pos hb ha
  constructor
  · calc (b-a)/b = 1-(b/a)⁻¹ := by field_simp
      _ ≤ _ := Real.one_sub_inv_le_log_of_pos hr
  · calc _ ≤ b/a-1 := Real.log_le_sub_one_of_pos hr
      _ = _ := by field_simp

theorem log_cell_bounds {p q ε U : ℝ} (hε : 0 < ε) (hp : p ∈ Set.Icc 0 U)
    (hq : q ∈ Set.Icc 0 U) :
    |p-q|/(U+ε) ≤ |Real.log (p+ε)-Real.log (q+ε)| ∧
    |Real.log (p+ε)-Real.log (q+ε)| ≤ |p-q|/ε := by
  have hord (p q : ℝ) (hP : p ∈ Set.Icc 0 U) (hQ : q ∈ Set.Icc 0 U) (h : q ≤ p) :
      |p-q|/(U+ε) ≤ |Real.log (p+ε)-Real.log (q+ε)| ∧
      |Real.log (p+ε)-Real.log (q+ε)| ≤ |p-q|/ε := by
    have hq' : 0 < q+ε := add_pos_of_nonneg_of_pos hQ.1 hε
    have hp' : 0 < p+ε := add_pos_of_nonneg_of_pos hP.1 hε
    have hh := log_ordered hq' (show q+ε ≤ p+ε by linarith)
    rw [show p+ε-(q+ε) = p-q by ring] at hh
    rw [abs_of_nonneg (sub_nonneg.mpr h), abs_of_nonneg
      (sub_nonneg.mpr (Real.log_le_log hq' (show q+ε ≤ p+ε by linarith)))]
    exact ⟨(div_le_div_of_nonneg_left (sub_nonneg.mpr h) hp'
        (show p+ε ≤ U+ε by linarith [hP.2])).trans hh.1,
      hh.2.trans (div_le_div_of_nonneg_left (sub_nonneg.mpr h) hε (by linarith [hQ.1]))⟩
  rcases le_total q p with h | h
  · exact hord p q hp hq h
  · simpa only [abs_sub_comm q p, abs_sub_comm (Real.log (q+ε)) (Real.log (p+ε))] using
      hord q p hq hp h

private theorem norm_domination (u v : J → ℝ) {b : ℝ} (hb : 0 ≤ b)
    (h : ∀ j, |u j| ≤ b * |v j|) : l2 u ≤ b * l2 v := by
  have hu : 0 ≤ l2 u := norm_nonneg _
  have hv : 0 ≤ l2 v := norm_nonneg _
  apply (sq_le_sq₀ hu (mul_nonneg hb hv)).mp
  rw [mul_pow,l2_sq,l2_sq,mul_sum]
  apply sum_le_sum
  intro j _
  simpa only [mul_pow, sq_abs] using pow_le_pow_left₀ (abs_nonneg _) (h j) 2

theorem log_vector_bounds (p q : J → ℝ) {ε U : ℝ} (hε : 0 < ε) (hU : 0 ≤ U)
    (hp : ∀ j, p j ∈ Set.Icc 0 U) (hq : ∀ j, q j ∈ Set.Icc 0 U) :
    let z := fun j ↦ Real.log (p j+ε)-Real.log (q j+ε)
    l2 (p-q)/(U+ε) ≤ l2 z ∧ l2 z ≤ l2 (p-q)/ε := by
  dsimp only
  have hu : 0 < U+ε := by positivity
  constructor
  · apply (div_le_iff₀ hu).mpr
    rw [mul_comm]
    apply norm_domination _ _ hu.le
    intro j
    simpa only [Pi.sub_apply, mul_comm (U+ε)] using
      (div_le_iff₀ hu).mp (log_cell_bounds hε (hp j) (hq j)).1
  · rw [div_eq_inv_mul]
    apply norm_domination _ _ (inv_nonneg.mpr hε.le)
    intro j
    simpa only [Pi.sub_apply, div_eq_inv_mul] using (log_cell_bounds hε (hp j) (hq j)).2

theorem residual_error (M : (T → ℝ) → J → ℝ) (C : Set (T → ℝ)) {k : ℝ} (hk : 0 < k)
    (h : ∀ x ∈ C, ∀ y ∈ C, k * distance x y ≤ l2 (M x - M y)) {x y : T → ℝ}
    (hx : x ∈ C) (hy : y ∈ C) (z : J → ℝ) {η ρ : ℝ}
    (hη : l2 (M x - z) ≤ η) (hρ : l2 (M y - z) ≤ ρ) : distance x y ≤ (η+ρ)/k := by
  apply (le_div_iff₀ hk).mpr
  rw [mul_comm]
  refine (h x hx y hy).trans ?_
  have ht := dist_triangle (WithLp.toLp 2 (M x) : EuclideanSpace ℝ J) (WithLp.toLp 2 z)
    (WithLp.toLp 2 (M y))
  rw [dist_comm (WithLp.toLp 2 z), dist_eq_norm, dist_eq_norm, dist_eq_norm] at ht
  exact ht.trans (add_le_add hη hρ)

#print axioms l2_scale
#print axioms cut_data_sq
#print axioms cut_distance_sq
#print axioms lowerConstant_spec
#print axioms cut_upper_bound
#print axioms extraction_bound
#print axioms fixedAmplitude_transfer
#print axioms pairExtract_exact
#print axioms pairExtract_bound
#print axioms dc_extraction_obstruction
omit [Fintype J] in
/-- A Hann endpoint absent from every cell prevents recovery of its product with a neighbor. -/
theorem blind_edge [DecidableEq T] (a : J → T → ℂ) {s t : T} (hst : s ≠ t) (ht : ∀ j, a j t = 0) :
    let u := Pi.single s (1 : ℝ)
    let v := Pi.single t (1 : ℝ)
    powers a (u+v) = powers a (u-v) ∧
      (u+v) s * (u+v) t ≠ (u-v) s * (u-v) t := by
  dsimp only
  constructor
  · apply missed_piece
    intro j
    right
    intro i
    by_cases hi : i = t
    · exact Or.inl (hi ▸ ht j)
    · exact Or.inr (Pi.single_eq_of_ne hi _)
  · norm_num [Pi.single_apply, hst, Ne.symm hst]

#print axioms blind_edge
#print axioms powers_scale
#print axioms no_positive_bound_on_ray
#print axioms scalar_lower
#print axioms edge_extraction
#print axioms fixedAmplitude_edgePowers
#print axioms log_cell_bounds
#print axioms log_vector_bounds
#print axioms residual_error

end TimeFrequency.MeasuredStability
