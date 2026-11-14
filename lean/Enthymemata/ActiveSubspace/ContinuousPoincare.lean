import Enthymemata.ActiveSubspace.MatrixPoincare
import Enthymemata.ActiveSubspace.ProductVariance
import Mathlib.Analysis.Fourier.AddCircle
import Mathlib.MeasureTheory.Integral.DivergenceTheorem

/-!
# Active subspaces — the matrix check on continuous laws (AS-2(d))

The statements of `Protolemmata/BASS/active-subspace-lemmas.md`, AS-2(d).

* **(i) Wirtinger on an interval.** `fourierCoeffOn_of_hasDerivAt_off_countable`: for a periodic
  `f` with a derivative off a countable set (a hinge's kinks), `cₙ(f) = T/(2πin) · cₙ(f′)`.
  `wirtinger_periodic`: so, by Parseval on `f − c₀`, `T⁻¹∫|f − c₀|² ≤ (T/2π)² T⁻¹∫|f′|²`.
  `wirtinger_interval`: on `[a, b]` with no end condition, through the even extension
  (`integral_reflect`), `Var_U(h) ≤ ((b − a)/π)² E_U[h′²]` under the uniform law, kinks allowed.
* **(iii) The glue, scalar form.** `uniformIcc` with `integral_uniformIcc` and `varM_uniformIcc`;
  `integral_integral_update`; `poincare_pi`: under a product of uniforms,
  `Var F ≤ Σᵢ ((bᵢ − aᵢ)/π)² E[Dᵢ²]`.
* **(iii) The glue, matrix form.** `covM`, `energyM`, `covM_form`, `energyM_form`;
  `matrix_poincare_pi`: `Cov(f) ≼ E[J D Jᵀ]`, `D = diag(((bᵢ − aᵢ)/π)²)`, on the product of
  uniforms itself. `matrix_poincare_box`: the same with every hypothesis asked only on the box,
  through `clampBox` and `ae_mem_box`, so unbounded hinges and kinks that move with the other
  coordinates (a rotated hinge) are covered.
-/

open MeasureTheory Set Complex
open scoped Real

namespace ActiveSubspace

section Wirtinger

variable {a b : ℝ}

/-- The zeroth Fourier coefficient on `[a, b]` is the mean. -/
theorem fourierCoeffOn_zero (hab : a < b) (f : ℝ → ℂ) :
    fourierCoeffOn hab f 0 = (1 / (b - a) : ℝ) • ∫ x in a..b, f x := by
  rw [fourierCoeffOn_eq_integral]
  simp

/-- **AS-2(d)(i).** Fourier coefficients of a function with equal ends and a derivative off a
countable set: `cₙ(f) = T / (2πin) · cₙ(f′)`, `T = b − a`. Mathlib's form asks a one-sided
derivative everywhere; a hinge's kinks need the countable exception. -/
theorem fourierCoeffOn_of_hasDerivAt_off_countable (hab : a < b) {f f' : ℝ → ℂ} {n : ℤ}
    (hn : n ≠ 0) {s : Set ℝ} (hs : s.Countable) (hf : ContinuousOn f (Icc a b))
    (hd : ∀ x ∈ Ioo a b \ s, HasDerivAt f (f' x) x) (hf' : IntervalIntegrable f' volume a b)
    (hper : f a = f b) :
    fourierCoeffOn hab f n = ((b - a : ℝ) : ℂ) / (2 * π * I * n) * fourierCoeffOn hab f' n := by
  haveI hT : Fact (0 < b - a) := ⟨by linarith⟩
  set e : ℝ → ℂ := fun y => fourier (-n) (y : AddCircle (b - a)) with he
  set φ : ℝ → ℂ := fun y => ((b - a : ℝ) : ℂ) / (-2 * π * I * n) * e y with hφdef
  have hec : Continuous e := (map_continuous (fourier (-n))).comp (AddCircle.continuous_mk' _)
  have hφc : Continuous φ := continuous_const.mul hec
  have hφ : ∀ x, HasDerivAt φ (e x) x := fun x => has_antideriv_at_fourier_neg hT hn x
  have huIcc : uIcc a b = Icc a b := uIcc_of_le hab.le
  have hI1 : IntervalIntegrable (fun y => e y * f y) volume a b :=
    (hec.continuousOn.mul hf).intervalIntegrable_of_Icc hab.le
  have hI2 : IntervalIntegrable (fun y => φ y * f' y) volume a b :=
    hf'.continuousOn_mul (huIcc ▸ hφc.continuousOn)
  have hftc := integral_eq_of_hasDerivAt_off_countable_of_le (fun y => φ y * f y)
    (fun y => e y * f y + φ y * f' y) hab.le hs (hφc.continuousOn.mul hf)
    (fun x hx => (hφ x).mul (hd x hx)) (hI1.add hI2)
  have hab' : (b : AddCircle (b - a)) = (a : AddCircle (b - a)) := by
    simpa using AddCircle.coe_add_period (b - a) a
  have hφab : φ b = φ a := by simp only [hφdef, he, hab']
  rw [hper, hφab, sub_self, intervalIntegral.integral_add hI1 hI2] at hftc
  have hsplit : ∫ y in a..b, φ y * f' y =
      ((b - a : ℝ) : ℂ) / (-2 * π * I * n) * ∫ y in a..b, e y * f' y := by
    rw [← intervalIntegral.integral_const_mul]
    exact intervalIntegral.integral_congr fun y _ => by simp only [hφdef]; ring
  rw [hsplit] at hftc
  rw [fourierCoeffOn_eq_integral, fourierCoeffOn_eq_integral]
  simp only [smul_eq_mul, Complex.real_smul]
  have h1 : ∫ y in a..b, e y * f y = -(((b - a : ℝ) : ℂ) / (-2 * π * I * n) *
      ∫ y in a..b, e y * f' y) := eq_neg_of_add_eq_zero_left hftc
  simp only [he] at h1
  rw [h1]
  have hT0 : ((b - a : ℝ) : ℂ) ≠ 0 := Complex.ofReal_ne_zero.2 hT.out.ne'
  have hn0 : (n : ℂ) ≠ 0 := Int.cast_ne_zero.2 hn
  have hpi : (π : ℂ) ≠ 0 := Complex.ofReal_ne_zero.2 Real.pi_ne_zero
  field_simp

/-- A function continuous on `[a, b]` is in `L²` of `(a, b]`. -/
theorem memLp_two_of_continuousOn {g : ℝ → ℂ} (hg : ContinuousOn g (Icc a b)) :
    MemLp g 2 (volume.restrict (Ioc a b)) := by
  obtain ⟨C, hC⟩ := isCompact_Icc.exists_bound_of_continuousOn hg
  refine MemLp.of_bound ((hg.mono Ioc_subset_Icc_self).aestronglyMeasurable measurableSet_Ioc) C ?_
  exact (ae_restrict_iff' measurableSet_Ioc).2
    (Filter.Eventually.of_forall fun x hx => hC x (Ioc_subset_Icc_self hx))

/-- **AS-2(d)(i), periodic Wirtinger.** For `f` with equal ends and a derivative off a countable
set, square-integrable: `T⁻¹ ∫ |f − c₀|² ≤ (T / 2π)² T⁻¹ ∫ |f′|²`, `c₀` the mean. -/
theorem wirtinger_periodic (hab : a < b) {f f' : ℝ → ℂ} {s : Set ℝ} (hs : s.Countable)
    (hf : ContinuousOn f (Icc a b)) (hd : ∀ x ∈ Ioo a b \ s, HasDerivAt f (f' x) x)
    (hper : f a = f b) (hf'2 : MemLp f' 2 (volume.restrict (Ioc a b))) :
    (b - a)⁻¹ * ∫ x in a..b, ‖f x - fourierCoeffOn hab f 0‖ ^ 2 ≤
      ((b - a) / (2 * π)) ^ 2 * ((b - a)⁻¹ * ∫ x in a..b, ‖f' x‖ ^ 2) := by
  set c₀ := fourierCoeffOn hab f 0 with hc₀
  set g : ℝ → ℂ := fun x => f x - c₀ with hgdef
  have hgc : ContinuousOn g (Icc a b) := hf.sub continuousOn_const
  have hgd : ∀ x ∈ Ioo a b \ s, HasDerivAt g (f' x) x := fun x hx => by
    simpa [hgdef] using (hd x hx).sub_const c₀
  have hgper : g a = g b := by simp only [hgdef, hper]
  have hf'i : IntervalIntegrable f' volume a b :=
    (intervalIntegrable_iff_integrableOn_Ioc_of_le hab.le).2 (hf'2.integrable one_le_two)
  have hPg := hasSum_sq_fourierCoeffOn hab (memLp_two_of_continuousOn hgc)
  have hPf' := hasSum_sq_fourierCoeffOn hab hf'2
  have hT : 0 < b - a := by linarith
  have hterm : ∀ n : ℤ, ‖fourierCoeffOn hab g n‖ ^ 2 ≤
      ((b - a) / (2 * π)) ^ 2 * ‖fourierCoeffOn hab f' n‖ ^ 2 := by
    intro n
    by_cases hn : n = 0
    · subst hn
      have h0 : fourierCoeffOn hab g 0 = 0 := by
        rw [fourierCoeffOn_zero, intervalIntegral.integral_sub (hf.intervalIntegrable_of_Icc hab.le)
          intervalIntegrable_const, intervalIntegral.integral_const, smul_sub, smul_smul, hc₀,
          fourierCoeffOn_zero, one_div, inv_mul_cancel₀ hT.ne', one_smul, sub_self]
      rw [h0, norm_zero, zero_pow two_ne_zero]
      positivity
    · rw [fourierCoeffOn_of_hasDerivAt_off_countable hab hn hs hgc hgd hf'i hgper, norm_mul,
        mul_pow]
      gcongr
      have hn1 : (1 : ℝ) ≤ |(n : ℝ)| := by
        rw [← Int.cast_abs]
        exact_mod_cast Int.one_le_abs hn
      rw [norm_div, Complex.norm_real, Real.norm_of_nonneg hT.le, norm_mul, norm_mul, norm_mul,
        Complex.norm_I, Complex.norm_real, Real.norm_of_nonneg Real.pi_pos.le, Complex.norm_intCast,
        Complex.norm_ofNat, mul_one]
      apply div_le_div_of_nonneg_left hT.le (by positivity)
      nlinarith [Real.pi_pos]
  have hle := hasSum_le hterm hPg (hPf'.mul_left (((b - a) / (2 * π)) ^ 2))
  simpa only [smul_eq_mul, hgdef] using hle

/-- The even reflection about `a` doubles an integral: `∫_{2a−b}^{b} φ(a + |x − a|) = 2 ∫_a^b φ`. -/
theorem integral_reflect (hab : a < b) {φ : ℝ → ℝ} (hφ : IntervalIntegrable φ volume a b) :
    ∫ x in (2 * a - b)..b, φ (a + |x - a|) = 2 * ∫ x in a..b, φ x := by
  have hleft : ∀ x ∈ uIcc (2 * a - b) a, φ (a + |x - a|) = φ (2 * a - x) := by
    intro x hx
    rw [uIcc_of_le (by linarith)] at hx
    rw [abs_of_nonpos (by linarith [hx.2])]
    congr 1
    ring
  have hright : ∀ x ∈ uIcc a b, φ (a + |x - a|) = φ x := by
    intro x hx
    rw [uIcc_of_le hab.le] at hx
    rw [abs_of_nonneg (by linarith [hx.1])]
    congr 1
    ring
  have hrefl : IntervalIntegrable (fun x => φ (2 * a - x)) volume (2 * a - b) a := by
    have := hφ.comp_sub_left (2 * a)
    simpa [show 2 * a - a = a by ring, show 2 * a - (2 * a - b) = b by ring] using this.symm
  have hI1 : IntervalIntegrable (fun x => φ (a + |x - a|)) volume (2 * a - b) a :=
    (intervalIntegrable_congr fun x hx => (hleft x (uIoc_subset_uIcc hx)).symm).1 hrefl
  have hI2 : IntervalIntegrable (fun x => φ (a + |x - a|)) volume a b :=
    (intervalIntegrable_congr fun x hx => (hright x (uIoc_subset_uIcc hx)).symm).1 hφ
  rw [← intervalIntegral.integral_add_adjacent_intervals hI1 hI2,
    intervalIntegral.integral_congr hleft, intervalIntegral.integral_congr hright,
    intervalIntegral.integral_comp_sub_left (fun x => φ x) (2 * a)]
  rw [show 2 * a - a = a by ring, show 2 * a - (2 * a - b) = b by ring]
  ring

/-- **AS-2(d)(i), Wirtinger on an interval.** For `h` continuous on `[a, b]` with a bounded
measurable derivative off a countable set (a hinge's kinks allowed), its variance under the uniform
law is at most `((b − a)/π)²` times the mean squared derivative:
`T⁻¹ ∫ (h − m)² ≤ (T/π)² T⁻¹ ∫ h′²`, `T = b − a`, `m` the mean. The proof applies
`wirtinger_periodic` to the even extension `h(a + |x − a|)` on `[2a − b, b]`. -/
theorem wirtinger_interval (hab : a < b) {h h' : ℝ → ℝ} {s : Set ℝ} (hs : s.Countable)
    (hc : ContinuousOn h (Icc a b)) (hd : ∀ x ∈ Ioo a b \ s, HasDerivAt h (h' x) x)
    (hm : Measurable h') {K : ℝ} (hK : ∀ x, |h' x| ≤ K) :
    (b - a)⁻¹ * ∫ x in a..b, (h x - (b - a)⁻¹ * ∫ y in a..b, h y) ^ 2 ≤
      ((b - a) / π) ^ 2 * ((b - a)⁻¹ * ∫ x in a..b, h' x ^ 2) := by
  have hT : 0 < b - a := by linarith
  set a' := 2 * a - b with ha'
  have ha'b : a' < b := by rw [ha']; linarith
  set m := (b - a)⁻¹ * ∫ y in a..b, h y with hmdef
  set F : ℝ → ℂ := fun x => ((h (a + |x - a|) : ℝ) : ℂ) with hFdef
  set F' : ℝ → ℂ := fun x => ((if a < x then h' x else -h' (2 * a - x) : ℝ) : ℂ) with hF'def
  have hmaps : ∀ x ∈ Icc a' b, a + |x - a| ∈ Icc a b := by
    intro x hx
    rw [ha'] at hx
    constructor
    · linarith [abs_nonneg (x - a)]
    · rcases le_total a x with h1 | h1
      · rw [abs_of_nonneg (by linarith)]; linarith [hx.2]
      · rw [abs_of_nonpos (by linarith)]; linarith [hx.1]
  have hFc : ContinuousOn F (Icc a' b) := by
    refine Complex.continuous_ofReal.comp_continuousOn ?_
    exact hc.comp (continuous_const.add (continuous_abs.comp (continuous_id.sub
      continuous_const))).continuousOn hmaps
  have hFper : F a' = F b := by
    simp only [hFdef, ha']
    rw [show 2 * a - b - a = -(b - a) by ring, abs_neg, abs_of_pos hT,
      show a + (b - a) = b by ring]
  set s' : Set ℝ := s ∪ (fun y => 2 * a - y) '' s ∪ {a} with hs'def
  have hs' : s'.Countable := (hs.union (hs.image _)).union (countable_singleton a)
  have hFd : ∀ x ∈ Ioo a' b \ s', HasDerivAt F (F' x) x := by
    intro x hx
    obtain ⟨⟨hx1, hx2⟩, hxs⟩ := hx
    simp only [hs'def, mem_union, mem_singleton_iff, mem_image, not_or] at hxs
    obtain ⟨⟨hxs1, hxs2⟩, hxa⟩ := hxs
    rcases lt_or_gt_of_ne hxa with hlt | hgt
    · -- left of `a`: `F = h(2a − ·)` near `x`
      have hin : 2 * a - x ∈ Ioo a b \ s := by
        refine ⟨⟨by linarith, by rw [ha'] at hx1; linarith⟩, fun hm' => hxs2 ⟨2 * a - x, hm', ?_⟩⟩
        ring
      have hder : HasDerivAt (fun y => h (2 * a - y)) (h' (2 * a - x) * (-1)) x := by
        have h1 : HasDerivAt (fun y => 2 * a - y) (-1) x := by
          simpa using (hasDerivAt_id x).const_sub (2 * a)
        have := (hd _ hin).comp x h1
        simp only [Function.comp_def] at this
        exact this
      have hev : (fun y => ((h (a + |y - a|) : ℝ) : ℂ)) =ᶠ[nhds x]
          fun y => ((h (2 * a - y) : ℝ) : ℂ) := by
        filter_upwards [Iio_mem_nhds hlt] with y hy
        rw [abs_of_neg (by simpa using hy)]
        congr 2
        ring
      have := (hder.ofReal_comp).congr_of_eventuallyEq hev
      simpa [hF'def, not_lt.2 hlt.le, mul_neg, mul_one] using this
    · -- right of `a`: `F = h` near `x`
      have hin : x ∈ Ioo a b \ s := ⟨⟨hgt, hx2⟩, hxs1⟩
      have hev : (fun y => ((h (a + |y - a|) : ℝ) : ℂ)) =ᶠ[nhds x] fun y => ((h y : ℝ) : ℂ) := by
        filter_upwards [Ioi_mem_nhds hgt] with y hy
        rw [abs_of_pos (by simpa using hy)]
        congr 2
        ring
      have := ((hd x hin).ofReal_comp).congr_of_eventuallyEq hev
      simpa [hF'def, hgt] using this
  have hF'm : Measurable F' := by
    refine Complex.measurable_ofReal.comp (Measurable.ite measurableSet_Ioi hm ?_)
    exact (hm.comp (measurable_const.sub measurable_id)).neg
  have hF'2 : MemLp F' 2 (volume.restrict (Ioc a' b)) := by
    refine MemLp.of_bound hF'm.aestronglyMeasurable K (Filter.Eventually.of_forall fun x => ?_)
    simp only [hF'def, Complex.norm_real, Real.norm_eq_abs]
    split_ifs
    · exact hK x
    · rw [abs_neg]; exact hK _
  have hW := wirtinger_periodic ha'b hs' hFc hFd hFper hF'2
  -- the integrals over `[a', b]` are twice those over `[a, b]`
  have hhi : IntervalIntegrable h volume a b := hc.intervalIntegrable_of_Icc hab.le
  have hmean : fourierCoeffOn ha'b F 0 = (m : ℂ) := by
    rw [fourierCoeffOn_zero, hFdef, intervalIntegral.integral_ofReal,
      integral_reflect hab hhi, Complex.real_smul, ← Complex.ofReal_mul]
    congr 1
    rw [hmdef, show b - a' = 2 * (b - a) by rw [ha']; ring]
    field_simp
  have hvar : ∫ x in a'..b, ‖F x - fourierCoeffOn ha'b F 0‖ ^ 2 =
      2 * ∫ x in a..b, (h x - m) ^ 2 := by
    rw [hmean]
    simp only [hFdef, ← Complex.ofReal_sub, Complex.norm_real, Real.norm_eq_abs, sq_abs]
    exact integral_reflect hab
      (((hc.sub continuousOn_const).pow 2).intervalIntegrable_of_Icc hab.le)
  have hen : ∫ x in a'..b, ‖F' x‖ ^ 2 = 2 * ∫ x in a..b, h' x ^ 2 := by
    have hpt : ∀ x, ‖F' x‖ ^ 2 = h' (a + |x - a|) ^ 2 := by
      intro x
      simp only [hF'def, Complex.norm_real, Real.norm_eq_abs, sq_abs]
      split_ifs with hx
      · rw [abs_of_pos (by linarith)]
        congr 2
        ring
      · rw [abs_of_nonpos (by linarith), neg_sq]
        congr 2
        ring
    simp only [hpt]
    have hi : IntervalIntegrable (fun x => h' x ^ 2) volume a b := by
      refine (intervalIntegrable_const (c := K ^ 2)).mono_fun'
        (hm.pow_const 2).aestronglyMeasurable (Filter.Eventually.of_forall fun x => ?_)
      simp only [Real.norm_eq_abs, abs_pow]
      exact pow_le_pow_left₀ (abs_nonneg _) (hK x) 2
    exact integral_reflect hab hi
  rw [hvar, hen, show b - a' = 2 * (b - a) by rw [ha']; ring] at hW
  have h2 : (2 * (b - a))⁻¹ * (2 * ∫ x in a..b, (h x - m) ^ 2) =
      (b - a)⁻¹ * ∫ x in a..b, (h x - m) ^ 2 := by field_simp
  have h3 : (2 * (b - a))⁻¹ * (2 * ∫ x in a..b, h' x ^ 2) =
      (b - a)⁻¹ * ∫ x in a..b, h' x ^ 2 := by field_simp
  rw [h2, h3, show 2 * (b - a) / (2 * π) = (b - a) / π by field_simp] at hW
  exact hW

end Wirtinger

section Uniform

/-- `|∫ f| ≤ C` under a probability measure when `|f| ≤ C`. -/
theorem abs_integral_le_of_abs_le {X : Type*} [MeasurableSpace X] (μ : Measure X)
    [IsProbabilityMeasure μ] {f : X → ℝ} {C : ℝ} (hC : ∀ x, |f x| ≤ C) : |∫ x, f x ∂μ| ≤ C := by
  refine (abs_integral_le_integral_abs).trans ?_
  calc ∫ x, |f x| ∂μ ≤ ∫ _x, C ∂μ :=
        integral_mono_of_nonneg (ae_of_all _ fun _ => abs_nonneg _) (integrable_const C)
          (ae_of_all _ hC)
    _ = C := by simp

/-- The uniform probability measure on `[a, b]`. -/
noncomputable def uniformIcc (a b : ℝ) : Measure ℝ :=
  (ENNReal.ofReal (b - a))⁻¹ • volume.restrict (Icc a b)

theorem uniformIcc_isProbabilityMeasure {a b : ℝ} (hab : a < b) :
    IsProbabilityMeasure (uniformIcc a b) := by
  constructor
  rw [uniformIcc, Measure.smul_apply, Measure.restrict_apply MeasurableSet.univ, univ_inter,
    Real.volume_Icc, smul_eq_mul, ENNReal.inv_mul_cancel
      (ENNReal.ofReal_pos.2 (sub_pos.2 hab)).ne' ENNReal.ofReal_ne_top]

/-- Under the uniform law, an integral is the interval's mean. -/
theorem integral_uniformIcc {a b : ℝ} (hab : a < b) (φ : ℝ → ℝ) :
    ∫ y, φ y ∂(uniformIcc a b) = (b - a)⁻¹ * ∫ y in a..b, φ y := by
  rw [uniformIcc, integral_smul_measure, ENNReal.toReal_inv,
    ENNReal.toReal_ofReal (by linarith), integral_Icc_eq_integral_Ioc,
    ← intervalIntegral.integral_of_le hab.le, smul_eq_mul]

theorem varM_uniformIcc {a b : ℝ} (hab : a < b) (h : ℝ → ℝ) :
    varM (uniformIcc a b) h =
      (b - a)⁻¹ * ∫ t in a..b, (h t - (b - a)⁻¹ * ∫ y in a..b, h y) ^ 2 := by
  simp only [varM, integral_uniformIcc hab]

end Uniform

section Glue

variable {α : Type*} [MeasurableSpace α]

/-- Integrating out coordinate `i` keeps a product law's integral. -/
theorem integral_integral_update {n : ℕ} (μ : Fin n → Measure α)
    [∀ i, IsProbabilityMeasure (μ i)] (i : Fin n) {G : (Fin n → α) → ℝ} (hG : Measurable G)
    {C : ℝ} (hGC : ∀ x, |G x| ≤ C) :
    ∫ x, ∫ t, G (Function.update x i t) ∂(μ i) ∂(Measure.pi μ) = ∫ x, G x ∂(Measure.pi μ) := by
  obtain ⟨m, rfl⟩ := Nat.exists_eq_succ_of_ne_zero (Fin.pos i).ne'
  set e := MeasurableEquiv.piFinSuccAbove (fun _ : Fin (m + 1) => α) i
  have hmp : MeasurePreserving e (Measure.pi μ)
      ((μ i).prod (Measure.pi fun j => μ (i.succAbove j))) :=
    measurePreserving_piFinSuccAbove μ i
  have hsymm : ∀ s z, e.symm (s, z) = Fin.insertNth (α := fun _ => α) i s z := fun _ _ => rfl
  have hΦ : Measurable fun x => ∫ t, G (Function.update x i t) ∂(μ i) :=
    (hG.comp measurable_update').stronglyMeasurable.integral_prod_right'.measurable
  have hΦC : ∀ x, |∫ t, G (Function.update x i t) ∂(μ i)| ≤ C := fun x =>
    abs_integral_le_of_abs_le _ fun t => hGC _
  rw [← (hmp.symm e).integral_comp' G, integral_prod_symm (fun p => G (e.symm p))
      (integrable_of_bdd (f := fun p => G (e.symm p)) _ (hG.comp e.symm.measurable)
        (fun p => hGC _)),
    ← (hmp.symm e).integral_comp' (fun x => ∫ t, G (Function.update x i t) ∂(μ i)),
    integral_prod (fun p => ∫ t, G (Function.update (e.symm p) i t) ∂(μ i))
      (integrable_of_bdd (f := fun p => ∫ t, G (Function.update (e.symm p) i t) ∂(μ i)) _
        (hΦ.comp e.symm.measurable) (fun p => hΦC _))]
  simp only [hsymm, Fin.update_insertNth]
  simp

/-- **AS-2(d)(iii), the scalar form.** Under a product of uniform laws on `[aᵢ, bᵢ]`, a bounded
measurable `F` whose sections along each coordinate are continuous, with a bounded measurable
derivative `Dᵢ` off a countable set, satisfies
`Var F ≤ Σᵢ ((bᵢ − aᵢ)/π)² E[Dᵢ²]`: Efron–Stein (`efron_stein_pi`), Wirtinger on each section
(`wirtinger_interval`), and the coordinate integrated back out. -/
theorem poincare_pi {n : ℕ} (a b : Fin n → ℝ) (hab : ∀ i, a i < b i)
    {F : (Fin n → ℝ) → ℝ} (hF : Measurable F) {C : ℝ} (hC : ∀ x, |F x| ≤ C)
    (D : Fin n → (Fin n → ℝ) → ℝ) (hD : ∀ i, Measurable (D i)) {K : ℝ}
    (hK : ∀ i x, |D i x| ≤ K) (s : Fin n → (Fin n → ℝ) → Set ℝ) (hs : ∀ i x, (s i x).Countable)
    (hcont : ∀ i x, ContinuousOn (fun t => F (Function.update x i t)) (Icc (a i) (b i)))
    (hder : ∀ i x, ∀ t ∈ Ioo (a i) (b i) \ s i x,
      HasDerivAt (fun t => F (Function.update x i t)) (D i (Function.update x i t)) t) :
    varM (Measure.pi fun i => uniformIcc (a i) (b i)) F ≤
      ∑ i, ((b i - a i) / π) ^ 2 *
        ∫ x, D i x ^ 2 ∂(Measure.pi fun i => uniformIcc (a i) (b i)) := by
  haveI : ∀ i, IsProbabilityMeasure (uniformIcc (a i) (b i)) := fun i =>
    uniformIcc_isProbabilityMeasure (hab i)
  refine (efron_stein_pi n (fun i => uniformIcc (a i) (b i)) F hF C hC).trans
    (Finset.sum_le_sum fun i _ => ?_)
  have hpt : ∀ x, varM (uniformIcc (a i) (b i)) (fun t => F (Function.update x i t)) ≤
      ((b i - a i) / π) ^ 2 * ∫ t, D i (Function.update x i t) ^ 2 ∂(uniformIcc (a i) (b i)) := by
    intro x
    rw [varM_uniformIcc (hab i), integral_uniformIcc (hab i)]
    exact wirtinger_interval (hab i) (hs i x) (hcont i x) (hder i x)
      ((hD i).comp (measurable_update x)) (fun t => hK i _)
  have hD2 : Measurable fun x => D i x ^ 2 := (hD i).pow_const 2
  have hD2C : ∀ x, |D i x ^ 2| ≤ K ^ 2 := fun x => by
    rw [abs_pow]
    exact pow_le_pow_left₀ (abs_nonneg _) (hK i x) 2
  have hint1 : Integrable (fun x => varM (uniformIcc (a i) (b i))
      (fun t => F (Function.update x i t))) (Measure.pi fun i => uniformIcc (a i) (b i)) :=
    integrable_of_bdd _ (measurable_coordVar (fun i => uniformIcc (a i) (b i)) i hF)
      (fun x => abs_varM_le _ fun t => hC _)
  have hΨ : Measurable fun x => ∫ t, D i (Function.update x i t) ^ 2 ∂(uniformIcc (a i) (b i)) :=
    (hD2.comp measurable_update').stronglyMeasurable.integral_prod_right'.measurable
  have hint2 : Integrable (fun x => ((b i - a i) / π) ^ 2 *
      ∫ t, D i (Function.update x i t) ^ 2 ∂(uniformIcc (a i) (b i)))
      (Measure.pi fun i => uniformIcc (a i) (b i)) :=
    (integrable_of_bdd _ hΨ (fun x => abs_integral_le_of_abs_le _ fun t => hD2C _)).const_mul _
  calc coordVarPi (fun i => uniformIcc (a i) (b i)) i F
      ≤ ∫ x, ((b i - a i) / π) ^ 2 *
          ∫ t, D i (Function.update x i t) ^ 2 ∂(uniformIcc (a i) (b i))
            ∂(Measure.pi fun i => uniformIcc (a i) (b i)) :=
        integral_mono hint1 hint2 hpt
    _ = ((b i - a i) / π) ^ 2 * ∫ x, D i x ^ 2 ∂(Measure.pi fun i => uniformIcc (a i) (b i)) := by
        rw [MeasureTheory.integral_const_mul,
          integral_integral_update (fun i => uniformIcc (a i) (b i)) i hD2 hD2C]

end Glue

section MatrixForm

open Matrix

variable {o : Type*} [Fintype o]

/-- `Var g = E[g²] − (E g)²` under a probability measure. -/
theorem varM_eq_sub {X : Type*} [MeasurableSpace X] (P : Measure X) [IsProbabilityMeasure P]
    {g : X → ℝ} (hg : Integrable g P) (hg2 : Integrable (fun x => g x ^ 2) P) :
    varM P g = ∫ x, g x ^ 2 ∂P - (∫ x, g x ∂P) ^ 2 := by
  set c := ∫ x, g x ∂P with hc
  have hexp : (fun x => (g x - c) ^ 2) = fun x => g x ^ 2 - 2 * c * g x + c ^ 2 := by
    funext x
    ring
  have i1 : Integrable (fun x => g x ^ 2 - 2 * c * g x) P := hg2.sub (hg.const_mul (2 * c))
  simp only [varM]
  rw [← hc, hexp, integral_add i1 (integrable_const _), integral_sub hg2 (hg.const_mul (2 * c)),
    MeasureTheory.integral_const_mul, MeasureTheory.integral_const, ← hc]
  simp only [probReal_univ, one_smul]
  ring

/-- The output covariance matrix under a law. -/
noncomputable def covM {X : Type*} [MeasurableSpace X] (P : Measure X) (f : X → o → ℝ) :
    Matrix o o ℝ :=
  fun p q => ∫ x, f x p * f x q ∂P - (∫ x, f x p ∂P) * ∫ x, f x q ∂P

/-- `E[J D Jᵀ]` with `D = diag c`. -/
noncomputable def energyM {X : Type*} [MeasurableSpace X] {n : ℕ} (P : Measure X)
    (J : X → Matrix o (Fin n) ℝ) (c : Fin n → ℝ) : Matrix o o ℝ :=
  fun p q => ∑ i, c i * ∫ x, J x p i * J x q i ∂P

omit [Fintype o] in
theorem covM_isHermitian {X : Type*} [MeasurableSpace X] (P : Measure X) (f : X → o → ℝ) :
    (covM P f).IsHermitian := by
  ext p q
  simp only [conjTranspose_apply, star_trivial, covM, mul_comm]

omit [Fintype o] in
theorem energyM_isHermitian {X : Type*} [MeasurableSpace X] {n : ℕ} (P : Measure X)
    (J : X → Matrix o (Fin n) ℝ) (c : Fin n → ℝ) : (energyM P J c).IsHermitian := by
  ext p q
  simp only [conjTranspose_apply, star_trivial, energyM, mul_comm (J _ q _)]

/-- `E[(Σₚ vₚ hₚ)²] = Σₚ Σ_q vₚ v_q E[hₚ h_q]` for bounded measurable `hₚ`. -/
theorem integral_sq_lin {X : Type*} [MeasurableSpace X] (P : Measure X) [IsProbabilityMeasure P]
    {h : o → X → ℝ} (hh : ∀ p, Measurable (h p)) {C : ℝ} (hC : ∀ p x, |h p x| ≤ C)
    (v : o → ℝ) :
    ∫ x, (∑ p, v p * h p x) ^ 2 ∂P = ∑ p, ∑ q, v p * v q * ∫ x, h p x * h q x ∂P := by
  have hii : ∀ p q, Integrable (fun x => h p x * h q x) P := fun p q =>
    integrable_of_bdd _ ((hh p).mul (hh q)) (C := C * C) fun x => by
      rw [abs_mul]
      exact mul_le_mul (hC p x) (hC q x) (abs_nonneg _) ((abs_nonneg _).trans (hC p x))
  have hsq : (fun x => (∑ p, v p * h p x) ^ 2) =
      fun x => ∑ p, ∑ q, v p * v q * (h p x * h q x) := by
    funext x
    rw [sq, Finset.sum_mul_sum]
    exact Finset.sum_congr rfl fun p _ => Finset.sum_congr rfl fun q _ => by ring
  rw [hsq, integral_finsetSum _ fun p _ => integrable_finsetSum _ fun q _ => (hii p q).const_mul _]
  refine Finset.sum_congr rfl fun p _ => ?_
  rw [integral_finsetSum _ fun q _ => (hii p q).const_mul _]
  exact Finset.sum_congr rfl fun q _ => integral_const_mul _ _

/-- `vᵀ Cov(f) v = Var(vᵀf)` for bounded measurable outputs. -/
theorem covM_form {X : Type*} [MeasurableSpace X] (P : Measure X) [IsProbabilityMeasure P]
    {f : X → o → ℝ} (hf : ∀ p, Measurable fun x => f x p) {C : ℝ} (hC : ∀ x p, |f x p| ≤ C)
    (v : o → ℝ) : v ⬝ᵥ (covM P f *ᵥ v) = varM P (fun x => v ⬝ᵥ f x) := by
  have hi : ∀ p, Integrable (fun x => f x p) P := fun p =>
    integrable_of_bdd _ (hf p) (fun x => hC x p)
  have hg : Integrable (fun x => v ⬝ᵥ f x) P := by
    simp only [dotProduct]
    exact integrable_finsetSum _ fun p _ => (hi p).const_mul (v p)
  have hgm : ∫ x, v ⬝ᵥ f x ∂P = ∑ p, v p * ∫ x, f x p ∂P := by
    simp only [dotProduct]
    rw [integral_finsetSum _ fun p _ => (hi p).const_mul _]
    exact Finset.sum_congr rfl fun p _ => integral_const_mul _ _
  have hg2 : ∫ x, (v ⬝ᵥ f x) ^ 2 ∂P = ∑ p, ∑ q, v p * v q * ∫ x, f x p * f x q ∂P :=
    integral_sq_lin P (h := fun p x => f x p) hf (fun p x => hC x p) v
  have hg2i : Integrable (fun x => (v ⬝ᵥ f x) ^ 2) P := by
    have hsq : (fun x => (v ⬝ᵥ f x) ^ 2) =
        fun x => ∑ p, ∑ q, v p * v q * (f x p * f x q) := by
      funext x
      simp only [dotProduct, sq, Finset.sum_mul_sum]
      exact Finset.sum_congr rfl fun p _ => Finset.sum_congr rfl fun q _ => by ring
    rw [hsq]
    exact integrable_finsetSum _ fun p _ => integrable_finsetSum _ fun q _ =>
      (integrable_of_bdd _ ((hf p).mul (hf q)) (C := C * C) fun x => by
        change |f x p * f x q| ≤ C * C
        rw [abs_mul]
        exact mul_le_mul (hC x p) (hC x q) (abs_nonneg _) ((abs_nonneg _).trans (hC x p))
        ).const_mul _
  rw [quadForm_eq, varM_eq_sub P hg hg2i, hg2, hgm, sq, Finset.sum_mul_sum,
    ← Finset.sum_sub_distrib]
  refine Finset.sum_congr rfl fun p _ => ?_
  rw [← Finset.sum_sub_distrib]
  refine Finset.sum_congr rfl fun q _ => ?_
  simp only [covM]
  ring

/-- `vᵀ E[J D Jᵀ] v = Σᵢ cᵢ E[(Σₚ vₚ Jₚᵢ)²]`. -/
theorem energyM_form {X : Type*} [MeasurableSpace X] {n : ℕ} (P : Measure X)
    [IsProbabilityMeasure P] {J : X → Matrix o (Fin n) ℝ} (hJ : ∀ p i, Measurable fun x => J x p i)
    {K : ℝ} (hK : ∀ x p i, |J x p i| ≤ K) (c : Fin n → ℝ) (v : o → ℝ) :
    v ⬝ᵥ (energyM P J c *ᵥ v) = ∑ i, c i * ∫ x, (∑ p, v p * J x p i) ^ 2 ∂P := by
  have hR : ∑ i, c i * ∫ x, (∑ p, v p * J x p i) ^ 2 ∂P =
      ∑ i, ∑ p, ∑ q, v p * v q * (c i * ∫ x, J x p i * J x q i ∂P) := by
    refine Finset.sum_congr rfl fun i _ => ?_
    rw [integral_sq_lin P (h := fun p x => J x p i) (fun p => hJ p i) (fun p x => hK x p i) v,
      Finset.mul_sum]
    refine Finset.sum_congr rfl fun p _ => ?_
    rw [Finset.mul_sum]
    exact Finset.sum_congr rfl fun q _ => by ring
  have hL : v ⬝ᵥ (energyM P J c *ᵥ v) =
      ∑ p, ∑ q, ∑ i, v p * v q * (c i * ∫ x, J x p i * J x q i ∂P) := by
    rw [quadForm_eq]
    refine Finset.sum_congr rfl fun p _ => Finset.sum_congr rfl fun q _ => ?_
    simp only [energyM, Finset.mul_sum]
  have hswap : ∀ p, ∑ q, ∑ i, v p * v q * (c i * ∫ x, J x p i * J x q i ∂P) =
      ∑ i, ∑ q, v p * v q * (c i * ∫ x, J x p i * J x q i ∂P) := fun p => Finset.sum_comm
  rw [hL, hR]
  simp only [hswap]
  exact Finset.sum_comm

omit [Fintype o] in
/-- **AS-2(d)(iii), the matrix check on continuous laws.** Under a product of uniform laws on
`[aᵢ, bᵢ]`, bounded measurable outputs `f` whose sections along each coordinate are continuous with
derivative `J` (bounded, measurable) off a countable set satisfy
`Cov(f) ≼ E[J D Jᵀ]`, `D = diag(((bᵢ − aᵢ)/π)²)`: the matrix check of E8's slice 5, proved on the
law it reads, with no hypothesis left but the regularity of the fit. -/
theorem matrix_poincare_pi [Finite o] {n : ℕ} (a b : Fin n → ℝ) (hab : ∀ i, a i < b i)
    {f : (Fin n → ℝ) → o → ℝ} (hf : ∀ p, Measurable fun x => f x p) {C : ℝ}
    (hC : ∀ x p, |f x p| ≤ C) (J : (Fin n → ℝ) → Matrix o (Fin n) ℝ)
    (hJ : ∀ p i, Measurable fun x => J x p i) {K : ℝ} (hK : ∀ x p i, |J x p i| ≤ K)
    (s : Fin n → (Fin n → ℝ) → Set ℝ) (hs : ∀ i x, (s i x).Countable)
    (hcont : ∀ p i x, ContinuousOn (fun t => f (Function.update x i t) p) (Icc (a i) (b i)))
    (hder : ∀ p i x, ∀ t ∈ Ioo (a i) (b i) \ s i x,
      HasDerivAt (fun t => f (Function.update x i t) p) (J (Function.update x i t) p i) t) :
    (energyM (Measure.pi fun i => uniformIcc (a i) (b i)) J (fun i => ((b i - a i) / π) ^ 2) -
      covM (Measure.pi fun i => uniformIcc (a i) (b i)) f).PosSemidef := by
  haveI : Fintype o := Fintype.ofFinite o
  haveI : ∀ i, IsProbabilityMeasure (uniformIcc (a i) (b i)) := fun i =>
    uniformIcc_isProbabilityMeasure (hab i)
  rw [loewner_iff _ _ (covM_isHermitian _ f) (energyM_isHermitian _ J _)]
  intro v
  rw [covM_form _ hf hC, energyM_form _ hJ hK]
  have hbd : ∀ x, |v ⬝ᵥ f x| ≤ ∑ p, |v p| * C := fun x => by
    refine (Finset.abs_sum_le_sum_abs _ _).trans (Finset.sum_le_sum fun p _ => ?_)
    rw [abs_mul]
    exact mul_le_mul_of_nonneg_left (hC x p) (abs_nonneg _)
  have hDbd : ∀ i x, |∑ p, v p * J x p i| ≤ ∑ p, |v p| * K := fun i x => by
    refine (Finset.abs_sum_le_sum_abs _ _).trans (Finset.sum_le_sum fun p _ => ?_)
    rw [abs_mul]
    exact mul_le_mul_of_nonneg_left (hK x p i) (abs_nonneg _)
  exact poincare_pi a b hab (F := fun x => v ⬝ᵥ f x)
    (Finset.measurable_sum _ fun p _ => (hf p).const_mul (v p)) hbd
    (fun i x => ∑ p, v p * J x p i)
    (fun i => Finset.measurable_sum _ fun p _ => (hJ p i).const_mul (v p)) hDbd s hs
    (fun i x => continuousOn_finsetSum _ fun p _ => continuousOn_const.mul (hcont p i x))
    (fun i x t ht => HasDerivAt.fun_sum fun p _ => (hder p i x t ht).const_mul (v p))


section Box

variable {n : ℕ}

/-- Each coordinate clamped into its interval `[aᵢ, bᵢ]`. -/
def clampBox (a b : Fin n → ℝ) (x : Fin n → ℝ) : Fin n → ℝ :=
  fun i => max (a i) (min (b i) (x i))

theorem clampBox_mem {a b : Fin n → ℝ} (hab : ∀ i, a i < b i) (x : Fin n → ℝ) (i : Fin n) :
    clampBox a b x i ∈ Icc (a i) (b i) :=
  ⟨le_max_left _ _, max_le (hab i).le (min_le_left _ _)⟩

theorem clampBox_of_mem {a b : Fin n → ℝ} {x : Fin n → ℝ} (hx : ∀ i, x i ∈ Icc (a i) (b i)) :
    clampBox a b x = x := by
  funext i
  simp only [clampBox, min_eq_right (hx i).2, max_eq_right (hx i).1]

theorem measurable_clampBox (a b : Fin n → ℝ) : Measurable (clampBox a b) :=
  measurable_pi_lambda _ fun i =>
    measurable_const.max (measurable_const.min (measurable_pi_apply i))

theorem clampBox_update {a b : Fin n → ℝ} (x : Fin n → ℝ) (i : Fin n) {t : ℝ}
    (ht : t ∈ Icc (a i) (b i)) :
    clampBox a b (Function.update x i t) = Function.update (clampBox a b x) i t := by
  funext j
  by_cases hj : j = i
  · subst hj
    simp only [clampBox, Function.update_self, min_eq_right ht.2, max_eq_right ht.1]
  · simp only [clampBox, Function.update_of_ne hj]

/-- The box has full measure under the product of uniform laws. -/
theorem ae_mem_box {a b : Fin n → ℝ} (hab : ∀ i, a i < b i) :
    ∀ᵐ x ∂(Measure.pi fun i => uniformIcc (a i) (b i)), ∀ i, x i ∈ Icc (a i) (b i) := by
  haveI : ∀ i, IsProbabilityMeasure (uniformIcc (a i) (b i)) := fun i =>
    uniformIcc_isProbabilityMeasure (hab i)
  refine ae_all_iff.2 fun i => ?_
  have hi : ∀ᵐ t ∂(uniformIcc (a i) (b i)), t ∈ Icc (a i) (b i) :=
    Measure.ae_smul_measure (ae_restrict_mem measurableSet_Icc) _
  exact Measure.tendsto_eval_ae_ae.eventually hi

omit [Fintype o] in
/-- **AS-2(d)(iii), on the box.** The matrix check with every hypothesis asked only on the box
`Πᵢ [aᵢ, bᵢ]`: bounds, continuity of the sections, and their derivatives off a countable set that
may move with the other coordinates. A fit need not be bounded off the box (a hinge is not): it is
read through `clampBox`, which changes nothing under the law. -/
theorem matrix_poincare_box [Finite o] (a b : Fin n → ℝ) (hab : ∀ i, a i < b i)
    {f : (Fin n → ℝ) → o → ℝ} (hf : ∀ p, Measurable fun x => f x p) {C : ℝ}
    (hC : ∀ x, (∀ i, x i ∈ Icc (a i) (b i)) → ∀ p, |f x p| ≤ C)
    (J : (Fin n → ℝ) → Matrix o (Fin n) ℝ) (hJ : ∀ p i, Measurable fun x => J x p i) {K : ℝ}
    (hK : ∀ x, (∀ i, x i ∈ Icc (a i) (b i)) → ∀ p i, |J x p i| ≤ K)
    (s : Fin n → (Fin n → ℝ) → Set ℝ) (hs : ∀ i x, (s i x).Countable)
    (hcont : ∀ p i x, (∀ j, x j ∈ Icc (a j) (b j)) →
      ContinuousOn (fun t => f (Function.update x i t) p) (Icc (a i) (b i)))
    (hder : ∀ p i x, (∀ j, x j ∈ Icc (a j) (b j)) → ∀ t ∈ Ioo (a i) (b i) \ s i x,
      HasDerivAt (fun t => f (Function.update x i t) p) (J (Function.update x i t) p i) t) :
    (energyM (Measure.pi fun i => uniformIcc (a i) (b i)) J (fun i => ((b i - a i) / π) ^ 2) -
      covM (Measure.pi fun i => uniformIcc (a i) (b i)) f).PosSemidef := by
  set κ := clampBox a b with hκ
  have hmem : ∀ x i, κ x i ∈ Icc (a i) (b i) := fun x i => clampBox_mem hab x i
  have hupd : ∀ x i, ∀ t ∈ Icc (a i) (b i),
      κ (Function.update x i t) = Function.update (κ x) i t := fun x i t ht =>
    clampBox_update x i ht
  have hmain := matrix_poincare_pi a b hab (f := fun x => f (κ x))
    (fun p => (hf p).comp (measurable_clampBox a b)) (fun x p => hC _ (hmem x) p)
    (fun x => J (κ x)) (fun p i => (hJ p i).comp (measurable_clampBox a b))
    (fun x p i => hK _ (hmem x) p i) (fun i x => s i (κ x)) (fun i x => hs i (κ x))
    (fun p i x => (hcont p i (κ x) (hmem x)).congr fun t ht => by
      simp only [hupd x i t ht])
    (fun p i x t ht => by
      have hin : t ∈ Icc (a i) (b i) := Ioo_subset_Icc_self ht.1
      have hd := hder p i (κ x) (hmem x) t ht
      have hev : (fun t' => f (κ (Function.update x i t')) p) =ᶠ[nhds t]
          fun t' => f (Function.update (κ x) i t') p := by
        filter_upwards [Ioo_mem_nhds ht.1.1 ht.1.2] with t' ht'
        rw [hupd x i t' (Ioo_subset_Icc_self ht')]
      rw [hupd x i t hin]
      exact hd.congr_of_eventuallyEq hev)
  have hae := ae_mem_box hab
  have hfix : ∀ᵐ x ∂(Measure.pi fun i => uniformIcc (a i) (b i)), κ x = x :=
    hae.mono fun x hx => clampBox_of_mem hx
  have hcov : covM (Measure.pi fun i => uniformIcc (a i) (b i)) (fun x => f (κ x)) =
      covM (Measure.pi fun i => uniformIcc (a i) (b i)) f := by
    ext p q
    have e1 : ∫ x, f (κ x) p * f (κ x) q ∂(Measure.pi fun i => uniformIcc (a i) (b i)) =
        ∫ x, f x p * f x q ∂(Measure.pi fun i => uniformIcc (a i) (b i)) :=
      integral_congr_ae (hfix.mono fun x hx => by simp only [hx])
    have e2 : ∀ r : o, ∫ x, f (κ x) r ∂(Measure.pi fun i => uniformIcc (a i) (b i)) =
        ∫ x, f x r ∂(Measure.pi fun i => uniformIcc (a i) (b i)) := fun r =>
      integral_congr_ae (hfix.mono fun x hx => by simp only [hx])
    simp only [covM, e1, e2]
  have hen : energyM (Measure.pi fun i => uniformIcc (a i) (b i)) (fun x => J (κ x))
        (fun i => ((b i - a i) / π) ^ 2) =
      energyM (Measure.pi fun i => uniformIcc (a i) (b i)) J (fun i => ((b i - a i) / π) ^ 2) := by
    ext p q
    have e3 : ∀ i, ∫ x, J (κ x) p i * J (κ x) q i ∂(Measure.pi fun i => uniformIcc (a i) (b i)) =
        ∫ x, J x p i * J x q i ∂(Measure.pi fun i => uniformIcc (a i) (b i)) := fun i =>
      integral_congr_ae (hfix.mono fun x hx => by simp only [hx])
    simp only [energyM, e3]
  rwa [hcov, hen] at hmain

end Box

end MatrixForm

end ActiveSubspace
