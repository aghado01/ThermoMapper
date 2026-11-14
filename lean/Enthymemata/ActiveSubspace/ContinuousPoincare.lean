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
-/

open MeasureTheory Set Complex intervalIntegral
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

end ActiveSubspace
