import Enthymemata.TimeFrequency.DilogKernel
import Mathlib.Probability.Moments.Covariance

/-! # The log-power covariance from the joint Laplace transform -/

open Real MeasureTheory Set Filter ProbabilityTheory
open scoped ProbabilityTheory

namespace TimeFrequency

variable {Ω : Type*} [MeasurableSpace Ω] {μ : Measure Ω} [IsProbabilityMeasure μ]

def PowerLaplaceLaw (κ : ℝ) (μ : Measure Ω) (X Y : Ω → ℝ) : Prop :=
  ∀ s t : ℝ, 0 ≤ s → 0 ≤ t →
    (∫ ω, exp (-X ω * s) * exp (-Y ω * t) ∂μ) =
      1 / ((1+s)*(1+t)-κ*s*t)

noncomputable def meanLogKernel (t : ℝ) := (exp (-t) - 1 / (1+t)) / t

theorem integrable_meanLogKernel : Integrable meanLogKernel positiveLine := by
  refine ((laplaceBasis_spec 0).1.const_mul 2).mono'
    (by apply Measurable.aestronglyMeasurable; unfold meanLogKernel; fun_prop) ?_
  filter_upwards [ae_restrict_mem measurableSet_Ioi] with t ht
  have ht' : 0 < t := ht
  have hd : 0 < 1+t := by linarith
  have he : exp (-t) ≤ 1 / (1+t) := by
    rw [exp_neg, one_div]
    exact inv_anti₀ hd (by linarith [add_one_le_exp t])
  rw [Real.norm_of_nonpos (show meanLogKernel t ≤ 0 from div_nonpos_of_nonpos_of_nonneg
    (sub_nonpos.mpr he) ht'.le)]
  unfold meanLogKernel laplaceBasis
  rw [pow_zero]
  rcases le_total t 1 with h | h
  · have hl : 1-t ≤ exp (-t) := by linarith [add_one_le_exp (-t)]
    have hh : -( (exp (-t) - 1/(1+t))/t) ≤ t/(1+t) := by
      field_simp
      nlinarith [mul_nonneg (sub_nonneg.mpr hl) hd.le]
    apply hh.trans
    field_simp
    nlinarith [mul_nonneg (sub_nonneg.mpr h) hd.le]
  · have hh : -((exp (-t) - 1/(1+t))/t) ≤ 1/(t*(1+t)) := by
      field_simp
      nlinarith [mul_nonneg (exp_pos (-t)).le hd.le]
    apply hh.trans
    field_simp
    nlinarith

theorem memLp_exp_test {X : Ω → ℝ} (hX : Measurable X) (hp : ∀ᵐ ω ∂μ, 0 ≤ X ω)
    {s : ℝ} (hs : 0 ≤ s) : MemLp (fun ω ↦ exp (-X ω * s)) 2 μ := by
  apply MemLp.of_bound (by fun_prop) 1
  filter_upwards [hp] with ω h
  rw [Real.norm_of_nonneg (exp_pos _).le, exp_le_one_iff]
  exact mul_nonpos_of_nonpos_of_nonneg (neg_nonpos.mpr h) hs

theorem memLp_logKernel_test {X : Ω → ℝ} (hX : Measurable X) (hp : ∀ᵐ ω ∂μ, 0 ≤ X ω)
    {s : ℝ} (hs : 0 ≤ s) : MemLp (fun ω ↦ logKernel (X ω) s) 2 μ := by
  simpa only [logKernel, div_eq_inv_mul, Pi.sub_apply] using
    ((memLp_const (exp (-s))).sub (memLp_exp_test hX hp hs)).const_mul s⁻¹

/-- A nonnegative power with the unit-exponential Laplace transform has no atom at zero. -/
theorem positive_of_unit_laplace {X : Ω → ℝ} (hX : Measurable X)
    (hp : ∀ᵐ ω ∂μ, 0 ≤ X ω) (law : ∀ s : ℝ, 0 ≤ s →
      (∫ ω, exp (-X ω*s) ∂μ) = 1/(1+s)) : ∀ᵐ ω ∂μ, 0 < X ω := by
  let A := {ω | X ω = 0}
  have hm : MeasurableSet A := measurableSet_eq_fun hX measurable_const
  have hzero : μ.real A = 0 := by
    by_contra! hn
    have ha : 0 < μ.real A := lt_of_le_of_ne (by positivity) (Ne.symm hn)
    let s := (μ.real A)⁻¹
    have hs : 0 ≤ s := inv_nonneg.mpr ha.le
    have hle : μ.real A ≤ 1/(1+s) := by
      rw [← integral_indicator_one hm, ← law _ hs]
      apply integral_mono ((integrable_const _).indicator hm)
        ((memLp_exp_test hX hp hs).integrable (by norm_num))
      intro ω
      by_cases h : ω ∈ A
      · rw [Set.indicator_of_mem h]
        change (1 : ℝ) ≤ exp (-X ω*s)
        rw [show X ω = 0 from h, neg_zero, zero_mul, exp_zero]
      · simp [Set.indicator_of_notMem h, (exp_pos _).le]
    have hh := (le_div_iff₀ (show 0 < 1+s by linarith)).mp hle
    rw [mul_add, mul_one, mul_inv_cancel₀ ha.ne'] at hh
    linarith
  have hne : ∀ᵐ ω ∂μ, X ω ≠ 0 := by
    rw [ae_iff]
    simpa [A] using (measureReal_eq_zero_iff (measure_ne_top μ A)).mp hzero
  exact hp.and hne |>.mono fun _ h ↦ lt_of_le_of_ne h.1 h.2.symm

theorem logKernel_moments {κ s t : ℝ} {X Y : Ω → ℝ}
    (law : PowerLaplaceLaw κ μ X Y) (hk : κ ∈ Icc 0 1)
    (hX : Measurable X) (hY : Measurable Y) (pX : ∀ᵐ ω ∂μ, 0 ≤ X ω) (pY : ∀ᵐ ω ∂μ, 0 ≤ Y ω)
    (hs : 0 < s) (ht : 0 < t) :
    (∫ ω, logKernel (X ω) s ∂μ) = meanLogKernel s ∧
    (∫ ω, logKernel (Y ω) t ∂μ) = meanLogKernel t ∧
    (∫ ω, logKernel (X ω) s * logKernel (Y ω) t ∂μ) =
      laplaceCovKernel κ s t + meanLogKernel s * meanLogKernel t := by
  let f := fun ω ↦ exp (-X ω * s)
  let g := fun ω ↦ exp (-Y ω * t)
  have hf := memLp_exp_test hX pX hs.le
  have hg := memLp_exp_test hY pY ht.le
  have hfi := hf.integrable (by norm_num)
  have hgi := hg.integrable (by norm_num)
  have hx : ∫ ω, f ω ∂μ = 1/(1+s) := by simpa [f] using law s 0 hs.le le_rfl
  have hy : ∫ ω, g ω ∂μ = 1/(1+t) := by simpa [g] using law 0 t le_rfl ht.le
  have hx' : ∫ ω, logKernel (X ω) s ∂μ = meanLogKernel s := by
    unfold logKernel meanLogKernel
    rw [integral_div, integral_sub (integrable_const _) hfi, integral_const]
    simp [f] at hx
    simp [hx]
  have hy' : ∫ ω, logKernel (Y ω) t ∂μ = meanLogKernel t := by
    unfold logKernel meanLogKernel
    rw [integral_div, integral_sub (integrable_const _) hgi, integral_const]
    simp [g] at hy
    simp [hy]
  refine ⟨hx',hy',?_⟩
  have hKX := memLp_logKernel_test hX pX hs.le
  have hKY := memLp_logKernel_test hY pY ht.le
  have hc := covariance_eq_sub hKX hKY
  simp only [Pi.mul_apply, hx', hy'] at hc
  suffices hval : cov[fun ω ↦ logKernel (X ω) s, fun ω ↦ logKernel (Y ω) t; μ] =
      laplaceCovKernel κ s t by linarith
  unfold logKernel
  rw [covariance_fun_div_left, covariance_fun_div_right, covariance_const_sub_left hfi,
    covariance_const_sub_right hgi, neg_neg, covariance_eq_sub hf hg]
  change (((∫ ω, f ω * g ω ∂μ) - (∫ ω, f ω ∂μ)*(∫ ω, g ω ∂μ)) / t) / s = _
  rw [hx, hy, show (∫ ω, f ω * g ω ∂μ) = _ from law s t hs.le ht.le]
  have hd : 0 < (1+s)*(1+t)-κ*s*t := by
    nlinarith [mul_nonneg (sub_nonneg.mpr hk.2) (mul_nonneg hs.le ht.le)]
  unfold laplaceCovKernel
  generalize hD : (1+s)*(1+t)-κ*s*t = D at hd ⊢
  field_simp [hd.ne']
  nlinarith only [hD]

omit [IsProbabilityMeasure μ] in
theorem diagonal_laplace {X : Ω → ℝ} (h : ∀ s : ℝ, 0 ≤ s →
    (∫ ω, exp (-X ω*s) ∂μ) = 1/(1+s)) : PowerLaplaceLaw 1 μ X X := by
  intro s t hs ht
  simp_rw [← exp_add, ← mul_add]
  rw [h _ (add_nonneg hs ht)]
  congr 1
  ring

theorem log_memLp_of_diagonal_laplace {X : Ω → ℝ} (hX : Measurable X)
    (hp : ∀ᵐ ω ∂μ, 0 < X ω) (law : PowerLaplaceLaw 1 μ X X) :
    MemLp (fun ω ↦ log (X ω)) 2 μ := by
  let F := fun z : (ℝ × ℝ) × Ω ↦ logKernel (X z.2) z.1.1 * logKernel (X z.2) z.1.2
  have hi : Integrable F ((positiveLine.prod positiveLine).prod μ) := by
    apply (integrable_prod_iff (by dsimp [F, logKernel]; fun_prop)).mpr
    constructor
    · filter_upwards [positive_scales] with p h
      have hK (s : ℝ) (hs : 0 < s) : MemLp (fun ω ↦ logKernel (X ω) s) 2 μ := by
        exact memLp_logKernel_test hX (hp.mono fun _ h ↦ h.le) hs.le
      exact (hK _ h.1).integrable_mul (hK _ h.2)
    · apply ((integral_laplaceCovKernel (κ := 1) (by constructor <;> norm_num)).1.add
        (integrable_meanLogKernel.mul_prod integrable_meanLogKernel)).congr
      filter_upwards [positive_scales] with p h
      obtain ⟨_,_,hm⟩ := logKernel_moments law (by constructor <;> norm_num) hX hX
        (hp.mono fun _ h ↦ h.le) (hp.mono fun _ h ↦ h.le) h.1 h.2
      refine hm.symm.trans (integral_congr_ae ?_)
      filter_upwards [hp] with ω _
      apply (Real.norm_of_nonneg _).symm
      dsimp [F,logKernel]
      rcases le_total 1 (X ω) with hx | hx
      · apply mul_nonneg <;> apply div_nonneg <;> first | exact h.1.le | exact h.2.le |
          apply sub_nonneg.mpr; apply exp_le_exp.mpr; nlinarith [h.1,h.2]
      · apply mul_nonneg_of_nonpos_of_nonpos <;> apply div_nonpos_of_nonpos_of_nonneg <;>
          first | exact h.1.le | exact h.2.le |
            apply sub_nonpos.mpr; apply exp_le_exp.mpr; nlinarith [h.1,h.2]
  apply (memLp_two_iff_integrable_sq (Real.measurable_log.comp hX).aestronglyMeasurable).mpr
  apply hi.integral_prod_right.congr
  filter_upwards [hp] with ω h
  dsimp [F]
  rw [integral_prod_mul (logKernel (X ω)) (logKernel (X ω)), (logKernel_spec h).2, sq]

theorem integral_log_from_kernel {X : Ω → ℝ} (hX : Measurable X)
    (hp : ∀ᵐ ω ∂μ, 0 < X ω) (hi : Integrable (fun ω ↦ log (X ω)) μ) :
    (∫ t, (∫ ω, logKernel (X ω) t ∂μ) ∂positiveLine) = ∫ ω, log (X ω) ∂μ := by
  let F := fun z : Ω × ℝ ↦ logKernel (X z.1) z.2
  have hf : Integrable F (μ.prod positiveLine) := by
    apply (integrable_prod_iff (by
      apply Measurable.aestronglyMeasurable
      dsimp [F,logKernel]
      fun_prop)).mpr
    exact ⟨hp.mono fun _ hx ↦ (logKernel_spec hx).1, hi.norm.congr
      (hp.mono fun _ hx ↦ (integral_norm_logKernel hx).symm)⟩
  refine (integral_prod_symm F hf).symm.trans ?_
  rw [integral_prod F hf]
  exact integral_congr_ae (hp.mono fun _ hx ↦ (logKernel_spec hx).2)

theorem logPower_series_of_laplace {κ : ℝ} (hk : κ ∈ Icc 0 1) {X Y : Ω → ℝ}
    (hX : Measurable X) (hY : Measurable Y) (hpX : ∀ᵐ ω ∂μ, 0 < X ω)
    (hpY : ∀ᵐ ω ∂μ, 0 < Y ω) (law : PowerLaplaceLaw κ μ X Y) :
    HasSum (fun n : ℕ ↦ κ^(n+1)/((n:ℝ)+1)^2)
      ((∫ ω, log (X ω)*log (Y ω) ∂μ) - (∫ ω, log (X ω) ∂μ)*(∫ ω, log (Y ω) ∂μ)) := by
  have hLX := log_memLp_of_diagonal_laplace hX hpX (diagonal_laplace fun s hs ↦
    by simpa using law s 0 hs le_rfl)
  have hLY := log_memLp_of_diagonal_laplace hY hpY (diagonal_laplace fun s hs ↦
    by simpa using law 0 s le_rfl hs)
  have hK {s t : ℝ} (hs : 0 < s) (ht : 0 < t) :=
    logKernel_moments law hk hX hY (hpX.mono fun _ h ↦ h.le) (hpY.mono fun _ h ↦ h.le) hs ht
  have hμX := integral_log_from_kernel hX hpX (hLX.integrable (by norm_num))
  have hμY := integral_log_from_kernel hY hpY (hLY.integrable (by norm_num))
  have hmeanX : (∫ s, meanLogKernel s ∂positiveLine) = ∫ ω, log (X ω) ∂μ :=
    (integral_congr_ae ((ae_restrict_mem measurableSet_Ioi).mono fun s hs ↦
      (hK (s := s) (t := 1) hs zero_lt_one).1.symm)).trans hμX
  have hmeanY : (∫ t, meanLogKernel t ∂positiveLine) = ∫ ω, log (Y ω) ∂μ :=
    (integral_congr_ae ((ae_restrict_mem measurableSet_Ioi).mono fun t ht ↦
      (hK (s := 1) (t := t) zero_lt_one ht).2.1.symm)).trans hμY
  have hxy := (log_product_laplace hX hY hpX hpY (hLX.integrable_mul hLY)).2
  have ⟨hi,hs⟩ := integral_laplaceCovKernel hk
  convert hs using 1
  rw [← hxy, ← hmeanX, ← hmeanY, ← integral_prod_mul, ← integral_sub]
  · apply integral_congr_ae
    filter_upwards [positive_scales] with p hp
    linarith only [(hK hp.1 hp.2).2.2]
  · exact (log_product_laplace hX hY hpX hpY (hLX.integrable_mul hLY)).1
  · exact integrable_meanLogKernel.mul_prod integrable_meanLogKernel

/-- Both logarithms are square-integrable and their covariance is the stated series,
including at the two singular parameter endpoints. -/
theorem covariance_log_of_laplace {κ : ℝ} (hk : κ ∈ Icc 0 1) {X Y : Ω → ℝ}
    (hX : Measurable X) (hY : Measurable Y) (pX : ∀ᵐ ω ∂μ, 0 ≤ X ω)
    (pY : ∀ᵐ ω ∂μ, 0 ≤ Y ω) (law : PowerLaplaceLaw κ μ X Y) :
    HasSum (fun n : ℕ ↦ κ^(n+1)/((n:ℝ)+1)^2)
      (cov[fun ω ↦ log (X ω), fun ω ↦ log (Y ω); μ]) := by
  have mx (s : ℝ) (hs : 0 ≤ s) : (∫ ω, exp (-X ω*s) ∂μ) = 1/(1+s) := by
    simpa using law s 0 hs le_rfl
  have my (s : ℝ) (hs : 0 ≤ s) : (∫ ω, exp (-Y ω*s) ∂μ) = 1/(1+s) := by
    simpa using law 0 s le_rfl hs
  have hx := positive_of_unit_laplace hX pX mx
  have hy := positive_of_unit_laplace hY pY my
  rw [covariance_eq_sub (log_memLp_of_diagonal_laplace hX hx (diagonal_laplace mx))
    (log_memLp_of_diagonal_laplace hY hy (diagonal_laplace my))]
  exact logPower_series_of_laplace hk hX hY hx hy law

#print axioms covariance_log_of_laplace
#print axioms integrable_meanLogKernel
#print axioms memLp_exp_test
#print axioms memLp_logKernel_test
#print axioms positive_of_unit_laplace
#print axioms logKernel_moments
#print axioms diagonal_laplace
#print axioms log_memLp_of_diagonal_laplace
#print axioms integral_log_from_kernel
#print axioms logPower_series_of_laplace

end TimeFrequency
