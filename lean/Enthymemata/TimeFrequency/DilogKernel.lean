import Enthymemata.TimeFrequency.LogPowerCovariance
import Enthymemata.TimeFrequency.LogLaplace
import Mathlib.MeasureTheory.Integral.DominatedConvergence

/-! # Evaluation of the bivariate log-power covariance kernel -/

open Real Set MeasureTheory Filter Topology

namespace TimeFrequency

noncomputable def laplaceBasis (n : ℕ) (t : ℝ) := (t / (1 + t)) ^ n / (1 + t) ^ 2

theorem laplaceBasis_spec (n : ℕ) :
    Integrable (laplaceBasis n) positiveLine ∧
      (∫ t, laplaceBasis n t ∂positiveLine) = 1 / ((n : ℝ) + 1) := by
  let g := fun t : ℝ ↦ (t / (1+t)) ^ (n+1) / ((n : ℝ)+1)
  have hd {t : ℝ} (ht : 0 ≤ t) : HasDerivAt g (laplaceBasis n t) t := by
    have h := (hasDerivAt_id t).div ((hasDerivAt_id t).const_add 1) (by dsimp; positivity)
    refine ((h.pow (n+1)).div_const _).congr_deriv ?_
    simp only [Pi.div_apply, id_eq, Nat.cast_succ, Nat.add_one_sub_one, one_mul, mul_one]
    unfold laplaceBasis
    field_simp
    ring
  have hp {t : ℝ} (ht : t ∈ Ioi 0) : 0 ≤ laplaceBasis n t := by
    have ht' : 0 < t := ht
    unfold laplaceBasis
    positivity
  have hq : Tendsto (fun t : ℝ ↦ t / (1+t)) atTop (𝓝 1) := by
    have h₀ : Tendsto (fun t : ℝ ↦ (1+t)⁻¹) atTop (𝓝 0) :=
      tendsto_inv_atTop_zero.comp (tendsto_atTop_add_const_left atTop (1 : ℝ) tendsto_id)
    have h := (tendsto_const_nhds (x := (1 : ℝ))).sub h₀
    apply Tendsto.congr' (f₁ := fun t : ℝ ↦ 1 - (1+t)⁻¹) _ (by simpa using h)
    filter_upwards [eventually_gt_atTop (0 : ℝ)] with t ht
    field_simp
    ring
  have hg : Tendsto g atTop (𝓝 (1 / ((n:ℝ)+1))) := by
    simpa [g] using (hq.pow (n+1)).div_const ((n:ℝ)+1)
  refine ⟨integrableOn_Ioi_deriv_of_nonneg' (fun t ht ↦ hd ht) (fun t ht ↦ hp ht) hg, ?_⟩
  simpa [g] using integral_Ioi_of_hasDerivAt_of_nonneg' (fun t ht ↦ hd ht) (fun t ht ↦ hp ht) hg

noncomputable def laplaceCovKernel (κ s t : ℝ) : ℝ :=
  κ / ((1+s)*(1+t)*((1+s)*(1+t)-κ*s*t))

theorem laplaceCovKernel_series {κ s t : ℝ} (hk : κ ∈ Icc 0 1) (hs : 0 < s) (ht : 0 < t) :
    HasSum (fun n ↦ κ^(n+1) * laplaceBasis n s * laplaceBasis n t) (laplaceCovKernel κ s t) := by
  let r := κ * (s/(1+s)) * (t/(1+t))
  have hs' : s / (1+s) < 1 := (div_lt_one (by positivity)).mpr (by linarith)
  have ht' : t / (1+t) < 1 := (div_lt_one (by positivity)).mpr (by linarith)
  have hr : |r| < 1 := by
    have hk0 := hk.1
    rw [abs_of_nonneg (show 0 ≤ r by dsimp [r]; positivity)]
    exact mul_lt_one_of_nonneg_of_lt_one_right
      (mul_lt_one_of_nonneg_of_lt_one_right hk.2 (by positivity) hs').le (by positivity) ht'
  have h := (hasSum_geometric_of_abs_lt_one hr).mul_left (κ/((1+s)^2*(1+t)^2))
  have he : κ / ((1+s)^2*(1+t)^2) * (1-r)⁻¹ = laplaceCovKernel κ s t := by
    have hr' := (abs_lt.mp hr).2
    have hn := sub_ne_zero.mpr (ne_of_gt (show 1 > r from hr'))
    have hd : (1+s)*(1+t)-κ*s*t ≠ 0 := by
      apply ne_of_gt
      nlinarith only [mul_nonneg (sub_nonneg.mpr hk.2) (mul_nonneg hs.le ht.le), hs, ht]
    dsimp [r] at hn ⊢
    unfold laplaceCovKernel
    field_simp
  rw [← he]
  apply h.congr
  intro s'
  apply Finset.sum_congr rfl
  intro n _
  simp only [r, laplaceBasis, mul_pow, pow_succ, div_eq_mul_inv, mul_inv_rev]
  ring

theorem positive_scales : ∀ᵐ p : ℝ × ℝ ∂positiveLine.prod positiveLine, 0 < p.1 ∧ 0 < p.2 := by
  rw [positiveLine, Measure.prod_restrict]
  exact ae_restrict_mem (measurableSet_Ioi.prod measurableSet_Ioi)

theorem integral_laplaceCovKernel {κ : ℝ} (hk : κ ∈ Icc 0 1) :
    Integrable (fun p : ℝ × ℝ ↦ laplaceCovKernel κ p.1 p.2) (positiveLine.prod positiveLine) ∧
    HasSum (fun n : ℕ ↦ κ^(n+1) / ((n : ℝ)+1)^2)
      (∫ p, laplaceCovKernel κ p.1 p.2 ∂positiveLine.prod positiveLine) := by
  let F := fun n (p : ℝ × ℝ) ↦ κ^(n+1) * (laplaceBasis n p.1 * laplaceBasis n p.2)
  have hI (n) : Integrable (F n) (positiveLine.prod positiveLine) :=
    ((laplaceBasis_spec n).1.mul_prod (laplaceBasis_spec n).1).const_mul _
  have hval (n) : ∫ p, F n p ∂positiveLine.prod positiveLine = κ^(n+1) / ((n:ℝ)+1)^2 := by
    simp only [F, integral_const_mul, integral_prod_mul, (laplaceBasis_spec n).2]
    field_simp
  have hpos (n) : 0 ≤ᵐ[positiveLine.prod positiveLine] F n := by
    filter_upwards [positive_scales] with p hp
    obtain ⟨hs, ht⟩ := hp
    have hk' := hk.1
    dsimp [F, laplaceBasis]
    positivity
  have hsum : Summable fun n ↦ ∫ p, ‖F n p‖ ∂positiveLine.prod positiveLine := by
    have hnorm (n) : (∫ p, ‖F n p‖ ∂positiveLine.prod positiveLine) = κ^(n+1)/((n:ℝ)+1)^2 := by
      rw [← hval]
      exact integral_congr_ae ((hpos n).mono fun p hp ↦ Real.norm_of_nonneg hp)
    simp_rw [hnorm]
    apply Summable.of_nonneg_of_le (fun n ↦ by have := hk.1; positivity) _ dilog_one.summable
    intro n
    exact div_le_div_of_nonneg_right (pow_le_pow_left₀ hk.1 hk.2 _) (sq_nonneg _)
  have he : (fun p ↦ ∑' n, F n p) =ᵐ[positiveLine.prod positiveLine]
      (fun p : ℝ × ℝ ↦ laplaceCovKernel κ p.1 p.2) := by
    filter_upwards [positive_scales] with p hp
    simpa only [F, mul_assoc] using (laplaceCovKernel_series hk hp.1 hp.2).tsum_eq
  have hresult := hasSum_integral_of_summable_integral_norm hI hsum
  simp_rw [hval] at hresult
  rw [integral_congr_ae he] at hresult
  refine ⟨?_, hresult⟩
  rcases eq_or_lt_of_le hk.1 with hk0 | hk0
  · simp [← hk0, laplaceCovKernel]
  · apply Integrable.of_integral_ne_zero
    rw [← hresult.tsum_eq]
    exact (hresult.summable.tsum_pos (fun n ↦ by positivity) (0 : ℕ) (by positivity)).ne'

#print axioms laplaceBasis_spec
#print axioms laplaceCovKernel_series
#print axioms positive_scales
#print axioms integral_laplaceCovKernel

end TimeFrequency
