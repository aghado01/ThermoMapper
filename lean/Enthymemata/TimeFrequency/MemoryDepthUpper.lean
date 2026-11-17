import Enthymemata.TimeFrequency.MemoryDepthTesting
import Mathlib.Probability.ProbabilityMassFunction.Constructions
import Mathlib.Probability.ProbabilityMassFunction.Integrals
import Mathlib.Probability.Moments.SubGaussian

/-! # Sufficient budgets for the finite taper-bit experiment -/

namespace TimeFrequency.ContextTree.DepthEvidence

open Finset MeasureTheory ProbabilityTheory
open scoped ENNReal NNReal

variable {X Y : Type*} [Fintype X] [Fintype Y]

noncomputable def Law.pmf (P : Law X) : PMF X :=
  PMF.ofFintype (fun x ↦ ENNReal.ofReal (P.mass x)) (by
    rw [← ENNReal.ofReal_sum_of_nonneg (fun x _ ↦ (P.positive x).le), P.total]
    norm_num)

noncomputable def Law.measure [MeasurableSpace X] (P : Law X) : Measure X := P.pmf.toMeasure

instance [MeasurableSpace X] (P : Law X) : IsProbabilityMeasure P.measure :=
  inferInstanceAs (IsProbabilityMeasure P.pmf.toMeasure)

theorem measure_singleton [MeasurableSpace X] [MeasurableSingletonClass X] (P : Law X) (x : X) :
    P.measure {x} = ENNReal.ofReal (P.mass x) := by
  exact P.pmf.toMeasure_apply_singleton x (measurableSet_singleton x)

theorem integral_expect [MeasurableSpace X] [MeasurableSingletonClass X]
    (P : Law X) (f : X → ℝ) : ∫ x, f x ∂P.measure = expect P f := by
  rw [Law.measure, PMF.integral_eq_sum]
  simp only [Law.pmf, PMF.ofFintype_apply, ENNReal.toReal_ofReal (P.positive _).le,
    smul_eq_mul, expect]

theorem expect_indicator [MeasurableSpace X] [MeasurableSingletonClass X]
    (P : Law X) (s : Set X) : expect P (s.indicator (fun _ ↦ (1 : ℝ))) = P.measure.real s := by
  rw [← integral_expect]
  exact integral_indicator_one (Set.toFinite s).measurableSet

theorem iid_measure [MeasurableSpace X] [MeasurableSingletonClass X] (n : ℕ) (P : Law X) :
    (iid n P).measure = Measure.pi (fun _ : Fin n ↦ P.measure) := by
  apply Measure.ext_of_singleton
  intro x
  rw [measure_singleton, Measure.pi_singleton]
  simp only [measure_singleton, iid]
  exact ENNReal.ofReal_prod_of_nonneg (fun i _ ↦ (P.positive (x i)).le)

noncomputable def bitValue (b : Bool) : ℝ := if b then 1 else 0
noncomputable def bitSum {n : ℕ} (x : Word n) := ∑ i, bitValue (x i)

private theorem bit_mean {p : ℝ} (hp : p ∈ Set.Ioo 0 1) :
    expect (bernoulli p hp) bitValue = p := by
  simp [expect, bernoulli, bitMass, bitValue]

private theorem centered_subgaussian {p : ℝ} (hp : p ∈ Set.Ioo 0 1) (n : ℕ) (i : Fin n) :
    HasSubgaussianMGF (fun x : Word n ↦ bitValue (x i) - p) (1 / 4)
      (iid n (bernoulli p hp)).measure := by
  have hb (x : Word n) : bitValue (x i) ∈ Set.Icc (0 : ℝ) 1 := by
    cases x i <;> norm_num [bitValue]
  have h := hasSubgaussianMGF_of_mem_Icc
    (μ := (iid n (bernoulli p hp)).measure) (measurable_of_finite _).aemeasurable
    (Filter.Eventually.of_forall hb)
  rw [integral_expect, expect_coordinate, bit_mean] at h
  convert h using 1
  norm_num

private theorem centered_independent {p : ℝ} (hp : p ∈ Set.Ioo 0 1) (n : ℕ) :
    iIndepFun (fun i : Fin n ↦ fun x : Word n ↦ bitValue (x i) - p)
      (iid n (bernoulli p hp)).measure := by
  rw [iid_measure]
  exact iIndepFun_pi (X := fun _ : Fin n ↦ fun b : Bool ↦ bitValue b - p)
    (fun _ ↦ (measurable_of_finite _).aemeasurable)

private theorem negative_independent {p : ℝ} (hp : p ∈ Set.Ioo 0 1) (n : ℕ) :
    iIndepFun (fun i : Fin n ↦ fun x : Word n ↦ p - bitValue (x i))
      (iid n (bernoulli p hp)).measure := by
  rw [iid_measure]
  exact iIndepFun_pi (X := fun _ : Fin n ↦ fun b : Bool ↦ p - bitValue b)
    (fun _ ↦ (measurable_of_finite _).aemeasurable)

theorem bit_upper_tail {p ε : ℝ} (hp : p ∈ Set.Ioo 0 1) (he : 0 ≤ ε) (n : ℕ) :
    (iid n (bernoulli p hp)).measure.real {x | n * (p + ε) ≤ bitSum x} ≤
      Real.exp (-2 * n * ε ^ 2) := by
  by_cases hn : n = 0
  · subst n
    simp [bitSum]
  have hnR : (n : ℝ) ≠ 0 := by exact_mod_cast hn
  have h := HasSubgaussianMGF.measure_sum_ge_le_of_iIndepFun
    (centered_independent hp n) (s := univ)
    (fun i _ ↦ centered_subgaussian hp n i)
    (show 0 ≤ (n : ℝ) * ε by positivity)
  have hs : {x : Word n | (n : ℝ) * ε ≤ ∑ i, (bitValue (x i) - p)} =
      {x | (n : ℝ) * (p + ε) ≤ bitSum x} := by
    ext x
    simp only [Set.mem_setOf_eq, sum_sub_distrib, sum_const, card_univ,
      Fintype.card_fin, nsmul_eq_mul, bitSum]
    constructor <;> intro h <;> linarith
  simp only [sum_const, card_univ, Fintype.card_fin, nsmul_eq_mul] at h
  rw [hs] at h
  convert h using 1
  norm_num
  field_simp
  ring

theorem bit_lower_tail {p ε : ℝ} (hp : p ∈ Set.Ioo 0 1) (he : 0 ≤ ε) (n : ℕ) :
    (iid n (bernoulli p hp)).measure.real {x | bitSum x ≤ n * (p - ε)} ≤
      Real.exp (-2 * n * ε ^ 2) := by
  by_cases hn : n = 0
  · subst n
    simp [bitSum]
  have hnR : (n : ℝ) ≠ 0 := by exact_mod_cast hn
  have hg (i : Fin n) :
      HasSubgaussianMGF (fun x : Word n ↦ p - bitValue (x i)) (1 / 4)
        (iid n (bernoulli p hp)).measure := by
    have he : (-(fun x : Word n ↦ bitValue (x i) - p)) =
        (fun x : Word n ↦ p - bitValue (x i)) := by funext x; simp
    rw [← he]
    exact (centered_subgaussian hp n i).neg
  have h := HasSubgaussianMGF.measure_sum_ge_le_of_iIndepFun
    (negative_independent hp n) (s := univ) (fun i _ ↦ hg i)
    (show 0 ≤ (n : ℝ) * ε by positivity)
  have hs : {x : Word n | (n : ℝ) * ε ≤ ∑ i, (p - bitValue (x i))} =
      {x | bitSum x ≤ (n : ℝ) * (p - ε)} := by
    ext x
    simp only [Set.mem_setOf_eq, sum_sub_distrib, sum_const, card_univ,
      Fintype.card_fin, nsmul_eq_mul, bitSum]
    constructor <;> intro h <;> linarith
  simp only [sum_const, card_univ, Fintype.card_fin, nsmul_eq_mul] at h
  rw [hs] at h
  convert h using 1
  norm_num
  field_simp
  ring

noncomputable def taperDecision (n : ℕ) (a : ℝ) : Word n → ℝ :=
  {x | (n : ℝ) * ((2 + a) / 16) ≤ bitSum x}.indicator (fun _ ↦ 1)

theorem taper_decision_range (n : ℕ) (a : ℝ) (x : Word n) :
    taperDecision n a x ∈ Set.Icc 0 1 := by
  classical
  unfold taperDecision
  by_cases h : (n : ℝ) * ((2 + a) / 16) ≤ bitSum x <;> simp [h]

theorem midpoint_errors {a : ℝ} (ha : 0 < a) (hu : a < 1 / 2) (n : ℕ) :
    expect (iid n (bernoulli (1 / 8) (taper_parameters ha hu).1)) (taperDecision n a) ≤
        Real.exp (-(n : ℝ) * a ^ 2 / 128) ∧
      expect (iid n (bernoulli ((1 + a) / 8) (taper_parameters ha hu).2))
        (fun x ↦ 1 - taperDecision n a x) ≤ Real.exp (-(n : ℝ) * a ^ 2 / 128) := by
  have he : -2 * (n : ℝ) * (a / 16) ^ 2 = -(n : ℝ) * a ^ 2 / 128 := by ring
  constructor
  · have h := bit_upper_tail (taper_parameters ha hu).1 (show 0 ≤ a / 16 by positivity) n
    have hm : (1 / 8 : ℝ) + a / 16 = (2 + a) / 16 := by ring
    rw [hm, he] at h
    simpa only [taperDecision, expect_indicator] using h
  · have h := bit_lower_tail (taper_parameters ha hu).2 (show 0 ≤ a / 16 by positivity) n
    have hm : (1 + a) / 8 - a / 16 = (2 + a) / 16 := by ring
    rw [hm, he] at h
    have hc : (fun x : Word n ↦ 1 - taperDecision n a x) =
        {x | bitSum x < (n : ℝ) * ((2 + a) / 16)}.indicator (fun _ ↦ (1 : ℝ)) := by
      classical
      funext x
      by_cases hx : (n : ℝ) * ((2 + a) / 16) ≤ bitSum x <;>
        simp_all [taperDecision, not_le]
    rw [hc, expect_indicator]
    have hmono :
        (iid n (bernoulli ((1 + a) / 8) (taper_parameters ha hu).2)).measure.real
          {x : Word n | bitSum x < (n : ℝ) * ((2 + a) / 16)} ≤
        (iid n (bernoulli ((1 + a) / 8) (taper_parameters ha hu).2)).measure.real
          {x : Word n | bitSum x ≤ (n : ℝ) * ((2 + a) / 16)} :=
      measureReal_mono (by
        intro x hx
        change bitSum x < (n : ℝ) * ((2 + a) / 16) at hx
        change bitSum x ≤ (n : ℝ) * ((2 + a) / 16)
        exact hx.le)
    exact hmono.trans h

theorem common_expect (R : Law Y) (P : Law X) (f : X → ℝ) :
    expect (product R P) (fun z ↦ f z.2) = expect P f := by
  simp only [expect, product, Fintype.sum_prod_type, mul_assoc, ← mul_sum, ← sum_mul,
    R.total, one_mul]

theorem sufficient_budget {a δ : ℝ} (ha : 0 < a) (hu : a < 1 / 2)
    (hδ : 0 < δ) (n : ℕ) (R : Law Y)
    (hn : 128 * Real.log (1 / δ) / a ^ 2 ≤ n) :
    expect (product R (iid n (bernoulli (1 / 8) (taper_parameters ha hu).1)))
        (fun z ↦ taperDecision n a z.2) ≤ δ ∧
      expect (product R (iid n (bernoulli ((1 + a) / 8) (taper_parameters ha hu).2)))
        (fun z ↦ 1 - taperDecision n a z.2) ≤ δ := by
  have h := (div_le_iff₀ (sq_pos_of_pos ha)).mp hn
  have he : Real.exp (-(n : ℝ) * a ^ 2 / 128) ≤ δ := by
    apply (Real.le_log_iff_exp_le hδ).mp
    simp only [one_div, Real.log_inv] at h
    nlinarith
  rw [common_expect R _ (taperDecision n a),
    common_expect R _ (fun x ↦ 1 - taperDecision n a x)]
  exact ⟨(midpoint_errors ha hu n).1.trans he, (midpoint_errors ha hu n).2.trans he⟩

theorem sufficient_budget_rounded {a δ : ℝ} (ha : 0 < a) (hu : a < 1 / 2)
    (hδ : 0 < δ) (n : ℕ) (R : Law Y)
    (hn : ⌈128 * Real.log (1 / δ) / a ^ 2⌉₊ ≤ n) :
    expect (product R (iid n (bernoulli (1 / 8) (taper_parameters ha hu).1)))
        (fun z ↦ taperDecision n a z.2) ≤ δ ∧
      expect (product R (iid n (bernoulli ((1 + a) / 8) (taper_parameters ha hu).2)))
        (fun z ↦ 1 - taperDecision n a z.2) ≤ δ := by
  apply sufficient_budget ha hu hδ n R
  exact (Nat.le_ceil _).trans (Nat.cast_le.mpr hn)

#print axioms measure_singleton
#print axioms integral_expect
#print axioms expect_indicator
#print axioms iid_measure
#print axioms bit_upper_tail
#print axioms bit_lower_tail
#print axioms taper_decision_range
#print axioms midpoint_errors
#print axioms common_expect
#print axioms sufficient_budget
#print axioms sufficient_budget_rounded

end TimeFrequency.ContextTree.DepthEvidence
