import Enthymemata.TimeFrequency.ContextTreeFourier
import Mathlib.Analysis.SpecialFunctions.Log.Basic
import Mathlib.Algebra.BigOperators.Field
import Mathlib.Algebra.BigOperators.Ring.Finset

/-! # Finite observed likelihoods for the all-depth source pair -/

namespace TimeFrequency.ContextTree.DepthEvidence

open Finset

/-- Strict positivity is imposed on the retained finite support. -/
structure Law (X : Type*) [Fintype X] where
  mass : X → ℝ
  positive : ∀ x, 0 < mass x
  total : ∑ x, mass x = 1

variable {X Y : Type*} [Fintype X] [Fintype Y]

noncomputable def expect (P : Law X) (f : X → ℝ) := ∑ x, P.mass x * f x
noncomputable def likelihood (P Q : Law X) (x : X) := P.mass x / Q.mass x
noncomputable def divergence (P Q : Law X) := expect P (fun x ↦ Real.log (likelihood P Q x))

noncomputable def product (P : Law X) (Q : Law Y) : Law (X × Y) where
  mass x := P.mass x.1 * Q.mass x.2
  positive x := mul_pos (P.positive _) (Q.positive _)
  total := by
    simp only [Fintype.sum_prod_type, ← mul_sum, Q.total, mul_one, P.total]

noncomputable def iid (n : ℕ) (P : Law X) : Law (Fin n → X) where
  mass x := ∏ i, P.mass (x i)
  positive x := prod_pos (fun i _ ↦ P.positive (x i))
  total := by
    rw [← Fintype.prod_sum]
    simp only [P.total, prod_const_one]

theorem expect_coordinate (n : ℕ) (P : Law X) (i : Fin n) (f : X → ℝ) :
    expect (iid n P) (fun x ↦ f (x i)) = expect P f := by
  let g (j : Fin n) (x : X) := P.mass x * if j = i then f x else 1
  have factor (x : Fin n → X) :
      (∏ j, g j (x j)) = (∏ j, P.mass (x j)) * f (x i) := by
    simp only [g, prod_mul_distrib, Fintype.prod_ite_eq']
  have terms (j : Fin n) : (∑ x, g j x) = if j = i then expect P f else 1 := by
    by_cases h : j = i <;> simp [g, h, expect, P.total]
  unfold expect
  change (∑ x : Fin n → X, (∏ j, P.mass (x j)) * f (x i)) = _
  simp_rw [← factor]
  rw [← Fintype.prod_sum]
  simp only [terms, Fintype.prod_ite_eq']
  rfl

theorem iid_log_likelihood (n : ℕ) (P Q : Law X) (x : Fin n → X) :
    Real.log (likelihood (iid n P) (iid n Q) x) =
      ∑ i, Real.log (likelihood P Q (x i)) := by
  unfold likelihood iid
  rw [← prod_div_distrib, Real.log_prod]
  intro i _
  exact ne_of_gt (div_pos (P.positive _) (Q.positive _))

theorem iid_divergence (n : ℕ) (P Q : Law X) :
    divergence (iid n P) (iid n Q) = n * divergence P Q := by
  unfold divergence
  simp only [expect, iid_log_likelihood, mul_sum]
  rw [sum_comm]
  change (∑ i : Fin n, expect (iid n P) (fun x ↦ Real.log (likelihood P Q (x i)))) = _
  have hi (i : Fin n) :
      expect (iid n P) (fun x ↦ Real.log (likelihood P Q (x i))) =
        expect P (fun x ↦ Real.log (likelihood P Q x)) :=
    expect_coordinate n P i (fun x ↦ Real.log (likelihood P Q x))
  simp_rw [hi]
  simp only [sum_const, card_univ, Fintype.card_fin, nsmul_eq_mul, expect, mul_sum]

theorem common_likelihood (R : Law Y) (P Q : Law X) (y : Y) (x : X) :
    likelihood (product R P) (product R Q) (y, x) = likelihood P Q x := by
  unfold likelihood product
  exact mul_div_mul_left _ _ (ne_of_gt (R.positive y))

theorem common_divergence (R : Law Y) (P Q : Law X) :
    divergence (product R P) (product R Q) = divergence P Q := by
  simp only [divergence, expect, Fintype.sum_prod_type, common_likelihood]
  change (∑ y, ∑ x, R.mass y * P.mass x * Real.log (likelihood P Q x)) = _
  simp only [mul_assoc, ← mul_sum, ← sum_mul, R.total, one_mul]

noncomputable def bitMass (p : ℝ) (b : Bool) : ℝ := if b then p else 1 - p

noncomputable def bernoulli (p : ℝ) (hp : p ∈ Set.Ioo 0 1) : Law Bool where
  mass := bitMass p
  positive b := by cases b <;> simp [bitMass] <;> linarith [hp.1, hp.2]
  total := by simp [bitMass ]

noncomputable def kl (p q : ℝ) :=
  p * Real.log (p / q) + (1 - p) * Real.log ((1 - p) / (1 - q))

theorem bernoulli_divergence {p q : ℝ} (hp : p ∈ Set.Ioo 0 1) (hq : q ∈ Set.Ioo 0 1) :
    divergence (bernoulli p hp) (bernoulli q hq) = kl p q := by
  simp [divergence, expect, likelihood, bernoulli, bitMass, kl ]

def successes {n : ℕ} (x : Word n) : ℕ := (univ.filter (fun i ↦ x i = true)).card

theorem iid_bit_mass {n : ℕ} {p : ℝ} (hp : p ∈ Set.Ioo 0 1) (x : Word n) :
    (iid n (bernoulli p hp)).mass x = p ^ successes x * (1 - p) ^ (n - successes x) := by
  have hc : (univ.filter (fun i : Fin n ↦ ¬ x i = true)).card = n - successes x := by
    have h := card_filter_add_card_filter_not (s := (univ : Finset (Fin n)))
      (fun i ↦ x i = true)
    simp only [card_univ, Fintype.card_fin] at h
    unfold successes
    omega
  simp only [iid, bernoulli, bitMass]
  rw [prod_ite]
  simp only [prod_const, hc, successes]

theorem taper_parameters {a : ℝ} (ha : 0 < a) (hu : a < 1 / 2) :
    (1 / 8 : ℝ) ∈ Set.Ioo 0 1 ∧ (1 + a) / 8 ∈ Set.Ioo 0 1 := by
  constructor <;> constructor <;> linarith

theorem taper_likelihood {a : ℝ} (ha : 0 < a) (hu : a < 1 / 2) (n : ℕ)
    (R : Law Y) (y : Y) (x : Word n) :
    likelihood
      (product R (iid n (bernoulli ((1 + a) / 8) (taper_parameters ha hu).2)))
      (product R (iid n (bernoulli (1 / 8) (taper_parameters ha hu).1))) (y, x) =
      (1 + a) ^ successes x * (1 - a / 7) ^ (n - successes x) := by
  rw [common_likelihood]
  simp only [likelihood, iid_bit_mass]
  rw [mul_div_mul_comm, ← div_pow, ← div_pow]
  congr 1 <;> congr 1 <;> ring

theorem zero_tapers {a : ℝ} (ha : 0 < a) (hu : a < 1 / 2)
    (R : Law Y) (y : Y) (x : Word 0) :
    likelihood
      (product R (iid 0 (bernoulli ((1 + a) / 8) (taper_parameters ha hu).2)))
      (product R (iid 0 (bernoulli (1 / 8) (taper_parameters ha hu).1))) (y, x) = 1 := by
  rw [taper_likelihood ha hu]
  simp [successes]

theorem expected_log_evidence {a : ℝ} (ha : 0 < a) (hu : a < 1 / 2)
    (n : ℕ) (R : Law Y) :
    divergence
      (product R (iid n (bernoulli ((1 + a) / 8) (taper_parameters ha hu).2)))
      (product R (iid n (bernoulli (1 / 8) (taper_parameters ha hu).1))) =
      n * kl ((1 + a) / 8) (1 / 8) := by
  rw [common_divergence, iid_divergence, bernoulli_divergence]

private theorem expect_div (P : Law X) (f : X → ℝ) (c : ℝ) :
    expect P (fun x ↦ f x / c) = expect P f / c := by
  simp only [expect, ← mul_div_assoc, ← sum_div]

theorem expected_bits {a : ℝ} (ha : 0 < a) (hu : a < 1 / 2) (n : ℕ) (R : Law Y) :
    let P := product R (iid n (bernoulli ((1 + a) / 8) (taper_parameters ha hu).2))
    let Q := product R (iid n (bernoulli (1 / 8) (taper_parameters ha hu).1))
    expect P (fun x ↦ Real.log (likelihood P Q x) / Real.log 2) =
      n * (kl ((1 + a) / 8) (1 / 8) / Real.log 2) := by
  dsimp only
  rw [expect_div]
  simpa only [divergence, mul_div_assoc] using
    congrArg (fun z : ℝ ↦ z / Real.log 2) (expected_log_evidence ha hu n R)

noncomputable def posterior (u v : ℝ) (P Q : Law X) (x : X) :=
  u * P.mass x / (u * P.mass x + v * Q.mass x)

theorem posterior_odds {u v : ℝ} (hu : 0 < u) (hv : 0 < v) (P Q : Law X) (x : X) :
    posterior u v P Q x / (1 - posterior u v P Q x) = u / v * likelihood P Q x := by
  have hd : u * P.mass x + v * Q.mass x ≠ 0 :=
    ne_of_gt (add_pos (mul_pos hu (P.positive x)) (mul_pos hv (Q.positive x)))
  have hc : 1 - posterior u v P Q x =
      v * Q.mass x / (u * P.mass x + v * Q.mass x) := by
    unfold posterior
    field_simp
    ring
  rw [hc]
  unfold posterior likelihood
  field_simp

#print axioms expect_coordinate
#print axioms iid_log_likelihood
#print axioms iid_divergence
#print axioms common_likelihood
#print axioms common_divergence
#print axioms bernoulli_divergence
#print axioms iid_bit_mass
#print axioms taper_parameters
#print axioms taper_likelihood
#print axioms zero_tapers
#print axioms expected_log_evidence
#print axioms expected_bits
#print axioms posterior_odds

end TimeFrequency.ContextTree.DepthEvidence
