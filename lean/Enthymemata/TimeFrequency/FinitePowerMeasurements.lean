import Enthymemata.TimeFrequency.TwoSampleRecovery
import Mathlib.Analysis.SpecialFunctions.Log.Basic
import Mathlib.LinearAlgebra.Matrix.DotProduct

/-! # Finite real waveforms observed through specified complex Fourier cells -/

namespace TimeFrequency.FinitePower

open Complex Matrix Finset
open scoped ComplexConjugate

variable {T J : Type*} [Fintype T]

noncomputable def amplitude (a : T → ℂ) (x : T → ℝ) := ∑ t, a t * (x t : ℂ)

noncomputable def powers (a : J → T → ℂ) (x : T → ℝ) (j : J) := normSq (amplitude (a j) x)

noncomputable def logFeatures (ε : ℝ) (a : J → T → ℂ) (x : T → ℝ) (j : J) :=
  Real.log (powers a x j + ε)

theorem logFeatures_eq_iff {ε : ℝ} (hε : 0 < ε) (a : J → T → ℂ) (x y : T → ℝ) :
    logFeatures ε a x = logFeatures ε a y ↔ powers a x = powers a y := by
  constructor
  · intro h
    funext j
    have he := congrArg Real.exp (congrFun h j)
    change Real.exp (Real.log (powers a x j + ε)) = Real.exp (Real.log (powers a y j + ε)) at he
    have hp (z : T → ℝ) : 0 < powers a z j + ε :=
      add_pos_of_nonneg_of_pos (normSq_nonneg _) hε
    rwa [Real.exp_log (hp x), Real.exp_log (hp y), add_left_inj] at he
  · intro h
    funext j
    simp [logFeatures,h]

def outer (x : T → ℝ) : Matrix T T ℝ := fun t s ↦ x t*x s

def cellMatrix (a : T → ℂ) : Matrix T T ℝ := fun t s ↦ (a t).re*(a s).re + (a t).im*(a s).im

noncomputable def lifted (a : J → T → ℂ) (H : Matrix T T ℝ) (j : J) :=
  ∑ t, ∑ s, cellMatrix (a j) t s * H t s

theorem powers_lifted (a : J → T → ℂ) (x : T → ℝ) : powers a x = lifted a (outer x) := by
  funext j
  simp only [powers, amplitude, lifted, cellMatrix, outer, Complex.normSq_apply, Complex.re_sum,
    Complex.im_sum, Complex.mul_re, Complex.mul_im, Complex.ofReal_re, Complex.ofReal_im,
    mul_zero, sub_zero, zero_add, sum_mul, mul_sum, ← sum_add_distrib]
  apply sum_congr rfl
  intro t _
  apply sum_congr rfl
  intro s _
  ring

theorem power_quadratic (a : T → ℂ) (x : T → ℝ) :
    normSq (amplitude a x) = x ⬝ᵥ (cellMatrix a *ᵥ x) := by
  have h := congrFun (powers_lifted (fun _ : Unit ↦ a) x) ()
  refine h.trans ?_
  simp only [lifted, outer, dotProduct, mulVec, mul_sum]
  apply sum_congr rfl
  intro t _
  apply sum_congr rfl
  intro s _
  ring

omit [Fintype T] in
theorem outer_eq_iff (x y : T → ℝ) : outer x = outer y ↔ y = x ∨ y = -x := by
  constructor
  · intro h
    have hxy (u v) : x u*x v = y u*y v := congrFun (congrFun h u) v
    by_cases hz : ∀ u, x u = 0
    · left
      funext v
      have hh := hxy v v
      rw [hz] at hh ⊢
      nlinarith
    · push Not at hz
      obtain ⟨u,hu⟩ := hz
      have hs : (y u - x u)*(y u+x u) = 0 := by nlinarith [hxy u u]
      rcases mul_eq_zero.mp hs with hp | hm
      · left
        funext v
        have hsame := sub_eq_zero.mp hp
        apply mul_left_cancel₀ hu
        simpa only [hsame] using (hxy u v).symm
      · right
        funext v
        have hneg : y u = -x u := by linarith
        have he := hxy u v
        rw [hneg] at he
        change y v = -x v
        apply mul_left_cancel₀ hu
        linarith
  · rintro (rfl | rfl)
    · rfl
    · ext i j
      simp [outer]

/-- Injectivity up to sign on the candidate class is equivalent to excluding its nonzero
differences of outer products from the lifted measurement kernel. -/
theorem kernel_criterion (a : J → T → ℂ) (C : Set (T → ℝ)) :
    (∀ x ∈ C, ∀ y ∈ C, powers a x = powers a y → y = x ∨ y = -x) ↔
      ∀ x ∈ C, ∀ y ∈ C, lifted a (outer x - outer y) = 0 → outer x - outer y = 0 := by
  have hl (x y : T → ℝ) :
      lifted a (outer x - outer y) = lifted a (outer x) - lifted a (outer y) := by
    funext j
    simp [lifted, mul_sub, sum_sub_distrib]
  simp_rw [hl, sub_eq_zero, ← powers_lifted, outer_eq_iff]

theorem amplitude_add (a : T → ℂ) (x y : T → ℝ) :
    amplitude a (x+y) = amplitude a x + amplitude a y := by
  simp [amplitude, mul_add, sum_add_distrib]

theorem amplitude_sub (a : T → ℂ) (x y : T → ℝ) :
    amplitude a (x-y) = amplitude a x - amplitude a y := by
  simp [amplitude, mul_sub, sum_sub_distrib]

theorem cut_difference (a : J → T → ℂ) (u v : T → ℝ) (j : J) :
    powers a (u+v) j - powers a (u-v) j = 4*((amplitude (a j) u)*conj (amplitude (a j) v)).re := by
  simp only [powers, amplitude_add, amplitude_sub, Complex.normSq_apply, Complex.add_re,
    Complex.add_im, Complex.sub_re, Complex.sub_im, Complex.mul_re, Complex.conj_re,
    Complex.conj_im]
  ring

theorem cut_visible_iff (a : J → T → ℂ) (u v : T → ℝ) :
    powers a (u+v) ≠ powers a (u-v) ↔
      ∃ j, ((amplitude (a j) u)*conj (amplitude (a j) v)).re ≠ 0 := by
  simp only [Function.ne_iff]
  apply exists_congr
  intro j
  have h := cut_difference a u v j
  constructor <;> intro hn he <;> apply hn <;> linarith

theorem missed_piece (a : J → T → ℂ) (u v : T → ℝ)
    (h : ∀ j, (∀ t, a j t = 0 ∨ u t = 0) ∨ (∀ t, a j t = 0 ∨ v t = 0)) :
    powers a (u+v) = powers a (u-v) := by
  have hz (b : T → ℂ) (z : T → ℝ) (hb : ∀ t, b t = 0 ∨ z t = 0) : amplitude b z = 0 := by
    apply sum_eq_zero
    intro t _
    rcases hb t with h | h <;> simp [h]
  funext j
  have hd := cut_difference a u v j
  rcases h j with hu | hv
  · simpa [hz _ _ hu, sub_eq_zero] using hd
  · simpa [hz _ _ hv, sub_eq_zero] using hd

/-- The real and imaginary components implement the pilot's negative-exponent Fourier row.
A window can encode its translation and boundary extension by its values on the recording. -/
noncomputable def fourierRow (g θ : T → ℝ) (t : T) : ℂ :=
  ⟨g t * Real.cos (θ t), -g t * Real.sin (θ t)⟩

theorem amplitude_impulse [DecidableEq T] (a : T → ℂ) (t : T) (b : ℝ) :
    amplitude a (Pi.single t b) = a t * b := by
  classical
  rw [amplitude, sum_eq_single t]
  · simp
  · intro u _ hu
    simp [Pi.single_eq_of_ne hu]
  · simp

theorem impulse_cross [DecidableEq T] (g θ : T → ℝ) (t s : T) (b c : ℝ) :
    ((amplitude (fourierRow g θ) (Pi.single t b)) *
      conj (amplitude (fourierRow g θ) (Pi.single s c))).re =
      b*c*g t*g s*Real.cos (θ t-θ s) := by
  rw [amplitude_impulse, amplitude_impulse]
  simp only [fourierRow, Complex.mul_re, Complex.mul_im, Complex.conj_re, Complex.conj_im,
    Complex.ofReal_re, Complex.ofReal_im, Real.cos_sub]
  ring

theorem impulse_visibility [DecidableEq T] (g θ : J → T → ℝ) (t s : T) (b c : ℝ) :
    powers (fun j ↦ fourierRow (g j) (θ j)) (Pi.single t b + Pi.single s c) ≠
      powers (fun j ↦ fourierRow (g j) (θ j)) (Pi.single t b - Pi.single s c) ↔
    ∃ j, b*c*g j t*g j s*Real.cos (θ j t-θ j s) ≠ 0 := by
  rw [cut_visible_iff]
  simp_rw [impulse_cross]

theorem quarter_cycle_and_dc :
    let u := Pi.single (0 : Fin 2) (1 : ℝ)
    let v := Pi.single (1 : Fin 2) (1 : ℝ)
    let a := fun j : Fin 2 ↦ fourierRow (fun _ : Fin 2 ↦ 1)
        (fun t ↦ (j.val : ℝ)*t.val*(Real.pi/2))
    powers a (u+v) 1 = powers a (u-v) 1 ∧
      powers a (u+v) 0 - powers a (u-v) 0 = 4 := by
  dsimp only
  constructor
  · apply sub_eq_zero.mp
    rw [cut_difference, impulse_cross]
    norm_num [Real.cos_neg]
  · rw [cut_difference, impulse_cross]
    norm_num

/-- Every frequency and every positive common log floor is blind to an impulse at
a zero window coefficient, even if the signal class includes the zero waveform. -/
theorem blind_endpoint [DecidableEq T] (a : J → T → ℂ) (t : T) (ha : ∀ j, a j t = 0)
    {ε : ℝ} (hε : 0 < ε) :
    logFeatures ε a (Pi.single t 1) = logFeatures ε a 0 ∧
      (Pi.single t (1 : ℝ) ≠ (0 : T → ℝ) ∧ Pi.single t (1 : ℝ) ≠ -(0 : T → ℝ)) := by
  refine ⟨(logFeatures_eq_iff hε _ _ _).mpr (funext fun j ↦ ?_), ?_⟩
  · unfold powers
    rw [amplitude_impulse, ha]
    simp [amplitude]
  · have h : Pi.single t (1 : ℝ) ≠ (0 : T → ℝ) := by
      intro he
      have := congrFun he t
      simp at this
    exact ⟨h, by rwa [neg_zero]⟩

theorem pilot_hann_blind {N : ℕ} (frequencies : J → ℝ) (scale : ℝ) {ε : ℝ} (hε : 0 < ε) :
    let a := fun j ↦ fourierRow (fun t : Fin (N + 1) ↦
        scale*(1-Real.cos (2*Real.pi*t.val/(N+1)))/2)
      (fun t ↦ frequencies j * t.val)
    logFeatures ε a (Pi.single 0 1) = logFeatures ε a 0 ∧
      (Pi.single 0 1 ≠ (0 : Fin (N+1) → ℝ) ∧ Pi.single 0 1 ≠ -(0 : Fin (N+1) → ℝ)) := by
  apply blind_endpoint _ _ _ hε
  intro j
  apply Complex.ext <;> simp [fourierRow]

example (f : Fin 64 → ℝ) := pilot_hann_blind (N := 511) f 1 (ε := 1) zero_lt_one
example (f : Fin 128 → ℝ) := pilot_hann_blind (N := 1023) f 1 (ε := 1) zero_lt_one

#print axioms logFeatures_eq_iff
#print axioms powers_lifted
#print axioms power_quadratic
#print axioms outer_eq_iff
#print axioms kernel_criterion
#print axioms amplitude_add
#print axioms amplitude_sub
#print axioms cut_difference
#print axioms cut_visible_iff
#print axioms missed_piece
#print axioms amplitude_impulse
#print axioms impulse_cross
#print axioms impulse_visibility
#print axioms quarter_cycle_and_dc
#print axioms blind_endpoint
#print axioms pilot_hann_blind

end TimeFrequency.FinitePower
