import Lemmas.TimeFrequency.PhasePropagation
import Mathlib.Combinatorics.SimpleGraph.Hasse
import Mathlib.Analysis.SpecialFunctions.Trigonometric.Basic

/-! # Finite real phase recovery from one- and two-sample views -/

namespace TimeFrequency

open Complex SimpleGraph

/-- Unnormalized DFT with the negative-exponent convention on one complete real frame. -/
noncomputable def frameDFT {n : ℕ} (x : Fin n → ℝ) (k : ℕ) : ℂ :=
  ∑ j, (x j : ℂ) * Complex.exp (-(2 * Real.pi * ((k : ℂ) * j.val / n) * I))

noncomputable def framePower {n : ℕ} (x : Fin n → ℝ) (k : ℕ) : ℝ :=
  Complex.normSq (frameDFT x k)

theorem framePower_dc {n : ℕ} (x : Fin n → ℝ) : framePower x 0 = (∑ j, x j) ^ 2 := by
  simp [framePower, frameDFT, ← Complex.ofReal_sum, Complex.normSq_ofReal, sq]

theorem twoSample_powers (a b : ℝ) :
    framePower ![a,b] 0 = (a+b)^2 ∧ framePower ![a,b] 1 = (a-b)^2 := by
  constructor
  · rw [framePower_dc, Fin.sum_univ_two]
    rfl
  · have h : frameDFT ![a,b] 1 = ((a-b : ℝ) : ℂ) := by
      simp [frameDFT, Fin.sum_univ_two]
      ring_nf
      simp only [← neg_mul]
      rw [neg_mul, Complex.exp_neg_pi_mul_I]
      ring
    rw [framePower, h, Complex.normSq_ofReal]
    ring

theorem twoSample_power_difference (a b : ℝ) :
    framePower ![a,b] 0 - framePower ![a,b] 1 = 4*a*b := by
  rw [(twoSample_powers a b).1, (twoSample_powers a b).2]
  ring

/-- Singleton magnitudes and both bins of every complete adjacent pair, without wrap or padding. -/
def SameTwoSampleViews {n : ℕ} (x y : Fin n → ℝ) : Prop :=
  (∀ i, |y i| = |x i|) ∧ ∀ u v : Fin n, u.val + 1 = v.val → ∀ k : Fin 2,
    framePower ![y u, y v] k.val = framePower ![x u, x v] k.val

theorem twoSampleViews_iff_products {n : ℕ} {x y : Fin n → ℝ} :
    SameTwoSampleViews x y ↔ (∀ i, |y i| = |x i|) ∧
      ∀ u v, (pathGraph n).Adj u v → y u * y v = x u * x v := by
  constructor
  · rintro ⟨hm, hp⟩
    have h (u v : Fin n) (ha : u.val + 1 = v.val) : y u * y v = x u * x v := by
      have he := congrArg₂ (· - ·) (hp u v ha 0) (hp u v ha 1)
      change framePower _ 0 - framePower _ 1 = framePower _ 0 - framePower _ 1 at he
      rw [twoSample_power_difference, twoSample_power_difference] at he
      linarith
    refine ⟨hm, fun u v ha ↦ ?_⟩
    rcases pathGraph_adj.mp ha with ha | ha
    · exact h u v ha
    · simpa [mul_comm] using h v u ha
  · rintro ⟨hm, hp⟩
    refine ⟨hm, fun u v huv k ↦ ?_⟩
    have h := hp u v (pathGraph_adj.mpr (.inl huv))
    have hu := (sq_eq_sq_iff_abs_eq_abs _ _).mpr (hm u)
    have hv := (sq_eq_sq_iff_abs_eq_abs _ _).mpr (hm v)
    fin_cases k
    · change framePower _ 0 = framePower _ 0
      rw [(twoSample_powers _ _).1, (twoSample_powers _ _).1]
      nlinarith only [hu,hv,h]
    · change framePower _ 1 = framePower _ 1
      rw [(twoSample_powers _ _).2, (twoSample_powers _ _).2]
      nlinarith only [hu,hv,h]

def supportPath {n : ℕ} (x : Fin n → ℝ) : SimpleGraph {i // x i ≠ 0} :=
  (pathGraph n).induce {i | x i ≠ 0}

/-- One sign per connected component, with the original zero samples fixed. -/
def ComponentSigns {n : ℕ} (x y : Fin n → ℝ) : Prop :=
  (∀ i, x i = 0 → y i = 0) ∧ ∃ σ : (supportPath x).ConnectedComponent → ℝ,
    (∀ C, σ C = 1 ∨ σ C = -1) ∧ ∀ v : {i // x i ≠ 0},
      y v.val = σ ((supportPath x).connectedComponentMk v) * x v.val

theorem twoSampleViews_iff_componentSigns {n : ℕ} {x y : Fin n → ℝ} :
    SameTwoSampleViews x y ↔ ComponentSigns x y := by
  rw [twoSampleViews_iff_products]
  constructor
  · rintro ⟨hm, hp⟩
    let r := fun v : {i // x i ≠ 0} ↦ y v.val / x v.val
    have hr (v) : |r v| = 1 := by rw [abs_div, hm, div_self (abs_ne_zero.mpr v.property)]
    have hstep {u v} (h : (supportPath x).Adj u v) : r u = r v := by
      have hu : y u.val = r u * x u.val := (div_mul_cancel₀ _ u.property).symm
      have he : y u.val * y v.val = x u.val * x v.val := hp _ _ h
      have hv := phase_step (c := (r u : ℂ)) (xu := (x u.val : ℂ)) (xw := (x v.val : ℂ))
        (yu := (y u.val : ℂ)) (yw := (y v.val : ℂ)) (by simpa using hr u)
        (by exact_mod_cast u.property) (by exact_mod_cast hu)
        (by simpa only [Complex.conj_ofReal, ← Complex.ofReal_mul] using congrArg Complex.ofReal he)
      exact ((eq_div_iff v.property).mpr (by exact_mod_cast hv.symm))
    have hw {u v} (p : (supportPath x).Walk u v) : r u = r v := by
      induction p with
      | nil => rfl
      | cons he _ ih => exact (hstep he).trans ih
    let σ := SimpleGraph.ConnectedComponent.lift (G := supportPath x) r fun _ _ p _ ↦ hw p
    refine ⟨fun i hi ↦ abs_eq_zero.mp (by simpa [hi] using hm i), σ, ?_, ?_⟩
    · intro C
      induction C using SimpleGraph.ConnectedComponent.ind with
      | h v =>
        change r v = 1 ∨ r v = -1
        exact abs_eq_abs.mp (by simpa using hr v)
    · intro v
      exact (div_mul_cancel₀ _ v.property).symm
  · rintro ⟨hz, σ, hσ, hv⟩
    constructor
    · intro i
      by_cases hi : x i = 0
      · simp [hi, hz i hi]
      · have h := hv ⟨i, hi⟩
        rcases hσ ((supportPath x).connectedComponentMk ⟨i,hi⟩) with hs | hs <;>
          simp [h, hs]
    · intro u v he
      by_cases hu : x u = 0
      · simp [hu, hz u hu]
      by_cases hh : x v = 0
      · simp [hh, hz v hh]
      have ha : (supportPath x).Adj ⟨u,hu⟩ ⟨v,hh⟩ := he
      have hc := SimpleGraph.ConnectedComponent.connectedComponentMk_eq_of_adj ha
      have h₁ := hv ⟨u,hu⟩
      have h₂ := hv ⟨v,hh⟩
      rw [hc] at h₁
      rcases hσ ((supportPath x).connectedComponentMk ⟨v,hh⟩) with hs | hs <;>
        simp [h₁,h₂,hs]

theorem global_sign_of_connected_support {n : ℕ} {x y : Fin n → ℝ}
    (hc : (supportPath x).Connected) (h : SameTwoSampleViews x y) : y = x ∨ y = -x := by
  obtain ⟨hz, σ, hσ, hv⟩ := twoSampleViews_iff_componentSigns.mp h
  obtain ⟨v⟩ := hc.nonempty
  let C := (supportPath x).connectedComponentMk v
  have he (u : {i // x i ≠ 0}) : (supportPath x).connectedComponentMk u = C :=
    SimpleGraph.ConnectedComponent.sound (hc.preconnected u v)
  have hall (i : Fin n) : y i = σ C * x i := by
    by_cases hi : x i = 0
    · simp [hi, hz i hi]
    · simpa only [he] using hv ⟨i, hi⟩
  rcases hσ C with hs | hs
  · exact Or.inl (funext fun i ↦ by simpa [hs] using hall i)
  · exact Or.inr (funext fun i ↦ by simpa [hs] using hall i)

theorem zero_signal_views {n : ℕ} {y : Fin n → ℝ} : SameTwoSampleViews 0 y ↔ y = 0 := by
  constructor
  · exact fun h ↦ funext fun i ↦ abs_eq_zero.mp (by simpa using h.1 i)
  · rintro rfl
    exact ⟨fun _ ↦ rfl, fun _ _ _ _ ↦ rfl⟩

theorem single_sample_views {x y : Fin 1 → ℝ} : SameTwoSampleViews x y ↔ y = x ∨ y = -x := by
  rw [twoSampleViews_iff_products]
  constructor
  · intro h
    rcases abs_eq_abs.mp (h.1 0) with hs | hs
    · exact Or.inl (funext fun i ↦ by fin_cases i; exact hs)
    · exact Or.inr (funext fun i ↦ by fin_cases i; exact hs)
  · rintro (rfl | rfl) <;> simp

theorem gap_witness (a b : ℝ) (hab : a * b ≠ 0) :
    SameTwoSampleViews ![a,0,b] ![a,0,-b] ∧
    ¬((![a,0,-b] : Fin 3 → ℝ) = ![a,0,b] ∨ (![a,0,-b] : Fin 3 → ℝ) = -![a,0,b]) ∧
    framePower ![a,0,b] 0 - framePower ![a,0,-b] 0 = 4*a*b := by
  refine ⟨twoSampleViews_iff_products.mpr ⟨?_, ?_⟩, ?_, ?_⟩
  · intro i
    fin_cases i <;> simp
  · intro u v huv
    rw [pathGraph_adj] at huv
    fin_cases u <;> fin_cases v <;> simp_all
  · rintro (h | h)
    · have he := congrFun h 2
      change -b = b at he
      exact hab (by rw [show b = 0 by linarith, mul_zero])
    · have he := congrFun h 0
      change a = -a at he
      exact hab (by rw [show a = 0 by linarith, zero_mul])
  · simp [framePower_dc, Fin.sum_univ_succ]
    ring

#print axioms framePower_dc
#print axioms twoSample_powers
#print axioms twoSample_power_difference
#print axioms twoSampleViews_iff_products
#print axioms twoSampleViews_iff_componentSigns
#print axioms global_sign_of_connected_support
#print axioms zero_signal_views
#print axioms single_sample_views
#print axioms gap_witness

end TimeFrequency
