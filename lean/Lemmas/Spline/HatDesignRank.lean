import Lemmas.Spline.DesignRank
import Mathlib.Tactic.FieldSimp
import Mathlib.Tactic.Ring

/-! # Schoenberg–Whitney sufficiency for order-two B-splines -/

open Matrix Finset

namespace Spline

/-- Right-continuous degree-one B-spline. A repeated knot contributes no zero-length span. -/
noncomputable def hat (y : ℕ → ℝ) (j : ℕ) (x : ℝ) : ℝ :=
  if y j ≤ x ∧ x < y (j + 1) then (x - y j) / (y (j + 1) - y j)
  else if y (j + 1) ≤ x ∧ x < y (j + 2) then (y (j + 2) - x) / (y (j + 2) - y (j + 1))
  else 0

def hatSupport {n : ℕ} (y : ℕ → ℝ) (i : Fin n) (x : ℝ) : Prop :=
  InSupport (fun j : Fin n ↦ y j = y (j + 1)) (fun j ↦ y j) (fun j ↦ y (j + 2)) i x

theorem hat_support {n : ℕ} {y : ℕ → ℝ} (hy : Monotone y) {i : Fin n} {x : ℝ}
    (hx : hat y i x ≠ 0) : hatSupport y i x := by
  unfold hat at hx
  unfold hatSupport InSupport
  split_ifs at hx with h h
  · exact ⟨Or.inl (lt_of_le_of_ne h.1 (fun he ↦ hx (by simp [he]))),
      h.2.trans_le (hy (by omega))⟩
  · have hle : y i ≤ y (i + 1) := hy (by omega)
    rcases lt_or_eq_of_le (hle.trans h.1) with he | he
    · exact ⟨Or.inl he, h.2⟩
    · exact ⟨Or.inr ⟨by linarith, he⟩, h.2⟩
  · contradiction

/-- Values at knot indices, padded with zero before the first and after the last basis peak. -/
noncomputable def hatCoefficient {n : ℕ} (c : Fin n → ℝ) (k : ℕ) : ℝ :=
  ∑ j : Fin n, if j.val + 1 = k then c j else 0

theorem hatCoefficient_at {n : ℕ} (c : Fin n → ℝ) (i : Fin n) :
    hatCoefficient c (i.val + 1) = c i := by
  simp [hatCoefficient, Fin.val_inj]

theorem hatCoefficient_boundary {n : ℕ} (c : Fin n → ℝ) :
    hatCoefficient c 0 = 0 ∧ hatCoefficient c (n + 1) = 0 := by
  constructor
  · simp [hatCoefficient]
  · apply sum_eq_zero
    intro i _
    rw [if_neg (by omega)]

/-- The two nonzero pieces on one knot span have the ordinary linear interpolation weights. -/
theorem hat_on_span {y : ℕ → ℝ} (hy : Monotone y) {k : ℕ} {x : ℝ}
    (hx : y k ≤ x ∧ x < y (k + 1)) (j : ℕ) :
    (y (k + 1) - y k) * hat y j x =
      (if j + 1 = k then y (k + 1) - x else 0) + (if j = k then x - y k else 0) := by
  have hspan (p : ℕ) : y p ≤ x ∧ x < y (p + 1) ↔ p = k := by
    constructor
    · rintro ⟨hl, hr⟩
      rcases lt_trichotomy p k with h | h | h
      · have := hy (show p + 1 ≤ k by omega)
        linarith [hx.1]
      · exact h
      · have := hy (show k + 1 ≤ p by omega)
        linarith [hx.2]
    · rintro rfl
      exact hx
  have hd : y (k + 1) - y k ≠ 0 := sub_ne_zero.mpr (hx.1.trans_lt hx.2).ne'
  unfold hat
  simp only [hspan j, hspan (j + 1)]
  by_cases h : j = k
  · subst k
    simp [hd, mul_div_cancel₀]
  · by_cases h' : j + 1 = k
    · subst k
      simp [hd, mul_div_cancel₀]
    · simp [h, h']

theorem hat_interpolation {n : ℕ} {y : ℕ → ℝ} (hy : Monotone y) (c : Fin n → ℝ)
    {k : ℕ} {x : ℝ} (hx : y k ≤ x ∧ x < y (k + 1)) :
    (y (k + 1) - y k) * (∑ j : Fin n, hat y j x * c j) =
      (y (k + 1) - x) * hatCoefficient c k + (x - y k) * hatCoefficient c (k + 1) := by
  simp only [hatCoefficient, mul_sum, ← sum_add_distrib]
  apply sum_congr rfl
  intro j _
  rw [← mul_assoc, hat_on_span hy hx, add_mul]
  by_cases h : j.val + 1 = k <;> by_cases h' : j.val = k <;> simp [h, h']

/-- Two distinct zeros on a nonempty affine span force both of its endpoint values to vanish. -/
theorem affine_two_zeros {l u a b p q : ℝ} (hlu : l < u) (hab : a < b)
    (ha : (u - a) * p + (a - l) * q = 0) (hb : (u - b) * p + (b - l) * q = 0) :
    p = 0 ∧ q = 0 := by
  have h : (b - a) * (p - q) = 0 := by linarith only [ha, hb]
  have he := sub_eq_zero.mp ((mul_eq_zero.mp h).resolve_left (sub_ne_zero.mpr hab.ne'))
  rw [he] at ha
  have hz : (u - l) * q = 0 := by linarith only [ha]
  have hq := (mul_eq_zero.mp hz).resolve_left (sub_ne_zero.mpr hlu.ne')
  exact ⟨he.trans hq, hq⟩

theorem affine_zero_propagates {l u x p q : ℝ} (hl : l ≤ x) (hu : x < u) (hp : p ≠ 0)
    (he : (u - x) * p + (x - l) * q = 0) : l < x ∧ q ≠ 0 := by
  have h : (x - l) * q ≠ 0 := fun hz ↦
    mul_ne_zero (sub_ne_zero.mpr hu.ne') hp (by rwa [hz, add_zero] at he)
  exact ⟨lt_of_le_of_ne hl fun h' ↦ h (by rw [h', sub_self, zero_mul]),
    fun h' ↦ h (by rw [h', mul_zero])⟩

theorem h461_hat {n : ℕ} {y : ℕ → ℝ} (hy : Monotone y) (s : Fin n → ℝ) (hs : StrictMono s)
    (hsite : ∀ i, hatSupport y i (s i)) :
    (Matrix.of fun i j : Fin n ↦ hat y j (s i)).det ≠ 0 := by
  cases n
  · simp
  rename_i m
  apply isUnit_iff_ne_zero.mp
  apply (isUnit_iff_isUnit_det _).mp
  apply mulVec_injective_iff_isUnit.mp
  apply (injective_iff_map_eq_zero (Matrix.mulVecLin _)).mpr
  intro c hc
  have hrow (i : Fin (m + 1)) {k : ℕ} (h : y k ≤ s i ∧ s i < y (k + 1)) :
      (y (k + 1) - s i) * hatCoefficient c k + (s i - y k) * hatCoefficient c (k + 1) = 0 := by
    rw [← hat_interpolation hy c h]
    have hi : ∑ j : Fin (m + 1), hat y j (s i) * c j = 0 := congrFun hc i
    rw [hi, mul_zero]
  have right (i : Fin (m + 1)) (hi : hatCoefficient c (i.val + 1) ≠ 0) :
      y (i + 1) < s i ∧ hatCoefficient c (i.val + 2) ≠ 0 := by
    induction i using Fin.induction with
    | zero =>
      have hz := (hatCoefficient_boundary c).1
      have ht := hsite 0
      change (y 0 < s 0 ∨ (y 0 = y 1 ∧ y 0 = s 0)) ∧ s 0 < y 2 at ht
      have hr : y 1 ≤ s 0 := by
        by_contra! hl
        have hlt : y 0 < s 0 := by rcases ht.1 with h | ⟨h, h'⟩ <;> linarith
        have he := hrow 0 ⟨hlt.le, hl⟩
        rw [hz, mul_zero, zero_add] at he
        exact hi ((mul_eq_zero.mp he).resolve_left (sub_ne_zero.mpr hlt.ne'))
      exact affine_zero_propagates hr ht.2 hi (hrow 0 (k := 1) ⟨hr, ht.2⟩)
    | succ i ih =>
      have ht := hsite i.succ
      change (y (i.val+1) < s i.succ ∨ _) ∧ s i.succ < y (i.val+3) at ht
      have hr : y (i.val+2) ≤ s i.succ := by
        by_contra! hl
        have hlt : y (i.val+1) < s i.succ := by
          rcases ht.1 with h | ⟨h, h'⟩
          · exact h
          · change y (i.val+1) = y (i.val+2) at h
            change y (i.val+1) = s i.succ at h'
            linarith
        have he := hrow i.succ ⟨hlt.le, hl⟩
        by_cases hp : hatCoefficient c (i.val+1) = 0
        · rw [hp, mul_zero, zero_add] at he
          exact hi ((mul_eq_zero.mp he).resolve_left (sub_ne_zero.mpr hlt.ne'))
        · obtain ⟨hprev, _⟩ := ih hp
          have hprev' := (hsite i.castSucc).2
          have he' := hrow i.castSucc (k := i.val+1) ⟨hprev.le, hprev'⟩
          exact hi (affine_two_zeros (hlt.trans hl) (hs Fin.castSucc_lt_succ) he' he).2
      exact affine_zero_propagates hr ht.2 hi (hrow i.succ (k := i.val+2) ⟨hr, ht.2⟩)
  funext i
  rw [← hatCoefficient_at c i]
  induction i using Fin.reverseInduction with
  | last =>
    by_contra hn
    exact (right _ hn).2 (hatCoefficient_boundary c).2
  | cast i ih =>
    by_contra hn
    exact (right _ hn).2 ih

open Classical in
theorem hat_fullRank_iff_greedy {n N : ℕ} {y : ℕ → ℝ} (hy : Monotone y)
    (t : Fin N → ℝ) (ht : StrictMono t) :
    (∀ c : Fin n → ℝ, (Matrix.of fun i j ↦ hat y j.val (t i)) *ᵥ c = 0 → c = 0) ↔
      ∀ k < n, (greedy (liftPred (n := n) fun j i ↦ hatSupport y j (t i)) N k).isSome := by
  apply fullRank_iff_greedy t ht _ _ _ _
  · exact fun _ _ h ↦ hy h
  · exact fun _ _ h ↦ hy (by omega)
  · intro i j hij he hclose
    by_cases h : i = j
    · simpa [h] using hclose
    · have h₁ : i.val + 1 ≤ j.val := by omega
      exact le_antisymm (hy (by omega)) ((hy h₁).trans he.ge)
  · exact fun _ _ ↦ hat_support hy
  · exact h461_hat hy

/-- Both selected sites are strictly inside the same span of the clamped linear basis. -/
theorem hat_two_sites_witness :
    (Matrix.of fun i j : Fin 2 ↦ hat (fun k ↦ if k < 2 then 0 else 1) j
      (![1 / 4, 3 / 4] i)).det = 1 / 2 := by
  norm_num [Matrix.det_fin_two, hat]

#print axioms hat_support
#print axioms hatCoefficient_at
#print axioms hatCoefficient_boundary
#print axioms hat_on_span
#print axioms hat_interpolation
#print axioms affine_two_zeros
#print axioms affine_zero_propagates
#print axioms h461_hat
#print axioms hat_fullRank_iff_greedy
#print axioms hat_two_sites_witness

end Spline
