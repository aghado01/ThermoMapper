import Enthymemata.TimeFrequency.ContextTreeFourier

/-! # Completing the depth-four binary observation design -/

namespace TimeFrequency.ContextTree.DepthFour

open Finset

/-- Excluding r indexes all five four-site indicator tapers. A excludes 3 and B excludes 2. -/
def quartic (r : Fin 5) (x : Word 5) : ℤ := ∏ i ∈ univ.erase r, sign (x i)
def taper (r : Fin 5) (x : Word 5) : ℤ := (∑ i ∈ univ.erase r, sign (x i)) ^ 2
def selectedPairs (r : Fin 5) (x : Word 5) : ℤ :=
  ∑ i, ∑ j, if i < j ∧ i ≠ r ∧ j ≠ r then sign (x i) * sign (x j) else 0

theorem taper_is_power (r : Fin 5) (x : Word 5) :
    framePower (fun i ↦ if i = r then 0 else (sign (x i) : ℝ)) 0 = (taper r x : ℝ) := by
  rw [framePower_dc]
  simp only [taper, Int.cast_pow, Int.cast_sum]
  congr 1
  have he (i : Fin 5) : (if i = r then 0 else (sign (x i) : ℝ)) =
      if i ∈ univ.erase r then (sign (x i) : ℝ) else 0 := by
    by_cases hi : i = r <;> simp [hi]
  simp_rw [he]
  rw [← sum_filter]
  simp only [filter_mem_eq_inter, univ_inter]

theorem tapered_fourth (r : Fin 5) (x : Word 5) :
    8 * hit (taper r x) 16 = 1 + selectedPairs r x + quartic r x := by decide +revert

theorem five_top (x : Word 5) :
    16 * hit (dc x) 25 = 1 + pairSum (fun _ ↦ true) x + ∑ r, quartic r x := by decide +revert

private theorem inversion_at_one (x : Word 5) :
    16 * (hit x (fun _ ↦ true) + hit x (fun _ ↦ false)) =
      1 + pairSum (fun _ ↦ true) x + ∑ r, quartic r x := by decide +revert

theorem five_inversion (t x : Word 5) :
    16 * (hit x t + hit x (flip t)) =
      1 + pairSum t x + ∑ r, quartic r t * quartic r x := by
  let y : Word 5 := fun i ↦ decide (x i = t i)
  have hy (i) : sign (y i) = sign (t i) * sign (x i) := by
    dsimp [y]
    cases ht : t i <;> cases hx : x i <;> decide
  have heq : y = (fun _ ↦ true) ↔ x = t := by
    constructor
    · intro he; funext i
      have hi := congrFun he i
      simpa [y] using hi
    · intro he; subst x; funext i; simp [y]
  have hflip : y = (fun _ ↦ false) ↔ x = flip t := by
    constructor
    · intro he; funext i
      have hi := congrFun he i
      dsimp [y] at hi
      change x i = !(t i)
      cases hx : x i <;> cases ht : t i <;> simp_all
    · intro he
      subst x
      funext i
      simp [y, flip]
  have hp : pairSum (fun _ ↦ true) y = pairSum t x := by
    unfold pairSum
    apply sum_congr rfl
    intro i _
    apply sum_congr rfl
    intro j _
    by_cases h : i < j
    · simp only [if_pos h, show sign true = 1 from rfl, one_mul]
      rw [hy, hy]
      ring
    · simp only [if_neg h]
  have hq (r) : quartic r y = quartic r t * quartic r x := by
    simp only [quartic, hy, prod_mul_distrib]
  have he := inversion_at_one y
  simpa only [hit, heq, hflip, hp, hq] using he

theorem endpoint_quartics (x : Word 5) :
    quartic 4 x = character (slice x 0 4 (by decide)) ∧
    quartic 0 x = character (slice x 1 4 (by decide)) := by decide +revert

private theorem mean_linear {n : ℕ} (p : Word n → ℝ) (a b : ℤ) (f g : Word n → ℤ) :
    mean p (fun x ↦ a * f x + b * g x) = a * mean p f + b * mean p g := by
  simp only [mean, Int.cast_add, Int.cast_mul, mul_add, sum_add_distrib, mul_sum]
  congr 1 <;> apply sum_congr rfl <;> intros <;> ring

private theorem agree_linear {n : ℕ} {P Q : Source} {f g : Word n → ℤ}
    (a b : ℤ) (hf : Agree P Q f) (hg : Agree P Q g) :
    Agree P Q (fun x ↦ a * f x + b * g x) := by
  unfold Agree at *
  rw [mean_linear, mean_linear, hf, hg]

private theorem agree_constant (P Q : Source) (n : ℕ) (c : ℤ) :
    Agree P Q (fun _ : Word n ↦ c) := by
  simp [Agree, mean, ← sum_mul, P.total, Q.total]

private theorem agree_sum {n : ℕ} {I : Type*} [Fintype I] {P Q : Source}
    (f : I → Word n → ℤ) (h : ∀ i, Agree P Q (f i)) :
    Agree P Q (fun x ↦ ∑ i, f i x) := by
  have hs (p : Word n → ℝ) : mean p (fun x ↦ ∑ i, f i x) = ∑ i, mean p (f i) := by
    simp only [mean, Int.cast_sum, mul_sum]
    rw [sum_comm]
  unfold Agree at *
  rw [hs, hs]
  exact sum_congr rfl (fun i _ ↦ h i)

private theorem mean_from_bins {n : ℕ} {p q : Word n → ℝ} (f : Word n → ℤ)
    (h : ∀ v, bin p f v = bin q f v) (g : ℤ → ℤ) :
    mean p (g ∘ f) = mean q (g ∘ f) := by
  rw [← bin_mean, ← bin_mean]
  exact sum_congr rfl (fun y _ ↦ by rw [h])

theorem lags_identified {P Q : Source} (h : SameDC 4 P Q) {k : ℕ} (hk : 1 ≤ k) (h4 : k ≤ 4) :
    lag P k = lag Q k := by
  have hm (w : ℕ) (hw : w ≤ 5) : powerMean P w = powerMean Q w := by
    cases w with
    | zero => simp [powerMean, mean, dc]
    | succ w => exact mean_from_bins dc (h _ (by omega) hw) id
  obtain ⟨j, rfl⟩ := Nat.exists_eq_succ_of_ne_zero (by omega : k ≠ 0)
  have hp := mean_second_difference P j
  have hq := mean_second_difference Q j
  rw [hm _ (by omega), hm _ (by omega), hm _ (by omega)] at hp
  linarith

private theorem pairs_identified {P Q : Source} (h : SameDC 4 P Q) (i j : Fin 5) (hij : i < j) :
    Agree P Q (fun x : Word 5 ↦ sign (x i) * sign (x j)) := by
  unfold Agree
  rw [mean_pair P i j hij, mean_pair Q i j hij]
  exact lags_identified h (by omega) (by omega)

private theorem selectedPairs_identified {P Q : Source} (h : SameDC 4 P Q) (r : Fin 5) :
    Agree P Q (selectedPairs r) := by
  apply agree_sum
  intro i
  apply agree_sum
  intro j
  by_cases he : i < j ∧ i ≠ r ∧ j ≠ r
  · simpa only [if_pos he] using pairs_identified h i j he.1
  · simpa only [if_neg he] using agree_constant P Q 5 0

private theorem pairSum_identified {P Q : Source} (h : SameDC 4 P Q) (t : Word 5) :
    Agree P Q (pairSum t) := by
  apply agree_sum
  intro i
  apply agree_sum
  intro j
  by_cases he : i < j
  · simp only [if_pos he]
    simpa only [zero_mul, add_zero] using agree_linear (sign (t i) * sign (t j)) 0
      (pairs_identified h i j he) (agree_constant P Q 5 0)
  · simpa only [if_neg he] using agree_constant P Q 5 0

private theorem added_quartic {P Q : Source} (h : SameDC 4 P Q) (r : Fin 5)
    (hr : ∀ v, bin (P.law 5) (taper r) v = bin (Q.law 5) (taper r) v) :
    Agree P Q (quartic r) := by
  have he : quartic r = fun x ↦ 1 * (8 * hit (taper r x) 16 - selectedPairs r x) + (-1) * 1 := by
    funext x
    linarith [tapered_fourth r x]
  rw [he]
  have hh := mean_from_bins _ hr (fun v ↦ hit v 16)
  have ha := agree_linear 8 (-1) hh (selectedPairs_identified h r)
  simpa only [Function.comp_apply, neg_one_mul, sub_eq_add_neg] using
      agree_linear 1 (-1) ha (agree_constant P Q 5 1)

private theorem old_quartic {P Q : Source} (h : SameDC 4 P Q) :
    Agree P Q (quartic 0) ∧ Agree P Q (quartic 4) := by
  have hm : P.law 4 = Q.law 4 :=
    small_block_identified h (by decide) (by decide) (by decide)
  have he : mean (P.law 4) character = mean (Q.law 4) character := by rw [hm]
  have ht (o : ℕ) (ho : o + 4 ≤ 5) :
      Agree P Q (fun x : Word 5 ↦ character (slice x o 4 ho)) := by
    unfold Agree
    rw [P.translation, Q.translation]
    exact he
  exact ⟨by simpa only [← (endpoint_quartics _).2] using ht 1 (by decide),
    by simpa only [← (endpoint_quartics _).1] using ht 0 (by decide)⟩

def SameAugmented (P Q : Source) : Prop :=
  SameDC 4 P Q ∧
    (∀ v, bin (P.law 5) (taper 3) v = bin (Q.law 5) (taper 3) v) ∧
    (∀ v, bin (P.law 5) (taper 2) v = bin (Q.law 5) (taper 2) v)

theorem quartics_identified {P Q : Source} (h : SameAugmented P Q) (r : Fin 5) :
    Agree P Q (quartic r) := by
  have h0 := (old_quartic h.1).1
  have h4 := (old_quartic h.1).2
  have h3 := added_quartic h.1 3 h.2.1
  have h2 := added_quartic h.1 2 h.2.2
  have ht := mean_from_bins dc (h.1 5 (by decide) (by decide)) (fun v ↦ hit v 25)
  have hsum : Agree P Q (fun x ↦ ∑ r, quartic r x) := by
    have he : (fun x : Word 5 ↦ ∑ r, quartic r x) =
        fun x ↦ 1 * (16 * hit (dc x) 25 - pairSum (fun _ ↦ true) x) + (-1) * 1 := by
      funext x; linarith [five_top x]
    rw [he]
    have ha := agree_linear 16 (-1) ht (pairSum_identified h.1 (fun _ ↦ true))
    simpa only [Function.comp_apply, neg_one_mul, sub_eq_add_neg] using
      agree_linear 1 (-1) ha (agree_constant P Q 5 1)
  have h1 : Agree P Q (quartic 1) := by
    have hu : (fun x : Word 5 ↦ ∑ r, quartic r x) =
        fun x ↦ quartic 0 x + quartic 1 x + quartic 2 x + quartic 3 x + quartic 4 x := by
      funext x; simp [Fin.sum_univ_succ]; ring
    unfold Agree at *
    rw [hu] at hsum
    have hadd (p : Word 5 → ℝ) (f g : Word 5 → ℤ) :
        mean p (fun x ↦ f x + g x) = mean p f + mean p g := by
      simpa only [one_mul, Int.cast_one] using mean_linear p 1 1 f g
    simp only [hadd] at hsum
    linarith
  fin_cases r <;> assumption

theorem completed_block {P Q : Source} (h : SameAugmented P Q) : P.law 5 = Q.law 5 := by
  funext t
  have hp := pairSum_identified h.1 t
  have hq : Agree P Q (fun x ↦ ∑ r, quartic r t * quartic r x) := by
    apply agree_sum
    intro r
    simpa only [zero_mul, add_zero] using
      agree_linear (quartic r t) 0 (quartics_identified h r) (agree_constant P Q 5 0)
  have hall := agree_linear 1 1 (agree_linear 1 1 (agree_constant P Q 5 1) hp) hq
  simp only [one_mul] at hall
  unfold Agree at hall
  simp_rw [← five_inversion t] at hall
  have he (p : Source) :
      mean (p.law 5) (fun x ↦ 16 * (hit x t + hit x (flip t))) = 32 * p.law 5 t := by
    simp [mean, hit, mul_add, mul_ite, sum_add_distrib, p.symmetry]
    ring
  rw [he, he] at hall
  linarith

theorem completed_kernel_tree {P Q : Source} (h : SameAugmented P Q) :
    nextKernel 4 (P.law 5) = nextKernel 4 (Q.law 5) ∧
    minimalTree (nextKernel 4 (P.law 5)) = minimalTree (nextKernel 4 (Q.law 5)) := by
  rw [completed_block h]
  exact ⟨rfl, rfl⟩

#print axioms taper_is_power
#print axioms tapered_fourth
#print axioms five_top
#print axioms five_inversion
#print axioms endpoint_quartics
#print axioms lags_identified
#print axioms quartics_identified
#print axioms completed_block
#print axioms completed_kernel_tree

end TimeFrequency.ContextTree.DepthFour
