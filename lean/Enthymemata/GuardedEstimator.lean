import Mathlib.Algebra.BigOperators.Field
import Mathlib.Algebra.Order.BigOperators.Group.Finset
import Mathlib.Algebra.Order.BigOperators.Ring.Finset
import Mathlib.Data.Fintype.BigOperators
import Mathlib.Data.Finset.Powerset
import Mathlib.Data.Nat.Choose.Basic
import Mathlib.Data.Real.Basic
import Mathlib.Tactic.FieldSimp
import Mathlib.Tactic.LinearCombination
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.Ring

/-!
# BASS — estimators that stay valid whatever the emulator does (GE-1, GE-2, GE-3)

The finite statements of `Protolemmata/BASS/guarded-estimator-lemmas.md`. `m` is a run's mean at a
configuration, `f` a fit's prediction, `r = m − f`. Expectations are written out as finite sums,
over simple random samples of a pool or over independent draws from a probability vector. A run's
noise enters only through its conditional mean `m`, so these are the statements after that first
expectation.

* **GE-1.** `crossfit_gap` is the algebra; `srs_mean` is the simple-random-sample expectation; and
  `crossfit_bias` is the conditional bias they give, over simple random samples of the pool outside
  the training set. `srs_difference_unbiased` is the case with no training set on the pool.
* **GE-2.** `guarded_fold_unbiased`: with independent draws and the integral in the first term, each
  fold is exactly unbiased, for any map from training draws to fits.
* **GE-3.** `telescope_srs` and `telescope_iid`, telescoping levels on a pool and on draws;
  `tunedVariance` and `tunedVariance_sub_tuned`, the tuned weight as the minimizer of the stated
  quadratic, with its corollaries; `horvitz_thompson_unbiased`, known inclusion probabilities.
-/

open Finset

namespace BASS

/-! ### Simple random samples -/

section Pool

variable {P : Type*} [DecidableEq P]

/-- The `(k+1)`-subsets of `U` containing `i` are as many as the `k`-subsets of `U` without it. -/
theorem card_powersetCard_filter_mem {U : Finset P} {i : P} (hi : i ∈ U) (k : ℕ) :
    ((U.powersetCard (k + 1)).filter (i ∈ ·)).card = (U.card - 1).choose k := by
  have himage : (U.powersetCard (k + 1)).filter (i ∈ ·) =
      ((U.erase i).powersetCard k).image (insert i) := by
    ext F
    simp only [mem_filter, mem_powersetCard, mem_image]
    constructor
    · rintro ⟨⟨hFU, hc⟩, hiF⟩
      refine ⟨F.erase i, ⟨erase_subset_erase i hFU, ?_⟩, insert_erase hiF⟩
      rw [card_erase_of_mem hiF, hc, Nat.add_sub_cancel]
    · rintro ⟨G, ⟨hGU, hc⟩, rfl⟩
      have hiG : i ∉ G := fun h => (mem_erase.1 (hGU h)).1 rfl
      exact ⟨⟨insert_subset hi (hGU.trans (erase_subset _ _)), by
        rw [card_insert_of_notMem hiG, hc]⟩, mem_insert_self _ _⟩
  rw [himage, card_image_of_injOn, card_powersetCard, card_erase_of_mem hi]
  intro G hG H hH hGH
  have hnot : ∀ {K : Finset P}, K ∈ ((U.erase i).powersetCard k : Set (Finset P)) → i ∉ K := by
    intro K hK h
    have := (mem_powersetCard.1 hK).1 h
    exact (mem_erase.1 this).1 rfl
  have hGH' : insert i G = insert i H := hGH
  calc G = (insert i G).erase i := (erase_insert (hnot hG)).symm
    _ = (insert i H).erase i := by rw [hGH']
    _ = H := erase_insert (hnot hH)

omit [DecidableEq P] in
/-- Each element's value is counted once per subset containing it. -/
theorem sum_powersetCard_sum {U : Finset P} (k : ℕ) (r : P → ℝ) :
    ∑ F ∈ U.powersetCard (k + 1), ∑ i ∈ F, r i =
      ((U.card - 1).choose k : ℝ) * ∑ i ∈ U, r i := by
  classical
  have hswap : ∑ F ∈ U.powersetCard (k + 1), ∑ i ∈ F, r i =
      ∑ i ∈ U, ∑ F ∈ (U.powersetCard (k + 1)).filter (i ∈ ·), r i := by
    rw [Finset.sum_comm' (s' := fun i => (U.powersetCard (k + 1)).filter (i ∈ ·)) (t' := U)]
    intro F i
    simp only [mem_filter, mem_powersetCard]
    constructor
    · rintro ⟨⟨hFU, hc⟩, hiF⟩
      exact ⟨⟨⟨hFU, hc⟩, hiF⟩, hFU hiF⟩
    · rintro ⟨⟨⟨hFU, hc⟩, hiF⟩, -⟩
      exact ⟨⟨hFU, hc⟩, hiF⟩
  rw [hswap, Finset.mul_sum]
  refine Finset.sum_congr rfl fun i hi => ?_
  rw [Finset.sum_const, card_powersetCard_filter_mem hi, nsmul_eq_mul]

omit [DecidableEq P] in
/-- **The simple-random-sample expectation.** Over `(k+1)`-subsets of `U`, the expected sample mean
is `U`'s mean. -/
theorem srs_mean (U : Finset P) (k : ℕ) (hk : k + 1 ≤ U.card) (r : P → ℝ) :
    (∑ F ∈ U.powersetCard (k + 1), (∑ i ∈ F, r i) / ((k : ℝ) + 1)) /
        (U.card.choose (k + 1) : ℝ) =
      (∑ i ∈ U, r i) / (U.card : ℝ) := by
  obtain ⟨m, hm⟩ : ∃ m, U.card = m + 1 := ⟨U.card - 1, by omega⟩
  have hkm : k ≤ m := by omega
  rw [← Finset.sum_div, sum_powersetCard_sum, hm, Nat.add_sub_cancel]
  have hc : ((m + 1 : ℕ) : ℝ) * (m.choose k : ℝ) =
      ((m + 1).choose (k + 1) : ℝ) * ((k + 1 : ℕ) : ℝ) := by
    exact_mod_cast Nat.add_one_mul_choose_eq m k
  have hpos : ((m + 1).choose (k + 1) : ℝ) ≠ 0 := by
    have : 0 < (m + 1).choose (k + 1) := Nat.choose_pos (by omega)
    exact_mod_cast this.ne'
  have hk0 : ((k : ℝ) + 1) ≠ 0 := by positivity
  have hm0 : (((m + 1 : ℕ) : ℝ)) ≠ 0 := by positivity
  push_cast at hc ⊢
  push_cast at hm0
  field_simp
  linear_combination (∑ i ∈ U, r i) * hc

omit [DecidableEq P] in
/-- The expected value of a constant over simple random samples is that constant. -/
theorem srs_const (U : Finset P) (k : ℕ) (hk : k ≤ U.card) (c : ℝ) :
    (∑ _F ∈ U.powersetCard k, c) / (U.card.choose k : ℝ) = c := by
  have hpos : (U.card.choose k : ℝ) ≠ 0 := by
    have : 0 < U.card.choose k := Nat.choose_pos hk
    exact_mod_cast this.ne'
  rw [Finset.sum_const, card_powersetCard, nsmul_eq_mul]
  field_simp

end Pool

/-! ### GE-1 -/

section CrossFit

variable {P : Type*} [Fintype P] [DecidableEq P]

/-- **GE-1, the algebra.** The bias of a fold is the training share times the gap between the fit's
mean error off and on its training points. -/
theorem crossfit_gap (m f : P → ℝ) (A : Finset P) (hA : A.Nonempty) (hAc : Aᶜ.Nonempty) :
    (∑ i, f i) / (Fintype.card P : ℝ) + (∑ i ∈ Aᶜ, (m i - f i)) / (Aᶜ.card : ℝ)
        - (∑ i, m i) / (Fintype.card P : ℝ) =
      ((A.card : ℝ) / (Fintype.card P : ℝ)) *
        ((∑ i ∈ Aᶜ, (m i - f i)) / (Aᶜ.card : ℝ) - (∑ i ∈ A, (m i - f i)) / (A.card : ℝ)) := by
  have hN : (Fintype.card P : ℝ) = (A.card : ℝ) + (Aᶜ.card : ℝ) := by
    exact_mod_cast (Finset.card_add_card_compl A).symm
  have hsplit : ∑ i, (m i - f i) = ∑ i ∈ A, (m i - f i) + ∑ i ∈ Aᶜ, (m i - f i) :=
    (Finset.sum_add_sum_compl A _).symm
  have hdiff : ∑ i, m i - ∑ i, f i = ∑ i, (m i - f i) := (Finset.sum_sub_distrib _ _).symm
  have ha : (A.card : ℝ) ≠ 0 := by exact_mod_cast hA.card_pos.ne'
  have hc : (Aᶜ.card : ℝ) ≠ 0 := by exact_mod_cast hAc.card_pos.ne'
  have hNpos : (Fintype.card P : ℝ) ≠ 0 := by
    rw [hN]
    have h₁ : (0 : ℝ) < A.card := by exact_mod_cast hA.card_pos
    have h₂ : (0 : ℝ) < Aᶜ.card := by exact_mod_cast hAc.card_pos
    positivity
  set Ra := ∑ i ∈ A, (m i - f i)
  set Rc := ∑ i ∈ Aᶜ, (m i - f i)
  have hfm : ∑ i, f i = ∑ i, m i - (Ra + Rc) := by rw [← hsplit, ← hdiff]; ring
  rw [hfm]
  field_simp
  rw [hN]
  ring

/-- **GE-1.** Given a training set `A` and a fit on it, with the fold's runs a simple random sample
of size `k + 1` of the pool outside `A`, the estimator's expectation misses the pool's mean by the
training share times the gap. -/
theorem crossfit_bias (m f : P → ℝ) (A : Finset P) (hA : A.Nonempty) (k : ℕ)
    (hk : k + 1 ≤ Aᶜ.card) :
    (∑ F ∈ Aᶜ.powersetCard (k + 1),
        ((∑ i, f i) / (Fintype.card P : ℝ) + (∑ i ∈ F, (m i - f i)) / ((k : ℝ) + 1)
          - (∑ i, m i) / (Fintype.card P : ℝ))) / (Aᶜ.card.choose (k + 1) : ℝ) =
      ((A.card : ℝ) / (Fintype.card P : ℝ)) *
        ((∑ i ∈ Aᶜ, (m i - f i)) / (Aᶜ.card : ℝ) - (∑ i ∈ A, (m i - f i)) / (A.card : ℝ)) := by
  have hAc : Aᶜ.Nonempty := Finset.card_pos.1 (by omega)
  rw [← crossfit_gap m f A hA hAc]
  have hpos : (Aᶜ.card.choose (k + 1) : ℝ) ≠ 0 := by
    have : 0 < Aᶜ.card.choose (k + 1) := Nat.choose_pos hk
    exact_mod_cast this.ne'
  have hmean := srs_mean Aᶜ k hk (fun i => m i - f i)
  have hconst := srs_const Aᶜ (k + 1) hk
    ((∑ i, f i) / (Fintype.card P : ℝ) - (∑ i, m i) / (Fintype.card P : ℝ))
  rw [← hmean]
  have hsum : ∑ F ∈ Aᶜ.powersetCard (k + 1),
      ((∑ i, f i) / (Fintype.card P : ℝ) + (∑ i ∈ F, (m i - f i)) / ((k : ℝ) + 1)
        - (∑ i, m i) / (Fintype.card P : ℝ)) =
      ∑ _F ∈ Aᶜ.powersetCard (k + 1),
        ((∑ i, f i) / (Fintype.card P : ℝ) - (∑ i, m i) / (Fintype.card P : ℝ)) +
      ∑ F ∈ Aᶜ.powersetCard (k + 1), (∑ i ∈ F, (m i - f i)) / ((k : ℝ) + 1) := by
    rw [← Finset.sum_add_distrib]
    exact Finset.sum_congr rfl fun _ _ => by ring
  rw [hsum, add_div, hconst]
  ring

omit [DecidableEq P] in
/-- With the fit made apart from the pool, the runs a simple random sample of it, the estimator is
exactly unbiased for the pool's mean, whatever `f` is. -/
theorem srs_difference_unbiased (m f : P → ℝ) (k : ℕ) (hk : k + 1 ≤ Fintype.card P) :
    (∑ F ∈ (univ : Finset P).powersetCard (k + 1),
        ((∑ i, f i) / (Fintype.card P : ℝ) + (∑ i ∈ F, (m i - f i)) / ((k : ℝ) + 1))) /
        ((Fintype.card P).choose (k + 1) : ℝ) =
      (∑ i, m i) / (Fintype.card P : ℝ) := by
  have hk' : k + 1 ≤ (univ : Finset P).card := by rwa [Finset.card_univ]
  have hmean := srs_mean (univ : Finset P) k hk' (fun i => m i - f i)
  have hconst := srs_const (univ : Finset P) (k + 1) hk' ((∑ i, f i) / (Fintype.card P : ℝ))
  rw [Finset.card_univ] at hmean hconst
  rw [Finset.sum_add_distrib, add_div, hconst, hmean, Finset.sum_sub_distrib, sub_div]
  ring

end CrossFit

/-! ### Independent draws: GE-2 and telescoping -/

section Draws

variable {X : Type*} [Fintype X]

/-- The integral of `g` under a probability vector `p`. -/
noncomputable def integralOf (p g : X → ℝ) : ℝ := ∑ y, p y * g y

/-- The expectation of `G` over `n` independent draws from `p`. -/
noncomputable def drawsExpect (p : X → ℝ) (n : ℕ) (G : (Fin n → X) → ℝ) : ℝ :=
  ∑ x : Fin n → X, (∏ i, p (x i)) * G x

/-- `n` independent draws have total weight 1. -/
theorem drawsExpect_const (p : X → ℝ) (hp : ∑ y, p y = 1) (n : ℕ) (c : ℝ) :
    drawsExpect p n (fun _ => c) = c := by
  have h : ∑ x : Fin n → X, ∏ i, p (x i) = ∏ _i : Fin n, ∑ y, p y := by
    rw [Finset.prod_univ_sum, Fintype.piFinset_univ]
  rw [drawsExpect, ← Finset.sum_mul, h, hp, Finset.prod_const_one, one_mul]

/-- One draw's expectation is the integral: the other draws sum out. -/
theorem drawsExpect_coord (p : X → ℝ) (hp : ∑ y, p y = 1) (n : ℕ) (j : Fin n) (g : X → ℝ) :
    drawsExpect p n (fun x => g (x j)) = integralOf p g := by
  have h := Finset.prod_univ_sum (t := fun _ : Fin n => (univ : Finset X))
    (f := fun i y => p y * if i = j then g y else 1)
  rw [Fintype.piFinset_univ] at h
  have hlhs : ∏ i : Fin n, ∑ y, p y * (if i = j then g y else 1) = integralOf p g := by
    rw [Finset.prod_eq_single j]
    · simp [integralOf]
    · intro i _ hij
      simp [hij, hp]
    · simp
  rw [← hlhs, h, drawsExpect]
  refine Finset.sum_congr rfl fun x _ => ?_
  rw [Finset.prod_mul_distrib]
  simp [Finset.prod_ite_eq']

theorem drawsExpect_add (p : X → ℝ) (n : ℕ) (G H : (Fin n → X) → ℝ) :
    drawsExpect p n (fun x => G x + H x) = drawsExpect p n G + drawsExpect p n H := by
  simp only [drawsExpect, mul_add, Finset.sum_add_distrib]

/-- The mean of the draws indexed by `J` has the integral as its expectation. -/
theorem drawsExpect_mean (p : X → ℝ) (hp : ∑ y, p y = 1) (n : ℕ) (J : Finset (Fin n))
    (hJ : J.Nonempty) (g : X → ℝ) :
    drawsExpect p n (fun x => (∑ i ∈ J, g (x i)) / (J.card : ℝ)) = integralOf p g := by
  have hc : (J.card : ℝ) ≠ 0 := by exact_mod_cast hJ.card_pos.ne'
  have h : drawsExpect p n (fun x => ∑ i ∈ J, g (x i)) = (J.card : ℝ) * integralOf p g := by
    simp only [drawsExpect, Finset.mul_sum]
    rw [Finset.sum_comm]
    have : ∀ i ∈ J, ∑ x : Fin n → X, (∏ k, p (x k)) * g (x i) = integralOf p g :=
      fun i _ => drawsExpect_coord p hp n i g
    rw [Finset.sum_congr rfl this, Finset.sum_const, nsmul_eq_mul]
  have h' : drawsExpect p n (fun x => (∑ i ∈ J, g (x i)) / (J.card : ℝ)) =
      drawsExpect p n (fun x => ∑ i ∈ J, g (x i)) / (J.card : ℝ) := by
    simp only [drawsExpect, Finset.sum_div, mul_div_assoc]
  rw [h', h]
  field_simp

theorem drawsExpect_sum (p : X → ℝ) (n : ℕ) {κ : Type*} (s : Finset κ)
    (G : κ → (Fin n → X) → ℝ) :
    drawsExpect p n (fun x => ∑ a ∈ s, G a x) = ∑ a ∈ s, drawsExpect p n (G a) := by
  simp only [drawsExpect, Finset.mul_sum]
  exact Finset.sum_comm

/-- Two distinct draws are independent. -/
theorem drawsExpect_pair (p : X → ℝ) (hp : ∑ y, p y = 1) (n : ℕ) {i j : Fin n} (hij : i ≠ j)
    (g h : X → ℝ) :
    drawsExpect p n (fun x => g (x i) * h (x j)) = integralOf p g * integralOf p h := by
  have H := Finset.prod_univ_sum (t := fun _ : Fin n => (univ : Finset X))
    (f := fun l y => p y * if l = i then g y else if l = j then h y else 1)
  rw [Fintype.piFinset_univ] at H
  have hlhs : ∏ l : Fin n, ∑ y, p y * (if l = i then g y else if l = j then h y else 1) =
      integralOf p g * integralOf p h := by
    rw [Finset.prod_eq_mul i j hij]
    · simp [integralOf, hij.symm]
    · intro l _ hl
      simp [hl.1, hl.2, hp]
    · simp
    · simp
  rw [← hlhs, H, drawsExpect]
  refine Finset.sum_congr rfl fun x _ => ?_
  rw [Finset.prod_mul_distrib]
  congr 1
  rw [Finset.prod_eq_mul i j hij]
  · simp [hij.symm]
  · intro l _ hl
    simp [hl.1, hl.2]
  · simp
  · simp

/-- The second moment of a sum of `n` draws of a mean-zero function. -/
theorem drawsExpect_sum_sq (p g : X → ℝ) (hp : ∑ y, p y = 1) (hg : integralOf p g = 0) (n : ℕ) :
    drawsExpect p n (fun x => (∑ i, g (x i)) ^ 2) = (n : ℝ) * integralOf p (fun y => g y ^ 2) := by
  have hexp : (fun x : Fin n → X => (∑ i, g (x i)) ^ 2) =
      fun x => ∑ i, ∑ j, g (x i) * g (x j) := by
    funext x
    rw [sq, Finset.sum_mul_sum]
  rw [hexp, drawsExpect_sum]
  have hinner : ∀ i : Fin n, drawsExpect p n (fun x => ∑ j, g (x i) * g (x j)) =
      integralOf p (fun y => g y ^ 2) := by
    intro i
    rw [drawsExpect_sum, Finset.sum_eq_single i]
    · have := drawsExpect_coord p hp n i (fun y => g y ^ 2)
      simpa only [sq] using this
    · intro j _ hj
      rw [drawsExpect_pair p hp n (Ne.symm hj), hg, zero_mul]
    · simp
  simp only [hinner]
  rw [Finset.sum_const, Finset.card_univ, Fintype.card_fin, nsmul_eq_mul]

/-- Some outcome of the draws does at least as well as the expectation. -/
theorem exists_le_drawsExpect (p : X → ℝ) (h₀ : ∀ y, 0 ≤ p y) (hp : ∑ y, p y = 1) (n : ℕ)
    (G : (Fin n → X) → ℝ) : ∃ x, G x ≤ drawsExpect p n G := by
  by_contra hne
  simp only [not_exists, not_le] at hne
  have hw : ∀ x : Fin n → X, 0 ≤ ∏ i, p (x i) := fun x => Finset.prod_nonneg fun i _ => h₀ _
  have htot : ∑ x : Fin n → X, ∏ i, p (x i) = 1 := by
    have h := drawsExpect_const p hp n 1
    simpa [drawsExpect] using h
  have hs : ∑ x : Fin n → X, (∏ i, p (x i)) * (G x - drawsExpect p n G) = 0 := by
    simp only [mul_sub]
    rw [Finset.sum_sub_distrib, ← Finset.sum_mul, htot, one_mul]
    simp [drawsExpect]
  have hz := (Finset.sum_eq_zero_iff_of_nonneg
    fun x _ => mul_nonneg (hw x) (sub_pos.2 (hne x)).le).1 hs
  have hall : ∑ x : Fin n → X, ∏ i, p (x i) = 0 :=
    Finset.sum_eq_zero fun x hx =>
      (mul_eq_zero.1 (hz x hx)).resolve_right (sub_pos.2 (hne x)).ne'
  rw [hall] at htot
  exact zero_ne_one htot

/-- **GE-2.** With `a` independent training draws, `k + 1` further independent draws for the fold,
and the first term the integral of the fit, the fold is exactly unbiased for `m`'s integral,
whatever
map `fit` turns training draws into a prediction. -/
theorem guarded_fold_unbiased (p m : X → ℝ) (hp : ∑ y, p y = 1) (a k : ℕ)
    (fit : (Fin a → X) → X → ℝ) :
    drawsExpect p a (fun tr => drawsExpect p (k + 1) (fun x =>
        integralOf p (fit tr) + (∑ i, (m (x i) - fit tr (x i))) / ((k + 1 : ℕ) : ℝ))) =
      integralOf p m := by
  have hfold : ∀ tr, drawsExpect p (k + 1) (fun x =>
      integralOf p (fit tr) + (∑ i, (m (x i) - fit tr (x i))) / ((k + 1 : ℕ) : ℝ)) =
      integralOf p m := by
    intro tr
    have hJ : (univ : Finset (Fin (k + 1))).Nonempty := Finset.univ_nonempty
    have hmean := drawsExpect_mean p hp (k + 1) univ hJ (fun y => m y - fit tr y)
    rw [Finset.card_univ, Fintype.card_fin] at hmean
    rw [drawsExpect_add, drawsExpect_const p hp, hmean]
    simp only [integralOf, mul_sub, Finset.sum_sub_distrib]
    ring
  simp only [hfold]
  exact drawsExpect_const p hp a _

/-- **GE-3, telescoping on independent draws.** For fits `f`, `Y₁` made apart from the draws, and
any nonempty index set `J` of them,
`∫f + mean_n (Y₁ − f) + mean_J (Y₂ − Y₁)` has expectation `∫Y₂`. -/
theorem telescope_iid (p f Y₁ Y₂ : X → ℝ) (hp : ∑ y, p y = 1) (n : ℕ) (hn : 0 < n)
    (J : Finset (Fin n)) (hJ : J.Nonempty) :
    drawsExpect p n (fun x =>
        (integralOf p f + (∑ i, (Y₁ (x i) - f (x i))) / (n : ℝ)) +
          (∑ i ∈ J, (Y₂ (x i) - Y₁ (x i))) / (J.card : ℝ)) = integralOf p Y₂ := by
  have huniv : (univ : Finset (Fin n)).Nonempty := ⟨⟨0, hn⟩, mem_univ _⟩
  have h₁ := drawsExpect_mean p hp n univ huniv (fun y => Y₁ y - f y)
  have h₂ := drawsExpect_mean p hp n J hJ (fun y => Y₂ y - Y₁ y)
  rw [Finset.card_univ, Fintype.card_fin] at h₁
  rw [drawsExpect_add, drawsExpect_add, drawsExpect_const p hp, h₁, h₂]
  simp only [integralOf, mul_sub, Finset.sum_sub_distrib]
  ring

end Draws

/-! ### GE-3 on a pool -/

section Telescope

variable {P : Type*} [Fintype P]

/-- **GE-3, telescoping on a pool.** The runs a simple random sample of the pool, the costly level's
runs a simple random sample of those: `mean_pool f + mean_F (Y₁ − f) + mean_G (Y₂ − Y₁)` has the
pool's mean of `Y₂` as its expectation, for fits made apart from the samples. -/
theorem telescope_srs (f Y₁ Y₂ : P → ℝ) (n k : ℕ) (hn : n + 1 ≤ Fintype.card P) (hk : k ≤ n) :
    (∑ F ∈ (univ : Finset P).powersetCard (n + 1),
        ((∑ G ∈ F.powersetCard (k + 1),
            ((∑ i, f i) / (Fintype.card P : ℝ) + (∑ i ∈ F, (Y₁ i - f i)) / ((n : ℝ) + 1) +
              (∑ i ∈ G, (Y₂ i - Y₁ i)) / ((k : ℝ) + 1))) /
          ((n + 1).choose (k + 1) : ℝ))) / ((Fintype.card P).choose (n + 1) : ℝ) =
      (∑ i, Y₂ i) / (Fintype.card P : ℝ) := by
  have hn' : n + 1 ≤ (univ : Finset P).card := by rwa [Finset.card_univ]
  have hinner : ∀ F ∈ (univ : Finset P).powersetCard (n + 1),
      (∑ G ∈ F.powersetCard (k + 1),
          ((∑ i, f i) / (Fintype.card P : ℝ) + (∑ i ∈ F, (Y₁ i - f i)) / ((n : ℝ) + 1) +
            (∑ i ∈ G, (Y₂ i - Y₁ i)) / ((k : ℝ) + 1))) / ((n + 1).choose (k + 1) : ℝ) =
        (∑ i, f i) / (Fintype.card P : ℝ) + (∑ i ∈ F, (Y₁ i - f i)) / ((n : ℝ) + 1) +
          (∑ i ∈ F, (Y₂ i - Y₁ i)) / ((n : ℝ) + 1) := by
    intro F hF
    have hcard : F.card = n + 1 := (mem_powersetCard.1 hF).2
    have hk' : k + 1 ≤ F.card := by omega
    have hmean := srs_mean F k hk' (fun i => Y₂ i - Y₁ i)
    have hconst := srs_const F (k + 1) hk'
      ((∑ i, f i) / (Fintype.card P : ℝ) + (∑ i ∈ F, (Y₁ i - f i)) / ((n : ℝ) + 1))
    rw [hcard] at hmean hconst
    rw [Finset.sum_add_distrib, add_div, hconst, hmean]
    push_cast
    ring
  rw [Finset.sum_congr rfl hinner]
  have hsample : ∀ F ∈ (univ : Finset P).powersetCard (n + 1),
      (∑ i, f i) / (Fintype.card P : ℝ) + (∑ i ∈ F, (Y₁ i - f i)) / ((n : ℝ) + 1) +
          (∑ i ∈ F, (Y₂ i - Y₁ i)) / ((n : ℝ) + 1) =
        (∑ i, f i) / (Fintype.card P : ℝ) + (∑ i ∈ F, (Y₂ i - f i)) / ((n : ℝ) + 1) := by
    intro F _
    have : ∑ i ∈ F, (Y₂ i - f i) = ∑ i ∈ F, (Y₁ i - f i) + ∑ i ∈ F, (Y₂ i - Y₁ i) := by
      rw [← Finset.sum_add_distrib]
      exact Finset.sum_congr rfl fun _ _ => by ring
    rw [this, add_div]
    ring
  rw [Finset.sum_congr rfl hsample]
  exact srs_difference_unbiased Y₂ f n hn

end Telescope

/-! ### GE-3: the tuned weight -/

section Tuned

/-- The stated quadratic: the variance of `λ·mean_pool f̂ + mean_runs (y − λf̂)` for a fixed `f̂`
and a
pool independent of the runs, with `v = Var f̂`, `c = Cov(y, f̂)`, `w = Var y`. -/
noncomputable def tunedVariance (v c w Np nr l : ℝ) : ℝ :=
  l ^ 2 * v / Np + (w - 2 * l * c + l ^ 2 * v) / nr

/-- The tuned weight. -/
noncomputable def tunedWeight (v c Np nr : ℝ) : ℝ := c / ((1 + nr / Np) * v)

/-- **GE-3.** The quadratic's excess over its value at `λ*` is a square, so `λ*` is its minimizer.
-/
theorem tunedVariance_sub_tuned {v c w Np nr : ℝ} (hv : v ≠ 0) (hNp : 0 < Np) (hnr : 0 < nr)
    (l : ℝ) :
    tunedVariance v c w Np nr l - tunedVariance v c w Np nr (tunedWeight v c Np nr) =
      v * (1 / Np + 1 / nr) * (l - tunedWeight v c Np nr) ^ 2 := by
  have hsum : Np + nr ≠ 0 := by positivity
  unfold tunedVariance tunedWeight
  field_simp
  ring

/-- `λ*` does at least as well as any weight, for `Var f̂ > 0`. -/
theorem tunedVariance_tuned_le {v c w Np nr : ℝ} (hv : 0 < v) (hNp : 0 < Np) (hnr : 0 < nr)
    (l : ℝ) :
    tunedVariance v c w Np nr (tunedWeight v c Np nr) ≤ tunedVariance v c w Np nr l := by
  have h := tunedVariance_sub_tuned (w := w) (c := c) hv.ne' hNp hnr l
  have hnn : 0 ≤ v * (1 / Np + 1 / nr) * (l - tunedWeight v c Np nr) ^ 2 := by positivity
  linarith

/-- At `λ = 0` the quadratic is the plain mean's variance. -/
theorem tunedVariance_zero (v c w Np nr : ℝ) : tunedVariance v c w Np nr 0 = w / nr := by
  simp [tunedVariance]

/-- `λ = 1`, the difference estimator, is worse than the plain mean exactly when the emulator's
covariance with the runs falls short of half its variance times `1 + n_r/N_p`. -/
theorem tunedVariance_one_gt_zero_iff {v c w Np nr : ℝ} (hNp : 0 < Np) (hnr : 0 < nr) :
    tunedVariance v c w Np nr 0 < tunedVariance v c w Np nr 1 ↔ 2 * c < (1 + nr / Np) * v := by
  have key : tunedVariance v c w Np nr 1 - tunedVariance v c w Np nr 0 =
      ((1 + nr / Np) * v - 2 * c) / nr := by
    unfold tunedVariance
    field_simp
    ring
  constructor
  · intro h
    have : 0 < ((1 + nr / Np) * v - 2 * c) / nr := by rw [← key]; linarith
    have := (div_pos_iff_of_pos_right hnr).1 this
    linarith
  · intro h
    have : 0 < ((1 + nr / Np) * v - 2 * c) / nr := div_pos (by linarith) hnr
    linarith

end Tuned

/-! ### GE-3: known inclusion probabilities -/

section Inclusion

variable {P : Type*} [Fintype P] [DecidableEq P]

/-- The probability of the sample `S` when each configuration is run independently with
probability `q t`. -/
noncomputable def inclusionWeight (q : P → ℝ) (S : Finset P) : ℝ :=
  (∏ t ∈ S, q t) * ∏ t ∈ Sᶜ, (1 - q t)

/-- The sample weights have marginal `q t` at `t`. -/
theorem sum_inclusionWeight_of_mem (q : P → ℝ) (t : P) :
    ∑ S : Finset P, (if t ∈ S then inclusionWeight q S else 0) = q t := by
  have h := Finset.prod_add q (fun i => if i = t then 0 else 1 - q i) (univ : Finset P)
  rw [Finset.powerset_univ] at h
  have hprod : ∏ i, (q i + if i = t then 0 else 1 - q i) = q t := by
    rw [Finset.prod_eq_single t]
    · simp
    · intro i _ hit
      simp [hit]
    · simp
  rw [hprod] at h
  rw [h]
  refine Finset.sum_congr rfl fun S _ => ?_
  by_cases hS : t ∈ S
  · rw [if_pos hS, inclusionWeight, Finset.compl_eq_univ_sdiff]
    congr 1
    refine Finset.prod_congr rfl fun i hi => ?_
    have : i ≠ t := by
      rintro rfl
      exact (Finset.mem_sdiff.1 hi).2 hS
    simp [this]
  · rw [if_neg hS]
    symm
    apply mul_eq_zero_of_right
    apply Finset.prod_eq_zero (i := t)
    · simp [hS]
    · simp

/-- The sample weights sum to 1. -/
theorem sum_inclusionWeight (q : P → ℝ) : ∑ S : Finset P, inclusionWeight q S = 1 := by
  have h := Finset.prod_add q (fun i => 1 - q i) (univ : Finset P)
  rw [Finset.powerset_univ] at h
  simp only [add_sub_cancel, Finset.prod_const_one] at h
  simpa [inclusionWeight, Finset.compl_eq_univ_sdiff] using h.symm

/-- **GE-3, Horvitz–Thompson.** Each configuration run independently with known `q t ≠ 0`, the fit
made apart: `mean_pool f + (1/N) Σ_{t run} (m t − f t) / q t` has the pool's mean as its
expectation.
Only `q t ≠ 0` is used in the algebra; `q` being probabilities is what makes the weights a law. -/
theorem horvitz_thompson_unbiased (q m f : P → ℝ) (hq : ∀ t, q t ≠ 0) (N : ℝ) :
    ∑ S : Finset P, inclusionWeight q S *
        ((∑ i, f i) / N + (∑ t ∈ S, (m t - f t) / q t) / N) = (∑ i, m i) / N := by
  have hsum : ∑ S : Finset P, inclusionWeight q S * ∑ t ∈ S, (m t - f t) / q t =
      ∑ t, (m t - f t) := by
    have hexpand : ∀ S : Finset P, inclusionWeight q S * ∑ t ∈ S, (m t - f t) / q t =
        ∑ t, (m t - f t) / q t * (if t ∈ S then inclusionWeight q S else 0) := by
      intro S
      rw [Finset.mul_sum]
      have : ∑ t, (m t - f t) / q t * (if t ∈ S then inclusionWeight q S else 0) =
          ∑ t ∈ S, (m t - f t) / q t * inclusionWeight q S := by
        rw [← Finset.sum_filter_add_sum_filter_not univ (· ∈ S)]
        have h₀ : ∑ t ∈ univ.filter (¬ · ∈ S),
            (m t - f t) / q t * (if t ∈ S then inclusionWeight q S else 0) = 0 :=
          Finset.sum_eq_zero fun t ht => by simp [(Finset.mem_filter.1 ht).2]
        rw [h₀, add_zero, Finset.filter_mem_eq_inter, Finset.univ_inter]
        exact Finset.sum_congr rfl fun t ht => by rw [if_pos ht]
      rw [this]
      exact Finset.sum_congr rfl fun _ _ => mul_comm _ _
    simp only [hexpand]
    rw [Finset.sum_comm]
    refine Finset.sum_congr rfl fun t _ => ?_
    rw [← Finset.mul_sum, sum_inclusionWeight_of_mem, div_mul_cancel₀ _ (hq t)]
  have hsplit : ∀ S : Finset P, inclusionWeight q S *
      ((∑ i, f i) / N + (∑ t ∈ S, (m t - f t) / q t) / N) =
      inclusionWeight q S * ((∑ i, f i) / N) +
        (inclusionWeight q S * ∑ t ∈ S, (m t - f t) / q t) / N := by
    intro S
    ring
  simp only [hsplit]
  rw [Finset.sum_add_distrib, ← Finset.sum_mul, sum_inclusionWeight, ← Finset.sum_div, hsum,
    Finset.sum_sub_distrib]
  ring

end Inclusion

end BASS
