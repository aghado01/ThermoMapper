import Enthymemata.TimeFrequency.TwoSampleRecovery
import Mathlib.Algebra.BigOperators.Fin

/-! # Finite block reconstruction for stationary binary sources -/

namespace TimeFrequency.ContextTree

open Finset

abbrev Word (n : ℕ) := Fin n → Bool

def sign (b : Bool) : ℤ := if b then 1 else -1
def flip {n : ℕ} (x : Word n) : Word n := fun i ↦ !(x i)
def reverse {n : ℕ} (x : Word n) : Word n := fun i ↦ x i.rev
def slice {n : ℕ} (x : Word n) (o w : ℕ) (h : o + w ≤ n) : Word w :=
  fun i ↦ x ⟨o + i.val, by omega⟩

def dc {n : ℕ} (x : Word n) : ℤ := (∑ i, sign (x i)) ^ 2

noncomputable def mean {X : Type*} [Fintype X] (p : X → ℝ) (f : X → ℤ) : ℝ :=
  ∑ x, p x * (f x : ℝ)

noncomputable def bin {X Y : Type*} [Fintype X] [DecidableEq Y]
    (p : X → ℝ) (f : X → Y) (y : Y) : ℝ := ∑ x, if f x = y then p x else 0

theorem bin_mean {X Y : Type*} [Fintype X] [DecidableEq Y] (p : X → ℝ) (f : X → Y) (g : Y → ℤ) :
    ∑ y ∈ univ.image f, bin p f y * (g y : ℝ) = mean p (g ∘ f) := by
  classical
  simp only [bin, sum_mul, mean]
  rw [sum_comm]
  apply sum_congr rfl
  intro x _
  rw [sum_eq_single (f x)]
  · simp
  · intro y _ hy
    simp [Ne.symm hy]
  · simp

/-- Finite cylinder laws of a stationary sign-symmetric source. Translation consistency
is stated against integer tests, which include every word indicator. -/
structure Source where
  law : (n : ℕ) → Word n → ℝ
  nonneg : ∀ n x, 0 ≤ law n x
  total : ∀ n, ∑ x, law n x = 1
  translation : ∀ n o w (h : o + w ≤ n) (f : Word w → ℤ),
    mean (law n) (fun x ↦ f (slice x o w h)) = mean (law w) f
  symmetry : ∀ n x, law n (flip x) = law n x

def SameDC (D : ℕ) (P Q : Source) : Prop :=
  ∀ w, 1 ≤ w → w ≤ D + 1 → ∀ v, bin (P.law w) dc v = bin (Q.law w) dc v

def Agree {n : ℕ} (P Q : Source) (f : Word n → ℤ) := mean (P.law n) f = mean (Q.law n) f

private theorem agree_const (P Q : Source) (n : ℕ) (c : ℤ) : Agree P Q (fun _ : Word n ↦ c) := by
  simp [Agree, mean, ← sum_mul, P.total, Q.total]

private theorem agree_linear {n : ℕ} {P Q : Source} {f g : Word n → ℤ}
    (a b : ℤ) (hf : Agree P Q f) (hg : Agree P Q g) :
    Agree P Q (fun x ↦ a * f x + b * g x) := by
  have lin (p : Word n → ℝ) : mean p (fun x ↦ a * f x + b * g x) =
      a * mean p f + b * mean p g := by
    simp only [mean, Int.cast_add, Int.cast_mul, mul_add, sum_add_distrib, mul_sum]
    congr 1 <;> apply sum_congr rfl <;> intros <;> ring
  exact (lin _).trans ((congrArg₂ (fun u v : ℝ ↦ a * u + b * v) hf hg).trans (lin _).symm)

private theorem agree_sum {n : ℕ} {I : Type*} [Fintype I] {P Q : Source} (f : I → Word n → ℤ)
    (h : ∀ i, Agree P Q (f i)) : Agree P Q (fun x ↦ ∑ i, f i x) := by
  have hs (p : Word n → ℝ) : mean p (fun x ↦ ∑ i, f i x) = ∑ i, mean p (f i) := by
    simp only [mean, Int.cast_sum, mul_sum]
    rw [sum_comm]
  exact (hs _).trans ((sum_congr rfl (fun i _ ↦ h i)).trans (hs _).symm)

private theorem agree_dc {D : ℕ} {P Q : Source} (h : SameDC D P Q) {w : ℕ} (hw : w ≤ D + 1)
    (g : ℤ → ℤ) : Agree P Q (fun x : Word w ↦ g (dc x)) := by
  cases w with
  | zero =>
    have he : (fun x : Word 0 ↦ g (dc x)) = fun _ ↦ g 0 := by funext x; simp [dc]
    rw [he]
    exact agree_const _ _ _ _
  | succ w =>
    change mean _ (g ∘ dc) = mean _ (g ∘ dc)
    rw [← bin_mean, ← bin_mean]
    exact sum_congr rfl (fun v _ ↦ by rw [h _ (by omega) hw v])

/-- Squared power admits an exact second difference on any one block. -/
theorem dc_second_difference (k : ℕ) (x : Word (k + 2)) :
    dc x - dc (slice x 0 (k + 1) (by omega)) - dc (slice x 1 (k + 1) (by omega)) +
      dc (slice x 1 k (by omega)) = 2 * sign (x 0) * sign (x (Fin.last (k + 1))) := by
  let s := ∑ i : Fin k, sign (x i.succ.castSucc)
  have ha : (∑ i, sign (x i)) = sign (x 0) + s + sign (x (Fin.last (k + 1))) := by
    rw [Fin.sum_univ_castSucc, Fin.sum_univ_succ]; rfl
  have hb : (∑ i, sign (slice x 0 (k + 1) (by omega) i)) = sign (x 0) + s := by
    rw [Fin.sum_univ_succ]
    simp [slice, s]
    apply sum_congr rfl; intro i _; rfl
  have hc : (∑ i, sign (slice x 1 (k + 1) (by omega) i)) = s + sign (x (Fin.last (k + 1))) := by
    rw [Fin.sum_univ_castSucc]
    congr 1
    · apply sum_congr rfl; intro i _; unfold slice; apply congrArg (sign ∘ x); apply Fin.ext; simp; omega
    · unfold slice; apply congrArg (sign ∘ x); apply Fin.ext; simp; omega
  have hd : (∑ i, sign (slice x 1 k (by omega) i)) = s := by
    apply sum_congr rfl; intro i _; unfold slice; apply congrArg (sign ∘ x); apply Fin.ext; simp; omega
  simp only [dc, ha, hb, hc, hd]
  ring

noncomputable def powerMean (P : Source) (w : ℕ) := mean (P.law w) dc
noncomputable def lag (P : Source) (k : ℕ) :=
  mean (P.law (k + 1)) (fun x ↦ sign (x 0) * sign (x (Fin.last k)))

theorem mean_second_difference (P : Source) (k : ℕ) :
    powerMean P (k + 2) - 2 * powerMean P (k + 1) + powerMean P k = 2 * lag P (k + 1) := by
  have h := congrArg (mean (P.law (k + 2))) (funext (dc_second_difference k))
  have lin {n : ℕ} (p : Word n → ℝ) (f g : Word n → ℤ) :
      mean p (fun x ↦ f x - g x) = mean p f - mean p g := by
    simp [mean, mul_sub, sum_sub_distrib]
  have add {n : ℕ} (p : Word n → ℝ) (f g : Word n → ℤ) :
      mean p (fun x ↦ f x + g x) = mean p f + mean p g := by
    simp [mean, mul_add, sum_add_distrib]
  have dbl {n : ℕ} (p : Word n → ℝ) (f : Word n → ℤ) :
      mean p (fun x ↦ 2 * f x) = 2 * mean p f := by
    simp only [mean, Int.cast_mul, Int.cast_ofNat, mul_sum]
    apply sum_congr rfl; intros; ring
  simp only [mul_assoc, add, lin, dbl] at h
  rw [P.translation _ 0 (k + 1) _ dc, P.translation _ 1 (k + 1) _ dc, P.translation _ 1 k _ dc] at h
  unfold powerMean lag
  linarith

theorem mean_pair (P : Source) {n : ℕ} (i j : Fin n) (hij : i < j) :
    mean (P.law n) (fun x ↦ sign (x i) * sign (x j)) = lag P (j.val - i.val) := by
  let d := j.val - i.val
  have hlen : i.val + (d + 1) ≤ n := by dsimp [d]; omega
  let f := fun y : Word (d + 1) ↦ sign (y 0) * sign (y (Fin.last d))
  have he : (fun x : Word n ↦ sign (x i) * sign (x j)) = fun x ↦ f (slice x i.val (d + 1) hlen) := by
    funext x
    have h := show i.val + d = j.val by dsimp [d]; omega
    simp [f, slice, h]
  rw [he, P.translation]
  rfl

private theorem agree_pair {D n : ℕ} {P Q : Source} (h : SameDC D P Q)
    (hn : n ≤ D + 1) (i j : Fin n) (hij : i < j) :
    Agree P Q (fun x : Word n ↦ sign (x i) * sign (x j)) := by
  unfold Agree
  rw [mean_pair P i j hij, mean_pair Q i j hij]
  let k := j.val - i.val - 1
  have hlag : j.val - i.val = k + 1 := by dsimp [k]; omega
  have hw : k + 2 ≤ n := by dsimp [k]; omega
  have H (w : ℕ) (hw : w ≤ D + 1) : powerMean P w = powerMean Q w := agree_dc h hw id
  have hp := mean_second_difference P k
  have hq := mean_second_difference Q k
  rw [H _ (by omega), H _ (by omega), H _ (by omega)] at hp
  rw [hlag]
  linarith

def pairSum {n : ℕ} (t x : Word n) : ℤ :=
  ∑ i, ∑ j, if i < j then (sign (t i) * sign (t j)) * (sign (x i) * sign (x j)) else 0

def character {n : ℕ} (x : Word n) : ℤ := ∏ i, sign (x i)
def hit {X : Type*} [DecidableEq X] (a b : X) : ℤ := if a = b then 1 else 0

private theorem agree_pairs {D n : ℕ} {P Q : Source} (h : SameDC D P Q)
    (hn : n ≤ D + 1) (t : Word n) : Agree P Q (pairSum t) := by
  apply agree_sum
  intro i
  apply agree_sum
  intro j
  by_cases hij : i < j
  · simp only [if_pos hij]
    simpa only [zero_mul, add_zero] using agree_linear (sign (t i) * sign (t j)) 0
      (agree_pair h hn i j hij) (agree_const P Q n 0)
  · simpa only [if_neg hij] using agree_const P Q n 0

/-- The top DC bin supplies the fourth moment, in addition to the six pair terms. -/
theorem fourth_identity (x : Word 4) :
    8 * hit (dc x) 16 = 1 + pairSum (fun _ ↦ true) x + character x := by decide +revert

private theorem agree_fourth {D : ℕ} {P Q : Source} (h : SameDC D P Q) (hn : 4 ≤ D + 1) :
    Agree P Q (@character 4) := by
  have he : @character 4 = fun x ↦ (8 * hit (dc x) 16 - pairSum (fun _ ↦ true) x) + (-1) * 1 := by
    funext x; linarith [fourth_identity x]
  rw [he]
  have h₁ := agree_dc h hn (fun v ↦ hit v 16)
  have h₂ := agree_pairs h hn (fun _ ↦ true)
  have h₃ := agree_linear 8 (-1) h₁ h₂
  simpa only [one_mul, neg_one_mul, ← sub_eq_add_neg] using
    agree_linear 1 (-1) h₃ (agree_const P Q 4 1)

/-- Even Walsh inversion through four signs, with the complementary atoms paired. -/
theorem small_inversion {n : ℕ} (hn : 1 ≤ n) (h4 : n ≤ 4) (t x : Word n) :
    2 ^ (n - 1) * (hit x t + hit x (flip t)) =
      1 + pairSum t x + (if n = 4 then character t * character x else 0) := by
  obtain rfl | rfl | rfl | rfl : n = 1 ∨ n = 2 ∨ n = 3 ∨ n = 4 := by omega
  all_goals decide +revert

private theorem mean_hit {n : ℕ} (p : Word n → ℝ) (t : Word n) :
    mean p (fun x ↦ hit x t) = p t := by
  classical
  simp [mean, hit, mul_ite]

theorem small_block_identified {D n : ℕ} {P Q : Source} (h : SameDC D P Q)
    (hn : 1 ≤ n) (hw : n ≤ D + 1) (h4 : n ≤ 4) : P.law n = Q.law n := by
  funext t
  have hpair := agree_pairs h hw t
  have hquart : Agree P Q (fun x : Word n ↦ if n = 4 then character t * character x else 0) := by
    by_cases he : n = 4
    · subst n
      simpa only [if_true, zero_mul, add_zero] using
        agree_linear (character t) 0 (agree_fourth h hw) (agree_const P Q 4 0)
    · simpa only [if_neg he] using agree_const P Q n 0
  have hall := agree_linear 1 1 (agree_linear 1 1 (agree_const P Q n 1) hpair) hquart
  simp only [one_mul] at hall
  unfold Agree at hall
  simp_rw [← small_inversion hn h4 t] at hall
  have ht (p : Source) :
      mean (p.law n) (fun x ↦ 2 ^ (n - 1) * (hit x t + hit x (flip t))) =
        ((2: ℝ) ^ (n - 1)) * (p.law n t + p.law n t) := by
    calc _ = (2: ℝ) ^ (n - 1) * (mean (p.law n) (fun x ↦ hit x t) +
          mean (p.law n) (fun x ↦ hit x (flip t))) := by
            simp only [mean, Int.cast_mul, Int.cast_pow, Int.cast_ofNat, Int.cast_add, ← sum_add_distrib, mul_sum]
            apply sum_congr rfl; intros; ring
      _ = _ := by rw [mean_hit, mean_hit, p.symmetry]
  rw [ht, ht] at hall
  have hp : (0: ℝ) < 2 ^ (n - 1) := by positivity
  nlinarith

noncomputable def nextKernel (D : ℕ) (p : Word (D + 1) → ℝ) (h : Word D) (b : Bool) :=
  p (Fin.snoc h b) / (p (Fin.snoc h false) + p (Fin.snoc h true))

/-- A suffix is a leaf when its next-sign law is constant over all full histories
that end in it, and no proper suffix has that property. -/
def minimalTree {D : ℕ} (K : Word D → Bool → ℝ) : Set (List Bool) :=
  let constant := fun s : List Bool ↦ ∀ h h' : Word D, s <:+ List.ofFn h → s <:+ List.ofFn h' → K h = K h'
  {s | (∃ h : Word D, s <:+ List.ofFn h) ∧ constant s ∧
    ∀ u : List Bool, u <:+ s → u.length < s.length → ¬ constant u}

theorem low_depth_recovery {D : ℕ} (hD : D ≤ 3) {P Q : Source} (h : SameDC D P Q) :
    nextKernel D (P.law (D + 1)) = nextKernel D (Q.law (D + 1)) ∧
    minimalTree (nextKernel D (P.law (D + 1))) = minimalTree (nextKernel D (Q.law (D + 1))) := by
  rw [small_block_identified h (by omega : 1 ≤ D + 1) le_rfl (by omega)]
  exact ⟨rfl, rfl⟩

theorem kernel_denominator_pos {D : ℕ} {p : Word (D + 1) → ℝ} (hp : ∀ x, 0 < p x) (h : Word D) :
    0 < p (Fin.snoc h false) + p (Fin.snoc h true) := add_pos (hp _) (hp _)

private theorem mean_add {X : Type*} [Fintype X] (p : X → ℝ) (f g : X → ℤ) :
    mean p (fun x ↦ f x + g x) = mean p f + mean p g := by
  simp [mean, mul_add, sum_add_distrib]

private theorem mean_scale {X : Type*} [Fintype X] (p : X → ℝ) (a : ℤ) (f : X → ℤ) :
    mean p (fun x ↦ a * f x) = a * mean p f := by
  simp only [mean, Int.cast_mul, mul_sum]
  apply sum_congr rfl; intros; ring

theorem power_mean_formula (P : Source) (w : ℕ) :
    powerMean P w = w + 2 * ∑ i ∈ range w, ((w : ℝ) - (i + 1)) * lag P (i + 1) := by
  let F := fun n : ℕ ↦ (n : ℝ) + 2 * ((n : ℝ) * (∑ i ∈ range n, lag P (i + 1)) -
    ∑ i ∈ range n, (i + 1 : ℝ) * lag P (i + 1))
  have hstep (k : ℕ) : F (k + 2) - 2 * F (k + 1) + F k = 2 * lag P (k + 1) := by
    simp only [F, sum_range_succ, Nat.cast_succ]
    ring
  have hm (n : ℕ) : powerMean P n = F n := by
    induction n using Nat.twoStepInduction with
    | zero => simp [powerMean, mean, dc, F]
    | one =>
      have he : (@dc 1) = fun _ ↦ 1 := by funext x; decide +revert
      simp [powerMean, he, mean, F, P.total]
    | more k h₁ h₂ =>
      have h := mean_second_difference P k
      rw [h₁, h₂] at h
      linarith [hstep k]
  rw [hm]
  simp only [F, mul_sum, ← sum_sub_distrib, sub_mul]

theorem fourth_from_dc (P : Source) : 8 * bin (P.law 4) dc 16 =
    1 + 3 * lag P 1 + 2 * lag P 2 + lag P 3 + mean (P.law 4) character := by
  let p := fun (i j : Fin 4) (x : Word 4) ↦ sign (x i) * sign (x j)
  have he : (fun x : Word 4 ↦ 8 * hit (dc x) 16) = fun x ↦
      1 + p 0 1 x + p 0 2 x + p 0 3 x + p 1 2 x + p 1 3 x + p 2 3 x + character x := by
    funext x
    change 8 * hit (dc x) 16 = _
    dsimp [p]
    decide +revert
  have h := congrArg (mean (P.law 4)) he
  have htop : mean (P.law 4) (fun x ↦ hit (dc x) 16) = bin (P.law 4) dc 16 := by
    simp [mean, hit, bin, mul_ite]
  have hone : mean (P.law 4) (fun _ ↦ 1) = 1 := by simp [mean, P.total]
  have hp (i j : Fin 4) (hij : i < j) : mean (P.law 4) (p i j) = lag P (j.val - i.val) :=
    mean_pair P i j hij
  simp only [mean_add, mean_scale, htop, hone, hp 0 1 (by decide), hp 0 2 (by decide),
    hp 0 3 (by decide), hp 1 2 (by decide), hp 1 3 (by decide), hp 2 3 (by decide)] at h
  norm_num at h
  linarith

theorem odd_mean_zero (P : Source) {n : ℕ} (f : Word n → ℤ) (hf : ∀ x, f (flip x) = -f x) :
    mean (P.law n) f = 0 := by
  have hflip : Function.Involutive (@flip n) := by intro x; funext i; simp [flip]
  have h := Equiv.sum_comp hflip.toPerm (fun x ↦ P.law n x * (f x : ℝ))
  change (∑ x, P.law n (flip x) * (f (flip x) : ℝ)) = _ at h
  simp only [P.symmetry, hf, Int.cast_neg, mul_neg, sum_neg_distrib] at h
  unfold mean
  linarith

theorem prefix_mass (P : Source) (D : ℕ) (h : Word D) :
    P.law (D + 1) (Fin.snoc h false) + P.law (D + 1) (Fin.snoc h true) = P.law D h := by
  let g := fun x : Word (D + 1) ↦ P.law (D + 1) x * ((hit (slice x 0 D (by omega)) h : ℤ) : ℝ)
  have he (w : Word D) (b : Bool) : slice (Fin.snoc w b : Word (D + 1)) 0 D (by omega) = w := by
    funext i
    have hi : (⟨0 + i.val, by omega⟩ : Fin (D + 1)) = i.castSucc := Fin.ext (by simp)
    simp only [slice, hi, Fin.snoc_castSucc]
  have hm := P.translation (D + 1) 0 D (by omega) (fun w ↦ hit w h)
  rw [mean_hit] at hm
  change (∑ x, g x) = _ at hm
  rw [← (Fin.snocEquiv (fun _ : Fin (D + 1) ↦ Bool)).sum_comp g, Fintype.sum_prod_type] at hm
  change (∑ b : Bool, ∑ w : Word D, P.law _ (Fin.snoc w b) *
    ((hit (slice (Fin.snoc w b : Word (D + 1)) 0 D (by omega)) h : ℤ) : ℝ)) = _ at hm
  simp only [he, Fintype.sum_bool, hit, Int.cast_ite, Int.cast_one, Int.cast_zero, mul_ite, mul_one, mul_zero,
    sum_ite_eq', mem_univ, if_true] at hm
  linarith

theorem constant_kernel_tree {D : ℕ} (K : Word D → Bool → ℝ) (hc : ∀ h h', K h = K h') :
    minimalTree K = {[]} := by
  ext s
  constructor
  · rintro ⟨_, _, hm⟩
    by_contra hne
    have hs : s ≠ [] := hne
    have hl : ([] : List Bool).length < s.length := by cases s <;> simp_all
    exact hm [] (by simp) hl (fun h h' _ _ ↦ hc h h')
  · rintro rfl
    exact ⟨⟨fun _ ↦ true, by simp⟩, fun h h' _ _ ↦ hc h h', by simp⟩

theorem empty_history_kernel (P : Source) (hp : ∀ x, 0 < P.law 1 x) (h : Word 0) (b : Bool) :
    nextKernel 0 (P.law 1) h b = 1 / 2 := by
  have he : flip (Fin.snoc h false : Word 1) = Fin.snoc h true := by
    funext i; fin_cases i; rfl
  have hs := P.symmetry _ (Fin.snoc h false : Word 1)
  rw [he] at hs
  have hpos := hp (Fin.snoc h false)
  unfold nextKernel
  cases b <;> rw [hs] <;> field_simp <;> ring

#print axioms constant_kernel_tree
#print axioms empty_history_kernel
#print axioms power_mean_formula
#print axioms mean_pair
#print axioms fourth_from_dc
#print axioms odd_mean_zero
#print axioms prefix_mass
#print axioms bin_mean
#print axioms dc_second_difference
#print axioms mean_second_difference
#print axioms fourth_identity
#print axioms small_inversion
#print axioms small_block_identified
#print axioms kernel_denominator_pos
#print axioms low_depth_recovery

end TimeFrequency.ContextTree
