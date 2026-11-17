import Enthymemata.TimeFrequency.DepthFourCompletion
import Mathlib.Algebra.BigOperators.Field

/-! # Positive stationary witnesses for the depth-four measurement design -/

namespace TimeFrequency.ContextTree.DepthFour

open Finset Collision

noncomputable def weight (a b : ℝ) (x : Word 5) :=
  (1 + a * (quartic 3 x : ℝ) + b * (quartic 2 x : ℝ) - (a + b) * (quartic 1 x : ℝ)) / 32

def Admissible (a b : ℝ) := |a| + |b| + |a + b| < 1

private theorem quartic_abs (r : Fin 5) (x : Word 5) : |quartic r x| = 1 := by decide +revert

theorem positive {a b : ℝ} (h : Admissible a b) (x : Word 5) : 0 < weight a b x := by
  have bound (c : ℝ) (r : Fin 5) : |c * (quartic r x : ℝ)| = |c| := by
    rw [abs_mul, ← Int.cast_abs, quartic_abs]; simp
  have hA := neg_abs_le (a * (quartic 3 x : ℝ))
  have hB := neg_abs_le (b * (quartic 2 x : ℝ))
  have hC := le_abs_self ((a + b) * (quartic 1 x : ℝ))
  rw [bound] at hA hB hC
  unfold weight Admissible at *
  linarith

theorem symmetry (a b : ℝ) (x : Word 5) : weight a b (flip x) = weight a b x := by
  have h (r : Fin 5) : quartic r (flip x) = quartic r x := by decide +revert
  simp [weight, h]

private theorem endpoint_sum (h : Word 4) (r : Fin 5) (hr : r ≠ 0 ∧ r ≠ 4) :
    quartic r (Fin.snoc h false) + quartic r (Fin.snoc h true) = 0 ∧
    quartic r (Fin.cons false h) + quartic r (Fin.cons true h) = 0 := by decide +revert

theorem endpoint_masses (a b : ℝ) (h : Word 4) :
    weight a b (Fin.snoc h false) + weight a b (Fin.snoc h true) = 1 / 16 ∧
    weight a b (Fin.cons false h) + weight a b (Fin.cons true h) = 1 / 16 := by
  have ha := endpoint_sum h 3 (by decide)
  have hb := endpoint_sum h 2 (by decide)
  have hc := endpoint_sum h 1 (by decide)
  have cast (u v : ℤ) (he : u + v = 0) : (u : ℝ) + (v : ℝ) = 0 := by exact_mod_cast he
  have ha₁ := cast _ _ ha.1
  have ha₂ := cast _ _ ha.2
  have hb₁ := cast _ _ hb.1
  have hb₂ := cast _ _ hb.2
  have hc₁ := cast _ _ hc.1
  have hc₂ := cast _ _ hc.2
  constructor <;> unfold weight
  · nlinarith [congrArg (fun z : ℝ ↦ a * z) ha₁, congrArg (fun z : ℝ ↦ b * z) hb₁,
      congrArg (fun z : ℝ ↦ (a + b) * z) hc₁]
  · nlinarith [congrArg (fun z : ℝ ↦ a * z) ha₂, congrArg (fun z : ℝ ↦ b * z) hb₂,
      congrArg (fun z : ℝ ↦ (a + b) * z) hc₂]

theorem total_mass (a b : ℝ) : ∑ x, weight a b x = 1 := by
  rw [← (Fin.snocEquiv (fun _ : Fin 5 ↦ Bool)).sum_comp (weight a b),
    Fintype.sum_prod_type, sum_comm]
  change (∑ h : Word 4, ∑ t : Bool, weight a b (Fin.snoc h t)) = 1
  simp only [Fintype.sum_bool]
  have he (h : Word 4) :
      weight a b (Fin.snoc h true) + weight a b (Fin.snoc h false) = 1 / 16 := by
    linarith [(endpoint_masses a b h).1]
  simp only [he, sum_const, card_univ, nsmul_eq_mul]
  norm_num [Fintype.card_fun]

noncomputable def kernel (a b : ℝ) (h : Word 4) (t : Bool) := 16 * weight a b (Fin.snoc h t)

theorem stochastic {a b : ℝ} (ha : Admissible a b) (h : Word 4) :
    (∀ t, 0 < kernel a b h t) ∧ ∑ t, kernel a b h t = 1 := by
  refine ⟨fun t ↦ mul_pos (by norm_num) (positive ha _), ?_⟩
  simp only [Fintype.sum_bool, kernel]
  linarith [(endpoint_masses a b h).1]

theorem kernel_from_block (a b : ℝ) : kernel a b = nextKernel 4 (weight a b) := by
  funext h t
  simp only [kernel, nextKernel, (endpoint_masses a b h).1]
  ring

theorem stationary (a b : ℝ) (f : Word 4 → ℝ) :
    ∑ h : Word 4, ∑ t : Bool, (1 / 16 : ℝ) * kernel a b h t * f (advance h t) =
      ∑ h, (1 / 16 : ℝ) * f h := by
  have flow (h : Word 4) (t : Bool) :
      (1 / 16 : ℝ) * kernel a b h t = weight a b (Fin.snoc h t) := by
    unfold kernel; ring
  simp_rw [flow]
  let F := fun x : Word 5 ↦ weight a b x * f (Fin.tail x)
  calc _ = ∑ x : Word 5, F x := by
        rw [sum_comm, ← (Fin.snocEquiv (fun _ : Fin 5 ↦ Bool)).sum_comp F, Fintype.sum_prod_type]
        rfl
    _ = _ := by
      rw [← (Fin.consEquiv (fun _ : Fin 5 ↦ Bool)).sum_comp F, Fintype.sum_prod_type, sum_comm]
      apply sum_congr rfl
      intro h _
      change (∑ t : Bool, weight a b (Fin.cons t h) * f h) = _
      rw [← sum_mul, Fintype.sum_bool]
      have he := (endpoint_masses a b h).2
      rw [add_comm] at he
      rw [he]

theorem sign_symmetric_kernel (a b : ℝ) (h : Word 4) (t : Bool) :
    kernel a b (flip h) (!t) = kernel a b h t := by
  have he : Fin.snoc (flip h) (!t) = flip (Fin.snoc h t : Word 5) := by
    funext i; fin_cases i <;> rfl
  simp only [kernel, he, symmetry]

private theorem endpoint_antisymmetric (x : Word 5) (r i : Fin 5)
    (hr : r ≠ 0 ∧ r ≠ 4) (hi : i = 0 ∨ i = 4) :
    quartic r (toggle i x) = -quartic r x := by decide +revert

private theorem cancel_subblock {w o : ℕ} (hw : w < 5) (ho : o + w ≤ 5)
    (r : Fin 5) (hr : r ≠ 0 ∧ r ≠ 4) {Y : Type*} [DecidableEq Y] (f : Word w → Y) (v : Y) :
    bin (fun x : Word 5 ↦ (quartic r x : ℝ)) (fun x ↦ f (slice x o w ho)) v = 0 := by
  let i : Fin 5 := if o = 0 then 4 else 0
  have hi : i = 0 ∨ i = 4 := by by_cases h : o = 0 <;> simp [i, h]
  have avoid (j : Fin w) : (⟨o + j.val, by omega⟩ : Fin 5) ≠ i := by
    intro he
    have hv := congrArg Fin.val he
    by_cases hz : o = 0 <;> simp [i, hz] at hv
    omega
  apply fiber_cancels (toggleWords i)
  · intro x
    change ((quartic r (toggle i x) : ℤ) : ℝ) = _
    rw [endpoint_antisymmetric _ _ _ hr hi, Int.cast_neg]
  · intro x
    apply congrArg f
    funext j
    simp [slice, toggleWords, Function.Involutive.toPerm, toggle, avoid]

theorem law_decomposition {Y : Type*} [DecidableEq Y] (a b : ℝ) (f : Word 5 → Y) (v : Y) :
    bin (weight a b) f v = bin (fun _ ↦ (1 / 32 : ℝ)) f v +
      (a * bin (fun x ↦ (quartic 3 x : ℝ)) f v +
       b * bin (fun x ↦ (quartic 2 x : ℝ)) f v -
       (a + b) * bin (fun x ↦ (quartic 1 x : ℝ)) f v) / 32 := by
  simp only [bin, mul_sum, ← sum_add_distrib, ← sum_sub_distrib, sum_div]
  apply sum_congr rfl
  intro x _
  by_cases h : f x = v <;> simp [h, weight]
  ring

theorem proper_laws {w o : ℕ} (hw : w < 5) (ho : o + w ≤ 5)
    {Y : Type*} [DecidableEq Y] (a b : ℝ) (f : Word w → Y) (v : Y) :
    bin (weight a b) (fun x ↦ f (slice x o w ho)) v = bin (fun _ ↦ uniform w) f v := by
  rw [law_decomposition, cancel_subblock hw ho 3 (by decide),
    cancel_subblock hw ho 2 (by decide), cancel_subblock hw ho 1 (by decide)]
  simp only [mul_zero, sub_zero, add_zero, zero_div]
  let e := fun j : Fin w ↦ (⟨o + j.val, by omega⟩ : Fin 5)
  have he : Function.Injective e := by
    intro j k h; apply Fin.ext; have := congrArg Fin.val h; dsimp [e] at this; omega
  have hu := uniform_restrict e he (fun x ↦ if f x = v then (1 : ℝ) else 0)
  unfold bin slice
  simpa only [e, Function.comp_def, mul_ite, mul_one, mul_zero, uniform, one_div,
    show (2 : ℝ) ^ 5 = 32 by norm_num] using hu

def triple (r : Fin 4) (h : Word 4) : ℤ := ∏ i ∈ univ.erase r, sign (h i)

theorem prediction_formula (a b : ℝ) (h : Word 4) :
    kernel a b h true =
      (1 + a * (triple 3 h : ℝ) + b * (triple 2 h : ℝ) - (a + b) * (triple 1 h : ℝ)) / 2 := by
  have he (r : Fin 4) : quartic r.castSucc (Fin.snoc h true) = triple r h := by decide +revert
  unfold kernel weight
  rw [show (3:Fin 5) = (3:Fin 4).castSucc from rfl, he,
      show (2:Fin 5) = (2:Fin 4).castSucc from rfl, he,
      show (1:Fin 5) = (1:Fin 4).castSucc from rfl, he]
  ring

theorem exact_four {a b : ℝ} (hn : a ≠ 0 ∨ b ≠ 0) :
    ∃ h h' : Word 4, Fin.tail h = Fin.tail h' ∧ kernel a b h true ≠ kernel a b h' true := by
  by_contra! h
  have hA := h (![false, true, false, true] : Word 4) ![true, true, false, true] rfl
  have hB := h (![false, true, true, false] : Word 4) ![true, true, true, false] rfl
  simp only [prediction_formula] at hA hB
  have ha0 : triple 3 (![false, true, false, true] : Word 4) = 1 ∧
      triple 2 (![false, true, false, true] : Word 4) = -1 ∧
      triple 1 (![false, true, false, true] : Word 4) = 1 := by decide
  have ha1 : triple 3 (![true, true, false, true] : Word 4) = -1 ∧
      triple 2 (![true, true, false, true] : Word 4) = 1 ∧
      triple 1 (![true, true, false, true] : Word 4) = -1 := by decide
  have hb0 : triple 3 (![false, true, true, false] : Word 4) = -1 ∧
      triple 2 (![false, true, true, false] : Word 4) = 1 ∧
      triple 1 (![false, true, true, false] : Word 4) = 1 := by decide
  have hb1 : triple 3 (![true, true, true, false] : Word 4) = 1 ∧
      triple 2 (![true, true, true, false] : Word 4) = -1 ∧
      triple 1 (![true, true, true, false] : Word 4) = -1 := by decide
  rw [ha0.1, ha0.2.1, ha0.2.2, ha1.1, ha1.2.1, ha1.2.2] at hA
  rw [hb0.1, hb0.2.1, hb0.2.2, hb1.1, hb1.2.1, hb1.2.2] at hB
  norm_num at hA hB
  rcases hn with hn | hn <;> apply hn <;> linarith

#print axioms positive
#print axioms symmetry
#print axioms endpoint_masses
#print axioms total_mass
#print axioms stochastic
#print axioms kernel_from_block
#print axioms stationary
#print axioms sign_symmetric_kernel
#print axioms law_decomposition
#print axioms proper_laws
#print axioms prediction_formula
#print axioms exact_four

end TimeFrequency.ContextTree.DepthFour
