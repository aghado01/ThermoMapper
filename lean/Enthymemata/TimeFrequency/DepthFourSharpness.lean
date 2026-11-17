import Enthymemata.TimeFrequency.DepthFourStationary

/-! # Two tapers are necessary within the four-site indicator class -/

namespace TimeFrequency.ContextTree.DepthFour

open Complex Finset Collision

def rotate (x : Word 5) : Word 5 := fun i ↦ x (i + 1)

noncomputable def root (k : ℕ) : ℂ := Complex.exp (-(2 * Real.pi * ((k : ℂ) / 5) * I))

private theorem dft_polynomial (x : Fin 5 → ℝ) (k : ℕ) :
    frameDFT x k = ∑ i : Fin 5, (x i : ℂ) * (root k) ^ i.val := by
  unfold frameDFT root
  apply sum_congr rfl
  intro i _
  rw [← Complex.exp_nat_mul]
  congr 2
  ring

private theorem fifth_root (k : ℕ) : root k ^ 5 = 1 := by
  unfold root
  rw [← Complex.exp_nat_mul]
  have he : (5 : ℂ) * -(2 * Real.pi * ((k : ℂ) / 5) * I) = -((k : ℂ) * (2 * Real.pi * I)) := by ring
  norm_num only [Nat.cast_ofNat]
  rw [he, Complex.exp_neg, Complex.exp_nat_mul_two_pi_mul_I, inv_one]

private theorem root_unit (k : ℕ) : Complex.normSq (root k) = 1 := by
  have he : root k = Complex.exp (((-(2 * Real.pi * ((k : ℝ) / 5)) : ℝ) : ℂ) * I) := by
    unfold root
    push_cast
    congr 1
    ring
  rw [he, Complex.normSq_eq_norm_sq, Complex.norm_exp_ofReal_mul_I]
  norm_num

theorem rotate_power (x : Fin 5 → ℝ) (k : ℕ) :
    framePower (fun i ↦ x (i + 1)) k = framePower x k := by
  have he : root k * frameDFT (fun i ↦ x (i + 1)) k = frameDFT x k := by
    simp only [dft_polynomial, Fin.sum_univ_succ]
    norm_num
    change root k * ((x 1 : ℂ) + ((x 2 : ℂ) * root k + ((x 3 : ℂ) * root k ^ 2 +
      ((x 4 : ℂ) * root k ^ 3 + (x 0 : ℂ) * root k ^ 4)))) =
      (x 0 : ℂ) + ((x 1 : ℂ) * root k + ((x 2 : ℂ) * root k ^ 2 + ((x 3 : ℂ) * root k ^ 3 +
      (x 4 : ℂ) * root k ^ 4)))
    calc _ = (x 0 : ℂ) * root k ^ 5 + (x 1 : ℂ) * root k + (x 2 : ℂ) * root k ^ 2 +
        (x 3 : ℂ) * root k ^ 3 + (x 4 : ℂ) * root k ^ 4 := by ring
      _ = _ := by rw [fifth_root]; ring
  have hn := congrArg Complex.normSq he
  simpa only [map_mul, root_unit, one_mul, framePower] using hn

theorem rotate_vector (x : Word 5) : powerVector (rotate x) = powerVector x := by
  funext k
  exact rotate_power (fun i ↦ (sign (x i) : ℝ)) k.val

def rotateWords : Word 5 ≃ Word 5 where
  toFun := rotate
  invFun := fun x i ↦ x (i + 4)
  left_inv x := by funext i; decide +revert
  right_inv x := by funext i; decide +revert

private theorem rotated_quartic (r : Fin 5) (x : Word 5) :
    quartic r (rotate x) = quartic (r + 1) x := by decide +revert

private theorem cyclic_fiber {Y : Type*} [DecidableEq Y] (f : Word 5 → Y)
    (hf : ∀ x, f (rotate x) = f x) (v : Y) (r : Fin 5) :
    bin (fun x ↦ (quartic (r + 1) x : ℝ)) f v = bin (fun x ↦ (quartic r x : ℝ)) f v := by
  have he := Equiv.sum_comp rotateWords (fun x ↦ if f x = v then (quartic r x : ℝ) else 0)
  change (∑ x, if f (rotate x) = v then (quartic r (rotate x) : ℝ) else 0) = _ at he
  simpa only [hf, rotated_quartic, bin] using he

theorem full_vector_law (a b : ℝ) (v : Fin 5 → ℝ) :
    bin (weight a b) powerVector v = bin (fun _ ↦ (1 / 32 : ℝ)) powerVector v := by
  have h1 := cyclic_fiber powerVector rotate_vector v 1
  have h2 := cyclic_fiber powerVector rotate_vector v 2
  change bin (fun x ↦ (quartic 2 x : ℝ)) powerVector v = _ at h1
  change bin (fun x ↦ (quartic 3 x : ℝ)) powerVector v = _ at h2
  rw [law_decomposition, h2, h1]
  ring

theorem all_fourier_laws (a b : ℝ) {w o : ℕ} (ho : o + w ≤ 5) (v : Fin w → ℝ) :
    bin (weight a b) (fun x ↦ powerVector (slice x o w ho)) v =
      bin (fun _ ↦ uniform w) powerVector v := by
  by_cases hw : w = 5
  · subst w
    have hz : o = 0 := by omega
    subst o
    have he (x : Word 5) : slice x 0 5 ho = x := by funext i; simp [slice]
    simp_rw [he]
    convert full_vector_law a b v using 1
    norm_num [uniform]
  · exact proper_laws (by omega) ho a b powerVector v

/-- Replacing the omitted sign preserves a taper and cancels every other quartic. -/
private theorem taper_character_bin (r s : Fin 5) (v : ℤ) (h : r ≠ s) :
    bin (fun x : Word 5 ↦ (quartic s x : ℝ)) (taper r) v = 0 := by
  apply fiber_cancels (toggleWords r)
  · intro x
    have he : quartic s (toggle r x) = -quartic s x := by decide +revert
    change ((quartic s (toggle r x) : ℤ) : ℝ) = _
    rw [he, Int.cast_neg]
  · intro x
    change taper r (toggle r x) = taper r x
    decide +revert

noncomputable def coefficient (a b : ℝ) (r : Fin 5) :=
  if r = 3 then a else if r = 2 then b else if r = 1 then -(a + b) else 0

theorem taper_law (a b : ℝ) (r : Fin 5) (v : ℤ) :
    bin (weight a b) (taper r) v = bin (fun _ ↦ (1 / 32 : ℝ)) (taper r) v +
      coefficient a b r /32 * bin (fun x ↦ (quartic r x : ℝ)) (taper r) v := by
  rw [law_decomposition]
  fin_cases r <;> simp [coefficient, taper_character_bin, Fin.ext_iff] <;> ring

theorem taper_maximum (a b : ℝ) (r : Fin 5) :
    bin (weight a b) (taper r) 16 = (1 + coefficient a b r) / 8 := by
  have hcount : ∑ x : Word 5, hit (taper r x) 16 = 4 := by decide +revert
  have hsum : (∑ x : Word 5, if taper r x = 16 then quartic r x else 0) = 4 := by decide +revert
  have hc : ∑ x : Word 5, ((hit (taper r x) 16 : ℤ) : ℝ) = 4 := by exact_mod_cast hcount
  have hs : (∑ x : Word 5, if taper r x = 16 then (quartic r x : ℝ) else 0) = 4 := by
    exact_mod_cast hsum
  rw [taper_law]
  have hu : bin (fun _ : Word 5 ↦ (1 / 32 : ℝ)) (taper r) 16 = 1 / 8 := by
    have he : bin (fun _ : Word 5 ↦ (1 / 32 : ℝ)) (taper r) 16 =
        (1 / 32 : ℝ)*∑ x : Word 5, ((hit (taper r x) 16 : ℤ) : ℝ) := by
      simp only [bin, hit, mul_sum, Int.cast_ite, Int.cast_one, Int.cast_zero,
        mul_ite, mul_one, mul_zero]
    rw [he, hc]; norm_num
  rw [hu, show bin (fun x ↦ (quartic r x : ℝ)) (taper r) 16 = 4 from hs]
  ring

theorem two_parameters_recovered {a b c d : ℝ}
    (hA : bin (weight a b) (taper 3) 16 = bin (weight c d) (taper 3) 16)
    (hB : bin (weight a b) (taper 2) 16 = bin (weight c d) (taper 2) 16) :
    a = c ∧ b = d := by
  simp [taper_maximum, coefficient] at hA hB
  constructor <;> linarith

def invisibleA (r : Fin 5) : ℝ := if r = 3 then 0 else 1
def invisibleB (r : Fin 5) : ℝ := if r = 3 then 1 else if r = 1 then -1 else 0

theorem admissible_neg {a b : ℝ} (h : Admissible a b) : Admissible (-a) (-b) := by
  simpa only [Admissible, ← neg_add, abs_neg] using h

theorem coefficient_neg (a b : ℝ) (r : Fin 5) :
    coefficient (-a) (-b) r = -coefficient a b r := by
  fin_cases r <;> simp [coefficient]
  ring

theorem one_taper_direction (r : Fin 5) {ε : ℝ} (hp : 0 < ε) (hu : ε < 1 / 2) :
    Admissible (ε*invisibleA r) (ε*invisibleB r) ∧
    (ε*invisibleA r ≠ 0 ∨ ε*invisibleB r ≠ 0) ∧
    coefficient (ε*invisibleA r) (ε*invisibleB r) r = 0 := by
  fin_cases r <;> simp [invisibleA, invisibleB, coefficient, Admissible, abs_of_pos hp,
    abs_neg, hp.ne'] <;> linarith

theorem one_taper_insufficient (r : Fin 5) {ε : ℝ} (hp : 0 < ε) (hu : ε < 1 / 2) :
    let a := ε*invisibleA r
    let b := ε*invisibleB r
    (∀ x, 0 < weight a b x) ∧
    (∀ v, bin (weight a b) (taper r) v = bin (fun _ ↦ (1 / 32 : ℝ)) (taper r) v) ∧
    (∃ h h' : Word 4, Fin.tail h = Fin.tail h' ∧ kernel a b h true ≠ kernel a b h' true) ∧
    (∀ w o (ho : o + w ≤ 5) (v : Fin w → ℝ),
      bin (weight a b) (fun x ↦ powerVector (slice x o w ho)) v =
        bin (fun _ ↦ uniform w) powerVector v) := by
  obtain ⟨ha, hn, hc⟩ := one_taper_direction r hp hu
  refine ⟨positive ha, ?_, exact_four hn, fun _ _ ho ↦ all_fourier_laws _ _ ho⟩
  intro v
  rw [taper_law, hc, zero_div, zero_mul, add_zero]

theorem four_site_shape (S : Finset (Fin 5)) (hS : S.card = 4) :
    ∃ r : Fin 5, S = univ.erase r := by
  have hc : Sᶜ.card = 1 := by rw [card_compl, hS]; decide
  obtain ⟨r, hr⟩ := card_eq_one.mp hc
  refine ⟨r, ?_⟩
  apply compl_injective
  rw [hr]
  ext i
  simp

theorem fair_tree : minimalTree (kernel 0 0) = {[]} := by
  apply constant_kernel_tree
  intro h h'
  funext t
  simp [kernel, weight]

theorem collision_tree {a b : ℝ} (hn : a ≠ 0 ∨ b ≠ 0) :
    [] ∉ minimalTree (kernel a b) := by
  rintro ⟨_, h, _⟩
  obtain ⟨s, t, _, hst⟩ := exact_four hn
  exact hst (congrFun (h s t (by simp) (by simp)) true)

#print axioms four_site_shape
#print axioms fair_tree
#print axioms collision_tree
#print axioms rotate_power
#print axioms rotate_vector
#print axioms full_vector_law
#print axioms all_fourier_laws
#print axioms taper_law
#print axioms taper_maximum
#print axioms two_parameters_recovered
#print axioms admissible_neg
#print axioms coefficient_neg
#print axioms one_taper_direction
#print axioms one_taper_insufficient

end TimeFrequency.ContextTree.DepthFour
