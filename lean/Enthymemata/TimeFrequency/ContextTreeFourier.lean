import Enthymemata.TimeFrequency.ContextTreeCollision

/-! # The full rectangular Fourier-power vector and its asymmetric control -/

namespace TimeFrequency.ContextTree

open Complex Finset
open scoped ComplexConjugate

private noncomputable def angle (n k j : ℕ) : ℝ := -(2 * Real.pi * (k * j / n))

private theorem frame_as_exp {n : ℕ} (x : Fin n → ℝ) (k : ℕ) :
    frameDFT x k = ∑ j, (x j : ℂ) * Complex.exp ((angle n k j.val : ℂ) * I) := by
  unfold frameDFT angle
  push_cast
  apply sum_congr rfl
  intro j _
  congr 2
  ring

/-- Time reversal conjugates the DFT and multiplies it by one common unit phase per bin. -/
theorem dft_reversal {n : ℕ} (x : Fin n → ℝ) (k : ℕ) :
    frameDFT (fun j ↦ x j.rev) k =
      Complex.exp ((angle n k (n - 1) : ℂ) * I) * conj (frameDFT x k) := by
  rw [frame_as_exp, frame_as_exp, ← Equiv.sum_comp (@Fin.revPerm n)
    (fun j : Fin n ↦ (x j.rev : ℂ) * Complex.exp ((angle n k j.val : ℂ) * I))]
  simp only [Fin.revPerm_apply, Fin.rev_rev, map_sum, map_mul, conj_ofReal, ← Complex.exp_conj,
    Complex.conj_I, Finset.mul_sum]
  apply sum_congr rfl
  intro j _
  have he : j.rev.val + j.val = n - 1 := by simp [Fin.rev]; omega
  have hc : (j.rev.val : ℝ) + j.val = ((n - 1 : ℕ) : ℝ) := by exact_mod_cast he
  have ha : angle n k j.rev.val = angle n k (n - 1) - angle n k j.val := by
    unfold angle
    rw [← hc]
    ring
  rw [ha, Complex.ofReal_sub, sub_mul, sub_eq_add_neg, Complex.exp_add]
  simp only [mul_neg]
  ring

theorem rectangular_reversal {n : ℕ} (x : Fin n → ℝ) (k : ℕ) :
    framePower (fun j ↦ x j.rev) k = framePower x k := by
  simp only [framePower, dft_reversal, map_mul, Complex.normSq_conj]
  rw [Complex.normSq_eq_norm_sq, Complex.norm_exp_ofReal_mul_I]
  simp

noncomputable def powerVector {n : ℕ} (x : Word n) (k : Fin n) :=
  framePower (fun i ↦ (sign (x i) : ℝ)) k.val

theorem powerVector_reverse {n : ℕ} (x : Word n) : powerVector (reverse x) = powerVector x := by
  funext k
  exact rectangular_reversal (fun i : Fin n ↦ (sign (x i) : ℝ)) k.val

theorem dc_is_rectangular {n : ℕ} (x : Word n) :
    framePower (fun i ↦ (sign (x i) : ℝ)) 0 = (dc x : ℝ) := by
  simp [framePower_dc, dc]

namespace Collision

theorem full_power_law (m : ℕ) (ε : ℝ) (v : Fin (L m) → ℝ) :
    bin (weight m ε) powerVector v = bin (fun _ ↦ uniform (L m)) powerVector v :=
  invariant_observable ε _ powerVector_reverse _

/-- All frequencies are delivered jointly within a width. The starts may be any
contiguous subblocks. No joint experiment across widths is asserted. -/
theorem all_width_power_laws {m w o : ℕ} (ho : o + w ≤ L m) (ε : ℝ) (v : Fin w → ℝ) :
    bin (weight m ε) (fun x ↦ powerVector (slice x o w ho)) v = bin (fun _ ↦ uniform w) powerVector v := by
  by_cases hw : w = L m
  · subst w
    have hz : o = 0 := by omega
    subst o
    have he (x : Word (L m)) : slice x 0 (L m) ho = x := by funext i; simp [slice]
    simp_rw [he]
    apply full_power_law
  · exact proper_observable (by omega) ho ε powerVector v

theorem all_width_dc_laws {m w o : ℕ} (ho : o + w ≤ L m) (ε : ℝ) (v : ℤ) :
    bin (weight m ε) (fun x ↦ dc (slice x o w ho)) v = bin (fairSource.law w) dc v := by
  by_cases hw : w = L m
  · subst w
    have hz : o = 0 := by omega
    subst o
    have he (x : Word (L m)) : slice x 0 (L m) ho = x := by funext i; simp [slice]
    simp_rw [he]
    apply invariant_observable
    intro x
    apply congrArg (fun r : ℤ ↦ r ^ 2)
    exact Equiv.sum_comp (@Fin.revPerm (L m)) (fun j ↦ sign (x j))
  · exact proper_observable (by omega) ho ε dc v

theorem asymmetric_taper_dc (m : ℕ) (x : Word (L m)) :
    framePower (fun i ↦ if i ∈ univ.image (selected m) then (sign (x i) : ℝ) else 0) 0 =
      (asymmetric x : ℝ) := by
  rw [framePower_dc]
  unfold asymmetric dc
  push_cast
  congr 1
  rw [← sum_filter]
  simp only [filter_mem_eq_inter, univ_inter]
  rw [sum_image (fun a _ b _ h ↦ selected_injective m h)]
  rfl

/-- One statement combines the entire depth-D witness with the asymmetric-cell separation. -/
theorem high_depth_collision {d : ℕ} (hd : 4 ≤ d) {ε : ℝ} (hε : 0 < ε) (hu : ε < 1 / 2) :
    ∃ (p : Word (d + 1) → ℝ) (K : Word d → Bool → ℝ) (a : Word (d + 1) → ℤ),
      (∀ x, 0 < p x) ∧ (∑ x, p x) = 1 ∧ (∀ x, p (flip x) = p x) ∧
      (∀ h, (∀ b, 0 < K h b) ∧ ∑ b, K h b = 1) ∧
      (∀ h b, uniform d * K h b = p (Fin.snoc h b)) ∧
      (∀ f : Word d → ℝ, ∑ h, ∑ b, uniform d * K h b * f (advance h b) =
        ∑ h, uniform d * f h) ∧
      (∃ h h' : Word d, (∀ i : Fin d, i.val ≠ 0 → h i = h' i) ∧ K h true ≠ K h' true) ∧
      (∀ w o (ho : o + w ≤ d + 1) (v : Fin w → ℝ),
        bin p (fun x ↦ powerVector (slice x o w ho)) v = bin (fun _ ↦ uniform w) powerVector v) ∧
      bin p a 16 = (1 + ε) / 8 ∧ bin (fun _ ↦ uniform (d + 1)) a 16 = 1 / 8 := by
  obtain ⟨m, rfl⟩ : ∃ m, d = m + 4 := ⟨d - 4, by omega⟩
  obtain ⟨h, h', ht, hdiff⟩ := exact_memory m hε
  have he (i : Fin (D m)) : i.val ≠ 0 → h i = h' i := by
    refine Fin.cases (by simp) (fun j _ ↦ congrFun ht j) i
  exact ⟨weight m ε, transition m ε, asymmetric, weight_positive hε hu, weight_total _ _,
    weight_symmetric _, transition_stochastic _ hε hu, block_flow _ _, stationary_history _ _,
    ⟨h, h', he, hdiff⟩, fun _ _ h ↦ all_width_power_laws h ε, perturbed_asymmetric _ _, fair_asymmetric _⟩

#print axioms dft_reversal
#print axioms rectangular_reversal
#print axioms powerVector_reverse
#print axioms dc_is_rectangular
#print axioms full_power_law
#print axioms all_width_power_laws
#print axioms all_width_dc_laws
#print axioms asymmetric_taper_dc
#print axioms high_depth_collision

end Collision
end TimeFrequency.ContextTree
