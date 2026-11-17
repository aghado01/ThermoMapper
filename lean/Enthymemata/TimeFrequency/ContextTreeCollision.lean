import Enthymemata.TimeFrequency.ContextTreeRecovery

/-! # Stationary binary sources indistinguishable by every rectangular power vector -/

namespace TimeFrequency.ContextTree

open Finset

noncomputable def uniform (n : ℕ) : ℝ := (2 ^ (n : ℕ) : ℝ)⁻¹

theorem uniform_total (n : ℕ) : ∑ _ : Word n, uniform n = 1 := by
  simp [uniform, Fintype.card_pi_const, Fintype.card_fun]

private theorem uniform_succ (n : ℕ) : 2 * uniform (n + 1) = uniform n := by
  unfold uniform
  rw [pow_succ]
  field_simp

/-- Every selection of distinct coordinates of a fair word is itself a fair word. -/
theorem uniform_restrict {n w : ℕ} (e : Fin w → Fin n) (hi : Function.Injective e)
    (g : Word w → ℝ) :
    ∑ x : Word n, uniform n * g (x ∘ e) = ∑ y : Word w, uniform w * g y := by
  classical
  let R := Set.range e
  let E := Equiv.ofInjective e hi
  let S : Word n ≃ Word w × ({i : Fin n // i ∉ R} → Bool) :=
    (Equiv.piEquivPiSubtypeProd (· ∈ R) (fun _ ↦ Bool)).trans
      (Equiv.prodCongr (Equiv.arrowCongr E.symm (Equiv.refl Bool)) (Equiv.refl _))
  have hS (x : Word n) : (S x).1 = x ∘ e := rfl
  let N := Fintype.card ({i : Fin n // i ∉ R} → Bool)
  have hcard := Fintype.card_congr S
  have hN : (2: ℝ) ^ n = 2 ^ w * N := by
    exact_mod_cast (show 2 ^ n = 2 ^ w * N by simpa [N, Fintype.card_fun] using hcard)
  have hfac : (N : ℝ) * uniform n = uniform w := by
    unfold uniform
    field_simp
    linarith
  calc ∑ x : Word n, uniform n * g (x ∘ e) =
        ∑ z : Word w × ({i : Fin n // i ∉ R} → Bool), uniform n * g z.1 :=
          S.sum_comp (fun z ↦ uniform n * g z.1)
    _ = _ := by
      rw [Fintype.sum_prod_type]
      apply sum_congr rfl
      intro y _
      simp only [sum_const, card_univ, nsmul_eq_mul]
      change (N : ℝ) * (uniform n * g y) = _
      rw [← mul_assoc, hfac]

/-- Independent fair signs give mass 2⁻ⁿ to each cylinder of width n. -/
noncomputable def fairSource : Source where
  law n _ := uniform n
  nonneg _ _ := by unfold uniform; positivity
  total := uniform_total
  translation n o w ho f := by
    let e := fun j : Fin w ↦ (⟨o + j.val, by omega⟩ : Fin n)
    have he : Function.Injective e := by
      intro a b hab; apply Fin.ext; have := congrArg Fin.val hab; dsimp [e] at this; omega
    exact uniform_restrict e he (fun x ↦ (f x : ℝ))
  symmetry _ _ := rfl

private theorem sign_not (b : Bool) : sign (!b) = -sign b := by cases b <;> decide
private theorem sign_abs (b : Bool) : |sign b| = 1 := by cases b <;> decide

@[simp] theorem reverse_reverse {n : ℕ} (x : Word n) : reverse (reverse x) = x := by
  funext i; simp [reverse]

def revWords (n : ℕ) : Word n ≃ Word n :=
  ⟨reverse, reverse, reverse_reverse, reverse_reverse⟩

def toggle {n : ℕ} (i : Fin n) (x : Word n) : Word n :=
  fun j ↦ if j = i then !(x j) else x j

def toggleWords {n : ℕ} (i : Fin n) : Word n ≃ Word n :=
  Function.Involutive.toPerm (f := toggle i) (by intro x; funext j; by_cases h : j = i <;> simp [toggle, h])

/-- An antisymmetric weight sums to zero on each invariant observable fiber. -/
theorem fiber_cancels {X Y : Type*} [Fintype X] [DecidableEq Y]
    (T : X ≃ X) (d : X → ℝ) (f : X → Y)
    (hd : ∀ x, d (T x) = -d x) (hf : ∀ x, f (T x) = f x) (y : Y) : bin d f y = 0 := by
  have he := T.sum_comp (fun x ↦ if f x = y then d x else 0)
  simp only [hf, hd] at he
  have hneg (x : X) : (if f x = y then -d x else 0) = -(if f x = y then d x else 0) := by
    split_ifs <;> simp
  simp_rw [hneg] at he
  rw [sum_neg_distrib] at he
  unfold bin
  linarith

namespace Collision

abbrev L (m : ℕ) := m + 5
abbrev D (m : ℕ) := m + 4

def one (m : ℕ) : Fin (L m) := ⟨1, by dsimp [L] <;> omega⟩
def two (m : ℕ) : Fin (L m) := ⟨2, by dsimp [L] <;> omega⟩

def chi {m : ℕ} (x : Word (L m)) : ℤ :=
  sign (x 0) * sign (x (one m)) * sign (x (two m)) * sign (x (Fin.last (D m)))

def skew {m : ℕ} (x : Word (L m)) : ℤ := chi x - chi (reverse x)

noncomputable def weight (m : ℕ) (ε : ℝ) (x : Word (L m)) : ℝ :=
  uniform (L m) * (1 + ε * skew x)

private theorem chi_abs {m : ℕ} (x : Word (L m)) : |chi x| = 1 := by
  simp [chi, abs_mul, sign_abs]

theorem weight_positive {m : ℕ} {ε : ℝ} (hε : 0 < ε) (hu : ε < 1 / 2)
    (x : Word (L m)) : 0 < weight m ε x := by
  have h₁ := (abs_le.mp (le_of_eq (chi_abs x))).1
  have h₂ := (abs_le.mp (le_of_eq (chi_abs (reverse x)))).2
  have hl : (-2: ℝ) ≤ (skew x : ℤ) := by exact_mod_cast (show -2 ≤ skew x by unfold skew <;> omega)
  unfold weight
  apply mul_pos (by unfold uniform; positivity)
  nlinarith [mul_le_mul_of_nonneg_left hl hε.le]

theorem weight_symmetric {m : ℕ} (ε : ℝ) (x : Word (L m)) : weight m ε (flip x) = weight m ε x := by
  simp [weight, skew, chi, flip, reverse, sign_not]

private theorem skew_reverse {m : ℕ} (x : Word (L m)) : skew (reverse x) = -skew x := by
  simp [skew]

private theorem tilted_bin {Y : Type*} [DecidableEq Y] {m : ℕ} (ε : ℝ) (f : Word (L m) → Y) (y : Y) :
    bin (weight m ε) f y = bin (fun _ ↦ uniform (L m)) f y +
      (uniform (L m) * ε) * bin (fun x ↦ (skew x : ℝ)) f y := by
  simp only [bin, mul_sum, ← sum_add_distrib]
  apply sum_congr rfl
  intro x _
  by_cases h : f x = y <;> simp [h, weight] <;> ring

theorem invariant_observable {Y : Type*} [DecidableEq Y] {m : ℕ} (ε : ℝ) (f : Word (L m) → Y)
    (hf : ∀ x, f (reverse x) = f x) (y : Y) :
    bin (weight m ε) f y = bin (fun _ ↦ uniform (L m)) f y := by
  rw [tilted_bin, fiber_cancels (revWords _) (fun x ↦ (skew x : ℝ)) f (fun x ↦ by
    change ((skew (reverse x) : ℤ) : ℝ) = _; rw [skew_reverse, Int.cast_neg]) hf y]
  ring

theorem weight_total (m : ℕ) (ε : ℝ) : ∑ x, weight m ε x = 1 := by
  have he := invariant_observable ε (fun _ : Word (L m) ↦ ()) (fun _ ↦ rfl) ()
  simpa only [bin, if_true, uniform_total] using he

private theorem endpoint_toggle {m : ℕ} (x : Word (L m)) (i : Fin (L m))
    (hi : i = 0 ∨ i = Fin.last (D m)) : skew (toggle i x) = -skew x := by
  have neq (a b : Fin (L m)) (h : a.val ≠ b.val) : a ≠ b := fun e ↦ h (congrArg Fin.val e)
  have h1 : (one m) ≠ 0 := neq _ _ (by simp [one])
  have h2 : (two m) ≠ 0 := neq _ _ (by simp [two])
  have h3 : (one m).rev ≠ 0 := neq _ _ (by simp [one, Fin.rev, L] <;> omega)
  have h4 : (two m).rev ≠ 0 := neq _ _ (by simp [two, Fin.rev, L] <;> omega)
  have h5 : (one m) ≠ Fin.last (D m) := neq _ _ (by simp [one, D] <;> omega)
  have h6 : (two m) ≠ Fin.last (D m) := neq _ _ (by simp [two, D] <;> omega)
  have h7 : (one m).rev ≠ Fin.last (D m) := neq _ _ (by simp [one, Fin.rev, L, D] <;> omega)
  have h8 : (two m).rev ≠ Fin.last (D m) := neq _ _ (by simp [two, Fin.rev, L, D] <;> omega)
  have hn : (Fin.last (D m)) ≠ 0 := neq _ _ (by simp [D])
  rcases hi with rfl | rfl <;>
    simp [skew, chi, reverse, toggle, h1, h2, h3, h4, h5, h6, h7, h8, hn, Ne.symm hn, sign_not] <;> ring

theorem proper_observable {Y : Type*} [DecidableEq Y] {m w o : ℕ} (hw : w < L m)
    (ho : o + w ≤ L m) (ε : ℝ) (f : Word w → Y) (y : Y) :
    bin (weight m ε) (fun x ↦ f (slice x o w ho)) y = bin (fun _ ↦ uniform w) f y := by
  let i : Fin (L m) := if o = 0 then Fin.last (D m) else 0
  have hi : i = 0 ∨ i = Fin.last (D m) := by by_cases h : o = 0 <;> simp [i, h]
  have avoid (j : Fin w) : (⟨o + j.val, by omega⟩ : Fin (L m)) ≠ i := by
    intro he
    have h := congrArg Fin.val he
    by_cases hz : o = 0 <;> simp [i, hz, L, D] at h <;> dsimp [L] at hw ho <;> omega
  have hf (x : Word (L m)) : slice (toggle i x) o w ho = slice x o w ho := by
    funext j
    simp [slice, toggle, avoid]
  rw [tilted_bin, fiber_cancels (toggleWords i) (fun x ↦ (skew x : ℝ)) _ (fun x ↦ by
    change ((skew (toggle i x) : ℤ) : ℝ) = _; rw [endpoint_toggle _ _ hi, Int.cast_neg])
    (fun x ↦ congrArg f (hf x)) y, mul_zero, add_zero]
  let e := fun j : Fin w ↦ (⟨o + j.val, by omega⟩ : Fin (L m))
  have he : Function.Injective e := by intro a b hab; apply Fin.ext; have := congrArg Fin.val hab; dsimp [e] at this <;> omega
  unfold bin slice
  simpa only [e, Function.comp_def, mul_ite, mul_one, mul_zero] using
    uniform_restrict e he (fun x ↦ if f x = y then 1 else 0)

private theorem toggle_snoc {n : ℕ} (h : Word n) :
    toggle (Fin.last n) (Fin.snoc h false) = Fin.snoc h true := by
  funext j
  refine Fin.lastCases ?_ (fun i ↦ ?_) j <;> simp [toggle, Fin.castSucc_ne_last]

private theorem toggle_cons {n : ℕ} (h : Word n) :
    toggle 0 (Fin.cons false h) = Fin.cons true h := by
  funext j
  refine Fin.cases ?_ (fun i ↦ ?_) j <;> simp [toggle]

theorem outgoing_mass (m : ℕ) (ε : ℝ) (h : Word (D m)) :
    weight m ε (Fin.snoc h false) + weight m ε (Fin.snoc h true) = uniform (D m) := by
  rw [weight, weight, ← toggle_snoc h, endpoint_toggle _ _ (Or.inr rfl), Int.cast_neg]
  have hu := uniform_succ (D m)
  change 2 * uniform (L m) = _ at hu
  linarith

theorem incoming_mass (m : ℕ) (ε : ℝ) (h : Word (D m)) :
    weight m ε (Fin.cons false h) + weight m ε (Fin.cons true h) = uniform (D m) := by
  rw [weight, weight, ← toggle_cons h, endpoint_toggle _ _ (Or.inl rfl), Int.cast_neg]
  have hu := uniform_succ (D m)
  change 2 * uniform (L m) = _ at hu
  linarith

noncomputable def transition (m : ℕ) (ε : ℝ) (h : Word (D m)) (b : Bool) :=
  weight m ε (Fin.snoc h b) / uniform (D m)

theorem transition_stochastic (m : ℕ) {ε : ℝ} (hε : 0 < ε) (hu : ε < 1 / 2) (h : Word (D m)) :
    (∀ b, 0 < transition m ε h b) ∧ ∑ b, transition m ε h b = 1 := by
  have hpos : 0 < uniform (D m) := by unfold uniform; positivity
  refine ⟨fun _ ↦ div_pos (weight_positive hε hu _) hpos, ?_⟩
  simp [transition, Fintype.sum_bool, ← add_div, add_comm, outgoing_mass, hpos.ne']

theorem transition_sign_flip (m : ℕ) (ε : ℝ) (h : Word (D m)) (b : Bool) :
    transition m ε (flip h) (!b) = transition m ε h b := by
  have he : Fin.snoc (flip h) (!b) = flip (Fin.snoc h b : Word (L m)) := by
    funext i; refine Fin.lastCases ?_ (fun j ↦ ?_) i <;> simp [flip]
  simp only [transition, he, weight_symmetric]

/-- Products of the next-sign probabilities assign the same mass to a path and its sign reversal. -/
noncomputable def pathLikelihood {n : ℕ} (K : Word n → Bool → ℝ) (h : Word n) : List Bool → ℝ
  | [] => 1
  | b :: bs => K h b * pathLikelihood K (Fin.tail (Fin.snoc h b : Word (n + 1))) bs

theorem path_sign_flip (m : ℕ) (ε : ℝ) (h : Word (D m)) (bs : List Bool) :
    pathLikelihood (transition m ε) (flip h) (bs.map Bool.not) =
      pathLikelihood (transition m ε) h bs := by
  induction bs generalizing h with
  | nil => rfl
  | cons b bs ih =>
    simp only [List.map_cons, pathLikelihood, transition_sign_flip]
    have he : Fin.snoc (flip h) (!b) = flip (Fin.snoc h b : Word (L m)) := by
      funext i; refine Fin.lastCases ?_ (fun j ↦ ?_) i <;> simp [flip]
    rw [he]
    exact congrArg (fun r ↦ transition m ε h b * r) (ih (Fin.tail (Fin.snoc h b : Word (L m))))

theorem block_flow (m : ℕ) (ε : ℝ) (h : Word (D m)) (b : Bool) :
    uniform (D m) * transition m ε h b = weight m ε (Fin.snoc h b) := by
  unfold transition uniform
  field_simp

def advance {n : ℕ} (h : Word n) (b : Bool) : Word n := Fin.tail (Fin.snoc h b : Word (n + 1))

/-- The shift-and-append pushforward of stationary history and next-sign pairs.
Reindexing the edge law by the prepended bit and destination history gives the invariant law. -/
theorem stationary_history (m : ℕ) (ε : ℝ) (f : Word (D m) → ℝ) :
    ∑ h : Word (D m), ∑ b : Bool, uniform (D m) * transition m ε h b * f (advance h b) =
      ∑ h, uniform (D m) * f h := by
  simp_rw [block_flow]
  let F := fun x : Word (L m) ↦ weight m ε x * f (Fin.tail x)
  calc _ = ∑ x : Word (L m), F x := by
          rw [sum_comm, ← (Fin.snocEquiv (fun _ : Fin (L m) ↦ Bool)).sum_comp F, Fintype.sum_prod_type]
          rfl
    _ = _ := by
      rw [← (Fin.consEquiv (fun _ : Fin (L m) ↦ Bool)).sum_comp F, Fintype.sum_prod_type, sum_comm]
      apply sum_congr rfl
      intro h _
      change (∑ b : Bool, weight m ε (Fin.cons b h) * f (Fin.tail (Fin.cons b h : Word (L m)))) = _
      simp only [Fin.tail_cons, ← sum_mul, Fintype.sum_bool]
      rw [add_comm, incoming_mass]

def marked {n : ℕ} (oldest : Bool) (i : Fin n) : Bool :=
  if i.val = 0 then oldest else decide (i.val ≠ 1)

private theorem snoc_marked (m : ℕ) (a : Bool) :
    Fin.snoc (marked a : Word (D m)) true = (marked a : Word (L m)) := by
  funext i
  refine Fin.lastCases ?_ (fun j ↦ ?_) i <;> simp [marked, D]

private theorem skew_marked (m : ℕ) (a : Bool) :
    skew (marked a : Word (L m)) = -2 * sign a := by
  have h₁ : (one m).rev.val = m + 3 := by simp [one, Fin.rev, L] <;> omega
  have h₂ : (two m).rev.val = m + 2 := by simp [two, Fin.rev, L] <;> omega
  simp [skew, chi, reverse, marked, h₁, h₂, one, two, D, L, sign]
  split <;> ring

theorem exact_memory (m : ℕ) {ε : ℝ} (hε : 0 < ε) :
    ∃ h h' : Word (D m), Fin.tail h = Fin.tail h' ∧ transition m ε h true ≠ transition m ε h' true := by
  refine ⟨marked false, marked true, ?_, ?_⟩
  · funext i; simp [Fin.tail, marked]
  · have hu : 0 < uniform (L m) := by unfold uniform; positivity
    have hd : uniform (D m) ≠ 0 := by unfold uniform; positivity
    simp only [transition, snoc_marked, weight, skew_marked, Int.cast_mul]
    intro he
    have he' := (div_left_inj' hd).mp he
    have he'' := mul_left_cancel₀ hu.ne' he'
    norm_num [sign] at he''
    linarith

theorem fair_kernel (m : ℕ) (h : Word (D m)) (b : Bool) :
    nextKernel (D m) (fun _ ↦ uniform (L m)) h b = 1 / 2 := by
  unfold nextKernel uniform
  field_simp
  <;> ring

theorem fair_tree (m : ℕ) : minimalTree (nextKernel (D m) (fun _ ↦ uniform (L m))) = {[]} := by
  apply constant_kernel_tree
  intro h h'
  funext b
  rw [fair_kernel, fair_kernel]

theorem different_trees (m : ℕ) {ε : ℝ} (hε : 0 < ε) :
    [] ∈ minimalTree (nextKernel (D m) (fun _ ↦ uniform (L m))) ∧
    [] ∉ minimalTree (transition m ε) := by
  constructor
  · rw [fair_tree]; rfl
  · rintro ⟨_, hk, _⟩
    obtain ⟨h, h', _, hne⟩ := exact_memory m hε
    exact hne (congrFun (hk h h' (by simp) (by simp)) true)

/-- The four coordinates of the asymmetric DC taper, in time order. -/
def selected (m : ℕ) (i : Fin 4) : Fin (L m) :=
  if i.val = 0 then 0 else if i.val = 1 then one m else if i.val = 2 then two m else Fin.last (D m)

theorem selected_injective (m : ℕ) : Function.Injective (selected m) := by
  intro i j he
  have hv := congrArg Fin.val he
  fin_cases i <;> fin_cases j <;> (first | rfl | (exfalso; norm_num [selected, one, two, D] at hv <;> omega))

def asymmetric {m : ℕ} (x : Word (L m)) := dc (x ∘ selected m)

private theorem selected_chi {m : ℕ} (x : Word (L m)) : character (x ∘ selected m) = chi x := by
  simp [character, Fin.prod_univ_four, selected, chi]

private theorem chi_at_top {m : ℕ} (x : Word (L m)) (h : asymmetric x = 16) : chi x = 1 := by
  have h₄ : ∀ y : Word 4, dc y = 16 → character y = 1 := by decide
  rw [← selected_chi]
  exact h₄ _ h

theorem fair_asymmetric (m : ℕ) : bin (fun _ : Word (L m) ↦ uniform (L m)) asymmetric 16 = 1 / 8 := by
  have hu := uniform_restrict (selected m) (selected_injective m) (fun x ↦ if dc x = 16 then 1 else 0)
  have hc : (∑ x : Word 4, hit (dc x) 16) = 2 := by decide
  have hr : (∑ x : Word 4, ((hit (dc x) 16 : ℤ) : ℝ)) = 2 := by exact_mod_cast hc
  calc _ = uniform 4 * ∑ x : Word 4, ((hit (dc x) 16 : ℤ) : ℝ) := by
          rw [mul_sum]
          simpa only [bin, asymmetric, hit, Int.cast_ite, Int.cast_one, Int.cast_zero, mul_ite, mul_one, mul_zero] using hu
    _ = _ := by rw [hr]; norm_num [uniform]

private theorem reverse_character_bin (m : ℕ) :
    bin (fun x : Word (L m) ↦ (chi (reverse x) : ℝ)) asymmetric 16 = 0 := by
  let i := (one m).rev
  have hv : i.val = m + 3 := by simp [i, one, Fin.rev, L] <;> omega
  have hl : Fin.last (D m) ≠ i := by intro he; have := congrArg Fin.val he; simp [hv, D] at this
  have hz : (0 : Fin (L m)) ≠ i := by intro he; have := congrArg Fin.val he; simp [hv] at this
  have htwo : (two m).rev ≠ i := by intro he; have := congrArg Fin.val he; simp [hv, two, Fin.rev, L] at this <;> omega
  apply fiber_cancels (toggleWords i)
  · intro x
    change ((chi (reverse (toggle i x)) : ℤ) : ℝ) = _
    have he : chi (reverse (toggle i x)) = -chi (reverse x) := by
      simp [chi, reverse, toggle, hl, hz, htwo, i, sign_not] <;> ring
    rw [he, Int.cast_neg]
  · intro x
    apply congrArg dc
    funext j
    change toggle i x (selected m j) = x (selected m j)
    have hj : selected m j ≠ i := by
      intro he
      have h := congrArg Fin.val he
      fin_cases j <;> simp [selected, hv, one, two, D] at h <;> omega
    simp [toggle, hj]

theorem perturbed_asymmetric (m : ℕ) (ε : ℝ) :
    bin (weight m ε) asymmetric 16 = (1 + ε) / 8 := by
  have he : uniform (L m) * bin (fun x : Word (L m) ↦ (skew x : ℝ)) asymmetric 16 = 1 / 8 := by
    rw [← fair_asymmetric m]
    have hs : bin (fun x : Word (L m) ↦ (skew x : ℝ)) asymmetric 16 =
        bin (fun _ : Word (L m) ↦ (1: ℝ)) asymmetric 16 -
        bin (fun x : Word (L m) ↦ (chi (reverse x) : ℝ)) asymmetric 16 := by
      simp only [bin, ← sum_sub_distrib]
      apply sum_congr rfl
      intro x _
      by_cases h : asymmetric x = 16
      · simp [h, skew, chi_at_top _ h]
      · simp [h]
    rw [hs, reverse_character_bin, sub_zero]
    simp only [bin, mul_sum, mul_ite, mul_one, mul_zero]
  rw [tilted_bin, fair_asymmetric]
  linarith [show uniform (L m) * ε * bin (fun x : Word (L m) ↦ (skew x : ℝ)) asymmetric 16 =
    ε * (1 / 8) by rw [← he]; ring]

#print axioms selected_injective
#print axioms exact_memory
#print axioms different_trees
#print axioms fair_kernel
#print axioms fair_tree
#print axioms fair_asymmetric
#print axioms perturbed_asymmetric
#print axioms uniform_total
#print axioms uniform_restrict
#print axioms reverse_reverse
#print axioms fiber_cancels
#print axioms weight_positive
#print axioms weight_symmetric
#print axioms invariant_observable
#print axioms weight_total
#print axioms proper_observable
#print axioms outgoing_mass
#print axioms incoming_mass
#print axioms transition_stochastic
#print axioms fairSource
#print axioms transition_sign_flip
#print axioms path_sign_flip
#print axioms block_flow
#print axioms stationary_history

end Collision
end TimeFrequency.ContextTree
