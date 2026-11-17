import Enthymemata.TimeFrequency.MemoryDepthTesting
import Mathlib.Data.Fintype.Sets

/-! # Actual separate-view laws for the all-depth Fourier collision -/

namespace TimeFrequency.ContextTree.DepthEvidence

open Finset Collision

variable {X Y : Type*} [Fintype X] [Fintype Y]

@[ext] theorem Law.ext {P Q : Law X} (h : P.mass = Q.mass) : P = Q := by
  cases P
  cases Q
  cases h
  rfl

noncomputable def observed {Z : Type*} [DecidableEq Z] (P : Law X) (f : X → Z) :
    @Law {z // z ∈ univ.image f} (Finset.Subtype.fintype (univ.image f)) where
  mass z := bin P.mass f z.val
  positive z := by
    obtain ⟨x, _, hx⟩ := mem_image.mp z.property
    have hn (y : X) (_ : y ∈ univ) : 0 ≤ if f y = z.val then P.mass y else 0 := by
      split_ifs
      · exact (P.positive y).le
      · rfl
    have h := single_le_sum hn (mem_univ x)
    rw [if_pos hx] at h
    exact lt_of_lt_of_le (P.positive x) h
  total := by
    rw [sum_coe_sort]
    have h := bin_mean P.mass f (fun _ ↦ 1)
    simpa [mean, P.total] using h

theorem observed_eq {Z : Type*} [DecidableEq Z] (P Q : Law X) (f : X → Z)
    (h : ∀ z, bin P.mass f z = bin Q.mass f z) : observed P f = observed Q f := by
  apply Law.ext
  funext z
  exact h z.val

noncomputable def independent {I : Type*} [Fintype I] [DecidableEq I]
    {Z : I → Type*} [∀ i, Fintype (Z i)] (P : ∀ i, Law (Z i)) : Law (∀ i, Z i) where
  mass x := ∏ i, (P i).mass (x i)
  positive x := prod_pos (fun i _ ↦ (P i).positive (x i))
  total := by
    rw [← Fintype.prod_sum]
    simp only [fun i ↦ (P i).total, prod_const_one]

noncomputable def fairLaw (n : ℕ) : Law (Word n) where
  mass _ := uniform n
  positive _ := by unfold uniform; positivity
  total := uniform_total n

noncomputable def sourceLaw (m : ℕ) {a : ℝ} (ha : 0 < a) (hu : a < 1 / 2) :
    Law (Word (L m)) where
  mass := weight m a
  positive := weight_positive ha hu
  total := weight_total m a

def taperBit {m : ℕ} (x : Word (L m)) := decide (asymmetric x = 16)

theorem binary_complement (P : Law X) (f : X → Bool) :
    bin P.mass f false = 1 - bin P.mass f true := by
  have he (x : X) :
      (if f x = false then P.mass x else 0) + (if f x = true then P.mass x else 0) =
        P.mass x := by cases f x <;> simp
  have h := congrArg (fun g : X → ℝ ↦ ∑ x, g x) (funext he)
  simp only [sum_add_distrib, P.total] at h
  unfold bin
  linarith

theorem fair_taper_bit (m : ℕ) (b : Bool) :
    bin (fairLaw (L m)).mass taperBit b = bitMass (1 / 8) b := by
  have ht : bin (fairLaw (L m)).mass taperBit true = 1 / 8 := by
    simpa only [fairLaw, taperBit, bin, decide_eq_true_eq] using fair_asymmetric m
  cases b
  · rw [binary_complement, ht]
    rfl
  · exact ht

theorem source_taper_bit (m : ℕ) {a : ℝ} (ha : 0 < a) (hu : a < 1 / 2) (b : Bool) :
    bin (sourceLaw m ha hu).mass taperBit b = bitMass ((1 + a) / 8) b := by
  have ht : bin (sourceLaw m ha hu).mass taperBit true = (1 + a) / 8 := by
    simpa only [sourceLaw, taperBit, bin, decide_eq_true_eq] using perturbed_asymmetric m a
  cases b
  · rw [binary_complement, ht]
    rfl
  · exact ht

open scoped Classical in
noncomputable def rectangularLaw {m : ℕ} (P : Law (Word (L m))) (w o : ℕ)
    (ho : o + w ≤ L m) :
    @Law {v : Fin w → ℝ // v ∈ univ.image (fun x ↦ powerVector (slice x o w ho))}
      (Finset.Subtype.fintype _) :=
  observed P (fun x ↦ powerVector (slice x o w ho))

theorem rectangular_law_eq {m : ℕ} {a : ℝ} (ha : 0 < a) (hu : a < 1 / 2)
    (w o : ℕ) (ho : o + w ≤ L m) :
    rectangularLaw (sourceLaw m ha hu) w o ho = rectangularLaw (fairLaw (L m)) w o ho := by
  classical
  apply observed_eq
  intro v
  have h0 := all_width_power_laws ho 0 v
  have ha0 : (fun _ : Word (L m) ↦ uniform (L m)) = weight m 0 := by
    funext x
    simp [weight]
  change bin (weight m a) (fun x ↦ powerVector (slice x o w ho)) v =
    bin (fun _ ↦ uniform (L m)) (fun x ↦ powerVector (slice x o w ho)) v
  rw [all_width_power_laws ho a, ha0, h0]

theorem rectangular_experiment_eq {m : ℕ} {a : ℝ} (ha : 0 < a) (hu : a < 1 / 2)
    (r : ℕ) (w o : Fin r → ℕ) (ho : ∀ i, o i + w i ≤ L m) :
    independent (fun i ↦ rectangularLaw (sourceLaw m ha hu) (w i) (o i) (ho i)) =
      independent (fun i ↦ rectangularLaw (fairLaw (L m)) (w i) (o i) (ho i)) := by
  have h : (fun i ↦ rectangularLaw (sourceLaw m ha hu) (w i) (o i) (ho i)) =
      (fun i ↦ rectangularLaw (fairLaw (L m)) (w i) (o i) (ho i)) :=
    funext (fun i ↦ rectangular_law_eq ha hu (w i) (o i) (ho i))
  rw [h]

theorem iid_observation_bins (n : ℕ) (P : Law X) {Z : Fin n → Type*}
    [∀ i, DecidableEq (Z i)] (f : ∀ i, X → Z i) (z : ∀ i, Z i) :
    bin (iid n P).mass (fun x i ↦ f i (x i)) z = ∏ i, bin P.mass (f i) (z i) := by
  classical
  unfold bin iid
  rw [Fintype.prod_sum]
  apply sum_congr rfl
  intro x _
  simp only [prod_ite_zero, mem_univ, forall_const, funext_iff]

theorem actual_taper_records (m n : ℕ) {a : ℝ} (ha : 0 < a) (hu : a < 1 / 2)
    (x : Word n) :
    (∏ i, bin (sourceLaw m ha hu).mass taperBit (x i)) =
      (iid n (bernoulli ((1 + a) / 8) (taper_parameters ha hu).2)).mass x ∧
    (∏ i, bin (fairLaw (L m)).mass taperBit (x i)) =
      (iid n (bernoulli (1 / 8) (taper_parameters ha hu).1)).mass x := by
  constructor <;> simp only [source_taper_bit, fair_taper_bit, iid, bernoulli]

theorem raw_taper_records (m n : ℕ) {a : ℝ} (ha : 0 < a) (hu : a < 1 / 2)
    (x : Word n) :
    bin (iid n (sourceLaw m ha hu)).mass (fun records i ↦ taperBit (records i)) x =
      (iid n (bernoulli ((1 + a) / 8) (taper_parameters ha hu).2)).mass x ∧
    bin (iid n (fairLaw (L m))).mass (fun records i ↦ taperBit (records i)) x =
      (iid n (bernoulli (1 / 8) (taper_parameters ha hu).1)).mass x := by
  rw [iid_observation_bins, iid_observation_bins]
  exact actual_taper_records m n ha hu x

theorem raw_rectangular_records {m : ℕ} {a : ℝ} (ha : 0 < a) (hu : a < 1 / 2)
    (r : ℕ) (w o : Fin r → ℕ) (ho : ∀ i, o i + w i ≤ L m)
    (v : ∀ i, Fin (w i) → ℝ) :
    bin (iid r (sourceLaw m ha hu)).mass
        (fun records i ↦ powerVector (slice (records i) (o i) (w i) (ho i))) v =
      bin (iid r (fairLaw (L m))).mass
        (fun records i ↦ powerVector (slice (records i) (o i) (w i) (ho i))) v := by
  classical
  rw [iid_observation_bins r _ (fun i x ↦ powerVector (slice x (o i) (w i) (ho i))),
    iid_observation_bins r _ (fun i x ↦ powerVector (slice x (o i) (w i) (ho i)))]
  apply prod_congr rfl
  intro i _
  have ha0 := all_width_power_laws (ho i) 0 (v i)
  have ha1 := all_width_power_laws (ho i) a (v i)
  have hz : weight m 0 = (fun _ ↦ uniform (L m)) := by
    funext x
    simp [weight]
  change bin (weight m a) (fun x ↦ powerVector (slice x (o i) (w i) (ho i))) (v i) =
    bin (fun _ ↦ uniform (L m)) (fun x ↦ powerVector (slice x (o i) (w i) (ho i))) (v i)
  rw [← hz]
  exact ha1.trans ha0.symm

theorem taper_sites_and_span (m : ℕ) :
    (selected m 0).val = 0 ∧ (selected m 1).val = 1 ∧
      (selected m 2).val = 2 ∧ (selected m 3).val = D m ∧ L m = D m + 1 := by
  exact ⟨rfl, rfl, rfl, rfl, rfl⟩

theorem every_depth {d : ℕ} (hd : 4 ≤ d) : ∃ m, d = D m := by
  exact ⟨d - 4, by dsimp [D]; omega⟩

#print axioms Law.ext
#print axioms observed_eq
#print axioms binary_complement
#print axioms fair_taper_bit
#print axioms source_taper_bit
#print axioms rectangular_law_eq
#print axioms rectangular_experiment_eq
#print axioms iid_observation_bins
#print axioms actual_taper_records
#print axioms raw_taper_records
#print axioms raw_rectangular_records
#print axioms taper_sites_and_span
#print axioms every_depth

end TimeFrequency.ContextTree.DepthEvidence
