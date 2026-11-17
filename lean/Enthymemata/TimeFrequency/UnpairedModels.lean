import Enthymemata.TimeFrequency.PairedViewLaws
import Mathlib.Data.Matrix.Mul

/-! # Identification from separate labeled views on a finite dictionary -/

namespace TimeFrequency.UnpairedModels

open Set Finset Matrix

variable {C J : Type*} {V : J → Type*} (M : (j : J) → C → V j)

def hull (S : Set C) : Set C := {c | ∀ j, M j c ∈ M j '' S}

def Private (S : Set C) (x : C) : Prop := ∃ j, ∀ y ∈ S, M j y = M j x → y = x

theorem subset_hull (S : Set C) : S ⊆ hull M S := fun x hx _ ↦ ⟨x,hx,rfl⟩

theorem hull_images (S : Set C) (j : J) : M j '' hull M S = M j '' S := by
  apply Subset.antisymm _ (image_mono (subset_hull M S))
  rintro _ ⟨x,hx,rfl⟩
  exact hx j

theorem compatible_subset_hull {S A : Set C} (h : ∀ j, M j '' A = M j '' S) : A ⊆ hull M S := by
  intro x hx j
  rw [← h]
  exact ⟨x,hx,rfl⟩

theorem insert_images {S : Set C} {x : C} (hx : x ∈ hull M S) (j : J) :
    M j '' (insert x S) = M j '' S := by
  apply Subset.antisymm _ (image_mono (Set.subset_insert x S))
  rintro _ ⟨y,hy,rfl⟩
  rcases Set.mem_insert_iff.mp hy with rfl | hy
  · exact hx j
  · exact ⟨y,hy,rfl⟩

theorem erase_images {S : Set C} {x : C} (h : ¬ Private M S x) (j : J) :
    M j '' (S \ {x}) = M j '' S := by
  classical
  apply Subset.antisymm (image_mono Set.sdiff_subset) _
  rintro _ ⟨y,hy,rfl⟩
  by_cases he : y = x
  · subst y
    have hn : ¬ ∀ y ∈ S, M j y = M j x → y = x := fun hp ↦ h ⟨j,hp⟩
    push Not at hn
    obtain ⟨y,hy,he,hd⟩ := hn
    exact ⟨y,⟨hy,hd⟩,he⟩
  · exact ⟨y,⟨hy,he⟩,rfl⟩

theorem support_unique_iff (S : Set C) :
    (∀ A : Set C, (∀ j, M j '' A = M j '' S) → A = S) ↔
      hull M S = S ∧ ∀ x ∈ S, Private M S x := by
  classical
  constructor
  · intro h
    refine ⟨h _ (hull_images M S), ?_⟩
    intro x hx
    by_contra hn
    have he := h _ (erase_images M hn)
    have : x ∈ S \ {x} := he.symm ▸ hx
    exact this.2 rfl
  · rintro ⟨hh,hp⟩ A ha
    have hAS : A ⊆ S := hh ▸ compatible_subset_hull M ha
    apply Subset.antisymm hAS
    intro x hx
    obtain ⟨j,hj⟩ := hp x hx
    have hv : M j x ∈ M j '' A := (ha j).symm ▸ ⟨x,hx,rfl⟩
    obtain ⟨y,hy,he⟩ := hv
    exact hj y (hAS hy) he ▸ hy

section Finite

variable [Fintype C] [∀ j, DecidableEq (V j)]

/-- Each block uses exactly the attained values of its operator. -/
def alphabet (j : J) := Finset.univ.image (M j)

abbrev Bin (j : J) := {v // v ∈ alphabet M j}

def incidence : Matrix ((j : J) × Bin M j) C ℝ :=
  fun row c ↦ if M row.1 c = row.2.val then 1 else 0

theorem incidence_marginal (p : C → ℝ) (j : J) (b : Bin M j) :
    (incidence M).mulVec p ⟨j,b⟩ = ∑ c, if M j c = b.val then p c else 0 := by
  simp only [Matrix.mulVec, dotProduct, incidence, ite_mul, one_mul, zero_mul]

theorem block_mass (p : C → ℝ) (j : J) :
    ∑ b : Bin M j, (incidence M).mulVec p ⟨j,b⟩ = ∑ c, p c := by
  simp_rw [incidence_marginal]
  rw [sum_comm]
  apply sum_congr rfl
  intro c _
  let b : Bin M j := ⟨M j c, mem_image.mpr ⟨c,mem_univ _,rfl⟩⟩
  rw [sum_eq_single b]
  · simp [b]
  · intro v _ hv
    exact if_neg (fun he ↦ hv (Subtype.ext he.symm))
  · simp

def Probability (p : C → ℝ) : Prop := (∀ c, 0 ≤ p c) ∧ ∑ c, p c = 1

/-- A zero-sum signed direction gives two probability laws by its Jordan decomposition. -/
theorem kernel_witness {R : Type*} (D : Matrix R C ℝ) {z : C → ℝ} (hz : z ≠ 0)
    (hm : ∑ c, z c = 0) (hk : D.mulVec z = 0) :
    ∃ p q, Probability p ∧ Probability q ∧ p ≠ q ∧ D.mulVec p = D.mulVec q := by
  classical
  let u := fun c ↦ max (z c) 0
  let v := fun c ↦ max (-z c) 0
  have hu (c) : 0 ≤ u c := le_max_right _ _
  have hv (c) : 0 ≤ v c := le_max_right _ _
  have he : z = u-v := by
    funext c
    change z c = max (z c) 0 - max (-z c) 0
    rcases le_total 0 (z c) with h | h <;> simp [h]
  have hsum : (∑ c, u c) = ∑ c, v c := by
    rw [he] at hm
    simpa only [Pi.sub_apply, sum_sub_distrib, sub_eq_zero] using hm
  have hn : 0 < ∑ c, u c := by
    by_contra! hl
    have h0 := le_antisymm hl (sum_nonneg (fun c _ ↦ hu c))
    have uz := (sum_eq_zero_iff_of_nonneg (fun c _ ↦ hu c)).mp h0
    have vz := (sum_eq_zero_iff_of_nonneg (fun c _ ↦ hv c)).mp (hsum ▸ h0)
    apply hz
    rw [he]
    funext c
    simp [uz c (mem_univ _), vz c (mem_univ _)]
  let s := (∑ c, u c)⁻¹
  have hs : 0 < s := inv_pos.mpr hn
  have hd : D.mulVec u = D.mulVec v := by
    rw [he, Matrix.mulVec_sub, sub_eq_zero] at hk
    exact hk
  refine ⟨s • u,s • v,⟨fun c ↦ mul_nonneg hs.le (hu c), ?_⟩,
    ⟨fun c ↦ mul_nonneg hs.le (hv c), ?_⟩, ?_, ?_⟩
  · simp only [Pi.smul_apply, smul_eq_mul, ← mul_sum, s, inv_mul_cancel₀ hn.ne']
  · simp only [Pi.smul_apply, smul_eq_mul, ← mul_sum, ← hsum, s, inv_mul_cancel₀ hn.ne']
  · intro huv
    apply hz
    rw [he,sub_eq_zero]
    funext c
    exact mul_left_cancel₀ hs.ne' (congrFun huv c)
  · rw [Matrix.mulVec_smul, Matrix.mulVec_smul, hd]

theorem simplex_injective_iff {R : Type*} (D : Matrix R C ℝ) :
    (∀ p q, Probability p → Probability q → D.mulVec p = D.mulVec q → p = q) ↔
    (∀ z, D.mulVec z = 0 → ∑ c, z c = 0 → z = 0) := by
  constructor
  · intro hi z hk hm
    by_contra hn
    obtain ⟨p,q,hp,hq,hd,he⟩ := kernel_witness D hn hm hk
    exact hd (hi p q hp hq he)
  · intro hk p q hp hq he
    apply sub_eq_zero.mp
    apply hk
    · rw [Matrix.mulVec_sub,he,sub_self]
    · simp [sum_sub_distrib,hp.2,hq.2]

theorem complete_block_kernel [Nonempty J] {z : C → ℝ} (hk : (incidence M).mulVec z = 0) :
    ∑ c, z c = 0 := by
  obtain ⟨j⟩ := ‹Nonempty J›
  rw [← block_mass M z j, hk]
  simp

end Finite

section FourierDictionary

open PairedViews

/-- The indices refer, in order, to the four signed representatives of slice 5. -/
def wave (i : Fin 4) : Wave :=
  if i.val = 0 then (1,2) else if i.val = 1 then (2,-1) else if i.val = 2 then (1,-2) else (2,1)

def viewCode (j : Bool) (i : Fin 4) : Bool :=
  decide (i.val = 1 ∨ i.val = (if j then 2 else 3))

def decode (j b : Bool) : Wave :=
  if j then (if b then (1,3) else (3,1)) else (if b then (2,1) else (1,2))

theorem decode_injective (j : Bool) : Function.Injective (decode j) := by
  intro b c h
  cases j <;> cases b <;> cases c <;> norm_num [decode] at *

theorem four_physical_views (j : Bool) (i : Fin 4) :
    (if j then frameView else singletonView) (wave i) = decode j (viewCode j i) := by
  fin_cases i <;> cases j <;> norm_num [wave,viewCode,decode,singletonView,frameView]

def fourDirection (i : Fin 4) : ℝ := if i.val < 2 then 1 else -1

theorem four_kernel :
    fourDirection ≠ 0 ∧ ∑ i, fourDirection i = 0 ∧
      (incidence viewCode).mulVec fourDirection = 0 := by
  refine ⟨fun hz ↦ ?_, ?_, ?_⟩
  · have := congrFun hz 0
    norm_num [fourDirection] at this
  · norm_num [Fin.sum_univ_four, fourDirection]
  · funext row
    obtain ⟨j,⟨b,hb⟩⟩ := row
    rw [incidence_marginal]
    cases j <;> cases b <;> norm_num [Fin.sum_univ_four,viewCode,fourDirection]

/-- The two occupied-bin lists are exactly the same for these disjoint supports. -/
theorem four_supports (j : Bool) :
    ({0,1} : Finset (Fin 4)).image (viewCode j) = ({2,3} : Finset (Fin 4)).image (viewCode j) ∧
    ({0,1} : Finset (Fin 4)).image (viewCode j) = Finset.univ.image (viewCode j) := by
  cases j <;> decide

def threeCode (j : Bool) (i : Fin 3) := viewCode j i.castSucc

theorem three_collisions (j : Bool) : ∃ x y, x ≠ y ∧ threeCode j x = threeCode j y := by
  cases j <;> decide

private def attained {A B : Type*} [Fintype A] [DecidableEq B] (f : Bool → A → B)
    (j : Bool) (x : A) : Bin f j := ⟨f j x, mem_image.mpr ⟨x,mem_univ _,rfl⟩⟩

/-- The P and B rows give the first and second mass. The A row gives the first plus the third. -/
theorem three_injective : Function.Injective (incidence threeCode).mulVec := by
  intro p q he
  have h0 := congrFun he ⟨true, attained threeCode true 0⟩
  have h1 := congrFun he ⟨false, attained threeCode false 1⟩
  have h02 := congrFun he ⟨false, attained threeCode false 0⟩
  simp only [incidence_marginal] at h0 h1 h02
  norm_num [attained, Fin.sum_univ_three, threeCode, viewCode] at h0 h1 h02
  funext i
  fin_cases i <;> first | exact h0 | exact h1 | (change p 2 = q 2; linarith)

/-- Complete enumeration of the occupied-bin experiment on the three-element dictionary. -/
theorem three_supports :
    (∀ A : Finset (Fin 3),
      (∀ j, A.image (threeCode j) = ({2} : Finset (Fin 3)).image (threeCode j)) → A = {2}) ∧
    ({0,1} : Finset (Fin 3)) ≠ Finset.univ ∧
    (∀ j, ({0,1} : Finset (Fin 3)).image (threeCode j) = Finset.univ.image (threeCode j)) := by
  decide

end FourierDictionary

#print axioms decode_injective
#print axioms four_physical_views
#print axioms four_kernel
#print axioms four_supports
#print axioms three_collisions
#print axioms three_injective
#print axioms three_supports
#print axioms subset_hull
#print axioms hull_images
#print axioms compatible_subset_hull
#print axioms insert_images
#print axioms erase_images
#print axioms support_unique_iff
#print axioms incidence_marginal
#print axioms block_mass
#print axioms kernel_witness
#print axioms simplex_injective_iff
#print axioms complete_block_kernel

end TimeFrequency.UnpairedModels
