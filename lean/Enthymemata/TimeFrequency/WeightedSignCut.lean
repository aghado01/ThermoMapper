import Enthymemata.TimeFrequency.TwoSampleRecovery
import Mathlib.Analysis.InnerProductSpace.PiL2

/-! # Sharp conditioning of edge products at fixed real amplitudes -/

namespace TimeFrequency.SignCut

open Finset SimpleGraph

variable {V : Type*} [Fintype V] [LinearOrder V]

noncomputable def l2 {J : Type*} [Fintype J] (z : J → ℝ) : ℝ :=
  ‖(WithLp.toLp 2 z : EuclideanSpace ℝ J)‖

theorem l2_sq {J : Type*} [Fintype J] (z : J → ℝ) : l2 z ^ 2 = ∑ j, z j ^ 2 :=
  EuclideanSpace.real_norm_sq_eq _

def flip (x : V → ℝ) (S : Finset V) (u : V) : ℝ := if u ∈ S then -x u else x u

noncomputable def distance (x y : V → ℝ) := min (l2 (x-y)) (l2 (x+y))

/-- Ordering only selects one copy of each unordered edge. Other coordinates are zero. -/
noncomputable def products (G : SimpleGraph V) (x : V → ℝ) (e : V × V) :=
  open Classical in if e.1 < e.2 ∧ G.Adj e.1 e.2 then x e.1*x e.2 else 0

noncomputable def mass (x : V → ℝ) (S : Finset V) := ∑ v ∈ S, x v ^ 2

noncomputable def light (x : V → ℝ) (S : Finset V) := min (mass x S) (mass x Sᶜ)

noncomputable def boundary (G : SimpleGraph V) (x : V → ℝ) (S : Finset V) :=
  ∑ e : V × V, if (e.1 ∈ S ↔ e.2 ∈ S) then 0 else products G x e ^ 2

omit [Fintype V] in
theorem flip_abs (x : V → ℝ) (S : Finset V) (u : V) : |flip x S u| = |x u| := by
  by_cases h : u ∈ S <;> simp [flip,h]

omit [Fintype V] in
theorem exists_flip [Finite V] {x y : V → ℝ} (h : ∀ u, |y u| = |x u|) : ∃ S, y = flip x S := by
  letI := Fintype.ofFinite V
  refine ⟨univ.filter fun u ↦ y u ≠ x u, funext fun u ↦ ?_⟩
  rcases abs_eq_abs.mp (h u) with he | he
  · simp [flip, he]
  · by_cases hsame : y u = x u <;> simp [flip, hsame, he]

theorem distance_flip_sq (x : V → ℝ) (S : Finset V) :
    distance x (flip x S) ^ 2 = 4 * light x S := by
  have he : flip x Sᶜ = -flip x S := by
    funext v
    by_cases h : v ∈ S <;> simp [flip,h]
  have hboth (T : Finset V) : l2 (x - flip x T)^2 = 4 * mass x T := by
    rw [l2_sq, mass, Finset.mul_sum]
    calc ∑ v, ((x - flip x T) v)^2 = ∑ v, if v ∈ T then 4*x v^2 else 0 := by
          apply Finset.sum_congr rfl
          intro v _
          by_cases h : v ∈ T <;> simp [flip,h]
          ring
      _ = _ := by simp
  have ha := hboth Sᶜ
  have hd := hboth S
  rw [he, sub_neg_eq_add] at ha
  unfold distance light
  rcases le_total (l2 (x-flip x S)) (l2 (x+flip x S)) with h | h
  · rw [min_eq_left h, hd, min_eq_left]
    have hs : l2 (x-flip x S)^2 ≤ l2 (x+flip x S)^2 := pow_le_pow_left₀ (norm_nonneg _) h 2
    rw [hd, ha] at hs
    linarith
  · rw [min_eq_right h, ha, min_eq_right]
    have hs : l2 (x+flip x S)^2 ≤ l2 (x-flip x S)^2 := pow_le_pow_left₀ (norm_nonneg _) h 2
    rw [ha, hd] at hs
    linarith

theorem products_flip_sq (G : SimpleGraph V) (x : V → ℝ) (S : Finset V) :
    l2 (products G x - products G (flip x S))^2 = 4 * boundary G x S := by
  rw [l2_sq, boundary, Finset.mul_sum]
  apply Finset.sum_congr rfl
  intro e _
  by_cases he : e.1 < e.2 ∧ G.Adj e.1 e.2 <;>
    by_cases hu : e.1 ∈ S <;> by_cases hv : e.2 ∈ S <;>
      simp [products, flip, he, hu, hv] <;> ring

theorem complementary_cut (G : SimpleGraph V) (x : V → ℝ) (S : Finset V) :
    boundary G x Sᶜ / light x Sᶜ = boundary G x S / light x S := by
  simp [boundary, light, min_comm, not_iff_not]

noncomputable def properCuts (x : V → ℝ) := univ.filter fun S : Finset V ↦ 0 < light x S

noncomputable def cutConstant (G : SimpleGraph V) (x : V → ℝ) (h : (properCuts x).Nonempty) :=
  ((properCuts x).image fun S ↦ boundary G x S / light x S).min' (h.image _)

theorem constant_attained (G : SimpleGraph V) (x : V → ℝ) (h : (properCuts x).Nonempty) :
    ∃ S, 0 < light x S ∧ boundary G x S / light x S = cutConstant G x h := by
  have hm := Finset.min'_mem _ (h.image (fun S ↦ boundary G x S / light x S))
  obtain ⟨S,hS,he⟩ := mem_image.mp hm
  exact ⟨S, (mem_filter.mp hS).2, he⟩

theorem sharp_constant_iff (G : SimpleGraph V) (x : V → ℝ) (h : (properCuts x).Nonempty) (k : ℝ) :
    (∀ y : V → ℝ, (∀ v, |y v| = |x v|) →
      k * distance x y ^ 2 ≤ l2 (products G x - products G y) ^ 2) ↔ k ≤ cutConstant G x h := by
  constructor
  · intro hk
    obtain ⟨S,hS,he⟩ := constant_attained G x h
    have hz := hk _ (flip_abs x S)
    rw [distance_flip_sq, products_flip_sq] at hz
    rw [← he, le_div_iff₀ hS]
    linarith
  · intro hk y hy
    obtain ⟨S,rfl⟩ := exists_flip hy
    rw [distance_flip_sq, products_flip_sq]
    have hm : 0 ≤ light x S :=
      le_min (sum_nonneg fun _ _ ↦ sq_nonneg _) (sum_nonneg fun _ _ ↦ sq_nonneg _)
    by_cases hs : 0 < light x S
    · have hl : cutConstant G x h ≤ boundary G x S / light x S :=
        Finset.min'_le _ _ (mem_image.mpr ⟨S, by simpa [properCuts] using hs, rfl⟩)
      linarith [(le_div_iff₀ hs).mp (hk.trans hl)]
    · have hz : light x S = 0 := le_antisymm (le_of_not_gt hs) hm
      have hb : 0 ≤ boundary G x S := sum_nonneg fun e _ ↦ by split <;> positivity
      rw [hz, mul_zero, mul_zero]
      positivity

omit [Fintype V] [LinearOrder V] in
theorem mass_pos (x : V → ℝ) (S : Finset V) (v : V) (hv : v ∈ S) (hx : x v ≠ 0) :
    0 < mass x S := sum_pos' (fun _ _ ↦ sq_nonneg _) ⟨v,hv,sq_pos_of_ne_zero hx⟩

theorem extremizing_flip (G : SimpleGraph V) (x : V → ℝ) (h : (properCuts x).Nonempty) :
    ∃ y : V → ℝ, (∀ v, |y v| = |x v|) ∧ 0 < distance x y ^ 2 ∧
      l2 (products G x - products G y)^2 = cutConstant G x h * distance x y ^ 2 := by
  obtain ⟨S,hs,he⟩ := constant_attained G x h
  refine ⟨flip x S, flip_abs x S, ?_, ?_⟩
  · rw [distance_flip_sq]
    positivity
  · rw [distance_flip_sq, products_flip_sq, (div_eq_iff hs.ne').mp he]
    ring

theorem properCuts_nonempty {x : V → ℝ} {u v : V} (huv : u ≠ v) (hu : x u ≠ 0) (hv : x v ≠ 0) :
    (properCuts x).Nonempty := by
  refine ⟨{u}, mem_filter.mpr ⟨mem_univ _, lt_min (mass_pos x _ u (by simp) hu) ?_⟩⟩
  exact mass_pos x _ v (by simp [huv.symm]) hv

omit [Fintype V] [LinearOrder V] in
theorem subsingleton_support_sign {x y : V → ℝ}
    (hs : ∀ u v, x u ≠ 0 → x v ≠ 0 → u = v) (hy : ∀ u, |y u| = |x u|) :
    y = x ∨ y = -x := by
  by_cases h : ∃ u, x u ≠ 0
  · obtain ⟨u,hu⟩ := h
    have hz (v) (hv : v ≠ u) : y v = 0 ∧ x v = 0 := by
      have hx : x v = 0 := by by_contra hn; exact hv (hs v u hn hu)
      exact ⟨abs_eq_zero.mp (by simpa [hx] using hy v), hx⟩
    rcases abs_eq_abs.mp (hy u) with h | h
    · left
      funext v
      by_cases hv : v = u
      · simpa [hv] using h
      · rw [(hz _ hv).1, (hz _ hv).2]
    · right
      funext v
      by_cases hv : v = u
      · simpa [hv] using h
      · simp [(hz _ hv).1, (hz _ hv).2]
  · push Not at h
    left
    funext v
    simp [h, abs_eq_zero.mp (by simpa [h] using hy v : |y v| = 0)]

theorem boundary_zero_iff (G : SimpleGraph V) {x : V → ℝ} (hx : ∀ v, x v ≠ 0) (S : Finset V) :
    boundary G x S = 0 ↔ ∀ u v, G.Adj u v → (u ∈ S ↔ v ∈ S) := by
  unfold boundary
  rw [sum_eq_zero_iff_of_nonneg (fun e _ ↦ by split <;> positivity)]
  constructor
  · intro h
    have hor (u v : V) (huv : u < v) (ha : G.Adj u v) : u ∈ S ↔ v ∈ S := by
      by_contra hn
      have he := h (u,v) (mem_univ _)
      simp [hn, products, huv, ha, hx u, hx v] at he
    intro u v ha
    rcases lt_trichotomy u v with huv | rfl | hvu
    · exact hor _ _ huv ha
    · rfl
    · exact (hor _ _ hvu ha.symm).symm
  · intro h e _
    by_cases he : e.1 < e.2 ∧ G.Adj e.1 e.2
    · simp [h _ _ he.2]
    · simp [products, he]

/-- Apply to the graph induced on the nonzero support. The positive-cut hypothesis supplies
both sides, so the minimum never needs a convention for an empty family. -/
theorem constant_zero_iff_disconnected (G : SimpleGraph V) {x : V → ℝ}
    (hx : ∀ v, x v ≠ 0) (h : (properCuts x).Nonempty) :
    cutConstant G x h = 0 ↔ ¬G.Preconnected := by
  classical
  have hnonneg : 0 ≤ cutConstant G x h := by
    obtain ⟨S,hs,he⟩ := constant_attained G x h
    rw [← he]
    exact div_nonneg (sum_nonneg fun e _ ↦ by split <;> positivity) hs.le
  constructor
  · intro hz hc
    obtain ⟨S,hs,he⟩ := constant_attained G x h
    rw [hz, div_eq_zero_iff] at he
    have hstep := (boundary_zero_iff G hx S).mp (he.resolve_right hs.ne')
    have hn (T : Finset V) (hT : 0 < mass x T) : T.Nonempty := by
      by_contra! he
      simp [he, mass] at hT
    obtain ⟨u,hu⟩ := hn S ((lt_min_iff.mp hs).1)
    obtain ⟨v,hv⟩ := hn Sᶜ ((lt_min_iff.mp hs).2)
    have hw {i j} (p : G.Walk i j) : i ∈ S ↔ j ∈ S := by
      induction p with
      | nil => rfl
      | cons ha _ ih => exact (hstep _ _ ha).trans ih
    exact (mem_compl.mp hv) ((hw (hc u v).some).mp hu)
  · intro hc
    change ¬∀ u v, G.Reachable u v at hc
    push Not at hc
    obtain ⟨u,v,huv⟩ := hc
    let S := univ.filter fun i ↦ G.Reachable u i
    have hS : 0 < light x S := lt_min
      (mass_pos x _ u (by simp [S]) (hx u))
      (mass_pos x _ v (by simp [S, huv]) (hx v))
    have hz : boundary G x S = 0 := (boundary_zero_iff G hx S).mpr (by
      intro i j ha
      simp only [S, mem_filter, mem_univ, true_and]
      exact ⟨fun hi ↦ hi.trans ha.reachable, fun hj ↦ hj.trans ha.symm.reachable⟩)
    have hle : cutConstant G x h ≤ boundary G x S / light x S :=
      Finset.min'_le _ _ (mem_image.mpr ⟨S, by simpa [properCuts] using hS, rfl⟩)
    exact le_antisymm (by simpa [hz] using hle) hnonneg

/-- The two unnormalized length-two powers change by opposite signed amounts.
At fixed amplitudes their squared Euclidean discrepancy is eight times the product discrepancy. -/
theorem raw_pair_factor {a b c d : ℝ} (ha : |a| = |c|) (hb : |b| = |d|) :
    (framePower ![a,b] 0 - framePower ![c,d] 0)^2 +
      (framePower ![a,b] 1 - framePower ![c,d] 1)^2 = 8*(a*b-c*d)^2 := by
  rw [(twoSample_powers a b).1, (twoSample_powers a b).2,
    (twoSample_powers c d).1, (twoSample_powers c d).2]
  rcases abs_eq_abs.mp ha with h | h <;> rcases abs_eq_abs.mp hb with h' | h' <;> rw [h,h'] <;> ring

/-- A single quiet sample on the path separates two unit-amplitude endpoints. -/
theorem quiet_bridge (δ : ℝ) :
    let x : Fin 3 → ℝ := ![1,δ,1]
    let y := flip x {2}
    distance x y ^ 2 = 4 ∧ l2 (products (pathGraph 3) x - products (pathGraph 3) y)^2 = 4*δ^2 := by
  dsimp only
  rw [distance_flip_sq, products_flip_sq]
  have hm : mass (![1,δ,1] : Fin 3 → ℝ) ({2}ᶜ : Finset (Fin 3)) = 1 + δ^2 := by
    have he : ({2}ᶜ : Finset (Fin 3)) = {0,1} := by decide
    rw [he]
    simp [mass]
  constructor
  · rw [light, hm]
    simp [mass, min_eq_left (by linarith [sq_nonneg δ] : (1 : ℝ) ≤ 1 + δ^2)]
  · simp [boundary, products, Fintype.sum_prod_type, Fin.sum_univ_three, pathGraph_adj]

theorem quiet_bridge_zero :
    SameTwoSampleViews (![1,0,1] : Fin 3 → ℝ) ![1,0,-1] ∧
    (![1,0,-1] : Fin 3 → ℝ) ≠ ![1,0,1] ∧ (![1,0,-1] : Fin 3 → ℝ) ≠ -![1,0,1] := by
  have h := gap_witness (1 : ℝ) 1 (by norm_num)
  exact ⟨h.1,fun he ↦ h.2.1 (.inl he), fun he ↦ h.2.1 (.inr he)⟩

/-- The first frame is the same in both waveforms. The second frame accounts for all
of the squared discrepancy in the four retained Fourier powers. -/
theorem quiet_bridge_raw (δ : ℝ) :
    (framePower ![δ,1] 0 - framePower ![δ,-1] 0)^2 +
      (framePower ![δ,1] 1 - framePower ![δ,-1] 1)^2 = 32*δ^2 := by
  rw [raw_pair_factor (a := δ) (c := δ) (b := 1) (d := -1) rfl (by simp)]
  ring

#print axioms l2_sq
#print axioms flip_abs
#print axioms exists_flip
#print axioms distance_flip_sq
#print axioms products_flip_sq
#print axioms complementary_cut
#print axioms constant_attained
#print axioms sharp_constant_iff
#print axioms mass_pos
#print axioms extremizing_flip
#print axioms properCuts_nonempty
#print axioms subsingleton_support_sign
#print axioms boundary_zero_iff
#print axioms constant_zero_iff_disconnected
#print axioms raw_pair_factor
#print axioms quiet_bridge
#print axioms quiet_bridge_zero
#print axioms quiet_bridge_raw

end TimeFrequency.SignCut
