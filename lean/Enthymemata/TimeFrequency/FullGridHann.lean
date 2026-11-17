import Mathlib.Analysis.Fourier.ZMod
import Mathlib.Combinatorics.SimpleGraph.Connectivity.Connected
import Mathlib.Analysis.SpecialFunctions.Trigonometric.Basic

/-! # Full-grid periodic Hann frames on a zero-extended real waveform -/

namespace TimeFrequency.FullGrid

open Finset Complex
open scoped ComplexConjugate

variable {F : ℕ} [NeZero F]

/-- Circular autocorrelation on the frequency grid's time domain. -/
noncomputable def correlation (f : ZMod F → ℂ) (d : ZMod F) :=
  ∑ j, f (j+d)*conj (f j)

/-- The finite Fourier transform of the circular autocorrelation is the full power spectrum. -/
theorem dft_correlation (f : ZMod F → ℂ) (k : ZMod F) :
    ZMod.dft (correlation f) k = ↑(normSq (ZMod.dft f k)) := by
  rw [← Complex.mul_conj]
  simp only [ZMod.dft_apply, smul_eq_mul, correlation, mul_sum, map_sum, map_mul, sum_mul]
  rw [Finset.sum_comm]
  apply sum_congr rfl
  intro j _
  refine Fintype.sum_equiv (Equiv.addLeft j) _ _ fun d ↦ ?_
  simp only [show (Equiv.addLeft j) d = j+d from rfl, ← AddChar.map_neg_eq_conj, neg_neg]
  have hc : ZMod.stdAddChar (-(d*k)) =
      ZMod.stdAddChar (-((j+d)*k)) * ZMod.stdAddChar (j*k) := by
    rw [← AddChar.map_add_eq_mul]
    congr 1
    ring
  rw [hc]
  ring

theorem powers_iff_correlation (f g : ZMod F → ℂ) :
    (∀ k, normSq (ZMod.dft f k) = normSq (ZMod.dft g k)) ↔ correlation f = correlation g := by
  constructor
  · intro h
    apply ZMod.dft.injective
    funext k
    rw [dft_correlation, dft_correlation, h]
  · intro h k
    apply Complex.ofReal_injective
    rw [← dft_correlation, ← dft_correlation, h]

/-- A finite waveform is zero off its support on ℤ. All integer starts are permitted.
Grid positions beyond the specified window are padded with zeros. -/
noncomputable def frame (W : ℕ) (g : ℕ → ℝ) (x : ℤ →₀ ℝ) (r : ℤ) (j : ZMod F) : ℂ :=
  if j.val < W then ↑(g j.val * x (r+j.val)) else 0

def SamePowers (F : ℕ) [NeZero F] (W : ℕ) (g : ℕ → ℝ) (x y : ℤ →₀ ℝ) : Prop :=
  ∀ (r : ℤ) (k : ZMod F), normSq (ZMod.dft (frame W g x r) k) =
    normSq (ZMod.dft (frame W g y r) k)

/-- For an N-sample recording only the starts 1-W through N-1 can contribute. -/
theorem finite_starts {W N : ℕ} (g : ℕ → ℝ) (x y : ℤ →₀ ℝ)
    (hx : ∀ t : ℤ, t < 0 ∨ N ≤ t → x t = 0) (hy : ∀ t : ℤ, t < 0 ∨ N ≤ t → y t = 0) :
    SamePowers F W g x y ↔ ∀ r : ℤ, r ∈ Set.Ico (1 - (W : ℤ)) (N : ℤ) →
      ∀ k : ZMod F, normSq (ZMod.dft (frame W g x r) k) = normSq (ZMod.dft (frame W g y r) k) := by
  refine ⟨fun h r _ ↦ h r, ?_⟩
  intro h r k
  by_cases hr : r ∈ Set.Ico (1 - (W : ℤ)) (N : ℤ)
  · exact h r hr k
  have hf : frame (F := F) W g x r = frame W g y r := by
    funext j
    by_cases hj : j.val < W
    · have ht : r + j.val < 0 ∨ (N : ℤ) ≤ r + j.val := by
        simp only [Set.mem_Ico] at hr
        omega
      simp only [frame, if_pos hj, hx _ ht, hy _ ht]
    · simp [frame,hj]
  rw [hf]

private theorem no_wrap {W d : ℕ} (hF : 2 * W - 1 ≤ F) (hd : d < W) {j : ZMod F} (hj : j.val < W) :
    (j+(d : ZMod F)).val = j.val+d := by
  have hlt : j.val+d < F := by omega
  rw [← ZMod.natCast_zmod_val j, ← Nat.cast_add, ZMod.val_natCast_of_lt hlt,
    ZMod.val_natCast_of_lt (ZMod.val_lt j)]

/-- For a nonnegative lag in the window, the power data give this unaliased linear row
in the unknown sample products. -/
theorem lag_row {W d : ℕ} (hF : 2 * W - 1 ≤ F) (hd : d < W) (g : ℕ → ℝ) (x : ℤ →₀ ℝ) (r : ℤ) :
    correlation (frame (F := F) W g x r) (d : ZMod F) =
      ↑(∑ j : ZMod F, if j.val+d < W then
        g j.val*g (j.val+d)*(x (r+j.val)*x (r+j.val+d)) else 0 : ℝ) := by
  rw [correlation, Complex.ofReal_sum]
  apply sum_congr rfl
  intro j _
  by_cases hj : j.val < W
  · simp only [frame, if_pos hj, no_wrap hF hd hj, Nat.cast_add]
    by_cases hl : j.val+d < W <;> simp [hl, add_assoc, mul_comm, mul_left_comm, mul_assoc]
  · have hl : ¬ j.val+d < W := by omega
    simp [frame, hj, hl]

/-- The last nonzero taper product isolates the first differing sample product in a finite
convolution. The taper may vanish at either end. -/
theorem lag_products {W d : ℕ} (hF : 2 * W - 1 ≤ F) (hd : d < W) (g : ℕ → ℝ)
    (hg : ∃ s, s + d < W ∧ g s * g (s + d) ≠ 0) {x y : ℤ →₀ ℝ} (h : SamePowers F W g x y) :
    ∀ t : ℤ, x t*x (t+d) = y t*y (t+d) := by
  classical
  let D := fun t : ℤ ↦ x t*x (t+d) - y t*y (t+d)
  let S := (x.support ∪ y.support).filter fun t ↦ D t ≠ 0
  have hmem (t) : t ∈ S ↔ D t ≠ 0 := by
    simp only [S, mem_filter, mem_union, Finsupp.mem_support_iff, and_iff_right_iff_imp]
    intro hn
    by_contra! hz
    simp [D,hz] at hn
  by_contra! hn
  obtain ⟨u,hu⟩ := hn
  have hnS : S.Nonempty := ⟨u,(hmem u).mpr (sub_ne_zero.mpr hu)⟩
  let t := S.min' hnS
  let K := (range W).filter fun s ↦ s + d < W ∧ g s*g (s + d) ≠ 0
  obtain ⟨s, hs, hgs⟩ := hg
  have hK : K.Nonempty := ⟨s, mem_filter.mpr ⟨mem_range.mpr (by omega), hs, hgs⟩⟩
  let b := K.max' hK
  have hb := (mem_filter.mp (max'_mem K hK)).2
  let r : ℤ := t-b
  have hbd : b+d < W := hb.1
  have hbF : b < F := by omega
  have hxrow := congrFun ((powers_iff_correlation _ _).mp (h r)) (d : ZMod F)
  rw [lag_row hF hd, lag_row hF hd, Complex.ofReal_inj] at hxrow
  have hrow : (∑ j : ZMod F, if j.val+d < W then g j.val*g (j.val+d)*D (r+j.val) else 0) = 0 := by
    rw [← sub_eq_zero, ← sum_sub_distrib] at hxrow
    convert hxrow using 1
    apply sum_congr rfl
    intro j _
    by_cases hj : j.val+d < W <;> simp [hj,D,mul_sub]
  rw [Finset.sum_eq_single (b : ZMod F)] at hrow
  · simp only [ZMod.val_natCast_of_lt hbF, if_pos hbd] at hrow
    have ht : r+(b:ℤ) = t := by dsimp [r]; ring
    rw [ht] at hrow
    exact ((hmem t).mp (Finset.min'_mem S hnS)) ((mul_eq_zero.mp hrow).resolve_left hb.2)
  · intro j _ hj
    by_cases hl : j.val+d < W
    · have hne : j.val ≠ b := fun he ↦ hj (by rw [← he, ZMod.natCast_zmod_val])
      by_cases hzK : g j.val*g (j.val+d) = 0
      · simp [hl, hzK]
      have hjb := le_max' K j.val (mem_filter.mpr ⟨mem_range.mpr (by omega), hl, hzK⟩)
      have hlt : r+(j.val:ℤ) < t := by change j.val ≤ b at hjb; dsimp [r]; omega
      have hz : D (r+j.val) = 0 := by
        by_contra hnD
        exact (not_le.mpr hlt) (Finset.min'_le S _ ((hmem _).mpr hnD))
      simp only [if_pos hl, hz, mul_zero]
    · simp [hl]
  · simp

/-- The scale can include the unit-energy normalization. -/
noncomputable def periodicHann (scale : ℝ) (W s : ℕ) : ℝ :=
  scale * (1 - Real.cos (2*Real.pi*s/W)) / 2

theorem hann_zero (scale : ℝ) (W : ℕ) : periodicHann scale W 0 = 0 := by simp [periodicHann]

theorem hann_pos {scale : ℝ} (ha : 0 < scale) {W s : ℕ} (hs : 0 < s) (hsW : s < W) :
    0 < periodicHann scale W s := by
  have hW : 0 < (W : ℝ) := by exact_mod_cast (by omega : 0 < W)
  have hs' : 0 < (s : ℝ) := by exact_mod_cast hs
  have hsW' : (s : ℝ) < W := by exact_mod_cast hsW
  have hlt : 2*Real.pi*s/W < 2*Real.pi := (div_lt_iff₀ hW).mpr
    (mul_lt_mul_of_pos_left hsW' (by positivity))
  have ht : 0 < 2*Real.pi*s/W := by positivity
  have hcos := lt_of_le_of_ne (Real.cos_le_one _) (fun he ↦ ht.ne'
    ((Real.cos_eq_one_iff_of_lt_of_lt (by linarith [Real.pi_pos]) hlt).mp he))
  exact div_pos (mul_pos ha (sub_pos.mpr hcos)) (by norm_num)

theorem hann_lags {W : ℕ} (hW : 3 ≤ W) (hF : 2 * W - 1 ≤ F) {scale : ℝ} (hs : 0 < scale)
    {x y : ℤ →₀ ℝ} (h : SamePowers F W (periodicHann scale W) x y) {d : ℕ} (hd : d ≤ W - 2) :
    ∀ t : ℤ, x t*x (t+d) = y t*y (t+d) :=
  lag_products hF (by omega) _ ⟨W-1-d, by omega,
    mul_ne_zero (hann_pos hs (by omega) (by omega)).ne' (hann_pos hs (by omega) (by omega)).ne'⟩ h

def supportGraph (W : ℕ) (x : ℤ →₀ ℝ) : SimpleGraph {t : ℤ // x t ≠ 0} where
  Adj u v := u ≠ v ∧ u.val ≤ v.val+(W - 2 : ℕ) ∧ v.val ≤ u.val+(W - 2 : ℕ)
  symm := ⟨fun _ _ h ↦ ⟨h.1.symm,h.2.2,h.2.1⟩⟩
  loopless := ⟨fun _ h ↦ h.1 rfl⟩

def ComponentSigns (W : ℕ) (x y : ℤ →₀ ℝ) : Prop :=
  (∀ t, x t = 0 → y t = 0) ∧
    ∃ σ : (supportGraph W x).ConnectedComponent → ℝ, (∀ C, σ C = 1 ∨ σ C = -1) ∧
      ∀ v : {t : ℤ // x t ≠ 0}, y v = σ ((supportGraph W x).connectedComponentMk v) * x v

theorem signs_of_lags {W : ℕ} {x y : ℤ →₀ ℝ}
    (h : ∀ d : ℕ, d ≤ W - 2 → ∀ t : ℤ, x t*x (t+d) = y t*y (t+d)) : ComponentSigns W x y := by
  have diag (t : ℤ) := h 0 (Nat.zero_le _) t
  simp only [Nat.cast_zero, add_zero] at diag
  have hs (t : ℤ) : y t = x t ∨ y t = -x t := by
    have he : (y t - x t)*(y t + x t) = 0 := by nlinarith [diag t]
    rcases mul_eq_zero.mp he with he | he <;> [left; right] <;> linarith
  let ratio := fun v : {t : ℤ // x t ≠ 0} ↦ y v / x v
  have hv (v) : y v.val = ratio v * x v.val := (div_mul_cancel₀ _ v.property).symm
  have hpm (v) : ratio v = 1 ∨ ratio v = -1 := by
    rcases hs v with h | h <;> simp [ratio,h,v.property]
  have hp (u v : ℤ) (huv : u ≤ v) (hd : v ≤ u+(W - 2 : ℕ)) :
      x u*x v = y u*y v := by
    have he := h (v-u).toNat (by omega) u
    simpa only [show u+((v-u).toNat : ℤ) = v by omega] using he
  have hed {u v} (ha : (supportGraph W x).Adj u v) : ratio u = ratio v := by
    obtain ⟨_,h1,h2⟩ := ha
    have he : x u*x v = y u*y v := by
      rcases le_total u.val v.val with huv | hvu
      · exact hp _ _ huv h2
      · simpa [mul_comm] using hp _ _ hvu h1
    have hh : ratio u*ratio v = 1 := by
      apply mul_right_cancel₀ (mul_ne_zero u.property v.property)
      calc ratio u * ratio v * (x u * x v) = (ratio u * x u) * (ratio v * x v) := by ring
        _ = y u*y v := by rw [← hv u, ← hv v]
        _ = 1*(x u*x v) := by rw [one_mul, he]
    rcases hpm u with hu | hu <;> rcases hpm v with hv | hv <;> simp_all
  have hw {u v} (p : (supportGraph W x).Walk u v) : ratio u = ratio v := by
    induction p with
    | nil => rfl
    | cons he _ ih => exact (hed he).trans ih
  let σ := SimpleGraph.ConnectedComponent.lift (G := supportGraph W x) ratio fun _ _ p _ ↦ hw p
  refine ⟨fun t ht ↦ ?_, σ, ?_, hv⟩
  · rcases hs t with he | he <;> simp [he,ht]
  · intro C
    induction C using SimpleGraph.ConnectedComponent.ind with
    | h v => exact hpm v

theorem powers_of_signs {W : ℕ} {g : ℕ → ℝ} (hg : g 0 = 0) {x y : ℤ →₀ ℝ}
    (h : ComponentSigns W x y) : SamePowers F W g x y := by
  classical
  obtain ⟨hz,σ,hσ,hy⟩ := h
  intro r k
  let live := fun j : ZMod F ↦ 0 < j.val ∧ j.val < W ∧ x (r+j.val) ≠ 0
  have hdead (j : ZMod F) (hj : ¬live j) :
      frame W g x r j = 0 ∧ frame W g y r j = 0 := by
    by_cases hW : j.val < W
    · by_cases hzero : j.val = 0
      · simp [frame,hzero,hg]
      · have hx : x (r+j.val) = 0 := by
          by_contra hne
          exact hj ⟨Nat.pos_of_ne_zero hzero,hW,hne⟩
        simp only [frame, if_pos hW, hx, hz _ hx, mul_zero, Complex.ofReal_zero, and_self]
    · simp [frame,hW]
  by_cases hn : ∃ j, live j
  · obtain ⟨j,hj⟩ := hn
    let u : {t : ℤ // x t ≠ 0} := ⟨r+j.val,hj.2.2⟩
    let C := (supportGraph W x).connectedComponentMk u
    have hc {l : ZMod F} (hl : live l) :
        (supportGraph W x).connectedComponentMk ⟨r+l.val,hl.2.2⟩ = C := by
      by_cases he : l.val = j.val
      · congr 1
        exact Subtype.ext (by simp [u,he])
      · apply SimpleGraph.ConnectedComponent.connectedComponentMk_eq_of_adj
        refine ⟨?_,?_,?_⟩
        · intro hh
          apply he
          have eq := congrArg Subtype.val hh
          dsimp [u] at eq
          omega
        · dsimp [u]; obtain ⟨_,h₁,_⟩ := hl; obtain ⟨h₂,_,_⟩ := hj; omega
        · dsimp [u]; obtain ⟨h₁,_,_⟩ := hl; obtain ⟨_,h₂,_⟩ := hj; omega
    have hf : frame W g y r = (σ C : ℂ) • (frame (F := F) W g x r) := by
      funext l
      by_cases hl : live l
      · have he := hy ⟨r+l.val,hl.2.2⟩
        rw [hc hl] at he
        simp only [frame, if_pos hl.2.1, Pi.smul_apply, smul_eq_mul, he, Complex.ofReal_mul]
        ring
      · simp [(hdead l hl).1,(hdead l hl).2]
    rw [hf, map_smul]
    rcases hσ C with he | he <;> simp [he]
  · push Not at hn
    have he : frame (F := F) W g x r = frame W g y r := funext fun j ↦
      (hdead j (hn j)).1.trans (hdead j (hn j)).2.symm
    rw [he]

/-- Unit hop, all translations, F ≥ 2W−1, all F bins, finite real waveforms.
One independent sign per support component is necessary and sufficient. -/
theorem hann_recovery {W : ℕ} (hW : 3 ≤ W) (hF : 2 * W - 1 ≤ F) {scale : ℝ} (hs : 0 < scale)
    (x y : ℤ →₀ ℝ) : SamePowers F W (periodicHann scale W) x y ↔ ComponentSigns W x y :=
  ⟨fun h ↦ signs_of_lags fun _ hd ↦ hann_lags hW hF hs h hd,
    powers_of_signs (hann_zero _ _)⟩

/-- Length four can see the sign across a one-sample gap with seven retained frequencies. -/
theorem hann_sees_gap :
    ¬ SamePowers 7 4 (periodicHann 1 4)
      (Finsupp.single (0:ℤ) 1 + Finsupp.single 2 1)
      (Finsupp.single (0:ℤ) 1 - Finsupp.single 2 1) := by
  intro h
  have he := hann_lags (by decide : 3 ≤ 4) (by decide : 2*4-1 ≤ 7) zero_lt_one h
    (d := 2) (by decide) 0
  norm_num [Finsupp.single_apply] at he

#print axioms dft_correlation
#print axioms powers_iff_correlation
#print axioms finite_starts
#print axioms lag_row
#print axioms lag_products
#print axioms hann_zero
#print axioms hann_pos
#print axioms hann_lags
#print axioms signs_of_lags
#print axioms powers_of_signs
#print axioms hann_recovery
#print axioms hann_sees_gap

end TimeFrequency.FullGrid
