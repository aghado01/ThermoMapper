import Mathlib.Algebra.Order.Floor.Semiring
import Mathlib.Analysis.Normed.Module.Basic
import Mathlib.Topology.MetricSpace.CoveringNumbers
import Mathlib.Tactic.FieldSimp
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.Ring

/-!
# BGP-4, UV-4a — covering a scaled family and an interval

The entropy step (61′) of `Protolemmata/BARS/bgp-unknown-variance-lemmas.md` (UV-4) covers the
amplitude-scaled sieve `{s • g | s ∈ [a, b], g ∈ B}` by products of two covers, using
`s•g − s′•g′ = s•(g − g′) + (s − s′)•g′`:

* `isCover_image2_smul`: covers of `B` at radius `ε₁` and of `[a, b]` at radius `ε₂`, with the
  centres of the first bounded by `R`, give a cover of the scaled family at any `ε ≥ b ε₁ + R ε₂`.
* `coveringNumber_image2_smul_le`: so the internal covering numbers multiply; with the brief's
  radii `ε/(2b)` and `ε/(2R)` this is `coveringNumber_image2_smul_le_half`.
  `externalCoveringNumber_image2_smul_le` is the external form, through the standard factor-2
  comparison between the two.
* `exists_isCover_Icc` and `coveringNumber_Icc_le`: `[a, b]` at radius `ε > 0` needs at most
  `⌈(b − a)/(2ε)⌉ + 1` centres, the points `a + 2εk` clipped at `b`.

Internal covers are the natural choice for the product: centres in `B` are bounded by `R`, and the
product of internal covers lies in the scaled family, so it is internal again. The space is any
seminormed space over `ℝ`, which takes the GP sieve's sup-norm space.
-/

open Metric Set
open scoped NNReal ENNReal

namespace BGP

variable {E : Type*} [SeminormedAddCommGroup E] [NormedSpace ℝ E]

/-- The cardinality of `image2 f s t` is at most the product of the cardinalities. -/
theorem encard_image2_le {α β γ : Type*} (f : α → β → γ) (s : Set α) (t : Set β) :
    (image2 f s t).encard ≤ s.encard * t.encard := by
  rw [← image_prod, ← encard_prod]
  exact encard_image_le _ _

/-- **UV-4a (the cover).** If `N₁` covers `B` at radius `ε₁`, its points have norm at most `R`,
and `N₂` covers `[a, b]` (with `0 ≤ a`) at radius `ε₂`, then `{s′ • g′}` covers the scaled family
`{s • g | s ∈ [a, b], g ∈ B}` at any radius `ε ≥ b ε₁ + R ε₂`: by
`s•g − s′•g′ = s•(g − g′) + (s − s′)•g′`. -/
theorem isCover_image2_smul {a b R : ℝ} {B N₁ : Set E} {N₂ : Set ℝ} {ε₁ ε₂ ε : ℝ≥0}
    (ha : 0 ≤ a) (hN₁R : ∀ g ∈ N₁, ‖g‖ ≤ R) (hε : b * ε₁ + R * ε₂ ≤ ε)
    (h₁ : IsCover ε₁ B N₁) (h₂ : IsCover ε₂ (Icc a b) N₂) :
    IsCover ε (image2 (· • ·) (Icc a b) B) (image2 (· • ·) N₂ N₁) := by
  rw [isCover_iff_subset_iUnion_closedBall] at h₁ h₂ ⊢
  rintro _ ⟨s, hs, g, hg, rfl⟩
  obtain ⟨g', hg'N, hgg'⟩ := mem_iUnion₂.1 (h₁ hg)
  obtain ⟨s', hs'N, hss'⟩ := mem_iUnion₂.1 (h₂ hs)
  rw [mem_closedBall, dist_eq_norm] at hgg'
  rw [mem_closedBall, Real.dist_eq] at hss'
  refine mem_iUnion₂.2 ⟨s' • g', mem_image2_of_mem hs'N hg'N, ?_⟩
  rw [mem_closedBall, dist_eq_norm]
  obtain ⟨hsa, hsb⟩ := hs
  have hs0 : 0 ≤ s := ha.trans hsa
  have hsplit : s • g - s' • g' = s • (g - g') + (s - s') • g' := by
    rw [smul_sub, sub_smul]
    abel
  have h1 : ‖s • (g - g')‖ ≤ b * ε₁ := by
    rw [norm_smul, Real.norm_of_nonneg hs0]
    exact mul_le_mul hsb hgg' (norm_nonneg _) (hs0.trans hsb)
  have h2 : ‖(s - s') • g'‖ ≤ R * ε₂ := by
    rw [norm_smul, Real.norm_eq_abs, mul_comm]
    exact mul_le_mul (hN₁R g' hg'N) hss' (abs_nonneg _) ((norm_nonneg _).trans (hN₁R g' hg'N))
  calc ‖s • g - s' • g'‖ = ‖s • (g - g') + (s - s') • g'‖ := by rw [hsplit]
    _ ≤ ‖s • (g - g')‖ + ‖(s - s') • g'‖ := norm_add_le _ _
    _ ≤ b * ε₁ + R * ε₂ := add_le_add h1 h2
    _ ≤ ε := hε

/-- **UV-4a (internal covering numbers).** For `B` bounded by `R` and `0 ≤ a`, the scaled family
has `N(ε) ≤ N(ε₂, [a, b]) · N(ε₁, B)` whenever `b ε₁ + R ε₂ ≤ ε`. -/
theorem coveringNumber_image2_smul_le {a b R : ℝ} {B : Set E} {ε₁ ε₂ ε : ℝ≥0}
    (ha : 0 ≤ a) (hR : ∀ g ∈ B, ‖g‖ ≤ R) (hε : b * ε₁ + R * ε₂ ≤ ε) :
    coveringNumber ε (image2 (· • ·) (Icc a b) B)
      ≤ coveringNumber ε₂ (Icc a b) * coveringNumber ε₁ B := by
  rcases eq_empty_or_nonempty B with rfl | hB
  · simp
  rcases eq_empty_or_nonempty (Icc a b) with hI | hI
  · rw [hI]
    simp
  have hB0 : coveringNumber ε₁ B ≠ 0 := (coveringNumber_pos_iff.2 hB).ne'
  have hI0 : coveringNumber ε₂ (Icc a b) ≠ 0 := (coveringNumber_pos_iff.2 hI).ne'
  by_cases h1 : coveringNumber ε₁ B = ⊤
  · rw [h1, ENat.mul_top hI0]
    exact le_top
  by_cases h2 : coveringNumber ε₂ (Icc a b) = ⊤
  · rw [h2, ENat.top_mul hB0]
    exact le_top
  obtain ⟨C₁, hC₁B, -, hC₁, hC₁card⟩ := exists_set_encard_eq_coveringNumber h1
  obtain ⟨C₂, hC₂I, -, hC₂, hC₂card⟩ := exists_set_encard_eq_coveringNumber h2
  calc coveringNumber ε (image2 (· • ·) (Icc a b) B)
      ≤ (image2 (· • ·) C₂ C₁).encard :=
        (isCover_image2_smul ha (fun g hg => hR g (hC₁B hg)) hε hC₁ hC₂).coveringNumber_le_encard
          (image2_subset hC₂I hC₁B)
    _ ≤ C₂.encard * C₁.encard := encard_image2_le _ _ _
    _ = coveringNumber ε₂ (Icc a b) * coveringNumber ε₁ B := by rw [hC₁card, hC₂card]

/-- UV-4a with the radii of the brief: covers of `B` at `ε/(2b)` and of `[a, b]` at `ε/(2R)`. -/
theorem coveringNumber_image2_smul_le_half {a b R : ℝ} {B : Set E} (ha : 0 ≤ a) (hb : 0 < b)
    (hR0 : 0 < R) (hR : ∀ g ∈ B, ‖g‖ ≤ R) (ε : ℝ≥0) :
    coveringNumber ε (image2 (· • ·) (Icc a b) B)
      ≤ coveringNumber (ε / (2 * R.toNNReal)) (Icc a b)
        * coveringNumber (ε / (2 * b.toNNReal)) B := by
  refine coveringNumber_image2_smul_le ha hR (le_of_eq ?_)
  simp only [NNReal.coe_div, NNReal.coe_mul, NNReal.coe_ofNat, Real.coe_toNNReal _ hb.le,
    Real.coe_toNNReal _ hR0.le]
  field_simp
  ring

/-- UV-4a for external covering numbers, through `coveringNumber (2δ) ≤ externalCoveringNumber δ`:
the product needs the factor `2` in both radii. -/
theorem externalCoveringNumber_image2_smul_le {a b R : ℝ} {B : Set E} {ε₁ ε₂ ε : ℝ≥0}
    (ha : 0 ≤ a) (hR : ∀ g ∈ B, ‖g‖ ≤ R) (hε : b * (2 * ε₁) + R * (2 * ε₂) ≤ ε) :
    externalCoveringNumber ε (image2 (· • ·) (Icc a b) B)
      ≤ externalCoveringNumber ε₂ (Icc a b) * externalCoveringNumber ε₁ B := by
  have h := coveringNumber_image2_smul_le ha hR (ε₁ := 2 * ε₁) (ε₂ := 2 * ε₂)
    (by exact_mod_cast hε)
  refine (externalCoveringNumber_le_coveringNumber ε _).trans (h.trans ?_)
  exact mul_le_mul' (coveringNumber_two_mul_le_externalCoveringNumber _ _)
    (coveringNumber_two_mul_le_externalCoveringNumber _ _)

/-- **UV-4a (the interval).** `[a, b]` has an internal cover at radius `ε > 0` by at most
`⌈(b − a)/(2ε)⌉ + 1` points: `min (a + 2εk) b` for `k ≤ ⌈(b − a)/(2ε)⌉`. A point `x` is within
`ε` of `a + 2εk` for `k = ⌊(x − a)/(2ε) + 1/2⌋`, and clipping at `b` only moves the centre
closer. -/
theorem exists_isCover_Icc (a b : ℝ) {ε : ℝ≥0} (hε : 0 < ε) :
    ∃ N : Finset ℝ, (N : Set ℝ) ⊆ Icc a b ∧ IsCover ε (Icc a b) N ∧
      N.card ≤ ⌈(b - a) / (2 * ε)⌉₊ + 1 := by
  rcases lt_or_ge b a with hba | hab
  · refine ⟨∅, by simp, ?_, by simp⟩
    rw [Icc_eq_empty (not_le.2 hba)]
    exact IsCover.empty
  have hε' : (0 : ℝ) < ε := hε
  set K := ⌈(b - a) / (2 * ε)⌉₊ with hK
  refine ⟨(Finset.range (K + 1)).image (fun k : ℕ => min (a + 2 * ε * k) b), ?_, ?_, ?_⟩
  · intro c hc
    simp only [Finset.coe_image, Finset.coe_range, mem_image, mem_Iio] at hc
    obtain ⟨k, -, rfl⟩ := hc
    have hk : (0 : ℝ) ≤ 2 * ε * k := by positivity
    exact ⟨le_min (by linarith) hab, min_le_right _ _⟩
  · rw [isCover_iff_subset_iUnion_closedBall]
    rintro x ⟨hxa, hxb⟩
    set y := (x - a) / (2 * ε) + 1 / 2 with hy
    have hy0 : 0 ≤ y := by
      have : 0 ≤ (x - a) / (2 * ε) := div_nonneg (by linarith) (by positivity)
      linarith
    set k := ⌊y⌋₊ with hk
    have hky : (k : ℝ) ≤ y := Nat.floor_le hy0
    have hyk : y < k + 1 := Nat.lt_floor_add_one y
    have hkK : k < K + 1 := by
      have h1 : (x - a) / (2 * ε) ≤ (b - a) / (2 * ε) :=
        div_le_div_of_nonneg_right (by linarith) (by positivity)
      have h2 : (b - a) / (2 * ε) ≤ K := Nat.le_ceil _
      have : (k : ℝ) < K + 1 := by linarith
      exact_mod_cast this
    have hdist : |x - (a + 2 * ε * k)| ≤ ε := by
      have e1 : (x - a) / (2 * ε) * (2 * ε) = x - a := by field_simp
      rw [abs_le]
      constructor
      · have : y - 1 < k := by linarith
        nlinarith
      · nlinarith
    refine mem_iUnion₂.2 ⟨min (a + 2 * ε * k) b, ?_, ?_⟩
    · simp only [Finset.coe_image, Finset.coe_range, mem_image, mem_Iio]
      exact ⟨k, hkK, rfl⟩
    · rw [mem_closedBall, Real.dist_eq]
      rcases le_total (a + 2 * ε * k) b with h | h
      · rw [min_eq_left h]
        exact hdist
      · rw [min_eq_right h, abs_of_nonpos (by linarith)]
        have := (abs_le.1 hdist).1
        linarith
  · exact (Finset.card_image_le).trans (by simp)

/-- **UV-4a (the interval's count).** `N(ε, [a, b]) ≤ ⌈(b − a)/(2ε)⌉ + 1` for `ε > 0`. -/
theorem coveringNumber_Icc_le (a b : ℝ) {ε : ℝ≥0} (hε : 0 < ε) :
    coveringNumber ε (Icc a b) ≤ (⌈(b - a) / (2 * ε)⌉₊ + 1 : ℕ) := by
  obtain ⟨N, hNs, hN, hcard⟩ := exists_isCover_Icc a b hε
  calc coveringNumber ε (Icc a b) ≤ (N : Set ℝ).encard := hN.coveringNumber_le_encard hNs
    _ = N.card := encard_coe_eq_coe_finsetCard N
    _ ≤ _ := by exact_mod_cast hcard

end BGP
