import Enthymemata.BassStructureChain
import Enthymemata.BgpCovering
import Enthymemata.GuardedEstimator
import Mathlib.Analysis.Normed.Module.FiniteDimension
import Mathlib.Analysis.SpecialFunctions.Log.Basic
import Mathlib.LinearAlgebra.FiniteDimensional.Lemmas
import Mathlib.MeasureTheory.Measure.Lebesgue.EqHaar

/-!
# BASS — a union of subspaces' balls has a local cover free of `r`

The central claim of `Protolemmata/BASS/contraction-walk.md` about local entropy, checkable with no
paper. Under a structure `S` a regression function lies in a finite-dimensional subspace, so the
part of the structure's subspace within `r` of the truth is covered at radius `ξr/4` by a count that
depends on the dimension and `ξ` alone, with no bound on any coefficient.

* `card_le_of_pairwise_dist`: the volume argument, in any finite-dimensional normed space, so with
  no appeal to the norm being Euclidean. A separated set in a ball of radius `R` has at most
  `(1 + 2R/ε)^dim` points.
* `exists_cover_subspace_ball`, `exists_cover_iUnion_subspace_balls`,
  `exists_local_cover_iUnion`: the cover of one subspace's ball, of a finite union, and the union at
  radius `ξr/4` with count `|κ| (1 + 16/ξ)^M`. The centres lie in the covered set.
* `exists_local_cover_with_scale`: with the noise scale, in `‖f − f′‖ + |σ − σ′|`, through
  `BGP.exists_isCover_Icc`: the count is `|κ| (1 + 32/ξ)^M (⌈8/ξ⌉ + 1)`.
* `card_finsets_le`, `exists_dictionary_local_cover`, `log_dictionary_count_le`: for structures of
at
  most `M` terms from a finite dictionary, `log count ≤ M log(|D| + 1) + M log(1 + 32/ξ) +
  log(⌈8/ξ⌉ + 1)`.

Beside the cover, three small pieces the walks use:
* `exists_mTerm_approx`: Maurey's empirical method in the design norm, the `M`-term rate that *MARS
via
  LASSO*'s class gives with no further structure: `‖f − f_M‖ₙ² ≤ V² R² / M`, for `V` the
    coefficients'
  `ℓ¹` norm and `R` a bound on the dictionary. It reads neither the dictionary's size nor its
  dimension, which is what lets `p` enter the prior-mass step only through the structure's prior.
* `sum_powersetCard_prod_le`: the unordered `M`-term collections' weights sum to at most
  `(Σ c)^M / M!`, which bounds the structure prior's normalizer by a constant.
* `jeffreysProfile_mono`, `rSquaredOf_mono`, `jeffreys_superset_ge`: D03's Jeffreys marginal is
  increasing in `R²`, and `R²` cannot fall when the column space grows, so SC-5's superset step
    holds
  for it with no hypothesis on the profile. `invGammaProfile_mono` and `invGamma_superset_ge` do
  the same under BASS's default proper prior IG(c, c·s²). They sit here for the logarithm's
  import and the design norm.

The space is any normed space over `ℝ`. BASS's has the design vectors in `ℝ^N` and `‖·‖ₙ`, the
Euclidean norm scaled; the dictionary is the terms' columns, the intercept among them.
-/

open Metric Set MeasureTheory Module
open scoped NNReal ENNReal

namespace BASS

section Volume

variable {V : Type*} [NormedAddCommGroup V] [NormedSpace ℝ V] [FiniteDimensional ℝ V]

/-- **The volume argument.** Points more than `ε` apart in a closed ball of radius `R` number at
most
`(1 + 2R/ε)^dim`: the open balls of radius `ε/2` about them are disjoint and sit in a ball of radius
`R + ε/2`. -/
theorem card_le_of_pairwise_dist (C : Finset V) (a : V) {ε R : ℝ} (hε : 0 < ε) (hR : 0 ≤ R)
    (hsep : ∀ x ∈ C, ∀ y ∈ C, x ≠ y → ε < dist x y) (hC : ∀ x ∈ C, dist x a ≤ R) :
    (C.card : ℝ) ≤ (1 + 2 * R / ε) ^ finrank ℝ V := by
  borelize V
  let μ : Measure V := Measure.addHaar
  set M := finrank ℝ V with hM
  have hε2 : 0 < ε / 2 := by positivity
  have hRε : 0 < R + ε / 2 := by positivity
  set m := μ (ball (0 : V) 1) with hm
  have hm0 : m ≠ 0 := (measure_ball_pos μ 0 one_pos).ne'
  have hmT : m ≠ ⊤ := measure_ball_lt_top.ne
  have hdisj : (C : Set V).PairwiseDisjoint fun x => ball x (ε / 2) := by
    intro x hx y hy hxy
    apply ball_disjoint_ball
    have := hsep x hx y hy hxy
    linarith
  have hsub : (⋃ x ∈ C, ball x (ε / 2)) ⊆ ball a (R + ε / 2) := by
    intro z hz
    simp only [mem_iUnion, mem_ball, exists_prop] at hz
    obtain ⟨x, hx, hzx⟩ := hz
    rw [mem_ball]
    have := hC x hx
    have htri := dist_triangle z x a
    linarith
  have hsum : μ (⋃ x ∈ C, ball x (ε / 2)) = (C.card : ℝ≥0∞) * (ENNReal.ofReal ((ε / 2) ^ M) * m) :=
    by
    rw [measure_biUnion_finset hdisj fun _ _ => measurableSet_ball]
    have : ∀ x ∈ C, μ (ball x (ε / 2)) = ENNReal.ofReal ((ε / 2) ^ M) * m :=
      fun x _ => Measure.addHaar_ball_of_pos μ x hε2
    rw [Finset.sum_congr rfl this, Finset.sum_const, nsmul_eq_mul]
  have hle : (C.card : ℝ≥0∞) * (ENNReal.ofReal ((ε / 2) ^ M) * m) ≤
      ENNReal.ofReal ((R + ε / 2) ^ M) * m := by
    rw [← hsum, ← Measure.addHaar_ball_of_pos μ a hRε]
    exact measure_mono hsub
  have hcancel : (C.card : ℝ≥0∞) * ENNReal.ofReal ((ε / 2) ^ M) ≤
      ENNReal.ofReal ((R + ε / 2) ^ M) := by
    rw [← mul_assoc] at hle
    exact (ENNReal.mul_le_mul_iff_left hm0 hmT).1 hle
  have hreal : (C.card : ℝ) * (ε / 2) ^ M ≤ (R + ε / 2) ^ M := by
    rw [← ENNReal.ofReal_natCast, ← ENNReal.ofReal_mul (by positivity)] at hcancel
    exact (ENNReal.ofReal_le_ofReal_iff (by positivity)).1 hcancel
  have hbase : 1 + 2 * R / ε = (R + ε / 2) / (ε / 2) := by
    field_simp
    ring
  rw [hbase, div_pow, le_div_iff₀ (by positivity)]
  exact hreal

end Volume

section Subspace

variable {E : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E]

/-- A finite separated set in a subspace's ball of radius `r` obeys the volume bound in the
subspace's own dimension, whatever its centre: two of its points are within `2r`. -/
theorem card_le_of_subspace_ball (W : Submodule ℝ E) [FiniteDimensional ℝ W] {M : ℕ}
    (hW : finrank ℝ W ≤ M) (θ₀ : E) {ε r : ℝ} (hε : 0 < ε) (hr : 0 ≤ r) (C : Finset E)
    (hCW : ∀ x ∈ C, x ∈ W ∧ dist x θ₀ ≤ r) (hsep : ∀ x ∈ C, ∀ y ∈ C, x ≠ y → ε < dist x y) :
    (C.card : ℝ) ≤ (1 + 4 * r / ε) ^ M := by
  classical
  rcases C.eq_empty_or_nonempty with rfl | ⟨a, ha⟩
  · simp only [Finset.card_empty, Nat.cast_zero]
    positivity
  have hsub : (C.subtype (· ∈ W)).card = C.card := by
    rw [Finset.card_subtype, Finset.filter_true_of_mem fun x hx => (hCW x hx).1]
  have hbound := card_le_of_pairwise_dist (C.subtype (· ∈ W)) (⟨a, (hCW a ha).1⟩ : W) hε
    (R := 2 * r) (by positivity)
    (fun x hx y hy hxy => by
      rw [Subtype.dist_eq]
      exact hsep _ (Finset.mem_subtype.1 hx) _ (Finset.mem_subtype.1 hy)
        fun h => hxy (Subtype.ext h))
    (fun x hx => by
      rw [Subtype.dist_eq]
      have h₁ := (hCW _ (Finset.mem_subtype.1 hx)).2
      have h₂ := (hCW a ha).2
      have htri := dist_triangle_right (x : E) a θ₀
      linarith)
  rw [hsub] at hbound
  have hbase : 2 * (2 * r) / ε = 4 * r / ε := by ring
  rw [hbase] at hbound
  exact hbound.trans (pow_le_pow_right₀ (by
    have : 0 ≤ 4 * r / ε := by positivity
    linarith) hW)

/-- **One subspace.** The part of a finite-dimensional subspace within `r` of any point has a cover
at radius `ε`, with centres in that part, by at most `(1 + 4r/ε)^M` points. -/
theorem exists_cover_subspace_ball (W : Submodule ℝ E) [FiniteDimensional ℝ W] {M : ℕ}
    (hW : finrank ℝ W ≤ M) (θ₀ : E) {ε r : ℝ} (hε : 0 < ε) (hr : 0 ≤ r) :
    ∃ N : Finset E, (N : Set E) ⊆ (W : Set E) ∩ closedBall θ₀ r ∧
      (W : Set E) ∩ closedBall θ₀ r ⊆ ⋃ c ∈ N, closedBall c ε ∧
      (N.card : ℝ) ≤ (1 + 4 * r / ε) ^ M := by
  set S : Set E := (W : Set E) ∩ closedBall θ₀ r with hS
  set ε' : ℝ≥0 := ε.toNNReal with hε'
  have hcoe : (ε' : ℝ) = ε := Real.coe_toNNReal _ hε.le
  set B : ℕ := ⌊(1 + 4 * r / ε) ^ M⌋₊ with hB
  have hsepdist : ∀ C : Set E, Metric.IsSeparated (ε' : ℝ≥0∞) C →
      ∀ x ∈ C, ∀ y ∈ C, x ≠ y → ε < dist x y := by
    intro C hC x hx y hy hxy
    have h : (ε' : ℝ≥0∞) < edist x y := hC hx hy hxy
    rw [edist_dist] at h
    have h' : ENNReal.ofReal ε < ENNReal.ofReal (dist x y) := h
    exact (ENNReal.ofReal_lt_ofReal_iff (dist_pos.2 hxy)).1 h'
  have hfinite : ∀ t : Finset E, (t : Set E) ⊆ S → Metric.IsSeparated ε' (t : Set E) →
      t.card ≤ B := by
    intro t htS hsep
    have := card_le_of_subspace_ball W hW θ₀ hε hr t
      (fun x hx => by
        have := htS hx
        exact ⟨this.1, by simpa [mem_closedBall] using this.2⟩)
      (fun x hx y hy hxy => hsepdist _ hsep x hx y hy hxy)
    exact Nat.le_floor this
  have hpack : packingNumber ε' S ≤ B := by
    refine iSup₂_le fun C hCS => iSup_le fun hsep => ?_
    by_contra hlt
    rw [not_le] at hlt
    obtain ⟨t, htC, htcard⟩ := Set.exists_subset_encard_eq (Order.add_one_le_of_lt hlt)
    have htfin : t.Finite := Set.finite_of_encard_eq_coe htcard
    have hle := hfinite htfin.toFinset (by simpa using htC.trans hCS)
      (by rw [Set.Finite.coe_toFinset]; exact Metric.IsSeparated.subset htC hsep)
    have : (htfin.toFinset.card : ℕ∞) = (B : ℕ∞) + 1 := by
      rw [← htcard, htfin.encard_eq_coe_toFinset_card]
    have hc : htfin.toFinset.card = B + 1 := by exact_mod_cast this
    omega
  have hne : packingNumber ε' S ≠ ⊤ := ne_top_of_le_ne_top (ENat.coe_ne_top B) hpack
  have hcov := isCover_maximalSeparatedSet hne
  have hsubset : maximalSeparatedSet ε' S ⊆ S := maximalSeparatedSet_subset
  have hencard : (maximalSeparatedSet ε' S).encard ≤ B := by
    rw [encard_maximalSeparatedSet hne]
    exact hpack
  have hfin : (maximalSeparatedSet ε' S).Finite := by
    apply Set.finite_of_encard_le_coe hencard
  refine ⟨hfin.toFinset, by simpa using hsubset, ?_, ?_⟩
  · rw [isCover_iff_subset_iUnion_closedBall] at hcov
    intro x hx
    obtain ⟨c, hc, hxc⟩ := mem_iUnion₂.1 (hcov hx)
    refine mem_iUnion₂.2 ⟨c, by simpa using hc, ?_⟩
    rwa [hcoe] at hxc
  · have hle := hfinite hfin.toFinset (by simpa using hsubset)
      (by rw [Set.Finite.coe_toFinset]; exact isSeparated_maximalSeparatedSet)
    calc (hfin.toFinset.card : ℝ) ≤ (B : ℝ) := by exact_mod_cast hle
      _ ≤ (1 + 4 * r / ε) ^ M := Nat.floor_le (by positivity)

variable {κ : Type*} [Fintype κ]

/-- **A finite union.** -/
theorem exists_cover_iUnion_subspace_balls (V : κ → Submodule ℝ E) [∀ k, FiniteDimensional ℝ (V k)]
    {M : ℕ} (hV : ∀ k, finrank ℝ (V k) ≤ M) (θ₀ : E) {ε r : ℝ} (hε : 0 < ε) (hr : 0 ≤ r) :
    ∃ N : Finset E, (N : Set E) ⊆ ⋃ k, (V k : Set E) ∩ closedBall θ₀ r ∧
      (⋃ k, (V k : Set E) ∩ closedBall θ₀ r) ⊆ ⋃ c ∈ N, closedBall c ε ∧
      (N.card : ℝ) ≤ Fintype.card κ * (1 + 4 * r / ε) ^ M := by
  classical
  choose C hCsub hCcov hCcard using
    fun k => exists_cover_subspace_ball (V k) (hV k) θ₀ hε hr
  refine ⟨Finset.univ.biUnion C, ?_, ?_, ?_⟩
  · intro x hx
    obtain ⟨k, -, hk⟩ := Finset.mem_biUnion.1 (Finset.mem_coe.1 hx)
    exact mem_iUnion.2 ⟨k, hCsub k hk⟩
  · intro x hx
    obtain ⟨k, hk⟩ := mem_iUnion.1 hx
    obtain ⟨c, hc, hxc⟩ := mem_iUnion₂.1 (hCcov k hk)
    exact mem_iUnion₂.2 ⟨c, Finset.mem_biUnion.2 ⟨k, Finset.mem_univ _, hc⟩, hxc⟩
  · calc ((Finset.univ.biUnion C).card : ℝ) ≤ ((∑ k, (C k).card : ℕ) : ℝ) := by
          exact_mod_cast Finset.card_biUnion_le
      _ = ∑ k, ((C k).card : ℝ) := by push_cast; rfl
      _ ≤ ∑ _k : κ, (1 + 4 * r / ε) ^ M := Finset.sum_le_sum fun k _ => hCcard k
      _ = Fintype.card κ * (1 + 4 * r / ε) ^ M := by
          rw [Finset.sum_const, Finset.card_univ, nsmul_eq_mul]

/-- **The walk's claim about local entropy.** A union of `|κ|` subspaces of dimension at most `M`:
the part within `r` of any point is covered at radius `ξr/4` by at most `|κ| (1 + 16/ξ)^M` points of
itself. The count is free of `r` and of any bound on the coefficients. -/
theorem exists_local_cover_iUnion (V : κ → Submodule ℝ E) [∀ k, FiniteDimensional ℝ (V k)]
    {M : ℕ} (hV : ∀ k, finrank ℝ (V k) ≤ M) (θ₀ : E) {ξ r : ℝ} (hξ : 0 < ξ) (hr : 0 < r) :
    ∃ N : Finset E, (N : Set E) ⊆ ⋃ k, (V k : Set E) ∩ closedBall θ₀ r ∧
      (⋃ k, (V k : Set E) ∩ closedBall θ₀ r) ⊆ ⋃ c ∈ N, closedBall c (ξ * r / 4) ∧
      (N.card : ℝ) ≤ Fintype.card κ * (1 + 16 / ξ) ^ M := by
  have hbase : 4 * r / (ξ * r / 4) = 16 / ξ := by
    field_simp
    ring
  have := exists_cover_iUnion_subspace_balls V hV θ₀ (ε := ξ * r / 4) (by positivity) hr.le
  rwa [hbase] at this

/-- **With the noise scale.** In `‖f − f′‖ + |σ − σ′|`, the part of the union times `ℝ` within `r`
of
`(f₀, σ₀)` is covered at radius `ξr/4` by at most `|κ| (1 + 32/ξ)^M (⌈8/ξ⌉ + 1)` balls. So the cover
reads no norm bound on `f` and no compact range for `σ`. -/
theorem exists_local_cover_with_scale (V : κ → Submodule ℝ E) [∀ k, FiniteDimensional ℝ (V k)]
    {M : ℕ} (hV : ∀ k, finrank ℝ (V k) ≤ M) (f₀ : E) (σ₀ : ℝ) {ξ r : ℝ} (hξ : 0 < ξ) (hr : 0 < r) :
    ∃ N : Finset (E × ℝ),
      {θ : E × ℝ | (∃ k, θ.1 ∈ V k) ∧ ‖θ.1 - f₀‖ + |θ.2 - σ₀| ≤ r} ⊆
        ⋃ c ∈ N, {θ : E × ℝ | ‖θ.1 - c.1‖ + |θ.2 - c.2| ≤ ξ * r / 4} ∧
      (N.card : ℝ) ≤ Fintype.card κ * (1 + 32 / ξ) ^ M * ((⌈8 / ξ⌉₊ : ℝ) + 1) := by
  have hε : 0 < ξ * r / 8 := by positivity
  obtain ⟨Nf, -, hNfcov, hNfcard⟩ :=
    exists_cover_iUnion_subspace_balls V hV f₀ (ε := ξ * r / 8) hε hr.le
  set ε' : ℝ≥0 := (ξ * r / 8).toNNReal with hε'
  have hcoe : (ε' : ℝ) = ξ * r / 8 := Real.coe_toNNReal _ hε.le
  have hε'pos : 0 < ε' := Real.toNNReal_pos.2 hε
  obtain ⟨Nσ, -, hNσcov, hNσcard⟩ := BGP.exists_isCover_Icc (σ₀ - r) (σ₀ + r) hε'pos
  rw [isCover_iff_subset_iUnion_closedBall] at hNσcov
  classical
  refine ⟨Nf ×ˢ Nσ, ?_, ?_⟩
  · rintro ⟨f, σ⟩ ⟨⟨k, hfk⟩, hle⟩
    simp only at hfk hle
    have hfnorm : ‖f - f₀‖ ≤ r := by
      have := abs_nonneg (σ - σ₀)
      linarith
    have hσabs : |σ - σ₀| ≤ r := by
      have := norm_nonneg (f - f₀)
      linarith
    have hfmem : f ∈ ⋃ k, (V k : Set E) ∩ closedBall f₀ r :=
      mem_iUnion.2 ⟨k, hfk, by rwa [mem_closedBall, dist_eq_norm]⟩
    obtain ⟨c, hc, hfc⟩ := mem_iUnion₂.1 (hNfcov hfmem)
    have hσmem : σ ∈ Icc (σ₀ - r) (σ₀ + r) := by
      rw [abs_le] at hσabs
      exact ⟨by linarith [hσabs.1], by linarith [hσabs.2]⟩
    obtain ⟨s, hs, hσs⟩ := mem_iUnion₂.1 (hNσcov hσmem)
    refine mem_iUnion₂.2 ⟨(c, s), Finset.mem_product.2 ⟨hc, hs⟩, ?_⟩
    rw [mem_closedBall, dist_eq_norm] at hfc
    rw [mem_closedBall, Real.dist_eq, hcoe] at hσs
    change ‖f - c‖ + |σ - s| ≤ ξ * r / 4
    linarith
  · have hratio₁ : 4 * r / (ξ * r / 8) = 32 / ξ := by
      field_simp
      ring
    have hratio₂ : (σ₀ + r - (σ₀ - r)) / (2 * (ε' : ℝ)) = 8 / ξ := by
      rw [hcoe]
      field_simp
      ring
    rw [hratio₁] at hNfcard
    rw [hratio₂] at hNσcard
    rw [Finset.card_product]
    push_cast
    have h₂ : (Nσ.card : ℝ) ≤ (⌈8 / ξ⌉₊ : ℝ) + 1 := by exact_mod_cast hNσcard
    exact mul_le_mul hNfcard h₂ (Nat.cast_nonneg _) (by positivity)

end Subspace

/-! ### Structures from a finite dictionary -/

section Dictionary

variable {D : Type*} [Fintype D]

/-- At most `(|D| + 1)^M` subsets of `D` have at most `M` elements: each is the image of a map from
`M` slots to `D` or nothing. -/
theorem card_finsets_le (M : ℕ) :
    Fintype.card {S : Finset D // S.card ≤ M} ≤ (Fintype.card D + 1) ^ M := by
  classical
  have hle : ∀ g : Fin M → Option D, ((Finset.univ.image g).eraseNone).card ≤ M := by
    intro g
    calc ((Finset.univ.image g).eraseNone).card ≤ (Finset.univ.image g).card :=
          Finset.card_eraseNone_le _
      _ ≤ (Finset.univ : Finset (Fin M)).card := Finset.card_image_le
      _ = M := by simp
  have hsurj : Function.Surjective
      (fun g : Fin M → Option D => (⟨(Finset.univ.image g).eraseNone, hle g⟩ :
        {S : Finset D // S.card ≤ M})) := by
    rintro ⟨S, hS⟩
    refine ⟨fun i => S.toList[(i : ℕ)]?, Subtype.ext ?_⟩
    ext x
    simp only [Finset.mem_eraseNone, Finset.mem_image, Finset.mem_univ, true_and]
    constructor
    · rintro ⟨i, hi⟩
      exact Finset.mem_toList.1 (List.mem_of_getElem? hi)
    · intro hx
      obtain ⟨i, hi⟩ := List.mem_iff_getElem?.1 (Finset.mem_toList.2 hx)
      have hlt : i < S.toList.length := (List.getElem?_eq_some_iff.1 hi).1
      have hltM : i < M := by
        rw [Finset.length_toList] at hlt
        omega
      exact ⟨⟨i, hltM⟩, hi⟩
  calc Fintype.card {S : Finset D // S.card ≤ M}
      ≤ Fintype.card (Fin M → Option D) := Fintype.card_le_of_surjective _ hsurj
    _ = (Fintype.card D + 1) ^ M := by simp

variable {E : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E] [DecidableEq E]

/-- The column space of a structure, its terms' span. -/
def structureSpace {D : Type*} (b : D → E) (S : Finset D) : Submodule ℝ E :=
  Submodule.span ℝ ((S.image b : Finset E) : Set E)

instance (b : D → E) (S : Finset D) : FiniteDimensional ℝ (structureSpace b S) :=
  FiniteDimensional.span_finset ℝ (S.image b)

omit [Fintype D] in
theorem finrank_structureSpace_le (b : D → E) (S : Finset D) :
    finrank ℝ (structureSpace b S) ≤ S.card :=
  (finrank_span_finset_le_card (S.image b)).trans Finset.card_image_le

/-- **Structures of at most `M` terms, with the noise scale.** The local cover in
`‖f − f′‖ + |σ − σ′|`, its count free of `r`, of any norm bound on the coefficients, and of any
range for `σ`. -/
theorem exists_dictionary_local_cover (b : D → E) (M : ℕ) (f₀ : E) (σ₀ : ℝ) {ξ r : ℝ}
    (hξ : 0 < ξ) (hr : 0 < r) :
    ∃ N : Finset (E × ℝ),
      {θ : E × ℝ | (∃ S : Finset D, S.card ≤ M ∧ θ.1 ∈ structureSpace b S) ∧
          ‖θ.1 - f₀‖ + |θ.2 - σ₀| ≤ r} ⊆
        ⋃ c ∈ N, {θ : E × ℝ | ‖θ.1 - c.1‖ + |θ.2 - c.2| ≤ ξ * r / 4} ∧
      (N.card : ℝ) ≤ ((Fintype.card D : ℝ) + 1) ^ M * (1 + 32 / ξ) ^ M * ((⌈8 / ξ⌉₊ : ℝ) + 1) := by
  obtain ⟨N, hcov, hcard⟩ := exists_local_cover_with_scale
    (fun S : {S : Finset D // S.card ≤ M} => structureSpace b S.1) (M := M)
    (fun S => (finrank_structureSpace_le b S.1).trans S.2) f₀ σ₀ hξ hr
  refine ⟨N, ?_, hcard.trans ?_⟩
  · rintro θ ⟨⟨S, hS, hθ⟩, hle⟩
    exact hcov ⟨⟨⟨S, hS⟩, hθ⟩, hle⟩
  · have hk : (Fintype.card {S : Finset D // S.card ≤ M} : ℝ) ≤ ((Fintype.card D : ℝ) + 1) ^ M := by
      exact_mod_cast card_finsets_le (D := D) M
    have h32 : (0 : ℝ) ≤ (1 + 32 / ξ) ^ M := by positivity
    have h8 : (0 : ℝ) ≤ (⌈8 / ξ⌉₊ : ℝ) + 1 := by positivity
    exact mul_le_mul_of_nonneg_right (mul_le_mul_of_nonneg_right hk h32) h8

/-- **The log count.** `log N ≤ M log(|D| + 1) + M log(1 + 32/ξ) + log(⌈8/ξ⌉ + 1)`: the walk's
`M̄ₙ (K log(pN))`, once `|D|` is at most a polynomial in `pN` of degree `K`. -/
theorem log_dictionary_count_le (card : ℕ) (M : ℕ) {ξ : ℝ} (hξ : 0 < ξ) :
    Real.log (((card : ℝ) + 1) ^ M * (1 + 32 / ξ) ^ M * ((⌈8 / ξ⌉₊ : ℝ) + 1)) =
      M * Real.log ((card : ℝ) + 1) + M * Real.log (1 + 32 / ξ) +
        Real.log ((⌈8 / ξ⌉₊ : ℝ) + 1) := by
  have h₁ : ((card : ℝ) + 1) ≠ 0 := by positivity
  have h₂ : (1 + 32 / ξ) ≠ 0 := by positivity
  have h₃ : ((⌈8 / ξ⌉₊ : ℝ) + 1) ≠ 0 := by positivity
  rw [Real.log_mul (mul_ne_zero (pow_ne_zero _ h₁) (pow_ne_zero _ h₂)) h₃,
    Real.log_mul (pow_ne_zero _ h₁) (pow_ne_zero _ h₂), Real.log_pow, Real.log_pow]

end Dictionary

/-! ### The structure prior's normalizer -/

section Normalizer

variable {D : Type*} [Fintype D]

/-- Choosing an `(M+1)`-set and one of its members is choosing an `M`-set and a member outside it.
-/
private theorem sum_powersetCard_succ_swap [DecidableEq D] (c : D → ℝ) (M : ℕ) :
    ∑ T ∈ (Finset.univ : Finset D).powersetCard (M + 1), ∑ _t ∈ T, ∏ i ∈ T, c i =
      ∑ S ∈ (Finset.univ : Finset D).powersetCard M, ∑ t ∈ Sᶜ, c t * ∏ i ∈ S, c i := by
  rw [Finset.sum_sigma', Finset.sum_sigma']
  refine Finset.sum_nbij' (fun x => ⟨x.1.erase x.2, x.2⟩) (fun y => ⟨insert y.2 y.1, y.2⟩)
    ?_ ?_ ?_ ?_ ?_
  · rintro ⟨T, t⟩ hx
    obtain ⟨hT, ht⟩ := Finset.mem_sigma.1 hx
    have hcard := (Finset.mem_powersetCard.1 hT).2
    refine Finset.mem_sigma.2 ⟨Finset.mem_powersetCard.2 ⟨Finset.subset_univ _, ?_⟩, ?_⟩
    · simp only at ht ⊢
      rw [Finset.card_erase_of_mem ht, hcard, Nat.add_sub_cancel]
    · simp
  · rintro ⟨S, t⟩ hy
    obtain ⟨hS, ht⟩ := Finset.mem_sigma.1 hy
    have hcard := (Finset.mem_powersetCard.1 hS).2
    have htS : t ∉ S := Finset.mem_compl.1 ht
    refine Finset.mem_sigma.2 ⟨Finset.mem_powersetCard.2 ⟨Finset.subset_univ _, ?_⟩, ?_⟩
    · simp only
      rw [Finset.card_insert_of_notMem htS, hcard]
    · simp
  · rintro ⟨T, t⟩ hx
    have ht := (Finset.mem_sigma.1 hx).2
    simp only at ht
    simp [Finset.insert_erase ht]
  · rintro ⟨S, t⟩ hy
    have htS : t ∉ S := Finset.mem_compl.1 (Finset.mem_sigma.1 hy).2
    simp [Finset.erase_insert htS]
  · rintro ⟨T, t⟩ hx
    have ht := (Finset.mem_sigma.1 hx).2
    simp only at ht ⊢
    rw [Finset.mul_prod_erase T c ht]

/-- **The structure prior's normalizer.** `M!` times the sum over `M`-sets of the product of their
weights is at most `(Σ c)^M`. With `Σ c ≤ 1` the collections of `M` terms weigh at most `1/M!`, so
`W_F`, the count weight's sum over `F`, lies between 1 and `((h₂+1)/h₂)^{h₁}`. -/
theorem sum_powersetCard_prod_le (c : D → ℝ) (hc : ∀ t, 0 ≤ c t) (M : ℕ) :
    (M.factorial : ℝ) * ∑ S ∈ (Finset.univ : Finset D).powersetCard M, ∏ t ∈ S, c t ≤
      (∑ t, c t) ^ M := by
  classical
  induction M with
  | zero => simp
  | succ M ih =>
    have hL : ((M : ℝ) + 1) * ∑ T ∈ (Finset.univ : Finset D).powersetCard (M + 1), ∏ i ∈ T, c i =
        ∑ T ∈ (Finset.univ : Finset D).powersetCard (M + 1), ∑ _t ∈ T, ∏ i ∈ T, c i := by
      rw [Finset.mul_sum]
      refine Finset.sum_congr rfl fun T hT => ?_
      rw [Finset.sum_const, (Finset.mem_powersetCard.1 hT).2, nsmul_eq_mul]
      push_cast
      ring
    have hS : 0 ≤ ∑ t, c t := Finset.sum_nonneg fun t _ => hc t
    have hstep : ((M : ℝ) + 1) * ∑ T ∈ (Finset.univ : Finset D).powersetCard (M + 1), ∏ i ∈ T, c i ≤
        (∑ S ∈ (Finset.univ : Finset D).powersetCard M, ∏ i ∈ S, c i) * ∑ t, c t := by
      rw [hL, sum_powersetCard_succ_swap, Finset.sum_mul]
      refine Finset.sum_le_sum fun S _ => ?_
      have hprod : 0 ≤ ∏ i ∈ S, c i := Finset.prod_nonneg fun i _ => hc i
      rw [mul_comm (∏ i ∈ S, c i), Finset.sum_mul]
      exact Finset.sum_le_sum_of_subset_of_nonneg (Finset.subset_univ _)
        fun t _ _ => mul_nonneg (hc t) hprod
    rw [Nat.factorial_succ]
    push_cast
    calc ((M : ℝ) + 1) * (M.factorial : ℝ) *
          ∑ S ∈ (Finset.univ : Finset D).powersetCard (M + 1), ∏ t ∈ S, c t
        = (M.factorial : ℝ) * (((M : ℝ) + 1) *
            ∑ S ∈ (Finset.univ : Finset D).powersetCard (M + 1), ∏ t ∈ S, c t) := by ring
      _ ≤ (M.factorial : ℝ) * ((∑ S ∈ (Finset.univ : Finset D).powersetCard M, ∏ i ∈ S, c i) *
            ∑ t, c t) := mul_le_mul_of_nonneg_left hstep (Nat.cast_nonneg _)
      _ = ((M.factorial : ℝ) * ∑ S ∈ (Finset.univ : Finset D).powersetCard M, ∏ i ∈ S, c i) *
            ∑ t, c t := by ring
      _ ≤ (∑ t, c t) ^ M * ∑ t, c t := mul_le_mul_of_nonneg_right ih hS
      _ = (∑ t, c t) ^ (M + 1) := (pow_succ _ _).symm

end Normalizer

/-! ### An `M`-term approximation in the design norm -/

section MTerm

variable {D : Type*} [Fintype D]

/-- `‖v‖ₙ²`, on design vectors. -/
noncomputable def designNormSq {N : ℕ} (v : Fin N → ℝ) : ℝ := (∑ k, v k ^ 2) / (N : ℝ)

theorem designNormSq_nonneg {N : ℕ} (v : Fin N → ℝ) : 0 ≤ designNormSq v :=
  div_nonneg (Finset.sum_nonneg fun _ _ => sq_nonneg _) (Nat.cast_nonneg _)

/-- **Maurey's empirical method, in the design norm.** A combination of a dictionary with
coefficient
norm `V = Σ |a_d|` and atoms of size at most `R` has an `M`-term approximation with
`‖f − f_M‖ₙ² ≤ V² R² / M`: draw `M` atoms with probability `|a_d|/V`, and some draw beats the
expected squared error. The bound reads neither the dictionary's size nor the design. -/
theorem exists_mTerm_approx {N M : ℕ} (hM : 0 < M) (b : D → Fin N → ℝ) (a : D → ℝ) {R : ℝ}
    (hb : ∀ d, designNormSq (b d) ≤ R ^ 2) :
    ∃ S : Finset D, S.card ≤ M ∧ ∃ g ∈ Submodule.span ℝ (b '' (S : Set D)),
      designNormSq (fun k => ∑ d, a d * b d k - g k) ≤ (∑ d, |a d|) ^ 2 * R ^ 2 / (M : ℝ) := by
  classical
  set V := ∑ d, |a d| with hV
  have hV0 : 0 ≤ V := Finset.sum_nonneg fun d _ => abs_nonneg _
  have hMr : (0 : ℝ) < M := by exact_mod_cast hM
  rcases hV0.eq_or_lt with hVz | hVpos
  · have ha : ∀ d, a d = 0 := fun d => by
      have := (Finset.sum_eq_zero_iff_of_nonneg fun d _ => abs_nonneg (a d)).1 hVz.symm d
        (Finset.mem_univ d)
      exact abs_eq_zero.1 this
    refine ⟨∅, by simp, 0, Submodule.zero_mem _, ?_⟩
    simp [designNormSq, ha, ← hVz]
  -- The draws: atom `d` with probability `|a d|/V`, carrying the sign and scale `V`.
  set p : D → ℝ := fun d => |a d| / V with hp
  have hp0 : ∀ d, 0 ≤ p d := fun d => div_nonneg (abs_nonneg _) hV0
  have hp1 : ∑ d, p d = 1 := by
    simp only [hp]
    rw [← Finset.sum_div, ← hV, div_self hVpos.ne']
  set s : D → ℝ := fun d => if 0 ≤ a d then V else -V with hs
  have hps : ∀ d, p d * s d = a d := by
    intro d
    simp only [hp, hs]
    split_ifs with h
    · rw [abs_of_nonneg h]
      field_simp
    · rw [abs_of_neg (not_le.1 h)]
      field_simp
  have hs2 : ∀ d, s d ^ 2 = V ^ 2 := by
    intro d
    simp only [hs]
    split_ifs <;> ring
  set c : D → Fin N → ℝ := fun d k => s d * b d k with hc
  set f : Fin N → ℝ := fun k => ∑ d, a d * b d k with hf
  have hfc : ∀ k, integralOf p (fun d => c d k) = f k := by
    intro k
    simp only [integralOf, hc, hf]
    exact Finset.sum_congr rfl fun d _ => by rw [← hps d]; ring
  set h : Fin N → D → ℝ := fun k d => f k - c d k with hh
  have hh0 : ∀ k, integralOf p (h k) = 0 := by
    intro k
    have : integralOf p (h k) = f k * ∑ d, p d - integralOf p (fun d => c d k) := by
      simp only [integralOf, hh, mul_sub, Finset.sum_sub_distrib, Finset.mul_sum]
      congr 1
      exact Finset.sum_congr rfl fun _ _ => by ring
    rw [this, hp1, hfc]
    ring
  have hh2 : ∀ k, integralOf p (fun d => h k d ^ 2) ≤ ∑ d, p d * c d k ^ 2 := by
    intro k
    have : integralOf p (fun d => h k d ^ 2) =
        ∑ d, p d * c d k ^ 2 - f k ^ 2 := by
      have hexp : ∀ d, p d * h k d ^ 2 =
          p d * c d k ^ 2 - 2 * f k * (p d * c d k) + f k ^ 2 * p d := fun d => by
        simp only [hh]
        ring
      simp only [integralOf, hexp]
      rw [Finset.sum_add_distrib, Finset.sum_sub_distrib, ← Finset.mul_sum, ← Finset.mul_sum, hp1]
      have := hfc k
      simp only [integralOf] at this
      rw [this]
      ring
    rw [this]
    nlinarith [sq_nonneg (f k)]
  -- The error of a draw, and its expectation.
  set Φ : (Fin M → D) → ℝ :=
    fun x => designNormSq (fun k => f k - (∑ i, c (x i) k) / (M : ℝ)) with hΦ
  have hrewrite : ∀ x : Fin M → D, Φ x =
      ∑ k, ((∑ i, h k (x i)) ^ 2 / ((M : ℝ) ^ 2 * (N : ℝ))) := by
    intro x
    simp only [hΦ, designNormSq]
    rw [Finset.sum_div]
    refine Finset.sum_congr rfl fun k _ => ?_
    have hsum : ∑ i : Fin M, h k (x i) = (M : ℝ) * f k - ∑ i, c (x i) k := by
      simp only [hh, Finset.sum_sub_distrib, Finset.sum_const, Finset.card_univ, Fintype.card_fin,
        nsmul_eq_mul]
    have hdiff : f k - (∑ i, c (x i) k) / (M : ℝ) = ((M : ℝ) * f k - ∑ i, c (x i) k) / (M : ℝ) := by
      field_simp
    rw [hsum, hdiff, div_pow, div_div]
  have hexpect : drawsExpect p M Φ ≤ V ^ 2 * R ^ 2 / (M : ℝ) := by
    have hEΦ : drawsExpect p M Φ =
        ∑ k, integralOf p (fun d => h k d ^ 2) / ((M : ℝ) * (N : ℝ)) := by
      have : Φ = fun x => ∑ k, (∑ i, h k (x i)) ^ 2 / ((M : ℝ) ^ 2 * (N : ℝ)) :=
        funext hrewrite
      rw [this, drawsExpect_sum]
      refine Finset.sum_congr rfl fun k _ => ?_
      have hk := drawsExpect_sum_sq p (h k) hp1 (hh0 k) M
      have hdiv : drawsExpect p M (fun x => (∑ i, h k (x i)) ^ 2 / ((M : ℝ) ^ 2 * (N : ℝ))) =
          drawsExpect p M (fun x => (∑ i, h k (x i)) ^ 2) / ((M : ℝ) ^ 2 * (N : ℝ)) := by
        simp only [drawsExpect, Finset.sum_div, mul_div_assoc]
      rw [hdiv, hk]
      by_cases hN0 : (N : ℝ) = 0
      · simp [hN0]
      · field_simp
    rw [hEΦ]
    have hstep : ∀ k, integralOf p (fun d => h k d ^ 2) / ((M : ℝ) * (N : ℝ)) ≤
        (∑ d, p d * c d k ^ 2) / ((M : ℝ) * (N : ℝ)) := fun k =>
      div_le_div_of_nonneg_right (hh2 k) (by positivity)
    refine (Finset.sum_le_sum fun k _ => hstep k).trans ?_
    rw [← Finset.sum_div, Finset.sum_comm]
    have hatom : ∀ d, ∑ k, p d * c d k ^ 2 = p d * V ^ 2 * (∑ k, b d k ^ 2) := by
      intro d
      rw [Finset.mul_sum]
      refine Finset.sum_congr rfl fun k _ => ?_
      simp only [hc, mul_pow, hs2]
      ring
    simp only [hatom]
    have hbound : ∀ d, p d * V ^ 2 * (∑ k, b d k ^ 2) / ((M : ℝ) * (N : ℝ)) ≤
        p d * (V ^ 2 * R ^ 2 / (M : ℝ)) := by
      intro d
      have hbd := hb d
      have hcalc : p d * V ^ 2 * (∑ k, b d k ^ 2) / ((M : ℝ) * (N : ℝ)) =
          p d * (V ^ 2 * designNormSq (b d) / (M : ℝ)) := by
        simp only [designNormSq]
        by_cases hN0 : (N : ℝ) = 0
        · simp [hN0]
        · field_simp
      rw [hcalc]
      refine mul_le_mul_of_nonneg_left ?_ (hp0 d)
      refine div_le_div_of_nonneg_right ?_ hMr.le
      exact mul_le_mul_of_nonneg_left hbd (sq_nonneg V)
    rw [Finset.sum_div]
    refine (Finset.sum_le_sum fun d _ => hbound d).trans ?_
    rw [← Finset.sum_mul, hp1, one_mul]
  obtain ⟨x, hx⟩ := exists_le_drawsExpect p hp0 hp1 M Φ
  refine ⟨Finset.univ.image x, ?_, fun k => (∑ i, c (x i) k) / (M : ℝ), ?_, hx.trans hexpect⟩
  · exact Finset.card_image_le.trans (by simp)
  · have heq : (fun k => (∑ i, c (x i) k) / (M : ℝ)) = ∑ i, (s (x i) / (M : ℝ)) • b (x i) := by
      funext k
      simp only [Finset.sum_apply, Pi.smul_apply, smul_eq_mul, hc, Finset.sum_div]
      exact Finset.sum_congr rfl fun i _ => by ring
    rw [heq]
    refine Submodule.sum_mem _ fun i _ => Submodule.smul_mem _ _ (Submodule.subset_span ?_)
    exact ⟨x i, by simp, rfl⟩

end MTerm

/-! ### SC-5's Jeffreys profile -/

/-- The profile of the Jeffreys marginal `(1+g)^{(ν−M)/2} (1+g(1−R²))^{−ν/2}`, with `ν = N − 1`
(`BayesianRegression.md`, D03), up to the Occam term. -/
noncomputable def jeffreysProfile (g ν R2 : ℝ) : ℝ :=
  (ν / 2) * (Real.log (1 + g) - Real.log (1 + g * (1 - R2)))

/-- **SC-5.** The Jeffreys profile increases in `R²`, for `g > 0` and `R² ≤ 1`, so
`collapsedLogMarginal_superset_ge` applies to it with Occam cost `½ log(1 + g)`. -/
theorem jeffreysProfile_mono {g ν R2 R2' : ℝ} (hg : 0 < g) (hν : 0 ≤ ν) (hR : R2 ≤ R2')
    (hR' : R2' ≤ 1) : jeffreysProfile g ν R2 ≤ jeffreysProfile g ν R2' := by
  unfold jeffreysProfile
  have hpos : 0 < 1 + g * (1 - R2') := by
    have : 0 ≤ g * (1 - R2') := mul_nonneg hg.le (by linarith)
    linarith
  have hle : 1 + g * (1 - R2') ≤ 1 + g * (1 - R2) := by nlinarith
  have := Real.log_le_log hpos hle
  exact mul_le_mul_of_nonneg_left (by linarith) (by positivity)

/-- The profile of the collapsed marginal under the proper σ² prior `IG(c, c·s²)`, BASS's default
(`Sigma2Prior.ScaledToResponse`), up to the Occam term and a constant common to every structure:
`−(ν/2 + c) log(κ + (1 + g(1 − R²))/(1 + g))`, with `κ = 2c·s²/TSS` the prior's sum on the total's
scale. The prior adds `c` to the exponent and `c·s²` to the sum (`invGammaProfile_zero`: at
`c = κ = 0` it is the Jeffreys profile). -/
noncomputable def invGammaProfile (g ν c κ R2 : ℝ) : ℝ :=
  -(ν / 2 + c) * Real.log (κ + (1 + g * (1 - R2)) / (1 + g))

theorem invGammaProfile_zero {g R2 : ℝ} (hg : 0 < g) (hR : R2 ≤ 1) (ν : ℝ) :
    invGammaProfile g ν 0 0 R2 = jeffreysProfile g ν R2 := by
  unfold invGammaProfile jeffreysProfile
  have h₁ : 1 + g * (1 - R2) ≠ 0 := by
    have : 0 ≤ g * (1 - R2) := mul_nonneg hg.le (by linarith)
    linarith
  have h₂ : 1 + g ≠ 0 := by linarith
  rw [zero_add, add_zero, Real.log_div h₁ h₂]
  ring

/-- **SC-5 under IG(c, c·s²).** The profile increases in `R²`, for `g > 0`, `κ ≥ 0`, `ν/2 + c ≥ 0`
and `R² ≤ 1`, so `collapsedLogMarginal_superset_ge` applies with the same Occam cost
`½ log(1 + g)`. -/
theorem invGammaProfile_mono {g ν c κ R2 R2' : ℝ} (hg : 0 < g) (hνc : 0 ≤ ν / 2 + c) (hκ : 0 ≤ κ)
    (hR : R2 ≤ R2') (hR' : R2' ≤ 1) : invGammaProfile g ν c κ R2 ≤ invGammaProfile g ν c κ R2' := by
  unfold invGammaProfile
  have hg1 : 0 < 1 + g := by linarith
  have hpos : 0 < κ + (1 + g * (1 - R2')) / (1 + g) := by
    have : 0 ≤ g * (1 - R2') := mul_nonneg hg.le (by linarith)
    have : 0 < (1 + g * (1 - R2')) / (1 + g) := div_pos (by linarith) hg1
    linarith
  have hle : κ + (1 + g * (1 - R2')) / (1 + g) ≤ κ + (1 + g * (1 - R2)) / (1 + g) := by
    have : 1 + g * (1 - R2') ≤ 1 + g * (1 - R2) := by nlinarith
    have := div_le_div_of_nonneg_right this hg1.le
    linarith
  have hlog := Real.log_le_log hpos hle
  have := mul_le_mul_of_nonneg_left hlog hνc
  linarith

section Nesting

variable {N : ℕ}

/-- The residual sum of squares against a set of fits, per case: the smallest `‖y − v‖ₙ²`. -/
noncomputable def rssOf (V : Set (Fin N → ℝ)) (y : Fin N → ℝ) : ℝ :=
  sInf ((fun v => designNormSq (y - v)) '' V)

theorem rssOf_nonneg (V : Set (Fin N → ℝ)) (y : Fin N → ℝ) : 0 ≤ rssOf V y :=
  Real.sInf_nonneg (by
    rintro _ ⟨v, -, rfl⟩
    exact designNormSq_nonneg _)

/-- A larger set of fits has a smaller residual. -/
theorem rssOf_anti {V V' : Set (Fin N → ℝ)} (hV : V.Nonempty) (hsub : V ⊆ V') (y : Fin N → ℝ) :
    rssOf V' y ≤ rssOf V y := by
  refine csInf_le_csInf ⟨0, ?_⟩ (hV.image _) (image_mono hsub)
  rintro _ ⟨v, -, rfl⟩
  exact designNormSq_nonneg _

/-- `R² = 1 − RSS/TSS`, against a set of fits, for a total `tss`. -/
noncomputable def rSquaredOf (V : Set (Fin N → ℝ)) (y : Fin N → ℝ) (tss : ℝ) : ℝ :=
  1 - rssOf V y / tss

theorem rSquaredOf_le_one (V : Set (Fin N → ℝ)) (y : Fin N → ℝ) {tss : ℝ} (ht : 0 < tss) :
    rSquaredOf V y tss ≤ 1 := by
  have := div_nonneg (rssOf_nonneg V y) ht.le
  unfold rSquaredOf
  linarith

/-- `R²` cannot fall when the set of fits grows. -/
theorem rSquaredOf_mono {V V' : Set (Fin N → ℝ)} (hV : V.Nonempty) (hsub : V ⊆ V')
    (y : Fin N → ℝ) {tss : ℝ} (ht : 0 < tss) : rSquaredOf V y tss ≤ rSquaredOf V' y tss := by
  have := div_le_div_of_nonneg_right (rssOf_anti hV hsub y) ht.le
  unfold rSquaredOf
  linarith

/-- **SC-5 for D03's Jeffreys marginal.** With `V ⊆ V'` the column spaces of a structure and of its
superset by `k` columns, the collapsed log marginal falls by at most `(k/2) log(1 + g)`, whatever
`N`
and whatever variance the new columns explain. -/
theorem jeffreys_superset_ge {V V' : Set (Fin N → ℝ)} (hV : V.Nonempty) (hsub : V ⊆ V')
    (y : Fin N → ℝ) {tss g ν : ℝ} (ht : 0 < tss) (hg : 0 < g) (hν : 0 ≤ ν) (M k : ℕ) :
    collapsedLogMarginal (Real.log (1 + g) / 2) (jeffreysProfile g ν) M (rSquaredOf V y tss) -
        (k : ℝ) * (Real.log (1 + g) / 2) ≤
      collapsedLogMarginal (Real.log (1 + g) / 2) (jeffreysProfile g ν) (M + k)
        (rSquaredOf V' y tss) :=
  collapsedLogMarginal_superset_ge M k
    (jeffreysProfile_mono hg hν (rSquaredOf_mono hV hsub y ht) (rSquaredOf_le_one V' y ht))

/-- **SC-5 under BASS's default σ² prior, IG(c, c·s²).** With `V ⊆ V'` the column spaces of a
structure and of its superset by `k` columns, the collapsed log marginal falls by at most
`(k/2) log(1 + g)`, as under the Jeffreys form, whatever `N` and whatever variance the new columns
explain. `s²` reads the response alone, so `κ` is the same for both structures. -/
theorem invGamma_superset_ge {V V' : Set (Fin N → ℝ)} (hV : V.Nonempty) (hsub : V ⊆ V')
    (y : Fin N → ℝ) {tss g ν c κ : ℝ} (ht : 0 < tss) (hg : 0 < g) (hνc : 0 ≤ ν / 2 + c)
    (hκ : 0 ≤ κ) (M k : ℕ) :
    collapsedLogMarginal (Real.log (1 + g) / 2) (invGammaProfile g ν c κ) M (rSquaredOf V y tss) -
        (k : ℝ) * (Real.log (1 + g) / 2) ≤
      collapsedLogMarginal (Real.log (1 + g) / 2) (invGammaProfile g ν c κ) (M + k)
        (rSquaredOf V' y tss) :=
  collapsedLogMarginal_superset_ge M k
    (invGammaProfile_mono hg hνc hκ (rSquaredOf_mono hV hsub y ht) (rSquaredOf_le_one V' y ht))

end Nesting

end BASS
