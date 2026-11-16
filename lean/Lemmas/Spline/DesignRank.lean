import Mathlib.Data.Finset.Sort
import Mathlib.Data.Real.Basic
import Mathlib.LinearAlgebra.Basis.VectorSpace
import Mathlib.LinearAlgebra.Dimension.Constructions
import Mathlib.LinearAlgebra.LinearIndependent.Lemmas
import Mathlib.LinearAlgebra.Matrix.NonsingularInverse
import Mathlib.LinearAlgebra.Pi
import Mathlib.Order.Interval.Finset.Fin
import Mathlib.Tactic.Linarith

/-!
# Ordered supports and full column rank (S7)

After Schumaker's Theorem 4.61. A design whose columns have interval supports with nondecreasing
ends, sampled at sorted sites, has full column rank exactly when an increasing selection of sites
puts the `k`-th selected site in the `k`-th support. A greedy pass over the sites decides that.

* `det_eq_zero_of_zero_block` (S7a): an `n × n` matrix with a zero block of `r` rows and `s`
  columns, `r + s > n`, is singular. Every permutation meets the block.
* `inSupport_of_det_ne_zero` (S7b): a nonsingular square collocation matrix has its `i`-th site in
  the `i`-th support. A site left of its support makes rows `≤ i` vanish on columns `≥ i`; a site
  right of it makes rows `≥ i` vanish on columns `≤ i`. The supports may close on the left at a
  full-multiplicity knot, (4.133).
* `exists_selection_of_fullRank` (S7c): full column rank gives `n` rows whose square submatrix is
  nonsingular (`exists_strictMono_det_ne_zero`). Sorted, they are the selection.
* `exists_isSelection_iff_greedy` (S7d): the first-fit pass succeeds exactly when a selection
  exists. Its picks lie at or before any selection's, for any membership predicate.
* `fullRank_of_selection` (S7e): with Theorem 4.61's sufficiency half as the named premise
  `h461`, a selection gives full column rank. `fullRank_iff_greedy` closes the loop.
  `det_indicator_eq_one` proves the premise for order 1, so it is satisfiable.
-/

open Matrix

namespace Spline

/-! ### S7a, the zero block -/

/-- **S7a.** A square matrix with every entry zero on `rows × cols`, where
`|rows| + |cols| > n`, has determinant zero. -/
theorem det_eq_zero_of_zero_block {ι R : Type*} [Fintype ι] [DecidableEq ι] [CommRing R]
    (M : Matrix ι ι R) (rows cols : Finset ι) (hcard : Fintype.card ι < rows.card + cols.card)
    (hz : ∀ i ∈ rows, ∀ j ∈ cols, M i j = 0) : M.det = 0 := by
  rw [det_apply]
  refine Finset.sum_eq_zero fun σ _ => ?_
  obtain ⟨j, hj, hσ⟩ : ∃ j ∈ cols, σ j ∈ rows := by
    by_contra h
    push Not at h
    have hmaps : Set.MapsTo σ.symm (rows : Set ι) ((colsᶜ : Finset ι) : Set ι) := by
      intro i hi
      rw [Finset.mem_coe] at hi
      rw [Finset.mem_coe, Finset.mem_compl]
      intro hc
      exact h _ hc (by simpa using hi)
    have h1 := Finset.card_le_card_of_injOn σ.symm hmaps σ.symm.injective.injOn
    have h2 := cols.card_le_univ
    rw [Finset.card_compl] at h1
    omega
  have hprod : ∏ i, M (σ i) i = 0 := Finset.prod_eq_zero (Finset.mem_univ j) (hz _ hσ _ hj)
  rw [hprod, smul_zero]

theorem det_eq_zero_of_left {R : Type*} [CommRing R] {n : ℕ} (M : Matrix (Fin n) (Fin n) R)
    (i : Fin n) (h : ∀ ν μ, ν ≤ i → i ≤ μ → M ν μ = 0) : M.det = 0 :=
  det_eq_zero_of_zero_block M (Finset.Iic i) (Finset.Ici i)
    (by have := i.isLt; rw [Fin.card_Iic, Fin.card_Ici, Fintype.card_fin]; omega)
    fun ν hν μ hμ => h ν μ (Finset.mem_Iic.1 hν) (Finset.mem_Ici.1 hμ)

theorem det_eq_zero_of_right {R : Type*} [CommRing R] {n : ℕ} (M : Matrix (Fin n) (Fin n) R)
    (i : Fin n) (h : ∀ ν μ, i ≤ ν → μ ≤ i → M ν μ = 0) : M.det = 0 :=
  det_eq_zero_of_zero_block M (Finset.Ici i) (Finset.Iic i)
    (by have := i.isLt; rw [Fin.card_Iic, Fin.card_Ici, Fintype.card_fin]; omega)
    fun ν hν μ hμ => h ν μ (Finset.mem_Ici.1 hν) (Finset.mem_Iic.1 hμ)

/-! ### S7b, necessity -/

/-- Membership in the `j`-th support: right of `ℓ j` (or at it, when `closed j`, the
full-multiplicity case of (4.133)) and left of `u j`. -/
def InSupport {n : ℕ} (closed : Fin n → Prop) (ℓ u : Fin n → ℝ) (j : Fin n) (x : ℝ) : Prop :=
  (ℓ j < x ∨ (closed j ∧ ℓ j = x)) ∧ x < u j

/-- **S7b.** Sorted sites `t`, functions `φ j` vanishing off the `j`-th support, ends `ℓ` and `u`
nondecreasing, and a support that closes at `ℓ μ` only if every earlier support with the same
left end closes too (B-splines satisfy this). If the collocation matrix `(φ j (t i))` is
nonsingular, then `t i` lies in the `i`-th support for every `i`. -/
theorem inSupport_of_det_ne_zero {n : ℕ} (t : Fin n → ℝ) (ht : Monotone t)
    (φ : Fin n → ℝ → ℝ) (closed : Fin n → Prop) (ℓ u : Fin n → ℝ) (hℓ : Monotone ℓ)
    (hu : Monotone u) (hclosed : ∀ i μ, i ≤ μ → ℓ i = ℓ μ → closed μ → closed i)
    (hφ : ∀ j x, φ j x ≠ 0 → InSupport closed ℓ u j x)
    (hdet : (Matrix.of fun i j => φ j (t i)).det ≠ 0) (i : Fin n) :
    InSupport closed ℓ u i (t i) := by
  by_contra hnot
  unfold InSupport at hnot
  rcases not_and_or.1 hnot with hl | hr
  · push Not at hl
    apply hdet
    refine det_eq_zero_of_left _ i fun ν μ hν hμ => ?_
    simp only [of_apply]
    by_contra hne
    obtain ⟨h1, -⟩ := hφ μ (t ν) hne
    have htν : t ν ≤ t i := ht hν
    have hℓμ : ℓ i ≤ ℓ μ := hℓ hμ
    rcases h1 with h1 | ⟨hc, he⟩
    · linarith [hl.1]
    · have e1 : ℓ i = t i := by linarith [hl.1]
      have e2 : ℓ i = ℓ μ := by linarith [hl.1]
      exact hl.2 (hclosed i μ hμ e2 hc) e1
  · push Not at hr
    apply hdet
    refine det_eq_zero_of_right _ i fun ν μ hν hμ => ?_
    simp only [of_apply]
    by_contra hne
    obtain ⟨-, h2⟩ := hφ μ (t ν) hne
    have h3 := ht hν
    have h4 := hu hμ
    linarith

/-- **S7b, open supports.** With every support open, `ℓ i < t i < u i`. -/
theorem mem_Ioo_of_det_ne_zero {n : ℕ} (t : Fin n → ℝ) (ht : Monotone t) (φ : Fin n → ℝ → ℝ)
    (ℓ u : Fin n → ℝ) (hℓ : Monotone ℓ) (hu : Monotone u)
    (hφ : ∀ j x, φ j x ≠ 0 → x ∈ Set.Ioo (ℓ j) (u j))
    (hdet : (Matrix.of fun i j => φ j (t i)).det ≠ 0) (i : Fin n) : t i ∈ Set.Ioo (ℓ i) (u i) := by
  obtain ⟨h1 | ⟨⟨⟩, -⟩, h2⟩ := inSupport_of_det_ne_zero t ht φ (fun _ => False) ℓ u hℓ hu
    (fun _ _ _ _ h => h) (fun j x hx => ⟨Or.inl (hφ j x hx).1, (hφ j x hx).2⟩) hdet i
  exact ⟨h1, h2⟩

/-! ### S7c, designs -/

/-- Full column rank gives `n` rows, in increasing order, whose square submatrix is nonsingular. -/
theorem exists_strictMono_det_ne_zero {N n : ℕ} (M : Matrix (Fin N) (Fin n) ℝ)
    (hrank : ∀ c : Fin n → ℝ, M *ᵥ c = 0 → c = 0) :
    ∃ e : Fin n → Fin N, StrictMono e ∧ (M.submatrix e id).det ≠ 0 := by
  -- the rows span: a functional vanishing on every row is `⟨·, c⟩` with `M c = 0`
  have hspan : Submodule.span ℝ (Set.range M.row) = ⊤ := by
    by_contra hne
    obtain ⟨f, hf0, hle⟩ := Submodule.exists_le_ker_of_lt_top _ (lt_top_iff_ne_top.2 hne)
    set c : Fin n → ℝ := fun j => f fun i => if j = i then 1 else 0 with hc_def
    have hfc : ∀ v, f v = v ⬝ᵥ c := fun v => by
      rw [LinearMap.pi_apply_eq_sum_univ f v]
      simp [dotProduct, hc_def, smul_eq_mul]
    have hc : M *ᵥ c = 0 := by
      funext i
      have hi : M.row i ∈ LinearMap.ker f := hle (Submodule.subset_span ⟨i, rfl⟩)
      rw [LinearMap.mem_ker, hfc] at hi
      simpa [mulVec, dotProduct, Matrix.row] using hi
    apply hf0
    refine LinearMap.ext fun v => ?_
    rw [hfc, hrank c hc, dotProduct_zero, LinearMap.zero_apply]
  -- a linearly independent subfamily of rows spans, so it has `n` members
  obtain ⟨κ, a, ha, hspa, hli⟩ := exists_linearIndependent' (K := ℝ) (v := M.row)
  have : Finite κ := Finite.of_injective a ha
  have : Fintype κ := Fintype.ofFinite κ
  have hsp : ⊤ ≤ Submodule.span ℝ (Set.range (M.row ∘ a)) := by rw [hspa, hspan]
  have hcard : Fintype.card κ = n := by
    have := Module.finrank_eq_card_basis (Module.Basis.mk hli hsp)
    rw [Module.finrank_fin_fun] at this
    exact this.symm
  -- sort the chosen rows
  set s : Finset (Fin N) := Finset.univ.image a with hs_def
  have hs : s.card = n := by
    rw [hs_def, Finset.card_image_of_injective _ ha, Finset.card_univ, hcard]
  set e := s.orderEmbOfFin hs
  have hmem : ∀ k, ∃ x : κ, a x = e k := fun k => by
    have h1 : e k ∈ Finset.univ.image a := s.orderEmbOfFin_mem hs k
    obtain ⟨x, -, hx⟩ := Finset.mem_image.1 h1
    exact ⟨x, hx⟩
  choose g hg using hmem
  have hginj : Function.Injective g := fun k₁ k₂ h => e.injective (by rw [← hg, ← hg, h])
  refine ⟨e, e.strictMono, ?_⟩
  have hli' : LinearIndependent ℝ (M.submatrix e id).row := by
    have := hli.comp g hginj
    convert this using 1
    · funext k
      simp only [Function.comp, hg]
      rfl
    all_goals rfl
  rw [← isUnit_iff_ne_zero, ← isUnit_iff_isUnit_det]
  exact linearIndependent_rows_iff_isUnit.1 hli'

/-- **S7c.** If the `N × n` design `(φ j (t i))` on sorted sites has full column rank, some
increasing selection `e` puts site `e k` in the `k`-th support for every `k`. -/
theorem exists_selection_of_fullRank {N n : ℕ} (t : Fin N → ℝ) (ht : Monotone t)
    (φ : Fin n → ℝ → ℝ) (closed : Fin n → Prop) (ℓ u : Fin n → ℝ) (hℓ : Monotone ℓ)
    (hu : Monotone u) (hclosed : ∀ i μ, i ≤ μ → ℓ i = ℓ μ → closed μ → closed i)
    (hφ : ∀ j x, φ j x ≠ 0 → InSupport closed ℓ u j x)
    (hrank : ∀ c : Fin n → ℝ, (Matrix.of fun i j => φ j (t i)) *ᵥ c = 0 → c = 0) :
    ∃ e : Fin n → Fin N, StrictMono e ∧ ∀ k, InSupport closed ℓ u k (t (e k)) := by
  obtain ⟨e, he, hdet⟩ := exists_strictMono_det_ne_zero _ hrank
  exact ⟨e, he, inSupport_of_det_ne_zero (t ∘ e) (ht.comp he.monotone) φ closed ℓ u hℓ hu hclosed
    hφ hdet⟩

/-! ### S7d, the greedy check -/

section Greedy

variable (P : ℕ → ℕ → Prop) [DecidableRel P] (N : ℕ)

/-- An increasing selection of `n` sites below `N`, site `e k` serving column `k`. -/
def IsSelection (n : ℕ) (e : ℕ → ℕ) : Prop :=
  (∀ k < n, e k < N ∧ P k (e k)) ∧ ∀ k, k + 1 < n → e k < e (k + 1)

/-- The first site at or after `lo` that serves column `k`. -/
noncomputable def nextPick (lo k : ℕ) : Option ℕ :=
  if h : ∃ i, i < N ∧ lo ≤ i ∧ P k i then some (Nat.find h) else none

/-- The greedy check: column `0` takes the first site that serves it, and column `k + 1` the first
site after column `k`'s pick. -/
noncomputable def greedy : ℕ → Option ℕ
  | 0 => nextPick P N 0 0
  | k + 1 => (greedy k).bind fun p => nextPick P N (p + 1) (k + 1)

variable {P N}

theorem nextPick_spec {lo k p : ℕ} (h : nextPick P N lo k = some p) :
    lo ≤ p ∧ p < N ∧ P k p := by
  unfold nextPick at h
  split_ifs at h with hex
  cases h
  obtain ⟨h1, h2, h3⟩ := Nat.find_spec hex
  exact ⟨h2, h1, h3⟩

theorem nextPick_le {lo k i : ℕ} (hi : lo ≤ i ∧ i < N ∧ P k i) :
    ∃ p, nextPick P N lo k = some p ∧ p ≤ i := by
  have hex : ∃ i, i < N ∧ lo ≤ i ∧ P k i := ⟨i, hi.2.1, hi.1, hi.2.2⟩
  exact ⟨Nat.find hex, by unfold nextPick; rw [dif_pos hex],
    Nat.find_min' hex ⟨hi.2.1, hi.1, hi.2.2⟩⟩

theorem greedy_spec {k p : ℕ} (h : greedy P N k = some p) : p < N ∧ P k p := by
  cases k with
  | zero => exact (nextPick_spec h).2
  | succ k =>
    obtain ⟨p', -, h'⟩ := Option.bind_eq_some_iff.1 h
    exact (nextPick_spec h').2

theorem greedy_lt_succ {k p q : ℕ} (hp : greedy P N k = some p)
    (hq : greedy P N (k + 1) = some q) : p < q := by
  simp only [greedy, hp, Option.bind_some] at hq
  exact (nextPick_spec hq).1

/-- The greedy picks lie at or before any selection's. -/
theorem greedy_le_of_isSelection {n : ℕ} {e : ℕ → ℕ} (he : IsSelection P N n e) :
    ∀ k < n, ∃ p, greedy P N k = some p ∧ p ≤ e k := by
  intro k
  induction k with
  | zero => exact fun hk => nextPick_le ⟨Nat.zero_le _, (he.1 0 hk).1, (he.1 0 hk).2⟩
  | succ k ih =>
    intro hk
    obtain ⟨p, hp, hpe⟩ := ih (by omega)
    have hlt := he.2 k hk
    obtain ⟨q, hq, hqe⟩ :=
      nextPick_le (P := P) (N := N) (lo := p + 1)
        ⟨by omega, (he.1 (k + 1) hk).1, (he.1 (k + 1) hk).2⟩
    exact ⟨q, by simp only [greedy, hp, Option.bind_some, hq], hqe⟩

theorem isSelection_greedy {n : ℕ} (h : ∀ k < n, (greedy P N k).isSome) :
    IsSelection P N n fun k => (greedy P N k).getD 0 := by
  have hsome : ∀ k < n, ∃ p, greedy P N k = some p := fun k hk =>
    Option.isSome_iff_exists.1 (h k hk)
  refine ⟨fun k hk => ?_, fun k hk => ?_⟩
  · obtain ⟨p, hp⟩ := hsome k hk
    simp only [hp, Option.getD_some]
    exact greedy_spec hp
  · obtain ⟨p, hp⟩ := hsome k (by omega)
    obtain ⟨q, hq⟩ := hsome (k + 1) hk
    simp only [hp, hq, Option.getD_some]
    exact greedy_lt_succ hp hq

/-- **S7d.** The greedy check succeeds exactly when a selection exists, whatever the predicate. -/
theorem exists_isSelection_iff_greedy {n : ℕ} :
    (∃ e, IsSelection P N n e) ↔ ∀ k < n, (greedy P N k).isSome := by
  constructor
  · rintro ⟨e, he⟩ k hk
    obtain ⟨p, hp, -⟩ := greedy_le_of_isSelection he k hk
    simp [hp]
  · exact fun h => ⟨_, isSelection_greedy h⟩

end Greedy

/-- A predicate on `Fin`-indexed columns and sites, read on `ℕ`. -/
def liftPred {n N : ℕ} (Q : Fin n → Fin N → Prop) (k i : ℕ) : Prop :=
  ∃ (hk : k < n) (hi : i < N), Q ⟨k, hk⟩ ⟨i, hi⟩

theorem isSelection_iff_strictMono {n N : ℕ} (Q : Fin n → Fin N → Prop) :
    (∃ e, IsSelection (liftPred Q) N n e) ↔ ∃ e : Fin n → Fin N, StrictMono e ∧ ∀ k, Q k (e k) := by
  constructor
  · rintro ⟨e, he⟩
    have hstep : ∀ d k, k + d + 1 < n → e k < e (k + d + 1) := by
      intro d
      induction d with
      | zero => exact fun k hk => he.2 k hk
      | succ d ih =>
        intro k hk
        exact (ih k (by omega)).trans (he.2 (k + d + 1) (by omega))
    refine ⟨fun k => ⟨e k, (he.1 k k.2).1⟩, fun a b hab => ?_, fun k => ?_⟩
    · have hab' : (a : ℕ) < b := hab
      obtain ⟨d, hd⟩ : ∃ d, (b : ℕ) = a + d + 1 := ⟨b - a - 1, by omega⟩
      change e a < e b
      rw [hd]
      exact hstep d a (by omega)
    · obtain ⟨hk, hi, hQ⟩ := (he.1 k k.2).2
      exact hQ
  · rintro ⟨e, he, hQ⟩
    refine ⟨fun k => if hk : k < n then e ⟨k, hk⟩ else 0, fun k hk => ?_, fun k hk => ?_⟩
    · simp only [dif_pos hk]
      exact ⟨(e ⟨k, hk⟩).isLt, hk, (e ⟨k, hk⟩).isLt, hQ _⟩
    · simp only [dif_pos (show k < n by omega), dif_pos hk]
      exact he (show (⟨k, by omega⟩ : Fin n) < ⟨k + 1, hk⟩ from Fin.mk_lt_mk.2 (by omega))

/-! ### S7e, sufficiency with Theorem 4.61 as premise -/

/-- **S7e.** Theorem 4.61's sufficiency half as the premise `h461`: sites strictly increasing and
each in its support give a nonsingular collocation matrix. Then a selection on strictly increasing
sites gives the design full column rank. -/
theorem fullRank_of_selection {N n : ℕ} (t : Fin N → ℝ) (ht : StrictMono t) (φ : Fin n → ℝ → ℝ)
    (closed : Fin n → Prop) (ℓ u : Fin n → ℝ)
    (h461 : ∀ s : Fin n → ℝ, StrictMono s → (∀ k, InSupport closed ℓ u k (s k)) →
      (Matrix.of fun i j => φ j (s i)).det ≠ 0)
    {e : Fin n → Fin N} (he : StrictMono e) (hsel : ∀ k, InSupport closed ℓ u k (t (e k))) :
    ∀ c : Fin n → ℝ, (Matrix.of fun i j => φ j (t i)) *ᵥ c = 0 → c = 0 := by
  intro c hc
  have hdet := h461 (t ∘ e) (ht.comp he) hsel
  have hunit : IsUnit (Matrix.of fun i j => φ j ((t ∘ e) i)) :=
    (isUnit_iff_isUnit_det _).2 (isUnit_iff_ne_zero.2 hdet)
  apply mulVec_injective_iff_isUnit.2 hunit
  rw [mulVec_zero]
  funext k
  have := congrFun hc (e k)
  simpa [mulVec, dotProduct] using this

open Classical in
/-- **S7, the loop closed.** Given Theorem 4.61's sufficiency half, a design on strictly increasing
sites has full column rank exactly when the greedy check succeeds. -/
theorem fullRank_iff_greedy {N n : ℕ} (t : Fin N → ℝ) (ht : StrictMono t) (φ : Fin n → ℝ → ℝ)
    (closed : Fin n → Prop) (ℓ u : Fin n → ℝ) (hℓ : Monotone ℓ) (hu : Monotone u)
    (hclosed : ∀ i μ, i ≤ μ → ℓ i = ℓ μ → closed μ → closed i)
    (hφ : ∀ j x, φ j x ≠ 0 → InSupport closed ℓ u j x)
    (h461 : ∀ s : Fin n → ℝ, StrictMono s → (∀ k, InSupport closed ℓ u k (s k)) →
      (Matrix.of fun i j => φ j (s i)).det ≠ 0) :
    (∀ c : Fin n → ℝ, (Matrix.of fun i j => φ j (t i)) *ᵥ c = 0 → c = 0) ↔
      ∀ k < n, (greedy (liftPred fun k i => InSupport closed ℓ u k (t i)) N k).isSome := by
  rw [← exists_isSelection_iff_greedy, isSelection_iff_strictMono]
  constructor
  · exact exists_selection_of_fullRank t ht.monotone φ closed ℓ u hℓ hu hclosed hφ
  · rintro ⟨e, he, hsel⟩
    exact fullRank_of_selection t ht φ closed ℓ u h461 he hsel

/-! ### The premise at order 1, and witnesses -/

/-- The order-1 B-splines on knots `y`: the indicators of the open intervals `(y j, y (j + 1))`. -/
noncomputable def indicatorBasis {n : ℕ} (y : ℕ → ℝ) (j : Fin n) (x : ℝ) : ℝ :=
  if y j < x ∧ x < y (j + 1) then 1 else 0

/-- Theorem 4.61's sufficiency half at order 1: at sites in their own supports, the collocation
matrix of the indicators is the identity. -/
theorem det_indicator_eq_one {n : ℕ} (y : ℕ → ℝ) (hy : StrictMono y) (s : Fin n → ℝ)
    (hs : ∀ k : Fin n, y k < s k ∧ s k < y (k + 1)) :
    (Matrix.of fun i j : Fin n => indicatorBasis y j (s i)).det = 1 := by
  have h1 : (Matrix.of fun i j : Fin n => indicatorBasis y j (s i)) = 1 := by
    ext i j
    simp only [of_apply, one_apply, indicatorBasis]
    by_cases hij : i = j
    · subst hij
      simp [hs i]
    · rw [if_neg hij, if_neg]
      rintro ⟨h1, h2⟩
      rcases lt_or_gt_of_ne (fun h => hij (Fin.ext h)) with hlt | hlt
      · have : y (i + 1) ≤ y j := hy.monotone (by omega)
        linarith [(hs i).2]
      · have : y (j + 1) ≤ y i := hy.monotone (by omega)
        linarith [(hs i).1]
  rw [h1, det_one]

/-- **The premise is satisfiable:** Theorem 4.61's sufficiency half holds at order 1. -/
theorem h461_indicator {n : ℕ} (y : ℕ → ℝ) (hy : StrictMono y) :
    ∀ s : Fin n → ℝ, StrictMono s →
      (∀ k, InSupport (fun _ => False) (fun j : Fin n => y j) (fun j => y (j + 1)) k (s k)) →
      (Matrix.of fun i j => indicatorBasis y j (s i)).det ≠ 0 := by
  intro s _ hsupp
  rw [det_indicator_eq_one y hy s fun k => ?_]
  · exact one_ne_zero
  · obtain ⟨h1 | ⟨⟨⟩, -⟩, h2⟩ := hsupp k
    exact ⟨h1, h2⟩

/-- The loop's hypotheses hold for order-1 splines on knots `0, 1, 2` with sites `½, 3/2`. -/
example : ∀ c : Fin 2 → ℝ,
    (Matrix.of fun i j => indicatorBasis (fun k : ℕ => (k : ℝ)) j (![1 / 2, 3 / 2] i)) *ᵥ c = 0 →
      c = 0 :=
  fullRank_of_selection ![1 / 2, 3 / 2]
    (Fin.strictMono_iff_lt_succ.2 fun i => by fin_cases i; norm_num)
    (indicatorBasis fun k : ℕ => (k : ℝ)) (fun _ => False) (fun j => ((j : ℕ) : ℝ))
    (fun j => (((j : ℕ) + 1 : ℕ) : ℝ))
    (h461_indicator (fun k : ℕ => (k : ℝ)) Nat.strictMono_cast) strictMono_id
    fun k => by fin_cases k <;> exact ⟨Or.inl (by norm_num), by norm_num⟩

end Spline
