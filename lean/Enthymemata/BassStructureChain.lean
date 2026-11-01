import Mathlib.Algebra.BigOperators.Field
import Mathlib.Algebra.Order.Star.Real
import Mathlib.Data.Fintype.BigOperators
import Mathlib.Data.Real.Basic
import Mathlib.LinearAlgebra.Matrix.Hadamard
import Mathlib.LinearAlgebra.Matrix.PosDef
import Mathlib.Tactic.FieldSimp
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.Ring

/-!
# BASS — the structure chain (SC-1, SC-2, SC-3)

The self-contained statements of `Protolemmata/BASS/structure-chain-lemmas.md`.

* **SC-1.** `mhRatio_mul_swap`: a Metropolis–Hastings ratio and its reverse multiply to 1.
  `mh_detailed_balance`: the kernel `q · min 1 r` is reversible. `not_isMhPair_of_mul_ne_one`: so a
  printed pair of ratios whose product is not 1 comes from no `π` and `q`.
* **SC-2.** `gap_le_cut`, the easy half of Cheeger: a Poincaré constant `γ` is at most any set's
  flow over its variance, with `varianceOf_indicator` and `dirichlet_indicator` the two
  computations. `dirichlet_indicator_of_stationary` and `gap_le_cut_of_stationary` take a stationary
  kernel in place of a reversible one, since a sweep's block scan need not be reversible.
* **SC-3.** `khaos_precision_entry`, its matrix form `khaos_precision_eq`, the conjugate form
  `khaos_precision_eq_add_inv` (`G ∘ A = A + S₀⁻¹`), and `khaos_precision_posDef`.
  `khaos_precision_of_const` and `khaos_precision_of_diagonal` are the two cases where it reduces
  to a form with a known fast path: a scalar times the Gram, and a diagonal.
* **SC-5,** posed by the walk of Ronen et al. `gap_le_separator`: a set's flow runs through any
  separator, which is the form that walk's Lemma 3 takes when no single state gates the cut.
  `collapsedLogMarginal_superset_ge`: why BASS has no such gate, since adding columns costs only
  the Occam term.
-/

open Finset

namespace BASS

/-! ### SC-1 -/

/-- The Metropolis–Hastings ratio, reverse over forward. -/
noncomputable def mhRatio {S : Type*} (π : S → ℝ) (q : S → S → ℝ) (s t : S) : ℝ :=
  π t * q t s / (π s * q s t)

/-- **SC-1, reciprocity.** With nonzero `π` and `q` both ways, a ratio and its reverse multiply
to 1. -/
theorem mhRatio_mul_swap {S : Type*} (π : S → ℝ) (q : S → S → ℝ) (s t : S)
    (hs : π s ≠ 0) (ht : π t ≠ 0) (hst : q s t ≠ 0) (hts : q t s ≠ 0) :
    mhRatio π q s t * mhRatio π q t s = 1 := by
  unfold mhRatio
  field_simp

/-- **SC-1, the contrapositive that reads a manuscript.** A pair of acceptance factors whose
product is not 1 is the ratio pair of no `π` and `q` that are nonzero on both states. -/
theorem not_isMhPair_of_mul_ne_one {S : Type*} (s t : S) (a b : ℝ) (hab : a * b ≠ 1) :
    ¬ ∃ (π : S → ℝ) (q : S → S → ℝ), π s ≠ 0 ∧ π t ≠ 0 ∧ q s t ≠ 0 ∧ q t s ≠ 0 ∧
      mhRatio π q s t = a ∧ mhRatio π q t s = b := by
  rintro ⟨π, q, hs, ht, hst, hts, rfl, rfl⟩
  exact hab (mhRatio_mul_swap π q s t hs ht hst hts)

/-- **SC-1, detailed balance.** With positive `π` and `q` both ways, the kernel
`q · min 1 r` satisfies detailed balance with respect to `π`: both sides equal
`min (π s · q s t) (π t · q t s)`. -/
theorem mh_detailed_balance {S : Type*} (π : S → ℝ) (q : S → S → ℝ) (s t : S)
    (hs : 0 < π s) (ht : 0 < π t) (hst : 0 < q s t) (hts : 0 < q t s) :
    π s * (q s t * min 1 (mhRatio π q s t)) = π t * (q t s * min 1 (mhRatio π q t s)) := by
  have ha : 0 < π s * q s t := mul_pos hs hst
  have hb : 0 < π t * q t s := mul_pos ht hts
  have key : ∀ {a b : ℝ}, 0 < a → a * min 1 (b / a) = min a b := by
    intro a b h
    rw [mul_min_of_nonneg _ _ h.le, mul_one, mul_div_cancel₀ _ h.ne']
  calc π s * (q s t * min 1 (mhRatio π q s t))
      = (π s * q s t) * min 1 ((π t * q t s) / (π s * q s t)) := by rw [mhRatio]; ring
    _ = min (π s * q s t) (π t * q t s) := key ha
    _ = min (π t * q t s) (π s * q s t) := min_comm _ _
    _ = (π t * q t s) * min 1 ((π s * q s t) / (π t * q t s)) := (key hb).symm
    _ = π t * (q t s * min 1 (mhRatio π q t s)) := by rw [mhRatio]; ring

/-! ### SC-2 -/

section Flow

variable {S : Type*} [Fintype S] [DecidableEq S]

/-- The Dirichlet form of a kernel `P` under weights `π`. -/
noncomputable def dirichlet (π : S → ℝ) (P : S → S → ℝ) (f : S → ℝ) : ℝ :=
  (1 / 2) * ∑ s, ∑ t, π s * P s t * (f s - f t) ^ 2

/-- The variance of `f` under weights `π`. -/
noncomputable def varianceOf (π : S → ℝ) (f : S → ℝ) : ℝ :=
  ∑ s, π s * f s ^ 2 - (∑ s, π s * f s) ^ 2

/-- The flow out of `A`. -/
noncomputable def cutFlow (π : S → ℝ) (P : S → S → ℝ) (A : Finset S) : ℝ :=
  ∑ s ∈ A, ∑ t ∈ Aᶜ, π s * P s t

/-- The indicator of `A`. -/
def indicatorOf (A : Finset S) : S → ℝ := fun s => if s ∈ A then 1 else 0

/-- `∑_A π · g` as a weighted sum against the indicator of `A`. -/
private theorem sum_mul_indicator (π : S → ℝ) (A : Finset S) (k : ℕ) (hk : k ≠ 0) :
    ∑ s, π s * indicatorOf A s ^ k = ∑ s ∈ A, π s := by
  rw [← Finset.sum_add_sum_compl A]
  have h₁ : ∑ s ∈ A, π s * indicatorOf A s ^ k = ∑ s ∈ A, π s :=
    Finset.sum_congr rfl fun s hs => by simp [indicatorOf, hs]
  have h₂ : ∑ s ∈ Aᶜ, π s * indicatorOf A s ^ k = 0 :=
    Finset.sum_eq_zero fun s hs => by
      have : s ∉ A := Finset.mem_compl.1 hs
      simp [indicatorOf, this, hk]
  rw [h₁, h₂, add_zero]

/-- **SC-2, first computation.** The indicator's variance is `π(A)(1 − π(A))`. The form needs no
normalization of `π`; it is a variance when `π` is a probability vector. -/
theorem varianceOf_indicator (π : S → ℝ) (A : Finset S) :
    varianceOf π (indicatorOf A) = (∑ s ∈ A, π s) * (1 - ∑ s ∈ A, π s) := by
  have h₁ := sum_mul_indicator π A 2 two_ne_zero
  have h₂ := sum_mul_indicator π A 1 one_ne_zero
  simp only [pow_one] at h₂
  rw [varianceOf, h₁, h₂]
  ring

/-- The Dirichlet form of an indicator is half the flow across the cut, counted both ways. -/
theorem dirichlet_indicator_eq_half (π : S → ℝ) (P : S → S → ℝ) (A : Finset S) :
    dirichlet π P (indicatorOf A) = (1 / 2) * (cutFlow π P A + cutFlow π P Aᶜ) := by
  have hin : ∀ s ∈ A, ∑ t, π s * P s t * (indicatorOf A s - indicatorOf A t) ^ 2 =
      ∑ t ∈ Aᶜ, π s * P s t := by
    intro s hs
    rw [← Finset.sum_add_sum_compl A]
    have h₁ : ∑ t ∈ A, π s * P s t * (indicatorOf A s - indicatorOf A t) ^ 2 = 0 :=
      Finset.sum_eq_zero fun t ht => by simp [indicatorOf, hs, ht]
    have h₂ : ∑ t ∈ Aᶜ, π s * P s t * (indicatorOf A s - indicatorOf A t) ^ 2 =
        ∑ t ∈ Aᶜ, π s * P s t :=
      Finset.sum_congr rfl fun t ht => by
        have : t ∉ A := Finset.mem_compl.1 ht
        simp [indicatorOf, hs, this]
    rw [h₁, h₂, zero_add]
  have hout : ∀ s ∈ Aᶜ, ∑ t, π s * P s t * (indicatorOf A s - indicatorOf A t) ^ 2 =
      ∑ t ∈ A, π s * P s t := by
    intro s hs
    have hs' : s ∉ A := Finset.mem_compl.1 hs
    rw [← Finset.sum_add_sum_compl A]
    have h₁ : ∑ t ∈ A, π s * P s t * (indicatorOf A s - indicatorOf A t) ^ 2 =
        ∑ t ∈ A, π s * P s t :=
      Finset.sum_congr rfl fun t ht => by simp [indicatorOf, hs', ht]
    have h₂ : ∑ t ∈ Aᶜ, π s * P s t * (indicatorOf A s - indicatorOf A t) ^ 2 = 0 :=
      Finset.sum_eq_zero fun t ht => by
        have : t ∉ A := Finset.mem_compl.1 ht
        simp [indicatorOf, hs', this]
    rw [h₁, h₂, add_zero]
  rw [dirichlet, ← Finset.sum_add_sum_compl A, Finset.sum_congr rfl hin, Finset.sum_congr rfl hout,
    cutFlow, cutFlow, compl_compl]

/-- Reversibility makes the flow the same both ways across any cut. -/
theorem cutFlow_compl_of_reversible (π : S → ℝ) (P : S → S → ℝ)
    (hrev : ∀ s t, π s * P s t = π t * P t s) (A : Finset S) :
    cutFlow π P Aᶜ = cutFlow π P A := by
  rw [cutFlow, cutFlow, compl_compl, Finset.sum_comm]
  exact Finset.sum_congr rfl fun s _ => Finset.sum_congr rfl fun t _ => hrev t s

/-- **SC-2, second computation.** Under reversibility the Dirichlet form of `1_A` is `A`'s flow. -/
theorem dirichlet_indicator (π : S → ℝ) (P : S → S → ℝ)
    (hrev : ∀ s t, π s * P s t = π t * P t s) (A : Finset S) :
    dirichlet π P (indicatorOf A) = cutFlow π P A := by
  rw [dirichlet_indicator_eq_half, cutFlow_compl_of_reversible π P hrev]
  ring

/-- **SC-2.** For a kernel reversible with respect to `π`, a Poincaré constant `γ` satisfies
`γ π(A)(1 − π(A)) ≤ flow(A)` for every `A`. The protolemma's normalization `∑ π = 1` plays no part
in the inequality, and is dropped; it is what makes the left side a variance. -/
theorem gap_le_cut (π : S → ℝ) (P : S → S → ℝ)
    (hrev : ∀ s t, π s * P s t = π t * P t s) (γ : ℝ)
    (hgap : ∀ f : S → ℝ, γ * varianceOf π f ≤ dirichlet π P f) (A : Finset S) :
    γ * ((∑ s ∈ A, π s) * (1 - ∑ s ∈ A, π s)) ≤ cutFlow π P A := by
  have := hgap (indicatorOf A)
  rwa [varianceOf_indicator, dirichlet_indicator π P hrev] at this

/-- Relaxation time from SC-2: with a positive gap and a positive flow, the relaxation time `1/γ`
is at least the set's variance over its flow. -/
theorem relaxation_ge_of_cut (π : S → ℝ) (P : S → S → ℝ)
    (hrev : ∀ s t, π s * P s t = π t * P t s) {γ : ℝ} (hγ : 0 < γ)
    (hgap : ∀ f : S → ℝ, γ * varianceOf π f ≤ dirichlet π P f) (A : Finset S)
    (hflow : 0 < cutFlow π P A) :
    (∑ s ∈ A, π s) * (1 - ∑ s ∈ A, π s) / cutFlow π P A ≤ 1 / γ := by
  rw [div_le_div_iff₀ hflow hγ, one_mul, mul_comm]
  exact gap_le_cut π P hrev γ hgap A

/-- A stochastic kernel with stationary `π` has the same flow both ways across any cut, reversible
or not. -/
theorem cutFlow_compl_of_stationary (π : S → ℝ) (P : S → S → ℝ) (hP : ∀ s, ∑ t, P s t = 1)
    (hstat : ∀ t, ∑ s, π s * P s t = π t) (A : Finset S) :
    cutFlow π P Aᶜ = cutFlow π P A := by
  have hout : ∀ s, ∑ t, π s * P s t = π s := fun s => by rw [← Finset.mul_sum, hP, mul_one]
  have hA : ∑ s ∈ A, (∑ t ∈ A, π s * P s t + ∑ t ∈ Aᶜ, π s * P s t) = ∑ s ∈ A, π s :=
    Finset.sum_congr rfl fun s _ => by rw [Finset.sum_add_sum_compl, hout]
  have hin : ∑ t ∈ A, (∑ s ∈ A, π s * P s t + ∑ s ∈ Aᶜ, π s * P s t) = ∑ t ∈ A, π t :=
    Finset.sum_congr rfl fun t _ => by rw [Finset.sum_add_sum_compl, hstat]
  rw [Finset.sum_add_distrib] at hA hin
  rw [cutFlow, cutFlow, compl_compl]
  have hswap : ∑ t ∈ A, ∑ s ∈ Aᶜ, π s * P s t = ∑ s ∈ Aᶜ, ∑ t ∈ A, π s * P s t :=
    Finset.sum_comm
  have hswap' : ∑ t ∈ A, ∑ s ∈ A, π s * P s t = ∑ s ∈ A, ∑ t ∈ A, π s * P s t :=
    Finset.sum_comm
  linarith

/-- **SC-2 without reversibility.** The Dirichlet form of a stationary stochastic kernel is the
form of its additive symmetrization, and `1_A` still reads `A`'s flow. So a sweep's composition,
which SC-0 notes need not be reversible, keeps the bound. -/
theorem dirichlet_indicator_of_stationary (π : S → ℝ) (P : S → S → ℝ) (hP : ∀ s, ∑ t, P s t = 1)
    (hstat : ∀ t, ∑ s, π s * P s t = π t) (A : Finset S) :
    dirichlet π P (indicatorOf A) = cutFlow π P A := by
  rw [dirichlet_indicator_eq_half, cutFlow_compl_of_stationary π P hP hstat]
  ring

/-- **SC-2 for a stationary kernel.** The same bound, with `γ` the constant of the Dirichlet form
of `P`, which is that of its additive symmetrization. -/
theorem gap_le_cut_of_stationary (π : S → ℝ) (P : S → S → ℝ) (hP : ∀ s, ∑ t, P s t = 1)
    (hstat : ∀ t, ∑ s, π s * P s t = π t) (γ : ℝ)
    (hgap : ∀ f : S → ℝ, γ * varianceOf π f ≤ dirichlet π P f) (A : Finset S) :
    γ * ((∑ s ∈ A, π s) * (1 - ∑ s ∈ A, π s)) ≤ cutFlow π P A := by
  have := hgap (indicatorOf A)
  rwa [varianceOf_indicator, dirichlet_indicator_of_stationary π P hP hstat] at this

/-- **SC-5, the separator.** If every step out of `A` lands in `W` or stays in `A`, reversibility
puts
`A`'s flow through `W`: it is at most `W`'s mass times the largest chance of stepping from `W` into
`A`. Ronen et al.'s Lemma 3 is the case `W = {T₀}`, the trivial tree. -/
theorem cutFlow_le_separator (π : S → ℝ) (P : S → S → ℝ) (hπ : ∀ s, 0 ≤ π s)
    (hrev : ∀ s t, π s * P s t = π t * P t s) (A W : Finset S) (hWA : Disjoint W A) {e : ℝ}
    (hA : ∀ s ∈ A, ∀ t, t ∉ A → t ∉ W → P s t = 0) (he : ∀ t ∈ W, ∑ s ∈ A, P t s ≤ e) :
    cutFlow π P A ≤ e * ∑ t ∈ W, π t := by
  have hWsub : W ⊆ Aᶜ := fun t ht => Finset.mem_compl.2 (Finset.disjoint_right.1 hWA.symm ht)
  have hflow : cutFlow π P A = ∑ s ∈ A, ∑ t ∈ W, π s * P s t := by
    refine Finset.sum_congr rfl fun s hs => ?_
    refine (Finset.sum_subset hWsub fun t ht htW => ?_).symm
    rw [hA s hs t (Finset.mem_compl.1 ht) htW, mul_zero]
  rw [hflow, Finset.sum_comm, Finset.mul_sum]
  refine Finset.sum_le_sum fun t ht => ?_
  have hswap : ∑ s ∈ A, π s * P s t = π t * ∑ s ∈ A, P t s := by
    rw [Finset.mul_sum]
    exact Finset.sum_congr rfl fun s _ => hrev s t
  rw [hswap, mul_comm]
  exact mul_le_mul_of_nonneg_right (he t ht) (hπ t)

omit [DecidableEq S] in
/-- **SC-5.** So a Poincaré constant is at most `W`'s mass, times the exit chance, over `A`'s
variance:
the form a bound for BASS has to take, since no single state gates it. -/
theorem gap_le_separator (π : S → ℝ) (P : S → S → ℝ) (hπ : ∀ s, 0 ≤ π s)
    (hrev : ∀ s t, π s * P s t = π t * P t s) (γ : ℝ)
    (hgap : ∀ f : S → ℝ, γ * varianceOf π f ≤ dirichlet π P f) (A W : Finset S)
    (hWA : Disjoint W A) {e : ℝ} (hA : ∀ s ∈ A, ∀ t, t ∉ A → t ∉ W → P s t = 0)
    (he : ∀ t ∈ W, ∑ s ∈ A, P t s ≤ e) :
    γ * ((∑ s ∈ A, π s) * (1 - ∑ s ∈ A, π s)) ≤ e * ∑ t ∈ W, π t := by
  classical
  exact (gap_le_cut π P hrev γ hgap A).trans (cutFlow_le_separator π P hπ hrev A W hWA hA he)

end Flow

/-! ### SC-5: the superset path -/

section Superset

/-- A collapsed log marginal: an Occam cost `c` per column (`c = ½ log(1 + g)` under the g-form),
plus
a profile in `R²`. -/
def collapsedLogMarginal (c : ℝ) (h : ℝ → ℝ) (M : ℕ) (R2 : ℝ) : ℝ := -((M : ℝ) * c) + h R2

/-- **SC-5, the superset step.** Adding `k` columns cannot lower `R²`, so for a profile increasing
in
`R²` the collapsed log marginal falls by at most `k c`, whatever `N` and whatever variance the
columns
explain. Removing the first structure's terms afterwards is the same step read from the second
structure. So the path `S → S ∪ T → T` has no likelihood valley: the depth a tree's root split
forces
does not arise here. -/
theorem collapsedLogMarginal_superset_ge {c : ℝ} {h : ℝ → ℝ} (M k : ℕ) {R2 R2' : ℝ}
    (hmono : h R2 ≤ h R2') :
    collapsedLogMarginal c h M R2 - (k : ℝ) * c ≤ collapsedLogMarginal c h (M + k) R2' := by
  unfold collapsedLogMarginal
  push_cast
  linarith

end Superset

/-! ### SC-3 -/

section Precision

variable {m : Type*}

/-- **SC-3, one entry.** `G_ij A_ij = A_ij + A_ij / (g₀ g_i g_j)`. Here `g₀` stands for the paper's
`g₀²`, which only ever enters as a nonzero scalar. -/
theorem khaos_precision_entry (A : Matrix m m ℝ) (g : m → ℝ) (g₀ : ℝ) (i j : m)
    (hg₀ : g₀ ≠ 0) (hi : g i ≠ 0) (hj : g j ≠ 0) :
    (g₀ * g i * g j + 1) / (g₀ * g i * g j) * A i j = A i j + A i j / (g₀ * g i * g j) := by
  field_simp

/-- KHAOS's weight matrix `G_ij = (g₀ g_i g_j + 1) / (g₀ g_i g_j)`. -/
noncomputable def khaosWeights (g : m → ℝ) (g₀ : ℝ) : Matrix m m ℝ :=
  Matrix.of fun i j => (g₀ * g i * g j + 1) / (g₀ * g i * g j)

open Matrix

variable [Fintype m] [DecidableEq m]

/-- **SC-3, matrix form.** `G ∘ A = A + g₀⁻¹ D⁻¹ A D⁻¹` with `D = diag g`. -/
theorem khaos_precision_eq (A : Matrix m m ℝ) (g : m → ℝ) (g₀ : ℝ)
    (hg₀ : g₀ ≠ 0) (hg : ∀ i, g i ≠ 0) :
    khaosWeights g g₀ ⊙ A =
      A + g₀⁻¹ • (diagonal (fun i => (g i)⁻¹) * A * diagonal (fun i => (g i)⁻¹)) := by
  ext i j
  simp only [hadamard_apply, khaosWeights, of_apply, Matrix.add_apply, Matrix.smul_apply,
    smul_eq_mul, mul_diagonal, diagonal_mul]
  rw [khaos_precision_entry A g g₀ i j hg₀ (hg i) (hg j)]
  have hi := hg i
  have hj := hg j
  field_simp

/-- **SC-3, the conjugate form.** With `S₀ = g₀ D A⁻¹ D` for an invertible Gram `A`,
`G ∘ A = A + S₀⁻¹`: the paper's precision is the ordinary conjugate one. -/
theorem khaos_precision_eq_add_inv (A : Matrix m m ℝ) (g : m → ℝ) (g₀ : ℝ)
    (hA : IsUnit A.det) (hg₀ : g₀ ≠ 0) (hg : ∀ i, g i ≠ 0) :
    khaosWeights g g₀ ⊙ A = A + (g₀ • (diagonal g * A⁻¹ * diagonal g))⁻¹ := by
  rw [khaos_precision_eq A g g₀ hg₀ hg]
  congr 1
  have hD : diagonal g * diagonal (fun i => (g i)⁻¹) = (1 : Matrix m m ℝ) := by
    rw [diagonal_mul_diagonal, ← diagonal_one]
    exact congrArg diagonal (funext fun i => mul_inv_cancel₀ (hg i))
  have hD' : diagonal (fun i => (g i)⁻¹) * diagonal g = (1 : Matrix m m ℝ) := by
    rw [diagonal_mul_diagonal, ← diagonal_one]
    exact congrArg diagonal (funext fun i => inv_mul_cancel₀ (hg i))
  symm
  apply inv_eq_right_inv
  have hAA : A⁻¹ * A = 1 := nonsing_inv_mul A hA
  calc g₀ • (diagonal g * A⁻¹ * diagonal g) *
        (g₀⁻¹ • (diagonal (fun i => (g i)⁻¹) * A * diagonal (fun i => (g i)⁻¹)))
      = (g₀ * g₀⁻¹) • (diagonal g * A⁻¹ * (diagonal g * diagonal (fun i => (g i)⁻¹)) * A *
          diagonal (fun i => (g i)⁻¹)) := by
        rw [smul_mul_smul_comm]
        simp only [Matrix.mul_assoc]
    _ = diagonal g * (A⁻¹ * A) * diagonal (fun i => (g i)⁻¹) := by
        rw [mul_inv_cancel₀ hg₀, one_smul, hD]
        simp only [Matrix.mul_one, Matrix.mul_assoc]
    _ = 1 := by rw [hAA, Matrix.mul_one, hD]

set_option linter.unusedFintypeInType false in
set_option linter.unusedDecidableInType false in
/-- **SC-3, positive definiteness.** For a positive definite Gram and `g₀ > 0`, the precision
`G ∘ A` is positive definite, as a sum of two positive definite matrices: no Schur product needed.
-/
theorem khaos_precision_posDef (A : Matrix m m ℝ) (g : m → ℝ) (g₀ : ℝ)
    (hA : A.PosDef) (hg₀ : 0 < g₀) (hg : ∀ i, g i ≠ 0) :
    (khaosWeights g g₀ ⊙ A).PosDef := by
  rw [khaos_precision_eq A g g₀ hg₀.ne' hg]
  refine hA.add ?_
  have hinj : Function.Injective (diagonal (fun i => (g i)⁻¹) : Matrix m m ℝ).mulVec := by
    intro x y hxy
    funext i
    have := congrFun hxy i
    simp only [mulVec_diagonal] at this
    exact mul_left_cancel₀ (inv_ne_zero (hg i)) this
  have hconj : (diagonal (fun i => (g i)⁻¹) * A * diagonal (fun i => (g i)⁻¹) :
      Matrix m m ℝ).PosDef := by
    have := hA.conjTranspose_mul_mul_same (B := diagonal (fun i => (g i)⁻¹)) hinj
    simpa [diagonal_conjTranspose] using this
  exact hconj.smul (inv_pos.2 hg₀)

omit [Fintype m] [DecidableEq m] in
/-- With `g` constant the precision is a scalar times the Gram: `GForm`'s case. -/
theorem khaos_precision_of_const (A : Matrix m m ℝ) (c g₀ : ℝ) :
    khaosWeights (fun _ => c) g₀ ⊙ A = ((g₀ * c * c + 1) / (g₀ * c * c)) • A := by
  ext i j
  simp [hadamard_apply, khaosWeights]

omit [Fintype m] in
/-- With a diagonal Gram the precision is diagonal, for any `g`. So SC-3's "a third form" is
about a Gram with off-diagonal entries: on orthogonal columns it is a Gram plus a diagonal. -/
theorem khaos_precision_of_diagonal (d g : m → ℝ) (g₀ : ℝ) :
    khaosWeights g g₀ ⊙ diagonal d =
      diagonal fun i => (g₀ * g i * g i + 1) / (g₀ * g i * g i) * d i := by
  ext i j
  by_cases hij : i = j
  · subst hij
    simp [hadamard_apply, khaosWeights]
  · simp [hadamard_apply, khaosWeights, hij]

end Precision

end BASS
