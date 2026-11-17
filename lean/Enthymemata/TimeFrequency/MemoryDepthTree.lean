import Enthymemata.TimeFrequency.MemoryDepthBounds
import Enthymemata.TimeFrequency.MemoryDepthPaired

/-! # Exact suffix-tree recovery from paired powers and separated transitions -/

namespace TimeFrequency.ContextTree.DepthEvidence

def constantSuffix {d : ℕ} (q : Word d → ℝ) (s : List Bool) : Prop :=
  ∀ h h' : Word d, s <:+ List.ofFn h → s <:+ List.ofFn h' → q h = q h'

def suffixTree {d : ℕ} (C : List Bool → Prop) : Set (List Bool) :=
  {s | (∃ h : Word d, s <:+ List.ofFn h) ∧ C s ∧
    ∀ u : List Bool, u <:+ s → u.length < s.length → ¬ C u}

def separated {d : ℕ} (q : Word d → ℝ) (γ : ℝ) : Prop :=
  ∀ s, ¬ constantSuffix q s → ∃ h h' : Word d,
    s <:+ List.ofFn h ∧ s <:+ List.ofFn h' ∧ γ ≤ |q h - q h'|

def empiricalConstant {d : ℕ} (qhat : Word d → ℝ) (γ : ℝ) (s : List Bool) : Prop :=
  ∀ h h' : Word d, s <:+ List.ofFn h → s <:+ List.ofFn h' →
    |qhat h - qhat h'| ≤ γ / 2

def binaryKernel {d : ℕ} (q : Word d → ℝ) (h : Word d) (b : Bool) : ℝ :=
  if b then q h else 1 - q h

theorem scalar_minimalTree {d : ℕ} (q : Word d → ℝ) :
    minimalTree (binaryKernel q) = suffixTree (d := d) (constantSuffix q) := by
  have he (s : List Bool) :
      (∀ h h' : Word d, s <:+ List.ofFn h → s <:+ List.ofFn h' →
        binaryKernel q h = binaryKernel q h') ↔ constantSuffix q s := by
    constructor
    · intro hc h h' hh hh'
      have hx := congrFun (hc h h' hh hh') true
      exact hx
    · intro hc h h' hh hh'
      funext b
      simp only [binaryKernel, hc h h' hh hh']
  ext s
  simp only [minimalTree, suffixTree, Set.mem_setOf_eq, he]

theorem conditional_kernel {d : ℕ} (P : Law (Word (d + 1))) :
    binaryKernel (conditional (historyResponseLaw P)) = nextKernel d P.mass := by
  funext h b
  have hp := historyMass_pos (historyResponseLaw P) h
  cases b
  · change 1 - P.mass (Fin.snoc h true) /
      (P.mass (Fin.snoc h false) + P.mass (Fin.snoc h true)) =
      P.mass (Fin.snoc h false) /
        (P.mass (Fin.snoc h false) + P.mass (Fin.snoc h true))
    change 0 < P.mass (Fin.snoc h false) + P.mass (Fin.snoc h true) at hp
    field_simp [hp.ne']
    ring
  · rfl

theorem threshold_contexts {d : ℕ} {q qhat : Word d → ℝ} {γ : ℝ}
    (hsep : separated q γ) (herr : ∀ h, |qhat h - q h| < γ / 4) (s : List Bool) :
    empiricalConstant qhat γ s ↔ constantSuffix q s := by
  constructor
  · intro hc
    by_contra hn
    obtain ⟨h, h', hh, hh', hgap⟩ := hsep s hn
    have he := hc h h' hh hh'
    have h1 := abs_sub_le (q h) (qhat h) (q h')
    have h2 := abs_sub_le (qhat h) (qhat h') (q h')
    have e1 := herr h
    have e2 := herr h'
    rw [abs_sub_comm (q h) (qhat h)] at h1
    linarith
  · intro hc h h' hh hh'
    have h1 := abs_sub_le (qhat h) (q h) (qhat h')
    rw [hc h h' hh hh', abs_sub_comm (q h') (qhat h')] at h1
    have e1 := herr h
    have e2 := herr h'
    rw [hc h h' hh hh'] at e1
    linarith

theorem threshold_tree {d : ℕ} {q qhat : Word d → ℝ} {γ : ℝ}
    (hsep : separated q γ) (herr : ∀ h, |qhat h - q h| < γ / 4) :
    suffixTree (d := d) (empiricalConstant qhat γ) = minimalTree (binaryKernel q) := by
  rw [scalar_minimalTree]
  have he : empiricalConstant qhat γ = constantSuffix q :=
    funext (fun s ↦ propext (threshold_contexts hsep herr s))
  rw [he]

noncomputable def estimatedTree {d n : ℕ} (x : Fin n → Word d × Bool) (γ : ℝ) :=
  suffixTree (d := d) (empiricalConstant (estimate x) γ)

theorem paired_tree_budget {d : ℕ} (P : Law (Word (d + 1))) (n : ℕ) {π γ δ : ℝ}
    (hπ : 0 < π) (hm : ∀ h, π ≤ historyMass (historyResponseLaw P) h)
    (hγ : 0 < γ) (hu : γ ≤ 1) (hδ : 0 < δ)
    (hsep : separated (conditional (historyResponseLaw P)) γ)
    (hn : ⌈16 * Real.log ((2 : ℝ) ^ (d + 1) / δ) / (π * γ ^ 2)⌉₊ ≤ n) :
    probability (iid n (historyResponseLaw P))
      (fun x ↦ estimatedTree x γ ≠ minimalTree (nextKernel d P.mass)) ≤ δ := by
  have hcount : Fintype.card (Word d) = 2 ^ d := by simp [Word]
  have he : Real.log (2 * Fintype.card (Word d) / δ) / (π * (γ / 4) ^ 2) =
      16 * Real.log ((2 : ℝ) ^ (d + 1) / δ) / (π * γ ^ 2) := by
    rw [hcount]
    push_cast
    rw [pow_succ, mul_comm (2 : ℝ) (2 ^ d)]
    ring
  have hb := transition_budget (historyResponseLaw P) n
    (show 0 < γ / 4 by positivity) (show γ / 4 ≤ 1 / 2 by linarith)
    hπ hm hδ (by rw [hcount]; positivity) (by rwa [he])
  apply (probability_mono _ (E := fun x ↦ estimatedTree x γ ≠
    minimalTree (nextKernel d P.mass)) (F := fun x ↦
      ∃ h, γ / 4 ≤ |estimate x h - conditional (historyResponseLaw P) h|) ?_).trans hb
  intro x hx
  by_contra hn
  have herr : ∀ h, |estimate x h - conditional (historyResponseLaw P) h| < γ / 4 := by
    intro h
    exact lt_of_not_ge (fun he ↦ hn ⟨h, he⟩)
  have ht := threshold_tree hsep herr
  rw [conditional_kernel] at ht
  exact hx ht

#print axioms scalar_minimalTree
#print axioms conditional_kernel
#print axioms threshold_contexts
#print axioms threshold_tree
#print axioms paired_tree_budget

end TimeFrequency.ContextTree.DepthEvidence
