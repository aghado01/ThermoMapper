import Mathlib.Analysis.Complex.Basic
import Mathlib.Combinatorics.SimpleGraph.Connectivity.Connected
import Mathlib.Data.ZMod.Basic
import Mathlib.Tactic.FinCases
import Mathlib.Tactic.LinearCombination
import Mathlib.Tactic.NormNum

/-!
# Phase propagation on a connected graph (S1)

The mechanism of Bodmann–Hammen–Voroninski's structured-pattern recovery (arXiv 1109.0573,
Theorem `teo:1D`) in general form. Signals `x, y : V → ℂ` live on the vertices of a simple graph,
`x` nowhere zero. Magnitudes on the vertices and the products `x u · conj (x v)` on the edges
determine `y` up to one global phase when the graph is connected, and only then.

* `exists_phase_of_vertex`: on a preconnected graph, the edge products and the magnitude at a
  single vertex give `y = c · x` with `‖c‖ = 1`.
* `exists_phase`: S1 as posed, on a connected graph with every magnitude.
* `exists_phase_iff_preconnected`: with every magnitude given, the edge products determine `y`
  up to a global phase for every `y` exactly when the graph is preconnected. A split into two
  parts with no edge between them admits `y = x` on one part and `y = −x` on the other.
* `mul_conj_eq_polarization` and `mul_conj_eq_polarization_sub`: each edge product is read from
  magnitudes, `‖a + b‖` and `‖a + i·b‖` (or the paper's `‖a − i·b‖`) beside `‖a‖` and `‖b‖`;
  `exists_phase_of_norms` states S1 on magnitudes alone.
* `shiftGraph_connected_iff`: on `ZMod n` with edges `k ~ k + s`, the graph is connected exactly
  when `s` is coprime to `n`, the paper's condition.
-/

open ComplexConjugate

namespace TimeFrequency

variable {V : Type*} {G : SimpleGraph V}

/-- One edge carries the phase: if `yᵤ = c·xᵤ` with `‖c‖ = 1` and `yᵤ·conj y_w = xᵤ·conj x_w`,
then `y_w = c·x_w`. -/
theorem phase_step {c xu xw yu yw : ℂ} (hc : ‖c‖ = 1) (hxu : xu ≠ 0) (hu : yu = c * xu)
    (he : yu * conj yw = xu * conj xw) : yw = c * xw := by
  rw [hu, mul_comm c, mul_assoc] at he
  have h1 : c * conj yw = conj xw := mul_left_cancel₀ hxu he
  have h2 : conj c * yw = xw := by simpa using congrArg conj h1
  have hcc : c * conj c = 1 := by
    rw [Complex.mul_conj', hc]
    simp
  calc yw = (c * conj c) * yw := by rw [hcc, one_mul]
    _ = c * (conj c * yw) := by ring
    _ = c * xw := by rw [h2]

/-- **S1, from one vertex.** On a preconnected graph, with `x` nowhere zero, the edge products
and the magnitude at one vertex `v₀` determine `y` up to a global phase. -/
theorem exists_phase_of_vertex (hG : G.Preconnected) {x y : V → ℂ} (hx : ∀ v, x v ≠ 0)
    (v₀ : V) (h₀ : ‖y v₀‖ = ‖x v₀‖)
    (hedge : ∀ u v, G.Adj u v → y u * conj (y v) = x u * conj (x v)) :
    ∃ c : ℂ, ‖c‖ = 1 ∧ y = c • x := by
  have hx₀ : ‖x v₀‖ ≠ 0 := norm_ne_zero_iff.2 (hx v₀)
  refine ⟨y v₀ / x v₀, by rw [norm_div, h₀, div_self hx₀], ?_⟩
  have key : ∀ {u v : V} (p : G.Walk u v),
      y u = y v₀ / x v₀ * x u → y v = y v₀ / x v₀ * x v := by
    intro u v p
    induction p with
    | nil => exact id
    | cons h _ ih =>
      exact fun hu => ih (phase_step (by rw [norm_div, h₀, div_self hx₀]) (hx _) hu
        (hedge _ _ h))
  funext v
  exact (hG v₀ v).elim fun p => key p (div_mul_cancel₀ _ (hx v₀)).symm

/-- **S1 as posed.** On a connected graph, with `x` nowhere zero, equal magnitudes at every
vertex and equal products on every edge give `y = c·x` with `‖c‖ = 1`. -/
theorem exists_phase (hG : G.Connected) {x y : V → ℂ} (hx : ∀ v, x v ≠ 0)
    (hnorm : ∀ v, ‖y v‖ = ‖x v‖)
    (hedge : ∀ u v, G.Adj u v → y u * conj (y v) = x u * conj (x v)) :
    ∃ c : ℂ, ‖c‖ = 1 ∧ y = c • x :=
  have ⟨v₀⟩ := hG.nonempty
  exists_phase_of_vertex hG.preconnected hx v₀ (hnorm v₀) hedge

/-- **S1, sharpness.** For `x` nowhere zero, the magnitudes and edge products determine every
`y` up to a global phase exactly when the graph is preconnected. -/
theorem exists_phase_iff_preconnected {x : V → ℂ} (hx : ∀ v, x v ≠ 0) :
    (∀ y : V → ℂ, (∀ v, ‖y v‖ = ‖x v‖) →
      (∀ u v, G.Adj u v → y u * conj (y v) = x u * conj (x v)) →
        ∃ c : ℂ, ‖c‖ = 1 ∧ y = c • x) ↔ G.Preconnected := by
  classical
  constructor
  · intro H
    by_contra hpre
    simp only [SimpleGraph.Preconnected, not_forall] at hpre
    obtain ⟨u, v, huv⟩ := hpre
    have hadj : ∀ a b, G.Adj a b → (G.Reachable u a ↔ G.Reachable u b) := fun a b h =>
      ⟨fun hr => hr.trans h.reachable, fun hr => hr.trans h.symm.reachable⟩
    obtain ⟨c, -, hc⟩ := H (fun w => if G.Reachable u w then x w else -x w)
      (fun w => by by_cases h : G.Reachable u w <;> simp [h])
      (fun a b h => by
        by_cases ha : G.Reachable u a
        · simp [ha, (hadj a b h).1 ha]
        · have hb : ¬G.Reachable u b := fun hb => ha ((hadj a b h).2 hb)
          simp [ha, hb])
    have hu := congrFun hc u
    have hv := congrFun hc v
    simp only [Pi.smul_apply, smul_eq_mul, if_pos (SimpleGraph.Reachable.refl u)] at hu
    simp only [Pi.smul_apply, smul_eq_mul, if_neg huv] at hv
    have hc1 : c = 1 := (mul_right_cancel₀ (hx u) (by rw [one_mul]; exact hu)).symm
    rw [hc1, one_mul] at hv
    exact hx v (by linear_combination (-1 / 2 : ℂ) * hv)
  · intro hG y hnorm hedge
    rcases isEmpty_or_nonempty V with hV | hV
    · exact ⟨1, by simp, funext fun v => isEmptyElim v⟩
    · obtain ⟨v₀⟩ := hV
      exact exists_phase_of_vertex hG hx v₀ (hnorm v₀) hedge

/-- **Polarization.** `a·conj b = ½[(‖a + b‖² − ‖a‖² − ‖b‖²) + i(‖a + i·b‖² − ‖a‖² − ‖b‖²)]`. -/
theorem mul_conj_eq_polarization (a b : ℂ) :
    a * conj b = (((‖a + b‖ ^ 2 - ‖a‖ ^ 2 - ‖b‖ ^ 2 : ℝ) : ℂ) +
      ((‖a + Complex.I * b‖ ^ 2 - ‖a‖ ^ 2 - ‖b‖ ^ 2 : ℝ) : ℂ) * Complex.I) / 2 := by
  apply Complex.ext <;> simp [Complex.sq_norm, Complex.normSq_apply] <;> ring

/-- **Polarization, the paper's sign.** With `a − i·b` in place of `a + i·b` the imaginary part
changes sign: `a·conj b = ½[(‖a + b‖² − ‖a‖² − ‖b‖²) − i(‖a − i·b‖² − ‖a‖² − ‖b‖²)]`. -/
theorem mul_conj_eq_polarization_sub (a b : ℂ) :
    a * conj b = (((‖a + b‖ ^ 2 - ‖a‖ ^ 2 - ‖b‖ ^ 2 : ℝ) : ℂ) -
      ((‖a - Complex.I * b‖ ^ 2 - ‖a‖ ^ 2 - ‖b‖ ^ 2 : ℝ) : ℂ) * Complex.I) / 2 := by
  apply Complex.ext <;> simp [Complex.sq_norm, Complex.normSq_apply] <;> ring

/-- **S1 on magnitudes alone.** On a connected graph, with `x` nowhere zero, equal magnitudes
at every vertex and of `y u + y v` and `y u + i·y v` on every edge give `y = c·x`, `‖c‖ = 1`. -/
theorem exists_phase_of_norms (hG : G.Connected) {x y : V → ℂ} (hx : ∀ v, x v ≠ 0)
    (hnorm : ∀ v, ‖y v‖ = ‖x v‖)
    (hsum : ∀ u v, G.Adj u v → ‖y u + y v‖ = ‖x u + x v‖)
    (hrot : ∀ u v, G.Adj u v → ‖y u + Complex.I * y v‖ = ‖x u + Complex.I * x v‖) :
    ∃ c : ℂ, ‖c‖ = 1 ∧ y = c • x :=
  exists_phase hG hx hnorm fun u v h => by
    rw [mul_conj_eq_polarization, mul_conj_eq_polarization, hsum u v h, hrot u v h, hnorm u,
      hnorm v]

/-! ### The shift graph on `ZMod n` -/

/-- The graph on `ZMod n` with edges `k ~ k + s`. -/
def shiftGraph (n : ℕ) (s : ZMod n) : SimpleGraph (ZMod n) :=
  SimpleGraph.fromRel fun a b => b = a + s

theorem shiftGraph_reachable_nsmul (n : ℕ) (s : ZMod n) (k : ℕ) :
    (shiftGraph n s).Reachable 0 (k • s) := by
  induction k with
  | zero =>
    rw [zero_nsmul]
  | succ k ih =>
    by_cases hs : s = 0
    · subst hs
      simp only [smul_zero] at ih ⊢
      exact ih
    · refine ih.trans (SimpleGraph.Adj.reachable ?_)
      rw [shiftGraph, SimpleGraph.fromRel_adj, succ_nsmul]
      exact ⟨fun h => hs (by simpa using h), Or.inl rfl⟩

/-- Every vertex reachable from `0` is an integer multiple of the shift. -/
theorem shiftGraph_reachable_zmultiple (n : ℕ) (s : ZMod n) {v : ZMod n}
    (h : (shiftGraph n s).Reachable 0 v) : ∃ k : ℤ, v = k • s := by
  have key : ∀ {u w : ZMod n} (p : (shiftGraph n s).Walk u w),
      (∃ k : ℤ, u = k • s) → ∃ k : ℤ, w = k • s := by
    intro u w p
    induction p with
    | nil => exact id
    | cons hadj _ ih =>
      rintro ⟨k, rfl⟩
      rw [shiftGraph, SimpleGraph.fromRel_adj] at hadj
      rcases hadj.2 with h | h
      · exact ih ⟨k + 1, by rw [h, add_smul, one_smul]⟩
      · exact ih ⟨k - 1, by rw [sub_smul, one_smul, h, add_sub_cancel_right]⟩
  exact h.elim fun p => key p ⟨0, by simp⟩

/-- **S1's corollary.** On `ZMod n`, `n ≥ 1`, the shift graph `k ~ k + s` is connected exactly
when `s` is coprime to `n`. -/
theorem shiftGraph_connected_iff (n : ℕ) [NeZero n] (s : ℕ) :
    (shiftGraph n s).Connected ↔ s.Coprime n := by
  rw [← ZMod.isUnit_iff_coprime]
  constructor
  · intro hG
    obtain ⟨k, hk⟩ := shiftGraph_reachable_zmultiple n s (hG.preconnected 0 1)
    rw [zsmul_eq_mul] at hk
    exact IsUnit.of_mul_eq_one_right (k : ZMod n) hk.symm
  · intro hs
    refine ⟨fun a b => ?_⟩
    have hreach : ∀ v, (shiftGraph n s).Reachable 0 v := fun v => by
      have h := shiftGraph_reachable_nsmul n (s : ZMod n) (v * ↑hs.unit⁻¹).val
      rwa [nsmul_eq_mul, ZMod.natCast_zmod_val, mul_assoc, hs.val_inv_mul, mul_one] at h
    exact (hreach a).symm.trans (hreach b)

/-! ### Witnesses -/

/-- The path `0 — 1 — 2` with `x = (1, 2, 3)` and `y = i·x`: S1's hypotheses hold and its
conclusion follows. -/
example : ∃ c : ℂ, ‖c‖ = 1 ∧ (fun v : ZMod 3 => Complex.I * ((v.val : ℂ) + 1)) =
    c • fun v : ZMod 3 => (v.val : ℂ) + 1 := by
  refine exists_phase ((shiftGraph_connected_iff 3 1).2 (Nat.coprime_one_left 3)) (fun v => ?_)
    (fun v => ?_) (fun u v _ => ?_)
  · exact_mod_cast Nat.succ_ne_zero v.val
  · simp
  · simp only [map_mul, Complex.conj_I]
    linear_combination (-(((u.val : ℂ) + 1) * conj ((v.val : ℂ) + 1))) * Complex.I_sq

/-- Sharpness on the edgeless graph on two vertices: `y = (1, −1)` and `x = (1, 1)` meet both
hypotheses, and no global phase relates them. -/
example : (∀ v, ‖(![1, -1] : Fin 2 → ℂ) v‖ = ‖(![1, 1] : Fin 2 → ℂ) v‖) ∧
    (∀ u v, (⊥ : SimpleGraph (Fin 2)).Adj u v →
      (![1, -1] : Fin 2 → ℂ) u * conj ((![1, -1] : Fin 2 → ℂ) v) =
        (![1, 1] : Fin 2 → ℂ) u * conj ((![1, 1] : Fin 2 → ℂ) v)) ∧
    ¬∃ c : ℂ, (![1, -1] : Fin 2 → ℂ) = c • ![1, 1] := by
  refine ⟨fun v => by fin_cases v <;> simp, fun u v h => absurd h (by simp), ?_⟩
  rintro ⟨c, hc⟩
  have h0 := congrFun hc 0
  have h1 := congrFun hc 1
  simp only [Pi.smul_apply, smul_eq_mul, Matrix.cons_val_zero, Matrix.cons_val_one,
    mul_one] at h0 h1
  rw [← h0] at h1
  norm_num at h1

end TimeFrequency
