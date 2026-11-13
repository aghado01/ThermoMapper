import Enthymemata.BassStructureChain
import Mathlib.Algebra.BigOperators.Fin
import Mathlib.Analysis.InnerProductSpace.PiL2
import Mathlib.Analysis.SpecialFunctions.Pow.Real
import Mathlib.Data.Fin.Tuple.Basic
import Mathlib.Tactic.FieldSimp
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.Ring

/-!
# Active subspaces — E36's direction moves (AS-3)

The statements of `Protolemmata/BASS/active-subspace-lemmas.md`, AS-3: the acceptance ratios that
projection-pursuit ridges on the BASS core (BPPR) would code. They reuse SC-1's `BASS.mhRatio` and
`BASS.mh_detailed_balance`.

* **(a) The change move.** `powerSpherical_symm`: the power spherical kernel reads `⟪μ, x⟫` alone,
  so on one sphere `q(θ′ | θ) = q(θ | θ′)`. `mhRatio_of_symm`: a symmetric proposal leaves the
  target ratio. `change_ratio`, `change_detailed_balance`: with the direction's prior uniform,
  the knot drawn from its prior given the proposed direction, the ratio is the likelihood's alone.
  `change_ratio_powerSpherical` instantiates it on the sphere of the active coordinates.
  `powerSpherical_antipode`: the law puts nothing on `−μ`, so on `S⁰` the move never flips a sign.
* **(b) Birth and death,** on ordered tuples with insertion and deletion at uniform positions, the
  paper's form. `birth_green_ratio`: Green's ratio is `L′/L · w(M+1)/w(M) · P_D/P_B · π_r/q_r`; the
  positions' `1/(M+1)` cancel. `poisson_count_ratio`: a Poisson count gives `λ/(M+1)`.
  `double_count_ratio`: a further `1/(M+1)` is the ratio of the count weight `λ^M/(M!)²`, so a
  chain that writes it targets that prior and not Poisson's. `prior_draw_cancel`: a direction and
  knot drawn from their priors leave only the `(a, 𝒥)` terms. `birth_death_balance`: each birth
  and its reverse death balance, through `mh_detailed_balance`.

The measure-theoretic form on `S^{a−1} × ℝ` is not stated: the finite statements carry the ratio's
algebra, and the numeric check AS-3d reads the continuous form.
-/

open Finset BASS

namespace ActiveSubspace

/-! ### (a) The change move -/

/-- With a symmetric proposal, Metropolis–Hastings' ratio is the target's. -/
theorem mhRatio_of_symm {S : Type*} (π : S → ℝ) (q : S → S → ℝ) (s t : S) (h : q s t = q t s)
    (hq : q s t ≠ 0) : mhRatio π q s t = π t / π s := by
  unfold mhRatio
  rw [h, mul_div_mul_right _ _ (h ▸ hq)]

/-- The power spherical kernel `(1 + ⟪μ, x⟫)^κ / N` (De Cao and Aziz), with `N = N(κ, a)` its
normalizer, which does not read `μ`. -/
noncomputable def powerSpherical {a : ℕ} (κ N : ℝ) (μ x : EuclideanSpace ℝ (Fin a)) : ℝ :=
  (1 + inner ℝ μ x) ^ κ / N

/-- **AS-3(a).** The kernel is symmetric in the pair: `q(x | μ) = q(μ | x)`. -/
theorem powerSpherical_symm {a : ℕ} (κ N : ℝ) (μ x : EuclideanSpace ℝ (Fin a)) :
    powerSpherical κ N μ x = powerSpherical κ N x μ := by
  simp only [powerSpherical, real_inner_comm]

/-- **AS-3(a), at `a = 1`.** The law puts nothing on the antipode: on `S⁰ = {μ, −μ}` the change move
keeps the direction's sign. -/
theorem powerSpherical_antipode {a : ℕ} (κ N : ℝ) (hκ : κ ≠ 0) (μ : EuclideanSpace ℝ (Fin a))
    (hμ : ‖μ‖ = 1) : powerSpherical κ N μ (-μ) = 0 := by
  simp only [powerSpherical, inner_neg_right, real_inner_self_eq_norm_sq, hμ, one_pow,
    add_neg_cancel, Real.zero_rpow hκ, zero_div]

/-- **AS-3(a), the change move's ratio.** States are a direction and a knot. The direction's prior
is uniform (the constant `c`), the knot's prior `πt θ` may read the direction, and `L` is the
collapsed likelihood. The proposal draws the direction from a symmetric kernel `k` and the knot from
its prior given the proposed direction. Then the ratio is the likelihood ratio alone. -/
theorem change_ratio {Θ T : Type*} (c : ℝ) (k : Θ → Θ → ℝ) (hk : ∀ a b, k a b = k b a)
    (πt L : Θ → T → ℝ) (s s' : Θ × T) (hc : c ≠ 0) (hπ : πt s.1 s.2 ≠ 0) (hπ' : πt s'.1 s'.2 ≠ 0)
    (hk0 : k s.1 s'.1 ≠ 0) (hL : L s.1 s.2 ≠ 0) :
    mhRatio (fun s => c * πt s.1 s.2 * L s.1 s.2) (fun s s' => k s.1 s'.1 * πt s'.1 s'.2) s s' =
      L s'.1 s'.2 / L s.1 s.2 := by
  unfold mhRatio
  dsimp only
  rw [hk s'.1 s.1]
  field_simp

/-- **AS-3(a), balance.** The change move with the likelihood ratio as its acceptance satisfies
detailed balance. -/
theorem change_detailed_balance {Θ T : Type*} (c : ℝ) (k : Θ → Θ → ℝ) (hk : ∀ a b, k a b = k b a)
    (πt L : Θ → T → ℝ) (s s' : Θ × T) (hc : 0 < c) (hπ : 0 < πt s.1 s.2) (hπ' : 0 < πt s'.1 s'.2)
    (hk0 : 0 < k s.1 s'.1) (hL : 0 < L s.1 s.2) (hL' : 0 < L s'.1 s'.2) :
    c * πt s.1 s.2 * L s.1 s.2 * (k s.1 s'.1 * πt s'.1 s'.2 * min 1 (L s'.1 s'.2 / L s.1 s.2)) =
      c * πt s'.1 s'.2 * L s'.1 s'.2 *
        (k s'.1 s.1 * πt s.1 s.2 * min 1 (L s.1 s.2 / L s'.1 s'.2)) := by
  have hk0' : 0 < k s'.1 s.1 := hk s.1 s'.1 ▸ hk0
  have h := mh_detailed_balance (fun s => c * πt s.1 s.2 * L s.1 s.2)
    (fun s s' => k s.1 s'.1 * πt s'.1 s'.2) s s' (by positivity) (by positivity) (by positivity)
    (by positivity)
  rw [change_ratio c k hk πt L s s' hc.ne' hπ.ne' hπ'.ne' hk0.ne' hL.ne',
    change_ratio c k hk πt L s' s hc.ne' hπ'.ne' hπ.ne' hk0'.ne' hL'.ne'] at h
  exact h

/-- **AS-3(a), with the power spherical proposal.** On the sphere of the active coordinates, the
change move that proposes the direction from the power spherical law about the current one, and the
knot from its prior given the proposed direction, accepts on the likelihood ratio alone. -/
theorem change_ratio_powerSpherical {a : ℕ} {T : Type*} (κ N c : ℝ)
    (πt L : Metric.sphere (0 : EuclideanSpace ℝ (Fin a)) 1 → T → ℝ)
    (s s' : Metric.sphere (0 : EuclideanSpace ℝ (Fin a)) 1 × T) (hc : c ≠ 0)
    (hπ : πt s.1 s.2 ≠ 0) (hπ' : πt s'.1 s'.2 ≠ 0)
    (hk0 : powerSpherical κ N (s.1 : EuclideanSpace ℝ (Fin a)) s'.1 ≠ 0) (hL : L s.1 s.2 ≠ 0) :
    mhRatio (fun s => c * πt s.1 s.2 * L s.1 s.2)
        (fun s s' => powerSpherical κ N (s.1 : EuclideanSpace ℝ (Fin a)) s'.1 * πt s'.1 s'.2) s s' =
      L s'.1 s'.2 / L s.1 s.2 :=
  change_ratio (Θ := Metric.sphere (0 : EuclideanSpace ℝ (Fin a)) 1) c
    (fun x y => powerSpherical κ N (x : EuclideanSpace ℝ (Fin a)) (y : EuclideanSpace ℝ (Fin a)))
    (fun x y => powerSpherical_symm κ N (x : EuclideanSpace ℝ (Fin a)) y) πt L s s' hc hπ hπ' hk0 hL

/-! ### (b) Birth and death -/

section BirthDeath

variable {R : Type*}

/-- The target on ordered tuples: a count weight `w M` (the count prior), independent ridges with
prior `pr`, and the collapsed likelihood `L`. -/
noncomputable def tupleTarget (w : ℕ → ℝ) (pr : R → ℝ) (L : (M : ℕ) → (Fin M → R) → ℝ) (M : ℕ)
    (c : Fin M → R) : ℝ :=
  w M * (∏ j, pr (c j)) * L M c

/-- The tuple with `r` inserted at position `i`. -/
def insertAt {M : ℕ} (i : Fin (M + 1)) (r : R) (c : Fin M → R) : Fin (M + 1) → R :=
  Fin.insertNth (α := fun _ => R) i r c

/-- Inserting `r` at `i` multiplies the ridges' prior by `pr r`. -/
theorem prod_insertAt (pr : R → ℝ) {M : ℕ} (c : Fin M → R) (i : Fin (M + 1)) (r : R) :
    ∏ j, pr (insertAt i r c j) = pr r * ∏ j, pr (c j) := by
  rw [Fin.prod_univ_succAbove _ i]
  simp only [insertAt, Fin.insertNth_apply_same, Fin.insertNth_apply_succAbove]

/-- **AS-3(b), Green's ratio for a birth.** A birth from the tuple `c` of `M` ridges chooses the
birth move with probability `P_B(M)`, a position among `M + 1` uniformly, and the ridge `r` from
`q`; its reverse death chooses the death move with probability `P_D(M + 1)` and that position
uniformly. The positions' `1/(M + 1)` cancel, and the ratio is
`L′/L · w(M+1)/w(M) · P_D(M+1)/P_B(M) · π_r(r)/q(r)`. -/
theorem birth_green_ratio (w : ℕ → ℝ) (pr q : R → ℝ) (L : (M : ℕ) → (Fin M → R) → ℝ)
    (PB PD : ℕ → ℝ) {M : ℕ} (c : Fin M → R) (i : Fin (M + 1)) (r : R) (hw : w M ≠ 0)
    (hpr : ∏ j, pr (c j) ≠ 0) (hL : L M c ≠ 0) (hq : q r ≠ 0) (hPB : PB M ≠ 0) :
    tupleTarget w pr L (M + 1) (insertAt i r c) * (PD (M + 1) / (M + 1)) /
        (tupleTarget w pr L M c * (PB M / (M + 1) * q r)) =
      L (M + 1) (insertAt i r c) / L M c * (w (M + 1) / w M) * (PD (M + 1) / PB M) *
        (pr r / q r) := by
  simp only [tupleTarget, prod_insertAt]
  have hM : ((M : ℝ) + 1) ≠ 0 := by positivity
  field_simp

/-- **AS-3(b).** A Poisson count, `w M = λ^M / M!`, gives the factor `λ/(M + 1)`. -/
theorem poisson_count_ratio (lam : ℝ) (hlam : lam ≠ 0) (M : ℕ) :
    (lam ^ (M + 1) / (M + 1).factorial) / (lam ^ M / M.factorial) = lam / (M + 1) := by
  rw [Nat.factorial_succ, Nat.cast_mul, Nat.cast_succ, pow_succ]
  have h1 : (M.factorial : ℝ) ≠ 0 := by positivity
  have h2 : ((M : ℝ) + 1) ≠ 0 := by positivity
  field_simp

/-- **AS-3(b), the count read twice.** `λ/(M + 1) · 1/(M + 1)` is the factor of the count weight
`λ^M / (M!)²`: a birth ratio that adds a further `1/(M + 1)` targets that count prior, not
Poisson's. -/
theorem double_count_ratio (lam : ℝ) (hlam : lam ≠ 0) (M : ℕ) :
    (lam ^ (M + 1) / ((M + 1).factorial : ℝ) ^ 2) / (lam ^ M / (M.factorial : ℝ) ^ 2) =
      lam / (M + 1) * (1 / (M + 1)) := by
  rw [Nat.factorial_succ, Nat.cast_mul, Nat.cast_succ, pow_succ]
  have h1 : (M.factorial : ℝ) ≠ 0 := by positivity
  have h2 : ((M : ℝ) + 1) ≠ 0 := by positivity
  field_simp

/-- **AS-3(b).** With the direction and knot drawn from their priors, their terms cancel and the
ratio keeps only the active count's and the index set's. -/
theorem prior_draw_cancel (pa pJ qa qJ pθ pt : ℝ) (hθ : pθ ≠ 0) (ht : pt ≠ 0) :
    (pa * pJ * pθ * pt) / (qa * qJ * pθ * pt) = (pa * pJ) / (qa * qJ) := by
  rw [mul_assoc (pa * pJ), mul_assoc (qa * qJ), mul_div_mul_right _ _ (mul_ne_zero hθ ht)]

/-- **AS-3(b), balance.** Each birth (at position `i`, ridge `r`) and its reverse death balance
with the acceptance `min 1 ρ` forward and `min 1 ρ⁻¹` back, `ρ` the ratio above. Summed over the
matched pairs, this is the detailed balance of the birth–death kernel. -/
theorem birth_death_balance (w : ℕ → ℝ) (pr q : R → ℝ) (L : (M : ℕ) → (Fin M → R) → ℝ)
    (PB PD : ℕ → ℝ) {M : ℕ} (c : Fin M → R) (i : Fin (M + 1)) (r : R) (hw : 0 < w M)
    (hw' : 0 < w (M + 1)) (hpr : ∀ x, 0 < pr x) (hL : 0 < L M c)
    (hL' : 0 < L (M + 1) (insertAt i r c)) (hq : 0 < q r) (hPB : 0 < PB M)
    (hPD : 0 < PD (M + 1)) :
    tupleTarget w pr L M c * (PB M / (M + 1) * q r) *
        min 1 (L (M + 1) (insertAt i r c) / L M c * (w (M + 1) / w M) *
          (PD (M + 1) / PB M) * (pr r / q r)) =
      tupleTarget w pr L (M + 1) (insertAt i r c) * (PD (M + 1) / (M + 1)) *
        min 1 (L (M + 1) (insertAt i r c) / L M c * (w (M + 1) / w M) *
          (PD (M + 1) / PB M) * (pr r / q r))⁻¹ := by
  have hprod : 0 < ∏ j, pr (c j) := Finset.prod_pos fun j _ => hpr _
  have hA : 0 < tupleTarget w pr L M c := by unfold tupleTarget; positivity
  have hB : 0 < tupleTarget w pr L (M + 1) (insertAt i r c) := by
    unfold tupleTarget
    have : 0 < ∏ j, pr (insertAt i r c j) := Finset.prod_pos fun j _ => hpr _
    positivity
  have hfwd : 0 < PB M / (M + 1) * q r := by positivity
  have hrev : 0 < PD (M + 1) / (M + 1) := by positivity
  rw [← birth_green_ratio w pr q L PB PD c i r hw.ne' hprod.ne' hL.ne' hq.ne' hPB.ne', inv_div]
  set A := tupleTarget w pr L M c * (PB M / (M + 1) * q r) with hAdef
  set B := tupleTarget w pr L (M + 1) (insertAt i r c) * (PD (M + 1) / (M + 1)) with hBdef
  have h := mh_detailed_balance (S := Bool) (fun b => if b then B else A)
    (fun _ _ => (1 : ℝ)) false true (by simp only [Bool.false_eq_true, ↓reduceIte]; positivity)
    (by simp only [↓reduceIte]; positivity) one_pos one_pos
  simp only [mhRatio, Bool.false_eq_true, if_false, if_true, mul_one, one_mul] at h
  exact h

end BirthDeath

end ActiveSubspace
