import Enthymemata.BassSensitivity
import Mathlib.Analysis.Calculus.Deriv.Polynomial
import Mathlib.MeasureTheory.Integral.IntervalIntegral.IntegrationByParts
import Mathlib.MeasureTheory.Integral.Pi
import Mathlib.Topology.Algebra.Polynomial

/-!
# BASS — the integral statements (SG-3 continuous, SG-4(b))

The parts of `Protolemmata/BASS/sensitivity-generator-lemmas.md` that need Mathlib's integrals. They
sit apart from `BassSensitivity` so that module's imports stay free of measure theory.

* **SG-3.** `mixture_product_integral_pi`: under a joint mixture of product measures, a product-form
  term's integral is the mixture of the per-component products of one-variable integrals.
* **SG-4(b).** With `⟨p, q⟩ = ∫₀¹ p q` and `J(p, q) = ∫₀¹ x(1−x) p′ q′`:
  * `l2Form_derivative_jacobi`: integration by parts, whose boundary term vanishes;
  * `jacobiForm_shiftedLegendre`: `J(P_m, P_n) = n(n+1) ⟨P_m, P_n⟩`;
  * `shiftedLegendre_orthogonal`: orthogonality with no separate proof;
  * `variance_le_half_jacobiForm`: `Var f ≤ ½ J(f, f)` for a finite expansion, with `½` the
    eigenvalue
    gap, and `variance_eq_half_jacobiForm_of_degree_one` showing it sharp at degree 1.

These stop short of normalizing: `∫P_n² = 1/(2n+1)` is not needed for any of them, and the
coefficients here are on `P_n` rather than on the orthonormal `ψ_n`.
-/

open MeasureTheory Polynomial intervalIntegral Finset
open scoped ENNReal NNReal

namespace BASS

/-! ### SG-3 -/

/-- **SG-3, on measures.** Under `Σ_c w_c ⊗_i μ_{c,i}`, a product-form term's integral is
`Σ_c w_c Π_i ∫ a_i dμ_{c,i}`. -/
theorem mixture_product_integral_pi {ι κ : Type*} [Fintype ι] [Fintype κ]
    {X : ι → Type*} [∀ i, MeasurableSpace (X i)]
    (w : κ → ℝ≥0) (μ : κ → ∀ i, Measure (X i)) [∀ c i, IsProbabilityMeasure (μ c i)]
    (a : ∀ i, X i → ℝ) (ha : ∀ c i, Integrable (a i) (μ c i)) :
    ∫ x, ∏ i, a i (x i) ∂(∑ c, (w c : ℝ≥0∞) • Measure.pi (μ c)) =
      ∑ c, (w c : ℝ) * ∏ i, ∫ y, a i y ∂(μ c i) := by
  have hint : ∀ c, Integrable (fun x : ∀ i, X i => ∏ i, a i (x i))
      ((w c : ℝ≥0∞) • Measure.pi (μ c)) := fun c =>
    (Integrable.fintype_prod_dep (fun i => ha c i)).smul_measure ENNReal.coe_ne_top
  rw [integral_finsetSum_measure fun c _ => hint c]
  refine Finset.sum_congr rfl fun c _ => ?_
  rw [MeasureTheory.integral_smul_measure, integral_fintype_prod_eq_prod, ENNReal.coe_toReal,
    smul_eq_mul]

/-! ### SG-4(b) -/

/-- `P_n` as a real polynomial. -/
noncomputable def legendreReal (n : ℕ) : ℝ[X] := (shiftedLegendre n).map (Int.castRingHom ℝ)

/-- **SG-4(a)** over `ℝ`. -/
theorem legendreReal_jacobi (n : ℕ) :
    derivative ((X * (1 - X)) * derivative (legendreReal n)) =
      -C ((n : ℝ) * (n + 1)) * legendreReal n := by
  have h := congrArg (Polynomial.map (Int.castRingHom ℝ)) (shiftedLegendre_jacobi n)
  simpa [legendreReal, Polynomial.map_mul, Polynomial.derivative_map] using h

/-- The `L²` form on `[0, 1]`. -/
noncomputable def l2Form (p q : ℝ[X]) : ℝ := ∫ x in (0 : ℝ)..1, p.eval x * q.eval x

/-- The weighted Dirichlet form `J(p, q) = ∫₀¹ x(1−x) p′ q′`. -/
noncomputable def jacobiForm (p q : ℝ[X]) : ℝ :=
  ∫ x in (0 : ℝ)..1, x * (1 - x) * (derivative p).eval x * (derivative q).eval x

theorem l2Form_comm (p q : ℝ[X]) : l2Form p q = l2Form q p := by
  simp only [l2Form, mul_comm]

theorem jacobiForm_comm (p q : ℝ[X]) : jacobiForm p q = jacobiForm q p := by
  simp only [jacobiForm]
  congr 1
  funext x
  ring

private theorem intervalIntegrable_eval (p : ℝ[X]) :
    IntervalIntegrable (fun x => p.eval x) volume 0 1 :=
  p.continuous.intervalIntegrable _ _

/-- **SG-4(b), integration by parts.** `∫₀¹ p · [x(1−x) q′]′ = −J(p, q)`: the boundary term
vanishes,
since `x(1−x)` is 0 at both ends. -/
theorem l2Form_derivative_jacobi (p q : ℝ[X]) :
    l2Form p (derivative ((X * (1 - X)) * derivative q)) = -jacobiForm p q := by
  set v : ℝ[X] := (X * (1 - X)) * derivative q with hv
  have hparts := integral_deriv_mul_eq_sub (a := 0) (b := 1)
    (u := fun x => p.eval x) (u' := fun x => (derivative p).eval x)
    (v := fun x => v.eval x) (v' := fun x => (derivative v).eval x)
    (fun x _ => p.hasDerivAt x) (fun x _ => v.hasDerivAt x)
    (intervalIntegrable_eval _) (intervalIntegrable_eval _)
  have hv0 : v.eval 0 = 0 := by simp [hv]
  have hv1 : v.eval 1 = 0 := by simp [hv]
  rw [hv0, hv1, mul_zero, mul_zero, sub_self] at hparts
  have hsplit := intervalIntegral.integral_add (μ := volume) (a := 0) (b := 1)
    (f := fun x => (derivative p).eval x * v.eval x) (g := fun x => p.eval x * (derivative v).eval
      x)
    (((derivative p).continuous.mul v.continuous).intervalIntegrable _ _)
    ((p.continuous.mul (derivative v).continuous).intervalIntegrable _ _)
  rw [hsplit] at hparts
  have hJ : jacobiForm p q = ∫ x in (0 : ℝ)..1, (derivative p).eval x * v.eval x := by
    simp only [jacobiForm, hv, eval_mul, eval_sub, eval_one, eval_X]
    congr 1
    funext x
    ring
  simp only [l2Form]
  linarith

/-- **SG-4(b).** `J(P_m, P_n) = n(n+1) ⟨P_m, P_n⟩`. -/
theorem jacobiForm_shiftedLegendre (m n : ℕ) :
    jacobiForm (legendreReal m) (legendreReal n) =
      (n : ℝ) * (n + 1) * l2Form (legendreReal m) (legendreReal n) := by
  have h := l2Form_derivative_jacobi (legendreReal m) (legendreReal n)
  rw [legendreReal_jacobi] at h
  have hscale : l2Form (legendreReal m) (-C ((n : ℝ) * (n + 1)) * legendreReal n) =
      -((n : ℝ) * (n + 1)) * l2Form (legendreReal m) (legendreReal n) := by
    simp only [l2Form, eval_mul, eval_neg, eval_C]
    rw [← intervalIntegral.integral_const_mul]
    congr 1
    funext x
    ring
  rw [hscale] at h
  linarith

/-- **SG-4(b).** Orthogonality, with no separate proof: `J` and `⟨·,·⟩` are symmetric and the
eigenvalues differ. -/
theorem shiftedLegendre_orthogonal {m n : ℕ} (hmn : m ≠ n) :
    l2Form (legendreReal m) (legendreReal n) = 0 := by
  have h₁ := jacobiForm_shiftedLegendre m n
  have h₂ := jacobiForm_shiftedLegendre n m
  rw [jacobiForm_comm, l2Form_comm] at h₂
  have hsub : ((n : ℝ) * (n + 1) - (m : ℝ) * (m + 1)) *
      l2Form (legendreReal m) (legendreReal n) = 0 := by
    linarith
  have hne : ((n : ℝ) * (n + 1) - (m : ℝ) * (m + 1)) ≠ 0 := by
    intro h
    have h' : (n : ℝ) * (n + 1) = (m : ℝ) * (m + 1) := by linarith
    have hnat : n * (n + 1) = m * (m + 1) := by exact_mod_cast h'
    apply hmn
    rcases lt_trichotomy m n with hlt | heq | hgt
    · nlinarith
    · exact heq
    · nlinarith
  exact (mul_eq_zero.1 hsub).resolve_left hne

/-- `P_0 = 1`. -/
theorem legendreReal_zero : legendreReal 0 = 1 := by
  ext k
  rcases k with _ | k <;> simp [legendreReal, coeff_shiftedLegendre, Polynomial.coeff_one]

section Expansion

variable (a : ℕ → ℝ) (d : ℕ)

/-- A finite expansion `f = Σ_{n ≤ d} a_n P_n`. -/
noncomputable def legendreExpansion : ℝ[X] := ∑ n ∈ range (d + 1), C (a n) * legendreReal n

private theorem l2Form_sum_left (s : Finset ℕ) (c : ℕ → ℝ) (g : ℕ → ℝ[X]) (q : ℝ[X]) :
    l2Form (∑ n ∈ s, C (c n) * g n) q = ∑ n ∈ s, c n * l2Form (g n) q := by
  simp only [l2Form, eval_finsetSum, eval_mul, eval_C, Finset.sum_mul]
  rw [intervalIntegral.integral_finsetSum]
  · refine Finset.sum_congr rfl fun n _ => ?_
    rw [← intervalIntegral.integral_const_mul]
    congr 1
    funext x
    ring
  · intro n _
    exact ((continuous_const.mul (g n).continuous).mul q.continuous).intervalIntegrable _ _

private theorem jacobiForm_sum_left (s : Finset ℕ) (c : ℕ → ℝ) (g : ℕ → ℝ[X]) (q : ℝ[X]) :
    jacobiForm (∑ n ∈ s, C (c n) * g n) q = ∑ n ∈ s, c n * jacobiForm (g n) q := by
  simp only [jacobiForm, derivative_sum, derivative_C_mul, eval_finsetSum, eval_mul, eval_C,
    Finset.mul_sum, Finset.sum_mul]
  rw [intervalIntegral.integral_finsetSum]
  · refine Finset.sum_congr rfl fun n _ => ?_
    rw [← intervalIntegral.integral_const_mul]
    congr 1
    funext x
    ring
  · intro n _
    exact (((continuous_id.mul (continuous_const.sub continuous_id)).mul
      (continuous_const.mul (derivative (g n)).continuous)).mul
      (derivative q).continuous).intervalIntegrable _ _

/-- Parseval for the expansion. -/
theorem l2Form_expansion_self :
    l2Form (legendreExpansion a d) (legendreExpansion a d) =
      ∑ n ∈ range (d + 1), a n ^ 2 * l2Form (legendreReal n) (legendreReal n) := by
  rw [legendreExpansion, l2Form_sum_left]
  refine Finset.sum_congr rfl fun n hn => ?_
  rw [l2Form_comm, l2Form_sum_left, Finset.sum_eq_single n]
  · rw [l2Form_comm]
    ring
  · intro k _ hk
    rw [shiftedLegendre_orthogonal hk, mul_zero]
  · intro h
    exact absurd hn h

/-- `J(f, f) = Σ_n a_n² n(n+1) ⟨P_n, P_n⟩`. -/
theorem jacobiForm_expansion_self :
    jacobiForm (legendreExpansion a d) (legendreExpansion a d) =
      ∑ n ∈ range (d + 1), a n ^ 2 * ((n : ℝ) * (n + 1) *
        l2Form (legendreReal n) (legendreReal n)) := by
  rw [legendreExpansion, jacobiForm_sum_left]
  refine Finset.sum_congr rfl fun n hn => ?_
  rw [jacobiForm_comm, jacobiForm_sum_left, Finset.sum_eq_single n]
  · rw [jacobiForm_shiftedLegendre]
    ring
  · intro k _ hk
    rw [jacobiForm_shiftedLegendre, shiftedLegendre_orthogonal hk, mul_zero, mul_zero]
  · intro h
    exact absurd hn h

/-- `∫₀¹ f` reads the expansion's constant coefficient. -/
theorem integral_expansion :
    l2Form 1 (legendreExpansion a d) = a 0 := by
  have h1 : l2Form 1 1 = 1 := by simp [l2Form]
  rw [l2Form_comm, legendreExpansion, l2Form_sum_left,
    Finset.sum_eq_single 0]
  · rw [legendreReal_zero, l2Form_comm, h1, mul_one]
  · intro k _ hk
    rw [← legendreReal_zero, shiftedLegendre_orthogonal hk, mul_zero]
  · simp

theorem l2Form_self_nonneg (p : ℝ[X]) : 0 ≤ l2Form p p := by
  refine intervalIntegral.integral_nonneg zero_le_one fun x _ => mul_self_nonneg _

/-- **SG-4(b), the Poincaré form.** `Var f ≤ ½ J(f, f)`: the Jacobi generator's gap is 2. This is
the
statement `total(v) ≤ ½ E[x_v(1−x_v)(∂_v f)²]` on one variable's fibre. -/
theorem variance_le_half_jacobiForm :
    l2Form (legendreExpansion a d) (legendreExpansion a d) - a 0 ^ 2 ≤
      (1 / 2) * jacobiForm (legendreExpansion a d) (legendreExpansion a d) := by
  have h1 : l2Form (legendreReal 0) (legendreReal 0) = 1 := by
    rw [legendreReal_zero]
    simp [l2Form]
  rw [l2Form_expansion_self, jacobiForm_expansion_self, Finset.mul_sum,
    Finset.sum_range_succ', Finset.sum_range_succ', h1]
  simp only [Nat.cast_zero, zero_add, mul_zero, pow_two, mul_one]
  have hterm : ∀ n ∈ range d, a (n + 1) * a (n + 1) *
      l2Form (legendreReal (n + 1)) (legendreReal (n + 1)) ≤
      1 / 2 * (a (n + 1) * a (n + 1) * (((n + 1 : ℕ) : ℝ) * (((n + 1 : ℕ) : ℝ) + 1) *
        l2Form (legendreReal (n + 1)) (legendreReal (n + 1)))) := by
    intro n _
    have hnn := l2Form_self_nonneg (legendreReal (n + 1))
    have hsq : 0 ≤ a (n + 1) * a (n + 1) := mul_self_nonneg _
    have hn : (2 : ℝ) ≤ ((n + 1 : ℕ) : ℝ) * (((n + 1 : ℕ) : ℝ) + 1) := by
      push_cast
      have : (0 : ℝ) ≤ n := Nat.cast_nonneg n
      nlinarith
    have : 0 ≤ a (n + 1) * a (n + 1) * l2Form (legendreReal (n + 1)) (legendreReal (n + 1)) :=
      mul_nonneg hsq hnn
    nlinarith
  have := Finset.sum_le_sum hterm
  linarith

/-- **SG-4(b), sharpness.** At degree 1 the bound is an equality. -/
theorem variance_eq_half_jacobiForm_of_degree_one (a₀ a₁ : ℝ) :
    let g : ℕ → ℝ := fun n => if n = 0 then a₀ else a₁
    l2Form (legendreExpansion g 1) (legendreExpansion g 1) - g 0 ^ 2 =
      (1 / 2) * jacobiForm (legendreExpansion g 1) (legendreExpansion g 1) := by
  intro g
  have h1 : l2Form (legendreReal 0) (legendreReal 0) = 1 := by
    rw [legendreReal_zero]
    simp [l2Form]
  rw [l2Form_expansion_self, jacobiForm_expansion_self]
  simp only [Finset.sum_range_succ, Finset.sum_range_zero, zero_add, h1]
  push_cast
  ring

end Expansion

end BASS
