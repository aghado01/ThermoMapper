import Enthymemata.BassSensitivity
import Mathlib.Analysis.Calculus.Deriv.Polynomial
import Mathlib.MeasureTheory.Integral.IntervalIntegral.IntegrationByParts
import Mathlib.MeasureTheory.Integral.Pi
import Mathlib.Topology.Algebra.Polynomial

/-!
# BASS — the integral statements (SG-3 continuous, SG-4(b) with its tensor form)

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

* **SG-4(b) on a tensor basis**, under the product uniform law on `[0, 1]^d`, for
  `f = Σ_{α ∈ S} a_α Π_v P_{α_v}` (`tensorExpansion`, with `tensorPartial` its derivative in `x_v`
  by `hasDerivAt_tensorExpansion`):
  * `integral_jacobi_tensorPartial`: `E[x_v(1−x_v)(∂_v f)²] = Σ_α a_α² α_v(α_v+1) N_α`;
  * `sum_integral_jacobi_tensorPartial`: summed over `v`, `Σ_α a_α² (Σ_v α_v(α_v+1)) N_α`;
  * `tensorTotal_eq`: the total numerator `E[(f − E_v f)²]` is `Σ_{α_v ≥ 1} a_α² N_α`;
  * `tensorTotal_le` and `tensorTotal_eq_half`: `total(v) ≤ ½·E[x_v(1−x_v)(∂_v f)²]`, an equality
    when every term has degree at most 1 in `x_v`.

These stop short of normalizing: `∫P_n² = 1/(2n+1)` is not needed for any of them, and the
coefficients here are on `P_n` rather than on the orthonormal `ψ_n`. On the tensor basis,
`a_α² N_α`, with `N_α = Π_v ⟨P_{α_v}, P_{α_v}⟩` (`tensorNorm`), is the squared orthonormal
coefficient, so the identities read as the protolemma states them.
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

/-! ### SG-4(b) on a tensor basis -/

section Tensor

variable {d : ℕ}

/-- The uniform law on `[0, 1]`, as the interval integrals read it. -/
noncomputable abbrev unitUniform : Measure ℝ := volume.restrict (Set.Ioc (0 : ℝ) 1)

/-- The product uniform law on `[0, 1]^d`. -/
noncomputable abbrev cubeUniform (d : ℕ) : Measure (Fin d → ℝ) :=
  Measure.pi fun _ : Fin d => unitUniform

theorem integral_unitUniform (g : ℝ → ℝ) : ∫ t, g t ∂unitUniform = ∫ t in (0 : ℝ)..1, g t :=
  (intervalIntegral.integral_of_le zero_le_one).symm

theorem integrable_unitUniform (p : ℝ[X]) : Integrable (fun t => p.eval t) unitUniform :=
  (p.continuous.integrableOn_Icc (a := 0) (b := 1)).mono_set Set.Ioc_subset_Icc_self

/-- Fubini on a product of polynomial factors: the cube's integral of `Π_w q_w(x_w)` is
`Π_w ∫₀¹ q_w`. -/
theorem integral_cube_prod (q : Fin d → ℝ[X]) :
    ∫ x, ∏ w, (q w).eval (x w) ∂cubeUniform d = ∏ w, ∫ t in (0 : ℝ)..1, (q w).eval t := by
  rw [integral_fintype_prod_eq_prod (f := fun w t => (q w).eval t)]
  exact Finset.prod_congr rfl fun w _ => integral_unitUniform _

/-- The cube's integral of a double sum of products of polynomial factors, term by term. -/
theorem integral_cube_sum_sum (T : Finset (Fin d → ℕ)) (c : (Fin d → ℕ) → (Fin d → ℕ) → ℝ)
    (q : (Fin d → ℕ) → (Fin d → ℕ) → Fin d → ℝ[X]) :
    ∫ x, ∑ α ∈ T, ∑ β ∈ T, c α β * ∏ w, (q α β w).eval (x w) ∂cubeUniform d =
      ∑ α ∈ T, ∑ β ∈ T, c α β * ∏ w, ∫ t in (0 : ℝ)..1, (q α β w).eval t := by
  have hint : ∀ α β, Integrable (fun x : Fin d → ℝ => c α β * ∏ w, (q α β w).eval (x w))
      (cubeUniform d) := fun α β =>
    (Integrable.fintype_prod (f := fun w t => (q α β w).eval t)
      fun w => integrable_unitUniform _).const_mul _
  rw [integral_finsetSum _ fun α _ => integrable_finsetSum _ fun β _ => hint α β]
  refine Finset.sum_congr rfl fun α _ => ?_
  rw [integral_finsetSum _ fun β _ => hint α β]
  refine Finset.sum_congr rfl fun β _ => ?_
  rw [MeasureTheory.integral_const_mul, integral_cube_prod]

/-- The tensor basis function `x ↦ Π_v P_{α_v}(x_v)`. -/
noncomputable def tensorEval (α : Fin d → ℕ) (x : Fin d → ℝ) : ℝ :=
  ∏ v, (legendreReal (α v)).eval (x v)

/-- A finite tensor expansion `f = Σ_{α ∈ S} a_α Π_v P_{α_v}`. -/
noncomputable def tensorExpansion (S : Finset (Fin d → ℕ)) (a : (Fin d → ℕ) → ℝ)
    (x : Fin d → ℝ) : ℝ :=
  ∑ α ∈ S, a α * tensorEval α x

/-- Its partial derivative in `x_v`: each term's `v`th factor differentiated
(`hasDerivAt_tensorExpansion`). -/
noncomputable def tensorPartial (S : Finset (Fin d → ℕ)) (a : (Fin d → ℕ) → ℝ) (v : Fin d)
    (x : Fin d → ℝ) : ℝ :=
  ∑ α ∈ S, a α * ((derivative (legendreReal (α v))).eval (x v) *
    ∏ w ∈ univ.erase v, (legendreReal (α w)).eval (x w))

/-- A basis function's squared norm, `N_α = Π_v ⟨P_{α_v}, P_{α_v}⟩`; `a_α² N_α` is the squared
coefficient on the orthonormal basis. -/
noncomputable def tensorNorm (α : Fin d → ℕ) : ℝ :=
  ∏ v, l2Form (legendreReal (α v)) (legendreReal (α v))

theorem tensorNorm_nonneg (α : Fin d → ℕ) : 0 ≤ tensorNorm α :=
  Finset.prod_nonneg fun _ _ => l2Form_self_nonneg _

theorem prod_update_eq (g : Fin d → ℝ → ℝ) (x : Fin d → ℝ) (v : Fin d) (t : ℝ) :
    ∏ w, g w (Function.update x v t w) = g v t * ∏ w ∈ univ.erase v, g w (x w) := by
  rw [← mul_prod_erase univ _ (mem_univ v), Function.update_self]
  congr 1
  exact Finset.prod_congr rfl fun w hw => by rw [Function.update_of_ne (ne_of_mem_erase hw)]

/-- `tensorPartial` is the partial derivative along `x_v`. -/
theorem hasDerivAt_tensorExpansion (S : Finset (Fin d → ℕ)) (a : (Fin d → ℕ) → ℝ) (v : Fin d)
    (x : Fin d → ℝ) :
    HasDerivAt (fun t => tensorExpansion S a (Function.update x v t)) (tensorPartial S a v x)
      (x v) := by
  have hfun : (fun t => tensorExpansion S a (Function.update x v t)) = fun t =>
      ∑ α ∈ S, a α * ((legendreReal (α v)).eval t *
        ∏ w ∈ univ.erase v, (legendreReal (α w)).eval (x w)) := by
    funext t
    simp only [tensorExpansion, tensorEval]
    refine Finset.sum_congr rfl fun α _ => ?_
    rw [prod_update_eq (fun w s => (legendreReal (α w)).eval s)]
  rw [hfun]
  exact HasDerivAt.fun_sum (u := S)
    (A := fun α t => a α * ((legendreReal (α v)).eval t *
      ∏ w ∈ univ.erase v, (legendreReal (α w)).eval (x w)))
    (A' := fun α => a α * ((derivative (legendreReal (α v))).eval (x v) *
      ∏ w ∈ univ.erase v, (legendreReal (α w)).eval (x w)))
    fun α _ => (((legendreReal (α v)).hasDerivAt (x v)).mul_const _).const_mul _

/-- **SG-4(b), per variable.** `E[x_v(1−x_v)(∂_v f)²] = Σ_α a_α² α_v(α_v+1) N_α` under the
product uniform law: only the diagonal survives, by `J` and `⟨·,·⟩` orthogonality. -/
theorem integral_jacobi_tensorPartial (S : Finset (Fin d → ℕ)) (a : (Fin d → ℕ) → ℝ)
    (v : Fin d) :
    ∫ x, x v * (1 - x v) * tensorPartial S a v x ^ 2 ∂cubeUniform d =
      ∑ α ∈ S, a α ^ 2 * ((α v : ℝ) * (α v + 1) * tensorNorm α) := by
  classical
  set q : (Fin d → ℕ) → (Fin d → ℕ) → Fin d → ℝ[X] := fun (α β : Fin d → ℕ) (w : Fin d) =>
    if w = v then X * (1 - X) * derivative (legendreReal (α v)) * derivative (legendreReal (β v))
    else legendreReal (α w) * legendreReal (β w) with hq
  have hpt : ∀ x : Fin d → ℝ, x v * (1 - x v) * tensorPartial S a v x ^ 2 =
      ∑ α ∈ S, ∑ β ∈ S, a α * a β * ∏ w, (q α β w).eval (x w) := by
    intro x
    rw [sq, tensorPartial, Finset.sum_mul_sum, Finset.mul_sum]
    refine Finset.sum_congr rfl fun α _ => ?_
    rw [Finset.mul_sum]
    refine Finset.sum_congr rfl fun β _ => ?_
    rw [← mul_prod_erase univ _ (mem_univ v)]
    have hv : (q α β v).eval (x v) = x v * (1 - x v) *
        (derivative (legendreReal (α v))).eval (x v) *
        (derivative (legendreReal (β v))).eval (x v) := by
      simp [hq]
    have hw : ∏ w ∈ univ.erase v, (q α β w).eval (x w) =
        (∏ w ∈ univ.erase v, (legendreReal (α w)).eval (x w)) *
          ∏ w ∈ univ.erase v, (legendreReal (β w)).eval (x w) := by
      rw [← Finset.prod_mul_distrib]
      exact Finset.prod_congr rfl fun w hw => by simp [hq, ne_of_mem_erase hw]
    rw [hv, hw]
    ring
  simp_rw [hpt]
  rw [integral_cube_sum_sum (c := fun α β => a α * a β)]
  have hone : ∀ α β, ∏ w, ∫ t in (0 : ℝ)..1, (q α β w).eval t =
      jacobiForm (legendreReal (α v)) (legendreReal (β v)) *
        ∏ w ∈ univ.erase v, l2Form (legendreReal (α w)) (legendreReal (β w)) := by
    intro α β
    rw [← mul_prod_erase univ _ (mem_univ v)]
    congr 1
    · simp only [hq, if_pos rfl, jacobiForm, eval_mul, eval_sub, eval_one, eval_X]
    · exact Finset.prod_congr rfl fun w hw => by
        simp [hq, ne_of_mem_erase hw, l2Form, eval_mul]
  simp only [hone]
  refine Finset.sum_congr rfl fun α hα => ?_
  rw [Finset.sum_eq_single α]
  · rw [jacobiForm_shiftedLegendre, tensorNorm,
      ← mul_prod_erase univ (fun w => l2Form (legendreReal (α w)) (legendreReal (α w)))
        (mem_univ v)]
    ring
  · intro β _ hβ
    obtain ⟨w, hw⟩ : ∃ w, α w ≠ β w := Function.ne_iff.1 (Ne.symm hβ)
    by_cases hwv : w = v
    · rw [hwv] at hw
      rw [jacobiForm_shiftedLegendre, shiftedLegendre_orthogonal hw]
      ring
    · rw [Finset.prod_eq_zero (i := w) (mem_erase.2 ⟨hwv, mem_univ w⟩)
        (shiftedLegendre_orthogonal hw)]
      ring
  · intro h
    exact absurd hα h

/-- **SG-4(b), summed over the variables.** `Σ_v E[x_v(1−x_v)(∂_v f)²] =
Σ_α a_α² (Σ_v α_v(α_v+1)) N_α`: on the orthonormal basis, `Σ_α b_α² Σ_v α_v(α_v+1)`, each
coefficient weighted by the tensor generator's eigenvalue. -/
theorem sum_integral_jacobi_tensorPartial (S : Finset (Fin d → ℕ)) (a : (Fin d → ℕ) → ℝ) :
    ∑ v, ∫ x, x v * (1 - x v) * tensorPartial S a v x ^ 2 ∂cubeUniform d =
      ∑ α ∈ S, a α ^ 2 * (∑ v, (α v : ℝ) * (α v + 1)) * tensorNorm α := by
  simp only [integral_jacobi_tensorPartial]
  rw [Finset.sum_comm]
  refine Finset.sum_congr rfl fun α _ => ?_
  rw [Finset.mul_sum, Finset.sum_mul]
  exact Finset.sum_congr rfl fun v _ => by ring

/-- Parseval on the tensor basis: `E[(Σ_α a_α Ψ_α)²] = Σ_α a_α² N_α`. -/
theorem integral_tensor_sq (T : Finset (Fin d → ℕ)) (a : (Fin d → ℕ) → ℝ) :
    ∫ x, (∑ α ∈ T, a α * tensorEval α x) ^ 2 ∂cubeUniform d = ∑ α ∈ T, a α ^ 2 * tensorNorm α := by
  have hpt : ∀ x : Fin d → ℝ, (∑ α ∈ T, a α * tensorEval α x) ^ 2 = ∑ α ∈ T, ∑ β ∈ T,
      a α * a β * ∏ w, (legendreReal (α w) * legendreReal (β w)).eval (x w) := by
    intro x
    rw [sq, Finset.sum_mul_sum]
    refine Finset.sum_congr rfl fun α _ => Finset.sum_congr rfl fun β _ => ?_
    simp only [tensorEval, eval_mul, Finset.prod_mul_distrib]
    ring
  simp_rw [hpt]
  rw [integral_cube_sum_sum (c := fun α β => a α * a β)
    (q := fun α β w => legendreReal (α w) * legendreReal (β w))]
  refine Finset.sum_congr rfl fun α hα => ?_
  rw [Finset.sum_eq_single α]
  · simp only [tensorNorm, l2Form, eval_mul]
    ring
  · intro β _ hβ
    obtain ⟨w, hw⟩ : ∃ w, α w ≠ β w := Function.ne_iff.1 (Ne.symm hβ)
    have h0 : ∫ t in (0 : ℝ)..1, (legendreReal (α w) * legendreReal (β w)).eval t = 0 := by
      simpa [l2Form, eval_mul] using shiftedLegendre_orthogonal hw
    rw [Finset.prod_eq_zero (mem_univ w) h0]
    ring
  · intro h
    exact absurd hα h

theorem integral_legendreReal (n : ℕ) :
    ∫ t in (0 : ℝ)..1, (legendreReal n).eval t = if n = 0 then 1 else 0 := by
  have h : ∫ t in (0 : ℝ)..1, (legendreReal n).eval t =
      l2Form (legendreReal 0) (legendreReal n) := by
    simp [l2Form, legendreReal_zero]
  rw [h]
  split_ifs with hn
  · subst hn
    rw [legendreReal_zero]
    simp [l2Form]
  · exact shiftedLegendre_orthogonal (Ne.symm hn)

/-- Integrating out `x_v`: `E_v f(x) = ∫₀¹ f(x with x_v = t) dt`, the mean given the other
variables. -/
noncomputable def tensorCondMean (S : Finset (Fin d → ℕ)) (a : (Fin d → ℕ) → ℝ) (v : Fin d)
    (x : Fin d → ℝ) : ℝ :=
  ∫ t in (0 : ℝ)..1, tensorExpansion S a (Function.update x v t)

/-- `E_v f` keeps exactly the terms constant in `x_v`. -/
theorem tensorCondMean_eq (S : Finset (Fin d → ℕ)) (a : (Fin d → ℕ) → ℝ) (v : Fin d)
    (x : Fin d → ℝ) :
    tensorCondMean S a v x = ∑ α ∈ S.filter (fun α => α v = 0), a α * tensorEval α x := by
  have hpt : ∀ t, tensorExpansion S a (Function.update x v t) = ∑ α ∈ S,
      (a α * ∏ w ∈ univ.erase v, (legendreReal (α w)).eval (x w)) *
        (legendreReal (α v)).eval t := by
    intro t
    simp only [tensorExpansion, tensorEval]
    refine Finset.sum_congr rfl fun α _ => ?_
    rw [prod_update_eq (fun w s => (legendreReal (α w)).eval s)]
    ring
  simp only [tensorCondMean, hpt]
  rw [intervalIntegral.integral_finsetSum fun α _ =>
    (continuous_const.mul (legendreReal (α v)).continuous).intervalIntegrable 0 1]
  simp only [intervalIntegral.integral_const_mul, integral_legendreReal]
  rw [Finset.sum_filter]
  refine Finset.sum_congr rfl fun α _ => ?_
  split_ifs with h
  · simp only [tensorEval]
    rw [← mul_prod_erase univ _ (mem_univ v), h, legendreReal_zero]
    simp
  · simp

/-- The Sobol total numerator of `x_v`, `E[(f − E_v f)²]`. -/
noncomputable def tensorTotal (S : Finset (Fin d → ℕ)) (a : (Fin d → ℕ) → ℝ) (v : Fin d) : ℝ :=
  ∫ x, (tensorExpansion S a x - tensorCondMean S a v x) ^ 2 ∂cubeUniform d

/-- The total numerator read from the coefficients: `Σ_{α_v ≥ 1} a_α² N_α`. -/
theorem tensorTotal_eq (S : Finset (Fin d → ℕ)) (a : (Fin d → ℕ) → ℝ) (v : Fin d) :
    tensorTotal S a v = ∑ α ∈ S.filter (fun α => α v ≠ 0), a α ^ 2 * tensorNorm α := by
  have hpt : ∀ x, tensorExpansion S a x - tensorCondMean S a v x =
      ∑ α ∈ S.filter (fun α => α v ≠ 0), a α * tensorEval α x := by
    intro x
    rw [tensorCondMean_eq, tensorExpansion,
      ← Finset.sum_filter_add_sum_filter_not S (fun α => α v = 0)]
    ring
  simp only [tensorTotal, hpt]
  exact integral_tensor_sq _ a

/-- **SG-4(b), the per-variable bound.** `total(v) ≤ ½·E[x_v(1−x_v)(∂_v f)²]`: the weighted
Poincaré inequality with the Jacobi generator's gap, `α_v(α_v+1) ≥ 2` on every term with
`α_v ≥ 1`. -/
theorem tensorTotal_le (S : Finset (Fin d → ℕ)) (a : (Fin d → ℕ) → ℝ) (v : Fin d) :
    tensorTotal S a v ≤
      (1 / 2) * ∫ x, x v * (1 - x v) * tensorPartial S a v x ^ 2 ∂cubeUniform d := by
  rw [tensorTotal_eq, integral_jacobi_tensorPartial, Finset.mul_sum, Finset.sum_filter]
  refine Finset.sum_le_sum fun α _ => ?_
  have hN := tensorNorm_nonneg α
  have hsq : 0 ≤ a α ^ 2 * tensorNorm α := mul_nonneg (sq_nonneg _) hN
  split_ifs with h
  · have h1 : (1 : ℝ) ≤ α v := by exact_mod_cast Nat.one_le_iff_ne_zero.2 h
    have hk : (2 : ℝ) ≤ (α v : ℝ) * (α v + 1) := by nlinarith
    nlinarith [mul_le_mul_of_nonneg_right hk hsq]
  · have h0 : (0 : ℝ) ≤ (α v : ℝ) * (α v + 1) := by positivity
    exact mul_nonneg (by norm_num) (mul_nonneg (sq_nonneg _) (mul_nonneg h0 hN))

/-- **SG-4(b), sharpness.** With every term of degree at most 1 in `x_v`, the bound is an
equality. -/
theorem tensorTotal_eq_half (S : Finset (Fin d → ℕ)) (a : (Fin d → ℕ) → ℝ) (v : Fin d)
    (h1 : ∀ α ∈ S, α v ≤ 1) :
    tensorTotal S a v =
      (1 / 2) * ∫ x, x v * (1 - x v) * tensorPartial S a v x ^ 2 ∂cubeUniform d := by
  rw [tensorTotal_eq, integral_jacobi_tensorPartial, Finset.mul_sum, Finset.sum_filter]
  refine Finset.sum_congr rfl fun α hα => ?_
  split_ifs with h
  · have : α v = 1 := by have := h1 α hα; omega
    rw [this]
    push_cast
    ring
  · have : α v = 0 := not_not.1 h
    rw [this]
    push_cast
    ring

end Tensor

end BASS
