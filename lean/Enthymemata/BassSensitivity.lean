import Mathlib.Algebra.BigOperators.Field
import Mathlib.Algebra.Polynomial.Derivative
import Mathlib.Analysis.Real.Sqrt
import Mathlib.Data.Fintype.BigOperators
import Mathlib.RingTheory.Adjoin.Basic
import Mathlib.RingTheory.Polynomial.ShiftedLegendre
import Mathlib.Tactic.FieldSimp
import Mathlib.Tactic.LinearCombination
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.Ring

/-!
# BASS — sensitivity as a Markov generator (SG-1, SG-2(a), SG-3, SG-4, SG-5)

The self-contained statements of `Protolemmata/BASS/sensitivity-generator-lemmas.md`.

**SG-1, abstractly.** `CommutingMeans` is a commuting family of idempotents, self-adjoint for a
nonnegative symmetric form. It takes a finite product and a continuous one alike; Fubini supplies
the
family in the second case. Its components `P_u = Π_{v∈u}(I − E_v)·Π_{v∉u}E_v` live in the
commutative
algebra the family generates (`hoeffdingOp` is the same product in any commutative ring).
* `sum_component`, `component_component_of_ne`, `component_idem`: a complete orthogonal family.
* `generator_eigen`, and `generator_component`: `L P_u = −|u| P_u`.
* `parseval`, `total_eq_sum`, `sum_total`, `dirichlet_eq_sum_total`: `Σ_v total(v) = −⟨f, Lf⟩ =
  Σ_u |u| V_u`.
* `efron_stein_form`, `efron_stein_form_eq_iff`: the spectral gap of 1, with equality exactly when
  `V_u = 0` for `|u| ≥ 2`; `efron_stein_eq_iff_additive` reads that as additivity on a finite
    product.
* `form_noise`: `⟨f, U_ρ f⟩ = Σ_u ρ^|u| V_u`, the polynomial whose coefficients are the variance by
  interaction order.

**SG-1, on a finite product.** `coordMean_idem`, `coordMean_comm`, `coordMean_selfAdjoint`;
`efron_stein`; and `noise_prod_form`, the factorization of `U_ρ` on product-form functions.

**SG-2(a).** `sqrt_compVar_sub_le`, `sqrt_compVar_sub_le_form`: `√V_u` is 1-Lipschitz.

**SG-3.** `mixture_product_integral`, the finite form.

**SG-4.** `shiftedLegendre_jacobi`, from the coefficients.

**SG-5.** `no_common_eigenvector`, `top_eigenspace_ne_subspace_mean`: with explicit 2 × 2 matrices,
`E[C]` and the mean of the draws' projectors share no eigenvector, so their Rayleigh maximizers
differ. `rayleigh_max_is_eigen` is the Rayleigh step.
-/

open Finset

namespace BASS

/-! ### SG-1: the components in a commutative ring -/

section Ring

variable {ι R : Type*} [Fintype ι] [DecidableEq ι] [CommRing R]

/-- `P_u = Π_{v∈u}(1 − e_v) · Π_{v∉u} e_v`. -/
def hoeffdingOp (e : ι → R) (u : Finset ι) : R := (∏ v ∈ u, (1 - e v)) * ∏ v ∈ uᶜ, e v

/-- **SG-1(b), completeness.** `Σ_u P_u = 1`. This needs nothing of `e`. -/
theorem sum_hoeffdingOp (e : ι → R) : ∑ u, hoeffdingOp e u = 1 := by
  have h := Finset.prod_add (fun v => 1 - e v) e (univ : Finset ι)
  simp only [sub_add_cancel, Finset.prod_const_one, Finset.powerset_univ] at h
  rw [h]
  exact Finset.sum_congr rfl fun u _ => by rw [hoeffdingOp, Finset.compl_eq_univ_sdiff]

/-- **SG-1(e), the operator form.** `Π_v (ρ + (1 − ρ) e_v) = Σ_u ρ^|u| P_u`. -/
theorem noise_eq_sum (e : ι → R) (r : R) :
    ∏ v, (r + (1 - r) * e v) = ∑ u, r ^ u.card * hoeffdingOp e u := by
  have h := Finset.prod_add (fun v => r * (1 - e v)) e (univ : Finset ι)
  rw [Finset.powerset_univ] at h
  have hl : ∀ v, r + (1 - r) * e v = r * (1 - e v) + e v := fun v => by ring
  simp only [hl]
  rw [h]
  refine Finset.sum_congr rfl fun u _ => ?_
  rw [Finset.prod_mul_distrib, Finset.prod_const, hoeffdingOp, Finset.compl_eq_univ_sdiff]
  ring

variable {e : ι → R}

/-- Inside `u`, `e_v` kills `P_u`. -/
theorem e_mul_hoeffdingOp_of_mem (he : ∀ v, e v * e v = e v) {u : Finset ι} {v : ι} (hv : v ∈ u) :
    e v * hoeffdingOp e u = 0 := by
  rw [hoeffdingOp, ← Finset.mul_prod_erase u _ hv]
  have : e v * (1 - e v) = 0 := by rw [mul_sub, he, mul_one, sub_self]
  linear_combination ((∏ x ∈ u.erase v, (1 - e x)) * ∏ x ∈ uᶜ, e x) * this

/-- Outside `u`, `e_v` fixes `P_u`. -/
theorem e_mul_hoeffdingOp_of_notMem (he : ∀ v, e v * e v = e v) {u : Finset ι} {v : ι}
    (hv : v ∉ u) : e v * hoeffdingOp e u = hoeffdingOp e u := by
  have hvc : v ∈ uᶜ := Finset.mem_compl.2 hv
  rw [hoeffdingOp, ← Finset.mul_prod_erase uᶜ _ hvc]
  linear_combination ((∏ x ∈ u, (1 - e x)) * ∏ x ∈ uᶜ.erase v, e x) * he v

/-- **SG-1(b), orthogonality.** Distinct components multiply to 0. -/
theorem hoeffdingOp_mul_of_ne (he : ∀ v, e v * e v = e v) {u w : Finset ι} (huw : u ≠ w) :
    hoeffdingOp e u * hoeffdingOp e w = 0 := by
  have : ∃ v, (v ∈ u ∧ v ∉ w) ∨ (v ∉ u ∧ v ∈ w) := by
    by_contra hne
    apply huw
    ext v
    have := not_exists.1 hne v
    tauto
  obtain ⟨v, ⟨hu, hw⟩ | ⟨hu, hw⟩⟩ := this
  · have h₁ := e_mul_hoeffdingOp_of_mem he hu
    have h₂ := e_mul_hoeffdingOp_of_notMem he hw
    linear_combination hoeffdingOp e u * (-h₂) + hoeffdingOp e w * h₁
  · have h₁ := e_mul_hoeffdingOp_of_notMem he hu
    have h₂ := e_mul_hoeffdingOp_of_mem he hw
    linear_combination hoeffdingOp e u * h₂ + hoeffdingOp e w * (-h₁)

/-- **SG-1(b), idempotence.** -/
theorem hoeffdingOp_idem (he : ∀ v, e v * e v = e v) (u : Finset ι) :
    hoeffdingOp e u * hoeffdingOp e u = hoeffdingOp e u := by
  have h : hoeffdingOp e u * ∑ w, hoeffdingOp e w = hoeffdingOp e u := by
    rw [sum_hoeffdingOp, mul_one]
  rw [Finset.mul_sum, Finset.sum_eq_single u] at h
  · exact h
  · intro w _ hw
    exact hoeffdingOp_mul_of_ne he (Ne.symm hw)
  · simp

end Ring

/-! ### SG-1(c): the abstract step -/

/-- **SG-1(c), the abstract step.** If the `P i` kill `f` on `u` and fix it off `u`, then
`Σ_i (P_i f − f) = −|u| f`. -/
theorem generator_eigen {E ι : Type*} [AddCommGroup E] [Module ℝ E] [Fintype ι]
    (P : ι → E →ₗ[ℝ] E) (u : Finset ι) (f : E)
    (h₀ : ∀ i ∈ u, P i f = 0) (h₁ : ∀ i ∉ u, P i f = f) :
    (∑ i, (P i f - f)) = -(u.card : ℝ) • f := by
  classical
  have hsplit : ∑ i, (P i f - f) = ∑ i ∈ u, (P i f - f) + ∑ i ∈ uᶜ, (P i f - f) :=
    (Finset.sum_add_sum_compl u _).symm
  have hin : ∑ i ∈ u, (P i f - f) = ∑ _i ∈ u, -f :=
    Finset.sum_congr rfl fun i hi => by rw [h₀ i hi, zero_sub]
  have hout : ∑ i ∈ uᶜ, (P i f - f) = 0 :=
    Finset.sum_eq_zero fun i hi => by rw [h₁ i (Finset.mem_compl.1 hi), sub_self]
  rw [hsplit, hin, hout, add_zero, Finset.sum_const, neg_smul, smul_neg, Nat.cast_smul_eq_nsmul]

/-! ### SG-1 over a commuting family -/

/-- A commuting family of idempotents, self-adjoint for a nonnegative symmetric form: conditional
expectations on a product. -/
structure CommutingMeans (ι E : Type*) [AddCommGroup E] [Module ℝ E] where
  /-- `E_v`. -/
  mean : ι → E →ₗ[ℝ] E
  /-- Each `E_v` is idempotent. -/
  idem : ∀ i f, mean i (mean i f) = mean i f
  /-- The family commutes: where independence enters. -/
  comm : ∀ i j f, mean i (mean j f) = mean j (mean i f)
  /-- The form `⟨f, g⟩`. -/
  form : E →ₗ[ℝ] E →ₗ[ℝ] ℝ
  /-- The form is symmetric. -/
  form_symm : ∀ f g, form f g = form g f
  /-- Each `E_v` is self-adjoint. -/
  form_adj : ∀ i f g, form (mean i f) g = form f (mean i g)
  /-- The form is nonnegative. -/
  form_nonneg : ∀ f, 0 ≤ form f f

namespace CommutingMeans

variable {ι E : Type*} [Fintype ι] [DecidableEq ι] [AddCommGroup E] [Module ℝ E]
variable (F : CommutingMeans ι E)

/-- The commutative algebra the family generates. -/
abbrev algebraOf : Subalgebra ℝ (E →ₗ[ℝ] E) := Algebra.adjoin ℝ (Set.range F.mean)

omit [Fintype ι] [DecidableEq ι] in
theorem mean_mem (i : ι) : F.mean i ∈ F.algebraOf := Algebra.subset_adjoin ⟨i, rfl⟩

omit [Fintype ι] [DecidableEq ι] in
theorem algebraOf_comm : IsMulCommutative F.algebraOf :=
  Algebra.isMulCommutative_adjoin ℝ (by
    rintro _ ⟨i, rfl⟩ _ ⟨j, rfl⟩
    exact LinearMap.ext fun f => F.comm i j f)

/-- The generated algebra is commutative. -/
@[reducible] noncomputable def commRing : CommRing F.algebraOf :=
  { (inferInstance : Ring F.algebraOf) with
    mul_comm := (F.algebraOf_comm).is_comm.comm }

/-- `E_v` in the generated algebra. -/
def gen (i : ι) : F.algebraOf := ⟨F.mean i, F.mean_mem i⟩

omit [Fintype ι] [DecidableEq ι] in
theorem gen_mul_gen (i : ι) : F.gen i * F.gen i = F.gen i :=
  Subtype.ext (LinearMap.ext fun f => F.idem i f)

/-- **SG-1(b).** The component `P_u`, as an operator. -/
noncomputable def component (u : Finset ι) : E →ₗ[ℝ] E :=
  letI := F.commRing
  ((hoeffdingOp F.gen u : F.algebraOf) : E →ₗ[ℝ] E)

/-- `V_u = ‖P_u f‖²`. -/
noncomputable def compVar (u : Finset ι) (f : E) : ℝ := F.form (F.component u f) (F.component u f)

/-- `total(v) = E[Var(f | x_{−v})] = ⟨f, f − E_v f⟩`. -/
noncomputable def total (v : ι) (f : E) : ℝ := F.form f (f - F.mean v f)

/-- The generator `L = Σ_v (E_v − I)`, applied. -/
noncomputable def generator (f : E) : E := ∑ v, (F.mean v f - f)

/-- `Π_{v∈s}(ρI + (1 − ρ)E_v)`, as an operator. -/
noncomputable def noisePart (ρ : ℝ) (s : Finset ι) : E →ₗ[ℝ] E :=
  letI := F.commRing
  ((∏ v ∈ s, (algebraMap ℝ F.algebraOf ρ + (1 - algebraMap ℝ F.algebraOf ρ) * F.gen v) :
    F.algebraOf) : E →ₗ[ℝ] E)

/-- `U_ρ = Π_v(ρI + (1 − ρ)E_v)`, as an operator. -/
noncomputable def noise (ρ : ℝ) : E →ₗ[ℝ] E := F.noisePart ρ univ

omit [Fintype ι] [DecidableEq ι] in
theorem noisePart_empty (ρ : ℝ) (f : E) : F.noisePart ρ ∅ f = f := by
  letI := F.commRing
  simp [noisePart]

omit [Fintype ι] in
theorem noisePart_insert (ρ : ℝ) {a : ι} {s : Finset ι} (ha : a ∉ s) (f : E) :
    F.noisePart ρ (insert a s) f =
      ρ • F.noisePart ρ s f + (1 - ρ) • F.mean a (F.noisePart ρ s f) := by
  letI := F.commRing
  simp only [noisePart]
  rw [Finset.prod_insert ha]
  simp [gen, Module.End.mul_apply, Module.algebraMap_end_apply, sub_smul]

/-- **SG-1(b), completeness.** `Σ_u P_u f = f`. -/
theorem sum_component (f : E) : ∑ u, F.component u f = f := by
  letI := F.commRing
  have h := congrArg (fun a : F.algebraOf => ((a : E →ₗ[ℝ] E) f)) (sum_hoeffdingOp F.gen)
  simp only [AddSubmonoidClass.coe_finsetSum, LinearMap.sum_apply] at h
  simpa [component] using h

/-- Inside `u`, `E_v` kills `P_u f`. -/
theorem mean_component_of_mem {u : Finset ι} {v : ι} (hv : v ∈ u) (f : E) :
    F.mean v (F.component u f) = 0 := by
  letI := F.commRing
  have h := congrArg (fun a : F.algebraOf => ((a : E →ₗ[ℝ] E) f))
    (e_mul_hoeffdingOp_of_mem F.gen_mul_gen hv)
  simpa [component, gen, Module.End.mul_apply] using h

/-- Outside `u`, `E_v` fixes `P_u f`. -/
theorem mean_component_of_notMem {u : Finset ι} {v : ι} (hv : v ∉ u) (f : E) :
    F.mean v (F.component u f) = F.component u f := by
  letI := F.commRing
  have h := congrArg (fun a : F.algebraOf => ((a : E →ₗ[ℝ] E) f))
    (e_mul_hoeffdingOp_of_notMem F.gen_mul_gen hv)
  simpa [component, gen, Module.End.mul_apply] using h

/-- `P_u` commutes with `E_v`. -/
theorem component_mean (u : Finset ι) (v : ι) (f : E) :
    F.component u (F.mean v f) = F.mean v (F.component u f) := by
  letI := F.commRing
  have h := congrArg (fun a : F.algebraOf => ((a : E →ₗ[ℝ] E) f))
    (mul_comm (hoeffdingOp F.gen u) (F.gen v))
  simpa [component, gen, Module.End.mul_apply] using h

/-- **SG-1(b), orthogonality as operators.** -/
theorem component_component_of_ne {u w : Finset ι} (huw : u ≠ w) (f : E) :
    F.component u (F.component w f) = 0 := by
  letI := F.commRing
  have h := congrArg (fun a : F.algebraOf => ((a : E →ₗ[ℝ] E) f))
    (hoeffdingOp_mul_of_ne F.gen_mul_gen huw)
  simpa [component, Module.End.mul_apply] using h

/-- **SG-1(b), idempotence as operators.** -/
theorem component_idem (u : Finset ι) (f : E) :
    F.component u (F.component u f) = F.component u f := by
  letI := F.commRing
  have h := congrArg (fun a : F.algebraOf => ((a : E →ₗ[ℝ] E) f))
    (hoeffdingOp_idem F.gen_mul_gen u)
  simpa [component, Module.End.mul_apply] using h

omit [Fintype ι] [DecidableEq ι] in
/-- Every element of the generated algebra is self-adjoint: the product of two commuting
self-adjoint
operators is self-adjoint, which is where the family's commuting enters. -/
theorem form_algebra_adj (a : F.algebraOf) (f g : E) :
    F.form ((a : E →ₗ[ℝ] E) f) g = F.form f ((a : E →ₗ[ℝ] E) g) := by
  obtain ⟨a, ha⟩ := a
  induction ha using Algebra.adjoin_induction generalizing f g with
  | mem x hx =>
    obtain ⟨i, rfl⟩ := hx
    exact F.form_adj i f g
  | algebraMap r =>
    simp [Module.algebraMap_end_apply]
  | add x y _ _ hx hy =>
    simpa only [LinearMap.add_apply, map_add] using congrArg₂ (· + ·) (hx f g) (hy f g)
  | mul x y hxm hym hx hy =>
    have hcomm : x * y = y * x := by
      have h := F.algebraOf_comm.is_comm.comm (⟨x, hxm⟩ : F.algebraOf) ⟨y, hym⟩
      exact congrArg Subtype.val h
    have hswap : x (y g) = y (x g) := by
      have h := congrArg (fun a : E →ₗ[ℝ] E => a g) hcomm
      simpa only [Module.End.mul_apply] using h
    simp only [Module.End.mul_apply]
    change F.form (x (y f)) g = F.form f (x (y g))
    have h₁ := hx (y f) g
    have h₂ := hy f (x g)
    rw [hswap]
    exact h₁.trans h₂

/-- `P_u` is self-adjoint. -/
theorem form_component_adj (u : Finset ι) (f g : E) :
    F.form (F.component u f) g = F.form f (F.component u g) := by
  letI := F.commRing
  exact F.form_algebra_adj (hoeffdingOp F.gen u) f g

/-- **SG-1(b), orthogonality in the form.** -/
theorem form_component_of_ne {u w : Finset ι} (huw : u ≠ w) (f g : E) :
    F.form (F.component u f) (F.component w g) = 0 := by
  rw [form_component_adj, component_component_of_ne F huw, map_zero]

/-- The form against a component reads the component's variance. -/
theorem form_component_self (u : Finset ι) (f : E) :
    F.form f (F.component u f) = F.compVar u f := by
  rw [compVar, form_component_adj, component_idem]

/-- **SG-1(c).** `L P_u = −|u| P_u`. -/
theorem generator_component (u : Finset ι) (f : E) :
    F.generator (F.component u f) = -(u.card : ℝ) • F.component u f :=
  generator_eigen F.mean u _ (fun _ hi => F.mean_component_of_mem hi f)
    (fun _ hi => F.mean_component_of_notMem hi f)

/-- **SG-1.** Parseval: `⟨f, f⟩ = Σ_u V_u`. -/
theorem parseval (f : E) : F.form f f = ∑ u, F.compVar u f := by
  calc F.form f f = F.form (∑ u, F.component u f) f := by rw [F.sum_component]
    _ = ∑ u, F.form (F.component u f) f := by rw [map_sum, LinearMap.sum_apply]
    _ = ∑ u, F.compVar u f :=
        Finset.sum_congr rfl fun u _ => by rw [form_component_adj, form_component_self]

theorem compVar_nonneg (u : Finset ι) (f : E) : 0 ≤ F.compVar u f := F.form_nonneg _

/-- `total(v) = Σ_{u ∋ v} V_u`. -/
theorem total_eq_sum (v : ι) (f : E) :
    F.total v f = ∑ u, if v ∈ u then F.compVar u f else 0 := by
  have hexp : f - F.mean v f = ∑ u, if v ∈ u then F.component u f else 0 := by
    have hmean : F.mean v f = ∑ u, if v ∈ u then 0 else F.component u f := by
      conv_lhs => rw [← F.sum_component f]
      rw [map_sum]
      refine Finset.sum_congr rfl fun u _ => ?_
      split_ifs with hu
      · exact F.mean_component_of_mem hu f
      · exact F.mean_component_of_notMem hu f
    conv_lhs => rw [hmean]; lhs; rw [← F.sum_component f]
    rw [← Finset.sum_sub_distrib]
    refine Finset.sum_congr rfl fun u _ => ?_
    split_ifs <;> simp
  rw [total, hexp, map_sum]
  refine Finset.sum_congr rfl fun u _ => ?_
  split_ifs
  · exact F.form_component_self u f
  · simp

/-- **SG-1(d).** `Σ_v total(v) = Σ_u |u| V_u`. -/
theorem sum_total (f : E) : ∑ v, F.total v f = ∑ u : Finset ι, (u.card : ℝ) * F.compVar u f := by
  simp only [total_eq_sum]
  rw [Finset.sum_comm]
  refine Finset.sum_congr rfl fun u _ => ?_
  rw [← Finset.sum_filter, Finset.sum_const, nsmul_eq_mul]
  congr 2
  rw [Finset.filter_mem_eq_inter, Finset.univ_inter]

omit [DecidableEq ι] in
/-- **SG-1(d).** The Dirichlet form: `Σ_v total(v) = −⟨f, Lf⟩`. -/
theorem dirichlet_eq_sum_total (f : E) : -F.form f (F.generator f) = ∑ v, F.total v f := by
  rw [generator, map_sum, ← Finset.sum_neg_distrib]
  refine Finset.sum_congr rfl fun v _ => ?_
  rw [total, map_sub, map_sub]
  ring

/-- `Var f = ⟨f, f⟩ − V_∅ = Σ_{u ≠ ∅} V_u`. -/
theorem variance_eq_sum (f : E) :
    F.form f f - F.compVar ∅ f = ∑ u ∈ univ.erase (∅ : Finset ι), F.compVar u f := by
  rw [F.parseval f, ← Finset.add_sum_erase univ _ (mem_univ (∅ : Finset ι))]
  ring

/-- The excess of the Dirichlet form over the variance, a nonnegative sum. -/
theorem sum_total_sub_variance (f : E) :
    ∑ v, F.total v f - (F.form f f - F.compVar ∅ f) =
      ∑ u ∈ univ.erase (∅ : Finset ι), ((u.card : ℝ) - 1) * F.compVar u f := by
  rw [F.variance_eq_sum, F.sum_total, ← Finset.add_sum_erase univ _ (mem_univ (∅ : Finset ι))]
  rw [Finset.card_empty, Nat.cast_zero, zero_mul, zero_add, ← Finset.sum_sub_distrib]
  exact Finset.sum_congr rfl fun u _ => by ring

/-- **SG-1(d), Efron–Stein.** `Var f ≤ Σ_v total(v)`: the spectral gap of 1. -/
theorem efron_stein_form (f : E) : F.form f f - F.compVar ∅ f ≤ ∑ v, F.total v f := by
  have h := F.sum_total_sub_variance f
  have hnn : 0 ≤ ∑ u ∈ univ.erase (∅ : Finset ι), ((u.card : ℝ) - 1) * F.compVar u f := by
    refine Finset.sum_nonneg fun u hu => mul_nonneg ?_ (F.compVar_nonneg u f)
    have : 0 < u.card := Finset.card_pos.2 (Finset.nonempty_iff_ne_empty.2 (Finset.ne_of_mem_erase
      hu))
    have : (1 : ℝ) ≤ u.card := by exact_mod_cast this
    linarith
  linarith

/-- **SG-1(d), equality.** Efron–Stein is an equality exactly when `V_u = 0` for every `|u| ≥ 2`. -/
theorem efron_stein_form_eq_iff (f : E) :
    F.form f f - F.compVar ∅ f = ∑ v, F.total v f ↔ ∀ u : Finset ι, 2 ≤ u.card → F.compVar u f = 0
      := by
  have hexc := F.sum_total_sub_variance f
  have hterm : ∀ u ∈ univ.erase (∅ : Finset ι), 0 ≤ ((u.card : ℝ) - 1) * F.compVar u f := by
    intro u hu
    refine mul_nonneg ?_ (F.compVar_nonneg u f)
    have : 0 < u.card := Finset.card_pos.2 (Finset.nonempty_iff_ne_empty.2 (Finset.ne_of_mem_erase
      hu))
    have : (1 : ℝ) ≤ u.card := by exact_mod_cast this
    linarith
  constructor
  · intro heq u hu
    have hzero : ∑ u ∈ univ.erase (∅ : Finset ι), ((u.card : ℝ) - 1) * F.compVar u f = 0 := by
      rw [← hexc, ← heq, sub_self]
    have hne : u ∈ univ.erase (∅ : Finset ι) := by
      refine Finset.mem_erase.2 ⟨?_, mem_univ _⟩
      rintro rfl
      simp at hu
    have h := (Finset.sum_eq_zero_iff_of_nonneg hterm).1 hzero u hne
    have hc : ((u.card : ℝ) - 1) ≠ 0 := by
      have : (2 : ℝ) ≤ u.card := by exact_mod_cast hu
      linarith
    exact (mul_eq_zero.1 h).resolve_left hc
  · intro hzero
    have : ∑ u ∈ univ.erase (∅ : Finset ι), ((u.card : ℝ) - 1) * F.compVar u f = 0 := by
      refine Finset.sum_eq_zero fun u hu => ?_
      by_cases h2 : 2 ≤ u.card
      · rw [hzero u h2, mul_zero]
      · have hpos : 0 < u.card :=
          Finset.card_pos.2 (Finset.nonempty_iff_ne_empty.2 (Finset.ne_of_mem_erase hu))
        have hone : u.card = 1 := by omega
        simp [hone]
    linarith

/-- **SG-1(e), the operator form.** `U_ρ f = Σ_u ρ^|u| P_u f`. -/
theorem noise_apply (ρ : ℝ) (f : E) :
    F.noise ρ f = ∑ u : Finset ι, ρ ^ u.card • F.component u f := by
  letI := F.commRing
  have h := congrArg (fun a : F.algebraOf => ((a : E →ₗ[ℝ] E) f))
    (noise_eq_sum F.gen (algebraMap ℝ F.algebraOf ρ))
  simp only [AddSubmonoidClass.coe_finsetSum, LinearMap.sum_apply] at h
  refine h.trans (Finset.sum_congr rfl fun u _ => ?_)
  simp only [component, Subalgebra.coe_mul, Module.End.mul_apply,
    Subalgebra.coe_algebraMap, ← map_pow, Module.algebraMap_end_apply]

/-- **SG-1(e).** `⟨f, U_ρ f⟩ = Σ_u ρ^|u| V_u`: the coefficient of `ρ^k` is the variance at
interaction order `k`. -/
theorem form_noise (ρ : ℝ) (f : E) :
    F.form f (F.noise ρ f) = ∑ u : Finset ι, ρ ^ u.card * F.compVar u f := by
  rw [noise_apply, map_sum]
  refine Finset.sum_congr rfl fun u _ => ?_
  rw [map_smul, smul_eq_mul, form_component_self]

/-! ### SG-2(a) -/

omit [Fintype ι] [DecidableEq ι] in
/-- Cauchy–Schwarz for the nonnegative form. -/
theorem form_sq_le (x y : E) : F.form x y ^ 2 ≤ F.form x x * F.form y y := by
  have hexp : ∀ s t : ℝ, F.form (s • x + t • y) (s • x + t • y) =
      s ^ 2 * F.form x x + 2 * s * t * F.form x y + t ^ 2 * F.form y y := by
    intro s t
    simp only [map_add, map_smul, LinearMap.add_apply, LinearMap.smul_apply, smul_eq_mul,
      F.form_symm y x]
    ring
  have hy := F.form_nonneg y
  rcases hy.eq_or_lt with hy0 | hy0
  · have hxy : F.form x y = 0 := by
      by_contra hne
      set t : ℝ := -(F.form x x + 1) / (2 * F.form x y) with ht
      have h := F.form_nonneg ((1 : ℝ) • x + t • y)
      rw [hexp, ← hy0] at h
      have hc : 2 * (1 : ℝ) * t * F.form x y = -(F.form x x + 1) := by
        rw [ht]
        field_simp
      have hsum : (1 : ℝ) ^ 2 * F.form x x + 2 * 1 * t * F.form x y + t ^ 2 * 0 = -1 := by
        rw [hc]
        ring
      linarith
    rw [hxy, ← hy0]
    simp
  · have h := F.form_nonneg (F.form y y • x + (-F.form x y) • y)
    rw [hexp] at h
    have h' : 0 ≤ F.form y y * (F.form x x * F.form y y - F.form x y ^ 2) := by
      have : F.form y y ^ 2 * F.form x x + 2 * F.form y y * (-F.form x y) * F.form x y +
          (-F.form x y) ^ 2 * F.form y y =
          F.form y y * (F.form x x * F.form y y - F.form x y ^ 2) := by ring
      linarith
    have := nonneg_of_mul_nonneg_right h' hy0
    linarith

omit [Fintype ι] [DecidableEq ι] in
/-- The seminorm's triangle inequality. -/
theorem sqrt_form_add_le (x y : E) :
    √(F.form (x + y) (x + y)) ≤ √(F.form x x) + √(F.form y y) := by
  have hx := F.form_nonneg x
  have hy := F.form_nonneg y
  have hcs : F.form x y ≤ √(F.form x x) * √(F.form y y) := by
    rw [← Real.sqrt_mul hx]
    exact Real.le_sqrt_of_sq_le (F.form_sq_le x y)
  have hexp : F.form (x + y) (x + y) = F.form x x + 2 * F.form x y + F.form y y := by
    simp only [map_add, LinearMap.add_apply, F.form_symm y x]
    ring
  rw [hexp]
  apply Real.sqrt_le_iff.2
  refine ⟨by positivity, ?_⟩
  have hx' := Real.sq_sqrt hx
  have hy' := Real.sq_sqrt hy
  nlinarith

/-- `V_u h ≤ ⟨h, h⟩`. -/
theorem compVar_le_form (u : Finset ι) (h : E) : F.compVar u h ≤ F.form h h := by
  rw [F.parseval h]
  exact Finset.single_le_sum (fun w _ => F.compVar_nonneg w h) (mem_univ u)

/-- **SG-2(a).** `|√V_u(f) − √V_u(g)| ≤ ‖P_u(f − g)‖`. -/
theorem sqrt_compVar_sub_le (u : Finset ι) (f g : E) :
    |√(F.compVar u f) - √(F.compVar u g)| ≤ √(F.compVar u (f - g)) := by
  have hmap : ∀ a b : E, F.component u (a + b) = F.component u a + F.component u b := map_add _
  have hside : ∀ a b : E, √(F.compVar u a) ≤ √(F.compVar u (a - b)) + √(F.compVar u b) := by
    intro a b
    have h := F.sqrt_form_add_le (F.component u (a - b)) (F.component u b)
    rwa [← hmap, sub_add_cancel] at h
  have hsym : F.compVar u (g - f) = F.compVar u (f - g) := by
    have hc : F.component u (g - f) = -F.component u (f - g) := by rw [← map_neg, neg_sub]
    simp only [compVar, hc, map_neg, LinearMap.neg_apply, neg_neg]
  rw [abs_le]
  constructor
  · have := hside g f
    rw [hsym] at this
    linarith
  · have := hside f g
    linarith

/-- **SG-2(a).** So `√V_u` is 1-Lipschitz in the form's seminorm: a rate for `f` passes to it. -/
theorem sqrt_compVar_sub_le_form (u : Finset ι) (f g : E) :
    |√(F.compVar u f) - √(F.compVar u g)| ≤ √(F.form (f - g) (f - g)) :=
  (F.sqrt_compVar_sub_le u f g).trans (Real.sqrt_le_sqrt (F.compVar_le_form u _))

end CommutingMeans

/-! ### SG-1 on a finite product -/

section Finite

variable {ι : Type*} [Fintype ι] [DecidableEq ι] {X : ι → Type*} [∀ i, Fintype (X i)]

/-- `E_i`, averaging out coordinate `i`. -/
noncomputable def coordMean (μ : ∀ i, X i → ℝ) (i : ι) (f : (∀ j, X j) → ℝ) :
    (∀ j, X j) → ℝ :=
  fun x => ∑ y : X i, μ i y * f (Function.update x i y)

/-- The mean under the product law. -/
noncomputable def meanOf (μ : ∀ i, X i → ℝ) (f : (∀ j, X j) → ℝ) : ℝ :=
  ∑ x : ∀ j, X j, (∏ i, μ i (x i)) * f x

omit [Fintype ι] in
/-- `E_i f` does not read coordinate `i`. -/
theorem coordMean_update (μ : ∀ i, X i → ℝ) (i : ι) (f : (∀ j, X j) → ℝ) (x : ∀ j, X j) (y : X i) :
    coordMean μ i f (Function.update x i y) = coordMean μ i f x := by
  simp only [coordMean, Function.update_idem]

omit [Fintype ι] in
/-- **SG-1(a).** `E_i` is idempotent, for a probability vector. -/
theorem coordMean_idem (μ : ∀ i, X i → ℝ) (hμ : ∀ i, ∑ y, μ i y = 1) (i : ι)
    (f : (∀ j, X j) → ℝ) : coordMean μ i (coordMean μ i f) = coordMean μ i f := by
  funext x
  have h : coordMean μ i (coordMean μ i f) x = ∑ y, μ i y * coordMean μ i f x :=
    Finset.sum_congr rfl fun y _ => by rw [coordMean_update]
  rw [h, ← Finset.sum_mul, hμ, one_mul]

omit [Fintype ι] in
/-- **SG-1(a).** The `E_i` commute: `Finset.sum_comm` with `Function.update_comm`. -/
theorem coordMean_comm (μ : ∀ i, X i → ℝ) (i j : ι) (f : (∀ j, X j) → ℝ) :
    coordMean μ i (coordMean μ j f) = coordMean μ j (coordMean μ i f) := by
  by_cases hij : i = j
  · subst hij
    rfl
  funext x
  simp only [coordMean, Finset.mul_sum]
  rw [Finset.sum_comm]
  refine Finset.sum_congr rfl fun b _ => Finset.sum_congr rfl fun a _ => ?_
  rw [Function.update_comm (Ne.symm hij)]
  ring

omit [∀ i, Fintype (X i)] in
/-- Changing one coordinate and its factor in the product law. -/
private theorem weight_update (μ : ∀ i, X i → ℝ) (i : ι) (x : ∀ j, X j) (y : X i) :
    (∏ j, μ j (Function.update x i y j)) * μ i (x i) = (∏ j, μ j (x j)) * μ i y := by
  rw [← Finset.mul_prod_erase univ _ (mem_univ i),
    ← Finset.mul_prod_erase univ (fun j => μ j (x j)) (mem_univ i)]
  have h : ∏ j ∈ univ.erase i, μ j (Function.update x i y j) = ∏ j ∈ univ.erase i, μ j (x j) :=
    Finset.prod_congr rfl fun j hj => by rw [Function.update_of_ne (Finset.ne_of_mem_erase hj)]
  rw [h, Function.update_self]
  ring

/-- The swap `(x, y) ↦ (x[i := y], x i)`, an involution. -/
private def swapCoord (i : ι) : ((∀ j, X j) × X i) ≃ ((∀ j, X j) × X i) where
  toFun p := (Function.update p.1 i p.2, p.1 i)
  invFun p := (Function.update p.1 i p.2, p.1 i)
  left_inv p := by
    simp [Function.update_idem, Function.update_eq_self]
  right_inv p := by
    simp [Function.update_idem, Function.update_eq_self]

/-- Averaging out a coordinate keeps the mean. -/
theorem meanOf_coordMean (μ : ∀ i, X i → ℝ) (hμ : ∀ i, ∑ y, μ i y = 1) (i : ι)
    (f : (∀ j, X j) → ℝ) : meanOf μ (coordMean μ i f) = meanOf μ f := by
  have hswap := Equiv.sum_comp (swapCoord (X := X) i)
    (fun p => (∏ j, μ j (p.1 j)) * (μ i p.2 * f (Function.update p.1 i p.2)))
  have hl : meanOf μ (coordMean μ i f) =
      ∑ p : (∀ j, X j) × X i, (∏ j, μ j (p.1 j)) * (μ i p.2 * f (Function.update p.1 i p.2)) := by
    rw [Fintype.sum_prod_type]
    exact Finset.sum_congr rfl fun x _ => by rw [coordMean, Finset.mul_sum]
  have hr : meanOf μ f =
      ∑ p : (∀ j, X j) × X i, (∏ j, μ j (Function.update p.1 i p.2 j)) *
        (μ i (p.1 i) * f p.1) := by
    rw [Fintype.sum_prod_type]
    refine Finset.sum_congr rfl fun x _ => ?_
    have : ∀ y, (∏ j, μ j (Function.update x i y j)) * (μ i (x i) * f x) =
        (∏ j, μ j (x j)) * f x * μ i y := by
      intro y
      have := weight_update μ i x y
      linear_combination f x * this
    simp only [this]
    rw [← Finset.mul_sum, hμ, mul_one]
  rw [hl, hr, ← hswap]
  refine Finset.sum_congr rfl fun p _ => ?_
  simp [swapCoord, Function.update_idem, Function.update_eq_self]

/-- **SG-1(a).** `E_i` is self-adjoint under the product law. -/
theorem coordMean_selfAdjoint (μ : ∀ i, X i → ℝ) (hμ : ∀ i, ∑ y, μ i y = 1) (i : ι)
    (f g : (∀ j, X j) → ℝ) :
    meanOf μ (fun x => f x * coordMean μ i g x) = meanOf μ (fun x => coordMean μ i f x * g x) := by
  have hpull : ∀ a b : (∀ j, X j) → ℝ, coordMean μ i (fun x => a x * coordMean μ i b x) =
      fun x => coordMean μ i a x * coordMean μ i b x := by
    intro a b
    funext x
    simp only [coordMean, Finset.sum_mul, mul_assoc]
    exact Finset.sum_congr rfl fun _ _ => by
      rw [show (∑ y, μ i y * b (Function.update (Function.update x i _) i y)) =
        ∑ y, μ i y * b (Function.update x i y) by
          simp only [Function.update_idem]]
  rw [← meanOf_coordMean μ hμ i (fun x => f x * coordMean μ i g x), hpull,
    ← meanOf_coordMean μ hμ i (fun x => coordMean μ i f x * g x)]
  have hsym : coordMean μ i (fun x => coordMean μ i f x * g x) =
      coordMean μ i (fun x => g x * coordMean μ i f x) := by
    congr 1
    funext x
    exact mul_comm _ _
  rw [hsym, hpull]
  simp only [mul_comm]

/-- **SG-1(a).** `E_i` as a linear map. -/
noncomputable def coordMeanLinear (μ : ∀ i, X i → ℝ) (i : ι) :
    ((∀ j, X j) → ℝ) →ₗ[ℝ] ((∀ j, X j) → ℝ) where
  toFun := coordMean μ i
  map_add' f g := by
    funext x
    simp only [coordMean, Pi.add_apply, mul_add, Finset.sum_add_distrib]
  map_smul' c f := by
    funext x
    simp only [coordMean, Pi.smul_apply, smul_eq_mul, RingHom.id_apply, Finset.mul_sum]
    exact Finset.sum_congr rfl fun _ _ => by ring

/-- The product law's form. -/
noncomputable def meanForm (μ : ∀ i, X i → ℝ) :
    ((∀ j, X j) → ℝ) →ₗ[ℝ] ((∀ j, X j) → ℝ) →ₗ[ℝ] ℝ :=
  LinearMap.mk₂ ℝ (fun f g => meanOf μ (fun x => f x * g x))
    (fun a b c => by
      simp only [meanOf, Pi.add_apply, add_mul, mul_add, Finset.sum_add_distrib])
    (fun r a b => by
      simp only [meanOf, Pi.smul_apply, smul_eq_mul, Finset.mul_sum]
      exact Finset.sum_congr rfl fun _ _ => by ring)
    (fun a b c => by
      simp only [meanOf, Pi.add_apply, mul_add, Finset.sum_add_distrib])
    (fun r a b => by
      simp only [meanOf, Pi.smul_apply, smul_eq_mul, Finset.mul_sum]
      exact Finset.sum_congr rfl fun _ _ => by ring)

/-- The finite product's commuting means. -/
noncomputable def productMeans (μ : ∀ i, X i → ℝ) (h₀ : ∀ i y, 0 ≤ μ i y)
    (hμ : ∀ i, ∑ y, μ i y = 1) : CommutingMeans ι ((∀ j, X j) → ℝ) where
  mean := coordMeanLinear μ
  idem i f := coordMean_idem μ hμ i f
  comm i j f := coordMean_comm μ i j f
  form := meanForm μ
  form_symm f g := by
    simp only [meanForm, LinearMap.mk₂_apply, mul_comm (f _)]
  form_adj i f g := (coordMean_selfAdjoint μ hμ i f g).symm
  form_nonneg f := by
    simp only [meanForm, LinearMap.mk₂_apply, meanOf]
    exact Finset.sum_nonneg fun x _ => mul_nonneg
      (Finset.prod_nonneg fun i _ => h₀ i _) (mul_self_nonneg _)

/-- The product law has total mass 1. -/
theorem sum_prod_mu (μ : ∀ i, X i → ℝ) (hμ : ∀ i, ∑ y, μ i y = 1) :
    ∑ x : ∀ j, X j, ∏ i, μ i (x i) = 1 := by
  have h := Finset.prod_univ_sum (t := fun i => (univ : Finset (X i))) (f := fun i y => μ i y)
  rw [Fintype.piFinset_univ] at h
  rw [← h]
  simp [hμ]

omit [∀ i, Fintype (X i)] in
set_option linter.unusedFintypeInType false in
/-- A function that no single coordinate moves is constant. -/
theorem eq_of_update_eq (g : (∀ j, X j) → ℝ)
    (h : ∀ (v : ι) (x : ∀ j, X j) (y : X v), g (Function.update x v y) = g x) (x x' : ∀ j, X j) :
    g x' = g x := by
  have hs : ∀ s : Finset ι, g (fun j => if j ∈ s then x' j else x j) = g x := by
    intro s
    induction s using Finset.induction_on with
    | empty => simp
    | insert a s _ ih =>
      have hup : (fun j => if j ∈ insert a s then x' j else x j) =
          Function.update (fun j => if j ∈ s then x' j else x j) a (x' a) := by
        funext j
        by_cases hj : j = a
        · subst hj
          simp
        · simp [Finset.mem_insert, hj]
      rw [hup, h, ih]
  simpa using hs univ

/-- `P_∅ f` is the constant `E f`, so `V_∅ = (E f)²`. -/
theorem compVar_empty_productMeans (μ : ∀ i, X i → ℝ) (h₀ : ∀ i y, 0 ≤ μ i y)
    (hμ : ∀ i, ∑ y, μ i y = 1) (f : (∀ j, X j) → ℝ) :
    (productMeans μ h₀ hμ).compVar ∅ f = meanOf μ f ^ 2 := by
  set F := productMeans μ h₀ hμ with hF
  set g := F.component ∅ f with hg
  have hone : ∀ i, coordMean μ i (fun _ : ∀ j, X j => (1 : ℝ)) = fun _ => 1 := by
    intro i
    funext x
    simp [coordMean, hμ]
  have hfix : ∀ (v : ι) (x : ∀ j, X j) (y : X v), g (Function.update x v y) = g x := by
    intro v x y
    have hv : coordMean μ v g = g :=
      F.mean_component_of_notMem (u := ∅) (Finset.notMem_empty v) f
    rw [← hv, coordMean_update]
  have hconst : ∀ x x', g x' = g x := eq_of_update_eq g hfix
  have hoff : ∀ u : Finset ι, u ≠ ∅ → F.component u (fun _ => (1 : ℝ)) = 0 := by
    intro u hu
    obtain ⟨v, hv⟩ := Finset.nonempty_iff_ne_empty.2 hu
    have h := F.mean_component_of_mem hv (fun _ => (1 : ℝ))
    rwa [← F.component_mean, show F.mean v (fun _ => (1 : ℝ)) = fun _ => 1 from hone v] at h
  have hone' : F.component ∅ (fun _ => (1 : ℝ)) = fun _ => 1 := by
    have h := F.sum_component (fun _ => (1 : ℝ))
    rwa [Finset.sum_eq_single ∅ (fun u _ hu => hoff u hu) (by simp)] at h
  have hmean : meanOf μ f = meanOf μ g := by
    have h := F.form_component_adj ∅ f (fun _ => (1 : ℝ))
    rw [hone'] at h
    simp only [hF, productMeans, meanForm, LinearMap.mk₂_apply, mul_one] at h
    exact h.symm
  by_cases hne : Nonempty (∀ j, X j)
  · obtain ⟨x₀⟩ := hne
    have hgc : g = fun _ => g x₀ := funext fun x => hconst x₀ x
    have hmg : meanOf μ g = g x₀ := by
      rw [hgc, meanOf, ← Finset.sum_mul, sum_prod_mu μ hμ, one_mul]
    rw [hmean, hmg, CommutingMeans.compVar, ← hg]
    simp only [hF, productMeans, meanForm, LinearMap.mk₂_apply]
    rw [hgc, meanOf, ← Finset.sum_mul, sum_prod_mu μ hμ, one_mul]
    ring
  · exfalso
    have : ∑ x : ∀ j, X j, ∏ i, μ i (x i) = 1 := sum_prod_mu μ hμ
    have hempty : IsEmpty (∀ j, X j) := not_nonempty_iff.1 hne
    simp at this

/-- **SG-1(d), Efron–Stein on a finite product.** -/
theorem efron_stein (μ : ∀ i, X i → ℝ) (h₀ : ∀ i y, 0 ≤ μ i y) (hμ : ∀ i, ∑ y, μ i y = 1)
    (f : (∀ j, X j) → ℝ) :
    meanOf μ (fun x => f x ^ 2) - meanOf μ f ^ 2 ≤
      ∑ i, meanOf μ (fun x => f x * (f x - coordMean μ i f x)) := by
  have h := (productMeans μ h₀ hμ).efron_stein_form f
  rw [compVar_empty_productMeans] at h
  have hl : (productMeans μ h₀ hμ).form f f = meanOf μ (fun x => f x ^ 2) := by
    simp only [productMeans, meanForm, LinearMap.mk₂_apply, sq]
  rw [hl] at h
  exact h

omit [Fintype ι] in
/-- `E_w` fixes a function of another coordinate. -/
theorem coordMean_of_coord_ne (μ : ∀ i, X i → ℝ) (hμ : ∀ i, ∑ y, μ i y = 1) {w v : ι} (hwv : w ≠ v)
    (h : X v → ℝ) : coordMean μ w (fun x => h (x v)) = fun x => h (x v) := by
  funext x
  have hup : ∀ y, h (Function.update x w y v) = h (x v) := fun y => by
    rw [Function.update_of_ne (Ne.symm hwv)]
  simp only [coordMean, hup]
  rw [← Finset.sum_mul, hμ, one_mul]

omit [Fintype ι] in
/-- `E_w` fixes a constant. -/
theorem coordMean_const (μ : ∀ i, X i → ℝ) (hμ : ∀ i, ∑ y, μ i y = 1) (w : ι) (c : ℝ) :
    coordMean μ w (fun _ : ∀ j, X j => c) = fun _ => c := by
  funext x
  simp only [coordMean]
  rw [← Finset.sum_mul, hμ, one_mul]

omit [∀ i, Fintype (X i)] in
set_option linter.unusedFintypeInType false in
/-- A function no coordinate but `v` moves reads only `x v`. -/
theorem eq_of_update_eq_of_ne (g : (∀ j, X j) → ℝ) (v : ι)
    (h : ∀ (w : ι), w ≠ v → ∀ (x : ∀ j, X j) (y : X w), g (Function.update x w y) = g x)
    {x x' : ∀ j, X j} (hv : x v = x' v) : g x' = g x := by
  have hs : ∀ s : Finset ι, g (fun j => if j ∈ s then x' j else x j) = g x := by
    intro s
    induction s using Finset.induction_on with
    | empty => simp
    | insert a s ha ih =>
      by_cases hav : a = v
      · subst hav
        have hsame : (fun j => if j ∈ insert a s then x' j else x j) =
            fun j => if j ∈ s then x' j else x j := by
          funext j
          by_cases hj : j = a
          · subst hj
            simp [ha, hv]
          · simp [Finset.mem_insert, hj]
        rw [hsame, ih]
      · have hup : (fun j => if j ∈ insert a s then x' j else x j) =
            Function.update (fun j => if j ∈ s then x' j else x j) a (x' a) := by
          funext j
          by_cases hj : j = a
          · subst hj
            simp
          · simp [Finset.mem_insert, hj]
        rw [hup, h a hav, ih]
  simpa using hs univ

omit [DecidableEq ι] in
set_option linter.unusedFintypeInType false in
/-- A probability vector on every coordinate makes the product inhabited. -/
theorem nonempty_of_prob (μ : ∀ i, X i → ℝ) (hμ : ∀ i, ∑ y, μ i y = 1) : Nonempty (∀ j, X j) := by
  classical
  by_contra h
  have htot := sum_prod_mu μ hμ
  rw [not_nonempty_iff] at h
  simp at htot

/-- **SG-1(d), the equality.** Under full support, Efron–Stein is an equality exactly when `f` is
additive, a constant plus functions of one coordinate each. -/
theorem efron_stein_eq_iff_additive (μ : ∀ i, X i → ℝ) (hpos : ∀ i y, 0 < μ i y)
    (hμ : ∀ i, ∑ y, μ i y = 1) (f : (∀ j, X j) → ℝ) :
    meanOf μ (fun x => f x ^ 2) - meanOf μ f ^ 2 =
        ∑ i, meanOf μ (fun x => f x * (f x - coordMean μ i f x)) ↔
      ∃ (c : ℝ) (h : ∀ v, X v → ℝ), ∀ x, f x = c + ∑ v, h v (x v) := by
  have h₀ : ∀ i y, 0 ≤ μ i y := fun i y => (hpos i y).le
  set F := productMeans μ h₀ hμ with hF
  have hl : F.form f f = meanOf μ (fun x => f x ^ 2) := by
    simp only [hF, productMeans, meanForm, LinearMap.mk₂_apply, sq]
  have hkey : (meanOf μ (fun x => f x ^ 2) - meanOf μ f ^ 2 =
        ∑ i, meanOf μ (fun x => f x * (f x - coordMean μ i f x))) ↔
      ∀ u : Finset ι, 2 ≤ u.card → F.compVar u f = 0 := by
    have hc : F.compVar ∅ f = meanOf μ f ^ 2 := compVar_empty_productMeans μ h₀ hμ f
    rw [← F.efron_stein_form_eq_iff f, hl, hc]
    rfl
  rw [hkey]
  have hzero : ∀ g : (∀ j, X j) → ℝ, F.form g g = 0 → g = 0 := by
    intro g hg
    have hg' : ∑ x, (∏ i, μ i (x i)) * (g x * g x) = 0 := hg
    have hterm := (Finset.sum_eq_zero_iff_of_nonneg fun x _ => mul_nonneg
      (Finset.prod_nonneg fun i _ => h₀ i _) (mul_self_nonneg _)).1 hg'
    funext x
    have hw : 0 < ∏ i, μ i (x i) := Finset.prod_pos fun i _ => hpos i _
    have := (mul_eq_zero.1 (hterm x (Finset.mem_univ x))).resolve_left hw.ne'
    exact mul_self_eq_zero.1 this
  constructor
  · intro hV
    have hcomp : ∀ u : Finset ι, 2 ≤ u.card → F.component u f = 0 :=
      fun u hu => hzero _ (hV u hu)
    -- The decomposition stops at singletons.
    have hsum : f = F.component ∅ f + ∑ v, F.component {v} f := by
      have hfull := F.sum_component f
      have hsub : (insert (∅ : Finset ι) (univ.image fun v => ({v} : Finset ι))) ⊆ univ :=
        Finset.subset_univ _
      rw [← Finset.sum_subset hsub] at hfull
      · rw [Finset.sum_insert, Finset.sum_image] at hfull
        · exact hfull.symm
        · intro a _ b _ hab
          exact Finset.singleton_injective hab
        · simp
      · intro u _ hu
        apply hcomp u
        by_contra hlt
        have hne : u ≠ ∅ := fun he => hu (by simp [he])
        obtain ⟨a, ha⟩ := Finset.nonempty_iff_ne_empty.2 hne
        have hle : u.card ≤ 1 := by omega
        have hsing : u = {a} := by
          refine Finset.eq_singleton_iff_unique_mem.2 ⟨ha, fun b hb => ?_⟩
          exact Finset.card_le_one.1 hle b hb a ha
        exact hu (Finset.mem_insert_of_mem (Finset.mem_image.2 ⟨a, Finset.mem_univ a, hsing.symm⟩))
    obtain ⟨x₀⟩ := nonempty_of_prob μ hμ
    -- Each singleton's component reads only its own coordinate, and `P_∅ f` none.
    have hE : ∀ (g : (∀ j, X j) → ℝ) (w : ι), coordMean μ w g = g →
        ∀ (x : ∀ j, X j) (y : X w), g (Function.update x w y) = g x := by
      intro g w hg x y
      rw [← hg, coordMean_update]
    have hconst : ∀ x, F.component ∅ f x = F.component ∅ f x₀ := fun x =>
      eq_of_update_eq (F.component ∅ f) (fun v =>
        hE _ v (F.mean_component_of_notMem (u := ∅) (Finset.notMem_empty v) f)) x₀ x
    have hsingle : ∀ (v : ι) (x : ∀ j, X j),
        F.component {v} f x = F.component {v} f (Function.update x₀ v (x v)) := by
      intro v x
      refine (eq_of_update_eq_of_ne (F.component {v} f) v (fun w hwv => ?_) (by simp)).symm
      exact hE _ w (F.mean_component_of_notMem (u := {v}) (by simpa using hwv) f)
    refine ⟨F.component ∅ f x₀, fun v y => F.component {v} f (Function.update x₀ v y), fun x => ?_⟩
    have hx := congrFun hsum x
    rw [Pi.add_apply, Finset.sum_apply, hconst x] at hx
    rw [hx]
    congr 1
    exact Finset.sum_congr rfl fun v _ => hsingle v x
  · rintro ⟨c, h, hf⟩ u hu
    have hconst : F.component u (fun _ => c) = 0 := by
      obtain ⟨w, hw⟩ : u.Nonempty := Finset.card_pos.1 (by omega)
      have h0 := F.mean_component_of_mem hw (fun _ => c)
      rw [← F.component_mean] at h0
      rwa [show F.mean w (fun _ => c) = fun _ => c from coordMean_const μ hμ w c] at h0
    have hcoord : ∀ v, F.component u (fun x => h v (x v)) = 0 := by
      intro v
      obtain ⟨w, hw, hwv⟩ := Finset.exists_mem_ne (by omega : 1 < u.card) v
      have h0 := F.mean_component_of_mem hw (fun x => h v (x v))
      rw [← F.component_mean] at h0
      rwa [show F.mean w (fun x => h v (x v)) = fun x => h v (x v) from
        coordMean_of_coord_ne μ hμ hwv (h v)] at h0
    have hfeq : f = (fun _ => c) + ∑ v, fun x => h v (x v) := by
      funext x
      rw [hf x]
      simp [Finset.sum_apply]
    have hcompf : F.component u f = 0 := by
      rw [hfeq, map_add, map_sum, hconst, Finset.sum_eq_zero fun v _ => hcoord v, add_zero]
    simp [CommutingMeans.compVar, hcompf]

/-- The mean of a product-form function factors. -/
theorem meanOf_prod (μ : ∀ i, X i → ℝ) (g : ∀ i, X i → ℝ) :
    meanOf μ (fun x => ∏ i, g i (x i)) = ∏ i, ∑ y, μ i y * g i y := by
  have h := Finset.prod_univ_sum (t := fun i => (univ : Finset (X i))) (f := fun i y => μ i y * g i
    y)
  rw [Fintype.piFinset_univ] at h
  rw [h, meanOf]
  exact Finset.sum_congr rfl fun x _ => by rw [Finset.prod_mul_distrib]

/-- `E_v` on a product: it replaces the `v`-th factor by its mean. -/
theorem coordMean_prod (μ : ∀ i, X i → ℝ) (v : ι) (b : ∀ i, X i → ℝ) (x : ∀ j, X j) :
    coordMean μ v (fun z => ∏ i, b i (z i)) x =
      (∑ y, μ v y * b v y) * ∏ i ∈ univ.erase v, b i (x i) := by
  rw [coordMean, Finset.sum_mul]
  refine Finset.sum_congr rfl fun y _ => ?_
  rw [← Finset.mul_prod_erase univ _ (mem_univ v), Function.update_self]
  have : ∏ i ∈ univ.erase v, b i (Function.update x v y i) = ∏ i ∈ univ.erase v, b i (x i) :=
    Finset.prod_congr rfl fun i hi => by rw [Function.update_of_ne (Finset.ne_of_mem_erase hi)]
  rw [this]
  ring

/-- The noise factor at `v` on a product-form function, as a product-form function. -/
theorem noiseFactor_prod (μ : ∀ i, X i → ℝ) (v : ι) (ρ : ℝ) (b : ∀ i, X i → ℝ) :
    (fun x => ρ * (∏ i, b i (x i)) + (1 - ρ) * coordMean μ v (fun z => ∏ i, b i (z i)) x) =
      fun x => ∏ i, Function.update b v
        (fun y => ρ * b v y + (1 - ρ) * ∑ t, μ v t * b v t) i (x i) := by
  funext x
  rw [coordMean_prod, ← Finset.mul_prod_erase univ _ (mem_univ v),
    ← Finset.mul_prod_erase univ (fun i => Function.update b v _ i (x i)) (mem_univ v)]
  simp only [Function.update_self]
  have : ∏ i ∈ univ.erase v, Function.update b v
      (fun y => ρ * b v y + (1 - ρ) * ∑ t, μ v t * b v t) i (x i) =
      ∏ i ∈ univ.erase v, b i (x i) :=
    Finset.prod_congr rfl fun i hi => by rw [Function.update_of_ne (Finset.ne_of_mem_erase hi)]
  rw [this]
  ring

/-- `Π_{v∈s}(ρI + (1 − ρ)E_v)` on a product-form function replaces each factor in `s`. -/
theorem noisePart_prod (μ : ∀ i, X i → ℝ) (h₀ : ∀ i y, 0 ≤ μ i y) (hμ : ∀ i, ∑ y, μ i y = 1)
    (ρ : ℝ) (b : ∀ i, X i → ℝ) (s : Finset ι) :
    (productMeans μ h₀ hμ).noisePart ρ s (fun x => ∏ i, b i (x i)) =
      fun x => ∏ i, (if i ∈ s then fun y => ρ * b i y + (1 - ρ) * ∑ t, μ i t * b i t
        else b i) (x i) := by
  induction s using Finset.induction_on with
  | empty =>
    rw [CommutingMeans.noisePart_empty]
    funext x
    simp
  | insert a s ha ih =>
    rw [CommutingMeans.noisePart_insert _ _ ha, ih]
    have hf := noiseFactor_prod μ a ρ
      (fun i => if i ∈ s then fun y => ρ * b i y + (1 - ρ) * ∑ t, μ i t * b i t else b i)
    have hb : (if a ∈ s then fun y => ρ * b a y + (1 - ρ) * ∑ t, μ a t * b a t else b a) = b a :=
      if_neg ha
    rw [hb] at hf
    refine Eq.trans ?_ (hf.trans ?_)
    · funext x
      rfl
    · funext x
      refine Finset.prod_congr rfl fun i _ => ?_
      by_cases hia : i = a
      · subst hia
        simp
      · simp [Finset.mem_insert, hia]

/-- **SG-1(e), the factorization.** For product-form `a` and `b`,
`⟨a, U_ρ b⟩ = Π_v (E a_v E b_v + ρ Cov(a_v, b_v))`. So the pair polynomial needs only each
variable's
`E a_v`, `E b_v` and `E[a_v b_v]`. -/
theorem noise_prod_form (μ : ∀ i, X i → ℝ) (h₀ : ∀ i y, 0 ≤ μ i y) (hμ : ∀ i, ∑ y, μ i y = 1)
    (ρ : ℝ) (a b : ∀ i, X i → ℝ) :
    (productMeans μ h₀ hμ).form (fun x => ∏ i, a i (x i))
        ((productMeans μ h₀ hμ).noise ρ (fun x => ∏ i, b i (x i))) =
      ∏ v, ((∑ y, μ v y * a v y) * (∑ y, μ v y * b v y) +
        ρ * (∑ y, μ v y * (a v y * b v y) - (∑ y, μ v y * a v y) * (∑ y, μ v y * b v y))) := by
  rw [CommutingMeans.noise, noisePart_prod]
  simp only [Finset.mem_univ, if_true, productMeans, meanForm, LinearMap.mk₂_apply,
    ← Finset.prod_mul_distrib]
  rw [meanOf_prod μ (fun i y => a i y * (ρ * b i y + (1 - ρ) * ∑ t, μ i t * b i t))]
  refine Finset.prod_congr rfl fun v _ => ?_
  set B := ∑ t, μ v t * b v t with hB
  have h : ∑ y, μ v y * (a v y * (ρ * b v y + (1 - ρ) * B)) =
      ρ * ∑ y, μ v y * (a v y * b v y) + (1 - ρ) * (∑ y, μ v y * a v y) * B := by
    calc ∑ y, μ v y * (a v y * (ρ * b v y + (1 - ρ) * B))
        = ∑ y, (ρ * (μ v y * (a v y * b v y)) + (1 - ρ) * B * (μ v y * a v y)) :=
          Finset.sum_congr rfl fun y _ => by ring
      _ = ρ * ∑ y, μ v y * (a v y * b v y) + (1 - ρ) * (∑ y, μ v y * a v y) * B := by
          rw [Finset.sum_add_distrib, ← Finset.mul_sum, ← Finset.mul_sum]
          ring
  rw [h]
  ring

end Finite

/-! ### SG-3 -/

/-- **SG-3, the finite form.** Under a joint mixture of products, a product-form term's integral is
the mixture of the per-component products of one-variable integrals. -/
theorem mixture_product_integral {ι κ : Type*} [Fintype ι] [DecidableEq ι] [Fintype κ]
    {X : ι → Type*} [∀ i, Fintype (X i)]
    (w : κ → ℝ) (μ : κ → ∀ i, X i → ℝ) (a : ∀ i, X i → ℝ) :
    ∑ x : ∀ j, X j, (∑ c, w c * ∏ i, μ c i (x i)) * ∏ i, a i (x i) =
      ∑ c, w c * ∏ i, ∑ y, μ c i y * a i y := by
  simp only [Finset.sum_mul]
  rw [Finset.sum_comm]
  refine Finset.sum_congr rfl fun c _ => ?_
  have h := Finset.prod_univ_sum (t := fun i => (univ : Finset (X i)))
    (f := fun i y => μ c i y * a i y)
  rw [Fintype.piFinset_univ] at h
  rw [h, Finset.mul_sum]
  refine Finset.sum_congr rfl fun x _ => ?_
  rw [Finset.prod_mul_distrib]
  ring

/-- **SG-3, the pair-moment corollary**, the case `a = φ_j · φ_l` the claim is stated for. -/
theorem mixture_pair_integral {ι κ : Type*} [Fintype ι] [DecidableEq ι] [Fintype κ]
    {X : ι → Type*} [∀ i, Fintype (X i)]
    (w : κ → ℝ) (μ : κ → ∀ i, X i → ℝ) (φ ψ : ∀ i, X i → ℝ) :
    ∑ x : ∀ j, X j, (∑ c, w c * ∏ i, μ c i (x i)) * ((∏ i, φ i (x i)) * ∏ i, ψ i (x i)) =
      ∑ c, w c * ∏ i, ∑ y, μ c i y * (φ i y * ψ i y) := by
  rw [← mixture_product_integral w μ (fun i y => φ i y * ψ i y)]
  exact Finset.sum_congr rfl fun x _ => by rw [Finset.prod_mul_distrib]

/-! ### SG-4 -/

section Legendre

open Polynomial

/-- The recurrence on the coefficients of `P_n`: `(k+1)² c_{k+1} = −(n−k)(n+k+1) c_k`. -/
theorem coeff_shiftedLegendre_succ (n k : ℕ) :
    ((k : ℤ) + 1) ^ 2 * (shiftedLegendre n).coeff (k + 1) =
      -(((n : ℤ) - k) * ((n : ℤ) + k + 1)) * (shiftedLegendre n).coeff k := by
  rw [coeff_shiftedLegendre, coeff_shiftedLegendre]
  rcases le_or_gt k n with hkn | hkn
  · have hsub : ((n - k : ℕ) : ℤ) = (n : ℤ) - k := by omega
    have h₁ : (n.choose (k + 1) : ℤ) * ((k : ℤ) + 1) = (n.choose k : ℤ) * ((n : ℤ) - k) := by
      rw [← hsub]
      exact_mod_cast Nat.choose_succ_right_eq n k
    have h₂ : (((n + k).choose n : ℕ) : ℤ) * ((n : ℤ) + k + 1) =
        ((n + (k + 1)).choose n : ℤ) * ((k : ℤ) + 1) := by
      have h := Nat.choose_mul_succ_eq (n + k) n
      have hk : n + k + 1 - n = k + 1 := by omega
      rw [hk, show n + k + 1 = n + (k + 1) by ring] at h
      have h' : (((n + k).choose n : ℕ) : ℤ) * ((n + (k + 1) : ℕ) : ℤ) =
          ((n + (k + 1)).choose n : ℤ) * ((k + 1 : ℕ) : ℤ) := by exact_mod_cast h
      push_cast at h'
      linarith
    linear_combination
      (-(-1 : ℤ) ^ k * ((k : ℤ) + 1) * ((n + (k + 1)).choose n : ℤ)) * h₁ -
        (-(-1 : ℤ) ^ k * ((n : ℤ) - k) * (n.choose k : ℤ)) * h₂
  · have h₁ : n.choose (k + 1) = 0 := Nat.choose_eq_zero_of_lt (by omega)
    have h₂ : n.choose k = 0 := Nat.choose_eq_zero_of_lt hkn
    simp [h₁, h₂]

/-- **SG-4(a).** `[x(1−x) P_n′]′ = −n(n+1) P_n`. The relation does not see Mathlib's sign. -/
theorem shiftedLegendre_jacobi (n : ℕ) :
    derivative ((X * (1 - X)) * derivative (shiftedLegendre n)) =
      -C ((n : ℤ) * (n + 1)) * shiftedLegendre n := by
  ext k
  have hrec := coeff_shiftedLegendre_succ n k
  have hexp : ((X * (1 - X)) * derivative (shiftedLegendre n)) =
      X * derivative (shiftedLegendre n) - X ^ 2 * derivative (shiftedLegendre n) := by ring
  have hc₁ : (X * derivative (shiftedLegendre n)).coeff (k + 1) =
      (shiftedLegendre n).coeff (k + 1) * ((k : ℤ) + 1) := by
    rw [coeff_X_mul, coeff_derivative]
  have hc₂ : (X ^ 2 * derivative (shiftedLegendre n)).coeff (k + 1) =
      (shiftedLegendre n).coeff k * (k : ℤ) := by
    cases k with
    | zero => simp [coeff_X_pow_mul']
    | succ j =>
      rw [show j + 1 + 1 = j + 2 by ring, coeff_X_pow_mul, coeff_derivative]
      push_cast
      ring
  rw [coeff_derivative, hexp, coeff_sub, hc₁, hc₂, neg_mul, coeff_neg, coeff_C_mul]
  linear_combination hrec

end Legendre

/-! ### SG-5 -/

section Rayleigh

/-- A symmetric 2 × 2 quadratic form `p x² + 2 q x y + r y²`. -/
def quad2 (p q r x y : ℝ) : ℝ := p * x ^ 2 + 2 * q * x * y + r * y ^ 2

/-- A unit vector maximizing the quadratic form over the unit circle. -/
def IsRayleighMax (p q r x y : ℝ) : Prop :=
  x ^ 2 + y ^ 2 = 1 ∧ ∀ s t : ℝ, s ^ 2 + t ^ 2 = 1 → quad2 p q r s t ≤ quad2 p q r x y

/-- **SG-5, the Rayleigh step.** A maximizer of the Rayleigh quotient is an eigenvector: the form
`A·I − M` is nonnegative and vanishes at it. Stated through the cross product `w × Mw = 0`. -/
theorem rayleigh_max_is_eigen {p q r x y : ℝ} (h : IsRayleighMax p q r x y) :
    x * (q * x + r * y) - y * (p * x + q * y) = 0 := by
  obtain ⟨hunit, hmax⟩ := h
  set A := quad2 p q r x y with hA
  -- Homogeneity: every vector's form is at most `A` times its squared norm.
  have hall : ∀ s t : ℝ, quad2 p q r s t ≤ A * (s ^ 2 + t ^ 2) := by
    intro s t
    by_cases h0 : s ^ 2 + t ^ 2 = 0
    · have hs : s = 0 := by nlinarith [sq_nonneg s, sq_nonneg t]
      have ht : t = 0 := by nlinarith [sq_nonneg s, sq_nonneg t]
      simp [quad2, hs, ht]
    · have hpos : 0 < s ^ 2 + t ^ 2 := lt_of_le_of_ne (by positivity) (Ne.symm h0)
      set d := √(s ^ 2 + t ^ 2) with hd
      have hdpos : 0 < d := Real.sqrt_pos.2 hpos
      have hdsq : d ^ 2 = s ^ 2 + t ^ 2 := Real.sq_sqrt hpos.le
      have := hmax (s / d) (t / d) (by field_simp; linarith)
      have hq : quad2 p q r (s / d) (t / d) = quad2 p q r s t / d ^ 2 := by
        simp only [quad2]
        field_simp
      rw [hq, div_le_iff₀ (by positivity), hdsq] at this
      exact this
  -- The perpendicular direction, perturbed by `u`.
  set B := -y * (p * x + q * y) + x * (q * x + r * y) with hB
  set C := quad2 p q r (-y) x with hC
  have hline : ∀ u : ℝ, 2 * u * B + u ^ 2 * (C - A) ≤ 0 := by
    intro u
    have h := hall (x - u * y) (y + u * x)
    have hn : (x - u * y) ^ 2 + (y + u * x) ^ 2 = (1 + u ^ 2) * (x ^ 2 + y ^ 2) := by ring
    have hform : quad2 p q r (x - u * y) (y + u * x) = A + 2 * u * B + u ^ 2 * C := by
      simp only [hA, hB, hC, quad2]
      ring
    rw [hform, hn, hunit] at h
    nlinarith
  have hB0 : B = 0 := by
    have hc := hline (B / (|C - A| + 1))
    set D := |C - A| + 1 with hD
    have hDpos : 0 < D := by positivity
    have hle : C - A ≤ |C - A| := le_abs_self _
    have hexp : 2 * (B / D) * B + (B / D) ^ 2 * (C - A) = B ^ 2 * (2 * D + (C - A)) / D ^ 2 := by
      field_simp
    rw [hexp] at hc
    have hneg : 2 * D + (C - A) > 0 := by
      have := neg_abs_le (C - A)
      nlinarith
    by_contra hne
    have hpos : 0 < B ^ 2 * (2 * D + (C - A)) / D ^ 2 := by positivity
    linarith
  simp only [hB] at hB0
  linarith

/-- With `C̄ = 2aaᵀ + 100bbᵀ` and `P̄ = 2aaᵀ + bbᵀ` for `a = (1, 0)` and `b = (3/5, 4/5)`, whose
angle
is about 53°: `(38, 48, 64)` and `(59/25, 12/25, 16/25)`, the entries of `3 E[C]` and `3 E[P_j]`.
They share no eigenvector. -/
theorem no_common_eigenvector {x y : ℝ} (h : x ^ 2 + y ^ 2 = 1)
    (hC : x * (48 * x + 64 * y) - y * (38 * x + 48 * y) = 0)
    (hP : x * (12 / 25 * x + 16 / 25 * y) - y * (59 / 25 * x + 12 / 25 * y) = 0) : False := by
  have hxy : x * y = 0 := by nlinarith
  have hxx : x ^ 2 = y ^ 2 := by nlinarith
  rcases mul_eq_zero.1 hxy with hx | hy
  · rw [hx] at hxx h
    nlinarith
  · rw [hy] at hxx h
    nlinarith

/-- **SG-5, r = 1.** The top eigenvector of `E[C]` is not the top eigenvector of the draws' mean
projector, for draws `(λ, v) = (1, a), (1, a), (100, b)`: each form's Rayleigh maximizer differs
from the other's, up to sign. -/
theorem top_eigenspace_ne_subspace_mean {x y s t : ℝ}
    (hC : IsRayleighMax 38 48 64 x y) (hP : IsRayleighMax (59 / 25) (12 / 25) (16 / 25) s t) :
    ¬ ((x = s ∧ y = t) ∨ (x = -s ∧ y = -t)) := by
  have heC := rayleigh_max_is_eigen hC
  have heP := rayleigh_max_is_eigen hP
  have hunit := hC.1
  rintro (⟨rfl, rfl⟩ | ⟨rfl, rfl⟩)
  · exact no_common_eigenvector hunit heC heP
  · refine no_common_eigenvector hunit heC ?_
    linear_combination heP

end Rayleigh

end BASS
