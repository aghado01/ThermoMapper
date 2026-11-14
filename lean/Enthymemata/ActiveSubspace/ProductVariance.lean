import Mathlib.MeasureTheory.Constructions.Pi
import Mathlib.MeasureTheory.Integral.Prod
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.Ring

/-!
# Active subspaces — Efron–Stein on a product of probability measures (AS-2(d)(ii))

The continuous case of `BASS.efron_stein`, which AS-2's finite form reaches only on finite products.
* `sq_integral_le`: `(∫ φ)² ≤ ∫ φ²` under a probability measure.
* `var_prod_le`, the two-block step: under `ν ⊗ ρ`,
  `Var F ≤ ∫ Var_ρ F(y, ·) dν + ∫ Var_ν F(·, z) dρ` — the law of total variance, then Jensen on the
  inner average and Fubini.
* `efron_stein_pi`: under `Measure.pi μ`, `Var F ≤ Σᵢ E[Varᵢ F]` (`coordVarPi`), by induction
  over coordinates, splitting off coordinate `0` with `measurePreserving_piFinSuccAbove`
  (`integral_pi_succ`).

Functions are bounded and measurable, so every integral is finite.
-/

open MeasureTheory

namespace ActiveSubspace

/-- The variance under a measure. -/
noncomputable def varM {X : Type*} [MeasurableSpace X] (μ : Measure X) (f : X → ℝ) : ℝ :=
  ∫ x, (f x - ∫ y, f y ∂μ) ^ 2 ∂μ

section Basic

variable {X : Type*} [MeasurableSpace X]

/-- A bounded measurable function is integrable under a finite measure. -/
theorem integrable_of_bdd (μ : Measure X) [IsFiniteMeasure μ] {f : X → ℝ} (hf : Measurable f)
    {C : ℝ} (hC : ∀ x, |f x| ≤ C) : Integrable f μ :=
  Integrable.of_bound hf.aestronglyMeasurable C
    (ae_of_all _ fun x => by simpa [Real.norm_eq_abs] using hC x)

/-- `(∫ φ)² ≤ ∫ φ²` under a probability measure. -/
theorem sq_integral_le (μ : Measure X) [IsProbabilityMeasure μ] {φ : X → ℝ}
    (hφ : Integrable φ μ) (hφ2 : Integrable (fun x => φ x ^ 2) μ) :
    (∫ x, φ x ∂μ) ^ 2 ≤ ∫ x, φ x ^ 2 ∂μ := by
  set c := ∫ x, φ x ∂μ with hc
  have h0 : 0 ≤ ∫ x, (φ x - c) ^ 2 ∂μ := integral_nonneg fun x => sq_nonneg _
  have hexp : (fun x => (φ x - c) ^ 2) = fun x => φ x ^ 2 - 2 * c * φ x + c ^ 2 := by
    funext x
    ring
  have i1 : Integrable (fun x => φ x ^ 2 - 2 * c * φ x) μ := hφ2.sub (hφ.const_mul (2 * c))
  rw [hexp, integral_add i1 (integrable_const _), integral_sub hφ2 (hφ.const_mul (2 * c)),
    integral_const_mul, integral_const] at h0
  simp only [probReal_univ, one_smul, ← hc] at h0
  nlinarith [h0]

end Basic

section TwoBlock

variable {Y Z : Type*} [MeasurableSpace Y] [MeasurableSpace Z]

/-- **AS-2(d)(ii), the two-block step.** Under a product of probability measures, the variance is
at most the mean conditional variance in each block, summed:
`Var F ≤ ∫ Var_ρ F(y, ·) dν + ∫ Var_ν F(·, z) dρ`. -/
theorem var_prod_le (ν : Measure Y) (ρ : Measure Z) [IsProbabilityMeasure ν]
    [IsProbabilityMeasure ρ] {F : Y × Z → ℝ} (hF : Measurable F) {C : ℝ} (hC : ∀ p, |F p| ≤ C) :
    varM (ν.prod ρ) F ≤
      ∫ y, varM ρ (fun z => F (y, z)) ∂ν + ∫ z, varM ν (fun y => F (y, z)) ∂ρ := by
  -- the block means and the grand mean
  set g : Y → ℝ := fun y => ∫ z, F (y, z) ∂ρ with hgdef
  set k : Z → ℝ := fun z => ∫ y, F (y, z) ∂ν with hkdef
  set m : ℝ := ∫ p, F p ∂(ν.prod ρ) with hmdef
  have hg : Measurable g := (hF.stronglyMeasurable.integral_prod_right').measurable
  have hk : Measurable k := (hF.stronglyMeasurable.integral_prod_left').measurable
  have hgC : ∀ y, |g y| ≤ C := fun y => by
    refine (abs_integral_le_integral_abs).trans ?_
    calc ∫ z, |F (y, z)| ∂ρ ≤ ∫ _z, C ∂ρ :=
          integral_mono_of_nonneg (ae_of_all _ fun _ => abs_nonneg _) (integrable_const C)
            (ae_of_all _ fun z => hC (y, z))
      _ = C := by simp
  have hkC : ∀ z, |k z| ≤ C := fun z => by
    refine (abs_integral_le_integral_abs).trans ?_
    calc ∫ y, |F (y, z)| ∂ν ≤ ∫ _y, C ∂ν :=
          integral_mono_of_nonneg (ae_of_all _ fun _ => abs_nonneg _) (integrable_const C)
            (ae_of_all _ fun y => hC (y, z))
      _ = C := by simp
  have hFi : Integrable F (ν.prod ρ) := integrable_of_bdd _ hF hC
  have hmg : m = ∫ y, g y ∂ν := integral_prod F hFi
  have hmk : m = ∫ z, k z ∂ρ := integral_prod_symm F hFi
  -- bounds for the squared deviations
  have hsqC : ∀ (u v : ℝ), |u| ≤ C → |v| ≤ C → |(u - v) ^ 2| ≤ (2 * C) ^ 2 := by
    intro u v hu hv
    rw [abs_pow]
    have : |u - v| ≤ 2 * C := (abs_sub _ _).trans (by linarith)
    exact pow_le_pow_left₀ (abs_nonneg _) this 2
  have hmC : |m| ≤ C := by rw [hmg]; exact (abs_integral_le_integral_abs).trans (by
    calc ∫ y, |g y| ∂ν ≤ ∫ _y, C ∂ν :=
          integral_mono_of_nonneg (ae_of_all _ fun _ => abs_nonneg _) (integrable_const C)
            (ae_of_all _ fun y => hgC y)
      _ = C := by simp)
  -- Step A: the variance as an iterated integral
  have hA : varM (ν.prod ρ) F = ∫ y, ∫ z, (F (y, z) - m) ^ 2 ∂ρ ∂ν := by
    rw [varM, ← hmdef]
    exact integral_prod (fun p => (F p - m) ^ 2)
      (integrable_of_bdd _ ((hF.sub measurable_const).pow_const 2)
        (fun p => hsqC _ _ (hC p) hmC))
  -- Step B: for each `y`, the law of total variance
  have hB : ∀ y, ∫ z, (F (y, z) - m) ^ 2 ∂ρ = varM ρ (fun z => F (y, z)) + (g y - m) ^ 2 := by
    intro y
    have hFy : Measurable fun z => F (y, z) := hF.comp measurable_prodMk_left
    have hi1 : Integrable (fun z => (F (y, z) - g y) ^ 2) ρ :=
      integrable_of_bdd _ ((hFy.sub measurable_const).pow_const 2)
        (fun z => hsqC _ _ (hC (y, z)) (hgC y))
    have hi2 : Integrable (fun z => 2 * (g y - m) * (F (y, z) - g y)) ρ :=
      (integrable_of_bdd _ (hFy.sub measurable_const) (C := 2 * C) (fun z =>
        (abs_sub _ _).trans (by linarith [hC (y, z), hgC y]))).const_mul _
    have hexp : (fun z => (F (y, z) - m) ^ 2) = fun z =>
        (F (y, z) - g y) ^ 2 + 2 * (g y - m) * (F (y, z) - g y) + (g y - m) ^ 2 := by
      funext z
      ring
    have hmid : ∫ z, 2 * (g y - m) * (F (y, z) - g y) ∂ρ = 0 := by
      have iFy : Integrable (fun z => F (y, z)) ρ :=
        integrable_of_bdd _ hFy (fun z => hC (y, z))
      rw [integral_const_mul, integral_sub iFy (integrable_const _), integral_const]
      simp [hgdef]
    have hi12 : Integrable (fun z => (F (y, z) - g y) ^ 2 +
        2 * (g y - m) * (F (y, z) - g y)) ρ := hi1.add hi2
    rw [hexp, integral_add hi12 (integrable_const _), integral_add hi1 hi2, hmid,
      integral_const]
    simp [varM, hgdef]
  -- Step C: Jensen on the inner average, then Fubini
  have hdev : ∀ y, g y - m = ∫ z, (F (y, z) - k z) ∂ρ := by
    intro y
    have iF : Integrable (fun z => F (y, z)) ρ :=
      integrable_of_bdd _ (hF.comp measurable_prodMk_left) (fun z => hC (y, z))
    rw [integral_sub iF (integrable_of_bdd _ hk hkC), ← hmk]
  have hFk : Measurable fun p : Y × Z => (F p - k p.2) ^ 2 :=
    (hF.sub (hk.comp measurable_snd)).pow_const 2
  have hC' : ∀ y, (g y - m) ^ 2 ≤ ∫ z, (F (y, z) - k z) ^ 2 ∂ρ := by
    intro y
    rw [hdev y]
    exact sq_integral_le ρ
      (integrable_of_bdd _ ((hF.comp measurable_prodMk_left).sub hk) (C := 2 * C) (fun z =>
        (abs_sub _ _).trans (by linarith [hC (y, z), hkC z])))
      (integrable_of_bdd _ (hFk.comp measurable_prodMk_left)
        (fun z => hsqC _ _ (hC (y, z)) (hkC z)))
  have hswap : ∫ y, ∫ z, (F (y, z) - k z) ^ 2 ∂ρ ∂ν = ∫ z, varM ν (fun y => F (y, z)) ∂ρ := by
    rw [integral_integral_swap (integrable_of_bdd _ hFk (fun p => hsqC _ _ (hC p) (hkC p.2)))]
    rfl
  -- assemble
  have hvarρ : Integrable (fun y => varM ρ (fun z => F (y, z))) ν := by
    have hm' : Measurable fun y => varM ρ (fun z => F (y, z)) :=
      (((hF.sub (hg.comp measurable_fst)).pow_const 2).stronglyMeasurable
        |>.integral_prod_right').measurable
    refine integrable_of_bdd _ hm' (C := (2 * C) ^ 2) fun y => ?_
    simp only [varM]
    rw [abs_of_nonneg (integral_nonneg fun _ => sq_nonneg _)]
    calc ∫ z, (F (y, z) - ∫ z', F (y, z') ∂ρ) ^ 2 ∂ρ ≤ ∫ _z, (2 * C) ^ 2 ∂ρ :=
          integral_mono_of_nonneg (ae_of_all _ fun _ => sq_nonneg _) (integrable_const _)
            (ae_of_all _ fun z => (le_abs_self _).trans (hsqC _ _ (hC (y, z)) (hgC y)))
      _ = (2 * C) ^ 2 := by simp
  have hgm : Integrable (fun y => (g y - m) ^ 2) ν :=
    integrable_of_bdd _ ((hg.sub measurable_const).pow_const 2) (fun y => hsqC _ _ (hgC y) hmC)
  have hinner : Integrable (fun y => ∫ z, (F (y, z) - k z) ^ 2 ∂ρ) ν :=
    (integrable_of_bdd _ hFk (fun p => hsqC _ _ (hC p) (hkC p.2))).integral_prod_left
  rw [hA, integral_congr_ae (ae_of_all _ hB), integral_add hvarρ hgm, ← hswap]
  have := integral_mono hgm hinner hC'
  linarith

end TwoBlock

section Pi

/-- A bounded function's variance is at most `(2C)²`. -/
theorem varM_le {X : Type*} [MeasurableSpace X] (μ : Measure X) [IsProbabilityMeasure μ]
    {f : X → ℝ} {C : ℝ} (hC : ∀ x, |f x| ≤ C) : varM μ f ≤ (2 * C) ^ 2 := by
  have hm : |∫ y, f y ∂μ| ≤ C := (abs_integral_le_integral_abs).trans (by
    calc ∫ y, |f y| ∂μ ≤ ∫ _y, C ∂μ :=
          integral_mono_of_nonneg (ae_of_all _ fun _ => abs_nonneg _) (integrable_const C)
            (ae_of_all _ hC)
      _ = C := by simp)
  calc varM μ f ≤ ∫ _x, (2 * C) ^ 2 ∂μ :=
        integral_mono_of_nonneg (ae_of_all _ fun _ => sq_nonneg _) (integrable_const _)
          (ae_of_all _ fun x => by
            change (f x - ∫ y, f y ∂μ) ^ 2 ≤ (2 * C) ^ 2
            have h : |f x - ∫ y, f y ∂μ| ≤ 2 * C := (abs_sub _ _).trans (by linarith [hC x])
            have h2 := pow_le_pow_left₀ (abs_nonneg _) h 2
            rwa [sq_abs] at h2)
    _ = (2 * C) ^ 2 := by simp

theorem varM_nonneg {X : Type*} [MeasurableSpace X] (μ : Measure X) (f : X → ℝ) :
    0 ≤ varM μ f :=
  integral_nonneg fun _ => sq_nonneg _

theorem abs_varM_le {X : Type*} [MeasurableSpace X] (μ : Measure X) [IsProbabilityMeasure μ]
    {f : X → ℝ} {C : ℝ} (hC : ∀ x, |f x| ≤ C) : |varM μ f| ≤ (2 * C) ^ 2 := by
  rw [abs_of_nonneg (varM_nonneg μ f)]
  exact varM_le μ hC

/-- The variance is carried by a measure-preserving equivalence. -/
theorem varM_comp_equiv {X Y : Type*} [MeasurableSpace X] [MeasurableSpace Y] {μ : Measure X}
    {ν : Measure Y} (e : X ≃ᵐ Y) (he : MeasurePreserving e μ ν) (f : Y → ℝ) :
    varM μ (fun x => f (e x)) = varM ν f := by
  simp only [varM]
  rw [he.integral_comp' f]
  exact he.integral_comp' (fun y => (f y - ∫ y', f y' ∂ν) ^ 2)

/-- A section's variance is measurable in the section's parameter. -/
theorem measurable_varM_section {X Y : Type*} [MeasurableSpace X] [MeasurableSpace Y]
    (ν : Measure Y) [SFinite ν] {G : X × Y → ℝ} (hG : Measurable G) :
    Measurable fun x => varM ν (fun y => G (x, y)) := by
  have hm : Measurable fun x => ∫ y, G (x, y) ∂ν :=
    hG.stronglyMeasurable.integral_prod_right'.measurable
  exact ((hG.sub (hm.comp measurable_fst)).pow_const 2).stronglyMeasurable.integral_prod_right'
    |>.measurable

variable {α : Type*} [MeasurableSpace α]

/-- `E[Varᵢ F]`: the mean variance of `F` along coordinate `i`, the others held. -/
noncomputable def coordVarPi {d : ℕ} (μ : Fin d → Measure α) (i : Fin d)
    (F : (Fin d → α) → ℝ) : ℝ :=
  ∫ x, varM (μ i) (fun y => F (Function.update x i y)) ∂(Measure.pi μ)

theorem measurable_coordVar {d : ℕ} (μ : Fin d → Measure α) [∀ i, IsProbabilityMeasure (μ i)]
    (i : Fin d) {F : (Fin d → α) → ℝ} (hF : Measurable F) :
    Measurable fun x => varM (μ i) (fun y => F (Function.update x i y)) :=
  measurable_varM_section (μ i) (G := fun p : (Fin d → α) × α => F (Function.update p.1 i p.2))
    (hF.comp measurable_update')

/-- The product splits as coordinate `0` times the rest. -/
theorem integral_pi_succ {d : ℕ} (μ : Fin (d + 1) → Measure α) [∀ i, IsProbabilityMeasure (μ i)]
    {K : (Fin (d + 1) → α) → ℝ} (hK : Measurable K) {C : ℝ} (hKC : ∀ x, |K x| ≤ C) :
    ∫ x, K x ∂(Measure.pi μ) =
      ∫ y, ∫ z, K (Fin.cons y z) ∂(Measure.pi fun j : Fin d => μ j.succ) ∂(μ 0) := by
  set e := MeasurableEquiv.piFinSuccAbove (fun _ : Fin (d + 1) => α) 0
  have hmp : MeasurePreserving e (Measure.pi μ) ((μ 0).prod (Measure.pi fun j => μ j.succ)) :=
    measurePreserving_piFinSuccAbove μ 0
  have hsymm : ∀ y z, e.symm (y, z) = Fin.cons y z := fun y z =>
    (rfl : e.symm (y, z) = Fin.insertNth (α := fun _ => α) 0 y z).trans (Fin.insertNth_zero' y z)
  rw [← (hmp.symm e).integral_comp' K, integral_prod (fun p => K (e.symm p))
    (integrable_of_bdd (f := fun p => K (e.symm p)) _ (hK.comp e.symm.measurable)
      (fun p => hKC _))]
  simp only [hsymm]

/-- **AS-2(d)(ii), Efron–Stein on a product of probability measures.** For bounded measurable
`F`, `Var F ≤ Σᵢ E[Varᵢ F]`. -/
theorem efron_stein_pi (d : ℕ) : ∀ (μ : Fin d → Measure α) [∀ i, IsProbabilityMeasure (μ i)]
    (F : (Fin d → α) → ℝ), Measurable F → ∀ C : ℝ, (∀ x, |F x| ≤ C) →
      varM (Measure.pi μ) F ≤ ∑ i, coordVarPi μ i F := by
  induction d with
  | zero =>
    intro μ _ F _ C _
    rw [Finset.univ_eq_empty, Finset.sum_empty]
    have hconst : ∀ x : Fin 0 → α, F x = F isEmptyElim := fun x =>
      congrArg F (Subsingleton.elim _ _)
    simp only [varM, hconst, integral_const, probReal_univ, one_smul, sub_self]
    simp
  | succ d ih =>
    intro μ _ F hF C hC
    set μ' : Fin d → Measure α := fun j => μ j.succ with hμ'
    set e := MeasurableEquiv.piFinSuccAbove (fun _ : Fin (d + 1) => α) 0
    have hmp : MeasurePreserving e (Measure.pi μ) ((μ 0).prod (Measure.pi μ')) :=
      measurePreserving_piFinSuccAbove μ 0
    set G : α × (Fin d → α) → ℝ := fun p => F (e.symm p) with hGdef
    have hG : Measurable G := hF.comp e.symm.measurable
    have hGC : ∀ p, |G p| ≤ C := fun p => hC _
    have hGcons : ∀ y z, G (y, z) = F (Fin.cons y z) := fun y z =>
      congrArg F ((rfl : e.symm (y, z) = Fin.insertNth (α := fun _ => α) 0 y z).trans
        (Fin.insertNth_zero' y z))
    have hvar : varM (Measure.pi μ) F = varM ((μ 0).prod (Measure.pi μ')) G := by
      have h := varM_comp_equiv e hmp G
      have hfe : (fun x => G (e x)) = F := funext fun x => by simp [hGdef]
      rwa [hfe] at h
    -- the block of coordinate `0`
    have hblock0 : coordVarPi μ 0 F =
        ∫ z, varM (μ 0) (fun y => G (y, z)) ∂(Measure.pi μ') := by
      rw [coordVarPi, integral_pi_succ μ (measurable_coordVar μ 0 hF)
        (fun x => abs_varM_le (μ 0) (fun y => hC _))]
      have hpt : ∀ y z, varM (μ 0) (fun t => F (Function.update (Fin.cons y z) 0 t)) =
          varM (μ 0) (fun t => G (t, z)) := fun y z => by
        simp only [Fin.update_cons_zero, hGcons]
      simp only [hpt, integral_const, probReal_univ, one_smul]
      rfl
    -- the other blocks, through the induction hypothesis
    have hblockj : ∀ j : Fin d, ∫ y, coordVarPi μ' j (fun z => G (y, z)) ∂(μ 0) =
        coordVarPi μ j.succ F := by
      intro j
      conv_rhs => rw [coordVarPi, integral_pi_succ μ (measurable_coordVar μ j.succ hF)
        (fun x => abs_varM_le _ (fun y => hC _))]
      refine integral_congr_ae (ae_of_all _ fun y => ?_)
      simp only [coordVarPi]
      refine integral_congr_ae (ae_of_all _ fun z => ?_)
      simp only [← Fin.cons_update, hGcons]
      rfl
    have hmeasj : ∀ j : Fin d, Measurable fun y => coordVarPi μ' j (fun z => G (y, z)) := by
      intro j
      have hΦ : Measurable fun p : α × (Fin d → α) =>
          varM (μ' j) (fun t => G (p.1, Function.update p.2 j t)) :=
        measurable_varM_section (μ' j) (G := fun q : (α × (Fin d → α)) × α =>
          G (q.1.1, Function.update q.1.2 j q.2))
          (hG.comp ((measurable_fst.comp measurable_fst).prodMk
            (measurable_update'.comp ((measurable_snd.comp measurable_fst).prodMk
              measurable_snd))))
      exact hΦ.stronglyMeasurable.integral_prod_right'.measurable
    have hbddj : ∀ j : Fin d, ∀ y, |coordVarPi μ' j (fun z => G (y, z))| ≤ (2 * C) ^ 2 := by
      intro j y
      refine (abs_integral_le_integral_abs).trans ?_
      calc ∫ z, |varM (μ' j) (fun t => G (y, Function.update z j t))| ∂(Measure.pi μ')
          ≤ ∫ _z, (2 * C) ^ 2 ∂(Measure.pi μ') :=
            integral_mono_of_nonneg (ae_of_all _ fun _ => abs_nonneg _) (integrable_const _)
              (ae_of_all _ fun z => abs_varM_le _ (fun t => hGC _))
        _ = (2 * C) ^ 2 := by simp
    have hint1 : Integrable (fun y => varM (Measure.pi μ') (fun z => G (y, z))) (μ 0) :=
      integrable_of_bdd _ (measurable_varM_section _ hG)
        (fun y => abs_varM_le _ (fun z => hGC _))
    have hint2 : Integrable (fun y => ∑ j, coordVarPi μ' j (fun z => G (y, z))) (μ 0) :=
      integrable_finsetSum _ fun j _ => integrable_of_bdd _ (hmeasj j) (hbddj j)
    have hIH : ∫ y, varM (Measure.pi μ') (fun z => G (y, z)) ∂(μ 0) ≤
        ∑ j : Fin d, coordVarPi μ j.succ F := by
      calc ∫ y, varM (Measure.pi μ') (fun z => G (y, z)) ∂(μ 0)
          ≤ ∫ y, ∑ j, coordVarPi μ' j (fun z => G (y, z)) ∂(μ 0) :=
            integral_mono hint1 hint2 fun y =>
              ih μ' (fun z => G (y, z)) (hG.comp measurable_prodMk_left) C (fun z => hGC _)
        _ = ∑ j : Fin d, coordVarPi μ j.succ F := by
            rw [integral_finsetSum _ fun j _ => integrable_of_bdd _ (hmeasj j) (hbddj j)]
            exact Finset.sum_congr rfl fun j _ => hblockj j
    rw [hvar, Fin.sum_univ_succ, hblock0]
    have h2 := var_prod_le (μ 0) (Measure.pi μ') hG hGC
    linarith

end Pi

end ActiveSubspace
