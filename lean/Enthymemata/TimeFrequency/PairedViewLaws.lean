import Enthymemata.TimeFrequency.FinitePowerMeasurements
import Enthymemata.TimeFrequency.TwoSampleRecovery
import Mathlib.MeasureTheory.Measure.Dirac
import Mathlib.Topology.Instances.Matrix
import Mathlib.MeasureTheory.Constructions.BorelSpace.Basic

/-! # Identification from the labeled views of the same real waveform -/

namespace TimeFrequency.PairedViews

open MeasureTheory Set FinitePower
open scoped ENNReal

variable {T J : Type*} [Fintype T]

/-- The compact outer-product image represents the waveform class modulo its real global sign. -/
def outerClass (C : Set (T → ℝ)) := outer '' C

omit [Fintype T] in
theorem outerClass_compact {C : Set (T → ℝ)} (hc : IsCompact C) : IsCompact (outerClass C) :=
  hc.image (show Continuous (outer (T := T)) by
    change Continuous (fun x : T → ℝ ↦ fun i j ↦ x i * x j)
    fun_prop)

theorem lifted_continuous (a : J → T → ℂ) : Continuous (lifted a) := by
  change Continuous (fun H : T → T → ℝ ↦ fun j ↦ ∑ t, ∑ s, cellMatrix (a j) t s * H t s)
  fun_prop

theorem lifted_injective_iff (a : J → T → ℂ) (C : Set (T → ℝ)) :
    Function.Injective (fun H : outerClass C ↦ lifted a H.val) ↔
    ∀ x ∈ C, ∀ y ∈ C, powers a x = powers a y → y = x ∨ y = -x := by
  constructor
  · intro hi x hx y hy he
    apply (outer_eq_iff x y).mp
    apply congrArg Subtype.val (hi (a₁ := ⟨outer x, ⟨x,hx,rfl⟩⟩) (a₂ := ⟨outer y,⟨y,hy,rfl⟩⟩) ?_)
    simpa only [← powers_lifted] using he
  · intro hs U V he
    obtain ⟨x,hx,hX⟩ := U.property
    obtain ⟨y,hy,hY⟩ := V.property
    apply Subtype.ext
    rw [← hX,← hY]
    apply (outer_eq_iff _ _).mpr (hs _ hx _ hy ?_)
    rwa [powers_lifted, powers_lifted, hX, hY]

section CompactModel

variable {K Y : Type*} [TopologicalSpace K] [CompactSpace K] [TopologicalSpace Y] [T2Space Y]
  {M : K → Y}

/-- The inverse is continuous on precisely the observed image. -/
noncomputable def imageHomeomorph (hc : Continuous M) (hi : Function.Injective M) :
    K ≃ₜ Set.range M := (hc.isClosedEmbedding hi).toIsEmbedding.toHomeomorph

theorem inverse_recovers (hc : Continuous M) (hi : Function.Injective M) (x : K) :
    (imageHomeomorph hc hi).symm ⟨M x, ⟨x,rfl⟩⟩ = x :=
  (imageHomeomorph hc hi).symm_apply_apply x

omit [TopologicalSpace K] [CompactSpace K] [TopologicalSpace Y] [T2Space Y] in
theorem image_sets_eq_iff (hi : Function.Injective M) (A B : Set K) :
    M '' A = M '' B ↔ A = B := hi.image_injective.eq_iff

theorem paired_laws_eq_iff [MeasurableSpace K] [BorelSpace K] [MeasurableSpace Y] [BorelSpace Y]
    (hc : Continuous M) (hi : Function.Injective M) (μ ν : Measure K) :
    μ.map M = ν.map M ↔ μ = ν := (hc.isClosedEmbedding hi).measurableEmbedding.map_injective.eq_iff

end CompactModel

theorem collision_diracs {K Y : Type*} [MeasurableSpace K] [MeasurableSingletonClass K]
    [MeasurableSpace Y] {M : K → Y} (hm : Measurable M) {x y : K} (hxy : x ≠ y) (he : M x = M y) :
    Measure.dirac x ≠ Measure.dirac y ∧ (Measure.dirac x).map M = (Measure.dirac y).map M :=
  ⟨dirac_ne_dirac hxy, by simp [Measure.map_dirac' hm,he]⟩

theorem signed_diracs [Countable J] (a : J → T → ℂ) {x : T → ℝ} (hx : x ≠ 0) :
    Measure.dirac x ≠ Measure.dirac (-x) ∧
    (Measure.dirac x).map (powers a) = (Measure.dirac (-x)).map (powers a) := by
  apply collision_diracs
  · have he : powers a = lifted a ∘ outer := funext (powers_lifted a)
    rw [he]
    apply Continuous.measurable
    apply (lifted_continuous a).comp
    change Continuous (fun x : T → ℝ ↦ fun i j ↦ x i * x j)
    fun_prop
  · intro h
    apply hx
    funext t
    have he := congrFun h t
    change x t = -x t at he
    change x t = 0
    linarith
  · simp_rw [powers_lifted]
    rw [(outer_eq_iff x (-x)).mpr (Or.inr rfl)]

abbrev Wave := ℝ × ℝ

def singletonView (x : Wave) : Wave := (|x.1|, |x.2|)
def frameView (x : Wave) : Wave := (|x.1+x.2|, |x.1-x.2|)
def jointView (x : Wave) := (singletonView x, frameView x)

/-- Squaring the four magnitudes gives the raw singleton, DC and Nyquist Fourier powers. -/
theorem views_squared (x : Wave) :
    (singletonView x).1 ^ 2 = framePower ![x.1] 0 ∧
    (singletonView x).2 ^ 2 = framePower ![x.2] 0 ∧
    (frameView x).1 ^ 2 = framePower ![x.1,x.2] 0 ∧
    (frameView x).2 ^ 2 = framePower ![x.1,x.2] 1 := by
  simp only [singletonView, frameView, sq_abs]
  refine ⟨?_, ?_, (twoSample_powers _ _).1.symm, (twoSample_powers _ _).2.symm⟩ <;>
    simp [framePower_dc]

noncomputable def twoPoint {K : Type*} [MeasurableSpace K] (x y : K) : Measure K :=
  (2 : ℝ≥0∞)⁻¹ • (Measure.dirac x + Measure.dirac y)

instance {K : Type*} [MeasurableSpace K] (x y : K) : IsProbabilityMeasure (twoPoint x y) := by
  constructor
  simp [twoPoint, smul_eq_mul, ENNReal.inv_two_add_inv_two]

theorem twoPoint_map {K Y : Type*} [MeasurableSpace K] [MeasurableSpace Y]
    (M : K → Y) (hm : Measurable M) (x y : K) :
    (twoPoint x y).map M = twoPoint (M x) (M y) := by
  simp [twoPoint, Measure.map_smul, Measure.map_add _ _ hm, Measure.map_dirac' hm]

noncomputable abbrev modelX := twoPoint ((1,2) : Wave) (2,-1)
noncomputable abbrev modelY := twoPoint ((1,-2) : Wave) (2,1)

theorem four_tuples :
    jointView (1,2) = ((1,2),(3,1)) ∧ jointView (2,-1) = ((2,1),(1,3)) ∧
    jointView (1,-2) = ((1,2),(1,3)) ∧ jointView (2,1) = ((2,1),(3,1)) := by
  norm_num [jointView,singletonView,frameView]

theorem equal_single_views :
    modelX.map singletonView = modelY.map singletonView ∧
      modelX.map frameView = modelY.map frameView := by
  have h₁ : Measurable singletonView := by unfold singletonView; fun_prop
  have h₂ : Measurable frameView := by unfold frameView; fun_prop
  constructor
  · rw [twoPoint_map _ h₁, twoPoint_map _ h₁]
    norm_num [singletonView]
  · rw [twoPoint_map _ h₂, twoPoint_map _ h₂]
    norm_num [frameView, twoPoint, add_comm]

theorem different_joint_laws :
    let event : Set (Wave × Wave) := {((1,2),(3,1))}
    modelX.map jointView event = (2 : ℝ≥0∞)⁻¹ ∧ modelY.map jointView event = 0 := by
  have hm : Measurable jointView := by unfold jointView singletonView frameView; fun_prop
  rw [twoPoint_map _ hm, twoPoint_map _ hm, four_tuples.1, four_tuples.2.1,
    four_tuples.2.2.1, four_tuples.2.2.2]
  norm_num [twoPoint, Measure.smul_apply, Measure.add_apply, Measure.dirac_apply,
    Set.indicator, smul_eq_mul]

theorem different_quotient_supports :
    ((1,2) : Wave) ∉ ({(1,-2),(2,1),-(1,-2),-(2,1)} : Set Wave) := by
  norm_num [Prod.ext_iff]

/-- Both laws live on the same compact outer-product class with its Borel structure. -/
theorem paired_power_laws [Countable J] (a : J → T → ℂ) {C : Set (T → ℝ)} (hc : IsCompact C)
    [MeasurableSpace (outerClass C)] [BorelSpace (outerClass C)]
    (hs : ∀ x ∈ C, ∀ y ∈ C, powers a x = powers a y → y = x ∨ y = -x)
    (μ ν : Measure (outerClass C)) :
    μ.map (fun H : outerClass C ↦ lifted a H.val) = ν.map (fun H : outerClass C ↦ lifted a H.val) ↔
      μ = ν := by
  letI : CompactSpace (outerClass C) := isCompact_iff_compactSpace.mp (outerClass_compact hc)
  exact paired_laws_eq_iff ((lifted_continuous a).comp continuous_subtype_val)
    ((lifted_injective_iff a C).mpr hs) μ ν

#print axioms paired_power_laws
#print axioms outerClass_compact
#print axioms lifted_continuous
#print axioms lifted_injective_iff
#print axioms inverse_recovers
#print axioms image_sets_eq_iff
#print axioms paired_laws_eq_iff
#print axioms collision_diracs
#print axioms signed_diracs
#print axioms twoPoint_map
#print axioms views_squared
#print axioms four_tuples
#print axioms equal_single_views
#print axioms different_joint_laws
#print axioms different_quotient_supports

end TimeFrequency.PairedViews
