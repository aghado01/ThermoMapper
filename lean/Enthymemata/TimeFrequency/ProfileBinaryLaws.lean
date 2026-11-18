import Enthymemata.TimeFrequency.ProfilePairMoments
import Mathlib.MeasureTheory.Measure.Dirac

/-! # Two attainable harmonic fractions realizing complete pair laws -/

namespace TimeFrequency.LinkedProfiles
open MeasureTheory Set
open scoped ENNReal NNReal
noncomputable section

variable {V : Type*}

def fraction (l h : ℝ) (hit : Bool) := if hit then h else l
def success (l h θ : ℝ) := (θ - l) / (h - l)
def bitMass (p : ℝ) (hit : Bool) := if hit then p else 1 - p
def sign (s : Bool) : ℝ := if s then 1 else -1

def bitLaw (p : ℝ) : Measure Bool :=
  p.toNNReal • Measure.dirac true + (1 - p).toNNReal • Measure.dirac false

instance (p : ℝ) : IsFiniteMeasure (bitLaw p) := by unfold bitLaw; infer_instance

theorem success_mem {l h θ : ℝ} (hlh : l < h) (hθ : θ ∈ Icc l h) :
    success l h θ ∈ Icc 0 1 :=
  ⟨div_nonneg (sub_nonneg.mpr hθ.1) (sub_nonneg.mpr hlh.le),
    (div_le_one (sub_pos.mpr hlh)).mpr (by linarith [hθ.2])⟩

theorem bitLaw_probability {p : ℝ} (hp : p ∈ Icc 0 1) : IsProbabilityMeasure (bitLaw p) := by
  constructor
  simp [bitLaw, ← ENNReal.coe_add,
    ← Real.toNNReal_add hp.1 (sub_nonneg.mpr hp.2)]

theorem bitLaw_mass {p : ℝ} (hp : p ∈ Icc 0 1) (i : Bool) :
    (bitLaw p).real {i} = bitMass p i := by
  cases i <;> simp [bitLaw, bitMass, measureReal_def,
    Real.toNNReal_of_nonneg hp.1, Real.toNNReal_of_nonneg (sub_nonneg.mpr hp.2)]

theorem fraction_mean {l h θ : ℝ} (hlh : l < h) (hθ : θ ∈ Icc l h) :
    ∫ b, fraction l h b ∂bitLaw (success l h θ) = θ := by
  have hp := success_mem hlh hθ
  rw [bitLaw, integral_add_measure]
  · calc
      _ = success l h θ * h + (1 - success l h θ) * l := by
        simp [integral_smul_nnreal_measure, fraction, Real.toNNReal_of_nonneg hp.1,
          Real.toNNReal_of_nonneg (sub_nonneg.mpr hp.2), NNReal.smul_def, smul_eq_mul]
      _ = θ := by
        unfold success
        field_simp [sub_ne_zero.mpr hlh.ne']
        ring
  all_goals apply Integrable.smul_measure <;> simp

def hitProb (l h : ℝ) (μ a : V → ℝ) (s : Bool) (v : V) :=
  success l h (μ v + sign s * a v)

/-- The same fair source is used for both conditionally independent emissions. -/
def bitPairs (l h : ℝ) (μ a : V → ℝ) (u v : V) : Measure (Bool × Bool) :=
  balanced ((bitLaw (hitProb l h μ a true u)).prod (bitLaw (hitProb l h μ a true v)))
    ((bitLaw (hitProb l h μ a false u)).prod (bitLaw (hitProb l h μ a false v)))

instance {X : Type*} [MeasurableSpace X] (p q : Measure X) [IsFiniteMeasure p] [IsFiniteMeasure q] :
    IsFiniteMeasure (balanced p q) := ⟨by
  simp only [balanced, Measure.smul_apply, smul_eq_mul]
  finiteness⟩

instance (l h : ℝ) (μ a : V → ℝ) (u v : V) : IsFiniteMeasure (bitPairs l h μ a u v) := by
  unfold bitPairs
  infer_instance

theorem bitPairs_cell {l h : ℝ} (hlh : l < h) (μ a : V → ℝ)
    (hadm : ∀ s v, μ v + sign s * a v ∈ Icc l h) (u v : V) (i j : Bool) :
    (bitPairs l h μ a u v).real {(i,j)} =
      bitMass (success l h (μ u)) i * bitMass (success l h (μ v)) j +
        sign i * sign j * (a u * a v) / (h - l)^2 := by
  have hp (s w) : hitProb l h μ a s w ∈ Icc 0 1 := success_mem hlh (hadm s w)
  unfold bitPairs balanced
  rw [measureReal_ennreal_smul_apply, measureReal_add_apply, ← singleton_prod_singleton,
    measureReal_prod_prod, measureReal_prod_prod]
  simp only [bitLaw_mass (hp _ _)]
  cases i <;> cases j <;> norm_num [hitProb, success, bitMass, sign] <;>
    field_simp [sub_ne_zero.mpr hlh.ne'] <;> ring

/-- Marginal means and one contrast product determine the entire two-by-two table. -/
theorem bitPairs_eq {l h : ℝ} (hlh : l < h) (μ a b : V → ℝ)
    (ha : ∀ s v, μ v + sign s * a v ∈ Icc l h)
    (hb : ∀ s v, μ v + sign s * b v ∈ Icc l h)
    {u v : V} (he : a u * a v = b u * b v) :
    bitPairs l h μ a u v = bitPairs l h μ b u v := by
  apply MeasureTheory.ext_iff_measureReal_singleton.mpr
  rintro ⟨i,j⟩
  rw [bitPairs_cell hlh μ a ha, bitPairs_cell hlh μ b hb, he]

/-- Push the full low/high draw through any fixed synthesis or measurement. -/
def readPairs {X : Type*} [MeasurableSpace X] (readout : Bool → X) (law : Measure (Bool × Bool)) :=
  law.map (Prod.map readout readout)

theorem readPairs_cell {X : Type*} [MeasurableSpace X] [MeasurableSingletonClass X]
    (f : Bool → X) (hf : Function.Injective f) (law : Measure (Bool × Bool)) (i j : Bool) :
    (readPairs f law).real {(f i,f j)} = law.real {(i,j)} := by
  have hm : Measurable (Prod.map f f) := Measurable.of_discrete
  have hi : Function.Injective (Prod.map f f) := hf.prodMap hf
  unfold readPairs
  rw [measureReal_def, Measure.map_apply hm (measurableSet_singleton _)]
  exact congrArg (fun S ↦ (law S).toReal) (hi.preimage_image (s := {(i,j)}) ▸ by simp)

/-- Population model obtained by mapping each Bernoulli draw to its harmonic fraction. -/
def binaryModel {l h : ℝ} (hlh : l < h) (μ a : V → ℝ)
    (hadm : ∀ s v, μ v + sign s * a v ∈ Icc l h) : ResponseModel V where
  law s v := (bitLaw (hitProb l h μ a s v)).map (fraction l h)
  probability s v := by
    letI : IsProbabilityMeasure (bitLaw (hitProb l h μ a s v)) :=
      bitLaw_probability (success_mem hlh (hadm s v))
    exact Measure.isProbabilityMeasure_map (Measurable.of_discrete.aemeasurable)
  integrable s v := by
    apply (integrable_map_measure measurable_id.aestronglyMeasurable
      Measurable.of_discrete.aemeasurable).mpr
    unfold bitLaw
    apply Integrable.add_measure <;> apply Integrable.smul_measure_nnreal <;> simp

theorem binaryModel_mean {l h : ℝ} (hlh : l < h) (μ a : V → ℝ)
    (hadm : ∀ s v, μ v + sign s * a v ∈ Icc l h) (s : Bool) (v : V) :
    (binaryModel hlh μ a hadm).mean s v = μ v + sign s * a v := by
  change (∫ t, id t ∂(bitLaw (hitProb l h μ a s v)).map (fraction l h)) = _
  rw [integral_map Measurable.of_discrete.aemeasurable measurable_id.aestronglyMeasurable]
  exact fraction_mean hlh (hadm s v)

theorem binaryModel_pairs {l h : ℝ} (hlh : l < h) (μ a : V → ℝ)
    (hadm : ∀ s v, μ v + sign s * a v ∈ Icc l h) (u v : V) :
    (binaryModel hlh μ a hadm).pairLaw u v = readPairs (fraction l h) (bitPairs l h μ a u v) := by
  have hmeas : Measurable (fraction l h) := Measurable.of_discrete
  unfold ResponseModel.pairLaw binaryModel readPairs bitPairs balanced
  rw [Measure.map_smul, Measure.map_add _ _ (hmeas.prodMap hmeas)]
  simp only [Measure.map_prod_map _ _ hmeas hmeas]

#print axioms binaryModel_mean
#print axioms binaryModel_pairs

#print axioms success_mem
#print axioms bitLaw_probability
#print axioms bitLaw_mass
#print axioms fraction_mean
#print axioms bitPairs_cell
#print axioms bitPairs_eq
#print axioms readPairs_cell

end
end TimeFrequency.LinkedProfiles
