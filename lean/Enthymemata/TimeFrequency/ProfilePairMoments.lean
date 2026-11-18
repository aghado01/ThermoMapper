import Enthymemata.TimeFrequency.ProfilePairGraph
import Mathlib.MeasureTheory.Integral.Prod

/-! # Balanced same-source pair laws and recovery of context-specific means -/

namespace TimeFrequency.LinkedProfiles
open MeasureTheory Set
open scoped ENNReal
noncomputable section

variable {V : Type*}

/-- Each source has its own response distribution in each observed context. -/
structure ResponseModel (V : Type*) where
  law : Bool → V → Measure ℝ
  probability : ∀ s v, IsProbabilityMeasure (law s v)
  integrable : ∀ s v, Integrable (fun t : ℝ ↦ t) (law s v)

attribute [instance] ResponseModel.probability

def balanced {X : Type*} [MeasurableSpace X] (p q : Measure X) : Measure X :=
  (2 : ℝ≥0∞)⁻¹ • (p + q)

instance balanced_probability {X : Type*} [MeasurableSpace X] (p q : Measure X)
    [IsProbabilityMeasure p] [IsProbabilityMeasure q] : IsProbabilityMeasure (balanced p q) := by
  constructor
  simp only [balanced, Measure.smul_apply, Measure.add_apply, measure_univ, smul_eq_mul,
    one_add_one_eq_two]
  exact ENNReal.inv_mul_cancel (by norm_num) (by norm_num)

theorem balanced_integral {X : Type*} [MeasurableSpace X] {p q : Measure X} {f : X → ℝ}
    (hp : Integrable f p) (hq : Integrable f q) :
    ∫ x, f x ∂balanced p q = ((∫ x, f x ∂p) + ∫ x, f x ∂q) / 2 := by
  rw [balanced, integral_smul_measure, integral_add_measure hp hq]
  norm_num
  ring

namespace ResponseModel
variable (P : ResponseModel V)

def mean (s : Bool) (v : V) := ∫ t, t ∂P.law s v
def center (v : V) := (P.mean true v + P.mean false v) / 2
def contrast (v : V) := (P.mean true v - P.mean false v) / 2

/-- Product measures express conditional independence, also when the contexts coincide. -/
def pairLaw (u v : V) := balanced ((P.law true u).prod (P.law true v))
  ((P.law false u).prod (P.law false v))

theorem paired_first (u v : V) : (∫ z, z.1 ∂P.pairLaw u v) = P.center u := by
  unfold pairLaw
  rw [balanced_integral ((P.integrable true u).comp_fst _) ((P.integrable false u).comp_fst _)]
  simp [integral_fun_fst (fun t : ℝ ↦ t), center, mean]

theorem paired_second (u v : V) : (∫ z, z.2 ∂P.pairLaw u v) = P.center v := by
  unfold pairLaw
  rw [balanced_integral ((P.integrable true v).comp_snd _) ((P.integrable false v).comp_snd _)]
  simp [integral_fun_snd (fun t : ℝ ↦ t), center, mean]

theorem paired_centered (u v : V) :
    ∫ z, (z.1 - P.center u) * (z.2 - P.center v) ∂P.pairLaw u v = P.contrast u * P.contrast v := by
  have hi (s w) : Integrable (fun t ↦ t - P.center w) (P.law s w) :=
    (P.integrable s w).sub (integrable_const _)
  have hc (s w) : (∫ t, t - P.center w ∂P.law s w) = P.mean s w - P.center w := by
    rw [integral_sub (P.integrable s w) (integrable_const _)]
    simp [mean]
  unfold pairLaw
  rw [balanced_integral ((hi true u).mul_prod (hi true v)) ((hi false u).mul_prod (hi false v))]
  simp only [integral_prod_mul (fun t ↦ t - P.center u) (fun t ↦ t - P.center v), hc]
  unfold center contrast
  ring

/-- Population variance conditional on the source, written in raw moments. -/
def withinVariance (s : Bool) (u : V) := (∫ t, t^2 ∂P.law s u) - (P.mean s u)^2

theorem single_square_excess (u : V) (h₂ : ∀ s, Integrable (fun t : ℝ ↦ t^2) (P.law s u)) :
    (∫ t, t^2 ∂balanced (P.law true u) (P.law false u)) - (P.center u)^2 =
      (P.contrast u)^2 + (P.withinVariance true u + P.withinVariance false u) / 2 := by
  rw [balanced_integral (h₂ true) (h₂ false)]
  unfold withinVariance center contrast
  ring

/-- On the sampled edges, equality of complete pair laws gives equality of both profiles up
to a single exchange of the source names. A root repeat or an odd closed walk anchors scale. -/
theorem identify_profiles (Q : ResponseModel V) (G : SimpleGraph V) (hG : G.Preconnected)
    (sampled : V → V → Prop) (covered : ∀ u, ∃ v, sampled u v)
    (edges : ∀ u v, G.Adj u v → sampled u v)
    (same : ∀ u v, sampled u v → Q.pairLaw u v = P.pairLaw u v)
    (ha : ∀ v, P.contrast v ≠ 0) (root : V)
    (anchor : sampled root root ∨ ∃ w : G.Walk root root, Odd w.length) :
    Q.mean = P.mean ∨ Q.mean = fun s ↦ P.mean (!s) := by
  have hμ (u : V) : Q.center u = P.center u := by
    obtain ⟨v,h⟩ := covered u
    rw [← Q.paired_first u v, same u v h, P.paired_first]
  have hm (u v : V) (h : sampled u v) :
      Q.contrast u * Q.contrast v = P.contrast u * P.contrast v := by
    rw [← Q.paired_centered, same u v h, hμ, hμ, P.paired_centered]
  have he := fun u v h ↦ hm u v (edges u v h)
  have ident : Q.contrast = P.contrast ∨ Q.contrast = -P.contrast := by
    rcases anchor with h | ⟨w,hw⟩
    · exact identify_from_square hG ha he root (by simpa [sq] using hm root root h)
    · exact identify_from_odd_walk hG ha he w hw
  rcases ident with plus | minus
  · left
    funext s v
    have h := congrFun plus v
    have k := hμ v
    cases s <;> unfold center at k <;> unfold contrast at h <;> linarith
  · right
    funext s v
    have h := congrFun minus v
    have k := hμ v
    simp only [Pi.neg_apply, contrast] at h
    cases s <;> simp only [Bool.not_true, Bool.not_false] <;> unfold center at k <;> linarith

end ResponseModel

#print axioms balanced_probability
#print axioms balanced_integral
#print axioms ResponseModel.paired_first
#print axioms ResponseModel.paired_second
#print axioms ResponseModel.paired_centered
#print axioms ResponseModel.single_square_excess
#print axioms ResponseModel.identify_profiles

end
end TimeFrequency.LinkedProfiles
