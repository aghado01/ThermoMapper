/-
  Stability.lean — the confidence pushforward (capstone of the BARS → barcode chain).

  Companion to:
    • Bifiltration.lean — the estimand is monotone (spin_agreement_mono ⇒ functorial filtration).
    • BARS.lean         — the per-draw reduction is exact (argmax / span are closed-form roots).
  This file joins them: a sup-norm credible band on the curve ⇒ a bottleneck confidence ball on
  the diagram, with no inflation.

  AXIOM  : diagram_stability   (function → diagram is 1-Lipschitz; REF-PH §5.2, BCKL2010 §2.2/§4).
  PROVED : stability_event_subset, confidence_pushforward   (measure monotonicity over the axiom).
  OBLIGATION : landscape_stable   (FH2024 §4 — the averageable surrogate; build with TDA.Ph).

  Proto-lemma: ../proto-lemmas/confidence-pushforward-lemmas.md
-/
import Mathlib.MeasureTheory.Measure.MeasureSpace
import Mathlib.Topology.MetricSpace.Basic

namespace Spc.Stability

/-! ### 0. Opaque stand-ins for the TDA-layer objects (as Bifiltration.lean does for the measure). -/

opaque PersistenceDiagram : Type
opaque bottleneck : PersistenceDiagram → PersistenceDiagram → ℝ
/-- Sup-norm of a curve perturbation, `‖f − g‖_∞`. -/
opaque supNorm : (ℝ → ℝ) → ℝ
/-- The filtration-then-diagram functor: a curve (filtration index → value) to its barcode. The
    monotonicity that makes this well-defined is `Bifiltration.lean`'s `spin_agreement_mono`. -/
opaque Dgm : (ℝ → ℝ) → PersistenceDiagram

/-! ### 1. The one cited input — stability for functions. -/

/-- **AXIOM — stability (REF-PH §5.2; BCKL2010 §2.2/§4).**
    The diagram is 1-Lipschitz in the sup-norm: perturbing the curve moves the barcode no further
    (bottleneck) than the curve moved (sup). The single deep input; everything below is plumbing.
    Cited, not formalized — the classical CEEH function-stability theorem; BCKL2010 gives the
    statistical version with minimax rates over nonparametric curve estimators. -/
axiom diagram_stability (f g : ℝ → ℝ) :
    bottleneck (Dgm f) (Dgm g) ≤ supNorm (fun x => f x - g x)

/-! ### 2. Provable plumbing — event inclusion, then the measure pushforward. -/

/-- A sup-norm ball around `fhat` sits inside the bottleneck ball around its diagram. Pure
    consequence of `diagram_stability`. -/
theorem stability_event_subset (fhat : ℝ → ℝ) (ε : ℝ) :
    {f : ℝ → ℝ | supNorm (fun x => f x - fhat x) ≤ ε}
      ⊆ {f : ℝ → ℝ | bottleneck (Dgm f) (Dgm fhat) ≤ ε} :=
  fun _ hf => le_trans (diagram_stability _ fhat) hf

/-- **L3 — the confidence pushforward.**
    If the BARS posterior `μ` places mass `≥ 1−α` on the sup-norm credible band around the fitted
    curve `fhat`, it places mass `≥ 1−α` on the bottleneck ball around the fitted diagram. Honest
    confidence on the barcode, pushed from honest confidence on the curve — `d_b ≤ ‖·‖_∞` (axiom)
    makes the band-event a subset of the ball-event, and `measure_mono` carries the mass across.
    The "no inflation" is `BARS.lean`'s exactness: the band-ε is *also* the bottleneck-ε. -/
theorem confidence_pushforward
    {Ω : Type*} [MeasurableSpace Ω] (μ : MeasureTheory.Measure Ω)
    (f : Ω → (ℝ → ℝ))                       -- a random curve = one posterior draw
    (fhat : ℝ → ℝ) (ε α : ℝ)
    (hband : ENNReal.ofReal (1 - α)
              ≤ μ {ω | supNorm (fun x => f ω x - fhat x) ≤ ε}) :
    ENNReal.ofReal (1 - α)
      ≤ μ {ω | bottleneck (Dgm (f ω)) (Dgm fhat) ≤ ε} :=
  le_trans hband <| MeasureTheory.measure_mono <|
    fun ω hω => le_trans (diagram_stability (f ω) fhat) hω

/-! ### 3. OBLIGATION — the averageable surrogate (diagrams are not a vector space). -/

/-- Persistence landscape: the Banach-valued, stable, averageable summary of a diagram (FH2024). -/
opaque PersistenceLandscape : Type
opaque landscape : PersistenceDiagram → PersistenceLandscape
opaque landscapeSupNorm : PersistenceLandscape → PersistenceLandscape → ℝ

/-- **L5 — landscape stability (OBLIGATION; FH2024 §4).**
    The landscape map is 1-Lipschitz from (bottleneck) to (sup), so the draw-wise landscapes admit a
    Bochner mean with a CLT — the surrogate that lets "posterior over diagrams" become "mean ± band
    over landscapes." OBLIGATION: discharged once `PersistenceLandscape` + its stability are defined
    in the `TDA.Ph` layer (no Mathlib persistence-landscape theory to cite). -/
theorem landscape_stable (D D' : PersistenceDiagram) :
    landscapeSupNorm (landscape D) (landscape D') ≤ bottleneck D D' := by
  sorry  -- [OBLIGATION] FH2024 §4: 1-Lipschitz landscape (build with TDA.Ph PersistenceLandscape)

end Spc.Stability
