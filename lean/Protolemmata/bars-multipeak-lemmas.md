# BARS — multi-peak lemmas

The continuation of [`fable-BARS-lemma.md`](fable-BARS-lemma.md) into the multi-peak regime.
That note certified the *single* readout — the global argmax is a closed-form root
(`argmax_in_closed_form_set`) and peak-of-mean ≠ mean-of-peaks (`argmax_expectation_noncommute`).
These four lift the same two pillars from the argmax to the **full set of transitions**, plus the two
new pillars that only appear once the peak count `K` can exceed 1.

**Scope: the general engine, ahead of any consumer.** BARS lives under `maths/**` as a standalone,
domain-agnostic curve/function-inference engine — it is **not** part of SPC and is not consumed by it.
The classical SPC instruments (`MagnetizationPeakDetector` and the Domany-lineage T-sweep) are a
separate, already-shipped track written from the original papers; BARS does not touch them. **SPC-BARS
is net-new downstream functionality — a deliberate *parallel* track to the classical approach** — to be
built once the engine is tight and general. These lemmas certify the engine *on its own terms*, ahead of
that adoption, precisely so the adoption path is unambiguous when it comes. Where they write `χ(T)` /
`T_c` it names the eventual application, never a current coupling, and the classical detector is *not* a
co-spec to reconcile against.

**What the engine already does** (so we certify what's built, not what's imagined). Per draw,
`BarsEnsemble.Advance` (`src/maths/regression/freeknot/BarsEnsemble.cs:276`) reduces to two peak
functionals via `SplineExtrema` (`src/maths/regression/freeknot/SplineExtrema.cs`): the global
`Argmax`, and a prominence-gated `SignificantPeakCount`. The pooling is **asymmetric** — the global
peak gets a full `PeakPosterior` (location, 95 % CI, R̂, ESS); the count gets only a *mean*
(`PeakCountMean`); the 2nd/3rd transitions get *no* location posterior. So multi-peak detection is
half-built: count-per-draw exists, a posterior over *where* the non-dominant transitions sit does not.

**Disambiguation (load-bearing).** Three distinct things get called "multi-peak"; a lemma for one says
nothing about the others:

1. **Curve multimodality** — several maxima within one draw `f̃(T)`. → `SignificantPeakCount`. **These
   lemmas target this.**
2. **Posterior multimodality** — the *single* peak's location posterior is multimodal. → what
   `temperLevels` actually targets (the `BarsEnsemble.cs:66` docstring calls this "multi-peak T_c
   posterior", but it is #2, not #1). MP-2′ below pins the boundary.
3. **Change-points** — piecewise-constant segment count, the *separate* exact-DP engine
   (`ExactChangepoint`), not this path.

Genre throughout (per Lemma A/B): a *surprising exact reduction that renders apparatus vestigial*,
with the honest guard named where it bites.

---

## MP-1 — the full local-max set IS the closed-form candidate set

The lift of `argmax_in_closed_form_set` from the global maximiser to **every** local maximum: the set
of local maxima of a spline draw is finite and contained in {span boundaries} ∪ {interior
derivative-roots}. So `SignificantPeakCount` enumerates over an exact finite set — the multi-peak scan
is as vestigial as the argmax scan, zero injected error. The extra fact the C² structure buys you: `f`
is strictly monotone *between* consecutive criticals, which is exactly what makes the code's discrete
height-comparison among critical points (`SplineExtrema.cs:46`) equal the true local-max test.

```lean
import Mathlib.Analysis.Calculus.LocalExtr.Basic
import Mathlib.Algebra.Polynomial.Derivative
import Mathlib.Topology.Algebra.Order.Compact

open Set Polynomial

variable (p : Polynomial ℝ) (a b : ℝ)

/-- Fermat, multi-point form: every interior local maximiser of a span is a derivative root. -/
theorem local_max_is_critical
    (t : ℝ) (ht : t ∈ Ioo a b) (hmax : IsLocalMax (fun x => p.eval x) t) :
    (p.derivative).eval t = 0 := by
  -- IsLocalMax.deriv_eq_zero ▸ Polynomial.deriv
  sorry

/-- On a non-constant span the local-max set is finite: it injects into the (finite) root set
    of `p.derivative ≠ 0`. The finiteness the per-draw count relies on. -/
theorem finite_local_max (hp : p.derivative ≠ 0) :
    {t ∈ Icc a b | IsLocalMax (fun x => p.eval x) t}.Finite := by
  sorry

/-- The discrete-neighbour test is faithful: between consecutive criticals `f` is strictly
    monotone, so a recorded critical `cᵢ` is a local max iff `f cᵢ` exceeds both neighbours.
    REQUIRES the alternation guard (no degenerate critical — see below). -/
theorem neighbour_test_correct
    (hp : p.derivative ≠ 0) (hnd : (p.derivative).roots.Nodup) :   -- no double roots on the span
    True := by   -- statement to be sharpened to the monotone-between-criticals equivalence
  sorry
```

**Fishable:** yes — same Fermat + finiteness-of-roots core as the single-peak lemma, directly liftable
from the closed `argmax_in_closed_form_set` proof.

**Guards, where they bite (both are latent issues MP-1 forces you to settle):**

- **Degree-3 hypothesis.** The `q1/q2/q3` derivative recovery from four evaluations
  (`SplineExtrema.cs:106`) is exact *only* for cubic spans, though `basis.Degree` is read elsewhere.
  For any non-cubic basis the roots are wrong. `Degree == 3` belongs in the statement (and arguably an
  assert in the code).
- **Boundary-maximum semantics are an engine decision (and probably a configurable axis).**
  `SignificantPeakCount` counts a *boundary* maximum (`i==0`/`i==n-1`, `SplineExtrema.cs:46`) as a peak.
  Whether a *general* curve-peak-counter should is a definitional choice for the engine on its own terms —
  and because different consumers will want different answers (a thermal sweep typically rejects a
  transition sitting at the very edge of the swept bracket; a generic change-detector may not), it is
  best exposed as a *count-interior-only vs count-all* axis rather than hard-coded. The MP-1 hypothesis is
  where that choice is stated explicitly. It is **not** a reconciliation against the classical SPC
  detector — that is a separate track with its own domain convention, and slaving the engine to it would
  be a category error.

**Open formalization question:** sharpen `neighbour_test_correct` to the monotone-between-criticals
equivalence, or take "criticals partition the span into monotone arcs" as the axiom and prove the
height-comparison faithfulness from there?

---

## MP-2 — peak *count* of the mean undercounts  (LEAD)

The sharp lift of `argmax_expectation_noncommute`. Not only does peak *location* fail to survive
averaging — the peak *count* fails too, and **directionally**: averaging fills valleys, so the count of
the pooled mean curve is a systematic **under**count of the per-draw transitions. This is the formal
teeth certifying that running any detector on the pooled `r.Fit` is *wrong*, and that the per-draw
`PeakCountSum` reduce (`BarsEnsemble.cs:283`) is forced, not stylistic.

```lean
import Mathlib.Topology.Order.LocalExtr

open Set

/-- Number of local maxima of `f` on `[a,b]` (well-defined when the set is finite — the spline case). -/
noncomputable def peakCount (f : ℝ → ℝ) (a b : ℝ) : ℕ :=
  {x ∈ Icc a b | IsLocalMax f x}.ncard

/-- Count-of-mean undercounts mean-of-count: a single explicit witness (two double-bump curves whose
    peaks interleave so the average flattens to one bump). Certifies multi-peak detection MUST be
    reduce-per-draw — `peakCount r.Fit` is a biased estimate of the count posterior. -/
theorem count_of_mean_undercounts :
    ∃ (f g : ℝ → ℝ) (a b : ℝ), a < b ∧
      peakCount f a b = 2 ∧ peakCount g a b = 2 ∧
      peakCount (fun x => (f x + g x) / 2) a b = 1 := by
  -- witness: two offset quadratic double-tents; valleys of one sit under peaks of the other,
  -- the mean has a single interior maximum. One-witness existence, no theory gap.
  sorry
```

**Fishable:** yes, cleaner than the location version — a single explicit witness, no infrastructure.
The only Mathlib dependency beyond the original is `Set.ncard` finiteness for the witness, which is
discharged by the witnesses being piecewise-quadratic (finitely many criticals).

**Why lead with this:** smallest true statement that certifies the original contribution
(reduce-per-draw for the count), exact "apparatus is vestigial" genre, and it protects a real
foot-gun — peak-detecting the pooled mean curve `r.Fit` silently undercounts transitions.

**Open formalization question:** state `peakCount` over `IsLocalMax` (slicker, needs a finiteness
side-goal), or over the explicit closed-form critical Finset from MP-1 (heavier type, finiteness free)?
The latter ties MP-2 to MP-1; the former keeps MP-2 standalone.

---

## MP-2′ — curve multimodality ⟂ posterior multimodality  (boundary guard)

A guard rather than a reduction, but load-bearing because the engine's own language blurs it (#1 vs #2
above). The two "multi" are independent: a single-peaked curve family can have a *bimodal* peak-location
posterior (chains place the one peak in two basins), and a multi-peaked curve can have a unimodal count
posterior. Consequence: **R̂ on peak location detects #2, never #1** — the count posterior is the only
instrument for #1, and `temperLevels` buys mixing for #2, not detection for #1.

```lean
/-- The two notions are orthogonal: a curve-posterior with per-draw count ≡ 1 but a bimodal
    location posterior. (Existence; states that `peakCount ≡ 1` does NOT imply a unimodal T_c
    posterior, so R̂(location) cannot be read as a transition-count diagnostic.) -/
theorem count_vs_location_multimodality_independent : True := by
  -- formalised as: ∃ a draw-family with peakCount ≡ 1 whose location pushforward is bimodal.
  -- sharpen once a posterior-as-measure carrier is chosen.
  sorry
```

**Fishable:** the C# / docs consequence is immediate; the Lean statement waits on a posterior-as-measure
carrier (see the engine-design carrier deferral). Keep as a *named guard* until then — its job is to stop
`temperLevels` being sold as multi-peak *detection*.

---

## MP-3 — the prominence count is a well-defined functional, fragile only on the catastrophe set

`SignificantPeakCount` is integer-valued and **piecewise-constant** in the spline coefficients; it jumps
exactly when (a) a critical pair annihilates (the per-span quadratic discriminant crosses 0) or (b) a
peak's prominence crosses `θ·range` (`SplineExtrema.cs:41,54`). So the count posterior is a genuine
categorical pushforward, and mass split between `k` and `k+1` near a bifurcation is *correct ambiguity*,
not numerical noise — the honest representation of a marginal transition.

```lean
import Mathlib.Topology.Basic

/-- Count as a function of the coefficient vector. -/
noncomputable def sigCount (n : ℕ) (coef : Fin n → ℝ) (θ : ℝ) : ℕ := sorry

/-- FISHABLE core: away from the catastrophe set the count is locally constant
    (so the pushforward onto ℕ is a well-defined categorical law). -/
theorem sigCount_locally_constant
    (n : ℕ) (θ : ℝ) (coef : Fin n → ℝ) (hreg : Regular n coef θ) :
    ∀ᶠ c in nhds coef, sigCount n c θ = sigCount n coef θ := by
  sorry

/-- NOT FISHABLE today (state it, don't prove it): the catastrophe set
    `{coef | ¬ Regular n coef θ}` is Lebesgue-null, so under an absolutely-continuous
    posterior the count is a.s. well-defined. Needs Sard / transversality + abs-cont prior —
    the same "real but deferred" status as the approximation-economy (spline ≈ dense-grid) pillar. -/
theorem catastrophe_set_null (n : ℕ) (θ : ℝ) :
    True /- volume {coef : Fin n → ℝ | ¬ Regular n coef θ} = 0 -/ := by
  sorry
```

**Fishable:** the well-definedness/local-constancy core, yes. The measure-zero of the catastrophe set,
no — scope it like the unfishable pillars in the parent note. The `Regular` predicate (no zero-discriminant
span, no on-threshold prominence) is the explicit guard.

**Open formalization question:** is `catastrophe_set_null` worth axiomatizing now (so MP-3's a.s.
well-definedness can be *used* by downstream lemmas), or left as an out-of-scope marker until Mathlib's
transversality story is usable?

---

## MP-4 — the peak SET posterior is a point process; its intensity is the poolable summary

Design-level (not yet fishable), but it's the one that lands multi-peak BARS at the **C ∩ B ∩ P** centre
and names the engine gap. Pooling multiple peaks per draw is not a fixed-`K` vector — it is an **unordered
random finite set** on `[0,1]`. The matching-free first moment is the **intensity** `λ(B) = E[#(peaks ∩
B)]`, poolable per-draw with no label-switching; any fixed-`K` vector readout `(peak₁,…,peak_K)` is only
well-defined when `K` is a.s. constant, and under varying `K` requires an arbitrary labelling.

```lean
/-- Per-draw peak set: an unordered finite subset of [0,1]. -/
def PeakSet (f : ℝ → ℝ) : Finset ℝ := sorry   -- the MP-1 closed-form local-max set, prominence-gated

/-- The intensity (first moment) of the peak point process is additive and matching-free:
    `λ(B ⊎ C) = λ(B) + λ(C)`. The correct poolable readout. -/
theorem peak_intensity_additive : True := by
  -- λ(B) := E_draw[ (PeakSet f̃ ∩ B).card ];  additivity is linearity of expectation over disjoint B.
  sorry

/-- Ill-posedness of the fixed-K vector readout: when K varies across draws there is no
    canonical labelling, so component-wise pooling is not invariant to the labelling. -/
theorem fixedK_readout_ill_posed : True := by
  sorry
```

**Connection (the compendium bridge):** this is exactly **MNO2019**'s marked-Poisson-process model of a
persistence diagram, instantiated on the T-axis — peaks ≅ a marked point process, birth = location,
mark = height/prominence. The conjugate-intensity machinery transfers wholesale, which is what carries the
edge work (C ∩ B) into the centre (C ∩ B ∩ P).

**Engine gap it names:** to give the 2nd/3rd transition a credible interval, accumulate a per-draw peak
**intensity histogram** over the grid, not just `PeakCountSum` — the same reduce-per-draw shape, one
field instead of one scalar. That is the missing readout, and MP-4 is its specification.

**Not fishable** until a point-process / random-measure carrier exists in-tree; this is the forward
marker, the genus the other three are species of.

---

## Fishability ladder & sequencing

| Lemma | Certifies | Fishable now? | Genre |
|---|---|---|---|
| **MP-2** | reduce-per-draw forced for the **count** (mean undercounts) | ✅ one witness | exact non-commutation |
| **MP-1** | full local-max set = exact finite closed-form set | ✅ Fermat + finiteness | apparatus vestigial |
| **MP-3** | count is a well-defined categorical pushforward | ◑ core yes, null-set no | stability / measurability |
| **MP-4** | peak intensity is the matching-free poolable readout | ✗ needs PP carrier | point-process specification |

**Recommendation:** land **MP-2** first (smallest true statement, protects a live foot-gun, pure
witness), **MP-1** right behind it (the enumeration certificate, and where the endpoint + degree specs get
pinned). MP-3 stages once its catastrophe-null is either axiomatized or out-scoped; MP-4 is the design
marker that should drive the peak-intensity readout the engine is currently missing. The two MP-1 guards
— boundary-max semantics (an engine-level definitional choice, likely a configurable axis) and the
cubic-span assumption — are engine-level findings worth acting on independent of the Lean work, and
neither references SPC: they live entirely inside the general engine.
