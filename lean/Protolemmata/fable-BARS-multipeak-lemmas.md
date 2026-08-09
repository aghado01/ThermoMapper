# BARS — the multi-peak slate

Companion to `fable-BARS-lemma.md`, which posed the *single*-peak pair (`argmax_in_closed_form_set`,
`argmax_expectation_noncommute`). This lifts them to the **multi-peak** readout — the per-draw
`SplineExtrema.SignificantPeakCount` that feeds `BarsResult.PeakCountMean`.

**Disambiguate "multi-peak" first** (a lemma for one says nothing about the others):

1. **Curve multimodality** — several maxima within *one* draw f̃(T). → `SignificantPeakCount`. **The slate
   below is about this.**
2. **Posterior multimodality** — the *single* peak's location posterior is multimodal (chains disagree where
   the one peak is). → what `temperLevels` targets; a different object.
3. **Change-points** — piecewise-constant segment count, the separate exact-DP engine (`ExactChangepoint`).

These are engine-level rigor for an agnostic capability; **no BARS consumer is built yet**, so nothing here
asks the classic Domany detectors to change. Lead with MP-2.

---

### MP-2 — peak count does not survive averaging (lead; the count-lift of `argmax_expectation_noncommute`)

The sharpest and highest-value statement, and it certifies a design decision already shipped: running *any*
peak detector on the **pooled** `r.Fit` curve is wrong, because averaging fills valleys and **systematically
undercounts** transitions. So the per-draw `PeakCountSum` reduce is *forced*, not stylistic.

It is a one-witness existence claim — two offset double-tent curves, each with two significant peaks, whose
average has one. Self-contained real analysis, no external dependency. **Lemma tier** (complete, no apology)
once the witnesses are formalized.

```lean
import Mathlib

/-- Number of significant local maxima of `f` on `[a,b]` at relative prominence `θ` (the functional behind
    `SignificantPeakCount`). Definition is shared with MP-1/MP-3; sorried here only because its home is MP-1. -/
noncomputable def peakCount (θ : ℝ) (f : ℝ → ℝ) : ℕ := sorry

/-- Count does not commute with averaging, and the failure has a sign: averaging can only *merge* peaks.
    Witness: two double-tents offset so each valley sits under the other's peak. -/
theorem peak_count_noncommute_undercounts :
    ∃ (f g : ℝ → ℝ) (θ : ℝ), 0 < θ ∧ θ < 1 ∧
      peakCount θ f = 2 ∧ peakCount θ g = 2 ∧ peakCount θ (fun x => (f x + g x) / 2) < 2 := by
  sorry
```

*Note:* this is the same genre as the PKWang "apparatus is vestigial" result — a small true statement that
certifies an engineering choice. It leans on `peakCount` from MP-1; until that def is real, this is stated but
unproved (an enthymeme *over a sorried def* — keep it stated, not promoted, until MP-1 lands).

---

### MP-1 — the local-max set is the closed-form candidate set (the count-lift of `argmax_in_closed_form_set`)

Every local maximum of a spline draw lies in {span boundaries} ∪ {interior derivative-roots}, a finite set, so
`SignificantPeakCount` enumerates an exact candidate set with zero scan error — exactly like the global argmax.
The extra fact the C² structure buys: `f` is monotone between consecutive criticals, so a discrete
neighbor-comparison on the recorded heights *is* the true local-max test.

```lean
import Mathlib

variable (config : KnotConfig) (coef : ℝ → ℝ) -- placeholder for the spline draw

/-- The finite candidate set: span boundaries together with interior roots of the per-span derivative. -/
def criticalCandidates : Finset ℝ := sorry

/-- Every local max is a candidate, and the set is finite.  Citation boundary: finiteness of polynomial roots
    (mathlib `Polynomial.setOf_isRoot` finiteness) — the per-span derivative is a polynomial. -/
theorem local_maxima_subset_candidates :
    ∀ x, IsLocalMax (fun t => eval config coef t) x → x ∈ criticalCandidates := by
  sorry
```

*Guard where it bites:* the monotone-between-criticals property holds only for **non-degenerate** criticals — a
double root of the per-span derivative (a horizontal inflection) breaks the max/min alternation. The hypothesis
must exclude it (companion to the existing `hp : p.derivative ≠ 0`). **Lemma tier** modulo that hypothesis +
the root-finiteness cite.

*Two code facts this lemma pins, both now settled outside Lean:*
- **Degree is no longer a hypothesis.** The old "exact only for cubic spans" caveat is *resolved in code*:
  `SplineExtrema` now reconstructs the degree-d derivative and root-finds generally (commit `afe9689b`), so the
  candidate set is exact at any degree. MP-1 is therefore degree-general; drop the `Degree == 3` hypothesis the
  earlier analysis wanted.
- **Endpoint inclusion is a consumer policy, not part of this lemma.** `SignificantPeakCount` currently treats a
  boundary maximum as a candidate; the classic `MagnetizationPeakDetector` says "endpoints are never peaks."
  MP-1 proves the candidate *set* is exact and finite — pure geometry. Whether a boundary rise *counts as a
  transition* is a downstream calibration the (unbuilt) SPC-consumes-BARS layer decides; the engine just exposes
  the candidate. Do **not** fold an endpoint convention into the lemma.

---

### MP-3 — the prominence count is a well-defined integer functional, jumping only on a fold/threshold set

`N(coef)` is piecewise-constant and integer-valued; it jumps exactly when a critical pair annihilates (a per-span
discriminant crosses 0) or a prominence crosses `θ·range`. So the count-posterior is a genuine categorical
pushforward, and mass split between `k` and `k+1` near a bifurcation is *correct ambiguity*, not noise.

```lean
/-- Well-definedness + the explicit jump set.  Fishable half. -/
theorem peak_count_piecewise_constant :
    -- N is locally constant off  D = {discriminant = 0} ∪ {prominence = θ·range}
    sorry

/-- D is measure-zero.  UNFISHABLE today: needs a Sard/transversality argument + an absolutely-continuous prior
    on coef.  State it, scope it like the spline≈dense-grid approximation pillar — do not grind it. -/
theorem jump_set_measure_zero : sorry
```

**Enthymeme** — the well-definedness is provable, but the measure-zero half is a forward boundary (Sard +
AC-prior). Real but not fishable today; keep it apologizing.

---

### MP-4 — the peak set is a point process on [0,1] (design-level; names the engine gap)

Pooling multiple peaks per draw is **not** a fixed-K vector — it is an unordered random *finite set*. The
matching-free summary is the **intensity** `λ(T) = E[peak density]`, poolable per-draw with no label-switching;
any fixed-K vector summary is ill-posed when K varies across draws.

Two payoffs, neither fishable yet:
- **It names a real engine-capability gap.** `PeakPosterior` gives the *global* peak a full posterior, but the
  2nd/3rd transitions get *no* location posterior and the count gets only `PeakCountMean`. To give a non-dominant
  transition a credible interval you accumulate a peak-**intensity histogram**, not just `PeakCountSum`. That is
  an agnostic engine readout we *could* build — a candidate future increment, not a consumer.
- **A marked-point-process treatment** (intensity conjugacy) would land this at the center of the slate. The
  earlier analysis cites "MNO2019" for peaks-as-a-marked-Poisson-process — **verify that reference is actually
  in-corpus before relying on it** (it is not in the BARS compendium index; may be misattributed).

```lean
/-- The per-draw peak set as a finite subset of [0,1]; its expected counting measure is the intensity λ. -/
def peakSet (config : KnotConfig) (coef : ℝ → ℝ) : Finset ℝ := sorry
-- Design note, not yet a theorem: λ(T) = E_draws[ #(peakSet ∩ dT) ] is the well-posed pooled summary.
```

**Enthymeme / design note** — not posed as a theorem yet; recorded because it both names the engine gap and
points at the conjugate-intensity machinery that would close it.

---

### BARS-S — the sufficiency premise (and its self-check)

The premise every lemma above silently assumes: that the per-draw curve f̃(T) was *faithfully recovered from
phase-1*. MP-1..4 all run on the candidate set / count of that recovered curve — so if phase-1 under-resolved the
structure, they are exact statements about the *wrong* curve. Sufficiency is logically prior to the whole slate.

It is also what makes the "clip is moot under BARS" reading hold (the arch thread's `75240f4c` clip-semantics
reframe): BARS's phase-1 is uniform over all of [0,1], so there is no interior sampling gap — sufficiency ⇒ the
peaks appear over the full domain ⇒ no *coverage* clip. The clip flag does not vanish; it migrates to a
*bracket-adequacy / edge-transition* signal, a domain premise the (unbuilt) consumer owns.

The elegant part — and why it belongs beside the span machinery — is that **the spans audit their own premise**. A
peak whose FWHM span is comparable to the phase-1 grid spacing was under-resolved, so the narrowest returned span
sets a Nyquist-ish floor that certifies, after the fact, whether the sparse grid was dense enough to trust the
placement. The FWHM extension both *uses* sufficiency (to place the deep grid) and *measures* it.

```lean
variable (f : ℝ → ℝ)  -- the true response on [0,1]

/-- Sufficiency: there is a phase-1 spacing Δ (set by the finest structural scale of f) at which the BARS fit on
    the uniform Δ-grid recovers the true landmark set.  FORWARD BOUNDARY, not fishable today: a Nyquist /
    approximation-theory bound on spline-fit recovery vs grid density — scope it like the spline≈dense-grid pillar. -/
theorem phase1_density_suffices : sorry

/-- The measurable shadow: the narrowest significant-peak span width lower-bounds the resolution the phase-1 grid
    achieved, so the returned spans diagnose first-pass sufficiency post hoc. -/
theorem spans_audit_sufficiency : sorry
```

**Enthymeme.** The sufficiency bound itself is the forward boundary (Nyquist-for-splines — real but not fishable
today); its concrete, auditable half is the span widths already computed in `SignificantPeakSpans`. This is the
load-bearing premise the FWHM-span work rests on, and the cleanest tie between the multi-peak readout and the
adaptive-schedule role BARS plays downstream.

---

**Promotion order:** MP-2 leads (smallest true statement, certifies the shipped per-draw reduce) but rides on
MP-1's `peakCount` def, so MP-1 must land first to give MP-2 something non-sorried to stand on. MP-3 and MP-4
each carry an unfishable half — stage them as enthymemata; let MP-4 be the one that motivates the peak-intensity
readout the engine is currently missing. **BARS-S sits *under* the slate** — the premise MP-1..4 quantify over —
and lands as an enthymema whose measurable half (spans-as-resolution) is already in code.
