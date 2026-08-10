# SPC × BARS readout — application obligations

Status: protolemma ledger, 2026-08-08.

Target module: `Spc.BarsReadout` (initially
`Enthymemata/Spc/BarsReadout.lean`).

This note owns the gap between the per-function facts drafted in
[`BARS.v2.lean.md`](BARS.v2.lean.md) and the readout implemented by
[`SplineExtrema`](../../src/maths/regression/spline/bars/SplineExtrema.cs).
It does **not** formalize the BARS sampler. BARS supplies posterior spline
draws and method-intrinsic candidate/crossing arithmetic; SPC supplies the
meaning of a thermal peak, the boundary policy, the prominence threshold, the
span level, and the treatment of clipped spans.

Provenance: O1 recovers the useful part of
[`fable-BARS-lemma.md`](fable-BARS-lemma.md); O2 sharpens the candidate and
neighbour-test obligations in
[`bars-multipeak-lemmas.md`](bars-multipeak-lemmas.md); O3 reconciles
[`bars-span-lemmas.md`](bars-span-lemmas.md) with the current implementation.
The larger count-instability and merge-tree direction remains recorded in
[`Fable5-BARS-revisions.md`](Fable5-BARS-revisions.md).

| Obligation | Already available | Missing deliverable |
| --- | --- | --- |
| Closed-form global peak | Interval-free polynomial Fermat lemma and finite derivative-root set | Boundary-aware argmax corollary and finite candidate set |
| Full-spline enumeration | Per-span derivative-root and level-crossing facts | Finite union over endpoints, knots, and per-span roots; faithfulness of the ordered-neighbour test |
| SPC readout contract | Raw local-max count, global level occupancy, and non-commutation witnesses | Prominence-gated peaks, boundary policy, branch-relative spans, clipping, and equivalence to `SplineExtrema` |

## O1 — boundary-aware global argmax

Recover the useful single-peak corollary referenced by both BARS drafts. The
membership result does not need a non-constant hypothesis: if the polynomial
is constant, every point satisfies the derivative-root disjunct. The
non-constant hypothesis belongs on the separate finiteness result.

Candidate statement:

```lean
/-- A maximizer on a closed polynomial span is an endpoint or a critical point. -/
theorem argmax_in_closed_form_set
    (p : Polynomial ℝ) (a b t : ℝ)
    (ht : t ∈ Set.Icc a b)
    (hmax : IsMaxOn (fun x => p.eval x) (Set.Icc a b) t) :
    t = a ∨ t = b ∨ (Polynomial.derivative p).eval t = 0 := by
  sorry
```

Then package `{a, b}` with the roots of `derivative p` and prove that this
candidate set is finite when `derivative p ≠ 0`. This is the direct
“no grid search” certificate. It certifies absence of scan or optimizer slop;
it does not claim that floating-point root evaluation is exact arithmetic.

The older `argmax_expectation_noncommute` sketch is **not** ready to promote:
as written, it puts the pooled maximum at `(a + b) / 2`, exactly the mean of
the two named endpoint locations, and therefore does not state its advertised
inequality. Replace it with a statement that represents unique selected
argmax locations, or let the count/branch-valued non-commutation results carry
the application claim.

## O2 — spline globalization and neighbour-test faithfulness

The implementation enumerates domain endpoints, every knot/break point, and
every interior derivative root of each polynomial span. Formalize that
implementation-faithful candidate set first:

```text
Candidates(draw)
  = {domain endpoints}
    ∪ {knots}
    ∪ ⋃ span, {interior roots of span.derivative}
```

Required results:

1. `Candidates(draw)` is finite when the knot set is finite and no span is
   flat. A flat span is an explicit degeneracy result, not an `ncard = 0`
   accident.
2. For a continuously glued curve, every local maximum restricted to the
   thermal window lies in `Candidates(draw)`. Including knots unconditionally
   keeps this theorem valid below `C¹`; a later `C¹` corollary may show when
   knot candidates are already critical from both adjacent pieces. If a
   discontinuous carrier is admitted later, its one-sided knot semantics need
   a separate statement.
3. Between consecutive, deduplicated candidates, a non-flat polynomial span
   has no derivative root and is strictly monotone.
4. Consequently, an interior candidate is a strict local maximum exactly
   when its value strictly exceeds the values at its two neighbours. Endpoint
   candidates are routed through the explicit SPC boundary policy rather than
   folded into this equivalence.
5. The ordered mathematical enumeration agrees with `SplineExtrema`'s
   `CriticalPoints`/neighbour comparison, subject to a separately stated
   numerical contract for root tolerance and deduplication.

The old degree-3 guard is stale and must not be revived: the implementation
now retains the cubic fast path but reconstructs and solves higher-degree span
polynomials as well.

## O3 — the actual SPC thermal-feature readout

The present `peakCount` is deliberately a raw witness functional. It is not
the specification of `SplineExtrema.SignificantPeakCount`: it has no
prominence gate, uses full-topology `IsLocalMax`, excludes restricted-window
boundary maxima, and maps an infinite weak-max set to `0` through `ncard`.

The application-level specification needs these ingredients:

- A boundary policy, at minimum `interiorOnly` versus `includeBoundary`.
- A finite peak candidate set from O2, followed by a strict local-maximum
  predicate rather than unrestricted weak `IsLocalMax`.
- Topographic prominence and the significance gate
  `prominence ≥ relativeProminence · (curveMax - curveMin)`.
- An explicit degenerate-range result when `curveMax = curveMin`.
- A branch/peak col level `c` and peak height `h`. For drop fraction `ρ`, the
  implemented crossing level is

  ```text
  ℓ = h - ρ (h - c).
  ```

  At `ρ = 1/2` this is a branch-relative half-prominence level. It equals the
  global `sSup(f)/2` convention only under additional baseline and
  single-peak assumptions.
- A span result carrying left/right crossings and independent clipping flags.
  A missing in-domain crossing is a clipped observation, not a fabricated
  root and not a failed theorem hypothesis hidden from the caller.
- Algorithm/specification equivalence for `SignificantPeaks`,
  `SignificantPeakCount`, and `SignificantPeakSpans`.

The merge/col and prominence definitions should come from the planned
carrier-generic `Pi0`/merge-tree foundation. `Spc.BarsReadout` should consume
that structure and add the thermal-domain policies; it should not define a
second, BARS-specific notion of prominence.

The current `levelWidth` remains useful under the honest name **global
occupancy above a level**. The current `halfMax`/`fwhm` pair remains a compact
non-commutation witness, but it is not the engine's per-branch span contract.

## Discharge sequence

1. Land the revised per-polynomial scaffolding and O1.
2. Introduce the minimal piecewise-polynomial spline carrier and prove O2.
3. Build the `Pi0` branch/merge-level carrier used to define prominence.
4. Specify the policy-bearing SPC peak and span readouts and prove their
   equivalence to the mathematical `SplineExtrema` algorithm.
5. Restrict the non-commutation witnesses to the admitted spline/draw carrier,
   rather than relying only on existential witnesses over arbitrary
   functions.
6. Add the posterior measure carrier and prove the peak-set/intensity and
   span-coverage pushforward statements.
