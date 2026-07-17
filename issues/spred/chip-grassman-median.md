Adding the probe to the pilot runner — one seed at I=1000 plus the warm-start baseline in the same run, so per-block improvement is measured directly against PCA:

The probe is running in the background (~2 min predicted: warm-start baseline + I=1000 at seed 211). It will also validate the extrapolation model at its far end (predicted 71.3 s for the annealed run). I'll read out the verdict — budget-bound vs proposal-scale-bound — when it lands.

Committed as `2c84cac`. The probe answered the question cleanly:

## Verdict: proposal-scale-bound, not budget-bound

At **I=1000**, seven of eight blocks end **bit-identical** to the PCA warm start — zero accepted proposals in 1000 iterations each. Block 5 accepted exactly _one_ step, and forensically it's the first iteration: its Grassmann move is 0.1000, which is precisely the annealer's `step = temp·0.1` at initial temperature. One accepted move in 8,000 block-iterations.

The mechanism, pinned to source: [SubspaceAnnealer.cs:81-82](src/maths/geometry/dim-reduction/SubspaceAnnealer.cs:81) proposes an _isotropic random_ horizontal tangent of schedule-fixed length under `0.99^iter` geometric cooling. In Gr(617, 30) the intrinsic dimension is 30×587 = **17,610** — from a PCA warm start, the fraction of improving random directions is vanishing, and by iteration ~500 cooling has extinguished uphill acceptance anyway. More budget provably doesn't help; this cylinder-scale proposal design just doesn't transfer to high-codimension Grassmannians.

Also recorded a faithfulness correction: same-run differencing at I=1000 gives a marginal of **0.092 s/iter** — the pilot's 0.025 was small-I noise against a ±4 s fixed-share variance. Corrected screen cost: ~34 min for full S0 at I=1000. Still a non-constraint.

## Where this leaves Phase 2

**S0 is blocked on annealer mobility** — an engine question, now with redesign candidates recorded in the brief:

1. **Acceptance-targeted adaptive step** — adapt step length toward ~25% acceptance instead of a fixed `temp·0.1` schedule; minimal change.
2. **Structured two-plane proposals** — Givens-style rotations mixing one retained with one discarded direction by a small angle; the classic subspace-search move, with a far higher improving fraction than isotropic tangents at this dimension.
3. Gradient-informed Riemannian search — the bigger lift (PH subgradients via matched pairs), probably not first.

My recommendation is 1+2 together — they're small, they stay faithful to the SA design, and they attack exactly the vanishing-improving-fraction failure the probe exposed. But this is engine redesign on `SubspaceAnnealer`, so it's yours to direct — say the word and I'll build it.

## Disposition (2026-07-16, parent thread)

Triage: **execute 1+2 as one bundled engine work item; track 3.** Independent check of the
recommendation concurs, with one sharpening recorded in the brief: adaptive step alone degenerates
(it buys acceptance by shrinking toward diffusion-scale moves), so the two-plane Givens proposal
carries the fix and the adaptive scale is its controller. Design settled and deposited in
[annealer-mobility-brief.md](annealer-mobility-brief.md); dev-sequence P3 re-pointed at it; chip
issued against the brief with report-back here. The two-plane move is itself a Grassmann geodesic
(rank-1 horizontal tangent), so the design.md deviation story is unchanged. Tracked, not executed:
gradient-informed Riemannian search; restart ensembles / plateau early-stop as later add-ons. The
S0 re-probe stays with the pilot lane (`SpredIsoletPilotTests`) to avoid colliding with `2c84cac`'s
owner.
