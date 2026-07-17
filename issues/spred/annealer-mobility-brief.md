# SubspaceAnnealer mobility — proposal redesign brief

**Status:** ready to execute
**Date:** 2026-07-16
**Home:** engine layer — `src/maths/geometry/dim-reduction/SubspaceAnnealer.cs` + `tests/geometry/dim-reduction/`
**Seed evidence:** [chip-grassman-median.md](chip-grassman-median.md) iteration-budget probe (`2c84cac`)

## The finding this executes against

At I=1000 on ISOLET S0 (Gr(617, 30), intrinsic dimension 30×587 = 17,610), seven of eight blocks end
bit-identical to the PCA warm start — one accepted move in 8,000 block-iterations, and that one is the
very first step at full temperature. Compute is a non-constraint (corrected marginal 0.092 s/iter,
~34 min full S0); the block is proposal design. Mechanism: an isotropic random horizontal tangent's
directional derivative shrinks like 1/√dim while the finite-step curvature penalty stays O(step²), so
the improving fraction of fixed-length isotropic proposals vanishes at high codimension, and the
`0.99^iter` cooling extinguishes uphill acceptance by ~iteration 500. More budget provably does not
help. The cylinder-scale proposal design does not transfer.

## Design (settled)

1. **Two-plane Givens proposals — the primary move.** Pick a retained column i ∈ [k] and a unit
   direction v orthogonal to span(Y) (draw a random d-vector, project off span(Y) — the existing
   `HorizontalTangent` projection specialized to one column — and normalize); rotate
   `y_i ← y_i·cosθ + v·sinθ`. This is exactly the Grassmann geodesic for the rank-1 horizontal
   tangent `θ·v·e_iᵀ` (its SVD has a single σ = θ), so it stays inside the design.md
   "geodesic SA on Grassmann" story — `ExpMap` may be reused as-is or specialized closed-form.
   Rationale: restores O(1)-dimensional moves whose improving fraction does not vanish with
   codimension; the classic subspace-search move.
2. **Acceptance-adaptive step scale.** Adapt the θ scale toward ~25% acceptance (multiplicative or
   windowed update) with a floor and ceiling, replacing the fixed `temp·0.1` magnitude. Cooling
   continues to govern the Metropolis temperature only. Adaptive-alone is known-insufficient: it
   buys acceptance by shrinking toward diffusion-scale steps; it is the controller inside the
   two-plane move, not the fix itself.
3. **Isotropic fraction.** A small isotropic-proposal mixture (default candidate 10%) is optional
   ergodicity insurance — decide during implementation; pure two-plane is acceptable if the cylinder
   validates.
4. **Config surface.** A declarative `SubspaceAnnealerOptions` record (proposal kind/mixture, target
   acceptance, step floor/ceiling, cooling), defaults reproducing sensible current behavior;
   `Compute` takes the options. One Xoshiro stream; adaptation is state-dependent but deterministic
   given the seed — same-seed bit-identical runs remain the contract.

## Validation

- **New engine-level mobility fact (the key test):** a synthetic quadratic-on-Grassmann objective at
  high codimension (e.g. d=200, k=5, minimize a distance/trace objective to a target subspace) —
  the current isotropic proposals stall from a deliberately offset start, the new proposals descend.
  Fast, deterministic, no PH dependency.
- Cylinder fixture still beats every axis-aligned view (`SpredCylinderTests`) — no small-d regression.
- Same-seed determinism facts updated for the options surface.
- The S0 acceptance/objective re-probe belongs to the pilot lane (`SpredIsoletPilotTests` — owned by
  the concurrent pilot session; coordinate, do not collide).

## Tracked, not executed

- **Gradient-informed Riemannian search** (PH subgradients via matched pairs) — the bigger lift;
  revisit only if two-plane + adaptive stall on S0.
- **Screened multi-proposal steps (best-of-M per iteration)** — phase-2 layer on the two-plane
  design, exploiting the oracle/proposal cost asymmetry (~0.092 s exact eval vs microsecond
  candidate generation): generate M two-plane candidates, screen with a cheap stage, spend the
  exact evaluation only on the screen's best. The two-plane structure makes screening O(n) per
  candidate after one O(nd) pass per fresh direction (a Givens move changes one projected
  coordinate), and the cheap stage already exists as the oracle-validated fidelity ladder
  (`DistanceKind` sliced/Sinkhorn + landmark subsample + coarse `MinPersistence`). Principled
  variant when reversibility matters: multiple-try Metropolis (Liu–Liang–Wong); greedy best-of-M
  is fine for optimization-mode SA with best-so-far tracking. **Gate before building:** measure
  surrogate fidelity — rank correlation between the cheap score and the true Δobjective on the
  cylinder plus one S0 block; a misranking screen only adds cost.
- **Fidelity-cascade restarts (successive halving)** — the restart-ensembles / plateau-early-stop
  remainder of dev-sequence P3, upgraded to a Hyperband shape: many seeds × short budget × cheap
  fidelity, keep top-K by objective, graduate survivors to exact Hungarian on full data. Pure
  orchestration over existing config knobs; no new math.
- **SPC-side transfer (parked here until an SPC issues home exists):** cascade top-K's native
  shape — score many candidate basins cheaply, prune with a persistence guardrail, refine
  survivors — belongs to the SPC resolver/hierarchy machinery, where it is nearly already true:
  selective clustering is a purity-coverage gate, and lineage persistence across the temperature
  filtration is a survival criterion. Candidate application: basin/lineage ranking in the ISOLET
  Track B resolver work, with persistence-above-threshold as a prune exemption.

Seed credit: AG's XBridge note "Topo Top K" (Perplexity ideation, captured 2026-07). The surviving
core is the cascade economics — never spend the expensive oracle on an unscreened candidate; the
note's basin-pruning framing is the SPC-side item above. Its citations are keyword noise; the
load-bearing anchors are multiple-try Metropolis and successive halving/Hyperband.

## Report-back

Append implementation and probe results here, per chip conventions.

### 2026-07-16 — implemented; geometry + TDA.DimReduction suites green

**What landed** (`SubspaceAnnealer.cs`, new `SubspaceAnnealerOptions.cs`, thin call-site updates in
`Spred.cs`/`DistributedSpred.cs`, rewritten `SubspaceAnnealerTests.cs`):

- **Two-plane Givens primary move**, closed-form rank-1 specialization of `ExpMap`:
  `y_i ← y_i·cosθ + v·sinθ` with v Gaussian-drawn and projected off span(Y), θ = scale × N(0,1).
  Isotropic kind retained behind the mixture knob.
- **`SubspaceAnnealerOptions`** — IsotropicFraction (default 0.1, kept as ergodicity insurance; the
  cylinder validates), TargetAcceptance 0.25, InitialStep 0.1, StepFloor 1e-3, StepCeiling π/2,
  InitialTemperature 1.0, CoolingRate 0.99. `Compute(data, targetDim, objective, maxIters, seed,
  options, ct)` — result surface unchanged. Same-seed bit-identical contract kept (one Xoshiro
  stream); determinism facts cover defaults and a non-default options record.
- **Deviation the validation forced: per-move-coordinate step scales.** A single pooled scale is
  structurally strangled in the refinement regime — acceptance is capped near
  (mobile columns / k)·50%, so when the 25% target is unreachable the controller drives the pooled
  step to the floor and kills the mobile column with it. Now one multiplicative controller per
  retained column plus one for the isotropic kind, sharing the floor/ceiling and the zero-drift
  update law (grow^p · shrink^(1−p) = 1 at p = target). Deterministic; no options-surface change.

**Mobility fact** (`Compute_HighCodimension_TwoPlaneDescends_WhereIsotropicStalls`): Gr(5, 200),
chordal-distance² objective to span{e_0..e_4}, exact-rank-5 data placing the PCA warm start on
span{e_0..e_3, e_100} — aligned in four columns, offset in one, the ISOLET situation in miniature.
2000 iters, seed 7, temperature 1e-3 (commensurate with per-move increments — melt controlled for):
two-plane at shipped defaults **0.642**; fixed-length isotropic (old kind, pinned 0.1) **1.000000
bit-frozen**; adaptive-step isotropic **1.000000 bit-frozen** — adaptive-alone confirmed
insufficient, exactly as the design section predicted. Runs in ~1 s, no PH dependency.

**Findings that matter for the S0 re-probe** (pilot lane):

1. **Degenerate-pair fragility.** When the PCA eigen-tail is near-degenerate, the offset direction
   smears across ≥2 retained columns; every single-column Givens then trades an O(q_mixed) loss
   against the 1/(d−k) sliver and descent collapses to a tiny-θ cross-term crawl (observed at
   fixture scale: a 0.49/0.59 mixed pair → ~1e-6/move). ISOLET's tail eigenvalues at k=30 are
   plausibly closely spaced — if S0 still crawls, suspect this before suspecting compute; the
   targeted extension is a paired (rank-2) column move, cheaper than gradient-informed search.
2. **Temperature melt.** InitialTemperature 1.0 with 0.99-cooling accepts wholesale damage to a
   good warm start whenever per-move objective increments ≪ 1. Default kept (schedule out of scope
   here), but the S0 probe should set the initial temperature commensurate with observed proposal
   increments.
3. Even at intrinsic dim 975, *fixed-scale* isotropic proposals descend on a smooth objective from
   a fully-orthogonal start (measured 1.70 at 4000 iters from f₀ = 5) — the hard stall needs the
   good-warm-start + sliver-offset geometry, which the engine fact encodes deliberately.

**Verification note:** TestHarness.Runner can't locate the repo root right now (expects
`changelog.md` beside `Directory.Build.props`; the file is absent), so suites ran via plain
`dotnet test`: Maths.Geometry.Tests 21/21; TDA.DimReduction.Tests 66/66 with
`--filter "Category!=Benchmark"` (pilot facts untouched, per coordination).
*(Root anchor fixed in `ec0e9d4`.)*

**Parent-thread addendum (2026-07-16):** finding 1's premise is already evidenced, not
hypothetical — the ISOLET benchmark brief's diagnostics record covariance effective rank ≈ 45.4,
participation ratio ≈ 12.1, and 95% variance at 173 dims, so the k = 30 cut sits inside a flat
stretch of the spectrum where near-degenerate tail pairs are the expected case. Cheap pre-check
before any S0 budget: read out the eigengaps λ₂₅…λ₃₅ from the existing PCA — if the gaps near the
cut are small relative to their magnitudes, treat the paired rank-2 column move as a prerequisite
rather than a fallback, and calibrate `InitialTemperature` from the I=1000 probe's recorded
per-move increments in the same pass (finding 2).
