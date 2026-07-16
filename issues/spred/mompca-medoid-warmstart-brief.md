# Brief: medoid warm start for the MoM-PCA Grassmann median

**Status:** closed — executed in the parent thread 2026-07-16 (see Report)
**Scope:** `src/maths/geometry/dim-reduction/MoMPCA.cs` + `tests/geometry/dim-reduction/` (K.You DR track)
**Pattern source:** `DistributedSpred.AggregateProjections` fix, 2026-07-16 (see Report section for the landed shape)

## Context

`DistributedSpred.AggregateProjections` had a robustness defect: the Grassmann geometric
median was warm-started at `frames[0]`. When the leading block is corrupted and sits at
principal angle exactly π/2 from the clean majority (YᵀZ singular — the Grassmann cut
locus), `GrassmannManifold.LogMap` degenerates (the π/2 component's U₁ column is zeroed by
the thin SVD), the Weiszfeld/IRLS subgradient vanishes, and the iteration never leaves the
corrupted initialization. Fixed by warm-starting at the **medoid** frame (min total
Grassmann distance to all frames, O(k²) distance calls, ties → lowest index). This keeps
the warm start inside the majority cluster (within LogMap's injectivity radius) and makes
the aggregate order-invariant. Regression tests:
`AggregateProjections_CorruptedFirstBlock_StillFindsCleanMajority`,
`ComputeWithDiagnostics_CorruptedFirstBlock_AggregatesCleanMajority`, and the in-domain
contrast `AggregateProjections_TiltedCorruptedFirstBlock_ConvergesToCleanMajority` in
`tests/tda/dim-reduction/DistributedSpredTests.cs`.

## Defect (sibling instance)

`MoMPCA.cs` line ~62:

```csharp
double[] framePrelim = (double[])blockFrames[0].Clone();
GeometricMedian.Compute(grass, blockFrames, weights, framePrelim);
```

Identical `blockFrames[0]` warm start. MoM-PCA is *specifically* a
median-of-means contamination-robustness construct, so a contaminated leading block is
in-scope by design, and `framePrelim` additionally feeds the scale calibration (`CalibrateScale`)
and the joint product-median warm start — a stalled preliminary median poisons both.

## Fix recipe

1. Compute the medoid index over `blockFrames` under `grass.Distance` (accumulate pairwise
   totals symmetrically, strict `<` so ties resolve to the lowest index — deterministic).
   Mirror the private `MedoidIndex(GrassmannManifold, double[][])` helper in
   `DistributedSpred.cs`, or hoist a shared helper next to `GeometricMedian` in
   `Maths.Geometry.Estimators.Intrinsic` and re-seat the DistributedSpred call site on it
   (preferred if naming is settled with Azriel first — don't coin API surface unilaterally).
2. Warm-start `framePrelim` at the medoid clone.
3. Consider whether `productMedian` (the joint product-manifold median warm start further
   down) needs the same treatment — it is built from `meanPrelim`/`framePrelim`, so it
   inherits the fix if step 2 lands, but verify.
4. Do NOT touch `GeometricMedian` itself — its contract is explicitly caller-supplied warm
   start (`ManifoldMedian.cs`), and the medoid selection is caller policy.

## Test recipe

Port the cut-locus fixture from `DistributedSpredTests.cs`: clean-majority blocks spanning
the xy plane (plus small in-plane rotations/tilts), corrupted leading block spanning xz —
every clean frame contains the y axis, xz is normal to it, so the old warm start sits at
principal angle exactly π/2 from the whole majority. Assert the MoM-PCA frame lands near xy
(Grassmann distance < 0.1) and strictly closer to xy than xz. Add the tilted (0.3 rad)
in-domain contrast variant. Run the geometry test project to verify.

## Report

**2026-07-16 — executed in the parent thread** (chip superseded; if a chip session picked
this up, stop — the work below is already in the working tree).

- **Hoisted the helper** (fix-recipe step 1, hoist variant): `GeometricMedian.MedoidIndex<TManifold>`
  now lives in `src/maths/geometry/estimators/intrinsic/ManifoldMedian.cs`, generic over
  `IRiemannianManifold` like `Compute`, ArrayPool scratch per the file's zero-alloc discipline,
  strict `<` ties → lowest index. Justified as shared surface: MoMPCA is the second consumer.
- **Re-seated `DistributedSpred.AggregateProjections`** on the shared helper; its private
  `MedoidIndex` deleted.
- **Fixed `MoMPCA.ComputeMoM`**: `framePrelim` warm-starts at the medoid frame. `productMedian`
  inherits the fix (seeded from `framePrelim`); `CalibrateScale` uses only `Distance` (total
  function, no LogMap) so its exposure was a poisoned α̂ from a stalled center — also resolved.
- **Tests** in `tests/geometry/dim-reduction/MoMPCATests.cs`:
  `ComputeMoM_CorruptedFirstBlock_RecoversCleanSubspace` (cut-locus guard, 5 blocks: xz circle
  leading, 3× xy, yz trailing) and `ComputeMoM_TiltedCorruptedFirstBlock_RecoversCleanSubspace`
  (0.3 rad in-domain contrast).
- **Mutation-verified**: reverting to `blockFrames[0]` fails the cut-locus guard at exactly π/2
  (hard stall) and the tilted contrast at ≈0.13 (no stall, but the α̂ calibration + joint median
  seeded off the imperfect preliminary miss the 0.1 tolerance) — so in MoMPCA even in-domain
  contamination in block 0 degraded the result under the old warm start; both effects are
  covered.
- **Suites**: Maths.Geometry.Tests 17/17, TDA.DimReduction.Tests 57/57 (Release). Uncommitted.
