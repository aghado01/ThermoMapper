# SPRED — Dev Sequence

**Updated:** 2026-07-03

Ordered build plan for the SPRED track. For the architecture and the rationale behind each layer,
see [design.md](design.md).

## Done

### 0 — Engine (`SubspaceAnnealer`) ✔ 2026-07-03
- Renamed `Spred` → `SubspaceAnnealer`, reserving the name *SPRED/Spred* for the driver.
- Rebuilt on `GrassmannManifold`: horizontal-tangent proposals retracted along Grassmann geodesics,
  replacing the hand-rolled perturb-then-Gram-Schmidt.
- Doc-comment corrected — the class states it is the engine, credits the SPRED seed, and records that
  the PH cost + driver belong in `TDA.Ph` / above.
- Tests green: same-seed bit-identical determinism + orthonormal projection rows.

Files:
- Engine: `src/maths/geometry/dim-reduction/SubspaceAnnealer.cs`
- Tests: `tests/geometry/dim-reduction/SubspaceAnnealerTests.cs`
- Cost-adapter ingredients (already present): `src/tda/ph/DiagramMetrics.cs`,
  `src/tda/ph/PersistenceBarcode.cs`, `src/tda/ph/RipsFiltration.cs`, `src/tda/ph/LazyRipsFiltration.cs`

## Critical path

### 1 — PH cost adapter (in `TDA.Ph`)  →  unblocks step 2
The function the driver injects as the `SubspaceCostFunction`.
- **Shape:** precompute the reference barcode `PH(X_ambient)` **once**; per proposal, build a Rips
  filtration on the projected cloud `project(X)` (k-dim), extract its barcode, and return a weighted
  sum over homological dimensions `Σ_p w_p · DiagramMetrics.Wasserstein(proj_p, ref_p, p)`.
- **Open choices:** which `H_p` to match (H0 alone is cheap and often enough; H1 adds loop
  fidelity); the Wasserstein order `p`; the essential-bar policy (`DiagramMetrics.EssentialPolicy` —
  `InfiniteOnMismatch` vs `FinitePenalty`).
- **Constraint:** `DiagramMetrics.Bottleneck` is `NotImplementedException` (marked P1) — gate on
  `Wasserstein` until it lands.
- **Verify first:** the `RipsFiltration` / `PersistenceBarcode` API surface (constructor inputs, how
  a barcode-per-dimension is returned) before wiring — the recent migration may have shifted
  signatures.

### 2 — SPRED driver (new consumer)  →  depends on step 1
Composes `SubspaceAnnealer` (geometry) with the step-1 cost (ph). Named `Spred` (the reserved name).
- **Placement decision (open):** a dedicated consumer project (e.g. `TDA.DimReduction`) referencing
  both `Maths.Geometry` and `TDA.Ph`, versus adding a `TDA.Ph → Maths.Geometry` edge and hosting it
  in `TDA.Ph`. Prefer the dedicated consumer — it keeps `TDA.Ph` as PH *primitives* rather than DR
  *applications*, and matches the split's intent. No cycle risk either way (Geometry does not
  reference TDA.Ph).
- Surface: `Spred.Compute(data, targetDim, phOptions, maxIters, seed)` → delegates cost construction
  to step 1, optimization to the engine.

## After the driver (parallelizable)

### 3a — Performance
Each cost evaluation is a full k-dim Rips-PH over n points; the anneal runs `maxIters` of them, so PH
dominates wall-clock. Options to weigh (measure first): `LazyRipsFiltration`, landmark/subsampling
the cloud, restricting to H0/H1, and tuning the cooling schedule / iteration budget against a
proposal that moves further per accepted step. Do not silently cap — record any subsampling.

### 3b — Validation (R-oracle)
Per `validation-independence`: ground truth must originate outside the estimator's own model.
Engine-level manifold-opt facts are green; owed is end-to-end SPRED parity against an external oracle
(Rdimtools / maotai are installed in the in-repo `r/` renv — see `project_kisungyou_dr_track`),
scale-calibration, and the eigengap block weighting from the paper.

### 3c — Application: ISOLET
SPRED is unsupervised + linear, so it is not blocked by validation-independence but is bounded by the
unsupervised ceiling (the PCA-front-end wall in `project_isolet_pca_wall`). Track SPRED vs raw-617-d
and vs PCA-14 once the driver exists.

## Independent track (start anytime — needs only the engine)

### Distributed SPRED (§3.2)
Per-block simulated annealing yields one subspace per block; aggregate them by **geometric median on
the Grassmann manifold** (Weiszfeld/IRLS), reusing the existing
`GeometricMedian.Compute<GrassmannManifold>` / `RobustDistributedPCA` infrastructure. This is where a
Stiefel/Grassmann *estimator* (not the annealer) is the right tool. Ties into the
`project_kisungyou_dr_track` robust-DR work.
