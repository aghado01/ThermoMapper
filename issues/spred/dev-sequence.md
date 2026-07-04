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
  the PH objective + driver belong in `TDA.Ph` / above.
- Tests green: same-seed bit-identical determinism + orthonormal projection rows.

Files:
- Engine: `src/maths/geometry/dim-reduction/SubspaceAnnealer.cs`
- Tests: `tests/geometry/dim-reduction/SubspaceAnnealerTests.cs`
- Objective ingredients (already present): `src/tda/ph/DiagramMetrics.cs`,
  `src/tda/ph/PersistenceBarcode.cs`, `src/tda/ph/RipsFiltration.cs`, `src/tda/ph/LazyRipsFiltration.cs`

## Critical path

### 1 — PH objective (`PersistenceObjective`)  →  unblocks step 2
The `SubspaceObjectiveFunction` the driver injects (as `PersistenceObjective.Evaluate`). **API verified 2026-07-03.** The pipeline is real and
SPRED-anticipated: `RawDistanceWeights` is documented "CSR edge weights as Rips filtration values
(SPRED SA path)", and `Barcode` names "SPRED cost" as a consumer.

**Graph construction is injected, not hard-wired (AG steer 2026-07-03).** The projected graph is
built by the full `GraphCompiler` from a caller-supplied `GraphCompilerConfig` recipe — kNN/ε-ball,
OrRule/MutualKnn, MST repair, PathNeighbor refinement, and pathology interrupts are all the caller's
to turn. Injection seam: `GraphMetric.FromFeatures(projectedPoints, metric?)` — the compiler is
metric-driven (`Distance: Func<int,int,double>`), so handing it the projected coordinates runs the
entire recipe in the projected space. Optional `IDistanceMetric` allows a non-Euclidean projected
space (Poincaré / spherical / Fisher-Rao) too.

**Pipeline (per proposal P):**
`features_P = P·X` → `GraphMetric.FromFeatures(features_P)` → `GraphCompiler.Build(config, n, metric).Graph`
→ `RipsFiltration.RipsFromGraph(g, filtrationWeights, maxDim)` → `PersistentHomology.Compute(g, maxDim)`
→ `Σ_dim w_dim · DiagramMetrics.Wasserstein(B_P, B_ref, dim, p, essential)`.
`B_ref` built once the same way on ambient `GraphMetric.FromFeatures(X)` under the same recipe.

**Skeleton model — RESOLVED: recompile per proposal.** The compiler's stages are distance-driven (the
kNN edge *set* moves with P), so full-recipe fidelity ⇒ rebuild the graph each proposal. Correctness-
first, and it honors the "leverage the compiler" steer. SA tolerates the combinatorial discontinuity
(global stochastic optimizer, no continuity needed). The earlier fixed-skeleton + reflow-weights idea
is demoted to a step-3a perf mode (valid only for RawDistance on a frozen ambient topology).

**Filtration-value axis (also injected):** `RipsFromGraph` takes a `FiltrationWeights` — `RawDistance`
(needs a DistanceProjection graph) or `EffectiveResistance` (Laplacian-derived, any graph). Expose it.

**Placement (supersedes design.md row 2).** Not a new `TDA.Ph` type — TDA.Ph holds every barcode/
distance primitive but does not reference `Graphs`. The objective (`PersistenceObjective`) is a **type
assembled in the consumer/driver layer** (step 2), which references `Graphs` + `TDA.Ph`. Steps 1 and 2 merge.

**Robustness — SA must always get a finite, comparable objective value:**
- Catch `GraphPathologyException` from the compiler → return a large finite `PathologyPenalty` (SA
  rejects the proposal) rather than crashing.
- Essential policy `FinitePenalty`, never `InfiniteOnMismatch` (H0 always carries an essential bar, so
  counts routinely mismatch and the `+∞` would break the SA acceptance ratio).

**Proposed surface (declarative DTO; fluent shell at the REPL, per strict-core-fluent-shell) —
paper-grounded 2106.02096:**
`PersistenceObjectiveConfig { required GraphCompilerConfig Graph; GraphCompilerConfig? ReferenceGraph;
FiltrationWeights Filtration; IDistanceMetric? ProjectedMetric; int MaxDimension=2;
(int Dim,double Weight)[] Dimensions=[(0,0.5),(1,0.5)]; double WassersteinOrder=2;
EssentialPolicy Essential; double PathologyPenalty=1e6; double VarianceRegularizer=0.0 }`.
`new PersistenceObjective(data, config)` — builds `B_ref` once; exposes `Barcode ReferenceBarcode` and
`double Evaluate(double[][] projection)` (a `SubspaceObjectiveFunction`, passed to the engine as a method
group). The type carries the state a closure would have hidden.
- `Dimensions` = the paper's `(λ, 1−λ)` multi-order weighting (§6); `ReferenceGraph` null → reuse
  `Graph` (paper-faithful — same construction both sides).
- `VarianceRegularizer` = §6 PCA-spirit term `w·tr(P Σ_X Pᵀ)`; **negative** rewards variance (PCA
  maximizes that trace, so a positive coefficient penalizes it). Default 0 (off).
- `WassersteinOrder` = p; the ground-metric **q is fixed at L∞** by `DiagramMetrics` (the paper allows
  general q) — extend `DiagramMetrics` if a q knob is ever wanted.

**Settled sub-choices:** `maxDim = 2` (H0+H1; loops need triangle fillers); Wasserstein `p = 2`.
`Bottleneck` stays out — still `NotImplementedException` (P1). Full Rips/Čech is recoverable as a dense
recipe; Čech proper deferred (nerves/ unwired here).

### 2 — SPRED driver (new consumer)  →  depends on step 1
Composes `SubspaceAnnealer` (geometry) with the step-1 objective (ph). Named `Spred` (the reserved name).
- **Placement decision (open):** a dedicated consumer project (e.g. `TDA.DimReduction`) referencing
  both `Maths.Geometry` and `TDA.Ph`, versus adding a `TDA.Ph → Maths.Geometry` edge and hosting it
  in `TDA.Ph`. Prefer the dedicated consumer — it keeps `TDA.Ph` as PH *primitives* rather than DR
  *applications*, and matches the split's intent. No cycle risk either way (Geometry does not
  reference TDA.Ph).
- Surface: `Spred.Compute(data, targetDim, phOptions, maxIters, seed)` → delegates objective construction
  to step 1, optimization to the engine.

## After the driver (parallelizable)

### 3a — Performance
Each objective evaluation is a full k-dim Rips-PH over n points; the anneal runs `maxIters` of them, so PH
dominates wall-clock. Options to weigh (measure first): `LazyRipsFiltration`, landmark/subsampling
the cloud, restricting to H0/H1, and tuning the cooling schedule / iteration budget against a
proposal that moves further per accepted step. Do not silently cap — record any subsampling.

### 3b — Validation (R-oracle) + topological-equivalence measures
Per `validation-independence`: ground truth must originate outside the estimator's own model.
Engine-level manifold-opt facts are green; owed:
- End-to-end SPRED parity against an external oracle (Rdimtools / maotai are installed in the in-repo
  `r/` renv — see `project_kisungyou_dr_track`). **Parity may require a Stiefel-QR PERTURB mode** on
  the engine (add i.i.d. Gaussian to every entry + QR, per paper §3.1) to match the reference
  apples-to-apples — our Grassmann-geodesic walk is a deliberate improvement (design.md), so the
  oracle comparison needs the paper-faithful walk as an option.
- **§4 topological-equivalence measures** `μ_quasi-iso` / `μ_equiv` — post-hoc quality metrics via the
  filtration homomorphism (not the SA objective). `μ_quasi-iso` is barcode-computable (Prop 4.3: a
  height-matching sweep over `B_X`, `B_Y` with the η shift); `μ_equiv` needs π₁ of a quotient complex
  (SageMath-grade, hard) — defer. These make a natural reporting layer over an optimized projection.
- Scale-calibration and the paper's eigengap block weighting (distributed track).

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
