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

### 1–2 — Objective + driver (`PersistenceObjective`, `Spred`) ✔ 2026-07-03
`TDA.DimReduction` project stood up (csproj + test project + sln wiring), referencing Maths.Geometry
+ TDA.Ph + Graphs.Proximity. `PersistenceObjective` (recompile-per-proposal via GraphCompiler → Rips →
PH → Wasserstein; `Evaluate` + `ReferenceBarcode`; auto `FinitePenalty` at diam/2; `GraphPathologyException`
→ `PathologyPenalty`; optional variance regularizer) and the `Spred` driver land. Smoke tests green:
the objective scores a loop-preserving projection strictly below a loop-collapsing one on a 3D circle
(β₁=1), and `Spred.Compute` runs end-to-end returning an orthonormal loop-preserving projection.

Files: `src/tda/dim-reduction/{PersistenceObjectiveConfig,PersistenceObjective,Spred}.cs`,
`tests/tda/dim-reduction/SpredSmokeTests.cs`.

## Critical path (landed — detailed specs below)

### 1 — PH objective (`PersistenceObjective`) ✔ landed 2026-07-03
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

### 2 — SPRED driver (`Spred`, new consumer) ✔ landed 2026-07-03
Composes `SubspaceAnnealer` (geometry) with the step-1 objective (ph). Named `Spred` (the reserved name).
- **Placement decision (open):** a dedicated consumer project (e.g. `TDA.DimReduction`) referencing
  both `Maths.Geometry` and `TDA.Ph`, versus adding a `TDA.Ph → Maths.Geometry` edge and hosting it
  in `TDA.Ph`. Prefer the dedicated consumer — it keeps `TDA.Ph` as PH *primitives* rather than DR
  *applications*, and matches the split's intent. No cycle risk either way (Geometry does not
  reference TDA.Ph).
- Surface: `Spred.Compute(data, targetDim, phOptions, maxIters, seed)` → delegates objective construction
  to step 1, optimization to the engine.

## Performance & Scale

The recompile-per-proposal loop is the wall: one graph-build + PH + Wasserstein, ×`maxIters`.
**Measure before optimizing** (repo discipline) — two very different terms can dominate.

Per-`Evaluate` hotspot inventory:

| step | complexity | note |
|---|---|---|
| projected kNN | ~O(n²·k) | brute-force `pairDistance` (compiler can't index an arbitrary metric); k = target dim, tiny |
| **H0 Wasserstein** | **O(n³) time, O(n²) mem** | dense (n+m)² cost matrix + Hungarian; #H0 bars ≈ n. H1 is cheap. The sleeper cost — and why the paper runs order-0 as its own mode. |
| PathNeighbor refinement | ~O(n·(E + n log n)) | bounded SSSP (Euclidean pass-through is retired) |
| PH reduction | superlinear in #simplices | ≈ O(n·k + triangles) |
| reference barcode | once | amortized |

Prime suspects are projected kNN and H0-Wasserstein; which dominates depends on n and whether H0 is
matched — hence P0.

### P0 — Profile first
A small harness reporting the per-`Evaluate` breakdown (kNN vs PH vs Wasserstein vs refinement) and
total wall-clock at n ≈ few-hundred → few-thousand (cylinder / iris / synthetic). Optimization order
follows the profile, not guesswork.

### P1 — Per-evaluation speed
- **Fixed-skeleton mode** (the demoted skeleton option): freeze the ambient topology, reflow only edge
  weights per proposal → skip per-eval kNN entirely (valid for `RawDistance`). Biggest structural win;
  trades neighborhood-fidelity for speed.
- **KD-tree kNN** in the tiny projected space (k = 2–3) → O(n log n); a graph-layer change inside
  `DirectedKnn` for Euclidean metrics — keeps recompile fidelity, benefits every consumer.
- **Kill the O(n³):** restrict to H1 (config flip), or approximate H0-Wasserstein (paper §6: Gaussian-
  mixture / entropic / sliced OT).
- **Incremental PH via vineyards** (`RuVineyard.cs` exists): under fixed-skeleton only weights change →
  the filtration reorders → vineyard update instead of full recompute.

### P2 — Problem size
Subsampling / landmarks / witness complexes — shrinks n, helps kNN and Wasserstein superlinearly; also
the conceptual bridge to P4. Record any cap (no silent truncation).

### P3 — Iteration efficiency
Adaptive cooling, early-stop on plateau, larger accepted steps, restart ensembles. Each avoided iter
is a whole eval saved.

### P4 — Scale-out: Distributed SPRED (§3.2) — a facet of scale, its own sub-track
Partition X into m blocks, run `Spred.Compute` per block (embarrassingly parallel), aggregate the
block subspaces by **geometric median on the Grassmann** (Weiszfeld/IRLS). Distinct from P1–P3 — it
changes the *decomposition*, not per-eval speed — composes on top (each block still wants a fast eval),
and carries a **robustness co-benefit** (median breakdown point, outlier-resistant), so it is not
purely a performance lever.

**You-lineage (why it's mostly wiring).** Distributed SPRED and You's distributed PCA are the same
robust-aggregation pattern — geometric median of per-partition manifold-valued estimates — differing
only in the estimator and the factor:
- Distributed SPRED — `2106.02096` §3.2 (cites Lin et al. 2020): subspace-only → median on **Gr(n,k)**.
- Distributed PCA — `2605.20681` "Scale-Calibrated MoM for Robust Distributed PCA": (mean, subspace) →
  scale-calibrated MoM on the product **ℝᵖ × Gr(r,p)**.
- Underpinnings: product-median theory `2505.18844`, scale selection `2605.08001`, constant metric
  scaling `2601.10992`; the persistence-comparison sibling is `2208.12435`.
  Corpus: `codex-scientiae/corpora/KisungYou/`.

**Reuse — verified present, but only the primitive.** The piece distributed SPRED needs is the
Grassmann geometric median, already wired and exercised by `MoMPCA.ComputeMoM`
(`src/maths/geometry/dim-reduction/MoMPCA.cs` — it calls `GeometricMedian.Compute` on `GrassmannManifold`,
and does the full scale-calibrated product median over `ScaledManifold`×`ProductManifold`). So
distributed SPRED ≈ `Spred.Compute` per block → `GeometricMedian.Compute<Grassmann>` — **not** `MoMPCA`
itself (that carries the PCA mean factor SPRED lacks), but the same underlying primitive.

**Integration of You's program is only partially complete** (AG — started before the PH-engine wall,
resuming now that the SPRED track arrives here naturally). In place: `GrassmannManifold`,
`GeometricMedian`/Weiszfeld, `ScaledManifold`, `ProductManifold`, `MoMPCA` (the 2605.20681 consumer —
implemented, **not yet oracle-validated**). Owed / uncertain: validation of the whole median/`MoMPCA`
stack (no unit or R-oracle parity tests yet — `project_kisungyou_dr_track`), and much of the wider
cluster (the general product-median theory `2505.18844` beyond the PCA use, the scale-selection
approaches of `2605.08001`, the Wasserstein-median line `2209.03318`) not yet in code. Treat the reuse
as "lean on the primitive, validate as part of this sub-track" — not "it's done".

## Validation & measures
Per `validation-independence`: ground truth must originate outside the estimator's own model.
Engine-level manifold-opt facts are green; owed:
- End-to-end SPRED parity against an external oracle (Rdimtools / maotai are installed in the in-repo
  `r/` renv — see `project_kisungyou_dr_track`). **Parity may require a Stiefel-QR PERTURB mode** on the
  engine (add i.i.d. Gaussian to every entry + QR, per paper §3.1) to match the reference apples-to-
  apples — the Grassmann-geodesic walk is a deliberate improvement (design.md), so the oracle comparison
  needs the paper-faithful walk as an option.
- **§4 topological-equivalence measures** `μ_quasi-iso` / `μ_equiv` — post-hoc quality metrics via the
  filtration homomorphism (not the SA objective). `μ_quasi-iso` is barcode-computable (Prop 4.3: a
  height-matching sweep over `B_X`, `B_Y` with the η shift); `μ_equiv` needs π₁ of a quotient complex
  (SageMath-grade, hard) — defer. A natural reporting layer over an optimized projection.

## Application: ISOLET
SPRED is unsupervised + linear, so it is not blocked by validation-independence but is bounded by the
unsupervised ceiling (the PCA-front-end wall in `project_isolet_pca_wall`). Track SPRED vs raw-617-d and
vs PCA-14 once the driver exists.
