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

### Cylinder validation (paper Fig.1) ✔ 2026-07-03 — search works + full-persistence finding
`tests/tda/dim-reduction/SpredCylinderTests.cs`. On the noisy cylinder S¹×[-2,2] (β₁=1) where PCA's
warm start flattens the loop, the anneal converges to a projection **far better than any axis-aligned
view** (objective ~0.84 vs (h,x)=4.14, (x,y)=4.91) — the SA search is functional. **Key finding:** the
H0+H1 objective preserves *full* persistence, so the naive circle (x,y) is **not** optimal — it
collapses every height onto the same θ, destroying the H0 merge structure, and scores *worse* than the
loop-flat (h,x); SPRED instead finds a faithful **oblique** view keeping loop + height. **Perf wall
confirmed:** ~4 min at n=100, ~23 s at n=70 — H0-Wasserstein O(n³) is the binding constraint (→ P1).

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
→ `RipsFiltration.GraphRips(g, filtrationWeights, maxDim)` → `PersistentHomology.Compute(g, maxDim)`
→ `Σ_dim w_dim · DiagramMetrics.Wasserstein(B_P, B_ref, dim, p, essential)`.
`B_ref` built once the same way on ambient `GraphMetric.FromFeatures(X)` under the same recipe.

**Skeleton model — RESOLVED: recompile per proposal.** The compiler's stages are distance-driven (the
kNN edge *set* moves with P), so full-recipe fidelity ⇒ rebuild the graph each proposal. Correctness-
first, and it honors the "leverage the compiler" steer. SA tolerates the combinatorial discontinuity
(global stochastic optimizer, no continuity needed). The earlier fixed-skeleton + reflow-weights idea
is demoted to a step-3a perf mode (valid only for RawDistance on a frozen ambient topology).

**Filtration-value axis (also injected):** `GraphRips` takes a `FiltrationWeights` — `RawDistance`
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

### P0 — Profile first ✔ 2026-07-03 (`tests/tda/dim-reduction/SpredProfileTests.cs`)
Per-`Evaluate` breakdown on the cylinder (JIT-warmed, 3 reps), ms:

| n | #H0 | #H1 | graph | rips | ph | W(H0) | W(H1) | total |
|---|---|---|---|---|---|---|---|---|
| 60 | 60 | 238 | 0.8 | 1.3 | 3.5 | 8.1 | 68.2 | 82 |
| 100 | 100 | 394 | 1.6 | 3.0 | 5.2 | 29.3 | 105.2 | 144 |
| 150 | 150 | 584 | 1.3 | 4.6 | 7.2 | 77.0 | 105.6 | 196 |
| 200 | 200 | 785 | 1.5 | 2.8 | 3.5 | 165.3 | 242.4 | 416 |

**Decisive finding — overturns the a-priori plan.** The bottleneck is **diagram-Wasserstein (Hungarian
matching), overwhelmingly H1 (54–83%)** — because the graph-restricted Rips emits **~4n H1 bars** (noise
loops in the loopy kNN graph) vs n H0 bars, so the H1 Hungarian O((4n)³) dwarfs everything. **graph /
Rips / PH are negligible (<15 ms).** So kNN / fixed-skeleton / KD-tree would save almost nothing — the
earlier "H0-Wasserstein / kNN is the wall" guess was wrong on every count.

### P1 — Per-evaluation speed (reordered by P0: the cost is Wasserstein, not graph/PH)
- **Prune low-persistence bars before matching ✔ 2026-07-03 — the top lever, landed.** `MinPersistence`
  on `PersistenceObjectiveConfig`, applied in `PersistenceObjective.BarcodeFor` to both barcodes.
  Measured (profiler): the ~4n finite H1 bars are **all near-zero persistence** (`maxPersH1 ≈ 0` — the
  real loop is an *essential*/∞ bar), so any τ>0 cleanly separates noise from signal. W(H1) collapses
  **~1774–3600×** (→ 0.04–0.18 ms), leaving exactly the 1 essential loop bar — **near-exact** (cylinder
  objective 0.834 pruned vs 0.841 exact) and a **denoiser**. Whole cylinder eval ~2–3× faster (31→11 s);
  **H0-Wasserstein is now the residual** (next sub-lever).
- **Approximate the Wasserstein** (paper §6: entropic/Sinkhorn or sliced OT — `T4transport`'s
  `dist_sinkhorn`/`dist_swdist`, both reference and oracle). Helps H0 and any residual H1. H0 still
  can't be dropped (drives the SA descent — cylinder converged *with* H0+H1).
- **Sparser / less-loopy filtration** — fewer H1 noise loops at the source (fill more triangles,
  mutual-kNN, or a persistence-aware graph). Secondary; pruning is simpler and downstream-agnostic.
- **Deprioritized (P0 says negligible):** fixed-skeleton, KD-tree kNN, vineyards — graph + PH are
  <15 ms/eval; optimizing them saves ~nothing until the Wasserstein matching is fixed.

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

**First wiring slice landed 2026-07-06:** `DistributedSpred.Compute` now performs contiguous block
splitting, runs `Spred.Compute` per block, converts each k×d projection to the d×k Grassmann frame
currency, and aggregates by `GeometricMedian.Compute<GrassmannManifold>`. Tests cover the median
aggregation path, a clean-majority/corrupted-block robustness fixture, and public single-/multi-block
facade smoke. `ComputeWithDiagnostics` exposes block ranges, derived block seeds, per-block
projections, and the aggregate projection. A deterministic end-to-end corrupted-block fixture now
uses that surface to show three clean circle blocks outvote two corrupted block subspaces; larger
stochastic distributed robustness/performance fixtures remain future P4 work.

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
implemented, only partially oracle-validated). `mom_oracle.R` now validates the Grassmann geometric
median primitive against `Riemann::riem.median`; owed / uncertain: validation of the whole
`MoMPCA` product/scale stack (`project_kisungyou_dr_track`), and much of the wider
cluster (the general product-median theory `2505.18844` beyond the PCA use, the scale-selection
approaches of `2605.08001`, the Wasserstein-median line `2209.03318`) not yet in code. Treat the reuse
as "lean on the primitive, validate as part of this sub-track" — not "it's done".

## Validation & measures
Per `validation-independence`: ground truth must originate outside the estimator's own model. The
`r/` sub-project is that channel — a recently-introduced R-oracle harness (C# `ROracle` bridge; PCA
parity vs base-R `prcomp` already green = precedent) validating the C# ports against **Kisung You's
own reference code**, not just paper transcriptions. Sources: the CRAN packages in the `r/` renv
**and** the vendored source clones at `D:\aghado01\packages\kisungyou` (deeper than the CRAN API —
unexported internals + MATLAB per-paper toolboxes). Coverage is **only partially mapped** (harness is
new; the K.You integration paused at the PH-engine wall and is resuming now — AG).

**Independent-eval tiers (by what runs), pick the highest available per dependency:**
1. **Live R oracle** via `r/` — `Riemann` / `SHT` / `Rdimtools` source that executes, for numerical
   parity (strongest, most independent).
2. **Transcribed fixtures** where a reference exists but can't run in-repo — notably the `papers/`
   **MATLAB** toolboxes (no MATLAB license in ThermoMapper): transcribe the algorithm / expected values
   into C# facts as a best-effort independent eval.
3. **Paper-derived conditions** where no runnable or transcribable reference exists (SPRED itself):
   glean testable invariants from the text — Betti numbers, monotonicity, analytic special cases, the
   cylinder/iris outcomes. The general paper-adaptation fallback, and often sufficient on its own.

**Oracle-source map** — vendored set is growing (~12 repos now, toward the ~37 K.You has; a **project-
wide** oracle asset). SPRED-relevant, confirmed by inspection:
- **`TDAkit/`** (R) — the key new find: `homology_diagRips.R` (Rips → persistence diagram) is the oracle
  for the C# **PH pipeline** (`GraphRips` / `FullRips` → `PersistentHomology` → `Barcode`); `summaries_dist*.R` for
  **diagram distances** cross-checks `DiagramMetrics`. SPRED's objective is now component-oracle-able.
- `Riemann/` (R) — `inference_median.R` + `special_grassmann_*`: the **Grassmann geometric median** →
  the MoM / distributed-SPRED aggregation primitive.
- `T4transport/` (R) — computational OT: `dist_wasserstein` / `dist_sinkhorn` / `dist_swdist` (exact +
  entropic + sliced) → the OT layer and the §6 Wasserstein *approximations*; `free_median_*` (`2209.03318`).
- `SHT/` (R) — hypothesis testing → MxPbf (`2112.02580`); `Rdimtools/` (R) — linear DR (PCA oracle green).
- `DirectMedian/` — **Python** (not R): free-support Wasserstein medians; peripheral to SPRED (OT/median
  track), live-via-python or transcribe.
- Other-track clones: `T4cluster`, `mclustcomp` (clustering), `dglearn` (diffusion geometry),
  `IntrinsicESS` (`2605.03266`), `GeoInfoDec`.
- `papers/` — **MATLAB** (76 `.m`; `01-SPDtoolbox`…`06-Wasserstein-Heterogeneity`, `03-TopLSM` ≈
  `2208.12435`): **transcribe-only** (no MATLAB license), Tier 2. Not fully mapped.

**SPRED-the-composition still has no direct reference** (a `spred` grep across the clones is empty), but
each dependency-half now has a live R oracle — **`TDAkit`** (PH + diagram distance) and **`Riemann`**
(Grassmann median) — so **component-wise Tier-1 validation is strong** without a SPRED-direct oracle.
The paper's qualitative examples (cylinder β₁=1, iris) stay the end-to-end sanity. Any *direct* SPRED
comparison later needs the **Stiefel-QR PERTURB mode** (paper §3.1) for apples-to-apples — the
Grassmann-geodesic walk is a deliberate improvement (design.md).

Owed R oracles (per `project_kisungyou_dr_track`): extend `mom_oracle.R` to the full MoMPCA §5.3 α
product/scale median; add `mxpbf_oracle.R` (SHT / transcribe §2.2/§3.1). Plus the
**§4 topological-equivalence measures**
`μ_quasi-iso` / `μ_equiv` — post-hoc quality metrics via the filtration homomorphism (not the SA
objective): `μ_quasi-iso` is barcode-computable (Prop 4.3, height-matching sweep with the η shift);
`μ_equiv` needs π₁ of a quotient (SageMath-grade) — defer.

### Oracle-harness status (2026-07-05) — R-env repaired; TDA PH oracle green
**Where this fits:** the first Tier-1 live oracle — validating SPRED's **PH half** (`RipsFiltration` +
`PersistentHomology`) against gold-standard **Ripser**. Built:
- `r/oracles/tda_oracle.R` — emits Ripser's full-Rips diagram (via `TDAstats::calculate_homology`, the
  engine `TDAkit::diagRips` wraps); essential/infinite bars → `-1` sentinel.
- `tests/oracle/TdaParityTests.cs` — C# builds a full Rips via `FullRips.Build` → PH (a true full Rips,
  apples-to-apples with Ripser) and matches finite H0/H1 bars, plus asserts the one
  essential H0 bar Ripser omits. Gated on `ROracle.IsAvailable`. `TDA.Ph` ref added to the oracle csproj.
- **R env repaired (2026-07-06):** wiped stale `r/renv/library` + sandbox, rebuilt the package library at
  the ThermoMapper path, upgraded PDenv R to **4.6.1**, added `r/DESCRIPTION` as the explicit oracle
  dependency manifest, and pinned the full oracle dependency closure in `r/renv.lock`.
  `r/scripts/oracle-ci.ps1` now resolves the user-level PDenv R toolchain without `$PORTABLE_ROOT`, runs
  Rscript with timeouts, and passes.

**Migration cause fixed.** The old package dirs had missing `DESCRIPTION` files after the move; the fresh
library now has zero missing `DESCRIPTION` files. `renv::load()` was also hanging in renv's Windows
sandbox/junction setup, so `r/.Renviron` disables `RENV_CONFIG_SANDBOX_ENABLED` for this isolated oracle
project. Plain activation, `renv::status()`, and the PCA smoke are green.

**TDA oracle green with a contained Windows quirk.** `TDAstats::calculate_homology` computes and writes the
expected diagram, then `Rscript.exe` exits with Windows access violation `-1073741819` ("memory could not
be read"). Upgrading to R 4.6.1 did not remove that upstream native-exit quirk, and a source rebuild with
Rtools45 reaches link but fails during package lazy-loading. The C# bridge now suppresses the Windows fault
dialog for child R processes and tolerates this exact exit code only for `tda_oracle.R` when parseable JSON
was produced. `TdaParityTests.FullRips_PH_Matches_Ripser_H0_H1` is unparked and green; the comparison
filters diagonal zero-persistence H1 intervals because Ripser omits them while the explicit reducer
materializes them.

### Full Vietoris–Rips → `src` ✔ P1 landed 2026-07-06 — TDA.Ph capability
Canonical roadmap: `issues/ph/full-rips-roadmap.md`. Emerged from the parity test:
`TdaParityTests.FullRipsBarcode` had built a *complete* distance graph and handed it to the existing
`GraphRips(maxDim=2)` — a full Rips only in **density** (all pairs vs a kNN skeleton), still the same
2-skeleton (H0/H1), and ad-hoc in the test. P1 moved that into `TDA.Ph.FullRips.Build`, a
threshold-bounded complete Euclidean density API that reuses `GraphRips`. The graph-restricted path
(`GraphRips` / `LazyRipsFiltration` / `FlagComplex`, all triangle-capped) remains load-bearing for
SPRED / `ConditionedFiltration` / `H1CycleEdges`.

**Two orthogonal axes** (name them so they don't conflate):
- *Density* — complete graph (all pairs, threshold-bounded) vs kNN-restricted skeleton. The immediate one.
- *Dimension* — 2-skeleton (H0/H1) vs higher (H2 voids …). Deferred; it's an **enumeration** gap only —
  `FlagComplex` does triangles, but `PersistentHomology.Compute` is already dimension-agnostic.

**Naming settled:** `FullRips` (dense/complete — the standard TDA term) and `GraphRips`
(the graph-restricted builder). Avoid `SparseRips` — it collides with Sheehy's sparse-Rips
*approximation*; ours is exact on the skeleton, not an approximation. `DenseRips` isn't a standard term.

**Additive, not a swap:** production stays graph-restricted (SPRED's recompile-a-sparse-graph design is
load-bearing). Full Rips is for validation, small-cloud exact topology, and capability.

**Phased plan:**
- *P0 (prereq) ✔ 2026-07-05/06:* fix the renv env, unskip `TdaParityTests`, get it green.
- *P1 (ball-rolling) ✔ 2026-07-06:* migrate the test's full-Rips minimally into `src` — a first-class,
  **threshold-bounded** `FullRips` reusing `GraphRips`; refactor the parity test onto it.
- *P2 (full integration):* metric/enclosing-radius options; the dimension axis (k-clique enumeration in
  `FlagComplex` + higher simplices — deferred); performance (threshold bounding + route dense reductions
  through `PersistentCohomology` / `PersistenceClearing`, not the naive standard reducer — the n=30
  dense-complex cost, which is separate from the parked R-env hang).

## Application: ISOLET
SPRED is unsupervised + linear, so it is not blocked by validation-independence but is bounded by the
unsupervised ceiling (the PCA-front-end wall in `project_isolet_pca_wall`). Track SPRED vs raw-617-d and
vs PCA-14 once the driver exists.
