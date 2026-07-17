# Audit graph/PH pipeline

The audit is complete — clean across every surface. Adding the stress fact now, right after the existing determinism test:

Now verifying the new fact compiles and passes within budget — one targeted test run:

**Clean bill — no thread-safety hazards found in the parallelized pipeline.** I added the stress fact instead of fixes; it passes in ~0.4 s. Nothing committed.

## What was checked

**Static mutable state (task 1).** A regex sweep for static field declarations across `src/graphs`, `src/tda`, and all of `src/maths`, plus a second sweep for static auto-properties with setters (which hide mutable backing fields). Every hit is benign: `static readonly` immutable configs ([EuclideanMetrics.cs:10](src/graphs/distance/wrappers/EuclideanMetrics.cs:10) and siblings, `LanczosG7` lookup table), get-only properties, or expression-bodied getters returning fresh instances. There is **no non-readonly static field** anywhere in the swept surfaces, no `[ThreadStatic]` in all of `src`, and no static caches/memoization. The two singletons that warranted a look are stateless: `MetricRegistry.BoxingVisitor`/`PropertiesVisitor` ([MetricRegistry.cs:156-166](src/graphs/distance/MetricRegistry.cs:156)) hold nothing but one pure method, and `RawDistanceWeights.Instance` ([RipsFiltration.cs:24](src/tda/ph/RipsFiltration.cs:24)) is an empty immutable record.

**Per-file verification of everything inside the `Parallel.For` region:** `SubspaceAnnealer` (manifold, RNG, buffers all call-local), `PersistenceObjective` (per-block instance; `Evaluate` only reads `_data` and the projection — relevant because phase 3 of `ComputeWithDiagnostics` concurrently evaluates the _shared_ aggregate array across blocks, which is read-only there), `Pca`, `MatrixOps.Orthonormalize` (in-place on caller buffer only), `GrassmannManifold` (readonly struct, pure helpers), `Xoshiro256PlusPlus` (instance state only), `GraphMetric.FromFeatures` (read-only closures), `GraphCompiler.Build` (fully stateless — fresh log/generator/filter/scaler per call, config never mutated), `RipsFiltration`/`PersistentHomology`/`DiagramMetrics` (all `Sort()` calls in `tda/ph` operate on call-local collections; `DiagramMetrics` never sorts its barcode inputs in place).

**ArrayPool discipline.** `EigenFast` (under every PCA warm start and `ExpMap` SVD) rents and returns in `finally` with results copied out before return — no double-return, no retained pooled buffer. `PathNeighborRefiner` ([PathNeighborRefiner.cs:50-97](src/graphs/pipeline/refinement/PathNeighborRefiner.cs:50)) uses the `Parallel.For` localInit/localFinally pattern correctly: one `ThreadScratch` per worker, disposed exactly once. The `Spectral.cs`/`GraphLaplacian.cs` pool sites are off this pipeline (SPRED uses `RawDistanceWeights`, not the effective-resistance/spectral path) — noted, not audited in depth.

**Inner parallelism (task 2).** All four inner `Parallel.For` sites are _correct_, not just oversubscribed: `DirectedKnn` writes only per-`i` slots, `EpsilonBall` phase 1 writes only `halfLists[i]` with sequential symmetrization after, `LMPGraph` writes disjoint CSR row ranges of a fresh weights array, `PathNeighborRefiner` as above. Nothing to flag beyond the documented tradeoff.

## Test added

`ComputeWithDiagnostics_RepeatedParallelRuns_StayDeterministic` in [DistributedSpredTests.cs:146](tests/tda/dim-reduction/DistributedSpredTests.cs:146) — five serial-vs-parallel repetitions with distinct seeds, 3 blocks × 14 points, `maxIters: 2`, comparing aggregate projection, full objective, and all per-block projections/objectives each rep. Passed on first run: **438 ms**, far under the 30 s budget.

One incidental observation, not a hazard: [Spred.cs:35](src/tda/dim-reduction/Spred.cs:35) calls `.Projection` on `SubspaceAnnealer.Compute(...)` while the overload I read returns `double[][]` directly — there's evidently a second overload in the working tree (the build is green, so it resolves); the test file had also drifted on disk since my read, consistent with your in-flight edits.

# Deduplicate DistributedSPRED diag/evals

I'll start by reading the relevant files to understand the current structure.

Now let me look at `Spred.cs` and the test files.

**Edit 2 — `Spred.Compute` adapts to the new return type (public surface unchanged):**

**Edit 3 — `DistributedSpred`: single-block path reuses the one evaluation; `RunBlock` stops re-evaluating:**

**Edit 4 — update `SubspaceAnnealerTests` to the new signature and add the determinism tripwire:**

Now let me locate the test projects and run them.

Geometry tests pass (15/15, including the new tripwire fact). Now the TDA.DimReduction project:

Done — both findings fixed, stretch taken, changes left uncommitted in the working tree. Test results: geometry suite **15/15 green**; TDA.DimReduction **47 passed, 2 failed — both failures are from another session's concurrent, uncommitted work, not this change** (details below).

## What changed

**Finding 1 — single-block triple evaluation** ([DistributedSpred.cs:124](src/tda/dim-reduction/DistributedSpred.cs:124)): the single-block `ComputeWithDiagnostics` path now constructs the block result directly, reusing `RunBlock`'s one `LocalObjective` for `FullDataObjective` and the block's `AggregateObjective` — sound because the aggregate is the block projection, the block data is the full input, and `PersistenceObjective.Evaluate` is deterministic (I verified it holds no per-call mutable state). Three full PH evaluations outside the anneal became zero.

**Finding 2 — stretch taken**: `SubspaceAnnealer.Compute` now returns a `SubspaceAnnealerResult(double[][] Projection, double Objective)` record ([SubspaceAnnealer.cs:14](src/maths/geometry/dim-reduction/SubspaceAnnealer.cs:14)) carrying the internally tracked best value. The ripple was exactly the three sites the brief predicted: [Spred.cs:35](src/tda/dim-reduction/Spred.cs:35) appends `.Projection` (its public surface is unchanged), `DistributedSpred.RunBlock` uses `annealed.Objective` instead of re-evaluating (one full PH evaluation saved per block), and the annealer tests were updated. The regression tripwire you asked for is a new fact, `Compute_ReturnedObjective_MatchesFreshEvaluation`, asserting the returned value bit-matches a fresh evaluation at the returned projection.

## The two test failures — not from this change

While I worked, another session added uncommitted edits to this tree: two new tests in `DistributedSpredTests.cs` (both named `*_CorruptedFirstBlock_*`, stressing Weiszfeld/Grassmann-median convergence away from a corrupted block-0 initialization), +46 lines in `PersistenceObjective.cs`, a tweak to `SpredCylinderTests.cs`, and a new untracked `PersistenceObjectiveValidationTests.cs`. The two failures are exactly those two new tests, and they're independent of my change:

- `AggregateProjections_CorruptedFirstBlock_StillFindsCleanMajority` calls `AggregateProjections` directly with hand-built plane literals — it touches no code I modified (no annealer, no objective).
- `ComputeWithDiagnostics_CorruptedFirstBlock_AggregatesCleanMajority` runs with `maxIters: 0` and asserts only projection geometry, which my change doesn't alter (projections are computed identically; I only removed redundant evaluations of them).

They read like red TDD tests for a `GeometricMedian` initialization fix in progress — `GeometricMedian` itself is unmodified in the tree. Every assertion the brief called out passed: `LocalObjective == AggregateObjective` at precision 10 on the clean fixtures, parallel-vs-serial parity, and `FullDataObjective` matching the fresh global-baseline evaluation at precision 10.

One caution: with two sessions editing this tree simultaneously, our files didn't collide this time (my four files are disjoint from theirs), but that was luck — worth coordinating before the next overlapping run.

# Reconciliation ledger — tie-off (2026-07-16)

Parent-thread review of the four reports above plus the fifth chip (median corrupted-first stress —
its report was not captured here, but its artifacts are identifiable in-tree). Every chip's work
product was re-landed by the reconciliation commits; nothing from the reports remains stranded.

| Chip | Disposition |
|---|---|
| Harden block-size / error surface | Landed in `10d7939` — internal `RunBlocks` + `AggregateException` unwrap, `InvalidOperationException` block wrap with index/row count, tiny-block + body-cancellation tests. |
| Median corrupted-first stress | Tests landed (`CorruptedFirstBlock` ×2 + `Tilted` variant), including the `maxIters: 0` clarifying comments on all four aggregation fixtures. The failure they exposed was root-caused — Grassmann **cut-locus stall** at the `frames[0]` warm start — and fixed via `GeometricMedian.MedoidIndex`; `MoMPCA` re-seated on the same primitive (own brief + tests). |
| Config / input validation | Landed in `10d7939` — constructor validation + ragged-row guard, `PersistenceObjectiveValidationTests`, cylinder comment fix. |
| Dedup diagnostic evaluations | Landed in `10d7939` — `SubspaceAnnealerResult(Projection, Objective)`, `Spred`/`DistributedSpred` re-seated, single-block value reuse, determinism tripwire fact. |
| Static-state audit | Clean bill + `RepeatedParallelRuns` stress fact landed. Its one residual — the `Spectral.cs`/`GraphLaplacian.cs` pool sites, reachable in parallel via `FiltrationWeights.EffectiveResistance` — closed by direct read 2026-07-16: rent/return strictly in `try/finally`, buffers call-local, no static state; `ArrayPool.Shared` is thread-safe. |

Beyond the chips, the same window closed the ISOLET Phase-2 gate end-to-end: sliced + Sinkhorn
approximate diagram distances with the config selector (`cc5eb25`), live T4transport oracle parity
(`f7c0432`), and the S0 dim-30 pilot with the wall-clock extrapolation the gate demands (`5dc3ceb`).

Open after tie-off:

- **S0 annealer mobility** (`2c84cac`): the iteration-budget probe says S0 is blocked on the SA's
  fixed cooling/step schedule, not compute — dev-sequence P3 territory (adaptive cooling, larger
  accepted steps, restart ensembles), actively being probed in a concurrent lane.
- **Block-seed derivation** (`seed + 1009·block` aliasing): still the deliberate deferral into the
  repo-wide RNG-hygiene audit.

Process note for the record: the collisions the reports describe came from five single-item chips
sharing one tree. Per the updated dispatch guidance, future dispatches bundle coherent sequences
into fewer chips on complementary real estate, with briefs like this file as the coordination
substrate and reports appended here so the parent thread can pick the progress back up.

# Validate PersistenceObjective config/inputs

I'll read the relevant files first.

Now the existing tests and DistributedSpred's validation for message-style reference:

I have the full picture. Now let me check the smoke tests (for placement conventions) and locate the test project csproj:

All edits in place. Running the TDA.DimReduction test project to verify:

Grep says `CorruptedFirstBlock` IS in that file now, though my earlier read didn't show it — the file may have changed on disk. Let me look:

Done — validation is in, tests are green for everything this change touches, and nothing is committed.

## What changed

**[PersistenceObjective.cs](src/tda/dim-reduction/PersistenceObjective.cs)** — the constructor now rejects ragged data (`"All data rows must have the same dimension."`, matching DistributedSpred's wording, so both surfaces throw the same `ArgumentException`) and calls a new `ValidateConfig` that guards every rule: `MaxDimension >= 1`, `Dimensions` non-empty, each entry `0 <= Dim < MaxDimension` (message explains H_k needs (k+1)-simplices while the Rips complex only builds up to `MaxDimension`, so the term compares degenerate diagrams and goes constant), weights finite and > 0, `WassersteinOrder` finite and > 0, `MinPersistence >= 0` (NaN caught by the `!(x >= 0)` form), and `PathologyPenalty` finite and > 0. Since every entry point — `Spred`, `DistributedSpred`, direct use — constructs through here, they all inherit it.

Two judgment calls, reasoning stated in the messages:

- **Zero/negative weights are rejected**, not just non-finite ones. A zero weight isn't a harmless no-op: it still pays the Wasserstein matching, and under `EssentialPolicy.InfiniteOnMismatch` it turns an infinite penalty into NaN (0·∞), poisoning the objective. A negative weight rewards barcode mismatch. (`VarianceRegularizer` stays deliberately signed and unconstrained — its doc says negative rewards variance.)
- **`WassersteinOrder` must also be finite**, slightly beyond the brief's `> 0`: p = ∞ would silently break the cost exponentiation rather than give bottleneck distance.

**[PersistenceObjectiveValidationTests.cs](tests/tda/dim-reduction/PersistenceObjectiveValidationTests.cs)** (new) — 22 rejection facts, one theory/fact per rule with `ParamName` asserts, plus `Compute_RaggedInput_Throws` on `Spred.Compute` mirroring the DistributedSpred test's shape.

**[SpredCylinderTests.cs:52](tests/tda/dim-reduction/SpredCylinderTests.cs:52)** — comment reworded to "paper's multi-order idea, weights deliberately unnormalized" so the (0.5, 1.0) weights no longer read as the paper's convex (λ, 1−λ).

## Verification, with one heads-up

The full `TDA.DimReduction.Tests` run is 47/49. The 2 failures are **not from this change**: they're in `AggregateProjections_CorruptedFirstBlock_StillFindsCleanMajority` and `ComputeWithDiagnostics_CorruptedFirstBlock_AggregatesCleanMajority` — tests that were added to `DistributedSpredTests.cs` by something else while I worked (the file grew ~80 lines between my read and the test run; git shows it modified in the working tree, and I never edited it). They fail on aggregation geometry (Grassmann distance π/2 from the clean plane), which validation can't influence — my change only ever surfaces as `ArgumentException`, and their configs construct fine. A filtered run of my validation tests plus the cylinder test is 23/23 green.

Best guess on those two, per your delegate-don't-debug preference: they assert the clean majority still wins when the _corrupted_ block comes first, and the π/2 result suggests `AggregateProjections`' median aggregation is order-sensitive (e.g., seeded from block 0) — one for whoever owns that in-flight work.

# Harden DistributedSPRE block-size/error-surface

I'll harden the distributed SPRED error surface. First, let me read the relevant files.

Now let me look at `PersistenceObjective` and `Spred` to understand the failure paths, and check how internals are exposed to tests.

Let me establish a baseline build before making changes.

Build is green at baseline. Quick check of how Maths.Geometry globs its sources (to be sure which annealer is compiled):

Odd — a forced rebuild passes despite an apparent type mismatch. Let me verify the on-disk bytes directly.

Now the `RunBlocks` normalization (internal + `AggregateException` unwrap) and the `BlockSetupFailure` helper:

Now the tests. Since my earlier read of the test file may also be stale, let me verify the anchor region on disk first:

The test file on disk has two extra tests vs my stale read. Reading the current insertion region:

The region matches. Inserting the new error-surface tests after `Compute_PreCanceledToken_Throws`:

Now build and run the test project:

Same two pre-existing failures; everything else green. Confirming my six new test cases pass explicitly:

Done — the distributed SPRED error surface is hardened, all six new test cases pass, and changes are left uncommitted in the working tree.

## Changes to [DistributedSpred.cs](src/tda/dim-reduction/DistributedSpred.cs)

- **Normalized exception surface** — [RunBlocks](src/tda/dim-reduction/DistributedSpred.cs:275) is now `internal` (the test project already has `InternalsVisibleTo`) and carries a summary doc stating the contract. The parallel path catches `AggregateException` from `Parallel.For` and rethrows the first inner exception via `ExceptionDispatchInfo.Capture(...).Throw()`, so serial and parallel throw identical types with original stacks. `OperationCanceledException` for the shared token passes through untouched (both `Parallel.For`'s own cancellation collapse and the unwrap preserve it).
- **Informative small-block failures** — in both `Compute` and `ComputeWithDiagnostics`, the per-block work (objective construction + annealing) is wrapped in a `catch (Exception failure) when (failure is not OperationCanceledException)` that throws `InvalidOperationException` via a new `BlockSetupFailure` helper: _"Distributed SPRED block 3 (1 rows) failed objective construction or annealing; the block may be too small for the objective's graph recipe."_ with the original as `InnerException`. The wrap applies when `blockCount > 1`; a single block delegates to plain SPRED as before, where block context adds nothing.
- **XML docs** on `Compute` and `ComputeWithDiagnostics` now document both contracts (block-index/row-count wrap, why a validation-time floor can't exist since the recipe's K isn't visible, and cancellation never surfacing as `AggregateException`).

## Tests added to [DistributedSpredTests.cs](tests/tda/dim-reduction/DistributedSpredTests.cs:102)

- Tiny-block theories for `Compute` and `ComputeWithDiagnostics` (`blockCount == data.Length` → 1-row blocks under kNN K=6), each run serial and parallel — asserting `InvalidOperationException` in both modes with block index and row count in the message (block 0 asserted exactly for the deterministic serial case). The actual deep failure turns out to be the annealer's PCA warm start (`Pca.Compute` requires ≥2 samples), not the barcode build — a 1-node kNN graph builds fine — but the wrap is agnostic either way.
- `RunBlocks_BodyCancellation_SurfacesOperationCanceledException` — deterministic, no wall-clock timing: the body cancels the shared `CancellationTokenSource` then throws OCE with that token; asserts exact `OperationCanceledException` (never `AggregateException`) with the matching token, serial and parallel.

## Verification and two things you should know

Suite result: **55 passed, 2 failed of 57**. The two failures (`AggregateProjections_CorruptedFirstBlock_StillFindsCleanMajority`, `ComputeWithDiagnostics_CorruptedFirstBlock_AggregatesCleanMajority`) are **not from this change**: the first calls `AggregateProjections` directly (untouched by my diff), the second runs `maxIters: 0` with no exception thrown, so none of my code engages. Both fail with the Grassmann median at distance π/2 from the clean plane — i.e., Weiszfeld stuck at its corrupted block-0 initialization.

That connects to the second thing: **this working tree is being edited concurrently** by another session. Mid-task, `git status` showed modifications the session-start snapshot didn't have — `SubspaceAnnealer.cs`, `Spred.cs`, `PersistenceObjective.cs`, the SPRED test files, plus a new untracked `PersistenceObjectiveValidationTests.cs` and a new `AggregateProjections_TiltedCorruptedFirstBlock_ConvergesToCleanMajority` test. The two failures look like part of that in-flight Weiszfeld-initialization work rather than a regression. I verified every file's on-disk content immediately before editing, and a brief stash/pop I used for baseline checking restored everything intact — but you may want to be aware two sessions are touching the same SPRED files right now.
