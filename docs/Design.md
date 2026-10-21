# Design Principles

ThermoMapper was not written against a style guide, but the same design choices recur across all of its modules. This page names ten of them. For each, it explains the rationale, quotes a representative signature verbatim from the code, lists where the pattern applies, and states the trade-offs. Use it to read unfamiliar code faster, and to make new code fit.

| # | Principle | In one line |
|---|---|---|
| 1 | [Declarative config, pure engine](#1-declarative-config-pure-engine) | Immutable records describe a run; stateless engines, mostly static functions, execute it |
| 2 | [Compile-time specialization](#2-compile-time-specialization) | Struct generics and static abstract members let the JIT erase dispatch |
| 3 | [Currencies as the narrow waist](#3-currencies-as-the-narrow-waist) | A few shared data types decouple producers from consumers |
| 4 | [Honest abstention](#4-honest-abstention) | "Don't know" is a first-class answer |
| 5 | [Caller-owned buffers](#5-caller-owned-buffers) | Hot paths write into spans the caller provides |
| 6 | [Reproducible randomness](#6-reproducible-randomness) | Seeds are explicit, derived by index and checkpointable |
| 7 | [Reference/fast pairs and oracles](#7-referencefast-pairs-and-oracles) | Every fast path has a slow twin to check it |
| 8 | [Diagnostics as output](#8-diagnostics-as-output) | Decisions and pathologies are reported, not hidden |
| 9 | [Kernel/shell split and decorators](#9-kernelshell-split-and-decorators) | Algorithms are small kernels; orchestration wraps them |
| 10 | [Capability = nullability](#10-capability--nullability) | An optional array's presence *is* the feature flag |

## 1. Declarative config, pure engine

**Rationale.** When a run is described by an immutable record, it can be serialized, compared, logged next to its results and replayed. When the engine is a static function of (config, data), there is no hidden state to reset between runs, and concurrent runs cannot interfere. Validation happens once, at the engine's entry, so every caller gets the same guarantees. Several config types note in their source comments that fluent builders belong at the command-line or REPL boundary, not in the engine.

```csharp
public static GraphBuildResult Build(GraphCompilerConfig config, int n, GraphMetric metric, ProtectedEdgeSource? protectedEdges = null)
public static BarsResult Run(BarsConfig config, double[] x, double[] y, double[] grid)
```

**Where it applies.**

- `GraphCompilerConfig` → `GraphCompiler.Build` ([Graph Construction](./GraphConstruction.md))
- `BarsConfig` → `Bars.Run` ([Bayesian Regression](./BayesianRegression.md))
- `IrlsOptions` → `Irls.Solve` ([Riemannian Geometry](./RiemannianGeometry.md))
- `PersistenceObjectiveConfig` → `PersistenceObjective` ([SPRED](./Spred.md))
- `AccumulationSpec` for Swendsen–Wang runs ([SPC Engine](./SpcEngine.md)), and `BlattPartitionConfig` → `BlattPartitionStrategy` ([SPC Analysis](./SpcAnalysis.md))
- `SpcRunLayout` → `SpcRunPaths` ([SPC Analysis](./SpcAnalysis.md))

**Trade-offs.**

- Configs resolved at run time from earlier stages' diagnostics, such as graph construction's `Auto` choices, fail late: an invalid combination can surface only after several stages have run.
- Records hold arrays by reference, so "immutable" is often shallow. Mutable collections inside init-only records can still change.

## 2. Compile-time specialization

**Rationale.** Distance metrics, manifolds, robust losses and sampler configs sit inside the innermost loops. Making them `struct` type parameters, with `static abstract` flags where a compile-time switch is needed (`IsFlat`, `IRobustLoss`, `ISwConfig`), lets the JIT compile one specialized loop per concrete type, with no virtual calls and no dead branches. A **visitor bridge** turns a runtime choice, such as a metric spec string, into a compile-time type exactly once.

```csharp
static abstract bool IsFlat { get; }
TResult Visit<TMetric>(TMetric metric) where TMetric : struct, IDistanceMetric
public static TResult Invoke<TResult>(MetricSpec spec, IMetricVisitor<TResult> visitor)
```

**Where it applies.**

- `IRiemannianManifold.IsFlat`: flat manifolds erase the curved-space branches.
- `IRobustLoss`: `IsClosedForm`, `IsSingularAtZero` and `Weight`.
- `ISwConfig`: which per-edge counts the Swendsen–Wang hot loop writes.
- `TMetric : struct, IDistanceMetric` in `Prim`, `CoreDistances` and `HdbscanRunner.Run<TMetric>`.
- `MetricRegistry.Invoke` → `IMetricVisitor<TResult>`, used by HDBSCAN's metric dispatch.

See [Distance Metrics](./DistanceMetrics.md#design) and [Riemannian Geometry](./RiemannianGeometry.md).

**Trade-offs.**

- Cold callers use a **boxed** path (`MetricRegistry.Create`) rather than a visitor.
- Adding a metric touches four places: the struct, an enum value, a registry row, and the parse and invoke cases.
- Generic signatures become long, and errors surface as constraint violations.

## 3. Currencies as the narrow waist

**Rationale.** A handful of shared types carry results between modules:

- the CSR graph and its slot-indexed per-edge fields;
- the three edge currencies;
- node landscapes;
- `Assignment` and `Groups`;
- `Dendrogram` and its cost axis;
- `Barcode`;
- `NeighborSelection`, `EigenPair` and `ZigzagFiltration`.

Producers and consumers depend only on these types, never on each other. That is what lets any sampler feed any cut, any tree feed any selector, and any persistence engine feed any diagram distance. The full definitions are in [Concepts](./Concepts.md).

```csharp
public interface IGraphSignal<in TCurrency>
double[] Compute(TCurrency currency, CsrGraph graph)
```

**Where it applies.** Every module boundary. `IGraphSignal` is contravariant, so a signal written for a currency accepts that currency from any producer. `TDA.Ph`'s nerve types store only graphs and member arrays, so persistence code does not depend on Mapper.

**Trade-offs.**

- Semantics that live in the data rather than the type are checked only at run time, or not at all:
  - whether a graph's weights are distances or couplings (`EdgeWeightKind` is informational);
  - which slot of an edge is canonical;
  - which cost-axis string a tree carries.
- The cost-axis law compares strings ordinally, so a typo disables a walk.

## 4. Honest abstention

**Rationale.** A clustering method that is unsure should say so. Forcing every point into a cluster hides uncertainty and inflates agreement scores. ThermoMapper generalizes HDBSCAN's "noise" label to every method.

```csharp
public const int Unassigned = -1
```

**Where it applies.**

- `Assignment.Unassigned`. Resolvers and landscape walks return −1 rather than guess.
- `Assignment.Coverage` reports the assigned fraction separately from quality.
- Evaluators score only the assigned subset.
- Periphery-completion policies may fill in unassigned points but never relabel assigned ones.
- Diagnostics use NaN or −1 sentinels instead of fabricating values.

See [Cluster Resolution](./ClusterResolution.md) and [Cluster Evaluation](./ClusterEvaluation.md).

**Trade-offs.** Every consumer must handle −1 and NaN. `Groups.Argmax()` is the explicit exception: it never abstains. Scores computed on the assigned subset must always be read together with coverage.

## 5. Caller-owned buffers

**Rationale.** In hot paths, allocation costs more than arithmetic. Kernels that write into caller-supplied spans can reuse per-thread buffers under `Parallel.For`, can be called in tight loops without GC pressure, and make the memory cost of each call explicit. Factorizations are computed once and then applied many times.

```csharp
public static void DecomposeSymmetricInPlace(ReadOnlySpan<double> flatMatrix, int n, Span<double> outEigenvalues, Span<double> outEigenvectors, int maxSweeps = 256, double tol = 1e-12, DenseEigenFastVariant fastVariant = DenseEigenFastVariant.Default)
public static void ComputeMutualReachabilityMst<TMetric>(ReadOnlySpan<double> data, int n, int dim, ReadOnlySpan<double> coreDistances, Span<bool> visited, Span<double> minWeight, Span<int> parent, TMetric metric) where TMetric : struct, IDistanceMetric
```

**Where it applies.**

- `EigenFast`, `Prim` and `Dijkstra` work in caller spans (`Dijkstra` also takes a caller-supplied priority queue). `CoreDistances` writes into a caller span but allocates one small buffer per call.
- Tiered scratch memory in the metrics: `stackalloc` for small inputs, a per-thread buffer for medium ones, `ArrayPool` for large ones.
- Band Cholesky/LDLᵀ and dense Cholesky: factor once, solve many times.
- HDBSCAN reuses its scratch across runs.
- The CSR `UndirectedEdges()` walk allocates nothing.

See [Linear Algebra](./LinearAlgebra.md) and [Graph Substrate](./GraphSubstrate.md).

**Trade-offs.**

- Signatures grow long.
- Buffer-size preconditions become part of the contract.
- Not every path follows the rule. For example, `CoreDistances` allocates one buffer per call.

## 6. Reproducible randomness

**Rationale.** Scientific runs must be repeatable, including runs parallelized across replicas and temperatures. So:

- generators are passed explicitly;
- seeds are derived by *index*, so results do not depend on worker count or scheduling;
- generator state can be checkpointed and resumed exactly.

Independent random processes get independent streams, so changing one never shifts the other.

```csharp
public static int[] Derive(int master, int count)
public static int? Derive(int? baseSeed, double temperature, int replica, int round = 0, int quantizationFactor = DefaultTemperatureQuantizationFactor)
```

The first line is `SeedTree.Derive`. The second is the SPC per-task seed, which folds in the temperature and replica.

**Where it applies.**

- `Xoshiro256PlusPlus` with SplitMix64 seed expansion, and `SeedTree` ([Probability and Information](./ProbabilityAndInformation.md)).
- The ensemble shell derives 2C seeds, split into kernel streams and readout streams, using `SeedTree`'s prefix property ([Samplers](./Samplers.md)).
- SPC tasks are seeded by (seed, T, replica, round), and checkpoints store the generator's four state words ([SPC Engine](./SpcEngine.md)).
- Distributed SPRED seeds each block with `SeedTree` ([SPRED](./Spred.md)).
- The eye dataset gives its noise field a separate stream ([Synthetic Datasets](./SyntheticDatasets.md)).

**Trade-offs and exceptions.**

- A `null` seed means OS entropy and is not reproducible.
- Derived 32-bit child seeds can collide (birthday bound).
- Temperature seeds are quantized at $10^{-3}$.
- `Synthetic` generators, k-means++, the GMM initializers, the LOBPCG starting block, FastICA's initialization and the dense Cholesky sampler use `System.Random` instead.

## 7. Reference/fast pairs and oracles

**Rationale.** Fast algorithms (Jacobi with FMA, HDT dynamic connectivity, coned zigzag, cohomology with clearing) are subtle. Keeping a slow, obviously correct twin in the same codebase makes the fast path testable. Where the mathematics admits several independent algorithms, the code keeps a portfolio of **oracles**.

```csharp
public static Barcode Compute(ZigzagFiltration f, int maxDimension = 0, bool representatives = false, bool strictRepresentatives = false, GraphZigzagAlgorithm algorithm = GraphZigzagAlgorithm.Fast)
```

**Where it applies.**

- **`EigenFast` vs `Eigen`.** Selected at compile time with `EIGEN_REFERENCE` ([Linear Algebra](./LinearAlgebra.md)).
- **`GraphZigzag` Reference vs Fast.** A dense barcode forest vs a sparse forest backed by HDT.
- **Zigzag oracles.** The naive generalized-rank oracle (`ZigzagBarcodeNaive`), `ZigzagMapBarcode` and the strong-collapse references cross-check `FastZigzag` and the graph engines; `ReflectionZigzag` is a partial alternative engine ([Zigzag Persistence](./ZigzagPersistence.md)).
- **Standard persistence.** `PersistentHomology`'s plain column reduction is the reference for the cohomology and involuted fast paths ([Persistent Homology](./PersistentHomology.md)).
- **SPC.** The deterministic PKWang direct solver sits beside Monte Carlo Swendsen–Wang as an alternative producer of `Affinities`. Its temperature axis is in cumulative-energy units, so it is not a numerical reference for Swendsen–Wang ([SPC Engine](./SpcEngine.md)).

**Trade-offs.**

- Oracles can be dangerously slow at default settings. The naive zigzag oracle, for example, needs an explicit dimension to terminate.
- Pairs can drift apart. Documentation of fast-path complexity sometimes overstates what the implementation delivers.
- Reference and fast engines may emit bars in different orders, so compare barcodes as multisets.

## 8. Diagnostics as output

**Rationale.** A pipeline that silently picks a bandwidth, repairs a disconnected graph or falls back to singletons is hard to trust. ThermoMapper reports what it decided, as data, and lets configured thresholds turn pathologies into errors.

```csharp
public sealed record DiagnosticMessage(DiagnosticLevel Level, string Stage, string Text)
public sealed class GraphPathologyException : Exception
```

**Where it applies.**

- The `DiagnosticsLog` returned by every graph build.
- Pathology interrupts that throw `GraphPathologyException`.
- `GraphHealth`, `Hubness` and `MstBridge` ([Graph Construction](./GraphConstruction.md)).
- Report records with a `Computed` flag or NaN/−1 sentinels instead of exceptions ([Graph Substrate](./GraphSubstrate.md)).
- Mapper's coverage, topology and warnings ([Mapper](./Mapper.md)).
- SPRED turns a pathology into a finite objective penalty ([SPRED](./Spred.md)).

**Trade-offs.**

- Sentinels overload meaning. A girth of −1 can mean "acyclic" or "too large to compute".
- Warnings are strings, not codes.
- A report must be read to be useful, and nothing forces the caller to read it.

## 9. Kernel/shell split and decorators

**Rationale.** A kernel is small and knows only its own mathematics: one Markov-chain step, one kernel evaluation, one scaler. Orchestration concerns live in a *shell*: replicas, burn-in, convergence checks, parallelism, retries. Behavior changes are layered on by decorators rather than by modifying kernels. The same kernel then serves many workflows.

```csharp
public interface IChain<out TDraw>
public (EnsembleRun Summary, IReadOnlyList<IChainRun<TDraw>> Chains) Run<TDraw>(IEnsembleModel<TDraw> model, int chains, int masterSeed, int burn, int samples, double rHatTarget = 0.0, int batchSize = 0)
public ScheduledMove(IRjMove<TState> inner, Func<TState, double> weight)
```

**Where it applies.**

- `IChain` (one replica) vs `ChainEnsemble` (the shell).
- `ScheduledMove` wraps a move ([Samplers](./Samplers.md)). Tempering is not a target decorator: a chain samples its target at a `Temperature`, a query on each state's evaluation, because a collapsed model's tempered likelihood is not β times its log density.
- In graph construction, local-mutual-proximity scalers decorate the global-bandwidth scaler ([Graph Construction](./GraphConstruction.md)).
- In SPC, the sweep strategy orchestrates per-temperature Swendsen–Wang kernels through an executor ([SPC Engine](./SpcEngine.md)).

**Trade-offs.**

- Behavior is spread across layers, so finding which wrapper changed a value takes navigation.
- Not every sampler uses the shell. BAPS, for example, runs its own loop.

## 10. Capability = nullability

**Rationale.** Optional statistics cost memory and time in hot loops. Instead of a capability enum plus arrays that might be stale, an optional result is a **nullable array**, and its presence is the capability. The declarative spec says what to collect, and the record carries exactly what was collected.

```csharp
public int[]? BondFormedCount { get; init; }
public int[]? SpinAgreementCount { get; init; }
public int[]? CoMembershipCount { get; init; }
```

**Where it applies.**

- The Swendsen–Wang `Accumulator`'s per-edge counts and per-node sums, driven by `AccumulationSpec`. The SPCX checkpoint format writes a presence byte before each optional array ([SPC Engine](./SpcEngine.md)).
- `SweepResult`'s chosen currencies.
- `Bar`'s optional `Generator`, `Cocycle` and `Cycle` payloads.
- `GraphBuildResult`'s optional pre-repair graph and diagnostics.

**Trade-offs.**

- Every consumer null-checks.
- The spec and the actual presence can disagree. `SweepResult`'s documentation, for example, ties the presence of `ChosenAlignments` and `ChosenCoMembership` to the sweep's accumulation flags, but the default equilibrium pass collects every currency regardless.
- Downstream resolvers state their required capabilities in comments, such as "requires `AccumulationSpec.CoMembership`", rather than in types.

## Cross-cutting trade-offs

- **Duplicated rules instead of shared helpers.** The dense and matrix-free Laplacian builders deliberately repeat the same weight-clamp rules so that both routes see the same matrix. Tabular output repeats the atomic-write protocol because its project does not reference `Archivory`. Duplication keeps dependencies narrow, but the copies must be kept in step by hand.
- **Code over comments.** Doc comments sometimes describe intended rather than implemented behavior. The module docs describe what the code does, and point out mismatches where they matter.

**See also:** [Concepts](./Concepts.md) · [Architecture](./Architecture.md) · [Glossary](./Glossary.md)
