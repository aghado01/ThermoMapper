# ThermoMapper HPC and `hashish` review from the Doccer excavation

Date: 2026-08-06
Status: auxiliary implementation report; no source change is implied
Origin: cross-project review performed while planning the next Doccer expansion

## Purpose and boundary

This report preserves ThermoMapper-facing findings that should not live in the
codex-scientiae/Doccer decision canon. The Doccer documents use the same source reading to decide
which capabilities and engineering patterns are worth lifting. This report instead records what
is already strong in ThermoMapper, what appears incomplete or misleading, and what should be
checked when ThermoMapper development returns to these areas.

The review is not a proposal to make ThermoMapper follow Doccer's contracts or sequencing. The two
projects own their own semantics. A donor implementation can simultaneously contain a valuable
algorithmic idea and need repair in its native repository.

No benchmark was run and no compatibility claim was tested in this review. Performance remarks
below identify plausible costs and benchmark targets, not measured rankings. Correctness concerns
marked **audit** should be confirmed with a minimal oracle or published test vectors before source
changes land.

## Coverage

The review read all 22 files under [`src/hashish`](../src/hashish/) and sampled the broader
performance repertoire through:

- [`CsrGraph`](../src/graphs/primitives/CsrGraph.cs),
  [`UndirectedEdgeWalk`](../src/graphs/primitives/UndirectedEdgeWalk.cs),
  [`UnionFind`](../src/graphs/primitives/UnionFind.cs),
  [`Dijkstra`](../src/graphs/primitives/traversal/Dijkstra.cs), and
  [`PathNeighborRefiner`](../src/graphs/pipeline/refinement/PathNeighborRefiner.cs);
- [`GraphLaplacian`](../src/graphs/spectral/GraphLaplacian.cs),
  [`CoherentField`](../src/graphs/spectral/CoherentField.cs),
  [`MatrixOps`](../src/maths/linalg/MatrixOps.cs), and
  [`EigenFast`](../src/maths/linalg/EigenFast.cs);
- [`EarthMover`](../src/maths/distance/EarthMover.cs),
  [`OnlineMahalanobis`](../src/maths/linalg/WelfordMahal.cs),
  [`ScatterAccumulator`](../src/maths/geometry/estimators/intrinsic/ScatterAccumulator.cs), and
  [`GaussianMixtureModel`](../src/clustering/statistical/gmm/GaussianMixtureModel.cs);
- [`SeedTree`](../src/maths/rng/SeedTree.cs) and
  [`Xoshiro256PlusPlus`](../src/maths/rng/Xoshiro256PlusPlus.cs); and
- the persistent-homology reference/fast split, including
  [`PersistenceClearing`](../src/tda/ph/PersistenceClearing.cs),
  [`LazyRipsFiltration`](../src/tda/ph/LazyRipsFiltration.cs),
  [`FastZigzag`](../src/tda/ph/FastZigzag.cs),
  [`GraphZigzagFast`](../src/tda/ph/GraphZigzagFast.cs),
  [`GraphZigzagH1Fast`](../src/tda/ph/GraphZigzagH1Fast.cs), and the structures under
  [`src/tda/ph/dynamic`](../src/tda/ph/dynamic/).

## Strong patterns worth retaining

| Pattern | Good examples | Why it is valuable | Maintenance posture |
|---|---|---|---|
| Span-first kernels with caller-owned destinations | `Histogram.Normalize`, `GraphLaplacian.BuildDenseColumnMajor`, `Dijkstra.ComputeBoundedDistances` | Makes allocation and layout explicit while allowing array-backed convenience APIs. | Prefer this as the primitive shape; keep allocating wrappers at the shell. |
| Count → prefix sum → exact allocation → fill | `CsrGraph.FromEdges`, `CsrGraph.InducedSubgraph` | Converts irregular input into compact contiguous storage deterministically. | Reuse the construction shape for sparse matrices and indexes; add overflow and input validation. |
| Flat, layout-declared buffers | CSR targets/weights, flat TF-IDF rows, column-major Laplacian, parallel-array top-K heap | Improves locality and makes downstream SIMD/library calls possible without copying. | Keep layout in the contract of the operation or artifact; do not leave row/column order implicit. |
| Stack/pool hybrid scratch | `Levenshtein`, `MinHash`, `GaussianMixtureModel`, several graph/maths kernels | Avoids heap churn for small work while bounding large temporary allocations. | Centralize thresholds only after measurement; dynamic `stackalloc` still needs a hard safe ceiling. |
| Per-worker reusable scratch | `PathNeighborRefiner.ThreadScratch` + `Dijkstra` | Amortizes arrays and priority queues across a parallel partition. | Preserve worker ownership; consider touched-index/stamped reset where full `O(N)` clearing dominates sparse work. |
| Work proportional to the result set | masked Dijkstra early exit, upper-triangle edge walks, sorted-row intersection, bounded top-K heap | Avoids dense scans or full sorts when the requested population is sparse or bounded. | Keep a clear reference path and define deterministic ties. |
| Hardware-tiered kernel plus scalar tail | `MatrixOps`, `CoherentField`, `EigenFast`; `TensorPrimitives` use in cosine/TF-IDF | Lets the runtime use wide lanes without making vector width part of semantics. | Differential-test every tier; benchmark manual intrinsics against `TensorPrimitives` and the JIT before keeping duplicate machinery. |
| Resettable state and output spans | `UnionFind.Reset`, `UnionFind.WriteRootSizesTo`, online accumulators | Reuses long-lived storage and lets a caller own result materialization. | Expose snapshots separately from hot operations; remove unused retained scratch. |
| Stable online numerics | Welford updates, log-sum-exp in the GMM E-step | Handles streams and high-dynamic-range arithmetic without a batch materialization or avoidable instability. | Add merge/combine forms where parallel reduction is expected and test extreme ranges. |
| Structural RNG derivation and checkpointable streams | `SeedTree`, `Xoshiro256PlusPlus`, ensemble fan-out | Makes parallel runs reproducible without alias-prone arithmetic seed offsets. | Keep algorithm/state identity in persisted manifests; use `SeedTree` consistently. |
| Independent reference and fast engines | graph-zigzag reference vs dynamic-connectivity/MSF paths; naive/reflection/coned zigzag checks | Optimized structures remain accountable to a mechanistically different oracle. | This is one of the strongest patterns in the tree; expand bounded differential/property coverage before trusting asymptotic claims. |
| Alternate dictionary lookup and frozen fitted maps | IDF/TF-IDF span lookup and `FrozenDictionary` models | Avoids a token allocation per lookup and separates fitting from repeated query use. | Complete the immutability boundary: arrays and option objects exposed by an “immutable” model must not remain mutable. |

## Cross-cutting implementation guidance

### 1. Make algorithm identity and input basis explicit

The hashing and text-feature files currently mix UTF-16 code-unit loops, UTF-8 encoding, low-byte
character truncation, native-endian memory reinterpretation, normalization, case folding, and
seeded variants. Those choices are not interchangeable. Any digest or model that leaves process
memory should identify at least:

- algorithm and project-variant version;
- input basis and encoding (`char` code units versus bytes, plus byte order where relevant);
- normalization, case, tokenization, and shingle policy;
- seed/domain parameters and signature dimensions; and
- merge/error semantics for sketches.

Names such as TLSH, CTPH/ssdeep, MinHash, and HyperLogLog imply recognizable algorithms. Either
match an external specification and validate with published/interoperability vectors, or name the
implementation as a ThermoMapper-specific variant. “Inspired by” is a legitimate result; accidental
wire-compatibility claims are not.

### 2. Separate exact oracles, lossy signatures, candidate indexes, and streaming estimates

`hashish` contains all four, but they need different tests:

- exact Levenshtein/Jaccard/cosine results need direct mathematical or brute-force oracles;
- SimHash/MinHash/TLSH/CTPH need calibration, invariance, and distributional tests;
- LSH and Bloom filters need false-positive/false-negative contracts at the candidate boundary;
- Count-Min and HyperLogLog need error, merge, saturation, and large-cardinality tests.

A candidate collision is not equality. An approximate estimate should carry enough configuration
to reproduce its error model. Exact methods should remain available as calibration or verification
paths.

### 3. Use pooling at an operation lifetime, not inside the innermost repeated unit

The repository has good pool ownership in graph and geometry operations. The weakest use is where
scratch is rented or text is re-encoded once per token × hash-function pair. Prefer:

1. preprocess or encode once;
2. rent one workspace for the whole operation or worker;
3. perform all inner iterations over spans in that workspace; and
4. return the workspace in one `finally`.

Pooling is not automatically faster for small arrays. Record allocated bytes and throughput by
input distribution, especially around stack/pool thresholds.

### 4. Finish “immutable model” boundaries

`FrozenDictionary` is useful, but `TfIdfModel`, `TokenizedCorpus`, and `CooccurrenceModel` expose
mutable arrays or mutable option objects. A caller can therefore invalidate vocabulary-to-column
identity after fitting. Use private arrays with read-only views/copies, immutable option snapshots,
and explicit serialization artifacts when model identity matters.

### 5. Define deterministic ties and parallel thresholds

Stable results require more than a deterministic seed. Frequency pruning, top-K results, equal
scores, and sorted candidate outputs need explicit secondary keys. `Parallel.For` should not be a
default for every corpus size: add a workload threshold or a caller-owned execution policy, and
support cancellation where transforms may be large.

## `hashish` file review

### Hashes, sketches, and candidate indexes

| File | Useful capability | ThermoMapper follow-up |
|---|---|---|
| [`seeded.cs`](../src/hashish/seeded.cs) | Allocation-free seeded span hashing and a reusable 64-bit mixer. | Freeze the distinct `char`, byte, and `uint` input encodings; add domain/version identifiers before persistence. The `char` overload hashes code-unit values, not the byte sequence used by the byte overload. |
| [`bloom.cs`](../src/hashish/bloom.cs) | Word-packed approximate membership, double hashing, `PopCount` fill diagnostics. | Add merge/serialization/config identity if this becomes an artifact; test requested false-positive rates, insertion-count overflow, and concurrent-use expectations. |
| [`countmin.cs`](../src/hashish/countmin.cs) | Streaming approximate frequency with explicit epsilon/delta sizing. | Replace `long[,]` with a flat row-major buffer for the hot loop; define overflow, merge, conservative-update, and concurrent accumulation semantics. Validate empirical error against an exact counter. |
| [`hyperloglog.cs`](../src/hashish/hyperloglog.cs) | Compact distinct-count estimate with merge and leading-zero rank. | Add standard-vector and range/bias tests, artifact identity, and serialization. “Register” here is the standard HLL bucket term and should remain local to that algorithm, not become a project-wide namespace. |
| [`minhash.cs`](../src/hashish/minhash.cs) | MinHash signatures plus a banded LSH candidate index. | The current path allocates one string per unique character shingle and re-encodes every shingle once per hash slot. Hash each shingle once into a stable base value, then apply a validated hash family/permutation. Define the too-short/empty-signature result: two all-`uint.MaxValue` signatures currently estimate similarity as 1.0. Stamp shingle/band parameters, dedupe repeated ids, and retain exact Jaccard verification. |
| [`simhash.cs`](../src/hashish/simhash.cs) | Weighted 64-bit signatures, stack accumulator, and hardware `PopCount` comparison. | `Regex.Matches`, `Match.Value`, and `ToLowerInvariant` allocate per document/token. Reuse the shared span enumerator/tokenization artifact. Treat BM25/IDF/tokenization parameters as signature identity and test tie-at-zero bit behavior. |
| [`ctph.cs`](../src/hashish/ctph.cs) | Content-triggered chunking concept, dual resolutions, and span/two-row edit comparison. | **Audit compatibility.** The trigger is a cumulative FNV prefix, not an evicting rolling window, and digest bytes come from native-endian `ulong` reinterpretation. The block-size rule, truncation, cross-resolution comparison, malformed/zero-block parsing, and published ssdeep behavior need oracle vectors. Conform or rename as a project variant. |
| [`tlsh.cs`](../src/hashish/tlsh.cs) | Sliding-window bucket histogram, quartile coding, compact bit packing. | **Audit compatibility.** Window buckets use only each UTF-16 code unit's low byte while the checksum uses UTF-8, and comparison is a simplified character mismatch count rather than the standard TLSH distance. The long-input checksum allocates a full string and byte array. Conform to TLSH vectors or rename the digest. |

### Exact and heuristic similarity measures

| File | Useful capability | ThermoMapper follow-up |
|---|---|---|
| [`levenshtein.cs`](../src/hashish/levenshtein.cs) | Strong implementation shape: common-affix trimming, shorter column, two rows, stack/pool threshold, row swapping. | Add a banded/max-distance form for screening and a separately named trace/edit-script form when correspondence is needed. Test lone-surrogate/code-unit semantics and threshold edges. |
| [`jaccard.cs`](../src/hashish/jaccard.cs) | Exact Jaccard, containment, overlap, and Dice; useful calibration oracle for MinHash. | Every call materializes/clones both sets, even when the inputs are already `HashSet<T>`. Add prebuilt-set or span/sorted-input kernels and keep the empty-set conventions explicit. |
| [`cos.cs`](../src/hashish/cos.cs) | Span/TensorPrimitives cosine kernels, in-place normalization, upper-triangle matrix fill. | **Correctness fix:** the comment says a zero-norm input has distance 1.0, but `Similarity` maps NaN to 0 and `acos(0)/π` returns 0.5. Decide and test the intended convention. Validate every row length; longer rows are currently truncated to the first row's dimension. Guard `n*d` overflow. |
| [`ncd.cs`](../src/hashish/ncd.cs) | Codec-parametric normalized compression comparison. | Clarify whether the measure is directed or symmetrized (`C(xy)` can differ from `C(yx)`), bound or document finite-sample values outside `[0,1]`, and avoid new byte arrays/streams for every pair in large searches. Validate codec headers and empty inputs against an oracle. |
| [`measure.cs`](../src/hashish/measure.cs) | A small common distance/similarity vocabulary with struct adapters. | No in-tree consumer currently exploits this abstraction. Either integrate a real generic/static-dispatch consumer and benchmark it, or keep direct functions and remove the unused layer; do not claim devirtualization from adapter shape alone. |

### Tokenization, feature models, and search

| File | Useful capability | ThermoMapper follow-up |
|---|---|---|
| [`tokenizer.cs`](../src/hashish/tokenizer.cs) | One explicit normalization/tokenization stage and a reusable compiled regex. | Add a span-enumeration/token-view path so downstream features do not require `MatchCollection` plus token strings. Document that disabling compatibility normalization still performs canonical normalization. |
| [`shingler.cs`](../src/hashish/shingler.cs) | Ordered and set-valued word n-grams. | `string.Join` materializes every shingle, then `BuildSet` materializes the array and another set. Add hashed/index-window shingles or a streaming enumerator for MinHash/containment workloads. |
| [`histogram.cs`](../src/hashish/histogram.cs) | Zero-allocation normalization into a caller span; pooled unigram counts; explicit smoothing. | Reject negative counts and negative/non-finite alpha, define integer overflow behavior, and consider a reusable workspace for repeated same-vocabulary queries. |
| [`idf.cs`](../src/hashish/idf.cs) | Excellent alternate `ReadOnlySpan<char>` dictionary lookup, within-document deduplication without a per-document set, frozen fitted maps. | Validate formula enums and model parameters, snapshot text-preprocessing identity, and make persistence/version identity explicit. Consider a mergeable document-frequency accumulator for streaming or parallel fitting. |
| [`bm25.cs`](../src/hashish/bm25.cs) | Convenient bridge from the richer IDF model to SimHash. | It is explicitly a legacy shim. Under the repository's clean-break policy, either make it a real named BM25 statistics artifact or migrate callers and delete the compatibility surface. |
| [`cooc.cs`](../src/hashish/cooc.cs) | Two-pass vocabulary build, flat row-major matrix, frozen token map, reusable tokenized corpus. | **Audit counting:** the center loop visits both directed positions and each visit increments both symmetric cells, apparently doubling every pair and marginal. The global factor cancels in some probabilities but not in exposed raw counts, overflow, or work. Process each unordered window pair once if symmetric counts are intended. Guard `vocabSize*vocabSize`, use a stable token tie-break for equal frequencies, consider sparse/blocked storage, and stop exposing mutable arrays from the “immutable” model. |
| [`cooc_stats.cs`](../src/hashish/cooc_stats.cs) | PMI/PPMI, conditional probability, contextual entropy, and neighbor inspection. | Reuse the bounded heap already implemented in `TfIdfSearch`; sorting a `topN+1` list on every candidate is avoidable. Add index validation, sparse-row forms, smoothing choices, and a flat/streaming PPMI output to avoid jagged `V²` allocation. |
| [`tfidf.cs`](../src/hashish/tfidf.cs) | Reusable tokenized corpus, deterministic vocabulary columns, dense/sparse transforms, alternate lookup, TensorPrimitives, parallel row independence. | Validate all options at construction; make fitted options/arrays truly immutable; add a sparse-corpus artifact rather than only per-document sparse output; choose parallel thresholds. Handle `Dimension == 0` throughout consumers. |
| [`tfidf_search.cs`](../src/hashish/tfidf_search.cs) | Sparse-query × dense-row scoring and a good bounded parallel-array min-heap (`O(N log K)`). | Refuse or correctly normalize when `L2Normalize == false`; currently dot product is called cosine regardless. Handle zero-dimensional models before modulo/division, define equal-score document-id ties, and avoid parallel overhead below a measured threshold. |

## Broader source review

### Graph layout and traversal

- [`CsrGraph`](../src/graphs/primitives/CsrGraph.cs) contains a sound reusable construction shape,
  but its public mutable arrays make “graph identity” fragile. `FromEdges` assumes valid endpoints,
  no duplicates, and whatever row ordering the input happens to induce. If consumers require sorted
  rows, validate or sort them once during construction. `BuildReverseSlotMap` linearly searches the
  opposite row per directed edge and can degrade toward the sum of squared degrees; build reverse
  slots during a paired fill or with a one-pass edge-key map. The binary form needs version,
  endianness, stronger pointer monotonicity/range validation, and a declared duplicate policy.
- [`UndirectedEdgeWalk`](../src/graphs/primitives/UndirectedEdgeWalk.cs) is a useful zero-allocation
  pattern enumerator and consolidates an error-prone idiom. Retain its own adoption warning for
  fused ThermoMapper loops: benchmark lowered code and branch behavior before replacing them. Add
  differential enumeration tests over empty, self-loop, asymmetric, duplicate, and unsorted input,
  even if invalid cases are rejected by the graph constructor.
- [`PathNeighborRefiner`](../src/graphs/pipeline/refinement/PathNeighborRefiner.cs) and
  [`Dijkstra`](../src/graphs/primitives/traversal/Dijkstra.cs) demonstrate good per-worker scratch,
  reusable queues, bounded search, and target early exit. The current reset and target count still
  scan/fill `O(N)` arrays for each source. For sparse rows, stamped visitation arrays and a touched
  target list are concrete benchmark candidates. Make maximum parallelism and cancellation
  caller-controlled.
- [`UnionFind`](../src/graphs/primitives/UnionFind.cs) has the right reset/output-span shape. Consider
  exposing a caller-owned buffer constructor only if repeated construction remains measured; keep
  path compression and union-by-size as the reference semantics.

### Numerical kernels and scratch ownership

- [`EarthMover.Distance1D`](../src/maths/distance/EarthMover.cs) correctly targets repeated-call GC
  pressure, but **the mathematical input contract needs an oracle**. If arrays are samples, 1-D
  Wasserstein-1 is the mean absolute difference of sorted samples; if arrays are histogram masses,
  cumulative differences are appropriate but sorting destroys bin order. The current code sorts
  and then compares cumulative sums. It also silently truncates a longer `b` while a shorter one
  throws. Define the domain, require equal lengths, and test against a reference implementation.
- [`GraphLaplacian`](../src/graphs/spectral/GraphLaplacian.cs) is a good example of layout-aware
  output and scoped scratch. Keep the flat column-major entry point as the primitive consumed by the
  solver; ensure its weight conventions agree with the rectangular builder.
- [`MatrixOps`](../src/maths/linalg/MatrixOps.cs),
  [`CoherentField`](../src/graphs/spectral/CoherentField.cs), and
  [`EigenFast`](../src/maths/linalg/EigenFast.cs) contain valuable hardware-tiering patterns. They
  also create a large unsafe/manual surface. Benchmark them against `TensorPrimitives`, verify
  scalar-tail and small-length cases on every supported architecture, and document reduction-order
  differences where bitwise reproducibility is not promised.
- [`ScatterAccumulator`](../src/maths/geometry/estimators/intrinsic/ScatterAccumulator.cs) combines
  a generic struct kernel, flat destination, and one operation-lifetime rental well. If symmetry is
  guaranteed, half-matrix accumulation plus mirroring is a benchmark candidate. Verify that the
  generic constraint actually devirtualizes on supported runtimes before treating it as a reason
  to reshape public APIs.
- [`GaussianMixtureModel`](../src/clustering/statistical/gmm/GaussianMixtureModel.cs) correctly uses
  log-sum-exp. Its E-step dynamically `stackalloc`s `K` doubles but only comments that a pool should
  be used for very large `K`; implement the threshold rather than relying on callers to stay small.
- [`OnlineMahalanobis`](../src/maths/linalg/WelfordMahal.cs) is a useful streaming-statistics
  primitive. `_scratch` is unused and should be removed. Add mergeable `(count, mean, M2)` state for
  parallel streams if needed, and define counter overflow/reset behavior.

### Persistence and dynamic structures

- The reference/fast separation around graph zigzag is architecturally strong. The fast paths keep
  integer-id dynamic structures below the domain layer and name independent ground truth. Retain
  that arrangement while benchmarking the object-heavy dictionary/hash-set/treap implementation;
  asymptotic comments are not a throughput measurement.
- [`PersistenceClearing.ComputeH0`](../src/tda/ph/PersistenceClearing.cs) is allocation-heavy and
  duplicates union-find logic, but correctness comes first: **audit the elder rule**. The branch
  named `elder` appears to choose the component with the later/larger birth (and larger tie id),
  then kills the earlier component. A two-vertex, one-edge filtration should settle whether the
  names or the pairing are reversed. Once correct, replace edge-list/lambda/dictionary hot paths
  only with oracle coverage in place.
- [`LazyRipsFiltration`](../src/tda/ph/LazyRipsFiltration.cs) avoids constructing a
  `SimplicialFiltration` object but still materializes all discovered simplices, an index, and all
  cofacet arrays; “lazy” should be read narrowly or renamed. More importantly, its merge-style
  `CommonNeighbors` requires sorted CSR rows, while `CsrGraph.FromEdges` does not establish sorted
  targets. **Audit with unsorted edge input** because triangles may be missed.
- [`FastZigzag`](../src/tda/ph/FastZigzag.cs) shows a valuable transform-to-a-known-engine pattern,
  but it creates many lists, dictionaries, arrays, and coned cells. Pre-size from the replay counts
  and compact incarnation/event data only after the existing independent oracles remain green.

## Suggested ThermoMapper sequence

This is a maintenance priority list, not a Doccer dependency graph.

### P0 — correctness and identity

1. Add minimal oracle tests for `EarthMover.Distance1D`, the persistence elder rule, and unsorted
   CSR rows in `LazyRipsFiltration`.
2. Fix or explicitly document cosine zero-vector distance and TF-IDF non-normalized scoring.
3. Decide whether CTPH and TLSH are interoperable implementations or ThermoMapper-specific
   variants; add published vectors or rename them.
4. Freeze digest/model input bases and parameters; handle zero-dimensional and empty-signature
   cases explicitly.
5. Audit co-occurrence double counting, integer overflow, deterministic vocabulary ties, and model
   mutability.

### P1 — high-payoff allocation and layout

1. Remove MinHash's string-per-shingle and encode-per-seed multiplication.
2. Add streaming/token-view paths shared by tokenizer, shingler, SimHash, IDF, TF-IDF, and
   co-occurrence.
3. Flatten Count-Min storage; provide sparse/blocked co-occurrence storage; reuse the bounded top-K
   heap for context neighbors.
4. Implement actual stack/pool thresholds in GMM and operation/worker workspaces where repeated
   rentals remain hot.
5. Benchmark stamped/touched Dijkstra reset, reverse-slot construction, and manual SIMD paths.

### P2 — artifact and execution maturity

1. Make fitted models deeply immutable and serializable with algorithm/preprocessing identity.
2. Add mergeable streaming states for IDF/counting/HLL/Count-Min/Welford where the math permits it.
3. Add deterministic score ties, parallel thresholds, maximum degree, and cancellation.
4. Establish benchmark families that record elapsed time and allocated bytes for dense, sparse,
   short, long, ASCII, BMP, surrogate, empty, and adversarial inputs.

## Verification matrix for future changes

| Family | Correctness evidence | Performance evidence |
|---|---|---|
| Exact distances/statistics | brute-force or published oracle; symmetry/bounds/empty cases; Unicode basis | short/long, near-equal/dissimilar, stack/pool thresholds, allocated bytes |
| Digests/signatures | published vectors where the name implies a standard; deterministic serialization; parameter identity | throughput by input length and entropy; preprocessing allocation; comparison cost |
| Sketches | exact-counter calibration; merge equivalence; empirical error/confidence; saturation | update/query throughput, cache behavior, flat vs multidimensional storage |
| Feature models | independent small-corpus hand calculations; stable vocabulary/ties; immutable round-trip | fit versus transform, dense versus sparse, parallel crossover, memory peak |
| Graph/dynamic kernels | mechanistically independent reference engine and bounded randomized differential tests | graph size/degree/update mix, allocation, branch/cache profiles, crossover point |
| SIMD/manual intrinsics | scalar differential tests including tails, NaN/Inf/zero and unsupported hardware | runtime tier matrix; JIT/TensorPrimitives comparison; cold and steady-state runs |

The useful conclusion is not that `hashish` is good or bad as a block. It is a broad prototype
surface containing several strong low-allocation techniques, several valuable numerical/indexing
capabilities, and a handful of correctness and identity questions that need focused oracles before
the module becomes a trusted shared substrate.
