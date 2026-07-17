# ISOLET Unsupervised Benchmark Brief

**Status:** proposed research and implementation plan  
**Date:** 2026-07-15  
**Primary comparison:** Blatt, Wiseman, and Domany's full-data ISOLET SPC result  
**Constraint:** preprocessing, model selection, and clustering are fully unsupervised

## Executive decision

ISOLET should proceed as two deliberately separate tracks:

1. **Forensic SPC parity:** reproduce the published result on the supplied 617-dimensional
   representation without preprocessing.
2. **Unsupervised improvement tournament:** compare carefully controlled preprocessing and graph
   variants, including SPRED, after the parity surface is credible.

The parity track is the control experiment. It must not absorb a modern scaling rule, projection,
resolver, or graph heuristic merely because that change improves the score. The improvement track
may depart from the 1997 method, but each departure must be named, serialized, and evaluated against
the frozen parity configuration.

The working hypothesis is not simply "PCA failed, therefore SPRED wins." The current evidence says
that the raw representation already contains substantial local letter structure and that graph and
result-resolution fidelity are likely the first bottleneck. The strongest representation hypothesis
is:

```text
stable unsupervised feature selection
    -> correspondence-conscious H0 SPRED
    -> parity-grade SPC
```

Plain SPRED remains an important arm, not an assumed winner.

## Questions this benchmark must answer

1. Can ThermoMapper reproduce Domany's raw-617 ISOLET graph, thermodynamic response, selective
   clustering result, and hierarchy closely enough to claim SPC parity?
2. Does an unsupervised front end improve SPC at a fixed purity or fixed coverage without selecting
   parameters from the letter labels?
3. Does SPRED preserve structure that is more useful to SPC than PCA or a dimension-matched random
   projection?
4. Is any SPRED gain due to persistent-homology preservation, generic dimensional denoising, changed
   graph density, or stochastic selection?
5. Can the selected representation generalize across speaker cohorts rather than organizing mainly
   by speaker, gender, recording cohort, or paired utterance?

## What ISOLET is, and is not

The UCI matrix contains 7,797 observations, 617 continuous features, and a 26-valued letter label.
The features were produced by a speech-specific pipeline and include spectral coefficients, contour
measurements, and measurements around sonorant, pre-sonorant, and post-sonorant segments. The
features are already scaled to `[-1, 1]`.

This is not a raw waveform matrix. The exact ordering and grouping of the 617 feature columns is not
known. Consequently:

- A DWT across the 617 columns has no justified time or frequency adjacency.
- Generic image-like or sequence-like convolutions over feature index are not justified.
- The Wave_clus analogy belongs at the **feature-selection** stage: identify measurements with stable
  multimodal or locally discriminating structure, then cluster in that selected representation.
- A new waveform pipeline using wavelets, scattering transforms, MFCCs, or modulation spectra would
  be a separate benchmark requiring the original audio. It would not be Domany parity.

## Correction of the local record

The previous "ISOLET PCA wall" interpretation is not an adequate benchmark specification.

The published SPC result is:

- all 7,797 observations;
- all 617 supplied features;
- no PCA front end;
- 93% purity;
- 35% unclassified, therefore approximately 65% coverage;
- a hierarchy across temperature, not a single forced 26-way partition;
- labels used for external interpretation and scoring, not to fit the clustering.

Therefore the primary target is **selective hierarchical clustering**, not flat accuracy and not
"number of letters captured" alone. Purity without coverage is gameable by emitting many tiny
clusters, and coverage without purity hides mixed letter clusters.

Existing comments and tests that call ISOLET an established unsupervised ceiling should be treated as
historical controls until the parity track below is complete. HDBSCAN failing on one PCA
representation does not establish that the dataset lacks unsupervised structure.

## Evidence already in hand

The following exploratory measurements were made on the committed `datasets/isolet.csv.gz`. Labels
were used only to understand the failure, so these values are diagnostic evidence and must not become
selection criteria in the sealed benchmark.

| Diagnostic | Observation |
|---|---:|
| Data matrix | 7,797 x 617, plus label |
| Constant features | 0 |
| Raw Euclidean 1-NN with same letter | 90.8% |
| Mean same-letter fraction among first 5 neighbors | 81.8% |
| Mean same-letter fraction among first 10 neighbors | 80.2% |
| PCA variance, 20 dimensions | 71.6% |
| PCA variance, 30 dimensions | 77.2% |
| PCA variance, 50 dimensions | 83.5% |
| PCA dimensions for 95% variance | 173 |
| Covariance effective rank | about 45.4 |
| Covariance participation ratio | about 12.1 |

The existing local SPC results were far weaker than the local neighbor evidence:

| Representation | Local result |
|---|---|
| raw617 | 12 letters, 0.752 purity, 0.13 coverage |
| PCA50 | 14 letters, 0.748 purity, 0.19 coverage |
| PCA30 | 14 letters, 0.834 purity, 0.15 coverage |

This discrepancy is the reason graph construction, coupling, temperature sampling, and result
resolution precede representation work.

The graph edge audit provides another useful fingerprint:

| Local graph recipe | Undirected edges |
|---|---:|
| mutual `k=5` plus MST repair | about 11,547 |
| mutual `k=10` plus MST repair | about 20,249 |
| mutual `k=11` plus MST repair | about 22,168 |
| mutual `k=12` plus MST repair | about 24,112 |
| published ISOLET interaction graph | 22,471 |

The `k=11` result is a forensic clue, not proof that the authors used `k=11` or MST repair.

## Track A: forensic Domany parity

### A1. Frozen input

- Use the supplied 617 features exactly as scaled.
- Do not z-score, whiten, center by class, project, filter, or impute.
- Preserve observation order and a stable row identifier.
- Hash the decompressed feature matrix and label vector.
- Split labels from features before any graph or SPC API is called.

Translation does not affect Euclidean distances, but even harmless transformations should remain out
of the parity manifest so that the provenance is unambiguous.

### A2. Graph reconstruction

The paper specifies mutual K-nearest-neighbor adjacency in high dimensions and reports 22,471 ISOLET
edges. It says K is normally chosen to connect the interaction graph. It also describes a
`K=10`-plus-MST variant for very high dimensions, while saying that variant was used only in two
earlier examples. The ISOLET section does not state K explicitly.

The parity investigation must therefore record, for each plausible interpretation:

- directed KNN count before reciprocity;
- mutual edge count;
- connected-component count and sizes;
- isolated vertices;
- MST or connectivity-repair edge count;
- final edge count;
- degree minimum, median, mean, maximum, and tail;
- edge-distance quantiles;
- deterministic tie handling;
- exact-versus-approximate neighbor search;
- floating-point and duplicate-distance behavior.

Test a narrow forensic bracket around the observed edge fingerprint. Do not search K using labels.
The winning parity interpretation is the one best supported jointly by the paper, connectivity,
22,471-edge fingerprint, thermodynamic response, and downstream hierarchy.

### A3. Coupling reconstruction

For neighboring observations `i` and `j`, reconstruct the published interaction:

```text
J_ij = (1 / K_bar) * exp(-d_ij^2 / (2 * a^2))
```

where `a` is the mean distance over neighboring pairs and `K_bar` is the average graph degree. Record
`a`, `K_bar`, and coupling quantiles in the run artifact.

The parity arm uses this global scale. Self-tuning or locally scaled couplings belong only to Track B.

### A4. Potts and temperature protocol

- `q = 20` Potts states.
- `M = 1000` Swendsen-Wang sweeps per temperature.
- Use a temperature scan capable of resolving all susceptibility transitions; the paper describes
  typical runs of about 20 temperatures.
- Record burn-in policy, estimator policy, RNG algorithm, and seed.
- Refine temperature intervals around susceptibility peaks without looking at labels.
- Emit susceptibility, magnetization, cluster-size trajectories, and edge-correlation summaries.
- Repeat enough seeds to distinguish stochastic spread from a systematic parity failure.

### A5. Published resolver

At each selected super-paramagnetic regime:

1. Create core links where neighboring spin correlation `G_ij > 0.5`.
2. Link every point to its neighbor of maximum correlation to capture cluster peripheries.
3. Take connected components of the resulting graph.

Any minimum cluster size, lineage persistence, unclassified rule, or cross-temperature consolidation
must be explicit. A modern resolver can be evaluated later, but it must not silently replace this
published extraction during parity.

### A6. Parity outputs

Parity is not one scalar assertion. Preserve:

- graph manifest and edge fingerprint;
- coupling summary;
- full temperature schedule;
- susceptibility and magnetization traces;
- per-temperature component sizes;
- correlation distributions;
- resolved clusters and unclassified observations;
- cross-temperature lineage or hierarchy;
- purity-coverage operating points;
- the letter-group hierarchy used only after labels are revealed.

The headline target remains approximately 93% purity at 65% coverage. A miss should be decomposed into
graph, simulation, resolver, and hierarchy differences rather than attributed immediately to the data.

## Track B: unsupervised improvement tournament

Track B begins only after Track A has a frozen, credible baseline. Each candidate changes one
scientific factor at a time during screening. Full Cartesian products are prohibited until the
funnel has selected a small number of plausible finalists.

### B1. Scaling controls

Evaluate:

1. Raw supplied scaling.
2. Robust diagonal scaling using median and MAD, with a declared floor for very small MAD.
3. Conventional z-scoring as a control.

Do not begin with PCA whitening or full covariance whitening. Both can amplify low-variance noise and
erase the feature weighting introduced by the original speech pipeline. The input is already bounded,
so winsorization is not a default; it requires an observed tail pathology.

### B2. Wave_clus-inspired stable modality selection

This is the preprocessing arm most closely aligned with the spike-sorting workflow.

For every feature:

1. Estimate its one-dimensional density over a declared bandwidth or smoothing family.
2. Measure evidence for multiple separated modes using a cluster-tree, excess-mass, dip, or equivalent
   mode-persistence score.
3. Repeat over bootstrap samples and source speaker cohorts.
4. Reward modes that persist across smoothing scales and samples.
5. Penalize unstable modes and tiny tail bumps.
6. Rank features without labels.
7. Prune redundant selected features by absolute correlation or another declared dependence measure.

The historical Wave_clus Lilliefors/non-Gaussianity ranking should be retained as a control, not the
primary selector. Non-Gaussianity can represent skew, kurtosis, outliers, gender, or speaker effects
rather than letter modes.

Candidate feature budgets are 25, 50, and 100. The final budget is selected from an unlabeled stability
plateau, not from purity, NMI, or the known number 26.

### B3. Geometry-preserving feature-selection control

Use an unsupervised Laplacian-score-style ranking over the frozen raw parity graph:

- no letter labels;
- no requested number of clusters;
- graph recipe fixed before feature scoring;
- redundancy pruning after ranking;
- the same 25, 50, and 100 feature budgets.

This asks whether preserving the local raw graph is more useful than selecting marginal multimodality.
Because the raw graph may contain speaker nuisance, it is a control rather than ground truth.

### B4. Linear projection controls

Evaluate the following at target dimensions 20, 30, and 50:

- non-whitened PCA;
- seeded Gaussian or orthogonal random projection;
- SPRED on raw features;
- SPRED on the winning unsupervised feature subset.

Random projection is mandatory. It distinguishes a SPRED or PCA-specific gain from generic
distance concentration relief and reduced dimension.

Do not start with PCA followed by SPRED. PCA can remove topology before SPRED sees it. If a compute
preconditioner is later necessary, a high-dimensional PCA truncation must be labeled as a separate arm
and compared with a matched random preconditioner.

### B5. Graph and coupling improvements

Only after the representation screen should the leading candidates be crossed with modern graph
variants:

- self-tuning Gaussian coupling with `sigma_i = d(i, k_scale)` and
  `exp(-d_ij^2 / (2 sigma_i sigma_j))`;
- shared-neighbor edge similarity;
- mutual-proximity or another declared hubness correction;
- graph recipes that vary K in a narrow, unlabeled stability interval.

Retain the `1 / K_bar` thermodynamic normalization unless its removal is itself the isolated
experimental factor.

For every graph, report reciprocal-neighbor rate, in-degree skew, hubness, connectivity, repair edges,
and perturbation stability. Do not tune graph parameters against letter purity.

## SPRED experiment design

### Why H0 first

H0 records connected components merging over scale and is more directly related to a cluster tree than
H1 loops. ISOLET's clustering target gives no prior reason to spend half of the objective on loops.
The equal H0/H1 default is a general-purpose choice, not an ISOLET result.

The initial SPRED configuration should therefore use:

```text
MaxDimension      = 1
Dimensions        = [(0, 1.0)]
WassersteinOrder  = 2
Filtration        = RawDistance
MinPersistence    = 0
VarianceRegularizer = 0
target dimensions = {20, 30, 50}
```

This also avoids triangle construction and the H1 Hungarian matching cost.

### H0 matching-cost gate

Avoiding H1 does not make evaluation cheap. With `MinPersistence = 0` every finite H0 bar survives,
so each objective evaluation performs an exact Hungarian match on roughly one bar per observation
per side, and H0 bars do not prune the way the near-diagonal H1 noise loops did. The committed P0
profile measured W(H0) at about 165 ms for n = 200 with roughly cubic growth. At the planned
8-block split (about 975 observations per block) that extrapolates to 15-20 seconds per objective
evaluation per block, multiplied by annealing iterations, seeds, and target dimensions. Exact
full-data H0 matching at 7,797 observations is not feasible at all.

Phase 2 therefore does not begin until at least one of the following is landed and recorded:

1. an approximate diagram distance (entropic Sinkhorn or sliced Wasserstein, validated against the
   `T4transport` oracle) adopted as the screening metric;
2. a declared subsampling or landmark reduction with the cap serialized in the run artifact;
3. a pilot-derived wall-clock budget showing exact block-level matching fits the available compute
   at the chosen iteration count.

The two pilot seeds at dimension 30 must produce the full-screen wall-clock extrapolation before
any promotion decision.

**Status 2026-07-16 — condition 1 landed, with residuals.** `DiagramMetrics.SlicedWasserstein`
(deterministic diagonal-augmented slices, O(L·n log n)) and `DiagramMetrics.SinkhornWasserstein`
(log-domain entropic on the exact cost geometry) are in, selectable per run via
`PersistenceObjectiveConfig.DiagramDistance`. Validation: Sinkhorn converges to the in-repo exact
Hungarian oracle (ε→0, small diagrams, ≤1% relative), sliced matches its analytic single-bar
integral (2√2/π), and all three backends preserve the clean-vs-collapsed projection ranking at
objective level. Measured at n = 200 bars/side: exact 32 ms, sliced 4 ms, Sinkhorn 1.1 s
(iteration-capped, entropically biased at that scale). **Sliced is the screening metric** — at the
8-block split its n·log n replaces the cubic 15–20 s/eval extrapolation with tens of milliseconds;
Sinkhorn is a small-diagram fidelity tool, not a block-scale screen.

**Status 2026-07-16 (later) — condition 1 CLOSED.** The `T4transport` oracle clause is discharged
by `r/oracles/transport_oracle.R` + `tests/oracle/DiagramMetricsTransportParityTests.cs` (live
R toolchain, skip-when-absent), on an s = 100 diagonal-augmented fixture reconstructed
independently in R from raw bars, at p = 1 and p = 2:

- sliced vs `T4transport::swdist`: ≤ 8% tolerance against a measured 3–4% residual that is
  T4transport's own interpolated-quantile smoothing (its per-slice 1-D transport linearly
  interpolates ecdf quantiles on a 1000-point grid) plus ~1% Monte Carlo error; its scalar
  `distance` is also mean(W_p), not the documented (mean W_p^p)^(1/p), so the oracle recombines
  from `projdist`;
- Sinkhorn vs `T4transport::sinkhornD` at matched smoothing (λ = ε·cMax, identical kernel): ≤ 1%;
- exact Hungarian vs `lpSolve::lp.assign`: ≤ 1e-8 — an external oracle for the exact path too
  (the ε→0 Sinkhorn limit is pinned through this externally closed chain plus the unit suite).

**Status 2026-07-17 — pilot recorded; the gate section is fully discharged.** The two S0 pilot
seeds (211, 223) at dimension 30 ran on the shuffled raw617 (permutation seed 41, hashes in the
artifact), 8 blocks, H0-only, sliced screening metric, parity-grade mutual-K10+MST recipe on both
sides, I = 100 iterations. Full record: `pilot/spred-pilot-s0-dim30.json`; runner:
`tests/tda/dim-reduction/SpredIsoletPilotTests.cs` (`Category=Benchmark`).

- **Cost**: run(I) ≈ 46.7 s + 0.025 s·I (parallelism 8). The fixed share is dominated by the
  full-data objective's 7,797-row ambient reference; the marginal share is the annealing loop —
  the sliced metric prices a full 8-block iteration at ~25 ms wall. Full S0 screen
  (3 dims × 5 seeds): **12.3 / 14.7 / 17.8 min at I = 100 / 500 / 1000** — versus ≈ 30+ hours
  under the pre-gate exact-matching extrapolation. Cost is no longer a screen-design constraint.
- **Health**: every objective finite (no pathology penalties observed); block locals 65.8–66.9,
  aggregate-on-block 66.9–67.7, full-data 185.86; aggregate-to-aggregate Grassmann distance
  across seeds **0.024** — the distributed pipeline is seed-stable.
- **Finding for the screen design**: at I = 100 the anneal barely leaves the PCA warm start —
  five of eight blocks end bit-identical to it across both seeds, and the accepted improvements
  are ≤ 0.07% of the local objective. Block-to-aggregate Grassmann angles run 1.51–1.94, so
  975-row blocks disperse substantially in Gr(617, 30) while the medoid-seeded aggregate stays
  put. Before the full screen, either raise I well past 100 or tune the proposal/cooling scale —
  otherwise S0 will not differentiate from the B6 PCA arm by construction.

**Probe 2026-07-17 — iteration budget is NOT the fix; the annealer proposal scale is.** One seed
(211) at I = 1000 vs the same-run PCA warm-start baseline
(`pilot/spred-probe-s0-dim30-i1000.json`): **seven of eight blocks end bit-identical to the warm
start (Grassmann distance 0.0000 — zero accepted proposals in 1000 iterations); block 5 accepts
exactly one first-iteration step** (improvement 0.006%, Grassmann move 0.1000 = the annealer's
step length `temp·0.1` at initial temperature). Mechanism: `SubspaceAnnealer` proposes isotropic
random horizontal tangents at schedule-fixed length with geometric cooling `0.99^iter`; in
Gr(617, 30) — intrinsic dimension 30·587 = 17,610 — the improving-direction fraction from a PCA
start is vanishing, and cooling extinguishes uphill acceptance by iter ≈ 500. Raising I cannot
fix this; the proposal mechanism can (dimension-aware / acceptance-targeted step adaptation,
structured two-plane rotations mixing one retained with one discarded direction, or a
gradient-informed Riemannian search). **S0 is blocked on annealer mobility, not compute.**

Cost-model correction (faithfulness): the probe's same-run differencing gives a marginal of
0.092 s/iteration — the pilot's 0.025 s was small-I differencing noise (the ~47 s fixed share
varies ±4 s run-to-run, swamping a 2.5 s marginal at I = 100). Corrected model:
run(I) ≈ 43 s + 0.092 s·I → full S0 screen at I = 1000 ≈ **34 min**. Still no constraint.

### Barcode limitation

For a connected weighted graph, the finite H0 death values are closely related to spanning-tree merge
scales. Wasserstein distance between H0 diagrams compares that multiset of scales, but does not say
which observations or components merged. A projection can therefore obtain a good H0 objective while
rearranging neighborhoods in a way that harms letter clustering.

Consequences:

- SPRED objective value is a diagnostic, not the model-selection verdict.
- Raw-neighborhood retention and SPC-lineage stability are required companion diagnostics.
- Plain SPRED must be compared with PCA and random projection at the same dimension.
- Failure of plain SPRED does not by itself refute a topology-aware projection; it may expose the
  missing correspondence term.

### Ambient and projected graph recipes

Use the parity-grade graph recipe for the ambient reference. Initially use the same construction rule
for each projected cloud, matching the paper-faithful SPRED idea of reconstructing topology after
projection.

Record both `Graph` and `ReferenceGraph` in the objective manifest. A frozen ambient topology used for
all proposals can be a performance ablation, but not the main scientific result because it changes
what SPRED is optimizing.

### Distributed execution

`DistributedSpred` currently partitions contiguous rows. ISOLET source order is not an exchangeable
random sample; its files and cohorts are organized around speakers. Contiguous block optimization can
therefore produce speaker- or cohort-specific subspaces.

Before distributed SPRED:

- deterministically shuffle feature rows without using labels;
- retain original row IDs and an invertible permutation;
- record the permutation seed and hash;
- begin with 8 blocks, roughly 975 observations per block;
- repeat with at least 5 shuffle/annealing seeds for finalists;
- test block-count sensitivity at 4, 8, and 16 only after an 8-block candidate is promising;
- report local objectives, aggregate-on-block objectives, aggregate objective, and principal angles
  among local and aggregate projections.

The existing corrupted-block diagnostic supports the robustness story, but it does not replace this
representativeness check.

### H1 escalation

H1 enters only if H0 SPRED is competitive and the ambient barcode contains stable loops above a
predeclared persistence threshold. The first mixed objective should be:

```text
Dimensions = [(0, 0.8), (1, 0.2)]
MaxDimension = 2
```

Use an unlabeled ambient-barcode noise elbow to choose `MinPersistence` and freeze it before downstream
label scoring. Keep equal `0.5 / 0.5` as a paper-general control, not the default.

### Correspondence-aware extension

If H0 barcode preservation is good but projected neighborhood or SPC stability is poor, the next SPRED
development candidate is a composite objective:

```text
F(P) =
    W2(D0(X), D0(PX))
    + gamma * LocalGraphDistortion(P)
```

`LocalGraphDistortion` should use only the frozen ambient graph and feature distances, for example a
weighted distortion of ambient edge lengths plus a small declared sample of non-edge constraints.
It must not consume labels or SPC cluster assignments.

The generic `SubspaceAnnealer` objective seam is already suitable for this experiment. The new
topology/geometry objective belongs above graph construction and receives graph information as data;
it must not introduce consumer semantics into `src/graphs`.

## Label-sealed protocol

The CSV physically contains labels, so software separation is mandatory.

1. The loader returns a feature artifact and a separately hashed label artifact.
2. Preprocessors, graph compilers, SPRED, SPC, and intrinsic selectors receive only features and row
   IDs.
3. The known class count `26` is not supplied to any algorithm or parameter rule.
4. Candidate selection uses unlabeled diagnostics only.
5. The candidate set and all configs are frozen before the evaluator opens the label artifact.
6. Labels are used once for the final comparative report, except for explicitly marked exploratory
   diagnostics already listed in this brief.

For a secondary generalization test, fit scaling, feature selection, and projection on official
ISOLET1-4 data, freeze the transform, and apply it to the held-out ISOLET5 speaker cohort. Preserve the
full 7,797-point transductive result separately because that is the Domany comparison. Do not infer
speaker identities from row order unless the mapping has been verified from source metadata.

## Unlabeled candidate selection

Use a Pareto screen rather than an arbitrary weighted score. Relevant diagnostics are:

- graph connectivity without excessive repair;
- reciprocal-neighbor fraction and hubness;
- graph overlap under row resampling and small feature perturbations;
- stability of susceptibility peaks;
- stability of cluster-size trajectories and lineages across SPC seeds;
- projection principal-angle stability across SPRED seeds and distributed permutations;
- local-versus-aggregate SPRED objective gaps;
- neighborhood trustworthiness/continuity relative to the declared ambient representation;
- runtime, pathology count, and convergence behavior.

Neighborhood retention is a guardrail, not truth: an improved representation is allowed to correct bad
raw neighbors. A candidate should nevertheless be rejected if its claimed topological preservation is
paired with unstable or wholesale neighborhood replacement.

## Final evaluation

### Primary selective metrics

Report a purity-coverage curve over every defensible operating point:

```text
purity   = sum over assigned clusters of majority-letter count / assigned count
coverage = assigned count / 7797
```

Primary readings:

- purity at approximately 65% coverage;
- coverage at approximately 93% purity;
- area under the purity-coverage curve over a common coverage interval;
- number and size distribution of assigned clusters;
- number and distribution of unclassified observations.

No result wins by reporting purity alone.

### Secondary flat metrics

At declared, label-independent cuts, report:

- adjusted Rand index;
- normalized mutual information;
- B-cubed precision, recall, and F1;
- per-letter recall and dominant confusions;
- treatment of unclassified observations stated explicitly.

These metrics are secondary because neither Domany SPC nor the intended ThermoMapper result is a
forced flat 26-cluster partition.

### Hierarchy metrics

Preserve and evaluate:

- temperature-indexed lineage tree;
- split order and persistence of large groups;
- dendrogram purity or another declared hierarchical agreement measure;
- qualitative recovery of major acoustic groups shown in the original paper;
- the label-derived class-centroid hierarchy as evaluator-only context.

Do not report the paper's 0.98 CPCC as an SPC score. That value describes the label-derived reference
hierarchy's fit under its own construction.

### Stochastic reporting

For finalist configurations:

- use at least 5 declared seeds;
- report every seed, median, range, and failures;
- do not publish only the best stochastic run;
- bootstrap final label metrics only after configs are frozen;
- distinguish variability from SPRED annealing, distributed partitioning, SPC simulation, and
  temperature-grid refinement.

## Provisional decision rules

### Parity is credible when

- the graph interpretation is paper-supported and the 22,471-edge discrepancy is either eliminated or
  explained;
- coupling and temperature diagnostics match the published scale and qualitative susceptibility
  structure;
- the published resolver is implemented independently of modern lineage heuristics;
- the result approaches the reported 93% purity and 65% coverage without label tuning;
- the recovered hierarchy has recognizable agreement with the published hierarchy;
- repeated seeds show that remaining differences are systematic rather than Monte Carlo noise.

### A preprocessing candidate is competitive when

- it was selected entirely from unlabeled diagnostics;
- it improves the median purity-coverage frontier over the frozen raw parity baseline;
- the improvement repeats across seeds and is not caused by tiny-cluster proliferation;
- it is at least as stable as the raw baseline;
- it is superior to a dimension-matched random projection if dimensionality reduction is claimed as
  the cause;
- its provenance and compute cost are reported.

A practical initial win threshold is either at least 2 percentage points more coverage at 93% purity,
or at least 1 percentage point more purity near 65% coverage, with no material hierarchy or stability
regression. The full curve and uncertainty remain authoritative.

### SPRED earns further development when

- it beats both PCA and random projection at a matched target dimension;
- the result is stable across distributed permutations and annealing seeds;
- objective improvement correlates with stable graph and SPC behavior;
- raw-to-SPRED or selected-features-to-SPRED gains justify the additional compute.

If SPRED only lowers its own objective, it has not yet won.

## Experiment funnel

### Phase 0: parity

| ID | Representation | Graph/coupling | Purpose |
|---|---|---|---|
| A0 | raw617 | paper reconstruction | primary parity |
| A1 | raw617 | narrow graph-interpretation bracket | resolve K/connectivity ambiguity |

### Phase 1: cheap representation screen

| ID | Representation | Dimensions/features |
|---|---|---|
| B0 | raw supplied scaling | 617 |
| B1 | robust diagonal scaling | 617 |
| B2 | z-score control | 617 |
| B3 | stable modality selection | 25, 50, 100 |
| B4 | Lilliefors control | 25, 50, 100 |
| B5 | Laplacian score | 25, 50, 100 |
| B6 | PCA, not whitened | 20, 30, 50 |
| B7 | matched random projection | 20, 30, 50 |

Use unlabeled graph and stability diagnostics to retain at most three non-raw finalists.

### Phase 2: SPRED screen

| ID | Input | Objective | Target dimensions |
|---|---|---|---|
| S0 | raw617 | H0 | 20, 30, 50 |
| S1 | winning selected features | H0 | feasible members of 20, 30, 50 |
| S2 | raw617 | H0 plus small variance reward | winning dimension only |
| S3 | winning input | H0/H1 = 0.8/0.2 | winning dimension only, conditional |

Begin with two pilot seeds at dimension 30 to expose pathology and cost. Promote only healthy
configurations to the full five-seed, three-dimension screen.

### Phase 3: graph improvement

Cross only raw parity and the best two representations with:

- original global coupling;
- one self-tuning coupling;
- one hubness/shared-neighbor candidate selected from unlabeled graph stability.

### Phase 4: finalists

Run full temperature resolution, published resolver, modern resolver as a separate output, five or
more seeds, final label evaluation, hierarchy comparison, and compute accounting.

## Required artifacts

Every run should emit immutable, serializable artifacts:

- dataset hash and row-ID manifest;
- preprocessing config, fitted parameters, selected feature indices, and transform hash;
- projection matrix and projection diagnostics;
- distributed permutation, block boundaries, block seeds, and aggregate diagnostics;
- graph compiler config and graph manifest;
- serialized `CsrGraph`;
- coupling config and summary;
- SPC config, RNG identity, seeds, temperature schedule, and sweep counts;
- per-temperature thermodynamic and cluster diagnostics;
- resolver config and result;
- evaluator config and final metrics;
- wall-clock and peak-memory measurements where available.

Aggregate reports should include:

- graph audit table;
- susceptibility and cluster-size plots;
- purity-coverage curves with seed uncertainty;
- hierarchy visualization;
- projection principal-angle matrix;
- SPRED local/aggregate objective plot;
- per-letter confusion and unclassified table;
- compute-versus-quality table.

## Layering and implementation constraints

- `src/graphs` remains a neutral constructor of weighted `CsrGraph` artifacts plus provenance.
- Local scaling, hubness correction, or topology-derived edge scores enter graph construction through
  typed data/config/delegates, never by importing SPC or PH consumers.
- SPRED remains under `tda/dim-reduction` and consumes graph and PH surfaces according to the existing
  layer order.
- Label-aware evaluation remains above all constructors and clustering engines.
- Configs are declarative and JSON-serializable. RNG instances, delegates, and fitted live objects do
  not enter persisted configs.
- The CLI/REPL may provide fluent ergonomics, but benchmark runners consume frozen DTOs.
- Superseded experimental surfaces are removed rather than preserved behind compatibility aliases.

## Risks and mitigations

| Risk | Mitigation |
|---|---|
| Edge-count matching overfits an implementation accident | Require joint paper, connectivity, thermodynamic, and hierarchy evidence |
| Purity is inflated by tiny clusters | Always pair purity with coverage and cluster-size distribution |
| Labels leak through known `k=26` or parameter search | Physically separate label artifact and freeze configs before evaluation |
| Feature multimodality captures speaker or gender | Require bootstrap and cross-cohort stability; inspect only after freeze |
| Z-scoring amplifies weak bounded features | Keep raw scaling primary and z-score as an isolated control |
| H0 SPRED preserves merge scales but not memberships | Require neighborhood diagnostics; add graph-distortion term only if indicated |
| H1 cost dominates with noise loops | Start H0-only; gate H1 on stable ambient persistence |
| Exact H0 Wasserstein matching is cubic in block size | Gate Phase 2 on an approximate distance, declared subsampling, or a measured pilot budget |
| Contiguous distributed blocks encode source order | Deterministic permutation with retained row IDs and repeated seeds |
| Flexible graph/projection grid becomes label-tuned | Use staged Pareto screening and cap finalists before labels |
| Literature comparisons use ISOLET subsets or label-tuned grids | Report dataset size, known-cluster assumptions, and tuning protocol beside every comparison |

## Deliverables

- [ ] Frozen raw617 parity config and dataset hashes.
- [ ] ISOLET graph-forensics report resolving the 22,471-edge fingerprint.
- [ ] Published correlation resolver surface and tests.
- [ ] Purity-coverage and hierarchy evaluator isolated from training surfaces.
- [ ] Deterministic shuffled-partition support for distributed SPRED, with row-ID recovery.
- [ ] Stable-modality feature selector and historical Lilliefors control.
- [ ] Laplacian-score selector and redundancy pruning.
- [ ] Matched PCA and random-projection controls.
- [ ] H0 SPRED benchmark at 20, 30, and 50 dimensions.
- [ ] Conditional H0/H1 and correspondence-aware SPRED experiments.
- [ ] Modern graph/coupling comparison on the small finalist set.
- [ ] Final label-sealed report with stochastic uncertainty and compute accounting.

## Immediate execution order

1. Correct the benchmark oracle and reporting language to 93% purity / 65% coverage plus hierarchy.
2. Audit the raw mutual-neighbor graph until the 22,471-edge discrepancy is understood.
3. Exercise the exact published correlation resolver independently of current lineage resolution.
4. Freeze and run the raw617 parity configuration.
5. Build the label-sealed purity-coverage/hierarchy report surface.
6. Run the cheap representation screen.
7. Run H0 distributed SPRED on raw and the winning selected-feature input.
8. Promote only unlabeled-stable candidates to graph improvements and final label scoring.

This order prevents SPRED from becoming an expensive explanation for a parity defect and gives any
eventual gain a clean causal interpretation.

## Primary references

- Blatt, Wiseman, and Domany, "Data clustering using a model granular magnet":
  https://arxiv.org/abs/cond-mat/9702072
- UCI ISOLET dataset:
  https://archive.ics.uci.edu/dataset/54/isolet
- Fanty and Cole, "Spoken Letter Recognition":
  https://proceedings.neurips.cc/paper/1990/hash/49182f81e6a13cf5eaa496d51fea6406-Abstract.html
- Yu and You, "Shape-Preserving Dimensionality Reduction":
  https://arxiv.org/abs/2106.02096
- Quiroga, Nadasdy, and Ben-Shaul, "Unsupervised spike detection and sorting with wavelets and
  superparamagnetic clustering":
  https://pubmed.ncbi.nlm.nih.gov/15228749/

