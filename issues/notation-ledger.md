# Notation ledger

Splits in the module's language (symbols, identifiers, terms), logged as sessions find them. They wait for one unification pass across docs and code, run after the BASS regression work and its extensions are implemented, so every name is judged against the whole vocabulary at once. Entries record a split and where each reading lives; renames wait for that pass.

**Adding an entry.** Prepend it under `## Entries`, newest first, with:

- **Kind:** one concept under several names, or one name for several concepts.
- **Readings:** each meaning or spelling, with anchors. Code takes `file:line`; docs take a section or phrase, since their line numbers drift.
- **Found:** the date and the work that turned it up.
- **Note** (optional): what the pass needs to know, such as the sources' own symbol.

## Entries

### The gauge factor: g, κ, 1/K̂ and the old couplingScale (2026-11-01)

- **Kind:** one concept under several names.
- **Readings:** `κ` ("Normalization κ") in `docs/Graphs/GraphConstruction.md` before the compiler revision, `g` there now, and `GaugeFactor` in `GaugeScaler`, `ScalerResult` and `CouplingResolution`; `1/K̂` in the BWD papers, `BwdPottsCriticalEstimate`'s comments and `issues/spc-BARS/graph-compiler-brief.md`; `couplingScale` (`n / totalNeighbors`) in `GlobalBandwidthScaler`, removed by the revision.
- **Found:** the graph-compiler revision.
- **Note:** K̂ is the mean degree of the returned graph, so g = 1/K̂ only for `GaugeKind.InverseMeanDegree`.

### Bridge, MST bridge and repair edge (2026-11-01)

- **Kind:** one concept under several names.
- **Readings:** "bridge" and "MST bridge" in the compiler brief, `MstBridge`, `MstBridgeReport`, `EdgeWeightSummary.NearZeroBridges` and `ConnectivityRepair`'s "minimal bridges"; "repair edge" in `RepairEdgeSet` and the build log.
- **Found:** the graph-compiler revision.
- **Note:** a repair edge is not always a bridge in the graph-theoretic sense (a cut edge), and a graph-theoretic bridge need not be one the repair stage added.

### Near-zero: three thresholds (2026-11-01)

- **Kind:** one name for several concepts.
- **Readings:** the coupling floor, 1e-10, under which an edge is dropped (`CouplingProjection.CouplingFloor`); the collapse level, 1e-8, which the auto-LMP trigger and `MaxNearZeroEdgeRatio` count against (`AutoThresholds.NearZeroCoupling`); `EdgeWeights.Summary`'s `nearZeroEpsilon`, 1e-10, behind `NearZeroBridges`.
- **Found:** the graph-compiler revision.
- **Note:** only the first two are config fields; the third is a parameter of the summary.

### Pair connectivity: τ, ⟨n⟩, P(i↔j) and C, beside δ̄ (2026-11-01)

- **Kind:** one concept under several names, with a second concept one affine map away.
- **Readings:** τ_ij and τ_e in the stopping primer (§12.21, §13.6) and `issues/spc/bayesian-clustering-design.md`; ⟨n_ij⟩ in `CoMembership` (`src/graphs/primitives/CoMembership.cs:4`); P(i↔j) in `issues/spc/opus-gvdv-ideas.md`. δ̄ = ((q−1)τ + 1)/q, the spin agreement, is the different quantity `DenseTStack` cuts, read off co-membership by `SweepEdgeCurves.CoMembershipDelta`.
- **Found:** reading θ as a cost ratio for the Bayesian clustering note. The ratio belongs to τ, so θ on δ̄ maps through (qθ − 1)/(q − 1).
- **Note:** τ also names `LineagePersistence`'s split threshold in `docs/Clustering/SpcAnalysis.md`, so that doc and the closure gap's docstrings write ⟨n_ij⟩, and the design note writes τ, as the primer does.

### Binder: the loss and the cumulant (2026-11-01)

- **Kind:** one name, two concepts by two authors.
- **Readings:** Binder's loss (D. A. Binder 1978), the expected pair-counting loss with costs a and b, in `issues/spc/opus-gvdv-ideas.md`; the Binder cumulant U = 1 − ⟨m⁴⟩/(3⟨m²⟩²) (K. Binder 1981), `BinderCumulant` (`src/graphs/models/potts/observables/BinderCumulant.cs:16`).
- **Found:** naming the loss for the Bayesian clustering note.
- **Note:** proposed, not decided: pair loss over all pairs, link loss over graph edges (`issues/spc/bayesian-clustering-design.md`). Each would keep "Binder's loss" as its citation only.

### E: the energy and the window's width (2026-10-31)

- **Kind:** one symbol, two concepts in one doc.
- **Readings:** the window's width in e-folds in the stopping primer's §3 ("over a window of E e-folds", n_min = 1 + ⌈2E/w_min⌉) and §13.4 (n ≈ c·E/w_min); the energy in its §5, §12 and §13 (⟨E⟩, Var(E), σ_E, D_J) and in `Accumulator.EnergyPerDraw` (`src/clustering/graphical/spc/runtime/core/sampler/Accumulator.cs`).
- **Found:** writing the primer's round 2 (`issues/spc-BARS/stopping-theory-primer.md`, §13).
- **Note:** §13 keeps §3's E for the window where it quotes §3's count, and the energy everywhere else.

### Repeated runs at one temperature: replica and replicate (2026-10-31)

- **Kind:** one concept, two names.
- **Readings:** "replica" in `SweepCurve.Replicas` (`src/clustering/graphical/spc/profiling/SweepCurves.cs:115`), in `MonteCarloError`'s "the replicas as the chains" (`src/clustering/graphical/spc/profiling/MonteCarloError.cs:10`) and in the stopping primer's §12 and §13 ("8 replicas"); "replicate" in the primer's §1–§2 ("replicate spread"), in `Pointwise.cs`'s "a replicate's observation" (`src/maths/regression/Pointwise.cs:13`) and in the BASS extensions brief's "replicate sweeps" (`issues/regression/bass/extensions-brief.md`, the conformal entries).
- **Found:** the primer's round 2, where the user's own word was "replicant".
- **Note:** physics' replica (an independent copy of the system, as in replica exchange) and statistics' replicate (a repeated measurement) both fit. §13 says "replica" for a run and keeps "replicate sweep" where it follows the extensions brief.

### One SW update: draw, cycle and sweep (2026-10-29)

- **Kind:** one concept under three names, one of which ("sweep") also names a second concept.
- **Readings:** "draw" in `SwendsenWang.Draw` (`src/clustering/graphical/spc/runtime/core/sampler/SwendsenWang.cs:246`), `Accumulator.DrawCount` and the per-draw record `Accumulator.EnergyPerDraw` and siblings (`Accumulator.cs:29`, `:56`); "cycle" in `RunBudget.Cycles` (`RunBudget.cs:9`), `BurnIn(int cycles)` (`SwendsenWang.cs:214`) and the engine's `Σ_cycles` comments; "sweep" in the SW retention briefs' "per-sweep record" and "measured sweep" (`issues/spc/sw-retention-brief.md`, S1). "Sweep" also names the temperature sweep: `FixedGridSweepStrategy`, `SweepProfile`, `SweepCurves`, `SweepResult`.
- **Found:** naming the per-draw record for the SW retention chip (S1).
- **Note:** the record took "draw" (`…PerDraw`), after `DrawCount`.

### Coefficients: `Beta` beside `Coefficients` (2026-10-29)

- **Kind:** one concept, two identifiers on one type.
- **Readings:** `BassDesignSnapshot.Beta` and the field `_beta` (`src/maths/regression/spline/bass/BassDesign.cs:30`, `:113`); the public `BassDesign.Coefficients` (`:165`).
- **Found:** listing β's readings for this ledger.

### b_k beside the inverse temperature b (2026-10-29)

- **Kind:** one symbol, two concepts in one doc.
- **Readings:** the move schedule's birth and death probabilities b_k, d_k (r_k the rest) in BayesianRegression.md's "DMGK schedule"; the inverse temperature b in the same doc's tempering sections ("raised to an inverse temperature", "Tempered marginal (step 3.1)").
- **Found:** checking the doc's b uses for the f5ba714 question.
- **Note:** b_k, d_k is Green's (1995) own notation, which DMGK keep.

### The inverse temperature: b, β and T (2026-10-29)

- **Kind:** one concept, several names.
- **Readings:**
  - b / `inverseTemperature`: the sampler contract (`src/maths/samplers/targets/SamplingTarget.cs:7`, Samplers.md's target and temperature bullets) and the regression code and docs (`CollapsedGPrior.cs:47`, `ITemperedMarginal.cs:19`, `BassTarget.cs:194`, BayesianRegression.md).
  - β / `Betas` / `betaMin` / `temperBetaMin`: the tempering ladder (`ParallelTempering`, `TemperingLadder`, Samplers.md's tempering section) and `BarsEnsemble`, until f5ba714 moved them to b / `InverseTemperatures` / `hotEnd` / `temperHotEnd`. That commit stands until this pass, which may revert it; `tests/regression/BassMixingExitLadderTests.cs` was written against the new names after it.
  - T, the temperature itself: SPC's temperature grid and the Potts/SW code, with β = 1/T in `src/clustering/dendrogram/PartitionHierarchy.cs:101` and the SPC notes (`issues/spc/opus-spin-glass-threads.md`, `issues/spc/sw-retention-brief.md`).
- **Found:** 2026-10-28, writing E13 in the extensions brief; f5ba714.
- **Note:** the sources write β (Watanabe 2013, Farhi et al. 2013, Friel and Pettitt 2008); Bissiri et al. 2016 write w for the same weight. The SPC↔regression bridge (temperature integration, frustrated clustering) puts SPC's T and the samplers' b in one derivation.

### β: four readings (2026-10-29)

- **Kind:** one symbol, several concepts.
- **Readings:**
  - Regression coefficients: BARS, BAPS, HPS and BASS code and BayesianRegression.md; `HierarchicalPSpline.cs:160` (β_r, β_0).
  - The inverse temperature (entry above).
  - A spectral exponent, PSD ∝ |k|^−β: `src/maths/samplers/SpectralNoiseField.cs:11`.
  - The Beta distribution, by name: `AdaptiveLocalBetaKernel`, the `"beta"` kernel strings in `CountLawTests`.
- **Found:** 2026-10-28, the f5ba714 survey.
- **Note:** BayesianRegression.md's "Tempered marginal (step 3.1)" has coefficient β and inverse temperature b in one sentence.

### η: six readings (2026-10-29)

- **Kind:** one symbol, several concepts.
- **Readings:**
  - The GLM linear predictor, η = Bβ: `src/maths/regression/spline/bars/PoissonModel.cs:13`, `WhittleModel.cs:12`, BayesianRegression.md's `PoissonModel` bullet.
  - BGP's noise-to-signal ratio, σ² = s²η: `src/maths/regression/bgp/GpRegression.cs:18`, `BgpSampler.cs:13`, PosteriorContraction.md's prior [P1].
  - BAPS's local smoothing multipliers ηᵢ ~ Gamma(φ, φ): `src/maths/regression/spline/baps/Baps.cs:32`, BayesianRegression.md's `LocallyAdaptiveGibbs`.
  - HPS's parameter-expansion deviations η_r: `HierarchicalPSpline.cs:160`, BayesianRegression.md's `ParameterExpanded`.
  - A step size η_k in the Riemannian subgradient solver: `src/maths/geometry/solver/irls.cs:138`, RiemannianGeometry.md's solver flowchart.
  - A centred factor column of the Khatri–Rao product: `src/maths/linalg/KhatriRao.cs:16`, `BassDesign.cs:449`.
- **Found:** 2026-10-28, checking whether η meant a learning rate anywhere (E13).
- **Note:** the step size is the nearest thing to a learning rate. E13 writes the generalized-Bayes loss weight as b, the inverse temperature it is.
