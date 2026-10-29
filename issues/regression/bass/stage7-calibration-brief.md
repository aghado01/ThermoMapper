# Stage 7: calibration (chip brief)

Date: 2026-10-29. Builds the calibration checklist's C-X2, C-X1, C1, C2 and C3 against D09, which the user decided on 2026-10-29: [port-calibration-brief.md](port-calibration-brief.md), section "D09 (decided 2026-10-29)". Read that section first. Its *(2026-10-29)* amendments are part of the decision, and they govern wherever the earlier text differs.
- *What it builds:*
  - the emulator behind an interface;
  - the declared measurement-error model;
  - the two covariance formulations;
  - variance groups, the discrepancy and categorical calibration parameters;
  - the conditional sampler for θ given one emulator draw, with tempering on D05's likelihood path;
  - the nested conditional reference sampler.
- *Not in scope:*
  - cut SMC and the fitting SMC adapter (stage 8: C4, C5);
  - extensions beyond D09;
  - `docs/Maths`, the checklists and D09's text, which the main session updates after verifying.

Read first: [AGENTS.md](../../../AGENTS.md) and [CONTRIBUTING.md](../../../CONTRIBUTING.md). Then:
- [port-calibration-brief.md](port-calibration-brief.md) whole: the SMC paragraph, D09 and the room for extensions;
- [port-calibration-checklist.md](port-calibration-checklist.md)'s reports on C-O1 and C-O2, which state what each fixture records and what it does not cover;
- [port-outputs-brief.md](port-outputs-brief.md)'s D07: response units, the truncation residuals and covariance in physical units;
- the hub's [placement section](../bass-full-port-brief.md#placement-and-dependency-boundaries): calibration lives in `bass/calibration/`, and fitting never depends on it.

In the code and beside it:
- *R's calibration,* `D:\aghado01\packages\BASS\R\calib.R` (`calibrate.bassBasis` and its helpers), at the pinned commit. Read it for its API, its θ prior and proposals, and its outputs. Its defects are C-O1's fixtures, corrected here and never ported (the no-oracle-bug-parity rule in AGENTS.md).
- *The fixtures:* `r/oracles/bass_calibration_defects_oracle.R` (C-O1) and `r/oracles/bass_covariance_oracle.R` (C-O2), run by `r/scripts/bass_stage1_checks.R`. `tests/oracle/BassOutputDefectsParityTests.cs` shows how an oracle's JSON reaches a C# parity test.
- *BASS's fits:*
  - `BassFit` (a scalar or direct functional fit: `Predict`, its draws);
  - `output/BassBasisFit.cs` (a basis fit: its joint draws, `y_m`, `y_s`, `K`, the component noise variances, D07's truncation residual set);
  - `posterior/` (the predictive pieces).
- *Shared pieces:*
  - `DiagonalPlusLowRank` and `DenseCholesky` in `src/maths/linalg/`;
  - `Maths.Samplers.Blocks`: `BlockScan`, `RandomWalkMetropolis`, `SliceSampler`, `GaussianConditional`;
  - `ParallelTempering` and `TemperingLadder`, for the tempered option;
  - `Maths.Distributions`, which has no half-Cauchy yet: add one beside the others, and name it in the report.

## Working rules (a shared checkout, with a parallel chip)

- **Work in `D:\aghado01\ThermoMapper` on `main`.** Worktrees fail on this repo's long paths. Switch no branch, push nothing.
- **The N14 chip works beside you,** building D10's convergence verdict in `src/maths/samplers/ensemble/` and in BASS's `BassRun.cs`, `BassDiagnostics.cs`, `BassConfig.cs` and `BassResolved.cs`, and BARS's files.
  - Leave those files alone. If a build fails in a file you do not own, wait and retry; never edit it.
  - `BassFit.cs` and `output/BassBasisFit.cs` belong to neither chip. Implement the emulator contract for them by adapters in `calibration/`. If an adapter truly needs a change there, make the smallest accessor that changes no behaviour, and name it in the report.
- **Build and test into your own artifacts root:** add `-p:ArtifactsRoot=D:\aghado01\ThermoMapper\artifacts\chip-stage7\` to every `dotnet build` and `dotnet test`. Do not run `scripts/fact-harness.ps1`. Delete `artifacts\chip-stage7\` when you finish.
- **Your files:**
  - new files under `src/maths/regression/spline/bass/calibration/`, namespace `Maths.Regression.Spline.Bass.Calibration`;
  - a half-Cauchy in `Maths.Distributions`, if it is missing;
  - new tests `tests/regression/BassCalibration*Tests.cs` and `tests/oracle/BassCalibration*ParityTests.cs`;
  - additions to the two R oracles and `bass_stage1_checks.R`, when a corrected reference needs R's side, with every existing check still passing;
  - this brief's Reports section.
- **Commit your own files only, by explicit path,** one commit per step, in the repo's message style, ending with `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`.
  - The user's working-copy edits to `issues/regression/argon-bass-qa.md` and `bass-questions.md` stay unstaged.
  - Before each commit, read `git diff --cached` and confirm every hunk is yours.
- **Line endings:** briefs, C# and R are LF.
- **Configuration is values.** Calibration takes a values-only, validated record that serializes whole, in the style of `BassConfig`: closed sets as records with a `Kind` discriminator, `ToJson` and `FromJson`. It is calibration's own record. `BassConfig` is N14's file this session.

## The steps

Each step lands with its checks and ends with a fresh-eyes review of its own code and tests, for correctness and for coverage gaps, fixed before the next step. Name the passes and their fixes in the report.

**Step 1. The emulator behind an interface (C-X2).**
- *The contract* is what the cut posterior needs:
  - the number of retained draws;
  - for draw s at calibration inputs θ (and any fixed inputs the observations carry), the predicted output in physical units;
  - the draw's emulator-noise covariance in physical units: for a basis fit, (y_s ∘ K)·diag(σ²_{e,k})·(y_s ∘ K)ᵀ as a low-rank piece; for a fit without a basis, σ²_e·I on the diagonal (D09's amendment).
- *Adapters* serve a `BassBasisFit` and a scalar or direct functional `BassFit`.
  - Read from `calib.R` which emulator inputs are θ and which are fixed.
  - Keep the coordinates explicit: θ sampled in [0, 1], reported in physical units, with the transform recorded (D09).
- *Checks:* an adapter's prediction equals the fit's own `Predict` for the same draw, at random inputs, for each fit kind.

**Step 2. The measurement-error model (C-X1),** a declared component: Gaussian with variance groups now, behind a seam that robust models can later fill. Groups are declared by output index.

**Step 3. The covariance formulations (C1).** Declared types, with shapes and positive definiteness validated at construction and a mismatch refused (C-O1's branch-by-shape fixture asserts the correction).
- *Functional:*
  - Σ = diag(σ²_{g(j)}) + the emulator's piece (Step 1) + (1/n)·Σᵢ (y_s ∘ eᵢ)(y_s ∘ eᵢ)ᵀ over D07's n truncation residual vectors, the covariance of D07's uniform draw over them (taken about zero, as that draw is): a low-rank term of rank at most n;
  - with D07's diagonal truncation variance as the option R uses;
  - solved and log-determined through `DiagonalPlusLowRank`, without an m × m inverse.
- *Multivariate:* Σ = Σ_trunc + diag(σ_g)·R·diag(σ_g) + Σ_δ + the emulator's piece, on the dense factor.
- *Checks:*
  - C-O2's dense references, both formulations, at its conditionings;
  - units: scaled and unscaled y give the same physical posterior (C-O2's units section);
  - the low-rank truncation term against its dense form;
  - the diagonal option reproducing R's in physical units.

**Step 4. Variance groups, the discrepancy, categorical parameters (C2).**
- *Variance groups:* each group declares its prior.
  - The default is IG(ν_g/2, (ν_g/2)·ŝ_g²), which needs ν_g > 0.
  - The heavy-tailed choice is a half-Cauchy on σ_g with scale ŝ_g, proper.
  - R's (1 + σ²)⁻¹ kernel is a defect fixture asserted corrected, not an option.
  - Moves on log σ²_g carry their Jacobian.
- *The discrepancy (functional):* δ = μ_δ + D·v, v ~ N(0, s_v²·I), with D, s_v and the optional μ_δ declared. s_v has no default and is refused when missing; s_v = 1 reproduces R. v is drawn from its Gaussian conditional. In the multivariate formulation, Σ_δ and the mean shift are fixed, as D09 has it.
- *Categorical calibration parameters:* the joint discrete conditional, enumerated whole up to 1,024 combinations by default, with each term's categorical part precomputed (Francom et al. 2019). Above the cap, single-site updates take over; measure their mixing on a fixture, without assuming it.
- *Cache consistency:* any change to δ, σ² or the emulator draw invalidates the cached likelihood before the next acceptance test (C-O1's stale-likelihood fixture, corrected).
- *Checks:*
  - C-O1's fixtures, each asserting R's defect as recorded and the correction on this side: the stale likelihood, the branch by shape, the coordinates, the units, and the `s2.df` switch;
  - every conditional against quadrature or enumeration on small fixtures;
  - the half-Cauchy's normalization and scale-equivariance;
  - the enumerated categorical conditional against brute force.

**Step 5. The conditional sampler for θ given one emulator draw,** p(θ, σ², v | y, e), a `BlockScan`:
- θ in [0, 1] coordinates, with the prior and proposal you read from `calib.R`, corrected where C-O1 says so;
- the categorical block;
- the variance groups on log σ², with Jacobians;
- the discrepancy's Gaussian draw.

Tempering is an option on D05's likelihood path, through the shared `ParallelTempering`. C-O1's tempering section is the corrected reference: the kernels reach the likelihood path, and the swap is (b₁ − b₂)(ℓ(x₂) − ℓ(x₁)).
- *Checks:*
  - the chain's laws against enumeration or quadrature on small fixtures, untempered and at each rung;
  - the swap's balance;
  - the propriety rule: a b = 0 rung needs every prior proper.

**Step 6. The nested conditional reference sampler (C3).** This step needs N14's verdict, so first check that N14 has landed: `git log --oneline -- src/maths/samplers/ensemble/`, and `EnsembleRun.Verdict` in the code. N14's brief fixes the names: `ConvergenceVerdict`, `ConvergenceOutcome`, `GroupVerdict`, `EnsembleRun.Verdict` and `ChainEnsemble.ReferenceRHatTarget`.
- *If N14 has not landed,* stop after Step 5 and write your report for Steps 1–5, saying Step 6 waits. The main session will tell you when to continue.
- *The sampler:* for s = 1, …, S, take e_s from the retained emulator draws, and run Step 5's sampler on p(θ | y, e_s) as an ensemble until N14's verdict reads `Converged`, within a budget. Keep its draws, and report each inner run's verdict. An inner run that ends `NotConverged` is reported, never silently pooled.
- *Its error:* vary the inner length, and measure how the pooled posterior's quantities move with it.
- *Recovery:* a synthetic simulator with a known θ*.
  - Fit BASS on its runs, once with a basis and once direct.
  - Simulate field data at θ*, calibrate, and check θ*'s coverage over seeds.
  - Repeat with the discrepancy declared and with two variance groups.
- *Outputs* (D09): θ in physical coordinates, the variances, the discrepancy and the emulator index for each draw, with burn-in, thinning and the verdicts explicit.

## Review before reporting

Take a fresh-eyes pass over the whole change and fix what it finds:
- *Correctness:*
  - every R defect corrected, none ported;
  - units physical throughout;
  - the cut kept: y never updates the emulator;
  - the Jacobians;
  - the cache invalidated on every change;
  - the propriety rule at b = 0.
- *Coverage:*
  - both formulations, both fit kinds, and both truncation options;
  - one and several variance groups, with both priors;
  - with and without the discrepancy;
  - categorical parameters below and above the cap;
  - tempered and untempered.

## Report

Under "Reports" below, newest first:
- what landed, with the commits;
- the API and the configuration;
- what you read from `calib.R` (its θ prior, proposals and outputs) and where this port differs from it, and why;
- each check with its tolerance;
- the recovery runs and the inner-length error;
- any accessor added outside `calibration/`;
- the suites with their counts (Regression, the BASS and calibration oracle classes, R's stage-1 checks);
- the review passes and their fixes;
- anything the main session must know, choices this brief left open included.

## Reports

Newest first.

### Stage 7 built, Steps 1–6 — 2026-10-29

From the stage-7 chip. N14 had landed (`e8547d3`, `EnsembleRun.Verdict`) before Step 6, so all six steps are built. Nothing is pushed.

**Commits.**
- `066a054`: Step 1, the emulator contract (C-X2).
- `3e935b5`: Step 2, the measurement-error model (C-X1).
- `24014cd`: Step 3, the covariance formulations (C1).
- `50565b9`: Step 4, the conditional model, discrepancy and categorical parameters (C2).
- `d3af6a3`: Step 5, the conditional sampler and its tempering.
- `ab14608`: Step 6, the nested reference sampler (C3).

**API** (`src/maths/regression/spline/bass/calibration/`, namespace `Maths.Regression.Spline.Bass.Calibration`).
- *Emulator* (`CalibrationEmulator.cs`, `BassEmulators.cs`):
  - `ICalibrationEmulator` holds the parameters, `DrawCount`, `DrawWeights`, `Outputs`, `Noise(s)` and `TruncationResiduals()`, all in physical units (y_s ∘ eᵢ).
  - `Predict(s, continuous, levels, output)` takes θ in physical coordinates. `PredictLevels` evaluates every level combination in one call.
  - `EmulatorNoise` is `ComponentNoise(Factor = y_s ∘ K, Variances)` or `IndependentNoise(σ²_e)`.
  - `BassFitEmulator(fit, data, parameters, fixedInputs?, outputLocations?)` serves a scalar or direct functional fit. `BassBasisEmulator(basisFit, parameters, fixedInputs?)` serves a basis fit.
  - Each predictor of the fit is either a parameter, matched by name, or a fixed input carried by the observations' cases.
  - Each term's continuous part is evaluated once per θ, and its categorical part once per combination.
- *Parameters* (`CalibrationParameter.cs`): `ContinuousParameter(Name, Lower, Upper)` records the transform: `ToPhysical` and `ToUnit`. `CategoricalParameter(Name, Levels)`. `ParameterLayout` splits them by kind and enumerates the combinations.
- *Measurement* (`MeasurementModel.cs`): `GaussianMeasurement(Groups)` with `VarianceGroup(Outputs, Prior)`. The prior is `InverseGammaVariance(ν, ŝ)` or `HalfCauchyDeviation(ŝ)`. `MeasurementModel` validates the groups as a partition of the outputs.
- *Covariance* (`CalibrationCovariance.cs`):
  - `FunctionalFormulation { Truncation, Discrepancy = FunctionalDiscrepancy(Basis, Scale, Mean) }`.
  - `MultivariateFormulation(Correlation, DiscrepancyCovariance, DiscrepancyMean)`.
  - `TruncationTerm { ResidualSet, Diagonal }`.
  - `CalibrationCovariance.Create(...)` validates, and `.Factor(σ², noise)` returns a `CovarianceFactor` (log-determinant, solve, quadratic form).
- *Model* (`CalibrationModel.cs`):
  - `CalibrationModel(config, y, emulator)` binds the config to the data. `.Target(s)` is a `CalibrationTarget`, an `ISamplingTarget` on the likelihood path, whose `CalibrationEvaluation` carries ℓ₁.
  - `CalibrationState` is immutable: (unit θ, levels, σ², v).
  - `Physical(state)` and `Delta(v)` give θ in physical coordinates and the discrepancy.
- *Blocks* (`CalibrationBlocks.cs`): `DiscrepancyBlock` (with `.Moments(state, b)`), `VarianceBlock`, `CategoricalBlock`.
- *Sampler* (`CalibrationSampler.cs`): `ThetaWalk` (the move), `ThetaAdaptation`, `ThetaDecorrelation`, `CalibrationChain` (an `IChain`), and `CalibrationSampler.Run`, which returns a `ConditionalRun` of `CalibrationDraw`s (emulator draw, physical θ, levels, σ², v, ℓ).
- *Reference* (`NestedCalibration.cs`): `NestedCalibrationSampler.Run(model)` returns a `NestedCalibration`:
  - `Runs` of `InnerRun` (emulator draw, `ConvergenceVerdict`, each chain's draws, burn-in, thinning, kept draws per chain);
  - `NotConverged`;
  - `Pooled(includeNotConverged = false)`, which weights each run's draws 1/(S′·n_s).
- *Configuration*: `CalibrationConfig` is values only, in `BassConfig`'s JSON style, with `ToJson`, `FromJson` and `Validate`. Its fields:
  - `Parameters`, `Measurement` and `Covariance`, all required;
  - `EnumerationCap` (1,024);
  - `Sampling`: seed, `BurnIn` 2,000, `Samples` 2,000, `Thin` 1, `DecorrelationEvery` 100, and `Tempering`, a ladder or null;
  - `Reference`: `OuterDraws` (null takes every draw), `Chains` 4, `MaxSamples` 4,000, `BatchSize` 250, `RHatTarget` 1.01, `MinEssPerChain` 100.
- *Outside `calibration/`*: `Maths.Distributions.HalfCauchy` (C⁺(0, s), with exact tails, quantiles and draws). No accessor was added to `BassFit.cs` or `output/BassBasisFit.cs`: the adapters read existing internals (`BassDrawFit.Pool`, `BassBasisFit.Model`, `PredictionCases.Prepare`, `BassSurface`'s accumulation).

**What `calib.R` does, and where this port differs.**
- *Its API:* `calibrate.bassBasis(y, mod, type, sd.est, s2.df, s2.ind, meas.error.cor, bounds, discrep.mean, discrep.mat, nmcmc = 10000, temperature.ladder = 1.05^(0:30), decor.step.every = 100)`.
- *Its outputs:* θ (iterations × rungs × p, in [0, 1]), s2, discrep.vars, ii (the emulator draws), count, count.decor, tau, curr and dat.
- *The emulator.* Every emulator input is θ: `func` predicts at the transformed θ with one emulator draw, without nugget or truncation error, and R draws that emulator draw afresh, uniformly, every iteration.
  - The port fixes e_s for each inner run of the nested sampler (D09).
  - Its adapters also take fixed inputs per observation case (the brief's "any fixed inputs the observations carry").
- *θ.*
  - R: a uniform prior on the bounds, sampled in [0, 1] from a uniform start. Adaptive Metropolis per rung: S = 10⁻⁶·I until iteration 300, then (2.4²/p)(cov + ε)·e^τ, with τ moved every 100 iterations by min(0.1, 5/√i) towards 23/100. A proposal outside the box is rejected. The decorrelation step is a uniform independence proposal on each coordinate.
  - Kept as R has it: the prior, the walk, its 23% target and scale schedule, the 300-scan switch and the decorrelation step.
  - Changed:
    - the walk starts at 0.1² on each coordinate, since R's 10⁻⁶ hardly moves for 300 iterations;
    - adaptation stops at the end of burn-in, so the kept draws come from one fixed kernel;
    - θ is reported in physical coordinates, through each parameter's recorded transform (C-O1).
- *σ².*
  - R: a joint adaptive walk on log σ² with its Jacobian, from a start at 1. The prior is the IG kernel, or −log(1 + σ²) for every group when any `s2.df` is 0. R tempers the whole acceptance ratio.
  - Here: a random-walk step per group on ln σ², its scale adapted towards 0.44 during burn-in. The Jacobian stays untempered (D05's likelihood path). Each group declares its own proper prior. The start is a prior draw.
- *v.* R draws v from a Gibbs conditional with precision 1 and covariance `discrep.S/itl`, the powered path's. Here v's conditional is the likelihood path's, I/s_v² + b·DᵀΣ⁻¹D, with s_v declared.
- *Tempering.*
  - R tempers by default on 31 rungs and swaps random pairs after iteration 1,000, with a swap ratio that leaves its own target (C-O1).
  - Here tempering is an option with no default ladder. It runs through `ParallelTempering` on the likelihood path, with its even-odd adjacent swaps, whose ratio is (b₁ − b₂)(ℓ(x₂) − ℓ(x₁)).
- *Corrected, not ported:*
  - the stale likelihood: states are immutable and every change is evaluated afresh;
  - the branch by shape: the formulation is declared, and a mismatch is refused;
  - normalized coordinates in the output;
  - R's units: the covariance is in physical units;
  - the all-groups `s2.df` switch and its improper kernel;
  - the tempering path.

**Checks and tolerances** (54 calibration tests in `Maths.Regression.Tests`, 9 parity tests in `Maths.Oracle.Tests`).
- *Step 1:*
  - the adapters' prediction equals each fit's own `Predict` bit for bit, at random θ and draws: a scalar fit with fixed inputs, a categorical θ, a fixed categorical input, a direct functional fit at declared locations, a basis fit, and a basis fit with a categorical θ;
  - `PredictLevels` equals `Predict` per level;
  - the noise and the residuals are in physical units;
  - weights pass through, and declarations that do not match the fit are refused.
- *Step 2:*
  - the IG density against its closed form, to 10⁻¹²;
  - both priors normalized by trapezoid in ln σ², to 10⁻⁹;
  - scale-equivariance, to 10⁻¹¹ (half-Cauchy) and 10⁻¹⁰ (IG);
  - the half-Cauchy's CDF and quantile, to 10⁻¹², with draws within 4 binomial SE;
  - the partition refusals.
- *Step 3, low rank against dense:*
  - functional residual set, ln|Σ| and solves to 10⁻¹¹, and to 10⁻¹⁰ with n > m;
  - the diagonal option, as R's centred sample variance, to 10⁻¹¹;
  - σ²_e·I in both formulations, to 10⁻¹¹;
  - multivariate dense, with both truncation options, to 10⁻¹¹;
  - units: c-scaled outputs move ℓ by −m ln c, to 10⁻⁹.
- *Step 3, C-O2:*
  - well case through the formulation: ln|Σ| 2.8·10⁻¹⁴, solve 6.6·10⁻¹⁵ (tolerance 10⁻¹²);
  - moderate (κ 9·10⁷): 1.0·10⁻⁹ and 6.0·10⁻¹⁰ (tolerance 10⁻⁸);
  - extreme: backward error 7.1·10⁻¹⁶ (tolerance 10⁻¹⁴);
  - multivariate: 5·10⁻¹⁴ and 8·10⁻¹⁵ (tolerance 10⁻¹⁰);
  - units: ℓ against R's physical one 4.6·10⁻¹⁴ (tolerance 10⁻⁹), and posteriors against the physical and standardized ones 1.6·10⁻¹⁶ (tolerance 10⁻¹²); R's units give TV 0.96.
- *Step 3, C-O1:*
  - the branch by shape: a functional formulation with a square basis keeps 20 coefficients, and a multivariate one with a two-column discrepancy covariance is refused;
  - units: ln|Σ| against the dense physical one ≤ 1.1·10⁻¹⁵ (tolerance 10⁻¹⁰), and ≥ 0.447 from R's.
- *Step 4, the target:* against dense Gaussians with every prior written out, to 10⁻¹⁰, in both formulations and with σ²_e·I. The caches keep and rebuild exactly by what changed, counted by emulator calls.
- *Step 4, the conditionals:*
  - v against the joint-Gaussian (Kalman) form at b = 1, 0.35 and 0: mean to 10⁻¹⁰, covariance to 10⁻¹¹, draws within 5 SE;
  - ln σ² against quadrature (5 MCSE plus the grid) at b = 1 with the IG, and at b = 1, 0.4 and 0 with the half-Cauchy. For example, the IG at b = 1 gives −4.795 against −4.800, MCSE 0.005;
  - the categorical enumeration against brute force, to 10⁻¹², with frequencies within 5 SE;
  - single-site updates above the cap reach the enumerated law (5 MCSE + 0.002). Their mixing: acceptance 0.15, integrated autocorrelation 12.4 scans, ESS 4,857 of 60,000.
- *Step 4, the propriety rule:* b = 0 refuses an improper prior.
- *Step 4, R's kernel as a defect fixture:*
  - its mass over [0, X] is ln(1 + X), within 2·10⁻³;
  - the posterior median of σ in the original units is 0.2691 against 0.2386 with the data in units c = 10⁻³ against 10³;
  - the half-Cauchy gives 0.2362 in both.
- *Step 4, C-O1:*
  - the stale likelihood: R's rung 5 is off by 0.674; here cached minus fresh is 0 across R's iteration on every rung;
  - the coordinates: the emulator's calls equal the reported physical θ, inside the bounds;
  - the `s2.df` switch: R's −log1p(σ²) is asserted (10⁻¹⁵). Group 2's IG matches R's corrected kernel to 10⁻¹⁰, and a group with ν = 0 is refused and pointed to the half-Cauchy.
- *Step 5, the chain's laws against quadrature with v integrated out* (5 MCSE plus the grid):
  - untempered with a categorical: E[θ] 0.5364 against 0.5349, MCSE 0.0016;
  - multivariate: 0.5017 against 0.5012;
  - a ladder 1, 0.6, 0.25, 0, every rung matching (b = 1: 0.4739 against 0.4740), with swap acceptance 0.48 and 2,205 round trips.
- *Step 5, the rest:*
  - the swap's ratio, to 10⁻⁹;
  - the propriety rule for a ladder that reaches 0;
  - outputs physical, reproducible and thinned.
- *Step 5, C-O1's tempering section:*
  - the variance move's ratio reaches the likelihood path to 8.5·10⁻¹⁴, and the swap to 1.2·10⁻¹³ (tolerance 10⁻⁹);
  - the discrepancy's covariance matches R's likelihood-path matrix, now emitted by the oracle, to 3.2·10⁻¹⁵ (tolerance 10⁻¹⁰);
  - R's swap gaps are 234 and 140.
- *Step 6:*
  - the pooled cut posterior against the exact equal mixture of the conditionals by quadrature: E[θ] 0.5386 against 0.5340, MCSE 0.0037 (tempered inner runs, all three converged);
  - NotConverged is reported and left out of the pool unless asked;
  - the outer draws, plain and weighted;
  - the config's JSON round trip.

**Recovery** (`BassCalibrationRecoveryTests`).
- *Setup:*
  - η(t; θ) = θ₁·sin(2π(t + θ₂/4)) + θ₂·t/2 on 20 locations, run on a 6 × 6 grid;
  - BASS fitted on a PCA basis (rank 3, 200 joint draws) and directly;
  - field data at θ* = (1.3, 0.4) with noise sd 0.05;
  - the nested sampler with 6 outer draws, 2 chains, R̂ 1.05, ESS 50 per chain and at most 1,500 draws;
  - five seeds per case.
- *Coverage of θ* by the 95% intervals:*
  - basis: 5/5 and 5/5;
  - direct: 4/5 and 5/5. Its θ₁ means are 1.41–1.47, with intervals about ±0.16 wide: the direct emulator is the less accurate one;
  - with a discrepancy 0.1 + 0.1(t − ½) in the field data, declared on [1, t − ½] with s_v = 0.2: 5/5 and 5/5;
  - two groups, with noise sd 0.03 and 0.08, under an IG and a half-Cauchy: 5/5 and 5/5.
- *Convergence:* 3 to 6 of the 6 inner runs converged per seed within that budget, and the rest were reported and not pooled.

**The inner-run error** (`InnerLength_TheErrorItLeaves`). Fixed inner lengths: the ESS gate is set unreachable, and every run is pooled explicitly. Each length runs 2 chains × 3 emulator draws × 3 seeds, against the exact mixture's E[θ] = 0.4911 (unit coordinates, posterior sd 0.057).

| Kept draws per chain | RMS error of the pooled E[θ] | Pooled sd |
|---|---|---|
| 25 | 0.0115 | 0.0489 |
| 100 | 0.0015 | 0.0566 |
| 400 | 0.0023 | 0.0566 |
| 1,600 | 0.0009 | 0.0573 |

At 25 draws per chain the pooled spread is 15% low and the mean off by 0.2 sd. From 100 on, both sit at Monte Carlo noise.

**Suites.**
- Regression (`Maths.Regression.Tests`, whole): 1,101 passed, 0 failed (5 min 8 s).
- BASS and calibration oracle classes (`Maths.Oracle.Tests`, filter `Bass`): 63 passed, 0 failed (4 min 8 s).
- R's stage-1 checks: 97 verdicts pass, none fail, after both oracle additions. The oracles now also emit their inputs at 17 significant digits, and C-O1 its likelihood-path discrepancy covariance.

**Review passes and their fixes.**
- *Step 1:* added a fixed categorical input and weighted draws.
- *Step 2:* nothing to fix.
- *Step 3:*
  - the cache of the concatenated factor became a record swapped whole, since a value tuple can tear under concurrent chains;
  - added the residual set wider than the outputs (n > m) and the multivariate diagonal option.
- *Step 4:*
  - widened the quadrature grid for the half-Cauchy's tails at b = 0;
  - moved the unit-dependence check to units where R's kernel bites;
  - added the target under the multivariate formulation and with σ²_e·I, and the config's JSON round trip;
  - fixed two test bugs of mine: an exact slope comparison, and a call count that missed the check's own evaluation.
- *Step 5:* added the multivariate chain law.
- *Step 6:*
  - the e₂ run of the pooled-posterior test missed the θ gate untempered (R̂ 1.033 within 6,000 draws on a θ–level-coupled fixture), and the test now runs its inner runs on a ladder, the gate unchanged;
  - added the weighted outer draws and the config's round trip.
- *Whole change:*
  - every R defect corrected, none ported;
  - units physical throughout;
  - y never updates the emulator;
  - the Jacobians untempered;
  - every change evaluated afresh;
  - b = 0 needs proper priors;
  - the coverage list met.

**For the main session.**
- *The cost of `DiagonalPlusLowRank` (finding).* Its factor is a Jacobi thin SVD of the m × (r + n) factor on every σ² proposal: 1.4 ms at m = 20 with 39 columns.
  - The residual-set term is now compressed once to its thin SVD's numerical rank (at most min(m, n) columns, Σ unchanged to rounding), which made it 33 µs on the recovery problem, whose residuals are rounding noise.
  - In general the cost per proposal stays O(m·(r + min(m, n))²) with Jacobi's constant.
  - If real sizes make that binding, C-O2's own low-rank route, the Cholesky factor of I + E^{½}UᵀD⁻¹UE^{½}, costs O(m·k² + k³) with a small constant. It could join `DiagonalPlusLowRank` as an option. I did not change the linalg file.
- *Choices the brief left open:*
  - The multivariate truncation term. It defaults to the residual set's second moment about zero, as the functional one does; R's type 2 uses the centred n − 1 covariance. The C-O2 parity test reproduces R's by feeding the centred residuals scaled by √(n/(n − 1)).
  - The diagonal option is R's centred n − 1 sample variance, so it reproduces R. D07's density code reads the 1/n centred one.
  - The multivariate correlation is over the outputs (m × m, R's `meas.error.cor`), where D09's wording says "over variance groups".
  - Above the cap, single-site updates are Metropolis steps to a uniformly chosen other level, as D09's text has it. A per-site Gibbs draw would need no rejection.
  - `Pooled` excludes Undetermined runs as well as NotConverged ones.
  - `ReferenceSampling` defaults to N14's gate (R̂ 1.01, ESS 100 per chain). The recovery tests declare 1.05 and 50 for runtime.
  - There is no default tempering ladder.
- *Not done:* `docs/Maths`, the checklists and D09's text, which are the main session's.
