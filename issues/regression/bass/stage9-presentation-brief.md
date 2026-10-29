# BASS stage 9: the presentation data (chip brief)

Date: 2026-10-29. Stages 4–6 are on `main` and verified: the outputs checklist's reports, the last being "O7 verified: stage 6 complete". This chip builds stage 9's maths side, the outputs checklist's O9: the data R's `plot`, `print` and `summary` methods show (the hub's P09–P15), computed in `Maths.Regression`, with no rendering dependency. It also closes two small stage-6 leftovers.
- *What exists.*
  - `BassFit`, with its draws, `BassDiagnostics` (R̂ and ESS, acceptance, refusals, the tempering readout), `BassEvidence` and `BassResolved`, and the views `Predict`, `PointwiseLogDensity`, `ResidualDependenceFor`, `Sensitivity` and `Effects`.
  - `BassBasisFit`, with `Predict`, `PointwiseLogDensity` and `Sensitivity`.
  - `SobolAnalysis` and `BasisSobolAnalysis`, whose indices have posterior means but no quantiles.
  - Two weighted quantiles already: `DrawSample.Quantile` in `Maths.Samplers.Ensemble` (`src/maths/samplers/ensemble/Readouts.cs`), the shared readout, and `PredictiveMixture.Quantile(values, weights, p)`, which stage 5 wrote in `bass/posterior/`.
- *This chip builds:*
  1. *Posterior summaries:* the quantiles R's `summary` and `plot` report, weighted, through one shared weighted quantile.
  2. *Diagnostic data,* R's `plot.bass` and `plot.bassBasis`: the traces of M and σ² (and τ, λ) per chain; fitted against observed with intervals; residuals; and for basis fits the same in physical units.
  3. *Model metadata,* R's `print` and `summary`: the data's shape, the config and its resolved auto-picks, the run's identity and lengths, the moves' acceptance and refusals, the evidence and the tempering readout.
  4. *Sensitivity summaries,* R's `plot.bassSob`: each index's posterior mean and quantiles, for direct and basis analyses, per location and integrated.
  5. *Effect functions for basis fits* (O7's report: nearly free, since an effect is linear and has no cross terms), the gap stage 6 left.
- *Not in scope.*
  - Rendering and text formatting. They belong above maths (the hub's boundary rules). No consumer above maths takes BASS yet, and viz v2's contract shapes are still open, so an adapter would bind them early. Name in your report what an adapter would need.
  - F-M2 (the local moves), whose Phase 1 a later chip builds.
  - Calibration (stage 7) and SMC (stage 8).

Read first: [AGENTS.md](../../../AGENTS.md) and [CONTRIBUTING.md](../../../CONTRIBUTING.md). Then:
- the outputs segment's "Presentation (stage 9)" and D07's "Intervals" bullet: quantiles of the predictive mixture, never of per-draw interval endpoints;
- the hub's scope table (P09–P15) and its stage 9 row;
- R's methods, for their purposes: `plot.bass`, `plot.bassBasis` and `plot.bassSob` in `D:/aghado01/packages/BASS/R/plot.R`; `print.bass`, `print.bassBasis`, `summary.bass` and `summary.bassBasis` in `generics.R`. Cover what each is for, not its layout;
- the three outputs-segment chips' reports (stage 5, stage 6, O7), for the views you summarize;
- `docs/Maths/BayesianRegression.md`, BASS's sections.

In the code:
- `BassFit.cs`, `BassDiagnostics.cs`, `BassEvidence.cs` and `BassResolved.cs`;
- `bass/posterior/` (`BassPrediction`, `PredictiveMixture`) and `bass/output/BassComponents.cs`;
- `bass/sensitivity/` (`SobolDecomposition.cs`, `BassBasisSensitivity.cs`, `EffectFunctions.cs`, `BassSensitivity.cs`);
- `src/maths/samplers/ensemble/Readouts.cs` (`DrawSample`, `ChainTrace`), and `src/maths/regression/PosteriorSummaries.cs`, BARS's summaries, for the house pattern.

## Working rules (the shared checkout)

- **Work in `D:\aghado01\ThermoMapper` on `main`.** Worktrees fail on this repo's long paths. Switch no branch, push nothing.
- **Build and test into your own artifacts root:** add `-p:ArtifactsRoot=D:\aghado01\ThermoMapper\artifacts\chip-stage9\` to every `dotnet build` and `dotnet test`. Do not run `scripts/fact-harness.ps1`. Delete `artifacts\chip-stage9\` when you finish.
- **Others in the checkout.** The main session owns `docs/Maths`, the briefs and checklists. A later F-M2 chip will own the moves: `BassSweep.cs`, `BassTermProposal.cs`, `BassConfig.cs`, `BassRun.cs` and `BassDiagnostics.cs`. Read those, but do not edit them. Expect others' commits between yours.
- **Your files.**
  - New code in `src/maths/regression/spline/bass/presentation/`.
  - View methods on `BassFit.cs` and `BassComponents.cs` (`BassBasisFit`).
  - Additions in `bass/sensitivity/` for the quantiles and basis effects. Keep existing results bit-identical, and rerun the sensitivity tests.
  - Tests `tests/regression/BassPresentation*Tests.cs`, and basis-effect tests `tests/regression/BassSensitivityBasisEffect*Tests.cs`.
  - This brief's Reports section.
- **One weighted quantile.**
  - Use `DrawSample`, the shared readout, for posterior quantiles.
  - Compare its convention with `PredictiveMixture.Quantile`'s weighted empirical quantile. If they differ, keep one: a change to `PredictiveMixture` is yours, a change to `DrawSample` is shared code.
  - Name the change, its effect on any existing result, and the suites you reran.
- **Shared code stays generic.** A change in `Maths.Samplers` is tested there and named in your report.
- **Commit your own files only, by explicit path** (`git add <paths>`, never `-A` or `.`), one commit per task, in the repo's message style, ending with `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`.
- **Line endings:** briefs and C# are LF.
- **Decisions are not yours to settle.** Where the decisions leave something open (which quantiles by default, a record's shape), implement a reversible recommendation and put the question in your report with your reasons.
- **Values, not formatting.** Every summary is a record of values with explicit axes and names. No `string` tables, and no `ToString` layouts standing in for an adapter.

## Checks

- *Summaries against direct computation:* each quantile and mean from the draws, under equal and unequal weights; the traces against the draw stores; fitted and residual data against `Predict` at the training inputs.
- *Intervals* come from the predictive mixture (D07), never per-draw endpoints. O-O1's interval fixture already shows R's difference; reuse it.
- *Basis effects:* a one-component basis fit reproduces the direct fit's effects scaled by y_s(ℓ)·K_ℓ1; tensor quadrature of the reconstructed output agrees.
- *Sensitivity quantiles* agree with `DrawSample` over the defined draws, with zero-variance draws left out.
- *Metadata* round-trips what the fit holds: config, resolved values and identity.

## Review before reporting

Before you write the report, take a fresh-eyes pass and fix what it finds:
- *Is it correct?* Weights everywhere, axes and names explicit, nothing recomputed differently from the views it summarizes.
- *What is missing?* An R method's purpose not covered, a summary without a test.

Name the passes and their fixes in the report.

## Report

Under "Reports" below, newest first:
- what landed, with the commits;
- the tests and suite counts;
- any shared-code change;
- the questions you left as reversible defaults, with your reasons;
- the review passes;
- anything the main session must know: the docs to write, and what a rendering adapter above maths would need.

## Reports

Newest first.

### O9 built: the presentation data — 2026-10-29

**What landed.** Four commits, the presentation records in namespace `Maths.Regression.Spline.Bass.Presentation` (`bass/presentation/`), the basis effects in `...Bass.Sensitivity`.
- `6d37bab` fix(samplers): `DrawSample`'s quantile is type 7 under any equal weights (shared code, below).
- `c564f56` refactor(bass): one weighted quantile.
  - `DrawSample.Quantile` is the one weighted quantile.
  - `PredictiveMixture.Quantile(values, weights, p)` is deleted. It was the step CDF's inverse (R's type 1, the least value whose cumulative weight reaches p).
  - `BassBasisPrediction.Quantile` now takes a `DrawSample` per case and location.
  - The Gaussian mixture's exact quantile, which direct fits use, is unchanged.
- `f6652b8` feat(bass): effect functions for basis fits, `BassBasisEffects.Of` and `BassBasisFit.Effects(spec, effects)`, returning a `BasisEffectAnalysis`.
  - *How it is built.* Each component's effects come from `BassSensitivity.Effects` on its coefficient response, with the shared set of effects named, and are read at the joint draws. `BassSensitivity.cs` is untouched.
  - *The state* is kept free of ℓ. The views are:
    - `ComponentEffects(k)`, in coefficient units;
    - `AtLocation(ℓ)` and `PerLocation()`, with effects Σₖ g_ℓk·eₖ,u and mean y_m(ℓ) + g_ℓᵀE[w], where g_ℓk = y_s(ℓ)·K_ℓk;
    - `Integrated()`, with the location as an input, equally weighted: input effects weighted by ḡ; the location's main effect y_m(ℓ) + g_ℓᵀE[w] − f₀, y_m included as in O7; and each one-input effect with the location, Σₖ (g_ℓk − ḡₖ)·eₖ,v, the location's index fastest.
  - *The default set* is every main effect, plus each pair that some component's joint-draw terms hold.
  - *Refused:* output coordinates and output measures, as basis sensitivity refuses them; a cell budget over r × joint draws × cells; the integrated form's own cells, checked on demand; and the direct path's own refusals.
- `4892ae2` feat(bass): the presentation data, with the views on the fits.

  | R's method and purpose | The data | The view |
  |---|---|---|
  | `plot.bass`: traces of M and σ² | `BassTraces`, with one `BassChainTrace` per chain: `Iterations` (each draw's sweep after burn-in, (i·k + 1)·thin), `Count`, `Sigma2`, `Tau`, `Lambda` | `BassFit.Traces()` |
  | `plot.bass`: fitted against observed, intervals, residual histogram with its normal curve | `BassTrainingFit`: `Observed`, `Mean`, `Quantiles` (probability × case × location), `Residuals`, `NoiseSd` (E[σ], R's `mean(s)`), `LocationCoordinates` for a direct functional fit, `Intervals` = `PredictiveMixture` | `BassFit.TrainingFit(data, probabilities)`, through `Predict` at the training cases |
  | `plot.bassBasis`: traces per component; fit and residuals in physical units | the same records: a list of `BassTraces`, and a `BassTrainingFit` with `Intervals` = `PredictiveDraws`, `Observed` reconstructed exactly from the decomposition | `BassBasisFit.Traces()`, `TrainingFit(probabilities, seed)` |
  | `print.bass`, `print.bassBasis` | `BassMetadata`: `BassDataShape` (cases, locations, predictors' `RangeScale`/`TrainingLevels`, output variables, the categorical count), and the fit's own `Config`, `Resolved` (with `AutoPicks`), `Identity`, `Diagnostics` (convergence, moves, refusals, support, tempering) and `Evidence`, beside `BassRunLengths`. `BassBasisMetadata` adds the PCA account, `BassTruncationSummary`, the master seed, the joint draws and each component's metadata | `Metadata()` on both |
  | `summary.bass`, `summary.bassBasis` | `BassSummary` / `BassBasisSummary`: the metadata plus `BassPosteriorSummary` (M, σ², τ when sampled, λ). `Count.Min`/`Max` is R's range; `Sigma2.Mean` is R's posterior mean error variance | `Summary(probabilities)` on both |
  | `plot.bassSob`: boxplots of S and T; functional panels of mean S and mean S.var | `SobolSummary`: each subset's `MainIndex` and each input's `TotalIndex` over the defined draws; `MainVariance`, `TotalVariance` and `Variance` over every draw; the draw, defined, clamped and negative counts. `SobolAnalysisSummary` holds one per decomposition (per coordinate, or integrated), with the coordinates and the output variables' names. `BasisSobolSummary` holds `PerLocation` and `Integrated` | `BassFit.SensitivitySummary(data, spec, probabilities)`, `BassBasisFit.SensitivitySummary(spec, probabilities)`; `SobolSummary.Of(decomposition)` and `BasisSobolSummary.Of(analysis)` for analyses already computed |

  - Every scalar summary is a `PosteriorSummary(Mean, Sd, Min, Max, Probabilities, Quantiles, Draws)`, under the weights, through `DrawSample`.
  - A box's whiskers at R's `range = 0` are `Min` and `Max`.

**Tests.**
- New in Maths.Regression.Tests, 23 cases:
  - `BassPresentationTests` (10);
  - `BassPresentationBasisTests` (4);
  - `BassPresentationSensitivityTests` (4);
  - `BassSensitivityBasisEffectTests` (5).
- New in Maths.Samplers.Tests: `ReadoutTests.DrawSample_EqualWeights_AtAnyScale_GiveTheType7Quantile`.
- Changed:
  - `ReadoutTests.DrawSample_Weighted_PoolsTiesAndDropsZeroMass`: (0.5, 0.5, 0) now gives 10.5 at p = 0.05, inside the positive mass as before.
  - `BassPredictTests.Quantile_IsThePredictiveMixtures` drops the deleted empirical form's four asserts.
  - `BassPredictBasisTests.Quantile_IsTheWeightedDrawsQuantile` asserts type 7.
- The checks, as the brief lists them:
  - *Summaries against direct computation.* Means, sds, extremes and quantiles from the draws, under equal, normalized-equal and unequal weights (with zeros and ties). The quantiles are checked against a by-hand type 7 and a by-hand weighted-midpoint ECDF, and the summaries of M, σ², τ and λ are pooled over the chains.
  - *Traces* against the draw stores, including a store keeping one in three.
  - *Fitted and residual data* against `Predict` at the training inputs (scalar, categorical, functional; equal and unequal weights), and against basis `Predict` in physical units. The observed values match the PCA's Y to 10⁻¹².
  - *Intervals.* O-O1's construction: 1000 intercept-only draws whose means, or whose sds, spread. The training fit gives the mixture's ±2.771 and ±2.644, equal to `PredictiveMixture.Quantile` to 10⁻¹², against R's per-draw endpoints of ±3.91 and ±5.20, the fixture's numbers. The R-run fixture itself stays in the oracle suite, which pins `PredictiveMixture` to R's mixture.
  - *Basis effects.*
    - One component reproduces its direct effects scaled by y_s(ℓ)·K_ℓ1 (10⁻¹³), and the mean y_m(ℓ) + g·E[w].
    - Tensor quadrature of the reconstructed latent output (160² midpoint nodes, through `BassBasisFit.Predict`) agrees per location and for the integrated form's six effects and mean to 1.2·10⁻¹⁴, on effects up to 1.82.
    - Resampled joint draws are followed.
  - *Sensitivity quantiles* agree with `DrawSample` over the defined draws. Two zero-variance draws are left out of the indices and kept in the variances. The means are the decomposition's own, bit for bit. This holds for direct, conditioned (three coordinates, one outside the range), per-location and integrated analyses.
  - *Metadata* round-trips the config, resolution, auto-picks, identity, diagnostics and evidence (a tempered run with a reference rung: ln p(y) = −40.4 ± 0.5) by reference, and the shapes and lengths by value.
- Suites, on the committed tree:
  - Maths.Regression.Tests: 943/943 (920 before, plus 23).
  - Maths.Samplers.Tests: 187/187 (186 before, plus 1).
  - The BASS oracle suites, Maths.Oracle.Tests filtered to `Bass`: 54/54, with R present.

**Shared-code change:** `DrawSample.Quantile` (`src/maths/samplers/ensemble/Readouts.cs`).
- *Before:* unit weights gave type 7, and every other weighting the midpoint rule, so equal weights of 1/n gave the Hazen-type quantile. That broke the class's own stated invariance to a common rescaling, and normalized MCMC weights would have moved a fit's quantiles off type 7.
- *Now:* zero weights are dropped first. Equal positive weights at any value give type 7, and unequal ones keep the midpoint rule.
- *Effect on existing results:* none. The one other caller, BARS's `SummarisePeak`, adds unit weights. The only expectation that moved is the shared test above.
- *Not shared, but existing results that moved:* basis predictive quantiles, from type 1 to type 7, within their Monte Carlo error. The sensitivity code is untouched, so its results are bit-identical; its tests ran in the suite.

**Questions left as reversible defaults.**
1. *Default probabilities* are 0.025, 0.25, 0.5, 0.75 and 0.975, for every summary and for the training fit's intervals. One set covers R's plot interval (0.025, 0.975), the 50% box and the median. R's boxplot extremes are `Min` and `Max`.
2. *Near-equal weights.* `DrawSample`'s equal-weight rule leaves a jump between equal and near-equal weights (type 7 against the midpoint rule), inherent to pairing the two. A weighted type-7 generalization, through Kish's effective size, is the alternative. It matters once SMC (stage 8) brings unequal weights.
3. *Which draws each summary reads.* Indices are summarized over the defined draws, variances over every draw, since a zero variance is a value. The alternative is both over the defined draws.
4. *Basis training-fit intervals* come from the predictive draws: Monte Carlo, one per joint draw, from `seed` (0 by default). An exact quantile of the draws × residuals mixture costs S·n per evaluation and is not built (stage 5's default 2).
5. *A basis fit's observed responses* are reconstructed from the decomposition: y_m + y_s∘(Kw + e), with the in-sample residual, or the supplied one when the residuals were supplied. The fit does not store Y.
6. *`NoiseSd` is E[σ],* R's code (`mean(s)`). R's docs say "posterior mean variance"; √E[σ²] is the alternative.
7. *The trace axis is sweeps after burn-in,* a ladder's scan counting as one. R's x-axis is the retained index; the adapter can use either.
8. *The metadata holds the fit's own records* (`Config`, `Resolved`, `Diagnostics`, `Evidence`), so it round-trips by identity. `Diagnostics`' tallies are mutable classes, so the metadata is as immutable as the fit, no more.
9. *Basis effects:*
   - the default set is the union over the components;
   - the integrated form gives (v, location) for each one-input effect, and no input-pair × location three-way;
   - mixed hinge degrees are allowed, since an effect has no cross-component terms, though basis sensitivity refuses them.
10. *The namespace* is a sub-namespace, `...Bass.Presentation`, as `Sensitivity` is, not the parent namespace that `posterior/` and `output/` use.
11. *Not built:*
    - summaries of the effect functions (mean and bands per cell). R plots none; an adapter can reduce each cell with `PosteriorSummary.Of`.
    - a latent credible band in the training fit. R's panel shows the predictive interval only.

**Review passes.**
- *Correctness.*
  - **Weights.** Each view reads the right weights: the fit's (direct), each component's own (basis summaries), and the joint draws' (basis prediction, sensitivity, effects). This pass found the `DrawSample` scale defect, fixed in shared code, rather than in each caller passing unit weights.
  - **Recomputation.** `SobolSummary` overwrote its mean with the decomposition's mean, which would have left the sd about another mean had they diverged. Now there is one computation whose sum is the decomposition's own, asserted bit for bit.
  - **Names and axes** are explicit on every record. Each record's doc names its axes (case × location, probability × case × location, draws × cells with the last input fastest).
- *Coverage.*
  - Basis effects under resampled joint draws were untested; the case was added.
  - Each R method's purpose is mapped in the table above, and none is uncovered.
  - Commits stand alone: the effect test had borrowed a helper from the presentation tests, and it is now inlined.
  - A `MainIndexAcrossLocations` convenience was dropped, a projection an adapter makes.

**For the main session.**
- *Docs to write:*
  - In `docs/Maths/BayesianRegression.md`:
    - a "BASS: presentation" section covering the views, the records and their axes, the intervals (D07) and the defaults above;
    - the basis effects' formulas in the sensitivity section;
    - the API and source-map rows;
    - line 1037, "quantiles wait for stage 9", is now done;
    - line 982, basis intervals, now type 7 through `DrawSample`.
  - In `docs/Maths/Samplers.md:351`, `DrawSample`'s rule: type 7 under equal weights at any scale.
- *Checklist:* O9 can be ticked (`4892ae2`), closing the two stage-6 leftovers: the indices' quantiles, and basis effects (`f6652b8`).
- *What a rendering adapter above maths needs:*
  - **Trace panels.** Line series from each `BassChainTrace` (x = `Iterations`), one panel set per component for a basis fit.
  - **Fit panel.**
    - Scatter `Observed` against `Mean` per case and location, with interval segments from the two `Quantiles` rows it picks by probability, plus the 1:1 line.
    - The residual histogram bins `Residuals` (binning is the adapter's), overlaying N(0, `NoiseSd`²) when it is not null.
    - `Intervals` names the method, for a label.
  - **Print and summary text.** These format `BassMetadata` and `BassSummary`. `BassConfig.ToJson()` stands in for R's call.
  - **Sensitivity panels.**
    - Boxes from `MainIndex[k]`/`TotalIndex[v]`: `Min`, the quantiles, `Max`. Labels come from `Inputs` and the `Subsets` index lists.
    - Functional panels stack the mean main indices and mean main variances in subset order, against the conditioned form's `Coordinates` (names in `OutputVariables`).
    - A basis fit's locations are indices only: the fit holds no output coordinates, so the adapter takes the caller's grid.
  - **Effects.** Each `Effect` is draws × cells over its `Grids`, the last input fastest; bands are the adapter's reduction.
  - **Binding.** All of these are arrays and records, with no viz types. Binding them to viz v2's contracts waits until those are decided.
