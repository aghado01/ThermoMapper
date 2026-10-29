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
