# BASS stage 6: basis sensitivity, O7 (chip brief)

Date: 2026-10-29. Follows stage 5 ([its brief and report](stage5-outputs-brief.md#reports)) and stage 6's first part ([its brief and report](stage6-sensitivity-brief.md#reports)), both on `main` and verified (the outputs checklist's report "Stage 5 and stage 6's first part verified"). D08, decided on 2026-10-29 ([the outputs segment](port-outputs-brief.md)), fixes the design. This chip builds the outputs checklist's O7, R's `sobolBasis`.
- *What exists.*
  - `BassBasisFit` (`bass/output/BassComponents.cs`): `Decomposition.Map`, a `LinearReconstruction` with y_m `Offsets`, y_s `Scales`, K through `Basis(l, k)` or `BasisMatrix()`, and `Physical(C)` = diag(y_s)·C·diag(y_s); `Components[k]`, a `BassFit` per component, each with its own term table over the same predictors and knots; `ComponentData(k)`; `JointDraws.Draw(k, j)`.
  - Stage 6's `BassMoments` (`bass/sensitivity/`): `Pair(a, b)` takes terms by value, so terms from two components' tables meet in one service, and returns a `TermPairMoments` with the covariance's share on every subset of the shared variables. Also the measures, `ProductMeasure`, `SensitivitySpec`, `SobolDecomposition` and its statuses, and `BassSensitivity.Of` for direct fits.
  - The fixture `r/oracles/bass_sensitivity_defects_oracle.R` (O-O2): its units section, R's `sobolBasis` on a mock two-component fit with y_s from 2 to 5, waits for this chip.
- *This chip builds* basis sensitivity:
  - the variance of the reconstructed output under the input measure, decomposed into Sobol main and total indices per joint posterior draw;
  - the cross-component covariances kept;
  - contributions in output units;
  - a subset budget.
- *Not in scope.*
  - Effect functions for basis fits: R's `sobolBasis` has none. Name in your report whether they come nearly free.
  - Calibration (stage 7).
  - Presentation (stage 9).

Read first: [AGENTS.md](../../../AGENTS.md) and [CONTRIBUTING.md](../../../CONTRIBUTING.md). Then:
- D08 in full, above all its "Basis sensitivity", "Functional responses" and "Zero total variance" bullets;
- D07's "Joint draws across components" and "Covariance in physical units";
- stage 6's report, for its algebra (a term is κ·Π_v φ_v; two terms' components covary only on subsets of their shared variables) and for "What O7 needs from stage 5", which sketches this chip's sum;
- the outputs checklist's report on O-O2, its units section;
- R's `sobolBasis` (`D:/aghado01/packages/BASS/R/bassPCA.R:307`), for its outputs: indices per output location, and `S.var` in decomposition units.

In the code: `bass/output/BassComponents.cs`, `BassDecomposition.cs` and `bass/posterior/BassJointDraws.cs`; `bass/sensitivity/BassSensitivity.cs`, `BassMoments.cs` and `SobolDecomposition.cs`; and stage 6's tests, `tests/regression/BassSensitivity*Tests.cs` and `tests/oracle/BassSensitivityDefectParityTests.cs`.

## Working rules (the shared checkout)

- **Work in `D:\aghado01\ThermoMapper` on `main`.** Worktrees fail on this repo's long paths. Switch no branch, push nothing.
- **Build and test into your own artifacts root:** add `-p:ArtifactsRoot=D:\aghado01\ThermoMapper\artifacts\chip-o7\` to every `dotnet build` and `dotnet test`. Do not run `scripts/fact-harness.ps1`. Delete `artifacts\chip-o7\` when you finish.
- **The main session** owns `docs/Maths`, the briefs and checklists. No other chip is running; expect the main session's commits between yours.
- **Your files.**
  - New code in `bass/sensitivity/` (for example `BassBasisSensitivity.cs`).
  - `BassComponents.cs`, for `BassBasisFit.Sensitivity` only. The basis fit carries its predictors, so the entry needs no data argument.
  - Additions to `SobolDecomposition.cs` and `BassMoments.cs` where the basis case needs them. Keep direct sensitivity's results bit-identical, and rerun its tests.
  - Tests `tests/regression/BassSensitivityBasis*Tests.cs`, and the O-O2 units section in `tests/oracle/BassSensitivityDefectParityTests.cs` or a new `BassSensitivityBasisParityTests.cs`.
  - This brief's Reports section.
  - Leave other production files alone. If you need a change in one, make the smallest behaviour-free change, name it in your report, and rerun the suites that cover it.
- **Commit your own files only, by explicit path** (`git add <paths>`, never `-A` or `.`), one commit per task, in the repo's message style, ending with `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`.
- **Line endings:** briefs and C# are LF.
- **Decisions are not yours to settle.** D08 fixes the design. Where it leaves something open, implement a reversible recommendation and put the question in your report with your reasons. Do not re-open decided choices.
- **No oracle-bug parity.** R's decomposition units are a defect fixture: assert the physical values, and document R's beside them.

## What to build

- *The surface.* At output location ℓ, the reconstructed latent output is y_m(ℓ) + y_s(ℓ)·Σ_k K_ℓk·w_k(x), where w_k(x) = Σ_j a_kj B_kj(x) is component k's surface under joint draw j's paired component draws. The truncation residual and the noise are outside the decomposition, as observation noise is for direct fits.
- *The cross-component covariance.*
  - Var_x of that output is y_s(ℓ)²·Σ_{k,k′} K_ℓk K_ℓk′·Cov(w_k, w_k′), and Cov(w_k, w_k′) = Σ_{j,l} a_kj a_k′l·Cov(B_kj, B_k′l).
  - Each subset's main variance and each input's total come from the same `TermPairMoments` shares.
  - The r × r matrices of subset covariances are free of ℓ, so form them once per draw and weight them by K_ℓk K_ℓk′ at each location.
  - Fitting the components independently does not remove the cross terms.
- *The forms,* as D08's functional responses:
  - *per output location* (R's `sobolBasis`): indices as functions of ℓ;
  - *integrated over the locations,* each equally weighted (the empirical grid), with the location as one more input: a total decomposition.
  - Both report main and total indices, and the law of total variance links them.
- *Units.* Indices are dimensionless, and variance contributions are in output units, with y_s applied. R reports decomposition units, a factor of y_s² too small.
- *Joint draws* come from `JointDraws`, with their weights.
- *Order and budget.*
  - The default order is the largest term over all components, counted in variables.
  - Only subsets inside some term are enumerated.
  - A subset budget refuses before computing; R's `nind` filters after.
- *Zero total variance,* at a location or integrated, is a status with NaN indices, never 0 or 1. The clamping and tolerance follow D08 as stage 6 built them.
- *Measures* are `SensitivitySpec`'s, read through the shared transform.

## Checks

- *A one-component basis fit* gives exactly the direct fit's decomposition of that component, scaled by y_s(ℓ)²·K_ℓ1².
- *Against a dense computation:* the cross terms against the variance of the reconstructed output by tensor quadrature or Monte Carlo under the measure, per location, on a two- or three-component fit.
- *O-O2's units section,* the corrected expectation: the physical variances, against R's decomposition-unit ones; and the indices, against its unit-free reference.
- *The Sobol identities* per location and integrated; the law of total variance between the two forms.
- *Draw pairing and weights, singleton locations and components, and zero variance.*

## Review before reporting

Before you write the report, take a fresh-eyes pass and fix what it finds:
- *Is it correct?* The K weights and y_s at each location, the cross terms' symmetry, and the pairing of draws across components.
- *What is missing?* A D08 basis bullet with no test, a refusal without its reason.

Name the passes and their fixes in the report.

## Report

Under "Reports" below, newest first:
- what landed, with the commits;
- the tests and suite counts;
- any change outside your files;
- the questions you left as reversible defaults, with your reasons;
- the review passes;
- anything the main session must know, including the docs to write.

## Reports

Newest first.
