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
