# BASS stage 6: the moment service and direct sensitivity (chip brief)

Date: 2026-10-29. The user decided D08 on 2026-10-29 ([the outputs segment](port-outputs-brief.md), its D08 marked *(2026-10-29)* where the decision amended the draft), and this chip builds the first part of stage 6 on it: the outputs checklist's O5, O6 and O8.
- *What exists.* Stage 4 is on `main`: `BassFit` with its interned terms (`BassTermTable`) and per-chain `DrawStore<BassDraw>`, each draw's uncentred coefficients, and `BassInputTransform`. The fixtures `r/oracles/bass_moments_oracle.R` (O-O3: 156 cases under 9 measures) and `bass_sensitivity_defects_oracle.R` (O-O2) are built, with no C# side yet.
- *This chip builds:*
  1. O5, one moment service for direct and basis sensitivity;
  2. O6, direct sensitivity: Sobol main and total indices per posterior draw, for scalar fits and both functional forms;
  3. O8, effect functions for the declared measure and degree.
- *Not in scope.*
  - O7, basis sensitivity, which needs stage 5's basis fits and follows them. The moment service must serve it all the same (below).
  - Prediction (stage 5, a chip beside you).
  - Presentation (stage 9).
  - Input dependence: the measure is a product.

Read first: [AGENTS.md](../../../AGENTS.md) and [CONTRIBUTING.md](../../../CONTRIBUTING.md). Then:
- the outputs segment's D08 in full, and "Room for extensions" (O-X2: the moment service is written over piecewise polynomials, not hinges only; O-X4: term derivatives later);
- the hub's boundary rules: Distributions owns a measure's densities, CDFs, truncation and interval moments ∫ₐᵇ xᵏ dF; BASS owns the hinge products and their ANOVA structure;
- the outputs checklist's reports on O-O2 and O-O3, including "R's Student moments", which say what each fixture checks and where R fails;
- the JSS paper's §3 on Sobol decomposition for BASS, and Francom et al. 2018. The bibliography is in the hub's references.

In the code:
- `BassTerm.cs` (hinge and categorical factors, the degree), `BassBasis.cs` (how a term is evaluated, for the reference checks), `BassFit.cs`, `BassState.cs` (`BassDraw`), `BassInputTransform.cs`;
- `src/maths/distributions/` (`Normal`, `StudentT`, `Truncated`, `SpecialFunctions`);
- for how an oracle fixture is read, `tests/oracle/BassTargetParityTests.cs`.

## Working rules (a shared checkout)

- **Work in `D:\aghado01\ThermoMapper` on `main`.** Worktrees fail on this repo's long paths. Switch no branch, push nothing.
- **Build and test into your own artifacts root:** add `-p:ArtifactsRoot=D:\aghado01\ThermoMapper\artifacts\chip-stage6\` to every `dotnet build` and `dotnet test`. Do not run `scripts/fact-harness.ps1`. Delete `artifacts\chip-stage6\` when you finish.
- **Two chips work beside you.**
  - The stage-5 chip ([its brief](stage5-outputs-brief.md)) owns `bass/output/`, `bass/posterior/`, `BassFit.cs` and its tests.
  - The F-M2 chip owns `tests/regression/BassMixingExit*Tests.cs`.
  - The main session owns `docs/Maths`, the briefs and checklists.
  - Expect their commits between yours. If a build breaks mid-edit in a file not yours, wait and retry.
- **Your files.**
  - New code in `src/maths/regression/spline/bass/sensitivity/`, in the namespace convention the existing `bass/` code uses.
  - Interval-moment additions in `src/maths/distributions/`, generic and tested in `tests/Maths.Distributions.Tests`.
  - Tests `tests/regression/BassSensitivity*Tests.cs` and `tests/oracle/BassSensitivity*ParityTests.cs` and `BassMoment*ParityTests.cs`.
  - This brief's Reports section.
  - Leave `BassFit.cs` alone: expose a static entry (say `BassSensitivity.Of(fit, spec)`), and the main session wires `BassFit.Sensitivity` after you report.
- **Shared code stays generic.** Distributions' additions know nothing of hinges or BASS; name each in your report. Rerun the Distributions suite.
- **Commit your own files only, by explicit path** (`git add <paths>`, never `-A` or `.`), one commit per task. Use the repo's message style, ending with `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`.
- **Line endings:** briefs and C# are LF.
- **Decisions are not yours to settle.** D08 fixes the design. Where it leaves something open (an API name, a quadrature node count, a tolerance), implement a reversible recommendation and put the question in your report with your reasons.
- **No oracle-bug parity.** R's defects are fixtures: assert the corrected value, and document R's beside it.

## What to build

**O5. The moment service.**
- *Measures, one per input, declarative records:*
  - continuous:
    - uniform on an interval (default: the training range, which is [0, 1] in transformed units);
    - a finite mixture of truncated normal and Student t components, with weights before truncation by default (R's convention) and weights conditional on truncation as the other named form. Student components take a location-scale `scale`, with a constructor from a marginal sd for ν > 2;
    - the empirical marginal: the input's training values, equally weighted;
  - categorical: level probabilities, uniform over the observed levels by default.
- *Interval moments* ∫ₐᵇ xᵏ dF, in Distributions:
  - uniform and (truncated) normal in closed form at every k;
  - Student by Gauss–Legendre, split at the bounds, its error reported by doubling the nodes, at every ν, ν = 1 included;
  - the empirical marginal as exact sums.
- *Hinge moments,* in BASS:
  - each term pair's moment E[B_j B_l] and each term's E[B_j] under the product measure: a product over variables of one-dimensional moments of hinge powers, split at the knots and clipped to the measure's support;
  - categorical factors as sums of the supplied probabilities, intersections included, never ratios that can give 0/0;
  - every hinge degree.
- *Written for O7 too:* the service takes any two terms, from one fit or from two component fits' term tables, so basis sensitivity's cross-component moments use it unchanged. Write it over piecewise polynomials where that costs nothing more (O-X2).
- *Cost.* Draws repeat their structures, so cache moments per interned term and term pair across draws.

**O6. Direct sensitivity.**
- *Per posterior draw,* decompose the latent surface's variance under the measure:
  - Sobol main indices for each subset of inputs up to a declared order (default: the fit's largest interaction), and total indices per input;
  - variance contributions in the response's units.
- *The posterior distribution of each index,* from the draws, taking their weights (equal for MCMC).
- *The decomposed surface is a declared input* (O-X1), the mean by default.
- *Functional responses, both forms:*
  - integrate the functional variable as one more input, with its measure declared (uniform over its range by default, or the empirical output grid);
  - condition on a functional coordinate, giving indices as functions of t at declared coordinates (default: the training grid).
  - Both report main and total indices.
- *Zero total variance:* the indices are undefined and reported with a status, never 0 or 1. A small negative component variance from rounding is clamped within a tolerance relative to the total, and reported beyond it.
- *A subset budget* in place of R's filter after the fact.

**O8. Effect functions.** Main and two-way effects for the declared measure and degree, or an explicit refusal. R hard-codes the uniform measure and degree 1 (the O-O2 defect).

## Checks

- *Moments against O-O3's quadrature references,* all 156 cases under the 9 measures: degrees 1–3, knots inside and outside the support, equal and distinct knots, every sign pair, truncated mixtures under both conventions, and Student ν = 1, 2, 3 and a mixture.
- *Categorical moments* against enumeration, including a zero-probability level.
- *O-O2's corrected expectations:* category weights read, truncation bounds clipped, effects under the declared measure, and physical units. The units check needs basis fits and waits for O7; note it in your report.
- *The Sobol identities* on a draw: indices over all orders sum to 1, and totals are at least the mains.
- *Against Monte Carlo:* each index, on a fitted draw, against a direct Monte Carlo estimate of the variance decomposition under the measure, within its error.
- *Draw order, weights and zero variance.*

## Review before reporting

Before you write the report, take a fresh-eyes pass and fix what it finds:
- *Is it correct?* Each measure's support and truncation, the clipping at the knots, the categorical intersections, and the subset enumeration's completeness.
- *What is missing?* A degree, measure or ν untested; a refusal without its reason; a D08 bullet with no test.

Name the passes and their fixes in the report.

## Report

Under "Reports" below, newest first:
- what landed, with the commits;
- the tests and suite counts;
- each shared-code change;
- the questions you left as reversible defaults, with your reasons;
- the review passes;
- anything the main session must know: the docs it needs to write, the wiring of `BassFit.Sensitivity`, and what O7 needs from stage 5.

## Reports

Newest first.
