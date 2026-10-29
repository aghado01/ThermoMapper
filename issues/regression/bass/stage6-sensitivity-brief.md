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

### Stage 6 first part built: O5, O6, O8 — 2026-10-29

**What landed** (new BASS code in `bass/sensitivity/`, namespace `Maths.Regression.Spline.Bass.Sensitivity`):
- `17dfe2d` Distributions' interval moments (shared code, below).
- `2fa7ad6` O5, the moment service.
  - `ProductMeasure.Resolve` carries the declared measures (`UniformMeasure`, `MixtureMeasure` of `NormalComponent` and `StudentComponent` under `MixtureWeighting.BeforeTruncation` or `ConditionalOnTruncation`, `EmpiricalMeasure`, `CategoricalMeasure`) into the fit's units, each a `ContinuousMeasure` (its interval moments and support) or level probabilities.
  - `PiecewisePolynomial` is O-X2's form: breakpoints, and pieces that expand about any centre. Its `Expectation` integrates piece by piece against the measure's interval moments, about the centre of the piece clipped to the support. `HingeProduct` builds hinge products, and with offsets their centred products, so a covariance is integrated directly, never as E[φφ] − E[φ]E[φ].
  - `BassMoments.On(inputs, degree, measure)` gives `Mean(term)`, `Cross(a, b)` and `Pair(a, b)`. `Pair` returns a `TermPairMoments`: the covariance's share on every nonempty subset of the shared variables (`SubsetCovariance`) and each shared variable's total share (`VariableCovariance`). Terms are taken by value, so two component fits on the same inputs meet in one service. Factor, factor-pair, term and term-pair results are cached.
- `835a310` O6 and O8, in one commit because `BassSensitivity.cs` holds both entries.
  - `BassSensitivity.Of(fit, data, spec)` returns a `SobolAnalysis`: one `SobolDecomposition` for a scalar fit or the integrated form, one per coordinate for the conditioned form. A decomposition holds per draw each subset's main variance, each input's total and the total variance, in response units, with a `DecompositionStatus`. Indices and their weighted posterior means are views on it.
  - `BassSensitivity.Effects(fit, data, spec, effects)` returns an `EffectAnalysis` of main and two-way effects (`Effect`, `EffectGrid`) and each draw's mean f₀.
  - The algebra: a term is κ·Πᵥφᵥ, so a pair's components covary only on subsets u of their shared variables, with Cov_u = κκ′·Π_{u}Cov(φ, φ′)·Π_{rest}E φ. The variance, the main variances and the totals are sums of aⱼaₗ times the cached pair shares. Draws with repeated structures cost one pass over their pairs.

**Tests.**
- New:
  - Maths.Distributions.Tests: `IntervalMomentsTests` (19).
  - Maths.Regression.Tests: `BassSensitivityMomentTests` (10), `BassSensitivityTests` (9), `BassSensitivityEffectTests` (3).
  - Maths.Oracle.Tests: `BassMomentParityTests` (2), `BassSensitivityDefectParityTests` (4).
- Suites:
  - Maths.Distributions.Tests 148/148.
  - Maths.Regression.Tests 910/910, built from the tree at `cf0cc3a` with this chip's files in place. The review's added assertions came after that run and were rerun on their own: 22/22 sensitivity tests and 6/6 parity tests.
  - The six new parity tests pass, with R present.
- Evidence:
  - *O-O3:* all 156 cases under the 9 measures match the quadrature references to within 6·10⁻¹⁵ relative. The largest gap is 5.7·10⁻¹⁵, under Student ν = 2. The Student Gauss–Legendre estimates are at most 7·10⁻¹⁵, and the closed forms report 0. R's repaired Student forms miss by more than half in 104 cross-moment cases at degrees 2–3, because they return the degree-1 value.
  - *Categorical:* the moments equal enumeration exactly, with a zero-probability level giving 0 where R gives 0.25. On O-O2's cases they give 0, 0.5, 0 and 0.9, 0.95, 0.9 where R gives ⅓, ⅔, ⅓, and the two-level case gives 0.9 where R gives 0.5.
  - *O-O2, truncation:* every case under U(0.25, 0.75) matches quadrature to 10⁻¹⁴. R's basis path is wrong in 5 of the 6.
  - *O-O2, effects:* the mean's term moments match quadrature under each measure. R's C1.temp agrees only under U(0, 1). A degree-2 term's main effect is the degree-2 hinge less its mean, to 10⁻¹⁴, and sits more than 0.1 from the degree-1 curve that R draws.
  - *Term moments at degrees 1–3* match a tensor quadrature of `BassBasis.CaseColumn` to 10⁻¹². The measure is a Student mixture conditional on a truncation inside the training range, a uniform on an inner interval, and probabilities with a zero level.
  - *The empirical marginal* equals the exact average of the columns over the training cases to 10⁻¹⁴, at every degree. Every atom there sits on a knot, which checks the half-open pieces.
  - *Sobol identities* hold on every draw of a fitted model to 10⁻¹²: the mains over all orders sum to the total, each total is the sum of the subsets holding its input, and totals are at least the mains.
  - *Monte Carlo:* on a fitted draw, under a truncated normal mixture, a uniform on [0.1, 0.9] and level probabilities, the variance, each single main, each total and each pair's closed variance fall within 5 standard errors of pick-freeze estimates at N = 200 000. The single mains use Saltelli 2010, the totals Jansen.
  - *Functional forms:* the integrated variance under the empirical output grid equals the conditioned form's mean variance over that grid plus the variance of its means, to 10⁻¹¹ on every draw (the law of total variance). The conditioned variance matches a tensor quadrature of the case and output columns at four coordinates, one of them outside the training range and flagged.
  - *Effects* at degree 2 under a normal–Student mixture match E[f | x_v] − f₀ and E[f | x_v, x_w] − f_v − f_w − f₀ by quadrature to 10⁻¹¹.
  - *Draw order, weights and zero variance:* draws come in chain order. Weighted means are taken over the defined draws. An intercept-only draw, and a draw whose term is zero on the measure's support, are `ZeroVariance` with NaN indices, and they are left out of the means rather than counted as 0. Clamping and flagging are checked on `Status` directly.
- *Not checked here:* O-O2's physical-units section needs basis fits and waits for O7.

**Shared-code changes** (Distributions, generic, `17dfe2d`):
- `IIntervalMoments` (new): e^(−logScale)·∫ (x − c)ᵏ dF over (lower, upper], with an error estimate. The log scale lets a truncation keep its digits far in a tail; the half-open interval counts an atom on a breakpoint once.
- `GaussLegendre` (new): cached n-point rules.
- `Uniform` (new): closed form.
- `Mixture` (new): any components with normalized weights. `Truncated.Of(mixture)` is R's convention, and a mixture of truncated components is the conditional one.
- `Empirical` (new): exact compensated sums over the atoms.
- `Normal`: analytic interval moments at every power, by the density's Taylor series on stretches with h·(1 + |z|) ≤ 1, summed to convergence.
- `StudentT`: Gauss–Legendre at every ν on panels graded geometrically from the location, 20 nodes against 40 for the error. A finite interval is required, so it is truncated first.
- `Truncated<T>`: clips the interval and passes its log mass on as the scale.

**Questions left as reversible defaults:**
1. *The normal's "closed form"* is the density's Taylor series, analytic and with no quadrature error. The Φ/φ recurrence was tried in analysis and rejected: on a narrow piece, such as knots 0.0015 apart under σ = 0.2, it loses a factor of about 1/h² per power, which means every digit by degree 2–3. The series is checked against quadrature on wide, narrow, tail and 40–50 sd intervals. If D08's "closed form" meant the recurrence, this is the place to say so.
2. *Measures are declared in each input's own units* (R's convention: `sobol` scales `trunc`, `mean` and `sd` by the training range), defaulting to the training range. Supports beyond it are allowed and flagged (`ContinuousMeasure.Extrapolates`, and `CoordinatesOutside` for coordinates). R warns for inputs and stops for functional variables.
3. *Entry points* take the training data, as `PointwiseLogDensity` and `Predict` do. The fit keeps neither knots nor the training values the empirical measure needs.
4. *Default order:* the posterior's largest term, counted in variables, so every order is reported and the mains sum to 1. The config's `MaxInteraction` would leave out categorical and functional factors.
5. *Subset budget:* 10 000. Exceeding it refuses the call before anything is computed, with the count and what to change. Only subsets inside some posterior term are enumerated, since every other subset's main is exactly 0.
6. *Tolerance:* one `RoundingTolerance`, 10⁻¹⁰. The total variance is zero when it is at most the tolerance times Σ|aⱼaₗ|(|E[BⱼBₗ]| + |E[Bⱼ]E[Bₗ]|). A negative contribution is clamped above −tolerance × total and flagged below it.
7. *Posterior summaries* are weighted means over the defined draws only. Quantiles are left to stage 9, since the repo has no weighted quantile.
8. *Student quadrature:* 20 nodes per panel, the panels graded so that each stays three half-widths or more from the poles at location ± i·scale·√ν. The 40-node values are returned, and their difference from the 20-node ones is the estimate. Each analysis reports the largest estimate as `MomentError`.
9. *Effects:* grids evenly spaced over the measure's support, `GridPoints` 100 as in R; categorical inputs by level. The default set is every main plus the two-way effect of every pair a posterior term holds (R's `combs`). A cell budget of 10⁷ refuses before computing. All three forms are supported, the conditioned form per coordinate.
10. *An empirical output grid over several output variables* is the product of each variable's marginal over the training locations, not the joint grid, since the measure stays a product.
11. *Categorical probabilities* are given by level name. Levels left out get 0, the values are normalized, and a name outside the training levels is refused.
12. *The normal's far tail:* stretches below e^−92 of the interval's peak density are dropped. An interval beyond about 10⁸ sd is refused, because its stretches fall below one ulp.

**Review passes** (both before this report):
- *Correctness:*
  - Each measure's support and truncation: uniform ends, both mixture weightings, the empirical [min, max].
  - Clipping at the knots: a zero-width clipped piece is kept for a possible atom, and a piece beyond the support is skipped.
  - Categorical intersections, and the completeness of the subset enumeration: every pair component's subset lies inside some term, and a pair yields 2ⁿ − 1 components.
  - It found no defect. The one test that failed in development was an expectation of mine: a normal 40–50 sd outside the interval still has a truncated law, because its log mass is finite. The test now asserts that law.
- *Coverage:* three gaps, each closed with a test.
  - Pair mains were not checked against Monte Carlo; now each pair's closed variance is.
  - The conditioned form's identities were unchecked; the identities now run on every coordinate.
  - Nothing showed that the Student error estimate is reported; it is now asserted nonzero for Student measures and 0 for the closed forms.

**For the main session:**
- *Wiring:* `BassFit.Sensitivity(BassData data, SensitivitySpec? spec = null) => BassSensitivity.Of(this, data, spec)`, and `Effects` likewise. Both refuse data of another fingerprint.
- *Docs to write (docs/Maths):*
  - the pair ANOVA formula and why totals and mains are sums over pairs;
  - the normal's series and why not the recurrence;
  - the Student panels and their error;
  - the half-open piece convention for atoms;
  - the zero-variance and clamping rules.
- *What O7 needs from stage 5:*
  - the components' term tables and per-draw coefficients, paired (`BassJointDraws`), on one `BassInputs`;
  - K and y_s.
  - Then the cross-component covariance of the reconstruction is Σ_{c,c′} K_c K_c′ Σ_{j,l} a_cj a_c′l `BassMoments.Pair(term_cj, term_c′l)`, and its subset and total shares come from the same `TermPairMoments`.
  - Variances are scaled by y_s² into output units.
  - O-O2's units section is then the parity check.
