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

### O7 built: basis sensitivity — 2026-10-29

**What landed:** `f1d53ab`, in namespace `Maths.Regression.Spline.Bass.Sensitivity`.
- `BassBasisSensitivity.Of(fit, spec)`, and `BassBasisFit.Sensitivity(spec)` forwarding to it, return a `BasisSobolAnalysis`. No data argument is needed: the components share the fit's predictors.
- *The state.* For each joint draw, the analysis keeps the components' covariances as r × r matrices (packed over k ≤ k′) in the coefficient responses' units, free of ℓ:
  - on every reported subset, Cov(wₖ,u, wₖ′,u) = Σⱼₗ aₖⱼaₖ′ₗ·Cov_u(Bₖⱼ, Bₖ′ₗ);
  - each input's total;
  - the whole covariance;
  - each component's mean E[wₖ].

  They come from `BassMoments.Pair` across components as within one, on one service over component 0's inputs. They are exposed as `ComponentCovariance`, `ComponentMainCovariance`, `ComponentTotalCovariance` and `ComponentMeans`.
- *Two views,* each a `SobolDecomposition` in output units, with g_ℓk = y_s(ℓ)·K_ℓk:
  - `AtLocation(ℓ)`, R's form, computed on demand: the variance is Σₖₖ′ gₖgₖ′·Cov(wₖ, wₖ′), and the mains and totals are weighted alike. `PerLocation()` returns every location.
  - `Integrated`: the location is one more input, `"location"`, after the predictors, with the training locations equally weighted. With ḡ the locations' mean of g and Σ_G their covariance:
    - a subset u of the inputs has main variance ḡᵀC_uḡ;
    - u with the location has tr(Σ_G·C_u);
    - the location alone has the variance over ℓ of the location means y_m(ℓ) + g_ℓᵀE[w];
    - an input's total is tr(S_G·T_v), the locations' mean of its totals;
    - the location's total is its main plus tr(Σ_G·C).
- *Draws:* `JointDraws` pairs them. MCMC joint draws carry their chain and index within it. Resampled ones carry chain 0 and their joint index. The joint draws weigh equally, so `Weights` is null.
- *Order, budget and status* follow stage 6:
  - The order counts variables, the location counted as one in the integrated form. Undeclared, the default is the largest term over all components' joint draws per location, and one more integrated.
  - `SubsetBudget` bounds the subsets of the inputs (the kept matrices) and refuses before any pair moment is computed.
  - `BassSensitivity.Status` clamps negatives within the tolerance and flags those beyond it, and marks zero variance with NaN indices.
- *Refused, with reasons:*
  - `Coordinates` (the basis has none);
  - `OutputMeasures` (the locations are integrated equally weighted);
  - components fitted at different hinge degrees;
  - the spec's own refusals.

**Tests.**
- New:
  - Maths.Regression.Tests: `BassSensitivityBasisTests` (9).
  - Maths.Oracle.Tests: `BassSensitivityDefectParityTests.DecompositionUnits_AreThePhysicalOnes`.
- Suites:
  - Maths.Regression.Tests 920/920, built from the tree with this chip's files in place. The review's two added assertions came after that run; with them, all 31 sensitivity tests (22 direct, 9 basis) were rerun and pass.
  - The BASS oracle suites, Maths.Oracle.Tests filtered to `Bass`: 54/54, with R present.
- Evidence:
  - *One component:* the kept state is the direct decomposition of that component to the bit, apart from the values the direct one clamps. Each location's variance, mains and totals are y_s(ℓ)²·K_ℓ1² times the direct ones, to 10⁻¹³, with the same statuses and indices. The fixture is rank 1 of the stage-5 PCA fit, whose component-0 fit carries over by fingerprint.
  - *Dense check:* on the stage-5 PCA fit, under a truncated normal mixture and a uniform on [0.1, 0.9], tensor quadrature of y_m + y_s·Kw on three joint draws agrees to 10⁻¹⁰ of the variance with:
    - every location's variance, two mains, interaction and totals;
    - the integrated variance, the mains of z0, z1 and the location, z0 and z1 each with the location, and all three totals.

    The cross terms carry more than 10⁻³ of some location's variance.
  - *Cross terms, three checks:*
    - off-diagonal covariances are present;
    - they satisfy Cauchy–Schwarz;
    - the subsets' matrices sum to the whole.

    Reversing the components and K's columns leaves every location's and the integrated decomposition unchanged.
  - *Identities, both forms, every draw:* the mains over all orders sum to the total, each total is the sum of the subsets holding its input, totals are at least the mains, and the posterior mean indices sum to 1.
  - *Law of total variance:* the integrated variance equals the locations' mean variance plus the variance of the location means. Each input's total is the locations' mean of its totals. Each subset's per-location mean splits into its integrated main and its main with the location.
  - *Pairing:* the components were thinned (one component kept half of each chain) and resampled (all weighted). Each joint draw's diagonal and mean equal the direct decomposition and f₀ of `JointDraws.Draw(k, d)`.
  - *Edge cases:*
    - A singleton location: the two forms agree, and the location's main and total are 0.
    - A zero basis row: that location is `ZeroVariance` with NaN indices and NaN means, never 0.
    - Intercept-only draws: every location is zero, and the integrated form is all location.
    - Constant offsets with zero intercepts: zero integrated.
    - A 10⁻³ variation on offsets of 10⁶ stays defined.
  - *Mutations:* the tests were run against three deliberately broken versions of the code, which was restored after each:
    - dropping the cross-component pairs fails the quadrature and cross-term tests;
    - dropping y_s fails four tests, among them the O-O2 parity test, which then returns R's decomposition-unit value, 0.0886 against 0.354;
    - pairing by raw index fails the pairing test.
  - *O-O2's units section* (corrected expectation): at all 12 locations, to 10⁻¹³:
    - the port's variances equal the physical reference, y_s² times R's `Var.tot`;
    - its mains and interaction equal y_s² times the reference, where R's `S.var` are the decomposition-unit values;
    - its indices equal the unit-free reference, as R's do.

    R's largest relative miss on the physical variance is 0.96.

**Changes outside my files:**
- `BassSensitivity.cs` (behaviour-free): `Subsets` takes the terms' variable lists, not the direct path's private `Setup`, and `SubsetComparer` is internal, so basis sensitivity shares one enumeration. The direct tests pass unchanged in the 920 and the 31.
- `r/oracles/bass_sensitivity_defects_oracle.R`: the units section also emits `reference_effects_decomposition`, `r_effects` and `r_indices`, arrays it already computed. The JSON held only maxima and per-location totals, so the effects and indices could not be checked per location. Existing fields are unchanged.

**Questions left as reversible defaults:**
1. *One analysis, two views,* not the direct path's `SobolAnalysis` with `Functional` choosing a form.
   - The kept r × r matrices are each draw's ℓ-free state. Both forms come from one computation, and per-location arrays (m × draws × subsets) are built only when asked for.
   - The same state could decompose any linear functional of the output, such as a regional average, later.
   - `Functional` is therefore ignored.
2. *y_m is in the integrated form,* as the location means' variation: the brief's surface is y_m + y_s·Kw. Where the mean curve varies strongly, the location's main index dominates the integrated form. The per-location form and the inputs' totals do not see y_m. Decomposing y_s·Kw alone is the other choice.
3. *The location means' zero test.* That part is a sum of squared deviations. Its rounding scales with ε²M², not with the εM² of a difference of moments. It is therefore zero when its standard deviation is at most the tolerance times the means' magnitude M, which is |y_m(ℓ)| + Σₖ|g_ℓk|·(|aₖ₀| + Σⱼ|aₖⱼE Bₖⱼ|). The direct rule, the tolerance times M², would zero a 10⁻³ variation on a 10⁶ offset; a test holds this.
4. *Equal location weights,* as the brief says. Weights for an uneven grid would be one more declared input; `OutputMeasures` is refused meanwhile.
5. *Mixed hinge degrees are refused.* `BassMoments` holds one degree. A mixed-degree pair needs hinge products with a degree per hinge, a small change to `PiecewisePolynomial.HingeProduct` and `BassMoments` if it is ever wanted.
6. *The order counts the location* in the integrated form, as the direct form counts t. With a declared q, u with the location is reported only for |u| ≤ q − 1. The totals do not depend on the order.
7. *The budget counts the subsets of the inputs.* The integrated form adds at most 2U + 1 subsets derived from them, at no enumeration.
8. *Memory:* the kept state is draws × subsets × r(r+1)/2 doubles, for example 44 MB at 1000 draws, 100 subsets and r = 10. Nothing bounds it beyond `SubsetBudget`.
9. *Names:* the location input is `BasisSobolAnalysis.LocationInput` = `"location"`. A predictor of that name would share it, though its index stays distinct.

**Review passes** (both before this report):
- *Correctness:*
  - K and y_s at each location: g = y_s(ℓ)·K_ℓk in both views.
  - The cross terms' symmetry: every ordered pair across components, each unordered pair once within one, doubled off the diagonal in the weights.
  - The pairing of draws: through `JointDraws.Draw(k, d)` for every component.

  It found no defect. The mutation runs above confirmed that the tests catch each of the three.
- *Coverage:* two gaps, each closed with a test.
  - The location means' zero rule was a decision with no test: now the 10⁶ offset case.
  - "`Functional` is ignored" was untested: now the condition form's integrated variance equals the default's.

  Every D08 basis bullet now has a test:
  - cross-component covariances;
  - output units against R;
  - the budget;
  - zero variance at a location and integrated;
  - both functional forms with mains and totals;
  - declared measures;
  - weights and pairing.

  Each refusal carries its reason: coordinates, output measures, degree, tolerance, order, unknown name, budget, and a location out of range.

**For the main session:**
- *Docs to write (docs/Maths):*
  - the basis algebra: C_u per draw, free of ℓ, weighted by gₖgₖ′;
  - the integrated form's formulas (ḡᵀC_uḡ, tr(Σ_G C_u), the location means, totals as the locations' mean) and the law of total variance between the forms;
  - why independent fits keep cross terms;
  - y_m in the location's main effect;
  - the location means' zero rule;
  - R's decomposition units as a corrected defect.
- *Effect functions for basis fits come nearly free.* An effect is linear in the surface, so it has no cross terms:
  - at ℓ, the effect on u is y_s(ℓ)·Σₖ K_ℓk times component k's effect;
  - integrated, ḡ weights them, and (g_ℓ − ḡ) gives u with the location.

  The direct path's per-term centred factor tables give each component's effect per draw. The work is reindexing them by joint draw and a cell budget that scales with m: a refactor of `BassSensitivity.Effects`' term loop to take a term table and coefficients, not a fit.
- *The checklist:* O7 can be ticked (`f1d53ab`). O-O2's units section is now asserted.
