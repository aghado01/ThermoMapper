# BASS port: response bases, prediction and sensitivity: checklist

Companion to [port-outputs-brief.md](port-outputs-brief.md). An item is ticked when its report lands. Reports go under "Reports", newest first, each with its date, commits, evidence and remaining limitations.

## Checklist

### Decisions (stage 1)

- [x] **O-D07** Response units and uncertainty: decided 2026-10-29, as recommended (the report below).
- [x] **O-D08** Sensitivity measures and moments: decided 2026-10-29, as recommended (the report below).

### Reference fixtures

- [x] **O-O1** R defect fixtures: the truncation sampler's indexing, recycled noise, the residual's sign and units, 100% retention, constant columns, draw order and dropped axes, interval quantiles: `r/oracles/bass_output_defects_oracle.R`.
- [x] **O-O2** Sensitivity defect fixtures: ignored category weights, unclipped truncation bounds, decomposition units, hard-coded effects: `r/oracles/bass_sensitivity_defects_oracle.R`.
- [x] **O-O3** Moments against quadrature and category enumeration, including higher degrees and truncated mixtures: `r/oracles/bass_moments_oracle.R`.
- [x] **O-O4** Truncation residuals: the in-sample residual is zero at full retention; leave-one-case-out equals explicit refits, and the computation in case coordinates equals the one in output space: `r/oracles/bass_truncation_oracle.R`.

### Build (stages 5, 6, 9)

- [x] **O1** PCA and supplied-basis decomposition on the thin SVD (stage 5 chip, `2fec1ee`).
- [x] **O2** Component execution and joint draws across components (`bb065c6`).
- [x] **O3** Prediction: latent and predictive draws, component noise, truncation residuals, physical units (`29e7e23`, `bb065c6`).
- [x] **O4** Pointwise density for supplied-basis fits (`502fb0f`).
- [x] **O5** One moment service for direct and basis sensitivity (stage 6 chip, `17dfe2d`, `2fa7ad6`).
- [x] **O6** Direct sensitivity, including both functional forms (`835a310`).
- [ ] **O7** Basis sensitivity with cross-component terms and a subset budget.
- [x] **O8** Effect functions for the declared measure and degree (`835a310`).
- [ ] **O9** Presentation data above maths (stage 9).

### Room for extensions

- [x] **O-X1** The decomposed surface is a declared input (`DecomposedSurface`, the mean built).
- [x] **O-X2** The moment service is written over piecewise polynomials (`PiecewisePolynomial`).
- [x] **O-X3** The basis contract admits robust bases, an uncertainty slot and a declared truncation model (`FromBasis`, `BasisUncertainty`, `TruncationResidual`).
- [ ] **O-X4** Term derivatives are available to sensitivity.
- [x] **O-X5** The response representation carries a declared reconstruction map (linear now; composition for elastic emulation) (`ReconstructionMap`, `ReconstructionDomain`).

## Reports

### Stage 5 and stage 6's first part verified — 2026-10-29

Two chips built them beside each other on the shared checkout ([stage 5's report](stage5-outputs-brief.md#reports), [stage 6's](stage6-sensitivity-brief.md#reports)). Read and verified here.
- *Commits.*
  - Stage 5: `2fec1ee`, `29e7e23`, `bb065c6`, `502fb0f`, `dd533be`, `ac4730e`, and its report `cf0cc3a`.
  - Stage 6: `17dfe2d`, `2fa7ad6`, `835a310`, and its report `853b61f`.
  - Each touched only its own files. Stage 5 added only `Predict` to `BassFit.cs`. Stage 6's Distributions additions are generic (`IIntervalMoments`, `Uniform`, `Mixture`, `Empirical`, `GaussLegendre`, and interval moments on `Normal`, `StudentT` and `Truncated<T>`).
- *The suites on the combined tree,* with the main session's wiring: Regression 911/911, Distributions 148/148, Samplers 186/186, LinAlg 120/120, and the BASS oracle tests (below).
- *Wired here.* `BassFit.Sensitivity(data, spec)` and `BassFit.Effects(data, spec, effects)` call stage 6's static entries. A views test checks that they match the entries and refuse data of another fingerprint.
- *O-O1's constant-input field* now records D07's refusal: `corrected = "refused"`, which the parity test asserts. The stage-1 checks say the same.
- *Docs:*
  - `docs/Maths/BayesianRegression.md` gains the sections "BASS: decomposition, components and prediction" and "BASS: sensitivity", the API, configuration, status and source-map rows.
  - `docs/Maths/ProbabilityAndInformation.md` gains the interval moments.
- *Left as reversible defaults,* 12 in each chip's report; none re-opens D07 or D08. The main ones:
  - `Predict`, `Sensitivity` and `Effects` take the training data, since the fit keeps neither its knots nor its inputs;
  - a basis fit's intervals are weighted draw quantiles;
  - the normal's interval moments come from the density's Taylor series, not the Φ/φ recurrence, which loses every digit on narrow pieces by degree 2–3;
  - the indices' posterior quantiles wait for stage 9's weighted quantile.
- *Next:* O7, basis sensitivity, on stage 5's `BassBasisFit` and stage 6's `BassMoments.Pair`. O-O2's units section is its parity check.

### D07 and D08 decided — 2026-10-29

The user took the recommendations of a review against the built stage-4 engine. The segment brief marks each change *(2026-10-29)*.
- *D07, the choice:* leave-one-case-out truncation residuals by default, R's in-sample residual the named option `InSample`, their RMS ratio reported as a basis-stability diagnostic.
- *D07, amended:* a constant training column is refused at entry (the contract shapes; S1 refuses it), not mapped to 0; R's W = V_rᵀ convention kept for comparability, since D03's scale-free defaults make a component's fit independent of W's scale.
- *D07, added:*
  - components whose draw counts differ under the R̂ stop are thinned evenly to the shortest before pairing;
  - direct functional prediction at new output coordinates, flagged outside the training range;
  - a Rao–Blackwellized mean surface as an option, never a draw;
  - every operation takes the draw weights.
- *D08, the choices:* R's uniform default kept, the empirical marginal added as an option; Student moments by quadrature only, at every ν and degree, uniform and normal in closed form at every degree.
- *D08, placement:* interval moments in Distributions, hinge products and ANOVA in BASS (the hub's boundary rules).
- *Launch:* stage 5 (O1–O4) and direct sensitivity (O5, O6, O8) as two chips, [stage5-outputs-brief.md](stage5-outputs-brief.md) and [stage6-sensitivity-brief.md](stage6-sensitivity-brief.md); basis sensitivity (O7) after stage 5's basis fits land.

### R's Student moments (O-O3) — 2026-10-27

From the reference-fixtures chip's follow-up. `976ac0c` installs `hypergeo` 1.2-15 (with `contfrac` and `elliptic`) in the oracle environment, and `f7acf2a` runs R BASS's Student moments in `r/oracles/bass_moments_oracle.R`. All checks pass.
- **As written, R's Student moments fail.** They give a value only where a hinge misses the interval, which `C1` and `C2` return as 0 without evaluating anything; every other case errors.
  - `robust2f1` calls `hypergeo::f15.3.8`, the Abramowitz–Stegun 15.3.8 transformation. For |x| < 1 it returns a complex value, which fails `C1`'s and `C2`'s sign checks ("invalid comparison with complex values").
  - For |x| ≥ 1, R passes it the Pfaff-transformed argument in (0, 1), where it returns NA. R's commented-out code used `gsl::hyperg_2F1` there.
  - `hypergeo`'s general `hypergeo()` does not converge at these arguments either.
- **With only the 2F1 evaluator repaired:** R's own branches, Pfaff transform included, evaluated by the Gauss series summed to relative 10⁻¹⁷.
  - At degree 1 R's forms match quadrature at ν = 2 and 3 (first moments 8·10⁻¹⁴, cross moments 3·10⁻¹³).
  - At degrees 2 and 3, `intabq2.student` has no degree check and returns the degree-1 cross moment (relative error 1). `intabq1.student` refuses degrees above 1.
  - At ν = 1 the forms divide by ν − 1 and give NaN, and so does the ν = (1, 3) mixture.
- **For D08:** with a working 2F1, R's closed form is usable at degree 1 for ν ≠ 1. ν = 1 needs D08's quadrature fallback, and degrees above 1 have no R form to port. The quadrature references cover all of them.

### Output, sensitivity and moment fixtures (O-O1, O-O2, O-O3) — 2026-10-27

From the reference-fixtures chip ([brief](reference-fixtures-brief.md)), `a6960f0`. Three oracles, run by `r/scripts/bass_stage1_checks.R`, all checks passing. Each runs the pinned R code (`2671014`) as written, on mocks where a fitted object is needed, beside the corrected expectation.
- **O-O1,** `r/oracles/bass_output_defects_oracle.R`, on 8 curves at 40 locations, with spikes that make the residuals asymmetric.
  - *Truncation sampler.* With 40 locations and 8 cases R stops ("subscript out of bounds"). With 3 locations it draws only cases 1–3 of the 8. The corrected sampler draws whole case vectors over all 8.
  - *Component noise.* R's spread across 5 cases within each of 200 draws is exactly 0 (correlation 1). The correction is independent noise per case and draw.
  - *Residual,* with `scale = TRUE`, r = 2 and y_s from 0.61 to 1.27. R stores exactly minus D07's e, and its training reconstruction misses Y by up to 3.16 (Y's scale is 2.96). D07's y = y_m + y_s∘(Kw + e) reproduces Y to 4·10⁻¹⁶. The residuals' skewness is 3.76 in D07's orientation and −3.76 in R's.
  - *`perc.var = 100`.* R stops ("NA/NaN argument"): its `n.pc` expression gives NA because the cumulative share is exactly 100. The corrected rank is 7 = n − 1.
  - *Constant columns.* Training maps a constant input to 0, while prediction maps it to NaN at the same value and Inf at a new one. A constant output column under `scale = TRUE` stops R (y_s = 0).
  - *Draw order.* Draws requested as (4, 1, 5, 2) over model groups (1, 1, 2, 2, 3) come back as (4, 1, 2, 5). Requests of 4 × 1 and 1 × 1 lose their axes. With `nugget = TRUE` the noise lands on draw 5 (s² = 0) instead of draw 2 (s² = 1).
  - *Intervals.* With means spread over 1000 draws, R's 95% interval is ±3.91 against the predictive mixture's ±2.77. With the sds spread, it is ±5.20 against ±2.64.
- **O-O2,** `r/oracles/bass_sensitivity_defects_oracle.R`.
  - *Category weights.* On both paths R's moments are subset size over level count: 1/3, 2/3 and 1/3 for {a}, {a, b} and their intersection over three levels, whatever the supplied probabilities. Enumeration gives 0, 0.5, 0 under (0, ½, ½) and 0.9, 0.95, 0.9 under (0.9, 0.05, 0.05). In the inventory's two-level case R gives 0.5 and enumeration 0.9. None of `add_tlCat`, `C2Cat`, `C1Basis` and `C2Basis` reads `prior.cat`.
  - *Truncation bounds* under U(0.25, 0.75). The basis path's `C22Basis` gives 0.45 for the inventory's 0.1105903, and is wrong in 5 of the 6 sign and knot cases: for example 0.167 against 0 for (−, −) at (0.2, 0.4). The direct path's `C2` clips to the interval and matches quadrature to 7·10⁻¹⁷.
  - *Units.* `sobolBasis`, on a mock two-component fit with y_s from 2 to 5, returns variances equal to the decomposition-unit reference to 9·10⁻¹⁶, a factor y_s² below the physical ones (up to 96% off). Its indices match the unit-free reference to 1.5·10⁻¹⁵.
  - *Effects.* `C1.temp`, which `sobol_des`'s effects and returned mean use, is the U(0, 1) moment (0.35, 0.30) whatever the measure. `C1`, which the variances use, matches quadrature: 0.289 and 0.204 under U(0.25, 0.75), and 0.336 and 0.249 under the normal mixture. All three `makeBasis` calls in `sobol_des` (the effect curves, lines 405, 436 and 437) pass degree 1, so a degree-2 term's effect curve is the degree-1 hinge.
- **O-O3,** `r/oracles/bass_moments_oracle.R`. Reference values for 156 cases under 9 measures.
  - *Cases:* first moments with knots inside and outside the support, and cross moments with equal and distinct knots and every sign pair, each at degrees 1–3.
  - *Measures:*
    - uniform on [0, 1] and on [0.25, 0.75];
    - a normal mixture under R's truncation convention and under the conditional one;
    - a narrowly truncated normal;
    - Student ν = 1, 2 and 3, and a ν = (1, 3) mixture.
  - *Agreement.* Adaptive quadrature and a composite Gauss–Legendre rule agree to 5·10⁻¹⁵. R's own moments match them where R has them: uniform at every degree to 1.4·10⁻¹³, normal at degree 1 to 8·10⁻¹⁵. R refuses normal degrees above 1, and its Student moments need `hypergeo`, which is not installed.
  - *Categorical:* P(D) and the intersections by enumeration, with a zero-probability level (0 where R gives 0.25).
- **Not covered:**
  - R's Student moments.
  - The functional sensitivity path (`sobol_des_func`).
  - Cross-component moments against Monte Carlo (D08's own check, for O5–O7).
  - `perc.var = 100` when rounding pushes the cumulative share above 100.
  - Constant output columns under `scale = FALSE`, which R handles.

### Truncation-residual fixture (O-O4) — 2026-10-26

`r/oracles/bass_truncation_oracle.R`, run by `r/scripts/bass_stage1_checks.R` (cwd = r/). It emits both residual sets (in-sample, and leave-one-case-out by explicit refits in output space, the reference), each case's leave-one-out rank, the RMS of each, and their ratio (JSON null when the in-sample residual is zero). Checks, all passing:
- The case-coordinate computation, from the one thin SVD, matches the reference for both residuals to about 1e-14 relative, at n = 2, 5, 10 and 40.
- At full retention (n = 5, r = 4), the in-sample residual's largest entry is 1.3e-14, inside the rank tolerance of 8.3e-13; the leave-one-out ranks are capped at 3.
- With two cases, the leave-one-out residual is their difference (to 9e-15).
- The basis-stability ratio is 2.85 at n = 5, 1.41 at n = 10 and 1.08 at n = 40 (99% retention).

Every numerical rank is judged against one tolerance, max(n, m)·ε·d₁ for the whole data. A tolerance relative to each leave-one-out matrix's own largest singular value would read rounding noise as rank 1 when the other cases leave a zero matrix (n = 2). `_common.R`'s `emit_json` now passes options through to `write_json`, so this oracle writes NA as null. Scaling (y_s ≠ 1) is outside this oracle. No C# side yet; parity tests come with O1.

### Truncation residual held out of the basis — 2026-10-26

From reading Andros et al. (arXiv 2406.18751, 2602.09351). D07's draft now proposes leave-one-case-out truncation residuals by default, with R's in-sample residual as an option and their ratio as a basis-stability diagnostic. The scratch check behind it is exploration (one toy family), not a decision; O-O4 is the fixture to build. O-X5's reconstruction map now also declares where it can be evaluated.

### Prep for elastic emulation — 2026-10-26

"Room for extensions" gains the representation with a declared reconstruction map (O-X5), from Francom, Tucker et al.'s elastic calibration (arXiv 2305.08834).

### D07 and D08 drafted — 2026-10-25

Proposals, not decisions: `0844bdc`; categorical moments' source added in `a32c701`.

