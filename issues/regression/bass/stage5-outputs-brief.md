# BASS stage 5: decomposition, components and prediction (chip brief)

Date: 2026-10-29. The user decided D07 on 2026-10-29 ([the outputs segment](port-outputs-brief.md), its D07 marked *(2026-10-29)* where the decision amended the draft), and this chip builds stage 5 on it: the outputs checklist's O1–O4.
- *What exists.* Stage 4 is on `main`: `Bass.Fit`, `BassFit` with its interned terms and per-chain `DrawStore<BassDraw>`, `BassInputTransform` (S1), and `BassPointwise` for direct fits. The fixtures `r/oracles/bass_output_defects_oracle.R` (O-O1) and `bass_truncation_oracle.R` (O-O4) are built, with no C# side yet.
- *This chip builds:*
  1. O1, PCA and supplied-basis decomposition on the thin SVD, with D07's truncation residuals;
  2. O2, component execution and joint draws across components;
  3. O3, prediction: direct fits (scalar and direct functional, new output coordinates included) and basis fits, latent and predictive draws;
  4. O4, the pointwise density for supplied-basis fits.
- *Not in scope.*
  - Sensitivity (stage 6), which a chip builds beside you.
  - Calibration (stage 7).
  - SMC (stage 8): take weights everywhere, but nothing produces them yet.
  - Presentation (stage 9).
  - Solving coefficients from a supplied basis: prepared coefficient responses are the input (the hub's scope).

Read first: [AGENTS.md](../../../AGENTS.md) and [CONTRIBUTING.md](../../../CONTRIBUTING.md). Then:
- the outputs segment in full: D07 above all, the response-basis conventions, and "Room for extensions" (O-X1, O-X3, O-X5 shape your contracts);
- the hub's scope table (P02, P03, P07, P08) and its placement and boundary rules;
- the fitting segment's "Posterior ownership and prediction", its contract shapes (`BassData`, `BassInputTransform`, `BassFit`, `BassDraw`), and "Per-observation log density per draw (OQ17)";
- the outputs checklist's reports on O-O1 and O-O4, which say what each fixture checks and how.

In the code:
- `BassFit.cs`, `BassState.cs` (`BassDraw`), `BassInputTransform.cs` (`RangeScale.Contains`, `TrainingLevels.IndexOf`, the unseen-level policy), `BassBasis.cs` (evaluating a term's case and output columns: `CaseColumn`, `OutputColumn`), `BassPointwise.cs` and `BassDrawFit`;
- `Bass.cs` (the run, `BassRunOptions.MaxDegreeOfParallelism`) and the seed tree it uses;
- `src/maths/linalg/ThinSvd.cs`;
- for how an oracle fixture is read, `tests/oracle/BassTargetParityTests.cs`.

## Working rules (a shared checkout)

- **Work in `D:\aghado01\ThermoMapper` on `main`.** Worktrees fail on this repo's long paths. Switch no branch, push nothing.
- **Build and test into your own artifacts root:** add `-p:ArtifactsRoot=D:\aghado01\ThermoMapper\artifacts\chip-stage5\` to every `dotnet build` and `dotnet test`. Do not run `scripts/fact-harness.ps1`. Delete `artifacts\chip-stage5\` when you finish.
- **Two chips work beside you.**
  - The stage-6 chip ([its brief](stage6-sensitivity-brief.md)) owns `src/maths/regression/spline/bass/sensitivity/`, interval moments in `src/maths/distributions/`, and their tests.
  - The F-M2 chip owns `tests/regression/BassMixingExit*Tests.cs`.
  - The main session owns `docs/Maths`, the briefs and checklists.
  - Expect their commits between yours. If a build breaks mid-edit in a file not yours, wait and retry.
- **Your files.**
  - New code in `src/maths/regression/spline/bass/output/` (decomposition, components, reconstruction, truncation residuals) and `src/maths/regression/spline/bass/posterior/` (prediction and draw views), in the namespace convention the existing `bass/` code uses.
  - `BassFit.cs`, for `Predict` only.
  - Tests `tests/regression/BassOutput*Tests.cs` and `BassPredict*Tests.cs`, and `tests/oracle/BassOutput*ParityTests.cs`.
  - This brief's Reports section.
  - Leave other production files alone. If you need a change in one, make the smallest behaviour-free change, name it in your report, and rerun the suites that cover it.
- **`BassBasis` is taken** by the term basis. Name the output basis otherwise, for example `BassOutputBasis` or `BassDecomposition`.
- **Shared code stays generic.** A change to LinAlg (`ThinSvd`) is tested in its own host and named in your report.
- **Commit your own files only, by explicit path** (`git add <paths>`, never `-A` or `.`), one commit per task. Use the repo's message style, ending with `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`.
- **Line endings:** briefs and C# are LF.
- **Decisions are not yours to settle.** D07 fixes the design. Where it leaves something open (an API name, a default the decision does not state), implement a reversible recommendation and put the question in your report with your reasons. Do not re-open decided choices.
- **No oracle-bug parity.** R's defects are fixtures: assert the corrected value, and document R's beside it. Never port R's behaviour where D07 corrects it.

## What to build

**O1. Decomposition** (`bassPCA`, `bassBasis`).
- *PCA:*
  - centre by y_m (default: the column means) and scale by y_s (default: 1; the column sds an option);
  - a thin SVD of the standardized responses, Y_cᵀ = U D Vᵀ, with K = U_r D_r and W = V_rᵀ, as D07 states;
  - the rank: the smallest r whose cumulative share reaches the target (≥, not R's strict comparison), capped at the numerical rank. Judge every rank against one tolerance, max(n, m)·ε·d₁ for the whole data, as O-O4's report explains. Full retention means the numerical rank and zero in-sample truncation;
  - a constant output column is centred and gets y_s = 1;
  - a sign convention for comparisons only, with reconstructions or subspaces compared where singular values repeat.
- *Supplied basis:* K (m × r, nonorthogonal allowed) with prepared coefficient responses W (r × n), offsets and scales.
- *Truncation residuals,* D07's decided default.
  - Leave-one-case-out, computed in case coordinates from the one thin SVD (O(n⁴), independent of m), with the rank capped at the other cases' numerical rank.
  - R's in-sample residual is the named option `InSample`.
  - Report the RMS ratio of the two as a basis-stability diagnostic.
  - The residual is observed minus reconstructed, in standardized units.
- *The representation contract* (O-X5): a representation with its reconstruction map, the linear map y = y_m + y_s ∘ (K w + e) as the only member. The map declares where it can be evaluated (a PCA basis only at the training locations). Keep slots for O-X3: a basis-uncertainty slot and the truncation model as a declared choice.

**O2. Component execution.**
- *The fits.* Fit each coefficient response as a `ScalarResponse` on the shared predictors, through `Bass.Fit`, with one config for every component by default.
  - Seeds come from the master seed by the seed tree, one branch per component.
  - Component, chain and replica parallelism share one bound (`BassRunOptions`).
- *Joint draws across components.*
  - Pair equal-weight draws by index.
  - When the components' counts differ (the R̂ stop), thin each component evenly to the shortest before pairing (D07, *(2026-10-29)*).
  - Weighted draws are resampled within each component before pairing. Nothing produces weights yet, but the path is coded and tested with supplied weights.

**O3. Prediction.**
- *Direct fits* (`BassFit.Predict`):
  - latent draws f(x), or f(x, t) for a direct functional fit, from each draw's composition coefficients;
  - predictive draws, which add noise drawn independently for each case, location and draw;
  - draw × case × output axes always explicit, singletons included, in the requested draw order;
  - inputs through the fit's transform, with extrapolation flagged per case and the unseen-level policy declared on the call;
  - new output coordinates for direct functional fits, flagged per location outside the training range;
  - the Rao–Blackwellized mean surface (E[f | S, τ, y] per draw, averaged) as an option, never a draw.
- *Basis fits:*
  - y = y_m + y_s ∘ (K(w(x) + ε) + e), with ε_k ~ N(0, σ_k²) drawn for each prediction case and draw, and e a whole training residual vector, drawn uniformly over the training cases independently for each prediction case and draw;
  - the latent draw omits ε and e;
  - covariance in physical units, diag(y_s)·Cov·diag(y_s).
- *Both:*
  - intervals are quantiles of the predictive mixture from draws;
  - every operation takes the draw weights;
  - a prediction call takes its own seed, and the same seed reproduces the same draws.

**O4. The pointwise density for supplied-basis fits** (OQ17).
- The unit is a case, and the content is Σ_k log N(w_ik | w_k(xᵢ), σ_k²) in coefficient space, on the same `PointwiseLogDensity` seam as direct fits.
- The diagonal Gaussian truncation model, for comparing across bases, is offered as a declared approximation, not the default.

## Checks

- *D07's independent checks:*
  - the training reconstruction y = y_m + y_s ∘ (K W + e) holds exactly;
  - the edge cases: 100% retention, full rank, constant columns;
  - the residual's sign, on an asymmetric fixture;
  - independence of noise across prediction cases;
  - draw order and singleton axes;
  - the covariance transform against a dense computation.
- *Parity against O-O1 and O-O4,* the corrected expectations: the truncation sampler over all training cases, per-case component noise, the residual's orientation and units, 100% retention, constant columns, draw order and axes, interval quantiles; and the leave-one-case-out residual, which O-O4 builds from explicit refits.
- *Prediction at the training inputs* reproduces the stage-4 fitted values that `BassDrawFit` computes, bit for bit or to rounding, for scalar and direct functional fits.
- *Unequal component lengths:* the pairing is checked on fits stopped by R̂.

## Review before reporting

Before you write the report, take a fresh-eyes pass and fix what it finds:
- *Is it correct?* Units and orientation at every step, the noise's independence, the residual drawn whole, the pairing and its thinning, and the axes.
- *What is missing?* An edge case untested (n = 2, r = n − 1, one output location, a singleton prediction), a refusal without its reason, or a D07 bullet with no test.

Name the passes and their fixes in the report.

## Report

Under "Reports" below, newest first:
- what landed, with the commits;
- the tests and suite counts;
- any change to a file not yours;
- the questions you left as reversible defaults, with your reasons;
- the review passes;
- anything the main session must know: the docs it needs to write, and what O7 (basis sensitivity) will need from your contracts.

## Reports

Newest first.

### Stage 5 built: O1–O4 — 2026-10-29

**What landed** (new code in `bass/output/` and `bass/posterior/`, namespace `Maths.Regression.Spline.Bass`):
- `2fec1ee` O1. `BassDecomposition.FromPca` and `FromBasis`, with `BassPcaSpec`, `BassPcaAccount`, `TruncationResiduals`, the O-X5 map (`ReconstructionMap`, `LinearReconstruction`, `ReconstructionDomain`) and the O-X3 slot (`BasisUncertainty`, always null).
  - The rank is the smallest reaching the share (≥), capped at the numerical rank at max(n, m)·ε·d₁; a share of 1 is the numerical rank itself, whatever rounding does to the last shares.
  - Leave-one-case-out residuals in case coordinates, under the mean centre (n/(n−1)·cᵢ and cⱼ + cᵢ/(n−1)) and under a fixed centre (cᵢ and cⱼ).
- `29e7e23` O3, direct prediction.
  - `BassFit.Predict(data, cases, options)` returns `BassPrediction`, with `BassPredictOptions`, the inputs `PredictionCases`, `ContinuousColumn`, `CategoricalColumn` and `UnseenLevel`, and `PredictiveMixture` (the Gaussian mixture's exact quantile, and the weighted empirical one).
  - Internal: `BassSurface` evaluates draws at any inputs with `BassDrawFit.Evaluate`'s accumulation, and `BassConditionalMeans` gives the Rao–Blackwellized coefficients.
- `bb065c6` O2, with O3's basis half.
  - `BassComponents.Fit` and `Assemble` return a `BassBasisFit`, and `BassJointDraws.Pair` pairs the draws.
  - `BassBasisFit.Predict` returns a `BassBasisPrediction`, with the covariance in physical units, quantiles and the mean surface.
- `502fb0f` O4: `BassBasisFit.PointwiseLogDensity(units, BasisDensitySpace)`. Coefficient space is the default; the diagonal Gaussian in physical units is the declared option, through `DiagonalPlusLowRank`.
- `dd533be`: the O-O1 parity tests. `2fec1ee` carries O-O4's.
- `ac4730e`: the review's fixes.

**Tests.**
- New:
  - `BassOutputDecompositionTests` (19);
  - `BassPredictTests` (13);
  - `BassOutputComponentsTests` (7);
  - `BassPredictBasisTests` (8);
  - `BassOutputDensityTests` (3);
  - oracle `BassOutputTruncationParityTests` (5);
  - oracle `BassOutputDefectsParityTests` (7).
- Suites:
  - Maths.Regression.Tests 897/898 at `ac4730e`. The one failure is `BassSensitivityMomentTests.Resolution_RefusesUnknownNamesAndWrongKinds`, in the stage-6 chip's uncommitted tests that the build picked up; every stage-5 test passes.
  - Maths.Oracle.Tests' BASS parity 47/47, with R present.
  - The oracle run built into a second root, `artifacts\chip-stage5-oracle\`, so as not to contend with the background build. It is deleted with `chip-stage5\`.
- Evidence:
  - O-O4: rank, ranks and both residual sets match the oracle to 2·10⁻¹⁴ at n = 2, 5, 10 and 40, and at full retention. The stability ratios match too: 2.810, 1.654 and 1.096, null at n = 2 and at full retention.
  - O-O1: R's reconstruction misses Y by 3.09, this side's by 2·10⁻¹⁶. The skewness is 3.785, against R's −3.785.
  - O-O1: the mixture intervals are ±2.771 and ±2.644, against R's ±3.912 and ±5.201. The residuals were drawn from all 8 cases. The component noise has correlation 0.006 across cases, where R's is 1.
  - Prediction at the training inputs equals `BassDrawFit`'s fitted values bit for bit (scalar, categorical and direct functional).
  - The R̂-stopped components ran 800, 2000 and 2000 draws per chain and were paired over 800.

**Changes outside my files:** none. `BassFit.cs` gained `Predict` only, and `ThinSvd` is unchanged.

**Questions left as reversible defaults:**
1. *`Predict` takes the training data,* as `PointwiseLogDensity` does: the term table holds knot and level indices, and the fit keeps neither the knots nor the inputs. If the fit carried its prepared inputs, the argument would go.
2. *Basis intervals are weighted quantiles of the predictive draws,* the mixture's quantile up to Monte Carlo error. Direct fits solve their Gaussian mixture exactly. A basis mixture runs over draws × training residuals, so an exact root costs S·n per evaluation.
3. *Leave-one-out keeps the full data's y_s* under the column-sd option, so every residual is in one set of standardized units. Under a supplied centre it does not recentre. D07 states neither choice.
4. *Sign convention:* each component's largest |K| entry is positive.
5. *Even thinning* keeps draw ⌊(2j + 1)L/(2L')⌋, which is the identity when L' = L. Equal-weight components must share their chain count, and are refused with the reason otherwise, so the joint draws keep chains for PSIS.
6. *Weighted pairing:*
   - multinomial resampling, `SmcSpec`'s default, as a parameter;
   - to the smallest component's count, then shuffled, so that the sorted ancestors do not pair particle orders across components;
   - an equal-weight component beside weighted ones is thinned over its pooled draws.
7. *Seeds:* component k takes child k of the master seed, and the pairing's stream takes child r. With per-component configs, the seed tree replaces each config's `MasterSeed`.
8. *One bound:* max(1, ⌊P/C⌋) components run at once, each with min(P, C) chains. Progress is reported per component.
9. *Supplied basis:* its residuals are the in-sample ones from the training responses, which a basis fixed in advance did not see, or are supplied as given. Leave-one-out is not offered, because there is no construction to refit.
10. *Diagonal Gaussian density:*
    - the residual set's mean ē enters the mean, and its uniform-draw variance the diagonal;
    - a unit must hold whole cases;
    - it is refused where v = 0 (a constant column, or in-sample residuals at full retention) and for supplied residuals.
11. *Prediction columns* are matched by name; categorical values are given by level name.
12. *O-O1's constant-input fixture* still records R's training mapping to 0 as the "corrected" value, which predates D07's amendment. The parity test asserts the refusal. The fixture's field is the main session's to update.

**Review pass** (fresh eyes, after O4):
- *Correctness.* I re-read units and orientation at each step:
  - e is in standardized units, added before y_s;
  - the y-space density subtracts Σ log y_s;
  - the covariance is diag(y_s)·C·diag(y_s);
  - noise is drawn per draw and case (and per location, direct), and e whole;
  - the pairing and its thinning are within chains;
  - the axes are kept.

  The leave-one-out cost was stated as O(n⁴) independent of m: the refits are, but the map back costs O(nmρ). The docstring is corrected.
- *Coverage.*
  - Repeated singular values were untested. A tie is now compared by K Kᵀ and K W.
  - Already covered: n = 2, r = n − 1, full rank with m < n, one output location (decomposition and basis prediction), singleton draw and case, and a refusal with its reason for each input.

**For the main session:**
- *Docs* (`docs/Maths`):
  - leave-one-out in case coordinates under both centres;
  - the rank rule;
  - the joint draws (thinning, resampling and shuffling);
  - the predictive draw and its covariance;
  - the Rao–Blackwellized mean (ȳ + (ŷ_OLS − ȳ)/(1 + τ) at the training inputs under the mean centre);
  - the two density spaces.

  Checklist: O1–O4, and the O-X3 and O-X5 slots.
- *For O7:*
  - `BassBasisFit.Decomposition.Map` (y_m, y_s, K, `Physical`);
  - `Components[k]`, whose term tables share one transform and one set of knots;
  - `JointDraws.Draw(k, j)`, with the joint draws weighing equally;
  - `ComponentData(k)`, and the internal `Model(k)`.

  A note is in the shared session notes (`stage5-outputs-handover.md`).
