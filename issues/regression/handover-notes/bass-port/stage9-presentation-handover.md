# Stage 9 chip (O9, the presentation data): handover (2026-10-29)

The full report is in `issues/regression/bass/stage9-presentation-brief.md` under "Reports". Commits:
- `6d37bab`: `DrawSample.Quantile` is type 7 under equal weights at any scale (shared, Samplers);
- `c564f56`: one weighted quantile; `PredictiveMixture`'s empirical (type-1) overload deleted, basis intervals through `DrawSample`;
- `f6652b8`: basis effect functions (`BassBasisFit.Effects`, `BasisEffectAnalysis`);
- `4892ae2`: the presentation records (`bass/presentation/`, namespace `...Bass.Presentation`) and the views on both fits.

- **Views:** `Traces()`, `TrainingFit(...)`, `Metadata()`, `Summary(probabilities)`, `SensitivitySummary(...)` on `BassFit` and `BassBasisFit`; `Effects(spec, effects)` on `BassBasisFit`.
- **Moved results:** basis predictive quantiles go from type 1 to type 7, within Monte Carlo error. Nothing else moved: BARS's peak interval feeds `DrawSample` unit weights, and the sensitivity code is untouched.
- **Docs owed:** `docs/Maths/BayesianRegression.md` needs the presentation section and the basis effects' formulas, and its lines 982 and 1037 updated. `docs/Maths/Samplers.md:351` needs `DrawSample`'s equal-weight rule.
- **Decisions to confirm:** the report lists 11. The main ones:
  - the default probabilities (0.025, 0.25, 0.5, 0.75, 0.975);
  - indices summarized over the defined draws, variances over all draws;
  - `NoiseSd` = E[σ] (R's code);
  - basis intervals from the predictive draws, by seed.
