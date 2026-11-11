# Stage 6 chip (O5, O6, O8): handover (2026-10-29)

The full report is in `issues/regression/bass/stage6-sensitivity-brief.md` under "Reports". Commits:
- `17dfe2d`: Distributions' interval moments (shared);
- `2fa7ad6`: O5;
- `835a310`: O6 and O8;
- `853b61f`: the report.

- **Entry points** (namespace `Maths.Regression.Spline.Bass.Sensitivity`):
  - `BassSensitivity.Of(fit, data, spec)` returns a `SobolAnalysis`;
  - `BassSensitivity.Effects(fit, data, spec, effects)` returns an `EffectAnalysis`;
  - both take the training data and refuse other data by fingerprint.
- **Wiring is left to the main session:** `BassFit.Sensitivity` and `BassFit.Effects` forward to them. `BassFit.cs` was not touched.
- **For O7:**
  - `BassMoments.On(inputs, degree, ProductMeasure.Resolve(...))` takes terms by value, so terms from two component fits on one `BassInputs` share a service;
  - `Pair(a, b)` gives the covariance and its subset and variable shares;
  - O-O2's physical-units section is the check still owed.
- **A decision to confirm:** the normal's interval moments are the density's Taylor series on short stretches, not the Φ/φ recurrence. The recurrence loses every digit on narrow pieces (close knots). It is question 1 in the report.
- The suites: Distributions 148/148 and Regression 910/910, with the new parity tests green against R. `artifacts\chip-stage6\` is deleted.
