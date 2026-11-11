# Stage 6 chip (O7, basis sensitivity): handover (2026-10-29)

The full report is in `issues/regression/bass/stage6-basis-sensitivity-brief.md` under "Reports". Commits:
- `f1d53ab`: O7's code, tests and the O-O2 parity test;
- `4f6b4af`: the report.

- **Entry points:** `BassBasisFit.Sensitivity(spec)`, which forwards to `BassBasisSensitivity.Of(fit, spec)`, returns a `BasisSobolAnalysis`.
  - *What it keeps:* each joint draw's r × r component covariances on every subset, with the totals and means, in coefficient units and free of ℓ.
  - *Its views:* `AtLocation(ℓ)` (R's form, built on demand), `PerLocation()`, and `Integrated` (the location as an input named "location"), each a `SobolDecomposition` in output units.
- **Outside my files:**
  - `BassSensitivity.Subsets` takes the terms' variable lists, and `SubsetComparer` is internal. Both are behaviour-free.
  - The O-O2 oracle emits three more arrays it already computed.
- **The suites:** Regression 920/920; the BASS oracle suites 54/54, with R present. `artifacts\chip-o7\` is deleted.
- **Decisions to confirm:** the report lists 9 questions. The main ones:
  - one analysis with two views, rather than `Functional` choosing a form;
  - y_m included in the integrated location effect;
  - the location means' zero rule (standard deviation within the tolerance of the magnitude);
  - mixed hinge degrees refused.
