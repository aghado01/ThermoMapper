# The GIG's ratio-of-uniforms draw: landed (chip, 2026-11-06)

Commits on `main`, nothing pushed: `72fb8e2` (the ratio-of-uniforms primitive and Hörmann and Leydold's GIG draw, its tests, the envelope pin rebaselined, qbass's g-prior bar at 1.5), `a0f0292` (docs and report). Report: `issues/regression/bass/gbass-noise-brief.md`, under "Reports", newest first.

For whoever comes next:
- `GeneralizedInverseGaussian.DrawsByRatio(p, a, b)` picks the draw: ratio-of-uniforms for untruncated draws with a, b > 0, |p| <= 1e3 and 1e-8 <= omega <= 1e8; the tangent hat for truncated draws and beyond. Widening the range needs its law checked there first (`GigRatioOfUniformsTests`).
- `RatioOfUniforms` (internal, `Maths.Distributions`) takes any `IRelativeLogDensity` with closed or searched sides (`Bounds`).
- Quantile BASS's `qbass` readout: the g-prior arm's worst surface gap spans 0.49 to 1.26 over engine seeds 43 to 48; it is a gross-defect bar at 1.5, the ridge arm the parity bar at 1.0.
- `Maths.Distributions` now grants internals to `Maths.Regression.Tests` too (for timing the two draws side by side).

`artifacts\chip-gbass-noise\` was deleted.
