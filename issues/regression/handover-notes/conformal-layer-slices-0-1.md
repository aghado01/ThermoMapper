# Conformal layer, slices 0 and 1: hand-over (from the conformal chip, 2026-11-06)

Done on main, nothing pushed: 4fb54fe, cb48ac7, afb92fb, 58bb992, eebe51e (slice 0); 1b51458, 86aefe1, 15c29ce, 4b5150d (slice 1). The record is the two reports under "Reports" in `issues/regression/conformal-layer-brief.md`; the page is `docs/Maths/ConformalPrediction.md`.

For the GBASS session:
- `BassPointwise.Latents` now holds the latent list, and a noise law's latent per observation is reported `Integrated`. It used to be `IntegratedGivenOthers`, which made WAIC refuse every noise-law fit and PSIS-LOO's p_loo NaN.
- `BassFit.CandidateLogDensity` (`BassCandidates.cs`) calls `law.LogDensity` draws × candidates times. When `DensityByQuadrature` lands, it needs the same guard and `quadrature` flag as the pointwise path.
- Your uncommitted hunks in `BassFit.cs` and `BassPointwise.cs` were left out of my commits. Your working copies differ from HEAD by your hunks only, plus your wording of `Latents`' docstring.

For anyone timing: a reused MSBuild node keeps the priority it started with, and the test host it spawns inherits it. Set `MSBUILDDISABLENODEREUSE=1` in a below-normal shell so test hosts really run below normal.
