# GBASS noise, slice 3 and the envelope's guard: landed (chip, 2026-11-06)

Commits on `main`, nothing pushed: `f2e1e58` (the tangent hat's design-point guard), `5e837c4` (slice 3: log Bessel K, the GIG's normalizer and moments, the generalized hyperbolic, GIG mixing and the skewed t, the Normal-Wald on the core with β and γ drawn, the corrected `nw_est_mom` and `nw_gamma_prior`, tests and oracle modes), `4c8b587` (docs, ledger, slice 3's report). Report: `issues/regression/bass/gbass-noise-brief.md`, under "Reports".

For whoever builds slice 4 next:
- The record holds v, β, σ² and γ (`ObservationLatents`); a law's own parameter rides it as γ does, and `VarianceMeanMixture.MixingAt(latents)` reads v's law there. A new drawn scalar joins `BassNoiseShapeBlock` (scale slot, after σ²).
- σ² is on the record wherever `BassNoise.Skewed` (a skew drawn, or held away from 0); integrated otherwise. Generalized beta prime has β = 0, so it integrates σ² as t does unless a skew is declared.
- The readers take any law whose `At` returns an `IContinuousDistribution`; `PredictiveMixture.Quantile`'s law path is Newton now, so a CDF that is an integral costs a few reads per root. `GeneralizedHyperbolic` shows the pattern: an exact integral plus a lazily built normal-score table.
- The envelope accepts any reach now; a family's own good reach only saves steps.
- Trap: restoring a file from a backup with PowerShell `Copy-Item` keeps the backup's timestamp, and MSBuild then keeps a DLL built from the broken file; `touch` after restoring.

`artifacts\chip-gbass-noise\` was deleted.
