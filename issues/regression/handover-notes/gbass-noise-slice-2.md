# GBASS noise, slice 2: landed (chip, 2026-11-05)

Commits on `main`, nothing pushed: `5beaa49` (slice 2, the asymmetric Laplace and quantile BASS, the GIG and modified half-normal samplers, the retention lever, code and tests), `bfd58a9` (the new page `docs/Maths/VarianceMeanMixtures.md`, the ledger, slice 2's report). Report: `issues/regression/bass/gbass-noise-brief.md`, under "Reports".

For whoever builds slice 3 next:
- The noise's docs live on `docs/Maths/VarianceMeanMixtures.md` now, at the user's request: `BayesianRegression.md` is being broken up and is not to grow. Slice 1's section moved there; `BayesianRegression.md` keeps a pointer.
- σ² rides the record under skew (`ObservationLatents.Sigma2`); `BassTarget.GivenAt(fit)` reads it, and the state's prior adds ln p(σ²). σ²'s block is `BassScaleBlock`, the modified half-normal on σ⁻¹, with Q from `GPriorConditional.PriorEnergy` through the pencil.
- A free skew β and γ (slice 3) enter `VarianceMeanMixture`, whose `At` and `Reading` already take β from the record; `ExponentialMixing.SampleConditional` takes a drift; the GIG's conditional comes from `GigParameters`.
- Retention: draws of a model with `Rebuilds` keep `BassDraw.Latents` only; `WorkingOn`/`LatentsOn` rebuild the view by the same function that built it.
- Geweke harness: `BassQuantileQualificationTests`, σ² ~ IG(3, 2), τ ~ Gamma(5) or held.
- The `qbass` oracle at 6,000 iterations overruns `ROracle`'s 120 s at below-normal priority; it runs at 3,000.

`artifacts\chip-gbass-noise\` was deleted.
