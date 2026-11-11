# GBASS noise, slice 4: the generalized beta prime and the horseshoe (2026-11-06)

Landed: `4662639` (code, tests, the readout oracle's corrected hbass), `d48286a` (docs, ledger, the brief's report). Stopped after slice 4, as asked; slice 5 waits for the user.

- **What's on the config:** `BassNoise.Horseshoe()`, `StrawdermanBerger()`, `GeneralizedBetaPrime(GbpParameters)`. hbass = Horseshoe() with `Sigma2 = Sigma2Prior.HalfCauchy()` (c = 1; the engine's σ² is hbass's w·var(y)).
- **Seam change:** `IMixingLaw.SampleConditional` → `StepConditional(current, …)`. Exact laws answer with their draw, so their numbers are unchanged. GBP takes one slice step on ln v.
- **Pointwise density under GBP is opt-in:** `fit.PointwiseLogDensity(data, quadrature: true)`. It costs 41 µs a read against t's 0.1. The conformal chip's candidate grid reads `law.LogDensity` per point with no opt-in; I told that session.
- **The mixture's integral is now shared** (`NormalVarianceMeanMixture`), with three upgrades that move GH's predictive tables within 3·10⁻⁴ of a log tail: the integral about its peak, the table in asinh(x/s), and the steps into a pole.
- **GBASS defect:** hbass's adapted v step (0.12 at acceptance 0.98 at N = 200) mixes v slowly; M reads 4.16 against 5.36 with the step held. The readouts compare against the held-step sampler (`prop_sigma=2.5`).
- **Bars:** hbass's readout bars are qbass's: ridge ≤ 1.0, g-prior ≤ 1.5, from a six-seed spread (0.52–0.82 and 0.55–1.20).
- **Coordination with the conformal chip:**
  - It committed `BassPointwise.Latents` and the `Integrated` label itself (1b51458). My commit carries only the `quadrature` hunks in BassPointwise.cs and BassFit.cs.
  - While its test file didn't compile, I ran my gates from a scratch test project under my artifacts folder (now deleted).

Follow-ups, each with its decision, are in the brief's slice 4 report:
- one standardized table across draws (the predictive's 2,000 tables);
- the horseshoe's E₁ closed density;
- auxiliary-gamma Gibbs for p = 1;
- a stated scale for the hierarchical σ² prior, so Geweke can check it.

**Addendum.** The conformal chip handed over the candidate seam's opt-in once DensityByQuadrature reached HEAD: `5ec3a9b` (BassCandidates.LogDensity refuses a quadrature law without quadrature: true; BassFit.CandidateLogDensity takes the flag), `9e57141` (docs). The conformal chip has finished and doesn't write to the checkout.
