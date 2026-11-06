# Variance-mean mixture noise

> GBASS's noise on the BASS core (E1): a normal variance-mean mixture, εᵢ = σ(β·vᵢ + √(c·vᵢ)·ζᵢ), Gaussian given a latent scale per observation, which the Gaussian core fits through a working response that moves with the latents. Student t (robust BASS) and the asymmetric Laplace (quantile regression) are built; the Normal-Wald, generalized beta prime mixing and the engine's samplers on the carrier follow.

## Overview

GBASS (Rumsey, Francom and Shen, arXiv 2306.01911; CRAN 2.0.1, the oracle) is BASS with a different noise, not a new engine. This page is the noise's: what it declares, how the sweep moves its latents, what the readers of a fit take, and what was checked. The core it rides is [BayesianRegression](./BayesianRegression.md)'s: the [extension seams](./BayesianRegression.md#bass-the-extension-seams) it fills, the [weighted core](./BayesianRegression.md#bass-the-weighted-core), the [priors that stop reading a moving response](./BayesianRegression.md#bass-priors-that-stop-reading-a-moving-response) and the [design whose working response moves](./BayesianRegression.md#bass-a-design-whose-working-response-moves). The laws it draws from are [ProbabilityAndInformation](./ProbabilityAndInformation.md)'s: the [mixing law](./ProbabilityAndInformation.md#normal-variance-mean-mixtures-the-mixing-law-e1) and [the GIG, the modified half-normal and the asymmetric Laplace](./ProbabilityAndInformation.md#the-gig-the-modified-half-normal-and-the-asymmetric-laplace-e1). The build is [the brief's](../../issues/regression/bass/gbass-noise-brief.md), slice by slice, with the oracle modes of slice 0 under `r/oracles/gbass_*_oracle.R`.

| Noise | Mixing | c | β | σ² | Built |
|---|---|---|---|---|---|
| Student t (GBASS's `tbass`) | IG(ν/2, ν/2) | 1, or (ν − 2)/ν | 0 | Integrated in the structure step, redrawn with a | Slice 1 |
| Asymmetric Laplace (`qbass`) | Exp(1) | 2/(q(1 − q)) | 1/q − 1/(1 − q), fixed | On the record, drawn by the modified half-normal | Slice 2 |
| Normal-Wald (`nwbass`) | IG(1/γ, 1), γ sampled | 1 | Sampled | On the record | Slice 3 |
| Beta prime (`hbass`) | GBP(1, ½, ½) | var(y) | 0 | Half-Cauchy on σ | Slice 4 |

## Student t, robust BASS (E1's slice 1)

Built 2026-11-05 ([the build brief](../../issues/regression/bass/gbass-noise-brief.md#reports), slice 1), GBASS's `tbass` on the BASS core, and the first member of E1's family: εᵢ = σ(β·vᵢ + √(c·vᵢ)·ζᵢ), ζᵢ ~ N(0, 1), with β = 0 and vᵢ ~ IG(ν/2, ν/2), so ε/(σ√c) is t_ν. σ² is the core's, GBASS's w. The core fits y under precision weights 1/(c·v) at the current latents, D03's robust row: the coefficient prior on the unweighted columns, so the prior reads no W and v's conditional stays closed form. No new distribution was needed: v's conditional is an inverse gamma.

**Declared on the config.** `BassConfig.Noise = BassNoise.StudentT(ν, c, reference)` builds the observation model (`VarianceMeanMixture`) inside `BassModel.Create`, so `Bass.Fit`, `Resume`, `Fingerprint` and every view of a fit read it. `Noise` is left out of the JSON when null, so a Gaussian config's form and fingerprint stand as they were.

| c | σ² is | Matches |
|---|---|---|
| 1, the default | The t's squared scale | BARS's `StudentTWeights` |
| (ν − 2)/ν | The noise variance, for ν > 2 | GBASS's `tbass` |

c is a units convention: under a σ² prior that reads no scale it moves no structure, and under the scaled prior it sets what s² is compared with.

**What the priors read** (`BassNoise.Reference`, frozen at the start on `BassPriorReference`).

| Reference | Reads | For |
|---|---|---|
| `Moments`, the default | ȳ and the mean square about it | D03's read: at β = 0 under the unweighted Gram only the weights move, which no prior reads ([the options](../../issues/regression/bass/moving-response-priors-options.md)), so this restates what Gaussian BASS reads and large ν meets Gaussian BASS |
| `Robust` | The median and (MAD/Φ⁻¹(¾))² | Heavy tails, where an outlier would move the level and the scale's read |
| `None` | Nothing | Priors that read no response (a zero or flat intercept with a σ² prior that reads no s²); refused otherwise, at validation |

**The sweep** (D04's robust row): the structure with σ² and a integrated out, given v; τ; then σ² and a from their conditional given the structure, τ and v (`BassCoefficientsBlock`, the coefficients' slot); then each vᵢ from IG((ν + 1)/2, (ν + rᵢ²/(σ²c))/2), rᵢ = yᵢ − fᵢ (`BassLatentScaleBlock`, the latent slot), and the chain follows the working response at the new v. The scale slot stays empty: under a pure scale mixture σ² is redrawn with a each sweep. This is a partially collapsed Gibbs scan (van Dyk and Park 2008): σ² and a, integrated out of the first two steps, are drawn right after them and before the step that reads them. A retained draw composes σ² and a afresh on its own working response, as Gaussian BASS's does, so the trajectory does not depend on how often draws are kept.
- *The latents ride the working response* (`ObservationLatents` on `WorkingResponse.Latents`; E1's decision 1): v, and β, which is 0 here. A retained draw keeps the working response it was composed on, so it keeps its v, the generative state.
- *σ² and a ride the chain's own moves within a sweep* (`BassStructureMoves.SweepFit`), from the coefficients' block to the last block that reads them. Decision 1 put them on the record; a working response is shared by every chain at the start and by retained draws, so a value that changes within a sweep cannot ride it without a copy per sweep.
- *F under a moving W* stays on the gate that reads W (decision 4): a v that puts the current structure outside F is refused and the chain keeps its state, the exact step. The scan counts each block's attempts and acceptances, and a fit reports them per chain (`BassChainAccount.Blocks`).
- *Refused until their slices:* a tempered rung, SMC and the flat histogram, whose particles and walkers would hold their own latents (slice 5); a direct functional response, a per-observation latent there (slice 5); a skewed t, whose drift needs the GIG with a > 0 in v's conditional and Bessel K in its marginal (slice 3). Each refusal names its slice.

**The mixing law** (`IMixingLaw`, `StudentTMixing`; Maths.Regression, shared with BARS). Per observation: v's prior draw, its conditional at a power b, ∝ p(v)·N(r | m·v, s²·v)^b with s² = σ²c and m = σβ, the marginal ln ∫ N(r | m·v, s²·v)p(v)dv, and E[v] and Var[v], +∞ where infinite. Under t the conditional is drawn as 1/λ with λ ~ Gamma((ν + b)/2, (ν + b·r²/s²)/2), by the same call `StudentTWeights` makes, so at b = 1 the two give the same numbers from the same stream. The marginal is t_ν(0, s), its ln Γ((ν + 1)/2) − ln Γ(ν/2) taken from Stirling's series term by term above ν/2 = 10³ and ln(1 + r²/(νs²)) by `log1p`, so neither loses its digits as ν grows. `StudentTWeights.MarginalLogDensity` takes the direct forms, which on the same sweep read up to 2.0·10⁻⁵ off, relative, at ν = 3·10¹⁰; slice 6 routes it through this law.

**The noise's law for the readers** (`IBassObservationModel.Law`, null for the Gaussian, whose readers keep their own code):
- *Pointwise* (OQ17's scope rule): t_ν(yᵢ − fᵢ; 0, σ√c), each vᵢ integrated out, reported as an observation-level latent integrated. The collapsed predictive is refused under a law until it has its own check.
- *Prediction's three readings:* the surface f; y's mean f + σβE[v] (`BassPrediction.MeanShifts`, `ObservationMean`), f itself under t; predictive draws with v from its prior, then ζ. The predictive mixture's quantiles solve Σₛ wₛFₛ(q − fₛ) = p with each draw's t (`PredictiveMixture.Quantile` over `IContinuousDistribution`s).
- *The training fit's noise sd* is the posterior mean of √(σ²(β²Var[v] + c·E[v])), null where that variance is infinite: t at ν ≤ 2, where GBASS reports s2 = −w or NaN ([the corrections](./BayesianRegression.md#corrections-to-the-sources)).

**Checked** (`BassStudentTTests`; behind `BASS_QUALIFICATION`, `BassStudentTQualificationTests`; against R, `GbassConditionalParityTests` and `GbassStudentTReadoutTests`).
- *The law.* Its marginal against quadrature of the mixture over ln v on 11 chosen cases and a seeded sweep of 400, ν from 0.05 to 10¹¹ and s² over 24 decades: worst relative gap 1.4·10⁻¹³. Its conditional's draws against the Gamma's mean and log mean at b = 1, 0.4 and 0 (the prior), within 5 standard errors. The log-gamma ratio across its switch to 11 places.
- *v's conditional is `StudentTWeights`' draw*, number for number from one stream, at ν = 0.5, 3, 4 and 30; and GBASS's own v conditional under t mixing is that Gamma at every observation of a frozen state, shape exactly, rate to 14 digits.
- *Geweke's joint-distribution test* on N = 12, one predictor, M ≤ 2, with proper priors that read no response (the zero centre, σ² ~ IG(3, 2), t(10)): 200,000 draws of the marginal-conditional simulator against 300,000 steps of the successive-conditional one, each step a sweep on the newly drawn y. Eight functions (M, ln τ, ln σ², ln v₁, mean ln v, arctan f at two inputs, the standardized squared residual) on four cases, c = 0.7 on two seeds, c = 1, and τ held: every |z| ≤ 1.78. With τ ~ Gamma(5, rate 5N) or held: under Zellner–Siow's shape ½ the successive chain sticks in ln τ's long left tail, its mean 0.30 above the prior's on one seed and 0.58 below on another, every other function in place, which is the test's reach rather than the sweep's. Its power: dropping c from the v block's unit variance reads z 75 on mean ln v and −12 on the squared residual at c = 0.7, and is invisible at c = 1, which is why c = 0.7 is in the set. Nothing was refused: no structure left F at a drawn v, and no latent draw of 300,000 was refused.
- *Large ν reaches Gaussian BASS.* ν = 10⁷ against Gaussian BASS, 16 chains each, N = 200: σ² 0.23295 against 0.23268 (z −0.47), M 4.275 against 4.239 (z −0.92), the fit at all 200 cases rms z 1.06.
- *`tbass`'s readouts* at ν = 3, c = ⅓, N = 200 ([the oracle mode](../../r/README.md)): w 0.230 [0.180, 0.289] under the g-prior and 0.235 under the ridge, against `tbass`'s 0.224 [0.174, 0.283]; the surface's means within 0.23 of the wider 90% half-band at worst over 15 grid points (0.34 under the ridge); M 4.2 and 4.6 against 5.3. The chains differ by design: continuous knots and whole-term mutation in GBASS, pseudo-counts 1 against 5, its raw ridge with the intercept inside, its floor on w, λ sampled.
- *Recovery from noise outside the model*, a contaminated normal 0.85·N(0, 0.25²) + 0.15·N(0, 2²) on sin(2πx₁) + x₂, N = 300, six data sets: t(3) BASS's 90% intervals cover 89.5% of 294 grid points of the true surface, RMSE 0.059, against Gaussian BASS's 86.4% and 0.142; its predictive 90% intervals cover 88.3% of 1,800 held-out responses, Gaussian's 92.7%.
- *The run.* A fit resumes bit for bit, the latents of every retained draw included; the pointwise density, the predictive draws' deciles and the mixture's quantile against t's own; the noise sd against its closed form and null at ν = 2. Gaussian BASS's eight pinned hashes and the seams' suites hold (316 tests).

**Cost** (`BassStudentTQualificationTests.CostRow_PerPartOfASweep_AndRetention`; the thread's cycles on a loaded machine, so indicative, ratios before values). µs per sweep, each part's mean over the timed sweeps, t(4) on the cost report's rich surface, beside a Gaussian chain's whole sweep at the same structure:

| Point | Form | M | Structure | τ | a | v's N draws | Follow | Sweep | Gaussian sweep |
|---|---|---|---|---|---|---|---|---|---|
| n 1,000, p 10, the base | g | 20 | 261 | 6.5 | 46 | 114 | 533 | 961 | 40 |
| | ridge | 20 | 245 | 6.9 | 40 | 56 | 476 | 824 | 71 |
| n 200 | g | 7.7 | 23 | 2.4 | 5.1 | 11 | 41 | 83 | 2.9 |
| | ridge | 7.7 | 26 | 2.8 | 5.4 | 11 | 45 | 91 | 15 |
| n 10,000, M_max 200 | g | 200 | 96,811 | 30 | 5,079 | 544 | 285,544 | 388,008 | 1,798 |
| | ridge | 200 | 103,365 | 24 | 4,442 | 560 | 285,784 | 394,176 | 3,716 |

- *Follow is where a sweep goes:* 55% at the base, 74% at M 200, refactoring the weighted Gram at O(NM²). Then the structure step, whose proposals each decompose the pencil at O(p³) under the unweighted prior, then the N draws of v.
- *So a t sweep is 24 Gaussian sweeps at the base and 29 at n 200, and 216 at M 200.* The direct functional point is refused in this slice (a per-observation latent there is slice 5), and the weighted prior Gram is no valid chain under a latent scale without the Metropolis correction of [the options' step 3](../../issues/regression/bass/moving-response-priors-options.md#what-to-build-in-order), so only the unweighted reading ran.
- *The lever after the row,* structure sweeps per latent move, stays a measured option; the default stays one each, D04's.

**Retention.** Read here, built with slice 2: a retained draw keeps its latents alone ([below](#retention-a-draw-keeps-its-latents)). Under t its working response held the weights and v, 2N doubles, 15.6 KB a draw at N = 1,000; v alone is half that.

## The asymmetric Laplace, quantile regression (E1's slice 2)

Built 2026-11-05 ([the build brief](../../issues/regression/bass/gbass-noise-brief.md#reports), slice 2), GBASS's `qbass` on the BASS core: v ~ Exp(1), c = 2/(q(1 − q)) and β = 1/q − 1/(1 − q), so ε = σ(βv + √(cv)·ζ) is ALD(0, σ, q) and the surface f is the q-quantile of y (Kozumi and Kobayashi 2011). β and c are the level's, fixed; a free skew is slice 3's. `BassConfig.Noise = BassNoise.AsymmetricLaplace(q, reference)`.

**σ² rides the record.** Under a skew the working response is y − σβv, so σ² sits in the offset and the structure step cannot integrate it. It is drawn apart and held on the latent record (`ObservationLatents.Sigma2`; E1's decision 1, w included), and the target reads it from the state as it reads the statistics (`BassTarget.GivenAt`): the structure's marginal given σ², with ln p(σ²) in the state's prior (the hierarchical prior's at its rate). `BassTarget.Given` stays for a σ² from outside (E32). A held σ² is the config's, on the record from the start, and no block moves it.

**The sweep** (D04's robust row): the structure with a integrated out, given (σ², v); τ; a given σ² (`BassCoefficientsBlock`); each vᵢ from GIG(½, 2 + β²/c, rᵢ²/(σ²c)), rᵢ = yᵢ − fᵢ (`BassLatentScaleBlock`, the drift σβ in the mixing law's conditional); then σ² (`BassScaleBlock`), and the chain follows the offset, the weights shared. A Gibbs scan: each block draws its full conditional, a and v from the latest of the others.
- *σ²'s conditional* (decision 2: v and the modified half-normal first). With the coefficients' prior variance scaling with σ², K of them inside the prior with energy Q = τ(a − a₀)ᵀD(a − a₀) (`GPriorConditional.PriorEnergy`, through the prior's pencil, which every latent noise's prior holds) and σ²'s prior IG(α₀, β₀), u = σ⁻¹ is MHN(N + K + 2α₀, ½(Σrᵢ²/(cvᵢ) + Q + 2β₀), βΣrᵢ/c). K is M + 1 with the intercept inside the prior and M under the flat intercept. Jeffreys is α₀ = β₀ = 0.
- *Three routes to σ², for the record* (E1, the user's derivation of 2026-10-26): this one, in v, scale-equivariant; the GIG on u = √σ in the z parameterization, exact under any GIG prior on u but with the coefficients' variance ∝ u, which is not equivariant, so τ carries y's units; and a slice step in z for the equivariant prior. Interweaving at variance ∝ σ² pairs the first with the third, slice 7's.
- *What GBASS does differently,* corrected here (slice 0's findings): its w step at β ≠ 0 drops w's prior a and its floor. The engine's σ² priors are inverse gammas, whose shape and scale both enter, and it keeps no floor.

**What the priors read:** `Quantile`, by default, the sample q-quantile and the squared mean check loss (`BassPriorReference.Quantile`), the asymmetric Laplace's level and scale by maximum likelihood, whose ȳ is not its level. `Moments`, `Robust` and `None` stay selectable.

**The noise's law for the readers:** pointwise, ALD(yᵢ − fᵢ; 0, σ, q), the exponential mixture's marginal at the drift (`ExponentialMixing.MarginalLogDensity`); prediction's three readings, the surface f, which is the q-quantile, y's mean f + σβ (`MeanShifts`), and predictive draws through v and ζ, the mixture's quantiles over each draw's ALD; the training fit's noise sd, √(σ²(β² + c)), and its `Offset`, the draws' mean σβv.

**The primitives** (`Maths.Distributions`; [their forms](./ProbabilityAndInformation.md#the-gig-the-modified-half-normal-and-the-asymmetric-laplace-e1)): the GIG's sampler, the tangent hat in ln x where every GIG is log-concave; the modified half-normal's three exact draws; the asymmetric Laplace in closed form.

**Retention:** built here, [below](#retention-a-draw-keeps-its-latents).

**Refused until their slices:** a free skew (slice 3), a tempered rung, SMC and the flat histogram, and a direct functional response (slice 5), and the collapsed predictive under a noise law.

**Checked** (`BassQuantileTests`; behind `BASS_QUALIFICATION`, `BassQuantileQualificationTests`; against R, `GbassPrimitiveParityTests` and `GbassQuantileReadoutTests`).
- *The primitives,* against their laws integrated here and against `gbass_primitives_oracle.R`'s: the GIG on 134 cases, 14 chosen and 120 swept with p over ±10^(−3…3) and a and b over 16 decades each, 4,000 draws a case, worst decile z 3.01, and on four intervals; the modified half-normal on 134 cases over its three branches, worst decile z 3.32, 2.35 and 2.57; against R's law, 11 GIG and 9 modified half-normal cases at 50,000 draws, worst |z| 2.56 and 2.55, means included. The sweep found what the chosen cases had not, each corrected: the GIG's first log ratio cancelling where ω ≪ |p|, its hat's reach landing where the far tangent swamped it, and the half-normal's left reach landing on x = 0.
- *The asymmetric Laplace* is the exponential mixture's marginal to 12 digits at q = 0.9, ½, 0.1 and 0.01; its tails sum to 1 and its quantiles round-trip; its draws carry the law, on intervals too.
- *Exponential mixing:* its marginal against quadrature of the mixture over 305 cases, any drift, worst relative gap 8.7·10⁻¹⁵; its conditional's draws against the GIG at b = 1 and ½.
- *Geweke's joint-distribution test* under the skew, on the slice-1 design with σ² ~ IG(3, 2): four cases (q = 0.3 on two seeds, q = 0.8 with τ held, and q = ½, where the drift vanishes and σ² is integrated), every |z| ≤ 2.95, nothing redrawn or refused. Its power: flipping the skew's sign in σ²'s block sends σ² off to infinity within seconds and fails every skewed case; q = ½, where the break cannot show, passes.
- *The skew's sign* (E1's fixture): on symmetric noise, N = 400, the held-out share of 2,000 responses below the surface is 0.903 at q = 0.9 and 0.106 at q = 0.1, which a flipped sign would swap.
- *The quantile recovered* on noise outside the model, f + (0.3 + 0.4x₂)(G − 2) with G ~ Gamma(2, 1), skewed and heteroscedastic, N = 400, four data sets: the held-out share below the surface is 0.518 at q = ½ and 0.884 at q = 0.9; the surface's RMSE against the true quantile on a 7 × 7 grid is 0.160 and 0.325; f's 90% intervals cover the true quantile at 76.5% and 63.8% of the points, the working likelihood's undercoverage (Yang, Wang and He 2016).
- *`qbass`'s readouts* (`GbassQuantileReadoutTests`, q = 0.9, N = 200, w's prior IG(0.1, 0.05) on both sides, R 3,000 iterations, the engine 4 chains × 1,500 draws): σ² 0.00768 [0.00587, 0.00990] under the g-prior and 0.00820 [0.00640, 0.0103] under the ridge, against `qbass`'s w 0.00800 [0.00624, 0.00991]; the surface's means within 0.49 of the wider 90% half-band at worst over 15 grid points (0.65 under the ridge); M 4.94 and 3.96 against 4.87.
- *The run:* the record's σ² is the target's, its marginal given σ² to 12 digits over 100 sweeps; the offset is σβv; a quantile fit resumes bit for bit, latents and σ² included; the readers give y's mean shift σβ, the q-quantile at the surface, the asymmetric Laplace's pointwise density and noise sd, and the draws' mean offset.

**Cost** (indicative, the thread's cycles on a loaded machine). `BassQuantileQualificationTests.CostRow_PerPartOfASweep_AndRetentionBeforeAndAfter`, q = 0.9 on the cost report's rich surface; µs per sweep, each part's mean, beside a Gaussian chain's whole sweep at the same structure. σ²'s draw is its block less a second, offset-only move timed apart.

| Point | Form | M | Structure | τ | a | v's N draws | Follow v | σ² | Follow σ² | Sweep | Gaussian sweep |
|---|---|---|---|---|---|---|---|---|---|---|---|
| n 1,000, p 10, the base | g | 20 | 234 | 4.9 | 45 | 952 | 529 | 9.0 | 36 | 1,810 | 45 |
| | ridge | 20 | 260 | 5.4 | 52 | 975 | 554 | 14 | 34 | 1,894 | 111 |
| n 200 | g | 8.9 | 42 | 3.7 | 7.3 | 219 | 77 | 2.5 | 8.7 | 360 | 5.1 |
| | ridge | 5.9 | 16 | 2.7 | 4.3 | 179 | 37 | 2.1 | 5.0 | 246 | 16 |
| n 10,000, M_max 200 | g | 199 | 92,345 | 22 | 6,580 | 9,094 | 268,684 | 1,362 | 3,134 | 381,220 | 1,651 |

- *v's N GIG draws are the largest part at n 1,000,* 53% of a sweep, about 1 µs a draw through the hat; Follow after them is 29%. At M 200 Follow is 70%, refactoring the weighted Gram at O(NM²).
- *σ²'s move is cheap:* the offset alone moves, so its Follow keeps the factor, O(NM), and the modified half-normal is one draw.
- *So a quantile sweep is 40 Gaussian sweeps at the base,* against t's 24, the GIG's draws being the difference.

## Retention: a draw keeps its latents

Read in slice 1, built in slice 2: the brief's lever for a working response that dominates a retained draw. A model that rebuilds its working response from its latents (`IBassObservationModel.Rebuild`, `VarianceMeanMixture`'s) has its retained draws keep `BassDraw.Latents`, v with β and σ², and rebuild the working response on demand (`WorkingOn`); the readers that need the latents alone take them without a rebuild (`LatentsOn`). The blocks and the rebuild build the working response by one function from the latents, so a view rebuilt from v is the one the chain held, bit for bit. The generative state stays: v is kept. Measured over 200 retained draws of one quantile chain (`BassQuantileQualificationTests.Retention_BeforeAndAfter`), held both ways: 313.4 KB a draw with its working response against 70.4 KB with its latents alone at N = 10,000 (4.45×; the weights, the offset, the shifted response and v against v, 4N doubles against N); 32.7 KB against 9.0 KB at N = 1,000 (3.6×), where the draw's own few hundred bytes and the collector's noise weigh more. Under t the working response held the weights and v, so the ratio there is 2.
- *The alternatives:* keeping the working response, which holds the weights, the offset and the shifted response beside v, up to 4N doubles under the skew; or keeping nothing and redrawing, which loses the generative state.
- *The cost:* a rebuild is O(N) per draw read, for the views that read the offset or the weights (the residual check, the training fit's offset, the conditional means).

## Status and caveats

- *Built:* slices 0 to 2 of [the brief](../../issues/regression/bass/gbass-noise-brief.md), each with its report there.
- *Open, by slice:* a free skew and the Normal-Wald's γ, with Bessel K of real order and the corrected `nw_est_mom` and `nw_gamma_prior` (3); generalized beta prime mixing and its slice steps (4); the latents on tempered rungs, SMC and the flat histogram, and a per-observation latent under a direct functional design (5); BARS through the shared law (6); the z parameterization and interweaving (7).
- *The asymmetric Laplace is a working likelihood* where the noise is not one: its posterior intervals for the q-quantile undercover under other noise (63.8% for nominal 90% at q = 0.9 on skewed, heteroscedastic noise), as Yang, Wang and He (2016) show in general. The surface is the q-quantile; its intervals are not yet calibrated for it.

## Source map

| Area | Files |
|---|---|
| The noise on the config, the model, its blocks | `src/maths/regression/spline/bass/BassNoise.cs` (`BassNoise`, `ObservationLatents`, `IBassNoiseLaw`, `VarianceMeanMixture`, `BassCoefficientsBlock`, `BassLatentScaleBlock`, `BassScaleBlock`) |
| The carrier and the seam | `BassObservationModel.cs` (`WorkingResponse.Latents`, `IBassObservationModel.Law`, `Rebuilds`, `Rebuild`) |
| σ² read from the state | `BassTarget.cs` (`GivenAt`) |
| The sweep's draw of a and σ² | `BassSweep.cs` (`BassStructureMoves.DrawSweepFit`, `SweepFit`) |
| The draw keeping its latents | `BassState.cs` (`BassDraw.Latents`, `WorkingOn`, `LatentsOn`) |
| The prior's energy for σ²'s block | `src/maths/regression/gprior/CollapsedGPrior.cs` (`GPriorConditional.PriorEnergy`), `PriorPencil.cs` |
| The mixing laws | `src/maths/regression/MixingLaw.cs` (`IMixingLaw`, `StudentTMixing`, `ExponentialMixing`) |
| The laws drawn from | `src/maths/distributions/GeneralizedInverseGaussian.cs`, `ModifiedHalfNormal.cs`, `AsymmetricLaplace.cs` |
| Tests | `tests/regression/BassStudentTTests.cs`, `BassStudentTQualificationTests.cs`, `BassQuantileTests.cs`, `BassQuantileQualificationTests.cs`; `tests/oracle/Gbass*.cs` |
