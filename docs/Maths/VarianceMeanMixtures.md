# Variance-mean mixture noise

> GBASS's noise on the BASS core (E1): a normal variance-mean mixture, εᵢ = σ(β·vᵢ + √(c·vᵢ)·ζᵢ), Gaussian given a latent scale per observation, which the Gaussian core fits through a working response that moves with the latents. Student t (robust BASS), the asymmetric Laplace (quantile regression), the Normal-Wald with its skew and tails drawn, the skewed t, the generalized hyperbolic and the generalized beta prime with the horseshoe are built; the engine's samplers on the carrier follow.

## Overview

GBASS (Rumsey, Francom and Shen, arXiv 2306.01911; CRAN 2.0.1, the oracle) is BASS with a different noise, not a new engine. This page is the noise's: what it declares, how the sweep moves its latents, what the readers of a fit take, and what was checked. The core it rides is [BayesianRegression](./BayesianRegression.md)'s: the [extension seams](./BayesianRegression.md#bass-the-extension-seams) it fills, the [weighted core](./BayesianRegression.md#bass-the-weighted-core), the [priors that stop reading a moving response](./BayesianRegression.md#bass-priors-that-stop-reading-a-moving-response) and the [design whose working response moves](./BayesianRegression.md#bass-a-design-whose-working-response-moves). The laws it draws from are [ProbabilityAndInformation](./ProbabilityAndInformation.md)'s: the [mixing law](./ProbabilityAndInformation.md#normal-variance-mean-mixtures-the-mixing-law-e1) and [the GIG, the modified half-normal and the asymmetric Laplace](./ProbabilityAndInformation.md#the-gig-the-modified-half-normal-and-the-asymmetric-laplace-e1). The build is [the brief's](../../issues/regression/bass/gbass-noise-brief.md), slice by slice, with the oracle modes of slice 0 under `r/oracles/gbass_*_oracle.R`.

| Noise | Mixing | c | β | σ² | Built |
|---|---|---|---|---|---|
| Student t (GBASS's `tbass`) | IG(ν/2, ν/2) | 1, or (ν − 2)/ν | 0 | Integrated in the structure step, redrawn with a | Slice 1 |
| Asymmetric Laplace (`qbass`) | Exp(1) | 2/(q(1 − q)) | 1/q − 1/(1 − q), fixed | On the record, drawn by the modified half-normal | Slice 2 |
| Normal-Wald (`nwbass`) | IG(1/γ, 1), γ drawn or held | c, GBASS's 1 | Drawn or held | On the record under a skew; integrated where β is held at 0 | Slice 3 |
| Skewed t | IG(ν/2, ν/2) | as t | Drawn or held away from 0 | On the record | Slice 3 |
| Generalized hyperbolic (`gbass`, GIG `v_prior`) | GIG(p, a, b) | c | Drawn or held | As the Normal-Wald | Slice 3 |
| Horseshoe (`hbass`) | GBP(1, ½, ½), √v half-Cauchy | 1 here; `hbass`'s var(y) | Drawn or held (0 by default) | Half-Cauchy on σ at the reference's scale, by the hierarchical route; integrated where β is held at 0 | Slice 4 |
| Generalized beta prime (`gbass`, GBP `v_prior`) | GBP(p, a, b) | c | Drawn or held | As the horseshoe | Slice 4 |

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

## The Normal-Wald, the skewed t and the generalized hyperbolic (E1's slice 3)

Built 2026-11-06 ([the build brief](../../issues/regression/bass/gbass-noise-brief.md#reports), slice 3), GBASS's `nwbass` on the BASS core and the noise's general case: v ~ GIG(p, a, b), whose mixture is the generalized hyperbolic (Barndorff-Nielsen 1977), with the skew β drawn or held, so the data can choose the noise's skew and tails.

| Kind | v's law | β | γ | On the config |
|---|---|---|---|---|
| Normal-Wald (`nwbass`) | GIG(−½, γ², 1), the inverse Gaussian of mean 1/γ and shape 1 | N(m, s²) drawn, or held (0 by default) | N(m, s²) on γ > 0 drawn, or held (1 by default) | `BassNoise.NormalWald(gamma, skew, c, reference)` |
| Skewed t | IG(ν/2, ν/2), t's | N(m, s²) drawn, or held away from 0 | none | `BassNoise.SkewedT(ν, skew, c, reference)` |
| Generalized hyperbolic (`gbass` with a GIG `v_prior`) | GIG(p, a, b) as declared | drawn or held (0 by default) | none | `BassNoise.GeneralizedHyperbolic(mixing, skew, c, reference)` |

The defaults hold β at 0 and γ at 1, as GBASS's `nwbass` does; a fit that is to choose them states their priors. `BassNoisePrior(mean, sd)` is the normal, sd = 0 holding the scalar.

**σ² and the sweep.** σ² sits in the offset σβv wherever the skew can be nonzero, drawn or held away from 0, so it is held on the record and drawn by the modified half-normal (slice 2's block); where β is held at 0 the noise is a pure scale mixture and σ² is integrated out of the structure step, as under t. After σ² the noise-shape block (`BassNoiseShapeBlock`, the scale slot) draws, each where it is drawn:
- *β | ·*, normal: rᵢ/σ = βvᵢ + √(c·vᵢ)ζᵢ with rᵢ = yᵢ − fᵢ of the sweep's coefficients, so precision Σvᵢ/c + 1/s² and mean (Σrᵢ/(σc) + m/s²) over it, GBASS's β step;
- *γ | v*, normal on γ > 0: GIG(−½, γ², 1)'s normalizer is e^γ/√(2π), so Πp(vᵢ | γ) ∝ exp(Nγ − γ²Σvᵢ/2), precision Σvᵢ + 1/s² and mean (N + m/s²) over it, drawn on its support by the normal's truncated draw. GBASS takes that normal on the whole line, which is γ's conditional only on γ > 0 (slice 0's finding, confirmed there at 0.132 of its mass below 0); this is the correction the brief records, not a parity variant.

Then the chain follows the offset at the new β. v's law is read at the record's γ, so v's conditional is GIG(p − ½, a + β²/c, b + r²/(σ²c)) under the Normal-Wald at the current γ. The target reads neither β nor γ, and neither prior enters it: no move of the structure changes them. The skewed t is t mixing with a drift: its conditional is the GIG with a > 0 and its marginal the generalized hyperbolic's a = 0 limit (`StudentTMixing` routes a drift to `GigMixing`), so at m = 0 it is slice 1's t, draw for draw.

**The readers.** The pointwise density is the generalized hyperbolic's, closed in log Bessel K ([its forms](./ProbabilityAndInformation.md#log-bessel-k-the-gigs-normalizer-and-the-generalized-hyperbolic-e1)); the mean's shift σβE[v] and the variance σ²(β²Var[v] + c·E[v]) are the GIG's moments; predictive draws take v from its law. The predictive quantiles read each draw's noise as the `GeneralizedHyperbolic` distribution, whose tails come from a table of its normal score built on the first read, and the predictive mixture's root is now found by Newton with the mixture's density as its slope, inside the bracket bisection kept, so a CDF that is an integral is read a few times per root rather than once per halving (t's and the asymmetric Laplace's quantiles move by rounding).

**Calibrations: GBASS's `nw_est_mom` and `nw_gamma_prior`, corrected** (`NormalWaldCalibration`).
- *The NIG by its moments.* GBASS fits mean, variance, skewness and kurtosis by `optim` with the kurtosis 3 + 3(1 + 4β²)/(α²δγ); the NIG's is 3 + 3(1 + 4β²/α²)/(δγ) (slice 0 confirmed the error against the mixture's own fourth moment). With it the fit is closed: ρ = β/α has ρ² = S²/(3K − 4S²) for skewness S and excess kurtosis K, δγ = 3(1 + 4ρ²)/K, γ² = δγ/((1 − ρ²)Var), and δ, α, β and μ follow. Moments with 3K ≤ 4S² are no NIG's and are refused. In the core's terms σ = δ, β′ = βδ and γ′ = γδ (`NormalInverseGaussian.CoreSkew`, `CoreGamma`).
- *The prior on γ.* GBASS calibrates a normal on γ to the steepness ξ = (1 + γ)^(−½) of the NIG's shape triangle, P(ξ < q₁) = p₁ and P(ξ > q₂) = p₂, with a stray −1 in the second probability, and its default targets (0.1, 0.9, 0.5, 0.05) ask for half the mass above γ = 99 and 5% below 0.235, which no normal holds. Here the prior is the normal's part on γ > 0, the one γ is drawn from, and the targets are solved exactly. A bound says which targets a normal reaches: at P(γ > g₁) = ½ the most it puts below g₂ is 1 − 2^(−g₂/g₁), the exponential limit as the mean goes to −∞; GBASS's defaults ask 0.05 where that is 0.0016. Targets past it are refused with the most that is reached. The defaults here, (q₁, q₂, p₁, p₂) = (0.5, 0.9, 0.5, 0.03), put half the mass at ξ < ½ (γ > 3) and 3% at ξ > 0.9 (γ < 0.235, the heaviest tails), inside the 0.053 their bound allows.

| Prior on γ | Targets it reaches | γ's draw | Trade-off |
|---|---|---|---|
| **Normal on γ > 0 (chosen)** | Those under 1 − 2^(−g₂/g₁) at p₁ = ½; not a far median with much mass near 0 | Exact: the truncated normal | GBASS's family and conjugacy kept; its reach is limited, and stated |
| Normal on ln γ | Any, calibrated in closed form from two normal quantiles | Not conjugate: a slice step on ln γ, cheap (it reads N and Σv) | A prior spanning heavy to near-Gaussian tails over decades, at an MCMC step for an exact one |
| Gamma on γ | A far median with mass near 0, by a small shape | Conjugate through the modified half-normal, MHN(k, Σv/2, N − 1/θ) | Its shape below 1 puts the draw on the half-normal's series branch, whose terms grow with N |

**`lag_beta`, read** (decision 3: adopted only if a free skew is measured to soak up misfit early). GBASS holds β at its prior's mean for the first 1,000 iterations. Read on NIG data (N = 600, σ 0.3, β 1, γ 1.5), eight chains each way from the empty structure: free from the start, or held at 0 for 500 sweeps and then freed. Free from the start, the median chain reached the reference run's median M (5) in 23 sweeps (6 to 49); held for 500 sweeps, in 26 (13 to 56). Free, β over the first 200 sweeps averaged 1.68 against the reference's 1.21: the skew moves early, while M is small, without slowing the structure. Over sweeps 1,500 to 2,000 the two arms read M 5.06 and 5.16, β 1.64 and 1.40, σ² 0.127 and 0.121, against the reference's 5, 1.21 and 0.105; both arms' β sit above it by less than one of its posterior sd's (about 0.6), over windows short beside its mixing. By the rule stated beside it, adopt the hold only if the free arm takes half as long again to reach M or ends elsewhere, the hold is not adopted.

**Checked** (`BassNormalWaldTests`, `GigMixingTests`, `NormalWaldCalibrationTests`; behind `BASS_QUALIFICATION`, `BassNormalWaldQualificationTests`; against R, `GbassNormalWaldReadoutTests`, `GbassPrimitiveParityTests`):
- *The law.* GIG mixing's marginal against quadrature of the mixture itself over 600 cases, the Normal-Wald, the skewed t, the variance gamma and general GIGs, with drifts of both signs and none, to 10⁻⁸; the skewed t falls out of t mixing and meets the closed t at no drift; the conditional's draws carry its GIG's mean. The noise's tails, CDF and table as [the distribution's checks](./ProbabilityAndInformation.md#log-bessel-k-the-gigs-normalizer-and-the-generalized-hyperbolic-e1) state.
- *Geweke's joint-distribution test*, the slice 2 design at twice its size (400,000 marginal-conditional draws against 600,000 sweeps): five cases, the Normal-Wald with β and γ drawn on two seeds, with γ alone drawn (σ² integrated), the skewed t with β drawn and the generalized hyperbolic GIG(0.7, 1.2, 0.8) with β drawn; every |z| ≤ 2.15 over ten functions each (M, ln τ, ln σ², ln v₀, mean ln v, f at two points, β, ln γ and the standardized residual), and no block draw refused. Its power: β's conditional with its data term's sign flipped fails every skew-drawn case (|z| 20.8 to 336) and leaves the symmetric one passing. At slice 2's size (200,000 against 300,000) one function of one case reached |z| 4.08 (f₀, the Normal-Wald at seed 20261106); the same seed at twice the size reads −0.31, so the gate runs at twice the size.
- *Recovery* on NIG noise outside any fixture (N = 600, f = sin 2πx₁ + x₂, σ 0.3, β 1, γ 1.5; β ~ N(0, 2²), γ under the default calibrated prior; four chains): β 1.29 [0.57, 3.80], γ 1.82 [0.93, 3.89] and σ² 0.113 [0.062, 0.226] at their 99% intervals, each covering the truth; the surface f within RMSE 0.056 of the truth on 400 new points, and the mean f + σβE[v] within 0.029 of the true mean (its shift 0.249 against 0.200). β, γ and σ² trade along the noise's variance, which the data fix more tightly than any one of them.
- *`nwbass`'s readouts* (N = 200, both sides β ~ N(0, 1), γ under N(2, 1), w's prior IG(0.1, 0.05), R with `lag_beta` = 0): the noise's variance, which the data fix, 0.1814 under the g-prior and 0.1809 under the ridge against `nwbass`'s 0.1821; w 0.189 [0.085, 0.370] against 0.155 [0.097, 0.244], β 0.205 [−0.057, 0.613] against 0.191 [−0.033, 0.500], γ 1.13 [0.43, 2.30] against 0.92 [0.53, 1.48], each mean inside the other's 90% interval; the surface's means within 0.55 of the wider 90% half-band at worst (0.76 under the ridge); M 4.32 against 4.44. w and γ trade along the noise's variance, so each is read for overlap. The differences stated beside it: GBASS's untruncated γ draw, and slice 2's (continuous knots, whole-term mutation, pseudo-counts 1 against 5, its raw ridge, λ sampled, priors reading no response).
- *The calibrations.* The NIG by its moments recovers three NIGs from their moments taken by quadrature of their own densities to 7 digits, and α, β and δ from 400,000 draws; the prior on γ puts its two masses where asked to 8 digits on three targets, read by quadrature of the truncated normal, and refuses GBASS's defaults with the bound's value.
- *The run.* A Normal-Wald fit resumes bit for bit, β, γ and σ² included; a retained draw keeps its latents, γ with them, and its view rebuilds bit for bit under the Normal-Wald, the skewed t and the generalized hyperbolic.

**Cost** (`CostRow_PerBlockOfASweep`; thread cycles on a loaded machine, so indicative; µs per sweep by block, the Normal-Wald with β and γ drawn, under the g-form).

| Point | M | Structure | τ | a | v's block | σ² | β and γ | Sweep | Gaussian sweep |
|---|---|---|---|---|---|---|---|---|---|
| n 1,000, p 10, the base | 20 | 193 | 3.6 | 31 | 987 | 31 | 26 | 1,271 | 28 |
| n 200 | 9.4 | 33 | 2.5 | 5.0 | 181 | 6.6 | 6.4 | 235 | 4.2 |
| n 10,000, M_max 200 | 199 | 80,634 | 18 | 4,560 | 229,628 | 2,456 | 2,455 | 319,750 | 1,638 |

- *A Normal-Wald sweep is 45 Gaussian sweeps at the base, 56 at n 200 and 195 at M 200.* A second reading, beside slices 1's and 2's gates on four threads, read the Normal-Wald's blocks about twice as long and the Gaussian sweep the same: the machine's contention, so the solo reading is the table.

- *v's N GIG draws and the Follow after them are most of a sweep,* as under the asymmetric Laplace; the skew and γ cost a block of O(N) sums, near σ²'s.
- *Retention* is slice 2's: a retained draw keeps v with its scalars, γ now among them.

## The generalized beta prime and the horseshoe (E1's slice 4)

Built 2026-11-06 ([the build brief](../../issues/regression/bass/gbass-noise-brief.md#reports), slice 4), GBASS's `hbass` on the BASS core and its general case, `gbass` with a GBP `v_prior`: v ~ GBP(p, a, b), density ∝ v^(pa−1)·(1 + v^p)^(−(a+b)) ([the law](./ProbabilityAndInformation.md#the-generalized-beta-prime-and-its-normal-mixture-e1)), heavy at both ends of v, so that an observation is fitted closely, v small, or let go, v large: "interpolate or ignore".

| Kind | v's law | β | On the config |
|---|---|---|---|
| Horseshoe (`hbass`) | GBP(1, ½, ½), √v ~ C⁺(0, 1) | N(m, s²) drawn, or held (0 by default) | `BassNoise.Horseshoe(skew, c, reference)` |
| Strawderman–Berger (`hbass`, `likelihood = "sb"`) | GBP(½, 1, ½) | as above | `BassNoise.StrawdermanBerger(skew, c, reference)` |
| Generalized beta prime (`gbass` with a GBP `v_prior`) | GBP(p, a, b) as declared | as above | `BassNoise.GeneralizedBetaPrime(mixing, skew, c, reference)` |

**`hbass` as configuration.** `hbass` takes c = var(y) and w ~ GBP(1, ½, ½), √w ~ C⁺(0, 1), so the noise's scale per unit v, √(w·c), is half-Cauchy at sd(y). Here c = 1 and σ² under `Sigma2Prior.HalfCauchy()`: σ ~ C⁺(0, s) at the reference's mean square s², by the hierarchical route σ² | r ~ IG(½, r), r ~ Gamma(½, rate 1/s²) (Wand et al. 2011), with r's slice step (`BassSigma2RateBlock`). It is the same prior on the same quantity, the engine's σ² being `hbass`'s w·var(y), up to the reference's divisor N against var's N − 1; GBASS steps w by a log random walk. Nothing was built for it: the hierarchical prior was the core's and the noise never refused it. Slice 4 runs it under the noise: σ² integrated at the rate's prior where β is held at 0, and on the record, read by the rate's density (`BassTarget.Evaluate`'s σ² term) and drawn against IG(½, r) at the state's rate, under a skew. Its scale is read from y, so Geweke's test, whose y moves, cannot carry it; the readouts below can.

**v's step.** Under GBP mixing v's conditional, ∝ p(v)·v^(−½)·exp(−(r − m·v)²/(2s²v)), has no exact draw, so `IMixingLaw` now steps: `StepConditional(current, residual, s², m, b)`, which an exact law answers with its draw, the current value unread, so t, the asymmetric Laplace and the GIG draw the numbers they drew; the latent block passes each vᵢ's current value. `GbpMixing` takes one slice step on u = ln v (Neal 2003; `SliceSampler`, width 2, at most 64 widths), whose log density (pa − ½)·u − (a + b)·ln(1 + e^(pu)) − r²/(2s²)·e⁻ᵘ − m²/(2s²)·eᵘ is concave, so the slice is one interval. Its two exponential terms are taken in logs, so a residual whose square is below the double range keeps its wall. At a residual of exactly 0 the conditional is improper where pa ≤ ½, the horseshoe's and the Strawderman–Berger's case, the noise's pole: refused, an event of probability 0 under continuous data. `SliceSampler` gained a struct overload, so N steps a sweep allocate no delegates; the delegate form calls it, its numbers unchanged.

| v's step | Exact | Tuning | Trade-off |
|---|---|---|---|
| **Slice on ln v (chosen)** | A step that leaves the conditional invariant | None: stepping out and shrinkage size the interval | 6.6 to 6.8 evaluations of a few exponentials a step, 0.41 µs; one step a sweep |
| GBASS's random walk on ln v (`rgbp`) | A step, invariant at a fixed size | A size adapted in burn-in from mean(ln v)'s spread, which slice 0 found 0.543 against the per-observation 3.96 | Its mixing rides that tuning |
| Auxiliary gamma: at p = 1, v = G_a/G_b, so v \| λ ~ Gamma(a, λ), λ ~ Gamma(b, 1) | Exact draws: v \| λ, r is GIG(a − ½, 2λ + m²/s², r²/s²), λ \| v Gamma(a + b, 1 + v) (Armagan, Dunson and Clyde 2011) | None | p = 1 only, not the Strawderman–Berger; a second latent per observation; Gibbs between v and λ slows where they correlate |
| Rejection from the prior under the likelihood's peak | Exact | None | Its acceptance falls as \|r\|/s near the pole, under 1% at r = 10⁻³·s |

**The readers.** The noise's density has no closed form under a drift or a general GBP (the horseshoe's at m = 0 is (2π³)^(−½)·e^c·E₁(c)/s with c = r²/(2s²), which the checks read): it is the integral over ln v (`NormalGeneralizedBetaPrime`), hundreds of nodes a read. A pointwise density reads it N × draws times, so it is opt-in for that path: `BassFit.PointwiseLogDensity(data, quadrature: true)` and `BassPointwise.LogDensity(…, quadrature: true)` take it, and the conformal layer's candidate density the same way (`BassFit.CandidateLogDensity(…, quadrature: true)`, draws × candidates reads); without the flag a noise whose `IBassNoiseLaw.DensityByQuadrature` holds is refused with its reason. The predictive's quantiles read each draw's noise as that distribution, its tails from the mixture's table. The horseshoe and the Strawderman–Berger have neither a mean nor a variance: the variance reads +∞, so the training fit states no noise sd; under a drawn skew the mean's shift σβE[v] is ±∞, a noise without a mean. The pointwise density's latent report now calls the per-observation latent integrated (`LatentTreatment.Integrated`) under every noise law: given the globals nothing outside its observation moves it, so the term is in-sample, and WAIC and PSIS-LOO read it.

| The pointwise density | A read | Trade-off |
|---|---|---|
| **The integral over ln v, opt-in (chosen)** | 41 µs | Every GBP and drift; N × draws reads, so asked for |
| The horseshoe's closed form at m = 0 | An exponential integral | The symmetric horseshoe alone; a second path for one case |
| Monte Carlo over v | A few hundred draws | Its noise enters the pointwise checks as the model's |

**The mixture's integral, upgraded.** The generalized hyperbolic's integral and table are now a base the GBP's mixture shares (`NormalVarianceMeanMixture`), and three changes the GBP's heavy tails and pole forced move the generalized hyperbolic's tables too ([the forms](./ProbabilityAndInformation.md#the-generalized-beta-prime-and-its-normal-mixture-e1)): the integral about its peak where x and the drift share a sign, without the cancellation that made the far drift side's integrand noise and a table's build run for minutes; the table's cubic in asinh(x/s), where a power tail's score is smooth and the horseshoe's table in x missed its integral by up to half the log tail at 10⁻⁹; and the table's steps into a pole, at most ¼ of |x|. Its draws do not move; its predictive quantiles move at the table's level, under 3·10⁻⁴ of a log tail.

**Checked** (`BassGeneralizedBetaPrimeTests`, `GbpMixingTests`, `GeneralizedBetaPrimeTests`; behind `BASS_QUALIFICATION`, `BassGeneralizedBetaPrimeQualificationTests`; against R, `GbassHorseshoeReadoutTests`):
- *The law.* The slice step from exact draws of v's conditional (by rejection from the prior) leaves them on it, Kolmogorov–Smirnov against the conditional's CDF by quadrature in ln v, √n·D at most 1.17 against the cut 2.6 on eight cases: the horseshoe at residuals 10⁻³, ½ and 3, with a drift, at power ½ and at the prior, the Strawderman–Berger with a drift, a general GBP; every draw moved. A chain of 400,000 steps at r = 0.01 carries the conditional's mean of ln v (−4.582 ± 0.005 against −4.575). The exact laws' steps equal their draws, number for number. The noise's density, tails and table as [the distribution's checks](./ProbabilityAndInformation.md#the-generalized-beta-prime-and-its-normal-mixture-e1) state.
- *Geweke's joint-distribution test* at slice 3's size (400,000 marginal-conditional draws against 600,000 sweeps), v moved by its slice step: the horseshoe with σ² integrated on two seeds, with β drawn and σ² on the record, and the Strawderman–Berger; every |z| under the bar of 4 over nine functions each (M, ln τ, ln σ², ln v₀, mean ln v, f at two points, β and the standardized residual), the worst 3.67, ln τ under the horseshoe at seed 20261106, which reads −0.19 at twice the size. Three more seeds at that size read |z| ≤ 1.91, and four times the size on a fifth, the horseshoe and its skewed form, ≤ 2.64; the skewed form's v functions, to 2.64, change sign between seeds, slow mixing along the skew and scale rather than a bias. Of 600,000 latent draws 0 to 9 were refused, each outside F and kept as the gate's exact step. `hbass`'s half-Cauchy cannot enter: its scale reads y, which the successive simulator moves.
- *Recovery* from a contaminated normal outside the family (0.85·N(0, 0.25²) + 0.15·N(0, 2²), f = sin 2πx₀ + x₁, N = 300, six data sets), the horseshoe under `hbass`'s half-Cauchy against Gaussian BASS: f's 90% intervals cover 80.6% of 294 grid points against Gaussian BASS's 86.4%, at an RMSE of 0.086 against 0.142; the 90% predictive intervals cover 92.9% of 1,800 held-out responses against 92.7%. The horseshoe lets the outliers go and fits the surface closer than Gaussian BASS does. t(3) on the same data sets covers 89.5% at an RMSE of 0.059 (slice 1): its v is kept from 0 by its prior, where the horseshoe's v near 0 lets a fit follow some observations closely, and the horseshoe's intervals for f undercover.
- *`hbass`'s readouts* (`GbassHorseshoeReadoutTests`, N = 200, R 6,000 iterations, the engine 4 chains × 1,500 draws). `hbass` as written adapts v's step from mean(ln v)'s spread, slice 0's defect, to 0.12 at an acceptance of 0.98 here (0.038 at 0.99 on a second set of N = 200), so v mixes slowly: against the same sampler with v's step held at 2.5 and its adaptation off (acceptance 0.55), the correction and the reference, it reads M 4.16 against 5.36 and w·var(y) 0.082 against 0.068. Against the reference, at engine seed 41: the surface within 0.55 of the wider 90% half-band under the g-prior and 0.82 under the ridge; σ² 0.0765 and 0.0802 against w·var(y) 0.0685 [0.043, 0.103]; mean ln v −0.177 and −0.192 against −0.162; v below 0.1 in 16.7% of draws under both against 17.2%; M 5.84 and 5.10 against 5.36. Over engine seeds 41 to 46 the ridge read 0.52 to 0.82 and the g-prior 0.55 to 1.20, so the bars are `qbass`'s: parity at 1.0 for the ridge, GBASS's own form, and 1.5 for the g-prior, a stated model difference. Against `hbass` as written the g-prior read 0.97 to 1.56: the gap was its sampler, not the engine.
- *The run.* A horseshoe fit with a drawn skew under the half-Cauchy resumes bit for bit, v, σ², β and σ²'s rate included; a retained draw's view rebuilds bit for bit under the horseshoe, its skewed form and the Strawderman–Berger; the sweep runs the rate's block with the noise's, under σ² integrated and on the record.

**Cost** (`CostRow_PerPartOfASweep_AndThePointwiseRead`; thread cycles on a loaded machine, so indicative; µs per sweep by part, the horseshoe under the g-form, t(4)'s exact draws of v at the same residuals beside).

| Point | M | Structure | τ | a | v's slice steps | Follow | Sweep | t(4)'s v draws | Gaussian sweep |
|---|---|---|---|---|---|---|---|---|---|
| n 1,000, p 10, the base | 14 | 106 | 3.2 | 24 | 413 | 221 | 767 | 74 | 49 |
| n 200 | 3.1 | 11 | 2.0 | 2.6 | 83 | 12 | 110 | 9.9 | 4.2 |
| n 10,000, M_max 200 | 200 | 86,295 | 18 | 3,995 | 3,807 | 225,347 | 319,461 | 471 | 1,583 |

- *v's slice steps* take 6.6 to 6.8 evaluations each, 0.41 µs at the base, 5.6 times t's exact draw: half the base's sweep, under the Follow at M 200.
- *A horseshoe sweep is 16 Gaussian sweeps at the base,* 26 at n 200 and 200 at M 200.
- *The pointwise read* (wall clock) is 41 µs against t's closed form's 0.11: at N 1,000 and 2,000 draws, 82 s, the opt-in's reason.
- *The predictive's quantiles* read a table per draw, 15 to 28 ms for the horseshoe and the Strawderman–Berger (wall clock; 9 ms for a law without a pole): the recovery's twelve predictions over 2,000 draws each took most of its 29 minutes.

## The GIG's draw: ratio-of-uniforms where it is faster

Built 2026-11-06 after slice 3, at the user's word that spline surfaces are expensive and every bit of performance counts. v's conditional is a GIG under the asymmetric Laplace, the Normal-Wald, the skewed t and the generalized hyperbolic, and its N draws were most of a sweep. An untruncated GIG draw is now Hörmann and Leydold's (2014), GIGrvg's ([the method](./ProbabilityAndInformation.md#the-gig-the-modified-half-normal-and-the-asymmetric-laplace-e1)), on a ratio-of-uniforms primitive other families can take (`RatioOfUniforms`); a truncated draw stays the tangent hat's. Both draw the same law exactly and independently, so a chain mixes as it did; the measure is time a draw.

**Measured** (thread cycles on a loaded machine, so indicative; µs a draw, each with its own GIG's setup as the latent block draws, the least of five passes after six warm-up passes past tier 0). v's GIG parameters recorded from burned quantile (q = 0.9) and Normal-Wald chains at the cost row's points, in the standard form's λ = |p| and ω = √(ab):

| Fit and point | λ | ω at the 10th, 50th, 90th percentile | Ratio-of-uniforms | Tangent hat | Ratio faster |
|---|---|---|---|---|---|
| Quantile, n 1,000, p 10, the base | ½ | 0.59, 3.7, 8.0 | 0.171 | 0.950 | 5.6× |
| Quantile, n 200 | ½ | 0.48, 3.6, 8.5 | 0.126 | 0.703 | 5.6× |
| Quantile, n 10,000, M_max 200 | ½ | 0.58, 3.3, 7.2 | 0.123 | 0.637 | 5.2× |
| Normal-Wald, the base | 1 | 6.3, 6.6, 7.7 | 0.127 | 0.568 | 4.5× |
| Normal-Wald, n 200 | 1 | 3.0, 3.5, 4.4 | 0.206 | 0.518 | 2.5× |
| Normal-Wald, n 10,000, M_max 200 | 1 | 4.3, 4.5, 5.5 | 0.126 | 0.579 | 4.6× |

On far ranges, λ from 0 to 1,000 by eleven steps and ω from 10⁻⁸ to 10⁸ by twelve, 132 points, ratio-of-uniforms took 0.08 to 0.36 µs a draw and the hat 0.28 to 0.92, the ratio the faster at every point.

**Which one runs** (`GeneralizedInverseGaussian.DrawsByRatio`, a function of (p, a, b) alone, so a resumed run draws as the whole one did):

| Draw | Sampler | Why |
|---|---|---|
| Untruncated, a, b > 0, \|p\| ≤ 10³ and 10⁻⁸ ≤ ω ≤ 10⁸ | Ratio-of-uniforms | Faster at every point measured; its law checked over this range (18 chosen cases across its three regions and 300 seeded far-range cases, against quadrature; R's `rgig` law in the oracle) |
| Untruncated beyond that range | Tangent hat | One method for every parameter; the ratio's rectangle, from a cubic's roots by Cardano's rule, is not checked past it |
| Truncated to an interval | Tangent hat | The ratio's rectangle is the whole law's; cut to a far tail it would reject nearly everything, where the hat cut to the interval tends to accept all |
| b = 0 or a = 0 | The gamma's or the inverse gamma's own draw | The limits, as before |

**The trade-off.** The ratio is three methods by region, with their rectangles in closed form, so more code and a range checked rather than every parameter; the hat is one method for every parameter and interval, at the cost of a Newton search for its design points every draw, which is where its time goes. Within the range the ratio is the default; beyond it the hat answers.

**What it changed.** On the quantile fit at the base point (n 1,000, p 10, the g-form) v's GIG draws went from 952 µs a sweep to 208, from 53% of the sweep to 24%, and the sweep from 1,810 µs to 850; on the Normal-Wald, v's block (its draws and the Follow after them) went from 987 to 681 µs and the sweep from 1,271 to 1,027 at the base, from 235 to 133 at n 200. At M_max 200 the Follow after v dominates and the sweep is as it was. (Thread cycles, each row run alone, the before and after on different days: indicative.) Every quantile, Normal-Wald, skewed-t and generalized-hyperbolic fit draws v anew, so their draws and their readings moved within their gates (the change's report gives the numbers); no pinned hash covers those fits. Gaussian BASS draws no GIG and t without a skew draws v by its own gamma step, so neither moved. The envelope's pin moved by exactly its 200 untruncated GIG cases.

## Retention: a draw keeps its latents

Read in slice 1, built in slice 2: the brief's lever for a working response that dominates a retained draw. A model that rebuilds its working response from its latents (`IBassObservationModel.Rebuild`, `VarianceMeanMixture`'s) has its retained draws keep `BassDraw.Latents`, v with β and σ², and rebuild the working response on demand (`WorkingOn`); the readers that need the latents alone take them without a rebuild (`LatentsOn`). The blocks and the rebuild build the working response by one function from the latents, so a view rebuilt from v is the one the chain held, bit for bit. The generative state stays: v is kept. Measured over 200 retained draws of one quantile chain (`BassQuantileQualificationTests.Retention_BeforeAndAfter`), held both ways: 313.4 KB a draw with its working response against 70.4 KB with its latents alone at N = 10,000 (4.45×; the weights, the offset, the shifted response and v against v, 4N doubles against N); 32.7 KB against 9.0 KB at N = 1,000 (3.6×), where the draw's own few hundred bytes and the collector's noise weigh more. Under t the working response held the weights and v, so the ratio there is 2.
- *The alternatives:* keeping the working response, which holds the weights, the offset and the shifted response beside v, up to 4N doubles under the skew; or keeping nothing and redrawing, which loses the generative state.
- *The cost:* a rebuild is O(N) per draw read, for the views that read the offset or the weights (the residual check, the training fit's offset, the conditional means).

## Status and caveats

- *Built:* slices 0 to 4 of [the brief](../../issues/regression/bass/gbass-noise-brief.md), each with its report there.
- *Open, by slice:* the latents on tempered rungs, SMC and the flat histogram, and a per-observation latent under a direct functional design (5); BARS through the shared law (6); the z parameterization and interweaving (7).
- *The asymmetric Laplace is a working likelihood* where the noise is not one: its posterior intervals for the q-quantile undercover under other noise (63.8% for nominal 90% at q = 0.9 on skewed, heteroscedastic noise), as Yang, Wang and He (2016) show in general. The surface is the q-quantile; its intervals are not yet calibrated for it.
- *The Normal-Wald's prior on γ* is a normal on γ > 0, which cannot hold a far median with much mass near 0 (the bound above); a fit meant to range from near-Gaussian to the heaviest tails needs a prior on ln γ, open.
- *A mixture's predictive tails* come from a table built per draw: the generalized hyperbolic's some 60 tail integrals, about 10 ms; the horseshoe's 15 to 28 ms, its pole and heavier tails taking more points; a predictive over thousands of draws spends seconds to minutes building them.
- *The horseshoe's pointwise density* is an integral per observation and draw, 41 µs each, so it is opt-in; its symmetric case has a closed form (the exponential integral) not taken.
- *`hbass`'s half-Cauchy has no Geweke check under the noise:* its scale reads y, which Geweke's successive simulator moves; the readouts and the run carry it.

## Source map

| Area | Files |
|---|---|
| The noise on the config, the model, its blocks | `src/maths/regression/spline/bass/BassNoise.cs` (`BassNoise`, `BassNoisePrior`, `ObservationLatents`, `IBassNoiseLaw`, `VarianceMeanMixture`, `BassCoefficientsBlock`, `BassLatentScaleBlock`, `BassScaleBlock`, `BassNoiseShapeBlock`) |
| The carrier and the seam | `BassObservationModel.cs` (`WorkingResponse.Latents`, `IBassObservationModel.Law`, `Rebuilds`, `Rebuild`) |
| σ² read from the state | `BassTarget.cs` (`GivenAt`) |
| The sweep's draw of a and σ² | `BassSweep.cs` (`BassStructureMoves.DrawSweepFit`, `SweepFit`) |
| The draw keeping its latents | `BassState.cs` (`BassDraw.Latents`, `WorkingOn`, `LatentsOn`) |
| The prior's energy for σ²'s block | `src/maths/regression/gprior/CollapsedGPrior.cs` (`GPriorConditional.PriorEnergy`), `PriorPencil.cs` |
| The mixing laws | `src/maths/regression/MixingLaw.cs` (`IMixingLaw`, `StudentTMixing`, `ExponentialMixing`, `GigMixing`, `GbpParameters`, `GbpMixing`) |
| v's slice step | `src/maths/samplers/blocks/SliceSampler.cs` (`IScalarLogDensity`) |
| The pointwise opt-in | `BassPointwise.cs` (`LogDensity`'s `quadrature`, `Latents`), `BassFit.cs` (`PointwiseLogDensity`) |
| The Normal-Wald's calibrations | `src/maths/regression/NormalWaldCalibration.cs` (`NormalWaldCalibration`, `NormalInverseGaussian`) |
| The predictive's quantile under a law | `src/maths/regression/spline/bass/posterior/PredictiveMixture.cs` |
| The laws drawn from and read | `src/maths/distributions/GeneralizedInverseGaussian.cs`, `ModifiedHalfNormal.cs`, `AsymmetricLaplace.cs`, `GeneralizedHyperbolic.cs`, `BesselK.cs`, `RatioOfUniforms.cs`, `NormalVarianceMeanMixture.cs`, `GeneralizedBetaPrime.cs`, `NormalGeneralizedBetaPrime.cs`, `Beta.cs` (`InverseLogit`) |
| Tests | `tests/regression/BassStudentTTests.cs`, `BassStudentTQualificationTests.cs`, `BassQuantileTests.cs`, `BassQuantileQualificationTests.cs`, `BassNormalWaldTests.cs`, `BassNormalWaldQualificationTests.cs`, `GigMixingTests.cs`, `NormalWaldCalibrationTests.cs`, `BassGeneralizedBetaPrimeTests.cs`, `BassGeneralizedBetaPrimeQualificationTests.cs`, `GbpMixingTests.cs`; `tests/distributions/BesselKTests.cs`, `GigNormalizerTests.cs`, `LogConcaveEnvelopeTests.cs`, `GeneralizedBetaPrimeTests.cs`; `tests/oracle/Gbass*.cs` |
