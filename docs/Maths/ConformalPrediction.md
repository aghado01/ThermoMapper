# Conformal prediction

> The conformal layer (E16): prediction sets whose coverage holds in finite samples under exchangeability, whatever the model, from an engine's retained draws and pointwise densities with no refit. The layer owns the guarantee and the score owns the width, so the layer is built once, engine-free, and each engine contributes densities. Built: the full step (Fong and Holmes's conformal Bayes) over a training density and a candidate density; BASS's candidate density for every noise law; and CB-3's reading of the scores' Monte Carlo error, with the margin, the inflated rule and the fallback signal. The split form, the predictive distribution and the other engines' candidate seams follow.

## Overview

Conformal prediction ranks a candidate observation among the data by a conformity score that treats every unit alike. Under exchangeability the candidate's rank is uniform, so a set that refuses only the lowest ranks covers a new observation at a level fixed in advance, whether or not the model is right. The model only decides which units score high, and so how wide the set is.

What it serves, by problem class:
- prediction sets for a new observation that hold whether or not the noise model is right;
- coverage within strata, for partially exchangeable data (Mondrian conformal prediction, Vovk; grouped conformal Bayes, Fong and Holmes);
- later, a calibrated predictive distribution at each x, and a guarantee that does not lean on the sampler (the split form).

It does not reach a curve's or a parameter's own intervals: those stay with [calibrated intervals and resampling](../../issues/regression/bass/extensions-brief.md#calibrated-intervals-for-loss-based-posteriors).

Where things are:
- *The contract*, CB-1 to CB-4, is proved in `lean/Enthymemata/ConformalValidity.lean`, with its statements, checks and decisions in [conformal-lemmas.md](../../lean/Protolemmata/Conformal/conformal-lemmas.md). Every contract test below names its declaration.
- *The build* is [the layer's brief](../../issues/regression/conformal-layer-brief.md), slice by slice, and E16 in the [extensions brief](../../issues/regression/bass/extensions-brief.md#conformal-bayes-the-consumer) holds the literature.
- *The seams it reads* are [BayesianRegression](./BayesianRegression.md)'s pointwise density (OQ17), the draws × units matrix of log p(y_unit | draw), and its sibling at new cases, the candidate density.
- *The product it runs on* is [LinearAlgebra](./LinearAlgebra.md#dense-product-denseproduct)'s `DenseProduct`.

## The seams

**The training density** (`PointwiseLogDensity`): log p(y_unit | θ_s) for each retained draw s and each training unit, with the unit partition and its stratum map, how each latent was treated (`LatentReport`), the chain lengths, and draw weights for a weighted population.

**The candidate density** (`CandidateLogDensity`): log p(y | x_t, θ_s) for the same draws in the same rows, at each test case t and each value y of t's grid. The columns run case by case, each case's grid ascending; each case carries a stratum. Its latents, chain lengths and weights must be the training density's. A value outside a draw's support takes −∞. An engine joins the layer by producing it; the layer knows no engine.

## The full step (Fong and Holmes 2021)

For a candidate y at case t:
- **The add-one-in weights.** w̃_s ∝ W_s·p(y | x_t, θ_s), with W_s the draw's weight (equal for chains). The draws, from the posterior given the n training units, are reweighted to the posterior given all n + 1.
- **The scores.** Each unit's conformity is σᵢ = Σ_s w̃_s p(zᵢ | θ_s), the self-normalized estimate of its in-sample predictive p(zᵢ | all n + 1 units), and the candidate's own is σ_c = Σ_s w̃_s p(y | x_t, θ_s).
- **The rank.** Within the candidate's stratum of m units (its n_h training units and the candidate), π(y) = (#{i : σᵢ ≤ σ_c} + 1)/m, ties counted for the candidate.
- **The set.** At miscoverage α, k = ⌊α·m⌋, and y is accepted when π(y) > k/m, equivalently π(y) > α. The set covers at least 1 − k/m ≥ 1 − α, and exactly 1 − k/m with distinct scores (CB-1).

The guarantee is exact when the weights are the exact posterior's. From Monte Carlo draws the scores are estimates, and their error costs at most the margin CB-3 bounds (slice 2).

## The choices

Each choice below is decided; its alternatives and what they would change are recorded with it.

**The score's content: the retained parameters, in sample.** The score is Fong and Holmes's, each unit's predictive given all n + 1 units, read from the retained draws.
- *The alternative:* Hoff's (2023) Bayes-optimal score, each unit's predictive given every other unit, the candidate included. It ranks differently: with no covariate, y = (2, 2, 2, 0) puts the lone 0 lowest in sample and highest by Hoff's (CB-4's finding).
- *Why not now:* from draws it needs the weight g(y_c | θ)/g(yᵢ | θ) per data unit, add-one-in and leave-one-out at once, with PSIS's instability. Both scores cover; only the width differs. It is a later option.
- *Refused:* the collapsed content (`PointwiseContent.CollapsedPredictive`), whose terms leave each unit out of a rest that excludes the candidate. That is neither score.

**Ties count for the candidate** (CB-1, `coverage_orderings`, `card_rankLe_le`). The rule is valid whatever the ties: coverage at least 1 − k/m, and at most 1 − k/m plus the share of units whose ties straddle the threshold (`coverage_orderings_le`). With every score tied it accepts every candidate.
- *Never the strict rule*, which accepts when at least k units score strictly below: it covers at most 1 − k/m (`coverage_strict_le`) and refuses every candidate when every score ties.
- *Random tie-breaking* (`coverage_tieBreak`) gives exactly 1 − k/m. It waits for the predictive distribution (slice 4), where a uniform p-value is wanted. The layer reports, for each grid value, how many units score strictly below (`Below`) and how many tie (`Tied`), so any rule reads from the step.
- *Ties are exact.* A candidate whose column equals a training unit's scores exactly as that unit does (see the numerics below), so the tie the rule counts is the one the model has, not one rounding made or broke.

**Ranks within the candidate's stratum** (`coverage_law` with the stratum as the ranked set, `coverage_stratum`). The score may read the other strata in any way, through the posterior. A stratum of m units attains only levels 1 − k/m, and coverage exceeds 1 − α by less than 1/m.
- *Small strata:* below m = ⌈1/α⌉, k is 0 and the set is the whole space. The layer reports that (`WholeSpace`) and computes nothing for the case. A stratum no training unit holds has m = 1.
- *Pooling strata* instead, ranking across all units at a level set for the pooled count, loses the per-stratum guarantee: in the Lean checks, a stratum of 2 beside one of 3 covers 0.472 against 0.6 so ranked.

**Every unit is read by one density** (CB-4). The default, and the only reading built, integrates every observation-level latent (OQ17's scope rule): the training terms and the candidate's are the noise law's marginal density.
- *Refused:* a training matrix whose observation-level latents are conditioned on or mixed, and a candidate whose latents are reported otherwise than the training units'. Scoring the data units at their retained latents and the candidate integrated undercovers: it never covers more than the symmetric latent score (`accepts_latent_of_mixed`), and covers 0 against 1/2 in the Lean counterexample (`mixed_counterexample`), which the layer's test reproduces.
- *Not built:* the latent-posterior score for every unit, CB-4's (c), the data units at their retained latents and the candidate's latent drawn from p(v | θ, y) for each draw and candidate. It is valid (`latentScore_equivariant`) and needs a second candidate matrix, the conditioned density beside the integrated one that sets the weights.

**Sets on a grid in y.** The layer evaluates π at the stated grid and reports the accepted values.
- *A union of intervals.* π need not be monotone in y, so a set can split, as a bimodal predictive's does. `Intervals` gives the maximal runs of accepted grid values.
- *Its resolution.* Each boundary is known to within the grid's largest gap (`Resolution`). `OpenBelow` and `OpenAbove` flag a set that reaches the grid's end and may continue past it. `Width` reads each accepted value as its cell, within one gap per boundary.
- *Not the default:* refining the boundaries between grid values. It is an option for later, at the cost of more candidates where the set changes.

**A candidate outside every draw's support** scores −∞: its add-one-in posterior does not exist. It takes π = 1/m, the limit of a vanishing density, and ESS 0. It is refused whenever k ≥ 1.

**Weighted draws work as they are.** An SMC population's or the flat histogram's weights multiply the add-one-in weights, and draws of weight 0 drop out (`PredictiveChecks`' rows, shared). With exact posterior weights over a finite parameter the step is exact, which is how the CB-1 checks are run.

**The Monte Carlo error is read by CB-3** ([below](#cb-3-the-scores-error-and-what-it-costs)). Each candidate's Kish ESS of its add-one-in weights, 1/Σ w̃² (`SmcWeights.EffectiveSampleSize`'s form), is reported beside it: a low ESS marks a candidate far in a tail, where few draws carry the weights.

## CB-3: the scores' error and what it costs

The scores are Monte Carlo estimates, so the step's guarantee holds up to their error (CB-3, `coverage_close`, `coverage_inflated`). Slice 2 reads that error and what it can cost, and reports it with every set (`ConformalSet.Margin`, a `ConformalMargin`). Each choice below is decided; its alternatives and what they would change are recorded with it.

**The error is each score's Monte Carlo standard error, on the log scale.**
- *The form* is PSIS-LOO's own (`PredictiveChecks`' influence, shared). For a self-normalized estimate p̂ = Σ_s w̃_s p_s, each draw's influence is ψ_s = w̃_s(p_s/p̂ − 1), and the error is √log(1 + Σ_s ψ_s²/r_eff).
- *r_eff* is the relative efficiency of the candidate's density along the chains, the add-one-in ratio's, as PSIS reads the left-out unit's likelihood. It is 1 for independent or weighted draws.
- *For every unit at once:* Σψ² = A₂/p̂² − 2A₁/p̂ + Σw̃², with A₁ = Σ w̃²p and A₂ = Σ w̃²p², so two more products, P·W̃² and P²·W̃², give it beside the scores. An entry below 10⁻¹⁰⁰ of its unit's largest likelihood is summed in logs, since squared likelihoods underflow past 10⁻¹⁵⁴.
- *Alternatives:* the influence's ESS along the chains for each unit and candidate, which is exact for chains but an FFT for every pair; a concentration bound, which the Lean follow-up leaves open.

**The bound is z standard errors, with z = Φ⁻¹(1 − ν/(2m)).** Every one of a stratum's m units then stays within its bound with chance about 1 − ν, by the union of their normal tails. ν (`ErrorLevel`, 0.01 by default) adds to the miscoverage, as CB-3 says a bound that is not sure must.
- *Indicative:* the standard error is itself estimated, and the normal approximation fails where the weights degenerate. That is what k̂ reads (below).
- *Alternatives:* a fixed multiple, 2 or 3, with no stated failure probability and narrower bounds; a concentration bound, which would be sure but wider.

**The margin** (`coverage_close` with `card_fragile_le_estimated`) is always reported. For each grid value taken as the truth, it counts the band of units whose decision estimates within their bounds could change (`Fragile`). m̂ is the largest count over the grid values inside some draw's support, since any of them may be the truth. The plain set's coverage under error is then 1 − (k + m̂)/m − ν (`CoverageUnderError`).
- *Read with the uniform bound by default* (`MarginBounds.Uniform`): the largest of the stratum's units' bounds and the candidate's, under this candidate's weights. That is the safe reading of `coverage_close`'s unit-carried bounds.
- *Or with carried bounds* (`MarginBounds.Carried`): each unit takes the larger of its bound here and its bound as the candidate, under weights of its own density. Every unit is the candidate in some ordering, so this is the Lean chip's per-unit finding. It costs each unit's own weighting and, for chains, its r_eff, once per stratum.
- *Both read the orderings the layer holds.* A bound sure over every ordering would need each unit's error under every other unit's weights, n² weightings per stratum. That is not built.

**The inflated rule** (`coverage_inflated`, Ndiaye's stable set) is an option (`Inflate`, off by default) and always reported (`InflatedRanks`). It accepts y when more than k units have an estimate lowered by its bound at most the candidate's raised by its own. The bounds are per position, each unit's under this candidate's weights and the candidate's own. The rule keeps 1 − k/m − ν at the cost of width.

**Trust** is PSIS's Pareto k̂ of the add-one-in log ratios (`PredictiveChecks`' k̂, shared). It is read where the inflated rule accepts, since there a decision rests on the bounds; above PSIS's threshold at S, the error estimate there cannot be trusted.

**The fallback signal** (`FallBack`) is the signal to take the split form (CB-2, slice 3). It reads the inflated set against the plain one, decided by the user on 2026-11-07.
- *What it reads.* The inflated set holds the plain one, since raising the candidate's estimate and lowering the units' only adds units at or below it. The values only the inflated rule accepts are those whose decision rests on the scores' error.
- *The rule:* it fires when they exceed a share of the plain set, (|I| − |A|)/|A| > share (`Gap`, `FallbackShare`, 0.5 by default). The sets are read on the grid's cells as `Width` reads them, and by count on a grid of one value, where it fires exactly when the decision rests on the error.
- *Or on trust:* it also fires when a candidate in the inflated set has an untrusted k̂.
- *The default,* 0.5, fires when half the plain set's extent again rests on Monte Carlo error. On the readings below it fires on about half the BASS Gaussian fits' grids, a third of the t fits', and almost never on exact draws through the gap.
- *Alternatives.* A lower share fires earlier: 0.25 caught a quarter of the exact-draw grids at S 8 and about 70% of the BASS grids. Reading the margin instead, the rule first built in this slice, fired on every case measured. A share in coverage units would weigh the disagreement by the model's own predictive, which the layer exists not to lean on.

**What the readings showed** (the checks below).
- *The margin reads far wider than the error costs.* With the uniform bound, m̂ falls from 12.8 to 7.7 of 20 on exact draws as S grows from 40 to 2,560, and is 98 of 100 (Gaussian) and 94 (t) on the BASS fits at S 1,000. Carried bounds narrow it to 2.5 to 4.6 of 20, and to 72 and 52 of 100.
- *Why:* the band counts every unit whose rank adversarial errors of that size could move. With bounds of a few tenths of a nat beside log scores spread over a few nats, that is most units, while random errors move few decisions. The margin stays as a reported worst case and no longer drives the signal.
- *The gap tracks the error.* On exact draws its median falls from 0.157 at S 8 to 0.016 at S 2,560. It is 0.60 on the BASS Gaussian fits and 0.36 on the t fits, whose chains' r_eff widen the bounds.
- *The plain set lost no coverage* on any fixture: 0.905 at every S on exact draws, and within its error on the BASS fits. The inflated set covers 0.912 to 0.958 at 2 to 73% more width.

## Numerics

**The product.** For G candidates and the n_h training units of a stratum, the scores are one product: P (n_h × S), each unit's likelihoods scaled by their largest, P_is = exp(ℓ_is − maxₛ ℓ_is), times W̃ (S × G), each candidate's normalized weights. P is built once per stratum and call, W̃ per batch of 256 candidates. The kernel is `DenseProduct`: every entry a sum of fused multiply-adds in draw order.

**Why the order matters.** The candidate's own score is the same sum, four candidates to a vector, each lane in draw order, over its own likelihoods scaled the same way. So a candidate whose column equals a unit's gets that unit's score bit for bit, and a tie in the model is a tie in the step. A product that split the sum, by blocks or a tree, would be more accurate by an ulp and would break the ties the tie rule counts.

**Far ranges.** Scaling by each unit's largest likelihood keeps every entry of P in (0, 1], and the weights are normalized about their own largest, so neither overflows. An entry below 10⁻²⁰⁰ of its unit's largest likelihood is summed again in logs: above that floor, the terms that underflow in the product weigh less than S·2⁻¹⁰²²/10⁻²⁰⁰ of the entry. The candidate's own score takes the same rule, so equal columns still agree.

## Checks (slice 0)

Run on `ConformalBayesTests`, and the product on `DenseProductTests`.
- **CB-1 by enumeration** (`coverage_orderings`, `_eq`, `_le`, `coverage_strict_le`, `coverage_tieBreak`, `card_rankLe_le`). The model is discrete: three parameter values under a uniform prior, a unit (x, y) of six values, the posterior held exactly as three weighted draws. Over every data set of 2 to 5 units, 455 sets, each unit in turn the candidate:
  - the 56 sets with distinct scores cover exactly 1 − k/m at every k;
  - the 399 with ties cover at least 1 − k/m and at most that plus the straddling share, 720 (set, k) pairs with straddling ties among them;
  - the strict rule covers at most 1 − k/m, the tie-broken rule exactly 1 − k/m;
  - equal units tie exactly.
- **Per stratum** (`coverage_law`, `coverage_stratum`). The same model with x as the stratum, every multiset of 2 to 5 units in stratum 0 beside a fixed stratum 1: coverage at least 1 − k/m₀ within stratum 0, exactly that with distinct scores.
- **Coverage from truth outside the model**, exact posterior draws (S 400) of a linear Gaussian model with known unit variance, y | x, β ~ N(βx, 1). The truth is y = sin(1.5x) + (1 + 0.5|x|)(E − 1), E ~ Exp(1): curved, skewed and heteroscedastic. With n 19, α 0.1, k 2 over 3,850 replications:

| Set | Coverage (nominal 0.900, se 0.0048) | Mean width |
|---|---|---|
| The conformal set | 0.902 | 5.25, one interval each |
| The model's own 90% predictive interval | 0.798 | 3.39 |

- **Within strata.** x < 0 and x ≥ 0, 12 training units in each, the truth's noise scaled by 0.2 and by 3, the model the same. Ranked within its stratum, k 2 of m 13, nominal 0.846:
  - the quiet stratum covers 0.854 and the noisy one 0.824, se 0.0096, with mean widths 2.66 and 9.71;
  - the noisy stratum's reading is the fixture's seed. At 20,000 replications on three seeds it reads 0.835, 0.848 and 0.851; seed 611 at 40,000 reads 0.839, and its second 20,000 0.843. Pooled, 0.844 (z −1.2). The quiet stratum reads 0.851, 0.843 and 0.847;
  - raising S from 300 to 3,000 left the 1,400 decisions unchanged.
- **The small stratum** is the whole space below m = ⌈1/α⌉ at α 0.1, 0.15 and 0.5, and a set from there.
- **CB-4's refusals:** the collapsed content; an observation-level latent conditioned on or mixed; a candidate reading its latents otherwise; draws whose count, chains or weights differ. The counterexample, labelled integrated to pass the refusal, scores the data unit 5/8 and the candidate 1/2 and covers 0 against 1/2; the integrated reading ties them and covers.
- **The scores against a reference in logs**, a seeded sweep of 120 cases. It covers weights spanning e⁻⁷⁰⁰ with zeros among them, chains, log terms near −10⁴ that move by 900 nats across draws, candidates far in a tail and outside every draw's support, columns copied from units, and strata with one unit or none. Every score is within 10⁻¹² relative of the reference, with 62 exact ties and 481 unit scores 400 nats or more below the candidate's. The counts and π follow from the scores exactly, and the training matrix is left unchanged. A unit 800 nats below its best under the candidate's weights is summed in logs, to the reference's digits.
- **A bimodal predictive** gives a set of two intervals, its boundaries within a grid step of the closed form.

## Cost (slice 0)

The weights, the product and the ranks for one grid of G = 200 candidates, with the training side's once-per-stratum likelihoods beside them (`ConformalCostTests`, `BASS_QUALIFICATION=1`, the cycle clock, the best of five, Release). Indicative: the machine was shared.

| n, S | Training side | Weights | Product | Ranks | The full step | Per candidate (weights, product, ranks) |
|---|---|---|---|---|---|---|
| 1,000, 1,000 (the cost harness's base) | 15.8 ms | 2.5 ms | 14.3 ms | 1.2 ms | 36.7 ms | 12, 72, 6 µs |
| 1,000, 4,000 | 68.4 ms | 12.6 ms | 62.0 ms | 1.7 ms | 149.8 ms | 63, 310, 9 µs |

- *The product* runs at 31 and 27 GFLOP/s, 33 and 40 times the scalar loop of fused multiply-adds on the same operands.
- *The training side* is paid once per stratum and call. Many test cases in one call share it; a call per case pays it each time.
- *The exponentials* dominate the training side and the weights: n·S of them once, S·G per grid. A vector exponential would cut both, at the cost of bits that differ from `Math.Exp` across hardware.

## BASS: the candidate density

`BassFit.CandidateLogDensity(data, cases, grids, strata, unseenLevels)` (slice 1) is log p(y | x_t, θ_s) for every retained draw and every value y of each new case's grid, in the rows of `BassFit.PointwiseLogDensity(data)`, with its chain lengths or weights and its latent report.
- *The surface* at the cases is prediction's (`BassSurface`, the training basis at new inputs).
- *The density* is the training terms' own: N(y | f, σ²), or the noise law's with each observation's latent integrated out (`IBassNoiseLaw.LogDensity`), at the draw's σ² and its latents' globals. That covers Gaussian, Student t, quantile (asymmetric Laplace), Normal-Wald, skewed t and generalized hyperbolic BASS, and the layer reads them all by one density (CB-4's default).
- *The latent report* is the pointwise density's (`BassPointwise.Latents`, shared), so the layer's CB-4 check passes by construction and fails only if the two readings part.
- *A candidate is one observation* of a scalar response, so it pairs with the default units.
- *Refused:* a direct functional fit, whose unit is a curve (E17's whole-curve bands); a weighted fit, whose terms read a weight a new case does not carry; and draws whose working response moved, whose terms condition on an offset and weights a candidate does not have.

**The noise latent's report, corrected.** Until 2026-11-06 BASS reported a noise law's latent per observation as `IntegratedGivenOthers`, the collapsed content's treatment. That made `LeavesUnitsOut` true for every noise-law fit: WAIC refused them, PSIS-LOO's p_loo was NaN, and the layer would have refused them as leave-unit-out terms. Each latent is integrated against its prior, which no other observation moves given the globals, so the term is in sample. It is now reported `Integrated`, as [VarianceMeanMixtures](./VarianceMeanMixtures.md) always described it.

**Checks** (`BassCandidateTests`; coverage and cost with `BASS_QUALIFICATION=1`).
- *At the training cases it is the pointwise column, bit for bit,* for all six laws (120 draws × 80 cases each, largest gap 0): the same draws, chains and latent report. An SMC fit's candidates carry its population weights. The layer takes each pair.
- *Held-out coverage on data from outside every model.* The truth is f = sin(2πx₀) + x₁ with noise 0.2ζ plus, one time in seven, 0.8·Exp(1): contaminated and skewed. Each fit has n 99, so m 100 and k 10 at α 0.1, nominal 0.900. There are 20 fits per law with 50 held-out cases each, and the se is the fits' spread (binomial alone 0.0095). The width is read on 20 cases per fit, on a grid of step 0.08, and the model's own 90% interval from its predictive draws:

| Fit | Coverage (se) | Mean width | The model's own 90% interval | Its width | Least ESS at the true value |
|---|---|---|---|---|---|
| Gaussian | 0.905 (0.010) | 1.20 | 0.944 | 1.52 | 1.5 (8 of 1,000 below 50) |
| Student t (ν 3, c 0.7) | 0.894 (0.012) | 1.06 | 0.895 | 1.01 | 347 |
| Quantile (q 0.5) | 0.895 (0.012) | 1.02 | 0.921 | 1.22 | 87 |
| Normal-Wald (γ and β drawn) | 0.890 (0.015) | 1.07 | 0.900 | 1.18 | 74 |

  - Every law's set covers within its error, and it is narrower than the model's own interval where that overcovers: the Gaussian fit's σ² absorbs the contamination and widens its interval, while the density score still ranks the outliers low.
  - Two of 1,600 sets split into two intervals; none reached the grid's ends.
  - A Gaussian fit's add-one-in weights fall to an ESS of 1.5 at an outlying candidate, where the heavy-tailed laws keep 74 or more. That is where CB-3's margin (slice 2) will read.
- *Its cost,* for one grid of 200 values at the cost harness's base point (n 1,000, the rich surface, p 10, M_max 20), indicative on the cycle clock. Per candidate means the candidate's S density evaluations; the layer's step includes half the training side, which the two timed cases share.

| Law | S | The candidate density per grid | Per candidate | The layer's step per grid |
|---|---|---|---|---|
| Gaussian | 1,000 | 3.8 ms | 19 µs | 24.1 ms |
| Gaussian | 4,000 | 5.2 ms | 26 µs | 89.3 ms |
| Quantile | 1,000 | 5.9 ms | 29 µs | 24.0 ms |
| Student t | 1,000 | 12.7 ms | 64 µs | 23.4 ms |
| Generalized hyperbolic | 1,000 | 82.7 ms | 414 µs | 25.3 ms |
| Generalized hyperbolic | 4,000 | 289.6 ms | 1,448 µs | 82.6 ms |
| Normal-Wald | 1,000 | 179.5 ms | 898 µs | 21.2 ms |

  - Under the laws whose density is a Bessel function, the law's density is the slice's cost: 0.4 to 0.9 µs per evaluation, against 20 to 60 ns for the others.
  - Each call recomputes the law's normalizer, which is constant over a draw's grid. A density over many residuals at one draw would take it once; the pointwise density would share it.

## Checks and cost (slice 2)

Run on `ConformalMarginTests`, and with `BASS_QUALIFICATION=1` on `BassConformalMarginTests` and `ConformalCostTests`.
- **The band is the Lean count.** Under a uniform bound it is the literal `FragileEstimated` count over 3,000 random cases, with ties and bounds of 0 and ∞. Under per-unit bounds it is that count on continuous scores, where no comparison sits at rounding distance.
- **The deterministic statements on perturbed scores.** Exact scores with ties, per-unit bounds and estimates anywhere within them, over 4,000 cases:
  - `accepts_or_fragile`: 1,771 decisions were lost, every one fragile.
  - `acceptsInflated_of_accepts`, `accepts_shift_of_accepts_close` and `fragileEstimated_of_fragile` hold throughout.
  - `card_fragile_le_estimated` holds, and so does `fragileEstimated_mono`: the uniform band counts at least the per-unit one.
  - The adversary, the unit down by its bound and every other up by its own, refuses exactly the fragile units.
- **The layer against a reference.** Over 80 cases with weights, chains and far ranges:
  - its errors match Σψ²/r_eff read directly in logs, to 1.2·10⁻¹⁰;
  - its inflated ranks, uniform bounds and bands, uniform and carried, match the literal definitions on its own scores and bounds.
- **Inflated coverage on exact draws.** The conjugate fixture with S 40, so the scores carry real error, under noise the model does not share; n 19, k 2, 3,850 replications:
  - the inflated set covers 0.935 against 1 − k/m − ν = 0.89, and the plain set 0.905;
  - widths are 5.72 against 5.09;
  - the margin's reading at the true value averaged 0.485, which the plain set's coverage exceeds, as `coverage_close` says it must.
- **The margin as the draws grow** (exact draws, 100 data sets each, m 20):

| S | m̂, uniform | Coverage under error | m̂, carried | Under error, carried | Width, plain / inflated |
|---|---|---|---|---|---|
| 40 | 12.8 | 0.253 | 4.6 | 0.662 | 5.34 / 5.90 |
| 160 | 9.9 | 0.397 | 3.5 | 0.717 | 5.34 / 5.72 |
| 640 | 8.8 | 0.450 | 2.9 | 0.746 | 5.33 / 5.54 |
| 2,560 | 7.7 | 0.506 | 2.5 | 0.763 | 5.34 / 5.46 |

- **The gap as the draws grow**, the same data sets, with both rules' coverage over 2,000 true values (se 0.0067):

| S | Gap, median | 90th percentile | Above 0.25 | Above 0.5 | Fallback at 0.5 | Coverage, plain / inflated | Decision resting on the error |
|---|---|---|---|---|---|---|---|
| 8 | 0.157 | 0.313 | 24% | 1% | 100 | 0.905 / 0.941 | 3.6% |
| 16 | 0.132 | 0.253 | 11% | 0% | 100 | 0.905 / 0.941 | 3.6% |
| 40 | 0.103 | 0.189 | 5% | 0% | 58 | 0.906 / 0.937 | 3.1% |
| 160 | 0.063 | 0.115 | 1% | 0% | 8 | 0.906 / 0.927 | 2.1% |
| 640 | 0.032 | 0.072 | 0% | 0% | 3 | 0.906 / 0.918 | 1.3% |
| 2,560 | 0.016 | 0.039 | 0% | 0% | 1 | 0.905 / 0.912 | 0.7% |

  At S 8 and 16 the signal fires on trust: the add-one-in weights' tail is too short for k̂. At S 40, PSIS's threshold of 0.38 for that S fires it on 58 of 100. From S 160 it fires on a few.
- **The degenerate weights** (slice 1's Gaussian BASS fits, n 99, k 10, S 1,000, 20 fits × 50 held-out values, beside t fits of the same data):
  - *Where ESS is below 50* (8 values, ESS 1.5 to 48): the candidate's bound is 0.7 to 2.6 nats against a median 0.08 elsewhere, and π is 0.01 to 0.02. Both rules refuse every one, and the band holds 96 to 98 units under the uniform bound, 34 to 98 carried.
  - *So the degeneracy moves no decision there,* and the signal, read from the inflated set, stays quiet at all eight. Its refusal survives its bounds. The margin read wide there, as it did everywhere.
  - *The signal elsewhere:* at the true values it fires where the decision rests on the error, 53 of 992 Gaussian and 43 of 1,000 t values. Over the width grids the gap has median 0.60 (90th percentile 1.54) for the Gaussian fits and 0.36 (0.78) for the t fits, and fires on 218 and 120 of 400.
  - *Elsewhere* the medians are m̂ 98 (carried 72) for the Gaussian fits and 94 (carried 52) for the t fits.
  - *Coverage* is 0.905 plain and 0.958 inflated for the Gaussian, 0.894 and 0.937 for the t.
  - *Widths* are 1.20 against 2.07 for the Gaussian, 1.06 against 1.49 for the t.
  - *k̂* was read at 958 Gaussian and 937 t values in the inflated set, and passed PSIS's threshold at all but one t value.
- **The cost row**, per grid of 200 at n 1,000 (cycle clock, best of five, indicative). The weights include the errors' sums and each candidate's r_eff; the three products are the scores and the errors' two; the ranks include the bands and k̂:

| S, draws | Training side | Weights | Three products | Ranks | Full step | Slice 0's step |
|---|---|---|---|---|---|---|
| 1,000, independent | 16.2 ms | 3.3 ms | 41.4 ms | 62.2 ms | 127.9 ms | 36.7 ms |
| 1,000, two chains | 16.3 ms | 24.6 ms | 42.6 ms | 55.7 ms | 146.4 ms | |
| 4,000, independent | 63.6 ms | 16.9 ms | 171.2 ms | 122.8 ms | 383.4 ms | 149.8 ms |
| 4,000, two chains | 75.2 ms | 121.7 ms | 171.3 ms | 121.8 ms | 502.2 ms | |

  - The margin's always-on reading costs 2.6 to 3.5 times slice 0's step.
  - Its largest parts are the two error products, each candidate's r_eff along the chains (an FFT per candidate) and k̂ (a sort of S ratios per candidate in the inflated set).

## API at a glance

| Area | Types and members |
|---|---|
| The seams | `PointwiseLogDensity` (training); `CandidateLogDensity(values, grids, latents, stratumOf, chainLengths)` with `Weights`, `Offset(t)`, `Cases`, `Draws` |
| BASS's candidate density | `BassFit.CandidateLogDensity(data, cases, grids, strata, unseenLevels)`; `BassPointwise.Latents` (internal, shared) |
| The full step | `ConformalBayes.FullStep(training, candidates, ConformalOptions { Alpha, ErrorLevel, Inflate, FallbackShare, MarginBounds })` → `ConformalPrediction(Sets, Alpha)`; `ConformalBayes.Scores(training, candidates, column)` → `ConformalScores(Units, Scores, Candidate, EffectiveSampleSize)` |
| A case's set | `ConformalSet`: `Stratum`, `Units` (m), `Level` (k), `Coverage`, `WholeSpace`, `Grid`, `PValues`, `Below`, `Tied`, `CandidateScores`, `EffectiveSampleSizes`, `Accepts(g)`, `Intervals` (`GridInterval(Lower, Upper)`), `Resolution`, `OpenBelow`, `OpenAbove`, `Width` |
| CB-3's reading | `ConformalSet.Gap`, `AcceptsPlain(g)`, `AcceptsInflated(g)`; `ConformalSet.Margin` → `ConformalMargin` (`Multiple`, `CandidateErrors`, `UniformErrors`, `Fragile`, `InflatedRanks`, `ParetoK`, `ParetoKThreshold`, `Bounds`, `Count`, `Untrusted(k)`); `ConformalSet.Inflated`, `ErrorLevel`, `CoverageUnderError`, `FallBack`; `MarginBounds` (`Uniform`, `Carried`); `ConformalScores.Errors`, `CandidateError` |
| The product | `DenseProduct.Multiply`, `Pack`, `FusedDot`; `PackedColumns` |

## Not built yet

- A margin with bounds sure over every ordering, which needs each unit's error under every other unit's weights.
- The split form, with density, CQR and residual scores (slice 3).
- The conformal predictive distribution, with PIT and CRPS, and random tie-breaking with it (slice 4).
- The candidate seams of BARS, BAPS, the hierarchical P-spline, BGP and the changepoint posterior (slice 5), and BASS's for a direct functional fit's curves (E17).
- CB-4's latent-posterior reading, and Hoff's leave-one-out score.
- Refining a set's boundaries between grid values.

**See also:** [BayesianRegression](./BayesianRegression.md) (the pointwise seam) · [VarianceMeanMixtures](./VarianceMeanMixtures.md) (the noise laws whose latents are integrated) · [LinearAlgebra](./LinearAlgebra.md) (the product)
