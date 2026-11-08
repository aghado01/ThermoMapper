# Conformal prediction

> The conformal layer (E16): prediction sets whose coverage holds in finite samples under exchangeability, whatever the model, from an engine's retained draws and pointwise densities with no refit. The layer owns the guarantee and the score owns the width, so the layer is built once, engine-free, and each engine contributes densities. Built: the full step (Fong and Holmes's conformal Bayes) over a training density and a candidate density; BASS's candidate density for every noise law; CB-3's reading of the scores' Monte Carlo error, with the margin, the inflated rule and two diagnostics; and the split form, with density, CQR and residual scores. The predictive distribution and the other engines' candidate seams follow.

## Overview

Conformal prediction ranks a candidate observation among the data by a conformity score that treats every unit alike. Under exchangeability the candidate's rank is uniform, so a set that refuses only the lowest ranks covers a new observation at a level fixed in advance, whether or not the model is right. The model only decides which units score high, and so how wide the set is.

What it serves, by problem class:
- prediction sets for a new observation that hold whether or not the noise model is right;
- coverage within strata, for partially exchangeable data (Mondrian conformal prediction, Vovk; grouped conformal Bayes, Fong and Holmes);
- a guarantee that does not lean on the sampler: the split form, valid for any fit, a chain stopped early or a wrong model included;
- intervals whose width follows the noise, from any fit's quantiles (CQR) or centre and scale;
- later, a calibrated predictive distribution at each x.

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

**Two diagnostics on every set, and no trigger** (decided by the user on 2026-11-07). The layer reports how far the scores' Monte Carlo error could move each set, and switches to nothing on its own.
- *The inflated-set gap* (`Gap`). The inflated set holds the plain one, since raising the candidate's estimate and lowering the units' only adds units at or below it. The values only the inflated rule accepts are those whose decision rests on the error, as a share of the plain set: (|I| − |A|)/|A|, read on the grid's cells as `Width` reads them, or by count on a grid of one value.
- *The split-half mass* (`SplitHalfMass`, with `SplitHalfDisagreements`), the main session's construction and not from a paper. The retained draws split into two halves: each chain's first and second halves, or alternate draws without chains. The layer scores each half on the same grid, by two half-size score products, and reads the mass, under the full draws' predictive at x, of the values where the two halves' plain sets disagree. It is the coverage the error could shift here, read at its typical size with no union over units, and conservative, since each half carries √2 the full run's error. It is NaN on a grid of one value.
- *The fallback rule is an opt-in* (`FallbackShare`, null by default). With a share stated, `FallBack` fires when the gap exceeds it, or when a candidate in the inflated set has an untrusted k̂. Nothing calls it by default.
- *Why no trigger.* On the BASS fits the plain set covers at nominal, yet the split-half mass at 0.1α fires on 61 to 74% of cases and the gap at 0.5 on 32 to 54%. A default trigger would send most BASS cases to the split form, giving up data for no coverage gain. A trigger can be adopted once a fixture with chains like BASS's shows a measured shortfall to set its threshold by. None does yet, so none is scheduled.
- *Open, no decision waiting:* the main session's reading of why the halves disagree on BASS while coverage holds. BASS's chains move between structures slowly, so the two halves can sit in different regions of the posterior, and a chain held in one region may still score every unit nearly alike, where thin draws do not. If so, the split-half mass reads the set's instability, which differs from miscoverage. It is a reading, not a result.

**What the readings showed** (the checks below).
- *The margin reads far wider than the error costs.* With the uniform bound, m̂ falls from 12.8 to 7.7 of 20 on exact draws as S grows from 40 to 2,560, and is 98 of 100 (Gaussian) and 94 (t) on the BASS fits at S 1,000. Carried bounds narrow it to 2.5 to 4.6 of 20, and to 72 and 52 of 100.
- *Why:* the band counts every unit whose rank adversarial errors of that size could move. With bounds of a few tenths of a nat beside log scores spread over a few nats, that is most units, while random errors move few decisions. The margin stays as a reported worst case.
- *The gap tracks the error.* On exact draws its median falls from 0.157 at S 8 to 0.016 at S 2,560. It is 0.60 on the BASS Gaussian fits and 0.36 on the t fits, whose chains' r_eff widen the bounds.
- *The plain set lost no coverage* on any fixture: 0.905 at every S on exact draws, and within its error on the BASS fits. The inflated set covers 0.912 to 0.958 at 2 to 73% more width.

**The two diagnostics measured side by side** before they were built in, the split-half then read from two public steps (2026-11-07, at the user's request).
- *How it reads.* It splits the retained draws into two halves: each chain's first and second halves where there are chains, alternate draws otherwise. It builds the plain set from each half on the same grid, then reads the mass, under the full draws' predictive at x, of the region where the two half-sets disagree.
- *What the mass means.* It estimates how much coverage the Monte Carlo error could shift for this case. It reads the error at its typical size, with no union over units and no worst case, and is conservative, since each half carries √2 the full run's error.
- *Read as signals,* the gap through the rule at 0.5 and the mass at 0.1α. In the table, "the layer's signal" is that rule.

Both rules were read case by case: the signal on a grid at x, and the plain set's decision at the true value at the same x. So a rule's firing can be set against the coverage the plain set achieved. A useful rule fires where the plain set falls short and stays quiet where it holds. Each entry below is the fire rate, then coverage where the rule fires against where it is quiet:

| Fixture (1,000 cases each) | Plain coverage | The layer's signal | Split-half mass above 0.1α |
|---|---|---|---|
| Exact draws, S 40 | 0.907 | 59%: 0.912 / 0.899 | 20%: 0.827 / 0.927 |
| Exact draws, S 160 | 0.904 | 8%: 0.924 / 0.902 | 10%: 0.800 / 0.916 |
| Exact draws, S 640 | 0.888 | 2%: 0.833 / 0.889 | 6%: 0.738 / 0.898 |
| Exact draws, S 2,560 | 0.905 | 1%: 1.000 / 0.904 | 3%: 0.923 / 0.905 |
| Two chains, ρ 0.9 | 0.888 | 5%: 0.902 / 0.887 | 16%: 0.811 / 0.903 |
| Two chains, ρ 0.99 | 0.886 | 15%: 0.759 / 0.908 | 34%: 0.794 / 0.933 |
| Two chains, ρ 0.995 | 0.889 | 13%: 0.858 / 0.894 | 34%: 0.854 / 0.907 |
| Two chains, ρ 0.999 | 0.887 | 5%: 0.635 / 0.901 | 24%: 0.753 / 0.929 |
| BASS Gaussian, S 1,000 | 0.905 | 54%: 0.907 / 0.903 | 74%: 0.911 / 0.889 |
| BASS t, S 1,000 | 0.894 | 32%: 0.896 / 0.893 | 61%: 0.893 / 0.895 |

The chains are AR(1) about the exact posterior, two chains of 500, so they mix more slowly as ρ grows.

- *On exact draws and slow chains* the split-half mass discriminates. Where it fires, the plain set covers 0.75 to 0.85; where it is quiet, 0.90 to 0.93.
- *The layer's signal discriminates only for the slowest chains,* and fires rarely there. At small S its trust trigger fires where the plain set holds: at S 40, 59% of cases, covering 0.912. The gap alone above 0.5 behaves as the signal does on the chains and almost never fires on exact draws.
- *On the BASS fits neither discriminates,* because the plain set loses nothing there. The split-half mass is large anyway (median 0.033 Gaussian, 0.014 t) and fires on 74% and 61% of cases at 0.1α. 0.1α is a starting point, not a calibrated share.
- *At slice 1's degenerate Gaussian weights* (ESS 1.5 to 48), all eight true values are misses, outlying truths the plain set refuses. The layer's signal fires at six of their grids, the split-half at seven, with masses 0 to 0.125.

## The split form (CB-2)

**What it is.** The data split in two. A fit is made on the training split, and held-out calibration units are scored by it. For a test case in stratum h with n_h calibration units, m = n_h + 1 and k = ⌊α·m⌋. A candidate y is accepted when its score reaches the k-th smallest calibration score: the full step's rule, ties counted for the candidate (`SplitCalibration`). The set covers at least 1 − k/m, and exactly 1 − k/m with distinct scores (`split_coverage`).

**What it buys: a guarantee that does not lean on the sampler.** The fit may be any function of the training split, its sampler's seeds and its whole chain included (`split_coverage` holds for every `fit`; `splitScore_equivariant`). A chain stopped early, a wrong model and draws that never mixed all cover; they only widen the set. The draws' Monte Carlo error is part of the fit, so it costs no coverage. A split set carries no margin, no ESS and no diagnostics, and its `CoverageUnderError` is its `Coverage`. The hypothesis is that the calibration and test units are exchangeable given the training split, which exchangeability of every unit gives (`exchangeable_given_training`). The layer cannot check that the fit never saw the calibration units; the caller's split must hold it.

**What it costs: the calibration units.** The fit does not see them, so its scores are weaker and its sets wider, and the level steps in 1/m of the calibration's m, not the data's. So the full step suits good draws and scarce data, and the split form suits a guarantee that must not lean on the sampler (E16's reading). Neither replaces the other, and nothing switches between them on its own (the user's decision of 2026-11-07).

**Three scores,** each higher conforming better:
- *The predictive density* (`SplitConformal.Density`). The score is ln p̂(y | x) = ln Σ_s W_s p(y | x, θ_s), read through the candidate seam: the calibration cases each on a grid of one value, its response, and the test cases on their grids. The arithmetic is the full step's predictive, bit for bit, so equal columns tie. The set is π on the grid, a union of intervals, reported as a `ConformalSet` under `ConformalForm.Split`. It is the predictive's highest-density region, its level set by the calibration rather than by the model.
- *CQR's* (`SplitConformal.Quantiles`; Romano, Patterson and Candès 2019). The score is min(y − q_lo(x), q_hi(x) − y), Lean's `cqrScore` and the negative of their max(q_lo − y, y − q_hi). Its set is the interval [q_lo − Q, q_hi + Q], where Q is the negated threshold, the k-th largest calibration max(q_lo − y, y − q_hi) (`cqr_ordConnected`). It is empty when the test case's band q_hi − q_lo is less than −2Q. That needs Q < 0, which happens when most calibration units fall well inside their bands, or quantiles that cross by more than 2Q. Any fit's quantiles serve: a quantile fit's surface (quantile BASS) or a predictive's quantiles (`BassPrediction.Quantile`).
- *The standardized residual* (`SplitConformal.Residuals`). The score is −|y − μ(x)|/ρ(x), and its set is μ ± Qρ, where Q is the k-th largest calibration |y − μ|/ρ. The scale is Lei et al.'s ρ, since σᵢ on this page are the full step's scores. With no scale it is the absolute residual. Lei et al.'s locally weighted form is this score with ρ a fitted spread. BASS's `ObservationMean` with the predictive's interquartile range, `Quantile(0.75) − Quantile(0.25)`, serves every noise law, since no law reports its variance.

**Closed forms where the set is an interval.** CQR's and the residual's sets are intervals with ends in closed form, so they need no grid. Their ends are exact to rounding, and each case costs one addition after the calibration's sort. The density's set is not an interval in general, so it is read on a grid.

**Strata** (Mondrian; `coverage_stratum`). Each test case is ranked among its stratum's calibration units. Below m = ⌈1/α⌉ the set is the whole line, and a stratum with no calibration unit has m = 1.

**One function of the unit.** For the density, the calibration's and the test cases' candidate densities must hold the same draws, chain lengths, weights and latent report, or the layer refuses them. It also refuses a calibration unit read on a grid of more than one value.

**Not built,** with what each would change:
- CQR's asymmetric form, a threshold per tail at α/2 each, which Romano, Patterson and Candès also give. It would cover each tail separately, at more width.
- Cross-conformal and the jackknife+ (Barber, Candès, Ramdas and Tibshirani, 2021). They use every unit in both roles at 1 − 2α, at the cost of K fits or n.
- The tournament correction, recorded in E16, which keeps all the data at 1 − 2α for any sampler, at n leave-one-out chains.

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

**The two diagnostics' cost** (2026-11-07), per grid of 200 at n 1,000, on the cycle clock, best of five. Indicative: the machine was loaded, and parts unchanged since slice 2 read about 40% slower this time.
- *The split-half mass* costs 21 to 25 ms per grid at S 1,000 and 82 to 97 ms at S 4,000, about one more score product. Its halves' likelihoods are built once per stratum, 15 to 65 ms.
- *The gap* reads the inflated ranks the margin's errors give, so it costs almost nothing beside them.
- *The full step* with both diagnostics and the margin:

| S, draws | Full step | Slice 0's plain step |
|---|---|---|
| 1,000, independent | 185.8 ms | 36.7 ms |
| 1,000, two chains | 203.4 ms | |
| 4,000, independent | 643.9 ms | 149.8 ms |
| 4,000, two chains | 700.1 ms | |

- *Most of it is still the margin:* the two error products, the bands, k̂ and each candidate's r_eff.
- *The gap needs only the per-position errors,* the two error products and r_eff, not the bands or k̂. So the margin's always-on question, slice 2's, would leave the gap's cost and remove the bands' and k̂'s, 60 to 160 ms per grid in these rows.

## Checks and cost (slice 3)

**Checks** (`SplitConformalTests`; `BassSplitConformalTests` and the cost row with `BASS_QUALIFICATION=1`).
- *The rule against the count* (`card_rankLe_le`): 2,000 seeded calibrations with tied and infinite scores, strata of one and strata with none. Every probe agrees.
- *CB-2 by exact enumeration* (`split_coverage`, `exchangeable_given_training`). There are 200 laws, each with three training outcomes, units on five values, m from 2 to 5, and a random fit, wrong ones among them. Coverage is at least 1 − k/m in every law, and exactly 1 − k/m for injective fits on distinct values (`coverage_orderings_eq`).
- *The density score* equals a reference in logs to 10⁻¹², and the full step's predictive bit for bit. This held over 300 sweeps of weighted draws, chains, terms near −800 and +600, −∞ columns and strata. A test column copied from a calibration unit ties with it.
- *CQR's set is an interval* (`cqr_ordConnected`), and so is the residual's. The sweep covered 1,500 trials (7,980 sets, 550 of them empty) with crossing quantiles, tied responses and scales from 10⁻¹⁵⁰ to 10¹⁵⁰. On a grid through and past each set, the rank rule read by the count runs unbroken, and it agrees with the closed form except within rounding of an end.
- *Coverage for any fit on the conjugate model,* under the slice 0 truth the model does not share. There are 30 training units and 19 calibration units (m 20, k 2, nominal 0.900), and 2,000 replications of 20 test cases. The se, 0.0021, comes from coverage's Beta law given the calibration. Each cell is coverage (mean width), the density's width read on 600 sets per fit, none reaching the grid's ends:

| Fit | Its own 90% interval | Density | CQR | Residual |
|---|---|---|---|---|
| The exact posterior | 0.791 | 0.900 (5.36) | 0.900 (5.31) | 0.900 (5.30) |
| A wrong fit | 0.237 | 0.899 (6.43) | 0.898 (6.51) | 0.896 (12.39) |
| A chain stopped early | 0.243 | 0.899 (17.62) | 0.900 (17.64) | 0.899 (17.54) |

  - *The wrong fit* flips the slope's sign under an overconfident noise sd of 0.25, crosses CQR's quantiles, and shapes the residual's scale against the noise. *The chain stopped early* is forty Metropolis steps from β = 5.
  - Both cover, while their own intervals cover about a quarter of cases. Wrongness costs width, not coverage, and the misshaped scale costs the most.
- *Coverage for any BASS fit* on slice 1's contaminated data. There are 99 training units and 199 calibration units (m 200, k 20), and 30 fits of 100 test cases each, se 0.0067. Each cell is coverage (mean width):

| Fit | Its own 90% interval | Density | CQR | Residual |
|---|---|---|---|---|
| A chain stopped early (one chain, no burn-in, 20 retained draws) | 0.951 (2.50) | 0.908 (1.87) | 0.912 (2.02) | 0.911 (1.92) |
| A wrong fit (it never sees x₀) | 0.946 (2.81) | 0.909 (2.50) | 0.909 (2.50) | 0.910 (2.50) |
| A healthy fit | 0.925 (1.60) | 0.907 (1.21) | 0.909 (1.34) | 0.909 (1.26) |
| Quantile BASS at 0.05 and 0.95, CQR on their surfaces | 0.853 (1.27) | | 0.906 (1.48) | |

  - Every reading lies within two se of 0.900, and all sit above it together: the four facts share their 30 data sets, so one calibration draw moves them alike.
  - CQR corrects the quantile pair's undercoverage at 17% more width. One of its 3,000 intervals is empty.
  - No density set split into two intervals or reached the grid's ends.

**Cost** at n 1,000 calibration units and G 200, on the cycle clock, best of five. Indicative, on a loaded machine:

| S | Calibration, once per stratum | Per grid | Per candidate | The step on one grid | CQR, calibration and 200 intervals | The residual's |
|---|---|---|---|---|---|---|
| 1,000 | 12.0 ms | 2.1 ms | 10.6 µs | 14.5 ms | 0.70 ms | 0.65 ms |
| 4,000 | 48.9 ms | 7.3 ms | 36.5 µs | 54.4 ms | 0.33 ms | 0.40 ms |

- *The density score costs one predictive per candidate.* It needs no product over the units and no margin, against the full step's 644 to 700 ms per grid at S 4,000.
- *On BASS the readings cost more than the fits.* Each Gaussian fact at S 1,000 reads four predictive quantiles (`BassPrediction.Quantile`) at 299 cases, and took 812 to 860 s for 30 replications. The chain stopped early, at S 20, took 17 s. The quantile pair took 27 s for 60 fits, since its readings are surface means.

## API at a glance

| Area | Types and members |
|---|---|
| The seams | `PointwiseLogDensity` (training); `CandidateLogDensity(values, grids, latents, stratumOf, chainLengths)` with `Weights`, `Offset(t)`, `Cases`, `Draws` |
| BASS's candidate density | `BassFit.CandidateLogDensity(data, cases, grids, strata, unseenLevels)`; `BassPointwise.Latents` (internal, shared) |
| The full step | `ConformalBayes.FullStep(training, candidates, ConformalOptions { Alpha, ErrorLevel, Inflate, FallbackShare (opt-in, null by default), MarginBounds })` → `ConformalPrediction(Sets, Alpha)`; `ConformalBayes.Scores(training, candidates, column)` → `ConformalScores(Units, Scores, Candidate, EffectiveSampleSize)` |
| A case's set | `ConformalSet`: `Stratum`, `Units` (m), `Level` (k), `Coverage`, `WholeSpace`, `Grid`, `PValues`, `Below`, `Tied`, `CandidateScores`, `EffectiveSampleSizes`, `Accepts(g)`, `Intervals` (`GridInterval(Lower, Upper)`), `Resolution`, `OpenBelow`, `OpenAbove`, `Width` |
| CB-3's reading | `ConformalSet.Gap`, `SplitHalfMass`, `SplitHalfDisagreements`, `AcceptsPlain(g)`, `AcceptsInflated(g)`; `ConformalSet.Margin` → `ConformalMargin` (`Multiple`, `CandidateErrors`, `UniformErrors`, `Fragile`, `InflatedRanks`, `ParetoK`, `ParetoKThreshold`, `Bounds`, `Count`, `Untrusted(k)`); `ConformalSet.Inflated`, `ErrorLevel`, `CoverageUnderError`, `FallBack`; `MarginBounds` (`Uniform`, `Carried`); `ConformalScores.Errors`, `CandidateError` |
| The split form | `SplitConformal.Density(calibration, candidates, alpha)` → `ConformalPrediction(Sets, Alpha, Form)`; `SplitConformal.Quantiles(QuantileReadings(Lower, Upper, StratumOf), response, test, alpha)` and `SplitConformal.Residuals(ResidualReadings(Center, Scale, StratumOf), response, test, alpha)` → `ConformalIntervals(Intervals, Alpha, Calibration)`; `ConformalInterval`: `Units`, `Level`, `Threshold`, `Lower`, `Upper`, `WholeSpace`, `Coverage`, `Empty`, `Width`, `Contains(y)`; `SplitConformal.Calibrate(scores, stratumOf, alpha)` → `SplitCalibration`: `Units(h)`, `Level(h)`, `Threshold(h)`, `Rank(h, s)`, `PValue(h, s)`, `Accepts(h, s)`; `ConformalForm` (`Full`, `Split`) |
| The product | `DenseProduct.Multiply`, `Pack`, `FusedDot`; `PackedColumns` |

## Not built yet

- A margin with bounds sure over every ordering, which needs each unit's error under every other unit's weights.
- CQR's asymmetric form, cross-conformal and the jackknife+ (above).
- The conformal predictive distribution, with PIT and CRPS, and random tie-breaking with it (slice 4).
- The candidate seams of BARS, BAPS, the hierarchical P-spline, BGP and the changepoint posterior (slice 5), and BASS's for a direct functional fit's curves (E17).
- CB-4's latent-posterior reading, and Hoff's leave-one-out score.
- Refining a set's boundaries between grid values.

**See also:** [BayesianRegression](./BayesianRegression.md) (the pointwise seam) · [VarianceMeanMixtures](./VarianceMeanMixtures.md) (the noise laws whose latents are integrated) · [LinearAlgebra](./LinearAlgebra.md) (the product)
