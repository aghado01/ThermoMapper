# Conformal prediction

> The conformal layer (E16): prediction sets whose coverage holds in finite samples under exchangeability, whatever the model, from an engine's retained draws and pointwise densities with no refit. The layer owns the guarantee and the score owns the width, so the layer is built once, engine-free, and each engine contributes densities. Built: the full step (Fong and Holmes's conformal Bayes) over a training density and a candidate density; BASS's candidate density for every noise law; CB-3's reading of the scores' Monte Carlo error, with the margin, the inflated rule and two diagnostics; and the split form, with density, CQR and residual scores. And the conformal predictive distribution, with PIT, quantiles and CRPS (slice 4): the split form's, on the fit's predictive CDF, by default, and the full step's as an opt-in that flags the cases where its score falls, at the user's decision of 2026-11-08. The other engines' candidate seams follow.

## Overview

Conformal prediction ranks a candidate observation among the data by a conformity score that treats every unit alike. Under exchangeability the candidate's rank is uniform, so a set that refuses only the lowest ranks covers a new observation at a level fixed in advance, whether or not the model is right. The model only decides which units score high, and so how wide the set is.

What it serves, by problem class:
- prediction sets for a new observation that hold whether or not the noise model is right;
- coverage within strata, for partially exchangeable data (Mondrian conformal prediction, Vovk; grouped conformal Bayes, Fong and Holmes);
- a guarantee that does not lean on the sampler: the split form, valid for any fit, a chain stopped early or a wrong model included;
- intervals whose width follows the noise, from any fit's quantiles (CQR) or centre and scale;
- a calibrated predictive distribution at each x, with its PIT, quantiles and CRPS.

It does not reach a curve's or a parameter's own intervals: those stay with [calibrated intervals and resampling](../../issues/regression/bass/extensions-brief.md#calibrated-intervals-for-loss-based-posteriors).

Where things are:
- *The contract*, CB-1 to CB-4, is proved in `lean/Enthymemata/ConformalValidity.lean`, with its statements, checks and decisions in [conformal-lemmas.md](../../lean/Protolemmata/Conformal/conformal-lemmas.md). Every contract test below names its declaration.
- *The build* is [the layer's brief](../../issues/regression/conformal-layer-brief.md), slice by slice, and E16 in the [extensions brief](../../issues/regression/bass/extensions-brief.md#conformal-bayes-the-consumer) holds the literature.
- *The seams it reads* are [BayesianRegression](./BayesianRegression.md)'s pointwise density (OQ17), the draws × units matrix of log p(y_unit | draw), and its sibling at new cases, the candidate density.
- *The product it runs on* is [LinearAlgebra](./LinearAlgebra.md#dense-product-denseproduct)'s `DenseProduct`.

## The seams

**The training density** (`PointwiseLogDensity`): log p(y_unit | θ_s) for each retained draw s and each training unit, with the unit partition and its stratum map, how each latent was treated (`LatentReport`), the chain lengths, and draw weights for a weighted population.

**The candidate density** (`CandidateLogDensity`): log p(y | x_t, θ_s) for the same draws in the same rows, at each test case t and each value y of t's grid. The columns run case by case, each case's grid ascending; each case carries a stratum. Its latents, chain lengths and weights must be the training density's. A value outside a draw's support takes −∞. An engine joins the layer by producing it; the layer knows no engine.

**Each draw's log CDF** (`LogCdf` on both, slice 4): ln F(y | x, θ_s) in the density's layout, the latents treated alike, for a score that orders one value, so the training units must be one observation each. The log keeps both tails: with ℓ = ln F, F = e^ℓ, and 1 − F = −expm1(ℓ) to full relative precision near 1, where a probability rounds to 1. Null where the engine gave none.

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
- *Random tie-breaking* (`coverage_tieBreak`) gives exactly 1 − k/m. The predictive distribution (slice 4) takes it, as τ in Q(y, τ), where a uniform p-value is wanted. The layer reports, for each grid value, how many units score strictly below (`Below`) and how many tie (`Tied`), so any rule reads from the step.
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

**The margin** (`coverage_close` with `card_fragile_le_estimated`) is an opt-in (`MarginBounds`, null by default; decided below). For each grid value taken as the truth, it counts the band of units whose decision estimates within their bounds could change (`Fragile`). m̂ is the largest count over the grid values inside some draw's support, since any of them may be the truth. The plain set's coverage under error is then 1 − (k + m̂)/m − ν (`CoverageUnderError`).
- *With carried bounds* (`MarginBounds.Carried`, at the user's order of 2026-11-08). `coverage_close` needs bounds carried by the units, each holding in every ordering. So each unit takes its largest bound over the orderings the draws represent: the held one, and each in which a unit of the stratum, itself included, is the candidate. The candidate takes the larger of its own bound and its bounds under each unit's weights.
- *Or with the uniform bound* (`MarginBounds.Uniform`): the largest carried bound, for every unit, the Lean's uniform ε. Its band holds the carried one (`fragileEstimated_mono`), so the carried margin is never the wider, and it falls back to the uniform one where no unit stands apart.
- *Changed on 2026-11-08* ([below](#the-per-unit-bounds-between-slices-4-and-5)). Slice 2 read the uniform bound from the held ordering alone and the carried bound from two orderings, the held one and the unit's own turn. Neither dominated the units' errors over the orderings, which the band's guarantee needs.
- *Not read:* each ordering's own draws. The other orderings' errors come from the held draws reweighted by each unit's density, not from a run on each ordering's data.

**Decided by the user on 2026-11-08: the margin is opt-in** ("Make it opt-in").
- *What the layer does.* By default it reads no margin and pays none of its cost: no band, no carried bounds and no per-stratum orderings. `MarginBounds.Carried` or `MarginBounds.Uniform` asks for it. Not asked for, `ConformalMargin.Computed` is false, the band, its count and the carried bounds are null, and the plain set's `CoverageUnderError` is NaN: not read, rather than read as zero or as clear.
- *What stays on every set:* the inflated rule's readings, per position by default with its uniform option, and both diagnostics below. Neither the fallback rule nor the split-half mass reads the margin's band; the fallback reads the gap and k̂, which are the inflated rule's.
- *Why.* On BASS the margin reads 96 to 100 units of 100 under either bound, so it fires on every fit and says little. It costs the orderings, 0.33 to 1.8 s per stratum at n 1,000. The inflated rule already gives the protection, with no per-stratum cost. It still pays where a unit stands apart (4.9 against 16.3 of 21 in the tail fixture), so it stays available.
- *What would change it:* a fit where the margin reads a useful bound, or a use that needs the plain set's coverage under error.

**The inflated rule** (`coverage_inflated`, Ndiaye's stable set) is an option (`Inflate`, off by default) and always reported (`InflatedRanks`). It accepts y when more than k units have an estimate lowered by its bound at most the candidate's raised by its own. The rule keeps 1 − k/m − ν at the cost of width.
- *Per-position bounds by default* (`InflationBounds.PerPosition`): each unit's under this candidate's weights, and the candidate's own. `coverage_inflated` takes bounds that differ by ordering and position, so the held ordering serves alone, and a candidate's wide error widens only its own comparisons.
- *Or one uniform bound* (`InflationBounds.Uniform`), the largest per-position one, kept for comparison. Its set holds the per-position set.

**How each bound is estimated, and why it dominates** (2026-11-08). A bound must be at least the error of the estimate it covers, in every ordering the guarantee averages over (`fragileEstimated_mono` lets any larger bound serve).
- *The estimate.* Each bound is z of an estimate's Monte Carlo standard errors, z = Φ⁻¹(1 − ν/(2m)), the standard error being √log(1 + Σψ²/r_eff) under the weights that ordering puts on the draws, with r_eff the weighting density's along the chains.
- *Per position.* In the held ordering, the m estimates each lie within their own bound with chance at least 1 − ν, by the union over the m positions under the normal approximation. That is what `coverage_inflated` asks of an ordering, and no position's bound reads another's error.
- *Uniform, for the inflated rule.* The largest per-position bound is at least each of them, so it dominates them by construction.
- *Carried.* A unit's bound must hold in every ordering. The layer reads each unit's error under the weights of every unit of the stratum as the candidate, the held draws reweighted by that unit's density: three (n × S)·(S × n) products once per stratum. The candidate's error under each unit's weights takes one more (n × S)·(S × G) product per batch, from sums the step already holds. A unit's carried bound at y is the largest of its held bound and these, so it dominates its error in every ordering the draws represent, by the same union in each.
- *Uniform, for the margin.* The largest carried bound is at least every carried bound, so it dominates them by construction.
- *Still indicative,* as slice 2's bounds were. The standard error is itself estimated, and the normal approximation fails where the weights degenerate, which k̂ reads. The other orderings are an importance sample from the held draws, not draws of each ordering's own posterior.

**Trust** is PSIS's Pareto k̂ of the add-one-in log ratios (`PredictiveChecks`' k̂, shared). It is read where the inflated rule accepts, since there a decision rests on the bounds; above PSIS's threshold at S, the error estimate there cannot be trusted.

**Two diagnostics on every set, and no trigger** (decided by the user on 2026-11-07). The layer reports how far the scores' Monte Carlo error could move each set, and switches to nothing on its own.
- *The inflated-set gap* (`Gap`). The inflated set holds the plain one, since raising the candidate's estimate and lowering the units' only adds units at or below it. The values only the inflated rule accepts are those whose decision rests on the error, as a share of the plain set: (|I| − |A|)/|A|, read on the grid's cells as `Width` reads them, or by count on a grid of one value.
- *The split-half mass* (`SplitHalfMass`, with `SplitHalfDisagreements`), the main session's construction and not from a paper. The retained draws split into two halves: each chain's first and second halves, or alternate draws without chains. The layer scores each half on the same grid, by two half-size score products, and reads the mass, under the full draws' predictive at x, of the values where the two halves' plain sets disagree. It is the coverage the error could shift here, read at its typical size with no union over units, and conservative, since each half carries √2 the full run's error. It is NaN on a grid of one value.
- *The fallback rule is an opt-in* (`FallbackShare`, null by default). With a share stated, `FallBack` fires when the gap exceeds it, or when a candidate in the inflated set has an untrusted k̂. Nothing calls it by default.
- *Why no trigger.* On the BASS fits the plain set covers at nominal, yet the split-half mass at 0.1α fires on 61 to 74% of cases and the gap at 0.5 on 32 to 54%. A default trigger would send most BASS cases to the split form, giving up data for no coverage gain. A trigger can be adopted once a fixture with chains like BASS's shows a measured shortfall to set its threshold by. None does yet, so none is scheduled.
- *Open, no decision waiting:* the main session's reading of why the halves disagree on BASS while coverage holds. BASS's chains move between structures slowly, so the two halves can sit in different regions of the posterior, and a chain held in one region may still score every unit nearly alike, where thin draws do not. If so, the split-half mass reads the set's instability, which differs from miscoverage. It is a reading, not a result.

**What the readings showed** (the checks below).
- *The margin reads far wider than the error costs.* With the uniform bound, m̂ falls from 14.1 to 8.4 of 20 on exact draws as S grows from 40 to 2,560, and its median is 100 of 100 (Gaussian) and 98 (t) on the BASS fits at S 1,000. Carried bounds narrow it to 5.0 to 2.6 of 20, but only to 97 and 96 of 100 on BASS, where a unit far out, as the candidate, widens every unit's error (read on 2026-11-08 with the bounds over every ordering; slice 2's, from fewer orderings, read 12.8 to 7.7 and 2.5 to 4.6 of 20, and 98, 94, 72 and 52 of 100).
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
- *The standardized residual* (`SplitConformal.Residuals`). The score is −|y − μ(x)|/ρ(x), and its set is μ ± Qρ, where Q is the k-th largest calibration |y − μ|/ρ. The scale is Lei et al.'s ρ, since σᵢ on this page are the full step's scores. With no scale it is the absolute residual. Lei et al.'s locally weighted form is this score with ρ a fitted spread. BASS's `ObservationMean` with the predictive's interquartile range, `Quantile(0.75) − Quantile(0.25)`, serves every noise law: a law's variance (`IBassNoiseLaw.Variance`) is infinite for t at ν ≤ 2, and prediction reports no predictive variance per case.

**Closed forms where the set is an interval.** CQR's and the residual's sets are intervals with ends in closed form, so they need no grid. Their ends are exact to rounding, and each case costs one addition after the calibration's sort. The density's set is not an interval in general, so it is read on a grid.

**Strata** (Mondrian; `coverage_stratum`). Each test case is ranked among its stratum's calibration units. Below m = ⌈1/α⌉ the set is the whole line, and a stratum with no calibration unit has m = 1.

**One function of the unit.** For the density, the calibration's and the test cases' candidate densities must hold the same draws, chain lengths, weights and latent report, or the layer refuses them. It also refuses a calibration unit read on a grid of more than one value.

**Not built,** with what each would change:
- CQR's asymmetric form, a threshold per tail at α/2 each, which Romano, Patterson and Candès also give. It would cover each tail separately, at more width.
- Cross-conformal and the jackknife+ (Barber, Candès, Ramdas and Tibshirani, 2021). They use every unit in both roles at 1 − 2α, at the cost of K fits or n.
- The tournament correction, recorded in E16, which keeps all the data at 1 − 2α for any sampler, at n leave-one-out chains.

## The predictive distribution (slice 4): the reading

**What was read.** Vovk et al.'s conformal predictive distribution reads π(y) as a CDF in y (Vovk, Shen, Manokhin and Xie 2019). That needs a score monotone in y: the candidate's score rising with y, and each unit's falling as the candidate rises. The full step's natural score is the in-sample predictive CDF. A unit scores F̃ᵢ = Σ_s w̃_s(y) F(y_i | x_i, θ_s), and the candidate F̃_c = Σ_s w̃_s(y) F(y | x_t, θ_s), under the add-one-in weights w̃_s(y) ∝ W_s p(y | x_t, θ_s). (F̃, since α on this page is the miscoverage level.) The brief made adopting it conditional on its staying monotone, so slice 4 read that first.

**Why it can fall.** The candidate's score has slope dF̃_c/dy = Σ_s w̃_s p(y | θ_s) + Cov_w̃(∂_y ln p(y | θ), F(y | θ)). The first term is the reweighted density, positive. The second is negative where the draws differ in location: a draw centred higher gains weight as y rises, and has the lower CDF at y. It outruns the first once the draws' predictives at x separate. Two equal draws N(±a, 1) give φ(a) − a(Φ(a) − ½) at their midpoint, negative past a ≈ 0.877, so locations 1.75σ apart already make the score fall. A unit's score moves with the candidate's too, by the covariance of its own CDF with the candidate's log-density slope. It falls where the posterior ties the unit's surface to the candidate's, and rises where it opposes them.

**What it reads.** Each reading sweeps the candidate along a grid, step 0.05 on the conjugate model and 0.01 on BASS, and compares the scores on the tail that holds their digits. The conjugate model is slice 0's, under the truth it does not share (n 19, m 20). BASS is slice 1's contaminated data (n 99, m 100), with five fits of ten test cases per row.

| Draws | Cases | The candidate's score falls | π's counts fall | Units whose score rises with the candidate |
|---|---|---|---|---|
| Conjugate, the exact posterior at S 50, 400 and 4,000 | 100 each | none | none | 961 of 1,900 |
| Conjugate, a chain stopped early (40 Metropolis steps from β = 5) | 100 | none | none | 961 of 1,900 |
| Conjugate, two chains stuck apart (the exact draws shifted by ∓1 in turn) | 100 | on 44 cases, by up to 0.034 a step, from F̃ 0.52 | on 53 cases, by up to 4 of 20 | 961 of 1,900 |
| BASS, a healthy Gaussian fit | 50 | none | none | 4,816 of 4,950 |
| BASS, a Student t fit | 50 | none | on 1 case, once, by 1 of 100 | 4,950 of 4,950 |
| BASS, a wrong fit (it never sees x₀) | 50 | none | none | 4,412 of 4,950 |
| BASS, a chain stopped early (one chain, no burn-in, 20 draws) | 50 | on 3 cases, by up to 0.0046 a step, from F̃ 0.66 | on 4 cases, by up to 3 of 100 | 4,888 of 4,950 |

**What it means.**
- *Not adopted.* The score is monotone where the draws sit in one mode, and falls where their predictives at x separate: two chains in different places, or a chain stopped early that still holds its first structures. Those are the cases the layer's guarantee without the sampler is for, so the brief's condition fails, and the slice stops at the reading.
- *Vovk's second condition fails almost everywhere.* Half the conjugate model's units and 89 to 100% of BASS's rise with the candidate. Yet π stayed monotone on draws in one mode, because those units rise more slowly than the candidate. Monotonicity there is a reading, not a theorem.
- *Its PIT would still be valid.* A conformal p-value needs no monotone score: π at the true y, ties broken at random, is uniform under exchangeability whatever the score. What a falling score loses is the CDF, and with it the quantiles, CRPS and a distribution to report.
- *A probability rounds the upper tail.* Held as F, a score within 10⁻¹⁶ of 1 rounds to 1, so units rising toward it tie with the candidate. π's counts fell there on 4 to 54 of the 100 conjugate cases in a first reading. Held as ln F and compared tail by tail, those falls vanish, which is why the seams carry the log.

**The alternatives,** as they were put to the user:
1. *The split form's predictive distribution*, Vovk et al.'s split conformal predictive system (Vovk, Petej, Nouretdinov, Manokhin and Gammerman 2020). The score is the training fit's own predictive CDF, F̂(y | x) = Σ_s W_s F(y | x, θ_s). The calibration scores do not move with the candidate, and the candidate's is a CDF in y, so both of Vovk's conditions hold by construction, and CB-2 makes it valid for any fit. It recalibrates the fit's predictive by its calibration PITs. It costs slice 3's split, plus one CDF read per draw and grid value per case, with no product over the units. This is the recommendation.
2. *The full step's score with a per-case check.* It reads π as a CDF where its counts rise on the grid and flags the cases where they fall; its PIT is valid regardless. It costs a product per grid like the density score's. It fails where the sampler is in doubt, the case the full step's guarantee does not cover anyway.
3. *A monotone envelope* of the full step's π, its running maximum. It is always a CDF, but its PIT is no longer uniform, since it moves mass up in one tail. Not recommended.
4. *Cross-conformal predictive distributions* (Vovk, Nouretdinov, Manokhin and Gammerman 2018) average K split distributions, so every unit serves. Their validity is approximate, at K fits.

**Decided on 2026-11-08:** the first as the default and the second as an opt-in, built in the next section.

**The log CDF seams,** which both forms read. `PointwiseLogDensity.LogCdf` and `CandidateLogDensity.LogCdf` carry each draw's ln F beside its ln p. BASS fills them on `cdf: true` (`BassFit.PointwiseLogDensity`, `BassFit.CandidateLogDensity`) from the law at the draw, `IBassNoiseLaw.At`, or N(f, σ²).
- *The tails.* The CDF is a closed form for the t and the asymmetric Laplace. For the Normal-Wald, the skewed t and the generalized hyperbolic it is one table of tails per draw, held to its 10⁻⁵ relative. Above each draw's median the seam reads ln F as ln(1 − S) from the law's survival. The mixtures' tables end at their 10⁻¹² tails, and beyond the upper one their `LogCdf` integrates the lower side, losing 1 − F's digits and reading up to 3·10⁻¹⁵ above 0. No noise law changed.

**Checks** (`ConformalDistributionTests`, `BassConformalDistributionTests`; the full readings with `BASS_QUALIFICATION=1`).
- *Units and candidates read one distribution* (CB-4). At the training cases the units' log CDF equals the candidates', bit for bit, under all six laws (120 draws by 80 units each).
- *Each draw's CDF is its density's integral.* At a new case, the rise of F over a grid of step 0.0005 matches Simpson's rule on the same draw's density, for one chain's 20 draws. The largest gaps are 1.5·10⁻¹⁴ (Gaussian), 2.2·10⁻¹³ (t), 7.2·10⁻⁷ (asymmetric Laplace, across its kink), 2.3·10⁻⁷ (Normal-Wald), 2.7·10⁻⁷ (generalized hyperbolic) and 3.2·10⁻⁶ (skewed t, its table).
- *The upper tail.* From 0.5 to 32 scales above the median, −expm1(ln F) gives back the law's survival to 10⁻¹² relative, ln F at or below 0, for the five laws.
- *Refusals:* a unit of two observations, the collapsed content, a matrix not of the density's shape, a value above 0 or NaN, and a weighted fit's units.
- *The reading's own pieces:* the tail comparison holds the digits a probability loses near 1, and the two-draw slope matches its closed form to 10⁻⁶ on both sides of a ≈ 0.877.
- *The quick reading* asserts the table's split on 20 conjugate cases a kind: the exact posterior and the chain stopped early never fall, two chains apart do.

## The predictive distribution, built (slice 4)

**The decision (2026-11-08).** On the reading above, the user chose the split form's CDF as the score ("Split-form CDF"). They then added the full step as an option beside it ("Yes, add it"). So the split form's predictive distribution is the default, and the full step's in-sample CDF is an opt-in that flags the cases where it falls. The full step stays an option because its trade-off depends on the regime, not because its design is flawed. Where the draws sit in one mode its score held monotone, and it spends no data on calibration. Where they separate it gives no distribution.

**What a case reports** (`ConformalDistribution`). Each test case's candidate y is ranked by a CDF score among its stratum's m units: Q(y, τ) = (B + τ(T + 1))/m, with B the units below the candidate and T those tied with it.
- *The distribution.* Q⁻ = B/m (`Lower`) and Q⁺ = (B + T + 1)/m (`Upper`, which is π(y) under the CDF score) bound it on the grid. `Mid` is their mean over τ.
- *The PIT.* Q(y, τ) at the true y with τ uniform (`At`) is uniform whatever the score (`coverage_tieBreak`, `rank_uniform`).
- *The set.* The equi-tailed set at α takes k = ⌊α·m/2⌋ in each tail, ties counted for the candidate in both (`Accepts`, `Intervals`, `Width`). It covers at least 1 − 2k/m, exactly with distinct scores (CB-1 per tail), and it is the whole line when k = 0.
- *The quantiles.* `Quantile(p)` is the first grid value where Q's mean over τ reaches p, read to the grid's step. It is −∞ when the first grid value already reaches p, and +∞ when none does.
- *The CRPS.* `Crps(y)` is E_τ ∫ (Q(z, τ) − 1{z ≥ y})² dz = ∫ (Q̄ − 1{z ≥ y})² + (Q⁺ − Q⁻)²/12 dz. Q̄ and the band are read linearly between grid values, and each piece's square is integrated exactly. The tails are cut at the grid's ends: the mass below the grid sits at its first value, the mass above at its last. Reading Q̄ linearly moves each of its m − 1 steps within a cell, which costs at most h(1 + 1/(3m)) in all.
- *A case with no distribution.* Where Q falls along the grid, `Fall` (in units of 1/m) and `ScoreFall` (the candidate's own score) say how far. Such a case keeps its p-values and its set, and its quantiles and CRPS are NaN.

**The split form, the default** (`SplitConformal.Distribution`; Vovk, Petej, Nouretdinov, Manokhin and Gammerman 2020).
- *The score* is the training fit's predictive CDF, F̂(y | x) = Σ_s W̄_s F(y | x, θ_s), read from the seams' log CDF by both tails.
- *Monotone by construction.* The calibration's scores do not move with the candidate, and the candidate's rises with y, so both of Vovk's conditions hold and every case has a distribution.
- *What it is.* The fit's predictive, recalibrated by its calibration PITs. The distribution's steps sit at F̂⁻¹ of the calibration PITs' order statistics, so a fit whose PITs pile into the tails, as a wrong one's do, is widened there.
- *Validity* is CB-2's for any fit, PIT included.
- *Cost:* slice 3's split, then per case one CDF read per draw and grid value, with no product over the units.

**The full step, an opt-in** (`ConformalBayes.Distribution`).
- *The score* is the in-sample predictive CDF under the add-one-in weights, F̃ᵢ for each unit and F̃_c for the candidate: two (n × S)·(S × G) products a batch, one per tail. The candidate's own score is summed as the products' entries are, by fused multiply-adds in draw order, so a unit whose column equals the candidate's ties with it exactly.
- *The flag.* Each case is read along its grid. One whose counts fall is flagged, with how far, and keeps its set but has no distribution. The flag is a readout in its own right: a falling score marks an x where the draws' predictives disagree, a mixing signal. π can fall while the candidate's own score rises, where a unit's score rises faster and overtakes it, so `Fall` and `ScoreFall` are read apart.
- *Its PIT* holds on every case, flagged or not, up to the draws' Monte Carlo error, as the full step's sets do, for draws from the posterior. Draws that are not, a wrong fit's or a chain stopped early, break the full step's symmetry, and the split form is the guarantee there. No margin is read for this score.

**Choices,** with what each changes:
- *Both tails held* (the seams' log CDF): a probability would round F to 1 in the upper tail and make units tie with the candidate there.
- *The equi-tailed set* is the distribution's own. The highest-density set stays slice 3's density score.
- *Quantiles and set ends to the grid's step.* The split form's ends have a closed form, F̂⁻¹(· | x) at the calibration PITs' k-th and (n + 1 − k)-th order statistics. It would need the engine's quantile function, so it is not built: it would remove the step's error at two quantile solves per stratum and level.
- *The tails cut at the grid's ends* for the CRPS. Q itself never reaches 0 or 1 (Q⁺ ≥ 1/m, Q⁻ ≤ n/m), so its CRPS over the whole line diverges, and a grid that holds the mass gives the integral that matters.

**Checks** (`ConformalPredictiveDistributionTests`; `BassConformalDistributionTests`, the measured default and the cost row with `BASS_QUALIFICATION=1`).
- *The counts against references* in seeded sweeps, with weighted draws, zero weights, chains, strata with no calibration unit, both tails far out and copied columns. The split form's counts equal the reference at all 1,968 grid values. The full step's equal it at all 3,122, and a quarter of its columns are copied from a unit and tie with it exactly.
- *The readouts* on counts given by hand: Q⁻, Q⁺, the set at k per tail, quantiles at and beyond the grid, and a case with no distribution. The CRPS matches a midpoint sum of 10⁶ pieces to 10⁻⁶.
- *The split form's PIT is uniform for any fit* (`coverage_tieBreak`, `split_coverage`), against the fits' own uncalibrated PITs:

| Fit | m | PITs | Split form's PIT: largest bin gap (se), KS | The fit's own PIT, KS |
|---|---|---|---|---|
| Conjugate, the exact posterior | 20 | 40,000 | 0.004 (0.0021), 0.005 | 0.143 |
| Conjugate, a wrong fit | 20 | 40,000 | 0.003 (0.0021), 0.005 | 0.363 |
| Conjugate, a chain stopped early | 20 | 40,000 | 0.003 (0.0021), 0.005 | 0.335 |
| BASS, a healthy fit | 200 | 2,000 | 0.015 (0.0082), 0.037 | 0.143 |
| BASS, the Student t fit | 200 | 2,000 | 0.022 (0.0082), 0.036 | 0.044 |
| BASS, a wrong fit | 200 | 2,000 | 0.014 (0.0082), 0.023 | 0.057 |
| BASS, a chain stopped early | 200 | 2,000 | 0.018 (0.0082), 0.020 | 0.056 |

- *Vovk's conditions hold by construction.* On every grid case of every fit, 100 a fit on the conjugate model and 200 on BASS, neither the score nor Q falls. Q⁻ is 0 below the grid and Q⁺ is 1 above it.
- *The full step's opt-in.* On the conjugate exact posterior (m 50, 2,000 cases) no case is flagged, and its PIT is uniform over every case and over the unflagged ones (KS 0.030). On two chains stuck apart, 287 of 500 cases are flagged, 267 of them with the candidate's own score falling. Each flagged case keeps its set, with NaN quantiles and CRPS.
- *Agreement with slice 3's sets.* The equi-tailed set is slice 3's split rule in each tail at α/2 on the CDF score, at every one of 61,127 grid values in the bulk. On a fit that is right (data from the model, m 2,000), its ends, the density set's and CQR's all meet the fit's own central 90% interval. The mean gaps run 0.001 to 0.020 over 100 cases, each within four of its standard errors and a grid step (0.02).
- *CRPS against a quadrature reference.* The exact step distribution, read at 16 times the grid's resolution, gives the reference. The grid's CRPS lies within 1.3·10⁻³ of it at m 20 and 4.6·10⁻⁴ at m 200 (means 4.9·10⁻⁴ and 1.3·10⁻⁴), against the bound h(1 + 1/(3m)) = 0.010.

**The measured default.** Both forms run on the same data. The split form fits the training split and calibrates on the rest; the full step fits all of it. On the conjugate model that is 30 + 19 units against 49, 500 replications of 10 test cases. On BASS it is 99 + 199 against 298, 20 fits of 100 test cases. Width and CRPS are read where both give a distribution, on 200 grid cases a fit. Cost is the layer's per case, the seams built before the clock, on the wall clock in one pass, so indicative. It was read again on slice 4's final code: the first reading ran before the review's cost fix and overstated both forms, the split form's most. Every other reading came out the same to its last digit.

| Fit | The full step flags | 90% width, split / full | CRPS, split / full | The full step's PIT, every case / unflagged, KS | Cost per case, split / full |
|---|---|---|---|---|---|
| Conjugate, the exact posterior | 0% | 6.22 / 5.97 | 0.902 / 0.814 | 0.019 / 0.019 | 1.4 / 4.8 ms |
| Conjugate, a wrong fit | 1.6% | 7.46 / 7.99 | 1.262 / 1.178 | 0.021 / 0.018, a bin 4.1 se out | 1.9 / 5.2 ms |
| Conjugate, a chain stopped early | 0% | 18.15 / 18.51 | 3.225 / 3.029 | 0.011 / 0.011 | 0.3 / 1.1 ms |
| Conjugate, two chains stuck apart | 66% | 4.65 / 5.35 (66 cases) | 0.688 / 0.649 | 0.171 / 0.028 | 1.3 / 4.7 ms |
| BASS, a healthy fit | 0.1% | 1.322 / 1.335 | 0.2044 / 0.1993 | 0.022 / 0.022 | 12 / 72 ms |
| BASS, the Student t fit | 0.1% | 1.303 / 1.317 | 0.1956 / 0.1884 | 0.018 / 0.018 | 19 / 97 ms |
| BASS, a wrong fit (it never sees x₀) | 0% | 2.478 / 2.448 | 0.4770 / 0.4632 | 0.018 / 0.018 | 20 / 118 ms |
| BASS, a chain stopped early (one chain, no burn-in, 20 draws) | 8.5% | 2.030 / 2.161 (179 cases) | 0.3448 / 0.3439 | 0.031 / 0.036, a bin 3.8 se out | 0.6 / 4.4 ms |

The conjugate costs are per grid of 401, and BASS's per grid of 901 at S 1,000.

- *The split form's PIT holds on every fit.* Its data cost is the calibration: 19 of 49 units here, 199 of 298 on BASS.
- *The full step is sharper where its draws are a posterior that mixed,* by more where the data are scarce: its CRPS is 10% lower at n 49 on the conjugate model and 2.5 to 3.7% lower at n 298 on BASS, where its width is the split form's. It flags almost no case there (0 to 0.1%). **That is the regime to opt in: draws a posterior that mixed, and data too scarce to spare a calibration split.**
- *The full step needs draws that are a posterior.* The two chains stuck apart are not one, and its PIT breaks over every case (KS 0.17). Two thirds of those cases are flagged, so the flag marks where it fails. The conjugate wrong fit is not a posterior either, and its PIT is off by a bin at 4.1 se. The conjugate chain stopped early is a narrow cluster: the add-one-in weights barely move it, so every unit is read alike and its PIT holds. BASS's chain stopped early is no posterior either: 8.5% of its cases are flagged, and the full step's PIT sits a bin 3.8 se out over every case. BASS's wrong fit is a posterior that mixed, of a model that never sees x₀, and the full step holds there (KS 0.018): the model may be wrong, but the draws must be its posterior.
- *So the default stays the split form,* as decided. The reading gives no reason to move it: the full step's gain is a few percent where it is valid, and where its premise fails the flag often catches it, though not always.

**Cost** at n 1,000 units and G 200, on the cycle clock, best of five; indicative. Every row was read in one run on slice 4's final code, with the density score's and the density full step's beside it. The density score reads as it did in slice 3 (2.1 and 6.7 ms against 2.1 and 7.3), so the machine was about as loaded as then.

| S | The split form: calibration PITs, once per stratum | The split form: per grid | Slice 3's density score: per grid | The full step's opt-in: per grid |
|---|---|---|---|---|
| 1,000 | 18 ms | 2.3 ms | 2.1 ms | 52 ms |
| 4,000 | 70 ms | 9.1 ms | 6.7 ms | 219 ms |

- *The split form reads two tails,* one exp a value and a log above the median, at 1.1 to 1.4 times the density score's cost in the same run, and needs no product.
- *The full step's opt-in costs two products a batch,* about a third of the density full step's (152 to 194 ms at S 1,000 and 538 to 599 ms at S 4,000 in the same run), which reads three products and the margin.

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

## The other engines' candidate densities (slice 5)

Each engine joins as BASS does: log p(y | x_t, θ_s) at new cases, in its pointwise density's rows with its chain lengths and latent report, and on `cdf: true` each draw's log CDF beside it. Which form takes a seam turns on two premises. The full step needs each draw's terms in sample, the units independent given the draw. Both forms need the candidate exchangeable with the units it is ranked against.

| Engine | Its candidate density | The full step | The split form |
|---|---|---|---|
| BAPS | `Baps.CandidateLogDensity(result, design, grids, strata, cdf)`: N(y; z_tβ, σ²) at the new cases' design rows | Takes it | Takes it |
| BARS | `BarsEnsemble.CandidateLogDensity(result, x, grids, strata, quadrature, cdf)`: each observation model's term at the draw's knots and coefficients, each observation's latent integrated as in training | Takes it | Takes it |
| The hierarchical P-spline | `HierarchicalPSpline.CandidateLogDensity(result, design, replicates, grids, strata, cdf)`: N(y; z_tβ_r, σ²), a new point on a seen curve, conditioned on its drawn curve | Its guarantee does not reach it | Takes it |
| BGP | `BgpSampler.CandidateLogDensity(result, xCases, grids, strata, cdf)`: the GP's predictive given the training data at each draw (t, η) | Refuses it | Takes it |
| The changepoint posterior | None | | |

- *BAPS.* A candidate is one observation; its design row is the fitted basis at the new case, evaluated as the training design was. The smoothing parameters enter only β's prior, so the term is the Gaussian at the draw's β and σ², as in training.
- *BARS*, under every observation model it scores in sample:
  - the Gaussian;
  - robust weights' marginals (Student t, the contaminated normal);
  - the variance-mean noises at the draw's σ² and skew (Student t, the asymmetric Laplace, GIG and GBP mixing, the last by quadrature on `quadrature: true`);
  - the Poisson mass;
  - the Whittle ordinate's exponential, whose candidate is an ordinate held out of the periodogram, the reading under which its units are exchangeable.

  A value outside the response's support takes −∞: a count that is not a whole number, a negative ordinate. The Gibbs-loss model is refused, as its pointwise density is: its likelihood is a pseudo-likelihood.
- *The hierarchical P-spline.* The candidate conditions on β_r, as a training unit that splits its replicate does, and a point on a new curve is refused: it would integrate that curve while the training units condition on theirs (CB-4).
  - Its full step lacks exchangeability. The model takes one design for every replicate, so each training point's input is held by every curve, and a new point's input by none. The candidate is not exchangeable with the training points, and no partition restores it.
  - The layer cannot tell this from the seam, whose latents match the training density's, so the docstring states it.
  - The split form holds: calibration and test points on the seen curves, each at an input drawn as the others are, are exchangeable with one another for any fit.
  - A whole new curve on the shared design is exchangeable with the training curves, but it is a vector: E17's whole-curve bands.
- *BGP.* Its content is the collapsed one: f, the constant mean and s² are integrated against their posterior given the training data, never drawn.
  - A known σ² gives N(f̂(x), v + σ²), with v the variance of f, the mean's uncertainty included. An unknown scale under 1/s² gives a Student t on n − 1 degrees of freedom with scale² Q(v + η)/(n − 1). Either is the evidence ratio ln p(Y, y | t, η) − ln p(Y | t, η) (`GpRegression.PredictiveLogDensities`).
  - The full step refuses it, as it refuses the training density: a draw (t, η) leaves the units dependent through f, so each term leaves its unit out (`LeavesUnitsOut`). The split form takes it.
  - Draws that repeat their predecessor's (t, η) reuse its row, one O(n³) refit per accepted state.
- *The changepoint posterior* has no seam. Its posterior is over segmentations of the n positions it was given, one observation at each. A candidate at another position, a held-out one included, falls outside every segmentation, and a second value at a seen position is not a unit of that sequence. No candidate is exchangeable with the units.

**The log CDFs.**
- *The candidates' and the training units'.* BAPS's and BARS's pointwise densities take `cdf` for the full step's predictive distribution, units of one observation each. The candidates' log CDF reads the same law, for both forms.
- *Each law's:*
  - the Gaussian's;
  - Student t's, closed, from both tails;
  - the contaminated normal's mixture of two, each tail summed in logs (`IWeightResampler.MarginalLogCdf`, a default member that says a resampler without one gives none);
  - exponential mixing's asymmetric Laplace, closed, its smaller rate taken as 2/(√(m² + 2s²) + |m|);
  - the GIG and GBP mixtures' own laws, tabled once per draw (`BarsNoise.LogCdfAt`);
  - the Whittle ordinate's, ln(1 − exp(−y·e^(−η))), η the draw's linear predictor there;
  - BGP's Gaussian or Student t.
- *Refused:* a count's CDF steps, so the Poisson model gives none. Its PIT is not uniform without a randomized tie, which the predictive distribution does not take.
- *Shared:* every engine's seam reads ln F through `NoiseLogCdf`, the split that reads the upper half from the survival, moved from BASS's seam to the shared seam file.

**Checks** (`EngineCandidateTests`; coverage and cost with `BASS_QUALIFICATION=1`).
- *At the training cases, the pointwise column, bit for bit.*
  - BAPS, and BARS under all nine of its models and noises: Gaussian, t weights, contaminated, t noise, quantile, skewed GIG, skewed GBP, Poisson and Whittle. The log CDF and the latent report match too, and the full step takes each pair.
  - The hierarchical P-spline's, against the column of a unit that splits its replicate. Beside whole-curve units, which integrate β_r, the full step refuses it.
- *BGP's, against its pointwise content:* a GP fit on the data less one case gives, at that case, the leave-one-out term of the fit on all of them. The largest relative gap over 36 cases is 1.7·10⁻¹⁴ with a known σ² and 1.6·10⁻¹⁴ with an unknown scale. The full step refuses the seam, and the split form takes it.
- *Each law's log CDF against its density's integral,* a seeded sweep over σ² from e⁻⁴ to e⁴: a slice of the body, and each tail where the density falls to 10⁻¹² of the body's.
  - The closed forms agree to 10⁻¹² or better in the body and 3·10⁻¹⁰ in the tails.
  - The tabled mixtures agree within their tables' stated 10⁻⁵: 6.8·10⁻⁶ at worst in the body (GBP without a drift), 6·10⁻⁷ in the tails.
- *Held-out coverage on data from outside every model.*
  - The truth is f = sin(2πx) + x/2, with noise 0.2ζ plus, one time in seven, 0.8·Exp(1).
  - There are 20 fits each, with m 100 and k 10, nominal 0.900. The full step's fits have n 99. The split form's have 99 calibration units, beside a training fit of 99 for the hierarchical P-spline (three curves on a shared design of 33) and 60 for BGP.
  - There are 50 held-out cases per fit, and the se is the fits' spread. The width is read on 20 cases per fit, on a grid of step 0.08:

| Engine and form | Coverage (se) | Mean width | Sets of two intervals |
|---|---|---|---|
| BAPS, the full step | 0.913 (0.012) | 1.11 | 0 of 400 |
| BARS (Gaussian), the full step | 0.917 (0.012) | 1.10 | 1 of 400 |
| The hierarchical P-spline, the split form | 0.902 (0.012) | 1.23 | 0 of 400 |
| BGP (unknown scale), the split form | 0.903 (0.014) | 1.12 | 0 of 400 |

  - BARS's one split set is at x = 0.002, the domain's edge: [−3.00, −2.12] ∪ [−0.52, 0.68], its lower piece reaching the grid's end. Its add-one-in weights fell to an ESS of 2.2, and the inflated-set gap, 0.82, flags it. That is CB-3's reading at work, not a miss: the gate requires every set that reaches the grid's end to carry the gap's flag.
- *Its cost,* for one grid of 200 values, S = 2 chains × 500, indicative on the cycle clock. The layer's step is the full step for BAPS and BARS and the split form's for the others, on a calibration of 99:

| Engine | n | The candidate density per grid | Per candidate | The layer's step per grid |
|---|---|---|---|---|
| BAPS | 1,000 | 2.1 ms | 10 µs | 204 ms, the full step |
| BARS (Gaussian) | 1,000 | 3.3 ms | 17 µs | 196 ms, the full step |
| The hierarchical P-spline (10 curves of 100) | 1,000 | 1.9 ms | 10 µs | 3.6 ms, the split form's |
| BGP (unknown scale) | 400 | 20.1 s per call, whatever its cases | | 3.2 ms, the split form's |

  - The full step's 200 ms is its default step at S 1,000, as the margin's opt-in report measured it, the training side included.
  - BGP's seam is its refits: one O(n³) refit per accepted state, as its pointwise density pays, shared by every case in the call. The run already holds each draw's variance at its own test points, a cheaper path recorded as a follow-up in the brief.

**Found and handed on, not fixed here.**
- *The defect.* Under GBP mixing with a drift, `NormalVarianceMeanMixture`'s table reads ln F about 10⁻³ off near x = 0, while the density there is right to 12 digits.
  - At m/s = −0.2 and GBP(1, 1, 1.5), ln F(0) is −0.549285 against the integral's −0.550321.
  - The table's slope runs 4% low across ±10⁻³·s, and the gap falls to 5·10⁻⁶ by 0.5·s.
  - Without a drift F(0) is ½ exactly, and the GIG mixtures are within their 10⁻⁵.
- *Its effect.* BARS's GBP noise reads its log CDF from that table, so its predictive distribution inherits the error near the median. Its density, which the sets read, does not.
- *Where it stands.* The file belongs to the noise laws' owner. `TheLogCdf_UnderGbpMixingWithADrift_IsItsDensitysIntegral` is skipped until the fix, with the correct references as its assertions.

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

## The per-unit bounds (between slices 4 and 5)

Built on 2026-11-08 at the user's order, from the Lean chip's per-unit CB-3 (80a1993): `coverage_inflated` takes bounds by ordering and position, `coverage_close` bounds carried by the units, and `fragileEstimated_mono` lets any dominating bounds serve.

**What changed.**
- *The margin reads carried bounds* when asked for: each unit's largest bound over every ordering the draws represent, the candidate's the larger of its own and its bounds under each unit's weights. The uniform bound, the other option, is the largest carried one. The margin itself became opt-in the same day (above).
- *The inflated rule keeps per-position bounds by default,* with a uniform option for comparison.
- *Slice 2's two margin bounds are replaced.* Its carried bound read two orderings, the held one and the unit's own turn; its uniform bound read the held ordering alone, so it could sit below a unit's own-turn bound. Neither dominated the units' errors over the orderings.

**Why two orderings were not enough** (`ConformalPerUnitBoundsTests`, a reading with `BASS_QUALIFICATION=1`). On the conjugate fixture (n 20, 20 data sets, 161 grid values each), every unit's error was read under every unit's weights from public steps. Some other unit's turn as the candidate gave a larger error than slice 2's bound at most (unit, grid value) pairs:

| S | A unit in a tail | Pairs above slice 2's bound | By, median (largest) | m̂: slice 2's | Over every ordering | Uniform |
|---|---|---|---|---|---|---|
| 40 | no | 77% | 1.34 (11.6) | 3.60 | 4.10 | 13.30 |
| 640 | no | 64% | 1.75 (38.4) | 2.80 | 3.50 | 10.90 |
| 40 | yes | 84% | 1.43 (8.7) | 3.75 | 4.25 | 14.15 |
| 640 | yes | 87% | 2.89 (67.9) | 3.00 | 3.85 | 10.25 |

A unit far out, or of high leverage, makes the weights degenerate when it is the candidate, and that widens every unit's error. So a unit's worst error comes mostly from another unit's turn, not its own. The band moves by 0.5 to 0.85 of 21.

**Checks** (`ConformalMarginTests`, `ConformalPerUnitBoundsTests`).
- *The layer against a reference,* 395 candidates over 80 seeded fixtures with weights, chains and far ranges. Each unit's bound over the orderings and the candidate's carried bound match the same errors read from public steps, each unit's column taken as the candidate's, to 10⁻⁷ relative. The carried and uniform bands are the literal `FragileEstimated` counts on them, and both inflated ranks match their definitions.
- *No set wider than the uniform rule's.* At every grid value of every fixture the per-position inflated rank is at most the uniform one, and the carried band at most the uniform band.
- *Coverage at 1 − α* (exact draws, S 40, n 19, k 2, 3,850 true values, se 0.0048):

| Set | Coverage | Width | The margin's reading at the true value |
|---|---|---|---|
| Plain | 0.905 | 5.09 | carried 0.732, uniform 0.358 (slice 2's uniform 0.485) |
| Inflated, per position | 0.935, against 0.89 | 5.72 | |
| Inflated, uniform | 0.945, against 0.89 | 6.19 | |

- *A unit in a tail* (n 20, the last unit at x 3.5 with a residual of 6): the carried band reads 4.92 of 21 against the uniform 16.32 at S 40, and 3.38 against 8.10 at S 640. The tail unit's bound is the widest in 37 and 32 of 40 data sets, 10 to 12 times the median unit's, so the carried margin gains. The plain set covers 0.950 and 0.943 over 2,000 true values, above the carried readings of 0.790 and 0.843.
- *No unit standing apart* (a Beta–Bernoulli posterior, ten ones and ten zeros, symmetric about ½, S 4,000): every carried bound equals the uniform one to 10⁻⁹, and the bands are equal at both candidates and every α (10, 10 and 0 at α 0.1, 0.3 and 0.5). The carried margin falls back to the uniform one.

**Slice 2's readings, side by side** (the same seeds; `BASS_QUALIFICATION=1`). Plain coverage, the per-position inflated set and the gap read exactly as slice 2's did; the margin reads differently, and the uniform inflated rule is new.

| Exact draws, m 20 | S 8 | S 40 | S 160 | S 2,560 |
|---|---|---|---|---|
| m̂ carried: slice 2 → now | – → 7.44 | 4.6 → 5.04 | – → 3.73 | 2.5 → 2.57 |
| m̂ uniform: slice 2 → now | – → 17.10 | 12.8 → 14.08 | – → 11.60 | 7.7 → 8.44 |
| Inflated coverage, per position / uniform (se 0.0067) | 0.941 / 0.955 | 0.937 / 0.946 | 0.927 / 0.936 | 0.912 / 0.916 |
| Inflated width, per position / uniform | 6.15 / 6.71 | 5.85 / 6.28 | 5.68 / 5.95 | 5.47 / 5.60 |

| BASS, n 99, k 10, S 1,000 | Gaussian | Student t |
|---|---|---|
| m̂ at the true values, median, carried: slice 2 → now | 72 → 97 | 52 → 96 |
| m̂ uniform: slice 2 → now | 98 → 100 | 94 → 98 |
| Coverage: plain / inflated per position / inflated uniform | 0.905 / 0.958 / 0.986 | 0.894 / 0.937 / 0.945 |
| Width: plain / inflated per position / inflated uniform | 1.20 / 2.07 / 3.44 | 1.06 / 1.49 / 1.83 |

- *The per-position inflated rule pays most on BASS.* The uniform bound, a median 4.2 nats on the Gaussian fits against each candidate's 0.08, widens the inflated set by two thirds there for 3 points of coverage the per-position set does not need.
- *The carried margin pays on exact draws, about a third of the uniform band, and little on BASS,* where the units' errors over the orderings are wide for every unit. Slice 2's carried reading there, 72 and 52, came from bounds that did not dominate. The margin stays a reported worst case.

**Cost** (n 1,000, G 200, cycle clock, best of five, indicative). The units' orderings cost 331 to 447 ms once per stratum at S 1,000 and 1,137 to 1,795 ms at S 4,000, against 16 and 65 to 78 ms for the stratum's likelihoods; with two chains, each unit's r_eff adds 35 to 58%. Per grid the fourth product, the candidates' errors under the units' weights, brings the products to 70 ms at S 1,000 and 261 ms at S 4,000. On a single grid the step costs 522 to 662 ms and 1,717 to 2,400 ms with the orderings built; the orderings are paid once per stratum, so many cases share them.

**The margin, opt-in** (2026-11-08; `ConformalMarginTests`).
- *Both states.* Over 60 seeded fixtures with weights, chains, far ranges and a candidate outside every draw's support, under both inflation bounds, 240 sets were read without the margin and with it under each bound. Without it, `Computed` is false, the band, count and carried bounds are null and the plain set's coverage under error is NaN. With it, they are read. Every other reading is the same bit for bit: p-values, counts, scores, ESS, errors, inflated ranks, k̂, the gap, the split-half mass, the fallback signal and the widths.
- *The margin's own readings,* asked for, are those above, unchanged.
- *Cost* (n 1,000, G 200, the full step on one grid, on a busier machine than the row above): 170 to 193 ms by default against 728 to 803 ms with the margin at S 1,000, and 579 to 864 ms against 1,550 to 3,222 ms at S 4,000. The margin's orderings were 0.42 to 0.46 s and 1.1 to 2.0 s of that.

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
| The seams | `PointwiseLogDensity` (training); `CandidateLogDensity(values, grids, latents, stratumOf, chainLengths)` with `Weights`, `Offset(t)`, `Cases`, `Draws`; `LogCdf` on both, each draw's ln F |
| BASS's candidate density | `BassFit.CandidateLogDensity(data, cases, grids, strata, unseenLevels, quadrature, cdf)`; `BassFit.PointwiseLogDensity(..., cdf)`; `BassCandidates.UnitLogCdf` and `BassPointwise.Latents` (internal, shared) |
| The other engines' (slice 5) | `Baps.CandidateLogDensity(result, design, grids, strata, cdf)` and `Baps.PointwiseLogDensity(..., cdf)`; `BarsEnsemble.CandidateLogDensity(result, x, grids, strata, quadrature, cdf)` and `BarsEnsemble.PointwiseLogDensity(..., cdf)`; `HierarchicalPSpline.CandidateLogDensity(result, design, replicates, grids, strata, cdf)`, for the split form; `BgpSampler.CandidateLogDensity(result, xCases, grids, strata, cdf)` and `GpRegression.PredictiveLogDensities(fit, xTest, grids)`, for the split form; `IWeightResampler.MarginalLogCdf`; `NoiseLogCdf`, `CandidateSeam` and `BarsNoise.LogCdfAt` (internal, shared) |
| The full step | `ConformalBayes.FullStep(training, candidates, ConformalOptions { Alpha, ErrorLevel, Inflate, FallbackShare (opt-in, null by default), MarginBounds (opt-in, null by default), InflationBounds })` → `ConformalPrediction(Sets, Alpha)`; `ConformalBayes.Scores(training, candidates, column)` → `ConformalScores(Units, Scores, Candidate, EffectiveSampleSize)` |
| A case's set | `ConformalSet`: `Stratum`, `Units` (m), `Level` (k), `Coverage`, `WholeSpace`, `Grid`, `PValues`, `Below`, `Tied`, `CandidateScores`, `EffectiveSampleSizes`, `Accepts(g)`, `Intervals` (`GridInterval(Lower, Upper)`), `Resolution`, `OpenBelow`, `OpenAbove`, `Width` |
| CB-3's reading | `ConformalSet.Gap`, `SplitHalfMass`, `SplitHalfDisagreements`, `AcceptsPlain(g)`, `AcceptsInflated(g)`; `ConformalSet.Margin` → `ConformalMargin` (`Multiple`, `CandidateErrors`, `UniformErrors`, `CandidateCarriedErrors`, `CarriedUniformErrors`, `UnitCarriedErrors`, `Fragile`, `InflatedRanks`, `ParetoK`, `ParetoKThreshold`, `Bounds`, `Inflation`, `Computed`, `Count`, `Untrusted(k)`); `ConformalSet.Inflated`, `ErrorLevel`, `CoverageUnderError`, `FallBack`; `MarginBounds` (`Uniform`, `Carried`); `InflationBounds` (`PerPosition`, the default, `Uniform`); `ConformalScores.Errors`, `CandidateError` |
| The split form | `SplitConformal.Density(calibration, candidates, alpha)` → `ConformalPrediction(Sets, Alpha, Form)`; `SplitConformal.Quantiles(QuantileReadings(Lower, Upper, StratumOf), response, test, alpha)` and `SplitConformal.Residuals(ResidualReadings(Center, Scale, StratumOf), response, test, alpha)` → `ConformalIntervals(Intervals, Alpha, Calibration)`; `ConformalInterval`: `Units`, `Level`, `Threshold`, `Lower`, `Upper`, `WholeSpace`, `Coverage`, `Empty`, `Width`, `Contains(y)`; `SplitConformal.Calibrate(scores, stratumOf, alpha)` → `SplitCalibration`: `Units(h)`, `Level(h)`, `Threshold(h)`, `Rank(h, s)`, `PValue(h, s)`, `Accepts(h, s)`; `ConformalForm` (`Full`, `Split`) |
| The predictive distribution | `SplitConformal.Distribution(calibration, candidates)` (the default) and `ConformalBayes.Distribution(training, candidates)` (the opt-in) → `ConformalDistributions(Distributions, Form)` with `Flagged`; `ConformalDistribution`: `Units`, `Grid`, `Below`, `Tied`, `Form`, `Fall`, `ScoreFall`, `Monotone`, `Lower(g)`, `Upper(g)`, `At(g, τ)`, `Mid(g)`, `Level(α)`, `Coverage(α)`, `WholeSpace(α)`, `Accepts(g, α)`, `Intervals(α)`, `Width(α)`, `OpenBelow(α)`, `OpenAbove(α)`, `Quantile(p)`, `Crps(y)` |
| The product | `DenseProduct.Multiply`, `Pack`, `FusedDot`; `PackedColumns` |

## Not built yet

- A margin whose bounds read each ordering's own posterior, rather than the held draws reweighted, which needs a run per ordering.
- CQR's asymmetric form, cross-conformal and the jackknife+ (above).
- The split form's set ends and quantiles in closed form, F̂⁻¹(· | x) at the calibration PITs' order statistics, which needs each engine's quantile function; and cross-conformal predictive distributions (slice 4's alternatives).
- BASS's candidate seam for a direct functional fit's curves, and a whole new curve for the hierarchical P-spline (E17).
- A BGP content with f drawn at the training and candidate inputs per retained state, which the full step would take (slice 5's follow-up).
- The changepoint posterior's seam, which no candidate exchangeable with its units admits (slice 5).
- CB-4's latent-posterior reading, and Hoff's leave-one-out score.
- Refining a set's boundaries between grid values.

**See also:** [BayesianRegression](./BayesianRegression.md) (the pointwise seam) · [VarianceMeanMixtures](./VarianceMeanMixtures.md) (the noise laws whose latents are integrated) · [LinearAlgebra](./LinearAlgebra.md) (the product)
