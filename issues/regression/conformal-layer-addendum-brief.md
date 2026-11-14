# The conformal layer: addendum (chip brief)

Opened 2026-11-14 at the user's word ("that brief was getting way too long, addendum brief"). The [conformal-layer brief](archive/bass-gbass-khaos-activespaces/briefs/conformal-layer-brief.md), now archived, is the record of E16 through slice 15 and slice 14's follow-ups: its decisions 1 to 12, the plan for slices 7 to 18, the gates common to every slice, the working rules and every report. They all stand. This addendum carries what is open and the reports from here on.

**Chip:** the conformal chip, `local_69bae1e5`. **Verified by:** the main session, `local_b6732445`, from clean archives (`scripts/verify-csharp.ps1`, `scripts/verify-qualification.ps1` for gated tests).

## Where E16 stands

- Slices 0 to 15 are built and verified; slice 14's follow-ups 1 to 3 and its ratio estimator too, on their Lean (decision 12).
- Follow-up 2 is completed by `3f68330`, the merged K-fold CV+ bound under likelihood ratios, on the Lean at `47acb3c`.
- Slice 16, batch predictive inference (`60ba2cb`), slice 17, scores beyond the model's own (`0b91db9`), and slice 18, the guarantee read out (`5b55a6f`), are built and verified. Their reports are below. Every slice in the plan is built.
- The defects chip owns `PredictiveChecks.cs` and `Pointwise.cs`'s weight check; a change there is agreed with it first (the archived brief's coordination note).

## Open, each with the decision it would change

- **Slice 17's three follow-ups** (its report): RoBAS-Full, the horseshoe score; CBMA's evidence from the engines' own fits; a grid-free interval for the residual scores.
- **Slice 18's three follow-ups** (its report): a smoothed cell law for the log score, which is +∞ in a cell the step-function predictive gives no mass; the Lean's calibration at n ≥ 2; a set's own miscoverage bound.
- **Slice 15's precision rule** reads the training split's ESS: with 1,000 training and 399 calibration points it chose m of 294 to 500, and 7% of the sets were the whole space. It decides whether the fit takes the calibration split's size and scales the target itself.
- **Slice 15's alternation speed** along the factorization's flat directions. It decides the cost at large n.
- **τ through richer features,** the ratio estimator's open follow-up. It decides whether the layer offers a nonparametric τ.
- **Two follow-ups agreed with the GBASS chip** for E21's slice 5 (`127d64d`): new-group coverage by group-level units, and the layer's own check that strata lie inside cells.
- **Slice 12a's k̂ diagnostic,** opt-in or on by default (the handoff's waiting list).

## Working rules

As in the archived brief: each slice is reported here, then waits for the main session's verification before the next starts; anything that changes a default or removes an option goes to the user; chips leave `docs/` alone and each report ends with **Doc changes owed**.

## Reports

Newest first.

### 2026-11-14 — Slice 18's follow-up: the smoothed cell law for the log score

Built at the user's pick (2026-11-14, "Smoothed cell law for log score"), the first of slice 18's follow-ups.

**Landed** on `main`, nothing pushed:
- `b4b7e22`:
  - `ConformalDistribution.Smoothed` (`SmoothedCdf`: G̃ and its complement at each grid value);
  - `CellLaw` (`Steps`, `Smoothed`), which `CellMasses`, `LogScore` and `QuadraticScore` now take;
  - `SplitPredictiveDistribution.SmoothedAt`;
  - `SmoothedCellLawTests`.
- *Nothing existing changes:* `CellLaw.Steps` is the default, so every existing call reads as before.

**What it computes.**
- *The law.* G̃(F̂(y | x)), G̃ the calibration PITs' empirical CDF interpolated linearly between its knots (0, 0), (PIT₍ᵢ₎, i/(n + 1)) and (1, 1). At a tie with the PITs it is (B + T)/(n + 1).
- *Where it sits.* It lies inside the conformal band [Q⁻, Q⁺] at every y: it is the randomized Q(y, τ) with τ moving linearly across each gap between neighbouring PITs.
- *What it fixes.* Wherever the fit's CDF rises, every cell has mass, so the log score is finite: `logLoss_proper`'s positivity hypothesis, which the step law breaks, is met.
- *Precision.* Each point is read on the tail that holds its digits, as the CDF scores are, its complement kept beside it.
- *What it gives up.* Exact uniformity: its PIT at the truth lies within 1/(n + 1) of uniform. The randomized Q keeps the exact guarantee; this law is for scoring.
- *Where it applies.* The split form's unweighted cases. The full step's score is not monotone in y, and a weighed case reads weighted ranks, so both have none.

**Gates.**
- *Band and mass* (49 calibration units, 20 cases, cells of width 0.1): G̃ lies in [Q⁻, Q⁺] and rises at every grid value, G̃ + (1 − G̃) = 1, and every cell has mass. In every case the step law leaves cells empty.
- *The log score* (199 calibration units, 400 cases, cells of width 1): the smoothed score is finite everywhere and ranks the right fit ahead of one ignoring x, 1.447 against 1.523. In this draw the step law happened to be finite everywhere too; its +∞, 2 of 400 in slice 18's draw, is the band-and-mass gate's empty cells.
- *Uniformity* (99 calibration units, 4,000 replications): the smoothed PIT departs from uniform by at most 0.0068 at the deciles, within 1/(n + 1) = 0.0100.
- *The grid-free readout:* `SmoothedAt` matches the grid's G̃ bit for bit.
- *The full step and a weighed case* have no smoothed law, and their scores on it are NaN. An undefined `CellLaw` is refused.
- *Existing classes,* slice 18's and the exact-ends tests, pass unchanged.
- *Per commit,* the selected classes and the pinned tier: 1,228 of 1,228.

**Defects found: none.**

**Doc changes owed:** the guarantee-read-out section's cell scores gain the smoothed law, with its band property, its 1/(n + 1) uniformity, where it applies, and the API rows for `SmoothedCdf`, `CellLaw` and `SmoothedAt`.

### 2026-11-14 — Slice 17, scores beyond the model's own

Built at the user's pick (2026-11-14, "Slice 17, scores beyond the model"), the archived plan's row 17, the last in the plan. Sources read in TeX:
- Bhagwat, Kong and Jiang, CBMA (arXiv 2511.16924v1), §3, Algorithm 1, the lemma and both theorems;
- Ashouritaklimi, Cortinovis and Caron, RoBAS (arXiv 2607.04236v2), §2.2–§3.4;
- Wu, Lohmanova, Kaski and Caprio (arXiv 2602.03331v2), §4 and Algorithm 1.

Each is a symmetric score or a rule already proved, so coverage stands on CB-1, and for Wu et al. on slice 18's `fail_le_fixed`. No new Lean was needed. The plan's Lean column asked for Wu et al.'s risk bound "if it differs from coverage": it is slice 18's high-probability rule, already proved and verified.

**Landed** on `main`, nothing pushed:
- `0b91db9`:
  - `ModelAveraging.FullStep`, giving an `AveragedPrediction`;
  - `ResidualShrinkage.Prediction` with `ResidualScore` (`DistanceToOrigin`, `DistanceToAverage`, `Shrinkage`);
  - `RiskControlledSets.Density`, giving `HighDensitySets`;
  - `ScoresBeyondTheModelTests`.
- *Nothing existing changes.*
- *Verified by the main session on 2026-11-14,* from a clean archive of `ca4ff81` with zero skips: the Commit selection, 48 of 48 ([receipts](receipts/verify-conformal-s17-0b91db9/)).

**What it computes.**
- *CBMA* (`ModelAveraging`).
  - From K models' full-step seams over the same units, strata and grids, and each model's log posterior probability given the n units (its log evidence plus log prior), the score is σᵢ = Σ_k q_k(y)·σᵢ^{M_k}, q_k(y) ∝ p(M_k | Z₁:ₙ)·p_{M_k}(y | x, Z₁:ₙ).
  - That is p(M_k | Z₁:ₙ₊₁), so σᵢ is the hierarchical model's posterior predictive, symmetric in the n + 1 units (their lemma). CB-1's rule on it covers at least 1 − α.
  - The per-model scores are `ConformalBayes.Scores`'s, mixed by log-sum-exp. The weights are reported at every candidate.
- *The residual scores* (`ResidualShrinkage`). Full conformal on a fixed predictor's residuals, each of the n + 1 scored against the other n:
  - distance to origin, |r|;
  - distance to average, |r − r̄|, the Bayes score of a flat prior on the mean (their proposition);
  - RoBAS-EB, |r − â·r̄| with â = υ̂²/(υ̂² + σ̂²/n) from their EB estimates. It shrinks to the origin when the residuals sit near zero and to their average under a mean shift.

  The others' mean and spread are read from residuals centred at their overall mean, so they keep their digits under a large shift.
- *Wu et al.'s high-density sets* (`RiskControlledSets`). The split form's density score −ln p̂(y | x) and the set {s ≤ λ*}, λ* the least threshold whose Gaffke bound on the calibration miscoverage is at most α (slice 18's `HighProbability`). The set's miscoverage exceeds α with chance at most 1 − β over calibration sets, a training-conditional guarantee. Their Bayesian quadrature only estimates the set's size, so it is not built.

**Gates.**
- *CBMA with one model* is the full step's set bit for bit, its weight 1.
- *CBMA's weights against the exact evidence.* Two noise laws, σ 1 and 1.5, on 8 units; the exact posterior model probabilities run from 0.015 to 0.669 over the grid. With 4,000 exact draws the weights match them to 0.0021.
- *CBMA across models wrong in their mean* (y = 0.7x + N(0, 1), models βx and βx², n 30, 500 replications, α 0.1):

  | score | coverage | mean width |
  |---|---|---|
  | CBMA | 0.912 | 3.52 |
  | βx alone | 0.912 | 3.44 |
  | βx² alone | 0.910 | 4.14 |

  The mean weight on βx at the truth is 0.958.
- *A model wrong only in its noise scale* (σ 3 against σ 1) barely moves the conformal set: 3.52 against 3.47 at 800 replications. A symmetric unimodal density ranks units by |residual| whatever its scale. So the plan's "wrong noise" bites only where the law's shape or centre differs, and the gate above uses a wrong mean.
- *The residual scores' definitions:* 300 random vectors, some shifted by up to 10⁴, against two-pass definitions.
- *The residual scores under a right and a shifted predictor* (n 15, α 0.1, 2,000 replications, se 0.0067):

  | shift | distance to origin | distance to average | RoBAS-EB |
  |---|---|---|---|
  | 0: coverage, width | 0.930, 4.11 | 0.926, 4.27 | 0.927, 4.17 |
  | 3: coverage, width | 0.946, 9.50 | 0.947, 4.31 | 0.948, 4.32 |

  RoBAS-EB sits between the references when the predictor is right and matches distance to average under the shift.
- *Wu et al.'s sets* (a fit to 50 training units, 60 calibration units, the miscoverage read on 1,000 fresh units, 300 calibration sets, α 0.1, β 0.9): the miscoverage is at most α in 0.957 of the calibration sets (se 0.017), with mean miscoverage 0.046.
- *The refusals:*
  - no models, weights of another count, models over other units, no model with positive weight;
  - a NaN residual, centres and grids of different counts, a grid out of order, an undefined score;
  - a confidence of 1.
- *Per commit,* the selected classes and the pinned tier: 48 of 48 (one class, the new one).

**Defects found: none.**

**Follow-ups, each with the decision it would change.**
- *RoBAS-Full,* the horseshoe score with ₁F₁, beside the EB score. Decides whether the layer offers the fully Bayesian robust score.
- *CBMA's evidence from the engines' own fits.* Here it is exact and given; BASS and BARS would supply estimates such as bridge sampling. Decides whether model averaging reads an engine's evidence directly.
- *A grid-free interval for the residual scores* (their bracketed root search), RoBAS-EB's sets being intervals. Decides whether `ResidualShrinkage` reports exact ends.

**Doc changes owed** (for the offline docs pass): `docs/Maths/Regression/ConformalPrediction.md`:
- a section on scores beyond the model's own: CBMA's score and weights, the three residual scores, and Wu et al.'s risk-controlled sets with their training-conditional guarantee;
- the finding that a wrong noise scale barely moves conformal sets;
- the gates' numbers;
- API rows for `ModelAveraging`, `AveragedPrediction`, `ResidualShrinkage`, `ResidualScore`, `RiskControlledSets` and `HighDensitySets`.

### 2026-11-14 — Slice 18, the guarantee read out

Built at the user's pick (2026-11-14, "Slice 18, the guarantee read out"), under the Lean's sections `LossPosterior` and `ProperScores` (proved 2026-11-10, checked against their sources 2026-11-13). It follows the archived plan's row 18 and its 2026-11-13 correction: b*_β is Gaffke's confidence bound, calibrated over replicated calibration sets, and never named a posterior.

**Landed** on `main`, nothing pushed:
- `5b55a6f`:
  - `LossBound`, with `Cdf`, `Survival`, `Mean` and `Gaffke`;
  - `RiskControl.Mean` and `RiskControl.HighProbability`, giving a `RiskThreshold`;
  - `ConformalDistribution.CellMasses`, `LogScore` and `QuadraticScore`;
  - `GuaranteeReadoutTests`.
- *Nothing existing changes.*
- *Verified by the main session on 2026-11-14,* from a clean archive of `07dc64e` with zero skips: the Commit selection, 1,218 of 1,218 ([receipts](receipts/verify-conformal-s18-5b55a6f/)).

**What it computes.**
- *The law of L⁺ = Σ Uᵢ ℓ₍ᵢ₎,* over n losses in [0, B] and ℓ₍ₙ₊₁₎ = B, U ~ Dir(1, …, 1).
  - P(L⁺ > x) is the divided difference of (t − x)₊ⁿ on the knots, confluent where they tie. It is read in exact rational arithmetic on the losses' binary values, so discrete and tied losses read exactly.
  - Losses of two values take an exact binomial sum, O(n) a reading: there L⁺ = a + (B − a)·Beta(c_B, c_a).
  - `Mean` is (Σ ℓ + B)/(n + 1).
- *`Gaffke(β)`,* b*_β, the least x with P(L⁺ ≤ x) ≥ β, by bisection on the exact law. It is Gaffke's upper confidence bound on the risk J: P(J ≤ b*_β) ≥ β over replicated calibration sets of independent losses at most B (Vlassis and Thomas's Theorem 1). At one loss it is ℓ + β(B − ℓ) (`gaffke_one`).
- *`RiskControl`,* on a loss matrix nonincreasing in the threshold:
  - `Mean`: the first threshold with (Σ ℓ + B)/(n + 1) ≤ α. The expected test loss is at most α (`crc_orderings`), given the losses at the last threshold at most α.
  - `HighProbability`: the first threshold with b*_β ≤ α, by bisection. Its risk exceeds α with chance at most 1 − β (`fail_le_fixed`).
  - Each reports whether any threshold met α.
- *Proper scores on the grid's cells.* `CellMasses` is the predictive's law on (−∞, g₀], (g₀, g₁], …, (g_{G−1}, ∞), from Q̄.
  - `LogScore(y)` is −ln of y's cell (`logLoss_proper`, `logLoss_strict`).
  - `QuadraticScore(y)` is Σ p² − 2p of y's cell, minus Gneiting and Raftery's quadratic score (`gneitingLoss_gap`).
  - Both are NaN for a case with no distribution.

**Gates.**
- *The law is exact:*
  - one loss gives ℓ + β(B − ℓ) on 50 cases;
  - 0/1 losses give Beta(k + 1, n − k), P(Bin(n, x) ≥ k + 1), to 10⁻¹² on 60 cases;
  - three loss sets, ties included, agree with 100,000 Dirichlet draws at five points within 4.5 se.
- *Calibration over replicated calibration sets* (2,000 each), P(J ≤ b*_β):

  | loss law | n 3 | n 8 |
  |---|---|---|
  | Bernoulli(0.3), β 0.75 | 1.000 | 0.941 |
  | Bernoulli(0.3), β 0.9 | 1.000 | 0.938 |
  | three points on {0, ½, 1}, β 0.9 | 1.000 | 0.973 |
  | U², β 0.9 | 1.000 | 0.997 |

  Every rate is at least β, conservative at these small n, discrete losses included.
- *The rules* (miscoverage of {s ≤ λ} on a grid of 101 thresholds, scores uniform, n 30, α 0.1, 1,500 calibration sets): conformal risk control's test loss is 0.093; the high-probability rule (β 0.9) fails 0.026, against the allowed 0.1.
- *Propriety, exactly,* on 500 pairs of laws on 2 to 8 cells: the log score's gap is KL(q ‖ p) > 0, and the quadratic score's gap is Σ (p − q)².
- *On the split predictive distribution* (199 calibration units, 400 test cases, cells of width 1): the right fit beats one ignoring x on the quadratic score, −0.254 against −0.237.
- *The refusals:*
  - losses outside [0, B] or NaN, an infinite B, β outside (0, 1);
  - losses that rise in the threshold, losses above B, no thresholds.
  - A rule no threshold satisfies returns the last, flagged.
- *Existing classes,* the predictive-distribution tests, pass unchanged.
- *Per commit,* the selected classes and the pinned tier: 1,218 of 1,218.

**Defects found: none in code.** One finding about the readout. The split predictive distribution is a step function in the calibration ranks, so a cell no calibration PIT rises across has no mass, and the log score is +∞ for an observation there: 2 of 400 cases for the right fit at cells of width 1. The log-loss's propriety needs the forecast positive where the truth is (`logLoss_proper`'s hypothesis), which such cells break. Over the finite cases the right fit scores 1.48 against 1.61. The quadratic score has no such gap and is the one that ranks fits on these distributions.

**Follow-ups, each with the decision it would change.**
- *A smoothed cell law for the log score:* the PIT ranks' randomized Q spread over the cell, or a floor. Decides whether the layer offers the log score on predictive distributions without +∞.
- *The calibration for n ≥ 2 in the Lean* (Vlassis and Thomas's chain-measure induction), formalized only at one unit here. Decides whether the Lean states what `HighProbability` reports.
- *A set's own miscoverage bound:* b*_β applies to losses of a fixed set rule, such as a validation split's misses. A set whose threshold is chosen on the same calibration units needs the training-conditional law instead. Decides whether conformal sets report a held-out bound.

**Doc changes owed** (for the offline docs pass): `docs/Maths/Regression/ConformalPrediction.md`:
- a section on the guarantee read out: Gaffke's bound as a confidence bound, never a posterior, with Snell and Griffiths's Theorem 4.3 failing as stated;
- the two risk rules and their guarantees;
- the cell scores beside the CRPS, and the log score's +∞ on step-function laws;
- the gates' numbers;
- API rows for `LossBound`, `RiskControl`, `RiskThreshold`, `CellMasses`, `LogScore` and `QuadraticScore`.

### 2026-11-14 — Slice 16, batch predictive inference

Built at the user's pick (2026-11-14, "Slice 16, batches"), under the Lean's section `Batch` (proved 2026-11-10) and Lee, Tchetgen Tchetgen and Dobriban (arXiv 2409.13990v5, read in TeX: §2, Condition 1, Algorithm 1, Theorem 1, the quantile and selection corollaries, the joint-bounds theorem and corollary).

**Landed** on `main`, nothing pushed:
- `60ba2cb`: `BatchConformal`, with
  - `RankLaw` and `Quantile`, giving a `BatchInterval`;
  - `BoxCoverage`, `Box` and `Simultaneous`, giving `BatchBounds`;
  - `Function` and `Select`, the latter giving a `BatchSelection`;
  - and `BatchConformalTests`.
- *Nothing existing changes.*
- *Verified by the main session on 2026-11-14,* from a clean archive of `eb0b064` with zero skips: the Commit selection, 49 of 49 ([receipts](receipts/verify-conformal-s16-60ba2cb/)).

**What it computes.** Batch sets for functions of m test units' scores from n calibration scores, exchangeable together. The scores' floor and ceiling, S_(0) and S_(n+1), are ±∞ unless given.
- *The basis.* The test units' sorted ranks are uniform on H = {1 ≤ r_1 ≤ … ≤ r_m ≤ n + 1}, |H| = C(n + m, m), and every sorted test score lies between S_(r_j − 1) and S_(r_j) at once (`card_testRank_mem`, `batch_sandwich`). Each set reads that sandwich on rank sequences holding at least 1 − α of H.
- *`Quantile`,* the ζ-th smallest score of the batch: [S_(q_L − 1), S_(q_U)], from the exact law of r_(ζ) (their quantile corollary). q_U is the least rank with at most γ of the law above it; q_L the greatest with at most β below. At m = 1 it is split conformal's set.
- *`Box` and `BoxCoverage`,* several order statistics at once: [S_(w_j − 1), S_(q_j)] with the box's share of H (their joint-bounds theorem, `batch_simultaneous`). The count runs by dynamic programming on L(ρ), with the product over consecutive orders, the Lean's correction of their closed form.
- *`Simultaneous`:* their centred search, each order's box centred at round(t_j·n/m) and widened by the least a whose share reaches 1 − α, by bisection.
- *`Function`,* any h monotone in each sorted score: their Algorithm 1, with the rank order h̃ given (h on the ranks by default). It enumerates H twice, refused above 2,000,000 sequences, since computing it is NP-hard in general.
- *`Select`:* scores μ̂·1{y ≤ c}, the threshold T̂ = S_(q_η), and the units with μ̂ above it declared above c, at most η of them falsely with chance at least 1 − α (their selection corollary; `batch_exceed`, `false_claim_iff`).
- *Exact levels.* The rank counts are integers (`BigInteger`). A share of H at most γ is the count at most γ·|H| for the given γ read exactly, so whole levels such as (1 − γ)(n + 1) = 18 land on their rank and do not round away from it.

**Gates.**
- *The counts are enumerations of H:*
  - the law of r_(ζ) at every k, for n ≤ 6 and m ≤ 4;
  - 300 random boxes' shares, to 10⁻¹².
- *At one test unit:* the quantile set is split conformal's, q_U = ⌈(1 − γ)(n + 1)⌉, across five (n, γ), whole levels included.
- *The batch median* (n 499, m 20, β = γ = 0.05, 4,000 batches): it covers 0.900, within its bound 0.9 + max p = 0.907, with mean width 0.94. Bonferroni over single sets covers 1.000 at mean width 6.07.
- *The quartiles at once* (n 200, m 40, α 0.1, 3,000 batches): the box holds 0.911 of H at ranks [21, 79] and [121, 179] of 201, and covers 0.914.
- *The batch mean* (n 12, m 3, scores in [0, 1], 4,000 batches): it covers 0.983 at mean width 0.78 of the range. With unbounded scores the mean's set is the whole line, since a sequence reaching S_(0) = −∞ is always held.
- *`Function` with h the second smallest and h̃ = r_(2)* is the quantile set exactly.
- *Selection* (n 200, m 50, μ̂(x) = x, c 0.7, η 2, α 0.1, 3,000 batches): false claims stay within 2 in 0.916 of the batches, with 4.7 units selected on average.
- *The refusals:*
  - a batch of none, an order outside 1, …, m, levels summing to 1 or more;
  - a NaN score, scores outside the floor or ceiling;
  - orders that do not rise, a box out of order;
  - an H above the enumeration cap;
  - negative predictions, η = m.
- *Per commit,* the selected classes and the pinned tier: 49 of 49 (one class, the new one).

**Defects found: none.**

**Follow-ups, each with the decision it would change.**
- *The rank order h̃ for the mean.* The default, h on the ranks, is conservative here (0.983 at nominal 0.9). Their remark on choosing h̃ suggests orders closer to h's own. Decides whether `Function` offers a better default for sums and means.
- *Sparse functions* (their §5: h reading a few order statistics), as `Quantile` reads one. It decides whether functions of a few order statistics avoid the enumeration.
- *Their upper bound 1 − α + ε* for the quantile, not formalized. The Lean states the lower bound only. Decides whether the layer reports the set's tightness.

**Doc changes owed** (for the offline docs pass): `docs/Maths/Regression/ConformalPrediction.md`, a batch section:
- the sandwich and the four forms, with the Lean names;
- exact levels;
- the gates' numbers;
- API rows for `BatchConformal`, `BatchInterval`, `BatchBounds` and `BatchSelection`.

### 2026-11-14 — Slice 14's follow-up 2, completed: weighted CV+'s merged bound under likelihood ratios

Built at the main session's word once its Lean (`47acb3c`; `cvPlus_merged_likelihoodRatio`, `cvPlus_merged_likelihoodRatio_refined`) was verified. It closes the open item in follow-up 2's report.

**Landed** on `main`, nothing pushed:
- `3f68330`: weighted CV+ under likelihood ratios reports a merged bound (`CvPlusBounds.Merged`), where it reported none. `WeightedCrossFitTests` gains its checks.
- *Verified by the main session on 2026-11-14,* from a clean archive of `1ea5d6b` with zero skips: the Commit selection, 1,213 of 1,213 ([receipts](receipts/verify-conformal-fu2e-3f68330/)).
- *Nothing else changes.* Position weights' 2c stands: its Lean was restated validly at `47acb3c`, with the same conclusion.

**What it computes,** per test point, from the bag's weights (the training units' and the test point's own, W their total):
- *simple:* min(1, 2c(w_max)), with c(γ) = (αW + (K − 1)γ)/(W + (K − 1)γ), K the folds and w_max the bag's largest weight;
- *refined:* γ the largest weight at most αW, L = αW + (K − 1)γ, N = W + (K − 1)γ and e = Σw²/W. The miss is 2L(W + (K − 1)e)/N² where 2L ≤ N, and 1 otherwise;
- where no weight is at most αW, no unit can be refused: miss 0;
- the merged coverage is 1 minus the lesser miss. It applies under `MergedRanks`, as slice 9's merged bound does.
- `Coverage` reports the largest of the three bounds that apply: one fold's, across folds', and merged.

**Gates.**
- *By hand* (three folds of weights 3, 2 and 4, the test's 2, α 0.2): the simple 2·8.2/17, the refined at γ = 2 (L = 6.2, N = 15, e = 21/11), the merged 1 minus their lesser, and the coverage the largest of the three bounds. Every weight above αW gives a merged bound of 1.
- *Covariate shift* (follow-up 2's fixture, 10 folds of 10 units): the merged bound averages 0.553 against one fold's 0.282 over 2,000 replications, while weighted CV+ covers 0.914; the across-folds bound, given w̄ = 2, is about 0.6 there (1 − 2α − 20/W with W near 101, computed rather than read from the run).
- *Existing classes,* the CV+ tests: unchanged and passing.
- *Per commit,* the selected classes and the pinned tier: 1,213 of 1,213.

**Defects found: none.**

**Doc changes owed:** the cross-fitted section's weighted CV+ bounds gain the merged bound under likelihood ratios, its two forms and where it serves (many units), with the Lean's table of which bound serves where.
