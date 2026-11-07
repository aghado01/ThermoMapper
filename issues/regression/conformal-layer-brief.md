# The conformal layer: E16's build (chip brief)

Date: 2026-11-06. Status: slices 0, 1 and 2 built 2026-11-06, the fallback signal read from the inflated set on 2026-11-07 at the user's call (reports below); slice 3 waits for the user. Written at the user's request from [E16](bass/extensions-brief.md#conformal-bayes-the-consumer) and the [conformal Lean contract](../../lean/Protolemmata/Conformal/conformal-lemmas.md).

## Why

Conformal prediction gives a prediction set coverage in finite samples under exchangeability, whatever the model. The model only sets the width. The user's call (2026-11-05) is to build conformal Bayes into GBASS and the regression engine generally, as one layer over retained draws and pointwise densities, serving every engine that has them. The layer owns the guarantee; the score owns the efficiency. So the layer is built once, engine-free, and each engine contributes densities.

What it serves, by problem class:
- prediction sets for a new observation that hold whether or not the noise model is right;
- a calibrated predictive distribution at each x (Vovk's conformal predictive distribution), with PIT and CRPS readouts;
- coverage within strata, for partially exchangeable data;
- a guarantee that does not lean on the sampler, where the chain's mixing is in doubt.

It does not reach a curve's or a parameter's own intervals: those stay with E15's calibration and the resampling entry ([E16](bass/extensions-brief.md#conformal-bayes-the-consumer), "What it does not reach").

## What is already built

- **The pointwise seam (OQ17).** `PointwiseLogDensity` (`src/maths/regression/Pointwise.cs`): the draws × units matrix of log p(y_unit | draw), with the unit partition, each latent's scope and treatment (`LatentReport`), chain lengths, and draw weights for SMC populations. `UnitPartition` carries the unit map and the stratum map, and its docstring already names grouped conformal Bayes as a consumer.
- **Who produces it.** BASS (`BassFit.PointwiseLogDensity`, `BassPointwise`, with each noise law's latent integrated per observation, OQ17's scope rule), BASS's basis density, BARS, BAPS, the hierarchical P-spline, BGP and the changepoint posterior.
- **The self-normalized importance estimate with its Monte Carlo error.** `PredictiveChecks.PsisLoo` computes per-unit self-normalized ratios, each draw's influence, and the relative efficiency from chain lengths. The add-one-in score is the same kind of estimate with other weights.
- **Kish ESS.** `SmcWeights.EffectiveSampleSize`.
- **Prediction at new cases.** `BassPrediction` holds each draw's surface at new cases, and `IBassNoiseLaw.LogDensity(residual, σ², latents)` is the noise law's density with v integrated, the same one the pointwise seam reads.
- **The contract, proved.** `lean/Enthymemata/ConformalValidity.lean`, CB-1 to CB-4, with the lines E16 can take in the [chip's report](../../lean/Protolemmata/Conformal/conformal-chip-brief.md#reports).
- **Quantile fits.** Quantile BASS (E1's slice 2), for CQR's score.

## The diff

### New
- **The candidate density seam.** For the same retained draws, log p(y_c | x_c, draw) at candidate units: a draws × candidates matrix, with the same `LatentReport` semantics and the same draw alignment (count, weights, chain lengths). It is what every engine adds to join the layer. BASS's comes from `BassPrediction` and its noise law, the latents integrated as in the training matrix.
- **The layer** (`Maths.Regression.Conformal`, or as the build finds fit), engine-free:
  - *the full step* (Fong and Holmes): add-one-in weights w_s ∝ draw weight × p(y_c | x_c, θ_s); each unit's conformity score, the self-normalized estimate of p(z_i | all n + 1 units), the candidate's included; the p-value π(y) within the candidate's stratum; the set on a candidate grid in y;
  - *the split form* (CB-2): calibration units scored by a fit on the training split, with scores from the predictive density, CQR's max(q_lo − y, y − q_hi), or a standardized residual;
  - *the conformal predictive distribution:* π(y) over the grid under a score monotone in y, Vovk et al.'s condition, which makes it a calibrated CDF, with PIT and CRPS. The density score is not monotone in y, so it gives sets and not a CDF;
  - *diagnostics:* each candidate's Kish ESS of the add-one-in weights, and CB-3's margin.
- **`docs/Maths/ConformalPrediction.md`,** a page of its own, beside `VarianceMeanMixtures.md`.

### Upgrades
- **`PredictiveChecks`' importance machinery** shared, not copied: the influence-based Monte Carlo error and the chain-aware relative efficiency serve the add-one-in estimate too.
- **A dense product for the score matrix.** The scores for G candidates are a (G × S) by (S × n) product of weights and likelihoods. No dense matrix product turned up in `Maths.LinAlg` on 2026-11-06. Check, and if there is none, build it there as a shared primitive, in lanes, with the efficiency pass's idioms. The user's standing note: spline surfaces are expensive, every bit of performance counts.
- **Each engine's prediction path** gains the candidate density, one adapter per engine.

### Left alone
- The engines' samplers and their pointwise matrices: the layer reads them and writes nothing, so every pinned hash and fingerprint stays.
- `PointwiseContent.CollapsedPredictive`: its terms leave each unit out of a rest that excludes the candidate, which is not the full step's in-sample score, so the layer refuses it.

## Decisions

**Recorded, and executed here.**
- *The layer is engine-wide,* built with GBASS (the user, 2026-11-05).
- *The contract* (CB-1 to CB-4, proved):
  - the accept rule counts ties for the candidate, π(y) > k/m; it covers at least 1 − k/m whatever the ties, and exactly 1 − k/m with distinct or tie-broken scores; the strict form is never used;
  - ranks are taken within the candidate's stratum; a stratum of m units attains only 1 − k/m, so below m = ⌈1/α⌉ the set is the whole space, and the layer says so;
  - every unit, the candidate included, is read by one density. The default is latents integrated; the alternative is the latent-posterior density with the candidate's latent drawn per draw. A matrix whose observation-scope latents are conditioned, or mixed, is refused, as is a candidate whose latents are treated otherwise than the training units';
  - the split form holds for any fit, the sampler's included;
  - an estimate within ε costs at most |M|/m, with |M| read from the estimates; ε must be a sure bound, or its failure probability adds to the miscoverage.

**Recommended, and accepted by the user on 2026-11-06 when the chip for slices 0 and 1 started.**
1. **One seam, not one type per engine.** The layer reads `PointwiseLogDensity` and the candidate matrix. An engine joins by adding the candidate density; nothing in the layer knows BASS.
2. **The full step's content is the retained parameters,** Fong and Holmes's in-sample score. Hoff's leave-one-out score, which carries his optimality, is a later option: it needs the weight g(y_c | θ)/g(y_i | θ) per data unit, add-one-in and leave-one-out at once, with PSIS's instability.
3. **Sets on a grid in y,** its resolution stated with the set. π(y) need not be monotone, so a set can be a union of intervals, and the layer reports it as one. Refinement of the boundaries is an option, not the default.
4. **CB-3's ε from the add-one-in estimate's Monte Carlo error,** through `PredictiveChecks`' influence machinery, reported as indicative: it is not a sure bound, and the Lean chip leaves ε's probabilistic form open. The margin |M| is always reported, inflation by 2ε is an option, and a wide margin is the signal to use the split form. ε's form is the decision the Lean follow-up would change.
5. **Weighted draws work as they are.** An SMC population's or the flat histogram's draw weights multiply the add-one-in weights.

## The slices

One commit per concern; tests with the code they test; each slice ends with its gates and a fresh-eyes review before the next starts.

| Slice | Delivers | First needs | Its own gates |
| --- | --- | --- | --- |
| 0. The layer's core | The full step over given matrices: scores, π(y), sets on a grid, strata, ties, ESS | `PointwiseLogDensity`; the candidate matrix as a plain input | CB-1 by enumeration on a discrete model with exact posterior draws; coverage within binomial error on a conjugate model with exact draws, under noise the model does not share; the small-stratum whole space; CB-4's refusals |
| 1. BASS's candidate seam | The candidate density for every noise law, GBASS's included | `BassPrediction`, `IBassNoiseLaw` | At training cases the candidate density equals `PointwiseLogDensity`'s column; coverage of held-out data on Gaussian, t, quantile and Normal-Wald BASS fits, with each one's width; the cost row |
| 2. The margin and inflation | CB-3's margin, the 2ε set, the fallback signal | Slice 0; `PredictiveChecks`' influence machinery | The margin against the Lean statement's count on perturbed scores; inflated coverage on exact draws |
| 3. The split form | CB-2 for any engine with a candidate seam: density, CQR and residual scores | Slice 1 | Coverage for any fit, a deliberately wrong one and a chain stopped early included; CQR's set an interval |
| 4. The predictive distribution | π(y) as a CDF under a score monotone in y, PIT and CRPS | Slice 0; each draw's CDF at the candidate (`IBassNoiseLaw.At`) | The in-sample predictive CDF is the natural score; whether it stays monotone in y under the add-one-in weights is read first, and the score is adopted only if it does. PIT uniform on held-out data |
| 5. The other engines | BARS's, BAPS's, the hierarchical P-spline's, BGP's and the changepoint posterior's candidate seams | Slice 1's shape | Each one's candidate density equal to its pointwise column at training cases; coverage on one fit each |

## Gates common to every slice

- **Nothing upstream moves.** The layer reads; every pinned hash and fingerprint holds.
- **The contract, in C#.** CB-1 to CB-4 as tests on enumerable or exact-draw fixtures, each test naming its Lean declaration.
- **Coverage from truth outside the model.** Data from a process the model does not share; coverage over replications against 1 − α within its binomial error, and the width beside it.
- **Primitives.** Chosen cases and a seeded random sweep, far ranges included: log-weights far in a tail, scores that tie, strata of one.
- **The cost row.** Per candidate and per grid: the weights, the score product and the ranks, at the cost report's base and at n 1,000 with S 4,000, on the cycle clock, marked indicative.
- **`docs/Maths`.** Decided statistical choices land on the new page with their alternatives and implications.

## Not in this brief

Each is recorded in E16, with the decision it would change:
- the tournament correction (Bhattacharyya, Zhang and Barber), full data at 1 − 2α whatever the draws, at n leave-one-out chains;
- Hoff's leave-one-out score;
- weighted conformal under shift, and conformal Bayes under label shift;
- batch predictive inference and E17's whole-curve bands;
- CBMA across noise laws, which needs each law's evidence, GBASS's slice 5;
- posterior conformal prediction's approximate conditional coverage, and highest-density sets by decision risk.

## Working rules

- On `main`, committing by pathspec with nothing left staged and nothing pushed. Commit messages describe the work, with no attribution lines ([CLAUDE.md](../../CLAUDE.md)).
- Builds set `ArtifactsRoot` inline into `artifacts\chip-conformal\`, deleted at the end.
- Heavy runs at below-normal priority with capped threads; the machine is shared.
- New beside old; a removal is its own question to the user, with what is lost.
- Findings framed by problem class, not by any one application.
- A follow-up goes in a report with the decision it would change and what each outcome would do. A question with no decision waiting on it is recorded as open, not proposed.

## Report

Under "Reports" below, newest first, per slice: what landed by commit; the gates' numbers; coverage and width by fit; the cost row; defects found; decisions taken beyond this brief; what was left unbuilt and why. E16 takes a line per slice that lands.

## Reports

Newest first.

### 2026-11-07 — Slice 2's fallback signal, read from the inflated set

At your call on slice 2's first follow-up, the signal now reads the inflated set against the plain one, not the margin.

**Landed** on `main`, nothing pushed:
- `74c871b` `ConformalSet.Gap`, `AcceptsPlain` and `AcceptsInflated`, the new rule, and its tests.
- `7a7bc3f` the rule, its default and alternatives, and the readings on `ConformalPrediction.md`.

**The rule.**
- *What it reads.* The inflated set holds the plain one, so the values only it accepts are those whose decision rests on the scores' Monte Carlo error.
- *When it fires:* when they exceed a share of the plain set, (|I| − |A|)/|A| > `FallbackShare`, read on the grid's cells, or by count on a grid of one value, where it fires exactly when the true value's decision rests on the error. It also still fires when a candidate in the inflated set has an untrusted k̂.
- *The margin* stays a reported worst case and no longer drives the signal.

**The default, 0.5,** fires when half the plain set's extent again rests on the error. It was set from these readings:

| Case | Gap, median | 90th percentile | Fires at 0.25 | Fires at 0.5 |
|---|---|---|---|---|
| Exact draws, S 8 | 0.157 | 0.313 | 24% | 1% |
| Exact draws, S 160 | 0.063 | 0.115 | 1% | 0% |
| Exact draws, S 2,560 | 0.016 | 0.039 | 0% | 0% |
| BASS Gaussian, S 1,000, two chains | 0.60 | 1.54 | 77% | 55% |
| BASS t, S 1,000, two chains | 0.36 | 0.78 | 69% | 30% |

- The exact-draw grids rarely fire on the gap at 0.5. At S 8 and 16 they fire on trust, since the weights' tail is too short for k̂, and at S 40 on 58 of 100, where PSIS's threshold for that S is 0.38.
- At the true values the signal fires where the decision rests on the error: 53 of 992 Gaussian and 43 of 1,000 t values.
- *At the degenerate Gaussian weights* (ESS 1.5 to 48) both rules refuse, so no decision rests on the error and the signal stays quiet at all eight.

**What it does not show.** The plain set lost no coverage on any fixture: 0.905 at every S on exact draws, and 0.905 and 0.894 on the BASS fits. So the signal reads how far the error could move the set, not a loss measured here. No fixture has yet shown the plain set losing coverage to Monte Carlo error.

**Gates.**
- Two slowly mixing chains (AR(1) with ρ 0.995) fire on the gap alone, at 0.884.
- Six draws fire on trust alone, with a gap of 0.
- 2,000 exact draws do not fire.
- Every set's signal is the rule's at shares 0, 0.01, 0.25 and 1.
- The inflated set holds the plain one throughout.
- 209 of 209 neighbouring tests pass (the eight pins, PSIS-LOO, WAIC, pointwise, candidate and conformal), and 59 of 59 after the review's fix.

**Defect found and fixed in review.** An error that does not exist, which happens when every weighted term of a unit is −∞, came out NaN. The inflated rule would then have skipped a unit the plain rule counts. It is now taken as unbounded.

**Follow-ups, each with the decision it would change.**
- *A fixture where the plain set does lose coverage to Monte Carlo error,* such as chains stopped before mixing. It decides whether the 0.5 default is set against measured harm rather than only against the gap's spread.
- *Whether the margin stays always-on* (2.6 to 3.5 times slice 0's step). Nothing reads it now but the report.

### 2026-11-06 — Slice 2, the margin, inflation and the fallback signal

**Landed** on `main`, nothing pushed:
- `fa2123c` PSIS-LOO's draw influence, log-scale error, tail length and k̂ became shared helpers in `PredictiveChecks`, with its relative efficiency internal. PSIS-LOO's arithmetic is unchanged.
- `e7f0049` CB-3 in the layer (`ConformalMargin`, `MarginBounds`, the new `ConformalOptions` and `ConformalSet` members), with `ConformalMarginTests`, `BassConformalMarginTests` and the cost row.
- `aab438d` the CB-3 section, slice 2's checks and cost, and the API on `ConformalPrediction.md`.

**What it does.**
- *ε:* each score's Monte Carlo standard error, √log(1 + Σψ²/r_eff) with ψ_s = w̃_s(p_s/p̂ − 1). That is `PredictiveChecks`' influence, shared; r_eff is the candidate's density's along the chains. Two more products give it for every unit.
- *The bound:* z standard errors, with z = Φ⁻¹(1 − ν/(2m)) and ν = 0.01, which adds to the miscoverage.
- *The margin, always reported:* the estimated band per grid value, with m̂ the largest over the grid. It is read under the uniform bound by default, or each unit's carried one, and gives the plain set's coverage under error, 1 − (k + m̂)/m − ν.
- *The inflated rule:* per-position bounds, an option, its ranks always reported.
- *Trust:* PSIS's k̂ where the inflated rule accepts.
- *The fallback signal:* m̂ > 0.1·α·m, or an untrusted candidate in the inflated set.

**Gates.** Each contract test names its declarations.
- *The margin against the Lean statement's count on perturbed scores.* The layer's band is the literal `FragileEstimated` count over 3,000 cases with ties and bounds of 0 and ∞, and over 3,000 with per-unit bounds on continuous scores. On 4,000 cases of exact scores, per-unit bounds and estimates within them:
  - `accepts_or_fragile`: 1,771 decisions lost, every one fragile;
  - `acceptsInflated_of_accepts`, `accepts_shift_of_accepts_close`, `fragileEstimated_of_fragile`, `card_fragile_le_estimated` and `fragileEstimated_mono` all hold;
  - the adversary refuses exactly the fragile units.
- *The layer against a reference* (80 cases with weights, chains and far ranges): its errors match Σψ²/r_eff in logs to 1.2·10⁻¹⁰, and its inflated ranks, bounds and bands, uniform and carried, match the definitions on its own scores.
- *Inflated coverage on exact draws* (S 40, n 19, k 2, 3,850 replications, noise outside the model, se 0.0048):

| Set | Coverage | Width |
|---|---|---|
| Inflated | 0.935, against 1 − k/m − ν = 0.89 | 5.72 |
| Plain | 0.905, above the margin's reading of 0.485 | 5.09 |

- *Nothing upstream moves:* 326 of 326 tests pass. They include the eight `BassBitIdentityTests` pins, the pointwise, PSIS-LOO, WAIC, diagnostics, SMC, candidate and conformal tests, on HEAD plus this slice's files. The gated readings ran with `BASS_QUALIFICATION=1`, and each printed its rows.

**The added gate: slice 1's degenerate Gaussian weights.** These are slice 1's held-out fits again (n 99, k 10, S 1,000, 20 fits × 50 values), beside t fits of the same data.
- *What the layer says at the 8 values with ESS below 50* (ESS 1.5 to 48):
  - the candidate's bound is 0.7 to 2.6 nats, against a median 0.08 elsewhere, and π is 0.01 to 0.02;
  - both the plain and the inflated rule refuse every one;
  - the band holds 96 to 98 units under the uniform bound and 34 to 98 carried, and the fallback signal fires.
- *So the degeneracy moves no decision.* The candidate's refusal survives its own bounds. What fires the signal is the error of every unit's score, not the candidate's weights.
- *Everywhere else:* m̂ has median 98 (carried 72) for the Gaussian fits and 94 (carried 52) for the t fits, and the signal fires on 992 of 992 Gaussian values and 969 of 1,000 t values.

| Fit | Plain coverage | Inflated coverage | Plain width | Inflated width |
|---|---|---|---|---|
| Gaussian | 0.905 | 0.958 | 1.20 | 2.07 |
| t | 0.894 | 0.937 | 1.06 | 1.49 |

  k̂ was read at 958 Gaussian and 937 t values in the inflated set, and passed PSIS's threshold at all but one t value.

**The margin as the draws grow** (exact draws, m 20). The fallback signal fired on all 100 data sets at every S, under both bounds.

| S | m̂, uniform | m̂, carried | Coverage under error, carried |
|---|---|---|---|
| 40 | 12.8 | 4.6 | 0.662 |
| 2,560 | 7.7 | 2.5 | 0.763 |

**The cost row,** per grid of 200 at n 1,000:
- the full step is 127.9 ms at S 1,000 and 383.4 ms at S 4,000 with independent draws, and 146.4 and 502.2 ms with two chains;
- slice 0's step was 36.7 and 149.8 ms, so the step now costs 2.6 to 3.5 times as much;
- the largest parts are the two error products, each candidate's r_eff (an FFT per candidate) and k̂ (a sort of S ratios per accepted candidate).

**Defects found.**
- *In the cost test,* not the layer: `Count` overwrites the score product with the scores, so a repeated timed call read scores as products and fell onto the exact log path. That showed seconds per grid. The test now refreshes the product before each timed call.
- *A stopped run's test host* held the Release DLL and failed the next build. The host was stopped and the readings rerun.

**Decisions taken beyond this brief.**
- *ε's bound:* z = Φ⁻¹(1 − ν/(2m)) with ν = 0.01, the union over the stratum's units. Its failure probability then adds to the miscoverage, as CB-3 says.
- *r_eff:* one per candidate, from its density along the chains, rather than the influence's per unit and candidate, which would be an FFT per pair.
- *The relative variance expanded* over two products, summed in logs below 10⁻¹⁰⁰.
- *The carried bounds, built as an option.* Each unit takes the larger of its bound here and as the candidate. A bound sure over every ordering would need n² weightings.
- *k̂ only where the inflated rule accepts;* candidates outside every draw's support left out of m̂.
- *The fallback rule and its share* (0.1 of α's budget).

**Left unbuilt:** bounds sure over every ordering; and a fallback signal that tells a degenerate case from a healthy one, which waits on the decision below.

**Follow-ups, each with the decision it would change.**
- *What the fallback signal reads.* As specified, a wide margin, it fires on every case measured, under both bounds and at every S tried, while the plain set covers within its error and the inflated set above nominal. Three options:
  - Read it from the inflated set instead: fire when the inflated set is wider than the plain one by more than a stated share, or when the two disagree on more than a share of the grid. That measures how far the error moves this case's decisions, and its discrimination would be measured on these fixtures.
  - Take ε at one standard error, without the union over m. The band narrows about 3.5 times, and the stated failure probability is given up.
  - Keep the signal as it is, as "not certifiable at this S", and make the inflated set the default guard.
- *Whether the margin stays always-on.* It costs 2.6 to 3.5 times the step. Opt-in would return the plain step to slice 0's cost, and the inflated rule needs only the per-position errors, not the band.
- *ε's probabilistic form* (open in the Lean follow-up) decides whether the bound is z standard errors or a concentration bound, and with it how wide every reading above is.

### 2026-11-06 — Slice 1, BASS's candidate seam

**Landed** on `main`, nothing pushed:
- `1b51458` fix: a noise law's latent per observation is reported `Integrated` (defect below), and the latent list moves into `BassPointwise.Latents`, which the seam shares.
- `86aefe1` `BassFit.CandidateLogDensity(data, cases, grids, strata, unseenLevels)`, in the new `BassCandidates.cs`, with `BassCandidateTests`.
  - It gives log p(y | x_t, draw) on each new case's grid, in the pointwise density's rows, with its chain lengths or weights and its latent report.
  - The surface is prediction's. The density is the training terms' own, N(y | f, σ²) or the law's with each latent integrated, for every law.
  - Refused: a direct functional fit, a weighted fit, and draws whose working response moved.
- `15c29ce` its section on `ConformalPrediction.md`, and the seam, the corrected report and the API entry on `BayesianRegression.md`.

**Gates.**
- *At training cases the candidate density is the pointwise column, bit for bit:* largest gap 0, over 120 draws × 80 cases for each of Gaussian, t, quantile, Normal-Wald, skewed t and generalized hyperbolic BASS. The latent report and chains are the same, an SMC fit's weights are carried, and the layer takes every pair.
- *Held-out coverage from truth outside every model.* The truth is f = sin(2πx₀) + x₁ with noise 0.2ζ plus, one time in seven, 0.8·Exp(1). There are 20 fits per law of n 99, so k 10 of m 100, nominal 0.900, with 50 held-out cases per fit.
  - The se is the fits' spread: each fit costs seconds, so the design is clustered, and the binomial se alone is 0.0095.
  - The width is read on 20 cases per fit at grid step 0.08. The model's own 90% interval comes from its predictive draws.

| Fit | Coverage (se) | Width | The model's own interval | Its width | Least ESS at the true value |
|---|---|---|---|---|---|
| Gaussian | 0.905 (0.010) | 1.20 | 0.944 | 1.52 | 1.5 (8 of 1,000 below 50) |
| Student t | 0.894 (0.012) | 1.06 | 0.895 | 1.01 | 347 |
| Quantile (q 0.5) | 0.895 (0.012) | 1.02 | 0.921 | 1.22 | 87 |
| Normal-Wald | 0.890 (0.015) | 1.07 | 0.900 | 1.18 | 74 |

  - Two of 1,600 sets split into two intervals, and none reached the grid's ends.
  - The Gaussian fit's own interval overcovers and is wider, its σ² absorbing the contamination. Its add-one-in weights degenerate at outlying candidates, which is CB-3's margin to read (slice 2).
- *The cost row* is per grid of 200 at the cost harness's base point (n 1,000), on the cycle clock, indicative:

| Law | S | The candidate density | Per candidate | The layer's step |
|---|---|---|---|---|
| Gaussian | 1,000 | 3.8 ms | 19 µs | 24.1 ms |
| Gaussian | 4,000 | 5.2 ms | 26 µs | 89.3 ms |
| Quantile | 1,000 | 5.9 ms | 29 µs | 24.0 ms |
| Student t | 1,000 | 12.7 ms | 64 µs | 23.4 ms |
| Generalized hyperbolic | 1,000 | 82.7 ms | 414 µs | 25.3 ms |
| Generalized hyperbolic | 4,000 | 289.6 ms | 1,448 µs | 82.6 ms |
| Normal-Wald | 1,000 | 179.5 ms | 898 µs | 21.2 ms |

- *Counts confirmed:* each of the four coverage facts ran as one test and printed its reading, and the cost row ran 7 of 7, all with `BASS_QUALIFICATION=1`. Without it they return at once, as the gate intends.
- *Nothing upstream moves.* 450 of 450 tests pass: the eight `BassBitIdentityTests` pins, the pointwise, PSIS-LOO, WAIC and diagnostics tests, the t, quantile, Normal-Wald, SMC, prediction and mixing suites, and the conformal and candidate tests. 153 of 153 pass on the exact committed state of the two shared files.

**Defects found.**
- *The noise latent's report* (fixed, `1b51458`). BASS reported each law's latent per observation as `IntegratedGivenOthers`, the collapsed content's treatment, so `LeavesUnitsOut` was true for every noise-law fit. WAIC refused t, quantile, Normal-Wald, skewed t and GH fits, PSIS-LOO's p_loo was NaN, and the layer would have refused them all. The term is in sample, as `VarianceMeanMixtures.md` says. The GBASS session made the same change in its working copy after my note.
- *A cost finding, not a defect.* `BassPrediction.Quantile` under the Normal-Wald takes 6 to 8 s per call at S 1,000 and 50 cases: each draw's generalized hyperbolic builds its tail table. The candidate density there takes 13 to 23 s per fit at 10,000 candidates.
  - The first coverage run, which read the exact mixture quantile and used grids of 201, ran past 100 minutes and was stopped without a reading.
  - The test now reads the model's interval from predictive draws, uses grids of 101 on 20 cases per fit, and runs each law as its own fact.

**Decisions taken beyond this brief.**
- *The seam's refusals.* A candidate is one observation of a scalar response; a weighted fit's terms read a weight a new case lacks; moved working responses condition on offsets a candidate lacks.
- *One latent list* for both densities, so CB-4's match holds by construction.
- *The clustered coverage design and its se,* above.
- *Quadrature laws are left to GBASS.* GBASS's generalized beta prime noise is in progress in the checkout, with `DensityByQuadrature` and a `quadrature` flag on the pointwise path, neither in HEAD. The seam needs the same guard and flag once they land; the GBASS session was told.
- *Commits from the index.* The GBASS session holds uncommitted edits in `BassFit.cs` and `BassPointwise.cs`. Those two files were committed as HEAD plus this slice's hunks (`git hash-object` and `update-index`), and `BassFit.cs`'s hunks were re-added to its working copy, which now differs from HEAD by GBASS's hunks only.
- *Priority.* A reused MSBuild node started at normal priority carried that priority to some of this chip's test hosts: slice 0's runs and the first coverage run, until it was caught and lowered. Node reuse is now off in the runner, so every child inherits below-normal.

**Left unbuilt:** the quadrature laws' guard, which is not in HEAD; candidates for a weighted fit, which would take each case's weight; and a direct functional fit's curves, which are E17's.

**Follow-ups, each with the decision it would change.**
- *A noise law's density over many residuals at one draw.* It decides whether `IBassNoiseLaw` gains a batched density that takes each draw's normalizer once.
  - Under the Normal-Wald and the generalized hyperbolic, the law's density is this slice's cost: 0.4 to 0.9 µs per evaluation, 83 to 180 ms per grid against the layer's 24 ms. The pointwise density would share it.
  - If built, those rows fall by about the normalizer's share. If not, the row stands. It touches `BassNoise.cs`, which the GBASS session is editing now.
- *Open, no decision waiting:* the exact mixture quantile's cost under GH-family laws, for prediction's readers.

### 2026-11-06 — Slice 0, the layer's core

**Landed** on `main`, nothing pushed:
- `4fb54fe` `DenseProduct` in `Maths.LinAlg`, with `DenseProductTests`. No dense product existed there; the search covered `MatrixOps`, `PairwiseSum`, `KhatriRao` and every `Multiply`/`Product` in `src/maths`.
- `cb48ac7` the layer: `CandidateLogDensity` (the seam), `ConformalBayes.FullStep` and `Scores`, `ConformalSet`, with `ConformalBayesTests` and `ConformalCostTests`. `PredictiveChecks`' draw rows and the pointwise density's weight and chain checks became internal and shared.
- `afb92fb` `docs/Maths/ConformalPrediction.md`, a `LinearAlgebra.md` section for the product, and `BayesianRegression.md`'s pointwise section pointing at the layer.

**Gates.** Each contract test names its Lean declarations.
- *CB-1 by enumeration.* The model is discrete: three parameter values under a uniform prior, a unit (x, y) of six values, the posterior held exactly as weighted draws. All 455 data sets of 2 to 5 units were run, each unit in turn the candidate, at every k:
  - the 56 with distinct scores cover exactly 1 − k/m;
  - the 399 with ties cover at least 1 − k/m and at most that plus the straddling share, with 720 (set, k) pairs straddling;
  - the strict rule covers at most 1 − k/m, the tie-broken rule exactly 1 − k/m;
  - `card_rankLe_le` holds, and equal units tie exactly.
- *Per stratum:* the same on every stratum-0 multiset of 2 to 5 beside a fixed stratum 1, ranked within stratum 0.
- *Coverage from truth outside the model.* Exact posterior draws (S 400) of y | x, β ~ N(βx, 1), against a truth sin(1.5x) + (1 + 0.5|x|)(E − 1) that is curved, skewed and heteroscedastic, with n 19 and α 0.1, so k 2.

| Fit | Coverage (nominal 0.900, se 0.0048, 3,850 replications) | Mean width |
|---|---|---|
| The conformal set | 0.902 | 5.25, one interval each |
| The model's own 90% predictive interval | 0.798 | 3.39 |

- *Within strata.* Two strata of 12 units, the noise ×0.2 and ×3 under one homoscedastic model, k 2 of m 13 (nominal 0.846, se 0.0096). They cover 0.854 and 0.824, with widths 2.66 and 9.71.
  - The noisy stratum's 0.824 (z −2.3) was read further. At 20,000 replications on seeds 611, 12345 and 777 it is 0.835, 0.848 and 0.851, and seed 611's next 20,000 give 0.843. Pooled it is 0.844 (z −1.2), so the low reading is seed 611's data.
  - S 3,000 left all 1,400 decisions as S 300 had them, so it is not Monte Carlo error.
  - The fixture keeps seed 611 and its 3.5-se band.
- *The small stratum* is the whole space below ⌈1/α⌉ at α 0.1, 0.15 and 0.5.
- *CB-4's refusals* all fire. `mixed_counterexample` is reproduced: labelled integrated to pass the refusal, the data unit scores 5/8, the candidate 1/2, and coverage is 0 against 1/2. The integrated reading ties them.
- *Primitives.* A seeded sweep of 120 cases against a reference in logs. It covers weights from e⁻⁷⁰⁰ with zeros, chains, log terms near −10⁴ moving 900 nats across draws, candidates in a tail and outside the support, copied columns, and strata of one or none. Every score is within 10⁻¹² relative, with 62 exact ties; counts and π follow exactly, and the inputs are unchanged. The product matches the scalar fused loop bit for bit on chosen shapes and a 60-shape sweep, half of it with subnormals, signed zeros and ±2⁵⁰⁰.
- *Nothing upstream moves.* The eight `BassBitIdentityTests` pins and the pointwise, PSIS-LOO, WAIC and diagnostics tests pass: 264 of 264, and `Maths.LinAlg.Tests` 152 of 152. These ran on a tree of HEAD plus this slice's files.

**The cost row** (`BASS_QUALIFICATION=1`, `BASS_COST_CLOCK=cycles`, Release, best of five, indicative on a shared machine). One grid is G = 200 candidates.

| n, S | Training side, once per stratum | Weights | Product | Ranks | Full step | Per candidate |
|---|---|---|---|---|---|---|
| 1,000, 1,000 | 15.8 ms | 2.5 ms | 14.3 ms | 1.2 ms | 36.7 ms | 12, 72, 6 µs |
| 1,000, 4,000 | 68.4 ms | 12.6 ms | 62.0 ms | 1.7 ms | 149.8 ms | 63, 310, 9 µs |

The product runs at 31 and 27 GFLOP/s, 33 and 40 times the scalar fused loop.

**Defects found.** None in the layer. Two of my own fixture bounds were wrong and were fixed: a count of distinct data sets that cannot exceed 56, and a width ratio above what the fixture gives. One defect upstream, in BASS's noise latent report, belongs to slice 1 and is in its report.

**Decisions taken beyond this brief.**
- *The product sums in order, not pairwise,* by fused multiply-adds. It is bit for bit across hardware, and a candidate whose column equals a unit's gets that unit's score exactly, so CB-1's ties are the model's. Its error on non-negative terms is at most S·u relative.
- *The log floor:* a product entry below 10⁻²⁰⁰ of its unit's largest likelihood is summed again in logs, by the same rule for the candidate.
- *Outside every draw's support,* a candidate takes π = 1/m and ESS 0.
- *A whole-space case computes nothing;* its p-values are null.
- *Random tie-breaking waits for slice 4,* where a uniform p-value is wanted. The step reports `Below` and `Tied` so any rule can be read from it.
- *Kish's ESS only.* The chain-aware relative efficiency and the influence-based error are what slice 2's ε reads, and they are shared there. Slice 0 reads only `PredictiveChecks`' draw rows, which it shares.
- *Isolated builds.* The GBASS session resumed during this slice and edits BASS and distribution files in the checkout. Every build here ran on `git archive HEAD` plus this chip's files, so neither session builds the other's half-made edits.

**Left unbuilt:** CB-4's latent-posterior reading, which needs a second candidate matrix; Hoff's score; refining boundaries between grid values. These are as the brief has them.

**Follow-ups, each with the decision it would change.**
- *A vector exponential* for the training side and the weights. Those are 81 of 150 ms for one case at S 4,000. It decides whether those exponentials may differ in their low bits from `Math.Exp` across hardware. The scores' ties would still hold within one machine, since every column takes the same function.
