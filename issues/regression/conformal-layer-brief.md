# The conformal layer: E16's build (chip brief)

Date: 2026-11-06. Status: open, not started. Written at the user's request from [E16](bass/extensions-brief.md#conformal-bayes-the-consumer) and the [conformal Lean contract](../../lean/Protolemmata/Conformal/conformal-lemmas.md).

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

**Recommended; a chip starts on these unless the user says otherwise.**
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
