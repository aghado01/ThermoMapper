# The conformal layer: E16's build (chip brief)

Date: 2026-11-06. Status: slices 0, 1 and 2 built 2026-11-06, the fallback signal read from the inflated set on 2026-11-07 at the user's call, with a split-half comparator measured beside it; on 2026-11-07 both became diagnostics on every set, with no trigger, at the user's decision; slice 3, the split form, built 2026-11-07 at the user's go-ahead (reports below); slice 4 read on 2026-11-07 (the in-sample predictive CDF is not monotone in y where the draws' predictives separate) and built on 2026-11-08 at the user's decision: the split form's predictive distribution by default, the full step's as an opt-in with its flag; CB-3's per-unit bounds built on 2026-11-08 at the user's order, a slice between 4 and 5, and the margin made opt-in the same day at the user's decision; the next wave approved by the user on 2026-11-08 ("Yes, start it", relayed by the main session): slice 5, the other engines' candidate seams, built the same day (report below) and verified by the main session; E35 as slice 6, its plan below, built the same day (6·0 to 6c). On 2026-11-08 the user asked for E16 to be completed after E35 (relayed by the main session): its plan, slices 7 onward, is below. E35 was verified by the main session on 2026-11-09. Slice 7 landed in e3071cf after those three calls. Slice 8 landed in 331f9e7 and received coordinating sign-off on the combined archive 69398e4. The first proved slice 9 increment landed in 3fbe41b following Lean 4084c55 and received coordinating sign-off on bccb661. The equal-fold CV+ follow-up landed in dd2532c following Lean 33665cf and 17ba5c9; it passed final qualification and awaits coordinating review. Slice 12a, Hoff's leave-one-out score as the full step's opt-in, landed in ab21104 at the user's decision 9; slice10 followed its verification (decision10). Slice10's density and CDF diagnostics landed in4395b18 and6fc1da2 and received coordinating sign-off on2026-11-11: independent67/67 gates,1,053/1,053 Commit and3,161/3,161 whole assembly, zero skips. Slices11 and12b remain held; no further slice is commissioned. Written at the user's request from [E16](bass/extensions-brief.md#conformal-bayes-the-consumer) and the [conformal Lean contract](../../lean/Protolemmata/Conformal/conformal-lemmas.md).

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

**Decided by the user on 2026-11-08: the margin is opt-in** ("Make it opt-in").
- *What the layer does.* The margin, its carried band, its uniform form and the per-stratum orderings only it needs, is off by default (`MarginBounds`, null). Not asked for, it says so: `ConformalMargin.Computed` is false, the band, its count and the carried bounds are null, and the plain set's coverage under error is NaN. The inflated rule stays on every set, per position by default with its uniform option.
- *Where it applies.* It changes the 2026-11-07 decision, diagnostics on every set, for the margin only. The fallback signal and the split-half comparator stay as decided; neither reads the margin's band.
- *Why* (the user's reasons, relayed). On BASS the margin reads 96 to 100 of 100 units, so it fires on every fit and says little. It costs 0.33 to 1.8 s per stratum at n 1,000. The inflated rule gives the protection at no per-stratum cost. It still pays where a unit stands apart, 4.9 against 16.3 in the tail fixture, so it stays available.

**Ordered by the user on 2026-11-08: CB-3's per-unit bounds.**
- *What was asked.* The inflated rule on per-position bounds; the margin on per-unit bounds with the uniform fallback; the uniform forms kept reachable, so the two can be compared; each bound's estimate, and why it dominates. It follows the Lean chip's per-unit CB-3 (80a1993): `coverage_inflated` with bounds by ordering and position, `coverage_close` with bounds carried by the units, `fragileEstimated_mono`.
- *Decided in building it, recorded in the report below.* A unit's carried bound is its largest over every ordering the draws represent, not the two slice 2 read, since two did not dominate. The uniform margin bound is the largest carried one.

**Decided by the user on 2026-11-08: the predictive distribution's score.**
- *What the layer does.* The default is the split form's score, the training fit's predictive CDF F̂(y | x) ("Split-form CDF"; Vovk, Petej, Nouretdinov, Manokhin and Gammerman 2020). The full step's in-sample predictive CDF is an opt-in ("Yes, add it"). It flags each case whose counts fall along its grid, reporting how far, and returns that case's set without a distribution. There is no monotone envelope, which would break the PIT.
- *Why.* Slice 4's reading (2026-11-07) found the full step's score monotone where the draws sit in one mode and falling where their predictives separate. The split form's is monotone by construction and valid for any fit. The full step's trade-off depends on the regime rather than on a flaw: it spends no data on calibration.
- *The measured default (2026-11-08)* keeps the split form. The full step is 2.5 to 10% sharper in CRPS where its draws are a posterior that mixed, more where the data are scarce. It breaks where they are not. The page names that regime as the place to opt in.

**Decided by the user on 2026-11-07: report the fallback diagnostics, with no trigger.**
- *What the layer does.* Every set reports the inflated-set gap and the split-half mass; the split-half is the main session's construction, not from a paper. Neither fires anything, and nothing switches to the split form on its own. The fallback rule (the gap with k̂'s trust) stays as an opt-in with a caller-stated threshold, off by default.
- *Why.* On the BASS fits the plain set covers at nominal, yet the split-half mass at 0.1α fires on 61 to 74% of cases and the gap at 0.5 on 32 to 54%. A default trigger would send most BASS cases to the split form, giving up data for no coverage gain.
- *What would change it.* A trigger can be adopted once a fixture with chains like BASS's shows a measured shortfall to set its threshold by. None does yet, so none is scheduled.
- *Open, no decision waiting:* the main session's reading of why. BASS's chains move between structures slowly, so the halves can disagree while coverage holds. A chain held in one region of the posterior may still score every unit nearly alike, where thin draws do not. If so, the split-half mass reads set instability, which differs from miscoverage. It is a reading, not a result.

**Recommended, and accepted by the user on 2026-11-06 when the chip for slices 0 and 1 started.**
1. **One seam, not one type per engine.** The layer reads `PointwiseLogDensity` and the candidate matrix. An engine joins by adding the candidate density; nothing in the layer knows BASS.
2. **The full step's content is the retained parameters,** Fong and Holmes's in-sample score. Hoff's leave-one-out score, which carries his optimality, is a later option: it needs the weight g(y_c | θ)/g(y_i | θ) per data unit, add-one-in and leave-one-out at once, with PSIS's instability.
3. **Sets on a grid in y,** its resolution stated with the set. π(y) need not be monotone, so a set can be a union of intervals, and the layer reports it as one. Refinement of the boundaries is an option, not the default.
4. **CB-3's ε from the add-one-in estimate's Monte Carlo error,** through `PredictiveChecks`' influence machinery, reported as indicative: it is not a sure bound, and the Lean chip leaves ε's probabilistic form open. The margin |M| is always reported (opt-in since 2026-11-08), inflation by 2ε is an option, and a wide margin is the signal to use the split form. ε's form is the decision the Lean follow-up would change.
5. **Weighted draws work as they are.** An SMC population's or the flat histogram's draw weights multiply the add-one-in weights.

**Decided by the user on 2026-11-09, for E16's completion plan (slices 7 to 18), relayed by the main session:**
6. **Slice 7's flag: the full step refuses the hierarchical P-spline's seam** ("ok"). This removes an input the full step accepted, and the removal is approved. That path never gave valid coverage, because a design shared by every replicate breaks the full step's exchangeability. The split form keeps taking the seam.
7. **Lean first** for the slices that bring new guarantees (9, 10, 12 to 16, 18). Each guarantee is stated, checked numerically and proved in `ConformalValidity.lean`, by the conformal-Lean chip, before its C#, as CB-1 to CB-4 and E35's M-1 and M-2 were.
8. **The order** ("ok"): slices 7 to 12 now, in order. Slices 13 to 18 are held as a queue that the user orders by need.

**Decided by the user on 2026-11-10, relayed by the main session:**
9. **Slice 12 builds Hoff's leave-one-out score as an opt-in** ("Yes, opt-in"). The full step's default score stays Fong and Holmes's in-sample score.
   - *Why opt-in:* the Lean measured a small gain in its finite model, sets 0.2 to 0.5% smaller than the in-sample score's at exact coverage. Each data unit needs its own weight g(z_new | θ)/g(zᵢ | θ), with PSIS's instability.
   - *What slice 12 measures:* the width against the in-sample score on BASS fits. Any question of a default comes back with those numbers.
   - *The Lean's finding:* under CB-1's tie rule on discrete data, the in-sample set can be the smaller one.

10. **The queue, on the conformal chip** (the user, 2026-11-10): slice 12a (Hoff's score, building), then slice 10, slice 11, and slice 12b (CB-4's latent-posterior reading). Each slice is reported and verified before the next starts. Slices 13 to 18 stay held, and the user orders them after slice 10's BASS reading ("After slice 10's reading"). If the margin still reads nothing on BASS fits with the degenerate turns left out, the tournament correction (13) is the candidate to come first. Slices 13 to 16 and 18 take their Lean first (decision 7).

11. **Slice 13 first among 13 to 18** (the user, 2026-11-11, "Yes, 13 first"), on slice 10's BASS reading.
    - *The reading:* the margin reads nothing on BASS fits with the degenerate turns left out. The Kept diagnostic is negative on the Gaussian and Student-t fits at every ESS floor, the carried one is near 0, and the own-turn proxy is not a guaranteed superset.
    - *The ordering:* the tournament correction gives the full-data guarantee that does not lean on mixing, and its Lean is verified (`f998131`). Slices 11 and 12b keep their place ahead of it (decision 10). Slices 14 to 18 follow by need, each with its Lean already verified or landed (14 to 16 verified, 18 landed in `4383916`).

**The Lean for slices 10 and 12 is done and verified** (`dfe5c99`, `a708eb2`, `f82ae97`, `f692e37`; verified on 3b61763). The rules the C# implements are in the [Lean chip's reports](../../lean/Protolemmata/Conformal/conformal-chip-brief.md#reports) and the ledger's [slice 10](../../lean/Protolemmata/Conformal/conformal-lemmas.md#slice-10--the-full-steps-error-completed) and slice 12 sections.
- *Slice 10, the excluded set:* it is read per grid value y over all m units, the candidate included. The data units' flags come from conservative proxies.
- *Slice 10, the reading:* only at kept y. The largest excluded share over the grid is charged, and 2ν in place of ν.
- *Slice 10, the grid* (`coverage_grid_excluding`, `bb5df90`). The maximum over the grid bounds coverage only if the grid holds the test unit's value. On a continuous response, a fixed grid's maximum is a reading at its points, not a bound.
  - Each count is constant between the y where an estimate comparison or a flag changes. A sound report puts a grid value in every such interval, using the breakpoints, as slice 8's refined boundaries do.
  - The PIT is linear between grid values (`pitCdf_interp`), so its largest departure sits at a grid value.
- *Slice 12:* the draw identity cᵢ = A/B̂ᵢ, with the candidate at A. The latent-posterior reading draws the candidate's latent per draw, or uses its closed-form conditional mean.

## The slices

One commit per concern; tests with the code they test; each slice ends with its gates and a fresh-eyes review before the next starts.

| Slice | Delivers | First needs | Its own gates |
| --- | --- | --- | --- |
| 0. The layer's core | The full step over given matrices: scores, π(y), sets on a grid, strata, ties, ESS | `PointwiseLogDensity`; the candidate matrix as a plain input | CB-1 by enumeration on a discrete model with exact posterior draws; coverage within binomial error on a conjugate model with exact draws, under noise the model does not share; the small-stratum whole space; CB-4's refusals |
| 1. BASS's candidate seam | The candidate density for every noise law, GBASS's included | `BassPrediction`, `IBassNoiseLaw` | At training cases the candidate density equals `PointwiseLogDensity`'s column; coverage of held-out data on Gaussian, t, quantile and Normal-Wald BASS fits, with each one's width; the cost row |
| 2. The margin and inflation | CB-3's margin, the 2ε set, the fallback signal | Slice 0; `PredictiveChecks`' influence machinery | The margin against the Lean statement's count on perturbed scores; inflated coverage on exact draws |
| 3. The split form | CB-2 for any engine with a candidate seam: density, CQR and residual scores | Slice 1 | Coverage for any fit, a deliberately wrong one and a chain stopped early included; CQR's set an interval |
| 4. The predictive distribution | π(y) as a CDF under a score monotone in y, PIT and CRPS | Slice 0; each draw's CDF at the candidate (`IBassNoiseLaw.At`) | The in-sample predictive CDF is the natural score; whether it stays monotone in y under the add-one-in weights is read first, and the score is adopted only if it does. PIT uniform on held-out data |
| 4b. CB-3's per-unit bounds | The inflated rule on per-position bounds; the margin on bounds carried by the units, falling back to the uniform bound where no unit stands apart; the uniform forms kept as options. Ordered by the user on 2026-11-08 | Slice 2; the Lean chip's per-unit CB-3 (80a1993) | Coverage at 1 − α on the existing fixtures; no set wider than the uniform rule's; a fixture with a unit in a tail, where the margin gains; slice 2's readings beside the new ones |
| 5. The other engines | BARS's, BAPS's, the hierarchical P-spline's, BGP's and the changepoint posterior's candidate seams | Slice 1's shape | Each one's candidate density equal to its pointwise column at training cases; coverage on one fit each |
| 6. Conformal test martingales | Exchangeability checked online, as [E35](bass/extensions-brief.md#conformal-test-martingales-exchangeability-checked-online-e35) sets out: each new score ranked against every earlier one with the randomized tie, a betting function, and an alarm with its guarantee. Ordered by the user on 2026-11-08 to follow slice 5 | Slice 3's fixed-fit scores; slice 4's log CDF and its rank code | [Slice 6's plan](#slice-6-the-plan-conformal-test-martingales-e35), below |
| 7 to 18. E16's completion | The rest of E16, as its plan below sets out. Asked for by the user on 2026-11-08, after E35 | Slices 0 to 6 | [Its plan](#e16s-completion-the-plan-slices-7-onward), below |

## Slice 6, the plan: conformal test martingales (E35)

Written on 2026-11-08 from [E35's entry](bass/extensions-brief.md#conformal-test-martingales-exchangeability-checked-online-e35) and its Lean (`82d3dc6`; M-1 and M-2 in `ConformalValidity.lean`, with the record in `conformal-lemmas.md`), after the main session verified slice 5.

**What it delivers.** Exchangeability checked as the data arrive. Each new observation's score is ranked among every score so far with a randomized tie. A betting function that reads only the past turns the p-values into capital, and an alarm fires on the capital with its guarantee: a false-alarm rate fixed in advance, with no stopping rule. Like the layer, it is engine-free. It reads scores, and a fixed fit's scores come through the seams.

**The six gates, as the code enforces them.**
1. *The monitored sequence is exchangeable.* That is the caller's premise, stated on the type. In the inductive form the score is fixed by a training set proper that is not in the stream (CB-2's conditioning). The tests break it with a change in the stream and require an alarm.
2. *Each new score is ranked against every score so far, itself included.* One order-statistic structure holds the whole stream's scores. There is no calibration-set entry point.
3. *A fresh θₙ breaks each step's ties,* pₙ = (L + θₙT)/n, with T counting the newest unit itself. θ is drawn from the monitor's own seeded stream, in (0, 1], so p never reads 0, and independently of the scores.
4. *Every score is recomputed at every time.* The inductive form's score is a fixed function, so the monitor takes each score once and keeps no score from an earlier fit. The full form is not in this plan (below).
5. *Each bet is ≥ 0, reads only the past and integrates to at most 1.* The betting function gives fₙ before pₙ is seen, then observes pₙ. Each member's integral is checked by quadrature, and the plug-in's is renormalized by quadrature and scaled down by its rounding.
6. *The capital starts at 1, and alarms read only the past.* ln S₀ = 0, and each alarm reads S₁, …, Sₙ.

**The pieces, by name.**
- `OrderStatisticTree<T>`: an insert-and-rank structure, O(log n) per arrival, returning (Below, Tied) as the split calibration's and the CDF scores' ranks do. It holds doubles or `CdfScore`s.
- `ConformalTestMartingale`: the monitor. It takes each score, ranks it, draws θ and bets, reports the step, and alarms.
- `IBettingFunction` and its members: power, the power mixture, Simple Jumper and the plug-in.
- The scores of a fixed fit, each higher conforming better:
  - the predictive density through the candidate seam;
  - the standardized residual;
  - the PIT's centrality, min(F̂, 1 − F̂), from the log CDF seams.

**The slices.**

| Slice | Delivers | Its own gates |
|---|---|---|
| 6·0. BGP's seam at the run's test points | The main session's follow-up from slice 5. The run keeps each draw's variance of f, unscaled, beside the draw (a shared reference, no copy) and its test inputs on the result. The seam reads them at those inputs instead of refitting | Bit for bit the refit path's; the cost row against slice 5's 20 s |
| 6a. The online p-values (M-1) | The order-statistic tree, the smoothed p-value with its θ, and the fixed fit's three scores | M-1's finite core by enumeration, ties included: each rank sequence 1/n! over the orderings. p-values i.i.d. uniform on exchangeable streams. Each premise broken, with the failure it causes: a fixed calibration set (dependence), CB-1's tie rule (not uniform), scores kept from earlier fits (drift). The tree's ranks against a brute-force count, a seeded sweep with ties and far values. The cost per arrival |
| 6b. Bets and the capital (M-2) | Power, the power mixture, Simple Jumper (Vovk et al. 2021, Algorithm 1) and the plug-in (a kernel density of the earlier p-values); the capital in logs; each step's ln f; the e-value | Each bet's integral by quadrature: 1, or at most 1 for the plug-in. E[Sₙ \| history] = Sₙ₋₁ exactly, on a grid enumeration. A bet that reads the current p, or integrates above 1, breaks it. E S_N = 1 within Monte Carlo error |
| 6c. Alarms and readouts | Ville's alarm, proved. CUSUM and Shiryaev–Roberts as options, their guarantees stated as posed, not proved. The path, the alarm time and procedure, and the guarantee in its terms | Ville's false alarms at most 1/c on exchangeable streams: uniform p-values, and a fixed BASS and BARS fit's scores. CUSUM's and Shiryaev–Roberts's mean time to a false alarm at least c, measured and marked unproven. The delay to an alarm after a change in the noise or the mean, per bet. The cost row |

**Not in this plan, each with the decision it would change.**
- *The full form.* Every score would be recomputed at every arrival under the posterior given the stream so far, the held draws reweighted by the whole stream. That costs n scores and a re-sort per arrival, and the weights degenerate as the stream grows. It decides whether the monitor gains a second score source; until then the inductive form serves.
- *Strata.* A Mondrian martingale is one monitor per stratum. It decides whether the monitor routes by stratum or the caller runs one per stratum.

## E16's completion, the plan: slices 7 onward

Written on 2026-11-08 at the user's request, relayed by the main session: complete E16 once E35 is done. It gathers what [E16's entry](bass/extensions-brief.md#conformal-bayes-the-consumer), this brief's reports and the page record as not built. E17, whole-curve bands, stays its own entry.

**How each slice runs.**
- *Reading first, where the method is not settled.* A slice whose papers have only been read from their abstracts (2026-11-05) reads them first and reports the reading, as slice 4 did, before building.
- *Each new guarantee is named,* with whether the Lean chip should pose it first, as CB-1 to CB-4, M-1 and M-2 were. The decision is the user's (below).
- *As before:* each slice reports, then waits for the main session's verification. Anything that changes a default or removes an option goes to the main session for the user.

| Slice | Delivers | The decision it serves | Its own gates | A new guarantee |
|---|---|---|---|---|
| 7. The seams completed | A flag on `CandidateLogDensity` naming the forms a seam serves, which the full step reads and refuses a split-only seam by; candidates for a weighted BASS fit, each case with its weight (slice 1's left unbuilt) | Whether the full step refuses the hierarchical P-spline's seam in code rather than by docstring. **The user's call:** it removes an input the full step now takes | Each engine's seam names its forms. The full step refuses the hierarchical seam and takes every other. A weighted fit's candidates are its pointwise column at training cases, bit for bit, with coverage on one weighted fit | None |
| 8. Exact ends and quantiles | The split form's set ends and its distribution's quantiles in closed form, F̂⁻¹(· \| x) at the calibration PITs' order statistics, through each engine's quantile function or a root of F̂; a set's boundaries refined between grid values by bisection on π(y) | Whether the split form's distribution and sets need a grid, and whether a set's ends are exact to a tolerance rather than the grid's step | The closed ends against fine grids; the PIT of the closed quantiles uniform on held-out data; refined ends within the tolerance, with π checked either side; the cost against the grid | None |
| 9. Cross-fitted forms | Cross-conformal and the jackknife+ or CV+ (Barber, Candès, Ramdas and Tibshirani 2021); cross-conformal predictive distributions (Vovk et al. 2018); CQR's asymmetric form, a threshold per tail at α/2 | Which form keeps more data than the split form without the full step's reliance on the sampler, at K fits or n | Jackknife+ coverage at least 1 − 2α; CV+'s assumption-specific finite fold bounds from `33665cf` and `17ba5c9`; coverage near 1 − α as measured; width against split at equal data; CQR per-tail coverage | The jackknife+'s 1 − 2α; CV+'s merged-rank, one-fold and fold-tournament bounds; the asymmetric CQR's per tail |
| 10. The full step's error, completed | A CB-3 margin for the CDF score, for the full step's opt-in distribution; the carried bound that leaves out the orderings whose weights degenerate, adding their share to the miscoverage as ν is added (CB-3's open follow-up) | Whether the full step's distribution reports its Monte Carlo error, and whether the margin can read anything on fits like BASS's | Dominance against every ordering the draws represent, as 4b's gates; coverage under error at 1 − α on exact draws; slice 2's and 4b's readings beside the new | The exclusion's added miscoverage, if the Lean chip poses it |
| 11. BGP's drawn f | A retained-parameters content for BGP: f drawn jointly at the training and candidate inputs for each retained state, so the full step takes BGP | Whether BGP gains the full step, at O((n + T)³) per state, beside the split form | The pointwise column at training cases bit for bit; the draws' f against the GP's conditional moments; coverage on one fit; the cost against the split form's | None |
| 12. Hoff's score and the latent posterior | Hoff's leave-one-out score cᵢ = p(yᵢ \| y₋ᵢ), with the add-one-in and leave-one-out weights at once and PSIS's diagnostics (Hoff 2023); CB-4's latent-posterior reading, the candidate's latent drawn per draw, with its second candidate matrix | Whether the full step offers Hoff's Bayes-optimal score, and reads a conditioned observation latent | CB-1 by enumeration with each score; coverage; width against the in-sample score under a right and a wrong prior; the ESS and k̂ where the weights degenerate | CB-4's latent-posterior reading, posed but not proved in the Lean record |
| 13. The tournament correction | Full data at 1 − 2α whatever the draws (Bhattacharyya, Zhang and Barber), from n leave-one-out chains that treat their data symmetrically and need not mix; its 2ε-inflated threshold. Influence-function and RKHS approximations of the full step (2202.01315, 2601.13102) weighed in its reading | Whether the layer offers a guarantee that keeps all the data without leaning on mixing, at n chains | Coverage at least 1 − 2α with chains stopped before mixing; the inflated form's 1 − α − 3√ν on a fixture where ν is known; the cost | 1 − 2α, and the inflated form's |
| 14. Shift | Weighted conformal under covariate shift with a stated or estimated likelihood ratio (Tibshirani, Barber, Candès and Ramdas 2019; weighted Bayesian conformal, Lou and Luo); label shift (Choi); Barber et al.'s weights beyond exchangeability. E35 names it as what comes first for adaptive designs (E30) | Whether the layer takes a weight per unit, and from what | Coverage under a known shift with the true ratio; the loss with an estimated one; exchangeable data unchanged at equal weights | The weighted guarantee, and its loss bound beyond exchangeability |
| 15. Conditional coverage, approximate | Posterior conformal prediction (Zhang and Candès): coverage near nominal per cluster of the score's mixture, beyond named strata | Whether sets adapt to x beyond the stratum map | Marginal coverage exact; per-cluster coverage near nominal on a heteroscedastic fixture; width against Mondrian strata | Its approximate per-cluster coverage |
| 16. Batches | Batch predictive inference (Lee, Tchetgen Tchetgen and Dobriban): sets for functions of a batch of scalar units, simultaneous sets, and selection with false claims controlled | Whether the layer covers a batch jointly | Each of its guarantees on exchangeable batches; width against Bonferroni over single sets | Its batch guarantees |
| 17. Scores beyond the model's own | CBMA across noise laws (Bhagwat, Kong and Jiang), which needs each law's evidence; RoBAS's robust-prior score (Ashouritaklimi et al.); highest-density sets by decision risk (Wu et al.) | Which robust scores the layer offers beside the predictive density | Coverage for each, which any symmetric score keeps; width against one law's score under a wrong noise and a wrong prior; CBMA's weights against each law's evidence | Wu et al.'s risk bound, if it differs from coverage |
| 18. The guarantee read out | Snell and Griffiths's posterior over the loss a calibration set implies, by Bayesian quadrature, for small m; proper scores for the predictive distributions (Gressmann et al.) beside the CRPS | Whether a set reports a posterior on its own coverage where m is small | The posterior's calibration over replicated small calibration sets; the scores' propriety by simulation | Its posterior's calibration |

**Coordination, not a slice here.** The GBASS chip will ask about the candidate seam under BASS's random-intercept block (E21's remaining work). Slice 5's reading applies: a candidate in a seen group conditions on the group's drawn effect, as the hierarchical P-spline's does on its curve. Whether the groups share a design decides which form takes it.

**Not in this plan, each with the decision it would change.**
- *E35's full form and one monitor per stratum,* recorded in slice 6's plan.
- *A noise law's density over many residuals at one draw* (slice 1's follow-up). It belongs to the noise laws' owner, and decides whether `IBassNoiseLaw` gains a batched density.
- *A trigger for the fallback signal,* which waits on a fixture with a measured shortfall (2026-11-07's decision).
- *Dobriban and Lin's joint coverage regions,* which need a conditional pivot a nonparametric fit lacks (E16's entry).

**Decisions for the user, before the slices they touch.**
- *Slice 7's flag* removes an input the full step now takes.
- *Lean first* for slices 9, 10, 12 to 16 and 18, as for CB and M, or C# with the Lean following.
- *The order.* As above, by what each completes and needs. Slices 13 to 18 are wider additions that could be ordered by need instead.

## Slice 12a, the plan: Hoff's score as an opt-in (2026-11-10)

Commissioned by the user on 2026-11-10 ("lets build Hoff's score as an opt-in"), under decision 9, ahead of slices 10 and 11. Slice 12's other half, CB-4's latent-posterior reading, stays queued.

**Own the files.** This slice holds `ConformalBayes.cs` and its tests until it lands. Slice 10's C# touches the same file and waits for it.

**What the Lean fixes** (in `ConformalValidity.lean` and the ledger's slice 12 section; each test names its declaration).
- *The score.* Hoff's leave-one-out score is cᵢ = p(zᵢ | every other unit, the candidate included), in CB-4's integrated reading (b), which the layer already reads for latent-noise models (`hoffJointScore_eq`).
- *From the n-unit draws* (`hoffScore_none`, `hoffScore_some`):
  - the candidate scores A = E_n[g(z_new | θ)], the plain n-unit predictive;
  - each data unit scores A/Bᵢ, with Bᵢ = E_n[g(z_new | θ)/g(zᵢ | θ)];
  - so the candidate is accepted iff #{i : Bᵢ ≥ 1} ≥ k (`hoff_accepts_iff`).
- *Validity.* It is equivariant, so CB-1 covers under the layer's tie rule (ties count for the candidate; `hoff_coverage`). Hoff's strict rule fails under ties and is not offered.
- *Its error* (`coverage_close_shift`, `coverage_inflated_shift`). Â is common to every unit of an ordering, so on the log scale the ranks ignore it. The bounds are log B̂ᵢ's Monte Carlo error for the data units and 0 for the candidate's position.
- *Latents.* Every unit, the candidate included, is read the same way. A mixed reading undercovers (the ledger's counterexample) and stays refused under the existing contract.

**Build.**
- *The option.* `ConformalOptions` gains the score choice. The default stays Fong and Holmes's in-sample score, and every existing set, hash and fingerprint stays unchanged with the option off.
  - Name it by what it does, for example `ConformalScore.InSample` and `ConformalScore.LeaveOneOut`, and credit Hoff (2023) in the docstring.
- *The full step under the option:*
  - the plain set and the inflated rule, the latter on per-position bounds by the shift form;
  - per candidate and data unit, the weights' PSIS k̂ and ESS, as the add-one-in path reports them.
- *What the option does not yet serve, refused with the reason or reported as not computed* (never read with the in-sample score's bounds):
  - the margin (`MarginBounds`), unless its carried bounds for B̂ᵢ come cleanly in this slice;
  - the fallback rule's k̂ trigger, the split-half comparator and the full step's predictive distribution, wherever they read the in-sample score.

  Each is recorded in the report with the decision it would change.
- *`Scores`,* the per-candidate diagnostic, reports the chosen score.

**Gates** (with the common gates below):
- `hoff_accepts_iff` and the draw identity on exact-draw fixtures. On the conjugate fixture the closed-form leave-one-out predictive is the reference.
- Equivariance: permuting the units permutes the scores.
- CB-1 coverage by enumeration with the option on.
- The inflated rule's coverage under the shift form's bounds.
- **Decision 9's measurement:** set width against the in-sample score at equal coverage, and coverage from truth outside the model. Run it on the conjugate fixture and on BASS fits (a healthy Gaussian fit, a Student t fit and a wrong fit, as slice 4's reading used), behind the qualification switch.
  - Optional: the finite model of `conformal-checks.py`, where Hoff's expected size was 0.2 to 0.5% below the in-sample score's.
- The cost row: the extra weight per data unit against the in-sample step.

**Report** under "Reports": what landed by commit, the gates' numbers, the measurement with no default proposed, the refusals and their reasons, and the "Doc changes owed". Message the main session tersely when it lands.

## Gates common to every slice

- **Nothing upstream moves.** The layer reads; every pinned hash and fingerprint holds: the pinned tier at each commit, the whole assembly at the slice's end ([targeted runs](../../CONTRIBUTING.md#testing-standards)).
- **The contract, in C#.** CB-1 to CB-4 as tests on enumerable or exact-draw fixtures, each test naming its Lean declaration.
- **Coverage from truth outside the model.** Data from a process the model does not share; coverage over replications against 1 − α within its binomial error, and the width beside it.
- **Primitives.** Chosen cases and a seeded random sweep, far ranges included: log-weights far in a tail, scores that tie, strata of one.
- **The cost row.** Per candidate and per grid: the weights, the score product and the ranks, at the cost report's base and at n 1,000 with S 4,000, on the cycle clock, marked indicative.
- **Docs, offline.** Chips leave `docs/` alone (the user's standing guidance of 2026-11-10, [AGENTS.md](../../AGENTS.md#docops--conventions-for-issues)). Each report lists under "Doc changes owed" the decided statistical choices, with their alternatives and implications, for the offline docs pass.

## Not in this brief

Each is recorded in E16, with the decision it would change:
- the tournament correction (Bhattacharyya, Zhang and Barber), full data at 1 − 2α whatever the draws, at n leave-one-out chains;
- Hoff's leave-one-out score, now slice 12a (above);
- weighted conformal under shift, and conformal Bayes under label shift;
- batch predictive inference and E17's whole-curve bands;
- CBMA across noise laws, which needs each law's evidence, GBASS's slice 5;
- posterior conformal prediction's approximate conditional coverage, and highest-density sets by decision risk.

## Working rules

- **Messages.** Coordinate shared files and interfaces directly with their owner before editing. For slice 10, the conformal worker owns the full-step code and its tests. E21 owns the grouped BASS candidate seam. Send the coordinator a short event with the report link when a slice lands or needs a decision. Results and receipts go in the report. Sessions do not read each other's transcripts. [Repository DocOps](../../AGENTS.md#docops--conventions-for-issues)
- **Verification.** Build from a snapshot of committed sources plus this worker's own edits. The coordinator independently checks the committed slice. For slice 10, wait for that check on 10a before starting 10b.
- On `main`, committing by pathspec with nothing left staged and nothing pushed. Commit messages describe the work, with no attribution lines ([CLAUDE.md](../../CLAUDE.md)).
- Builds set `ArtifactsRoot` inline into `artifacts\chip-conformal\`, deleted at the end.
- Heavy runs at below-normal priority with capped threads; the machine is shared. Tests run by [targeted runs](../../CONTRIBUTING.md#testing-standards) (`scripts/test-select.ps1`): the selected classes per edit, with the pinned tier per commit, the whole assembly once at the slice's end, and the qualification gates as each slice states.
- New beside old; a removal is its own question to the user, with what is lost.
- Findings framed by problem class, not by any one application.
- A follow-up goes in a report with the decision it would change and what each outcome would do. A question with no decision waiting on it is recorded as open, not proposed.

## Report

Maintain an **In progress** entry under "Reports" with uncommitted work, running checks and the next step. Keep probe code and data in the repository's `artifacts/`. Fold the entry into the slice report when the work lands.

Each report gives the commits, gate results, coverage and width by fit, cost, defects, decisions and open work. Include a **Doc changes owed** heading for the offline docs pass. E16 takes a line per slice that lands.

## Reports

Newest first.

### 2026-11-11 — Slice 12b, CB-4's latent-posterior reading as the full step's opt-in

Built at the main session's word of 2026-11-11, after slice 11's verification, against the Lean's CB-4(c) and slice 12 ("CB-4(c) read from draws").

**Landed** on `main`, nothing pushed:
- `c5a3f1f`, with its tests.
  - *The layer:* `ConformalOptions.ObservationLatents` (`LatentReading.Integrated`, the default; `Conditioned`, the opt-in), and `ConformalBayes.Scores(…, latents)`.
  - *The seam:* `CandidateLogDensity.ConditionedValues`, the candidate's term at its latent drawn from its conditional given its value. Its `Values` stay the integrated terms, which the add-one-in weights read.
  - *BARS:* `observationLatents` (and `seed` for the candidates) on `PointwiseLogDensity` and `CandidateLogDensity`, covering the robust weights and the noise laws. The draws keep W = 1/(c·v).
  - *BASS:* the same on `BassFit.PointwiseLogDensity` and `CandidateLogDensity`, through `BassPointwise` and `BassCandidates`. The draws keep their `ObservationLatents`, and under random effects the residual is net of the draw's group effects. The GBASS chip confirmed both files were free; it gave the recipe and declined to build this side.
  - `LatentPosteriorReadingTests`.
- *Nothing moves with the option off.* The default path's arithmetic is unchanged, and the engines' new parameters default to the integrated reading. The default's refusal of conditioned latents now names the option.

**What it computes.**
- *The training units* read their retained latents: ln N(eᵢ | m·vᵢ, s²·vᵢ), with s² = σ²c and m = σβ, at each draw's own vᵢ.
- *The candidate* draws v′ per draw and per grid value from its exact conditional given its residual, p(v′ | θ_s, y): the t, exponential and GIG mixings, and BARS's two robust-weight resamplers. Its term at that v′ is `ConditionedValues`.
- *The weights* are the integrated g(y | θ_s), as by default. So every unit's weighted reading is CB-4(c)'s latent-posterior score (`latentScore_some_draws`, `latentScore_none_draws`).
- *The candidate's own score and its bound* read the conditioned terms under the same weights. Its error therefore carries the drawn latent's Σ w̃² Var(f | θ, y) (`drawnLatent_mean_var`). The weights' ESS, efficiency and k̂ stay the integrated weights'.
- *A density with no observation latent reads alike either way.* A Gaussian BARS fit gives no conditioned matrix, and its sets are the default's, bit for bit.

**Refused, each with its reason.**
- Under the option:
  - the margin: its carried bounds read the training terms as each unit's turn's weights, which are now f, not g;
  - Hoff's score: the conditioned Hoff reading, `hoffCondScore`, is not built;
  - a training density whose observation latents are integrated: the reading would be the default's under another name;
  - a seam without `ConditionedValues`;
  - latent reports that differ;
  - the CDF options overload.
- In the engines:
  - the CDF under the conditioned reading;
  - a noise without an exact conditional draw (GBP's slice step is a step from the current value, not a draw);
  - BASS's new groups, whose integrated effect makes the latent's conditional given y not the latent's alone;
  - the collapsed content.

**Gates.**
- *The layer reads CB-4(c) exactly.* On a finite model (θ ∈ (0, 1), a latent bit per unit, y ∈ (0, 1, 2)), the exact joint posterior is given as weighted draws, every (θ, v) combination, with the candidate's conditional-mean reading. Every unit's score and the candidate's equal the closed form Σ_θ p(θ | all) Σ_v q f²/g to 10⁻¹² (`latentScore_eq`).
- *CB-1 by enumeration,* every data set of three units and a candidate, k 1 of m 4, bound 0.75.
  - The latent-posterior reading covers 0.785 under an iid law outside the model and 0.799 under the model's own.
  - The integrated reading covers 0.872 and 0.871.
- *The drawn reading against the conditional mean,* 40,000 latents at one draw and residual, against quadrature over ln v from the latent's prior:
  - BARS's t weights: 0.142497 against 0.142573, se 0.000214;
  - BARS's t noise: 0.139793 against 0.139892, se 0.000213;
  - BARS's skewed GIG: 0.257875 against 0.257964, se 0.000310.
- *BARS's seams,* t weights and t noise. Each training term is its draw's Normal at its retained latent, written out. The candidates' values are the integrated terms, and one seed repeats the draws. The full step takes the conditioned pair under the option, and each mix is refused.
- *BASS's seams,* Student t noise. Each training term is N(e | σβv, σ²c·v) at the draw's retained latent and fit, written out. Each candidate's conditioned term is reproduced from the seam's own stream in draw, case and grid order. The values stay the integrated terms. Mixes and the CDF are refused.
- *Per commit,* the selected classes and the pinned tier: 2,659 of 2,659. *At the slice's end,* the whole assembly: 3,214 of 3,214 at `c5a3f1f`.

**The measurement** (behind `BASS_QUALIFICATION=1`; a reading, no default proposed): 20 fits of 99 points (m 100, k 10), 50 held-out cases each at the truth and 20 on a grid. Both readings run on the same fits and cases, with the width difference paired.

| Fit and truth | Coverage, integrated / latent-posterior (se) | Width, integrated / latent-posterior | Latent-posterior less integrated (se) | Least ESS at the truth |
|---|---|---|---|---|
| BARS, t₄ weights; truth t₄ noise, the model's own | 0.886 / 0.885 (0.014) | 0.869 / 0.864 | −0.6% (0.2%) | 249 |
| BARS, t₄ weights; truth slice 5's contaminated, outside the model | 0.895 / 0.895 (0.009) | 0.967 / 0.960 | −0.6% (0.3%) | 125 |
| BASS, Student t noise (slice 4's fit); truth outside every model | 0.895 / 0.896 (0.013) | 1.124 / 1.127 | +0.3% (0.1%) | 161 |

- *Both readings cover within error on every fit.* The width gap is small, and its sign varies with the fit: on BARS the latent-posterior sets are 0.6% narrower (two to three se); on BASS they are 0.3% wider (three se).
- *Hoff's per-unit weights degenerated on this BASS fit* (least ESS 14.4, slice 12a). The latent-posterior reading's weights are the integrated add-one-in weights, least ESS 161.
- *In the finite model,* where y moves strongly with the latent, the Lean measured the latent-posterior sets 1.2 to 4.3% larger. The sign of the gap is the model's.
- *No k̂ above its threshold at the truth* on BARS, under either reading.

**The cost row** (BARS, n 400, S 1,000, 20 cases of 200, cycle clock, indicative, the machine shared):
- the training density: 306 ms integrated, 124 ms conditioned;
- the candidates: 815 ms integrated, 659 ms with the drawn latents;
- the full step: 1,617 ms integrated, 1,375 ms latent-posterior.

The conditioned terms are Normal, cheaper than t's marginal. The drawn latents cost less than that saving.

**Defects found: none.** One gap was closed in review, before any commit: the distribution's options overload would have ignored the reading silently, and it now refuses it.

**Decided in building it.**
- *The option names the reading,* `LatentReading`, and the engines take `LatentTreatment`. This keeps the engines free of the conformal layer's types.
- *The drawn reading only,* the plan's construct. The conditional-mean reading, which removes the drawn latent's variance, is a follow-up.
- *Drawn per grid value,* each candidate value with its own latent, on one stream per call in draw, case and grid order.
- *The split-half diagnostic reads the conditioned terms* as the full step does.

**Follow-ups, each with the decision it would change.**
- *The conditional-mean reading,* ln E[f | θ, y], in closed form for the GIG family (Bessel K ratios), the t and the robust weights, and by quadrature for GBP. It would remove Σ w̃² Var from the candidate's bound and admit GBP. It decides whether the mixing laws gain a conditional-mean member, which is shared by BASS and BARS.
- *The conditioned Hoff reading* (`hoffCondScore`). It decides whether Hoff's score offers the latent-posterior reading. The Lean measured it 1.9 to 3.1% above Hoff's integrated size in the finite model.
- *The margin under the reading.* Carried bounds would need each unit's turn's weights g beside its conditioned terms f, a second training matrix. It decides whether the reading reports a margin.
- *BASS's new groups under the reading:* the joint conditional of the new effect and the latent given y.

**Doc changes owed** (for the offline docs pass; chips leave `docs/` alone).
- *`docs/Maths/Regression/ConformalPrediction.md`, CB-4's section:*
  - the latent-posterior reading as an opt-in: `ConformalOptions.ObservationLatents`, `CandidateLogDensity.ConditionedValues`;
  - what each unit reads, and the weights;
  - the refusals and the drawn latent's variance;
  - the enumeration and the measurement tables.
- *The engine pages,* BARS and BASS: the `observationLatents` and `seed` parameters on the training and candidate seams, the exact-draw laws, and GBP's refusal.
- *The API rows:* `LatentReading`, `ConformalBayes.Scores(…, latents)`.

### 2026-11-11 — Slice 11, BGP's drawn f: the retained-parameters content the full step takes

Built at the main session's word of 2026-11-11 (the resume order 11, 12b, 13), from E16's completion plan.

**Landed** on `main`, nothing pushed:
- `5c25b48`, with its tests.
  - `BgpSampler.DrawFunction(result, xCases, seed)` returns a `BgpDrawnFunction`, whose `PointwiseLogDensity(units, cdf)` and `CandidateLogDensity(grids, strata, cdf, cases)` serve both forms. It also exposes `AtObservation`, `AtCase` and `NoiseVariance`.
  - `GpRegression.JointPosteriorInFrame`, internal; the GP also keeps its responses as given.
  - The collapsed content's docstrings, and its seam's refusal reason, now point to the drawn function.
  - `BgpDrawnFunctionTests`.
- *Nothing moves for existing calls.* The collapsed content, its seam and the run are unchanged; only the refusal's wording changed.

**What it computes.**
- *The draw.* For each retained draw (t, η, s²), g = m + f is drawn at the distinct mapped inputs, the training inputs first and then the cases'.
  - Its law is the kriging mean and s² times the universal-kriging covariance, the constant mean's flat prior included.
  - s² is the run's own exact draw given (t, η), so (t, η, s², g) is a draw from the joint posterior.
  - A draw that repeats its predecessor's state shares the mean and the factor, and takes its own normals. The normals come from a stream per chain, derived from the seed.
- *The terms:* in sample, N(y | g(x), σ²) with σ² = s²η, for the training units and the candidates alike. Every latent is reported conditioned, none left out.
- *A case at a training input* reads that input's g, so its term is the training term, bit for bit. Two units at one input share one value.
- *Pairing.* The training density and the candidates come from one drawn function. A second call draws g afresh, and its candidates would not be jointly drawn with the first call's training terms. The docstring says so; the layer cannot tell.

**Gates.**
- *Bit for bit at the training inputs,* known σ² and unknown scale (n 40, S 300, every draw and unit), the log CDF included. Duplicate cases share g. The full step and its predictive distribution take the pair, and the collapsed content stays refused. One seed repeats the draws, and another draws afresh.
- *The draws against the GP's conditional moments* at one fixed state, 6,000 draws, six points, three of them new, one outside the data's range. The reference is the bordered kriging system [[R, 1], [1ᵀ, 0]] solved by elimination, an algebra the engine does not use.
  - Means fall within 4.5 standard errors.
  - Variances fall within 4.5·√(2/S), with ratios 0.98 to 1.02.
  - Correlations hold to their tolerance, from −0.03 to 0.91.
- *The pivoted factor* reproduces a near-singular covariance (30 inputs, a bandwidth of four spacings, three cases): rank 24 of 33, largest gap 5.4·10⁻¹⁴ against the largest variance 0.128.
- *The refusals:* inputs of another dimension; a case not drawn; a grid count unlike the cases'; a CDF over units of two observations. With no cases, the training density stands alone.
- *Per commit,* the selected classes and the pinned tier: 149 of 149. *At the slice's end,* the whole assembly: 3,175 of 3,175 at `5c25b48`.

**Held-out coverage** (behind `BASS_QUALIFICATION=1`). Slice 5's truth outside the model: f = sin(2πx) + x/2, with noise 0.2ζ plus, one time in seven, 0.8·Exp(1). 20 fits of 99 points, unknown scale, 50 held-out cases per fit at their true value, and 20 on a grid for the width; the se is the fits' spread.

| Form | Coverage (se) | Mean width | Sets in two pieces | Open at the grid's ends |
|---|---|---|---|---|
| The full step on the drawn function (m 100, k 10) | 0.889 (0.012) | 1.017 | 0 of 400 | 0 |
| The split form on the same 99 points, 50 to fit and 49 to calibrate (m 50, k 5) | 0.873 (0.017) | 1.084 | — | 0 |
| Slice 5's split form, 60 to fit and 99 more to calibrate (m 100), for reference | 0.903 (0.014) | 1.12 | 0 of 400 | 0 |

- *Both cover within error.* At equal data the full step's sets are 6% narrower than the split form's, which spends half the points on calibration.
- *The add-one-in weights degenerate at some held-out cases:* the least ESS at the truth is 1.0, at the contamination's outliers, as on BASS fits.

**The cost row** (n 400, S 1,000 with 240 distinct states, a grid of 200 per case, cycle clock, indicative, the machine shared):

| Cases (points) | The drawn function (per state) | Its training density | Its candidates | The full step | The split form's seam | Its step |
|---|---|---|---|---|---|---|
| 1 (401) | 40.4 s (168 ms) | 3.0 ms | 1.4 ms | 74 ms | 13.9 s | 4.2 ms |
| 50 (450) | 57.3 s (239 ms) | 7.9 ms | 63.6 ms | 3.41 s | 16.6 s | 119 ms |

- *The drawn function costs about three times the split form's seam.* Each state pays one O(n³) refit, as the seam does, plus a solve per point and the pivoted factor.
- *It is paid once per call,* for every case it was drawn at.

**Defects found: none in this slice's code.** Two findings during the build, before any commit.
- *The shared `CholeskyDecomposition` diverges on this covariance.* My first draw used it, and g came back NaN.
  - Its floor sets a pivot at or below 10⁻¹² to 10⁻¹². On a covariance whose spectrum falls below rounding, the rows under a floored pivot divide their rounding by 10⁻⁶. The next pivots then went to −3.7·10¹⁵⁷, with 55 entries of the factor infinite.
  - `LogDet` stayed finite (−683.6). Only `Status.Finite` reported it, so a caller that checks the log determinant misses it.
  - The drawn function now uses Cholesky with diagonal pivoting, stopped at 10⁻¹² of the largest variance (Harbrecht, Peters and Schneider 2012).
  - The shared factor is not changed: it belongs to the linear algebra's owner (below).
- *A docstring was wrong after the build:* the collapsed content said BGP's "one content" never draws f. It now points to the drawn function.

**Decided in building it.**
- *g is drawn jointly at every input,* not g(X) first and then g(x*) given g(X). The noise-free conditional of a smooth kernel is the near-singular solve the pivoted factor avoids.
- *The inputs are deduplicated on their mapped values,* so a repeated input is one coordinate, not a singular pair.
- *The training density depends on the cases drawn,* though its law does not, so the pair comes from one object.
- *The terms read the responses as given,* kept by the GP, so a training value read as a candidate reproduces its term exactly.

**Follow-ups, each with the decision it would change.**
- *Drawing g during the run,* while each state's factor exists, as the run keeps its fit variance at its test points (slice 6·0). It would save the drawn function's O(n³) refit per state, nearly all of its cost at small T. It decides whether the run takes the cases in advance for the full step's content.
- *The shared factor's floor.*
  - A caller that reads `LogDet` and not `Status` can sample from an infinite factor. It decides whether `CholeskyDecomposition` refuses a non-finite factor, or offers a pivoted one. If it does, callers that rely on the soft floor change.
  - This belongs to the linear algebra's owner, not this layer.

**Doc changes owed** (for the offline docs pass; chips leave `docs/` alone).
- *`docs/Maths/Regression/ConformalPrediction.md`:*
  - line 35: BGP now serves both forms through its drawn function;
  - the seam table's BGP row (line 444): add `DrawFunction` and its seams, which the full step takes;
  - line 461: two contents, the collapsed and the drawn;
  - the coverage and cost tables: the rows above;
  - the "not built" line (852): built.
- *`docs/Maths/Regression/BGP.md`, "Terms" (line 9):* the retained-parameters content.
  - `DrawFunction`, its law and its pairing rule;
  - the deduplicated inputs and the pivoted factor's tolerance;
  - the bit-for-bit candidate, and the cost.
- *`docs/Maths/LinearAlgebra.md`, the Cholesky section (line 197):* a floored pivot can make the factor's later rows diverge, with `LogDet` finite. `Status.Finite` is the check. The pivoted factor in BGP's drawn function is the alternative for semidefinite covariances.

### 2026-11-11 — Slice 10, density and CDF error diagnostics

Completed after slice 12a's verification, under decision 10's rules and the grid rule (1d76a22). The coordinator verified 10a before releasing 10b, then signed off the completed slice on `6fc1da2`. Defaults are unchanged. The numerical coverage and PIT readouts remain conditional diagnostics, not guaranteed bounds.

**Delivered in two parts, one code/test commit each.**
- *10a, the density score's margin with the degenerate turns left out.* A new `MarginBounds.Kept`.
  - *D(y), one rule for every turn:* a turn degenerates when its weights' ESS falls below `ConformalOptions.DegenerateEss` (100 by default). The ESS is Kish's times the relative efficiency along the chains.
  - *The candidate's own turn* is read directly from its held add-one-in weights.
  - *A data unit's own turn* is read through a proxy, the union of three flags:
    - its ESS without the candidate: Kish(W) · rel · efficiency, with rel = 1/(E[p]·E[1/p]) under the held draws;
    - its PSIS leave-one-out k̂ above the threshold;
    - at kept y, its ESS with the candidate: rel = 1/(Σ w̃ p · Σ w̃/p), from one more product.
  - *The bounds* each unit carries are taken over the turns outside the y-free flags. ĉ(y) = k + the band under those bounds + D's data units that the inflated rule accepts, read at kept y only.
  - *The report:* 1 − (max over kept y of ĉ + max over every y of |D(y)|)/m − 2ν. The inflated set's is 1 − (k + max |D|)/m − ν.
  - *The grid rule:* an opt-in probe (`ConformalOptions.MarginProbe`) bisects every cell whose readings differ at its ends, reusing the stratum's scorer, so a value sits in each interval found. Any margin mode can take it.
  - *Gates:*
    - the proxy against exact own turns, by refits on the conjugate model;
    - dominance and coverage under error on exact draws;
    - the BASS reading beside the carried margin's (decision 10's question), at several ESS floors;
    - whether the margin can serve Hoff's score, recorded.
- *10b, the CDF score's margin* on the full step's distribution (an options overload):
  - the equi-tailed set's union band and inflated tails;
  - the PIT's counts per level, as bounds on P(PIT ≤ k/m);
  - read with D and the grid as in 10a.

**State (November 11, sole conformal lane).** 10a code and tests landed in **`4395b18d2cac7b236a42520bf115542c3e8a48c5`** (engine Git blob `a60871c055854d33234ba2e5556a0cbf128f8b7f`), with nothing left staged and nothing pushed. The coordinator independently verified its exact archive: qualification **15/15**, Commit **2,558/2,559**, sole failure the reproduced pre-E21 assertion ([receipts](../../artifacts/verify-wip-nov11-10a-focused-4395b18/)); **10b released**, 11 and 12b held. The three 10a paths are `src/maths/regression/conformal/ConformalBayes.cs`, `tests/regression/ConformalKeptMarginTests.cs`, and `tests/regression/BassConformalKeptTests.cs`. Other concurrent BASS edits belonged to E21 and were excluded from both 10a worker snapshots. The resumed worker first checked all three live/snapshot SHA-256 hashes and their compiled portable-PDB source checksums: all matched ([original PDB receipt](../../artifacts/worker-conformal-s10/receipts/source-pdb-hashes.json), [original input manifest](../../artifacts/worker-conformal-s10/receipts/input-manifest.json)). The inherited Commit run was preserved and collected; no duplicate heavy selection was started.

- Snapshot: `artifacts/worker-conformal-s10/source`, `git archive 9335d4f` plus only the worker's three C# paths; `base.tar` retains the archive. Latest source includes the analytic exact-own-turn enhancement. Build outputs: `artifacts/worker-conformal-s10/build/`. Settings `two.runsettings` cap xUnit at two threads. Direct builds need `-m:1 -p:BuildInParallel=false`: default parallel restore failed with exit 1 and no diagnostic errors (`restore.log`, `edit.log`); this is a tooling candidate, not a test failure. NU1900 vulnerability-query access warning also occurred.
- Completed receipts in `artifacts/worker-conformal-s10/receipts/`: `focused.trx` retains one incorrect fixture expectation (negative lower readings are permitted); `focused2.trx` passes 7; `qualification-final.trx` passes 12; `qualification-complete.trx` passes 14 including BASS and cost; `analytic-refit.trx` passes the enhanced exact reference measurement. The last analytic enhancement changes tests only relative to qualification-complete. Full logs share these names under the worker artifact root.
- **Commit selection collected:** [Commit.trx](../../artifacts/worker-conformal-s10/build/test-select/20261111_121336-commit/Commit.trx), **1,002 passed / one failed**, 508 seconds, qualification switch off. The sole failure is the known baseline `RandomInterceptTests.Pointwise_ReadsEachDrawsGroupEffects_AndRefusesWhatIsNotBuilt`: expected `NotSupportedException`, got `ArgumentException` for non-predictor column g, before any full-step conformal call. The snapshot test blob `d55567d72973fd36936ab31e769da3cf69c5a199` matches `9335d4f`. This assertion and the candidate seam remain E21's. Outer log `artifacts/worker-conformal-s10/commit.log`.
- **Earlier-blob whole-assembly Slice interrupted:** `scripts/test-select.ps1 -SourceRoot .../source -Changed <the same three paths> -ArtifactsRoot .../build/ -Scope Slice -Threads 1 -NoBuild`, qualification switch off, below-normal priority. At the coordinator's explicit release, stopped only the verified process tree belonging to `build/test-select/20261111_122300-slice`, after 1,077 seconds. Its partial [slice.log](../../artifacts/worker-conformal-s10/slice.log), settings and [interruption receipt](../../artifacts/worker-conformal-s10/receipts/slice-interrupted.json) are retained. It emitted the same known E21 assertion failure and produced no completed TRX; no pass count or whole-assembly success is claimed. The coordinator owns the superseding final combined assembly. This run used the pre-midpoint sources/binaries, untouched by the corrected revision below.
- **Corrected midpoint revision:** coordinator review found that `lo + (hi-lo)/2` skipped a valid finite interval near ±1e308 when the width overflowed. `4395b18` uses the same opposite-sign calculation as `ConformalBisection`, `(lo+hi)*.5`, and tests a probe at zero with differing endpoint flags. The original Commit/Slice and saved qualification receipts belong to the earlier engine SHA-256 `1b102b3b73d9982cb85c34f4dff1469d7e0a7bc731413b45d240c05e3c3beaa8`; they do not test this new branch. New engine SHA-256 **`477553acfaa77cecc425920e13dcf14cdda1037616ee18dbbba4b9501a788ff1`**; only the midpoint expression changes in production relative to that earlier blob. A separate archive-plus-own-edits snapshot `source-midpoint` built serially into `artifacts/chip-conformal/midpoint/` and passed [midpoint-focused.trx](../../artifacts/worker-conformal-s10/receipts/midpoint-focused.trx): **12/12**, qualification switch off (seven deterministic cases including the new regression; five qualification-only cases return without measuring). [Corrected input manifest](../../artifacts/worker-conformal-s10/receipts/midpoint-input-manifest.json) and [PDB receipt](../../artifacts/worker-conformal-s10/receipts/midpoint-pdb-hashes.json) match the committed sources. The coordinator owns the corrected combined Commit and final whole-assembly verification; no duplicate heavy selection is commissioned here.
- Code connects DegenerateEss and the returned flag arrays, computes candidate-conditioned inverse products and kept-only third flags with the **base draw ESS**, excludes y-free flagged turns from carried candidate bounds, reads D's inflated-accepted count, and implements MarginRefinement using the existing scorer. Uniform inflation at flagged y remains unchanged; far-candidate/constant-training tests cover both inflation modes. Refinement probes only cells with differing count/flag readings; equal-end hidden transitions and tails remain unobserved. It is not a continuous-response certificate.
- **Proxy is not sound as a guaranteed flag superset.** Exact analytic Gaussian own-turn ESS, S512, 24 fixtures, floors 10/100/400: at kept y 2261/1957/1254 turns; proxy flags 5/31/181; actual own flags 8/27/173; misses **3/1/3**, extra flags 0/5/11. At flagged y, y-free flags miss **0/1/29** actual own flags, with actual flags 0/4/132 among 19/323/1026 turns. Independent exact posterior refits S16384 disagree with analytic classifications at 5/1/1 kept turns; their own misses are 2/1/3. Largest proxy/analytic ESS ratio 3.8918. Earlier union-reference counts in qualification-final are superseded by `analytic-refit.trx`. `MarginBounds.Kept` and `CoverageUnderError` docstrings explicitly call the numeric result an empirical fixed/refined-grid proxy reading conditional on flag containment and grid coverage, **not a guaranteed lower bound**.
- Exact-draw coverage (1500 replications, n19, alpha .1, seed-tree root 2026111102) at S64/S512: plain .9007/.9060, inflated .9327/.9193 (binomial se .0077). S512 singleton true-value readings at floors10/100/400: plain .8263/.8212/.6652, inflated .8780/.8711/.7594; these are not response-uniform grid bounds. Dominance of computed errors over represented kept orderings is tested separately from actual error containment.
- BASS reading: 5 fits per Gaussian/t/wrong kind, n99, burn-in1000, samples500, 30 truths and four fixed grids per fit, seed-tree root2026111103. Plain/inflated coverages .9067/.9600 (Gaussian), .9200/.9400 (t), .8800/.9000 (wrong). At floors10/100/400, mean kept readings are -.4955/-.4115/-.2495 versus carried -.0945 Gaussian; -.3985/-.6355/-.4150 versus carried .0010 t; .2570/.6270/.3455 versus carried .2025 wrong. Widths plain/inflated 1.205/2.060, 1.100/1.675, 2.480/2.635 respectively. No default or order recommendation. Hoff remains refused because per-unit ratio carried bounds and matching own-turn flags are unbuilt; density bounds cannot be substituted.
- Indicative cycle cost, n1000/G200, best of three: S1000 once/stratum372.76ms, weights/products/counts2.86/73.74/58.44ms per grid (14.31/368.68/292.20us per candidate); S4000 1485.57ms,15.72/259.09/219.66ms (78.61/1295.45/1098.30us).
- **Closed at verified slice10.** Coordinator sign-off on `6fc1da2`: independent gated **67/67**, exact Commit **1,053/1,053**, final combined whole regression assembly **3,161/3,161**, **zero skips**. [Independent gates and three public-API oracles](../../artifacts/verify-wip-nov11-10b-oracle-6fc1da2/qualification.trx) · [Commit](../../artifacts/verify-wip-nov11-10b-6fc1da2/commit-selection.trx) · [Whole assembly](../../artifacts/verify-wip-nov11-10b-6fc1da2/slice.trx). These supersede the interrupted worker runs for final qualification. No further slice is commissioned; 11 and12b remain held.

**Qualification interpretation.** Passing the measurement test is not passing the proxy-containment hypothesis: `analytic-refit.trx` disproves that hypothesis on these fixtures. Here "actual own flags" means the population Kish ESS, scaled to S512, under the exact Gaussian posterior for the data without unit i plus the candidate; independent S16384 refits provide a separate Monte Carlo check. The numerical coverage readouts require containment of those flags, adequate grid coverage, and valid simultaneous error bounds. Only the arithmetic and dominance over represented orderings are tested as properties. Refinement remains a diagnostic search: equal-end transitions, transitions inside the terminal tolerance and the response tails can be missed.

The coordinator separately reported a public-API conjugate count check with the returned flag count below the analytic own-turn count at 3/240 grid points (floor100) and 31/240 (floor400), worst shortfall three units, plus four original reproductions and a grouped/full-step integration check passing ([independent receipt](../../artifacts/verify-e21-conformal-nov11/review.trx)). This corroborates the containment limitation. Final source verification is recorded separately above; sign-off does not promote the proxy/grid diagnostics into guarantees.

The actual BASS comparisons are distinct from those numerical diagnostics. Each model kind has five fits, 150 true responses and 20 fixed grids (161 points from −6 to 10); the fits use n99, burn-in1000 and samples500. Student-t here changes the fitted noise model; `wrong` omits predictor x0. These small fixture measurements establish no model-wide ordering:

| Fitted kind | Actual plain / inflated coverage | Mean grid width plain / inflated | Kept diagnostic, ESS floors 10 / 100 / 400 | Carried diagnostic |
|---|---:|---:|---:|---:|
| Gaussian | .9067 / .9600 | 1.205 / 2.060 | −.4955 / −.4115 / −.2495 | −.0945 |
| Student-t | .9200 / .9400 | 1.100 / 1.675 | −.3985 / −.6355 / −.4150 | .0010 |
| Omitted predictor | .8800 / .9000 | 2.480 / 2.635 | .2570 / .6270 / .3455 | .2025 |

The coverage and widths are unchanged across these floors; the ESS floor changes the diagnostic's flags and charged band, not the underlying plain or inflated rule. Negative diagnostic readings are permitted and remain visible. Defaults stay unchanged: in-sample density score, no requested margin, per-position inflation, ESS floor100 when Kept is selected, and refinement off. Grouped BASS candidates and their tests remain E21's responsibility.

**10b landed and independently verified.** **`6fc1da2bef901188247b15f35410970e2e2a2b85`**, four owned paths: new `ConformalCdfMargin.cs`, edits in `ConformalBayes.cs` and `ConformalDistribution.cs`, and new `ConformalCdfMarginTests.cs`. Snapshot `artifacts/worker-conformal-s10b/source` is archive `f7cee05` (includes E21 `40158dc`) plus these four paths. Corrected qualification passed **39/39**, no skips, `BASS_QUALIFICATION=1`, one test thread, serial build ([TRX](../../artifacts/worker-conformal-s10b/receipts/corrected-qualification.trx)). The independent gates, Commit and assembly are recorded in the sign-off above. [Input manifest](../../artifacts/worker-conformal-s10b/receipts/final-input-manifest.json) and [PDB/snapshot/live checksums](../../artifacts/worker-conformal-s10b/receipts/final-source-pdb-hashes.json) match the committed files. Engine Git blobs: Bayes `83eaf11aa38c3247b4a131591f7cfb49d37346e6`, CDF margin `81a87f8fe9ba0cc11e4bd8a8b422df62f8e2216e`, distribution `6dc93ea7bde08cd292971b121fd0b286215594b7`; test blob `e45ed937bade7174629b73e2dbc6f162cbe2322f`.

Earlier focused2 passed **61/61** with qualification off; first qualification passed **36/36**. The initial focused receipt retains three failed expectations (constant-field summation noise, fixed by centered influences, and the fixture's omitted second tail, corrected to the two-tail formula). A superseded broader qualification including the old `ConformalDistributionTests` gates was stopped only at the coordinator's explicit release: no completed result is claimed; its partial `qualification-final.log`, `broad-inputs/`, PDB hash receipt and `receipts/broad-interrupted.json` remain. No duplicate Commit/Slice selection was started; coordinator owns the integrated selection and final assembly. Index clear after the path-scoped code commit; nothing pushed.

Disposable `artifacts/chip-conformal/midpoint/` and `s10b/` builds were removed after the PDB/source evidence was saved (96,280,146 bytes; [cleanup receipt](../../artifacts/worker-conformal-s10b/receipts/build-cleanup.json)). Worker archives, snapshots, TRX files, partial logs and hash manifests remain in the repository's artifacts. No worker test process remains active.

The options overload retains the original CDF ranks and monotonicity flags. Errors are probability-scale centered self-normalized influences, read on the smaller tail; carried maxima use blocks of at most256 turns, not six N² moment matrices. Weighted training turns normalize `log(weight)+logDensity` before exponentiation (the subnormal `[1,double.Epsilon]`, `[-745,0]` regression agrees with an independent carried-error calculation). Near cancellation/underflow the centered influences are anchored at the largest-weight draw and accumulated by a scaled norm, preserving tiny deviations even when their squares underflow. The chain correction follows the existing likelihood-ratio efficiency convention; no influence-specific MCMC efficiency qualification is claimed. The density-weight proxy is reused only for flags. Counts are read at every rank level: the two tails' fragile sets enter by their union; D's inflated-accepted units and flags are charged separately, and PIT counts use doubled shifts. Refinement reuses the stratum scorer, with the same conditional grid limitations as10a. A stable log-lower-tail normal quantile helper now serves both density and CDF paths, including ν so small that subtracting its tail probability from1 would round to1. Margin-null inflation counts only the candidate's two ranks, avoiding all-unit pair loops.

**10b API and interpretation.** `ConformalBayes.Distribution(training, candidates, options)` attaches `ConformalCdfMargin`. The original overload performs no margin work. `AcceptsPlain`, `AcceptsInflated`, `Inflated` and `CoverageUnderError(alpha)` read the two-tail set; `Lower`, `Upper`, `At`, the monotonicity flag, quantiles and CRPS retain the original ranks. Alpha remains a query on the returned distribution: the margin collects every rank level, rather than freezing one alpha from the options. `Margin.Pit(k)`, `PitAt(u)` and `PitDeparture` read the PIT probability diagnostic, its interpolation and its largest departure from uniformity.

For k per tail, Kept's plain reading is `1 − (2k + max_kept(band_union + D_accepted) + max_all D)/m − 2ν`; the inflated reading is `1 − (2k + max_all D)/m − ν`. Without exclusion the plain reading pays one ν. PIT lower/upper counts are `#{u not flagged: Uhat(u) ≤ k}` and `#{u flagged or Lhat(u) < k}`, with doubled shifts. Over kept response values they take the minimum lower / maximum upper, charge `max D/m + 2ν` at both ends under Kept (ν without exclusion), and clip probabilities to [0,1]. Endpoint probabilities at0 and1 are exact. Interpolation in **PIT rank level** u is the proved linear rule; coverage of the **response grid** y remains conditional. A flagged candidate turn supplies only flags, no band/PIT reading. Neither the CDF field nor the two-tail set cures the proxy-containment failure from10a.

**10b measurements (corrected exact revision, 39/39).** Exact Gaussian posterior draws, n19, 1,000 replications per draw budget, alpha .1, root2026111112: S64 plain/inflated coverage **.9020/.9270**, S512 **.9010/.9130** (nominal-.9 binomial SE .0095). At S512 and floors10/100/400, singleton plain diagnostics .7642/.7641/.6342 and inflated .8783/.8716/.7593; mean PIT departure diagnostics .1894/.2042/.4218. PIT bins held88–111 of1,000 at S512 and90–113 at S64. These singleton response grids do not certify response-uniform coverage, and passing a measurement method is not proving containment.

BASS CDF uses a separate seed tree (root2026111113), five fits per kind, n99, burn-in1000, samples500, 150 true responses and20 fixed grids of161 points from−6 to10. These are not paired density-versus-CDF comparisons with10a's different seed tree. The actual set measurements do not vary across the ESS floors; the numeric diagnostics do:

| Fitted kind | Actual CDF-set coverage plain / inflated | Mean width plain / inflated | Kept diagnostic, floors10 /100 /400 | Carried diagnostic | PIT departure diagnostic, floors10 /100 /400 |
|---|---:|---:|---:|---:|---:|
| Gaussian | .9333 / .9800 | 1.520 / 3.020 | −.3150 / −.3670 / −.2320 | −.0545 | .8885 / .9595 / .9900 |
| Student-t | .9333 / .9733 | 1.320 / 2.155 | −.5425 / −.6755 / −.3455 | −.1025 | .8645 / .9900 / .9900 |
| Omitted predictor | .9067 / .9400 | 2.565 / 2.685 | .3390 / .5715 / .4690 | .2810 | .3680 / .2625 / .3615 |

None of these20 grids per kind was nonmonotone; this fixture observation does not establish general monotonicity. PIT-bin counts and all readouts are retained in `artifacts/worker-conformal-s10b/receipts/corrected-qualification.trx`. The Gaussian/t diagnostics remain largely uninformative despite reasonable measured coverage, while the omitted-predictor fixture reads higher; no default or work-order change follows. Turn blocking bounds temporary moments by N×256, but the opt-in still retains N×S fields/weights and uses quadratic unit-pair counting when a margin is requested; no large-N throughput claim is made.

**Parked shared weighted-row follow-up (coordinator decision; outside the new turn normalizer).** Static review of [PredictiveChecks.Rows](../../src/maths/regression/PredictiveChecks.cs), lines123–145, found that it takes `log(weight/total)` before the conformal paths receive their log weights. Input1, finite raw **`[1e300,1e-300]`**, loses the second log weight to−∞, while subtracting logs would retain approximately−1381.55. Input2, **`[1e308,1e308]`**, overflows the sum despite individually finite weights; [Pointwise.CheckWeights](../../src/maths/regression/Pointwise.cs), lines210–222, accepts it. This predates10b and affects shared predictive consumers; it is not repaired by normalizing the new training turns' combined logs. The coordinator explicitly parked it to keep this change on the qualified `6fc1da2`.

**Follow-up verification criterion.** A separately qualified shared-reader repair must retain the finite second normalized log weight for input1 and return two weights of1/2 for input2, without NaN or loss of rows. Compare normalized log weights and resulting weighted predictive/CDF and carried-error readings against an independent high-precision/log-sum-exp reference. Include a finite log-likelihood tilt that makes input1's initially tiny row material, plus positive-rescaling invariance and zero-weight controls, then run the shared consumers' selection and assembly. That decision would extend extreme raw-weight support across readers beyond E21/slice10; current receipts establish only the tested ranges, including `[1,double.Epsilon]`, rather than arbitrary raw-weight dynamic range.

#### Doc changes owed

- In the offline conformal documentation, show `MarginBounds.Kept`, `DegenerateEss`, `Flagged`, `CandidateFlagged`, `DegenerateAccepted` and `MarginRefinement` together. Label fixed/refined-grid numbers as diagnostics conditional on proxy containment and grid coverage, with the analytic missed-flag evidence above; do not describe them as guaranteed coverage lower bounds.
- Separate actual exact-draw/BASS coverage and width measurements from the diagnostic values. Retain the small-fit scope, the ESS-floor table, unchanged defaults and the refusal of Hoff margins.
- Explain that optional bisection reuses the stratum scorer and can only inspect cells whose endpoint readings differ. It certifies neither interval coverage nor response tails. Public XML comments on `MarginBounds.Kept` and `CoverageUnderError` already carry the conditional interpretation; `docs/` stays with the offline documentation pass.
- Add the CDF options overload and `ConformalCdfMargin`: probability-scale errors, the two-tail union, inflated-tail set queries and PIT counts for every rank level. Distinguish exact interpolation over PIT rank levels from conditional coverage of the response grid. Inflation changes set queries, while distribution ranks, monotonicity flags, quantiles and CRPS retain the original form. Show the separate CDF BASS table and its different seed tree. Note blockwise carried maxima, stable combined-log normalization, the log-tail normal quantile and centered-influence fallback without promoting these numerical repairs into coverage guarantees.

### 2026-11-10 — Slice 12a, Hoff's leave-one-out score as the full step's opt-in

Built at the user's decision 9 and commission of 2026-11-10, relayed by the main session, against the Lean's slice 12.

**Landed** on `main`, nothing pushed:
- `ab21104`, with its tests:
  - `FullStepScore` (`InSample`, the default; `LeaveOneOut`) on `ConformalOptions.Score`;
  - `ConformalOptions.LeaveOneOutParetoK`, off by default;
  - `LeaveOneOutScorer` (internal), beside the in-sample `StratumScorer`, which is untouched;
  - `ConformalBayes.Scores(…, score)`;
  - `HoffScoreTests`.
  - `ConformalBayesTests`' discrete-model helpers became internal, for reuse.
- *Nothing moves with the option off.* The in-sample path's code is unchanged, so every hash and fingerprint holds. The one shared change is that `ConformalMargin.Untrusted` counts a k̂ not read (NaN) as none; the in-sample path never leaves a NaN where the inflated rule accepts.

**What it computes** (each test names its declaration).
- *The candidate* scores ln A = ln Σ_s W̄_s p(y | θ_s) (`hoffScore_none`).
- *A data unit* scores ln A − ln B̂ᵢ, with B̂ᵢ = Σ_s W̄_s p(y | θ_s)/p(zᵢ | θ_s) (`hoffScore_some`).
  - As ln cᵢ = −ln Σ_s w̃_s/p(zᵢ | θ_s), it is one (n × S)·(S × G) product of the units' inverse likelihoods, each scaled by its smallest log term, against the add-one-in weights.
  - An entry below 10⁻¹⁰⁰ is summed again in logs.
- *Ties.* A unit within 10⁻⁹ of the candidate is read again as Bᵢ = Σ_s W_s e^(ℓ_new − ℓᵢ)/Σ_s W_s in draw order, so a unit whose column equals the candidate's ties exactly (`hoff_accepts_iff`, `hoff_coverage`).
- *The error, by the shift form* (`coverage_close_shift`, `coverage_inflated_shift`). A data unit's bound is ln B̂ᵢ's Monte Carlo error, its relative variance G₂/G₁² − 2G₃/G₁ + Σ W̄² from two more products. The candidate's position bound is 0, since Â is common to every unit.
- *Diagnostics per grid value:*
  - the least Kish ESS of the data units' weights, G₁²/G₂;
  - asked for, the largest k̂ among the deciding units, those within their own bound of the candidate, where the inflated rule accepts; −∞ where no unit decides.

**Refused, or not computed, with the reason.**
- *The margin* (`MarginBounds`) is refused: its carried bounds are built for the in-sample score's weights, and each data unit here has its own.
- *The fallback rule* (`FallbackShare`) is refused: it reads the margin's trust in the add-one-in weights' k̂.
- *The split-half reading* is not computed (NaN): it reads the in-sample score's halves.
- *The full step's predictive distribution* is unaffected: it takes no options, and its score is the in-sample predictive CDF.
- *A unit with zero density under a retained draw* is refused: its weight there is unbounded, and the n-unit posterior gives that draw no mass.

**Gates.**
- *The discrete model with exact posterior weights* (504 data set and candidate pairs of four units). Each score is Hoff's ratio p(z)/p(z₋ᵢ) to 10⁻¹², and the set accepts exactly when at least k units have Bᵢ ≥ 1.
- *CB-1 by enumeration, the option on:* 53 distinct data sets, 398 with ties, 4 skipped as within 10⁻⁹. Coverage is at least 1 − k/m on every set, and exactly that with distinct scores; equal units tie exactly.
- *The conjugate model, 40,000 exact draws.* The draw identity holds to 10⁻⁹ against a direct B̂ᵢ. Every unit meets the closed-form leave-one-out predictive within 1.5 of its own reported error, and the candidate meets the closed-form predictive within 0.01.
- *Equivariance:* permuting the units permutes the scores bit for bit (`hoffScore_equivariant`).
- *A seeded sweep* (60 trials). Log terms reach −700, draw weights reach far into a tail, and strata hold one unit. Every score matches a log-sum-exp reference to 10⁻⁹, and 60 of 60 exact ties hold.
- *The inflated rule under the shift form,* S 40, 4,000 replications: plain coverage 0.895, inflated 0.917 per position and 0.947 uniform, against 1 − k/m − ν = 0.89.
- *The option's refusals,* and k̂ NaN unless asked for.
- *Per commit,* the selected classes and the pinned tier: 959 of 959. *At the slice's end,* the whole assembly: 3,073 of 3,073 at `ab21104`.

**Decision 9's measurement** (behind `BASS_QUALIFICATION=1`; a reading, no default proposed).

*The conjugate model,* n 19, m 20, k 2. Coverage is over 4,000 replications and widths over 1,000, with the paired difference's standard error:

| Truth | S | Coverage, in-sample / Hoff | Width, in-sample / Hoff | Hoff − in-sample (se) |
|---|---|---|---|---|
| The model's own prior predictive | 400 | 0.899 / 0.896 | 3.593 / 3.588 | −0.0045 (0.0072), −0.12% |
| The model's own prior predictive | 4,000 | 0.900 / 0.897 | 3.594 / 3.590 | −0.0035 (0.0073), −0.10% |
| Outside the model | 400 | 0.909 / 0.910 | 5.393 / 5.461 | +0.068 (0.012), +1.27% |
| Outside the model | 4,000 | 0.910 / 0.910 | 5.398 / 5.465 | +0.067 (0.012), +1.24% |

- *Under the model's own law,* `hoff_expectedSize_le`'s setting, the two scores are indistinguishable at n 19. The Lean's finite model found 0.2 to 0.5% at m 3 and 4, and expected the gap to shrink with n.
- *Under a truth outside the model,* Hoff's sets are 1.2% wider, about five standard errors.
- *The widths agree at S 400 and S 4,000,* so the gap is the score's, not Monte Carlo error.

*BASS fits,* slice 4's: n 99, m 100, k 10, 20 fits, 50 held-out cases each, widths on 20 cases a fit. k̂ was read, over the deciding units:

| Fit | Coverage, in-sample / Hoff (se) | Width, in-sample / Hoff | Hoff's weights: least ESS, below 50 | Full step's time, in-sample / Hoff, over 20 fits |
|---|---|---|---|---|
| Healthy Gaussian | 0.889 / 0.896 (0.012) | 1.291 / 1.302, +0.9% | 1.0; 994 of 1,000 | 27 / 202 s |
| Student t | 0.895 / 0.900 (0.010) | 1.124 / 1.148, +2.1% | 14.4; 6 of 1,000 | 25 / 135 s |
| Wrong fit (blind to x₀) | 0.885 / 0.887 (0.009) | 2.520 / 2.521, +0.05% | 4.0; 317 of 1,000 | 32 / 55 s |

- *Both scores cover within error on every fit.* Hoff's sets are 0.05 to 2.1% wider.
- *The healthy Gaussian fit's per-unit weights degenerate* at nearly every true value. Its add-one-in weights already did, in slice 1's reading.
- *k̂ above its threshold* came at 322, 213 and 104 of the inflated sets holding the truth, for the healthy, Student t and wrong fit.

**The cost row** (n 1,000, G 200, cycle clock, indicative, the machine loaded by other sessions):

| S | The full step, in-sample | Hoff's, k̂ off | Hoff's, k̂ on (one call) | Hoff's training side | Hoff's readings per grid |
|---|---|---|---|---|---|
| 1,000 | 412 ms | 176 ms | 50.6 s | 49 ms | 138 ms |
| 4,000 | 636 ms | 371 ms | 35.3 s | 98 ms | 244 ms |

- *Hoff's score itself costs less than the in-sample step,* which also reads the split-half halves.
- *k̂ over the deciding units is the payload.* In the cost row's synthetic likelihoods, which are no posterior, 47,000 to 109,000 units decide across the accepted values.

**Defects found: none in the layer.** Two of my readings changed during the slice, before any commit.
- *The first k̂ read every data unit,* fifteen times the in-sample step on the BASS fits and hours at n 1,000. It now reads the deciding units only.
- *k̂ became an opt-in,* by the user's standing guidance that a payload making a path expensive becomes opt-in on that path, keeping its output. **This departs from the plan's wording,** "the weights' PSIS k̂ … as the add-one-in path reports them". It is for the main session and the user: k̂ can be default-on again by one line, at the cost above.

**Decided in building it.**
- *A scorer of its own beside the in-sample one,* so the default path is untouched bit for bit.
- *The reported ESS* is the least over the data units' weights.
- *The error's relative efficiency* is the candidate's density's along the chains, as the in-sample path's is; each unit's own ratio's efficiency would cost S per unit and candidate.
- *Zero density under a draw is refused,* not read as an infinite weight.

**Follow-ups, each with the decision it would change.**
- *The margin for Hoff's score:* carried bounds of ln B̂ᵢ over the orderings. Slice 10 records whether the margin can now serve it. If it can, the refusal lifts; if not, it stays.
- *Each unit's own relative efficiency* in its error. It would decide whether the bounds read each ratio's chain efficiency, at S per unit and candidate.
- *k̂ by default under Hoff's score:* the one-line change above, if the user wants it on despite the cost.
- *Open, no decision waiting:* the gap's dependence on n under the model's own law.

**Doc changes owed** (for the offline docs pass; chips leave `docs/` alone).
- *`docs/Maths/ConformalPrediction.md`, a subsection on the full step's score choice:*
  - Fong and Holmes's in-sample score, the default, against Hoff's leave-one-out score, the opt-in;
  - the rule A/B̂ᵢ and `hoff_accepts_iff`; the shift form's bounds; ties exact;
  - the diagnostics and the k̂ opt-in with its cost;
  - the refusals;
  - decision 9's measurement, both tables, with the implication: alike under the model's own law at n 19, up to 2% wider outside it, and dearer only through k̂.
- *Its API rows:* `ConformalOptions.Score`, `FullStepScore`, `ConformalOptions.LeaveOneOutParetoK`, `ConformalBayes.Scores(…, score)`, and the `ConformalSet` fields' meaning under the option (`CandidateScores` is ln A, `EffectiveSampleSizes` the least unit's).

### 2026-11-10 — Weighted BASS candidates: `Bass.Fit`'s case-weights entry, then the seam's weighted reading

This closes slice 7's open item, "The weighted BASS fit's candidates: not built, blocked on a missing entry point". Both parts are on `main`: the entry point in `7343bc8`, and the candidates in `c1bbf51`. Each was committed through a private index, so the other sessions' uncommitted hunks in `BassModel.cs` stayed out (agreed with the GBASS chip).

**Part 1, the entry point (`7343bc8`).** `BassData` takes optional `CaseWeights`: one positive, finite weight per case, a `WeightMeaning` stated with them, and the coefficient prior's `WeightedPriorGram`. `BassModel.Create` hands them to the working response. `BassResolved`, the design, the target and the draws then read them as the weighted core does, so `Bass.Fit`, `Resume` and every view on a fit run the weighted core.
- *Where they sit.* The weights are on the data, per case and copied, like a predictor. A weight is a statement about its row, and `BassData` is what the fingerprint and `BassFit.ModelOn` already read. The meaning has no default: precision and replication differ in σ²'s ν and in the marginal, and only the caller knows which one holds. The prior's Gram defaults to `Weighted`. Case weights are data fixed through the run, so that prior is the GLS g-prior at the unit block's cost. Under replication, only `Weighted` makes integer weights duplicated rows through the marginal (BASS's page). `Unweighted` keeps D03's prior. `Frozen` is refused, since it freezes a W that moves and case weights do not.
- *The fingerprint* takes the weights, their meaning and their prior, last, and only when declared, so unweighted data keep their fingerprint and a resume cannot change the weights.
- *Refused, each with its reason:* case weights beside a declared noise or another observation model (a noise's latents hand the core weights of their own, and their product is not built); random effects (the block's whitening reads precision weights alone, and a drawn scale's frozen reference would have to read them); a direct functional response (a curve's weight stands at each of its locations, which the direct design holds stacked, and a replicated curve's count would have to reach the support); several outputs. Also a weight that is not positive and finite, a count that is not one per case, and an empty set.
- *Gates* (`BassCaseWeightsTests`, 10 tests). Integer replication weights fit as the duplicated rows through `Bass.Fit`: two chains' retained draws agree state by state, the terms exactly and τ, σ² and the coefficients to 10⁻⁸, 10⁻⁸ and 10⁻⁷ (60 scattered cases, 50 draws a chain, the structure moving 33 times). Each weighted row's pointwise term is its copies' sum to 7.9·10⁻¹³. The tempered marginal and its slope through the data agree with the duplicates at three b, three σ² priors, both forms and three intercepts. Precision weights on the data are the core's weighted target bit for bit, under either Gram: the resolved record, S, the determinant term, the marginal and ℓ_b. A fit's pointwise terms are N(y | f, σ²/w), written out, to 10⁻¹². A resume continues on its own weights exactly as the uninterrupted run does, and refuses other weights or none.
- *Unweighted fits are unchanged bit for bit.* The commit run (selection and pinned tier) ran on a snapshot of HEAD with only this work's files: 2,208 of 2,212 passed. The four failures were all part 2's, fixed before its commit.

**Part 2, the candidates (`c1bbf51`).** `BassFit.CandidateLogDensity(..., caseWeights)` takes each new case's weight. On a weighted fit, each candidate reads the training row's own term in the same arithmetic: `normalizer + ½ln w − ½w·e²/σ²` under precision, `w·(normalizer − ½e²/σ²)` under replication. The two terms and the log CDF now live in `BassPointwise` (`NormalTerm`, `WeightedNormalTerm`, `WeightedNormalLogCdf`), which the training terms, the candidates and `UnitLogCdf` all call, so a training case's candidate is its column bit for bit. The log CDF is Φ(√w·e/σ) under precision, and Φ(e/σ) under replication, since each copy is one observation of N(f, σ²). `UnitLogCdf` no longer refuses weighted fits. Working responses that move per draw stay refused.
- *Refused,* each an `ArgumentException`: no weights for a weighted fit (missing); a weight that is not positive and finite; a count that is not one per case; weights given to a fit whose terms read none, Gaussian or under a noise (extra).
- *A defect the gates found, fixed where the reading is read.* A draw keeps its working response whenever that response is not the response itself. So a fixed-weight draw holds the fit-time model's instance, and `BassFit.ModelOn` builds a fresh model for each view. The seam's "moved" check compared instances, so it read every weighted draw as moved and refused. The pointwise latent report did the same, and the full step refuses that report. `BassPointwise.Moved` now reads a draw as moved only when its reading differs from the model's in offset, weights, meaning, prior Gram, latents or whitening. Each distinct instance is compared once. The latent report and both candidate checks use it, and readings that differ by draw stay conditioned on and refused.
  - *The route not taken:* storing no working response on a draw composed on the model's own reading (`BassSweep`, agreed with the GBASS chip, then withdrawn). It broke nine `BassMovingResponseTests`, which pin each draw keeping its reading: the draw's generative state, which a test observation model cannot rebuild from the data.
- *Gates* (`BassWeightedCandidateTests`). At the training cases, each with its own weight, the candidates are the pointwise columns bit for bit (`Assert.Equal` on every draw and case), log CDF included, under both meanings. The same holds on a tempered ladder's b = 1 replica and on SMC's replicate populations. The full step and its predictive distribution take the pair. The refusals for a missing, non-positive, extra or miscounted weight hold. The commit run on a snapshot of HEAD with only part 2's files: 2,214 of 2,214. *At the end,* the whole assembly on a snapshot of `c1bbf51`: 3,055 of 3,055.
- *Two existing tests changed* (cleared with the conformal chip): `BassCandidateTests.Refuses_...` now expects the missing-weight `ArgumentException`, and `BassConformalDistributionTests.LogCdf_Refuses_...` drops its four lines that asserted `UnitLogCdf`'s weighted refusal, the refusal this brief lifts.
- *Held-out coverage* (`HeldOutCoverage_PrecisionWeights`, `BASS_QUALIFICATION=1`, about a minute; rerun on the committed code's build with the same numbers): f = sin(2πx₀) + x₁, with noise at each case's precision w, log-uniform on [¼, 4]: (0.2ζ plus, one time in seven, 0.8·Exp(1))/√w, a shape outside the model's Normal terms. Twenty fits of n 99, fifty held-out cases each with its weight, at α 0.1. **Coverage 0.9020** (nominal 0.900, se 0.0090 over fits). Mean width 1.407; 1 of 400 sets in two pieces; none open at the grid's ends.

**A measurement, no default proposed: coverage tilts with the weight.** By the weight's half, cases with w ≥ 1 cover 0.952 and cases with w < 1 cover 0.852, at nearly one width (1.401 against 1.412). The cases are exchangeable as (x, w, y) triples, so the guarantee is marginal, and it holds. The full step's score is the predictive density, which reads a case's height, not its tails: a set accepts y where N(y; f, σ²/w) clears a threshold shared across precisions. Its half-width in the case's own standard deviations, √(2(ln √w + c)), grows with w. So precise cases over-cover and noisy ones under-cover. A scale-free score would read each case at its own scale: the full step's predictive distribution (the PIT, which reads Φ(√w·e/σ)) or a split form. Whether weighted fits should prefer one is the coordinating session's call; measuring the distribution form's coverage by weight half on these fits would decide it.

**Left unbuilt, and why.**
- Case weights beside a noise, random effects, a curve or several outputs (refused above). Each needs its own reading: the product with a noise's latent weights, the block's whitening under replication, the stacked direct design and its support.
- Prediction's noise at a weighted new case: `Predict` reads σ² as a unit-weight observation's (BASS's page). A new case's precision weight could enter its predictive draws as the candidates' does.
- A flat-histogram run under case weights is not exercised at fit level. Chains, tempered ladders and SMC are.

**Found in passing, not fixed (outside this brief).** Under random intercepts, `BassCandidates.UnitLogCdf` reads `Working.Values`, the whitened response, against `BassDrawFit`'s f + u_g, while the retained pointwise terms read the response itself. So `PointwiseLogDensity(cdf: true)` on a random-intercept fit pairs the whitened response with an unwhitened fit. This predates this brief; the E21 owner should check it.

**Doc changes owed.**
- `docs/Maths/Regression/BASS.md`, "BASS: the weighted core": "Nothing fills the seam yet, and `Bass.Fit` stays Gaussian" becomes: case weights on the data (`CaseWeights`) fill it for Gaussian BASS on a scalar response. Add `CaseWeights` to the "Who hands it" row, beside aggregated rows and resampling, with the defaults and why (the meaning stated, `Weighted` by default, `Frozen` refused) and the refusals. "What reads the design" gains the candidate density: each new case's weight, read at the fit's meaning, with the replication CDF taken per copy. "Checked" gains the gates above.
- `docs/Maths/Regression/BASS.md`, the moving-response section's "`Bass.Fit` stays Gaussian": a latent block is still the consumer's, but a fixed weighting now reaches `Bass.Fit`.
- `docs/Maths/Regression/ConformalPrediction.md`, the "Not built yet" line "Candidates for a weighted BASS fit, which wait for such a fit to be reachable": built. Add each case's weight to the seam's inputs, and record the tilt measurement (marginal 0.902; 0.952 and 0.852 by weight half) beside the density score's scale reading.
- Issue records, for their owners: `issues/regression/bass/extensions-checklist.md` ("`Bass.Fit` takes no observation model, so a weighted run reaches a model and a chain, not a run") and `extensions-brief.md`'s slice 7 line: a weighted run now reaches `Bass.Fit` through `CaseWeights`.

### Slice 9 follow-up — CV+ with equal-sized folds

Built on `main` on 2026-11-09 in `dd2532c`, following the CV+ proofs in `33665cf` and the one-fold and order-reading update in `17ba5c9`. The coordinating review authorized this follow-up after signing off the first slice 9 increment. [API and statistical contract](../../docs/Maths/ConformalPrediction.md#cv-with-equal-sized-folds) · [Lean handoff](../../lean/Protolemmata/Conformal/conformal-chip-brief.md#reports)

`CrossConformal.CvPlus` consumes `FoldReadings(Response, HeldOut, FoldOf, Predictions)`: n training units in K positive equal-sized folds of m = n/K, a held-out prediction for every training unit, and K prediction rows at the test cases. Every residual and test prediction comes from that unit's own fold fit. Fold labels are dense from zero, and the layer checks the dimensions, equal counts and finite readings. Fitting, preprocessing and tuning remain the caller's responsibility.

The closed ends are the ℓ-th smallest lower and ℓ-th largest upper values, at ℓ = ⌊α(n + 1)⌋. The absolute-residual and endpoint-sort kernels are shared with jackknife+; its existing grids and defaults remain available. Ties are accepted, empty intervals are represented, and zero ℓ gives the whole line with coverage one. K = 1 is supported with an empty-training fit. K = n reduces to jackknife+ endpoints, including permuted singleton fold labels.

**Each bound has its own hypotheses.** `CvPlusBounds` returns all three numerical bounds, clipped at zero, and each `CvPlusBound` carries `Value`, `Theorem`, `RequiredAssumptions` and `Applies`.

| Bound | Coverage lower bound | Caller assertions |
|---|---|---|
| `Merged` | 1 − 2(ℓ + K − 1)/(n + K) | `ExcludesFoldAndTest` and `FoldTestExchangeability`, combined as `MergedRanks` |
| `OneFold` | 1 − ℓ/(m + 1) | The same assertions suffice; the theorem needs just one fold's exclusion and exchangeability with the test unit |
| `Tournament` | 1 − max(2ℓ + m − 2, 0)/(n + m) | `ExtendedFoldExchangeability` and `SymmetricFit`, combined as `FoldTournament` |

`Coverage` takes the largest applicable value, or null when no positive-level claim applies. `Both` supplies both assertion sets. The merged and one-fold bounds admit fits reading the other training units' order or positions. Each fold with the test unit must be exchangeable given the rest and the fitting randomness; its fit, preprocessing and tuning read neither. The tournament requires a joint law on n + m units invariant under fold and slot permutations with the matching action on fitting randomness, and a fit symmetric under that action. I.i.d. data with independent fitting randomness supply the extension. Exchangeability of n + 1 units alone supplies no such extension. The API's flags are caller assertions, not runtime certification of a data law. The one-fold API reuses the sufficient all-fold assertion instead of adding a separate chosen-fold contract.

**Finite controls.** There are 134,648 strict-count comparisons across tied endpoints, empty sets and dyadic scales near 10⁻¹⁵⁰ to 10¹⁵⁰. The merged and one-fold implications hold over 360 and 15,120 complete orderings at K = 2 and K = 3, respectively, using symmetric, first-unit and compact-order position-weighted fits. All 564 and 39,432 refusals obey both implications. Another 120 one-fold orderings check empty-training fits, ties and split endpoint/bound parity. The tournament covers above its bound over 2,880 orderings of n + m units with mean, median, nearest-neighbor and shared-priority fits; the unused test-fold filler never enters the public readings. All 16 assumption-flag combinations are checked. A K = 2 case selects the strictly stronger one-fold bound; at ℓ = 0 no assumptions are needed. K = n endpoint and applicable tournament-bound parity holds through n = 50 at four levels.

**The order-reading correction.** The public jackknife+ and n-fold docs now include the sufficient compact ordered-row contract proved by `jackknifePlus_orderReading`: the held-out place is closed up, and original positions, indices, masks and row-keyed randomness are absent. Ordered sweeps sharing a fixed seed can qualify. This supersedes the first increment's blanket exclusion of fixed index-order sweeps below. The symmetric-fit route remains available. Position-aware fits receive no such jackknife+ guarantee: the proof ledger's parity construction covers exactly one half at positive levels. CV+'s merged and one-fold hypotheses are separate and still permit position-reading fits outside their held-out fold. Tournament assumptions are unchanged.

**Opt-in coverage and width.** The committed qualification uses 500 independent replications and 2,500 test units per fit at n = 40, α = 0.1 and ℓ = 4. Its truth is skewed, heteroscedastic and curved. The symmetric least-squares line omits that curvature; a second fit returns the first training response. Standard errors use replication means to account for the five test units sharing a training set. The line supplies both assertion sets; the first-unit fit supplies only `MergedRanks`, enabling the merged and one-fold bounds.

| Fit | K | Measured coverage | SE | Mean width | Selected proved lower bound |
|---|---:|---:|---:|---:|---:|
| Symmetric line | 2 | 0.91160 | 0.00601 | 5.2643 | 0.80952, one fold |
| Symmetric line | 5 | 0.90720 | 0.00597 | 5.0626 | 0.70833, tournament |
| Symmetric line | 10 | 0.90640 | 0.00605 | 4.9846 | 0.77273, tournament |
| Symmetric line | 40 | 0.90600 | 0.00608 | 4.9431 | 0.82927, tournament |
| First-unit fit | 2 | 0.93000 | 0.00546 | 7.9900 | 0.80952, one fold |
| First-unit fit | 5 | 0.93080 | 0.00554 | 7.7873 | 0.64444, merged |
| First-unit fit | 10 | 0.92600 | 0.00538 | 7.4015 | 0.48000, merged |
| First-unit fit | 40 | 0.91320 | 0.00590 | 6.9858 | 0.00000, merged/one fold |

At the same 40-unit budget, split residual conformal uses 20 training and 20 calibration units. The line gives coverage 0.90200 (SE 0.00644), width 5.5234; the first-unit fit gives 0.91040 (SE 0.00622), width 7.1782. Its finite lower bound is 19/21 = 0.90476. These near-nominal measurements are fixture readings, not upgraded guarantees. The first-unit CV+ check explicitly excludes the tournament even where its numerical value is higher. Its supplied `MergedRanks` assertion permits position-aware fits and does not assert the separate compact order-reading jackknife+ contract; at K = n it therefore keeps the conservative bound selected by those flags.

**Cost.** Release calibrated thread-cycle clock, best of three after a warmup, n = 1,000 and 200 test cases. Endpoint inputs are built before that clock. The split comparator times residual calibration and endpoints around a zero center, with its inputs likewise prepared first.

| K | Linear fold fits and predictions | CV+ endpoint reduction | Per test case | Split calibration/endpoints |
|---:|---:|---:|---:|---:|
| 5 | 1.806 ms | 71.045 ms | 355.226 µs | 0.568 ms |
| 10 | 1.287 ms | 21.540 ms | 107.702 µs | 0.215 ms |
| 1,000 | 437.417 ms | 20.711 ms | 103.556 µs | 0.175 ms |

These timings are indicative on a shared machine, with runtime tiering effects, and price simple linear refits rather than BASS sampler refits. The final qualification also reran all 11 shared `ConformalCostTests` cases.

**Validation.** The selected precommit Edit gate passed 859/859. On a fresh `git archive dd2532c`, all four owned source/test blobs were matched to their Git objects before a serial Release build. Focused checks passed 16/16. With `BASS_QUALIFICATION=1` and `BASS_COST_CLOCK=cycles`, the combined focused/qualification/common-cost filter passed 32/32. All 16 gated bodies printed output: two CV+ coverage cases, three CV+ cost cases and eleven common-cost cases. Qualification is disabled for the Edit, Commit and whole-assembly runs. The commit-scope gate passed 874/874, including pins. The whole regression assembly passed 2,932/2,932 with zero skips on the same archive (14 m 14 s in the test runner).

Receipts: [source manifest](receipts/chip-e16-cvplus/results/source-dd2532c.json), [precommit Edit](../../artifacts/chip-e16-cvplus/results/edit-precommit.trx), [focused](../../artifacts/chip-e16-cvplus/results/focused-dd2532c.trx), [qualified](../../artifacts/chip-e16-cvplus/results/qualification-dd2532c.trx), [Commit](../../artifacts/chip-e16-cvplus/results/commit-dd2532c.trx), [whole assembly](../../artifacts/chip-e16-cvplus/results/slice-dd2532c.trx), [qualification settings](receipts/chip-e16-cvplus/results/qualification.runsettings). The final archive uses a distinct build directory; preliminary incremental-development receipts are superseded.

Reproduction from a fresh `dd2532c` archive in `$snapshot`, with output in `$build` and `$receipts` (below-normal process priority, `MSBUILDDISABLENODEREUSE=1`, `DOTNET_PROCESSOR_COUNT=2`):

```powershell
dotnet build "$snapshot/tests/Maths.Regression.Tests/Maths.Regression.Tests.csproj" -c Release "-p:ArtifactsRoot=$build/" --disable-build-servers -m:1
$env:BASS_QUALIFICATION='1'
$env:BASS_COST_CLOCK='cycles'
dotnet test "$snapshot/tests/Maths.Regression.Tests/Maths.Regression.Tests.csproj" -c Release "-p:ArtifactsRoot=$build/" --no-build --no-restore --disable-build-servers --settings "$receipts/qualification.runsettings" --filter 'FullyQualifiedName~CvPlusTests|FullyQualifiedName~CvPlusQualificationTests|FullyQualifiedName~CrossConformalTests|FullyQualifiedName~ConformalCostTests'
$env:BASS_QUALIFICATION='0'
./scripts/test-select.ps1 -SourceRoot $snapshot -Commit dd2532c -Scope Commit -ArtifactsRoot $build -Threads 2 -NoBuild
$env:DOTNET_PROCESSOR_COUNT='3'
./scripts/test-select.ps1 -SourceRoot $snapshot -Scope Slice -ArtifactsRoot $build -Threads 3 -NoBuild
```

**Scope.** Equal-sized folds only. The cost qualification uses simple linear fits, with fitting/prediction time separated from endpoint reduction and split calibration. Cross-conformal predictive distributions have no claim here. E21 grouping and BASS prediction, E34 readouts, sampler trajectories, and Lean source remain outside this change.

### Slice 9 — The proved cross-fitted forms and asymmetric CQR

Built on `main` on 2026-11-09 in `3fbe41b`, following Claude's `4084c55` proof contract. The implemented rules are `JackRefuses`, `AcceptsPairwise`, `jackknifePlus_coverage`, `crossFit_coverage` and `cqrAsym_coverage_half`. [Lean handoff](../../lean/Protolemmata/Conformal/conformal-chip-brief.md) · [Statistical choices and API](../../docs/Maths/ConformalPrediction.md#cross-fitted-sets-and-asymmetric-cqr-slice-9)

- `CrossConformal.JackknifePlus` reduces leave-one-out prediction rows to closed endpoints. With n training units and m = n + 1, it takes k = ⌊αm⌋ and the k-th smallest lower and k-th largest upper endpoints. At k = 0 it returns the whole line.
- `CrossConformal.Pairwise` accepts when the candidate is strictly stranger than fewer than m − k held-out units, under each unit's own fit. Ties count for neither side. Its `ConformalSet.Form` is `CrossFitted`, and the result reports the applicable finite bound. `CrossConformal.Density` reduces paired held-out and candidate seams one fit at a time. Draws align within each fit and can differ between fits.
- For both forms the reported bound is 1 − (2k − 1)/m at positive k, clipped at zero, and 1 at zero k. It requires exchangeable units and fitting symmetry under relabeling, including a relabeling of the fit's randomness. The library reads supplied fits and cannot certify their symmetry. Fixed index-order sampling with a fixed seed falls outside the supplied contract.
- `SplitConformal.AsymmetricQuantiles` calibrates y − q_lo and q_hi − y separately at α/2. It returns both tail thresholds, levels and miscoverage bounds, with total bound 1 − (k_lo + k_hi)/m. A zero tail level gives an infinite endpoint.

**Finite controls.** The public pairwise rule is checked on all 760 comparison patterns through m = 4, with ties and infinite scores, and 400 seeded matrices through m = 30. Public jackknife+ intervals agree with 74,068 strict endpoint comparisons across dyadic response scales near 10⁻¹⁵⁰ to 10¹⁵⁰. The 960 permutation cases exercise a stable fit, a sample-size-unstable fit and a shared priority order carried with the units. The pairwise residual set stays inside the jackknife+ interval. A symmetric regular-tournament fixture attains the finite bound, including coverage 17/20 at n = 19, α = 0.1. Asymmetric CQR agrees with 23,554 count comparisons, including 178 empty intervals. Its controls cover 360 exact-tail orderings with distinct scores, ties, crossing quantiles, unseen strata and small calibration sets.

**Opt-in coverage and width.** On 500 independent replications of the skewed, heteroscedastic truth used by the earlier conformal checks, 2,500 held-out units gave the readings below. The regression model is a symmetric least-squares line and omits that truth's curvature and skew. Standard errors use replication means, accounting for the five test units sharing their training data.

| Form | Measured coverage | SE | Mean width | Proved lower bound |
|---|---:|---:|---:|---:|
| Jackknife+, n = 39 | 0.89960 | 0.00644 | 4.9382 | 0.825 |
| Pairwise residual, n = 39 | 0.89960 | 0.00644 | 4.9381 | 0.825 |
| Split residual, 20 training + 19 calibration | 0.90280 | 0.00621 | 5.6187 | 0.900 |

Each uses the same total budget of 39 observed units. Jackknife+ uses 39 leave-one-out fits. The pairwise width is read on a grid of spacing 0.025. No set reached a grid edge. These near-nominal coverages are measurements of this fixture.

On a separate 500-replication fixture with both quantile fits shifted below the truth and 39 calibration units, asymmetric CQR missed below at 0.05480 (SE 0.00463) and above at 0.05000 (SE 0.00463). Coverage was 0.89520 (SE 0.00619), with width 5.4924. Its bounds are 0.05 per tail and 0.90 total. Symmetric CQR on the same fits missed above at 0.10200 (SE 0.00622), with width 9.3126.

**Cost.** Indicative Release calibrated thread-cycle readings, best of three, with density matrices built before the clock. At n = 1,000 and G = 200, the paired density reductions retain one fold's matrices at a time. The table separates the measured cost of simple linear leave-one-out fitting from the reductions.

| Draws per fit | n-fold density reductions and ranks | Per candidate | Split density | 1,000 linear refits and predictions | Jackknife+ endpoint sort |
|---:|---:|---:|---:|---:|---:|
| 1,000 | 1,902.32 ms | 9,511.58 µs | 13.77 ms | 195.99 ms | 0.7253 ms |
| 4,000 | 10,820.69 ms | 54,103.44 µs | 91.62 ms | 349.16 ms | 0.3646 ms |

A separate opt-in cost case (`bccb661`) measured calibration plus 200 test intervals at n = 1,000. Symmetric CQR took 1.3005 ms and asymmetric CQR 2.0863 ms, best of three on the same cycle clock, with the quantile readings built before the clock. [CQR cost receipt](../../artifacts/chip-e16-s9/results/cqr-cost.trx)

**Qualification flags and receipts.** The main run uses a clean archive of `3fbe41b`, including the committed E34 runtime work in `41727ac`. The focused contract filter passed 11/11. The selected classes plus pinned tier passed 843/843, both with `BASS_QUALIFICATION=0`. The opt-in filter passed 26/26 with `BASS_QUALIFICATION=1` and `BASS_COST_CLOCK=cycles`.

The whole regression assembly passed 2,846/2,846 on that archive with qualification off, including the slow tier, at three xUnit workers. The test duration was 26 minutes. The other main runs used two workers.

All 15 gated cases in that filter printed their measurement bodies. These were `HeldOutCoverageAndWidth_JackknifePairwiseAndSplit` (1), `AsymmetricCqr_EachTailUnderWrongQuantileFits` (1), `DensityCost_ThousandFitsAgainstTheSplit` (2), and the existing `ConformalCostTests` cost rows (11). The supplementary CQR cost body ran and printed separately, with the same opt-in flags. In flag-off runs those cases return early and still count as passed.

The supplemental test commits keep the production source of `3fbe41b`. The `a983e1a` CQR filter passed 4/4, then its selected class plus pinned tier passed 40/40. The `bccb661` opt-in CQR cost filter passed 1/1, then its selected class plus pinned tier passed 41/41. All reported counters have zero failures and zero skipped.

| Receipt | File |
|---|---|
| Focused contracts on `3fbe41b` | [11/11](../../artifacts/chip-e16-s9/results/focused-3fbe41b.trx) |
| Commit selection on `3fbe41b` | [843/843](../../artifacts/chip-e16-s9/results/commit-3fbe41b.trx) |
| Opt-in measurements on `3fbe41b` | [26/26](../../artifacts/chip-e16-s9/results/qualification-3fbe41b.trx) |
| Whole regression assembly on `3fbe41b` | [2,846/2,846](../../artifacts/chip-e16-s9/results/slice-3fbe41b.trx) |
| CQR rounding regression | [4/4](../../artifacts/chip-e16-s9/results/rounding.trx) |
| Rounding commit selection | [40/40](../../artifacts/chip-e16-s9/results/commit-a983e1a.trx) |
| Supplementary CQR cost | [1/1](../../artifacts/chip-e16-s9/results/cqr-cost.trx) |
| CQR cost commit selection | [41/41](../../artifacts/chip-e16-s9/results/commit-bccb661.trx) |

The commit receipt also shows new methods missing from the older slow profile. The two E21 drawn-scale ladder and SMC methods sum to 310.0 s and 191.4 s respectively. E34's multi-output SMC enumeration method adds 53.8 s. These observations went to the coordinator for the next slow-tier refresh. The older listed warm-start method ran through the selector's intentional pinned-tier override.

**Commands.** Builds and tests used Release, below-normal priority and disabled build-server reuse. The main archive excludes `r`, following the handoff convention.

```powershell
$run = Join-Path (Get-Location) 'artifacts/chip-e16-s9'
$src = Join-Path $run 'snap-3fbe41b'
$build = Join-Path $run 'head-build'
$results = Join-Path $run 'results'
$project = Join-Path $src 'tests/Maths.Regression.Tests/Maths.Regression.Tests.csproj'
(Get-Process -Id $PID).PriorityClass = 'BelowNormal'
$env:MSBUILDDISABLENODEREUSE = '1'
$env:DOTNET_PROCESSOR_COUNT = '2'
$env:BASS_QUALIFICATION = '0'
$env:BASS_COST_CLOCK = 'cycles'
dotnet build $project -c Release "-p:ArtifactsRoot=$build/" --disable-build-servers -m:1
$selection = @{ SourceRoot = $src; Commit = '3fbe41b'; ArtifactsRoot = $build; NoBuild = $true }
dotnet test $project -c Release "-p:ArtifactsRoot=$build/" --no-build --disable-build-servers -m:1 --settings "$results/qualification.runsettings" --filter 'FullyQualifiedName~CrossConformalTests|FullyQualifiedName~AsymmetricCqrTests' --logger 'trx;LogFileName=focused-3fbe41b.trx' --results-directory $results
./scripts/test-select.ps1 @selection -Scope Commit -Threads 2
$env:BASS_QUALIFICATION = '1'
$filter = 'FullyQualifiedName~CrossConformalTests|FullyQualifiedName~AsymmetricCqrTests|FullyQualifiedName~CrossConformalQualificationTests|FullyQualifiedName~ConformalCostTests'
dotnet test $project -c Release "-p:ArtifactsRoot=$build/" --no-build --disable-build-servers -m:1 --settings "$results/qualification.runsettings" --filter $filter --logger 'trx;LogFileName=qualification-3fbe41b.trx' --results-directory $results
$env:BASS_QUALIFICATION = '0'
$env:DOTNET_PROCESSOR_COUNT = '3'
./scripts/test-select.ps1 @selection -Scope Slice -Threads 3
```

**Numerical choices.** At subnormal α, α/2 can round to zero. Its ranks remain zero at every supported sample size, so the implementation uses the smallest positive representable tail level there. If finite responses and quantiles produce an overflowing CQR conformity difference, the new asymmetric path asks for a common rescaling instead of constructing misleading endpoints.

**Handoff wording corrected.** The implemented bounds are B_asym = 1 − 2⌊αm/2⌋/m and B_sym = 1 − ⌊αm⌋/m. Their relation is B_asym − (1 − α) ∈ [0, 2/m), while B_asym − B_sym is 0 or 1/m. At α = 0.1 and m = 31 the bounds are 29/31 and 28/31. The asymmetric bound's excess above nominal 0.9 is 0.035484, greater than 1/31. The handoff's prose claimed an excess of at most 1/m above nominal. The coordinator has the correction for the Lean notes' owner. The exact level statements and the proved lower bound are unaffected.

The noninteger αm case is covered by `Rounding_NonintegerAlphaTimesM_AgainstBothFiniteBounds` in `a983e1a`. Its snapshot retains the production source of `3fbe41b`. [Supplemental CQR receipt](../../artifacts/chip-e16-s9/results/rounding.trx)

**Next contract.** Claude committed the CV+ proof in `33665cf` while this increment's whole-assembly run was in progress. The [handoff](../../lean/Protolemmata/Conformal/conformal-chip-brief.md#2026-11-09--slice-9s-follow-up-cv-proved-both-bounds) provides both merged-rank and fold-tournament bounds. The merged-rank bound admits a fit without symmetry that excludes its held-out fold and test unit. The tournament bound needs a symmetric fit and the larger exchangeability structure supplied by i.i.d. data. A bound is available only when its own assumptions apply. CV+ is the separate implementation after coordinating verification of this increment. Cross-conformal predictive distributions remain a separate question.

### Slice 8 — Grid-free split quantiles and refined set boundaries

Built on `main` on 2026-11-09 in `331f9e7`. `SplitConformal.CalibrateDistribution` sorts each stratum's calibration PITs once. The grid API and the new grid-free API use the same weighted two-tail scores. At a fixed test input the layer accepts the fit's inverse CDF or brackets a root through its candidate seam. Every probe checks the split form, the draw metadata, latent treatment, stratum and requested response. The statistical choices and API are on [the conformal page](../../docs/Maths/ConformalPrediction.md).

For m units including the candidate, the equi-tailed set inverts the k-th and (m − k)-th PITs, with k = ⌊αm/2⌋. It covers the whole line at k = 0. Quantiles invert the mean-over-τ rank used by the grid. When the test CDF is flat at a calibration PIT, the tied rank decides which end of the plateau to return. Numerical set ends round outward by at most the requested absolute tolerance. `ConformalSet.RefineBoundaries` returns an accepted and a rejected witness in each grid cell whose endpoint decisions differ.

**Numerical and statistical gates**

- On a weighted mixture with a zero-weight draw and three strata, the largest difference across the closed quantiles and set ends was 0.00190685 against a public grid of step 0.002. The BASS predictive inverse agreed with the root of its fitted mixture CDF to 10⁻⁸. BAPS, Gaussian BARS, the hierarchical P-spline and both BGP scale models put the median crossing between the public grid ranks on either side. The split-only restrictions are exercised alongside those reads.
- On a disconnected density set, each of the split plain, full plain and full inflated rules bracketed all four analytic crossings. The widest bracket was 7.451 × 10⁻⁹ at a requested 10⁻⁸. Each accepted witness exceeded the level and each rejected witness was at or below it under a fresh singleton evaluation of the applicable rule.
- On 500 independent training and calibration samples from skewed, heteroscedastic truth, 3,602 of 4,000 held-out responses fell inside the nominal 90% interval (90.05%, mean width 5.9698). The randomized PIT had KS distance 0.0159. The ten bins induced by the closed deciles held 9.20% to 11.00% each, with a per-bin SE of 0.55 percentage points accounting for shared calibration. The ungated controls also cover duplicate PITs, upper-tail ranks beyond ordinary CDF precision, unseen strata, a saturated PIT, plateaus and a tolerance below floating-point spacing.

**Cost on the commit.** Indicative calibrated thread-cycle times in Release, best of five. The median and two 90% ends share one calibration of n = 1,000. The test input has a standard Gaussian predictive with a known native inverse. Both paths include candidate-seam construction. Calibration time is listed separately.

| Draws | Calibrate PITs | Grid, G = 200 | Fitted inverse, 7 probes | CDF roots, 89 probes |
|---:|---:|---:|---:|---:|
| 1,000 | 18.46 ms | 53.40 ms | 1.14 ms | 21.72 ms |
| 4,000 | 73.06 ms | 233.43 ms | 5.74 ms | 91.66 ms |

For the split density rule at n = 39, S = 40 and tolerance 0.001, a 17-point grid plus 18 singleton probes cost 1.16 ms. An 8,001-point grid at spacing 0.001 cost 95.14 ms. The qualification receipt also includes the full step's cost breakdown at both draw budgets for independent draws and two chains, the split density cost rows and both distribution cost rows.

**Limits of the readings.** Refinement locates a crossing in each observed transition cell. Equal endpoint decisions can conceal further components, and an accepted grid edge gives no bound outside the grid. The cost of a root includes all calls to the supplying engine. BGP can refit away from its retained test inputs, and a split density callback using the public grid method rereads calibration on each probe.

**Verification on `331f9e7`.** All runs used Release. The code commit includes E34's committed `ac9726b` in its ancestry. The TRX counters for each run show zero failures and zero skipped.

| Check | Passed | Qualification | Receipt |
|---|---:|---|---|
| Edit selection, 44 classes selected | 754/754 | Off | [Edit TRX](../../artifacts/chip-e16-s8/results/edit-331f9e7.trx) |
| Commit selection with the pinned tier | 764/764 | Off | [Commit TRX](../../artifacts/chip-e16-s8/results/commit-331f9e7.trx) |
| Slice qualification filter | 34/34 | On, cycle clock | [Qualification TRX](../../artifacts/chip-e16-s8/results/qualification-331f9e7.trx) |
| Whole regression assembly | 2,755/2,755 | Off | [Slice TRX](../../artifacts/chip-e16-s8/results/slice-331f9e7.trx) |

Twelve gated cases executed their measurement bodies in the qualification run and printed readings in the TRX. They were `HeldOutPitAndClosedDeciles_AreUniformOnSkewedTruth` (1), `ExactEndsCostRow` (2), `RefinedBoundariesCostRow` (1), `CostRow` (4), `SplitCostRow` (2) and `DistributionCostRow` (2). The off runs count successful early returns for those gates. The whole-assembly TRX also records exact agreement with all 1,968 public split grid ranks and 3,122 full grid ranks in the existing seeded reference sweeps.

**Commands for the committed snapshot.** These use the repository's PowerShell selector. The source came from `git archive 331f9e7`, extracted with `--exclude=r`. All seven changed source and test blobs matched the commit. The run used below-normal priority and `MSBUILDDISABLENODEREUSE=1`. On the initial development snapshot, a build with multiple MSBuild workers exited before tests started. Building explicitly with one worker and testing that DLL supplied the recorded runs.

```powershell
$root = Join-Path (Get-Location) 'artifacts/chip-e16-s8'
$src = Join-Path $root 'snap-331f9e7'
$build = Join-Path $root 'head-build'
$results = Join-Path $root 'results'
$project = Join-Path $src 'tests/Maths.Regression.Tests/Maths.Regression.Tests.csproj'
(Get-Process -Id $PID).PriorityClass = 'BelowNormal'
$env:MSBUILDDISABLENODEREUSE = '1'
$env:DOTNET_PROCESSOR_COUNT = '2'
$env:BASS_QUALIFICATION = '0'
$env:BASS_COST_CLOCK = 'cycles'
dotnet build $project -c Release "-p:ArtifactsRoot=$build/" --disable-build-servers -m:1
$selection = @{ SourceRoot = $src; Commit = '331f9e7'; ArtifactsRoot = $build; NoBuild = $true }
./scripts/test-select.ps1 @selection -Scope Edit -Threads 2
./scripts/test-select.ps1 @selection -Scope Commit -Threads 2
$env:BASS_QUALIFICATION = '1'
$filter = 'FullyQualifiedName~ConformalExactEndsTests|FullyQualifiedName~ConformalCostTests|FullyQualifiedName~EngineCandidateTests.Baps_AtTrainingCases|FullyQualifiedName~EngineCandidateTests.Bars_AtTrainingCases|FullyQualifiedName~EngineCandidateTests.Hierarchical_AtTrainingCases|FullyQualifiedName~EngineCandidateTests.Bgp_TheSplitForm'
dotnet test $project -c Release "-p:ArtifactsRoot=$build/" --no-build --disable-build-servers -m:1 --settings "$results/qualification.runsettings" --filter $filter --logger 'trx;LogFileName=qualification-331f9e7.trx' --results-directory $results
$env:BASS_QUALIFICATION = '0'
$env:DOTNET_PROCESSOR_COUNT = '3'
./scripts/test-select.ps1 @selection -Scope Slice -Threads 3
```

### 2026-11-09 — Slice 7, the seams completed: each names the forms it serves

Built at the user's decisions of 2026-11-09 (6 to 8 above), relayed by the main session.

**Landed** on `main`, nothing pushed (the docs and E16's line are in the log beside this report):
- `e3071cf`, with its tests in `EngineCandidateTests`:
  - `ConformalForms` (`Full`, `Split`, `Both`) and `CandidateLogDensity.Forms` and `FormsReason`, both by default;
  - every entry point holds a seam to its forms, quoting the engine's reason;
  - the hierarchical P-spline's and BGP's seams name the split form alone.
- The page: the forms on "The seams", slice 5's table updated, and the API row.

**Gates.**
- *Each engine's forms:* BAPS's and BARS's seams serve both, and the full step takes them as before. The hierarchical P-spline's serves the split form, and the full step now refuses it beside point units and whole-curve units alike, quoting the shared design. BGP's serves the split form, and its refusal says why its terms are leave-unit-out.
- *The layer holds every entry point to the flag:* the full step, its distribution and its scores; the split form and its distribution; an inductive stream's density and centrality scores. Out-of-range flags are refused.
- *The conformal tests:* 78 of 78. *Per commit,* the selected classes and the pinned tier: 2,141 of 2,141. *At the slice's end,* the whole assembly at e3071cf: 2,736 of 2,736.

**The weighted BASS fit's candidates: not built, blocked on a missing entry point** (corrected 2026-11-10 at the user's word: "it seems like some things are missing that are expected, not that there's nothing to build").
- *The core holds weighted fits.* `WorkingResponse` takes weights with their meaning, precision or replication, and `BassWeightedCoreTests` gates them at `BassModel`'s level.
- *The path to them is missing.* `Bass.Fit` builds its model from the config alone and takes no weights. BASS's page (`WeightMeaning`) names who hands a weight: aggregated rows and resampling (E18) for replication; a latent scale, Pólya–Gamma or a mixture for precision. Aggregated rows have no entry point, so no fit a user can make carries weights.
- *Once a weighted fit is reachable, the candidates take each case's weight.* Each candidate reads the training row's term: N(y; f, σ²/w) under precision weights, w copies of the Gaussian term under replication, with the CDF read alike. Until then the seam's refusal stays as a guard.
- *My first record called this "nothing to build".* That mistook a missing entry point for an absence of need.
- *Doc changes owed:* the page's "Not built yet" line on these candidates should name the missing `Bass.Fit` entry point, not a fit "to be reachable".
- *Built 2026-11-10* (`7343bc8`, `c1bbf51`): `CaseWeights` on `BassData`, and each new case's weight on the candidates. See the report of that date above.

**Decided in building it.**
- *The flag is three values.* A seam that leaves a form out gives its reason, which the refusal quotes; a refusal without one says the engine gave none.
- *An inductive stream's scores need the split form,* since they score as its calibration does.

**Next:** slice 8, exact ends and quantiles. Slices 9, 10 and 12 wait for their Lean.

### 2026-11-08 — Slice 6c, the alarms (E35's Ville, CUSUM and Shiryaev–Roberts); E35 built

Built after the main session verified 6b. With it E35 is built as planned: the p-values (6a), the bets and capital (6b), and the alarms.

**Landed** on `main`, nothing pushed (the docs, E16's and E35's lines are in the log beside this report):
- `befb449` the following, with `MartingaleAlarmTests`:
  - `AlarmProcedure`, `MartingaleAlarm` and `AlarmGuarantee`;
  - the martingale's alarms, with the procedures fired at each step.
- The subsection on `ConformalPrediction.md`, with which procedure to read.

**The gates, as enforced.**
- *The alarms read only the past:* each statistic is a recursion on the capital so far.
- *Ville's guarantee is stated as proved,* CUSUM's and Shiryaev–Roberts's as posed, in each alarm's `Guarantee`.
- *Each alarm restarts after it fires,* so its guarantee is per run.
- *Refused:* a threshold at or below 1, an undefined procedure, and a bet whose log is NaN or +∞.

**Gates** (numbers on the page).
- *Each statistic against its definition from the path,* within 10⁻⁹, restarts included.
- *Ville's false alarms at most 1/c* over 10,000 streams, for five bets and two c: 0.021 to 0.077 at c 10, and 0.184 at c 4.
- *CUSUM and Shiryaev–Roberts* meet their posed bounds in every reading. Shiryaev–Roberts with Simple Jumper sits near its bound: a mean time of 21.7 at c 20.
- *Delays after a change,* per bet and procedure. Simple Jumper detects a 1-sd shift in the mean, or the sd raised to 1.6, in 98 to 99% of streams under Ville's alarm (median 110 to 119), and CUSUM in 97 to 98% (75 to 81).
- *A fitted engine's scores:* a BARS fit and a BASS fit, whose noise rose from 0.2 to 0.5. Ville's false alarms are 0.02 and 0.01 at c 20, within 1/20. The change is detected in 98 and 99% of streams, median delay 36 and 38.
- *Cost:* 46 ns per arrival for the three alarms' recursions. The arrival, rank included, costs 1.5 to 4.3 µs with the load.
- *Per commit,* the selected classes and the pinned tier: 77 of 77. *At the slice's end,* the whole assembly at befb449: 2,635 of 2,635.

**Defects found:** none. The cost row's first form timed the bare martingale after the alarmed one and read a 1.9 µs difference. A probe with the order swapped put the alarms at 0.12 µs and the rest down to load, so the row now times the alarms' recursions alone.

**Decided in building it.**
- *Every procedure restarts after an alarm.* That is Vovk's repeated form for CUSUM and Shiryaev–Roberts, and for Ville's a fresh test of the stream since.
- *A monitor takes several alarms on one capital.*
- *Shiryaev–Roberts's threshold is documented as a mean run length.* Over 500 arrivals at c 100 it alarms falsely in every stream with Simple Jumper, as its guarantee allows.

**Not built, from E35's plan:** the full form and one monitor per stratum, each with the decision it would change, recorded in the plan.

### 2026-11-08 — Slice 6b, the bets and the capital (E35's M-2)

Built after the main session verified 6·0 and 6a.

**Landed** on `main`, nothing pushed (the docs and E16's line are in the log beside this report):
- `206eca9` the following, with `BettingCapitalTests`:
  - `IBettingFunction` and its four members: power, the power mixture, Simple Jumper and the plug-in;
  - `BettingCapital`;
  - `ConformalTestMartingale`, which joins the ranks and the capital, with each step's ln f, the path and the e-value.
- The subsection on `ConformalPrediction.md`.

**The gates M-2 sets, as enforced.**
- *Each bet is placed before its p-value is seen:* `LogBet`, then `Observe`. The tests check that each member's curve is a function of the past alone.
- *Each bet integrates to at most 1:* to 1 for power, the mixture and Simple Jumper, and to at most 1 for the plug-in, scaled down by 10⁻¹² past its rounding.
- *The capital starts at 1:* ln S₀ = 0.
- *Refused:* a p-value outside (0, 1], and a bet that is NaN or +∞.

**Gates** (numbers on the page).
- *Integrals,* after random histories, uniform and skewed: within 2.3·10⁻¹⁴ of 1, the plug-in at most 1 within 1.1·10⁻¹⁰.
- *Past only:* each curve before pₙ is a fresh bet's fed p₁, …, pₙ₋₁.
- *Against references:* Simple Jumper against Algorithm 1 to 10⁻¹² at four J; the mixture against its integral over ε to 10⁻¹¹; the plug-in against a rebuild from its counts.
- *E S₂₀ = 1 end to end:* 0.9987 (se 0.0017) for power 0.9, 0.9997 (0.0025) for Simple Jumper and 0.898 (0.068) for the plug-in, over 100,000 streams.
- *Broken bets:* a peeking bet and one integrating to 1.1 reach 10 in every one of 2,000 streams, where Simple Jumper reaches it in 0.0075 of them.
- *After a change,* Simple Jumper and the plug-in rise from about −2 to 36 and 37 in ln S. The mixture and power 0.2 decay before the change, to −63 and −165, and do not recover: the reason for 6c's CUSUM and Shiryaev–Roberts.
- *Cost:* 1.8 to 3.3 µs per arrival with its rank, under load.
- *Per commit,* the selected classes and the pinned tier: 57 of 57. *At the slice's end,* the whole assembly at 206eca9: 2,576 of 2,576, nothing skipped.

**Defects found:** none. Two of my test bounds were wrong and were corrected.
- xUnit's 10-decimal rounding split two values that agree to 2·10⁻¹⁴; the comparison is now relative.
- A rise of 100 in capital after the change, asked of every bet, failed for the bets that had decayed under the null. Only the adaptive bets are held to it now, and the fixed bets' decay is recorded instead.

**Decided in building it.**
- *Bets give ln f,* so a power bet at a tiny p stays exact.
- *The plug-in's density is binned,* 100 bins by default, with its rows cached between refits, O(bins) per arrival. Fedorova et al.'s kernel density over every p-value would cost O(n).
- *The plug-in's uniform floor is an option, 0 by default,* the paper's pure kernel bet.

**Next:** 6c, the alarms and their readouts.

### 2026-11-08 — Slice 6a, the smoothed p-values of a stream (E35's M-1)

Built at the user's go-ahead for E35, the main session having verified slice 5.

**Landed** on `main`, nothing pushed (the docs and E16's line are in the log beside this report):
- `415c08c` the following, with `SmoothedPValueTests`:
  - `SmoothedPValues`, the p-values;
  - `OrderStatisticTree`, its ranks;
  - `StreamScores`, a fixed fit's density, standardized residual and PIT centrality.
- The section on `ConformalPrediction.md`, "Conformal test martingales".

**The gates M-1 sets, as enforced.**
- *Ranked against the whole stream:* one tree holds every score, and there is no calibration-set entry point.
- *A fresh θₙ per step,* in (0, 1], from the stream's own seeded generator.
- *Each score taken once, never revised:* the inductive form's fixed fit, scored as the split form scores its calibration.
- *The stream's exchangeability* is the caller's premise, stated on the type.

**Gates** (numbers on the page).
- *M-1's finite core, by enumeration:* seven data sets of four and five units, ties included. Every rank sequence has probability 1 over the orderings, to 10⁻¹².
- *i.i.d. uniform:* 400 streams of 50, with continuous, three-valued and all-tied scores, pooled and per position, and the lag-one correlation. Also uniform on a fixed N(μ, 1) fit's density and centrality scores of a t₃ stream.
- *Each premise broken, beside the type:*
  - a fixed calibration set: P(both lowest) is 1/10 against the stream's 1/20, exact;
  - CB-1's tie rule: 1 at every step of an all-tied stream;
  - scores kept from earlier fits: KS 0.057 against a fixed score's 0.0067.
- *The tree against a brute-force count,* ties, zeros and infinities included, its height within the AVL bound.
- *Cost:* about 0.4 µs per arrival warmed at n 10⁵, on par with `SortedSet`. The loaded cost row read 0.9 to 1.6 µs, and three-valued scores cost 25 to 30 ns.
- *Per commit,* the selected classes and the pinned tier: 52 of 52. *At the slice's end,* the whole assembly at 415c08c: 2,551 of 2,552, the one skipped the GBP test.

**Defects found: one of mine, before the commit.** The tree's first form wrote a child back with `_left[node] = Insert(...)`. C# reads the array before the call, so when an insert below grew the arrays, the write went to the old one: ranks ran one short, and p rose above 1. The tests caught it. The child is now written back after the call, and the comment there says why.

**Decided in building it.**
- *One walk ranks and inserts.* Ranking first and then inserting walked the tree twice.
- *The nodes are one array of structs,* not six parallel arrays. That halved the cost: each level's fields sit in one cache line.
- *θ lies in (0, 1],* so p never reads 0 and a power bet's p^(ε − 1) stays finite.
- *The centrality is min(F̂, 1 − F̂)* from the two tails, which orders units as the PIT's distance from ½ does, reversed. It is the two-sided choice; the PIT itself would test one direction.

**Next:** 6b, the bets and the capital.

### 2026-11-08 — Slice 6·0, BGP's seam at the run's test points

Built at the main session's note of 2026-11-08, slice 5's first follow-up, as slice 6's prelude.

**Landed** on `main`, nothing pushed:
- `4dd1dfb`: the run keeps each draw's fit variance of f, unscaled, beside the draw, as a reference shared by the draws of one fit with no copy (`BgpDraw.FitVariance`). It keeps its test inputs on the result (`BgpResult.TestInputs`).
- At those inputs, any subset in any order, the seam reads the kept variance and the mean instead of refitting. Elsewhere it refits as before.
- The page's BGP lines and cost row are updated in the docs commit beside this report.

**Gates** (`EngineCandidateTests`).
- *Bit for bit the refit path's,* for a known σ² and an unknown scale, the log CDF included, on a reordered subset of the test points with one point twice.
- *It reads what was kept:* doubling the kept variance moves the density. A case off the test points refits.
- *Cost* at n 400, S 1,000, one grid of 200, cycle clock, indicative: 3.5 ms per grid against 16.0 s per call by refitting in the same run (20.1 s in slice 5's row).
- *Per commit,* the selected classes and the pinned tier: 155 of 156, the one skipped the GBP test.

**Decided in building it.**
- *An unscaled copy is kept,* not the scaled one divided by s², whose quotient would differ from the refit in the last bit.
- *Cases match the test points on their values exactly.* A resumed run cannot change its test points: they are in its fingerprint.

### 2026-11-08 — Slice 5, the other engines' candidate seams

Built at the user's go-ahead of 2026-11-08 ("Yes, start it", relayed by the main session).

**Landed** on `main`, nothing pushed (the docs and E16's line are in the log beside this report):
- `eef6f47` the seams, with `EngineCandidateTests`:
  - BAPS's and BARS's, which both forms take, with the training units' log CDF on `cdf: true` for the full step's predictive distribution;
  - the hierarchical P-spline's and BGP's, which the split form takes;
  - shared pieces: `NoiseLogCdf`, moved from BASS's seam; `CandidateSeam`; `IWeightResampler.MarginalLogCdf`, a default member; `BarsNoise.LogCdfAt`; and `GpRegression.PredictiveLogDensities`.
- The section on `ConformalPrediction.md`, with its table of which form takes each seam and why.

**Per engine.**
- *BAPS:* N(y; z_tβ, σ²) at the new cases' design rows.
- *BARS:* every observation model it scores in sample, each observation's latent integrated as in training. A value outside the support takes −∞. The Gibbs loss is refused, as its pointwise density is.
- *The hierarchical P-spline:* a new point on a seen curve, conditioned on its drawn curve, as a unit that splits its replicate is.
- *BGP:* the GP's predictive given the training data at each draw (t, η), f, the mean and s² integrated.
- *The changepoint posterior:* no seam (below).

**Gates** (`EngineCandidateTests`; coverage and cost with `BASS_QUALIFICATION=1`).
- *At the training cases, the pointwise column, bit for bit, largest gap 0.* BAPS, and BARS under all nine of its models and noises: Gaussian, t weights, contaminated, t noise, quantile, skewed GIG, skewed GBP, Poisson and Whittle. Each log CDF and latent report matches too, and the full step takes each pair.
- *The hierarchical P-spline's,* against the column of a unit that splits its replicate. Beside whole-curve units the full step refuses it.
- *BGP's, against its pointwise content:* a fit on the data less one case gives, at that case, the leave-one-out term of the fit on all of them, to 1.7·10⁻¹⁴ (known σ²) and 1.6·10⁻¹⁴ (unknown scale). The full step refuses the seam, and the split form takes it.
- *Each law's log CDF against its density's integral,* a seeded sweep over σ² from e⁻⁴ to e⁴, the body and both tails out to where the density falls 10⁻¹². The closed forms agree to 10⁻¹² in the body and 3·10⁻¹⁰ in the tails. The tabled mixtures agree within their tables' stated 10⁻⁵: 6.8·10⁻⁶ at worst in the body (GBP without a drift), 6·10⁻⁷ in the tails.
- *Held-out coverage from truth outside every model:* f = sin(2πx) + x/2 with noise 0.2ζ plus, one time in seven, 0.8·Exp(1). There are 20 fits each, m 100 and k 10, 50 held-out cases per fit, and the se is the fits' spread.

| Engine and form | Coverage (se) | Mean width | Sets of two intervals |
|---|---|---|---|
| BAPS, the full step | 0.913 (0.012) | 1.11 | 0 of 400 |
| BARS (Gaussian), the full step | 0.917 (0.012) | 1.10 | 1 of 400 |
| The hierarchical P-spline, the split form | 0.902 (0.012) | 1.23 | 0 of 400 |
| BGP (unknown scale), the split form | 0.903 (0.014) | 1.12 | 0 of 400 |

  - BARS's one split set, at x = 0.002, reaches the grid's end: [−3.00, −2.12] ∪ [−0.52, 0.68], with ESS 2.2 and the inflated-set gap at 0.82.
  - My first gate asserted that no set reaches the grid's end, slice 1's gate for BASS. At BARS's edge that overstates what the draws give, so the gate now asserts that every such set carries the gap's flag, 0.5 or more.
- *The cost row* (S 1,000, one grid of 200, cycle clock, indicative; run alone after a first run under load was discarded):

| Engine | n | The candidate density per grid | Per candidate | The layer's step per grid |
|---|---|---|---|---|
| BAPS | 1,000 | 2.1 ms | 10 µs | 204 ms, the full step |
| BARS (Gaussian) | 1,000 | 3.3 ms | 17 µs | 196 ms, the full step |
| The hierarchical P-spline (10 curves of 100) | 1,000 | 1.9 ms | 10 µs | 3.6 ms, the split form's |
| BGP (unknown scale) | 400 | 20.1 s per call, whatever its cases | | 3.2 ms, the split form's |

  - The full step's 200 ms is its default step at S 1,000, as the margin's opt-in report measured it, the training side included.
  - BGP's seam is its refits: one O(n³) refit per accepted state, as its pointwise density pays, shared by every case in the call. The run already holds each draw's variance at its own test points, a cheaper path recorded as a follow-up in the brief.

- *Per commit,* the selected classes and the pinned tier: 2,056 of 2,057, the one skipped the GBP test below. *At the slice's end,* the whole assembly at d2e153b: 2,496 of 2,497, the same one skipped.

**Defects found.**
- *None in the seams.*
- *Handed on, not fixed here.* Under GBP mixing with a drift, `NormalVarianceMeanMixture`'s table reads ln F about 10⁻³ off near x = 0, while the density there is right to 12 digits.
  - At m/s = −0.2 and GBP(1, 1, 1.5), ln F(0) is −0.549285 against the integral's −0.550321. The table's slope runs 4% low across ±10⁻³·s, and the gap falls to 5·10⁻⁶ by 0.5·s. Without a drift F(0) is ½ exactly.
  - BARS's GBP noise reads its log CDF from that table, so its predictive distribution inherits the error near the median. Its sets, which read the density, do not.
  - The file belongs to the noise laws' owner. `TheLogCdf_UnderGbpMixingWithADrift_IsItsDensitysIntegral` is skipped until the fix, with the correct references as its assertions.

**Decisions taken beyond this brief.**
- *Which form takes each seam,* by two premises. The full step needs each draw's terms in sample, the units independent given the draw. Both forms need the candidate exchangeable with the units it is ranked against.
  - *The hierarchical P-spline serves the split form.* Its model takes one design for every replicate, so each training point's input is held by every curve and a new point's by none: no candidate point is exchangeable with the training points, and no partition restores it. Its latents match the training density's, so the layer cannot tell; the docstring and the page state it.
  - *BGP serves the split form.* f is integrated, never drawn, so the full step refuses its content by its own check (`LeavesUnitsOut`).
  - *The changepoint posterior has no seam.* Its posterior is over segmentations of the n positions it was given, one observation at each. A candidate at another position, a held-out one included, falls outside every segmentation, and a second value at a seen position is not a unit of that sequence. The refusal is in the class's docstring, since there is no method to refuse from.
- *BARS's support and CDFs.* A value outside the response's support takes −∞, as the seam allows. The Poisson model gives no CDF: a count's CDF steps, and its PIT is not uniform without a randomized tie. Whittle's candidate is read as an ordinate held out of the periodogram, under which its units are exchangeable.
- *The training units' log CDF* for BAPS and BARS, slice 4's shape for BASS. The hierarchical P-spline does not get one, since the full step's predictive distribution would read it and the shared design denies the full step.
- *`BayesianRegression.md` is not touched.* Another session holds uncommitted edits there. This page holds the API.

**Left unbuilt:** BGP's content with f drawn, a whole new curve for the hierarchical P-spline (E17), and the changepoint posterior's seam.

**Follow-ups, each with the decision it would change.**
- *BGP's seam at the run's own test points.*
  - The run computes each draw's conditional variance at its test points while the factor exists (`BgpDraw.ConditionalVariance`). A seam read there would take that variance and the mean instead of refitting: O(n²) per case instead of an O(n³) refit per accepted state, which the cost row shows is the seam's whole cost.
  - It decides whether the seam gains that path. If it does, a split form whose test points are the run's costs a fraction of a second. If not, the row stands.
- *BGP with f drawn.*
  - A content with f drawn jointly at the training and candidate inputs per retained state would put BGP's terms in sample, and the full step would take it, at O((n + T)³) per state.
  - It decides whether BGP gains a retained-parameters content. If not, BGP stays on the split form.
- *The seam stating the forms it serves.*
  - The hierarchical P-spline's seam passes the full step's checks though the guarantee does not reach it. A flag on `CandidateLogDensity` that the full step refuses would make that a refusal, not a docstring.
  - It decides whether the seam carries that flag. If not, the docstring stands.
- *Open, no decision waiting:* BASS's random-intercept fit. A candidate in a seen group would condition on the group's drawn effect, as the hierarchical P-spline's does on its curve. Whether BASS's groups share a design, which decides the full step's premise there, is E21's reading.

### 2026-11-08 — The margin made opt-in

Built at the user's decision of 2026-11-08, "Make it opt-in", relayed by the main session after it read 3fc80f9 and ran the Commit scope on HEAD (1,685 of 1,685).

**Landed** on `main`, nothing pushed (the commits are in the log beside this report):
- *`ConformalOptions.MarginBounds` is null by default.* Asked for, `Carried` or `Uniform` reads the margin as before. Not asked for, the scorer builds no orderings, runs no fourth product and counts no bands.
- *The readout says when it was not computed.* `ConformalMargin.Computed` is false; `Fragile`, `Count`, `Bounds` and the carried bounds are null; the plain set's `CoverageUnderError` is NaN. The inflated rule's readings are on every set.
- *Tests:* `TheMargin_IsOptIn_AndEveryOtherReadingStandsWithoutIt` (new). Slice 2's and the per-unit bounds' tests ask for the margin, and the coverage test checks that the default's decisions are the margin's. The fallback test checks that asking for the margin moves no fallback signal. The cost row times the default step and the margin's apart.
- *The page and this brief* record the decision with its date and words; E16's line too.

**Read before changing anything:** the fallback rule reads the inflated-set gap and k̂, and the split-half comparator reads two half-step plain sets. Neither reads the margin's band, so both stay as decided on 2026-11-07. Nothing serializes `ConformalOptions`, so no fingerprint moves.

**Gates** (`ConformalMarginTests`, `ConformalPerUnitBoundsTests`, `BassConformalMarginTests`; the readings with `BASS_QUALIFICATION=1`).
- *Both states,* over 60 seeded fixtures with weights, chains, far ranges and a candidate outside every draw's support, under both inflation bounds: 240 sets read without the margin and with it under each bound. Without it, the not-computed readout holds. With it, the margin is read. Every other reading is the same bit for bit: p-values, counts, scores, ESS, errors, inflated ranks, k̂, the gap, the split-half mass, the fallback signal and the widths.
- *The margin's readings, asked for,* are the per-unit bounds' report's, every printed line the same to its last digit as at 3fc80f9: on exact draws m̂ 5.04 and 2.57 carried at S 40 and 2,560, on BASS 97 and 96, the tail fixture 4.92 against 16.32, and the coverage test's readings 0.732 and 0.358.
- *Cost* (n 1,000, G 200, one grid, on a busier machine): the default step 170 to 193 ms at S 1,000 and 579 to 864 ms at S 4,000, against 728 to 803 and 1,550 to 3,222 ms with the margin.
- *Per commit,* the selected classes and the pinned tier: 673 of 673. *At the slice's end* (`scripts/test-select.ps1 -Scope Slice` on a snapshot of adb8344): the whole `Maths.Regression.Tests` assembly, 2,367 of 2,367, the eight pinned hashes among them, in 10 m 31 s. Nothing broke, so no follow-up commit was needed.

**Closed:** the per-unit bounds' follow-up, whether the margin stays always-on. **Open, no decision waiting:** a carried bound that leaves out the orderings whose weights degenerate, recorded in the per-unit bounds' report.

### 2026-11-08 — CB-3's per-unit bounds: the inflated rule per position, the margin on bounds carried over every ordering

Built at the user's order of 2026-11-08, relayed by the main session, from the Lean chip's per-unit CB-3 (80a1993).

**Found first: slice 2 had built most of the shape.** Its inflated rule read per-position bounds, and its margin had a carried option. Two of its bounds did not dominate what the Lean statements need:
- its carried bound read two orderings, the held one and the unit's own turn as the candidate;
- its uniform bound, the margin's default, read the held ordering alone, so it could sit below a unit's own-turn bound.

A probe settled the first before any build. On the conjugate fixture, some other unit's turn as the candidate gave a unit a larger error than slice 2's bound at 64 to 87% of (unit, grid value) pairs, by a median 1.3 to 2.9 times and up to 68. A unit far out, or of high leverage, makes the weights degenerate when it is the candidate, and every unit's error widens.

**Landed** on `main`, nothing pushed (the commits are in the log beside this report):
- *The margin on carried bounds, the default* (`MarginBounds.Carried`). Each unit's bound is the largest of its held bound and its bounds under the weights of every unit of the stratum as the candidate: three (n × S)·(S × n) products once per stratum. The candidate's is the larger of its own and its bounds under each unit's weights: one more (n × S)·(S × G) product per batch, from sums the step holds. `MarginBounds.Uniform` is the largest carried bound.
- *The inflated rule on per-position bounds, the default* (`InflationBounds.PerPosition`), with `InflationBounds.Uniform`, the largest of them, for comparison.
- *`ConformalMargin`* gains `CandidateCarriedErrors`, `CarriedUniformErrors`, `UnitCarriedErrors` and `Inflation`. `UniformErrors` keeps the held ordering's largest bound, now the inflated rule's uniform one.
- *Tests:* `ConformalPerUnitBoundsTests` (new); `ConformalMarginTests`' reference sweep and coverage test, for both forms; `BassConformalMarginTests`, reading both; the cost row, timing the orderings.
- *The page:* the margin's and the inflated rule's paragraphs, how each bound is estimated and why it dominates, and a section with the checks, slice 2's readings side by side, and the cost. The ledger gains the two uniform bounds and ε.

**How each bound is estimated, and why it dominates.** Every bound is z = Φ⁻¹(1 − ν/(2m)) of an estimate's Monte Carlo standard errors, under the weights its ordering puts on the draws.
- *Per position:* each estimate in the held ordering against its own error, so the union over the m positions keeps every one within its bound with chance 1 − ν, which `coverage_inflated` asks of an ordering.
- *Uniform, for the inflated rule:* the largest per-position bound, at least each by construction.
- *Carried:* a unit's largest bound over every ordering the draws represent, so it is at least its error in each of them, by the same union in each.
- *Uniform, for the margin:* the largest carried bound, at least each by construction.
- *Still indicative:* the standard errors are estimated, the normal approximation fails where weights degenerate (k̂ reads that), and the other orderings are the held draws reweighted, not draws of each ordering's own posterior.

**Gates** (`ConformalMarginTests`, `ConformalPerUnitBoundsTests`; the readings with `BASS_QUALIFICATION=1`). The full numbers are on [the page](../../docs/Maths/ConformalPrediction.md#the-per-unit-bounds-between-slices-4-and-5).
- *Against a reference:* 395 candidates over 80 seeded fixtures with weights, chains and far ranges. Every bound matches the errors read from public steps to 10⁻⁷ relative, and every band is the literal `FragileEstimated` count.
- *No set wider than the uniform rule's,* and no carried band above the uniform one, at every grid value of every fixture.
- *Coverage at 1 − α* (exact draws, S 40, k 2 of 20, se 0.0048): plain 0.905; inflated 0.935 per position and 0.945 uniform, against 0.89. The plain set covers above the margin's readings, 0.732 carried and 0.358 uniform.
- *A unit in a tail,* where the margin gains: the carried band reads 4.9 of 21 against the uniform 16.3 at S 40, and 3.4 against 8.1 at S 640. The plain set covers 0.950 and 0.943, above the carried readings of 0.790 and 0.843.
- *No unit standing apart* (a symmetric Beta–Bernoulli posterior): every carried bound equals the uniform one to 10⁻⁹, and the bands are equal.
- *Per commit,* the selected classes and the pinned tier: 136 of 136. *At the slice's end* (`scripts/test-select.ps1 -Scope Slice` on a snapshot of 477834c): the whole `Maths.Regression.Tests` assembly, 2,336 of 2,336, the eight pinned hashes among them, in 9 m 15 s. Nothing broke, so no follow-up commit was needed.

**Which readings change, and by how much.** The seeds are slice 2's.
- *Unchanged:* plain coverage, the per-position inflated set and the gap, exactly as slice 2 read them, on exact draws and on BASS.
- *The margin's default reading* moves from slice 2's uniform band to the carried one: at S 40 on exact draws, m̂ from 12.8 to 5.04 of 20, and the coverage it reads from 0.485 to 0.732.
- *The carried band, over every ordering:* on exact draws from 4.6 to 5.04 at S 40 and 2.5 to 2.57 at S 2,560. On BASS's median at the true values, from 72 to 97 of 100 (Gaussian) and 52 to 96 (t).
- *The uniform band, now the largest carried bound:* from 12.8 to 14.08 and 7.7 to 8.44 of 20; on BASS, from 98 to 100 and 94 to 98.
- *New, the uniform inflated rule:* on exact draws 0.4 to 1.4 points more coverage at 2.5 to 9% more width. On BASS it covers 0.986 and 0.945 at widths 3.44 and 1.83, against the per-position set's 0.958 and 0.937 at 2.07 and 1.49.

**What the readings show.**
- *The per-position bounds pay most on BASS.* The uniform bound there, a median 4.2 nats against each candidate's 0.08, widens the inflated set by two thirds.
- *The carried margin pays on exact draws, a third of the uniform band, and little on BASS.* Each unit's error over the orderings is wide there, since a unit far out degenerates the weights in its turn. Slice 2's carried reading on BASS, 72 and 52, came from bounds that did not dominate.
- *The tail unit carries the widest bound* in 37 and 32 of 40 data sets, 10 to 12 times the median unit's. Not in all: another unit can carry a wider one.

**Cost** (n 1,000, G 200, cycle clock, indicative). The units' orderings cost 0.33 to 0.45 s once per stratum at S 1,000 and 1.1 to 1.8 s at S 4,000, against 16 and 65 to 78 ms for the stratum's likelihoods. The fourth product brings the products to 70 and 261 ms per grid. On one grid the step costs 0.52 to 0.66 s and 1.7 to 2.4 s with the orderings built, which many cases share.

**Found in review and fixed.** The units' carried bounds were rebuilt into a new array for every case; they are now built once per stratum and shared.

**Follow-ups, each with the decision it would change.**
- *Whether the margin stays always-on* (slice 2's question, still open). It now costs the orderings once per stratum on top of slice 2's 2.6 to 3.5 times the plain step, and on BASS it reads 96 to 100 units of 100 under either bound. The inflated rule needs only the per-position bounds. Opt-in would skip the orderings, the fourth product and the bands.
- *A carried bound that leaves out the orderings whose weights degenerate,* adding their share of m to the miscoverage as ν is added. It would decide whether the margin can read anything on fits like BASS's, where a few units far out widen every bound.

### 2026-11-08 — Slice 4, the predictive distribution built: the split form's by default, the full step's as an opt-in

Built at the user's decisions of 2026-11-08, relayed by the main session: "Split-form CDF" for the default score, then "Yes, add it" for the full step as an opt-in beside it. They follow the reading of 2026-11-07, below.

**Landed** on `main`, nothing pushed (the commits are in the log beside this report):
- *`ConformalDistribution`,* per test case on its grid.
  - Q(y, τ) = (B + τ(T + 1))/m, with Q⁻ and Q⁺ (π under the CDF score) and the PIT at the true y. `coverage_tieBreak` and `rank_uniform` make the PIT uniform whatever the score.
  - The equi-tailed set takes k = ⌊α·m/2⌋ per tail, ties counted for the candidate in each, and covers 1 − 2k/m (CB-1 per tail).
  - Quantiles are read to the grid's step.
  - CRPS is τ-averaged, linear between grid values, its tails cut at the grid's ends, and within h(1 + 1/(3m)) of the exact value.
  - `Fall` and `ScoreFall` say how far a case's counts and its own score fall. A case whose counts fall keeps its set but has no distribution.
- *The default, `SplitConformal.Distribution`* (Vovk, Petej, Nouretdinov, Manokhin and Gammerman 2020). The score is the training fit's predictive CDF F̂(y | x), read by both tails from the seams' log CDF. Both of Vovk's conditions hold by construction, and CB-2 makes it valid for any fit.
- *The opt-in, `ConformalBayes.Distribution`.* The score is the in-sample predictive CDF F̃ under the add-one-in weights, two products a batch, the candidate's own score summed as the products' entries are so equal columns tie. Each case is read along its grid and flagged where its counts fall. The flag is a readout of its own, a mixing signal.
- *Shared pieces:* `CdfScore`, which compares on the tail that holds the digits; `CdfScores`; and `CdfScorer`. `SplitConformal.Density` and `Distribution` share one alignment check.
- *The page* gains the built section, with the decision beside the reading, its checks, the measured default and the cost row. The ledger gains Q beside Q(y, τ).

**Gates** (`ConformalPredictiveDistributionTests`; `BassConformalDistributionTests`, the measured default and the cost row with `BASS_QUALIFICATION=1`). The full numbers are on [the page](../../docs/Maths/ConformalPrediction.md#the-predictive-distribution-built-slice-4).
- *The counts against references* in seeded sweeps: every count equal, the split form's at 1,968 grid values and the full step's at 3,122, copied columns tying exactly.
- *The PIT is uniform on held-out data for any fit,* under the split form. On the conjugate model (m 20, 40,000 PITs a fit) the exact posterior, a wrong fit and a chain stopped early all read KS 0.005, with every bin within 0.004 of a tenth (se 0.0021). The fits' own PITs read KS 0.14 to 0.36. On BASS (m 200, 2,000 PITs a fit) the healthy, t, wrong and stopped-early fits are all within four se in every bin.
- *Vovk's conditions hold by construction* on every grid case of every fit, conjugate and BASS: no fall of the score or of Q, and Q⁻ = 0 below the grid, Q⁺ = 1 above.
- *The opt-in's PIT, over every case and the unflagged:*
  - uniform on the conjugate exact posterior and on BASS's healthy, t and wrong fits, the last a posterior that mixed of a wrong model;
  - off where the draws are no posterior: the conjugate wrong fit (a bin 4.1 se out), two chains stuck apart (KS 0.17 over every case) and BASS's chain stopped early (a bin 3.8 se out);
  - where the draws are nearly a point mass, as with the conjugate chain stopped early, every unit is read alike and it holds.
- *The flag.* The conjugate exact posterior flags no case. Two chains stuck apart flag 287 of 500, 267 of them with the candidate's own score falling; each keeps its set, with NaN quantiles and CRPS.
- *Agreement with slice 3's sets.* The equi-tailed set is slice 3's split rule in each tail at α/2 on the CDF score at all 61,127 grid values in the bulk. On a fit that is right (m 2,000), its ends, the density set's and CQR's meet the fit's own central interval, each mean gap within four se and a grid step.
- *CRPS against a quadrature reference:* within 1.3·10⁻³ at m 20 and 4.6·10⁻⁴ at m 200, against the bound 0.010.
- *Per commit,* the selected classes and the pinned tier: 129 of 129. *At the slice's end* (`scripts/test-select.ps1 -Scope Slice` on a snapshot of 77d0c24): the whole `Maths.Regression.Tests` assembly, 2,297 of 2,297, the eight pinned hashes among them, in 12 m 1 s. Nothing broke, so no follow-up commit was needed.

**The measured default.** Both forms ran on the same data, the full step fitting all of it.
- *Where the draws are a posterior that mixed,* the full step is sharper by more where data are scarce. Its CRPS is 10% lower at n 49 on the conjugate model, and 2.5 to 3.7% lower at n 298 on BASS, at the split form's width. It flags 0 to 0.1% of cases there.
- *Where they are not,* it breaks: two chains stuck apart flag 66% of cases.
- *The default stays the split form, as decided.* The page writes down the regime to opt in: draws a posterior that mixed, and data too scarce to spare a calibration split.
- *Cost.* The layer's own cost per case is 1.4 against 4.8 ms on the conjugate exact posterior at a grid of 401, and 12 against 72 ms on the healthy BASS fit at a grid of 901 and S 1,000, read on the final code (below). The split form's data cost is 19 of 49 units and 199 of 298.

**Cost row** (cycle clock, indicative), every row in one run on the final code. The density score, read in the same run, reads as in slice 3, so the machine was about as loaded as then. Per grid of 200 at n 1,000:
- *the split form:* 2.3 ms at S 1,000 and 9.1 ms at S 4,000, after its calibration PITs once per stratum (18 and 70 ms), against the density score's 2.1 and 6.7 ms;
- *the full step's opt-in:* 52 and 219 ms, against the density full step's 152 to 194 and 538 to 599 ms.

**Found in review and fixed.**
- *The split form ran ten times the density score's cost.* It read 1 − F by expm1 below the median, where 1 − F ≥ ½ loses nothing, and walked the matrices down their columns. It now takes one exp a value and a log only above the median, and walks rows in the matrices' order with each column's sums unchanged. The per-grid cost went from about 20 ms to 7.0 ms at S 1,000.
  - The BASS gates and the measured default first ran on the build before this change, after the tie fix below. The change moves the upper sums by an ulp where ℓ ≤ −ln 2 and keeps every sum's order, so it could not move their readings beyond rounding. It did move their cost per case, read on the wall clock.
  - All five runs were read again on slice 4's final code (b48158b), at the user's go-ahead later on 2026-11-08. Every reading came out the same to its last printed digit. The cost per case fell for both forms, the split form's most: the healthy BASS fit reads 12 against 72 ms, where it read 80 against 308. Part of the fall is a lighter load on the machine. The page and this report carry the new costs, and the cost row was read again with them, after the full step's setup was tiled.
- *The full step's own score was summed in another order than the products,* so a copied column did not tie with its unit. It now sums by fused multiply-adds in draw order, and the sweep's 296 near-ties became exact ones.
- *The flag test asserted that a flagged case's own score falls.* π can fall while it rises, when a unit's score overtakes it: 20 of the 287 flagged cases do. `Fall` and `ScoreFall` are now read apart.
- *Three ungated tests ran for one to two minutes.* They were shrunk to 13 to 49 s, each keeping its gate.

**Notation:** Q, the split form's threshold, beside Q(y, τ), the distribution; F̂, the split form's score, beside F̃. Both are in the ledger.

**Not built,** recorded on the page:
- the split form's set ends and quantiles in closed form, F̂⁻¹ at the calibration PITs' order statistics, which needs each engine's quantile function;
- cross-conformal predictive distributions;
- a CB-3 margin for the CDF score.

**Follow-ups:** none with a decision waiting. Slice 5, the other engines' candidate seams, waits for the user.

### 2026-11-07 — Slice 4, the predictive distribution: the reading, and the score not adopted

Started at the user's go-ahead of 2026-11-07 ("Yes, start both", relayed by the main session). The row reads the in-sample predictive CDF's monotonicity first and adopts it as the score only if it holds. It does not hold in general, so the slice stops at the reading, as the go-ahead asked, and the score waits for the user's choice.

**The reading.** The full numbers are on [the page](../../docs/Maths/ConformalPrediction.md#the-predictive-distribution-slice-4-the-reading). The score is F̃(y) = Σ_s w̃_s(y) F(y | x, θ_s) under the add-one-in weights, swept along a grid and compared on the tail that holds its digits.
- *It holds where the draws sit in one mode.* That covers the conjugate model's exact posterior at S 50, 400 and 4,000 and a chain stopped early there, and BASS's healthy, Student t and wrong fits. The t fit's π fell once, by 1 of 100.
- *It falls where the draws' predictives at x separate.* With two conjugate chains stuck apart, the score falls on 44 of 100 cases and π by up to 4 of 20. With BASS's chain stopped early (one chain, no burn-in, 20 draws), it falls on 3 of 50 cases and π by up to 3 of 100.
- *Why.* The slope is the reweighted density plus Cov_w̃(∂_y ln p, F), which is negative where draws differ in location. Two equal draws N(±a, 1) give φ(a) − a(Φ(a) − ½) at their midpoint, which falls once they stand 1.75σ apart (a ≈ 0.877).
- *Vovk's second condition fails almost everywhere.* Half the conjugate units and 89 to 100% of BASS's rise with the candidate. Where π stays monotone it is a reading, not a theorem.
- *What a falling score loses.* Its PIT would still be uniform, since a conformal p-value needs no monotone score. What it loses is the CDF: quantiles, CRPS and a distribution to report.

**The alternatives,** recorded on the page for the user's choice:
1. *The split form's predictive distribution* (Vovk, Petej, Nouretdinov, Manokhin and Gammerman 2020). Its score is the training fit's predictive CDF F̂(y | x). Both conditions hold by construction, CB-2 makes it valid for any fit, and it needs no product over the units. This is the recommendation: slice 3's split and these seams are all it needs.
2. *The full step's score with a per-case check.* It is a CDF where π rises and flagged where it falls. It fails where the sampler is in doubt.
3. *A monotone envelope.* Its PIT is not uniform; not recommended.
4. *Cross-conformal predictive distributions.* Approximate validity, at K fits.

**Landed** on `main`, nothing pushed (the commits are in the log beside this report):
- *The log CDF seams.* `PointwiseLogDensity.LogCdf` (draws × units, each unit one observation) and `CandidateLogDensity.LogCdf` carry each draw's ln F beside its ln p. A probability rounds F to 1 in the upper tail. Held that way in a first reading, the units rising toward the candidate tied with it, and π fell spuriously on 4 to 54 of 100 conjugate cases. The log keeps 1 − F, and those falls vanish.
- *BASS fills them on `cdf: true`* (`BassFit.PointwiseLogDensity`, `BassFit.CandidateLogDensity`; `BassCandidates.UnitLogCdf`), from the law at the draw (`IBassNoiseLaw.At`) or N(f, σ²). Above each draw's median it reads ln(1 − S) from the law's survival (`BassCandidates.LogCdf`). No noise law changed.
- *The reading as tests,* the quick form ungated and the full form with `BASS_QUALIFICATION=1`, and the page's section.
- *Outside the listed files:* `BassFit.cs`, the seam's two entry points that slice 1 added. Its diff was clean before the edit.

**Gates.**
- *Units and candidates read one distribution* (CB-4). At the training cases the units' log CDF equals the candidates', bit for bit, under all six laws.
- *Each draw's CDF is its density's integral* by Simpson's rule, step 0.0005. The largest gaps are 1.5·10⁻¹⁴ (Gaussian), 2.2·10⁻¹³ (t), 7.2·10⁻⁷ (asymmetric Laplace, at its kink), and 2.3·10⁻⁷ to 3.2·10⁻⁶ for the tabled laws.
- *The upper tail* holds the law's survival to 10⁻¹², out to 32 scales.
- *The refusals.* A unit of two observations, the collapsed content, a bad shape, a value above 0 or NaN, and a weighted fit's units are refused.
- *The reading's pieces:* the tail comparison, and the two-draw slope against its closed form.
- *Nothing upstream moves.* All 235 neighbouring tests pass: the conformal layer, the pointwise and candidate seams, PSIS-LOO and the predictive checks. No sampler or pointwise matrix changed, so no pinned hash can move.
- *The slice-end run* (the rule of 2026-11-07, `scripts/test-select.ps1 -Scope Slice` on a snapshot of 298e569): the whole `Maths.Regression.Tests` assembly, 2,222 of 2,222, the eight pinned hashes among them, in 8 m 29 s. Nothing broke, so no follow-up commit was needed.
- *Not run.* PIT on held-out data and the cost row wait for an adopted score. The seam's own cost is the law's CDF read beside its density read: a closed form for the Gaussian, t and asymmetric Laplace, and for the tabled laws one table built per draw.

**Found and handed on, not fixed here.** The Normal-Wald's, skewed t's and generalized hyperbolic's `LogCdf` comes from `NormalVarianceMeanMixture`, whose tail table ends at the 10⁻¹² tails. Beyond the upper end (`TailTable.Read`, past `_at[^1]`), it integrates the lower side, so it loses 1 − F's digits and can read above 0: +2.9·10⁻¹⁵ for the generalized hyperbolic at 32 scales, and above 0 too on a Normal-Wald fit's skewed draws. Reading that end as ln(1 − S) from the upper integral would keep them. The file belongs to the noise laws' chip. The seam reads the survival above the median, so it does not wait on the fix.

**Found in review and fixed.**
- Slice 3's report and the page said no noise law reports its variance. `IBassNoiseLaw.Variance` exists, but it is infinite for t at ν ≤ 2, and prediction gives none per case. Both lines are corrected, and the IQR stays the spread that serves every law.
- The page wrote the CDF scores as α, which on that page is the miscoverage level. They are F̃ now.
- The tabled laws' integral tolerance was tightened from 2·10⁻⁴ to their tables' 2·10⁻⁵.

**Notation:** the ledger now holds slice 3's ρ(x) beside σ and σᵢ, and slice 4's F̃ beside α and Vovk's α.

**Follow-ups:** the score for the predictive distribution, the user's choice among the alternatives above. The mixtures' far upper tail is for the noise laws' chip.

### 2026-11-07 — Slice 3, the split form

Built at the user's go-ahead of 2026-11-07 ("Slice 1 + conformal 3", relayed by the main session).

**Landed** on `main`, nothing pushed (the commits are in the log beside this report):
- *`SplitConformal`,* engine-free, gives CB-2 (`split_coverage`) for any fit. There are three scores, each ranked within the test case's stratum among the calibration units by the full step's rule, ties counted for the candidate:
  - *the predictive density* (`Density`), ln Σ W_s p(y | x, θ_s), read through the candidate seam at each calibration case's response. Its arithmetic is the full step's predictive, bit for bit. Its sets are `ConformalSet`s on the grid under `ConformalForm.Split`, with no margin, ESS or diagnostics;
  - *CQR's* (`Quantiles`), min(y − q_lo, q_hi − y), Lean's `cqrScore`, whose interval [q_lo − Q, q_hi + Q] comes in closed form;
  - *the standardized residual's* (`Residuals`), −|y − μ|/ρ, whose interval is μ ± Qρ.
- *`SplitCalibration`* serves any score with per-stratum thresholds, ranks, π and the accept rule.
- *Upgrades:*
  - a `ConformalSet` with no margin reads its inflated rule as its plain one, and its coverage under error as its coverage;
  - `ConformalPrediction` carries its `Form`;
  - `PredictiveChecks.Rows` reads a weights list, so the candidate density's weights take the same rows.
- *No engine code changed.* BASS joins through the seam it already has (slice 1), with `BassPrediction.Quantile` for CQR and `ObservationMean` with the predictive's interquartile range for the residual. A law's variance is infinite for t at ν ≤ 2 and prediction reports none per case, so the IQR is the spread that serves every law (corrected in slice 4's review: this line first said no law reports its variance).
- *The page* gains a split form section with its choices and what is not built, its checks and cost, and the API rows.

**Gates** (`SplitConformalTests`; `BassSplitConformalTests` and the cost row with `BASS_QUALIFICATION=1`). The full numbers are on [the page](../../docs/Maths/ConformalPrediction.md#checks-and-cost-slice-3).
- *The rule against the count* (`card_rankLe_le`): 2,000 seeded calibrations with ties, infinite scores and empty strata, every probe agreeing.
- *CB-2 by exact enumeration* (`split_coverage`, `exchangeable_given_training`, `coverage_orderings_eq`) over 200 laws with random fits, wrong ones among them. Coverage is at least 1 − k/m in every law, and exactly 1 − k/m for injective fits on distinct values.
- *The density score* matches a reference in logs to 10⁻¹² and the full step's predictive bit for bit. A copied column ties.
- *CQR's set is an interval* (`cqr_ordConnected`), and so is the residual's. Over 7,980 sets, 550 of them empty, the rank rule runs unbroken and agrees with the closed form except within rounding of an end.
- *Coverage for any fit.* On the conjugate model (m 20, nominal 0.900, se 0.0021), the exact posterior, a wrong fit and a chain stopped early cover 0.896 to 0.900 under all three scores. The two bad fits' own intervals cover 0.24. On BASS (m 200, se 0.0067), a chain stopped early, a fit that never sees x₀, and a healthy fit cover 0.907 to 0.912. The four BASS facts share their data, so their readings move together. CQR on a quantile pair covers 0.906 where its own intervals cover 0.853.
- *Nothing upstream moves:* 100 of 100 neighbouring tests pass, PSIS-LOO, the pointwise density and slices 0 to 2 among them. No engine file changed, so no pinned hash can move.

**Cost** (cycle clock, indicative, loaded machine), per grid of 200 at n 1,000 calibration units:
- the density score costs 2.1 ms at S 1,000 and 7.3 ms at S 4,000, against the full step's 644 to 700 ms, plus the calibration's predictive once per stratum (12 and 49 ms);
- CQR's and the residual's calibration with 200 intervals each cost under 1 ms.
- *A cost finding, not a defect.* On BASS the readings cost more than the fits. Each Gaussian fact at S 1,000 read four predictive quantiles at 299 cases and took 812 to 860 s for 30 replications. The quantile pair, read from surface means, took 27 s for 60 fits. If the split form on BASS runs routinely at large S, `PredictiveMixture.Quantile`'s speed sets its cost. Nothing waits on it now.

**Found in review and fixed.**
- The coverage test's width grid was too narrow: 7 of 1,800 density sets reached its ends, so their widths were read short. It now spans [−50, 50], none reach it, and the chain stopped early reads 17.62 where it read 15.36.
- Two `ConformalSet` parameter docs described only the full step; they now give the split form's readings.
- E16's slice 2 line still said the fallback signal fires; it now records the decision of 2026-11-07.

**Notation:** the residual's scale is ρ(x), Lei et al.'s symbol, since σᵢ on the page are the full step's scores. This is for `issues/notation-ledger.md` once the GBASS chip commits it.

**Not built,** recorded on the page with what each would change: CQR's asymmetric form, cross-conformal and the jackknife+, and the tournament correction (E16).

**Follow-ups:** none with a decision waiting. Slice 4, the predictive distribution, waits for the user.

### 2026-11-07 — The fallback diagnostics on every set, with no trigger

This executes your decision of 2026-11-07, recorded under "Decisions" above with its reasons and the open reading.

**Landed** on `main`, nothing pushed (the commits are in the log beside this report):
- *Two diagnostics on every set.* `ConformalSet.Gap`, the inflated-set gap, and `SplitHalfMass` with `SplitHalfDisagreements`, the main session's split-half construction.
  - The split-half is built into the layer. Two plain-mode scorers take the draws' halves (each chain's first and second halves, or alternate rows without chains) and run the score product and the ranks only.
  - Its mass is read under the full draws' predictive, which the full step already holds.
- *No trigger.* `FallbackShare` is now null by default and `FallBack` stays false. With a share stated, the rule fires as before (gap or untrusted k̂). The rule was kept, not removed.
- *The page and the brief* record the decision, why, and the open reading of why the halves disagree on BASS while coverage holds.

**Gates.**
- *The layer's split-half mass is the reference's* (two public steps on the half densities) to rounding, over 40 cases with and without chains, 44 of them with halves that disagree.
- *The edge cases are NaN:* a grid of one value, and a half that gives some unit zero density under every draw. The step goes on.
- *The rule is off by default and the rule when asked,* at shares 0 to 1, on the strained, slow-chain and healthy cases.
- *Nothing upstream moves:* 63 of 63 neighbouring tests pass, among them the eight pins and PSIS-LOO. The comparison test now reads the layer's own diagnostic.

**Defect found and fixed in review.** A half in which some unit has zero density under every draw would have made the half scorer throw, failing a step the full draws can take. The diagnostic is now NaN there.

**Cost** (indicative, the machine loaded, about 40% above slice 2's readings for unchanged parts), per grid of 200 at n 1,000:
- the split-half costs 21 to 25 ms at S 1,000 and 82 to 97 ms at S 4,000, plus its halves' likelihoods once per stratum (15 to 65 ms);
- the gap costs almost nothing beside the inflated ranks;
- the full step is 186 to 203 ms at S 1,000 against slice 0's plain 36.7, and 644 to 700 ms at S 4,000 against 149.8.

**Follow-ups, each with the decision it would change.**
- *Whether the margin stays always-on* (slice 2's question, still open). The gap needs only the per-position errors, the two error products and r_eff. Turning the margin's bands and k̂ off would save 60 to 160 ms per grid in these rows and keep both diagnostics.
- *A trigger,* once a fixture with chains like BASS's shows a measured shortfall to set its threshold by. None is scheduled.

### 2026-11-07 — The fallback signal: the inflated-set rule and a split-half comparator, side by side

At your request through the main session, a split-half comparator was measured beside the layer's signal on the same fixtures, with the same discrimination test. It is the main session's construction, not from a paper, and the page says so.
- *Neither was adopted.* The layer's `FallBack` still reads the inflated set, as you chose directly, and the comparator lives only in the test.
- *How the comparator reads.* The retained draws are split into two halves: each chain's first and second halves, or alternate draws without chains. It builds the plain set from each half on the same grid, then reads the mass, under the full draws' predictive at x, of the region where the two half-sets disagree. It fires at a mass above 0.1α.

**Landed** on `main`, nothing pushed: `ConformalFallbackComparisonTests` (gated, `BASS_QUALIFICATION=1`), the side-by-side section on `ConformalPrediction.md`, and this report. Both rules were read case by case: the signal on a grid at x, and the plain set's decision at the true value at the same x.

| Fixture (1,000 cases each) | Plain coverage | The layer's signal: fires, covered on / off | Split-half above 0.1α: fires, covered on / off |
|---|---|---|---|
| Exact draws, S 40 | 0.907 | 59%, 0.912 / 0.899 | 20%, 0.827 / 0.927 |
| Exact draws, S 160 | 0.904 | 8%, 0.924 / 0.902 | 10%, 0.800 / 0.916 |
| Exact draws, S 640 | 0.888 | 2%, 0.833 / 0.889 | 6%, 0.738 / 0.898 |
| Exact draws, S 2,560 | 0.905 | 1%, 1.000 / 0.904 | 3%, 0.923 / 0.905 |
| Two chains, ρ 0.9 | 0.888 | 5%, 0.902 / 0.887 | 16%, 0.811 / 0.903 |
| Two chains, ρ 0.99 | 0.886 | 15%, 0.759 / 0.908 | 34%, 0.794 / 0.933 |
| Two chains, ρ 0.995 | 0.889 | 13%, 0.858 / 0.894 | 34%, 0.854 / 0.907 |
| Two chains, ρ 0.999 | 0.887 | 5%, 0.635 / 0.901 | 24%, 0.753 / 0.929 |
| BASS Gaussian, S 1,000 | 0.905 | 54%, 0.907 / 0.903 | 74%, 0.911 / 0.889 |
| BASS t, S 1,000 | 0.894 | 32%, 0.896 / 0.893 | 61%, 0.893 / 0.895 |

**What the measurement says.**
- *Where the plain set falls short, the split-half mass finds it.* On exact draws and slowly mixing chains, the cases it flags cover 0.75 to 0.85, and the rest cover 0.90 to 0.93. It is the rule that fires where the plain set falls short and stays quiet where it holds.
- *The layer's signal discriminates only for the slowest chains,* and fires rarely there. At small S its trust trigger (k̂ undiagnosable or above PSIS's threshold) fires where the plain set holds: at S 40, 59% of cases, covering 0.912.
- *On the BASS fits neither rule discriminates.* The plain set loses nothing there, yet the split-half mass is large (median 0.033 Gaussian, 0.014 t) and fires on 74% and 61% of cases. So 0.1α, the stated starting point, is not a calibrated share for these chains.
- *At slice 1's degenerate Gaussian weights* all eight true values are misses: outlying truths the conformal rule refuses, part of the 10% by design. The layer's signal fires at six of their grids, the split-half at seven.
- *The slow-chain fixture is the one this measurement added.* It is the first where the plain set's coverage varies with how well the draws represent the posterior, which the discrimination test needs.

**Cost.** The comparator is two half-size full steps per call, about one extra score product, measured here through the public step. Built into the layer, it would need only the two half-size score products, not their margins.

**Follow-ups, each with the decision it would change.**
- *Which rule is the signal.* Each outcome changes this:
  - The split-half mass: the signal would discriminate where the error costs coverage, but its share needs calibrating on chains like BASS's, where 0.1α fires on most cases that cover.
  - The inflated-set gap, as now: it fires less often, but at small S on cases that cover, through the trust trigger.
  - Both, the split-half mass firing and the gap reported: the signal reads typical error and the inflated set stays the width-for-safety option.
- *The split-half share,* if adopted, set on a fixture with chains like BASS's against measured shortfall. The BASS fits here showed no shortfall to set it against.

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
