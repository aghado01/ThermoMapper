# The conformal layer: E16's build (chip brief)

Date: 2026-11-06. Status: slices 0, 1 and 2 built 2026-11-06, the fallback signal read from the inflated set on 2026-11-07 at the user's call, with a split-half comparator measured beside it; on 2026-11-07 both became diagnostics on every set, with no trigger, at the user's decision; slice 3, the split form, built 2026-11-07 at the user's go-ahead (reports below); slice 4 read on 2026-11-07 (the in-sample predictive CDF is not monotone in y where the draws' predictives separate) and built on 2026-11-08 at the user's decision: the split form's predictive distribution by default, the full step's as an opt-in with its flag; CB-3's per-unit bounds built on 2026-11-08 at the user's order, a slice between 4 and 5, and the margin made opt-in the same day at the user's decision; slice 5 waits for the user. Written at the user's request from [E16](bass/extensions-brief.md#conformal-bayes-the-consumer) and the [conformal Lean contract](../../lean/Protolemmata/Conformal/conformal-lemmas.md).

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
| 6. Conformal test martingales | Exchangeability checked online, as [E35](bass/extensions-brief.md#conformal-test-martingales-exchangeability-checked-online-e35) sets out: each new score ranked against every earlier one with the randomized tie, a betting function, and an alarm with its guarantee. Ordered by the user on 2026-11-08 to follow slice 5 | Slice 3's fixed-fit scores; slice 4's log CDF and its rank code | E35's, set when its brief is written |

## Gates common to every slice

- **Nothing upstream moves.** The layer reads; every pinned hash and fingerprint holds: the pinned tier at each commit, the whole assembly at the slice's end ([targeted runs](../../CONTRIBUTING.md#testing-standards)).
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
- Heavy runs at below-normal priority with capped threads; the machine is shared. Tests run by [targeted runs](../../CONTRIBUTING.md#testing-standards) (`scripts/test-select.ps1`): the selected classes per edit, with the pinned tier per commit, the whole assembly once at the slice's end, and the qualification gates as each slice states.
- New beside old; a removal is its own question to the user, with what is lost.
- Findings framed by problem class, not by any one application.
- A follow-up goes in a report with the decision it would change and what each outcome would do. A question with no decision waiting on it is recorded as open, not proposed.

## Report

Under "Reports" below, newest first, per slice: what landed by commit; the gates' numbers; coverage and width by fit; the cost row; defects found; decisions taken beyond this brief; what was left unbuilt and why. E16 takes a line per slice that lands.

## Reports

Newest first.

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
