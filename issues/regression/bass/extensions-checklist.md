# BASS extensions: checklist

Companion to [extensions-brief.md](extensions-brief.md). An item is ticked when its report lands. Reports go under "Reports", newest first, each with its date, commits, evidence and remaining limitations.

## Checklist

### Entries

- [ ] **E1** Generalized-hyperbolic noise (GBMARS), with OQ20's BARS side. Building ([brief](gbass-noise-brief.md)): slices 0–4 and the GIG's ratio-of-uniforms draw landed by 2026-11-06; slices 5 and 6, the engine's reach and OQ20's BARS side, landed 2026-11-07; slice 7 building.
- [ ] **E2** Robust centring and bases (OQ16).
- [ ] **E3** Free-knot tensor-product surfaces (ABMRS).
- [ ] **E4** Additional proposal kernels. The Bernoulli-inclusion birth is briefed with KHAOS ([brief](khaos-brief.md), 2026-11-07).
- [ ] **E5** KHAOS's modified g-prior. Briefed with KHAOS ([brief](khaos-brief.md), 2026-11-07).
- [ ] **E6** Basis-estimation uncertainty.
- [ ] **E7** Moments for piecewise polynomials.
- [ ] **E8** Active-subspace analysis.
- [ ] **E9** Cross-engine evidence.
- [ ] **E10** Direct-functional BASS for families of SPC sweeps.
- [ ] **E11** Elastic emulation of misaligned functional outputs.
- [ ] **E12** Residuals correlated along the output index.
- [ ] **E13** The free-energy curve over the ladder.
- [ ] **E14** Semi-modular calibration.
- [ ] **E15** Calibrated intervals for loss-based posteriors (BARS's Gibbs-loss route first).
- [ ] **E16** Conformal Bayes, the consumer, with conformal predictive distributions. Building as an engine-wide layer ([brief](../conformal-layer-brief.md)): slices 0–2 landed 2026-11-06.
- [ ] **E17** Conformal bands for whole curves.
- [ ] **E18** Posteriors for readouts by resampling units.
- [ ] **E19** Polynomial-chaos terms. Briefed ([brief](khaos-brief.md), 2026-11-07), with CRAN's `khaos` 2.1.1 pinned.
- [ ] **E20** Data-augmented observation models.
- [ ] **E21** A prior-precision builder, with ridge blocks. Briefed with KHAOS ([brief](khaos-brief.md), 2026-11-07).
- [ ] **E22** Correlated outputs and coupling within a run.
- [ ] **E23** Multi-fidelity emulation through the offset.
- [ ] **E24** Guarded estimators.
- [ ] **E25** Sensitivity under dependent inputs.
- [ ] **E26** Stochastic emulation through latent inputs.
- [ ] **E27** External benchmarking.
- [ ] **E28** Bases that mix term types (lead).
- [ ] **E29** Categorical levels: learned positions and a similarity readout (lead).
- [ ] **E30** Choosing the next runs.
- [ ] **E31** Tree ensembles, with BCF's multiplier.
- [ ] **E32** Additive models with a GP component.
- [ ] **E33** Regression-family integration.

## Reports

### E1's build brief, and CRAN's GBASS pinned — 2026-11-04

The user's reading: GBASS is BASS with a different noise model, a diff and not a new engine. [The build brief](gbass-noise-brief.md) states it: four new distributions and functions, the mixing law, one observation model and the upgrades, in slices from Student t on. CRAN's GBASS 2.0.1 is in `r/renv.lock` and ran for the first time; its code matches the clone E1 read but for four small differences, in the brief. Order by the user: GBMARS, then KHAOS. Nothing in the engine is built, and E1's checkbox stays open.

### GBASS read — 2026-11-04

GBMARS's reference implementation, read against the paper, E1 and G1; nothing built. The detail is in [E1](extensions-brief.md#generalized-hyperbolic-noise-gbmars).
- *E1 gains* the sampled γ and β it had not named, the MHN as an exact draw at variance ∝ w in the v parameterization, and five defects to correct. The text's two defects are slips in the text alone.
- *E20 gains* logistic noise by Stefanski's mixing, in GBASS since July 2026 and not yet in its sampler. Two reading leads go to the [lineage review](lineage-review.md#reading-leads).
- *G1 gains* candidate reading (h), a whole-term redraw.
- *Not checked:* the package was not run; `rMHN`'s series branch; `lag_beta`'s reason; the GBP step's scale, read and not measured.

### The tiny-shape evidence fixed — 2026-11-03

The defect [the priors' check](#the-moving-response-priors-checked--2026-11-02) found, on the user's go-ahead. No default changes and the hashes hold. The forms and the alternative are in [the docs](../../../docs/Maths/BayesianRegression.md#bass-priors-that-stop-reading-a-moving-response).
- *Built.*
  - `Gamma.SampleLog`: ln X on `Sample`'s stream, which stands where X would underflow.
  - `ResolvedSigma2Rate` holds to v: `LogDensity` at any v, `SampleLogRate` through `SampleLog`, `PriorAt` and `RateAt` at the smallest normal double r_min or above, and `LogRemainder` for the rest of rᵃ, which `BassEvaluation` adds.
  - `BassState.StartLogRate`: a chain starts no lower than r_min, while a prior rung and SMC's particles keep the draw.
  - `Gamma.LogCdf` and `LogSurvival` where Rate·x falls under 10⁻³⁰⁰, the defect noted beside this one.
  - The two U rows the [inventory](../capability-inventory.md) added, in [the brief's table](extensions-brief.md#upgrades-to-existing-code).
- *Evidence.* All on the final tree at below-normal priority: Distributions 287/287, Regression 1,821/1,821 with `BassPathProfile` left out and the bit-identity tests in, and the oracle's BASS and quantile filters 112/112. That also reruns `ab6b45b`'s tests, whose docs I read and whose code I did not. ln Z, reported error in brackets, against each enumeration with r integrated whole:

  | a, h | Enumeration | Ladder | Linear ladder | SMC |
  |---|---|---|---|---|
  | ½, ½ | −16.294 | −16.296 | −16.281 | −16.302 |
  | ½, 10⁻³ | −21.000 | −21.130 (0.211) | −21.019 (0.069) | −21.010 (0.044) |
  | ½, 10⁻⁴ | −23.294 | −23.248 (0.131), 30,000 samples | not run | −23.471 (0.140) |
  | 10⁻³, 10⁻³ | −21.686 | −21.668 (0.034) | not run | −21.685 (0.009) |

  The last row is where r's own posterior reaches below the doubles, so it is the check on `LogRemainder`.
- *The review found.*
  - **The fix as first proposed was not enough.** With a under about 0.01 it would have left a factor on the slice below the doubles, hence `LogRemainder`. A draw's rate could have read 0, which the pointwise densities refuse, hence `RateAt`. And a chain from a prior draw far below the doubles could not have climbed within burn-in, hence `StartLogRate`.
  - **A tiny h costs precision, and that is no new bias.** The first run at h = 10⁻⁴ read −23.75 on a ladder and −23.47 on SMC, low by 0.45 and 0.18, since of order h/a of the prior's draws land where the likelihood is, and cutting the prior had hidden that by a factor of 14. A budget probe on the same seeds: the ladder at 30,000 samples read −23.25 (0.13), and at 16 levels −23.19 (0.13); SMC at 16,000 particles read −23.30 (0.07). The reported error covers each gap, and a ladder's largest pair's rejection read 0.998 and above, so `Saturated`. So the check allows three reported errors up to a ceiling of 0.4 on the error, in place of a flat 0.15.
  - Two of the chip's tests asserted the old behaviour, `LogDensity(−800) = −∞` and a state below the doubles off the support, and now assert the new.
- *Limits.*
  - The integrated form's slope at b = 0 is unbounded for h ≤ 1, since E[1/r] diverges. That form reports no path sampling, so nothing reads it.
  - `LogRemainder` needs S far above r_min, so a response whose squares sit near 10⁻³⁰⁸ would read wrongly, as it would elsewhere.
  - A draw retains r_min for an r below it, since a double cannot hold less. ln r stays on the state.
  - The alternative's size, about 0.25 nats at a = h = 10⁻³, is derived, not run.

### A design whose working response moves — 2026-11-02

`ab6b45b`. A chain takes a new working response between sweeps: the design refits the current structure, and the state's fit, the target's statistics, pointwise and the training fit follow. No default changes. The forms, options and trade-offs are in [the docs](../../../docs/Maths/BayesianRegression.md#bass-a-design-whose-working-response-moves).
- *Built.*
  - `WorkingResponse.Moved(offset, weights)`, the model's reading at another latent state.
  - `BassStructureMoves.ProposeFollow`, a candidate for a block to weigh, and `Follow` and `BassChain.Follow`, which move the chain. A block reaches them by `BassStructureMoves.Of(chain)`.
  - An offset alone keeps the columns, factor and gate, about 10–16 µs per move on Friedman. Weights whiten and factor again, about 80–140 µs, by a built design's arithmetic to the bit.
  - `BassDraw.Working`, read by pointwise, the conditional means and the residual check. `BassTrainingFit.Offset`, with `Residuals` on y − offset − f̄.
- *One choice beyond the brief.* The statistics went on the state, not on a target per snapshot: `BassStructureFit.Working`, read by the one target. σ²'s rate went the same way, for the same reason: a swap carries an evaluation, an SMC stage holds one target, and the flat histogram wraps one.
- *For the user's review: F under a moving W.* Built with the gate reading W: a move that leaves F is refused, and a block keeps the old state, which is exact. **BARS's review rejected this reading:** its robust block admits on G alone, since a joint gate makes the support depend on latent outlier weights. BARS could, because it factors P_a. BASS's design holds a factor of G, so the other reading needs a path with no such factor, on the pencil alone, and `Weighted` cannot take it at all. It adds to this one with no seam changed.
  - *Measured:* 0 of 7,200 draws refused on Friedman and direct functional, over t(3), t(1) and a Laplace's scales, with the smallest angle at 10⁻³, nearly seven orders above s_min. At N = 30 with a support floor of one row, a Laplace's refused 16 of 1,200, or 1.3%, and the t's none.
  - *The call:* keep this, or build the other path for `Unweighted` and `Frozen`, at about the size of the pencil's own. My recommendation is to keep it until E1 reads its own refusals, which a block's accept counter gives.
- *Decided, with the alternatives in the docs.*
  - A model whose priors read the snapshot throws at its first move (`BassModel.PriorsReadResponse`). The routes that read nothing are wider than the brief's: a flat or zero intercept with Jeffreys, a declared inverse gamma or a held σ², so E32's zero centre qualifies.
  - A ladder's and a run's composition read the target their sampling chain holds. Having a ladder read that target exposed something the row had not named: a swap left a replica's target behind its state, so a σ² drawn apart was read against the wrong state. A replica now takes the target with the state. SMC and the flat histogram keep the model's.
- *Evidence.*
  - Regression 1,820 of 1,820 on the final tree, 119 of them new, with the eight hashes held. The BASS oracle filter 90 of 90, at below-normal priority.
  - `AMovingResponse_WalksTheSameStructures_OnKYPlusM` holds through the new path on all seven routes. A chain moved in place walks a chain rebuilt per snapshot for 320 sweeps on fifteen routes and forms, draws included.
  - At each snapshot: detailed balance on all 2,263 structures, worst 3.6·10⁻¹⁴, dense algebra to 10⁻⁹, and every structure through the moved path to 10⁻⁹.
  - The joint over (structure, latent) against enumeration on seven readings and a ladder. A move refused exactly where F leaves the structure, 720 of 6,543.
  - 28 deliberate breaks, one at a time in an exported tree, all caught.
- *Cost, under load.* After/before over 72 paired blocks of 200,000 sweeps per fixture, M identical in all. Friedman: wall 1.005 (IQR 0.93–1.13), CPU 1.017 (IQR 0.93–1.11). Direct functional: wall 1.009 (IQR 0.95–1.09), CPU 1.028 (IQR 0.93–1.11). On minima: wall 0.974 and 1.025, CPU 0.978 and 1.018. Others kept about 2 cores busy, and 5 of 144 pairs met the quiet rule. The signs flip by fixture, so no slowdown is resolved, and a change of a few percent could hide in it; the statistical half waits on a quiet window, as the others do.
- *The review found.*
  - **A test with no power.** The first gate case passed with z = 0 holding π = 1.0000, so nothing was tested off the support. A prior that evens z's posterior now gives every value mass, and the gate case must put over 3% where F moves.
  - **A real gap the first mutants showed.** Two of the first 24 survived. The conditional's draw reads Σw where its mean does not, and the chain test had no route on the weighted Gram; and no test held a state no design holds. Each has one now.
  - **A wrong reading of a draw.** A regression test reads a draw on another dataset, and the view took the draw's own n-row response for N rows. A draw's working response now applies only on the response it was composed on (`BassDraw.WorkingOn`). The mean offset bypassed that rule too.
  - **Two scares that were not defects.** A lone chain's z-law read 0.25 against ⅓ under extreme weights; over 5,000,000 sweeps it read 0.327, so that was mixing. A chain held a structure "outside F" at s_min = 0.6, an exact rational tie on the grid's prior Gram, which rounding decides by column order, so the test's floor moved to 0.61.
- *Seam shapes changed, with no recorded ruling.*
  - `BassStructureFit` and `BassDraw` gain `Working`, and `BassTrainingFit` gains `Offset`.
  - `BassTarget`'s constructor takes the working response it read. `NewDesign` and `BassDesign.Load` take one.
  - `WorkingResponse` gains `Moved` and `IsResponseItself`, and its mean is read once.
  - `ReversibleJumpChain` exposes `Moves`. `RidgeSpectrum` and `PriorPencil` gain `WithCross`.
- *Limits.*
  - Nothing fills it yet: latent blocks and `Bass.Fit` are out of scope, so a run's composition under a handed target is unreached.
  - A latent's draw at b = 0 is unbuilt, so a reference rung and SMC's start sit at the model's working response. A ladder with a reference rung is checked for its mechanics only, and SMC and the flat histogram not at all.
  - A conjugate draw under `Weighted` needs [the correction of step 3](moving-response-priors-options.md#what-to-build-in-order). A block that weighs the marginal is exact, as checked.
  - Weights that count rows do not move, and a model with no weights takes none.
  - The training fit's quantiles stay the surface's, and prediction adds no offset back.
  - A σ² drawn apart is not in a snapshot, so a restart under one is not bit for bit.
  - Under a latent scale a retained draw holds its own working response, about 24 MB at N = 1,000 over 1,000 draws, and the collapsed content and conditional means take a design per draw.
  - `BassPriorReferenceTests` is committed, which frees the `SampleLogRate` fix in the report below.
  - The [capability inventory](../capability-inventory.md) stays as read at `c1eda6a`, its U4 included.

### The moving-response priors checked — 2026-11-02

By the main session, on a clean export of `b5ff1e6`: Regression 1,693/1,693 with `BassPathProfile` left out, the BASS oracle filter 90/90 at below-normal priority. Reread and held: the level's energy on all three forms, κ read as τ/κ and its link values, the frozen Gram, the family's mean and both limits, half-Cauchy at a = h = ½, ln IG on v. No reader found that drops a state's rate. The docs of `09774a2` and `d30ab9e` hold against the same derivations.
- *One defect, fixed since ([the report](#the-tiny-shape-evidence-fixed--2026-11-03)).* `ResolvedSigma2Rate.SampleLogRate` draws again when a Gamma draw underflows, so the prior rung and SMC's start draw the Gamma cut to the doubles, while `LogDensity` keeps the whole Gamma's constant. ln Z then reads high by −ln P(r in the doubles), with P(below) about e^(−750h), silently. On the enumeration fixture at a = ½:

  | h | Predicted | Ladder | SMC |
  |---|---|---|---|
  | 2 | 0 | −0.01 | −0.01 |
  | 10⁻² | 0.0006 | not run | −0.02 |
  | 10⁻³ | 0.64 | +0.60 | +0.62 |
  | 10⁻⁴ | 2.62 | +2.56 | +2.61 |

  So it is under 0.01 for h ≥ 0.006. The chip's evidence checks ran at h = ½ and its small-h checks on the target, so neither could see it. By the argument, not run: a ratio between models sharing h and s² cancels it, and chains are unaffected unless a is small too.
- *The fix proposed, not built or tried.* Draw v exactly: on an underflow, v = ln ε + ln(U)/h, the left tail to double precision, with `LogDensity` finite to the left and `PriorAt` reading r at ε. It waits because `BassPriorReferenceTests` asserts `LogDensity(−800) = −∞` and holds another session's uncommitted edits. With a under about 0.01 as well, the marginal below ε would still read r at ε, so that corner stays approximate, within a factor (a + h)/h on that slice.
- *Beside it.* `Gamma.LogCdf` reads −∞ where `Rate * x` underflows with x still in the doubles, at shape 10⁻³ about −0.75 in truth. The window is x within a factor 1/Rate of the doubles' end. Not fixed.

### Priors that stop reading a moving response, steps 1 and 2 — 2026-11-02

`950daca`. Steps 1 and 2 of [the options](moving-response-priors-options.md#what-to-build-in-order), built on the user's go-ahead as a family of selectable forms on one seam. No default changes. The forms, alternatives and trade-offs are in [the docs](../../../docs/Maths/BayesianRegression.md#bass-priors-that-stop-reading-a-moving-response).
- *Built.*
  - `BassPriorReference` on `WorkingResponse`: a frozen level and s², with a rule named per noise (`Moments`, `About`, `Robust`, `Quantile`, `FrozenGram`, `Link`). D03's refusals and the exact-fit floor read it.
  - The level's energy Σwᵢ(ȳ_w − c)², on all three coefficient forms.
  - The coefficient prior's constant κ, read as τ/κ, with Sabanés Bové and Held's link values.
  - `WeightedPriorGram.Frozen`, a prior on a fixed W₀.
  - `Sigma2Prior.Hierarchical` and `HalfCauchy`, with `BassSigma2RateBlock` last in the scale slot.
  - `BassTauBlock` and a chain's composition read the target the chain holds, which closes row 2's seam note.
- *Choices made in building, each argued in the docs.*
  - σ²'s rate went on the state, not the target as the note implied, so swaps, particles, the reference rung and the flat histogram carry it. BARS refuses the prior.
  - A slice step, not Wand's two-block form, since it is exact on every path and tempered form.
  - κ reads through τ/κ rather than moving τ₀, so it reaches a held τ and hyper-g.
  - τ₀ stays 1/N under W₀, so W₀ = (1/κ)·1 is κ under the g-form.
- *Evidence.*
  - Regression 1,701 of 1,701 on the final tree, 76 of them new, with the bit-identity hashes held. The BASS oracle filter 90 of 90, at below-normal priority.
  - `BassPriorReferenceTests`: [the list](../../../docs/Maths/BayesianRegression.md#bass-priors-that-stop-reading-a-moving-response). Dense algebra and the joint normal to 10⁻⁹. Equivariance on the targets to 10⁻⁸ and through chains, with a declared prior failing it as the trilemma says. The family against `HalfCauchy`'s own density and both limits. Enumeration of all 2,263 structures, worst balance 3.6·10⁻¹⁴.
  - The evidence against the enumeration's −16.294: −16.280 on a ladder, −16.327 on a linear ladder, −16.303 on SMC.
- *Cost, under load.* After/before over 72 paired blocks of 200,000 sweeps, each pair's M identical. CPU 0.99 on Friedman (IQR 0.89–1.03) and 1.00 on direct functional (IQR 0.96–1.08); wall 0.96 and 0.99. No pair met the quiet rule, with about 2.4 other cores busy. On minima over four runs: wall 0.98 and 1.02, CPU 0.98 and 1.06. The signs flip by fixture, so no slowdown is resolved at this load, and a change of a few percent could still hide in it; the statistical half waits on a quiet window, as the others do.
- *The review found.*
  - **A real defect a loose tolerance had passed.** `BassParticles.Draw` rebuilt its state from terms and τ, so every SMC particle started at the rate's mean. ln Z read −16.00 against −16.29, inside my first tolerance of 0.5. The offset was the same at 600 × 2 and 4,000 × 4 particles, which is what showed it was no noise. Fixed, and the tolerance is now 0.15.
  - A default run paid 24 bytes per evaluation for the option; now one double.
  - The timing harness could not read this change: CPU time ticks every 15.6 ms, 8% of a 0.2 s block, and the first round's ratios disagreed in sign. `BASS_WEIGHTS_BLOCK` sets a longer block.
  - The note's "chains through `ReplaceTarget`" could not move a response, since it swaps only the target. The tests use a model per snapshot, and the note is corrected.
  - c was already σ²'s shape in this module, so the constant is κ, and the splits are in the [notation ledger](../../notation-ledger.md).
  - Two stale statements about `BassTauBlock`, in the docs and the brief, and an overstatement in the docs' own list of checks.
- *Seam shapes changed, approved with the go-ahead.*
  - `WorkingResponse` gains `Reference`, and `WeightedPriorGram` gains `Frozen`.
  - `BassState` gains `LogRate`, `BassDraw` gains `Sigma2Rate`, and `BassChain`, `BassLadder` and `BassStructureMoves` take a starting rate.
  - `BassResolved` gains `PriorLevel`, `CoefficientConstant` and `Sigma2Rate`, and `BassEvidencePriors` gains `Sigma2Rate`.
  - `ReversibleJumpChain` exposes `Target`.
  - `BassTauBlock` takes no target, and `BassSweepSlots.Blocks` is `(model, first)`.
- *Limits.*
  - A design still cannot follow a response, so a moving response here is a model per snapshot. With no reference the priors read the snapshot, so the row beside this one has to require a reference or a route that reads nothing.
  - `Bass.Fit` takes no observation model. So a reference reaches a model and a chain, not a run; the hierarchical prior, being config, reaches a run.
  - A ladder's and a run's composition read the model's target; only a lone chain reads the one it holds. A σ² drawn apart on a ladder will need that.
  - σ²'s rate is not monitored: R̂ reads σ².
  - Which route each consumer takes is unruled, and steps 3 and 4 are open.
  - The [capability inventory](../capability-inventory.md) stays as read at `c1eda6a`.

### The weighted core and the marginal given σ² — 2026-11-02

The first two rows of [Upgrades to existing code](extensions-brief.md#upgrades-to-existing-code) are built. The third, priors for a moving response, is [written out as options](moving-response-priors-options.md) for the user's call and not built. Commits: `a772fd9` (the core and the form), `ce15a49` (the unit path's support count, the cost block). The forms, alternatives and trade-offs are in [BayesianRegression](../../../docs/Maths/BayesianRegression.md#bass-the-weighted-core).
- *A check that could not hold as written.* "Integer weights equal duplicated rows, through the marginal" conflicts with D03's robust row, the recorded default under weights, in two ways. Its g-prior stays on the unweighted Gram, while a duplicated row's own g-prior reads its copies. And precision weights count N in σ²'s ν, where duplicated rows count Σwᵢ. So both axes are selectable (`WeightedPriorGram`, `WeightMeaning`), with D03's prior and precision weights as the contract's defaults. The check holds in full on `Replication` with `Weighted`, and for the likelihood's statistics under all four readings. Under precision, given σ², it holds up to a constant free of the structure.
- *What D03's prior added.* `FromWeighted`, BARS's block, is built at one τ and b = 1, which τ's slice step and tempered swaps cannot read, so the upgrade needed `PriorPencil`: any (τ, b) in O(p), O(p³) per admitted proposal.
- *What the form added.* `Sigma2Prior.Held`, a fourth kind that leaves the other kinds' JSON and fingerprints as they were, and `BassTarget.Given`, for a σ² drawn apart through `ReplaceTarget`.
- *Beyond the rows.* The pointwise densities now read the working response and its weights, so they read a fixed offset too.
- *Evidence.*
  - Regression 1,625 of 1,625 on the final tree, with the eight hashes held.
  - The BASS oracle filter 90 of 90 on `a772fd9`: 85 serially, then 5 that had hit the R oracle's 120 s limit at below-normal priority. None was a mismatch, and they pass at normal priority in 17–46 s. After `ce15a49`, which only restores the unit path's old loop, the admission parity tests pass 3 of 3; the full filter was not rerun.
  - The checks are listed with [the forms](../../../docs/Maths/BayesianRegression.md#bass-the-weighted-core): dense algebra, σ²'s integral to 10⁻⁹, enumeration of all 2,263 structures, duplicated rows state by state for 400 sweeps, and a held σ²'s evidence at ln Z −14.860 against the enumeration's −14.870.
- *Cost, under load.* After/before over 72 paired blocks: Friedman 0.99 (IQR 0.94–1.05), direct functional 1.07 (IQR 0.98–1.26), or 0.99 on the pairs where others kept at most 4 cores busy. No block met the quiet rule: one other test host throughout, others' cores at a median of 3. So a difference under about 5% cannot be read. The pencil costs 3–4× the closed forms at M ≈ 9, growing as M³.
- *The review found.*
  - A real slowdown of +6 to +14% on Friedman, from a support count that took a nullable lookup and a floating-point add per nonzero element. A first timing block, short and unpaired, could not see it, and hid it behind 2–6× noise; the user's question about a busy machine prompted the paired one that did. Fixed in `ce15a49`.
  - One name for the record (`Weighting`, with `Weights` left as the raw W), and a test name that had gone false.
  - Coverage gaps, now tests: the strict factorization near the gate, the exact-fit floor under W, a ladder's hand-over, and the evidence under a held σ².
- *Seam shapes to report, since they have no recorded ruling.*
  - `WorkingResponse` gains `Meaning`, `PriorGram`, `Weighting` and `ResponseMean`.
  - `BassDesign`'s constructor takes the intercept, which only the pencil reads.
  - `Residual` reads W^½r.
  - `BassEvidencePriors.Sigma2` is nullable, for a held σ².
  - `BassTauBlock` keeps the target it was built with, so a sampled τ beside a moving σ² reads a stale σ². The fix is for it to read the chain's current target, which the chain does not yet expose.
  - Option B for the third row would add a reference to `WorkingResponse`.
- *Limits.*
  - Nothing moves yet: weights and the working response are fixed per model, and the row beside the third, a design that follows them, is not built.
  - `Bass.Fit` takes no observation model, so a weighted run reaches a model and a chain, not a run.
  - Weights with a Q other than the identity are still refused.
  - A direct design under weights costs O(N) per inner product.
  - The gate reads W, so F moves with a latent scale (flagged in the options).
  - `RidgeColumns` and its kin now also hold a pencil, logged in the [notation ledger](../../notation-ledger.md).

### The three chips checked — 2026-11-02

By the main session, on a clean export of `608964c`: Distributions 272/272, R parity 22/22, Regression 1,617/1,617 with `BassPathProfile` left out and `BassBitIdentityTests` in.
- *The inventory (`1a38994`).* Nine sampled citations hold at `c1eda6a`. One correction, to its marks: N1, N2, N4 and part of N3 said the model forbids what is a choice of prior or an approximation of the method, by the same reasoning its U15 gives BGP. They are now U20 (a proper null-space prior, for BAPS's evidence and full prior draw) and U21 (exact coefficient draws for Poisson, Whittle and Gibbs-loss BARS), so it marks 21 upgrades and 5 gaps true of the model as stated. Not checked: the claims it marks as its own readings. Its commit is attributed to Sonnet 5.5, which is what that session ran on.
- *Truncated draws and quantiles (`1f382ab`).* Reread and held: each family's coordinate, peak, spans and offset, the hat's tangents and cut, the bound, the mixture's pick. But a fuzz at random parameters and intervals found four defects its seventy chosen cases and R parity had not, all in the far range and all silent:
  - `HalfCauchy.SampleOn` returned the lower bound on every draw past 10¹⁶ scales, where π/2 − atan a cancels to 0, and on a finite interval past 10¹⁵⁴, where ab overflows;
  - `HalfCauchy.LogPdf` read −∞ past 10¹⁵⁴ scales, as z² overflowed;
  - `Weibull.InverseLogCdf` and `LogCdf` underflowed below ln p ≈ −708 with the root still in range;
  - `Mixture` returned a finite number for a root beyond the double range.
  Corrected in `0acb3f6`, whose message overstates one thing: no draw or density that was right changes, but a quantile can move by rounding, through the solver's midpoint and the Weibull's log form. The fuzz is now `DistributionSweepTests`, seeded. After them: Distributions 280/280, R parity 22/22, bit identity 8/8. The Regression suite was not rerun on the fixes, which touch quantiles and those far-range values only.
- *The weighted core and the marginal given σ² (`a772fd9`).* Reread and held: `PriorPencil`'s whitening, solve, draw and determinant; `LogMarginalGiven` and its slope; the two weight meanings. `TheBlock_MatchesDenseAlgebra` is a real check, building S and the determinant term from the columns with a dense Cholesky and never touching the pencil. Still that chip's: its docs are uncommitted, its before-and-after timing is unresolved, and the [moving-response options](moving-response-priors-options.md) wait on the user.

### Truncated draws and the rest of the quantiles — 2026-11-02

The "Truncated draws" and "Quantiles for the rest" rows of [the upgrades](extensions-brief.md#upgrades-to-existing-code), built ahead of E20 and F-R1's `s2.lower` arm (`1f382ab`). The alternatives and trade-offs are in [ProbabilityAndInformation.md](../../../docs/Maths/ProbabilityAndInformation.md#truncated-draws-and-quantiles).
- *The shape.* `IContinuousDistribution` gains `Sample`, `SampleOn(rng, lower, upper)`, `Quantile(p)`, `InverseLogCdf(logP)` and `InverseLogSurvival(logS)`, on every family, `Uniform`, `Mixture` and `Truncated<T>` included. `Truncated<T>.Sample` is its parent's `SampleOn`. The brief expected per-family samplers; they are, behind one interface member, so no family can lack a draw.
- *What made one sampler enough.* Each family is log-concave in one coordinate at every parameter value: z, ln x, logit x, asinh(t/√ν). So a gamma below shape 1, an inverse gamma's upper tail and a t's polynomial tail need no branch, and the same hat and the same solver serve all of them.
- *Draws.* A hat of tangents, against h = ln f − ln f(anchor). It reads no CDF or normalizer. Weibull and half-Cauchy invert in closed form.
  - *Compared:* Robert 1995 and Botev's `trandn` in one dimension, measured. Minimax tilting is multivariate, and reduces to `trandn` here, so it is left for E22. Inversion costs 2–9× the hat and loses a narrow interval.
  - *Chosen:* the hat. It proves over half accepted and measures 87% or more, where Robert drops to 63% and `trandn` to 31%.
- *Quantiles.* Newton on the smaller log tail in that coordinate, where the tail is concave, from the tail's asymptote or the bulk's approximation.
- *Evidence.* The distributions project passes, 272 tests. The oracle project's 22 new parity cases pass (`r/oracles/distribution_quantile_oracle.R`, base R only). `BassBitIdentityTests` holds on every hash, so no existing draw changed.
  - *Draws:* a chi-square against the density and the interval moments on 70-odd cases, N(0, 1) on [10, ∞) and [10, 10.0001] among them; acceptance counted from the generator with no hook in the sampler; and exact truncated means and variances from R's log tails, which never read a draw.
  - *Quantiles:* closed forms to ln p = −10³⁰⁰; 418 round trips through R's `p*`; 361 gaps from R's `q*`, the worst 1.4×10⁻¹⁰, at ν = 10⁵, where this side's ln B cancels.
- *R's flaws, corrected rather than ported.* `qt` stops 1% short at ν = 1.5 and tail −700, and overflows to ∞ at ν = ½ and tail −230, where t = 1.03×10¹⁹⁹. `qbeta` floors at 2⁻¹⁰²³ where the root is below the double range.
- *What the review found and fixed.* Below ln p ≈ −10¹⁵ the solver's slope is all rounding, so the normal and gamma now return their asymptotes there. `Truncated<T>.Quantile(0)` missed the bound by an ulp. Test seeds from `HashCode.Combine` were random per process.
- *A contract moved.* `Normal.Quantile` now refuses a p outside [0, 1], as the interface says, where it returned NaN; no caller relied on it. `Normal.StandardQuantile` keeps its NaN, for its callers.
- *Limits.*
  - A hat draw costs about 3× Robert's on a normal half-line, probit's case: about 100 ns against 30, the 5th percentile of short rounds on a shared machine, so rough. A closed half-line path could sit behind `Normal.SampleOn` if E20's latent block profiles there; nothing has measured that.
  - Quantiles inherit the tails' accuracy, and throw above a shape near 10⁸, by the series' iteration cap; the draws read no tail, and hold at shape 10¹² and ν = 10¹⁴.
  - `Empirical` has no draw or quantile.
  - Dagpunar 1978 and Philippe 1997 are unmeasured, with both citations checked by abstract only. Gilks and Wild 1992, Prékopa 1973, AS 91, AS 109 and Hill 1970 are from memory.
  - F-R1's `s2.lower` still needs its derivation and enumeration checks. E20 still needs its latent block.

### Complete values declared — 2026-11-02

`4e3ea95`, by the main session, at the user's word.
- *Built:* `ISamplingTarget.HasCompleteValues`, off by default, with `completeValues` on the two factories. `LadderEvidence.FromLadders`, `FromSmc` and `Reweighting.LogEvidence` refuse without it, which `ParallelTempering`, `SmcResult` and `MulticanonicalDraws` carry. `BassTarget`, `BassParticleTarget` and `SplineTarget` declare it, the last through a new `IObservationModel.HasCompleteMarginal`.
- *No run's numbers move:* BASS's three copies of the condition (`BassRun`, `BassSmcRun`, `BassFlatHistogramRun`) now read the one property, with τ's prior kept apart, since a flat-histogram ratio from b_lo > 0 under an improper τ is legitimate. The refusal only reaches targets that do not declare, and the toys that read an evidence now do.
- *Suites:* on the shared checkout, with other chips' uncommitted work in it: Samplers 285/285, Regression 1,546/1,546, BASS oracle 90/90. On the commit alone, in a clean export: `LadderEvidenceTests` with `ReweightingTests` 21/21, `BarsEvidenceTests` with `BassTemperingTests` 97/97. The full suites were not run on the export.
- *The commit:* the BASS core chip had already extended the property in its uncommitted `BassTarget.cs`, for a held σ². The commit takes only the one property from that file and leaves its working copy alone, so the chip's version lands with its own commit.
- *Not covered:* `SmcResult.LogNormalizerRatio`, and the two estimators that take draws directly. The inventory's slope flag (U8) is not built, and the slope is still probed by catching an exception.

### The capability inventory written — 2026-11-02

[The inventory](../capability-inventory.md) is E33's first part. BARS, BAPS, the hierarchical P-spline, BGP and BASS and the pieces under them are read against nine questions (marginal likelihood, weights and offset, tempering, prior draws, conditional draws, derivatives and smoothness, predictive draws and pointwise densities, resume, evidence), each answer with file and line. Nineteen gaps are marked upgrades and six true of the model itself. E33 stays unticked, since its declared capabilities and four verdicts are not built; the inventory's last section says what earns a declaration.
- *The table.* Its seven rows were checked: every claim holds, the weighted-Gram row is corrected (the shared block takes weights only at one τ and b = 1, with no slope), and five gain facts the reading found. Eleven rows are added: a moving working response in BASS's design, BARS's integrated-form slope, path sampling on BASS's integrated form, predictive draws outside BASS, derivatives of a fit, declared smoothness, a power on the likelihood, evidence outside Normal BARS and BASS, exposed prior draws, weights and offset in BARS, and the same for BAPS and the hierarchical model.
- *A correction outside the table.* The free-energy entry said BASS does not produce the integrated form's slope. The target does (`BassTarget.cs:259-266`); the ladder declines to trace it (`BassLadder.cs:218`). The gap in the other direction is BARS's: its integrated `NormalMarginal` supplies none, though its τ is fixed and the slope at b = 0 would be finite.
- *Found on the way.* Poisson, Whittle and Gibbs loss ignore `weights` without a word. BARS's resume test omits `Evidence` and has no setup with the integrated form and a likelihood-path reference rung, the default that reports one. BGP is the one family with a marginal given σ², at s² = 1. The brief's "BARS's guard was not read" is answered: it tests the model type and the prior kind, and needs no normalization test.
- *Limits.* No code changed and nothing ran. Other sessions' uncommitted work (a declared-completeness change and a weighted-core timing test) was in the tree and is not reflected, so the inventory's U8 and U1 are the rows to re-read when it lands. Not read: BASS's structure prior and proposals, the calibration folder, the ladder, SMC and ensemble internals, the penalties beyond the cited lines, the tests beyond the resume ones' names and `BarsResumeTests.cs`, and `docs/Maths`. Claims that are readings, such as the weighted block's eigen-factor and the tempered Gibbs conditionals, are marked in the file.

### The post-BASS notes read against the code — 2026-11-02

[Gemini's digest](../discussion/gemini-post-bass-feedback.md) of the Perplexity threads, checked against the code and the threads by the main session. New: [Upgrades to existing code](extensions-brief.md#upgrades-to-existing-code), E31, E32, E33 and a paragraph under cross-engine evidence, each tied to the idea it came from. What the lineage review already held is not repeated.
- *Corrections to the digest.* BASS is not smooth by default, since `Degree = 1` has kinks. `BlockScan` snapshots do not copy state in BASS or BARS, which pass an identity and save state themselves. Loss-based updating is not for later, since `GibbsLossModel` is already a Gibbs posterior. BCF's `aᵢ = zᵢ` is one coding. The forest-tempering caveat needs no forest-integrated form, by a reading untested.
- *What it dropped.* The threads' order, with trees last, and their scope question.
- *Limits.* No code changed, and the from-memory citations are marked where they appear.

### The lineage review landed — 2026-11-02

A discussion of the Francom lineage, recorded after a sweep of the thread found most of it recorded nowhere. New: the [lineage review](lineage-review.md) (the five axes, separability, F-X1 as a data-augmentation seam, the simulator workflow, a proposed build order, reading leads, the notes against it) and entries E19–E30. E1, E8 and E11 gain what the first reading found. Nothing built, nothing scheduled, and most of the sources beyond the five the review names are unread.
- *A correction from checking the code.* Probit was called free on unit weights. It also needs σ² held at 1, the marginal given σ², and priors that stop reading the response (E20).
- *Open.* The build order is proposed and undecided, and the exact-against-sampled consumer shape is to settle before E8.

### Conformal Bayes and its curve-level relatives read — 2026-10-29

Raised by the user: is conformal Bayes already in, and how would it reach the full posterior and combine with BARS?
- *Already in:* OQ17's seam, but not the conformal step, which was deferred to the SPC × BARS work.
- *New entries:*
  - E16, the consumer, with conformal predictive distributions;
  - E17, simultaneous bands for whole curves;
  - E18, posteriors for readouts by resampling units (BayesBag, martingale posteriors), with E15 as the calibrated variant.
- *The references' details* were checked on 2026-10-29. Fong and Holmes have a local copy; the others are still to be read.
- *Not here:* Lei, Rinaldo and Wasserman's conformal clustering trees, a pseudo-density hierarchy with finite-sample guarantees, belong to SPC, beside the spin-glass threads' Thread 7.

### Bissiri, Holmes and Walker read — 2026-10-29

*A general framework for updating belief distributions* (2016, *JRSS B* 78; arXiv 1306.6430), read against the port.
- *Already absorbed:* the tempered rungs as power-likelihood updates (E13), and BARS's `GibbsLossModel` on a unit-information scale.
- *New entries:*
  - E14, semi-modular calibration, from §4.2's partial-information update, which is also the cut posterior's coherence argument;
  - E15, calibrated intervals for loss-based posteriors, from §3.4, for BARS's Gibbs-loss route first.
- *A framing for F-R1:* §3.3's hierarchical loss reads a deterministic simulator's σ² as a learned loss weight.
- *The follow-ups* the entries name (Carmona and Nicholls 2020; Holmes and Walker 2017; Lyddon, Holmes and Walker 2019; Syring and Martin 2019) were named from memory at first; their details were checked the same day (the brief's entries give each title, venue, DOI and arXiv number), and they are to be read before either entry starts.
- *Not here:* §5.2's model-free clustering, which belongs to SPC.

### Temperature integration read — 2026-10-28

From Farhi, Hed, Bon, Caticha, Mak and Domany (2013, *Physica A* 392), read for SPC (`issues/spc/opus-spin-glass-threads.md`, Thread 5). BASS's evidence already has the paper's shape, and its reference rung carried back to SPC as the closure at β = 0. New entry E13 (the free-energy curve over the ladder: tempered posteriors, WBIC, the real log canonical threshold). E9 gains why the engines compare at all and the Bayes-factor route under the improper options. Nothing built. Liang et al. 2008 and Watanabe 2013 are checked against the library PDFs (`supellex-private/PDFs/Bayesian`: `Liang2008_StatAssc103.pdf`, `watanabe13a.pdf`); the check added E9's base-model condition and E13's exact regime (b = b₀/ln N). Bissiri, Holmes and Walker 2016 is checked too (`bissiri2016.pdf`, the arXiv 1306.6430v2 copy): their loss weight (eq. 8, §3.1) is the module's inverse temperature b, and an uncited SafeBayes aside is dropped. The entries use the module's symbols (b, θ, ℓ), with the authors' symbols named once where a reader needs them for the source.

### GBMARS's global scale derived — 2026-10-26

From the user's derivation. With u = √w and zᵢ = u·vᵢ, u has a GIG full conditional for any GIG mixing. The prior must go on u (Gamma, inverse-gamma or 1/u) and the coefficient prior's variance must scale as u, or not at all. Variance ∝ w leaves the family, and variance ∝ u puts y's units into Σ₀ or τ (OQ19). The scale-equivariant alternative is log-concave in log u and is left to a slice step. A numerical check is in E1's brief entry. E1 stays deferred.

### Andros et al. read — 2026-10-26

From arXiv 2406.18751 and 2602.09351: new entry E12 (correlated residuals along the output index); E9 gains held-out proper scores (the energy and interval scores); E10 notes the functional GP as a comparison point. Sketching and the functional GP are not adopted, with reasons in the hub's register.

### Elastic emulation added — 2026-10-26

New entry from Francom, Tucker, Huerta, Shuler and Ries (arXiv 2305.08834), with its prep in the outputs segment (O-X5).
