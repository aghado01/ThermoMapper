# BASS extensions: checklist

Companion to [extensions-brief.md](extensions-brief.md). An item is ticked when its report lands. Reports go under "Reports", newest first, each with its date, commits, evidence and remaining limitations.

## Checklist

### Entries

- [ ] **E1** Generalized-hyperbolic noise (GBMARS), with OQ20's BARS side.
- [ ] **E2** Robust centring and bases (OQ16).
- [ ] **E3** Free-knot tensor-product surfaces (ABMRS).
- [ ] **E4** Additional proposal kernels.
- [ ] **E5** KHAOS's modified g-prior.
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
- [ ] **E16** Conformal Bayes, the consumer, with conformal predictive distributions.
- [ ] **E17** Conformal bands for whole curves.
- [ ] **E18** Posteriors for readouts by resampling units.
- [ ] **E19** Polynomial-chaos terms.
- [ ] **E20** Data-augmented observation models.
- [ ] **E21** A prior-precision builder, with ridge blocks.
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
- *Cost, under load.* After/before over 72 paired blocks of 200,000 sweeps, each pair's M identical. CPU 0.99 on Friedman (IQR 0.89–1.03) and 1.00 on direct functional (IQR 0.96–1.08); wall 0.96 and 0.99. No pair met the quiet rule, with about 2.4 other cores busy. On minima over four runs, the signs flip by fixture, within about ±2%, so no slowdown is resolved at this load.
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
