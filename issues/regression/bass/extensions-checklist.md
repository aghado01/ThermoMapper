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

### The three chips checked — 2026-11-02

By the main session, on a clean export of `608964c`: Distributions 272/272, R parity 22/22, Regression 1,617/1,617 with `BassPathProfile` left out and `BassBitIdentityTests` in.
- *The inventory (`1a38994`).* Nine sampled citations hold at `c1eda6a`. One correction, to its marks: N1, N2, N4 and part of N3 said the model forbids what is a choice of prior or an approximation of the method, by the same reasoning its U15 gives BGP. They are now U20 (a proper null-space prior, for BAPS's evidence and full prior draw) and U21 (exact coefficient draws for Poisson, Whittle and Gibbs-loss BARS), so it marks 21 upgrades and 5 gaps true of the model as stated. Not checked: the claims it marks as its own readings. Its commit is attributed to Sonnet 5.5, which is what that session ran on.
- *Truncated draws and quantiles (`1f382ab`).* Reread and held: each family's coordinate, peak, spans and offset, the hat's tangents and cut, the bound, the mixture's pick. But a fuzz at random parameters and intervals found four defects its seventy chosen cases and R parity had not, all in the far range and all silent:
  - `HalfCauchy.SampleOn` returned the lower bound on every draw past 10¹⁶ scales, where π/2 − atan a cancels to 0, and on a finite interval past 10¹⁵⁴, where ab overflows;
  - `HalfCauchy.LogPdf` read −∞ past 10¹⁵⁴ scales, as z² overflowed;
  - `Weibull.InverseLogCdf` and `LogCdf` underflowed below ln p ≈ −708 with the root still in range;
  - `Mixture` returned a finite number for a root beyond the double range.
  Corrected in `0acb3f6`, none changing a value that was right, and the fuzz is now `DistributionSweepTests`, seeded. After them: Distributions 280/280, R parity 22/22, bit identity 8/8. The Regression suite was not rerun on the fixes, which touch quantiles and those far-range values only.
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
