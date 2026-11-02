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
