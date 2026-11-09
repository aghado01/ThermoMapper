# GBASS's noise on the BASS core: E1's build (chip brief)

Date: 2026-11-04. Status: written by the main session at the user's request; not yet commissioned. It turns [E1](extensions-brief.md#generalized-hyperbolic-noise-gbmars) and OQ20 into a diff against the engine and an order to build it in. E1 keeps the derivations and the defects; this brief keeps the pieces, the slices and the gates. Reports go under "Reports" below, newest first.

## Why

The user's reading (2026-11-04): GBASS is BASS with a different noise model, not a new engine. Its noise is a normal variance-mean mixture, εᵢ = √w(βvᵢ + √(c·vᵢ)·ζᵢ), and F-X1 and F-X2 were built for it. The structure chain, the term table, the design and the structure priors stay as they are. What arrives is a working response at a latent state, the blocks that move that state, and the noise's own law wherever a reader assumed a Normal. The core already fits y − β√w·v under weights 1/(c·v), with GBASS's w as σ² and GBASS's τ as the reciprocal of the core's.

**The order is the user's call (2026-11-04):** GBASS first, with CRAN's package in the oracle harness before any of it, then KHAOS.

## The oracle: CRAN's GBASS 2.0.1

Pinned in `r/renv.lock` through `r/DESCRIPTION` with its imports `BASS` 1.3.1 and `GIGrvg` 0.8. `library(GBASS)` from `r/` reaches it. The source oracles of R BASS are apart: they read the checkout at `2671014` through `bass_source`, and the CRAN `BASS` does not touch them.

- **CRAN's package is the reference; the clone was the first reading.** The clone at `D:\aghado01\packages\bayesian-related\GBASS` (`46c9afd`, 2026-07-27) postdates CRAN's tarball (packaged 2026-04-22). With comments and roxygen removed, the R code differs in four places: CRAN exports `rmhn` with `rMHN` as its alias; the clone adds a default `w_prior$prop_sigma` for the skewed case, which the GIG path does not read; the clone moved `gbass2bass` to its own file and added knot strategies; and `stefanski.R` is in the clone alone. So E1's reading and its defects hold for CRAN's code, and E20's Stefanski route has no CRAN oracle.
- **Licence.** CRAN registers BSD-3 with its `LICENSE`, which settles E1's question against the README's badge. `GIGrvg` is GPL (≥ 2) and R BASS GPL-3, so nothing here copies their code: samplers come from the papers (Hörmann and Leydold 2014; Sun, Kong and Pal 2021), and R is run for comparison.
- **Three kinds of parity, each only where the source is right** (CONTRIBUTING's "Sources are checked, not trusted"):
  - *Primitives.* `GIGrvg::rgig` and `rgig2` against the GIG's exact moments and ours; `rmhn` against ours and quadrature; base `besselK(x, nu, expon.scaled = TRUE)` against the log Bessel K.
  - *Conditionals at a frozen state.* At a given (B, a, v, w, β, γ, τ): v's GIG parameters, β's and γ's normals, τ's inverse gamma, w's quadratic term and modified half-normal parameters. `_common.R` already evaluates a source's expressions as written; `gbass`'s body is one function, so the same walk applies.
  - *Posterior readouts, never traces.* The chains differ by design, so `tbass`, `qbass`, `nwbass` and `hbass` meet the engine on the fitted surface over a grid, w, β, γ and M, each difference below stated beside the comparison. τ is held equal where the ridge's scaling would otherwise separate them (a raw uncentred ridge with its intercept inside, against centred columns at unit root mean square).
- **The read defects get run.** E1's list came from reading. The chip confirms each against the installed package, or strikes it.

## What is already built

| GBASS needs | Here |
| --- | --- |
| The structure chain | `Change`, `PruneOrGrow`, `VariableChange`, the term table, D01's grid |
| y − β√w·v and its weights, moving each sweep | `WorkingResponse.Moved`, `BassStructureMoves.Follow`, `BassChain.Follow`, the weighted core |
| The structure's marginal given w | `GPriorMarginal.LogMarginalGiven`, `BassTarget.Given`, `Sigma2Prior.Held` |
| The sweep's order | `BassSweepSlot`, D04's robust row: structure, τ, a, the latent scales, then σ² and the skew |
| Priors that stop reading the working response | `BassPriorReference` (`Moments(y)`), `RidgeForm`, `GPriorCentre.Zero`, `WeightedPriorGram` |
| Conjugate updates and exact draws | `GigParameters`, `ConjugateUpdates`, `GaussianConditional`, `InverseGamma`, `SampleOn`, `SliceSampler` |
| `hbass`'s half-Cauchy on √w, conjugately | `Sigma2Prior.HalfCauchy`, `BassSigma2RateBlock` |
| A draw's latent reading, kept | `BassDraw.Working`, `BassDrawFit.MeanOffset`, `BassTrainingFit.Offset` |
| Student t's conditional and marginal | `StudentTWeights` on BARS's side, the pin for slice 1 |

## The diff

### New

| Piece | Where | What it takes |
| --- | --- | --- |
| GIG distribution and sampler | `Maths.Distributions` | Hörmann and Leydold 2014, with `GigParameters` for the algebra and its limits. A residual of exactly zero gives b = 0 at p > 0, the Gamma limit |
| Modified half-normal | `Maths.Distributions` | Kernel x^{α−1}exp(−βx² + γx), Sun, Kong and Pal 2021. Every branch checked, including the 0 < α ≤ 1 series branch that E1's scratch left unchecked |
| Log Bessel K of real order | `SpecialFunctions` | In logs, over far ranges. Only Normal-Wald (K₁) and general GIG mixing read it: t, Laplace and asymmetric Laplace have closed densities |
| Generalized beta prime | `Maths.Distributions` | Density and draw. v's conditional is not conjugate, so `SliceSampler` on ln v, one step per observation |
| The mixing law | `Maths.Regression`, shared with BARS | Per observation: the prior draw, the conditional draw given (r, w, β) at rung b, the marginal log density, E[v] and Var[v]. `IWeightResampler`'s `Resample` and `MarginalLogDensity` are this seam at β = 0, so OQ20 takes the same members |
| The variance-mean observation model | Beside `GaussianObservation` | `Working` from v, β, w and c, and the blocks: a (with σ² under a pure scale mixture), v, then σ² and the skew. A name true across the family: GBP mixing is not generalized hyperbolic |

### Upgrades

| Piece | Today | Upgrade |
| --- | --- | --- |
| `IBassObservationModel` | `Working` and `Blocks` | Also the noise law: pointwise density, predictive noise, the mean's shift. Gaussian is the default, so its members do not change |
| `WorkingResponse.Moved` | Offset and weights | Also carries the latent record (decision 1 below) |
| `BassPointwise.Retained` | Hard-coded N(y \| f, σ²/w) | Reads the noise law, integrating vᵢ (OQ17). The collapsed predictive refuses a non-Gaussian law, with the reason, until it has its own check |
| `BassPrediction.Predict` | No offset | Three readings per law: the surface f, which is the q-quantile under `qbass`; the mean f + √w·β·E[v]; predictive draws with v and ζ. E20 reads the same |
| `BassTrainingFit.NoiseSd` | Null under weights | The law states it: √(w(β²Var[v] + c·E[v])), and null where the variance is infinite (t at ν ≤ 2) |
| `Bass.Fit`, `Resume`, `Fingerprint` | No observation model: a run through them is Gaussian | A declarative noise on `BassConfig`, left out of the JSON when absent as `Sigma2Prior`'s `Value` is, so every existing fingerprint holds |
| `BassSmc`, `BassFlatHistogramRun` | One model target (`BassSmc.cs:151` refuses another) | Read the carrier, so a particle holds its own latents and w (slice 5) |
| BARS `StudentTWeights` | Its own conditional | Through the shared law, numbers unchanged, then skew and quantile BARS with an offset BARS does not yet take (OQ20, the checklist's row) |

### Left alone

The moves, the proposals, the term table, the design kernels, D01's knots, the composition's exact draws (which read the target the chain holds), and `UseTables`, whose rule already turns prefix tables off under a latent block (`651cb90`).

## GBASS's wrappers as configuration

| GBASS | Mixing | c | β | w's prior and draw here |
| --- | --- | --- | --- | --- |
| `tbass(ν)` | GIG(−ν/2, 0, ν) | (ν−2)/ν, so w is the variance; c = 1 is `StudentTWeights`' t | 0 | 1/w; integrated in the structure step, then conjugate (D04) |
| `qbass(q)` | GIG(1, 2, 0) | 2/(q(1−q)) | 1/q − 1/(1−q), fixed | IG(0.1, 0.05); modified half-normal |
| `nwbass` | GIG(−½, γ², 1) | 1 | Sampled, with γ | IG(0.1, 0.05); modified half-normal |
| `hbass` | GBP(1, ½, ½) | var(y), read once through the reference | 0 | Half-Cauchy on √w by the hierarchical route, where GBASS takes a log random walk |

## Decisions

**Recorded, and executed here.**
- *D03's robust row:* the unweighted g-prior is the default; the ridge with a sampled τ is the GBASS-parity option.
- *D04's robust row:* the order above; a each sweep; σ² integrated out of the structure step under a pure scale mixture and held otherwise, since under skew w sits in the offset.
- *OQ17's scope rule:* the pointwise density integrates vᵢ.
- *Corrections, not parity variants:* E1's defects, and the γ > 0 piece drawn with a truncated draw.

**Recommended, and accepted by the user on 2026-11-05 when the chip for slices 0 and 1 started.**
1. **The latent record rides `WorkingResponse`; this is the decision on `BassState` that the seams left to E1.** v, β, γ and, when held, w go on a record that `Moved` takes.
   - *Why there.* `BassStructureFit.Working` and `BassDraw.Working` already carry the working response through snapshots, swaps, SMC's copies and retention, and the moving-response tests already find a chain's latent by it (`LatentBlock`). It is by definition the reading of the response at a latent state.
   - *w too.* Under skew a move of w moves the offset, so `ReplaceTarget(Given(w))` and `Follow` are two updates that must stay in step, with a wrong density between them. On the record they are one move, and `BassTarget.Evaluate` reads w from the state, as it reads the statistics and σ²'s rate. SMC and the flat histogram, which hold one target, then work without a target per particle. `Given` stays for E32's outside σ².
   - *Within a sweep.* a, and σ² under a pure scale mixture, sit on the record from their draw to the last block that reads them. Being redrawn each sweep, they are not state.
   - *Alternatives.* A field on `BassState` changes every snapshot and duplicates what the working response carries. State held by the observation model breaks on swaps and forks, since one model serves every chain.
2. **v and the modified half-normal first.** This matches the reference and stays scale-equivariant under variance ∝ w. α = N + K − 2p_w, with K the coefficients whose prior variance scales with w: M + 1 with the intercept inside the prior, M under the flat intercept (GBASS's form is the first). The user's u and z parameterization (E1, 2026-10-26) joins at slice 7 for interweaving. The three routes and their trade-offs go in `docs/Maths`.
3. **Take GBASS's models, not its chain.** D01's grid, not continuous knots. Pseudo-counts stay 5, not 1, and any mixing comparison states them. The whole-term mutation only if [G1's reading (h)](lineage-gaps-register.md#g1-does-the-sampler-reach-the-posterior) earns it. `npart` symmetric. `lag_beta` is read at slice 3, and adopted only if a free skew is measured to soak up misfit early.
4. **F under a moving W stays on the gate that reads W,** the checklist's open call, until this build reads its own refusals. Every slice reports them.

## The slices

One commit per concern; tests with the code they test; each slice ends with its gates and a fresh-eyes review before the next starts.

| Slice | Delivers | First needs | Its own gates |
| --- | --- | --- | --- |
| 0. Oracle modes | `gbass_*_oracle.R` and their parity tests, the read defects run | None; the package is pinned | The modes run from `r/`; each defect confirmed or struck |
| 1. Student t | Robust BASS | No new distribution: v's conditional is inverse gamma. The carrier, the noise law, pointwise, prediction, the run and fingerprint, the a and v blocks | v's conditional equal to `StudentTWeights`' number for number; large ν reaches Gaussian BASS's readouts; `tbass` readouts |
| 2. Asymmetric Laplace | Quantile regression | The GIG sampler, the modified half-normal, w on the carrier | The skew's sign (E1's fixture); the q-quantile recovered on data with a known quantile; `qbass` readouts |
| 3. Normal-Wald | Skew and tails chosen by the data | Log Bessel K, sampled β and γ, the truncated γ, general GIG mixing | Density against quadrature of the mixture; corrected `nw_est_mom` and `nw_gamma_prior`; `nwbass` readouts; `lag_beta` read |
| 4. Beta prime | "Interpolate or ignore" noise, the Horseshoe | GBP, slice steps on ln v, half-Cauchy w | Pointwise by per-observation quadrature, which costs N × draws × nodes, so opt-in for that path; `hbass` readouts |
| 5. The engine's reach | What GBASS lacks: categorical and functional parts, ladders, SMC, the flat histogram | The tempered v conditional, GIG(p − b/2, a + bβ²/c, b_v + b·r²/(w·c)) in the integrated form (derived, not run), which at b = 0 is the prior draw `docs/Maths` asks for; SMC and the flat histogram on the carrier | The law at each rung by the gates below; a per-observation latent under a direct functional design; evidence where the values are complete |
| 6. BARS (OQ20) | Skewed and quantile BARS | `StudentTWeights` through the law, numbers unchanged, then an offset in BARS | `StudentTWeights`' pinned values; the same oracle modes |
| 7. Interweaving | The z parameterization beside v | The GIG on u, the slice step for the equivariant prior | Same law as v; R̂ and effective sizes reported per unit of work |

The carrier and the noise law in slice 1 also serve E20's Pólya–Gamma and probit blocks, which is a reason to keep their names free of this family.

## Gates common to every slice

- **Gaussian untouched.** The eight pinned hashes unmoved and every existing fingerprint the same: the pinned tier at each commit, the whole assembly at the slice's end ([targeted runs](../../../CONTRIBUTING.md#testing-standards)).
- **The law.**
  - A discrete mixing law on an enumerable design is exact by enumeration per value, the moving-response tests' `Joint`.
  - For continuous v, Geweke's joint-distribution test: the marginal-conditional and successive-conditional simulators, on a design small enough to run long, with proper priors.
- **Recovery from truth outside the model.** Data generated by a process the estimator does not share: heavy tails, known skew, known quantiles. Interval coverage reported.
- **Primitives.** Chosen cases and a seeded random sweep, far ranges included.
- **Refusals and acceptance.** Each block's counters reported at each point, for decision 4.
- **The cost row,** below.
- **`docs/Maths`.** Decided statistical choices land with their alternatives and implications; new distributions in `ProbabilityAndInformation.md`. New writing takes the module's symbols, with the w, σ² and τ crossings logged in `issues/notation-ledger.md`.

## The moving-weight cost row

No cost report prices this chain, and [the efficiency pass](efficiency-pass-brief.md#reports) kept to plain chains. What is known:
- An offset move costs about 10–16 µs and a weight move about 80–140 µs on Friedman, both read before the pass. A quiet base sweep was 22 µs. So at these sizes `Follow` alone would be several plain sweeps: the row is about where a sweep's time goes, not about small factors.
- Prefix tables are off by rule under a latent block. Levers 2–5 sit in code this chain runs, including whitening in lanes, so the 80–140 µs wants re-reading.
- Under `WeightedPriorGram.Unweighted` each admitted proposal costs its own O(p³) (`PriorPencil`). `BassWeightsCostTests` already times the combinations of reading, moves and form.
- R's `tbass` took about 11 ms an iteration at N = 300 and `qbass` the same (2,000 iterations, M about 6, loaded machine, one reading): a reference, not a target.

**To read, at the cost report's base, n 200, M 200 at n = 10,000 and a direct functional point, under both Gram readings, g-form and ridge:** the structure step, τ, a's draw, the N draws of v, w and β, and `Follow`, per sweep and per part. The record's pairs want a quiet window; until one comes, the thread-cycle clock's ratios, marked indicative.

**Retention.** A retained draw shares its `WorkingResponse`, which holds at least the shifted values, the offset, the weights and their logs: about 32 KB a draw at N = 1,000 and 3.2 MB at N = 100,000. Read it in slice 1. If it dominates, keep v and the scalars on the draw and rebuild the rest on demand, since they are functions of those. The generative state stays either way.

**A lever to read after the row, not before:** structure sweeps per latent move. Any fixed cadence is a valid scan, so it is a measured option, and the default stays one each, as D04's sweep has it.

## Not in this brief

- **G1's exact machinery.** With v moving, the structure is no longer a chain on an enumerable space. Whether heavy tails raise or lower the barriers between basins is unread, and it runs only if a decision comes to depend on it.
- **Stefanski's logistic mixing (E20).** In the clone only; its route to logit is E20's.
- **KHAOS.** Next ([the build order](lineage-review.md#build-order)), in a brief of its own, with CRAN's `khaos` 2.1.1 pinned the same way first. It shares this build's lineage and its method, not its pieces: it works on the basis, the proposals and the coefficient prior (E19, E4, E5, E21), which this build leaves alone, and it moves no latent. Once both exist their product costs nothing more: a polynomial-chaos basis under quantile or robust noise.

## Working rules

- On `main`, committing by pathspec with nothing left staged and nothing pushed. Commit messages describe the work, with no attribution lines ([CLAUDE.md](../../../CLAUDE.md)).
- Builds set `ArtifactsRoot` inline into `artifacts\chip-gbass-noise\`, deleted at the end.
- Heavy runs at below-normal priority with capped threads; the machine is shared, and basin-hold's chains were running when this was written. Tests run by [targeted runs](../../../CONTRIBUTING.md#testing-standards) (`scripts/test-select.ps1`): the selected classes per edit, with the pinned tier per commit, the whole assembly once at the slice's end, and the qualification gates as each slice states.
- New beside old; a removal is its own question to the user, with what is lost.
- Findings framed by problem class, not by any one application.
- A follow-up goes in a report with the decision it would change and what each outcome would do. A question with no decision waiting on it is recorded as open, not proposed.

## Report

Under "Reports" below, newest first, per slice: what landed by commit; the gates' numbers; the parity comparisons with their stated differences; refusals and acceptance; the cost row; defects confirmed, struck or newly found; decisions taken beyond this brief; what was left unbuilt and why. `docs/Maths/BayesianRegression.md`'s table of what each extension still writes shrinks as slices land.

## Reports

Newest first.

### The mixture's tail table across a cusp at 0 (chip, 2026-11-08)

**Landed** in `019738d` (the fix, the sweep, the conformal chip's test unskipped) and the commit that carries this report (`ProbabilityAndInformation.md`'s pole entry, the page's entry beside the tail table's). The conformal chip found it in its slice 5 ([its brief](../conformal-layer-brief.md), "Handed on, not fixed here"), and the main session handed it on.
- *The defect.* Under GBP mixing with a drift, `NormalVarianceMeanMixture`'s table read ln F about 10⁻³ off near x = 0, the density there right to 12 digits: at GBP(1, 1, 1.5) and m/s = −0.2, ln F(0) read −0.549285 against the integral's −0.550321, the table's slope 4% low across ±10⁻³·s. BARS's GBP noise reads its log CDF from the table, so its predictive inherited the error near the median; the conformal sets read the density and did not.
- *The cause: the table's grid near 0,* not its interpolation's form or the drift's centring. Where v's density near 0 goes as v^(α−1), the mixture's density goes as |x|^(2α−1) at 0, so the normal score is not smooth there. At α ≤ ½, a pole, the table already stepped into 0 on cells short against |x|. Above it, at a cusp (α < 1) or a kink (α = 1, as here), it took ordinary steps: without a drift the median is 0 and a node sat on it, while with one a cubic cell spanned it. A probe showed the gap vanishing at the nodes and largest across the cell holding 0 (1.0·10⁻³ at 0, 1.8·10⁻³ at −0.03·s), and 0 at x = 0 without the drift.
- *The fix,* as the tail table's was, in the table alone. Every finite left exponent takes the pole's walk: within s of 0 steps of at most ¼ of |x|, ½ within 10⁻²·s, down to 10⁻⁶·s, and 0 as a node. Only a pole's two cells are read by the integral, its slope at 0 infinite; a cusp's slope there is finite.
- *What moves.* The tables of the GBP at pa > ½ and of the variance gamma at p > ½: their nodes within s of 0, and past 0 on the far side. Away from 0 their readings move within the table's stated accuracy; near it they lose the old error. Unchanged bit for bit: the horseshoe and the Strawderman–Berger, poles that already took this walk, so BASS's default noise; every GIG mixture with b > 0, the Normal-Wald, the skewed t and the generalized hyperbolic, whose densities are smooth at 0. A table's build time did not move measurably under load: four affected laws 41 to 73 ms after against 50 to 84 ms before, while the unchanged horseshoe read 86 and 163 ms in the same runs.
- *Checked.* `NearZero_UnderADrift_TheTableHoldsItsIntegral`, a seeded sweep over 40 laws, 30 GBPs with α = pa from 0.52 to 3 and 10 variance gammas with p from 0.55 to 3, drifts m/s within ±1 and scales e^(±2): it reads both log tails at 0, at ±10⁻⁹·s to ±s and at the median against the integral. The worst relative gap was 5.8·10⁻² before the fix (a GBP at α = 0.55, just above the pole), and is 7.9·10⁻⁶ after. The conformal chip's `TheLogCdf_UnderGbpMixingWithADrift_IsItsDensitysIntegral` is unskipped and passes. Per commit: distributions 32 of 32, regression 1,952 of 1,952.

**Follow-ups.** None new.

### The move along β·v, the skewed horseshoe's default (chip, 2026-11-08)

**Landed** in `bf21578` (the switch, its log and the tests) and the commit that carries this report (the page, `BayesianRegression.md`'s GBMARS row, the ledger's anchors). The user approved the move as the default under the skewed horseshoe only ("Default there", 2026-11-08, through the main session); it stays opt-in under every other skewed law. The user's reasons: under the skewed horseshoe the sweep without it fails Geweke's test at b = 1 (4.13 at four times the size) and the move passes (1.43), and doing the job outweighs a run's 1.2 to 1.6 times; under the generalized hyperbolic, the Normal-Wald and the skewed t the sweep passes without it, so the cost buys nothing. b = ½ stays this noise's limit on hot rungs. This replaces the previous report's "stays opt-in" for the skewed horseshoe.

| Piece | Where |
|---|---|
| The switch | `BassNoise.SkewRescale`, now `bool?`: null takes the noise's default, true or false declares it, true under a held skew refused; null left out of the serialized form |
| The default | `BassNoise.SkewRescaleByDefault`, 1.99(1, ½, ½) with β drawn; `RescalesSkew` resolves the two, and the sweep reads it |
| The log | `BassResolved.AutoPicks` gains `Noise.SkewRescale`, on or off with its reason, wherever a drawn skew leaves it to the default |

**What changes, and what does not.**
- *Serialized forms and fingerprints:* a config written before serializes and fingerprints as it did, since null is left out and true is written as before. A skewed-horseshoe config that declares nothing now runs the move, so its draws change under the same fingerprint. Declaring `false` keeps the sweep without it, under a fingerprint of its own.
- *Resume: neither a refusal nor a restore of the old setting, as nothing can reach either.* A BASS checkpoint is the run's live restart state: `EnsembleCheckpoint` holds `BassRun`'s private chain state, live objects, with no serializer, and nothing outside `src/maths` runs BASS. So no checkpoint taken before this change reaches a build with it. Within a build the default resolves from the config, which the fingerprint holds, so a resume runs the sweep its checkpoint ran; the test resumes a tempered fit under the default bit for bit. A checkpoint ever written out would carry the resolved setting.
- *Pinned hashes:* none moved. The bit-identity pins are Gaussian BASS's, and nothing was re-baselined.
- *Harnesses whose readings are the sweep without the move* declare it off, so those readings stand: the reach Geweke's plain forms (slices 5 and 7) and interweaving's measurement (slice 7). The generalized beta prime's sweep test lists the move among the skewed horseshoe's blocks.
- *The generalized beta prime's Geweke gate* under the skewed horseshoe (seed 20261106) now runs the default sweep: worst |z| 1.99.

**Gates.**
- *Its unit tests:* the default resolves on under the skewed horseshoe and off under eight other laws: the horseshoe symmetric and with its skew held, and the skewed Strawderman–Berger, generalized hyperbolic, Normal-Wald and t, the t and the asymmetric Laplace. A declaration overrides it either way, true needs a drawn skew, and the log holds the pick exactly where a drawn skew leaves it to the default. The switch's three forms serialize apart and fingerprint apart. A tempered fit under the default resumes bit for bit with the move running.
- *Per commit,* the selected classes and the pinned tier from a snapshot: 1,627 of 1,627.
- *At the slice's end,* the whole assembly: 2,329 of 2,329.

**Follow-ups, with the decision each would change.** None new.

### The move along β·v, corrected: the coefficients drawn after an accepted move (chip, 2026-11-08)

**Landed** in `1a8459c` (the move, its tests and the gates) and the commit that carries this report (the page's section, the ledger's anchors). The main session's review of `3be066e` found that after an accepted move the blocks that follow read coefficients drawn before it. It asked for a fix, choosing between two exact forms by mixing per second, then a rerun of the skewed horseshoe's Geweke rows and the effective sizes, and a correction of the page and this report. The previous report's readings with the move come from the version with the defect; these replace them. [The page's section](../../../docs/Maths/VarianceMeanMixtures.md#the-skews-rescaling-a-move-along-βv-slices-5s-and-7s-follow-up) holds the forms, the tables and the decisions.

**The defect, and a second one beside it.**
- *Stale coefficients.* The move accepts on ℓ_b, the coefficients a integrated, so it is a partially collapsed step. σ²'s block, the shape's and interweaving's, each drawing given this sweep's a, then read an a drawn at the old (β, v).
- *A proposal read through a.* κ's scale took the likelihood's curvature as b·R/(2σ²), R from the sweep's fit. On the a-marginal that is a proposal depending on a variable the step does not condition on, which is not exact either.

**Two exact forms, built and measured side by side.**

| Form | The fix | On its orbit |
|---|---|---|
| The coefficients integrated | a drawn at the new state after an accepted move (`DrawSweepFit`); the scale's likelihood term b·N/2, from x alone | \|z\| ≤ 1.54, 44 to 46% accepted, a a draw at the new state after every acceptance |
| Given this sweep's a | ln A on the likelihood given a, −(bN/2)·t − (b·R/(2σ²))·(e⁻ᵗ − 1), O(N); the weights refitted only when accepted | \|z\| ≤ 1.97, 44% accepted |

**The choice: the coefficients integrated.**
- *Mixing per second.* On v's mean, which the move moves, it mixes 1.7 times faster under the skewed horseshoe (177 and 226 a second against 98 and 136 at seeds 11 and 23) and alike under the generalized hyperbolic and the Normal-Wald. ln σ², β and M swing by seed in both forms, following the structure's R̂ of 1.3 to 3.8 (given a under the Normal-Wald, 2.14 against the plain sweep on one seed and 0.09 on the other), so they do not separate the forms. The form given a runs 4 to 10% faster.
- *Geweke under the skewed horseshoe,* worst |z| at b = 1, ¾ and ½: the coefficients integrated, 1.93, 2.15 and 9.71 on seed 20261107; 1.43, 2.69 and 4.29 at four times the size; 2.41, 2.59 and 5.18 on seed 20261108. Given a: 4.01, 3.00 and 4.54; 2.23, 1.48 and 5.23; 4.21, 1.96 and 5.70. As first committed: 2.36, 5.70 and 11.13; 4.04, 4.69 and 4.26; 2.02, 5.47 and 3.46. The fix clears b = ¾, which the stale a failed.
- *Under the other skewed laws,* the kept form: the generalized hyperbolic, the Normal-Wald and the skewed t at b = 1, 1.65, 2.20 and 1.48; the generalized hyperbolic and the Normal-Wald at b = ½, 2.58 and 1.29.
- The form given a and its knob are removed. Its derivation stays on the page beside the kept one.

**The default, against the decision rule.** The move was to become the default under a skewed heavy-tailed law if it passed where the scan failed. It passes at b = 1, where the scan alone failed at four times the size (1.43 against 4.13), and at ¾. At b = ½ it fails on every seed and size, as the scan does, though less at four times the size (4.29 against 9.16). Per second it buys v's mean, 3.5 times under the skewed horseshoe, and not ln σ², β or M (0.41 to 1.70), at 1.2 to 1.6 times a run's time. **So it stays opt-in:** it clears the full-power reading the scan failed, but not the hot rungs', and its gain is in v alone. The skewed horseshoe's failures at b = ½ stand as its mixing limit on hot rungs.

**Interweaving reads a conditionally, confirmed.** The scale slot runs σ²'s block, the shape's, then interweaving, each from its conditional given this sweep's a, none integrating it. The only step between the coefficients' draw and them that integrates a is this move, which now draws a afresh whenever it changes x. v's block draws v given a, so a stands after it.

**Gates.**
- *Per commit,* the selected classes and the pinned tier from a snapshot (`scripts/test-select.ps1 -Scope Commit`): 1,618 of 1,618.
- *Geweke gates now* the kept form at b = 1 and ¾ under the skewed horseshoe on two seeds and four times the size, beside the other skewed laws' rows; b = ½ stays a reading.
- *A ladder's swap rates* under the skewed horseshoe (6 rungs to 0.05): the cold pair accepts no swap with or without the move, the others 0.007 to 0.040 without it and 0.009 to 0.042 with it; no round trip either way.
- *The cost row:* a sweep at n 200, both chains at M 2, 275 µs with the move against 91 without (3.0×); at the base point 976 against 448, at M 7 against M 1. A whole run takes 1.2 to 1.6 times as long.
- *At the slice's end,* the whole assembly (`-Scope Slice`): 2,320 of 2,320.

**Follow-ups, with the decision each would change.** None new. Slice 5's reading for the scan alone did not involve the move, and stands as recorded.

### The move along β·v, the skew's rescaling (chip, 2026-11-08)

**Landed** in `3be066e` (the code and its tests) and the commit that carries this report (`docs/Maths`, the ledger, the checklist). The user's go-ahead came through the main session on 2026-11-08 ("Yes, build it"), as slice 5's and slice 7's follow-up. [The page's section](../../../docs/Maths/VarianceMeanMixtures.md#the-skews-rescaling-a-move-along-βv-slices-5s-and-7s-follow-up) holds the derivation, κ's proposal against its alternatives, and the measurement.

| Piece | Where |
|---|---|
| The switch | `BassNoise.SkewRescale`, serialized only when on, refused without a drawn skew; on `BassNoise`, since `BassConfig.cs` was KHAOS's |
| The move | `BassSkewRescaleBlock`, in the latent slot after v's block: (β, v) → (β/κ, κ·v), the weights refitted once, accepted on the sweep's joint with the Jacobian κ^(N−1) |
| κ's proposal | ln κ ~ N(0, h(x)²), h = 2.4/√C(x) from the log target's curvature along ln κ at the state, the Hastings ratio of the state-dependent scale in A |

**Gates.**
- *Gaussian untouched:* the whole regression assembly once at the end (`scripts/test-select.ps1 -Scope Slice` on a snapshot of `3be066e`), 2,317 of 2,317, the pinned tier among them. Before the commit, the selected classes and the pinned tier (`scripts/test-select.ps1 -Scope Commit` from a snapshot): 1,558 of 1,558.
- *Its detailed balance by its statement on a frozen state:* from one state the move stays on its orbit, whose exact law ∝ π(Tₜ(x₀))·e^((N−1)t) the engine's target gives on a grid. 60,000 moves meet its mean and mean square, |z| ≤ 0.87, under the skewed horseshoe and the generalized hyperbolic at b = 1 and ½, at 44% acceptance; the forward and reverse acceptances cancel to 2·10⁻¹³.
- *Geweke with the move under the other skewed laws:* the generalized hyperbolic, the Normal-Wald and the skewed t at b = 1, worst |z| 2.19, 2.61 and 2.03; the generalized hyperbolic and the Normal-Wald at b = ½, 2.96 and 1.38.
- *The run:* a tempered fit with the move and a reference rung resumes bit for bit.
- *The cost row:* a sweep at the base point 2,794 µs with the move against 1,097 without (2.5×), 481 against 307 at n 200 (1.6×).

**The measurement, and the decision it made.** The move does not pass where the scan failed:
- *The skewed horseshoe's Geweke with the move,* worst |z|: on seed 20261107, 2.36, 5.70 and 11.13 at b = 1, ¾ and ½, and 4.04, 4.69 and 4.26 at four times the size; on the control seed 20261108, 2.02, 5.47 and 3.46. The scan alone on seed 20261107 at four times the size: 4.13 at b = 1, 3.13 at ¾, and slice 5's 9.16 at ½.
- *Effective sizes per second,* with the move over without it at seeds 11 and 23: ln σ² 1.03 and 0.76 (skewed horseshoe), 0.53 and 0.35 (skewed GH), 0.23 and 0.73 (skewed Normal-Wald); β 0.22 to 0.82; mean ln v 0.49 to 3.33.
- *A ladder's swap rates* under the skewed horseshoe (6 rungs to 0.05): the cold pair accepts no swap with the move or without it, the others 0.007 to 0.040 without and 0.001 to 0.035 with; no round trip either way.
- **So the move stays opt-in, and the skewed horseshoe's failures stand as that noise's mixing limit under the sweep, recorded with these numbers.** The readings change size and side with seed and size, the scan alone fails at b = 1 at four times the size, and the move holds its law on its orbit and under the lighter-tailed laws, on rungs too.

**Decisions taken beyond the brief.**
- *κ's scale from the state's curvature, not adapted in burn-in.* D10's adaptation needs a scale per replica in its restart state and a burn-in signal to the block, state that would live in `BassStructureMoves` and the run, which other chips hold, and a ladder's rungs would each want their own. The curvature scale needs neither and holds 44% acceptance across laws, rungs and data. A slice step on ln κ was the other route, at some six refits a move.
- *The move needs a drawn skew,* which it moves; a held skew refuses it.

**Left unbuilt, and why.** A default for the move: the measurement decided against it.

**Follow-ups, with the decision each would change.** None new. The skewed horseshoe's slow mixing has had the three moves this brief proposed: the scan, interweaving and the move along β·v. A further move would come with its own measurement on these readings.

### The mixture's tail table past its ends (chip, 2026-11-08)

**Landed** in the commit that carries this entry, after slice 7 and apart from it, at the main session's ask (the defect the conformal chip found in its slice 4).

- *The defect.* Past either end of `NormalVarianceMeanMixture`'s tail table (beyond the 10⁻¹² quantiles) `TailTable.Read` returned the integral on the side asked for, so the side whose mass is near 1 there lost 1 − F's digits and could read above 0: +3.0·10⁻¹⁵ for a generalized hyperbolic, +6.3·10⁻¹⁵ for a Normal-Wald, +2.2·10⁻¹⁶ for a skewed horseshoe. It reached the Normal-Wald, the skewed t, the generalized hyperbolic and the GBP mixtures.
- *The fix.* Past an end, that side is ln(1 − e^(far tail)), the far tail integrated on its own side. Inside the table both sides already came from one normal score, read as ln Φ(±z) with `Log1mExp`.
- *Checked* (`GeneralizedBetaPrimeTests.PastTheTablesEnds_TheSideNearOne_KeepsTheFarTailsDigits`): at 8 to 1,000 scales on both sides, under the four families, no log tail above 0, and the near side equals ln(1 − e^(far tail)) exactly. Before the commit, the distributions' selected class (31 of 31) and the regression selection with the pinned tier (1,724 of 1,724).
- `docs/Maths/ProbabilityAndInformation.md`'s table sentence says how the ends are read.

### Slice 7: interweaving, the z parameterization beside v (chip, 2026-11-08)

**Landed** in `211d768` (the code and its tests) and the commit that carries this report (`docs/Maths`, the checklist). [The page's section](../../../docs/Maths/VarianceMeanMixtures.md#interweaving-the-z-parameterization-beside-v-e1s-slice-7) holds σ's conditional given z, its draw against the exact GIG route, and what it buys.

| Piece | Where |
|---|---|
| The switch | `BassNoise.Interweave`, serialized only when on; `BassConfig.Validate` refuses it with a held σ² |
| The step | `BassInterweavingBlock`, last in the scale slot: z = σ·v held, σ by one slice step on ln σ, then v = z/σ |
| v's prior density per law | `IMixingLaw.LogPrior` on t, the exponential, the GIG (its cached normalizer) and the GBP |

**Slice 7's gates.**
- *Gaussian untouched:* the whole regression assembly once at the slice's end (`scripts/test-select.ps1 -Scope Slice` on a snapshot of `211d768`), 2,241 of 2,241 in 13 minutes, the pinned tier among them. Before the commit, the selected classes and the pinned tier (`scripts/test-select.ps1 -Scope Commit`, from a snapshot of HEAD and this slice's files): 1,689 of 1,689.
- *Same law as v:* the step alone, σ moved by it only (β held at 0.5, σ²'s block out), Geweke at b = 1: the skewed horseshoe 1.75, the generalized hyperbolic 2.31. The interwoven sweep: t(5) at b = 1 and 0.4, 1.42 and 1.90; the generalized hyperbolic with a drawn skew, 1.40; the quantile noise, 1.42; the skewed horseshoe at seed 20261108, 3.03, 2.42 and 3.07 at b = 1, ¾ and ½. A rung at b = 0 with the step leaves every block on its prior, √(n/2)·D at most 1.63.
- *The main session's added gate, the skewed horseshoe at b = ½* (seed 20261107, slice 5's failing reading): with the step it fails, 4.96 at 1× and 10.12 at four times the size; at seed 20261108 it passes at 3.07. It fails on that seed at b = 1 and ¾ too (4.10 and 5.11), and without the step on the same design, σ moved only by its own block, at 4.27; the step alone passes there. **So interweaving does not cure the reading, and slice 5's joint move along β·v stays a follow-up with these readings as its measurement.**
- *R̂ and effective sizes per unit of work* (4 noises × 2 seeds, 4 chains × 8,000 draws): bulk ESS per second interwoven over plain, ln σ² 0.04 to 1.15, β 0.40 to 0.93, mean ln v 0.14 to 5.84 (three- to sixfold under the skewed generalized hyperbolic and horseshoe), M 0.11 to 3.75; under the quantile noise σ²'s ESS falls about twentyfold. R̂ moves with the ESS: the skewed horseshoe's σ² from 1.51 to 1.20 and 1.70 to 1.25, the quantile noise's from 1.01 to 1.04.
- *Refusals and acceptance:* the step moved every time it ran (8,000 of 8,000 per chain); no Geweke case refused a step for F.
- *The cost:* a fit with the step takes 1.3 to 2.5 times the work of one without.
- *The run:* a tempered fit with the step and a reference rung resumes bit for bit.

**Defects found and fixed in this slice.**
- *The GIG's prior density recomputed its Bessel normalizer on every read* (`GigMixing.LogPrior` through `LogPdf`), N reads per evaluation: an interwoven generalized-hyperbolic fit took 4.7 times the plain fit's work. It now reads the law's cached normalizer, 1.8 times. Found by the ESS row before the commit.
- *The Geweke harness read t's interwoven sweep with the sweep's own σ²,* which the step moves without keeping under a pure scale mixture: |z| 11 to 14 on the standardized residual. It now draws the coefficients and σ² afresh given the structure and the new latents. Harness only.

**Decisions taken beyond the brief.**
- *σ given z by one slice step on ln σ for every law and prior.* The derivation's exact GIG holds for a prior on σ in the GIG family with the coefficients' variance scaling as σ or not at all; BASS's priors are equivariant (the g-prior's variance scales with σ², σ²'s inverse gamma sits on w), so the slice step is the route, as the derivation itself names for the equivariant prior.
- *Opt-in, off by default,* so every fit's draws stand; given what it buys, it stays so.
- *Under a pure scale mixture the step moves v only:* σ² is the sweep's draw with the coefficients and is integrated again by the next structure step.
- *A held σ² refuses the switch,* which has no scale to move.

**Left unbuilt, and why.** The exact GIG route (no BASS prior takes it). E1's brief is built through slice 7: what stays open is the flat histogram under a noise (the user's decision, slice 5's report) and the follow-ups below and in slices 4 to 6.

**Follow-ups, with the decision each would change.**
- *z′ = σ²·v under a pure scale mixture.* Decision: whether interweaving earns a default for t and the other scale mixtures. Under z′ the likelihood is free of σ, so σ given z′ reads only its priors and moves along the ridge σ²·v holds, where z = σ·v reads every observation's variance and moves σ little. Measured by this slice's ESS-per-work row: if σ²'s ESS per second rises above the plain sweep's on t and the quantile noise's symmetric case, the switch defaults on for pure scale mixtures; if not, interweaving stays opt-in.
- *Slice 5's joint move along β·v,* with this slice's readings: the skewed horseshoe fails on seed 20261107 with either block moving σ.

### Slice 6: BARS under the noise, skewed and quantile BARS (chip, 2026-11-07)

**Landed** in `957629b` (the code and its tests) and the commit that carries this report (`docs/Maths`, the checklist). The user's go-ahead for slices 6 and 7 came through the main session on 2026-11-07 ("Yes, start both"). [The page's section](../../../docs/Maths/VarianceMeanMixtures.md#bars-under-the-noise-e1s-slice-6-oq20) holds σ²'s route under a skew against its alternatives, the priors read once, and the readers.

| Piece | Where |
|---|---|
| `StudentTWeights` through the law | `StudentTMixing.SamplePrecision`, the precision draw both take; `StudentTWeights` holds the law |
| BARS's noise | `BarsNoise` (`Quantile`, `Skewed`, `HeldSkew`, `StudentT`, `ScaleMixture`), on `BarsConfig.Noise` and the ensemble |
| The noise block | `BarsEnsemble`'s noise block on the weight block's cadence, `NoiseState` on the chain, its restart state and fingerprint |
| The Normal model under a noise | `WeightedNormalModel.ReadingFrom` (the priors read once), `Given` (σ² held), `SampleNoiseFit` (the coefficients, f and the prior's energy) |
| The shared conditionals | `VarianceMeanConditionals.SampleInverseScale` and `SampleSkew`, which BASS's σ² and shape blocks now call |
| The conditional variance at a given σ² | `GPriorConditional.ConditionalVarianceGiven` |
| A chain for tests and cost readings | `BarsEnsemble.Probe` (Geweke's successive simulator), `ChainOn` (one chain on the caller's thread) |

**Slice 6's gates.**
- *Gaussian untouched, and BASS's noises:* the whole regression assembly on `957629b`, 2,189 of 2,189 in 15 minutes, the eight pinned hashes and every fingerprint unmoved. Before the commit, the 624 tests the change reaches (BARS's, BASS's noise classes, the g-prior's readers and `BassBitIdentityTests`' eight pinned hashes) and the GBASS oracles' 27 passed; BASS's σ² and skew blocks call the shared conditionals with their arithmetic unchanged.
- *`StudentTWeights`' pinned values:* its weights equal the step as written before, number for number (ν = 0.5, 3, 4, 30; 500 residuals each); its marginal is the law's, which the BARS suites' checks to 9 and 12 places still meet.
- *The law, Geweke's test* on 16 points with priors that read no response (cubic B-splines, Poisson(1) knots, the centre at zero, σ² ~ IG(4, 0.6)), 400,000 against 600,000, eight functions: t(5) worst |z| 1.53, AL(0.8) 1.66, skewed GIG 2.41, skewed GBP 3.80 at seed 20261108, 1.29 at 20261109 and 0.85 at four times the size.
- *The given marginal* against a dense normal density of the working response: within 6·10⁻¹⁴.
- *Recovery,* the brief's and OQ20's: three quantile curves under heteroscedastic exponential noise, 10.5%, 50.2% and 89.5% of the data at or below the q = 0.1, 0.5 and 0.9 fits, RMSE 0.023, 0.051 and 0.162 against the true curves; a known skew of 0.8 covered by the 90% interval in three of four data sets, posterior means 0.635 to 0.865.
- *Refusals and acceptance:* no noise block refused in any Geweke case (BARS's support reads the unweighted Gram, so a weight cannot leave it); the ensemble refuses a noise with another model, with a resampler, tempered or on an adaptive kernel, each with its reason.
- *The cost row* (thread cycles, loaded, indicative; one chain on the test thread): a step at n = 1,000 costs 257 to 313 µs under a noise against 211 plain and 291 for `StudentTWeights`; at n = 200 70 to 139 µs, the load and the knot count swinging more than the noise.
- *The run:* quantile and skewed fits resume bit for bit, the skew and the weights included; the pointwise density is the law's at each draw's σ² and skew.
- *"The same oracle modes":* BARS's σ² and β steps are now the code GBASS's conditional oracle checks through BASS (`VarianceMeanConditionals`), and its v step is the shared law's; the oracle's 27 pass.

**Defects found and fixed in this slice.**
- *The first cost harness read nothing true.* It timed whole runs on the test thread's cycle clock while the chains ran on the pool (robust BARS read 0.6 µs a step), and then on the process's CPU time, which the pool's spinning inflates (plain BARS read slower at n = 200 than at 1,000). The row now steps one chain on the test thread.
- *A snapshot build kept a stale assembly:* `git archive` stamps files with the commit's time, older than the previous build's outputs, so MSBuild skipped recompiling the source under KHAOS's new tests. A fresh snapshot now drops its build output. Harness only; no commit carried it.

**Decisions taken beyond the brief.**
- *σ² is held on BARS's chain under a skew,* with the knot step reading the marginal given σ²: σ sits in the offset, so it cannot stay integrated. The alternative, the drift in the response's units, would have kept σ² integrated but cannot carry the quantile noise, whose skew is fixed in σ's units.
- *BARS's priors under a noise read the data once* (`ReadingFrom`), as BASS's do under a moving working response; a run without a noise reads as it did.
- *The σ² and β conditionals are shared code* (`VarianceMeanConditionals`), so the two engines cannot drift apart.
- *`BarsNoise` beside `IWeightResampler`,* not in it: the resampler's seam has no offset, no σ² state and no current v for a stepped law; the contaminated normal and `StudentTWeights` keep it.
- *No tempering, no adaptive kernel, no drawn γ under BARS's noise,* as its resampler has none; each is refused with its reason.
- *The ensemble's run context is set in one place* (`SetContext`), shared by a run, the probe and `ChainOn`.

**Left unbuilt, and why.** Tempered BARS under a noise: BARS's weighted marginal is built at b = 1 only, so robust BARS has no tempered form either, and a tempered noise would start there. The Normal-Wald's drawn γ on BARS: `BarsNoise` holds its mixing law. Both are recorded as open, with no decision waiting on them. Slice 7 follows.

**Follow-ups, with the decision each would change.** None from this slice.

### Slice 5: the engine's reach, ladders, SMC and a direct functional response (chip, 2026-11-07)

**Landed** in `37137f2` (the code and its tests) and the commit that carries this report (`docs/Maths`, the notation ledger). [The page's section](../../../docs/Maths/VarianceMeanMixtures.md#the-engines-reach-ladders-smc-and-a-direct-functional-response-e1s-slice-5) holds the tempered target, the prior's draw with the latents, the evidence and the flat histogram's refusal, each against its alternatives.

| Piece | Where |
|---|---|
| Every noise block at the rung's power b | `BassCoefficientsBlock.Power`, which `BassCoefficientsBlock`, `BassLatentScaleBlock`, `BassScaleBlock` and `BassNoiseShapeBlock` read (`BassNoise.cs`) |
| The latents from their prior | `VarianceMeanMixture.PriorLatents` |
| The prior's draw with the latents, the b = 0 member | `BassPriorSampler.DrawPrior` (the reference sampler's `Draw`), `LastDesign`, `Holds`; `DrawState`, a chain's start, on the model's working response |
| SMC's particles under a noise | `BassParticles.Draw`: the particle adopts the sampler's design (`BassSmc.cs`) |
| The ladder's hand-over of the reference's draw | `BassLadder.HandOver`: under a noise the replica adopts the sampler's design |
| A latent per observation under a direct functional design | `VarianceMeanMixture.Working`, its refusal lifted |
| The config's rules | `BassConfig.Validate`: the likelihood path in the integrated form for ladders and SMC, the flat histogram refused, a Jeffreys σ² under a skew refused where b = 0 draws it |

**Slice 5's gates.**
- *Gaussian untouched:* the whole regression assembly 2,121 of 2,121 on `e43f209` with this slice and 2,160 of 2,160 on `b408bf4` (KHAOS slice 1) with it, the eight pinned hashes and every fingerprint unmoved; the distributions' 347 of 347; the samplers' 285 of 285; the GBASS oracles 27 of 27. A Gaussian path takes nothing new: a Gaussian ladder reloads the reference's draw as before, and `DrawState`'s load on the model's working response is its own design's.
- *The law at each rung, Geweke's test on the tilted model* (π_b is the posterior of a model whose noise variance is s²/b, whose σ² prior is IG(α₀ − N(1 − b)/2, β₀) and whose v law's exponent is raised by (1 − b)/2; 400,000 marginal-conditional draws against 600,000 sweeps at b, N = 12, nine functions): seven cases, every |z| under 4, the worst per case t(5) at b = 0.4 2.24 and at 0.75 2.28, the generalized hyperbolic with a drawn skew at 0.4 2.06 and at 0.75 1.88, GBP(1, ½, 1.5) with a drawn skew at ½ 2.23, the skewed horseshoe at ¾ 2.62 and the symmetric horseshoe at ½ 2.17.
- *The skewed horseshoe at b = ½ fails it,* recorded as a reading: |z| 5.33 (seed 20261107, mean ln v), 5.41 (seed 20261108, β), 9.16 at four times the size (mean ln v), the successive chain's ln v low each time. Its tilted v law, GBP(1, ¾, ¼), falls as v^(−5/4); under a drift a large vᵢ makes rᵢ ≈ σβvᵢ, so β and vᵢ pin each other along a ridge spanning decades of v. Beside it everything passes: the symmetric horseshoe at b = ½ (worst 2.17 and 2.64 on two seeds, 0.72 at four times the size), the skewed horseshoe at ¾ and the lighter GBP at ½ in the gate, and v's step under a drift at power ½ from exact draws (√n·D 0.81 and 1.30, cut 2.6). So the step's conditional is right and the scan crawls along the ridge; the follow-up below names the move that would test it.
- *b = 0 is the whole prior:* 3,000 prior draws of a skewed Normal-Wald's latents (β, γ drawn) against their priors, √n·D 0.59 (β), 0.83 (γ), 1.27 (σ²), 0.62 (v₀ through its inverse Gaussian CDF), in 3,352 attempts; SMC's kernel at b = 0 on exact prior draws against fresh ones, √(n/2)·D at most 1.20 over M, ln τ, ln v₀, mean ln v, β and ln σ², under t(4), the skewed Normal-Wald and the skewed GBP.
- *Evidence where the values are complete,* against plain Monte Carlo over 2,000,000 prior draws (N = 12; the ladder at 10 levels, hot end 0.02 and a reference rung, over 4 chains; SMC at 2,000 particles over 4 replicates): t(5), ln p(y) −17.239 ± 0.016 against the ladder's −17.247 ± 0.011 and SMC's −17.248 ± 0.010 (z 0.43, 0.47); the generalized hyperbolic with a drawn skew, −18.050 ± 0.011 against −18.026 ± 0.012 and −18.013 ± 0.007 (z 1.48, 2.90); the bar 4.
- *A per-observation latent under a direct functional design:* 60 cases × 20 locations under 0.85·N(0, 0.2²) + 0.15·N(0, 1.5²), two data sets, 4 chains × 500 draws: t(3)'s 90% intervals for f cover 92.6% of 840 case-location points at an RMSE of 0.044, Gaussian BASS's 96.1% at 0.074. The unit tests: N·M latents a draw, the pointwise unit the curve, its observations' law densities summed to 10 digits.
- *Categorical predictors* run a skewed Normal-Wald with nothing new.
- *Refusals and acceptance* (decision 4): no structure move refused in any Geweke case; the latent block refused 5 of 600,000 draws under the symmetric horseshoe at b = ½ and none elsewhere, each outside F and kept as the gate's exact step; the prior's draws with their latents took some 2,300,000 attempts for 2,000,000 in the evidence gate, the rest outside F or drawing an intercept term.
- *The cost row:* (thread cycles, loaded, indicative; whole sweeps after an untimed warm-up pass, M in brackets) t(4) at the base 2,092 µs at b = 1 (M 20), 1,310 at ½ and 240 at 0.05 (M 0), against Gaussian BASS's 32, 30 and 28; at n 200 190, 81 and 84 against 16, 7.4 and 16; on the direct functional point (100 cases × 20 locations, the stacked design under per-row weights) 2,722, 2,787 and 496 against 30, 32 and 21, its Gaussian readings swinging to 100 between runs. Slice 1's per-part harness on this code reads its base sweep at 962 µs against the 961 it documented.
- *The run:* a tempered fit under a skewed Normal-Wald with a reference rung and the evidence resumes bit for bit, v, σ², β and γ and ln p(y) included; tempered fits under t, the horseshoe and the skewed Normal-Wald, and SMC fits under t and the skewed Normal-Wald, report ln p(y), their ln σ² within 0.25 of the chains'.
- *`docs/Maths`:* [the page's section](../../../docs/Maths/VarianceMeanMixtures.md#the-engines-reach-ladders-smc-and-a-direct-functional-response-e1s-slice-5); the ledger takes the rung's b beside the GIG's and the GBP's.

**Defects found and fixed in this slice.**
- *A resumed ladder with a reference rung refused the reference's first draw under a noise.* The restored replicas' designs read the response as the checkpointed model had; the new reference sampler built its latents on the new model's reading; a design takes only working responses of its own reading, so the replica's load threw. The hand-over above fixes it, and the resume test pins it.
- *`DrawState` after a `DrawPrior`* would have loaded a chain's start on the latents the last prior draw left in the sampler's design. It now loads on the model's working response. Found while building, before any commit.

**Decisions taken beyond the brief.**
- *The flat histogram is refused under a noise.* The brief listed it in this slice; its walkers' π_w and WHAM read the marginal over a and σ², b·ℓ₁, whose weights couple every vᵢ, so v has no step of its own there. Two routes would carry it, each a sampler of its own (the page has both); the follow-up below asks which, if either.
- *Rungs temper the likelihood given the latents, on the likelihood path in the integrated form only;* the linear form and the powered path are refused for the same coupling. The integrated form is the configs' default, and its b = 0 member is the whole prior.
- *A chain's start keeps v at 1* (`DrawState` unchanged): only a reference rung and SMC draw the latents from their prior, so every MCMC fit's draws stand as they were.
- *A Jeffreys σ² under a skew is refused with a reference rung or SMC,* whose b = 0 draws σ² from its prior.
- *The reference's draw is adopted with its design under a noise only.* A Gaussian ladder reloads it as before, its pinned bits untouched. Under a noise the evidence gates read the same digits before and after the change, so adopting and reloading agree there.
- *The reach cost row warms the JIT* with an untimed first pass: its first cell read 3,157 µs cold against about 2,000 warm (a probe: `BassChain` 4,024 cold, 2,024 warm; a bare chain at b = 1, 2,108), and slice 1's per-part harness, rerun on this code, reads 962 µs at the base against its documented 961, so the noise's sweep did not move. A bare chain at b = 1 and an MCMC chain cost the same: a ladder's cold rung pays nothing for its temperature.

**Left unbuilt, and why.** The flat histogram under a noise (above). Slices 6 (BARS through the shared law) and 7 (interweaving) wait for the user, beside KHAOS's build. A direct functional fit's candidate density is the conformal layer's and refuses every direct functional fit, Gaussian included (its candidate is a scalar observation; whole curves are E17's), so nothing here waits on it. `docs/Maths/BayesianRegression.md`'s row for GBMARS in what each extension still writes is stale for slices 4 and 5; that file held KHAOS's uncommitted edits, so its update went to that session to carry.

**Follow-ups, with the decision each would change.**
- *The flat histogram under a noise.* Decision: whether a noise runs on the flat histogram, and by which route. (a) The walker holds a and σ² beside v, so ℓ is the conditional likelihood, a sum whose change per vᵢ is O(1); a and σ² then need Metropolis steps under w(ℓ) that the collapsed structure moves do not, and the multicanonical weights span a wider range of ℓ. (b) Each vᵢ moves by Metropolis under w(ℓ₁), the marginal refitted by a rank-one weight update, O(M²) a vᵢ and O(N·M²) a sweep, the structure moves as they are. (c) Leave it refused: ladders and SMC carry every noise, with the evidence.
- *A joint move along the drift's ridge.* Decision: whether a skewed heavy-tailed noise on a ladder's hot rungs needs a move beyond the one-at-a-time scan. A Metropolis move of β and every vᵢ along β·vᵢ held (β → β/κ, v → κv, the weights refitted once) would be measured by the skewed horseshoe's Geweke at b = ½ and the swap rates of a ladder under it: if it passes, the move joins the latent slot; if not, the reading stands as that noise's mixing limit below b ≈ ¾.

### Slice 4: the generalized beta prime and the horseshoe (chip, 2026-11-06)

**Landed** in `4662639` (the code and its tests, the readout oracle's corrected `hbass`) and the commit that carries this report (`docs/Maths`, the notation ledger). [The page's section](../../../docs/Maths/VarianceMeanMixtures.md#the-generalized-beta-prime-and-the-horseshoe-e1s-slice-4) holds the choices with their alternatives; [the distributions'](../../../docs/Maths/ProbabilityAndInformation.md#the-generalized-beta-prime-and-its-normal-mixture-e1) the law and the mixture's integral.

| Piece | Where |
|---|---|
| The GBP | `GeneralizedBetaPrime` (`src/maths/distributions/GeneralizedBetaPrime.cs`); its quantile through `Beta.InverseLogit` |
| The normal mixture over any log-concave mixing law | `NormalVarianceMeanMixture`, the generalized hyperbolic's integral and table made a base; `NormalGeneralizedBetaPrime`; `GeneralizedHyperbolic` on it |
| GBP mixing, v's slice step | `GbpParameters`, `GbpMixing`; `IMixingLaw.StepConditional`; `SliceSampler`'s struct overload (`IScalarLogDensity`) |
| The kinds on the config | `BassNoise.GeneralizedBetaPrime`, `Horseshoe`, `StrawdermanBerger`; `BassNoiseKind.GeneralizedBetaPrime`, `BassNoise.BetaPrimePrior` |
| The latent block from the current v | `BassLatentScaleBlock` |
| The pointwise density, opt-in | `IBassNoiseLaw.DensityByQuadrature`, `IMixingLaw.MarginalByQuadrature`; `BassPointwise.LogDensity`'s and `BassFit.PointwiseLogDensity`'s `quadrature`; the conformal layer's candidate density the same (`BassCandidates.LogDensity`, `BassFit.CandidateLogDensity`, in its own commit after that layer's seam landed) |
| `hbass`'s half-Cauchy | Nothing new: `Sigma2Prior.HalfCauchy()` and `BassSigma2RateBlock`, run under the noise |
| `hbass` corrected in the oracle | `gbass_readouts_oracle.R`'s `prop_sigma` for `hbass`, v's posterior and step emitted |

**Slice 4's gates.**
- *Gaussian untouched:* the eight pinned hashes and every fingerprint unmoved, the whole regression assembly 2,060 of 2,060, the distributions' 347 of 347, the samplers' 285 of 285 (the slice sampler's struct form draws the delegate form's numbers), the GBASS oracles 27 of 27; the generalized hyperbolic's tails and tables hold their slice 3 bars.
- *The law, the step:* the slice step from exact draws of v's conditional (rejection from the prior) leaves them on it, √n·D at most 1.17 against the cut 2.6 at 20,000 draws, on eight cases (the horseshoe at r = 10⁻³, ½ and 3, with a drift, at power ½ and at the prior, the Strawderman–Berger with a drift, GBP(2, 1.5, 3)), every draw moved; a 400,000-step chain at r = 0.01 carries the conditional's mean of ln v, −4.582 ± 0.005 against −4.575. The exact laws' steps equal their draws, number for number.
- *The law, Geweke's test* at slice 3's size (400,000 against 600,000): the horseshoe with σ² integrated on two seeds, with β drawn and σ² on the record, and the Strawderman–Berger, every |z| under 4 over nine functions; the worst, 3.67 (ln τ, the horseshoe at seed 20261106), reads −0.19 at twice the size; three more seeds read ≤ 1.91, and four times the size on a fifth ≤ 2.64 (the skewed form's v functions, whose sign changes between seeds). `hbass`'s half-Cauchy has no Geweke case: its scale reads y.
- *Primitives:* the GBP against the half-Cauchy on √v to 10⁻¹⁴ from 10⁻¹² to 10¹²; its tails and levels against Simpson's rule to 6·10⁻¹¹ on seven laws; a seeded sweep of 400 laws over far ranges (quantiles back through the tails from 10⁻²⁶¹, the complement to 10⁻¹²); draws by Kolmogorov–Smirnov on 65 laws. The mixture's density by the integral against the generalized hyperbolic's closed form to 7·10⁻¹⁵ and the horseshoe's exponential integral to 10⁻¹⁵ from |x| = 10⁻⁸·s to 10⁸·s; density and tails against Simpson's rule to 3·10⁻¹³ on six laws, poles and drifts included; the table against the integral out to the 10⁻⁹ quantiles within 2.5·10⁻⁴, 10⁻⁵ in the body. The distributions' whole assembly 347 of 347; the samplers' 285 of 285.
- *Recovery* from a contaminated normal (N = 300, six data sets): the horseshoe under `hbass`'s half-Cauchy covers the true surface at 80.6% of 294 grid points for 90% intervals, RMSE 0.086, against Gaussian BASS's 86.4% and 0.142 and slice 1's t(3) 89.5% and 0.059; held-out responses 92.9% against 92.7%. Its v near 0 lets a fit follow some observations closely; its intervals for f undercover.
- *`hbass`'s readouts* (N = 200, R 6,000 iterations), against `hbass` with v's step held at 2.5 (the defect below): at engine seed 41 the surface within 0.55 of the wider 90% half-band under the g-prior and 0.82 under the ridge; σ² 0.0765 and 0.0802 against w·var(y) 0.0685 [0.043, 0.103]; mean ln v −0.177 and −0.192 against −0.162 and v below 0.1 at 16.7% against 17.2%; M 5.84 and 5.10 against 5.36. Over seeds 41 to 46 the ridge read 0.52 to 0.82 and the g-prior 0.55 to 1.20, so the bars are `qbass`'s, the ridge at parity 1.0 and the g-prior at 1.5; against `hbass` as written the g-prior read 0.97 to 1.56.
- *Refusals and acceptance* (decision 4): v's slice step always moves (every one of 20,000 draws in each of eight cases); the latent block refused 0 to 9 of 600,000 draws per Geweke case, each outside F; the unit tests' chains refused none in 150 sweeps.
- *The cost row* (thread cycles, loaded, indicative): a horseshoe sweep is 767 µs at the base under the g-form against a Gaussian sweep's 49 (16×), 110 against 4.2 at n 200 and 319 ms against 1.6 at M 200; v's slice steps 413 µs at the base, 6.8 evaluations and 0.41 µs a step, 5.6 times t's exact draw; the pointwise read 41 µs (wall) against t's 0.11, 82 s at N 1,000 and 2,000 draws; a predictive's table per draw 15 to 28 ms (wall), most of the recovery's 29 minutes.
- *The run:* a horseshoe fit with a drawn skew under the half-Cauchy resumes bit for bit, v, σ², β and σ²'s rate included; a retained draw's view rebuilds bit for bit under the horseshoe, its skewed form and the Strawderman–Berger; the sweep runs σ²'s rate's block beside the noise's.
- *`docs/Maths`:* the page's section, with v's step, the pointwise density and the mixture's upgrades each against their alternatives; `ProbabilityAndInformation.md` takes the GBP, the shared mixture and GBP mixing; the ledger takes the GBP's p, a, b beside the GIG's and `hbass`'s c and w.

**Defects.**
- *GBASS's GBP v step, worse than slice 0 read it.* Slice 0 found `rgbp`'s step adapted from mean(ln v)'s spread (0.543 at N = 60, acceptance 0.90) and called the chain valid but slow. At N = 200 the step adapts to 0.12 (acceptance 0.98) on the readouts' data and 0.038 (0.99) on a second set, so in 6,000 iterations v barely leaves its start at 1: `hbass` reads M 4.16 where the same sampler with its step held at 2.5 reads 5.36, and v below 0.1 in 15.9% of draws against 17.2%. The readouts compare against that correction, not against the frozen chain.
- *Found and fixed in the mixture's integral* (the generalized hyperbolic's, now shared):
  - Far out on the drift's side z's two terms cancel, and once the peak is narrower than u's rounding the integrand is noise Kronrod's tolerance cannot meet: a skewed horseshoe's table took 107 s, a generalized hyperbolic's with a strong drift 20 s. Fixed by taking the integral about its peak.
  - The density's narrow peak was stepped over by the scan and missed by the panels' clamp at 10⁻⁶. Fixed by Newton about the peak and a golden section over the span.
  - The table's cubic in x missed a heavy tail's integral by up to half the log tail at 10⁻⁹. Fixed by the cubic in asinh(x/s).
  - A pole at the median gave the table a point of infinite slope, read as NaN. Fixed by the pole's steps.
- *The GBP's quantile through the smaller tail's y/(1 − y)* lost up to 3·10⁻³ where y sat near 1; it now takes the beta's logit.

**Decisions taken beyond the brief.**
- *`IMixingLaw` steps rather than draws:* `StepConditional(current, …)` replaces `SampleConditional` on the seam, an exact law answering with its draw (the classes keep `SampleConditional`); the alternative, a draw that throws for GBP, would have made the seam lie for one family.
- *The mixture's integral shared and upgraded,* moving the generalized hyperbolic's tables by up to 3·10⁻⁴ of a log tail in the far tails (its draws and densities unmoved); the alternative, a second integral for the GBP, would have left the generalized hyperbolic's far-drift tables at 20 s and its heavy tails' cubic in x.
- *`SliceSampler` takes a struct density,* the delegate form calling it with its numbers unchanged: N steps a sweep allocate no closures.
- *The pointwise latent read as integrated in-sample* (`LatentTreatment.Integrated`, at the conformal chip's ask, which edited `BassStudentTTests` to expect it): given the globals nothing outside its observation moves a vᵢ, so WAIC and PSIS-LOO read the term. And `BassPointwise.Latents`, the report extracted for its `BassCandidates`.
- *The readouts' bars are `qbass`'s:* the ridge at parity ≤ 1.0, the g-prior at 1.5.

**Left unbuilt, and why.** The engine's reach (5), BARS (6) and interweaving (7) are later slices. The horseshoe's closed density at m = 0, the auxiliary-gamma Gibbs for p = 1 and a Geweke check of the half-Cauchy under the noise wait on the follow-ups below.

**Follow-ups, with the decision each would change.**
- *One table for every draw whose noise has the same shape.* Decision: whether a predictive under a mixture noise builds a table per draw. Without a drawn skew or γ the noise at each draw is one law scaled by √(σ²c) (m/s = β/√c is fixed), so one standardized table, read at x/s, would serve them all: today the horseshoe's predictive quantiles over 2,000 draws build 2,000 tables, most of the recovery's 29 minutes. The table's reads move by rounding.
- *The horseshoe's closed density at m = 0.* Decision: whether the symmetric horseshoe's pointwise density needs the opt-in. (2π³)^(−½)·e^c·E₁(c)/s is one exponential integral, against 41 µs of quadrature; a second path for the case `hbass` declares.
- *Auxiliary-gamma Gibbs for p = 1.* Decision: v's step under the horseshoe. v | λ ~ Gamma(a, λ), λ ~ Gamma(b, 1) makes v's conditional a GIG, drawn by ratio-of-uniforms in some 0.15 µs, and λ | v a gamma: exact draws, perhaps half the slice step's 0.41 µs, at a second latent per observation to retain, and Gibbs between v and λ that slows where they correlate. A measurement of the effective draws per µs against the slice step would decide it.
- *A stated scale for the hierarchical σ² prior under a noise.* Decision: whether `hbass`'s half-Cauchy gets a joint-law gate. Its scale reads y through the reference, which Geweke's successive simulator moves; a reference stated on the config would let the test hold it, and the readouts carry it until then.

### The GIG's ratio-of-uniforms draw (chip, 2026-11-06)

**Landed** in `72fb8e2` (the code, its tests and the qbass gate's bars) and the commit that carries this report (`docs/Maths`). At the user's word that spline surfaces are expensive and every bit of performance counts, after slice 3. [The page's section](../../../docs/Maths/VarianceMeanMixtures.md#the-gigs-draw-ratio-of-uniforms-where-it-is-faster) holds the measurement, the rule and the trade-off.

| Piece | Where |
|---|---|
| Ratio-of-uniforms, with and without a shift by the mode | `RatioOfUniforms.Draw`; `RatioOfUniforms.Bounds`, the rectangle by search for a family without closed sides (`src/maths/distributions/RatioOfUniforms.cs`) |
| Hörmann and Leydold's GIG draw | `GeneralizedInverseGaussian.SampleByRatio`: the shifted ratio where λ > 2 or ω > 3, the unshifted where T₋½-concave or ω > 0.2, a three-piece hat below both |
| Which one runs | `GeneralizedInverseGaussian.DrawsByRatio(p, a, b)` |

**Measured** (thread cycles, indicative; µs a draw with its own setup, the least of five passes after six warm-up passes): on v's GIG parameters recorded from burned quantile and Normal-Wald chains at the cost row's points, ratio-of-uniforms 0.12 to 0.21 against the tangent hat's 0.52 to 0.95, 2.5 to 5.6 times faster; over 132 far-range points, λ to 1,000 and ω from 10⁻⁸ to 10⁸, 0.08 to 0.36 against 0.28 to 0.92, faster at every point. Each sampler's share of a sweep, each row run alone: on the quantile fit at the base v's draws went from 952 µs to 208 (53% of the sweep to 24%) and the sweep from 1,810 to 850; on the Normal-Wald v's block from 987 to 681 and the sweep from 1,271 to 1,027, at n 200 from 235 to 133. At M 200 the Follow after v dominates and nothing moved.

**Which one runs.** Ratio-of-uniforms for every untruncated draw with a, b > 0, |p| ≤ 10³ and 10⁻⁸ ≤ ω ≤ 10⁸, the range its law was checked over; the tangent hat beyond it and for every truncated draw, where the ratio's rectangle would reject nearly everything in a far tail; the gamma and inverse gamma at the limits. A function of (p, a, b) alone, so a resumed run draws as the whole one did.

**The gates.**
- *The law against quadrature at the deciles and by Kolmogorov–Smirnov:* 18 chosen cases across the three regions and 300 seeded far-range cases (|p| to 10³, a and b over 16 decades), worst √n·D 1.76 against the cut 2.6; at 1,000,000 draws on nine of the fits' own parameter sets, means, variances and deciles within sampling error (the one case at mean z −2.48 read 1.16 and −0.25 at 4,000,000 draws on two other seeds). The primitive finds N(0, 1)'s rectangle at ±√(2/e) to 7 digits and draws it and a gamma's by search on their laws.
- *Against `gbass_primitives_oracle.R`:* the default draw, now this one, carries R's GIG law at its deciles (`EngineSamplers_CarryTheOraclesLaw`, 26 of 26 GBASS oracle tests).
- *Gaussian BASS untouched:* it draws no GIG; the whole regression assembly 2,022 of 2,022, the distributions' 317 of 317, the samplers' 285 of 285. t without a skew draws v by its own gamma step and did not move.
- *The law under the sweep:* slices 2's and 3's Geweke cases pass with the new draw (23 of 23 qualification tests).
- *What moved, and why.* Every quantile, Normal-Wald, skewed-t and generalized-hyperbolic fit draws v anew, so its draws and readings moved; no pinned hash covers those fits. The quantile recovery's held-out shares read 0.517 and 0.888 (0.518 and 0.884 before), the skew's sign 0.898 and 0.103 (0.903 and 0.106), `qbass`'s σ² 0.00723 [0.00567, 0.00909] under the g-prior (0.00768 before) and 0.00803 under the ridge, `nwbass`'s noise variance 0.1825 and 0.1804 against 0.1821. The envelope's pin moved by exactly its 200 untruncated GIG cases and was rebaselined, ce13801abaff15d3 to 624f8b056215181a.
- *`qbass`'s surface gate, re-set* (approved 2026-11-06). Its g-prior arm read 1.26 of the wider half-band at seed 43 against the 1.0 bar. Over engine seeds 43 to 48 that arm read 0.49 to 1.23 (mean 0.85) under the tangent hat and 0.61 to 1.26 (mean 0.92) under the new draw, two of six above 1.0 under either: the bar measured seed 43's stream more than the sampler. The g-prior arm carries a stated model difference from `qbass`, so its gap is not expected to be zero; it now takes worst ≤ 1.5, a bar for gross defects, not parity. The ridge arm, GBASS's own form, keeps the parity bar ≤ 1.0 and read 0.51 to 0.81 on all six seeds.

**Decisions taken beyond the brief.** The ratio's three regions take GIGrvg's thresholds, read from Hörmann and Leydold's algorithms, as GIGrvg 0.8's C source is not installed; its law was checked at those thresholds, each region on chosen cases and the sweep. Its range is the one checked, not GIGrvg's whole domain.

**Left unbuilt.** No other family takes the primitive yet: none has a draw it would speed; the gamma's and the inverse gamma's untruncated draws are their own.

### Slice 3: the Normal-Wald, the skewed t and the generalized hyperbolic, with the envelope's guard (chip, 2026-11-06)

**Landed** in `f2e1e58` (the envelope's guard, its own commit as asked), `5e837c4` (slice 3's code and tests) and the commit that carries this report (`docs/Maths`, the notation ledger). [The page's section](../../../docs/Maths/VarianceMeanMixtures.md#the-normal-wald-the-skewed-t-and-the-generalized-hyperbolic-e1s-slice-3) holds the model, the alternatives and the checks; [the distributions' section](../../../docs/Maths/ProbabilityAndInformation.md#log-bessel-k-the-gigs-normalizer-and-the-generalized-hyperbolic-e1) the log Bessel K, the GIG's normalizer and the generalized hyperbolic; this is the account.

| Piece | Where |
|---|---|
| The envelope's guard | `LogConcaveEnvelope.GuardedPoint`, `Guarded`; `LogConcaveEnvelopeTests` |
| Log Bessel K of real order | `SpecialFunctions.LogBesselK`, `LogBesselKScaled` (`src/maths/distributions/BesselK.cs`) |
| The GIG's normalizer, density and moments | `GeneralizedInverseGaussian.LogNormalizer`, `LogPdf`, `Mean`, `Variance` |
| The generalized hyperbolic | `GeneralizedHyperbolic` (density closed, tails by the integral over v, a reader's table) |
| General GIG mixing, the skewed t | `GigMixing`; `StudentTMixing`'s drift routed to it |
| The kinds on the config | `BassNoise.NormalWald`, `SkewedT`, `GeneralizedHyperbolic`; `BassNoisePrior` |
| β and γ on the record, their block | `ObservationLatents.Gamma`, `WithShape`; `BassNoiseShapeBlock` |
| The predictive's root by Newton | `PredictiveMixture.Quantile` (the law path) |
| `nw_est_mom` and `nw_gamma_prior`, corrected | `NormalWaldCalibration.Moments`, `GammaPrior`; `NormalInverseGaussian` |

**The guard.**
- *What it does.* Newton's search for the design point runs as it did; where it runs out of steps short of its point, a guarded search runs from the reach: Newton where h ≥ −2²⁰, and below it Newton on ln(−h) against ln t, exact for a power of t and close for an exponential, or a halving in ln t where h has overflowed. So every search that converged before is unchanged, and any reach is valid.
- *Every existing family's draws:* the eight pinned hashes, every fingerprint and each family's truncated-draw tests pass unchanged (the whole regression assembly 2,022 of 2,022, and the GBASS and distribution oracles 48 of 48, the distributions' 297 of 297, the samplers' 285 of 285). A new pin over 1,311 cases of the families' own draws, each case on its own stream, moved 31 cases and no family but one: all the modified half-normal's, 25 in its branch in x at α from just above 1 to 2.2, where the search met x = 0, and 6 in its ln x branch at α near 10⁻³, from a reach on the exponential side. Each is a case whose own reach ran Newton out of steps, the guard's target; no other draw moved.
- *A seeded sweep with bad reaches:* every family's own log ratio (normal, t, beta, gamma, GIG at and off the peak, the half-normal's three), 30 cases each over far ranges, drawn with reaches from the root of h = −1 to 10¹⁰⁰ times it and +∞, against its law by quadrature. Before the guard, from ×10 or ×1,000 on, the GIG's, gamma's and the half-normal's hats threw or accepted nearly nothing, every family threw at +∞, and at ×10¹⁰⁰ beta's and t's draws failed Kolmogorov–Smirnov silently, the hat sunk below h (acceptance read above 1). After it, every law holds at every reach.
- *Acceptance before and after* on those far ranges: the worst case per family at reaches ×10 to +∞ went from 0 (and no finite area) to 0.867 or more; at the root's own reach the half-normal's branch in x, which threw in 7 of 30 cases, draws, and the GIG's worst rose from 0.771 to 0.872; elsewhere it is unchanged.
- *Defects it fixed outside the sweep:* the GIG's own truncated draw on far upper tails at small |p|, a and b, which threw in 4 of 4,000 cases (`Gig_OnAFarUpperTail_DrawsItsLaw`).
- *The GIG's own reaches* (`ExponentialReach`, `LinearReach`) are no longer needed for correctness: with them replaced by the curvature's reach, slice 2's GIG law tests pass 17 of 17. They stay as first guesses, and save the guarded search on the ranges that need it, a draw four to seven times faster there (2.3 against 15 µs, 2.1 against 8.4).

**Slice 3's gates.**
- *The density against quadrature of the mixture:* GIG mixing's marginal against ∫ N(r | m·v, s²·v)p(v) dv taken over ln v with the GIG normalized by its own quadrature, 600 cases across the Normal-Wald, the skewed t, the variance gamma and general GIGs, drifts of both signs and none, to 10⁻⁸. The generalized hyperbolic's tails by the integral against its density's mass at its quantiles 0.001 to 0.999, to 1.3·10⁻⁹ on six laws.
- *Log Bessel K against `gbass_primitives_oracle.R`'s `besselK`:* 4.2·10⁻¹⁵ at its nine points; against its integral on 2,000 cases, ν to 10⁴ and x from 10⁻³⁰⁰ to 10³⁰⁰, 10⁻¹².
- *`nw_est_mom` corrected:* the NIG by its moments, in closed form with the NIG's kurtosis, recovers three NIGs from their moments taken by quadrature of their own densities to 7 digits.
- *`nw_gamma_prior` corrected,* with targets a normal can reach: the normal on γ > 0, solved exactly, puts its two masses where asked to 8 digits; at P(γ > g₁) = ½ the most it puts below g₂ is 1 − 2^(−g₂/g₁), so GBASS's defaults, asking 0.05 where that is 0.0016, are refused with that value; the defaults here are (0.5, 0.9, 0.5, 0.03). The prior on ln γ and the gamma prior are on the page with their trade-offs.
- *`nwbass`'s readouts* (β ~ N(0, 1), γ under N(2, 1), w's prior IG(0.1, 0.05), R at `lag_beta` = 0): the noise's variance 0.1814 under the g-prior and 0.1809 under the ridge against `nwbass`'s 0.1821; w, β and γ each inside the other's 90% interval (w 0.189 against 0.155, β 0.205 against 0.191, γ 1.13 against 0.92); the surface within 0.55 of the wider half-band at worst (0.76 under the ridge). w and γ trade along the noise's variance, which the data fix; the readouts oracle now reports GBASS's own `s2` for it. The stated differences: GBASS's untruncated γ draw and slice 2's list.
- *`lag_beta`, read:* eight chains each way from the empty structure on NIG data, free against held at 0 for 500 sweeps: the median chain reached the reference's median M in 23 sweeps free and 26 held; free, β ran at 1.68 over the first 200 sweeps against the reference's 1.21, without slowing the structure; over sweeps 1,500 to 2,000 the arms read M 5.06 and 5.16, β 1.64 and 1.40 and σ² 0.127 and 0.121, against 5, 1.21 and 0.105. By the rule stated beside it, the hold is not adopted.
- *The law, by Geweke's test* at twice slice 2's size (400,000 marginal-conditional draws against 600,000 sweeps): the Normal-Wald with β and γ drawn on two seeds, with γ alone (σ² integrated), the skewed t and the generalized hyperbolic with β drawn, every |z| ≤ 2.15 over ten functions each, β and ln γ among them. Its power: β's conditional with its data term's sign flipped fails every skew-drawn case (|z| 20.8, 34.1, 108.7 and 335.9) and leaves the case without a skew passing. At slice 2's size one function of one case reached |z| 4.08 (f₀, the Normal-Wald at seed 20261106); the same seed at twice the size gives −0.31.
- *Recovery* on NIG noise (N = 600, σ 0.3, β 1, γ 1.5): β 1.29 [0.57, 3.80], γ 1.82 [0.93, 3.89], σ² 0.113 [0.062, 0.226] at 99%, each covering the truth; the surface within RMSE 0.056 of the truth, the mean within 0.029.
- *Refusals and acceptance* (decision 4): no block draw refused in any Geweke case, 0 of 600,000 each.
- *The cost row* (thread cycles, loaded, indicative): a Normal-Wald sweep is 1,271 µs at the base under the g-form against a Gaussian sweep's 28 (45×), 235 against 4.2 at n 200 and 320 ms against 1.6 at M 200; v's block, its N GIG draws and the Follow after them, is 78% of it at the base and 72% at M 200, and β's and γ's block 26 µs, near σ²'s. The table is on the page.
- *Retention* is slice 2's: a retained draw keeps v with its scalars, γ now among them, rebuilt bit for bit under all three kinds; a Normal-Wald fit resumes bit for bit, β, γ and σ² included.
- *`docs/Maths`:* the page's section; `ProbabilityAndInformation.md` takes the guard, log Bessel K, the GIG's normalizer and the generalized hyperbolic; the ledger takes γ beside the half-normal's and the NIG's, and repoints four of its references this slice shifted.

**Defects.** *Corrected here* (slice 0's findings): GBASS's γ drawn from its γ > 0 piece's normal on the whole line, now drawn on γ > 0; `nw_est_mom`'s kurtosis; `nw_gamma_prior`'s stray −1, and its defaults out of reach. *Found and fixed:* the envelope's design point (above), which slice 2's GIG truncated draw met on far upper tails. *In this slice's own first forms:* a noise CDF from a mixture of normals at nodes of v, 10⁻⁵ in the body and wrong far out under a strong drift, replaced by the generalized hyperbolic's integral and table; the tail integral's Gauss–Kronrod recursing without end on a NaN integrand (0·∞ far out at no drift), made finite and NaN-proof; a table on an x grid stretched to the 10⁻¹² quantile, which left a heavy tail's body unresolved, replaced by the normal score walked out from the median.

**Decisions taken beyond the brief.**
- *The guard as a second search after Newton's,* not a change to Newton's own steps: a first form did that, and moved 38 of the pinned cases whose searches had recovered on their own.
- *σ² on the record wherever the skew can be nonzero,* drawn or held away from 0; integrated where it is held at 0, so a symmetric Normal-Wald is a pure scale mixture as t is.
- *β's and γ's priors stay out of the target,* entering their own block's conditional only.
- *The predictive's noise is the exact generalized hyperbolic, read through a table:* the integral costs 110 to 230 µs a tail, and a predictive quantile reads hundreds per draw; the table holds the log tails to 10⁻⁵ in the body and 3·10⁻⁴ at a power tail's 10⁻⁹, within what the Monte Carlo over draws resolves.
- *The predictive mixture's root by Newton* with the mixture's density, inside the bracket: a training fit's quantiles under the Normal-Wald went from minutes to a second. t's and the asymmetric Laplace's quantiles move by rounding; their readers' tests pass.
- *The NIG by moments in closed form,* no optimizer; *the prior on γ on γ > 0,* the conjugate truncated normal, with its reach stated.
- *The skewed t and the generalized hyperbolic exposed on the config,* both falling out of GIG mixing at no extra cost.
- *`lag_beta` not adopted,* by the read and its rule.

**Left unbuilt, and why.** Generalized beta prime (4), the engine's reach (5), BARS (6) and interweaving (7) are later slices; a prior on ln γ waits on the follow-up below.

**Follow-ups, with the decision each would change.**
- *A prior on ln γ for the Normal-Wald.* Decision: the default prior on γ. A normal on γ > 0 cannot hold a far median with mass near 0, so a fit meant to range from near-Gaussian to the heaviest tails needs targets past its reach; a prior on ln γ reaches them in closed form at the cost of a slice step for γ.
- *A faster normal log tail in the generalized hyperbolic's integral.* Decision: whether its readers keep the table or read the integral everywhere; today a draw's table costs some 60 tail integrals, about 10 ms, so a predictive over thousands of draws spends seconds building them.

### Slice 2: the asymmetric Laplace, quantile regression (chip, 2026-11-05)

**Landed** in `5beaa49` (the code and its tests) and the commit that carries this report (the docs, the notation ledger). At the user's request (2026-11-05) the noise's documentation moved to its own page, [Variance-mean mixture noise](../../../docs/Maths/VarianceMeanMixtures.md), which takes slice 1's section from `BayesianRegression.md` and this slice's; `BayesianRegression.md` keeps a pointer and is 65 lines shorter. The page holds the forms, the alternatives and the checks; this is the account.

| Piece | Where |
|---|---|
| The GIG's sampler, the modified half-normal, the asymmetric Laplace | `GeneralizedInverseGaussian`, `ModifiedHalfNormal`, `AsymmetricLaplace` (`Maths.Distributions`); `GigParameters.ToDistribution` |
| Exponential mixing | `ExponentialMixing` (`src/maths/regression/MixingLaw.cs`) |
| The noise on the config | `BassNoise.AsymmetricLaplace(q, reference)`, `BassNoiseReference.Quantile` |
| σ² on the record, the target reading it (decision 1, w included) | `ObservationLatents.Sigma2`; `BassTarget.GivenAt`, with ln p(σ²) in the state's prior |
| σ²'s block (decision 2) | `BassScaleBlock`: the modified half-normal on σ⁻¹, K and Q from `GPriorConditional.PriorEnergy` through the prior's pencil |
| v's block with the drift, the offset following | `BassLatentScaleBlock`, `VarianceMeanMixture.At` |
| The retention lever | `BassDraw.Latents`, `WorkingOn`, `LatentsOn`; `IBassObservationModel.Rebuilds`, `Rebuild` |
| The readers | the asymmetric Laplace's density, mean shift, quantiles, noise sd, and the draws' mean offset |

**The gates.**
- *The skew's sign* (E1's fixture): on symmetric noise, N = 400, the held-out share of 2,000 responses below the surface is 0.903 at q = 0.9 and 0.106 at q = 0.1, which a flipped sign would swap.
- *The q-quantile recovered,* on f + (0.3 + 0.4x₂)(G − 2), G ~ Gamma(2, 1), skewed and heteroscedastic, N = 400, four data sets: held-out share below the surface 0.518 at q = ½ and 0.884 at q = 0.9, within the ±0.04 asserted. Readings beside it: the surface's RMSE against the true quantile, 0.160 and 0.325, and f's 90% intervals covering the true quantile at 76.5% and 63.8% of a 7 × 7 grid.
- *`qbass`'s readouts* (`GbassQuantileReadoutTests`, q = 0.9, N = 200, w's prior IG(0.1, 0.05) on both sides): σ² 0.00768 [0.00587, 0.00990] under the g-prior and 0.00820 [0.00640, 0.0103] under the ridge against `qbass`'s w 0.00800 [0.00624, 0.00991]; the surface's means within 0.49 of the wider 90% half-band at worst (0.65 under the ridge); M 4.94 and 3.96 against 4.87. The differences stated beside it: GBASS's continuous knots and whole-term mutation, pseudo-counts 1 against 5, its raw ridge with the intercept inside, λ sampled, and the engine's priors reading the sample 0.9-quantile. R's run is 3,000 iterations, since 6,000 overran the oracle's 120 s at below-normal priority on the loaded machine.
- *The engine's GIG and modified half-normal against `gbass_primitives_oracle.R`'s law:* 11 and 9 cases, 50,000 draws each, worst |z| 2.56 and 2.55 at R's deciles and means.
- *Gaussian untouched:* the eight pinned hashes and the seams' suites inside the whole regression assembly, 1,984 of 1,984 (slice 1's 1,959 and this slice's 25), and the GBASS oracle classes 24 of 24 with R present. A Gaussian config's JSON still has no `Noise` key. The gated tests return at once in that run; theirs above ran with `BASS_QUALIFICATION` = 1.
- *The law, by Geweke's test* under the skew: four cases (q = 0.3 on two seeds, q = 0.8 with τ held, q = ½ where the drift vanishes), 200,000 marginal-conditional draws against 300,000 sweeps each, every |z| ≤ 2.95. Its power: flipping the skew's sign in σ²'s block sends σ² to infinity within seconds and fails every skewed case.
- *Primitives, chosen cases and seeded sweeps:* the GIG on 134 cases (worst decile z 3.01) and four intervals; the modified half-normal on 134 across its three branches (3.32, 2.35, 2.57); the asymmetric Laplace against the exponential mixture to 12 digits; exponential mixing's marginal over 305 cases, 8.7·10⁻¹⁵.
- *Refusals and acceptance* (decision 4): the latent and scale blocks refused nothing: 0 of 300,000 each in every Geweke case, 120 of 120 per block and chain in the resume fit.
- *The cost row* (thread cycles, loaded, indicative): a quantile sweep is 1,810 µs at the base under the g-form against a Gaussian sweep's 45 (40×; t's was 24×), 360 against 5.1 at n 200, 381 ms against 1.7 at M 200. v's N GIG draws are 53% of it at the base, about 1 µs a draw; Follow after them 29%, and 70% at M 200. σ²'s move is cheap: the offset alone moves. The table is on the page.
- *Retention, before and after:* 313.4 KB a retained draw with its working response against 70.4 KB with its latents alone at N = 10,000 (4.45×), 32.7 against 9.0 KB at N = 1,000. Views rebuilt from v equal the stored ones bit for bit (values, offset, weights, ȳ_w and the weighting, 40 draws each under the asymmetric Laplace and t), and a fit resumes bit for bit, latents and σ² included.
- *`docs/Maths`:* the new page; `ProbabilityAndInformation.md`'s section on the GIG, the half-normal and the asymmetric Laplace; the notation ledger takes q and level, and the half-normal's α, β, γ.

**Defects.** *Corrected here* (slice 0's finding): GBASS's w step at β ≠ 0 drops w's prior a and its floor; the engine's σ² priors are inverse gammas whose shape and scale both enter σ²'s modified half-normal, and it keeps no floor. Recorded on the page; the corrections register in `BayesianRegression.md` takes no new row, at the user's word not to grow that file. *Found and fixed in this slice's own first forms,* by the seeded sweeps: the GIG's log ratio cancelling where ω ≪ |p|; a hat reach that put the envelope's far tangent where h is −4·10¹⁶, whose rounding swamped the hat; the half-normal's left reach landing on x = 0.

**Decisions taken beyond the brief.**
- *The GIG is drawn by the tangent hat in ln x,* where every GIG is log-concave, in place of Hörmann and Leydold's three regions: one exact method across all parameters and far ranges, and a truncated draw for free. It costs about 1 µs a draw.
- *The modified half-normal by exact routes of its own* (the hat in x for α ≥ 1, in ln x for γ ≤ 0, the series mixture otherwise) in place of Sun, Kong and Pal's rejection schemes, each branch checked against quadrature and R.
- *The quantile read by default:* the sample q-quantile and the squared mean check loss; `Moments`, `Robust` and `None` selectable.
- *The draw keeps its latents for any model that says it rebuilds* (`Rebuilds`), so t's draws keep v alone too, halving theirs.
- *The quantile recovery's assertion is the held-out share,* with the RMSE and the coverage as readings: no prior standard exists for them.

**Left unbuilt, and why.** A free skew and γ (slice 3), generalized beta prime (4), the engine's reach (5), BARS (6) and interweaving (7) are later slices. A skewed t too: its conditional is a GIG the new sampler draws, but its marginal needs Bessel K of real order, so both wait for slice 3's general GIG mixing; t mixing's refusal now points there. The shared envelope behind every truncated draw is unchanged: the GIG's own reaches avoid the overshoot.

**Follow-ups, with the decision each would change.**
- *Calibrated intervals for the quantile.* Decision: whether quantile BASS reports its intervals as credible for the q-quantile. Under other noise they undercover (63.8% for nominal 90% at q = 0.9); Yang, Wang and He's (2016) adjustment, or the learning rate of [the calibrated-intervals entry](extensions-brief.md#calibrated-intervals-for-loss-based-posteriors), would widen them toward nominal. Left, the page states the undercoverage.
- *A faster GIG sampler as a measured option.* Decision: the GIG's default sampler under quantile BASS. Hörmann and Leydold's ratio-of-uniforms would cut v's draws, 53% of a sweep at n 1,000; it changes the default only if effective sample size per second rises.
- *A guard on the envelope's design point.* Decision: whether the shared sampler behind every truncated draw bounds its Newton search where h falls astronomically. It changes draws only where today's hat would accept nothing, which no existing family's reach reaches.

### Slice 1: Student t, robust BASS (chip, 2026-11-05)

**Landed** in `c6e0baf` (the code and its tests) and the commit that carries this report (`docs/Maths`, the notation ledger). [The docs' section](../../../docs/Maths/VarianceMeanMixtures.md#student-t-robust-bass-e1s-slice-1) holds the forms, the alternatives and the checks; this is the account.

| Piece | Where |
|---|---|
| The mixing law, shared with BARS: prior draw, conditional at power b, marginal, E[v], Var[v] | `IMixingLaw`, `StudentTMixing` (`src/maths/regression/MixingLaw.cs`); `SpecialFunctions.Log1p` made public for it |
| The noise on the config, left out of the JSON when null | `BassConfig.Noise`, `BassNoise` (`StudentT(ν, c, reference)`) |
| The carrier: the latents on the working response (decision 1) | `ObservationLatents`; `WorkingResponse.Latents`, `Moved(offset, weights, latents)` |
| The observation model and its blocks | `VarianceMeanMixture`, `BassCoefficientsBlock`, `BassLatentScaleBlock`; `BassStructureMoves.DrawSweepFit` |
| The noise law for the readers | `IBassObservationModel.Law` (null for the Gaussian), `IBassNoiseLaw`; pointwise, prediction (`MeanShifts`, `ObservationMean`, the mixture's quantile over `IContinuousDistribution`s), the training fit's `NoiseSd` |
| The run and the fingerprint | `BassModel.Create` reads `Noise`; `Bass.Fit`, `Resume` and `Fingerprint` take it through the config |
| The blocks' counters per chain | `BassChainAccount.Blocks` |

**The gates.**
- *v's conditional equals `StudentTWeights`' number for number:* from one stream, the law's draw is 1/λ of `StudentTWeights`' weight exactly, 500 residuals at each of ν = 0.5, 3, 4, 30, and the marginals agree to 11 places there. Against R, GBASS's own v conditional under t mixing at a frozen state is that Gamma at every observation (slice 0's `GbassConditionalParityTests`).
- *Large ν reaches Gaussian BASS's readouts:* ν = 10⁷ against Gaussian BASS, 16 chains of 1,000 draws each, N = 200: σ² 0.23295 against 0.23268 (z −0.47), M 4.275 against 4.239 (z −0.92), the fit at 200 cases rms z 1.06, max 2.54.
- *`tbass`'s readouts* (`GbassStudentTReadoutTests`, ν = 3, c = ⅓, N = 200, R 6,000 iterations, the engine 4 chains × 1,500 draws): w 0.230 [0.180, 0.289] under the g-prior, 0.235 [0.181, 0.297] under the ridge, against 0.224 [0.174, 0.283]; the surface's means within 0.23 of the wider 90% half-band at worst over 15 grid points, rms 0.13 (under the ridge, worst 0.34, rms 0.19); M 4.17 and 4.63 against 5.29. The differences stated beside it: GBASS's continuous knots and whole-term mutation, pseudo-counts 1 against 5, its raw ridge with the intercept inside and τ ~ IG(½, N/2), its floor on w, λ sampled under Gamma(1, 1) where D02 integrates it.
- *Gaussian untouched:* the eight pinned hashes and the seams' suites (bit identity, observation model, sweep slots, weighted core, moving response, prior reference, contract, resume) 316 of 316 after the change; the whole regression assembly 1,959 of 1,959. A Gaussian config's JSON has no `Noise` key, so every fingerprint holds.
- *The law, by Geweke's test:* N = 12, M ≤ 2, proper priors that read no response, 200,000 marginal-conditional draws against 300,000 successive-conditional sweeps, eight functions on four cases (c = 0.7 on two seeds, c = 1, τ held): every |z| ≤ 1.78. Its power, by a deliberate break: dropping c from v's unit variance reads z 75 on mean ln v. The cases take τ ~ Gamma(5, rate 5N) or τ held; under Zellner–Siow's shape ½ the successive chain sticks in ln τ's long left tail (its mean 0.30 above the prior's on one seed, 0.58 below on another, the other seven functions in place), the test's reach rather than the sweep's.
- *Recovery from truth outside the model:* contaminated normal noise, 0.85·N(0, 0.25²) + 0.15·N(0, 2²), on sin(2πx₁) + x₂, N = 300, six data sets. t(3) BASS's 90% intervals cover 89.5% of 294 grid points of the true surface, RMSE 0.059; Gaussian BASS's 86.4%, RMSE 0.142. Held-out responses inside the 90% predictive intervals: 88.3% under t, 92.7% under the Gaussian.
- *Primitives, chosen cases and a seeded sweep:* the law's marginal against quadrature of the mixture, 11 chosen and 400 swept cases, ν from 0.05 to 10¹¹, s² over 24 decades: worst relative gap 1.4·10⁻¹³. The conditional's draws against the Gamma's moments at b = 1, 0.4 and 0. The log-gamma ratio across its switch. The sweep found the quadrature reference's own cancellation at ν = 10⁶ (2.5·10⁻¹⁰), fixed in the reference, and BARS's direct forms (below).
- *Refusals and acceptance* (decision 4): the latent block refused no draw anywhere: 0 of 300,000 in each Geweke case, 0 of 120 per chain in the resume fit, 0 in 50 sweeps on the 80-case fixture. The structure's own moves are counted as before. So F under a moving W stays on the gate that reads W, with nothing yet to weigh against it.
- *The cost row* (thread cycles, loaded, indicative): a t sweep is 961 µs at the base under the g-form against a Gaussian sweep's 40 (24×), 83 against 2.9 at n 200 (29×), 388 ms against 1.8 at M 200, n 10,000 (216×). Follow, the weighted Gram's refactoring at O(NM²), is 55% of it at the base and 74% at M 200; then the structure step, the pencil's O(p³) per proposal; then v's N draws. The table is in the docs.
- *Retention:* 15.6 KB a retained draw at N = 1,000 and 1.5 MB at N = 100,000, its working response's weights and v, against about 0.2 KB of its own at M ≈ 10. It dominates a draw. The weights are a function of v, so keeping v alone halves it under t; under skew, where the offset and the shifted response join, it would quarter it.
- *`docs/Maths`:* the section above, the seams' table, what GBMARS still writes, the corrections table (GBASS's moments where infinite), `ProbabilityAndInformation.md`'s mixing law, the source map. The notation ledger takes w, β, τ, c and v between GBASS and the core.

**Defects.** *Found:* BARS's `StudentTWeights.MarginalLogDensity` takes ln Γ((ν + 1)/2) − ln Γ(ν/2) and ln(1 + x) directly, which lose digits as ν grows: on the law's sweep it reads up to 2.0·10⁻⁵ off, relative, at ν = 3·10¹⁰, where the law holds 1.4·10⁻¹³. No ν BARS takes by default comes near it. Slice 6 routes it through the law. *Corrected:* GBASS's s2 where the moments are infinite (slice 0's finding): the training fit states no noise sd at ν ≤ 2.

**Decisions taken beyond the brief.**
- *σ² and a ride the chain's own moves within a sweep* (`BassStructureMoves.SweepFit`), not the latent record. Decision 1 put them on the record; but a working response is shared by every chain at the start and by retained draws, so a value that changes within a sweep could ride it only as a copy per sweep, or as shared mutable state between chains. The latents themselves, v and β, are on the record as decided.
- *What the priors read:* `Moments` by default, which restates D03's read and what Gaussian BASS reads (the moving-response note settles this row at β = 0 under the unweighted Gram); `Robust` and `None` selectable.
- *c* defaults to 1, `StudentTWeights`' t; (ν − 2)/ν, `tbass`'s, is a value of the same field.
- *Refused with their slice named:* a tempered rung, SMC and the flat histogram (slice 5), a direct functional response (slice 5), a drift (slices 2 and 3), and the collapsed predictive under a noise law, until it has its own check.
- *GBASS's floor on w is not ported:* σ² is integrated in the structure step and redrawn with a, with no floor, as D04 has it.

**Left unbuilt, and why.** The retention lever, v alone on a draw with the rest rebuilt on demand: under t it halves the 15.6 KB, and slice 2, where the offset joins, is where it quarters it, so it belongs with the skew. The cost row's direct functional point and the weighted prior Gram's reading: the first is refused in this slice, and the second is no valid chain under a latent scale without the options' step 3.

**Follow-ups, with the decision each would change.**
- *Retained draws keep v alone.* Decision: what a retained draw holds under a latent noise. Built, retention halves under t and quarters under skew, and the views rebuild the working response from v, O(N) per draw read. Left, long runs at N ≥ 10⁵ hold 1.5 MB a draw. Slice 2's to take.
- *Structure sweeps per latent move,* the brief's lever after the row. Decision: the default cadence, one each (D04). Follow is most of a t sweep, so k structure sweeps per latent move cut its share by k; it changes the default only if effective sample size per second rises, read across fixtures, with the default staying one otherwise.

**Open, with no decision waiting.** Whether the Geweke test reaches Zellner–Siow's τ at a length a session can run.

### Slice 0: the oracle modes (chip, 2026-11-05)

**Landed** in `ee2cce6`: four modes under `r/oracles` that run from `r/` against the pinned GBASS 2.0.1, with their parity classes in `tests/oracle` and rows in `r/README.md`'s oracle map. 21 tests, each run against R and each printing its readings, 21 passed (`Maths.Oracle.Tests`, filter `Gbass`, about 85 s together on the loaded machine). A class returns early only when R is absent; a missing package fails it.

| Mode | What it evaluates | Parity class, and what the C# side computes |
| --- | --- | --- |
| `gbass_primitives_oracle.R` | `GIGrvg::rgig`, `rgig2`, `mu_gig`, `var_gig`; `rmhn`; base `besselK` | `GbassPrimitiveParityTests`: each law's integral in logs from `GigParameters.LogKernel` and the kernels, so the reference is computed on both sides |
| `gbass_conditionals_oracle.R` | `gbass`'s and `nwbass`'s conditionals for v, β, γ, τ, w, a and a birth's likelihood, as written, against the joint GBASS states | `GbassConditionalParityTests`: `StudentTWeights`' step, `ConjugateUpdates` for τ and w, β's and γ's closed forms |
| `gbass_defects_oracle.R` | E1's defects, each beside a reference that shares nothing with the code | `GbassDefectParityTests` |
| `gbass_readouts_oracle.R` | One wrapper's surface over a grid, w, β, γ, τ, M, the share below the fit, move counts | `GbassReadoutParityTests`; the engine meets it from slice 1 on |

**The gates.** The modes run from `r/`. Each defect is confirmed or struck below.

**Primitives** (`ours` here is a quadrature, since the engine's GIG sampler, modified half-normal and log Bessel K are slices 2 and 3; they meet the same mode then):
- *The GIG is one parameterization.* `GigParameters`' kernel integrated here gives R's log mass to 9 places and its moments to 9 relative digits on 11 cases, the limits b = 0 (the qbass prior and a residual of exactly zero) and a = 0 (t mixing) included, with the Bessel K moments where a, b > 0. GIGrvg's draws and `rgig2`'s carry the law's probabilities at their deciles, worst |z| 2.23 over 50,000 draws a case.
- *`rmhn`* matches quadrature on all four branches: worst decile |z| 2.09, the 0 < α ≤ 1 series branch at α = 0.3, 0.6 and 1 included, which E1's scratch had left unchecked. Its moments there also equal the series Σₖ (γ/√β)ᵏΓ((α + k)/2)/k! read in R to 15 digits. No defect.
- *`besselK(x, ν, expon.scaled = TRUE)`* equals K's integral representation to 10⁻¹¹ relative on nine points from (10⁻⁶, 3) to (800, 40) and (20, 250), with no warning: a reference for slice 3's log Bessel K.

**Conditionals at a frozen state** (N = 40, M = 2, two states: skewed with general GIG mixing; t(5) at β = 0). The joint less each conditional's log kernel is constant to 1.1·10⁻¹³ nats over eight points for v, β, τ, w's GIG at β = 0, and γ on γ > 0. a's mean and covariance as written equal dense algebra to 10⁻¹⁰, and a birth's likelihood part equals the difference of log N(z | 0, wcV + wτBBᵀ) to 9 places. On the C# side: under t mixing λ = 1/v is `StudentTWeights`' Gamma((ν + 1)/2, (ν + r²/(wc))/2) for every observation, shape exactly and rate to 14 digits; τ is `ConjugateUpdates.NormalVariance` over M + 1 coefficients, and w at β = 0 its GIG form (N-X3) over N + M + 1 terms, to 13 digits.

**Defects.** Each is a correction the engine takes, not a parity variant.

| Defect | Run | Verdict |
| --- | --- | --- |
| `nw_est_mom`'s kurtosis | The mixture's fourth moment by quadrature over the Wald, on both sides: 6.4641 at (α, β, δ) = (2, 1, 1), where the code gives 5.1651; 10.4536 for 8.0684 at (3, −2, ½); 3.5920 for 3.3085 at (1.5, 0.3, 4). The corrected form equals the mixture to 9 digits | Confirmed |
| `nw_gamma_prior`'s stray −1 | P(ξ > 0.9) at its start (100, 100): 0.001135 by definition, 0.003537 as written; P(ξ < 0.1) right | Confirmed |
| The w step at β ≠ 0 drops w's prior a and the floor | The modified half-normal's parameters are the same at a = 0 and a = 5; the joint less its kernel spreads 0.038 nats at a = 0.4, 10⁻¹³ at a = 0. With the floor at the conditional's median, 4,000 runs of the w block as written: at β = 0 none below it (49.7% keep the old w), at β = 0.4 59.2% below it | Confirmed |
| Birth admits npart rows, mutation needs more | The two conditions as written on a column of exactly npart positive rows: birth admits, mutation refuses | Confirmed |
| The GBP v step scaled from mean(ln v) | `hbass` at N = 60, 3,000 iterations: adapted step 0.543 against 2.38 × sd(mean ln v) = 0.641 and 2.38 × the median per-observation sd of ln vᵢ = 3.96; acceptance 0.90. A tuning defect: adaptation stops at burn-in, so the chain stays valid | Confirmed |
| γ's draw untruncated | At N = 5, Σv = 20, (m, s) = (0.1, 3), the code's normal puts 0.132 below 0, where it is no conditional: the joint less its kernel spreads 58 nats on γ < 0 and 10⁻¹³ on γ > 0. Under a prior on γ > 0 that mass is the draw's total-variation error | Confirmed |
| The text's r = ŷ − y | `r <- y - yhat` as written; the β conditional is the joint's for r = y − ŷ | Struck for the code |
| Table 1's t row, b = ν/2 | `tbass`'s v prior as written is GIG(−ν/2, 0, ν), and √v·ζ's density is dt to 9 digits; b = ν/2 gives √2·t(√2x) | Struck for the code |

**Found here.**
- *`mu_gig` and `var_gig` where a moment is infinite.* At a = 0 they return −b/(2(p + 1)) and −b²/(4(p + 1)²(p + 2)) for any p < 0: −3 and −7.2 at tbass's ν = 1.5, −18 for t(3)'s variance. So `tbass` reports s2 = −w at ν = 1, −3w at ν = 1.5, NaN at ν = 2 (and bias NaN, 0·∞). The engine's `NoiseSd` is null there, as this brief already said.
- *`nw_gamma_prior`'s default targets are out of reach.* (q₁, q₂, p₁, p₂) = (0.1, 0.9, 0.5, 0.05) asks for half of γ's mass above 99 and 5% within 0.235 of 0, which no normal holds. As written the fit lands at (−22.3, 90.4) with probabilities 0.288 and 0.0020; corrected at (−1.42, 82.7) with 0.231 and 0.0023. A calibration on the triangle should state targets a normal can reach, or a prior on log γ; slice 3's question.

**Readouts** (N = 200, f = sin(2πx₁) + x₂, noise 0.3·t(3), 2,000 iterations, 1,000 burn-in): each wrapper runs in 9–11 s on the loaded machine. `qbass` at q = 0.9 has 90.5% of the responses below its surface in sample, with β held at 1/q − 1/(1 − q). `tbass` at ν = 3 reads w = 0.232 [0.181, 0.294] against the noise variance 0.27. Mean M 4.5–5.1, acceptance 0.10–0.18 per move.

**Refusals and acceptance; the cost row.** Nothing of the engine runs in this slice, so neither applies; R's own move acceptances are above.

**Left for later.** The engine's GIG sampler, modified half-normal and log Bessel K meet `gbass_primitives_oracle.R` at slices 2 and 3; `qbass`, `nwbass` and `hbass` meet the engine there and at slice 4. Readings, not proposals: `lower_bound` is a valid independence step that kept the old w half the time with the floor at the median.

### CRAN's GBASS pinned (main session, 2026-11-04)

Slice 0's precondition. `r/DESCRIPTION` takes `GBASS`; `renv::install("GBASS@2.0.1")` and `renv::snapshot` add exactly three packages to `r/renv.lock`: GBASS 2.0.1, BASS 1.3.1, GIGrvg 0.8. `renv::status()` is consistent.

First run of the package, from `r/` on a loaded machine, N = 300, 2,000 iterations: `tbass` on t(3) noise 21.7 s (mean M 6.3, mean w 0.83); `qbass` at q = 0.9 on Normal noise 23.1 s (mean M 5.8), with 91.0% of the responses below the fitted surface, in sample; `rmhn` runs. No comparison with the engine yet, and no defect has been run.

**The harness, on the way.** R BASS's checkout had moved to `D:/aghado01/packages/bayesian-related/BASS`, and `bass_source` and five parity classes (`BassBasisParityTests`, `BassMomentParityTests`, `BassSensitivityDefectParityTests`, `BassTargetParityTests`, `WalleniusParityTests`) still named the old default. Those five `return` when the folder is missing, so they had been passing without running. Their default now names the new folder and they run: 33 of 33 in 5 m 24 s, loaded. The other oracle classes were not rerun. The silent `return` stands.
