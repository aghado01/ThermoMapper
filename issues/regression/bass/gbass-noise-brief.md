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

**Recommended; a chip starts on these unless the user says otherwise.**
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

- **Gaussian untouched.** The eight pinned hashes unmoved and every existing fingerprint the same.
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

- **G1's exact machinery.** With v moving, the structure is no longer a chain on an enumerable space. Whether heavy tails raise or lower the barriers between basins is unread, and worth a chip once slice 1 runs.
- **Stefanski's logistic mixing (E20).** In the clone only; its route to logit is E20's.
- **KHAOS.** Next ([the build order](lineage-review.md#build-order)), in a brief of its own, with CRAN's `khaos` 2.1.1 pinned the same way first. It shares this build's lineage and its method, not its pieces: it works on the basis, the proposals and the coefficient prior (E19, E4, E5, E21), which this build leaves alone, and it moves no latent. Once both exist their product costs nothing more: a polynomial-chaos basis under quantile or robust noise.

## Working rules

- On `main`, committing by pathspec with nothing left staged and nothing pushed. Commit messages describe the work, with no attribution lines ([CLAUDE.md](../../../CLAUDE.md)).
- Builds set `ArtifactsRoot` inline into `artifacts\chip-gbass-noise\`, deleted at the end.
- Heavy runs at below-normal priority with capped threads; the machine is shared, and basin-hold's chains were running when this was written.
- New beside old; a removal is its own question to the user, with what is lost.
- Findings framed by problem class, not by any one application.

## Report

Under "Reports" below, newest first, per slice: what landed by commit; the gates' numbers; the parity comparisons with their stated differences; refusals and acceptance; the cost row; defects confirmed, struck or newly found; decisions taken beyond this brief; what was left unbuilt and why. `docs/Maths/BayesianRegression.md`'s table of what each extension still writes shrinks as slices land.

## Reports

Newest first.

### CRAN's GBASS pinned (main session, 2026-11-04)

Slice 0's precondition. `r/DESCRIPTION` takes `GBASS`; `renv::install("GBASS@2.0.1")` and `renv::snapshot` add exactly three packages to `r/renv.lock`: GBASS 2.0.1, BASS 1.3.1, GIGrvg 0.8. `renv::status()` is consistent.

First run of the package, from `r/` on a loaded machine, N = 300, 2,000 iterations: `tbass` on t(3) noise 21.7 s (mean M 6.3, mean w 0.83); `qbass` at q = 0.9 on Normal noise 23.1 s (mean M 5.8), with 91.0% of the responses below the fitted surface, in sample; `rmhn` runs. No comparison with the engine yet, and no defect has been run.

**The harness, on the way.** R BASS's checkout had moved to `D:/aghado01/packages/bayesian-related/BASS`, and `bass_source` and five parity classes (`BassBasisParityTests`, `BassMomentParityTests`, `BassSensitivityDefectParityTests`, `BassTargetParityTests`, `WalleniusParityTests`) still named the old default. Those five `return` when the folder is missing, so they had been passing without running. Their default now names the new folder and they run: 33 of 33 in 5 m 24 s, loaded. The other oracle classes were not rerun. The silent `return` stands.
