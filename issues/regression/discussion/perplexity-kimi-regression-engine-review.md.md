
# so i've been working on my bayesian regression engine, heres some snapshot. i still need to add a few more files next turn

The back half fills in the BASS execution engines and the downstream analysis surface:

**Execution engines.** Three distinct sampling paths now visible: the MCMC `BassEnsemble` with optional tempering ladders (`BassLadder`) and evidence reporting gated on a reference rung at b = 1; the SMC path with replicate-level outer parallelism, adaptive inverse-temperature stepping via conditional ESS, per-replicate weight pooling, and replicate-spread diagnostics standing in for R-hat; and the multicanonical/flat-histogram path (`BassMulticanonicalChain`) with a learning-phase tally reset and reweighting against learned histogram weights.[^2_2][^2_3]

**Noise as scale mixtures.** `BassNoise` is a genuinely deep family — Student-t, asymmetric Laplace, Normal-Wald, generalized hyperbolic, and generalized beta prime (with horseshoe and Strawderman–Berger as GBP specializations) — all realized through a single `VarianceMeanMixture` observation model with per-observation latent scales $v_i$, skewness drawn or held, and precision weights $1/(c\,v_i)$. The sweep slots (Structure → Tau → Coefficients → Latent → Scale) enforce that latent-scale and scale blocks read the current sweep's coefficient fit, with explicit `InvalidOperationException` guards on block ordering.[^2_2]

**Structure prior.** Per-part priors (continuous / categorical / functional) with uniform order, uniform variable subsets, categorical subsets either size-uniform or measure-uniform over $2^n - 2$ proper subsets, and a negative-binomial count law whose log-collection weights are lazily cached in a 1024-entry NaN-initialized table — the comment noting τ's slice step reads it six times per sweep.[^2_2]

## Two things I noticed

- There are roadmap markers embedded in production guards: `VarianceMeanMixture.Working` and `BassCoefficientsBlock.Update` both throw `NotSupportedException` citing "E1's slice 5" — per-observation latents under functional designs and latent-noise blocks on tempered rungs are deferred work. If tempering + robust noise is on your critical path, those are the seams.[^2_2]
- The SMC convergence story is honest but structurally weaker than MCMC's: convergence diagnostics are all `Unavailable` and the verdict is `Undetermined`, with `BassSmcSpread` (between-replicate SD / pooled SD) as the surrogate. Worth being deliberate about when you surface that readout as "converged."


## What stands out

**The flat-histogram module is the most novel machinery here.** It's a proper two-phase SAMC design: a learning phase with diminishing-gain stochastic approximation (`SamcLearning`, gain scale t₀ defaulting to 10 × node count), WHAM warm-start from a tempering pilot, runtime top-window extension (`ExtendTopUntil`, with the occupancy histogram extended alongside), and then a production phase that freezes weights and reuses the learned chains — checkpointable at both boundaries. R-hat-style guarantees are replaced by explicit flatness diagnostics: per-chain round-trip counts across the window ends, `ChainsWithoutRoundTrip`, and max flatness deviation against the trapezoidal `FlatOccupancy` target. That's a much more honest convergence story than the SMC path's `Undetermined` verdicts.[^3_3]

**The reweighting self-consistency is well enforced.** `MulticanonicalEvaluation.LogLikelihood(b)` throws for any $b \neq 1$ with the exact right reason — $\pi_w$ is not on a tempering path — and `MulticanonicalChain` rejects any chain not at `Temperature.Posterior`. The `Reweighting.LogEvidence` curve additionally gates on `HasCompleteValues` so a tilted log-normalizer can't be silently misread as evidence, mirroring the same `IsProper`/`HasCompleteValues` guards in `BassEnsemble.ReportsEvidence`. This is a consistent invariant propagated across three modules, which is the kind of thing that's easy to get wrong in exactly one place.[^3_3]

**Proposal abstraction is minimal and clean** — `Proposal<TState>(Candidate, LogProposalRatio, LogJacobian)` as a single record struct, with the RJ accept logic presumably centralized in `ReversibleJumpChain` rather than smeared across moves.[^3_2]

## What's notable in this layer

**The factorization story is richer than a default.** `ISpdFactor` has four concrete SPD paths: dense, banded (`BandCholesky`), wide-band (`WideBandCholesky`, row-packed with $p+1$ entries per row and a half-bandwidth gate at order ≥ 64), and `BorderedCholesky` — the bordered one is presumably what services g-prior rank-one structure updates in `BassDesign`. Every factor carries `FactorStatus` + a pivot tally + `EquilibratedCondition` estimation rather than throwing on the first non-PD pivot, which matches the `BassRefusal` tally design: numerical degeneracy is a first-class, counted outcome, not an exception. And `SpdFactors.ForBand` auto-selects representation by the 2(p+1) ≤ n rule.[^4_2]

**`ThinSvd` is a serious implementation**, not a wrapper: one-sided Jacobi with Householder pre-reduction, scaled reflectors to keep vᵀv in range, Drmač-Veselić-style rotation angle with the min/max ratio trick, rank-tolerance-based null-space completion via re-orthogonalized pivot search, and overflow-safe norms with a 1e±280 rescale branch. AVX2-vectorized throughout.[^4_2]

**The RNG layer is exactly what reproducibility needs and no more.** `SeedTree.Derive` is SplitMix64-based hierarchical seeding (the same SplitMix64 that warms xoshiro's state), and `Xoshiro256PlusPlus` has `SaveState`/`Restore` with an explicit all-zero-state guard — that's what makes the per-chain checkpointing in `BassChainRun.SaveState` actually sound, since the composition RNG, move RNGs, and readout RNGs all go through it. Rejection-free `NextInt` via the multiply-shift trick, too.[^4_2]

**`PairwiseSum`** — recursive pairwise blocked at 128 with Vector256 lanes — is the quiet workhorse: every Gram, every log-density accumulation, `WorkingResponse.ResponseSquares`, the prior-reference moments all route through it, so the statistical layer gets O(log n) error growth instead of O(n) for free.[^4_2]

One curiosity: `linalg` also carries PCA/ICA/LOBPCG/FFT — none of which anything in the regression tree references. That suggests either an older feature extraction layer (pre-`BassInputTransform`?) or machinery reserved for the ThermoMapper application layer itself.[^4_1]


## Review 

Here's the review. I've now read deep enough into the load-bearing files — `ParallelTempering`, `SamcLearning`/`Wham`, `BassStructureMoves`/`BassTarget`, the conformal scorer, GP/changepoint/BARS, and the primitive blocks — to go beyond structure into judgment calls. Overall: this is research-grade work with an unusually consistent discipline (refusal tallies instead of exceptions, fingerprints everywhere, checkpointability down to per-RNG-word state). The notes below are where that discipline has uneven coverage or where the substrate is leaving capability on the table.

## Gaps and weak seams

**SMC is the soft underbelly, and I don't think readout alone fixes it.** Beyond the missing convergence diagnostics I flagged earlier, `BassSmcKernel.Move` runs exactly one RJ step + static blocks per stage with no within-stage adaptation — mutation intensity is set by `MutationSteps` in the spec, not by anything the kernel learns. For the RJ dimension moves this is usually fine, but the failure mode is a replicate that drifts into a high-b dimension regime where the single-step birth/death palette can't traverse, and replicate spread won't catch it because a *stuck* replicate converges internally. Concrete addition: per-stage acceptance-rate-adaptive mutation counts (double mutation when stage acceptance drops below a floor, capped by budget), plus a post-hoc per-replicate ESS trace from the stage reports you already store — you have the data, just not the readout.[^5_1]

**The E1-slice-5 frontier is now the sharpest capability gap.** Three separate guards — functional-design latents, latent-noise blocks on tempered rungs, working-response moves under response-reading priors — all fence off the same region: non-Gaussian observation models inside tempering/evidence workflows. `ProposeFollow` at least gives a principled error message instead of letting it happen silently. But this means skewed-t / Normal-Wald noise is currently MCMC-only, unevidenced, which matters for model comparison across noise laws (arguably the main reason to have five mixing laws). The pragmatic path to evidence for mixture models is not lifting the tempering restriction — it's the SMC path, where `BassParticleTarget` has no b=1 block constraint in principle; the latent blocks just need stage-aware guards, not outright refusal.[^5_2]

**BassTauBlock's fixed slice width is the one un-adapted knob in an otherwise adaptive system.** Width 1.0 / 32 expansions against a posterior that shifts by order-of-magnitude in τ across the run's life (it also drives the collection-weight cache 6×/sweep) is a classic place for a slow mixing tail. Cheap fix with big mixing payoff: adapt the width online from the accepted interval span (e.g., shrink-to-half on >8 expansions, grow on zero), reset at burn end like `ResetTallies`. Same for `BassSigma2RateBlock`, which inherits the same constants verbatim.[^5_3][^5_4]

**`TemperingPilot.FromLadders` correctly refuses powered-path ladders, but the powered path is what an improper-priority fit needs.** `BassTarget.Defines` shows you're already reasoning carefully about where the integrated form is undefined (flat intercept, Jeffreys σ² at b→0) — the *linear* form is the escape hatch at b=0, and the ladder layer supports `TemperedForm.Linear`, but the pilot hard-requires `TemperingPath.Likelihood` with the (b_c − b)·ℓ₁ linearity check. A powered-path pilot would need a generalized slope functional instead of raw ℓ₁. If GBASS or KHAOS intends flat-histogram evidence under non-Gaussian noise (where likelihood-path b=0 may be improper), this check will be the wall you hit.[^5_5][^5_3]

## Performance opportunities

**ParallelTempering stepping is serial — the largest single-threaded cost left in the engine.** `Step()` walks replicas in a for-loop and then swaps; with DEO/schedules you have at most two synchronization points per scan, so replica steps are embarrassingly parallel. For BASS this is the whole ladder's cost — a 12-rung ladder runs 12 full BASS sweeps serially per scan on one core while sibling chains idle on others (BassEnsemble parallelizes *across* ladders, not within). A `Parallel.For` over replica steps with per-replica RNGs (already isolated) is a contained change since swaps are pure post-step bookkeeping.[^5_5]

**GP `FitCore` refactors from scratch per bandwidth evaluation.** The O(n³) Cholesky per (t, η) grid point is unavoidable for a dense kernel, but `KernelAffinity`'s O(n²) distance-pair sum is recomputed identically across t — no wait, it isn't, it depends on t. The real waste: `EmpiricalBayesBandwidthPrior.LogDensity` needs only `KernelAffinity(t)`, but any marginal-likelihood optimization over t pays full refactorization. Given `Refit(fit)` already reuses structure, a one-sided finite-difference-friendly memo of the last (t, η) factorization would cut Nelder–Mead over bandwidth roughly in half. At scale, the honest escape is a pivoted/incomplete Cholesky or FIC approximation — `DiagonalPlusLowRank` is sitting in `linalg` unused by anything in the tree, which makes me think it was built for exactly this.[^5_6][^5_7]

**The BASS τ-evaluation hot path duplicates marginal construction.** `BassEvaluation.LogValueAt(b)` calls `Marginal(fit, tau, b)` per tempering evaluation during ladder sweeps, and `Marginal` rebuilds the full `GPriorMarginal` (including CentredMarginalAt weighting/mean plumbing) for states identical up to logits. Since `BassStructureMoves` already caches column means and sweep columns keyed by term, a `(fit, tau, b)`-keyed last-value memo on the evaluation or in the collapse layer would eliminate redundant marginal work in swap-heavy scans, where each state gets evaluated at its own b *and* its neighbor's.[^5_3]

**`Wham.Estimate` is fine; `SamcLearning.SolveSpd` is a private hand-rolled Cholesky with a scalar ridge.** It duplicates `linalg`'s dense Cholesky functionality and, more importantly, ridge-floors silently at 1e-14·max — a Hessian that's drifting singular (thin pilot tails) gets silently regularized *inside* the solve rather than surfacing as a `FactorStatus` you could route into the refusal/diagnostics vocabulary you use everywhere else. Routing it through `DenseCholesky` + `FactorStatus` would align SAMC failure modes with the BARS design-admission gates — same philosophy, currently diverging.[^5_7][^5_8][^5_5]

## Capability runway

**For GBASS: the substrate is closer than the code layout suggests.** BARS already carries exactly the machinery BASS lacks for exponential-family responses: `PoissonModel` with IRLS + banded Laplace + g-prior terms + independence-MH conditional coefficient sampling against the Laplace factor, and `GibbsLossModel`'s robust IRLS. The seam is `IBassObservationModel`, which today speaks only "weights + offsets of a Gaussianized working response." A Laplace-family observation model would make `WorkingResponse` a second-order expansion and reuse `PoissonConditional.Coefficients` verbatim. The categorical-factor + hinge basis needs no changes — the whole GLM story lives below the basis.[^5_8][^5_2]

**A `TemperingPilot` from SMC stage traces would fuse your two evidence engines.** Wham only needs (b_i, draws of the slope functional, inefficiency) — `BassSmcReplicate` already captures stage-level temperatures, log-normalizer ratios, and stored trace stages. Folding SMC stages into a pseudo-pilot gives SAMC a warm start without ever running a tempering ladder, which is precisely the workflow for models where the reference rung is expensive or non-existent. Today SAMC↔WHAM↔ladder is one closed triangle; SMC sits outside it.[^5_1][^5_5]

**The conformal layer's stratum machinery is genuinely novel and one abstraction short.** The full-step Bayesian p-value computation with Pareto-k importance diagnostics, carried vs. uniform margin bounds, ESS-corrected error inflation, and the CB-4 latent-treatment consistency checks is well beyond anything in the standard literature — but it's singleton-grid only. Joint/multi-response conformal (for `FunctionalResponse` outputs — the Sobol/sensitivity surfaces) needs the same scorer running over output units with the latent-treatment checks generalized from `Observation` scope to grouped scopes. `LatentReport`/`LatentScope` already has the vocabulary (Group/Global scope appears in the changepoint posterior's pointwise density); the conformal checks just don't consume it yet.[^5_9][^5_6]

**The unused spectral/reduction stack (FFT, ICA, LOBPCG, Sparse PCA-ish `ThinSvd`, `SpectralNoiseField`) is a latent input-layer capability.** `SpectralNoiseField` (fractional Gaussian fields by FFT filtering) is clearly synthetic-surface generation for validation — if it isn't already wired into a test harness, `ExactChangepoint` + `SpectralNoiseField` + BASS recovery benchmarks would be the cheapest calibration evidence you can generate for the RJ structure moves.[^5_4][^5_7]

---

# GBASS aka Generalized Bayesian MARS and KHAOS are both Lancome lineage extensions of the original BASS concept, ABMRS is Adaptive Bayesian Multivariate Regression Splines, not technically Kass lineage, but sort of a generalization of it to multivariate, introduces some neat things and was only published in its final form recently. In all cases, the repo code associated with the papers tends to sell short of the theory of the papers and i'm aiming for the papers or better in my synthesis

That tracks — "the repos undersell the papers" is the real story of this lineage: BMARS/RJMCMC reference implementations are usually single-chain birth/death/change with a Gaussian collapse and no evidence story, and the papers' actual content (tempering, model comparison, robustness) lives in the text, not the code. Your engine already exceeds every published repo I know of on the sampler side — so the question becomes where the *paper-level* theory lives relative to what you have. Mapping the three targets onto the substrate:

## GBASS: the seam is already cut

My earlier read holds and sharpens — the generalized axis (exponential-family responses) is a below-the-basis problem, and your architecture agrees with that: `IBassObservationModel` owns the working-response construction, and BARS's `PoissonModel` is a complete worked example of the exact machinery the papers use: IRLS to the mode, banded LDLᵀ Laplace factor with g-prior Gram terms injected into the Hessian, and independence-MH coefficient draws against that Laplace proposal.[^6_1][^6_2]

The parts the papers have that BARS currently handles better than most repos, worth porting deliberately:

- **Laplace-quality discipline.** `PoissonModel.Score` returns the Laplace approximation *with the g-prior log normalizer separated* (`logNormalizer = 0.5·ν·log scale − 0.5·log|Gram|`), so the marginal is comparable across knot counts — that's the quantity GBASS's RJ acceptance actually needs, not the raw deviance. Keep that separation when you move it under `IBassObservationModel`. [^6_1]
- **Tempering honesty.** BARS's GLM marginals are `LinearMarginal` — evaluate the Laplace value once at b=1, then scale linearly in b, with `CanTemperToPrior` gating the reference rung. That's the defensible choice: re-running IRLS per rung is both expensive and approximation-inconsistent across b. Whatever GBASS does should inherit this, and your ladder/SAMC evidence chain will tolerate it because the linear-tempered slope functional is exactly what `TemperingPilot` validates.[^6_3][^6_1]
- **The one design decision to make now:** whether the GLM working response is still a `WorkingResponse` (offsets + diagonal weights) or a generalization of it. Laplace working responses are naturally (linear shift, second-order metric) — your `Moved`/`Rebuild`/`Rebuilds` triad is already carrying two semantic roles, and GLM adds a third. If you widen the seam once into "quadratic working response with optional non-diagonal metric" now, ABMRS inherits it for free later instead of paying a second refactor.

One caution on scope: the null-model reference rung (`KnotReference` in BARS samples count-and-knots under the tempered prior) is what makes GLM evidence *absolute* rather than comparative. Standard repos skip it and report unnormalized Laplace values that can't be compared across models; you already own the solution pattern via `IReferenceSampler`.[^6_1]

## ABMRS: you already built its skeleton

The striking thing from the tree is that the functional-BASS design is structurally *isomorphic* to the multivariate-response problem: `BassTerm` has a functional part with output-location knots, `BassBasis` produces tensor columns $a_i \otimes h_r$ (`CaseColumn` × `OutputColumn`), `FunctionalResponse` carries output locations plus a metric (`IInnerProduct`), and `BassPointInputs` already evaluates monitored functionals through that separable design. An ABMRS synthesis on this substrate doesn't need a new basis layer at all.[^6_2][^6_4]

The real work is the noise model, and this is where "papers or better" has genuine room:

- **The collapse survives the multivariate generalization.** With separable design $X = A \otimes H$ and matrix-normal slab $\mathrm{vec}(\beta)\mid\Sigma,\tau \sim N(0,\, \tau G^{-1} \otimes \Sigma)$ with $\Sigma \sim \mathrm{IW}$, the marginal is matrix-t and stays closed-form — the g-prior machinery (`CollapsedGPrior`, weighting, centered intercept) generalizes to Kronecker form rather than being replaced. Nobody's repo implements this cleanly; it's the difference between a multivariate model and "univariate BASS run q times."
- **Your existing `WeightedNormalModel`/`WhittleModel` in BARS give a second, cheaper axis** — correlated noise across output locations in the spectral domain, for the case where outputs are functional (temporal/spectral) rather than unstructured. Worth deciding ABMRS's stance on separable $\Sigma \otimes W$ vs. full matrix-normal vs. Whittle-style structure before the refactor.[^6_1]
- **Known missing piece, already on your own roadmap:** E1 slice 5's "per-observation latent under a direct functional design" guard is exactly the place where ABMRS's robust-noise variants will land. Fixing the seam once (previous section) unblocks both.[^6_2]
- **Watch the RJ palette against the paper.** ABMRS's "neat things" in its final form mostly live in proposals; your palette already exceeds the standard BMARS birth/death/change baseline with prune/grow/variable-change factor moves, `BassFactorShares`, local windowed knot steps, and Wallenius-biased variable selection conditioned on current appearances — you'd be adding the multivariate paper's moves *onto* a strictly stronger base.[^6_4]


## The synthesis differentiator

If the aim is papers-or-better, the edge over every reference implementation I know of isn't any single model — it's the audit layer wrapped around identical math: evidence with Monte Carlo error on three independent estimators (ladder + reference rung, SAMC/WHAM, SMC log normalizer) that can be cross-checked against each other, refusal/support tallies making numerical failure first-class, exact-reference validation available off the shelf (`ExactChangepoint` + `SpectralNoiseField` for recovery experiments the papers did by hand), and conformal on top when the version space settles. I'd strongly suggest making that three-way evidence cross-check an explicit diagnostic — disagreement between estimators is how you catch the subtle RJ/prior bugs that the published repos never detect.[^6_5][^6_3]

One thing I still can't map: what does KHAOS contribute that GBASS's generalized-response axis doesn't? Knowing whether its axis is dynamical/streaming, calibration-flavored, or structure-prior-side would let me say concretely which pieces of the current tree it reuses versus what it makes you build.

---

# KHAOS has some things about polynomials i think, its very much adjacent to GBASS

Found it — and your memory was right. KHAOS is Rumsey, Francom, Gibson, Tucker \& Huerta's **Bayesian Adaptive Polynomial Chaos Expansions** (Stat, 2026): the BASS RJMCMC machinery re-based on products of shifted standardized Legendre polynomials instead of hinge splines. It's genuinely adjacent — same skeleton, different basis atoms — which means for your synthesis it's less a new model than a new `BassBasis` family with a couple of prior wrinkles worth understanding.[^7_1][^7_2]

## What the paper actually does

- **Atom:** $f(x) = \sum_m a_m \prod_i \psi_{d_{im}}(x_i)$, $\psi_d$ shifted-standardized Legendre, orthonormal under uniform input measure on. Terms are degree-partitions across variables rather than knot/sign partitions.[^7_0][^7_3][^7_1]
- **Palette:** four moves — birth, death, mutation(degree), mutation(variable). Notice: your `BassMoveKind` set (Birth, Death, Redraw, LocalStep, Toggle, Prune, Grow, VariableChange) is already a superset in spirit; degree mutation is your `LocalStep`'s analog (bounded integer walk on the degree) and variable mutation is your `VariableChange`.[^7_4][^7_1]
- **Two innovations worth lifting:** (1) a *modified g-prior* whose shrinkage scales with basis-function complexity (degree/interaction order), implemented via Laplace approximation for tuning-free inference; (2) a global regularizer $g_0^2$ sampled by Metropolis–Hastings with an inverse-gamma Laplace-fit proposal. They also offer a fixed ridge alternative, which wins noise-free and loses high-noise — the g-prior's data-driven shrinkage is the robustness mechanism.[^7_1]
- **Neighbor knowledge discovery:** for many-p birth proposals they cite Nott, Kuk \& Duc (NKD) — but your Wallenius-appearance-weighted variable draw is already a strictly adaptive version of that idea (learns during the run rather than fixed a priori).[^7_3][^7_4]


## What changes in your engine — and what becomes free

**The basis swap is small; the payoff lands in sensitivity.** A `PolynomialFactor(Variable, Degree)` next to `HingeFactor(Variable, Sign, Knot)`, `CaseColumn` evaluating per-variable Legendre recurrences and multiplying across factors — contained work. But orthonormality under the input measure means Sobol indices become sums of squared coefficients over term subsets: your entire `BassBasisSensitivity`/`SobolDecomposition`/`BassMoments` apparatus, which currently has to integrate hinge products against the measure, collapses to closed-form bookkeeping over interned terms. That's the honest "better than the package" angle: they get Sobol cheaply *and approximately*; you'd get it exactly from `BassDraw` coefficients you already store per draw.[^7_4]

**The complexity-shrinkage g-prior already has a home in your design.** Their modified g-prior is per-column ridge structure by another name — you have `RidgeForm` with per-column ridge terms carried through `CollapsedRidgePrior.CentredRidgeMarginalAt`, plus `PriorPencil` for per-term precision structure. Complexity-based shrinkage = construct the pencil/ridge terms from (degree partition, interaction order) at proposal time. No new marginal math, just a new ridge-term generator.[^7_4]

**Their $g_0^2$ MH-with-Laplace-proposal is your weakest-link candidate.** Your equivalent global scale is τ on a slice step reading the cached negative-binomial collection weights — already better-integrated. If you implement their exact form for replication purposes, an independence-MH against Inv-Gamma fits inside one `IBlock` using your existing `NewtonMetropolis` and `InverseGamma`; but I'd keep the slice as the default.[^7_5][^7_4]

**Two things PCE gives up vs. hinges, worth surfacing in the synthesis:** no knots means no knot-placement adaptation at all — locality comes only from degree, so discontinuous or sharply-localized surfaces need high degree against the complexity shrinkage (this is exactly their high-noise result in reverse); and high-degree univariate factors have Runge-ish edge behavior that hinge bases don't. The interesting synthesis move available to you and to nobody with the reference code: **mixed bases** — hinge factors and polynomial factors as separate `BassPart`s under one structure prior, letting the RKHS-free RJ palette choose atoms per term. Your part system (Continuous/Categorical/Functional with per-part caps and minimum orders) is architecturally one array slot away from that.[^7_4]

## Sequencing implication

~~KHAOS being basis-local changes the ordering calculus slightly: it doesn't need the `IBassObservationModel` seam widening at all (Gaussian response, g-prior collapse as-is, Laplace only for the $g_0^2$ conditional), whereas GBASS is entirely gated on that seam. If you want a fast win after conformal — one new basis family plus a ridge-term generator, both with closed-form test surfaces available — KHAOS before GBASS de-risks the observation-model refactor by not being coupled to it. It also validates your sensitivity layer against an identical functional (hinge and Legendre Sobol indices should agree on smooth test functions), which is a validation route even for the mixed-basis ambition later.~~


[^7_1]: https://arxiv.org/html/2510.25036

[^7_2]: https://cran.r-project.org/web/packages/khaos/khaos.pdf

[^7_3]: https://github.com/knrumsey/khaos

[^7_4]: src_s013_maths.txt

[^7_5]: src_s003_maths.txt

[^7_6]: src_s004_maths.txt

[^7_7]: https://cran.rstudio.com/web/packages/BASS/BASS.pdf

[^7_8]: https://onlinelibrary.wiley.com/doi/pdf/10.1002/sta4.70151

[^7_9]: https://steamcommunity.com/

[^7_10]: https://dfrancom.r-universe.dev/BASS

[^7_11]: https://steamcommunity.com/sharedfiles/filedetails/?id=209000244

[^7_12]: https://publications.artsci.wustl.edu/web/packages/BASS/BASS.pdf

[^7_13]: https://github.com/lanl/BASS

[^7_14]: https://search.r-project.org/CRAN/refmans/BASS/html/bass.html

[^7_15]: https://steamcommunity.com/login/home

[^7_16]: https://steamcommunity.com/id/officials1mple

[^7_17]: https://steamcommunity.com/id/karrigancsgo

[^7_18]: https://steamcommunity.com/dev/managegameservers

