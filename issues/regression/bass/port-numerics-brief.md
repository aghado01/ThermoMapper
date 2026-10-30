# BASS port: shared numerics and sampler contracts

A segment of the [full port brief](../bass-full-port-brief.md), which holds the scope and completion contract, the dependency boundaries, the stage table and the decision register. This segment covers the BASS-specific parts of stages 2 and 3, on top of the slices the [preparation brief](../bass-preparation-brief.md) landed. The generic SMC core is specified with its consumers in the [calibration and SMC segment](port-calibration-brief.md). Its checklist, where the reports go as work lands, is [port-numerics-checklist.md](port-numerics-checklist.md).

The decision drafts here are proposals for the user's design pass, not decisions. The register, with each decision's status, and the table of sources the drafts were drawn from are in the hub's [decisions section](../bass-full-port-brief.md#statistical-decisions-to-settle).

## Reuse and shared renovations

| Existing surface | Treatment | Concrete work | Status (2026-10-25) |
|---|---|---|---|
| [ReversibleJumpChain and proposal machinery](../../../src/maths/samplers/rjmcmc/ReversibleJumpChain.cs) | Reuse with targeted expansion | Preserve proposal/Jacobian/reverse-selection accounting and telemetry. BASS supplies its own state/moves. Make cache refresh and state replacement safe after hyperparameter or target-stage changes. | Done: steps 0.2 (target refresh in place) and 3.1 (`SetState`, `ReplaceTarget`). |
| `IRjTarget` | Generalize | Move its model-neutral log-density responsibility into a neutral sampler target contract usable by RJ, tempering, and SMC. Migrate callers cleanly. | Done: step 3.1 replaced it with [`ISamplingTarget<TState>`](../../../src/maths/samplers/targets/SamplingTarget.cs). |
| [ParallelTempering](../../../src/maths/samplers/rjmcmc/ParallelTempering.cs) | Expand and relocate | Replace concrete RJ-chain construction with per-replica kernel factories and explicit state exchange. Put general replica exchange in `samplers/tempering`; retain RJ mechanics in `samplers/rjmcmc`. | Expanded: steps 3.1–3.3 (per-replica move factories, declared tempering paths, the non-reversible swap schedule). Not relocated; it remains in `samplers/rjmcmc/`. |
| [IChain](../../../src/maths/samplers/ensemble/IChain.cs), [IEnsembleModel](../../../src/maths/samplers/ensemble/IEnsembleModel.cs), ensemble execution | Reuse after qualification | A complete BASS sweep presents as a chain. BASS supplies invariant readouts and its rich posterior accumulators. Construct private mutable kernels/workspaces and preserve cumulative statistics across state or target replacement. | Qualified: steps 0.2–0.3 and 3.4. BAPS, BGP and `HierarchicalPSpline` also run on it. |
| [ChainEnsemble](../../../src/maths/samplers/ensemble/ChainEnsemble.cs), [ChainDiagnostics](../../../src/maths/samplers/ensemble/ChainDiagnostics.cs) | Repair in stage 0, then expand | Undefined/infinite R-hat must not count as convergence or be displayed as zero. Qualify minimum-data and constant-functional policies first. Add trace-based rank-normalized split/folded diagnostics and appropriate bulk/tail ESS in the shared renovation stage. Particle diagnostics remain separate. | Done: steps 0.1 and 3.4 (rank-normalized diagnostics, bulk and tail ESS). The external reference is pending (OQ9). |
| [CholeskyDecomposition](../../../src/maths/linalg/Cholesky.cs) | Add a strict factorization path | Explicit status, rank/conditioning information, solves into buffers, bordered insertion, deletion/replacement, and rebuild checks. Do not silently inherit pivot flooring or eager inverse construction. | Strict reference done: step 2.2, [`DenseCholesky`](../../../src/maths/linalg/DenseCholesky.cs). Incremental factor done 2026-10-26: [`BorderedCholesky`](../../../src/maths/linalg/BorderedCholesky.cs), with bordered births, Givens deaths, staged changes, strict rebuilds and diag(A⁻¹) (N1–N4). |
| [MatrixOps](../../../src/maths/linalg/MatrixOps.cs) | Reuse and expand | Reuse dot products and column-oriented Gram operations; add general contractions as needed. Term evaluation and cache semantics stay in BASS. | Pairwise cross-products added 2026-10-26: [`PairwiseSum`](../../../src/maths/linalg/PairwiseSum.cs) (N10). The direct functional contractions too: [`KhatriRao`](../../../src/maths/linalg/KhatriRao.cs) (N7). |
| [PCA](../../../src/maths/linalg/PCA.cs) | Expand | Add stable thin SVD and numerical-rank handling in LinAlg. BASS owns response scaling and its coefficient convention. SVD implementation/provider selection is a numerical design task. | Thin SVD done 2026-10-26: [`ThinSvd`](../../../src/maths/linalg/ThinSvd.cs), Householder QR then one-sided Jacobi, with numerical rank (N5). `Pca` and the Stiefel and Grassmann manifolds moved onto it 2026-10-27. |
| Covariance solve/log-determinant operations | Add in LinAlg | Diagonal-plus-low-rank operations for functional calibration, plus the strict dense path for full covariance. Calibration supplies the covariance components. | Diagonal-plus-low-rank path done 2026-10-26: [`DiagonalPlusLowRank`](../../../src/maths/linalg/DiagonalPlusLowRank.cs) (N6). The dense path is `DenseCholesky`. |
| [Normal](../../../src/maths/distributions/Normal.cs), [Gamma](../../../src/maths/distributions/Gamma.cs), [Beta](../../../src/maths/distributions/Beta.cs), [SpecialFunctions](../../../src/maths/distributions/SpecialFunctions.cs) | Reuse and expand | Add inverse-gamma/truncation support and stable normal/Student interval probabilities and moments. Make rate versus scale explicit. | Done: step 2.1 (`InverseGamma`, `StudentT`, `Truncated`, CDFs and interval probabilities; rate and scale named at construction). |
| [Xoshiro256PlusPlus](../../../src/maths/rng/Xoshiro256PlusPlus.cs), [SeedTree](../../../src/maths/rng/SeedTree.cs) | Reuse and expand where needed | Preserve RNG capture/restore. Supply reproducible independent offspring streams for resampling and stable component/replica identities. | Unchanged; offspring streams for resampling wait for the SMC core. |

The generalized sampling contracts affect existing `SplineTarget`, `BarsEnsemble`, and sampler tests. First record the intended statistical corrections in the existing-engine repair track. Then migrate callers against that qualified baseline without further incidental model changes. Do not preserve a known incorrect distribution merely to keep old expectations green. During refactoring, follow the repository's staged replacement and final cleanup discipline; the end state has no compatibility aliases.

Several existing components provide useful patterns but cannot be substituted unchanged:

- [IBasis](../../../src/maths/regression/spline/IBasis.cs) and `KnotConfig` describe one-dimensional ordered knots. Preserve that meaning; BASS needs multivariate product terms.
- [BAPS TensorProductDesign](../../../src/maths/regression/spline/baps/TensorProductDesign.cs) constructs a fixed bivariate B-spline basis. It is not BASS's changing interaction basis.
- [BandCholesky](../../../src/maths/linalg/BandCholesky.cs) provides sampling/factorization patterns, but BASS's active-term Gram is generally dense.
- [WeightedNormalModel](../../../src/maths/regression/spline/bars/WeightedNormalModel.cs) has a different coefficient prior and marginal target. Its numerical operations are reusable; its statistical formula is not the BASS formula.
- [IObservationModel.SampleConditional](../../../src/maths/regression/spline/bars/IObservationModel.cs) defaults to the conditional mean with a NaN σ², which is not a full posterior draw. Weighted Normal draws (σ², β) exactly (step 0.5c), Poisson and Whittle approximately, and Gibbs loss returns the mean; the model's `CoefficientUncertainty` label says which. Reuse its extrema intervals as an uncertainty reference only where the label is `Exact`.
- BARS's fixed complexity-prior helpers do not implement BASS's sampled count-intensity hierarchy.

Existing covariance callers may legitimately require regularization or inverse output. A new strict factor path must not silently change their semantics. No blanket conversion of the existing engines from one dimension to two dimensions is required.

### Shared Gibbs building blocks (OQ18)

**Settled: shared Gibbs building blocks, and separating samplers from model classes** (opened 2026-10-23; the plan completed 2026-10-24, phases A–D below). Tracked as the preparation brief's OQ18. The user's direction: this kind of unification is the ideal, and samplers should be separate from model classes where that is not unwieldy.
- **The observation.** The engines hand-roll the same conjugate steps: 21 sites across 7 files in `Maths.Regression`.
  - **Banded Gaussian conditional draws:** factor a precision band, solve for the mean, add a `SampleInnovation` draw. They appear in BAPS's Gibbs and adaptive sweeps and its marginal-mode β, in `HierarchicalPSpline`'s β_r and β₀, and in the g-prior block's coefficient draw.
  - **Inverse-gamma updates from a sum of squares:** σ², τ², τ_u², τ₀², and the σ² behind BARS's robust weights.
  - **Gamma updates for local multipliers:** BAPS's η, and the Student-t weights.
  - **Random-walk Metropolis on a scalar:** BAPS's log λ, and BGP's log t.

  BASS needs all of them, for its coefficient, variance, hyperparameter and count-intensity blocks.
- **The shape, to be designed (not decided):**
  - **Primitives in `Maths.Samplers`:** a Gaussian conditional draw from a precision factor and right-hand side, the inverse-gamma and gamma conjugate updates, and a scalar random-walk Metropolis step that counts its acceptances. `Maths.Samplers` already references LinAlg and Distributions, so they fit there.
  - **A generic scan in `Maths.Samplers`:** a chain over a model's state that runs its blocks in order each step, counts acceptance per block, and records a draw. It saves its streams and state for checkpoints, and the ensemble shell drives it like any other chain.
  - **Composition with the reversible-jump chain:** BARS's step is already a reversible-jump step plus the weight block, assembled ad hoc in `StepOnce`, and BASS is reversible-jump moves plus conjugate blocks. One composite kernel would serve both.
  - **What stays with the model in `Maths.Regression`:** the state type, each block's sufficient statistics (residual sum of squares, roughness, the penalty band), and the result assembly.
- **Why `BapsSampler` cannot simply move to `Maths.Samplers`.** Its conditionals use `PenalizedSpline`, `DifferencePenalty` and `BandedDesign`, all in `Maths.Regression`, and Regression references Samplers, so the move would make a cycle. The split above is the version that works: the generic scan and primitives move down, and the model's blocks stay. It matches the reversible-jump engine's existing split, with `ReversibleJumpChain` in Samplers and `SplineTarget` and the knot moves in Regression.
- **The naming plan (decided by the user on 2026-10-23; applied 2026-10-24 in phase D, `21f6c04`, `6880dcf`, `e9e021b`).** The user: "baps is the method, something that sits above just the sampler."
  - **BARS, BAPS and BASS are method names.** They name the folder, the namespace, the entry point and the result, and they stay.
  - **The samplers inside a method are named for their algorithm in the field.** `BapsSampler` is really the BAPS method with two samplers embedded:
    - **block Gibbs**, for the `Gibbs` and `AdaptiveGibbs` schemes: a conjugate block draw of β, then inverse-gamma variance updates, plus η in the adaptive scheme;
    - **collapsed Metropolis**, for the `MarginalEvidence` scheme: random-walk Metropolis on log λ with β and σ² integrated out, then the exact conditional draws of σ² and β.

  | Level | Today | Planned | Note |
  |---|---|---|---|
  | Method | `BapsSampler` | `Baps` | The entry point: data, priors, scheme, result assembly; parallel to `Bars` |
  | Method | `BapsResult`, `BapsDraw` | unchanged | |
  | Method | `BapsLambdaUpdate { Gibbs, MarginalEvidence, AdaptiveGibbs }` | `BapsScheme { Gibbs, Collapsed, LocallyAdaptiveGibbs }` | "λ update" does not describe the adaptive scheme, which changes the model (local weights on the penalty), not the update |
  | Sampler | inside the private `BapsChainRun` | the shared block Gibbs scan in `Maths.Samplers` | BAPS's blocks (β, σ², τ², η) stay in `baps/` |
  | Sampler | inside the private `BapsChainRun` | collapsed Metropolis | A scalar random-walk Metropolis block (BGP's log t uses the same step), followed by the conditional draws |

  - **Canonical references,** for the docs when the names land: Bayesian P-splines (Lang & Brezger 2004), with block updates of the coefficients as banded Gaussian Markov random field draws (Rue 2001). For locally adaptive smoothing: Lang & Brezger 2004, and Baladandayuthapani, Mallick & Carroll 2005. For the uniform prior on effective degrees of freedom: Hodges & Sargent 2001.
  - **`HierarchicalPSpline`** is already named for what it is. Its sweep becomes a block Gibbs scan in the same way.
- **The check when it lands:** each engine moved onto the blocks must be bit-identical to its current sampler, with the same streams and the same order, as BAPS's move onto the ensemble shell was. The renames ride along in the same change.
- **Sequencing:** `HierarchicalPSpline`'s move onto the ensemble shell waits for this, so it lands on the shared scan rather than another hand-rolled chain class. Its draw retention does not depend on it.
- **Phase A landed (2026-10-24, `c525cd4`, `fb077f3`):** `GaussianConditional` (over any `ISpdFactor`, returning the conditional mean beside the draw), `ConjugateUpdates` (inverse-gamma and gamma posteriors of a normal scale) and `BlockScan<T>` (ordered blocks, per-block cadence and counters, `IChain`, exact snapshots) in `Maths.Samplers.Blocks`; `HierarchicalPSpline`'s sweep runs on them as five blocks, bit-identical to its former loop. The scalar random-walk Metropolis step waits for BAPS's λ and BGP's t, whose acceptance draws consume the stream differently (BAPS always draws the uniform; BGP skips it when the ratio is at least 1). Next: `HierarchicalPSpline` onto the ensemble shell (phase B), the mixing study (C), then BAPS and BARS's weight block (D).
- **Phase B landed (2026-10-24, `d2db820`):** `HierarchicalPSpline` runs on the ensemble shell through the block scan: overdispersed starts by default (the original fixed start kept), rank-normalized diagnostics on the three log variances, a classic-R̂ consensus map on the population curve, the early stop and exact checkpoints. Its diagnostics showed τ_u² mixing slowly at once, which phase C takes up.
- **Phase C landed (2026-10-24, `d59a0c4`):** `HierarchicalScheme { Centred, CollapsedPopulation, Interweaved, ParameterExpanded }` and `HalfTDeviationPrior`, each a different block list on the same scan. The collapsed block uses the dense path of `GaussianConditional`, and the interweaving's non-centred step is a Metropolis block. The block scan composed all of them without change.
- **Phase D landed (2026-10-24, `21f6c04`, `6880dcf`, `e9e021b`):** `RandomWalkMetropolis` in `Maths.Samplers.Blocks`; BAPS renamed (`Baps`, `BapsScheme { Gibbs, Collapsed, LocallyAdaptiveGibbs }`) with its Gibbs schemes as a block Gibbs scan and its collapsed scheme as a Metropolis block; BARS's step as a scan of the knot step and the weight block on its cadence. All bit-identical. BGP's walk keeps its own step, since its acceptance skips the uniform when the ratio is at least 1. OQ18's plan is complete; BASS's blocks are the next clients.

## Dense factor maintenance

Establish a strict full-factorization implementation as the reference. Then add bordered births, Givens-based deletions, and initially a staged delete/insert for a changed term. A border extension changes dimension; it is not a same-size rank-one update. A replaced Gram row and column has a symmetric change of rank at most two.

Use relative finite/rank/conditioning checks and periodic or triggered full rebuilds. A fixed absolute `1e-12` threshold alone is not a scale-aware admission policy. When an incremental update fails because of accumulated numerical error, check a strict rebuild before deciding that the candidate is inadmissible. Singular g-prior designs cannot be made mathematically defined merely by hiding a ridge inside the factorization. A different proper coefficient prior may admit a rank-deficient design; LinAlg must not impose a universal statistical rejection policy.

Data cross-products, the Gram's entries, are summed pairwise or with compensation, not left to right. On a near-constant column (a categorical indicator true for all but a few cases), left-to-right summation adds many equal centred values after one large one, and they round the same way every time. In D06's fixture F-O7 that cost 10⁻³ relative in the admission angle near the gate, where pairwise summation reached ε·κ. Under D03's Liang convention the strict rebuild also forms the centred Gram without the rank-one downdate (D06).

Candidate numerical work remains private until acceptance. Rejection leaves the accepted state unchanged. Coefficients are drawn using precision-factor triangular solves, not a materialized inverse. Factor maintenance can cost `O(M²)` while term evaluation and data cross-products retain separate data-dependent costs.

## D10 (note, adopted 2026-10-27): a policy for adaptive kernels

- **BASS's own proposals are not adaptive.** I(K | M) and the variable weights follow the current model's terms, not the chain's history, so they are ordinary Metropolis–Hastings proposals that depend on the state. They need exact, positive reverse probabilities (D01), not adaptation theory.
- **Adaptive kernels follow one policy.** These include mutation-type weights tuned from acceptance rates (as KHAOS does), tempering ladders, and a history-based variable proposal. The policy follows Nott and Kohn (2005, *Biometrika* 92, Theorems 1–2, for finite and compact state spaces) and Roberts and Rosenthal (2007, for general spaces):
  - tune only from the past, so that each step's kernel leaves the posterior invariant;
  - keep a floor on every component, a uniform minorization such as KHAOS's 10% floor on each mutation type;
  - let the adaptation diminish (changes of order 1/n, as running averages give) or freeze it after burn-in, the simplest sufficient route and the default;
  - store the adaptation state in checkpoints and replay it on restart.

  Adaptation that continues past burn-in leaves the diagnostics' Markov assumption only approximately true, which is a further reason freezing is the default.
- **A candidate kernel:** a birth proposal that weights variables by inclusion frequencies accumulated over the run, aimed at Nott and Kohn's setting of many inert predictors. It is an option to measure under this policy, not a default ([extensions brief](extensions-brief.md#additional-proposal-kernels)).

## D10 (decided 2026-10-29): stopping and convergence semantics

Drafted and decided on 2026-10-29: the user took the draft as written, with its four points. Built the same day by N14 (`e8547d3`, `e88b33e`, `1c26dcd`; the numerics checklist's report). Sources:
- Vehtari, Gelman, Simpson, Carpenter and Bürkner (2021), on rank-normalized R̂ and bulk and tail ESS;
- B01's repaired stopping mechanics (`ChainEnsemble`'s R̂ stop);
- F-M2's pass (2026-10-29, the fitting checklist's reports), whose untempered chains never agreed at S10's lengths;
- argon's notes (`issues/regression/argon-bass-qa.md`, 2026-10-29), weighed in the main session.

- **The gap.** `ChainEnsemble` checks its R̂ stop between rounds, never after the last.
  - `StoppedEarly` is false both when no target is set and when a target runs out of budget, so a fit never says it did not converge.
  - On untempered Friedman, M's R̂ near 2 shows only to a reader of the per-quantity diagnostics.
- **A verdict on every run, in the shared ensemble.** `EnsembleRun` gains a convergence verdict, read at the end of every run whether or not a target was set, by the rules the stop already uses. BARS and BASS read it alike, which is why D10 sits in this segment.
  - *The target:* `RHatTarget` when set; otherwise 1.01, Vehtari et al.'s threshold for rank-normalized R̂.
  - *The ESS gate:* `MinEssPerChain`, 100 per chain in the bulk and the tails, as now.
  - *B01's rules, unchanged:* a frozen or unavailable reading fails, and a constant functional passes but cannot certify the verdict alone.
  - *Three outcomes:* converged; not converged, naming each failing quantity with its R̂ and ESS; and undetermined, when nothing could be estimated (one chain, or every functional constant).
  - *A fit is always returned.* Not converging is a reported state, never an exception or a null, and the checkpoint lets the caller extend the run (`Bass.Resume`).
- **Three agreement groups.** The engine declares the grouping of its monitored quantities; the verdict is generic. For BASS:
  - *predictive:* σ² and the fitted values at the monitored points. Do the chains predict the same surface and noise?
  - *posterior level:* each retained draw's log posterior, the collapsed log π(structure, u | y) at b = 1, as F-M2's pass computed it. Do the chains sit in regions of the same mass? This draft adds it to the monitored quantities: the pass found it set the minimum ESS in many fits.
  - *structural:* M, and τ in the g mode. Do the chains agree on the representation? In an overcomplete hinge basis M is not identified by the surface. So structural disagreement beside predictive and level agreement means non-identifiability, and the verdict reports it rather than hiding it.
  - Each group gets its own verdict, and the run's verdict is their conjunction by default.
  - *The pass, read this way.*
    - On untempered Friedman every group fails: the fitted values reach R̂ 1.84, and the log posterior and M sit near 2.
    - On the mixed and functional fits σ² agrees (R̂ 1.00–1.01), but one monitored point does not (1.45 and 1.48), and the log posterior and M do not. That is genuine disagreement, not only non-identifiability.
    - Tempered, every group's median R̂ is at most 1.38 and no fit reads above 1.61, still short of 1.01.
- **The stop's scope** (`StoppingScope` on `BassSampling`; in the ensemble, a set of groups): which groups an R̂ target gates.
  - The default is every group.
  - The option is predictive only, for an emulator user who needs the surface and the noise and accepts a reported structural disagreement. The verdict still reports every group.
  - A predictive-only stop can end a run whose chains sit at different posterior levels, so it is a declared choice, not the default.
- **Run lengths are a budget.** `BurnIn` and `Samples` stay R's 9,000 + 1,000, provisional, as the most a run spends. An R̂ target may stop sooner, and the verdict says whether the budget sufficed.
  - F-M2's pass gives no basis for sizing them from ESS per sweep: no arm cleared, and every untempered fit read the ESS floor.
  - So they stay R's until a sizing has a converged reading to start from.
- **Mode-separation readings,** reported beside R̂. Once chains sit in separate modes, R̂ and bulk ESS saturate, with the ESS clamping near the split chains' count over two. These readings measure how far apart the modes are:
  - the gap between the chains' mean log posteriors, the largest minus the smallest (12–42 nats in Phase 0's M4);
  - the largest distance between two chains' posterior-mean surfaces at the monitored points.
  - The verdict reports them. No stop reads them in this draft, since they have no calibrated threshold yet.
- **Friedman's worse modes and stage 8.** Move work on them is parked.
  - A mode 12–42 nats below the best carries about e⁻¹² of its mass or less, so a correct sampler all but never sits there. *(Corrected 2026-10-29.)* That is a gap in density per structure, not in mass. A basin's mass sums over the structures it holds, and basins differ in how many they hold, so the mass ratio needs basin fractions or per-basin evidence. The Opus chat's read in `issues/regression/argon-friedman-etc.md` caught it.
  - SMC's resampling along D05's path weights particles by those ratios, which MCMC chains without round trips cannot. Stage 8's fitting adapter (C4) is the test, recorded as a qualification question: does a population of a given size find the better representation, and keep it?
  - *The caveat:* births are bounded by the Occam cost, so some particles must reach that representation at low b. The particle count is the variable to read.
  - Any later move work reads the agreement time from dispersed starts: the sweeps or seconds until every chain's mean log posterior lies within Δ nats of the leader's. It does not read ESS per second.
  - *Answered 2026-10-29* (C4, then the Friedman measurement, [its report](friedman-band-brief.md#reports)): no population read, up to 16,000 particles, finds and keeps one answer. The band is a first-order transition, its two phases coexisting at fixed b, and the lineages collapse through it at every N.
  - *Then the reference path* (2026-10-29, [its report](friedman-reference-brief.md#reports)): a population started on the pilots' structures agrees within a nat at 4,000 particles, but on those structures' own evidence, 22–23 nats below the 99.9% lower bound that R2's from-prior normalizers set on ln p(y). Agreement read from such a start measures the reference's coverage, not the evidence.
- **Not adopted:**
  - argon's 1.05 as the reference target, which is looser than the threshold the shared R̂ was built against; a caller can still declare 1.05;
  - a stop on the mode-separation readings, until they have a threshold;
  - an exception on non-convergence: a fit that can be inspected and extended serves better.
- **What it changes:**
  - `EnsembleRun` gains the verdict, and `BassDiagnostics` carries it by group; BARS's diagnostics display it under B01's display policy.
  - Adding the log posterior to BASS's monitored quantities enters the R̂ stop's gate, so a run with a target may stop later. A run without a target is unaffected, and so are the draws.
- **Independent checks:**
  - a synthetic ensemble whose chains sit at two levels: the verdict fails the right groups and names the failing quantities, and constant and frozen functionals follow B01's rules;
  - the end-of-run verdict agrees with the stop's rule applied to the same traces;
  - the log posterior against the chains' own density, as F-M2's pass checked it (to 10⁻¹¹);
  - bit-identity: the verdict reads the run and changes nothing in it, so the existing hashes hold.
- **Decided with it** (the user, 2026-10-29):
  - the reference target, 1.01;
  - the log posterior added to the monitored quantities;
  - every group gating the stop by default;
  - R's budget as the provisional default.

## Room for extensions

Prep work in this segment, so that the [extensions](extensions-brief.md) plug in without reworking the port:
- **One face-splitting contraction.** The direct functional algebra (fitting segment) needs the column-wise Khatri–Rao product A ⊙ H, with the Hadamard Gram G = (AᵀA)∘(HᵀH) when the samples form a Cartesian grid. Its columns are a_j ⊗ h_j, its rows the Hadamard products of a case's and a location's factor rows, and its transpose a face-splitting product. ABMRS's rows are face-splitting (row-wise Kronecker) products themselves. Built as [`KhatriRao`](../../../src/maths/linalg/KhatriRao.cs) (N7, 2026-10-26). Build it once in LinAlg with that shortcut and a sparse accumulation for scattered samples. ABMRS's tensor designs reuse it: a Kronecker Gram on gridded inputs, and accumulation at O(n(p+1)^(2d)) on scattered ones.
- **A factor contract open to a wide-band factor.** BASS uses the dense path. ABMRS's Gram is block-banded, so `ISpdFactor` should admit a wide-band factor beside `DenseCholesky`, chosen by size.
- **Conjugate updates as limits of one form.** Keep `ConjugateUpdates`' Gamma and inverse-gamma updates expressible as limits of a GIG update, so GBMARS's latent-scale block generalizes them rather than sitting beside them. No GIG code is needed now.
- **The variance-inflation vector.** The diag(G⁻¹) that D06's order-free gate maintains also gives every coefficient's variance inflation, which an optional `MaxVarianceInflation` cap reads.
- **D10's policy** covers the extension kernels: the history-based variable proposal and adaptive mutation weights.
