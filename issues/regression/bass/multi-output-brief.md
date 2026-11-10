# Several outputs on one basis: E34's build (chip brief)

Date: 2026-11-08. Status: commissioned by the user on 2026-11-08 ("Yes, start E34"), after E8 was completed; slice 0 starts. Written by the main session from [E34's entry](extensions-brief.md#several-outputs-on-one-basis-with-σ-integrated-e34), which holds the model, its pieces and the checks. Reports go under "Reports" below, newest first.

## Why

K continuous outputs per case that share one structure: several responses of one experiment, unordered or few enough that a free K × K covariance is affordable. One fit selects terms for all outputs at once, reads their residual correlation, and gives joint predictions:
- joint bands;
- the cross-output correlation of predictions;
- the probability that every output meets its specification.

It serves any multi-response surrogate or regression, and E22 (coupled outcomes) builds on it. Ordered outputs on a grid stay with the direct and PCA modes.

## Sources

- **The model** is multivariate Bayesian regression with a matrix-normal prior scaled by Σ and an inverse-Wishart on Σ. Brown, Vannucci and Fearn (JRSS B, 1998) select variables in it with inclusion shared across responses; Huang and Wand (*Bayesian Analysis*, 2013) give the scale-mixture alternative to the inverse-Wishart. Both are in `supellex-private\PDFs\Bayesian\` (`jrssb98.pdf`, `HuangWand2013.pdf`), added by the user on 2026-11-08 and read in slice 0 (its addendum below). BVF write Σ ∼ IW(δ; Q) in Dawid's shape, δ = ν₀ − K + 1.
- **No oracle** is known in CRAN's lineage. The checks are internal and against closed forms, as the entry lists.

## What is already built

- **The collapsed core:** `CollapsedGPrior`, `GPriorConditional`, `PriorPencil`. Its coefficient prior runs through the prior-precision builder (`BassPriorPrecision`: the g-form, the normalized ridge, KHAOS's raw ridge and modified g-prior), and σ² is integrated under `Sigma2Prior` (`ScaledToResponse` by default, D03).
- **Responses:** `BassData` and `BassResponse` (`ScalarResponse`, and a functional response over shared output locations), plus basis fits through `BassDecomposition` and `BassComponents` (D07: PCA, each component fitted apart).
- **The structure chain and its samplers:** the ladder with `RungWorkers`, SMC with its two flags, the flat histogram, resume.
- **GBASS's noises:** the latent record per observation (`ObservationLatents`), weights and offset on the working response, and the move along β·v as a partially collapsed step (`1a8459c`: the coefficients redrawn after an accepted move).
- **Random intercepts** (E21, `0a2168c`): the rows whitened by V = I + s²ZZᵀ, with the effects drawn from their conditional.
- **Readouts:** the moment service and Sobol; the gradient readout; E8's co-active posterior for two outputs of one fit (`CoActivePosterior`).

## Decisions

These are the main session's, taken in commissioning. Each goes on the page with its alternatives and trade-offs. The user may reverse any of them.

1. **The model.** Y (n × K) = 1μᵀ + BA + E, with each row of E ~ N(0, Σ).
   - *The priors:* A | Σ, τ ~ MN(0, S₀(τ), Σ) through the prior-precision builder, and Σ ~ IW(ν₀, Ψ₀). Inclusion is shared across outputs.
   - *The collapse:* with A and Σ integrated, the structure step reads a K × K quadratic form and its log-determinant, and the marginal is matrix-t. The intercept is level-centred per output, inside the prior as scalar BASS's is, so Σ's posterior counts n observations, nb at rung b; the flat intercept's counts n − 1, nb − 1 at rung b. *Corrected 2026-11-08* from "n − 1 degrees of freedom" by slice 0's reading: the K = 1 reduction to scalar BASS and Brown, Vannucci and Fearn's eq. (20) both count n.
2. **Σ's default prior: D03's `ScaledToResponse`, carried to K outputs.**
   - *The form:* ν₀ = K − 1 + 2c and Ψ₀ = 2c·S_Y, with c = 10⁻³, stated as an empirical prior. S_Y is scalar BASS's s² carried to K outputs: the outputs' mean squares and cross-products about the intercept's centre levels, with its divisor n. Those levels are ȳ under level centring and the flat intercept, and 0 under the zero centre. *Corrected 2026-11-08* from "sample covariance about their means" by slice 0.
   - *At K = 1* it is IG(c, c·s²), the scalar default exactly. That reduction is a gate.
   - *Units:* because S_Y → LᵀS_YL when Y → YL, the structure posterior is invariant to any invertible rescaling or mixing of the outputs. This settles the entry's output-units choice.
   - *Its correlations:* under a diagonal S_Y each prior correlation has density ∝ (1 − ρ²)^(c − 1), Huang and Wand's Property 3 at their ν = 2c, so at c = 10⁻³ nearly all its mass sits at ±1. The data overrule it at moderate n/K. It stays the default provisionally, for the reduction and the invariance; slice 3 measures it at small n/K.
   - *Options:* a declared Ψ₀ in the outputs' units; standardized outputs (scale-invariant only); and Huang and Wand's prior at ν = 2.
     - *Huang and Wand's:* Σ | a ∼ IW(ν + K − 1, 2ν diag(1/a)) with a_k ∼ IG(½, 1/A_k²), A_k each output's root mean square about its level. It is the engine's hierarchical σ² prior carried to K rates, each moved by a slice step as τ is, and the collapse holds given them. It gives uniform correlations, and it is equivariant under rescaling and permuting the outputs, not under mixing them.
     - *When:* slice 2, since its rates are state coordinates that every sampler and resume carry, as the scalar rate is (added 2026-11-08 by the main session; the slice the chip's choice).
   - *Refused:* an S_Y that is singular or nearly so, with the reason. Collinear outputs need fewer outputs.
3. **Properness under tempering, correcting the entry.** At rung b the tempered target integrates p(Y | A, Σ)^b against the priors (D05's integrated form), so Σ's conditional is IW(ν₀ + ν_b, Ψ₀ + R_b). Here ν_b = nb under level centring and the zero centre, and nb − 1 under the flat intercept (*corrected 2026-11-08* by slice 0 from "b·ν_data with ν_data = n − 1").
   - *When it is proper:* whenever ν₀ + ν_b > K − 1. Under a proper prior (ν₀ > K − 1) with the intercept inside the prior, that holds at every b ≥ 0, including the reference rung. Under the flat intercept it fails at nb ≤ K − ν₀ under a proper prior too (nb ≤ 1 − 2c under decision 2), and b = 0 diverges, as scalar BASS's flat pole does.
   - *What changes from the entry:* it asked the tilted prior alone to be proper (ν₀ − n(1 − b) > K − 1), which is stricter than needed. Refuse only a rung, or a reference rung, where the condition above fails: under an improper ν₀, or at the flat intercept's hottest rungs.
4. **Shared inclusion only.** Inclusion per output needs an indicator per (term, output) and loses the matrix-t.
5. **Missing outputs are refused** at first, with the reason.
6. **GBASS's noise by row.** A per-row scale vᵢ is shared across a case's outputs, with a skew vector γ ∈ ℝᴷ, giving the multivariate t and the multivariate generalized hyperbolic.
   - *The γ·v move* is the β·v move's analogue: a partially collapsed step, with A and Σ drawn afresh after an accepted move, and a proposal that reads the state alone.
   - *What carries over:* the skewed horseshoe's default move and its recorded hot-rung limit.
7. **Random intercepts with one s shared across outputs.** Group effects u_g ∈ ℝᴷ ~ N(0, s²Σ), so the row covariance is V ⊗ Σ, the matrix-normal form holds, and the effects' correlation is Σ's. An s per output breaks the form and is not built.
8. **The PCA mode's joint option.** D07's retained scores are fitted as K outputs here, with truncation residuals as now. It depends only on the scores' span (decision 2), and trades per-component sparsity for that.

## The slices

One commit per concern; tests with the code they test; each slice ends with its gates and a fresh-eyes review before the next starts.

| Slice | Delivers | Its own gates |
| --- | --- | --- |
| 0. References | The matrix-t marginal in closed form; the sources read once supplied | Against the dense vec form with an asymmetric Σ at K = 2, so Σ ⊗ S₀ and S₀ ⊗ Σ cannot pass for each other under the code's vec order; against quadrature over Σ at K = 2 and Monte Carlo at K = 3; at K = 1 against the scalar marginal under decision 2's reduction |
| 1. The core | A multi-output response; Σ's prior (decision 2) and its validation; the collapse in the structure step and τ's step through `PriorPencil`; per draw Σ | Y and A | Σ; level-centred intercepts | K = 1 equal to scalar BASS to rounding, with scalar BASS's pinned hashes unmoved; the law by enumeration on small structure spaces; Geweke; equivariance Y → YL; recovery of a known Σ and of surfaces from outside the model |
| 2. Samplers and evidence | The ladder under decision 3, SMC with its flags, the flat histogram, resume; Huang and Wand's prior (decision 2's option), its K rates through the sweep and every sampler; the normalized ridge and KHAOS's two coefficient forms for several outputs, their spectra and pencils carrying K cross products (moved from slice 1, 2026-11-09) | Each rung's law and the evidence against the enumeration; resume bit for bit; a refusal only where decision 3's condition fails; Huang and Wand's prior by enumeration with its rates integrated, and Geweke; each coefficient form's target against the matrix-t reference under its prior covariance, and at K = 1 against scalar BASS |
| 3. Readouts | The matrix-t predictive; joint bands; the predictions' cross-output correlation; P(every output in its specification); E[Σ] and its correlation matrix; Sobol per output and for declared linear combinations; co-active pairs across outputs through E8 | Against Monte Carlo over the draws; the joint probability against counting; Sobol per output against single-output fits on the same structure. *Decision 2's measurement:* the coverage of the correlations and of the joint predictive at n/K from 3 to 30, decision 2 against Huang and Wand at ν = 2. If decision 2 undercovers where Huang and Wand does not, the default goes to the user with the numbers; otherwise decision 2 stands |
| 4. Noise and random intercepts | Decisions 6 and 7 | Geweke under each noise and with the block; the γ·v move on a frozen state against its orbit's law; recovery from outside the model |
| 5. The PCA mode's joint option | Decision 8 | Against separate component fits on data whose scores are independent; the equivariance gate on the scores |

## Gates common to every slice

- **Scalar BASS untouched.** Its pinned hashes and fingerprints hold: the pinned tier at each commit, the whole assembly at the slice's end ([targeted runs](../../../CONTRIBUTING.md#testing-standards)).
- **Truth from outside the model,** on surfaces and covariances not drawn from the engine's priors.
- **Primitives:** chosen cases and a seeded random sweep, far ranges and near-singular Σ included.
- **The cost row:** the structure step at K = 1, 2, 5 and 10 against scalar BASS, and the readouts. Cycle clock, marked indicative.
- **`docs/Maths`:** the decisions on the regression page or a new page, with their alternatives and implications. New writing takes the module's symbols, and crossings go in `issues/notation-ledger.md`.

## Not in this brief

- **E22's coupled outcomes** and its parameter expansion for a correlation matrix.
- **Inclusion per output;** missing outputs as latents.
- **The response envelope** (E34's entry), and an output basis from E[J W_x Jᵀ].
- **Conformal prediction for several outputs.** The candidate seam refuses a multi-output fit, with its reason.

## Working rules

- On `main`, committing by pathspec with nothing left staged and nothing pushed. Commit messages describe the work, with no attribution lines ([CLAUDE.md](../../../CLAUDE.md)). Leave `basin-hold-validation.csv` and `path-profile-validation.csv` out of every commit.
- Other sessions may be writing in the checkout: the conformal chip on E35, the GBASS chip on a fix to `NormalVarianceMeanMixture`. Before editing a file, check `git diff -- <file>` for changes you did not make; build from a snapshot (`git archive HEAD` plus your files); commit only compiling states. Shared docs (`docs/Maths/BayesianRegression.md`, the notation ledger, `extensions-brief.md`) are edited and committed in one step.
- Builds set `ArtifactsRoot` inline into `artifacts\chip-e34\`, deleted at the end.
- Heavy runs at below-normal priority with capped threads, MSBuild node reuse off; the machine is shared. Tests run by [targeted runs](../../../CONTRIBUTING.md#testing-standards) (`scripts/test-select.ps1`, `-SourceRoot` at your snapshot): the selected classes per edit, with the pinned tier per commit, the whole assembly once at the slice's end.
- **Sources:** read a paper before stating what it says, and mark anything cited from memory. When a source is not on disk, ask the main session, which asks the user. Fetch nothing yourself.
- New beside old; a removal is its own question to the user, with what is lost.
- Findings framed by problem class, not by any one application.
- A follow-up goes in a report with the decision it would change and what each outcome would do.
- After each slice, report to the main session and wait for its verification note before starting the next.

## Report

Under "Reports" below, newest first, per slice: what landed by commit; the gates' numbers; decisions taken beyond this brief; what was left unbuilt and why; follow-ups with the decision each would change.

## Reports

### Slice 2: samplers, covariance rates and coefficient forms (chip, 2026-11-09)

**Landed:** `4f3f210` adds the normalized test-side rung and evidence references and the properness boundaries. `41727ac` adds the ladder, SMC, flat histogram and exact resume for several outputs; optional Huang and Wand rates in the state, sweep and every sampler; and the normalized ridge, raw ridge and modified g-form through a shared `OutputPencil`. The math and implementation choices are on [the several-output page](../../../docs/Maths/MultiOutput.md).

**Reference and law gates.** The enumeration normalizes the structural prior on the admitted space F and integrates on ln τ with its Jacobian. The Huang–Wand fixture holds τ = 0.1 and instead integrates both log rates. Its two-dimensional quadrature halves its panel width from 2 to 1 on [−30, 8]; ln Z changes by 2.4·10⁻¹⁴. Neither calculation uses the production likelihood's collapse.

| Gate | Result |
| --- | --- |
| Integrated and linear rung laws under a declared IW, with τ integrated | Count law and mean ln τ agree at four rungs, including the normalized prior rung. Exact ln Z = −6.25269935; ladder estimates −6.23311960 ± 0.00687 and −6.24719222 ± 0.00785 (Monte Carlo errors) |
| SMC under both forms, declared IW and Huang–Wand | All four evidence and final count comparisons pass. IW estimates −6.27682162 ± 0.0157 and −6.26675872 ± 0.00946; Huang–Wand estimates −6.14495274 ± 0.0288 and −6.13978315 ± 0.0137 |
| Huang–Wand count and rates, independent dense quadrature | ln Z = −6.1574765353; P(M = 1) = 0.32862039; E r₀ = 0.36002664, E r₁ = 0.043636699. The plain chain and three ladder rungs meet their integrated references; ladder ln Z = −6.157579796 ± 0.00996 |
| Flat histogram, reweighted count and evidence | IW P(M = 1) = 0.11543239 versus 0.11392457, ln Z = −6.260378266 ± 0.0121; Huang–Wand P(M = 1) = 0.32966483 versus 0.32862039, ln Z = −6.147181961 ± 0.0102 |
| Tempered properness | Accepted and refused sides of ν₀ + Nb − f > K − 1 agree with the dense reference, including the flat boundary exactly; the linear form and requested reference rung keep their own meanings |
| Coefficient forms, dense matrix-t and K = 1 reduction | 108 dense cases per form, K = 1, 2, 3 and all three intercept conventions, τ from 10⁻⁸ to 10⁴ and b from 10⁻⁶ to 1. Maximum relative errors 1.08·10⁻¹⁴ ridge, 1.07·10⁻¹⁴ raw ridge, 9.28·10⁻¹⁵ modified g; K = 1 agrees with scalar BASS |
| Proposed generalized forms and accepted projection reuse | 1,345 admitted proposal rows meet fresh dense fits, maximum relative error 3.03·10⁻¹³. All 1,287 changed pending acceptances reuse the scored pencil: 441 ridge, 441 raw, 405 modified |
| Conditional coefficient draws | 20,000 draws per form, standardized by the independently built dense precision and each sampled Σ. Worst moment absolute z = 2.77 |
| Huang–Wand's prior | JSON and explicit-scale independence, the scalar conditional hierarchy, and 100,000 prior draws checking the half-t scales and uniform-correlation moments; worst absolute z = 0.64 |
| Resume and parallel execution | Twelve resume cases cover all four samplers, both covariance priors, and generalized coefficient forms; structures, τ, Σ, coefficients, rates, weights and evidence resume exactly. SMC's overdispersion and collapse flags agree. Raw-ridge Huang–Wand ladder and SMC draws agree with one and two workers |

The flagged qualification bodies passed 46/46, including 36 pins, four declared-IW and two Huang–Wand Geweke seeds, the two histograms, the integrated-rate ladder and the g-form cost row. Each Geweke case compares 400,000 independent draws with 12,000 independent paths of 50 consecutive sweeps, using path means for its standard error. The largest absolute z was 1.86 for the declared IW and 2.33 for Huang–Wand, whose two log rates are also checked. This qualification snapshot predates the generalized proposal's prior-column handoff fix and accepted-pencil cache retention. These g-form bodies have no coefficient-prior pencil; neither fix changes their numerical path. The source difference is retained with the receipts. The final fast run passed 49/49 with the flag enabled, including the pins, generalized cache assertions, parallel draws and cycle-clock costs. The original uninterrupted seed-20261123 tail finding remains open as reported in the slice 1 recheck below.

[Sampler and generalized-form results](../../../artifacts/e34-slice2/results/generalized-entrypoints.trx) · [Flagged qualification](../../../artifacts/e34-slice2/results/qualification.trx) · [Final cache, draw and cost checks](../../../artifacts/e34-slice2/results/cache-and-cost.trx) · [Qualification provenance](../../../artifacts/e34-slice2/results/provenance.json)

**Clean-commit verification.** A clean archive of `41727ac` builds in Release. The Commit selection passes 2,225/2,225, zero skips, with `BASS_QUALIFICATION` unset (8 min 57 s test duration). The whole regression assembly, Slice scope with the slow tier, passes 2,831/2,831, zero skips, on those same binaries (32 min 11 s test duration, flag unset). All runs use below-normal priority, `DOTNET_PROCESSOR_COUNT=2`, MSBuild node reuse off and at most two xUnit threads. The receipts, settings, source hashes and qualification difference are retained under `artifacts/e34-slice2/results/`; the temporary source snapshots and build trees are removed after verification. [Commit receipt](../../../artifacts/e34-slice2/results/commit.trx) · [Commit settings](../../../artifacts/e34-slice2/results/commit.runsettings) · [Full assembly](../../../artifacts/e34-slice2/results/slice.trx) · [Slice settings](../../../artifacts/e34-slice2/results/slice.runsettings)

**Cost and implementation choice.** The accepted baseline uses one fresh full-spectrum generalized decomposition per proposed structure, O(M³), and O(KM²) for the K cross projections, with the same production `DenseEigen` selector. Repeated evaluations use O(MK²) shrinkage followed by O(K³) to factor the posterior output scale; the slope also solves K systems. An accepted proposal keeps its already-scored immutable `OutputPencil`. Scalar ridge still uses its bordered spectrum. A bordered multi-output generalized update is a possible later optimization, requiring the same dense-value, acceptance-cache and resume gates.

Cycle-clock medians below are indicative µs per call on N = 256 and four predictors, over five interleaved batches of 200 calls. They separate the decomposition and K projections from a full design load and proposal. No other test host was reported. The full marginal includes the output-scale determinant and normalization; the shrinkage row does not. IQRs and the clock's load record are in the TRX.

| M, K | Form | Load and setup | Change proposal | Decomposition | K projections | Cached shrinkage | Full marginal | Full slope |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 8, 2 | g | 43.3 | 5.31 | — | — | 0.320 | 1.16 | 1.48 |
| 8, 2 | normalized ridge | 105.4 | 56.5 | 48.4 | 0.846 | 1.60 | 2.39 | 2.93 |
| 24, 5 | g | 170.6 | 14.4 | — | — | 0.496 | 1.35 | 4.85 |
| 24, 5 | normalized ridge | 460.7 | 331.9 | 313.5 | 5.96 | 4.05 | 4.90 | 5.10 |

Reading accepted ridge forms costs 0.260 and 0.267 µs at the two sizes, with object identity gated separately. The whole-sweep row, N = 1000, p = 4, M_max = 20 on separate chains, reads scalar 34.2 µs; K = 1, 2, 5, 10 read 41.4, 49.7, 86.4, 127.0 µs (1.21, 1.45, 2.53, 3.71×). These measurements price this implementation on these fixtures; they do not establish an asymptotic crossover.

**Remaining scope.** Covariance-default coverage and joint readouts stay in slice 3. The optional Huang–Wand prior does not settle that comparison. Row noise, moving working responses and random effects stay in slice 4; the PCA option stays in slice 5. All shared runtime files were released after `41727ac`; E21's pending named-grouping work was excluded from these snapshots and this commit.

### Slice 1 Geweke recheck (2026-11-09)

The original uninterrupted run reproduced on `b360316`. With 400,000 independent draws, 600,000 successive sweeps, 50 batches and the original bound of 4, seed 20261123 gave mean ln τ of −3.75481 on the independent side and −3.46080 on the successive side. Its 1% quantiles were −11.254 and −8.475. The original probe agreed with the dense matrix-t marginal on all 600,000 visited states (maximum relative error 2.84·10⁻¹³). The qualification change is test-side in `ac9726b`.

The shorter and longer diagnostics drew their successive initial states from the same RNG **after** consuming 100,000 or 400,000 independent draws. Their successive trajectories therefore start at different states and stream positions. The change from z = −3.38 to −7.25 cannot be read as the effect of extending one realization by a factor of four.

**Fixed responses and structures.** The declared IW(6, Ψ), zero centre and Gamma(½, rate 6) scale prior are independent of the response. The independent simulator samples admissible structures from the prior and draws Σ and the centred matrix-normal coefficients. All eight original statistics depend only on these parameters, so that simulator can integrate out Y. The successive simulator generates Gaussian observations, then takes one structure/τ sweep and one posterior draw of Σ and the coefficients for each new response. The original probe used the same random draws and added the reference calculation after recording the state.

The two simulators share the basis evaluator and inverse-Wishart sampler. At a prior draw, each output starts with an independent vector with covariance G⁻¹ across its basis coefficients. The Cholesky factor of Σ mixes the vectors across outputs, and division by √τ gives coefficient covariance Σ ⊗ (τG)⁻¹. The separately drawn intercept has covariance Σ/(12τ). The dense matrix-normal and inverse-Wishart references have their own checks in slices 0 and 1.

At each of 101 equally spaced ln τ values from −25 through 0, the production density on **ln τ** was compared, up to one fixed-structure constant, with the dense matrix-t log marginal plus `u/2 − 6·exp(u)`. All 404 evaluations were finite and on support. The nonempty structure has one hinge and one two-way interaction. The shifted fixture adds (400, −200) to every response, putting appreciable conditional mass in the disputed tail. Each row used 500 warm-up updates and 60,000 updates of the production `BassTauBlock`, with the structure and response held fixed. The reference used Gauss–Legendre quadrature, split at the probability cutoff. Halving the panel width and widening [−50, 6] to [−60, 8] changed the mean and probability by less than 10⁻⁹.

| Terms | Shift | Max log-density error (nats) | Mean ln τ, quadrature / block | Lower-tail probability, quadrature / block |
| ---: | ---: | ---: | ---: | ---: |
| 0 | 0 | 2.84·10⁻¹⁴ | −1.96785 / −1.96644 | 0.035845 / 0.035683 below −4 |
| 2 | 0 | 4.26·10⁻¹⁴ | −1.22178 / −1.21984 | 0.000293 / 0.000217 below −4 |
| 0 | 400 | 1.49·10⁻¹⁰ | −14.40515 / −14.40224 | 0.624001 / 0.624733 below −14 |
| 2 | 400 | 2.94·10⁻¹⁰ | −13.05732 / −13.05695 | 0.097510 / 0.096983 below −14 |

Across the eight readings the largest |z| was 1.46. A separate proposal check on the 12-case interacting design compared 24,452 proposed fits, before any acceptance decision, with dense refits at ln τ = −25, −15, −7 and 0. The maximum absolute marginal difference was 3.53·10⁻¹¹ nats, the smallest admitted angle statistic was 0.000719, and 6,509 refusals agreed with explicit design loads. This checks candidate scoring and support on the sampled proposals. [Fixed structures and replicated paths](../../../artifacts/e34-geweke-recheck/investigation.trx) · [Proposal checks](../../../artifacts/e34-geweke-recheck/proposals.trx)

**Why the successive simulator can miss the tail.** Under the zero centre the intercept conditional mean is ȳ/(1+τ). At fixed τ, observing new Y around the previous fitted level and drawing the intercept again retains a fraction 1/(1+τ) of that level in expectation. At the analytic prior 1% quantile (ln τ ≈ −11.2436), that fraction is about 0.999987. This calculation holds τ fixed. The actual diagnostic samples it. Empirically, 337 of 488 paths starting below the 1% cutoff ended below it 50 sweeps later. Dividing a single path into batches cannot account for excursions that path never observed.

Two controls used the same seeds and the same 400,000/600,000 budgets:

- The committed joint-invariance test uses 12,000 independently initialized paths, 50 consecutive sweeps per path. Only the posterior sweep and composition advance a path after initialization. Its standard error uses the independent path means. The eight original statistics keep the same bound, with a ninth for the analytic lower 1% event (`n·τ` has the law of a squared standard normal). The independent simulator is also checked against the analytic mean ln τ and the tail probability.
- An uninterrupted control kept structure, τ and Σ through all 600,000 sweeps. Before generating every 50th response it redrew only the fitted coefficients from their declared conditional prior given those three quantities. This preserves the joint law at the simulation stage where the previous Y has been discarded. The next production sweep and posterior composition ran as before. It passed all eight original statistics on all four seeds. On seed 20261123 its 1% quantile moved from −8.475 to −10.665. Its below-cutoff frequencies still ranged from 0.00481 to 0.01050 across the four runs.

| Seed | Original single path, ln τ z / worst absolute z | Coefficient-refresh path, ln τ z / worst absolute z | Independent paths, ln τ z / worst absolute z | Independent paths, fraction below prior q₀.₀₁ |
| --- | ---: | ---: | ---: | ---: |
| 20261122 | −0.89 / 1.13 | −0.57 / 1.31 | +0.45 / 1.84 | 0.01103 |
| 20261123 | **−7.25 / 7.25** | −0.18 / 1.30 | −0.52 / 1.86 | 0.00943 |
| 20261124 | −0.53 / 2.38 | +0.10 / 0.96 | −0.44 / 1.49 | 0.00941 |
| 20261125 | −1.36 / 1.92 | −0.83 / 1.13 | −0.56 / 0.85 | 0.01030 |

The controls identify slow exploration of the simulated fitted levels as a contributor to the single-path discrepancy. The original seed's continuous-path tail estimate remains unqualified at this budget. The fixed-structure target, conditional update and replicated joint-invariance gates pass for the fixtures above. The production sampler and priors required no correction on this evidence. [Original four-seed rerun](../../../artifacts/e34-geweke-recheck/original.trx) · [Uninterrupted control](../../../artifacts/e34-geweke-recheck/coefficient-refresh.trx)

The reproduction inputs sit with the local receipts. On a snapshot of `b360316`, the first patch adds the pending seeds and sets the original length to 400,000/600,000. Applying the second patch on top adds the conditional coefficient redraw and tail readout. Run `dotnet test` on the regression test project with `--filter FullyQualifiedName~MultiOutputLawTests.Geweke` and `BASS_QUALIFICATION=1` for either comparison. To reproduce the proposal check, add its source beside the test, mark `MultiOutputLawTests` partial in the snapshot, and select `MultiOutputLawTests.ProposalFits`. [Original overlay](../../../artifacts/e34-geweke-recheck/original-length-seeds.patch) · [Control overlay](../../../artifacts/e34-geweke-recheck/coefficient-refresh.patch) · [Proposal probe](../../../artifacts/e34-geweke-recheck/proposal-check.cs)

**Verification on a clean archive of `ac9726b`.** With `BASS_QUALIFICATION=1`, the four seeds, frozen-structure tests, multi-output references, core and fit tests, and the complete pinned tier passed 94/94. The 33 pins include all ten scalar bit-identity cases and the qualification-only WHAM hash `86f45b016b3f3f47`. The largest independent-simulator |z| against the analytic mean or tail probability was 2.14. The Commit selection passed 45/45 with the flag unset. The runs used Release, below-normal priority, `MSBUILDDISABLENODEREUSE=1`, `DOTNET_PROCESSOR_COUNT=2`, and the dedicated `artifacts/e34-geweke-recheck/` root. [Qualification results](../../../artifacts/e34-geweke-recheck/qualified.trx) · [Commit selection](../../../artifacts/e34-geweke-recheck/commit.trx)

The full regression assembly passed 2,742/2,742 with zero skips and `BASS_QUALIFICATION` unset (Slice scope, two test threads, 18 min 12 s). The initial Commit invocation stopped during restore before producing a TRX. The serial Release build succeeded, then both test-select scopes ran with `-NoBuild` on those binaries. [Full assembly](../../../artifacts/e34-geweke-recheck/slice.trx)

**Remaining scope.** Slice 2 covers the sampler and evidence extensions in the table above. The Σ default decision remains with slice 3's small-n coverage comparison against Huang and Wand. This recheck covers the slice 1 g-form with a declared IW and a sampled τ.

### Slice 1: the core (chip, 2026-11-09)

**Landed:**
- `ad72e2f`: `InverseWishart` (Maths.Distributions) and `OutputMarginal` (gprior), the K × K closed form.
- `19dda31`: the K-output core. `MultiOutputResponse`; the outputs' means on `WorkingResponse`; S_Y and Σ's prior in `BassResolved`; the design's K-output side (`DesignOutputs`, `BorderedCholesky.ProposedWhiten`); the forms on `BassStructureFit`; the target's branch; the draw (`BassOutputDraws`).
- `2b907aa`: fits of several outputs (`Compose`, the run's functionals, the draws' rows across chains, the fingerprint, the readers' refusal) and Σ's prior options (`BassConfig.OutputCovariance`, `OutputCovariancePrior`).
- `878fd19`: the laws, Geweke, recovery and the cost row.
- `b03d926`: the refusal gate names KHAOS's two forms beside the ridge (they move to slice 2, below).
- The docs commit with this report: [Several outputs](../../../docs/Maths/MultiOutput.md), a pointer from BayesianRegression.md, and the notation ledger's entry.

**Gates:**

| Gate | Result |
| --- | --- |
| The inverse Wishart's density against slice 0's reference; its draws against the inverse gamma's at K = 1 on one stream | 8·10⁻¹⁵; 10⁻¹⁴ |
| `OutputMarginal` against the matrix-t reference: 548 cases, K to 4, the three conventions and a frozen level, τ over seven decades, b to 0, units 10^±6 | 3.4·10⁻¹³; 52 refused by both, none by one |
| The design against one scalar design per output and per pair of outputs, by polarization: 232 proposals, 149 acceptances | 7.5·10⁻¹⁵ |
| The target against the matrix-t reference at random structures, τ and b, each convention | 1.3·10⁻¹⁵ |
| K = 1 against scalar BASS | S_Y and the forms to the bit; the target to 3.6·10⁻¹⁶; the draw to 6·10⁻¹⁶; a fit draw for draw, 500 draws, structures and τ equal, σ², coefficients and λ to 2.5·10⁻¹⁴ |
| Scalar BASS untouched | each commit's selection with the pinned tier passed (1,896 of 1,896 at `19dda31`, 1,892 of 1,892 at `2b907aa`); the whole assembly below |
| The law by enumeration: count law, E[ln τ] and six likeliest structures, small F, 200,000 steps | worst \|z\| 3.01 (K = 2), 2.64 (K = 3) |
| Geweke (qualification): zero centre, declared IW(6, Ψ₀), τ sampled, n = 12, K = 2, 100,000 and 150,000 | worst \|z\| 3.38 (ln τ) and 1.78 on two seeds |
| Equivariance Y → YL under decision 2, every convention | −(N − f) ln\|det L\| to 10⁻¹⁵ |
| Recovery: three surfaces outside the model, Σ with correlations 0.3, −0.2, 0, n = 300 | each E[Σ] entry within 1.3 posterior sd; surfaces' RMSE 0.04 to 0.06 against noise sds 0.15 to 0.3 |
| Σ's prior options | each kind resolved, round-tripped through JSON, refused where it does not fit; the target under a declared prior against the reference to 10⁻¹² |
| The refusals | a noise, the ladder, SMC, the normalized ridge, KHAOS's raw ridge and modified g-prior, Sigma2's other kinds, a constant output and collinear outputs, each refused at `BassModel.Create` with its reason |
| The cost row (wall clock, indicative): a sweep on n = 1000, p = 4, M_max = 20, each on its own chain | scalar 30.3 µs; K = 1, 2, 5, 10 at 1.06, 1.26, 1.62, 2.79× |

The whole regression assembly (Slice scope, slow tier included) on a snapshot of `878fd19`: 2,699 of 2,699 in 11 min 21 s.

**Decisions beyond the brief, with their alternatives:**
- *Where K enters the design: a K-output side beside the scalar statistics,* which every K-output path is gated on, so scalar BASS's arithmetic is unchanged by construction. Alternatives: generalize the scalar fields to K columns in place, which reshapes snapshots, pencils and the follow records and puts scalar bits at risk; or a separate design class, which duplicates admission, the gate and the tables.
- *The g-form alone in slice 1.* The main session's call (2026-11-09): the g-form first, the normalized ridge and KHAOS's two forms as a second commit if they come cleanly, else in slice 2. They do not come cleanly, so they move to slice 2. They need the spectrum's bordered candidates and the pencil to carry K cross vectors, and the K × K shrinkage Σⱼ h_j w_j w_jᵀ/λ_j from them; the g-form's τ step is the closed form S(τ) = b(R + τ/(τ + b)F), a pencil whose eigenvalues are equal. Until slice 2 each is refused with its reason, a gate of this slice's.
- *Σ's prior in config:* `BassConfig.OutputCovariance`, null carrying `Sigma2`'s ScaledToResponse(c) and refusing Sigma2's other kinds. Its kinds are decision 2's, Standardized (Ψ₀ = 2c·diag S_Y) and a declared IW. The alternative, carrying every kind from `Sigma2`, cannot state a K × K Ψ₀. The field is left out of the JSON when null, so every existing fingerprint holds.
- *Huang and Wand in slice 2,* with the samplers its K rates thread through (the main session left the slice to the chip).
- *The run's functionals:* ln|Σ| in σ²'s place, which reads the correlations and whose R̂ no rescaling or mixing moves; the monitored points (case, output) pairs. The alternative is Σ's lower triangle, K(K + 1)/2 functionals and another shape for the diagnostics.
- *The draw's stream order* follows the scalar draw's, so K = 1 draws equal scalar BASS's to rounding. `BassDraw` holds `Sigma` and `OutputCoefficients`, its σ² NaN and its coefficients empty, and the readers of one response refuse it with the reason.
- *The hinge prefix tables stand aside for several outputs* (exact borders), a cost the row above includes.
- *The collinearity refusal* reads the correlation matrix's Cholesky pivot 1 − R² against 10⁻¹⁰, under decision 2's form alone; a constant output is refused under every form, as scalar BASS refuses a constant response.
- *Noise and random intercepts are refused at `BassModel.Create`,* with the reason, before the observation model builds a working response that several outputs refuse.

**Unbuilt, and why:** the normalized ridge and KHAOS's priors for several outputs (slice 2, above); the improper |Σ|^(−(K + 1)/2) in config, which `OutputMarginal` supports but whose exact-fit support rule for K outputs is not built; the slope in b (slice 2's evidence); the samplers (slice 2); the readouts (slice 3); noise and random intercepts (slice 4).

**Follow-ups, with the decision each would change:**
- *Σ's lower triangle in the run's diagnostics,* if a run's correlations converge more slowly than ln|Σ| shows. That would change the diagnostics record's shape; nothing has shown the need.
- The four-seed Geweke follow-up and the limits of its original single-path diagnostic are recorded in the slice 1 recheck above.

### Slice 0, addendum: the sources read (chip, 2026-11-08)

**Landed:** `b55a827`, two gates against the sources and the reference's citations. The selected classes with the pinned tier: 59 of 59.

**Brown, Vannucci and Fearn (1998)**, read through a decoding of the PDF's math fonts, since the text layer drops the Greek.
- *Their model*, (4)–(7): Y − 1αᵀ − XB ∼ N(I_n, Σ) on centred X and Y, with B ∼ N(H_γ, Σ) by rows (vec(B) ∼ N(·, Σ ⊗ H_γ), the reference's order) and Σ ∼ IW(δ; Q). γ indexes B's rows, so inclusion is shared across the responses, as the brief cites them.
- *Their shape:* δ follows Dawid (1981), with δ = ν − q + 1 against the conventional degrees of freedom. Their default δ = 3 is ν₀ = q + 2, with Q = kI and E(Σ) = Q/(δ − 2).
- *Their intercept:* α ∼ N(α₀, hΣ) with h → ∞. Its |Σ|^(−1/2) cancels in α's integral, so their (20) reads |Q_γ|^(−(n + δ + q − 1)/2): a count of n, as under level centring. This confirms item 1 of slice 0's findings from the source: n − 1 is neither BVF's count nor scalar BASS's.
- *Gate:* over every subset of five columns, under their H = v₁I and their g-prior (22), the reference differs from their ln g(γ) by one constant, to 2.8·10⁻¹⁴. Their Q_γ (21) is Ψ₀ + R to 2·10⁻¹⁵, and (22)'s shortcuts hold: |H_γ||K_γ| = (c + 1)^(p_γ), and Q_γ = Q + c/(c + 1)·RSS + YᵀY/(c + 1).
- *Read for later slices:* their g-prior with c from 10 to 100 (Smith and Kohn) is the g-form at τ = 1/c. Their sampler is Gibbs on γ componentwise, and Metropolis with add, delete and swap moves in Brown et al. (1998). Their §8's matrix-T predictive is slice 3's.

**Huang and Wand (2013).**
- *Their prior*, (2): Σ | a ∼ IW(ν + K − 1, 2ν diag(1/a)) with a_k ∼ IG(½, 1/A_k²), and their IW(κ, B) conventional.
- *Its properties:* σ_k ∼ Half-t(ν, A_k) (Property 2), and each correlation has density ∝ (1 − ρ²)^(ν/2 − 1), uniform at ν = 2 (Properties 3 and 4). Σ and the a_k are conditionally conjugate: a_k | rest ∼ IG((ν + K)/2, ν(Σ⁻¹)_kk + 1/A_k²).
- *Gate:* inverse Wishart draws with a diagonal scale give that correlation law at K = 2 and 3, for ν 0.5, 2 and 6, within 1.7 se on E[ρ²] and E[ρ⁴]. It holds at each a, since a correlation is unchanged by a diagonal rescaling.
- *At K = 1* their prior is the engine's `Sigma2Prior.Hierarchical(ν/2, ½)` with A = s: σ² | r ∼ IG(ν/2, r), with r = ν/a ∼ Gamma(½, rate 1/(νA²)). The engine's half-Cauchy is their ν = 1. Built for E34, their prior would be the hierarchical prior carried to K outputs: the state carries K rates, each moved as r is, and given them the collapse stays matrix-t. That is what the entry's "its scales would take a slice step as τ does" comes to.

**Found, for the main session (decision 2):** decision 2's ν₀ = K − 1 + 2c is Huang and Wand's inverse Wishart at their ν = 2c. Under a diagonal S_Y, each prior correlation then has density ∝ (1 − ρ²)^(c − 1); at c = 10⁻³ nearly all its mass sits at ±1. This is Gelman's (2006) objection to IG(ε, ε), which Huang and Wand take up. The data overrule it once n is large against K, as they overrule the scalar IG(c, c·s²).
- *The decision it bears on:* decision 2's default against Huang and Wand's prior at ν = 2, which the brief leaves out. For few cases per output, or for correlations read as a result, their prior gives uniform correlations at the cost of K more coordinates in the state.
- *A measurement that would decide it:* the posterior correlation's coverage on a small-n fixture (n ≈ 2K to 4K) under both priors. If the default's coverage holds, nothing changes; if it fails, Huang and Wand's prior becomes the default candidate.

**Notation:** BVF's IW(δ; Q) and the conventional IW(ν, Ψ) of the reference and of Huang and Wand differ by δ = ν − K + 1. That crossing goes in `issues/notation-ledger.md` with slice 1's docs, whose ledger another session is editing now.

### Slice 0: the references (chip, 2026-11-08)

**Landed:** `d5b5dc6`, test-side as E8's slice 0 put its references; no engine code.
- `tests/regression/MatrixTReferences.cs`. `MatrixTReference` is the collapsed marginal by dense algebra. A integrated against MN(A₀, Λ₀⁻¹, Σ) leaves ν = nb − f (f = 1 under the flat intercept), the K × K form R = b(Y − XĀ)ᵀ(Y − XĀ) + (Ā − A₀)ᵀΛ₀(Ā − A₀), a sum of semidefinite terms, and (K/2)(ln|Λ₀| − ln|P|). Σ against IW(ν₀, Ψ₀) then gives ln p_b(Y) through Γ_K, and Σ | Y ~ IW(ν₀ + ν, Ψ₀ + R). Proper and improper priors; exactly 0 at b = 0 under a proper prior; the matrix-t density in the data's n × n row scale at b = 1.
- The routes it is checked by: `DenseVecReference` (vec(Y) as one Gaussian, the Kronecker products entry by entry, the flat intercept by Helmert rotation), `ReferenceInverseWishart` (density, Bartlett draws), `SigmaQuadrature` (K = 2), `ReferenceSigmaPrior.ScaledToResponse` (decision 2).
- `MatrixTReferenceTests.cs`: 21 tests, 6 s.

**Gates:**

| Gate | Result |
| --- | --- |
| A's integral against the dense vec form, K = p = 2, Σ with unequal variances and a correlation, b 1 and 0.4, intercept inside and flat | 3·10⁻¹⁶ relative. Λ₀⁻¹ ⊗ Σ in Σ ⊗ Λ₀⁻¹'s slot misses by 0.27 to 4.0 nats |
| The same over a seeded sweep: K ≤ 4, p ≤ 4, n ≤ 9, b to 10⁻⁶, units 10^±6 | 250 moderate cases 1.9·10⁻¹². 250 far cases (cond Σ to 10⁸, Λ₀ over eight decades) within 1.3·10⁻¹⁴·κ̂, κ̂ to 7·10¹² being the dense Cholesky's conditioning |
| The matrix-t density, untempered, 200 cases | 1.7·10⁻¹² |
| Σ's integral against quadrature at K = 2: informative prior, b = ½, decision 2 under level centring, flat intercept | 1.4·10⁻¹³ nats; a coarser grid within 8·10⁻¹⁴; faces 44 nats below the peak. E[Σ] against (Ψ₀ + R)/(ν₀ + ν − 3) to about 10⁻¹³ |
| The inverse Wishart's normalization by the same rule, ν 4.5 and 11 | ln ∫ = −2.9·10⁻¹³ and −1.1·10⁻¹³ |
| Monte Carlo over Bartlett draws from the prior, K = 3, 3·10⁵ draws | z 0.63 (b = 1, se 0.008) and −1.44 (b = ½, se 0.004) |
| K = 1 against the scalar block under decision 2's IW(2c, 2c·s²) = IG(c, c·s²): g-form, normalized ridge, pencil; level centring, zero centre, flat; Jeffreys; τ 10⁻⁶ to 10⁴; units 10^±6; b to 0 | 574 compared, 4.1·10⁻¹⁴; 26 refused by both, none by one alone |
| Equivariance under decision 2: Y → YL with scales 10^±3, mixing and a sign | every structure shifts by −(nb − f) ln\|det L\|, M = 1 and 3, to 7·10⁻¹⁴ relative |

Whole regression assembly on `15eb151` and the slice: 2,573 of 2,573.

**Found, for the main session (decisions 1 to 3):**
1. *Level centring's count is nb, not b(n − 1).* Scalar BASS's level centring keeps the intercept inside the g-prior, centred on ȳ, k = M + 1 ([BayesianRegression.md](../../../docs/Maths/BayesianRegression.md), "Why level centring"). Integrating it spends no observation, so Σ | Y ~ IW(ν₀ + nb, Ψ₀ + R), and slice 1's gate (K = 1 equal to scalar BASS) holds only with nb.
   - n − 1 is the flat intercept's count at b = 1, and nb − 1 is its tempered count.
   - b(n − 1) would be a third path, the restricted likelihood (Y projected off 1) tempered. It is neither BASS's integrated form nor equal to scalar BASS at K = 1.
   - *Proposed wording:* decision 1, "level-centred per output, inside the prior as scalar BASS's is, so ν = nb"; decision 3, "ν₀ + nb − f > K − 1".
2. *Under the flat intercept, decision 3's "only under an improper ν₀" does not hold.* There ν₀ + nb − 1 > K − 1 fails at nb ≤ K − ν₀ under a proper prior too: at nb ≤ 1 − 2c under decision 2's ν₀. b = 0 diverges, as scalar BASS's flat pole does. Level centring and the zero centre keep decision 3 as written. The reference refuses where the condition fails.
3. *Decision 2's S_Y.* For K = 1 to be the scalar default exactly, S_Y takes divisor n and is read about the centre levels, as the scalar s² is (`BassResolved`). Those are ȳ under level centring and the flat intercept, and 0 under the zero centre, where it is the second moment.
4. *Where the entry's first condition came from.* Integrating Σ against the tilted prior IW(ν₀ + n(b − 1)), the matrix-t's own form below b = 1, needs that prior proper. Integrating the tempered likelihood directly needs only ν₀ + nb − f > K − 1. This confirms decision 3's correction.

**Decisions beyond the brief:**
- *References test-side, no engine code.* The alternative was an engine inverse Wishart and the K × K marginal now, gated against the dense forms. Test-side keeps slice 1's engine block gated against routes that share none of its code.
- *Vec order:* vec stacks A's columns, one per output, so Cov(vec A) = Σ ⊗ Λ₀⁻¹. Slice 1's draws of A | Σ follow it.
- *The quadrature in Bartlett coordinates of Σ⁻¹, with a second pass on the first pass's moments.* A box over the Cholesky factor of Σ missed the tails: its off-diagonal coordinate spreads as e^(ln L₁₁). The Jacobian is gated by the normalization.
- *Monte Carlo through the collapse at each draw.* The dense form checks the first 100 draws to 10⁻¹¹; at 3·10⁵ draws it cost more than the test's budget. Nothing of Γ₃ enters the draws' side.

**Unbuilt:** the reading of Brown, Vannucci and Fearn (1998) and Huang and Wand (2013), not on disk; asked of the main session on 2026-11-08 (read since; the addendum above). Nothing here depends on them. The matrix-t is cited from Dickey (1967) and Gupta and Nagar (2000), from memory, and derived in the reference's comments.

**Follow-ups:**
- *Items 1 to 3:* slice 1 builds level centring with ν = nb, S_Y about the centre levels with divisor n, and decision 3's condition with f, unless the main session or the user decides otherwise. Another choice changes the count, and with it the K = 1 gate.
- *Brown, Vannucci and Fearn, once supplied:* their marginal against this reference. If their intercept convention differs, the difference goes in docs/Maths with this one; otherwise nothing changes.
