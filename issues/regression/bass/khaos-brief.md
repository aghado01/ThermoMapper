# KHAOS on the BASS core: E19's build (chip brief)

Date: 2026-11-07. Status: slice 0 landed 2026-11-07 (`31362ad`); slice 1 landed 2026-11-07 (`fa05048`, `b408bf4`, `baf69f9`), with declared input domains (`3eb56e1`, the user's decision of 2026-11-07); the check on slice 1's input-box finding done 2026-11-08 (real, not a trap); slice 2 landed 2026-11-08 (`f289af3`); the rest in the user's order of 2026-11-08 ("Yes, 5 first"): slice 5, then 4, then 3, each after the coordinating session's verification of the one before; slice 5 landed 2026-11-08 (`98b0e7e`); slice 4 landed 2026-11-08 (`8f57487`); slice 3 landed 2026-11-08 (`0a2168c`), random-effect intercepts as a capability of the core, which completes KHAOS's slices, verified by the coordinating session 2026-11-08, and E8 the build order's next item, waiting for the user to commission it; and SMC's trap is flagged (`42739ac`, the coordinating session's call of 2026-11-08), its remedies recorded, not built. Written at the user's request after GBASS's slices 0–4, from the KHAOS paper's TeX, the `khaos` 2.1.1 source and the engine.

## Why

KHAOS (Rumsey, Francom et al., *Bayesian Adaptive Polynomial Chaos Expansions*, arXiv 2510.25036, in `gauntlet/bayesian/2510.25036v1`) is an adaptive Bayesian polynomial chaos expansion:
- **the basis:** reversible jump over sparse products of standardized shifted-Legendre polynomials;
- **the birth:** a Bernoulli-inclusion proposal with rejection of empty and over-order draws;
- **the coefficient prior:** a modified, n-component g-prior that shrinks each term by its complexity, or a fixed ridge.

It serves smooth responses, computer experiments above all, where polynomial chaos is the standard surrogate and Sobol indices read straight off the coefficients.

It is the same lineage as BASS on other axes: the basis (E19), the proposals (E4, which builds F-X4) and the coefficient prior (E5 through E21). It moves no latent, and it is deterministic and conjugate. So, as GBASS was, it is a diff against the engine, not a new engine. The build order is the user's: GBMARS, then KHAOS, then E8 ([the build order](lineage-review.md#build-order)). Once both exist, their product costs nothing more: a polynomial-chaos basis under GBASS's noise.

## The oracle: CRAN's `khaos` 2.1.1

Pinned in `r/` (`fb38e44`), with its imports glmnet, shape and survival; mvtnorm was already locked.
- **What it holds:** `adaptive_khaos` with `prior = "g"` (`adaptive_khaos_gprior`) or the ridge (`adaptive_khaos_ridge`); `sparse_khaos`, Shao et al.'s forward-selection SBPCE with five enrichment rules; `ordinal_khaos`; `sobol_khaos`; and their `predict` methods.
- **Reading its source.** It installs as a binary, so its source is read by deparsing the namespace. The main session did that on 2026-11-07 into its scratchpad; regenerate it with `deparse` over `ls(asNamespace("khaos"))`.

## What is already built

- **Term types (F-X3).** `IDesignKernel<TPoints>` per part, with `DesignKernels` holding the continuous hinge product, the categorical subsets and the functional hinge product. Each part supplies `Multiply`, `MultiplyDerivative` and its separable factors, and `BassMoments` serves every type's moments by one rule, a term being κ·Πᵥ φᵥ(xᵥ). `PiecewisePolynomial` integrates one-piece polynomials, so Legendre factors need no new moment code (E19).
- **The structure prior and proposals (D01).** `BassStructurePrior` holds the Poisson count with λ integrated, times a term measure. `BassTermProposal` holds D01's birth (Wallenius), death and change, and the fixed-M family.
- **The move contract.** `IRjMove<TState>` splits Green's factors between engine and move. F-X4, the proposal-kernel contract with exact forward and reverse probabilities, is not built: it was decided to be built with KHAOS's Bernoulli-inclusion proposal.
- **The collapsed coefficient blocks.**
  - `CollapsedGPrior`, the g-form.
  - `CollapsedRidgePrior` with `RidgeSpectrum`, the normalized ridge.
  - `PriorPencil`: a posterior precision bA + τD with a prior precision D apart from the likelihood's Gram A. It is decomposed once per structure, after which any (τ, b) costs O(p). KHAOS's prior precision is D_g⁻¹ΨᵀΨD_g⁻¹/g₀², so this is the pencil with τ = 1/g₀².
- **The Lean groundwork.**
  - SC-3, proved (`khaos_precision_entry`): G∘ΨᵀΨ = ΨᵀΨ + S₀⁻¹, so the modified g-prior stays conjugate.
  - SC-1: KHAOS's printed ratios ([the register](lineage-gaps-register.md#manuscript-readings-to-confirm)).
  - SG-4: `shiftedLegendre_jacobi`, the Jacobi generator's eigenvalues, an option beside the complexity penalty.
- **The exact law on small spaces.** The exact-gap work's `ReversibleChain` and `EliminatedChain` (`tests/regression/BassKernelSpectrum.cs`) enumerate a structure chain's kernel. A Legendre term space at small p, d and q is enumerable, so the law is checkable exactly.
- **What else is reusable.** `SliceSampler`, with its struct overload from GBASS slice 4. `GaussLegendre`. GBASS's noise laws, for the product.

## The diff

### New
- **The Legendre part (E19).** A case-side `IDesignKernel` whose part is a multi-index over the term's variables: Πⱼ ψ_{αⱼ}(xⱼ), with ψ_α(x) = √(2α + 1)·P_α(2x − 1).
  - Values come from the three-term recurrence, and derivatives from the recurrence for P′.
  - The factors are one-piece polynomials for `PiecewisePolynomial`.
  - Under the uniform measure the ψ are orthonormal, so Sobol indices are sums of squared coefficients: an independent check on `BassMoments`.
- **The degree-split term space.** A term is its variables with their degrees.
  - Order q = |{j : αⱼ > 0}| and degree d = Σαⱼ.
  - The measure is uniform over 𝒜_{p, d_max, q_max}, whose size is Σ_{i≤q} Σ_{j≤d} C(p, i)·C(j − 1, i − 1).
  - The count prior is `BassStructurePrior`'s, with λ integrated.
- **F-X4, the proposal-kernel contract,** built with E4's Bernoulli-inclusion birth:
  - q₀ is drawn from {1, …, q_max} ∝ q₀^(−s_q);
  - each variable's inclusion probability ηⱼ comes from the current term counts, summing to q₀ (`make_weights`);
  - inclusion is χⱼ ~ Bern(ηⱼ), redrawn while q = 0 or q > q_max;
  - the degree is drawn from {q, …, d_max} ∝ d^(−s_d), and split by a uniform composition into q positive parts.
  - The forward probability is exact: Σ_{q₀} J(q₀)·P(χ | q₀)/P(1 ≤ q ≤ q_max | q₀), each validity taken from the Poisson-binomial. The death's reverse is the same expression.
  - The two mutations, degree and variable swap, carry exact ratios.
  - D01's proposals become one member of the contract, their draws unchanged.
- **The Poisson-binomial** in `Maths.Distributions`: pmf and cdf by dynamic programming over the p coins, in O(p·q_max).
- **The modified g-prior (E5) as a prior-precision builder (E21).**
  - The prior precision is S₀⁻¹ = D_g⁻¹ΨᵀΨD_g⁻¹/g₀², with g_m = (1 + q_m(d_m + q_m − 2))^(−ζ/2), through `PriorPencil` with τ = 1/g₀².
  - g₀² ~ IG(a_g, b_g) is moved by one `SliceSampler` step on ln g₀², under its conditional derived with β and σ² integrated: the IG prior, |Σ_n|^½|S₀|^(−½) and the data term.
  - ζ = 0 is the g-prior with g₀² drawn.
- **Ridge blocks (E21).** A prior-precision block for random-effect intercepts on a grouping variable, with its centring or sum constraint declared beside it, so the block does not overlap the level-centred intercept. Every grouped-data problem then uses it.
- **`docs/Maths/PolynomialChaos.md`,** a page of its own, beside `VarianceMeanMixtures.md` and `ConformalPrediction.md`.

### Upgrades
- **`BassTerm` and the design kernels.** The contract was written so later members join without changing it, but `BassTerm` is closed (port-fitting F-X3). A Legendre part may need it opened, or a part kind added. The build decides how, and reports it.
- **`BassStructurePrior`** gets a term measure per part kind. A hinge fit's measure is unchanged.
- **The collapsed core** takes a prior-precision builder, not one more enum case (E21). The g-form and the normalized ridge become two builders whose draws are unchanged.
- **`PriorPencil`** serves g₀²'s steps: one decomposition per structure, then O(p) per g₀².
- **The Legendre columns, for speed.** They are evaluated in lanes, with new `Vector256` loops marked `AggressiveOptimization`. Each variable's ψ_α columns up to d_max may be cached, chosen by size, never by hit history, since the cached values equal the recomputed ones bit for bit.

### Left alone
- **Hinge BASS and GBASS's noise laws.** Their draws are unchanged.
- **The product with GBASS's noise.** It needs nothing new and is exercised by one test.

## Decisions

**Recorded, and executed here.**
- *The order:* GBMARS, then KHAOS, then E8 (the user, 2026-11-04).
- *F-X4:* built with KHAOS's Bernoulli-inclusion proposal (port-fitting checklist).
- *g₀²:* one `SliceSampler` step on ln g₀² under its derived conditional, in place of the paper's Laplace-approximation independence proposal (E19).
- *Adaptive mutation probabilities:* under D10's policy.
- *Corrections, not parity variants:* each defect below that is confirmed is corrected, and the tests assert the corrected reference, as GBASS's were.

**Recommended, and accepted by the user on 2026-11-07 when the chip for slice 0 started.**
1. **Legendre as a part type behind F-X3, not a separate engine.** The tempering ladders, SMC, the flat histogram, the pointwise density, the conformal layer and GBASS's noise then serve it as they serve hinges. A basis that mixes term types stays an option, not built.
2. **β and σ² integrated in the structure step,** as BASS does and as KHAOS's own ratio already is. β is drawn for readouts only.
3. **Uniform compositions:** q − 1 distinct cuts from {1, …, d − 1}, with probability 1/C(d − 1, q − 1). This corrects `random_partition` and the C(d, q) in the ratios (below).
4. **Exact Poisson-binomial validity on both sides,** in place of the birth's product over attempts and the death's normal approximation.
5. **The mutation-type probability fixed at ½ by default,** degree only when p ≤ 3. KHAOS's rule adapted to acceptance rates (floor 0.1) is available under D10, as diminishing adaptation or burn-in only.
6. **Inputs on [0, 1]** by the fit's declared ranges, recorded on the fit, as BASS's are. Moments are exact under the uniform measure. `MixtureMeasure` and `EmpiricalMeasure` come through the moment service: E19's "past the paper".

## The defects to run (read 2026-11-07; each to confirm or strike)

| Defect | Where | What slice 0 runs |
| --- | --- | --- |
| The g-prior's structure ratio leaves out the prior's normalizer. The birth, death and degree-change ratios read ½Δln\|Σ_n\| + (a_σ + n/2)Δln Q, without ½Δln\|S₀⁻¹\| = −Δ[(M + 1)ln g₀ + Σₘ ln g_m] + ½Δln\|ΨᵀΨ\|. The ridge path does carry its ½ln(1/τ²). | `adaptive_khaos_gprior`; the same omission in `log_evidence_fast_core` and `log_evidence_full_core`, which `sparse_khaos`'s Bayes factors read | The closed-form marginal ln p(y \| structure) at frozen structures against the code's terms. If confirmed: on an orthogonal design the code's penalty per term is about ½ ln n whatever g₀² and g_m are, where the stated model's is ½ ln(1 + g₀²g_m²), so the complexity factor barely reaches the structure moves |
| g₀²'s conditional has no data term (E5). `dgsq_full` is θ^(−a−M/2)·e^(−b/θ)·\|Σ_n\|^½. With β and σ² integrated, Q moves with g₀²; conditioned on them, βᵀS₀⁻¹β/σ² enters. | `dgsq_full`, `rg0sq_laplace_*` | The joint less each candidate conditional's log kernel at frozen states, as GBASS's slice 0 did |
| The birth's forward probability multiplies every attempt's: `delayed_reject_term` adds ln P(χ) for each redraw, invalid ones included. The death's reverse divides the q₀-marginal by the q₀-marginal validity. Since q₀ is drawn once before the redraws, the exact forward is Σ_{q₀} J(q₀)·P(χ \| q₀)/P(valid \| q₀), and both sides differ from it whenever a birth redraws. | `adaptive_khaos_*`'s birth and death | The forward probability by enumeration at small p, against the three forms |
| `random_partition` is not uniform over compositions: it sorts q − 1 cuts drawn with replacement, which weights compositions unequally. The ratios use C(d, q), which is neither the uniform 1/C(d − 1, q − 1) nor the draw's probability. | `random_partition`; the birth, death and Mutate1 ratios | The draw's distribution by enumeration at small d and q, against both counts |
| The variable swap can pick a variable already active in the term: only the swapped-out variable's weight is zeroed. That gives a term with a repeated variable, which is outside the PCE basis. | Mutate (variable) | The count on a seeded chain, and a constructed state |
| The death's reverse takes P(q > q_max) by a normal approximation. | `adaptive_khaos_*`'s death | Exact Poisson-binomial against it, at small p |
| σ²'s floor `s2_lower`, as GBASS's floor on w was. | All variants | Not ported; recorded |
| The paper's typos (E19's list): P₁(x) = 1; the order sum over M; the degree from {q, …, p}; the partition probability binom(d − 1, d − 1)⁻¹. | The TeX against the code | Each confirmed in the code's own reading; the code draws the degree from {q, …, d_max} |
| SC-1's printed ratios against the code's | The paper's birth, death and Mutate1 | Each factor's direction, by enumeration |

## The slices

One commit per concern; tests with the code they test; each slice ends with its gates and a fresh-eyes review before the next starts.

| Slice | Delivers | First needs | Its own gates |
| --- | --- | --- | --- |
| 0. Oracle modes | `khaos_*_oracle.R` and their parity tests; the defects above run | None; the package is pinned | The modes run from `r/`; each defect confirmed or struck |
| 1. The Legendre basis under the ridge | The Legendre part, the term space and its measure, F-X4 with the Bernoulli-inclusion birth and both mutations, the Poisson-binomial, under the normalized ridge and KHAOS's raw ridge | Slice 0's corrected references | The exact law by enumeration on small term spaces (`ReversibleChain`); Geweke; `adaptive_khaos_ridge` readouts; the hinge fit's draws unchanged |
| 2. The modified g-prior | The prior-precision builder, the complexity factors, g₀² by slice steps through `PriorPencil` | Slice 1 | g₀²'s conditional against the joint at frozen states; Geweke with g₀² drawn; `adaptive_khaos_gprior` readouts under the corrected marginal; ζ = 0 against `CollapsedGPrior` with g₀² drawn |
| 3. Ridge blocks | E21's random-effect intercepts, their convention declared | Slice 2's builder | Recovery of known group effects on data from outside the model; the identifiability convention under level centring |
| 4. Readouts past the paper | Sobol from coefficients beside `BassMoments`; moments under `MixtureMeasure` and `EmpiricalMeasure`; derivatives for E8 | Slice 1 | Sobol against coefficients under the uniform measure, and against `sobol_khaos`; moments against quadrature; derivatives against central differences |
| 5. The engine's reach | Tempering, SMC, the flat histogram, the pointwise density, the conformal layer's candidate density, and GBASS's noise on Legendre terms | Slices 1–2; GBASS slice 5 | The law at each rung; one Student t and one quantile KHAOS fit; conformal coverage on one Legendre fit |

## Gates common to every slice

- **Hinge BASS untouched.** The eight pinned hashes and every existing fingerprint hold: the pinned tier at each commit, the whole assembly at the slice's end ([targeted runs](../../../CONTRIBUTING.md#testing-standards)).
- **The law.**
  - Exact by enumeration on term spaces small enough to list.
  - Geweke's joint-distribution test on a design small enough to run long, with proper priors.
- **Recovery from truth outside the model,** on smooth test functions with known Sobol indices: Ishigami, Sobol's g-function, the borehole. Coverage of the true surface and of held-out responses, reported beside the width. These are the problem class's standard benchmarks; none is a fixture of any one application.
- **Primitives.** Chosen cases and a seeded random sweep, far ranges included:
  - the recurrence to degree 30 and beyond, with orthonormality by Gauss–Legendre quadrature;
  - the Poisson-binomial against brute-force enumeration;
  - |𝒜| against counting.
- **Refusals and acceptance.** Each move's counters at each point.
- **The cost row.** A Legendre column against a hinge column; a structure sweep at the cost report's base and at p 20, d_max 15, q_max 5; g₀²'s step. All on the cycle clock, marked indicative.
- **`docs/Maths`.** Decided statistical choices land on the new page with their alternatives and implications. New writing takes the module's symbols, and crossings go in `issues/notation-ledger.md`. KHAOS's d, q, α, g_m and g₀ meet BASS's degree, interaction order and g.

## Not in this brief

- **E8.** It follows KHAOS and reads slice 4's derivatives.
- **`sparse_khaos`,** Shao et al.'s forward selection. It is a point estimate, not a posterior. It may serve as an initializer or a benchmark (E27), and is recorded with the decision it would change if one comes.
- **`ordinal_khaos`.** An ordinal response through a latent normal with thresholds belongs to E20's data-augmented observation models, on the carrier GBASS built. With this basis it is a later product.
- **SG-4's Jacobi penalty,** an option beside g_m, posed in Lean.
- **Bases that mix term types.** An option the contract allows, not built.

## Working rules

- On `main`, committing by pathspec with nothing left staged and nothing pushed. Commit messages describe the work, with no attribution lines ([CLAUDE.md](../../../CLAUDE.md)).
- Other chips may be writing in the checkout. Before editing a file, check `git diff -- <file>` for changes you did not make; build from a snapshot (`git archive HEAD` plus your files); commit only compiling states.
- Builds set `ArtifactsRoot` inline into `artifacts\chip-khaos\`, deleted at the end.
- Heavy runs at below-normal priority with capped threads, MSBuild node reuse off; the machine is shared. Tests run by [targeted runs](../../../CONTRIBUTING.md#testing-standards) (`scripts/test-select.ps1`, `-SourceRoot` at your snapshot): the selected classes per edit, with the pinned tier per commit, the whole assembly once at the slice's end, and the qualification gates as each slice states.
- New beside old; a removal is its own question to the user, with what is lost.
- Findings framed by problem class, not by any one application.
- A follow-up goes in a report with the decision it would change and what each outcome would do. A question with no decision waiting on it is recorded as open, not proposed.

## Report

Under "Reports" below, newest first, per slice: what landed by commit; the gates' numbers; the parity comparisons with their stated differences; refusals and acceptance; the cost row; defects confirmed, struck or newly found; decisions taken beyond this brief; what was left unbuilt and why.

## Reports

Newest first.

### Slice 3: random-effect intercepts (chip, 2026-11-08)

At the user's order (5, then 4, then 3), after the coordinating session verified slice 4. With it, KHAOS's slices are complete, and the build order's next item is E8.

**Landed** on `main`, nothing pushed:
- `0a2168c`: the block and its gates.
  - `BassRandomIntercept` and `RandomInterceptBlock` (`BassRandomIntercepts.cs`).
  - The grouping on `BassData`, the config's field, the fingerprint and the refusals.
  - The whitening, through `WorkingResponse`'s log-Jacobian and `BassDesign`.
  - The group effects' draw in the sweep's composition (`BassDraw.GroupEffects`).
  - Prediction, the mean surface, the training fit and the retained pointwise content, each reading f + u_g.
  - The gates: `RandomInterceptTests` and `RandomInterceptQualificationTests`. `BassRunTests.AssertSameFit` now compares the group effects too.
- This commit: `BayesianRegression.md`'s section, its config row and source map; `PolynomialChaos.md`'s row; the notation ledger's entry; the status line and this report.

**What is built,** as a capability of the core, not a KHAOS feature. Any grouping and any coefficient form, on hinge and Legendre terms alike, as [the section](../../../docs/Maths/BayesianRegression.md#bass-random-effect-intercepts-e21) states it.
- *The model:* y = α + Σⱼ βⱼBⱼ(x) + u_g + ε, with u | σ² ~ N(0, σ²s²Π) on the observed levels. s is the effects' sd over σ, and it is held.
- *The grouping* is a `CategoricalPredictor` that `BassData` holds beside the predictors. No term reads it, so the hinge's categorical terms leave it alone, and the Legendre basis's refusal of categorical inputs does not fire.
- *The block is integrated out exactly,* by whitening the response and every column by V^(−½), V = I + s²(ZH)(ZH)ᵀ.
  - The price −½ln|V| rides on the working response as the weights' Σ ln w does, so every marginal, rung and particle, and the evidence, carries it.
  - The effects are drawn last in each composed draw, from their exact conditional given α, β and σ².

**The convention under level centring,** declared as `RandomInterceptConvention.ObservationCentred`.
- Σ_g n_g u_g = 0, so Zu has mean zero over the cases and the intercept alone carries the level.
- V1 = 1, so the whitening keeps the intercept's column, every column's mean and ȳ.
- The unweighted sum-to-zero is not built. On unbalanced levels it leaves Zu a nonzero mean, which overlaps the level-centred intercept.

**The engine's reach,** as the coordinating session asked. The section's table holds the detail.
- *Run and gated,* on a hinge fit (the g-form, τ drawn) and a Legendre fit (the modified g-prior, g₀² drawn):
  - a chain, and the ladder with its reference rung: each rung's law against the enumeration;
  - SMC and its flags: ln p(y) against the enumeration, with the flags clear;
  - the ladder's evidence, against the enumeration;
  - the flat histogram, against the enumeration (Legendre);
  - resume, bit for bit;
  - the retained pointwise content.
- *Refused, each with its reason in the message:*
  - the collapsed pointwise content: leaving a unit out would leave it out of the block too, whose whitening couples its level's cases;
  - the candidate density, and so the conformal layer: a new case's group effect, and the exchangeability across groups that the layer reads, need a reading of their own;
  - GBASS's noise: the block whitens under unit weights, and the latent weights move every sweep;
  - a direct functional response, and a decomposition's components, which share no grouping.

**The gates.**

| Gate | Reading |
| --- | --- |
| The whitening against dense algebra | W·V·W = I to 2.6·10⁻¹⁵, W symmetric, W·1 = 1 to 2·10⁻¹⁶, and ln\|V\| against a dense Cholesky to 10⁻¹², on 30 cases over six unbalanced levels |
| The effects' conditional against dense algebra | The mean s²ΠZᵀV⁻¹r to 8·10⁻¹⁶. Over 200,000 draws, the covariance σ²(s²Π − s⁴ΠZᵀV⁻¹ZΠ) to 0.4% of the largest variance and the mean's worst z 1.6. Σ n_g u_g = 0 to 5·10⁻¹⁵ in every draw |
| The target, absolutely | ln p(y \| S, τ) against the dense multivariate t with scale V + ΨS₀Ψᵀ under σ²'s inverse gamma: to 7·10⁻¹⁴ for the raw ridge and 4·10⁻¹⁴ for the modified g-prior with the block, and to 10⁻¹³ without it |
| Each rung's law | Ladders of 3 and 4 rungs, the reference rung included, against the enumeration of every structure's tempered law, on a hinge fit and a Legendre fit: within 5 standard errors at every M and for E[ln τ] |
| The evidence | SMC (2,000 particles, 4 replicates): hinge −16.550 ± 0.015 and Legendre −74.124 ± 0.018, against the enumeration's −16.539 and −74.106, the flags clear. The ladder: −16.472 ± 0.030 and −74.113 ± 0.038 |
| The flat histogram | Against the enumeration on the Legendre fit with the block, as slice 5's gate reads it |
| The pointwise density | The retained content against N(yᵢ \| fᵢ + u_g, σ²) to 6·10⁻¹⁶ over 60 draws. The group effects are reported as conditioned on, at the group's scope. The collapsed content and the candidate density are refused |
| Prediction and the mean surface | A case of an observed level reads f + u_g to 10⁻¹². A population case's predictive spread over σ² is 2.67 against 1 + s² = 2.69, and an observed level's is 1.006. An unseen level is refused, or placed at the population. The mean surface against the dense joint posterior mean, E[Ψa + Zu \| S, y] with C = V + ΨΨᵀ/τ, to 2·10⁻¹⁴ |
| The convention | Each draw's mean fit over the cases is the intercept's and the terms', to 10⁻¹². The fits' mean level is 2.1341 against ȳ 2.1333 |
| Recovery from outside the model | Twelve unbalanced levels (4 to 37 cases), and effects drawn as 0.8·(E − 1) with E ~ Exp(1), sorted onto the levels by size, which the block's centred Normal does not hold. Four datasets, s = 2. The hinge fit's 90% intervals cover 44 of 48 effects and the Legendre fit's 43, with standardized errors' mean squares of 1.03 and 1.04. σ comes out at 0.29 to 0.33 against 0.3, and at 0.38 to 1.0 without the block. At held-out cases of the observed levels, the RMSE is 0.11 to 0.15, and 0.39 to 1.0 without the block. The population surface is no worse in any dataset, and 3× better where the effects are large |
| The convention, from outside | The fitted level lies within 2 posterior sd of the surface's mean plus Σ n_g u*_g / n in all eight fits. Its squared gaps sum to 0.0018 (hinge) and 0.0017 (Legendre), against 0.35 from the levels' unweighted mean |
| Hinge BASS | The whole regression assembly at the slice's end, once, on HEAD `76f3781` with the slice: 2,396 of 2,396, the pinned hashes and the SMC baseline among them. The commit's selection with the pinned tier: 1,710 of 1,710, after the defect below was fixed (17 had failed). Every fit without the block draws as before: its effects are drawn last, and the pinned tier holds. The slow tier was rewritten from that profile (48 methods), with three of the new gates among them |

**The cost.** Wall clock, one run on the shared machine, indicative.
- *The qualification's fits* (n 246; 2 chains of 2,000 + 1,000 sweeps on Legendre terms, 3,000 + 1,000 on hinges): the four Legendre fits took 9.5 s with the block and 6.9 s without, and the four hinge fits 0.3 s and 0.2 s.
- *What that measures.* Without the block, σ is larger and the structures smaller, so the pair bounds the whitening's cost rather than isolating it.
- *Where the cost sits.* Each column the design evaluates adds O(n + L²) for the whitening. The hinge's prefix tables, which build raw columns, are off under the block. The block's setup is one eigendecomposition of size L − 1.

**Found here.**
- *My own defect, caught by the commit's selection.* The log-Jacobian's line went in between an `if` and its `else if` in `WorkingResponse`'s constructor. The frozen prior Gram's unit weights then replaced a weighted response's own weights, and 17 frozen-route tests failed, `BassMovingResponseTests` and `BassPriorReferenceTests` among them. The line now follows the chain, and all 196 tests in those classes pass.
- *Two readers took the whitened response for the observations.* The retained pointwise content read the working response, Wy, and the training fit would have put the whitening back as an offset. Both were found in review and fixed, and the dense check holds the first.
- *The first qualification gates were unsound,* and they failed before they were rewritten.
  - The effects' RMSE was gated against their realised spread. Hinge dataset 0 drew a spread of 0.32, so the ratio read 0.27 against 0.25, at the posterior's own spread.
  - The population surface was gated as better with the block than without. Where the effects are small, approximation bias dominates both fits, and Legendre dataset 2 read 0.0796 against 0.0783.
  - The rewrite gates standardized errors pooled over the 48 levels, and the predictive mean at held-out cases of the observed levels, where the block's gain is the effects themselves.
  - The surface is now gated as no worse, by at most 10%.
  - The generator now sorts the effects onto the levels by size, so that the convention's target and the unweighted alternative differ in every dataset.

**Decisions beyond the brief.**
- *Whitening, not columns in slice 2's builder.* The brief placed the block on `BassPriorPrecision`.
  - As columns there, the effects would join every structure's factor and residual sum of squares, through every move, rung and particle.
  - Whitening keeps the factor, the samplers and their caches as they are, and its marginal is exact.
  - The cost is a held s: a drawn s re-whitens every column, O(nM) per value.
- *Edits in files this chip does not own.*
  - The conformal session owns `posterior/BassPrediction.cs`, and its tree was clean. My edit adds the group effects to the latent, a population case's own effect to the predictive draw, and E[u | S, τ, y] to the mean surface.
  - `presentation/BassTrainingFit.cs` puts back no offset under the block.
  - `output/BassComponents.cs` refuses a component config that declares the block.

**Left unbuilt, and why.**
- *A drawn s,* a variance component: each value re-whitens, at O(nM). It is E21's next step if a consumer needs one.
- *More than one grouping, crossed or nested, and random slopes.* The whitening generalizes, but V loses its two-level structure.
- *The collapsed pointwise content, the candidate density and conformal prediction,* for the reasons in their refusals.
- *GBASS's noise, weights, and a functional response under the block.* Each needs V^(−½) under moving or non-unit weights, or across locations.

### Slice 4: readouts past the paper (chip, 2026-11-08)

At the user's order (5, then 4, then 3), after the coordinating session verified slice 5 and SMC's flags.

**Landed** on `main`, nothing pushed:
- `8f57487`: the coefficient route (`BassSensitivity.FromCoefficients` → `PolynomialChaosSensitivity`), the gradient readout (`BassFit.Gradient` → `BassGradient`, through `BassSurface.EvaluateDerivative`), and their gates (`KhaosReadoutTests`; `KhaosSobolParityTests` with `khaos_sobol_oracle.R`).
- This commit: `PolynomialChaos.md`'s section and rows, `BayesianRegression.md`'s rows, the status line and this report.

**What is built.**
- *Sobol and the Jacobi measure from the coefficients,* under the product uniform law on each input's box: its declared domain, else its training range.
  - The decomposition is the moment route's own type: Var f = Σ b_α², each subset's main variance, each input's total Σ_{α_v ≥ 1} b_α², read from the same draws at O(M) per draw.
  - The Jacobi measure is E[x_v(1 − x_v)(∂_v f)²] = Σ b_α² α_v(α_v + 1). It is invariant under the affine map to raw units.
  - These are Lean SG-4(b)'s identities (`faca090`): `tensorTotal_eq`, `integral_jacobi_tensorPartial`, `tensorTotal_le` and `tensorTotal_eq_half`.
  - Refused for a hinge or categorical term and for a functional response.
  - It takes no measure, and its readout states the law it holds under.
- *The gradient readout:* ∂f/∂z per draw, case and location, in each predictor's own units. Each term's derivative comes through the design kernels, and the chain rule through the box. It holds pointwise, under no measure, for hinge and Legendre fits alike.
- *Moments under the other measures:* no code changed. The moment route already takes Legendre factors through their polynomial expansions and each measure's interval moments, so the slice gates it.

**Which measure each readout holds under,** as the coordinating session asked:

| Readout | Measure | Route |
| --- | --- | --- |
| `PolynomialChaosSensitivity.Sobol` | The product uniform law on the boxes | The coefficients, exact |
| `PolynomialChaosSensitivity.Jacobi` | The same law | The coefficients, exact |
| `BassFit.Sensitivity` on a Legendre fit | Any product measure, `MixtureMeasure` and `EmpiricalMeasure` included | The moment service, exact (a Student component by Gauss–Legendre, its error reported) |
| `BassFit.Gradient` | None: pointwise | The design kernels |

E8's C under a product measure is the moment service with derivative factors, and that is E8's own build. Under any measure, a dependent one included, the gradient averaged over the measure's points is the route, so E8's first consumer has an input that does not assume a product measure.

**The gates.**

| Gate | Reading |
| --- | --- |
| Sobol against the coefficients under the uniform measure | The moment route against the coefficient route, draw for draw, on an Ishigami fit (degree 8, order 3; 80 draws, M up to 18, 7 subsets): 6·10⁻¹⁴ of Var f, the statuses equal |
| Against `sobol_khaos` | Its own code on the engine's 80 draws, packed as khaos 2.1.1's `adaptive_khaos` object (multi-indices, coefficients, σ²): every index, its division by Var f + σ² included, to 2.2·10⁻¹⁶ |
| SG-4(b)'s identities against numeric integration | A 10³ Gauss–Legendre rule, exact at these degrees, on 8 draws each at degree 8, order 3 and at degree 2, order 2, from the draws' own surfaces and gradients: Var f, the main and total variances and the Jacobi measure to 4·10⁻¹⁵. Total ≤ ½·Jacobi at every draw and input, the least ratio 2.79 at degree 8. Equality to 10⁻¹² at the 4 (draw, input) pairs held to degree 1 |
| Moments against quadrature | The moment route on a Legendre fit (degree 4, order 2) against brute-force tensor sums over each measure, the density written out in the test. Under the empirical marginals, all 27,000 triples of training values: 3·10⁻¹⁴ of Var f. Under truncated mixtures (normal and Student components, both weightings), 40 density-weighted Gauss–Legendre nodes per input: 3·10⁻¹⁴, the Student moments' own quadrature error 7·10⁻¹² |
| Derivatives against central differences | The gradient readout against the predicted surface's central differences at 10⁻⁶ of each range, 40 cases: a Legendre fit on [−π, π]³ with that domain declared, so the chain rule shows, 3·10⁻¹⁰ relative; a hinge fit, 2·10⁻¹⁰ |
| Hinge BASS | The whole regression assembly at the slice's end, once, on HEAD `6bfb497` with the slice: 2,374 of 2,374, the pinned hashes and the SMC baseline among them. The commit's selection with the pinned tier: 1,678 of 1,678. The draws of every existing fit are unchanged: `BassSurface.Evaluate` keeps its arithmetic, and the pinned tier holds |

**The cost** (wall clock, one run, indicative): on the comparison's 80 draws, the moment route took 18.1 ms and the coefficients 0.7 ms, each with its setup.

**Found here.**
- *Degree 1 everywhere cannot hold interactions.* `LegendreBasis.MaxDegree` is a term's total degree, which must be at least its order. So SG-4(b)'s equality case is checked wherever a draw holds an input to degree 1 alone, at degree 2 and order 2.
- *The reference's own cancellation.* A near-empty draw's variance as E[f²] − (E f)² loses every digit. The test's references are two-pass sums about the mean.

**Decisions beyond the brief.**
- *The coefficient route returns the moment route's type,* so the two compare field for field and a consumer reads either. The moment route stays `BassFit.Sensitivity`'s default, since it serves every measure and every term type.
- *The gradient readout is general,* hinge and Legendre alike, and pointwise. It is the input E8 and E25 can take under a dependent measure.
- *`BassFit.Sensitivity`'s doc* now says the default law is on each input's box, its declared domain, else its training range. Slice 2 had noted this as left.

**Left unbuilt, and why.**
- *E8's derivative moments* (C = E[∇f∇fᵀ] by the moment service with derivative factors): E8's own build.
- *The gradient of a basis fit* (`BassBasisFit`), its components' gradients combined through the output basis. E8's functional reading, co-active blocks contracted with the output Gram, decides that shape, so it is recorded for E8.
- *Slice 3,* waiting on the coordinating session's note.

### Slice 5 addendum: SMC's flags (chip, 2026-11-08)

**The decision.** Flag SMC, option 2 of the choice the trap's measurement put. The coordinating session made the call on 2026-11-08: the user was unsure and left it to that session, and can still reverse it. Nothing else of the choice is built (below).

**Landed** on `main`, nothing pushed: `42739ac`, the flags with their gates and the trap harness's flag gate; and this commit, the docs and this report.

**What is built.** Two flags on BASS's SMC fit generally, read from what the run already holds. They change no draw, weight or normalizer.
- *A jump over a collapse* (`BassSmcReplicate.JumpedOverACollapse`, counted in `BassSmcReadout.JumpedReplicates`). The last step to b = 1 left from b below ½ (`JumpBelow`), over a population of at most 5% of N distinct structures (`CollapsedShare`).
  - The run counts the population's distinct structures once, before the step that reaches b = 1, and a checkpoint carries the count. A fit resumed after that step still reads the flag.
  - The flag is null short of b = 1, and where the count was not taken.
- *Disagreeing replicates* (`BassSmcReadout.Overdispersed`, mirrored in `BassEvidence.Overdispersed`). The replicates' ln Ẑ between-to-within ratio is above 3 (`BassEvidence.OverdispersionThreshold`), pooled as the evidence pools them.
  - The evidence's error is then the replicates' spread already, by the pooling rule's larger error. The flag says the genealogy errors are not to be read.
  - It is null with one replicate, or where a genealogy error is NaN (a resampling scheme other than multinomial). The readout then says so in words rather than reading clear.
- *In words* (`BassSmcReadout.Flags`): one line per jumped replicate, one for the disagreement, or one saying it could not be read. A summary or adapter prints them; presentation stays values only.
- *The thresholds were fixed before any healthy fit was read,* from the trap's measurement.
  - Trapped replicates jumped from b 0.003 to 0.35 over 1 to 9 structures of 500. Those that reached the posterior left from b 0.98 and up.
  - If each genealogy error is right, the ratio squared is about χ²_{R−1}/(R − 1). That passes 9 with probability 0.27% at R = 2, and below 10⁻⁵ at R = 4.

**The gates.**

| Gate | Reading |
| --- | --- |
| Draws unchanged | A new pinned SMC baseline (`BassBitIdentityTests.Smc_KeepsTheBaseline`) hashes the draws, the weights, each replicate's ln Ẑ and genealogy error, and its stages' b and sweeps. Its values were read on `fc39acc` before the flags were built (hinge `3440c6398ce3c0d5`, KHAOS `a51939ff36591784`), and they hold with the flags |
| The trap | Production fits on the trap fixture, N 500 at 5 sweeps a stage (`KhaosSmcTrapQualificationTests.Flags_OnTheTrap_AndNoneOnHealthyFits`). On the raw ridge, 32 replicates: 2 reached the posterior and neither jumped. 24 of the 30 trapped replicates jumped over a collapse; the other six grew part of the structure (ln Ẑ up to −334.8) and left from b 0.99 and up. The disagreement reads those six: the whole fit's ratio is 30.0, and all 8 groups of four are flagged one way or the other. On the modified g-prior, 16 replicates: all 16 reached σ² within twice the chains' and none jumped. Their ln Ẑ disagree at 9.2 times their errors, and all 4 groups of four are flagged |
| Healthy fits | No flag fires. GBASS's t(5) and skewed GH evidence fixtures at the adaptive-mutation measurement's size (32 replicates of 500): ratios 0.96 and 0.85, no jump, no group of four flagged; at their evidence gate's size (4 of 2,000): 1.41 and 0.97. The SMC fits against R's enumeration (six normalizers and two posteriors): 0.23 to 1.20. KHAOS's exact-law fixture, three fits: 0.38 to 1.44. `BassSmcTests`' fit: 0.19. Every last step left from b 0.948 and up. The collapse share alone does not separate: these small spaces held 2.1% to 8.8% of N distinct before their last step, against the trapped runs' 0.2% to 1.8%, so the jump carries the separation in every reading. A small space whose data barely inform could leave from b < ½ over so few structures and read as jumped. No reading here did, so the threshold was not moved |
| The rule's edges | Constructed replicates (`BassSmcTests`): from b exactly ½, no jump; 25 of 500 structures is collapsed, 26 is not; short of b = 1, or uncounted, null. Two replicates 60 nats apart against errors 0.1 and 0.8 read 74.4 times; two trapped alike read neither ratio nor flag, only the jumps; a NaN genealogy error leaves the ratio unread. A one-replicate fit reads null in both places and says so in words. The count survives every resume, a complete fit's resume included |
| Hinge BASS | The whole regression assembly on HEAD `fc39acc` with the flags, once: 2,366 of 2,366, the pinned hashes and the new SMC baseline among them. The SMC parity fits against R: 12 of 12 |

**What the flags see beyond KHAOS.** On Friedman's band (stage 8's qualification) the replicates' ln Ẑ spread 6.5 to 8.4 nats against genealogy errors of at most 0.83. The disagreement flag reads that, as it should.

**What they cannot see.** A wrong region that every replicate reaches alike, by a schedule that never jumps. Its replicates agree, and no step leaves a collapse. Only a reference outside the run reads it: another sampler, or an enumeration where one exists.

**Recorded, not built.** Each would change what an SMC fit draws, and each waits on its own decision.

| Remedy | What the measurement showed | What it would cost |
| --- | --- | --- |
| A cap on the step, b′ ≤ r·b | It reaches the basin but not the evidence. At r = 1.25, 11 of 16 replicates reach σ² within twice the chains', against 3 of 32 uncapped, at twice the time. At r = 1.1, 8 of 16, at 3.5 times. They end at M 7 to 10 against the posterior's 10.15, and their ln Ẑ spreads over 49 nats with errors under 1 | An `SmcSpec` field and one bound in the schedule's step, opt-in. Its draws would differ from the uncapped fit's wherever it binds. *Built opt-in 2026-11-09 at the user's word (`SmcSpec.StepCap`, `0ffc0e0`): with E4's residual birth, 16 of 16 replicates reached the posterior ([the proposal-kernels brief](proposal-kernels-brief.md#reports)).* |
| A birth proposed from the data's residual, through F-X4 | Unmeasured. Births accept about 2% on the trap fixture, which limits the capped runs and MCMC's climb from empty alike | A new `IBassProposalKernel` member that weights each term of 𝒜 by its fit to the current residual, with exact forward and reverse probabilities: O(\|𝒜\|·n) per birth to score 𝒜, or a score kept and updated per accepted move. Its gates are the exact law by enumeration and Geweke, then this measurement again. Slice-sized |
| A reference nearer the posterior | Unmeasured | A pilot MCMC run, and a reference over structures with an exact density and full support (the pilot's structures mixed with the prior). A tempering path from it to the posterior in the SMC core and BASS's target, with ln Z = ln Z_ref + the ratio. Larger than the birth |

### Slice 5: the engine's reach (chip, 2026-11-08)

At the user's order of 2026-11-08 ("Yes, 5 first"): slice 5, then 4, then 3.

**Landed** on `main`, nothing pushed, the commit messages without attribution lines:
- `98b0e7e`: the gates below. Three small lifts in the hinge tests let them run on Legendre terms:
  - the ladder's per-rung batches (`BassTemperingTests.RungBatches`);
  - the flat histogram's check against the enumeration (`BassFlatHistogramTests.MatchesTheEnumeration`, its node count a parameter);
  - the pointwise refit with its prior given (`BassPointwiseTests.Refit`).
  The term enumerations and the structure key also read Legendre terms. GBASS's per-rung Geweke takes Legendre rows (`legendre`). The SMC trap's harness (`KhaosSmcTrapQualificationTests`) has its table (`khaos-smc-trap-replicates.csv`).
- This commit: `PolynomialChaos.md`'s section and rows, `BayesianRegression.md`'s KHAOS row, the status line and this report.

**No engine code changed.** A Legendre term is a part type behind F-X4 and F-X3, so the ladder, SMC, the flat histogram, the pointwise density, the candidate density and GBASS's noise serve it through the design and the target. The slice is its gates.

**The gates.** The fixture is slice 0's exact-law fixture (p 4, d_max 3, q_max 2, M ≤ 2, n 40, 466 structures) with σ² ~ IG(2, 0.5), under the raw ridge at τ held 0.1 and the modified g-prior with g₀² ~ IG(2, 80) drawn.

| Gate | Reading |
| --- | --- |
| The law at each rung | Four ladders of 40,000 scans each, from the ladder's own replicas, against the enumeration's tempered law: the count law at every rung within 5 SE + 0.004, and E[ln τ] within 5 SE + 0.01. Raw ridge at b = 1, 0.4 and 0, and with the reference rung at 1, 0.4, 0.1 and 0. Modified g-prior at 1, 0.4 and 0, and with the reference rung at 1, 0.3, 0.05 and 0, E[ln τ] −3.5119 against −3.5095 (± 0.0015) through −3.9585 against the prior's −3.9592 |
| SMC and the ladder, ln p(y) | Against the enumeration's ln p(y): SMC (2,000 particles, 4 replicates) −72.2965 ± 0.0144 against −72.2910 on the raw ridge, and −69.6168 ± 0.0167 against −69.5869 on the modified g-prior. The ladder (reference rung, self-sizing) −72.2579 ± 0.0359 and −69.6057 ± 0.0426. SMC under the adaptive count: −69.599 ± 0.023 on the modified g-prior. At 32 replicates of 8,000 particles, run to check a 2.75 SE gap in P(M = 0) at 8 replicates, the count law sits within 1.2 SE of the posterior with g₀² held and 0.72 SE with it drawn, and ln p(y) within 0.8 SE. The gap was the 8 replicates' error. Every count law within 5 SE + 0.01 |
| Rungs at once | A tempered KHAOS fit with its rungs stepped at once (`RungWorkers` −1) is the fit stepped in order, hash for hash: e5b1758628567579 both |
| The flat histogram | `MatchesTheEnumeration` under each prior: on the modified g-prior (12 nodes) the learned ln w is 0.039 from SAMC's limit, P(M) within 0.31 SE and every check within its allowance; on the raw ridge (6 nodes) 0.020, P(M = 0) at 1.87 SE. At a held τ the density of states is F's 466 point masses, which leave SAMC's limit undefined on 12 nodes |
| The pointwise density | Both contents, on draws from 13 structures with g₀² drawn and 4 at τ held: the retained content against the Normal at the explicit uncentred design (4·10⁻¹⁶), and the collapsed content against explicit refits under KHAOS's prior written out, τI on the uncentred coefficients or τ·D_g⁻¹ΨᵀΨD_g⁻¹ with Ψ = [1, B] (1.4·10⁻¹⁴). Under g₀² drawn a structure is met at several τ, so each τ takes its own leave-out quadratics |
| The candidate density | At the training cases it is the pointwise column, under the Gaussian and the Student t noise: 0 relative |
| GBASS's noise on Legendre terms | GBASS's per-rung Geweke (`BassNoiseReachQualificationTests`), 400,000 marginal-conditional draws against 600,000 successive steps, on Legendre terms of degree up to 3 under the modified g-prior with g₀² ~ IG(2, 2) drawn: the Student t at b = 1 and 0.4, worst \|z\| 2.92 and 2.07; the quantile at b = 1, interwoven, 1.74; the skewed horseshoe with its move along β·v, its default, at b = 1 and ¾, 0.96 and 2.82 (mean ln v, its slow coordinate, as on hinges). No move refused |
| Conformal coverage on a Legendre fit | The full conformal Bayes set on a Gaussian KHAOS fit (KHAOS's defaults, degree 8, order 2) under the contaminated noise of `BassCandidateTests`: 20 fits of 99 cases, 50 held-out cases each. Coverage 0.878 against the nominal 0.900 (se 0.011, 1.9 SE low); mean width 1.08; no set open at the grid's ends |
| A Student t KHAOS fit | A t(3) KHAOS fit against a Gaussian one under 0.25·t₂ noise, n 150, the posterior mean surface against the truth at 200 held-out cases. Over 16 sets the t fit is better on all 16, mean RMS error 0.088 against 0.188. A first run of 6 sets read 5 of 6 and a paired gap of 1.3 SE, short of the 2 SE it was gated on. t₂'s tails make the gap heavy-tailed: one set reads 0.82 against 0.09. So the gate is now the sign, at least 12 of 16, on fresh seeds |
| A quantile KHAOS fit | q = 0.8 on heteroscedastic noise (0.1 + 0.4x₁)ζ, 4 sets of 300. 0.768 of 1,200 held-out cases fall below the fitted surface (± 0.012, within 3 SE + 0.02). Its RMS distance from the true 0.8-quantile surface is 0.079, against the Gaussian mean surface's 0.289 |
| Hinge BASS | The whole regression assembly at the slice's end, once, on HEAD `477834c` with the slice's tests: 2,361 of 2,361, the pinned hashes and every fingerprint among them. The hinge classes whose helpers were lifted, slow methods included: 310 of 310. The slow tier is rewritten from that profile; it now holds the three KHAOS reach methods above 10 s, and timing under load moved a few hinge methods across the line both ways |

**SMC on KHAOS: the trap, measured. Its remedy is the user's open choice.** The harness is `KhaosSmcTrapQualificationTests`, on the adaptive-mutation measurement's fixture and seeds. Its engine arms reproduce that measurement's replicates ln Ẑ for ln Ẑ.

| Arm, N 500 | Reached (σ² < 0.449) | ln Ẑ, reached | ln Ẑ, trapped | Seconds, reached and trapped |
| --- | --- | --- | --- | --- |
| Engine schedule, 5 sweeps | 3/32 | −316.5 ± 4.8 (−319.6 to −311.0) | −349.3 ± 3.0 | 8.4, 4.7 |
| Engine schedule, 20 sweeps | 2/16 | −291.2, −283.1 | −349.86 ± 0.09 | 30.5, 14.6 |
| b′ ≤ 1.25b, 5 sweeps | 11/16 | −313.4 ± 16.6 (−331.1 to −282.1) | −331.0 ± 4.0 | 10.2, 9.4 |
| b′ ≤ 1.1b, 5 sweeps | 8/16 | −318.0 ± 5.0 (−327.7 to −310.9) | −327.6 ± 4.8 | 16.5, 15.9 |
| Modified g-prior, engine, 5 sweeps | 16/16 | −221.1 ± 11.1 (−239.1 to −196.3) | none | 30.0 |

- *The mechanism is the schedule's jump.*
  - Below b ≈ 0.005 a term costs more than it gains, about 14 nats under this raw ridge, so the tempered law sits at near-empty structures. The population collapses there: at the stage before the jump it holds 1 to 9 distinct structures of 500 particles, at M ≈ 0.
  - They all reweight alike, so the conditional ESS sees no spread, and the next b is 1, from b between 0.003 and 0.35.
  - ln Ẑ is read at that step. The final stage's sweeps start growing terms (mean M 0.2 to 1.0 after it and two more stages), but the evidence is already spent.
  - At b = 1 itself there is no trap: greedy forward selection from M = 0 climbs monotonely, and the ladder's cold chains read the long chains' M and σ².
- *The cap reaches the basin but not the posterior.* The capped runs end at M 7 to 10, against the posterior's 10.15 and σ² 0.2245. Their ln Ẑ spreads over 49 nats while each replicate reports an error under 1. Births accept about 2% of the time, and that limits them.
- *KHAOS's own g-prior is not trapped,* at about 9 nats a term. All 16 replicates reach σ² 0.19 to 0.21 at M 13.5 to 14.9. Its ln Ẑ still spreads over 43 nats, sd 11, each replicate reporting an error under 1.
- *The ladder has the same transition and flags it.*
  - Four chains, a reference rung, self-sizing, 5,000 scans after 5,000.
  - The cold chains read M 10.13 to 10.16 and σ² 0.2246 to 0.2267.
  - One pair in every chain rejects every swap, and there are no round trips. The chains' stepping stones read −296.2, −266.8, −246.8 and −205.6.
  - The pooled −253.9 ± 19.0 is flagged saturated, between-to-within 28.
- *What sees the trap.*
  - The last step left from b < ½: 41 of 43 trapped replicates, 0 of 5 reached.
  - `LeastMoved` under 0.1 over the last three stages: 40 of 43, and 1 of 5 reached.
  - The mean log posterior's drift over two extra stages at b = 1: 0 of 29 at 5 sweeps.
  - The pooled evidence's between-to-within ratio, which SMC fits already compute over their replicates: it sees the dispersed evidence the jump does not.
- *The choice,* sent 2026-11-08 to the coordinating session for the user. Nothing is built until the user picks.
  1. Refuse SMC on a Legendre basis. Builds nothing, and loses the g-prior's working posterior.
  2. Flag it on the jump and on the between-to-within ratio, recommended. A small build that leaves the draws unchanged, and it covers any SMC fit.
  3. Fix it. A step cap is cheap but measured insufficient for the evidence. A residual-informed birth through F-X4 targets the 2% acceptance; it is slice-sized and unmeasured. A reference nearer the posterior means a pilot and a new tempering path; it is larger and unmeasured.

**Found here.**
- *A noise's reference beside the zero centre fails at the fit, not at validation.* `BassNoiseReference.Moments`, the default, reads a level from the response, and `BassResolved` refuses it against the zero centre once the fit runs. KHAOS's fits take the zero centre, so under a noise they state `BassNoiseReference.None` with Jeffreys σ². What the reference should read under the zero centre, or whether validation should refuse it with the config's path, is GBASS's call. Recorded, not changed.
- *SMC's reported evidence error under-covers on KHAOS fixtures of this size,* under either prior (above).

**Decisions beyond the brief.**
- *The noise and the modified g-prior.* The prior reads Ψ as evaluated, not the latents' weights, as `khaos`'s S₀ does. The noise's unweighted prior Gram already serves it.
- *The test helpers were lifted, not copied.* The hinge tests keep their bodies and seeds.

**Left unbuilt, and why.**
- *The trap's remedy:* the user's choice.
- *Slices 4 and 3,* in the user's order, each after the coordinating session's verification.

### Slice 2: the modified g-prior (chip, 2026-11-08)

**Landed** on `main`, nothing pushed, the commit messages without attribution lines:
- `f289af3`: the prior on τ stated on τ itself (`TauPrior.Gamma`); the prior-precision builder (`BassPriorPrecision`, E21) with slice 1's raw ridge folded into it; KHAOS's modified g-prior (`ModifiedGForm`) through the builder; its config and resolution; two overflow fixes in the shared g-prior code (below); and the gates.
- This commit: `PolynomialChaos.md`, `BayesianRegression.md`, the notation ledger and this report.

GBASS committed `BassConfig.cs` in `211d768` while the slice was in its gates, so the config's piece went into the working tree before the slice commit. Until then it lived in the build snapshot only.

**What it is.**
- *The model.* `adaptive_khaos_gprior`'s as written: (β₀, β) | σ², g₀² ~ N(0, σ²S₀) on the columns as evaluated, the intercept's column of ones among them, with S₀ = g₀²·D_g(ΨᵀΨ)⁻¹D_g, D_g = diag(1, g₁, …, g_M) and g_m = (1 + q_m(d_m + q_m − 2))^(−ζ/2). `khaos`'s defaults are ζ = 1 and g₀² ~ IG(0.001, 1,000).
- *The builder.* The design held a flag per form, the normalized ridge's and slice 1's raw ridge's. It now holds a `BassPriorPrecision`, which says:
  - whether the prior stands apart from the Gram;
  - what each term's prior-side column is;
  - the pencil's D_α and D_β;
  - whether the angle gate reads the prior's Gram.

  It has four members: the g-form, the normalized ridge, the raw ridge (D = I) and the modified g-prior.
- *The modified g-prior is the g-form's pencil on scaled columns.* In the design's coordinates it is the g-form's prior-apart pencil on the prior-side columns scaled by 1/g_m: D_α = N, δ = z̄_u/g − z̄_w and D_β = D_g⁻¹G_uD_g⁻¹. No new pencil code was needed, and the angle gate reads it as it reads the g-form's.
- *g₀²'s step is τ's.* τ = 1/g₀², and IG(a, b) on g₀² is Gamma(a, rate b) on τ. τ's slice step on ln τ (D04) is then g₀²'s on ln g₀², reflected, against the conditional with β and σ² integrated, at O(p) per g₀² after one decomposition per structure. The brief's "one `SliceSampler` step on ln g₀²" needed no block of its own.
- *ζ = 0* is the g-form with the intercept inside under the zero centre.

**The gates.**

| Gate | Reading |
| --- | --- |
| g₀²'s conditional against the joint | At six structures and g₀² from 0.5 to 10⁶, for ζ = 0, 1 and 2.5: the target's log marginal against a dense multivariate t from the columns alone, y \| g₀² ~ t with scale I + B S₀ Bᵀ, 1.3·10⁻¹⁰ after one constant; the conditional τ's step reads against the joint, per structure, 1.3·10⁻¹⁰ |
| ζ = 0 against `CollapsedGPrior` | The pencil against the g-form with the intercept inside, both with g₀² drawn: 466 enumerated structures at five τ each, 4·10⁻¹⁶ relative |
| The exact law, g₀² integrated | Each of the 466 structures' weight by Simpson's rule over ln τ on the oracle's nodes, against `khaos_exactlaw_oracle.R`: 1.6·10⁻¹⁵ in P(M) and the 30 inclusions, 1.7·10⁻¹⁵ in g₀²'s posterior CDF at 300, 1,000 and 10⁴, on both boxes. The production chain with g₀² ~ IG(0.001, 1,000) drawn against it over 400,000 steps: worst \|z\| 1.46 |
| Geweke, g₀² drawn | 12 cases, M ≤ 3, σ² ~ IG(3, 2), g₀² ~ IG(2, 2): worst \|z\| 2.37 and 1.52 (ζ = 1, two seeds), 2.14 (ζ = 0), over nine functions, ln g₀² among them |
| `adaptive_khaos_gprior` readouts | Corrected, degree 8, order 3; R 6,000 iterations with 3,000 kept, the engine two chains of 10,000 after 10,000. On the raw inputs with [0, 1] declared: M 14.66 against 13.87, σ² 0.201 against 0.203, g₀²'s median 3,137 against 3,256, the surfaces 0.068 apart, total indices within 0.003. On the rescaled inputs: M 20.23 against 19.08, σ² 0.221 against 0.236, the surfaces 0.074 apart, the indices within 0.003. Long runs on the declared box settle the gap in M: eight engine chains of 40,000 read M 15.06 (sd 0.36), σ² 0.197, g₀²'s median about 3,130; two `khaos` runs of 40,000 read M 15.06 and 15.03, σ² 0.198 and g₀²'s median 3,125 and 3,115 |
| Slice 1 under the builder | The raw ridge's exact-law chain and Geweke read as before to the printed digit (\|z\| 1.79, 0.71; 2.15, 1.39, 1.31) |
| Hinge BASS | The whole regression assembly at the slice's end, 2,253 of 2,253, the eight pinned hashes and every fingerprint among them; the commit's selection on the later HEAD, 1,849 of 1,849 |

**Parity with `khaos`.** The engine's model is `khaos`'s as written. Two parts differ, each a slice 0 defect:
- *g₀²'s kernel.* The engine takes the brief's slice step on the derived conditional. `khaos` takes an independence step from a Laplace approximation to a conditional without the data term.
- *The structure ratio.* The engine's carries ½Δln|S₀⁻¹|, which `khaos`'s leaves out.

**Refusals.** Each is refused at resolution with its reason:
- the modified g-prior on a hinge basis, since its factor reads a Legendre term's order and degree;
- the modified g-prior under weights whose prior Gram is the weighted one, which moves with W.

**The cost row** (`KhaosCostTests`, the cycle clock, indicative): at the base a modified g-prior sweep is 343 µs at M 20, beside the raw ridge's 363 µs, the normalized ridge's 35 µs and the hinge fit's 21 µs in the same run. Each candidate is decomposed afresh, as under the raw ridge. At p 20, q_max 5 it is 58 µs at M 6.

**Found here.**
- *Two overflows in the shared g-prior code, fixed.* `khaos`'s IG(0.001, 1,000) on g₀² puts about half its mass on τ below 10⁻³⁰⁰, where the density is finite. Two computations broke there:
  - The spectral terms' ln(1 + bλ/τ) (`RidgeColumns`, every pencil and spectrum) read +∞ once bλ/τ overflowed.
  - Goldberg's ln(1 + x) (`LogFunctions.Log1p`) overflowed in its product ln(u)·x above x ≈ 2.5·10³⁰⁵.

  The log marginal read −∞ at such a τ, and a band of −∞ for u between about −709.8 and −703.2 stopped τ's slice step below it. On the cost row the chain sat at u ≈ −711 and M 0 for 4,000 sweeps, every birth refused. Each computation now takes ln bλ − ln τ or ln u where it overflowed, and keeps its operations everywhere else, so every finite reading stands (the whole assembly, with its pinned hashes). From u = −721.7 the chain now reaches u ≈ −31 by sweep 100 and M 20 by sweep 1,000. A test pins the marginal's ½-per-direction slope through the old wall, and the climb.
- *The box matters less under the g-prior.* On Ishigami the g-prior reads M about 15 and σ² 0.20 on the declared box, and M about 20 and σ² 0.22 on the rescaled one; the raw ridge read 10 and 0.23 against 6 and 0.95. A term costs about ½ln(1 + g₀²g_m²) + ln|𝒜| ≈ 3 + 6 nats here against the raw ridge's 14, so the rescaled box's denser expansion still pays.
- *The readout gate's M gap was Monte Carlo.* R's single chain of 3,000 kept draws read M 0.8 below the engine's. Long runs of both agree to 0.03.

**Decisions beyond the brief.**
- *g₀² is τ inverted, through τ's existing block,* with a Gamma on τ itself, rather than a g₀² block. The draws record τ, and g₀² = 1/τ.
- *The builder holds the prior side, not the target's dispatch.* The target already reads a structure's fit, whether a factor, a spectrum or a pencil, and needed no change. The g-form's and the normalized ridge's fast forms stand wherever the prior is not apart, so their draws are unchanged.
- *The raw ridge's prior-apart pencil folded into the builder* as its unit member, its draws unchanged.
- *ζ is finite and not negative,* and the factor is the paper's for every ζ, so ζ = 0 needs no special case.

**Left unbuilt, and why.**
- *The modified g-prior on hinges.* Its factor reads a Legendre term's order and degree, and a hinge analogue is not KHAOS's.
- *The weighted prior Gram for the modified g-prior* (above).
- *Tempering and SMC under KHAOS's priors* are not gated, as for slice 1's raw ridge. The pencil serves them as it serves weighted g-form fits.
- *Slices 3 and on,* waiting for the user.

### Slice 1 addendum: the input box's effect is the posterior's, not a trap (chip, 2026-11-08)

**The question.** The user's go-ahead of 2026-11-08 ("Check, then slice 2") asked whether the readout fixture's chains hold the posterior or a trap both samplers share. On the inputs rescaled to [0, 1] both read M about 6.2 and σ² about 0.95, against M about 10.1 and σ² about 0.23 on the raw inputs. Both accept about 2% of births, so agreement alone could not tell.

**Landed** on `main`, nothing pushed: this commit, with `KhaosRichStartTests` behind `BASS_QUALIFICATION`.

**The method.** On slice 0's Ishigami readout fixture under the raw ridge (degree 8, order 3, khaos's τ² = 10⁵, σ² Jeffreys), in two boxes: the raw inputs with [0, 1] declared, and the inputs rescaled to [0, 1] exactly with no domain.
- The collapsed log posterior, the log marginal plus the structure prior with λ integrated, is scored for three states: greedy forward selection from M = 0; its best-improvement descent by deaths, swaps and births; and, in the rescaled box, the raw structure's image under the affine map. The image is the closure of the raw structure's multi-indices, and it spans the raw columns exactly.
- Four chains run from each start, 10,000 steps after 10,000, traced in log π, M and E[σ² | S].

**The reading: real.**

| Box | State | M | log π | E[σ² \| S] |
| --- | --- | --- | --- | --- |
| Raw, [0, 1] declared | Greedy and its descent | 10 | −214.27 | 0.232 |
| | Chains from M = 0 and from the descent | 10.14 | −214.27 (median) | 0.2275 |
| Rescaled | Greedy and its descent | 6 | −252.11 | 0.958 |
| | The raw structure's image | 27 | −385.05 | 0.156 |
| | The image's descent | 6, the same structure | −252.11 | 0.958 |
| | Chains from M = 0 and from it | 6.42 and 6.37 | −252.11 (median) | 0.902 and 0.908 |

- *The rich state fits better and scores far lower.* The image beats the mode by about 109 nats of likelihood (E[σ²] 0.16 against 0.96 at n = 120), yet sits 133 nats below it. Each term costs about 14 nats under the raw ridge, ln|𝒜| ≈ 6 and ½ ln(nτ²) ≈ 8, not the 9 the go-ahead estimated.
- *Every route lands on the mode.* The image's descent returns to the sparse mode, and chains from either start agree to 0.05 in M.
- *The cause is Ishigami's parity about ½.* The raw structure is parity-pure: odd degrees in x₀ (sin z₁), even in x₁ (sin² z₂) and in x₂ (z₃⁴). A box centred on ½ keeps that parity, since ψ_α(½ + c(x − ½)) spans only the ψ_β with β ≤ α of α's parity, and a shifted box breaks it. Two declared boxes of one width, 1.04, separate the two:

| Declared box | Chains from M = 0 | E[σ²] |
| --- | --- | --- |
| Centred on ½, [−0.02, 1.02] | M 10.16, log π −212.55 | 0.227 |
| Shifted by 0.01, [−0.03, 1.01] | M 7.69, log π −246.7 | 0.646 |

**What it changes.** The readout gate stands as a posterior reference in both boxes, and slice 2's readouts take the same fixtures. A Legendre fit's box is part of its model: on a function with a symmetry about the box's centre, a box that is off-centre by 1% of its width costs the expansion its sparsity. That is the case for declaring the domain.

**Found on the way.** Greedy forward selection is the weaker search on these surfaces. On the centred box it stops at M 7 and log π −221.59, below the chains' −212.55. The check's states are the descent's local optima, and the check reads them only beside the chains.

### Slice 1 addendum: declared input domains (chip, 2026-11-07)

**The decision.** The slice 1 report left open whether a continuous predictor may declare its domain. The user leaned toward it ("predictor able to declare its domain") and chose to build it, on 2026-11-07, which the main session confirmed ("i approved the declared domains thing").

**Landed** on `main`, nothing pushed, the commit messages without attribution lines:
- `3eb56e1`: `InputDomain` on `ContinuousPredictor`, the box in `RangeScale`, the fingerprint, the sensitivity measures' default and flag, the resolution's log, and the tests.
- This commit: `PolynomialChaos.md`, `BayesianRegression.md` (with GBASS's GBMARS row for its slice 6, at its request) and this addendum.

**What it is.**
- *The declaration.* `InputDomain(Lower, Upper)` on `ContinuousPredictor`, output variables included. Refused: an end that is not finite, a lower end not below the upper, and a value outside the domain.
- *The box.* `RangeScale` scales by the domain when one is declared and by the training range otherwise (`Min`, `Max`), and keeps the training range (`TrainingMin`, `TrainingMax`). `Contains`, which flags prediction cases and calibration bounds, reads the training range, since past the data a fit extrapolates whatever the box.
- *The fingerprint* adds a domain after its column's values, and only when one is declared, so data that declare none keep their fingerprint.
- *Sensitivity.* A measure left undeclared is uniform over the box, and `ContinuousMeasure.Extrapolates` reads the training range in the box's units.
- *The log.* Resolution logs each domain, and a Legendre fit's line says when a predictor declares none.

**The gates.**

| Gate | Reading |
| --- | --- |
| The exact law | With [0, 1] declared on slice 0's raw rows, the engine's π against `khaos_exactlaw_oracle.R`, 3.4·10⁻¹⁵. With the inputs drawn into [0.05, 0.95] and [0, 1] declared, π against a dense marginal on ψ at the raw values, 3.4·10⁻¹⁵; the training range's law differs from it by up to 0.019 |
| Readouts | Ishigami's raw inputs with [0, 1] declared, against corrected `adaptive_khaos_ridge`: M 10.15 against 10.15, σ² 0.2275 against 0.2278, the surfaces 0.005 apart, total indices within 0.001. The rescaled case, kept beside it: M 6.22 against 6.13, σ² 0.929 against 1.005, the surfaces 0.035 apart, the indices within 0.004 |
| Hinge BASS | 2,263 enumerated structures under the g-form with a domain of [−0.5, 2] around data on [0, 1]: the same structures admitted, π unchanged to 7.8·10⁻¹⁶. The whole regression assembly, 2,179 of 2,179 |
| The rest | Refusals, the box and the flags, the fingerprint, the measure and the log (`InputDomainTests`) |

**Found here.** The readout gate's engine chains ran 3,000 steps after 3,000. On the raw inputs they sat at M 9.78 and σ² 0.263, a σ² ratio of 1.155 against the gate's 1.15. Eight engine chains of 30,000 read M 10.13 and σ² 0.2282, and two `khaos` runs of 30,000 read M 10.14 and 10.15 and σ² 0.2274 and 0.2275. With about 1.4% of births accepted, the short chains had not finished growing from M = 0. The engine's chains now run 10,000 after 10,000, about two seconds.

**Decisions beyond the brief.**
- *The domain is the data's, not the config's.* It is a fact about the input, it travels with the predictor to every fit, and the fingerprint already covers the data.
- *The flags read the training range.* Past the data a hinge fit extrapolates beyond its last knot and a Legendre fit has no observations, whatever the box. A case outside the declared domain is flagged the same way, not refused.
- *Hinge fits take a domain too.* Their law is unchanged under the g-form and the normalized ridge, since the knots are the data's and the prior reads each column up to its scale. Two readings move: D06's column mass, which reads hinges normalized over the box, and an undeclared sensitivity measure, which spans it.

**Left as it was.** `BassFit.Sensitivity`'s doc still says "the training range's uniform by default". `BassFit.cs` holds another session's uncommitted edits, so the line waits for whoever commits that file.

### Slice 1: the Legendre basis under the ridge (chip, 2026-11-07)

**Landed** on `main`, nothing pushed, the commit messages without attribution lines:
- ``fa05048``: `ShiftedLegendre` and `PoissonBinomial` in `Maths.Distributions`, with their tests. The engine's primitives also join slice 0's R comparisons.
- ``b408bf4``: the Legendre part, KHAOS's term space and its measure, F-X4 (`IBassProposalKernel`) and KHAOS's proposal, the sweep reading its moves through the contract, and their tests, with the engine's `make_weights` and forward law against `khaos`.
- `baf69f9`: KHAOS's raw ridge in the design, the config (`BassStructureSpec.Legendre`, `BassCoefficientPrior.RawRidge`) and its resolution, and the engine-level gates: the exact law, the chain against it, Geweke, the config's refusals, the cost row and the readouts against corrected `adaptive_khaos_ridge`.
- This commit: `docs/Maths/PolynomialChaos.md`, pointers from `BayesianRegression.md`, and this report. At GBASS's request the same commit brings `BayesianRegression.md`'s GBMARS row to its slices 4 and 5.

**The waiting.** GBASS slice 5 held `BassConfig.cs` uncommitted until `37137f2`. Until then the config's part was applied to HEAD's copy inside the chip's build snapshot only, and the working tree's `BassConfig.cs` was left alone. `BassPriorSampler.cs` needed no change: its prior draws go through `BassStructurePrior.DrawTerm`, which serves the Legendre part.

**How the Legendre part entered `BassTerm`** (the brief left it to the build). `BassTerm` takes a fourth part, `LegendreFactor(Variable, Degree)`, canonical by variable, with `BassPart.Legendre` and `BassStructurePrior.PartCount` = 4. A hinge term's hash, equality and columns are unchanged: the hash adds the part only when it is nonempty, and `LegendreProductKernel` sits last among the case-side kernels and leaves the column alone for a term without the part. The contract gained one member, `IDesignKernel.ReadsVariable`, with a default, so a kernel written outside it (the tests' directional hinge) still compiles. `BassBasis.CaseDerivative` now takes the slope of the part that reads the variable. Its old rule took every differentiable part's derivative, which would have zeroed a Legendre term's through the empty hinge part.

**The gates** (each number from this chip's runs; the docs page carries the table):

| Gate | Reading |
| --- | --- |
| Primitives | ψ_α orthonormal by quadrature to degree 40 (10⁻¹²), zero at the Gauss–Legendre nodes to 64, endpoints to 60, Bonnet's identity on a 4,000-case sweep to 60 (10⁻¹²), the lanes bit for bit; against R's `ss_legendre_poly` (10⁻¹²). The Poisson-binomial against brute force over 400 sweep cases across twelve decades (10⁻¹² relative) and the binomial at 2,000 coins (10⁻⁹); against slice 0's enumeration and R's tail (10⁻¹⁴). \|𝒜\| three ways at six sizes, and at p 20, d_max 15, q_max 5 exactly |
| The law, exactly | p 4, d_max 3, q_max 2, M ≤ 2: 466 structures. π against a dense marginal from the columns alone, 6·10⁻¹⁵ under either ridge; against `khaos_exactlaw_oracle.R`, 3·10⁻¹⁵ (raw ridge). Each move's kernel from its statement, in detailed balance and stationary to 4·10⁻¹⁴ (1,800 birth/death, 1,740 degree and 3,132 swap transitions) |
| The chain against the law | 400,000 steps: P(M), four inclusions and the six likeliest structures, worst \|z\| 1.79 (raw ridge) and 0.71 (normalized) |
| Geweke | 12 cases, p 4, M ≤ 3, σ² ~ IG(3, 2), τ held; 200,000 marginal-conditional draws against 300,000 successive steps: worst \|z\| 2.15 and 1.39 (raw ridge, two seeds), 1.31 (normalized) |
| `adaptive_khaos_ridge` readouts | Corrected, on Ishigami rescaled to [0, 1] (below), degree 8, order 3, 6,000 iterations, 3,000 kept: M 6.32 against 6.13, σ² 0.959 against 1.005, the surfaces 0.049 apart (function sd 3.72), total indices within 0.002, no repeated variable |
| Hinge BASS | The whole regression assembly on HEAD `8e42f54` (GBASS slice 5) with this slice, 2,172 of 2,172, the eight pinned hashes and every fingerprint among them. A hinge config's JSON holds neither `Legendre` nor `RawRidge` |

**Parity with `khaos`.** `make_weights` agrees to 3·10⁻¹⁶ and the forward law to 7·10⁻¹⁶ at slice 0's frozen states. The engine's ratios are the corrected ones, so they are not `khaos`'s as written. Each difference is a confirmed slice 0 defect: the forward probability over redraws, the death's normal approximation, C(d, q) for C(d − 1, q − 1), `random_partition`'s law, and the swap onto a held variable.

**Refusals and acceptance.**
- *Refused:* a Legendre basis with categorical predictors, a functional response, d_max below q_max, a hinge degree other than 1, the column-mass support rule, local steps, the fixed-M family, or a raw ridge without its τ. Each is refused with its reason, at validation or at resolution.
- *A Legendre fit's `BassModel.Proposal`* (D01's) is refused, and `Kernel` is KHAOS's.
- *Acceptance,* exact-law fixture (raw ridge): birth 14%, death 17%, degree 40%, swap 8%. On Ishigami: birth 2.3%, death 2.0%, degree 11%, swap 0.1%.

**The cost row** (`KhaosCostTests`, the thread's cycle clock, round-robin in one process, indicative, the machine loaded):
- *A column* at n 1,000, in two runs, the second on a quieter machine: ψ₈·ψ₇ 11.9 and 8.6 µs, ψ₁₅ 10.4 and 7.3 µs, against 66 and 46 µs for the scalar per-point dividing form and 2.6 and 1.4 µs for a two-hinge product: ψ₈·ψ₇ at 4.5× and 6.2× the hinges.
- *A sweep* at the cost report's base with d_max 15 and M 20: 37.7 and 18.9 µs under the normalized ridge, 356 and 190 µs under the raw ridge, 22.3 and 10.2 µs for the hinge fit.
- *A sweep at p 20, d_max 15, q_max 5* (M 2 to 4): 38.2 and 16.3 µs raw, 48.5 and 20.0 µs normalized.

The first column form, four lanes per point-block, a serial chain of dividing steps, read 6.4× a hinge product in an earlier window. It was rewritten step-major and division-free before landing. The raw ridge pays for a fresh pencil decomposition per candidate, the brief's accepted design for a prior apart from the Gram, which slice 2's g-prior shares.

**Found here.**
- *The posterior is not invariant to the input box.* On slice 0's Ishigami fixture, `khaos` on the raw inputs (about [0.008, 0.992]) reads M 10.2 and σ² 0.23 against the noise's 0.25. On the same data rescaled to [0, 1] exactly, which is the box the engine's range scaling reads, `khaos` and the engine alike read M about 6.2 and σ² about 0.95. Eight engine chains of 35,000 steps agree. The likely cause is Ishigami's exact symmetry about ½, which keeps its raw Legendre expansion sparse and which a shift of the centre breaks. That is unverified. The readout gate therefore compares the two on the same rescaled inputs.
- *Two raw-ridge slips in resolution, found by the review and fixed before landing.* The held τ's halving point, which places a ladder's automatic hot end, read the raw ridge as the g-form: τ in place of τ/N, though a Legendre column carries about N observations' information as a normalized one does. The τ prior's log line named the g-form's unit. A test pins the halving point.
- *Slow mixing on a rich surface.* On Ishigami at degree 8 and order 3, `khaos`'s corrected birth and the engine's both accept about 2% of births and almost no swaps. Each chain still matches the other and its own long run, so this is a reading, not a defect.

**Decisions beyond the brief.**
- *The Legendre part as a fourth `BassTerm` part* (above), with the derivative rule and `ReadsVariable`.
- *KHAOS's raw ridge through the design's prior-apart pencil* (D_α = 1, δ = −z̄_w, D_β = I), not through slice 2's prior-precision builder, which does not exist yet. Slice 2's builder subsumes it.
- *q_max is `MaxInteraction`,* a term's continuous factors, and d_max is `LegendreBasis.MaxDegree`.
- *The swap is tallied as `VariableChange`* (the same kind of move), and the degree mutation as a new `DegreeChange`.
- *No column cache* (the brief allowed one by size). A degree-15 column costs 4.5× to 6× a two-hinge product's, by run, about a third of a normalized-ridge sweep (the cost row).
- *The swap's pseudo-count is a setting of its own,* `SwapPseudoCount`, where `khaos` reuses `make_weights`' α (`coin_pars[[2]]`) plus 10⁻⁹.
- *The emulator's term projection* (`BassEmulators`) carries the Legendre part, which no Legendre fit reaches today, since those fits have no categorical predictors.

**Notation crossings** went into `issues/notation-ledger.md` once GBASS committed it: KHAOS's q (a term's order) is BASS's K, BASS's q is the hinge power, and KHAOS's degree d = Σα has no BASS counterpart.

**Left unbuilt, and why.**
- *Slice 2,* waiting for the user.
- *Declared input ranges* (the open question below).
- *A bordered pencil for the raw ridge,* which the cost row does not ask for at these sizes.
- *The adaptive mutation probability,* an option under D10.
- *Moments under `MixtureMeasure` and `EmpiricalMeasure` past degree 15,* slice 4's.
- *A gate for a Legendre basis or the raw ridge under a declared noise (E1).* The config allows both, and the core fits the noise's working response under its weights through the same design and pencil, so nothing in the noise is the basis's; no test runs the pair.

**Recorded as open, with the decision each would change.**
- *Declared input ranges.* Decision 6 reads "inputs on [0, 1] by the fit's declared ranges"; the engine scales by the training range, as BASS does, and the posterior moves with the box (above). The decision: whether a continuous predictor may declare its domain, which `BassInputTransform` would then scale by in place of the training range. If yes, a Legendre fit reads the user's coordinates, as `khaos` does, and its fingerprint takes the ranges. If no, a user who knows the domain rescales before fitting, and `PolynomialChaos.md` says so. The user leans toward a declared domain (2026-11-07).

### Slice 0 addendum: the corrected khaos samplers (chip, 2026-11-07)

**The decision.** For slice 2's readout gate the user chose, on 2026-11-07, to patch khaos's g-prior sampler in R rather than drop the gate ("Patch khaos in R"). The defects shared by both samplers sit in the ridge's too, so the same patches correct both. Slice 1's ridge readouts then meet a corrected reference as well.

**Landed** in `86873ee`: `khaos_corrected(prior)` in `khaos_helpers.R`, the new mode `khaos_exactlaw_oracle.R` with `KhaosExactLawTests`, `corrected=1` in the readouts mode, and two more defect tests. With `BASS_QUALIFICATION=1`, the 36 KHAOS tests pass, and all 36 printed their readings, the two gated runs included. They took about 11 minutes, one class at a time (`xUnit.ParallelizeTestCollections=false`). In parallel on the loaded machine, five R processes at once pushed calls past the bridge's 120 s.

**How the patch is made.** Each patch replaces one expression of the pinned body, matched by `identical()`, and states how many times it must match, so a change in the pinned source stops the oracle rather than leaving a patch unapplied. Both samplers take:
- *The birth's and death's proposal terms:* the exact forward probability Σ_{q₀} J(q₀)·P(χ \| q₀)/P(valid \| q₀) on both sides, by a Poisson-binomial over the coins; the uniform composition with its count C(d − 1, q − 1); and the move probabilities the draw itself gives at M = 0 and at `max_basis`.
- *The degree mutation's count,* C(d − 1, q − 1).
- *The degree draw* indexed through its support (see below).
- *The death's prior ratio,* −ln M.
- *The variable swap* onto inactive variables only.
- *The mutation type* held at ½ (decision 5).

The g-prior also takes:
- ½Δln\|S₀⁻¹\| in all four structure ratios. The swap changes ΨᵀΨ, and the degree mutation changes both g and ΨᵀΨ.
- The current state's Σ_n, Q and ln\|Σ_n\| refreshed at the g₀² the iteration holds.
- g₀² by one slice step on ln g₀² against its derived conditional (below).
- σ² given the structure, then β given σ².

The ridge also takes σ² given the structure, IG(g1 + n/2, d/2), the marginal its ratios read, then β given σ².

**g₀²'s kernel.** The g₀² step is replaced, not only re-targeted. khaos's independence proposal is a Laplace fit to `dgsq_full`, whose tails are lighter than the corrected conditional's. With it re-targeted, the exact-law check below read P(g₀² < 10⁴) as 0.869 against 0.826, with an effective sample size of 196 for ln g₀². That chain is valid but slow in the tail. One slice step on ln g₀², the kernel E19 decided for the engine, reads 0.833 with an effective sample size near 2,500. The step reads the conditional through a pencil: one Cholesky and one eigendecomposition per iteration, then O(p) per g₀².

**The check: the exact structure law** (`khaos_exactlaw_oracle.R`, `KhaosExactLawTests`). The space: four inputs, degree 3, order 2, so 30 multi-indices and the 466 structures of at most two terms, on 40 cases weak enough to spread the posterior.
- *The exact law* integrates λ against its Gamma prior. Under the g-prior it integrates g₀² by Simpson's rule through a pencil, written apart from the chain's; that pencil equals the dense multivariate t to 8·10⁻¹³.
- *The test.* Each reading, P(M), each term's inclusion and under the g-prior P(g₀² < 300, 1000, 10⁴), is set against the chain's batch-means error. E ln g₀² is not used: with no term held, ln g₀² has a long right tail under the vague IG(0.001, 1000), and its mean's batch-means error is unreliable.

| 8,000 iterations | P(M = 0, 1, 2) | Exact | Worst z |
| --- | --- | --- | --- |
| Ridge, corrected | 0.414 / 0.389 / 0.197 | 0.427 / 0.374 / 0.199 | 2.1 over 33 readings |
| G-prior, corrected | 0.377 / 0.357 / 0.266; P(g₀² < 300, 1000, 10⁴) 0.036 / 0.338 / 0.833 against 0.035 / 0.334 / 0.826 | 0.389 / 0.350 / 0.261 | 2.9 over 36 |
| Ridge, as written | 0.047 / 0.449 / 0.503; 5.8% of draws hold a term outside 𝒜 | 0.427 / 0.374 / 0.199 | 49 over 33 |
| G-prior, as written | 0.002 / 0.179 / 0.819; P(g₀² < 300, 1000, 10⁴) 0.015 / 0.217 / 0.925; 16.7% outside 𝒜 | 0.389 / 0.350 / 0.261 | 67 over 36 |

The as-written runs are long, so they sit behind `BASS_QUALIFICATION`.

**The corrected readouts** (Ishigami with one inert input, 120 cases, degree 8, order 3, 6,000 iterations):
- *Ridge:* M 10.2, s2 0.228 against the noise's 0.25, and the surface 0.21 RMSE from the function. Its total indices are 0.571 / 0.429 / 0.248 / 0.000 and its partial indices x₁ 0.317, x₂ 0.422 and x₁:x₃ 0.244, against Ishigami's 0.558 / 0.442 / 0.244 / 0 and 0.314, 0.442 and 0.244.
- *G-prior:* M 13.9, s2 0.203, the surface 0.24 RMSE from the function, and g₀² 3,684 (90% [1,781, 6,965]). Its total indices are 0.577 / 0.423 / 0.251 / 0.001 and its partial indices x₁ 0.321, x₂ 0.417 and x₁:x₃ 0.246. As written, the same data read M 72, g₀² 43 and x₄ 0.153.
- Neither keeps a term that repeats a variable. These, not the samplers as written, are what slices 1 and 2 meet.

**Found here.**
- *The g-prior's current state goes stale when g₀² moves.* Σ_n, Q and ln\|Σ_n\| are updated only when a structure move is accepted, while g₀² moves every iteration. Every structure ratio therefore sets a candidate at the current g₀² against a current state left at the g₀² of its last acceptance; β is drawn from that stale Σ_n too. A probe that draws nothing (the chain is the one written) reads the gap in the current state's log marginal before every move: on the defects' seeded chain (100 cases, 1,500 iterations) the gap is 7.4 nats on average and 32.1 at most, and it passes 1 nat on 84.7% of iterations. On a second data set, generated in R, it reads 39.5, 86.7 and 74%. Confirmed. Together with the missing normalizer, it plausibly accounts for most of the as-written g-prior's M ≈ 72.
- *The degree draw fails where the order is the largest degree.* `sample(q:d_max, 1, prob = d_probs)` at q = d_max hands `sample` a single number, which R reads as 1:q, and it stops with "incorrect number of probabilities". So khaos as written cannot run with order equal to degree. At the defaults, degree 15 and order 5, it never meets this. Confirmed.

**The follow-up is resolved.** Slice 2's g-prior readout gate meets `corrected=1`, and slice 1's ridge gate does too.

### Slice 0: the oracle modes (chip, 2026-11-07)

**Landed** in `31362ad`: four modes under `r/oracles` that run from `r/` against the pinned khaos 2.1.1, sharing `khaos_helpers.R` (the namespace, khaos's expressions evaluated as written, and references that share nothing with its code), with their parity classes in `tests/oracle` and four rows in `r/README.md`'s oracle map. 26 tests (`Maths.Oracle.Tests`, filter `Khaos`), 26 passed, and all 26 printed their readings, so each ran against R; about 106 s together on the loaded machine. None is gated. The C# side reuses `GbassOracle`'s runner and computes its references itself: the collapsed marginal through `PriorPencil`, the Legendre recurrence, |𝒜|, compositions, `random_partition`'s law by counting, and the Poisson-binomial.

| Mode | What it evaluates | Parity class |
| --- | --- | --- |
| `khaos_primitives_oracle.R` | `ss_legendre_poly`, `A_size` and `generate_A`, `make_weights`, `random_partition` | `KhaosPrimitiveParityTests` |
| `khaos_conditionals_oracle.R` | At frozen states (60 cases, four Legendre terms): ln p(y \| structure) under the ridge and the modified g-prior; the samplers' structure ratios as written; g₀²'s conditional integrated and given β; σ² given β; the birth's exact forward law against khaos's own draw, the death's reverse, both mutations, the prior ratios and the move probabilities at M = 0 | `KhaosConditionalParityTests`, and the frozen-state rows of `KhaosDefectParityTests` |
| `khaos_defects_oracle.R` | The other rows: an orthogonal design's per-term penalty, the log-evidence cores and their callers, `random_partition` against both counts, a seeded g-prior chain, the death's normal approximation, the paper's slips, each ratio factor's direction, σ²'s floor | `KhaosDefectParityTests` |
| `khaos_readouts_oracle.R` | `adaptive_khaos` under each prior and `sobol_khaos`: the surface over a grid, M, s2, λ, g₀², the move counts, the indices and the terms that repeat a variable | `KhaosReadoutParityTests`; the engine meets it from slice 1 on |

**The gates.** The modes run from `r/`. Each defect is confirmed or struck below.

**Parity, and what the engine meets** (corrected references throughout):
- *The collapsed marginal.* `PriorPencil` with prior precision D_g⁻¹ΨᵀΨD_g⁻¹ at τ = 1/g₀², and I at 1/τ² for the raw ridge, equals the dense multivariate t over n × n to 10⁻⁹ at five nested structures under each prior. So E21's builder through the pencil already serves both. The ridge's structure ratio as written is that marginal's difference to 10⁻⁹, and its death is the birth's negative.
- *g₀²'s conditional.* With β and σ² integrated, read through the pencil at each g₀², it is the joint less a constant to 6·10⁻¹³ over g₀² ∈ [0.5, 5000]. Given β and σ², IG(a_g + (M + 1)/2, b_g + βᵀD_g⁻¹ΨᵀΨD_g⁻¹β/(2σ²)) holds to 8·10⁻¹³.
- *σ² given β* is IG(a + n/2 + (M + 1)/2, b + RSS/2 + βᵀS₀⁻¹β/2) under either prior, to 10⁻¹⁴.
- *The birth's forward probability* of an inclusion vector is Σ_{q₀} J(q₀)·P(χ \| q₀)/P(valid \| q₀): khaos's own draws carry it (worst z 1.3 and 2.0 over 20,000 births, at order < p and at order = p). The variable swap's ratio is exact for its own kernel, to 10⁻¹².
- *The Legendre polynomials.* `ss_legendre_poly` equals Bonnet's recurrence here bit for bit, since both run the same recurrence. The independent checks: orthonormality to 5·10⁻¹⁵ under the engine's `GaussLegendre` and 1.5·10⁻¹⁴ under R's, and ψ_J vanishing at the J-point nodes to 6·10⁻¹⁴, degrees 8 and 24. R's Golub–Welsch nodes are the engine's to 13 digits.
- *|𝒜|.* `A_size` equals the sum in exact integers and the multi-indices walked here, up to 53,710,887 at (p, d, q) = (20, 15, 5); the sum alone at (30, 15, 5), 467,246,768.
- *`random_partition`.* The law found by enumerating its draws equals the sorted-cuts count to 10⁻¹⁴, and 40,000 of its draws carry it (worst z 3.1).
- *Readouts, ridge* (Ishigami's function on x₁..x₃ and an inert x₄, 120 cases, noise 0.5, degree 8, order 3, 6,000 iterations):
  - The fit: M 10.8, s2 0.232 against the noise's 0.25, the surface 0.17 RMSE from the function, and the 90% bands cover 15 of 16 grid points.
  - The total indices: 0.573, 0.427, 0.248 and 0.000, against Ishigami's 0.558, 0.442, 0.244 and 0. The partial indices: x₁ 0.319, x₂ 0.420 and x₁:x₃ 0.244, against 0.314, 0.442 and 0.244.
  - *A stated difference:* `sobol_khaos` divides by Var f + s2, with the noise taking a 0.016 share. The engine's `BassMoments` divides by Var f, so the mode emits both.
- *Readouts, g-prior as written:* checked for form only (below).

**Refusals and acceptance; the cost row.** Nothing of the engine runs in this slice, so neither applies.
- *R's own acceptance* (birth/death/degree/swap): 0.039/0.033/0.113/0.000 under the ridge; 0.221/0.137/0.132/0.233 under the g-prior.
- *R's costs on the loaded machine.* The ridge readout's 6,000 iterations take 14–43 s and the g-prior's 2,500 take 26–38 s, since it holds M ≈ 72. One `legendre_poly` call at degree 24 costs 0.95 s, because the recursion calls itself on j − 1 and j − 2 with no memory, about φ^j calls. The default degree 15 keeps that cheap; it is no cost the engine's lanes share.

**Defects.** Each confirmed one is a correction the engine takes, not a parity variant.

| Defect | Run | Verdict |
| --- | --- | --- |
| The g-prior's structure ratio leaves out ½Δln\|S₀⁻¹\| | Frozen state: as written equals the marginals' difference less ½Δln\|S₀⁻¹\| to 10⁻¹⁰ per birth, the missing part 0.17–1.28 nats a step. On a 64-case Hadamard design, the code's per-term penalty is ½ln 64 + ½ln(1 + 1/(g₀²g²)), so 2.08–3.94 nats; the model's is ½ln(1 + g₀²g²), 0.012–4.61 nats, over g₀² ∈ [1, 10⁴] and g ∈ {1, 0.30, 0.16}. At g₀² = 1 a q = 3, d = 12 term costs 3.94 nats as written and 0.012 under the model; at 10⁴, 2.08 and 2.76. `log_evidence_fast_core` leaves out the same term (to 10⁻¹⁰). `log_evidence_full_core` also reads the least-squares residual in place of the posterior's quadratic form, 2.5–6.2 nats a step beyond. Nothing in 2.1.1 calls either core: `sparse_khaos` ranks by KIC | Confirmed. Struck: its reach into `sparse_khaos`'s Bayes factors |
| g₀²'s conditional has no data term | `dgsq_full` against the joint with β and σ² integrated: it spreads 80.1 nats as written; 9.21 nats, which is ln 10⁴, with the data term −(a + n/2)ln Q added; and 8·10⁻¹³ with one more factor 1/g₀². Against the joint given β and σ², 2,042 nats | Confirmed, with a second part: the code's power θ^(−a−M/2) (M counting the intercept) leaves out the inverse gamma's −1 |
| The birth's forward probability | Against the exact law above, over every valid χ at p = 4 and p = 3. With no redraw the birth's term is the q₀-marginal P̄(χ), up to 0.33 nats off. Each redraw past q_max multiplies its P̄ in, down to −15.5 nats. An all-zero redraw adds ln 1, since R's `x[-integer(0)]` is empty, so at order = p the term never leaves P̄(χ), 0.29 nats off. The death's reverse is P̄(χ)/(1 − P̄(0) − P̄(q > q_max) by its normal), 0.033 nats off at order < p and 0.148 at order = p | Confirmed |
| `random_partition` not uniform; the ratios' C(d, q) | Uniform for q ≤ 2 (one cut or none). At q ≥ 3 the total variation to uniform is 0.12–0.175 at (5, 3), (6, 3), (7, 4) and (8, 3), with tied cuts weighted down: at (7, 4) the law runs from 1/64 to 6/64 against the uniform 1/20. 1/C(d, q) matches neither law anywhere: it is 1/6 where q = 1's probability is 1, and the degree mutation as written reads 1.204 nats at q = 2, 3 → 5, where uniform (here also the draw's own law) gives 0.693 | Confirmed |
| The swap can pick an active variable | Constructed state, the term {x₁, x₃} at η (2, 1, 1, 0): a repeated variable is proposed with probability 0.40 when x₁ is swapped out and 0.50 when x₃ is. Seeded g-prior chain (100 cases, 1,500 iterations): 55.6% of the kept terms repeat a variable, in all 750 kept draws. In the g-prior readout, 29% of terms, and `sobol_khaos` gives x₂:x₂ 0.388 of the variance. Under the ridge the swap is never accepted | Confirmed |
| The death's normal approximation | Four (η, q_max) cases and nine q₀, against the Poisson-binomial by enumeration on both sides (agreeing to 10⁻¹³): up to 64% off in small tails (0.0018 for 0.0036, 0.0029 for 0.0082) and 7% in the body (0.346 for 0.322) | Confirmed |
| σ²'s floor `s2_lower` | Its default is 0 in both samplers, applied as `max(s2_lower, ·)` | Not ported; recorded |
| The paper's slips | In the code: `legendre_poly(x, 1)` returns x; the coins are `rbinom(p, 1, wts)`; the degree is drawn from {q, …, d_max} ({2, …, 6} at p = 3); the partition count is C(d, q), which is the row above | Struck for the code |
| SC-1's printed ratios | Each factor's exponent in the code's birth, death and degree mutation, as written, is its Metropolis–Hastings exponent; 8 of the 12 are inverted in the print. The code also marginalizes q₀ inside the coins' term, where the print conditions on it | Struck for the code; the manuscript's inversions stand as SC-1 records them |

**Found here.**
- *The death's prior ratio is not the birth's reciprocal.* Birth plus death log prior ratio is ln(K/(K + 1)) to 10⁻¹², at K = 1..4. The birth carries the multiset prior's +ln(M + 1); the death writes −ln(M + 1) where −ln M is meant, so deaths pass K/(K + 1) as often as they should. The corrected pair is λ/|𝒜| and |𝒜|/λ (λ drawn), or its λ-integrated form.
- *The move probabilities at M = 0.* The move is forced to birth there (share 1.000 over 20,000 draws), yet the birth from M = 0 and the death to M = 0 read P_B = P_D = 1/3. Their move parts are 0 as written where the exact parts are −ln 3 and +ln 3. Both errors, with the previous one, push mass away from the empty structure.
- *σ²'s draw given β leaves out β's prior* (the paper's display too). As written it is IG(a + n/2, b + RSS/2); the joint less it spreads 63.6 nats under the g-prior and 20.5 under the ridge over σ² ∈ [0.03, 0.3]. Two more problems sit beside it:
  - The ridge's draw reads its prior rate as g2 where its marginal reads g2/2, which matters only at g2 ≠ 0, and the default is 0.
  - β is drawn given the previous sweep's σ², so (β, σ²) is not a draw from its joint given the structure.

  Decision 2 already corrects all three: σ² given the structure, then β given σ².
- *The degree weights index from q:* (d − q + 1)^(−s_d), where the paper writes d^(−s_d) with default 1. The code's default s_d = 0 is uniform, where the two agree. Either is a valid proposal, since both directions read the same weights. The engine follows the brief's d^(−s_d).
- *`make_weights` sums to q₀ only to its three passes,* up to 3.5·10⁻⁵ off, where the exact root of Σv^β = q₀ meets it to 10⁻¹⁰. Either is a valid proposal; the brief's "summing to q₀" names the root.
- *The g-prior chain as written is no reference.* It holds M ≈ 72 on the readouts' 120 cases. Its indices move with the seed: x₁'s total is 0.69 on one seed and 0.19 on another at 4,000 iterations, and 0.104 on the test's data, where Ishigami's is 0.558, with the inert x₄ at 0.153. Its surface is 1.21 RMSE from the function, and its bands cover 6 of 16 points. These are the defects above acting together.

**Decisions beyond the brief.**
- *The readouts' configuration:* Ishigami with one inert input (p = 4, so the swap runs), 120 cases, degree 8, order 3; the ridge at 6,000 iterations, the g-prior at 2,500 thinned by 2. The defaults fit neither budget:
  - At degree 15, order 5, 200 cases and 5,000 iterations, the ridge never found x₁:x₃ (T₃ 0.004).
  - At 150 cases and 10,000 iterations it did on three seeds (T₃ 0.245–0.255), but took 66–100 s, past the bridge's 120 s on this machine. Under the chosen configuration it does on every seed tried (four).

  A mixing reading of khaos's birth on this function, not a defect.
- *`ROracle`'s timeout* is left at 120 s; the g-prior's chain is sized to it.

**Left unbuilt.** Nothing of the engine: slice 1 waits for the user. The readouts under the g-prior are checked for form only, for the reason above.

**Follow-up, with the decision it changes.** Slice 2's gate reads "`adaptive_khaos_gprior` readouts under the corrected marginal", and the g-prior as written is no such reference. The choice is between two things:
- *Patch khaos's g-prior sampler in the readouts mode,* body deparsed and replaced as GBASS's `nw_gamma_prior` was. The patches: the normalizer in the three ratios, g₀²'s conditional, σ²'s draw, the swap restricted to inactive variables, the death's prior ratio and the move probabilities at M = 0. Slice 2 then has an R readout to meet.
- *Rest the g-prior's gates on the exact law, Geweke and recovery from truth,* with the ridge's readouts as the only R comparison. Slice 2 then drops that gate.

Resolved the same day: the user chose the patch (the addendum above).
