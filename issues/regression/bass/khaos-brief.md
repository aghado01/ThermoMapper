# KHAOS on the BASS core: E19's build (chip brief)

Date: 2026-11-07. Status: written, not started. Written at the user's request after GBASS's slices 0–4, from the KHAOS paper's TeX, the `khaos` 2.1.1 source and the engine.

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

**Recommended; a chip starts on these unless the user says otherwise.**
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

- **Hinge BASS untouched.** The eight pinned hashes and every existing fingerprint hold.
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
- Heavy runs at below-normal priority with capped threads, MSBuild node reuse off; the machine is shared.
- New beside old; a removal is its own question to the user, with what is lost.
- Findings framed by problem class, not by any one application.
- A follow-up goes in a report with the decision it would change and what each outcome would do. A question with no decision waiting on it is recorded as open, not proposed.

## Report

Under "Reports" below, newest first, per slice: what landed by commit; the gates' numbers; the parity comparisons with their stated differences; refusals and acceptance; the cost row; defects confirmed, struck or newly found; decisions taken beyond this brief; what was left unbuilt and why.

## Reports

Newest first.
