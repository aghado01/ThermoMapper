# KHAOS on the BASS core: E19's build (chip brief)

Date: 2026-11-07. Status: slice 0 landed 2026-11-07 (`31362ad`); slice 1 landed 2026-11-07 (`fa05048`, `b408bf4`, `baf69f9`), with declared input domains (`3eb56e1`, the user's decision of 2026-11-07); the check on slice 1's input-box finding, then slice 2, under way (the user's go-ahead of 2026-11-08). Written at the user's request after GBASS's slices 0–4, from the KHAOS paper's TeX, the `khaos` 2.1.1 source and the engine.

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
