# Projection-pursuit ridges on the BASS core: E36's build

Date: 2026-11-10; approved 2026-11-14. Status: **approved by the user on 2026-11-14,** who asked for a chip to build it and gave the go-ahead to install `BayesPPR` (decision 8). The decisions still marked **open** go at their recommendations. That is the coordinating session's reading of the approval, told to the user the same day; the user can turn any of them before the slice where it first binds (listed under Decisions). Drafted by the active-subspace chip from [E36's entry](extensions-brief.md#projection-pursuit-ridges-on-the-bass-core-bppr-e36), the [Grassmann discussion](../discussion/grassmann-stiefel-subspaces-brief.md)'s settlement of the direction prior, the paper (read for this draft) and the read-only `gcbeck/bppr` clone (ideas only, under its GPL-3.0). Reports go under "Reports" below, newest first.

## Why

BASS's terms are axis-aligned. A response that bends along a diagonal, such as one that depends on a few linear combinations of many inputs, needs many interactions to follow it. A ridge g(θᵀx) follows it with one direction.
- *The lineage's third basis,* after hinges and Legendre polynomials: E8 finds a posterior active subspace, rotate-and-refit refits in it (E8's slice 7, a plan with the user), and BPPR learns the directions while it fits.
- *Its natural readout is E8:* C = Σ_{m,l} θ_mθ_lᵀE[g_m′g_l′] lies in span{θ_m}. E8's subspace summaries avoid both label switching across ridges and a varying M, which comparing draws' directions directly would meet.

## Sources

- **Collins, Francom and Rumsey, *Bayesian Projection Pursuit Regression*, arXiv 2210.09181 v1.** Held in `D:\aghado01\codex-scientiae\supellex-private\gauntlet\bayesian\2210.09181v1`, TeX source. §2.3 (the model), §2.4 (the sampler) and §3.1 (Friedman's function) were read for this draft; the rest is for slice 0. What it states:
  - *The model:* f(x) = β₀ + Σ_m β_mᵀb_m(xᵀθ_m | t_m), with M ~ Poisson(λ), λ = 10 suggested.
  - *The ridge:* a modified natural cubic spline with K = 4 basis functions, b₁ = (u − t₀)₊, so the ridge is zero below t₀. The others are divided differences of (u − t_ℓ)₊³, linear beyond the outer knot. t₀ ~ U(L, U), with U the q-quantile of the projections and L = U − (U − Q₀)/p₀, p₀ = 2/3, and q chosen so that at least about 20 projections lie above U. The remaining knots sit at equally spaced quantiles of the projections above t₀.
  - *The directions:* a_m ~ U{1, …, A}, A = min(3, p); J_m uniform among the p-choose-a_m sets; θ_m uniform on the coordinate subsphere of J_m. It does not constrain orthogonality.
  - *The inputs* are standardized to mean 0 and variance 1, "important" for the uniform prior.
  - *The coefficients:* Zellner–Siow, β | B, σ², τ ~ N(0, τσ²(BᵀB)⁻¹), π(σ²) ∝ 1/σ², τ ~ IG(1/2, n/2).
  - *The sampler:* birth, death or change, each with probability one third, then Gibbs steps for β, σ² and τ.
    - Births draw θ and t₀ from their priors and a and J from adaptive proposals: weights ω₀ + counts and Wallenius' noncentral hypergeometric, ω₀ = υ₀ = 1. These read the current state only, so the move is valid.
    - Changes keep a and J, and propose the active coordinates of θ from a power spherical law centred on the current direction, κ = 1000, and t₀ from its prior.
    - It found Bayesian back-fitting slower to mix.
  - *Friedman's function,* the example to meet: n = 300, noise sd 1, 20,000 iterations with the last 2,000 kept. Out-of-sample RMSE 1.05 for the noisy response and 0.38 for f; 95% intervals cover 94.9% and 91.1%.
  - *A slip in the text:* Algorithm 4 writes the change step's bounds with Q and q where eq. (UL) uses q and p₀ (a symbol crossing; the formulas agree).
- **R's `BayesPPR` 0.2.0.** Its manual is in `PDFs/Bayesian/BayesPPR.pdf`. It is not installed here (checked 2026-11-10); pinning it is decision 8, settled on 2026-11-14.
- **`gcbeck/bppr`** (C++20, GPL-3.0), read-only at `D:\aghado01\packages\bayesian-related\bppr-gcbeck` (`1b0d47a`). E36's entry records its reading: a repulsive direction prior det(ΘᵀΘ) whose normalizer c_M is missing as written; reshaped quantile knots; online updates. Ideas only; its claims count as a paper's, with no code to run as an oracle.
- **The power spherical law** (De Cao and Aziz 2020, arXiv 2006.04437, the paper's reference for it) is on disk since 2026-11-11: `D:\aghado01\codex-scientiae\supellex-private\gauntlet\bayesian\2006.04437v2\` (PDF and TeX). Its sampler and density are cited from memory until slice 0 reads it.

## What is already built

- **The design-kernel contract** (F-X3, `IDesignKernel`): value, derivative, admission and the factors a moment service reads. `BassDesignKernelTests` already carries a test-side directional hinge (wᵀx − t)₊ that fits the contract as it stands.
- **The structure chain,** its reversible jump, the collapsed coefficients (Zellner–Siow and the g-form), D06's admission (an order-free angle gate on columns), tempering and SMC with their flags, GBASS's noises, the gradient readout, prediction, the pointwise density and the conformal layer.
- **E8,** whose sampled route reads any fit's C from the gradient readout, under a product law or a joint sample.
- **Basis fits** (D07) and E34's joint fits, which take a term type through the component fits' configuration.

## Decisions

Those marked **settled** are the user's; the rest are **open**, each with a recommendation that stands unless the user turns it. Where each first binds: 8 at slice 0; 1, 3 and 7 at slice 1; 4, 5 and 6 at slice 2. Each lands on a page with its alternatives and trade-offs when built.

1. **The term type, open.**
   - *Recommended, staged:* first the directional hinge (θᵀx − t)₊^q as BASS's fourth term type behind `IDesignKernel`, one column per term, so every carry-over above holds unchanged. Then the paper's ridge as a block of K columns over one direction, each column a fixed combination of directional hinges: degree one at t₀, cubic at the quantile knots.
   - *Alternatives:* the ridge block alone (parity with the paper, a new structural unit in the chain from the start); the hinge alone (BASS-native, with no parity with the paper's numbers).
2. **The direction prior, settled** (the user, 2026-11-08, on the Grassmann discussion's reading): the paper's independent directions, a_m ≤ A uniform, J_m uniform, θ_m uniform on its coordinate subsphere. The determinant det(ΘᵀΘ) stays a labelled option for many inputs and few sparse ridges, with c_M stated and its tilt on M reported. D06's gate applies to the ridge columns, and E21's ridge blocks handle collinearity on the coefficient side.
3. **The metric the direction prior is uniform in, open.**
   - *Recommended:* the paper's standardization, mean 0 and variance 1, recorded in the fit's transform beside the box, so that Friedman's comparison is like for like.
   - *Alternative:* the box, BASS's units and E8's default metric. Under a uniform design the two differ by one common factor, so they give the same prior; under a skewed design they do not.
4. **The knot along a projection, open.**
   - *Recommended:* the paper's continuous t₀ ~ U(L, U) with p₀ and q, for both term types. A knot indexed to a data point, as BASS's are, would jump whenever θ moves the projections.
   - *Alternatives:* gcbeck's reshaped quantile levels (E36's entry), or a discrete knot at a projected case, redrawn on a change.
5. **The moves, open in one respect.** Birth, death and change as the paper's, with the change's power spherical proposal on the active coordinates (κ = 1000) and t₀ from its prior.
   - *Recommended:* J from BASS's own variable proposal, already validated. *Alternative:* the paper's Wallenius form.
   - A prior adapted from the chain's history is not taken (E36's entry); κ adapts only under D10's diminishing adaptation, if at all.
   - *The ratios, from Lean* (`lean/Enthymemata/ActiveSubspace/DirectionMoves.lean`, verified from clean archives at `42c980d` and `36bd90e`; the [Lean chip's brief](../../../lean/Protolemmata/BASS/active-subspace-lean-chip-brief.md) holds the reports):
     - *The change move* accepts on the collapsed likelihood ratio alone (`change_ratio_powerSpherical`), with t₀ drawn from its prior given the proposed direction. At a = 1 it keeps the direction's sign (`powerSpherical_antipode`), so the sampler treats a = 1 apart.
     - *A birth* accepts on L′/L · λ/(M + 1) · P_D(M + 1)/P_B(M) · π(a)π(𝒥|a)/(q(a)q(𝒥|a)) (`birth_green_ratio`, `poisson_count_ratio`, `prior_draw_cancel`); a death on its reciprocal.
     - *No further 1/(M + 1):* with it the chain targets λ^M/(M!)² (`double_count_ratio`).
     - *The paper's slip:* its eq. (mhb) divides by W_{a\*} at a\* = 1, where its Algorithm 2 draws 𝒥 uniformly. The engine codes the corrected ratio, and a parity read against `BayesPPR` is made against the corrected ratio too.
     - *The moves balance as measures* (`mh_flux_symm`, `birth_death_flux`), and inserting a ridge preserves measure (`measurePreserving_piFinSuccAbove`). The code takes prior and proposal densities against one base measure, which slice 2's Geweke gate exercises.
     - *Not on disk:* the paper's supplement, which holds its own change ratio. The derived ratio stands either way.
6. **A basis mixing ridges with BASS's hinge terms, open.**
   - *Recommended:* ridges only first. The contract allows mixing; the moment service has no exact route for a ridge, so a mixed fit's Sobol reads by sampling.
   - *Alternative:* mixed from the start.
7. **Categorical inputs, open.**
   - *Recommended:* refused in the first slices, then BASS's subset factors multiplying a ridge, separable in the category.
   - *Alternative:* the paper's dummy coding with its "in at least one category" indicator ridges.
8. **The oracle, settled** (the user, 2026-11-14): pin `BayesPPR` 0.2.0 from CRAN in `r/renv.lock`, as `coactivity` and BASS's R were pinned. Read its parity on fixed draws: the basis matrix, the collapsed likelihood ratio, the birth's acceptance. Correct its flaws rather than port them.
   - *Not taken:* the paper's numbers and the engine's internal checks alone.
9. **The coefficients, settled by the engine:** Zellner–Siow by default, as the paper's; the g-form as an option; GBASS's noises carried over.
10. **The readouts, settled by the contract:** prediction, the pointwise density, the gradient readout through the derivative, the conformal layer, and E8 on its sampled route. The moment service refuses a ridge with the reason, so Sobol and effects read by sampling (E25).

## The slices

One commit per concern; tests with the code they test; each slice ends with its gates and a fresh-eyes review before the next starts.

| Slice | Delivers | Its own gates |
| --- | --- | --- |
| 0. Sources, references and the power spherical law | The paper read in full; `BayesPPR`'s manual; the oracle if pinned; reference functions with known directions (a single index, sums of ridges, Friedman's function); the power spherical sampler and density | The law's draws against its density, and the cosine's law against its closed form, on chosen κ and a seeded sweep to κ = 10⁵; the proposal's symmetry; each reference's C in closed form |
| 1. The directional hinge | The term type: sparse θ, t, degree q; value and derivative behind `IDesignKernel`; admission; prediction and the gradient readout; the moment service's refusal | The contract test promoted from the test-side kernel; the derivative against central differences; columns against a hand computation; D06 admits ±θ and several t on one θ and refuses a near-duplicate; existing fits bit for bit with the type off |
| 2. Directions in the structure chain | Birth, death and change with θ and t (decisions 4 and 5); tempering and SMC | Geweke's test, prior only, since directions are continuous; the collapsed law given frozen directions against the dense marginal; each move's ratio against a direct evaluation; single-index and sums-of-ridges functions from outside the model, recovered through E8 (principal angles to the truth) |
| 3. The ridge block | The paper's K-column ridge over one direction as one structural unit, with its knots | The basis against eq. (spline_expansion) by hand; Geweke's test; the oracle on fixed draws if pinned; Friedman's function against the paper's 1.05 and 0.38 and its 94.9% and 91.1% (exploration, not validation) |
| 4. The labelled options | The determinant prior with c_M (a Monte Carlo table for sparse directions, the closed form for dense ones) and its tilt on M reported; gcbeck's quantile knots if decision 4 takes them | c_M against p!/((p − M)!p^M) for dense directions; the tilt against its definition; Geweke's test under the option |

## Gates common to every slice

- **Fitting untouched.** With the term type off, every existing fit is bit for bit as before: the pinned tier at each commit, the whole assembly at the slice's end.
- **Truth from outside the model.** Recovery functions are not drawn from the prior, and their directions are known in closed form.
- **Primitives.** Chosen cases and a seeded random sweep, far ranges included.
- **The cost row.** Per sweep and per draw against BASS on the same data: a birth's, a change's and the readout's. Cycle clock, marked indicative.
- **Docs.** Chips leave `docs/` alone; each report lists the doc changes owed.

## Not in this brief

- **Mixed bases and categorical inputs** past decisions 6 and 7's first steps.
- **The multivariate mode.** It comes through basis fits and E34 once the term type is in a configuration. The paper's own mode drops the truncation error, which D07 already reads.
- **gcbeck's online updates.** The engine's route to the same end is SMC along a data-addition path (E36's entry).
- **The Grassmann discussion's multi-index model,** θ_m = Wφ_m with W on Gr(k, p): a different model, worth building only if pooling across ridges matters.
- **Binary classification.** The paper comments it out.

## Working rules

- On `main`, committing by pathspec with nothing left staged and nothing pushed. Commit messages describe the work, with no attribution lines ([CLAUDE.md](../../../CLAUDE.md)). Leave `basin-hold-validation.csv` and `path-profile-validation.csv` out of every commit.
- Other sessions may be writing in the checkout. Before editing a file, check `git diff -- <file>` for changes you did not make; on a shared file, message its other writer first and commit only your own hunk. Build from a snapshot (`git archive HEAD` plus your files) and commit only compiling states.
- Builds set `ArtifactsRoot` inline into `artifacts\chip-bppr\`, deleted at the end. Heavy runs at below-normal priority with capped threads, MSBuild node reuse off; the machine is shared. Tests run by [targeted runs](../../../CONTRIBUTING.md#testing-standards).
- **Sources:** read a paper before stating what it says, and mark anything cited from memory. When a source is not on disk, ask the main session; fetch nothing yourself. The one fetch allowed is decision 8's: `BayesPPR` 0.2.0 and what it depends on, from CRAN through `r/`'s renv. `gcbeck/bppr` is read for ideas only.
- New beside old; a removal is its own question to the user, with what is lost.
- No notes outside the repository: hand-over state goes in this brief's report as an "In progress" entry, probes and their data go in `artifacts/` (AGENTS.md, 2026-11-10).
- A commit that touches shared code names in its report every test project that `scripts/test-select.ps1` selects for it, read from disk, with each one's count.
- After each slice, report to the main session and wait for its verification note before starting the next. For this chip "the main session" is the second coordinating session (the [handoff](../coordination-handoff-2026-11-11.md#sessions)'s Sessions table), which verifies each slice from a clean archive. `artifacts/verify-main` is the coordinators' workspace; leave it alone.

## Report

Under "Reports" below, newest first, per slice: what landed by commit; the gates' numbers; decisions taken beyond this brief; what was left unbuilt and why; follow-ups with the decision each would change; the doc changes owed.

## Reports

### Slice 0: sources, references and the power spherical law (2026-11-14)

**Landed.**
- `a7a7667` pins CRAN's `BayesPPR` 0.2.0 in `r/renv.lock` (DESCRIPTION import, README line). It imports only base R's `methods` and `parallel`, so nothing else entered the lock. CRAN binary, published 2026-09-18, maintainer J. Derek Tucker. renv's description hash `d20a627318d483a1d721a89c3b2376fb`; the lock is sanitized of hashes. A short fit of exp(0.6x₁ + 0.8x₂) on [0, 1]³ settles on one ridge at (0.596, 0.803, 0.003).
- `63b5a63` fixes `SpecialFunctions.Log1p`: Goldberg's ln(u)·x/(u − 1) overflowed to +∞ above x ≈ 2.5·10³⁰⁵. It now returns ln u there. Only outputs that were +∞ change, so no finite result moves.
- `8be8b1d` adds `SpecialFunctions.LogGammaRatio(x, s)` = ln Γ(x + s) − ln Γ(x). It raises x to 16 by the recurrence, then differences Stirling's series term by term. The power spherical normalizer at κ = 10⁵ needs it: there each ln Γ is about 10⁶, and their difference would lose about five digits.
- `738630b` adds `PowerSpherical` (`Maths.Distributions`): density, log normalizer, `CosineLaw` (Beta(α, β) of (1 + μᵀx)/2) and draws, on S^(d−1) for d ≥ 1.
  - A draw is two gamma variates, so 1 − t and √(1 − t²) keep their digits as t → 1.
  - The map onto μ is the reflection along e₁ − μ, or minus the one along e₁ + μ when μ₁ > 0, so the map is defined at μ = e₁ and never subtracts nearly equal numbers.
  - The density reads 1 + μᵀx as ‖x + μ‖²/2. It is exact near the antipode and symmetric in μ and x bit for bit.
  - On S⁰ it is the point mass at μ for κ > 0 (Lean's `powerSpherical_antipode`), and uniform on ±1 at κ = 0.
- `741d5ec` adds E36's references (`tests/regression/ProjectionPursuitReferences.cs`) beside E8's `ReferenceFunction` family. Each carries its law and the unit directions that span C's range:
  - a single index exp(0.6x₁ + 0.8x₂) on [0, 1]⁴;
  - a dense single index sin(π(x₁ − x₂ + ½x₃ + ½x₄)) on [0, 1]⁵, with four active coordinates, more than one ridge's three;
  - a single index cos(2(x₁ + x₃)) on [0, 1] × [0, 4] × [−1, 1], where standardizing turns (1, 0, 1)/√2 into (1, 0, 2)/√5, for decision 3's metric;
  - two ridges sharing x₂, at cosine 0.36;
  - Friedman's function on [0, 1]⁶ (`FriedmanReference`), whose C is in closed form through Cin(2π). ∂₄f : ∂₅f is 2 : 1 everywhere, so 10x₄ + 5x₅ is one ridge on (2, 1)/√5 and C has rank 4, spanned by e₁, e₂, e₃ and (2e₄ + e₅)/√5. The paper's posterior agrees: its ridges sit on {3}, {4} and three times on {1, 2, 5}.

**Gates.**
- *The law's draws against its density:* 12 cases, d ∈ {2, 3} by κ ∈ {0, ½, 3, 30, 10³, 10⁵}, 20,000 draws each. Cells sit off-centre from μ (arcs on the circle, polar bands by azimuthal sectors on the sphere), and each cell's probability is quadrature of `LogPdf` itself. Smallest χ² p-value 0.083. The cells sum to 1 within 1.7·10⁻¹¹, except 4.1·10⁻⁸ at d = 3, κ = ½, where the antipode is a cone point of the density inside a panel.
- *The cosine's law against its closed form, and the tangent's against the uniform's:* Kolmogorov–Smirnov on (1 − μᵀx)/2 ~ Beta(β, α) and on the tangent along a fixed direction orthogonal to μ, over all draws and within each half split at the cosine's median.
  - 11 chosen cases, d up to 16 and κ up to 10⁵: worst √n·D 1.43.
  - A seeded sweep of 300 cases, d ∈ [2, 16], κ ∈ {0, 1000} ∪ [10⁻³, 10⁵]. A quarter of the means are at the reflection's branch (μ₁ = ±10⁻¹⁷…10⁻¹³) and a quarter next to ±e₁ (1 − |μ₁| from 10⁻¹⁴ to 10⁻⁶). Worst √n·D 1.81, worst sign |z| 1.83, worst |‖x‖² − 1| 2.7·10⁻¹⁵. The cut is 2.6.
- *The normalizer:* against the beta's, through the tangent-normal decomposition, for d = 2…40, nine κ to 10⁵ and seven cosines down to w = 10⁻¹². Worst error is 0.13 of a rounding tolerance scaled by the ln Γ terms. At κ = 0 the density equals 1/|S^(d−1)|.
- *The proposal's symmetry:*
  - The density is symmetric in μ and x bit for bit, over 4,000 pairs with d from 1 to 12.
  - One step from the uniform law: 5 cases, 40,000 pairs each. The law of θ* stays uniform, with √n·D at most 1.26. Three statistics that change sign under swapping the pair, h(θ) − h(θ*) for two choices of h and θ₁θ*₂ − θ*₁θ₂, are positive about half the time: sign |z| at most 1.83.
  - On S⁰, 1,000 draws all keep the sign.
- *`LogGammaRatio`:*
  - integer steps against Σ ln(x + i), for x from 10⁻⁶ to 10⁹;
  - s = ½ against the central binomial product to 10⁻¹⁴, for n up to 10⁶;
  - 2,000 moderate cases against the ln Γ difference;
  - a 20,000-case additivity sweep over x ∈ [10⁻³⁰⁰, 10¹⁵] and s, r ∈ [10⁻⁶, 10⁶]. This sweep found the `Log1p` overflow.
- *Each reference's C in closed form:*
  - against tensor Gauss–Legendre quadrature of its gradient: at most 1.8·10⁻¹⁵ of C's largest entry;
  - gradients against fourth-order differences: at most 7.9·10⁻¹²;
  - the declared directions against C: ‖C − PCP‖/‖C‖ at most 4.6·10⁻¹⁶, and C is positive definite on the span (smallest pivot over largest diagonal 0.10 for the two ridges, 0.55 for Friedman's), so the span is C's range, rank included;
  - Cin(2π) = 2.437653393057224, summed in decimal, against compensated quadrature within 1.5·10⁻¹⁵.
- *Fitting untouched:* every test project that test-select selects for these commits, run on a snapshot of `HEAD` holding the four commits' files (`artifacts/chip-bppr/snap`).
  - *Selected, per commit (classes, read from disk):*
    - `63b5a63` (`Log1p`) and `8be8b1d` (`LogGammaRatio`): Distributions 19 of 20, Regression 324 of 329 (325 of 330 for `8be8b1d`), Samplers 26 of 31, Oracle 36 of 39, Synthetic 6 of 8, DimReduction 8 of 8, Mapper 1 of 10, Ph 1 of 43, Geometry 1 of 9, VizCore 50 of 85, CumulativeField 3 of 3.
    - `738630b`: Distributions 1 of 20, nothing else. `741d5ec`: Regression 1 of 330, nothing else.
    - `Clustering.Graphical.SPC.Tests` holds the changes in its closure and selects 0 of 1; it has no pinned test. HdbScanSmoke, the Spc projects and VizCoreSmoke have no test classes. No other project holds the changes in its closure, so none can move.
  - *Whole assemblies:* Regression 3497 of 3497, Distributions 429 of 429 (and 416 of 416 at commit scope before it).
  - *Commit scope (selected classes, slow tier left out, pinned tier in):* Samplers 259 of 259, Synthetic 29 of 29, DimReduction 83 of 83, Mapper 1 of 1, Ph 1 of 1, Geometry 4 of 4, VizCore 384 of 386 (2 skipped by their own `Skip`), CumulativeField 15 of 15, SPC no test.
  - *Oracle at commit scope:* 204 of 207. The three failures were R oracles timing out at 120 s under load: `transport_oracle.R` (`DiagramDistances_MatchT4TransportAndLpAssign`, p = 1 and 2) and `gbass_readouts_oracle.R` (`HBass_MeetsCorrectedHbass_…`). Rerun alone, all three pass, in 115 s, 114 s and 79 s. The transport pair sits 5 s inside its limit even alone.
- *The cost row:* nothing to cost yet. Slice 0 has no move or readout; the first row comes with slice 1's readout and slice 2's moves.

**The paper, read in full (with its supplement) and De Cao and Aziz.**
- *The supplement is on disk,* as `supplemental_materials.tex` beside `main.tex`; the Sources section above says it is not. It holds the death and change ratios. Its change ratio is the collapsed likelihood ratio alone, as `change_ratio_powerSpherical` derives.
- *The paper's slips beyond the two the brief records* (eq. (mhb) at a\* = 1, and Algorithm 4's bounds):
  - Algorithm 1 draws β and then σ² ~ IG(n/2, ‖y − Bβ‖²/2). That is not σ²'s full conditional, since β's prior scales with σ²; the full conditional is IG((n + 1 + KM)/2, (‖y − Bβ‖² + βᵀBᵀBβ/τ)/2). `BayesPPR` draws σ² from its marginal given τ, IG(n/2, (yᵀy − τ/(1 + τ)·qf)/2), before β, which is correct.
  - The supplement's death ratio reads υ_j at s − 1, with the dying ridge still counted. The reverse birth proposes from the reduced state, so the reduced counts are right, as they already are for ω_a. `BayesPPR` uses the reduced counts for both.
- *Conventions to note for parity:*
  - The intercept sits inside the g-prior, shrunk toward 0 by τ/(1 + τ). This is R BASS's zero centre and is not location-invariant. A parity read takes the engine's zero-centre option.
  - Bounds and knots are R's type-7 quantiles.
  - q = 1 − min(20, n/10)/n.
- *Not on disk:* the published version (Statistics and Computing 2023, DOI 10.1007/s11222-023-10334-z, as the package's DESCRIPTION cites it). v1 was read.
- *De Cao and Aziz:* the law, its Beta marginal and Algorithm 1 agree with the code above. Their Householder vector e₁ − μ is undefined at μ = e₁; the engine's sign choice removes that.

**`BayesPPR` 0.2.0 against the paper and the Lean ratios.** Read from the installed namespace, deparsed to `artifacts/chip-bppr/oracle-source/`.
- *Where it meets Lean:*
  - A birth accepts on L′/L · λ/(M + 1) · P_D(M + 1)/P_B(M) · π(a)π(𝒥|a)/(q(a)q(𝒥|a)), with no further 1/(M + 1). Births append and deaths pick uniformly; that is the point-process accounting, where the paper's random insertion position cancels the same factor.
  - P_D/P_B follows the move menu: birth alone at M = 0; death or change at the cap; birth or death when no ridge has a direction. The paper's uniform choice of move omits this factor, which is right only if impossible moves are proposed and rejected.
  - At a = 1, 𝒥 is uniform and no Wallenius term enters, which corrects eq. (mhb).
  - A death's reverse birth uses the reduced counts, which corrects the supplement.
  - A change accepts on the likelihood ratio alone.
  - BASS's own variable proposal is the same law: uniform for one variable, Wallenius on counts plus one for several. So decision 5's recommendation coincides with the corrected paper.
- *Departures:*
  1. **At a = 1 a change draws the sign by a fair coin,** θ* = ±1 independent of θ, with a fresh t₀. This is valid: a symmetric independence proposal under the uniform prior on ±1, accepted on the likelihood ratio. The paper's law at a = 1 would keep θ. The brief says only that the sampler treats a = 1 apart; the form is a choice for slice 2, below.
  2. M is capped at min(150, ⌊n/K⌋ − 2), a truncated Poisson.
  3. Inputs with two distinct values become categorical and are not standardized. Inputs with three to K values become "disc", with a linear ridge (one projection column). Neither is in the paper.
  4. Admission refuses a design whose Cholesky diagonal ratio, max(d₂, …)/min(d), exceeds 1000. That is R BASS's `getQf` rule, order-dependent, and D06 replaces it.
  5. `adapt_act_feat = FALSE` drops π/q for a and 𝒥 from the ratio, so `w_n_act` and `w_feat` become the prior. With TRUE they are the proposal's starting counts under a uniform prior. A sampler flag thus changes the target; the manual documents it.
  6. **`bppr_resume` restarts the move-menu factor at 0** whatever M is resumed. Births and deaths after a resume are off by ln 3 or ln 2 until one is accepted. This is a defect in resume only.
  7. `dwallenius` is R BASS's inclusion–exclusion (`dmwnchBass`), which loses digits under skewed weights.
  8. `rps` divides by zero at μ = e₁, a measure-zero event that ends in a rejection. Through t = 2z − 1 it loses about log₁₀κ digits of √(1 − t²); this is immaterial at κ = 1000.
  9. `prior_coefs = "flat"` uses c = 1 with a ln 10⁻⁵ penalty per ridge, which is no coherent posterior over M. It is not ported.
  10. `n_adapt` iterations skip the σ², β and τ steps. That affects burn-in only; the default is 0.
  11. `scale_proj_dir_prop = s` maps to κ = (1/s − 1) + √((1/s)(1/s − 1)), so s ≈ 0.002 gives κ = 1000. This map is not in the paper.
  12. The manual's argument table is shifted one row (for example, `n_ridge_max` is described as the Poisson mean).
- *The recommendations, read against the sources:* nothing argues against decisions 1, 3, 4, 5, 6 or 7.
  - On 3: `BayesPPR` standardizes with the sample sd (n − 1). The engine's transform should match it for a like-for-like Friedman comparison.
  - On 7: since `BayesPPR` makes two-valued inputs categorical by itself, a parity read must use inputs with more than K distinct values.

**Decisions taken beyond the brief.**
- The law is a class in `Maths.Distributions`, defined for d ≥ 1, S⁰ included.
- `LogGammaRatio` lives in `SpecialFunctions`. Because `SpecialFunctions` is a hub, test-select selects most classes in ten projects for these commits.
- The `Log1p` fix was made in place, since it moves only +∞ results.
- The references are test fixtures beside E8's.

**Unbuilt.**
- The oracle script, `r/oracles/bppr_oracle.R`: the fixed-draw parity read (the basis matrix, the collapsed likelihood ratio, the birth's acceptance) binds at slice 3, where the brief's table puts the oracle.
- No term type yet; that is slice 1.

**Follow-ups.**
- *The change move at a = 1* (binds at slice 2):
  - (i) `BayesPPR`'s fair coin on the sign, with t₀ from its prior: a one-coordinate ridge can flip ends in one step.
  - (ii) Keep the sign and redraw t₀ only: a sign then changes only through death and rebirth.
  - Recommended: (i). Both accept on the likelihood ratio.
- *`LogBeta` at large shapes:* `LogBeta` has the cancellation `LogGammaRatio` removes, about 10⁻¹⁰ absolute at shapes near 10⁵. Moving `LogBeta` onto `LogGammaRatio` would change the bits of fits that read beta densities (BARS's locality proposal), so it is a question of its own. It is not done.
- *The Sources section above* still says the supplement is not on disk. Its text is left as approved.
- *Not E36's, found by its gates:*
  - `GeneralizedBetaPrimeTests.TheIntegratedDensity_MeetsTheGeneralizedHyperbolicsClosedForm` at (−0.5, 0.09, 1.0, 1.2, 0.5) ran 22 min at commit scope and 15 min in the whole assembly; its comment speaks of 20 s. `Maths.Distributions.Tests` has no `slow-tests.txt`, so every Distributions scope runs it. Every `Log1p` call on its path takes |x| ≤ 1, where `63b5a63` is bit for bit, so these commits do not move it.
  - The test host of a test-select run is not below normal. Read during the Regression run, the chain was testhost (Normal) ← dotnet (Normal) ← dotnet (Normal) ← test-select's dotnet (BelowNormal). So test-select runs do not keep the shared-machine rule.
  - `transport_oracle.R`'s cases take 115 s alone against the oracle's 120 s timeout, so any load fails them.

**Doc changes owed.**
- `docs/Maths/ProbabilityAndInformation.md`: the power spherical law (density, cosine law, draw, symmetry, S⁰), `LogGammaRatio`, and `Log1p`'s range.
- `docs/Maths/Regression/Sources&Changes.md`: the paper's σ² step and the supplement's death weights, `BayesPPR`'s resume defect, the `adapt_act_feat` target switch, and the a = 1 coin.
- `docs/Maths/Regression/ActiveSubspaces.md`: E36's references beside E8's, and Friedman's rank-4 C.
