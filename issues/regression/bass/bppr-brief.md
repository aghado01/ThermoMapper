# Projection-pursuit ridges on the BASS core: E36's build (draft brief)

Date: 2026-11-10. Status: draft for the user's review. Nothing is built until the user approves it, with the decisions marked **open** below. Drafted by the active-subspace chip from [E36's entry](extensions-brief.md#projection-pursuit-ridges-on-the-bass-core-bppr-e36), the [Grassmann discussion](../discussion/grassmann-stiefel-subspaces-brief.md)'s settlement of the direction prior, the paper (read for this draft) and the read-only `gcbeck/bppr` clone (ideas only, under its GPL-3.0). Reports go under "Reports" below, newest first.

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
- **R's `BayesPPR` 0.2.0.** Its manual is in `PDFs/Bayesian/BayesPPR.pdf`. It is not installed here (checked 2026-11-10); pinning it is decision 8.
- **`gcbeck/bppr`** (C++20, GPL-3.0), read-only at `D:\aghado01\packages\bayesian-related\bppr-gcbeck` (`1b0d47a`). E36's entry records its reading: a repulsive direction prior det(ΘᵀΘ) whose normalizer c_M is missing as written; reshaped quantile knots; online updates. Ideas only; its claims count as a paper's, with no code to run as an oracle.
- **The power spherical law** (De Cao and Aziz 2020, arXiv 2006.04437, the paper's reference for it) is **not held**. Its sampler and density are cited from memory until slice 0 reads it; the main session asks the user if it is wanted on disk.

## What is already built

- **The design-kernel contract** (F-X3, `IDesignKernel`): value, derivative, admission and the factors a moment service reads. `BassDesignKernelTests` already carries a test-side directional hinge (wᵀx − t)₊ that fits the contract as it stands.
- **The structure chain,** its reversible jump, the collapsed coefficients (Zellner–Siow and the g-form), D06's admission (an order-free angle gate on columns), tempering and SMC with their flags, GBASS's noises, the gradient readout, prediction, the pointwise density and the conformal layer.
- **E8,** whose sampled route reads any fit's C from the gradient readout, under a product law or a joint sample.
- **Basis fits** (D07) and E34's joint fits, which take a term type through the component fits' configuration.

## Decisions

Those marked **settled** are the user's already; the rest are **open**, each with a recommendation. Each lands on a page with its alternatives and trade-offs when built.

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
6. **A basis mixing ridges with BASS's hinge terms, open.**
   - *Recommended:* ridges only first. The contract allows mixing; the moment service has no exact route for a ridge, so a mixed fit's Sobol reads by sampling.
   - *Alternative:* mixed from the start.
7. **Categorical inputs, open.**
   - *Recommended:* refused in the first slices, then BASS's subset factors multiplying a ridge, separable in the category.
   - *Alternative:* the paper's dummy coding with its "in at least one category" indicator ridges.
8. **The oracle, open.**
   - *Recommended:* pin `BayesPPR` 0.2.0 from CRAN, the user's go-ahead for the install, as `coactivity` and BASS's R were pinned. Read its parity on fixed draws: the basis matrix, the collapsed likelihood ratio, the birth's acceptance. Correct its flaws rather than port them.
   - *Alternative:* the paper's numbers and the engine's internal checks alone.
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
- **Sources:** read a paper before stating what it says, and mark anything cited from memory. When a source is not on disk, ask the main session; fetch nothing yourself. `gcbeck/bppr` is read for ideas only.
- New beside old; a removal is its own question to the user, with what is lost.
- No notes outside the repository: hand-over state goes in this brief's report as an "In progress" entry, probes and their data go in `artifacts/` (AGENTS.md, 2026-11-10).
- After each slice, report to the main session and wait for its verification note before starting the next.

## Report

Under "Reports" below, newest first, per slice: what landed by commit; the gates' numbers; decisions taken beyond this brief; what was left unbuilt and why; follow-ups with the decision each would change; the doc changes owed.

## Reports
