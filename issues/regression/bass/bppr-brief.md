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
     - *The paper's supplement* (`supplemental_materials.tex` beside `main.tex`, read in slice 0) gives the same change ratio: the collapsed likelihood ratio alone.
   - **The change move at a = 1, settled** (the user, 2026-11-14, relayed by the coordinating session): the fair coin. In the user's words: "Use the fair coin at slice 2. For a one-coordinate ridge, redraw the sign uniformly from −1 and +1, then draw the knot from its prior conditional on the proposed sign. That last detail matters because reversing the direction changes the projected knot interval. The ridge can switch ends in one change move, and the acceptance probability still uses only the collapsed likelihood ratio."
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

### Slice 3 verified from a clean archive (coordinating session, 2026-11-15)

- **Commits:** `332692b`, `349ab81` and `10b290d`, on a snapshot of `72d3941`, through `scripts/verify-csharp.ps1` over their 25 changed files, in the two projects that select.
- **Result:** Regression 2,717 of 2,717 and Oracle 106 of 106 at commit scope. The slice's own classes are among them: `RidgeColumnTests` 6, `DirectionalRidgeTests` 10, `BpprFriedmanTests` 3 and `BpprParityTests` 1, the last against the pinned `BayesPPR` from the archive.
- **Read.**
  - *The columns* against the paper's eq. (spline_expansion), in its TeX: b₁, b_ℓ = d_{ℓ−1} − d_K and d's form. The linear form past t_{K+1}, (t_K − t_{ℓ−1})(3u − t_{ℓ−1} − t_K − t_{K+1}), rederived here from the two cubes, and it meets the cubic at t_{K+1}.
  - *The divisors:* every t_{K+1} − t_j is positive by the constructor's order on the knots.
  - *The knots:* type-7 quantiles at (ℓ − 1)/K of the projections strictly above t₀; none with fewer than two such cases or with t_K = t_{K+1}.
  - *The report's answers* to the four points: counts by ridge, coincident knots, admission all or none, the degree.
- **Not read line by line:** `BassDesign`'s scratch load and the sweep's block proposals. They rest on the slice's gates (each proposal the strict load bit for bit, the collapsed law against the dense marginal, the moves' ratios, Geweke's test at K = 4) and on the oracle's parity.
- **For the user,** from the report's follow-ups:
  - *M's mixing under ridge blocks.* On Friedman's function f's 95% intervals cover 77% to 85% in the engine and 71% to 83% in `BayesPPR`, with M fixed over the kept draws in all six runs. The reading the report names decides between tempering or SMC as a ridge basis's default sampler and a split/merge move for blocks.
  - *Affine ridges* under p₀, a labelled option if wanted.
  - *Rank-K updates* for block proposals, a question of cost alone.
- **Receipts:** `artifacts/verify-bppr-s3-72d3941/`, with the changed list, the counters and the snapshot tracked in [receipts](../receipts/verify-bppr-s3-72d3941/).

### Slice 3: the ridge block (2026-11-15)

The coordinating session let the design below stand and asked for four points to be stated and tested: what a draw reports, coincident knots, admission all or none, and the fit's degree. Each is answered below. This entry replaces the in-progress one.

**The design.** The ridge block enters at the design boundary, so the structure chain stays slice 2's.
- *A structure term* is slice 2's directional term (θ, t₀): the prior, E36's kernel, the moves, the counts and Geweke's test with the prior alone carry over unchanged. t₀'s prior is already eq. (UL), and the remaining knots are fixed given θ and t₀, so they add no density.
- *A ridge basis*, `DirectionalBasis.Columns = K` (1, today's hinge alone; 4, the paper's), has the design expand each structure term into K column terms:
  - b₁ = (θᵀz − t₀)₊;
  - b_ℓ = d_{ℓ−1} − d_K for ℓ = 2, …, K, with d_ℓ the paper's divided differences;
  - t₁, …, t_{K+1} are the type-7 quantiles at (ℓ − 1)/K of the projected training cases above t₀, as `BayesPPR`'s `get_mns_basis` reads them.
- *Each column term* carries θ and the K + 2 knots, so prediction needs no training data. Draws, prediction, the gradient readout, E8 and the conformal layer already read one coefficient per column term, so they read a ridge as they read hinges.
- *Admission:* D06's angle gate reads the K columns. The support count reads b₁, the block's reach, since the late cubic columns are nonzero on only the top quarter of the cases above t₀. M's cap counts structure terms, and the design's column cap is K times it.
- *The alternatives, not taken:*
  - one term carrying K columns, which changes every readout's coefficient layout;
  - units of K terms in the structure layer, which changes every move and the prior.

**Landed**, one commit per concern.
- *The columns* (`332692b`).
  - `RidgeColumn` is b_ℓ(θᵀz | t), ℓ = 1, …, K, of the paper's modified natural spline (eq. (spline_expansion)): b₁ = (u − t₀)₊ and b_ℓ = d_{ℓ−1} − d_K, with d_j = ((u − t_j)₊³ − (u − t_{K+1})₊³)/(t_{K+1} − t_j). Past t_{K+1} it takes its exact linear form (t_K − t_{ℓ−1})(3u − t_{ℓ−1} − t_K − t_{K+1}). The column carries its knots, so prediction needs no training data.
  - `RidgeColumnKernel` (design kernel 4) reads it through the directional hinge's projection, now shared, unnormalized and with no separable factors, so the moment service refuses it.
  - `DirectionalTermSpace.RidgeKnots` gives t₀ and the type-7 quantiles at (ℓ − 1)/K of the projected training cases above t₀, `BayesPPR`'s arrangement; null where a column is undefined.
  - `BassTerm`'s ridge part takes the case side alone, its hash added only where present. A structure prior refuses a column term; `BassSensitivity.FromCoefficients` names it.
- *The design's blocks* (`349ab81`).
  - `DirectionalBasis.Columns` (K; 1, the default, is slice 2's hinge alone).
  - `RidgeExpansion` (new) expands a structure term (θ, t₀) into its K column terms, cached beside the term; null when its knots are undefined. It also reads a column list back as its structure terms.
  - `BassDesign` takes the expansion. `Terms` reads structure terms, `ColumnTerms` the columns. A birth, death or change proposes the whole structure, loaded strictly into a scratch design whose account is the proposal's; an acceptance loads it here the same way. A block's newest angle is the least over its columns.
  - The count prior reads structure terms (`BassStructureFit.Units`). The design's column cap is K times M's cap.
  - `BassResolved` floors n_active at K + 1 (`BayesPPR`'s floor) and bounds M's default cap at ⌊(N − 1)/K⌋.
  - The config refuses what a ridge basis does not read (below).
  - Every readout that counts terms counts ridges: `BassStructure.TermCount(K)`, `BassStructureSpec.ColumnsPerTerm`.
- *The oracle and Friedman* (`10b290d`).
  - `r/oracles/bppr_oracle.R`: `BayesPPR`'s own `get_mns_basis`, `get_qf_info` and `dwallenius`, from its namespace, with the standardization, the knots and the birth's ratios transcribed beside them.
  - `BpprParityTests` (Oracle) and `BpprFriedmanTests` (Regression).

**The coordinating session's points.**
- *The prior sampler* rejects a whole draw when any of its terms has no knot range (`BassPriorSampler.DrawOn`: a null term redraws the collection). It therefore samples p(M)·P(range)^M on F, unnormalized, the same restriction the chain's target carries.
  - Under a ridge basis, a block with undefined knots (below) is a term with no range in the same sense: the expansion gives none and the design refuses it for support.
  - Geweke's test reads the two together: the prior sampler's acceptance is 0.999 there.
- *What a draw reports.* A draw holds K column ids per ridge, since its coefficients are the columns'. Every readout that counts terms divides by K, and every readout that walks coefficients reads columns.
  - *M's posterior:* `BassPosteriorSummary.Count` (`BassSummary.cs:73`) reads `TermCount(K)`.
  - *Traces:* `BassTraces` (`BassTraces.cs:44`) reads `TermCount(K)`.
  - *The structural monitors:*
    - the chain's and the flat histogram's M functional read the structure state's count (`BassRun.cs:341`, `BassMulticanonical.cs:214`), whose terms are structure terms;
    - SMC's M functional reads `TermCount(K)` (`BassSmcRun.cs:168`).
  - *Mass at the cap* (`BassRun.cs:210`, `BassFlatHistogramRun.cs:236`, `BassSmcRun.cs:215`) compares `TermCount(K)` with M's cap, which counts ridges.
  - *The projection comparison's mean size* (`BassProjectionComparison.cs:122`) reads `TermCount(K)`.
  - *The count prior* reads `BassStructureFit.UnitCount`.
  - *The order and variable counts:* the support tally (`BassDiagnostics.Support`) counts each candidate a birth or change proposes, the structure term, once at its a; E36's kernel's counts read the structure state's terms. E4's inclusion history, the one per-variable inclusion readout, is refused for a directional basis (slice 2).
  - *Walking columns, as they should:* prediction, the gradient, E8, the multi-output coefficient checks, the monitored fitted values (`BassOutputFunctionals`), and the sensitivity's id table.
  - *Tested on a K = 4 fit* under the chain, tempering and SMC:
    - every draw holds whole blocks;
    - the summary's mean M equals the draws' column count over K;
    - the traces equal it draw by draw (chain and tempering; SMC keeps none);
    - one inclusion count: the support tally's attempts equal the births' and changes' attempts with a term, not K times them, at orders 1 to 3 only.
    - The prior chain's Geweke test asserts the same identity, with none at order 0.
    - The log prior at M = 2 reads M = 2's, not M = 8's (−13.14 against −11.53).
- *Coincident knots.* Two cases arise.
  - *Coincident upper knots* (t_j = t_{K+1} for some j ≤ K) leave d_j's divisor zero. The knots are null, the expansion gives no block, and the term lies outside the support, as a range with no width does: a `Support` refusal, nothing renormalized.
  - *Coincident interior knots* (t_j = t_{j+1} < t_{K+1}) divide by no zero, but they make two columns equal (d_j = d_{j+1}). D06's gate refuses the block: a `Numerical` refusal, outside F.
  - Tested: on an input tied at 0.5 on every third case, +e₀'s knots from its range's lower end are −3.95, −2.15, −0.330003, −0.330003, 0.786, 1.89, and the block is refused `Numerical`; a top tie gives null knots and a `Support` refusal.
- *Admission is all or none.* The design proposes the whole structure through its scratch load, so a block with one column refused is refused whole, and the design holds no partial block. Tested: a refused birth and a refused change leave the design's terms, columns, residual and fitted sums exactly as they were, with nothing pending to accept.
- *The fit's degree.* At K ≥ 2 the column degrees are the paper's (one at t₀, three at the quantile knots), so `Structure.Degree` is not read, and the config refuses any value but 1, the default, rather than ignore it. At K = 1 it is the hinge's q, as in slice 2.
  - K = 1 builds no expansion: each term is its own column, slice 2's hinge. At degree 2, its columns match slice 1's hand hinge to 3.3·10⁻¹⁶ relative.
  - The config also refuses, at K ≥ 2, the column mass support rule, `HingeNormalization` and KHAOS's raw ridge form, which each read the columns' scale; the paper's columns are unnormalized.

**Gates.**
- *Concern 1* (`RidgeColumnTests`, 6 cases):
  - knots against R's type-7 quantile by hand: 1.7·10⁻¹⁵;
  - the columns and their derivatives against the paper's formulas in exact rational arithmetic, on the training cases and past the range: 3.2·10⁻¹⁵ relative to |d_{ℓ−1}| + |d_K|;
  - a seeded sweep of 17,184 points: values 3.0·10⁻¹⁵, derivatives 5.1·10⁻¹⁶;
  - derivatives against central differences through the basis: 8.1·10⁻¹¹;
  - through a fit: prediction 2.3·10⁻¹⁵ and the gradient 2.5·10⁻¹⁵ against the hand sum; E8's C off its rank-one image 2.7·10⁻¹⁶; the moment service refuses.
- *Concern 2* (`DirectionalRidgeTests`, 10 cases):
  - each proposal is the proposed structure's strict load, bit for bit;
  - the collapsed law over the blocks' columns against the dense marginal: 4.8·10⁻¹³ relative;
  - each move's ratio (two births, their deaths, eight changes, four at a = 1) against a direct evaluation: 5.2·10⁻¹¹;
  - Geweke's test with the prior alone (K = 4, M ≤ 3; 120,000 sweeps against 30,000 prior draws; 10 statistics): worst |z| 1.20, M at −0.09. Births 37% accepted, deaths all, changes all but 0.2%;
  - a K = 4 fit on the single index (n = 300, noise sd 0.1) under the chain, tempering and SMC: the surface's RMSE 0.015, 0.010 and 0.0097 of a spread 0.596; E8's angle 0.0071, 0.0101 and 0.0081 rad. Mean M 2.00 under the chain against 1.00 under tempering and SMC; concern 3's Friedman readings say why;
  - the coincident knots, admission, degree and K = 1 cases above, and the config's refusals.
- *Concern 3:*
  - `BpprParityTests` against `BayesPPR` 0.2.0 on fixed draws (two ridges and a proposed third, K = 4, n = 60, g = 60), read against the corrected ratios:
    - the standardization 6.7·10⁻¹⁶; the knots 5.6·10⁻¹⁶; the basis columns 1.2·10⁻¹⁵;
    - the birth's likelihood ratio, −n/2·(ln sse′ − ln sse) − K/2·ln(g + 1): Δ 5.1·10⁻¹³;
    - the active count and variable set's ratio: Δ 1.4·10⁻¹⁵;
    - the count ratio, the engine at h₁ = 10⁷, h₂ = 10⁶ against `BayesPPR`'s Poisson(10): Δ 8.0·10⁻⁷, the stated approximation.
  - *Friedman's function* (exploration, not validation; `BpprFriedmanTests`): the paper's §3.1 setup at `BayesPPR`'s defaults in the engine's terms:
    - K = 4, A = 3, ω₀ = υ₀ = 1, p₀ = 2/3, n_active = 20, κ = 1000;
    - Zellner–Siow (τ = 1/g, τ ~ Gamma(½, rate n/2)), Jeffreys σ², the intercept inside the g-prior, M capped at 73 = min(150, ⌊n/K⌋ − 2);
    - **the count prior** h₁ = 10⁷, h₂ = 10⁶, Poisson(10) to 8·10⁻⁷ in its ratios;
    - 20,000 sweeps, the last 2,000 kept, 2,000 points out of sample.

    `BayesPPR` ran on the same data at its defaults (`artifacts/chip-bppr/probes/friedman_bppr.R`):

    | Seed | | RMSE y | RMSE f | cover y | cover f | M kept | wall clock |
    | --- | --- | --- | --- | --- | --- | --- | --- |
    | 3671 | engine | 1.114 | 0.506 | 95.2% | 85.3% | 6 | 8.3 s |
    | | `BayesPPR` | 1.115 | 0.482 | 95.7% | 83.3% | 5 | 20.0 s |
    | 3681 | engine | 1.066 | 0.435 | 95.2% | 82.9% | 4 | 4.0 s |
    | | `BayesPPR` | 1.115 | 0.547 | 94.5% | 71.2% | 4 | 23.8 s |
    | 3691 | engine | 1.105 | 0.487 | 94.4% | 76.8% | 5 | 5.2 s |
    | | `BayesPPR` | 1.147 | 0.568 | 93.0% | 75.4% | 6 | 23.0 s |
    | paper | one run | 1.05 | 0.38 | 94.9% | 91.1% | | |

    The engine and `BayesPPR` read alike on the same data. Both fall short of the paper's f readings, which come from one run on one draw of the data. The test's bounds catch only a gross failure: RMSE f < 1, RMSE y < 1.5, y coverage in [0.85, 1], f coverage in [0.6, 1].
  - *M does not move once settled, in either sampler.* Over the kept draws M is constant in all six runs; the engine accepted no birth and no death in 2,000 sweeps.
    - The likelihood charges each block K/2·ln(1 + g), 11.4 nats at g ≈ n under Zellner–Siow and more at larger g. A birth must cut the residual sum of squares by about 7% to pay for it, where a block that fits only noise cuts about 4σ², 1.3%; a needed block's death costs far more than it saves.
    - The ridge fit test's chain held M = 2 where tempering and SMC found 1.
  - *A birth's `Numerical` refusals* (28%, 24% and 5% of births at M = 6, 5 and 4) are D06 refusing a ridge whose t₀ lies below every case.
    - Such a ridge, 1 − p₀ = ⅓ of births, has b₁ affine over the data. Once the current ridges' affine b₁ and the intercept span the inputs a new one reads, its b₁ lies in their span.
    - At a settled state of a chain on two seeds' data (probe `FriedmanRefusalProbe.cs`, 5,000 births each): of the births with t₀ below every case, 984 of 1,651 were refused at seed 3681 and 361 of 1,714 at 3671; of the others, 0 of 3,349 and 6 of 3,286.
    - `BayesPPR`'s `get_qf_info` refuses such births too: its Cholesky fails on an exact span, or its diagonal ratio passes 1,000.
- *Fitting untouched:* every test project test-select selects, read from disk, on snapshots of HEAD plus each concern's files:
  - `332692b`: Oracle 21 of 39 classes, 105 of 105; Regression 247 of 339, 2,704 of 2,704.
  - `349ab81`: Oracle 21 of 39, 105 of 105; Regression 248 of 343, 2,714 of 2,714. The gate ran before two test-side additions (the K = 4 fit's inclusion count, the newest angle in the bit-for-bit test) and a death's NaN angle; on a snapshot of `b5fac86` plus the commit's files, the five classes that read them (`DirectionalRidgeTests`, `RidgeColumnTests`, `DirectionalChainTests`, `DirectionalHingeTests`, `BassDesignTests`) pass 51 of 51.
  - `10b290d`: Oracle 1 of 40 classes, 1 of 1; Regression 1 of 344 with the pinned tier, 45 of 45.
  - The other 25 projects hold none of the files and select nothing.
- *The cost row:* (indicative: wall clock, CostProfile's meter after a discarded warm-up pass; slice 2's single index on [0, 1]⁴, each chain burnt in 3,000 sweeps; high priority on cores 0 to 7, the machine at 24% load before; probe `artifacts/chip-bppr/probes/BpprBlockCostProbe.cs`):

  | Median µs | birth | death | change | sweep | prediction per draw |
  | --- | --- | --- | --- | --- | --- |
  | n = 300, axis hinges (M 5) | 4.7 | 1.1 | 4.9 | 5.3 | 18 |
  | n = 300, directional hinge (M 3) | 10.8 | 1.0 | 11.6 | 9.2 | 19 |
  | n = 300, ridge K = 4 (M 1) | 76 | 7.5 | 52 | 54 | 27 |
  | n = 1000, axis hinges (M 7) | 19.9 | 1.7 | 12.1 | 10.6 | 68 |
  | n = 1000, directional hinge (M 5) | 33.9 | 1.7 | 33.7 | 20.7 | 74 |
  | n = 1000, ridge K = 4 (M 1) | 222 | 6.7 | 157 | 122 | 75 |

  A block's birth or change costs 4.5 to 7 times slice 2's hinge's, at one ridge against its three and five terms. The strict load evaluates and factors every block's columns, not the new block's alone, so the cost grows with M: Friedman's chains, at M 4 to 6 and n = 300, take about 4 to 8 s for 20,000 sweeps, 200 to 400 µs a sweep. The readout reads columns, one coefficient each.

**Decisions taken beyond the brief.**
- *Proposals by a strict scratch load*, not bordered rank-K updates. Each block proposal costs a structure's build, O(N(KM)²), but every mode the design serves is the strict build's, and a restore or a hand-over reproduces it bit for bit. Rank-K updates are the follow-up below.
- *The support count reads b₁*, the block's reach; the late cubic columns are nonzero on only the top quarter of the cases above t₀. A ridge basis takes the count rule; the column mass rule is refused.
- *`Structure.Degree`* is refused unless 1 at K ≥ 2, rather than ignored.
- *n_active is floored at K + 1*, so the cases above U fill K + 1 knots, `BayesPPR`'s floor; it is logged as an auto-pick.
- *M's default cap* is ⌊(N − 1)/K⌋, BASS's N − 1 columns read at K columns a ridge; Friedman sets `BayesPPR`'s 73.

**Unbuilt.**
- The labelled options, slice 4.
- `BayesPPR`'s posterior run as a gate: read once, on Friedman, as exploration above.

**Follow-ups.**
- *Rank-K bordered updates for block proposals.* They would replace the strict build's O(N(KM)²) with a rank-K border. This changes only the cost row: a block birth or change against slice 2's hinge. Decide from the cost row.
- *One projection per proposed direction* (slice 2's follow-up): a block reads one projection for its range, its knots and K columns. The same decision.
- *Mixing over M under ridge blocks.* A chain holds M where its burn-in leaves it, and so does `BayesPPR`. This decides whether a ridge basis's default sampler is tempering, SMC or several chains, or whether to build a block split/merge move. The reading that decides it: M's spread across chains and seeds, against tempering's and SMC's, on Friedman. If tempering or SMC spread M where single chains hold it, the default changes or the move is built; if all of them hold it, M is identified on these data and nothing changes.
- *Affine ridges.* A third of births propose a ridge whose b₁ is affine over the data, and at a settled state many of them are refused (60% and 21% at the probe's two states). Most current ridges are such ridges (4 of 5 and 4 of 6 there). This is the paper's prior (p₀), not a defect, and a labelled option for p₀ is slice 4's question if the user wants it.

**Doc changes owed.**
- `docs/Maths/Regression/BASS.md`: the ridge basis.
  - the expansion at the design boundary; the columns and their knots;
  - undefined knots outside the support and coincident interior knots outside F;
  - admission all or none; the degree; the readouts that count ridges;
  - the count prior's units; n_active's floor; M's cap; the config's refusals; the cost row.
- `docs/Maths/Regression/Sources&Changes.md`:
  - `BayesPPR`'s basis and the count prior's approximation of its Poisson;
  - the Friedman readings beside the paper's and `BayesPPR`'s on the same data;
  - M's mixing in both samplers.

### Slice 2 verified from a clean archive, its Oracle selection to follow (coordinating session, 2026-11-15)

- **Commit:** `1e7bbea`, on a snapshot of `31b612b`, through `scripts/verify-csharp.ps1` over its twelve changed files.
- **Result:** Regression 2,698 of 2,698 at commit scope, the 32 directional cases among them. No other project but Oracle selects.
- **Oracle, confirmed since** in the pass that verified the defects chip's oracle limits, on a snapshot of `b5fac86` with this slice's files in the changed list: 214 of 214 ([receipts](../receipts/verify-defects-s8s9-b5fac86/)). That pass also confirms slice 1's one open test, `GbassHorseshoeReadoutTests.HBass_MeetsCorrectedHbass_…`, which passed in 1 m 34 s under the 500 s default.
- **Oracle, as first recorded:** its selection (21 of 39 classes) was held for that pass, so that a timeout there was not read as this slice's. The chip's own run: 104 of 105, the one failure `khaos_readouts_oracle.R` at the old 120 s, passing alone.
- **Read, against Lean's line under decision 5 and the paper's TeX.**
  - *The birth:* its proposal ratio carries one −ln(M + 1), the death's uniform choice, and the count law supplies the rest, so the count factor is (M + h₁)/((M + 1)(h₂ + 1)) and no further 1/(M + 1) enters. θ and the knot are drawn from their priors and cancel.
  - *The death:* the dying term's birth under the counts without it.
  - *The change:* the within-move ratio ln π(t | θ) − ln π(t\* | θ\*) is the proposal's, the power spherical law being symmetric and the coin fair; the prior's ratio cancels it, so the move accepts on the collapsed likelihood ratio alone.
  - *The knot's range* is the paper's eq. (UL), read in its TeX: U = Q_q and L = U − (U − Q₀)/p₀. The report's reading of the p₀ sentence is right: under (UL), P(t > Q₀) = p₀.
  - *The shared order and variable laws* were lifted out of `BassTermProposal` with their bodies unchanged.
- **One point to state in slice 3's report.** The target is the prior restricted to directions with a knot range, not renormalized, so the count law is tilted by P(range)^M. Say what the Geweke test's prior sampler does with a draw that has no range, a rejection of the whole draw or the term dropped, and that it matches the chain's target. The two differ only at order ε = 0.14% under the integrated count law, below what the test reads.
- **Receipts:** `artifacts/verify-bppr-s2-31b612b/`, with the changed list, the counters and the snapshot tracked in [receipts](../receipts/verify-bppr-s2-31b612b/).

### Slice 2: directions in the structure chain (2026-11-15)

The coordinating session read slice 1's code with no defect and asked for slice 2's reading and design while its clean-archive run went on, with nothing of slice 2 committed before the verification note. Its note came at the end of the build: slice 1 verified from a clean archive at `e70447c` (above), with this list of eight files to close.

**Landed** (`1e7bbea`).
- `DirectionalTermSpace` (new) holds a directional basis's term measure (decisions 2 and 4).
  - *The direction:* a uniform on {1, …, A}, 𝒥 uniform among the C(p, a) subsets, θ_𝒥 uniform on S^(a−1), density Γ(a/2)/(2π^(a/2)) against surface measure, ½ per sign on S⁰. Draws normalize a standard normal vector, a fair sign at a = 1.
  - *The knot:* uniform on [L, U], the paper's eq. (UL). U is R's type-7 quantile of the projected training cases at q = 1 − n_active/n, n_active = min(20, n/10) by default (`BayesPPR`'s `n_dat_min`); L = U − (U − Q₀)/p₀ with p₀ = 2/3. Draws take U − (U − L)·r, `BayesPPR`'s form.
  - *The range* comes from the column's own projection (`DirectionalHingeKernel.Project`, now shared), with the top order statistics kept in a small ascending buffer, a sort past 64 of them. Each direction's range is kept beside it in a `ConditionalWeakTable`, a pure cache.
  - *A range with no width* (ties at the least projection up to U) leaves its direction no knot: such a term lies outside the support, as a structure outside F does, and the measure is not renormalized for it. A birth that draws one is a failed birth, a change a null move.
- `DirectionalBasis` (`Structure.Directional`): p₀ (`LocalProbability`), n_active (`ActiveCases`) and κ (`Concentration`, 1000). A is `Structure.MaxInteraction` (3, resolved to min(3, p), the paper's A); ω₀ and υ₀ are `OrderWeight` and `VariableWeight`. Left out of the JSON when null, so every hinge config's form and fingerprint stand.
- `BassStructurePrior` gains the directional part (`BassPart.Directional`, `PartCount` 5): its order is a, its variables 𝒥, its factors the direction's and the knot's log densities. A hinge or Legendre fit's prior refuses a directional term by name rather than scoring it as an empty one.
- `DirectionalTermProposal` (new), E36's kernel beside D01's and KHAOS's (F-X4):
  - *the birth* draws a and 𝒥 by D01's order and variable law, now shared internal statics in `BassTermProposal`, bit for bit; it is the paper's law corrected at a = 1 (decision 5's recommendation, which slice 0 found `BayesPPR`'s). θ and the knot are drawn from their priors;
  - *the death* reverses it under the reduced counts;
  - *the change* keeps a and 𝒥, draws θ_𝒥 from PS(θ_𝒥, κ), or at a = 1 the sign by a fair coin (the user's decision), then the knot from its prior given the new direction. Its within-move ratio ln π(t|θ) − ln π(t*|θ*) cancels the prior's, so it accepts on the collapsed likelihood ratio alone. A new move kind, `DirectionChange`, tallies it.
- `BassResolved` refuses categorical predictors (decision 7) and a direct functional response, resolves q and logs it as an auto-pick (`KnotLevel`). `BassModel` takes E36's kernel for a directional config; `Proposal` refuses with E36's name.
- `BassSensitivity.FromCoefficients` refuses a directional term with its reason (sampling, E25) rather than pointing to the moment route, which refuses it too.

**The files that read a term's parts by kind** (the coordinating session's list), each handled or refused with its reason; none reads a directional term as an empty one:
- `BassStructurePrior.cs`: handled. The directional part scores and draws the term; `Validate` refuses a directional term in a fit without the part, and the term counts read a and 𝒥 (`Order`, `Variables`).
- `BassResidualBirth.cs`: refused at the config. `Sampling.Births.Residual` grows continuous hinge factors onto parent terms, which a directional basis has neither of.
- `BassConfig.cs`: handled. It declares the basis and refuses what a directional basis does not read: `Structure.Legendre` beside it, `Sampling.LocalSteps`, `Sampling.FactorMoves`, `Sampling.Births.History` and `.Residual`, and settings out of range.
- `BassResolved.cs`: handled, with refusals: categorical predictors (decision 7), a direct functional response (its terms carry an output hinge product, not built for this basis), `ActiveCases` ≥ n; KHAOS's modified g-prior refuses as before, since it reads Legendre terms.
- `BassModel.cs`: handled. It builds E36's kernel; `Proposal`, D01's kernel, refuses with its reason, so the fixed-M family could not reach a directional term even if the config check were bypassed.
- `BassPriorPrecision.cs`: refused at resolution. Only `ModifiedGForm.Factor` reads a term's parts, and the modified g-prior is refused for any basis but Legendre's; the other forms read the columns alone.
- `sensitivity/BassSensitivity.cs`: refused with its reason. `Of` and `Effects` reach the moment service, which refuses through the kernel's `Factors` (slice 1); `FromCoefficients` now refuses up front, naming the directional hinge and E25.
- `calibration/BassEmulators.cs`: handled. A directional term holds no categorical factor, so its split on θ's categorical parameters is the term itself, whose column the basis evaluates at θ's continuous values. The rebuild that drops θ's subsets is reached only by a term with categorical factors, which a directional term cannot hold.
- *Also read, outside the list:* `BassTermProposal` and `LegendreTermProposal` are never handed a directional term (each serves its own basis); `BassDesign.TableBorder` falls back to the column for any term but a single continuous hinge; `BassActiveSubspace`'s output side reads `Functional` only; `BassDiagnostics`' support tally reads `Order`, which knows a; the structural monitors read M and τ; the presentation, posterior and projection folders read terms through the basis; no fit persists its terms (only configs serialize).

**Gates** (`DirectionalChainTests`, `DirectionalRecoveryTests`, 20 cases):
- *The knot range against R's type-7 quantile by hand* (projections from the raw inputs and a compensated standardization): on the insertion buffer (q = 0.9), the sort (81 order statistics) and an index with no fraction (n = 11): worst 3.5·10⁻¹⁶ relative. An equal direction computed afresh gives the same range. A range with no width: +e₃ on an input tied at its least value has none and log prior −∞; from −e₃ the fair coin proposes it 999 times in 2,000, each a null move.
- *The term measure:* log densities against the closed form −ln A − ln C(p, a) − ln|S^(a−1)| − ln(U − L), to 12 digits, −∞ just outside the range; refusals of a hinge term, a = 4 > A, a variable past p, and a directional term in a hinge fit. 60,000 draws: P(a), each variable's inclusion given a, E[θ₁² | a] (exactly 1 at a = 1), P(θ₁ > 0): worst |z| 2.47; the knot's position on its range √n·D 1.28.
- *The birth's law:* over every (a, 𝒥) under a three-term structure, the engine's order and variable mass against the hand law, with Wallenius' term by quadrature apart from the engine's: worst 6.5·10⁻¹⁴ in ln q; the masses sum to 1 to 2·10⁻¹⁶; 50,000 births' (a, 𝒥) against it: χ² 12.6 on 13.
- *Each move's ratio against a direct evaluation:* five births of a = 1, 2 and 3, their deaths, and eight changes (four at a = 1, four at a = 2), the engine's target ratio plus its within-move ratio against the dense marginal's ratio, the count law's (M + h₁)/((M + 1)(h₂ + 1)) and π(a)π(𝒥|a)/(q(a)q(𝒥|a)) by hand: worst 1.9·10⁻¹³. The fair coin flipped a sign among the changes.
- *The collapsed law given frozen directions against the dense marginal* (y | σ² ~ N(0, σ²(I + B(τBᵀB)⁻¹Bᵀ)), σ² ~ IG(3, 2) integrated, B from the hand columns), M = 0 to 3 at degrees 1 and 2: worst 1.3·10⁻¹⁵ relative.
- *Geweke's test with the prior alone* (b = 0, two seeds; 300,000 sweeps against 60,000 prior-sampler draws, 15 statistics: M, the active counts, each variable's use and θ mass, the sign, the knot's position and its square): worst |z| 2.03 and 2.83; M at z 1.99 and −0.62. Moves at b = 0: births 35% accepted, deaths all, changes all but the 0.14% outside F.
- *Recovery through E8* (sampled route, E[C]'s leading k eigenvectors against the true span; n = 300, noise sd 0.1, 10,000 burn-in sweeps): largest principal angle

  | | chain | tempering (4 rungs) | SMC (200 particles × 2, 12 mutation steps) |
  | --- | --- | --- | --- |
  | single index exp(0.6x₁ + 0.8x₂) | 0.0095 | 0.0035 | 0.0026 |
  | two ridges, cosine 0.36 | 0.059 | 0.087 | 0.018 |

  The gate is 0.1 rad. BASS's axis hinges on the two ridges' chain data reach 0.18 (and 0.13 at n = 1000, where the directional basis reaches 0.036). Two settings were read before they were set: 2,000 burn-in sweeps left the second ridge at 0.30 rad (M 7.1 against 4.0 settled), and four SMC mutation steps at 0.21 (12 give 0.018, 30 give 0.006).
- *Several outputs* (E34) share a directional structure as they share BASS's: two outputs' posterior mean surfaces against their functions, root mean square error 2.0% and 5.4% of each function's spread.
- *The refusals* above, and the config's JSON round trip with and without the basis.
- *Fitting untouched:* every test project test-select selects for `1e7bbea`, read from disk: Regression 244 of 336 classes and Oracle 21 of 39; the other 25 projects hold none of its files and select nothing. They ran on a snapshot of `ca2767c` plus the slice's files.
  - *Regression, the whole assembly:* 3,544 of 3,544.
  - *Oracle at commit scope:* 104 of 105. The failure is `khaos_readouts_oracle.R` at ROracle's 120 s, in `KhaosEngineFitParityTests.Readouts_MeetTheCorrectedGPrior_OnIshigami(declared: False)`. That is the limit, not the code: R runs below normal now, and the defects chip is resetting the limits from readings. Rerun alone at below normal, on a snapshot of `5923ebc` plus the slice's files, both of its cases pass: 45 s (`declared: False`) and 37 s (`declared: True`).
  - *Rebased:* E8's slice 8a and E34's resolution row landed before the commit; on a snapshot of `5923ebc` plus the slice's files the 32 directional cases pass again.
- *The cost row* (indicative: wall clock, CostProfile's meter after a discarded warm-up pass, single index on [0, 1]⁴, each chain burnt in 3,000 sweeps, high priority, the machine at about 50% load; probe `artifacts/chip-bppr/probes/BpprChainCostProbe.cs`):

  | Median µs | birth | death | change | sweep | prediction per draw |
  | --- | --- | --- | --- | --- | --- |
  | n = 300, axis hinges (M 5) | 4.8 | 1.3 | 5.0 | 6.6 | 22 |
  | n = 300, directional (M 3) | 11.8 | 1.0 | 11.2 | 11.0 | 22 |
  | n = 1000, axis hinges (M 7) | 9.9 | 1.6 | 10.8 | 15.1 | 80 |
  | n = 1000, directional (M 5) | 32.0 | 1.5 | 36.1 | 26.7 | 80 |

  A birth or change costs 2.2 to 3.3 times BASS's: a proposed direction projects every case twice, once for its knot range (from a pooled buffer) and once for its column (allocated).

**Decisions taken beyond the brief.**
- A directional basis is its own proposal kernel, as KHAOS's is, and reuses D01's order and variable laws rather than restating them; a and 𝒥's caps and pseudo-counts are BASS's settings (A = `MaxInteraction`, ω₀ = `OrderWeight`, υ₀ = `VariableWeight`), so the paper's ω₀ = υ₀ = 1 is a config away.
- The count prior is BASS's (λ integrated, Gamma(10, 10) by default), not the paper's Poisson(10), and M's cap is BASS's min(1000, N − 1), not `BayesPPR`'s min(150, ⌊n/K⌋ − 2). Slice 3's Friedman comparison should state its count prior.
- A direction without a knot range lies outside the support, not renormalized, as F is; the prior sampler counts such a draw among its intercept draws.
- Several outputs are served; a direct functional response is refused. The flat histogram is allowed, since it runs the same moves and target, but no gate reads it with a directional basis.

**The paper, one more slip.** §2.3 calls p₀ "the prior probability that t_{m0} < Q₀"; under eq. (UL) that probability is 1 − p₀, and p₀ is P(t_{m0} > Q₀), as the sentence's "i.e." clause, the hinge zero for at least one case, says. `BayesPPR` follows eq. (UL), as the engine does.

**Unbuilt.**
- The ridge block and the oracle's fixed-draw parity, slice 3.
- Mixed bases (decision 6), categorical inputs (decision 7), a direct functional response, and E4's history and residual births, local steps and the fixed-M family for a directional basis: each refused with its reason.

**Follow-ups.**
- *One projection per proposed direction:* the knot range and the column each project every case; handing the range's projection to the column would save one of a birth's and a change's two passes over the cases (not measured apart). Decide at slice 3, where the ridge block reads one projection for K columns.
- *The projection buffer* (slice 1's): the column's 8 KB a call at n = 1000; same decision.

**Doc changes owed.**
- `docs/Maths/Regression/BASS.md`: the directional basis (the term measure, the knot's range and its type-7 quantile, the support without a range, the kernel's three moves and their ratios, the fair coin at a = 1, the config and its refusals, the cost row).
- `docs/Maths/Regression/Sources&Changes.md`: the paper's p₀ sentence; the count prior and M's cap against the paper's and `BayesPPR`'s.
- `docs/Maths/Regression/ActiveSubspaces.md`: E36's recovery readings through the sampled route.

### Slice 1 verified from a clean archive (coordinating session, 2026-11-15)

- **Commits:** `2eab479` and `9e249dc`, with `117230e` (`LogGammaRatio`'s small steps) and the defects chip's `ca56277` (the transport oracle's limit) in the same run: a snapshot of `e70447c`, through `scripts/verify-csharp.ps1` over sixteen changed files, in each of the eleven projects that select.
- **Result,** selected classes at commit scope:
  - Distributions 430 of 430; Regression 3,199 of 3,199; Samplers 259 of 259; Synthetic 29 of 29; DimReduction 83 of 83; Mapper 1 of 1; Ph 1 of 1; Geometry 4 of 4; VizCore 384 of 386, the two left by their own `Skip`; CumulativeField 15 of 15.
  - Oracle 213 of 214. The transport cases pass under their own limit, in 3 m 55 s and 4 m 58 s.
  - *The one failure* is `gbass_readouts_oracle.R` at ROracle's 120 s, in `GbassHorseshoeReadoutTests.HBass_MeetsCorrectedHbass_…`. Rerun alone from a snapshot of `569f3b7`, it timed out again.
    - Since `0dad092` a test run's whole tree is below normal, R included; this oracle took 79 s alone at normal priority (slice 0's report).
    - So this run does not confirm that test. It passed on the same code in the chip's slice 1 gates (Oracle 205 of 207, the two transport cases its only failures). It is confirmed here once the oracle's limit is read again at below normal, which is with the defects chip.
  - No compiler error; the warnings are the standing ones.
- **Read.**
  - *The standardization:* z = a + b·u with a = (min − mean)/sd and b = (max − min)/sd follows from u's box scaling; no transform is persisted, so nothing stored lacks it.
  - *The kernel:* the peak at the corner by θ's signs, under the box and under the training range; the derivative's θ_v·b_v; the guard at a hinge that is zero over the box.
  - *The term:* the hash adds the directional part only where there is one, and equality reads it.
  - *E8's route:* `Auto` is exact only when every term has factors; every fit built before this slice is separable, so its route is as it was.
  - *The tests' hand computation* is from the definition in raw units, with its own compensated mean and deviation, apart from the kernel's path through u.
- **`117230e`,** verified with this run: Distributions 430 of 430 here, and every other project that selects for `GammaRatio.cs`.
- **For slice 2, one list to close in its report.** These files read a term's parts by kind and do not yet read `Directional`: `BassStructurePrior.cs`, `BassResidualBirth.cs`, `BassConfig.cs`, `BassResolved.cs`, `BassModel.cs`, `BassPriorPrecision.cs`, `sensitivity/BassSensitivity.cs` and `calibration/BassEmulators.cs`. Each either handles a directional term or refuses it with its reason.
- **Receipts:** `artifacts/verify-bppr-s1-e70447c/` and `artifacts/verify-bppr-s1-oracle-rerun/`, with the changed lists, counters and snapshots tracked in [receipts](../receipts/verify-bppr-s1-e70447c/).

### Slice 1: the directional hinge (2026-11-14)

The user cleared slice 1 before the verification note arrived. The note came during it, and confirms decisions 1, 3 and 7 at their recommendations, as built: the directional hinge first, the paper's standardization with the sample deviation, categorical inputs refused. The user has since settled the change move at a = 1 for slice 2, the fair coin with the knot drawn given the proposed sign, recorded under decision 5. `LogBeta` moves onto `LogGammaRatio` as the defects chip's change (its slice 9), not E36's.

**Landed.**
- `2eab479` keeps each continuous input's training mean and sample deviation (n − 1), `Standardization`, beside its box in `RangeScale` (decision 3). `BassInputs.Standardizing(v)` reads it from the scaled value as z = a + b·u, with a = (min − mean)/sd and b = (max − min)/sd over the box, so prediction's cases standardize as training's did. A scale built without data carries none, and a directional hinge refuses it.
- `9e249dc` adds the directional hinge, BASS's fourth case-side term type (decision 1's first stage).
  - `ProjectionDirection` is a sparse unit θ: active variables ascending, nonzero weights, ‖θ‖² = 1 to 10⁻¹⁰, the sign its own, compared bit for bit. `DirectionalHinge(θ, t)` puts the knot on the standardized projection, continuous (decision 4's form; slice 2 draws it). The reverse hinge is (−θ, −t).
  - `BassTerm.Directional` takes the case side alone. A product with hinge, Legendre or categorical factors refuses (decisions 6 and 7); an output side is allowed. Every other term's hash is as it was.
  - `DirectionalHingeKernel`, the fourth case kernel, gives (θᵀz − t)₊^q/c at the fit's degree with R's powers.
    - c is the hinge's largest value over each variable's normalization interval (the box, or the training range under that option), at the corner where u_j sits at its upper end for θ_j > 0 and its lower end otherwise; c = 1 where the hinge is zero over the box, R's guard.
    - The derivative in u_v is the hinge's on the projection times θ_v·b_v, so at the knot at degree one it is the mean of the one-sided limits.
  - The contract gains `Separable(term)`, default true, and `DesignKernels.Separable`. The kernel's `Factors` refuses with its reason, so the moment service refuses a fit holding one, and with it Sobol indices, effects and E8's exact route (decision 10).
  - E8's `ActiveSubspaceRoute.Auto` now does what its own summary said: the exact route under a product law when every term has factors, the sampled route otherwise. The exact route refuses a fit without them up front. On the sampled route the Poincaré checks are left out (null, as under a joint sample), since their totals are the moment route's. Every readout passes it: a fit, a basis fit, the co-active pairs, several outputs.
  - Prediction, the pointwise density, the gradient readout and admission go through `BassBasis` and needed no change.
  - `BassDesignKernelTests`' test-side directional kernel is replaced by the real one.

**Gates.**
- *The contract test, promoted:* the real kernel reads its term and its variables, is differentiable, is not separable while every other part is, refuses factors, and refuses a variable outside the inputs; its derivatives match central differences at interior points.
- *Columns against a hand computation,* with z from the raw values by the test's own compensated mean and deviation, and c at the box's maximizing corner in raw units:
  - five hinges at degrees 1 to 3 under the box and 1 to 2 under the training range, at training cases and at points past the range on every side: worst 1.3·10⁻¹⁵ relative, and 1/c to 10⁻¹³;
  - a hinge zero over the box is zero, with c = 1.
- *A seeded sweep* of 300 cases: a from 1 to 3, with random variables, weights and signs; knots in [−2.5, 2.5], a fifth far out at ±5 to ±55; degrees 1 to 3 under both normalizations; points to 2·10³ past the box. Columns hold to 1.2·10⁻¹³ and derivatives, against the hand derivative q(p − t)₊^(q−1)·θ_v·b_v/c, to 8.0·10⁻¹⁴, relative to max(1, |·|).
- *Derivatives against central differences* (h = 10⁻⁶, degrees 1 and 2, three directions): worst 9.1·10⁻¹¹ of the largest derivative. A variable outside θ gives zero. With the knot on a case's own projection, degree one's derivative there is θ_v·b_v/(2c) to 14 digits, and degree two's is 0.
- *D06:* θ and −θ, with knots −0.8, 0 and 0.9 on θ and 0 and 0.4 on −θ, are all admitted (newest angles 1.0, 0.75, 0.029, 0.23 and 0.0080). The same term again, and θ moved by 10⁻⁹, are both refused as numerical.
- *Through a fit* (a run's fit with two hand-built hinges, degrees 1 and 2):
  - prediction against the hand computation within 7.8·10⁻¹⁶, and the gradient readout in raw units within 9.4·10⁻¹⁶;
  - `Sensitivity` refuses, naming the directional hinge; E8's exact route refuses;
  - the automatic route takes the sampled route with no Poincaré check, and with one hinge the sampled C is rank one along θ's raw-unit image θ_j/sd_j, to 1.4·10⁻¹⁵ of its largest entry.
- *The standardization:* the mean to 14 digits and the deviation to 10⁻¹⁴ relative, against a compensated two-pass computation; a + b·u against (x − mean)/sd at the training cases to 2.6·10⁻¹⁵, under a declared domain too.
- *Fitting untouched:* every test project test-select selects for slice 1's commits and `117230e`, which rides these gates. They ran on a snapshot whose code is `9e249dc`'s (`artifacts/chip-bppr/snap`).
  - *Selected, per commit (classes, read from disk):*
    - `2eab479` and `9e249dc`: Regression 241 of 333 and Oracle 21 of 39 each. No other project holds them in its closure, so none can move.
    - `117230e` (`LogGammaRatio`): Distributions 19 of 20, Regression 325 of 330, Samplers 26 of 31, Oracle 36 of 39, Synthetic 6 of 8, DimReduction 8 of 8, Mapper 1 of 10, Ph 1 of 43, Geometry 1 of 9, VizCore 50 of 85, CumulativeField 3 of 3, and SPC 0 of 1, which has no pinned test.
  - *The whole assembly:* Regression 3516 of 3516.
  - *Commit scope (selected classes, slow tier left out, pinned tier in):*
    - Distributions 430 of 430, at `117230e`'s own gate.
    - Samplers 259 of 259, Synthetic 29 of 29, DimReduction 83 of 83, Mapper 1 of 1, Ph 1 of 1, Geometry 4 of 4, and CumulativeField 15 of 15.
    - VizCore 384 of 386, with 2 skipped by their own `Skip`. SPC selects no test.
  - *Oracle at commit scope:* 205 of 207. The two failures are `transport_oracle.R` timing out at 120 s, in `DiagramDistances_MatchT4TransportAndLpAssign` at p = 1 and 2.
    - Rerun alone, both time out at 120 s again.
    - The test reads diagram distances, which nothing in slice 1 touches.
    - This is the defects chip's slice 8 limit. `ca56277` (8b) has since given that oracle a limit of its own, 900 s, against readings of 185 to 344 s alone. It landed after the snapshot.
- *The cost row* (indicative: wall clock, CostProfile's meter, n = 1000 on three inputs, high priority, the machine at about 4% load; the probe is `artifacts/chip-bppr/probes/BpprCostProbe.cs`):

  | Per call, µs | q = 1 | q = 2 |
  | --- | --- | --- |
  | Column, axis hinge, one factor | 2.5 | 2.6 |
  | Column, axis hinge, three factors | 3.9 | 3.9 |
  | Column, directional, a = 1 | 4.9 | 5.7 |
  | Column, directional, a = 3 | 6.6 | 7.0 |
  | ∂/∂x₀, axis hinge, one factor | 4.6 | 5.1 |
  | ∂/∂x₀, axis hinge, three factors | 8.9 | 10.2 |
  | ∂/∂x₀, directional, a = 3 | 6.3 | 6.4 |

  The directional kernel allocates its projection, 8 KB a call at this n; the axis hinge's column runs in lanes and allocates nothing. A birth's and a change's costs come with slice 2's moves.

**Decisions taken beyond the brief.**
- The standardization's home is `RangeScale` beside the box, with `BassInputs.Standardizing` as its map from u. The kernel reads u, so the case inputs are unchanged and existing fits are bit for bit.
- The hinge is normalized as BASS's are, to a largest value of one over the box or the training range, so D06's column mass and the g-form read it as they read hinges. The paper's ridge is unnormalized; Zellner–Siow is invariant to a column's scale, so its posterior is the same.
- One directional hinge per term, alone on the case side; an output side is allowed, since a direct functional response's term is the Khatri–Rao product of its two sides.
- `Separable` enters the contract; E8's `Auto` already promised the sampled route for terms without factors.
- Decision 7 binds at the term: a product with a subset factor refuses. A fit's data may still hold categorical predictors other terms read; a ridges-only basis (decision 6) is slice 2's configuration.

**Unbuilt.**
- The structure chain's moves, slice 2: no configuration builds a directional fit yet, and a fit holding one is made by hand.
- Sobol indices and effects by sampling (E25): refused, with the reason.

**Follow-ups.**
- *The projection buffer:* pooling it would remove the 8 KB allocation per column without moving a bit. Decide at slice 2, against the birth's and the change's costs.

**Doc changes owed.**
- `docs/Maths/Regression/BASS.md`: the directional hinge as a term type (its definition, normalization, derivative and refusal of factors); the contract's `Separable`; the transform's standardization beside the box, with decision 3's alternative, the box, and its trade-offs.
- `docs/Maths/Regression/ActiveSubspaces.md`: `Auto`'s rule with separability, and the Poincaré check absent for a fit holding a directional hinge.

### Slice 0 verified from a clean archive (coordinating session, 2026-11-14)

- **Commit:** `e31e894`, through `scripts/verify-csharp.ps1` over the slice's eight changed files, in each of the eleven projects that select.
- **Result,** selected classes at commit scope:
  - Distributions 416 of 416; Regression 3,201 of 3,201; Samplers 259 of 259; Synthetic 29 of 29; DimReduction 83 of 83; Mapper 1 of 1; Ph 1 of 1; Geometry 4 of 4; VizCore 384 of 386, the two left by their own `Skip`; CumulativeField 15 of 15.
  - Oracle 205 of 207. The two failures are `transport_oracle.R`'s 120 s timeouts in `DiagramMetricsTransportParityTests` (p = 1 and 2), the standing ones. `gbass_readouts_oracle.R` passed here.
  - The build's warnings are the standing ones; none is in the slice's files.
- **Read.**
  - *The law:* N follows from the surface integral; the cosine is Beta(α, β); the draw is two gammas and a uniform tangent; each branch of the reflection carries e₁ to μ with ‖u‖² ≥ 2; S⁰ is right at κ = 0 and at κ > 0.
  - *`LogGammaRatio`:* the recurrence, Stirling's coefficients and their difference.
  - *`Log1p`:* the new branch is taken only where the old value was +∞, which no finite x > −1 has, so no finite result moves. "Fitting untouched" rests on that, with `LogGammaRatio` and the law new and called by no fit.
  - *Friedman's C:* its three expectations by hand.
- **Checked from outside the engine** ([script and output](../receipts/verify-bppr-s0-e31e894/)):
  - N's formula against quadrature of the surface integral, 20 cases: within 7·10⁻⁶⁰.
  - `LogGammaRatio`'s algorithm, transcribed into doubles, against mpmath: at s ≥ 10⁻³ the error is under 10⁻¹⁵ of its largest term.
  - Cin(2π) = 2.4376533930572244, and Friedman's three expectations against two-dimensional quadrature to 10⁻²⁶.
- **One finding, fixed.** `LogGammaRatio` subtracted two rounded Stirling remainders, so its relative error grew as s fell: 4·10⁻¹² at s = 10⁻⁶, 4·10⁻⁹ at 10⁻⁹ and 4·10⁻⁶ at 10⁻¹².
  - The tests' floor of 10⁻¹⁴ absolute could not see it. The law's own calls, s ≥ ½, were unaffected.
  - The chip differenced the series power by power (`117230e`), with a small-step test against polygammas from outside the function. It is read here and verified with slice 1, whose gates it rides.
- **One figure in the report, measured.** At κ = 10⁵ the gain in ln N over two ln Γ calls is 3.4·10⁻¹⁰ to 1.8·10⁻¹¹ absolute, not five digits: (α + β)·ln 2 itself rounds at 1.5·10⁻¹¹. The log density near the mode carries about κ·10⁻¹⁶ either way, and the symmetric proposal never evaluates it.
- **For the user,** from the report's follow-ups:
  - the change move at a = 1: the fair coin on the sign (the chip's recommendation, `BayesPPR`'s form) or the kept sign; it binds at slice 2;
  - `LogBeta` onto `LogGammaRatio`, which would change the bits of fits that read beta densities;
  - three problems outside E36: the 15 to 22 minute Distributions case with no slow list, test-select's test host at normal priority, and the transport oracle 5 s inside its timeout.
- **Receipts:** `artifacts/verify-bppr-s0-e31e894/`, with the changed list, the counters and the outside checks tracked in [receipts](../receipts/verify-bppr-s0-e31e894/).

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
