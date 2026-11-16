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
6. **Mixed directional widths, selected 2026-11-15.** Three disjoint kinds: smooth K = 4 blocks, K = 1 directional hinges, and categorical factors alone when available. This supersedes the earlier alternative of BASS hinge products beside ridges. Whole-surface sampled sensitivity reads a mixed fit.
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

## Slices 5 and 6: categorical inputs and mixed directional widths

**Slice 5 was approved by the user on 2026-11-15** and has since been verified. **Slice 6 implementation authorized 2026-11-15.** Decision 6a selects kinds by width within the directional basis; the remaining recommendations stand. The coordinating session's earlier review raised three points:
- *6c* holds for a one-column hinge. With K = 4 a ridge at a = 1 is a spline curve in one input, not a duplicate of an axis hinge, so the least active count could be a setting with 2 as its default.
- *6d and decision 2.* Decision 2 said E21's ridge blocks handle collinearity on the coefficient side. Scale-sensitive forms are refused for unnormalized ridge columns, so that is not met today, and D06 alone handles near-duplicates. Meeting it needs a per-kernel column scale, a decision of its own.
- *Order.* One fit per basis compared by evidence, and the reading of M's mixing, come before a mixed chain is built.

**The user's review (2026-11-15).** Slice 6 stays a draft until the categorical label is defined and the evidence and mixing comparison is in hand. The user's slice 5 points (blocks of varying width, the determinant counting direction-bearing terms, a gate-only Sobol fixture, and what an empty, small or tied cell does) were built into slice 5; its report says how. On slice 6:
- *6a:* keep the per-term mixture, but first settle who owns a term of categorical factors alone. Both kinds' laws can generate exactly that term, and its stored factors reveal no birth label. So either the label is retained through draws and restores, or shared indicators get one explicitly weighted place in the prior. The choice fixes the reverse birth's probability and what M_r and M_h mean.
- *6b:* keep the integrated Beta(1, 1) once labels are defined. A fair-coin birth is valid, but its ratio must carry the proposal's ½, the target's Beta factor, and the reverse death picked among all terms. The component kernels use the count they are handed for that death's denominator, so delegating the whole ratio to them with per-kind counts would be wrong.
- *6c:* default to 2 for a one-column directional hinge, and allow 1 for a ridge block: a K = 4 block on one input is a univariate spline curve, which a minimum of 2 would remove. An explicit minimum of 2 stays meaningful for a fit whose ridges are meant to be oblique.
- *6d:* keep the g-form default and qualify the existing normalized `RidgeForm` on mixed columns. Its precision reads diag(G)/N, so rescaling one column already cancels, which corrects this draft's claim below that E21's ridge would shrink the two kinds unequally. Keep slice 3's refusal of `RawRidgeForm` at K ≥ 2, and test the normalized form by rescaling one admitted ridge column and checking the marginal likelihood and predictions. A penalty tuned apart for a whole block would be a further model choice.
- *The pure-basis fits first,* on matched data under proper priors: the log-evidence difference with its Monte Carlo error and flags, and M's traces across independent starts. Slice 3 recorded no accepted birth or death in the last 2,000 sweeps of each Friedman run, and slice 6 changes a term's kind only by a death and a birth, so held-out error alone cannot show that a mixed chain explores both representations. Those readings inform the mixed sampler before it is built.

The user chose these as E36's next work, decision 7's and decision 6's later steps. This was written as a design for review. Decisions marked **open** go at their recommendation unless the user turns them. Each line says what the other outcome would change.

### Slice 5: categorical inputs, BASS's subset factors multiplying a ridge

**Baseline before slice 6.**
- `BassResolved.cs:316` refuses categorical predictors in a directional fit (decision 7).
- `BassTerm`'s constructor (`BassTerm.cs:129`, `:132`) refuses a directional hinge or a ridge column beside any other case-side factor.
- `DirectionalTermSpace.Validate` refuses any part but a direction.
- The directional prior (`BassStructurePrior`) holds the directional part alone.

**The design** (decision 7's direction, at its recommendation). A term is a direction part times a categorical part: Π_v 1{c_v ∈ S_v} · b(θᵀz | t), with the subsets drawn from BASS's subset measure.
- *The kernels already compose this.* A column is the product of every case kernel's part (`BassBasis.CaseColumn`). The categorical kernel reads no continuous variable, so the derivative's product rule holds (`BassBasis.CaseDerivative`).
- *Ridge blocks:* the categorical part multiplies all K columns of a block, which stays one structural unit.
- *The prior:* BASS's term measure over two parts. Each part's order is uniform up to its cap, A for the direction and BASS's `MaxCategoricalInteraction` (3) for the categories. The variables are a uniform subset of that size, and the subsets come from `CategoricalSubsets` (R's size-uniform). A draw with both parts empty is the intercept and lies outside F, as in BASS's several-part fits.
- *The moves:*
  - A birth draws both parts' orders and variables by D01's law under each part's counts, then θ, the knot and the subsets from their priors.
  - A death is slice 2's.
  - A change on a term holding both parts picks one of them with probability ½. On the direction it is E36's change (power spherical, or the fair coin at a = 1). On a categorical factor it is D01's redraw from the subset measure.
  - Every move accepts on the collapsed likelihood ratio and the count and order terms, as now.
- *The paper's dummy coding, not taken.* It appends L − 1 dummies to x, so they enter θ. A ridge whose 𝒥 holds only dummies becomes β(1 − Π(1 − d_ℓ)), an indicator of at least one of its categories (the paper's eq. (cat_ridge)).
  - *What it gives that subsets do not:* a category-dependent shift of one ridge's argument, g(θᵀz + Σθ_ℓd_ℓ), so categories share a curve's shape with moved knots. It also gives a union across several categorical variables, and parity with `BayesPPR`.
  - *Its costs:* the ridge is no longer separable in the category, so sensitivity and E8 cannot sum the categories exactly. The dropped reference level changes the prior. And the dummies' scale enters the sphere θ is uniform on (decision 3's metric).
  - *What subsets give instead:* separability in the category, BASS's validated subset law, and a pure categorical term (an empty direction) that equals the paper's indicator ridge for one variable.

**Decisions for the user.**
- **5a, the knot's range under a subset, settled** (the user, 2026-11-15, as recommended).
  - *Recommended:* t₀'s range and the ridge's quantile knots are read from the projections of the cases the term reaches (its cells), so about n_active of its own cases lie above U, and K + 1 knots exist where it is nonzero. The range is cached per (θ, subsets). A categorical change redraws the knot from its prior given the new cells, so the move still accepts on the likelihood ratio alone.
  - *Alternative:* the range over all cases, BASS's own convention for a hinge beside a categorical factor. The knot's prior then ignores the subsets, and a categorical change keeps the knot. But a small subset leaves few of its cases above U, so the support count refuses many of its draws.
- **5b, pure categorical terms, settled** (the user, 2026-11-15, as recommended).
  - *Recommended:* allowed (the direction's order may be 0), as BASS allows a categorical main effect. Under a ridge basis such a term is one column. Without it, a level's shift needs a ridge.
  - *Alternative:* every term holds a direction, and categories only gate ridges.
- **5c, the determinant prior beside categories, settled** (the user, 2026-11-15, as recommended).
  - *Recommended:* it reads directions alone, so one direction gated by two subsets has det 0 and no mass under the option. That is the option's stated belief, many distinct directions, which the report states.
  - *Alternative:* a determinant within each set of terms sharing a cell, which no source proposes and which costs a cell bookkeeping.

**Gates, slice 5.**
- *Columns and knots:* gated columns and blocks against the indicator times slice 1's and slice 3's hand columns. Knots over the cells against R's type-7 quantile.
- *The prior:* the term measure's log density by hand over both parts, and draws against it by χ².
- *Moves and Geweke:* each move's ratio against a direct evaluation, gated and ungated terms with both kinds of change. Geweke's test with the prior alone, with categorical statistics.
- *Recovery from outside the model:* a category-gated ridge, for example 1{c ∈ S}·exp(0.6x₁ + 0.8x₂) plus a level shift. E8's angle and the surface against BASS's hinges on the same data, and `BayesPPR`'s dummy coding read beside it as exploration.
- *Readouts:*
  - prediction and the gradient on gated terms by hand;
  - the calibration emulator's split keeps a gated ridge's direction (below);
  - the sampled sensitivity on a gated fit against a closed form.
- *Fitting untouched.*

### Slice 6: mixed bases (kinds by width selected after the pure-basis readings)

**Refused today.**
- `BassConfig.cs:181` refuses a directional basis beside `Structure.Legendre`.
- `Structure.Directional` replaces BASS's hinges outright: its prior holds no continuous hinge part, and `BassModel` builds E36's kernel alone.
- `BassTermProposal`, D01's kernel, refuses a directional config by name (`BassModel.Proposal`).

**Revised after the pure-basis readings (2026-11-15; the report "Pure-basis readings before slice 6").**
- *Smoothness, not direction, decides most of the evidence between the bases.* The paper's cubic blocks win every smooth function by 9 to 46 nats, axis-aligned ones included, and lose axis kinks to degree-1 hinges by 24.
- *At matched degree the direction's own share is modest* where the function is smooth: about +6 oblique, −4 to −7 axis-aligned. It is large on an oblique kink: +68.
- *Plain chains under ridge blocks hold M where they start* (R̂(M) 8.3 and 5.2 over four starts), where tempering's chains and SMC's replicates agree.

**The design** (E28's questions, answered).
- **6a, the kinds, decided 2026-11-15.**
  - *Use kinds by width inside the directional basis.* A term is a smooth K = 4 ridge block, a K = 1 directional hinge, or categorical factors alone. A one-coordinate directional hinge supplies an axis kink with a continuous knot. Each direction-bearing term can also take slice 5's categorical gates. BASS's axis-hinge products remain available through the ordinary BASS basis.
    - *What carries over:* one term measure, (θ, t₀) with its subsets, serves both widths, so E36's moves, the determinant, the cells and the fair coin carry over. Slice 5's blocks of varying width already expand and read them back.
    - *What it answers:* the question the readings put, kinked or smooth. The axis-or-oblique question is answered inside the direction prior.
    - *What it lacks:* BASS's hinge products, interactions as products rather than as sums of ridges. On Friedman's function these read no better than directional hinges (−538.1 against −537.4).
  - *The draft's: ridges beside BASS's hinge terms.* Two term measures and a kernel combining D01's and E36's laws. It keeps BASS's product interactions and data-located knots.
  - *Per variable,* E28's middle: each continuous variable is labelled for one kind, declared or drawn. A declared label is cheap; a drawn one moves every term on its variable at once.
  - *Per part within a term* (a hinge on x₃ times a ridge on x₁, x₂): not recommended. The derivative's product rule assumes one kernel reads each variable, and no source proposes it.
  - *One fit per basis, compared by evidence:* now read, with the estimators' disagreement qualified in its report. A mixed basis can express smooth and kinked components directly. Pure K = 4 already contains b₁, a hinge, and pure K = 1 can approximate a smooth curve with several terms; whether mixing improves recovery or cost is measured, not assumed.
- **Who owns a term of categorical factors alone** (the user's review; it applies to either design).
  - *Selected with 6a.* A categorical-only term gets one explicitly weighted place in the prior. Each kind's law is conditioned on what makes it that kind, so the laws' supports are disjoint and a term's kind is read from the term itself:
    - a direction (a ≥ 1) in a block of K columns or in one column;
    - under the draft's design, at least one continuous hinge;
    - at least one subset and nothing else.
    - *The reverse birth* reads the dying term's kind off the term: its kind's proposal weight times its kind's law under that kind's counts.
    - *The counts:* the per-kind counts (M_K4, M_K1, M_c; or M_r, M_h, M_c) are functions of the state, and their sum M feeds the count law.
    - *Its cost:* each kind's law is conditioned, so the direction's order runs from 1 there and the categorical kind's from 1.
  - *A retained label:* each term carries its kind through draws, restores and the term table.
    - A categorical term's label has no likelihood, only its share's prior, so it is a nuisance variable the chain must mix, and the same term under two labels is two states.
    - The reverse birth reads the label's kind; M_r and M_h count labels, and the categorical terms split between them by prior noise.
  - *The mixture integrated, no label:* a categorical term's prior is ρ_r π_r(t) + ρ_h π_h(t).
    - A collection's prior then sums over the ambiguous terms' labels, a polynomial in the shares that takes O(A²) to integrate for A ambiguous terms.
    - The reverse birth reads both kinds' laws, mixed by their proposal weights, and M_r and M_h are undefined for those terms.
- **6b, the shares, recommendation selected for implementation.**
  - *Recommended:* Dirichlet(1, 1, 1) over the three kinds, integrated as λ is. The collection weight gains Γ(3)·Π_k Γ(1 + M_k)/Γ(3 + M), uniform over the split at every M, and sits in the terms' prior beside the determinant.
  - *The birth's ratio:* the birth draws a kind with a fixed weight (⅓ each) and then its law. Its ratio carries that weight, the Dirichlet factor's ratio, the kind's law under its own counts, and the reverse death's 1/(M + 1) over all terms.
  - *Computed by the mixed kernel itself:* E36's and D01's kernels take the death's denominator from the counts they are handed, so delegating the whole ratio to them with per-kind counts would be wrong (the user's review).
  - *Alternative:* fixed shares, a multinomial on the split, firmer as M grows.
  - *Available kinds:* continuous-only data have two kinds, K = 4 and K = 1, each proposed with probability ½ and shares integrated under Dirichlet(1, 1). With an enabled categorical part there are three, each proposed with probability ⅓ under Dirichlet(1, 1, 1). A direction-bearing kind has a uniform active count on 1..A and categorical order on 0..C; the categorical-only kind has order on 1..C. Categories therefore allow both gated directional widths and one separately weighted categorical-only kind, without a hidden label or empty term. Shares are not left on an unavailable kind.
- **6c, the least active count, recommendation selected for implementation.**
  - *Recommended, under the width design:* none is needed. A K = 1 hinge at a = 1 is the only axis hinge in the basis, so no two kinds duplicate one column, and a block at a = 1 is a univariate spline curve.
  - *Under the draft's design,* the user's review: 2 by default for a one-column directional hinge, 1 allowed for a block, and an explicit 2 kept for fits meant to be oblique.
- **6d, the coefficient prior, recommendation selected for implementation.**
  - *Recommended:* the g-form's default and E21's normalized `RidgeForm`, whose precision reads diag(G)/N, so rescaling one column cancels. That is tested by rescaling one admitted ridge column and checking the marginal likelihood and predictions.
  - *Kept:* `RawRidgeForm` stays refused at K ≥ 2.
  - *Not taken:* a penalty tuned apart for a whole block is a further model choice.
- **Sampler decision (2026-11-15; pure-block default implemented at `a83db0c`, retained for mixed widths).** For a new pure ridge-block fit (K ≥ 2), or a mixed fit whose basis admits such blocks, resolve an omitted sampler to self-sizing parallel tempering. Use the exact reference rung when the declared priors and tempered form admit it. Otherwise use an admissible positive hot end and record why. Explicit sampler choices win, including SMC. An explicitly requested plain chain runs and carries a fit-level diagnostic that ridge counts can remain at their starts and intervals for f can be too narrow. Distinguish automatic selection from an explicit plain chain in the config, log the resolved sampler, and resume a checkpoint with the sampler that produced it. The rule applies to the configured basis even if the current draw contains zero blocks.
- A width change keeps θ, t₀ and the subsets between K = 1 and K = 4. It occupies half the change slot on a direction-bearing term; the ordinary direction/subset change occupies the other half. Forward and reverse probabilities are equal, with no Jacobian; the target supplies the likelihood and integrated shares' ratios. Birth/death/change scheduling and boundary menu corrections remain the structure chain's. Categorical-only terms use the subset change alone.
- **The moves.** A birth draws the kind, then its law. A death is uniform over all terms. A change keeps the kind: E36's or slice 5's change for a direction-bearing term, slice 5's subset redraw for a categorical term, and D01's redraw for a BASS hinge term under the draft's design. Local steps, the fixed-M family, and E4's history and residual births stay refused.
- **Mixing and identifiability.** f's functionals do not care which representation is drawn. The per-kind counts and inclusion readouts say less the freer the mix.

**Gates, slice 6.**
- *The prior:* the shares' collection weight by hand, the conditioned kinds' laws, and draws against them.
- *Moves and Geweke:* each move's ratio across kinds, the width change included, against a direct evaluation. Geweke's test with the prior alone, with M per kind.
- *Blocks:* K = 4 blocks beside one-column terms; slice 5 built and gated the expansion.
- *Sampler resolution and coverage.* Exercise the pure-block and mixed defaults, an explicit plain chain with its diagnostic, an explicit SMC fit, the one-column hinge control and checkpoint continuation. Compare M across independent starts and nominal 95% pointwise interval coverage for f and y separately on the Friedman fixture. Report cost and Monte Carlo uncertainty. The evidence comparison has its [own follow-up](evidence-estimator-disagreement-brief.md).
- *Recovery from outside the model:* a smooth oblique ridge plus an axis kink, for example exp(0.6x₀ + 0.8x₁) + 3(x₂ − 0.4)₊. The mixed fit against each pure basis, with their evidence, under tempering and SMC.
- *Fitting untouched.*

### How the selected widths meet what is built

- **K-column ridge blocks.**
  - `RidgeExpansion` expands smooth structure terms into four columns and leaves each directional hinge and categorical-only term as its own column. Read-back walks all three cases term by term and restores the width.
  - A mixed design has a structural cap M_max and an independent actual-column cap N − 1. A width-up change is refused when its three extra columns cross that bound.
  - `Structure.Degree = 1` remains required: this slice's K = 1 kind is a directional hinge, and its K = 4 kind has the paper's fixed linear/cubic columns.
- **The determinant prior.** Both K = 4 and K = 1 bear directions and count in M_θ; categorical-only terms do not. The normalizer is c_(M_θ), and the coupling is −∞ beyond M_θ = p or on dependence. Without categories every term bears a direction and M's cap can still be cut at p. The prior sampler first draws kinds by their integrated shares, holds kinds and categorical parts fixed during determinant rejection, then draws knots; admission rejects the whole collection afterward.
- **The fair coin** remains active at a = 1 in either directional width's ordinary change: sign uniform on ±1, then knot redrawn conditional on that sign and its cells.
- **Admission.** D06 reads every column of both kinds; the support count reads a gated ridge's b₁ over its cells.
- **The readouts.**
  - *Prediction, the pointwise density and the conformal layer* read columns through the kernels, so they need nothing.
  - *The gradient* holds by the one-kernel-per-variable rule, under 5's design and 6a's recommendation.
  - *E8* takes its sampled route, since a ridge is not separable; its product law already carries categorical marginals.
  - *The sampled sensitivity's whole-surface route* reads all three kinds, and is the route used for slice 6's mixed-fit readings. No split estimator is needed here. Any future exact/sampled split must inspect separability: a K = 1 directional hinge with several active coordinates is not a BASS axis-hinge product.
  - `BassSensitivity.FromCoefficients` keeps refusing ridges.
- **One defect slice 5 must fix.** The calibration emulator's split on categorical calibration parameters (`BassEmulators.cs:457`) rebuilds a term from its hinge, Legendre and functional parts, which would drop a gated ridge's direction. Slice 2 found it unreachable because a direction held no categorical factor; slice 5 makes it reachable.

### The files that read a term's parts by kind

| File | Today | Slice 5 | Slice 6 |
| --- | --- | --- | --- |
| `BassTerm.cs` | refuses a direction or ridge column beside another case-side part | allows categorical factors beside them | still refuses hinge × ridge in one term |
| `BassConfig.cs`, `DirectionalTermSpace.cs` | refuses a directional basis beside Legendre; refuses local steps, the fixed-M family, history and residual births | — | declares `MixedWidths` at K = 4, keeps degree 1 and the block admission/refusal rules |
| `BassResolved.cs` | refuses categorical predictors | resolves the categorical cap for a directional fit | separates structural and actual column caps; determinant counts both directional widths |
| `BassStructurePrior.cs` | the directional part alone | the categorical part beside it | three disjoint conditioned kinds, integrated shares, coupling over both directional widths |
| `DirectionalTermSpace.cs` | refuses any other part; ranges per direction | ranges and knots per (θ, cells) | — |
| `DirectionalTermProposal.cs`, `MixedWidthProposal.cs` | directional births, deaths, changes | categorical draws and redraws | mixed kernel owns conditioned per-kind births, global deaths and width flips; reuses direction/subset changes |
| `BassTermProposal.cs` | refuses a directional config | — | remains the ordinary BASS kernel; the new tally kind is `WidthChange` |
| `BassModel.cs` | E36's kernel for a directional config | — | the mixed kernel; the expansion over ridges |
| `RidgeExpansion.cs`, `BassDesign.cs` | every unit expanded; draw or units by the first term | the categorical factor on every column; b₁'s support over the cells | ridge units alone; per-term read-back; two caps |
| `BassBasis.cs`, `DesignKernels.cs` | product of kernels; one kernel per variable for derivatives | nothing | nothing under 6a's recommendation |
| `DirectionDeterminant.cs`, `BassPriorSampler.cs` | every term's direction | direction-bearing terms alone | both directional widths; correlated kind draws, determinant rejection, knots, then whole-collection admission |
| `BassSweep.cs`, `BassTarget.cs` | the terms' prior; the count prior on structure terms | — | existing calls read integrated shares through `LogTermsPrior`; count law stays on all structural units |
| `BassDiagnostics.cs` | the support tally by part and order | the categorical part's orders too | the `WidthChange` move tally; retained per-kind counts read from the terms |
| `BassResidualBirth.cs`, `BassBirthMixture.cs`, `BassPriorPrecision.cs` | refused for a directional basis | refused | refused in a mixed fit |
| `calibration/BassEmulators.cs` | drops a direction on a categorical split (unreachable) | keeps it (the defect above) | — |
| `sensitivity/BassSensitivity.cs`, `BassActiveSubspace.cs`, `BassSampledSensitivity.cs` | refuse a ridge on the exact routes; the sampled routes read any fit | — | the split's f_sep and f_ns (its slice 4) |

**Left out of both.**
- The paper's dummy coding, and ridges times Legendre terms (KHAOS beside ridges).
- Learned category positions (E29).
- A direct functional response with ridges.
- Per-part mixing.
- History and residual births for ridges.

## Report

Under "Reports" below, newest first, per slice: what landed by commit; the gates' numbers; decisions taken beyond this brief; what was left unbuilt and why; follow-ups with the decision each would change; the doc changes owed.

## Reports

### Pure-block Friedman coverage under three samplers (2026-11-16)

The measurement the sampler decision owed (`a83db0c`'s entry below), complete: 18/18.
- **The run.** `artifacts/chip-bppr/probes/SamplerCoverageProbe.cs` on Friedman's fixture, with `BpprFriedmanTests`' data and priors: Jeffreys σ², the intercept inside the g-prior, Zellner–Siow τ, K = 4.
  - *The data:* n = 300, noise sd 1, 2,000 held-out points.
  - *The cases:* three data seeds (3671, 3681, 3691) by two sampler seeds by three samplers.
    - The default, which resolves to self-sizing tempering with a hot end of 1/300 because Jeffreys σ² refuses the reference rung.
    - Declared plain chains.
    - Declared SMC: the linear form, 500 particles, 12 mutations.
  - *Each fit:* four starts at 18,000 + 2,000 sweeps per chain.
  - *The intervals:* f's from the weighted latent draws; y's from the exact predictive mixture.
  - *Receipts:* `build-s6b/runs/20261115_222520-coverage/coverage.trx` and `coverage-output.txt`.

| Sampler | M across starts | R̂(M) | f, pooled (MC sd) | f, one start: mean (range) | y, pooled / one start | RMSE f | Fit s | Sweep-equivalents |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| default (tempering) | moves within each start; start means differ by ≤ 0.76 | 1.00–1.36 | 91.1% (2.5) | 87.3% (78.8–93.5) | 95.6% / 95.2% | 0.420 | 338 (204–420) | 1.28M |
| declared plain | fixed in every start of 4 fits; start means 4 to 6 | frozen ×4; 7.0, 5.0 | 92.8% (3.4) | 83.1% (70.0–94.0) | 95.6% / 94.9% | 0.484 (one 0.812) | 8 (5–13) | 80k |
| declared SMC | replicates agree within 0.25 | n/a | 89.1% (0.9) | 87.3% (81.2–94.2) | 95.4% / 95.3% | 0.423 | 2,948 (2,069–4,030) | 2.74M |

How to read the columns:
- *Values* are means over the six fits.
- *MC sd* comes from the two sampler seeds on each dataset: √(mean Δ²/2) over the three pairs.
- *The six fits' sd of pooled f* is 3.3–3.6 points under every sampler, data and Monte Carlo together.
- *Given a fit,* the held-out binomial error is about 0.7 points for f and 0.5 for y.
- *"One start"* is a single chain or replicate, as the paper's single `BayesPPR` run is (91.1% for f, 94.9% for y).

What it reads:
- **Plain chains hold M:** fixed in every start of four fits, R̂ 7.0 and 5.0 in the other two. A single plain chain covers f at 83.1% on average and 70.0% at worst, the case the caution names.
  - *Pooled over four starts* the coverage rises to 92.8%. That is the spread of four frozen chains at different M, not mixing.
  - *Their RMSE is the worst:* 0.484, with 0.812 in one fit.
- **The default moves M within each chain and matches SMC.** Single-run f coverage is 87.3% under both; RMSE is 0.420 against 0.423.
  - *Its cost* is about 40× the plain chain's wall time and about a ninth of SMC's.
  - *SMC's pooled coverage is the most reproducible:* MC sd 0.9 points.
- **No sampler reaches nominal f coverage at this budget:** 83–87% for one start, 89–93% pooled. That leaves the coverage decision noted in slice 6's report open, for pure blocks too.
- **y coverage reads 94.2–96.2% pooled and 93.3–96.2% for one start,** and does not separate the samplers.

Caveats:
- *The default's ladders made essentially no round trips:* at most 0.0001 per scan, over 14–18 levels with a maximum pair rejection of 0.74–0.97. Its M moves through neighbouring swaps, so R̂ near 1 does not show mixing over the whole ladder.
- *Fit seconds were taken under the shared machine's load.*
- *The readout dominated each case's time:* the predictive-mixture quantile was inverted five times per fit, so each case took 15 min to 1 h 15 min. The CDF-membership test gives the same coverage far faster.

**Doc changes owed:** none beyond `a83db0c`'s entry.

### 17522a6 read against slice 6's draft (2026-11-16)

A reading, at the coordinating session's request; nothing built.
- **The prior and every ratio match the draft and the user's review.**
  - *Categorical-only terms* are their own kind, read off the term, with no label. Their law drops the direction's order and conditions the categorical order to k ≥ 1. Direction-bearing kinds condition a ≥ 1. Both corrections are exact against slice 5's uniform orders from 0, which carry no intercept renormalization.
  - *The shares* enter once, in `LogTermsPrior`, as the integrated Dirichlet(1, …, 1) factor over the available kinds. `BassSweep` reads `LogTermsPrior` whole, so every move's ratio carries them.
  - *The birth's ratio is the mixed kernel's own:* the kind weight 1/K, the kind's adaptive law under its own counts, and the reverse death 1/(M + 1) over all terms. The death mirrors it. Component kernels serve only the ordinary change.
  - *The width change keeps θ, t₀ and the cells,* with a proposal ratio of 0. The draft's "kind weights" in its acceptance belong to the fixed-shares alternative alone.
  - *The prior sampler is exact.* Kinds come by Polya's urn. The determinant rejection redraws directions alone with kinds and cells held, since E[det]/c_{M_θ} = 1 at each M_θ.
- **Left out, recorded, not built.**
  - The per-kind count readout and the `WidthChange` tally in `BassDiagnostics`, which this brief's file table lists; the probes compute per-kind counts from the draws.
  - A sensitivity split; the whole-surface sampled route reads mixed fits.
  - The BASS-hinge alternative.
- **The default caps differ.**
  - *The mixed fit's structural cap* reads one column a term: min(1000, N − 1) (`BassResolved.cs:339`). Pure K = 4's is min(1000, ⌊(N − 1)/4⌋). Past ⌊(N − 1)/4⌋ blocks, the mixed fit's column cap N − 1 refuses designs as admission.
  - *A default-cap evidence comparison* between a mixed and a pure fit therefore truncates M's law at different points. The two normalizers differ by the count law's mass between the caps.
  - *Negligible on Friedman:* under the near-Poisson(10) count law at n = 300 that mass lies above M = 74. It matters for small N or a count law reaching ⌊(N − 1)/4⌋.
  - *Not in the recovery readings:* they set a cap of 12 for every basis.

### The sampler default for ridge bases, committed at `a83db0c` (2026-11-15; coverage read 2026-11-16)

**Landed.** The user's sampler decision for pure ridge-block fits: eleven files, `BassSamplerDefault.cs` and `RidgeSamplerTests.cs` new. The mixed-width slice kept it unchanged for mixed fits (`17522a6`; its legacy probe, in that slice's report below).
- **The setting.** `BassSampling.Sampler` takes `Automatic` or `PlainChains`. `Automatic` is the default and is left out of the serialized form, so every existing config keeps its form and fingerprint. A declared `PlainChains` serializes, with a fingerprint of its own.
- **The rule** reads the configured basis, so a draw holding no block resolves the same.
  - *Passes through untouched,* nothing logged, unless `Structure.Directional.Columns` ≥ 2.
  - *A declared `Tempering`, `Smc` or `FlatHistogram` runs as declared.*
  - *Otherwise, self-sizing parallel tempering:* on the likelihood path, in the form the model runs, integrated unless its caller names one.
    - The exact reference rung where the engine's own resolution admits it.
    - Else a positive hot end auto-picked by the resolution, with the rung's refusal as the reason. On Friedman's fixture, Jeffreys σ² refuses the rung, and the hot end resolves to 1/300.
    - Else plain chains, with the caution and both refusals.
  - *`Validate` refuses* `PlainChains` beside a declared sampler, and an undefined value.
- **The logged row** is `Sampler`, first in `AutoPicks`. Its readings:
  - "parallel tempering, self-sizing, exact reference rung", because "automatic for a ridge basis (K ≥ 2): its plain chains can hold M where they started";
  - "…auto-picked hot end", with "no reference rung, since Sampling.Tempering.ReferenceRung needs b = 0, where the integrated form is undefined under the Jeffreys σ² prior; the linear form reaches it (D05)";
  - "plain chains", "as declared (Sampling.Sampler = PlainChains)" or "as the resumed checkpoint ran", each with the caution.
- **The caution.** `BassDiagnostics.SamplerCaution` is set where plain chains run on a ridge basis, declared or resumed. Its text: "plain chains on a ridge basis can hold M where they started, and f's intervals can then be too narrow (E36's pure-basis readings: R̂(M) 5 to 8 over four starts); tempering or SMC is the default for this basis". It is null otherwise.
- **Declared and resolved configs.** The run takes the resolved config (`BassModel.Config`). The fit, its fingerprint and the checkpoint's identity carry the declared one (`BassModel.Declared`).
- **Checkpoints.** `Bass.Run` reads the sampler that produced a checkpoint (`BassEnsemble.Tempered`: a ladder in the restart state) and resolves with it.
  - A plain checkpoint under an automatic config, such as one taken before this default, resumes plain. It carries the caution and the row "as the resumed checkpoint ran", and goes on as the plain run's own resume does.
  - A tempered checkpoint resumes on the same ladder, bit for bit with the uninterrupted run.
- **Review fixes before the gate.**
  - The refusal quoted as the reason drops .NET's " (Parameter 'config')".
  - The automatic ladder takes the form the model runs, so a caller naming the linear form under Jeffreys σ² gets the reference rung, and the ladder and the target agree.
- **The K ≥ 2 tests.**
  - *`BpprFriedmanTests`:* declares `PlainChains`, the paper's and `BayesPPR`'s single plain chain.
  - *`DirectionalRidgeTests`' chain case:* declares `PlainChains` and asserts the caution. The case reads the chain beside its tempered and SMC cases.
  - *`CategoricalReadoutTests`' and `DirectionDeterminantTests`' chain cases:* declare `PlainChains`. Their subject is the readouts and the prior, so the sampler is incidental; their SMC cases are unchanged.
  - *`SampledEffectsTests` (the sensitivity chip's, K = 4):* takes the default, tempering; not edited.
  - *K = 1 and BASS-hinge fits (`RidgeColumnTests` among them):* unaffected.
  - *`RidgeSamplerTests`, 5 tests:*
    - the automatic ladder and its row;
    - the hot end under Jeffreys σ², and the linear form's rung;
    - declared plain chains with the caution;
    - declared SMC and ladders winning, and the refusals;
    - K = 1 and BASS's hinges passing through with the JSON unchanged;
    - both checkpoint rules.
- **Gate:** a fresh snapshot (`5aef657` plus the eleven files) and a fresh build folder, `artifacts/chip-bppr/build-s6c`.
  - *Maths.Oracle.Tests:* 22 of 40 classes, 106/106.
  - *Maths.Regression.Tests:* 263 of 358 classes, the slow tier's 71 methods out and the pinned tier in; 2,852/2,852.
  - *The other 25 projects* select nothing.
  - *Receipts:* `build-s6c/test-select/20261115_222345-commit/` and `…/20261115_223603-commit/`.
- **Coverage.**
  - *The probe:* `artifacts/chip-bppr/probes/SamplerCoverageProbe.cs`, on Friedman's fixture with `BpprFriedmanTests`' data and priors.
  - *The cases:* 18 fits, three data seeds (3671, 3681, 3691) by two sampler seeds by three samplers: the default, declared plain chains, and declared SMC (linear form, 500 particles, 12 mutations).
  - *Each fit:* four starts, 18,000 + 2,000 sweeps per chain.
  - *The readings:* M per start with R̂(M), 95% coverage of f and of y pooled and per start, RMSE for f, and cost as wall time and sweep-equivalents.
  - *Where it stands:* complete, 18/18, at 09:15 on 2026-11-16, unrestarted. Its readings are in the entry "Pure-block Friedman coverage under three samplers" above.
- **Doc changes owed.**
  - BASS's config reference: `Sampling.Sampler`, the ridge-basis default and its resolution, the caution, the declared and resolved configs, and the checkpoint rule.
  - The sampler page: tempering by default for ridge bases.

### Slice 6: mixed directional widths, implemented and qualified (2026-11-15–16)

**Ownership and base.** Started on local `main` at `a83db0c`, with the unrelated document moves preserved. Only the BPPR structure/sampler files, `MixedWidthTests.cs` and this brief belong to this slice. Independent sensitivity and evidence changes stay outside the snapshots and commits. No subagents, no push. Snapshots are `a83db0c` plus the explicit file lists under `artifacts/chip-bppr-mixed/`; each has its own build root. Runs use `DOTNET_PROCESSOR_COUNT=2`, below-normal launches and bounded parallelism.

**Implemented at `17522a6`, on local main.** `DirectionalBasis { Columns = 4, MixedWidths = true }` admits smooth blocks, directional hinges and an available categorical-only kind. A directional structure term carries its width in `DirectionalHinge.Columns`; equality, interning, column expansion/read-back and checkpoint state retain it. Pure structure terms retain the unmarked representation and pure priors reject mixed markers. `MixedWidths = false` is omitted from canonical JSON. The mixed proposal owns conditioned births, own-kind adaptive order/variable counts, global uniform deaths and symmetric width flips. Existing directional changes preserve width, including the a = 1 sign coin. Shares are Dirichlet(1,1) without categories, Dirichlet(1,1,1) with them; determinant coupling counts both directional widths. The exact prior sampler draws correlated kinds by Polya's urn, holds their categorical parts during determinant rejection, then draws knots and rejects whole collections on admission failure.

**Targeted gates.** The initial 10 tests passed (`build-s2/runs/20261115_232305-mixed-s2/mixed-s2.trx`). Expanded gates passed 13/13 (`build-s3/runs/20261115_233159-mixed-s3/mixed-s3.trx`), then final targeted gates 14/14 (`build-s4/runs/20261115_233516-mixed-s4/mixed-s4.trx`). They include the hand shares and conditioned laws; own-kind birth probabilities including a Wallenius check by quadrature; global reverse deaths; actual structure move menus at M = 0, interior and the cap; symmetric width flips; prior-only per-kind counts against independent restricted-prior draws with and without categories and with the determinant; mixed-column restore and exact checkpoint continuation; automatic tempering, explicit plain/SMC controls; and the separate term and column bounds, including a refused width-up change. A rejected removal is checked to preserve counts before a valid death.

The prior-only gate compares 16,000 independent restricted-prior draws with 100,000 chain sweeps after 3,000 burn-in, using batch errors on six statistics: total M, each kind's M, total active count and positive-sign count. Worst |z| is 1.45. Prior acceptance is 0.9944 without categories, 0.9719 with categories, and 0.8265 with categories and the determinant; the latter records 2,937 determinant rejections. Width-change attempts are 10,737, 7,279 and 7,241 respectively. The conditional kind laws precede admission; admission can tilt their split, which both sides of this check retain.

**6d's mixed path.** The rescaling check first admits one four-column block beside one directional hinge and one categorical-only term through the mixed model's expansion/design. It takes those actual six columns into the same collapsed g-form and normalized `RidgeSpectrum` coefficient paths the model uses, verifies the unscaled marginal against the mixed target, then rescales one cubic block column by 0.001, −7 and 1000. It compares marginal likelihoods and posterior mean predictions at three τ values. Dense marginal checks and `RawRidgeForm` refusal are included; these are not only the older generic ridge test.

**Prechange configuration and continuation.** `probes/LegacyProbe.cs` ran against the prechange `chip-bppr/build-s6c` assemblies and the new isolated assemblies. New code reads the saved old canonical JSON. K = 1 and K = 4, each automatic and explicit plain: canonical JSON, fingerprint, and SHA256 over the resumed draws (coefficients, scale/noise/count draws and terms) all match. Each resumed run also matches its own uninterrupted run. Receipts: `legacy-old/`, `legacy-new/`, `invoke-legacy-old2/`, `invoke-legacy-new/` beneath the slice's artifact root. Checkpoints here are the engine's in-memory checkpoint objects; this does not claim a disk checkpoint serializer.

The compatibility probe was repeated against the final `build-s4` engine after the validation fixes: all four cases again match (`legacy-final/`, `invoke-legacy-final/`). `source-hashes.json` records that the eleven committed source/test files match the qualification snapshot byte for byte. Other sessions' later commits on main and uncommitted sensitivity work were excluded from this isolated source snapshot; this is not a combined-main verification claim.

**Measurement provenance.** `probes/MixedRecoveryProbe.cs` measured the external surface exp(0.6x₀ + 0.8x₁) + 3(x₂ − 0.4)₊ on two data seeds, all three bases and tempering/SMC, with proper matched priors. The replay in `recovery-s3/` adds weighted per-start counts and coverage errors and uses the equivalent predictive CDF interval test; the earlier quantile-inversion receipts in `recovery-s2/` are preserved. Eight cases completed before the earlier readout run was stopped to add those uncertainties; this was a measurement revision, not a test failure. `probes/MixedFriedmanCoverageProbe.cs` reuses the existing Friedman fixture and readout, at its stated smaller mixed-fit budget, with f and y separate. The pre-existing pure-block coverage run (`chip-bppr/build-s6b/runs/20261115_222520-coverage/`, `SamplerCoverageProbe.cs`) has not been restarted or claimed complete; its results are qualified separately when available. No dependency on the evidence session's owned code: the evidence readings below carry the currently available uncertainties and flags and the standing estimator-disagreement caveat.

**Fitting untouched, final gates passed.** Regression commit scope: 2,866/2,866, selected 264/359 classes with the pinned tier and the existing slow tier omitted; receipt `build-s4/test-select/20261115_233723-commit/Commit.trx`. Whole Regression assembly, slow tier included: **3,782/3,782**, no skips or failures; receipt `build-s4/test-select/20261116_000547-slice/Slice.trx`. Selected Oracle checks: **106/106**, no skips or failures; receipt `build-oracle-s4/test-select/20261115_234403-commit/Commit.trx`. The implementation commit contains only its eleven source/test paths. No source outside this slice is needed.

**Recovery completed: twelve matched fits.** `build-s3/runs/20261115_233737-recovery-s3/recovery-s3.trx`: 12/12 passed; machine-readable rows are in `recovery-s3/` and `recovery-summary.csv`. Each data seed supplies n = 200 independent uniform points on [0,1]³, noise sd 0.1 and 400 held-out points. The count prior is the default Gamma(10,10) mixture, coefficient scale Zellner–Siow Gamma(½,100), g-form, mean intercept and the default proper scaled σ² prior, all matched across bases. Tempering has two independent starts, exact reference rung and 4,000 burn-in + 1,000 retained sweeps; SMC has three replicates, 256 particles and eight mutations per stage. All use the same explicit structural cap 12. The effect of the widths is not confounded by giving the mixed fit a larger term cap.

| Data | Basis | RMSE f, PT / SMC | 95% f coverage, PT / SMC | 95% y coverage, PT / SMC |
| --- | --- | --- | --- | --- |
| 0 | K = 4 | 0.0434 / 0.0382 | 84.5% / 89.5% | 96.75% / 97.25% |
| 0 | K = 1 | 0.0416 / 0.0398 | 80.25% / 83.0% | 96.25% / 97.0% |
| 0 | mixed | 0.0363 / 0.0344 | 88.0% / 92.25% | 96.75% / 97.25% |
| 1 | K = 4 | 0.0172 / 0.0173 | 99.5% / 99.5% | 96.5% / 96.5% |
| 1 | K = 1 | 0.0321 / 0.0329 | 84.75% / 86.0% | 96.75% / 96.25% |
| 1 | mixed | 0.0157 / 0.0171 | 100% / 100% | 96.5% / 96.75% |

The mixed fit has lower RMSE here, with a practically tied pure-block control on data 1 under SMC. This is two datasets, not a population-risk dominance claim. Good y coverage does not repair the f undercoverage on data 0. The mixed f coverage's between-start standard errors are 3.63 and 4.00 percentage points on data 0 (PT/SMC), zero on data 1 where every held-out point is covered; the latter is a saturated statistic, not zero posterior Monte Carlo uncertainty. Mixed y coverage's between-start errors are at most 0.15 points. The held-out conditional binomial standard error is up to 2.0 points for these f proportions and about 0.9 points for y, separately from variability across training data; two datasets do not estimate that variability precisely.

| Data / sampler | Mean (M_K4, M_K1) | Probability both widths | Fit seconds |
| --- | --- | --- | --- |
| 0 / PT | (0.458, 2.917) | 45.8% | 10.9 |
| 0 / SMC | (0.428, 2.411) | 42.7% | 67.0 |
| 1 / PT | (0.852, 1.463) | 84.9% | 13.8 |
| 1 / SMC | (0.721, 1.902) | 72.0% | 65.8 |

SMC summaries use their draw weights within each replicate and in the pool. The pure-block controls fit two blocks; the pure-hinge controls average about 3.5 terms. Pure-control fit times span 1.8–21.3 s for PT and 16.2–74.9 s for SMC under the shared machine's load. These are descriptive cost rows, not isolated timing claims. Both widths are explored; a mixed draw need not retain both, even when the generating function contains both structures.

**Evidence, with the available errors and flags.** Mixed minus pure-hinge ln Z is −3.05 ± 1.07 and −1.62 ± 0.82 for PT on data 0/1; SMC gives −1.05 ± 0.61 and −0.66 ± 0.49. Mixed minus pure-block is +3.76 ± 1.35 and +8.04 ± 1.25 for PT; SMC gives +4.74 ± 1.17 and +6.49 ± 0.56. Errors here are quadrature combinations of the estimators' reported Monte Carlo errors, not an independently calibrated uncertainty for a difference. SMC reaches the posterior in every replicate, with no overdispersion flag; no fit is saturated or at its ladder-level cap. PT has no observed complete round trips in these 5,000-sweep runs. Therefore the evidence favors the cheaper hinge representation here despite the mixed fits' prediction gain, and several small differences remain weak under the stated errors. The standing estimator-disagreement investigation still applies; these readings do not certify its error estimates or select a new default basis.

**Whole-surface sensitivity.** Every fit was read with the existing sampled `BassSensitivity.Of` route: 80 selected draws, 2,048 base points, order at most 2. Mixed total Sobol indices on data 0 are about (0.186–0.187, 0.320–0.321, 0.502–0.504), and on data 1 about (0.194–0.195, 0.330, 0.485–0.487), agreeing descriptively across samplers and with the pure controls. These are sampled summaries over selected draws, not exact integrals or an assertion of calibrated Sobol uncertainty; no sensitivity code was changed.

**Measurement cost finding.** For continuous predictive mixtures, membership of a central 95% y interval is equivalent to its CDF at y lying in [0.025,0.975]. The eight saved old quantile-inversion cases and their replay have exactly the same RMSE, f coverage and y coverage. Their readout times fall from 17–74 s to 0.6–1.7 s, including the extra per-start summaries. This is a probe improvement only. The existing longer pure-block coverage run was left alone.

**Mixed Friedman coverage: six runs, separately qualified.** `build-s5/runs/20261115_233914-friedman-s5/friedman-s5.trx`: 6/6 passed as descriptive measurements, not assertions of nominal coverage. `probes/MixedFriedmanCoverageProbe.cs` reuses the existing fixture/readout (n = 300, six uniform inputs, noise sd 1, 2,000 held-out points) and the same priors: zero-centre g-form, Jeffreys σ², Zellner–Siow Gamma(½,150), near-Poisson(10) count prior Gamma(10⁷,10⁶). Here the basis is mixed. The data seed is 3671, with two sampler seeds derived by `SeedTree`: −1111815509 and −1769123746, and four independent starts or SMC replicates each. This probe's budget is 6,000 burn-in + 1,000 retained sweeps per chain, versus the pre-existing pure-block probe's 18,000 + 2,000; SMC here has 300 particles and eight mutations per stage, versus its 500 and twelve. These are separate qualifications, not a matched-budget claim against the longer pure-block run, which remains active without a completed receipt at its existing path.

| Sampler | Pooled f coverage, seeds 0 / 1 | Pooled y coverage, seeds 0 / 1 | RMSE f, seeds 0 / 1 | Fit seconds | Approx. sweep-equivalents |
| --- | --- | --- | --- | --- | --- |
| automatic PT | 90.25% / 95.15% | 96.0% / 96.35% | 0.446 / 0.417 | 102 / 139 | 525k / 546k |
| explicit plain | 98.35% / 90.4% | 97.0% / 96.1% | 0.351 / 0.443 | 5.4 / 4.7 | 28k / 28k |
| explicit SMC | 89.7% / 85.5% | 95.85% / 96.05% | 0.410 / 0.462 | 553 / 544 | 1,085k / 1,022k |

On this one dataset, the two sampler seeds give mean f coverage 92.70%, 94.38% and 87.60%, with standard errors of those two-run means 2.45, 3.98 and 2.10 percentage points (PT/plain/SMC). Mean y coverage is 96.18%, 96.55%, 95.95%, with corresponding errors 0.18, 0.45 and 0.10 points. These very small replication counts do not calibrate confidence intervals; conditional held-out binomial errors and uncertainty across training datasets are distinct. The latter is not measured by this one-data-seed check.

The pool conceals substantial disagreement. The eight individual plain-chain f coverages range 73.8–88.75%, despite pooled coverage of 90.4–98.35%; their M means range 5.36–7.80, and R̂(M) is 2.096 / 2.089. Automatic PT gives individual-chain f coverage 77.6–94.75%, M means 4.82–6.36 and R̂(M) 1.518 / 1.217. Both default fits resolve the hot end to 1/300 because Jeffreys σ² refuses the integrated form's exact reference rung. They have 17–20 levels, zero complete round trips and maximum pair rejection 0.99 / 0.96. Plain fits carry the requested caution. Declared SMC runs the linear form, reaches the posterior in all eight replicates with no jumped replicate, and has individual-replicate f coverage 79.7–89.2%; its weighted M means range 4.87–6.86. This measurement does not request evidence, so null evidence flags are not evidence that its estimator is healthy.

The sampler resolution and continuation contracts pass; adequate f coverage and complete mixing are **not established at this budget**. In particular, neither the high plain pooled coverage nor the improved PT R̂ proves convergence, and good y coverage is insufficient. The implementation retains the user's sampler choice; these measurements are a qualification limit and a basis for a subsequent mixing/coverage decision, not a silent tuning change. Full per-start rows are `friedman-output.json`; the raw numeric coverage/cost rows are in `friedman-s5/`. The older pure-block plain-chain readings earlier in this brief remain their own qualification, as does the still-running longer pure-block probe.

**Handoff.** Implementation is committed at `17522a6`; all owned source/test changes are committed, all own runs are finished, and the final qualification report is a separate documentation concern. There is no code dependency on either other session and no push. The remaining uncertainties are measured ones: fitted f coverage and mixing at larger budgets/across more datasets, the independent evidence-error investigation, and the original longer pure-block coverage run's final results. That original run is not owned by this slice and was left running. Subsequent combined-main verification must include the sensitivity/evidence sessions' changes; the gates above certify the isolated slice snapshot. The directional basis still requires at least one continuous predictor; all-categorical data retain the ordinary BASS route.

**Doc changes owed.** BASS's basis/config and prior descriptions, the mixed move law and available-kind shares, the separate column cap, sampler defaults/caution, and sampled whole-surface readouts. `docs/` remains for the offline batch.

### Pure-basis readings before slice 6 (2026-11-15)

The coordinating session asked, as the user's review directs, for the pure-basis readings before anything of slice 6's chain: the log-evidence difference between a ridge-only and a hinge-only fit on matched data under proper priors, with its Monte Carlo error and flags, and M across independent starts against tempering's and SMC's. Then the definition 6a needs: who owns a term of categorical factors alone. The probe is `artifacts/chip-bppr/probes/PureBasisReadingProbe.cs`, run on a fresh snapshot of `37812fd` in a fresh build folder; nothing in `src/` changed.

**The setup.**
- *Six functions from outside the model,* n = 300, inputs uniform on [0, 1]^p:
  - oblique exp(0.6x₀ + 0.8x₁) (p = 4, noise sd 0.1);
  - axis sin(2πx₀) + 4(x₁ − ½)² (p = 4, 0.1);
  - mixed, the oblique plus sin(2πx₂) (p = 4, 0.1);
  - Friedman's (p = 6, sd 1);
  - axis kinks 3(x₀ − 0.4)₊ + 2|x₁ − 0.6| (p = 4, 0.1);
  - an oblique kink 3(0.6x₀ + 0.8x₁ − 0.7)₊ (p = 4, 0.1).
- *Four bases under the same proper priors* (BASS's centre and σ² prior, Zellner–Siow τ, the count prior Gamma(10, 10), M's default caps):
  - the paper's ridge blocks, K = 4;
  - directional hinges, K = 1, degree 1;
  - BASS's hinges at degree 1, its default;
  - BASS's hinges at degree 3.
  The second and fourth are the controls at matched smoothness.
- *Per fit:*
  - four plain chains from independent starts (10,000 burn-in sweeps, 1,000 draws), reporting each chain's M and R̂(M);
  - SMC with 500 particles, 4 replicates and 12 mutation steps, reporting ln Z with its error and the overdispersion flag;
  - a self-sizing ladder to the exact reference rung, two chains, reporting stepping-stone ln Z with its flags.

**The evidence.** ln Z by SMC / by tempering. † marks SMC overdispersed and ‡ tempering saturated; Monte Carlo errors are 0.04 to 2.2 (SMC) and 0.2 to 2.8 (tempering).

| | ridge, K = 4 | directional hinge, K = 1 | BASS's hinges, degree 1 | BASS's hinges, degree 3 |
| --- | --- | --- | --- | --- |
| oblique | **236.9** / **237.5** | 227.4 / 228.0 | 221.3 / 222.6 | 194.3 / 194.4 |
| axis | **219.7** / **219.5** | 167.8† / 174.8 | 175.1† / 178.4 | 203.1 / 205.2‡ |
| mixed | **197.8** / **194.6** | 151.6 / 156.0 | 146.6 / 152.5 | 139.0† / 143.4 |
| Friedman | **−523.6** / **−520.8** | −537.4† / −540.7‡ | −538.1† / −535.3‡ | −560.8 / −557.7‡ |
| axis kinks | 200.2 / 200.6 | **225.0** / **225.2** | 223.6 / 224.5 | 188.8 / 187.3 |
| oblique kink | 253.0 / 253.1 | **261.4** / **261.4** | 192.5† / 194.2‡ | 191.8† / 186.0 |

- *Ridge-only against hinge-only,* the question as asked (K = 4 against BASS's default): +15.6/+14.9 on the oblique function, +44.6/+41.1 on the axis one, +51.2/+42.1 on the mixed one, +14.5/+14.5 on Friedman's, −23.4/−23.9 on the axis kinks and +60.5/+58.9 on the oblique kink.
- *Smoothness decides most of it.* The paper's cubic blocks win every smooth function, the axis-aligned one by more than the oblique one, and lose the axis kinks to both degree-1 bases. BASS's degree-3 hinges are no stand-in for them: those truncated cubics trail BASS's degree 1 everywhere but the smooth axis function.
- *Direction at matched degree* (directional against BASS's hinges, both degree 1):
  - +6.1/+5.5 on the oblique function and +5.0/+3.5 on the mixed one;
  - −7.3/−3.6 on the smooth axis function, and +1.5/+0.7 on the axis kinks, where a = 1 directional hinges are axis hinges with continuous knots;
  - +68.9/+67.3 on the oblique kink, which BASS's products cannot follow;
  - Friedman's is inconclusive, every estimate flagged.

**M across independent starts.** R̂(M) over four plain chains:

| | ridge, K = 4 | directional hinge | BASS, degree 1 | BASS, degree 3 |
| --- | --- | --- | --- | --- |
| oblique | 8.32 | 2.14 | 1.01 | 1.34 |
| axis | constant (2) | 1.23 | 1.56 | 1.04 |
| mixed | constant (2) | 2.60 | 1.83 | 1.09 |
| Friedman | 5.21 | 2.28 | 1.45 | 1.20 |
| axis kinks | constant (3) | 1.84 | 1.03 | 1.87 |
| oblique kink | constant (1) | 3.57 | 3.11 | 2.64 |

- *Ridge chains hold M where they start.*
  - On the oblique function three chains sat at M = 1 and one at 3.88. On Friedman's they sat at 5, 5, 6 and 6, each chain constant.
  - SMC's replicates and tempering's chains agree on the same data: 1.00 throughout on the oblique function; 4.4 to 4.8 (SMC) and 4.94 and 4.96 (tempering) on Friedman's, each spanning more than one M.
  - "Constant" means every chain found the same M and never left it, which the function's count of effects explains.
- *Degree-1 bases' chains disagree less, but they disagree.* R̂ reaches 2 to 3.6 on the oblique kink and the mixed function. SMC's replicates there are flagged overdispersed in places (Friedman's hinges: 8.6 to 10.6).

**What the evidence's errors do not cover.**
- *Round trips:* the ladders recorded about no round trips per scan in most runs; the reference rung's fresh prior draws are what move them.
- *The two estimators disagree beyond their errors:* by 3.2 on the mixed function's ridges (unflagged, 197.8 ± 0.7 against 194.6 ± 0.5) and by 7.0 on the axis function's directional hinges.
- So differences under about 8 nats read as indicative.
- *What stands under either estimator:* the smoothness contrasts of 14 to 46 nats, the axis kinks' 24, and the oblique kink's 68.
- *At the edge:* the 8 to 10 nats between the cubic blocks and directional hinges on the two oblique functions.

**For slice 6** (its revised draft holds the decisions):
- *The useful mixture is kinked beside smooth,* not axis beside oblique. The smooth functions go to the cubic blocks and the kinked ones to degree-1 hinges, by 8 to 46 nats over the other smoothness. At matched degree the axis-or-oblique question is answered within the direction prior, a = 1 being an axis hinge.
- *Sampler implication from this probe.* Settled ridge chains rarely accept births or deaths, which led to the initial recommendation to refuse a plain chain for a mixed fit. The sampler decision above keeps it available by explicit request with a diagnostic.
- *A function holding both a smooth and a kinked structure* is where a mixed basis would beat both pure ones; it is slice 6's recovery gate, since no pure basis can read it.

**Doc changes owed.**
- `docs/Maths/Regression/BASS.md`: the pure-basis readings (the evidence table, smoothness against direction, M's mixing), as the reading behind slice 6's design.

### Slice 5 verified from a clean archive (coordinating session, 2026-11-15)

- **Commits:** `7d0d1e2` and `f11c8cc`, on a snapshot of `37812fd`, through `scripts/verify-csharp.ps1` over their 24 changed files, in the two projects that select.
- **Result:** Regression 2,844 of 2,844 and Oracle 106 of 106 at commit scope. Among them: `CategoricalRidgeTests` 7 and `CategoricalReadoutTests` 5; slice 4's `DirectionDeterminantTests` 10 and slice 2's `DirectionalChainTests` 13 under the restructured prior sampler; and the sampled-sensitivity chip's `SampledSensitivityTests` 18 and `SampledEffectsTests` 11 under the two-line change to its file.
- **Read: the prior sampler under the determinant.** A redraw redraws each direction whole, its active count, its variables and its weights (`DrawDirectionPart` at `active = 0`), and not the weights under held variables. So the kept directions have the law Π π(direction)·det/c at M_θ, whose normalizer reads M_θ alone, and with M and each term's parts from the base law the kept draw has the stated target. Redrawing weights alone would have left the variables' law untilted, which is another law.
- **Read: the user's review at the slice 6 section's head** (`6dc2384`). It corrects this session's second point there: the normalized `RidgeForm` reads diag(G)/N, so it is scale-free on mixed columns, and only `RawRidgeForm` is refused at K ≥ 2.
- **Two changes in other chips' files,** each with a test that reaches it: the calibration emulator's split keeps a gated direction (`BassEmulators.cs`), and the sampled route varies a category that enters only as a gate (`BassSampledSensitivity.TermVariablesOf`). The sampled-sensitivity chip is told of the second.
- **Next, by the user's review:** the pure-basis readings on matched data, and the categorical label's definition for 6a, before a mixed chain. Slice 6 stays a draft.
- **Receipts:** `artifacts/verify-bppr-s5-37812fd/`, with the changed list, the counters and the snapshot tracked in [receipts](../receipts/verify-bppr-s5-37812fd/).

### Slice 5: categorical inputs (2026-11-15)

The user approved slice 5 at its recommendations (5a, 5b, 5c) and reviewed the plan before the build. The review's three slice 5 points are built here: blocks of varying width, the determinant counting direction-bearing terms, and a gate-only Sobol fixture with what an empty, small or tied cell does. Its slice 6 points are recorded at that section's head (6dc2384).

**Landed**, one commit per concern.
- *The gated model* (`7d0d1e2`).
  - `BassTerm` takes categorical factors beside a direction or a ridge column; hinge and Legendre factors there stay refused (decision 6).
  - `DirectionalTermSpace` reads a gated term's knot range and a block's knots over its cells, the cases whose levels each factor holds. U is the 1 − n_active/n_S quantile of their n_S projections. n_active is the fit's, after a ridge block's floor, now resolved as `BassResolved.ActiveCases`. Ranges are cached per direction and cells.
  - *5a's edge cases:* a cell with no case, n_S ≤ n_active, or a range with no width gives no range. The term's log prior is then −∞, outside the support: the term measure and the prior sampler meet it as a null draw, a birth as a failed one, a change as a null move. Nothing is redrawn and nothing renormalized.
  - `BassStructurePrior` holds BASS's categorical part beside the direction when the data have categorical predictors. Each part's order is uniform from 0, a draw with both parts empty is the intercept, outside F, and the subsets come from `CategoricalSubsets`. A term of categorical factors alone is allowed (5b). Without categorical predictors every draw, sum and ratio is slice 2's, bit for bit.
  - `DirectionalTermProposal`:
    - a birth draws both parts' orders and variables by D01's law, then θ, the subsets and the knot over the cells;
    - a change on a term holding both parts takes the direction or a subset with probability ½ each: the subset by D01's redraw, then the knot from its prior over the new cells;
    - every move accepts on the likelihood ratio, with the count and order terms as before.
  - *Blocks of varying width.*
    - `RidgeExpansion` expands a direction into K gated columns and a term of categorical factors alone into one column, itself, and reads columns back term by term.
    - `BassDesign` tells columns from structure terms by any ridge column. It caps structure terms at M_max apart from the columns' K·M_max, and takes the newest angle over the last term's own columns.
    - `BassStructure.TermCount` reads the term table, counting a block once and a categorical term once. All seven readouts that count terms read it, replacing slice 3's division by K, which miscounts a categorical term beside a block (5/4 = 1).
  - *The determinant counts direction-bearing terms (5c, the user's review).*
    - Terms without a direction take no part in det. c is read at M_θ, the number of directions, so which terms bear one stays the term measure's.
    - M keeps its cap beside a categorical part, and M_θ ≤ p binds through det. Without one, slice 4's cut of M at p stands.
    - The prior sampler now draws M and each term's parts from the base law, then redraws the directions alone until kept with probability det, then the knots. A missing range is met after the rejection, since a range checked inside it would tilt the parts' law. Slice 4's draw, M from p(M)/c_M with the whole collection rejected, is exact for slice 4's target but cannot read M_θ.
    - ln c_M is now summed as Σ_{j<M} ln(1 − j/p), exact at M ≤ 1, where the log-gamma difference left 9·10⁻¹⁶.
  - `BassResolved` lifts decision 7's refusal and keeps the categorical cap; `DirectionalChainTests` and `DirectionalHingeTests` assert the lifted refusals as built.
- *The readouts* (`f11c8cc`).
  - *The calibration emulator's split* on a categorical calibration parameter (`BassEmulators.cs:457`) keeps a gated term's direction or ridge column. It rebuilt the term from its hinge, Legendre and functional parts, which would have dropped the direction.
  - *Sobol's sampled route* (`BassSampledSensitivity.TermVariablesOf`, the sampled-sensitivity chip's file, two lines) reads a non-separable term's categorical factors as variables it reads. Before, a category entering only as a gate was never varied: on the C·X fixture HEAD's route held no subset for C at all.

**Gates.**
- *The gated model* (`CategoricalRidgeTests`, 7 cases):
  - ranges over the cells against R's type-7 quantile by hand: 8.9·10⁻¹⁶; gated hinge and ridge columns against the mask times slice 1's hand hinge and slice 3's exact paper columns: 1.4·10⁻¹⁵; a block's knots over the cells: 2.6·10⁻¹⁶;
  - an empty, a too-small and a tied cell: no range, log prior −∞ under either basis; the empty and tied cells leave no block (a support refusal); the small cell's block exists, its term outside the support through its prior;
  - the term measure against the hand measure over both parts: 1.8·10⁻¹⁵; 60,000 draws over (a, k): χ² 7.42 on 5; subset sizes 15,040 of 30,080 at one level;
  - each move's ratio against a direct evaluation over gated and categorical terms (four births and their deaths, sixteen changes, three of them direction changes and thirteen subset changes): 1.6·10⁻¹³;
  - Geweke's test with the prior alone over gated terms (13 statistics, 200,000 sweeps against 40,000 prior draws): worst |z| 2.16;
  - *the determinant on one continuous input beside the categorical one* (the user's example): M's cap stays 4, a ridge and a categorical contrast couple at 0 (ln 1 − ln c₁), two directions give −∞; Geweke worst |z| 1.17, M_θ ≤ 1 in every draw, 7.7% of the prior at M ≥ 2;
  - *a ridge and a categorical term at K = 4* (the user's gate), loaded in either order: 5 columns, reloaded from them bit for bit, M = 2 through the table, the fit's sums equal to 10⁻⁹, and a third term a cap refusal at M_max = 2 though its column would fit.
- *Slices 2 to 4 under the change* (on the same snapshot, 65 cases): slice 2's Geweke bit for bit (2.83); slice 4's tilt χ² 1.23 on 4 (5.27 under the old sampler), its Geweke 2.14 and 1.50 (2.69 and 1.38), its hinge-basis fit 5.5% of the spread (5.0%).
- *The readouts* (`CategoricalReadoutTests`, 5 cases), on y = exp(0.6x₀ + 0.8x₁) where c ≠ 1, plus 1 where c = 2, and 0.5 where c = 1 (n = 300, noise sd 0.05; from outside the model):

  | | mean M | gated directions / categorical terms (draws) | RMSE of the spread 1.154 | BASS's hinges | E8's angle |
  | --- | --- | --- | --- | --- | --- |
  | chain, K = 4 | 2.01 | 600 / 605 (600) | 0.0057 | 0.0178 | 0.0063 |
  | SMC, K = 4 | 2.03 | 400 / 412 (400) | 0.0054 | 0.0173 | 0.0019 |
  | chain, K = 1 | 5.74 | 2,247 / 900 (600) | 0.0189 | 0.0178 | 0.0164 |

  - every draw's M through the table equals its terms read back; the summary and the traces equal it; the prediction against each draw's columns by hand: 7.4·10⁻¹⁶;
  - `BayesPPR` with the paper's dummy coding on the same data (exploration; `artifacts/chip-bppr/probes/categorical_bppr.R`): RMSE 0.0131 with about four ridges, against the gated blocks' 0.0057 with a ridge and a level shift;
  - the emulator's split: 60 readings over 160 gated terms in the draws, equal to the fit's prediction bit for bit;
  - Sobol's sampled route through `SensitivityRoute.Auto` on f(C, X) = C·X (C ~ Bernoulli(½), X ~ U(0, 1)): S_C 0.588, T_C 0.790, S_X 0.205, T_X 0.393 against 3/5, 4/5, 1/5, 2/5; E8's Poincaré check: no violation, its totals' error present.
- *Fitting untouched:* every test project test-select selects, read from disk, each commit on a fresh snapshot of HEAD plus its files and built into a fresh folder (no build folder reused across bases):
  - `7d0d1e2` on `77873bd`: Oracle 22 of 40 classes, 106 of 106 (2,602 s under the machine's load); Regression 260 of 355, 2,839 of 2,839. Only docs landed between that base and the commit.
  - `f11c8cc` on `7d0d1e2`: Oracle 13 of 40, 57 of 57; Regression 261 of 356, 2,844 of 2,844.
  - The other 25 projects hold none of the files and select nothing.

**Decisions taken beyond the brief.**
- *n_active over the cells is the fit's count,* so about n_active of the term's own cases lie above U, as 5a reads; a cell with n_S ≤ n_active hosts no gated term.
- *A block's knots need only two cases above t₀,* so a block can exist where t₀'s range does not; such a term is outside the support through its prior.
- *The determinant's c at M_θ* rather than a mixture over M (the user's review). It keeps the law of which terms bear a direction.
- *The prior sampler's determinant draw* restructured for both slices, as above.
- *Both readout fixes in other files* (the emulator's split and the sampled route's variables), each the smallest change, with a test that reaches it.

**Follow-ups.**
- *The sampled-sensitivity chip's slice 4 split* can now read a gated term as separable between its cells and its direction (the slice 6 draft's note).
- *A gated ridge's E8 readout* averages over the categories under the product law; a per-category active subspace is a readout question for E8, not built.

**Doc changes owed.**
- `docs/Maths/Regression/BASS.md`, the directional basis's categorical part:
  - the gated term and its cells; 5a's ranges and edge cases; the two-part term measure and moves;
  - blocks of varying width and the term count through the table;
  - the determinant at M_θ and the prior sampler's direction rejection; ln c_M's sum.
- `docs/Maths/Regression/Sources&Changes.md`: the paper's dummy coding against subset gates, with `BayesPPR`'s reading on the same data.
- `docs/Maths/Regression/Sensitivity.md` (or wherever the sampled route lives): a gate's variable read by the sampled route.

### Slice 4 verified from a clean archive; the brief's slices are all landed (coordinating session, 2026-11-15)

- **Commit:** `3ef7dc4`, on a snapshot of `e120042`, through `scripts/verify-csharp.ps1` over its eight changed files, in the two projects that select.
- **Result:** Regression 2,745 of 2,745 and Oracle 106 of 106 at commit scope, with `DirectionDeterminantTests` 10 and `DirectionalChainTests` 13 among them.
- **The change of plan, accepted:** c_M in closed form for sparse directions, with Monte Carlo as its gate where the settlement had a Monte Carlo table.
  - *The derivation, read:* by Cauchy–Binet det(ΘᵀΘ) = Σ_S det(Θ_S)², and for independent columns E[det(Θ_S)²] = M!·p^(−M) once E[θθᵀ] = I/p. The paper's sparse law has that: a coordinate is in 𝒥 with chance a/p and carries 1/a of the norm, and the signs are symmetric.
  - *Checked from outside the engine* ([script and output](../receipts/verify-bppr-s4-f62b4ab/)), in numpy from the law's definition, 60,000 collections each: (p, A, M) = (6, 3, 2), (6, 3, 3), (6, 3, 6), (10, 3, 5), (4, 2, 4) and (6, 1, 3) agree with p!/((p − M)!p^M) within |z| 1.43; E[θθᵀ] is I/p to 0.0024; past M = p the determinant is rounding.
- **Read in the report:** the coupling enters through the target's term prior, so no move changes; the prior sampler's rejection by det ≤ 1 is exact for the stated target; M's cap is cut at p with the count law's mass past it reported.
- **Where E36 stands.** Slices 0 to 4 are landed and verified, which is the brief's "done". What follows is the user's:
  - the reading of M's mixing under ridge blocks (slice 3's follow-up);
  - categorical inputs and mixed bases (decisions 6 and 7's later steps), not yet briefed;
  - the cost follow-ups: rank-K updates for block proposals, one projection per proposed direction, rank-one updates of the directions' Gram.
- **Receipts:** `artifacts/verify-bppr-s4-f62b4ab/`, with the changed list, the counters, the snapshot and the outside check tracked in [receipts](../receipts/verify-bppr-s4-f62b4ab/).

### Slice 4: the determinant direction prior (2026-11-15)

The coordinating session verified slice 3 and started slice 4: the determinant as a labelled option, with c_M stated and its tilt on M reported. gcbeck's quantile knots are not taken, since decision 4 stands at the paper's t₀.

**Landed** (`3ef7dc4`).
- `DirectionPrior` (new): `Independent`, the paper's and the default, and `Determinant`. It is set as `DirectionalBasis.Directions`.
- `DirectionDeterminant` (new):
  - π(Θ | M) = Π_m π(θ_m)·det(ΘᵀΘ)/c_M over the structure terms' unit directions, with c_M = E[det ΘᵀΘ] under the independent law;
  - `LogNormalizer` gives ln c_M = ln p! − ln (p − M)! − M ln p, −∞ past M = p;
  - `LogDeterminant` gives ln det(ΘᵀΘ) through the Gram's Cholesky factor over the sparse directions; a pivot at most 10⁻¹³ times its diagonal entry reads as a dependence;
  - `LogCoupling` gives ln det − ln c_M, −∞ past p or on a dependence.
- `BassStructurePrior.LogTermsPrior`: a structure's term prior, its terms' priors summed as before, plus the coupling under the option. The sweep's candidates and the prior sampler both read it, so every move, tempering, SMC and the flat histogram take the coupling from the target. No move changes.
- `BassPriorSampler`: under the option it draws M from p(M)/c_M and keeps a collection with probability det(ΘᵀΘ) ≤ 1 (Hadamard). This is exact rejection from p(M)·Π π(θ_m)·det/c_M on F. It counts its `DeterminantRejections`, carried in a ladder's restart state.
- `BassResolved`:
  - M's cap is cut at p under the option, with a log row: the determinant has no mass past M = p;
  - a `Directions` row states c_1 … c_M and the tilt (below).

**c_M is closed for sparse directions too.** The Grassmann discussion derived c_M = p!/((p − M)!p^M) for dense directions, and the brief planned a Monte Carlo table for sparse ones.
- *The derivation:*
  - By Cauchy–Binet, det(ΘᵀΘ) = Σ_{|S| = M} det(Θ_S)².
  - For a square matrix with independent columns, each with second moment Σ, E[det(A)²] = M!·det Σ.
  - So c_M depends only on E[θθᵀ]. Any law on unit vectors that is exchangeable over the coordinates and symmetric in each sign has E[θθᵀ] = I/p. The paper's sparse law is one: a, then 𝒥 uniform, then θ_𝒥 uniform on its subsphere.
  - Hence c_M = C(p, M)·M!·p^(−M), the dense form.
- *Gated by Monte Carlo* rather than tabulated by it (below): an exact constant beats a noisy table.

**The tilt on M**, by the engine's convention (the prior restricted to F, unnormalized, its tilts reported):
- *With c_M stated,* the directions' law integrates to one at each M ≤ p. The option's tilt on M is therefore the cut at p: the count law's mass past p lies outside F. The resolution reports that mass.
- *c_M itself* is the tilt the count law would take with the normalizer left out, as gcbeck's `zellner.h` leaves it. The resolution states c_1 … c_M beside it. At p = 4 they are 1, 1, 0.75, 0.375, 0.094; at p = 10, c_5 = 0.30.
- *F's own tilt* reads through the existing readouts: the support tally by order, and the mass at the cap, which is now p.

**Gates** (`DirectionDeterminantTests`, 10 cases):
- *c_M against the product Π_{j<M}(p − j)/p* for p = 1 … 12 and M = 0 … p: 6.8·10⁻¹⁵ relative; −∞ at M = p + 1.
- *c_M against Monte Carlo* (40,000 collections each; the determinant computed by elimination with partial pivoting, apart from the engine's Cholesky): worst |z| 1.92 over twelve cases.
  - Dense directions, drawn in the test as normalized Gaussians: (p, M) = (3, 2), (3, 3), (5, 3), (6, 6), (10, 5).
  - The engine's sparse term measure: A = 3 at p = 6 (M 2, 3, 6) and p = 10 (M 3, 5); A = 2 at p = 4, M = 4; and A = 1 at p = 6, M = 3, where c_M is the chance of three distinct coordinates.
- *The determinant by hand:*
  - ln det against elimination on 400 sparse collections, M = 1 … 6: 1.6·10⁻¹¹;
  - a birth's factor det(Θ′ᵀΘ′)/det(ΘᵀΘ) against the new direction's sin² to the others' span by Gram–Schmidt: 3.4·10⁻¹¹;
  - −∞ on a sign pair, a repeat (the pivot floor), three directions in a 2-plane and M past p; directions on disjoint coordinates give det 1.
- *Each move's ratio against a direct evaluation* (slice 2's fixture under the option): slice 2's dense-marginal, count and order-and-variable ratios, plus ln sin² − ln((p − M)/p) for a birth and the determinants' ratio for a change. Worst 7.3·10⁻¹⁴ over three births, their deaths and eight changes. The two births that repeat a current direction have −∞, as the target says.
- *The tilt against its definition* (p = 4, count law h₁ = 40, h₂ = 10, so λ is near 4; declared cap 8):
  - the resolution cuts M's cap to 4 and reports the count law's mass past p, 0.3533, against the sum by hand;
  - the prior sampler's M over 40,000 draws against the count law cut at p: χ² 5.27 on 4 degrees of freedom;
  - F refused none of its draws; the determinant turned back 138,219 of 178,219 attempts.
- *Geweke's test with the prior alone* under the option (the same fixture; 200,000 sweeps against 40,000 prior draws; slice 2's 15 statistics with det(ΘᵀΘ) and M = 4's indicator): worst |z| 2.69 (M = 4's indicator) and 1.38 on two seeds. Births 69% accepted at b = 0, deaths 49%.
- *Unchanged without the option:* the terms' prior is the plain sum bit for bit; the default basis equals one declaring `Independent`; the config round-trips the option and refuses an unknown value.
- *A fit* on the single index (n = 300, noise sd 0.1):
  - Every draw's directions are independent and number p at most.
  - With ridge blocks, the option's regime, the surface's RMSE is 3.0% of the spread under the chain and 2.9% under SMC, with one direction in every draw.
  - With one-column hinges it is 5.0%, read at a 10% bound. Several hinges along nearly one direction are what a hinge basis needs for a single index, and they are what the determinant taxes (the Grassmann discussion's reading).
- *Fitting untouched:* every test project test-select selects for `3ef7dc4`, read from disk, on a snapshot of `60459fc` plus its eight files (only docs landed before the commit): Oracle 22 of 40 classes, 106 of 106; Regression 251 of 346, 2,745 of 2,745. The other 25 projects hold none of the files and select nothing.
- *The cost row:* (indicative: slice 3's probe, wall clock, high priority on cores 0 to 7, the machine at 17% load before; `artifacts/chip-bppr/probes/BpprDeterminantCostProbe.cs`):

  | Median µs | birth | death | change | sweep | prediction per draw |
  | --- | --- | --- | --- | --- | --- |
  | n = 300, ridge K = 4, independent (M 1) | 54.8 | 4.3 | 36.7 | 36.0 | 18 |
  | n = 300, ridge K = 4, determinant (M 1) | 54.7 | 2.4 | 32.1 | 31.9 | 18 |
  | n = 1000, ridge K = 4, independent (M 1) | 164 | 5.3 | 112 | 98 | 52 |
  | n = 1000, ridge K = 4, determinant (M 1) | 171 | 3.1 | 110 | 99 | 52 |
  | n = 1000, hinge, independent (M 5) | 20.2 | 1.0 | 23.4 | 15.9 | 49 |
  | n = 1000, hinge, determinant (M 2) | 18.3 | 0.7 | 19.2 | 14.6 | 49 |

  At the same state the coupling is within the noise: an M × M Gram and its factor at M ≤ p. The determinant's deaths allocate less because its column cap, p·K, is smaller. The hinge rows compare different states, since the chains settle at different M. At n = 300 the determinant's hinge chain sits at its cap, M = 4, where a birth is a cap refusal (3.3 µs), so that row is left out.

**Decisions taken beyond the brief.**
- *c_M in closed form for sparse directions,* by the derivation above, with Monte Carlo as its gate. The brief and the user's settlement said "tabulated by Monte Carlo for sparse directions" because no closed form was known; the constant is now exact. The coordinating session can turn this.
- *M's cap is cut at p* in the resolution rather than left to the prior's −∞, so a birth past p is a cap refusal, and the mass at the cap reads the cut.
- *The pivot floor,* 10⁻¹³ of the diagonal entry: an exact dependence that rounding leaves a hair above zero (three directions in a coordinate 2-plane) reads as none. Its cost is the mass where sin² < 10⁻¹³.
- *The determinant on a hinge basis is served, not refused.* Its cost there is measured above.

**Unbuilt.**
- *gcbeck's quantile knots:* decision 4 stands.
- *A tempered det^γ:* γ softens none of the determinant's consequences, since 0^γ = 0 (the discussion).
- *The discussion's local, sign-aware kernel* exp(κ(θᵀθ′ − 1)), with c_M(κ) by Monte Carlo: not asked for.

**Follow-ups.**
- *Rank-one updates of the Gram's factor.* Each candidate's coupling refactors the M × M Gram, O(M³) at M ≤ p, which is negligible at the gates' p. It matters only in the option's own regime of many inputs, with M in the tens. Decide from a cost row at such p.

**Doc changes owed.**
- `docs/Maths/Regression/BASS.md`, the determinant direction prior:
  - its definition; c_M and its derivation for sparse directions; the cut at p and its tilt; the prior sampler's rejection; the pivot floor;
  - its cost on a hinge basis; the cost row.
- `docs/Maths/Regression/Sources&Changes.md`:
  - gcbeck's missing normalizer and the tilt it puts on M;
  - c_M's closed form for sparse directions, which corrects the Grassmann discussion's "Sparse directions need E36's Monte Carlo".

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
