# Ridge blocks on the BASS core: E21's remaining scope (chip brief)

Date: 2026-11-08. Status: written by the GBASS chip at the user's request (relayed by the main session, 2026-11-08), for the main session's verification before slice 0. It turns what [E21](extensions-brief.md#a-prior-precision-builder-with-ridge-blocks) leaves open into a diff against the engine and an order to build it in. Reports go under "Reports" below, newest first.

## Why

KHAOS built E21's first parts: the prior-precision builder (`f289af3`, `BassPriorPrecision`) and random-effect intercepts on one grouping with s held, integrated out by whitening the rows by V = I + s²(ZH)(ZH)ᵀ, the effects drawn from their conditional under the observation-weighted sum-to-zero convention (`0a2168c`, `BassRandomIntercepts.cs`; [KHAOS brief, slice 3](khaos-brief.md#slice-3-random-effect-intercepts-chip-2026-11-08)). What is left, as the main session listed it:
1. **s drawn.**
2. **Several groupings,** crossed and nested, each with its identifiability convention declared.
3. **Random slopes** on a declared continuous input, as ridge blocks, under a declared standardization: ridge priors change under rescaling.
4. **The block under GBASS's noises,** where the weights and latents move V each sweep.
5. **The block under conformal prediction:** the candidate seam for a case in a known group and in a new one, its shape agreed with the conformal chip.

## What is already built

| Needed | Here |
| --- | --- |
| A block integrated out exactly | `RandomInterceptBlock`: V^(−½)x = x + Z·M·Zᵀx at O(n + L²) a column, ln\|V\| = Σ ln(1 + s²λ) on the working response's log-Jacobian, so every marginal, rung, particle and the evidence carry it |
| The effects' exact conditional | `RandomInterceptBlock.DrawEffects`, drawn last in each composed draw (`BassSweep`), `BassDraw.GroupEffects` |
| The grouping and the declaration | `BassData.Grouping`, `BassConfig.RandomIntercept` (`BassRandomIntercept(Scale, Convention)`), their refusals and the fingerprint |
| Prediction and the readers | `BassPrediction` (f + u_g at an observed level; a population case's own effect), the retained pointwise content, the training fit |
| A working response that moves each sweep, with its refit | GBASS's path: `WorkingResponse` at a latent state, `BassStructureMoves.ProposeFollow` and `Follow`, the design's `Followed` move, which re-evaluates the columns and the factor when the weights move |
| A Metropolis step on a marginal, by a refit at the proposal | `BassSkewRescaleBlock`: ProposeFollow, the target's evaluation at the chain's temperature, the coefficients redrawn after an accepted move |
| The engine's reach with the block | Chains, ladders with the reference rung, SMC and its flags, the flat histogram, resume: KHAOS's slice 3 gated each against the enumeration |

Refused today, each with its reason in the message: GBASS's noise with the block; the candidate density (so the conformal layer); the collapsed pointwise content; a direct functional response and a decomposition's components.

## The model, as built and as extended

y = α + Σⱼ βⱼBⱼ(x) + Σₖ Z̃ₖuₖ + ε, ε ~ N(0, σ²W⁻¹) (W = I for Gaussian BASS, the noise's 1/(c·vᵢ) under GBASS's, with its offset σβv), with each block k a grouping's intercepts or its slopes on a declared input: Z̃ₖ the cases' level indicators, times the standardized input for a slope. uₖ | σ² ~ N(0, σ²sₖ²Πₖ), Πₖ the projector its convention declares. With Z̃ = [Z̃₁Hₖ…] over orthonormal bases Hₖ of the conventions' subspaces and S = diag(sₖ) by block, y's covariance is σ²V, V = W⁻¹ + Z̃S²Z̃ᵀ. The structure's marginal integrates a, σ² (where the structure step integrates it) and every uₖ by whitening the rows by V^(−½); the effects are drawn from their exact conditional given a and σ². The sₖ are drawn or held.

## Decisions

1. **s drawn by a Metropolis step on ln s on the marginal.** It is a step on the marginal with a, σ² and u integrated, as τ's slice step is, made by a refit at the proposal (the β·v move's pattern, `ProposeFollow` and the target at the chain's temperature). The effects are then drawn last, given the draw's s, as now.
   - *Validity.* An exact Metropolis–Hastings step on the marginal, with no auxiliary variable read by the proposal. Nothing reads u between that step and u's draw, so no partially collapsed ordering is at stake under Gaussian BASS. Under a noise (decision 7) u is drawn right after a, before any block that conditions on it.
   - *The proposal.* t = ln(s′/s) ~ N(0, h²), h = 2.4/√C, C = (L − 1)/2 plus the prior's curvature at s, the information L levels carry about ln s. C reads the state alone, so the step is exact with C's Hastings ratio in A. With several blocks each sₖ takes its own step, C from its own qₖ.
   - *Its place.* τ's slot, after τ's slice step, so tempering replicas and SMC particles run it. Under Gaussian BASS the coefficients and u are drawn at composition, after it; under a noise the coefficients' slot follows τ's, so a and u are drawn after it every sweep and nothing reads a stale draw.
   - *The alternatives.*
     - **The effects as latents, s | u conjugate,** the partially collapsed scheme: (S, a, σ²) on the u-marginal at s, then u | a, σ², s, then s | u. It is valid in that order, but the structure step on the u-marginal still needs V at the new s, so it re-whitens every sweep as the Metropolis step does. It adds the funnel: s | u sticks near 0 when the groups hold few cases. The marginal step dominates it.
     - **The structure moves conditional on u,** u an offset like GBASS's σβv. It needs no whitening at all and leaves the design's caches as they are, but gives up the collapse: u and the level and terms mix by Gibbs, slowly where they are correlated. It is the alternative to measure if slice 1's cost row shows the refit dominating a sweep; the measurement then decides between the two, by effective sizes per second of ln s, u's spread and M.
     - **A slice step on ln s:** some six refits a step against one.
   - *The cost.* One refit at the proposal a sweep: the current structure's columns re-whitened at s′ (O(n + L²) each) and the factor rebuilt (O(nM² + M³)), where held s costs nothing. Slice 1's cost row measures it. The sufficient-statistics route makes a new s cost O(LM² + L²M + M³), free of n: the whitened Gram is XᵀX − (ZᵀX)ᵀK(s)(ZᵀX) with K(s) by Woodbury. It lives in the collapsed core's factor, E34's ground, so it is a follow-up with that row as its decision.
2. **s's prior: half-Cauchy on s, C⁺(0, A),** with A declared and 1 by default, since s is in units of σ: a group effect about as large as the noise (Gelman 2006; Polson and Scott 2012). It is proper, so the reference rung and SMC's start draw s from it. Held s stays the default, so every config written before declares it as before.
   - *Alternatives:* an inverse gamma on s², conjugate only for the latent scheme, whose small shape and scale put mass near 0 that Gelman shows the posterior inherits; a half-normal, with a lighter tail that lets the data move s less when the groups disagree with the prior.
3. **The whitening rides the working response,** not the design. V depends on s and, under a noise, on W, both of which move with the state. So the block's whitening (V^(−½) and ln|V|) is carried by the working response at its state, as the weights are, and the design reads it from the response it is on, through the `Followed` move that already re-evaluates columns when the weights move. Held s and unit weights give the arithmetic of today, bit for bit.
4. **Several groupings: the general whitening by a thin factorization, with one grouping's closed form kept.** With Z̃ (n × q, q = Σₖ qₖ) factored once as QR (rank-revealing, Q n × r), V = I + Q·(RS²Rᵀ)·Qᵀ, so V^(−½) = I + Q·E((1 + Λ)^(−½) − I)Eᵀ·Qᵀ from an r × r eigendecomposition per value of S. That is O(r³) per s and O(nr) a column.
   - One intercept grouping keeps its closed form, O(n + L²) a column, so its fits stay bit for bit; the general form agrees with it to 10⁻¹² (gate).
   - *Alternative:* the per-grouping closed forms composed, which is exact only when the groupings' spans are orthogonal in the case metric, as nested groupings under the within-parent convention are, and wrong for crossed ones.
5. **Each grouping's convention, declared.**
   - *ObservationCentred* (built): Σ_g n_g u_g = 0, the grouping's contribution has mean zero over the cases, and the intercept keeps the level. The default for intercepts, crossed or not.
   - *WithinParent* for a nested grouping, which names its parent grouping: Σ_{g ∈ p} n_g u_g = 0 for each parent level p, so the inner effects are orthogonal to the parent's in the case metric and the parent carries each parent level's mean. A grouping declared nested whose level falls under two parents is refused.
   - Crossed groupings whose levels alias (one level of each holding the same cases) are identified by their variance components alone; the fit's log reports the alias, and no constraint is added.
6. **Random slopes: a declared input, a declared standardization, and a convention against the level and the slope.**
   - The input is a continuous predictor by name, and its scale is declared, with no default: `Unit`, the [0, 1] the input transform gives the terms, or `Standardized`, mean 0 and sd 1 over the training cases. s_b is then in units of σ per unit of that scale.
   - The slopes' convention, *ObservationCentred*: their contribution has mean zero over the cases and zero inner product with the standardized input. So the intercept keeps the level, and the population slope stays with the terms.
   - Slopes are independent of the intercepts' effects at first. A slope correlated with its grouping's intercept, a 2 × 2 covariance per grouping, is a follow-up whose decision is whether recovery from outside the model shows the independent form missing it.
7. **Under GBASS's noises, the block whitens under the weights, and u is drawn each sweep.**
   - *The whitening at W.* For one intercept grouping it keeps its closed form with the counts weighted, n_g^w = Σ_{i∈g} wᵢ, and the indicator vectors scaled by W^½. In general Z̃_w = W^½Z̃ is factored afresh when the weights move, O(nr²) a sweep.
   - *The draw of u.* u is drawn in the coefficients' slot right after a, given a, σ², S and the latents, so v's block, σ²'s, the shape's and interweaving's condition on it.
   - *What each block then reads.* Every block reads the residual less Z̃u. Beyond that:
     - σ²'s block adds u's prior to its conditional: q more in the count and |γ|²/s² in the energy, since u's variance scales with σ² as the coefficients' does.
     - Interweaving adds u's term, as it adds the coefficients' prior.
     - The β·v move becomes a step on the (a, u)-marginal, with a and then u drawn afresh after an accepted move.
     - s's step reads the marginal at the current weights.
8. **The candidate seam: a case in a known group and in a new one,** its shape agreed with the conformal chip before slice 5 opens. The proposal:
   - *A known group:* log p(y | x_t, θ_s) = the training terms' density at f + Σₖ u_{k,g(t)}, conditioned on the draw's effects, as the hierarchical P-spline's new point on a seen curve is. The groups are offered as strata for Mondrian sets.
   - *A new group:* the draw's population predictive with its effects integrated: N(y | f, σ²(1 + Σₖ sₖ²·z̃ₖ²)) under Gaussian BASS. Under a noise, y | v ~ N(f + σβv, σ²(cv + Σ s²z̃²)) is a mixture over v again, read by quadrature on v or refused at first.
   - *Units are cases* for the split form. The full step stays refused, because the whitening couples a level's cases, so leaving a case out moves its group's effects. Groups as the units for new-group coverage (calibration groups held out whole) is the chip's call.
9. **The docs get a page of their own,** `docs/Maths/RandomEffects.md`, by the user's rule for GBASS (2026-11-08: no growing the regression monolith). `BayesianRegression.md`'s KHAOS section keeps its place with one pointer line, edited and committed in one step.

## The slices

| Slice | What | Its gates |
| --- | --- | --- |
| 0 | The whitening onto the working response (decision 3): the block's V^(−½) and ln\|V\| carried by the response at its state, the design reading them through `Followed` | KHAOS's slice 3 gates unchanged. A held-s block fit's draws pinned by hash before the change and met bit for bit after. The pinned tier |
| 1 | s drawn (decisions 1 and 2): `BassRandomIntercept`'s scale prior, the step on ln s in τ's slot, s on the state, the draw and the retained draw (`BassDraw` carries s), prediction's population spread at the draw's s | The step's orbit law on a frozen structure: its exact law in ln s on a grid of the marginal, met by the step's chain. The law by enumeration: every structure's posterior with τ and s integrated by quadrature, on a chain and on each rung of a ladder with its reference rung. SMC's ln p(y) against the enumeration. Geweke with s drawn: ln s, the effects' spread, M, ln τ, ln σ². The flat histogram against the enumeration, as KHAOS's slice 3 gated it, or refused with its reason if its weights read a marginal that s moves. Recovery (below). Resume. The cost row: a sweep with s drawn against held |
| 2 | Several groupings, crossed and nested (decisions 4 and 5): `BassData`'s groupings, the config's list of blocks with the one-intercept form kept, the general whitening, the conventions | The general whitening against dense algebra, and against the closed form to 10⁻¹² for one grouping. Each convention's identifiability check. The law by enumeration and Geweke with two crossed groupings and with a nested pair, s's drawn. Recovery. Resume. The cost row by r |
| 3 | Random slopes (decision 6) | Dense algebra. The slopes' convention. The standardization's invariance: a fit with x rescaled by c under the declared scale draws as the fit on x, bit for bit where the arithmetic allows and to 10⁻¹² otherwise. The law by enumeration and Geweke. Recovery of slopes. Resume |
| 4 | The block under GBASS's noises (decision 7) | Geweke with the block under t, the asymmetric Laplace, the Normal-Wald, the skewed horseshoe with its β·v move, and interweaving. The whitening at W against dense algebra. The β·v move's orbit law on a frozen state with the block. Recovery under heavy-tailed noise. Resume. Every noise fit without a block bit for bit |
| 5 | The candidate seam (decision 8), its shape agreed with the conformal chip | At the training cases the candidate is the pointwise column, bit for bit. The new group's density against the population predictive's closed form, and under a noise against quadrature. Coverage by simulation: within known groups, marginally and by group stratum, and for new groups, split over groups. The refusals with their reasons |

Slices 1 to 3 run under Gaussian BASS; slice 4 carries every block to the noises; slice 5 reads all of them.

## Gates common to every slice

- **Recovery from outside the model.** Effects (and slopes) generated from laws the block's centred normal does not hold, such as a skewed or a two-point law, on unbalanced groups, with s unknown. The gates are: the effects' 90% intervals and their standardized errors' mean square, pooled over the levels as KHAOS's slice 3 gated them; s's interval; σ's; the RMSE at held-out cases of observed levels; the population surface no worse than without the block, by at most 10%.
- **The law, by enumeration and by Geweke,** with the drawn scales among the functions.
- **Each convention's identifiability check:** its constraint in every draw to 10⁻¹², and the fitted level (and slope) where the convention puts it, within 2 posterior sd of the generator's.
- **Resume, bit for bit,** the scales and the effects included.
- **Fits without a block, bit for bit:** the pinned tier at each commit (`scripts/test-select.ps1 -Scope Commit`), the whole assembly once at the slice's end (`-Scope Slice`), and the qualification gates behind `BASS_QUALIFICATION`.
- **A fresh-eyes review** of each slice before its report, fixed before reporting.

## Coordination

- **E34's chip** works in the collapsed core and `BassData`. Its decision 7 puts random intercepts with one shared s under K outputs, with row covariance V ⊗ Σ, so decision 3's whitening on the working response is the form it reuses. Before slice 0's edit of `BassDesign` and slice 2's of `BassData` I send it the change's shape, check `git diff -- <file>`, and commit each edit by itself.
- **The conformal chip** owns the layer and `posterior/BassPrediction.cs`. Slice 5's shape is agreed with it first (decision 8), and any edit to its files waits for its tree to be clean.
- **The KHAOS chip** takes E4's proposal kernels and leaves `BassRandomIntercepts.cs` to this brief.
- **Shared docs** (`BayesianRegression.md`, the notation ledger): checked, edited and committed in one step. The ledger takes sₖ, Z̃, Πₖ and the slopes' scale beside KHAOS's s.

## Not in this brief

- A slope correlated with its intercept (decision 6's follow-up).
- The sufficient-statistics route for a new s (decision 1's follow-up, E34's ground).
- Random effects under a direct functional response or a decomposition, and E34's K-output form, which its own brief builds on decision 3.
- BARS's and the other engines' blocks.

## Working rules

The GBASS chip's rules carry over:
- Work on `main`. Commit by pathspec, adding new files first, and leave nothing staged. Push nothing, and put no attribution lines in commits.
- Builds go to `artifacts\chip-ridge-blocks\`, which is deleted at the end. Heavy runs run below normal priority with capped threads, one heavy class at a time.
- Tests run through `scripts/test-select.ps1`: the Commit scope per commit and the Slice scope once at a slice's end, from a snapshot of HEAD plus this chip's files while other sessions hold edits.
- No subagents without permission.
- Each slice's report goes to the main session, and the next slice waits for its verification note.

## Report

Each slice reports here and to the main session:
- what landed, with commits;
- each gate's reading;
- the decisions taken beyond this brief;
- what was found;
- the follow-ups, each with the decision it would change.

## Reports

Newest first.
