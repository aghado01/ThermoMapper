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
4. **Several groupings: the general whitening by a thin factorization, with one grouping's closed form kept.** With Z̃ (n × q, q = Σₖ qₖ) factored once into Q·R (rank-revealing, Q n × r), V = I + Q·(RS²Rᵀ)·Qᵀ, so V^(−½) = I + Q·E((1 + Λ)^(−½) − I)Eᵀ·Qᵀ. The scale-free factor uses the existing thin SVD. A reduced SVD of R·S obtains E and Λ at each scale vector without explicitly forming the covariance Gram: the latter lost the small component at scales (10⁶, 0.6) in an independent two-contrast reference. Work is O(qr²+r³) per scale vector and O(nr+r²) per column, independent of n for the scale factorization. Effects use the reduced factor's right vectors in their conditional gain and retain all q latent coordinates under aliases.
   - One intercept grouping keeps its closed form, O(n + L²) a column, so its fits stay bit for bit; the general form agrees with it to 10⁻¹² (gate).
   - *Alternative:* the per-grouping closed forms composed, which is exact only when the groupings' spans are orthogonal in the case metric, as nested groupings under the within-parent convention are, and wrong for crossed ones.
5. **Each grouping's convention, declared.**
   - **Public API decision, 2026-11-09.** The user chose “Use named collections and migrate the singular API.” Slice 2 replaces `BassData.Grouping` and `BassConfig.RandomIntercept` with named groupings and intercept blocks, and moves each draw's `GroupScale` and `GroupEffects` into named per-group records. A single group is a one-entry collection and uses the existing arithmetic and random draw order. Migrate callers and checked-in examples and show the one-entry conversion for saved configs.
   - **Run identity.** Config JSON contributes to the [run fingerprint](../../../src/maths/regression/spline/bass/Bass.cs). The renamed grouped configuration changes that fingerprint. Under the existing resume check, converting a config alone cannot resume an old grouped checkpoint. Document this when reporting the migration. Check that ungrouped fingerprints and the existing one-group draw hashes still hold, and that a checkpoint written with the new representation resumes.
   - *ObservationCentred* (built): Σ_g n_g u_g = 0, the grouping's contribution has mean zero over the cases, and the intercept keeps the level. The default for intercepts, crossed or not.
   - *WithinParent* for a nested grouping, which names its parent grouping: Σ_{g ∈ p} n_g u_g = 0 for each parent level p, so the inner effects are orthogonal to the parent's in the case metric and the parent carries each parent level's mean. A grouping declared nested whose level falls under two parents is refused.
   - Crossed groupings whose levels alias (one level of each holding the same cases) are identified by their variance components alone; the fit's log reports the alias, and no constraint is added.
6. **Random slopes: a declared input, a declared standardization, and a convention against the level and the slope.**
   - The input is a continuous predictor by name, and its scale is declared, with no default: `Unit`, the [0, 1] the input transform gives the terms, or `Standardized`, mean 0 and sd 1 over the training cases. s_b is then in units of σ per unit of that scale.
   - The slopes' convention, *ObservationCentred*: their contribution has mean zero over the cases and zero inner product with the standardized input. So the intercept keeps the level, and the population slope stays with the terms.
   - Slopes are independent of the intercepts' effects at first. A slope correlated with its grouping's intercept, a 2 × 2 covariance per grouping, is a follow-up whose decision is whether recovery from outside the model shows the independent form missing it.
   - **Rescaling contract (user decision, 2026-11-09).** Rescale the continuous input and its domain together for both `Unit` and `Standardized`. At ×4, require bit-identical full sampled-scale chains on the fixture. This power-of-two control preserves the transformed inputs exactly. Investigate any mismatch as a scale-dependence defect. At ×3.7, compare fixed-state whitening, marginals, conditionals and held-scale fits to 10⁻¹² relative to their magnitude. For full sampled-scale chains at ×3.7, require identical structure paths and every accept/reject decision, with continuous draws within 10⁻⁹ relative. Decimal rescaling can change the last bit of the transformed inputs, and the adaptive updates propagate those differences. Record the relative comparison scale and its treatment of zero alongside each gate.
7. **Under GBASS's noises, the block whitens under the weights, and u is drawn each sweep.**
   - *The whitening at W.* For one intercept grouping it keeps its closed form with the counts weighted, n_g^w = Σ_{i∈g} wᵢ, and the indicator vectors scaled by W^½. In general Z̃_w = W^½Z̃ is factored afresh when the weights move, O(nr²) a sweep.
   - *The draw of u.* u is drawn in the coefficients' slot right after a, given a, σ², S and the latents, so v's block, σ²'s, the shape's and interweaving's condition on it.
   - *What each block then reads.* Every block reads the residual less Z̃u. Beyond that:
     - σ²'s block adds u's prior to its conditional: q more in the count and |γ|²/s² in the energy, since u's variance scales with σ² as the coefficients' does.
     - Interweaving adds u's term, as it adds the coefficients' prior.
     - The β·v move becomes a step on the (a, u)-marginal, with a and then u drawn afresh after an accepted move.
     - s's step reads the marginal at the current weights.
   - **Amended by the user on 2026-11-10 ("B with your refinement"), after slice 4's finding.**
     - *The finding.* Integrating u out and tempering that likelihood is exact only at b = 1. Under a noise, the blocks for v, σ², the shape and interweaving condition on u, and their conditionals at b < 1 belong to tempering the likelihood given u.
     - *The form.* Under a noise, every rung carries u as a latent that the structure reads as an offset Z̃u, with no whitening. That makes ladders, SMC and the reference rung exact, and s steps on N(u; 0, σ²s²) with no refit.
     - *At b = 1 only* (plain chains, the cold rung, SMC's final stage), the structure step reads the u-integrated likelihood at the current weights, and u is drawn from its conditional right after, before any block reads it. This is the partially collapsed order of the β·v fix, so each rung's kernel leaves its own target invariant, and swaps read the shared state (structure, latents, u).
     - *The evidence is unchanged:* both tempering paths end at p(y).
     - *The Gaussian collapse stays as built.*
     - *Slice 4's finding 2:* under weights, the centring convention uses the GLS intercept direction c = L·1, not unweighted counts.
   - **Amended again by the user on 2026-11-10 ("Make it opt-in"), on slice 4's measurement.**
     - *What was measured.* The b = 1 collapse reaches the same posterior but runs at 0.14 to 0.71 of plain form B's ESS per second.
     - *The default.* Plain form B at every rung.
     - *The collapse* becomes a sampling option, off by default. Its gain is expected where u and the basis compete for the same variation, so it is measured on a confounded dataset (a grouping aligned with an input, or a few large groups) before any default changes.
     - *Still with the user:* the u–s funnel's fix (ln s's integrated autocorrelation time is about 19,000 sweeps), by a step on s with u integrated, exact at b = 1 only, or by interweaving u's centred and non-centred forms, exact at every rung.
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
| 2 | Several groupings, crossed and nested (decisions 4 and 5): named collections for the data, config and retained effects, the existing one-group calculation, the general whitening and the conventions | The general whitening against dense algebra, and against the closed form to 10⁻¹² for one grouping. Each convention's identifiability check. The law by enumeration and Geweke with two crossed groupings and with a nested pair, s's drawn. Recovery. Resume. The cost row by r |
| 3 | Random slopes (decision 6) | Dense algebra. The slopes' convention. Decision 6's three rescaling gates on both declared scales (×4 bit identity, ×3.7 primitives and held fits at 10⁻¹² relative, ×3.7 sampled chains with identical decisions and continuous values at 10⁻⁹ relative). The law by enumeration and Geweke. Recovery of slopes. Resume |
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

## Slices 4 and 5, commissioned (2026-11-10)

Commissioned by the user on 2026-11-10: the GBASS chip takes slices 4 and 5, which the GPT workers that built slices 1 to 3 left.
- **Read first** slice 2's [report](ridge-multiple-groupings-report.md) and slice 3's [report](ridge-random-slopes-report.md). They landed after this chip's slice 0, so read the current code before editing.
- **Order.** Slice 4, then slice 5, each reported to the main session and verified before the next starts.
- **Who else is writing:**
  - *The E34 chip* builds its slice 5, the PCA mode's joint option, now. Its slice 4, noise by row and random intercepts under K outputs, waits for this slice 4's verification, because it reuses the block under the noises.
  - *The conformal chip* holds `ConformalBayes.cs` and its tests through slices 12a, 10, 11 and 12b, and owns `posterior/BassPrediction.cs`.
- **Shared files** (`BassConfig.cs`, `BassModel.cs`, `BassResolved.cs`, `BassSweep.cs` and the like):
  - message the other chip with the change's shape before editing;
  - check `git diff -- <file>` for changes you did not make;
  - commit each such edit by itself.
- **Slice 5's seam** is agreed with the conformal chip by message first (decision 8). Edits to the conformal chip's files wait until its tree is clean.
- **Docs.** Chips leave `docs/` alone (AGENTS.md, the user's standing guidance of 2026-11-10). This replaces the "Shared docs" line under Coordination: list the doc changes owed in the report. Notation-ledger entries stay in `issues/`, each edited and committed in one step.
- **Messages** are terse: events and coordination, with a pointer to the report. No session reads another's transcript.

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

### In progress: slice 5, the candidate seam under random effects (GBASS chip, updated 2026-11-10)

Hand-over state, updated in place until the slice's report replaces it.
- *State.* The shape below went to the conformal chip on 2026-11-10. That chip is paused for the user's quota. Nothing is built until the shape is agreed.
- *Proposed shape* (decision 8, refined against the code):
  1. *`BassFit.CandidateLogDensity` takes random-effect fits*, which it refuses today. Each candidate names its groups as prediction does (`BassPredictionInputs.GroupsOf`).
     - *A level seen in training* reads the draw's effects. Under Gaussian noise that is N(y | f + Σₖ z̃ₖ·u_{k,g}, σ²), with z̃ the layout's multiplier: 1 for an intercept, the scaled input for a slope. Under a noise it is the law's density with v integrated. These are the training pointwise terms (`BassPointwise.Retained` with `BassDrawFit` adding Z̃u), so a training case's candidate is its column bit for bit.
     - *A new group* takes that grouping's effect integrated from its prior. That covers a level unseen in training under `UnseenLevel.OutsideEverySubset`, or the grouping's column left out. Under Gaussian noise it is N(y | f + the known effects, σ²(1 + Σ_new s_k²·z̃_k²)). Under a noise it is a mixture over v again, y | v ~ N(f + known + σβv, σ²(c·v + Σ_new s_k²·z̃_k²)). That density is taken by quadrature on v when `quadrature: true` is passed, and refused otherwise. It needs a law density with an added variance beside `IBassNoiseLaw.LogDensity`, in `BassNoise.cs`, a shared file. Its CDF needs the same, on `cdf: true`.
  2. *The latent report:* the group effects are conditioned for known groups and integrated for new ones. Under a noise, v is also reported as integrated. `BassPointwise.Latents` drops that line today when random effects are present, a slice-4 gap this slice fixes.
  3. *Open, the conformal chip's call:* may one call mix known-group and new-group cases (with a per-case flag on `CandidateLogDensity`), or is a mixed call refused, so each call is read by one density (CB-4)?
  4. *Open:* does BASS give a helper with each training unit's and each candidate's level of grouping k, as strata for Mondrian sets, or does the conformal layer build strata from the data?
  5. *Units* are cases, for the split step. The full step and Hoff's leave-one-out score stay refused for random-effect fits, because u couples a level's cases, so leaving a case out moves its group's draw. Coverage for new groups, with whole groups held out, is the conformal chip's design.
- *Files:* `BassCandidates.cs`, `BassFit.cs` (the candidate entry), `BassPointwise.cs` (the latent report), `BassNoise.cs` (the widened law density, shared). Ownership was asked of the conformal chip.
- *Gates, per the slices table:*
  - at the training cases the candidate is the pointwise column, bit for bit;
  - the new group's density against the population predictive's closed form, and under a noise against independent quadrature;
  - coverage by simulation within known groups, marginally and by group stratum, and for new groups split over groups;
  - the refusals with their reasons.

### Slice 4: random effects under GBASS's noises (GBASS chip, 2026-11-10)

**Landed** in `4f5bf3f` (the block under the noises, its tests), `94133f8` (the notation ledger), `1940cc0` (the exact-start Geweke rows, committed as WIP by the main session at the user's word), `397ca88` (the mixing measurement lengthened to 40,000 sweeps per chain), and the commit that carries this report. Decision 7 as the user amended it on 2026-11-10 ("B, with the refinement").

**The form.**
- *u is a latent on every rung.* Under a noise the reading carries u as the offset Z̃u beside σβv (`RandomEffectsLatent`, `WorkingResponse.Effects`). σ² stays on the record, as a skew already keeps it (`VarianceMeanMixture.RecordsSigma2`), because u's prior scales with σ². Every rung tempers the likelihood given u and the latents, so ladders, SMC and the reference rung are exact.
- *u's block* (`BassGroupEffectsBlock`, right after the coefficients) draws u from N(K·Z̃ᵀ(bW)e, σ²K) at the rung's power, with K = (S⁻² + Z̃ᵀ(bW)Z̃)⁻¹, through the weighted whitening at bW. The sweep's fit is handed on with Z̃u in it, so v's, σ²'s, the shape's and interweaving's blocks read the residual net of u. σ²'s and interweaving's priors add u's energy Σ‖u_k‖²/s_k² on Σq_k coordinates.
- *s's step* (`BassLatentGroupScaleBlock`, in τ's slot) is a Metropolis step on ln s_k against its half-Cauchy and N(u_k; 0, σ²s_k²), with no refit and h = 2.4/√(2q_k + the prior's curvature).
- *At b = 1 only* (`BassCollapsedStructureBlock`), the chain follows the u-integrated reading at the current weights (`VarianceMeanMixture.Integrated`). It takes one reversible-jump move on that marginal, draws a and then u, and follows the reading given u again. A plain chain takes this step in place of its move step, and a replica or a particle takes it after its own moves.
- *Retained draws* keep their latents, their u and their scales (`BassDraw.Groups`), and a draw's reading is rebuilt with its u. Prediction now adds a population case's own effect under a noise law too, and the mean surface's E[u | ·] reads the draw's weights.

**Found.**
1. *The collapse alone is exact only at b = 1* (reported before the amendment). The core tempers a state's evaluation analytically in b, with u integrated out of the likelihood. The noise's blocks condition on u, and at b < 1 their conditionals belong to tempering the likelihood given u. Hence form B, with the collapse kept at b = 1.
2. *Under weights, the GLS intercept direction is c = L·1, not W^½·1.* The conventions centre u with unweighted counts, so L = V_w^(−½)W^½ moves the constant. The design's weighted path takes L and centres against c, with the GLS mean and mass cᵀc (`WeightedWhitening`); at W = I this is the unweighted block's. On the dense checks, c departs from W^½·1 by 0.27 to 0.68.
3. *A stale prior Gram across readings* (fixed). A design keeps its prior Gram when it follows a reweighting, which holds only while the prior's side does not change with the reading. The integrated reading's prior first read V^(−½)z and the reading given u read the raw z, so a chain moving between them carried a Gram that depended on its history. A resumed chain parted from the continuing one at its first draw. Under B the coefficient prior belongs to the model given u, so both readings now read the raw columns, and resume is bit for bit.
4. *The u–s funnel under form B.* With u a latent, s moves only through u, and where the data pin u, at large s, s barely moves. Geweke's successive-conditional chain, which redraws y each step, has to cross that funnel; a fit on fixed data stays near where its data put u. See the Geweke gate below. A step on s with u integrated (at b = 1), or an interweaving move on (s, u), would address it. Both are design changes, so neither is built.
5. *The design's layout rule* refused a reading that gains or loses a whitening. Under weights every column is evaluated afresh on a reweighting, so a weighted design now passes between the two.
6. *The column-mean cache* was keyed by the weight array alone. Under a weighted whitening the GLS centring moves with the scales too, so the cache is now keyed by the whitening.

**Gates.**
- *Dense algebra:* with and without a slope, at b = 1 and 0.3: xᵀC⁻¹y, the mass, the mean, ln|V_w| and u's conditional mean agree with dense algebra to 1.4·10⁻¹⁴ at worst. c departs from W^½·1 by 0.27 to 0.68.
- *The design on the integrated reading is the GLS fit:* on three terms the residual against a dense GLS fit is ‖r‖² = 1.1·10⁻¹⁵, and the mean, the mass and Σ ln w − ln|V_w| agree to 7·10⁻¹⁵.
- *The sweep's order and a fit:* a chain runs the collapsed moves, τ, the group scale, the coefficients, u, the latent scales, then σ². A replica runs the same order after its own moves. 42 of 300 collapsed moves were accepted in a fit.
- *Resume, bit for bit:* a tempered fit with a drawn scale and the reference rung. A probe covered plain chains, ladders with and without the reference rung, and held and drawn scales.
- *Noise fits without a block, bit for bit:* three pins taken on HEAD before the change (t `5471c4bcad7d0b09`, skewed horseshoe `afd1aa3f399ea35f`, tempered Normal-Wald `98b9820e14563b71`).
- *Geweke, first design (successive-conditional), unqualified and kept on record.* Each row ran 300,000 marginal-conditional draws against 450,000 steps of one chain. The z for each function was the mean difference over √(variance/300,000 + the square of the 50-batch se). The functions were M, ln τ, ln σ², mean ln v, β, ln s, u0, Σu², atan f0, atan f6 and z0², with bar |z| < 4.

  | row | worst \|z\| | ln s: draws, chain, z | other \|z\| ≥ 3 |
  |---|---|---|---|
  | t, b = 1, collapse | 1.94 (pass) | −0.0001, −0.031, 0.52 | — |
  | t + interweave, b = 1, collapse | **4.51** | −0.0001, −0.127, 4.51 | — |
  | quantile, b = 1, collapse | 3.53 (pass) | −0.0005, −0.081, 3.53 | — |
  | Normal-Wald, b = 1, collapse | **5.12** | −0.0000, −0.149, 5.12 | ln σ² 3.33, Σu² 3.14 |
  | t, b = 1, plain | **4.30** | −0.0001, −0.130, 4.30 | — |
  | t, b = 0.5, replica | **4.02** | −0.0020, −0.111, 4.02 | — |

  The Normal-Wald at 0.5 and both horseshoe rows were stopped unrun under this design. Every chain's mean of ln s fell below the prior's 0. Σu² and u0 have no finite variance under the half-Cauchy s, so their z values cannot be read.
- *Why, on a probe of the t plain row's chain run four times as long* (1.8 million steps, another start):
  - ln s's autocorrelation is 0.965 at lag 1, 0.48 at lag 100, 0.28 at lag 10,000 and 0.04 at lag 50,000.
  - The batch-means IAT keeps growing with the batch: 444 at 1,000, 3,210 at 9,000 (the test's batches), 8,523 at 30,000, and 19,374 and 18,547 at 90,000 and 180,000. The ESS over 1.8 million steps is about 95, so the test's se was too small by a factor of about 2.4.
  - z shrinks as the chain grows: 3.09, 1.77 and 1.09 after 450,000, 900,000 and 1.8 million steps on 50 batches (1.68, 0.97 and 0.78 on 10). The mean falls from 0.83 to 0.35 to 0.14: the sign follows the start, and the size follows the chain's length.
  - The chain spends 0.150 of its steps at s > 5, against the prior's 0.126, in 22,584 excursions, the longest 63,563 steps. At s < 0.2 it spends 0.122 against 0.126, the longest excursion 205 steps. The large-s tail is the funnel's.
- *Geweke, exact starts (the gate):* 10 sweeps at b from each of 40,000 independent draws of the tilted joint. Each function's change has mean 0 under invariance, and the z reads the paired changes, so mixing cannot bias it. The functions are the first design's, with u0 read as atan u0 and Σu² as ln Σu² and as u's energy Σu²/(σ²s²). All nine rows pass |z| < 4:

  | row | worst \|z\| (function) | ln s z |
  |---|---|---|
  | t, b = 1, collapse | 0.91 | −0.06 |
  | t, b = 1, plain | 1.98 (atan f0) | 0.83 |
  | t, b = 0.5, replica | 1.32 | −1.07 |
  | quantile, b = 1, collapse | 2.58 (atan f0) | −0.28 |
  | Normal-Wald, b = 1, collapse | 1.13 | −0.94 |
  | Normal-Wald, b = 0.5, replica | 2.64 (u energy) | 0.51 |
  | skewed horseshoe + β·v move, b = 1, collapse | 2.14 (β) | −1.33 |
  | skewed horseshoe + β·v move, b = 0.75, replica | 3.17 (ln Σu²) | −2.05 |
  | t + interweave, b = 1, collapse | 2.81 (atan f0, f6 and u0 together: one effect) | 0.05 |

  The two horseshoe rows also read the β·v orbit with the block.
- *The exact-start gate's power, on a probe.* s's step was replaced by a mutant whose target carries an extra s^(−ε); under the half-Cauchy this moves E[ln s] by −ε·π²/4. On the t plain row:
  - At ε = 0.053, the −0.13 shift behind the first design's z 5.1, ln s's change is −0.0270 ± 0.0054, z −5.00.
  - At ε = 0.02 (−0.049) it is −0.0074 ± 0.0054, z −1.36.
  - Unmutated it is 0.0043 ± 0.0054, z 0.80.

  Ten sweeps carry about a fifth of a stationary shift. So the gate sees the failure's size, and misses a bias of a third of it.
- *A ladder and SMC end to end* (n = 120, t₄, twelve unbalanced levels):
  - ln p(y): ladder −114.98 ± 1.08, SMC −114.41 ± 0.09.
  - ln s: ladder 0.956 ± 0.007 and SMC 0.960 ± 0.003, against chains' 0.952 ± 0.002.
  - M: ladder 4.218 ± 0.021 and SMC 4.189 ± 0.003, against chains' 4.201 ± 0.016.
- *Recovery under heavy-tailed noise* (t₃, effects generated outside the model, four sets): f + u RMSE 0.098, 0.107, 0.141 and 0.098 with the block, against 0.365, 0.643, 1.086 and 0.493 without. The effects' 90% intervals cover 46 of 48, and the standardized errors' mean square is 0.76.
- *Mixing, the collapse against form B without it:* n = 300, t₄, twelve unbalanced levels, with four chains of 2,000 + 40,000 sweeps per arm on two seeds, plus a probe on other chain seeds.
  - *Both arms reach the same posterior.* M reads 4.365 ± 0.022 collapsed against 4.299 ± 0.020 plain (seed 11), 4.307 ± 0.048 against 4.262 ± 0.017 (seed 23), and 4.290 ± 0.026 against 4.303 ± 0.051 (probe). ln s and ln σ² agree within 0.003.
  - *Bulk ESS per sweep, collapsed over plain:* ln σ² 1.17, 1.15 and 1.60; ln s 1.02, 1.02 and 1.04; u0 0.75 and 1.06; M 0.55, 0.61 and 1.20.
  - *Time per sweep, collapsed over plain:* 1.76, 4.40 and 2.24, on the cycle clock with the machine under load.
  - *ESS per second, collapsed over plain:* ln σ² 0.67 and 0.26, ln s 0.58 and 0.23, u0 0.43 and 0.24, M 0.31 and 0.14.
  - The first reading kept 8,000 sweeps per chain and gave M ESS 28 against 90 to 99. That was too short for M, whose IAT is about a thousand sweeps, so the test now keeps 40,000.
- *Per commit:* 2230 of 2230 (`4f5bf3f`'s files on its parent). *Slice end:* 3095 of 3095, on a snapshot of `397ca88`.

**For the user: keep or drop the b = 1 collapse?** Both arms reach the same posterior, so this is a question of cost. On this data, a grouping unrelated to the inputs, the collapse gives:
- ln σ² 15 to 60% more ESS per sweep;
- nothing on ln s;
- an unsettled change on M, 0.55 to 1.2 per sweep;
- 1.8 to 4.4 times the time per sweep.

Net, it reaches 0.14 to 0.71 of form B's ESS per second.
- *Keep* (decision 7's amendment as it stands): this keeps the structure step on the u-marginal at b = 1. Its gain is expected where u and the basis compete for the same variation, such as a grouping aligned with an input or a few large groups; that case is not measured here. It costs this slowdown on every fit with random effects under a noise.
- *Drop:* plain chains take form B's move step, and a ladder's cold rung and SMC's final stage lose their extra step. The untested gain is lost, and fits run at the plain sweep's cost. The step stays reachable through `BassSweepSlots.Blocks(…, collapse: true)` until a cleanup.
- *Opt-in,* the middle way: a sampling option, off by default, measured on a confounded dataset before any default changes.

**The user's call, 2026-11-10: opt-in** (decision 7). Landed in `aa0d0ee`: `BassSampling.CollapseEffectsInMoves`, off by default. Every chain, ladder replica and SMC particle builds its sweep from `BassSweepSlots.Blocks`, which takes the collapsed step only when asked. Without the option every sampler's sweep is form B's, block by block. A chain, a ladder and SMC fit with the option stated off are their default fits bit for bit, and with it on they differ. The tempered resume runs bit for bit both ways. The option is refused under Gaussian noise and without random effects, and is left out of the serialized form when off. Per commit: 2280 of 2280.

**Decisions taken beyond this brief.**
- The Geweke gate reads exact starts, not a successive-conditional chain, for the reasons and with the power above. u0 and Σu² are read bounded and in logs.
- The coefficient prior reads the raw columns under a noise, on both readings (finding 3).
- σ² is kept on the record under random effects without a skew, and an improper σ² prior is refused with random effects, since u's prior scales with σ².

**Coordination.**
- E34's chip had the change's shape before the edits.
- The case-weights chip (`local_57bad71b`) committed its `BassModel.cs` hunk by private index (`7343bc8`) and withdrew its `BassSweep.cs` line.
- E8's chip edits none of these files.
- `posterior/BassPrediction.cs` is the conformal chip's. Its tree was clean, the edit is two noise-path lines, and the chip was told.

**Doc changes owed** (chips leave `docs/` alone):
- `docs/Maths/RandomEffects.md`, a new section "Under GBASS's noises":
  - form B and the b = 1 collapse;
  - u's block, and s's step on u's prior;
  - σ² on the record;
  - the GLS intercept direction c = L·1, with the dense checks;
  - the prior on raw columns under a noise, against the whitened design under Gaussian noise;
  - the readers;
  - the u–s funnel;
  - the exact-start Geweke design and its power;
  - the mixing measurement.
- `docs/Maths/VarianceMeanMixtures.md`: one line in the status list saying the noises take random effects, pointing to RandomEffects.md.
- `docs/Maths/regression/…` (the split BayesianRegression pages): the GBMARS row's "random effects under the noises".

**Follow-ups, with the decision each would change.**
- *The u–s funnel* (finding 4): a step on s with u integrated at b = 1, or an interweaving move on (s, u). It would change whether s mixes where the data pin u at large s. Fits on fixed data mixed ln s well in the measurement (ESS about 30,000 from 160,000 sweeps), so the funnel is the Geweke chain's and weak data's. A design change for the user.
- *Case weights with random effects*: `BassData.cs` refuses them because the block's whitening read precision weights alone. `WeightedWhitening` at precision weights is that GLS, so the refusal could lift for precision weights. That is the weighted-fits chip's call.
- *The collapse's cost*: on the integrated reading every column is evaluated afresh, with tables off under a whitening. That is likely most of the 2.2× per sweep; it was not measured apart. If the collapse is kept, tables on weighted whitenings would change its ESS per second.

### Slice 2 — named crossed and nested Gaussian intercepts (2026-11-09)

The scoped commit `feat(bass): support named crossed and nested intercept blocks` carries the named data, configuration and retained-effect collections, both centring conventions, joint whitening and per-group scale steps. The [report](ridge-multiple-groupings-report.md) records the E34 reconciliation, API and checkpoint migration, dense references, nonlinear recovery and rank-cost rows. The independent unequal-scale review led to the reduced weighted-factor SVD and stable gain in decision 4. Combined qualification passed 66/66, including deterministic covariance and actual draws at scales differing by 10⁶, one-group hash pins, grouped resume and E34 output-rate resume. Commit passed 2,037/2,037 and full Slice passed 2,874/2,874, with no skips. The current-main test assembly also compiled with zero errors. [Receipts](../../../artifacts/e21-ridge-s2-qualification/) ([tracked inputs](../receipts/e21-ridge-s2-qualification/)) retain the baseline, source hashes, settings and TRXs.

### Slice 1 — drawn scale on one Gaussian intercept grouping (2026-11-09)

The scoped commit `feat(bass): draw the Gaussian random-intercept scale on the collapsed target` carries the implementation, tests and [user documentation](../../../docs/Maths/RandomEffects.md). Verification used a snapshot of `541b772` plus this slice's files. The dedicated [receipts](../../../artifacts/e21-ridge-s1-qualification/) ([tracked inputs](../receipts/e21-ridge-s1-qualification/)) include the run settings, baseline and saved TRXs.

- **Sampler and restart.** The optional half-Cauchy prior draws the scale in log coordinates. The shared target carries its density and Jacobian for moves, diagnostics and powered swaps. The proposal rebuilds the whitening from the current design's raw response and reference, including after checkpoint restoration and replica hand-over. Ordinary and tempered fits resumed bit for bit on separately constructed data with identical values, including s and the effects. All three held-scale fingerprints matched.
- **Quadrature.** Frozen two-term hinge and Legendre structures used trapezoids on [−24, 24] at spacing 0.05 and [−28, 28] at 0.025. The largest changes in the first and second log moments were 2.30 × 10⁻⁹ and 6.13 × 10⁻⁸. At 60,000 scale moves, the two moment z scores were −1.46 and 1.04 for hinges, 0.20 and 0.02 for Legendre. The moves accepted 46.8% and 47.0%.
- **Structure and scale law.** The oracle enumerates every zero-, one- and two-term candidate and checks its support at each s. It normalizes by the integrated prior mass on that support. Extending the log-scale grid from [−16, 16] to [−22, 22] and halving the step on both axes moved log evidence, log means and individual structure probabilities by at most 5.08 × 10⁻⁶. The unit half-Cauchy's reference E[ln s] is zero. Observed reference means were −0.0075 ± 0.0054 and −0.0001 ± 0.0053. All structure probabilities, count probabilities, E[ln τ] and E[ln s] passed at b = 1, 0.4 and 0. Four ordinary chains per basis also passed after tempered burn-in, using 200,000 production sweeps each with no swaps. The comparison retained the original batch-error tolerances.
- **Evidence and prediction.** Four SMC replicates of 2,000 particles gave log evidence −17.0414 ± 0.0120 against −17.0540 for hinges and −71.1691 ± 0.0257 against −71.1340 for Legendre. Both returned zero SMC flags. Population predictions use σ²(1+s²) at the selected draw. The config refuses a drawn s with flat-histogram sampling because its weights and WHAM reconstruction have not been qualified on this moving marginal.
- **Joint law and recovery.** Each of two Geweke seeds used 24,000 independent joint draws paired with four posterior sweeps at simulated y. The worst |z| was 2.18 across ln s, log RMS effect, M, ln τ and ln σ². The four skewed, unbalanced recovery datasets had twelve levels each. Both bases' 90% intervals covered 41/48 effects, 3/4 generating scales and 4/4 generating noise deviations. Mean squared standardized effect errors were 1.158 and 1.142. The largest |Σ n_g u_g| was 7.8 × 10⁻¹⁴. Every fitted level lay within two posterior standard deviations of its generating value. Grouped holdout RMSE was less than half the no-block fit's on every dataset, and population surface RMSE stayed within the 10% bound.
- **Cost.** Medians over three 20,000-sweep runs were 4.02 µs held and 10.46 µs drawn on the 18-case hinge fixture, 9.04 µs held and 18.72 µs drawn on the 40-case Legendre fixture. The draft's proposal width uses the grouping's Gaussian information, 2 Σ (s²λ/(1+s²λ))², plus the log-prior curvature with a floor of one. It includes the Hastings correction. Larger n, M and L would determine whether a sufficient-statistics refit improves effective samples per second enough to justify changing the factor.
- **Cold-start finding.** At 20,000 sweeps four ordinary chains begun from the structure prior missed a hinge structure's probability (0.3427 against 0.4650). At 200,000 sweeps a rare pair had pooled occupancy 0.2481 against 0.0002. The fourth chain, with move/block seeds 303/403, occupied the pair `x0−k2 · x1+k0 | x0−k2 · x1−k2` at 418 of 420 observations spaced 500 sweeps apart. On all four chains the incremental and freshly refitted log marginals agreed to 10⁻⁷ at those observations. Tempered burn-in resolved the fixture for the subsequent ordinary-chain comparison. Reliable convergence from a structure-prior start on this fixture remains open.

**Regression receipts.** Release builds ran below normal priority with `MSBUILDDISABLENODEREUSE=1` and `DOTNET_PROCESSOR_COUNT=2`. The focused filter `FullyQualifiedName~RandomEffectsTests|FullyQualifiedName~RandomInterceptTests` passed 39/39. The four opt-in rows in that count returned early. With `BASS_QUALIFICATION=1`, the Geweke, cost and `RandomInterceptQualificationTests` filter ran all eight bodies and passed 8/8. With the flag unset, `scripts/test-select.ps1 -Scope Commit` passed 1,869/1,869 and `-Scope Slice` passed 2,772/2,772. Both used `-SourceRoot` for the snapshot, a private `-ArtifactsRoot` and the explicit changed source/test paths. The three held-scale hashes ran in both. Review then added their missing pinned-tier trait. On the final test source `Tier=Pinned&FullyQualifiedName~RandomEffectsTests.AHeldScaleFit_KeepsItsDraws` selected and passed all three. The four new ladder and SMC rows accounted for 259 seconds in the Commit run and are candidates for the next measured slow-tier refresh.

### Slice 0: the whitening on the working response (chip, 2026-11-08)

**Landed** in `7cf1761` (the change and its gates) and the commit that carries this report (`docs/Maths/RandomEffects.md`, begun, with one pointer line in `BayesianRegression.md`). The main session approved the plan on 2026-11-08; E34's chip had the change's shape before `BassDesign` was edited (no reply yet, its files clean).

**What changed.**
- `WorkingResponse.Whitening`: the block whose V^(−½) the reading's offset applies, its ln|V| in the log-Jacobian as before. `RandomInterceptBlock.Whitened` sets it; a reading moved to another offset keeps it.
- `BassDesign` drops its own block and its constructor's `block` parameter, which `BassModel.NewDesign` alone passed. It whitens each column by the reading it is on; whether its columns are whitened is its layout, fixed when it is built.
- `ProposeFollow` treats a reading at another scale as it treats moved weights: every column evaluated and whitened afresh, the Gram factored strictly. Gaining or losing a whitening is refused as another layout.

**A defect found and fixed with it: the hinge prefix tables under a block.**
- They give a one-way hinge's border by summing each held column against the raw hinge. Under a block the held columns are whitened, so they read xⱼᵀV^(−½)h where the border is xⱼᵀV⁻¹h.
- `WarmTables` already refused to build them under a block; `TableBorder` read them anyway. They switch on at 500 cases and up to 16 continuous predictors, past KHAOS's gates' n 246, so nothing caught it.
- On HEAD the new gate failed: with the tables asked for, a proposal the dots admitted came back refused as numerical. A design with whitened columns now reads every border by the dots.

**Gates.**
- *Held-s fits bit for bit:* a hinge chain, a Legendre chain and a tempered hinge fit with its reference rung, pinned by hash on HEAD before the change (`b9b95aee57ef459c`, `a1f4528886afa5f7`, `05122582e8f56f75`) and met after.
- *Following the whitening:* a design loaded at s = 1.3 and followed to a reading at s = 0.4, against a design built at 0.4 with the same terms: the residual and fitted squares and the coefficients equal to the bit (gated at 10⁻¹²). Followed back to 1.3, its first fit within 10⁻¹². A reading without a whitening is refused.
- *The tables under a block:* 600 cases, 1,500 proposals with the tables asked for, each equal to the dots' to the bit, the tables reading 0 entries.
- *KHAOS's gates* (`RandomInterceptTests`) and the prefix-table tests pass unchanged.
- *Fingerprints* do not move: nothing serialized changed.
- *Per commit,* the selection with the pinned tier: 1,845 of 1,845. *At the slice's end,* the whole assembly: 2,640 of 2,640.

**Decisions beyond the brief.**
- *The sweep's draw of the effects still reads the model's block.* Under a held s it is the reading's own object, and `BassSweep.cs` held the KHAOS chip's uncommitted birth-and-death edit, which I left alone. Slice 1 moves that read to the reading's block, which a drawn s needs.

**Noted for slice 1.** Under level centring the model's priors read the working response (`BassModel.PriorsReadResponse`), and `BassStructureMoves.ProposeFollow` refuses a chain of such a model to move. ȳ is unmoved by the whitening (V1 = 1), but s², which a scaled σ² prior reads, is the whitened response's. A drawn s therefore needs the priors to read the data, not the reading at s: a frozen reference of the unwhitened response's level and mean square, so the σ² prior does not move with s.

**Follow-ups.** None new.
