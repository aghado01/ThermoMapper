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

### Slice 1 — drawn scale on one Gaussian intercept grouping (2026-11-09)

The scoped commit `feat(bass): draw the Gaussian random-intercept scale on the collapsed target` carries the implementation, tests and [user documentation](../../../docs/Maths/RandomEffects.md). Verification used a snapshot of `541b772` plus this slice's files. The dedicated [receipts](../../../artifacts/e21-ridge-s1-qualification/) include the run settings, baseline and saved TRXs.

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
