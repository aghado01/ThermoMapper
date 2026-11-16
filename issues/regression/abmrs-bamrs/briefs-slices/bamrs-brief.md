# BAMRS: free-knot tensor-product spline surfaces (preparation brief)

Date: 2026-11-12. Status: preparation; nothing is built. This is E3 of the [extensions checklist](../../bass/extensions-checklist.md), queued by the user on 2026-11-08. It replaces the design notes under E3 in the [extensions brief](../../bass/extensions-brief.md#free-knot-tensor-product-surfaces-abmrs) as the working document. Written by the main session at the user's request, after a re-read of the sources and a survey of the engine. Its decisions are proposals until the user takes them. Revised 2026-11-15 after an external review the user relayed: [the record](#review-of-2026-11-15) lists what changed. The rollout was re-planned the same day at the user's direction: BAMRS is built as a new engine beside BARS and BAPS ([Decision 1](#decisions), [its record](#rollout-re-plan-of-2026-11-15)).

## Why

BAMRS fits a surface over two or three inputs with a tensor-product B-spline whose knots are free on every axis: reversible-jump births, deaths and relocations per axis, with the coefficients and the noise variance integrated out. In method it is BARS on a product carrier. In code it is its own engine, as BASS is. He, Yang and Kang's ABMRS is the published instance. BAMRS is this engine's version, on the corrected forms the regression family already holds.

What it adds to the family:

- **Smooth surfaces with exact per-draw geometry.** BASS's hinges default to degree 1, and BAPS's tensor splines keep fixed knots. A cubic tensor spline with free knots has exact derivatives, restrictions and extrema for each posterior draw, with the knot uncertainty carried.
- **Jump-like transitions along an axis,** from knots that cluster within the data's resolution. Each draw stays continuous, and its change across a gap in the data is free. A literal jump needs more ([item 4](#what-the-product-structure-changes)).
- **Scope.** Low dimension: d = 2 for cubics, and d = 3 at degree 1 or with many observations (the support gate, below). The engine is application-agnostic.

**The name.** BAMRS is Bayesian Adaptive Multivariate Regression Splines, the source's ABMRS in the family's word order. [BAMRS.md](../../../../docs/Maths/Regression/BAMRS.md) and [Regression.md](../../../../docs/Maths/Regression/Regression.md) expand it as "Multiresolution", which is a correction owed. A multiresolution engine is a different model and a separate candidate ([its section](#multiresolution-a-separate-candidate)).

## Sources

- **The preprint:** He, Yang and Kang, arXiv 2405.13353v1 (May 2024), where the method is called EBARS. Markdown at [papers/HYK2024.md](../papers/HYK2024.md); the source is held under `codex-scientiae/supellex-private/gauntlet/bayesian/2405.13353v1`. arXiv has no later version (checked 2026-11-12).
- **The article:** *Biometrics* 82(3), 2026, ujag121, doi 10.1093/biomtc/ujag121. Its title is *Adaptive Bayesian multivariate spline knot inference with prior specifications on model complexity*. **The main text is not held** and the journal page refuses automated fetch. Crossref holds the publisher's abstract and reference list for the DOI (`api.crossref.org/works/10.1093/biomtc/ujag121`, read 2026-11-15). Every statement here about the published version rests on that record or on the open supplement.
- **The supplement:** [ABMRS_Revision_V2_Supplementary_Materials](../ABMRS_Revision_V2_Supplementary_Materials), with the web appendix as [text](../ABMRS_Revision_V2_Supplementary_Materials/web_appendix_extracted.txt). Web Appendix A proves Lemma 1 and Lemma 3 only, with the preprint's statements, which implies the article carries no theorem needing a proof beyond the preprint's lemmas.
- **The author's repositories,** cloned at `D:\aghado01\packages\bayesian-related`: `ABMRS` (its `src` and `R` are identical to the supplement's copy) and `exampleABMRS`, the scripts behind the article's simulations and its real-data example.
- **Chen and Chen (2008),** *Extended Bayesian information criteria for model selection with large model spaces*, Biometrika. Not held. Its consistency condition is quoted below from memory and wants checking against the paper before it enters `docs/`.
- **Checked earlier:** the report half of [abmrs-verification-questions.md](../discussion/abmrs-verification-questions.md), which holds against the supplement.
- **Unreliable, and not to be cited:** [HYK-2024-2026-differences.md](../discussion/HYK-2024-2026-differences.md) and [HYK-2024-2026-clarifications.md](../discussion/HYK-2024-2026-clarifications.md). They are external summaries, and their detail is unverified. The earlier pass found three of their claims without a counterpart in the supplement's code: a knot-multiplicity kernel, a gradient and Hessian peak finder, and a changed death formula. [argon-abmrs-verification-debrief.md](../discussion/argon-abmrs-verification-debrief.md) retracts part of it.
  - *Withdrawn 2026-11-15.* This brief first listed their fMRI and cardiac applications as invented. The published abstract names an analysis of functional MRI data from the Human Connectome Project, so that application is real. The cardiac one is unverified: the reference list cites a cardiopulmonary-bypass study (De Somer and co-authors 2011), and nothing held says how the article uses it. What the summaries say of either application's detail still waits for the main text.

## The source's theory, re-derived

Notation is the module's (BARS.md): m observations, d axes indexed by a, degree pₐ, kₐ interior knots ξₐ on axis a, νₐ = kₐ + pₐ + 1 columns per axis, ν = Πₐ νₐ columns in all, Z the m × ν design, G = ZᵀZ, C the candidate count (the source's n).

- **Model.** y = Zβ + ε with ε ~ N(0, σ²I). Row i of Z is the row-wise Kronecker product of the axes' B-spline rows at xᵢ.
- **Structure prior** (their equation 3). Over a grid of Cₐ candidate sites per axis, π(k) ∝ τ(k)^(1−γ) and π(ξ | k) = 1/τ(k), with τ(k) = Πₐ binom(Cₐ, kₐ). The joint is π(k, ξ) ∝ τ(k)^(−γ). This is Chen and Chen's prior: uniform within a model size, and a size's total mass ∝ τ^(1−γ).
- **Coefficient prior** (equation 4). β | σ ~ N(0, mσ²G⁻¹), the unit-information g-prior, with π(σ) = 1/σ.
- **Lemma 1** holds, with a model-free constant. Integrating β gives y | σ ~ N(0, σ²(I + m·P_Z)), whose determinant is (m+1)^ν and whose inverse is I − (m/(m+1))·P_Z. The σ integral then gives p(y | k, ξ) ∝ (m+1)^(−ν/2)·a^(−m/2), with a = yᵀ(I − (m/(m+1))P_Z)y. It needs Z of full column rank.
- **The conditional draws are closed-form:** σ² ~ IG(m/2, a/2), then β | σ² ~ N((m/(m+1))β̂, (m/(m+1))σ²G⁻¹). The preprint's "Gibbs sampler" for them is unnecessary.
- **Lemma 2** is the BIC-type approximation: p̂ ∝ m^(−(ν+1)/2)·(σ̂²)^(−m/2)·τ^(−γ), with σ̂² the maximum-likelihood variance.
- **Lemma 3** holds for every γ ≥ 0. The schedule is b_k = c·min(1, ((C−k)/(k+1))^(1−γ)) and d_k = c·min(1, (k/(C−k+1))^(1−γ)). One of the two minima is always 1, so b_k/d_(k+1) = π(k+1)/π(k) exactly, the structure prior and the proposal cancel, and the acceptance is the marginal-likelihood ratio alone. It holds per axis because the count law factorizes and the axis is picked uniformly.
- **Where [0.5, 2] likely comes from.** This is my derivation; their reasoning is unread. Chen and Chen's condition is γ > 1 − 1/(2κ) when the candidates number O(m^κ). Their runs set C = 2m or 3m, so κ = 1 and γ > ½. γ > 1 lies outside Chen and Chen's range, and the supplement runs it (1.2, 1.5).
- **What the article adds to the preprint,** as far as the supplement shows: γ past 1, c = 0.3, the column filter, coefficient draws with prediction intervals, a peak-detection simulation, and a real-data example (triceps skinfold). The abstract adds a second real-data example, functional MRI data from the Human Connectome Project, and names three tasks: fitting discontinuous multivariate functions, detecting univariate change points, and locating peaks in multivariate regression. `exampleABMRS` holds no code for the fMRI analysis. The preprint's manifold-denoising application has no trace in `exampleABMRS`.

## Corrections to the source

Rows 2–5 are in BAMRS.md today. The rest are new from this pass.

| # | Source | What the check found | The engine's design |
|---|---|---|---|
| 1 | Lemma 1's prior is centred at zero (both versions) | a = RSS + ‖Zβ̂‖²/(m+1), so the marginal depends on the response's level, and a response far from zero swamps the residual. It is the defect BARS fixed on 2026-10-22 ([BARS.md](../../../../docs/Maths/Regression/BARS.md), "Centring the Normal prior"). Their code does not show it only because it runs the BIC-type ratio | The centred block (`GPriorCentre.Median` by default). A tensor B-spline basis sums to one, so the centre is c·1 with no solve |
| 2 | The code scores m^((ν−ν′)/2)·(σ̂/σ̂′)^m (`mebars.cpp:185`) | An approximation where Lemma 1 is exact and cheap. Its ν counts the columns it has dropped | The exact marginal of the centred block |
| 3 | Columns with at most `npart` = min(20, ⌊0.1m⌋) nonzero entries are dropped, and the Gram takes a 10⁻⁸ ridge (`utils.cpp:234`–`257`) | The model changes silently, and the basis stops summing to one | Refusal: the support gate and factor admission |
| 4 | β is drawn at the plug-in σ̂, and σ is never drawn (`mebars.cpp:199`–`202`) | The draws and the prediction intervals omit σ's uncertainty | The exact joint draw |
| 5 | At k = 0 a birth is forced, with b₀ in its acceptance; the reverse death omits 1/b₀ (`mebars.cpp:132`, `188`) | k = 0 is undersampled wherever a knot is favoured, by up to 1/b₀, and b₀ = c·C^(1−γ) is small for γ > 1. With d axes an axis without knots is an ordinary posterior state, so this matters more than in one dimension | The chain owns the selection ratio; an unavailable move stays put |
| 6 | The paper's schedule at k = 0 leaves relocation the weight 1 − b₀ − d₀ with no knot to move, and d₀ is 0⁰ at γ = 1 | Undefined as printed | As row 5 |
| 7 | The preprint's Algorithm 1, on a rejection, returns to the proposal step and repeats | That records accepted states only: the jump chain, whose law is the posterior reweighted by each state's acceptance rate. The code records the current state again, which is correct. Whether the article fixed the listing is unknown | Standard Metropolis–Hastings |
| 8 | γ ∈ [0, 1] in the preprint | The law is proper for every γ ≥ 0. For γ > 1 it is symmetric in k ↔ C − k | `EbicPrior(γ, C)` with C at least twice the largest count the data support (BARS.md) |
| 9 | γ sets the prior and the schedule together | Their γ comparisons do not separate the prior from mixing | Prior and schedule are separate here |
| 10 | Peak detection (`exampleABMRS/scripts/surface_peak.r`): degree 1, counts fixed at (5, 5), an OLS refit per knot draw, interior knot-lattice points tested against a ±0.02 neighbourhood, the per-draw global extremum averaged over draws | The idea is sound at degree 1, where extrema sit on the knot lattice (below). The script omits coefficient uncertainty and the boundary lattice, and there is no method for cubics | [Readouts](#readouts-how-the-bars-contract-extends) |
| 11 | Knot intervals come from fixed-count runs, as quantiles of the ordered knots; the triceps tests refit a `glm` at the posterior-mean knots | Ordered-knot summaries are valid only at a fixed count, and the plug-in ignores knot uncertainty | A fixed-count mode, and readouts taken per draw |
| 12 | Non-normal responses: a BIC-type plug-in, described and not coded | — | The engine's Laplace marginals, in a later slice |

## The model as BAMRS takes it

Every form below is settled for BARS and carries over. The product structure is the only new part.

- **Carrier.** One `KnotConfig` per axis: continuous knots on the open ordered simplex.
- **Structure prior.** A product across axes, restricted to the admissible set. Its log density there is Σₐ [ln p(kₐ) + ln kₐ!], the count law times the order-statistics density. The restriction couples the axes, so their counts are not independent under the prior (item 6 of the next section). Any `IComplexityPrior` serves an axis. `EbicPrior(γ, Cₐ)` reproduces the source's law, since a uniform draw of kₐ sites from a grid tends to the order statistics of kₐ uniforms as the grid refines; C stays as a nominal count.
- **Coefficients and noise.** The centred collapsed g-prior block, with the σ² priors BARS has (`Sigma2Prior`). Under the Jeffreys prior its marginal is Lemma 1 on ỹ = y − c.
- **Admission.** In order of cost: the simplex and the count priors; the size bound m > ν, read from the counts alone, with ν = Πₐ νₐ computed under overflow protection; the support gate on product columns; then factor admission (`DesignAdmission`). No design is allocated before the size bound passes (item 6 of the next section has the arithmetic). A refused configuration scores −∞. Truncating the prior to the admissible set keeps every acceptance ratio valid, since an inadmissible proposal is simply rejected.
  - The order is BAMRS's own. BARS's target keeps its order, in which the models refuse `TooFewObservations` after the design is built. So a configuration refused on two counts can carry a different reason in the two engines' tallies, while the admitted sets agree.
  - The gate needs no dense design. A product column's mass Σᵢ Πₐ Bₐ(xᵢₐ) accumulates from each row's Πₐ (pₐ + 1) nonzero entries. B-spline values lie in [0, 1], so that mass is at most each factor column's own mass, and an axis column below the minimum refuses before any product is formed.
- **Moves.** Birth, death and relocation on one axis at a time, with the chain owning the selection ratio. Constant weights are BARS's default and the prior-cancelling schedule its option (`DmgkSchedule`); the source's schedule is that option at c = 0.3 on the EBIC law.

## What the product structure changes

These follow from the model. None is in the source.

1. **A knot on one axis costs a whole slice of columns.** A birth on axis a adds Δν = Π over a′ ≠ a of ν_a′ columns, and the marginal charges (m+1)^(−Δν/2) for them.
   - At m = 500 with 8 columns on the other axis, that is 25 nats: the knot must cut the residual sum by about 9.5% to break even. In one dimension it is 3.1 nats, about 1.2%.
   - At d = 3 with 6 columns on each other axis it is 112 nats.
   - *Consequence.* Per-axis count posteriors will be sharp and count changes rarely accepted, so relocation carries the mixing. BARS found that tempering rarely pays (BARS.md, "When to temper"); that finding was for one dimension and may not hold here.
   - *The measurement and its decision.* Chains can agree on each axis's count and still hold different count pairs or knot sites. So a screen comes first, and two readings decide. The screen: start chains at different count vectors and read R̂ of each kₐ. Disagreement after the burn-in BARS uses shows the plain chain failing. Agreement settles nothing.
     - *The small case, for correctness.* With `EbicPrior(γ, 1)` on each axis the support is the four count vectors of [Validation](#validation)'s quadrature. The plain chain, a ladder and the flat-histogram sampler over counts are each read against it: the joint count occupancy, and the posterior mean of a surface functional (the fit at a point, or an integral over part of the box).
     - *The benchmark, for efficiency.* From the dispersed starts, compare the three samplers on effective samples per unit time, for the count-vector indicators and for a surface functional, and on whether the starts reach the same joint occupancy.
     - *Decision.* A tempered sampler becomes BAMRS's default where it reproduces the small case and, on the benchmark, either the plain chain's starts disagree or the tempered sampler gives more effective samples per unit time. Otherwise tempering stays optional as in BARS. The default is stated for the surfaces and sizes measured.
2. **The support gate bites at the corners.** A clamped boundary B-spline of degree p on a first interval of width h has integral h/(p+1). For inputs roughly uniform on the box, a corner column's mass is about m·Πₐ hₐ/(pₐ+1).
   - Under the default of one point's worth, cubics at d = 2 need h₁h₂ ≥ 16/m at each corner: at m = 500, first knots at 0.1 on both axes are refused (mass 0.31), and 0.2 by 0.4 is admitted (2.5).
   - At d = 3 the divisor is 64 for cubics and 8 for degree 1. Cubics with m = 1000 and intervals of 0.3 are refused (0.42).
   - The source's filter is harsher in count and silent in effect: it drops those columns and keeps the knots.
   - *The measurement and its decision* is under Decision 2.
3. **Size.** ν grows as a product and the block needs m > ν. In lexicographic order the Gram's half-bandwidth is about p·ν_inner + p, so the factor is cheapest with the axis of fewer columns innermost. That choice is a pure function of the state, so it does not disturb bit-for-bit resumption.
   - *Axis identity.* A birth can change which axis has fewer columns, so the cheapest order flips during a run. BAPS flattens coefficients as j·ν_y + k, and `TensorPenalty` strides by ν_y (`TensorPenalty.cs:101`–`105`). Its caller fixes the order once, which fixed knots allow. A flip that reached the coefficients would swap the axes under every consumer, with no error raised.
   - *Recommended.* BAMRS's assembler is its own, and BAPS's `TensorProductDesign` and `TensorPenalty` stay as they are. Outside the assembler a draw's coefficients are a d-way array addressed by multi-index (j₁, …, j_d), with one meaning whatever the counts. The physical order belongs to the assembler and its factor, and nothing else reads it. Three things can then change the physical form without reaching a consumer:
     - the size rule (at d = 3 the half-bandwidth is p times the product of every axis but the outermost, so the rule is the largest axis outermost);
     - a factor that is not a flat band;
     - the full-grid shortcut, which produces the array from per-axis factors with no band at all.
   - *Sequence.* The first build keeps the declared order, and reordering arrives with BAMRS's own Gram path. A penalty, if Decision 2 adds one, is per-axis data that the assembler places through `DifferencePenalty.AccumulateStrided`, which already takes a stride.
   - *Its test.* A rectangular design (ν₁ ≠ ν₂) with different penalties on the two axes, against dense Kronecker algebra, in both size orders so that the layout flips between them.
4. **A steep transition's columns span the box.** On this carrier the knots are distinct, so every draw is continuous with p − 1 continuous derivatives, and f(c⁺) − f(c⁻) = 0 at every c however tightly the knots cluster. What p + 1 knots inside one gap in the data give is a free change across the gap: the polynomial pieces on its two sides no longer constrain each other, and no observation sits in the transition. To the data that is a jump. Every column keeps its support: no B-spline has all p + 2 of its knots inside a cluster of p + 1, so each reaches data on one side. I argued this from B-spline support and did not run it against the gate.
   - *What spans the box is the cost.* The cluster's knots lie on a whole line (or plane) across the box, and each costs a slice of columns (item 1). The contrast across the gap, f(g⁺, θ) − f(g⁻, θ), is itself a spline in θ on the other axis's knots. Its amplitude can vary along the line and vanish on part of it. So a transition confined to part of the domain is representable, and it is paid for across the whole box.
   - *A literal jump* needs a discontinuous basis on that axis or a knot of multiplicity p + 1. `StepBasis` is a discontinuous basis the product already admits, piecewise constant along its axis. Multiplicity lies outside the open simplex and would be a new carrier. The source stands in the same place: its sampler draws distinct sites from an even grid (`mebars.cpp:22`–`33`, `:142`), while two of its test surfaces carry fourfold knots (`exampleABMRS/scripts/surface_fitting.r:61`, `:68`).
5. **Expect the mean function's intervals to under-cover.** The supplement reports 66–96% coverage at a nominal 95% (Web Tables 2c and 3), below 90% in most cells, against 93–95% for prediction intervals. My hypothesis is that sharp count posteriors (item 1) leave approximation bias outside the band, with their plug-in σ̂ a smaller part. It is a hypothesis; BAMRS's own coverage is a reported number, not a gate.
6. **Admission couples the axes' counts.** The prior is a product across axes only before admission. Under a constant likelihood the mass of a count vector k is ∝ Πₐ pₐ(kₐ)·Pr(admitted | k), where the second factor is the share of knot placements (uniform order statistics on each axis) that the target admits. It depends on all the counts at once, through m > Πₐ νₐ and through the product columns' masses.
   - *Size alone shows it.* At m = 500 with bicubics and `EbicPrior(1, 1000)` on each axis, 1,359 of the 1,002,001 count pairs satisfy (k₁ + 4)(k₂ + 4) < 500, by enumeration. Each kₐ is flat before admission. After it, an axis's count carries mass in proportion to its admissible partners: 121 of the 1,359 pairs at kₐ = 0, and 17 at kₐ = 20. "γ = 1 is flat in k" does not survive the product, and the location gates bend the law further.
   - *Qualification reads the joint law.* The chain and the reference sampler are tested on the law of the count vector, against the quadrature oracle of [Validation](#validation). A test of each axis's marginal against pₐ would fail a correct sampler.
   - *The reference sampler rejects whole states.* It draws every axis, then keeps or discards the state as one. Redrawing only a refused axis would sample a different law. The count law follows the tempering path, as in BARS (`SplineTarget.Reference`): Poisson(1) on each axis on the powered path, whose b = 0 member is the base measure; each axis's configured count prior on the likelihood path. Both take uniform order statistics for the sites.
   - *Size first, in the sampler too.* On the likelihood path the flat law above fails the size bound 99.86% of the time, and a draw at (1000, 1000) would ask for a 500 × 1,008,016 dense design, 3.76 GiB. The sampler tabulates the count vectors that pass the size bound and draws k from the product law restricted to them. That is exact, since the bound is a function of k alone and every admitted state passes it. The location gates then act by whole-state rejection.
   - *Lemma 3 is untouched.* The schedule cancels the ratio of pₐ, and admission enters each acceptance as an indicator.

## What is already built

Read from the source on 2026-11-12. "As is" means by reading, not by running. Nothing in this table is edited in place: BAMRS calls it as it stands, or adds an entry beside it (Decision 1's rule).

| Piece | Where | Use for BAMRS |
|---|---|---|
| Per-axis carrier and basis | `KnotConfig`, `SplineBasis`, `BSpline`, `StepBasis`, `IBasis` | As is, one per axis. `IBasis` needs no change: the product is taken over the axes' own designs, so axes may mix bases (a step axis with a cubic axis) |
| Tensor design | `Baps.TensorProductDesign` | Left as it is: bivariate, over `SplineBasis`, in the order its caller picks ("put the smaller dimension as y"). BAMRS builds its own d-axis assembler over `IBasis` in its folder (item 3) |
| Face-splitting contraction (N-X1) | `KhatriRao.FaceSplittingGram`, `FaceSplittingCrossProducts`, `FaceSplittingHalfBandwidth`, `KroneckerGram`, `FaceSplittingRow` | As is. Built for this; only tests call it |
| Wide-band factor (N-X2) | `WideBandCholesky`, `SpdFactors.ForBand` | As is. Built for this; nothing selects it |
| Observation models | `IObservationModel`: `WeightedNormalModel`, `PoissonModel`, `WhittleModel`, `GibbsLossModel` | As is in correctness: each takes a dense design and knows nothing of the basis. `BandedDesign` reads each row's first and last nonzero column, so a tensor row's span is p·ν_inner + p wide and the accumulation costs O(m·bw²). At the source's benchmark sizes (m ≤ 1000, ν ≲ 150) that is small. When sizes ask for more, BAMRS forms its own Gram by `FaceSplittingRow` with a factor from `SpdFactors.ForBand`, behind its own seam, and the models stay as they are |
| Collapsed block | `CollapsedGPrior`, `BandedGPrior`, `GPriorCentre`, `Sigma2Prior` | As is: any `ISpdFactor`, the centre conventions, tempered marginals, leave-out densities |
| Target | `SplineTarget` | Left as it is. BAMRS's target is new, in its folder, with the same composition: per-axis count priors, the size bound ahead of any design, the gate on product columns, the model's marginal, and a reference sampler that draws every axis and rejects the whole state (item 6) |
| Moves | `KnotBirthMove`, `KnotDeathMove`, `KnotRelocateMove`, `KnotBirthDensity`, the kernels, `DmgkSchedule`, `ScheduledMove` | As is, wrapped: BAMRS's adaptor applies an `IRjMove<KnotConfig>` to one axis of the product carrier. The chain normalizes selection over the 3d moves, and a null proposal at kₐ = 0 stays put, which is the balance row 5 asks for |
| Chains and ensembles | `ReversibleJumpChain<TState>`, `ParallelTempering<TState>`, `TemperingLadder`, `ChainEnsemble`, `IEnsembleModel<TState>`, `LadderEvidence` | As is: generic in the state |
| Other samplers | `SmcRun`, `MulticanonicalRun` and the flat-histogram family | Candidates for item 1's count barriers. Their contracts are unread |
| Orchestration | `BarsEnsemble : IEnsembleModel<KnotConfig>` (82 KB): chains, robust weights, noise, checkpoints, evidence, the conformal hooks | Left as it is, and read as a pattern. BAMRS has its own ensemble on the generic shell, as BASS (`BassEnsemble`) and BAPS (`BapsEnsemble`) have theirs (Decision 1) |
| Curve readouts | `SplineExtrema`, `BarsDraw`, `PeakDraw`, `PeakSpan`, λ, π, `PeakModes` | `SplineExtrema` is called as it stands on a restriction. The draw and peak records are BARS's, and BAMRS has its own. [Readouts](#readouts-how-the-bars-contract-extends) |
| Moments over piecewise polynomials | `Bass.Sensitivity.PiecewisePolynomial`, `BassMoments` | Sobol indices and effects for tensor B-splines, as the extensions brief notes |
| Declared coefficient priors | `PriorPencil`, `CollapsedPencilPrior`, `TauPrior`, BAPS's `TensorPenalty` | Decision 2's ridge option |
| Term contract (F-X3) | `Bass.IDesignKernel<TPoints>`, `DesignFactor` | The other route: a tensor B-spline part as a BASS term type. Not recommended as the first build (Decision 1) |
| Lean | [BARS.lean](../../../../lean/Enthymemata/BARS.lean): MP-1, MP-2, SP-1 to SP-3 stated, proofs open; protolemmata MP-3, MP-4, O1–O3; [Stability.lean](../../../../lean/Enthymemata/Stability.lean) | [Lean](#lean) |

## Decisions

Each is the main session's recommendation, with its alternatives. The user may take or reverse any.

1. **A new engine beside BARS and BAPS.** Decided by the user on 2026-11-15. It reverses this brief's first recommendation, which lifted `BarsEnsemble` over the product carrier and made BARS the one-axis member.
   - *Why.* BARS has grown well past the original Kass method, and BAPS is fitted work. Unifying an engine while building a new one puts both at risk for the new one's sake. BASS was built the other way, as its own engine on shared primitives, and BAMRS follows it.
   - *The shape,* after BASS and BARS. Its own folder, `src/maths/regression/spline/bamrs/`, in `Maths.Regression`, with tests under `tests/regression/`. Its own state, target, move adaptor, ensemble, and draw and result records. A declarative config with a static entry point. The names follow the family's pattern: `BamrsConfig`, `BamrsState`, `BamrsTarget`, `BamrsEnsemble`, `Bamrs.Fit`.
   - *The reuse rule,* in order of preference:
     1. Call a primitive as it stands.
     2. Add to a shared primitive beside what is there: a new member, overload or type, with every existing signature and path left as it is. The commit gates on every test project from disk, with BARS's and BAPS's pinned tiers unchanged.
     3. Where neither serves, BAMRS takes its own variant in its folder, and the split is logged as a follow-up. Whether to unify is a separate decision for the user, after BAMRS stands.
   - *What this withdraws.* Lifting `BarsEnsemble`; the tensor builder shared with BAPS; BARS's pinned tier as an acceptance through a lifted engine; and opening `SplineExtrema` in place.
   - *What it gains.* BARS becomes an independent oracle for BAMRS at d = 1 ([Validation](#validation)). The first plan made the two one code path, which could not check itself.
   - *Its cost.* BAMRS's ensemble repeats the shape of BARS's. The first build carries only what a Normal fit needs: the knot block, retained draws, the fit and its spread, checkpoint and resume. Robust weights, the noise block, evidence and the conformal hooks each arrive in their own slice when asked for, reusing BARS's primitives where those take a design and not a carrier.
   - *Still not recommended: a BASS term type.* BASS sums many sparse terms under its term measure; ABMRS is one full tensor term with moving knots. The fit is poor, though a tensor part inside BASS remains a later option through F-X3.
   - *Open until slice 0:* which of BARS's later blocks (weights, noise, the conformal hooks) take a design and not a carrier. I read `BarsEnsemble`'s surface, not its body.
2. **Thin product columns: the g-prior with the gate, and a ridge option only on evidence.**
   - *Recommended.* Keep the g-prior and the gate, the source's model in its corrected form, as the default.
   - *The measurement and its decision.* The question is how much posterior mass the gate excludes. The share of evaluations refused as `ThinSupport` does not answer it, since that share follows where the proposals go. So on the source's benchmark surfaces at d = 2, run the same g-prior model with the gate off and factor admission kept, and record each retained draw's least product-column mass. The posterior share of draws below the default minimum is the mass the gate removes, and the record answers for every other minimum too. If that share is material, add a declared ridge or penalty prior through `PriorPencil`, under which thin columns are regularized and the gate can relax. If it is not, the option waits for a caller. The refusal share and the counts the source reports stay as context. What share counts as material is fixed in the slice brief before the first reading, with the run's size and seeds.
   - *Alternative: lower the gate's default with d.* It trades refusals for coefficients the g-prior leaves almost unconstrained.
3. **The count law's default is open.** `EbicPrior(1, 2m)` matches the source's runs and gives comparable tables. BARS's `PoissonPrior(5)` per axis implies about 81 columns a priori for a bicubic. Both figures are before admission, which makes the counts' prior joint (item 6). Take the default after the benchmark, not before.
4. **Fixed per-axis counts are a mode.** A relocation-only palette at given counts. The article's knot intervals and its peak simulation both use it, and it is cheap.
5. **Readouts arrive in tiers** (next section): affine and restriction readouts first, degree-1 lattice extrema second, joint cubic extrema last.
6. **The full-grid shortcut is its own slice.** On a complete grid Z = Z₁ ⊗ Z₂ and ỹᵀP ỹ = ‖Q₁ᵀ Ỹ Q₂‖²_F, so a move on one axis is a one-dimensional refit. Missing cells and robust weights that do not factor across axes both break it, and the scattered path then serves.
7. **Non-normal responses wait** until the Normal model qualifies. `PoissonModel` takes a tensor design as it takes any other.

## Readouts: how the BARS contract extends

**What BARS has.** Each retained draw holds its knots and one joint draw of (β, σ²). `SplineExtrema` reduces that draw exactly: the candidate set is the endpoints, the knots and each span's derivative roots, in closed form for degree ≤ 3. Peaks pass a prominence gate, and each peak's span is the level crossing at a drop fraction of its prominence, again as a closed-form root. The pooled views are λ(T), π(T) and `PeakModes`. This is a cubic method, and it is a method for one variable.

**What carries over unchanged.**
- The rule: the posterior mean of a readout affine in the surface pools each draw's conditional mean, and nonlinear readouts use the joint draw.
- The order: reduce each draw exactly, then pool. MP-2 and SP-3 are statements about averaging and do not depend on the dimension.
- The affine readouts: the fit, partial derivatives and integrals over an axis.
- *Spread needs the coefficients.* `FitSd` is not affine. BARS takes it from the law of total variance: the variance over draws of the conditional mean, plus the mean over draws of the conditional variance (`BarsEnsemble.cs:573`–`577`). The spread of every affine readout has the same two terms. A fixed-knot fixture isolates the second: every draw has the same conditional mean, and the fit's posterior spread is still positive.
- *Sobol outputs are functionals of a draw.* A variance component is quadratic in the coefficients, and an index is a ratio of two. Neither is a functional of the conditional mean. Each is computed on the draw's joint (β, σ²) and summarized over draws, as BASS does (`SobolDecomposition`). A variance component βᵀQβ also has a closed-form conditional expectation, β̄ᵀQβ̄ + tr(Q·Cov(β | ξ)); an index has none. The decomposition of the posterior mean surface is a different quantity, which leaves out the posterior spread, and it takes its own name if it is ever reported.

**Tier 1: restrictions reuse the curve engine as it stands.** Fix every axis but one. f(·, θ) is a spline in T on axis 1's knots, with coefficients cⱼ(θ) = Σₗ βⱼₗ Bₗ(θ). `SplineExtrema` applies to it unchanged, so each draw gives the ridge T*(θ), its height, its significant peaks and their spans along T, exactly, at any θ. Pooled over a θ grid these are λ(T | θ) and π(T | θ), and a posterior band for the ridge.
- *What needs no change, and what needs an entry beside.* `SplineExtrema` takes (`KnotConfig`, coefficients, `SplineBasis`).
  - An axis-aligned restriction is a spline on axis 1's own `KnotConfig` with coefficients cⱼ(θ), so `SplineExtrema` reads it as it stands.
  - A derivative spline (below) is a spline of degree p − 1 on the same knots, which the same entry reads with a `SplineBasis` of that degree.
  - Only an oblique line asks for more. There a bicubic restricts to a piecewise polynomial of degree 6, which is no spline on either axis's knots. Its reader takes breakpoints, an evaluator and a degree, and the general-degree path already handles degree 6. By Decision 1's rule it is an entry beside `SplineExtrema`'s own or a type in BAMRS's folder, and BARS's calls are untouched.
- *The ridge is a peak set at each θ.* T*(θ) jumps where two ridges exchange dominance. So the pooled readout is the peak set along T at each θ, λ(T | θ), and the argmax is one view of it.
- *Ridge births and ends belong to this tier.* The interior ones are found among the roots of one polynomial in θ, and the rest are a branch leaving the box (see the fold on a surface, below).

**Tier 2: degree-1 surfaces have a closed form.** On each cell a degree-1 tensor spline is affine in each coordinate. A function affine in each coordinate attains its maximum over a box at a vertex, so every strict local extremum of the surface sits on the knot lattice, the boundary lattice included. A weak extremum can fill a line, as a ridge of constant height along a cell edge does. A bilinear cell's one interior critical point is a saddle. So at degree 1 peak inference is knot inference, which is the source's simulation, here with coefficient draws and the boundary.

**Tier 3: joint extrema of a cubic surface are new mathematics.** Fermat's candidate set stratifies over the cells: lattice vertices; critical points of the cubic restricted to each cell edge, which the existing closed form finds; and interior solutions of ∂₁f = ∂₂f = 0. The interior system has bidegrees (2, 3) and (3, 2), so at most 2·2 + 3·3 = 13 isolated solutions per cell. A cell can also hold a curve of critical points, as −(u − v)² does. That is the surface's counterpart of a flat span, and it needs its own result.
- *Recommended solver.* Convert each cell's patch to Bernstein form and subdivide. The control net bounds the patch, so cells are excluded with a certificate, and the global argmax needs only branch and bound on that bound.
- *Alternatives.* A resultant reduces the system to one polynomial of degree 13, which is ill-conditioned in the monomial basis. Newton from seeds gives no completeness.
- *What changes in the contract.* "Exact, in closed form" becomes "enclosed to a stated tolerance". The guarantee is one-sided: no critical point lies outside the returned boxes, and that a box holds exactly one is a separate test. The `argmax_in_closed_form_set` specification weakens accordingly for cubic surfaces, and holds as stated at degree 1.
- The same solver finds saddles, so it is one tool for two readouts. The fold points of ridges need only one-variable roots (below).

**Significance in two variables.** Prominence generalizes as the persistence of a maximum in the superlevel-set merge tree, and a span as the superlevel component at the drop level: a region, with an area and axis extents. Both need saddles. Obligation O3 already assigns the merge and col definitions to the planned carrier-generic π₀ and merge-tree foundation; the surface readout should consume that foundation too, and wait for it.

## Features on spline curves beyond peaks

This answers the user's wish to characterize more kinds of feature. It is a proposal. The readouts are functions of one curve (knots, coefficients, basis), so they serve every Tier 1 restriction of a surface and would serve a BARS draw as well. They are built with BAMRS. Whether BARS takes them up is a separate decision, outside this rollout.

**The organizing fact.** The derivative of a degree-p spline is a degree-(p−1) spline on the same knots, with coefficients from differences. So one extrema-and-crossings engine, run on f, f′ and f″, reads features of every order:

| Order | Zero set | Features | For a cubic |
|---|---|---|---|
| 0 | f = ℓ | Level crossings, spans, occupancy above a level | Cubic roots; built |
| 1 | f′ = 0 | Peaks and valleys; monotone stretches between them | Quadratic roots; built |
| 2 | f″ = 0 | Inflections: steepest ascents and descents, a peak's flanks, shoulders | Linear per span; new and trivial |
| 3 | jumps of f‴ | The knots | The knot posterior |

- The last row is why "knot inference" and "peak inference" coincide in the source only at degree 1: there the order-p features, the knots, are the order-1 features.
- **Steepest slopes** are the peaks and valleys of f′: the curve engine on the derivative spline.
- **Plateaus** are level crossings of f′ at ±ε: the span arithmetic on the derivative spline.
- **Inflection width.** The distance between the two inflections flanking a peak needs no baseline. The half-prominence span does, and O3 records that its level is a convention.
- **Apex curvature** f″(x*) and **left and right half-widths** (already in `PeakSpan`) give sharpness and asymmetry. The area above the drop level is a closed-form integral.
- **Monotonicity on [a, b]** has a posterior probability for each direction: the share of draws with f′ ≥ 0 throughout [a, b], and likewise with f′ ≤ 0. The test is the sign of f′'s least (or greatest) value on the interval, from the curve engine on the derivative spline. Counting critical points is the wrong test: (T − ½)³ is strictly increasing through an interior stationary point.
- **Transitions and kinks.** With distinct knots f(c⁺) − f(c⁻) = 0 on every draw, so a cluster of p + 1 knots inside one data gap has no jump size to read. Its readout is the contrast across a named interval, f(g⁺) − f(g⁻), with [g⁻, g⁺] a gap in the data or a window the caller gives. The contrast is affine in the coefficients, so its mean and spread follow the fit's. On a surface it is a spline along the other axis ([item 4](#what-the-product-structure-changes)). The steepest-slope readout above finds where to place the interval. At degree 1 each knot is a slope change, which is the quantity the source's triceps tests address by plug-in.
- **The topological summary** is the curve's 0-dimensional persistence diagram, exact from the candidate set. Prominence is persistence under the elder rule, up to the boundary convention of O1, and `Stability.lean` already pushes a sup-norm band forward to it.

**Shoulders and the fold.** A shoulder is a local minimum of |f′| at a nonzero value: a root of f″ where f′ keeps its sign. On a cubic span with f′(u) = Au² + Bu + C and discriminant Δ = B² − 4AC, taken where the inflection u* = −B/(2A) lies inside the span:
- f′(u*) = −Δ/(4A), the shoulder's residual slope;
- for Δ > 0 a maximum and a minimum straddle the inflection, separated by √Δ/|A|, with height difference Δ^(3/2)/(6A²);
- for Δ < 0 there is no critical pair, and |Δ/(4A)| measures how far the span is from having one.

So one closed-form scalar per span orders a strong peak-and-valley pair, a marginal pair, a shoulder and a plain monotone stretch. Its zero set is the fold in MP-3's catastrophe set ("a critical pair annihilates: the per-span quadratic discriminant crosses 0"). A shoulder and a marginal peak are the two sides of that fold. The newborn pair's prominence is the height difference above, which starts at zero, so a positive prominence gate is passed beyond the fold, where Δ^(3/2)/(6A²) reaches the gate. One feature then has two posterior probabilities: that the pair exists (Δ > 0), and that its peak is significant (the gate is passed). The fold moves the raw count and leaves the gated count alone. In general, where a pair straddles a knot, the definition goes through the derivative spline and not the per-span discriminant.

**On a surface the fold becomes a point.** In a family f(T, θ), a ridge is born or ends where ∂_T f = 0 and ∂²_T f = 0 together: the parameter value at which a shoulder becomes a peak. On a bicubic cell ∂_T f is a quadratic in the span coordinate whose coefficients A(θ), B(θ) and C(θ) are cubics in θ. So the ridge branches are the quadratic formula in θ, and the roots of Δ(θ) = B(θ)² − 4A(θ)C(θ), one polynomial of degree at most 6, are the candidate folds.
- *A root is a fold where Δ changes sign,* with A ≠ 0 and the double root u* = −B/(2A) inside the span. A pair of branches then exists on the side where Δ > 0 and on no other. At the double root Δ′ = −4A·∂_θ∂_T f, so a simple root of Δ is the transversal fold, the generic case.
- *A root of even order is no birth.* ∂_T f = (T − ½)² − (θ − ½)² has Δ = 4(θ − ½)², which vanishes at θ = ½, and it has two stationary points on both sides. The branches cross there.
- *A cell's branch count also changes with no fold,* when a root passes through the span's end into the next span or out of the box. A ridge's ends are its folds and its exits from the box.

The root finder behind `SplineExtrema`'s general-degree path covers degree 6, so ridge births and ends are a Tier 1 readout.

## Lean

Everything here is at the protolemma stage: candidate statements, none compiled. Their home, once worked, is a note under `lean/Protolemmata/` beside the BARS ones.

**Where the estate stands.**
- `lean/Enthymemata/BARS.lean` holds the first statements of MP-1, MP-2 and SP-1 to SP-3, with six open proofs.
- `lean/Protolemmata/BARS/BARS.v2.lean.md` repairs them:
  - MP-1 loses its interval, which blocked `finite_local_max` at the endpoints.
  - SP-2 gains its antitone dual, and a peak's width becomes the difference of the two bracket roots. `levelWidth` is read as occupancy above a level.
  - SP-3 splits into `levelWidth_not_affine`, which two constants witness, and a version whose level is tied to each curve.
  - MP-2′ joins MP-2 as the overcount witness, so the count has no signed bound under averaging.
  - Its header names what keeps it from compiling: a renamed root-finiteness lemma, now `Polynomial.finite_setOf_isRoot`, and two measure rewrites.
- `lean/Protolemmata/BARS/spc-bars-readout-obligations.md` owns the gap to `SplineExtrema` as O1 to O3. Its discharge sequence is: the revised scaffolding and O1; a piecewise-polynomial carrier and O2; the branch and merge-level carrier; the readout specification; the witnesses restricted to the carrier; the posterior measure carrier.
- `lean/Protolemmata/BARS/Fable5-BARS-revisions.md` adds a general pair: the raw count is unstable under a small perturbation, and the thresholded count is stable. It places the stable half downstream of `Stability.lean`, and moves MP-3 into merge-tree definitions on a connected-components foundation. The active stages hold no such module yet, and `Stability.lean` still takes curves on ℝ with no tameness hypothesis.
- O1 restates `argmax_in_closed_form_set`, and records that the old `argmax_expectation_noncommute` sketch does not state its inequality. A witness that does: f = −a(x − p)² and g = −b(x − q)² have a mean whose argmax is (ap + bq)/(a + b), which differs from (p + q)/2 whenever a ≠ b and p ≠ q.
- On the sampler side, `lean/Enthymemata/BassStructureChain.lean` has SC-1 with no open proof: `mhRatio_mul_swap`, `not_isMhPair_of_mul_ne_one` and `mh_detailed_balance`.
- `lean/Enthymemata/BassIntegrals.lean` has a d-variable `tensorExpansion` with `tensorPartial` and `hasDerivAt_tensorExpansion`, fixed to Legendre factors.

**Order.** Landing the revised draft comes first, and new statements are written against it. They join the ledger's sequence and add no step ahead of it. BAMRS changes what the second step's carrier has to be.

**The carrier.** Polynomial pieces on a finite break set, with continuity as a hypothesis.
- It is closed under each operation the readouts use: differentiation, the restriction of a tensor product to a line, and the difference of two draws on their merged breaks.
- It matches Tier 1's piecewise-polynomial reader, which needs breakpoints, an evaluator and a degree.
- The B-spline coefficient form is a later bridge onto it. Partition of unity and the numerical contract belong to that bridge.
- A surface is one polynomial per cell of a box partition. `tensorExpansion` shows the pattern over `Fin d → ℝ`. Opening its factor family past Legendre would give BASS and BAMRS one tensor carrier; that is a candidate, not a decision.

**Four statements that shape a contract.** These are worth stating before the code exists, and proving after.
1. *Maxima are strict.* A function affine in each coordinate attains its maximum over a box at a vertex, by induction on the coordinates. So the strict local maxima of a degree-1 surface are lattice vertices. Weak maxima can fill a line, which is why the predicate is strict, as O3 already asks for curves.
2. *The ridge readout is a peak set at each θ.* The argmax ridge jumps where two ridges exchange dominance, so the pooled object is MP-4's peak set, read at each θ.
3. *Tier 3's enclosure is one-sided.* No critical point lies outside the returned boxes. That a box holds exactly one is a separate claim.
4. *Gate and fold are two readouts.* The fold statements below give each its coordinate.

**Curve features.**
- *The derivative ladder is by instantiation.* MP-1 and SP-1 hold for any polynomial, so steepest slopes and inflections are those statements at `derivative p`. What is owed is that the carrier is closed under differentiation. O2 already takes every knot as a candidate, which covers the lower smoothness of f′ and f″. An inflection is a sign change of f″, and a span with f″ ≡ 0 is the degenerate case, as the flat span is for peaks.
- *The fold on a cubic span.* For f′(u) = Au² + Bu + C with A ≠ 0: a critical pair exists exactly when Δ > 0; its two values differ by Δ^(3/2)/(6A²); and for Δ < 0 the least |f′| is |Δ|/(4|A|). Mathlib's `discrim` and `quadratic_eq_zero_iff` carry the first part.
- *The fold leaves the gated count alone.* The newborn pair's prominence is that difference of values, which starts at zero. So the pair fails any positive gate until Δ^(3/2)/(6A²) reaches it. MP-3's event (a) is then a jump of the raw count, and the gated count jumps only on the threshold set. In general this is the thresholded-count stability that the revisions note places downstream of `Stability.lean`; here it has a closed form. Read against `SignificantPeaks` and `ValleyFloor`, it holds for a pair born on a monotone stretch. Peaks of equal height are a tie convention there.
- *Monotone is a sign condition.* f′ ≥ 0 on [a, b] gives a nondecreasing f, by Mathlib's `monotoneOn_of_deriv_nonneg`. (T − ½)³ witnesses that having no interior critical point is sufficient and not necessary.

**Surfaces.**
- *Restriction.* Σⱼ Σₗ βⱼₗ Bⱼ(x) Cₗ(θ) = Σⱼ (Σₗ βⱼₗ Cₗ(θ)) Bⱼ(x) for any two families of functions. It is an exchange of sums and needs no spline theory. It reduces Tier 1's exactness to O1 and O2.
- *The ridge in closed form.* On a bicubic cell ∂_T f = A(θ)u² + B(θ)u + C(θ) with cubic coefficients. The ridge branches are the quadratic formula in θ, and the candidate folds are the roots of Δ(θ), of degree at most 6. Where Δ > 0 the ridge's slope is −∂_θ∂_T f / ∂²_T f. The regularity that slope needs is Δ ≠ 0.
- *Which roots are folds.* With A ≠ 0 and Δ = 0, Δ′ = −4A·∂_θ∂_T f at the double root: an identity in A, B, C and their derivatives. So a simple root of Δ is a transversal fold, and the pair exists on exactly the side where Δ > 0. The witness that a root need not be a fold is (T − ½)² − (θ − ½)²: its Δ = 4(θ − ½)² has a double root, and two stationary points stand on both sides of it.
- *The bilinear saddle.* The interior critical point of a + bx + cy + dxy with d ≠ 0 is a saddle.
- *The convex-hull bound.* Nonnegative weights that sum to one keep a combination between its least and its greatest coefficient. It has two uses:
  - With Mathlib's `bernstein_nonneg` and `bernstein.probability` it is the exclusion certificate for statement 3.
  - With the B-splines' partition of unity it bounds a curve's sup-norm by its largest coefficient. That bridges coefficient perturbations to the sup-norm `Stability.lean` consumes.
- *The stratified candidate set.* Lattice vertices, edge critical points and interior critical points exhaust the local maxima of a continuous piecewise-polynomial surface: O2's statement in two variables. Finiteness needs a regularity hypothesis, since a cell can hold a curve of critical points, and that degeneracy gets its own result.
- *Stability for surfaces.* Give `Dgm` a domain type and thread the tameness hypothesis; both are cheap while the types are opaque. A surface band's ε needs only an upper bound on the sup of a difference of two draws. The convex-hull bound gives one on merged cells, and it errs on the safe side.

**Sampler and model.** These sit beside SC-1, since they serve every reversible-jump engine.
- *The source's k = 0 pair.* By `mebars.cpp:132` and `:188` the printed ratio from k = 0 is b₀·L, and the one back from k = 1 is 1/L, with L the likelihood ratio. Their product is b₀ ≠ 1, so `not_isMhPair_of_mul_ne_one` applies as it stands: an instance, with no new proof.
- *The schedule cancels the count prior* (Lemma 3): r·min(1, 1/r) = min(1, r) for r > 0, at r = ((C − k)/(k + 1))^(1−γ).
- *A mixture of balanced move pairs is balanced,* with weights that depend on the state and a null proposal as a stay step. This is the bookkeeping `ReversibleJumpChain` owns.
- *Repeating until acceptance changes the target* (Algorithm 1). If a kernel is reversible for π with holding probability h < 1, its accept-only chain is reversible for π(x)(1 − h(x)).
- *The corner mass.* ∫₀ʰ (1 − x/h)ᵖ dx = h/(p + 1), and the product over axes by Mathlib's `integral_fintype_prod_eq_prod`, as `integral_cube_prod` uses it for polynomial factors. It checks item 2's divisor, which is arithmetic so far.
- *The readout rule.* Its affine half is linearity of expectation. Its other half is the witness family: MP-2, MP-2′, SP-3 and the argmax witness above.

**Left as citation or arithmetic:** the dimension cost of an axis knot, Chen and Chen's condition, the bound of 13 (the solver does not use it), and Lemma 1's Gaussian integral. The Kronecker shortcut's identity waits for its slice, since the harness works back from built code.

## Multiresolution: a separate candidate

The docs' mistaken expansion points at a real successor. A multiresolution spline refines locally: hierarchical B-splines and their truncated form, locally refined splines, T-splines, and multilevel B-spline approximation. I name these from memory (Forsey and Bartels 1988; Giannelli, Jüttler and Speleers 2012; Dokken and co-authors 2013; Sederberg and co-authors 2003; Lee, Wolberg and Shin 1997). None is held, and a Bayesian free-refinement version wants a prior-work survey before any design.

- **What it would fix.** Three of the items above: a refinement adds a few columns where the data ask, not a slice (item 1); refinement stays where the data are, which eases the corner rule (item 2); and a steep transition is paid for where it is, not across the box (item 4). The preprint's own closing remark points the same way: the domain must be a box, and "a flexible multivariate spline" in place of the tensor product is its suggested remedy.
- **What it would reuse.** Each basis function is still a product of one-variable B-splines, so the face-splitting rows, the collapsed block, the reversible-jump chain and the readouts all serve. Its moves refine and coarsen cells, closer to BASS's term births and deaths than to knot moves.
- **What it would cost.** The Gram loses its band, the basis needs a linear-independence rule (truncation supplies one), and the structure prior becomes a prior over refinement trees.
- **A one-variable member exists too:** nested dyadic knots on a curve.
- **Naming.** Both models would abbreviate with M and R. BAMRS stays with the multivariate engine, to match its source; the multiresolution one takes a name when it is real.

## Slices

Proposed order; engine first, with no consumer wired in. Every slice lands in `spline/bamrs/` and its tests. A slice that touches a shared file does so by an entry beside what is there, and its commit gates on every test project from disk with BARS's and BAPS's pinned tiers unchanged (Decision 1).

0. **Read and measure.** Read `BarsEnsemble` and `BassRun` as patterns for the ensemble's shape. List the shared types BAMRS will call, each marked as it stands or as needing an entry beside. Run a bicubic design through `WeightedNormalModel` as it stands, and check the marginal against dense algebra.
1. **State, design, target, moves.** The product state (one `KnotConfig` per axis), the d-axis assembler with its coefficient form (item 3), the target with size-first admission and the whole-state reference sampler (item 6), and the axis adaptor over BARS's knot moves. At d = 1 the target's evaluations are read against BARS's `SplineTarget`. Qualify on the generic chain under a constant likelihood: the joint count law after admission (item 6), kₐ = 0 included, for the chain and for the reference sampler on both tempering paths, against the quadrature oracle.
2. **The gate's measurement** (Decision 2) and **the count-mixing measurement** (item 1), with their decisions. Each reading's size, seeds and rerun rule are fixed in the slice brief before it is taken.
3. **The ensemble.** `BamrsEnsemble` on the generic shell, with the config and the entry point: the knot block, retained draws with their joint (β, σ²), the fit and its spread, checkpoint and resume. Normal responses with fixed weights. Acceptance: resume bit for bit, and at d = 1 agreement with a BARS run on the same data inside a Monte Carlo gate stated beforehand.
4. **Tier 1 readouts,** with the order-2 curve features and ridge births and ends. `SplineExtrema` is called as it stands, and the piecewise-polynomial reader for oblique lines is the one new entry.
5. **Fixed counts and Tier 2.**
6. **BAMRS's own Gram path:** accumulation by face-splitting rows, the band factor and the size-ordered layout, when a measured size asks for them.
7. **The full-grid shortcut.**
8. **What BARS grew, one slice each when asked for:** robust weights, the noise block, evidence, the conformal hooks.
9. **Tier 3 and the merge tree,** after the π₀ foundation.
10. **Non-normal responses.**

## Validation

Ground truth comes from outside the estimator, and the source's code is a sanity band, never a parity target.

- **Dense algebra** for the marginal, the centre and the conditional draws on small tensor designs.
- **Quadrature** when each axis holds at most one knot. Each of the count vectors (0, 0), (1, 0), (0, 1) and (1, 1) needs an integral over at most two knot positions.
  - *With the data* the integrand is the marginal likelihood. It gives the posterior over the four and the posterior mean of a surface functional. With `EbicPrior(γ, 1)` on each axis the four are the whole support.
  - *Under a constant likelihood* the integrand is the admission indicator, and the integral is Pr(admitted | k). The law of k, given that it is one of the four, is Πₐ pₐ(kₐ)·Pr(admitted | k) normalized over the four. That holds under any count law, so it checks the chain and the reference sampler on the powered path (Poisson(1) counts) and on the likelihood path (the configured law). The indicator is discontinuous, so the rule locates the admitted region's edge and does not integrate across it.
- **BARS as the one-axis oracle.** The two are separate code paths, so their agreement is evidence.
  - *Evaluations.* For the same knots, BAMRS's target at d = 1 returns BARS's log prior and marginal: to rounding, and bit for bit where both call the same model on the same design.
  - *Posteriors.* On shared data the count law, the fit and `FitSd` agree with a BARS run within Monte Carlo error, under a gate whose size, seeds and rerun rule are fixed first.
- **Nothing built is disturbed.** Every commit that touches a shared file passes every test project from disk, with BARS's and BAPS's pinned tiers unchanged.
- **A full grid by both paths:** the Kronecker shortcut and the scattered path agree to rounding.
- **Axis exchange:** transposing the inputs leaves every posterior summary unchanged.
- **A Geweke joint test** of the wrapped moves, as E34's.
- **Seeded random sweeps** of the numerical primitives: the design builder, the Bernstein conversion and the subdivision solver.
- **Column order:** a rectangular design with different penalties on its two axes, against dense Kronecker algebra, in both size orders (item 3).
- **A fixed-knot fixture** for the spread: with the knots held, every draw's conditional mean is the same, and `FitSd` equals the closed-form posterior sd of the fit.
- **External bands:** Web Table 2, pinned to its script and units.
  - *The script* is `exampleABMRS/scripts/surface_fitting_appendix.r`: five smooth cases (`g1` to `g5`, lines 54–75), m = 250 and 1000, noise sds 0.25, 0.22, 0.22, 0.25 and 0.25. By the appendix's text the runs are cubic, with γ ∈ {0.8, 1, 1.2} and 50 replicates.
  - *The units.* The script stores each replicate's MSE against the noise-free function at m fresh uniform test points (line 163). `summary_output.qmd` takes the square root of each replicate's value (line 298) and reports the mean and sd over replicates (line 282). The table prints those in units of 10⁻¹. So the comparable number is a per-replicate RMSE, averaged.
  - *`surface_fitting.r` is a different experiment:* three cases at m = 500 and 1000 with noise sds 0.8, 2 and 2, two of them tensor splines with fourfold knots. By Web Table 3's three surfaces at m = 500 it is the main text's surface study. Web Table 3 holds its coverage, and its RMSE table is in the main text, which is not held. Its functions are recoverable from the script and serve as fixtures without a published band.
- **Coverage** of the mean function is reported, with the source's 66–96% as context.

## Not read, not verified

- The article's main text. The prose reasons for c = 0.3 and for [0.5, 2], and whether rows 6 and 7 were fixed in print, stay unknown. None blocks the build.
- The bodies of `BarsEnsemble` and `CollapsedGPrior`, and the contracts of `SmcRun` and the flat-histogram family.
- Chen and Chen's condition, quoted from memory.
- The references under Multiresolution, and the contour-tree literature for bilinear cells.
- Nothing in the Lean section was compiled. Its Mathlib names were checked for existence in the local checkout, and the fold's effect on the gated count was read against `SignificantPeaks` and `ValleyFloor`, not tested.
- Nothing here was run. The numbers under "What the product structure changes" are arithmetic, and item 6's count pairs are an enumeration.
- The fMRI analysis, and the article's use of the cardiopulmonary-bypass reference, are known from the abstract and the reference list alone.
- Whether `SplineBasis` and `SplineExtrema` serve a derivative spline at degree p − 1, as their signatures suggest. And which of BARS's later blocks take a design and not a carrier.
- The 2026-11-15 revision ran no sampler. The fold identity Δ′ = −4A·∂_θ∂_T f was derived by hand and checked on two examples. That p + 1 knots free the two sides of a gap is a dimension count, not a run.

## Doc changes owed

- **BAMRS.md:** the expansion to "Multivariate Regression"; the overview's "multi-scale hierarchical spline bases"; row 1's "exact marginal" to name the centred block; rows 1, 6, 7, 10 and 11 of the table above.
- **Regression.md,** line 14: the expansion and its gloss.
- **Sources&Changes.md,** the He, Yang and Kang row: the same new rows.
- **The two HYK-2024-2026 discussion files:** a header line marking them unverified and not to be cited, and noting that the fMRI application they mention is real.
- **The extensions brief's E3 entry:** superseded by this brief for design; its queue line stays.

## Review of 2026-11-15

An external read-only review, relayed by the user. Its main finding is that admission couples the axes, while the reference sampler and its acceptance test were written for independent counts. Each point was checked against the source before it was taken; all nine held.

| # | The review's point | Checked against | Where it changed this brief |
|---|---|---|---|
| 1 | The count law after admission is joint. The test and the reference sampler assumed independent counts, and the tempering path was unstated | `SplineTarget.cs:85`–`106`, `:148`–`207` | Item 6; the structure prior; the Target row; slice 1; Validation's quadrature |
| 2 | The tensor dimension is checked before any design is allocated | The enumeration: 1,359 of 1,002,001 pairs, and 3.76 GiB at (1000, 1000) | Admission's order; item 6 |
| 3 | A refusal share and a per-axis R̂ do not answer the gate's and the mixing question | A point of method | Decision 2; item 1 |
| 4 | A storage order that flips with the counts swaps the axes silently | `TensorPenalty.cs:101`–`105`; `TensorProductDesign.cs:13`–`17` | Item 3; the Tensor design row; Validation |
| 5 | Distinct knots give a continuous transition, so a jump has no size; the feature is a contrast across a gap, and its amplitude varies along the other axis | `mebars.cpp:22`–`33`, `:142`; `surface_fitting.r:61`, `:68` | Why; item 4; "Transitions and kinks"; Multiresolution |
| 6 | `FitSd` and the Sobol outputs need coefficient uncertainty | `BarsEnsemble.cs:573`–`577`; `SobolDecomposition` | "What carries over unchanged"; Validation |
| 7 | A discriminant root is a candidate fold, and monotonicity is a sign test | Both counterexamples by hand | "On a surface the fold becomes a point"; Monotonicity; Lean |
| 8 | The fMRI application is real, and four links resolved one directory too high | Crossref's record for the DOI; the folder listing | Sources |
| 9 | Web Table 2 belongs to `surface_fitting_appendix.r`, which stores MSE where the table reports RMSE | Both scripts; `summary_output.qmd:282`, `:298` | Validation's external bands |

Added in the same pass, beyond the review's text:
- An axis's marginal count mass after the size bound (121 pairs against 17), as the concrete form of point 1.
- The gate without a dense design, and the per-axis mass as a bound on the product column's.
- The tally label a reordered admission can change. The re-plan below made this a difference between two engines' tallies, where the first text had a one-axis path to protect.
- For point 4, the choice between the review's two remedies: a fixed coefficient order, with the permutation confined to the band.
- For point 5, that the source's own sampler draws distinct sites while its test surfaces carry fourfold knots, and that `StepBasis` is the discontinuous basis already at hand.
- For point 7, the identity Δ′ = −4A·∂_θ∂_T f, which makes "simple root" and "transversal fold" one condition.
- For point 8, the cardiopulmonary-bypass reference, which leaves the cardiac application unverified where the brief had called it invented.
- For point 9, the aggregation in `summary_output.qmd`: a square root per replicate, then the mean and sd.

## Rollout re-plan of 2026-11-15

The user's direction, given after the review above: unifying the engine while building BAMRS was a mistake. BARS has gone beyond the original Kass work and BAPS is fitted, and neither is to be dismantled or reshaped for BAMRS's sake. BAMRS is built as a new engine, as BASS was, reusing primitives and upgrading them along the way only where that breaks nothing already built.

What the survey behind the re-plan found (2026-11-15, by reading):
- Every engine already has its own ensemble on the generic shell in `src/maths/samplers/`: `BarsEnsemble`, `BapsEnsemble`, `BassEnsemble`, `BgpSampler`. BASS takes two primitives from BARS's folder (`TemperedForm`, `StudentTWeights`) and none of its target, moves, models or ensemble. BAMRS will call far more of BARS's public surface than that, as a consumer, which is why the reuse rule never asks those signatures to change.
- `src/maths/regression/spline/bamrs/` exists and is empty.
- The spline root already holds the neutral pieces: `KnotConfig`, `IBasis`, `SplineBasis`, `BSpline`, `BandedDesign`, `BandedGPrior`.
- BARS's models, moves, count priors and `SplineExtrema` are public and take a design or a `KnotConfig`, so BAMRS calls them without an edit. The whole tree is one assembly, `Maths.Regression`.

Where the brief changed:
- **Decision 1** is the user's: a new engine, with the reuse rule in three steps.
- **"What is already built"**: the Tensor design, Observation models, Target, Moves, Orchestration and Curve readouts rows now say what is left as it is and what BAMRS builds for itself.
- **Item 3**: the assembler and the coefficient form are BAMRS's own, and BAPS's builder and penalty are untouched. The clause that draws "pool without a permutation" is gone; draws at different counts never pooled.
- **Tier 1**: restrictions and derivative splines need no change to `SplineExtrema`. Only the oblique-line reader is new.
- **The curve features** are built with BAMRS, and their adoption by BARS is a separate decision.
- **Slices**: the ensemble is BAMRS's own (slice 3), the Gram path is BAMRS's own (slice 6), and what BARS grew arrives one slice at a time (slice 8).
- **Validation**: BARS is the one-axis oracle, and every commit on a shared file leaves BARS's and BAPS's pinned tiers unchanged.
