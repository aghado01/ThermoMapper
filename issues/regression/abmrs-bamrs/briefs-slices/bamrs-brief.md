# BAMRS: free-knot tensor-product spline surfaces (preparation brief)

Date: 2026-11-12. Status: preparation; nothing is built. This is E3 of the [extensions checklist](../../bass/extensions-checklist.md), queued by the user on 2026-11-08. It replaces the design notes under E3 in the [extensions brief](../../bass/extensions-brief.md#free-knot-tensor-product-surfaces-abmrs) as the working document. Written by the main session at the user's request, after a re-read of the sources and a survey of the engine. Its decisions are proposals until the user takes them.

## Why

BAMRS fits a surface over two or three inputs with a tensor-product B-spline whose knots are free on every axis: reversible-jump births, deaths and relocations per axis, with the coefficients and the noise variance integrated out. It is BARS on a product carrier. He, Yang and Kang's ABMRS is the published instance. BAMRS is this engine's version, on the corrected forms the regression family already holds.

What it adds to the family:

- **Smooth surfaces with exact per-draw geometry.** BASS's hinges default to degree 1, and BAPS's tensor splines keep fixed knots. A cubic tensor spline with free knots has exact derivatives, restrictions and extrema for each posterior draw, with the knot uncertainty carried.
- **Discontinuities along an axis,** from knots that cluster within the data's resolution.
- **Scope.** Low dimension: d = 2 for cubics, and d = 3 at degree 1 or with many observations (the support gate, below). The engine is application-agnostic.

**The name.** BAMRS is Bayesian Adaptive Multivariate Regression Splines, the source's ABMRS in the family's word order. [BAMRS.md](../../../../docs/Maths/Regression/BAMRS.md) and [Regression.md](../../../../docs/Maths/Regression/Regression.md) expand it as "Multiresolution", which is a correction owed. A multiresolution engine is a different model and a separate candidate ([its section](#multiresolution-a-separate-candidate)).

## Sources

- **The preprint:** He, Yang and Kang, arXiv 2405.13353v1 (May 2024), where the method is called EBARS. Markdown at [papers/HYK2024.md](../papers/HYK2024.md); the source is held under `codex-scientiae/supellex-private/gauntlet/bayesian/2405.13353v1`. arXiv has no later version (checked 2026-11-12).
- **The article:** *Biometrics* 82(3), 2026, ujag121, doi 10.1093/biomtc/ujag121. **The main text is not held** and the journal page refuses automated fetch. Every statement here about the published version is inferred from its open supplement.
- **The supplement:** [ABMRS_Revision_V2_Supplementary_Materials](../ABMRS_Revision_V2_Supplementary_Materials), with the web appendix as [text](../ABMRS_Revision_V2_Supplementary_Materials/web_appendix_extracted.txt). Web Appendix A proves Lemma 1 and Lemma 3 only, with the preprint's statements, which implies the article carries no theorem needing a proof beyond the preprint's lemmas.
- **The author's repositories,** cloned at `D:\aghado01\packages\bayesian-related`: `ABMRS` (its `src` and `R` are identical to the supplement's copy) and `exampleABMRS`, the scripts behind the article's simulations and its real-data example.
- **Chen and Chen (2008),** *Extended Bayesian information criteria for model selection with large model spaces*, Biometrika. Not held. Its consistency condition is quoted below from memory and wants checking against the paper before it enters `docs/`.
- **Checked earlier:** the report half of [abmrs-verification-questions.md](../../discussion/abmrs-verification-questions.md), which holds against the supplement.
- **Unreliable, and not to be cited:** [HYK-2024-2026-differences.md](../../discussion/HYK-2024-2026-differences.md) and [HYK-2024-2026-clarifications.md](../../discussion/HYK-2024-2026-clarifications.md). They are external summaries with invented content (a knot-multiplicity kernel, a gradient and Hessian peak finder, a changed death formula, fMRI and cardiac applications). [argon-abmrs-verification-debrief.md](../../discussion/argon-abmrs-verification-debrief.md) retracts part of it.

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
- **What the article adds to the preprint,** as far as the supplement shows: γ past 1, c = 0.3, the column filter, coefficient draws with prediction intervals, a peak-detection simulation, and a real-data example (triceps skinfold). The preprint's manifold-denoising application has no trace in `exampleABMRS`.

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
- **Structure prior.** Independent across axes. Its log density is Σₐ [ln p(kₐ) + ln kₐ!], the count law times the order-statistics density, restricted to the admissible set. Any `IComplexityPrior` serves an axis. `EbicPrior(γ, Cₐ)` reproduces the source's law, since a uniform draw of kₐ sites from a grid tends to the order statistics of kₐ uniforms as the grid refines; C stays as a nominal count.
- **Coefficients and noise.** The centred collapsed g-prior block, with the σ² priors BARS has (`Sigma2Prior`). Under the Jeffreys prior its marginal is Lemma 1 on ỹ = y − c.
- **Admission.** The support gate on product columns, then factor admission (`DesignAdmission`), then m > ν. A refused configuration scores −∞. Truncating the prior to the admissible set keeps every acceptance ratio valid, since an inadmissible proposal is simply rejected.
- **Moves.** Birth, death and relocation on one axis at a time, with the chain owning the selection ratio. Constant weights are BARS's default and the prior-cancelling schedule its option (`DmgkSchedule`); the source's schedule is that option at c = 0.3 on the EBIC law.

## What the product structure changes

These follow from the model. None is in the source.

1. **A knot on one axis costs a whole slice of columns.** A birth on axis a adds Δν = Π over a′ ≠ a of ν_a′ columns, and the marginal charges (m+1)^(−Δν/2) for them.
   - At m = 500 with 8 columns on the other axis, that is 25 nats: the knot must cut the residual sum by about 9.5% to break even. In one dimension it is 3.1 nats, about 1.2%.
   - At d = 3 with 6 columns on each other axis it is 112 nats.
   - *Consequence.* Per-axis count posteriors will be sharp and count changes rarely accepted, so relocation carries the mixing. BARS found that tempering rarely pays (BARS.md, "When to temper"); that finding was for one dimension and may not hold here.
   - *The measurement and its decision.* Start chains at different per-axis counts and read R̂ of each kₐ. If they disagree after the burn-in BARS uses, a ladder or the flat-histogram sampler over counts becomes BAMRS's default. If they agree, tempering stays optional as in BARS.
2. **The support gate bites at the corners.** A clamped boundary B-spline of degree p on a first interval of width h has integral h/(p+1). For inputs roughly uniform on the box, a corner column's mass is about m·Πₐ hₐ/(pₐ+1).
   - Under the default of one point's worth, cubics at d = 2 need h₁h₂ ≥ 16/m at each corner: at m = 500, first knots at 0.1 on both axes are refused (mass 0.31), and 0.2 by 0.4 is admitted (2.5).
   - At d = 3 the divisor is 64 for cubics and 8 for degree 1. Cubics with m = 1000 and intervals of 0.3 are refused (0.42).
   - The source's filter is harsher in count and silent in effect: it drops those columns and keeps the knots.
   - *The measurement and its decision* is under Decision 2.
3. **Size.** ν grows as a product and the block needs m > ν. In lexicographic order the Gram's half-bandwidth is about p·ν_inner + p, so the axis with fewer columns goes innermost. That choice is a pure function of the state, so it does not disturb bit-for-bit resumption.
4. **A discontinuity spans the box.** p + 1 knots clustered inside one gap in the data make a jump, and every column keeps its support: no B-spline has all p + 2 of its knots inside a cluster of p + 1, so each reaches data on one side. I argued this from B-spline support and did not run it against the gate. In a tensor product the jump is a whole line (or plane) across the box. A discontinuity confined to part of the domain has no compact representation in this model.
5. **Expect the mean function's intervals to under-cover.** The supplement reports 66–96% coverage at a nominal 95% (Web Tables 2c and 3), below 90% in most cells, against 93–95% for prediction intervals. My hypothesis is that sharp count posteriors (item 1) leave approximation bias outside the band, with their plug-in σ̂ a smaller part. It is a hypothesis; BAMRS's own coverage is a reported number, not a gate.

## What is already built

Read from the source on 2026-11-12. "As is" means by reading, not by running.

| Piece | Where | Use for BAMRS |
|---|---|---|
| Per-axis carrier and basis | `KnotConfig`, `SplineBasis`, `BSpline`, `StepBasis`, `IBasis` | As is, one per axis. `IBasis` needs no change: the product is taken over the axes' own designs, so axes may mix bases (a step axis with a cubic axis) |
| Tensor design | `Baps.TensorProductDesign` | Upgrade: it is bivariate and takes `SplineBasis`. Lift it to d axes over `IBasis`, and BAPS's `AnisotropicTensorSpline` and BAMRS share one builder |
| Face-splitting contraction (N-X1) | `KhatriRao.FaceSplittingGram`, `FaceSplittingCrossProducts`, `FaceSplittingHalfBandwidth`, `KroneckerGram`, `FaceSplittingRow` | As is. Built for this; only tests call it |
| Wide-band factor (N-X2) | `WideBandCholesky`, `SpdFactors.ForBand` | As is. Built for this; nothing selects it |
| Observation models | `IObservationModel`: `WeightedNormalModel`, `PoissonModel`, `WhittleModel`, `GibbsLossModel` | As is in correctness: each takes a dense design and knows nothing of the basis. `BandedDesign` reads each row's first and last nonzero column, so a tensor row's span is p·ν_inner + p wide and the accumulation costs O(m·bw²). At the source's benchmark sizes (m ≤ 1000, ν ≲ 150) that is small. The upgrade, when sizes ask for it, is accumulation by `FaceSplittingRow` and a factor from `SpdFactors.ForBand` |
| Collapsed block | `CollapsedGPrior`, `BandedGPrior`, `GPriorCentre`, `Sigma2Prior` | As is: any `ISpdFactor`, the centre conventions, tempered marginals, leave-out densities |
| Target | `SplineTarget` | New beside it, same composition: per-axis count priors, the gate on product columns, the model's marginal, and a reference sampler that draws each axis independently |
| Moves | `KnotBirthMove`, `KnotDeathMove`, `KnotRelocateMove`, `KnotBirthDensity`, the kernels, `DmgkSchedule`, `ScheduledMove` | Lift: an adaptor applies an `IRjMove<KnotConfig>` to one axis of the product carrier. The chain normalizes selection over the 3d moves, and a null proposal at kₐ = 0 stays put, which is the balance row 5 asks for |
| Chains and ensembles | `ReversibleJumpChain<TState>`, `ParallelTempering<TState>`, `TemperingLadder`, `ChainEnsemble`, `IEnsembleModel<TState>`, `LadderEvidence` | As is: generic in the state |
| Other samplers | `SmcRun`, `MulticanonicalRun` and the flat-histogram family | Candidates for item 1's count barriers. Their contracts are unread |
| Orchestration | `BarsEnsemble : IEnsembleModel<KnotConfig>` (82 KB): chains, robust weights, noise, checkpoints, evidence, the conformal hooks | Decision 1. Its carrier and its peak readout are one-dimensional |
| Curve readouts | `SplineExtrema`, `BarsDraw`, `PeakDraw`, `PeakSpan`, λ, π, `PeakModes` | [Readouts](#readouts-how-the-bars-contract-extends) |
| Moments over piecewise polynomials | `Bass.Sensitivity.PiecewisePolynomial`, `BassMoments` | Sobol indices and effects for tensor B-splines, as the extensions brief notes |
| Declared coefficient priors | `PriorPencil`, `CollapsedPencilPrior`, `TauPrior`, BAPS's `TensorPenalty` | Decision 2's ridge option |
| Term contract (F-X3) | `Bass.IDesignKernel<TPoints>`, `DesignFactor` | The other route: a tensor B-spline part as a BASS term type. Not recommended as the first build (Decision 1) |
| Lean | [BARS.lean](../../../../lean/Enthymemata/BARS.lean): MP-1, MP-2, SP-1 to SP-3 stated, proofs open; protolemmata MP-3, MP-4, O1–O3; [Stability.lean](../../../../lean/Enthymemata/Stability.lean) | [Lean](#lean) |

## Decisions

Each is the main session's recommendation, with its alternatives. The user may take or reverse any.

1. **One engine over the carrier, with BARS as its d = 1 member.**
   - *Recommended.* Build the product carrier, the lifted moves and the tensor target beside the existing types, and qualify them on the generic chain and ensemble. Then lift `BarsEnsemble`'s orchestration over the carrier, with the per-draw readout behind a seam. BARS becomes the one-axis case. Deletions wait for the final sweep.
   - *Acceptance.* BARS's pinned tier passes unchanged through the lifted engine, so the one-axis path must consume its streams exactly as today.
   - *Alternative: a sibling ensemble.* Faster to a first surface, at the cost of a second copy of the orchestration (robust weights, noise, checkpoints, evidence, conformal hooks) to keep in step.
   - *Alternative: a BASS term type.* BASS sums many sparse terms under its term measure; ABMRS is one full tensor term with moving knots. The fit is poor, though a tensor part inside BASS remains a later option through F-X3.
   - *Open until slice 0:* how much of `BarsEnsemble` depends on the carrier. I read its surface, not its body.
2. **Thin product columns: the g-prior with the gate, and a ridge option only on evidence.**
   - *Recommended.* Keep the g-prior and the gate, the source's model in its corrected form, as the default.
   - *The measurement and its decision.* On the source's benchmark surfaces at d = 2, count the share of evaluations refused as `ThinSupport` and compare the admitted counts with the counts the source reports. If the gate removes configurations the posterior would hold, add a declared ridge or penalty prior through `PriorPencil`, under which thin columns are regularized and the gate can relax. If it does not, the option waits for a caller.
   - *Alternative: lower the gate's default with d.* It trades refusals for coefficients the g-prior leaves almost unconstrained.
3. **The count law's default is open.** `EbicPrior(1, 2m)` matches the source's runs and gives comparable tables. BARS's `PoissonPrior(5)` per axis implies about 81 columns a priori for a bicubic. Take the default after the benchmark, not before.
4. **Fixed per-axis counts are a mode.** A relocation-only palette at given counts. The article's knot intervals and its peak simulation both use it, and it is cheap.
5. **Readouts arrive in tiers** (next section): affine and restriction readouts first, degree-1 lattice extrema second, joint cubic extrema last.
6. **The full-grid shortcut is its own slice.** On a complete grid Z = Z₁ ⊗ Z₂ and ỹᵀP ỹ = ‖Q₁ᵀ Ỹ Q₂‖²_F, so a move on one axis is a one-dimensional refit. Missing cells and robust weights that do not factor across axes both break it, and the scattered path then serves.
7. **Non-normal responses wait** until the Normal model qualifies. `PoissonModel` takes a tensor design as it takes any other.

## Readouts: how the BARS contract extends

**What BARS has.** Each retained draw holds its knots and one joint draw of (β, σ²). `SplineExtrema` reduces that draw exactly: the candidate set is the endpoints, the knots and each span's derivative roots, in closed form for degree ≤ 3. Peaks pass a prominence gate, and each peak's span is the level crossing at a drop fraction of its prominence, again as a closed-form root. The pooled views are λ(T), π(T) and `PeakModes`. This is a cubic method, and it is a method for one variable.

**What carries over unchanged.**
- The rule: readouts affine in the surface use the conditional mean, and nonlinear ones use the joint draw.
- The order: reduce each draw exactly, then pool. MP-2 and SP-3 are statements about averaging and do not depend on the dimension.
- The affine readouts: the fit, `FitSd`, partial derivatives, integrals over an axis, and the Sobol moments.

**Tier 1: restrictions reuse the curve engine as it stands.** Fix every axis but one. f(·, θ) is a spline in T on axis 1's knots, with coefficients cⱼ(θ) = Σₗ βⱼₗ Bₗ(θ). `SplineExtrema` applies to it unchanged, so each draw gives the ridge T*(θ), its height, its significant peaks and their spans along T, exactly, at any θ. Pooled over a θ grid these are λ(T | θ) and π(T | θ), and a posterior band for the ridge.
- *One upgrade that serves both engines.* `SplineExtrema` takes (`KnotConfig`, coefficients, `SplineBasis`), though inside it needs only breakpoints, an evaluator and a degree. Opening it to a piecewise-polynomial input lets it read restrictions, derivative splines (below) and oblique lines. Along an oblique line a bicubic restricts to a piecewise polynomial of degree 6, which the general-degree path already handles.
- *The ridge is a peak set at each θ.* T*(θ) jumps where two ridges exchange dominance. So the pooled readout is the peak set along T at each θ, λ(T | θ), and the argmax is one view of it.
- *Ridge births and ends belong to this tier.* They are roots of one polynomial in θ (see the fold on a surface, below).

**Tier 2: degree-1 surfaces have a closed form.** On each cell a degree-1 tensor spline is affine in each coordinate. A function affine in each coordinate attains its maximum over a box at a vertex, so every strict local extremum of the surface sits on the knot lattice, the boundary lattice included. A weak extremum can fill a line, as a ridge of constant height along a cell edge does. A bilinear cell's one interior critical point is a saddle. So at degree 1 peak inference is knot inference, which is the source's simulation, here with coefficient draws and the boundary.

**Tier 3: joint extrema of a cubic surface are new mathematics.** Fermat's candidate set stratifies over the cells: lattice vertices; critical points of the cubic restricted to each cell edge, which the existing closed form finds; and interior solutions of ∂₁f = ∂₂f = 0. The interior system has bidegrees (2, 3) and (3, 2), so at most 2·2 + 3·3 = 13 isolated solutions per cell. A cell can also hold a curve of critical points, as −(u − v)² does. That is the surface's counterpart of a flat span, and it needs its own result.
- *Recommended solver.* Convert each cell's patch to Bernstein form and subdivide. The control net bounds the patch, so cells are excluded with a certificate, and the global argmax needs only branch and bound on that bound.
- *Alternatives.* A resultant reduces the system to one polynomial of degree 13, which is ill-conditioned in the monomial basis. Newton from seeds gives no completeness.
- *What changes in the contract.* "Exact, in closed form" becomes "enclosed to a stated tolerance". The guarantee is one-sided: no critical point lies outside the returned boxes, and that a box holds exactly one is a separate test. The `argmax_in_closed_form_set` specification weakens accordingly for cubic surfaces, and holds as stated at degree 1.
- The same solver finds saddles, so it is one tool for two readouts. The fold points of ridges need only one-variable roots (below).

**Significance in two variables.** Prominence generalizes as the persistence of a maximum in the superlevel-set merge tree, and a span as the superlevel component at the drop level: a region, with an area and axis extents. Both need saddles. Obligation O3 already assigns the merge and col definitions to the planned carrier-generic π₀ and merge-tree foundation; the surface readout should consume that foundation too, and wait for it.

## Features on spline curves beyond peaks

This answers the user's wish to characterize more kinds of feature. It is a proposal, and it applies to BARS today and to every Tier 1 restriction of a surface.

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
- **Monotonicity on [a, b]** has a posterior probability: the share of draws with no interior critical point there.
- **Jumps and kinks.** A cluster of p + 1 knots inside one data gap is a jump, with size f(c⁺) − f(c⁻) per draw. At degree 1 each knot is a slope change, which is the quantity the source's triceps tests address by plug-in.
- **The topological summary** is the curve's 0-dimensional persistence diagram, exact from the candidate set. Prominence is persistence under the elder rule, up to the boundary convention of O1, and `Stability.lean` already pushes a sup-norm band forward to it.

**Shoulders and the fold.** A shoulder is a local minimum of |f′| at a nonzero value: a root of f″ where f′ keeps its sign. On a cubic span with f′(u) = Au² + Bu + C and discriminant Δ = B² − 4AC, taken where the inflection u* = −B/(2A) lies inside the span:
- f′(u*) = −Δ/(4A), the shoulder's residual slope;
- for Δ > 0 a maximum and a minimum straddle the inflection, separated by √Δ/|A|, with height difference Δ^(3/2)/(6A²);
- for Δ < 0 there is no critical pair, and |Δ/(4A)| measures how far the span is from having one.

So one closed-form scalar per span orders a strong peak-and-valley pair, a marginal pair, a shoulder and a plain monotone stretch. Its zero set is the fold in MP-3's catastrophe set ("a critical pair annihilates: the per-span quadratic discriminant crosses 0"). A shoulder and a marginal peak are the two sides of that fold. The newborn pair's prominence is the height difference above, which starts at zero, so a positive prominence gate is passed beyond the fold, where Δ^(3/2)/(6A²) reaches the gate. One feature then has two posterior probabilities: that the pair exists (Δ > 0), and that its peak is significant (the gate is passed). The fold moves the raw count and leaves the gated count alone. In general, where a pair straddles a knot, the definition goes through the derivative spline and not the per-span discriminant.

**On a surface the fold becomes a point.** In a family f(T, θ), a ridge is born or ends where ∂_T f = 0 and ∂²_T f = 0 together: the parameter value at which a shoulder becomes a peak. On a bicubic cell ∂_T f is a quadratic in the span coordinate whose coefficients A(θ), B(θ) and C(θ) are cubics in θ. So the ridge branches are the quadratic formula in θ, and the fold points are the roots of Δ(θ) = B(θ)² − 4A(θ)C(θ), one polynomial of degree at most 6, kept where the inflection lies inside the span. The root finder behind `SplineExtrema`'s general-degree path covers that degree, so ridge births and ends are a Tier 1 readout.

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
- It matches the upgrade proposed for `SplineExtrema`, which needs breakpoints, an evaluator and a degree.
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

**Surfaces.**
- *Restriction.* Σⱼ Σₗ βⱼₗ Bⱼ(x) Cₗ(θ) = Σⱼ (Σₗ βⱼₗ Cₗ(θ)) Bⱼ(x) for any two families of functions. It is an exchange of sums and needs no spline theory. It reduces Tier 1's exactness to O1 and O2.
- *The ridge in closed form.* On a bicubic cell ∂_T f = A(θ)u² + B(θ)u + C(θ) with cubic coefficients. The ridge branches are the quadratic formula in θ, and the folds are the roots of Δ(θ), of degree at most 6. Where Δ > 0 the ridge's slope is −∂_θ∂_T f / ∂²_T f. The regularity that slope needs is Δ ≠ 0, the fold's own predicate.
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

- **What it would fix.** Three of the five items above: a refinement adds a few columns where the data ask, not a slice (item 1); refinement stays where the data are, which eases the corner rule (item 2); and a discontinuity can be local (item 4). The preprint's own closing remark points the same way: the domain must be a box, and "a flexible multivariate spline" in place of the tensor product is its suggested remedy.
- **What it would reuse.** Each basis function is still a product of one-variable B-splines, so the face-splitting rows, the collapsed block, the reversible-jump chain and the readouts all serve. Its moves refine and coarsen cells, closer to BASS's term births and deaths than to knot moves.
- **What it would cost.** The Gram loses its band, the basis needs a linear-independence rule (truncation supplies one), and the structure prior becomes a prior over refinement trees.
- **A one-variable member exists too:** nested dyadic knots on a curve.
- **Naming.** Both models would abbreviate with M and R. BAMRS stays with the multivariate engine, to match its source; the multiresolution one takes a name when it is real.

## Slices

Proposed order; engine first, with no consumer wired in.

0. **Read and measure.** Read `BarsEnsemble`'s body for Decision 1. Run a bicubic design through `WeightedNormalModel` as it stands, and check the marginal against dense algebra.
1. **Carrier, design, target, moves.** The product carrier, the d-axis design builder, the tensor target and the axis adaptor. Qualify on the generic chain: each axis's count law under a constant likelihood, kₐ = 0 included, and the reference sampler.
2. **The gate's measurement** (Decision 2) and **the count-mixing measurement** (item 1), with their decisions.
3. **The ensemble** (Decision 1), with BARS's pinned tier as the acceptance.
4. **Tier 1 readouts,** with `SplineExtrema` opened to piecewise polynomials, the order-2 curve features, and ridge births and ends.
5. **Fixed counts and Tier 2.**
6. **Accumulation and factor upgrades,** when a measured size asks for them.
7. **The full-grid shortcut.**
8. **Tier 3 and the merge tree,** after the π₀ foundation.
9. **Non-normal responses.**

## Validation

Ground truth comes from outside the estimator, and the source's code is a sanity band, never a parity target.

- **Dense algebra** for the marginal, the centre and the conditional draws on small tensor designs.
- **Quadrature** for the posterior over counts when each axis holds at most one knot: the evidence of (0, 0), (1, 0), (0, 1) and (1, 1) is an integral over at most two knot positions.
- **The one-axis case equals BARS** bit for bit.
- **A full grid by both paths:** the Kronecker shortcut and the scattered path agree to rounding.
- **Axis exchange:** transposing the inputs leaves every posterior summary unchanged.
- **A Geweke joint test** of the lifted moves, as E34's.
- **Seeded random sweeps** of the numerical primitives: the design builder, the Bernstein conversion and the subdivision solver.
- **External bands:** the RMSEs of Web Table 2, and the surfaces in `exampleABMRS/scripts/surface_fitting.r`, whose functions are recoverable from the script.
- **Coverage** of the mean function is reported, with the source's 66–96% as context.

## Not read, not verified

- The article's main text. The prose reasons for c = 0.3 and for [0.5, 2], and whether rows 6 and 7 were fixed in print, stay unknown. None blocks the build.
- The bodies of `BarsEnsemble` and `CollapsedGPrior`, and the contracts of `SmcRun` and the flat-histogram family.
- Chen and Chen's condition, quoted from memory.
- The references under Multiresolution, and the contour-tree literature for bilinear cells.
- Nothing in the Lean section was compiled. Its Mathlib names were checked for existence in the local checkout, and the fold's effect on the gated count was read against `SignificantPeaks` and `ValleyFloor`, not tested.
- Nothing here was run. The numbers under "What the product structure changes" are arithmetic.

## Doc changes owed

- **BAMRS.md:** the expansion to "Multivariate Regression"; the overview's "multi-scale hierarchical spline bases"; row 1's "exact marginal" to name the centred block; rows 1, 6, 7, 10 and 11 of the table above.
- **Regression.md,** line 14: the expansion and its gloss.
- **Sources&Changes.md,** the He, Yang and Kang row: the same new rows.
- **The two HYK-2024-2026 discussion files:** a header line marking them retracted.
- **The extensions brief's E3 entry:** superseded by this brief for design; its queue line stays.
