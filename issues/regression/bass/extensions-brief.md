# BASS extensions brief

Extensions beyond the port's completion contract, collected from the [full port brief](../bass-full-port-brief.md) and the reviews behind it. Each entry says what it adds, its source, what it needs from the port (each segment's "Room for extensions" section), and its status. Nothing here is scheduled: an entry starts on the user's go-ahead, with its own specification and acceptance criteria. Its checklist, where the reports go, is [extensions-checklist.md](extensions-checklist.md). Gaps past the lineage's own work, some of which bear on these entries, are in the [lineage-gaps register](lineage-gaps-register.md). The [lineage review](lineage-review.md) is the map behind entries E19–E30: the lineage's axes, what its closed forms rest on, and a proposed build order.

| Extension | Source | Needs from the port | Status |
|---|---|---|---|
| Generalized-hyperbolic noise (GBMARS) | Rumsey, Francom and Shen, arXiv 2306.01911 | Fitting: the offset seam, the latent slot, per-observation densities. Numerics: conjugate updates as limits. Outputs: the declared surface | Deferred |
| Robust centring and bases (OQ16) | Gervini 2008; Locantore et al. 1999 | Outputs: the basis contract | Deferred |
| Free-knot tensor-product surfaces (ABMRS) | He, Yang and Kang, *Biometrics* 2026 | Numerics: the face-splitting contraction, a wide-band factor. Fitting: the design-kernel contract | Parked |
| Additional proposal kernels | KHAOS, arXiv 2510.25036; Nott and Kohn 2005 | Fitting: the proposal-kernel contract. Numerics: D10's policy | Candidate |
| KHAOS's modified g-prior | Rumsey, Francom et al., arXiv 2510.25036 | Fitting: the coefficient-prior set | Candidate |
| Basis-estimation uncertainty | Sartini et al., arXiv 2412.11340 | Outputs: the basis contract's uncertainty slot | Candidate |
| Moments for piecewise polynomials | Derivations review §7, item 7 | Outputs: the moment service | Candidate |
| Active-subspace analysis | arXiv 2307.11241, via the PCA-extensions discussion | Outputs: term derivatives | Candidate |
| Cross-engine evidence | Derivations review §7, item 8; Gneiting and Raftery 2007 for held-out scores; Liang et al. 2008 and Farhi et al. 2013 for Bayes factors under the improper options | D03's proper defaults and D05's ladder evidence; D03's evidence rule. Outputs: D07's joint draws | Candidate |
| The free-energy curve over the ladder | Farhi et al. 2013; Watanabe 2013; Bissiri, Holmes and Walker 2016 | D05's ladder with the reference rung, its rungs' draws kept; the integrated form's per-draw slope | Candidate |
| Elastic emulation of misaligned functional outputs | Francom, Tucker, Huerta, Shuler and Ries, arXiv 2305.08834 | Outputs: a representation with a declared reconstruction map. Calibration: D09's likelihood over both parts | Candidate |
| Residuals correlated along the output index | Andros, Guhaniyogi, Francom and Pasqualini, arXiv 2406.18751, 2602.09351 | Fitting: the output-side inner product; the residual-dependence check | Candidate |
| Semi-modular calibration | Carmona and Nicholls 2020, arXiv 2003.06804; Bissiri, Holmes and Walker 2016, §4.2 | Calibration: D09's cut posterior and its emulator interface. Stage 8: the cut SMC's induced-density bridges | Candidate |
| Calibrated intervals for loss-based posteriors | Bissiri, Holmes and Walker 2016, §3.4; Holmes and Walker 2017; Lyddon, Holmes and Walker 2019; Syring and Martin 2019 | BARS: `GibbsLossModel`. Any engine that weights a likelihood or a loss | Candidate |
| Conformal Bayes, the consumer | Fong and Holmes 2021, arXiv 2106.06137; Vovk, Shen, Manokhin and Xie 2019 | The pointwise seam (OQ17, built): retained draws, the n × S matrix, the unit and stratum maps | Candidate |
| Conformal bands for whole curves | Lei, Rinaldo and Wasserman 2015; Diquigiovanni, Fontana and Vantini, arXiv 2106.01792 | The unit map's whole-curve unit; each draw's pointwise predictive sd | Candidate |
| Posteriors for readouts by resampling units | Huggins and Miller 2023; Fong, Holmes and Walker 2023 | Retained draws with their readouts (BARS's peaks); replicate units | Candidate |
| Polynomial-chaos terms (E19) | Rumsey, Francom et al., arXiv 2510.25036 | Fitting: the design-kernel and proposal-kernel contracts. E4, E5 | Candidate |
| Data-augmented observation models (E20) | Albert and Chib 1993; Polson, Scott and Windle 2013; Lawrence et al. 2008 | Fitting: the offset seam, the latent slot, a weighted Gram, σ² held. Distributions: a truncated draw (built) | Candidate |
| A prior-precision builder, with ridge blocks (E21) | KHAOS, via SC-3 | Fitting: the coefficient-prior set | Candidate |
| Correlated outputs and coupling within a run (E22) | Chib and Greenberg 1998; Ghosh and Dunson 2009; Hoff 2007 | Fitting: the output-side inner product, its Q-centred forms. E20. Samplers: a covariance sampler | Candidate |
| Multi-fidelity emulation through the offset (E23) | Kennedy and O'Hagan 2000; Tucker et al., arXiv 2609.27174 | Fitting: the offset seam. Stage 8's cut | Candidate |
| Guarded estimators (E24) | [GE-1 to GE-3](../../../lean/Protolemmata/BASS/guarded-estimator-lemmas.md); Angelopoulos, Duchi and Zrnic 2023; Giles 2008; Zrnic and Candès 2024 | Per-draw integrals under a declared measure, predictions, fitting on folds | Candidate |
| Sensitivity under dependent inputs (E25) | [G4](lineage-gaps-register.md#g4-one-moment-service-across-dependent-inputs); Francom and Nachtsheim, arXiv 2506.11471; Apley and Zhu 2020 | Outputs: `InputMeasures`, the moment service, term derivatives | Candidate |
| Stochastic emulation through latent inputs (E26) | Zhu and Sudret, arXiv 2202.03344 | Fitting: the latent slot, the design-kernel contract | Candidate |
| External benchmarking (E27) | Rumsey, Gibson, Francom and Morris, arXiv 2512.09060 | Predictive draws; E9's proper scores | Candidate |
| Bases that mix term types (E28) | None known | Fitting: the design-kernel contract, E21. E19 | Lead |
| Categorical levels: learned positions and a similarity readout (E29) | Zhang, Tao, Chen and Apley 2020 | Fitting: the design-kernel contract. Retained draws | Lead |
| Choosing the next runs (E30) | Binois et al. 2019; Ranjan, Bingham and Michailidis 2008; Watson and Pelli 1983 | Predictive draws, a noise model that varies with x | Candidate |
| Tree ensembles, with BCF's multiplier (E31) | Chipman, George and McCulloch 2010; Hahn, Murray and Carvalho 2020; [the notes](../perplexity-bart-etc.md) | Samplers: `BlockScan`, `IRjMove`, the target's separable case. A leaf contract beside `IObservationModel`. Readouts: declared smoothness | Candidate |
| Additive models with a GP component (E32) | [The notes](../perplexity-bart-etc.md) | BGP: `y` per fit, σ² from outside with s² free, a declared mean. Fitting: the marginal given σ² | Candidate |
| Regression-family integration (E33) | [The notes](../perplexity-bart-etc.md) | Every family: an inventory, declared capabilities, four verdicts | Candidate |
| Direct-functional BASS for families of SPC sweeps | Derivations review §7, item 9 | The completed port | Application |

## Upgrades to existing code

Existing code that entries below would otherwise each change on their own, with the ideas that call for it, read against the code on 2026-11-02. These can be built whenever the user chooses, without waiting for an entry to start (the user, 2026-11-02). An analysis of an entry counts them as part of its work, not as reasons it cannot be done. The [capability inventory](../capability-inventory.md) holds the evidence for each row, family by family, and the gaps that are true of the model itself.

| Upgrade | Today | Called for by |
|---|---|---|
| Weighted Gram, centring and marginal in the BASS core | **Built 2026-11-02.** `BassDesign` carries `WorkingResponse.Weights`, with both priors: D03's on the unweighted columns (`PriorPencil`, any (τ, b) in O(p)) and one on the weighted columns, where integer weights equal duplicated rows ([the forms](../../../docs/Maths/BayesianRegression.md#bass-the-weighted-core), [the report](extensions-checklist.md#reports)). BARS's own weighted block is still built at one τ and b = 1; its tempering could read the pencil | E1, E20, E22 through E20, E31's multiplier; BARS's weighted tempering |
| The structure step's marginal given σ² | **Built 2026-11-02.** `GPriorMarginal.LogMarginalGiven`, selected by `Sigma2Prior.Held` or `BassTarget.Given`. BGP's known-σ² mode has one too, at s² = 1. `BassTauBlock` and a chain's composition read the target the chain holds, so a sampled τ can sit beside a moving σ² | E1's mode with σ² held, E20's probit, E32 with BASS as a component |
| Priors that stop reading a moving response | **Steps 1 and 2 of [the options](moving-response-priors-options.md) built 2026-11-02.** `BassPriorReference` on `WorkingResponse` freezes the level, s², the coefficient prior's constant and a Gram on W₀ (`WeightedPriorGram.Frozen`), with a named rule per noise; `Sigma2Prior.Hierarchical` and `HalfCauchy` put σ²'s rate on the state ([the forms](../../../docs/Maths/BayesianRegression.md#bass-priors-that-stop-reading-a-moving-response), [the report](extensions-checklist.md#reports)). Steps 3 and 4 wait on the latent blocks and the run. BARS's `WeightedNormalModel` re-reads at every evaluation, over a fixed y. BAPS's component prior and the hierarchical half-t scale read it once | E20, E32, E23 |
| A BASS design whose working response moves | **Built 2026-11-02.** `WorkingResponse.Moved` and `BassStructureMoves.Follow` move a chain between sweeps: the design refits the current structure, and the state's fit names the working response it was fitted on, which the target reads, so swaps and particles carry it. A move that leaves F is refused, a model whose priors read the snapshot refuses to move, and a ladder's and a run's composition read the target they hold. Pointwise, the training fit and the residual check read each draw's own ([the forms and the options](../../../docs/Maths/BayesianRegression.md#bass-a-design-whose-working-response-moves), [the report](extensions-checklist.md#reports)). F reading W is open to the user's review, since BARS took the other reading. Latent blocks, a latent's draw at b = 0 and `Bass.Fit` are the consumers' | E20, E23, E32, E1 |
| Truncated draws | **Built 2026-11-02.** `SampleOn(rng, lower, upper)` on every `IContinuousDistribution`, with `Truncated<T>.Sample`, by one [tangent hat](../../../docs/Maths/ProbabilityAndInformation.md#truncated-draws-and-quantiles) in a coordinate where each family is log-concave, reading no tail. A closed half-line path for the normal is the open cost question, about 3× Robert's per draw, if E20's latent block profiles there | E20's probit and discrete-time survival; F-R1's [`s2.lower` arm](deterministic-profile-brief.md#for-a-rerun) |
| Quantiles for the rest | **Built 2026-11-02.** `Quantile(p)`, `InverseLogCdf(logP)` and `InverseLogSurvival(logS)` on every `IContinuousDistribution`, Newton on the smaller log tail, checked against base R. The draws did not need them: [inversion](../../../docs/Maths/ProbabilityAndInformation.md#truncated-draws-and-quantiles) costs 2–9× the hat and loses a narrow interval | A complete family |
| `GpRegression` on a changing response | `y` is fixed at construction, and a known σ² forces s² = 1. The squared distances do not read `y`, so a `y` per fit is cheap. The nugget is one scalar, with no weights | E32 |
| Declared complete values | **Built,** `4e3ea95`: `ISamplingTarget.HasCompleteValues`, which `LadderEvidence` and `Reweighting.LogEvidence` read. The slope is still probed by catching an exception | [Cross-engine evidence](#cross-engine-evidence), E33 |
| BARS's integrated-form slope in b | `NormalMarginal` overrides no slope, so BARS's path sampling is linear-form only, though `GPriorMarginal.LogMarginalSlope` exists and BARS's τ is fixed, which keeps the slope at b = 0 finite | E13, E9 |
| Path sampling on BASS's integrated form | The target supplies the slope, and `BassLadder` traces it on the linear form only, since the slope at b = 0 carries E[1/τ], infinite under Zellner–Siow | E13 |
| Predictive draws outside BASS | Only BASS draws new responses. BARS, BAPS, the hierarchical model and BGP keep retained draws or a mean and variance, and BGP never draws f | E27, E9's held-out scores, E16, E17, E30 |
| Derivatives of a fit | BASS has term-level derivatives and moments. `IBasis`, `BSpline` and `IGpKernel` have none, and no family reads out a fit's derivative | E8, E25 |
| Declared smoothness | No contract states a fit's smoothness. BASS's default `Degree = 1` has kinks; BARS's cubic fit is C² at simple knots; a squared-exponential GP is C^∞ | E31, E8, E25, E33 |
| A power on the likelihood | BAPS, the hierarchical model and BGP take no inverse temperature | E15; E9 through BGP's evidence |
| Evidence outside Normal BARS and BASS | BGP: a flat mean prior, 1/s², an unnormalised bandwidth prior, no tempering. BARS's Poisson and Whittle drop a known Σ ln yᵢ!. BAPS and the hierarchical model cannot: their penalty is improper on its null space | E9 |
| Exposed prior draws of the full state | BARS draws its knots and BASS its (structure, τ) exactly. No method composes (β, σ²) at b = 0 onto them, and the repo has no Geweke test | E31, E33 |
| Weights and offset in BARS | Poisson, Whittle and Gibbs loss ignore `weights` without a word, and `Bars.Run` takes none | E1's BARS side, E20, E33 |
| Weights, offset and a response per sweep in BAPS and the hierarchical model | y, ZᵀZ and Zᵀy are fixed at construction; `BandedDesign.Accumulate` already takes weights | E15, E20's seam applied to BAPS, E33 |

## Generalized-hyperbolic noise (GBMARS)

**Deferred extension (2026-10-25): generalized-hyperbolic noise (GBMARS).** Rumsey, Francom and Shen (arXiv 2306.01911) extend BASS's noise to a normal variance-mean mixture, εᵢ = √w(βvᵢ + √vᵢ·ζᵢ) with ζᵢ ~ N(0, c) and vᵢ ~ GIG(p, a, b). That covers Student t, Laplace, asymmetric Laplace (quantile regression), normal-inverse Gaussian and variance-gamma noise. It lies outside the completion contract, and generalizes the derivations review's §7 item 4 (robust BASS from BARS's resamplers).
- **BARS already has the β = 0 half.** [`IWeightResampler`](../../../src/maths/regression/spline/bars/WeightResampling.cs) is a scale mixture with Gamma precision weights, and `StudentTWeights`' conditional is GBMARS's GIG conditional at β = 0 and a = 0. Its sweep is GBMARS's.
- **What the seam needs:**
  - a GIG distribution and sampler in `Maths.Distributions` (Hörmann and Leydold 2014, the bounded rejection sampler GBMARS uses; Devroye 2014), and a GIG conjugate update in `Maths.Samplers.Blocks` whose limits b → 0 and a → 0 are `ConjugateUpdates`' Gamma and inverse-gamma updates;
  - a mean offset: a variance-mean mixture hands the Gaussian core the working response y − β√w·v along with its weights, where a scale mixture hands only the weights;
  - log Bessel K of real order, for the pointwise density, which integrates each vᵢ because it belongs to one observation (OQ17's scope rule). Nothing in the repository computes it yet.

  `ContaminatedNormalWeights` is a discrete mixture and stays its own implementation.
- **Sweep and prior:** D04's robust row and D03's unweighted g-prior, in the [fitting segment](port-fitting-brief.md). GBMARS's own choice, a ridge with a sampled τ, stays selectable. GBMARS calls its λ update approximate under its stated prior; it is exact under D02's.
- **The global scale (derived by the user, 2026-10-26).** When β ≠ 0, GBMARS samples w by a modified half-normal step, because the drift puts w in the mean. Put the scale on u = √w instead and set zᵢ = u·vᵢ, so that εᵢ = βzᵢ + √(c·u·zᵢ)·ξᵢ with ξᵢ ~ N(0, 1). Write the mixing as GIG(p, ψ₀, χ₀), so that a stays the coefficients. Given z, the factors involving u are:
  - the likelihood, u^{−n/2}·exp(−Σ(rᵢ − βzᵢ)²/(2czᵢu));
  - the mixing density, u^{−np}·exp(−(ψ₀Σzᵢ/u + χ₀u·Σ1/zᵢ)/2);
  - a coefficient prior whose variance scales as u, u^{−M/2}·exp(−aᵀΣ₀⁻¹a/(2u)).

  Their product is a GIG kernel in u for any GIG mixing, so u has a GIG full conditional under any GIG-family prior on u: Gamma, inverse-gamma, or 1/u (Jeffreys on w). Two conditions:
  - *The prior goes on u, not on w.* An inverse-gamma prior on w adds exp(−b_w/u²), which falls outside the GIG family.
  - *The coefficient prior's variance scales as u, or not at all; never as w.* The g-prior's scaling, variance ∝ w = u², adds exp(−aᵀΣ₀⁻¹a/(2u²)) and also falls outside. Variance ∝ u is not scale-equivariant: rescaling y by k scales the prior variance by k, not k². Σ₀, or τ under the g-prior's Σ₀ = G⁻¹/τ, then carries y's units, so D03's Zellner–Siow default acts in raw units. That is OQ19's problem, settled by standardizing y. The scale-equivariant alternative keeps variance ∝ w. u's conditional is then u^{λ−1}·exp(−(A/u + Bu)/2 − Q/(2u²)), which is log-concave in log u, so one step of the `SliceSampler` block samples it.

  This matches Kozumi and Kobayashi (2011), whose σ multiplies the noise linearly and gets an inverse-gamma conditional. Theirs is the case with an unscaled coefficient prior and exponential mixing (χ₀ = 0). Checked numerically (scratch, 2026-10-26):
  - in the GIG cases, log p(u | ·) fits [1, log u, 1/u, u] to 10⁻¹⁶;
  - it misses by 2·10⁻² under an inverse-gamma prior on w, and by 3·10⁻³ with variance ∝ u².

  The z and v parameterizations are a centred against non-centred pair, so the interweaving schemes apply, with the same priors in both.
- **Defect fixtures from the text** (TeX source, checked 2026-10-25):
  - The residual is defined as r = ŷ − y, but the β and w conditionals hold only for r = y − ŷ. Implemented literally, the skew's sign flips. Fixture: recover a known β from skewed synthetic data.
  - Table 1's t(ν) row, GIG(−ν/2, 0, ν/2), gives under the paper's own GIG density a t_ν scaled by 1/√2. A standard t needs b = ν, which is what `StudentTWeights` implements, so the seam's t specialization must reproduce `StudentTWeights` exactly.
- **Sensitivity:** D08's draft, which declares the decomposed surface.
- **Sensitivity by quantile.** GBMARS fits one emulator per quantile level and tracks the Sobol indices across levels (its stochastic SIR example), showing which inputs matter in each tail. Two routes, the second a reading here, and a trade-off for `docs/Maths` rather than a default:
  - one fit per level, which can cross. The paper finds crossing seldom in practice, and less with replicates and with a flexible count prior;
  - quantiles read from one skewed fit, which cannot cross but rests on the generalized-hyperbolic family being right.
- **A generalized beta-prime prior on the per-observation variances** (the paper's supplement, SM2). It has mass near both ends of the shrinkage scale, so that a point is either fitted or ignored ("interpolate or ignore"), with the Horseshoe as a special case. Its full conditional is not tractable, so on BARS's side it would be a non-conjugate `IWeightResampler`, sampled by a Metropolis–Hastings or slice step. The paper reports "some success" with it in certain problems.
- **The weighted Gram came first.** Built 2026-11-02, with D03's unweighted prior and the marginal given σ² ([the upgrades](#upgrades-to-existing-code)). A design whose working response moves followed, the same day. What remains here is the GIG and GH pieces and the latent state.
- **The BARS side:** quantile and skewed BARS come from the same seam (the preparation brief's OQ20).

## Robust centring and bases for functional responses (OQ16)

**Deferred extension (2026-10-22): robust centring and bases for functional responses.** An extension beyond the parity contract, which uses R's mean centring and SVD. It is the functional counterpart of the preparation brief's OQ15 (the Normal prior's median centre), and is tracked there as OQ16.
- **A robust curve centre (D03).**
  - **Why.** In direct functional BASS the level to centre on is a curve. The pointwise median of the curves need not resemble any case, and it is not rotation-equivariant. The geometric (spatial) median in the space of curves is the standard robust functional location (Gervini 2008; Cardot, Cénac and Zitt 2013).
  - **What it keeps** of the scalar median's properties: it depends only on the data, so a Gibbs step over robust weights stays exact; it is translation-equivariant; and a minority of contaminated curves cannot drag it far (breakdown point 1/2).
  - **What exists.** `Maths.Geometry`'s `EuclideanMedian` computes it, by weighted Weiszfeld iteration from a warm start.
- **Robust bases for `bassPCA` (P02, D07).**
  - **Why.** `bassPCA` centres at the mean and takes an SVD, so one gross curve contaminates both the centre and the leading components.
  - **The functional-data answer** is spherical PCA (Locantore et al. 1999): centre at the spatial median, then take the eigenvectors of the spatial-sign covariance, Σ vᵢvᵢᵀ/‖vᵢ‖². Its eigenvectors are consistent for the covariance's under elliptical symmetry; its eigenvalues are not, so variance-based rank selection needs its own rule.
  - **`ComputeWithScatter`'s Weiszfeld scatter** weights each curve by 1/‖vᵢ‖ rather than 1/‖vᵢ‖², so a gross curve's influence still grows linearly with its distance. The spatial-sign form bounds it, so it is the safer basis.
  - **`MoMPCA`** (median-of-means aggregation of block PCAs on Rᵖ × Gr(r, p)) is a third route for large n.
- **The metric.** Distances between curves must be the L²(t) norm over the output locations, so each coordinate is scaled by the square root of its quadrature weight (√Δt on a non-uniform grid) before the median or the scatter is computed. `EuclideanMedian`'s weights are per case, not per coordinate.
- **The dependency to check.** `Maths.Regression` does not reference `Maths.Geometry` today. Geometry depends on LinAlg, Distance and Distributions, so the edge makes no cycle, but it must be checked against the repository's layering rules, or the routines placed where both can reach them.
- **What it does not do.** A robust centre and basis leave BASS's Normal likelihood non-robust. Downweighting whole contaminated curves needs case-level robust weights, a separate extension.

## Free-knot tensor-product surfaces (ABMRS)

**Parked extension (2026-10-20): free-knot tensor-product surfaces.**
- **What.** An expensive variant that carries BARS's continuous knot moves (birth, death and relocation per input dimension) to tensor-product spline surfaces, reserved for low-dimensional problems where its cost is tractable.
- **Pieces it would reuse,** once the BASS upgrades are done: the shared g-prior block, the `IBasis` carrier seam, the reversible-jump engine and the ensemble shell.
- **Cost.** For scattered data the Gram matrix is block-banded rather than a Kronecker product, so factor cost grows quickly with the per-dimension basis sizes.
- **Before designing it,** survey prior work on Bayesian free-knot surface splines. This is outside the completion contract.
- **Prior work found (2026-10-20).** He, Yang and Kang's EBARS (arXiv 2405.13353, 2024), published in *Biometrics* (2026) as ABMRS. The published details were verified against the article's open supplement on 2026-10-21 (preparation brief, step 3.8 checklist). It is the multivariate BARS this extension describes:
  - a discrete candidate grid on each axis, with tensor-product B-splines;
  - the same collapsed g-prior marginal (their Lemma 1), though the released code computes the BIC-type ratio;
  - Green's prior-aware schedule on the EBIC count law;
  - EBIC plug-ins for non-normal responses, described but not in the released package;
  - basis columns with at most `npart` = min(20, ⌊0.1·m⌋) supported points dropped from the fit, rather than the configuration refused.

  The rest of the survey stays open.
- **Design notes from reviewing it:**
  - **Support.** Step 3.7's column-mass support gate carries over unchanged, as Σᵢ Πₐ Bⱼₐ(xᵢₐ) over tensor columns. Near-empty corner cells are common there.
  - **Diagnostics.** Track knot counts per axis separately from the basis-column count. A birth on axis j adds Πₐ≠ⱼ(kₐ + pₐ + 1) columns at once, though the move itself is one-dimensional.
  - **Assembly.** Z's rows are row-wise Kronecker (Khatri–Rao) products with (p+1)^d nonzeros each. Accumulate G sparsely, at O(n·(p+1)^(2d)) rather than O(nν²); a move changes only its axis's one-dimensional basis matrix.
  - **Factor.** In lexicographic order the Gram is block-banded with half-bandwidth about p·ν₂ + p. Use a wide-band factor or step 2.2's `DenseCholesky`, chosen by size.
- **Its value here, and its sequencing (2026-10-25):**
  - *Value: smooth extrema.* BASS's hinges default to degree 1, which Francom et al. 2018 chose because "we often encounter instabilities when q > 1 because of the erratic tail behavior of higher order polynomials", so BASS's maxima sit at knots. ABMRS's cubic tensor B-splines would carry BARS's peak machinery (`SplineExtrema`, peak modes over full draws) to two or three inputs: a peak ridge over (T, θ) from families of SPC sweeps, for example. For peak and ridge readouts ABMRS is the better fit; for many inputs, BASS is (derivations review §7, item 9).
  - *Sequencing.* Building ABMRS first would not de-risk BASS. BASS's main risks are D01's term measure, the bordered factor maintenance, categorical terms and the functional contractions, and ABMRS touches none of them. It is scheduled on its own merits.
  - *A shared primitive.* Both designs are face-splitting (row-wise Khatri–Rao) products. Direct-functional BASS gets G = (AᵀA)∘(HᵀH) because cases × output locations form a Cartesian sample. ABMRS on gridded inputs gets a Kronecker Gram, and on scattered inputs accumulates G sparsely. One LinAlg contraction with a Cartesian shortcut serves both.

## Additional proposal kernels

These preserve the target when their forward and reverse probabilities are exact, and are validated before any becomes a default (the hub's modernization table).
- **Bernoulli-inclusion variable proposal** (KHAOS). Each input's inclusion is an independent coin with probabilities summing to an expected interaction order q₀; the reverse probability sums over q₀, and delayed rejection needs Poisson-binomial densities. It avoids Wallenius' distribution and any hard cap on the interaction order. KHAOS pairs it with a variable-swap mutation whose ratio is stated explicitly.
- **History-based variable proposal** (Nott and Kohn 2005, *Biometrika* 92). Variables weighted by inclusion frequencies accumulated over the run, aimed at many inert predictors. It is adaptive, so it follows D10's policy.
- **Residual-informed proposal mixtures** (the hub's modernization table), with exploratory support kept, as BASS's own births keep it.

## KHAOS's modified g-prior

An n-component g-prior that shrinks each term by its complexity: per-term factors from the interaction order and degree, and a global g₀² with an inverse-gamma prior (Rumsey, Francom et al., arXiv 2510.25036). KHAOS found its fixed ridge best on noise-free simulators and its g-prior better on most noisy ones. Its displayed conditional for g₀² has no data term, so the conditional is derived rather than copied. An option beside its penalty, built from the Jacobi generator's eigenvalues, is posed as [SG-4](../../../lean/Protolemmata/BASS/sensitivity-generator-lemmas.md).

## Basis-estimation uncertainty

D07 keeps the output basis fixed, so its estimation uncertainty is not integrated out. Fast Bayesian FPCA (Sartini, Zhou, Selvin, Zeger and Crainiceanu, arXiv 2412.11340) is one route to modelling it, through the outputs segment's basis contract.

## Moments for piecewise polynomials

A moment service written over piecewise polynomials rather than hinges only gives BAPS's tensor and anisotropic splines, and ABMRS's B-splines, variance decompositions from the same code as BASS's (derivations review §7, item 7).

## Active-subspace analysis

Posterior active-subspace analysis beside the Sobol stage (arXiv 2307.11241, *Discovering Active Subspaces for High-Dimensional Computer Models*, by overlapping authors; [PCA-extensions discussion](../discussion/astra-grassman-stiefel-pca-extensions.md)). For tensor-product representations, including MARS, the gradient integrals reduce to one-dimensional ones. It needs term derivatives, hinge reference cases, numerical integration checks and gradient-accuracy checks before its directions reduce another model, and it comes before jointly learning projection directions during fitting. [Gap G3](lineage-gaps-register.md#g3-sobol-indices-and-active-subspaces-have-different-standing) records why: unlike Sobol indices, C inherits no rate from a contraction result for f, E[C] is not the mean of the draws' subspaces, and under an empirical measure it reads the kink convention.

This entry is larger than its first form.
- **It is the moment service with derivative factors.** C = E[∇f ∇fᵀ] splits into three one-variable integrals per factor pair, ∫φ′ψ, ∫φψ and ∫φ′ψ′, which is what the moment service computes with `FactorDerivative` in place of the values.
- **Pairs of models.** Co-active subspaces generalize C to pairs of functions, C_kl = E[∇f_k ∇f_lᵀ], for adjacent simulator models (arXiv 2311.18146). That asks whether two fits, or two versions of one simulator, respond to the same perturbations in the same direction, which matching accuracy does not show.
- **Functional outputs, a reading here.** For a basis fit, f(x, t) = Σ φ_k(t) f_k(x), so the vector-valued active subspace of Zahm et al. (arXiv 1801.07922, unread) is Σ_kl ⟨φ_k, φ_l⟩_Q C_kl. So functional active subspaces are co-active blocks contracted with the output Gram, and F-X6's declared Q is the output norm their result depends on.
- **A consistency check per draw.** The diagonal of C is the derivative-based sensitivity measure, and a Poincaré inequality bounds each input's total Sobol variance by (b − a)²/π² times that entry, for a uniform input on [a, b]. That ties value moments to derivative moments. It is a consistency check, not ground truth.
- **Across engines.** Wycoff, Binois and Wild (arXiv 1907.11572, unread) give closed-form C for Gaussian processes, so one consumer covers BGP.
- **The sequence the lineage implies.** A posterior active subspace, then rotating the inputs and refitting, then BPPR (arXiv 2210.09181, unread), which learns directions during fitting. F-X3's test already implements a directional hinge, and Zhang et al. 2003 is the frequentist precursor.
- **Dependent inputs.** The first consumer that cannot assume a product measure, so the [lineage review's open question](lineage-review.md#separability-is-what-the-closed-forms-rest-on) is settled here, and [E25](#sensitivity-under-dependent-inputs) follows it.

## Cross-engine evidence

Evidence for BARS, BASS and BGP by path sampling (derivations review §7, item 8). It needs proper priors in every engine, stated in each engine's units. BASS has them by default since D03 was decided (2026-10-27: σ² scaled to the response, the intercept level-centred), as Normal BARS has. Under D03's evidence rule, evidence is reported from any run that can give one, empirical-Bayes in the level and scale, and the predictive checks already built for OQ17 (WAIC, PSIS-LOO) are the default comparison across engines.

**Held-out runs (added 2026-10-26).** Where runs are held out, proper scores compare any engines that produce predictive draws, with no evidence and no shared prior: the energy score over whole-output draws, and the interval score for marginal intervals (Gneiting and Raftery 2007, as used by Andros et al., arXiv 2406.18751). The energy score rewards the cross-output dependence that D07's joint draws keep, and it is where a missing truncation residual or an independent-noise model shows up. The scorers read draws only and know no model, like the OQ17 consumers.

**Why the engines compare at all (added 2026-10-28).** Farhi, Hed, Bon, Caticha, Mak and Domany (2013, *Physica A* 392, *Temperature Integration*) compute the free-energy difference between two systems by integrating each one's energy over temperature from a common high-temperature limit. That needs the two to share a phase space, so that Z agrees there. Evidence needs less: a normalized prior gives Z(0) = 1 in any dimension, so BARS, BASS and BGP compare without sharing a parameter space. BASS's evidence already has the paper's shape: D05's ladder doubles as the integration points, and its reference rung closes the integral at b = 0 exactly. SPC's notes adopt that closure (`issues/spc/opus-spin-glass-threads.md`, Thread 5).

**Bayes factors under the improper options (added 2026-10-28).** D03 reports a run's evidence only under proper priors, and that stands. A Bayes factor between two runs is a narrower quantity. Liang, Paulo, Molina, Clyde and Berger (2008, *JASA* 103) exclude improper priors on model-specific parameters, whose arbitrary constants stay in the marginal likelihoods (their §1), and put a flat prior on the intercept and 1/σ² on the variance once the columns are centred, so that both are common to every model (§2.1, eq. 3, written there on the precision). On common parameters the constants cancel in Bayes factors. Farhi et al.'s energy cutoff is the computational form of that argument: the same truncation in both systems, wide enough not to matter at the temperature of interest, and checked for insensitivity to where it sits.
- Here: truncate the improper priors identically in both runs, run each ladder to its reference rung, take the difference of the evidences, and widen the truncation until the difference stops moving.
- It holds only between models whose improper parameters mean the same thing: the same response, the same intercept convention.
- It also needs the same base model. With improper priors on the common parameters, Bayes factors depend on whether they are null-based or full-based (Liang et al. §2); BASS's centred marginal (`CollapsedGPrior.CentredMarginalAt`) is null-based.
- Cost: a truncated σ² or intercept changes the collapsed integrated form (incomplete-gamma and normal-CDF terms in place of the closed form).
- It is an option, not a gap: D03's proper defaults already give evidence.

**Complete values, declared (from [the notes'](../perplexity-bart-etc.md) capability list; built 2026-11-02, `4e3ea95`).** Most missing capabilities fail loudly, but a target with incomplete values still returned a number, and the condition was prose on `LadderEvidence`. A target now declares it: `ISamplingTarget.HasCompleteValues`, off by default. `FromLadders`, `FromSmc` and `Reweighting.LogEvidence` refuse without it. BASS declares it under proper σ² and intercept priors, and its three run-level copies of that condition read the one property. BARS reads `IObservationModel.HasCompleteMarginal`. A proper prior stays a separate condition, for an evidence. Not covered: `SmcResult.LogNormalizerRatio`, and the two estimators that take draws directly, which have no target to ask.

## The free-energy curve over the ladder

The ladder's evidence is one point, b = 1, of a curve: ln Z(b) at every rung, the free energy of each tempered posterior π(F)·p_b. Farhi et al. read ΔF at any temperature of interest, and three readouts sit on the curve here.
- **Tempered posteriors.** Each rung is Bissiri, Holmes and Walker's general-Bayes update π(θ)·exp{−b·loss} (2016, *JRSS B* 78, eq. 8, where b is their loss weight w) with the self-information loss −ℓ: the power likelihood of their §3.1. There 0 < b < 1 gives the data less influence and b = 0 keeps the prior. Their framework is motivated by misspecification, and their §3 offers ways to set b (annealing, unit information loss, operational characteristics) without prescribing one; ln Z(b) is not among them. What the ladder gives without new sampling is each rung's tempered posterior, so the sensitivity of any readout to b is a view over the ladder once its rungs' draws are kept. This is unrelated to GBMARS, whose "generalized" is the noise family (first entry).
- **WBIC** (Watanabe 2013, *JMLR* 14, eq. 4). The free energy from a single rung, E_b[−ℓ] at b = 1/ln N with ℓ the log-likelihood (−ℓ = nL_n in Watanabe's notation). It is valid for singular models, where BIC is not: its first two terms, −ℓ(θ₀) + λ ln N with θ₀ the parameter nearest the truth in Kullback–Leibler divergence, match the Bayes free energy's (Theorem 4 with b₀ = 1). With the hot end at 1/N that rung lies inside the ladder (at N = 200, 1/N = 0.005 and 1/ln N ≈ 0.19), and the ladder's own ln Z(1) checks it.
- **The real log canonical threshold λ.** At b = b₀/ln N, with b₀ a constant, E_b[−ℓ] = −ℓ(θ₀) + λ/b + U_N·√(λ/(2b)) + O_p(1), where U_N has mean zero and is asymptotically Gaussian (Watanabe 2013, Theorem 4). So the slope of E_b[−ℓ] against 1/b between two such rungs converges to λ (Corollary 3), and the ladder's other rungs in that range check it. λ is the Bayesian effective dimension of a singular model. Adaptive splines are singular: a term whose coefficient is zero leaves its knots and variables unidentified, so λ, not the basis count, measures their complexity. Which decision it would inform (calibrating the basis-count prior, comparing engines) is open.

**Needs.**
- D05's ladder with the reference rung, and its rungs' draws kept, beyond the scalar traces evidence reads.
- The integrated form's per-draw slope ∂_b ln p_b at a rung. BASS's target supplies it (`BassEvaluation.LogLikelihoodSlope`), but its ladder traces slopes on the linear form only, since on the integrated form the slope at b = 0 carries E[1/τ], infinite under Zellner–Siow, though finite at 1/ln N. Tracing it from a b_lo above 0, or with τ held, is an [upgrade](#upgrades-to-existing-code).
- A rung pinned at 1/ln N, or a rule for reading between rungs, since an adaptive ladder places its interior rungs by swap rejections.

## Elastic emulation of misaligned functional outputs

From *Elastic Bayesian Model Calibration* (Francom, Tucker, Huerta, Shuler and Ries, arXiv 2305.08834).
- **The problem.** When functional outputs differ in phase, so that their features move along t, PCA on the raw curves mixes phase with amplitude. It needs more components, and its averages blur features that each run has sharply.
- **Their method.** They align the runs in the square-root-velocity framework, which splits each curve into an aligned curve and a warping function. They represent each warping function by its shooting vector: the inverse exponential map of its square-root-slope function on the unit Hilbert sphere. They then emulate the aligned curves and the shooting vectors with separate BMARS emulators, and calibrate modularly on both.
- **What it needs here.**
  - The outputs segment's representation with a declared, nonlinear reconstruction map (composition with the warping).
  - Two ordinary BASS basis fits, one for the aligned curves and one for the shooting vectors.
  - D09's cut likelihood over both parts.
  - The alignment itself, which is new code; the sphere and its exponential map relate to `Maths.Geometry`.
- **A simpler route to compare first: landmark warping** (Francom, Sansó and Kupresanin, *SIAM/ASA JUQ* 10, 2022; unread). The predecessor to elastic warping, for outputs with identifiable features such as arrival times, spectral lines or peaks, where a landmark can stand in for a square-root-velocity alignment.
- **Why it matters here.** A susceptibility peak whose temperature moves with a control parameter is exactly this situation. It is a second route to families of SPC sweeps, beside ABMRS's peak readouts (next entry).

## Residuals correlated along the output index

From Andros, Guhaniyogi, Francom and Pasqualini (arXiv 2406.18751, 2602.09351), whose models give each run a Gaussian-process residual over the output domain.
- **The problem.** Direct BASS treats noise as independent across output locations. When the misfit is smooth along t, a dense grid then counts n·m independent observations, so intervals are too narrow and births absorb smooth residual. The fitting segment's residual-dependence check detects this; this entry models it.
- **The model.** Residual covariance σ²(I ⊗ R_θ): independent across cases, correlated along the index. Because R acts only on the output side, the direct contractions keep their separable form with Q = R⁻¹ (the fitting segment's output-side inner product). The collapsed marginal gains a log|R| term per case, and the g-prior is taken in the same inner product.
- **Cost.** For a Markov residual on a one-dimensional index, Q is tridiagonal and every output-side product is linear in m. Sampling θ recomputes the output-side Gram for the current terms, O(M² m) with a banded Q.
- **Not covered.** The PCA mode's residual is the empirical truncation residual, whole-vector, so it already keeps within-output dependence (D07); this entry concerns direct fits.

## Semi-modular calibration

D09 takes the cut posterior: the observations update θ given each emulator draw, and never feed back into the emulator. Bissiri, Holmes and Walker's partial-information update (2016, *JRSS B* 78, §4.2) is its coherence argument. Updating θ on only the part of the information that bears on it is a valid belief update, though it is no joint model's posterior.
- **The extension.** Semi-modular inference (Carmona and Nicholls, "Semi-Modular Inference: enhanced learning in multi-modular models by tempering the influence of components", *AISTATS* 2020, PMLR 108:4226–4235; arXiv 2003.06804; not yet read here) puts a power η on the feedback term. Its follow-up, "Scalable Semi-Modular Inference with Variational Meta-Posteriors" (arXiv 2204.00296), fits the whole η family at once. η = 0 is the cut and η = 1 the joint posterior, with η chosen by a predictive criterion.
- **Why it is cheap here.** Stage 8's cut SMC already bridges successive conditionals through induced densities, and a path in η is one more bridge sequence on the same machinery.
- **Needs.** D09's emulator interface (C-X2) and stage 8's cut SMC (C5). The joint end, η = 1, also needs the emulator's likelihood under the observations, which the cut never forms.

## Calibrated intervals for loss-based posteriors

BARS's `GibbsLossModel` is a Gibbs posterior on Tukey's bisquare loss, and its docstring says its credible intervals have no calibrated interpretation.
- *Its scale is already set.* Its normalization, ρ(u) ≈ u²/2 on residuals scaled by a robust σ̂, is effectively Bissiri et al.'s unit-information weight (§3.2), which for squared loss gives w = 1/(2σ²), the Gaussian case.
- *Its coverage is not.* A unit-information weight fixes the loss's scale, not the intervals' coverage.
- **The extension.** Set the loss weight w so that credible intervals reach frequentist coverage: Bissiri et al.'s operational calibration (§3.4). The algorithms that do it are its follow-ups, their details checked on 2026-10-29 and not yet read here:
  - Holmes and Walker, "Assigning a value to a power likelihood in a general Bayesian model", *Biometrika* 104(2):497–503 (2017), doi 10.1093/biomet/asx010, arXiv 1701.08515: the power set by matching Fisher information;
  - Lyddon, Holmes and Walker, "General Bayesian updating and the loss-likelihood bootstrap", *Biometrika* 106(2):465–478 (2019), doi 10.1093/biomet/asz006, arXiv 1709.07616: a bootstrap sampler for loss minimizers, which also calibrates the weight;
  - Syring and Martin, "Calibrating general posterior credible regions", *Biometrika* 106(2):479–486 (2019), doi 10.1093/biomet/asy054, arXiv 1509.00922: general posterior calibration (GPC), which tunes the weight by bootstrap until the regions reach nominal coverage;
  - Wu and Martin, "A comparison of learning rate selection methods in generalized Bayesian inference", arXiv 2012.11349: a survey of these and others.
- **Where it applies.** BARS's Gibbs-loss route first, then any engine that weights a likelihood for robustness, such as a power posterior under misspecification. BASS's tempered rungs are outside it: they are a path to b = 1, not a posterior to report.
- **When.** Only if the Gibbs-loss route is used for intervals, rather than the scale-mixture resamplers BARS prefers.

## General-Bayes readings for the port

Bissiri, Holmes and Walker, "A general framework for updating belief distributions", *JRSS B* 78(5):1103–1130 (2016), arXiv 1306.6430, were read against the port on 2026-10-29. Most of the paper is already absorbed:
- each tempered rung as its power-likelihood update (§3.1, eq. 8), in the free-energy entry above;
- BARS's `GibbsLossModel`, a loss-based update on a unit-information scale (§3.2).

Two readings frame decisions without changing them.
- **D09's cut posterior** is §4.2's partial-information update; see semi-modular calibration, above.
- **The deterministic-simulator profile (F-R1).**
  - *The weight is a precision.* §3.3's hierarchical loss with l(w) = log w gives w^ξ·exp{−w·l}·π: a Gamma prior on a precision.
  - *So σ² is a loss weight.* For a noise-free simulator, BASS's sampled σ² is a learned weight on the squared-error loss, not a noise variance.
  - *And the update is valid.* §4.1 makes updating on non-stochastic data a valid update.
  - *What it means for F-R1:* it supports keeping σ² sampled in the profile, and it reads the exact-fit impropriety under Jeffreys as w → ∞.

Not taken:
- §4.1's loss-based expert statements, with no current need;
- §3.5's conjugate loss prior, since BASS's priors are already conjugate;
- the supplement's joint quantiles, which would become a quantile-emulation entry only if wanted.

§5.2's model-free clustering belongs to SPC.

## Conformal Bayes, the consumer

The seam is built and the consumer is not. OQ17 (2026-10-23, built 2026-10-24; BASS's terms in S9) gives:
- every draw's state, retained;
- the n × S matrix of each unit's log density under each draw;
- a unit map, under which a whole curve can be one unit;
- a stratum map for the grouped (partially exchangeable, or Mondrian) form.

The conformal step itself was deferred to the SPC × BARS work (`docs/Maths/BayesianRegression.md`, the pointwise section).
- **The step** (Fong and Holmes, *Conformal Bayesian Computation*, NeurIPS 2021; arXiv 2106.06137, a local copy in `codex-scientiae/supellex-private/gauntlet/bayesian/2106.06137v2`). Add-one-in importance weights, w_s ∝ p(y | x, θ_s), turn the retained draws into conformal prediction sets for a candidate y, with no refits. Coverage holds in finite samples under exchangeability, up to Monte Carlo error, whether or not the noise model is right.
- **A full predictive distribution at each x.** Sweeping the candidate y gives the conformal p-value as a function of y: a calibrated predictive CDF, the conformal predictive distribution of Vovk, Shen, Manokhin and Xie (*Machine Learning* 108, 445–474, 2019; arXiv 1710.08894). It comes from the same draws and weights, at the cost of a grid in y, and gives PIT-style calibration and a calibrated CRPS.
- **The grouped form** is the noise-model calibration check: replicates at one design point are exchangeable, so their coverage tests the posterior's predictive width there.
- **What it cannot reach** (recorded with OQ17). On SPC × BARS's fixed temperature grid, only replicates at one temperature are exchangeable. So coverage holds for a new replicate at a measured temperature, in steps of 1/(k + 1), and says nothing about unmeasured temperatures, the curve or the readout's features. The next two entries address that.

## Conformal bands for whole curves

Simultaneous prediction bands for a new unit's whole curve: a new simulator run's output for BASS, a new replicate sweep for SPC × BARS.
- **The sources.**
  - Lei, Rinaldo and Wasserman, *A conformal prediction approach to explore functional data* (*Annals of Mathematics and Artificial Intelligence*, 2015; arXiv 1302.6452): bands from a finite-dimensional projection and the inductive conformal predictor.
  - Diquigiovanni, Fontana and Vantini, *Conformal prediction bands for multivariate functional data* (arXiv 2106.01792): exact bands around any estimator, in closed form, whose width can follow the data's local behaviour. Its univariate companion is *The Importance of Being a Band* (*Statistica Sinica* 35, 853–871, 2025).
- **Here.** The nonconformity score is the curve's largest scaled residual, the maximum over locations of |y − f̂| over the posterior's pointwise predictive sd, so the band widens where the posterior is uncertain. The unit map already takes a whole curve as one unit, and conformal Bayes's add-one-in step applies to it.
- **Why it matters for BASS.** Andros et al. 2024 found R BASS's 95% intervals covering 0.59 on held-out runs. D07's leave-one-out truncation residuals improve the intervals; conformal bands would guarantee their coverage for a new run, under exchangeable runs.
- **Needs** exchangeable whole-curve units: the runs of a simulator design, or replicate SPC sweeps over one temperature grid.

## Posteriors for readouts by resampling units

The posterior of a population feature, such as the peak T* of SPC × BARS's curve, robust to the noise model being wrong. Conformal methods cover observables only, so this needs another route. Three relatives share one idea, that the uncertainty comes from resampling exchangeable units:
- **BayesBag** (Huggins and Miller, *Reproducible model selection using bagged posteriors*, *Bayesian Analysis* 18(1), 79–104, 2023; arXiv 2007.14845): average the posteriors over bootstrap resamples, typically 50–100 refits. It is more stable than the standard posterior under misspecification, and slightly more conservative under a correct model.
- **Martingale posteriors** (Fong, Holmes and Walker, *JRSS B* 85(5), 1357–1391, 2023): uncertainty as missing future data.
  - Impute the rest of the population by predictive resampling, and compute the feature on the completed data.
  - With the Bayesian predictive this recovers the usual posterior (Doob's theorem).
  - With a predictive that resamples units (a Bayesian bootstrap), it gives a posterior that does not rest on the likelihood being right, which is BayesBag's construction in another reading.
- **Syring and Martin's calibration** (E15) is the coverage-calibrated variant.
- **The synthesis with BARS** (SPC × BARS): refit BARS on resamples of replicate sweeps, and pool the peak readouts its draws already carry. Together with the two entries above, that gives one layer over BARS's retained draws and the unit map, at three levels:
  - a new replicate at a measured temperature (conformal Bayes);
  - a new sweep's whole curve (the bands);
  - the readouts (this entry).
- **The constraint.** Both curve-level levels need replicate sweeps as the exchangeable unit. On a single sweep, nothing distribution-free exists for curve features, so replicate sweeps are worth budgeting for in SPC runs.

## Polynomial-chaos terms

KHAOS's basis (Rumsey, Francom et al., arXiv 2510.25036), for smooth responses, where polynomial chaos is the standard surrogate. A term is a product of one-variable shifted-Legendre polynomials, with a structure defined by a degree split across variables rather than by a knot. E4 and E5 hold its proposals and prior; this entry is the term type.
- **It fits F-X3 as a case-side kernel.** Its factors are one-piece polynomials, so `PiecewisePolynomial` integrates them with no new moment code.
- **A cross-check it brings.** Under a uniform measure, Sobol indices reduce to sums of squared coefficients, an independent check on `BassMoments`.
- **Past the paper.** KHAOS assumes uniform, independent inputs. The moment service extends it to `MixtureMeasure` and `EmpiricalMeasure`.
- **New.**
  - The degree-split structure and its birth proposal.
  - Poisson-binomial densities for delayed rejection, with E4's Bernoulli-inclusion proposal.
  - Adaptive mutation probabilities, under D10.
  - One step of `SliceSampler` on ln g₀², in place of the paper's Laplace-approximation independence proposal, with its conditional derived (E5).
- **Defect fixtures from the TeX, to confirm against the `khaos` source.** The ratios are in [SC-1](../../../lean/Protolemmata/BASS/structure-chain-lemmas.md). The rest are probable typos: P₁(x) = 1 where x is meant; the order sum runs over M where p is meant; the degree is drawn from {q, …, p}, probably meaning d_max; the partition probability is binom(d−1, d−1)⁻¹ where binom(d−1, q−1)⁻¹ is meant.

## Data-augmented observation models

F-X1 hands the core a working response and weights, so any observation model that turns Gaussian given a latent variable is a member. The [lineage review](lineage-review.md#f-x1-is-a-data-augmentation-seam) tabulates them: probit, ordinal, logit through Pólya–Gamma, a scale per declared group, a mixture of normals.
- **The weighted Gram, centring and marginal.** Built 2026-11-02 ([the upgrades](#upgrades-to-existing-code)). A latent's precision weights read `WeightMeaning.Precision` and D03's `WeightedPriorGram.Unweighted`, which keeps its conditional in closed form.
- **Probit needs more than unit weights.** An earlier reading called it free. It also needs:
  - σ² held at 1, for identification: `Sigma2Prior.Held(1)` (built 2026-11-02), with the marginal given σ² from the same S and determinant term. D04's generalized-hyperbolic mode reads that form too when σ² is not integrated, through `BassTarget.Given`;
  - priors that stop reading the response. D03's σ² scale and level-centred intercept read it through `BassResolved.Resolve`, and here the working response moves every sweep, so this entry names its route among the ones built ([the options](moving-response-priors-options.md));
  - the latent state and its block in F-X2's slot, as for GBMARS, with the cross-products recomputed each time the latents move;
  - a truncated-normal draw, tail-stable and qualified apart from density evaluation: `Normal.MeanSd(μᵢ, 1).SampleOn(rng, 0, ∞)`, built 2026-11-02 ([the upgrade](#upgrades-to-existing-code)).
- **Putting the offset back.** Prediction and the pointwise densities read y today, not y − offset (the F-X1 report).
- **Groups.** A heavy-tailed scale per declared group, rather than per observation, is GBMARS's latent indexed differently, and it gives units sharing a group joint extremes (E22).
- **Mixtures.** A Dirichlet-process mixture of normals hands the core (y − μᵢ, 1/σᵢ²). Its mixing law is shared across x, so it gives an irregular residual, not heteroscedasticity, and its components are a density estimate rather than discovered populations.
- **Discrete-time survival** (Sparapani et al. 2016, from memory, via the [notes](../perplexity-literature-deep-dive.md)) is probit on a data set expanded per time point, so it arrives with this entry plus its own expansion and readouts. What a run's ending means, censoring or a competing outcome, is a modelling choice to declare.
- **Not covered:** a scale that varies with x, which needs a second modelled surface; and likelihoods with no Gaussian augmentation, for which Linero's generalized BART (arXiv 2202.09924, unread, via the [notes](../perplexity-bayesian-trees.md)) uses nonconjugate reversible-jump updates.

## A prior-precision builder, with ridge blocks

KHAOS's modified g-prior needs the collapsed core to accept a prior precision other than ΨᵀΨ/g. [SC-3](../../../lean/Protolemmata/BASS/structure-chain-lemmas.md) shows G∘ΨᵀΨ = ΨᵀΨ + S₀⁻¹, so it stays conjugate. Declaring that as a builder, rather than adding one more case to the coefficient-prior set, gives a second kind of prior: a ridge block.
- **Random-effect intercepts** on any grouping variable with many levels are ridge blocks, so E5's work pays off for every grouped-data problem.
- **Scale-sensitive priors.** Ridge priors, unlike Zellner's, change when columns are rescaled, so they need the per-kernel standardization E28 raises.
- **Identifiability.** A random-effect block overlaps the level-centred intercept unless it is centred or sum-constrained, so its convention is declared with it.

## Correlated outputs and coupling within a run

For several outputs from one run, or several outcomes for one case. Means do not reach how outcomes combine within a run: P(all succeed | x), conditional failure rates and the distribution of the number of successes all depend on it, however good the joint posterior of the means.
- **Two reasons outcomes co-occur.** A hard configuration raises both, which the mean model captures. Within-run coupling is one run going bad across the board, which only R captures. Telling them apart needs repeated runs at one configuration.
- **The model.** Latent z ~ N(μ(x), R) per run, with F-X6's Q = R⁻¹ for the mean. R is the correlation matrix, identified as one (Chib and Greenberg 1998). Binary and ordinal outcomes come through E20.
- **A third inner product.** Beside the identity and E12's tridiagonal, Q = (ΛΛᵀ + D)⁻¹ by Woodbury, with Q-centred forms (documented on `KhatriRao`; the design refuses Q ≠ I today).
- **Sampling R, in `Maths.Samplers`, apart from any model.** Parameter expansion for a free R: one conjugate inverse-Wishart draw, then a rescaling to a correlation matrix (Liu and Wu 1999; Lawrence et al. 2008). Factor loadings under parameter expansion for a low-rank R (Ghosh and Dunson 2009). All from memory.
- **Where a Gaussian latent stops.**
  - *Joint extremes.* A Gaussian copula has no tail dependence, so it underweights runs where everything fails together. A heavy-tailed scale shared across a run gives positive tail dependence (E20's groups), as does a heavy-tailed factor.
  - *Weak identification.* When failures are rare, the information about R is in the co-failures, which are rarer still. A low-rank R with shrinkage holds up where a free K × K one does not.
- **Alternatives, with different readouts.**
  - *A chain factorization,* P(Y₁ | x)·P(Y₂ | Y₁, x)·…, reusing the binary machinery with earlier outcomes as categorical inputs. Directed and conditional, where the joint latent is symmetric; it needs an order.
  - *A rank likelihood* (Hoff 2007) for mixed outcome types, estimating the copula without specifying the marginals.
- **A proposed order:** a factor latent with one run-level factor, then the shared heavy-tailed scale, then ordinal cutpoints.
- **Parked:** R(x), covariance regression (Hoff and Niu 2012), until a fixed R proves too stiff.

## Multi-fidelity emulation through the offset

Kennedy and O'Hagan's autoregressive model, y_hi = ρ·f_lo + δ, is close to a configuration here. F-X1's offset, built for calibration discrepancy, carries it: fit δ against y_hi − ρ·f_lo, with stage 8's cut deciding whether the high-fidelity data feed back into f_lo.
- **Fixed against drawn ρ.** A fixed ρ, or one draw of f_lo at a time, is a plain offset. Drawing ρ moves the working response every sweep, so it meets E20's question of priors that read the response.
- **The same structure, read across.** Any cheap and costly pair of a code, or a surrogate against the real system. δ is estimable only where both are observed.
- **Prior work.** Tucker, Collins, Huerta and Brown (arXiv 2609.27174) use BASS emulators for multi-fidelity mapping and fusion, with elastic alignment. Read it to the end, then decide what this entry adds.

## Guarded estimators

Estimators that stay valid whatever the emulator does: the design supplies validity, the emulator supplies efficiency. They suit any problem where generating a configuration is cheap and running it is not. The mathematics is [GE-1 to GE-3](../../../lean/Protolemmata/BASS/guarded-estimator-lemmas.md), proved in `GuardedEstimator`. This entry is the layer beside the engine, taking predictions from any emulator.
- **What it reads:** a per-draw integral under a declared measure (the moment service; the Monte Carlo fallback elsewhere), predictions at points, and a way to fit on folds. Under draws from the generator with the integral, the cross-fit bias is exactly zero (GE-2), which is a concrete return on the lineage's closed forms.
- **The tuned weight, not a fixed difference.** With a poor emulator, Var(y − f̂) can exceed Var(y). Survey sampling's regression estimator and PPI++ tune it, and the tuned weight is asymptotically never worse than either end.
- **Levels.** Several costs stack as a telescoping sum: multilevel and multifidelity Monte Carlo, and approximate control variates with a tuned weight per level (Gorodetsky et al. 2020, from memory).
- **Adaptive runs.** Unbiasedness needs a random sample, which targeted designs (E30) break. Two options, kept apart until one is chosen, with the trade-off for `docs/Maths`:
  - *An audit stream.* A fixed share of runs is drawn at random for the design-based estimates, and the emulator fits on all. Simple, but those runs do less for the design.
  - *Known probabilities.* Each configuration is drawn with a probability the emulator sets, above a floor, and its residual is weighted by the inverse. Every run serves both purposes; variance grows as the floor drops.
- **Limits.** The guarantees cover linear quantities, reweighted means and, through PPI, regression coefficients and quantiles. Sobol shares, boundaries and subspaces still depend on the emulator, and are checked by held-out prediction. After reweighting to another measure, report the Kish effective sample size of the correction term, which says whether the design covered that measure.

## Sensitivity under dependent inputs

Observational covariates, constrained parameter spaces and any calibrated posterior all break the product measure, and D09's θ posterior is never one. Francom and Nachtsheim (arXiv 2506.11471) sort methods by "enough runs for an emulator?" and "independent inputs?". In those terms the engine covers one quadrant.
- **The measure.** [G4](lineage-gaps-register.md#g4-one-moment-service-across-dependent-inputs)'s joint mixture of products keeps every pair integral closed, from one component (`ProductMeasure`) to N with delta marginals (a joint sample). The G4 limits apply: moments and C are exact, and Hoeffding is lost.
- **The readouts.**
  - *Accumulated local effects* (Apley and Zhu 2020, from memory) need only derivatives and the empirical joint, and O-X4 already reaches the derivatives, so they are the first consumer.
  - *Shapley effects* cost more, dropping under G4 from a double loop to a single one.
  - *Monte Carlo over θ draws with O-X4's analytic gradients,* the fallback for any measure and any kernel.
- **The consumer's shape** is the [review's open question](lineage-review.md#separability-is-what-the-closed-forms-rest-on), settled at E8.

## Stochastic emulation through latent inputs

Simulators whose output law at one input can be multimodal: agent-based, epidemiological and Monte Carlo codes. GBMARS's latent scales handle only unimodal laws. Zhu and Sudret (*IJUQ* 2023, arXiv 2202.03344, unread) add a latent input and noise to a polynomial-chaos model, which captures multimodal output laws without replicate runs.
- **Here.** F-X2's latent slot holds a per-observation latent, and a kernel reads it as one more input: E19's polynomials, or a hinge. It is a different model from GBMARS's, with the latent in the basis rather than the noise.
- **What is new:** the latent's update enters the design, so it moves columns rather than the working response, which F-X1 does not cover. The cost and the mixing of that update are the open questions.
- **Before starting:** read the paper for how it fits and selects terms, since it is not known here to be an MCMC fit.

## External benchmarking

*All Emulators are Wrong…* (Rumsey, Gibson, Francom and Morris, arXiv 2512.09060, unread) compares about 30 emulators on 60 test functions and 40 real datasets, scored by CRPS, through the `duqling` R package.
- **Why.** The port's oracle tests check parity with R BASS. This checks the engines against the field across many domains, with ground truth from outside any estimator.
- **Needs:** predictive draws from each engine and E9's proper scorers, which read draws only. The test functions and datasets, and whether they can be reached without running R in the loop, are to be read first.

## Bases that mix term types

An option, not a recommendation. F-X3's contract already allows a term to mix kernels across variables (hinge × categorical exists), so a Legendre part on x₂ times a hinge part on x₁ fits the rule, and collapsing makes kernel-changing moves discrete Metropolis–Hastings steps with no Jacobian. No mixed hinge-and-polynomial basis is known in the lineage, which is an impression rather than a search result.
- **The structure prior.** Each kernel has its own term measure, so each variable's part gets a kernel label, and that label's prior has to be chosen.
- **The coefficient prior.** Zellner's g-prior is unchanged by rescaling columns. Ridge priors are not, so they need per-kernel standardization (E21). KHAOS's complexity factor is defined on a multi-index, so each kernel would report its own complexity.
- **Mixing and identifiability.** One smooth effect can be one polynomial term or several hinges, and the path between them passes through poor states, which is the multimodality the tempering machinery is for. f's functionals do not care which representation the sampler is in. Structural summaries do: inclusion counts, term counts and the adaptive proposals that feed on them (KHAOS's variable weights, Nott and Kohn) say less the more freely the basis mixes.
- **The alternative and a middle.** One fit per basis, compared by evidence or held-out scores (E9). Or one label per variable rather than per part: fewer labels and a simpler prior, which still gives "kinked in x₁, smooth in x₂".

## Categorical levels: learned positions and a similarity readout

Three ways of sharing strength across a categorical input's levels. Latent-variable GPs (Zhang, Tao, Chen and Apley 2020) place each level at a point in a low-dimensional space. Multidimensional item response theory gives each a latent vector. BASS's categorical factors group levels into subsets, term by term.
- **A readout, cheap and untested.** How often two levels share a subset across the posterior, weighted by each term's variance, gives a similarity matrix from any fit with categorical factors, with no model change.
- **A term type, the lead.** A level at a learned continuous position, with a hinge on that position: LVGP's idea as a BASS factor. It stays separable, so under a categorical measure the moments are a finite sum. It is not conjugate in the positions, and it is identified only up to rotation and translation, so its cost is statistical as well as in code.

## Choosing the next runs

Step five of the simulator workflow, for which no seam exists: where to run next, given a fitted emulator.
- **Replicate or explore** (Binois, Huang, Gramacy and Ludkovski, *Technometrics* 2019), with noise that varies across the space, so it wants a noise model that varies with x.
- **Contours** (Ranjan, Bingham and Michailidis 2008): place runs near a level set of the surface.
- **From psychophysics, which has done this for decades.** A 50% threshold is a psychometric threshold: QUEST (Watson and Pelli 1983), the Psi method (Kontsevich and Tyler 1999), QUEST+ (Watson 2017) for several dimensions, and qCSF (Lesmes et al. 2010) for a threshold curve along a second dimension. Its psychometric function also has a lapse rate, success topping out below 1, kept apart from the slope. Contour designs from computer experiments usually have no such term.
- **With E24.** Targeted runs break a design-based guarantee unless their probabilities are known, so the two are designed together.
- **Unread, closest on the lineage's side:** Hutchings et al. (arXiv 2405.16298), active learning for functional output.

All of these are from memory or from the [evals note](../opus-bayesian-evals.md), to be read before this entry starts.

## Tree ensembles, with BCF's multiplier (E31)

BART (Chipman, George and McCulloch 2010) and its relatives, as a regression family of their own beside BARS and BASS. The source is [the notes](../perplexity-bart-etc.md) and [Gemini's digest of them](../discussion/gemini-post-bass-feedback.md), whose code claims were checked on 2026-11-02. The notes read the repo at `440d8a7`, before the extension seams. The reading list (SBART, DART, GP-BART, survival and NFT BART, XBART) stays in them, and the [lineage review](lineage-review.md#against-the-notes-in-issuesregression) has the engine boundary.
- **A leaf contract beside `IObservationModel`, from BART's backfitting.** The update integrates one tree's leaves given partial residuals and a shared σ²: `LogIntegratedLeafLikelihood(stats, noiseVariance, leafPrior)` and `SampleLeaf`, over a leaf's count, Σr and Σr². Indicator columns in a design would lose those O(n) statistics.
- **The component multiplier, from BCF.** A tree's contribution is aᵢ·g(xᵢ), so the count and Σr become Σa² and Σa·r. aᵢ = zᵢ is only one coding, under which controls say nothing about the treatment forest. Hahn, Murray and Carvalho's (2020) invariant coding, b_{zᵢ}, lets them inform it (from memory, unread here). So the contract takes any known aᵢ.
- **The scan, from forests' mutable state.** Grow and prune fit `IRjMove` as BASS's collapsed moves do, with no coefficient map. A forest holds cached predictions and leaf assignments, so `BlockScan`'s `copy` has to be real. BASS and BARS pass an identity, since each saves its state itself (`BassSweep.cs:512`, `BarsEnsemble.cs:452`). A resume bit-identity test catches a wrong one.
- **Block order, from the digest's fixed pipeline.** It puts shrinkage last, while D04 puts τ second because τ's block integrates the coefficients out. A partially collapsed sampler can lose its target when reordered (van Dyk and Park 2008, from memory). So each family declares its order and checks it by enumeration, or by Geweke's joint test (2004, from memory) where the space is too large to enumerate. The repo has no Geweke test.
- **Tempering a forest (a reading here, untested).** The notes warn that one tree's integrated likelihood under b is not the forest's. That is true, and no forest-integrated form is needed. Temper the joint target, prior × p(y | trees, leaves, σ²)^b. Integrating one tree's leaves out, given the rest, is a valid conditional of that target. Swaps read b·ℓ₁ on the full state once the leaves are drawn, `ITargetEvaluation`'s separable case. The joint's margin on the structures is the integrated form, not the linear one.
- **Readouts, from SoftBART and ridgeBART.** A standard BART draw is piecewise constant, so derivatives and peak counts need a smooth variant. No generic posterior-function contract exists for one to fall through today: `SplineExtrema` is BARS's own. The digest's "smooth" overstates BASS as well, whose default `Degree = 1` has kinks at its knots. So a shared posterior-output contract declares smoothness, which BASS needs as much as trees do.

## Additive models with a GP component (E32)

y = g(x) + f(x) + ε, with g a GP and f a forest or BASS, and treed GPs. The source is [the notes](../perplexity-bart-etc.md), and the question they left open was checked on 2026-11-02.
- **`GpRegression` today.** It fixes `y`, its shift and the squared distances at construction. A known σ² sets s² = 1, and the unknown-scale mode learns both.
- **What it needs.** `y` per fit, which the distances do not read; σ² given from outside with s² free; and a declared mean, in place of the training mean it subtracts today. `BgpSampler` follows.
- **The other component.** BASS needs the marginal given σ² ([above](#upgrades-to-existing-code)). The lineage review's warning stands: such a sum is a new model, with centring to declare and component recovery to test.

## Regression-family integration (E33)

[The notes](../perplexity-bart-etc.md) put this after the spline families, with trees last as a test that the engine is basis-agnostic. They also pose a scope question, a library of qualified families or a general model-definition system, and favour the library.
- **A capability inventory.** For BARS, BAPS, the hierarchical P-spline, BGP and BASS: what each can evaluate, integrate, draw conditionally, predict, resume, temper and report as evidence. Each gap is marked as an upgrade to make or as true by the nature of the model. The [completion matrix](completion-audit-matrix.md) covers BASS against R, and the lineage review's axes seed the rest. Written 2026-11-02: [the capability inventory](../capability-inventory.md), with each gap marked an upgrade or true of the model, and the rows it added to the table above.
- **Declared capabilities.** Partly there: `CanTemperLikelihood` and `CanTemperToPrior` on BARS's model, and `BassConfig.Validate`'s refusals. `ITargetEvaluation` refuses a b or a slope by exception at first evaluation, so after a run starts. The inventory decides what else earns a declaration. Evidence's is [above](#cross-engine-evidence).
- **Four verdicts.** Target correctness, computational adequacy, predictive adequacy and evidence qualification, each certified apart. F-R1 is the example: the oracle and enumeration checks pass while no arm clears the R̂ gate. `ConvergenceVerdict` reads computational adequacy alone. A new family's qualification report reads by these four.
- **No composition layer.** `BlockScan` and `IBlock` already compose blocks, for collapsed BASS and explicitly blocked BAPS alike. What the notes' own sequence asks for is contracts on block invariance, tempering semantics, normalization and restart.

## Direct-functional BASS for families of SPC sweeps

Fitting f(x, T) across families of SPC temperature sweeps (derivations review §7, item 9); BARS remains the engine for a single sweep. For peak and ridge readouts over (T, θ), ABMRS is the better fit (its entry above), and elastic emulation handles peaks that move with θ. An application, outside the completion contract.

Andros et al.'s functional Gaussian process (arXiv 2602.09351) addresses the same problem with the factorization the other way round: a fixed spline basis in the inputs, whose coefficients are GPs over the output domain. BASS-PCA fixes the output side from the data and adapts the input side; direct BASS adapts both, with sparse hinge products. With one or two control parameters a fixed input basis is cheap, but the model's marginal likelihood is a GP over all n·m outputs, so it is a comparison point here, not a route the port prepares for.
