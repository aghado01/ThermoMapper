# BASS extensions brief

Extensions beyond the port's completion contract, collected from the [full port brief](../bass-full-port-brief.md) and the reviews behind it. Each entry says what it adds, its source, what it needs from the port (each segment's "Room for extensions" section), and its status. Nothing here is scheduled: an entry starts on the user's go-ahead, with its own specification and acceptance criteria. Its checklist, where the reports go, is [extensions-checklist.md](extensions-checklist.md).

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
| Direct-functional BASS for families of SPC sweeps | Derivations review §7, item 9 | The completed port | Application |

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

An n-component g-prior that shrinks each term by its complexity: per-term factors from the interaction order and degree, and a global g₀² with an inverse-gamma prior (Rumsey, Francom et al., arXiv 2510.25036). KHAOS found its fixed ridge best on noise-free simulators and its g-prior better on most noisy ones. Its displayed conditional for g₀² has no data term, so the conditional is derived rather than copied.

## Basis-estimation uncertainty

D07 keeps the output basis fixed, so its estimation uncertainty is not integrated out. Fast Bayesian FPCA (Sartini, Zhou, Selvin, Zeger and Crainiceanu, arXiv 2412.11340) is one route to modelling it, through the outputs segment's basis contract.

## Moments for piecewise polynomials

A moment service written over piecewise polynomials rather than hinges only gives BAPS's tensor and anisotropic splines, and ABMRS's B-splines, variance decompositions from the same code as BASS's (derivations review §7, item 7).

## Active-subspace analysis

Posterior active-subspace analysis beside the Sobol stage (arXiv 2307.11241, *Discovering Active Subspaces for High-Dimensional Computer Models*, by overlapping authors; [PCA-extensions discussion](../discussion/astra-grassman-stiefel-pca-extensions.md)). For tensor-product representations, including MARS, the gradient integrals reduce to one-dimensional ones. It needs term derivatives, hinge reference cases, numerical integration checks and gradient-accuracy checks before its directions reduce another model, and it comes before jointly learning projection directions during fitting.

## Cross-engine evidence

Evidence for BARS, BASS and BGP by path sampling (derivations review §7, item 8). It needs proper priors in every engine, stated in each engine's units. BASS has them by default since D03 was decided (2026-10-27: σ² scaled to the response, the intercept level-centred), as Normal BARS has. Under D03's evidence rule, evidence is reported from any run that can give one, empirical-Bayes in the level and scale, and the predictive checks already built for OQ17 (WAIC, PSIS-LOO) are the default comparison across engines.

**Held-out runs (added 2026-10-26).** Where runs are held out, proper scores compare any engines that produce predictive draws, with no evidence and no shared prior: the energy score over whole-output draws, and the interval score for marginal intervals (Gneiting and Raftery 2007, as used by Andros et al., arXiv 2406.18751). The energy score rewards the cross-output dependence that D07's joint draws keep, and it is where a missing truncation residual or an independent-noise model shows up. The scorers read draws only and know no model, like the OQ17 consumers.

**Why the engines compare at all (added 2026-10-28).** Farhi, Hed, Bon, Caticha, Mak and Domany (2013, *Physica A* 392, *Temperature Integration*) compute the free-energy difference between two systems by integrating each one's energy over temperature from a common high-temperature limit. That needs the two to share a phase space, so that Z agrees there. Evidence needs less: a normalized prior gives Z(0) = 1 in any dimension, so BARS, BASS and BGP compare without sharing a parameter space. BASS's evidence already has the paper's shape: D05's ladder doubles as the integration points, and its reference rung closes the integral at b = 0 exactly. SPC's notes adopt that closure (`issues/spc/opus-spin-glass-threads.md`, Thread 5).

**Bayes factors under the improper options (added 2026-10-28).** D03 reports a run's evidence only under proper priors, and that stands. A Bayes factor between two runs is a narrower quantity. Liang, Paulo, Molina, Clyde and Berger (2008, *JASA* 103) exclude improper priors on model-specific parameters, whose arbitrary constants stay in the marginal likelihoods (their §1), and put p(α, φ) ∝ 1/φ on the intercept and the precision (1/σ² on the variance) once the columns are centred, so that both are common to every model (§2.1, eq. 3). On common parameters the constants cancel in Bayes factors. Farhi et al.'s energy cutoff is the computational form of that argument: the same truncation in both systems, wide enough not to matter at the temperature of interest, and checked for insensitivity to where it sits.
- Here: truncate the improper priors identically in both runs, run each ladder to its reference rung, take the difference of the evidences, and widen the truncation until the difference stops moving.
- It holds only between models whose improper parameters mean the same thing: the same response, the same intercept convention.
- It also needs the same base model. With improper priors on the common parameters, Bayes factors depend on whether they are null-based or full-based (Liang et al. §2); BASS's centred marginal (`CollapsedGPrior.CentredMarginalAt`) is null-based.
- Cost: a truncated σ² or intercept changes the collapsed integrated form (incomplete-gamma and normal-CDF terms in place of the closed form).
- It is an option, not a gap: D03's proper defaults already give evidence.

## The free-energy curve over the ladder

The ladder's evidence is one point, b = 1, of a curve: ln Z(b) at every rung, the free energy of each tempered posterior π(F)·p_b. Farhi et al. read ΔF at any temperature of interest, and three readouts sit on the curve here.
- **Tempered posteriors.** A posterior with likelihood weight b < 1 is a generalized-Bayes posterior with learning rate η = b (Bissiri, Holmes and Walker 2016, *JRSS B* 78), used for robustness to misspecification. ln Z(η) is not by itself a rule for choosing η; SafeBayes, for one, chooses by sequential predictive loss. What the ladder gives without new sampling is each rung's tempered posterior, so the sensitivity of any readout to η is a view over the ladder once its rungs' draws are kept. This is unrelated to GBMARS, whose "generalized" is the noise family (first entry).
- **WBIC** (Watanabe 2013, *JMLR* 14, eq. 4). The free energy from a single rung, E_b[−ℓ] at b = 1/ln N with ℓ the log-likelihood (−ℓ = nL_n in Watanabe's notation). It is valid for singular models, where BIC is not: its first two terms, −ℓ(w₀) + λ ln N, match the Bayes free energy's (Theorem 4 with β₀ = 1). With the hot end at 1/N that rung lies inside the ladder (at N = 200, 1/N = 0.005 and 1/ln N ≈ 0.19), and the ladder's own ln Z(1) checks it.
- **The real log canonical threshold λ.** At b = β₀/ln N, with β₀ a constant, E_b[−ℓ] = −ℓ(w₀) + λ/b + U_N·√(λ/(2b)) + O_p(1), where U_N has mean zero and is asymptotically Gaussian (Watanabe 2013, Theorem 4). So the slope of E_b[−ℓ] against 1/b between two such rungs converges to λ (Corollary 3), and the ladder's other rungs in that range check it. λ is the Bayesian effective dimension of a singular model. Adaptive splines are singular: a term whose coefficient is zero leaves its knots and variables unidentified, so λ, not the basis count, measures their complexity. Which decision it would inform (calibrating the basis-count prior, comparing engines) is open.

**Needs.**
- D05's ladder with the reference rung, and its rungs' draws kept, beyond the scalar traces evidence reads.
- The integrated form's per-draw slope ∂_b ln p_b at a rung. BASS does not produce it today: in `BassEvidence` only the linear form's slope comes free, and on the integrated form it is infinite at b = 0 under Zellner–Siow, though finite at 1/ln N.
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
- **Why it matters here.** A susceptibility peak whose temperature moves with a control parameter is exactly this situation. It is a second route to families of SPC sweeps, beside ABMRS's peak readouts (next entry).

## Residuals correlated along the output index

From Andros, Guhaniyogi, Francom and Pasqualini (arXiv 2406.18751, 2602.09351), whose models give each run a Gaussian-process residual over the output domain.
- **The problem.** Direct BASS treats noise as independent across output locations. When the misfit is smooth along t, a dense grid then counts n·m independent observations, so intervals are too narrow and births absorb smooth residual. The fitting segment's residual-dependence check detects this; this entry models it.
- **The model.** Residual covariance σ²(I ⊗ R_θ): independent across cases, correlated along the index. Because R acts only on the output side, the direct contractions keep their separable form with Q = R⁻¹ (the fitting segment's output-side inner product). The collapsed marginal gains a log|R| term per case, and the g-prior is taken in the same inner product.
- **Cost.** For a Markov residual on a one-dimensional index, Q is tridiagonal and every output-side product is linear in m. Sampling θ recomputes the output-side Gram for the current terms, O(M² m) with a banded Q.
- **Not covered.** The PCA mode's residual is the empirical truncation residual, whole-vector, so it already keeps within-output dependence (D07); this entry concerns direct fits.

## Direct-functional BASS for families of SPC sweeps

Fitting f(x, T) across families of SPC temperature sweeps (derivations review §7, item 9); BARS remains the engine for a single sweep. For peak and ridge readouts over (T, θ), ABMRS is the better fit (its entry above), and elastic emulation handles peaks that move with θ. An application, outside the completion contract.

Andros et al.'s functional Gaussian process (arXiv 2602.09351) addresses the same problem with the factorization the other way round: a fixed spline basis in the inputs, whose coefficients are GPs over the output domain. BASS-PCA fixes the output side from the data and adapts the input side; direct BASS adapts both, with sparse hinge products. With one or two control parameters a fixed input basis is cheap, but the model's marginal likelihood is a GP over all n·m outputs, so it is a comparison point here, not a route the port prepares for.
