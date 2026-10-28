# BASS port: fitting

A segment of the [full port brief](../bass-full-port-brief.md), which holds the scope and completion contract, the dependency boundaries, the stage table and the decision register. This segment covers stage 4, complete mixed-predictor scalar and direct-functional fitting, and the stage-1 decisions it consumes, D01–D06. It builds on the [numerics segment](port-numerics-brief.md). Its checklist, where the reports go as work lands, is [port-fitting-checklist.md](port-fitting-checklist.md).

D01–D06 were decided in the user's design pass (2026-10-27), and the contract shapes on 2026-10-28; each section says so. The stage-4 plan was taken on 2026-10-28 too. The register, with each decision's status, and the table of sources the drafts were drawn from are in the hub's [decisions section](../bass-full-port-brief.md#statistical-decisions-to-settle).

## Direct functional algebra

Let `A` contain case-side factors and `H` output-location factors. Preserve the R implementation's separability:

\[
G=(A^\top A)\circ(H^\top H),\qquad h_j=A_j^\top YH_j.
\]

The BASS design layer owns construction and caching of `A` and `H`. General numerical contractions may live in LinAlg. Avoid materializing the full `nm × (M+1)` design. Dense contraction with `Y` still depends on the response-grid size.

**Covariates known at the output locations (2026-10-26).** The output coordinates `T: m × d` may carry, beside the coordinates themselves, covariates known at each output location and shared by every case, such as elevation in Andros et al.'s SLOSH studies (arXiv 2406.18751, 2602.09351). They enter as further functional variables, and direct prediction at new coordinates then needs their values there. R's `xx.func` takes several columns the same way, so the papers' statement that BASS cannot use such covariates holds for the PCA mode only, whose basis absorbs them and cannot carry them to new coordinates.

**Residual dependence along the output index (proposal, 2026-10-26).** D07's noise is independent across output locations. When the misfit is smooth along t, as it is for a deterministic simulator, a dense grid then counts n·m independent observations: intervals come out too narrow, and births are rewarded for absorbing smooth residual. Each direct fit reports a residual-dependence check along the output coordinates (for a one-dimensional index, the lag-one autocorrelation of each case's residual under the posterior mean). Modelling the dependence is an [extension](extensions-brief.md#residuals-correlated-along-the-output-index); its prep is in "Room for extensions" below. Andros et al.'s models carry such a residual, a Gaussian process per run.

## Posterior ownership and prediction

- Keep immutable retained structures/draws separate from mutable per-chain or per-particle workspaces. Pooled buffers must never back retained samples after reuse.
- Store optional training predictions separately from posterior state. Intern repeated structures where useful without introducing mutable aliases.
- Preserve requested draw order and explicit draw × case × output axes, including singleton cases.
- Separate conditional-mean estimation from full latent-function draws. Retained uncertainty and nonlinear functionals must include conditional coefficient variation, using valid conditional draws or an independently justified integration method.
- Separate latent function draws, per-case observation noise, and empirical truncation residual vectors. Preserve within-output correlation and apply response scale and residual orientation consistently.
- Fit one input transform and reuse it everywhere. Specify behavior for constant columns, unseen categories, and extrapolation.
- Keep posterior storage separate from checkpoints. Checkpoints include actual RNG state, iteration/stage, adaptation state, replica/particle mapping, and data/config identity.
- Bound component, chain, replica, particle, and numerical-worker parallelism together. Read-only training data may be shared; mutable scratch and RNG streams may not.

Prediction, sensitivity, and summary operations must support weighted posterior draws. Reconstructed output from independent component posteriors must form valid joint draws rather than assuming incompatible component weights align by index.

## Per-observation log density per draw (OQ17)

**Open (2026-10-23): per-observation log density per draw.** Tracked as the preparation brief's OQ17.
- **What.** The observation-model seam that BARS and BASS will share should return log p(yᵢ | xᵢ, draw) for each observation and each retained draw: the n × S pointwise matrix. For functional responses the unit is a case (a whole curve), and its density keeps the within-case correlation.
- **Who uses it.**
  - **Conformal Bayes** (Fong and Holmes 2021). Add-one-in importance weights, wₛ ∝ p(y | x, θₛ), turn the existing draws into conformal prediction sets. Their coverage holds in finite samples under exchangeability, up to Monte Carlo error, with no refits and whether or not the noise model is right. The weights are bounded wherever the density is, unlike leave-one-out's 1/p. Its consumer is the SPC × BARS work, and it is deferred to there.
  - **WAIC and PSIS-LOO** (Watanabe 2010; Vehtari, Gelman and Gabry 2017), for comparing models across BARS, BAPS and BASS.
  - **The noise-model calibration check.** Replicates at one design point are exchangeable, so the grouped (partially exchangeable) form of conformal Bayes tests whether the posterior's predictive width is right there.
- **What exists.** `IObservationModel` returns the collapsed `LogMarginalLikelihood` (coefficients and scale integrated out), conditional means and coefficient draws, and the conditional variance. The chain's state is the knot configuration. No pointwise term exists.
- **Decided 2026-10-23 by the user: both, as configurable options.**
  - **Retained parameters.** Each structural draw keeps its coefficients and noise scale, so the pointwise density is the plain likelihood.
  - **Collapsed predictive.** Each structural draw supplies the collapsed posterior predictive (Student t for the Normal g-prior).
  - **Defaults.** Each engine picks its default at its port (D04, D07), and the other option stays selectable. The seam's default `SampleConditional` returns the conditional mean, which is not a full draw, so a density evaluated there would understate the uncertainty.
- **The prerequisite: a retained posterior sample.**
  - **What BARS does today.** Since 2026-10-24 it keeps every draw's state in `BarsResult.Draws` (`BarsDraw`): the knot configuration, the coefficients and σ² from one joint conditional draw, the robust weights in force at the draw, and the draw's peak readout. The seam returns σ² with β (`SampleConditional`, which replaced `SampleCoefficients`). Before that, `BarsEnsemble` reduced each draw into accumulators fixed at configuration time and dropped the knots, the coefficient draw and its σ², so every new readout needed a rerun (the preparation brief's report "BARS: draw state retained").
  - **What the retained sample holds.** Each draw's generative state: the structure, coefficients, noise scale and, where present, the robust weights.
  - **What is derived, not stored.** Per-draw curves are recomputed from that state. The streaming accumulators stay as online views (the R̂ stop needs them), not as the only output.
  - **Size.** At BARS's sizes this is small: 4 chains × 10,000 draws × about 40 numbers is 1.6 million doubles, about 13 MB.
- **Models without a plain density.** The Gibbs-loss route's term is a pseudo-likelihood, exp(−η·loss), and Whittle's observations are periodogram ordinates. Each states what its pointwise term means, or declines the capability.
- **The unit is the caller's (decided 2026-10-24 by the user).** For SPC × BARS the unit of interest starts as a replicate at a temperature, but the code must treat that as one application, "because the engine shouldn't be built just to do this one thing."
- **Design for the seam (2026-10-24, agreed in direction; built the same day for BARS, `HierarchicalPSpline` and the changepoint draws, `c1adf92`, `7bae205`).**
  - **Two caller-supplied maps, kept apart.**
    - The *unit map* sends each observation to a unit: what is left out or added in. The default is each observation its own unit.
    - The *stratum map* sends each unit to an exchangeability stratum. The default is one stratum. Only the consumers read it; the density does not. Stratified coverage is Mondrian conformal prediction (Vovk); Fong and Holmes's grouped conformal Bayes is its Bayesian form.
    - SPC × BARS sets unit = one replicate's observation at T and stratum = T. Functional data would set unit = a curve; plain regression, unit = an observation with one stratum or covariate bins; a time series, unit = a block.
  - **Which latents are integrated follows from the unit.** A latent is integrated out when everything that depends on it lies inside the unit, and conditioned on otherwise. So the unit chooses the predictive question:
    - `HierarchicalPSpline` with curve units integrates β_r, which asks about a new curve; with observation units it conditions on β_r, which asks about a new point on a seen curve. Both are valid.
    - BARS's robust weight wᵢ belongs to one observation, so it is integrated (the marginal Student t) under any unit containing that observation. Conditioning a unit's density on its own fitted latent would bias leave-one-out.
  - **Models declare their latents' scope:** per observation (wᵢ), per group (β_r, given the group map), or global (β, σ², knots). The result carries the draws × units matrix, the unit map it was computed for, and which latents were integrated and which conditioned on, so a consumer can check the question it answers. A unit map that splits a group latent's dependents is allowed and reported as conditioning, not refused.
  - **Placement.** The density lives with each model: the shared observation-model seam for BARS and BASS, and `HierarchicalPSpline` and `ExactChangepoint` beside it. The contract types (the two maps and the pointwise result) are shared in `Maths.Regression`. The consumers (WAIC, PSIS-LOO, conformal Bayes, the calibration check) read the matrix and the strata and know no model. Only application code, such as SPC × BARS, picks the maps.
  - **The two draw contents fit inside it:** retained parameters integrate only the unit-private latents, and the collapsed predictive also integrates the structure's coefficients and scale.
  - **Built, 2026-10-24.**
    - `UnitPartition`, `PointwiseLogDensity` and `LatentReport` in `Maths.Regression`.
    - The retained-parameters content for BARS (Normal with robust weights integrated through each resampler's marginal, Poisson, Whittle; Gibbs loss declines), `HierarchicalPSpline` (curve units integrate β_r) and the changepoint draws (whole-segment units integrate the level).
    - The model-free consumers `PredictiveChecks.Waic` and `PsisLoo`.
    - BASS's observation model implements the same contract at its port.
  - **Built later the same day:** BAPS's retained terms (`78810a7`); the collapsed content for Normal BARS through `IConditionalFit.LeaveOutLogDensities` (`a6e8e5c`); BGP's collapsed terms, and σ²'s prior on `GPriorConditional.LeaveOutLogDensity`, which BASS's collapsed content will read (`a1f56b1`).
  - **PSIS's relative efficiency** (`e4f4813`): the draws' chains ride on the pointwise result, and PSIS-LOO reads each unit's r_eff for its tail and Monte Carlo error, as loo does. BASS's producer sets its chains the same way.
  - **Not built:** the conformal consumer (SPC × BARS), and the collapsed content for BAPS, `HierarchicalPSpline` and the changepoint draws.

## Starting on the support

- **Chains start on the support (2026-10-24, `6bc5134`).** BARS's overdispersed start sheds knots until the posterior admits it, refuses a run whose zero-knot design is refused, and throws on a retained draw off the support. A chain started off the support can be stranded there for the whole run, so BASS's reversible-jump start over basis functions needs the same rule.

## Exploratory support in births

**Exploratory support in BASS's own births** (added 2026-10-21).
- **What BASS does.** Its births draw knots from observed values. They choose the interaction order with weight `w1` added to the current model's counts. A main-effect birth picks its variable uniformly; an interaction birth weights variables by `w2` added to their counts. Both weights default to 5.
- **Why the pseudo-counts matter.** A death's acceptance carries the reverse birth's probability of the removed term, given the reduced model. With `w1` at 0, an order used only by the removed term would have proposal probability 0 there; with `w2` at 0, so would a variable used only by a removed interaction term. Either death could then never be accepted.
- **What the port does.** It keeps both weights positive and validates them as such, never offering 0.
- **The evidence.** BARS lacked such a component: its births were local only. A 10% uniform component cut its knot count's autocorrelation time 12× and about 20× on two fits (preparation brief, OQ11 and step 3.9).

## Decisions: D01–D06

### D01 (decided 2026-10-27): term representation and prior measure

Decided by the user on 2026-10-27, as drafted. The categorical subset measure follows R, with the paper's uniform over subsets selectable; the prior-cancelling move schedule stays an option for F-M1.

- **A term** is an interaction order K ∈ {1, …, K_max}, a set of K distinct variables, and a sign and a knot for each variable. A model holds an unordered collection of M terms.
- **The measure on the collection** is the prior on labelled terms times the M! orderings. Nott, Kuk and Duc's prior ratio carries (m+1)!/m!, "a term which accounts for the different ways of ordering the basis functions", and a death picks its term uniformly.
  - Net, a birth from M carries p(M+1)/p(M) once and no further factor of M+1: λ/(M+1) under the Poisson prior. R's [birth](D:/aghado01/packages/BASS/R/rjmcmc_des.R:36) (`log(lam) − log(nc)`) and GBMARS's S_B agree.
  - The KHAOS manuscript keeps both the death's 1/(M+1) and the labelled prior's, which puts an extra 1/M! on the count prior: the same family as BARS's R01. That ratio becomes a defect fixture.
- **Repeated terms lie outside the support.** An identical term duplicates a column, so the strict factor refuses it and D06's admission angle is zero. Every admitted collection then has distinct members, which keeps the M! count exact. Interning terms (derivations review §7, item 5) is a storage choice and leaves the measure unchanged.
- **Knot measure: uniform over each variable's distinct observed values.** Francom et al. 2018 build the prior constant from "the number of unique marginal values", and R draws from them ([`unique.ind`](D:/aghado01/packages/BASS/R/bass.R:218), `vars.len`). Each functional variable uses the distinct values in its column of the output coordinates T (the n_r output locations for a one-dimensional index).
  - *Not adopted:* Nott, Kuk and Duc's uniform over the n observations (the 1/n factors in their ratios), which weights tied values by multiplicity; and the continuous Unif(0, 1) written in the JSS and GBMARS texts, which describes the range rather than the measure the code samples.
- **Prior on a term:** K uniform on {1, …, K_max}, the variable set uniform over the C(p, K) subsets, signs uniform, knots as above. Unconstrained, this is Francom et al.'s constant c_K = (1/2)^K · C(p, K)⁻¹ · Πₖ 1/n_{vₖ}, times 1/K_max for the order.
- **Support: the term prior is restricted globally to the feasible set F of D02, without per-order normalizers.**
  - *The sources.* Francom et al. state the per-order version, where c_K is one over the number of admissible order-K terms, and recommend the unconstrained constant "as a conservative proxy" where counting them is infeasible. R uses the proxy ([`logProbChangeMod`](D:/aghado01/packages/BASS/R/mcmc_funcs.R:22)). The global reading is what the implementation samples, stated as a prior.
  - *The tilt.* The order's marginal becomes π(K) ∝ P(admissible | K)/K_max under the unconstrained measure. Products of hinges have more zeros and fail the support count more often, so high orders get less mass than under the per-order version. Report the support-rejection rate per order so the tilt stays visible.
  - *The support rule* is the nonzero-count rule b_m ≥ b of every source (R's `npart`). For direct functional responses the design part and the functional part are counted separately (`npart`, `npart.func`), as Francom et al. suggest. D06's draft keeps the count as the default and offers Q6's column-mass form as an option.
  - *Not adopted as the default:* the per-order normalization. No implementation samples it, and its constant needs a data-dependent count of admissible terms. It stays computable on enumerable fixtures, to measure the tilt.
- **Move-type probabilities depend on the state:** equal probability over the moves available at M (Nott, Kuk and Duc). Birth only at M = 0; death and change ½ each at M_max; ⅓ each in between. A birth's ratio uses d_{M+1}/b_M and a death's b_{M−1}/d_M.
  - R's move *selection* already follows this rule ([selection](D:/aghado01/packages/BASS/R/mcmc_update.R:16)), but its *ratios* use fixed thirds ([constants](D:/aghado01/packages/BASS/R/bass.R:250)). Under a constant likelihood that makes π(1)/π(0) three times the correct value and π(M_max)/π(M_max − 1) two-thirds of it (derivations review §4). It is a defect fixture, not an option.
  - *Option: a prior-cancelling schedule* (Green 1995; used by BARS's DMGK schedule and by He, Yang and Kang). It sets b_M ∝ min(1, p(M+1)/p(M)) and d_{M+1} ∝ min(1, p(M)/p(M+1)), so the count prior cancels from births and deaths; with λ integrated out (D02), p(M+1)/p(M) = (M + h₁)/((M + 1)(h₂ + 1)). BARS found that the best schedule depends on the regime (preparation brief, OQ14), so it is an option to measure in a pre-registered comparison, not a default.
- **Birth proposal** (Nott, Kuk and Duc, extended by BASS to K ≤ 3):
  - the order from I(K | M) ∝ w₁ + (number of order-K terms);
  - for K = 1, the variable uniformly over the p predictors ([R](D:/aghado01/packages/BASS/R/rjmcmc_all.R:21));
  - for K ≥ 2, the variables by weighted sampling without replacement, with weights z_u ∝ w₂ + (appearances of u);
  - signs uniformly, knots from the knot measure.

  A death's reverse probability is that birth's probability for the removed term under the reduced model's counts. w₁ > 0 keeps it positive for every order, and w₂ > 0 for interaction terms.
- **The variable-set probability** is Wallenius' noncentral hypergeometric with one ball per variable: P(S) = ∫₀¹ Π_{i∈S}(1 − t^{αᵢ}) dt, with αᵢ = zᵢ / Σ_{k∉S} zₖ. R's [`dmwnchBass`](D:/aghado01/packages/BASS/R/mcmc_funcs.R:12) evaluates it exactly by inclusion–exclusion over the 2^K subsets; Nott, Kuk and Duc's pair formula is the case K = 2. The port uses inclusion–exclusion up to K = 3 and quadrature of the positive integrand above that, where the alternating sum loses digits.
- **Change move:** pick a term uniformly and one of its factors uniformly, then redraw that factor's sign and knot from their prior. The proposal is symmetric, so the ratio is the likelihood ratio alone (R's [change](D:/aghado01/packages/BASS/R/rjmcmc_des.R:95)).
- **Additional kernel, not the default:** KHAOS's Bernoulli-inclusion proposal for the variable set, with its variable-swap mutation, offered beside the Wallenius proposal once validated ([extensions brief](extensions-brief.md#additional-proposal-kernels)).
- **Considered, not adopted** (the hub's lineage register):
  - *Delayed rejection of structurally invalid births* (KHAOS). KHAOS can condition on its rejected set because that set's probability, a Poisson-binomial, is computable. BASS's births fail on data-dependent support, whose probability is not, so a failed birth counts as a rejection.
  - *A uniform prior over basis functions* (KHAOS), with complexity penalized through per-term shrinkage instead. BASS's prior, uniform over the order and then over the terms of that order, already corrects for how many terms each order has; KHAOS's route is the modified g-prior extension.
- **Independent checks** (F-O3 and F-O4, built 2026-10-27):
  - on a tiny enumerable model (p = 2, a few distinct values, K_max = 2, small M_max) under a constant likelihood, the count law and the order law against enumeration of the globally truncated measure;
  - the two boundary ratios;
  - Wallenius probabilities against brute-force sums over draw orders;
  - defect fixtures for R's fixed thirds and KHAOS's double count factor.
- **Categorical terms** (Francom, Sansó, Bulaevskaya, Lucas and Simpson 2019, §3.3; R's [`genCandBasisCat`](D:/aghado01/packages/BASS/R/rjmcmc_all.R:39) and [`logProbChangeModCat`](D:/aghado01/packages/BASS/R/mcmc_funcs.R:43)):
  - *The term.* A term's categorical part is Π_{l=1}^{L} 1(z_{u_l} ∈ D_l) over L distinct categorical variables, where each D_l is a proper subset of its variable's levels; L = 0 gives 1. A term multiplies a continuous part of order J and a categorical part of order L, each with its own cap (R's `maxInt` and `maxInt.cat`), and J + L ≥ 1.
  - *Variables.* As for continuous variables: uniform when L = 1, and Wallenius with the z weights when L ≥ 2.
  - *Subset measure.* The subset's size is uniform on {1, …, n_v − 1}, then the subset is uniform among those of that size. This is how R samples and scores it, with the proposal and prior factors cancelling, and it mirrors the continuous rule: the order uniform, then the variable set uniform given the order.
    - The paper states that "each combination of parameters conditional on interaction degree has equal prior probability", which is uniform over the 2^{n_v} − 2 proper subsets.
    - The two agree up to three levels, which covers every factor in the paper's application (its 162 combinations are 2 · 3⁴). From four levels on they differ, the paper's version putting more mass on mid-sized subsets.
    - R's measure is the default (decided 2026-10-27), and the paper's uniform is selectable. Both are coherent, so this is a divergence between paper and code, not a correction (the hub's register).
  - *Support.* The combined term needs at least b nonzero points: "each resulting basis function must have at least 20 non-zero points".
  - *Change move.* Redraws one categorical factor's subset from its prior, which is symmetric (R's [`genBasisChangeCat`](D:/aghado01/packages/BASS/R/rjmcmc_all.R:77)).
  - *Levels.* The training factor's levels. A factor with one level carries no information and is refused at setup.
  - *Complement pairs.* With the intercept in every model (D03), a categorical main effect on D and one on its complement span the same space together with it, so they define the same model, and the measure counts each such model twice. This is a tilt to report, not a defect.
- **The functional variables** share their own interaction cap (R's `maxInt.func`), as the hub's scope section requires, and each one's knots are its distinct values.

### D02 (decided 2026-10-27): count prior at the capacity limit

Decided by the user on 2026-10-27, as drafted.

- **Prior:** the derivations review's global truncation (§4) with D01's term measure: π(λ, M, terms) ∝ Gamma(λ; h₁, rate h₂) · Pois(M; λ) · M! · Πₘ π(termₘ) · 1_F, where F is the cap, the support constraints and numerical admission (D06).
- **λ is integrated out**, as in KHAOS. The normalizer is one global constant, so λ's integral stays closed-form: p(M) ∝ Γ(h₁ + M)/(M! (h₂ + 1)^M), and a birth's count factor is (M + h₁)/((M + 1)(h₂ + 1)), restricted to F. λ leaves the chain and the sweep loses a block.
  - For each retained draw, λ | M ~ Gamma(h₁ + M, rate h₂ + 1), which is exact under global truncation. R and GBMARS update λ in this form; GBMARS calls it approximate because its stated prior is the conditionally normalized one.
- **This changes the model the papers state.** JSS and Francom et al. 2018 write the Poisson truncated and renormalized for each λ. Under that prior the Gamma update needs the cap's normalizer, a regularized incomplete gamma function: tractable for the cap alone, intractable once support and admission join F. The draft takes the global reading because it is what R and GBMARS sample, and one mechanism then covers the cap, support and admission.
- **The tilt:** λ's marginal is its Gamma prior tilted by P_λ(F). At R's defaults (h₁ = h₂ = 10, `maxBasis` = min(1000, n), [bass.R:267](D:/aghado01/packages/BASS/R/bass.R:267)) the cap rarely binds. Report the posterior mass at the cap and the support and admission rejection rates.
- **Defaults:** h₁ = h₂ = 10, from R and Nott, Kuk and Duc's Gamma(10, 10). Francom et al. 2018 name λ's prior as the main control on overfitting, and the docs should say so.
- **Independent check** (F-O3, built 2026-10-27): under a constant likelihood on D01's fixture, the count law against the negative binomial restricted to F by enumeration, and the retained λ draws against the Gamma mixture.

### D03 (decided 2026-10-27): coefficient-prior modes and the τ hyperprior

Decided by the user on 2026-10-27, taking the proper-defaults alternative: σ² scaled to the response and the intercept level-centred, so every default prior is proper. The record of that choice, its cost and the options kept selectable is [the subsection below](#d03-the-proper-defaults-decision-2026-10-27).

- **g mode:** a | σ², τ, B ~ N(μ₀, σ² (BᵀB)⁻¹/τ), with τ the precision multiplier as in R, JSS and Francom et al., and μ₀ the intercept's centre (below). σ² ~ IG(c, c·s²) by default (`Sigma2Prior.ScaledToResponse`, c = 10⁻³, s² the response's mean square about its level; c, since a names the coefficients and τ's shape); the Jeffreys p(σ²) ∝ 1/σ² (g₁ = g₂ = 0, R's default) is selectable. R's selectable modes are `g` and `jeffreys`; it also ships a ridge sampler (`rjmcmc_des_ridge.R`) that `bass()` never reaches.
- **τ hyperprior: the unit-information Gamma family,** τ ~ Gamma(a, rate a·n), whose mean is 1/n for every shape a. The default is a = ½, Zellner–Siow: τ ~ Gamma(1/2, rate n/2), equivalently g = 1/τ ~ IG(1/2, n/2) (Liang et al. 2008). The shape is a declared constant, never a function of n: E[ln τ] − ln(1/n) = ψ(a) − ln a at every n, −1.27 at a = ½, and a larger a tightens τ towards unit information, with the fixed-g mode below as its limit a → ∞.

  | Source | Prior as written | Stated intent | E[τ] as written |
  |---|---|---|---|
  | Francom et al. 2018 | Gamma(1/N, 1) | "centered over the unit information prior" | 1/N |
  | R documentation ([b.tau](D:/aghado01/packages/BASS/R/bass.R:34)) | rate "one over the number of observations", shape 1/2 | "centers … on the unit information prior" | n/2 |
  | JSS 2020 | Gamma(1/2, 2/n), "(shape and rate)" | "the Zellner–Siow Cauchy prior" | n/4 |
  | GBMARS | g ~ IG(1/2, 2/n) on the variance multiplier, under its own conjugate update | "leads to the Zellner–Siow Cauchy prior" | n/4 |
  | R code ([default](D:/aghado01/packages/BASS/R/bass.R:287)) | Gamma(1/2, rate 2/n) | none | n/4 |

  Every stated intent centres τ on unit information, and two sources name Zellner–Siow, which does. The written 2/n gives Zellner–Siow only when read as a scale. The code's rate contradicts every stated intent, so the derivations review's rate/scale reading now rests on documented intent rather than inference.
  - *Francom et al. 2018's Gamma(1/N, 1) is the member a = 1/N, corrected* (the hub's [register](../bass-full-port-brief.md#what-the-lineage-taught), 2026-10-27). A shape that shrinks with N keeps the mean at 1/N but spreads ln τ over about N nats: E[ln τ] = ψ(1/N) ≈ −N and sd(ln τ) = √ψ′(1/N) ≈ N, so the median is about 2⁻ᴺ (e⁻²⁸ at N = 40). The prior is centred on unit information in its mean only, which for a precision is carried by a thinning tail. Its effect on inference: under the Liang convention, leaving the intercept-only structure costs 1.4–5.1 nats more than under Zellner–Siow, growing with N, while between two structures that both carry τ the two priors differ by under 0.3 nats. The correction keeps the paper's intent and its family and fixes the shape; fixture F-O8.
  - The family's rate is stated explicitly, since C#'s Gamma takes a scale (derivations review §6).
  - *R's code rate is not ported* ([the parity rule](../bass-full-port-brief.md#validation-and-performance-qualification), 2026-10-27). Gamma(1/2, rate 2/n) contradicts every stated intent, R's documentation included, so it is a defect rather than a choice. It would buy no parity either: τ's prior does not enter the deterministic statistics, which are compared at a fixed τ, and stochastic updates are not matched across languages.
- **n counts observations.** For direct functional responses n = n_x · n_r, as in Francom et al.'s N and R's [`data$n`](D:/aghado01/packages/BASS/R/bass.R:248). Each fit to a supplied-basis or PCA component has n = n_x coefficient responses.
- **Intercept: level centring (`GPriorCentre.Mean`) by default; the Liang et al. convention and R's zero centre selectable.** Under level centring the intercept stays a column inside the g-prior (k = M + 1 in the τ/(1+τ) factor) and the prior is centred on the level, μ₀ = ȳ·e₀: R's convention on y − ȳ. Under the Liang convention (`GPriorCentre.Flat`) the intercept gets a flat prior, the remaining columns are centred and get the g-prior, and k = M. R's zero centre, μ₀ = 0, is DiMatteo, Genovese and Kass's prior, selectable on its own terms (`GPriorCentre.Zero`).
  - *Why.* R keeps the intercept inside the g-prior (k = M + 1). The prior then shrinks the fitted level towards zero, and τ's conditional reads the uncentred fit, whose intercept share n·ȳ² grows with y's level, so τ's posterior and model-size selection depend on where y sits. Under level centring and under the Liang convention τ's posterior is invariant to shifting y (F-O9). (Bayes factors across models are defined under both conventions, since each puts its improper priors only on parameters every model shares; invariance is the difference.)
  - *How much it matters.* A y far from zero relative to its noise pushes τ small under R's convention, and a small τ makes every term's Occam factor, about (τ/(1+τ))^{1/2}, harsher, so fits come out sparser than the data warrant. R's own advice, a commented-out stop in its σ² update, is to centre y. Level centring builds that advice into a proper prior; the Liang convention integrates the level instead, under a flat prior. For supplied-basis and PCA fits the coefficient responses already have zero mean across cases, so there the two conventions nearly coincide; the choice bites for scalar and direct-functional fits.
  - *Computation.* Under every convention the marginal depends only on the column space of [1, Z], so it is computed in the centred basis [1, Z_c], whose Gram is n ⊕ Z_cᵀZ_c: the intercept is orthogonal to every centred column. The centred Gram is the uncentred one minus a rank-one term, n·b̄b̄ᵀ, so the bordered updates and the functional Hadamard form carry over; D06 says where that downdate loses accuracy, and which cancellation-free form the strict rebuild uses instead. The uncentred Gram, R's, serves comparisons with R's deterministic statistics.
  - *It changes the model relative to R,* whose default is the zero centre: decided by the user, 2026-10-27.
- **Propriety by default.** Every default prior is proper: σ² (scaled), the intercept (level-centred), τ (Zellner–Siow) and D02's count prior. So the β = 0 rung of either tempered form, SMC started from the prior, and evidence need no separate mode (D05), and no raw-unit prior remains among the defaults (OQ19 settled for BASS). The improper options, Jeffreys σ² and the Liang flat intercept, give up the absolute marginal and, on the integrated form, the β = 0 rung.
- **Other coefficient priors**, recorded so the contract leaves room for them:
  - *A ridge with a sampled τ* (GBMARS), an option.
  - *A ridge with a fixed scale,* an option only with the scale stated as part of the model: every added term then carries a factor in that scale, so the penalty on model size is set by it (Nott, Kuk and Duc's diffuse N(0, σ²/ε)). KHAOS found its fixed ridge best on noise-free simulators and its g-prior better on most noisy ones, which supports keeping both selectable. Defaults per regime come from a pre-registered comparison of our own.
  - *KHAOS's modified g-prior,* a candidate extension ([extensions brief](extensions-brief.md#khaoss-modified-g-prior)): an n-component g-prior that shrinks each term by its complexity, with per-term factors from its interaction order and degree and a global g₀² under an inverse-gamma prior. Its displayed conditional for g₀² has no data term, so derive it rather than copy it.
  - *Under robust or generalized-hyperbolic noise* (outside the completion contract): the unweighted g-prior as the default, which reduces to Gaussian BASS at unit weights and keeps the latent scales' conditionals in closed form, at one extra factorization per move as in BARS's Q2; the ridge with a sampled τ as the GBMARS-parity option.
- **Independent checks:** τ's marginal for a fixed structure against one-dimensional quadrature; Bayes factors that hold under shifts of y for level centring and the Liang convention and move under R's zero centre, and the level-against-Liang identity (F-O9); invariance to rescaling y under the scaled σ² prior; the complete marginal against quadrature over σ² and τ; the marginal at a fixed τ, under R's intercept convention, against R's deterministic statistics.
- **R's Jeffreys mode has no single target** (derived from R's source, 2026-10-25).
  - *What R does.* `beta.prior = "jeffreys"` holds τ at 0 ([start](D:/aghado01/packages/BASS/R/bass.R:334), [update](D:/aghado01/packages/BASS/R/mcmc_update.R:100)), so the coefficient prior is flat, and a is drawn from N(â, σ²G⁻¹).
  - *Its structure moves* add, for each term, (2π)^{1/2} · σ^{−1} · (|G′|/|G|)^{1/4} ([birth](D:/aghado01/packages/BASS/R/rjmcmc_des.R:39), [death](D:/aghado01/packages/BASS/R/rjmcmc_des.R:80); R's Σ log diag R is ½ log|G|). Its change move adds (|G′|/|G|)^{1/4}.
  - *Its σ² draw,* with a integrated out, is IG(g₁ + (N − 1)/2, g₂ + (yᵀy − Q)/2) ([σ² shape](D:/aghado01/packages/BASS/R/mcmc_update.R:74)).
  - *What coherent priors would give.* A flat prior a ∝ c^k gives c · (2πσ²)^{1/2} · (|G′|/|G|)^{−1/2} per term and a σ² shape of g₁ + (N − k)/2. The Jeffreys prior p(a | σ²) ∝ |G|^{1/2} σ^{−k} gives (2π)^{1/2} per term and no |G| factor. R's factors match neither, and its moves and its σ² draw imply different coefficient priors (flat times σ^{−2k}|G|^{3/4} for the moves, σ^{−(k−1)} for the σ² draw).
  - *So it is not ported as is.* It is an unqualified transition in the sense of the hub's scope section, and becomes a defect fixture.
  - *The root cause is general.* Any improper coefficient prior leaves the posterior over structures undefined, because its normalizing constant changes with the dimension (Bartlett's paradox), so an implementation has to pick an arbitrary constant per term.
- **The second coefficient mode, for the completion contract: fixed g.** The g-prior with τ held at a declared value, unit information τ = 1/N by default, which is BARS's convention and the unit-information family's limit a → ∞. It keeps the Jeffreys mode's practical intent: coefficients are essentially unshrunk (by the factor N/(N+1)), and structure selection is set by the unit-information Occam factor, about ½·log(N + 1) per term, the BIC penalty. It is coherent: §1's collapsed marginal with no τ block.
  - *Not adopted:* a flat prior with a declared per-term constant, coherent but arbitrary, which is an information-criterion penalty in disguise.
- **Posterior propriety** (derived, and checked by quadrature on 2026-10-25).
  - *Reduction.* The structure space is finite: K ≤ K_max, knots on distinct values, finitely many categorical subsets, M ≤ M_max. So the posterior is proper exactly when every admitted structure's marginal p(y | S) is finite.
  - *g mode, R's intercept convention, p(σ²) ∝ 1/σ², Zellner–Siow τ.* Integrating a and σ² leaves ∫ p(τ) (τ/(1+τ))^{k/2} (S(τ)/2)^{−N/2} dτ, with S(τ) = ‖r‖² + τQ/(1+τ). It is finite whenever ‖r‖² > 0. When the structure fits y exactly (‖r‖² = 0), the integrand behaves as τ^{(k−N−1)/2} near zero, so it diverges for k < N and is finite at k = N.
  - *Liang convention.* The same in g = 1/τ: (1+g)^{(N−1−M)/2} (1 + g(1−R²))^{−(N−1)/2} with g ~ IG(1/2, N/2). It is finite for R² < 1 and, when R² = 1, diverges unless M ≥ N − 1.
  - *A proper σ² prior* (g₂ > 0) replaces S(τ) by 2g₂ + S(τ) > 0, and every marginal is finite.
  - *The quadrature.* A 40-point fixture with k = 3 hinge columns spanning 1 and x, the integral taken on log τ with its lower limit moved from 10⁻⁴ to 10⁻⁶⁴:
    - noisy y: converged (log integral 22.77);
    - exact y: grew by 18 per unit of log(1/τ_lower), the predicted rate (N − k − 1)/2;
    - exact y with g₂ = 10⁻³: converged (110.47);
    - k = N: converged;
    - the Liang convention matched, with a slope of 18 in log g.

    The fixture itself hit the cancellation D06 warns about: yᵀy − Q came out slightly negative for the exact fit, so the residual has to be computed directly. These are stage 1's first reference fixtures: [`bass_marginal_oracle.R`](../../../r/oracles/bass_marginal_oracle.R), with every number above reproduced by [`bass_stage1_checks.R`](../../../r/scripts/bass_stage1_checks.R).
  - *Consequence.* Under the Jeffreys σ² option, any admitted structure with k < N that fits y exactly makes the whole posterior improper. Deterministic simulators with piecewise-linear structure can produce one: two opposite hinges at the same knot span x exactly.
  - *Rule.* Under the Jeffreys σ² option, F excludes structures whose residual is zero to working precision (‖r‖ ≤ c·ε·‖y‖, from the residual itself; the oracle uses c = 10³), and reports the refusal. The default's proper σ² makes every marginal finite, so it needs no such rule. A zero response (R's convention) or a constant one (Liang's and level centring's) is refused, as BARS refuses an identically zero response (R13).
  - *The ridge with a sampled τ* has the same structure, and is checked by the same quadrature when it is built.
- **Evidence: where it is reported** (decided 2026-10-27).
  - *From any run that can give one:* the SMC normalizer (stage 8) or a ladder estimate over the tempering rungs (D05's evidence bullet). It is a result field, carried with its Monte Carlo error and with the priors and transforms it depends on, since the marginal likelihood is sensitive to prior scales in a way the posterior is not. The defaults read the level and the scale from y, so it is empirical-Bayes there and comparable only on the same data.
  - *Runs under the improper options report none.* Their priors leave the absolute marginal undefined. Within one fit the posterior over structures already carries the Bayes factors between its models, which is what model-size readouts report.
  - *Comparison across configurations and engines* defaults to the predictive checks built for OQ17 (WAIC and PSIS-LOO), which need no proper prior and no shared prior scales.
  - *Cross-engine evidence* is an extension ([extensions brief](extensions-brief.md#cross-engine-evidence)), since it needs proper priors in every engine.

#### D03: the proper-defaults decision (2026-10-27)

Drafted beside D03 at the user's request and taken by the user in the design pass, 2026-10-27; D03 above now states it. It changed two defaults, σ²'s prior and the intercept's, and kept the rest: τ's family and default, n's meaning, the fixed-g mode, the Jeffreys-mode defect fixture and the other coefficient priors.
- **The choice.** For the intercept, shift-invariance, a proper prior and no double use of the data cannot all hold, since a proper prior on a location is shift-invariant only if it reads the data:

  | Convention | Shift-invariant | Proper | No double use of the data |
  |---|---|---|---|
  | Liang (flat intercept) | yes | no | yes |
  | Level centring (mean), the default | yes | yes | no: about 0.23 nats of extra credit per added term |
  | R's zero centre | no: up to 9 nats across offsets on F-O9's design | yes | yes |

  Level centring was taken for propriety, for one convention across BARS and BASS through the shared block, and because it is R with y centred, as R's own code advises. Liang with the scaled σ² prior was the coherent runner-up.

- **σ²: `Sigma2Prior.ScaledToResponse`, c = 10⁻³,** IG(c, c·s²) with s² the response's mean square about the centre's level. It is Normal BARS's default since 2026-10-27, in the g-prior namespace BASS shares.
  - *Scale-free.* Rescaling y moves every structure's marginal by the same amount, which is OQ19's route (b); BARS checks this from y × 10⁻⁶ to × 10⁶.
  - *Close to Jeffreys.* At c = 10⁻³ the differences between structures sit within 0.05 nats of the Jeffreys prior's on BARS's sine fit.
  - *`Jeffreys` stays selectable,* as the reference prior that reads nothing from y.
- **Intercept: level centring, `GPriorCentre.Mean`.** This is R's convention applied to ỹ = y − ȳ: the intercept stays a column inside the g-prior (k = M + 1), and the prior is centred on the level, μ₀ = ȳ·G⁻¹Zᵀ1, which is ȳ·e₀ for an explicit intercept column. The shared block already takes a prior mean (`CollapsedGPrior`, `GPriorCentre`).
  - *Proper* given σ² and τ.
  - *Shift-equivariant.* ỹ does not move when y shifts, so neither does the posterior over structures and τ.
  - *Mean rather than median under Gaussian noise.* With the mean, the level's own energy is exactly zero (n·mean(ỹ)² = 0), so τ's conditional reads no level at all. `Median`, BARS's default (OQ15), is the robust choice once the noise extensions bring weights.
  - *The Liang convention stays selectable* as the form that is exact in the level but improper, and the zero centre stays selectable on its own terms.
- **What it costs against Liang.** Level centring reads the level from the data at no charge, so σ² gets N degrees of freedom, not N − 1. Under the Jeffreys limit its marginal is Liang's with −N/2 in place of −(N − 1)/2 on (1 + g(1 − R²)), and the intercept-only structure gains the constant ∫p(g)(1 + g)^{−1/2}. Fixture F-O9 (`bass_intercept_oracle.R`) checks that identity on an explicit design to 10⁻⁸, and on F-O8's grid (N = 40, 100, 1000), against Liang's Bayes factors:
  - against the intercept-only structure, level centring favours a structure by 0.43–0.47 nats more at M = 1, 0.67–0.69 at M = 2 and 1.06–1.08 at M = 5;
  - between two structures one term apart (M = 2 at R² = 0.25 against M = 1 at R² = 0.2), by 0.23–0.24 nats more, at every N.

  A term's Occam factor is about ½·log(1 + g), 1.8–3.5 nats at these N, so the extra credit is a small fraction of it, but systematic and not shrinking with N.
- **What it buys.** Every default prior is proper: σ², the intercept, τ under Zellner–Siow, and D02's count prior.
  - *Evidence by default,* wherever a run can give one (a ladder with a reference rung, the SMC normalizer). It carries the priors it read from y (D05's priors record), so it is empirical-Bayes in the level and the scale, and comparable only on the same data.
  - *Exact fits are proper.* The propriety bullet's exclusion of zero-residual structures is needed only under the improper options, and the deterministic-simulator profile no longer has to steer to a proper mode.
  - *The exact reference rung and SMC from the prior* on either tempered form, with no separate proper mode.
  - *OQ19 settled for BASS:* no raw-unit prior remains among the defaults.
  - *Comparisons with R.* Level centring is R's convention on y − ȳ, so R's deterministic statistics compare directly on the centred response. The draft also claimed it spared the centred Gram's downdate. It need not: the marginal depends only on the column space of [1, Z], so the centred basis computes it exactly under every convention (the design pass's correction).
- **Applied in the design pass** (2026-10-27):
  - *D03:* the g mode's σ² default; the intercept bullet (level centring by default, Liang and the zero centre as options); the proper-prior-mode bullet goes, its content now the defaults; the propriety bullet's exclusion applies to the improper options only; the evidence bullet becomes "default runs report evidence where the run can give one", with the predictive checks still the default for comparisons across engines; the independent checks gain shift and scale invariance under the defaults, and the complete marginal against quadrature over σ² and τ.
  - *D04:* σ²'s conjugate draw becomes IG(c + N/2, c·s² + S/2); the chain state is unchanged.
  - *D05:* "Which priors each use needs" collapses, since every use has proper priors by default. The linear form's propriety advantage holds for the improper options only, so the choice between the forms is efficiency alone, the ground on which BARS kept the integrated form.
  - *D06:* unchanged in substance. The draft said the gate would read R's uncentred Gram, stricter by ρ_j = ‖z_j − z̄_j‖²/‖z_j‖²; the design pass corrected that: the marginal depends only on the column space of [1, Z], so the gate reads the centred angles under every convention, and D06 no longer depends on the intercept choice.
  - *The calibration brief:* its paragraph on the improper default goes, since the default then starts a prior-to-posterior SMC run.
  - *The regime note:* its propriety bullet holds for the improper options only.
  - *Fixtures:* F-O9 is built (2026-10-27). On its 40-point design, level centring and Liang hold at offsets 0, 10 and −250, while R's zero centre moves the two-hinge structure's Bayes factor from 5.56 to 2.16 and −3.87 nats.

### D04 (decided 2026-10-27): collapsed and conditional update order

Decided by the user on 2026-10-27, as drafted, with D03's level centring.

- **The rule.** Integrate the coefficients a, and σ² where it is conjugate, out of every block up to their own draw, and draw them immediately before the first block that conditions on them. Every block then samples a conditional of the same joint posterior, and no block reads a statistic from a stale draw. R's τ update is the pattern this rules out: it reads ‖Rβ‖² from a coefficient draw made before σ² was redrawn (the inventory's partially collapsed ordering, derivations review §2).
- **Gaussian BASS.** The chain state is (structure, τ).
  1. Birth, death and change on p(y | structure, τ), with a and σ² integrated out (derivations review §1.2).
  2. τ by slice sampling on log τ, with its Jacobian, each evaluation O(1) from the cached (‖r‖², Q, k, N) (§2).

  No block conditions on a, σ² or λ, so each retained draw gets them by composition: σ² from its collapsed conditional, a given σ², and λ | M (D02). These are exact draws given the chain state, so no argument about scan order is needed.
  - *What changed from §2:* its third block, λ from the conjugate conditional, leaves the chain, because D02's draft integrates λ out. It returns as a draw at retention.
  - Under D03's level centring, k = M + 1 and the statistics are those of y − ȳ, computed in the centred basis. σ²'s draw at retention is IG(c + N/2, c·s² + S/2) under the default scaled prior (D03's c).
- **Retention.** Each retained draw keeps the structure, τ, σ², a and λ: every draw's generative state.
- **Pointwise density default (OQ17).** OQ17's scope rule conditions on latents whose dependents span several units, and a and σ² are global, so the default is the retained-parameters content: the Normal density at each draw's (a, σ²). The collapsed predictive, Student t with a and σ² integrated through the g-prior block's `LeaveOutLogDensity`, stays selectable. BARS's brute-force leave-one-out found its collapsed k̂ pessimistic under clumping (preparation brief), so BASS's collapsed content gets the same check before a comparison relies on it.
- **Robust and generalized-hyperbolic modes** (outside the completion contract; see the [GBMARS extension](extensions-brief.md#generalized-hyperbolic-noise-gbmars)). Per-observation latent scales read the residuals, so a is drawn every sweep, just before them:
  1. the structure given the latents, with a integrated out (and σ² too when the noise is a pure scale mixture);
  2. τ;
  3. a;
  4. the latent scales;
  5. σ², when it is not integrated, and the skew.

  τ stays before a so that it can stay collapsed. Under weights its collapsed density needs both the weighted and the unweighted Gram (D03's unweighted g-prior), so τ given a is the cheap alternative, at the cost of τ–a coupling; the interweaving of `HierarchicalPSpline`'s phase C is the in-house remedy. This is the shape of BARS's robust sweep ([`IWeightResampler`](../../../src/maths/regression/spline/bars/WeightResampling.cs)): knots given the weights with β and σ² integrated, then (σ², β), then the weights.
- **Independent checks:** §2's, for a fixed structure τ's marginal against quadrature and, on a tiny enumerable model, composition draws against a long joint-Gibbs reference; and the retained λ draws against D02's check.

### D05 (decided 2026-10-27): tempering and SMC target paths

Decided by the user on 2026-10-27, as drafted: the likelihood path, untempered by default; the integrated form as the default tempered form, the linear form selectable; the exact reference rung opt-in; SMC on multinomial resampling; stepping-stone evidence.

Here b ∈ (0, 1] is the inverse temperature, as in `Temperature(b, path)`.
- **Default path: likelihood-tempered** (derivations review §3; also BARS's default under Q5). π_b ∝ prior(structure, τ) · p_b(y | structure, τ), where p_b is §3's closed-form tempered marginal with a and σ² integrated out. The count prior, the term measure and τ's prior are untempered at every rung.
- **Why this path:**
  - *Swaps are cheap.* The prior factors are identical across rungs and cancel, so a swap needs only ℓ_b = log p_b, O(1) from cached statistics (§3).
  - *Hot replicas stay meaningful.* They approach the actual prior over structures, where a powered posterior flattens it.
  - *No reference-measure dependence.* A powered density on a space of varying dimension depends on the reference measure of its coordinates. With D01's discrete knot measure the structures are countable, so they carry counting measure on either path, but τ is continuous, and powering its density depends on whether τ or log τ is the coordinate. The likelihood path's b = 0 member is the prior under any parameterization.
  - *An exact reference rung (OQ2).* At b = 0 the target is the prior, which can be drawn exactly: M and the terms from the unconstrained measure, accepted when they lie in F (rejection from the globally truncated measure is exact), and τ from its prior. D03's proper defaults supply it.
- **The powered path stays selectable** (Q5), with τ powered in log τ coordinates and that choice recorded. Its b → 0 limit is its reference measure, which is not normalizable here, so it has no exact reference rung.
- **A second tempered form on the likelihood path, selectable** (decided 2026-10-27; the integrated form stays the default). Beside p_b, the linear form ℓ_b = b·log p₁(y | structure, τ): the collapsed marginal powered as a whole.
  - *What it buys.* Its b = 0 member is the prior over (structure, τ), which is proper even under D03's improper options, since those sit inside p₁. So under those options it keeps an exact reference rung and an SMC start; only the evidence is lost.
  - *BARS evidence* (preparation brief, OQ2's rerun). On Normal BARS the linear form had the smaller communication barrier on two of three fits, by up to 37% at a deep hot end, and no consistent difference on the third. Its efficiency edge, about 1.2×, fell short of the pre-registered threshold, so BARS keeps the integrated default.
  - *BASS.* BASS's knots sit at data values, so its powered path has no useful reference (the counting measure concentrates at the cap), and under D03's improper options only the linear form reaches b = 0.
- **Which priors each use needs.** Under D03's proper defaults every use has them: parallel tempering at any b_min, the b = 0 rung on either form, SMC started from the prior, and ladder evidence. Under the improper options parallel tempering with b_min > 0 still runs, since each rung's collapsed marginal is finite for b > 0 and the improper prior's arbitrary constant appears with both signs in every within-rung ratio and every swap; the b = 0 rung and SMC from the prior then need the linear form, and there is no evidence.
- **Moves at rung b** use ℓ_b(new) − ℓ_b(old), plus the untempered prior and proposal ratios. τ is slice-sampled on log τ from p(τ)·(τ/(τ+b))^{k/2}·(g₂ + b·S_b/2)^{−(g₁ + Nb/2)} (§3).
- **R's tempering is a defect fixture, not an option.**
  - R multiplies the prior and proposal factors of birth and death by the inverse temperature ([birth](D:/aghado01/packages/BASS/R/rjmcmc_des.R:34)).
  - Its swap density ([`lp`](D:/aghado01/packages/BASS/R/mcmc_funcs.R:64)) leaves out the term-prior blocks.
  
  Together these have no single target. Fixture: swap detailed balance on an enumerable structure space.
- **Schedule and ladder.** The non-reversible DEO schedule (step 3.3) and the adaptive ladder are reused; the ladder's adaptation follows D10's policy and freezes after burn-in.
- **Default: untempered.** R defaults to a single temperature (`temp.ladder` = 1). Francom et al. 2018 (§3.3) describe BMARS posteriors as often highly multimodal, but BARS found that tempering helped only where chains were stuck (OQ4). So BASS starts untempered, with tempering a declared option, until its own pre-registered measurement says otherwise.
- **SMC (stage 8)** uses the same sequence, b₀ = 0 < … < b_S = 1, started from the prior.
  - Each particle's incremental weight is ℓ_{b′} − ℓ_b, which costs O(1).
  - Mutations are D04's sweep at the stage's b.
  - Stages are placed adaptively by conditional ESS (Zhou, Johansen and Aston 2016).
  - Resampling uses independent offspring streams.
  
  The normalizer estimate is the evidence under D03's proper defaults, starting from the normalized prior; under the improper options it is not.
- **Evidence (revised 2026-10-27 from the built seam; decided with D05).** `LadderEvidence` (N11) estimates it over the same rungs or SMC stages: stepping-stone by default, and path sampling where the slope is supplied. Four points from reading that seam against this draft:
  1. *The estimator follows the tempered form.* On the integrated form p_b, ∂_b ℓ_b has a closed form (§3), but ℓ_b is not linear in b, so Friel, Hurn and Wyse's discretization diagnostic misses E_b[∂²ℓ_b]: on the seam's nonlinear test target it read 0.094 against a bias of 0.669. Stepping-stone reads ℓ_b only and has no discretization bias, so it is that form's estimator.
  2. *On the linear form, path sampling is free.* Its slope is ℓ₁, already computed for every move, so path sampling runs beside stepping-stone over the same draws, the diagnostic is valid there, and the two estimates' difference reads the trapezoid's error. Rehearsed on BARS (`1fdaacf`): one knot, 10 levels, stepping-stone −1.6313 ± 0.0047 and path sampling −1.6322 ± 0.0049, against a quadrature of −1.6231, with a diagnostic of 0.0010.
  3. *The priors record* (D03's evidence field). BARS's `BarsEvidencePriors` carries the count prior, the support gate, the centre and its resolved level, and the σ² prior as specified and as resolved. BASS's would add τ's prior (D03's shape a and its rate), the term measure (D01), M's count prior (D02) and the interaction cap. Where a prior is resolved from the response (the centre's level, `ScaledToResponse`'s s²), the evidence is empirical-Bayes, and a comparison is licensed only on the same data.
  4. *SMC's reported error under-covers at one sweep per stage.* The Kish-based error covered 0.935 (mean z² 1.14) at five moves per stage and 0.795 (z² 2.15) at one (200 seeds, 1,000 particles): shared ancestry the weights miss. **Built since** (2026-10-27): the genealogy estimate (Lee and Whiteley 2018), calibrated under multinomial resampling at one and five moves (coverage 0.950 and 0.965), and NaN under the other schemes, where it reads too small. Multinomial is `SmcOptions`' default since, so BASS's stage-8 runs report the error unless they choose another scheme (the numerics checklist's report).
- **Independent checks:**
  - the endpoints: b = 1 gives D04's target, and b = 0 gives the prior;
  - swap detailed balance on an enumerable structure space;
  - the SMC normalizer against enumeration with proper priors;
  - exact reference draws against D02's count law;
  - the R defect fixtures above.

### D06 (decided 2026-10-27): numerical model admission

Decided by the user on 2026-10-27, as drafted, with the count as the support default; the gate reads the centred angles under every intercept convention (D03's computation).

- **Two gates, as Q6 decided for BARS** (step 3.7):
  - *Numerical:* the derivations review's angle rule (§5), in an order-free form. A structure is admitted only when every column's angle to the span of all the others satisfies s_j = 1/(G_jj (G⁻¹)_jj) ≥ s_min, read from the centred Gram (below). It is scale-free, and the bordered step's s is its value for the new column.
    - *Which Gram (2026-10-26; settled 2026-10-27).* The marginal depends only on the column space of [1, Z] under every intercept convention, so the gate reads the centred basis [1, Z_c], whose Gram is n ⊕ Z_cᵀZ_c: s_j is the reciprocal of column j's usual, centred, variance inflation factor, and the intercept's own s is 1. R's uncentred Gram, with the intercept as a column, gives the uncentred factor, which serves comparisons with R.
    - *How the two compare.* Both start from the same residual: b_j projected off the span of the other columns and 1. The centred s_j divides it by ‖b_j − b̄_j·1‖², the uncentred one by ‖b_j‖². So uncentred s_j = centred s_j · ρ_j, with ρ_j = ‖b_j − b̄_j·1‖²/‖b_j‖² ≤ 1.
    - *Which is stricter.* The uncentred gate is, most of all for columns whose mean is large against their spread. A full hinge on uniform data has ρ = 1/4; a categorical indicator true for all but a few cases has ρ near zero. Both gates are downward closed, and near s_min they admit slightly different sets. F-O7 has an example: an indicator true for all but one case in 10⁵, nearly collinear with a second column once both are centred. Its centred s is 10⁻⁹, which the Liang gate admits; its uncentred s is 10⁻¹⁴, which R's convention refuses.
    - *R's `getQf` rule,* which reads the uncentred factor's pivots rather than s_j, is not ported; it stays a defect fixture (below).
    - *Checked* (F-O7, 2026-10-26): the identity holds to 3·10⁻¹⁵ on a hinge design.
  - *Support:* the nonzero-count rule b_m ≥ b of D01, counted separately for the design and functional parts (`npart`, `npart.func`). It is a prior support restriction inside F, not a numerical test.

  Exact duplicates and exact collinearity have s = 0 and fail the numerical gate, which D01's rule on repeated terms relies on.
- **Why order-free** (checked 2026-10-25).
  - *The counterexample.* The new column's angle depends on insertion order. With b₁ = e₁, b₂ = e₁ + 10⁻³e₂ and b₃ = e₂ + 10⁻³e₃, the leave-one-out values are (10⁻¹², 10⁻¹², 10⁻⁶). Added last, b₃ shows 10⁻⁶, while b₁ added last shows 10⁻¹². Under a newest-column rule at 10⁻¹⁰, the same set would be admitted or refused depending on its history. F would then not be a set, and the target restricted to it would be undefined.
  - *The gate is downward closed.* Removing a column only widens the others' angles, so a death never fails the gate. A birth or a change needs the whole vector.
  - *Updates in O(k²).* On a birth, with w = G⁻¹u, each old entry of diag(G⁻¹) gains w_j²/δ and the new entry is 1/δ (the bordered inverse). On a death of column m, each entry loses (G⁻¹)_{jm}²/(G⁻¹)_{mm}.
  - *Relation to BARS.* The bordered s is `FactorStatus.MinRelativePivot` (d_j/A_jj) for the new column. Each leave-one-out s_j is at most column j's pivot in any order, so the order-free gate is at least as strict. BARS's columns come in a canonical order (sorted knots), so its `MinRelativePivot` is already a property of the set; BASS's terms are unordered, which is why it needs the leave-one-out form.
- **s_min = 10⁻¹⁰,** the same threshold as BARS's `MinRelativePivot` gate (step 3.7).
  - *Numerics.* A computed s_j carries a relative error of order u_G·κ(E), where E is the equilibrated Gram and u_G is the relative error made in forming it. That is about ε·log₂N when the cross-products over the N rows are summed pairwise, and up to N·ε when they are summed left to right. For designs that pass the gate this stays small at the term counts BASS reaches. `EquilibratedCondition()` reports κ₁(E), and a decision within a factor of two of s_min is confirmed by a strict rebuild.
    - *Measured* (F-O7, 2026-10-26, against QR of the centred design). This replaces the scratch reading recorded earlier the same day, which put the error down to accumulation over the rows in general.
      - Summed left to right, the Gram of a design with a near-constant column gave s ≈ 10⁻⁹ with relative errors of 7·10⁻⁴ to 2·10⁻³, against ε·κ(E) ≈ 2·10⁻⁷. The cause is the summation order: the column's many equal centred values, added after one large one, round the same way every time.
      - Summed pairwise, the same Gram gave 2·10⁻⁷ to 3·10⁻⁷.
      - So the design layer sums its cross-products pairwise or with compensation (the [numerics segment](port-numerics-brief.md#dense-factor-maintenance)). Either way the error stays inside the factor-of-two band.
    - *The centred Gram's accuracy.* D03's rank-one downdate, BᵀB − N·b̄b̄ᵀ, cancels for columns with small ρ_j: their diagonal entry loses about ε/ρ_j relative, however the cross-products are summed. F-O7's numbers, both forms summed pairwise:
      - an indicator true for all but one case in 10⁵ (ρ = 10⁻⁵): relative error 3·10⁻³ from the downdate, against 3·10⁻⁷ from explicitly centred columns, so still inside the band;
      - one true for all but a thousand (ρ = 10⁻²): the two agree, at 4·10⁻⁸ and 2·10⁻⁷.
    - *What uses which.* Incremental updates may downdate. The strict rebuild forms the centred Gram without cancellation, under every convention:
      - scalar fits use explicitly centred columns;
      - direct fits use the centred Hadamard form G_c = (A_cᵀA_c)∘(H_cᵀH_c) + n·(āāᵀ)∘(H_cᵀH_c) + m·(h̄h̄ᵀ)∘(A_cᵀA_c). A_c and H_c are the centred case-side and output-side factors, and ā and h̄ their means. It is a sum of three positive semidefinite terms, so its diagonal forms without cancellation, and it keeps the separable structure. In F-O7, with factors whose means are large against their spread, plus one design-only and one functional-only term, it matched explicit centring to 2·10⁻¹⁶ in the Gram and 3·10⁻¹⁶ in s. The downdate gave 4·10⁻¹³ and 1·10⁻¹³.
  - *Statistics.* 1/s_min = 10¹⁰ bounds every coefficient's variance inflation. It is a numerical guard, not a modelling restriction; the marginal's Occam factor already discourages redundant terms.
  - *Option:* a declared cap on variance inflation (`MaxVarianceInflation`), a statistical support restriction for stable extrapolation. Off by default.
- **R's rule: a defect fixture, not ported.** R's `getQf` refuses a design when max_{j≥2} r_jj / min_j r_jj > 10³ ([`getQf`](D:/aghado01/packages/BASS/R/mcmc_funcs.R:118), commented "a hack, otherwise we get huge variance inflation in beta"). That is a pivot-ratio gate at 10⁻⁶ on the squared pivots.
  - It depends on column norms. Like BARS's `PivotRatio` (R15; step 3.7 counted 1,590 refusals against 76 under the scale-free gate), it refuses well-conditioned local terms whose columns are small, and does so more often as n grows.
  - Fixture: a term orthogonal to the rest, with b nonzero points, refused by R's rule at large n and admitted by the angle gate. At n = 40,000, a 20-point local term with values 0.01 gives R's ratio 1291, so R refuses it, while the smallest leave-one-out angle in the design is 0.25.
- **Reference fixture:** [`bass_admission_oracle.R`](../../../r/oracles/bass_admission_oracle.R) emits the leave-one-out angles, the pivots in the given order and R's `getQf` decision. Its default mode reads the uncentred Gram, with the intercept as a column. Its other modes (F-O7) cover the centred angles the gate reads:
  - `centred` gives the centred angles by QR, ρ and the uncentred angles, and the centred Gram formed explicitly (summed pairwise and left to right) and by the downdate;
  - `hadamard` gives the same for direct fits from A and H, with the three-term form. The order counterexample and the large-n fixture are reproduced by [`bass_stage1_checks.R`](../../../r/scripts/bass_stage1_checks.R).
  - It is not kept as an option for comparisons with R's posterior ([the parity rule](../bass-full-port-brief.md#validation-and-performance-qualification)). Those comparisons could decide nothing, since R's own chain is off its target: the fixed ⅓ move probabilities alone put it 0.020 off in total variation on F-O6's enumerable space.
- **Count against column mass.** For BARS, Q6 replaced the count with column mass Σᵢ Bⱼ(xᵢ). The count is BASS's default (decided 2026-10-27): every BASS source uses it, and Francom et al. 2018 justify it by the edge instability of hinge terms with few nonzero values. Column mass, Σᵢ B_m(xᵢ) with hinges normalized to a maximum of one, is an option, and a fixture compares the predictive variance of edge terms under the two.
- **A strict rebuild before refusal.** When an incremental update fails the numerical gate, a strict rebuild decides before the candidate is refused (the [numerics segment's](port-numerics-brief.md#dense-factor-maintenance) dense factor maintenance). Refusals are reported by reason: numerical, support, cap.
- **A ridge is a prior, never a numerical repair.** Under the g-prior a rank-deficient design lies outside F; a ridge enters only as one of D03's declared coefficient priors.
- **Shared state.** The residual r and δ serve the factor border, the numerical gate and informed proposals (§5). r is updated on each accepted birth, recomputed after deaths and on a cadence, and S is computed from r rather than as yᵀy − Q/(1+τ), whose cancellation bites when τ is small (§5).
- **Independent checks:** §5's, randomized insert and delete sequences against fresh factorizations, residual drift against recomputation and s against the SVD's principal angle; plus one constructed case for each refusal reason.

## Contract shapes (decided 2026-10-28)

Stage 1's last deliverable: the data, config and result contracts for fitting, drafted after D01–D06 were decided and consuming them. The user took the recommendations on 2026-10-28; the choices are recorded at the end of the section. They follow CONTRIBUTING's strict core: one immutable config record, a static entry point, and delegates and live objects on the call, never in the config. Placement is the hub's: `Maths.Regression.Spline.Bass`, under `src/maths/regression/spline/bass/`.

**Entry points.**

```csharp
public static class Bass
{
    public static BassFit Fit(BassConfig config, BassData data, BassRunOptions? run = null);
    public static BassFit Resume(BassCheckpoint checkpoint, BassConfig config, BassData data, BassRunOptions? run = null);
}
```
- `BassRunOptions` carries what is not the model: one bound on parallelism over chains, replicas, particles and workers (the posterior-ownership rules above), a cancellation token, a checkpoint request and progress callbacks. It is neither serialized nor fingerprinted.
- `Resume` refuses a checkpoint whose fingerprint differs from that of `config` and `data`.

**Data.**

```csharp
public sealed record BassData(IReadOnlyList<BassPredictor> Predictors, BassResponse Response);

public abstract record BassPredictor(string Name);
public sealed record ContinuousPredictor(string Name, double[] Values) : BassPredictor(Name);
public sealed record CategoricalPredictor(string Name, int[] Codes, IReadOnlyList<string> Levels) : BassPredictor(Name);

public abstract record BassResponse;
public sealed record ScalarResponse(double[] Y) : BassResponse;
public sealed record FunctionalResponse(double[,] Y, IReadOnlyList<ContinuousPredictor> OutputVariables) : BassResponse;   // Y is n × m
```
- *Predictors are named columns* in the caller's order, which the fit keeps, so results, sensitivity indices and predictions refer to them by name and position. Mixed types are first-class, and categorical codes index `Levels`.
- *A functional response* shares one set of m output locations across cases (the hub's representation table). The output variables describe them: the coordinates, and any covariates known at each location (the direct functional algebra section), each a further functional variable.
- *Refused at entry, with a reason:* a constant continuous column, which has no knot; a factor with one level; a non-finite response; a zero or constant response (D03); mismatched lengths.
- *Supplied-basis and PCA responses are not a fitting input.* The outputs segment fits each coefficient response as a `ScalarResponse` and keeps the basis beside the component fits (stage 5).
- *Room for the noise extensions:* an offset and observation weights join `BassData` with them (F-X1), and are absent now.

**The input transform** (`BassInputTransform`), fitted once from the training data and kept in the fit:
- continuous predictors and output variables: R's range scaling to [0, 1], stored as each column's (min, max); knots are the scaled distinct values (D01);
- categorical predictors: the training levels, by name.

At prediction, a value outside the training range extrapolates the hinges linearly and is flagged per case. An unseen level is refused by default, or placed outside every subset (`UnseenLevel.OutsideEverySubset`), declared on the prediction call.

**Config.** Every field is a value, so the record serializes to JSON and is its own manifest.

```csharp
public sealed record BassConfig
{
    public BassStructureSpec Structure { get; init; } = new();                       // D01, D06
    public BassCountPrior Count { get; init; } = new();                              // D02
    public BassCoefficientPrior Coefficients { get; init; } = new SampledG();        // D03
    public Sigma2Prior Sigma2 { get; init; } = Sigma2Prior.ScaledToResponse(1e-3);  // D03
    public GPriorCentre Intercept { get; init; } = GPriorCentre.Mean;                // D03
    public BassSampling Sampling { get; init; } = new();                             // D04, D05, D10
}

public sealed record BassStructureSpec
{
    public int Degree { get; init; } = 1;
    public int MaxInteraction { get; init; } = 3;
    public int MaxCategoricalInteraction { get; init; } = 3;
    public int MaxFunctionalInteraction { get; init; } = 3;
    public int? MaxBasis { get; init; }                                  // null: min(1000, N − 1)
    public SupportRule Support { get; init; } = new();
    public double AngleFloor { get; init; } = 1e-10;                     // D06's s_min
    public double? MaxVarianceInflation { get; init; }                   // off
    public CategoricalSubsetMeasure CategoricalSubsets { get; init; } = CategoricalSubsetMeasure.SizeUniform;
    public double OrderWeight { get; init; } = 5;                        // w1, kept positive
    public double VariableWeight { get; init; } = 5;                     // w2, kept positive
}
public sealed record SupportRule(SupportKind Kind = SupportKind.Count, double? Minimum = null, double? MinimumFunctional = null);
                                                                         // null: min(20, 0.1·n) per part, R's default
public sealed record BassCountPrior(double Shape = 10, double Rate = 10);  // h1, h2; λ integrated out

public abstract record BassCoefficientPrior;
public sealed record SampledG(double TauShape = 0.5) : BassCoefficientPrior; // τ ~ Gamma(a, rate a·n); a = ½ is Zellner–Siow
public sealed record FixedG(double? Tau = null) : BassCoefficientPrior;      // null: unit information, 1/n
// Ridge (GBMARS's sampled-τ ridge) and FixedRidge (its scale stated as part of the model) join when built (stage 4's choice 1)

public sealed record BassSampling
{
    public int Chains { get; init; } = 4;
    public int MasterSeed { get; init; }
    public int BurnIn { get; init; } = 9000;                 // R's split of 10,000; provisional
    public int Samples { get; init; } = 1000;
    public int Thin { get; init; } = 1;
    public double RHatTarget { get; init; }                  // 0: no early stop
    public int BatchSize { get; init; }
    public double MinEssPerChain { get; init; } = 100;
    public MoveSchedule Moves { get; init; } = new EqualOverAvailable();   // or PriorCancelling(C = 0.4)
    public TemperingSpec? Tempering { get; init; }           // null: untempered (D05); joins with S8
    public SmcSpec? Smc { get; init; }                       // null: MCMC; joins with stage 8
    public bool ReportEvidence { get; init; } = true;
}
```
- `TemperingSpec` holds D05's choices: the path (likelihood by default, powered selectable), the tempered form (integrated by default, linear selectable), the ladder (levels, hot end, swap schedule, adaptation frozen after burn-in, D10) and the exact reference rung (opt-in). `SmcSpec` holds the particle count, the resampling scheme (multinomial by default), the adaptive schedule's conditional-ESS target and the sweeps per stage.
- *Why values, not component instances.* `BarsConfig` holds instances (`IObservationModel`, priors, kernels), so it neither serializes nor fingerprints whole; that is how BARS's checkpoint identity came to cover component types and not their parameters (OQ10). BASS's choices are closed sets that the decisions enumerate, so values cover them, the manifest is the config, and the fingerprint is a hash of its serialization and of the data. A new prior or kernel joins its closed set with its decision.
- *Auto-picks* (the basis cap, the support counts, n, s² and the centre's level) are resolved at entry and reported in the fit's `Resolved` record, the automation log CONTRIBUTING asks for.
- *Not carried:* R's `s2.lower` (a truncated σ² prior, which would have to be stated as one) and `beta.prior = "jeffreys"` (a defect fixture, D03). R's `save.yhat` and `small` become a view and the draw store's thinning.

**Result.** The retained posterior and its account, immutable.

```csharp
public sealed record BassFit(
    BassConfig Config, BassResolved Resolved, BassInputTransform Transform, RunIdentity Identity,
    BassTermTable Terms,                                 // interned terms
    IReadOnlyList<DrawStore<BassDraw>> Draws,            // per chain, in draw order
    double[]? DrawWeights,                               // SMC; null: equal
    BassDiagnostics Diagnostics, BassEvidence? Evidence, ResidualDependence? ResidualDependence,
    EnsembleCheckpoint? Checkpoint);

public sealed record BassDraw(BassStructure Structure, double Tau, double Sigma2, double[] Coefficients, double Lambda);
public readonly record struct BassStructure(int[] TermIds);   // sorted ids into Terms: an unordered collection (D01)
```
- *Every draw's generative state* (D04): the structure, τ, σ², the coefficients and λ. Fitted values, posterior means and intervals are views over it, computed on demand. Training fits are a cached view only on request (R's `save.yhat`), stored apart from the draws.
- *Interned terms.* A `BassTermTable` shared by every draw holds each term's continuous factors (variable, sign, knot), categorical factors (variable, level subset) and functional factors (output variable, sign, knot), with the fit's degree. Structures hold ids, so a retained draw never aliases mutable state.
- *Weighted draws.* SMC's draws carry their weights, and every view (prediction, sensitivity, summaries) takes weights, equal ones for MCMC.
- *Diagnostics* (`BassDiagnostics`):
  - rank-normalized R̂ and ESS for the monitored functionals: M, τ, σ², and fitted values at declared points;
  - acceptance per move type;
  - refusals by reason (`Tally<BassRefusal>`: cap, support, numerical; D06);
  - the support-rejection rate per order (D01's tilt) and the posterior mass at the cap (D02's);
  - the tempering readout when tempered, and SMC's stage reports with the genealogy error when SMC ran.
- *Evidence* (`BassEvidence`): stepping-stone, and path sampling where the slope comes free, with their Monte Carlo errors, and `BassEvidencePriors` (D05's point 3): the coefficient prior and τ's shape, the count prior, the term measure and interaction caps, the support rule, the centre and its resolved level, and σ²'s prior as specified and as resolved. Null under the improper options.
- *The checkpoint* rides beside the posterior, never inside it.
- *Derived operations* on the fit, each a function of the draws: `PointwiseLogDensity(UnitPartition, PointwiseContent)`, the shared OQ17 seam, with retained parameters by default (D04); `Predict` (stage 5, the outputs segment); `Sensitivity` (stage 6).

**The choices** (decided by the user, 2026-10-28, as recommended).
1. *Closed sets as records, values throughout,* not `BarsConfig`'s instances, so the fingerprint covers the whole config.
2. *The intercept enum.* The shared `GPriorCentre` gains `Flat` for the Liang convention, one shape per pattern; BARS refuses it until the shared block supports the Liang form.
3. *σ²'s prior as a value.* The shared `Sigma2Prior` becomes a serializable record, which BARS's config then shares.
4. *Sampling defaults.* R's 9,000 burn-in and 1,000 kept draws, provisional until BASS's own pre-registered measurement, as BARS's are (OQ14).
5. *Unseen levels at prediction* are refused by default; `UnseenLevel.OutsideEverySubset` is declared on the prediction call.

**Refined in S1** (2026-10-28, building the contracts):
- `GPrior` is `SampledG`: the name collides with the namespace `Maths.Regression.GPrior`, which code in `Maths.Regression` resolves first. The pair reads as the g-prior with τ sampled or fixed.
- The support thresholds are `Minimum` and `MinimumFunctional`, doubles: R's own default, min(20, 0.1·n), is fractional, and under the column-mass option the threshold is a mass.
- The data's arrays are copied at construction and exposed read-only (`ReadOnlyMemory<T>`), so a fit, its fingerprint and its draws never see a value change; the constructors still take arrays.
- A setting joins with the code that reads it: `Tempering` with S8, `Smc` with stage 8, the ridge priors when built.
- The median centre is refused with its reason: D03 offers level centring, Liang and the zero centre, and keeps the median for the noise extensions' weights.
- The basis cap defaults to min(1000, N − 1), since the intercept counts toward the dimension bound (the hub's scope); R's min(maxBasis, n) allows one column too many.
- A declared level no training case takes is left out of the training levels and logged; one-level factors are judged on the levels the cases take.
- The unions' discriminator is `Kind`, as `Sigma2Prior`'s field is, so the manifest reads one way throughout.

Choices 2 and 3 changed shared code before stage 4's first fit (F-C2, landed 2026-10-28): `Sigma2Prior` is a record of its kind and parameters, which BARS's Normal model holds, and `GPriorCentre.Flat` exists, refused by BARS with its reason.

## Stage 4 plan (decided 2026-10-28)

Drafted after the baseline was requalified at `cf37584`, and taken by the user the same day with its four recommendations (below). It sequences F1–F11 into ten steps, each landing with its tests and a report in the checklist, and names what each reuses, what it adds to shared code and which fixtures check it. The user's go-ahead authorizes its source code, S1 first.

**What it builds on.**
- *Factors and contractions:* `BorderedCholesky` (N1–N4), an `ISpdFactor`, so `CollapsedGPrior.FromLeastSquares` takes it and `MarginalAt` (N12) gives the marginal at any (τ, b) in O(1); `KhatriRao` (N7); `PairwiseSum` (N10).
- *Priors:* `Sigma2Prior` and `GPriorCentre` (F-C2); `Gamma`.
- *Sampling:* `ReversibleJumpChain`, which owns the state-dependent selection ratio, so D01's schedule is move weights (`ScheduledMove`); `BlockScan` and `SliceSampler`; `ParallelTempering` with the reference rung and `TraceRungs`; `LadderEvidence` (N11).
- *The run:* the chain ensemble, `DrawStore`, `RunIdentity` and checkpoints; `UnitPartition` and `PointwiseLogDensity` (OQ17).

**Shared-code additions** (each generic, tested in its own host, landed before its first BASS use):
- *Wallenius' multivariate probability* with one ball per variable, in Distributions: inclusion–exclusion up to K = 3 and quadrature of the positive integrand above (D01), checked against F-O5.
- *A per-replica block in `ParallelTempering`:* a kernel that leaves π_b invariant, run on each replica at its own b between swap rounds. A replica takes reversible-jump moves only today, so τ's slice step (D04) has no seat under tempering. Tested on an analytic target with a continuous coordinate.
- *The flat-intercept form in `CollapsedGPrior`* (Liang: k = M, N − 1 degrees of freedom), for `GPriorCentre.Flat`. BARS keeps refusing it, for its basis.

**Steps.**

| Step | Delivers | Checklist | Checks |
|---|---|---|---|
| S1 | Contracts: `BassData` and its refusals, `BassInputTransform`, `BassConfig` serialized with its closed sets as polymorphic records, the fingerprint, `BassResolved` | (F-C1's shapes) | JSON round trip; each refusal; the transform at the training range's edges |
| S2 | Terms: hinge factors of any positive degree, categorical subsets, functional factors; the interned term table; column evaluation, case side A and output side H; support counts per part | F1, F2 (in part) | basis evaluation against R (a new oracle mode); a tiny functional design's explicit nm × k columns against the separable form |
| S3 | D01's term measure and proposals with D02's count factor: births (w₁, w₂, Wallenius), deaths with the reverse birth's probability, the change move; R's categorical subset measure, the paper's uniform selectable | F1, F4 (in part) | proposal probabilities summing to one over enumerable spaces; Wallenius against F-O5 |
| S4 | The design layer and admission: the centred Gram (explicitly centred columns for scalar fits, the three-term Hadamard form for direct ones), pairwise cross-products, bordered updates, the order-free angle gate with a strict rebuild near s_min, the support gate, refusals by reason, the residual on its cadence; column mass and `MaxVarianceInflation` as options | F2, F7 | randomized insert and delete sequences against fresh factorizations; F-O7; F-O2's order counterexample; one constructed case per refusal reason; residual drift |
| S5 | The collapsed target: level centring (default), the zero centre and Liang; the g and fixed-g modes; σ² scaled (default), Jeffreys with the exact-fit exclusion, and a fixed IG; the integrated and linear tempered forms on the likelihood path, and the powered path | F3 | F-O1's quadrature; F-O9's identities and shift invariance; scale invariance; F-O8; R's deterministic statistics at fixed τ under its zero centre (parity) |
| S6 | D04's sweep: the reversible-jump palette on D01's schedule (prior-cancelling selectable), τ's slice block on log τ, composition at retention (σ², a, λ), `BassDraw` | F4, F5 | F-O3's and F-O4's laws from the chain under a constant likelihood, with R's thirds asserted corrected; an enumerable posterior with data (a new oracle mode); composition draws against a long joint-Gibbs reference; λ against the Gamma mixture |
| S7 | The run: `Bass.Fit` and `Bass.Resume` on the chain ensemble, each chain's draws in a `DrawStore`, starts on the support, the diagnostics (R̂ and ESS for M, τ, σ² and declared fitted values; acceptance per move; refusal tallies; support rejection per order; mass at the cap), the R̂ stop | F8, F10 | resume bit-identical; results independent of parallelism; every start in F |
| S8 | Tempering: the per-replica block, the likelihood-path ladder (integrated form by default, linear selectable), the exact reference rung (opt-in), DEO swaps and the adaptive ladder frozen after burn-in, the evidence with `BassEvidencePriors` | F6 | the endpoints b = 0 and 1; swap detailed balance on an enumerable space, with R's swap density asserted corrected (F-O6); the ladder's evidence against enumeration under proper priors |
| S9 | The pointwise density (retained parameters by default, the collapsed predictive selectable; case units for direct fits; chains for r_eff) and the residual-dependence readout | F9, F11 | against dense computation; WAIC and PSIS-LOO read it |
| S10 | Qualification: controlled recovery on simulated surfaces (scalar continuous, mixed categorical, direct functional), a cost baseline by n, p and M, the docs section with the decisions and their trade-offs, the suites | | recovery within the posterior's intervals; R's posteriors as exploration only (the parity rule) |

The order is S1 → S2 → S3 and S4 → S5 → S6 → S7 → S8 → S9 → S10. The shared additions land with S3 (Wallenius), S5 (the flat-intercept form) and S8 (the per-replica block).

**Starts.** Each chain starts from an exact draw of the prior restricted to F, the sampler the reference rung uses (D05). The starts are then dispersed, which R̂ needs, and on the support by construction (F10). R starts every chain at M = 0.

**New oracle modes** (R, under `r/oracles`, reproduced by a stage-4 checks script):
- basis evaluation and the deterministic statistics at fixed τ, against R's own functions on the pinned checkout (parity where R is right);
- `posterior`: the exact posterior over an enumerable structure space with data, τ integrated by quadrature;
- the ladder's evidence on that space, against its normalizer.

**The choices** (taken by the user, 2026-10-28, as recommended):
1. *The ridge priors are deferred.* Stage 4 builds the g and fixed-g modes, which the completion contract asks for. `Ridge` and `FixedRidge` join the closed set when they are built with their own marginal, which needs |G + τI| rather than the g-form and a refactoring per τ.
2. *τ under tempering takes the per-replica block,* which keeps D04's slice step and serves any hybrid replica. The alternative, τ as a Metropolis move inside the palette, would have changed D04's block.
3. *Chains start from exact prior draws,* not R's M = 0.
4. *One sequence runs the steps,* in the session that planned them.

**Not in stage 4.** F-M1 and F-R1, the pre-registered comparisons, need the running engine and follow it. Supplied-basis and PCA fits, and prediction at new inputs, are stage 5.

## Regime note: deterministic simulators (proposal, 2026-10-26)

Three findings from the lineage and the propriety check point the same way.
- **KHAOS:** its fixed ridge did best on noise-free simulators, its g-prior on most noisy ones.
- **GBMARS's review of BMARS:** a_λ and b_λ smaller than the usual 10 suit deterministic surrogates, which need more terms to interpolate.
- **D03's propriety result:** under the improper σ² prior, a structure with k < N that fits y exactly makes the whole posterior improper, and deterministic data are where exact fits happen. D03's default σ² prior is proper, so this binds only under the Jeffreys option.

**Proposal:** a named profile, `DeterministicSimulator`, which is not a default:
- D03's default proper σ² prior, scaled to the response, which needs no standardization;
- a wider λ prior (h₁ and h₂ well below 10, their values recorded);
- the ridge with a sampled τ as the alternative to the g-prior;
- in direct functional fits, the residual-dependence check read before any interval is reported, since a deterministic simulator's residual is all misfit and so is smooth along t.

Which coefficient prior the profile uses is left to the pre-registered comparison, not set from KHAOS's study. Under the Jeffreys option the exact-fit refusal stays, and points users to this profile.

## Room for extensions

Prep work in this segment, so that the [extensions](extensions-brief.md) plug in without reworking the port:
- **An observation-model seam with an offset.** BASS's observation model hands the Gaussian core a working response and weights, (y − offset, W), not weights alone. Gaussian BASS uses a zero offset and unit weights; generalized-hyperbolic noise (GBMARS) then needs only its latent-scale block.
- **A latent slot in the sweep.** D04's robust order (structure, τ, a, latent scales, then σ² and the skew) is the sweep's general form; Gaussian BASS leaves the latent slots empty.
- **Term types behind one design-kernel contract** (derivations review §7, item 6): hinge products, categorical subsets and the functional variable now, and later tensor-product B-splines (ABMRS) and directional hinges (learned projections). Each term type supplies its column evaluation, its support count, its moments and its derivatives.
- **Proposal kernels behind one contract,** each with exact forward and reverse probabilities. The Wallenius variable proposal is the first member; KHAOS's Bernoulli-inclusion proposal, the variable swap and a history-based proposal (D10's policy) join it later.
- **Coefficient priors as a declared set:** g with a sampled τ, fixed g and a ridge with a sampled τ now, so KHAOS's modified g-prior is one more member.
- **Pointwise densities per observation,** with per-observation latents integrated under OQ17's scope rule, so the generalized-hyperbolic density is one more case.
- **A declared inner product on the output side.** The direct contractions are written G = (AᵀA)∘(HᵀQH) and h_j = A_jᵀYQH_j, with Q = I now. A residual correlated along the index, with covariance σ²(I ⊗ R) and Q = R⁻¹, then changes only the output-side factors, plus Σᵢ yᵢᵀQyᵢ and a log|R| term per case in the collapsed marginal. The case side and the separable algebra stay as they are. A Markov residual on a one-dimensional index has a tridiagonal Q, which keeps the cost linear in m.
