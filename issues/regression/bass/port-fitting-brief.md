# BASS port: fitting

A segment of the [full port brief](../bass-full-port-brief.md), which holds the scope and completion contract, the dependency boundaries, the stage table and the decision register. This segment covers stage 4, complete mixed-predictor scalar and direct-functional fitting, and the stage-1 decisions it consumes, D01–D06. It builds on the [numerics segment](port-numerics-brief.md). Its checklist, where the reports go as work lands, is [port-fitting-checklist.md](port-fitting-checklist.md).

The decision drafts here are proposals for the user's design pass, not decisions. The register, with each decision's status, and the table of sources the drafts were drawn from are in the hub's [decisions section](../bass-full-port-brief.md#statistical-decisions-to-settle).

## Direct functional algebra

Let `A` contain case-side factors and `H` output-location factors. Preserve the R implementation's separability:

\[
G=(A^\top A)\circ(H^\top H),\qquad h_j=A_j^\top YH_j.
\]

The BASS design layer owns construction and caching of `A` and `H`. General numerical contractions may live in LinAlg. Avoid materializing the full `nm × (M+1)` design. Dense contraction with `Y` still depends on the response-grid size.

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

## Decision drafts: D01–D06

### D01 (draft): term representation and prior measure

- **A term** is an interaction order K ∈ {1, …, K_max}, a set of K distinct variables, and a sign and a knot for each variable. A model holds an unordered collection of M terms.
- **The measure on the collection** is the prior on labelled terms times the M! orderings. Nott, Kuk and Duc's prior ratio carries (m+1)!/m!, "a term which accounts for the different ways of ordering the basis functions", and a death picks its term uniformly.
  - Net, a birth from M carries p(M+1)/p(M) once and no further factor of M+1: λ/(M+1) under the Poisson prior. R's [birth](D:/aghado01/packages/BASS/R/rjmcmc_des.R:36) (`log(lam) − log(nc)`) and GBMARS's S_B agree.
  - The KHAOS manuscript keeps both the death's 1/(M+1) and the labelled prior's, which puts an extra 1/M! on the count prior: the same family as BARS's R01. That ratio becomes a defect fixture.
- **Repeated terms lie outside the support.** An identical term duplicates a column, so the strict factor refuses it and D06's admission angle is zero. Every admitted collection then has distinct members, which keeps the M! count exact. Interning terms (derivations review §7, item 5) is a storage choice and leaves the measure unchanged.
- **Knot measure: uniform over each variable's distinct observed values.** Francom et al. 2018 build the prior constant from "the number of unique marginal values", and R draws from them ([`unique.ind`](D:/aghado01/packages/BASS/R/bass.R:218), `vars.len`). The functional variable uses its n_r distinct output locations.
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
- **Independent checks:**
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
    - The draft follows R, and the paper's uniform is the named alternative.
  - *Support.* The combined term needs at least b nonzero points: "each resulting basis function must have at least 20 non-zero points".
  - *Change move.* Redraws one categorical factor's subset from its prior, which is symmetric (R's [`genBasisChangeCat`](D:/aghado01/packages/BASS/R/rjmcmc_all.R:77)).
  - *Levels.* The training factor's levels. A factor with one level carries no information and is refused at setup.
  - *Complement pairs.* Under D03's Liang convention, a categorical main effect on D and one on its complement have the same centred column up to sign, so they define the same model, and the measure counts each such model twice. This is a tilt to report, not a defect.
- **The functional variable** uses its own interaction cap (R's `maxInt.func`), as the hub's scope section requires, and its knots are the n_r distinct output locations.

### D02 (draft): count prior at the capacity limit

- **Prior:** the derivations review's global truncation (§4) with D01's term measure: π(λ, M, terms) ∝ Gamma(λ; h₁, rate h₂) · Pois(M; λ) · M! · Πₘ π(termₘ) · 1_F, where F is the cap, the support constraints and numerical admission (D06).
- **λ is integrated out**, as in KHAOS. The normalizer is one global constant, so λ's integral stays closed-form: p(M) ∝ Γ(h₁ + M)/(M! (h₂ + 1)^M), and a birth's count factor is (M + h₁)/((M + 1)(h₂ + 1)), restricted to F. λ leaves the chain and the sweep loses a block.
  - For each retained draw, λ | M ~ Gamma(h₁ + M, rate h₂ + 1), which is exact under global truncation. R and GBMARS update λ in this form; GBMARS calls it approximate because its stated prior is the conditionally normalized one.
- **This changes the model the papers state.** JSS and Francom et al. 2018 write the Poisson truncated and renormalized for each λ. Under that prior the Gamma update needs the cap's normalizer, a regularized incomplete gamma function: tractable for the cap alone, intractable once support and admission join F. The draft takes the global reading because it is what R and GBMARS sample, and one mechanism then covers the cap, support and admission.
- **The tilt:** λ's marginal is its Gamma prior tilted by P_λ(F). At R's defaults (h₁ = h₂ = 10, `maxBasis` = min(1000, n), [bass.R:267](D:/aghado01/packages/BASS/R/bass.R:267)) the cap rarely binds. Report the posterior mass at the cap and the support and admission rejection rates.
- **Defaults:** h₁ = h₂ = 10, from R and Nott, Kuk and Duc's Gamma(10, 10). Francom et al. 2018 name λ's prior as the main control on overfitting, and the docs should say so.
- **Independent check:** under a constant likelihood on D01's fixture, the count law against the negative binomial restricted to F by enumeration, and the retained λ draws against the Gamma mixture.

### D03 (draft): coefficient-prior modes and the τ hyperprior

- **g mode:** a | σ², τ, B ~ N(0, σ² (BᵀB)⁻¹/τ), with τ the precision multiplier as in R, JSS and Francom et al., and p(σ²) ∝ 1/σ² by default (g₁ = g₂ = 0). R's selectable modes are `g` and `jeffreys`; it also ships a ridge sampler (`rjmcmc_des_ridge.R`) that `bass()` never reaches.
- **τ hyperprior default: Zellner–Siow,** τ ~ Gamma(1/2, rate n/2), equivalently g = 1/τ ~ IG(1/2, n/2) (Liang et al. 2008), with E[τ] = 1/n.

  | Source | Prior as written | Stated intent | E[τ] as written |
  |---|---|---|---|
  | Francom et al. 2018 | Gamma(1/N, 1) | "centered over the unit information prior" | 1/N |
  | R documentation ([b.tau](D:/aghado01/packages/BASS/R/bass.R:34)) | rate "one over the number of observations", shape 1/2 | "centers … on the unit information prior" | n/2 |
  | JSS 2020 | Gamma(1/2, 2/n), "(shape and rate)" | "the Zellner–Siow Cauchy prior" | n/4 |
  | GBMARS | g ~ IG(1/2, 2/n) on the variance multiplier, under its own conjugate update | "leads to the Zellner–Siow Cauchy prior" | n/4 |
  | R code ([default](D:/aghado01/packages/BASS/R/bass.R:287)) | Gamma(1/2, rate 2/n) | none | n/4 |

  Every stated intent centres τ on unit information, and two sources name Zellner–Siow, which does. The written 2/n gives Zellner–Siow only when read as a scale. The code's rate contradicts every stated intent, so the derivations review's rate/scale reading now rests on documented intent rather than inference.
  - Named options: `BassLiteral`, Gamma(1/2, rate 2/n), for parity tests against R; and `Francom2018`, Gamma(1/N, 1). Each states its rate explicitly, since C#'s Gamma takes a scale (derivations review §6).
- **n counts observations.** For direct functional responses n = n_x · n_r, as in Francom et al.'s N and R's [`data$n`](D:/aghado01/packages/BASS/R/bass.R:248). Each fit to a supplied-basis or PCA component has n = n_x coefficient responses.
- **Intercept: the Liang et al. convention as the default, R's as a parity option.** The intercept gets a flat prior, the remaining columns are centred and get the g-prior, and k = M in the τ/(1+τ) factor.
  - *Why.* R keeps the intercept inside the g-prior (k = M + 1). The prior then shrinks the fitted level towards zero, and τ's conditional reads the uncentred fit, whose intercept share n·ȳ² grows with y's level, so τ's posterior and model-size selection depend on where y sits. Under the Liang convention τ's posterior is invariant to shifting y. (Bayes factors across models are defined under both conventions, since each puts its improper priors only on parameters every model shares; invariance is the difference.)
  - *How much it matters.* A y far from zero relative to its noise pushes τ small under R's convention, and a small τ makes every term's Occam factor, about (τ/(1+τ))^{1/2}, harsher, so fits come out sparser than the data warrant. R's own advice, a commented-out stop in its σ² update, is to centre y. The Liang convention is the Bayesian version of centring: it integrates the level instead of subtracting an estimate. For supplied-basis and PCA fits the coefficient responses already have zero mean across cases, so there the two conventions nearly coincide; the choice bites for scalar and direct-functional fits.
  - *Cost.* The centred Gram is the uncentred one minus a rank-one term, n·b̄b̄ᵀ, so the bordered updates and the functional Hadamard form carry over. In the proper mode the flat intercept needs a proper prior too, which the same standardization that settles σ²'s units supplies (next bullet); R's convention needs only σ²'s.
  - *It changes the model relative to R,* so it is the item in D03 that most needs the user's judgment.
- **Proper-prior mode, for D05.** The β = 0 rung of the likelihood-tempered path and SMC initialization need proper priors on σ² and, under the Liang convention, on the intercept. A proper IG(g₁, g₂) with g₂ > 0 acts in y's raw units (OQ19's problem), so this mode either standardizes y and records the transform or uses a scale-free proper prior. D05 chooses.
- **Other coefficient priors**, recorded so the contract leaves room for them:
  - *A ridge with a sampled τ* (GBMARS), an option.
  - *A ridge with a fixed scale,* an option only with the scale stated as part of the model: every added term then carries a factor in that scale, so the penalty on model size is set by it (Nott, Kuk and Duc's diffuse N(0, σ²/ε)). KHAOS found its fixed ridge best on noise-free simulators and its g-prior better on most noisy ones, which supports keeping both selectable. Defaults per regime come from a pre-registered comparison of our own.
  - *KHAOS's modified g-prior,* a candidate extension ([extensions brief](extensions-brief.md#khaoss-modified-g-prior)): an n-component g-prior that shrinks each term by its complexity, with per-term factors from its interaction order and degree and a global g₀² under an inverse-gamma prior. Its displayed conditional for g₀² has no data term, so derive it rather than copy it.
  - *Under robust or generalized-hyperbolic noise* (outside the completion contract): the unweighted g-prior as the default, which reduces to Gaussian BASS at unit weights and keeps the latent scales' conditionals in closed form, at one extra factorization per move as in BARS's Q2; the ridge with a sampled τ as the GBMARS-parity option.
- **Independent checks:** τ's marginal for a fixed structure against one-dimensional quadrature; shifting y leaves τ's posterior unchanged under the Liang convention and moves it under R's (a fixture); `BassLiteral` against R on deterministic statistics.
- **R's Jeffreys mode has no single target** (derived from R's source, 2026-10-25).
  - *What R does.* `beta.prior = "jeffreys"` holds τ at 0 ([start](D:/aghado01/packages/BASS/R/bass.R:334), [update](D:/aghado01/packages/BASS/R/mcmc_update.R:100)), so the coefficient prior is flat, and a is drawn from N(â, σ²G⁻¹).
  - *Its structure moves* add, for each term, (2π)^{1/2} · σ^{−1} · (|G′|/|G|)^{1/4} ([birth](D:/aghado01/packages/BASS/R/rjmcmc_des.R:39), [death](D:/aghado01/packages/BASS/R/rjmcmc_des.R:80); R's Σ log diag R is ½ log|G|). Its change move adds (|G′|/|G|)^{1/4}.
  - *Its σ² draw,* with a integrated out, is IG(g₁ + (N − 1)/2, g₂ + (yᵀy − Q)/2) ([σ² shape](D:/aghado01/packages/BASS/R/mcmc_update.R:74)).
  - *What coherent priors would give.* A flat prior a ∝ c^k gives c · (2πσ²)^{1/2} · (|G′|/|G|)^{−1/2} per term and a σ² shape of g₁ + (N − k)/2. The Jeffreys prior p(a | σ²) ∝ |G|^{1/2} σ^{−k} gives (2π)^{1/2} per term and no |G| factor. R's factors match neither, and its moves and its σ² draw imply different coefficient priors (flat times σ^{−2k}|G|^{3/4} for the moves, σ^{−(k−1)} for the σ² draw).
  - *So it is not ported as is.* It is an unqualified transition in the sense of the hub's scope section, and becomes a defect fixture.
  - *The root cause is general.* Any improper coefficient prior leaves the posterior over structures undefined, because its normalizing constant changes with the dimension (Bartlett's paradox), so an implementation has to pick an arbitrary constant per term.
- **The second coefficient mode, for the completion contract: fixed g.** The g-prior with τ held at a declared value, unit information τ = 1/N by default, which is BARS's convention. It keeps the Jeffreys mode's practical intent: coefficients are essentially unshrunk (by the factor N/(N+1)), and structure selection is set by the unit-information Occam factor, about ½·log(N + 1) per term, the BIC penalty. It is coherent: §1's collapsed marginal with no τ block.
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
  - *Consequence.* Under the improper σ² default, any admitted structure with k < N that fits y exactly makes the whole posterior improper. Deterministic simulators with piecewise-linear structure can produce one: two opposite hinges at the same knot span x exactly.
  - *Proposal.* In the improper mode, F excludes structures whose residual is zero to working precision (‖r‖ ≤ c·ε·‖y‖, from the residual itself; the oracle uses c = 10³), and reports the refusal. Deterministic simulators are steered to the proper mode. A zero response (R's convention) or a constant one (Liang's) is refused, as BARS refuses an identically zero response (R13).
  - *The ridge with a sampled τ* has the same structure, and is checked by the same quadrature when it is built.
- **Evidence: where it is reported** (proposal, 2026-10-25).
  - *Only from proper-mode runs:* the SMC normalizer (stage 8) or path sampling over a tempering ladder (D05). It is an optional result field, carried with its Monte Carlo error and with the priors and transforms it depends on, since the marginal likelihood is sensitive to prior scales in a way the posterior is not.
  - *Default runs report none.* Their improper priors leave the absolute marginal undefined. Within one fit the posterior over structures already carries the Bayes factors between its models, which is what model-size readouts report.
  - *Comparison across configurations and engines* defaults to the predictive checks built for OQ17 (WAIC and PSIS-LOO), which need no proper prior and no shared prior scales.
  - *Cross-engine evidence* is an extension ([extensions brief](extensions-brief.md#cross-engine-evidence)), since it needs a proper mode in every engine.

### D04 (draft): collapsed and conditional update order

- **The rule.** Integrate the coefficients a, and σ² where it is conjugate, out of every block up to their own draw, and draw them immediately before the first block that conditions on them. Every block then samples a conditional of the same joint posterior, and no block reads a statistic from a stale draw. R's τ update is the pattern this rules out: it reads ‖Rβ‖² from a coefficient draw made before σ² was redrawn (the inventory's partially collapsed ordering, derivations review §2).
- **Gaussian BASS.** The chain state is (structure, τ).
  1. Birth, death and change on p(y | structure, τ), with a and σ² integrated out (derivations review §1.2).
  2. τ by slice sampling on log τ, with its Jacobian, each evaluation O(1) from the cached (‖r‖², Q, k, N) (§2).

  No block conditions on a, σ² or λ, so each retained draw gets them by composition: σ² from its collapsed conditional, a given σ², and λ | M (D02). These are exact draws given the chain state, so no argument about scan order is needed.
  - *What changed from §2:* its third block, λ from the conjugate conditional, leaves the chain, because D02's draft integrates λ out. It returns as a draw at retention.
  - Under D03's intercept convention, k = M and the statistics are centred.
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

### D05 (draft): tempering and SMC target paths

Here b ∈ (0, 1] is the inverse temperature, as in `Temperature(b, path)`.
- **Default path: likelihood-tempered** (derivations review §3; also BARS's default under Q5). π_b ∝ prior(structure, τ) · p_b(y | structure, τ), where p_b is §3's closed-form tempered marginal with a and σ² integrated out. The count prior, the term measure and τ's prior are untempered at every rung.
- **Why this path:**
  - *Swaps are cheap.* The prior factors are identical across rungs and cancel, so a swap needs only ℓ_b = log p_b, O(1) from cached statistics (§3).
  - *Hot replicas stay meaningful.* They approach the actual prior over structures, where a powered posterior flattens it.
  - *No reference-measure dependence.* A powered density on a space of varying dimension depends on the reference measure of its coordinates. With D01's discrete knot measure the structures are countable, so they carry counting measure on either path, but τ is continuous, and powering its density depends on whether τ or log τ is the coordinate. The likelihood path's b = 0 member is the prior under any parameterization.
  - *An exact reference rung (OQ2).* At b = 0 the target is the prior, which can be drawn exactly: M and the terms from the unconstrained measure, accepted when they lie in F (rejection from the globally truncated measure is exact), and τ from its prior. This needs D03's proper mode.
- **The powered path stays selectable** (Q5), with τ powered in log τ coordinates and that choice recorded. Its b → 0 limit is its reference measure, which is not normalizable here, so it has no exact reference rung.
- **Which priors each use needs.** Parallel tempering with b_min > 0 runs under the default improper priors on σ² and the intercept. Each rung's collapsed marginal is finite for b > 0, and the improper prior's arbitrary constant appears with both signs in every within-rung ratio and every swap, so it cancels. The b = 0 rung, SMC started from the prior and path-sampling evidence need D03's proper mode, whose units question is OQ19's.
- **Moves at rung b** use ℓ_b(new) − ℓ_b(old), plus the untempered prior and proposal ratios. τ is slice-sampled on log τ from p(τ)·(τ/(τ+b))^{k/2}·(g₂ + b·S_b/2)^{−(g₁ + Nb/2)} (§3).
- **R's tempering is a defect fixture, not an option.**
  - R multiplies the prior and proposal factors of birth and death by the inverse temperature ([birth](D:/aghado01/packages/BASS/R/rjmcmc_des.R:34)).
  - Its swap density ([`lp`](D:/aghado01/packages/BASS/R/mcmc_funcs.R:64)) leaves out the term-prior blocks.
  
  Together these have no single target. Fixture: swap detailed balance on an enumerable structure space.
- **Schedule and ladder.** The non-reversible DEO schedule (step 3.3) and the adaptive ladder are reused; the ladder's adaptation follows D10's policy and freezes after burn-in.
- **Default: untempered.** R defaults to a single temperature (`temp.ladder` = 1). Francom et al. 2018 (§3.3) describe BMARS posteriors as often highly multimodal, but BARS found that tempering helped only where chains were stuck (OQ4). So BASS starts untempered, with tempering a declared option, until its own pre-registered measurement says otherwise.
- **SMC (stage 8)** uses the same sequence, b₀ = 0 < … < b_S = 1, started from the prior in the proper mode.
  - Each particle's incremental weight is ℓ_{b′} − ℓ_b, which costs O(1).
  - Mutations are D04's sweep at the stage's b.
  - Stages are placed adaptively by conditional ESS (Zhou, Johansen and Aston 2016).
  - Resampling uses independent offspring streams.
  
  The normalizer estimate is the evidence only in the proper mode, starting from a normalized prior.
- **Evidence:** path sampling over the same rungs or stages, with ∂_b ℓ_b in closed form (§3).
- **Independent checks:**
  - the endpoints: b = 1 gives D04's target, and b = 0 gives the prior;
  - swap detailed balance on an enumerable structure space;
  - the SMC normalizer against enumeration with proper priors;
  - exact reference draws against D02's count law;
  - the R defect fixtures above.

### D06 (draft): numerical model admission

- **Two gates, as Q6 decided for BARS** (step 3.7):
  - *Numerical:* the derivations review's angle rule (§5), in an order-free form. A structure is admitted only when every column's angle to the span of all the others satisfies s_j = 1/(G_jj (G⁻¹)_jj) ≥ s_min. s_j is the reciprocal of column j's uncentred variance inflation factor. It is scale-free, and the bordered step's s is its value for the new column.
  - *Support:* the nonzero-count rule b_m ≥ b of D01, counted separately for the design and functional parts (`npart`, `npart.func`). It is a prior support restriction inside F, not a numerical test.

  Exact duplicates and exact collinearity have s = 0 and fail the numerical gate, which D01's rule on repeated terms relies on.
- **Why order-free** (checked 2026-10-25).
  - *The counterexample.* The new column's angle depends on insertion order. With b₁ = e₁, b₂ = e₁ + 10⁻³e₂ and b₃ = e₂ + 10⁻³e₃, the leave-one-out values are (10⁻¹², 10⁻¹², 10⁻⁶). Added last, b₃ shows 10⁻⁶, while b₁ added last shows 10⁻¹². Under a newest-column rule at 10⁻¹⁰, the same set would be admitted or refused depending on its history. F would then not be a set, and the target restricted to it would be undefined.
  - *The gate is downward closed.* Removing a column only widens the others' angles, so a death never fails the gate. A birth or a change needs the whole vector.
  - *Updates in O(k²).* On a birth, with w = G⁻¹u, each old entry of diag(G⁻¹) gains w_j²/δ and the new entry is 1/δ (the bordered inverse). On a death of column m, each entry loses (G⁻¹)_{jm}²/(G⁻¹)_{mm}.
  - *Relation to BARS.* The bordered s is `FactorStatus.MinRelativePivot` (d_j/A_jj) for the new column. Each leave-one-out s_j is at most column j's pivot in any order, so the order-free gate is at least as strict. BARS's columns come in a canonical order (sorted knots), so its `MinRelativePivot` is already a property of the set; BASS's terms are unordered, which is why it needs the leave-one-out form.
- **s_min = 10⁻¹⁰,** the same threshold as BARS's `MinRelativePivot` gate (step 3.7).
  - *Numerics.* A computed s_j carries a relative error of order ε·κ(E), where E is the equilibrated Gram. For designs that pass the gate this stays small at the term counts BASS reaches. `EquilibratedCondition()` reports κ₁(E), and a decision within a factor of two of s_min is confirmed by a strict rebuild.
  - *Statistics.* 1/s_min = 10¹⁰ bounds every coefficient's variance inflation. It is a numerical guard, not a modelling restriction; the marginal's Occam factor already discourages redundant terms.
  - *Option:* a declared cap on variance inflation (`MaxVarianceInflation`), a statistical support restriction for stable extrapolation. Off by default.
- **R's rule: a defect fixture and a parity option.** R's `getQf` refuses a design when max_{j≥2} r_jj / min_j r_jj > 10³ ([`getQf`](D:/aghado01/packages/BASS/R/mcmc_funcs.R:118), commented "a hack, otherwise we get huge variance inflation in beta"). That is a pivot-ratio gate at 10⁻⁶ on the squared pivots.
  - It depends on column norms. Like BARS's `PivotRatio` (R15; step 3.7 counted 1,590 refusals against 76 under the scale-free gate), it refuses well-conditioned local terms whose columns are small, and does so more often as n grows.
  - Fixture: a term orthogonal to the rest, with b nonzero points, refused by R's rule at large n and admitted by the angle gate. At n = 40,000, a 20-point local term with values 0.01 gives R's ratio 1291, so R refuses it, while the smallest leave-one-out angle in the design is 0.25.
- **Reference fixture:** [`bass_admission_oracle.R`](../../../r/oracles/bass_admission_oracle.R) emits the leave-one-out angles, the pivots in the given order and R's `getQf` decision. The order counterexample and the large-n fixture are reproduced by [`bass_stage1_checks.R`](../../../r/scripts/bass_stage1_checks.R).
  - `BassLiteral` admission reproduces R's rule for posterior comparisons against R, and is never a default.
- **Count against column mass.** For BARS, Q6 replaced the count with column mass Σᵢ Bⱼ(xᵢ). The draft keeps the count as BASS's default: every BASS source uses it, and Francom et al. 2018 justify it by the edge instability of hinge terms with few nonzero values. Column mass, Σᵢ B_m(xᵢ) with hinges normalized to a maximum of one, is an option, and a fixture compares the predictive variance of edge terms under the two.
- **A strict rebuild before refusal.** When an incremental update fails the numerical gate, a strict rebuild decides before the candidate is refused (the [numerics segment's](port-numerics-brief.md#dense-factor-maintenance) dense factor maintenance). Refusals are reported by reason: numerical, support, cap.
- **A ridge is a prior, never a numerical repair.** Under the g-prior a rank-deficient design lies outside F; a ridge enters only as one of D03's declared coefficient priors.
- **Shared state.** The residual r and δ serve the factor border, the numerical gate and informed proposals (§5). r is updated on each accepted birth, recomputed after deaths and on a cadence, and S is computed from r rather than as yᵀy − Q/(1+τ), whose cancellation bites when τ is small (§5).
- **Independent checks:** §5's, randomized insert and delete sequences against fresh factorizations, residual drift against recomputation and s against the SVD's principal angle; plus one constructed case for each refusal reason.

## Regime note: deterministic simulators (proposal, 2026-10-26)

Three findings from the lineage and the propriety check point the same way.
- **KHAOS:** its fixed ridge did best on noise-free simulators, its g-prior on most noisy ones.
- **GBMARS's review of BMARS:** a_λ and b_λ smaller than the usual 10 suit deterministic surrogates, which need more terms to interpolate.
- **D03's propriety result:** under the improper σ² prior, a structure with k < N that fits y exactly makes the whole posterior improper, and deterministic data are where exact fits happen.

**Proposal:** a named profile, `DeterministicSimulator`, which is not a default:
- D03's proper σ² mode, with y standardized and the transform recorded;
- a wider λ prior (h₁ and h₂ well below 10, their values recorded);
- the ridge with a sampled τ as the alternative to the g-prior.

Which coefficient prior the profile uses is left to the pre-registered comparison, not set from KHAOS's study. The improper mode's exact-fit refusal stays, and points users to this profile.

## Room for extensions

Prep work in this segment, so that the [extensions](extensions-brief.md) plug in without reworking the port:
- **An observation-model seam with an offset.** BASS's observation model hands the Gaussian core a working response and weights, (y − offset, W), not weights alone. Gaussian BASS uses a zero offset and unit weights; generalized-hyperbolic noise (GBMARS) then needs only its latent-scale block.
- **A latent slot in the sweep.** D04's robust order (structure, τ, a, latent scales, then σ² and the skew) is the sweep's general form; Gaussian BASS leaves the latent slots empty.
- **Term types behind one design-kernel contract** (derivations review §7, item 6): hinge products, categorical subsets and the functional variable now, and later tensor-product B-splines (ABMRS) and directional hinges (learned projections). Each term type supplies its column evaluation, its support count, its moments and its derivatives.
- **Proposal kernels behind one contract,** each with exact forward and reverse probabilities. The Wallenius variable proposal is the first member; KHAOS's Bernoulli-inclusion proposal, the variable swap and a history-based proposal (D10's policy) join it later.
- **Coefficient priors as a declared set:** g with a sampled τ, fixed g, `BassLiteral` and a ridge with a sampled τ now, so KHAOS's modified g-prior is one more member.
- **Pointwise densities per observation,** with per-observation latents integrated under OQ17's scope rule, so the generalized-hyperbolic density is one more case.
