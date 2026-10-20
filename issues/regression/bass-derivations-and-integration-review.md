# BASS derivations and integration review

Date: 2026-10-19. Status: review. The resolutions below are proposed entries for the brief's decision register; each still needs its independent check before adoption. Source-code implementation is not authorized by this document.

This note cross-examines the [full port brief](bass-full-port-brief.md) and the [public inventory and backend review](bass-public-inventory-and-backend-review.md). It works the derivations those documents leave open for D02, D04, D05 and D06, records a D03 source finding on the τ hyperprior, shows that the existing BARS target belongs to the same collapsed g-prior family as BASS, and ranks the integration opportunities that follow. The [adaptation prerequisites review](bass-adaptation-prerequisites-review.md) was not re-examined here.

## Notation

Responses `y` (length `N`), design `B` (`N × k`, intercept included, full column rank), `G = BᵀB`, `h = Bᵀy`, `Q = hᵀG⁻¹h`, projector `P = BG⁻¹Bᵀ`, OLS residual `r = y − Py`, so `‖r‖² = yᵀy − Q`. Weighted variants use `BᵀWB`, `BᵀWy` and `yᵀWy`.

Priors, in R's parameter names:

\[
a\mid\sigma^2,\tau\sim N\!\left(0,\ \sigma^2G^{-1}/\tau\right),\quad
\sigma^2\sim \mathrm{IG}(g_1,g_2),\quad
\tau\sim\mathrm{Gamma}(a_\tau,\ \text{rate } b_\tau),\quad
M\mid\lambda\sim\mathrm{Pois}(\lambda),\quad
\lambda\sim\mathrm{Gamma}(h_1,\ \text{rate } h_2).
\]

R's defaults are `g1 = g2 = 0` (the improper `1/σ²` limit), `a.tau = 0.5`, `b.tau = 2/n` and `h1 = h2 = 10` ([signature](D:/aghado01/packages/BASS/R/bass.R:54), [τ default](D:/aghado01/packages/BASS/R/bass.R:285)).

## 1. One collapsed g-prior family

### 1.1 Coefficient-marginal target

Integrating out `a` gives the inventory's target. Written with nonnegative terms:

\[
\log p(y\mid B,\sigma^2,\tau)=-\frac N2\log(2\pi\sigma^2)+\frac k2\log\frac{\tau}{1+\tau}-\frac{S}{2\sigma^2},
\qquad S=\lVert r\rVert^2+\frac{\tau}{1+\tau}\,Q .
\]

`S` equals `yᵀy − Q/(1+τ)` algebraically; §5 explains why the form above is the one to compute.

The target depends on `B` only through `P` and `k`, so it is invariant to rescaling or any invertible recombination of the active columns:

- Hinge normalization constants do not enter the g-prior target. They matter for coefficient interpretation, sensitivity moments and numerics, and for flat or Jeffreys coefficient priors, which bring `|G|` into the marginal. That is one more reason the `jeffreys` mode needs its own derivation (D03).
- Near-collinearity between terms is a numerical and extrapolation concern, not a defect of the target. §5 turns it into the admission rule.

### 1.2 Variance-collapsed target

Integrating `σ² ~ IG(g1, g2)` as well:

\[
p(y\mid B,\tau)=(2\pi)^{-N/2}\left(\frac{\tau}{1+\tau}\right)^{k/2}
\frac{g_2^{\,g_1}\,\Gamma(g_1+N/2)}{\Gamma(g_1)\,(g_2+S/2)^{g_1+N/2}} .
\]

In the improper limit `g1 = g2 = 0` this is proportional to `(τ/(1+τ))^{k/2} S^{−N/2}`, with a constant shared by every model. The exact conditionals for composition sampling are

\[
\sigma^2\mid y,B,\tau\sim\mathrm{IG}\!\left(g_1+\tfrac N2,\ g_2+\tfrac S2\right),\qquad
a\mid y,B,\sigma^2,\tau\sim N\!\left(\frac{G^{-1}h}{1+\tau},\ \frac{\sigma^2G^{-1}}{1+\tau}\right).
\]

R already samples σ² from this conditional, with the comment "with beta marginalized out (maybe better stability)"; the earlier uncollapsed update is still commented out above it ([variance update](D:/aghado01/packages/BASS/R/mcmc_update.R:66)).

### 1.3 BARS is a member of the family

[WeightedNormalModel](../../src/maths/regression/spline/bars/WeightedNormalModel.cs:9) returns `p(y|Z) ∝ (m+1)^(−ν/2) · a^(−m/2)` with `a = yᵀWy − (m/(m+1))·Q`, for `m` observations and `ν` basis columns. Setting `τ = 1/N` in §1.2, with the improper variance prior, gives `τ/(1+τ) = 1/(N+1)` and `S = yᵀy − (N/(N+1))Q`: the BARS target exactly. The unit-information prior of DiMatteo, Genovese and Kass is the g-prior with `g = N`.

| Setting | BARS today | R BASS | Proposed BASS |
|---|---|---|---|
| τ | Fixed at `1/N` | Sampled under a Gamma hyperprior | Sampled under a declared hyperprior (§6) |
| σ² | Integrated out (improper prior) | Sampled from the collapsed conditional | Integrated out; drawn by composition for retained samples (§2) |
| Weights `W` | Scale-mixture robustness | None | Optional, from the BARS resamplers (§7) |
| Factor | Banded LDLᵀ | Dense, refactored per proposal | Dense incremental; Hadamard Gram for direct functional responses (§5) |

The brief is right that `WeightedNormalModel` cannot be reused unchanged: its τ is fixed and its factor is banded. What both engines share sits one level down. A collapsed g-prior block holding `(‖r‖², Q, k, N)` and a factor backend serves both through the settings above, and every closed form in §§2–5 is written once. BAPS uses a difference-penalty precision rather than one proportional to `G`. It adds a `log|Λ+G| − log|Λ|` term and can join the block when it migrates.

### 1.4 Why the identity was hard to see, and how to sort engine differences

The two implementations hide the shared target in different ways. BARS writes it in closed form, with τ folded into the constant `m+1` and σ² already integrated out. R's BASS never writes the full target. It computes per-move changes inside acceptance ratios (`½log τ − ½log(1+τ)` per added term, a change in `Q` divided by `1+τ`), spread across six data-mode files, and it samples σ². The identity appears once BASS's target is written whole, σ² is integrated out, and BARS's constant is read as a parameter: `m+1 = (1+τ)/τ` at `τ = 1/m`. Both engines descend from Denison, Mallick and Smith's 1998 pair of papers on reversible-jump curve fitting and Bayesian MARS.

Because the target depends only on the spanned space (§1.1), every difference between two engines falls into one of three bins:

| Bin | Meaning | BARS vs BASS |
|---|---|---|
| Representation | Same space, same target | B-splines or hinges for a given space; hinge normalization; banded LDLᵀ, dense Cholesky or QR; code organization |
| Sampler | Same posterior, different chain | σ² integrated out or sampled; move design; any valid tempering scheme |
| Model | Different posterior | Reachable spaces; the structure prior; τ fixed or sampled |

Code is shared freely across the first two bins; only the model bin holds genuine engine differences. The banded-versus-dense contrast, for instance, comes from basis support (B-splines are local, hinges one-sided) and is invisible to the statistics.

Sorted by the model bin, the spline engines are points in one space:

| Engine | Knots | Reachable space | Coefficient prior |
|---|---|---|---|
| BARS | Adaptive, continuous | The full one-dimensional spline space for its knots | g-prior, τ = 1/N |
| BAPS | Fixed, dense | The full tensor-product spline space | Difference penalty |
| BASS | Adaptive, at data values | Sparse subspaces of tensor-product spline spaces, spanned by hinge products | g-prior, τ sampled |

In one dimension with matching degree, each BASS model spans a subspace of the spline space BARS uses for the same knots.

What stays per-engine is the basis family, the moves and the structure prior: which subspaces an engine reaches, how it proposes them, and how its prior weights them. The ensemble shell and the collapsed observation model are shared. Transfers run both ways. BARS offers BASS its robustness weights, its ensemble shell, and its structure priors as named options. BASS offers BARS a sampled τ: the fixed unit-information value gives a BIC-like penalty of `(k/2)·log(N+1)`, while a mixture-of-g hyperprior such as Zellner–Siow adapts the shrinkage to the data (Liang et al. 2008).

## 2. D04 proposal: full collapse with composition sampling

The chain state is `(structure, τ, λ)`. One sweep:

1. Structural birth, death and change moves on `p(y|B,τ)` from §1.2. Neither `σ²` nor `a` is conditioned on.
2. τ by slice sampling on `log τ` (with its Jacobian). The density is `p(τ) · (τ/(1+τ))^{k/2} · (g2 + S(τ)/2)^{−(g1+N/2)}`, and each evaluation is O(1) from the cached `(‖r‖², Q, k, N)`.
3. λ from its exact conjugate conditional (§4).

For each retained sample, draw `σ²` and then `a` from the §1.2 conditionals. These are exact conditionals given the chain state. Because `σ²` and `a` are not part of the chain, no scan-order invariance argument is needed, which settles D04's partial-collapse question by construction. Collapsing typically improves mixing as well (Liu, Wong and Kong 1994).

Source context: as traced in the inventory, R draws coefficients, then σ² from the collapsed conditional, then updates τ using `qf2 = ‖Rβ‖²` from the earlier coefficient draw together with the new σ² ([qf2](D:/aghado01/packages/BASS/R/mcmc_update.R:67), [τ update](D:/aghado01/packages/BASS/R/mcmc_update.R:97)). That is the partially collapsed ordering the inventory flags.

Independent check: for a fixed structure, compare the τ chain's marginal to quadrature of the one-dimensional density; on a tiny enumerable model, compare composition draws with a long joint-Gibbs reference.

## 3. D05 proposal: the likelihood-tempered path

Temper the likelihood only: `π_β ∝ prior(structure, λ, τ, σ², a) · L(a, σ²)^β`. This is the power posterior of Friel and Pettitt (2008). It differs from [TemperedTarget](../../src/maths/samplers/rjmcmc/TemperedTarget.cs), which powers the whole posterior. Integrating `a` stays closed-form:

\[
\log p_\beta(y\mid B,\sigma^2,\tau)=-\frac{N\beta}2\log(2\pi\sigma^2)+\frac k2\log\frac{\tau}{\tau+\beta}-\frac{\beta S_\beta}{2\sigma^2},
\qquad S_\beta=\lVert r\rVert^2+\frac{\tau}{\tau+\beta}\,Q,
\]

\[
a\mid\cdot\sim N\!\left(\frac{\beta\,G^{-1}h}{\tau+\beta},\ \frac{\sigma^2G^{-1}}{\tau+\beta}\right),\qquad
\sigma^2\mid\cdot\sim\mathrm{IG}\!\left(g_1+\tfrac{N\beta}2,\ g_2+\tfrac{\beta S_\beta}2\right).
\]

With σ² integrated out as well:

\[
\log p_\beta(y\mid B,\tau)=-\frac{N\beta}2\log 2\pi+\frac k2\log\frac{\tau}{\tau+\beta}
+g_1\log g_2+\log\Gamma\!\left(g_1+\tfrac{N\beta}2\right)-\log\Gamma(g_1)-\left(g_1+\tfrac{N\beta}2\right)\log\!\left(g_2+\tfrac{\beta S_\beta}2\right).
\]

Properties:

- **Endpoints.** β = 1 recovers §1. β = 0 gives zero, which is the prior, and that endpoint is finite only with a proper variance prior (`g1, g2 > 0`). This is the brief's SMC initialization requirement in concrete form.
- **Prior blocks stay untempered.** The λ and τ conditionals are the same at every rung, and hot replicas approach the actual structure prior as β falls, where `TemperedTarget`'s whole-posterior power flattens it. R instead powers them ([rgammaTemper](D:/aghado01/packages/BASS/R/misc.R:13), used for [λ](D:/aghado01/packages/BASS/R/mcmc_update.R:40) and [τ](D:/aghado01/packages/BASS/R/mcmc_update.R:100)) while also tempering the proposal factors in birth and death. That combination has no single target.
- **Swaps.** The prior factors are identical across rungs and cancel. With `ℓ_β = log p_β(y|·)`, the swap is `log α = ℓ_{β_i}(x_j) + ℓ_{β_j}(x_i) − ℓ_{β_i}(x_i) − ℓ_{β_j}(x_j)`, each term O(1) from cached statistics. `ℓ_β` is nonlinear in β, so the swap at [ParallelTempering.cs:71](../../src/maths/samplers/rjmcmc/ParallelTempering.cs:71), `(β_i − β_{i+1})(L_{i+1} − L_i)`, applies to the powered-posterior path only.
- **Reference measure.** On a variable-dimension space, a powered density depends on the reference measure chosen for each model's coordinates. The likelihood-tempered path does not: its β = 0 member is the prior under any parameterization.
- **SMC.** A particle's incremental weight is `ℓ_{β′} − ℓ_β`, which is O(1). Only the mutations do data work.
- **Evidence.** `log p(y) = ∫₀¹ E_{π_β}[∂_β ℓ_β] dβ` (path sampling; Gelman and Meng 1998), and `∂_β ℓ_β` is closed-form. PT rungs and SMC stages both supply the integrand.

Naming: in the literature, "power posterior" means this likelihood-tempered path. Code names should say "likelihood-tempered" to stay distinct from `TemperedTarget`.

Independent check: the β endpoints; swap detailed balance on an enumerable structure space; the SMC normalizer against enumeration with proper priors.

## 4. D02 proposal: a globally truncated joint prior

Let `F` be the feasible set: the basis-count cap, the support constraints, and numerical admission (§5). Define

\[
\pi(\lambda,M,\text{terms})\propto\mathrm{Gamma}(\lambda;h_1,h_2)\,\mathrm{Pois}(M;\lambda)\,p(\text{terms}\mid M)\,\mathbf 1_F .
\]

The normalizer is one global constant, so `λ | M ~ Gamma(h1 + M, h2 + 1)` exactly. That is R's update ([λ update](D:/aghado01/packages/BASS/R/mcmc_update.R:38)), untempered under §3. The inventory's concern about the simple Gamma update is therefore a question of which prior is meant, and under this reading the update is exact.

The price is that λ's marginal prior becomes the stated Gamma tilted by the feasibility probability `P_λ(F)`. The conditionally normalized alternative is tractable for the cap alone, through the regularized incomplete gamma function. It is intractable for admission, because the per-M normalizer over admissible term sets depends on the training inputs. Global truncation handles the cap, support constraints and admission with one mechanism.

At R's defaults (`h1 = h2 = 10`, so λ has prior mean 1, and `maxBasis = 1000`) the cap rarely binds, and admission is the part that matters. Report the posterior mass at the cap and the admission rejection rate so the tilt stays visible.

Independent check: under a constant likelihood, the count law on a tiny model matches the tilted, truncated Poisson computed by enumeration.

**D01 fixture expectation.** For the inventory's boundary-selection issue: under a constant likelihood, and with the term-level proposal factors otherwise correct, R's fixed one-third move probabilities make `π_R(1)/π_R(0)` three times the correct ratio and `π_R(M_max)/π_R(M_max − 1)` two-thirds of it. The empty model is under-weighted by a factor of three, which matters most for null and weak-signal data.

## 5. D06 proposal: admission by angle, with the residual as shared state

For a candidate column `b`, the bordered step computes `u = Bᵀb`, `v = bᵀb`, `Lz = u` and `δ = v − zᵀz = ‖(I − P)b‖²`. Three quantities follow:

\[
s=\frac{\delta}{v}=\sin^2\theta,\qquad
\Delta Q=\frac{(b^\top r)^2}{\delta},\qquad
r'=r-\frac{b^\top r}{\delta}(I-P)b,\quad \lVert r'\rVert^2=\lVert r\rVert^2-\Delta Q,
\]

where θ is the angle between `b` and the current span.

- **Admission.** Admit a term when `s ≥ s_min`, with `s_min` declared as part of `F` (§4) so that the rule is a stated prior support restriction. It is scale-free, unlike an absolute `δ ≤ 1e-12`. Statistically, it bounds extrapolation variance: a column nearly inside the training span produces large opposing coefficients that behave badly at new inputs. Numerically, `δ` carries an absolute error of order `ε·v` times a modest growth factor, plus accumulated drift, so the decision is reliable when `s_min` is well above `k·ε`.
- **Shared state.** The same `(r, δ)` serves the factor border, the admission test and the informed-proposal score.
- **Cancellation.** Computing `S` as `yᵀy − Q/(1+τ)` loses about `log₁₀(yᵀy/S)` digits. Because `S ≥ τQ/(1+τ)`, the loss is at most `log₁₀((1+τ)/τ)`; it is large only when τ is small *and* the fit is close (`‖r‖²` small). The ratio is uncentered, so a response offset adds to it.
  - **BARS:** τ = 1/N caps the loss at `log₁₀(N+1)` digits, which is negligible, so BARS's current form is adequate. Its `a ≤ 0` guard rejects only an identically zero response (R13 in the prerequisites register; corrected 2026-10-19, after the first version of this note overstated the BARS case).
  - **BASS:** a sampled τ tracks the signal-to-noise ratio, roughly `k/(N·SNR)`, so it is small in the emulator regime, where fits are also close. That is where the loss becomes real. R shows the symptom: it floors a nonpositive `s2.rate` to machine epsilon, next to a commented-out stop asking users to center or rescale `y` ([floor](D:/aghado01/packages/BASS/R/mcmc_update.R:69)).

  For the shared block, keep `r` explicitly and compute `S` as in §1.1:
  - **Scalar responses:** update `r` in O(N) per accepted birth, and recompute it at O(Nk) after deaths and on a fixed cadence.
  - **Direct functional responses:** `r` is `n × m`; refresh it per accepted move or on a cadence, at O(nmk).
- **Factor reference.** The target depends only on the span, so an updatable QR, an orthonormal basis for the span, is a natural strict reference for scalar responses and avoids squaring the condition number. Direct functional responses need the Gram path because of their Khatri–Rao structure. Both can satisfy one factor contract.

Independent check: randomized insert/delete sequences against fresh factorizations; residual drift against recomputation; `s` against the SVD-computed principal angle.

## 6. D03 finding: R's default τ hyperprior

R sets `τ ~ Gamma(0.5, rate 2/n)` ([τ default](D:/aghado01/packages/BASS/R/bass.R:285)). Two places confirm the rate semantics: the conjugate update ([τ update](D:/aghado01/packages/BASS/R/mcmc_update.R:97)) and the prior term in the log posterior used for swaps ([log posterior](D:/aghado01/packages/BASS/R/mcmc_funcs.R:75)). The parameter documentation says the default is `1/n` and that it centers the coefficient prior on the unit-information prior ([b.tau](D:/aghado01/packages/BASS/R/bass.R:34)).

Under the rate reading, `E[τ] = n/4`, whereas unit information is `τ = 1/n`. Read as a scale, `2/n` gives `Gamma(1/2, rate n/2)`, which is exactly the Zellner–Siow prior (`g ~ IG(1/2, n/2)`; Liang et al. 2008). That is consistent with a rate/scale slip, though the source cannot prove intent.

The impact depends on signal strength. The τ conditional's rate is `b_τ + aᵀGa/(2σ²)`, and the data term is about `n·SNR/2`. The two priors agree when the signal-to-noise ratio is large and diverge for weak signals, which is where model-size selection is most sensitive.

Proposed D03 content:

- **Named hyperpriors.** Name each τ hyperprior and state its rate explicitly. Zellner–Siow is the default candidate; R's actual default remains available as a named option for parity tests. The C# Gamma takes a scale, so this is a concrete case of the inventory's rate/scale warning.
- **Intercept convention.** R keeps the intercept inside the g-prior (`k = M + 1`). The Liang et al. convention puts a flat prior on the intercept and the g-prior on the centered remainder, so improper priors fall only on parameters every model shares. That makes Bayes factors across models well-defined, which §7 item 8 needs. It changes the model relative to R and needs its own decision.

## 7. Integration opportunities

Ranked by cost, given the code that already exists.

| # | Opportunity | Builds on | Lands with |
|---|---|---|---|
| 1 | Collapsed g-prior block shared by BARS and BASS (§1.3); BARS gains a sampled τ through it (§1.4). | `WeightedNormalModel`, `BandCholesky` | Stage 2, after D04 and D06 |
| 2 | Likelihood-tempered PT for both engines, with even/odd non-reversible swaps and an adaptive ladder (Syed et al. 2022). Under that scheme the rate of full round trips across the ladder stays bounded as rungs are added; under reversible PT it decays. | `ParallelTempering`, `GeometricLadder` | Stage 3, after D05 |
| 3 | SMC with O(1) reweighting (§3). | Item 1 and a proper σ² prior | Stages 3 and 8 |
| 4 | Robust BASS from the BARS scale-mixture resamplers (`StudentTWeights`, `ContaminatedNormalWeights`). Each weight sweep rebuilds `G` and `h`, so run it at a slow cadence. The largest benefit is in the calibration measurement model. | BARS robustness seam | After stage 4 |
| 5 | Immutable, interned terms. A training-column cache keyed by term id serves PT swaps and a term's death and rebirth. Posterior prediction becomes one product: unique-term columns at the new inputs times a sparse term × draw coefficient matrix. Sensitivity moments are computed per unique term pair. | New posterior layer | Stages 4–6 |
| 6 | One shared target and move-accounting layer over generic struct design kernels for the six data modes. Static abstract interface members let the JIT specialize each mode. | — | Stage 4 |
| 7 | A moment service for piecewise polynomials rather than hinges only, so BAPS tensor and anisotropic splines get variance decompositions from the same code. | Distributions | Stage 6 |
| 8 | Cross-engine evidence (BARS, BASS, BGP) by path sampling. | Item 2, the §6 intercept convention, proper priors | After stage 8 |
| 9 | Direct-functional BASS `f(x, T)` for families of SPC temperature sweeps; BARS remains the engine for a single sweep. | Application layer | Outside the completion contract |

## 8. Document notes

- The inventory's seven-stage sequence predates the brief's stages 0–9; the inventory now points to the brief's stage table.
- **Stage-0 gate scope.** The brief justifies gating BASS integration on all of B01–B06 by confounding ([brief](bass-full-port-brief.md#existing-engine-baseline-before-integration)). B01, B02 and B06 sit on BASS's code path: shared diagnostics, mutable sampler state and LinAlg factors. B03 partly does, because generalizing `IRjTarget` migrates the BARS callers; if the shared block of §1.3 is adopted, B03's singular-design and conditional-mean items join that path as well. B04 (BAPS) and B05 (BGP) have no code path into BASS and confound only cross-engine comparisons. Either gate B04 and B05 on the first cross-engine comparison, or state the gate as a discipline rule (reviewed defects close before a new engine lands) instead of resting it on confounding.
- An empty `src/maths/regression/bass/` directory exists, while the brief places BASS at `regression/spline/bass/`.
- Brief line 274 states its modernization constraints as prohibitions. A positive restatement: dispatch SIMD by hardware with a scalar fallback; represent categories with a bitset that scales past 64 levels; score candidate terms by fitted gain so interaction-only signals stay visible; let dimension estimates inform proposals while the interaction cap remains a declared control.
- The brief's `WeightedNormalModel` reuse row can reference §1.3 once the shared block is adopted.

## References

- Denison, Mallick and Smith (1998). Automatic Bayesian curve fitting. *JRSS-B* 60.
- Denison, Mallick and Smith (1998). Bayesian MARS. *Statistics and Computing* 8.
- DiMatteo, Genovese and Kass (2001). Bayesian curve-fitting with free-knot splines. *Biometrika* 88.
- Friel and Pettitt (2008). Marginal likelihood estimation via power posteriors. *JRSS-B* 70.
- Gelman and Meng (1998). Simulating normalizing constants: from importance sampling to bridge sampling to path sampling. *Statistical Science* 13.
- Liang, Paulo, Molina, Clyde and Berger (2008). Mixtures of g priors for Bayesian variable selection. *JASA* 103.
- Liu, Wong and Kong (1994). Covariance structure of the Gibbs sampler with applications to the comparisons of estimators and augmentation schemes. *Biometrika* 81.
- Syed, Bouchard-Côté, Deligiannidis and Doucet (2022). Non-reversible parallel tempering: a scalable highly parallel MCMC scheme. *JRSS-B* 84. [arXiv:1905.02939](https://arxiv.org/abs/1905.02939).
