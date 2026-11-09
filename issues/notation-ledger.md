# Notation ledger

Splits in the module's language (symbols, identifiers, terms), logged as sessions find them. They wait for one unification pass across docs and code, run after the BASS regression work and its extensions are implemented, so every name is judged against the whole vocabulary at once. Entries record a split and where each reading lives; renames wait for that pass.

**Adding an entry.** Prepend it under `## Entries`, newest first, with:

- **Kind:** one concept under several names, or one name for several concepts.
- **Readings:** each meaning or spelling, with anchors. Code takes `file:line`; docs take a section or phrase, since their line numbers drift.
- **Found:** the date and the work that turned it up.
- **Note** (optional): what the pass needs to know, such as the sources' own symbol.

## Entries

### ω_H beside the local steps' ω; δ for the history birth's floor beside the factor border's δ; ĥ_u beside h₁, h₂ and ĉ (2026-11-08)

- **Kind:** one name for several concepts.
- **Readings:**
  - ω_H, the history birth's share of the births (`HistoryBirth.Share`; `docs/Maths/BayesianRegression.md`, "BASS: the birth mixture and the history birth (E4)"), beside ω, the local steps' window fraction (`LocalSteps.WindowFraction`, "BASS: the local-move options (F-M2)").
  - δ, the history weights' floor (`HistoryBirth.Floor`), which is Nott and Kohn's δ, beside δ, the factor border's residual norm in the design (the port-fitting brief's shared state; `BassDesign`).
  - ĥ_u, a variable's inclusion frequency over the chain's history (`BassInclusionHistory.Frequencies`), which is Nott and Kohn's ĉ, beside h₁ and h₂, the count prior's hyperparameters (`BassCountPrior`), and h in the priors for a moving response (the 2026-11-02 entry).
- **Found:** 2026-11-08, the proposal-kernels brief's slice 1.

### G for the reconstruction's map and for the output side; k for a fixed and a per-draw rank (2026-11-08)

- **Kind:** one name for several concepts.
- **Readings:**
  - G = diag(y_s)K, a basis fit's reconstruction map, in slice 3's text (`docs/Maths/ActiveSubspaces.md`, "Functional outputs and co-active pairs"; `BassActiveSubspaces.ReconstructionMap`, the XML on `OutputMetric` and `BassMoments.GradientCross`).
  - G = E[J W_x Jᵀ], the output side, in the Grassmann discussion, E8's brief (slice 5) and the page's slice-5 section (`OutputPosterior`); slice 5 writes the map out as diag(y_s)K.
  - k, the readout's fixed dimension (`ActiveSubspaceSpec.Dimension`), beside k_s, a draw's own rank by a rule, and k*, the rank summary's (`RankPosterior.SummaryRank`); the projection distance's k in the discussion is whichever subspace a draw is compared with.
- **Found:** 2026-11-08, E8's slice 5.

### Whitening's map L = Σ^{1/2} beside 2307.11241's Σ^{−1/2} and its L for mixture components; "Whitened" and "Standardized" (2026-11-08)

- **Kind:** one concept under several names, and one name for several concepts.
- **Readings:**
  - The module's L maps metric to box coordinates, x = μ + Lu, so whitening takes L = Σ_ρ^{1/2} and C_u = LᵀCL (`ActiveSubspaceAnalysis.Map`, `Whitening.SymmetricRoot`; `docs/Maths/ActiveSubspaces.md`, "Dependent inputs"). 2307.11241 §3.3.2 maps the other way, z = Σ^{−1/2}(x − μ) with A = Σ^{−1/2}, and writes C_x = (Σ^{−1/2})ᵀC_zΣ^{−1/2}.
  - L also counts mixture components in 2307.11241 (§3.3.3, §4.1.1's L uniform boxes) and in its supplement (L = 478 boxes).
  - `InputMetric.Whitened` scaled each input by its sd in slice 2, which is whitening only because every law then was a product. From slice 4 it is L = Σ_ρ^{1/2} under any law, the same values under a product law. `InputMetric.Standardized` names the diagonal form.
- **Found:** 2026-11-08, E8's slice 4.

### W for an output metric and for eigenvectors; the co-active matrices' indices; the output functionals (2026-11-08)

- **Kind:** one name for several concepts, and one concept under several names.
- **Readings:**
  - The output metric of C = E[JᵀW_yJ]: W_y in the active-subspace brief and `docs/Maths/ActiveSubspaces.md` ("Functional outputs and co-active pairs"), `OutputMetric` and `ActiveSubspaceSpec.OutputMetric` in `src/maths/regression/spline/bass/sensitivity/ActiveSubspaces.cs`, the local `w` in `BassActiveSubspaces.Of`, and Q in F-X6.
  - W for eigenvectors: W_s, a draw's top-k eigenvectors in decision 4's P̄ = Σ w_s W_s W_sᵀ (`DrawSubspace.Basis`), and 2311.18146's W_kℓ, V_kℓ's eigenvectors (eq. Vdecomp), here `CoActivePosterior.Directions`; the paper's columns w_i^(kℓ) sit beside the module's w_s for draw weights.
  - The co-active matrices: the paper's C_kℓ, V_kℓ and t_kℓ, indexed by model, against the module's C_fg, V and t_fg (`CoActivePair`, `CoActivePosterior` in `src/maths/regression/subspaces/CoActivePosterior.cs`); coactivity calls C its "Constantine matrix".
  - The output functionals: a_f and a_g in the docs and XML comments, `CoActiveSpec.Output` and `OtherOutput`, the locals `a`, `b` (two fits) and `a`, `c` (one fit) in `BassActiveSubspaces.CoActive`.
- **Found:** 2026-11-08, E8's slice 3.
- **Note:** b stays a draw's coefficients, so the output functionals are not a and b in writing; W_y is the brief's own symbol, taken because the page already had W_s.

### "Slope" and "derivative" for φ′, "gradient" for ∇f and for C, and the paper's h, γ, I₁–I₃ (2026-11-08)

- **Kind:** one concept under several names, and one name for several concepts.
- **Readings:**
  - φ′, a factor's derivative: `DesignFactor.Derivative` at a point, `FactorDerivative`, `BassBasis.CaseDerivative`, against its moments `DesignFactor.SlopeMean` and `SlopePair` (`src/maths/regression/spline/bass/DesignKernel.cs`) and the slope flags of `PiecewisePolynomial.HingeProduct`; `ShiftedLegendre` names it `Derivative`.
  - "Gradient": ∇f pointwise (`BassGradient`, `BassFit.Gradient`), against C = E[∇f∇fᵀ] (`BassMoments.GradientOuterProduct`) and a term pair's E[∂ᵢB_a∂ⱼB_b] (`BassMoments.GradientPair`, `TermPairGradientMoments`), `src/maths/regression/spline/bass/sensitivity/BassMoments.cs`.
  - Rumsey et al.'s h_{im}, γ_m and I₁, I₂, I₃ (arXiv 2307.11241, eq. Cij), against the module's φ, b and the factor-pair moments E[φ′ψ], E[φψ], E[φ′ψ′] (`docs/Maths/ActiveSubspaces.md`, "C on the exact route").
- **Found:** 2026-11-08, E8's slice 1.
- **Note:** "slope" was taken for the moments because `Derivative` already names the pointwise value on the same type.

### V, D and s in the random-intercept block beside the ridge's spectrum and σ²'s scaled prior (2026-11-08)

- **Kind:** one name for several concepts.
- **Readings:**
  - V: the block's covariance factor, σ²V = Cov(y), V = I + s²(ZH)(ZH)ᵀ (`RandomInterceptBlock`, `src/maths/regression/spline/bass/BassRandomIntercepts.cs`; `docs/Maths/BayesianRegression.md`, "BASS: random-effect intercepts (E21)"), against V, the eigenvectors of G̃ = VΛVᵀ and of L⁻¹AL⁻ᵀ in the ridge's and the pencil's spectra (the same page, "the coefficient prior's two axes, and the ridge").
  - D: the levels' counts ZᵀZ in the block, against D = diag(G_c)/N in the normalized ridge, and the prior precision τD of the pencil with its D_α, D_β.
  - s: the block's scale, the effects' sd over σ (`BassRandomIntercept.Scale`), against s² in `Sigma2Prior.ScaledToResponse`'s IG(c, c·s²), the response's variance, and s = ±1, a hinge's sign.
- **Found:** 2026-11-08, KHAOS slice 3 (E21).

### "Uniform bound" for the inflated rule and for the margin, and ε for CB-3's bounds beside Vovk's level (2026-11-08)

- **Kind:** one name for several concepts.
- **Readings:**
  - The inflated rule's uniform bound, the held ordering's largest per-position bound (`InflationBounds.Uniform`, `ConformalMargin.UniformErrors`, `src/maths/regression/conformal/ConformalBayes.cs`), against the margin's uniform bound, the largest bound carried over the orderings (`MarginBounds.Uniform`, `ConformalMargin.CarriedUniformErrors`). `docs/Maths/ConformalPrediction.md`, "CB-3: the scores' error and what it costs".
  - ε, CB-3's error bound, the Lean's εⱼ per unit and "the Lean's uniform ε" (`MarginBounds.Uniform`'s doc; the brief's contract line), against ε, the level in Vovk et al.'s conformal predictive distributions (the page's slice 4 notes, beside α_i).
- **Found:** 2026-11-08, the conformal layer's per-unit bounds.

### κ, the skew's rescaling factor, beside κ in the priors for a moving response (2026-11-08)

- **Kind:** one name for several concepts.
- **Readings:**
  - κ, the factor of the move along β·v, (β, v) → (β/κ, κ·v), drawn as t = ln κ (`BassSkewRescaleBlock`, `src/maths/regression/spline/bass/BassNoise.cs:713`; the move's section in `docs/Maths/VarianceMeanMixtures.md`).
  - κ, the coefficient prior's constant in σ²κ(τG)⁻¹ and a term's normalizing constant (the 2026-11-02 entry).
- **Found:** 2026-11-08, the move along β·v.

### g in KHAOS's modified g-prior beside BASS's g-form and τ (2026-11-08)

- **Kind:** one concept under two names, and one name for several concepts.
- **Readings:**
  - KHAOS's g₀², the modified g-prior's overall scale, is the engine's τ inverted, τ = 1/g₀² (`BassState.LogTau`; `ModifiedGForm`, `src/maths/regression/spline/bass/BassPriorPrecision.cs:14`), and `khaos`'s IG(a_g, b_g) on it is `TauPrior.Gamma(a, b)` on τ (`src/maths/regression/gprior/TauPrior.cs`). BASS's g-form writes the same scale as 1/τ in σ²G_c⁻¹/τ (`docs/Maths/BayesianRegression.md`, D03), so g₀² is the g-form's g at ζ = 0.
  - g_m, a term's complexity factor (1 + q(d + q − 2))^(−ζ/2) (`ModifiedGForm.Factor`), scales one coefficient, not the prior; BASS has no counterpart.
  - `khaos`'s g₁ and g₂ are its ridge's σ² hyperparameters (`adaptive_khaos_ridge`), and its `g2` field is the draws of g₀² (`fit$g2`), so `khaos` spells two quantities g2.
  - TauPrior's Gamma kinds: `UnitInformationGamma` is Gamma(a, rate a/τ₀) on the form's unit information, and `Gamma` is Gamma(a, rate b) on τ itself.
- **Found:** 2026-11-08, KHAOS slice 2 (the modified g-prior).
- **Note:** the paper writes g₀² and g_m; `khaos`'s code writes `g2` for g₀² and `g_vec` for the g_m.

### Q, the split form's threshold, beside Q(y, τ), the conformal predictive distribution (2026-11-08)

- **Kind:** one name for several concepts.
- **Readings:**
  - Q, the split form's negated threshold, the k-th largest calibration nonconformity, in CQR's set [q_lo − Q, q_hi + Q] and the residual's μ ± Qρ (`docs/Maths/ConformalPrediction.md`, "The split form (CB-2)"; `SplitConformal.Quantiles`, `SplitConformal.Residuals`).
  - Q(y, τ) = (B + τ(T + 1))/m, the conformal predictive distribution, with Q⁻ and Q⁺ its ends at τ = 0 and 1 (`ConformalDistribution.Lower`, `.Upper`, `.At`, `src/maths/regression/conformal/ConformalDistribution.cs`; `ConformalPrediction.md`, "The predictive distribution, built (slice 4)").
  - F̂(y | x), the split form's score, the training fit's predictive CDF, beside F̃, the full step's in-sample one (entry of 2026-11-07 below).
- **Found:** 2026-11-08, the conformal layer's slice 4 (the predictive distribution, built).
- **Note:** Vovk et al. write Q for the distribution and the CQR paper writes Q for the threshold, so each is its source's own; the page keeps both, each in its own section.

### u and z of the z parameterization, beside u = ln v and z in a moving response, and BARS's w beside GBASS's (2026-11-08)

- **Kind:** one name for several concepts.
- **Readings:**
  - u = σ = √w and z = σ·v, the noise's scale and its latent under the z parameterization (the user's derivation, E1, 2026-10-26; `BassInterweavingBlock`, `src/maths/regression/spline/bass/BassNoise.cs:931`); slice 7's section in `docs/Maths/VarianceMeanMixtures.md` writes σ and z.
  - u = ln v, the coordinate of the GBP's slice step (`GbpMixing.StepConditional`, `src/maths/regression/MixingLaw.cs:331`) and of the mixture's integral (`NormalVarianceMeanMixture`).
  - z, a term's raw column and a latent that hands a working response (the 2026-11-02 entry).
  - w, BARS's per-point weights, the precisions 1/(c·vᵢ) under a noise (`BarsDraw.Weights`, `src/maths/regression/spline/bars/BarsEnsemble.cs:62`), beside GBASS's w, the global variance multiplier the core writes σ² (the 2026-11-06 entry).
- **Found:** 2026-11-07 and 08, E1's slices 6 and 7.
- **Note:** Yu and Meng (2011) name the two augmentations ancillary and sufficient; the page does too, rather than by their letters.

### α, the miscoverage level, beside Vovk's conformity scores α, and the in-sample predictive CDF F̃ (2026-11-07)

- **Kind:** one name for several concepts, and a symbol taken to keep them apart.
- **Readings:**
  - α, the miscoverage level, with k = ⌊α·m⌋ (`ConformalOptions.Alpha`, `src/maths/regression/conformal/ConformalBayes.cs:20`; `docs/Maths/ConformalPrediction.md` throughout).
  - α_i, a unit's conformity score in Vovk et al.'s conformal predictive distributions, whose level is ε.
  - F̃ᵢ and F̃_c, the in-sample predictive CDF as a score, Σ_s w̃_s F(y | x, θ_s) (`ConformalPrediction.md`, "The predictive distribution (slice 4): the reading"; `tests/regression/ConformalDistributionTests.cs`), taken so it reads as neither α. The page's density scores stay σᵢ.
- **Found:** 2026-11-07, the conformal layer's slice 4 (the predictive distribution's reading).
- **Note:** if the split form's predictive distribution is built, its score F̂(y | x), the training fit's predictive CDF, takes the hat; Vovk's Q(y, τ) would meet the page's Q (next entry).

### σ, the noise's scale, beside the conformity scores σᵢ, and the residual's scale ρ (2026-11-07)

- **Kind:** one name for several concepts, and a third symbol taken to keep them apart.
- **Readings:**
  - σ, σ², the noise's scale and variance across the regression pages and code (`BassDraw.Sigma2`; `docs/Maths/BayesianRegression.md`).
  - σᵢ and σ_c, the full step's conformity scores, Fong and Holmes's symbol (`docs/Maths/ConformalPrediction.md`, "The scores").
  - ρ(x), the standardized residual's scale in the split form, −|y − μ(x)|/ρ(x) (`ResidualReadings.Scale`, `src/maths/regression/conformal/SplitConformal.cs:132`; `ConformalPrediction.md`, "The split form (CB-2)"), taken so it reads as neither σ.
- **Found:** 2026-11-07, the conformal layer's slice 3 (the split form).
- **Note:** Lei et al. (2018) write ρ̂ for the locally weighted form's spread, so ρ is the source's own; the page's Q, the split form's threshold on |y − μ|/ρ and CQR's, is not Vovk's Q(y, τ) for the predictive distribution.

### q and degree in KHAOS beside BASS's K, q and Degree (2026-11-07)

- **Kind:** one name for several concepts, and one concept under two names.
- **Readings:**
  - A term's order, the variables it reads: KHAOS's q, capped at q_max (`LegendreTermSpace.MaxOrder`, `src/maths/regression/spline/bass/LegendreTermSpace.cs:66`; `docs/Maths/PolynomialChaos.md`, "The term space"), is BASS's K (`docs/Maths/BayesianRegression.md`, the order's marginal π(K)), and both are capped by `BassStructureSpec.MaxInteraction` (`src/maths/regression/spline/bass/BassConfig.cs:306`).
  - q in BASS is the hinge power (`BassStructureSpec.Degree`, `BassConfig.cs:303`; `BassBasis.Degree`, `src/maths/regression/spline/bass/BassBasis.cs:41`; "power q" in `BayesianRegression.md`), so KHAOS's q and BASS's q name different things.
  - Degree in KHAOS is a term's total polynomial degree d = Σαⱼ, capped at d_max (`BassTerm.LegendreDegree`, `src/maths/regression/spline/bass/BassTerm.cs:179`; `LegendreBasis.MaxDegree`, `src/maths/regression/spline/bass/LegendreBasis.cs:25`), with each factor's αⱼ (`LegendreFactor.Degree`, `BassTerm.cs:18`). BASS's `Degree` is the hinge power, which a Legendre fit refuses at anything but 1.
- **Found:** 2026-11-07, KHAOS slice 1 (the Legendre basis under the ridge).
- **Note:** the KHAOS paper and `khaos` write q for the order and d for the degree; slice 2 adds KHAOS's g_m and g₀ beside BASS's g-form.

### b, the rung's inverse temperature, beside the GIG's and the GBP's b (2026-11-07)

- **Kind:** one name for several concepts.
- **Readings:**
  - b, the likelihood's power at a rung (`Temperature.InverseTemperature`, `src/maths/samplers/targets/SamplingTarget.cs:67`; `BassCoefficientsBlock.Power`, `src/maths/regression/spline/bass/BassNoise.cs:606`), which the noise's blocks now read: v's conditional p(v)·v^(−b/2)·exp(−b(r − m·v)²/(2s²v)), σ²'s MHN(bN + K + 2α₀, …) and β's precision bΣvᵢ/c + 1/s².
  - b, the GIG's rate on 1/v, GIG(p, a, b) (`GigParameters.B`, `src/maths/samplers/blocks/GigParameters.cs:19`), and the GBP's second shape, GBP(p, a, b) (`GbpParameters.B`, `src/maths/regression/MixingLaw.cs:258`).
  - The brief writes the GIG's as b_v where the two meet, GIG(p − b/2, a + bβ²/c, b_v + b·r²/(w·c)); `docs/Maths/VarianceMeanMixtures.md` does the same in slice 5's section.
- **Found:** 2026-11-07, E1's slice 5 (the engine's reach under a noise).
- **Note:** Friel and Pettitt (2008) write the power posterior's exponent t; D05 and the samplers write b.

### p, a, b of the GBP beside the GIG's, and hbass's c and w against the core's σ² (2026-11-06)

- **Kind:** one name for several concepts, and one concept in two parameterizations.
- **Readings:**
  - p, a, b, the generalized beta prime's power and shapes, GBP(p, a, b) ∝ v^(pa−1)·(1 + v^p)^(−(a+b)) (`GbpParameters`, `src/maths/regression/MixingLaw.cs`; `GeneralizedBetaPrime.P`, `A`, `B`, `src/maths/distributions/GeneralizedBetaPrime.cs`), GBASS's `p`, `a`, `b` of a GBP `v_prior`; against the GIG's p, a, b (`GigParameters`, `src/maths/samplers/blocks/GigParameters.cs`), whose roles differ (a and b there are the rates on v and 1/v), and a, the coefficients.
  - c and w under `hbass`: GBASS takes c = var(y) and w the global scale with √w ~ C⁺(0, 1); the core takes c = 1 (`BassNoise.Horseshoe`, `src/maths/regression/spline/bass/BassNoise.cs`) and σ² under `Sigma2Prior.HalfCauchy()`, so the core's σ² is `hbass`'s w·var(y), the noise's scale per unit v on both sides.
  - The pole's rate, how fast ln g(u) rises from u = −∞ (`NormalVarianceMeanMixture.LeftExponent`, `src/maths/distributions/NormalVarianceMeanMixture.cs`), pa under the GBP and p under the variance gamma; written in words, not α, which names the intercept and the modified half-normal's first parameter.
- **Found:** 2026-11-06, E1's slice 4 (the generalized beta prime on the BASS core).
- **Note:** McDonald (1984) and Armagan, Dunson and Clyde (2011) write the GBP's shapes a and b with the power; `docs/Maths/VarianceMeanMixtures.md` writes GBP(p, a, b) in full where it meets the GIG.

### γ, the Normal-Wald's tail, beside the half-normal's γ and the NIG's α, β, δ (2026-11-06)

- **Kind:** one name for several concepts, and one concept in two parameterizations.
- **Readings:**
  - γ, the Normal-Wald's tail, v ~ GIG(−½, γ², 1) (`ObservationLatents.Gamma`, `src/maths/regression/spline/bass/BassNoise.cs:275`; `BassNoise.GammaPrior`, `BassNoise.cs:197`), GBASS's `gam`, `m_gamma` and `s_gamma`; against the modified half-normal's γ (`ModifiedHalfNormal.Gamma`), which σ²'s block sets to βΣr/c.
  - The NIG's γ = √(α² − β²) (`NormalInverseGaussian.Gamma`, `src/maths/regression/NormalWaldCalibration.cs:17`), Barndorff-Nielsen's, is the core's γ over δ; its β is the core's skew over δ, and its δ the core's σ. `NormalInverseGaussian.CoreSkew` and `CoreGamma` give the core's.
  - ξ = (1 + δγ)^(−½), the steepness of the NIG's shape triangle (`NormalInverseGaussian.Steepness`, `NormalWaldCalibration.cs:29`), GBASS's `nw_triangle`, which `NormalWaldCalibration.GammaPrior`'s targets read.
- **Found:** 2026-11-06, E1's slice 3 (the Normal-Wald on the BASS core).
- **Note:** GBASS writes `gam` and `bet` in the core's parameterization and `alpha`, `beta`, `delta` in `nw_est_mom`'s, the NIG's. `docs/Maths/VarianceMeanMixtures.md` writes the core's γ and β and names the NIG's where they meet.

### q, level, and the modified half-normal's α, β, γ (2026-11-05)

- **Kind:** one concept under several names, and one name for several concepts.
- **Readings:**
  - q, the quantile level, is `BassNoise.Quantile` (`src/maths/regression/spline/bass/BassNoise.cs:189`) and `AsymmetricLaplace.Level` (`src/maths/distributions/AsymmetricLaplace.cs:39`), where `AsymmetricLaplace.Quantile(p)` is the distribution's quantile function, as on every family.
  - Level is also the prior's centre, `BassPriorReference.Level`, which under the quantile read is the sample q-quantile.
  - α, β, γ, the modified half-normal's parameters (`ModifiedHalfNormal.Alpha`, `Beta`, `Gamma`, `src/maths/distributions/ModifiedHalfNormal.cs:42`), against α the intercept and β the centred coefficients ("BASS: the weighted core"), β the skew (`ObservationLatents.Skew`) and γ the Normal-Wald's tail (GBASS's `nwbass`, slice 3). In σ²'s block (`BassScaleBlock`) the half-normal's β is ½(Σrᵢ²/(c·vᵢ) + Q + 2β₀) and its γ is the skew β times Σrᵢ/c.
- **Found:** 2026-11-05, E1's slice 2 (the asymmetric Laplace on the BASS core).
- **Note:** Sun, Kong and Pal (2021) write MHN(α, β, γ); Yu and Moyeed (2001) and GBASS write q. `docs/Maths/VarianceMeanMixtures.md` writes the half-normal's parameters in full where they meet BASS's.

### w, β, τ, c and v between GBASS and the BASS core (2026-11-05)

- **Kind:** one name for several concepts, and one concept under several names.
- **Readings:**
  - w: GBASS's global scale, the core's σ² (`BassDraw.Sigma2`; `gbass-noise-brief.md`, "Why"), against wᵢ and W, the observation weights (`WorkingResponse.Weights`, "BASS: the weighted core").
  - β: GBASS's skew, the drift per unit v (`ObservationLatents.Skew`, `src/maths/regression/spline/bass/BassNoise.cs:284`), against β and β_c, the centred coefficients in (α, β) ("BASS: the weighted core") and β̂ (`BassDesign.Coefficients`).
  - τ: GBASS's τ multiplies the coefficients' variance, a ~ N(0, wτI), so it is the reciprocal of the core's τ, which multiplies their precision.
  - c: GBASS's `scale`, the Gaussian part's variance per unit w·v (`BassNoise.VarianceScale`, `BassNoise.cs:182`), against c, σ²'s shape in `BassResolved`'s (1 − 2c)/N and `BassEvidencePriors` (logged on 2026-11-02 below).
  - v: GBASS's latent scale (`ObservationLatents.Values`, `BassNoise.cs:281`; `IMixingLaw`), against v = ln r, the hierarchical rate's coordinate ("BASS: priors that stop reading a moving response"; `BassState.LogRate`).
- **Found:** 2026-11-05, E1's slices 0 and 1 (GBASS's noise on the BASS core).
- **Note:** the paper and GBASS's code write w, β, τ, c and v as here; the slice's own text writes σ² for w, and m = σβ, s² = σ²c for the mixing law's unit drift and variance.

### Lens and bracket (2026-11-04)

- **Kind:** one concept under several names.
- **Readings:** "lens" and "span lens" in the stopping primer (§13.19, §13.21) and `SpcAnalysis.md`'s open point. "Bracket" in code (`SpcScheduleHelpers.AutoBracket`, `SpanBracket`) and `GraphConstruction.md`'s decision table. "Grid" for the same ends with their levels, in `--temperatures` and `LogSpaceGrid`.
- **Found:** 2026-11-04, building `SpanBracket` from the primer's probe.
- **Note:** a bracket is two ends and a grid adds the levels; the primer's "lens" means both.

### z in a moving response (2026-11-02)

- **Kind:** one name for several concepts.
- **Readings:**
  - z, a term's raw column (`BassDesign.Whitened`, `docs/Maths/BayesianRegression.md`, "BASS: the weighted core": W^½(z − z̄_w)).
  - z, a latent state or value that hands a working response (`moving-response-priors-options.md`: "the latent z (probit)"; "BASS: a design whose working response moves": π(S, z)).
- **Found:** 2026-11-02, a design whose working response moves.
- **Note:** Albert and Chib's own z is the probit's latent, so the column may be the one to rename.

### κ, r and h in the priors for a moving response (2026-11-02)

- **Kind:** one name for several concepts.
- **Readings:**
  - κ, the coefficient prior's constant in σ²κ(τG)⁻¹ (`BassPriorReference.CoefficientConstant`, `docs/Maths/BayesianRegression.md`, "BASS: priors that stop reading a moving response"), against a term's normalizing constant κ ("BASS: sensitivity", "BASS: the extension seams") and the condition number κ(E) ("BASS: the design and admission").
  - r, σ²'s rate under the hierarchical prior (`ResolvedSigma2Rate`, `BassState.LogRate`), which is b in `Sigma2Prior`'s IG(a, b) and g₂ in `GPriorMarginal`'s IG(g₁, g₂), against the residual r throughout `CollapsedGPrior`.
  - h, the rate's shape (`Sigma2Prior.HyperShape`), against the link h(·) in `BassPriorReference.Link`'s text and the count prior's h₁, h₂.
  - c stayed clear of the constant: it is already σ²'s shape in `BassResolved`'s (1 − 2c)/N and `BassEvidencePriors`, where `Sigma2Prior` says a.
- **Found:** 2026-11-02, steps 1 and 2 of the moving-response priors.
- **Note:** the sources' own symbols are c for κ (Sabanés Bové and Held 2011) and 1/ξ for r (Wand et al. 2011). b was not taken for r because b is the inverse temperature in BASS.

### The ridge's names on a prior that is no ridge (2026-11-02)

- **Kind:** one name for several concepts.
- **Readings:** `RidgeColumns`, `RidgeTerms`, `DesignProposal.Ridge` and `BassStructureFit.Ridge` (`src/maths/regression/gprior/CollapsedGPrior.cs`, `src/maths/regression/spline/bass/BassDesign.cs`, `BassTarget.cs`) name the normalized ridge's terms, and now also hold `PriorPencil`'s, under the g-form on the unweighted columns as much as under the ridge. `BassDesign.PriorTerms` is the general spelling.
- **Found:** 2026-11-02, the weighted core.
- **Note:** what they hold is a decomposed pencil's log determinant and shrinkage, so a name for that fits both. `RidgeSpectrum` stays the ridge's own.

### The gauge factor: g, κ, 1/K̂ and the old couplingScale (2026-11-01)

- **Kind:** one concept under several names.
- **Readings:** `κ` ("Normalization κ") in `docs/Graphs/GraphConstruction.md` before the compiler revision, `g` there now, and `GaugeFactor` in `GaugeScaler`, `ScalerResult` and `CouplingResolution`; `1/K̂` in the BWD papers, `BwdPottsCriticalEstimate`'s comments and `issues/spc-BARS/graph-compiler-brief.md`; `couplingScale` (`n / totalNeighbors`) in `GlobalBandwidthScaler`, removed by the revision.
- **Found:** the graph-compiler revision.
- **Note:** K̂ is the mean degree of the returned graph, so g = 1/K̂ only for `GaugeKind.InverseMeanDegree`.

### Bridge, MST bridge and repair edge (2026-11-01)

- **Kind:** one concept under several names.
- **Readings:** "bridge" and "MST bridge" in the compiler brief, `MstBridge`, `MstBridgeReport`, `EdgeWeightSummary.NearZeroBridges` and `ConnectivityRepair`'s "minimal bridges"; "repair edge" in `RepairEdgeSet` and the build log.
- **Found:** the graph-compiler revision.
- **Note:** a repair edge is not always a bridge in the graph-theoretic sense (a cut edge), and a graph-theoretic bridge need not be one the repair stage added.

### Near-zero: three thresholds (2026-11-01)

- **Kind:** one name for several concepts.
- **Readings:** the coupling floor, 1e-10, under which an edge is dropped (`CouplingProjection.CouplingFloor`); the collapse level, 1e-8, which the auto-LMP trigger and `MaxNearZeroEdgeRatio` count against (`AutoThresholds.NearZeroCoupling`); `EdgeWeights.Summary`'s `nearZeroEpsilon`, 1e-10, behind `NearZeroBridges`.
- **Found:** the graph-compiler revision.
- **Note:** only the first two are config fields; the third is a parameter of the summary.

### Pair connectivity: τ, ⟨n⟩, P(i↔j) and C, beside δ̄ (2026-11-01)

- **Kind:** one concept under several names, with a second concept one affine map away.
- **Readings:** τ_ij and τ_e in the stopping primer (§12.21, §13.6) and `issues/spc/bayesian-clustering-design.md`; ⟨n_ij⟩ in `CoMembership` (`src/graphs/primitives/CoMembership.cs:4`); P(i↔j) in `issues/spc/opus-gvdv-ideas.md`. δ̄ = ((q−1)τ + 1)/q, the spin agreement, is the different quantity `DenseTStack` cuts, read off co-membership by `SweepEdgeCurves.CoMembershipDelta`.
- **Found:** reading θ as a cost ratio for the Bayesian clustering note. The ratio belongs to τ, so θ on δ̄ maps through (qθ − 1)/(q − 1).
- **Note:** τ also names `LineagePersistence`'s split threshold in `docs/Clustering/SpcAnalysis.md`, so that doc and the closure gap's docstrings write ⟨n_ij⟩, and the design note writes τ, as the primer does.

### Binder: the loss and the cumulant (2026-11-01)

- **Kind:** one name, two concepts by two authors.
- **Readings:** Binder's loss (D. A. Binder 1978), the expected pair-counting loss with costs a and b, in `issues/spc/opus-gvdv-ideas.md`; the Binder cumulant U = 1 − ⟨m⁴⟩/(3⟨m²⟩²) (K. Binder 1981), `BinderCumulant` (`src/graphs/models/potts/observables/BinderCumulant.cs:16`).
- **Found:** naming the loss for the Bayesian clustering note.
- **Note:** proposed, not decided: pair loss over all pairs, link loss over graph edges (`issues/spc/bayesian-clustering-design.md`). Each would keep "Binder's loss" as its citation only.

### E: the energy and the window's width (2026-10-31)

- **Kind:** one symbol, two concepts in one doc.
- **Readings:** the window's width in e-folds in the stopping primer's §3 ("over a window of E e-folds", n_min = 1 + ⌈2E/w_min⌉) and §13.4 (n ≈ c·E/w_min); the energy in its §5, §12 and §13 (⟨E⟩, Var(E), σ_E, D_J) and in `Accumulator.EnergyPerDraw` (`src/clustering/graphical/spc/runtime/core/sampler/Accumulator.cs`).
- **Found:** writing the primer's round 2 (`issues/spc-BARS/stopping-theory-primer.md`, §13).
- **Note:** §13 keeps §3's E for the window where it quotes §3's count, and the energy everywhere else.

### Repeated runs at one temperature: replica and replicate (2026-10-31)

- **Kind:** one concept, two names.
- **Readings:** "replica" in `SweepCurve.Replicas` (`src/clustering/graphical/spc/profiling/SweepCurves.cs:115`), in `MonteCarloError`'s "the replicas as the chains" (`src/clustering/graphical/spc/profiling/MonteCarloError.cs:10`) and in the stopping primer's §12 and §13 ("8 replicas"); "replicate" in the primer's §1–§2 ("replicate spread"), in `Pointwise.cs`'s "a replicate's observation" (`src/maths/regression/Pointwise.cs:13`) and in the BASS extensions brief's "replicate sweeps" (`issues/regression/bass/extensions-brief.md`, the conformal entries).
- **Found:** the primer's round 2, where the user's own word was "replicant".
- **Note:** physics' replica (an independent copy of the system, as in replica exchange) and statistics' replicate (a repeated measurement) both fit. §13 says "replica" for a run and keeps "replicate sweep" where it follows the extensions brief.

### One SW update: draw, cycle and sweep (2026-10-29)

- **Kind:** one concept under three names, one of which ("sweep") also names a second concept.
- **Readings:** "draw" in `SwendsenWang.Draw` (`src/clustering/graphical/spc/runtime/core/sampler/SwendsenWang.cs:246`), `Accumulator.DrawCount` and the per-draw record `Accumulator.EnergyPerDraw` and siblings (`Accumulator.cs:29`, `:56`); "cycle" in `RunBudget.Cycles` (`RunBudget.cs:9`), `BurnIn(int cycles)` (`SwendsenWang.cs:214`) and the engine's `Σ_cycles` comments; "sweep" in the SW retention briefs' "per-sweep record" and "measured sweep" (`issues/spc/sw-retention-brief.md`, S1). "Sweep" also names the temperature sweep: `FixedGridSweepStrategy`, `SweepProfile`, `SweepCurves`, `SweepResult`.
- **Found:** naming the per-draw record for the SW retention chip (S1).
- **Note:** the record took "draw" (`…PerDraw`), after `DrawCount`.

### Coefficients: `Beta` beside `Coefficients` (2026-10-29)

- **Kind:** one concept, two identifiers on one type.
- **Readings:** `BassDesignSnapshot.Beta` and the field `_beta` (`src/maths/regression/spline/bass/BassDesign.cs:30`, `:113`); the public `BassDesign.Coefficients` (`:165`).
- **Found:** listing β's readings for this ledger.

### b_k beside the inverse temperature b (2026-10-29)

- **Kind:** one symbol, two concepts in one doc.
- **Readings:** the move schedule's birth and death probabilities b_k, d_k (r_k the rest) in BayesianRegression.md's "DMGK schedule"; the inverse temperature b in the same doc's tempering sections ("raised to an inverse temperature", "Tempered marginal (step 3.1)").
- **Found:** checking the doc's b uses for the f5ba714 question.
- **Note:** b_k, d_k is Green's (1995) own notation, which DMGK keep.

### The inverse temperature: b, β and T (2026-10-29)

- **Kind:** one concept, several names.
- **Readings:**
  - b / `inverseTemperature`: the sampler contract (`src/maths/samplers/targets/SamplingTarget.cs:7`, Samplers.md's target and temperature bullets) and the regression code and docs (`CollapsedGPrior.cs:47`, `ITemperedMarginal.cs:19`, `BassTarget.cs:194`, BayesianRegression.md).
  - β / `Betas` / `betaMin` / `temperBetaMin`: the tempering ladder (`ParallelTempering`, `TemperingLadder`, Samplers.md's tempering section) and `BarsEnsemble`, until f5ba714 moved them to b / `InverseTemperatures` / `hotEnd` / `temperHotEnd`. That commit stands until this pass, which may revert it; `tests/regression/BassMixingExitLadderTests.cs` was written against the new names after it.
  - T, the temperature itself: SPC's temperature grid and the Potts/SW code, with β = 1/T in `src/clustering/dendrogram/PartitionHierarchy.cs:101` and the SPC notes (`issues/spc/opus-spin-glass-threads.md`, `issues/spc/sw-retention-brief.md`).
- **Found:** 2026-10-28, writing E13 in the extensions brief; f5ba714.
- **Note:** the sources write β (Watanabe 2013, Farhi et al. 2013, Friel and Pettitt 2008); Bissiri et al. 2016 write w for the same weight. The SPC↔regression bridge (temperature integration, frustrated clustering) puts SPC's T and the samplers' b in one derivation.

### β: four readings (2026-10-29)

- **Kind:** one symbol, several concepts.
- **Readings:**
  - Regression coefficients: BARS, BAPS, HPS and BASS code and BayesianRegression.md; `HierarchicalPSpline.cs:160` (β_r, β_0).
  - The inverse temperature (entry above).
  - A spectral exponent, PSD ∝ |k|^−β: `src/maths/samplers/SpectralNoiseField.cs:11`.
  - The Beta distribution, by name: `AdaptiveLocalBetaKernel`, the `"beta"` kernel strings in `CountLawTests`.
- **Found:** 2026-10-28, the f5ba714 survey.
- **Note:** BayesianRegression.md's "Tempered marginal (step 3.1)" has coefficient β and inverse temperature b in one sentence.

### η: six readings (2026-10-29)

- **Kind:** one symbol, several concepts.
- **Readings:**
  - The GLM linear predictor, η = Bβ: `src/maths/regression/spline/bars/PoissonModel.cs:13`, `WhittleModel.cs:12`, BayesianRegression.md's `PoissonModel` bullet.
  - BGP's noise-to-signal ratio, σ² = s²η: `src/maths/regression/bgp/GpRegression.cs:18`, `BgpSampler.cs:13`, PosteriorContraction.md's prior [P1].
  - BAPS's local smoothing multipliers ηᵢ ~ Gamma(φ, φ): `src/maths/regression/spline/baps/Baps.cs:32`, BayesianRegression.md's `LocallyAdaptiveGibbs`.
  - HPS's parameter-expansion deviations η_r: `HierarchicalPSpline.cs:160`, BayesianRegression.md's `ParameterExpanded`.
  - A step size η_k in the Riemannian subgradient solver: `src/maths/geometry/solver/irls.cs:138`, RiemannianGeometry.md's solver flowchart.
  - A centred factor column of the Khatri–Rao product: `src/maths/linalg/KhatriRao.cs:16`, `BassDesign.cs:449`.
- **Found:** 2026-10-28, checking whether η meant a learning rate anywhere (E13).
- **Note:** the step size is the nearest thing to a learning rate. E13 writes the generalized-Bayes loss weight as b, the inverse temperature it is.
