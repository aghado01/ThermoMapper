# Several outputs on one basis (E34)

> K continuous outputs per case that share one structure, with their K × K noise covariance Σ integrated out: one fit selects terms for every output at once and reads their residual correlation. Built slice by slice from [the brief](../../issues/regression/bass/multi-output-brief.md), which holds the decisions; [E34's entry](../../issues/regression/bass/extensions-brief.md#several-outputs-on-one-basis-with-σ-integrated-e34) holds the model's place among the extensions. This page extends BASS's collapsed core ([BayesianRegression.md](./BayesianRegression.md)) and does not repeat it. Ordered outputs on a grid stay with the direct and PCA modes.

## The model and its collapse (decision 1; slice 1)

Y (N × K) = 1μᵀ + Z_cA + E, each row of E ∼ N(0, Σ), with the intercept μ ∈ ℝᴷ placed by BASS's convention per output (D03).
- *The priors.* A | Σ, τ ∼ MN(0, (τG_c)⁻¹, Σ) on the centred columns, and Σ ∼ IW(ν₀, Ψ₀) in the conventional degrees of freedom. vec stacks A's columns, one per output, so Cov(vec A) = Σ ⊗ (τG_c)⁻¹. Inclusion is shared across the outputs: the structure is one (decision 4).
- *The collapse.* With A, the intercept and Σ integrated against the likelihood raised to b, the structure step reads the K × K scale S = b(R + τ/(τ + b)·(F + N(ȳ − ℓ)(ȳ − ℓ)ᵀ)) and K times the scalar determinant term, where R = Σᵢ rᵢrᵢᵀ and F = B_cᵀG_c⁻¹B_c are the centred least-squares fit's residual and fitted forms (the scalar ‖r‖² and Q_c carried to K outputs) and ℓ the intercept's centre levels, so the last term is the zero centre's Nȳȳᵀ and none under level centring. S is a sum of semidefinite terms, so a small τ on a close fit keeps its digits as the scalar S does. Σ | Y ∼ IW(ν₀ + ν, Ψ₀ + S), and
  ln p_b(Y | structure, τ) = −(νK/2) ln π + K·(scalar determinant term) + (ν₀/2) ln|Ψ₀| − ln Γ_K(ν₀/2) + ln Γ_K((ν₀ + ν)/2) − ((ν₀ + ν)/2) ln|Ψ₀ + S|,
  exactly 0 at b = 0 under a proper prior. Untempered it is the matrix-t marginal, Brown, Vannucci and Fearn's (1998) eqs. (19)–(21); at K = 1 it is the scalar block's marginal under IG(ν₀/2, Ψ₀/2) (`OutputMarginal`).
- *The count* ν is Nb under level centring and the zero centre, whose intercept sits inside the prior, and Nb − 1 under the flat intercept, whose integral spends an observation. *Corrected 2026-11-08* from "n − 1" in the brief's first draft, at slice 0: the K = 1 reduction to scalar BASS requires Nb, and Brown, Vannucci and Fearn's vague intercept, α ∼ N(α₀, hΣ) with h → ∞, counts N too, its |Σ|^(−½) cancelling in α's integral. A count b(N − 1) would be a third path, the likelihood of Y projected off 1 tempered, which is neither BASS's integrated form nor equal to scalar BASS at K = 1.
- *Alternatives.* Separate scalar fits per output read no residual correlation and select K structures. Inclusion per output needs an indicator per (term, output) and loses the matrix-t (decision 4, not built). The response envelope and an output basis from E[J W_x Jᵀ] are E34's later options.

## Σ's prior (decision 2; slice 1)

**Decision 2, the default:** ν₀ = K − 1 + 2c and Ψ₀ = 2c·S_Y, c = 10⁻³, stated as an empirical prior. S_Y is scalar BASS's s² carried to K outputs: the outputs' mean squares and cross-products about the intercept's centre levels with divisor N, each entry by s²'s own arithmetic (ȳ under level centring and the flat intercept, 0 under the zero centre; *corrected 2026-11-08* from "sample covariance about their means").
- *At K = 1* it is IG(c, c·s²), the scalar default: S_Y is s² to the bit, and the resolved prior's density is the inverse gamma's.
- *Units.* S_Y → LᵀS_YL when Y → YL, so ln p(YL) = ln p(Y) − (Nb − f) ln|det L| at every structure and τ, f = 1 under the flat intercept: the structure and τ posteriors are invariant to any invertible rescaling or mixing of the outputs.
- *Its correlations.* Under a diagonal S_Y each prior correlation has density ∝ (1 − ρ²)^(c − 1): Huang and Wand's (2013) Property 3 at their ν = 2c, read through a correlation's invariance to a diagonal rescaling. At c = 10⁻³ nearly all of its mass sits at ±1, Gelman's (2006) objection to IG(ε, ε) carried to correlations. The data overrule it at moderate N/K; with few cases per output they may not. It stays the default provisionally, for the reduction and the invariance, and slice 3 measures the coverage of the correlations and of the joint predictive at N/K from 3 to 30 against Huang and Wand's prior at ν = 2.
- *Refused:* an output with no spread about its level, as scalar BASS refuses a constant response; and under this form, outputs whose S_Y is singular or nearly so, read on the correlation scale as a Cholesky pivot 1 − R² under 10⁻¹⁰. Collinear outputs need fewer outputs.

**The options** (`BassConfig.OutputCovariance`, `OutputCovariancePrior`; null carries `Sigma2`'s ScaledToResponse(c), and Sigma2's other kinds are refused with the reason):
- *Standardized(c):* Ψ₀ = 2c·diag S_Y, the default on outputs standardized to unit mean square: invariant to rescaling the outputs, not to mixing them, and centred on uncorrelated outputs. Collinear outputs pass.
- *A declared IW(ν₀, Ψ₀)* in the outputs' own units: free of the data, as Geweke's test needs.
- *Huang and Wand's prior at ν = 2* (slice 2): Σ | a ∼ IW(ν + K − 1, 2ν diag(1/a)), a_k ∼ IG(½, 1/A_k²) with A_k each output's root mean square about its level. It gives uniform correlations and half-t standard deviations; at K = 1 it is the engine's hierarchical σ² prior at hyper-shape ½ (`Sigma2Prior.Hierarchical(ν/2, ½)`, A = s), so built it is that prior carried to K rates, each moved by a slice step as τ is, the collapse holding given them. Equivariant under rescaling and permuting the outputs, not under mixing them.

## Properness under tempering (decision 3)

At rung b, Σ's conditional is IW(ν₀ + ν_b, Ψ₀ + S_b), proper where ν₀ + ν_b > K − 1. With the intercept inside a proper prior that holds at every b ≥ 0, the reference rung included. Under the flat intercept it fails at Nb ≤ K − ν₀ under a proper prior too (Nb ≤ 1 − 2c under decision 2), and b = 0 diverges, as scalar BASS's flat pole does (`BassTarget.Defines`). The entry's first condition, ν₀ − N(1 − b) > K − 1, asked the tilted prior IW(ν₀ + N(b − 1)) alone to be proper, which integrating Σ against that tilted prior needs but integrating the tempered likelihood directly does not. The samplers that read it arrive in slice 2.

## The design's K-output side (slice 1)

A design of several outputs keeps the scalar design's structure, factor, admission and gate, which read no response, and holds the response's statistics in a K-output side (`DesignOutputs`) in place of the scalar ones: the outputs centred on their means, each term's K cross products, β̂ (M × K), the explicit residuals and the forms R and F.
- *Proposals* carry the scalar bordered updates to K × K. A birth or a change reads R′_kl = R̃_kl − z̃_k z̃_l/(newest·δ) on the reduced model's R̃ and z̃, a change's gain β̂_mβ̂_mᵀ/(G⁻¹)_mm and β̂_m·uᵀG⁻¹e_m/(G⁻¹)_mm; a death's R′ = R + β̂_mβ̂_mᵀ/(G⁻¹)_mm. Where a move explains nearly all of one output's residual the candidate's residuals are formed explicitly, as the scalar design forms its own. F′ is read from the candidate factor's whitening of each output's crosses, F′_kl = (L′⁻¹b′_k)·(L′⁻¹b′_l) (`BorderedCholesky.ProposedWhiten`, beside `ProposedInverseQuadraticForm`).
- *Acceptances* solve β̂ and recompute every residual explicitly, so nothing drifts; snapshots carry the side for a restart bit for bit.
- *Scalar BASS is unchanged by construction:* every K-output path is gated on the side, and the pinned tier passes at each commit.
- *What it costs:* a proposal adds O(NK) for the crosses and O(M²K + MK²) for F′ over the scalar step, an acceptance O(NMK) for the residuals; the hinge prefix tables stand aside for several outputs in slice 1, since their borders are approximate and the K-output side reads exact ones. The cost row is below.

## The draw (slice 1)

A retained draw is exact given (structure, τ): Σ | Y from its inverse Wishart by Bartlett's decomposition through Ψ's own factor (`InverseWishart`), then A | Σ ∼ MN(β̂/(τ + 1), (τ + 1)⁻¹G_c⁻¹, Σ), and the intercepts α | Σ ∼ N(m, Σ·c₀/N) apart from them: m = ℓ + (ȳ − ℓ)/(τ + 1) and c₀ = 1/(τ + 1) inside the prior, m = ȳ and c₀ = 1 under the flat intercept (`BassOutputDraws`). The stream is read as the scalar draw reads it, Σ's gamma variate and normals, then each output's M normals through the factor's innovation, mixed across the outputs by Σ's factor, then the intercepts' K normals, so at K = 1 a fit is scalar BASS's draw for draw on one seed. `BassDraw` holds Σ and the coefficients as one row per term and one column per output (`OutputCoefficients`), uncentred in the ids' order.

## The run (slice 1)

A run of several outputs monitors ln|Σ| in σ²'s place, and each monitored case's fitted value at every output: the monitored points are (case, output) pairs, a declared case without an output filling in every output. ln|Σ| reads the correlations as well as the scales, and its rank-normalized R̂ is moved by no rescaling or mixing of the outputs. The alternative, Σ's whole lower triangle, reads each entry apart at the cost of K(K + 1)/2 functionals and another shape for the diagnostics; it is a follow-up if a run's correlations converge more slowly than ln|Σ| shows.

## What slice 1 refuses, with its reason

Noises and random intercepts (decisions 6 and 7, slice 4); a working response that moves; the ladder, SMC and the flat histogram (slice 2); the normalized ridge and KHAOS's priors, whose spectra and pencils would carry K cross products (slice 2); a direct functional response; missing outputs (decision 5); and the readers of one response (prediction, the pointwise density, sensitivity, the conformal candidates), which refuse a draw of several outputs until the readouts (slice 3).

## Checked

- **Slice 0's references** (`MatrixTReferenceTests`): the closed form against the dense vec form with the Kronecker products built entry by entry (Λ₀⁻¹ ⊗ Σ in Σ ⊗ Λ₀⁻¹'s slot misses by 0.27 to 4 nats), quadrature over Σ at K = 2 (10⁻¹³ nats), Monte Carlo over the prior at K = 3, the matrix-t density, the scalar block at K = 1 (574 cases, 4·10⁻¹⁴), Brown, Vannucci and Fearn's eq. (20) over 64 structures (2.8·10⁻¹⁴) and Huang and Wand's Property 3.
- **The pieces** (`OutputMarginalTests`): the inverse Wishart's density against the reference (8·10⁻¹⁵) and its draws against the inverse gamma's at K = 1 on one stream (10⁻¹⁴); the K-output marginal against the reference over 548 cases with K to 4, the three conventions, a frozen level, τ over seven decades, b to 0 and units 10^±6 (3.4·10⁻¹³, refusals agreeing).
- **The core** (`MultiOutputCoreTests`): the design against one scalar design per output and per pair of outputs, by polarization, over 232 proposals and 149 acceptances (7.5·10⁻¹⁵); the target against the matrix-t reference under each convention (1.3·10⁻¹⁵); K = 1 against scalar BASS, the forms to the bit, the target to 3.6·10⁻¹⁶ and the draw to 6·10⁻¹⁶; equivariance to 10⁻¹⁵.
- **Fits** (`MultiOutputFitTests`, `OutputCovariancePriorTests`): at K = 1 a fit is scalar BASS's, 500 draws, structures and τ equal, σ², coefficients and λ to 2.5·10⁻¹⁴; each prior kind resolved and round-tripped through JSON, the target under a declared prior against the reference (10⁻¹²).
- **Laws** (`MultiOutputLawTests`): the chain's count law, E[ln τ] and six likeliest structures against the enumeration of F with τ integrated, K = 2 and 3, 200,000 steps (worst |z| 3.0 and 2.6); the recovery of a known, correlated Σ and of three surfaces outside the model (each Σ entry within 1.3 posterior sd; the surfaces' RMSE 0.04 to 0.06 against noise sds 0.15 to 0.3). Geweke's joint test and the cost row run under qualification; their numbers are in the brief's report.

## Notation

Σ is the K × K noise covariance and ν, ν₀ the inverse Wishart's conventional degrees of freedom; Brown, Vannucci and Fearn write IW(δ; Q) in Dawid's shape, δ = ν − K + 1, and their q is this page's K. R and F are the K × K forms of the scalar ‖r‖² and Q_c, and S the K × K scale of the scalar S. Where K, F, R, S_Y and ν cross the module's own names is logged in [the notation ledger](../../issues/notation-ledger.md) (2026-11-09).
