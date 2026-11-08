# Polynomial chaos

> KHAOS on the BASS core (E19, with E4 and E21): an adaptive Bayesian polynomial chaos expansion, every term a product of standardized shifted Legendre polynomials over the continuous predictors, with reversible jump over KHAOS's degree-split term space and its Bernoulli-inclusion proposal, under the normalized ridge or KHAOS's raw ridge. Slice 1 was built on 2026-11-07. The modified g-prior (slice 2), ridge blocks (3), readouts past the paper (4) and the engine's reach (5) follow.

## Overview

KHAOS (Rumsey, Francom et al., *Bayesian Adaptive Polynomial Chaos Expansions*, arXiv 2510.25036; CRAN `khaos` 2.1.1, the oracle) is BASS's lineage on three axes: the basis, the proposals and the coefficient prior. It moves no latent and is conjugate, so it is a diff against the engine and not a new engine. The basis is a part type behind [the design-kernel contract](./BayesianRegression.md) (F-X3). The proposal fills F-X4, the proposal-kernel contract. The ridges are the core's collapsed blocks. So the sweep, the target, the design's admission, the run, the fingerprint and the readers serve a Legendre fit as they serve a hinge fit. The build is [the brief's](../../issues/regression/bass/khaos-brief.md), slice by slice. Slice 0's oracle modes sit under `r/oracles/khaos_*_oracle.R`, with `khaos`'s samplers corrected for the defects slice 0 confirmed (`khaos_corrected`), and those corrected samplers are the reference the engine meets.

| Piece | Where | Slice |
|---|---|---|
| ψ_α, values, slopes, lanes, expansions | `Maths.Distributions.ShiftedLegendre` | 1 |
| The Poisson-binomial | `Maths.Distributions.PoissonBinomial` | 1 |
| The Legendre part | `LegendreFactor` on `BassTerm`; `LegendreProductKernel` among `DesignKernels.Case`; `LegendreDesignFactor` | 1 |
| The term space and its measure | `LegendreTermSpace`; `BassStructurePrior`'s Legendre part (`BassPart.Legendre`) | 1 |
| F-X4 and KHAOS's proposal | `IBassProposalKernel`; `LegendreTermProposal` (D01's `BassTermProposal` the hinge member) | 1 |
| KHAOS's raw ridge | `RawRidgeForm`, through `BassDesign`'s prior-apart pencil | 1 |
| The config | `BassStructureSpec.Legendre` (`LegendreBasis`), `BassCoefficientPrior.RawRidge` | 1 |
| Declared input domains | `InputDomain` on `ContinuousPredictor`; `RangeScale` scales by it | 1 |
| The modified g-prior, g₀² by τ's slice step | `ModifiedGForm`; `BassPriorPrecision`, the design's prior-precision builder (E21); `TauPrior.Gamma` | 2 |
| Sobol and the Jacobi measure from the coefficients | `BassSensitivity.FromCoefficients` → `PolynomialChaosSensitivity` | 4 |
| The gradient readout | `BassFit.Gradient` → `BassGradient`; `BassSurface.EvaluateDerivative` | 4 |
| The engine's reach | The ladder, SMC, the flat histogram, the pointwise and candidate densities and GBASS's noise laws, unchanged; gated on Legendre terms | 5 |
| Random-effect intercepts (E21) | `BassRandomIntercept`, `RandomInterceptBlock`: a capability of the core, gated on a hinge fit and a Legendre fit ([BayesianRegression](./BayesianRegression.md#bass-random-effect-intercepts-e21)) | 3 |

## The basis (E19)

**ψ_α(x) = √(2α + 1)·P_α(2x − 1)**, orthonormal on [0, 1] under the uniform measure. A term is Πⱼ ψ_{αⱼ}(xⱼ) over its variables, with no normalizer: unlike a hinge product, which is scaled to peak at one, it already has unit mean square under the uniform measure.
- *Values* come from Bonnet's recurrence, (k + 1)P_{k+1} = (2k + 1)uP_k − kP_{k−1}.
- *Slopes* come from P′_{k+1} = P′_{k−1} + (2k + 1)P_k. That recurrence holds at u = ±1, where the closed form k(uP_k − P_{k−1})/(u² − 1) divides by zero, so a derivative at an end of the training range needs no special case.
- *A column is evaluated in lanes,* four points at a time (`ShiftedLegendre.Multiply`, `Vector256`, `AggressiveOptimization`). The recurrence's coefficients aₖ = (2k + 1)/(k + 1) and bₖ = k/(k + 1) are tabled, so a step multiplies and never divides. The column runs one step at a time across all its points, so the points' steps overlap, where one point's chain of α steps would wait on each. Each element sees the scalar form's operations in its order, so the bits are the scalar form's on any hardware. That is about 6× the scalar per-point dividing form at n = 1,000 (the cost row).
- *No column cache.* The brief allowed a per-variable cache of ψ_α up to d_max, chosen by size. A degree-15 column still costs about 4.5× to 6× a two-hinge product's, by run, about a third of a normalized-ridge sweep at the cost report's base. A cache would cut that share by evaluating each variable's ψ once per fit. It is left as an option, with the reading that would decide it.

**Inputs** are range-scaled to [0, 1] as BASS's are, (x − min)/(max − min), recorded on the fit, over each predictor's declared domain when it declares one (`InputDomain`) and over its training range otherwise: the brief's sixth decision, "inputs on [0, 1] by the fit's declared ranges". `khaos` expects its inputs on (0, 1) already, so a user who knows the domain declares it and the engine reads the values `khaos` reads. **The posterior is not invariant to the box.** On slice 0's Ishigami fixture, `khaos` on the raw inputs, which span about [0.008, 0.992], reads M 10.2 and σ² 0.23 (noise 0.25). On the same data rescaled to [0, 1] exactly, `khaos` and the engine alike read M about 6.2 and σ² about 0.95, and eight engine chains of 35,000 steps agree. A polynomial space of total degree d is invariant to an affine map, but a sparse selection from it is not: under a shift and a stretch, each ψ_α spreads over the ψ_β with β ≤ α. The likely cause is Ishigami's exact symmetry about ½, which keeps its raw expansion sparse (odd degrees in x₁, even in x₂ and x₃) and which a shift of its centre breaks. That is unverified.
- *The domain is the data's.* `InputDomain` sits on `ContinuousPredictor` (an output variable's too), and every value must lie in it. `RangeScale` keeps the training range beside the box, so prediction still flags a case outside the training range, where the fit extrapolates. The fingerprint takes a domain only when one is declared, so data that declare none keep theirs. A sensitivity measure left undeclared is uniform over the box, and `ContinuousMeasure.Extrapolates` reads the training range in the box's units.
- *Hinge fits under a domain:* the knots are the data's, and the g-form and the normalized ridge read each column up to its scale, so the law is the training range's (2,263 enumerated structures agree to 8·10⁻¹⁶). Two readings do move: D06's column mass, which reads hinges normalized over the box, and an undeclared measure, which spans the box.
- *The alternative:* the training range alone, with a user who knows the domain rescaling before the fit. It was the build until the user leaned toward declaring the domain (2026-11-07): the rescaling would sit outside the fit, so its record would not say which box the basis is orthonormal on.
- *Implication:* a Legendre fit without a declared domain reads the training range, which moves with the sample, and its resolution log says so.

**Moments** go through the moment service, as every term type's do. A factor is one piece of a `PiecewisePolynomial`, expanded about the measure's centre (`ShiftedLegendre.Coefficients`), and the uniform measure gives a mean of 0 and orthonormal pairs to 10⁻¹⁰ through degree 10. The expansion's coefficients grow like 2^(2α)/√α, so a moment past degree 15 or so loses digits to cancellation. Under the uniform measure the moments are exact by orthonormality without the expansion: slice 4's Sobol-from-coefficients path reads them that way, and also builds `MixtureMeasure` and `EmpiricalMeasure`.

**A derivative** (`BassBasis.CaseDerivative`) takes the slope of the part that reads the variable, the new `IDesignKernel.ReadsVariable`, and the other parts' values. The old rule took every differentiable part's derivative, which zeroed a Legendre term's through the empty hinge part. A hinge term's derivative is the same bits.

## The term space and its measure

**𝒜_{p, d_max, q_max}** is every multi-index over the p continuous predictors with order q = |{j : αⱼ > 0}| ≤ q_max and degree d = Σαⱼ ≤ d_max. A term is its variable set with a composition of d over it, so

  |𝒜| = Σ_{q ≤ q_max} Σ_{d ≤ d_max} C(p, q)·C(d − 1, q − 1) = Σ_{q ≤ q_max} C(p, q)·C(d_max, q),

the inner sum by the hockey stick. The term measure is uniform over 𝒜, π(term) = 1/|𝒜|, and the count prior is BASS's, λ ~ Gamma(h₁, h₂) integrated. `khaos`'s defaults are h₁ = 4 and h₂ = 40/n.
- *A uniform draw* (`LegendreTermSpace.Draw`, the prior's `DrawTerm`): q with probability C(p, q)·C(d_max, q)/|𝒜|, exact in integers while |𝒜| fits a long; the variables a uniform q-subset; the degrees a uniform q-subset s₁ < … < s_q of {1, …, d_max}, read as α₁ = s₁ and αₖ = sₖ − sₖ₋₁. That is uniform over the compositions of every d ≤ d_max.
- *Uniform compositions:* q − 1 distinct cuts from {1, …, d − 1}, each composition with probability 1/C(d − 1, q − 1) (the brief's third decision). `khaos`'s `random_partition` sorts cuts drawn with replacement and weighs compositions unequally, and its ratios use C(d, q). Slice 0 confirmed both, and both are corrected here.
- *The alternative measure* was BASS's for hinges: the order uniform, then the terms of that order uniformly. KHAOS's uniform measure weighs each order by the terms it holds, which at d_max 15 and q_max 5 puts most of the prior mass on high orders. The count prior's Occam factor and the marginal's penalty then carry the sparsity. Uniform over 𝒜 is the paper's model and the oracle's, so the engine takes it.
- *Refusals:* a Legendre fit holds Legendre terms alone. Categorical predictors, a functional response, a d_max below q_max, the hinge degree, the column-mass support rule, D01's local steps and the fixed-M family are each refused with the reason. Bases that mix term types are an option the contract allows, not built.

## The proposal: F-X4 and KHAOS's Bernoulli inclusion (E4)

**F-X4** (`IBassProposalKernel`) is the contract by which a structure proposal states its forward and reverse probabilities exactly, given the current structure's counts:
- a birth's term, and ln q_b(term | counts);
- a death's ratio, the removed term's birth under the reduced counts against the death's 1/M;
- a change of one term, with its within-move ratio.

The sweep reads its moves through it. D01's hinge proposal is one member, its draws unchanged, and KHAOS's is the other.

**The birth.**
1. q₀ is drawn from {1, …, q_max} ∝ q₀^(−s_q).
2. Each variable's inclusion probability is η = `make_weights`(appearance counts, q₀), ported as `khaos` writes it: vⱼ = (aⱼ^α + ε)/Σ(a^α + ε)·q₀/p, then three passes of β = mean((ln(q₀ + δ) − ln p)/ln v), δ += q₀ − Σv^β, giving ηⱼ = vⱼ^β.
3. χⱼ ~ Bern(ηⱼ), redrawn while q = Σχ is 0 or above q_max.
4. d is drawn from {q, …, d_max} ∝ (d − q + 1)^(−s_d).
5. A uniform composition of d is spread over the chosen variables.

q₀ is drawn once, before the redraws, so the inclusion's exact law is

  P(χ) = Σ_{q₀} J(q₀)·P(χ | q₀)/P(1 ≤ q ≤ q_max | q₀),

each validity the Poisson-binomial's, by dynamic programming over the p coins in O(p·q_max) (`PoissonBinomial`). The forward probability is P(χ)·D(d | q)/C(d − 1, q − 1).
- *The alternatives* were `khaos`'s two forms. Its birth multiplies every attempt's probability in, invalid draws included. Its death approximates the validity by a normal, with up to 64% error in the tail at small p. Both miss the law their own draw follows, which slice 0 confirmed by enumeration. Exact validity on both sides was the brief's fourth decision.
- *`make_weights` as written:* its sum is q₀ only approximately (to about 10⁻⁴). The proposal reads the same η both ways, so the ratio stays exact whatever the sum.

**The death** picks a term uniformly. Its reverse is the same expression under the counts without it.

**The change** is one of KHAOS's two mutations of a uniformly picked term:
- *The degree,* with probability ½ (`LegendreBasis.DegreeMutation`), or always when p ≤ 3 (KHAOS's rule). The degree is redrawn from the birth's law at the term's order, with a new uniform composition; the ratio is each side's D(d | q)/C(d − 1, q − 1).
- *A variable swap* otherwise. A factor is picked uniformly, and its variable is swapped for one the term does not use, ∝ w + its appearances, the factor's degree kept. The reverse draws among the variables the swapped term leaves out, under the counts after the swap. `khaos` zeroes only the swapped-out variable's weight, so it can swap onto a variable the term holds, giving a repeated variable outside the basis. Slice 0 counted 55.6% of a seeded chain's kept terms repeating one. A term over every variable has no swap: a null move.
- *The alternative* was KHAOS's adaptive mutation probability, a rule on the two mutations' acceptance rates with floor 0.1. It is available under D10 as diminishing adaptation or burn-in only, and not built.

## The coefficient priors

**The normalized ridge** (`RidgeForm`) is the core's: σ²I/τ on the centred columns scaled to unit root mean square. It is decomposed once per structure, and each proposal is a bordered candidate.

**KHAOS's raw ridge** (`RawRidgeForm`, `BassCoefficientPrior.RawRidge`): (β₀, β) | σ², τ ~ N(0, σ²I/τ) on the columns as evaluated, uncentred, with the intercept inside under the zero centre. This is `adaptive_khaos_ridge`'s model, with `khaos`'s τ² = 1/τ. In the design's coordinates, with the columns centred, β₀ = α − z̄ᵀβ. So it is a prior apart from the Gram, `PriorPencil`'s D with D_α = 1, δ = −z̄ and D_β = I. The design carries it as it carries the robust row's weighted prior, at any weights, and admission is unchanged. Under level centring β₀'s prior centres on the level, and a flat intercept leaves β₀ out.
- *Cost:* a pencil re-triangulates per proposal, with no bordering, O(M³) per candidate against the normalized ridge's O(M²) (the cost row).
- *The alternative* was a bordered pencil, a rank-one update of L⁻¹AL⁻ᵀ per birth. Not built.
- *When each serves:* on a Legendre basis the columns sit near unit root mean square, so the two ridges nearly agree. The raw ridge is the oracle's model, for parity. The normalized ridge is scale-free across term types.

**KHAOS's modified g-prior** (`ModifiedGForm`, `BassCoefficientPrior.ModifiedG`): (β₀, β) | σ², g₀² ~ N(0, σ²S₀) on the columns as evaluated, the intercept's column of ones among them, with S₀ = g₀²·D_g(ΨᵀΨ)⁻¹D_g, D_g = diag(1, g₁, …, g_M) and g_m = (1 + q_m(d_m + q_m − 2))^(−ζ/2). A richer term's coefficient is shrunk harder. This is `adaptive_khaos_gprior`'s model, with its defaults ζ = 1 and g₀² ~ IG(0.001, 1,000).
- *Through the builder (E21).* The design holds a prior-precision builder, `BassPriorPrecision`, one per coefficient form, in place of a flag per form. It says whether the prior stands apart from the Gram, what each term's prior-side column is, the pencil's D_α and D_β, and whether the angle gate reads the prior's Gram. In the design's coordinates the modified g-prior is the g-form's prior-apart pencil on the prior-side columns scaled by 1/g_m: D_α = N, δ = z̄_u/g − z̄_w and D_β = D_g⁻¹G_uD_g⁻¹. Scaling a column leaves its angles alone, so F is the g-form's. The raw ridge is the builder's unit member, D = I. The g-form and the normalized ridge keep their own forms, a factor and a spectrum, wherever the prior does not stand apart, so their draws are unchanged.
- *g₀²'s step is τ's.* τ = 1/g₀², and IG(a, b) on g₀² is Gamma(a, rate b) on τ: `TauPrior.Gamma(a, b)`, stated on τ itself whatever the form's τ₀. τ's slice step on ln τ (D04) is then g₀²'s on ln g₀², reflected, against the conditional with β and σ² integrated: the IG prior, |Σ_n|^½|S₀|^(−½) and the data term. The pencil is decomposed once per structure, after which each g₀² costs O(p). The draws record τ, and g₀² = 1/τ.
- *ζ = 0* is the g-form with the intercept inside under the zero centre. The two agree at every enumerated structure and τ (the checks).
- *The alternative* was `khaos`'s own step: an independence proposal from a Laplace approximation fitted to a conditional without the data term, which slice 0 confirmed as a defect. The brief chose the slice step on the derived conditional (E19).
- *Refused:* a hinge basis, since the factor reads a Legendre term's order and degree; and, under weights, the weighted prior Gram, which moves the prior with W. The unweighted and frozen prior Grams serve it.
- *Where τ is subnormal.* `khaos`'s IG(0.001, 1,000) on g₀² puts about half its mass on τ below 10⁻³⁰⁰, where the density is finite. Two products overflowed there: the spectral terms' ln(1 + bλ/τ), and Goldberg's ln(1 + x) above x ≈ 2.5·10³⁰⁵. Each now takes ln bλ − ln τ or ln(1 + x) where it overflowed, and keeps its operations elsewhere. A chain started from such a draw climbs to the posterior's τ in about a hundred sweeps, where before it could not move.
- *The box under the g-prior.* On Ishigami the g-prior reads M about 15 and σ² 0.20 on the declared box, and M about 20 and σ² 0.22 on the rescaled one; the raw ridge read 10 and 0.23 against 6 and 0.95. A term costs about ½ln(1 + g₀²g_m²) + ln|𝒜| ≈ 3 + 6 nats there, against the raw ridge's 14, so the rescaled box's denser expansion still pays.

σ²'s prior is the core's; `khaos` defaults to g₁ = g₂ = 0 for the ridge's σ² and a_σ = b_σ = 0 for the g-prior's, Jeffreys both.

## The config

`BassStructureSpec.Legendre` is a `LegendreBasis`; null, the default, is BASS's hinges. It holds d_max (`MaxDegree`, `khaos`'s 15), s_q (`OrderPenalty`, 1), s_d (`DegreePenalty`, 0), `make_weights`' α, ε and passes (1, 2, 3), the swap's pseudo-count (1) and the mutation probability (½). q_max is `MaxInteraction`, a term's continuous factors, resolved to at most p. `Legendre` and `RawRidge` are left out of the JSON when unused, so a hinge config's form and fingerprint stand as they were. A Legendre fit's `BassModel.Proposal` (D01's) is refused, and `BassModel.Kernel` is its proposal.

## Readouts past the paper (slice 4)

Each readout holds under one measure, read by one route:

| Readout | Measure | Route | Exact |
| --- | --- | --- | --- |
| `BassSensitivity.FromCoefficients` → `PolynomialChaosSensitivity.Sobol`, the moment route's own type | The product uniform law on each input's box: its declared domain, else its training range | Σ b_α² over each subset's terms, O(M) per draw | Yes, on a Legendre fit alone |
| `PolynomialChaosSensitivity.Jacobi`, E[x_v(1 − x_v)(∂_v f)²] | The same law | Σ b_α² α_v(α_v + 1) | Yes |
| `BassFit.Sensitivity` (`BassSensitivity.Of`) on a Legendre fit | Any product measure: a uniform on any interval, `MixtureMeasure`, `EmpiricalMeasure` | The moment service: each factor's moments through its polynomial expansion and the measure's interval moments, O(M²) per draw | Yes; a Student component's moments are taken by Gauss–Legendre, their error reported (`MomentError`) |
| `BassFit.Gradient` → `BassGradient`, ∂f/∂z per draw, case and location in each predictor's own units | None: pointwise | Each term's derivative through the design kernels (O-X4), and the chain rule through the box | Yes; at a hinge's knot, the mean of the one-sided limits |

**From the coefficients** (Lean SG-4(b), `BassIntegrals.lean`, `faca090`). On the product uniform law the ψ_α are orthonormal with mean 0 for α ≠ 0, so for f = b₀ + Σ_α b_α ψ_α:
- Var f = Σ b_α². A subset's main variance is the sum over the terms whose variables are exactly that subset. An input's total is Σ_{α_v≥1} b_α², which is E[(f − E_v f)²] itself (`tensorTotal_eq`).
- The Jacobi measure is E[x_v(1 − x_v)(∂_v f)²] = Σ b_α² α_v(α_v + 1) (`integral_jacobi_tensorPartial`).
  - Under the affine map to raw units it reads E[(z − a)(b − z)(∂_z f)²]: the same number.
- Each total is at most half its Jacobi measure, with equality where every term holding the input has degree 1 in it (`tensorTotal_le`, `tensorTotal_eq_half`).
  - This is the Poincaré-type check E8's entry asks for, read per draw. On the Ishigami fixture at degree 8, ½·Jacobi over the total was 2.79 or more.
- The decomposition is the moment route's on the same draws, at O(M) per draw in place of O(M²). The moment route stays the default (`BassFit.Sensitivity`), since it serves every measure and every term type.

**Under any other measure the coefficient formulas do not hold.** There the ψ_α are neither orthonormal nor of mean 0, so `FromCoefficients` takes no measure. Its readout says the law it holds under.
- *The decomposition* is the moment route's, exact through the measure's interval moments. KHAOS's "past the paper" is this route on Legendre factors, and no code changed for it.
- *A derivative measure* is the gradient readout averaged over the measure's points. That route serves any measure, a dependent one included.
- *Refused:* `FromCoefficients` on a fit with a hinge or categorical term, whose coefficients are not on an orthonormal basis, and on a functional response.

**For E8.** E8's C = E[∇f∇fᵀ] is the moment service with derivative factors, under a product measure. That is not built here: it is E8's own.
- What E8 can take from this slice:
  - the gradient readout, pointwise and so good under any measure, which is the fallback for the first consumer that cannot assume a product measure;
  - the Jacobi diagonal in closed form under the uniform law;
  - its bound against the totals.
- The lineage review's open question on the consumer's shape is still settled at E8.

**The cost** (wall clock, one run, indicative): on the comparison's 80 draws with M up to 18, the moment route took 18.1 ms and the coefficients 0.7 ms, each with its setup.

## The engine's reach (slice 5)

A Legendre term is a part type behind F-X3, so every sampler and reader serves it as it serves hinges, through the design and the target. Slice 5 changed no engine code; it gated each one on KHAOS's priors (the checks below):
- the ladder, with and without its reference rung, and with its rungs stepped at once (`BassRunOptions.RungWorkers`);
- SMC, under the fixed and the adaptive count of sweeps;
- the flat histogram;
- the pointwise density, both contents, and the conformal layer's candidate density;
- GBASS's noise laws, the skewed horseshoe with its move along β·v among them.

Two things to know when fitting:
- *Under a noise the modified g-prior reads Ψ as evaluated,* not the latents' weights, as `khaos`'s S₀ does. The noise's prior Gram is the unweighted one, which the builder serves.
- *A noise's reference and the zero centre.* KHAOS's fits take the zero centre. A noise's reference (`BassNoiseReference.Moments`, the default) reads a level from the response, and the zero centre refuses it at the fit. A KHAOS fit under a noise states `BassNoiseReference.None`, with a σ² prior that reads nothing (Jeffreys, `khaos`'s).

**SMC on KHAOS: the trap, measured; its remedy is the user's open decision.** On Ishigami under the raw ridge (τ held at 10⁻⁵, degree 8, order 3, n 120, N 500; `KhaosSmcTrapQualificationTests`):
- *The mechanism is the schedule's jump, not the kernel at b = 1.*
  - Below b ≈ 0.005 a term costs more than it gains, so the tempered law sits at near-empty structures. The population collapses there to one to nine distinct structures of 500 particles, which all reweight alike.
  - The conditional ESS then sees no spread, and the schedule jumps to b = 1 in one stage, from b between 0.003 and 0.35.
  - ln Ẑ is read at that jump, about 60 nats low, with an error under 1.
  - At b = 1 itself the posterior is reachable: greedy forward selection from M = 0 climbs monotonely, and the ladder's cold chains read the long chains' M and σ².
- *How often.*
  - The engine's schedule reached the posterior's basin, σ² within twice the long chains', in 3 of 32 replicates at 5 sweeps a stage, and 2 of 16 at 20.
  - Capping each step at b′ ≤ 1.25b reached the basin in 11 of 16 at twice the time, and b′ ≤ 1.1b in 8 of 16 at 3.5 times. The capped runs still end short of the posterior (M 7 to 10 against 10.15). Their ln Ẑ spreads over 49 nats, each replicate reporting an error under 1.
- *KHAOS's own g-prior is not trapped there.* A term costs about 9 nats against the raw ridge's 14, and 16 of 16 replicates reach σ² 0.19 to 0.21. Its ln Ẑ still spreads over 43 nats (sd 11), each replicate reporting an error under 1.
- *The ladder has the same transition and says so.*
  - Its cold chains read M 10.13 to 10.16 and σ² 0.225, the long chains' values.
  - One pair in every chain rejects every swap, and there are no round trips. The four chains' stepping stones read −296 to −206.
  - The evidence reports itself saturated, with a between-to-within error ratio of 28.
- *What sees it.*
  - "The last step left from b < ½": 41 of 43 trapped replicates and none of the 5 found.
  - `LeastMoved` under 0.1 at the last three stages: 40 of 43 trapped, and 1 of 5 found.
  - A drift in the mean log posterior over two extra stages at b = 1: none.
- *Decided 2026-11-08: flag it.* The coordinating session made the call, the user unsure and leaving it to that session; the user can reverse it.
  - BASS's SMC fit now flags a replicate whose last step jumped over a collapsed population, and replicates whose ln Ẑ disagree beyond their own errors ([Bayesian regression](./BayesianRegression.md#bass-smc), "The flags"). Neither changes a draw.
  - On this fixture, as production fits:
    - *Raw ridge, 32 replicates:* 24 of the 30 trapped replicates jumped, and neither of the 2 that reached the posterior did. The other six read as disagreement (between-to-within 30), and every group of four was flagged.
    - *The g-prior, 16 replicates:* all reached the posterior's σ², none jumped, and they disagreed at 9.2 times their errors.
- *Recorded, not built,* each waiting on its own decision ([the brief](../../issues/regression/bass/khaos-brief.md#reports)):
  - a cap on the step: measured above, it reaches the basin but not the evidence;
  - a birth proposed from the data's residual through F-X4, aimed at the 2% birth acceptance;
  - a reference nearer the posterior, a pilot run and a tempering path from it.

## What was checked (slices 1, 2, 4 and 5)

| Gate | Reading |
|---|---|
| ψ_α | Orthonormal by Gauss–Legendre quadrature to degree 40 (10⁻¹²); zero at the Gauss–Legendre nodes to degree 64; endpoint values and slopes to degree 60; Bonnet's identity on a seeded sweep to degree 60 (10⁻¹²); lanes and batch bit for bit against the scalar |
| The Poisson-binomial | Against brute force over 2ᵖ inclusion vectors at chosen coins and a 400-case sweep across twelve decades, zero and one included (10⁻¹² relative); the binomial at 2,000 coins (10⁻⁹); the truncated recursion bit for bit |
| \|𝒜\| | The brief's sum, a direct count of multi-indices and the enumeration agree at six sizes; the cost report's p 20, d_max 15, q_max 5 exactly; p 400, q_max 12 in logs (10⁻¹²) |
| The draws | The uniform term draw and the uniform composition by χ²; the birth's draws against its stated law (400,000 draws); each mutation's draws against its statement |
| The forward law | The inclusion against brute force (10⁻¹²); the birth's probabilities sum to one over 𝒜 (10⁻¹²); `make_weights` and the forward law against `khaos` at slice 0's frozen states (3·10⁻¹⁶ and 7·10⁻¹⁶) |
| The exact law | p 4, d_max 3, q_max 2, M ≤ 2: 466 structures. π from the target against a dense marginal from the columns alone (6·10⁻¹⁵ under either ridge) and against `khaos_exactlaw_oracle.R` (3·10⁻¹⁵); every move's kernel from its statement in detailed balance and stationary (4·10⁻¹⁴; 1,800 birth/death, 1,740 degree, 3,132 swap transitions); the production chain's P(M), inclusions and six likeliest structures against π over 400,000 steps (worst \|z\| 1.8 raw, 0.7 normalized) |
| Geweke | 12 cases, M ≤ 3, σ² ~ IG(3, 2), τ held, p 4: worst \|z\| 2.15 and 1.39 (raw ridge, two seeds), 1.31 (normalized) over M, σ², the surface at two points, two inclusions, the degree sum and a scaled residual |
| Readouts | Against corrected `adaptive_khaos_ridge` on Ishigami (degree 8, order 3; R 6,000 iterations with 3,000 kept, the engine two chains of 10,000 after 10,000): on the raw inputs with [0, 1] declared, M 10.15 against 10.15, σ² 0.2275 against 0.2278, the surfaces 0.005 apart (root mean square, function sd 3.72), total Sobol indices within 0.001; on the inputs rescaled to [0, 1], M 6.22 against 6.13, σ² 0.929 against 1.005, the surfaces 0.035 apart, the indices within 0.004. Eight engine chains of 30,000 on the raw inputs read M 10.13 and σ² 0.2282, and two `khaos` runs of 30,000 read 10.14 to 10.15 and 0.2274 to 0.2275 |
| Declared domain | The exact law with [0, 1] declared on the raw rows against `khaos_exactlaw_oracle.R` (3·10⁻¹⁵); with inputs in [0.05, 0.95] and [0, 1] declared, π against a dense marginal on ψ at the raw values (3·10⁻¹⁵), where the training range's law differs by up to 0.019; refusals, the box and the flags, the fingerprint, the measure and the log (`InputDomainTests`) |
| The modified g-prior's target | At six structures and g₀² from 0.5 to 10⁶, for ζ = 0, 1 and 2.5: the log marginal against a dense multivariate t built from the columns alone, y \| g₀² ~ t with scale I + B S₀ Bᵀ (1.3·10⁻¹⁰ after one constant), so g₀²'s conditional, which τ's step reads, is the joint's |
| ζ = 0 | The pencil against the g-form's `CollapsedGPrior` with the intercept inside, both with g₀² drawn: 466 enumerated structures at five τ each, 4·10⁻¹⁶ relative |
| The exact law, g₀² integrated | Each of the 466 structures' weight by Simpson's rule over ln τ on the oracle's nodes, against `khaos_exactlaw_oracle.R`: 1.6·10⁻¹⁵ in P(M) and the inclusions and 1.7·10⁻¹⁵ in g₀²'s posterior CDF at 300, 1,000 and 10⁴, on both boxes. The production chain with g₀² ~ IG(0.001, 1,000) drawn, against it over 400,000 steps: worst \|z\| 1.46 |
| Geweke, g₀² drawn | 12 cases, p 4, M ≤ 3, σ² ~ IG(3, 2), g₀² ~ IG(2, 2) drawn by τ's step; 200,000 marginal-conditional draws against 300,000 successive steps over M, σ², ln g₀², the surface at two points, two inclusions, the degree sum and a scaled residual: worst \|z\| 2.37 and 1.52 (ζ = 1, two seeds), 2.14 (ζ = 0) |
| `adaptive_khaos_gprior` readouts | Corrected, degree 8, order 3, R 6,000 iterations with 3,000 kept against the engine's two chains of 10,000 after 10,000. On the raw inputs with [0, 1] declared: M 14.66 against 13.87, σ² 0.201 against 0.203, g₀²'s median 3,137 against 3,256, the surfaces 0.068 apart, total indices within 0.003. On the rescaled inputs: M 20.23 against 19.08, σ² 0.221 against 0.236, the surfaces 0.074 apart, the indices within 0.003. Long runs on the declared box: eight engine chains of 40,000 read M 15.06 (sd 0.36), σ² 0.197 and g₀²'s median about 3,130; two `khaos` runs of 40,000 read M 15.06 and 15.03, σ² 0.198 and g₀²'s median 3,125 and 3,115 |
| Sobol from the coefficients (slice 4) | Against the moment route, draw for draw, on an Ishigami fit (degree 8, order 3; 80 draws, M up to 18, 7 subsets): 6·10⁻¹⁴ of Var f. Against `sobol_khaos` on the same 80 draws, packed as khaos's own object (`khaos_sobol_oracle.R`, `KhaosSobolParityTests`): 2.2·10⁻¹⁶ in every index, its division by Var f + σ² included |
| SG-4(b)'s identities | Against a 10³ Gauss–Legendre rule, exact at these degrees, on 8 draws each of fits at degree 8, order 3 and at degree 2, order 2, from the draws' own surfaces and gradients: Var f, the main and total variances and the Jacobi measure to 4·10⁻¹⁵. Total ≤ ½·Jacobi at every draw and input, the least ratio 2.79 at degree 8; equality to 10⁻¹² at the 4 (draw, input) pairs held to degree 1 |
| The gradient | Against central differences of the predicted surface at 10⁻⁶ of each range, on 40 cases: a Legendre fit on [−π, π]³ with that domain declared, so the chain rule through the box shows, 3·10⁻¹⁰ relative; a hinge fit, 2·10⁻¹⁰ |
| Moments under other measures | The moment route's decomposition of a Legendre fit (degree 4, order 2) against brute-force tensor sums over each measure, its density written out in the test: under the empirical marginals, all 27,000 triples of training values, 3·10⁻¹⁴ of Var f; under truncated mixtures (normal and Student components, both weightings) on 40 density-weighted Gauss–Legendre nodes per input, 3·10⁻¹⁴, the Student moments' own quadrature error 7·10⁻¹² |
| The law at each rung (slice 5) | Slice 0's exact-law fixture with σ² ~ IG(2, 0.5), under the raw ridge at τ held and the modified g-prior with g₀² ~ IG(2, 80) drawn: four ladders of 40,000 scans, with and without the reference rung, read from their own replicas against the enumeration's tempered law. The count law at every rung within 5 SE + 0.004, and E[ln τ] at every rung within 5 SE + 0.01, from b = 1 (−3.512 against −3.510) to the prior at b = 0 (−3.959 against −3.959) |
| SMC and the ladder | ln p(y) against the enumeration's: SMC (2,000 particles, 4 replicates) −72.297 ± 0.014 against −72.291 on the raw ridge, −69.617 ± 0.017 and, under the adaptive count, −69.599 ± 0.023 against −69.587 on the modified g-prior; the ladder (reference rung, self-sizing) −72.258 ± 0.036 and −69.606 ± 0.043. At 32 replicates of 8,000 particles SMC's count law sits within 1.2 SE of the posterior with g₀² held, and within 0.72 SE with it drawn, its ln p(y) within 0.8 SE with it drawn. A tempered fit with its rungs stepped at once is the fit stepped in order, hash for hash |
| The flat histogram | Against the enumeration (`BassFlatHistogramTests.MatchesTheEnumeration`): the learned weights 0.039 (modified g-prior, 12 nodes) and 0.020 (raw ridge, 6 nodes) from SAMC's limit, the reweighted count law, inclusions, likeliest structures and E[τ] within their allowances. At a held τ the density of states is F's 466 point masses, which leave SAMC's limit undefined on 12 nodes |
| The pointwise and candidate densities | Both contents on draws over 13 structures with g₀² drawn and 4 at τ held: the retained content against the Normal at the explicit uncentred design (4·10⁻¹⁶); the collapsed content against explicit refits under KHAOS's prior written out, τI on the uncentred coefficients or τ·D_g⁻¹ΨᵀΨD_g⁻¹ (1.4·10⁻¹⁴), each τ its own leave-out quadratics. The candidate density at the training cases is the pointwise column under the Gaussian and Student t noise, exactly |
| GBASS's noise on Legendre terms | GBASS's per-rung Geweke (`BassNoiseReachQualificationTests`) on Legendre terms of degree up to 3 under the modified g-prior with g₀² drawn: the Student t at b = 1 and 0.4 (worst \|z\| 2.92 and 2.07), the quantile at b = 1 interwoven (1.74), and the skewed horseshoe with its move along β·v at b = 1 and ¾ (0.96 and 2.82), no move refused |
| Conformal coverage | The full conformal Bayes set on a Gaussian KHAOS fit under contaminated noise no model here holds, twenty fits of 99 with fifty held-out cases each: coverage 0.878 against the nominal 0.900 (se 0.011), mean width 1.08, no set open at the grid's ends |
| A Student t and a quantile fit | A Student t(3) fit under 0.25·t₂ noise, sixteen sets of 150: its surface beats the Gaussian fit's on all sixteen, mean RMS error 0.088 against 0.188. A quantile fit at q = 0.8 on heteroscedastic noise, four sets of 300: 0.768 of 1,200 held-out cases below its surface (± 0.012), the fitted surface 0.079 from the true 0.8-quantile surface in root mean square against the Gaussian mean surface's 0.289 |
| SMC's trap | See the section above: the engine's schedule jumps over a near-empty population on the raw-ridge Ishigami fixture, and the ladder flags the same transition (`KhaosSmcTrapQualificationTests`). SMC's flags, as production fits there: 24 of 30 trapped replicates jumped and no replicate that reached the posterior did; between-to-within 30 (raw ridge) and 9.2 (g-prior); every group of four flagged |
| Hinge BASS | The whole regression assembly at the slice's end, 2,253 of 2,253, the eight pinned hashes and every fingerprint among them, and the commit's selection on the later HEAD, 1,849 of 1,849; a hinge config's JSON holds neither `Legendre`, `RawRidge` nor `ModifiedG`, and a τ prior other than the Gamma on τ writes no rate; at slice 5's end, on HEAD `477834c` with its tests, 2,361 of 2,361, the pinned hashes among them |

**The cost row** (`KhaosCostTests`, the thread's cycle clock, read round-robin in one process; indicative, the machine loaded):
- *A column* at n = 1,000, in two runs, the second on a quieter machine: ψ₈·ψ₇ 11.9 and 8.6 µs and ψ₁₅ 10.4 and 7.3 µs, against 66 and 46 µs for the scalar per-point form that divides at each step, and 2.6 and 1.4 µs for a product of two hinges. ψ₈·ψ₇ reads 4.5× and 6.2× the hinges.
- *A sweep* at the cost report's base (n 1,000, p 10, M_max 20, q_max 2, d_max 15), at M 20: 37.7 and 18.9 µs under the normalized ridge and 356 and 190 µs under KHAOS's raw ridge, against 22.3 and 10.2 µs for the hinge fit there. The raw ridge's pencil decomposes each candidate anew, the brief's accepted design for a prior apart from the Gram, which slice 2's g-prior shares.
- *A sweep at p 20, d_max 15, q_max 5,* where the chains sat at M 2 to 4: 38.2 and 16.3 µs (raw) and 48.5 and 20.0 µs (normalized).
- *The modified g-prior* (slice 2, a third run): at the base, a sweep 343 µs at M 20, beside the raw ridge's 363 and the normalized ridge's 35 in the same run, since each candidate is decomposed afresh as under the raw ridge; at p 20, q_max 5, 58 µs at M 6

## Notation

This page takes KHAOS's symbols. α is a multi-index, q a term's order (its variables), d its degree (Σα), q_max and d_max their caps, and 𝒜 the term space. Where they cross BASS's:
- *q:* KHAOS's q is a term's interaction order, BASS's K. BASS's q is the hinges' power, the config's `Structure.Degree`.
- *d / degree:* KHAOS's degree is a sum over the term's factors. BASS's `Degree` is the hinge power, so a Legendre fit refuses it and caps its degree by `Legendre.MaxDegree`.
- *g:* KHAOS's g₀² is the engine's 1/τ, as the g-form's g is in σ²G_c⁻¹/τ, so the two meet at ζ = 0. g_m is a term's complexity factor, with no BASS counterpart. `khaos`'s g₁ and g₂ are the ridge's σ² hyperparameters.

The code follows the page: `LegendreFactor.Degree` is α, `LegendreDegree` is d, `MaxInteraction` is q_max. These crossings are in [the notation ledger](../../issues/notation-ledger.md).
