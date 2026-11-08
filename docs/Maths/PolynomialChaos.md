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
| The modified g-prior, g₀² by slice steps | `PriorPencil` as a prior-precision builder | 2 |

## The basis (E19)

**ψ_α(x) = √(2α + 1)·P_α(2x − 1)**, orthonormal on [0, 1] under the uniform measure. A term is Πⱼ ψ_{αⱼ}(xⱼ) over its variables, with no normalizer: unlike a hinge product, which is scaled to peak at one, it already has unit mean square under the uniform measure.
- *Values* come from Bonnet's recurrence, (k + 1)P_{k+1} = (2k + 1)uP_k − kP_{k−1}.
- *Slopes* come from P′_{k+1} = P′_{k−1} + (2k + 1)P_k. That recurrence holds at u = ±1, where the closed form k(uP_k − P_{k−1})/(u² − 1) divides by zero, so a derivative at an end of the training range needs no special case.
- *A column is evaluated in lanes,* four points at a time (`ShiftedLegendre.Multiply`, `Vector256`, `AggressiveOptimization`). The recurrence's coefficients aₖ = (2k + 1)/(k + 1) and bₖ = k/(k + 1) are tabled, so a step multiplies and never divides. The column runs one step at a time across all its points, so the points' steps overlap, where one point's chain of α steps would wait on each. Each element sees the scalar form's operations in its order, so the bits are the scalar form's on any hardware. That is about 6× the scalar per-point dividing form at n = 1,000 (the cost row).
- *No column cache.* The brief allowed a per-variable cache of ψ_α up to d_max, chosen by size. A degree-15 column still costs about 4.5× to 6× a two-hinge product's, by run, about a third of a normalized-ridge sweep at the cost report's base. A cache would cut that share by evaluating each variable's ψ once per fit. It is left as an option, with the reading that would decide it.

**Inputs** are range-scaled to [0, 1] as BASS's are, (x − min)/(max − min) on the training data, recorded on the fit (the brief's sixth decision). `khaos` expects its inputs on (0, 1) already. **The posterior is not invariant to that map.** On slice 0's Ishigami fixture, `khaos` on the raw inputs, which span about [0.008, 0.992], reads M 10.2 and σ² 0.23 (noise 0.25). On the same data rescaled to [0, 1] exactly, `khaos` and the engine alike read M about 6.2 and σ² about 0.95, and eight engine chains of 35,000 steps agree. A polynomial space of total degree d is invariant to an affine map, but a sparse selection from it is not: under a shift and a stretch, each ψ_α spreads over the ψ_β with β ≤ α. The likely cause is Ishigami's exact symmetry about ½, which keeps its raw expansion sparse (odd degrees in x₁, even in x₂ and x₃) and which a shift of its centre breaks. That is unverified.
- *The alternative:* declared input ranges, the domain's box beside the training range's. Not built, and recorded as an open question in the slice 1 report.
- *Implication:* a fit is not invariant to its input box, so a domain known to the user is worth stating. Until declared ranges exist, a Legendre fit reads the training range.

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

σ²'s prior is the core's; `khaos` defaults to g₁ = g₂ = 0, Jeffreys. The modified g-prior, KHAOS's other prior, is slice 2's.

## The config

`BassStructureSpec.Legendre` is a `LegendreBasis`; null, the default, is BASS's hinges. It holds d_max (`MaxDegree`, `khaos`'s 15), s_q (`OrderPenalty`, 1), s_d (`DegreePenalty`, 0), `make_weights`' α, ε and passes (1, 2, 3), the swap's pseudo-count (1) and the mutation probability (½). q_max is `MaxInteraction`, a term's continuous factors, resolved to at most p. `Legendre` and `RawRidge` are left out of the JSON when unused, so a hinge config's form and fingerprint stand as they were. A Legendre fit's `BassModel.Proposal` (D01's) is refused, and `BassModel.Kernel` is its proposal.

## What was checked (slice 1)

| Gate | Reading |
|---|---|
| ψ_α | Orthonormal by Gauss–Legendre quadrature to degree 40 (10⁻¹²); zero at the Gauss–Legendre nodes to degree 64; endpoint values and slopes to degree 60; Bonnet's identity on a seeded sweep to degree 60 (10⁻¹²); lanes and batch bit for bit against the scalar |
| The Poisson-binomial | Against brute force over 2ᵖ inclusion vectors at chosen coins and a 400-case sweep across twelve decades, zero and one included (10⁻¹² relative); the binomial at 2,000 coins (10⁻⁹); the truncated recursion bit for bit |
| \|𝒜\| | The brief's sum, a direct count of multi-indices and the enumeration agree at six sizes; the cost report's p 20, d_max 15, q_max 5 exactly; p 400, q_max 12 in logs (10⁻¹²) |
| The draws | The uniform term draw and the uniform composition by χ²; the birth's draws against its stated law (400,000 draws); each mutation's draws against its statement |
| The forward law | The inclusion against brute force (10⁻¹²); the birth's probabilities sum to one over 𝒜 (10⁻¹²); `make_weights` and the forward law against `khaos` at slice 0's frozen states (3·10⁻¹⁶ and 7·10⁻¹⁶) |
| The exact law | p 4, d_max 3, q_max 2, M ≤ 2: 466 structures. π from the target against a dense marginal from the columns alone (6·10⁻¹⁵ under either ridge) and against `khaos_exactlaw_oracle.R` (3·10⁻¹⁵); every move's kernel from its statement in detailed balance and stationary (4·10⁻¹⁴; 1,800 birth/death, 1,740 degree, 3,132 swap transitions); the production chain's P(M), inclusions and six likeliest structures against π over 400,000 steps (worst \|z\| 1.8 raw, 0.7 normalized) |
| Geweke | 12 cases, M ≤ 3, σ² ~ IG(3, 2), τ held, p 4: worst \|z\| 2.15 and 1.39 (raw ridge, two seeds), 1.31 (normalized) over M, σ², the surface at two points, two inclusions, the degree sum and a scaled residual |
| Readouts | Against corrected `adaptive_khaos_ridge` on rescaled Ishigami (degree 8, order 3, 6,000 iterations, 3,000 kept): M 6.32 against 6.13, σ² 0.959 against 1.005, the surfaces 0.049 apart (root mean square, function sd 3.72), total Sobol indices within 0.002 |
| Hinge BASS | The whole regression assembly, 2,172 of 2,172, the eight pinned hashes and every fingerprint among them; a hinge config's JSON holds neither `Legendre` nor `RawRidge` |

**The cost row** (`KhaosCostTests`, the thread's cycle clock, read round-robin in one process; indicative, the machine loaded):
- *A column* at n = 1,000, in two runs, the second on a quieter machine: ψ₈·ψ₇ 11.9 and 8.6 µs and ψ₁₅ 10.4 and 7.3 µs, against 66 and 46 µs for the scalar per-point form that divides at each step, and 2.6 and 1.4 µs for a product of two hinges. ψ₈·ψ₇ reads 4.5× and 6.2× the hinges.
- *A sweep* at the cost report's base (n 1,000, p 10, M_max 20, q_max 2, d_max 15), at M 20: 37.7 and 18.9 µs under the normalized ridge and 356 and 190 µs under KHAOS's raw ridge, against 22.3 and 10.2 µs for the hinge fit there. The raw ridge's pencil decomposes each candidate anew, the brief's accepted design for a prior apart from the Gram, which slice 2's g-prior shares.
- *A sweep at p 20, d_max 15, q_max 5,* where the chains sat at M 2 to 4: 38.2 and 16.3 µs (raw) and 48.5 and 20.0 µs (normalized).

## Notation

This page takes KHAOS's symbols. α is a multi-index, q a term's order (its variables), d its degree (Σα), q_max and d_max their caps, and 𝒜 the term space. Where they cross BASS's:
- *q:* KHAOS's q is a term's interaction order, BASS's K. BASS's q is the hinges' power, the config's `Structure.Degree`.
- *d / degree:* KHAOS's degree is a sum over the term's factors. BASS's `Degree` is the hinge power, so a Legendre fit refuses it and caps its degree by `Legendre.MaxDegree`.
- *g:* KHAOS's g_m and g₀ (slice 2) meet BASS's g-form, σ²G_c⁻¹/τ.

The code follows the page: `LegendreFactor.Degree` is α, `LegendreDegree` is d, `MaxInteraction` is q_max. These crossings go to [the notation ledger](../../issues/notation-ledger.md) once its open edit is committed.
