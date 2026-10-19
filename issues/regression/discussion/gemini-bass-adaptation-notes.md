
BASS


regression


samplers


Incorporating BASS (Bayesian Adaptive Spline Surfaces
, Francom & Sansó / LANL) into ThermoMapper as a new pillar of Maths.Regression is very feasible and architecturally well-suited.

Because ThermoMapper has already generalized its trans-dimensional MCMC engine (ReversibleJumpChain,ParallelTempering, and ChainEnsemble), approximately 50–60% of the complex MCMC machinery needed by BASS is already implemented and verified.

Here is an architectural assessment of where BASS fits, which primitives are shared, what new components must be built, and an estimated implementation roadmap.

1. Architectural Positioning: Where BASS Fits in ThermoMapper
ThermoMapper’s src/maths/regression
 currently has three primary pillars:

- spline/bars: 1D free-knot B-spline RJMCMC (DiMatteo, Genovese & Kass) — engineered for high-resolution 1D curve recovery, non-local peak extraction ($T_c$ transition detection, $\lambda(T)$ intensity, $\pi(T)$ span coverage) with banded linear algebra.
- spline/baps: 1D/tensor penalized B-splines with adaptive difference penalty $\lambda$ (and local adaptive $\eta_i$) via Gibbs/REML marginal evidence.
- bgp: Nonparametric stationary Gaussian process with empirical-Bayes bandwidth sampling.

## The Role of BASS
BASS represents Multivariate Nonparametric Regression (Bayesian MARS): $$f(\mathbf{x}) = a_0 + \sum_{m=1}^M a_m B_m(\mathbf{x}), \quad B_m(\mathbf{x}) = \prod_{k=1}^{K_m} \frac{\left[ s_{m,k} (x_{v_{m,k}} - t_{m,k}) \right]_+^q}{c_m}$$

It addresses limitations in the existing pillars:

- **Curse of Dimensionality:** Full Cartesian tensor products (in BAPS) scale exponentially with input dimension $p$. BASS adaptively selects sparse interaction sets ($K_m \le \text{maxInt}$) over continuous and categorical inputs.
- **Global Sensitivity Analysis (GSA):** Because products of truncated power hinges integrate in closed form against uniform/mixture priors, BASS provides exact analytical Sobol indices ($S_j, S_{j,k}, T_j$) without Monte Carlo integration.
- **Functional Output Emulation:** High-dimensional curve/surface outputs $y_i(t)$ are modeled either by design-functional tensor products or via PCA projection (bassPCA).

Recommended placement: src/maths/regression/spline/bass/
 (or src/maths/regression/bass/).

### Primitives Shared With Existing Infrastructure

| Subsystem | Existing ThermoMapper Primitive | How BASS Reuses It | Effort Saved |
| --- | --- | --- | --- |
| RJMCMC Engine | 

ReversibleJumpChain<TState>
The core Green 1995 4-factor acceptance engine is completely model-agnostic. BASS just defines its state carrier BassState and implements 

IRjMove<BassState>
	High (engine already handles move-selection ratios, reverse keys, and proposal stats)
Parallel Tempering	

ParallelTempering<TState>
BASS in R relies heavily on parallel tempering (temp.ladder) to traverse multimodal basis spaces. ThermoMapper already implements replica-exchange parallel tempering conforming to IChain<TState>.	High (instant multi-temperature swap support)
Multi-Chain Ensemble	

ChainEnsemble
Provides multi-core fan-out (Parallel.For), round-based burn-in, early stopping, and Gelman–Rubin $\hat{R}$ / ESS convergence diagnostics across invariant functionals.	High (all convergence diagnostics & seed trees are solved)
Reproducibility	

SeedTree
, 

Xoshiro256PlusPlus
Hierarchical deterministic seeding across replicas and tempering chains.	High
Linear Algebra	

CholeskyDecomposition
In-place lower-triangular Cholesky factorizations with log-determinant and soft floor clamping.	Medium
Functional PCA	

PCA
SVD/PCA basis decomposition for multivariate/functional outputs (bassPCA).	High
Distributions	Maths.Distributions	Gamma, Normal, Inverse-Gamma, Beta, and Truncated distributions for hyperparameter priors ($\sigma^2$, $\lambda$, $\tau$).	Medium

### 3. Critical Linear Algebra Nuance: Dense vs. Banded Cholesky
In univariate BARS and BAPS, B-spline basis functions have local compact support ($p+1$ overlapping neighbors), producing banded Gram matrices $Z^T W Z$ that use 

BandCholesky
 with $O(n \cdot \text{deg}^2)$ LDLᵀ solves.

In BASS, basis functions are multivariate tensor-product hinges with global support. Therefore:

$X^T X$ is dense, not banded.
However, the active basis count $M$ is typically small ($M \approx 10$ to $60$).
A $60 \times 60$ dense Cholesky solve in C# takes less than a microsecond.
Optimization opportunity:
**Birth move**: Appends 1 basis function $\implies$ bordered Cholesky rank-1 update ($O(M^2)$ triangular solve).
**Death move**: Removes 1 basis function $\implies$ standard Cholesky column deletion / Givens downdate ($O(M^2)$).
**Change move**: Modifies 1 basis function $\implies$ rank-1 downdate + rank-1 update.

## What Needs to Be Built (The Gaps)
**Carrier Representation** (BassBasis, BassState):

**Continuous hinge**: $[s(x_v - t)]_+^q / c$ with degree $q \ge 1$ and normalization constant $c$.
**Categorical indicator**: $I(x_v \in S)$ with subset $S \subset \text{Levels}$.
**Tensor product term** combining $K \in [1, \text{maxInt}]$ continuous and categorical factors.
**BassState**: immutable or copy-on-write set of $M$ active basis terms.
**Reversible Jump Moves** (IRjMove<BassState>):

**BassBirthMove**: Samples interaction order $K$, variable subset (Wallenius / weighted without replacement), signs $s_k \in {-1, +1}$, and knot $t_k$ from candidate data points. Computes proposal density $q(m)$ and reversibility factor logProbChangeMod.
**BassDeathMove**: Selects an active basis uniformly at random; calculates reverse birth probability.
**BassChangeMove**: Perturbs a knot location, sign, or variable in an existing term.
Target / Marginal Likelihood (BassTarget : IRjTarget<BassState>):

Evaluates $X^T X$ and $X^T y$.
Computes quadratic form $q_f = y^T X (X^T X)^{-1} X^T y$ with Zellner's $g$-prior or unit-information prior, plus the Poisson prior on basis count $M \sim \text{Poisson}(\lambda)$.
Closed-Form Sobol Decomposition (SobolDecomposition):

Closed-form 1D polynomial overlap integrals: $$\int_0^1 [s_1(x - t_1)]+^q [s_2(x - t_2)]+^q , dx$$
Analytical expansion to compute total variance $V$, main effect variances $V_j$, and interaction variances $V_{j,k}$, outputting exact posterior distributions of Sobol indices ($S_j, T_j$).
Strict Core API & DTOs:

**BassConfig**: Declarative record (degree, max interaction, priors $g_1, g_2, h_1, h_2$, tempering ladder, chains, seed).
**BassResult**: Posterior mean surface, credible intervals, active variable frequencies, Sobol indices, and convergence diagnostics ($\hat{R}$ and ESS).
**Bass.Run(BassConfig config, double[,] x, double[] y, ...)**: Strict-core entry point.

5. Implementation Effort & Roadmap
Phase	Scope	Estimated Effort
**Phase 1:** Basis & State	BassTerm, BassState, design evaluation $X$, Gram matrix and bordered Cholesky updates	~1.5 days
**Phase 2:** RJMCMC Moves & Target	BassTarget, BassBirthMove, BassDeathMove, BassChangeMove, hyperparameter priors	~2 days
**Phase 3:** Ensemble Shell Integration	BassEnsemble : IEnsembleModel<BassState>, BassConfig, BassResult, wiring into ChainEnsemble & ParallelTempering	~1 day
**Phase 4:** Analytical Sobol Decomposition	Exact 1D hinge integral table, ANOVA tensor contraction, $S_j$ / $T_j$ posterior summaries	~1.5 days
**Phase 5:** Functional Outputs & Categorical	Categorical subset factors, functional response tensor basis, and BassPca via Maths.LinAlg.PCA	~1.5 days
**Phase 6:** Oracle Parity Tests	Validation tests in tests/regression asserting parameter recovery and numerical parity against LANL's R package 


6. Architectural Discipline
When implementing BASS in ThermoMapper, preserve the invariants in AGENTS.md and CONTRIBUTING.md:

## Layer ordering:
Keep it strictly in Maths.Regression, depending only down on Maths.Samplers, Maths.Distributions, Maths.LinAlg, and Maths.Rng.

## Strict core, fluent shell:
Keep BassConfig as a pure, serializable DTO. Live delegates or interactive REPL helpers belong at the CLI/REPL boundary.

## Oracle validation:
Use the existing R setup in r/oracles to validate C# BASS outputs against LANL's packages/BASS reference implementation.

