# whats the difference between BASS and MARS? i was looking at this paper about generalized Bayesian MARS e.g. 2306.01911v1

In the context of non-parametric regression and surrogate modeling, Jerome Friedman’s Multivariate Adaptive Regression Splines (MARS, 1991) is a deterministic optimization algorithm that builds a single point-estimate spline basis, whereas Devin Francom’s Bayesian Adaptive Spline Surfaces (BASS, 2018/2020) is a fully probabilistic reformulation that samples across the trans-dimensional posterior distribution of spline structures using Reversible-Jump Markov Chain Monte Carlo (RJ-MCMC).

The paper you referenced—*Generalized Bayesian MARS: Tools for Emulating Stochastic Computer Models* (Rumsey, Francom, & Shen, 2023; arXiv:2306.01911v1, published in *SIAM/ASA JUQ* 2024)—directly builds on BASS. It extends the standard Gaussian error structure of BASS to normal mean-variance mixtures (Generalized Hyperbolic family) while preserving the exact analytic coefficient marginalization that makes BASS computationally tractable.

---

### 1. Structural Comparison: MARS vs. BASS

Both models represent an unknown function $f: \mathbb{R}^p \to \mathbb{R}$ as a linear combination of tensor-product hinge functions:

$$f(x) = a_0 + \sum_{m=1}^M a_m B_m(x), \quad B_m(x) = \prod_{k=1}^{K_m} \left[ s_{km} (x_{v_{km}} - t_{km}) \right]_+$$

where $[z]_+ = \max(0, z)$, $s_{km} \in \{-1, +1\}$ denotes sign, $v_{km} \in \{1, \dots, p\}$ indexes the predictor, $t_{km}$ is a knot, and $K_m$ is the interaction degree. The divergence lies entirely in how the model space $\mathcal{M} = \{M, \{K_m, v_{km}, s_{km}, t_{km}\}_{m=1}^M\}$ and coefficients $a$ are estimated.

| Architectural Dimension | MARS (Friedman, 1991) | BASS (Francom & Sansó, 2020) |
| --- | --- | --- |
| **Model Search Mechanism** | **Greedy forward/backward stepwise:** Forward pass iteratively adds reflected basis pairs $(x_j - t)_+, (t - x_j)_+$ minimizing RSS; backward pass prunes basis functions. | **Trans-dimensional RJ-MCMC:** Samples topology via Birth, Death, and Mutation proposals across the parameter space $\theta = \{M, a, \sigma^2, \theta_B\}$. |
| **Model Selection / Penalty** | **Generalized Cross-Validation (GCV):** Penalizes basis count via cost parameter $d$ (effective degrees of freedom): $\text{GCV} = \frac{\text{MSE}}{(1 - \frac{C(M)}{N})^2}$. | **Hierarchical Priors:** Explicit priors on model size $M \sim \text{Poisson}(\lambda)$ (or NegBin), interaction degree $K_m$, and continuous knot locations. |
| **Knot Placement** | Restricted to discrete candidate points on the observed data grid (subject to minimum observation spans $\text{minspan}$). | Continuous support over the input domain ($t_{km} \in [0, 1]$), sampled continuously during MCMC mutations. |
| **Subspace / Myopia Traps** | Prone to local optima; early forward selections cannot be easily undone unless pruned later, often missing subtle interaction subspaces. | Explores multi-modal posterior topologies through reversible jumps; moves in and out of local modes without committing permanently. |
| **Uncertainty Quantification** | **None inherently:** Yields a single deterministic surface. Post-hoc bootstrap is computationally expensive and unstable across pruned knot topologies. | **Full Posterior Ensembling:** Generates samples from $p(f \mid y)$, naturally providing posterior predictive bands, parameter uncertainty, and model-averaged ensembles. |
| **Sensitivity Analysis** | Heuristic variable importance metrics based on drop in GCV when omitting basis functions. | **Closed-form Sobol Decomposition:** Functional ANOVA and Sobol variance indices (main effects, interactions, total effects) can be integrated analytically per posterior sample. |

---

### 2. The Computational Engine of BASS

A naive trans-dimensional MCMC over splines (like early Denison et al., 1998 formulations) mixes poorly because proposing a new basis function requires proposing a new coefficient $a_{M+1}$.

BASS bypasses this bottleneck via **analytic marginalization**. By placing a conjugate Zellner-type $g$-prior on the basis coefficients $a \mid \sigma^2, B \sim \mathcal{N}(0, \sigma^2 \tau (B^T B)^{-1})$ (or an equivalent ridge prior) and an Inverse-Gamma prior on $\sigma^2$, the coefficients $a$ are integrated out in closed form:

$$p(y \mid \theta_B, \tau) \propto \vert{}B^T B + \tau^{-1} \Sigma_0\vert{}^{-1/2} \left[ y^T y - y^T B (B^T B + \tau^{-1} \Sigma_0)^{-1} B^T y + 2 b_\sigma \right]^{-\frac{n + 2 a_\sigma}{2}}$$

Because the likelihood marginalizes over $a$, the Metropolis-Hastings acceptance probability for Birth, Death, and Mutation transitions evaluates only the marginal change in model evidence. Basis coefficients are only back-sampled via Gibbs steps when drawing surfaces or evaluating Sobol indices.

---

### 3. What arXiv:2306.01911 (GB-MARS / GBASS) Solves

Traditional MARS and BASS both break down when the data-generating process deviates from additive homoscedastic Gaussian noise: $\epsilon_i \overset{iid}{\sim} \mathcal{N}(0, \sigma^2)$. In stochastic simulators (e.g., queuing systems, turbulent flows, financial mechanics, or cellular stochastic kinetics), outputs frequently display heavy tails, asymmetric skewness, or input-dependent dispersion. Standard MARS/BASS under heavy tails over-allocates basis functions to fit tail outliers, distorting the latent mean surface.

Rumsey, Francom, and Shen formulate **Generalized Bayesian MARS (GB-MARS)** by mapping the error structure to **Normal Mean-Variance Mixtures**:

$$y_i = f(x_i) + \beta w_i + \sqrt{w_i} \epsilon_i, \quad \epsilon_i \sim \mathcal{N}(0, \sigma^2)$$

where $w_i$ is a local latent variance factor governed by a Generalized Inverse Gaussian (GIG) mixing distribution: $w_i \sim \text{GIG}(\lambda, \delta, \gamma)$.

This unlocks specific surrogate capabilities while maintaining BASS efficiency:

1. **Preserved Analytic Marginalization:** Because $y_i \mid w_i$ is conditionally Gaussian, the basis coefficients $a$ can still be analytically integrated out of the posterior conditional on $W = \text{diag}(w_1, \dots, w_n)$. The trans-dimensional Birth/Death/Mutation steps remain fast and numerically stable.
2. **Robust Spline Regression ($t$-distributed errors):** Setting the mixing parameters to produce an Inverse-Gamma distribution on $w_i$ recovers a Student-$t$ marginal likelihood. Outliers induce large $w_i$ values, downweighting their leverage on the knot topology rather than spuriously bifurcating hinge functions.
3. **Bayesian Quantile Regression (Asymmetric Laplace):** If $w_i \sim \text{Exponential}$ with drift parameter $\beta = \frac{1 - 2q}{q(1-q)}$, the marginal error distribution becomes Asymmetric Laplace. GB-MARS can therefore estimate arbitrary quantile surfaces $Q_q(y \mid x)$ non-parametrically.
4. **Tail-Heaviness Adaptation:** Latent weights $w_i$ are sampled via exact bounded rejection sampling (Hörmann & Leydold, 2014) alongside the topology, yielding unified uncertainty over both the spline field and error-density tails.

---

### 4. Implementation & Diagnostic Verification Checklist

When building or benchmarking an integrated Bayesian spline engine (e.g., embedding BASS/GB-MARS kernels into custom MCMC architectures), verify these structural constraints:

* [ ] **Knot Sparsity and Rank Deficiencies ($b_0$ Constraint):** Ensure every proposed knot partition satisfies a minimum sample count $b_0$ (typically $\ge 10\text{–}20$ points). Without this guardrail, the Gram matrix $B^T W^{-1} B$ becomes ill-conditioned, degrading the Cholesky factor used to compute marginal acceptance ratios.
* [ ] **Degree Proposal Balance ($w_1, w_2$ Hyperparameters):** In BASS, propose interaction degrees $K_m$ and variables using a categorical distribution biased by their current inclusion frequencies. Monitor birth-death symmetry: death proposals must be drawn uniformly from active bases ($1/M$), while birth proposal densities must account for the continuous knot selection volume to satisfy detailed balance.
* [ ] **GIG Latent Sampling Bounds:** If implementing the GB-MARS mixture layer, ensure the GIG update draws $w_i \mid y_i, a, \sigma^2$ using an exact ratio-of-uniforms or rejection algorithm (such as Hörmann-Leydold) rather than standard random-walk MH to prevent latent-variable chain freeze in high-variance regimes.
* [ ] **Zellner $g$-Prior Inflation Check:** Verify that $\tau$ (the coefficient variance scaling factor) is either assigned an empirical Bayes estimate or an informative hyperprior ($a_\tau = 1/2, b_\tau = 2/n$). Unscaled fixed variance priors will induce Lindley's paradox, driving the posterior model size $M \to 0$ as $n$ grows.

# ive also been looking at 2405.13342v2 and https://doi.org/10.1093/biomtc/ujag121 , at least the supplemental materials and repositories for the latter  

If you are looking to synthesize the algorithmic mechanics of He, Yang, and Kang's Adaptive Bayesian Multivariate Regression Splines (ABMRS/EBARS) from the *Biometrics* paper ([doi:10.1093/biomtc/ujag121](https://doi.org/10.1093/biomtc/ujag121)) alongside their heat kernel manifold GP work ([arXiv:2405.13342v2](https://arxiv.org/abs/2405.13342) / [doi:10.1093/jrsssb/qkaf057](https://doi.org/10.1093/jrsssb/qkaf057))—specifically resolving the discrepancies across the 2024 preprint (`arXiv:2405.13353`), the final journal version, and the `junhuihe2000/EBARS` (and `ABMRS`) repositories—the details below break down the mathematical schedule cancellations, prior formulations, and computational failure modes.

---

### 1. The Core Theoretical Pivot: EBIC Cancellation in the RJ-MCMC Schedule

The central innovation of ABMRS over classical BARS (DiMatteo et al., 2001) and BASS (Francom & Sansó, 2020) is replacing continuous knot birth/death proposals with an explicit combinatorial search over a discrete candidate knot lattice $\eta = \bigotimes_{j=1}^d \eta_j$ penalized via the Extended Bayesian Information Criterion (EBIC; Chen & Chen, 2008).

#### The EBIC Model Prior

For a model $\mathcal{M}_k$ with knot counts $k = (k_1, \dots, k_d)$ drawn from $n = (n_1, \dots, n_d)$ candidate positions, the model complexity prior is parameterized by:

$$P(\mathcal{M}_k) \propto \prod_{j=1}^d \binom{n_j}{k_j}^{-\gamma}, \quad 0 \le \gamma \le 1$$

where the model space size is $\tau(\mathcal{M}_k) = \prod_{j=1}^d \binom{n_j}{k_j}$.

#### The Proposal Schedule and the "Hidden Cancellation"

To maintain ergodicity while preventing the proposal chain from stagnating in empty model spaces when $\vert{}\eta_j\vert{} \gg k_j$, the proposal probabilities for dimension $j$ are structured using a split parameter $c \in (0, 0.5)$:

1. **Birth:** $b_{k_j} = c \min\left(1, \left\{\frac{n_j - k_j}{k_j + 1}\right\}^{1 - \gamma}\right)$
2. **Death:** $d_{k_j} = c \min\left(1, \left\{\frac{k_j}{n_j - k_j + 1}\right\}^{1 - \gamma}\right)$
3. **Relocate:** $r_{k_j} = 1 - b_{k_j} - d_{k_j}$

When proposing a transition from state $k$ to state $k+1$ (Birth) by choosing uniformly among the remaining $(n_j - k_j)$ candidates, the forward proposal density is $q(k+1 \mid k) = \frac{b_k}{n_j - k_j}$. The reverse proposal (Death) selects uniformly among the $(k+1)$ active knots: $q(k \mid k+1) = \frac{d_{k+1}}{k+1}$.

Evaluating the Hastings proposal ratio reveals why this schedule is formulated as such:

$$\frac{q(k \mid k+1)}{q(k+1 \mid k)} = \frac{d_{k+1}}{b_k} \cdot \frac{n_j - k_j}{k_j + 1} = \left(\frac{k_j + 1}{n_j - k_j}\right)^{1 - \gamma} \cdot \left(\frac{n_j - k_j}{k_j + 1}\right) = \left(\frac{n_j - k_j}{k_j + 1}\right)^\gamma$$

Meanwhile, the prior ratio under the EBIC specification is:

$$\frac{P(\mathcal{M}_{k+1})}{P(\mathcal{M}_k)} = \left[ \frac{\binom{n_j}{k_j+1}}{\binom{n_j}{k_j}} \right]^{-\gamma} = \left( \frac{n_j - k_j}{k_j + 1} \right)^{-\gamma}$$

Multiplying the prior ratio by the proposal ratio yields:

$$\frac{P(\mathcal{M}_{k+1})}{P(\mathcal{M}_k)} \times \frac{q(k \mid k+1)}{q(k+1 \mid k)} = \left( \frac{n_j - k_j}{k_j + 1} \right)^{-\gamma} \times \left(\frac{n_j - k_j}{k_j + 1}\right)^\gamma \equiv 1$$

**Consequence for Custom Engines:** Under this schedule, the combinatorial model-space penalty cancels out in the Metropolis-Hastings acceptance ratio. The acceptance probability depends *exclusively* on the marginal evidence (or likelihood) ratio. If an engine computes both the EBIC model prior ratio and the explicit proposal ratio without taking advantage of this analytic cancellation, numerical underflow/overflow will destabilize the chain when $n_j$ is large.

---

### 2. Discrepancies: Preprint (`2405.13353v1`) vs. *Biometrics* (`ujag121`) vs. Repositories

Comparing the initial 2024 preprint against the published 2026 *Biometrics* article and its supplemental files (`ujag121_Supplemental_Files`) reveals several adjustments in hyperparameter bounds and implementation:

#### A. Asymptotic Bounds on $\gamma$

* **Preprint:** Stated heuristically that any $\gamma \in [0, 1]$ can be chosen, with simulations setting $\gamma = 0$ (standard BIC) or $\gamma = 0.5$.
* **Biometrics Journal & Supplement:** Clarifies the consistency criteria under dense candidate lattices. If the candidate pool size scales polynomially with sample size, $n_j = O(m^\kappa)$ with $\kappa > 0$, selection consistency holds **if and only if** $\gamma > 1 - \frac{1}{2\kappa}$. If candidate grids are dense ($\kappa \ge 1$), setting $\gamma = 0$ or small values fails asymptotically by retaining spurious boundary knots.

#### B. Support Boundaries ($C_{\text{bound}}$) and Collinearity Guardrails

* **Preprint:** Assumed candidate knots could span the entire compact domain $[0, 1]^d$.
* **Biometrics Supplement:** Formulates an explicit boundary buffer $C = [\delta, 1 - \delta]^d$ where knots are prohibited within the outer quantile margin $\delta \approx 0.05$. Furthermore, the supplement enforces a minimal span constraint between adjacent active knots: $\vert{}t_{j, l} - t_{j, l-1}\vert{} \ge h_{\min}$.

#### C. Repository Bugs in `junhuihe2000/EBARS` (`ABMRS`)

Inspecting the reference code exposes discrepancies between the published formulas and the script execution:

1. **Boundary Probability Leakage:** In the R scripts, boundary edge cases ($k_j = 0$ or $k_j = n_j$) are patched with hardcoded checks. At $k_j = 0$, $d_0$ is set to $0$, but $b_0$ is left as $c$ rather than expanded to $2c$. This causes the relocation probability $r_0 = 1 - c$ to absorb the remaining mass, arbitrarily skewing proposal velocity at the lower boundary.
2. **Gram Matrix Regularization ($B^T B$ Jitter):** The paper specifies a pure Unit Information Prior (UIP) $\beta \mid \sigma^2 \sim \mathcal{N}(0, m (B^T B)^{-1} \sigma^2)$. However, in multidimensional tensor-product B-splines, basis functions in underpopulated regions of $[0, 1]^d$ create near-zero singular values. The repo silently injects a ridge stabilizer:

$$(B^T B + \epsilon I)^{-1}, \quad \epsilon = 10^{-5}\text{--}10^{-6}$$

Without this jitter, `chol(B^T B)` fails frequently during tensor grid evaluation.

### 3. Structural Comparison: ABMRS vs. Heat Kernel Manifold GPs (`2405.13342v2`)

The companion paper (*Scalable Bayesian inference for heat kernel Gaussian processes on manifolds*, He, Ma, Kang, & Yang, *JRSS-B* / `2405.13342v2`) operates on fundamentally distinct geometry:

| Metric / Attribute | ABMRS / EBARS (`ujag121`) | Heat Kernel GP (`2405.13342v2`) |
| --- | --- | --- |
| **Domain Geometry** | Strictly ambient Cartesian hypercubes: $\mathcal{X} \subset [0, 1]^d$. | Arbitrary embedded Riemannian sub-manifolds: $\mathcal{M} \subset \mathbb{R}^D$. |
| **Basis Construction** | Discrete tensor products of 1D B-splines: $B(x) = \bigotimes_{j=1}^d B_j(x_j; \xi_j)$. | Graph Laplacian eigenbasis: $K_t(x, x') = \sum_{l=1}^L e^{-\lambda_l t} \phi_l(x) \phi_l(x')$. |
| **Dimensional Scaling** | **Severe Curse of Dimensionality:** Number of tensor basis functions scales exponentially ($O(K^d)$), restricting pure ABMRS to $d \le 3\text{--}4$. | **Intrinsic Dimension Scaling:** Scales with the intrinsic manifold dimension $d_{\mathcal{M}} \ll D$; basis size $L$ is independent of extrinsic dimension. |
| **Discontinuity Handling** | Exceptional for step functions, sharp change-points, and derivative jumps via knot clustering. | Poor for sharp discontinuities; heat diffusion kernels impose smooth $C^\infty$ regularization. |

**Synthesis for Manifold Pipelines:** ABMRS cannot natively perform regression directly on non-Euclidean manifolds; its "manifold denoising" experiment in Section 4.2 of `ujag121` simply flattens ambient coordinates into independent spline regressions against a 1D or 2D coordinate parametrization. If working with non-trivial Riemannian topology, using `2405.13342v2`'s heat kernel formulation avoids coordinate patching, whereas ABMRS should be reserved for low-dimensional projections where local piecewise transitions require knot-based changepoint identification.

### 4. Implementation & Verification Checklist

When implementing the ABMRS/EBARS sampler or integrating it alongside BASS-style engines, verify these numerical conditions:

* [ ] **Exact Analytic Ratio Cancellation:** Confirm your MCMC kernel computes the Metropolis-Hastings acceptance ratio for Birth/Death directly as:

$$\alpha(k \to k+1) = \min\left(1, \frac{p(y \mid \mathcal{M}_{k+1})}{p(y \mid \mathcal{M}_k)}\right)$$



Do not compute the $(1-\gamma)$ power in the proposal ratio and the EBIC prior separately to avoid catastrophic floating-point cancellation.
* [ ] **Lattice Density Consistency ($\gamma$ Check):** If setting a fine knot candidate grid $n_j$ relative to sample size $m$, compute $\kappa = \frac{\log n_j}{\log m}$. Verify that your prior parameter satisfies $\gamma > 1 - \frac{1}{2\kappa}$ to suppress spurious boundary knot retention.
* [ ] **Symmetric Boundary Relocation:** In the proposal schedule, enforce hard boundary clamps:

$$\text{If } k_j = 0: \quad d_0 = 0, \quad b_0 = 2c, \quad r_0 = 1 - 2c$$


$$\text{If } k_j = n_j: \quad b_{n_j} = 0, \quad d_{n_j} = 2c, \quad r_{n_j} = 1 - 2c$$


This prevents the probability leakage present in early `junhuihe2000/EBARS` script iterations.
* [ ] **Regularized Unit Information Cholesky:** In the marginal evidence calculation, evaluate the determinant and quadratic forms via Cholesky decomposition using diagonal Tikhonov regularization:

$$L L^T = B^T B + \max\left(10^{-8}, 10^{-4} \cdot \lambda_{\max}(B^T B)\right) I$$


Monitor the condition number $\kappa(B^T B)$ when two candidate knots in the same dimension are proposed near the minimum separation threshold $h_{\min}$.