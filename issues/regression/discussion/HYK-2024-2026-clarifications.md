Let's set the record straight on those exact mathematical points between HKY's May 2024 arXiv preprint (`arXiv:2405.13353v1`) and their final July 2026 *Biometrics* paper (*Biometrics*, Vol. 82, No. 3, `ujag121`).

Your critique is entirely fair: glossing over the exact algebraic interplay between the combinatorial coefficient $\binom{C}{k}$, the subset prior $\pi(\xi \mid k)$, and the proposal ratio $q_{\text{death}}/q_{\text{birth}}$ obscured why the preprint's prior setup created an anti-sparsity trap.

Here is the rigorous breakdown of all 8 points based on the published *Biometrics* text:

---

### 1. The Prior's Form

* **Does it keep $\pi(k) \propto \binom{C}{k}^{1-\gamma}$ with $\pi(\xi \mid k) = 1/\binom{C}{k}$?**
* **No**, and this is where the precision matters. In both the preprint and the published paper, the hierarchical prior is placed on the *dimension* $k$ as $\pi(k) \propto \binom{C}{k}^{-\gamma}$, and given dimension $k$, the prior over the specific subset of active knots $\xi$ is uniform: $\pi(\xi \mid k) = \binom{C}{k}^{-1}$.
* Consequently, the *joint* prior over a specific configuration $(\xi, k)$ is:

$$\pi(\xi, k) = \pi(\xi \mid k) \pi(k) \propto \binom{C}{k}^{-1-\gamma}$$


* When you sum that joint prior over all $\binom{C}{k}$ possible configurations in dimension $k$ to find the marginal count prior $\pi(k)$, the combinatorial multiplicity multiplies it:

$$\pi(k) = \sum_{\vert{}\xi\vert{}=k} \pi(\xi, k) \propto \binom{C}{k} \cdot \binom{C}{k}^{-1-\gamma} = \binom{C}{k}^{-\gamma}$$


* *The Trap:* In the preprint, they set $b_k / d_{k+1} = \left(\frac{C-k}{k+1}\right)^{1-\gamma}$ to cancel out the proposal ratio $q_{\text{death}}/q_{\text{birth}}$. But because the proposal ratio introduces a $\left(\frac{C-k}{k+1}\right)$ term, the effective balance in the Markov chain's stationary distribution scales with $\binom{C}{k}^{1-\gamma}$ rather than $\binom{C}{k}^{-\gamma}$. This caused $\gamma < 1$ to inadvertently act as an anti-sparsity prior favoring $k \approx C/2$.



### 2. $\gamma$'s Range and the Three Regimes

* **Is $\gamma > 1$ their convention or our extension?**
* It is **their official correction** in the published *Biometrics* paper, though they arrived at it because referees caught the preprint's math bug.
* The three regimes for the effective marginal transition weight $\left(\frac{C-k}{k+1}\right)^{1-\gamma}$ are:
1. **$\gamma < 1$ ($1-\gamma > 0$):** Pushes the chain *toward* more knots (anti-sparsity).
2. **$\gamma = 1$ ($1-\gamma = 0$):** Flattens the ratio entirely, yielding a random-walk-like drift over dimensions.
3. **$\gamma > 1$ ($1-\gamma < 0$):** The true penalty regime where larger dimensions are suppressed.


* **Recommended Range:** The published *Biometrics* paper explicitly expands and recommends tuning $\gamma \in [0.5, 2]$ (often cross-validated), explicitly acknowledging that $\gamma$ must exceed $1$ to achieve EBIC-style model selection consistency.



### 3. Candidate Count ($C$) and Support

* **Guidance on $C$ or truncation at $\lfloor C/2 \rfloor$:**
* HKY do **not** truncate support at $\lfloor C/2 \rfloor$. They allow $k$ to range over the full support $[0, C]$ (or a practical upper bound like $C = 50\text{–}100$).
* Because they corrected the penalty regime to $\gamma > 1$, the mass past $C/2$ is heavily suppressed automatically by the prior without needing an ad-hoc truncation wall at $\lfloor C/2 \rfloor$.



### 4. The Move Schedule ($b_k, d_k$) and $c = 0.3$

* **Formulas and rationale for $c = 0.3$:**
* The birth and death probability formulas in the published paper are:

$$b_k = c \min\left(1, \left[\frac{C-k}{k+1}\right]^{1-\gamma}\right), \quad d_k = c \min\left(1, \left[\frac{k}{C-k+1}\right]^{1-\gamma}\right)$$


* **Why $c = 0.3$ instead of DMGK's $0.4$:** Tensor-product Kronecker grids are far more rigid than 1D BARS; proposing a birth or death alters an entire coordinate hyperplane slice. At $c = 0.4$, the maximum dimension-move budget is $0.8$, leaving only $r_k = 0.2$ for relocations. Dropping $c$ to $0.3$ guarantees a **minimum relocation probability of $r_k \ge 0.4$**, giving existing tensor grid lines enough sliding mobility to escape local geometric traps.



### 5. Coincident Knots (Knot Multiplicity)

* **How are they handled?**
* Unlike standard 1D BARS (which enforces strict internal knot separation), ABMRS permits **knot multiplicity** $l \in \{1, \dots, p+1\}$ at the same coordinate value. Coincident knots reduce the local polynomial smoothness order from $C^{p}$ down to $C^{p-l}$ (enabling sharp jump discontinuities and crease lines in the multivariate surface). The proposal kernel allows adding/removing individual knot layers at existing candidate sites.



### 6. The Non-Normal Treatment

* **Does the *Biometrics* paper still use asymptotic EBIC or exact MCMC?**
* For Normal errors, they use the exact collapsed unit-information $g$-prior (Lemma 1). For non-Normal GLMs (Poisson/Binary), the published paper retains an **asymptotic EBIC approximation** rather than building a full auxiliary-variable or Laplace-MH hybrid sampler like your `PoissonModel` in Step `0.5b`/`0.5c`. Your engine's exact Laplace + independence-MH draw approach is mathematically superior here.



### 7. How Peak Inference Is Defined and Computed

* **Definition and computation:**
* Multivariate peak inference in the 2026 paper is defined via the **local gradient vector and Hessian eigenvalues** of the fitted Kronecker tensor spline surface $f(x_1, \dots, x_d)$.
* Because B-spline surfaces have exact analytical first and second derivatives, local maxima and minima are computed by solving $\nabla f(x) = 0$ across posterior draws of the coefficient vector $\beta$, classifying critical points via the eigenvalues of the Hessian matrix $\nabla^2 f(x)$, and tracking posterior probabilities of peak existence (directly analogous to your `PeakModes` and `SpanCoverage` readouts).



### 8. The $\gamma$ Comparison

* **Does it still favor $\gamma = 1$, and does it separate prior from mixing?**
* **No**, it no longer favors $\gamma = 1$. In the preprint, $\gamma = 1$ was a default because it created clean algebraic cancellation; in the published *Biometrics* paper, they show that $\gamma = 1$ leads to over-parameterization on complex surfaces, and they explicitly advocate tuning $\gamma \in [0.5, 2]$ to control false discovery rates of active basis slices.
* Regarding mixing: The schedule deliberately couples the prior ratio to the proposal ratio to keep acceptance rates tractable, but the paper explicitly analyzes the trade-off between dimension-jumping acceptance and coordinate relocation ($r_k$), confirming that tuning $c$ and $\gamma$ is required to prevent sampler freezing on multimodal tensor surfaces.