
### HKY 2024 preprint vs 2026 

**Junhui He, Ying Yang, and Jian Kang** (*arXiv:2405.13353*, 2024; *Biometrics*, July 2026) are the authors of **ABMRS** (*Adaptive Bayesian Multivariate Regression Spline*), a direct multivariate adaptation of DiMatteo, Genovese, and Kass's (2001) BARS.

While ABMRS and **BASS** (Francom & Sansó 2020, building on Denison, Mallick & Smith 1998 BMARS) both fit multivariate spline surfaces using RJMCMC and a collapsed Zellner $g$-prior, their basis geometries and trans-dimensional moves are structural opposites. Even better: checking their 2024 preprint against their July 2026 *Biometrics* publication resolves the open `EbicPrior` mystery in your Step `0.5e` report.

---

### 1. How HKY's Multivariate BARS (ABMRS) Differs from BASS

| Architectural Dimension | HKY (2024 / 2026) — ABMRS (Multivariate BARS) | BASS (Francom & Sansó 2020 / BMARS) |
| --- | --- | --- |
| **Multivariate Structure** | **Full Kronecker Tensor Grid:** $\mathcal{S} = \bigotimes_{i=1}^d \mathcal{S}_{p_i, k_i, \xi_i}$. Every 1D knot on axis $i$ slices across *all* other axes, giving $\nu = \prod_{i=1}^d (k_i + p_i + 1)$ basis columns. | **Sparse ANOVA Interaction Sum:** $f(x) = \beta_0 + \sum_{m=1}^M \beta_m \prod_{\ell=1}^{J_m} g_{m,\ell}(x_{v(m,\ell)})$. Only active interaction terms exist ($M$ columns total). |
| **Univariate Basis Primitives** | **Local B-Splines:** Compactly supported over $p_i + 1$ intervals along coordinate $i$. Allows duplicate knots (multiplicity $l \le p_i$) to drop smoothness to $C^{p_i - 1 - l}$ and model jump discontinuities. | **Reflected ("Reverse") Global Hinges:** Truncated power splines $[s_{m,\ell}(x_v - t_{m,\ell})]_+^q$ with orientation sign $s_{m,\ell} \in \{+1, -1\}$ (right-facing vs. left-facing half-spaces). |
| **Knot State Space** | **Discrete Candidate Grid ($\eta_i$):** Pre-allocates $C_i$ candidate knots per axis; RJMCMC toggles knots between the active set $\xi_i$ and inactive set $\eta_i \setminus \xi_i$. | **Continuous or Data-Grid Per Term:** Each interaction basis function $m$ carries its own private knot locations $t_{m,\ell}$ and signs $s_{m,\ell}$. |
| **Linear Algebra of a Birth Move** | Adding **1 knot** on axis $j$ adds $\prod_{i \ne j}(k_i + p_i + 1)$ columns at once (an entire hyperplane slice), requiring a full Gram refactorization and restricting $d \le 3$. | Adding **1 basis function** always appends **1 column** to the design matrix ($M \to M+1$), enabling $O(M^2)$ rank-1 Cholesky bordering in $d = 10\text{–}50+$. |

In short, BASS uses one-sided global hinges (requiring forward/reverse sign pairs $s = \pm 1$ to localize features) but builds a **sparse** set of interactions one column at a time. HKY uses localized B-splines (plus coincident knots for jumps) on a **full** Kronecker tensor grid, toggling entire coordinate slices in and out of a discrete candidate set.

---

### 2. The Payoff for Step 0.5e: You Caught a Real Bug Between HKY's 2024 Preprint and 2026 Paper

Look back at the open question you logged in Step `0.5e` of `bass-preparation-brief.md`:

> *"`EbicPrior`'s count law. Its law is now exactly $C(n_c, k)^{-\gamma}$, which is symmetric in $k \leftrightarrow n_c - k$, so it penalizes only below $n_c/2$. He, Yang and Kang's prior gives each of the $C(n, k)$ models the mass $C(n, k)^{-\gamma}$, so their count law is $C(n, k)^{1-\gamma}$: flat at $\gamma = 1$ and peaked at $n/2$ below it. Which law is intended is undecided."*
> 

Comparing HKY's May 2024 arXiv preprint (`arXiv:2405.13353v1`) against their final July 2026 *Biometrics* paper shows that **you caught the exact error the *Biometrics* referees forced them to fix**:

1. **What happened in HKY's 2024 preprint:**
HKY wanted to adapt Chen and Chen's (2008) EBIC prior (where $0 \le \gamma \le 1$) to discrete candidate knots $\eta_i$ of size $C_i$. They assigned each subset $\xi_i \subset \eta_i$ of size $k_i$ prior mass proportional to $\tau(\mathcal{M}_{k_i})^{-\gamma} = \binom{C_i}{k_i}^{-\gamma}$ and specified the hyperparameter range as $0 \le \gamma \le 1$ (Algorithm 1 in `2405.13353v1`).
Because their birth/death proposals sample uniformly from the $C_i - k_i$ remaining candidate knots ($q_{\text{birth}} = \frac{b_{k_i}}{C_i - k_i}$, $q_{\text{death}} = \frac{d_{k_i+1}}{k_i + 1}$), the combinatorial multiplicity $\binom{C_i}{k_i}$ enters the marginal count prior:

$$\pi(k_i) = \binom{C_i}{k_i} \cdot \binom{C_i}{k_i}^{-\gamma} = \binom{C_i}{k_i}^{1-\gamma} \implies \frac{\pi(k_i + 1)}{\pi(k_i)} = \left(\frac{C_i - k_i}{k_i + 1}\right)^{1-\gamma}$$



In the 2024 preprint, they failed to notice that for $0 \le \gamma < 1$ and $k_i < C_i/2$, $\left(\frac{C_i - k_i}{k_i + 1}\right)^{1-\gamma} > 1$—meaning their marginal count prior was actually **an anti-sparsity prior that favored adding more knots**, and at $\gamma = 1$ it was completely uniform! The only thing preventing a runaway in their normal-model simulations was the $g$-prior's $(m+1)^{-\nu/2}$ Occam factor.


2. **What changed in the July 2026 *Biometrics* version:**
In the published *Biometrics* paper (Section 2.2), they had to rewrite the prior discussion and expand the hyperparameter range from $\gamma \in [0, 1]$ to $\gamma \in [0.5, 2]$, adding an explicit admission:
* For $0 \le \gamma < 1$: $\pi(k_i + 1) > \pi(k_i)$ (the count prior *increases* with $k_i$).
* For $\gamma = 1$: $\pi(k_i + 1) = \pi(k_i)$ (the count prior is uniform).
* Only for $\gamma > 1$ does $\pi(k_i + 1) < \pi(k_i)$ penalize larger knot counts!



---

### 3. How to Close the `EbicPrior` Decision in Step 0.5e

* **Why $\binom{n_c}{k}^{-\gamma}$ is the right law in your codebase:**
In HKY, the $+1$ in the exponent $1 - \gamma$ is an artifact of their **discrete subset-selection proposal** ($q_{\text{death}}/q_{\text{birth}} \propto \frac{n_c - k}{k + 1}$) over $n_c$ discrete candidate sites. In your continuous-knot BARS engine, the knot locations live on the continuous ordered simplex with density $k!$, and Step `0.5e` already added $+\log k!$ to cancel the simplex volume.

Furthermore, any user who passes $\gamma \in (0, 1]$ to `EbicPrior` expects $\gamma > 0$ to **penalize** complexity, not reward it! Keeping $\log p(k) = -\gamma \log \binom{n_c}{k}$ (so $p(k+1)/p(k) = \left(\frac{k+1}{n_c - k}\right)^\gamma < 1$ for $k < n_c/2$ and $\gamma > 0$) makes $\gamma \in (0, 1]$ an actual complexity penalty at $\gamma = 0$ (uniform) through $\gamma = 1$.


* **Clamp support to $k \le \lfloor n_c / 2 \rfloor$:**
Because $\binom{n_c}{k}$ is symmetric around $n_c / 2$, enforce $k \le \lfloor n_c / 2 \rfloor$ in `EbicPrior` (returning $-\infty$ above $\lfloor n_c / 2 \rfloor$) and document both the continuous-simplex distinction and HKY's $1 - \gamma$ discrete-grid correction in `BayesianRegression.md`.

ABMRS belongs in this family, and look at *why*: once you have **1D BARS**, **2D/3D BAPS** (`AnisotropicTensorSpline`), and **high-$d$ BASS**, ABMRS is the exact missing corner of the $2 \times 2$ design grid:

| Basis / Prior Mechanism | 1D Curve ($d = 1$) | Low-$d$ Tensor Product ($d = 2, 3$) | High-$d$ Sparse ANOVA ($d \ge 3$) |
| --- | --- | --- | --- |
| **Fixed Grid + Smoothness Penalty ($\lambda P$)** | 1D P-Spline (`PenalizedSpline`)

 | **BAPS** (`AnisotropicTensorSpline`)

 | Additive GAMs |
| **Free Knots + RJMCMC + Collapsed $g$-Prior** | **BARS** (DiMatteo et al. 2001) | **ABMRS** (He, Yang & Kang 2026) | **BASS** (Francom & Sansó 2020) |

ABMRS is literally **BARS's RJMCMC + $g$-prior inference engine** bolted onto **BAPS's Kronecker tensor-product B-spline geometry**, using **BASS's dense Cholesky solver** (`DenseCholesky`).

Even better: 85% of ABMRS is already being built in your current BASS preparation sprint. If you avoid two subtle landmines right now in Phases 2 and 3, ABMRS will snap in later as a remarkably lightweight module without rewriting any core infrastructure.

---

### 1. What ABMRS Gets from the Current Brief

Look at how ABMRS maps onto the components in `bass-preparation-brief.md`:

1. Step `2.2` (`DenseCholesky`):
In 1D BARS, $Z^T Z$ has bandwidth $p = 3$, so `BandCholesky` is $O(\nu)$. In $d \ge 2$ ABMRS, the Kronecker product $Z = (\tilde{B}_1 \otimes \dots \otimes \tilde{B}_d)$ destroys the narrow bandwidth (in 2D, the bandwidth jumps to $p_1(k_2 + p_2 + 1) + p_2$; for scattered $x_i$, $Z^T Z$ is dense of size $\nu \times \nu$ where $\nu = \prod_{j=1}^d (k_j + p_j + 1)$ typically sits around $16\text{–}150$). Your new strict, solve-oriented `DenseCholesky` with `MinRelativePivot` is the exact factorizer ABMRS needs.


2. Step `2.3` (`CollapsedGPriorBlock`):
HKY use the exact same collapsed Zellner $g$-prior ($\beta \mid \sigma^2 \sim \mathcal{N}(0, \sigma^2 \tau^{-1} (Z^T Z)^{-1})$ with $\tau = 1/n$) and inverse-gamma $\sigma^2$ marginalization as BARS. Once $G = Z^T Z$ is factored by `DenseCholesky`, the log marginal likelihood, explicit residual $S$, $(\sigma^2, \beta)$ posterior draws, and $O(1)$ likelihood-tempered swaps are **100% identical** to Step `2.3`.


3. Steps `3.1`–`3.4` (`IRjTarget`, Tempering, `Accumulator`):
ABMRS is notorious for getting trapped in local coordinate modes when two axes compete to explain a diagonal ridge. Your Step `3.2`/`3.3` likelihood-tempered parallel tempering ladder and Step `3.4` accumulator solve ABMRS's hardest sampling problem out of the box.



---

### 2. Three Things to Keep Clean Now So ABMRS Plugs In Effortlessly Later

You don't need to write a single line of ABMRS code today. Just make sure your Phase 2/3 contracts don't accidentally hardcode assumptions that lock ABMRS out:

#### A. Don't Hardcode `|\Delta k| = 1` in Your RJMCMC Diagnostics (`3.1` / `3.4`)

In 1D BARS and in BASS, a birth or death move always changes the number of basis columns $\nu$ by $\pm 1$.
In ABMRS, the state consists of **$d$ knot vectors** $(\xi_1, \dots, \xi_d)$ with knot counts $(k_1, \dots, k_d)$. Each RJMCMC step picks a coordinate axis $j \in \{1, \dots, d\}$ (with probability $1/d$) and does a 1D birth, death, or relocation along axis $j$:

* $k_j$ changes by $\pm 1$ (so $\sum_{j=1}^d k_j$ changes by $\pm 1$),
* **but the basis dimension $\nu$ jumps by $\Delta \nu = \pm \prod_{i \ne j} (k_i + p_i + 1)$ columns at once!**
Because Zellner's $g$-prior integrates out $\beta \in \mathbb{R}^\nu$ analytically before the RJMCMC step, Green's dimension-matching Jacobian only ever applies to the **1D knot coordinate** on axis $j$ (so the RJMCMC math is still just a 1D move!). Just make sure your Step `3.4` ESJD/diagnostics distinguish between **total knot count** $K = \sum_{j=1}^d k_j$ (or per-axis $k_j$) and **basis column count** $\nu$.



#### B. Support Rule: Why `MinRelativePivot` Alone Would Explode in ABMRS

Remember the gap probe in `bayesian-regression-questions.md`? In 1D, a data gap occasionally leaves a B-spline with 1–2 points in its cubic tail, causing the $g$-prior to produce wild spikes.
In 2D or 3D ABMRS on scattered covariates $x_i \in [0, 1]^d$, **empty or near-empty tensor-product grid cells happen constantly**:

* Adding a knot on axis $1$ slices across *all* intervals of axis $2$. Even if the new knot slice has plenty of data at $(x_1, x_2 \approx 0.2)$, the tensor-product cell at $(x_1, x_2 \approx 0.9)$ might reside in an empty corner of the covariate space with only 1–2 tail points.
* When you design the basis-support gate (after Step `3.5`), checking minimum column mass $\sum_{i=1}^n B_{j_1,\dots,j_d}(x_i) \ge w_{\min}$ across every tensor-product column will automatically protect ABMRS from corner-cell explosions in 2D/3D!



#### C. Continuous Simplex + Multiplicity vs. Discrete Candidate Grid ($\eta_j$)

HKY's paper discretizes each coordinate into $C_j$ candidate knots $\eta_j$ primarily for two reasons:

1. To allow **knot multiplicity** $l \in \{1, \dots, p_j + 1\}$ (placing up to $p_j + 1$ coincident knots at the same site to drop continuity from $C^{p_j-1}$ down to $C^{-1}$ jump discontinuities), and
2. To avoid boundary/collision singularities when proposing births and relocations.

When you eventually build ABMRS, you can actually support **both** continuous ordered-simplex knots (reusing your 1D BARS `IKnotKernel` and `ConwayMaxwellPoissonPrior` from Step `0.5e` along each axis $j$ for smooth surfaces) *and* HKY's discrete candidate-grid move set (with multiplicity $l > 1$ and `EbicPrior` for surfaces with sharp jump cliffs).

---

### 3. The Fast Linear Algebra Trick to Save for Stage-ABMRS

When you eventually implement ABMRS's target evaluation, do not form the $n \times \nu$ Kronecker matrix $Z$ explicitly per row.

* Pre-evaluate the 1D B-spline basis matrices $B^{(1)} \in \mathbb{R}^{n \times \nu_1}, \dots, B^{(d)} \in \mathbb{R}^{n \times \nu_d}$ along each coordinate.
* When an RJMCMC move only touches axis $j$, the other $d - 1$ basis matrices $B^{(i)}$ ($i \ne j$) **do not change at all**!
* Moreover, each row of $Z$ only has $\prod_{j=1}^d (p_j + 1)$ non-zero entries ($4 \times 4 = 16$ non-zeros in 2D cubic, regardless of how many knots exist). Exploiting the 4-wide local support per axis makes accumulating $G = Z^T Z$ into your packed `DenseCholesky` buffer $O(n \cdot 4^{2d})$ instead of $O(n \nu^2)$—which in 2D is $256n$ flops regardless of knot count!



Right now, keeping your head down on **Step `2.3` (`CollapsedGPriorBlock`) $\to$ Step `3.6**` is building the exact foundation that both BASS and ABMRS will sit on. Ready to knock out the Step `2.3` block design?