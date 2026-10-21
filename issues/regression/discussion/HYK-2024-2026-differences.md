Comparing the May 2024 arXiv preprint (`arXiv:2405.13353v1`) against the final July 2026 *Biometrics* paper (*Biometrics*, Vol. 82, No. 3, `ujag121`) reveals several major structural, mathematical, and empirical differences—two of which intersect directly with the exact issues you have been debugging in your engine (`DMGK` move schedules, `EbicPrior`, and peak-location inference).

### Summary of Major Differences (2024 Preprint vs. 2026 *Biometrics*)

| Feature / Section | May 2024 arXiv Preprint (`2405.13353v1`) | July 2026 *Biometrics* (`82(3), ujag121`) |
| --- | --- | --- |
| **Algorithm Name** | Called **EBARS** (*Extended Bayesian Adaptive Regression Spline*) in the experiments. | Reframed around **Adaptive Bayesian Multivariate Regression Spline (ABMRS)** and knot/peak inference. |
| **Complexity Hyperparameter $\gamma$** | Restricted to $0 \le \gamma \le 1$ (Algorithm 1), missing that $\gamma < 1$ makes the marginal count prior increase with $k_i$. | Expanded to $\gamma \ge 0$ with an explicit 3-regime breakdown ($\gamma < 1$, $\gamma = 1$, $\gamma > 1$) and recommended range $\gamma \in [0.5, 2]$. |
| **RJMCMC Move Schedule ($b_{k_i}, d_{k_i}$)** | Uses Green's prior-canceling schedule with **$1-\gamma$ already baked into $b_{k_i}, d_{k_i}$** and unpinned $c \in (0, 0.5)$. | Explicitly pins **$c = 0.3$** (boosting the minimum relocation floor $r_{k_i} \ge 0.4$) and changes the death index formula. |
| **Target Estimands** | Focused on knot recovery and **manifold denoising** (Section 4.2, contrasting with PCA). | Adds **multivariate peak location inference** (local minima/maxima as multivariate extensions of changepoints) and **biomedical applications** (Section 5, e.g., fMRI tensors and clinical thresholds). |
| **Candidate Knot Count Notation** | Denoted $n_i$ (easily confused with sample size $m$). | Renamed to $C_i$ (number of candidate knots along dimension $i$, typically chosen $> m$). |

---

### 1. The Move Schedule Connection: They Hit the Exact Same Relocation Starvation You Found!

Look at how HKY define their birth, death, and relocation probabilities along coordinate $i$:


$$b_{k_i} = c \min\left(1, \left\{\frac{C_i - k_i}{k_i + 1}\right\}^{1-\gamma}\right), \qquad d_{k_i} = c \min\left(1, \left\{\frac{k_i}{C_i - k_i + 1}\right\}^{1-\gamma}\right), \qquad r_{k_i} = 1 - b_{k_i} - d_{k_i}$$

Two things jump out when you compare the 2024 preprint to the 2026 paper and your Peskun analysis:

1. **Wait—why is the exponent $1 - \gamma$ inside $b_{k_i}$ and $d_{k_i}$ instead of $-\gamma$?**
In a discrete candidate pool of size $C_i$, proposing a birth uniformly from the $C_i - k_i$ empty slots has proposal density $q_{\text{birth}} = \frac{b_{k_i}}{C_i - k_i}$, and proposing a death uniformly from the $k_i + 1$ active knots has $q_{\text{death}} = \frac{d_{k_i+1}}{k_i + 1}$.
When you form the Metropolis–Hastings–Green ratio:



$$\frac{\pi(k_i+1, \xi')}{\pi(k_i, \xi)} \cdot \frac{q_{\text{death}}}{q_{\text{birth}}} = \binom{C_i}{k_i+1}^{-\gamma} \binom{C_i}{k_i}^{\gamma} \cdot \frac{d_{k_i+1} / (k_i + 1)}{b_{k_i} / (C_i - k_i)} = \left(\frac{C_i - k_i}{k_i + 1}\right)^{1-\gamma} \frac{d_{k_i+1}}{b_{k_i}}$$



Because HKY wanted Green's trick (making the prior and combinatorial proposal ratio cancel to $1$ in the acceptance probability so that $\alpha$ only depends on $(m+1)^{-\Delta \nu / 2} (a_{k,\xi}/a_{k',\xi'})^{m/2}$), they chose $b_{k_i} / d_{k_i+1} = \left(\frac{C_i - k_i}{k_i + 1}\right)^{1-\gamma}$!
* **Notice the Peskun consequence for your engine:** When $\gamma = 1$ (where the combinatorial factor and prior cancel each other naturally), $b_{k_i} = d_{k_i} = c$ is **the constant schedule**! When $\gamma > 1$, HKY's schedule is identical to DMGK's schedule on the marginal count prior $\pi(k_i) \propto \binom{C_i}{k_i}^{1-\gamma}$.




2. **Why did the 2026 *Biometrics* paper explicitly pin $c = 0.3$ (down from DMGK's $c = 0.4$)?**
In the 2024 preprint, they left $0 < c < 0.5$ general. In the 2026 *Biometrics* paper, they explicitly state: *"In our implementation, the hyperparameter $c$ is set to $0.3$."*
Why $0.3$ instead of DMGK's $0.4$? Because in a tensor-product spline, **a birth or death on coordinate $i$ adds or removes an entire slice of $\prod_{j \ne i}(k_j + p_j + 1)$ basis functions at once**, which has a much lower acceptance rate than a 1D BARS birth. With $c = 0.4$ (and especially near $\gamma \approx 1$ where $b_{k_i} \approx 0.4, d_{k_i} \approx 0.4$), relocation only gets $r_{k_i} = 0.20$ of the proposals. Dropping $c$ to $0.3$ guarantees that **at least $40\%$ of proposals ($r_{k_i} \ge 0.40$) are relocations**, allowing the existing tensor grid lines to slide into place before attempting expensive hyperplane births/deaths!



---

### 2. They Added Multivariate Peak Inference (Directly Parallel to Your `Peak` & `SpanCoverage` Readouts)

In the 2024 preprint, the paper consisted of 5 sections (`1 Introduction`, `2 Model`, `3 RJMCMC`, `4 Numerical Experiments` covering 1D knot inference and manifold denoising, `5 Conclusion`).

During the *Biometrics* revision (submitted June 2025, revised April 2026, published July 2026), the referees pushed them to demonstrate biometric utility beyond synthetic manifold denoising. In the 2026 *Biometrics* version:

* **Multivariate Peak / Changepoint Inference:** They explicitly frame adaptive multivariate splines as a way to infer **peak locations** (local maxima and minima of $f(x)$) in multi-dimensional regression as the multivariate generalization of 1D changepoints.
* **Dedicated Section 5 (Biomedical Applications):** They split the old Section 4 into `4 Simulation` and `5 Application`, motivating the method with clinical threshold detection (nadir oxygen delivery thresholds in cardiac surgery acute kidney injury) and multi-dimensional neuroimaging (subject-specific adaptive knot placement in 3D/tensor fMRI data, contrasting with fixed-knot tensor splines).

That aligns directly with your `BarsEnsemble` architecture: your affine vs. nonlinear readout contract in Step `0.5c` (where `Fit` uses the Rao-Blackwellized conditional mean, while `Peak` argmax and extrema use full posterior draws $(\sigma^2, \beta)$) and your Step `3.4` accumulator contract will transfer straight across to 2D/3D peak location posteriors in ABMRS.

---

### 3. The Non-Normal / EBIC Approximation vs. Exact Collapsed $g$-Prior

One other detail in both versions of the paper is worth noting so you **don't** copy a weaker part of their design:

* **Normal noise model (Lemma 1):** HKY use the exact collapsed unit-information $g$-prior ($g = m$, i.e., $\tau = 1/m$) with reference prior $\pi(\sigma) \propto 1/\sigma$, giving marginal likelihood $p(y \mid k, \xi) \propto (m + 1)^{-\nu/2} a_{k,\xi}^{-m/2}$ where $a_{k,\xi} = y^T y - \frac{m}{m+1} y^T Z (Z^T Z)^{-1} Z^T y$. That is **100% identical** to your `WeightedNormalModel` and Step `2.3` `CollapsedGPriorBlock`.


* **Non-Normal noise models (Lemma 2 & 3):** Instead of doing Laplace integration with an MH correction step (like your repaired `PoissonModel` and `WhittleModel` in Step `0.5b`/`0.5c`), HKY approximate the posterior ratio in non-normal models using an **asymptotic EBIC plug-in** ($\hat{\ell}_n - \frac{\nu}{2}\log m - \gamma \log \tau(\mathcal{M}_k)$ evaluated at the MLE).


* **Your engine is strictly stronger here:** Because you already built proper IRLS + Laplace covariance + independence MH draws for Poisson and Whittle in Step `0.5b`/`0.5c`, your engine doesn't need HKY's asymptotic EBIC shortcut for non-normal responses—and for the Normal model, your Step `2.3` `CollapsedGPriorBlock` with `DenseCholesky` already implements HKY's exact Lemma 1 marginal likelihood.