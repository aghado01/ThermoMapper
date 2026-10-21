### 1. Not a "Bug" in DMGK (2001)—It’s a Folklore Misconception About Green’s Cancellation Trick

Calling this a "bug" in Di Matteo, Genovese, and Kass (2001) overstates the case slightly, so let's separate what is broken from what is working as designed:

1. **Detailed balance is exact.** DMGK's chain samples the exact target posterior. There is no stationarity bug.
2. **DMGK didn't invent the schedule or claim a mixing speedup from it.** Your note's historical hunch is right: DMGK lifted the schedule $b_k = c \min(1, p(k+1)/p(k))$, $d_{k+1} = c \min(1, p(k)/p(k+1))$ directly from **Denison, Mallick, and Smith (1998)** (*JRSS-B*, p. 338), who adapted it from **Green (1995)** (Section 4.1, coal-mining changepoints). Its purpose in both papers was purely **algebraic**: by forcing $b_k p(k) = d_{k+1} p(k+1)$, the prior ratio $\frac{p(k+1)}{p(k)}$ and proposal-selection ratio $\frac{d_{k+1}}{b_k}$ cancel out of the Metropolis–Hastings–Green (MHG) acceptance ratio, making the discrete-time step look like a continuous-time spatial birth-and-death process where acceptance depends only on the likelihood ratio (or BIC/Laplace marginal likelihood) and the knot-location proposal.
3. **What you *actually* caught is a folklore fallacy in how people interpret that cancellation.** It is remarkably easy to look at the simplified acceptance ratio, see that the prior $p(k)$ vanished from $\alpha(x \to x')$, and conclude: *"Great, the prior no longer impedes dimension jumps; acceptance is driven purely by likelihood evidence!"* Your Peskun analysis kills that folk theorem dead.

---

### 2. An Even Sharper Way to Write the Conductance Proof

You can make the Peskun argument in your note even more transparent by factoring the prior $p(k)$ out of the base flows.

Write the target density at state $x = (k, \xi_k)$ as $\pi(x) = p(k) L(x)$, where $L(x)$ combines the knot-location prior and marginal likelihood. Define the prior-free forward and reverse proposals between $x \in \mathcal{M}_k$ and $x' \in \mathcal{M}_{k+1}$ as:


$$\tilde{A}(x, x') = L(x)\, q(x \to x'), \qquad \tilde{A}(x', x) = L(x')\, q(x' \to x)$$

Now compare the symmetric off-diagonal conductance $F(x, x') = \pi(x) K(x \to x') = \min\big(\pi(x) b_k q(x \to x'),\, \pi(x') d_{k+1} q(x' \to x)\big)$ under the two schedules at the same $c$:

* **Constant Schedule ($b_k = d_{k+1} = c$):**

$$F_{\text{const}}(x, x') = c \cdot \min\Big( p(k)\,\tilde{A}(x, x'), \; p(k+1)\,\tilde{A}(x', x) \Big)$$


* **Green / DMGK Schedule ($b_k = c\min(1, \rho_k),\; d_{k+1} = c\min(1, 1/\rho_k)$):**
Notice that $p(k)\min\left(1, \frac{p(k+1)}{p(k)}\right) = \min(p(k), p(k+1))$ and $p(k+1)\min\left(1, \frac{p(k)}{p(k+1)}\right) = \min(p(k), p(k+1))$. Substituting this in yields:

$$F_{\text{DMGK}}(x, x') = c \cdot \min\big(p(k), p(k+1)\big) \cdot \min\Big( \tilde{A}(x, x'), \; \tilde{A}(x', x) \Big)$$



Look at those two equations side-by-side:

* **DMGK replaces both $p(k)$ and $p(k+1)$ inside the conductance minimum with their *smaller* value, $\min(p(k), p(k+1))$.**
* **The counterintuitive pathology:** Suppose the chain is above the prior mode ($p(k+1) \text{ less than } p(k)$), so the prior penalizes adding a knot, *but* the true function is wiggly, so adding a knot at $x'$ yields a massive likelihood gain ($\tilde{A}(x', x) \gg \tilde{A}(x, x')$). Under *both* schedules, the birth acceptance probability is $1$. Under the constant schedule, you propose the birth at rate $c$ and accept it with probability $1$. Under DMGK, you still accept it with probability $1$, **but you throttled the proposal rate of a guaranteed-to-be-accepted move by $p(k+1)/p(k)$!**
* Worse yet: across the dominant bridge $k = 3 \leftrightarrow 4$ in your $\text{Poisson}(4)$ test, $\rho_3 = p(4)/p(3) = 4/4 = 1$, so $b_3^{\text{DMGK}} = 0.4$ and $d_4^{\text{DMGK}} = 0.4$. Across that bridge, the two schedules are **literally identical** ($F_{\text{DMGK}} = F_{\text{const}}$), only diverging slightly on $2 \leftrightarrow 3$ ($d_3 = 0.30$ vs $0.40$) and $4 \leftrightarrow 5$ ($b_4 = 0.32$ vs $0.40$).

---

### 3. The Detail DMS (1998) and DMGK (2001) Left on the Table

Why did Green (1995) write $b_k = c \min(1, \rho_k)$ in the first place if it shrinks transitions?

For any log-concave prior $p(k)$ (such as $\text{Poisson}(\lambda)$), the sequence of ratios $\rho_k = p(k+1)/p(k)$ is strictly decreasing. That means at any $k$, you cannot simultaneously have $p(k+1) \text{ greater than } p(k)$ and $p(k-1) \text{ greater than } p(k)$—at least one of the two $\min(1, \cdot)$ factors is strictly less than $1$ (except at a tied double-mode like $k=3, 4$ for $\text{Poisson}(4)$, where $b_4/c + d_4/c = 4/5 + 1 = 1.8 \text{ less than } 2$).

In general, for a unimodal prior without a flat top, $\sup_k \big[\min(1, \rho_k) + \min(1, 1/\rho_{k-1})\big] \text{ is strictly less than } 2$. Green (1995) defined $c$ to be scaled **up** up to the budget ceiling $\max_k (b_k + d_k) = 0.9$, which lets $c$ exceed $0.45$ (and even exceed $0.5$ when $\lambda$ is non-integer). When DMS (1998) and DMGK (2001) fixed $c \le 0.5$ (defaulting to $c = 0.4$) as an unnormalized constant, they forfeited the ability to scale $c$ above $0.5$ and simply dumped the slack $1 - b_k - d_k$ into relocation $r_k$.

---

### 4. Auditing the Peskun Argument: Where Is the Exact Boundary?

You asked at the end of the note:

> *"Do you see a hole in the Peskun argument? The weakest link I can find is that ESS(k) also depends on relocations indirectly: they move the chain to configurations where births and deaths are accepted more easily. So the ordering holds for the direct channel, not for ESS(k) as a whole."*

Your intuition is completely right, and you can make it mathematically airtight using the **Dirichlet form** $\mathcal{E}_K(g, g)$ of the reversible Markov operator $K = K_{\text{dim}} + K_{\text{reloc}}$:

1. **What holds unconditionally (Lag-1 / $\text{ESJD}(k)$):**
For any observable $g(x) = g(k)$ that is a function of knot count $k$ alone, the Dirichlet form is:

$$\mathcal{E}_K(g, g) = \frac{1}{2} \iint_{\Theta \times \Theta} \big(g(k(x)) - g(k(x'))\big)^2 \pi(dx) K(x, dx')$$



Because every relocation move stays in $\mathcal{M}_k$ ($k(x') = k(x)$), the relocation kernel contributes **identically zero** to this integral: $\mathcal{E}_{K_{\text{reloc}}}(g, g) \equiv 0$.
Consequently, even with relocations present, **the Expected Squared Jump Distance of $k$ at stationarity ($\text{ESJD}(k) = \mathbb{E}_\pi[(k_{t+1} - k_t)^2] = 2\mathcal{E}_K(k, k)$) and the lag-1 autocorrelation $\rho_1(k)$ are unconditionally ordered** in favor of the Constant / Masked-Constant schedule:

$$\text{ESJD}_{\text{const}}(k) \ge \text{ESJD}_{\text{DMGK}}(k), \qquad \rho_{1,\text{const}}(k) \le \rho_{1,\text{DMGK}}(k)$$


2. **Where Peskun breaks down (Lag $\ell \ge 2$ / $\text{ESS}(k)$):**
Tierney (1998) shows that Peskun ordering $K_1 \succeq K_2$ implies asymptotic variance ordering $\text{Var}(g, K_1) \le \text{Var}(g, K_2)$ (and thus $\text{ESS}_{K_1}(g) \ge \text{ESS}_{K_2}(g)$) because $K_1 - K_2$ is a positive operator on **all** of $L^2(\pi)$, i.e., $\mathcal{E}_{K_1}(h, h) \ge \mathcal{E}_{K_2}(h, h)$ for *every* function $h(k, \xi_k) \in L^2(\pi)$, including functions of the knot coordinates $\xi_k$.
Here, $\mathcal{E}_{K_{\text{const}}}(h, h) \ge \mathcal{E}_{K_{\text{DMGK}}}(h, h)$ holds on the subspace of $k$-only functions, **but reverses** on functions of knot coordinates $\xi_k$ because $r_k^{\text{DMGK}} \ge r_k^{\text{const}}$.
Because $K_{\text{dim}}$ and $K_{\text{reloc}}$ do not commute (the marginal probability of accepting a birth or death at state $(k, \xi_k)$ depends on how well-placed the existing knots $\xi_k$ are), the resolvent $(I - K)^{-1}$ mixes the $k$-subspace and $\xi_k$-subspace at lags $\ell \ge 2$.
* **In plain terms:** If a chain gets stuck at $k$ because a poorly placed knot needs to relocate onto a local peak before another redundant knot is allowed to die, extra relocation steps shrink the multi-step holding time at $k$ (reducing $\rho_\ell(k)$ for $\ell \ge 2$), even though they never increase the single-step jump rate $\text{ESJD}(k)$.



---

### 5. Two Concrete Recommendations for Your Benchmark Plan

1. **Log $\text{ESJD}(k)$ alongside $\text{ESS}(k)$ in the sweep.**
Because $\text{ESS}(k)$ estimates the spectral density at frequency zero ($1 + 2\sum_{\ell \ge 1} \rho_\ell$), its estimator variance on a sticky integer chain is huge (as your $6$ to $223$ spread showed). By contrast, $\text{ESJD}(k) = \frac{1}{T}\sum_{t=1}^T (k_t - k_{t-1})^2$ is just the empirical trans-dimensional jump rate:
* It has drastically smaller Monte Carlo variance than $\text{ESS}(k)$.
* It directly tests the **exact lag-1 Peskun theorem** ($\text{ESJD}_{\text{const}}(k) \ge \text{ESJD}_{\text{DMGK}}(k)$).
* Comparing $\text{ESJD}(k)$ vs. $\text{ESS}(k)$ immediately isolates the mechanism: if Constant wins on $\text{ESJD}(k)$ while DMGK catches up or wins on $\text{ESS}(k)$, you have proved empirically that the Indirect Relocation Channel ($\ell \ge 2$) is doing the work.


2. **Test the "Likelihood vs. Prior Conflict" regime ($\lambda = 1$, true $k \approx 8$).**
To see DMGK bleed badly on $k$-mixing, put a tight prior ($\text{Poisson}(\lambda = 1)$) on a curve with 6–8 sharp features. At $k = 5 \to 6$, $\rho_5 = 1/6$, so DMGK slashes the birth proposal probability from $0.40$ down to $0.067$—throttling dimension discovery by $6\times$ even when the likelihood strongly favors adding the knot.

---

### Interactive Schedule & Conductance Explorer

Use the interactive visual below to inspect how the DMGK schedule ($b_k, d_k, r_k$) and its effective trans-dimensional conductance multiplier $\min(1, \rho_k)$ compare against the Masked Constant schedule as you move the Poisson prior $\lambda$ away from the posterior knot count $k$: