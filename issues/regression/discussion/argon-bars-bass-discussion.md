Look at the last two problems you uncovered:

1. **In the RJMCMC kernel:** A move-selection schedule ($b_k, d_k$) designed purely to cancel the prior out of Green's acceptance ratio was misattributed as a mixing optimization.
2. **In the Gram matrix solver:** A scale-dependent Cholesky heuristic (`PivotRatio` on $G = Z^T Z$) was silently doing double duty as a **statistical prior truncation** to stop Zellner's $g$-prior from exploding on data-starved B-splines.



When you pull a thread to unify BARS, BAPS, and BASS cleanly, those hidden couplings snap. That feels like getting bogged down, but you are actually uncovering the exact seam lines required to make the unified engine coherent.

Here is a sharp critique of the analysis in your note—including where the diagnosis is spot-on, where I call **BS** on one of its conclusions, and how to decide the two open questions.

---

### 1. What the Note Gets Right (Accept the Split)

Accept the two-part split (`MinRelativePivot` + explicit support restriction) and keep it sequenced *after* Step 2.3.

* **Never let a linear algebra solver silently define the statistical prior.** Using `PivotRatio` on an unscaled Gram matrix mixes floating-point stability with basis column norms. As the note correctly identifies, Cholesky accuracy is governed by the condition number of the diagonally equilibrated matrix $D^{-1/2} G D^{-1/2}$, which `MinRelativePivot` bounds directly.


* **The GLM / Robust weight trap is real.** In Poisson BARS (where the Fisher information weights by $\hat{\mu}_i$) and robust IRLS (where weights downweight outliers), column norms of $W^{1/2} Z$ shrink in low-count or high-outlier regions. Under `PivotRatio`, the sampler silently rejects valid knot placements in low-firing-rate regions purely because $\hat{\mu}_i$ scaled the column down. Switching the numerical gate to `MinRelativePivot` $\ge 10^{-10}$ with a strict, floor-free Cholesky fixes this permanently.


* **Sequencing after Step 2.3 is mandatory.** Step 2.3 is a structural refactor verified against the Step 0.6 baseline. If you change model admission rules at the same time as the refactor, you lose your regression invariant. Lock in Step 2.3 first, then land the admission split as a standalone, benchmarked diff.



---

### 2. Why 3–7 Points Still Blew Up (and Why BASS's `npart` Default Will Cripple BARS)

There are two critical mathematical subtleties with the `npart` recommendation:

#### Subtlety A: Why did configurations with 3–7 non-zero points *still* blow up to 226, even with `MinRelativePivot` $\in [0.1, 0.2]$?

configurations with 1–2 points reach 30–1,300, configurations with 3–7 points still hit 226 (at both 3 and 6 points), and blames the $g$-prior covariance $m\sigma^2(Z^T Z)^{-1}$ having "no constraint on a direction with no data".

Let's check the exact algebra of why that happens:

1. Under Zellner's $g$-prior ($\beta \mid \sigma^2 \sim \mathcal{N}(0, m\sigma^2 (Z^T Z)^{-1})$), the posterior mean is:

$$\mathbb{E}[\beta \mid y, Z] = \frac{m}{m+1} (Z^T Z)^{-1} Z^T y = \frac{m}{m+1} \hat{\beta}_{\text{OLS}}$$



Because the $g$-prior's precision matrix is proportional to $Z^T Z$, **its shrinkage factor $\frac{m}{m+1}$ is completely isotropic across all eigenvectors of $Z^T Z$**. When data in a direction is scarce, $(Z^T Z)$ shrinks, so the prior precision shrinks right along with the likelihood precision! It provides zero extra ridge regularization in data-starved directions.
2. Why did 6 non-zero points under a cubic B-spline still produce a fit of 226 while `MinRelativePivot` was $0.1$–$0.2$?
Because a cubic B-spline $B_j(x)$ supported on $[\xi_j, \xi_{j+4}]$ decays **cubically** at its boundaries:



$$B_j(x) \propto (x - \xi_j)^3 \quad \text{as } x \downarrow \xi_j$$



Suppose a knot interval spans a data gap, and 6 data points sit just inside the edge of $B_j$'s support where $B_j(x_i) \approx 0.003$.
* What does `sum(B_j != 0)` see? It counts **6 non-zero points** (`npart = 6`).


* What does `MinRelativePivot` see? Diagonal equilibration divides column $z_j$ by its Euclidean norm $\Vert{}z_j\Vert{}_2 = \sqrt{\sum_i B_j(x_i)^2} \approx 0.007$. Once normalized to unit length, $z_j / \Vert{}z_j\Vert{}_2$ is nearly orthogonal to its neighbors, so its relative pivot is healthy ($0.1$–$0.2$)!


* **What does $\hat{\beta}_{\text{OLS}}$ do?** To fit a tiny residual $r_i \approx 1.5$ at those edge points where $B_j(x_i) \approx 0.003$, least squares sets $\hat{\beta}_j \approx r_i / B_j(x_i) \approx 500$. Then, in the middle of the data gap where there are no data points, the B-spline $B_j(x)$ rises to its peak value of $\approx 2/3$, producing a phantom spike of $\hat{f}(x_{\text{gap}}) \approx 500 \times \frac{2}{3} \approx 333$!



That is why `PivotRatio` caught some of them: `PivotRatio` was dividing the unnormalized pivot $\Vert{}z_j^\perp\Vert{}_2^2$ by the largest pivot $\max_k \Vert{}z_k\Vert{}_2^2$, which was implicitly thresholding the **column energy** $\Vert{}z_j\Vert{}_2^2 = \sum_i B_j(x_i)^2$, not the unweighted count of non-zero entries $\sum_i \mathbf{1}(B_j(x_i) > 0)$

#### Subtlety B: Do NOT use BASS's `min(20, 0.1 * n)` as the BARS default!

At the end of the note, the recommendation is to adopt BASS's default `npart = min(20, 0.1 * n)` (which is 15 points per basis function at $n = 150$) for BARS.

Look at the basis geometry difference between BASS and BARS:

* **BASS** uses tensor products of MARS-style truncated linear hinges $[s(x_d - t)]_+^1$ supported on a **global half-space** $[t, 1]$ or $[0, t]$. Requiring $10\%$ of the data ($\min(20, 0.1n)$) to be non-zero in a half-space only keeps knots out of the extreme outer $10\%$ tails of the domain.


* **BARS** uses **compactly supported local cubic B-splines** spanning 4 knot intervals $[\xi_j, \xi_{j+4}]$ in the interior—and at a clamped boundary knot of multiplicity 4, the boundary B-spline $B_0(x)$ is supported on **a single knot interval $[\xi_0, \xi_1]$**!


* If you set `npart = 15` on $n = 150$ in BARS:


1. You forbid the first interior knot $\xi_1$ from ever being placed within the first $10\%$ of the data ($15$ points), because the boundary B-spline $B_0(x)$ is zero outside $[\xi_0, \xi_1)$.
2. You forbid placing closely spaced knots across a sharp local transient (like a neural PSTH onset or a critical temperature $T_c$ susceptibility spike), which is the entire reason to use BARS over fixed-knot splines.

---

### 3. Recommended Unified Support Gate

Instead of forcing BARS to use an unweighted non-zero indicator count with BASS's half-space threshold (`min(20, 0.1*n)`), define the support gate using **Effective Basis Weight (or Effective Sample Size per basis column)**:

| Basis Property | BASS (Global Half-Space Hinges) | BARS / BAPS (Local B-Splines) |
| --- | --- | --- |
| **Support geometry** | Half-space $[t, 1]^d$ | Local 4-interval window $[\xi_j, \xi_{j+4}]$ (1 interval at boundary)

 |
| **Why raw `count(B_j != 0)` fails on B-splines** | Works fine on half-spaces

 | Counts points in the $(x-\xi)^3 \sim 10^{-3}$ cubic tail as full points, allowing 226+ spikes at 6 points

 |
| **Scale-free Effective Support Metric** | $N_{\text{eff}}(j) = \sum_{i=1}^n \mathbf{1}(B_j(x_i) > 0)$<br> | **Column $L^1$ or $L^2$ mass:** $w_j = \sum_{i=1}^n B_j(x_i)$ or $N_{\text{eff}}(j) = \frac{(\sum_i B_j(x_i))^2}{\sum_i B_j(x_i)^2}$ (or counting points where $B_j(x_i) \ge \tau_{\text{tail}}$, e.g. $0.05$) |
| **Sensible Default Threshold** | `min(20, 0.1 * n)`<br> | `min_effective_weight` $\approx 2.0\text{–}3.0$ (or `npart = 4–6` points with $B_j(x_i) \ge 0.05$, or $\ge 2$ points per knot interval) |

Notice what requiring $\sum_{i=1}^n B_j(x_i) \ge w_{\min}$ (with $w_{\min} \approx 2.5$ or $3.0$) does for BARS:

* Because cubic B-splines form a partition of unity ($\sum_j B_j(x) = 1$), $\sum_{i=1}^n B_j(x_i)$ literally equals the **effective number of data points belonging to basis function $j$**.
* Those 6 tail points where $B_j(x_i) = 0.01$ have total weight $\sum_i B_j(x_i) = 0.06 \ll 2.5$, so they are **immediately rejected**—eliminating the 3-to-7-point blow-ups (`226`) without needing an artificially huge `npart = 15` that ruins sharp peaks!


* Meanwhile, 4 points sitting squarely under a narrow sharp peak where $B_j(x_i) \approx 0.65$ have $\sum_i B_j(x_i) \approx 2.6 \ge 2.5$, so sharp local peaks are **admitted**.

---

### 4.The 4-Gate Architecture

To stop feeling like every unification step opens a new rabbit hole, lock down the boundary between your modules so each concern lives in exactly one place:

1. **Gate 1 — Prior Support / Geometry Gate (`IBasisSupportRule`, evaluated *before* Cholesky):**
* Checks whether candidate knots give every basis function sufficient data support.


* Plug-in implementation: `NonZeroCountSupport(npart)` for BASS (`min(20, 0.1*n)`), and `EffectiveMassSupport(min_mass)` (or `ThresholdedCountSupport(npart, min_b_val)`) for BARS/BAPS.

2. **Gate 2 — Numerical Conditioning Gate (`IGramSolver`, evaluated *during* Cholesky):**
* Strictly scale-invariant: no diagonal floor (`FlooredPivots` removed), requires positive-definiteness and `MinRelativePivot >= 1e-10` on the equilibrated matrix.

3. **Gate 3 — Coefficient Prior (`ICoefficientPrior`):**
* Keeps Zellner's $g$-prior ($\tau = 1/m$ on $Z^T Z$) for exact BARS/BASS parity right now, while leaving a clean slot for a ridge-stabilized or mixtures-of-$g$ prior later so you don't have to redesign the engine when you want proper gap shrinkage.

4. **Gate 4 — RJMCMC Move Schedule (`IMoveSchedule`):**
* Pluggable `MaskedConstantSchedule` vs. `GreenDMGKSchedule`, evaluated on both $\text{ESJD}(k)$ and $\text{ESS}$.

### Immediate recommendations 

1. **Approve the split and the timing:** Switch the numerical check to strict `MinRelativePivot >= 1e-10` (no floor) + an explicit support check, implemented immediately *after* Step 2.3 lands against the 0.6 baseline.

2. Reject BASS's `min(20, 0.1*n)` as the BARS default: the cubic B-spline tail effect ($B_j(x) \propto (x-\xi)^3$) and boundary clamping. 

Have the probe log $\sum_i B_j(x_i)$ (column $L^1$ mass) and $\Vert{}z_j\Vert{}_2^2$ (column $L^2$ energy) alongside integer `npart` on those exact 3–7 point gap configurations—you will almost certainly see that a tiny threshold on effective mass ($\sum_i B_j(x_i) \ge 2.5$) cleanly separates every wild fit from sane fits without sacrificing sharp-peak resolution.



Want to inspect what the column $L^1$/$L^2$ mass cutoff looks like on clamped boundary splines, or sketch the interface connecting BARS, BAPS, and BASS ? 

Seeing `bass-preparation-brief.md` clarifies the whole picture. The brief's architecture is remarkably disciplined. 

* **Phase 0 was defect repair:** Every step fixed an unambiguous mathematical bug with a failing reproduction.

* **Phases 2 and 3 are seam extraction under a baseline lock (`0.6` $\to$ `3.5`):** The goal of Steps `2.3` through `3.4` is to extract shared infrastructure (`CollapsedGPriorBlock`, `IRjTarget`, `ParallelTempering`, `Accumulator`) so BASS can plug right in at Stage 4, while BARS's statistical behavior stays locked to the `0.6` baseline.



The friction is happening because inspecting the seams between BARS and BASS keeps exposing places where BARS's legacy conventions (`PivotRatio`, the `DMGK` schedule, `EbicPrior`) differ from BASS's conventions (`MinRelativePivot` + `npart`), tempting a redesign of BARS's statistical model in the middle of an infrastructure refactor.

Here is how to keep the boundaries sharp across your remaining steps (`2.3` and `3.1`–`3.6`) and close out the open questions without losing momentum.

---

### 1. Where the `PivotRatio` / `MinRelativePivot` / `npart` Decision Actually Belongs

In `bayesian-regression-questions.md`, the probe tested swapping `PivotRatio` for `MinRelativePivot` in BARS and noticed that `PivotRatio` was accidentally acting as a rough basis-support filter in data gaps.

Look at how that maps onto your architecture in `bass-preparation-brief.md`:

1. **`npart` does NOT belong in Step 2.3's collapsed $g$-prior block.**
Even in BASS (`rjmcmc_des_*.R`), `npart` (`sum(cand.des$basis != 0) < npart.des`) is checked in the **basis/proposal layer before any Gram matrix is factored**. The $g$-prior math block should know nothing about B-spline intervals or MARS half-spaces; its only job is to take a factored Gram system (whose `FactorStatus` passed the numerical gate) and compute marginal likelihoods, explicit residuals, tempered densities, and $(\sigma^2, \beta)$ draws.


2. Do not force BASS's `npart = min(20, 0.1*n)` onto BARS.
Minimizing redundancy means sharing *mechanisms*, not forcing a 1D local cubic B-spline model (where a clamped boundary spline spans 1 knot interval and interior tails decay as $(x-\xi)^3$) to use the same integer point-count heuristic as a multi-dimensional global half-space hinge model.


* In BASS, `npart` counts non-zero indicator hits in a tensor-product half-space.


* In BARS, as your probe showed, 3–7 points in the cubic $(x-\xi)^3$ tail of a B-spline near a gap edge still produce wild fits of $226$ because their basis values $B_j(x_i)$ are tiny while `MinRelativePivot` equilibrates the column to unit norm. If BARS adds an explicit support rule in `SplineTarget`, thresholding the **column mass** $\sum_i B_j(x_i) \ge w_{\min}$ (e.g., $w_{\min} \approx 2.0\text{–}3.0$) or points above a tail cutoff ($B_j(x_i) \ge 0.05$) matches B-spline geometry without blinding BARS to sharp narrow peaks.




3. Keep Step 2.3 and Step 3.5 completely clean first.
Leave BARS's `DesignAdmission` untouched during Steps `2.3`–`3.5` so your Step `2.3` refactor and Step `3.5` baseline rerun match Step `0.6` to machine rounding. Once `3.5` passes, land the BARS `MinRelativePivot` + spline-support rule as its own isolated commit with before/after rejection-reason counts from the Step `3.4` accumulator.



---

### 2. Designing Step 2.3 (Collapsed $g$-Prior Block) for Minimal Redundancy

Step `2.3` is the most important shared math seam between BARS and BASS. Notice what is genuinely shared versus what is engine-specific:

* **What differs:**
* **BARS** factors a bandwidth-4 matrix via `BandCholesky`, fixes $\tau = 1/m$, integrates out $\sigma^2$, and in robust mode (`0.5d`) factors *two* banded matrices ($G = Z^T Z$ and $P_a = Z^T(W + I/m)Z$) because $W \ne I$.


* **BASS** factors a dense matrix $G = B^T B$ via `DenseCholesky` (updated by bordering/Givens in Stage 2) and samples $\tau$.




* **What is 100% shared once the linear solve is done:**
In the unweighted normal case ($W = I$, which covers standard BARS and BASS) with prior $\beta \mid \sigma^2 \sim \mathcal{N}(0, \sigma^2 \tau^{-1} G^{-1})$ and likelihood inverse-temperature $\beta_{\text{temp}} \in (0, 1]$ (with $\beta_{\text{temp}} = 1$ for untempered):
1. **Only one factorization ($G$) is ever needed when $W = I$.** Because $P_a = (\beta_{\text{temp}} + \tau) G$, the determinant ratio between the prior normalizing constant $\vert{}\tau G\vert{}^{1/2}$ and the posterior Laplace integral $\vert{}(\beta_{\text{temp}} + \tau)G\vert{}^{-1/2}$ simplifies analytically:

$$\frac{\vert{}\tau G\vert{}^{1/2}}{\vert{}(\beta_{\text{temp}} + \tau) G\vert{}^{1/2}} = \left(\frac{\tau}{\beta_{\text{temp}} + \tau}\right)^{p/2}$$



Notice that $\log\vert{}G\vert{}$ cancels out completely when $W = I$! You only need $\log\vert{}G\vert{}$ and $\log\vert{}P_a\vert{}$ separately on the robust weighted path ($W \ne I$) from Step `0.5d`.


2. **The cancellation-free explicit residual $S$ works identically for both:**
With $\hat{\beta}_{\text{OLS}} = G^{-1} Z^T y$, $r_{\text{OLS}} = y - Z\hat{\beta}_{\text{OLS}}$, and $Q = y^T Z \hat{\beta}_{\text{OLS}} = \Vert{}L^{-1} Z^T y\Vert{}_2^2$:

$$S(\tau, \beta_{\text{temp}}) = \beta_{\text{temp}} \left( \Vert{}r_{\text{OLS}}\Vert{}_2^2 + \frac{\tau}{\beta_{\text{temp}} + \tau} Q \right)$$



At $\beta_{\text{temp}} = 1$ and $\tau = 1/m$, this is algebraically identical to $y^T y - \frac{m}{m+1} Q$, matches Step `0.6` to rounding for BARS, and eliminates catastrophic cancellation when BASS samples $\tau \ll 1$ on a near-interpolating fit.


3. O(1) Likelihood-Tempered Swaps (Step 3.2, Q5) fall right out:
Because $\Vert{}r_{\text{OLS}}\Vert{}_2^2$, $Q$, and $p$ do not depend on $\beta_{\text{temp}}$, caching those three numbers inside the block makes evaluating the log-marginal likelihood at *any* replica temperature $\beta_{\text{temp}}$ in Step `3.2` an $O(1)$ scalar calculation with zero matrix work.


---

### 3. Resolving the Open `EbicPrior` Question from Step 0.5e

In the Step `0.5e` report, you flagged an open question:

> *"`EbicPrior`'s count law. Its law is now exactly $C(n_c, k)^{-\gamma}$, which is symmetric in $k \leftrightarrow n_c - k$, so it penalizes only below $n_c/2$. He, Yang and Kang's prior gives each of the $C(n, k)$ models the mass $C(n, k)^{-\gamma}$, so their count law is $C(n, k)^{1-\gamma}$: flat at $\gamma = 1$ and peaked at $n/2$ below it. Which law is intended is undecided."*
> 

Here is the exact mathematical resolution:

* **Why He, Yang, and Kang have $p(k) \propto \binom{n_c}{k}^{1-\gamma}$:** Their formulation lives on a **discrete grid of $n_c$ candidate knots**. There are $\binom{n_c}{k}$ distinct models of dimension $k$. When you assign each individual model prior weight $\binom{n_c}{k}^{-\gamma}$ and sum over all $\binom{n_c}{k}$ models in dimension $k$, the discrete multiplicity $\binom{n_c}{k}$ multiplies the per-model weight, giving marginal count mass $\binom{n_c}{k} \cdot \binom{n_c}{k}^{-\gamma} = \binom{n_c}{k}^{1-\gamma}$.


* **Why continuous-knot BARS is different after Step 0.5e:** In continuous BARS, the within-dimension knot prior on the open ordered simplex $\Delta^k = \{0 < \xi_1 < \dots < \xi_k < 1\}$ has density $p(\xi \mid k) = k!$, which **already integrates to $1$** ($\int_{\Delta^k} k! \, d\xi = 1$). There is no $\binom{n_c}{k}$ multiplicity of discrete models inside dimension $k$—the continuous volume factor $1/k!$ was already canceled by the `+ log k!` term you added in Step `0.5e`!


* **Therefore:**
1. For `IKnotCountPrior` to honor the contract established in Step `0.5e`—namely, that `p(k)` **is** the marginal prior law on $k$ verified by `CountLaw.Sample`—`EbicPrior` in a continuous-knot engine should represent the marginal count penalty $p(k) \propto \binom{n_c}{k}^{-\gamma}$, **not** $\binom{n_c}{k}^{1-\gamma}$ (which at $\gamma < 1$ would actively *reward* adding knots up to $n_c / 2$!).


2. To prevent the $\binom{n_c}{k} = \binom{n_c}{n_c - k}$ symmetry from rewarding $k > n_c / 2$, `EbicPrior` should explicitly restrict support to $k \le \lfloor n_c / 2 \rfloor$ (returning $-\infty$ above $\lfloor n_c / 2 \rfloor$), since EBIC's $\log \binom{n_c}{k}$ penalty is only a sparsity penalty in the sub-half-space $k \le n_c / 2$.





---

### 4. Recommended Execution Path from Here

You don't have a structural problem—your DAG in `bass-preparation-brief.md` is tight and nearly at the finish line:

1. **Park the BARS admission change (`MinRelativePivot` + B-spline support threshold) until after Step `3.5**` so Steps `2.3`–`3.5` stay a pure, zero-drift refactor against the `0.6` baseline.


2. **Execute Step `2.3` (`CollapsedGPriorBlock`)** with the unweighted ($W=I$, single factor, cached $\Vert{}r_{\text{OLS}}\Vert{}_2^2$ and $Q$ for $O(1)$ tempering) and weighted ($W \ne I$, dual factor $G$ and $P_a$) paths cleanly separated.


3. **Run straight through `3.1` $\to$ `3.6**` (`IRjTarget`, likelihood-tempered swaps, non-reversible even/odd schedule, and the merged `Accumulator`), where every piece has a clear mechanical spec.



Want to sketch the exact C# struct/interface surface for Step `2.3`'s `CollapsedGPriorBlock` so it serves both `WeightedNormalModel` (`BandCholesky`) and Stage-4 BASS (`DenseCholesky`) without extra allocations?
[User did not ask, but user is always interested in avoiding extra allocations]
