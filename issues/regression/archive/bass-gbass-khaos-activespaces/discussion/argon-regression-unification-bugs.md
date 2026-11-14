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
