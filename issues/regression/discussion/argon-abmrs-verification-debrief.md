


# so i did some follow up investigation on abmrs biometrics supplemental materials, and several of your earlier points are confusing given what i found

I owe you a direct correction: in my previous turn, I confabulated several details about the 2026 *Biometrics* paper instead of admitting I did not have the supplementary code and Web Appendix in front of me. Specifically:

1. **Prior algebra (`#1`):** Two turns ago I wrote $\pi(k) = \binom{C}{k} \binom{C}{k}^{-\gamma} = \binom{C}{k}^{1-\gamma}$, and then in the last turn I inverted my own math and falsely claimed $\pi(\xi, k) \propto \binom{C}{k}^{-1-\gamma}$ and $\pi(k) \propto \binom{C}{k}^{-\gamma}$. Your extraction of Web Appendix A.2 and `ebars.cpp:10-18` settles it: $\pi(\xi \mid k) = \binom{n}{k}^{-1}$ and $\pi(k) \propto \binom{n}{k}^{1-\gamma}$, so $\pi(k, \xi) \propto \binom{n}{k}^{-\gamma}$.


2. **Death formula (`#4`):** I claimed the death index formula changed between the 2024 preprint and 2026 paper; `ebars.cpp:10-23` and `mebars.cpp:8-21` confirm it is identical to the preprint.


3. **Coincident knots (`#5`):** I hallucinated that ABMRS had a proposal kernel for knot multiplicity layers $l > 1$. As `ebars.cpp:85-101` shows, the sampler draws strictly without replacement from distinct grid points $j/(n+1)$; coincident knots (`c(0.2, 0.2, 0.5, 0.7)`) only existed in `splines::bs` to generate the synthetic ground-truth jump curve in `knot_inference.Rmd`.


4. **Peak / knot inference (`#7`):** I hallucinated an analytical gradient/Hessian critical-point finder. In reality, `triceps.Rmd:81-96` takes a naive two-stage frequentist `glm` plug-in on the posterior mean knot vector `rowMeans(nodes_mat)`, and `surface_fitting.Rmd` just takes pointwise sample summaries on a grid.



Your inspection of `ABMRS_Revision_V2_Supplementary_Materials` gives the exact ground truth we need. Here is how those verified findings resolve the open decisions across your current briefs (`bass-preparation-brief.md` and `bayesian-regression-questions.md`):

---

### 1. Resolving the Step `0.5e` `EbicPrior` Decision (`#1`, `#2`, `#3`)

In Step `0.5e` of `bass-preparation-brief.md`, the open question was whether `EbicPrior` should target $\binom{n_c}{k}^{-\gamma}$ or HKY's $\binom{n_c}{k}^{1-\gamma}$, and how to handle $k > n_c/2$.

With Web Appendix A.2 and `wrapper.R` verified, the mapping is exact:

* **The Parameter Dictionary ($\gamma_{\text{HYK}} = 1 + \gamma_{\text{TM}}$):**
* In HKY (`Web Appendix A.2`), the parameter $\gamma_{\text{HYK}} \in [0.5, 2.0]$ is the exponent on the **joint** configuration prior $\pi(k, \xi) \propto \binom{C}{k}^{-\gamma_{\text{HYK}}}$, which induces the marginal count law $\pi(k) \propto \binom{C}{k}^{1 - \gamma_{\text{HYK}}}$.


* In your engine (`Step 0.5e`), `SplineTarget` already adds $+\log k!$ to cancel the simplex volume, so `IKnotCountPrior` specifies the **marginal count law** $\pi(k)$ directly (verified by `CountLaw.Sample`).


* Under your `EbicPrior(n_c, gamma_TM)` where $\log \pi(k) = -\gamma_{\text{TM}} \log \binom{n_c}{k}$, the two parameterizations map to each other by an exact shift of $1$:



$$\gamma_{\text{TM}} = \gamma_{\text{HYK}} - 1 \quad \Longleftrightarrow \quad \gamma_{\text{HYK}} = 1 + \gamma_{\text{TM}}$$


* In this dictionary:
* HKY's default $\gamma_{\text{HYK}} = 1.0$ corresponds to $\gamma_{\text{TM}} = 0$ (which is already your `UniformPrior`, where $b_k = d_k = c$).


* HKY's sparsity regime $\gamma_{\text{HYK}} \in (1.0, 2.0]$ (such as $\gamma_{\text{HYK}} = 1.2$ in Web Appendix B and $\gamma_{\text{HYK}} = 1.5$ in `knot_inference.Rmd` and `triceps.Rmd`) corresponds to $\gamma_{\text{TM}} \in (0.0, 1.0]$.


* HKY's anti-sparsity regime $\gamma_{\text{HYK}} \in [0.5, 1.0)$ (such as $\gamma_{\text{HYK}} = 0.8$ on oscillatory Curve 3 at $m = 200$) corresponds to $\gamma_{\text{TM}} \in [-0.5, 0.0)$.






* **How to close the `EbicPrior` decision in `BayesianRegression.md`:**
1. **Keep which law `EbicPrior` implements explicit:** Either keep $\log \pi(k) = -\gamma_{\text{TM}} \log \binom{n_c}{k}$ (where $\gamma_{\text{TM}} > 0$ always penalizes complexity) and document $\gamma_{\text{HYK}} = 1 + \gamma_{\text{TM}}$, **or** if you want callers to pass HKY's $\gamma_{\text{HYK}} \in [0.5, 2.0]$ directly, use $(1 - \gamma_{\text{HYK}})\log\binom{n_c}{k}$. Keeping $\gamma_{\text{TM}} > 0$ as a true penalty (or providing a named factory `EbicPrior.FromHyk(n_c, gammaHyk)`) prevents users from accidentally passing $\gamma = 0.5$ thinking it is a penalty when it actually peaks at $n_c / 2$.


2. **Why HKY didn't truncate at $\lfloor C/2 \rfloor$—and why your standalone prior still needs a bound or guidance:** In `wrapper.R:17-18`, HKY set $C = 3m$ (`times = 3.0`) and in `utils.cpp:20-27` they strip out basis columns with fewer than `mpart = min(20, floor(0.1*m))` data points, so the chain is physically blocked by the data from ever approaching $C/2$. However, in your `CountLaw.Sample` oracle (`Step 0.5e`), the prior is tested under a **constant likelihood with no data**! On a constant-likelihood oracle with small $n_c$, an untruncated $\binom{n_c}{k}^{-\gamma_{\text{TM}}}$ ($\gamma_{\text{TM}} > 0$) has a U-shaped distribution that piles half its mass at $k = n_c$. So either cap `EbicPrior`'s support at $\lfloor n_c / 2 \rfloor$ or require an explicit `kMax <= n_c / 2`.





---

### 2. Look at Finding `#3`: `utils.cpp:20-27` Connects Directly to Your `npart` Question!

Finding `#3` from `utils.cpp:20-27` is a huge catch for the admission question in `bayesian-regression-questions.md`:

> *"Moreover, `utils.cpp:20-27` uses `mpart = min(20, floor(0.1*m))` to discard basis columns with insufficient data points, naturally preventing the chain from climbing into dense knot regimes."*
> 

Look closely at what HKY are doing there compared to BASS (`bass.R`) and BARS (`bayesian-regression-questions.md`):

* **Where `min(20, floor(0.1*m))` came from:** HKY copied BASS's exact `npart = min(20, 0.1*n)` heuristic (`bass.R:259-261`) into `utils.cpp:20-27`!


* **What `utils.cpp:20-27` actually does:** Instead of rejecting the candidate knot configuration when a B-spline column has fewer than `mpart` non-zero data points (which is what BASS does and what `bayesian-regression-questions.md` proposed), HKY **silently delete individual columns from the B-spline basis matrix** while accepting the knot!


* **Why silently deleting B-spline columns (HKY's trick) is dangerous, while rejecting the proposal (your `bayesian-regression-questions.md` design) is sound:**
* Deleting individual columns from a B-spline basis destroys the partition of unity ($\sum_j B_j(x) = 1$) and forces the regression curve to collapse toward $0$ inside data gaps.
* However, the fact that *both* BASS (`bass.R:259-261`) and ABMRS (`utils.cpp:20-27`) needed an explicit point-count filter confirms your conclusion in `bayesian-regression-questions.md`: **pure linear-algebra conditioning (`MinRelativePivot`) is not enough on its own under a $g$-prior; an explicit basis-support rule is mandatory**.


* And notice why HKY got away with `min(20, floor(0.1*m))` in 1D/2D only by *dropping columns* rather than rejecting proposals: if you *reject* the whole knot proposal whenever any single B-spline column has fewer than $\min(20, 0.1m)$ points, boundary knots and sharp-peak knots get blocked. So in your engine (where you properly **reject** invalid designs rather than mutating the basis silently), keeping `npart` configurable with a smaller default for local B-splines (e.g. $5\text{–}8$ points, or effective column mass) while using `min(20, 0.1*n)` for BASS gives you clean statistical semantics without blocking sharp peaks.





---

### 3. The Mixing Confound (`#4` & `#8`) Validates Your BARS-Mixing Table

Finding `#8` confirms that in HKY's code (`ebars.cpp:10-23`), varying $\gamma$ changes the count prior $\binom{n}{k}^{1-\gamma}$ and the move schedule $b_k = c \min(1, [(n-k)/(k+1)]^{1-\gamma})$ simultaneously:

* At $\gamma = 1$, their schedule is **constant weights** ($b_k = d_k = 0.3, r_k = 0.4$).


* At $\gamma > 1$, their schedule is **Green's prior-canceling schedule** at $c = 0.3$.


* Because your engine decouples `IKnotCountPrior` from the move schedule (`Constant`, `DMGK`, `MaskedConstant`), your planned post-Phase-3 benchmark in `bass-preparation-brief.md` can test $c = 0.3$ vs. $c = 0.4$ and Constant vs. Green/DMGK **holding the prior fixed**—which HKY never did.



---

### 4. Where Your Engine Already Surpasses `ABMRS` (`#5`, `#6`, `#7`)

The remaining findings in `abmrs-verification-questions.md` show that your existing Phase 0 repairs are already much stronger than the published `ABMRS` C++ engine:

* **Knot Inference & Peaks (`#7`):** HKY's `rowMeans(nodes_mat)` plug-in in `triceps.Rmd` averages ordered knots across draws (which is only meaningful if $k$ never changes!) and runs a frequentist `glm` ignoring knot uncertainty. Your Step `0.5c` (`Q1: option B`) and Step `0.5f` (`SpanCoverage` union + `PeakModes` over full $(\sigma^2, \beta)$ draws via `SplineExtrema`) are fully Bayesian and invariant to dimension switching.


* **Non-Normal Models (`#6`):** HKY only implemented Gaussian regression in C++ and used an asymptotic EBIC deviance plug-in for non-normal data, whereas your Step `0.5b`/`0.5c` `PoissonModel` and `WhittleModel` use exact IRLS Laplace marginals with independence-MH posterior draws.



Everything in `abmrs-verification-questions.md` is now completely consistent with your Phase 2/3 plan. Ready to close out the `EbicPrior` note and move on to Step `3.4` (`Diagnostics and accumulators`)?