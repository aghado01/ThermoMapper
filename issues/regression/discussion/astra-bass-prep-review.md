The dependency layout is coherent, but I found **two reproducible correctness gaps** behind the completed “exact robust weights” work. Three planned contracts also need tightening before BASS and SMC reuse.

I reviewed the [preparation brief](/D:/aghado01/ThermoMapper/issues/regression/bass-preparation-brief.md) against its supporting code, tests, and recorded results.

1. **Extreme outliers can receive full weight in the contaminated-normal model.**

   [WeightResampling.cs:85](/D:/aghado01/ThermoMapper/src/maths/regression/spline/bars/WeightResampling.cs:85) computes both mixture densities directly. When both underflow, their ratio becomes `0/0 → NaN`; the subsequent comparison is false, so the observation receives weight **1**.

   Calling the actual implementation with the defaults, σ² = 1, and residuals `[8, 100, 380, 400, 1000]` produced weights `[0.01, 0.01, 0.01, 1, 1]`. The largest residuals should overwhelmingly receive the outlier weight.

   Compute the responsibility using stable log odds and add an extreme-tail fixture. The existing conditional test covers moderate residuals and repeats the same ordinary-density arithmetic in its reference calculation.

2. **A weight update can move the current spline outside the target’s admitted support.**

   [BandedGPrior.cs:64](/D:/aghado01/ThermoMapper/src/maths/regression/spline/BandedGPrior.cs:64) applies admission to both the unweighted Gram matrix and the weighted posterior precision. But [BarsEnsemble.cs:348](/D:/aghado01/ThermoMapper/src/maths/regression/spline/bars/BarsEnsemble.cs:348) redraws weights without enforcing that support, then replaces the target while retaining the current knots.

   In a 40-observation cubic-spline probe, changing one valid weight from `1` to `0.01` changed the current configuration’s log marginal from **−4.71417 to −∞**, and coefficient readouts became **NaN**.

   This is more than missing error handling. If weighted-factor admission defines the target’s support, the weight conditional must respect that restriction. Otherwise, the supposedly exact sweep combines incompatible conditionals.

   Settle whether admission is a structural restriction based on the unweighted design or a restriction on the joint knot/weight state. Then make every transition preserve the declared support. This needs a full-sweep boundary test, beyond the existing fixed-design conditional checks.

3. **`MinRelativePivot` does not establish the proposed six-digit accuracy guarantee.**

   [Step 3.7:236](/D:/aghado01/ThermoMapper/issues/regression/bass-preparation-brief.md:236) retains that rationale for `MinRelativePivot ≥ 1e-10`. However, [FactorStatus.cs:14](/D:/aghado01/ThermoMapper/src/maths/linalg/FactorStatus.cs:14) correctly documents that its reciprocal is a **lower bound** on the equilibrated condition number—not an upper bound.

   A concrete counterexample is \(A=LL^\mathsf{T}\), where the 25×25 lower-bidiagonal \(L\) has diagonal 1 and subdiagonal 2. The implemented factor reports `PivotRatio = 1` and `MinRelativePivot = 0.2`, yet the equilibrated matrix’s condition number is at least **3.75 × 10¹⁴**.

   The statistic remains useful as a cheap rejection heuristic. It cannot certify accuracy. If accuracy qualification is the requirement, add a condition estimate and appropriate solve/residual checks. This distinction matters particularly for BASS’s broader dense designs.

4. **The stated `b = 0` endpoint is not necessarily the original prior.**

   The [linear-tempering report:1065](/D:/aghado01/ThermoMapper/issues/regression/bass-preparation-brief.md:1065) says the path starts at the prior while preserving `−∞` marginals at zero temperature. [SamplingTarget.cs:151](/D:/aghado01/ThermoMapper/src/maths/samplers/targets/SamplingTarget.cs:151) implements precisely that behavior.

   Consequently, states with finite prior density but zero likelihood remain excluded at `b = 0`. The reference is then the prior restricted to that support.

   This need not invalidate the current positive-temperature ladders. Before SMC, however, explicitly distinguish prior exclusions, zero likelihood, and numerical failure. Initialization and evidence accounting must use the reference distribution actually represented. A proper variance prior alone does not resolve this issue.

5. **Step 3.4 should separate accumulators, diagnostic traces, and execution checkpoints.**

   The [accumulator proposal:211](/D:/aghado01/ThermoMapper/issues/regression/bass-preparation-brief.md:211) combines commutative statistics, weighted particle readouts, and restart state.

   Mergeable readout statistics are appropriate shared infrastructure. Rank-normalized diagnostics and autocorrelation ESS additionally need ordered chain information; keeping only separate per-chain moments is insufficient. Restart state includes sampler configuration, RNGs, adaptation, replica assignments and swap-schedule state, which cannot be pooled like statistics.

   Use a checkpoint container that owns these distinct components. Keep MCMC diagnostic ESS distinct from particle-weight ESS. This preserves the dependency boundaries already established in the full-port brief.

Several parts of the preparation are particularly solid: the unweighted g-prior formulation, separating factor mechanics from model admission, neutral sampling-target contracts, and qualifying stochastic performance claims instead of asserting improvements from one seed. The recorded benchmark discussion is appropriately cautious.

I would reopen the robust-baseline qualification for findings 1–2, correct the numerical and reference-distribution claims, and clarify the accumulator boundary before proceeding through the remaining preparation stages. The brief’s opening **“not started”** status also needs reconciliation with its completion reports.

I verified the stored **207-test baseline**, plus the later **57/57 sampler** and **193/193 regression** results. My fresh suite attempt stopped before executing tests because the testhost dependency could not be resolved; the targeted compiled probes did execute. The review is anchored to `64d98a4`; concurrent diagnostics commits did not change the files supporting these findings.

No code, tests, or briefs were edited. I added only the durable review journal under `issues/regression/.doc-dive/`.