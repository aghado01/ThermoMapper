**The preparation is now substantially complete.** The earlier five findings are addressed in code or explicitly deferred to the SMC stage, and the dependency boundaries remain coherent. I found two new correctness gaps, one statistical decision worth revisiting, and a smaller LinAlg issue.

1. **Weighted quantiles are not invariant to merge order.**

   [Readouts.cs:144](/D:/aghado01/ThermoMapper/src/maths/samplers/ensemble/Readouts.cs:144) interpolates individual weighted observations without consolidating ties or excluding zero weights.

   Reproductions against the implementation:

   - Samples `(value, weight) = (0,1), (0,3), (1,1)` give a 70th percentile of **0.5 or 0**, depending on merge order.
   - `(10,0.5), (20,0.5), (−100,0)` gives a 5th percentile of **−78**, despite all positive mass being between 10 and 20.
   - A sample containing only `(123,0)` returns median **123** despite having zero total mass.

   This breaks the accumulator contract and matters for SMC, where duplicate particles and zero weights are ordinary cases. Exclude zero mass, handle ties consistently, and test permutation, merge-order, and weight-rescaling invariance.

2. **The new Normal centering path admits constant responses through rounding error.**

   [BandedGPrior.cs:116](/D:/aghado01/ThermoMapper/src/maths/regression/spline/BandedGPrior.cs:116) obtains the constant component through a numerical projection. Subtracting that reconstructed component leaves tiny residuals even when the original response is exactly constant.

   In a 40-observation cubic-spline fixture:

   | Constant response | Admission | Log marginal |
   |---|---|---:|
   | 0 | `ZeroResponse` | −∞ |
   | 1 | Accepted | 1378.208 |
   | 7 | Accepted | 1300.982 |

   The constant-one case reports conditional variance around **7 × 10⁻³³**. These numbers reflect projection roundoff. They contradict the [declared constant-response rejection](/D:/aghado01/ThermoMapper/issues/regression/bass-preparation-brief.md:1858) under the improper variance prior.

   Handle this degeneracy explicitly and preserve the constant component accurately when centering. Add nonzero constant-response fixtures alongside the existing noisy location-invariance tests.

3. **The plain-mean prior centre reintroduces outlier influence through the prior.**

   Keeping the centre independent of the sampled weights preserves the stated Gibbs conditional. However, using the ordinary response mean creates another route for an outlier to affect the fit.

   In an intercept-only fixture with 40 observations—39 zeros and one outlier assigned weight `1e-12`—the new centre produces:

   | Outlier | Conditional posterior mean |
   |---:|---:|
   | 1,000 | 0.625 |
   | 10,000 | 6.25 |
   | 1,000,000 | 625 |

   As the outlier’s likelihood weight tends to zero, the mean tends to \(R/n^2\), because the outlier still determines the prior centre. The conditional variance also grows substantially.

   **This is a statistical trade-off, not an incorrect Gibbs formula.** Shift invariance alone does not settle robust centering. Before carrying this choice into BASS’s D03, compare an explicit intercept treatment or a robust, weight-independent centre, and document the data-dependent prior interpretation. The shared block retaining `Zero` as its default keeps this decision appropriately local.

4. **The condition estimator gives invalid results for floored factors.**

   [BandCholesky.cs:326](/D:/aghado01/ThermoMapper/src/maths/linalg/BandCholesky.cs:326) combines normalization from the original matrix with solves from the floored factor.

   Factoring the one-dimensional matrix `[1e-20]` in floored mode returns an equilibrated condition estimate of **1e-8**. Every positive one-dimensional matrix has condition number **1**.

   Either refuse this estimate after flooring or compute it consistently for the modified matrix. **BARS’s strict factors avoid this path**, so this is a separate LinAlg contract correction.

The preparation has made substantial progress beyond fixing individual bugs. The full-support birth work identifies and removes a concrete mixing bottleneck; the subsequent tempering studies correctly remeasure performance after that kernel change. The Poisson work states its approximation and tests it against an independent integration reference. Readouts, ordered traces, and execution checkpoints now have distinct ownership, and Samplers acquires only downward mathematical dependencies.

I agree with keeping the remaining BARS tuning and SPC stopping-policy decisions parked. They need not delay BASS. The remaining BASS prerequisites are now specific: its statistical specification, dense incremental factor operations, basis/support rules, and reference fixtures. OQ10’s parameter-complete checkpoint identity belongs before dependable BASS resume; normalized reference distributions and evidence accounting belong before SMC.

At `1783921`, **82 sampler tests, 64 LinAlg tests, and 26 targeted regression tests passed** against the existing Release builds. I also verified the recorded 241-test regression result. The reproductions above exercise cases those tests omit.