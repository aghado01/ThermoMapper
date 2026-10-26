# BASS port: response bases, prediction and sensitivity

A segment of the [full port brief](../bass-full-port-brief.md), which holds the scope and completion contract, the dependency boundaries, the stage table and the decision register. This segment covers stages 5, 6 and 9: supplied-basis and PCA fits, prediction and uncertainty, sensitivity, and the presentation data above maths, with the stage-1 decisions D07 and D08. Its checklist, where the reports go as work lands, is [port-outputs-checklist.md](port-outputs-checklist.md).

The decision drafts here are proposals for the user's design pass, not decisions. The register, with each decision's status, and the table of sources the drafts were drawn from are in the hub's [decisions section](../bass-full-port-brief.md#statistical-decisions-to-settle).

## Response basis conventions

R's PCA construction uses `Y_scaledᵀ = U D Vᵀ`, basis `K = U_r D_r`, and coefficient responses `W = V_rᵀ`. Passing ordinary PCA scores directly into the component regressions would change response scale, priors, and residual-noise interpretation.

Handle rank, constant output columns, full variance retention, component signs, and repeated singular values explicitly. Compare reconstructions or subspaces when individual components are not identifiable. Supplied bases can be nonorthogonal. Independent coefficient fits still require cross-component terms when sensitivity is reconstructed under the input measure.

## Decision drafts: D07, D08

### D07 (draft): response units and uncertainty

Sources: the inventory's PCA, prediction and defect sections, and R's [`bassPCA` setup](D:/aghado01/packages/BASS/R/bassPCA.R:88) and [basis prediction](D:/aghado01/packages/BASS/R/bassPCA.R:218).
- **Scalar and direct functional responses.** The noise is Gaussian, independent across cases and output locations, with one variance. A predictive draw is f(x, t) under the draw plus noise drawn independently for each case and location. The pointwise unit defaults to a case (a whole curve), with the product of its locations' densities given the draw (D04's retained-parameters content).
- **Decomposition** (`bassPCA`):
  - Centre by y_m (default: the mean) and scale by y_s (default: 1).
  - Take a thin SVD of the standardized responses, Y_cᵀ = U D Vᵀ. The basis is K = U_r D_r, and the coefficient responses are W = V_rᵀ, with rows of unit norm. This is R's convention, kept because the priors and the meaning of each component's noise depend on it.
  - The rank is the smallest r whose cumulative variance share is at least the target, capped at the numerical rank. R's strict comparison fails at 100%; that becomes a defect fixture. Full retention means the numerical rank and zero truncation.
  - A constant output column is centred and gets y_s = 1, never a division by zero.
  - Each component's sign is fixed by a convention for comparisons only. Where singular values repeat, comparisons use reconstructions or subspaces.
- **Truncation residuals.**
  - A training case's residual is e_i = y_{c,i} − K wᵢ, observed minus reconstructed, in standardized units.
  - A predictive draw adds y_s ∘ e, in physical units. e is a whole training-case vector, which keeps its within-output correlation, drawn uniformly over the training cases independently for each prediction case and draw.
  - R stores reconstructed minus observed, adds it after restoring physical units, and indexes columns with the row count: three defect fixtures.
- **Component noise.** Each component's ε_k ~ N(0, σ_k²) is drawn for each prediction case and posterior draw. R draws one per posterior draw and reuses it across cases (defect fixture).
- **The predictive draw** for a case is y = y_m + y_s ∘ (K(w(x) + ε) + e). The latent draw omits ε and e. Latent draws include coefficient variation through D04's composition draws, so a mixture of conditional means is never reported as the full posterior.
- **Joint draws across components.** The components are fitted independently on separate coefficient responses, so their posteriors are independent and equal-weight MCMC draws pair validly by index. Weighted draws, from SMC, are resampled within each component before pairing; pairing them by index is invalid.
- **Covariance in physical units** is diag(y_s)·Cov·diag(y_s). Calibration (D09) and basis sensitivity (D08) apply the same transform. R's calibration uses the unscaled basis and truncation covariance (defect fixture).
- **Pointwise density for supplied-basis fits (OQ17).** The default unit is a case and the content is Σ_k log N(w_ik | w_k(xᵢ), σ_k²) in coefficient space. For a fixed basis this differs from the density in y-space by a constant and a residual term that are the same for every fit, so it is exact for comparing fits that share a basis. Comparing across bases or engines needs a model of the truncation residual. A Gaussian with diagonal empirical truncation variance (calibration's functional form) is offered as a declared approximation, not as the default.
- **Inputs.** One fitted transform is reused everywhere.
  - A constant training column maps to 0. R's training helper allows it, while its prediction helper divides by the zero span (defect fixture).
  - An unseen category is refused by default.
  - Continuous extrapolation is allowed, with each case flagged when it lies outside the training range.
- **Axes and order.** Draw × case × output axes are always explicit, including singletons, and the requested draw order is kept. R reorders by model group and drops axes (defect fixture).
- **Intervals** are quantiles of the predictive mixture from draws, not quantiles of each draw's Gaussian interval endpoints, which R's diagnostic plot uses.
- **Not in the baseline:** basis-estimation uncertainty (Fast Bayesian FPCA, arXiv 2412.11340, is one route if it is ever modelled) and robust bases (OQ16). Both are in the [extensions brief](extensions-brief.md).
- **Independent checks:**
  - the training reconstruction y = y_m + y_s ∘ (K W + e) holds exactly;
  - the edge cases: 100% retention, full rank, constant columns;
  - the residual's sign, on an asymmetric fixture;
  - independence of noise across prediction cases;
  - draw order and singleton axes;
  - the covariance transform against a dense computation.

### D08 (draft): sensitivity measures and moments

Sources: the JSS paper (§3), Francom et al. 2018, the inventory's Sobol section and its defect table, and R's [`sobol`](D:/aghado01/packages/BASS/R/sobol_normTmix.R:194).
- **What is decomposed.** The latent surface's variance under a declared independent product measure over the inputs, computed for each posterior draw, which gives a posterior distribution of indices. Observation noise is not in the denominator.
- **Which surface:**
  - *Gaussian BASS:* the indices decompose Ba, which is the mean surface.
  - *Skewed noise* (GBMARS with β ≠ 0): E[y | x] = Ba + √w·β·E[v]. The offset is constant in x, so the mean's indices equal Ba's.
  - *Fixed-β asymmetric Laplace:* Ba is the q-quantile surface, and indices per quantile (GBMARS's SIR example) decompose a different functional for each q.

  The contract takes the decomposed surface as a declared input, the mean by default, rather than assuming Ba is the mean.
- **Input measures,** one per input:
  - *continuous:* uniform on an interval (by default, the training range in transformed units), or a finite mixture of truncated normal and Student t components;
  - *categorical:* supplied level probabilities, uniform over the observed levels by default. R ignores supplied weights and always uses uniform moments (defect fixture: a level with probability 0.9 gives a second moment of 0.5).
- **Mixture conventions:**
  - *Truncation.* By default the weights are those before truncation and the whole mixture is truncated and renormalized, as in R. Weights conditional on truncation are the other named form.
  - *Student components* take a location-scale `scale`, not a marginal standard deviation, as in R. A constructor from a marginal standard deviation converts it, for ν > 2.
- **Moments: one service for direct and basis sensitivity.**
  - Each term pair's moment is a product of one-dimensional interval moments of hinge powers under each input's measure, split at the knots and clipped to the measure's support. R's basis path ignores the bounds (defect fixture: 0.45 against the integral, 0.1106).
  - Categorical moments are sums of the supplied probabilities, including intersections, never ratios that can give 0/0. Francom et al. 2019 give the analytic categorical solutions, in EOF space too, in their supplement, which is not in hand.
  - Every hinge degree is supported. Uniform measures have closed forms. Normal and Student mixtures use closed forms where they are stable, and otherwise Gauss–Legendre quadrature split at the knots with a stated error bound. R restricts these to degree 1. Truncated Student moments exist for every ν, so the closed form's factors in ν − 1 fall back to quadrature at the special values.
- **Functional responses** have two forms (JSS §3.1):
  - integrate the functional variable as one more input, for a total decomposition, with its measure declared: uniform over its range by default, or the empirical output grid;
  - condition on one functional coordinate, for indices as functions of t.

  Both report main and total indices. R's functional path omits the totals.
- **Basis sensitivity.**
  - The variance of the reconstructed output includes cross-component covariances under the input measure; fitting the components independently does not remove them.
  - Indices are dimensionless, and variance contributions are in output units (y_s applied). R reports decomposition units (defect fixture).
  - Enumerating subsets up to `int.order` drives the cost. R's `nind` filters after computation; the port exposes a budget instead.
- **Zero total variance:** the indices are undefined and reported with a status, never 0 or 1. A small negative component variance from rounding is clamped only within a tolerance relative to the total variance, and reported beyond it.
- **Effect functions:** main and two-way effects for the declared measure and degree, or an explicit refusal. R hard-codes uniform expressions and degree 1 (defect fixture).
- **Independent checks:**
  - moments against quadrature and category enumeration, including higher degrees, truncated mixtures and Student components with small ν;
  - the Sobol identities: indices over all orders sum to 1, and totals are at least the main effects;
  - cross-component moments against Monte Carlo;
  - draw order and zero variance;
  - the R defect fixtures above.

## Presentation (stage 9)

The capabilities P09–P15 (scope section of the hub): trace, fit and residual, interval and sensitivity data computed in `Maths.Regression`, with rendering in adapters above maths; inspectable model metadata and posterior summaries, with formatting above the engine. Intervals follow D07: quantiles of the predictive mixture, not of per-draw interval endpoints. No rendering dependency enters maths.

## Room for extensions

Prep work in this segment, so that the [extensions](extensions-brief.md) plug in without reworking the port:
- **The decomposed surface is a declared input** (D08): the mean now, and quantile surfaces under asymmetric-Laplace noise (GBMARS) later.
- **One moment service over piecewise polynomials,** not hinges only (derivations review §7, item 7). BAPS's tensor and anisotropic splines and ABMRS's B-splines then get variance decompositions from the same code.
- **A basis contract** that admits nonorthogonal and robustly estimated bases (OQ16's spherical PCA), keeps a slot for basis-estimation uncertainty, and takes the truncation-residual model as a declared choice (D07's empirical default, the diagonal Gaussian option).
- **Term derivatives** from the design-kernel contract (fitting segment), which active-subspace analysis needs.
- **A response representation with a declared reconstruction map.** D07's predictive draw assumes reconstruction is linear, y = y_m + y_s ∘ (K w + e). Elastic emulation (extensions brief) represents a curve by an aligned amplitude basis and a warping function, and reconstructs by composition, which is nonlinear. The basis contract therefore takes a representation together with its reconstruction map, with the linear map as the port's only member.
