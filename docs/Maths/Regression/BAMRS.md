# BAMRS — Bayesian Adaptive Multiresolution Splines

> Multiresolution regression spline basis under extended Bayesian information criterion (EBIC) count laws, developed from the ABMRS framework (He, Yang and Kang 2026). Part of `Maths.Regression`.

## Overview & Development Arc

BAMRS (Adaptive Bayesian Multiresolution Splines) extends the regression family in `Maths.Regression` to multi-scale hierarchical spline bases.

### Source Analysis & Corrections (He, Yang and Kang 2026)

From the audit of the ABMRS reference paper and code:

| Defect / Deviation in Source | What Our Check Found | The Engine's Design |
|---|---|---|
| BIC-type Normal marginal approximation | Approximates the marginal likelihood where Lemma 1 admits an exact closed form. | Keep the exact marginal likelihood without asymptotic BIC approximations. |
| Thin basis column pruning with 10⁻⁸ ridge | Thin basis columns dropped silently, breaking the partition of unity and altering the model. | Refuse unsupported configurations explicitly rather than silently altering basis columns. |
| Coefficient sampling at OLS $\hat{\sigma}$ | Fixes residual variance to point estimate, omitting $\sigma$'s posterior uncertainty. | Full conjugate draw: sample $\sigma^2$ conditionally, then $\beta \sim \mathcal{N}$. |
| Forced birth at $k=0$ | Carries proposal factor $b_0$ while reverse death omits $1/b_0$, undersampling $k=0$ when a knot is favoured. | Maintain proper reversible jump count balance at every knot configuration $k \ge 0$. |

### Knot Candidacy & Prior Representation

- **Count law:** `EbicPrior` reproduces the EBIC prior count law of He, Yang and Kang ($2\gamma \log C$).
- **Continuous vs Discrete Candidates:** While their original implementation restricts knots to subsets of a discrete candidate grid, `Maths.Regression` handles continuous knot locations on $(0, 1)$, preserving geometric invariance and avoiding grid artifact traps (keep $C \ge 2 k_{\max}$).
