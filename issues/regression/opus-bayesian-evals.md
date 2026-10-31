
## The primitives

1. **Exact readouts are roots and integrals of piecewise polynomials.** BARS's peak and span readouts are one case of something more general. Any feature defined by a polynomial equation in x, f, f′ and f″ has an exact posterior per draw. So do features defined between two fits, since f − g on merged knots is still piecewise polynomial. That covers:
   - peaks, onsets and inflection points (f″ = 0);
   - time above a threshold and area under the curve;
   - the point where two curves cross;
   - break-even points;
   - any condition that clears denominators into a polynomial. For example, revenue p·D(p) peaks where elasticity equals −1, which is the polynomial condition D + pD′ = 0.
2. **Variance decomposition with uncertainty.** Sobol indices and active subspaces work on any tensor-product representation, including a surrogate fitted to a black box.
3. **A noise family.** GBMARS gives heavy tails, skew and quantiles.
4. **Calibration and comparison.** Conformal sets over retained draws, LOO, evidence and proper scores sit on top of every engine. That makes "does the conclusion survive a change of prior" a routine check.
5. **Curves as data.** Functional outputs, elastic alignment, whole-curve bands, and replicates as exchangeable units.
6. **Sequential design aimed at features.** You sample where the feature is uncertain, rather than where the function is. Bayesian optimization is one special case.

