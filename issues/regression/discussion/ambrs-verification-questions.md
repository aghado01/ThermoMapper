Verify in the Biometrics text first. The paper is not yet in the compendium; the preprint is, as HYK2024.md. The external summary of the differences was wrong on one checkable point: it said the death formula changed, but the preprint already has the same formula. So check each item against the paper:

#	What to verify	Where it matters
1	Whether the prior keeps the preprint's form: π(k) ∝ C(C,k)^(1−γ) with π(ξ | k) = 1/C(C,k) (equation 3)	The γ_HYK = 1 + γ mapping, and so this step
2	γ's range (γ ≥ 0?), the three regimes and the recommended [0.5, 2]	Whether γ > 1 is their convention, to cite, or our extension, to label; constructor bounds; docs
3	Guidance on C ("typically > m"), or anything on k past C/2	Support [0, C] with guidance, or truncation
4	The b_k and d_k formulas; c = 0.3 and any stated reason or mixing claim	The mixing benchmark's arms (BARS mixing)
5	Coincident knots (multiplicity above 1)	The tensor-product extension's handling of jumps (port brief)
6	The non-normal treatment (still an EBIC plug-in?)	The docs' comparison with our Laplace marginals
7	How peak inference is defined and computed (draws or a plug-in)	The readout contract for the tensor-product extension
8	Whether the γ comparison still favors γ = 1, and whether it separates the prior from mixing	The confound noted under BARS mixing


9	Title, method name (ABMRS?), volume, issue, article ID, DOI	The citation
Answer: 

Junhui He, Ying Yang, Jian Kang, Adaptive Bayesian multivariate spline knot inference with prior specifications on model complexity, Biometrics, Volume 82, Issue 3, September 2026, ujag121, https://doi.org/10.1093/biomtc/ujag121
Method name: ABMRS (Adaptive Bayesian Multivariate Regression Splines)