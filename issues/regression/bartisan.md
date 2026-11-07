
Most BART software handles only Gaussian and binary outcomes. Using the sampler of Linero (2025), which goes beyond the conditionally conjugate case, bartisan supports a much wider range of response distributions: Poisson, negative binomial, gamma, beta, ordered beta, ordinal, multinomial, zero-inflated counts, Tweedie, AFT and proportional hazards models for survival data, location-scale regression, a Dirichlet process mixture for the errors, and custom families you write yourself. A model is specified just as in glm(), with a formula, data, and a family.

My goal with bartisan was to replace parametric models, such as those implemented in rstanarm, by combining the smoothness and flexibility of BART with the modern framework for computing interpretable quantities described by Julia Rohrer & Vincent Arel-Bundock (2026). No longer do you need to fit unjustifiably linear models, and no longer do you need to sacrifice inference or interpretability to fit a flexible, data-driven one.

bartisan works directly with bayesplot, posterior, loo, and marginaleffects for model checking, inference, and interpretation, and has its own tools for convergence diagnostics, variable importance, and causal effect estimation. It also includes several recent extensions of BART:

• Varying coefficient models, including Bayesian causal forests
• Sparsity-inducing priors for variable selection
• Smooth trees (as in SoftBART), on by default
• Random intercepts

The package site has documentation and tutorials to help you make the transition from parametric models to BART, including best practices for checking convergence and fit and choosing the right outcome distribution: https://lnkd.in/gVvgigXu
