# Several outputs on one basis: E34's build (chip brief)

Date: 2026-11-08. Status: commissioned by the user on 2026-11-08 ("Yes, start E34"), after E8 was completed; slice 0 starts. Written by the main session from [E34's entry](extensions-brief.md#several-outputs-on-one-basis-with-σ-integrated-e34), which holds the model, its pieces and the checks. Reports go under "Reports" below, newest first.

## Why

K continuous outputs per case that share one structure: several responses of one experiment, unordered or few enough that a free K × K covariance is affordable. One fit selects terms for all outputs at once, reads their residual correlation, and gives joint predictions:
- joint bands;
- the cross-output correlation of predictions;
- the probability that every output meets its specification.

It serves any multi-response surrogate or regression, and E22 (coupled outcomes) builds on it. Ordered outputs on a grid stay with the direct and PCA modes.

## Sources

- **The model** is multivariate Bayesian regression with a matrix-normal prior scaled by Σ and an inverse-Wishart on Σ. Brown, Vannucci and Fearn (JRSS B, 1998) select variables in it with inclusion shared across responses; Huang and Wand (*Bayesian Analysis*, 2013) give the scale-mixture alternative to the inverse-Wishart. **Neither is on disk.** Ask the main session for them in slice 0. Until then, cite them from memory and say so.
- **No oracle** is known in CRAN's lineage. The checks are internal and against closed forms, as the entry lists.

## What is already built

- **The collapsed core:** `CollapsedGPrior`, `GPriorConditional`, `PriorPencil`. Its coefficient prior runs through the prior-precision builder (`BassPriorPrecision`: the g-form, the normalized ridge, KHAOS's raw ridge and modified g-prior), and σ² is integrated under `Sigma2Prior` (`ScaledToResponse` by default, D03).
- **Responses:** `BassData` and `BassResponse` (`ScalarResponse`, and a functional response over shared output locations), plus basis fits through `BassDecomposition` and `BassComponents` (D07: PCA, each component fitted apart).
- **The structure chain and its samplers:** the ladder with `RungWorkers`, SMC with its two flags, the flat histogram, resume.
- **GBASS's noises:** the latent record per observation (`ObservationLatents`), weights and offset on the working response, and the move along β·v as a partially collapsed step (`1a8459c`: the coefficients redrawn after an accepted move).
- **Random intercepts** (E21, `0a2168c`): the rows whitened by V = I + s²ZZᵀ, with the effects drawn from their conditional.
- **Readouts:** the moment service and Sobol; the gradient readout; E8's co-active posterior for two outputs of one fit (`CoActivePosterior`).

## Decisions

These are the main session's, taken in commissioning. Each goes on the page with its alternatives and trade-offs. The user may reverse any of them.

1. **The model.** Y (n × K) = 1μᵀ + BA + E, with each row of E ~ N(0, Σ).
   - *The priors:* A | Σ, τ ~ MN(0, S₀(τ), Σ) through the prior-precision builder, and Σ ~ IW(ν₀, Ψ₀). Inclusion is shared across outputs.
   - *The collapse:* with A and Σ integrated, the structure step reads a K × K quadratic form and its log-determinant, and the marginal is matrix-t. The intercept is level-centred per output, with n − 1 degrees of freedom.
2. **Σ's default prior: D03's `ScaledToResponse`, carried to K outputs.**
   - *The form:* ν₀ = K − 1 + 2c and Ψ₀ = 2c·S_Y, with c = 10⁻³ and S_Y the outputs' sample covariance about their means, stated as an empirical prior.
   - *At K = 1* it is IG(c, c·s²), the scalar default exactly. That reduction is a gate.
   - *Units:* because S_Y → LᵀS_YL when Y → YL, the structure posterior is invariant to any invertible rescaling or mixing of the outputs. This settles the entry's output-units choice.
   - *Options:* a declared Ψ₀ in the outputs' units, and standardized outputs (scale-invariant only).
   - *Refused:* an S_Y that is singular or nearly so, with the reason. Collinear outputs need fewer outputs.
3. **Properness under tempering, correcting the entry.** At rung b the tempered target integrates p(Y | A, Σ)^b against the priors (D05's integrated form), so Σ's conditional is IW(ν₀ + b·ν_data, Ψ₀ + b·R), with ν_data = n − 1 under level centring.
   - *When it is proper:* whenever ν₀ + b·ν_data > K − 1. Under a proper prior (ν₀ > K − 1) that holds at every b ≥ 0, including the reference rung.
   - *What changes from the entry:* it asked the tilted prior alone to be proper (ν₀ − n(1 − b) > K − 1), which is stricter than needed. Refuse only a rung, or a reference rung, where the condition above fails, which happens only under an improper ν₀.
4. **Shared inclusion only.** Inclusion per output needs an indicator per (term, output) and loses the matrix-t.
5. **Missing outputs are refused** at first, with the reason.
6. **GBASS's noise by row.** A per-row scale vᵢ is shared across a case's outputs, with a skew vector γ ∈ ℝᴷ, giving the multivariate t and the multivariate generalized hyperbolic.
   - *The γ·v move* is the β·v move's analogue: a partially collapsed step, with A and Σ drawn afresh after an accepted move, and a proposal that reads the state alone.
   - *What carries over:* the skewed horseshoe's default move and its recorded hot-rung limit.
7. **Random intercepts with one s shared across outputs.** Group effects u_g ∈ ℝᴷ ~ N(0, s²Σ), so the row covariance is V ⊗ Σ, the matrix-normal form holds, and the effects' correlation is Σ's. An s per output breaks the form and is not built.
8. **The PCA mode's joint option.** D07's retained scores are fitted as K outputs here, with truncation residuals as now. It depends only on the scores' span (decision 2), and trades per-component sparsity for that.

## The slices

One commit per concern; tests with the code they test; each slice ends with its gates and a fresh-eyes review before the next starts.

| Slice | Delivers | Its own gates |
| --- | --- | --- |
| 0. References | The matrix-t marginal in closed form; the sources read once supplied | Against the dense vec form with an asymmetric Σ at K = 2, so Σ ⊗ S₀ and S₀ ⊗ Σ cannot pass for each other under the code's vec order; against quadrature over Σ at K = 2 and Monte Carlo at K = 3; at K = 1 against the scalar marginal under decision 2's reduction |
| 1. The core | A multi-output response; Σ's prior (decision 2) and its validation; the collapse in the structure step and τ's step through `PriorPencil`; per draw Σ | Y and A | Σ; level-centred intercepts | K = 1 equal to scalar BASS to rounding, with scalar BASS's pinned hashes unmoved; the law by enumeration on small structure spaces; Geweke; equivariance Y → YL; recovery of a known Σ and of surfaces from outside the model |
| 2. Samplers and evidence | The ladder under decision 3, SMC with its flags, the flat histogram, resume | Each rung's law and the evidence against the enumeration; resume bit for bit; a refusal only where decision 3's condition fails |
| 3. Readouts | The matrix-t predictive; joint bands; the predictions' cross-output correlation; P(every output in its specification); E[Σ] and its correlation matrix; Sobol per output and for declared linear combinations; co-active pairs across outputs through E8 | Against Monte Carlo over the draws; the joint probability against counting; Sobol per output against single-output fits on the same structure |
| 4. Noise and random intercepts | Decisions 6 and 7 | Geweke under each noise and with the block; the γ·v move on a frozen state against its orbit's law; recovery from outside the model |
| 5. The PCA mode's joint option | Decision 8 | Against separate component fits on data whose scores are independent; the equivariance gate on the scores |

## Gates common to every slice

- **Scalar BASS untouched.** Its pinned hashes and fingerprints hold: the pinned tier at each commit, the whole assembly at the slice's end ([targeted runs](../../../CONTRIBUTING.md#testing-standards)).
- **Truth from outside the model,** on surfaces and covariances not drawn from the engine's priors.
- **Primitives:** chosen cases and a seeded random sweep, far ranges and near-singular Σ included.
- **The cost row:** the structure step at K = 1, 2, 5 and 10 against scalar BASS, and the readouts. Cycle clock, marked indicative.
- **`docs/Maths`:** the decisions on the regression page or a new page, with their alternatives and implications. New writing takes the module's symbols, and crossings go in `issues/notation-ledger.md`.

## Not in this brief

- **E22's coupled outcomes** and its parameter expansion for a correlation matrix.
- **Huang and Wand's prior** on Σ; inclusion per output; missing outputs as latents.
- **The response envelope** (E34's entry), and an output basis from E[J W_x Jᵀ].
- **Conformal prediction for several outputs.** The candidate seam refuses a multi-output fit, with its reason.

## Working rules

- On `main`, committing by pathspec with nothing left staged and nothing pushed. Commit messages describe the work, with no attribution lines ([CLAUDE.md](../../../CLAUDE.md)). Leave `basin-hold-validation.csv` and `path-profile-validation.csv` out of every commit.
- Other sessions may be writing in the checkout: the conformal chip on E35, the GBASS chip on a fix to `NormalVarianceMeanMixture`. Before editing a file, check `git diff -- <file>` for changes you did not make; build from a snapshot (`git archive HEAD` plus your files); commit only compiling states. Shared docs (`docs/Maths/BayesianRegression.md`, the notation ledger, `extensions-brief.md`) are edited and committed in one step.
- Builds set `ArtifactsRoot` inline into `artifacts\chip-e34\`, deleted at the end.
- Heavy runs at below-normal priority with capped threads, MSBuild node reuse off; the machine is shared. Tests run by [targeted runs](../../../CONTRIBUTING.md#testing-standards) (`scripts/test-select.ps1`, `-SourceRoot` at your snapshot): the selected classes per edit, with the pinned tier per commit, the whole assembly once at the slice's end.
- **Sources:** read a paper before stating what it says, and mark anything cited from memory. When a source is not on disk, ask the main session, which asks the user. Fetch nothing yourself.
- New beside old; a removal is its own question to the user, with what is lost.
- Findings framed by problem class, not by any one application.
- A follow-up goes in a report with the decision it would change and what each outcome would do.
- After each slice, report to the main session and wait for its verification note before starting the next.

## Report

Under "Reports" below, newest first, per slice: what landed by commit; the gates' numbers; decisions taken beyond this brief; what was left unbuilt and why; follow-ups with the decision each would change.

## Reports
