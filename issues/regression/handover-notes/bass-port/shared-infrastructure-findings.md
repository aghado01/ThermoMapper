# Shared-infrastructure chip: findings for the D03/D04/D05 drafts (2026-10-27)

From the chip running `issues/regression/bass/shared-infrastructure-brief.md`. Full report: `issues/regression/bass/port-numerics-checklist.md`, "Ladder evidence (N11) and the g-prior's mean-free marginal (N12)". Commits: `0afe33d` (N12), `84c173f` (N11), `a315bd9` (checklist).

- **D04's τ slice and D05's swaps:** `block.MarginalAt(τ, b).LogMarginal(prior)`. It returns a `GPriorMarginal` value and allocates nothing. It gives the same double as `At(τ, b).LogMarginal(prior)`, because `At` is built on it.
- **D03/D05 evidence:** `Maths.Samplers.Evidence.LadderEvidence`. The default is stepping-stone, which reads only ℓ_b at the rungs.
  - **Path sampling on a collapsed (b-nonlinear) path.** The integrand is ∂_b ℓ_b: `ITargetEvaluation.LogLikelihoodSlope(b)`, whose default refuses. Friel–Hurn–Wyse's discretization diagnostic reads Var_b(slope) as the mean slope's b-derivative, and so misses E_b[∂²ℓ_b/∂b²].
    - On the test's nonlinear enumerable target it read 0.094 against a true trapezoid bias of 0.669, while stepping-stone's bias was −0.004.
    - So if D05 keeps path sampling with a diagnostic, BASS must also supply the curvature, which has no seam yet. For the collapsed σ² it is Var over σ² | x, b of ln p(y | x, σ²), a closed form with trigamma terms.
  - **Normal BARS** could report a stepping-stone evidence today under its proper σ² default, at O(1) per rung via `MarginalAt`. Whether it should is a consumer decision, left open.
- **S3 (N13, OQ9) is done.** The user approved the download, and `posterior` 1.7.0 and `loo` 2.10.1 are in `r/renv.lock` (`5eecbc1`).
  - `diagnostics_oracle.R` and `DiagnosticsParityTests` agree to 10⁻⁹ relative.
  - One convention differs and is recorded: the tail median and quantiles over all draws versus split chains, which matters at odd lengths only.
  - OQ9 is closed. D10's external-reference gate is therefore met for R̂, ESS, k̂, r_eff and the MCSE.
- **Follow-ups the user asked for (2026-10-27).** Commits `128fd1f`, `ae82374`, `00eb9da` and `3cd4e98`.
  - *Pareto guard:* posterior's `gpdfit` guard is in the C#.
  - *Tail states:* `PredictiveCheck.ParetoTails` reports each unit's treatment. The tied tail keeps −∞, a divergence from loo that is documented and asserted.
  - *Odd lengths:* aligned to posterior.
  - *SMC stage traces:* `SmcOptions.TraceStages` and `LadderEvidence.FromSmc`; the Kish error under-covers when mutations mix poorly.
  - *BARS evidence:* `BarsResult.Evidence` is on by default for tempered runs that end at a reference rung under a proper σ² prior, and off for the default untempered run.
  - *For D05:* path sampling for BARS and BASS still needs ∂_b ℓ_b, and its diagnostic needs E_b[∂²ℓ_b].

## D05 against the seam: inputs for the user's design pass (2026-10-27)

The user asked for these to be recorded here; the chip does not edit the segment briefs. Read against `port-fitting-brief.md`, D05 and D03's evidence rule.

- **Where they agree.**
  - The reporting rule is already the seam's gate. Default runs are untempered and report none; evidence comes only from proper-mode runs with an exact b = 0 rung. BARS's `ReportsEvidence` is the same rule.
  - D05's O(1) ℓ_b for moves, swaps and the τ slice is `CollapsedGPrior.MarginalAt(τ, b)`.
  - Stage 8's normalizer is `LadderEvidence.FromSmc`, which reproduces `LogNormalizerRatio` to 10⁻¹⁰.
- **1. The estimator: D05 says "Evidence: path sampling … with ∂_b ℓ_b in closed form (§3)".**
  - On the integrated form p_b, the seam's measurements favour stepping-stone. It reads only ℓ_b, which BASS computes for every swap anyway.
  - Path sampling there needs two closed forms: ∂_b ℓ_b, and E_b[∂²ℓ_b] for its discretization diagnostic. Without the second, the diagnostic is blind: on the nonlinear test target it read 0.094 against a true trapezoid bias of 0.669, where stepping-stone's bias was −0.004.
  - Suggested wording: stepping-stone on the integrated form by default; path sampling there only with the curvature term.
- **2. The linear form (D05's 2026-10-27 proposal) reverses that.**
  - ℓ_b = b·log p₁ is linear in b, so the slope is log p₁ itself, already computed. The curvature is zero, so the Friel–Hurn–Wyse diagnostic reads the integrand's true derivative.
  - Path sampling is then free, its diagnostic holds, and it cross-checks stepping-stone over the same draws.
  - So the estimator follows the form: stepping-stone on p_b, both on the linear form.
  - The linear form still needs the proper mode for the evidence itself, since the improper constants sit inside log p₁. D05 says so already.
  - *Rehearsed on BARS* (`TemperedForm.Linear`, commit in the checklist's next report). One knot, 10 levels: stepping-stone −1.6313 ± 0.0047 and path sampling −1.6322 ± 0.0049, against a quadrature of −1.6231, with a diagnostic of 0.0010.
- **3. The evidence carries its priors.** D05 asks for the evidence to be "carried with its Monte Carlo error and with the priors and transforms it depends on".
  - `BarsEvidence.Priors` (`BarsEvidencePriors`) does that for BARS: the count prior, the support gate, the centre and its resolved level, and the σ² prior as specified and as resolved.
  - BASS's equivalent would add τ's prior, the term measure, the count prior on M and the interaction cap.
  - Both of BARS's resolved pieces read the response: the median centre, and s² under the scaled σ² prior. So the default evidence is empirical-Bayes. Comparisons are licensed only between configurations that share those pieces, on the same data.
- **4. SMC error at D05's mutation rate.** D05 mutates with one D04 sweep per stage. That is the regime where the Kish-based error under-covered: 0.795 at one move per stage, against 0.935 at five.
  - Before stage 8 reports an error with its normalizer, it needs either several sweeps per stage or the genealogy variance of Lee and Whiteley (2018), which is not built.
  - It isn't urgent until stage 8.
