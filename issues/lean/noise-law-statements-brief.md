# Noise-law statements, open links, and the Q extension's preparation: the Lean chip's third brief

Date: 2026-11-16. Status: **posed at the user's request on 2026-11-16** ("do we have follow up work for lean chip. im also ok with additive engine work so long as its coherent"). Drafted by the coordinating session. It follows the [time–frequency](time-frequency-statements-brief.md) and [spline-space](spline-space-statements-brief.md) briefs, both closed.

## Why

The numerical chip's sweeps confirmed two formulas by measurement and corrected two of the coordinating session's claims; each now has an exact statement worth proving. The first two reports also left three links assumed. And the whitened marginal (S4) is the specification of an engine extension the working-response contract already declares, so the chip that proved it is the one best placed to prepare that extension.

## Sources

- **The numerical report:** [latent-recovery-pilot-brief.md](../regression/latent-recovery-pilot-brief.md), its Readings for B1 and B2 and its Primitives proposed.
- **Your own reports,** for the links left assumed: the Gaussian identity behind S3's reading, the marginal covariance S4 takes as given, and S7e's premise.
- **Schumaker,** as before, with `errata.md` for the corrected forms: § 4.8 for S16.
- **The engine,** for Part C: `src/maths/regression/spline/bass/BassObservationModel.cs` (the contract and its comment on Q), `BassDesign.cs`, `BassTarget.cs`, `BassPointwise.cs`, `posterior/`, `calibration/CalibrationCovariance.cs`, and `src/maths/regression/gprior/`.

## Part A: statements

Numbering continues. Each is the coordinating session's posing; the formulas from memory are marked, and a statement false as posed is a result.

**S10. The dilogarithm covariance** (B1 read it to within 1.6 standard errors at fifteen κ).
- *S10a.* For n ≥ 1, Σ_{k=0}^{n} (−1)ᵏ·C(n, k)·H_k = −1/n, with H_k the harmonic number and H₀ = 0.
- *S10b.* With the Laguerre polynomial L_n(x) = Σ_k C(n, k)·(−x)ᵏ/k!, ∫₀^∞ ln x · L_n(x)·e^{−x} dx = −1/n for n ≥ 1. It follows from S10a, Σ_k (−1)ᵏ·C(n, k) = 0, and ∫₀^∞ xᵏ·ln x·e^{−x} dx = k!·(H_k − γ).
- *S10c, with a named premise.* Premise, from memory: the pair (|a|², |b|²) of unit-variance jointly circular complex Gaussians with squared coherence κ < 1 has joint density e^{−x−y}·Σ_{n≥0} κⁿ·L_n(x)·L_n(y). Then the covariance of ln|a|² and ln|b|² is Σ_{n≥1} κⁿ/n². Close with the endpoint: the series at κ = 1 is π²/6.

**S11. Powers of a non-circular pair** (B1's near-DC reading). For zero-mean complex a and b whose fourth moments obey the Gaussian pairing rule, Cov(|a|², |b|²) = |E[a·conj b]|² + |E[a·b]|². The pairing rule (Isserlis) enters as a named premise unless Mathlib supplies it. The second term is the pseudo-Gram entry the numerical chip found necessary within two bins of DC.

**S12. The peak of n independent unit exponentials** (B2: a single frame's largest bin is not flat for noise). ∫₀^∞ (1 − (1 − e^{−x})ⁿ) dx = H_n. If it comes cheaply, also the probabilistic form: the expected maximum is H_n.

**S13. The rectangular window's turnover** (B2 read 1.255·τ). φ(u) = (1 − e^{−u})²/u on u > 0 increases to a unique maximizer u⋆ and decreases after it; u⋆ is the unique positive root of eᵘ = 1 + 2u; and 1.25 < u⋆ < 1.26. A wider bracket is acceptable if the tight one costs too much.

**S14. The marginal covariance behind S4.** For β with mean zero and covariance σ²τ⁻¹(ZᵀΩ⁻¹Z)⁻¹, and ε uncorrelated with β with covariance σ²Ω, the covariance of Zβ + ε is σ²(Ω + τ⁻¹Z(ZᵀΩ⁻¹Z)⁻¹Zᵀ). State it on second moments; that the sum is Gaussian is a premise unless Mathlib's multivariate Gaussian closes it.

**S15. The variance of a Gaussian quadratic form** (closing S3's reading). For a real zero-mean vector x with covariance Σ whose fourth moments obey the pairing rule, and symmetric H, Var(xᵀHx) = 2·tr((HΣ)²).

**S16, optional and last. Theorem 4.61 at order 2.** Discharge S7e's premise for hat functions, so the rank-from-greedy equivalence holds unconditionally for degree-1 splines.

## Part B: promotion review

Run the README's promotion checks on the nine closed modules and report each one's eligibility and anything it still needs. Promote where the README lets a session do so; list the rest for the user.

## Part C: prepare the output-inner-product extension

The user has said additive engine work is acceptable so long as it is coherent. This part prepares that work and stops before building.

- **The extension.** BASS's core serves the output inner product Q that the working-response contract declares for a direct functional response: residuals correlated along the output index, covariance σ²(I ⊗ R), Q = R⁻¹. Today the core refuses Q ≠ I.
- **Its specification** is yours: `gPrior_logMarginal_whiten`, the Khatri–Rao Gram under I ⊗ Q, the level under Q, and S14 once proved.
- **Write the preparation brief** at `issues/regression/bass/output-inner-product-brief.md`, in the form the other BASS briefs take:
  - what the contract's comment says is missing, and every place in the core, the target, the pointwise contents, prediction and calibration that reads the independent form today;
  - the decisions, each with a recommendation and the slice where it first binds: where the coefficient prior's Gram reads under Q, how Q combines with observation weights and with several outputs, what tempering and the evidence need, what the residual-dependence check reports once Q is declared;
  - the slices, each with its gate;
  - the oracle: a fit under Q must equal the fit of the whitened response and design under Q = I, which S4 proves;
  - the invariant: with Q = I every existing result stays bit for bit.
- **Read the numerical chip's primitives 1 and 2** first; they are the two callers that asked for it.
- **Stop there.** Building starts on the user's go-ahead, which the coordinating session relays.

## Scope and order

- **Order:** S10a, S12, S13, S14 and S15 first; then S10b, S10c and S11; Part B; Part C; S16 last.
- **A witness for every statement** and `#print axioms` for every closed result, as before. A named premise gets an instance that meets it wherever one is cheap.
- **Part C writes one brief** under `issues/regression/bass/` and nothing under `src/`, `tests/` or `docs/`.

## Done

- S10a, S12, S13, S14 and S15 proved without `sorry`; S10b, S10c and S11 proved or compiling with their premises named.
- Part B's eligibility table, with promotions made where allowed.
- The preparation brief written, self-reviewed and committed.
- A report appended to this brief in the earlier briefs' form.

## Working rules

Those of the [first brief](time-frequency-statements-brief.md#working-rules).

## Report

### In progress

*Updated 2026-11-16.*
- **Landed:** Part A's S10 to S15, each with a witness and standard axioms only.
  - S10a, S10b, S12, S13, S14 and S15 are proved.
  - S10b needs no premise: Γ′(k + 1) = k!(H_k − γ) comes from Mathlib, and integrability from its Mellin derivative.
  - S10c and S11 are proved from their named premises, each met by an instance.
  - S13 holds with the tight bracket 1.25 < u⋆ < 1.26.
- **Next:** Part B (promotion review), Part C (the preparation brief), then S16.
- **Uncommitted:** nothing.
