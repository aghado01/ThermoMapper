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
- A report appended to this brief: per statement, its module, status, axioms and each change from the posed form, then a **Methods** section that documents the logic of the proofs concisely, with exposition and the mathematical key points, like a dense methods paper. It takes the place of a "Doc changes owed" heading (the user's guidance of 2026-11-16, for every report on this work).

## Working rules

Those of the [first brief](time-frequency-statements-brief.md#working-rules).

## Report

### 2026-11-16 — S10 to S15 proved, nine modules promoted, the Q extension prepared

**Landed:**
- `5b08d17`: Part A, S10 to S15.
- `f0cc4e2`: Part B, the promotions.
- `9afb503`: Part C, the preparation brief.
- The commit that adds this report.

S16 was not attempted.

**Part A, per statement.** Every declaration named below prints only `propext`, `Classical.choice` and `Quot.sound`. The modules are in `lean/Enthymemata/TimeFrequency/` except S14's.

| | Module | Status | Changes from the posed form, and why |
|---|---|---|---|
| S10a | `LogPowerCovariance.lean` | proved | Stated at n + 1 over ℚ (`alternating_choose_harmonic`), with an ℝ copy. |
| S10b | `LogPowerCovariance.lean` | proved, no premise | **The moment formula** ∫ xᵏ ln x e^{−x} = k!(H_k − γ) is proved, not assumed (`integral_pow_mul_log_mul_exp`): it is Mathlib's Γ′(k + 1) through the Gamma integral's derivative. Integrability comes from Mathlib's Mellin derivative lemma (`integrableOn_pow_mul_log_mul_exp`). **The Laguerre polynomial** is as posed (`laguerre`). |
| S10c | `LogPowerCovariance.lean` | proved from a named premise | **The premise** `LaguerreProductPremise` is the Hille–Hardy expansion integrated term by term: E[ln X·ln Y] = Σₙ κⁿ·(∫ ln x·L_n·e^{−x})². The interchange of sum and integral belongs to the premise. **The marginals' means** enter as hypotheses: each is ∫ ln x e^{−x}, which `integral_log_mul_exp` shows is −γ. **The endpoint** π²/6 is `dilog_one`, from Mathlib's ζ(2). **The instance** meeting the premise is a degenerate law at κ = 0, not the Gaussian pair: it shows the premise is no contradiction. |
| S11 | `LogPowerCovariance.lean` | proved from a named premise | **The pairing rule** is a hypothesis on the one fourth moment used (`cov_normSq_pairing`). **Zero means** are not needed. **The instance** is a six-point real law with E x² = 1 and E x⁴ = 3, read as a complex pair. |
| S12 | `PeakOfExponentials.lean` | proved | **As posed** (`integral_one_sub_pow_one_sub_exp`). **The probabilistic form** (the maximum's mean) is not formalized: it needs the layer-cake formula and a product law. |
| S13 | `WindowTurnover.lean` | proved | **The tight bracket** 1.25 < u⋆ < 1.26 holds (`turnover_spec`, with `turnover_le` for the maximizer). |
| S14 | `Lemmas/BassWhitenedMarginal.lean` | proved | **Stated on second moments** with integrable products (`integral_secondMoment_response`, `marginal_covariance`). That the sum is Gaussian stays a premise of S4's reading. **The instance** is four outcomes with β ⟂ ε. |
| S15 | `QuadraticFormVariance.lean` | proved from a named premise | **The premise** `Pairing` is Isserlis on the fourth-moment tensor. **The results:** `pairing_variance` (algebra) and `variance_quadForm` (expectation). **The instance** is the six-point law. |
| S16 | — | not attempted | Optional and last; see the follow-ups. |

**Part B, the promotion review.** Every check below ran on each of the nine modules:
- no `sorry` and no `axiom`;
- no `import Enthymemata.*`;
- only barrels refer to it;
- a witness for each main result;
- a fresh read of each statement against its brief.

All nine pass, and all are moved with namespaces and declaration names unchanged.

| Module | Now at | Review note |
|---|---|---|
| S1 `PhasePropagation` | `Lemmas/TimeFrequency/` | none |
| S2 `WindowScaling` | `Lemmas/TimeFrequency/` | none |
| S3 `FlatnessBound` | `Lemmas/TimeFrequency/` | The docstring now points to S15 for the Gaussian identity it cites. |
| S6 `Contraction` | `Lemmas/TimeFrequency/` | none |
| S4 + Khatri–Rao + S14 `BassWhitenedMarginal` | `Lemmas/` | `gaussLogDensity` is the density's formula with log det, a density only for positive-definite C; the docstring says so. |
| S5 `BamrsKroneckerGram` | `Lemmas/` | `BAMRS.proj` duplicates `BASS.proj`, in another namespace. |
| S7 `DesignRank` | `Lemmas/Spline/` | `fullRank_iff_greedy` is conditional on `h461`, proved at order 1. |
| S8 `TrigonometricBSpline` | `Lemmas/Spline/` | none |
| S9 `FreeKnotSteps` | `Lemmas/Spline/` | none |

- *Left in `Enthymemata`:* the four noise-law modules, which landed today and await the coordinating session's clean-archive check.
- *Barrels:* `Lemmas/TimeFrequency.lean` and `Lemmas/Spline.lean` are new. `Enthymemata/Spline.lean` is removed, and `Enthymemata/TimeFrequency.lean` now holds the noise-law modules.
- *The review is the author's.* The coordinating session's clean-archive verification is the independent check.

**Part C.** [output-inner-product-brief.md](../regression/bass/output-inner-product-brief.md) is prepared, not built. It covers:
- the proved specification;
- what is already built;
- a table of every place that reads the independent form;
- seven decisions with recommendations and the slice where each binds;
- five slices gated by the whitening oracle, with Q = I bit for bit as the invariant.

Building waits on the user's go-ahead.

**Declarations to print axioms for:**
- `TimeFrequency.alternating_choose_harmonic`, `alternating_choose_harmonic_real`.
- `TimeFrequency.integral_pow_mul_log_mul_exp`, `integrableOn_pow_mul_log_mul_exp`, `integral_log_mul_laguerre`, `integral_log_mul_exp`.
- `TimeFrequency.logPower_covariance`, `dilog_one`, `cov_normSq_pairing`.
- `TimeFrequency.integral_one_sub_pow_one_sub_exp`, `turnover_spec`, `turnover_le`.
- `TimeFrequency.pairing_variance`, `variance_quadForm`.
- `BASS.integral_mulVec_mul_mulVec`, `integral_secondMoment_response`, `marginal_covariance`.

**Gates:**
- `lean` on each file during drafting.
- `lake build` with the linter set on `Enthymemata.TimeFrequency`, `Enthymemata.BassWhitenedMarginal` (before its move) and `Lemmas`: no warnings. One lint is disabled locally on `integral_secondMoment_response`, whose statement needs no finiteness of ι while its proof does.
- `scripts/meta-ci.ps1 -NoBuild` after the moves: PASS.
- Cache: scoped fetches of `Mathlib.NumberTheory.Harmonic.GammaDeriv`, `Mathlib.NumberTheory.ZetaValues`, `Mathlib.Analysis.Complex.ExponentialBounds` (95 files) and `Mathlib.Data.Matrix.ColumnRowPartitioned` (1 file). D: had 22 GB free at the end.

**Follow-ups:**
- *S16,* Theorem 4.61 at order 2. It would discharge `h461` for hat functions and make `fullRank_iff_greedy` unconditional for degree-1 splines. The route: a nonzero spline vanishing at n sites interleaved with its supports has too many zeros, Schumaker's Theorem 4.56 at m = 2.
- *The two premises with Gaussian content,* S10c's Hille–Hardy expansion with its interchange and S15's Isserlis rule. Both await Mathlib's multivariate Gaussian.
- *Clean-archive verification* of the four noise-law modules and of the promoted tier.

#### Methods

**S10a.**
- *Pascal's rule* splits C(n+1, j+1) into C(n, j) + C(n, j+1), and H_{j+1} = H_j + 1/(j + 1).
- *The pieces:* the C(n, j+1)·H_{j+1} sum is the row-n sum A_n read with its index shifted. The C(n, j)·H_j sum is −A_n. They cancel, leaving −Σ_j (−1)ʲ C(n, j)/(j + 1).
- *The leftover sum:* that one is 1/(n + 1), since C(n, j)/(j + 1) = C(n+1, j+1)/(n + 1) and the alternating row sum of C(n+1, ·) is 0. Hence A_{n+1} = −1/(n + 1), with no induction hypothesis on A_n.

**S10b.**
- *The moments.* Mathlib gives Γ′(k + 1) = k!(H_k − γ) (`Complex.hasDerivAt_Gamma_nat`). It also gives the Gamma integral's derivative as ∫ t^{s−1} ln t e^{−t} dt (`hasDerivAt_GammaIntegral`). Γ equals that integral on Re s > 0, so the two derivatives at k + 1 coincide. The complex integral is the real one, since t^{(k+1)−1} = tᵏ for t > 0.
- *Integrability:* the first half of `mellin_hasDerivAt_of_isBigO_rpow`, with its side conditions copied from Mathlib's own proof.
- *Then:* expanding L_n gives Σ_k C(n, k)(−1)ᵏ(H_k − γ) = −1/n − γ·0.

**S10c.**
- *The n = 0 term* of the premise's series is (∫ ln x e^{−x})², which is the product of the two means. Removing it (`hasSum_nat_add_iff'`) leaves Σ_{n≥1} κⁿ·(−1/n)².
- *No γ appears,* so the result holds for any common value of the marginal means.

**S11.**
- *Lift to ℂ:* ‖a‖² = a·conj a, so each real integral is the complex integral of a product (`integral_complex_ofReal`).
- *Conjugates:* ∫ conj a·b = conj ∫ a·conj b and ∫ conj a·conj b = conj ∫ a·b (`integral_conj`).
- *Then:* the pairing rule's two cross terms are \|E a·conj b\|² and \|E ab\|². No integrability is needed, because the premise supplies the values.

**S12.**
- *The geometric sum,* with q = 1 − e^{−x}: 1 − qⁿ = (1 − q)·Σ_{k<n} qᵏ = e^{−x}·Σ qᵏ. That is the derivative of Σ_{k<n} q^{k+1}/(k + 1).
- *The improper integral:* Mathlib's `integral_Ioi_of_hasDerivAt_of_nonneg` gives limit minus initial value, H_n − 0. It needs only a nonnegative integrand.

**S13.**
- *The derivative:* φ′(u) = (1 − e^{−u})·(e^{−u}(1 + 2u) − 1)/u², and the last factor is e^{−u}·ψ(u) with ψ = 1 + 2u − eᵘ.
- *ψ's shape:* ψ(0) = 0, ψ′ = 2 − eᵘ, so ψ rises on [0, ln 2] and falls on [ln 2, ∞). Hence ψ > 0 on (0, u⋆), ψ < 0 after, and the positive root is unique.
- *The bracket:* e^{5/4} = e·e^{1/4} < 2.7182818286 × 1.2841 < 3.5, and e^{63/50} = e·e^{13/50} > 2.7182818283 × 1.2969 > 3.52. Both use `Real.exp_bound` with five terms. ψ therefore changes sign in (5/4, 63/50), and the intermediate value theorem places u⋆ there.
- *Monotonicity:* from the derivative's sign on each side (`strictMonoOn_of_deriv_pos`).

**S14.**
- *A general lemma:* E[(AX)_i (AX)_j] = (AMAᵀ)_ij for any second-moment matrix M, by expanding the double sum and integrating term by term.
- *Stacking:* X = (β, ε) over κ ⊕ ι, with A = [Z \| I] and M block-diagonal (Σ_β, Σ_ε), because E[βεᵀ] = 0.
- *Then:* `fromCols_mul_fromBlocks` and `fromCols_mul_fromRows` give ZΣ_βZᵀ + Σ_ε. With S4's Σ_β and Σ_ε, the `module` tactic finishes σ²(Ω + τ⁻¹ZG⁻¹Zᵀ).

**S15.**
- *Expansion:* E[(xᵀHx)²] = Σ_pqrs H_pq H_rs m4_pqrs.
- *The three pairings:*
  - the direct pairing S_pq S_rs gives (Σ H_pq S_pq)² = (E xᵀHx)²;
  - the two cross pairings each give tr((HS)²) once the indices are relabelled and H, S used symmetric.
- *Reindexing:* tr((HS)²) expands as Σ_p Σ_a Σ_c Σ_b H_pb S_ba H_ac S_cp. Matching it to each cross pairing takes two or three swaps of nested sums, done under `conv`.
- *S symmetric* follows from S's definition as second moments.

**Part B.** The README's bar is a module that stops apologizing and whose statement and dependencies pass semantic review. Meta-CI enforces the first half; this report records the second, by module.
