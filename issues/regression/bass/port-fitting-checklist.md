# BASS port: fitting: checklist

Companion to [port-fitting-brief.md](port-fitting-brief.md). An item is ticked when its report lands. Reports go under "Reports", newest first, each with its date, commits, evidence and remaining limitations.

## Checklist

### Decisions (stage 1)

- [ ] **F-D01** Term representation and prior measure: drafted 2026-10-25 as a proposal; settles at the user's design pass.
- [ ] **F-D02** Count prior at the capacity limit: drafted 2026-10-25 as a proposal; settles at the user's design pass.
- [ ] **F-D03** Coefficient-prior modes and the τ hyperprior: drafted 2026-10-25 as a proposal; settles at the user's design pass. The intercept convention is the user's call.
- [ ] **F-D04** Collapsed and conditional update order: drafted 2026-10-25 as a proposal; settles at the user's design pass.
- [ ] **F-D05** Tempering and SMC target paths: drafted 2026-10-25 as a proposal; settles at the user's design pass.
- [ ] **F-D06** Numerical model admission: drafted 2026-10-25 as a proposal; settles at the user's design pass.
- [ ] **F-R1** Deterministic-simulator profile: drafted 2026-10-26 as a proposal; its coefficient prior is set by a pre-registered comparison.
- [ ] **F-M1** Move schedule: equal over available moves (D01) against the prior-cancelling option, in a pre-registered comparison.

### Reference fixtures (stage 1)

- [x] **F-O1** Propriety of the collapsed marginal over τ: `r/oracles/bass_marginal_oracle.R`.
- [x] **F-O2** Order-free admission and R's `getQf` rule: `r/oracles/bass_admission_oracle.R`.
- [ ] **F-O3** Count and order laws under a constant likelihood, by enumeration (D01, D02).
- [ ] **F-O4** Boundary ratios π(1)/π(0) and π(M_max)/π(M_max − 1), and R's fixed thirds as a defect fixture.
- [x] **F-O5** Wallenius probabilities against brute-force sums over draw orders: `r/oracles/bass_wallenius_oracle.R`.
- [x] **F-O6** Defect fixtures: R's tempering and swap density, R's Jeffreys mode, KHAOS's double count factor: `r/oracles/bass_enumeration_oracle.R` (modes `swap` and `count`) and `r/oracles/bass_jeffreys_oracle.R`.
- [x] **F-O7** Centred admission angles for the Liang default, as the `centred` and `hadamard` modes of `r/oracles/bass_admission_oracle.R`. Checks: uncentred s_j = centred s_j · ρ_j; the downdate's accuracy on a near-constant column; and the centred Hadamard form against explicit centring.

### Build (stage 4)

- [ ] **F1** Term representation: continuous, categorical and functional parts, with interned terms.
- [ ] **F2** Design layer and the direct functional contractions.
- [ ] **F3** Collapsed g-prior target under both intercept conventions, with the fixed-g mode.
- [ ] **F4** Birth, death and change moves with exact boundary probabilities.
- [ ] **F5** D04's sweep, with composition draws at retention.
- [ ] **F6** Parallel tempering on the likelihood path, with the exact reference rung in the proper mode.
- [ ] **F7** D06's two gates and refusal reasons.
- [ ] **F8** Retention, checkpoints and exact restart.
- [ ] **F9** Pointwise density producer on the OQ17 seam, with its chains for PSIS's relative efficiency.
- [ ] **F10** Chains start on the support.
- [ ] **F11** Direct fits report the residual-dependence check along the output coordinates (proposal).

### Room for extensions

- [ ] **F-X1** Observation-model seam carries an offset with the weights.
- [ ] **F-X2** The sweep has a latent slot.
- [ ] **F-X3** Term types behind one design-kernel contract, with derivatives.
- [ ] **F-X4** Proposal kernels behind one contract with exact forward and reverse probabilities.
- [ ] **F-X5** Coefficient priors as a declared set.
- [ ] **F-X6** The output-side inner product is declared (Q = I now).

## Reports

### The parity rule — 2026-10-27

The user's decision: the C# preserves no R defect to align with a fixture (the hub's validation section). The `BassLiteral` options are gone from D03 (R's τ rate 2/n) and D06 (R's `getQf` admission), and from the declared set of coefficient priors. R's zero-centred intercept stays selectable as DiMatteo, Genovese and Kass's prior, `GPriorCentre.Zero`, not as a parity option. Earlier reports that mention `BassLiteral` record the drafts as they stood.

### F-O6 at three predictors, and reversibility — 2026-10-27

From the reference-fixtures chip's follow-up, `f7acf2a`. The enumeration oracle now reports each kernel's detailed balance with its own stationary law. The checks add p = 3 to the swap and count fixtures, so the Wallenius factor of R's proposal varies with the state. All checks pass.
- **Swap at p = 3.** N = 27 on a 3 × 3 × 3 grid, a support rule of 4 nonzero values: 126 terms, 48 supported, 1174 states, 14,160 moving transitions.
  - *Validity.* R's `getQf` agrees with the strict rank on every state, and R's `alpha` is b times its b = 1 value. At b = 1 R's kernel, with the state's move probabilities, is the posterior to 7·10⁻¹⁵. With its fixed thirds it is 0.100 off in total variation (F-O4's defect, larger here than at p = 2).
  - *Targets.* At b = 0.5 and 0.2, R's kernel is 0.077 and 0.143 in total variation from D05's likelihood path, 0.055 and 0.081 from the prior times p_1^b, and 0.41 to 0.65 from the two powered paths.
  - *Swaps.* R's swap misses balance against its own rung laws by spreads of 1.88 and 1.14, and against the candidate families by 0.88 to 6.18. D05's kernels reach their target to 7·10⁻¹⁶, and its swap balances on them to 1·10⁻¹⁴.
- **R's tempered kernel is not reversible at p = 3.** Its gap from detailed balance with its own stationary law is 5·10⁻⁷ at b = 0.5 and 1·10⁻³ at b = 0.2, against 2·10⁻¹⁶ at b = 1 and at most 7·10⁻¹⁶ for D05's kernels. At p = 2, where every interaction draws the one pair, R's kernel stays reversible at every rung (at most 5·10⁻¹⁶).
  - *Why.* R's acceptance is the untempered ratio to the power b, so the proposal ratio enters as (q_rev/q_fwd)^b. The kernel then has a detailed-balance target only if that factor's (1 − b) power is a ratio of state functions. At p = 3 the Wallenius factor's is not.
  - *So below b = 1* R's within-rung moves balance against no law at all, before any swap. That sharpens D05's "no single target".
- **Count law at p = 3.** A constant likelihood, two knots per variable, M_max = 2: 60 terms and 1831 states. D01's ratio gives the stated law (0.202, 0.404, 0.394) to 6·10⁻¹⁴. KHAOS's gives the stated law times 1/M! (0.252, 0.503, 0.245) to 4·10⁻¹⁴, 0.149 from the stated law.
- **Not covered:** p ≥ 4, and K_max = 3, where the Wallenius factor takes three-variable sets. Both grow the state space past dense GTH.

### Wallenius and defect fixtures (F-O5, F-O6) — 2026-10-27

From the reference-fixtures chip ([brief](reference-fixtures-brief.md)), `a6960f0`. Three oracles, run by `r/scripts/bass_stage1_checks.R`, all checks passing. They source the pinned R BASS checkout (`2671014`) through `_common.R`'s `bass_source`, and evaluate R's expressions as written, extracted from the function bodies, rather than transcriptions of them.
- **F-O5,** `r/oracles/bass_wallenius_oracle.R`. The reference is brute force over the K! draw orders. It sums to 1 over every enumerated order and equals 1/C(p, K) under equal weights, to 3·10⁻¹⁶. Its denominators sum the weights left in the urn: subtracting the drawn ones from the total, as a first version did, cost 9·10⁻¹³ with weights spanning five decades. R's own `sample(prob = z, replace = FALSE)` draws that law (largest |z| 2.8 over 20,000 draws per order).
  - *R's range.* `dmwnchBass` is exact at K = 2, 3: 5·10⁻¹⁵ with equal weights (p = 6), 1.4·10⁻¹³ with R-like counts (p = 8), 1.6·10⁻¹² with p = 50.
  - *Beyond it* the alternating sum loses digits: 5·10⁻¹³ at K = 4, 5 (p = 8) and 4·10⁻⁹ at K = 5 (p = 50). With weights spanning five decades it loses them from K = 2: 1·10⁻⁷ at K = 2, 3·10⁻⁵ at K = 3, 3·10⁻⁴ at K = 4, 8·10⁻³ at K = 5. Quadrature of the positive integrand (on t = e^{−u}) stays within 2·10⁻¹⁵ of brute force everywhere.
  - *At K = 1* `dmwnchBass` returns twice the probability, because its loop runs over `1:0`. R never calls it there (`logProbChangeMod` scores the uniform 1/p it samples), so this is a trap for a literal port, not an R defect.
- **F-O6 (a),** `r/oracles/bass_enumeration_oracle.R`, mode `swap`. The structure space is p = 2 with N = 24 on a 4 × 3 grid, K_max = 2 and M_max = 2: 62 terms, 30 with at least 3 nonzero values, 464 states and 6930 moving transitions. The continuous parameters are fixed at σ² = 0.1, τ = 1/N, λ = 2, the slice on which R's structure moves act. Stationary laws come from each kernel's transition matrix by GTH elimination, which keeps every component's relative accuracy (a dense solve left the corrected swap's spread at 10⁻⁶).
  - *Validity.* R's `getQf` agrees with the strict rank in both column orders on every state. R's `alpha` is exactly b times its b = 1 value. At b = 1, R's kernel with the state's move probabilities in place of its fixed thirds is the posterior to 2·10⁻¹⁴. With the thirds it is 0.020 off in total variation: F-O4's defect, measured in passing.
  - *R's tempered kernel targets none of the candidate laws.* At b = 0.5 and 0.2 its total variation is 0.086 and 0.145 from D05's likelihood path (prior times p_b), 0.036 and 0.078 from the prior times p_1^b, 0.26 and 0.51 from the powered collapsed path, and 0.32 and 0.58 from the joint over (structure, a) powered and then integrated.
  - *Swaps.* R's swap fails balance against its own rung laws: the spread of log s₁ − log s₂ − (b₁ − b₂)·lp is 1.78 (1 ↔ 0.5) and 1.07 (0.5 ↔ 0.2), where balance needs 0. Against the four candidate families it is 0.82 to 5.26. D05's corrected kernels reach their target to 4·10⁻¹⁶, and D05's swap balances on them to 5·10⁻¹⁵.
- **F-O6 (b),** `r/oracles/bass_jeffreys_oracle.R`. N = 40, an intercept and two hinges, σ² = 0.04. A tensor Gauss–Legendre rule for the integral over a matches the Gaussian closed form to 2·10⁻¹⁴. That confirms D03's coherent factors per birth, c(2πσ²)^{1/2}(|G′|/|G|)^{−1/2} (flat) and (2π)^{1/2} (Jeffreys), and its σ² shapes, g₁ + (N − k)/2 and g₁ + N/2.
  - *R's factors,* evaluated as written: the births give 2.667 and 2.272 on the log scale, against −0.967 and −0.179 (flat) and 0.919 (Jeffreys). The change move gives −0.047 against 0.093 and 0. The σ² shape is 19.5 against 19 and 20.
  - *The implied priors are confirmed by quadrature.* σ^{−2k}|G|^{3/4} reproduces R's birth, death and change factors, and σ^{−(k−1)} reproduces its σ² shape, each to 1·10⁻¹⁴. Two different priors, so no single target, as D03 says.
- **F-O6 (c),** `bass_enumeration_oracle.R`, mode `count`. A constant likelihood, p = 2 with two knots per variable, K_max = 2, M_max = 3, λ = 2: 24 terms and 2325 states. D01's birth ratio gives the stated count law (0.165, 0.330, 0.315, 0.190 for M = 0–3) to 6·10⁻¹⁴ in total variation. KHAOS's double count gives the stated law times 1/M! (0.241, 0.482, 0.230, 0.046) to 6·10⁻¹⁴, 0.229 from the stated one.
- **Findings for the drafts** (the briefs are not edited):
  - *D01, the inclusion–exclusion cut.* The cut at K ≤ 3 holds at R's default w₂ = 5, where the weights stay within a small factor of each other. With a small w₂ or skewed counts, the sum loses digits from K = 2 (above). Choosing by an error estimate, about ε·2^K/P(S), rather than by K alone would cover both. A port that routes K = 1 through the Wallenius formula must not copy R's loop.
  - *D05, R's lp.* `lp` departs from the collapsed log posterior by more than the term priors. It is the log joint density at a's conditional mode, so it also carries ½ log|G| − (k/2) log(2πσ²/(1 + τ)), which depends on the structure and on R's column scaling (`const`). The spreads above include both.
- **Not covered:**
  - The enumeration has p = 2, so every interaction draws the one pair and the Wallenius factor is 1.
  - R's tempered draws of λ, σ² and τ are held fixed on the slice.
  - Categorical and functional terms, and the robust modes.
  - The Jeffreys fixture is one design at one σ².

### Centred admission fixture (F-O7) — 2026-10-26

`r/oracles/bass_admission_oracle.R` gains two modes:
- `centred` reads the same design file as the default mode, with the intercept first, and drops the intercept.
- `hadamard` reads A and H for direct fits and builds the materialized design itself, for fixtures only.

The default mode's values are unchanged, apart from a new `mode` field. The QR of the explicitly centred design is the reference; it uses tol = 0, because R's LINPACK default drops nearly dependent columns. Cross-products are summed in plain double, pairwise or left to right, so that the numbers do not depend on the platform's long double. All checks are run by `r/scripts/bass_stage1_checks.R` and pass:
- *The identity.* On a hinge design, uncentred s = centred s · ρ to 3·10⁻¹⁵. The default mode agrees with the QR angles to 1·10⁻¹².
- *Near-constant column* (N = 10⁵, an indicator true for all but one case, ρ = 10⁻⁵).
  - Centred s is 1.0·10⁻⁹, which the Liang gate admits; uncentred s is 1.0·10⁻¹⁴, which R's convention refuses.
  - Relative errors in s: 3.4·10⁻³ by the downdate and 3.2·10⁻⁷ explicitly centred, both summed pairwise (the check asks for at least 100× apart, inside the band). Explicitly centred but summed left to right, 6.9·10⁻⁴.
  - With a thousand zeros (ρ = 10⁻²) the downdate and explicit forms agree, at 4·10⁻⁸ and 2·10⁻⁷. Left to right gives 1.8·10⁻³.
- *Centred Hadamard form* (n = 60, m = 80, means large against spread, with one design-only and one functional-only term). It matches explicit centring to 2·10⁻¹⁶ in the Gram and 3·10⁻¹⁶ in s; the downdate gives 4·10⁻¹³ and 1·10⁻¹³.

This corrects the scratch reading in the report below. Most of that error came from the summation order on near-constant columns, not from accumulating over the rows in general. D06's numerics bullet is revised to match, and the numerics segment gains N10 (pairwise or compensated cross-products). No C# side yet.

### D06's Gram under the Liang convention — 2026-10-26

From the user's review. D06 now says the gate reads the Gram that D03's intercept convention maintains. Under the Liang default that is the centred Gram, with s_j the reciprocal of the usual VIF. Under R's convention it is the uncentred Gram, with the intercept as a column. The two are related by uncentred s_j = centred s_j · ρ_j. `BassLiteral` reads `getQf`'s pivots rather than s_j.

A scratch check (not a fixture) confirmed the identity to 3·10⁻¹⁵ and found two numerical points, now in D06's numerics bullet:
- Accumulating the Gram over N rows, not ε alone, sets the error of s near the gate: up to 2·10⁻² at N = 10⁵.
- D03's rank-one downdate loses about ε/ρ_j on near-constant columns: 4·10⁻² against 8·10⁻⁴ at ρ = 10⁻⁵.

Both stay inside the factor-of-two confirmation band. The strict rebuild under the Liang convention uses explicitly centred columns, or the three-term centred Hadamard form in direct fits. New fixture item: F-O7.

### Andros et al. read — 2026-10-26

From arXiv 2406.18751 and 2602.09351. The direct functional algebra records that covariates known at the output locations enter as further functional variables (R's `xx.func` already allows it), and D01's knot wording now covers several functional variables. A residual-dependence check for direct fits is proposed (F11), and the output-side inner product is declared as prep for correlated residuals (F-X6). The deterministic-simulator profile reads that check.

### Lineage insights folded in — 2026-10-26

D01 gains the prior-cancelling schedule as an option and two "considered, not adopted" notes (KHAOS's delayed rejection, a uniform prior over basis functions). The deterministic-simulator profile is drafted. The hub's register ("What the lineage taught") records every insight and where it landed.

### First reference fixtures — 2026-10-25

`33b99bc`: `bass_marginal_oracle.R` and `bass_admission_oracle.R`, with `r/scripts/bass_stage1_checks.R` reproducing the D03 propriety and D06 admission numbers. No C# side yet; parity tests come with F3 and F7.

### D01–D06 drafted — 2026-10-25

Proposals, not decisions: `097e25b` (D01–D03), `53a4876` (D04, D06), `0844bdc` (D05), `a32c701` (categorical terms), `77274fb` (Jeffreys mode, propriety, order-free gate, s_min).

