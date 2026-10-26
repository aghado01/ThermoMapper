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
- [ ] **F-O5** Wallenius probabilities against brute-force sums over draw orders.
- [ ] **F-O6** Defect fixtures: R's tempering and swap density, R's Jeffreys mode, KHAOS's double count factor.

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

### Andros et al. read — 2026-10-26

From arXiv 2406.18751 and 2602.09351. The direct functional algebra records that covariates known at the output locations enter as further functional variables (R's `xx.func` already allows it), and D01's knot wording now covers several functional variables. A residual-dependence check for direct fits is proposed (F11), and the output-side inner product is declared as prep for correlated residuals (F-X6). The deterministic-simulator profile reads that check.

### Lineage insights folded in — 2026-10-26

D01 gains the prior-cancelling schedule as an option and two "considered, not adopted" notes (KHAOS's delayed rejection, a uniform prior over basis functions). The deterministic-simulator profile is drafted. The hub's register ("What the lineage taught") records every insight and where it landed.

### First reference fixtures — 2026-10-25

`33b99bc`: `bass_marginal_oracle.R` and `bass_admission_oracle.R`, with `r/scripts/bass_stage1_checks.R` reproducing the D03 propriety and D06 admission numbers. No C# side yet; parity tests come with F3 and F7.

### D01–D06 drafted — 2026-10-25

Proposals, not decisions: `097e25b` (D01–D03), `53a4876` (D04, D06), `0844bdc` (D05), `a32c701` (categorical terms), `77274fb` (Jeffreys mode, propriety, order-free gate, s_min).

