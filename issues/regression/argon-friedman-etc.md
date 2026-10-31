Analyzed bass-questions.md:L1-122
Analyzed stage8-fitting-smc-brief.md:L140-158
Ran check recent git commits
Ran check files in C4 and C5 step 4 commits
Analyzed stage8-cut-smc-brief.md:L130-149
Analyzed port-fitting-brief.md:L160-320
Analyzed port-fitting-brief.md:L321-480
Analyzed port-fitting-brief.md:L481-590
Analyzed stage8-fitting-smc-brief.md:L145-158
Analyzed stage8-cut-smc-brief.md:L135-149
Analyzed BassSmcQualificationTests.cs:L1-209
Ran check git status and recent commits
Analyzed local-move-measurements-brief.md:L1-160
Analyzed stage4-tempering-brief.md:L1-160
Analyzed CollapsedGPrior.cs:L1-180

Here is how each of those three open items stands and what to expect from them:

---

### 1. Friedman’s worse modes & C4’s SMC qualification
- **Why MCMC chains stall 24–42 nats below the best mode on Friedman:**
  - In [`BassMixingModeTests`](file:///D:/aghado01/ThermoMapper/issues/regression/bass/local-move-measurements-brief.md#L95-L105) (M4), chains that early on approximate $10\sin(\pi x_1 x_2)$ or $20(x_3 - 0.5)^2$ with suboptimal hinge placements cannot transition to the cleaner structure via single-term moves: deleting a load-bearing term drops $\ell_1$ by tens of nats, while adding the replacement first either trips [`AngleFloor`](file:///D:/aghado01/ThermoMapper/src/maths/regression/spline/bass/BassConfig.cs#L436) or pays the Occam penalty.
  - Under parallel tempering (M5), the structural transition concentrates in a steep band around $b \approx 0.16\text{–}0.32$ where pair rejections stay near `1.0` even at $\Delta b = 0.01$, choking round-trip swaps.
- **Why SMC (`C4`, [`BassSmcQualificationTests.cs`](file:///D:/aghado01/ThermoMapper/tests/regression/BassSmcQualificationTests.cs)) has a structural advantage here—and what to check in its readout:**
  - Unlike parallel tempering, SMC walks **one-way** from $b = 0 \to 1$ with $N = 1{,}000$ particles in lockstep (`TargetConditionalEss = 0.99`, `MutationSteps = 5`). It does not need round-trip swaps across $b \in [0.16, 0.32]$: if even a small fraction of particles assemble the better $(x_1, x_2)$ and $x_3^2$ hinges while crossing that band, incremental reweighting $\exp((b' - b)\,\Delta \ell)$ amplifies them across the population.
  - **The two numbers in C4's report that will decide the next step:**
    1. **Per-replicate log-posterior means/medians vs. `-892.4`:** Do all 4 replicates consistently find the `-892` basin at $N = 1{,}000$ (or $4{,}000$)?
    2. **Surviving lineages (`Families`):** If only 1–2 particles out of $1{,}000$ find the `-892` basin before $b = 0.32$, the $20+$ nat advantage will collapse `Families` to single digits during resampling.
  - **Decision rule after C4 reports:**
    - If SMC reliably lands all replicates in the top mode with healthy `Families` (or small between-replicate spread on predictions and $\sigma$), a pre-registered pass weighing SMC as a candidate default (or recommended profile) for structured fits is well earned.
    - If even $N = 4{,}000$ shows lineage collapse (`Families` $\ll 20$) or split replicates on Friedman, then the bottleneck is local mutation through the transition band—pointing to a **compound replacement move** (swapping one load-bearing term for a newly proposed term on the same variables in a single RJMCMC step) rather than a sampler-shell swap alone.

---

### 2. F-R1, the deterministic-simulator profile & the ridge prior
- **Why F-R1 is non-blocking right now:**
  - Under D03’s proper defaults ([`Sigma2Prior.ScaledToResponse(1e-3)`](file:///D:/aghado01/ThermoMapper/src/maths/regression/spline/bass/BassConfig.cs#L423) and [`GPriorCentre.Mean`](file:///D:/aghado01/ThermoMapper/src/maths/regression/gprior/CollapsedGPrior.cs#L21)), exact or near-exact interpolation on noise-free simulators is already proper ($2g_2 + S(\tau) > 0$). So deterministic simulators already work under `SampledG`.
- **What building `SampledRidge` and `FixedRidge` requires:**
  - In [`CollapsedGPrior`](file:///D:/aghado01/ThermoMapper/src/maths/regression/gprior/CollapsedGPrior.cs#L66), the prior precision $\tau G_c$ is proportional to $G_c = Z_c^\top Z_c$, so factoring $G_c$ once makes every $(\tau, b)$ evaluation $O(1)$ via [`CentredMarginalAt`](file:///D:/aghado01/ThermoMapper/src/maths/regression/gprior/CollapsedGPrior.cs#L159-L167).
  - Under a ridge prior $\beta_c \mid \sigma^2 \sim \mathcal{N}(0, \sigma^2 (\tau I_M)^{-1})$, the posterior precision at inverse temperature $b$ is $b G_c + \tau I_M$, and the marginal needs $\ln |I_M + \frac{b}{\tau} G_c|$ and $\tilde{y}^\top Z_c (b G_c + \tau I_M)^{-1} Z_c^\top \tilde{y}$.
  - **The key numerics design when you build it (`N-X1`):**
    - During **birth/death/change proposals** at fixed $(\tau, b = 1)$, [`BorderedCholesky`](file:///D:/aghado01/ThermoMapper/src/maths/linalg/BorderedCholesky.cs) can maintain $G_c + \tau I_M$ in $O(M^2)$ per proposal, just as it maintains $G_c$ today.
    - For **$\tau$ slice-sampling and tempered swaps ($b \ne 1$)** on the current structure, refactoring $b G_c + \tau I_M$ from scratch at every slice evaluation would cost $O(M^3)$. Instead, caching the eigendecomposition $G_c = V \Lambda V^\top$ (and projected right-hand side $w = V^\top Z_c^\top \tilde{y}$) once when a structure is accepted reduces every $(\tau, b)$ evaluation to **$O(M)$ scalar arithmetic**:
      $$\ln\left|I_M + \frac{b}{\tau}G_c\right| = \sum_{j=1}^M \ln\left(1 + \frac{b\lambda_j}{\tau}\right), \qquad Q_{\text{ridge}}(\tau, b) = \sum_{j=1}^M \frac{b\,\lambda_j w_j^2}{b\lambda_j + \tau}$$
  - **Recommendation:** Keep F-R1 parked until Stage 8 closes. When you schedule it, build `CollapsedRidgePrior` with that $O(M)$ spectral evaluation for $(\tau, b)$ and run the pre-registered comparison (`SampledG` vs. `SampledRidge` vs. `FixedRidge`, plus wider count prior $h_1, h_2 < 10$) on noise-free simulator benchmarks.

---

### 3. OQ2, exact reference rungs on long ladders
- **Why the 4- and 8-rung study showed no cold-chain ESS payoff, and why long ladders rarely change that:**
  - The exact reference rung ([`TemperingSpec.ReferenceRung`](file:///D:/aghado01/ThermoMapper/src/maths/regression/spline/bass/BassConfig.cs#L371), powered by [`BassPriorSampler`](file:///D:/aghado01/ThermoMapper/src/maths/regression/spline/bass/BassPriorSampler.cs) at $b = 0$) replaces the hottest MCMC replica with i.i.d. draws from the prior on $\mathcal{F}$.
  - On **mild surfaces**, the hot MCMC rung at $b_{\min} = 1/N$ already mixes near the prior in a few birth/death steps, so zeroing out its autocorrelation at $b = 0$ barely moves the cold chain's ESS.
  - On **hard multimodal surfaces** (like Friedman), the bottleneck is **not** slow mixing at the hot end ($b \to 0$)—it is the transport barrier across the intermediate phase-transition band ($b \approx 0.16\text{–}0.32$, M5). An exact prior draw injected at $b = 0$ starts with random knots and $M \sim \text{NegBin}(10, 10)$; as it swaps up toward $b = 0.16$, it still has to be sculpted by local moves through that same bottleneck.
- **What `ReferenceRung` *does* buy on any ladder length:**
  1. **Normalized evidence:** It anchors $b = 0$ at exact prior normalization ($\ln Z_0 = 0$) so [`BassEvidence`](file:///D:/aghado01/ThermoMapper/src/maths/regression/spline/bass/BassEvidence.cs) reports the full model evidence $\ln p(y)$.
  2. **Guaranteed prior regeneration (`AnnealedRestarts`):** Every round trip to $b = 0$ is a true renewal event, independent of how many terms $M$ the hot chain was carrying.
- **Recommendation:** Treat OQ2 as **settled in practice**—keep `ReferenceRung` for its evidence anchor and exact renewal diagnostic, rather than spending a separate long-ladder study hunting for an ESS speedup that is bottlenecked at intermediate $b$, not at $b = 0$.

---

# Opus chat's read 

I agree with the recommendations on F-R1 and OQ2. On Friedman the diagnosis is close, but I read the M5 numbers differently, and that changes what C4 should report. F-R1 also has a formula to fix before it's coded.

### Friedman

**The band is probably a phase transition, not the barrier between the two fitted modes.**
- To first order in Δb, a swap's log acceptance is −Δb times the difference in ∂_b ln p_b between the two states.
- So rejection near 0.95 at Δb = 0.01 means adjacent replicas differ by about 300 nats in that slope.
- The stalled modes are only 24–42 nats apart at b = 1, so the gap between neighbours is much larger than the gap between the modes.
- A gap that size looks like coexistence of two phases: a prior-dominated one (few terms, σ² near var(y)) and a fitted one. That's a first-order transition, the Potts q > 4 situation.
- To check it, histogram ℓ at fixed b inside the band. Two peaks means coexistence.
- If that's what it is, the cold chain's structure is settled when it condenses through the band, and each PT round trip is one more draw from that lottery.

Syed et al.'s communication barrier makes this quantitative. Non-reversible PT's round-trip rate is bounded by 1/(2 + 2Λ), with Λ = ∫λ(b) db. If the band dominates Λ, no ladder length gets around it. Your ladder adaptation's rejection estimates already give λ(b).

**That changes the account of SMC's advantage, and what C4 should read.** SMC wins here because it runs N independent condensations in one pass and weights them against each other. Avoiding round trips is secondary. Three consequences:
- *Families will usually be low in a correct run.* After a sharp transition, every surviving particle descends from the few that condensed into the best basin. Low Families tells you how many independent condensations the answer rests on, not whether it's right. I wouldn't make Families ≪ 20 a failure criterion on its own.
- *Decide on replicate agreement.*
  - SMC gives ln Ẑ for free. A spread above about a nat across the four replicates means at least one is missing mass.
  - Also report weighted basin fractions at b = 1 under a structural key, such as the multiset of variable sets in the terms.
  - Comparing SMC's ln Ẑ with `BassEvidence`'s ladder estimate on the same data measures how much mass PT misses.
- *Compare masses, not peak heights.* Each replicate's mean log posterior against −892.4 compares density heights. If knots are continuous and the modes differ in M, heights don't compare masses. An extra knot whose posterior width is about 0.01 inflates its basin's peak by about 4.6 nats relative to its mass. Basin fractions and ln Ẑ don't have this problem.

**Before building the compound replacement move, test `AngleFloor`.** If adding the replacement trips the floor, the intermediate state (old term plus replacement) is outside the support. Part of the barrier would then come from the floor rather than the posterior. Rerun M4 with a lower floor. If the gap shrinks, a swap that moves directly between in-support endpoints, checking the floor only at the endpoint, removes that part. Its reverse probability has to be exact, under the same conditions as E4's kernels.

### F-R1

Agree it should stay parked. Two things for when it's scheduled:

- **The spectral formula.**
  - With w = VᵀZ_cᵀỹ as argon defines it, ỹᵀZ_c(bG_c + τI)⁻¹Z_cᵀỹ = Σ w_j²/(bλ_j + τ), and the reduction in the exponent is b² times that.
  - Argon's Q = Σ bλ_j w_j²/(bλ_j + τ) is right only under two conditions. First, w_j must be the projection onto Z_c's left singular vectors, i.e. the j-th entry of VᵀZ_cᵀỹ divided by √λ_j. Second, the result needs one more factor of b outside.
  - Sanity check: as τ → 0 the correct form gives b·ỹᵀZ_cG_c⁻¹Z_cᵀỹ, the OLS explained sum of squares.
  - The unwhitened form is also safer when G_c is near singular, since it never divides by √λ_j.
- **Column scale.** A common τ on raw hinge columns depends on how the columns are scaled, and hinge columns vary a lot in scale.
  - For uniform inputs, a hinge (x − t)₊ has squared norm about n(1 − t)³/3, and products of hinges compound the variation.
  - So the ridge shrinks terms by where their knots sit, which the g-prior doesn't do.
  - KHAOS found a fixed ridge best on noise-free simulators, but with an orthonormal polynomial basis, where this doesn't come up.
  - Unless the port already normalizes columns, pre-register the normalization as part of the comparison. Otherwise the comparison measures the normalization as much as the prior.

### OQ2

Agree it's settled. The reference rung can't help with a barrier in the middle of the ladder, and its evidence anchor is reason enough to keep it. Record the Λ estimate alongside the decision, since it also says whether more rungs could ever help.