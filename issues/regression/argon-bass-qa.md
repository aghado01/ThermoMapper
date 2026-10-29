 Because untempered chains sit in distinct bands of $M$, Vehtari et al.’s split-chain bulk ESS is algebraically clamped near $C/2 \approx 4$ whenever between-chain variance $B$ dominates within-chain variance $W$—so `min-ESS/s` at a fixed 5,000 + 2,000 budget was measuring wall-clock speed against a clamped numerator of ~5.

Here is how I would lock in **Run Lengths + D10**, frame **Friedman’s worse modes**, and close **D09** to launch Stage 7.

---

### 1. Run Lengths & D10: What a Fit Reports When Chains Never Agree

Keeping `BurnIn = 9000, Samples = 1000` as the default budget and folding convergence accounting into **D10** is the right call. Right now, [`BassDiagnostics`](file:///D:/aghado01/ThermoMapper/src/maths/regression/spline/bass/BassDiagnostics.cs#L26-L40) only records `StoppedEarly` (`SamplesUsed < Samples`), which is `false` both when `RHatTarget == 0` (the default) and when `RHatTarget > 1` runs out of budget without converging—and in [`ChainEnsemble.Run`](file:///D:/aghado01/ThermoMapper/src/maths/samplers/ensemble/ChainEnsemble.cs#L136), `MayStop` isn't even evaluated on the final batch (`collected == total`).

For **D10**, I recommend three concrete rules when chains do not agree by the end of `Samples`:

1. **Always return the fit (never throw or null it), with an explicit convergence verdict on [`BassDiagnostics`](file:///D:/aghado01/ThermoMapper/src/maths/regression/spline/bass/BassDiagnostics.cs):**
   - Checkpoints (`fit.Checkpoint`) and `Bass.Resume` exist so a caller can inspect or extend an unconverged run.
   - Add explicit summary properties to `BassDiagnostics` (or a small `BassConvergenceSummary` record):
     - `MaxRHat` and `MinEss` across the monitored quantities;
     - `Converged` (whether every monitored functional has `RHat.MeetsTarget(target)` and `Ess.MeetsMinimum(C * MinEssPerChain)` at the end of the run, evaluated against `RHatTarget` when $> 1$, or a standard reference threshold of $1.05$ when `RHatTarget == 0`).
2. **Separate *predictive* agreement from *structural* agreement:**
   - Look at Phase 2’s $\hat{R}$ table ([local-move-pass-brief.md](file:///D:/aghado01/ThermoMapper/issues/regression/bass/local-move-pass-brief.md#L224-L258)): on `mixed` and `functional`, $\hat{R}_M$ is $2.15\text{–}2.50$, **yet $\hat{R}_{\sigma^2}$ is $1.00\text{–}1.01$, $95\%$ surface coverage is $98\text{–}99\%$, and $\sigma$ is covered in $5/5$ seeds.**
   - Why? Because an overcomplete hinge basis is not structurally identifiable: one chain can represent the same surface $f(x)$ with $M = 5$ terms while another uses $M = 7$ co-adapted terms.
   - If D10 reports both **Predictive $\hat{R}$** ($\max \hat{R}$ over $\sigma^2$ and `Fitted` at the monitored points) and **Structural $\hat{R}$** ($\max \hat{R}$ over $M$ and $\tau$), a user can immediately distinguish:
     - *Representation non-identifiability* (`PredictiveRHat` $\le 1.05$, `StructuralRHat` $> 1.1$): chains agree on the surface $f(x)$ and noise $\sigma^2$, and only differ in how they factored $f$ into hinges.
     - *Genuine mode disagreement* (`PredictiveRHat` $> 1.05$, as on untempered Friedman where `Fitted` $\hat{R}$ is $1.32\text{–}1.73$ and $\sigma^2$ $\hat{R}$ is $1.08$): chains are predicting different surfaces.
3. **Allow choosing which scope `RHatTarget` gates:**
   - A small enum or option on `BassSampling` (e.g., `StoppingScope.AllMonitored` by default vs `StoppingScope.PredictiveOnly` for $\sigma^2$ + `Fitted`) lets callers who care about surface predictions stop when $f$ and $\sigma^2$ agree without waiting for $M$ to mix across equivalent representations.

---

### 2. Friedman’s Worse Modes: The Right Metric (and Stage 8’s Role)

Park further MCMC move work on Friedman’s worse modes for now, and record two points for when you revisit it:
1. **The right mode-crossing metrics:**
   - **Between-chain log-posterior spread** $\operatorname{sd}_c(\overline{\ln\pi}_c)$ (or the gap $\max_c \overline{\ln\pi}_c - \min_c \overline{\ln\pi}_c$, which was 12–42 nats in M4) and **between-chain surface $L_2$ distance** $\max_{c, c'} \|\bar{f}_c - \bar{f}_{c'}\|_2$. Unlike split-$\hat{R}$ / bulk ESS (which saturates at $\text{ESS} \approx C/2$ as soon as chains sit in separate modes), both measure how far apart the modes are in function and density space.
   - **Hitting / agreement time** from dispersed prior starts: sweeps or seconds until all chains enter within $\Delta$ nats of the leading chain (or until `PredictiveRHat` drops below $1.10$).
2. **Stage 8’s Fitting SMC (`C4`) will test this directly without new local moves:**
   - In MCMC without round trips, a chain that settles into a worse co-adapted representation of $x_1 x_2$ during burn-in is trapped there forever unless it simultaneously tears down 3–5 load-bearing terms.
   - In tempered SMC along D05’s path $b \in [0, 1]$, particles that fall 12–42 nats behind the leading representation are **eliminated by resampling weights** as $b$ increases—so Stage 8 will tell you whether population resampling solves Friedman’s mode trapping on the existing local kernel.

---

### 3. Ready-to-Decide Summary for D09 (Unblocking Stage 7)

To kick off Stage 7 right now, here is the complete D09 decision set ([port-calibration-brief.md](file:///D:/aghado01/ThermoMapper/issues/regression/bass/port-calibration-brief.md#L38-L83)):
- **Target:** Cut posterior $p_{\text{cut}}(e, \theta \mid D_{\text{sim}}, y) = p(e \mid D_{\text{sim}})\,p(\theta \mid y, e)$, with emulator draws behind an interface (**C-X2** built in Stage 7).
- **Formulations:** Declared `Functional` (routed through [`DiagonalPlusLowRank`](file:///D:/aghado01/ThermoMapper/src/maths/linalg/DiagonalPlusLowRank.cs)) and `Multivariate` (routed through [`DenseCholesky`](file:///D:/aghado01/ThermoMapper/src/maths/linalg/DenseCholesky.cs)), all in physical units (`C-O2`).
- **Variance-group priors:** Declared per group; default $\text{IG}(\nu_g/2, (\nu_g/2)\hat{s}_g^2)$, with the heavy-tailed $(1 + \sigma_g^2 / s_{0,g}^2)^{-1}$ kernel selectable per group (marked `IsProper = false` so $b = 0$ on the likelihood path is refused).
- **Discrepancy prior (`Functional`):** $\delta = \mu_\delta + D v$ with $v \sim \mathcal{N}(0, \sigma_v^2 I_q)$, declared with default $\sigma_v = 1.0$ ($\tau_v = 1$, so $D$ defaults to the physical square-root factor $\Sigma_\delta = D D^\top$, while `PriorSd` $\sigma_v$ scales a dimensionless basis $D$).
- **Categorical calibration parameters:** Exact joint discrete conditional by enumeration up to a combination cap (e.g. $1{,}024$), single-site Metropolis/Gibbs above it.
- **Stage 7 reference sampler:** Nested conditional MCMC ($S$ outer draws $e_s$, inner chain on $p(\theta \mid y, e_s)$ to stationarity, inner-length error measured).