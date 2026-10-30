Analyzed bass-questions.md:L1-122
Analyzed stage8-fitting-smc-brief.md:L1-150
Analyzed stage8-cut-smc-brief.md:L1-138
Analyzed Readouts.cs:L1-283
Analyzing SmcSampler.cs
Analyzing SmcStage.cs
Analyzing SmcSampler.cs:L1-300
Ran list SMC and calibration files
Analyzed SmcContracts.cs:L1-53
Analyzed SmcRun.cs:L1-463
Analyzed SmcWeights.cs:L1-120
Analyzing SamplingTarget.cs
Analyzed CalibrationModel.cs:L1-357
Analyzed CalibrationSampler.cs:L1-406
Analyzed CalibrationBlocks.cs:L1-348
Analyzed NestedCalibration.cs:L1-262
Ran list files in src/maths/samplers
Ran list all files in src/maths/samplers
Analyzed SamplingTarget.cs:L1-180
Analyzed ConvergenceVerdict.cs:L1-126
Canceled check for Akinshin paper locally
Killed task Check for Akinshin paper locally
Ran list gauntlet paper directory
Used tool: manage_task

## 1. Readouts on N14 (`2d1ba7d`) and Stage 7 (`2232dc3`)

Both landings are clean:
- **N14 ([`ConvergenceVerdict.cs`](file:///D:/aghado01/ThermoMapper/src/maths/samplers/ensemble/ConvergenceVerdict.cs#L81-L125)):** Having `ConvergenceVerdict.Of` drive both the inter-round stop in [`ChainEnsemble`](file:///D:/aghado01/ThermoMapper/src/maths/samplers/ensemble/ChainEnsemble.cs) and the end-of-run verdict eliminates any drift between stopping and reporting, and failing on frozen readings matches D10. The S10 log-posterior gaps (`21.8` to `42.8` nats across chains, and `27.0` nats on tempered Friedman) immediately validate monitoring `LogPosterior`: the chains really are sitting in different basins.
- **Stage 7 (`calibration/`):** The implementation in [`CalibrationModel.cs`](file:///D:/aghado01/ThermoMapper/src/maths/regression/spline/bass/calibration/CalibrationModel.cs), [`CalibrationBlocks.cs`](file:///D:/aghado01/ThermoMapper/src/maths/regression/spline/bass/calibration/CalibrationBlocks.cs), [`CalibrationSampler.cs`](file:///D:/aghado01/ThermoMapper/src/maths/regression/spline/bass/calibration/CalibrationSampler.cs), and [`NestedCalibration.cs`](file:///D:/aghado01/ThermoMapper/src/maths/regression/spline/bass/calibration/NestedCalibration.cs) honours every D09 invariant (untempered $\ln \sigma^2$ Jacobian in [`VarianceBlock.Update`](file:///D:/aghado01/ThermoMapper/src/maths/regression/spline/bass/calibration/CalibrationBlocks.cs#L192-L203), likelihood-path precision $I/s_v^2 + b D^\top \Sigma^{-1} D$ in [`DiscrepancyBlock.Conditional`](file:///D:/aghado01/ThermoMapper/src/maths/regression/spline/bass/calibration/CalibrationBlocks.cs#L119-L133), and per-target [`PredictionCache`](file:///D:/aghado01/ThermoMapper/src/maths/regression/spline/bass/calibration/CalibrationModel.cs#L239) / [`FactorCache`](file:///D:/aghado01/ThermoMapper/src/maths/regression/spline/bass/calibration/CalibrationModel.cs#L240) ownership).
- **On Stage 7's three "Open, for later" items:**
  1. *Inner-run convergence rate under the default `1.01 / 100` gate:* Worth measuring once alongside C5's comparison pass, since [`NestedCalibration.NotConverged`](file:///D:/aghado01/ThermoMapper/src/maths/regression/spline/bass/calibration/NestedCalibration.cs#L68) excludes non-converged inner runs by default (`Pooled(includeNotConverged = false)`).
  2. *N15 (`DiagonalPlusLowRank` Woodbury/capacitance Cholesky vs. thin SVD):* Good call parking it as N15 until profile sizes warrant it.
  3. *Above-cap categorical single-site Gibbs:* Parked cleanly; most calibration problems have $\le 64$ level combinations and hit [`CategoricalBlock.Enumerate`](file:///D:/aghado01/ThermoMapper/src/maths/regression/spline/bass/calibration/CalibrationBlocks.cs#L295-L322).

---

## 2. Stage 8 Plan & the Four Revisions

**Endorse all four revisions (`1`, `3`, `6`, `10`) plus the ancestry-error note (`4`).** Together with the rest of the plan (`2`, `5`, `7`, `8`, `9`, `11`), the design is sound. Before launching the two chips in parallel, there are **two concrete technical catches** in the briefs worth tightening:

### Catch A — Item 5 ([`DrawSample.Quantile`](file:///D:/aghado01/ThermoMapper/src/maths/samplers/ensemble/Readouts.cs#L147-L199)): Akinshin's §6 formula vs. endpoint-anchored midpoint Type 7
In [stage8-fitting-smc-brief.md](file:///D:/aghado01/ThermoMapper/issues/regression/bass/stage8-fitting-smc-brief.md#L48-L51), Step 1 asks for:
1. Exact Type 7 at equal weights at any scale;
2. Continuity as weights move from equal;
3. $Q(0) = \min$ and $Q(1) = \max$ under unequal weights;
4. Tied values pooled and insertion-order independence.

- **Watch out:** Akinshin's (2023, arXiv:2304.07265, §6) $n_{\text{eff}}$-window Type-7 construction integrates the weighted empirical quantile step function $F_W^{-1}(u)$ over the window $\bigl[p(1 - 1/n_{\text{eff}}),\, p(1 - 1/n_{\text{eff}}) + 1/n_{\text{eff}}\bigr]$ of width $1/n_{\text{eff}} = \sum_i W_i^2$. Whenever an endpoint weight is smaller than $1/n_{\text{eff}}$, that window spills past the endpoint step—for example, on weights $(0.45, 0.45, 0.10)$, $1/n_{\text{eff}} = 0.415 > 0.10$, so at $p = 1$ the window $[0.585, 1.0]$ averages $x_{(2)}$ and $x_{(3)}$ instead of returning $Q(1) = x_{(3)} = \max$!
- **The exact formula that satisfies all four requirements:** Look at [`DrawSample.Quantile`](file:///D:/aghado01/ThermoMapper/src/maths/samplers/ensemble/Readouts.cs#L173-L199) as already written. It pools positive-weight tied values into `distinct` pairs $(x_i, W_i)$ with cumulative total $W = \text{total}$ and step midpoints $m_i = \sum_{j < i} W_j + \frac{1}{2}W_i$. The **only** reason it currently jumps at equal weights is [line 185](file:///D:/aghado01/ThermoMapper/src/maths/samplers/ensemble/Readouts.cs#L185):
  ```csharp
  double target = p * total; // maps [0, 1] -> [0, W] (Type 4/5 at equal weights)
  ```
  If line 185 instead maps $p \in [0, 1]$ affinely onto the midpoint span $[m_0, m_{k-1}] = [\tfrac{1}{2}W_0,\, W - \tfrac{1}{2}W_{k-1}]$:
  $$\text{target}(p) = \tfrac{1}{2}W_0 + p\cdot\bigl(W - \tfrac{1}{2}(W_0 + W_{k-1})\bigr)$$
  then:
  - At equal weights $W_i = c$ on $n$ distinct values, $\text{target}(p) = c\bigl(\tfrac{1}{2} + p(n - 1)\bigr)$ and $m_i = c(i + \tfrac{1}{2})$, so linear interpolation between $m_{i-1}$ and $m_i$ has fraction $t = p(n - 1) - (i - 1)$, **identically matching `Type7`**;
  - At $p = 0$, $\text{target}(0) = m_0 \implies Q(0) = x_0 = \min$; at $p = 1$, $\text{target}(1) = m_{k-1} \implies Q(1) = x_{k-1} = \max$;
  - It is a smooth rational function of the weights on the simplex, pools tied values, and keeps the equal-weights fast path bit-identical for BARS.

---

### Catch B — Items 1 & 11 (C5 Bridges on [`SmcRun<TState>`](file:///D:/aghado01/ThermoMapper/src/maths/samplers/smc/SmcRun.cs#L111) without editing `smc/`)
Keeping `src/maths/samplers/smc/` untouched (Revision 1) works cleanly, with two details C5 should know up front:
1. **First step of each bridge (`!ReferenceEquals(next.Target, Stage.Target)`):**
   - At the end of leg 1 (or bridge $s - 1$), `run.Stage` sits at `InverseTemperature == 1.0` on the previous target. [`SmcSampler.NextInverseTemperature`](file:///D:/aghado01/ThermoMapper/src/maths/samplers/smc/SmcRun.cs#L445-L462) only bisects `InverseTemperature` on `run.Stage.Target` from $b < 1$ to $1$.
   - Moreover, in [`SmcRun.EvaluationsFor`](file:///D:/aghado01/ThermoMapper/src/maths/samplers/smc/SmcRun.cs#L378-L385), calling `run.ConditionalEss(next)` when `next.Target` is a new bridge target (`!ReferenceEquals(next.Target, Stage.Target)`) re-evaluates all $N$ particles on every call. A 60-step bisection calling `run.ConditionalEss` directly for the first step $t_1 > 0$ of a bridge would evaluate the new target $60N$ times.
   - **Fix inside C5 (no `smc/` edits needed):** At the start of bridge $s$ ($e_{s-1} \to e_s$), evaluate $e_s$ on the $N$ particles once, bisect the first step $t_1 \in (0, 1]$ directly via [`SmcWeights.ConditionalEss`](file:///D:/aghado01/ThermoMapper/src/maths/samplers/smc/SmcWeights.cs#L53-L69) on the cached increments $t \cdot (\ell_s - \ell_{s-1})$, call `run.Advance(new SmcStage(bridgeTarget, new Temperature(t1, Likelihood)))`, and then call [`SmcSampler.TemperToPosterior(run, targetConditionalEss)`](file:///D:/aghado01/ThermoMapper/src/maths/samplers/smc/SmcRun.cs#L430-L438) for the remaining steps $t_1 < t_2 < \dots \le 1$ (where `Stage.Target` is now `bridgeTarget`, so `EvaluationsFor` reuses `_evaluations` at $O(1)$ cost).
2. **`SmcOptions.MaxStages` across a single continuous tour:**
   - Carrying a single [`SmcRun<CalibrationState>`](file:///D:/aghado01/ThermoMapper/src/maths/samplers/smc/SmcRun.cs#L111) across the first leg and all $S - 1$ bridges is what lets `Families` track surviving initial lineages across the whole tour. Because [`SmcRun`](file:///D:/aghado01/ThermoMapper/src/maths/samplers/smc/SmcRun.cs#L138) pre-derives `MaxStages + 1` stage seeds at `Start` (default `MaxStages = 1000`), C5 should scale `SmcOptions.MaxStages` with $S$ (e.g., `MaxStagesPerBridge * S`).

---

### Quick Verdict on the Remaining Revisions
- **Revision 1 (No edits to `smc/`):** Endorse. With the bridge pattern above, neither C4 nor C5 needs any changes in `src/maths/samplers/smc/`.
- **Revision 3 (`SmcSpec.Replicates = 4`, refuse `Chains` when set):** Endorse. Cleanly separates independent SMC populations from Markov chains.
- **Revision 4 (Report `Families` alongside weighted density / ESS readouts):** Endorse. After resampling, equal weights $1/N$ hide lineage collapse; reporting `Families` and between-replicate spread prevents over-interpreting within-population Kish ESS.
- **Revision 6 (Between-replicate spread of means in posterior-SD units instead of split-$\hat{R}$):** Endorse. Exchangeable particles have no time index for split-$\hat{R}$, and $\text{SD}_{\text{rep}}(\hat{\mu}_r) / \hat{\sigma}_{\text{post}}$ directly measures between-run agreement on the posterior's own scale.
- **Revision 8 (`Bass.Resume` for SMC):** Agree that it is the most deferrable item if wrapping [`SmcCheckpoint<BassState>`](file:///D:/aghado01/ThermoMapper/src/maths/samplers/smc/SmcRun.cs#L94-L97) alongside [`EnsembleCheckpoint`](file:///D:/aghado01/ThermoMapper/src/maths/samplers/ensemble/Checkpoint.cs) complicates `BassFit.Checkpoint`, since an adaptive SMC run finishes at $b = 1$ in one pass anyway.
- **Revision 10 (Tour distance on $\operatorname{SD}_i(\ell(\theta_i; e_s) - \ell(\theta_i; e_{s'}))$ over the first leg's particles):** Strongly endorse. Conditional ESS on a geometric bridge is a decreasing function of $\operatorname{Var}_W(\ell_s - \ell_{s'})$ under the current posterior population, whereas $L_2$ distance on a uniform box is dominated by high-error corners where $p(\theta \mid y, e_1)$ has zero mass. Note that the first leg's particles only need to be evaluated once per emulator draw ($S \times N$ evaluations), and those evaluations can be reused if cached or summarized into an $S \times N$ log-likelihood matrix before running nearest-neighbour + 2-opt.