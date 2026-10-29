Analyzed bass-questions.md:L1-65
Analyzed port-fitting-checklist.md:L1-250
Analyzed BassQualificationRecoveryTests.cs:L1-171
Analyzed BassQualificationLooTests.cs:L1-194

All five recommendations in S10’s report ([port-fitting-checklist.md](file:///D:/aghado01/ThermoMapper/issues/regression/bass/port-fitting-checklist.md#L92-L97)) are sound—here is my take on each, plus two technical findings from reading [`BassQualificationRecoveryTests.cs`](file:///D:/aghado01/ThermoMapper/tests/regression/BassQualificationRecoveryTests.cs) and [`BassQualificationLooTests.cs`](file:///D:/aghado01/ThermoMapper/tests/regression/BassQualificationLooTests.cs).

---

### 1. The Five S10 Decisions

1. **Friedman’s failed $\sigma$ criterion under untempered defaults:** **Keep as a recorded failure of the untempered defaults.**
   - Untempered BASS gives $\sigma \in [1.013, 1.161]$ (missing $1.0$ by $0.013$) with $\hat{R}_M = 2.35$, because at $M \approx 9.6$ piecewise-linear hinges on Friedman’s smooth surface leave $\text{RMSE} \approx 0.335$ of structural bias, which the residual sum of squares absorbs ($\sqrt{1^2 + 0.335^2} \approx 1.055$).
   - Tempered BASS (`new TemperingSpec()`) drops $\text{RMSE}$ to $0.282$, cuts $\hat{R}_M$ to $1.45$, raises $f$-coverage from $91.6\%$ to $97.2\%$, and covers $\sigma \in [0.991, 1.139]$. Keeping the pre-registered criterion untouched and showing that R `BASS` 1.3.1 fails the same way (while S8/Chip E tempering fixes $\sigma$ coverage) is the cleanest possible record.

2. **Slow local moves (1.5% birth/death acceptance, 5% change acceptance):** **Open as its own design item (`F-M1` / efficiency pass) rather than patching in Stage 4.**
   - Stage 4’s job was a faithful, verified foundation with proper priors, exact boundaries, and tempering.
   - When you open that design item, the three highest-leverage levers are:
     1. **Multiple RJ structure steps per sweep** before running $\tau$ / swaps / composition (since proposing and rejecting a birth/death takes only a rank-1 bordered Cholesky test once a column is evaluated, whereas swaps and composition pay per-sweep overhead).
     2. **Local knot/sign proposals inside `Change`** alongside full prior redraws, so existing terms can slide their knots without having to draw from $[0, 1]$ blind.
     3. **Data-informed / active-variable birth proposals** (with exact forward/reverse density ratios via `F-X4`'s proposal-kernel seam), since on Friedman 5 of 10 variables are inert and uniform/Wallenius births waste half their proposals on inert coordinates.

3. **Provisional run lengths (`BurnIn = 9000`, `Samples = 1000`):** **Leave to the efficiency pass.** Lengthening the defaults before improving local move efficiency only scales wall time against the same 1.5% birth/death bottleneck.

4. **Chip E’s six constants:** **Accept as built.**

5. **Self-sizing ladder settling counts:** **Leave to the efficiency pass** (see Finding A below for why `FunctionalData` over-settles slightly while `Friedman` stops at 20–23 rungs with one saturated pair).

---

### 2. Two Technical Findings from the S10 Runs

#### A. Why Tempered Friedman Stopped at 20–23 Rungs with One Saturated Pair ($r_{i^*} = 1.000$)
In [bass-questions.md](file:///D:/aghado01/ThermoMapper/issues/regression/bass-questions.md#L24), the chip noted:
> *"The tempered run settles at 20–23 rungs, but one pair still rejects every swap, and no complete round trip happened."*

Why didn't [`TemperingLadder.Size`](file:///D:/aghado01/ThermoMapper/src/maths/samplers/rjmcmc/TemperingLadder.cs#L88-L96) grow the ladder past 20–23 rungs toward `MaxLevels = 64` when a pair was still rejecting $100\%$ of swaps?
- **Right-censoring of $r_i$ at $1.0$ in $\hat{\Lambda} = \sum_i r_i$:**
  When *every* pair is saturated ($r_i \approx 1$), $\hat{\Lambda} \approx L - 1$, so $2\hat{\Lambda} + 1 \approx 2L - 1$ doubles the ladder.
  **But when the barrier is localized in one narrow temperature interval**—so **one** pair has $r_{i^*} = 1.000$ while the other $\sim 19$ pairs have already equalized to $r_i \approx 0.48$—that single saturated pair is capped at contributing at most $+1.0$ to $\hat{\Lambda} = \sum_i r_i$:
  $$\hat{\Lambda} = 1.0 + 19 \times 0.48 \approx 10.1 \implies \text{RecommendedLevels}(\hat{\Lambda}) = \lceil 2 \times 10.1 \rceil + 1 = 22$$
  Because $r_{i^*}$ cannot exceed $1.0$, $\hat{\Lambda}$ stops growing at $\sim 10$ and [`TemperingLadder.Place`](file:///D:/aghado01/ThermoMapper/src/maths/samplers/rjmcmc/TemperingLadder.cs#L119-L133) only assigns $1.0 / 10.1 \approx 10\%$ of the rungs (~2 intervals) across that cliff!
- Worth noting for the efficiency pass: if `max(r_i) > 0.95` while $L < \text{MaxLevels}$, either `Size` can add extra rungs for any right-censored pair ($r_i \to 1$), or `Place` can inflate the effective barrier increment of pairs near $1.0$ (e.g., via the Gaussian inverse $-2\Phi^{-1}((1 - r_i)/2)$) so new rungs concentrate where the cliff actually is.

#### B. Why `Psis_UnderSampledG_IsNoFalseNegative` Took So Long to Run
Look at [`BassQualificationLooTests.ExactLoo`](file:///D:/aghado01/ThermoMapper/tests/regression/BassQualificationLooTests.cs#L106-L153):
- For each of the 3 datasets (`"surface"`, `"outlier"`, `"spike"`) and each of the $|\mathcal{F}| = 2{,}263$ enumerated structures, it builds an array of $2{,}112$ quadrature nodes in $u = \ln\tau$ (all sharing the same structure $k$) and calls:
  [`BassPointwise.LogDensity(model, table, [draws], content: PointwiseContent.CollapsedPredictive)`](file:///D:/aghado01/ThermoMapper/tests/regression/BassQualificationLooTests.cs#L137)
- Inside [`BassPointwise.Collapsed`](file:///D:/aghado01/ThermoMapper/src/maths/regression/spline/bass/BassPointwise.cs#L135-L143), even though all $2{,}112$ draws share the same structure $k$, the inner loop calls [`conditional.LeaveOutLogDensity(...)`](file:///D:/aghado01/ThermoMapper/src/maths/regression/gprior/CollapsedGPrior.cs#L505-L600) from scratch for every node $i \in \{1 \dots 2112\}$ and unit $u \in \{1 \dots 27\}$.
- That executes $3 \times 2{,}263 \times 2{,}112 \times 27 \approx 3.87 \times 10^8$ calls to `LeaveOutLogDensity`, each repeating the $\tau$-independent Cholesky solve $G_c^{-1} Z_U^\top$ and heap-allocating scratch arrays. Hoisting the structure-only solve $G_c^{-1} Z_U^\top$ and OLS projections out of the per-$\tau$ loop in `CollapsedGPrior` / `BassPointwise.Collapsed` will cut both `ExactLoo` and normal `SampledG` `CollapsedPredictive` calls by more than an order of magnitude when draws share structures.