
Here is a read on the S7 slice in [bass-questions.md](file:///D:/aghado01/ThermoMapper/issues/regression/bass-questions.md) and its full report in [stage4-run-brief.md](file:///D:/aghado01/ThermoMapper/issues/regression/bass/stage4-run-brief.md#L82-L179).

---

### 1. The Three Open Questions at the Bottom of S7

All three were already adopted on `main` between S7 and S8, and the choices hold up well:

1. **`MonitoredPoints` ([`BassSampling.MonitoredPoints`](file:///D:/aghado01/ThermoMapper/src/maths/regression/spline/bass/BassConfig.cs#L297), [`BassDiagnostics.DefaultPoints`](file:///D:/aghado01/ThermoMapper/src/maths/regression/spline/bass/BassDiagnostics.cs#L55-L66)):**
   - **Accept as built.** Monitoring $M$, $\tau$, $\sigma^2$, plus the fitted surface $f(x_i, t_r)$ at quartile case indices (`0.25, 0.5, 0.75` of `cases - 1`, at the middle output location $\lfloor(m-1)/2\rfloor$ for functional data) gives rank-normalized $\hat{R}$ and ESS on both global hyperparameters and local function values without BARS's $O(\text{grid})$ functional vector.
   - Allowing a declared `IReadOnlyList<BassMonitoredPoint>?` on `BassSampling` (with `Location = null` resolved to the middle output location in [`BassResolved.ResolvePoints`](file:///D:/aghado01/ThermoMapper/src/maths/regression/spline/bass/BassResolved.cs#L200-L227)) is clean and keeps the config self-contained.
2. **`Resume`'s checkpoint type ([`Bass.Resume`](file:///D:/aghado01/ThermoMapper/src/maths/regression/spline/bass/Bass.cs#L31)):**
   - **Accept as built.** Taking [`EnsembleCheckpoint`](file:///D:/aghado01/ThermoMapper/src/maths/samplers/ensemble/Checkpoint.cs) directly matches [`BassFit.Checkpoint`](file:///D:/aghado01/ThermoMapper/src/maths/regression/spline/bass/BassFit.cs) and avoids a wrapper type that would only forward to `EnsembleCheckpoint`.
3. **`BassFit.Evidence`:**
   - Landed in S8 ([`BassEvidence`](file:///D:/aghado01/ThermoMapper/src/maths/regression/spline/bass/BassEvidence.cs)).

---

### 2. A Hidden Cost in S7's Design-Snapshot Departure (Now That S8 Uses It on Every Swap)

In S7, Chip A rightly noticed that reloading a structure from scratch via [`BassDesign.Load`](file:///D:/aghado01/ThermoMapper/src/maths/regression/spline/bass/BassDesign.cs#L263-L294) differs by a few ulps from an incrementally updated [`BorderedCholesky`](file:///D:/aghado01/ThermoMapper/src/maths/linalg/BorderedCholesky.cs) (because deaths and changes apply Givens rotations between periodic rebuilds). To keep `Resume` bit-identical, Chip A added [`BassDesign.Snapshot`](file:///D:/aghado01/ThermoMapper/src/maths/regression/spline/bass/BassDesign.cs#L301-L302) and [`BassDesign.Restore`](file:///D:/aghado01/ThermoMapper/src/maths/regression/spline/bass/BassDesign.cs#L306-L329).

Look closely at what `Restore` does, and how S8's [`BassLadder.HandOver`](file:///D:/aghado01/ThermoMapper/src/maths/regression/spline/bass/BassLadder.cs#L230-L250) now uses it:

```csharp
internal BassDesignSnapshot Snapshot()
    => new([.. _terms], _factor.Copy(), (double[])_beta.Clone(), (double[])_residual.Clone(), ResidualSquares, FittedSquares);

internal void Restore(BassDesignSnapshot snapshot)
{
    ...
    foreach (BassTerm term in snapshot.Terms)
    {
        if (Evaluate(term) is not var (column, output)) throw ...
        _terms.Add(term);
        _a.Add(column);
        if (_functional) _h.Add(output!);
        _b.Add(CrossResponse(column, output));
    }
    _factor = snapshot.Factor.Copy();
    _beta = (double[])snapshot.Beta.Clone();
    snapshot.Residual.CopyTo(_residual, 0);
    ResidualSquares = snapshot.ResidualSquares;
    FittedSquares = snapshot.FittedSquares;
}
```

Two things jump out when S7 and S8 are read together:

1. **In S8, [`BassLadder.HandOver`](file:///D:/aghado01/ThermoMapper/src/maths/regression/spline/bass/BassLadder.cs#L230-L250) calls `DesignFor` $\to$ `Snapshot()` and `Adopt` $\to$ `Restore()` on every accepted replica swap, every scan:**
   - Because `Restore` was written for checkpoint deserialization (where `_a`, `_h`, and `_b` are not saved), calling `Restore` during an in-memory swap **re-evaluates every basis column `Evaluate(term)` and its inner product `CrossResponse` from scratch** ($O(n_x M K + N M)$) on both the giver's and receiver's side of every swap!
   - Yet the giver replica's `BassDesign` *already holds* `_a` (`double[]` per term), `_h` (`double[]` per term), and `_b` (`double` per term), and those per-term column arrays are never mutated once created!
   - If `BassDesign` (or `BassDesignSnapshot`) carries or shares the immutable `_a` and `_h` column arrays and `_b` entries on in-memory hand-over, a replica swap drops all basis-column evaluations and `CrossResponse` dot products completely—copying only the $M$ array references and the $M \times M$ factor. **This will directly show up in S10's tempered cost baseline.**

2. **Conversely, for checkpoints, `_beta`, `_residual`, `ResidualSquares`, and `FittedSquares` are redundant with `_factor`:**
   - Look at [`BassDesign.Accept`](file:///D:/aghado01/ThermoMapper/src/maths/regression/spline/bass/BassDesign.cs#L235-L257) and [`BassDesign.Refit`](file:///D:/aghado01/ThermoMapper/src/maths/regression/spline/bass/BassDesign.cs#L531-L538): on every accepted move, `Accept()` calls `Refit()`, which computes `_beta = _factor.Solve(b)`, `_residual = Fit(_a, _h, _beta)`, `ResidualSquares = PairwiseSum.Dot(_residual, _residual)`, and `FittedSquares = PairwiseSum.Dot(b, _beta)` from scratch.
   - Only `_factor` has path-dependent state between rebuilds; `_beta` and `_residual` are deterministic functions of `(_factor, _a, _h, _b)`.

---

### 3. Two Small Cleanups Left Open in Chip A's Report

In [stage4-run-brief.md](file:///D:/aghado01/ThermoMapper/issues/regression/bass/stage4-run-brief.md#L175-L179), Chip A noted two small cross-file items it left untouched because they lived in S1/S3 files:

1. **Duplicate `DrawSubset` in [`BassTermProposal`](file:///D:/aghado01/ThermoMapper/src/maths/regression/spline/bass/BassTermProposal.cs#L244-L273):**
   [`BassTermProposal.DrawSubset`](file:///D:/aghado01/ThermoMapper/src/maths/regression/spline/bass/BassTermProposal.cs#L244-L273) is a verbatim 30-line copy of [`BassStructurePrior.Part.DrawSubset`](file:///D:/aghado01/ThermoMapper/src/maths/regression/spline/bass/BassStructurePrior.cs#L333-L363). Having `BassTermProposal` call `part.DrawSubset(rng, variable, Prior.Subsets)` removes the duplicate RNG consumer.
2. **Moving the two-chain check for `RHatTarget` into [`BassConfig.Validate`](file:///D:/aghado01/ThermoMapper/src/maths/regression/spline/bass/BassConfig.cs#L131-L132):**
   Currently, `if (sampling.RHatTarget > 0.0 && sampling.Chains < 2)` is checked in [`Bass.Run`](file:///D:/aghado01/ThermoMapper/src/maths/regression/spline/bass/Bass.cs#L47-L48) *after* `BassModel.Create(config, data)` has already resolved the data and built the basis. Moving that check into [`BassConfig.Validate`](file:///D:/aghado01/ThermoMapper/src/maths/regression/spline/bass/BassConfig.cs#L131-L132) lets `BassConfig.FromJson` and `BassResolved.Resolve` reject `RHatTarget > 0` with `Chains = 1` immediately alongside the other `Sampling` checks.

--- 

Here is how S9 ([stage4-density-brief.md](file:///D:/aghado01/ThermoMapper/issues/regression/bass/stage4-density-brief.md), commit `ab22bc6`) stands in the codebase and my recommendations on its three questions, plus two performance catches for S10.

---

### 1. Several Output Variables in [`BassResidualDependence`](file:///D:/aghado01/ThermoMapper/src/maths/regression/spline/bass/BassResidualDependence.cs#L121-L156)

Currently [`BassResidualDependence.Of`](file:///D:/aghado01/ThermoMapper/src/maths/regression/spline/bass/BassResidualDependence.cs#L125) returns `null` whenever `model.Inputs.Outputs.Count != 1`.

The two-pronged design recorded in [port-fitting-brief.md](file:///D:/aghado01/ThermoMapper/issues/regression/bass/port-fitting-brief.md#L547) is the right call:
1. **1D ordering declared on a multi-column output grid** (e.g., $(t, \text{elevation})$ where only $t$ is the index coordinate and elevation is a location covariate): let the caller pass either the output variable index/name or an explicit location permutation into [`BassResidualDependence.Of`](file:///D:/aghado01/ThermoMapper/src/maths/regression/spline/bass/BassResidualDependence.cs#L31-L34) and compute lag-1 autocorrelation $r_1$ along that ordering.
2. **$d \ge 2$ spatial coordinates** (e.g., $(s_1, s_2)$ in a 2D field): compute per-case Moran's $I$ over each location's $k$ nearest neighbours in the scaled coordinates $[0, 1]^d$:
   $$I_i = \frac{m}{\sum_{r,s} w_{rs}} \frac{\sum_{r,s} w_{rs}(e_{ir} - \bar{e}_i)(e_{is} - \bar{e}_i)}{\sum_r (e_{ir} - \bar{e}_i)^2}$$

**Two implementation details worth settling when S10 builds this:**
- **Row-standardize the $k$NN weights** ($w_{rs} = 1/k$ for $s \in \mathcal{N}_k(r)$, so $\sum_s w_{rs} = 1$ and $\sum_{r,s} w_{rs} = m$). Then the prefactor $m / \sum_{r,s} w_{rs}$ is identically $1$, $I_i \in [-1, 1]$, and its exact null expectation under independent white noise is $-1/(m-1)$—matching lag-1 $r_1$'s null behaviour ($\approx -1/m$) on the exact same scale.
- **Record shape and automatic assembly:** [`ResidualDependence`](file:///D:/aghado01/ThermoMapper/src/maths/regression/spline/bass/BassResidualDependence.cs#L19) currently names its per-case vector `LagOne`. If `ResidualDependence` also holds Moran's $I$, either add a `Neighbours` / `Kind` discriminator (or rename/alias `LagOne` before external callers lock onto it). Also decide what [`BassEnsemble.Assemble`](file:///D:/aghado01/ThermoMapper/src/maths/regression/spline/bass/BassRun.cs#L170) should do automatically when `Outputs.Count > 1`: either compute $k$NN Moran's $I$ over all output columns by default (assuming all `Outputs` are coordinates unless told otherwise) or leave `fit.ResidualDependence` as `null` when `Outputs.Count > 1` and let callers invoke `BassResidualDependence.Of(fit, ...)` with their coordinate subset. Running $k$NN Moran's $I$ ($k = \min(4, m-1)$) automatically over all `Outputs` is very cheap ($O(m^2 d)$ once for the shared location grid) and gives every multi-coordinate fit an out-of-the-box diagnostic.

---

### 2. Centring of $r_1$ in [`BassResidualDependence.Of`](file:///D:/aghado01/ThermoMapper/src/maths/regression/spline/bass/BassResidualDependence.cs#L53-L63)

**Recommendation: Accept as built (centred by each case's own mean residual $\bar{e}_i$).**

In a direct functional regression $Y_{ir} = f(x_i, t_r) + \varepsilon_{ir}$, any predictor-space misfit $f(x_i, \cdot) - \hat{f}(x_i, \cdot)$ that shifts curve $i$ up or down across all $m$ locations is **between-case ($x$-space) bias**, not within-curve serial correlation along $t$.
- If $r_1$ were computed without subtracting $\bar{e}_i$, a constant offset $c_i$ on case $i$ with white location noise $\varepsilon_{ir} \sim \mathcal{N}(0, \sigma^2)$ would produce
  $$r_{1,i}^{\text{uncentred}} \approx \frac{c_i^2}{c_i^2 + \sigma^2} \to 1$$
  falsely flagging pure $x$-space misfit as within-curve autocorrelation along $t$.
- Centring by $\bar{e}_i$ isolates shape misfit / serial correlation along $t$, and gives the standard textbook null expectation $-1/m$ and variance $\approx 1/m$.

---

### 3. D04's Brute-Force Leave-One-Out Check of `CollapsedPredictive`

In [`BassPointwise.Collapsed`](file:///D:/aghado01/ThermoMapper/src/maths/regression/spline/bass/BassPointwise.cs#L98-L147), `PointwiseContent.CollapsedPredictive` evaluates
$$\ell_{s,u} = \ln p(y_u \mid y_{-u}, S^{(s)}, \tau^{(s)})$$
via [`GPriorConditional.LeaveOutLogDensity`](file:///D:/aghado01/ThermoMapper/src/maths/regression/gprior/CollapsedGPrior.cs#L505-L600).

**Why `PsisLoo` on `CollapsedPredictive` is mathematically valid:**
Even though $y_u$ is not conditionally independent of $y_{-u}$ once $(a, \sigma^2)$ are integrated out, Bayes' rule on the collapsed target gives:
$$p(S, \tau \mid y) = p(S, \tau \mid y_{-u}) \, \frac{p(y_u \mid y_{-u}, S, \tau)}{p(y_u \mid y_{-u})} \implies \frac{p(S, \tau \mid y_{-u})}{p(S, \tau \mid y)} \propto \frac{1}{p(y_u \mid y_{-u}, S, \tau)} = \exp(-\ell_{s,u})$$
So feeding $\ell_{s,u} = \ln p(y_u \mid y_{-u}, S^{(s)}, \tau^{(s)})$ into `PredictiveChecks.PsisLoo` forms the exact Rao-Blackwellized importance ratios between the full collapsed posterior $p(S, \tau \mid y)$ and the leave-unit-out collapsed posterior $p(S, \tau \mid y_{-u})$.

**Why Chip B saw high $\hat{k}$ on scalar hinge data anyway (9 of 60 units with $\hat{k} > 0.62$):**
Integrating out $(a, \sigma^2)$ removes coefficient/scale variance, **not** knot-support sensitivity. In a 1D scalar hinge model with `SupportMinimum = 3`, a hinge near the boundary of $[0, 1]$ has only a handful of active observations on its non-zero side. Dropping one of those observations changes the local leverage $h_{uu} = z_u^\top G_c^{-1} z_u$ and the Woodbury matrix $I - (\tau / (\tau + 1)) Z_U G_c^{-1} Z_U^\top$ sharply for draws carrying that boundary knot. By contrast, on functional data ($m = 12$ locations per case), each predictor case is replicated across $m$ locations and the $m \times m$ block update is much smoother ($\hat{k} \le 0.45$ in Chip B's test).

When S10 adds the brute-force $K$-unit refit check (using [`BassChain`](file:///D:/aghado01/ThermoMapper/src/maths/regression/spline/bass/BassChain.cs#L43-L61) or `Bass.Fit` on $y_{-u}$), expect:
- Low-$\hat{k}$ units ($\hat{k} < 0.5$) under `CollapsedPredictive` to agree closely with brute-force refits on $y_{-u}$, and with smaller Monte Carlo variance than `RetainedParameters`.
- High-$\hat{k}$ boundary units to show why $\hat{k}$ warned: dropping $y_u$ actually shifts which boundary knots survive `SupportMinimum` in the true leave-one-out posterior $p(S \mid y_{-u})$.

---

### 4. Two Performance Catches in S9 Code

1. **Per-term posterior mean in [`BassResidualDependence.Of`](file:///D:/aghado01/ThermoMapper/src/maths/regression/spline/bass/BassResidualDependence.cs#L133-L145) (already noted for S10 in [port-fitting-brief.md](file:///D:/aghado01/ThermoMapper/issues/regression/bass/port-fitting-brief.md#L547)):**
   `BassResidualDependence.Of` runs automatically inside [`BassEnsemble.Assemble`](file:///D:/aghado01/ThermoMapper/src/maths/regression/spline/bass/BassRun.cs#L170) at the end of every functional `Bass.Fit`. Right now it loops over all $S$ retained draws and calls [`fit.Evaluate(pooled[s], f)`](file:///D:/aghado01/ThermoMapper/src/maths/regression/spline/bass/BassResidualDependence.cs#L138), evaluating the $n_x \times n_r$ outer product for every draw ($O(S \cdot M \cdot n_x \cdot n_r)$).
   Because the surface is linear in the coefficients $a$ and every basis function is interned in [`BassFit.Terms`](file:///D:/aghado01/ThermoMapper/src/maths/regression/spline/bass/BassFit.cs#L65), you can accumulate the mean intercept $\bar{a}_0$ and a `double[fit.Terms.Count]` array of posterior-mean coefficients $\bar{a}_{\text{id}} = \sum_s w_s a_{\text{id}}^{(s)}$ in $O(S \cdot M)$ scalar additions first, and then evaluate each distinct term in `fit.Terms` **once** ($O(|\mathcal{T}| \cdot n_x \cdot n_r)$).

2. **Redundant Cholesky solves across draws sharing a structure in [`BassPointwise.Collapsed`](file:///D:/aghado01/ThermoMapper/src/maths/regression/spline/bass/BassPointwise.cs#L118-L144):**
   `BassPointwise.Collapsed` groups draws by `Structure` and builds `design` and `x = fit.CentredDesign(ids)` once per distinct structure. However, inside the per-draw loop (`foreach (int s in rows)`), it calls [`conditional.LeaveOutLogDensity(...)`](file:///D:/aghado01/ThermoMapper/src/maths/regression/gprior/CollapsedGPrior.cs#L505-L600) for every draw $s$ and unit $u$.
   Inside [`GPriorConditional.LeaveOutLogDensity`](file:///D:/aghado01/ThermoMapper/src/maths/regression/gprior/CollapsedGPrior.cs#L518-L536):
   - `_precision` is the unweighted centred Gram Cholesky factor $G_c = Z^\top Z$, which **depends only on the structure $S$, not on $\tau^{(s)}$**.
   - Yet for every draw $s$ in `rows`, `LeaveOutLogDensity` re-runs `_precision.Solve` on every held-out row $Z_U$ and re-computes the $m \times m$ matrix $H_U = Z_U G_c^{-1} Z_U^\top$ and the OLS projections $Z_U \hat{\beta}_{\text{OLS}}$ from scratch.
   - When multiple draws share a structure (or when `FixedG` / `EBLocal` makes $\tau$ identical across all draws in `rows`!), $H_U = Z_U G_c^{-1} Z_U^\top$ and $r_U = (y_U - \bar{y}) - \text{shrink} \cdot Z_U \hat{\beta}_{\text{OLS}}$ only depend on $\tau$ through the scalars $\text{shrink} = \tau / (\tau + 1)$ and $S(\tau)$ (and if $\tau$ is identical across two draws with the same structure, the entire vector of $U$ leave-out log densities $\ell_{s,1\dots U}$ is identical and can simply be copied!).

   --- 


Chip E’s work (`ec565f6`, `d64ee4e`, `66c66be`, `abf16cf`, `36ce26a`, `c388461`) is very clean and resolves the S8 functional evidence discrepancy completely.

---

### 1. Code & Empirical Verification

- **The S8 functional fixture mystery is solved:** S8’s 16- and 32-rung readings of $180\text{–}186$ were short-run stepping-stone downward bias from under-sampled importance weights. At a 4,000 burn-in / 4,000 sample budget on `FunctionalData` ($N = 175$), all four estimators now agree:
  - Self-sized integrated (11–15 rungs): $192.59 \pm 0.52$
  - Self-sized linear (12–16 rungs): $193.45 \pm 0.41$ (path sampling: $193.42 \pm 0.53$)
  - Fixed 32 rungs, integrated: $193.40 \pm 0.37$
  - Fixed 64 rungs, integrated: $193.03 \pm 0.30$
- **Shared sampler invariants hold:**
  - In [`ParallelTempering`](file:///D:/aghado01/ThermoMapper/src/maths/samplers/rjmcmc/ParallelTempering.cs#L137-L143), `_growthRng` is appended as `seeds[^1]`, so any ladder that never calls [`Resize`](file:///D:/aghado01/ThermoMapper/src/maths/samplers/rjmcmc/ParallelTempering.cs#L482-L552) is bit-identical to pre-Chip-E (`BARS` and the 7 pinned fixed-`Levels` BASS runs are unchanged bit for bit).
  - [`LadderEvidence.Pool`](file:///D:/aghado01/ThermoMapper/src/maths/samplers/evidence/LadderEvidence.cs#L355-L379) reports $\max(\text{within}, \text{between})$ for $C \ge 2$ and exposes `WithinChainError`, `BetweenChainError`, and `BetweenToWithin` on both [`PooledLadderEstimate`](file:///D:/aghado01/ThermoMapper/src/maths/samplers/evidence/LadderEvidence.cs#L208-L215) and [`BassEvidenceEstimate`](file:///D:/aghado01/ThermoMapper/src/maths/regression/spline/bass/BassEvidence.cs#L87-L94).
- **Clone-safe replicas in [`BassLadder.Resize`](file:///D:/aghado01/ThermoMapper/src/maths/regression/spline/bass/BassLadder.cs#L147-L186):**
  - [`Fork`](file:///D:/aghado01/ThermoMapper/src/maths/regression/spline/bass/BassLadder.cs#L304) clones `TermArray` so no two replicas share a mid-run array reference, preventing a false `SharesStructure` match or a mid-run whole reload that could disagree with incremental Cholesky updates near D06’s angle floor.
  - If a resize ever forks the reference rung’s state (which is drawn via [`BassPriorSampler`](file:///D:/aghado01/ThermoMapper/src/maths/regression/spline/bass/BassPriorSampler.cs) without loading a design into `_structures[^1]` until needed), [`DesignFor`](file:///D:/aghado01/ThermoMapper/src/maths/regression/spline/bass/BassSweep.cs#L125-L129) syncs and loads it once on demand, covered by a dedicated test.

---

### 2. Recommendations on Chip E’s Open Questions

All six defaults Chip E chose can be **accepted as-is**:

1. **Initial count (`InitialLevels = 8` in [`TemperingSpec`](file:///D:/aghado01/ThermoMapper/src/maths/regression/spline/bass/BassConfig.cs#L329)):** **Accept.** Starts near where small scalar/categorical problems settle (4–7 rungs) and grows smoothly on larger $N$.
2. **Minimum count (`MinimumLevels` = 2, or 3 with reference rung) and the Round-1 dip:** **Accept.**
   - *Why Round 1 dips:* In [`BassChainRun.Burn`](file:///D:/aghado01/ThermoMapper/src/maths/regression/spline/bass/BassRun.cs#L335-L342), the first round is `round = 1` scan. Because all non-reference replicas start at the same state `cold.Start` and cold RJ moves accept only ~5–15% of the time, most adjacent pairs are still at the exact same state after 1 scan ($\log \alpha = 0 \implies r_i = 0$), so $\hat{\Lambda}$ under-reads on scan 1 and can shrink the ladder before regrowing on rounds $2, 4, 8, \dots$.
   - With normal burn-in lengths ($\ge 500$ sweeps, i.e., 9+ doubling rounds), the ladder regrows completely with negligible cost. Keeping it as landed preserves all of Chip E’s pinned sizing tests. *(If you ever want short burn-ins like 10–50 sweeps not to dip on scan 1, the one-line tweak is to only allow shrinking once `round >= 4`, without changing normal runs.)*
3. **One-rung deadband (`ShrinkDeadband = 1` in [`TemperingLadder.Size`](file:///D:/aghado01/ThermoMapper/src/maths/samplers/rjmcmc/TemperingLadder.cs#L46)):** **Accept.** Prevents oscillation when $2\hat{\Lambda} + 1$ hovers near an integer boundary while still allowing genuine downward adjustment after transient burn-in spikes.
4. **Saturation flag threshold (`RejectionThreshold = 0.9` in [`BassEvidence`](file:///D:/aghado01/ThermoMapper/src/maths/regression/spline/bass/BassEvidence.cs#L43)):** **Accept.** Under the Gaussian log-ratio model, relative weight variance is $e^{\sigma^2/2} - 1 \approx 220$ at $r = 0.9$ and $\approx 2{,}200$ at $r = 0.95$. Since `Saturated` is a non-destructive diagnostic flag (evidence is never nulled), $0.9$ is the right place to warn that stepping-stone weights are getting heavy-tailed.
5. **New rung’s source (`TemperingLadder.Sources` in [`TemperingLadder.cs`](file:///D:/aghado01/ThermoMapper/src/maths/samplers/rjmcmc/TemperingLadder.cs#L142-L173)):** **Accept.** Picking the nearest current rung on the cumulative barrier $\Lambda(\beta)$ (tie-broken by $|\Delta\beta|$, then colder) matches thermodynamic overlap rather than raw $\beta$-distance, which matters near the steep $b \to 0$ drop.
6. **Kept-replica assignment in [`ParallelTempering.Resize`](file:///D:/aghado01/ThermoMapper/src/maths/samplers/rjmcmc/ParallelTempering.cs#L500-L513):** **Accept.** Assigning each named source replica $s$ to the naming rung $k$ closest to $\beta_s$ (colder on a tie) automatically keeps the cold replica at index $0$ ($\Delta\beta = 0$) and the reference replica at `levels - 1` ($\Delta\beta = 0$) without special-casing, while minimizing retempering jumps on kept interior replicas.
7. **Burn-in barrier over-reading ($\hat{\Lambda} = 5\text{–}7.5$ in burn-in vs $3.8$ in sampling) & low round-trip rate on `FunctionalData`:**
   - **Why burn-in reads higher:** Right after a resize, newly forked replicas share their source rung's state, and while the cold chain is climbing toward the cap $M_{\max} = 3$ from the prior draw, the transient likelihood gap across rungs is wider than at stationarity. Settling a few rungs above $2\hat{\Lambda}_{\text{sampling}} + 1$ (11–17 rungs vs ~9) is conservative and keeps max rejection around $0.84\text{–}0.90$.
   - **Why round trips per scan are $10^{-4}\text{–}10^{-3}$ vs Theorem 1's ideal $\sim 0.03$:** That gap is Syed et al.’s diagnostic working as designed: once the ladder is well-spaced, the round-trip bottleneck on `FunctionalData` at $M = M_{\max} = 3$ is **local structure persistence** (at the cap, birth moves are unavailable, so a state can only change basis functions via death or change proposals), not the temperature ladder.

---

### 3. Remaining Follow-Ups Noted by Chip E

Only the two documentation updates in your S10 scope remain:
1. Update `docs/Maths/Samplers.md` with the new shared members (`ParallelTempering.Resize`, `RoundRejection`, `PairScans`, `TraceRungs(bool traceSlopes = true)`; `TemperingLadder.Place`, `Sources`, `Size`, `ShrinkDeadband`; `LadderEvidence.Pool`, `PooledLadderEstimate`).
2. Update the BASS tempering/evidence docs and S10 cost baseline with the self-sizing ladder, auto-picked `HotEnd`, per-chain `BassTemperingReadout`, and `BassEvidence` flags.

---

