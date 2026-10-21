# Samplers

> A model-agnostic Markov chain Monte Carlo (MCMC) toolkit covering a model-neutral target contract with declared tempering paths, reversible-jump MCMC with pluggable moves, parallel tempering, and a replica-ensemble shell that runs independent chains in parallel and diagnoses convergence on scalar functionals. The project also contains inverse-transform sampling and a power-law Gaussian noise-field generator. · **Project:** `Maths.Samplers` (namespaces `Maths.Samplers`, `Maths.Samplers.Rjmcmc`, `Maths.Samplers.Ensemble`, and others) · **Maturity:** core

## Overview

The samplers sit in layer 2 of the [architecture](./Architecture.md). They depend on `Maths.Rng`, and on `Maths.LinAlg` for the FFT. Their consumers are:

- [BayesianRegression](./BayesianRegression.md). BARS drives the reversible-jump engine and the ensemble shell. The Bayesian GP bandwidth sampler (BGP) also uses the ensemble shell. BAPS uses only the R̂ diagnostic.
- `InverseTransform`, which an SPC demonstration (`WangMonteCarloDemo`) uses for exponential energy draws. The Potts bond kernel `FkKernel` evaluates the matching CDF $1-e^{-x}$ but does not call the sampler.
- `SpectralNoiseField`, which [synthetic datasets](./SyntheticDatasets.md) use for colored-noise backgrounds.

```mermaid
flowchart LR
  subgraph Targets
    TG["ISamplingTarget&lt;T&gt;"] -->|Evaluate| EV["ITargetEvaluation"]
    TT["Temperature"] -.reads.-> EV
  end
  subgraph Rjmcmc
    RJ["ReversibleJumpChain&lt;T&gt;"] --> MV["IRjMove&lt;T&gt; / Proposal&lt;T&gt;"]
    SM["ScheduledMove&lt;T&gt;"] -.decorates.-> MV
    PT["ParallelTempering&lt;T&gt;"] -->|L replicas| RJ
  end
  RJ --> TG
  RJ -->|samples at| TT
  subgraph Ensemble
    CE[ChainEnsemble] -->|drives| IC["IChain&lt;T&gt;"]
    CE -->|StartChain / Burn / Accumulate| EM["IEnsembleModel&lt;T&gt; / IChainRun&lt;T&gt;"]
    CE --> CD[ChainDiagnostics]
  end
  RJ -.is an.-> IC
  PT -.is an.-> IC
```

## Key concepts

- **Carrier.** The state a chain walks. For BARS this is a knot configuration of varying dimension.
- **Move and palette.** A move is one proposal type, such as birth, death or relocate. The palette is the set of moves a chain chooses from. Each move names its **reverse** move by key.
- **Target and evaluation.** A target evaluates a carrier state once into an `ITargetEvaluation`: an untempered log prior and a log likelihood $\ell_b$ at any inverse temperature $b\in[0,1]$. A log prior of −∞ rejects a state without any likelihood work.
- **Temperature.** An inverse temperature and the path it is taken along: the **powered** posterior $b(\log p+\ell_1)$ or the **likelihood-tempered** path $\log p+\ell_b$. Both are the log posterior at $b=1$. A chain samples its target at one temperature.
- **Kernel vs shell.** The *kernel* (`IChain<T>`) advances one replica. The *shell* (`ChainEnsemble`) fans out replicas, runs rounds and pools diagnostics.
- **Functional.** A fixed-width scalar reduction of each draw, such as the knot count or fitted grid values. Convergence diagnostics are computed on functionals, not on raw states. R̂ is meaningless on a state whose dimension varies or whose labels are symmetric.

## Design

**Green's acceptance ratio, split by responsibility.** Only the engine sees the whole palette and each state's total move weight, so it owns two factors: the **target ratio** and the **move-selection ratio**. Only a move knows its bijection, so it owns the other two: its **auxiliary-density ratio** and its **Jacobian**. `Proposal<T>` carries exactly those two move-side terms. This split has two consequences:

- State-dependent selection weights can be a pure decorator (`ScheduledMove`, used by the DMGK schedule in BARS).
- Tempering acts on the target term alone: the chain's `Temperature` turns each evaluation into a log density. Tempering then leaves the Hastings terms untempered, as it should.

**The target contract (step 3.1).** A target returns an evaluation object, not a number.
- **Why the split.** The likelihood-tempered path keeps the prior whole, so it needs the prior and the likelihood apart. A collapsed model, such as BARS's Normal marginal with β and σ² integrated out, has an $\ell_b$ that is not $b\,\ell_1$, so "β times the log posterior" cannot express it.
- **Why an object.** The state's expensive work (for BARS, the design and its factorization) happens once, in `Evaluate`. Another temperature is then a query on the result, which is what a swap between rungs on the likelihood path, an SMC reweighting or a thermodynamic-integration estimate needs. A swap hands each replica the other's evaluation, and nothing is re-evaluated.
- **Constants.** Values are needed only up to an additive constant that may depend on $b$ and the data but never on the state. It cancels in Metropolis ratios, in swaps (each rung's constant enters once with each sign) and in normalized SMC weights. An evidence estimate needs complete values.
- **Refusals.** A model may refuse $b\ne1$ (`NotSupportedException`); the powered path asks only for $\ell_1$. A model may refuse $b=0$ when a nuisance prior is improper (`ArgumentOutOfRangeException`).
- **Alternatives set aside.**
  - A target that returns $\log\pi$ alone, with tempering as a decorator ($\beta\log\pi$, the `IRjTarget`/`TemperedTarget` pair it replaces): it cannot express the likelihood path or a nonlinear $\ell_b$.
  - Methods `LogPrior(state)` and `LogLikelihood(state, b)` on the target: each temperature would re-evaluate the state.
  - A closed-form tempered density only, `LogDensity(state, b)`: it hides the prior/likelihood split that SMC and evidence estimates read.
- **Costs.** One evaluation object per proposal, and each chain keeps its current state's evaluation alive (for BARS, the factor of $G$). Both are small beside the design and factor the model already builds per proposal.

**Conventions.**

- A `null` proposal means the move is unavailable here. The chain stays and does not count an attempt.
- A non-positive weight removes a move from the palette at that state.
- An infeasible state has a log prior of −∞ (no likelihood work is done) or a log likelihood of −∞ at every temperature.
- Log densities are defined only up to constants that do not depend on the state; a likelihood's may depend on $b$ and the data.

**Kernel/shell split.**

- `IChain<T>` is the kernel seam.
- `IEnsembleModel<T>` and `IChainRun<T>` form the client seam: the start policy, the burn hook, and the reduction from draw to functional.
- `ChainEnsemble` owns seeding, the parallel fan-out, rounds, early stopping and pooled diagnostics.
- Parallelism is across replicas only (`Parallel.For` over chains).

**Reproducible seeding.** One master seed is expanded once by `SeedTree.Derive` into 2C seeds: C *kernel* streams and C separate *readout* streams (see [ProbabilityAndInformation](./ProbabilityAndInformation.md)). A tempered chain re-derives L+1 sub-streams from its own chain seed.

## API at a glance

| Area | Entry points |
|---|---|
| Targets | `ISamplingTarget<T>` (`Evaluate`); `ITargetEvaluation` (`LogPrior`, `LogLikelihood(b)`); `Temperature(b, path)` with `LogDensity(evaluation)`, static `Posterior` and `LogPosterior(evaluation)`; `TemperingPath` (`Powered`, `Likelihood`); `SamplingTarget.FromLogDensity`, `SamplingTarget.Separable`, and the extension `target.LogPosterior(state)`; `SeparableEvaluation` |
| Reversible jump | `ReversibleJumpChain<T>(moves, target, start, rng, temperature = null)` with `Step()`, `Current`, `Temperature`, `CurrentEvaluation`, `CurrentLogDensity`, `Attempts`, `Accepted`, `MoveStats()`, `SetState(state, evaluation)`, `ReplaceTarget(target)`; `IRjMove<T>` (`Key`, `ReverseKey`, `Weight`, `Propose`); `Proposal<T>`; `MoveStat`; `ScheduledMove<T>` |
| Tempering | `ParallelTempering<T>(movesForReplica, target, start, betas, path, masterSeed)` with `Step()`, `Path`, `Cold`, `ColdChain`, `Replicas`, `SwapAttempts`, `SwapAccepts`, and static `GeometricLadder(levels, betaMin)` |
| Ensemble | `new ChainEnsemble().Run<T>(model, chains, masterSeed, burn, samples, rHatTarget = 0, batchSize = 0)` (an instance method on a stateless class) → `(EnsembleRun, chain runs)`; `ChainDiagnostics.RHat`, `ChainDiagnostics.Ess`; `RunningMoments` (`Add`, `Merge`, `Of`); `IChain<T>`, `IEnsembleModel<T>`, `IChainRun<T>`, `EnsembleRun` |
| i.i.d. and fields | `InverseTransform.Sample(rng, quantile, count)`, `InverseTransform.Exponential(rng, scale, count)`; `SpectralNoiseField.Generate(rng, gridSize, beta, boxMin, boxMax)` with `Sample`, `Min`, `Max` and presets `White` / `Pink` / `Brown` |

## Algorithms

### Reversible-jump step

At state $x$, let $w_i(x)$ be the move weights and $W(x)=\sum_{i:\,w_i(x)>0}w_i(x)$.

1. If $W(x)\le0$, stay. Otherwise select move m with probability $w_m(x)/W(x)$.
2. If `Propose` returns null, stay. This is not counted as an attempt. Otherwise count an attempt.
3. If the reverse move's weight at $x'$ is non-positive, or $W(x')\le0$, reject.
4. Otherwise accept with

$$\log\alpha=\big[\log\pi_b(x')-\log\pi_b(x)\big]+\underbrace{\log\tfrac{q_{m'}(u')}{q_m(u)}}_{\texttt{LogProposalRatio}}+\underbrace{\log|J|}_{\texttt{LogJacobian}}+\log\frac{w_{m'}(x')/W(x')}{w_m(x)/W(x)}.$$

Here $\log\pi_b$ is the chain's `Temperature.LogDensity` of a state's evaluation: $b(\log p+\ell_1)$ on the powered path and $\log p+\ell_b$ on the likelihood path, the log posterior at $b=1$, the default.

The chain accepts if $\log\alpha\ge0$ or $U<e^{\log\alpha}$. A candidate whose density is −∞ or NaN is rejected whatever its proposal ratio, so $-\infty+\infty$ cannot decide a move. The current evaluation and its density are cached, so each step makes at most one evaluation (none when the move is unavailable or the reverse move is blocked).

**What the constructor checks.** There must be at least one move, the keys must be unique, and every `ReverseKey` must exist in the palette.

**What guarantees detailed balance.** The engine does not check that a move's reverse points back to it. Detailed balance relies on three things:

- the moves being paired symmetrically;
- each move reporting correct auxiliary ratios and Jacobians;
- the selection ratio, which the engine computes.

**Cost.** Each step costs $O(M)$ weight evaluations for M moves, plus at most one proposal and one target evaluation.

### Tempering

**Paths (step 3.2).** A ladder declares one `TemperingPath` for all its replicas, and replica $i$ samples `Temperature(β_i, path)`.
- **Powered:** $\beta_i(\log p+\ell_1)$. Prior and likelihood are flattened together.
- **Likelihood-tempered:** $\log p+\ell_{\beta_i}$, the power posterior of Friel and Pettitt (2008). The prior stays whole, and $\beta=0$ is the prior itself.
- **What separates them.** On a space of varying dimension, a powered density depends on the reference measure of each subspace, so its hot rungs leave the prior in a direction that measure sets, not the model. For BARS's ordered knots under a Poisson(λ) count prior and a constant likelihood, the powered rung β samples a Poisson(λ^β) knot count, fewer knots as β falls; every likelihood-path rung samples the prior (`CountLawTests`, per rung). The likelihood path needs $\ell_b$ at every rung; the powered path needs only $\ell_1$, so it serves every target.
- **Choice.** The caller declares the path; no default lives in the sampler. BARS takes the likelihood path ([BayesianRegression](./BayesianRegression.md), Q5).

**The ladder.** In `ParallelTempering`, every replica starts from the same state, and each builds its own move palette from the `movesForReplica` factory, so mutable move state stays private to one replica. Each replica has its own RNG, and one further RNG drives swaps. `GeometricLadder` gives $\beta_i=\beta_{\min}^{\,i/(L-1)}$. `Replicas` exposes every rung, cold first, for per-rung diagnostics.

**Ladder validation.** At least two rungs, strictly descending, cold rung exactly $\beta_0=1$. The hot end must be positive on the powered path, where $\beta=0$ would be the bare reference measure. It may be 0 on the likelihood path, where $\beta=0$ is the prior; a target whose nuisance prior is improper refuses that rung when it is evaluated (BARS does).

**One step.**

1. Every replica steps once.
2. A deterministic sweep from cold to hot attempts a swap for each adjacent pair $(i,i+1)$, with

$$\log A=d_i(x_{i+1})+d_{i+1}(x_i)-d_i(x_i)-d_{i+1}(x_{i+1}),$$

where $d_i$ is replica $i$'s log density, its `Temperature.LogDensity`. The two current densities are cached, and the two cross terms query the other replica's evaluation at this replica's temperature, so a swap evaluates no state. An accepted swap hands each replica the other's state and evaluation, and each recomputes its density at its own β.
- **On the powered path** this is $(\beta_i-\beta_{i+1})(L_{i+1}-L_i)$ with $L=\log p+\ell_1$, the formula used until step 3.2; BARS's powered runs are bit-identical across the change.
- **On the likelihood path** the prior cancels, leaving $\ell_{\beta_i}(x_{i+1})+\ell_{\beta_{i+1}}(x_i)-\ell_{\beta_i}(x_i)-\ell_{\beta_{i+1}}(x_{i+1})$. The β-linear shortcut is wrong here, because $\ell_\beta$ need not be linear in β. On a six-state test with a nonlinear $\ell_\beta$, the shortcut biases even the cold rung (a state's frequency 0.0114 against the exact 0.0223). The general form matches every rung's law on both paths to total variation ≤ 0.005. Later pairs in the same sweep see earlier swaps. The chain interface reports the cold replica's move acceptance. Swaps are counted separately.

**The swap uniform is drawn on every attempt**, accepted or not.
- Between two states whose $L$ agree to rounding, $\log A$ is of rounding size and its sign is rounding noise. Until step 3.1, $L_i$ was recovered as `CurrentLogTarget`$/\beta_i$ after an earlier swap stored $\beta_i L$, so replicas holding the same state produced such ties routinely. $L_i$ is now read from the evaluation, so the same state gives the same $L$ exactly, but distinct states can still tie. The regression test builds two targets that break such a tie in opposite directions.
- Drawing every time keeps the swap stream independent of it, so a last-bit change in the target cannot shift later swaps.
- Before 2026-10-20 the draw was skipped when $\log A\ge0$, and step 2.3's cancellation-free scale moved a tempered BARS run through exactly this tie.
- The within-replica Metropolis step still skips the draw when $\log\alpha\ge0$. A tie there needs a candidate whose target equals the current one to rounding, which continuous knots do not produce.

### Ensemble shell

`ChainEnsemble.Run`:

1. Validate the arguments, then derive 2C seeds: kernel seeds are children $[0,C)$ and readout seeds are children $[C,2C)$.
2. In parallel, start each chain and run its burn-in.
3. Run sampling rounds of `batchSize` steps. The default is a single round of `samples` steps. Each draw is reduced to a functional vector, whose streaming moments (Welford's update) are accumulated per chain. Full sequences are kept for the leading `EssDim` functionals.
4. If `rHatTarget > 0`, stop once every functional's R̂ meets the target (`RHatDiagnostic.MeetsTarget`). The check runs at round boundaries once each chain holds at least `ChainEnsemble.MinDrawsForStop` (100) draws. An R̂ target with a single chain is rejected.
5. Pool the results:
   - means over $C\cdot n$ draws;
   - population variances over all $C\cdot n$ draws (`FunctionalVariance`), from the per-chain streaming moments merged by Chan's parallel update (`RunningMoments.Merge`, which is associative and commutative);
   - a typed R̂ (`RHatDiagnostic`) per functional;
   - ESS on the leading functionals;
   - acceptance $\sum\text{Accepted}/\sum\text{Attempts}$ from the chains' cumulative counters. These counters include burn-in.

### Convergence diagnostics

**R̂.** The classic Gelman–Rubin statistic, computed from per-chain means and sample variances. It is **not** split and **not** rank-normalized.

$$W=\tfrac1C\sum_i s_i^2,\quad B=\tfrac{n}{C-1}\sum_i(\bar x_i-\bar x)^2,\quad \hat V=\tfrac{n-1}{n}W+\tfrac{B}{n},\quad \hat R=\sqrt{\hat V/W}.$$

Each reading carries a state (`RHatState`):

| State | When | Meets a stopping target |
|---|---|---|
| `Estimated` | $W>0$ | When $\hat R\le$ target |
| `Constant` | $W=0$ and every chain mean is equal | Yes, since there is nothing to converge; but an early stop also needs at least one `Estimated` functional, so an all-constant set runs to the end |
| `Frozen` | $W=0$ and the chain means differ ($\hat R=+\infty$) | Never |
| `Unavailable` | $C<2$, $n<2$, or a non-finite draw | Never |

Streaming moments keep $W$ exactly zero for constant chains; running sums and sums of squares cancel to rounding noise and turn a constant functional into a spurious estimate. `RHatDiagnostic.Worst` summarizes a set in the order constant, estimates by value, frozen, unavailable.

**ESS.** The BDA3 variogram estimator:

$$V_t=\tfrac{1}{C(n-t)}\sum_j\sum_{i\ge t}(x_{j,i}-x_{j,i-t})^2,\qquad \hat\rho_t=1-\frac{V_t}{2\hat V}.$$

- The estimator sums the lag pairs $(\hat\rho_1+\hat\rho_2),(\hat\rho_3+\hat\rho_4),\dots$ and stops at the first *negative* pair. This pairing is offset by one lag from Geyer's $\hat\rho_{2m}+\hat\rho_{2m+1}$, and there is no monotone-sequence step.
- The result is $\tau=1+2\sum$ and $\text{ESS}=\min(Cn/\max(\tau,10^{-6}),\,Cn)$.
- It is NaN for $n<4$ and for a single chain ($C=1$).
- All lags are computed before truncation, so the cost is $O(C n^2)$.

### Inverse transform and spectral noise fields

**`InverseTransform.Sample`** returns $Q(U_i)$ with $U_i\in[0,1)$. `Exponential(scale)` uses $Q(u)=-\text{scale}\cdot\ln(1-u)$.

**`SpectralNoiseField.Generate`** builds a Gaussian random field on a cubic $N^d$ grid, where N must be a power of two:

1. Fill the grid with i.i.d. $N(0,1)$ values.
2. Apply the forward N-D FFT.
3. Scale each frequency bin by $|k|^{-\beta/2}$ (integer frequencies folded about N/2; DC set to 0).
4. Apply the inverse FFT, keep the real part, and standardize to mean 0 and population standard deviation 1.

The power spectrum is $\propto|k|^{-\beta}$ in **grid-index** units: $\beta=0$ is white noise, 1 is pink and 2 is brown.

**`SpectralNoiseField.Sample`** maps the query box onto grid nodes $0..N-1$ and interpolates multilinearly over the $2^d$ corners of the enclosing cell. The field is therefore **not periodic over the box**, even though the grid is periodic.

## Configuration

| Setting | Default | Meaning |
|---|---|---|
| `ChainEnsemble.Run` rHatTarget / batchSize | 0 / 0 | No early stop / one round of `samples` steps |
| `ChainEnsemble.MinDrawsForStop` | 100 | Draws per chain before an R̂ target may stop sampling |
| `GeometricLadder` levels / betaMin | (required) | L ≥ 2; $\beta_{\min}\in(0,1)$ |
| `SpectralNoiseField` presets | White 0 / Pink 1 / Brown 2 | Spectral exponent β |

## Status and caveats

- **ESS of a constant functional depends on rounding.** It comes out as about $C/2$ or as $Cn$ depending on whether the chain means round exactly to the constant. Rank-normalized bulk and tail ESS are planned.
- **Acceptance rates include burn-in.** The counters are never reset.
- **The ESS pairing is non-standard.** Treat ESS as indicative when you compare it with other libraries. It is NaN for single-chain runs.
- **`SpectralNoiseField` is not periodic.** Its documentation describes the grid as periodic, but queries map the box to nodes 0..N−1, so the sampled field is not periodic over the box.

## Source map

| Role | Files (`ThermoMapper/src/maths/samplers/`) |
|---|---|
| Targets | `targets/SamplingTarget.cs` |
| Reversible jump | `rjmcmc/ReversibleJumpChain.cs`, `IRjMove.cs`, `Proposal.cs`, `MoveStat.cs`, `ScheduledMove.cs`, `ParallelTempering.cs`; `IRjTarget.cs` and `TemperedTarget.cs` are unused since step 3.1 and go in step 3.6 |
| Ensemble | `ensemble/ChainEnsemble.cs`, `ChainDiagnostics.cs`, `EnsembleRun.cs`, `IChain.cs`, `IEnsembleModel.cs` |
| i.i.d. and fields | `iid/InverseTransform.cs`; `SpectralNoiseField.cs` (at the `samplers/` root) |

**Tests:** `tests/Maths.Samplers.Tests` (sources in `tests/samplers/`) owns the samplers' own checks, with targets built from `Maths.Samplers` alone: chain contracts, the reversible-jump generality check, the target contract, parallel tempering (each rung's law on both paths, ladder validation) and the diagnostics policy.

**See also:** [BayesianRegression](./BayesianRegression.md) (BARS, BGP) · [ProbabilityAndInformation](./ProbabilityAndInformation.md) · [DesignPrinciples](./DesignPrinciples.md)
