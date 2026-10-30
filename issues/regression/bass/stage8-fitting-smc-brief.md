# Stage 8, C4: BASS fitting by SMC (chip brief)

Date: 2026-10-29. Builds the calibration checklist's C4, the fitting SMC adapter on D05's sequence. D05 (in [port-fitting-brief.md](port-fitting-brief.md)), decided 2026-10-27, fixes the sequence:
- from the prior, b₀ = 0 < … < b_S = 1, on the likelihood path and the integrated form;
- each particle's incremental weight ℓ_{b′} − ℓ_b;
- D04's sweep at the stage's b as the mutation;
- stages placed by conditional ESS (Zhou, Johansen and Aston 2016);
- multinomial resampling, with independent offspring streams;
- the normalizer as the evidence under D03's proper defaults.

The generic SMC core (N8) is built. This chip writes BASS's adapter, and answers D10's stage 8 qualification question. The choices this brief makes beyond D05 are in "Decided for stage 8" below.
- *Not in scope:*
  - cut-calibration SMC, the parallel C5 chip's ([its brief](stage8-cut-smc-brief.md));
  - any change of default: MCMC stays the default, and SMC is an option;
  - `docs/Maths` and the checklists, which the main session updates after verifying.

Read first: [AGENTS.md](../../../AGENTS.md) and [CONTRIBUTING.md](../../../CONTRIBUTING.md). Then:
- D05 and D04 in [port-fitting-brief.md](port-fitting-brief.md), and its contract-shapes section (`BassSampling.Smc`, a `SmcSpec?`, "joins with stage 8");
- the calibration segment's "SMC integration" ([port-calibration-brief.md](port-calibration-brief.md));
- D10 in [port-numerics-brief.md](port-numerics-brief.md), its stage 8 bullet included;
- the numerics checklist's reports on N8 and on the normalizer's genealogy error (Lee and Whiteley 2018).

In the code:
- *The SMC core,* `src/maths/samplers/smc/`: `SmcStage`, `ISmcParticles`, `ISmcKernel`, `SmcRun`, `SmcSampler.RunTempered` and `NextInverseTemperature`, `SmcOptions`, `SmcStageReport`, `SmcStageTrace`, `SmcCheckpoint`, and `SmcWeights`, which holds the genealogy error.
- *BASS's pieces for particles:*
  - `BassPriorSampler`, the exact prior draws on F;
  - `BassTarget`'s tempered forms;
  - `BassStructureMoves` and `BassTauBlock`, a sweep at b;
  - `BassLadder`, whose hand-over copies a design and shares its columns, the pattern a particle's fork follows;
  - `BassEvidence` and `LadderEvidence`.
- *The fit:* `Bass.cs` (`Fit`, `Resume`), `BassRun.cs`, `BassFit.cs` (`DrawWeights`, already threaded through `Predict` and the other views), `BassPointwise.cs`, `BassDiagnostics.cs`.
- *The weighted quantile:* `DrawSample` in `src/maths/samplers/ensemble/Readouts.cs`, and the stage 9 docs' note on its jump between equal and near-equal weights.

## Decided for stage 8

These are the main session's recommendations, taken by the user on starting this chip.
- **`SmcSpec`** on `BassSampling.Smc`, null meaning MCMC. Its values are provisional per-engine defaults, read in the qualification and not tuned:
  - `Particles` 1,000;
  - `Replicates` 4;
  - `TargetConditionalEss` 0.99;
  - `MutationSteps` 5 sweeps per stage;
  - `ResampleBelow` 0.5;
  - multinomial resampling (D05).
- **`Replicates` counts independent SMC runs,** run in parallel within `MaxDegreeOfParallelism`. The name is true of SMC; "chains" is not (names follow methods).
  - *Refused with `Smc` set:* `Tempering`, since SMC is its own tempering, and `RHatTarget` > 0, since SMC ends when b reaches 1.
  - `Chains`, `BurnIn`, `Samples` and `Thin` are not read under SMC, and the resolved configuration says so.
- **Weights are kept.** Each replicate's final population, with its weights, is the fit's draws, through `BassFit.DrawWeights` (retain, don't discard).
  - `PointwiseLogDensity` takes the weights inside its log ratios, exactly. Check that against an equal-weight resample.
  - *Particles share ancestry,* so errors computed as though the draws were independent read optimistic: `PointwiseLogDensity`'s Monte Carlo errors, and any Kish ESS. Report the families beside them, and never present them as calibrated. The replicates carry the honest error.
- **`DrawSample` becomes one continuous weighted quantile:** Akinshin's weighted type 7 (*Weighted quantile estimators*, arXiv 2304.07265, §6). It averages the weighted quantile function over the window [(h − 1)/n*, h/n*], with h = (n* − 1)p + 1 and n* Kish's effective size.
  - It equals type 7 exactly at equal weights at any scale.
  - It is continuous in the weights, including as a weight goes to 0, so a negligible particle moves no quantile.
  - It does not anchor Q(1) at the maximum under unequal weights. No rule can do both: the maximum over the positive weights jumps when the top atom's weight reaches 0, so anchoring and continuity are incompatible. Continuity wins, since SMC's populations carry many negligible weights.
  - *The rejected alternative,* argon's note: today's midpoint rule with p mapped onto the span of the midpoints. It is exact at equal weights and anchored at both ends, but a negligible weight acts as a full node.
    - Atoms at 0 and 1 (weight ½ each) plus one at 100 with weight 10⁻⁹: Q(0.99) reads 97.03, against 0.99 without it.
    - Atoms at 0 and 10 with one at 1 of weight 10⁻⁹: Q(0.4) reads 0.80, against 4.00.
    - Akinshin's form gives 0.99 and 4.00 (the main session's check, 2026-10-29).
    - Today's midpoint rule has the same flaw. It never bit, since every weight so far is equal.
  - BARS's unit-weight calls stay bit-identical, through the equal-weights path.
- **The between-run reading.** For each monitored quantity, the spread of the replicates' posterior means, in posterior-sd units: SD over replicates of the weighted mean, over the pooled posterior sd. Split R̂ is meaningless here, because a population has no order to split, and an unsplit R̂ would invite comparison with D10's 1.01, which assumes autocorrelated chains. Keep each replicate's ESS and families beside it, with no stop rule.
- **The evidence.**
  - Each replicate's normalizer, with its genealogy error.
  - Stepping-stone over the stage traces, when traced.
  - Pooled across replicates by `LadderEvidence.Pool`'s rule, the larger of the within-run and between-run errors.
  - Reported under D03's proper defaults only, as the ladder's is.
- **Resume:** an SMC fit resumes from its checkpoint (`SmcCheckpoint`), bit-identically, if that fits `BassFit.Checkpoint` without reworking its type. Otherwise defer it and say so. It is the most droppable item, since an SMC run reaches b = 1 in one pass.

## Working rules (a shared checkout, with a parallel chip)

- **Work in `D:\aghado01\ThermoMapper` on `main`.** Worktrees fail on this repo's long paths. Switch no branch, push nothing.
- **The C5 chip works beside you** on `src/maths/regression/spline/bass/calibration/` and its tests. Leave its files alone. If a build fails in a file you do not own, wait and retry; never edit it. Calibration's adapters (`BassEmulators.cs`) read `BassFit`'s and `BassBasisFit`'s internals, so keep those shapes, or name any change in your report.
- **Build and test into your own artifacts root:** add `-p:ArtifactsRoot=D:\aghado01\ThermoMapper\artifacts\chip-c4\` to every `dotnet build` and `dotnet test`. Do not run `scripts/fact-harness.ps1`. Delete `artifacts\chip-c4\` when you finish.
- **Your files:**
  - in `src/maths/regression/spline/bass/`: `Bass.cs`, `BassConfig.cs`, `BassResolved.cs`, `BassRun.cs`, `BassDiagnostics.cs`, `BassFit.cs`, `BassPointwise.cs`, `BassEvidence.cs`, and new files for the adapter;
  - `src/maths/samplers/ensemble/Readouts.cs` (`DrawSample`);
  - `src/maths/samplers/smc/` stays as it is: its public API (`SmcRun.Start`, `Stage`, `ConditionalEss`, `Advance`, `SmcSampler.RunTempered`, `TemperToPosterior`, `NextInverseTemperature`, `SmcWeights`) is enough. Name any further need in the report rather than edit it;
  - their tests and new ones;
  - this brief's Reports section.
- **Commit your own files only, by explicit path,** one commit per step, in the repo's message style, ending with `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`.
  - The user's working-copy edits (`issues/regression/argon-bass-qa.md`, `bass-questions.md`, `issues/spc/opus-spin-glass-threads.md`) stay unstaged.
  - Before each commit, read `git diff --cached` and confirm every hunk is yours.
- **Line endings:** briefs and C# are LF.
- **Off means bit-identical:** with `Smc` null, every MCMC fit's draws and the eight `BassBitIdentityTests` hashes hold.

## The steps

Each step lands with its checks and ends with a fresh-eyes review of its own code and tests, for correctness and for coverage gaps, fixed before the next step. Name the passes and their fixes in the report.

**Step 1. The weighted quantile** (`DrawSample`). Build it as decided above, and check:
- type 7 exactly at equal weights at three scales;
- continuity as one weight moves from equal;
- every value within [min, max] under unequal weights, and continuity as a weight goes to 0, the negligible-atom cases above included;
- tied values pooled, zero weights dropped, and independence from insertion order;
- BARS's peak summaries bit-identical.

**Step 2. The particles and the kernel.**
- *A particle* is a `BassState` with its own design and moves.
- *`Draw`* is `BassPriorSampler`'s exact draw on F, with τ from its prior.
- *`Fork`* copies the design as the ladder's hand-over does: the factor copied, the columns shared.
- *The kernel* is `MutationSteps` of D04's sweep at the stage's b, the structure move then the τ block, with no adaptation (D10's policy).
- *Checks:*
  - each stage's law against enumeration on fixture A, at three b, through the kernel alone;
  - a fork is independent of its parent;
  - parallel moves give the same result as serial ones.

**Step 3. The adapter and the fit.** `Bass.Fit` with `Smc` set runs `Chains` replicates of `SmcSampler.RunTempered` from the prior to b = 1, with the stage seeds from the seed tree.
- Compose each final particle as a retained draw is composed (σ², the coefficients, λ), with its weight.
- Carry the population's diagnostics on `BassDiagnostics`: per replicate, the stage count, each stage's ESS, conditional ESS, families and acceptance; across replicates, the spread of their posterior means.
- Evidence as decided, with `BassEvidencePriors`; resume from the checkpoint.
- *Checks:*
  - on the oracle's enumerable fixtures, with proper priors: the normalizer against the exact ln Z on fixtures B and C (the ladder evidence's references: −28.589 and −30.903), within its genealogy error; and on fixture A, the posterior count law, inclusions and E[τ] against the exact posterior (`BassChainParityTests`'s references);
  - the stepping-stone estimate over the stages against the same;
  - weights summing to 1, and the views (`Predict`, `Sensitivity`, `TrainingFit`) taking them;
  - `PointwiseLogDensity` with weights against its equal-weight resample;
  - resume bit for bit;
  - every refusal in "Decided for stage 8".

**Step 4. The qualification** (D10's stage 8 question). This is exploration: it reads, and changes no default. Behind `BASS_QUALIFICATION` = 1, on S10's Friedman, mixed and functional fixtures, at N = 250, 1,000 and 4,000, with four replicates and S10's master seeds. Per fixture and N:
- the stage count and time;
- the final population's mean log posterior per replicate, against S10's best MCMC chain's median (Friedman −892.4) and its worse chains. The question is whether the population finds the better representation and keeps it;
- the replicates' spread of the posterior means of M, σ², the log posterior and the fitted values, in posterior-sd units;
- the families left;
- the evidence and its pooled error, against the tempered ladder's evidence on the same fixture where S10 or the pass recorded one;
- σ's coverage on Friedman.

Report what the readings show, kept apart from the numbers. Recommend no default.

## Review before reporting

Take a fresh-eyes pass over the whole change and fix what it finds:
- *Correctness:*
  - the incremental weight is ℓ_{b′} − ℓ_b on the integrated form;
  - the kernel leaves each stage invariant;
  - forks are independent;
  - the normalizer's error is the genealogy one, under multinomial resampling;
  - weights reach every view;
  - with `Smc` null, everything is bit-identical.
- *Coverage:*
  - fixed and sampled g;
  - a direct functional fit and a mixed one;
  - D03's improper options (the normalizer is then no evidence, and the run starts from the prior over (structure, τ) on the linear form, per D05);
  - one replicate;
  - a run that hits `MaxStages`.

## Report

Under "Reports" below, newest first:
- what landed, with the commits;
- the configuration and the API;
- each check with its tolerance;
- the qualification's tables, then what they show;
- the suites with their counts (Regression, Samplers, the BASS oracle classes);
- the review passes and their fixes;
- anything the main session must know.

## Reports

Newest first.

### C4's report: BASS fitting by SMC — 2026-10-29

**What landed** (on `main`, by explicit path, nothing pushed):
- `d0a0fd8` step 1: `DrawSample.Quantile` is Akinshin's weighted type 7.
- `9e8c349` step 2: the particles and the kernel (`BassSmc.cs`: `BassParticle`, `BassParticles`, `BassParticleTarget`, `BassSmcKernel`, `BassSmcAccount`).
- `b0978d5` step 3: the adapter and the fit (`BassSmcRun.cs`'s `BassSmcFit`; `SmcSpec` on `BassSampling.Smc`; `BassSmcReadout` on `BassDiagnostics.Smc`; `BassEvidence.OfSmc`; weights through `PointwiseLogDensity.Weights`; the stage budget and resume).
- `73352f0` step 4: the qualification, `BassSmcQualificationTests`, behind `BASS_QUALIFICATION` = 1.

**The configuration and the API.**
- `SmcSpec` (in `BassConfig.cs`), null meaning MCMC: `Particles` 1,000, `Replicates` 4, `TargetConditionalEss` 0.99, `MutationSteps` 5, `ResampleBelow` 0.5, `Resampling` multinomial, as decided; and three the build needed:
  - `Form`, integrated by default, the likelihood path's tempered form, one shape with `TemperingSpec.Form`. Under D03's improper options the integrated form has no b = 0, so `BassResolved` refuses it there with the reason "the linear form reaches it (D05)", as it refuses a tempering ladder's reference rung.
  - `MaxStages`, 1,000, each replicate's stage budget. A replicate that spends it before b = 1 returns flagged (`BassSmcReadout.ReachedPosterior` false), with no evidence or residual check, and `BassFit`'s posterior views (`Predict`, `PointwiseLogDensity`, `Sensitivity`, `Effects`, `TrainingFit`, `Summary`, `SensitivitySummary`, `ResidualDependenceFor`) refuse it; `Metadata` and `Traces` read it. Its checkpoint continues it. This is what "a run that hits `MaxStages`" does: a reported state, as D10's verdict is, not an exception.
  - `TraceStages`, off, for stepping-stone over the stages (and path sampling on the linear form) beside the normalizer.
- *Refused with `Smc` set* (`BassConfig.Validate`): `Tempering`, an `RHatTarget` above 0, `Particles` < 2, `Replicates` < 1, `TargetConditionalEss` outside (0, 1), negative `MutationSteps`, `ResampleBelow` outside [0, 1], `MaxStages` < 1. `Chains`, `BurnIn`, `Samples`, `Thin`, `BatchSize`, `MinEssPerChain` and `StoppingScope` are logged in `Resolved.AutoPicks` as "not read".
- *The run* (`BassSmcFit`): the master seed derives one seed per replicate (`SeedTree`), each fanned into the run's seed and the composition's. Replicates run at once within `BassRunOptions.MaxDegreeOfParallelism`, their particles sharing what is left. Each replicate is `SmcRun.Start` at b = 0, then `NextInverseTemperature` and `Advance` until b = 1 or the budget: `RunTempered`'s own loop, with the budget a return where `TemperToPosterior` throws. `src/maths/samplers/smc/` is untouched.
- *A particle* is a `BassState` with its own `BassStructureMoves`. `Draw` is `BassPriorSampler`'s exact draw on F (a fresh sampler per draw, so draws run concurrently), loaded into the particle's design; `Fork` copies the design as the ladder's hand-over does (`DesignFor`: the factor and fit copied, the columns shared) and gives the fork its own term array.
- *The kernel* (`BassSmcKernel`): each move is one D04 sweep at the stage's b, the structure move on D01's schedule then τ's slice step, through a reversible-jump chain built on SMC's stream over the particle's own moves; no adaptation. After each sweep the particle's tallies, its accepted candidate committed first, merge into the kernel's account under a lock and restart, so resampling's dropped particles lose none; the counts are integer sums, independent of order.
- *The fit*: `Draws` holds one store per replicate, its final population composed as a retained draw is (σ², coefficients, λ; D04); `DrawWeights` pools each particle's normalized weight over R, summing to 1; `Identity.Chains` is R. `Diagnostics.Smc` (`BassSmcReadout`): per replicate (`BassSmcReplicate`) its seed, stage reports (ESS, conditional ESS, resampling, distinct ancestors, families, acceptance, normalizer increments, genealogy variance), final b, ln Ẑ with its genealogy error, final Kish ESS and families, prior attempts, kernel account and mass at the cap; and `Spread`, per monitored quantity (M, τ, the log posterior, σ², each fitted point), the replicates' weighted means, the pooled posterior's mean and sd, and the SD of the replicates' means over that sd. The chain readings are unavailable under SMC (each `ConvergenceDiagnostic` `Unavailable`, the verdict `Undetermined` in every group, not gating, `Chains` empty), and `ModeSeparation` reads the replicates' weighted means.
- *The evidence* (`BassEvidence.OfSmc`), under the proper defaults when every replicate reached b = 1: each replicate's normalizer with its genealogy error as a stepping-stone `LadderEstimate` (with rungs, by `LadderEvidence.FromSmc`, when traced), pooled by `LadderEvidence.Pool`; path sampling joins on the linear form when traced. The ladder's flags read NaN and false.
- *Resume*: it fits `BassFit.Checkpoint` without reworking `EnsembleCheckpoint`. Each replicate's `SmcCheckpoint`, kernel account and prior attempts ride the engine slot (`BassSmcState : IChainState`), with the moments and trace slots empty. `Bass.Resume` retargets the checkpoint to the new model (each particle re-evaluated, which gives the saved evaluations) and runs up to `MaxStages` more stages, with `SmcOptions.MaxStages` = stages done + budget. The stage streams are derived by position, so the result is the uninterrupted run with the larger budget, bit for bit. A complete run's checkpoint resumes to the same fit; an MCMC checkpoint is refused.
- *Shared code outside the BASS files* (named, as the brief asks):
  - `src/maths/regression/Pointwise.cs` (CRLF kept): `PointwiseLogDensity.Weights`, an init-only property, validated (one per draw, finite, non-negative, positive sum, no chain lengths beside them).
  - `src/maths/regression/PredictiveChecks.cs`: WAIC's lppd as log Σ W p and p_WAIC as Σ W (ℓ − ℓ̄)² / (1 − Σ W²), which is the sample variance at equal weights; PSIS-LOO's log ratios log W − ℓ, its lppd weighted, zero-weight draws dropped. Null weights take the old code path, unchanged.
  - `BassModel.Create` reads `SmcSpec.Form` after `TemperingSpec.Form` (one line).
  - `BassRun.cs`: the monitored points' evaluation moved into `BassPointInputs`, shared with SMC's spread, same arithmetic; `BassEnsemble.Remap` internal.
  - `tests/regression/BassMovePassTests.cs`: S10's Friedman data hash re-pinned, `e38c2c6c05ff8dbc` → `a2946a493a466d54`, because the default config's JSON gained `"Smc":null`; the data unchanged (as N14's re-pin before it).
  - `tests/regression/BassPresentationTests.cs`: the independent weighted-quantile reference is now Akinshin's (interval overlap over unpooled draws), not the midpoint rule.

**The checks, with their tolerances.**
- *Step 1* (`ReadoutTests`, 4 new or rewritten): type 7 exactly (bit-equal) at equal weights at scales 1/n, 0.2 and 3.7; continuity as one weight moves from equal, |ΔQ| ≤ n·|ε|·range over ε from ±10⁻² to 10⁻¹¹, for the least, an inner and the largest draw; every quantile within [min, max] of the weighted draws and nondecreasing in p over 50 random weighted samples; as a weight goes to 0 (10⁻³ … 10⁻³⁰⁰) every quantile approaches the one without it within 10·w·range; the brief's cases, Q(0.99) = 0.99 and Q(0.4) = 4.00 to 10⁻⁶; ties pooled (0.8/11 by hand), zero weights dropped, insertion-order and merge independence (the existing test). BARS's unit weights take the untouched `Type7` path; its peak summaries are bit-identical by construction and its tests pass.
- *Step 2* (`BassSmcTests`, `BassSmcParityTests`): the kernel alone, run as a chain at fixed b from exact prior draws, against the enumeration oracle's tempered law on fixture A at b = 0.05, 0.3 and 1, on both forms and under fixed g: the count law within 5 batch-means SE + 0.004, E[τ] within 5 SE, every term's inclusion within 5 SE + 0.005 (4 × 30,000 sweeps per b). A fork moves exactly as its parent would from the same stream, and its 200 sweeps leave the parent's design (terms, ‖r‖², Q, β̂) and later path untouched, bit for bit. A tempered run of 96 particles with 6-way parallel moves gives the serial run's normalizer, weights, stage reports, particles and account exactly. The kernel refuses another model's stage.
- *Step 3*:
  - The normalizer against the oracle's exact ln Z (2,000 particles, 4 replicates, defaults otherwise): each replicate within 4 of its genealogy errors + 0.01, and the pooled within 4 of its error + 0.01. Fixture B, integrated: −28.597 ± 0.027 against −28.589 (replicates −28.586, −28.612, −28.588, −28.603, each ± 0.05); B, linear: −28.593 ± 0.032; C, categorical: −30.888 ± 0.019 against −30.903.
  - Stepping-stone over the traced stages equals each normalizer to 10⁻⁹, one rung per stage; path sampling on the linear form −28.592 ± 0.032, diagnostic 0.002.
  - Fixture A's exact posterior (8 replicates of 1,000; 5 SE of the replicates' spread + 0.003 on the count law, + 0.005 on inclusions, + 0 on E[τ]): defaults at M_max = 3, M = 1–3 at 0.176, 0.750, 0.074 against 0.176, 0.751, 0.073, E[τ] 0.015526 against 0.015476 (± 0.00014), all 48 inclusions; fixed g at M_max = 2, M = 1, 2 at 0.287, 0.713 against 0.280, 0.720.
  - The views take the weights: weights sum to 1 (10⁻¹²), 1/R per replicate; `TrainingFit`'s mean against Σ w f by hand (10⁻⁹); `Summary`'s σ² mean against Σ w σ² (10⁻¹²); `Sensitivity`'s decomposition carries `DrawWeights`; `Predict` runs under every option.
  - `PointwiseLogDensity` with weights: equal weights at scales 1/S and 3 give the unweighted WAIC and PSIS-LOO to 10⁻⁹; integer weights give a repeated sample's WAIC lppd exactly (10⁻⁹); against an equal-weight multinomial resample of 40,000 draws, WAIC −24.692 against −24.698 and PSIS-LOO −24.719 against −24.736 (allowance 0.1).
  - Resume: a fit with `MaxStages` 3 is flagged and refused by the posterior views, and resumes, three stages a time, to the uninterrupted run's draws (terms, τ, σ², λ, coefficients), weights, stage reports, normalizers, prior attempts and accounts, exactly; a complete run's checkpoint resumes to the same fit; another config's checkpoint and an MCMC one are refused.
  - Every refusal above, the config's JSON round trip, and the fingerprint moving with `Smc`; cancellation stops a run with `OperationCanceledException`, and progress reports the sweeps against the budget.
  - Options, each reaching b = 1 with weights summing to 1: fixed g (no τ in the spread, every τ the held one), a direct functional fit, a mixed one, Jeffreys σ² and the flat intercept on the linear form (no evidence), the linear form traced (path sampling present), one replicate (spreads and mode separation NaN, the evidence's error the within-run one), systematic resampling (genealogy error and the evidence's error NaN), and `MutationSteps` 0 (no sweeps, acceptance NaN).
- *Off means bit-identical*: the eight `BassBitIdentityTests` hashes hold.

**The qualification** (`BassSmcQualificationTests`, `BASS_QUALIFICATION` = 1, Release; SMC's defaults, four replicates of 1,000, S10's master seeds 53, 57 and 61). Each fixture ran at N ≈ 250, 1,000 and 4,000, and at S10's own size, where F-M2's Phase 0 (M4) recorded the MCMC chains. The functional fixture keeps its 20 locations, so its N is n_x · 20 (13, 50, 60 and 200 cases). The log posterior is M4's: the log marginal plus the log structure prior, τ's density left out. There is one row per run; replicates are in seed order. Times are wall-clock while C5's tests shared the machine, so they are rough.

| fixture | N | stages | s | mean log post per replicate | spread of the means (nats) | families | Kish ESS |
|---|---|---|---|---|---|---|---|
| Friedman | 250 | 157–167 | 71 | −545.2, −549.9, −527.4, −526.8 | 23.1 | 2, 2, 1, 1 | 625–999 |
| Friedman | 500 (S10) | 179–198 | 81 | −904.0, −914.3, −917.5, −914.1 | 13.5 | 1, 1, 1, 2 | 776–938 |
| Friedman | 1,000 | 203–228 | 198 | −1664.0, −1692.8, −1659.0, −1668.9 | 33.8 | 1, 1, 1, 1 | 703–996 |
| Friedman | 4,000 | 298–309 | 294 | −6096.6, −6031.0, −6079.6, −6076.6 | 65.6 | 1, 1, 1, 1 | 599–1000 |
| mixed | 250 | 122–128 | 29 | −271.1, −268.7, −269.5, −270.4 | 2.4 | 3, 4, 6, 4 | 729–1000 |
| mixed | 300 (S10) | 126–128 | 50 | −312.7, −313.1, −311.8, −312.1 | 1.3 | 9, 5, 3, 2 | 809–997 |
| mixed | 1,000 | 157–162 | 93 | −853.2, −854.6, −854.8, −850.1 | 4.7 | 1, 4, 3, 5 | 605–978 |
| mixed | 4,000 | 178–201 | 159 | −3134.2, −3124.2, −3132.2, −3124.2 | 10.0 | 1, 2, 3, 2 | 539–960 |
| functional | 260 | 125–128 | 29 | −23.3, −25.4, −22.1, −21.6 | 3.8 | 11, 9, 14, 10 | 862–1000 |
| functional | 1,000 | 155–168 | 100 | 135.2, 133.4, 133.2, 137.5 | 4.3 | 7, 3, 3, 6 | 526–982 |
| functional | 1,200 (S10) | 152–168 | 40 | 176.9, 174.5, 176.1, 176.4 | 2.4 | 6, 3, 3, 4 | 537–1000 |
| functional | 4,000 | 180–205 | 141 | 754.3, 767.2, 760.5, 755.0 | 12.9 | 2, 1, 1, 2 | 833–972 |

*Against M4's chains at S10's sizes* (medians of the log posterior):

| fixture | SMC replicates' weighted medians | M4's chains | what each replicate holds most |
|---|---|---|---|
| Friedman | −897.0, −910.7, −914.5, −907.3 | −892.4, −916.0, −926.2, −933.9 | mostly 3 × x1·x2, M 8.7–10.3; replicate 0's third structure is the best chain's 2 × x1·x2 (1.6%, −893.7) |
| mixed | −311.6, −311.8, −310.8, −311.3 | −309.9, −309.9, −312.5, −331.9 | x1 + x2·z + 3 × x3·z + z in all four, the best chains' structure; nothing near chain 2's −331.9 |
| functional | 177.6, 175.4, 177.0, 176.9 | 181.1, 168.6, 152.7, 151.7 | 2 × x0 + 2 × x0·t + x1·t (M = 5) in all four, chain 0's; nothing near chains 1 and 3's M = 7 at ≈152 |

*The replicates' spread of the posterior means, in posterior-sd units* (M, σ², the log posterior, the largest over the fitted points):

| fixture | N ≈ 250 | S10's size | 1,000 | 4,000 |
|---|---|---|---|---|
| Friedman | 0.64, 0.39, 0.66, 0.56 | 0.63, 0.40, 0.41, 0.82 | 0.56, 0.14, 0.74, 0.71 | 0.95, 0.17, 0.98, 0.85 |
| mixed | 0.24, 0.09, 0.20, 0.58 | 0.29, 0.04, 0.11, 0.17 | 0.37, 0.11, 0.32, 0.59 | 0.51, 0.19, 0.46, 0.58 |
| functional | 0.15, 0.10, 0.38, 0.23 | 0.12, 0.08, 0.27, 0.19 | 0.13, 0.10, 0.51, 0.27 | 0.71, 0.05, 0.82, 0.70 |

*The evidence* (the pooled normalizer ± its error, within-run and between-run errors; each replicate's genealogy error 0.42–0.83), and *σ's 95% interval*:

| fixture | N | evidence | within, between | σ's interval (truth) |
|---|---|---|---|---|
| Friedman | 250 | −502.56 ± 4.94 | 0.40, 4.94 | [1.020, 1.262] (1): missed |
| Friedman | 500 | −875.76 ± 4.21 | 0.42, 4.21 | [1.003, 1.151]: missed |
| Friedman | 1,000 | −1634.30 ± 1.41 | 0.42, 1.41 | [1.016, 1.116]: missed |
| Friedman | 4,000 | −5981.97 ± 8.72 | 0.42, 8.72 | [0.997, 1.043]: covered |
| mixed | 250 | −253.52 ± 1.85 | 0.37, 1.85 | [0.464, 0.557] (0.5): covered |
| mixed | 300 | −296.67 ± 1.09 | 0.36, 1.09 | [0.474, 0.559]: covered |
| mixed | 1,000 | −828.43 ± 0.39 | 0.39, 0.38 | [0.488, 0.534]: covered |
| mixed | 4,000 | −3096.94 ± 2.51 | 0.40, 2.51 | [0.5008, 0.5236]: missed |
| functional | 260 | −12.86 ± 0.49 | 0.23, 0.49 | [0.187, 0.223] (0.2): covered |
| functional | 1,000 | 139.24 ± 0.50 | 0.33, 0.50 | [0.187, 0.204]: covered |
| functional | 1,200 | 180.97 ± 1.57 | 0.38, 1.57 | [0.188, 0.203]: covered |
| functional | 4,000 | 768.96 ± 1.24 | 0.41, 1.24 | [0.1911, 0.1995]: missed |

Neither S10 nor F-M2's pass recorded a tempered ladder's evidence on these fixtures (S10's tempered runs have no reference rung), so there is none to set beside these.

**What the readings show.**
- *Mixed and functional: yes, at the representation's level.* At S10's sizes every replicate holds the structure the best MCMC chains hold, where one mixed chain and two functional chains sat 22–29 nats below. The replicates agree within 1.3–2.4 nats. Within that representation the replicates' medians sit 1–2 nats (mixed) and 3–6 nats (functional) below the best chain's, a gap in knots, not in the terms' variables.
- *Friedman: no.* At S10's size the replicates span the same range as the chains: −897 to −915 against the best chain's −892.4 and the worse chains' −916 to −934. Their populations mostly hold three x1·x2 terms where the best chain holds two. One replicate comes within 5 nats; it carries the best chain's structure at 1.6%. The replicates disagree by 13–66 nats, more at larger N.
- *Lineages collapse.* At the end every Friedman replicate descends from one or two initial particles, and from one at N ≥ 1,000. The other fixtures thin with N, the functional from 9–14 families at N = 260 to 1–2 at 4,000. So a final population is one lineage's local exploration. The comparison between modes that resampling makes happens among the few lineages alive at low b, not near b = 1, where the modes' mass ratios (e⁻¹² and below) would decide it. D10's caveat named the particle count as the variable to read; what binds here is the lineages, which these settings do not keep.
- *The genealogy error is uninformative here.* Once a replicate's families reach one, Lee and Whiteley's estimate reads a relative variance of 1, so its error is √ln 2 = 0.83 whatever the run; most replicates above read 0.83 or near it. The between-run error is 1–21× the within-run one, largest on Friedman, so the pooled error rests on the replicates' spread. With four replicates that spread is itself rough. On the enumerable fixtures, with 650–1,050 families surviving, the genealogy error was calibrated.
- *Prediction agrees better than structure,* as in D10's reading of F-M2's pass: σ²'s replicate spread is 0.04–0.4 posterior sd, while M's and the log posterior's reach 0.6–1.0 on Friedman and at N = 4,000 on every fixture.
- *σ on Friedman* is missed at N = 250, 500 and 1,000, every interval lying above the truth (lower ends 1.003–1.020), and covered at 4,000. That is the direction of untempered MCMC's registered failure. The mixed and functional intervals cover except at N = 4,000, where both miss narrowly (0.5008 and 0.1995).
- *Cost:* stages grow slowly with N, about 125 at N ≈ 250 to about 300 on Friedman at 4,000; a fit took 0.5–5 minutes here.
- No default is recommended.

**The suites.** Samplers 206/206. Regression 1169/1169, the full suite after the last commit. An earlier full run read 1170/1171; its one failure was the data hash above, since re-pinned, and C5's own test changes moved the count between the runs. The BASS oracle classes (`Maths.Oracle.Tests.Bass*`, which run here with R and the BASS checkout) 71/71, the 8 new `BassSmcParityTests` cases included.

**Review passes and their fixes.**
- *Step 1:* the clamp to [min, max] guards the window sum's rounding; tests for continuity at the equal ↔ unequal boundary and for monotonicity in p were added.
- *Step 2:* fixed g was missing from the kernel's oracle law; added. The account merges each sweep's tallies after committing the accepted candidate, so no acceptance is lost to `ResetTallies`.
- *Step 3:* the weighted WAIC's variance takes 1 − Σ W² so equal weights reproduce the S − 1 form; weighted draws refuse chain lengths, since a population has no order for r_eff; the resample check's allowance was tightened from 0.5 to 0.1 nats after reading the differences (≤ 0.02). Coverage gaps closed: systematic resampling, zero mutation steps, cancellation and progress.
- *The whole change*, against the brief's list: ℓ_{b′} − ℓ_b is `Temperature.LogDensity` on the likelihood path, the prior cancelling in the increment; the kernel's invariance and the forks' independence are checked above; the genealogy error is the multinomial one, NaN otherwise; every view takes the weights; `Smc` null is bit-identical. The fixed and sampled g, the functional and mixed fits, the improper options, one replicate and the stage budget are each covered.

**What the main session should know.**
- *The core's one further need:* `SmcSampler.TemperToPosterior` throws at `MaxStages`, so the adapter runs its loop itself. A variant that returns at the budget would let the adapter call the core directly.
- *Presentation:* `BassSummary`'s run lengths show `BurnIn`, `Samples` and `Thin` for an SMC fit, which reads none of them, and `BassTraces` reads each replicate's population in particle order, which is no trace. Both are stage 9's files; I left them.
- *Docs* (the main session's): `DrawSample`'s rule changed, so the stage 9 docs' note on the jump between equal and near-equal weights is obsolete; `PointwiseLogDensity.Weights` and the weighted checks; the SMC fit.
- *Cost:* a kernel move costs what a chain sweep costs (83 against 70 µs at N = 27, 221 against 216 µs on Friedman at N = 500, both under load); the schedule's bisection, 60 conditional ESS evaluations over the population, about 0.1 s a stage at 1,000 particles, runs single-threaded in the core.
- *C5's adapters:* `BassFit`'s and `BassBasisFit`'s shapes are unchanged; `BassDiagnostics` gained a trailing optional `Smc`.
- *`artifacts\chip-c4\`* is deleted.
