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
  - `TargetConditionalEss` 0.99;
  - `MutationSteps` 5 sweeps per stage;
  - `ResampleBelow` 0.5;
  - multinomial resampling (D05).
- **`Chains` counts independent SMC runs** (replicates), 4 by default, run in parallel within `MaxDegreeOfParallelism`.
  - *Refused with `Smc` set:* `Tempering`, since SMC is its own tempering, and `RHatTarget` > 0, since SMC ends when b reaches 1.
  - `BurnIn`, `Samples` and `Thin` are not read under SMC, and the resolved configuration says so.
- **Weights are kept.** Each replicate's final population, with its weights, is the fit's draws, through `BassFit.DrawWeights` (retain, don't discard).
  - `PointwiseLogDensity` takes the weights inside its log ratios, exactly. Check that against an equal-weight resample.
- **`DrawSample` becomes one continuous weighted quantile.** It is a weighted type 7 that equals type 7 at equal weights at any scale, varies continuously with the weights, and gives the minimum at p = 0 and the maximum at p = 1.
  - Read Akinshin's weighted type 7 (*Weighted quantile estimators*, arXiv 2304.07265, §3's effective sample size and §6) and build it, or a form you show equivalent.
  - Argon's node formula is not it: under weights (0.45, 0.45, 0.1) its last node lands at 1.27, so Q(1) is not the maximum.
  - BARS's unit-weight calls stay bit-identical, through the equal-weights path.
- **The between-run reading.** Across replicates, report R̂ of the monitored quantities over their final populations, resampled to equal weights, as a between-run reading. Keep each replicate's ESS and families beside it, never merged with chain R̂ (the SMC integration section), and with no stop rule.
- **The evidence.**
  - Each replicate's normalizer, with its genealogy error.
  - Stepping-stone over the stage traces, when traced.
  - Pooled across replicates by `LadderEvidence.Pool`'s rule, the larger of the within-run and between-run errors.
  - Reported under D03's proper defaults only, as the ladder's is.
- **Resume:** an SMC fit resumes from its checkpoint (`SmcCheckpoint`), bit-identically, like an MCMC fit.

## Working rules (a shared checkout, with a parallel chip)

- **Work in `D:\aghado01\ThermoMapper` on `main`.** Worktrees fail on this repo's long paths. Switch no branch, push nothing.
- **The C5 chip works beside you** on `src/maths/regression/spline/bass/calibration/` and its tests. Leave its files alone. If a build fails in a file you do not own, wait and retry; never edit it. Calibration's adapters (`BassEmulators.cs`) read `BassFit`'s and `BassBasisFit`'s internals, so keep those shapes, or name any change in your report.
- **Build and test into your own artifacts root:** add `-p:ArtifactsRoot=D:\aghado01\ThermoMapper\artifacts\chip-c4\` to every `dotnet build` and `dotnet test`. Do not run `scripts/fact-harness.ps1`. Delete `artifacts\chip-c4\` when you finish.
- **Your files:**
  - in `src/maths/regression/spline/bass/`: `Bass.cs`, `BassConfig.cs`, `BassResolved.cs`, `BassRun.cs`, `BassDiagnostics.cs`, `BassFit.cs`, `BassPointwise.cs`, `BassEvidence.cs`, and new files for the adapter;
  - `src/maths/samplers/ensemble/Readouts.cs` (`DrawSample`);
  - `src/maths/samplers/smc/`, for accessors only; name each one;
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
- Q(0) = min and Q(1) = max under unequal weights;
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
- Carry the population's diagnostics on `BassDiagnostics`: per replicate, the stage count, each stage's ESS, conditional ESS, families and acceptance; across replicates, the between-run R̂.
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
- the between-run R̂ of M, σ², the log posterior and the fitted values;
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
