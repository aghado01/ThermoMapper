# BASS stage 4: tempering and evidence, S8 (chip D brief)

Date: 2026-10-28. A chip of the [fitting segment's stage-4 plan](port-fitting-brief.md#stage-4-plan-decided-2026-10-28), which the user decided on 2026-10-28.
- *What exists.* S1–S7 and S9 are built, reviewed and on `main`. That includes the run (chip A's [report](stage4-run-brief.md#reports)): `Bass.Fit` and `Resume`, the exact prior sampler on F, the diagnostics, and checkpoints with bit-identical resume.
- *This chip builds S8:* parallel tempering on D05's likelihood path inside the run, the exact reference rung, and the evidence. It is the last engine step; S10's qualification follows in the main session.
- *Not in scope.* SMC is stage 8 of the port, not S8. The config's `Smc` field joins with it.

Read first: [AGENTS.md](../../../AGENTS.md) and [CONTRIBUTING.md](../../../CONTRIBUTING.md).

The decisions, in the [fitting segment](port-fitting-brief.md):
- D05 in full: the path, the two tempered forms, the reference rung, the schedule and ladder, the evidence and its priors record, and the R defects;
- D03's propriety, the evidence bullet and the proper-defaults subsection;
- the contract shapes: `TemperingSpec` in `BassSampling`, and `BassEvidence` with `BassEvidencePriors` in `BassFit`;
- the stage-4 plan's S8 row and choice 2 (τ under tempering through a per-replica block).

D10's adaptation note is in the [numerics segment](port-numerics-brief.md#d10-note-adopted-2026-10-27-a-policy-for-adaptive-kernels): the ladder adapts during burn-in and freezes after it.

The checklist's reports:
- S5, on the target's forms, `Defines(b)` and the slopes;
- S6, on τ's chain block, which a replica runs unchanged;
- chip A's notes for chip D, at the end of its report.

The stage-4 oracle chip's notes for S8 ([its brief](stage4-oracle-modes-brief.md), "For the main session"), on the integrated form near b = 0:
- its slope is singular there;
- log Z_b falls steeply, so rungs need geometric spacing down to about c/N;
- the form is undefined at b = 0 under Jeffreys, and below about 1/N under the flat intercept.

In the code:
- `ParallelTempering` with its per-replica blocks, and `TemperingLadder` with its schedules (DEO) and adaptation;
- `LadderEvidence` (N11) for stepping-stone and path sampling, and `IReferenceSampler`;
- BARS's ladder and evidence as the template: `BarsEnsemble.cs`, `BarsEvidencePriors`;
- BASS: `BassTarget` (forms, `Defines`, slopes), `BassTauBlock`, `BassStructureMoves` (its snapshot and restore, which are ready for replicas), `BassPriorSampler` (an `IReferenceSampler<BassState>`), and `BassRun.cs`, whose `StartChain` is where a chain becomes a ladder.

## Working rules (a shared checkout)

- **Work in `D:\aghado01\ThermoMapper` on `main`.** Worktrees fail on this repo's long paths. Switch no branch, push nothing.
- **Build and test into your own artifacts root:** add `-p:ArtifactsRoot=D:\aghado01\ThermoMapper\artifacts\chip-tempering\` to every `dotnet build` and `dotnet test`. Do not run `scripts/fact-harness.ps1`. Delete `artifacts\chip-tempering\` when you finish.
- **Run oracle tests filtered** to BASS (`--filter "FullyQualifiedName~Bass"`): the TDA oracle tests time out when several builds run at once.
- **Your files.** BASS's files under `src/maths/regression/spline/bass/` as the work needs them: `BassConfig.cs` (`TemperingSpec`), `BassFit.cs` (`Evidence`, between `Diagnostics` and `ResidualDependence`), `BassRun.cs`, `BassDiagnostics.cs`, `BassSweep.cs`, and new files for the ladder and the evidence. Tests go in new `tests/regression/BassTempering*Tests.cs` and `tests/oracle/BassTempering*ParityTests.cs`.
- **Shared code.** A change to `Maths.Samplers` stays generic, is tested in its own host, and is named in your report. The main session owns the fitting segment, its checklist and `docs/Maths`; no other chip is running.
- **Commit your own files only, by explicit path** (`git add <paths>`, never `-A` or `.`), one commit per task, in the repo's message style, ending with `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`.
- **Line endings:** briefs and C# are LF; `docs/Maths/*.md` are CRLF.
- **Decisions are not yours to settle.** Where D05, D10 or the contract shapes leave a default open (a ladder's default levels and hot end, the swap cadence, the evidence's estimator on each form), implement a reversible recommendation and put the question in your report with your reasons. Do not re-open decided choices.
- **Report** under "Reports" at the end of this brief, newest first: what landed, the commits, the tests and suite counts, the review passes (below) and anything the main session must know.

## Task 1: the tempered run

**The config.** `TemperingSpec` in `BassSampling` (null is untempered, the default, D05), a values-only closed set, serialized and fingerprinted like the rest. It holds:
- the path: likelihood by default, powered selectable;
- the tempered form: integrated by default, linear selectable;
- the ladder: its levels and hot end, or an adaptive ladder;
- the swap schedule: DEO;
- the exact reference rung, opt-in.

`Validate` and `BassResolved` refuse what cannot run:
- a rung where `BassTarget.Defines(b)` is false on the chosen form;
- the reference rung on a form or prior without b = 0, or on the powered path.

**The ladder.**
- One `ParallelTempering` per chain over the replicas. Each replica runs its own `BassStructureMoves` and `BassTauBlock` (the per-replica blocks) at its temperature.
- Swaps use the likelihood path's ℓ_b alone, since the prior cancels (D05).
- The ladder adapts during burn-in and is frozen after it (D10).
- The exact reference rung, when asked for, draws from `BassPriorSampler`.

**Retention and the rest of the run.**
- Draws are retained from the b = 1 replica only, composed as now.
- The diagnostics gain the tempering readout, following BARS's: swap acceptance per pair, round trips, the frozen ladder.
- Checkpoints carry the replicas' states, streams, the ladder and its adaptation state. Resume stays bit-identical.

## Task 2: the evidence

- **`BassEvidence`.** Stepping-stone, the integrated form's estimator. Path sampling runs beside it where the slope comes free: on the linear form it is ℓ₁. The integrated form's slope also exists (S5), but its b = 0 end is singular under a sampled τ, so state how path sampling treats it, or why it is left out. Each estimate carries its Monte Carlo error, pooled over chains as BARS pools.
- **`BassEvidencePriors`** (D05's point 3) records what the evidence depends on:
  - the coefficient prior and τ's shape;
  - the count prior;
  - the term measure and the interaction caps;
  - the support rule;
  - the centre and its resolved level;
  - σ²'s prior as specified and as resolved.
- **Null under the improper options:** the flat intercept or Jeffreys σ². Null too when `ReportEvidence` is false or the run is untempered.

## Checks

- *The endpoints.*
  - A b = 1 replica samples the posterior: the oracle's fixture A, as `BassFitParityTests` checks it.
  - A b = 0 replica samples the prior on F: the enumerated F of `BassSweepTests.EnumerateF`, with τ's prior.
  - On both forms.
- *Swap detailed balance* on an enumerable structure space, against the enumeration oracle's `swap` mode (F-O6). R's swap density is asserted corrected, a defect fixture, not an option.
- *The ladder's evidence* against the oracle's exact log Z on fixture B, under the proper defaults, on both forms. `BassTargetParityTests.Rungs_MatchTheOraclesLadder` already gives each rung's exact log Z_b and slope mean to compare against.
- *Bit-identity and independence.* Resume is bit-identical with tempering on, and results do not depend on the parallelism bound.
- *Refusals.* Undefined rungs are refused, and a reference rung where b = 0 is undefined is refused.

## Review before reporting

A task is not done when its tests pass. Before the report, review the work as a separate pass, as the main session reviewed S5 and S6 and the earlier chips reviewed their own (every one of those reviews found real defects or untested paths):

1. **Reread the whole change** as a reviewer seeing it for the first time: `git diff <the commit before your first>..HEAD -- <your paths>`. For each change, ask:
   - Is it correct? Check swap ratios (the prior's cancellation, ℓ_b on each form), the reference rung's law, the estimators and their errors, edge cases (a two-rung ladder, b_min near the form's limit, a replica refusing a swap-in state), the exception paths, the parameter names in exceptions, thread safety, and the contracts it implements.
   - Does it read like the surrounding code?
   - Do its docstrings say what it does?
2. **Hunt the coverage gaps.** Look for each option and convention the code accepts that no test runs:
   - both forms and both paths;
   - the reference rung on and off;
   - an adaptive or fixed ladder;
   - the intercept conventions and σ² priors (evidence null under the improper ones);
   - sampled or fixed g;
   - scalar, categorical and direct functional data.

   Look also for each refusal path no test reaches. Prefer independent references (the R oracle's modes, enumeration) to self-consistency.
3. **Fix what you find,** add the tests and rerun. Commit the fixes separately (`fix(bass): …` or `test(bass): …`).
4. **Review again** if the pass found anything substantive. Stop when a pass finds nothing.
5. **Record it** in your report: what each pass found, what was fixed, and anything left, with the reason.

## Reports
