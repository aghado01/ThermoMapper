# N14: the convergence verdict (chip brief)

Date: 2026-10-29. Builds D10's stopping semantics, decided by the user on 2026-10-29: [port-numerics-brief.md](port-numerics-brief.md), section "D10 (decided 2026-10-29)". Read that section first; this brief fixes the API and the order of work.
- *What it builds:*
  - a convergence verdict read at the end of every `ChainEnsemble` run;
  - engine-declared groups of functionals, and a declared stopping scope;
  - BASS's three groups, with each draw's log posterior among BASS's monitored quantities;
  - the mode-separation readings for BASS;
  - the verdict on BARS's diagnostics.
- *Not in scope:*
  - any change to a chain's draws;
  - a stop that reads the mode-separation readings;
  - `docs/Maths`, the checklists and D10's text, which the main session updates after verifying.

Read first: [AGENTS.md](../../../AGENTS.md) and [CONTRIBUTING.md](../../../CONTRIBUTING.md). Then:
- D10's section;
- the fitting checklist's report "F-M2 Phase 2 verified: no arm clears", for why the verdict matters: untempered chains that never agree;
- [local-move-pass-brief.md](local-move-pass-brief.md)'s report, whose harness computed each draw's log posterior and checked it against the chain's own density.

In the code:
- *The shared ensemble:* `src/maths/samplers/ensemble/` (`ChainEnsemble.Run` and its private `MayStop`, `EnsembleRun`, `IEnsembleModel<TDraw>`, `ChainDiagnostics`, and the `ConvergenceDiagnostic`, `RHatDiagnostic`, `EssDiagnostic` and `DiagnosticState` types).
- *BASS:* `BassRun.cs` (its functional reduction: `FunctionalDim`, `MonitoredDim` and the order of M, τ, σ² and the fitted values), `BassDiagnostics.cs`, and `BassConfig.cs` (`BassSampling`).
- *BARS:* `Bars.cs` and `BarsEnsemble.cs`, and how its diagnostics are presented under B01's display policy.
- *The other clients,* `BgpSampler.cs`, `Baps.cs` and `HierarchicalPSpline.cs`, which get the verdict through `EnsembleRun`.

## The API, fixed here

Stage 7's chip is building calibration in parallel. Its nested reference sampler (C3) will read these names, so build them as written. Shapes not named here are yours; state them in the report.

- **`ConvergenceOutcome`**, an enum in `Maths.Samplers.Ensemble`, beside `ChainEnsemble` with three values:
  - `Converged`;
  - `NotConverged`;
  - `Undetermined`, when nothing could be estimated: one chain, or every gating functional constant or frozen.
- **`ConvergenceVerdict`**, a record:
  - `Outcome`;
  - `Target`: the R̂ target read, which is the run's target when set, otherwise `ChainEnsemble.ReferenceRHatTarget`;
  - `MinEssPerChain`;
  - `Groups`, one `GroupVerdict` per declared group.
- **`GroupVerdict`**, a record:
  - `Name`;
  - `Outcome`;
  - `GatesStop`;
  - the failing functionals, each with its index, its name, its R̂ and its ESS.
- **`ChainEnsemble.ReferenceRHatTarget`** = 1.01 (Vehtari et al. 2021), a constant beside `MinEssPerChain`.
- **`EnsembleRun.Verdict`**, of type `ConvergenceVerdict`, on every run.
- **The model declares its groups:**
  - `IEnsembleModel<TDraw>.FunctionalNames` defaults to "f0", "f1", ….
  - `IEnsembleModel<TDraw>.Groups` defaults to one group named "all" holding every functional, so a model that declares nothing gets today's rule.
  - A group lists functional indices, and every functional belongs to exactly one group.
- **The stopping scope:** `ChainEnsemble.Run` takes the names of the groups that gate an R̂ stop, with null meaning every group.
  - The verdict reports every group, whatever the scope.
  - The run's outcome is the conjunction over the gating groups.
  - An unknown name is refused.

## The rules

They are the stop's rules (B01), applied to each group at the end of every run, and also between rounds when a target is set:
- *A monitored functional* passes when its rank-normalized R̂ meets the target and its bulk and tail ESS reach `MinEssPerChain` per chain.
- *An unmonitored functional* passes on its classic R̂.
- *A frozen or unavailable reading fails.*
- *A constant functional* passes but cannot certify on its own. A group with no estimate is `Undetermined`.
- *Combining outcomes:* any failing gating group makes the run `NotConverged`; otherwise, if every gating group is `Undetermined`, the run is `Undetermined`; otherwise it is `Converged`.
- *The stop between rounds* reads the same rule over the gating groups, so a run that stops early always ends `Converged`. Keep one implementation that both the stop and the verdict call.

## BASS

- **The log posterior joins the monitored quantities.** It is each retained draw's collapsed log π(structure, u | y) at b = 1, the chain's own density at its state; the pass's harness checked that value to 10⁻¹¹. Under tempering it is the b = 1 replica's.
  - Place it so the functional order stays readable, and update `MonitoredDim`.
  - Name every functional: "M", "tau", "sigma2", "log_posterior", and the fitted values by their monitored points.
- **The groups:**
  - "predictive": σ² and the fitted values at the monitored points;
  - "level": the log posterior;
  - "structural": M, and τ when it is sampled.
- **`BassSampling.StoppingScope`**, an enum `BassStoppingScope { All, Predictive }`, defaulting to `All`, validated, and echoed in the resolved configuration and its JSON. `Predictive` gates the stop on "predictive" alone.
- **`BassDiagnostics`** carries the run's `ConvergenceVerdict` (with `StoppedEarly` and `SamplesUsed` as they are) and the mode-separation readings:
  - the gap between the chains' mean log posteriors, the largest minus the smallest;
  - the largest Euclidean distance between two chains' mean fitted values at the monitored points.
  - The verdict reports them, and nothing stops on them.
- **Unchanged:** the draws, and every fit without an R̂ target. A fit with a target may stop later, since the log posterior now enters the gate; report which existing tests change, and why.

## BARS and the other clients

- *BARS* presents the verdict beside its R̂ readouts, under B01's display policy, and declares its groups if they read better than "all". State the choice.
- *BGP, BAPS and HierarchicalPSpline:* where a result already carries the ensemble's diagnostics, it carries the verdict too. Otherwise, list the client in the report.

## Working rules (a shared checkout, with a parallel chip)

- **Work in `D:\aghado01\ThermoMapper` on `main`.** Worktrees fail on this repo's long paths. Switch no branch, push nothing.
- **The stage 7 chip works beside you.** Its files are new ones under `src/maths/regression/spline/bass/calibration/`, its tests, and two R oracles. It leaves yours alone, and you leave its alone.
  - If a build fails in a file you do not own, wait and retry; never edit it.
  - `BassFit.cs` and `output/BassBasisFit.cs` belong to neither chip. If you need a change there, name it in the report instead.
- **Build and test into your own artifacts root:** add `-p:ArtifactsRoot=D:\aghado01\ThermoMapper\artifacts\chip-n14\` to every `dotnet build` and `dotnet test`. Do not run `scripts/fact-harness.ps1`. Delete `artifacts\chip-n14\` when you finish.
- **Your files:**
  - `src/maths/samplers/ensemble/`, and a new file there for the verdict's types;
  - in BASS, `BassRun.cs`, `BassDiagnostics.cs`, `BassConfig.cs` and `BassResolved.cs`;
  - BARS's `Bars.cs` and `BarsEnsemble.cs`;
  - the other clients' result records, where the verdict surfaces;
  - their tests, and new ones;
  - this brief's Reports section.
- **Commit your own files only, by explicit path,** one commit per step or a coherent pair, in the repo's message style, ending with `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`.
  - The user's working-copy edits to `issues/regression/argon-bass-qa.md` and `bass-questions.md` stay unstaged.
  - Before each commit, read `git diff --cached` and confirm every hunk is yours.
- **Line endings:** briefs and C# are LF.

## Order

Each step ends with a fresh-eyes review of its own code and tests, for correctness and for coverage gaps, fixed before the next step. Name the passes and their fixes in the report.
1. **The bit-identity baseline:** confirm `BassBitIdentityTests`' eight hashes hold before you start, and recheck them after every commit.
2. **The shared verdict:** the types, the groups and the scope, with one rule used by both the stop and the verdict. Commit this first, so stage 7's C3 can build on it.
3. **BASS:** the log posterior, the names and groups, `StoppingScope`, the verdict and the mode-separation readings on `BassDiagnostics`.
4. **BARS, and the other clients.**

## Checks

- *Synthetic ensembles* (`tests/samplers`):
  - chains at two levels in one group fail that group alone, and name its failing functionals;
  - constant and frozen functionals follow B01's rules;
  - one chain gives `Undetermined`;
  - an unknown scope name is refused;
  - a scope excluding a failing group can stop, while the verdict still reports that group.
- *The rule agrees with itself:* the stop's decision on a round equals the verdict's on the same traces, across randomized ensembles.
- *A run that stops early ends `Converged`.* A run whose target is not met by the end of its budget ends `NotConverged`, a case today's `StoppedEarly` cannot express.
- *BASS:*
  - the log posterior against the chain's own density at three draws per chain, untempered and tempered;
  - S10's untempered Friedman, mixed and functional fits at their seeds, where F-M2's pass read the R̂ values each group should report;
  - one of them tempered;
  - `StoppingScope.Predictive` against `All` on a fit with a target.
- *Bit-identity:* the draws, and the fits without a target, unchanged.
- *Suites:* Regression, Samplers, and the BASS oracle classes, all passing.

## Report

Under "Reports" below, newest first:
- what landed, with the commits;
- the API as built, against the fixed names above;
- the groups BARS declared, and the other clients' coverage;
- the verdicts read on S10's fits, untempered and tempered;
- the tests that changed, and why;
- the suites with their counts;
- the review passes and their fixes;
- anything the main session or the stage 7 chip must know.

## Reports

Newest first.

### N14 report (2026-10-29): the verdict landed, shared, on BASS, BARS and the other clients

*What landed,* on `main`, nothing pushed:
- `e8547d3` the shared verdict (step 2), committed first for stage 7's C3;
- `e88b33e` BASS: the log posterior, the three groups, `StoppingScope`, the verdict and the mode-separation readings (step 3);
- `1c26dcd` BARS's two groups and the verdict on BARS, BGP, BAPS and `HierarchicalPSpline` (step 4).

`BassBitIdentityTests`' eight hashes held at `27b7cf3` before the work and after every commit. No draw changed.

*The API as built,* against the fixed names (all in `Maths.Samplers.Ensemble`):
- **As fixed:**
  - `ConvergenceOutcome { Converged, NotConverged, Undetermined }`;
  - `ConvergenceVerdict(Outcome, Target, MinEssPerChain, Groups)`;
  - `GroupVerdict(Name, Outcome, GatesStop, Failing)`;
  - `ChainEnsemble.ReferenceRHatTarget = 1.01`, beside `MinEssPerChain`;
  - `EnsembleRun.Verdict`;
  - `IEnsembleModel<TDraw>.FunctionalNames` ("f0", "f1", … by default) and `.Groups` (one group "all" by default).
- **Shapes of my own:**
  - `FunctionalReading(Index, Name, RHat, Ess)`, one per failing functional. `Ess` is null for an unmonitored functional, which has classic R̂ only.
  - `FunctionalGroup(Name, Functionals)`, with `FunctionalGroup.AllName` and `FunctionalGroup.All(n)`. A model with no functionals declares no group, and its verdict is `Undetermined`.
  - `ConvergenceVerdict.Group(name)`, a lookup that throws `KeyNotFoundException`.
  - `Run(…, ParallelOptions? parallel = null, IReadOnlyCollection<string>? stoppingScope = null)`. The scope is the last parameter, so no positional caller broke; pass it by name. Null takes every group. An empty scope or an unknown name is refused (`ArgumentException`, parameter "stoppingScope"). A target needs a monitored functional in a gating group, the old check widened to the scope.
  - The declarations are checked at entry: one name per functional, and distinct, non-blank group names. Each group is non-empty, and each functional sits in exactly one group.
- **One rule:** `ConvergenceVerdict.Of`, read by the stop between rounds and by the verdict at the end, on the same readings (classic R̂ for every functional, the rank-normalized diagnostics for the monitored ones). Under the default scope, the conjunction over groups is the old `MayStop` exactly, which is deleted. So a model that declares nothing stops as before.
- **Frozen fails; the brief disagrees with itself here.** The API bullet gives `Undetermined` for "every gating functional constant *or frozen*". The rules section and D10 say a frozen reading fails. I built the rules: frozen and unavailable readings fail (`NotConverged`), and only one chain, or every gating functional constant, gives `Undetermined`. One chain is undetermined even though split R̂ exists for it: there is no between-chain reading.

*BASS:*
- **The log posterior:** `Temperature.LogPosterior(Target.Evaluate(state))` at b = 1, a ladder's b = 1 replica's state. It equals the chain's own density (`Moves.CurrentLogDensity` of the chain replayed as the run takes it) **bit for bit** at three draws per chain, untempered and tempered (`LogPosterior_IsTheChainsOwnDensity`). It costs one `Evaluate` per retained draw.
- **The order and names:** M, τ (g mode), log posterior, σ², the fitted values; named "M", "tau", "log_posterior", "sigma2", "fitted[i]" or "fitted[i,r]". So the groups run contiguously: structural, then level, then predictive. `MonitoredDim` is 4 + points (3 + points under fixed g).
- **The groups:** `BassDiagnostics.PredictiveGroup` = "predictive", `LevelGroup` = "level" and `StructuralGroup` = "structural", declared in that order.
- **`BassStoppingScope { All, Predictive }`:** `BassSampling.StoppingScope`, `All` by default. It is validated ("Sampling.StoppingScope"), in the config's JSON ("StoppingScope":"All"), and echoed resolved as `BassResolved.StoppingGroups`, the group names the run passes to the ensemble.
- **`BassDiagnostics`** gains `LogPosterior` (its diagnostic), `Verdict` and `ModeSeparation` = `BassModeSeparation(LogPosteriorGap, FittedDistance)`, read from the chains' traces. Both readings are NaN with one chain, and nothing stops on them. `StoppedEarly` and `SamplesUsed` are as they were.
- **`Bass.cs` needed one line,** passing `stoppingScope: model.Resolved.StoppingGroups`. The brief's file list does not name it, and neither chip owns it.

*BARS and the other clients:*
- **BARS declares two groups:** "predictive" (the fitted grid values, `fit[g]`, on classic R̂) and "structural" (`knots`), as D10 groups BASS's. They read better than "all": a knot count the curve does not identify shows as structural disagreement beside predictive agreement.
  - With both gating, their conjunction is the old rule, so BARS's stop is unchanged (`Bars_Tempered` keeps its hash).
  - `BarsResult.Verdict` sits beside `Knots` and `RHatFitMax`. Its readings keep their `DiagnosticState`, B01's display policy.
  - No stopping scope is exposed on `BarsConfig`.
- **BGP, BAPS and `HierarchicalPSpline`:** all three results already carried the ensemble's diagnostics, so all three now carry `Verdict` (`BgpResult`, `BapsResult`, `HierarchicalResult`). Each keeps the one group "all", with its functionals named:
  - BGP: "bandwidth", "noise_ratio";
  - BAPS: "log_lambda", "log_sigma2", "log_lambda[j]";
  - `HierarchicalPSpline`: "log_sigma2", "log_tau_u2", "log_tau_02", "log_tau_02[c]", "population[j]".

  No client is left without the verdict.

*The verdicts on S10's fits* (F-M2's A0 at k = 0 untempered, and A6 for Friedman tempered; `S10_Fits_FailEveryGroup`, gated). Every quantity's R̂ equals the pass's CSV value. Every fit is `NotConverged`, with every group `NotConverged`:

| Fit | Predictive failing (R̂, ESS) | Level | Structural | Gap (nats) | Fitted distance |
|---|---|---|---|---|---|
| Friedman | σ² 1.083 (30), f1 1.747 (6), f2 1.317 (10), f3 1.962 (5) | 2.050 (5) | M 2.348 (5), τ 1.013 (417) | 42.8 | 1.282 |
| mixed | f1 1.299 (17), f2 1.012 (388), f3 1.087 (29); σ² passes | 1.646 (7) | M 1.537 (6) | 21.8 | 0.164 |
| functional | σ² 1.023 (153), f1 1.134 (20), f2 1.497 (8), f3 1.140 (20) | 2.393 (5) | M 2.548 (5), τ 1.042 (61) | 30.9 | 0.050 |
| Friedman, tempered (A6) | σ² 1.077 (33), f1 1.235 (12), f2 1.189 (15), f3 1.049 (55) | 1.264 (11) | M 1.490 (7), τ 1.011 (2005) | 27.0 | 0.552 |

The mixed fit's σ² passes (R̂ 1.007, ESS 1,136) inside a failing predictive group, as the pass read it. Tempered τ fails on R̂ alone, 1.011 against 1.01.

*`StoppingScope.Predictive` against `All`,* on a target 1.05 fit (`BassSweepTests.Data(large)`, 3 chains, `MinEssPerChain` 50):
- predictive stops at 300 draws per chain, against 700 for every group;
- its verdict is `Converged`, and it still reports the non-gating groups: level fails (log_posterior R̂ 1.005, ESS 119), and structural fails (M R̂ 1.011, ESS 65);
- the predictive run's draws are the all-groups run's first ones.

*The tests that changed, and why:*
- **`BassOutputComponentsTests.RHatStoppedComponents_PairOverTheShortest`:** the log posterior in the gate kept its first component from stopping.
  - Under the old gate the components stopped at 800, 2000 and 2000; now all three run to 2000. The test still passed, but trivially: pairing over the shortest was no longer exercised.
  - It now sets `StoppingScope.Predictive` (700, 2000, 2000) and asserts the lengths differ.
- **`BassMovePassTests.Friedman_AtN500_IsS10s`:** its pin moved from `6b52ce789967ee0c` to `e38c2c6c05ff8dbc`, a deliberate regeneration.
  - The data hash reads `Bass.Fingerprint(new BassConfig(), data)`, which includes the config's JSON, now longer by `,"StoppingScope":"All"`.
  - With those 22 characters stripped, a probe reproduced `6b52ce789967ee0c` exactly, so the data are unchanged. The constant's comment says so.
- **Unchanged in outcome:**
  - `BassRunDiagnosticsTests.RHatStop_EndsTheRunEarly` stops at 300 draws under both gates.
  - `BassTemperingTests`' "early stop" stops at 150 of 200 under both.
  - No existing `ChainEnsemble` test changed.
- **Helpers:** `BassRunTests.AssertSameDiagnostics` and the resume helpers of BARS, BGP, BAPS and `HierarchicalPSpline` now compare the verdicts too, through `VerdictAssert.Same`. So every resume test checks that a resumed run's verdict, log posterior and mode separation are the uninterrupted run's.

*New tests:*
- `tests/samplers/ConvergenceVerdictTests.cs` (10). Its rule-agreement test runs 24 randomized ensembles (AR(1) noise, chain offsets, random scopes and targets). Each is replayed for j rounds in one, unchecked. Its verdict is `Converged` at the stop's round and not before; the trials held both early stops and exhausted budgets.
- `tests/regression/BassVerdictTests.cs` (7, the S10 test gated);
- `tests/regression/ClientVerdictTests.cs` (5);
- `tests/regression/VerdictAssert.cs`, the shared checks.

*The suites:*
- Samplers 204/204, after step 2 and again on the final commit;
- Regression 1,068/1,068 with steps 3 and 4 in place, then the bit-identity and verdict classes 20/20 on the final commit;
- the BASS oracle classes (`Maths.Oracle.Tests`, `~.Bass`): 54/54 after step 3, and 59/59 on the final commit, whose filter also picks up stage 7's new calibration parity classes;
- the gated S10 verdict test 1/1 (`BASS_QUALIFICATION=1`).

*The review passes and their fixes:*
- **Step 2:**
  - the rule against the old `MayStop`: identical under one group;
  - the stop's and the verdict's readings come from one helper (`Readings`);
  - the brief's frozen discrepancy, resolved by the rules;
  - in the tests, a frozen-and-unavailable case meant to read an unmonitored functional had it monitored (fixed), and a local function shadowed the test's chain type (renamed).
- **Step 3:**
  - a coverage gap: nothing checked that a resumed fit keeps its verdict and log posterior. Fixed in the shared resume helper.
  - probing the existing target-bearing tests under the old and new gates found the components test gone trivial. Fixed as above.
  - an over-long initializer line in `BassResolved` was split.
- **Step 4:**
  - `VerdictAssert` moved to its own file, so each commit's tests stand alone;
  - every client's resume helper now compares verdicts.

*For the main session and the stage 7 chip:*
- **Every BASS fingerprint changed:** the config's JSON now holds `StoppingScope`. A checkpoint taken before `e88b33e` will not resume, and any other pin that hashes `Bass.Fingerprint` moves the same way; I found only the one above.
- **For C3:**
  - the verdict reads 1.01 when a run sets no target;
  - one chain is always `Undetermined`;
  - pass `stoppingScope:` by name;
  - the API summary is in the shared notes folder (`n14-verdict-api-handover.md`).
- **Left for the main session, as the brief says:** `docs/Maths`, the checklists and D10's text.
- **Not built:** a stopping scope for BARS (`BarsConfig`), and any stop on the mode-separation readings.
