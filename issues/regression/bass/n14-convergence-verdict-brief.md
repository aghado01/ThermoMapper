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
