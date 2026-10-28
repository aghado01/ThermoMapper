# BASS stage 4: the run, S7 (chip A brief)

Date: 2026-10-28. A chip of the [fitting segment's stage-4 plan](port-fitting-brief.md#stage-4-plan-decided-2026-10-28), which the user decided on 2026-10-28. S1–S6 are built and reviewed (the [checklist's reports](port-fitting-checklist.md#reports)). This chip builds S7, the run: `Bass.Fit` and `Bass.Resume` on the chain ensemble, the exact prior sampler for starts, the diagnostics and checkpoints. Two other chips run beside it: chip B builds S9's pointwise density, chip C writes the docs section. The main session supervises, reviews your commits and wires B's density into your `BassFit` once both land. Chip D (tempering, S8) follows you and builds on your run.

Read first: [AGENTS.md](../../../AGENTS.md) and [CONTRIBUTING.md](../../../CONTRIBUTING.md) (zero-dependency C#, net10.0, lean docstrings, no compatibility shims, tests assert corrected references). Then, in the [fitting segment](port-fitting-brief.md): D02, D04, D05's exact reference rung, "Posterior ownership and prediction", "Starting on the support", the contract shapes (entry points, `BassFit`, `BassDraw`, the diagnostics bullet) and the stage-4 plan's S7 row and "Starts". In the checklist, S6's report and its "For S7" note.

## Working rules (a shared checkout)

- **Work in `D:\aghado01\ThermoMapper` on `main`.** Worktrees fail on this repo's long paths. Switch no branch, push nothing.
- **Build and test into your own artifacts root:** add `-p:ArtifactsRoot=D:\aghado01\ThermoMapper\artifacts\chip-run\` to every `dotnet build` and `dotnet test`. Do not run `scripts/fact-harness.ps1`, which builds the default root. Delete `artifacts\chip-run\` when you finish.
- **Run oracle tests filtered** to BASS (`--filter "FullyQualifiedName~Bass"`): the TDA oracle tests time out when several builds run at once.
- **Your files:** `src/maths/regression/spline/bass/Bass.cs`; new files for the run, the fit, the diagnostics, the run options and the prior sampler; `BassSweep.cs`, `BassState.cs` and `BassStructurePrior.cs` where the run needs hooks (per-order support tallies, a term sampler, snapshot state); new tests `tests/regression/BassRun*Tests.cs`, `tests/regression/BassPriorSampler*Tests.cs` and `tests/oracle/BassFitParityTests.cs`; and `tests/regression/BassSweepTests.cs` only to share its `EnumerateF`. Chip B owns new `BassPointwise*`/`BassResidualDependence*` files, chip C owns `docs/Maths/BayesianRegression.md`, the main session owns the fitting segment and its checklist. A change to shared code (`Maths.Samplers`, `Maths.Regression.GPrior`) stays generic, is tested in its own host, and is named in your report.
- **Commit your own files only, by explicit path** (`git add <paths>`, never `-A` or `.`), one commit per task, in the repo's message style, ending with `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`.
- **Line endings:** briefs and C# are LF; `docs/Maths/*.md` are CRLF. Check `git diff --stat` after a scripted edit.
- **Decisions are not yours to settle.** Where the decided docs leave a statistical default or a contract shape open, implement the recommended option only if it is reversible and needs no new contract field, and put the question in your report with your recommendation. Do not re-open decided choices (D01–D06, the contract shapes, the stage-4 choices).
- **Report** under "Reports" at the end of this brief, newest first: what landed, the commits, the tests and suite counts, the review passes (below) and anything the main session must know.

## Task 1: exact draws from the prior restricted to F

The sampler the starts use (the stage-4 plan's "Starts") and S8's exact reference rung (D05).
- **The law.** The globally truncated prior of D02: the unconstrained measure over collections, Γ(h₁ + M)/(h₂ + 1)^M · Πₘ π(termₘ) with the M! orderings, restricted to F. Draw λ ~ Gamma(h₁, rate h₂) and M ~ Poisson(λ), which gives the count law with its orderings; then M terms independently from D01's term measure π, including the mixed order rule decided 2026-10-28 (each part's order uniform on {0, …, cap}, all-empty outside F). In the g mode, τ from its prior. Accept the collection only when it lies in F: distinct terms, the cap, support, D06's gates, and under Jeffreys no exact fit (`BassTarget.Refusal`). Rejection from the globally truncated measure is exact.
- **What it needs.** A sampler for one term from π beside `BassStructurePrior.LogPrior`: the same law, drawn rather than scored.
- **Its account.** Report the acceptance rate; refuse with a clear message after a declared attempt limit.
- **Checks.** Draws against the prior on F enumerated exactly: S6's constant-likelihood fixture (`BassSweepTests.EnumerateF`, 2,263 states; make it reusable rather than copying it). Check the count law, the order law and τ's prior, and the categorical part on a small mixed fixture.

## Task 2: the run

- **Entry points** as the contract shapes give them: `Bass.Fit(config, data, run?)` and `Bass.Resume(checkpoint, config, data, run?)`, with `BassRunOptions` carrying what is not the model (a parallelism bound, cancellation, a checkpoint request, progress). Neither the run options nor parallelism changes a result.
- **The ensemble.** One `BassChain` per chain on `ChainEnsemble`'s shell (`IEnsembleModel`), with BARS's run as the template (`BarsEnsemble.cs`, `Bars.cs`).
  - Each chain gets three streams from the seed tree: the moves, τ's block and composition. The trajectory then does not depend on how often draws are retained (S6's note).
  - Starts come from Task 1.
  - Burn-in, samples and thin follow `BassSampling`.
  - Each chain's draws go into a `DrawStore<BassDraw>`.
- **Term tables.** Each chain interns into its own `BassTermTable`. The fit's table merges them in chain order, with each chain's ids remapped, so ids never depend on scheduling.
- **`BassFit`** per the contract shapes. `Evidence` and `ResidualDependence` are null until S8 and S9 land; `DrawWeights` is null.

## Task 3: diagnostics (`BassDiagnostics`)

- Rank-normalized R̂ and ESS (`ChainDiagnostics`) for M, τ (none under fixed g), σ², and fitted values at monitored points.
- **The monitored points are an open choice.** The config has no field for them. Implement a deterministic default that needs none, such as a few training cases spread through the data's order; for a direct fit, at a middle output location. Propose a declared list (`BassSampling`) in your report.
- Acceptance per move type, from `MoveStats`.
- Refusals by reason, the chains' `Tally<BassRefusal>` merged.
- Empty births.
- The support-rejection rate per order (D01's tilt): tally support refusals and attempts by the candidate's orders per part.
- The posterior mass at the cap (D02).
- The R̂ stop (`RHatTarget`, `BatchSize`, `MinEssPerChain`), as the shell provides it.

## Task 4: checkpoints and resume

- **The restart state** holds each chain's `BassState` (terms and u), its three streams, its counters, move statistics, tallies and term table. The design reloads from the state; S6's sync loads a foreign structure whole.
- **The checkpoint** is `EnsembleCheckpoint` with `RunIdentity`. It rides beside the posterior, never inside it.
- **Resume** refuses a fingerprint that differs from the config's and data's (`Bass.Fingerprint`).

**Checks for Tasks 2–4** (the plan's):
- a resumed run is bit-identical to an uninterrupted one;
- results do not depend on the parallelism bound;
- every start lies in F;
- an end-to-end `Bass.Fit` on the oracle's fixture A matches its exact posterior count law and τ, as `BassChainParityTests` does for the chain;
- the diagnostics read right on a fit you can reason about: mass at the cap 0 when the cap cannot bind, the refusal tallies equal to the chains' sum.

## Review before reporting

A task is not done when its tests pass. Before the report, review the work as a separate pass, as the main session reviewed S5 and S6 (each of those reviews found real defects or untested paths):

1. **Reread the whole change** as a reviewer seeing it for the first time: `git diff <the commit before your first>..HEAD -- <your paths>`. For each change, ask:
   - Is it correct? Check formulas, ratios, signs, edge cases, the off-support and exception paths, the parameter names in exceptions, thread safety and the contract it implements.
   - Does it read like the surrounding code?
   - Do its docstrings say what it does?
2. **Hunt the coverage gaps.** Look for each option and convention the code accepts that no test runs:
   - the intercept conventions (level centring, flat, zero);
   - the σ² priors (scaled, Jeffreys, a declared IG);
   - sampled or fixed g, and both move schedules;
   - scalar, categorical and direct functional data.

   Look also for each refusal or failure path no test reaches. Prefer independent references (the R oracle, enumeration, dense computation) to self-consistency.
3. **Fix what you find,** add the tests and rerun. Commit the fixes separately (`fix(bass): …` or `test(bass): …`).
4. **Review again** if the pass found anything substantive. Stop when a pass finds nothing.
5. **Record it** in your report: what each pass found, what was fixed, and anything left, with the reason.

## Reports
