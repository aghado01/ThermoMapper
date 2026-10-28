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

### S7: the run — 2026-10-28 (chip A)

All four tasks landed, each with its tests, then three review passes. Commits:
- `6609985` the prior sampler (Task 1);
- `b5bca5c` the run (Task 2);
- `87b23bc` the diagnostics (Task 3);
- `759136d` checkpoints and resume (Task 4);
- `b14c84d` and `7ac76ae` pass 1's fixes and gaps, `a402b14` pass 2's fix.

**Task 1: exact draws from the prior on F** (`BassPriorSampler.cs`; `BassStructurePrior.DrawTerm`).
- *The law.* M from the count law with its orderings, p(M) ∝ Γ(h₁ + M)/(M!(h₂ + 1)^M), drawn by inversion and truncated at M_max. That is the law of λ ~ Gamma(h₁, rate h₂), M | λ ~ Poisson(λ), which the brief describes, cut where F's cap rule would refuse anyway, so no Poisson sampler is needed and the sampler never refuses for the cap.
- *The rest of the draw.* Then M terms from D01's term measure (each part's order uniform, on {0, …, cap} with several parts), and τ from its prior in the g mode. A collection is kept only in F: no intercept term, distinct terms, the design's support and D06's gates, and no exact fit under Jeffreys.
- *Its seat.* It implements `IReferenceSampler<BassState>` (the likelihood path), so S8's reference rung can take it as it is.
- *Its account.* Attempts, acceptance, refusals by reason and intercept draws. A draw refuses after a declared limit (10⁶ by default), and its message reports that draw's attempts.
- *Checks.*
  - On S6's fixture, against `BassSweepTests.EnumerateF` (now shared, 2,263 states), 150,000 draws:
    - the count law 0.6021, 0.3072, 0.0772, 0.0135 against 0.6037, 0.3049, 0.0781, 0.0132;
    - the order law;
    - χ² 357 on 339 over every state;
    - E[ln τ], E[τ], and τ apart from M;
    - the acceptance rate 0.6540 against F's exact share of the proposal, 0.6550 ± 0.0010.
  - A general term-space enumerator in the tests (`BassTermSpace`) reproduces the same F, and checks:
    - a mixed fixture with a three-level factor (1,593 states, χ² 688 and 745 on 732, the joint (J, L) order law, the subset sizes);
    - a direct functional fixture (6,785 states, χ² 335 on 305);
    - Jeffreys on exact-fit data, which excludes 342 exact fits.
  - The term draws against the measure over every term of three spaces: both subset measures on a four-level factor, where they differ (χ² 697 and 750 on 734), and a functional fit (148 on 164).

**Task 2: the run** (`Bass.Fit`, `BassRun.cs`, `BassFit.cs`, `BassRunOptions.cs`).
- *Chains.* One `BassChain` per chain on `ChainEnsemble` (`BassEnsemble : IEnsembleModel<BassState>`). Burn-in counts sweeps, and a retained draw follows every `Thin` sweeps.
- *Retention.* Each draw is composed into the chain's `DrawStore<BassDraw>`. Each chain interns into its own term table, and the fit merges them in chain order, remapping each draw's ids and re-sorting its coefficients with them.
- *Streams.* The chain seed fans out to three streams: the start, the moves and τ's block. The shell's readout seed drives composition. That makes four streams where the brief named three; the start's own stream keeps the trajectory independent of how many attempts the start took.
- *`BassRunOptions`.* The parallelism bound, cancellation (checked each sweep), the checkpoint request, and progress (sweeps over every chain, reported from the chains' threads).
- *Checks.*
  - Each chain rebuilt by hand from the seed tree reproduces the fit's draws exactly, through the merge (thin 1 and 2).
  - The fit is identical under parallelism bounds 1, 2, 3 and none, and with or without the checkpoint and progress options.
  - A thinned run's states are every k-th state of the unthinned one.
  - Thirteen options run end to end, each draw in F: level centring, flat and zero intercepts; scaled, Jeffreys and declared-IG σ²; sampled, unit and declared fixed g; both schedules; categorical data under both subset measures; direct functional data; column mass; degree 2.
  - The options are checked; cancellation stops a run, before it and from a chain's thread; progress counts every sweep.
  - `BassFitParityTests` against the oracle's exact posterior, four chains of 40,000 draws:
    - fixture A at M_max = 3: the count law 0.1717, 0.7498, 0.0785 against 0.1759, 0.7508, 0.0733; E[τ] 0.015434 against 0.015476 ± 0.000059; all 48 inclusions, the worst at 0.16 of its allowance;
    - fixture C, with its categorical predictor: E[τ] 0.05006 against 0.04999 ± 0.00016; all 90 inclusions, the worst at 0.25.
- *Shared.* `ChainEnsemble.Run` takes `ParallelOptions` (the bound and cancellation), tested in `Maths.Samplers.Tests` (`ChainEnsembleParallelTests`).

**Task 3: diagnostics** (`BassDiagnostics.cs`).
- *What it reports.*
  - Rank-normalized R̂ and bulk and tail ESS for M, τ (null under fixed g), σ², and each draw's fitted value f = a₀ + Σ aⱼBⱼ at the monitored points.
  - Acceptance per move type; refusals by reason; empty births.
  - `BassSupportTally`: each birth's and change's candidate term at each part's order, with its support outcome.
  - The mass at the cap; `SamplesUsed` and `StoppedEarly` for the R̂ stop.
  - Each chain's own account, with its start.
- *The sampling phase only.* Chains reset their tallies at the end of burn-in, and the move counts lose their burn-in baseline, as the shell's acceptance rate does.
- *Checks.*
  - The totals are the chains' sums. A sweep is a move attempt or an empty birth, and the support tally's attempts are the births and changes.
  - Support refusals rise with the order, 0.197 at order 1 against 0.436 at order 2.
  - The mass at the cap is 0 when the cap cannot bind, and the draws' share there when it binds (above 0.9 at a cap of one).
  - Every diagnostic matches `ChainDiagnostics` recomputed from the retained draws.
  - The R̂ stop ends a run early on its batch boundary, with the full run's draws up to there.
  - Under Jeffreys on the exact-fit fixture, twelve chains' starts all lie in the enumerated F.

**Task 4: checkpoints and resume** (`Bass.Resume`).
- *The restart state.* Each chain's state (structure and u), the three chain streams and the composition stream, the counters, burn-in move counts, tallies, term table, draws and start. `Resume` refuses a fingerprint other than the config's and the data's.
- **Deviation from the brief: the design does not reload from the state.**
  - A strict reload of a structure differs from the chain's incrementally updated design in the last bits after deaths and changes: 10⁻¹⁶ in ‖r‖² and 2·10⁻¹⁵ in a coefficient on a direct functional chain.
  - That breaks bit-identity. A mutation check, restoring by reload, fails the functional resume and chain-restore tests.
  - So the restart state carries the design's snapshot: an exact copy of its factor (`BorderedCholesky.Copy`, shared LinAlg, tested in `BorderedCholeskyCopyTests`), and its coefficients and residual as held. A capacity growth replaces the factor without a refit, so a refit on restore could differ too.
  - This needed `BassDesign.Snapshot` and `Restore`, in a BASS file outside my list (internal, no change to existing behaviour).
- *Checks.*
  - A resumed run is bit-identical to the uninterrupted one under eight settings: defaults, functional, categorical, fixed g, prior-cancelling, Jeffreys with the flat intercept, thinned, and batched. That covers every draw, the merged table and every diagnostic.
  - Resuming twice, and resuming a checkpoint a second time, leave the fit it came with unchanged.
  - A config, seed or data change is refused, and a resume reports only its own sweeps.
  - A restored chain continues bit for bit, statistics, compositions and tallies included.

**Suites** (in `artifacts\chip-run\`, at `a402b14`): Regression 683/683; the BASS oracle tests 25/25, `BassFitParityTests`' two included; LinAlg 120/120; Samplers 146/146.

**Review passes.**
1. *Pass 1,* a reread of `b41072d..759136d` on my paths. While building, a test had already caught one defect: under fixed g, a chain's start reported e^{ln τ} where the held τ belongs, fixed before `87b23bc`.
   - *Found.*
     - The sweep's docstring still said a restored structure loads whole.
     - The support tally listed order 0 in fits of one part, where no term takes it.
     - The sampler's refusal message quoted its lifetime tallies.
     - *Gaps:* no check of the sampler through the functional design's two support counts; the default points at n = 1, 2 and 3; restore's failure paths (a foreign state, a state the snapshot does not hold, a snapshot over other data); the tally's merge; the checkpoint and progress options leaving a fit unchanged.
   - *Fixed* in `b14c84d` and `7ac76ae`.
2. *Pass 2* found that the ensemble's `Parallel.For` wraps a chain's exception. A start the sampler refuses then reached callers as an `AggregateException`, the message repeated per chain. `Fit` and `Resume` now rethrow the shared failure as itself (`a402b14`).
3. *Pass 3,* over that delta, found nothing.

**For the main session.**
- *Open: the monitored points.* By default, up to three training cases at the quartile positions of the data's order, rounded away from zero; for a direct fit, at the middle output location ⌊(m − 1)/2⌋. I recommend a declared `BassSampling.MonitoredPoints` (case indices with an optional location index; null takes the defaults), fingerprinted with the config. Points at new inputs would wait for stage 5's prediction.
- *Open: `Resume`'s type.* It takes `EnsembleCheckpoint`, the fit's own field, where the contract shapes write `BassCheckpoint`, a type that does not exist. I recommend amending the contract's text, so a fit's checkpoint resumes as it is.
- *`BassFit` has no `Evidence` yet.* `BassEvidence` is S8's; chip D adds the field between `Diagnostics` and `ResidualDependence`. `ResidualDependence` is there (chip B's type), left null.
- *For chip D.*
  - `BassPriorSampler` is the reference rung's sampler.
  - A ladder needs per-replica restart state: `BassStructureMoves.Snapshot(state)`, its restore constructor and the design's snapshot are internal and ready.
  - The run's `StartChain` is where a chain becomes a ladder.
- *Smaller items, left with their reasons.*
  - `BassStructurePrior`'s subset draw repeats `BassTermProposal`'s private one, with the same algorithm and stream use. The proposal could delegate to it; that file is not mine.
  - The R̂ stop's two-chain rule is checked in `Bass.Fit`, naming `Sampling.Chains`. It could live in `BassConfig.Validate`, S1's file.
  - A start's rejection loop does not check cancellation; the attempt limit bounds it.
  - The tail R̂ of a fitted value is sensitive to the last bits of the summation order, near ties about the median. To reproduce it from the draws, sum from a₀ in the ids' order.
