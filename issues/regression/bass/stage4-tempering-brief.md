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

### Chip D: S8, tempering and evidence — 2026-10-28

**Landed.** Both tasks, then the review cycle, all on `main`: `58fdeb2` (task 1, the tempered run), `8c1bb66` (task 2, the evidence), `fb22f73` and `f2d7149` (review pass 1), `aaf46b6` (pass 2).

**Task 1: the tempered run.**
- *`TemperingSpec`* in `BassSampling`, null by default (untempered, D05), a values-only record, serialized and fingerprinted with the config:
  - `Path`: likelihood by default, powered selectable;
  - `Form`: integrated by default, linear selectable;
  - `Levels` (8) and `HotEnd` (10⁻³): geometric from 1, then b = 0 with the reference rung;
  - `Adaptive` (true): re-spaced during burn-in, frozen after it;
  - `Schedule`: DEO;
  - `ReferenceRung` (false).
- *Refusals.*
  - `Validate` refuses an unknown path or form, the linear form on the powered path, the reference rung off the likelihood path, too few levels (2, or 3 with the rung), a hot end outside (0, 1), and any schedule but DEO.
  - `BassResolved` refuses a hot end where the integrated form is undefined (under the flat intercept, b ≤ (1 − 2c)/N), and the reference rung where the integrated form has no b = 0 (Jeffreys or the flat intercept). It reads a static `BassTarget.Defines(form, intercept, σ², N, b)`, which the instance method now calls.
  - `BassModel.Create` takes the config's form when none is passed.
- *`BassLadder`* (new) turns a chain into a ladder:
  - one `ParallelTempering` per chain, and per replica its own `BassStructureMoves` and `BassTauBlock`;
  - swaps on ℓ_b alone, the prior cancelling;
  - the reference rung drawing from its own `BassPriorSampler`;
  - an adaptive ladder re-spaced after rounds of 1, 2, 4, … sweeps (BARS's schedule), then frozen: `Adapt` refuses after `EndBurnIn`;
  - its seed from a fourth `SeedTree` child of the chain seed. The children are prefix-stable, so untempered runs are unchanged bit for bit.
- *One departure from S6's plan: designs travel with swapped states.* S6 loads a swapped-in state whole.
  - A whole load differs in the last bits from the giver's incremental design. Near D06's angle floor it could even refuse a structure the giver admitted, which would throw mid-run.
  - So the receiving replica adopts an exact copy of the giver's design (internal `BassStructureMoves.Holds`, `DesignFor`, `Adopt`), at O(N·M + M²) against a load's O(N·M² + M³).
  - Only reference draws are loaded whole, by the same `Load` that admitted them in the sampler.
- *Retention and diagnostics.*
  - Draws are composed from the b = 1 replica, and the moves, refusals and support account are that replica's.
  - `BassDiagnostics.Tempering` is a `BassTemperingReadout`, following BARS's: path, form, schedule, the reference flag, round trips and annealed restarts per scan, each pair's mean rejection (Rao–Blackwellized) and realized swap acceptance, the barrier Λ̂, Theorem 1's ideal round-trip rate, each chain's frozen ladder, and the reference sampler's attempts per draw.
- *Checkpoints* carry the shared `LadderSnapshot` (each replica's state, rung, streams and counters, the swap and reference streams, the adaptation's state and the traces), every replica's moves (its design as it is), the reference sampler's account and the burn-in marks. Resume is bit-identical.

**Task 2: the evidence.**
- *When.* `BassFit.Evidence`, between `Diagnostics` and `ResidualDependence`, is reported when `ReportEvidence` is set, the ladder ends at the reference rung, and the priors are proper. It is null under the flat intercept or Jeffreys σ², without the rung, on the powered path, and untempered.
- *The estimators.*
  - Stepping-stone gives `LogEvidence`. Each chain's ladder gives an estimate over the rungs it froze, and the chains pool as BARS pools them: the mean, with error √Σ err_c² / C.
  - Path sampling runs beside it on the linear form, where the slope is ℓ₁ and the discretization diagnostic holds.
  - *On the integrated form path sampling is left out:* its slope at b = 0 carries −(k/2)·E[1/τ], infinite under Zellner–Siow, so the trapezoid's end has no finite mean; elsewhere its diagnostic misses E_b[∂²ℓ_b/∂b²] (D05's first point). Docstring on `BassEvidence`.
- *`BassEvidencePriors`* records the coefficient prior as declared and as resolved (τ's Gamma or the held τ), the count prior with the M cap, the subset measure and the resolved interaction caps, the support rule with its resolved thresholds and D06's gates (they define F, where the prior is normalized), the centre and its level, and σ²'s prior as declared and as resolved.

**Tests.**
- *`BassTemperingTests`* (52 cases):
  - the spec's contract and every refusal;
  - adoption bit for bit, and every replica holding its state's design;
  - each rung's law against enumeration, over 8 options: both forms and paths, the reference rung, fixed g, flat + Jeffreys on the linear form and on the powered path, and a rung at the flat intercept's limit. Each b = 0 rung's τ matches the prior's closed form;
  - the readout over 10 options: two and three rungs, a fixed ladder, functional, categorical, the R̂ stop;
  - a ladder's scan as one sweep;
  - the evidence against enumeration's exact ln Z under 5 proper options and on a direct functional response, both forms;
  - the evidence null where it is undefined or not asked for.
- *`BassTemperingResumeTests`* (15): resume over 11 options, resume twice, the parallelism bound, and a ladder's snapshot before and after burn-in.
- *`BassTemperingParityTests`* (10, oracle):
  - *The b = 1 end on fixture A, both forms:* tempered fits against the exact posterior, all 48 inclusions within allowance.
  - *Each rung on fixture B (b = 1, 0.316, 0.1, 0), both forms:* the count law, E[τ] and every inclusion against the oracle's tempered laws, within five batch errors.
  - *The swap's exact balance (fixed g, fixture A, 299 states, both forms):* ln π_b₁ − ln π_b₂ − (ℓ_b₁ − ℓ_b₂) spreads at most 3.2·10⁻¹⁴ over F against the oracle's laws.
  - *R's swap density, a defect fixture:* on the `swap` slice it misses its own rung laws by 1.88 (1 ↔ 0.5) and 1.13 (0.5 ↔ 0.2), and the four candidate families by 1.10–5.08, where D05's misses its own by 1.4·10⁻¹⁴.
  - *The evidence:* fixture B, integrated −28.551 ± 0.058 and linear −28.596 ± 0.060 (path sampling −28.681 ± 0.066, diagnostic 0.078) against −28.589; fixture C, categorical, −30.935 ± 0.040 against −30.903.
- *Suites,* in `artifacts\chip-tempering\` (deleted since): Regression 756/756 (67 new), Samplers 146/146, the BASS oracle tests 35/35 (10 new).

**Review passes.**
1. *Pass 1* found four defects and several gaps.
   - *The defects* (`fb22f73`):
     - the reference rung's replica was handed no design, so a checkpoint reloaded a state swapped up to it whole;
     - while replicas share their start's structure, a giver that handed over its pending candidate lost its own state's design (now a second round);
     - `Adopt` replaced the design before checking it;
     - `BassLadder` refused a bad spec in `ParallelTempering`'s parameter names.

     The docs now say that the hot end moves with the reference rung, and that an unreached pair reads NaN.
   - *The gaps* (`f2d7149`): two- and three-rung ladders, the R̂ stop, progress counting, a rung at the flat intercept's limit, the evidence on a direct functional response, and categorical evidence against the oracle.
   - *A finding, not a defect.* With the default 8 levels on the N = 175 functional fixture (`BassRunTests.FunctionalData`, at its cap of 3), adjacent pairs rejected 0.96–1.00 of swaps.
     - The chains' stepping-stone estimates spread from 122 to 183 against errors of 2–5, and the two forms differed by 130 nats in a short run.
     - At 16 and 32 levels both forms agreed near 180–186.
     - So the reported Monte Carlo error under-covers when rungs barely overlap, and a fixed level count does not scale with the data (open question 1).
2. *Pass 2* found no defect. Its two gaps, the powered path under the improper options and the b = 0 τ against the prior's closed form, are closed in `aaf46b6`.
3. *Pass 3,* over that delta and the whole change, found nothing.

**Outside my listed files.** No change to `Maths.Samplers`. Edited:
- `BassTarget.cs`: the static `Defines`;
- `BassResolved.cs`: the rung refusals;
- `BassModel.cs`: the form from the config;
- `BassPriorSampler.cs`: an internal save and restore of its account;
- `BassSweep.cs`: the hand-over's three internal methods;
- `tests/oracle/BassTargetParityTests.cs`: `FixtureB` private → internal.

The config's JSON now carries `"Tempering":null`, so every fingerprint moved. That is pre-release, and no test pinned one.

**Open questions for the main session** (each implemented as a reversible default):
1. *The default levels and hot end,* 8 and 10⁻³.
   - On the small fixtures Λ̂ was 1.5–1.9, where Syed et al.'s N ≈ 2Λ asks about 5 rungs, and stepping-stone was within its error at 6–16 levels.
   - The functional finding above shows a fixed count fails as N grows.
   - I recommend an auto-pick from the burn-in's barrier, `TemperingLadder.RecommendedLevels(Λ̂)`, logged in `Resolved`. It needs a decision on how burn-in grows the ladder.
   - As a cheaper guard, the evidence could be flagged when a pair's rejection nears 1.
2. *Swap cadence:* after every sweep. It is `ParallelTempering`'s only cadence and BARS's, and the spec has no field for it.
3. *Path sampling on the integrated form:* left out, for the reasons above. Restoring it needs a curvature seam, E_b[∂²ℓ_b], for its diagnostic, and a finite end in place of b = 0.
4. *The schedule's closed set* holds DEO only, as D05 decided; SEO and the sequential schedule are refused, though BARS exposes them.
5. *The linear form on the powered path is refused,* since the form is the likelihood path's. The default integrated form is accepted there and ignored.
6. *With the reference rung, an adaptive ladder moves the hot end* (only b = 1 and b = 0 are fixed), as BARS's does. `HotEnd` is then the start.
7. *The integrated form's slope is traced at every rung every scan,* since `ParallelTempering` traces slopes whenever the target supplies them, though the evidence does not read it there. That costs one marginal per rung per scan; a trace switch in the shared ladder would save it.
8. *The diagnostics' defensive check* that a checkpoint's chain ran tempered as the config does is unreachable through `Bass.Resume`, since the fingerprint covers the config. It stays untested.
