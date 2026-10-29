# BASS stage 4: a self-sizing ladder and an honest evidence error (chip E brief)

Date: 2026-10-28. Follows chip D's S8 ([its report](stage4-tempering-brief.md#reports)). The user decided its first open question on 2026-10-28; the decision is the fitting checklist's report "S8 verified; chip D's questions decided", and this brief builds it.
- *What exists.* S8 is on `main`: `TemperingSpec`, `BassLadder`, `BassEvidence`, the tempering readout, and bit-identical resume with tempering on.
- *The finding this chip answers.* With the default 8 levels on the functional fixture (`BassRunTests.FunctionalData`, N = 175, at its cap of 3), adjacent pairs rejected 0.96–1.00 of swaps, and the chains' stepping-stone estimates spread from 122 to 183 against errors of 2–5. Two things fail:
  - *A fixed count cannot serve.* For a power posterior the barrier grows as Λ ≈ √(d/2π)·ln N (a Bernstein–von Mises heuristic). BASS's swaps read the collapsed ℓ_b, so its effective d is smaller, but the shape stands: the count must follow the data.
  - *The reported error under-covers* when rungs do not communicate, because it is pooled from within-chain errors alone.
- *This chip builds:*
  1. a ladder resize in the shared `ParallelTempering`, the placement from the caller;
  2. BASS's self-sizing ladder and its hot-end auto-pick;
  3. the evidence's between-chain guard, shared with BARS;
  4. the slope-trace switch (chip D's Q7).
- *Not in scope.* SMC (stage 8). BARS's ladder options and defaults: BARS gains only the pooled error. S10, which the main session runs beside you.

Read first: [AGENTS.md](../../../AGENTS.md) and [CONTRIBUTING.md](../../../CONTRIBUTING.md). Then:
- chip D's brief and report, in full, for what S8 built and why designs travel with swapped states;
- the fitting checklist's report "S8 verified; chip D's questions decided", which is this chip's decision;
- the fitting segment's D05, and the `TemperingSpec` bullet in its contract shapes (amended 2026-10-28);
- D10's adaptation note in the [numerics segment](port-numerics-brief.md#d10-note-adopted-2026-10-27-a-policy-for-adaptive-kernels): adapt during burn-in, freeze after it;
- `docs/Maths/Samplers.md` on tempering (rung count, `Equalize`, the barrier);
- [`issues/spc/opus-spin-glass-threads.md`](../../spc/opus-spin-glass-threads.md), threads 5 and 6, for why the placement is the caller's. SPC on parallel tempering will need this resize, and its rungs may be placed by acquisition, for a readout's resolution, and not by mixing.

In the code:
- *Shared:* `ParallelTempering.cs` (its `SeedTree` layout at the constructor, `AdaptLadder`, `Snapshot` and `Restore`, `TraceRungs`), `TemperingLadder.cs` (`Barrier`, `Equalize`, `RecommendedLevels`), `LadderEvidence.cs`, and BARS's pooling in `BarsEnsemble.cs` (`Pooled`, near line 558).
- *BASS:* `BassLadder.cs` (the hand-over: `Holds`, `DesignFor`, `Adopt` in `BassSweep.cs`), `BassRun.cs` (`BassChainRun.Burn`, the doubling rounds), `BassDiagnostics.cs` (`BassTemperingReadout.Of`, which pools by pair index at line 88), `BassEvidence.cs`, `BassResolved.cs` (the ladder's refusals) and `BassConfig.cs` (`TemperingSpec`).

## Working rules (a shared checkout)

- **Work in `D:\aghado01\ThermoMapper` on `main`.** Worktrees fail on this repo's long paths. Switch no branch, push nothing.
- **Build and test into your own artifacts root:** add `-p:ArtifactsRoot=D:\aghado01\ThermoMapper\artifacts\chip-ladder\` to every `dotnet build` and `dotnet test`. Do not run `scripts/fact-harness.ps1`. Delete `artifacts\chip-ladder\` when you finish.
- **Run oracle tests filtered** to BASS (`--filter "FullyQualifiedName~Bass"`): the TDA oracle tests time out when several builds run at once.
- **The main session works beside you, on S10.** It owns `BassResidualDependence.cs`, `BassPointwise.cs`, `BassFit.cs`, `docs/Maths`, the fitting brief and checklist, and new S10 test files (`tests/regression/BassRecovery*`, `BassQualification*`, `BassCost*`, and S10's oracle fixtures). Leave those alone. If you need a change in one, name it in your report. Expect its commits between yours; they touch none of your files.
- **Your files.**
  - *Shared:* `src/maths/samplers/rjmcmc/ParallelTempering.cs`, `TemperingLadder.cs`, `src/maths/samplers/evidence/LadderEvidence.cs`, `src/maths/regression/spline/bars/BarsEnsemble.cs` (its evidence pooling only), and their tests in `tests/samplers/` and BARS's evidence tests.
  - *BASS:* `BassConfig.cs` (`TemperingSpec`), `BassResolved.cs`, `BassLadder.cs`, `BassRun.cs`, `BassDiagnostics.cs`, `BassEvidence.cs`, `BassSweep.cs` (the hand-over's internals), and `tests/regression/BassTempering*Tests.cs`, `tests/oracle/BassTempering*ParityTests.cs`.
- **Shared code stays generic.** It is tested in its own host, and every change to it is named in your report. A ladder that never resizes must stay bit-identical, BARS's included: rerun BARS's tempered tests and the samplers' snapshot tests unchanged.
- **Commit your own files only, by explicit path** (`git add <paths>`, never `-A` or `.`), one commit per task, in the repo's message style, ending with `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`.
- **Line endings:** briefs and C# are LF; `docs/Maths/*.md` are CRLF.
- **Decisions are not yours to settle.** The checklist's report fixes the design. Where it leaves a constant open (the deadband's width, the flag's rejection threshold, the initial count, the minimum count, how a new rung picks its source state), implement a reversible recommendation and put the question in your report with your reasons. Do not re-open decided choices.
- **Report** under "Reports" at the end of this brief, newest first: what landed, the commits, the tests and suite counts, the review passes (below) and anything the main session must know.

## Task 1: the shared resize (`Maths.Samplers`)

**The resize.** `ParallelTempering` gains a burn-in resize: the caller gives the new rungs (a descending β array, the cold end and, with the reference rung, b = 0 kept) and, for each new rung, the replica whose state it starts from.
- Kept replicas keep their states and streams.
- A new replica starts from a clone of its source's state, made by a caller-supplied clone, since a model's state may need a deep copy. BASS needs one: see task 2.
- The swap statistics, the round-trip labels and the adaptation round's sums restart on the new ladder, since the old pairs do not map onto the new ones. The sampling-phase readout is unaffected: its marks are taken when burn-in ends.
- It is refused once the ladder is frozen. Name the freeze where `ParallelTempering` can see it, or keep the refusal in the callers, and say which in your report.

**Streams.** New replicas' move and block streams are drawn from one growth stream, appended as the last child of the existing `SeedTree` layout. `SeedTree` is prefix-stable, so every ladder that never resizes keeps its streams bit for bit. Snapshots carry the growth stream and the resized ladder, and a ladder restored after a resize continues bit for bit.

**The policy, in `TemperingLadder`.**
- *Placement.* K rungs placed on the cumulative barrier Λ(β): `Equalize`'s monotone cubic, inverted at K equal steps in place of the current count. `Equalize` is the case K = the current count.
- *The sizing rule:*
  - grow at once to `RecommendedLevels(Λ̂)` when it exceeds the current count. A saturated ladder reads Λ̂ ≤ pairs, so this roughly doubles it, and it keeps growing until the reading stops asking;
  - shrink only past a deadband (the proposal is a one-rung deadband: shrink when ⌈2Λ̂⌉ + 1 ≤ current − 2), since early rounds over-read;
  - clamp to a minimum and to the caller's maximum.

  Expose it as a pure function of (Λ̂, current, minimum, maximum), so BARS and SPC can call it.

**The slope switch.** `TraceRungs(bool traceSlopes = true)`. False skips the slope per rung per scan where no estimator reads it.

**The pooled evidence.** A helper in `Maths.Samplers.Evidence` pools per-chain `LadderEstimate`s:
- the mean;
- the within-chain error √Σs_c² / C;
- the between-chain error, the standard error of the chains' mean, √(Σ(ẑ_c − z̄)² / (C(C − 1))) for C ≥ 2;
- the reported error, the larger of the two, or the within-chain error alone at C = 1;
- the mean discretization diagnostic.

BARS's `Pooled` uses it, and its estimate record gains both components. That changes BARS's reported error and nothing else.

## Task 2: BASS's self-sizing ladder

**The spec.**
- `TemperingSpec.Levels` becomes `int?`. Null, the default, self-sizes from an initial count (8, open), and an explicit count fixes it as now.
- `MaxLevels` (64) is a compute budget.
- `HotEnd` becomes `double?`. Null, the default, auto-picks.
- `Validate` refuses:
  - a self-sizing ladder that does not adapt (`Adaptive` false with `Levels` null);
  - `MaxLevels` below the minimum, or below a fixed `Levels`;
  - and keeps its current refusals.

**The hot end's auto-pick,** resolved in `BassResolved` and logged as an auto-pick:
- 1/N, where b = E[τ] makes the shrinkage τ/(τ + b) one half: the tempered likelihood weighs one observation, as the unit-information prior does;
- 2/N under the flat intercept on the integrated form, whose pole is at (1 − 2c)/N, so that ν = Nb − 1 = 1;
- under fixed g, the smaller of the held τ and the above.

The resolved hot end is what the existing refusal (`BassTarget.Defines`) checks. A hot end that does not land in (0, 1) on tiny data is refused with a message to set it. With the reference rung the hot end stays the adaptation's start (chip D's Q6, accepted).

**The burn-in** (`BassChainRun.Burn`).
- At each doubling round except the last, the ladder reads Λ̂, sizes by the rule within `MaxLevels`, and places its rungs.
- The last round, the longest, equalizes at the settled count. With R's 9000-sweep burn-in it spans 4096 sweeps, and only 809 follow it.
- A fixed `Levels` runs exactly as S8 does: pin it, bit for bit, against values you record before the change.

**Clone-safe replicas.** A replica's structure identity is its term array, by reference (`BassState.SharesStructure`).
- A new replica takes its own copy of its source's term array and adopts its source's design (`DesignFor`, then `Adopt`, which compares content, not reference).
- *Why:* if two replicas held one mid-run array and both left it, the structure would be reloaded whole. A whole load differs from the incrementally updated design in the last bits, and near D06's angle floor it could refuse a structure the source admitted.
- The start shares one array across all replicas safely only because the start was itself loaded whole. State that invariant in the code.
- After any resize, no two replicas share a mid-run term array.

**The readout** (`BassTemperingReadout`).
- Each chain's pair rejections and acceptances, since chains freeze at different counts, beside its frozen ladder, which is already per chain.
- The pooled barrier, the chains' mean Λ̂.
- Each chain's level count.
- The largest pair rejection.
- A flag when any chain's ladder reached `MaxLevels`.

**The evidence** (`BassEvidence`).
- `BassEvidenceEstimate` gains the within- and between-chain errors, and its `MonteCarloError` is the larger of the two, through task 1's helper.
- The evidence also carries the flags: the largest pair rejection, against a threshold you recommend; the cap's binding; and the between-to-within ratio.
- A flagged evidence is still reported, with its flags, never nulled.

**The slope switch:** `BassLadder.EndBurnIn` traces slopes on the linear form only.

## Checks

- *The resize, in `tests/samplers/`:*
  - kept replicas keep their states and streams;
  - new replicas take the growth stream;
  - a ladder that never resizes is bit-identical to before;
  - snapshot and restore across a resize are bit-identical.
- *Placement and sizing:*
  - K rungs on a synthetic barrier profile give equal increments, and K = the current count reproduces `Equalize`;
  - the sizing rule's grow, deadband and clamp cases, as a table.
- *The pooled error:* on synthetic chain estimates, between against within, the larger reported, and C = 1 falling back to within. BARS's evidence tests still pass, with the new components checked.
- *BASS:*
  - On the functional fixture with self-sizing, no pair saturates, the two forms agree within their errors, and the chains agree. The S8 report's 16- and 32-level runs agreed near 180–186, a cross-check.
  - Fixture B's evidence against the oracle's exact ln Z with self-sizing, on both forms.
  - A fixed `Levels` run is unchanged from S8.
  - Resume is bit-identical across a resize: checkpoints before and after the resizing rounds, and after burn-in.
  - The hot end's auto-pick under each centre and under fixed g, logged.
  - Every refusal.
  - The clone invariant: after resizes, every replica holds its state's design, and no two replicas share a mid-run term array.

## Review before reporting

A task is not done when its tests pass. Before the report, review the work as a separate pass, as the main session reviewed S5 and S6 and the earlier chips reviewed their own; every one of those reviews found real defects or untested paths.

1. **Reread the whole change** as a reviewer seeing it for the first time: `git diff <the commit before your first>..HEAD -- <your paths>`. For each change, ask:
   - Is it correct? Check the resize's state and stream bookkeeping, the clone's independence, the placement at the ends (b = 1, the hot end, b = 0 with the reference rung), the sizing rule at its boundaries, the pooled errors at C = 1 and 2, the exception paths and their parameter names, thread safety, and the contracts it implements.
   - Does it read like the surrounding code?
   - Do its docstrings say what it does?
2. **Hunt the coverage gaps.** Look for each option and convention the code accepts that no test runs:
   - self-sizing and fixed ladders, with and without the reference rung;
   - both forms and both paths;
   - a cap that binds;
   - a short burn-in with few rounds;
   - sampled and fixed g;
   - each intercept convention for the hot end;
   - scalar, categorical and direct functional data.

   Look also for each refusal path no test reaches. Prefer independent references (the oracle's exact ln Z, enumeration) to self-consistency.
3. **Fix what you find,** add the tests and rerun. Commit the fixes separately (`fix(...)` or `test(...)`).
4. **Review again** if the pass found anything substantive. Stop when a pass finds nothing.
5. **Record it** in your report: what each pass found, what was fixed, and anything left, with the reason.

For the main session's S10: report each self-sized run's settled count with its N and the cold chain's mean M. The cost baseline will check the count's growth against the heuristic.

## Reports

### Chip E: a self-sizing ladder and an honest evidence error — 2026-10-28

**Landed.** Both tasks, then the review cycle, all on `main`:
- `ec565f6`: task 1, the shared resize, placement, sizing, slope switch and pooled error;
- `d64ee4e`: task 2, BASS's self-sizing ladder, hot end and flagged evidence;
- `66c66be` and `abf16cf`: review pass 1;
- `36ce26a`: pass 2.

**Task 1: the shared changes** (`Maths.Samplers`, and BARS's pooling), each generic and tested in its own host.
- *`ParallelTempering.Resize(betas, sources, fork, movesForReplica, blocksForReplica)`*, burn-in only. The caller places the rungs and names each new rung's source.
  - *Kept replicas.* A replica a rung names is kept, with its state, streams, moves, blocks and counters, at the naming rung nearest its own β (the colder on a tie). A replica no rung names is dropped. The cold rung must keep the cold replica, and with the reference rung the last rung must keep the reference replica.
  - *New replicas.* Each further rung naming a replica starts a new one from `fork` of its state, an independent copy as SMC's `Fork` is, evaluated once against the ladder's target, which the ladder now keeps. Its move stream, then its block stream, take the next draws of the growth stream, cold first.
  - *The growth stream* is the `SeedTree`'s last child, after the replica, swap, reference and block streams, so a ladder that never resizes keeps every stream. Pinned at `2bdc22a` with and without blocks and the reference rung; BARS's tempered tests and the samplers' snapshot tests pass unchanged.
  - *Statistics.* The pairs' statistics, the round-trip labels and the adaptation round restart; the totals carry on. `PairStats` averages over the new `PairScans` (the scans since the count was set, which is `Scans` for a ladder that never resizes, so its values are unchanged). `RoundRejection` reads the round's rates.
  - It returns, per new rung, the replica kept there or −1. A refused resize leaves the ladder as it was, its growth stream included.
  - *The freeze stays the callers'* (`BassLadder.Frozen`). The shared ladder refuses `Resize` only while it traces its rungs, as it refuses `AdaptLadder`.
  - `LadderSnapshot` gains `PairScans` and `GrowthRng`.
- *`TemperingLadder`*:
  - `Place(betas, rejection, K)` puts K rungs on the cumulative barrier, `Equalize`'s monotone cubic inverted at K equal steps; `Equalize` is now `Place` at the current count, bit for bit. Where the rates cannot place another count (no rejection, or too coarse), it places it as though every pair rejected alike, by the rungs' index.
  - `Sources(betas, rejection, placed)` gives each placed rung's nearest current rung on the barrier, then in β, then the colder.
  - `Size(Λ̂, current, minimum, maximum)` grows at once to `RecommendedLevels(Λ̂)`, shrinks to it only when it asks for at least two fewer (`ShrinkDeadband` = 1), and clamps.
- *`TraceRungs(bool traceSlopes = true)`*: false traces NaN slopes without querying them.
- *`LadderEvidence.Pool(chains)` → `PooledLadderEstimate`*:
  - the mean, the within-chain error √Σs_c²/C, and the between-chain error √(Σ(ẑ_c − z̄)²/(C(C − 1))), NaN at C = 1;
  - the larger of the two reported, or the within-chain error alone at C = 1 (NaN propagates);
  - the mean discretization diagnostic, and `BetweenToWithin`.
- *BARS*: `Pooled` uses it, and `BarsEvidenceEstimate` gains `WithinChainError` and `BetweenChainError`. Its reported error is the only value that moves.

**Task 2: BASS's self-sizing ladder.**
- *`TemperingSpec`.*
  - `Levels` (`int?`, null by default) self-sizes from `InitialLevels` (8) within `MaxLevels` (64). `HotEnd` (`double?`, null by default) auto-picks. `StartLevels()` and `MinimumLevels(referenceRung)` (2, or 3 with the rung) are new.
  - `Validate` adds three refusals: a self-sizing ladder that does not adapt (`Sampling.Tempering.Adaptive`); `MaxLevels` below the minimum; and `MaxLevels` below a fixed `Levels`.
  - The config's JSON now carries `"MaxLevels"`, so every tempered fingerprint moved; no test pinned one.
- *`BassResolved.HotEnd`*, null untempered, is the declared hot end or the auto-pick (`AutoHotEnd`), logged as `Tempering.HotEnd` with its reason. A self-sizing ladder's start is logged as `Tempering.Levels`.
  - The pick: 1/N; 2/N on the likelihood path's integrated form under the flat intercept; under fixed g, the smaller of the held τ and those.
  - Refused: an auto-pick outside (0, 1), with "set it"; and, by the existing `Defines` check, now reading the resolved value, an auto-pick the integrated form cannot take (a held τ below the flat intercept's pole), with "set it above".
  - The tempering refusals moved after τ's resolution, which the pick reads.
  - `BassLadder` makes the same pick when its model was resolved without the spec, as tests build ladders, and refuses it in the spec's name.
- *The burn-in* (`BassChainRun.Burn`). An adaptive self-sizing ladder calls `BassLadder.Resize` after every round but the last (`done + 2·round ≤ steps`), and `Adapt` after the last, which equalizes the settled count. A fixed `Levels` runs `Adapt` every round, as S8 did: seven fixed-ladder runs (both forms and paths, fixed g, functional, categorical, a non-adaptive ladder) are pinned bit for bit at `2bdc22a`.
- *`BassLadder.Resize`* reads Λ̂ from `RoundRejection` and sizes by `TemperingLadder.Size` within [`MinimumLevels`, `MaxLevels`]. At an unchanged count it is `AdaptLadder`, bit for bit. Otherwise it places the count (`Place`), picks the sources (`Sources`) and resizes the shared ladder.
- *Clone-safe replicas.*
  - `Fork` gives a new replica its own copy of its source's term array, with the source's fit and τ.
  - Its moves are built on the source's state and design (`DesignFor`, restored with empty tallies), then `Adopt` the fork, which compares content.
  - A source holding a reference draw it never loaded loads it whole first, by the same admission, as a swap's receiver would. The forked replica's design then equals a whole load bit for bit (tested).
  - The ladder's constructor states the start's invariant: sharing its array is safe only because the start was loaded whole.
- *`EndBurnIn`* traces slopes on the linear form only, and marks the pairs' sums over `PairScans`.
- *`BassTemperingReadout`* now has:
  - `PairRejection` and `PairAcceptance` per chain (`double[][]`);
  - `Barrier`, the chains' mean Λ̂, and `MaxPairRejection`;
  - `IdealRoundTripRate`, the chains' mean of each chain's;
  - `Levels`, each chain's count, and `AtMaxLevels`, set when a self-sizing chain's count reached its budget.
- *`BassEvidence`* gains `MaxPairRejection`, `AtMaxLevels` and `Saturated` (above `RejectionThreshold` = 0.9). Its estimates gain `WithinChainError`, `BetweenChainError` and `BetweenToWithin`, the error from `LadderEvidence.Pool`. A flagged evidence is still reported.

**The functional fixture** (N = 175, M_max 3, the reference rung, 4 chains of 4,000 burn-in and 4,000 draws):

| ladder | settled levels | max rejection | ln Z (within, between) |
| --- | --- | --- | --- |
| self-sized, integrated | 14, 15, 11, 12 | 0.902 | 192.59 ± 0.52 (0.35, 0.52) |
| self-sized, linear | 12, 16, 13, 13 | 0.837 | 193.45 ± 0.41 (0.41, 0.17); path sampling 193.42 ± 0.53 |
| fixed 32, integrated | 32 | 0.634 | 193.40 ± 0.37 |
| fixed 64, integrated | 64 | 0.348 | 193.03 ± 0.30 |
| fixed 16, integrated | 16 | 0.955 | 192.64 ± 0.59 |
| fixed 8, integrated (S8's default) | 8 | 1.000 | 189.68 ± 1.51 (0.43, 1.51) |

- The two forms agree, the chains agree (between/within 1.5 and 0.4), and both agree with the fixed 32 and 64 rungs.
- *S8's cross-check moves.* S8's 16- and 32-rung runs read 180–186; at this budget 32 and 64 rungs read 193.0–193.4. Those were short runs, low as stepping-stone reads where rungs barely overlap.
- At eight rungs the between-chain error now reports ±1.51 where the within-chain error gave ±0.43. It still misses by 2.3 of its errors, but is no longer 8 off.

**Tests.**
- *Samplers*, 40 new:
  - `ParallelTemperingResizeTests` (14): the never-resized pins; kept and new replicas' states and streams against the growth stream; the kept-replica rule; every rung's law after a resize against quadrature, on both paths and with the rung; snapshots before and after a resize, a shrink and the traces; every refusal, the ladder unchanged after a factory throws mid-resize; the slope switch.
  - `TemperingLadderTests` (24): placement on a known barrier at 5–33 rungs, equal to within 0.6%; `Place` at the current count is `Equalize`; the index fallback; the sizing table (15 cases) and its refusals; the sources against a brute-force reading.
  - `LadderEvidenceTests` (2 new): the pooled errors at C = 1, 2 and 3, and the refusals; and a calibration. Four chains spreading with sd 1 while each reports 0.1: the within-chain error covers 0.164 of the time, the reported one 0.849 (a t on 3 degrees of freedom gives 0.855).
- *BARS*: `BarsEvidenceTests` checks both new components and the larger reported.
- *Regression*, 44 new:
  - `BassTemperingSizingTests` (23): the S8 pins; the functional fixture above; the hot end's pick under each centre, the flat pole, both forms and paths, Jeffreys and fixed g, logged; its refusals; the ladder's refusals; the clone invariant after resizes; a resize that forks the reference rung's draw; a cap that binds (`MaxLevels` 4: max rejection 1.000, 67 ± 49, between/within 47, flagged, still reported); short burn-ins of 0–100 sweeps, replayed by hand from the chain's seeds, bit for bit; slopes on the linear form only.
  - `BassTemperingTests`: 7 self-sizing run options in the readout (scalar, linear, powered, fixed g, functional, categorical, flat) and 4 new refusals. The exact-ln Z evidence now runs self-sized under 5 options, against enumeration.
  - `BassTemperingResumeTests`: the 7 self-sizing options resumed after burn-in, and the ladder restored across its resizes (3: with and without the rung, fixed g), before and after a sizing round and after burn-in.
- *Oracle* (`BassTemperingParityTests`, self-sized by default now): fixture B against the oracle's exact −28.589, linear −28.554 ± 0.070 and integrated −28.609 ± 0.086; fixture C, categorical, −30.882 ± 0.075 against −30.903. None flagged.
- *Suites*, in `artifacts\chip-ladder\` (deleted since): Samplers 186/186 (40 new) and Regression 812/812 (44 new; the count includes the main session's S10 tests present in the checkout) at `36ce26a`; the BASS oracle tests 35/35 at `d64ee4e`, since when only docs, tests outside the oracle project and a bound no run reaches changed.

**Review passes.**
1. *Pass 1,* partly before task 2's commit:
   - *Defects:* the self-sizing start was clamped in two places and its log ignored the clamp (`StartLevels`); `Burn`'s `done + 2·round` overflowed past 2³⁰ sweeps; and stale docs on `AdaptLadder`'s round and `SwapPairStat`'s scope after a resize.
   - *Gaps:* no run reached a resize that forks the reference rung's draw, whose design its replica has not loaded (a dedicated test now forces one); no self-sized run under the flat intercept; the ladder's own auto-pick refusal; a ladder restored across resizes under fixed g; the clamped start's log.
2. *Pass 2* found two gaps: a resize's rung laws on the powered path, and sources for a rung placed beyond the ends.
3. *Pass 3,* over that delta and the whole change, found nothing.

**Settled counts for S10** (each self-sized run, with N and the cold chain's mean M):

| fixture | N | M_max | burn-in | settled levels | mean M |
| --- | --- | --- | --- | --- | --- |
| `BassSweepTests.Data` | 18 | 2 | 1,000 | 4–6 (integrated 5, 5, 6, 6; linear 5, 5, 5, 5; zero centre 5, 4, 5, 6; fixed g 4, 4, 4, 4; declared σ² 5, 5, 5, 6) | 1.07–1.17 |
| fixture C, categorical | 27 | 2 | 1,000 | 5, 6, 5, 5 | 1.84 |
| fixture B | 32 | 2 | 1,000 | 6, 6, 7, 6 (linear); 7, 7, 7, 7 (integrated) | 1.95–1.96 |
| `FunctionalData` | 175 | 3 | 4,000 | 14, 15, 11, 12 (integrated); 12, 16, 13, 13 (linear); 13, 13, 15, 15 (fixed g) | 3.000, at the cap (sampled g) |

The ladder reads the barrier high during burn-in: on the functional fixture the rounds read Λ̂ 5–7.5 where the sampling phase reads 3.8. So it settles above what sampling needs (11–17 rungs where 2Λ̂ + 1 is about 9). That is safe, but it costs compute, and the cost baseline should see it.

**Outside my listed files.** `tests/samplers/ParallelTemperingBlockTests.cs`: its fixtures went private → internal, for the resize tests. `BassSweep.cs`: one docstring sentence.

**For the main session.**
- *`docs/Maths/Samplers.md`'s* API table and tempering section need the new members: `Resize`, `RoundRejection`, `PairScans`, `TraceRungs(traceSlopes)`; `TemperingLadder.Place`, `Sources`, `Size` and `ShrinkDeadband`; `LadderEvidence.Pool` and `PooledLadderEstimate`.
- *S8's docs paragraph* can now describe the self-sizing ladder, the auto-picked hot end, the per-chain readout and the evidence's flags.
- Nothing in your files needs changing.

**Open questions** (each implemented as a reversible default):
1. *The initial count, 8.* Kept. The sizing rule moves away from it within a few rounds either way.
2. *The minimum count: the structural 2, or 3 with the reference rung.* Round 1's single scan reads low, since every replica starts from one state (Λ̂ 0.7–2.6 on the functional fixture, against the 11–17 it settles at). So the first sizing can drop to the minimum, then regrow by about a third each round (3 → 4 → 5 → 7 → 9 → 12 → 14). With R's 9,000-sweep burn-in there are 12 sizing rounds, and the dip costs nothing measurable. A short burn-in keeps a small ladder (7 sweeps: 4 and 6 rungs). I recommend keeping it. If short burn-ins matter, the cheap fix is a floor at `InitialLevels / 2`, or skipping the first round's shrink.
3. *The deadband: one rung,* as proposed. On the functional fixture it held 15 rungs through two rounds that asked for 14, then shrank when a round asked for 13.
4. *The rejection threshold, 0.9.* Under a Gaussian reading of the swap, a pair rejecting r has stepping-stone weights of relative variance e^{σ²/2} − 1, with r = 1 − 2Φ(−σ/2): about 26 at 0.8, 220 at 0.9 and 2,200 at 0.95, against a rung's few hundred to few thousand effective draws.
   - Measured: 0.955 (fixed 16) still read well, while 1.000 (fixed 8, the cap of 4) did not.
   - 0.9 flags before the weights outrun the draws. It flagged one self-sized pair at 0.902 whose estimate agreed, so it is a "look at the chains" flag, not a verdict.
5. *A new rung's source: its nearest rung on the barrier, then in β, then the colder.* That is the state a swap would hand it most readily. Ends map to themselves.
6. *Which naming rung keeps a replica: the one nearest its own β, the colder on a tie.* It moves least, and the ends stay kept without a special case.
7. *For S10's mixing work, not this chip's.* On the functional fixture the round trips are 0.0001–0.0008 per scan, against Theorem 1's ideal of a few hundredths at the settled counts. The cold rungs hold one structure for long stretches: the local moves, not the ladder, bound the mixing there.
