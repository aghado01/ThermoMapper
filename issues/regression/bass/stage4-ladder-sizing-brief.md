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
