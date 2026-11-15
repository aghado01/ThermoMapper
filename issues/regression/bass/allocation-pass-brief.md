# BASS's allocation pass (chip brief)

Date: 2026-11-15. Status: **a draft for the user's approval.** The decisions below carry recommendations, and approving the brief takes them. It was drafted by the main session at the user's request ("Draft the brief"), after a survey of what was left for BASS's efficiency. The pass's report goes under "Reports" below, newest first. The chip leaves `docs/` alone; each report ends with **Doc changes owed**.

## Why

The rung-parallelism work ([its report](rung-parallelism-report.md), 2026-11-08) found that garbage collection caps BASS's wall time, and that the cap grows with the threads that allocate. The fixture was Friedman, p = 10, on a fixed ladder of 12 rungs with the reference rung.
- *One chain:* workstation GC's pauses held 13–17% of the wall.
- *Four chains in order:* the pauses held 42–49%. Each collection stops every thread, so four allocating chains pay for each other's garbage.
- *Server GC:* the same 4-chain fit took 0.93 s against 2.46 s, at the same load. These are indicative rows: no core was free.

The two readings differ, and both stand:
- The cost report's ceiling of 6.5% of a sweep was read on one untempered chain, at its own points.
- The rung report's 13–17% is one tempered chain of 12 rungs, and its 42–49% is four of them.

The sweep allocates per proposal. The efficiency pass ([its report](efficiency-pass-brief.md#the-efficiency-pass-built-chip-2026-11-04)) left allocation alone except where a lever touched it, and its lever 5 left τ's state and evaluation objects under a time threshold. The SMC kernel builds a reversible-jump chain on every move and takes a lock to add its account (`BassSmcKernel.Move`, `BassSmc.cs`). The cost profile priced that as a lever that keeps bits ([its levers' table](../performance-cost-regression-samplers-brief.md#the-levers-table)).

**What the pass is for.** A steady-state sweep and an SMC move whose allocation does not grow with n, every draw bit for bit. Then a reading that decides how BASS fits are hosted (decision 4). Server GC treats the symptom in a process setting. This pass treats the cause, so every host gains, including those that cannot change their GC.

## The decisions (recommended; approving the brief takes them)

1. **Bits unchanged, everywhere.** Allocation changes where a value lives, never its arithmetic. Every lever leaves every draw, evaluation and reported number bit for bit, and the checks below say so. A change that would move a bit is outside this pass.
2. **Pure scratch.** Reused memory is sized from the owner's dimensions and configuration: n, the rows, the locations, M_max, the factor's order. It is never sized from history: the rule the efficiency pass found when it dropped lever 1's hit-count gate ([its report](efficiency-pass-brief.md#the-efficiency-pass-built-chip-2026-11-04)).
   - It belongs to the object that runs the proposal: a design, a chain, or an SMC worker. It is never shared across threads.
   - A restored or handed-over object builds its scratch at first use and continues bit for bit.
3. **What escapes keeps its own memory.** An array reachable from a state, a draw, a snapshot, an evaluation a caller holds, or a returned record is allocated as today, since states keep their generative state.
   - Only memory whose life ends inside a proposal, a commit or a τ step is reused.
   - An accepted proposal's buffers become the design's by swap, not by copy. A rejected proposal's buffers return to the workspace.
4. **Hosting is decided after the pass, from its reading, by the user.** The pass changes no process setting: no `ServerGarbageCollector`, no `DOTNET_gcServer`, no `GCHeapCount`. Its final reading puts workstation and server GC side by side, before and after, so the decision is made on the remaining gap.
5. **Inventory first, then one lever at a time.** The levers below are candidates. The inventory decides which are built, by the rule fixed before it runs. Each lever is its own commit, with its tests, its checks and a fresh-eyes review. The main session verifies it from a clean archive before the next starts.

## The inventory (counts: no quiet window needed)

Bytes are the same under load, so the inventory runs at once. It reads:
- **Bytes per sweep, and per part,** at the efficiency pass's ten points ([the table](efficiency-pass-brief.md#the-points)), through the cost harness's bytes meter (`BassCost.Meter.Bytes`, `BASS_COST_POINTS`).
  - The parts: birth, death and change proposals, the commits, τ's block, the composition, and the other blocks.
  - Beside them, the n-scaling at the base point: n 200, 1,000, 10,000 and 50,000.
- **By site:** one allocation trace at the base point and at n 10,000, through `dotnet-trace` with the GC allocation events (`--profile gc-verbose`), or an in-process `EventListener` on the same events if `dotnet-trace` is not installed. It gives the types and the call stacks, as a table of the top sites: bytes per sweep, and whether each grows with n or with M.
- **GC counts:**
  - generation 0, 1 and 2 collections per 1,000 sweeps (`GC.CollectionCount`);
  - the pause share (`GC.GetTotalPauseDuration` over the run's wall);
  - on the rung fixture at one chain and at four.
- **SMC:** bytes per move and per stage on an SMC fit of the Friedman fixture.

**The rule, fixed before the inventory runs.** A site is taken if it is at least 5% of bytes per sweep at any point, or if its bytes grow with n or with M_max. A smaller site is named in the report and left.

**What it decides:** which levers below are built, and each built lever's target in bytes per sweep.

## The levers (candidates; the inventory decides)

### 1. The column's scratch

A proposal evaluates its column and passes over it to centre it, whiten it, stack it for a direct design and form the prior side's column. `BassDesign.cs` allocates n-length arrays for these on each evaluation (the `new double[_cases]` and `new double[_rows]` around the column's evaluation).
- **The change:** a per-design workspace of the buffers a proposal needs, sized by the rows and the locations. An accepted proposal's column moves into the design by swap; a rejected proposal's buffers stay for the next.
- **Bits:** identical. The same operations, in the same order, write to other memory.

### 2. The pending's vectors and factors

The border, b, the diagonal, the bordered factor's row and the candidate pencil's arrays are sized by M.
- **The change:** a workspace sized by M_max + 1, swapped on acceptance as in lever 1.
- **Bits:** identical.

### 3. τ's step

Each slice evaluation builds a state (`WithLogTau`) and an evaluation object. The efficiency pass's lever 5 left both under its 25% time threshold; here they are judged by bytes.
- **The change:** if the inventory takes them, an evaluation path that writes into reused storage. The draw's own state, which escapes, stays a fresh object (decision 3).
- **Bits:** identical: the marginal's arithmetic and its order are untouched.

### 4. The SMC move

`BassSmcKernel.Move` constructs a `ReversibleJumpChain<BassState>` per move and adds its account under a lock.
- **The change:**
  - one chain per worker, reset per move to (state, rng, temperature);
  - one account per worker, merged between stages in worker order.
- **Bits:** identical draws. Integer tallies merge to today's totals in any order. A floating-point field, if the account has one, merges in a fixed order, and the report names it.

### 5. The shared reversible-jump chain

`ReversibleJumpChain.Step` (`src/maths/samplers/rjmcmc/`) evaluates per proposal, sums `TotalWeight` twice, and looks up the reverse move by string. The cost profile priced a precomputed reverse index and a single pass as levers that keep bits.
- **Shared code:** BARS and the other engines run this chain too, so this lever is its own commit, and their tests run with it.
- **Bits:** identical.

### 6. The ladder's scan

Per-scan allocations in the swap proposals and the fan-out's closure, built only if the inventory takes them.
- **Bits:** identical, and rungs at once stay bit for bit with rungs in order.

## Checks common to every lever

- **The eight `BassBitIdentityTests` hashes,** after every lever and with hardware intrinsics off (`DOTNET_EnableHWIntrinsic=0`).
- **Resume and hand-over:** `BassRunResumeTests`, `BassTemperingResumeTests` and the cost harness's `Untouched` check.
- **Parallel, bit for bit:** `BassRungParallelismTests` (rungs at once against in order, and the resumes across degrees), chains in parallel, and SMC particles across worker counts.
- **Aliasing.** A test holds a draw, a state and an evaluation across many further sweeps and moves, and finds every array it reads unchanged. Reused scratch must never reach anything that escapes.
- **An allocation gate, new.** Steady-state bytes per sweep at the base point after burn-in sit under a bound taken from the inventory's after-value. The built levers' bytes hold at n 1,000 and n 10,000 alike. Bytes do not depend on load, so this gate runs in the ordinary suite.
- **Fresh-eyes review at the end of each lever,** fixed before the lever is reported: correctness (a weight form, a direct design, a restored design, an unaccelerated path, a path the chosen cases miss) and coverage (a design form the tests do not drive).
- **The BASS fixtures** pass in Release after every lever: `BassDesignTests`, `BassBasisTests`, `BassDesignKernelTests`, `BassMovingResponseTests` and the others. A failure that predates the change is named, not fixed.

## The reading (timed)

- **Fixture:** the rung report's Friedman fixture at n 500 and 4,000.
  - Arms: one chain and four chains in order, each under workstation and server GC, before and after the pass.
  - Beside them, µs per sweep at the base point.
- **Protocol.** Timed sections never run at below-normal priority.
  - Build at below-normal priority, and start the test host at normal.
  - The timed section alone runs at process High, thread Highest, `ProcessorAffinity` 0xFF.
  - Record the machine's busy cores over each block (`GetSystemTimes`), and warm until tiering settles.
- **Pairing.** The baseline is a clean `git archive` build of the brief's commit, frozen and never rebuilt. It is paired ABBA with each lever's build. The draws are bit-identical, so each block does the same work in either binary.
- **Counts beside times:** bytes per sweep, collections per 1,000 sweeps, and the pause share.
- **Without a quiet window:** paired, pinned blocks across the two builds, labelled with their load. The counts stand either way; a timed row under load is marked indicative.

**What the reading decides:**
- **Hosting (decision 4).** If server GC still takes 25% or more off the 4-chain wall after the pass, the report recommends it for the hosts that run fits. It states the memory and heap-count cost of server GC on a shared machine beside the saving, and the user decides. If the gap falls under 10%, hosting stays as it is, and the report says so. Between the two, the report gives the rows and no recommendation.
- **Not decided here:** `RungWorkers`'s default still needs free cores. The pass reports the rungs-at-once arms, labelled, for the record only.

## Sequencing and coordination

- **Other chips in BASS's files, as of 2026-11-15:**
  - the BPPR chip (Regression expansion's): `BassInputTransform.cs`, `BassInputs.cs`, `BassTerm.cs`, `DesignKernel.cs`, `DesignKernels.cs`, `DirectionalHinge*.cs` and `sensitivity/`;
  - E8's slice 8: `sensitivity/`, and possibly `DesignKernels.cs` and `BassBasis.cs`;
  - the conformal chip: `BassCandidates.cs`, `BassPointwise.cs` and `BassFit.GroupConformal`.
- **This pass's files:** `BassDesign.cs`, `BassSweep.cs`, `BassTarget.cs`, `BassState.cs`, `BassSmc.cs`, `BassLadder.cs`, `src/maths/samplers/rjmcmc/ReversibleJumpChain.cs` and `src/maths/samplers/smc/SmcRun.cs`, with the tests it adds and the cost harness's meters.
- **Before each lever,** run `git status` and `git log` on its files. A file another chip holds is agreed through the main session first, and through Regression expansion for BPPR's.
- **Heads-up:** the defects chip's slice 9 reroutes `LogBeta` through `LogGammaRatio`, which may move pinned values that read Beta laws in their last bits. A pinned hash that moves during the pass is checked against that commit before it is read as this pass's.

## Not in this brief

Each is its own question, for later:
- prefix tables on tempered and SMC chains, and lever 1's skip until acceptance;
- the efficiency pass's other unbuilt tails;
- the cost profile's levers for BARS, BAPS, BGP and the changepoint model, and its paused timed grid;
- any process GC setting (decision 4);
- `RungWorkers`'s default.

## Working rules

- **`main`,** no worktree, no subagent, nothing pushed. Commit by explicit path (`git commit -- paths`) and leave nothing staged: the index is shared. Commit messages describe the work and carry no attribution lines.
- **Build and test** with `-p:ArtifactsRoot=D:\aghado01\ThermoMapper\artifacts\chip-allocation-pass\` inline on each `dotnet` command. Delete that folder at the end.
- **Machine:** shared.
  - Builds at below-normal priority with `-m:2`, tests with `-- xUnit.MaxParallelThreads=4`.
  - Timed sections as in the protocol above.
  - Clear `BASS_QUALIFICATION` and `BASS_COST_CLOCK` after a gated run.
- **Other sessions' uncommitted work stays untouched.**
- **Line endings:** briefs and C# are LF.

## Report

Under "Reports" below, newest first:
- **the inventory:** bytes per sweep and per part at the ten points, the n-scaling, the site table, GC counts and pause shares, SMC's bytes per move, and the rule's verdict per site;
- **per lever:** what changed, the checks with their numbers, bytes per sweep before and after, and what was left unbuilt and why;
- **the reading:** times, counts and pause shares, before and after, under both GCs, with the machine's load per block, and the hosting recommendation or its absence by the rule above;
- **deviations** from this brief, in the order met;
- **Doc changes owed.**

## Reports

Newest first.
