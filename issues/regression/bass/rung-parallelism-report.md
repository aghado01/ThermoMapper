# Rungs at once: a tempering ladder steps its rungs in parallel, every draw bit for bit (report, 2026-11-08)

Taken up from the lineage review's *parallel rungs* item ([lineage-review.md](lineage-review.md), "Taken up on 2026-11-08"). `ParallelTempering.Step` stepped a ladder's rungs one after another, and `BassEnsemble` ran chains in parallel only, so a tempered fit used at most one core per chain. A rung's step can now run `MaxDegreeOfParallelism` at a time, and a BASS run asks for that with `BassRunOptions.RungWorkers`, within what its chains leave of the bound. The draws are unchanged bit for bit at every setting. The default stays in order, because on the shared machine under load rungs at once were slower (see the wall time below).

## What landed

On `main`, nothing pushed: `454f2ef` the count-law table; `4776ff5` the ladder's option and its sampler tests; `efd32a1` BASS's run option, its wiring, tests and docs; this report after them (`0d72a37`); then `451cd49` the fit's paragraph in `BayesianRegression.md`, and `81cf586` SMC's split through `Split` (regression commit scope 1,615 of 1,615).

- **The ladder's option.** `ParallelTempering.MaxDegreeOfParallelism`, 1 by default: in order, on the caller's thread.
  - Above 1, the scan's rung steps (moves, blocks, or the reference rung's exact draw) fan out through `Parallel.For`. The swaps, the round-trip bookkeeping and the traces run after every rung has stepped, as before.
  - A rung's exception is held by its rung, and the coldest one is rethrown as the in-order loop throws it, not as the fan-out's aggregate.
  - The degree is not restart state: `LadderSnapshot` is unchanged, a restored ladder steps in order until set, and a resize keeps it.
- **BASS's run option, nested in its one bound.**
  - `BassRunOptions.Split(chains)` gives min(chains, P) chains at once and ⌊P / that⌋ workers left for each chain, at least 1. P is the bound, or the processor count when unbounded. This is the rule SMC already splits replicates and particles by.
  - `RungWorkers` (1 by default; −1 takes all the chains leave; n caps it at n) sets each ladder's degree through `RungWorkersFor(chains)`. `BassEnsemble` applies it at start and on restore, and the flat histogram's pilot does the same for its ladders.
  - `Bass.Fit` refuses a value that is neither −1 nor at least 1. All of this is a pure function of the options and the chain count, never of load.
- **The count-law table precomputed.** `BassStructurePrior.LogCollectionWeight` read a 1,024-entry table filled on first read (NaN until then), shared by every chain and rung of a fit. It is now filled when the prior is made, by the same expression, so nothing writes to it after.
- **Docs.** `docs/Maths/Samplers.md` (tempering: rungs at once; the API row; the configuration table) and `BassLadder`'s, `BassEnsemble`'s and the pilot's docstrings.

`BassSmcRun` kept its own inline split (the same arithmetic) while the adaptive-mutation session worked beside that file; it calls `Split` since `81cf586`.

## The audit: what a rung's step reaches

The model-wide state below was already shared by `BassEnsemble`'s parallel chains before this change. Stepping rungs at once adds sharing only within a ladder: the start state every replica begins on, the blocks, and what a hand-over passes between replicas' designs.

| State a rung's step reaches | Who shares it | Finding | Done |
|---|---|---|---|
| The replica's chain: state, cached evaluation, move stream, move counters | its rung | private | nothing |
| The replica's block stream and block counters (`ParallelTempering`) | its rung | private, one array per replica | nothing |
| The blocks: `MoveStep`, `BassTauBlock`, `BassSigma2RateBlock` (new per replica); the noise law's coefficient, latent, scale, shape and interweaving blocks (one instance per model, shared by every replica and chain) | every rung | stateless: readonly references to the model | nothing |
| `BassStructureMoves`: its design, counts, pending candidate, the column-means and sweep-column caches, its tallies | its rung | built once per replica | nothing |
| `BassDesign`'s columns, factor, β, residual; its prefix tables and hinge orders | its rung, columns shared after a hand-over (`DesignFor` shares them) | the columns are written once, when evaluated, and only read after; a ladder's designs build no prefix tables | nothing |
| `RidgeSpectrum` (and its bordered views), `PriorPencil` | a design and its snapshots: replicas after a hand-over | immutable after construction | nothing |
| `BassState` and its fit's statistics; the start state and term array every replica begins on | every rung | immutable | nothing |
| `BassEvaluation`'s untempered marginal, a `double?` filled on read (not an atomic store) | its rung | each chain evaluates its own start; a swap exchanges evaluations, never copies one; a reference draw and a resize's fork are evaluated afresh | nothing (per rung by construction) |
| **`BassStructurePrior`'s count-law table**, 1,024 NaN-initialized doubles filled on read | every chain and rung | a pure function of M, so a race wrote equal values; its safety rested on an aligned 8-byte store not tearing | **precomputed at construction** (size by constant, never by hit history) |
| `BassTarget`'s `Sigma2PriorTerms`, 8 slots keyed by the observation count | every chain and rung | immutable entries published by reference: a racing write only evicts, and every value equals recomputation | nothing |
| `WorkingResponse`'s ‖y − ȳ‖²_W, a NaN-sentinel double filled on read | every chain and rung (the model's), or a state's own | a pure value; its fill is one aligned 8-byte store on the 64-bit runtime | nothing (its own comment already says two readers compute the one value) |
| `StudentTMixing`'s skewed law, made on first use (`??=`) | the model's | two racing readers may each make one; `GigMixing` is immutable, so either serves | nothing (GBASS's file) |
| `GaussLegendre`'s rules | the process | `ConcurrentDictionary` | nothing |
| `BassPriorSampler` (the reference rung): its design, last draw, counts | the hot rung during the step; the hand-over reads it after | private to its rung during the step | nothing |
| `ParallelTempering`'s own: `Scans` (read by each rung's blocks), the swap stream, statistics, labels, traces | the ladder | `Scans` is only read during the step; the rest is touched after it | nothing |

**BARS** (`BarsEnsemble`'s ladders). One `SplineTarget` per ladder adds every evaluation to the chain's `BarsReadout.Evaluations`, a `Tally<DesignRejection>` whose `Add` is a plain `++`: rungs at once would lose counts. Every replica's moves wrap the chain's one kernel; `LocalBetaKernel` is immutable and `AdaptiveLocalBetaKernel` changes only in `Adapt`, between steps. BARS's ladders stay in order: `BarsEnsemble.Run` takes no bound, so nothing sets their degree, and the pinned BARS hash holds. See the open items.

**BGP** has no ladder (`GpRegression`, `BgpSampler`). **Calibration** (`CalibrationSampler`) gives each replica its own walk and blocks; its target and its adaptation lock were not audited, and its ladder stays in order.

## Gates

Every build was from a `git archive` snapshot plus this session's files, with `ArtifactsRoot` set under `artifacts/chip-rung-parallel/`, below normal, with MSBuild node reuse off.

| Gate | What ran | Result |
|---|---|---|
| Pinned tier, option on | `test-select -Scope Commit`, Maths.Regression.Tests, at 8b8c063. The split was then the default, so every tempered pinned fit stepped its rungs at once: 2-chain fixtures on 9 workers, 3-chain on 6, the flat histogram's pilots on 9. | 1,799 of 1,799, the 31 pinned among them |
| Pinned tier, option off | `Tier=Pinned` on a build whose ladder ignores the degree | 31 of 31 |
| Commit scope, final | `test-select -Scope Commit` at a670046 with the final files (rungs in order by default) | Maths.Regression.Tests 1,800 of 1,800; Maths.Samplers.Tests 174 of 174 |
| On against off, per fixture | `BassRungParallelismTests`: the self-sizing baseline's hash (`0d067c924599ed93` both ways, its pinned value), the seven fixed-level S8 pins, four flat-histogram pilot seeds, each in order and at once at one bound | all equal |
| Stress | `ParallelTemperingRungParallelismTests`: 32 rungs resized to 40, on 2, 7, 32 and 64 workers, every schedule, blocks, a reference rung, adaptation, traces. `BassRungParallelismTests`: 24-rung ladders on 16 workers over the integrated, linear and powered forms, the functional and categorical bases, the g-form held, the ridge, Student t, and a skewed noise with interweaving; a self-sizing ladder. Each compares every scan's cold state and every replica's final state, density, moves, blocks, refusals and traces. | all bit for bit |
| Stress at scale | the same at 48 rungs on 48 workers, 1,023 + 2,000 scans, behind `BASS_QUALIFICATION` = 1: reference, ridge, functional, skewed noise with interweaving | 4 of 4 bit for bit (3,167 entries each) |
| Resume mid-ladder | a self-sizing ladder snapshotted mid-burn-in at 8 workers and restored at 3; a fit checkpointed with rungs in order and resumed at once | bit for bit with the uninterrupted run |
| Exceptions | two rungs throwing in one scan, at 1, 4 and 32 workers | the coldest rung's own exception, at every degree |

## Wall time (indicative: no free cores)

`BassRungParallelismTests.WallTime_OneAndFourChainsOfTwelveRungs` behind `BASS_QUALIFICATION` = 1: the Friedman fixture (p = 10) on a fixed ladder of 12 rungs with the reference rung, each arm three times in turn after a warm-up, the median kept. Release, .NET 10, 18 logical processors, below-normal priority. No quiet window came: other processes held 6–7 cores throughout and the machine read 96–100% busy, so this process got about 1.7 cores on average. The single-threaded in-order arm at n = 4,000 got 0.82 of a core. These rows say what rungs at once do on the machine as it is shared. They do not say what they do with free cores. The rows are in `artifacts/chip-rung-parallel/cost-4/` (workstation GC) and `cost-servergc/`.

| Arm (n = 500, 1,000 + 1,000 sweeps) | Workstation GC: median s | CPU / wall | GC pause | Server GC: median s | GC pause |
|---|---|---|---|---|---|
| 1 chain, rungs in order | 0.97 | 1.27 | 13.5% | 0.60 | 4.4% |
| 1 chain, 2 rung workers | 1.32 | 1.50 | 15.2% | 0.92 | 7.8% |
| 1 chain, 4 rung workers | 1.48 | 1.29 | 16.8% | 1.03 | 13.4% |
| 1 chain, 6 rung workers | 1.33 | 1.19 | 15.8% | 0.95 | 7.9% |
| 1 chain, 12 rung workers | 1.45 | 1.08 | 16.9% | 0.99 | 13.2% |
| 4 chains, rungs in order | 2.46 | 1.56 | 49.3% | 0.93 | 35.6% |
| 4 chains, 2 rung workers each | 2.90 | 1.82 | 44.1% | 1.54 | 38.0% |
| 4 chains, 4 rung workers each | 2.89 | 2.20 | 45.8% | 1.85 | 38.8% |

| Arm (n = 4,000, 300 + 300 sweeps, workstation GC) | Median s | CPU / wall | GC pause |
|---|---|---|---|
| 1 chain, rungs in order | 0.61 | 0.82 | 13.2% |
| 1 chain, 4 rung workers | 0.71 | 1.46 | 14.6% |
| 1 chain, 12 rung workers | 0.71 | 1.83 | 18.2% |
| 4 chains, rungs in order | 1.59 | 2.02 | 41.8% |
| 4 chains, 4 rung workers each | 1.97 | 2.41 | 43.8% |

The 48-rung stress run read the same way: 3.0 s in order against 3.6 s at once (reference), 7.0 against 7.2 (functional), 5.4 against 3.7 (ridge), and 16.8 against 11.9 (skewed noise with interweaving).

**Reading.**

- **Rungs at once lost here.** A BASS rung step is short: about 25 µs at n = 500 and 85 µs at n = 4,000 on this fixture's structures. Each scan fans out and waits for its last rung. With the cores held elsewhere, a below-normal worker that is preempted holds the whole scan, so the wait costs more than the parallel work saves.
- **Hence the default.** That is why `RungWorkers` defaults to 1, in order, rather than to what the chains leave of the bound.
- **What a free-core reading must settle.** Whether the fan-out pays when cores are free, and from what step length.
- **GC is the larger cap, and it caps the chains already.** Workstation GC's pauses held 42–49% of a 4-chain fit's wall and 13–17% of one chain's. The sweep allocates per proposal, and every GC stops every thread. Under server GC the same 4-chain fit in order took 0.93 s against 2.46 s, at the same load.

## Open, each with the decision it would change

- **`RungWorkers`'s default: a free-core reading.** The wall-time test above, run in a quiet window, decides it.
  - If rungs at once win there by a margin, `RungWorkers` defaults to −1, what the chains leave.
  - If they lose or tie, it stays at 1, and the option is for runs whose rung steps are long (large n, many terms) on a machine with cores to spare.
  - Either way the draws are the same.
- **GC caps the parallel scaling the engine already has.**
  - Workstation GC's stop-the-world pauses held 42–49% of a 4-chain tempered fit's wall. Server GC took the same fit from 2.46 s to 0.93 s under the same load.
  - Two remedies, neither touching a draw: server GC in the hosts that run fits (a process setting: the test projects, the REPL), or an allocation pass on the sweep (per-proposal arrays, states and evaluations).
  - The decision either changes is how BASS fits are hosted, and whether the next efficiency pass targets allocation.
- **BARS's rungs at once.**
  - Needs per-rung evaluation tallies (one per replica, merged at the readout, or an atomic add), and a bound on `BarsEnsemble.Run` to split.
  - Its rung steps are shorter still (a 1-D spline's design and factor), so the same free-core reading decides whether to wire it at all.
- **The fan-out's grain.** One `Parallel.For` per scan, waiting on its last rung.
  - If the free-core reading shows the wait dominating at BASS's step lengths, steps could be grouped to fewer, longer tasks.
  - The ladder could also step in order below a rung count set in the config.
  - Neither may be chosen from measured load or timing at run time.
