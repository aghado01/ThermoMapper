# Rungs at once (from the rung-parallelism session, 2026-11-08)

For anyone running tempered BASS fits, timing them, or touching `ParallelTempering`, `BassRunOptions` or what a rung's step reaches. Report: `issues/regression/bass/rung-parallelism-report.md`.

- `ParallelTempering.MaxDegreeOfParallelism` (default 1) steps a scan's rungs at once, swaps after. Draws are bit for bit at every degree (stress-tested to 48 rungs on 48 workers; pinned tier holds with it on and off).
- BASS asks for it with `BassRunOptions.RungWorkers`: 1 by default (in order, as before), −1 for what the chains leave of the bound (`Split`), n for at most n. Nothing changes for a run that does not set it.
- On this machine under load, rungs at once were slower (BASS rung steps are 25–85 µs; a preempted below-normal worker holds the scan). A quiet-window run of `BassRungParallelismTests.WallTime_*` (BASS_QUALIFICATION=1) decides whether the default becomes −1.
- Finding for anyone timing fits: workstation GC pauses were 42–49% of a 4-chain tempered fit's wall; `DOTNET_gcServer=1` took it from 2.46 s to 0.93 s. Chain-level scaling is GC-bound today.
- Anything shared that a rung's step reaches must be read-only or a pure cache whose fill cannot tear. `BassStructurePrior`'s count-law table is now precomputed for that reason. If you add a lazily filled cache to the model, the target, a block or a design snapshot, keep it per rung or make its fill atomic.
- BARS's ladders stay in order: its `SplineTarget` adds every rung's evaluation to one `Tally` (plain `++`).
- `BassSmcRun` now splits its bound through `BassRunOptions.Split` (81cf586, same arithmetic); `BayesianRegression.md`'s fit paragraph names `RungWorkers` (451cd49).
