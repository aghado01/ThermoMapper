# Stage-4 ladder sizing (chip E): hand-over for S10 and the main session (2026-10-28)

From the chip running `issues/regression/bass/stage4-ladder-sizing-brief.md`. Full report: that brief, "Reports". Commits: `ec565f6` (task 1, shared), `d64ee4e` (task 2, BASS), `66c66be` + `abf16cf` (review pass 1), `36ce26a` (pass 2).

- **Entry points.** `TemperingSpec.Levels` null (default) self-sizes within `MaxLevels` (64) from 8; `HotEnd` null (default) auto-picks, resolved in `BassResolved.HotEnd` and logged. `BassLadder.Resize`. Readout: per-chain `PairRejection`/`PairAcceptance`, `Levels`, `MaxPairRejection`, `AtMaxLevels`. Evidence: `WithinChainError`, `BetweenChainError`, `BetweenToWithin`, `MaxPairRejection`, `Saturated` (> 0.9), `AtMaxLevels`.
- **Shared.** `ParallelTempering.Resize`, `RoundRejection`, `PairScans`, `TraceRungs(traceSlopes)`; `TemperingLadder.Place`, `Sources`, `Size`; `LadderEvidence.Pool`. A ladder that never resizes is bit-identical (pinned). BARS's evidence error is now the larger of within and between.
- **For S10's cost baseline.** Functional fixture (N = 175, M_max 3): settled 11–17 rungs, mean M 3.0; burn-in over-reads the barrier (5–7.5 against 3.8 in sampling), so counts settle above what sampling needs. N = 18: 4–6; fixture B (N = 32): 6–7; fixture C (N = 27): 5–6.
- **S8's 180–186 cross-check is superseded:** at 4,000 + 4,000 sweeps, fixed 32/64 rungs read 193.0–193.4 and the self-sized forms 192.6 and 193.4.
- **Docs owed** (main session's `docs/Maths`): the new shared members, and S8's paragraph with the self-sizing ladder.
