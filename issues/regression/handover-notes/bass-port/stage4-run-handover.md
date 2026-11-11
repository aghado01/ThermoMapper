# Stage-4 run (S7, chip A): hand-over for S8 and the main session (2026-10-28)

From the chip running `issues/regression/bass/stage4-run-brief.md`. Full report: that brief, "Reports". Commits: `6609985` (prior sampler), `b5bca5c` (run), `87b23bc` (diagnostics), `759136d` (checkpoints), `b14c84d` and `7ac76ae` (review pass 1), `a402b14` (pass 2).

- **Reference rung.** `BassPriorSampler : IReferenceSampler<BassState>` (likelihood path) draws exactly from the prior on F and evaluates under `model.Target`. It keeps its own design, so a ladder needs one per replica.
- **Restart state needs the design's snapshot, not a reload.**
  - A strict reload of a structure differs in the last bits from a chain's incrementally updated design after deaths and changes. On a direct functional chain: 10⁻¹⁶ in ‖r‖², 2·10⁻¹⁵ in a coefficient.
  - That breaks bit-identical resume; a mutation check confirms it.
  - The pieces are internal and ready for per-replica use:
    - `BassDesign.Snapshot()` / `Restore(...)`, which copy the factor via `BorderedCholesky.Copy()` and keep β and the residual as held;
    - `BassStructureMoves.Snapshot(state)` and its restore constructor;
    - `BassChain.Snapshot()` / `BassChain.Restore(model, snapshot, target?)`.
- **Where a chain becomes a ladder:** `BassEnsemble.StartChain` (`BassRun.cs`).
  - Each chain seed fans out with `SeedTree.Derive(chainSeed, 3)` into the start, the moves and τ's block. The readout seed drives composition.
  - SeedTree's children are prefix-stable, so a fourth child for a sub-ladder leaves the first three unchanged.
- **Diagnostics** count the sampling phase only: `BassStructureMoves.ResetTallies()` runs at the end of burn-in. `BassSupportTally` counts each candidate term by its orders.
- **`BassFit`** has `ResidualDependence` (null) and no `Evidence` yet; S8 adds `BassEvidence? Evidence` between `Diagnostics` and `ResidualDependence`.
