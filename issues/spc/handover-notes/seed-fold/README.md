# The SPC seed fold takes T's bits (2026-11-01)

Commits: `2d817d3` (fix, test, docs), `cd10044` (primer note). Full record: `issues/spc-BARS/stopping-theory-primer.md`, §13.14.

What other sessions should know:
- **Every seeded SPC stream moved, once.** Anything reading a seeded SPC number from before `2d817d3` was on the old fold. A harness that rebuilds partway through a sweep puts its rows on two folds.
- **`SpcExecutionOptions.TemperatureQuantizationFactor` and `Derive`'s `quantizationFactor` are gone.** Nothing set them.
- **A T reached by a different float path is now a different temperature,** with its own stream. No caller relied on that being otherwise.
- **The toy's lineage parity test is a coin flip:** 7 of 12 seeds on the old fold, 6 of 12 on the new, and seed 42 passes on both. §13.15 counts 5 of 12 on its run.
- **Checkpoint names now carry T's round-trip form (`4cd7b13`).** The 5-decimal names gave 117 of the 289 points on a [1e-6, 0.3] log grid a shared file: silently wrong frames on one worker, an `IOException` on the default workers. A checkpoint directory written before `4cd7b13` reruns rather than resumes, and the extract command's glob would read both forms, so delete old directories.
- **The Welch readings hold on one graph.** It is named as resolved in §13.14 and in `results/graphs-new.txt`. "BlattHierarchy, Cauchy K = 5" on today's defaults builds a different graph.
- **A harness note: `Copy-Item` keeps the source's modified time.** Copying an edit made during a build into a tree that build read can leave the edit uncompiled, with 0 errors. Touch copied sources.

What the cold-grid probe was, and wasn't:
- BlattHierarchy at the generator's default N = 6000, not §12.28's N = 1200. That cost is why it ran for hours.
- The Cauchy K = 5 arm finished under both folds. The Gaussian K = 10 arm was killed by the reboot, so it has no data.
- The probe wrote each file only when an arm ended. A rerun should append per row, per the shared-machine memory.

Contents:
- `scripts/`: `SbProbe.cs` (the probe, a Fact that lived only in an exported tree), `patch_tree.py` (the env switch back to the old fold, plus per-test readings), `run_parity.sh`, `summarize_parity.py`, `sweep_table.py`, `ranksum.py`, `analyze_se.py`, `analyze_shared.py`, `derive_pin.py` (the independent evaluation behind the test's pinned literals). Paths in them are this machine's, and they point at `artifacts/sb`, which has since been deleted.
- `results/probe/`: per-replica S2 on the 289-point schedule, 12 sweeps × 8 replicas per fold.
- `results/readings-{old,new}-{0..11}.tsv`: per-test parity readings by fold and seed offset, on the tree at `63589b3` (before the compiler revision).
- `results/final-vizcore.trx`, `final-spc.trx`: the clean run on `72763fb` plus this change, the unpatched code that landed.
- `results/graphs-old.txt`, `graphs-new.txt`: the graphs behind the readings, built at `63589b3` and at `e1eea5d`. The hashes cover the arrays, so equal hashes mean the same graph; `graphs-new.txt` carries `Resolved`.
- `results/checkpoint-*`: the checkpoint probe and final runs on `e1eea5d` plus the fix. The two timeouts in `checkpoint-final-vizcore.trx` came from a 56-minute run on a saturated machine; both pass alone in `checkpoint-timeouts-rerun.trx`.
- `scripts/GraphRead.cs`, `CpProbe.cs`: the Facts behind those, which lived only in exported trees.
