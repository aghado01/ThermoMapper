# LineagePersistence across the remaining graph settings (2026-11-02)

**Done:** the write-up is primer §13.18 (C59), scoped in `docs/Clustering/SpcAnalysis.md`. `harness/` holds the code, `results/` the resolved dump, the tables and the per-seed runs. The sweeps themselves (6.7 GB) stayed in the session's scratchpad.

Added after reading the results, so outside what's written below: Landsat K = 5 base, OR and q = 5 at 4 × (500, 1000) (`run-settings-b4.sh`), θ read on the q = 5 sweeps with no new sampling (`results/runs-theta/`), and the K = 20 graph's K̂.

Extends primer §13.16's probe (K and the coupling only) to the filter, the repair, the bandwidth rule, `ScaleEdges`, the gauge,
q and the stack's θ. Exploration, not validation: one pass, read against what's written here before it ran.

## Written before any sweep ran

Build: an export of `aeebe0d`, read through `LineagePersistence.Resolve`. 12 seeds (42, 1–11).

Base graphs, the recipe (mutual, MST, mean edge length after repair, gauge 1/K̂, no LMP, q = 20):
iris K = 5 and 10 Gaussian (w96, 1 × (200, 400) and 4 × (500, 1000)); toy K = 10 Gaussian (p48, low budget);
Landsat K = 5 and 10 Gaussian (p48, low budget); BlattHierarchy K = 5 Cauchy (w161, low budget).

One change at a time from each base: `or`, `norepair`, `mad`, `before`, `nogauge`, `q5`, `q40`. Then `def`, K and the kernel alone, so
everything else is the compiler's. θ = 0.3, 0.4, 0.6, 0.7 on the base sweeps, with no new sampling.
Grids: the dataset's lens at each graph's q. `nogauge` also on that lens × K̂, since the gauge is exactly a rescaling of T.
Graphs not in BWD's units also on their own `auto48`.

Each run reads: longest-lived on the full, T_fs and [T_fs, T_ps] windows; largest and midpoint; shares 0.10, and at the T_fs window the size
rule, 0.15, 0.25; `ColdCensored` on what's selected; a T\* ladder on iris and BlattHierarchy.

Expected before running:
- `nogauge` on the × K̂ lens matches `base` within seed noise, since it's the same model. If it doesn't, something's wrong with the harness or the seed fold.
- `nogauge` on the same lens loses the hot end, by a factor of K̂.
- Setosa's edge count and J move with `or`, `mad` and `before`, so "the bridge is shut" can stop holding. `norepair` removes it.
- q moves the bridge's opening only through ln(1 + q) and the floor 1/q, so T\* should move a little and the readings not much.
- θ sits between the floor and 1, and 0.5 isn't special to the resolver.

What counts as moved: the modal reading changes, its seeds fall by 3 or more of 12, or a median ARI or purity moves by 0.05 or more.

What would bring a default back to the user:
- *The largest level:* midpoint beats it (by 2 or more seeds passing, or 0.05 ARI or purity) on most graphs outside Landsat at K = 5, or on every arm of one dataset.
- *The 0.10 share:* another share beats it over all Landsat arms by 10% of seeds passing without losing median purity on most arms, or passes 9 of 12 on an arm where 0.10 passes under 6.
- *The defects:* any defect's fix making a reading worse than the old behaviour on some arm is reported, whichever way the defaults go.

Added after the graph dump and before any sweep: iris K = 10 Cauchy as a seventh base. On every Gaussian arm the one edge joining setosa stays far below the grid (J from 1.7·10⁻⁷ to 2·10⁻⁴ at best, and `norepair` has none), so those arms can't test that split. Before-sweep check: one base sweep at this build is byte-identical to the §13.16 probe's cache.

Budget follow-up, only if a largest-midpoint flip appears on Landsat at the low budget: that arm at 4 × (500, 1000), since §13.16 found the cost closes with sampling.

Not varied: distance, the ε-ball topology, LMP (bar what `def` resolves), the floor, `ScaleRank`, generator seeds, the kernel beyond the bases.
