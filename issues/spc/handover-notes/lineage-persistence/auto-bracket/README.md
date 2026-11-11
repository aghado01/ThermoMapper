# `auto`'s bracket under narrow bandwidths (2026-11-02)

Follows primer §13.18, "Past the brief: `auto` under a narrow bandwidth". Exploration, one pass, read against what's written here before it ran. No default changes: any change to the estimator goes to the user.

## Written before any sweep ran

Build: the same export of `aeebe0d` as the settings probe, so its results carry over. Graphs: its 63, as `settings-probe/results/settings-resolved.tsv` reports them.

**A correction to what §13.18 compared against.** Its table set `auto`'s top against T_fs and T_ps read on the T_ps lens. T_fs is the global argmax of χ on whatever grid it's given, and χ carries a factor of 1/T, so on a wide grid it leans cold, and on a narrow one it can sit at the grid's edge (it does on Landsat K = 5 MAD, the toy MAD and Landsat K = 10 defaults). So those ratios may not measure a transition. They get rechecked here, and this may overturn "12–23 times too cold".

**What's measured on a wide lens.** [10⁻⁷, 4] with 180 log levels (about 1.10 apart), every graph, 6 seeds (42, 1–5), 1 × (200, 400).
- T_fs and T_ps from the detector, as read on that lens.
- T_melt: the hottest level with median S₂·N ≥ 2, where a node's cluster holds on average at least one other. No χ, no 1/T, no reference to q. This is the primary target for the hot end.
- Peaks of T·χ (N·Var m, without the 1/T): local maxima of the 3-point smoothed median at 25% or more of its largest. The T·χ peak list is what the multi-stage graphs get read by.

**The statistics, all from the graph alone.** s_i is node strength, s̄ its mean, g the gauge.
- A: today's. median(J)·s̄/ln q, over g. Top 5×, bottom 0.05×.
- B: the strength scale. top = s̄/ln(1 + √q), bottom = top/80. Primer §13.17: T_ps's prefactor is e^(−1/2)/4, so on the recipe this returns BWD's bracket with no new constant.
- C: A with an upper quantile: p90(J)·s̄/ln q, over g. Same multipliers as A.
- D: top = p90(s)/ln(1 + √q), bottom = top/80. Primer §13.17: local melts scale with local strength.
- E: the percolation top (§13.13, C43): the T where the lower envelope Bernoulli(p̃)'s S₂·N falls to 1.13, the value the recipe reads at its top. Bottom = top/80.
- F: the connecting coupling. J_b, the coupling at which the strongest edges first join half the nodes. top = 4·J_b/ln(q − 1) is a guess at a constant, so F is read as a predictor first, and as a bracket only if its ratio to T_melt is stable.

**The explanation under test.** §13.18 says the median J falls far below the couplings that hold the transition.
- It holds if, over the 33 graphs on edge statistics, ln(T_melt/A's top) tracks ln(p90 J/median J) with Spearman ≥ 0.7, and C's spread of ln(T_melt/top) is at most a third of A's.
- Otherwise it's reported as not holding, whichever candidate wins.

**A candidate passes if all of these hold.**
1. *Equivariant.* J → gJ moves it to T → gT, checked on every base and no-gauge pair to 10⁻⁹. Report whether it needs the gauge's declaration to do so.
2. *The hot end.* T_melt ≤ top ≤ 10·T_melt on at least 60 of the 63 graphs. Failures are listed.
3. *The readings.* 12 seeds at 1 × (200, 400), 48 levels, T_fs window, largest, 0.10, on the 33 graphs that use edge statistics.
   - On Landsat, MAD at K = 5 and 10 and compiler defaults at K = 5 pass on 9 or more of 12 (they pass on 12 on the T_ps lens, and on 0 on `auto`).
   - No graph falls against today's `auto` by §13.18's rule (the modal reading changes, its seeds fall by 3 or more of 12, a median ARI or purity falls by 0.05 or more), except a move between two planted levels, which is reported.
4. *Cost.* Seconds on Landsat, from the graph.

Only candidates passing 1 and 2 get step 3's sweeps, at most two. Ties go to the simplest, in the order B, C, D, F, E. If none passes, that's the result.

**What would change the framing.** If the wide lens shows `auto` does contain T_melt on Landsat's MAD graphs, then "too cold" is wrong, and the reading's failure lies elsewhere: T_fs on a too-narrow window, say. That gets reported as the finding.

**Not varied:** the 48 levels, budgets above the low one, distance, LMP, the floor, generator seeds, the T_ps path for graphs in BWD's units (read for comparison, not proposed for change).

## Results (primer §13.19, C60)

No candidate passed, so nothing goes to the estimator.
- Criterion 1: A, B, C, D, F to about 10⁻¹⁵. E at 2.4%, sampling noise. A and C need the gauge's declaration.
- Criterion 2, of 63: A 49, B 59, C 57, D 63, E 60, F 7.
- Criterion 3, D: Landsat 8, 0, 0 of 12, against the 9 required, and falls on the toy and iris under MAD. B was swept past the criteria and fails too.
- The explanation: Spearman 0.93, C's spread 1.22 against A's 1.99, so it doesn't hold as written. A's top over B's is 5·(median J/g)·ln(1 + √q)/ln q exactly, which is the reason.
- The framing: `auto` doesn't contain T_melt on Landsat's MAD graphs, so "too cold" stands at the hot end. "12–23 times" doesn't: `auto` holds the first break, and T_fs on the T_ps lens wasn't measuring a transition.

## Corrections to what's written above

- χ doesn't carry a 1/T here. The detector reads the reduced FK sum, ⟨Σ|c|² − |c_max|²⟩/N, which goes to 1 when hot. So the reason given for T_fs leaning cold is wrong: it sits cold because the first break has the largest swing. T·χ just runs to T, so the peak list above is useless. Stages are read on χ instead (`auto-table.py`), with a rule chosen after seeing the curves.

## Past what's written above

- B swept although it missed criterion 2 by one graph.
- The scan of window starts (`windowscan`, `run-auto-scan.sh`), added after D and B failed criterion 3. It re-reads cached sweeps, with no new sampling.

## Files

- `harness/`: the probe with `stats`, `profile` and `windowscan`, the run scripts and every reader.
- `results/`: `auto-stats.tsv` (statistics and brackets, 63 graphs), `auto-wide-table.tsv` (the wide lens per graph), `auto-table33.md`, `readings-against-auto.txt`, `window-starts.txt`, and the raw `auto-wide/`, `auto-out/`, `auto-scan/`.
- The sweep caches (3 GB here, 7 GB in the settings probe) stayed in the scratchpad and regenerate from `run-auto-cand.sh`.
