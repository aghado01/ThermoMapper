# A lens spanning the graph, read per stage (2026-11-02)

Follows primer §13.19's "A lens spanning the graph, read per stage", which the user asked to have run. Exploration, one pass, read against what's written here before it ran. No default changes: anything for the estimator or the detector goes to the user.

## Written before any sweep ran

What I'd already looked at, all without labels: §13.19's wide-lens curves (χ and S₂·N, per seed), and each graph's coupling statistics on HEAD. Those fixed the cold end's margin and the shelf rule's constants below. The readings against labels are from §13.19 and earlier, on other lenses, so I know Landsat's upper χ stage passes and that the toy's stripes are cold.

**Build.** An export of HEAD, `394e3ef`, not §13.19's `aeebe0d`, since the coupling floor now keeps connectivity. 6 of the 63 graphs change: toy under MAD and compiler defaults, Landsat K = 5 and 10 under MAD and compiler defaults (`span-resolved.tsv`, against `settings-resolved.tsv`). `auto` is rerun on those 6 as the baseline, plus 2 unchanged graphs to check that they reproduce the old runs exactly; if they do, the other 57 baselines carry over.

**Graphs.** All 63, as `span-resolved.tsv` reports them. The 33 on edge statistics run first; the 30 on T_ps's bracket run after, as the "no cost" check.

**The lens, from the graph alone.** Rungs of 80^(1/47) = 1.098, the spacing of the 48-level lens, so the level count follows the span.
- Top: D's, p90(s)/ln(1 + √q). It held T_melt on 63 of 63.
- Bottom: J_first/(4·ln(1 + q)). J_first is the coupling of the coldest split a reading could keep: going down the edges by J, the last edge to join two components of m or more nodes, with m the reading's `MinClusterSize`. Edges below the floor, tested before the gauge, don't count: they're the graph's connectivity (repair and the floor's restored edges), and counting them drops the bottom by up to 33 decades. A bridge of coupling J opens at even odds at T = J/ln(1 + q). The 4 is a margin, set from the table, where the first break (χ ≥ m²/N on the wide lens) sits at 0.3 times that or above on connected graphs.
- If no such edge exists, the bottom is top/80.

**The stage rules, per seed, from that seed's own curve.** Both are in `harness/span-rules.py`.
- *χ's stages:* §13.19's rule, unchanged: local maxima of the 3-point smoothed χ, at 10% or more of its largest, above 2, standing 1.5 times above the deepest point to the next higher one on each side. On the wide lens it agrees between seeds on 31 of 33 graphs (the modal count on 5 or more of 6, positions within 2 times).
- *Shelves of S₂·N:* χ's height is in nodes, so a balanced hierarchy's later stages sit on shelves that never make a local maximum. BlattHierarchy's χ shows one stage where S₂·N shows 600, 200 and 50. So: g = −d ln(S₂·N)/d ln T on the 3-point smoothed curve, a central difference over ±4 rungs. A shelf is a minimum of g below the last melt (S₂·N ≥ 2) where g rises by 0.25 or more on both sides before reaching a lower value. Its cold and hot steps are the peaks of g on each side. Half-width 4 and depth 0.25 keep BlattHierarchy's three and agree on 23 of 33. Depth 0.5 agrees on 31 but loses them.
- A seed's stages are aligned to positions: the medians over the seeds with the modal count, each stage going to the nearest position within a factor of 2. A seed with no stage at a position counts as not passing there.

**The readings,** largest, share 0.10, θ = 0.5, m as before. 12 seeds (42, 1–11) at 1 × (200, 400), as in §13.19. K is a seed's count of stages.
- `chi up`, **the primary**: the whole lens, and [T_k, ∞) for each χ stage. K + 1 readings.
- `chi band`: [bottom, T_1], [T_1, T_2], ..., [T_K, ∞). `chi slice`: the whole lens, and the slice at each T_k.
- `shelf up`, `shelf band`, `shelf slice`: the whole lens, and for each shelf [cold step, ∞), [cold step, hot step], and the slice at its flat point.
- `null`: the whole lens, and [T, ∞) at K evenly spaced points in ln T across the lens, K the modal count. The same for shelves. This asks whether the stages do anything a blind scan wouldn't.
- `tfs`: today's detector on this lens, [T_fs, ∞) when it finds a cliff, else the whole lens. A single reading.
- `hot`: [T_K, ∞) for the hottest χ stage, else the whole lens. A single reading, and the one way this could end without handing anyone a list.
- `scan`: [T, ∞) at every third level, as a map, on the 33.

**A list "reads" a target** if at one position it does so on 9 or more of 12 seeds. "Holds `auto`" uses §13.18's rule against today's `auto` on the same seeds: the modal reading kept without falling by 3 or more of 12, and no median ARI or purity falling by 0.05 or more. A move between two planted levels (iris's S+merged and S+species, BlattHierarchy's coarse, medium and fine) is counted separately and reported. On Landsat's four MAD and compiler-default rows, a position reading the target isn't charged for ARI, since §13.19 showed ARI prefers the coarse level there; it's reported.

**An arm passes if all of these hold.**
1. *The lens,* one test for all arms. J → gJ moves it to T → gT, to 10⁻⁹, on every base and no-gauge pair. Median S₂·N at the top is below 2, and median χ at the bottom below m²/N, on 60 or more of 63. At most 240 levels (5 × 48) on every graph. Sweep time per graph is reported against `auto`'s.
2. *Seeds agree.* On 51 or more of 63, the modal count of stages holds on 9 or more of 12 seeds and every position stays within 2 times of its median. Applies to `chi` and `shelf`.
3. *Landsat.* Each of the four MAD and compiler-default rows reads 4 or more big clusters at purity 0.85 or more on 9 or more of 12. Only the upper stage reached this on all four in §13.19, read after the fact.
4. *No cost.* It holds `auto` on 60 or more of 63, counting moves between planted levels as held, with the strict count reported beside it. Failures are listed.
5. *Better than blind.* A list arm has to beat `null` of the same length: `null` fails criterion 3 or 4, or holds at least 5 fewer graphs. Otherwise the finding is "read several starts", and the stages add nothing.
6. *Short.* The modal list has at most 4 readings on every graph.

The arms get reported in this order: `chi up`, `hot`, `tfs`, then the rest as comparators. If more than one passes, the simplest goes forward, in the order `hot`, `chi up`, `shelf up`. `band` and `slice` come with a caveat: they get 3 chances beside `up`, so one passing alone needs its own pass before anyone relies on it.

**Reported, not criteria.** Which position each dataset wants. BlattHierarchy's planted levels reached at ARI ≥ 0.8. Between-seed ARI at each position. Cold-censored selections. How much of the lens a blind start passes on Landsat.

**What would change the framing.** If `tfs` passes, the lens alone was the fix and stages aren't needed. If `null` matches the lists, χ isn't locating anything. If Landsat's band and the toy's shelf both sit at one stage across their graphs, a rule like `hot` might replace T_fs, and that goes to the user as its own proposal.

**Not varied:** the budget, distance, LMP, the floor, generator seeds, θ, q beyond the arms named. No data outside these four sets, so this can't validate a rule. The constants above were fixed on the same graphs' curves.

## Amended after the sweep started, before any result was read

Two things in the criteria as first written would have decided a graph by construction. Only the plumbing test's one graph (iris, K = 10, Cauchy, base) had been looked at, one seed's rows.
- *Criterion 1's bottom.* A graph split below the floor or with no repair has χ above m²/N at every cold T, so no bottom could pass. The test becomes: median χ at the bottom is below χ₀ + m²/N, with χ₀ what the graph's own components above the floor give, (Σc² − c_max²)/N. The count against plain m²/N is reported beside it.
- *Criterion 4 on Landsat's four rows.* `auto` reads 0 of 12 there, and the README only excused ARI, but a position reading the target also changes the modal reading by construction. A position reading the target holds `auto` on those rows, and the count is reported with that said.

## Results (primer §13.21, C64)

No arm passes as written, so nothing goes to a default.
- Criterion 1, the lens: passes. Equivariant to 2·10⁻¹⁶, top 63 of 63, bottom 62 of 63 (53 against m²/N alone), 48–233 levels.
- Criterion 2: χ 52 of 63, shelves 41 of 63.
- Criterion 3, Landsat's four rows: `chi up`, `hot` and `chi band` 12, 12, 12, 12. `tfs` 0, 0, 0, 0. `chi slice` 4, 10, 12, 12. Shelves 0, 12, 0, 0. `null` 0, 12, 12, 12.
- Criterion 4, held of 63 against 60: `chi up` 47, `hot` 46, `chi band` 44, `tfs` 43, `null` 40, `chi slice` 38, shelves 35–37.
- Criterion 5: `chi up` beats `null`, 47 against 40, and `null` fails Landsat.
- Criterion 6: every list is 4 or shorter.
- So `chi up` passes all but 4, and `hot` all that apply but 4.

The bar in criterion 4 was wrong: §13.18's rule flags `auto` against itself on 18 of 63 (first six seeds, last six).

## Past what's written above

All after seeing the results, none with new sampling.
- `span-checks.py`: the rule on `auto` against itself, each list in its dataset's terms, what `hot`'s falls are, between-seed ARI, cost.
- `span-direction.py`: the verdict in each dataset's terms. `hot` 56 of 63 not a loss.
- `blattmix` in `lh-span2`: BlattHierarchy's clusters against the planted blocks. Under MAD and compiler defaults it's a per-branch mix.
- `span-blind.py`: blind slices and bands on BlattHierarchy, since `null` only tried `up`. Bands reach all three levels on 8 of 9, shelf bands on 6, shelf slices on 7, blind slices on none.
- `span-blind4.py`: four blind bands on all 63. Landsat 0, 12, 12, 12, not a loss on 46.

One thing was missed before the run: 12 iris graphs had no `auto` baseline on this grid, so they ran after the sweep (`run-span-base2.sh`).

## Built (2026-11-04, commit 6fa4298)

The user chose to build the lens before any held-out pass. `SpcScheduleHelpers.SpanBracket` and `--temperatures span`, beside `auto`, with no reading rule.
- `lenscheck` in `harness/Program-check.cs` (generated by `span-check-gen.py`, read by `span-check-read.py`) runs the shipped function beside the probe's `Span.Of` on all 63 graphs: lo, hi and levels are identical (`results/span-check.tsv`, `.out`).
- At `MinClusterSize` 1, the REPL's default and never swept: the bottom falls on 35 of 63, and the 63 take 55–240 levels, median 144, 8 at the cap. So the REPL warns, and SPC presets gained `min_cluster_size`.

## Files

- `harness/`: the probe (`Program.cs`, its project, and `Program-lens2.cs`, the copy with `blattmix` and χ₀), `span-rules.py` (the stage rules), `span-windows.py`, `span-read.py` (the criteria), the checks, and the run scripts. `span-pilot*.py` set the shelf rule's constants, with no labels.
- `results/`: `span-lens.tsv`, `span-resolved.tsv`, `span-stats.tsv` (each graph), `span-read.out` (the criteria and every position), `span-direction.out`, `span-checks.out`, `span-blind.out`, `span-blind4.out`, `span-tables.md`, and the raw `span-out/`, `span-base/`, `span-win/`, `curves/`, `span-blind/`, `span-blind4/`.
- The sweeps' co-membership caches (8.5 GB) stayed in the scratchpad and regenerate from `run-span.sh`.
