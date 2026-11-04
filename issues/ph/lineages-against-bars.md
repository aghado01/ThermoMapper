# Lineages against H₀ bars on a nested stack (2026-11-04)

The user asked for this comparison. It's one exploratory pass on existing sweeps, on the four datasets that fixed the
readings' constants, so it validates nothing and moves no default.

**The result.** A lineage isn't a bar. The two partition the same area of the merge tree and cut it differently, and
only about a third of lineages equal a bar. Given the same number of clusters, they still read the same partition on
most stacks. The tight clump I predicted, a small piece keeping the long bar while the big cluster's bar dies, never
happened at a sliver.

## What was compared

- **Graphs and sweeps.** The span lens's 63 graphs × 12 seeds, read from its caches, no new sampling: 756 stacks. Each
  graph is as resolved in `span-lens/results/span-resolved.tsv`, and none has a protect set. Budget 1 × (200, 400), the
  low one. θ = 0.5, `SplitShare` 0.10, representative `Largest`, m as the lens set it (iris 5, BlattHierarchy 10, toy and
  Landsat 50). Built from an export of `28ba17e`. `src/clustering`, `src/graphs` and the generators are unchanged since
  the sweeps' build, `394e3ef`.
- **Two windows.** The whole lens, as written before the run. Then `hot`, from each seed's hottest χ stage up, added
  after the first read because the whole lens gives one cluster on 268 of 756 stacks and so says nothing there.
- **Nesting.** Almost no sampled stack nests, so each edge's δ̄ is made non-increasing in T by unit-weight PAVA, and the
  stack and the lineages are read from those columns. Under `hot`, PAVA runs inside the window.
- **The filtration, gated by size.** A plain graph filtration has every vertex at the hot end, so all H₀ births tie and
  the elder rule has nothing to decide. So a vertex enters when its component first reaches m points, and an edge at the
  colder of its own entry and its endpoints' gate. Its components at any level are exactly the clusters of size m or
  more, and cooling only adds cells, so it is a filtration and `exists_intervalDecomposition` applies.
- **The bar reading.** With k the number of selected lineages, keep the k longest bars and absorb every other bar into the
  elder it merged with, as in ToMATo. k is borrowed from the lineages, so this reading has no way yet of choosing its own.

## Results

| | Whole lens | `hot` |
|---|---|---|
| Sampled stacks that nest | 0 of 756 | 1 of 756 |
| Lineages equal to a bar | 6201 of 17424, 36% | 5949 of 16155, 37% |
| Selected lineages equal to a bar | 450 of 2402, 19% | 1527 of 3809, 40% |
| Graphs with 90% or more of selected lineages equal | 0 of 63 | 5 of 63 |
| Sliver events where the elder isn't the largest child | 0 of 777 | 0 of 163 |
| Splits where the elder isn't the largest child | 1226 of 7319, 17% | 900 of 6463, 14% |
| Stacks reading more than one cluster | 488 | 756 |
| ...where the two readings agree at ARI 0.9 or more | 441 | 661 |
| ...where they agree exactly | 433 | 622 |

**Lineages and bars.**
- **The shape of it.** A bar follows the elder child through each split, while a lineage ends there. So a kept bar runs
  through a median of 3–5 lineages on the whole lens and 1–3.5 under `hot`, and the longest bar is 1.1–2.5 times the
  longest lineage.
- **Where they coincide.** On a leaf branch that loses to the elder where it joins its parent, or that never joins one. So
  narrowing the window raises the share: under `hot`, selected lineages equal a bar 64% of the time on stacks with 0–2
  splits in the window, and 19% with 11 or more. All 5 graphs at 90% have a median of at most 0.5 splits per stack. By
  dataset under `hot`: Landsat 73%, BlattHierarchy 25%, iris 16%, toy 1%.

**Who keeps the long bar.**
- **At a sliver, the largest child, every time.** Age ties (18 on the whole lens, 8 under `hot`, all Landsat) went to the
  larger by construction, so a finer grid could flip those. On this grid that bounds the rate at 2.3% and 4.9%.
- **It isn't a near miss, except on Landsat.** The largest child's reach beats the sliver's by a median of 0.93 e-folds on
  iris and 3.15 on toy. On Landsat it's 0.37, with the bottom decile at one rung (0.093). My reading, not tested: a big
  cluster's age is the hottest of its many sub-clumps, so age and mass go together.
- **So the two cut differently only at splits.** At a sliver the bar and the lineage both follow the largest child. At a
  split the lineage ends and the bar goes on.
- **At a real split, a comparable child.** The elder is smaller at 14–17% of splits. It's no clump: the median ratio of
  largest to elder is 1.2–2.0, and the largest is 5 or more times the elder on 1 of 1226 (whole) and 1 of 900 (`hot`).
- **Half of those are close to the grid.** The elder's median margin is one rung on toy and Landsat and 0.19–0.28 e-folds on
  iris and BlattHierarchy. 525 of 1226 and 464 of 900 are within one rung. Another 19% of splits are age ties.

**The readings.**
- **Mostly the same partition.** Median agreement is 1.000 on iris, toy and BlattHierarchy under both windows.
- **Where they differ.** Landsat (61 of 216 stacks below 0.9 under `hot`, 18 of 32 on the whole lens) and BlattHierarchy
  under MAD and compiler defaults (median 0.71–0.74).
- **The labels don't favour bars.** Absorption keeps more points. Under `hot`, on Landsat's graphs other than MAD and
  compiler defaults, the bars cover 3603 points against 2612 at purity 0.77 against 0.89. ARI counts an unassigned point a
  singleton, so it rewards that coverage: under `hot` bars score higher on 191 stacks and lineages on 80. BlattHierarchy
  goes both ways, with bars 0.88 against 0.79 under compiler defaults and 0.70 against 0.78 under MAD.

**Controls.**
- **Monotonizing moves little.** PAVA takes a median of 0.08–0.25% of (edge, level) pairs across θ, 0.77% at most. The
  lineage reading on raw columns against monotone ones has median ARI 1.000, or 0.990 on Landsat under `hot`.
- **It does sometimes move a reading.** ARI falls below 0.9 on 24 of 756 stacks (whole) and 62 of 756 (`hot`, 36 of them
  Landsat), and the worst case goes from 13 clusters to 3. For scale, the same Landsat reading's between-seed ARI has a
  median of 0.63 over its graphs, so monotonizing typically moves it less than changing the seed. These results carry
  back to raw stacks about as far as those numbers say.
- **The plain filtration doesn't read clusters.** Among its k longest bars, the share belonging to a side smaller than m
  is 47% on Landsat under `hot` and 36% on the whole lens, 31% on toy (whole). It's small on iris and BlattHierarchy.

## What this changes

- **The primer's sentence needs a qualifier.** §12.21 says a lineage's lifetime "is the persistence of a one-parameter
  filtration", and it missed the 90% I wrote down beforehand on every graph, whole lens. This session doesn't edit the
  primer, so the wording is offered here: *on a nested stack, a lineage is a branch of the merge tree between splits, and
  an H₀ bar of the size-gated filtration follows the elder child through each split, so a long bar runs through several
  lineages. Both partition the same area of the tree. They coincide only on a leaf branch that loses to the elder where it
  joins its parent, or that never joins one.* `issues/spc/bayesian-clustering-design.md:74` also writes "a lineage's bar
  endpoints", which reads the same way.
- **Proofs about bars reach lineages only through the tree.** `exists_intervalDecomposition` and `barcode_unique` cover
  the gated filtration, and no theorem covers the lineage's chains. The area identity below is the one exact link.
- **Almost no raw stack is a filtration at this budget.** Anything stated for nested stacks needs monotone columns or a
  restoration first. With `restoreNesting` off, `HierarchyEom` would abstain on all but one of these stacks. The columns
  estimate a connection probability proved monotone (`connProb_mono`, `lean/Lemmas/RandomCluster.lean`), so at equilibrium
  this non-nesting is estimation error and should shrink with budget. That shrinkage wasn't measured here.
- **Nothing here argues for reading bars instead of lineages.** They differ mostly where the data is hardest, and the bars
  give up purity there.

## Checks

All pass on all 1512 stacks, and on 1500 seeded random nested stacks in `selftest`.
- Monotone columns are monotone, and every monotone stack nests.
- **The area identity, exact:** Σ lineage lifetimes = Σ bar lengths = Σ over levels of cell width × clusters of size m or
  more.
- `PersistentHomology` on the gated `CellFiltration` equals the direct walk's bars as multisets in level units.
- The direct chain walk reproduces `Resolve`'s `AllLineages`, by each lineage's ends.
- PAVA matches the min-max formula on 2000 random vectors.
- **The checks can fail.** Breaking the elder rule trips 758 of them in `selftest`, and doubling the share in the chain rule
  365.

## Not covered

- One budget, the low one, and no held-out data. The window was chosen after seeing the first read.
- Only m, θ and `SplitShare` as the lens set them. A smaller m could let more slivers count.
- k is borrowed, so agreement says nothing about how a bar reading would choose its own number of clusters.
- Monotonizing is unit-weight PAVA, one way among several. BARS-monotonized columns or `HierarchyEom`'s restoration could
  read differently.
- No zigzag: the raw stacks were not read as they stand.
- The HDBSCAN condensed-tree resemblance is at the level of the rule, from reading the code. It wasn't measured here.

## Files

Harness, questions as written before the run, rows and full output, outside the repo:
`side-session-notes/local_c9eea969-0273-4b1b-a5ca-ed19101be413/lineage-bars/`. Start with `README.md`. The per-graph
tables are in `results/read-whole.out` and `results/read-hot.out`. The sweeps' caches stay with the lens, in the scratchpad
its `run-span.sh` writes to.
