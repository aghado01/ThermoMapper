# Lineages against H0 bars on a nested stack (2026-11-04)

The user asked for this comparison. Exploration, one pass, no default moves, nothing here validates a reading.
Written before any number was read.

## The claim under test

Primer §12.21: "the T-stack is then a dendrogram, and a lineage's lifetime is the persistence of a one-parameter
filtration". `LineagePersistence` ends a lineage at a split and starts both children fresh, and carries on through a
sliver below `SplitShare`. The elder rule decides by age alone. So I expect the two to cut the same tree differently,
and I told the user a tight clump can take the long bar. Neither has been measured.

## The object

- **Graphs and sweeps:** the span lens's 63 graphs and 12 seeds, read from its caches, no new sampling. Graphs as
  reported in `span-lens/results/span-resolved.tsv`. Budget 1 x (200, 400), the low one. θ = 0.5, q and `MinClusterSize`
  m per the lens (iris 5, BlattHierarchy 10, toy and Landsat 50), `SplitShare` 0.10, representative `Largest`,
  whole lens, default longest-lived reading.
- **Nesting:** sampled stacks needn't nest, so I record how often they do. Then each edge's δ̄ column is made
  non-increasing in T by unit-weight PAVA, the stack is cut from that, and `PartitionNesting.Holds` must be true.
  Lineages are resolved on those same monotone columns.
- **The filtration, gated by size.** Cooling adds edges. A plain graph filtration has every vertex at the hot end, so all
  H0 births tie and the elder rule has nothing to decide. So: a vertex enters when its component first reaches m
  points, an edge at the colder of its own entry and its endpoints' gate. That is a real filtration, since a component
  only grows on cooling and faces enter no later than cofaces. Bars come from the engine (`PersistentHomology` on a
  `CellFiltration`) and independently from a direct walk of the stack, which must agree.

## Questions

1. **Do lineages equal bars?** Share of all lineages, and of selected lineages, whose level span equals a bar's. How many
   lineages the K longest bars run through. How long the longest bar is against the longest lineage.
2. **Who keeps the long bar?** At each cold cluster with two or more children of size m or more:
   - a *sliver event* (second child's share below 0.10), where the lineage follows the largest child: how often is the
     elder child not the largest? This is the tight clump.
   - a *split* (share 0.10 or more), where the lineage ends: same count.
   - Age ties go to the larger child, so they can't count as disagreement. Reported separately.
3. **Does the reading change?** K = the number of selected lineages. The bar reading keeps the K longest bars and
   absorbs every other bar into the elder it merged with (persistence-based clustering, as in ToMATo). Agreement
   between the two readings on points both assign, each one's coverage, and each one's ARI and purity against the
   dataset's labels, as a report only.
4. **Controls.** Share of sampled stacks that nest as sampled. Share of (edge, level) pairs PAVA moves across θ. ARI of
   the lineage reading on raw columns against monotone columns: if this is low, nothing here carries back to raw stacks.
   In the plain filtration, how many of the K longest bars belong to a side smaller than m.

## How I'd read it

- **The primer's sentence** stands as written if 90% or more of selected lineages equal a bar on most graphs. Below
  that it needs a qualifier, which goes to the user: this session doesn't edit the primer.
- **The tight clump** is real if sliver events show the elder not largest at any appreciable rate. I set no threshold
  on this, and the rate is the result.
- **The reading** doesn't care which decomposition is used if agreement is near 1. If it's low, the choice matters, and
  the labels only say which way on these four datasets.

## Checks, each must pass or the run is wrong

- Monotone columns are monotone, and every monotone stack nests.
- Area identity: Σ lineage lifetimes = Σ bar lengths = Σ over levels of cell width × clusters of size m or more. Both
  decompositions partition the same area, so this holds exactly or a walk is wrong.
- Engine bars equal direct bars as multisets, in level units. `barcode_unique` says they must, whatever the ties.
- The direct chain walk reproduces `Resolve`'s `AllLineages`.
- A hand-built stack with known bars and lineages passes (`selftest`).

## Not varied

Budget, θ, m, `SplitShare`, the representative, generator seeds. No held-out data. Peripheral capture is off, as it is in
the lens.

## Past what's written above

Added after the first read of the whole lens, and labelled that way in the note.
- **The `hot` window.** The whole lens reads one cluster on 268 of 756 stacks, where Q3 says nothing. So the same code ran on
  each seed's window from its hottest χ stage up (`span-win/*.stages.json`), PAVA inside the window. This is the reading the
  lineage work is looking at.
- **Margins at events.** An exact 0 of 777 needed showing to be real, so events now carry the next oldest child's reach
  and size, and `read.py` recounts from the raw rows.
- **Scale,** in `scale.py`: between-seed ARI from the lens's `pairs.tsv`, Landsat's purity beside its ARI, and the selected
  share by splits in the window.
- **Reading on the 90%.** The note applies it to the whole lens, as written; `hot` is reported beside it.

## Results

The note is `issues/ph/lineages-against-bars.md` in the repo. No check failed on any of the 1512 stacks.

## Files

- `harness/`: `Program.cs`, `lineage-bars.csproj` (references an export of `28ba17e` at `head/`, via `git archive HEAD src
  projects Directory.Build.props project.toml presets`), `read.py`, `scale.py`.
- `results/whole/`, `results/hot/`: one `<graph>.tsv` and `<graph>.events.tsv` per graph. `read-whole.out`, `read-hot.out`,
  `scale.out` are the printed tables.
- To rerun: `lineage-bars.exe selftest`; `lineage-bars.exe run <span-cache> <outDir> 3 - <graphs>` for the whole lens, or
  with the `span-win` folder in place of `-` for `hot`. `LB_MUTATE=younger` or `chain` breaks a rule, to show the checks fail.
