# Identity across temperatures (2026-11-04)

A record of an open discussion with the user about using persistence ideas to study how clusters evolve across
temperature, and of the reading it prompted. Nothing here has been run on data, nothing is built and no default moves.
The derivations marked below were checked by enumeration on small graphs.

**Where it stands.** Colors carry no identity, even between consecutive draws, but every draw is an overlay on the same
points, so identity can be resolved afterwards from sets, with no coupling between temperatures. For clusters between
two partitions the overlap table is already the correspondence. For competing states across T, shares can be matched by
unbalanced transport. A canonical barcode needs one line of partitions, so the most common state at each T carries the
barcode and the others stay as shares beside it. q becomes a second axis, with T × q as a surface.

## The frame

These are the user's, and they set the direction:

- **Same points, different overlays.** The points are identical at every T, and each draw is one overlay of colors or
  other measurements on them. This is a data analysis problem that should use that registration. The table has points as
  rows and overlays as columns, one per (T, replica, draw).
- **Replicates are data in an experimental design.** There's a literature on combining them, depending on context.
- **Identity is different from a color assignment.** The latent grouping evolves over T while the color names don't
  carry over. Names can be resolved post hoc, after colors are drawn, so temperatures stay independent.
- **q as a knob.** Runs at different q are variants, like replicates with one setting changed, which makes each
  temperature a small factorial. T × q is a thermal surface, by analogy with a thermal curve at fixed q.
- **Sorting for identity.** Sort geometrically so that names stay with points, roughly.
- **Long chains.** One long run per (T, q) cell, without n replicas per q.

## What the engine gives today

Sampler files are under `src/clustering/graphical/spc/runtime/core/sampler/`, scheduling and execution beside it under
`runtime/`, and the sweep readouts under `src/clustering/graphical/spc/profiling/`.

- **Colors.** One fresh uniform color per cluster, every draw (`SwendsenWang.cs:317`). A cluster keeps its color only by
  chance, 1/q. Bonds form only between same-colored neighbors (`SwendsenWang.cs:293`), so a new cluster lies inside one
  color region: each cluster re-forms on its own edges, and two clusters fuse only after drawing the same color. Chains
  start from random spins (`SwendsenWang.cs:208`).
- **Names already exist.** `KeptPartitions` labels each node with the smallest node index in its cluster, comparable
  across draws and replicas (`Accumulator.cs:99`). I found no code reading these across T, though I searched only the
  sampler and the kept-partition docs.
- **Replicas are nested in T.** The seed hashes (base seed, T, replica, round) (`SpcSeedHelper.cs`), so replica 3 at two
  temperatures shares only a label. q isn't in the hash, so variants at one (T, replica) would share a seed under one
  base seed.
- **Point × T and edge × T tables exist.** Mean cluster size and giant participation per node (`SweepLandscapes.cs`),
  and co-membership per edge (`SweepEdgeCurves.cs`). Both are pooled across replicas, while the runs keep each replica's
  values.
- **Mixed q.** Each task carries `Q` (`SpcTaskSpec.cs:13`) and each frame records it. `SweepEdgeCurves` refuses mixed q.
  A search finds no such check in `SweepCurves` or `SweepLandscapes`, so frames should be split by q first. The sampler
  requires q ≥ 2 (`SwendsenWang.cs:137`).
- **Zigzag.** `GraphZigzagFast` and `GraphFrameZigzag` (not read for this note).

## Ways to name a set

Colors drop out first, since within a draw they only say which points share a set, and permuting them can't help: a
color is a random union of clusters.

- **Overlap table.** Two partitions a and b on the same points give a contingency table, a coupling with cluster sizes
  as marginals. No optimization is needed. Joining blocks that overlap gives a bipartite graph with k(a∧b) + k(a∨b) −
  k(a) − k(b) loops, the gap in the supermodularity proved in `lean/Lemmas/RandomCluster.lean`. It's zero when they
  nest. Weighted by overlap, it separates a hub moving between arms from a ring's rotating cuts, which ARI lumps
  together.
- **Leader.** Rank the points once and name each cluster after its top-ranked member (the user's sort). Node index is
  the rank today, and any other can be applied to kept partitions afterwards. A random rank shares a name between
  clusters with overlap J with probability J (MinHash). A density rank should keep names on dense cores, near ToMATo
  (from memory). The sort key matters: see McLachlan below.
- **Recurring sets.** Compare all clusters from all overlays by overlap, whichever T or draw they came from, and treat
  strongly overlapping ones as one identity (meta-clustering, Strehl & Ghosh 2002, from memory). Persistence is then the
  share of overlays at each T containing it. Chaining is the main risk, since a peeling sequence can link a big set to
  much smaller versions of itself, and overlap alone hides containment.
- **Tracking along a run.** Consecutive draws are linked by the chain. Two clusters sharing more than half of each pass
  the label on. The rest are real splits or merges: maximum-weight assignment, or a min-cost flow over the whole run
  (Zhang et al. 2008, from memory), with links that skip a draw. H₀ zigzag gives a canonical barcode instead, with a
  split starting a bar and a merge ending one, and union or intersection between draws is a real choice. It needs
  consecutive draws, a label can drift onto a different set, and track length depends on the sampler's dynamics, while
  shares at T depend only on the model.
- **Subspaces, no names at all.** A partition defines V_P ⊂ ℝ^N, the functions constant on its clusters, which is H⁰
  with real coefficients and doesn't change under renaming. The average projection P̄_T has entries ⟨1[i, j
  together]/|C_i|⟩ and trace E[k]. Between temperatures, the singular values of U_T^⊤ U_T′ are cosines of principal
  angles, and on a nested stack with exact agreement their count above a level is the H₀ persistent Betti number.
  Singletons have eigenvalue 1 and need a size gate, and equal eigenvalues let eigenvectors rotate, so compare
  subspaces, not vectors.
- **Modes.** The user asked why average rather than take a mode. Averaging blurs competition: with ab|c in 60% of draws
  and a|bc in 40%, the leading direction is a blend no draw contains. An exact mode rarely exists at large N, a
  per-point mode needs names first, and a per-pair mode is the existing cut, whose closure can build clusters no draw
  contains. So group draws at each T by variation of information, giving modes with shares, and average only inside a
  mode.

## Correspondence across temperatures

1. **Clusters:** the overlap table between two partitions.
2. **States:** match modes at adjacent T by unbalanced transport, with shares as mass and variation of information as
   cost. Any such pairing is a modeling choice, since nothing links independent draws.
3. **Barcode:** one line through each T's most common mode gives canonical bars, with other modes reported as shares.
4. **Windows and q:** treat them as Mapper lenses (see reading notes).

Two overlays can be compared directly at any pair of temperatures. They share no noise, so what they have in common
comes from the graph. Only following a single draw through T is undefined. A T × T matrix of mean agreement between
draws would show regimes as blocks, since the diagonal is a temperature's agreement with itself and shuffling draws
gives a test.

## q as a knob

**What q is.** Q limits colors, not clusters: at q = 2 a draw can still have up to N clusters. It enters the model only
as the factor q^k(ω) on each bond configuration, so T sets how strong contact is and q sets how often it can happen.

- **What can move.** At fixed T, d⟨f⟩/d ln q = Cov(f, k), at most sd(f)·sd(k), so q can't move much where k barely
  fluctuates. For an event A, d logit P(A)/d ln q = E[k | A] − E[k | not A]. For a pair, its size says how many fewer
  clusters there are in draws where they're together: a cooperativity at each T.
- **Topology in the weight.** With uniform J, the weight is v^{b₁}(q/v)^{b₀} up to a constant, with b₀, b₁ the Betti
  numbers of the bond graph (the Tutte polynomial, Sokal 2005, from memory). T × q is an exponential family in those
  statistics. T's statistic per draw is Σ_{e∈ω} J_e/p_e, whose expectation is −⟨H⟩, so energy can stand in for it.
  Locally, q acts as a rescaling of T when that statistic and k are perfectly correlated across draws.
- **A filtration axis.** At fixed T, connection probability only falls as q rises (Grimmett Thm 3.21). Raising q by a
  factor c lowers it no further than dividing every edge's odds by c at the old q. That bound covers increasing events
  only. A group that stays separate down to q = 2 survives constant exposure, while bond strengths stay fixed.

**Frozen blocks.** For two rigid blocks with coupling J: ⟨n⟩ odds are (e^{βJ} − 1)/q, so logit ⟨n⟩ is f(T) − ln q and T
and q are additive with slope −1. The ½ crossing is T = J/ln(q + 1) (the δ̄ form is in
`docs/Clustering/SpcAnalysis.md:225`). So q rescales frozen merges and keeps their order, one reason it often seems not
to matter. At the crossing, the chain between joined and apart is symmetric, flipping with probability 1/(q + 1) per
step.

| q | merge T, relative to q = 2 | 10–90% width (e-folds) | mean stay (draws) | autocorrelation time (draws) |
|---|---|---|---|---|
| 2 | 1.00 | 2.7 | 3 | 1 |
| 3 | 0.79 | 2.4 | 4 | 1.5 |
| 5 | 0.61 | 2.2 | 6 | 2.5 |
| 20 | 0.36 | 1.5 | 21 | 10 |
| 40 | 0.30 | 1.2 | 41 | 20 |

**Factorial design.** T × q on one common temperature grid, wider by about 1.2 e-folds for q from 2 to 40, with q dense
at the low end and an anchor at 20, where every existing result sits: 2, 3, 4, 5, 6, 8, 10, 20. One deep chain per cell,
with error bars from batch means and from kept partitions thinned to about one per autocorrelation time. Neighbors in T
and q act as partial replicates, and the T × q interaction is the signal, so it can't serve as the error term. Add
chains only at flagged cells. Those should include a start from the ordered side, because random starts all approach
from the disordered side and can share a bias across cells; I haven't seen an ordered start option in the sampler. On
the logit scale the slope in ln q can be estimated two ways, across variants and from draws within a cell, and
disagreement flags poor mixing.

**Thermal surface.** Each fixed-q curve is a slice, and each edge's and point's curves become surfaces. A peak becomes a
ridge T*(q). BASS could fit it, with `BassFit.Sensitivity` splitting T, q and T × q effects. On the logit scale, adding
ln q collapses frozen-block slices, leaving a residual surface. Frozen merges narrow as q rises while their gaps in ln T
stay fixed, which resembles scale-space analysis (Witkin 1983, from memory), with ridges forking as q rises. Cautions:
few q levels, ridges slant on a common grid so fits suit warped coordinates like T·ln(q + 1), and real jumps at high q
would be blurred by a smooth fit.

**Two disagreements, both testable.**

- *The user:* q's range should be bounded by a realistic number of clusters, and q = 20 is ceremonial on iris. *Me:* the
  low end should matter most, because the change in character sits there (first-order above q = 2 on dense graphs, from
  memory, and above q = 4 on the 2D lattice), but the frozen-block law has no K in it, and 5 → 20 narrows a merge more
  than 2 → 5. Synthetic graphs with K = 3, 6 and 12 blocks separate these: does the point where the residual surface
  flattens move with K?
- *The user:* fewer colors should make identity tracking easier. *Me:* palette size doesn't enter when tracking sets,
  and high q holds states longer, with low q giving more independent information per draw and fewer stuck chains. Iris
  may not fit: versicolor and virginica share a thick boundary, so those blocks aren't rigid when they merge. Comparing
  q = 4 and q = 20 at each q's own peak decides it.

**Set aside.** Recoloring each cluster with probability ρ and otherwise keeping its color preserves each draw's law
exactly (checked), and makes colors carry identity natively, at autocorrelation growing roughly as 1/ρ. The user prefers
post-hoc names with temperatures uncoupled, so this isn't pursued.

## Other readings worth keeping

- **Melting curves.** Points or pairs as subjects and T as dose. Thermal proteome profiling (Savitski et al. 2014; Tan
  et al. 2018; Childs et al. 2019, from memory) fits a curve per protein, and co-melting reveals complexes. An edge's ½
  crossing is its melting temperature, today's stack is H₀ filtered by that number, and curve width goes unused.
  Grouping points by whole curves could show clusters that melt alike without sharing a branch.
- **Support curves.** S(B, T), the share of partitions at T containing block B. A block with S > ½ is a component of the
  ½ cut on ⟨n⟩, for any sample of partitions with connected blocks. At the exact measure, majority blocks at different T
  can't cross, so they form a tree inside the pooled tree. A pooled node with low support exists only in the average,
  like a ring joined by its neighbor pairs. The catch: exact blocks rarely reach ½ at large N, and tolerant matching
  (the transfer bootstrap, Lemoine et al. 2018, from memory) may not keep this argument. Clade support in phylogenetics
  is the closest literature.
- **A metric at each T.** P_T(i, j apart) satisfies the triangle inequality, so each T is a finite metric space, and at
  the exact measure it shrinks on cooling. Across a stack, ∫ P_T(i, j apart) d ln T reduces to cophenetic distance when
  nested.
- **Name statistics.** A name's survival, P(its point connects to no higher-ranked point), falls on cooling at the exact
  measure, and survivals sum to E[k] for any ranking. Mean 1/|C_i| is the ranking-free version.
- **Replicate stages.** Pool before the filtration (today), average summaries, compare decisions (IDR, stability
  selection), or split so one half selects and the other measures. Disjoint stacks of one replica per T are i.i.d., so
  ranks among them give conformal coverage of the next run's barcode, with no stability theorem.
- **Fast versus dense.** Peaks in fast runs' curves are replicate peak calling, and reproducible peaks could set windows
  for dense runs.
- **Single-cell time courses.** Waddington-OT (from memory) faces the same one-measurement-per-unit problem.

## Reading notes

Sections read: Dey & Wang 2021 §4.3–4.4 (openings), §4.6, §7.1, §7.3 (opening), §8.2, §9.3 (opening), and definitions
11.18 and Fact 11.3 from search; Peyré & Cuturi 2020 §10.1–10.2 and §10.6 (opening); McLachlan & Peel 2000 §4.9.
Stability in Dey & Wang §9.4 and matching distance in §12.3 are unread.

- **Uniqueness needs a line (Dey & Wang Thm 4.10).** The interval decomposition is unique for a linear zigzag. Plural
  modes branch off that line, so I don't expect one canonical barcode for them (my reading). T × q is also off it, and
  the rank invariant is complete for one parameter only (Fact 11.3).
- **No bond graphs needed (§4.4).** Arrows may be simplicial maps, and collapses reduce to inclusions (Prop 4.12), so
  P_T → join ← P_T′ fits post hoc.
- **Bar types (§4.3).** Four end types, and on a fixed point set ends are splits, merges or the ends of the range.
  Diagrams match only same-type bars, and interleaving is awkward for zigzag, which matters for comparing barcodes
  across q.
- **Windows (Remark 9.1, §9.3.1).** Sublevel covers of T give merge trees and overlapping windows give a relaxed Reeb
  graph. Widening windows is a tower, with each component inside exactly one coarser one (Prop 9.13). With identity
  links between copies of a point, a window's component is a block of the join across that window (my reading).
- **Reeb space (end of §7.1).** A map to ℝ² gives a Reeb space, so a (T, q) lens is the identity version of the surface,
  with Mapper as its discretization. Reeb graph distances are strongly equivalent, and equal on merge trees (§7.3).
- **Per-point diagrams (§8.2).** Persistence distortion gives each point its own diagram, and needs Hausdorff only
  because its points aren't matched. Its stability result is for metric graphs and doesn't transfer automatically.
- **Unbalanced transport (Peyré & Cuturi §10.2).** Classical transport forces mass to split when shares differ, and a
  relaxation parameter trades transport against mass change. It has been used for cell dynamics. The parameter is a
  modeling choice.
- **Gromov–Wasserstein (§10.6).** Built for measures on different spaces. These are the same points, so it isn't needed.
- **Multimarginal (§10.1).** One tensor axis per temperature, so adjacent pairs look like the practical version.
- **Relabel afterwards (McLachlan §4.9).** Imposing an ordering during sampling can bias results, and the recommendation
  it reports is to simulate freely and relabel afterwards (Celeux; Stephens). In Richardson & Green's example, ordering
  by means gives clearer marginals and ordering by variances leaves them bimodal, so bimodal statistics per name would
  warn of a bad ranking.

## Status of the claims

**Checked by enumeration** on three small graphs (6 or 5 nodes), q from 1 to 20 and 8 values of β, all passing; the
script is `identity-across-temperatures/check_identity_note.py`, beside the `lineage-bars` harness in the side-session
notes, outside the repo:

- the logit slope in ln q, and survivals summing to E[k]
- leader survival and mean 1/|C_i| falling on cooling
- the triangle inequality for P_T(i, j apart)
- majority blocks as components of the ½ cut, and not crossing across T
- connection falling with q, and the bound with edge odds divided by c, both with a parameter per edge
- frozen-block odds against the two-state chain's stationary odds, and the table's formula values
- lazy recoloring preserving the Edwards–Sokal joint law, on one 4-node graph
- dim(V_a ∩ V_b) = k(a∨b), and the loop count being nonnegative, on random partitions

**By hand only:** the rank-function reduction to persistent Betti numbers, the cophenetic reduction, the local rescaling
criterion, and the reach of conformal coverage.

**Not formal.** None of this is in Lean. Extending the comparison inequality to q should follow the same Holley argument
as `rcProb_mono`, and the majority-block lemma should follow from `connCutGraph_reachable_mono`.

**Citations marked from memory** weren't checked against the sources.

## First looks, none run

1. Points × draws at one temperature, raw colors beside tracked labels.
2. A T × T matrix beside a heatmap of points × T colored by each point's most common name, on iris and Landsat.
3. q = 4 against q = 20 on iris, at each q's own peak, keeping every draw: track lengths and surviving sets.
4. K = 3, 6, 12 block graphs across q, for the bound on q's range.
5. A toy ring against a blob, reading support on the pooled tree.

Before any T × q run, q needs to enter the seed hash, or each variant needs its own base seed.
