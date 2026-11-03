# Bayesian clustering ideas for SPC: design note

Status: started 2026-10-30. The closure gap is built (2026-11-01), measured on data (2026-11-02; see "What the measurement found"), probed on Landsat (2026-11-02; see "Results of the probe"), tested against closure at a higher θ (2026-11-02; see "Results of the θ test") and against a closure that won't cross a thin cut (2026-11-02; see "Results of the thin-cut test"); the rest is discussion. Decisions are the user's and are marked as such; the rest is proposal.

Sources: `opus-gvdv-ideas.md`, Threads 1, 5 and 7 of `opus-spin-glass-threads.md`, `sw-retention-brief.md`, and the stopping primer's §13.3 and §13.6 (`issues/spc-BARS/stopping-theory-primer.md`).

Symbols are the module's: τ_ij = P(i↔j), the pair connectivity; p_e = 1 − e^(−J_e/T), the bond probability; δ̄ = ((q−1)τ + 1)/q, the spin agreement that `DenseTStack` cuts.

## The user's direction

2026-10-30, stated as a leaning and not as settled: pursue what fits a uniform temperature grid read by BARS. The primer's §13.3 records the grid decision; its decision 4 keeps BARS out of the stopping rule, so anything read through BARS here is on the measurement side.

2026-11-01, the user's decision: build the closure gap first, as idea 1's first step (see "Built"). No other idea is chosen yet.

## Status

| # | Idea | Status | Fit with the uniform grid and BARS |
|---|---|---|---|
| 1 | Loss cuts, against closure | The gap is built and measured, and a search over pair loss is tested and left with the tests. On Landsat at K = 5 it read worth building on that graph, on K = 10 and 20 no gain. Closure at a higher θ doesn't replace it. A closure that won't cross a thin cut does, on K = 5, and does so on the graph alone with no sampling, so the search isn't worth building on this evidence | A column kernel for `DenseTStack`; the closure gap below is a curve over T. The search's partitions don't nest, so a T-stack can't read them as they stand |
| 2 | VI cuts, credible balls, the robust window | Unblocked by S5 (`a777534`, SPCX v12) | Per-T. r(T) and each block's window are curves |
| 3 | Free energy by temperature integration | Unblocked by S1 | ⟨E⟩(T) on the grid. BARS can fit it, and integrating its draws gives ln Z(T) with a band |
| 4 | Calibrating T by posterior odds | Research-scale. The §5.2 rule does not carry over: the single cluster has the lowest H, and the colouring count q(q−1)…(q−k+1) favours splitting, the reverse of Bissiri et al. | Competes with the T* readout. Lowest priority |
| 5 | Split temperatures | Buildable once lineages are read as trajectories | BARS with `StepBasis` on the trajectory, so no separate changepoint engine |

Ways of sampling, by fit:
- **Fit:** an ordered start for half the replicas (the primer's C46, E16). Each T stays an independent task.
- **Fit at a cost:** parallel tempering on an unadapted grid. It runs the temperatures in step, which gives up per-task parallelism and the per-task SPCX cache.
- **No fit:** flat-histogram reweighting, which replaces the grid, and adaptive placement, which §13.3 rejects.

## Results on closure, the loss and the window

All are exact and hold for J ≥ 0, q ≥ 1 on any graph. "Enumerated" means a scratch check over every bond configuration of a small graph: a path of 6, two triangles joined by a bridge, and 30 random connected graphs (n = 6, 7–10 edges, J uniform on [0.15, 1]), at q = 2 and 20, 10 values of β from 0.3 to 20, and θ = 0.3, 0.5, 0.7. No SW and no dataset. Its counts check the algebra and aren't evidence about data.

**The θ-graph is the Bayes act; its components aren't.** With costs a for a false merge and b for a false split, the rule for each pair alone links at τ_ij > θ = a/(a + b). So the θ-graph minimizes the loss when transitivity isn't required. Closure then imposes transitivity by honouring every θ-edge and paying for every pair that drags in. The cost-ratio reading of θ belongs to the graph.

The code's θ cuts δ̄, not τ, so the cost ratio is a/(a + b) = (qθ − 1)/(q − 1): 0.474 at q = 20 and θ = ½.

**The closure gap.** With a = θ and b = 1 − θ, over edges (link loss) or all pairs (pair loss):

L(θ-graph) ≤ L(optimal partition) ≤ L(closure), and Δ(T, θ) = L(closure) − L(θ-graph) = Σ |τ_ij − θ| over the pairs where closure contradicts the pair's own rule.

- Δ = 0 makes closure exactly optimal, and Δ bounds what any search can gain. It takes one pass per column, with no search and no nesting question.
- Over edges only false merges can occur. Over all pairs both directions can: on K₂,₄ at q = 2 and β = 1, two poles with no edge between them read 0.700 through four parallel routes, while every edge reads 0.647. Between those levels the closure parts a pair whose own rule links it.
- It differs from the primer's pair indecision Σ_e τ_e(1 − τ_e), which is large wherever edges are undecided, even if closure is right.
- *Enumerated:* no sandwich violations in 3,840 cells; the formula holds to 1e-15.

**The gap integrates to a cophenetic distortion.** ∫₀¹ Δ dθ = ½ Σ (m_ij − τ_ij)², with m_ij the level at which closure first joins the pair, the max-min of τ_e over graph paths. Each pair contributes over the θ between τ_ij and m_ij, and ∫(θ − τ)dθ = (m − τ)²/2. Thread 1 says single linkage outputs the subdominant ultrametric, so this is that test on the connectivity field: a single number per T, and a curve over the grid. *Enumerated on the path:* 0.694027 both ways.

**Connectivity is a metric.** τ_ij ≥ τ_ik·τ_kj, by positive association (Grimmett 2006, Theorem 3.8; {i↔k} and {k↔j} are increasing, and together they imply {i↔j}). So −ln τ is a metric at every T, shrinking as T falls by C17. L hops at level θ guarantee only θ^L. On a tree it is an equality. *Enumerated:* no violation over all triples beyond rounding (−2e-16).

**Link loss can't see a pure chain; pair loss can.** A path at q = 20 with every edge at τ_e = 0.6 has τ(0, j) = 0.6^j exactly. At θ = ½ closure joins the whole path. Link loss agrees (1.00 for both). Pair loss drops from 4.79 to 2.41, by cutting the path into three pairs. This replaces the earlier statement that the edge form "reduces chaining": on a chain it does nothing.

**Nesting.** Measured, not proved, on these tiny graphs:

| Optimum | T-steps that fail to nest | θ-steps that fail to nest |
|---|---|---|
| Link loss | 0 of 1,728 | 0 of 1,280 |
| Pair loss, path | 2 of 54 | 3 of 40 |
| Pair loss, the rest | 0 of 1,674 | 0 of 1,240 |

An earlier message in this discussion said that loss optima carry no nesting guarantee and that a constrained search is needed. The guarantee is still missing, but the failures seen are confined to the chain. The likely cause, not inspected, is that segments of different lengths on a path can't refine each other. A loss-optimal family isn't a filtration in general, so it has no persistence module of its own without a constraint.

How often closure was already optimal: link 1,882 of 1,920 cells, pair 1,810 of 1,920. Most cells are crisp, so these rates are inflated, and graphs this small say nothing about ISOLET. The largest share of the loss a search recovered was 18% for link and 62% for pair.

**The window sits inside the bar.** Fix a block B with G[B] connected, and let W_B(α) be the temperatures where B is exactly one FK block with probability ≥ 1 − α. For θ ≥ ½ and α = 1 − θ, every T in W_B has B as a component of the θ-graph, ties aside:
- P(B a block) ≥ θ gives τ_e ≥ θ on B's internal edges, and τ_e ≤ 1 − θ ≤ θ on its boundary edges.
- So W_B lies inside B's bar in the T-stack, and bar minus window is where closure asserts a cluster the draws don't show whole. That is a per-cluster chaining readout.
- *Enumerated:* 9 cases, 0 violations. That is thin, so the proof carries it. The low count is itself a warning: exact-block events are rare.

**The window's ends are monotone (derived, not run).** Exact blocks are too rigid for real clusters, whose periphery flickers. The soft events keep the structure:
- "A cluster holds at least (1 − ε)|B| of B" is increasing in the bonds, so by C17's comparison its probability rises as T falls. It sets the hot end.
- "No cluster holds that much of B and more than ε|B| outside" is decreasing, so its probability falls as T falls. It sets the cold end.

Each end is a monotone curve, which BARS can fit with a band, giving a lineage's bar endpoints credible intervals. The exact containment above doesn't survive softening except on B's core.

A caveat on anything read through BARS here, these curves and ⟨E⟩(T) alike. One shared σ² against heteroscedastic SW noise makes BARS spend knots on noise and turns location intervals overconfident (primer §13.8). An event probability's noise is largest near ½, close to the crossings the window reads. Whether BARS's observation families can carry it is unchecked.

## The Persistent Homology thread's overlaps (forwarded by the user, 2026-11-01)

1. **θ as a cost ratio, against the θ-graph's components being a closure.** Resolved by the first two results: the ratio holds for the graph, the components are a closure, and Δ is the distance between the optimum's bound and the closure. m_ij is the H0 filtration's merge height, so the integrated gap is where the two readings meet.
2. **The robust window against H0 bar length.** Containment, for a fixed block, as above. The window as first proposed (`opus-gvdv-ideas.md`: the largest interval on which the optimal partition stays in its credible ball) is the wrong grain for a T-stack: clusters live at different temperatures, so a whole-partition window is empty where the resolvers are needed. The per-block window is primary; r(T) stays a per-T curve.
3. **Ultrametricity on sampled states.** Same data as VI cuts, credible balls and P(Q): kept partitions across replicas at one T, plus a distance between partitions. Two tests share the name:
   - *On states.* A triplet index over sampled partitions. Thread 1 notes that triplets inside one state are equilateral and uninformative. Ferromagnetic SPC is expected to have one state over most of T and two at a melt, which gives no hierarchy to test, so P(Q) gates it. Its natural home is frustrated SPC (Thread 4). Expected, not measured.
   - *On the field.* τ_ij against m_ij, the integrated gap above. Non-vacuous at every T.

## Built

**The closure gap (2026-11-01, `76e5dec`).**
- `ClosureGap` (`src/clustering/dendrogram/`): `OverEdges` and `OverPairs` at one θ on the co-membership scale, each taking an optional independent estimate to read the costs from, and `IntegratedOverEdges` and `IntegratedOverPairs` for ½Σ(m − τ)².
- `DenseCoMembership` (same folder): every pair's co-membership from kept partitions, N(N−1)/2 counts.
- `SweepClosureGap` (`src/clustering/graphical/spc/profiling/`): curves over a sweep's temperatures, on the cut `SweepDenseTStack` makes, with θ taken on the δ̄ scale and moved by (qθ − 1)/(q − 1).
- The curve carries the plug-in result from everything pooled, and a cross-fit gap with alternate replicas deciding and evaluating, averaged both ways.

**Why the cross-fit.** The question is whether the gap is small, and the plug-in form can't say so near a cut: a pair near θ whose noise lands it on the wrong side always adds. With K₆ and every true co-membership on the cut, the true gap is zero. The cross-fit gap evaluates each decision on data that didn't make it, so its mean sits at zero, while the plug-in mean sits near 0.26 (`Noise_InflatesTheSameEstimate_AndIndependentEvaluationCancelsIt`). The cross-fit gap measures decisions made on half the data and can come out negative.

**Tests** (`tests/VizCore.Tests/ClosureGapTests.cs`). Ground truth from outside the estimator, tolerances fixed before the first run and listed in the file's header:
- Co-membership by enumerating every bond configuration of the FK measure, the best partition by searching every partition, and the join level by a widest-path recursion.
- Six fixtures, including the chain and the poles, at q = 2 and 20, 12 values of β and 9 of θ. `Fixtures_ReachEveryKindOfContradiction` keeps the sandwich checks from passing vacuously; it failed on the first, coarser grid (8 and 3), which is why the grid is this wide.
- Swendsen–Wang against the enumeration, so the measure the checks assume is the one SW samples; `DenseCoMembership` against the accumulated edge counts, exact; and the curve against `SweepDenseTStack`'s own levels.

**Cost.** Filling `DenseCoMembership` costs N²/2 increments per kept partition under a giant cluster, about 5 × 10¹⁰ per temperature at N = 7,800 and 1,600 kept partitions. ISOLET's coldest rows took about a minute each, sampling included.

## Measurement of the closure gap

Asked for by the user on 2026-11-01. The protocol and readings below were written before the first run. The thresholds are judgement calls, fixed so the reading isn't fitted to the results.

**The question.** Does closure cost anything on real data at the levels that matter, and is that cost chaining?

**Protocol.** `tests/VizCore.Tests/ClosureGapHarness.cs`, built from an export of one commit, so the graph compiler's in-flight edits stay out of it.

| Set | N | Graph (mutual kNN + MST repair, Gaussian mean-edge kernel, 1/K̂ gauge) | Labels |
|---|---|---|---|
| `blobs` | 400 | k = 10: two Gaussians, σ = 1, centres 24 apart | The two blobs |
| `dumbbell` | 424 | k = 10: the same points, plus 24 evenly spaced on the line between them | Blob by side |
| `iris` | 150 | k = 5 | Species |
| `toy` | 3,200 | k = 10, `Bwd1995Toy` seed 42 | The generator's |
| `landsat` | 6,435 | k = 5 | Land cover |
| `isolet` | 7,797 | k = 10, raw 617 features | Letter |

- q = 20, and θ = ½ on δ̄, which is 0.474 on the co-membership scale.
- 65 temperatures evenly spaced in ln T over `BwdPottsCriticalEstimate.Bracket(20)`.
- 8 replicas × (1,000 burn-in, 2,000 cycles), keeping every 10th partition: 1,600 per temperature, pooled.
- One seed. This is exploration, so nothing is asserted.
- Each row goes to `artifacts/closure-gap/closure-gap.tsv` as it finishes, and finished rows are skipped on a rerun.

Below, *share* is the cross-fit pair gap over the closure's pair loss. A negative cross-fit gap counts as a share of zero.

**The instrument check.** There must be a temperature where `dumbbell`'s share is at least 20% and `blobs`' is at most 5%. Without one the instrument can't see a chain, and the readings below are void. Expected on the same levels, not required: `dumbbell`'s edge form stays under 5%.

**The readings, per real set.**
- *Best levels:* the temperatures whose closure has an ARI within 0.02 of the set's highest. Labels choose only where to read.
- *No room:* share ≤ 5% at every best level. A loss-optimal search isn't worth building on that set.
- *Room:* otherwise. This is necessary, not sufficient: the gap is an upper bound, and it doesn't show that a search would gain it.
- *Chaining shows:* room, and at those levels the labels side with the pair rule on more than half of the contradicted pairs.

Expected everywhere, `blobs` included: the gap peaks at melts, where a partition can't represent intermediate co-membership. That is indecision, not chaining, so `blobs` is the reference for what a plain melt costs.

## What the measurement found

Run 2026-11-01 to 11-02 from `4ed2867`, about two hours in all, at low priority on a shared machine. The rows are in `artifacts/closure-gap/closure-gap.tsv`, which git ignores.

**The instrument check passes,** at 6 temperatures from T = 0.0045 to 0.0063: `dumbbell`'s share is 25–79% with `blobs` at 0%. The expectation on the edge form missed: it stayed under 5% at only 2 of the 6, rising to 43%. The cause isn't inspected. A guess is that the bridge isn't a pure chain, since k = 10 gives it edges spanning several hops.

**The readings, as pre-registered.**

| Set | Best ARI | Best levels | Highest share there | Edge share there | Labels side with the rule | Reading |
|---|---|---|---|---|---|---|
| `iris` | 0.638 | 1, T = 0.0490 | 9.6% | 0.4% | 0% | Room |
| `toy` | 0.707 | 13, T = 0.0215–0.0490 | 4.4% | 1.6% | 71% | No room |
| `landsat` | 0.248 | 1, T = 0.0791 | 67.7% | 0.0% | 73% | Room, chaining shows |
| `isolet` | 0.148 | 3, T = 0.0907–0.1041 | 3.0% | 0.4% | 88% | No room |

**A flaw in the chaining reading.** Its threshold of one half ignores the class base rate. Pairs from different classes are 81% of Landsat's pairs and 96% of ISOLET's, so a contradicted pair drawn blindly sides with the rule about that often, and Landsat's 73% is below its base rate. The reading stands as written because it was fixed in advance, but it says little alone. The fairer comparison is the pair rule's precision against closure's at the same level. At Landsat's best level that is 0.809 against 0.408, so chaining holds up there. ISOLET's 0.96–0.99 at T ≤ 0.037 comes with precision at most 0.053 on either side, which is the base rate and no evidence.

**After the fact.** None of the rest was pre-registered. It is arithmetic on the same rows, with no new sampling.

*Peaks of the same height, of two kinds.* `toy` reads 42.3% at T = 0.0143 and 43.2% at T = 0.0847.
- On the cold peak the labels side with the rule on 97%, against a 75% base rate, and precision rises from 0.263 to 0.464. That is chaining.
- On the hot peak they side with it on 17%, precision is 0.829 either way, and only recall moves, 0.369 to 0.204. That is a melt, where pair loss would pay to part pairs that belong together. `blobs` shows only this kind, with precision 1.000 throughout.

So a large gap doesn't say which kind it is, and a search reading the same loss can't tell either.

*Chain peaks sit just before closure's own split.* The toy's pair rule reaches ARI 0.68 at T = 0.0153 and closure reaches 0.68 at T = 0.0164, one step of this grid, or 7% in T. `iris` has the same shape: at T = 0.0458 closure's pair precision is 0.630 and the rule's 0.937, and one step hotter closure reaches 0.956 by itself.

*Across the grid, the pair rule beats closure only on Landsat.* The rule is not a partition, so a search isn't sure to reach it, and this comparison favours the rule. ARI is by pair counts, to put both on one scale. Recall is the best each reaches over T at a given precision.

| Set | Best ARI, closure | Best ARI, pair rule | Rule's recall over closure's, at matched precision |
|---|---|---|---|
| `iris` | 0.638 | 0.577 | 0.90 for precision 0.69–0.91 |
| `toy` | 0.707 | 0.706 | 1.00 for 0.52–0.72 |
| `landsat` | 0.248 | 0.289 | 1.47 for 0.45–0.72, 1.16 at 0.80, 0.87 below 0.37; higher only where recall is under 0.05 |
| `isolet` | 0.149 | 0.149 | 0.76–0.90 |

*Landsat's best level is fragile.* Over T = 0.0739–0.0907 closure's ARI by pair counts runs 0.083, 0.248, 0.084, 0.211, while the pair rule's falls smoothly: 0.289, 0.241, 0.186, 0.130. Closure's joined pairs rise from 3.44 to 3.85 million between T = 0.0791 and 0.0847, though the exact closure can only lose pairs as T rises (C17). So this is noise at the cut, probably a few bridge edges, not inspected. The level the reading rests on may not repeat under another seed. The chaining signature itself covers T = 0.069–0.085.

*Smaller points.*
- The edge form stays at or under 2.6% at every level of every real set, against a pair form up to 79%. Link loss can't see chains on data either.
- False splits do occur on data: ISOLET at 11 levels, up to 2,681 pairs, and `dumbbell` at 3. Their cost is under 0.01% of the false merges'.
- The cross-fit changed the reading in places: Landsat at T = 0.0847 from 83.6% to 60.4%, the toy at T = 0.0644 from 16.0% to 8.1%.
- ISOLET's 10–18% at T = 0.0077–0.0095 is a ratio of small numbers. Its false-merge counts fit one or two dangling nodes under a single cluster.

**What this says (proposal, the user's call; for Landsat at K = 5, "Results of the probe" replaces it).**
- *Leave the search unbuilt for now.* On `toy` and ISOLET there is no room where it matters, and on `iris` the room at the best level points away from the labels: all 327 contradicted pairs are same-species pairs that pair loss would part. Only Landsat has a case: about 1.5 times the recall at matched precision, at a best ARI of 0.25, on one seed.
- *Keep the gap as a diagnostic,* with the limit above: alone it doesn't separate chaining from a melt. A label-free way to tell them apart is an open question.
- *Possibly why Landsat differs.* A T-stack reads each block at its own T, so a split delayed by a step or two costs nothing unless the block melts first. Landsat's lower recall at matched precision would fit that. Untested.

**Limits.** One seed and one θ. The gap is an upper bound on a search's gain, and the pair rule's scores against the labels aren't what a search would attain. Global ARI and pair counts don't see a T-stack's per-block reading. The integrated gap was recorded and not read. The graph is the one at `4ed2867`, without the graph compiler's later edits. Each set was read on one graph only, so every reading above belongs to that graph.

**The graphs, as resolved.** The protocol above named a recipe, which hides settings that could explain a result. These are as `GraphBuildResult.Resolved` reports them, read at the probe's commit; their edge counts match the measurement's rows. All are Euclidean kNN, mutual filter, `MstMin` repair, no LMP, a floor of 10⁻¹⁰, no repair edge dropped or zeroed, one component, and q = 20, with the mean-edge bandwidth over the repaired edges and gauge 1/K̂. `Auto` refinement resolves to `PathNeighbor`, which has no effect under Euclidean distance. The bracket is [0.004461, 0.3569], except `cauchy`'s [0.004778, 0.4778], whose units aren't BWD's.

| Graph | K | n, edges | Largest component before repair | Repair edges | Edges the floor cut | Kernel | Bandwidth | Gauge 1/K̂ (K̂) |
|---|---|---|---|---|---|---|---|---|
| `blobs` | 10 | 400, 1,503 | 50% | 1 | 0 | Gaussian | 0.3216 | 0.1331 (7.52) |
| `dumbbell` | 10 | 424, 1,574 | 53% | 1 | 22 | Gaussian | 0.4019 | 0.1347 (7.42) |
| `iris` | 5 | 150, 262 | 53% | 15 | 0 | Gaussian | 0.3065 | 0.2863 (3.49) |
| `toy` | 10 | 3,200, 12,655 | 99% | 12 | 0 | Gaussian | 0.1622 | 0.1264 (7.91) |
| `landsat`, arm `k05` | 5 | 6,435, 9,048 | 78% | 983 | 0 | Gaussian | 22.71 | 0.3556 (2.81) |
| `isolet` | 10 | 7,797, 19,454 | 91% | 297 | 0 | Gaussian | 7.365 | 0.2004 (4.99) |
| Arm `k10` | 10 | 6,435, 17,109 | 95% | 282 | 0 | Gaussian | 23.84 | 0.1881 (5.32) |
| Arm `k20` | 20 | 6,435, 34,234 | 99% | 52 | 0 | Gaussian | 25.27 | 0.0940 (10.64) |
| Arm `cauchy` | 5 | 6,435, 9,048 | 78% | 983 | 0 | Cauchy | 22.71 | 0.3556 (2.81) |

What these settings could explain:
- `dumbbell`, the instrument's positive case, is two pieces joined by one repair edge, and is the only graph where the floor cut edges. So the chain the check saw isn't a clean path of kNN edges.
- `iris` at K = 5 and `landsat` at K = 5 have the most repair (53% and 78% before it) and the lowest K̂, the conditions most likely to chain. They are the two sets that read *room*. `toy` and `isolet` read *no room* with K̂ of 7.91 and 4.99. Two sets each isn't evidence of a cause, but it is why the probe varies K.

## Probe of Landsat's case

Asked for by the user on 2026-11-02 ("Probe Landsat first"). The protocol and readings below were written before the first run. The thresholds are judgement calls, fixed so the reading isn't fitted to the results.

**The question.** The measurement left Landsat as the one case for a search, on one seed and one graph. Can a partition reach the gap there? Does what it reaches recover the land-cover labels better than closure? Does the best level hold on other seeds? And is any of this the graph's rather than closure's?

**Why the graph is in the question.** Default settings aren't correct settings, and these settings could explain the result. Landsat's graph at K = 5 has a mean degree of 2.81, with 983 of its 9,048 edges from repair and 78% of its points in the largest component before repair. That is close to a tree, where τ_ij is a product along the path, and closure's joins reach furthest past what the pairs hold. `docs/Clustering/SpcAnalysis.md` already records some of Landsat's readings as K = 5's. This chaining may be too.

**Arms.** One seed per arm except the measurement's own graph, whose seeds carry the best-level question.

| Arm | Graph | Master seeds |
|---|---|---|
| `k05` | K = 5, Gaussian: the measurement's graph | 20261102 (the measurement's), 20261103, 20261104 |
| `k10` | K = 10, Gaussian | 20261102 |
| `k20` | K = 20, Gaussian | 20261102 |
| `cauchy` | K = 5, Cauchy | 20261102 |

Each arm's graph is in the "as resolved" table above. Fixed across arms: those settings other than K and the kernel, all 6,435 points and 36 raw features, q = 20, θ = ½ on δ̄ (0.474 on the co-membership scale), and 8 replicas × (1,000 burn-in, 2,000 cycles), keeping every 10th partition.

**Grid.** 65 temperatures evenly spaced in ln T over `AutoBracket`'s range for that graph, worked from the hot end down. An arm stops going colder once closure's ARI has reached 0.05 and then stayed under 0.02 for three levels in a row. Labels choose only where to stop.

**The search.** `PairLossSearch` (`tests/VizCore.Tests/`), left with the tests since nothing has decided it should be more. It maximizes Σ_joined (τ_ij − θ), which is the pair loss up to a constant. Seven starts, each polished by single-node moves and positive merges, keeping the best score on the deciding field: the closure, agglomeration by summed and by mean weight, and four sequential allocations in random order. What it finds is a lower bound on what a partition can reach.

**Measures, per level.**
- *Saving:* the search's fall in pair loss below closure. *Cross-fit saving:* each replica half decides both partitions and the other half evaluates, averaged both ways, as the cross-fit gap does.
- *Attained fraction:* cross-fit saving over cross-fit gap, at levels where the gap's share is over 5%.
- *Attained share:* cross-fit saving over closure's pair loss.
- *ARI:* the pooled partition of each, closure and search, against the labels.

**The readings.**
- *Reproduction.* `k05` on its original seed must match the measurement's `ari` and pair loss exactly on every level they share. If not, comparison with the measurement is void, and the rest stands alone.
- *The best level.* It holds if closure's ARI at T = 0.0791 is within 0.05 of 0.248 on both fresh seeds of `k05`. Otherwise it doesn't.
- *Attainable.* The median attained fraction over `k05`'s levels T = 0.0690–0.0847 and its three seeds. *Attainable* at ½ or more, *little* at ¼ or less, *partly* between.
- *Better against the labels.* Gain is the search's best ARI over an arm's levels minus closure's best, for each seed. On `k05`: *worth building on this graph* if the gain is at least 0.05 on at least two of three seeds and the attained share is at least 5% at the search's best level; *no gain* if the gain is at most 0.02 on at least two of three; *small* otherwise. On the other arms, the same bands on their one seed.
- *The graph.* Each arm gets the measurement's *room* reading at its own best levels. *The chain is K = 5's* if `k10` and `k20` both read no room. *A graph is the cheaper lever* if some other arm's closure reaches a best ARI at or above the median over `k05`'s seeds of the search's best.

Reported with no threshold: ARI's range across seeds at each `k05` level, for closure and search; pairs joined at a hotter level and not at the next colder one, which breaks nesting, for both; recall at matched precision; cluster counts; run time. Partitions are kept in `artifacts/closure-gap/landsat-probe-labels.tsv` for these and for later reading.

### Results of the probe

Run 2026-11-02, 7.5 hours in all at low priority on a shared machine, every arm finishing. The rows are in `artifacts/closure-gap/landsat-probe.tsv` and the partitions beside them, which git ignores. Every finding below belongs to the graphs in the "as resolved" table, at θ = ½ and this budget.

**The readings, as pre-registered.**

| Reading | Result | Numbers |
|---|---|---|
| Reproduction | Holds | `k05` on its original seed matches the measurement exactly on all 43 levels: `ari`, pair loss, rule loss, cross gap, false merges, clusters, joined pairs |
| The best level | Doesn't hold | Closure's ARI at T = 0.0791 is 0.091 and 0.090 on the fresh seeds, against 0.248 |
| Attainable | Attainable | Median attained fraction 0.992 over 12 cells, range 0.919–0.995 |
| Better against the labels, `k05` | **Worth building on this graph** | Gain on all three seeds, attained share 48–78% at the search's best level |
| Better against the labels, `cauchy` | **Worth building on this graph** | One seed, attained share 65% |
| Better against the labels, `k10`, `k20` | No gain | One seed each |
| The chain is K = 5's | Not met | `k10` reads room at two of its three best levels, `k20` at none |
| A graph is the cheaper lever | Not met | The median of the search's best on `k05` is 0.317, and no other arm's closure reaches it |

| Arm | Master seed | Closure's best ARI (T) | Search's best ARI (T) | Gain | Highest share at closure's best levels |
|---|---|---|---|---|---|
| `k05` | 20261102 | 0.248 (0.0791) | 0.317 (0.0690) | +0.069 | 67.7% |
| `k05` | 20261103 | 0.226 (0.0847) | 0.284 (0.0739) | +0.058 | 63.7% |
| `k05` | 20261104 | 0.221 (0.0907) | 0.321 (0.0690) | +0.099 | 56.6% |
| `cauchy` | 20261102 | 0.259 (0.0736) | 0.312 (0.0593) | +0.053 | 37.0% |
| `k10` | 20261102 | 0.236 (0.0907) | 0.244 (0.0972) | +0.008 | 17.1% |
| `k20` | 20261102 | 0.277 (0.0972) | 0.281 (0.0791) | +0.005 | 1.5% |

Search's best level isn't near where any arm stopped: it sits 11–19 levels from the cold end of what ran, and the search's ARI falls well below its peak on the way there.

**With no threshold.**
- *ARI across `k05`'s seeds.* Away from the flips below, the range at any level with an ARI above 0.03 is at most 0.013, for both, apart from the search's 0.045 at T = 0.0231, where the arms end. Closure's is 0.159 at T = 0.0791 and 0.143 at T = 0.0847. The search's is 0.115 at T = 0.0690, and 0.002 at T = 0.0739, where it reads 0.281–0.284. So the search isn't steadier in general. Each flips at its own cliff, and the search's cliff is colder.
- *At the search's steady level.* Taking T = 0.0739 in place of each seed's peak, the gains are +0.034, +0.058 and +0.060, still at least 0.05 on two of three.
- *Nesting.* Closure breaks none on `k10` and `k20`, and 13,129 pairs on `cauchy`, at 4 of 40 levels. On `k05`, 16,980 and 26,995 at 6 and 7 of 42, and 1,590,691 on the original seed at 3 of 42, nearly all of it one flip at T = 0.0847, 41% of that level's joined pairs. The search breaks more levels at smaller sizes: 30–31 of 42 on `k05` (33,785–141,051 pairs), 23 of 40 on `cauchy` (91,319), 15 of 31 on `k10` (27,208), 10 of 34 on `k20` (3,749). Its worst level is 8.6% of the joined pairs, on `k05`. The search's partitions don't form a hierarchy, so a T-stack can't read them as they stand.
- *Recall at matched precision, search over closure.* On `k05`, 1.41–1.90 for precision 0.50–0.75 on all three seeds, 0.67–1.10 at 0.80, and 0.82–1.15 at 0.30–0.40. On `k10` and `k20`, 0.77–1.11. On `cauchy`, 1.32 for 0.45–0.70. The search reaches about the pair rule's recall on `k05`, within 0.02 for precision 0.65–0.75, so the rule's advantage in the measurement was attainable by a partition.
- *Clusters.* At T = 0.0739 on `k05`, closure has 1,652–1,674 and the search 2,347–2,360. Most are singletons on either side, 1,515–1,527 and 1,958–1,981.
- *Run time.* Per arm, for three searches a level: search 1,827–5,616 s, sampling 143–468 s, the dense field and gap 311–2,106 s. One search has a median of 10–42 s by arm, and the longest took 255 s.
- *Shares over 100%.* `cauchy` reads 111% at T = 0.0685. The cross gap uses each half's closure and the denominator the pooled one's, so a level where a half's closure chains and the pooled one doesn't can exceed 1.
- *Attained fraction elsewhere.* It stays at 0.83–1.12 on `k05`'s other levels and at 0.84–1.00 on `cauchy`'s, and falls where the gap is small: 0.55–1.00 on `k10`, and 0.38 at `k20`'s one level over 5%.

**After the fact.** None of this was pre-registered. It is arithmetic on the stored partitions, with no new sampling. Purity is the majority class's share over the clusters read, the definition of `Clustering.Evaluation.External.Purity` and `BwdParityTests.BigClusterPurity`, computed here by the reader. Those tests read clusters of at least 50. Leading with 20 below was my choice, so both sizes are given.

*What the cliff is.* On `k05`, seed 20261103, closure at T = 0.0907 holds class 1 (713 points, purity 1.00), class 3 (680, 0.82), class 7 (678, 0.82) and class 2 (133). By T = 0.0847, 3 and 7 have joined with half of class 4 (1,832, 0.42). By T = 0.0791 class 1 has joined too (3,360, 0.29). An ARI near 0.09 against 0.22 is that one merge, class 1 in the giant or out of it: over T = 0.0739–0.0907 and all three seeds, closure reads 0.083–0.093 with class 1 in and 0.209–0.248 with it out, with no exception. The fragility across seeds is which side of the merge a level falls on. The original seed has it out of order, in at T = 0.0847 and out at T = 0.0791. Class 4 splits between 3's and 7's clusters before they join, which would fit a path through 4. That isn't inspected.

*What the search does to it.* It keeps the same blocks apart further down. At T = 0.0739, on all three seeds, it reads class 3 (934–935, 0.78), class 7 (823–834, 0.81), class 1 (804–821, 1.00), class 2 (158) and class 5 (103–123), where closure has a giant of about 4,000 at 0.27. It reaches its own giant at T = 0.0644, two to three steps of this grid below closure's, or 14–22% in T. At T = 0.0690 two seeds read 3 and 7 apart and one reads them joined, which is the 0.206–0.321 flip. The attained fraction near 1 is dominated by cutting these giants, and says little about the fine structure.

*The same search parts what belongs together on the hot side.* At T = 0.0907 closure's class 1 of 713 points is the search's 371 and 258. At equal T the search's ARI is below closure's at every level from T = 0.0907 up to where both read 0 (0.122 against 0.21 at T = 0.0907). This is the melt kind of gap from the measurement: pair loss pays for that split too, and can't tell the two apart. The search does better only at a colder T than closure's best.

*The points covered at a given purity,* over each arm's levels, in clusters of at least 20 points:

| Arm | Purity ≥ 0.90: closure, search | Purity ≥ 0.85: closure, search | Purity ≥ 0.80: closure, search |
|---|---|---|---|
| `k05`, three seeds | 2,254–2,656, 2,481–2,784 | 2,648–2,658, 3,150–3,171 | 2,648–2,658, 3,150–3,171 |
| `cauchy` | 2,065, 2,633 | 2,957, 2,976 | 2,957, 3,439 |
| `k10` | 1,667, 1,562 | 2,118, 1,980 | 2,118, 1,980 |
| `k20` | 2,000, 1,938 | 2,409, 2,341 | 2,409, 2,341 |

So on `k05` the search covers 19% more points at purity 0.85, on all three seeds, and 19–22% at a minimum size of 50 (2,399–2,469 against 2,921–2,926). At 0.90 it depends on the seed, 0.93–1.24 times closure.

*K = 20's ARI is a different partition.* Its best of 0.277 reads 3 and 7 joined, with 4 (1,933 points: 783, 757, 362), at purity 0.664 over 3,429 points. By points covered at purity 0.85, the order changes: `k05`'s closure 2,650, then `k20`'s 2,409, then `k10`'s 2,118. So on the flat reading a higher K doesn't buy what the search buys, and 0.277 and 0.317 aren't to be read as close. This is one seed, K up to 20, and the flat reading only. `docs/Clustering/SpcAnalysis.md` records the T-stack's largest-level reading at K = 10 and 20 at median purity 0.940 and 0.961, so it says nothing about K for the reading SPC uses.

**What the settings could explain.**
- *K.* The share at the best levels falls with K, 57–68%, 17%, 1.5%, and the gain goes with it. The gain is the sparse graph's: closure pays for K̂ = 2.81 and the search pays it back. `k10` still reads room, so it isn't K = 5 alone.
- *The kernel.* `cauchy` behaves as `k05` does, on one seed, with closure alone covering 11% more at purity 0.85.
- *θ, the most open.* The probe held it at ½. A higher θ on closure may cut the same giants, with no dense field and no search. `docs/Clustering/SpcAnalysis.md` records θ = 0.6 passing on Landsat where θ = 0.5 reads a giant, on a T-stack at K = 5 and q = 5, which is another q and another reading. The levers tested were K and the kernel only, so *a graph is the cheaper lever* not being met says nothing about θ.
- *Repair.* K = 5 has 983 repair edges, which are candidates for what bridges the blocks. The rows don't keep τ on edges, so this isn't inspected. (Read since, in "Results of the θ test": they aren't what bridges them.)
- *Fixed and not varied:* Euclidean distance on the 36 raw features, the mutual filter, no LMP, mean-edge bandwidth, gauge 1/K̂, q = 20, the budget.

**What this says (proposal, the user's call).**
- *The measurement's proposal no longer stands for K = 5.* It rested on one seed and asked whether a partition could reach the gap. Now three seeds show one can, and it recovers the labels better: 0.03–0.10 in ARI, and 19% more points at purity 0.85. The reading, fixed in advance, is worth building on this graph.
- *Test the cheaper closure before building.* The search costs a dense N²/2 field, kept partitions, half a minute or so a level, and a lost hierarchy. Closure at a higher θ needs none of these, and the probe didn't try it. If it covers about 3,150 points at purity 0.85 on `k05` (about 2,920 in clusters of at least 50), the search isn't needed. It may well not: raising θ at one T may only retrace raising T at one θ, which closure has already read. Then the case for building holds, with that as its evidence. (Tested since, in "Results of the θ test": it doesn't reach it.)
- *The flat cut isn't SPC's reading.* The T-stack on this K = 5 graph covers 2,368–2,610 points at purity at least 0.90 at the higher budget (`docs/Clustering/SpcAnalysis.md`), in clusters of at least 50. At that size and purity, closure's flat cut covers 2,084–2,469 and the search's 2,163–2,550, at this budget, so the T-stack's reading is at least level with both. Whether the search adds anything to a per-block reading is unmeasured, and needs an answer on nesting first.

**Limits.** Labels choose the level on both sides, and no label-free way to pick the search's T is tested. Three seeds on `k05` and one on every other arm. One θ, one budget. The search is a lower bound on what a partition can reach. ARI and purity are global. The probe's binaries predate `8af6199`, which landed during the run. Rebuilt at `fca555a`, all ten graphs resolve the same, plus its new `RestoredEdges = 0`, and `ClosureGapTests` pass, 43 of 43.

## Closure at a higher θ

Asked for by the user on 2026-11-02 ("Test higher θ first"). The protocol and readings below were written before the first run. The thresholds are judgement calls, fixed so the reading isn't fitted to the results.

**The question.** The probe varied K and the kernel and held θ at ½. Can closure alone, at a higher θ, reach what the search reached on Landsat at K = 5: about 2,920 points at purity 0.85 in clusters of at least 50, and a best ARI of 0.28–0.32? If so, the search isn't needed there, and its costs go with it.

**Two predictions, either of which could hold.** Two frozen blocks cross θ at T = J/ln(θ(q−1)/(1−θ)) (`docs/Clustering/SpcAnalysis.md`). So going from θ = 0.5 to 0.9 at q = 20 lets a split survive to 0.57 of the temperature, where the blocks hold together better, which predicts a gain. Against that, raising θ at one T may only retrace raising T at one θ, which closure has already read, and that predicts none.

**Arms.** The probe's six: `k05` on master seeds 20261102, 20261103 and 20261104, and `k10`, `k20` and `cauchy` on 20261102. Their graphs are in the "as resolved" table. They are built from the current tree, which includes `8af6199`; those graphs resolve the same, and the reproduction reading below tests it.

**What runs.** The probe's own chains, from the same seeds at the same budget, 8 replicas × (1,000 burn-in, 2,000 cycles). All 65 temperatures of each arm's grid run, with no stopping rule, since a higher θ's best T may sit colder than where the probe stopped. Closure is `AffinityThreshold.Connect` on the edge column over all draws, at each θ on the δ̄ scale in 0.50, 0.55, …, 0.95, 0.99. Each level's edge column is kept in `artifacts/closure-gap/landsat-theta-columns/`, with the edge list and which edges came from repair, so later questions about edges or θ need no new sampling.

**Measures, per arm and θ, over T.**
- *Best ARI:* the highest ARI against the labels.
- *Points:* the most points covered at purity at least 0.85, in clusters of at least 50. Purity is `Clustering.Evaluation.External.Purity` with smaller clusters unassigned, which equals `BwdParityTests.BigClusterPurity`. A size of 50 is what the parity tests and `SpcAnalysis`'s Landsat table read, so it replaces the 20 the probe's results led with. Clusters of at least 20 are reported beside it.
- *The search's numbers,* from the probe on the same seed. Best ARI: 0.317, 0.284, 0.321 on `k05`, and 0.312 on `cauchy`. Points: 2,921, 2,923, 2,926 on `k05`, and 2,737 on `cauchy`.

**The readings.**
- *Reproduction.* At θ = ½, each arm's ARI must match the probe's `ari` exactly on every level they share. If not, comparison with the probe is void, and the rest stands alone.
- *Does θ help closure?* A seed *gains* at a θ if, against θ = ½ on that seed, its best ARI rises by at least 0.05 or its points by at least 10%. On `k05`, θ *helps* if one θ gains on at least two of three seeds. Otherwise it doesn't.
- *Does it replace the search?* A seed *reaches* the search at a θ if its best ARI is no more than 0.02 below the search's on that seed and its points are at least 95% of the search's. It *falls short* if its best ARI is more than 0.05 below, or its points are under 90%. On `k05`: *search not needed on this graph* if one θ reaches on at least two of three seeds; *search still ahead* if at every θ at least two seeds fall short; *in between* otherwise.
- *Other arms.* On `cauchy`, both readings, on its one seed. On `k10` and `k20`, whether θ helps, on their one seed. They have no search target, since the search showed no gain there.

These favour closure in two ways, both accepted so a *search still ahead* means something. Closure gets to choose θ, while the search ran only at ½. And closure reads all 16,000 draws per level, while the search read every 10th.

Reported with no threshold, by θ: the T of the best ARI; ARI's range across `k05`'s seeds at each T, which says whether the flip persists; pairs breaking nesting across T; recall at matched precision, against the search; cluster counts and singletons; points at purity 0.80 and 0.90; T × ln(θ(q−1)/(1−θ)) at the level where class 1 joins the giant, which tests the first prediction. And from the kept columns: at the level before class 1 joins, the edges across its block's boundary, how many there are, their δ̄, and how many came from repair. That is the note's open question on what bridges the blocks.

### Results of the θ test

Run 2026-11-02 from `b664d7d`, 46 minutes in all at low priority on a shared machine, every arm finishing. The rows are in `artifacts/closure-gap/landsat-theta.tsv` and the columns in `landsat-theta-columns/`, which git ignores. Every finding belongs to the graphs in the "as resolved" table, at this budget.

**The readings, as pre-registered.**

| Reading | Result | Numbers |
|---|---|---|
| Reproduction | Holds | ARI at θ = ½ matches the probe's `ari` exactly on every shared level of every arm (32–43 levels each), so `8af6199` changed nothing here |
| Does θ help closure, `k05` | Doesn't | No seed gains at any θ. Nearest: +0.049 in ARI (20261103, θ = 0.65) and 1.09 times the points (20261104, θ = 0.60) |
| Does θ help closure, `cauchy`, `k10`, `k20` | Doesn't | No θ gains. Best ARI rises at most 0.012, and points never rise |
| Does it replace the search, `k05` | **Search still ahead** | All three seeds fall short at every θ. Nearest: 0.026 below in ARI with 82% of the points (20261103, θ = 0.65); 92% of the points, 0.081 below in ARI (20261104, θ = 0.60) |
| Does it replace the search, `cauchy` | Search still ahead | Short at every θ. Nearest: 0.057 below in ARI with 89% of the points (θ = 0.70) |

On `k05`, across its three seeds:

| θ | Best ARI | At T | Points at purity ≥ 0.85, clusters ≥ 50 | Over the search's |
|---|---|---|---|---|
| 0.50 | 0.209–0.248 | 0.0791–0.0907 | 2,392–2,460 | 0.82–0.84 |
| 0.55 | 0.213–0.231 | 0.0791–0.0847 | 2,143–2,562 | 0.73–0.88 |
| 0.60 | 0.204–0.240 | 0.0739–0.0791 | 2,270–2,689 | 0.78–0.92 |
| 0.65 | 0.246–0.257 | 0.0739 | 2,408 | 0.82 |
| 0.70 | 0.216–0.221 | 0.0739 | 2,509–2,540 | 0.86–0.87 |
| 0.75 | 0.195–0.196 | 0.0739 | 2,385–2,412 | 0.82–0.83 |
| 0.80 | 0.199–0.200 | 0.0690 | 2,294–2,319 | 0.78–0.79 |
| 0.85 | 0.197–0.205 | 0.0690 | 2,065–2,462 | 0.71–0.84 |
| 0.90 | 0.192–0.197 | 0.0690 | 2,327–2,383 | 0.80–0.81 |
| 0.95 | 0.174–0.183 | 0.0644 | 1,960–2,064 | 0.67–0.71 |
| 0.99 | 0.136–0.146 | 0.0525–0.0602 | 1,693–1,760 | 0.58–0.60 |
| The search, θ = ½ | 0.284–0.321 | 0.0690–0.0739 | 2,921–2,926 | 1 |

**The two predictions.** The second holds and the first doesn't. Up to θ = 0.90 the best reading moves colder and its points barely change, so raising θ retraces raising T. Under the frozen-block form, T × ln(θ(q−1)/(1−θ)) at class 1's join would be constant. It rises from 0.218–0.233 at θ = ½ to 0.309 at θ = 0.90 and 0.369–0.396 at θ = 0.99. At θ = 0.90 class 1 joins at T ≈ 0.060, where the form predicts 0.042–0.045.

**With no threshold.**
- *The flip.* ARI's largest range across `k05`'s seeds at any level is 0.158 at θ = ½, at most 0.049 at θ = 0.55–0.65, 0.124 at 0.70, 0.107 at 0.75, 0.006 at 0.80, 0.047 at 0.85, 0.005 at 0.90, 0.048 at 0.95 and 0.097 at 0.99. So θ doesn't remove it. It depends on whether a grid level falls near the join.
- *Nesting.* On `k05`, 14,999–36,867 pairs at θ = ½, plus 1,586,144 on the original seed, nearly all one flip. At θ = 0.55–0.75, 7,595–111,559 at 3–9 levels; at 0.90, 84,940–184,060 at 14–15 levels; at 0.99, over 2 million at 47. Exact closure can't break nesting, so this is noise near the cut, and it grows as more edges sit near a high θ. The search's was 33,785–141,051 at 30–31 levels.
- *Recall at matched precision, `k05`.* Against closure's best over all θ and T, the search leads for precision up to 0.75 (0.221–0.292 against 0.165–0.208). Closure leads at 0.85 (0.156–0.167 against 0.083–0.105), and at 0.80 on two seeds of three.
- *Clusters at the best ARI.* 2,052–3,083 at θ = ½ and about 3,400 at θ = 0.90, mostly singletons (1,875–2,866 and about 3,250).
- *What crosses class 1's block* at the level before it joins, at θ = ½. 189–222 edges cross, 84–105 of them from repair, against repair's 11% of all edges. Their highest δ̄ is 0.499–0.500. Nearly all go to small clusters and singletons: exactly one goes to a cluster of at least 50, and it isn't a repair edge.
- *Run time.* 6–9 minutes an arm, with no search and no dense field. The columns take 4–17 MB an arm.

**After the fact.** None of this was pre-registered. It is arithmetic on the stored rows, partitions and columns, with no new sampling, and I went looking after seeing the readings, over ten θ and several purities. Labels chose which blocks the cuts below are between. Purity from the columns matches the evaluator's rows exactly on 1,430 closures.

*The search's lead depends on where it's read.* Points in clusters of at least 50, on `k05`:

| Purity ≥ | The search | Closure, θ = ½ | θ = 0.65 | θ = 0.70 |
|---|---|---|---|---|
| 0.95 | 1,025–1,034 | 641–981 | 725–1,012 | 631–950 |
| 0.90 | 2,163–2,550 | 2,083–2,460 | 2,075–2,408 | 2,509–2,540 |
| 0.85 | 2,921–2,926 | 2,392–2,460 | 2,408 | 2,509–2,540 |
| 0.80 | 2,921–2,926 | 2,392–2,460 | 2,811–2,820 | 2,509–2,540 |
| 0.75 | 2,923–3,478 | 2,392–2,460 | 2,811–2,820 | 2,509–2,540 |

So at purity 0.90, closure at θ = 0.70 covers about as much as the search or more (0.99–1.17 times, by seed) and with less spread. That is about the T-stack's 2,368–2,610 at the same purity and size, at another budget. At purity 0.80, θ = 0.65 reaches 96% of the search. The search's lead sits at purity 0.85, below 0.80, and in ARI. At T = 0.0739, θ = 0.70 holds all six classes apart at purity 0.90: 3 (about 745), 1 (about 720), 7 (about 630), 5 (119–206), 2 (130) and 4 (104). The search at that T reads more points in classes 1, 2, 3 and 7, at purity 0.87, with class 5 in two pieces and no block for class 4.

*What bridges the blocks: very few edges, and they are the graph's own.* The minimum cut between the blocks of classes 1, 3 and 7, over the whole graph at any T:

| Graph, and the level the blocks are from | 1 and 7 | 3 and 7 | 1 and 3 |
|---|---|---|---|
| `k05`, three seeds, T = 0.0907 | 3 | 6 | 3 |
| `cauchy`, T = 0.0736 | 3 | 6 | 3 |
| `k10`, T = 0.1278 | 25 | 59 | 25 |
| `k20`, T = 0.1278 | 92 | 178 | 92 |

K = 10 and 20 hold blocks that size apart only at a hotter level, so their blocks are smaller, which if anything understates the difference. Taking the repair edges out changes no cut on any graph, and no cut at any join holds a repair edge. So my guess that repair edges bridge the blocks was wrong. Over edges above θ = ½ on `k05`, the cut is 0–2 at the first levels below T = 0.0907 and stops at 2 for classes 1 and 7 and at 4 for 3 and 7, down to T = 0.043. On `k20` it climbs to 46 and 152.

At every join on `k05` at θ = ½, 0.65 and 0.90, the two blocks meet through a cut of 1 or 2 edges, and δ̄ on the widest path sits just above θ (0.50–0.56 at θ = ½). `k20`'s join of 3 and 7 goes through 15. In 9 of the 10 edges cut at `k05`'s joins of 3 and 7 both ends are class 4, which supports the path through 4 guessed above. So on K = 5 one edge can decide between an ARI of 0.09 and 0.22, and noise on that edge's δ̄ would explain the flip. This belongs to K = 5's graph: the mutual filter leaves these classes almost apart, which is also why it reads them at all.

*What the search cuts.* At T = 0.0739 it cuts 739–761 edges above θ, about 10.5% of them. Only 5–6 run between clusters of at least 50, and 1 between clusters of at least 300; at T = 0.0690 it's 2–7. So parting the big blocks costs it a handful of edges, and the rest shave nodes off the blocks' edges.

**What this says (proposal, the user's call).**
- *Closure at a higher θ doesn't replace the search,* by the reading fixed in advance, in closure's favour. The case for building holds on K = 5.
- *The cause points at something cheaper, untested.* If closure's failure on K = 5 is one or two edges, a closure that won't cross a thin cut between two large sides may reach what the search reaches, reading only the edge column: no dense field, no kept partitions. That is edge connectivity's ground, as in highly connected subgraphs (Hartuv & Shamir 2000, from memory). Being k-edge-connected is an equivalence that only coarsens as edges are added, so those classes would nest as far as closure's do, which the search's don't. It can be tried on the kept columns with no new sampling, against the same targets.
- *What could sink that.* With a mean degree of 2.81, plain k-edge connectivity strands most nodes, so it needs a size rule or peripheral capture, and those are settings fitted on one data set. The search gains by shaving too, which a cut wouldn't. And the T-stack already reads each block before it joins, so the gain over SPC's real reading may be small.
- *θ isn't the lever alone, but it isn't idle.* At purity 0.90 a flat cut at θ = 0.70 reads about what the T-stack reads at θ = ½. What the T-stack itself reads across θ on this graph isn't measured here.

**Limits.** Labels choose T on both sides and θ on closure's. The search ran only at θ = ½, so it might also move at a higher one. Three seeds on `k05` and one on every other arm. The grid steps 7% in T, and joins are sharper than that. The flat cut isn't the T-stack's reading. Minimum cuts were read between blocks chosen by label, so a label-free procedure has to find them first.

## A closure that won't cross a thin cut

Asked for by the user on 2026-11-02 ("Test thin cuts first"). The protocol and readings below were written before the first run, once `ThinCutClosureTests` passed. The thresholds are judgement calls, fixed so the reading isn't fitted to the results, and they are the θ test's wherever a reading carries over.

**The question.** The θ test found that on K = 5 the big blocks join through a cut of one or two edges. Can a closure that refuses such cuts between two large sides, reading only the edge column, reach what the search reached there: about 2,920 points at purity 0.85 in clusters of at least 50, and a best ARI of 0.28–0.32?

**The rule.** `ThinCutClosure`, which stays with the tests. At a T and a θ, take the edges closure reads, those above the cut. They build a cut tree (Gomory & Hu 1961, by Gusfield 1990's flows, both from memory), which holds a minimum cut for every pair of nodes. A cut's weight counts its edges first and sums their co-membership second, so among cuts with equally few edges the tree holds the weakest. `Partition(k, m)` then goes over the tree's cuts of at most k edges, lightest first, and removes each one whose two sides within the current piece both hold at least m points. A cluster is a piece's part connected above the cut, so it is connected, as closure's are. A cut's size is counted over all edges above the cut, not within the piece.

- k = 0 is closure.
- m = 1 is the plain form with no size rule: the connected parts of the (k + 1)-edge-connected classes, on the ground of Matula's k-components (1972, from memory). They only coarsen as edges are added.
- With m > 1 that isn't guaranteed, since a side can grow past m and part from what it had joined.

**Settings.** k in 1, 2, 3, 4, 6 and m in 1, 20, 50, 100, 200, 400, so 30 of them. For k: cuts above θ on `k05` stop at 2 and 4, and whole-graph cuts at 3–6. For m: from no rule, through the 50 the measures read, to the 300 and more that the big blocks hold.

**What runs.** No new sampling. The rule reads `landsat-theta-columns/` at all 65 levels. `k05` on its three seeds and `cauchy` run at all 11 θ. `k10` and `k20` run at θ = ½ only, since they have no search target and the densest trees. The graphs are those in the "as resolved" table, which these columns came from.

**Measures, per arm and setting (θ, k, m), over T.** As in the θ test. *Best ARI* against the labels. *Points:* the most points covered at purity at least 0.85, in clusters of at least 50, with purity from `Clustering.Evaluation.External.Purity`. *The search's numbers,* best ARI 0.317, 0.284, 0.321 on `k05` and 0.312 on `cauchy`; points 2,921, 2,923, 2,926 on `k05` and 2,737 on `cauchy`. Closure's numbers are the k = 0 rows.

**The readings.**
- *Reproduction.* The k = 0 rows must match `landsat-theta.tsv` exactly in every column but timing, on every arm, level and θ run. If not, the run is void.
- *Does the rule help closure?* A seed *gains* at a setting if, against closure at the same θ on that seed, its best ARI rises by at least 0.05 or its points by at least 10%. On `k05` the rule *helps* if one setting gains on at least two of three seeds. Otherwise it doesn't.
- *Does it replace the search?* A seed *reaches* the search at a setting if its best ARI is no more than 0.02 below the search's on that seed and its points are at least 95% of the search's. It *falls short* if its best ARI is more than 0.05 below, or its points are under 90%. On `k05`: *search not needed on this graph* if one setting reaches on at least two of three seeds; *search still ahead* if at every setting at least two seeds fall short; *in between* otherwise.
- *Other arms.* On `cauchy`, both readings on its one seed, over every setting and at each setting that reaches on `k05`, since that is the only check here that the setting wasn't chosen on it. On `k10` and `k20`, whether the rule helps at θ = ½, on their one seed.

These favour the rule. Labels choose T, it has 330 settings against the search's one, and it reads all 16,000 draws per level where the search read every 10th. So *search not needed* would say only that this family reaches the target on this data set. Judging one fixed setting takes a fresh data set, which takes new sampling and the user's go-ahead.

**Predictions, which decide nothing.**
1. *The plain form strands points.* At k = 2, m = 1, θ = ½ on `k05`, points fall under half of closure's 2,392–2,460.
2. *The blocks hold colder.* At k of at least 2 and m of at least 100, at θ = ½, the T of the best ARI on `k05` falls below closure's 0.0791–0.0907.
3. *The flip shrinks.* At k = 2, m = 100, θ = ½, ARI's largest range across `k05`'s seeds at any level is at most 0.05, against closure's 0.158.
4. *Shaving still counts.* ARI gets nearer the search than points do, since a cut doesn't take nodes off a block's edge.
5. *It is K = 5's graph.* `k10` and `k20` don't gain.

Reported with no threshold, by k and m: the T of the best ARI; clusters of at least 50 and how many classes lead one; ARI's range across `k05`'s seeds at each level; pairs breaking nesting across T; points at purity 0.80 and 0.90; time per tree.

**Limits, known in advance.** The tree holds one minimum cut for each pair, not every cut. So at k of at least 2, a thin cut between two large sides can be missed if it is the lightest between no pair. At k = 1 none can (`ThinCutClosureTests`). With a size rule the order of removal matters, and lightest first is a choice. The whole-graph cuts behind the settings were read with blocks chosen by label, so the settings aren't free of the labels.

### Results of the thin-cut test

Run 2026-11-02 from `2a92bc1`, 26 minutes in all at low priority on a shared machine, every arm finishing, with no new sampling. The rows are in `artifacts/closure-gap/landsat-thincut.tsv`, which git ignores. Every finding belongs to the graphs in the "as resolved" table.

**The readings, as pre-registered.**

| Reading | Result | Numbers |
|---|---|---|
| Reproduction | Holds | All 2,990 k = 0 rows match `landsat-theta.tsv` exactly in every column but timing |
| Does the rule help closure, `k05` | **Helps** | 263 of 330 settings gain on at least two seeds |
| Does it replace the search, `k05` | **Search not needed on this graph** | 157 settings reach on at least two seeds and 146 on all three. At θ = ½ alone, 18 of 30 reach on all three |
| `cauchy` | Helps, and reaches | 266 of 330 gain. Of `k05`'s 157, 154 reach here and 1 falls short |
| `k10` | Helps, by the letter | 14 of 30 gain: ARI up to +0.052, points up to 1.13 times |
| `k20` | Doesn't | None of 30 gain: ARI up to +0.006, points unchanged |

On `k05` at θ = ½, across its three seeds, in clusters of at least 50. The graph-alone rows are the control described below:

| Reading | Best ARI | At T | Points, purity ≥ 0.85 | ≥ 0.90 | ≥ 0.95 |
|---|---|---|---|---|---|
| Closure | 0.209–0.248 | 0.0791–0.0907 | 2,392–2,460 | 2,083–2,460 | 641–981 |
| The search | 0.284–0.321 | 0.0690–0.0739 | 2,921–2,926 | 2,163–2,550 | 1,025–1,034 |
| k = 2, m = 100 | 0.398–0.407 | 0.0739 | 4,187–4,221 | 3,050–3,066 | 641–981 |
| k = 4, m = 100 | 0.510–0.516 | 0.0399–0.0490 | 5,450–5,629 | 3,050–3,066 | 641–981 |
| k = 4, m = 200 | 0.530–0.548 | 0.0373–0.0427 | 5,141–5,450 | 2,083–2,460 | 641–981 |
| Graph alone, closure | 0.313 | | 1,618 | 1,530 | 1,057 |
| Graph alone, k = 2, m = 100 | 0.549 | | 4,701 | 2,959 | 1,057 |
| Graph alone, k = 4, m = 100 | 0.510 | | 5,500 | 3,133 | 1,057 |

**The predictions.**
1. *The plain form strands points.* Fails on its threshold, holds in substance. At k = 2, m = 1, points are 1,676–1,695, which is 69–71% of closure's and not under half, but ARI falls to 0.094 with about 4,020 clusters. At k = 3 points fall to 444, and at k = 4 to none. So no setting with m = 1 reaches, and only 1 of the 157 has m = 400: the size rule is what makes this work.
2. *The blocks hold colder.* Holds. The T of the best ARI runs 0.0739 at k = 2, 0.0690 at k = 3, 0.0373–0.0490 at k = 4, and 0.0283–0.0399 at k = 6. A larger k reads colder.
3. *The flip shrinks.* Holds. At k = 2, m = 100, ARI's largest range across seeds is 0.023, against closure's 0.158. That supports the reading that the flip was one or two edges' noise. It doesn't hold everywhere: at k = 1, m > 1 the range is 0.324–0.381, which fits one seed having one edge across where another has two, though I didn't inspect it.
4. *Shaving still counts.* Fails, since points don't lag. At k = 4, m = 100 and 200, they run 1.76–1.93 times the search's, with ARI 0.20–0.25 above it. Going colder takes in the periphery, and purity stays near 0.85 because the blocks stay apart, so the rule needs no shaving to pass this measure.
5. *It is K = 5's graph.* Half. `k20` doesn't gain, and `k10` gains on 14 settings. Their cuts of 25 and more are out of reach at k of at most 6, so this says little about larger k.

**With no threshold.**
- *Why points rise.* At a given T the rule covers about the same points as closure, since parted pieces mostly stay above the size read. What changes is purity. At T = 0.0739 on seed 20261102, closure reads 4,221 points at purity 0.311, and k = 2, m = 100 reads the same 4,221 at 0.878.
- *The cliff.* The rule's reading ends sharply, at the same level on every seed. At k = 4, m = 200 ARI climbs from 0.15–0.18 at T = 0.0907 to 0.53–0.55 at T = 0.0373, then falls to 0.35–0.37 at T = 0.0348. Since all three seeds fall together, this isn't a seed's flip. It looks like the cut growing past k, though I didn't read the cut there.
- *Nesting.* At m = 1 no pairs break it at θ up to 0.65 on any seed. Above that noise near a high cut does: up to 5,919 pairs at θ = 0.7, 8,207 at 0.95 and 377,609 at 0.99, worst at k = 1. At m > 1 it breaks more than closure: 55,699–8,928,522 pairs at θ = ½, against closure's 14,999–1,586,144.
- *Blocks at k = 4, m = 200, T = 0.0373, seed 20261102.* 1,610 (class 3 at 0.74), 1,466 (7 at 0.82), 1,115 and 428 (1 at 0.92 and 1.00), 390 (5 at 0.90), 370 (2 at 0.96), 250 (5 at 0.46), 244 (4 at 0.61). Two-thirds of class 4, 418 of its 626 points, go into the blocks of 3 and 7. At m = 100 and T = 0.0399 it gets a block of 154 at purity 0.91, holding 140 of its points. Outside clusters of at least 50 sit 562 points in 370 small clusters.
- *Time.* A tree takes a median of 0.45 s on K = 5, and 8 s and 10 s on K = 10 and 20, up to 20 s, with flows stopped at seven edges.

**After the fact.** None of this was pre-registered. I ran it because a win this large needed checking before anyone leaned on it.

*The numbers aren't the code's.* 426 kept partitions, recomputed in separate Python, match the rows in ARI, purity and points with no difference at all.

*The gain isn't SPC's.* The best readings at k = 4 sit where 92–95% of the edges are above θ, and the edges crossing between big clusters reach a co-membership of 1.000 at T = 0.0373–0.0525. So at those levels the sampler holds at least one such pair of sides together in nearly every draw, and can't see that split. To test what the sampling adds, the rule ran on the graph alone, with no sampling: each of 65 levels keeps the strongest share of edges by coupling, from 1/65 to all (`MeasureThinCutOnTheGraph`).
- On `k05` it reaches the search at 25 of 30 settings, within the reading's margins, even against the best seed's ARI and points.
- At or above SPC's best seed: ARI at 16 of 30 settings, points at 21. SPC's largest lead is 0.069 in ARI (k = 6, m = 20), and the graph's is 0.156 (k = 2, m = 400).
- It has no seeds, so it has no flip.
- It falls off a cliff at about the same share of edges. At k = 4, m = 200 that's between 95.4% and 96.9%, where SPC's falls between 95.3% and 96.0%.
- `cauchy` is no independent check here, since both kernels rank the edges alike: 63 of 2,015 rows differ.
- Where SPC's column does lead is plain closure on points: 2,392–2,460 at purity 0.85 against 1,618. So the sampling adds something a threshold on coupling doesn't, but the rule doesn't draw on it.

*The scale.* K-means on the raw features, from a stdlib script over three seeds, reads ARI 0.528–0.530 with 6 clusters at purity 0.74 over all 6,435 points, and 0.490–0.548 with 8 at 0.76–0.80. So ARI 0.55 is an ordinary level on this data. Closure's and the search's were low mostly because they stop hot and leave points out. HDBSCAN on the raw features (`ThinCutComparators`, with noise as singletons) reads ARI at most 0.332. With leaf selection and a minimum size of 50 it reads 2,700 points at purity 0.963, against 1,057 at most from anything above at 0.95. So the rule moves the middle of the frontier and not its top.

**What this says (proposal, the user's call).**
- *The search isn't worth building on this evidence.* By the reading fixed in advance, a rule reading only edges goes past it on K = 5. On K = 10 and 20 the search showed no gain to begin with.
- *The chain on K = 5 is a property of the graph.* The cut tree finds it without sampling. So this doesn't show a better way to read SPC. It shows that on this graph a large part of what closure loses, a graph partitioner keeps.
- *What's open.* Whether a fixed setting holds on other data, which on the graph alone needs no sampling at all. How to choose T, k and m without labels. Whether the T-stack, SPC's real reading, gains anything from it, which isn't measured.
- *T and k as two axes (an option, not a finding).* A larger k reads colder, and at m = 1 the classes coarsen as edges are added and refine as k rises. So (T, k) is a grid over nesting classes, in the spirit of the multiparameter frame, which isn't settled. The size rule that makes the numbers work is what breaks the nesting.

**Limits.** Labels choose T on every side, and the setting on the rule's. Landsat is one data set, and the settings came from cuts read with labels. Three seeds on `k05` and one on every other arm. The control's share of edges and SPC's T are different axes, each with its best level chosen by labels. HDBSCAN and k-means read the raw features and not this graph, and I tried few of their settings. The flat cut isn't the T-stack. The tree can miss a thin cut at k of at least 2, as said above, and I didn't measure how often on these columns.

## Proposed order

Not decided.

1. **The closure gap's measurement, probe, θ test and thin-cut test.** Read, above. On this evidence the search isn't worth building, and it stays with the tests. What's open is whether a thin-cut rule at a fixed setting holds on other data, which on the graph alone needs no sampling.
2. **Idea 3.** Stepping-stone and path sampling, validated on exact ln Z (a tree for any q, Kaufman at q = 2), then BARS on ⟨E⟩ checked against the same truth. Independent of 1.
3. **The partition distance table,** then VI cuts, credible balls and P(Q) as its reductions, and the per-block window through BARS.
4. **Idea 5.** BARS with `StepBasis` on lineage trajectories.

A placement question for step 3: partition comparisons are observables (`LinkOverlap` in `Graphs.Observables`), but `ContingencyTable` sits in `Clustering.Evaluation`, above the layer that would need it for VI.

## What bears on stopping, for the user to forward

- The per-block window is a measurement-side readout under §13.3's decision 4. Its one-sided events are monotone, so they would be endpoint certificates in the family of §13.5 for a block known in advance. Blocks aren't known at stopping time, so this is a note, not a proposal.
- The closure gap, over edges, could have flagged intervals where the θ-partition misleads even though C38's sandwich certifies it unchanged. Measured 2026-11-02, it wouldn't: on the real sets it stays at or under 2.6% at every level, and reads 0.0% at Landsat's best level, where the pair form reads 67.7%. The pair form sees more but costs N²/2 per kept partition, needs kept partitions, and doesn't separate chaining from a melt.
- Closure's reading can flip across seeds at a budget that reads steady elsewhere. On Landsat at K = 5, with 8 replicas × 2,000 cycles, one merge of about 900 points decides T = 0.0791–0.0847, and the ARI there runs 0.09 to 0.25 by seed, against a range of at most 0.013 at other levels. The original seed breaks nesting there by 41% of the joined pairs. Exact closure can't do that (C17), so it is noise at the cut, in a band a stopping rule would want flagged. Read since: the merge goes through a cut of one or two edges, none from repair, with δ̄ on the widest path at 0.50–0.56, so the flip may be one edge's noise on this sparse graph. Refusing cuts of up to two edges between sides of at least 100 points takes the range to 0.023 at θ = ½, which supports that. How it relates to C38's sandwich isn't examined.
- Budgets per temperature, set from measured ESS and R̂, were raised here (2026-10-30) as a way to spend draws. They clash with BARS as it stands: its run API takes no per-point weights (primer, "Equal weights"), so uneven budgets give it uneven noise it can't be told about. The same objection §13.3 raises against uneven spacing applies. They stay a diagnostic until BARS takes weights.

## Naming

Binder's loss (D. A. Binder 1978) must not read as `BinderCumulant` (K. Binder 1981); logged in `issues/notation-ledger.md`. Proposed, for the user's call: **pair loss** over all pairs, **link loss** over graph edges, the second matching `LinkOverlap`.

## Citations

Bissiri et al., Binder 1978, Wade & Ghahramani and Dahl et al. are from memory, as `opus-gvdv-ideas.md` says. So are Gomory & Hu 1961, Gusfield 1990, Matula 1972 and Hartuv & Shamir 2000. Grimmett's Theorem 3.8 is checked against the local copy.
