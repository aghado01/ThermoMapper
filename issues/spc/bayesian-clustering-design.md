# Bayesian clustering ideas for SPC: design note

Status: started 2026-10-30. The closure gap is built (2026-11-01) and measured on data (2026-11-02; see "What the measurement found"); the rest is discussion. Decisions are the user's and are marked as such; the rest is proposal.

Sources: `opus-gvdv-ideas.md`, Threads 1, 5 and 7 of `opus-spin-glass-threads.md`, `sw-retention-brief.md`, and the stopping primer's §13.3 and §13.6 (`issues/spc-BARS/stopping-theory-primer.md`).

Symbols are the module's: τ_ij = P(i↔j), the pair connectivity; p_e = 1 − e^(−J_e/T), the bond probability; δ̄ = ((q−1)τ + 1)/q, the spin agreement that `DenseTStack` cuts.

## The user's direction

2026-10-30, stated as a leaning and not as settled: pursue what fits a uniform temperature grid read by BARS. The primer's §13.3 records the grid decision; its decision 4 keeps BARS out of the stopping rule, so anything read through BARS here is on the measurement side.

2026-11-01, the user's decision: build the closure gap first, as idea 1's first step (see "Built"). No other idea is chosen yet.

## Status

| # | Idea | Status | Fit with the uniform grid and BARS |
|---|---|---|---|
| 1 | Loss cuts, against closure | The gap is built and measured. The search is buildable (link loss over graph edges from co-membership; pair loss over all pairs from kept partitions), and the measurement argues against building it now | A column kernel for `DenseTStack`; the closure gap below is a curve over T |
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

**What this says (proposal, the user's call).**
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

## Proposed order

Not decided.

1. **The closure gap's measurement.** Read, above. The probe of Landsat's case is what could change the proposal.
2. **Idea 3.** Stepping-stone and path sampling, validated on exact ln Z (a tree for any q, Kaufman at q = 2), then BARS on ⟨E⟩ checked against the same truth. Independent of 1.
3. **The partition distance table,** then VI cuts, credible balls and P(Q) as its reductions, and the per-block window through BARS.
4. **Idea 5.** BARS with `StepBasis` on lineage trajectories.

A placement question for step 3: partition comparisons are observables (`LinkOverlap` in `Graphs.Observables`), but `ContingencyTable` sits in `Clustering.Evaluation`, above the layer that would need it for VI.

## What bears on stopping, for the user to forward

- The per-block window is a measurement-side readout under §13.3's decision 4. Its one-sided events are monotone, so they would be endpoint certificates in the family of §13.5 for a block known in advance. Blocks aren't known at stopping time, so this is a note, not a proposal.
- The closure gap, over edges, could have flagged intervals where the θ-partition misleads even though C38's sandwich certifies it unchanged. Measured 2026-11-02, it wouldn't: on the real sets it stays at or under 2.6% at every level, and reads 0.0% at Landsat's best level, where the pair form reads 67.7%. The pair form sees more but costs N²/2 per kept partition, needs kept partitions, and doesn't separate chaining from a melt.
- Budgets per temperature, set from measured ESS and R̂, were raised here (2026-10-30) as a way to spend draws. They clash with BARS as it stands: its run API takes no per-point weights (primer, "Equal weights"), so uneven budgets give it uneven noise it can't be told about. The same objection §13.3 raises against uneven spacing applies. They stay a diagnostic until BARS takes weights.

## Naming

Binder's loss (D. A. Binder 1978) must not read as `BinderCumulant` (K. Binder 1981); logged in `issues/notation-ledger.md`. Proposed, for the user's call: **pair loss** over all pairs, **link loss** over graph edges, the second matching `LinkOverlap`.

## Citations

Bissiri et al., Binder 1978, Wade & Ghahramani and Dahl et al. are from memory, as `opus-gvdv-ideas.md` says. Grimmett's Theorem 3.8 is checked against the local copy.
