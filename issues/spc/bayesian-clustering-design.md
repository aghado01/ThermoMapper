# Bayesian clustering ideas for SPC: design note

Status: started 2026-10-30. The closure gap is built (2026-11-01) and not yet measured on data; the rest is discussion. Decisions are the user's and are marked as such; the rest is proposal.

Sources: `opus-gvdv-ideas.md`, Threads 1, 5 and 7 of `opus-spin-glass-threads.md`, `sw-retention-brief.md`, and the stopping primer's §13.3 and §13.6 (`issues/spc-BARS/stopping-theory-primer.md`).

Symbols are the module's: τ_ij = P(i↔j), the pair connectivity; p_e = 1 − e^(−J_e/T), the bond probability; δ̄ = ((q−1)τ + 1)/q, the spin agreement that `DenseTStack` cuts.

## The user's direction

2026-10-30, stated as a leaning and not as settled: pursue what fits a uniform temperature grid read by BARS. The primer's §13.3 records the grid decision; its decision 4 keeps BARS out of the stopping rule, so anything read through BARS here is on the measurement side.

2026-11-01, the user's decision: build the closure gap first, as idea 1's first step (see "Built"). No other idea is chosen yet.

## Status

| # | Idea | Status | Fit with the uniform grid and BARS |
|---|---|---|---|
| 1 | Loss cuts, against closure | Buildable. Link loss over graph edges from co-membership; pair loss over all pairs from kept partitions | A column kernel for `DenseTStack`; the closure gap below is a curve over T |
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

**Not done.** It hasn't been measured on the Domany battery or ISOLET. That is an experiment, for the user to ask for. Filling `DenseCoMembership` costs N²/2 increments per kept partition under a giant cluster, about 6 × 10⁹ per temperature at N = 7,800 and 200 kept partitions.

## Proposed order

Not decided.

1. **Measure the closure gap.** The edge form on the existing co-membership, and the pair form on kept partitions, scored as a curve over T on the Domany battery and ISOLET. If it is small everywhere, a loss-optimal search isn't worth building. If it is large, it says where.
2. **Idea 3.** Stepping-stone and path sampling, validated on exact ln Z (a tree for any q, Kaufman at q = 2), then BARS on ⟨E⟩ checked against the same truth. Independent of 1.
3. **The partition distance table,** then VI cuts, credible balls and P(Q) as its reductions, and the per-block window through BARS.
4. **Idea 5.** BARS with `StepBasis` on lineage trajectories.

A placement question for step 3: partition comparisons are observables (`LinkOverlap` in `Graphs.Observables`), but `ContingencyTable` sits in `Clustering.Evaluation`, above the layer that would need it for VI.

## What bears on stopping, for the user to forward

- The per-block window is a measurement-side readout under §13.3's decision 4. Its one-sided events are monotone, so they would be endpoint certificates in the family of §13.5 for a block known in advance. Blocks aren't known at stopping time, so this is a note, not a proposal.
- The closure gap, over edges, could flag intervals where the θ-partition misleads even though C38's sandwich certifies it unchanged. Unmeasured.
- Budgets per temperature, set from measured ESS and R̂, were raised here (2026-10-30) as a way to spend draws. They clash with BARS as it stands: its run API takes no per-point weights (primer, "Equal weights"), so uneven budgets give it uneven noise it can't be told about. The same objection §13.3 raises against uneven spacing applies. They stay a diagnostic until BARS takes weights.

## Naming

Binder's loss (D. A. Binder 1978) must not read as `BinderCumulant` (K. Binder 1981); logged in `issues/notation-ledger.md`. Proposed, for the user's call: **pair loss** over all pairs, **link loss** over graph edges, the second matching `LinkOverlap`.

## Citations

Bissiri et al., Binder 1978, Wade & Ghahramani and Dahl et al. are from memory, as `opus-gvdv-ideas.md` says. Grimmett's Theorem 3.8 is checked against the local copy.
