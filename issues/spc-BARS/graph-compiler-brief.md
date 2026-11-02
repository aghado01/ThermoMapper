# Graph compiler revision: brief

**Status: built 2026-11-01** (§2.1–2.3, §2.7, §2.9 and all of §4.1; §2.4–2.6 and §2.8 deferred, see the Report). Raised 2026-10-25 by the stopping-theory exploration (`stopping-theory-primer.md` §12.27–12.29,
claims C30–C34). The user has more compiler changes in mind than these. They go in §4 before this brief is handed to a
session.

## 1. Why

The graph compiler (`src/graphs/GraphCompiler.cs`) decides what structure SPC can see at all. Its stages are topology,
filter, repair, refinement, then projection and scaling.

On `BlattHierarchy`, a Gaussian graph shows only the finest level. A Cauchy graph with the same topology shows all three
levels, with stage widths a rigid-block model predicts to within a grid step. A union of KNN scales puts every level
into one T-stack.

Along the way the compiler turned out to cut edges it had just added. It also ties topology to coupling scale, and it
lacks couplings and topologies the probe had to build by hand. This brief collects those findings and what to change.

## 2. Findings and what to change

### 2.1 The coupling floor cuts repair edges (C34): the fix

- **Mechanism.**
  - Stage 3, repair, adds the minimal MST bridges: `pipeline/repair/MstMinRepair.cs` calls
    `ConnectivityRepair.EnsureConnected`.
  - Stage 5 then scales every edge. The scaler keeps an edge only if `coupling > 1e-10`, at
    `pipeline/scalers/GlobalBandwidthScaler.cs:174` (single kernel) and `:282` (mixture kernel).
  - A Gaussian bridge across a wide gap has coupling e^(−d²/2a²), where a is the bandwidth: about 10⁻³¹ on
    BlattHierarchy's coarser gaps. So the bridge is cut after the repair stage has counted it.
- **Reporting.** Coverage is checked after repair (`GraphCompiler.cs:138`) but not after scaling. The compiler reports a
  connected graph and returns a disconnected one.
- **The collapse detectors can't see it either.** They count edges under 1e-8 (`GraphCompiler.ComputeNearZeroRatio`,
  `AutoLmpScaler.ComputeNearZeroRatio`), and the edges under 1e-10 are already gone.
- **Evidence.** Mutual KNN, `MstMin`, `Refinement.Auto`, mean-edge-distance bandwidth:

  | Dataset | Gaussian | Cauchy, same topology |
  |---|---|---|
  | `BlattHierarchy.Generate(pointsPerFine: 50, leafSpread: 0.15, seed: 42)`, K = 5 | 0 edges between medium clusters, 0 between the coarse pair | 4 and 1 (the MST bridges) |
  | Same, K = 10, 20, 40 | likewise 0 and 0 | 4 and 1 |
  | `AnisotropicGaussian.Generate(pointsPerCluster: 400, meanSeparation: 4, anisotropyRatio: 1, seed: 55)`, K = 5 | S₂ = 0.333 at T = 10⁻⁶: three disconnected classes | S₂ = 1 |

- **Consequence.** A stage that begins at the grid's cold end on such a graph measures the grid, not a coupling, so
  stage widths stop being comparable across kernels. The shared-neighbour and mutual-reachability graphs the probe built
  on the compiled edge set inherit the cut (§2.6).
- **Options** (the user chooses):
  - **a. No floor:** keep every edge the topology stages produce, at its kernel coupling.
    - Connectivity and the cuts' statistics stay honest.
    - A 10⁻³¹ bridge is effectively broken at any reachable temperature, but it's still in the graph, and S₂ at T → 0
      reflects it.
    - SW pays for edges that never activate.
  - **b. A floor for repair edges only:** keep repair edges at a documented minimum coupling, and cut other edges below
    the threshold as now.
    - Repair then means something under any kernel.
    - But the floor is a coupling the kernel didn't produce, and the floor's value sets where the coarser stages land.
  - **c. Keep the cut, and report it:** after scaling, log the component count, and warn when repair edges are cut.
    - Honest, but it leaves the disconnection in place.
  - **Any option needs connectivity checked and logged after scaling** (§2.2).

### 2.2 Report connectivity of the final graph

- Log the component count (and largest-component coverage) of the graph the compiler returns, not only of the
  topology.
- Make the near-zero detectors count edges before the cut, or count the cut edges directly.

### 2.3 K and the coupling scale move together (C31)

- **The coupling.** The bandwidth is the mean edge length over the KNN edges, so raising K both adds edges and widens
  the kernel. On BlattHierarchy the mean edge length goes from 0.07 at K = 5 to 2.1 at K = 320.
- **The topology.** K selects by point count: edges cross a gap only once K exceeds the cluster's size.
- **The result:** one K shows one hierarchy level, because the finer gaps fill with edges as the bandwidth grows.
- **Change.** Make topology and coupling scale separate choices. For example, take the bandwidth from a fixed
  neighbour rank regardless of K, or set it explicitly, with the current coupling kept as a selectable option.

### 2.4 Multi-scale topology (C32)

- **The construction.** Take the union of the edge sets of several KNN scales. Couple each scale by its own kernel at
  its own bandwidth a_s and its own 1/K̂_s, weight it by w_s, and sum per edge.
- **The probe version** (`MultiScale` in the scratch probe) builds each scale with the Cauchy kernel, so that no repair
  edge is cut, then re-couples them with Gaussians.
- **Weights.**
  - Equal weights flood: the largest scale's edges hold everything at the coarse level.
  - w_s = (a_min/a_s)² makes the sum a discrete Gaussian scale mixture whose envelope falls as d⁻². It showed coarse
    6.5, medium 2.0 and fine 2.45 e-folds on BlattHierarchy.
- **What the weights decide.** They place each scale's stages along T. The widths within a scale come from the data.
  So the weight rule is a modelling choice and should be documented as one.
- **Related:** `coupling/Mixture.cs` mixes bandwidths on one edge set. The union differs because each scale brings its
  own edges.
- **Cost:** 155,412 edges at N = 1,200 for scales 5, 40 and 320. Sweeps are about 70× slower than K = 5.

### 2.5 An exponent for the power-law kernel

- `coupling/Cauchy.cs` fixes J = 1/(1 + d²/δ²), whose tail falls as d⁻².
- A general exponent α, either J = 1/(1 + (d/δ)^α) or a Student-t family, makes the tail d^(−α).
- **Prediction.** On a single bridge, stage width ≈ α·ln(the ratio of gap lengths). Cauchy (α = 2) matched: ln 70 =
  4.25 e-folds for the medium stage, where SW gave 4.26.
- Larger α separates the plateaus more cleanly, against the fixed width of each split (about 0.75 e-folds at q = 20).

### 2.6 Shared-neighbour and mutual-reachability couplings

- **Neither is in the compiler.** The probe re-weights the compiled edge set, so these graphs inherit §2.1's cut.
- **Add them as couplings over the topology's own edge set:**
  - shared neighbours: J = max(|N(i) ∩ N(j)|/K, floor)/K̂, with the neighbour count and the floor as parameters
    (Jaccard is a variant);
  - mutual reachability: a kernel on max(core distance of i, core distance of j, d).
- **Reuse:** core distances and mutual-reachability MSTs already exist for HDBSCAN (`primitives/mst/Prim.cs`).
- **On iris** at K = 5, shared neighbours were the one lever that widened the three-species stage (§12.25).

### 2.7 LMP at small K fragments every dataset

- At K = 5, LMP rescaling leaves S₂ at T = 10⁻⁶ at 0.09 (iris), 0.02 (AnisotropicGaussian), 0.012 (BlattHierarchy)
  and 0.010 (Landsat). At K = 10 it behaves.
- This looks like a problem with its local scale at small K, not with the data.
- **Check:** how many edges §2.1's cut removes under LMP at K = 5, and what `LocalMutualProximityScaler` does with a
  5th-neighbour scale.

### 2.8 `MstAll` (deferred)

- `GraphCompilerConfig.cs:106` defers adding the whole minimum spanning tree, on the grounds that there's no use case.
- **A possible use case:** from memory, BWD superimposed the full MST on the mutual-KNN graph. This is unverified, so
  check the paper.
- It would put a single-linkage path under every level of a hierarchy, not only across the gaps that need repair.

### 2.9 How a gap is coupled, as an explicit choice

- **The closest pair:** a repair edge kept at its kernel coupling. This is single linkage. Under a power-law kernel,
  the stage widths are the single-linkage merge heights on a log-temperature axis. It's exact on clean gaps and chains
  through noise.
- **A surface:** many crossing edges, as in the union or at a larger K. Closer to how a magnet couples, and robust to
  stray points, but the edge counts enter the widths.
- **Today the choice is an accident** of K and the 1e-10 cut. Make it a stated option.

## 3. Constraints

- **Pre-release:** replace wholesale, with no shims, `[Obsolete]` or migration windows. Write the new code beside the
  old, and delete what it supersedes in a final cleanup pass.
- **Defaults stay per engine, with the current behaviour kept as a selectable option** beside any new default.
- **The BlattHierarchy results are a probe, not a calibration.** No default or API shape is set from one dataset.
- **Document decided choices,** with their alternatives and trade-offs, in `docs/Graphs/GraphConstruction.md` (or a
  maths doc where the choice is statistical).
- **Searching:** use `rg` from the shell; the Grep tool is hooked off in this repo.
- **Commits:** targeted commits on main, each message ending with
  `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`. Worktrees fail in this repo (long paths).

## 4. Further changes (the user's, to add)

_The user has more compiler changes planned. List them here, with any ordering against §2, before this brief goes to a
session._

### 4.1 Order of operations: where the gauge and the scale are decided

Raised by the user on 2026-10-30: the gauge is decided at the wrong point relative to the other stages, and the split
with the cumulative-field solver (formerly PKWang inside SPC) is part of the same problem. Traced in code at `f52697e`.

**The order today, on the SW path.**
1. Topology, then filter, which symmetrizes the graph (OrRule or MutualKnn).
2. Repair adds the MST bridges.
3. Refinement.
4. Scaling, in `pipeline/scalers/GlobalBandwidthScaler.cs`:
   - estimate the bandwidth from the post-repair neighbour lists;
   - apply the kernel;
   - multiply by 1/K̂, only under `BandwidthStrategy.MeanEdgeDistance`;
   - drop every edge under the 1e-10 floor.
5. The caller then picks the temperature bracket.

**Findings.**
1. **The gauge is a side effect of the bandwidth choice.** `couplingScale = n / totalNeighbors` exists only in the
   `MeanEdgeDistance` branch of `ScaleSingle` (`GlobalBandwidthScaler.cs:147–153`). The mixture path (`ScaleMixture`)
   and both LMP scalers (`AutoLmpScaler`, `LocalMutualProximityScaler`) never apply it. Whether a graph has a gauge is
   decided by how its bandwidth is estimated.
2. **The gauge comes before the floor, so it changes the edge set.** `coupling *= couplingScale; if (coupling > 1e-10)`
   (`:172–174`): a larger K̂ (a larger K) shrinks every coupling and pushes more edges under the cut. The gauge is
   equivalent to rescaling T, as the parked coupling-normalization arc treats it, only when it comes after every
   threshold.
3. **K̂ and the bandwidth are measured on the wrong edge set.** Both read `refined.AllNeighbors`: post-repair and
   pre-cut.
   - The repair bridges, the longest edges by construction, enter the mean-edge-distance bandwidth, so a topology fix
     rescales every edge's coupling.
   - They also count in K̂.
   - K̂ is not the mean degree of the graph the compiler returns, since the cut comes after it.
4. **The pre-repair comparison graph changes two things.** `scaler.Scale(filtered, n)` (`GraphCompiler.cs:284`), built
   for `MstBridge.Compare`, re-estimates its own bandwidth and K̂ on the filtered edges. Comparing it with the returned
   graph mixes a topology change with a change of units.
5. **Temperature units are chosen downstream, by two rules that disagree.**
   - The REPL's `auto` bracket reads the bandwidth option: T_ps(q) under `MeanEdgeDistance`, the data heuristic
     otherwise (`user-repl/commands/SpcCommand.cs:262–280`).
   - The library's `AutoGridFixedSweep` always uses the heuristic (`AutoGridFixedSweep.cs:45`). The REPL's comment says
     the heuristic lands about K̂× too cold under the gauge; that factor is not measured here.
   - So one graph config gets different brackets from the REPL and the library. The parked arc's seam, "the coupling
     construction declares its T-estimate", is implemented as a caller reading a compiler option, not as the
     construction declaring it.
6. **The cumulative-field split.**
   - `CumulativeFieldSolver` sets its units after compile, from its own cumulative ladder: a cut sits at
     T = Hcum/ln 2, and `EstimateBracket` reads that ladder.
   - Its `Field.Local` (a per-node ladder) is a density-adaptive normalization applied after compile, the same family
     as LMP, which acts during compile. Nothing stops the two stacking, for example an LMP-scaled graph fed to
     `Field.Local`. The parked arc flagged this double counting.
   - Before the move (`e059df0`), `SpcCommand` called grid resolution "the LAST step shared by both siblings", yet
     `auto` gave each its own bracket.
   - SW's T and the cumulative-field T are in different units, and no artifact records either.

**A proposed order, with each stage owning one concern** (gauge versus measure as the parked coupling-normalization
note distinguishes them; for discussion, nothing decided):
1. **Structure:** topology, filter and symmetrization, repair. It ends with the edge set E.
2. **Scale:** estimated from a declared edge set, by default the data's own neighbourhoods before repair, and recorded.
3. **Shape:** the kernel gives raw couplings on all of E, with no floor.
4. **Measure, optional, one place:** LMP, a degree normalization, the α family, or a per-node field. These change the
   structure, so exactly one is chosen and stated; the cumulative-field `Field.Local` is a member of this stage, not a
   second, post-compile one.
5. **Thinning, optional and explicit:** edges dropped for cost are dropped on raw couplings, relative to the kernel's own
   scale, with repair edges protected or reported, and never after the gauge.
6. **Gauge, last:** one declared global factor, for example 1/K̂ with K̂ the mean degree of the final graph, or none.
   That makes it exactly a T rescaling.
   - The manifest records it, with the temperature estimate it implies.
   - SW, the cumulative-field solver and every bracket read their units from the graph, not from a CLI option, so the
     REPL and the library agree.

**Against §2.**
- §2.1's floor question is partly settled by placement: a floor in stage 5, before the gauge, with repair edges
  protected, is close to option b.
- §2.2's connectivity log belongs after stage 5 and again after stage 6, which should change nothing if the gauge is
  pure.
- §2.3 is stage 2's choice of edge set and estimator; separating K from the bandwidth falls out of it.
- §2.4's multi-scale union needs its per-scale 1/K̂_s to become stage 6's single gauge, or be stated as part of stage 4.
- The current behaviour stays selectable beside any new default (§3).

### 4.2 Decisions, 2026-11-01

The session took the decisions to the user, who answered "Build as revised". Measured first, on `7a35ffd` (`GraphCompilerReturnedGraphTests` and the
REPL): the BWD toy, Iris, and Landsat parity configs (mutual kNN, `MstMin`, Gaussian, mean-edge) cut **no** edge, so a floor change leaves Domany parity
and `BwdParityTests` untouched. On BlattHierarchy (Gaussian, mean-edge, K = 5, 10, 20) the floor cuts 5 repair edges and the graph goes from 1 component to
6; under MAD at K = 5 it cuts 25 repair edges and goes to 27. A Cauchy kernel on the same topology keeps one component.

1. **Floor and LMP.**
   - Repair edges are never cut. They keep their kernel coupling, clamped to 1e-300 where the kernel underflows. The clamp is applied before the gauge, so
     the gauge stays a pure rescale and 1e-300/K̂ stays a normal number. (`double.Epsilon` would round to 0 under 1/K̂.)
   - Every other edge keeps the 1e-10 thinning, now a recorded config field `CouplingFloor` (0 = none).
   - A kernel-honest caveat: this restores the graph's structure (components, what PH, Mapper and the cumulative-field ladder see). It does not change what SW
     sees on a Gaussian graph: a bridge at 1e-31 activates near T = 1e-31, so S₂ at T = 10⁻⁶ still reads the separate components. No floor value fixes that
     short of a large modelling constant, and none is added.
   - Under LMP a repair bridge is the weakest edge at both ends, so the strict rank writes it exactly 0 (`GraphCompilerReturnedGraphTests`). Decided:
     (a) repair edges join the LMP protect set, keeping their kernel coupling, as the H1 cycle edges do; (b) an inclusive rank rule (the empirical CDF, which
     counts the edge itself) is added as an option beside the strict one, default strict, so no existing LMP graph moves from it. The strict rank also writes 0 on
     every node's weakest edge, so a degree-1 node loses its only edge: that is the diagnosis of §2.7's fragmentation at small K.
2. **Order.** §4.1's stages, with thinning before the measure step: structure, scale, shape (kernel), thinning (repair edges protected), measure (LMP), gauge
   last. The floor is defined on kernel couplings, and LMP outputs rank products. The gauge is its own field (`Gauge`), K̂ is the returned graph's mean degree, and
   `MeanEdgeDistance` stops implying it.
3. **Edge set.** The mean-edge bandwidth reads the pre-repair edges by default (the bridges are the longest edges by construction); the post-repair reading stays
   selectable (`ScaleEdges`), and the BWD presets and the parity config select it, since the papers' mean includes the superimposed MST. The other strategies
   already read the pre-repair k-th sample. `PreRepairGraph` becomes the returned graph minus the repair edges: same units, topology change only, and the second
   scaling pass goes.
4. **The auto-LMP trigger is held.** Commit `7a35ffd` made the trigger count edges the floor dropped. LMP fragments sparse graphs (above), so a graph that loses
   5% of its edges to the floor must not switch to it. The trigger returns to counting the kept edges; the `MaxNearZeroEdgeRatio` interrupt and the log keep
   counting the cut edges. Revisit when the rank rule is settled.
5. **Temperature units.** The build result declares its resolved construction (kernel, bandwidth estimator, gauge factor, resolved LMP). One bracket function reads
   it, and the REPL's `auto` and `AutoGridFixedSweep` both call it. The q-only T_ps(q) bracket applies exactly when the declaration is the classic one (Gaussian,
   mean-pair bandwidth, 1/K̂ gauge, no LMP); otherwise the data heuristic divided by the gauge (the heuristic is ∝ J², the gauge ∝ J). **This moves `auto` grids
   for Cauchy or LMP graphs with mean-edge bandwidth**, which get T_ps today. Mapper patches carry no declaration and keep the heuristic. The T_ps correction
   itself (C44) is `681bf5c`, a separate commit by another session. The cumulative-field command warns when LMP meets `Field.Local`, and both T-unit systems are recorded.
6. **Resolved config and knobs.** `GraphConstructionManifest` gains a structured `Resolved` record, beside the requested config: what each `Auto` picked, the
   statistic and threshold behind it, the bandwidths, the gauge, the floor, the cut count, and the returned graph's components. The four buried Auto thresholds
   (skewness 3.0, coverage 0.95, near-zero 1e-8, near-zero ratio 0.05) become one `AutoThresholds` record with those defaults. The knob register is a table in
   `docs/Graphs/GraphConstruction.md` (concern, default, entanglements, invariance test), and the invariance tests are written with it.
7. **`ScaleRank`** (§2.3). Optional r ≤ K. For the mean-edge bandwidth the sample is the unique undirected pairs among each node's first r directed neighbours,
   before repair; for the MAD family, the r-th neighbour distance. Null keeps today's behaviour.
8. **Sequencing.** Now: 1–7 above. §2.9 folds into 1 (the closest pair is the default; the surface is K). Deferred, each landing on the restructured stages:
   §2.4 multi-scale union, §2.5 α-kernel, §2.6 shared-neighbour and mutual-reachability couplings, §2.8 `MstAll` (check the BWD paper first). The user has not
   yet said whether §4 holds more changes.

## 5. Done when

- A test reproduces §2.1's cut: a repair edge whose Gaussian coupling falls under the floor, and the compiled graph's
  component count. The user's chosen option makes it pass.
- The final graph's connectivity is logged (§2.2).
- Each other item the user keeps from §2 and §4 is implemented with tests, or explicitly deferred with a reason.
- The graphs and SPC tests are green, and the decided choices are documented.
- The work is reported below, including whether stage widths on Gaussian graphs change for far-apart clusters.

## Report

### 2026-11-01: the compiler revision

**What landed** (all on `main`, one concern per commit):

| Commit | Concern |
|---|---|
| `7a35ffd`, `b0b3186` | §2.2: the returned graph's component count and largest-component coverage, over all edges and over edges of positive coupling; the floor's dropped edges; what became of the repair bridges. The near-zero interrupt counts the dropped edges. |
| `63589b3` | The auto-LMP trigger reads the kept, non-repair edges only (`7a35ffd` had let dropped edges switch it on). |
| `be04f35` | §2.1: `CouplingProjection.CouplingFloor` (default 1e-10; 0 = none); repair edges are exempt, kept at their kernel coupling or 1e-300 on underflow. The §2.1 test reads connected. |
| `7178dae` | LMP keeps the repair edges at their kernel coupling; `LmpRankRule.Inclusive` beside `Strict` (§2.7). |
| `7506bf6` | §4.1 finding 4: `PreRepairGraph` is the returned graph without the repair edges; the second scaling pass is gone. |
| `55f97f9` | §4.1 findings 1 and 2: `Gauge` is its own field, applied last; the mean-edge bandwidth no longer implies it. |
| `5f60984` | §4.1 finding 3 and §2.3: `ScaleEdges` (pre-repair by default) and `ScaleRank`. |
| `fe32dd7` | `Resolved` in the result and the manifest; `GraphCompilerConfig.Auto` (`AutoThresholds`). |
| `f7da1f2` | §4.1 findings 5 and 6: the construction declares its units; `SpcScheduleHelpers.AutoBracket` is the one bracket function; the cumulative-field warning and the basis in both run manifests. |
| `5b6d357`, `9741b4c` | `docs/Graphs/GraphConstruction.md` (order, decisions with alternatives, knob register), the notation ledger, REPL wiring tests. |

**Checked against §5.** The §2.1 test (`GraphCompilerReturnedGraphTests.GaussianRepairBridge_UnderTheFloor_StaysInTheReturnedGraph`) passes under the chosen option. The final graph's connectivity is logged, over all edges and over coupled edges. Tests: `VizCore.Tests` 467 passed, 2 skipped (394 before the work; 73 new in `GraphCompilerReturnedGraphTests`, `GraphCompilerGaugeTests`, `GraphCompilerScaleTests`, `GraphCompilerResolvedTests`, `GraphUnitsBracketTests`, `GraphCompilerConfigWiringTests`, `LmpRankRuleTests`), `Clustering.Graphical.SPC.Tests` 13, `Clustering.Graphical.CumulativeField.Tests` 15, `TDA.Mapper.Tests` 28, `Graphs.Primitives.Tests` 3, all passing; `BwdParityTests` is unchanged and green. `TDA.DimReduction.Tests` (83) passed at `b0b3186` and again at `5b6d357` (1 h 25 m under load).

**Deferred, with reasons.** §2.4 (multi-scale union), §2.5 (an exponent for the power-law kernel), §2.6 (shared-neighbour and mutual-reachability couplings) and §2.8 (`MstAll`) each add a topology or coupling that lands on the restructured stages, and nothing else waited on them; §2.8 needs a read of the BWD paper first. §2.9 is the floor decision: the closest pair is the default, the surface is K.

**What moved.**
- *Graphs gaining edges:* those whose repair edges the floor used to cut. BlattHierarchy, Gaussian, mean-edge, K = 5, 10, 20 goes from 6 components to 1; MAD at K = 5 keeps 25 repair edges it lost. The BWD toy, Iris and Landsat parity configs cut no edge, so they are unchanged; `BwdParityTests` passes with its expected values untouched.
- *Changed defaults:* the mean-edge bandwidth reads the pre-repair edges (the BWD presets and the parity config select `AfterRepair`); the 1/K̂ factor is no longer implied by the mean-edge bandwidth (the presets and the parity config declare `Gauge`); the gauge's K̂ is the returned graph's mean degree; LMP graphs with repair edges keep them (they were written 0).
- *`auto` temperature grids* move for Cauchy or LMP graphs with the mean-edge bandwidth, which took T_ps before; the heuristic is now divided by the gauge. Graphs without the gauge are unchanged.
- *Graph caches:* every cached graph rebuilds (the config fingerprint covers the new projection fields). Older run manifests that used `mean-edge-distance` re-materialize without the gauge in `extract`.
- `docs/Clustering/SpcEngine.md` (line 163) still says the REPL's `auto` reads `MeanEdgeDistance`; it carries another session's uncommitted edits, so it was left.

**Stage widths on Gaussian graphs for far-apart clusters.** Before, the cut left no connection: the 1,200-point BlattHierarchy was six components, and there were no merge heights above the finest level. Now the bridges are in the graph, so the heights exist, but they are not where SW can reach them.
- *Gaussian, mean-edge, K = 5* (single-linkage merge heights, log₁₀ J): the finest merge at −13.5, then five merges at −300. The medium and coarse gaps underflow, so their repair edges sit on the 1e-300 clamp and tie: the medium and coarse stages are zero e-folds apart and some 640 e-folds below the fine one. At K = 10 the finest merge is at −7.3 and the rest tie at −300.
- *Cauchy, same topology:* −1.8, −3.7, −3.7, −3.8, −4.5: data-determined, adjacent merge heights 0.2 to 4.4 e-folds apart.
- *AnisotropicGaussian* (meanSeparation 4, K = 5), Gaussian: merges at −13.9, −10.7, −7.1, −4.4, so the widths are data-determined here, and the weakest bridge is 1.2·10⁻¹⁴. Its bond probability at T = 10⁻⁶ is about 10⁻⁸, so S₂ at T = 10⁻⁶ still reads three classes, as before: the fix restores the structure, and not the dynamics at a reachable temperature.
- A stage across a far Gaussian gap therefore lies at T ≲ 10⁻¹⁴ and, once the kernel underflows, at the clamp. Storing couplings in log space would give the underflowed heights their data again; that is not done.

**Findings for the user.**
1. *LMP fragmentation, §2.7.* The strict rank writes 0 on each node's weakest edge. Mutual kNN + MST-min, Cauchy, K = 5, no injected protect set: Iris 111 of 262 edges are zero and the coupled graph has 30 components; Landsat 3,754 of 9,048 and 1,862; BlattHierarchy 828 of 2,169 and 213. The inclusive rank has no zero edge and one coupled component in all of them, and K = 10 behaves the same way (Iris 6 components strict, 1 inclusive). The default stays strict; flipping it moves every LMP graph, and is the user's call.
2. *The heuristic against the corrected T_ps(20) = 0.0892.* The data heuristic on the gauged parity graphs lands 0.15 to 0.41 times it; divided by the gauge it lands 1.17 to 1.6 times it (toy 1.6, Iris 1.4, Landsat 1.5, Cauchy toy 1.17). The "K̂× too cold" in the REPL's old comment was measured against the old, wrong T_ps. Dividing by the gauge is exact (`g·h(raw)` and `h(gauged)/g` agree to four digits); the residual 1.2 to 1.6 is the heuristic's own calibration.
3. *The floor still cuts topology edges.* Only repair edges are exempt. In Blatt, MAD, K = 5, mutual kNN + MST-min, one non-repair edge dropped under the floor splits off a 100-node region (the log says so). Protecting graph bridges (cut edges), not only repair edges, would close it.
4. *Untested:* the spherical local-parametrix path for a repair edge beyond 0.9π (it keeps the 1e-300 clamp by code, with no test); the `Mixture` kernel is tested for the bridge only.
5. *Not done, §4:* the user had not said whether §4 holds further changes.


