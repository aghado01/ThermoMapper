# Graph compiler revision: brief

**Status: open.** Raised 2026-10-25 by the stopping-theory exploration (`stopping-theory-primer.md` §12.27–12.29,
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

## 5. Done when

- A test reproduces §2.1's cut: a repair edge whose Gaussian coupling falls under the floor, and the compiled graph's
  component count. The user's chosen option makes it pass.
- The final graph's connectivity is logged (§2.2).
- Each other item the user keeps from §2 and §4 is implemented with tests, or explicitly deferred with a reason.
- The graphs and SPC tests are green, and the decided choices are documented.
- The work is reported below, including whether stage widths on Gaussian graphs change for far-apart clusters.

## Report

_Appended by the session that does the work._
