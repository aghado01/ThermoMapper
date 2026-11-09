# The IRLS solver's scale constants (chip brief)

Date: 2026-11-08. Status: commissioned by the user ("chip off for that fix"), from the Grassmann discussion's probe. The user may add other ambient geometry fixes in this session; each gets its own commit and its own line in the report. Reports go under "Report" below, newest first.

## Why

`Irls.Solve` (`src/maths/geometry/solver/irls.cs`) is not scale-equivariant. `IrlsOptions.Default` (`irls-options.cs`) sets four constants in the data's absolute units:
- `SubgradientThreshold = 1e-4`: within this distance of a data point, the median switches to subgradient steps;
- `Eta0 = 1`: those steps are η₀/√(k+1) times an unnormalized sum, absolute. [RiemannianGeometry.md](../../docs/Maths/RiemannianGeometry.md#the-irls-engine-irlssolve) already records that "the step is not normalized by Σw or by the data scale";
- `Tolerance = 1e-8`, absolute by default (`ConvergenceCriterion.Absolute`);
- `Epsilon = 1e-10`, the coincidence radius.

**The probe** ([the discussion's Report](../regression/discussion/grassmann-stiefel-subspaces-brief.md#report), "probe" entry) gives the evidence. Its code and output are in the shared notes folder's `grassmann-probe/`. It used clusters of 21 points with spread s on Gr(2, 10) and Gr(3, 20), plus a Euclidean control, with truth from outside the solver, warm-started at the medoid as `MoMPCA` does:
- **With the defaults, below s = 1e-4,** the median runs to `MaxIterations` (200) and returns an estimate about 5% of s from the median. This happens on every manifold, the Euclidean one included. At s ≥ 1e-3 it converges in 6 to 9 iterations, to about 1e-6 of s.
- **With `HybridMode.WeiszfeldOnly`,** the absolute tolerance stops the iteration early as s nears it: 6e-6 of s at 1e-4, 3e-4 at 1e-6, 4e-2 at 1e-8.
- **The Karcher mean (L²)** is unaffected down to the Grassmann arccos floor, about 1e-7. That floor is a separate item, below.

**The callers:**
- `MoMPCA`'s four medians. Its Euclidean factor is in the data's units, so data in small units (SI lengths at millimetre scale, say) is exposed.
- `DistributedSpred`.
- `EuclideanMedian`/`EuclideanMean` (`estimators/ambient/`) have their own loop with absolute `tolerance = 1e-8` and `epsilon = 1e-10` defaults. `GaussianMixtureModel` calls `EuclideanMedian.ComputeWithScatter`. They were not probed. Check them, and fix them in the same pass if they share the defect.

## What done looks like

- **Equivariance is the contract.** On a flat manifold, the median and mean of c + s·u equal c + s·(the estimate for u), to a relative accuracy that does not depend on s, across at least 1e-12 ≤ s ≤ 1e6.
  - On a curved manifold the same holds down to the manifold's own resolution. For `GrassmannManifold` today that is about 1e-7, the arccos floor; the sphere's is to be measured.
  - Iteration counts stay bounded, and none ends at `MaxIterations` on these cases.
- **The decision is the chip's, with its alternatives on the page.** How the scale is defined and when it is computed, for example:
  - the warm start's median distance to the data, computed once;
  - recomputed per iteration;
  - the data's median pairwise distance.

  The scale should be robust to outliers, since the median exists to resist them. The threshold, the subgradient step and the tolerance are then read relative to it, and so is the coincidence radius if it should scale.
- **Nothing is removed without a question.** Today's absolute semantics stay reachable as an option, and the default becomes the scaled one. Removing an option or a criterion goes to the user as its own question, with what is lost. If the change moves low bits of existing outputs, report which tests and oracles moved and why. Do not loosen a tolerance to absorb the move.
- **Docs.** `docs/Maths/RiemannianGeometry.md`'s IRLS section states the new scale, its alternatives and their trade-offs.

## Gates

- **Truth from outside the solver.**
  - The Euclidean median by equivariance: the median of the u_i by a plain Weiszfeld at scale 1.
  - Curved manifolds by the tangent-space construction the probe uses, exact to O(s³).
- **Chosen cases plus a seeded random sweep.** Dimensions, n, spreads from 1e-12 to 1e6, centres far from the origin, weights unequal, a coincident point and a collinear configuration (where the median is a data point and the optimality check must return it).
- **The callers keep their results at ordinary scales.** `MoMPCA`'s tests, `GrassmannMedianParityTests` (oracle), `DistributedSpred`'s tests and GMM's tests pass unchanged, or any move is reported.
- **A scaled `MoMPCA` case.** The same data times 1e-6 gives the mean times 1e-6 and the same subspace.
- Tests run by [targeted runs](../../CONTRIBUTING.md#testing-standards) (`scripts/test-select.ps1`): the selected classes per edit, the commit scope per commit, and the affected assemblies once at the end.

## Not in this brief

- **The Grassmann arccos floor.** It needs a stable principal-angle routine (atan2 of cosines and residual sines) and `LogMap` by Edelman, Arias and Smith's atan form. The Grassmann discussion has proposed that this routine live in LinAlg beside `ThinSvd`, shared by Geometry and Regression. That placement waits for the user. The probe's `SineGrassmann` is a working reference for it, and the user may add it to this session.
- **`StiefelManifold`**, which ignores the frame and has no consumer.
- **E8 and anything under `src/maths/regression`.**

## Working rules

- **Commits.** On `main`, committing by pathspec (`git commit -- <paths>`), with nothing left staged and nothing pushed. Commit messages describe the work, with no attribution or co-author lines ([CLAUDE.md](../../CLAUDE.md)). Leave `issues/regression/bass/basin-hold-validation.csv` and `path-profile-validation.csv` out of every commit.
- **Shared checkout.** Other sessions are writing in it (E8's slices among them). Touch only the solver, the estimators, their tests, the geometry docs and this brief. Before editing a file, check `git diff -- <file>` for changes you did not make. Build from a snapshot (`git archive HEAD` plus your files), with `ArtifactsRoot` set inline into `artifacts\chip-irls\` and deleted at the end.
- **Machine load.** The machine is shared, so run at below-normal priority with capped threads and MSBuild node reuse off.
- **Review.** Each step ends with a fresh-eyes review of the diff, fixed before reporting.
- **Follow-ups.** Each goes in the report with the decision it would change.

## Report

### 2026-11-08: MoMPCA's calibration, scale-equivariant (finding 1, at the user's go-ahead)

**Commit.** `fa8776d`. The entry below this one records finding 1: the brief's 1e-6 `MoMPCA` case failed on `MoMPCA`'s own α̂ clamp, not on the solver. The user said to fix the clamp here.

**What changed.** `MoMPCA` had two constants in absolute units. Both are now relative:
- **The α̂ clamp.** α̂ was clamped to [1e-6, 2 − 1e-6], and the subspace weight was taken as 2 − α̂. The clamp bounds τ̂_U/τ̂_μ, a ratio in data units. Now α̂ = 2τ̂_U/(τ̂_μ + τ̂_U) and 2 − α̂ = 2τ̂_μ/(τ̂_μ + τ̂_U) are each computed directly. Each τ̂ is floored at its factor's resolution, which still keeps both scales positive when a majority of blocks agree exactly:
  - (ε·max‖xᵢ‖)²/p, the rounding of the block means, which rescales with the data;
  - 2ε/[r(p−r)], the arccos resolution of a Grassmann distance, in radians.
- **The Jacobi tolerance** of the projected-covariance eigendecomposition was an absolute 1e-12. At s ≤ 1e-6 it skipped every rotation and returned the in-subspace axes unrotated; the widened test found it. It is now 1e-12 of the matrix's largest entry.

**Gate.** `ComputeMoM_RescaledData_GivesTheRescaledMeanAndTheSameSubspace` runs s = 1e-12, 1e-9, 1e-6, 1e-3, 1e2, 1e3 and 1e6 on two fixtures: the tilted blocks (α̂₁ ≈ 0.15) and the corrupted 2-D blocks (α̂₁ ≈ 0.05). The mean, basis and variances hold within 1e-10 of scale 1, and α̂ within 1e-9 of 2/(1 + s²(2 − α̂₁)/α̂₁). The brief's 1e-6 case passes. With the clamp, the two fixtures moved at 1e-6 and 1e3 by 0.6% to 2.6% in the mean and 0.002 to 0.03 rad in the subspace.

**Moves.** α̂ changes where the clamp bound at scale 1. On the circle and ellipse fixtures the clean frames coincide (τ̂_U = 0) and the block means agree to rounding, so both τ̂ sit at their floors. α̂ goes from the clamp's 1e-6 to 2 within rounding, with a subspace weight near 1e-16. The joint median still follows the frames, since the means differ only by rounding. The calibrated-scale test asserted the clamp's range [1e-6, 2 − 1e-6]; it now asserts [0, 2]. Every other `MoMPCA` test passes unchanged, and the geometry assembly passes whole (68 of 68). No other code calls `MoMPCA`.

### 2026-11-08: the spread unit and the optimality check at the nearest data point (chip-irls)

**Commits.** `61a7fac` reads the solver's constants in units of the data's spread. `e70afaf` tests the optimality condition at the nearest data point and reads a radius cut short at the cut locus as the distance. This report is its own commit. Nothing was removed; the private `CheckOptimality` became `VertexIsOptimal`, which covers what it tested.

**The scale.** `IrlsOptions.Scale` picks the unit λ of `Tolerance`, `Epsilon`, `SubgradientThreshold` and `Eta0`:
- `IrlsScale.Spread`, the default and the enum's zero value. λ = σ, the weighted lower median of the radii from the current iterate to the data, recomputed each iteration. When half the weight or more sits on the iterate, σ is the weighted mean radius; when all of it does, the iterate is returned.
- `IrlsScale.Data`. λ = 1, the old semantics. `DataScale_ReproducesTheSolverInDataUnits_BitForBit` (pinned) matches afc814d's outputs bit for bit on a Euclidean median, an interior Grassmann median, a Grassmann mean, and `EuclideanMedian` with its scatter.
- **Why per iteration.** At convergence σ is the median distance from the median (50% breakdown), whatever the warm start. A σ fixed at the warm start inherits a contaminated start, such as `MoMPCA`'s arithmetic-mean start, and loses its inflation factor in relative accuracy. The median pairwise distance costs O(n²). [RiemannianGeometry.md](../../docs/Maths/RiemannianGeometry.md#the-unit-of-the-constants-irlsoptionsscale) has the table.

**What else `Spread` needed.** The probe found each of these, and each is stated in the docs:
- **Radii from ‖Log_p xᵢ‖.** The arccos `Distance` has relative error about ε/2θ², which stalls the median below s = 1e-4 even in spread units. The n log vectors of an iteration are kept (n·D doubles) and reused by the step.
- **The distance at the cut locus.** The Grassmann Log drops a principal angle at π/2, and the sphere's drops an antipode, so the Log comes out short of the distance, even zero. On curved manifolds the solver also computes `Distance`, and takes it as the radius where the Log is more than 0.1% short of a distance above 1e-6, the arccos noise of zero being about 2e-8. A curved iteration still costs n distances and n Log maps. `DistributedSpredTests.AggregateProjections_CorruptedFirstBlock_StillFindsCleanMajority` found this once the optimality check was in.
- **The shift is the tangent step's norm, plus a stall stop.** The sphere's and Grassmann's Exp drop tangent steps below 1e-12 (absolute). An iterate far from the origin cannot move by less than its ulps. Both kept the step above τσ while the iterate stopped moving. The solver now also stops when a step moves the iterate by at most 4ε‖p‖.
- **Flat steps through the tangent round trip,** p + Σω(xᵢ − p)/Σω, so a far centre costs only the rounding of p. The flat mean takes the same form from the warm start.
- **The subgradient step divided by Σw.**

**The ambient loop.** `EuclideanMedian` shared the defect: its floor ε and its stop τ(1 + ‖p‖) were absolute, so below a spread of about 1e-4 it stopped short of the median, by 1e-3 of the spread at s = 1e-4 and by half the spread at 1e-9. It takes a `scale` parameter (`Spread` by default) with the same σ, a floor εσ, a stop at max(τσ, 4ε‖p‖), and steps of xᵢ − p. With σ = 0 its final weights are the outer weights, not w/0. `EuclideanMean` is closed form with no constants.

**The optimality check (`e70afaf`).** The old check fired only at exact coincidence and used one coincident point's weight. Weiszfeld reaches a median that is a data point only in the limit, and the `Hybrid` switch then throws the iterate off again. From any warm start but that point, such a median ran to 200 iterations and stopped up to 6e-4 of the spread off, at every scale and in either unit. The collinear, duplicated and 5-dimensional cases showed it. The solver now tests the nearest data point once it coincides with the iterate, or once it holds half the IRLS weight together with its duplicates. It tests each point once, sums the weight sitting on it, and returns it exactly. A point on the candidate's cut locus enters the bound at its full weight. This applies in both units.

**Gates** (`IrlsScaleTests`, `MoMPCATests`):
- **Flat equivariance.** 12 seeded draws: d from 2 to 20, n from 5 to 40, equal and unequal weights, an occasional ×50 outlier. Centres at 0 and at 100s, spreads 1e-12 to 1e6, warm starts at the medoid and the mean, truth from a plain Weiszfeld at scale 1. The median holds 1e-6 of the spread (the probe: at most 8e-8), the mean 1e-12 (at most 1e-14), and `EuclideanMedian` 1e-6 (2.2e-7). No s costs more iterations than s = 1.
- **A far centre.** At |c| = 10 the error stays within 1000·ε‖c‖/s, the representation limit. At s = 1e-12 the data themselves hold only 2e-3 of the spread.
- **Curved.** S⁹, Gr(2, 10) and Gr(3, 20) for s from 1e-1 to 1e-6, in 7 to 12 iterations. The median holds about 1e-8 down to s = 1e-4, then about 1e-12/s from the Exp guard. **The sphere's resolution, measured:** the same as the Grassmannian's. At s = 1e-7 the Log maps' arccos drops angles below about 1e-8. S² misses by 7%, and one Gr(2, 10) draw in three runs to `MaxIterations`. The mean holds 1e-8 down to 1e-7 on both.
- **Chosen cases.** Each median is returned exactly:
  - collinear points with unequal weights, at every s, from the mean and from another data point;
  - a duplicated point holding half the weight, at every s, from the mean and from the medoid;
  - identical points, with a zero scatter;
  - frames on each other's cut locus, from the medoid.
- **Scaled `MoMPCA`.** At s = 1e-3 and 1e2 a balanced fixture (α̂ ≈ 0.15) reproduces scale 1 to 1e-14 in the mean, basis and variances. afc814d missed the 1e-3 case by 0.012 in the basis. **The brief's s = 1e-6 is not tested.** `MoMPCA` clamps α̂ to [1e-6, 2 − 1e-6]. The clamp binds once τ̂_U/τ̂_μ leaves about [5e-7, 2e6], and there the mean moves by about 1% and the subspace by about 0.01 rad. That non-equivariance is `MoMPCA`'s, not the solver's (finding 1).

**Callers.** `MoMPCATests` (6 old and 2 new), `DistributedSpred`/`Spred` (39) and `GrassmannMedianParityTests` (R ran) all pass unchanged. The selected classes ran at each commit, and the three assemblies ran whole at the end: Geometry 63 of 63, DimReduction 83 of 83, Oracle 209 of 211. The two failures are `DiagramMetricsTransportParityTests` (p = 1, 2), where the R process `transport_oracle.R` times out at 120 s before any comparison, alone as well as in the full run. That test uses only `TDA.Ph`, and the change does not reach it (finding 6). No test reaches GMM's `EuclideanMedian` call: only `MapperGMM` calls `RobustInitialize`, and it has no test. Low bits move because the stop is now τσ. For example, `MoMPCA`'s mean on the balanced fixture moves 1e-9 relative at scale 1.
- **Moved and called out:** the parity test's seven frames have a data point for their median, frame 0, where ‖G‖ ≈ 0.90 ≤ 1. The solver now returns frame 0 exactly. Before, it stopped about 2e-4 from it, and the oracle's 2e-3 tolerance hid that. The pinned test carries this case with the new expectation, and pins an interior six-frame case bit for bit against afc814d.

**Findings, each with the decision it would change:**
1. **`MoMPCA`'s α̂ clamp** makes it scale-dependent beyond about 1e3 either way from τ̂_μ ≈ τ̂_U, so the brief's 1e-6 case fails. The fix would compute the two factor weights 2τ̂_U/(τ̂_μ + τ̂_U) and 2τ̂_μ/(τ̂_μ + τ̂_U) directly, avoiding the cancellation in 2 − α̂, and guard only τ̂ = 0. That is a `MoMPCA` change outside this brief, and it decides whether the 1e-6 gate can pass.
2. **Weiszfeld's own rate.** Near a data point that barely fails the optimality condition (‖G‖/w → 1⁺), Weiszfeld contracts at w/‖G‖. One sweep draw (d = 3, n = 5, ‖G‖/w = 1.07) needs about 205 iterations at τ = 1e-8, at every s and in either unit. It ends at `MaxIterations`, still within 1e-6 of the spread. This decides between accelerating (a line search along the Weiszfeld step, or a Newton step after a rejected check with ‖G‖/W near 1) and a larger `MaxIterations`.
3. **`EuclideanMedian` stalls when warm-started on a data point.** That point's weight w/(εσ) holds the iterate, so it returns the start, in both units. GMM warm-starts at the cluster mean, so nothing reaches it today. This decides whether `EuclideanMedian` gets the optimality check (delegating to `Irls.Solve`), which would move GMM's robust initialization.
4. **The manifolds' absolute guards.** The sphere's and Grassmann's Exp drop tangent steps below 1e-12, which bounds the median's relative accuracy at about 1e-12/s. Their Log maps drop angles below about 1e-8. The cut-locus fallback depends on `Distance` being right near π/2. All of this belongs with the arccos-floor item and its LinAlg placement, which waits for the user.
5. **GMM's robust initialization has no test.**
6. **`DiagramMetricsTransportParityTests` times out in R** (`transport_oracle.R`, 120 s) on this machine, unrelated to this change. This decides between a longer oracle timeout and a look at the R environment.

The probe, its rows and the bit-for-bit pin probe are in the shared notes folder `local_58551eda-672f-47a3-8439-28c0db7e8f81/irls-scale-probe/`.
