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
