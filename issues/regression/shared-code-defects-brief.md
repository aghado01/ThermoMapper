# Two defects in shared numerics: the weight reader and the legacy Cholesky factor (chip brief)

Commissioned by the user on 2026-11-13 through the coordinating session. Both defects were found by conformal slices and parked for their owners ([handoff](coordination-handoff-2026-11-11.md#open-defects-in-shared-code)). Each returns wrong numbers with no error, to every caller that reads through it.

**Read first:** [AGENTS.md](../../AGENTS.md), [CONTRIBUTING.md](../../CONTRIBUTING.md), then the two records linked below.

## Slice 1: the weight reader at extreme weights

**The site.** `PredictiveChecks.Rows` (`src/maths/regression/PredictiveChecks.cs`) takes `log(weight / total)`. `Pointwise.CheckWeights` (`src/maths/regression/Pointwise.cs`) admits the weights it then misreads.
- `[1e300, 1e-300]`: the second ratio underflows and its log weight reads −∞. Subtracting logs keeps about −1381.55.
- `[1e308, 1e308]`: the sum overflows, `CheckWeights` accepts it, and both log weights read −∞.

**The record.** Conformal slice 10's report, ["Parked shared weighted-row follow-up" and "Follow-up verification criterion"](conformal-layer-brief.md#2026-11-11--slice-10-density-and-cdf-error-diagnostics). The criterion there is this slice's acceptance:
- the first input keeps its finite second log weight; the second returns two weights of 1/2; no NaN and no lost row;
- the normalized log weights, and the weighted predictive, CDF and carried-error readings built on them, agree with an independent high-precision or log-sum-exp reference;
- a finite log-likelihood tilt that makes the first input's tiny row material;
- invariance to a positive rescaling of the weights, and zero-weight controls;
- the shared consumers' selection and the whole assembly pass.

**Its readers.** The predictive checks (two sites in `PredictiveChecks.cs`) and nine conformal sites: `ConformalBayes` (three), `SmoothedPValues` (two), `CrossConformal`, `SplitConformal`, `SplitPitCalibration`, `Tournament`. Search for any reader that normalizes weights on its own; the repair is one shared routine that all of them read.

**The constraint.** Weights in the ordinary range keep their bits, so every pinned hash and fingerprint holds. Scaling by a power of two is exact, which makes this reachable. If no repair keeps the ordinary range bit for bit, stop and report the pins that would move: regenerating a pin is the user's call.

**Gates.** The criterion above; chosen cases and a seeded random sweep over the weights' whole dynamic range, subnormals and values near `double.MaxValue` included; the pinned tier at each commit; the whole regression assembly once at the end.

## Slice 2: `CholeskyDecomposition` on a near-singular covariance

**The site.** `src/maths/linalg/Cholesky.cs`. A pivot at or below 10⁻¹² is set to 10⁻¹², whatever the matrix's scale.

**The record.** Conformal slice 11's report, ["Defects found" and "The shared factor's floor"](conformal-layer-brief.md#2026-11-11--slice-11-bgps-drawn-f-the-retained-parameters-content-the-full-step-takes).
- On a smooth kernel's covariance, whose spectrum falls below rounding, the rows under a floored pivot divide their rounding error by 10⁻⁶. Later pivots went to −3.7·10¹⁵⁷ and were floored again; 55 entries of the factor were infinite.
- `LogDet` stayed finite (−683.6). Only `Status.Finite` reported the failure.
- BGP left for `DenseCholesky` (`5db2e40`), and its drawn function carries its own diagonally pivoted factor (`BgpDrawnFunction`).

**Remaining callers.** GMM: `GaussianComponent` and `BhattacharyyaCoefficient` (two sites), in `src/clustering/statistical/gmm/`. Viz-core: `src/viz-core/serializer.cs`. Tests and cost profiles under `tests/linalg/` and `tests/regression/` also construct it.

**Reading first.** Reproduce the failure in a test on that class of covariance. Then measure, on a ladder of condition numbers and scales, where the floored factor's solves and `LogDet` leave a reference's tolerance while every entry is still finite. The overflow is the visible end of the defect; wrong finite numbers come before it. Read how `DenseCholesky` treats the same inputs (`Status`, `FailedAt`, a NaN `LogDet`, `RelativePivot`), and what each remaining caller does with `Status` today.

**Build.** A factor with a non-finite entry is a failed factorization, and the class refuses it the way `DenseCholesky` does: no finite `LogDet`, and no solve, sample, inverse or factor read from it. Each remaining caller reads the outcome and handles it. The two factors then fail one way.

**For the user, reported with the measurements and a recommendation, and left unbuilt.** Each changes GMM's numbers or removes a surface:
- *The floor.* An absolute 10⁻¹² against a floor relative to the matrix's scale, and whether a pivot negative beyond rounding is floored or reported.
- *A pivoted factor in the linear algebra.* `BgpDrawnFunction`'s diagonal pivoting (Harbrecht, Peters and Schneider 2012) as a shared primitive that near-singular callers take.
- *Retiring `CholeskyDecomposition` into `DenseCholesky`.* What is lost (the soft floor, `WriteInverseTo`'s L⁻¹, `Sample`), and which callers, tests and cost profiles move. A removal is its own question.

**Gates.** The reproduction fails before the build and passes after. Every existing test of the class holds bit for bit on inputs the floor never reached. The pinned tier at each commit; the linear algebra, clustering and viz-core test projects that reach the class, and the regression assembly, once at the end.

## Working rules

- **Order.** Slice 1, then slice 2. They touch different code, so slice 2's reading starts when slice 1 lands. The coordinating session verifies each landed slice from a clean archive ([scripts/verify-csharp.ps1](../../scripts/verify-csharp.ps1)).
- **Files.** This chip owns `PredictiveChecks.cs`, the weight check in `Pointwise.cs`, `Cholesky.cs`, the three caller files, and their tests. Three other chips are at work:
  - the conformal chip (`local_69bae1e5`), on the remaining slices 14 to 18, owns `src/maths/regression/conformal/**` and reads `Rows`; a change inside that folder is agreed with it by message first;
  - the GBASS chip owns BASS's sampler files for the u–s funnel fix;
  - two Lean chips work in `lean/`.
- **Git.** On `main`, one commit per concern with its tests, committed by pathspec, nothing left staged, nothing pushed, no attribution lines.
- **Builds and runs.** `ArtifactsRoot` set inline to `artifacts\chip-defects\`, deleted at the end. Heavy runs at below-normal priority with capped threads; the machine is shared. Tests by [targeted runs](../../CONTRIBUTING.md#testing-standards) (`scripts/test-select.ps1`).
- **Review.** Each slice ends with a fresh-eyes pass for correctness and coverage gaps, fixed before the report.
- **Scope.** New beside old. Findings are framed by problem class. A follow-up is reported with the decision it would change and what each outcome would do. `docs/` stays untouched; the report lists "Doc changes owed". Sub-agents need the user's permission.
- **Notes.** An "In progress" entry under Reports holds uncommitted work, running checks and the next step. Probe code and data go to `artifacts/`.
- **Events.** End a turn with one line and the report's link when a slice lands or needs a decision. Results go in the report.

## Reports

Newest first. Each gives the commits, the gates' results, defects found, decisions for the user and "Doc changes owed".
