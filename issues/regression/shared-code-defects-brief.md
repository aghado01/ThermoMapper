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

### In progress (2026-11-13)

- *Slice 5, `BandCholesky` strict by default:* decided by the user on 2026-11-13. `PivotedCholesky`'s solve is still undecided and untouched.
- *The user's correction to slice 4's note on rollout.*
  - BARS's banded g-prior already asks for strict mode.
  - The BAPS penalized and tensor splines use the floored default but reject any factor that reports a floored pivot.
  - So the switch can admit a valid small-scale fit that BAPS rejects today. That direction is qualified explicitly, beside singular-design refusal and the public λx = λy = 1e20 tensor cases.
- *Next:* read `BandCholesky`'s floor, strict mode and every caller; then the switch, the three qualifications, the gates and the report.

### 2026-11-13 — Slice 4 verified from clean archives (coordinating session)

Over the nineteen files of `f448931`, `43bd61f`, `ff5d635`, `e6dc220` and the review fix `78603a2`.

- **The regression project:** `scripts/verify-csharp.ps1` on a clean archive of HEAD `a3ed9c0`, the Commit selection: 308 of 314 classes, 3,049 executed, 3,049 passed, 0 not executed.
- **The other projects,** from a clean archive of HEAD `8db288f`, after a sweep of every test project for the classes the files select (`closure.txt`):
  - `Maths.LinAlg.Tests`, whole: 190 of 190.
  - `TDA.Mapper.Tests`, whole: 28 of 28.
  - `VizCore.Tests`, whole: 703 passed and 2 skipped of 705. The skips are the two Möbius golden drafts.
  - `Maths.Samplers.Tests`, by selection: 8 of 8.
  - `Maths.Oracle.Tests`, by selection with the pinned tier: 146 of 147.
- **The one failure** is `DiagramMetricsTransportParityTests.DiagramDistances_MatchT4TransportAndLpAssign(p: 1)`: `transport_oracle.R` timed out after 120 s, inside the R call. The test names none of the changed types, and its p = 2 case passed. It is the report's finding for the TDA owner, and no part of this slice.
- **The other nine projects** that compile a changed file have no selected class and did not run.
- **The source was read** at `e6dc220`, and the two notes it drew landed in `78603a2`: the scene export's ε relative at every scale, and the remainder pass's cost in `PivotedCholesky`'s doc. Nothing else was found.
- **Receipts:** `artifacts/verify-defects-s4-a3ed9c0/` and `artifacts/verify-defects-s4-projects-8db288f/`, tracked in [receipts](receipts/) with the project script as it ran.

### 2026-11-13 — Slice 4: slice 2's three decisions, built

The user chose on 2026-11-13: the review's recommendations (`discussion/astra-shared-code-cleanup.md`), with this chip's two adjustments. The scene export draws a singular block instead of refusing it, and the pivoted factor's remainder check reads the off-diagonals too. The review's order: the floor first, the pivoted factor apart, the callers migrated and qualified, the legacy class deleted in its own commit.

**Commits.**
- `f448931` fix(linalg): CholeskyDecomposition drops its pivot floor.
- `43bd61f` feat(linalg): PivotedCholesky, the semidefinite factor BGP's drawn function reads.
- `ff5d635` refactor(gmm, viz): a Gaussian's covariance factor is DenseCholesky.
- `e6dc220` chore(linalg): retire CholeskyDecomposition and BGP's own pivoted routine.
- `78603a2` fix(viz): a singular block's ε is relative below unit scale too (after the coordinating session's review).

**1. The floor dropped (`f448931`).**
- *The class is strict.* The first pivot that is not finite and positive stops it, at the row `DenseCholesky` stops at, and an L⁻¹ past `double.MaxValue` fails it. On every input of a ladder the two agree on the failure row, `Status`, factor and `LogDet`. The ladder: 60 SPD matrices at scales 2⁻⁸⁰ to 2³⁹, 16 squared-exponential kernels, 10⁻¹³·I, and indefinite and non-finite inputs. 10⁻¹³·I now factors exactly, where the floor read every variance ten times too large.
- *The earlier class, kept verbatim and pinned*, reads bit for bit on 60 inputs the floor never reached.
- *The callers.*
  - GMM's `UpdateCache` and the Bhattacharyya distance refuse a covariance that is not positive definite. A warm start with a singular covariance is refused, an acceptance case.
  - The scene export draws a singular block as a thin ellipsoid. It adds ε = 10⁻⁸ of its largest variance, and refuses only a block that is still indefinite. The viewer loads all nine entries as a general 3×3 map (`viewer.html`, `buildGaussianGroup`), so a rank-deficient factor would make its normal matrix singular, and the ε keeps it invertible.
  - Rank one at 10²⁰⁰, which overflowed the floored factor, now exports a finite needle. `[[1, 2, 0], [2, 1, 0], [0, 0, 1]]` is refused, naming its layer and component.

**2. The pivoted factor promoted (`43bd61f`).**
- *The contract is narrow and sampling-first:* the rank, the loadings in the matrix's own row order, and `Map` from Rank standard normals to a draw of N(0, LLᵀ).
- *BGP's loop, verbatim.* The loadings match the earlier routine bit for bit (pinned, 40 kernels and low-rank Gram matrices). BGP's drawn function keeps its seeded draws: its tests, including the bit-for-bit candidate, and the conformal tests that read it pass.
- *The remainder is checked.* The earlier routine read `[[1, 2], [2, 1]]` as rank 1, truncating −3, and `[[0, 1], [1, 0]]` as rank 0. A semidefinite remainder whose diagonals are at most the stop has every entry within it, so the check reads the off-diagonals as well as the diagonal.
- *Measured: two tolerances.* BGP's joint posterior covariances leave remainders up to 1.2·10⁻¹¹ of their largest variance. That is rounding of the computed covariance, past the 10⁻¹² stop: the first form, which judged the remainder at the stop, refused two of BGP's tests.
  - The rank still stops at 10⁻¹².
  - Indefiniteness is judged at its own tolerance, √ε (1.5·10⁻⁸) by default, far below the refused cases.
  - A remainder `[[0, δ], [δ, 0]]` at δ = 5·10⁻¹¹ is kept at the default and refused at 10⁻¹².
- *The cost.* The factor costs O(p·k²) at rank k, and the remainder's check O((p − k)²·k), which is above the factor when p is far above k. For BGP both cost less than forming the covariance. The class's doc says so (`78603a2`).
- *The solve waits*, as the review asked. On a singular matrix it must refuse or project a right-hand side outside L's range, and BGP's sampling needs neither.

**3. The callers moved to `DenseCholesky` (`ff5d635`).**
- *`TryWriteInverseTo`* forms A⁻¹ on demand. It ports the retired L⁻¹ arithmetic, so GMM's Σ⁻¹ keeps its bits (pinned against the arithmetic kept verbatim). It refuses an inverse that overflows: `[1e-320]` factors, and its reciprocal does not fit, and a power-of-two bidiagonal L⁻¹ overflows at row 51. GMM's component is then left unusable (NaN caches, an exception) rather than holding a partial precision matrix.
- *`AddSample`* accumulates L·z into the mean in the retired order. GMM's seeded Box–Muller draws match bit for bit (pinned).
- *The exporter* reads `Factor(i, j)`: the same factor, the same bits.
- *Bhattacharyya* reads its mean term as ‖L⁻¹δ‖² through `InverseQuadraticForm`, so its distances move in their last bits. They agree with the explicit-inverse formula within 10⁻¹³ relative.

**4. The legacy class retired (`e6dc220`).**
- *Deleted:* `CholeskyDecomposition` with its own tests, and BGP's own pivoted routine.
- *Moved to `DenseCholesky`:* `FactorReferenceTests`, `DenseCholeskyTests` and `BandCholeskyTests`.
- *The cost profiles* time `DenseCholesky` as production. Their two rows that timed the retired path's `Decompose` are gone; the replica kernels of its two passes stay.

**The gates.**
- *Each commit*, from a snapshot of HEAD plus its files, by selection with the pinned tier:
  - A: linear algebra 90 of 90, VizCore 9 of 9, regression 80 of 80;
  - B: linear algebra 14 of 14, regression 145 of 145;
  - C: linear algebra 93 of 93, VizCore 12 of 12.
- *The end*, from a snapshot with D:
  - `Maths.LinAlg.Tests` 190 of 190;
  - `TDA.Mapper.Tests` 28 of 28;
  - `VizCore.Tests` 703 passed, with the 2 Möbius skips that stood before;
  - `Maths.Regression.Tests` 3,370 of 3,370 (18 min);
  - `Maths.Oracle.Tests` by selection, 145 of 147. The 2 failures are `transport_oracle.R` timing out at 120 s, alone as well as under load. Their test reads `TDA.Ph` and R only, and no TDA code reaches a changed factor, so the timeout stands apart from this slice (below).
- Receipts are in `artifacts/shared-code-defects/receipts/` (`slice4-*`).

**Defects found.**
- *One in the slice's code, found in the coordinating session's review and fixed in `78603a2`.* The export's ε started its largest variance at 1, so below unit scale it was absolute: a singular block with variances near 10⁻⁸ drew as thick as it was long.
  - It is now 10⁻⁸ of the block's largest variance at every scale, which is what the display needs: a flat ellipsoid stays flat.
  - An all-zero block keeps the missing dimensions' absolute 10⁻⁸, a sphere of radius 10⁻⁴.
  - Gates: VizCore 7 of 7 and linear algebra 14 of 14 by selection with the pinned tier.
- My own first remainder bound refused BGP's covariances, as above. One test used an assertion precision the analyzer refuses, and another rested on noise that could move the rank; both were rewritten before any run counted.
- *Finding: the transport oracle no longer fits its timeout.* `DiagramMetricsTransportParityTests` (both p) times out in `transport_oracle.R` at 120 s, even run alone. It has nothing to do with this slice, and is for the TDA owner: either the timeout or the fixture's size (s = 100 slices).
- *Finding: `BandCholesky`'s default floor is the same rule.* It "mirrored" the retired one: an absolute 10⁻¹², on by default, with strict mode a choice. By the same arithmetic it can introduce order-one error at small scale; this slice did not measure it.

**Decisions for the user.**
1. *`BandCholesky`'s default.* If it becomes strict, the banded marginals refuse a near-interpolating design that the floor now carries, and their numbers move only where the floor fires. If it stays floored, it is the one floored factor left, and its docs should say so. Recommendation: strict by default, as here, with the over-knotted design refused by its caller's prior.
2. *`PivotedCholesky`'s solve*, when a caller needs one: refuse a right-hand side outside the retained range, or project it to the minimum-norm least-squares answer.

**Doc changes owed.**
- *`docs/Maths/LinearAlgebra.md`:*
  - "Dense Cholesky (`CholeskyDecomposition`)" (line 193) and its rows in the API list (line 95), the module table (line 12), the diagram (line 33), the allocation note (line 73) and the storage table (line 81): retired.
  - `DenseCholesky` gains `TryWriteInverseTo` and `AddSample`, and serves the Gaussians.
  - A section for `PivotedCholesky`: its contract, the two tolerances, the remainder rule, its cost (the check's O((p − k)²·k) above the factor's O(p·k²) when p is far above k), and the measurements above.
  - `BandCholesky`'s floor no longer mirrors anything.
- *`docs/Clustering/GaussianMixtures.md`:*
  - line 11: `DenseCholesky`;
  - line 45: `UpdateCache` refuses a covariance not positive definite, or one whose inverse overflows;
  - the warm start (line 89): a singular covariance is refused;
  - "Bhattacharyya distance" (line 162): refusal, and the mean term through the factor.
- *`docs/Maths/Regression/BGP.md`*, the drawn function: it reads `PivotedCholesky`, whose remainder is checked.
- *Slice 2's "Doc changes owed"* for the floor are superseded: there is no floor left to describe.

### 2026-11-13 — Slice 3 verified from a clean archive (coordinating session)

- **Run:** `scripts/verify-csharp.ps1` on a clean archive of HEAD `badc0d9`, the Commit selection over slice 3's fourteen files (`c9f9a84`). All fourteen lie in the regression project's compile closure.
- **Selection:** 237 of 313 classes, with the pinned tier and without the slow tier's 71 methods.
- **Result:** 2,659 executed, 2,659 passed, 0 failed, 0 not executed (29 min). `DrawWeightRangeTests` is among them.
- **The source change was read:** nothing found. Each reader forms a ratio, so the power-of-two scaling cancels and ordinary weights keep their bits. The positivity checks still read the raw weights.
- **The two shared files:** the commit's diff of `SubspacePosterior.cs` and `CoActivePosterior.cs` holds this slice's lines only, and E8 slice 7e's orientation stands in the archive.
- **Not rerun here:** the `Maths.Oracle.Tests` selection, whose R oracles time out under this machine's load. It stands on the chip's run (112 of 116, the four timeouts passing when rerun singly).
- **Receipts:** `artifacts/verify-defects-s3-f1d70d5/` (named for the HEAD when the run was queued; its snapshot line records `badc0d9`), tracked in [receipts](receipts/verify-defects-s3-f1d70d5/).

### 2026-11-13 — Slice 2 verified from clean archives (coordinating session)

- **The standard run:** `scripts/verify-csharp.ps1` on a clean archive of HEAD `78d5077`, the Commit selection over slice 2's six files (`37eb63c`). It selects within the regression test project only: 2 of 308 classes, 79 executed, 79 passed, 0 not executed.
- **The other projects.** Five of the six files lie outside that project, so each test project that reaches them ran whole, from a clean archive of HEAD `ed7df5d`:
  - `Maths.LinAlg.Tests`: 179 of 179, `CholeskyDivergenceTests` among them.
  - `TDA.Mapper.Tests`: 28 of 28.
  - `VizCore.Tests`: 699 passed and 2 skipped of 701, `GaussianFactorRefusalTests` among them. The skips are the two Möbius golden drafts, skipped before this slice.
- **The sweep behind that choice** (`closure-78d5077.txt`): `Cholesky.cs` compiles into 14 test projects.
  - The name graph selects classes in three: linear algebra (7 of 19), VizCore (2 of 85) and regression (2 of 308).
  - Mapper holds three of the files and reaches the component without naming it.
  - `Clustering.Graphical.CumulativeField.Tests` also ran whole in the sweep's first pass: 15 of 15.
  - The other nine have no selected class and did not run.
- **The whole regression assembly** less its slow tier passed with the refusal in place, under slice 1's run (2,922 of 2,922 at `e39924c`).
- **The source change was read:** nothing found. A non-finite entry of L shows in its own row's pivot, and every read is gated on one flag.
- **Receipts:** `artifacts/verify-defects-s2-e429b29/` and `artifacts/verify-defects-s2-projects-ed7df5d/`. The lists, the counters and the one-off project script are tracked in [receipts](receipts/).

### 2026-11-13 — Slice 1 verified from a clean archive (coordinating session)

- **Run:** `scripts/verify-csharp.ps1` on a clean archive of HEAD `e39924c`, the Commit selection over slice 1's three files (`f50753f`).
- **Selection:** 288 of 308 classes of the regression assembly, with the pinned tier and without the slow tier's 71 methods.
- **Result:** 2,922 executed, 2,922 passed, 0 failed, 0 not executed (26 min). `WeightReaderTests` is among them.
- **The source change was read:** ordinary weights take the earlier expression, so their bits hold; nothing was found.
- **Slice 2** was verified later the same day (the entry above). This archive holds its source (`37eb63c`), so the regression assembly passed with the refusal in place.
- **Receipts:** `artifacts/verify-defects-s1-e39924c/`, with the changed-file list, counters and snapshot line tracked in [receipts](receipts/verify-defects-s1-e39924c/).

### 2026-11-13 — Slice 3: BASS's draw weights over their whole range

Added by the user to this chip on 2026-11-13, after the review in `discussion/astra-shared-code-cleanup.md`. It is slice 1's decision for BASS's own readers.

**Commit.** `c9f9a84` fix(regression): BASS's draw-weight readers keep every finite weight's share.

**The readers.** Eleven sites normalized linear weights through a plain sum or summed raw weights times values. All of them read `BassFit.DrawWeights`, an emulator's weights, or a subspace posterior's weights:
- *the selection:* `Draws3.Select`, which prediction, gradients, sensitivity, the active subspaces and the multi-output readouts read;
- *the summaries:* `PosteriorSummary.Over` (the fit's summary, the Sobol summaries) and `SobolDecomposition`'s means;
- *the resamplers:* `NestedCalibrationSampler.OuterDraws` and `BassJointDraws`' resampling;
- *the residual check:* its posterior mean, and `BassDrawFit.MeanOffset`;
- *the projection report:* its mean term count;
- *the subspaces:* `Whitening.Normalized`, `SubspacePosterior.Of` and `CoActivePosterior.Of`.

`DrawSample.Quantile` already reads weights relative to the largest, and `PredictiveChecks.Rows` was repaired in slice 1. The pointwise and candidate densities pass the weights to `Rows`.

**The repair.**
- *One routine, `WeightScale`* (`src/maths/regression/WeightScale.cs`). Each reader takes its weights scaled by 2^−e, the power of two that brings the largest into [1, 2), then does its own arithmetic unchanged.
- *Why ordinary weights keep their bits.* Each reader forms a ratio: a normalized weight, a weighted mean, a cumulative pick against u·Σ. The scaling is exact wherever the scaled weight is normal, and it cancels in the ratio, so ordinary weights give the bits they gave.
- *Why the extremes now read.* A sum of n scaled weights is below 2n, and a weighted sum of values stays below twice the values' own absolute sum.
- *Where the exponent comes from.* The selection takes it over the requested draws, repeats included; the Sobol means over the defined draws; every other reader over all its weights.
- *What stays as it was.* Zero-weight behavior and every refusal message at each entry point. The summary still drops a zero-weight draw from its count and range, testing the raw weight.
- *`PredictiveChecks.Rows`' fallback* now calls the same routine, with the same arithmetic.

**The gates.** Acceptance from the review: both extreme rows reproduce the ordinary row; equal weights resample the same indices at a fixed seed; the selection normalizes over the requested draws, repeats included; zero-weight behavior and the ordinary-range pins hold.
- *The reproduction*, `tests/regression/DrawWeightRangeTests.cs`, with the rows `[1, 1]`, `[1e308, 1e308]` (the sum overflows) and `[5e307, 5e307]` (the sum is finite but the weighted sums overflow). It also uses exact power-of-two multiples (2^±1000) and subnormal weights where their ratios are exact.
  - On HEAD's code, 9 of the 14 tests fail:
    - the selection read `[0, 0]`;
    - the summary's mean read NaN at `[1e308, 1e308]` and ∞ at `[5e307, 5e307]`;
    - every outer calibration draw fell to the last draw;
    - the joint draws resampled other ancestors;
    - the subspace weights normalized to 0;
    - the projection's mean term count read NaN;
    - the residual check threw on a non-finite residual.
  - The 5 that pass on HEAD are the cases that never overflowed, and the pinned comparison.
  - On the build, all 14 pass.
- *What the tests check.*
  - The two summary values read mean 3, SD 1 and median 3 at every row.
  - The selection, the subspace posteriors, the outer draws, the joint draws, and the fit's summary, prediction and sensitivity read the ordinary row's values. So do the residual check and the projection report. Power-of-two multiples match bit for bit; the extreme rows match within 10⁻¹¹ relative.
  - The pinned test compares 2,000 SMC-like weight vectors with the earlier selection and summary, kept verbatim, bit for bit.
- *The whole regression assembly*, from a snapshot of `77f80ba` plus the fourteen files: 3,341 of 3,341 (50 min).
- *`Maths.Oracle.Tests`*, which reaches the readers through 22 classes, by selection with its pinned tier: 112 of 116.
  - The 4 failures are R oracles timing out at 120 s on the loaded machine: `bass_enumeration_oracle.R` three times and `gbass_readouts_oracle.R` once. Each failed inside the R call, before any comparison.
  - Rerun one at a time, all 21 cases of those 4 tests pass (17 min).
- *The merged state.* E8 slice 7e (`199abd7`, the eigenvectors' signs) and conformal slice 15 (`55f4b93`) landed during the run, slice 7e in two of this slice's files. So the commit gate ran again from `048bbc4` plus the fourteen files, by selection (236 of 312 classes) with the pinned tier: 2,646 of 2,646 (12 min).
- Receipts are in `artifacts/shared-code-defects/receipts/` (`slice3-*`).

**Shared files.** E8 slice 7e committed its eigenvector orientation to `SubspacePosterior.cs` and `CoActivePosterior.cs` while this slice's gates ran.
- This slice's three lines in each sit beside that work, and its commit is the working tree, which is that HEAD plus them.
- A first staging, from the snapshot taken before slice 7e landed, would have reverted slice 7e's lines. Its staged diff showed it before any commit, and it was unstaged.

**Defects found.** None in the slice's code. My own test first built `2·1e308/3`, an infinite weight; it was rewritten to divide first before any run counted.

**Decisions for the user.** None for this slice. Slice 2's three decisions stand.

**Doc changes owed.**
- *`docs/Maths/Regression/BASS.md`*, `DrawWeights` (line 699): every view reads any finite weights, over their whole range, by `WeightScale`'s rule. Ordinary weights read as before.
- *`docs/Maths/Regression/Regression.md`*, the weighted draws (line 214): one rule for linear and log weights, `WeightScale` and `Rows`.
- *`docs/Maths/Regression/ActiveSubspaces.md`*, the summaries (line 65) and the posterior (line 105): `SubspacePosterior` and `CoActivePosterior` read their weights over the same range.

### 2026-11-13 — Slice 2: `CholeskyDecomposition` on a near-singular covariance

**Commit.** `37eb63c` fix(linalg): CholeskyDecomposition refuses a factor that is not finite, and its callers read it.

**The reading.** The probe is `artifacts/shared-code-defects/probe/cholesky-floor.cs`, with its tables in `results-before.md` and `results-edges.md` beside it.
- *The reproduction: squared-exponential covariances*, variance 1, on 50 to 401 points of [0, 1], lengths 0.05 to 1, no jitter.
  - 19 of the 20 configurations floor pivots, and 15 of those overflow: 351 to 73,536 entries of the factor are not finite, the first at rows 17 to 57.
  - `LogDet` still reads a finite −1,156 to −11,019 on every one.
  - `DenseCholesky` stops on all 19, at rows 5 to 16.
  - The other 4 floor 40 to 398 pivots and stay finite.
- *Where the floored factor goes wrong while finite.* A = s·Q·diag(κ^(−i/39))·Qᵀ with n = 40, s from 10⁻¹⁶ to 10⁶ and κ from 1 to 10²⁰. The reference is a double-double factor of the stored A.
  - *The floor fires on scale alone.* The pivots lie in [s/κ, s], so the floor fires when s/κ reaches 10⁻¹². At κ = 1, 10⁻¹³·I floors all 40 pivots.
  - *Once a pivot is floored, the outputs leave any tolerance while every entry stays finite.*
    - With 4 to 8 pivots floored, `LogDet` is off by 3.6 to 13 nats and the solves and inverse by 73 to 93 % (s = 1 at κ = 10¹⁴, 10⁻⁶ at 10⁸, 10⁻⁸ at 10⁶).
    - Across every cell where flooring changed a pivot, `LogDet` is off by 3.6 to 1,200 nats and the solves by 46 % up to 2·10¹¹ relative. At s = 10⁻¹², κ = 1, the floor replaces 10⁻¹² by itself, and the error there is 0.
    - The law `Sample` draws from, LLᵀ, is up to 5.9·10⁴ from A in relative terms, among the cells that have a reference.
  - *`DenseCholesky` on the same matrices is accurate to its κ·ε*: 2.4·10⁻¹⁰ nats and 1.4·10⁻¹⁰ at s = 10⁻⁶, κ = 10⁸.
    - It refuses only where the stored A has no factor in double: κ = 10²⁰ at every scale, and κ = 10¹⁸ at most.
    - There the floored factor returns finite numbers for a matrix that is not positive definite. At s = 10⁶, κ = 10²⁰, its LLᵀ is 5.9·10⁴⁴ from A.
  - *Where no pivot is floored, the two factors agree bit for bit.*
- *Two edges.*
  - A finite matrix can overflow a pivot while L stays finite: rank one at 10²⁰⁰ with one entry 2ε off reads `LogDet` 405.3, and its Σ⁻¹ has four non-finite entries.
  - L⁻¹ can overflow while L and every pivot are finite. The tests construct one.
  - `Status.Finite` catches the first edge and not the second.
- *How `DenseCholesky` treats the same inputs.* It stops at the first pivot that is not finite and positive and names it in `FailedAt`. `LogDet` is NaN, and every read, `RelativePivot` included, throws.
- *What the callers did with `Status`: none read it.*
  - `GaussianComponent` wrote its inverse and log-normalizer from any factor.
  - `BhattacharyyaCoefficient` read `LogDet` and the inverse.
  - The scene export wrote `WriteLTo`'s factor as it came: the edge above exported a finite factor whose LLᵀ is not the covariance.

**The build.**
- *The refusal.* `Decompose` fails when a pivot before flooring, an entry of L or an entry of L⁻¹ is not finite.
  - `FailedAt` (new) names the first such row, and `LogDet` is NaN.
  - `Solve`, `WriteInverseTo`, `WriteLTo` and `Sample` throw, and they also throw before any `Decompose`.
  - A floored factor that stays finite is read as before.
- *The callers.*
  - `GaussianComponent.UpdateCache` throws and leaves `CovarianceInverse` and `LogNormalizationFactor` NaN, so its density reads NaN and `Sample` is refused.
  - `BhattacharyyaCoefficient` refuses by name (`cov1`, `cov2`, `a`, `b`, or the averaged covariance).
  - The scene export refuses the component, naming its layer and index.
- *The two factors now fail one way on a factor that is not finite.* They still differ on a floored factor that stays finite (decision 1).

**The gates.**
- *The reproduction.*
  - On HEAD's code, 8 of the 11 new linear-algebra tests fail, every one that asserts a refusal. All 3 caller tests fail too: the component's cache, the distance (which read a finite −1,156 for `LogDet`), and the export of a meaningless factor.
  - The 3 that pass on HEAD hold what the build leaves alone: floored factors that stay finite are read, and the pinned comparison.
  - On the build, all of them pass.
- *Bit for bit.* The class's existing tests hold, `DenseCholeskyTests`' exact agreement included. The new pinned test compares 63 matrices with the earlier class kept verbatim: the factor, `LogDet`, the inverse, a solve and a sample, bit for bit.
  - 60 are random SPD matrices at scales 2⁻⁸⁰ to 2³⁹; the test requires at least 10 of them to floor.
  - The other 3 are two floored kernels and 10⁻¹³·I.
- *Whole projects*, from a snapshot of HEAD plus the six files:
  - `Maths.LinAlg.Tests`: 179 of 179.
  - `TDA.Mapper.Tests`: 28 of 28 (`MapperGMM` reaches the component).
  - `VizCore.Tests`: 699 passed, 2 skipped. The skips are the Möbius golden drafts, skipped before this slice.
  - `Maths.Regression.Tests`: 3,313 of 3,313 (27 min).
  - The clustering projects (SPC, cumulative field) do not reach the class: their selection is empty and they have no pinned tier. Both build.
- Receipts are in `artifacts/shared-code-defects/receipts/` (`slice2-*`).

**Defects found.**
- None in the slice's code.
- *The defect left by design* (decision 1): a floored factor that stays finite is read, and its `LogDet` and solves can carry order-one error. `Status.FlooredPivots` is its only signal, and no caller reads it.
- *Where GMM leans on the floor.*
  - `CovarianceRegularization` may be 0.
  - `InitializeWithParameters` takes covariances as given.
  - By the arithmetic, not measured: the 10⁻⁶ ridge falls below the pivots' rounding at large covariance scales.
  - In each, a singular covariance is floored, stays finite, and is wrong.

**Decisions for the user.** All three are unbuilt.
1. *The floor's rule.*
   - *Measured.* The absolute 10⁻¹² fires on scale alone, and flooring can introduce order-one error. A pivot at the floor is unchanged. Wherever flooring raised a pivot by orders of magnitude, the outputs carry the errors above: at κ = 1, 10⁻¹³·I has every variance multiplied by ten.
     - A floor relative to the largest diagonal would make the class scale-invariant: each κ column of the ladder would read the same at every s. Once it fires, it is just as wrong.
     - Pivots negative beyond rounding (the stored matrices at κ ≥ 10¹⁸) are floored silently today.
   - *Recommendation: no floor.* Report a pivot that is not positive, as `DenseCholesky` does, and let callers regularize. If a floor stays, make it relative to the largest diagonal, and refuse pivots below −n·ε·max A_jj.
     - Corrected 2026-11-13 (review in `discussion/astra-shared-code-cleanup.md`): not every caller regularizes.
     - GMM's EM adds its ridge, but its warm start and a regularization of 0 take covariances as given.
     - The export's ε fills only the display dimensions a layer lacks; it copies a supplied block as given.
     - Under no floor, a singular covariance in any of these is refused. Those refusals are acceptance cases.
   - *What changes.* GMM's numbers move only where the floor fires today, in the cases above. There `UpdateCache` would throw instead of returning a floored fit.
2. *A pivoted factor in the linear algebra.*
   - *Measured.* Diagonal pivoting, stopped at 10⁻¹² of the largest diagonal, factors every kernel of the reproduction at rank 8 to 57. Its ‖A − LLᵀ‖_F/‖A‖_F is at most 1.3·10⁻¹³, on the same kernels where the floored factor overflowed or floored up to 398 pivots.
   - *Recommendation.* Promote `BgpDrawnFunction`'s factor (Harbrecht, Peters and Schneider 2012) into `Maths.LinAlg` as the semidefinite factor, with a low-rank solve and sample. BGP's drawn function would be its first reader. It is the tool for smooth kernels' covariances; GMM does not need it.
3. *Retiring `CholeskyDecomposition` into `DenseCholesky`.* This is a removal, so it is its own question.
   - *What is lost:*
     - the soft floor (decision 1);
     - `WriteInverseTo`'s Σ⁻¹ through L⁻¹, which GMM's `CovarianceInverse`, its Mahalanobis distances and its mode ascent read;
     - `Sample`'s μ + Lz (`DenseCholesky` offers L⁻ᵀz);
     - `WriteLTo` (`DenseCholesky` offers `Factor(i, j)`).
   - *What moves:*
     - the callers: `GaussianComponent`, `BhattacharyyaCoefficient` (both overloads) and the scene export;
     - the tests: `DenseCholeskyTests` (2 comparisons), `FactorReferenceTests` (3 sites), `BandCholeskyTests` (3) and `CholeskyDivergenceTests`;
     - the cost profiles `CostProfileKernels`, `CostProfileLevers` and `CostProfileEngines`, which time the class as BGP's earlier path.
   - *Bits.* `DenseCholesky`'s factor, `LogDet` and `Status` equal this class's bit for bit where no pivot is floored.
     - An inverse writer that ports this class's L⁻¹ path keeps GMM's Σ⁻¹ bit for bit, where one rebuilt from solves would not.
     - It must refuse an inverse that overflows: `[1e-320]` has a finite factor and `LogDet`, but its reciprocal overflows.
     - Reading Bhattacharyya's mean term through `InverseQuadraticForm` would move its distances in their last bits.
   - *Recommendation.* Retire it once decision 1 removes the floor, and give `DenseCholesky` an inverse writer and an Lz sample so GMM loses nothing. Keep it if the floor stays: the floor is then the whole difference between the two.

**Doc changes owed.**
- *`docs/Maths/LinearAlgebra.md`:*
  - the API list (line 95): `CholeskyDecomposition` gains `FailedAt`;
  - "Dense Cholesky (`CholeskyDecomposition`)" (line 193): the refusal, its three triggers, and the reads that throw, before a first `Decompose` too. Also the floored factor that stays finite: it is read, and flooring can introduce order-one error, with the ladder's numbers; and the floor fires on scale alone.
- *`docs/Clustering/GaussianMixtures.md`:*
  - line 45: `UpdateCache` throws on a covariance with no finite factor, leaving Σ⁻¹ and the log-normalizer NaN;
  - the warm start (line 89): a singular warm-start covariance is floored and read, or refused when its factor is not finite;
  - "Bhattacharyya distance" (line 162): it refuses a covariance with no finite factor, by name.
- *The scene export's `cholesky_l`:* no doc in `docs/` describes it.

### 2026-11-13 — Slice 1: the weight reader at extreme weights

**Commit.** `f50753f` fix(regression): the weight reader keeps every finite weight's normalized log.

**What changed.**
- *`PredictiveChecks.Rows`* reads every finite weight to full precision.
  - Where the sum is finite and W_d/ΣW a normal double, it takes ln(W_d/ΣW), the expression it took before. Ordinary weights keep their bits by construction.
  - Where the sum passes `double.MaxValue`, it takes the same ratio of the weights scaled by the power of two that brings the largest into [1, 2). The scaling is exact.
  - Where the ratio is still below 2⁻¹⁰²², it takes ln W_d − e ln 2 − ln(the scaled sum), a difference of logs. These are the ratios the old form read as −∞ or as an imprecise subnormal.
- *`PointwiseLogDensity.CheckWeights`* tests for one positive weight instead of forming a sum that can overflow. It admits and refuses exactly what it did.
- *One shared routine.* WAIC, PSIS-LOO and the conformal layer's nine sites all read `Rows`. No other reader of the pointwise or candidate densities' weights normalizes them, so the conformal folder is unchanged.

**The gates.**
- *Slice 10's criterion*, in `WeightReaderTests`:
  - `[1e300, 1e-300]` reads rows [0, 1] with log weights [0, −1381.551055796427]; `[1e308, 1e308]` reads two weights of exactly ½. No NaN, no lost row.
  - The log weights against an exact reference: the weights summed as integers in units of 2⁻¹⁰⁷⁴, and the ratio taken to at least 62 bits by integer division. The tolerance is (n + 4)·2⁻⁵²·max(1, |reference|).
  - Nine chosen cases, and a seeded sweep of 4,000 vectors: 2 to 41 weights, binary exponents over random bands of [−1074, 1024), with zeros, the smallest subnormals and values within 10⁻³ of `double.MaxValue`. The worst error is 0.15 of the tolerance. 2,694 of the vectors hold a ratio the old form lost.
  - The readings, against log-sum-exp references on the exact log weights: the weighted predictive (`SplitConformal.LogPredictive`, WAIC's per-unit elpd), the CDF scores (`CdfScores.Fixed`), and slice 10's carried CDF error.
  - A finite tilt of about −1381 nats on the dominant draw makes the first input's tiny row material. It moves WAIC's lppd by 0.35 and 0.61 nats, and the carried bound's weights to .62 and .38.
  - Positive rescaling: 4,000 vectors scaled by a power of two. The 1,353 whose weights and ratios all stay normal agree bit for bit, and the rest within the tolerance. Scales that are not powers of two agree within three tolerances, and the carried error holds under 2^±20.
  - Zero weights drop out and leave the others' log weights bit for bit (2,000 vectors).
- *The reproduction.* On HEAD's code 12 of the 18 new tests fail, each extreme input reading −∞. After the repair, 18 of 18 pass.
- *The ordinary range.* The new pinned test compares 4,000 ordinary vectors with the old form, kept verbatim, bit for bit. The `[1, double.Epsilon]` receipts also keep their bits: that ratio is exact, and the difference of logs gives the same ln ε.
- *The whole regression assembly*, pinned tier included, from a `git archive HEAD` snapshot plus the three files: 3,306 passed, 0 failed (29 min). The selection by change was the whole assembly less its slow tier, since `Pointwise.cs` is a hub, so the slice's run covered it.
- Receipts are in `artifacts/shared-code-defects/receipts/` (`slice1-*`).

**Defects found.**
- None in the slice's code.
- *Finding, the same class outside the brief's files.* BASS's draw-weight readers normalize linear weights by a sum formed in linear space: `Draws3.Select` (`spline/bass/posterior/BassPrediction.cs`) and `PosteriorSummary.Over` (`spline/bass/presentation/`). Under `[1e308, 1e308]` the sum is +∞, and `Draws3.Select` returns weights of 0. They read `BassFit.DrawWeights`, not the pointwise densities. The samplers never produce such weights (SMC's and the flat histogram's are self-normalized), but `fit with { DrawWeights = … }` can.

**Decisions for the user.**
- *Whether BASS's draw weights get the same full-range reading.* If yes, a linear-space twin of `Rows`' rule serves `Draws3.Select` and `PosteriorSummary`, and prediction, sensitivity and the summaries read any finite weights. If no, `BassFit.DrawWeights` states that its sum must stay below `double.MaxValue`, and those readers refuse one that does not.

**Doc changes owed.**
- *`docs/Maths/Regression/Regression.md`*, the weighted draws (line 214): `PointwiseLogDensity.Weights` takes any finite, non-negative weights with one positive, over the whole finite range. The shared reader's rule: the ratio's log where it is normal, so ordinary weights read as before; the power-of-two-scaled ratio past `double.MaxValue`; a difference of logs below the normal range.
- *`docs/Maths/Regression/ConformalPrediction.md`*, "Weighted draws work as they are" (line 76): the shared rows read extreme weights to full precision. This closes slice 10's parked weighted-row follow-up.
