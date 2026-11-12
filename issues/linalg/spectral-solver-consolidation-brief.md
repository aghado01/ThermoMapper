# Spectral solver consolidation: LinAlg to regression

Queued 2026-11-12 for the next available LinAlg and regression owners. This brief brings together three sources:
- the [eigensolver performance brief](eigensolver-performance-brief.md) of 2026-11-09, whose witnesses and items it carries;
- the [Gemini notes](../regression/discussion/gemini-linalg-regression-notes.md) checked in its November 11 follow-ups;
- the 2026-11-12 review of every path from `Maths.LinAlg` into `Maths.Regression`, with two rounds of adversarial review.

This is now the working brief, and the earlier one stays as its record.

Most of the redundancy regression pays for does not come from having several eigensolvers. It falls into four kinds:
- the same matrix decomposed several times within one BASS sensitivity readout;
- per-draw solves whose input depends only on the draw's model structure;
- callers rescaling their matrices or choosing their own tolerances because `EigenFast`'s stopping rule is absolute;
- copies and reference solves in the LinAlg interface.

A fifth finding is a correctness defect in the matrix Poincaré check, whose repair also removes most of its cost.

## Where regression reaches LinAlg

| Route | Production entry | Regression callers |
|---|---|---|
| Dense symmetric, full spectrum | [DenseEigen](../../src/maths/linalg/DenseEigen.cs) → [EigenFast](../../src/maths/linalg/EigenFast.cs) | `RidgeSpectrum` ([CollapsedGPrior.cs:1033](../../src/maths/regression/gprior/CollapsedGPrior.cs)), [PriorPencil.cs:82](../../src/maths/regression/gprior/PriorPencil.cs), [PenalizedSpline.cs:185](../../src/maths/regression/spline/baps/PenalizedSpline.cs), [TensorPenalty.cs:160](../../src/maths/regression/spline/baps/TensorPenalty.cs), [BassRandomIntercepts.cs:190](../../src/maths/regression/spline/bass/BassRandomIntercepts.cs), and every caller of [ScaledEigen](../../src/maths/regression/subspaces/ScaledEigen.cs): E8's `SubspacePosterior`, `RankPosterior`, `CoActivePosterior`, `Whitening`, `MatrixPoincare`, BASS's `BassActiveSubspace` and `BassProjection` |
| Reference Jacobi | [Eigen](../../src/maths/linalg/Eigen.cs), meant for `EIGEN_REFERENCE` builds and oracles | [HierarchicalPSpline.cs:290](../../src/maths/regression/spline/baps/HierarchicalPSpline.cs), once per model; LOBPCG's Rayleigh–Ritz step ([Lobpcg.cs:304](../../src/maths/linalg/Lobpcg.cs)), twice per iteration |
| Thin SVD | [ThinSvd](../../src/maths/linalg/ThinSvd.cs), one implementation | `DrawSubspace.FromFactor` ([SubspacePosterior.cs:51](../../src/maths/regression/subspaces/SubspacePosterior.cs)), `OutputPosterior` ([OutputSide.cs:99](../../src/maths/regression/subspaces/OutputSide.cs)), [BassDecomposition.cs:389](../../src/maths/regression/spline/bass/output/BassDecomposition.cs) and its leave-one-out loop (line 525), [BassRandomEffects.cs](../../src/maths/regression/spline/bass/BassRandomEffects.cs) lines 493, 597 and 685, [CalibrationCovariance.cs:243](../../src/maths/regression/spline/bass/calibration/CalibrationCovariance.cs) |
| LOBPCG | [LOBPCG.Solve](../../src/maths/linalg/Lobpcg.cs) | none; graph code only ([GraphSpectralOperators.cs:200](../../src/graphs/spectral/GraphSpectralOperators.cs), [MagneticLaplacianOperator.cs:273](../../src/graphs/spectral/MagneticLaplacianOperator.cs)) |

Changes to `DenseEigen` also reach callers outside regression. Their results move with workstream B:
- [SpectralMath.cs](../../src/maths/linalg/SpectralMath.cs) lines 35 and 49, which graph spectral code uses;
- [ICA.cs:163](../../src/maths/linalg/ICA.cs) and [MoMPCA.cs:133](../../src/maths/geometry/dim-reduction/MoMPCA.cs);
- [NewtonMetropolis.cs:84](../../src/maths/samplers/blocks/NewtonMetropolis.cs);
- [AlignmentEigenCentrality.cs:114](../../src/graphs/observables/AlignmentEigenCentrality.cs) and [AlgebraicConnectivity.cs:45](../../src/graphs/observables/AlgebraicConnectivity.cs).

`ThinSvd` also serves `Pca`, `PrincipalAngles`, `DiagonalPlusLowRank` and the Stiefel and Grassmann manifolds.

The split between thin SVD and Gram eigendecomposition is sound everywhere it occurs. Callers that hold a factor use the SVD. Callers whose Gram is maintained incrementally (`RidgeSpectrum` under `BassDesign`, and `PriorPencil`) use the eigendecomposition. BassDecomposition's per-case leave-one-out SVDs could be rank-one updates, but they compare small singular values against the rank tolerance, and a Gram route would lose that accuracy. They stay as they are.

## Working rules

- **Solver selection.** One production solver serves each call. It is chosen by the matrix's representation, the part of the spectrum requested and the accuracy needed:
  - a dense symmetric matrix that needs every eigenpair goes to the dense solver;
  - a matrix-free operator with a request for a few extremal pairs goes to LOBPCG;
  - a factor goes to the thin SVD.

  Timings settle cost crossovers among suitable methods. The saved timing for `rand-32` is **0.4166 ms** for `EigenFast` ([receipt](../../artifacts/test-runs/EigenFastTimingHarness/20260602_001327954__CompareReferenceAndFastTimings__abc492ffaf3443b8818427cced15c74f/CompareReferenceAndFastTimings.direct-comparison.json)). The Gemini notes' estimate of under 10 µs through order 50 is not a dispatch threshold.
- **Checks.** Routine checks read residuals and the solver's convergence information. A scalar reference solve is kept for a targeted qualification, or for a specific numerical discrepancy run on the fixture that raises it.
- **Tests follow the right design.** A test that states a contract stays. Examples:
  - residual and orthogonality bounds;
  - the infinite matrix Poincaré ratio when covariance falls outside the bound's range;
  - its equality case at exactly 1.

  A test that records incidental behavior is rewritten when a slice moves it. Examples:
  - bits from one summation order;
  - the order of eigenvectors for a repeated eigenvalue;
  - parity with the reference solver's stopping rule.

  Regenerating a pinned expectation is a deliberate act ([CONTRIBUTING.md](../../CONTRIBUTING.md#testing-standards)). Each commit that moves a pin names it and says why.
- **Bit-identical results are the acceptance test for pure reuse only.** Pure reuse means the same call on the same input, or a cache. A change in how a root, a Gram or a posterior is computed is compared at stated tolerances.
- **New code lands beside old.** Superseded code (`ScaledEigen`, PriorPencil's rescaling, the two-matrix `MatrixPoincare.Ratio`) is deleted in a final cleanup commit of its own ([AGENTS.md](../../AGENTS.md)).
- **Removals are owner decisions.** The items listed under [Open decisions](#open-decisions) as removals stay in place until decided.
- **Machine load.** Heavy runs use `scripts/test-select.ps1` at below-normal priority with capped threads. Timings wait for a quiet machine.
- **Each slice ends with a fresh-eyes review**, and its findings are fixed before the report. Reports go under [Reports](#reports) with:
  - the commit, fixture, command and receipt;
  - the pins moved;
  - **Doc changes owed**.

## A. LinAlg correctness

Both fixes are independent of everything else here. Each lands as its own concern.

### A1. LOBPCG residual at the iteration cap

The loop measures the residual on the current block ([Lobpcg.cs:183–189](../../src/maths/linalg/Lobpcg.cs)). It then updates the block (lines 203–223) and can return the updated block with that earlier measurement (line 239).

**Witness:** the symmetric tridiagonal matrix with diagonal `[4, 3, 2]` and adjacent entries 1, largest pair, seed 173, tolerance `1e-12`, one iteration. The result reported 0.76200174. Its returned pair gives 0.34933228.

**Fix:** `ResidualNorm` and `Converged` describe the eigenpairs returned.

**Gate:** apply the operator independently to the returned pairs, on runs that hit the cap and runs that converge. Starting gates: [LobpcgDenseTests](../../tests/VizCore.Tests/LobpcgDenseTests.cs) and the graph solver parity checks in [EigenFastTimingHarness](../../tests/VizCore.Tests/EigenFastTimingHarness.cs).

### A2. EigenFast sweep with no rotations, and convergence information

A pivot is skipped when |a_pq| < tol ([EigenFast.cs:62](../../src/maths/linalg/EigenFast.cs)). The sweep exits when the Frobenius norm of the off-diagonal falls below tol (line 214). Every pivot can be skipped while their combined norm stays above the exit threshold, and then the solver runs to `maxSweeps`.

**Witness:** 32 × 32, diagonal 1 through 32, every off-diagonal `5e-13`, tolerance `1e-12`. Caps of 1 and 256 give identical eigenpairs.

**Fix:** a sweep that performs no rotations ends the solve. Skipping a pivot and stopping the solve then depend on the same test.

**Convergence information.** The dense result reports how the solve ended, as LOBPCG's `Result` already does ([Lobpcg.cs:112](../../src/maths/linalg/Lobpcg.cs)). Both quantities below are measured on the returned state, after the last sweep. A sample taken during a sweep precedes that sweep's rotations, so it cannot describe a sweep that rotated.
- **Sweeps and rotations.** These are the "deterministic check of work performed" the fix needs.
- **off(A) = (Σ_{i≠j} a_ij²)^½ and the largest off-diagonal entry** of the returned diagonalization.
- **`Converged`** holds when every off-diagonal entry passes the pivot test: below the tolerance for EigenFast, at or below it for the reference. It can hold while off(A) exceeds the tolerance. On the witness every pivot passes while off(A) ≈ 1.57·10⁻¹¹.

**Counting work.** Equal outputs also pass while duplicate work remains. `DecompositionLog` records each decomposition at its executing kernel (reference Jacobi, EigenFast, `ThinSvd`), with the matrix's shape, inside a scope on the async flow that opened it. Counting at the kernel rather than at `DenseEigen` records each solve once, and separate kernels make "zero reference solves" checkable. C's fixture counts the decompositions its reuse removes. B3 removes the reference solves from production code, and B4 makes the compiler enforce it.

**Tests.** [EigenFastTests](../../tests/VizCore.Tests/EigenFastTests.cs) currently compares `EigenFast` against the reference in all seven of its tests. They become property tests:
- residual ‖AV − VΛ‖ and orthogonality ‖VᵀV − I‖ on the existing fixtures;
- the convergence information on the witness, and on inputs that need rotations.

One parity test against the reference remains as a qualification.

## B. One production dense solver

### B1. A scale-relative stopping rule as `DenseEigen`'s only rule

`EigenFast` stops on an absolute tolerance, and callers work around it in different ways:
- [ScaledEigen](../../src/maths/regression/subspaces/ScaledEigen.cs) divides by the largest entry and solves at `1e-14`. It serves all of E8 and BASS sensitivity. [ActiveSubspaces.md](../../docs/Maths/Regression/ActiveSubspaces.md) (the "Eigendecompositions at unit scale" paragraph) records the evidence behind it: QᵀQ − I at 3.8·10⁻³ at small scales before the change, 1.8·10⁻¹¹ after. That paragraph names LinAlg's other consumers as a follow-up.
- [PriorPencil.cs:61–82](../../src/maths/regression/gprior/PriorPencil.cs) divides the whitened matrix by its largest diagonal and solves at `1e-15`.
- [MoMPCA.cs:125–133](../../src/maths/geometry/dim-reduction/MoMPCA.cs) passes `1e-12` times the largest absolute entry.
- `RidgeSpectrum` and `BassRandomIntercepts` pass a raw `1e-15`.

Applied inside the solver, a relative rule retires `ScaledEigen`, PriorPencil's division, MoMPCA's computed tolerance and the raw constants. RidgeSpectrum's column standardization ([CollapsedGPrior.cs:1029–1032](../../src/maths/regression/gprior/CollapsedGPrior.cs)) defines the normalized ridge prior. It is statistics, so it stays with that caller, as do each caller's input checks.

**Scope is an owner decision ([Open decisions](#open-decisions) item 3).** Whether the relative rule replaces the absolute one for every caller, or callers migrate to it one at a time, is open. Pinned tests are no reason to keep an option, since they are rewritten where a slice moves them. But replacing the rule changes the meaning of exposed tolerances. `AlignmentEigenCentrality.Tolerance` ([AlignmentEigenCentrality.cs:68](../../src/graphs/observables/AlignmentEigenCentrality.cs)) is public and documented as an absolute bound.

**The numerical decision is which relative rule ([Open decisions](#open-decisions) item 4):**
- **Relative to the matrix norm**, |a_pq| ≤ tol·‖A‖_F, as LAPACK does. Eigenvalue errors are about tol·‖A‖ in absolute terms.
- **Element-wise**, |a_pq| ≤ tol·√|a_pp·a_qq| ([Demmel and Veselić 1992](https://epubs.siam.org/doi/10.1137/0613074)). Small eigenvalues of graded positive definite matrices keep their relative accuracy, which is what E8's spectra and the prior pencil read. Indefinite matrices with zeros on the diagonal still need a way to stop.

**The floor cannot be a plain maximum with ε‖A‖.** Take a unit diagonal entry beside the independent block 10⁻²⁰·[[1, 0.5], [0.5, 1]]. A pivot floor of ε‖A‖ skips the block's only rotation and returns two eigenvalues of 10⁻²⁰, where the correct values are 0.5·10⁻²⁰ and 1.5·10⁻²⁰. The whole-matrix residual still looks excellent. Whatever stopping guard the rule takes for indefinite inputs must pass this witness. The recommendation is the element-wise rule with such a guard.

**Qualification:**
- The witness above. Gate on the relative error of the weak eigenvalues as well as on residuals.
- Rescale fixed matrices across many orders of magnitude, very small to very large, and check residuals in the original units.
- Rerun E8's equivariance sweep.
- Compare pencil evidence and ridge shrinkage with the current paths.
- Check the graph and geometry callers listed above at their tests' tolerances.

### B2. The interface: caller-owned buffers, eigenvalues only, one permutation

- **Copies.** A `double[,]` input is copied four times:
  1. flattened into a pooled column-major buffer ([EigenFast.cs:291–318](../../src/maths/linalg/EigenFast.cs));
  2. copied again into a pooled working matrix (line 43);
  3. its eigenvectors written into pooled storage;
  4. copied out to jagged rows (lines 276–281).

  Then callers transpose the rows back into `double[,]` columns: `SubspacePosterior` twice (lines 17–31 and 223–231), [RankPosterior.cs:211–214](../../src/maths/regression/subspaces/RankPosterior.cs) and [CoActivePosterior.cs](../../src/maths/regression/subspaces/CoActivePosterior.cs) lines 163–170 and 196–208. `ScaledEigen` adds a scaled copy and a LINQ array. `EigenFast.DecomposeSymmetricInPlace` (line 17) already writes into caller-owned buffers, but no regression caller uses it. `DenseEigen` gains that entry, writing eigenvectors in the layout callers read, and the transposing callers migrate to it.
- **Eigenvalues only.** [PenalizedSpline.cs:185](../../src/maths/regression/spline/baps/PenalizedSpline.cs) reads only eigenvalues, and so does D's reduced solve. A mode that skips the eigenvector rotations saves about half of each rotation's work.
- **Sort.** The final exchange sort ([EigenFast.cs:227–247](../../src/maths/linalg/EigenFast.cs)) swaps whole eigenvector columns, up to O(n²) swaps of length n. Sorting indices first and permuting once, during the copy the solver already makes, replaces it. The order of eigenvectors for an exactly repeated eigenvalue is not a contract. Tests that care compare subspaces.

The end-to-end savings from the eigenvalues-only mode and the sort are unmeasured. They ride with this slice rather than standing alone.

### B3. A production Rayleigh–Ritz solve in LOBPCG

The Rayleigh–Ritz step ([Lobpcg.cs:292–313](../../src/maths/linalg/Lobpcg.cs)) calls the reference `Eigen.DecomposeSymmetric` on the projection of order k, and again on the trial projection of order up to 3k, every iteration. It moves to the production kernel through B2's buffer entry, with workspace reused across iterations.

Measure after A1, and record:
- time and allocations for the projection and the whole solve;
- operator applications, iterations and convergence;
- independently recomputed residuals;
- eigenvalues, orthogonality and subspace agreement for clustered eigenvalues;
- runtime, build, tiering and PGO, warm-up and hardware SIMD support.

Include the sparse graph fixtures at k = 8 and a dense control. The earlier review's probes gave the reference solve 22% and 31% of whole-solve time. Calling the allocating `EigenFast` overload directly made both slightly slower, which is why the buffer entry comes first.

[HierarchicalPSpline.cs:287–292](../../src/maths/regression/spline/baps/HierarchicalPSpline.cs) also calls the reference directly, once per model, for the design Gram. It migrates in the same slice and is qualified against [HierarchicalPSplineConditionalTests](../../tests/regression/HierarchicalPSplineConditionalTests.cs).

### B4. The reference solver as a test oracle

Once B3 lands, the only production reference to `Eigen` is `DenseEigen`'s `#if EIGEN_REFERENCE` branch ([DenseEigen.cs:25, 43](../../src/maths/linalg/DenseEigen.cs)). Tests already use it as an oracle in several files:
- [FixtureGoldenHelpers.cs:157](../../tests/VizCore.Tests/FixtureGoldenHelpers.cs);
- [RandomEffectsWhiteningTests.cs:96](../../tests/regression/RandomEffectsWhiteningTests.cs);
- [SphericalBandwidthValidationHarness.cs:168](../../tests/VizCore.Tests/SphericalBandwidthValidationHarness.cs);
- `EigenFastTests`, `EigenDecomposeSymmetricTests`, `LobpcgDenseTests`.

Moving it into shared test code leaves production one dense kernel, and the compiler then guarantees that production makes zero reference solves. This depends on [Open decisions](#open-decisions) item 1. Until that is decided, B3's migration makes zero reference solves hold in default builds, and code review confirms it.

### B5. Reusing operator images across LOBPCG iterations

The trial update already rotates both X and AX ([Lobpcg.cs:216](../../src/maths/linalg/Lobpcg.cs)), and the next iteration applies the operator to X again (line 178). The experiment compares carrying AX across that boundary with applying the operator fresh, using the same projected solver on both sides. It counts operator applications and checks residual drift independently, including on constrained and preconditioned runs. Reuse is adopted if the whole solve benefits and the convergence checks hold.

### B6. One helper for symmetric matrix functions

Five places compute a function of a symmetric matrix through its eigensystem by hand, each with its own O(p³) triple loop:

| Site | Shape | Cutoff policy (stays with the caller) |
|---|---|---|
| `Whitening.SymmetricRoot` ([Whitening.cs:59–93](../../src/maths/regression/subspaces/Whitening.cs)) | V Λ^½ Vᵀ | refused unless the smallest eigenvalue exceeds 10⁻¹² of the largest |
| `PsdRoot` ([BassActiveSubspace.cs:456](../../src/maths/regression/spline/bass/sensitivity/BassActiveSubspace.cs)) | V Λ^½ Vᵀ | negative eigenvalues clamped to 0; every positive one kept |
| `Root` ([BassActiveSubspace.cs:1076](../../src/maths/regression/spline/bass/sensitivity/BassActiveSubspace.cs)) | rows √λ vᵀ, rectangular | eigenvalues at or below 10⁻¹⁴ of the scale dropped; below −10⁻¹² refused |
| `OutputFactor` (BassActiveSubspace.cs:296) | columns √λ v | eigenvalues above 10⁻¹⁴ of the largest kept |
| `InverseTranspose` ([BassProjection.cs:273](../../src/maths/regression/spline/bass/projection/BassProjection.cs)) | V Λ⁻¹ Vᵀ | none |

A LinAlg helper provides the symmetric form and the factor form from one eigensystem. It takes each caller's cutoff policy as input.

## C. BASS sensitivity: reuse what one readout has already computed

This work is independent of A and B and can run beside them. All sites are in [BassActiveSubspace.cs](../../src/maths/regression/spline/bass/sensitivity/BassActiveSubspace.cs) unless named. The line numbers here predate C's commits.

Two acceptance checks apply to each item:
- **Outputs:** a bit-level comparison of readouts before and after, or a stated tolerance where an item says so.
- **Work:** a fixture that counts the decompositions removed, through `DecompositionLog`. Equal outputs alone also pass while the duplicate work is still present.

### C1. The adequacy posterior is built twice

`Adequacy` (line 213) calls `Read` (line 65). On the direct fit, `Read` builds `OutputPosterior.Of(factors, weights, k)` (line 140): an SVD of each draw's factor, then one of the combined weighted factors. [BasisAdequacy.Of](../../src/maths/regression/subspaces/OutputSide.cs) (line 201) then builds it again from the same factors and the same k (line 206).

**Fix:** pass `Read`'s posterior to `BasisAdequacy.Of`, which removes a whole second set of decompositions. This is likely the largest single saving in this workstream.

The two calls currently normalize their weights by different routes: `Read`'s from the draw selection, and the second from `SubspacePosterior`'s normalized weights. So compare patterns, traces, shares and containment at a stated tolerance.

### C2. The output metric W_y is decomposed up to four times per readout

For a declared m × m output metric, one readout decomposes it at:
- `Metric` (line 1070), which calls `Root(w)` only to validate it and discards the result;
- `PsdRoot(w)` (lines 77 and 155);
- `Root(w)` on the sampled route (line 112);
- `MatrixCheckIf` (line 323), which calls `Root(w)` again.

`Adequacy` reruns `Metric` and `PsdRoot` (line 229) after `Read` has already done both, which makes six decompositions. The default metric is the identity: `PsdRoot` and `Root` still decompose it, and `PsdRoot` rebuilds the root in O(m³).

**Fix:**
- Decompose W_y once per readout. Read both roots from that eigensystem, each keeping its own rule. `Root` drops eigenvalues at or below 10⁻¹⁴ of the scale and produces a rectangular factor. `PsdRoot` keeps every positive eigenvalue and produces a symmetric root in output coordinates.
- Add a fast path for identity and diagonal metrics that reproduces the general path's rows and row order (descending eigenvalue).
- Pure reuse: bit-identical.

`Root(ComponentMetric(fit, w))` (line 184) decomposes a different matrix and stays.

### C3. Cache OutputFactor's term factor per structure

`OutputFactor` (lines 296–314) eigendecomposes K = `GradientGram(terms, W_x)` for every draw. It also rebuilds W_x = LLᵀ (lines 299–302) every time. K depends only on the draw's ordered terms and on the frame, and draws repeat structures.

**Fix:**
- Hoist W_x out of the per-draw call.
- Cache the factor per structure under the existing ordered `Structure` comparer (line 390), which matches terms in sequence by reference.
- Follow `MatrixCheck`'s cache (lines 342–358): a pure cache, sized by a fixed budget and never by hit history, whose entries are what a miss computes.
- Bit-identical.

### C4. Each sampled draw's Gram is formed twice

On the sampled route, `WeightedSums` accumulates C together with its Monte Carlo error (lines 274 and 290, class at line 1206). [DrawSubspace.FromFactor](../../src/maths/regression/subspaces/SubspacePosterior.cs) then forms FᵀF again (lines 43–50) before its thin SVD. That is O(N·p²) twice per draw.

**Fix:**
- `FromFactor` takes the accumulated C.
- The thin SVD stays, as does the fallback to C's own eigendecomposition when there are fewer right singular vectors than the requested dimension (line 54). Exercise that fallback explicitly.
- The summation order of `DrawSubspace.C` changes, so compare at a tolerance.

### C5. The whitened metric's inverse

`Whitening.SymmetricRoot` eigendecomposes Σ to build L = Σ^½ (line 969). [BassProjection.InverseTranspose](../../src/maths/regression/spline/bass/projection/BassProjection.cs) (line 273) then eigendecomposes L to get L⁻¹ (line 286), from the same eigenvectors. The matrix is p × p, so this is about consistency more than cost. Do it with B6, where the helper returns the root and its inverse from one eigensystem.

## D. The matrix Poincaré check in term space

**Contract.** `MatrixPoincare.Ratio` ([OutputSide.cs:135](../../src/maths/regression/subspaces/OutputSide.cs)) returns max_v vᵀAv / vᵀBv. That is the smallest c with Cov ≼ c·Bound, and it is infinite where the covariance has mass outside the bound's range.

**Current path.** BASS calls it once per draw (BassActiveSubspace.cs:361), with:
- covariance = M Σ_B Mᵀ and bound = M K_D Mᵀ (lines 359–360);
- M = R·H, of size rows × t, where R is the output metric's root and H maps the t terms to outputs;
- Σ_B and K_D cached per structure (lines 349–358).

Each call:
1. eigendecomposes the bound in output space (line 146) and keeps eigenvalues above 10⁻¹² of the largest (line 150);
2. declares the ratio infinite if covariance mass outside that set exceeds 10⁻¹² of the scale (line 158);
3. builds the reduced matrix from r² quadratic forms of O(m²) each, with r ≤ t (lines 161–165);
4. solves again (line 168).

The check is off by default because "it about doubles the readout" ([RankAndOutputSideTests.cs:391](../../tests/regression/RankAndOutputSideTests.cs)).

**Defect.** The cutoff acts on eigenvalues of M K_D Mᵀ, which scale with M squared. At 10⁻¹² of the largest eigenvalue, it drops any term direction the draw weights below about 10⁻⁶ of its leading one. The second round of review ran these witnesses against the built E8 assembly, with Σ_B = diag(0.25, 1) and K_D = I unless stated:

| M | `ThinSvd` numerical rank | Current ratio | Correct ratio |
|---|---:|---:|---:|
| diag(1, 1) | 2 | 1 | 1 |
| diag(1, 1e-7) | 2 | 0.25 | **1** |
| diag(1, 0) | 1 | 0.25 | 0.25 |
| diag(1, 1e-7), Σ_B = diag(0.25, 2) | 2 | 0.25 | **2** |

For the second row, 0.25 fails the contract along e₂: the bound would need 10⁻¹⁴ ≤ 0.25·10⁻¹⁴. The last row shows the cutoff can hide a violation. The expected values in the table are settled.

**Reduction.** Substituting u = Mᵀv turns the ratio into the t × t pencil (Σ_B, K_D) over range(Mᵀ):
- **M of full column rank** (decided by `ThinSvd`'s numerical rank, on M itself rather than its square): the ratio is that pencil's largest eigenvalue. It depends only on the structure, so it is cached beside Σ_B and K_D.
- **Rank-deficient M:** Q, the right singular vectors of M over its numerical rank, restricts both forms to (QᵀΣ_BQ, QᵀK_DQ) for that draw, at O(rows·t² + t³).

**Cutoffs move to term space; they do not disappear.**
- A direction null in K_D that carries Σ_B mass gives +∞.
- A direction null in both forms drops out, so the ratio stays finite. Σ_B can carry rounding noise in such a direction, because `CovarianceGram` computes E[ab] − E[a]E[b] ([BassMoments.cs:296–303](../../src/maths/regression/spline/bass/sensitivity/BassMoments.cs)).
- Each form is judged against its own largest eigenvalue, near rounding level, and no squared quantity is involved. The acceptance cases fix the constants.

**The output metric enters through its range, not its root.** `MatrixCheckIf` builds M from `Root(W_y)` (BassActiveSubspace.cs:323), which drops eigenvalues at or below 10⁻¹⁴ of the scale (line 1084). So a direction can disappear before any solver runs. Congruence by a factor R of full column rank leaves the ratio unchanged, because range(Rᵀ) = range(W_y). The check therefore depends on W_y only through that range.
- C's probe agrees: the ratios under the identity metric and under a random full metric match to about 10⁻¹³ on every draw.
- The second review's witness, through the BASS caller: W_y = diag(1, 10⁻¹⁶), H = I, K_D = I and Σ_B = diag(0.25, 2). The declared metric gives a ratio of 2. The truncated root keeps one row, [1, 0], and gives 0.25 even with an exact term-space solver.

The check takes M = UᵀH, with U an orthonormal basis of W_y's resolved range, or H itself when W_y is positive definite, and no square roots of small eigenvalues. A diagonal metric's range is exact: the locations with positive weight. That is what lets the witness return 2. For a general formed W_y, resolving its range is a rank decision at the eigensolver's resolution, about m·ε·‖W_y‖.

**Public entry.** `Ratio(termCovariance, termBound, map)` replaces the two-matrix signature. Formed matrices pass map = I, which is how the cosine cases in `MatrixPoincare_EqualityOnCosines_AndAMisScaledBoundFails` (line 322) migrate. Its contract states the resolution for formed inputs: a caller holding a factor passes the factor, since forming the product squares it.

**Acceptance cases:**
1. the existing equality on cosines (1), halved bound (2), mixed field (< 1) and outside-range (∞) cases;
2. the four witnesses in the table;
3. K_D null with Σ_B mass there (∞);
4. both forms null, with rounding noise in Σ_B (finite);
5. through the BASS caller, W_y = diag(1, 10⁻¹⁶) with H = I, K_D = I and Σ_B = diag(0.25, 2): ratio 2. A direction must not disappear in the metric's root;
6. through map = I, the formed pair Σ = diag(0.25, 10⁻¹⁴) and K = diag(1, 10⁻¹⁴): ratio 1. A direction must not disappear at the term-space cutoff, which bounds that cutoff's constant;
7. `MatrixPoincare_HoldsOnFittedDraws` (line 375), no violations, with the per-structure cache hit and missed in the same run.

**What the check measures in BASS.** For a product law the theorem gives Σ_B ≼ K_D exactly for every structure. So with exact moments, every ratio is at most 1, and a ratio above 1 can come only from wrong moments or a wrong Poincaré constant. In BASS the check therefore verifies the moment service and the constants. When M has full column rank, a draw's ratio is its structure's. Two rank-deficient draws can share terms and still see different ratios, so results are cached per structure only for full-rank draws, and reporting stays per draw ([Open decisions](#open-decisions) item 5).

## E. E8 leading output patterns: the one regression candidate for LOBPCG

`OutputPosterior.Of` already avoids forming the dense output covariance. Its `patterns` argument limits the pattern vectors it returns, while `MeanEigenvalues` exposes the full nonzero spectrum and `MeanTrace` sums it.

Benchmark a separately requested leading-patterns readout for large output dimension and small requested rank:
- Apply the mean operator as v ↦ Σₛ wₛ Fₛ(Fₛᵀv).
- Accumulate its trace as Σₛ wₛ ‖Fₛ‖²_F.
- Witness: two equally weighted draws with covariances diag(4, 3, 0) and diag(0, 3, 4). Their mean is diag(2, 3, 2), so the second coordinate leads even though it ranks second in each draw.
- Compare with the current SVD on this case and on clustered leading eigenvalues. Check the trace, Ritz residuals and leading subspace angles.
- Record end-to-end time, allocations and operator applications as factor width, output dimension and requested rank vary.

The analogous input-side experiment applies Gᵀ(Gv)/M for M equally weighted gradient rows, when only the first few directions are requested. E8's thin SVD still serves requests that resolve weak directions. Exact g-prior evidence and its spectral sums need the full pencil spectrum. Adoption is measurement-gated. [Output-side tests](../../tests/regression/RankAndOutputSideTests.cs) · [Prior pencil](../../src/maths/regression/gprior/PriorPencil.cs)

## F. Factorization follow-ups carried from the November 11 review

- **GP regression.** Give the ordinary fit a dense factor that solves rather than inverts.
  - The legacy `CholeskyDecomposition` constructs L⁻¹ on every fit. The GLS mean and the marginal evidence need only two covariance solves and a log determinant.
  - Compute the precision entries needed by the leave-out path on demand.
  - Preserve the GP's rejection of pivots at or below 1e-12, including positive pivots that the legacy factor floors.
  - Check fits, predictions, and singleton and grouped leave-out log densities, with a near-singular fixture for the rejection boundary.
  - Measure ordinary fitting and leave-out evaluation separately.

  [GpRegression](../../src/maths/regression/bgp/GpRegression.cs) · [Cholesky.cs](../../src/maths/linalg/Cholesky.cs) · [DenseCholesky](../../src/maths/linalg/DenseCholesky.cs)
- **Anisotropic tensor splines.** Benchmark against the caller's current compact, diagonal-major banded LDLᵀ factor.
  - The proposed 2–4× gain from `SpdFactors.ForBand` is unmeasured.
  - Compare end-to-end REML time and allocations across tensor sizes and bandwidths, including any conversion from compact storage.
  - Check coefficients, log determinants, REML optima and the handling of infeasible smoothing points before changing dispatch.

  [AnisotropicTensorSpline](../../src/maths/regression/spline/baps/AnisotropicTensorSpline.cs) · [WideBandCholesky](../../src/maths/linalg/WideBandCholesky.cs) · [WideBandCholeskyTests](../../tests/linalg/WideBandCholeskyTests.cs) · [AnisotropicRemlTests](../../tests/regression/AnisotropicRemlTests.cs)

## G. Kernel follow-ups: FMA, layout and the thin SVD

### G1. Fused multiply-add

**Current state.**
- `DenseEigenOptions.FastVariant = Fma` turns on the x86 `Fma` intrinsics in `EigenFast`, gated on `Fma.IsSupported` ([EigenFast.cs:35](../../src/maths/linalg/EigenFast.cs)). They apply in the vectorized eigenvector update (lines 94–109) and in the vectorized segment of A above row p (lines 138–153).
- The scalar tails and the two strided segments (lines 125–131 and 170–205) stay unfused.
- Only the timing harness ([EigenFastTimingHarness.cs](../../tests/VizCore.Tests/EigenFastTimingHarness.cs) lines 103, 130 and 353–370) and one dispatch test ([SpectralDenseDispatchTests.cs:85](../../tests/VizCore.Tests/SpectralDenseDispatchTests.cs)) select it.
- It was never enabled, pending validation, and no saved timing for it exists.

**Assessment.**
- **Speed:** limited in the current layout. Fusing turns each c·a − s·b into 2 operations instead of 3, but only in the vectorized segments. The two scalar segments that walk along rows p and q read one element from each column, so most reads land on a different cache line. They cover about two-thirds of A's elements per rotation on average, and they probably dominate the time.
- **Accuracy:** one rounding instead of two per element leaves Jacobi's error bounds where they are. The relative accuracy B1 is after comes from the stopping rule, not from fusing.
- **Reproducibility:** the variant breaks the project's rule in two ways.
  - On a CPU without FMA it silently falls back, so the same input gives different bits on different machines.
  - Within one rotation, some elements are fused and others are not, depending on n mod 4 and p.

  `ThinSvd` avoids FMA explicitly so that "vector width does not change the bits" ([ThinSvd.cs:212](../../src/maths/linalg/ThinSvd.cs)).

**The project's convention.** `DenseProduct` fuses every operation, in the vectorized body and the scalar tail alike, with `Vector256.FusedMultiplyAdd` and `Math.FusedMultiplyAdd` and no switch ([DenseProduct.cs](../../src/maths/linalg/DenseProduct.cs) lines 108 and 139–167). Because `Vector256`'s software fallback also fuses, every entry matches the scalar loop bit for bit on any machine ([LinearAlgebra.md](../../docs/Maths/LinearAlgebra.md), DenseProduct's arithmetic). The conformal layer follows the same convention. FMA enters a kernel this way, as its normal arithmetic, not as an option.

**Plan.**
- Remove the `Fma` option ([Open decisions](#open-decisions) item 2).
- Adopt fused arithmetic, under the convention, in kernels whose loops are contiguous:
  - `ThinSvd`'s rotations now;
  - `EigenFast` after G2.
- For each kernel, one timing at regression sizes decides: fused against unfused, with residual and orthogonality checks. Regression sizes are input-space Grams up to about 100 and output dimensions up to the output-location counts. Adoption moves bits, so pins are regenerated deliberately.

### G2. EigenFast's strided segments

Storing only the upper triangle makes the updates along rows p and q (EigenFast.cs:182–205) strided scalar loops. Two directions are possible, both gated by measurement:
- A layout that keeps every update on a contiguous column.
- For positive definite inputs, a Cholesky factor followed by `ThinSvd`'s one-sided Jacobi. Every rotation is then a contiguous column operation, and small eigenvalues keep their relative accuracy (Veselić and Hari; Demmel and Veselić; cited from memory). Positive definite inputs cover the Grams, covariances and the whitened pencil, which are most regression callers.

Inputs that are indefinite or singular keep two-sided Jacobi: the symmetrized co-active matrix, singular positive semidefinite matrices such as the Poincaré bounds and metric roots, and the Hessians in [NewtonMetropolis](../../src/maths/samplers/blocks/NewtonMetropolis.cs).

### G3. Column norms in ThinSvd

The sweep recomputes both column norms for every pair ([ThinSvd.cs:110](../../src/maths/linalg/ThinSvd.cs)), but a norm changes only when its column is rotated. Near convergence, most pairs are skipped, and those two norms are two of the three dot products each pair costs. Updating the norms through each rotation, and recomputing them at the start of each sweep to control drift, removes that cost. Bits move.

## Sequencing

1. **In parallel from the start:**
   - A1 and A2, with A2 adding convergence information and rewriting `EigenFastTests`;
   - C1, C2 and C3, the pure reuse;
   - D, whose expected values are settled.
2. **After A2:** B1 and B2, which need the convergence information for qualification. C4 can land any time; C5 waits for B6.
3. **After A1 and B2:** B3, then B5. B4 follows [Open decisions](#open-decisions) item 1.
4. **Measurement-gated:** E, F, G1–G3. G1's EigenFast timing waits for G2.
5. **Final cleanup commit:** delete `ScaledEigen`, PriorPencil's division, the two-matrix `Ratio` and the decided removals.

## Open decisions

Owner: Azriel. Slices leave these in place until decided.

1. **Remove the `EIGEN_REFERENCE` build?** What's lost: running the whole system on the reference solver to compare results on a numerical discrepancy. The test oracle covers targeted qualification either way. Affects [DenseEigen.cs](../../src/maths/linalg/DenseEigen.cs), [SpectralDenseDispatchTests.cs:65](../../tests/VizCore.Tests/SpectralDenseDispatchTests.cs), the harness's `#if` blocks, and the docs listed below.
   Recommended by both reviews: retire it once B3 has migrated both production callers, keeping the reference as the test oracle.
2. **Remove `DenseEigenFastVariant.Fma`?** What's lost: the FMA rotation path, which production never selects. Recommended: retire the switch, and let G1's timings decide whether fused arithmetic belongs in either kernel.
3. **B1's scope.** Should the relative rule replace the absolute one for every caller, or should callers migrate to it one at a time? What changes: the meaning of exposed tolerances such as `AlignmentEigenCentrality.Tolerance`.
4. **Which relative stopping rule for B1?** Recommended: element-wise, with a stopping guard for indefinite inputs that passes the 10⁻²⁰ witness.
5. **Matrix Poincaré reporting and default.** Recommended: cache full-rank results per structure; keep reporting per draw, since rank-deficient draws that share terms can differ; keep the check opt-in until whole-readout timings justify changing the default.

Settled in the 2026-11-12 reviews:
- D's expected values, the table's and the two added witnesses';
- the test rule in [Working rules](#working-rules);
- the scope of bit-identical acceptance.

## Doc changes owed

`docs/` is updated in batches once slices land:
- [LinearAlgebra.md](../../docs/Maths/LinearAlgebra.md):
  - the dense-route diagram and "Two dense routes" (stopping rule, reference status);
  - the table of which triangle each routine reads;
  - the notes on EigenFast's segments;
  - the SIMD notes;
  - the `FastVariant` option row.
- [ActiveSubspaces.md](../../docs/Maths/Regression/ActiveSubspaces.md): the unit-scale paragraph, which B1 retires; the matrix Poincaré section's contract, factor form and resolution.
- If item 1 is decided as removal: [Architecture.md](../../docs/Architecture.md) (`ExtraDefineConstants`), [Design.md](../../docs/Design.md) (the EigenFast versus Eigen entry), [README.md](../../README.md) (the reference build command).

## Reports

Each slice is reported here, newest first, with its commit, fixtures, commands, receipts, pins moved and doc changes owed.

### D: the matrix Poincaré check in term space (2026-11-12)

**Commit:** `56ac387`.

**What changed:**
- `MatrixPoincare.Ratio(termCovariance, termBound, map)` replaces the two-matrix signature. Its contract is the least c with MΣMᵀ ≼ c·MKMᵀ.
  - `RowSpace` decides M's rank by its thin SVD, on M itself.
  - A map of full column rank reads the pencil over the whole term space. Otherwise both forms are restricted to range(Mᵀ).
- `PencilRatio` judges each form on its own scale, at 8nε. An eigenvalue of K at or below 8nε of its largest counts as zero, and so does Σ's mass at or below 8nε of ‖Σ‖_F. Formed matrices pass map = I.
- `MetricRange` gives the BASS check W_y's range as rows. It returns null when W_y is positive definite, and it is exact for a diagonal metric. A dense metric's eigenvalues at or below m·ε of the largest count as zero.
- The check maps each draw by UᵀH.
- `TermPencil` keeps each structure's Σ_B and K_D with its whole-space ratio, read once, for draws of full column rank. A rank-deficient draw is read over its own range. The unused `Sandwich` helper is gone.

**Acceptance cases:**
1. The four existing cases (cosines 1, halved 2, mixed field below 1, outside the range ∞): unchanged, now through map = I.
2. The four witnesses in the table: 1, 1, 0.25 and 2.
3. K null with Σ mass there: ∞.
4. Both forms null, with Σ = 3·10⁻¹⁷ there: 1.
5. The metric-range witness, W_y = diag(1, 10⁻¹⁶): `MetricRange` returns null, so the map is H. The ratio is 2, where the old root's single row reads 0.25. It is tested through the metric-to-map step the BASS caller uses, not a full readout, since real moments cannot produce Σ_B = diag(0.25, 2).
6. The formed pair at 10⁻¹⁴ through map = I: 1.
7. Fitted draws: no violations.

Also tested:
- a rank-deficient map off the coordinate axes, M = [1, 1]: 0.625;
- `MetricRange` on a zero weight, a rank-one dense metric and a positive definite one;
- the same draw read twice in one readout, and read alone: the same bits.

**Effect on fitted draws:** the worst ratio is 0.8670 over 280 draws, the same as before to the digits printed. The existing check that ratios are the same read in a batch or alone still passes.

**Gate:** `scripts/test-select.ps1` at edit scope over D's paths, 214 classes: 2,489 passed, none failed. No pinned-tier test reaches an active-subspace readout.

**Not done:** whole-readout timings for [Open decisions](#open-decisions) item 5.

**Pins moved:** none. The four quadrature cases pass map = I.

**Doc changes owed:** [ActiveSubspaces.md](../../docs/Maths/Regression/ActiveSubspaces.md), the matrix Poincaré section. It needs:
- the contract, as the least c;
- the term-space reading and `Ratio`'s factor form;
- the metric read through its range;
- the two resolutions, 8nε for the forms and m·ε for a dense metric;
- the correction itself, with the witness table.

### C1–C3: BASS readout reuse (2026-11-12)

**Commit:** `23f599f`.

**What changed:**
- **C1.** `Adequacy` passes `Read`'s output posterior to `BasisAdequacy.Of`, whose signature is now `Of(factors, outputs, span, dimension)`. It refuses a posterior with a different draw count or more patterns. The unit tests build the posterior through a helper.
- **C2.** `MetricRoots` decomposes W_y at most once per readout, lazily, so a readout that needs no root runs no solve. `Root` and `PsdRoot` read that one eigensystem, each under its own rule and each built once. `Read` returns the metric, and `Adequacy` reads it from there. `PsdRoot`'s outer products skip zero eigenvector components, which leaves every sum unchanged, so an identity or diagonal metric costs O(m²) there rather than O(m³). This replaces a separate fast path.
- **C3.** `TermFactors` caches F_K per ordered structure within a readout, under the matrix check's budget (`GramCacheValues`), and forms W_x once.

**Outputs:** bit-identical. A probe hashed the bits of five readouts before and after:
- a direct exact readout with a full W_y;
- the same with the identity metric;
- a direct sampled readout;
- a basis exact readout;
- adequacy.

Every line matched, adequacy included: its weights, normalized by a different route, rounded the same here. The probe's source and both outputs are in `artifacts/spectral-c/`.

**Work:** [ReadoutDecompositionTests](../../tests/regression/ReadoutDecompositionTests.cs) counts decompositions with `DecompositionLog`:
- W_y is decomposed once on the sampled route with a declared metric. Before: 3.
- Three copies of one draw add two symmetric solves when patterns are asked. Before: 4.
- Adequacy adds only the span's SVD to the direct readout's decompositions. Before: 12 symmetric solves against 10.

Run against the pre-C sources, all three failed with exactly those counts.

**Gate:** `scripts/test-select.ps1 -Scope Commit` over the A and C paths, 286 of 301 classes: 2,920 passed, none failed or skipped ([trx](../../artifacts/spectral-c/test-select/20261112_112227-commit/Commit.trx)).

**Not done here:** C4 (each sampled draw's Gram formed twice) and C5 (the whitening inverse), as sequenced.

**Pins moved:** none.

**Doc changes owed:** [ActiveSubspaces.md](../../docs/Maths/Regression/ActiveSubspaces.md), `BasisAdequacy.Of`'s signature.

### A2: Jacobi stopping, convergence report and DecompositionLog (2026-11-12)

**Commit:** `3661e81`.

**What changed:**
- A sweep that rotates nothing ends the solve, in EigenFast and in the reference alike. The eigenpairs are the same bits as before, since the extra sweeps changed nothing.
- `EigenResult.Convergence` (`DenseEigenConvergence`) reports, measured on the returned state after the last sweep:
  - sweeps and rotations;
  - off(A) and the largest off-diagonal entry;
  - `Converged`, when every off-diagonal entry passes the pivot test.
- `EigenFast.DecomposeSymmetricInPlace` returns the same report. `ScaledEigen` carries it back to the caller's scale.
- `DecompositionLog` records decompositions at the three kernels.

**Tests:** [EigenFastTests](../../tests/VizCore.Tests/EigenFastTests.cs) became property tests:
- residuals, orthogonality and known spectra;
- the no-rotation witness, with one sweep, no rotation, and the same eigenpairs under caps of 1 and 256;
- a capped solve whose reported off(A) matches VᵀAV recomputed from the result;
- one agreement check against the reference.

[DecompositionLogTests](../../tests/VizCore.Tests/DecompositionLogTests.cs) covers scopes, nesting, and work on tasks the scope's flow starts.

**Gate:** the `VizCore.Tests` selection for the A paths, 34 classes: 233 passed, 2 skipped (the Möbius drafts' standing skips). The regression reach through `ScaledEigen` and the kernels ran in C's gate.

**Pins moved:** none. The parity-with-reference assertions were incidental behavior and are now property tests.

**Doc changes owed:** [LinearAlgebra.md](../../docs/Maths/LinearAlgebra.md), the stopping rule, `DenseEigenConvergence` and `DecompositionLog`. `AlignmentEigenCentrality.Tolerance` is documented as a bound on the off-diagonal infinity norm. With the pivot test now ending the solve, that wording describes what EigenFast does. Before, its exit read the Frobenius norm.

### A1: LOBPCG residual at the cap (2026-11-12)

**Commit:** `c2163a8`.

**What changed:** when the loop ends without converging, the operator is applied to the returned block and its residual measured, so `ResidualNorm` and `Converged` describe the pairs returned. That costs one extra operator application, at the cap only.

**Tests:** `LobpcgDenseTests.Solve_ResidualDescribesReturnedPairs` checks the witness (0.3493 reported, not 0.7620), a converged run, and a block of four capped at three iterations, each against the operator applied independently. The fix sits after the loop, whatever the options. Constrained and preconditioned capped runs have no test of their own.

**Gate:** as A2's.

**Pins moved:** none.

**Doc changes owed:** none.
