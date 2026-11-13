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
- **Removals are owner decisions,** recorded under [Decisions](#decisions). A slice removes nothing the owner has not decided to remove.
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

**Scope ([Decisions](#decisions) item 3): callers migrate in stages.** Pinned tests are no reason to keep an option, since they are rewritten where a slice moves them. But replacing the rule changes the meaning of exposed tolerances. `AlignmentEigenCentrality.Tolerance` ([AlignmentEigenCentrality.cs:68](../../src/graphs/observables/AlignmentEigenCentrality.cs)) is public and documented as an absolute bound.

**The rule ([Decisions](#decisions) item 4): one element-wise rule, under three acceptance cases.** The candidates were:
- **Relative to the matrix norm**, |a_pq| ≤ tol·‖A‖_F, as LAPACK does. Eigenvalue errors are about tol·‖A‖ in absolute terms.
- **Element-wise**, |a_pq| ≤ tol·√|a_pp·a_qq| ([Demmel and Veselić 1992](https://epubs.siam.org/doi/10.1137/0613074)). Small eigenvalues of graded positive definite matrices keep their relative accuracy, which is what E8's spectra and the prior pencil read. Indefinite matrices with zeros on the diagonal still need a way to stop.

**The floor cannot be a plain maximum with ε‖A‖.** Take a unit diagonal entry beside the independent block 10⁻²⁰·[[1, 0.5], [0.5, 1]]. A pivot floor of ε‖A‖ skips the block's only rotation and returns two eigenvalues of 10⁻²⁰, where the correct values are 0.5·10⁻²⁰ and 1.5·10⁻²⁰. The whole-matrix residual still looks excellent. Whatever stopping guard the rule takes for indefinite inputs must pass this witness. Decided: the element-wise rule, qualified on the three acceptance cases in [Decisions](#decisions) item 4.

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

Moving it into shared test code leaves production one dense kernel, and the compiler then guarantees that production makes zero reference solves. [Decisions](#decisions) item 1 keeps the build, so B4 is not pursued. B3's migration makes zero reference solves hold in default builds, and `LobpcgDenseTests` checks it.

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

**What the check measures in BASS.** For a product law the theorem gives Σ_B ≼ K_D exactly for every structure. So with exact moments, every ratio is at most 1, and a ratio above 1 can come only from wrong moments or a wrong Poincaré constant. In BASS the check therefore verifies the moment service and the constants. When M has full column rank, a draw's ratio is its structure's. Two rank-deficient draws can share terms and still see different ratios, so results are cached per structure only for full-rank draws, and reporting stays per draw ([Decisions](#decisions) item 5).

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

**Before the switch removal.** Commit `ec86071` preserves both the kernel and its timing harness.
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
- Remove the `Fma` option ([Decisions](#decisions) item 2).
- G1 stays open after the switch removal. Benchmark `ThinSvd`'s contiguous rotations now. When G2 settles the `EigenFast` layout, its report should point back here with the layout to benchmark.
- For each kernel, one timing at regression sizes decides: fused against unfused, with residual and orthogonality checks. Regression sizes are input-space Grams up to about 100 and output dimensions up to the output-location counts. Adoption moves bits, so pins are regenerated deliberately.
- Record the fused and unfused timings, CPU/SIMD path, input sizes and numerical checks in the G1 report, including a result that favors keeping the kernel unfused. A later layout or hardware change that alters the cost of the measured loop gives a concrete reason to rerun that comparison.

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
2. **After A2:** B1 and B2, which need the convergence information for qualification. A1, A2 and all of C and D have landed (see [Reports](#reports)).
3. **After A1 and B2:** B3, then B5. B4 is not pursued ([Decisions](#decisions) item 1).
4. **Measurement-gated:** E, F, G1–G3. G1's EigenFast timing waits for G2.
5. **Final cleanup commit:** delete `ScaledEigen`, PriorPencil's division, the two-matrix `Ratio` and the decided removals.

## Decisions

Decided by the owner on 2026-11-12, after the fourth review.

1. **The `EIGEN_REFERENCE` build stays.** The reference Jacobi remains in the production assembly behind the build flag, so the whole application can still run on it for comparison. B4 is not pursued. In default builds, production makes no reference solve: B3 migrated both callers, and `LobpcgDenseTests` checks LOBPCG's projections by kernel.
2. **The `Fma` switch is removed.** Fused arithmetic returns only through G1, unconditional and bit-stable per kernel, if its timings favor it.
3. **B1 migrates callers in stages.** One kernel supports both stopping rules while callers move, and each stage qualifies its own callers:
   1. `ScaledEigen` and `PriorPencil`;
   2. `RidgeSpectrum`, `BassRandomIntercepts` and `MoMPCA`;
   3. the graph callers, where `AlignmentEigenCentrality.Tolerance` changes its units deliberately.

   The absolute rule is deleted when the last caller has moved.
4. **B1 qualifies one element-wise rule**, |a_pq| ≤ tol·√|a_pp·a_qq|, on positive definite and indefinite inputs alike. It splits by definiteness only if a concrete case defeats it. The exact stopping guard is left to implementation, under three acceptance cases:
   - **A graded positive definite matrix.** A unit entry beside the block 10⁻²⁰·[[1, 0.5], [0.5, 1]] returns 0.5·10⁻²⁰ and 1.5·10⁻²⁰ to small relative error.
   - **An indefinite matrix with zero diagonals.** For [[0, 1], [1, 0]] the element-wise threshold is zero, so the nonzero pivot must rotate, and the solve returns −1 and +1.
   - **Unequally scaled entries in one pivot pair.** For [[1, 5·10⁻²¹], [5·10⁻²¹, 10⁻⁴⁰]] the small eigenvalue is about 7.5·10⁻⁴¹. A guard of ε·max(|a_pp|, |a_qq|) would skip the rotation and return 10⁻⁴⁰, about 33% high. The guard must not.

   When a pivot fails the accuracy test and floating-point arithmetic can make no further progress, the convergence report says so.
5. **The matrix Poincaré check reports one ratio per draw and stays opt-in.** Full-rank results are cached per structure; rank-deficient draws that share terms can differ. The question that would reopen the default is answered by timing representative BASS readouts with the check on and off.

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
  - the SIMD notes, including when to revisit FMA in another kernel. A profile showing substantial time in a multiply-add loop is a reason to compare fused and unfused arithmetic at caller sizes. Use `DenseProduct` as the worked example for matching scalar tails to the vector body, and link the G1 results for rotations when available;
  - the `FastVariant` option row.
- [ActiveSubspaces.md](../../docs/Maths/Regression/ActiveSubspaces.md): the unit-scale paragraph, which B1 retires; the matrix Poincaré section's contract, factor form and resolution.
- If item 1 is decided as removal: [Architecture.md](../../docs/Architecture.md) (`ExtraDefineConstants`), [Design.md](../../docs/Design.md) (the EigenFast versus Eigen entry), [README.md](../../README.md) (the reference build command).

## Reports

Each slice is reported here, newest first, with its commit, fixtures, commands, receipts, pins moved and doc changes owed.

### Supervised Sol sessions

| Slice | Worker chat | Report |
|---|---|---|
| B1 relative stopping rule and first caller migration | [B1 dense Jacobi and first caller migration](codex://threads/01a1d201-30f7-7941-aeb4-130bc14f7ae5) | [B1 stage one](spectral-b1-stage-one-report.md) |
| B1 stage two, migrating RidgeSpectrum, BassRandomIntercepts and MoMPCA in both builds | [Same B1 worker](codex://threads/01a1d201-30f7-7941-aeb4-130bc14f7ae5) | [B1 stage two](spectral-b1-stage-two-report.md) |
| F ordinary GP fit and leave-out factorization | [F GP factorization and leave-out costs](codex://threads/01a1d201-5ab1-74a2-81cd-0aa671ed7905) | [GP factorization](spectral-gp-factorization-report.md) |
| F anisotropic tensor-spline factorization and end-to-end REML comparison | [Same GP and factorization worker](codex://threads/01a1d201-5ab1-74a2-81cd-0aa671ed7905) | [Tensor-spline comparison](spectral-tensor-factorization-report.md) |

The [coordinating chat](codex://threads/01a1cc28-f4d3-7951-882a-65aeb6945f8c) reviews the reported commits and receipts. Each worker maintains its own report beside this brief.

At the stage-two handoff the coordinator checked the saved test counters. B1 stage one passed 66 focused regression tests in each build and 2,901 in the rebuilt default commit gate after five pin updates. The GP worker's isolated snapshot passed 152 in the commit gate and 3,280 in the full regression assembly. Stage two begins at `cb94b17`, which contains both code changes. Graph callers and the absolute default remain assigned to the later migration.

**Stage-two review.** MoMPCA landed in `f48613f` and the two regression callers in `672d1c2`. All six owned source and test files match the qualified hashes. The focused regression, geometry and LinAlg gates passed in both builds. The full default run passed 3,281 of 3,284 with the former held-group pins, and the rebuilt commit gate passed all 2,778 after the three pin updates. The hinge and tempered-hinge paths changed term sequences in 129 and 184 of 400 draws respectively. The report records the independent law checks and the limits of the trajectory comparison. [B1 stage-two qualification](spectral-b1-stage-two-report.md)

**Tensor factorization review.** Retain the compact banded LDLᵀ factor. Across two warmed processes and nine tensor shapes, baseline time divided by candidate time for the complete REML fit ranged from 0.84 to 1.55. The candidate allocated 1.28 to 5.74 times as many bytes. At fixed smoothing of 10²⁰ on both axes, feasibility changed on three shapes. All 225 search-box grid points and 990 distinct optimizer trial points agreed on feasibility. The full regression gate of 3,282 passes covers the production baseline. The probe supplies the candidate comparison. [Tensor-spline qualification](spectral-tensor-factorization-report.md)

### LOBPCG closeout: clustered eigenspaces, the capped constrained run and the Rayleigh–Ritz share (2026-11-12)

**Commit:** `e8d6b89` (tests only).

**Tests,** in `LobpcgDenseTests`, 12 of 12 passing in the default and `EIGEN_REFERENCE` builds:
- `Solve_ClusteredEigenvalues_AgreeOnTheClustersEigenspace` closes B3's clustered-eigenspace gap. It has six cases:
  - A = QΛQᵀ with Q random orthogonal places four eigenvalues at 10 above a spectrum on [0, 8], or three at 0 below one on [2, 10].
  - Each cluster is exact or spread by 10⁻¹².
  - The blocks are as wide as the cluster or narrower.

  The checks:
  - the returned vectors lie in the cluster's eigenspace to the Davis–Kahan bound ‖R‖_F/δ plus 10⁻¹²;
  - the values lie inside the cluster's interval widened by their residuals;
  - the reported residual is the independent one.

  Exact ties are where the production kernel's sort decides which Ritz vectors are kept.
- `Solve_CappedConstrainedPreconditioned_DescribesReturnedPairs` closes A1's gap. It runs the constrained, preconditioned fixture capped at three iterations. Each check holds to 10⁻¹²:
  - the reported residual is the operator applied to the returned pairs;
  - each value is its vector's Rayleigh quotient;
  - the vectors are orthogonal to the constraint.

  The converged carried-image test now shares its fixture.

**Measurements:**
- Setup: a `git archive` snapshot of `898a531` instrumented by `artifacts/spectral-closeout/patch_lobpcg.py`.
  - It times each Rayleigh–Ritz solve and switches its stopping rule.
  - It re-decomposes each projection under both rules, timed adjacent.
- Probe: `LobpcgRitzShareProbe.cs.txt`.
- Receipts: `ritz-share-{1,2}.txt`.
- Two rounds. Counts and iterations are identical across rounds, and the share ranges span both.

| Fixture | Rayleigh–Ritz share, absolute → relative | Iterations | Operator applications |
|---|---|---|---|
| Crescent, 500 nodes, k = 8, Jacobi preconditioner | 6.6–6.7% → 7.5–7.6% | 35 → 38 | 663 → 692 |
| Same, null space deflated | 6.6–7.5% → 7.3–7.4% | 33 → 32 | 641 → 632 |
| Möbius, 500 nodes, k = 8, Jacobi preconditioner | 5.9–7.0% → 7.2–7.3% | 42 → 41 | 756 → 748 |
| Same, null space deflated | 7.3% → 7.6–8.0% | 35 → 36 | 691 → 699 |
| Dense 256, k = 8 | 3.7–4.0% → 5.0–5.5% | 111 → 95 | 2,125 → 2,071 |
| Output operator G(Gᵀv), m = 2000, 200 factor columns, top 4 | 0.30–0.32% → 0.34–0.35% | 98 → 93 | 1,070 → 1,057 |
| Same, top 8 | 0.60–0.62% → 0.73–0.78% | 96 → 89 | 1,819 → 1,907 |
| Same, m = 8000, top 4 | 0.06% → 0.07% | 331 → 242 | 2,636 → 2,449 |

**Findings:**
- **The Rayleigh–Ritz solve is 4–8% of a solve.** That holds on the graphs and the dense control, at projected orders up to 24. On the output operators the operator dominates and the share is under 1%. That operator applies column by column, unblocked; a blocked apply would raise the share. So B1's change of default cannot move LOBPCG's cost much through this solve, and E's two arms can straddle it.
- **On identical projections the relative rule never stalls or fails to converge.** It takes the same number of sweeps or fewer, and the kept Ritz values agree to 1.5·10⁻¹⁵ of the projection's scale. The undeflated graph fixtures carry the null vector's Ritz value, about 10⁻¹⁵, on the diagonal, and it caused no extra sweeps.
- **Each relative solve still costs more.** It was slower in 31 of 32 adjacent pairs, typically by 10–25% and at most 36%. Its pivot test takes two square roots for every pivot it examines, and its tangent adds three divisions. At orders 12–24 the test is a large part of the work.

  For B1's final kernel or G2: a diagonal entry changes only when a rotation touches it, so its root could be carried rather than recomputed per pivot.
- **Rounding-level differences change LOBPCG's path.** Iteration counts move by −27% to +9% between rules whose Ritz values agree to 10⁻¹⁵. A single LOBPCG timing is one sample of a trajectory, so E times over several seeds.
- **Iteration counts depend on the spectral gap.** The output fixtures' random factors have no planted gap, which is why they take 89–331 iterations. E's fixtures state their gaps, since iteration counts scale with them.

**Not run:** B3's before/after count of operator applications. B3's receipts already record iterations (37 → 36, 42 → 42, 110 → 110), and B5 replaced the iteration they describe. The count would describe code that no longer runs and would change no decision. The table gives the current code's counts under each rule.

**Gates:** `LobpcgDenseTests` 12 of 12 in each build. The probe ran on the snapshot only; production is unchanged.

**Pins moved:** none.

**Doc changes owed:** [LinearAlgebra.md](../../docs/Maths/LinearAlgebra.md), LOBPCG:
- clustered spectra are compared by eigenspace;
- the Rayleigh–Ritz share;
- the relative rule's per-pivot cost at small orders.

### DenseEigenOptions and DenseEigenFastVariant removed (2026-11-12)

**Commit:** `75f76c1`. Decided by the owner after decision 2 left the enum with one value and the record carrying nothing.

**Consequences:**
- **Public API.**
  - Removed parameters:
    - the options parameter of every `DenseEigen` entry point;
    - the `fastVariant` parameter of `EigenFast`'s;
    - the options parameter of both `SpectralMath.BottomK` overloads;
    - `denseOptions` on `Spectral.ComputeBottomK` and `CoherentField.Build`.
  - Also removed: the reference build's compatibility check.
  - Already-compiled consumers must rebuild, even if their source omitted the optional argument, and callers that passed it explicitly must edit the call.
- **Saved configs.** Removing `AlignmentEigenCentrality.EigenOptions` changes the class's serialized shape. No tracked JSON, YAML or TOML config contains the field, and no in-repo loader reads this class. An external config that sets `EigenOptions` migrates by deleting that property and keeping `MaxSweeps` and `Tolerance`.
- **Benchmark output.** The timing harness drops `FastVariant` from its text summary and from newly generated JSON. `Materialization` stays, since it still distinguishes rectangular from flat column-major input. Historical receipts keep the variant that was actually measured.

B1 can introduce an explicitly named stopping rule when callers need to choose one.

**Gates:**
- `ThermoMapper.sln` builds, Release, with no errors.
- The focused dense and graph dispatch tests pass in both builds, 46 of 46 each. They are `EigenFastTests`, `EigenDecomposeSymmetricTests`, `DecompositionLogTests`, `LobpcgDenseTests`, `SpectralDenseDispatchTests`, `SpectralSolverPolicyTests`, `SpcDiagnosticsTests`, `IcaTests` and `PcaTests`. The two builds' `Maths.LinAlg.dll` differ in size, 110,080 against 109,568 bytes, so the reference define applied.

**Pins moved:** none.

**Doc changes owed:** [LinearAlgebra.md](../../docs/Maths/LinearAlgebra.md), the only page that names the removed types: `DenseEigenOptions`, the `FastVariant` row and the reference build's exception behavior.

### Decision 2: the FMA variant removed (2026-11-12)

**Commit:** `9635dca`.

**What changed:**
- `DenseEigenFastVariant.Fma` and EigenFast's fused branches are gone. The default arithmetic is unchanged, so outputs are the same bits.
- The timing harness drops its FMA comparison and spectral cases. `SpectralDenseDispatchTests` drops its FMA-against-default test.
- `DenseEigenFastVariant` keeps one value, and `DenseEigenOptions` now carries nothing. Deleting both would reach `Spectral.ComputeBottomK`'s and `CoherentField`'s `denseOptions` parameters, `SpectralMath.BottomK`, and the public `AlignmentEigenCentrality.EigenOptions`. That is a removal of its own, for the owner.

**Gate:** `VizCore.Tests` selection, 33 classes: 236 passed, 2 standing skips.

**Doc changes owed:** [LinearAlgebra.md](../../docs/Maths/LinearAlgebra.md), the `FastVariant` option row and the notes on EigenFast's segments.

### B5–B6: carried operator images and the matrix-function helper (2026-11-12)

**Commits:** `f253ec6` (B6), `9988a7f` (B5).

**B6, what changed:** `MatrixFunctions` in LinAlg provides:
- `Symmetric(eigen, f)`: f(A) = V·f(Λ)·Vᵀ, symmetrized, skipping zero weights and zero components, which leaves every sum unchanged;
- `FactorColumns` and `FactorRows`: V·f(Λ) over the eigenpairs a caller selects.

The five hand-written loops now call it:
- the whitening root and its inverse;
- the output metric's symmetric root (`PsdRootOf`) and row root (`RootOf`);
- OutputFactor's term factor;
- the analysis's derived inverse map.

Each caller keeps its own cutoff and refusal rules. The readout probe after B6 is identical to the probe after B2. The one change in value: the derived inverse of a hand-built analysis is now symmetrized, a change at rounding level. `MatrixFunctionsTests` checks the root, the inverse and the factors against a matrix of known spectrum.

**B6, gate:** run from a `git archive` snapshot of `f253ec6`, edit scope, 214 classes: 2,491 passed.

**B5, what changed:** after the first iteration, LOBPCG carries X and AX from the trial step. The operator and the k-dimensional Rayleigh–Ritz step no longer run again on them. A residual below the tolerance on carried images is confirmed on a fresh application before the solve stops, so the reported residual is the returned pairs'.

**B5, measurements:** reuse probe `artifacts/spectral-b/LobpcgReuseProbe.cs.txt`, receipts `reuse-{fresh,carried}-{1,2}.txt`. Two rounds each, 12 warm-up solves; the time column is the better minimum of the two rounds.

| Fixture | Applications, fresh → carried | Iterations | Minimum time, fresh → carried |
|---|---|---|---|
| Crescent, 500 nodes, k = 8, Jacobi preconditioner | 945 → 663 | 36 → 35 | 34.0 → 25.8 ms |
| Same, null space deflated | 907 → 641 | 34 → 33 | 32.8 → 24.8 ms |
| Möbius, 500 nodes, k = 8, Jacobi preconditioner | 1,076 → 756 | 42 → 42 | 36.7 → 29.7 ms |
| Same, null space deflated | 955 → 691 | 35 → 35 | 31.8 → 25.3 ms |
| Dense 256, k = 8 | 2,980 → 2,125 | 110 → 111 | 184 → 138 ms |

- **Effect:** applications fall 28–30% and minimum times 13–32%. Reported residuals match independently recomputed ones in every run.
- **One case is slower in iterations.** On a constrained, preconditioned dense case (n = 80, k = 5, the lowest eigenvector as the constraint), carried takes 98 iterations to fresh's 87, at about 15 operator applications per iteration against 20. Both stall at the same residual floor when asked for 10⁻¹⁰: 7.0·10⁻⁹ carried, 8.7·10⁻⁹ fresh. Both lean on the constraint alike: 3.8·10⁻⁹ and 4.6·10⁻⁹. The floor comes from the constraint being an approximate eigenvector, not from carrying.
- **Adopted:** the whole solve benefits and the convergence checks hold.

**B5, tests:**
- `Solve_RayleighRitzRunsOnTheDenseDispatch` now counts I projection solves for a solve converged at iteration I, and 4 for one capped at three iterations.
- `Solve_CarriedImages_ReportTheReturnedPairsResidual` checks a constrained, preconditioned solve at a tolerance of 10⁻⁷:
  - the reported residual equals the operator applied to its pairs, to 10⁻¹²;
  - the eigenvalues match the dense spectrum to 10⁻⁹;
  - orthogonality to the constraint holds to 10⁻⁸.

**B5, gate:** `VizCore.Tests` selection, 29 classes: 204 passed, 2 standing skips. `TDA.Mapper.Tests`: 14 passed.

**Pins moved:** none.

**Doc changes owed:** [LinearAlgebra.md](../../docs/Maths/LinearAlgebra.md):
- `MatrixFunctions`;
- LOBPCG's iteration with carried images and the confirmation on a fresh application before stopping.

### B2–B3: caller-owned dense solves and the production Rayleigh–Ritz step (2026-11-12)

**Commits:** `84f001e` (B2), `1adf27d` (B3).

**B2, what changed:**
- `EigenFast.DecomposeSymmetricInPlace` takes an empty eigenvector buffer for eigenvalues alone and skips the eigenvector rotations. The eigenvalues and the convergence report stay the same bits.
- The final sort runs the same exchanges on an index permutation and moves each column once along its cycles: O(n²) where it was up to O(n³). A repeated eigenvalue's vectors keep their order, so no output moved.
- The convergence report is read without a closure. A solve into caller-owned buffers allocates nothing once the pools are warm.
- `DenseEigen` gains `DecomposeSymmetric(flat, n, values, vectors)` and `DecomposeSymmetricInto(matrix or flat, values, columns)`. The second writes the eigenvectors once, as the columns of a `double[,]`.
- `ScaledEigen.Into` writes the same bits as `ScaledEigen.Of`.
- `SubspacePosterior` (twice), `RankPosterior` and `CoActivePosterior` (twice) read columns directly instead of transposing rows.
- `PenalizedSpline` and the Poincaré pencil's last solve read eigenvalues alone.

**B2, outputs:** the BASS readout probe after B2 is identical to the probe after C5, line for line. `EigenFastTests` checks:
- the caller-owned entries against the allocating overload, bit for bit;
- that eigenvalues alone match;
- the exchange order on diag(1, 1, 2);
- zero bytes allocated when warm.

**B3, what changed:**
- LOBPCG projects into a workspace sized for 3k and allocated once per solve, solves it through `DenseEigen`'s caller-owned entry, and reads the extremal Ritz vectors from the eigenvector columns. `BuildProjectedMatrix` and `BuildExtremalEigenvectorMatrix` are gone.
- `HierarchicalPSpline` decomposes its Gram through `DenseEigen` at the tolerance the reference used, 10⁻¹². Its conditional tests pass in the commit gate.
- The reference Jacobi's only production references are now `DenseEigen`'s `#if EIGEN_REFERENCE` branches.

**B3, work:** `LobpcgDenseTests.Solve_RayleighRitzRunsOnTheDenseDispatch` checks the projection solves:
- a solve that converges at iteration I runs 2I − 1;
- a solve capped at three iterations runs six;
- all run on the build's dense dispatch, none on the reference.

**B3, measurements:** probe `artifacts/spectral-b/LobpcgTimingProbe.cs.txt`, receipts `artifacts/spectral-b/lobpcg-*.txt`. Setup: Release, 18 cores, AVX2 and FMA, default tiering and PGO, below-normal priority on a shared machine.

| Fixture | Before, median / min | After, median / min | Allocation per solve | Iterations |
|---|---|---|---|---|
| Crescent graph, 500 nodes, k = 8 | 52 / 45 ms | 40 / 34 ms | 1,216 → 550 KB | 37 → 36 |
| Möbius graph, 500 nodes, k = 8 | 52 / 47 ms | 43 / 39 ms | 1,287 → 550 KB | 42 → 42 |
| Dense 256, k = 8 | 282 / 263 ms | 227 / 205 ms | 2,902 → 795 KB | 110 → 110 |

- **Warm-up:** these rows are the third round, with 12 warm-up solves on each side. With three warm-up solves, the post-change crescent medians read 106–115 ms while their minima held at 38–44 ms. That was a warm-up effect, gone at 12.
- **Load:** load moved timings by up to about 25% between rounds. So: no regression, the size of the gain unsettled.
- **Agreement:** reported and independently recomputed residuals match, all below 10⁻⁹, and the eigenvalues agree with the earlier solver's to about 10⁻¹⁵.
- **The kernels at the projected sizes**, reference against production:

  | Order | Reference | Production |
  |---|---|---|
  | 16 | 73–81 µs | 30–49 µs |
  | 24 | 265–331 µs | 97–112 µs |

  The production kernel allocates nothing, against 2.7–16 KB. At order 8 the readings disagree across rounds, 9.5–12 µs against 16–20 µs in three of them, with two rounds reversed by load.
- **Not measured:** the projection's share inside the whole solve, operator applications counted on their own (iteration counts stand in), and subspace agreement for clustered eigenvalues.

**Gates:**
- `VizCore.Tests`, 33 classes: 236 passed, 2 standing skips.
- `Maths.Regression.Tests` at commit scope, 284 classes with the pinned tier: 2,899 passed.
- `Maths.LinAlg.Tests` 22, `Maths.Geometry.Tests` 13, `TDA.Mapper.Tests` 14: all passed.
- `Maths.Oracle.Tests`: 135 of 136. `DiagramMetricsTransportParityTests(p: 2)` hit its R oracle's 120 s timeout under load. It reaches no eigensolver, and it passed when rerun alone.

**Pins moved:** none.

**Doc changes owed:** [LinearAlgebra.md](../../docs/Maths/LinearAlgebra.md):
- the caller-owned entries and eigenvalues alone;
- the sort;
- the diagram's LOBPCG Rayleigh–Ritz arrow, which now points to `DenseEigen`, not the reference;
- `DenseEigenOptions`' scope, now including LOBPCG's projections.

### C4–C5: the sampled Gram and the whitening inverse (2026-11-12)

**Commit:** `72c5838`.

**What changed:**
- **C4.** `DrawSubspace.FromFactor(factor, c, dimension)` takes C as the sampled route's running sums accumulated it and keeps it as the draw's C. It no longer forms FᵀF again. The thin SVD still reads the subspace, and fewer factor rows than the dimension still fall back to `FromC`.
- **C5.** This landed ahead of B6:
  - `Whitening.SymmetricRoots` returns Σ^{1/2} and Σ^{-1/2} from one eigensystem. The root's bits are unchanged, and `SymmetricRoot` reads it.
  - The readout's frame carries L⁻¹: the diagonal's reciprocals, or the inverse root.
  - `ActiveSubspaceAnalysis.MapInverse` (internal) keeps the inverse with the map it inverts. A copy with another `Map`, or an analysis built by hand, derives it from its own map, as the projection did before.
  - `BassProjection.InverseTranspose` is gone.
  - B6's helper will absorb `SymmetricRoots`' loop with the other four.

**Outputs:** compared against the earlier sources with the probe extended by numeric lines. Only the sampled readout moves, at rounding:
- each draw's C, by at most 1.0·10⁻¹⁴ relative;
- the mean C, by 5·10⁻¹⁶;
- the eigenvalues, by 1.5·10⁻¹⁵;
- the whitened projection, by 4.6·10⁻¹⁶.

Every other probe line is the same bits: the exact, identity-metric and basis readouts, adequacy, the matrix ratios, the sampled route's Monte Carlo error and output side, and the whitening map. The probe's outputs are `artifacts/spectral-c/probe-c45-*.txt`.

**Work:** two new tests in [ReadoutDecompositionTests](../../tests/regression/ReadoutDecompositionTests.cs).
- `FromFactor_KeepsTheAccumulatedGram` checks that:
  - the draw's C is the matrix passed in;
  - the only decomposition is one thin SVD;
  - the eigenvalues agree with C's own;
  - the fallback matches `FromC`;
  - a mis-shaped C is refused.
- `WhitenedProjection_ReadsTheReadoutsInverseRoot` checks that:
  - a whitened projection runs no symmetric decomposition (before: 1);
  - L·L⁻¹ = I to 10⁻¹³;
  - an analysis built by hand runs one, and projects alike to 10⁻¹²;
  - a copy with a new map derives its own inverse.

**Gate:** edit scope over the C4 and C5 paths, 214 classes: 2,491 passed, none failed. The one pinned test matching these names reads random effects' own `Whitening`, which this work does not touch.

**Pins moved:** none.

**Doc changes owed:** [ActiveSubspaces.md](../../docs/Maths/Regression/ActiveSubspaces.md), `Whitening.SymmetricRoots` and `DrawSubspace.FromFactor`'s signature.

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

**Not done:** whole-readout timings for [Decisions](#decisions) item 5.

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
