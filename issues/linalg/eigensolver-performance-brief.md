# LinAlg eigensolvers — residuals, Jacobi stopping and the LOBPCG projected solve

Queued 2026-11-09 for the next available LinAlg owner. The first two fixes can start independently of E21, E34 and conformal slice 8. Schedule comparative timings when the regression qualification runs leave the machine quiet. The [regression handoff](../regression/coordination-handoff-2026-11-09.md) tracks the separate E8 decision on a relative tolerance for small-scale matrices.

## 1. Make the two failure cases reproducible and fix them

- **LOBPCG residual at the iteration cap.** The loop measures the residual on the current block, then can return the updated block with that earlier measurement. Reproduce with the symmetric tridiagonal matrix with diagonal `[4, 3, 2]` and adjacent entries 1, largest pair, seed 173, tolerance `1e-12` and one iteration. In the review, the result reported 0.76200174 while its returned pair gave 0.34933228. Make `ResidualNorm` and `Converged` describe the returned eigenpairs. Verify by independently applying the operator to those pairs on capped and converged runs. [Iteration and return](../../src/maths/linalg/Lobpcg.cs)
- **EigenFast sweep with no rotations.** Use a 32 × 32 matrix with diagonal 1 through 32 and every off-diagonal `5e-13`, at the existing `1e-12` tolerance. Every pivot is skipped while their aggregate norm exceeds the sweep exit threshold. Caps of 1 and 256 produced identical eigenpairs in the review. Make a sweep that performs no rotations terminate. Gate this with a deterministic check of work performed, plus residual and orthogonality checks on the witness and on inputs that require rotations. Keep the existing pivot tolerance for this fix. [Pivot and sweep logic](../../src/maths/linalg/EigenFast.cs)

Land each fix as its own concern. The [dense LOBPCG tests](../../tests/VizCore.Tests/LobpcgDenseTests.cs), [EigenFast tests](../../tests/VizCore.Tests/EigenFastTests.cs) and [graph solver parity checks](../../tests/VizCore.Tests/EigenFastTimingHarness.cs) are the starting gates.

## 2. Give the projected solve a production implementation

LOBPCG applies the original operator to blocks. Its Rayleigh–Ritz step currently calls the scalar reference `Eigen.DecomposeSymmetric` on a projection of order `k`, and again on a trial projection of order up to `3k`. Keep that reference available as a numerical oracle. The production path should solve the projected problem with reusable workspace while retaining matrix-free application of the original operator. [Ritz step](../../src/maths/linalg/Lobpcg.cs) · [Dense solver interfaces](../../docs/Maths/LinearAlgebra.md)

Measure after the stopping fix. The review's in-memory probes gave the reference solver 22% and 31% of whole-solve time on two synthetic fixtures. Calling the allocating `EigenFast` overload directly made both runs slightly slower. Reproduce the baseline and compare alternatives at the projected sizes actually encountered. Include the sparse graph fixtures at `k = 8` and a dense control. Record time and allocations for the projection and whole solve, operator applications, iterations, convergence and independently recomputed residuals. Check eigenvalues, orthogonality and subspace agreement for clustered eigenvalues. Record runtime, build, tiering/PGO, warm-up and hardware SIMD support. The [timing harness](../../tests/VizCore.Tests/EigenFastTimingHarness.cs) provides the existing fixtures.

## 3. Measure reuse of the operator images

The trial update already rotates both `X` and `AX`. The next iteration reapplies the operator to `X`. Compare carrying `AX` across that boundary with fresh application as a separate experiment, using the same projected solver on both sides. Count operator applications and check residual drift independently, including for constrained and preconditioned runs. Adopt reuse if the complete solve benefits and the convergence checks hold.

Report each slice here with its commit, fixture, command and receipt. Decisions about the global dense tolerance, graph solver selection thresholds, and E8's full-spectrum and sampled-SVD readouts stay with their respective owners.
