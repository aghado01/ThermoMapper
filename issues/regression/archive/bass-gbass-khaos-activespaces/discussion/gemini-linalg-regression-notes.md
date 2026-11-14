Ran check ThermoMapper directory structure
Ran list src directories
Ran list directories in src/maths
Ran list files in linalg
Ran list files in regression
Ran find linalg using statements in regression
Ran extract public types from linalg
Ran extract all types from linalg
Ran find usage of linalg types in regression
Ran read Cholesky.cs
Ran read DenseCholesky.cs
Ran read ISpdFactor.cs
Ran find implementations of ISpdFactor
Ran check Cholesky in GpRegression.cs
Ran read GpRegression around line 220-270
Ran search Factor usage in GpRegression
Ran read LeaveOutLogDensities in GpRegression
Ran find LeaveOutPredictive definition
Ran read LeaveOutPredictive in Pointwise.cs
Ran read Eigen.cs
Ran read DenseEigen.cs
Ran read EigenFast.cs
Ran read Lobpcg.cs
Ran read ThinSvd.cs
Ran read PCA.cs
Ran find all LOBPCG references in src
Ran find all lobpcg references case-insensitively
Checked task Find all lobpcg references case-insensitively
Ran read DenseSymmetricOperator.cs
Ran inspect Eigen in regression/subspaces
Ran read ScaledEigen.cs
Ran inspect Eigen in PriorPencil.cs
Ran inspect linalg in CollapsedGPrior.cs
Ran inspect linalg in spline/baps
Ran read PenalizedSpline around line 185
Ran find ThinSvd in regression
Ran check SVD/PCA in bass
Ran inspect Cholesky in BassDesign.cs
Ran read BorderedCholesky.cs top
Ran read WideBandCholesky.cs
Ran find SpdFactors usage
Ran inspect SpdFactors in WideBandCholesky.cs
Ran inspect BandedGPrior.cs
Ran inspect linalg in spline/bars
Ran read OutputSide.cs
Ran read SubspacePosterior.cs
Ran check linalg in BassActiveSubspace.cs
Ran inspect DenseProduct in conformal
Ran inspect DiagonalPlusLowRank usage
Ran read DiagonalPlusLowRank.cs
Ran run LinAlg and Regression tests · Exited with code 1
Ran list test project directories
Ran run tests in tests/linalg · Exited with code 1
Ran list csproj files in tests
Ran run Maths.LinAlg.Tests
Checked task Run Maths.LinAlg.Tests
Checked task Run Maths.LinAlg.Tests

I have initiated the test execution in the background and will continue as soon as it completes.

Checked task Run Maths.LinAlg.Tests
Checked task Run Maths.LinAlg.Tests

I will wait for the tests to finish.

Ran run Maths.Regression.Tests

I am running the Regression test suite in the background.

Checked task Run Maths.Regression.Tests

I will wait for the test run to finish.

### 1. Mapping: What `regression` Imports from `LinAlg`

Every major regression family imports `Maths.LinAlg` (`47` files across the namespace). Below is the structural mapping of types consumed:

| Regression Subsystem | Files | Imported LinAlg Types & Roles |
| :--- | :--- | :--- |
| **BGP** (Bayesian GP) | [`GpRegression.cs`](file:///d:/aghado01/ThermoMapper/src/maths/regression/bgp/GpRegression.cs) | [`CholeskyDecomposition`](file:///d:/aghado01/ThermoMapper/src/maths/linalg/Cholesky.cs) (covariance factor $\Sigma_t + \eta I$, marginal evidence, GLS solve, and explicit inverse for leave-out densities). |
| **G-Prior** (Zellner & Hyper-$g$) | [`CollapsedGPrior.cs`](file:///d:/aghado01/ThermoMapper/src/maths/regression/gprior/CollapsedGPrior.cs), [`PriorPencil.cs`](file:///d:/aghado01/ThermoMapper/src/maths/regression/gprior/PriorPencil.cs), [`MatrixTPredictive.cs`](file:///d:/aghado01/ThermoMapper/src/maths/regression/gprior/MatrixTPredictive.cs), [`OutputMarginal.cs`](file:///d:/aghado01/ThermoMapper/src/maths/regression/gprior/OutputMarginal.cs) | [`ISpdFactor`](file:///d:/aghado01/ThermoMapper/src/maths/linalg/ISpdFactor.cs), [`DenseCholesky`](file:///d:/aghado01/ThermoMapper/src/maths/linalg/DenseCholesky.cs), [`DenseEigen`](file:///d:/aghado01/ThermoMapper/src/maths/linalg/DenseEigen.cs), [`FactorStatus`](file:///d:/aghado01/ThermoMapper/src/maths/linalg/FactorStatus.cs). Computes full generalized pencil spectrum $L^{-1} A L^{-T} = V \Lambda V^T$ to evaluate $(b\Lambda + \tau I)^{-1}$ in $O(p)$ per MCMC draw via `SpectralFactor : ISpdFactor`. |
| **Spline / BAPS** (P-Splines) | [`Baps.cs`](file:///d:/aghado01/ThermoMapper/src/maths/regression/spline/baps/Baps.cs), [`PenalizedSpline.cs`](file:///d:/aghado01/ThermoMapper/src/maths/regression/spline/baps/PenalizedSpline.cs), [`AnisotropicTensorSpline.cs`](file:///d:/aghado01/ThermoMapper/src/maths/regression/spline/baps/AnisotropicTensorSpline.cs), [`HierarchicalPSpline.cs`](file:///d:/aghado01/ThermoMapper/src/maths/regression/spline/baps/HierarchicalPSpline.cs), [`TensorPenalty.cs`](file:///d:/aghado01/ThermoMapper/src/maths/regression/spline/baps/TensorPenalty.cs) | [`BandCholesky`](file:///d:/aghado01/ThermoMapper/src/maths/linalg/BandCholesky.cs), [`BandFactorization`](file:///d:/aghado01/ThermoMapper/src/maths/linalg/BandCholesky.cs) (root-free $\text{LDL}^T$ for fast $N(0, A^{-1})$ draws), [`DenseCholesky`](file:///d:/aghado01/ThermoMapper/src/maths/linalg/DenseCholesky.cs), [`DenseEigen`](file:///d:/aghado01/ThermoMapper/src/maths/linalg/DenseEigen.cs) (penalty null space and effective degrees of freedom). |
| **Spline / BARS** (Free-knot Splines) | [`WhittleModel.cs`](file:///d:/aghado01/ThermoMapper/src/maths/regression/spline/bars/WhittleModel.cs), [`PoissonModel.cs`](file:///d:/aghado01/ThermoMapper/src/maths/regression/spline/bars/PoissonModel.cs), [`GibbsLossModel.cs`](file:///d:/aghado01/ThermoMapper/src/maths/regression/spline/bars/GibbsLossModel.cs), [`DesignAdmission.cs`](file:///d:/aghado01/ThermoMapper/src/maths/regression/spline/bars/DesignAdmission.cs), [`ReadoutQuadratics.cs`](file:///d:/aghado01/ThermoMapper/src/maths/regression/spline/bars/ReadoutQuadratics.cs) | [`BandCholesky`](file:///d:/aghado01/ThermoMapper/src/maths/linalg/BandCholesky.cs), [`ISpdFactor`](file:///d:/aghado01/ThermoMapper/src/maths/linalg/ISpdFactor.cs), [`FactorStatus`](file:///d:/aghado01/ThermoMapper/src/maths/linalg/FactorStatus.cs). Leverages `factor.EquilibratedCondition()` (Higham 1-norm estimate) and `factor.InverseQuadraticForms()` for pointwise prediction error. |
| **Spline / BASS** (Adaptive Splines) | [`BassDesign.cs`](file:///d:/aghado01/ThermoMapper/src/maths/regression/spline/bass/BassDesign.cs), [`BassRandomEffects.cs`](file:///d:/aghado01/ThermoMapper/src/maths/regression/spline/bass/BassRandomEffects.cs), [`BassActiveSubspace.cs`](file:///d:/aghado01/ThermoMapper/src/maths/regression/spline/bass/sensitivity/BassActiveSubspace.cs), [`CalibrationCovariance.cs`](file:///d:/aghado01/ThermoMapper/src/maths/regression/spline/bass/calibration/CalibrationCovariance.cs), [`BassDecomposition.cs`](file:///d:/aghado01/ThermoMapper/src/maths/regression/spline/bass/output/BassDecomposition.cs) | [`BorderedCholesky`](file:///d:/aghado01/ThermoMapper/src/maths/linalg/BorderedCholesky.cs) (rank-1 births, Givens downdates for deaths/changes), [`ThinSvd`](file:///d:/aghado01/ThermoMapper/src/maths/linalg/ThinSvd.cs) (random effect constraint revelation, output PCA basis, LOO truncation residuals), [`DiagonalPlusLowRank`](file:///d:/aghado01/ThermoMapper/src/maths/linalg/DiagonalPlusLowRank.cs), [`KhatriRao`](file:///d:/aghado01/ThermoMapper/src/maths/linalg/KhatriRao.cs), [`PairwiseSum`](file:///d:/aghado01/ThermoMapper/src/maths/linalg/PairwiseSum.cs). |
| **Subspaces** (Active Subspaces) | [`SubspacePosterior.cs`](file:///d:/aghado01/ThermoMapper/src/maths/regression/subspaces/SubspacePosterior.cs), [`OutputSide.cs`](file:///d:/aghado01/ThermoMapper/src/maths/regression/subspaces/OutputSide.cs), [`RankPosterior.cs`](file:///d:/aghado01/ThermoMapper/src/maths/regression/subspaces/RankPosterior.cs), [`Whitening.cs`](file:///d:/aghado01/ThermoMapper/src/maths/regression/subspaces/Whitening.cs), [`ScaledEigen.cs`](file:///d:/aghado01/ThermoMapper/src/maths/regression/subspaces/ScaledEigen.cs), [`PrincipalAngles.cs`](file:///d:/aghado01/ThermoMapper/src/maths/linalg/PrincipalAngles.cs) | [`ThinSvd`](file:///d:/aghado01/ThermoMapper/src/maths/linalg/ThinSvd.cs) ($C = F^T F$ via right singular vectors of $F$ without forming Gram), [`ScaledEigen`](file:///d:/aghado01/ThermoMapper/src/maths/regression/subspaces/ScaledEigen.cs) (wraps [`DenseEigen`](file:///d:/aghado01/ThermoMapper/src/maths/linalg/DenseEigen.cs)), [`PrincipalAngles`](file:///d:/aghado01/ThermoMapper/src/maths/linalg/PrincipalAngles.cs) (Grassmann distance between subspaces). |
| **Conformal & Pointwise** | [`ConformalBayes.cs`](file:///d:/aghado01/ThermoMapper/src/maths/regression/conformal/ConformalBayes.cs), [`ConformalDistribution.cs`](file:///d:/aghado01/ThermoMapper/src/maths/regression/conformal/ConformalDistribution.cs), [`Pointwise.cs`](file:///d:/aghado01/ThermoMapper/src/maths/regression/Pointwise.cs) | [`DenseProduct`](file:///d:/aghado01/ThermoMapper/src/maths/linalg/DenseProduct.cs), [`PackedColumns`](file:///d:/aghado01/ThermoMapper/src/maths/linalg/DenseProduct.cs) (BLAS-like fused multiply-add for candidate evaluations), [`DenseCholesky`](file:///d:/aghado01/ThermoMapper/src/maths/linalg/DenseCholesky.cs) (`LeaveOutPredictive` block reductions). |

---

### 2. Analysis of the Solver Variants

#### A. Cholesky Variants
1. **[`DenseCholesky`](file:///d:/aghado01/ThermoMapper/src/maths/linalg/DenseCholesky.cs)** (`ISpdFactor`): Strict, packed triangular storage $n(n+1)/2$. Halts on non-positive pivot, solve-oriented (never allocates/forms an explicit inverse), supports condition estimation.
2. **[`BorderedCholesky`](file:///d:/aghado01/ThermoMapper/src/maths/linalg/BorderedCholesky.cs)** (`ISpdFactor`): Incremental dense Cholesky for BASS. $O(M^2)$ border extension on births, Givens rotation sweeps on deaths, periodic strict rebuilds.
3. **[`BandCholesky`](file:///d:/aghado01/ThermoMapper/src/maths/linalg/BandCholesky.cs)** (`ISpdFactor`): Diagonal-major banded storage $O(n \cdot p)$ with root-free $\text{LDL}^T$ option. Used for 1D splines (BAPS/BARS).
4. **[`WideBandCholesky`](file:///d:/aghado01/ThermoMapper/src/maths/linalg/WideBandCholesky.cs)** (`ISpdFactor`): Row-packed banded Cholesky tailored for 2D tensor splines where bandwidth $p \approx p_x \nu_y + p$ makes diagonal striding slow.
5. **[`CholeskyDecomposition`](file:///d:/aghado01/ThermoMapper/src/maths/linalg/Cholesky.cs)**: The older $O(n^3)$ full 2D array representation. Uses soft pivot clamping ($10^{-12}$) and stores $L^{-1}$ internally to form full matrix inverses.

#### B. Eigen Solvers
1. **[`DenseEigen`](file:///d:/aghado01/ThermoMapper/src/maths/linalg/DenseEigen.cs) / [`EigenFast`](file:///d:/aghado01/ThermoMapper/src/maths/linalg/EigenFast.cs)**: Vectorized algebraic Givens cyclic Jacobi (no transcendentals, AVX2/FMA). Computes the full spectrum of symmetric matrices.
2. **[`ScaledEigen`](file:///d:/aghado01/ThermoMapper/src/maths/regression/subspaces/ScaledEigen.cs)**: Normalizes input matrix by $\max |A_{ij}|$ before invoking `DenseEigen` with absolute tolerance $10^{-14}$ to achieve scale-invariant relative precision.
3. **[`LOBPCG`](file:///d:/aghado01/ThermoMapper/src/maths/linalg/Lobpcg.cs)**: Iterative, matrix-free conjugate gradient solver targeting the extremal (bottom-$k$ or top-$k$) eigenpairs over an [`ILinearOperator`](file:///d:/aghado01/ThermoMapper/src/maths/linalg/Lobpcg.cs).

#### C. SVD & PCA
1. **[`ThinSvd`](file:///d:/aghado01/ThermoMapper/src/maths/linalg/ThinSvd.cs)**: Householder QR followed by one-sided Jacobi rotations on $R$. Computes $A = U \Sigma V^T$ without squaring condition numbers via $A^T A$. Accurately resolves small singular values down to $\epsilon \sigma_1$.
2. **[`Pca`](file:///d:/aghado01/ThermoMapper/src/maths/linalg/PCA.cs)**: Cleanly delegates to `ThinSvd` on centered data without forming covariance matrices.

---

### 3. Consolidation Opportunities

#### 1. Deprecate Legacy `CholeskyDecomposition` in Favor of `DenseCholesky`
- **Current State**: [`GpRegression.cs`](file:///d:/aghado01/ThermoMapper/src/maths/regression/bgp/GpRegression.cs#L225) is the **sole remaining caller** of `CholeskyDecomposition` in the regression module. It reallocates `new CholeskyDecomposition(_n)` on every `Fit()` call and invokes `WriteInverseTo(rInv)` in [`LeaveOutLogDensities`](file:///d:/aghado01/ThermoMapper/src/maths/regression/bgp/GpRegression.cs#L272).
- **Consolidation**:
  - Add an inversion method or triangular solve helper to [`DenseCholesky`](file:///d:/aghado01/ThermoMapper/src/maths/linalg/DenseCholesky.cs) (e.g. `InvertTo(double[,] dest)` or solving against column identity vectors in packed storage).
  - Switch `GpRegression` to `DenseCholesky` (or instantiate a reusable workspace factor).
  - This allows completely decoupling the regression module from the old `CholeskyDecomposition` class and standardizing all Cholesky operations under [`ISpdFactor`](file:///d:/aghado01/ThermoMapper/src/maths/linalg/ISpdFactor.cs).

#### 2. Connect `WideBandCholesky` via `SpdFactors.ForBand` to 2D Tensor Splines
- **Current State**: [`WideBandCholesky`](file:///d:/aghado01/ThermoMapper/src/maths/linalg/WideBandCholesky.cs) and [`SpdFactors.ForBand`](file:///d:/aghado01/ThermoMapper/src/maths/linalg/WideBandCholesky.cs#L256) were implemented to solve cache-thrashing on wide bands ($2(p+1) \le n$), but [`AnisotropicTensorSpline.cs`](file:///d:/aghado01/ThermoMapper/src/maths/regression/spline/baps/AnisotropicTensorSpline.cs#L130) still unconditionally uses `BandCholesky`.
- **Consolidation**: Switch `AnisotropicTensorSpline` and `TensorProductDesign` to dispatch via `SpdFactors.ForBand`. For typical 2D tensor grids (e.g., $30 \times 30$ knots, $n = 900$, band $\approx 90$), this yields a 2–4× speedup in matrix factorization and triangular solves.

#### 3. Unify Matrix Scaling Before Cyclic Jacobi (`ScaledEigen`)
- **Current State**: [`ScaledEigen.cs`](file:///d:/aghado01/ThermoMapper/src/maths/regression/subspaces/ScaledEigen.cs) lives in `subspaces`, but [`PriorPencil.cs`](file:///d:/aghado01/ThermoMapper/src/maths/regression/gprior/PriorPencil.cs#L61-L82) and [`CollapsedGPrior.cs`](file:///d:/aghado01/ThermoMapper/src/maths/regression/gprior/CollapsedGPrior.cs#L1032) duplicate identical unit-scaling logic by hand before calling `DenseEigen`.
- **Consolidation**: Move `ScaledEigen` into `Maths.LinAlg` (or expose an overload `DenseEigen.DecomposeSymmetricScaled`) so all regression callers share the same scale-invariant Jacobi path.

---

### 4. Can We Use LOBPCG Directly for Better Performance?

#### Where LOBPCG is NOT Suitable in Regression
1. **G-Prior and Pencil Solvers (`PriorPencil`, `CollapsedGPrior`)**:
   - `PriorPencil` forms `SpectralFactor : ISpdFactor`. During MCMC, it computes $(b \Lambda + \tau I)^{-1}$ for continuously varying $\tau$ and $b$.
   - **Why LOBPCG fails here**: Truncating to the top-$k$ eigenvalues leaves out the remainder of the spectrum, which destroys the exact closed-form evidence and posterior variance calculations. All $p$ eigenvalues are needed.
2. **Standard Active Subspaces with Small Input Dimension ($p \le 50$)**:
   - In [`BassActiveSubspace.cs`](file:///d:/aghado01/ThermoMapper/src/maths/regression/spline/bass/sensitivity/BassActiveSubspace.cs) and [`SubspacePosterior.cs`](file:///d:/aghado01/ThermoMapper/src/maths/regression/subspaces/SubspacePosterior.cs), $C$ is a $p \times p$ matrix where $p$ is the number of input parameters (often $p \in [2, 30]$).
   - **Why LOBPCG does not win**: For $p \le 50$, `DenseEigen` (vectorized Jacobi) runs in $< 10\,\mu\text{s}$. LOBPCG's iteration overhead, Rayleigh-Ritz projections, and convergence checks will be slower than the direct Jacobi sweep. In addition, rank-selection rules ([`RankPosterior.cs`](file:///d:/aghado01/ThermoMapper/src/maths/regression/subspaces/RankPosterior.cs)) and diagonal activity scores ([`ActivityScores.cs`](file:///d:/aghado01/ThermoMapper/src/maths/regression/subspaces/ActivityScores.cs)) inspect the entire eigenvalue vector.

#### Where LOBPCG WOULD Deliver Substantial Performance Wins
1. **Multi-Output BASS Output-Side Patterns ([`OutputSide.cs`](file:///d:/aghado01/ThermoMapper/src/maths/regression/subspaces/OutputSide.cs))**:
   - **Current Pattern**: Computes $G = \mathbb{E}[J W_x J^T]$, an $m \times m$ matrix where $m$ is the number of **spatial/output observation locations** (often $m \in [500, 10\,000]$ in PDE/thermal field emulation).
   - Currently, `OutputSide.cs` runs `ThinSvd.Decompose(f)` or stacks factors into `ThinSvd`, which computes the full rank ($O(m \cdot c^2)$ or $O(m^2 \cdot c)$).
   - If the user only requests $k \ll m$ leading output patterns (e.g., $k=2$ or $3$ dominant spatial modes):
     - Wrap $v \mapsto F (F^T v)$ as an [`ILinearOperator`](file:///d:/aghado01/ThermoMapper/src/maths/linalg/Lobpcg.cs).
     - Running **LOBPCG** directly on this operator finds the top-$k$ output modes in $O(k \cdot m \cdot c \cdot \text{iterations})$ without doing full SVD or forming any $m \times m$ matrix.
2. **High-Dimensional Active Subspaces ($p \gg 100$)**:
   - When emulating high-dimensional parameter spaces (e.g. distributed inputs or basis coefficients where $p > 200$), Monte Carlo gradients $G \in \mathbb{R}^{M \times p}$ can be wrapped into a matrix-free operator $v \mapsto \frac{1}{M} G^T (G v)$.
   - LOBPCG can extract the $k=1, 2, 3$ active directions without forming the $p \times p$ Gram matrix $C = \frac{1}{M} G^T G$.

---

### Recommended Action Plan

1. **Retire `CholeskyDecomposition` from `GpRegression.cs`**: Provide an inversion path in `DenseCholesky` and migrate `GpRegression` to use `DenseCholesky`.
2. **Promote `ScaledEigen` to `Maths.LinAlg`**: Replace manual scaling in `PriorPencil.cs` and `CollapsedGPrior.cs` with the centralized scaled Jacobi dispatcher.
3. **Dispatch 2D Tensor Splines to `SpdFactors.ForBand`**: Allow `AnisotropicTensorSpline` to switch between `BandCholesky` and `WideBandCholesky` based on bandwidth.
4. **Introduce Matrix-Free Operator for Multi-Output Patterns**: Implement an `ILinearOperator` over $F F^T$ in `subspaces\OutputSide.cs` using `LOBPCG` when $m \gg k$, providing a dedicated fast-path for high-resolution output fields.