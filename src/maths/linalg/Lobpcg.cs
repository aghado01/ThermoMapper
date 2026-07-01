// src/maths/linalg/Lobpcg.cs
#nullable enable
using System;
using System.Collections.Generic;

namespace Maths.LinAlg;

/// <summary>
/// Matrix-free symmetric linear operator <c>A</c>, applied to a block of column
/// vectors at once. Buffers are flat <b>column-major</b>: a block of
/// <c>columns</c> vectors occupies <see cref="Dimension"/> ×
/// <c>columns</c> doubles, with column <c>j</c> at offset <c>j * Dimension</c>.
/// </summary>
/// <remarks>
/// This is the dispatch seam that keeps <see cref="LOBPCG"/> a pure numerical
/// primitive: the solver never sees a graph, a Laplacian type, or a sparse
/// format. Concrete operators (graph Laplacian, dense matrix, Gram matrix of a
/// data block, …) live next to the structures they wrap and implement this
/// interface.
/// </remarks>
public interface ILinearOperator
{
    /// <summary>Row/column count <c>n</c> of the (square, symmetric) operator.</summary>
    int Dimension { get; }

    /// <summary>
    /// Computes <c>AX</c> for the column-major block <paramref name="block"/>
    /// (length <c>Dimension * columns</c>) into <paramref name="result"/>
    /// (same length). Must not alias <paramref name="block"/> with
    /// <paramref name="result"/>.
    /// </summary>
    void Apply(ReadOnlySpan<double> block, Span<double> result, int columns);
}

/// <summary>
/// Optional preconditioner — an approximation of <c>A⁻¹</c> applied to the
/// residual block to accelerate convergence. A <see langword="null"/>
/// preconditioner is treated as the identity. Same column-major block layout as
/// <see cref="ILinearOperator"/>.
/// </summary>
public interface IPreconditioner
{
    /// <summary>
    /// Applies the preconditioner to <paramref name="residualBlock"/>
    /// (length <c>n * columns</c>) into <paramref name="resultBlock"/>.
    /// </summary>
    void Apply(ReadOnlySpan<double> residualBlock, Span<double> resultBlock, int columns);
}

/// <summary>
/// Locally Optimal Block Preconditioned Conjugate Gradient eigensolver for the
/// extremal (smallest- or largest-magnitude) eigenpairs of a symmetric operator.
/// </summary>
/// <remarks>
/// <para>
/// Pure numerical primitive: it consumes an <see cref="ILinearOperator"/>, an
/// optional <see cref="IPreconditioner"/>, and an optional <em>constraint
/// block</em> <c>Y</c> the solution is kept orthogonal to. It has no knowledge of
/// graphs or Laplacians — callers wire those in via the operator interface
/// (see <c>Graphs.Spectral.GraphLaplacianOperator</c>).
/// </para>
/// <para>
/// <b>Constraint handling (deflation).</b> When known eigenvectors are supplied
/// as <see cref="Options.Constraints"/> (e.g. the constant null vector of a
/// combinatorial Laplacian), the iterate is deflated against <c>Y</c> on
/// <em>every</em> iteration, not merely at startup. A non-constant
/// preconditioner such as Jacobi (<c>D⁻¹</c>) reintroduces a component along
/// <c>Y</c> into the preconditioned residual even when the raw residual is
/// orthogonal to <c>Y</c>; re-projecting each iteration is what keeps the solver
/// from collapsing back onto the deflated mode.
/// </para>
/// </remarks>
public static class LOBPCG
{
    /// <summary>Solver inputs beyond the operator and the count <c>k</c>.</summary>
    public sealed record Options
    {
        /// <summary>Maximum outer iterations.</summary>
        public int MaxIterations { get; init; } = 300;

        /// <summary>Convergence threshold on the max per-column residual norm.</summary>
        public double Tolerance { get; init; } = 1e-9;

        /// <summary>
        /// When <see langword="true"/>, target the <c>k</c> <em>largest</em>
        /// eigenvalues; otherwise the <c>k</c> smallest.
        /// </summary>
        public bool WantLargest { get; init; } = false;

        /// <summary>
        /// Optional constraint block <c>Y</c> (flat column-major, length
        /// <c>Dimension * ConstraintColumns</c>) the solution is deflated against
        /// every iteration. Need not be orthonormal — the solver re-orthonormalizes
        /// it internally. <see langword="null"/> ⇒ no deflation.
        /// </summary>
        public double[]? Constraints { get; init; }

        /// <summary>Number of columns in <see cref="Constraints"/>.</summary>
        public int ConstraintColumns { get; init; }

        /// <summary>Optional preconditioner; <see langword="null"/> ⇒ identity.</summary>
        public IPreconditioner? Preconditioner { get; init; }

        /// <summary>
        /// Seed for the random orthonormal initial block. Fixed by default so a
        /// given operator yields reproducible eigenpairs run-to-run.
        /// </summary>
        public int Seed { get; init; } = 12345;
    }

    /// <summary>Computed eigenpairs plus convergence telemetry.</summary>
    public sealed record Result(
        IReadOnlyList<EigenPair> Eigenpairs,
        int Iterations,
        double ResidualNorm,
        bool Converged);

    /// <summary>
    /// Computes the <c>k</c> extremal eigenpairs of <paramref name="op"/>.
    /// </summary>
    public static Result Solve(ILinearOperator op, int k, Options? options = null)
    {
        if (op is null) throw new ArgumentNullException(nameof(op));
        options ??= new Options();

        int n = op.Dimension;
        if (n == 0 || k <= 0)
            return new Result(Array.Empty<EigenPair>(), 0, 0.0, true);

        k = Math.Min(k, n - 1);
        if (k <= 0)
            return new Result(Array.Empty<EigenPair>(), 0, 0.0, true);

        bool wantLargest = options.WantLargest;
        double tol = options.Tolerance;
        int maxIter = options.MaxIterations;
        IPreconditioner? preconditioner = options.Preconditioner;

        double[] x = MatrixOps.RandomOrthonormal(n, k, options.Seed);
        double[] xWork = new double[n * k];
        double[] ax = new double[n * k];
        double[] axWork = new double[n * k];
        double[] residual = new double[n * k];
        double[] preconditioned = new double[n * k];
        double[] conjugate = new double[n * k];
        double[] trial = new double[n * (3 * k)];
        double[] trialBasis = new double[n * (3 * k)];
        double[] appliedTrialBasis = new double[n * (3 * k)];
        double[] theta = new double[k];

        // Constraint block Y: re-orthonormalized once, then re-applied every
        // iteration. This is the deflation that keeps the iterate off known modes.
        double[]? constraints = null;
        int constraintCols = 0;
        if (options.Constraints is { } rawConstraints && options.ConstraintColumns > 0)
        {
            int requested = options.ConstraintColumns;
            var normalized = new double[n * requested];
            constraintCols = CompactOrthonormalize(
                rawConstraints.AsSpan(0, n * requested), n, requested, normalized);
            if (constraintCols > 0)
            {
                constraints = normalized;
                ProjectBlockAgainst(x, n, k, constraints, constraintCols);
                MatrixOps.Orthonormalize(x, n, k);
            }
        }

        double lastResidual = double.MaxValue;
        bool converged = false;
        int conjugateCols = 0;
        int completedIterations = 0;

        for (int iter = 0; iter < maxIter; iter++)
        {
            completedIterations = iter + 1;

            op.Apply(x, ax, k);
            RotateToExtremalRitzBasis(x, ax, n, k, k, theta, xWork, axWork, wantLargest);
            Swap(ref x, ref xWork);
            Swap(ref ax, ref axWork);

            ComputeResidualBlock(ax, x, theta, residual, n, k);
            lastResidual = ComputeMaxColumnNorm(residual, n, k);
            if (lastResidual < tol)
            {
                converged = true;
                break;
            }

            ApplyPreconditioner(preconditioner, residual, preconditioned, n, k);

            // Deflate the search direction. Order matters only for numerical
            // conditioning; the constraint projection is the bug-critical one —
            // without it a non-constant preconditioner leaks the deflated mode
            // back into the iterate.
            if (constraintCols > 0)
                ProjectBlockAgainst(preconditioned, n, k, constraints!, constraintCols);
            ProjectBlockAgainst(preconditioned, n, k, x, k);
            if (conjugateCols > 0)
                ProjectBlockAgainst(preconditioned, n, k, conjugate, conjugateCols);

            int trialCols = 0;
            trialCols = AppendBlock(trial, trialCols, x, n, k);
            trialCols = AppendBlock(trial, trialCols, preconditioned, n, k);
            if (conjugateCols > 0)
                trialCols = AppendBlock(trial, trialCols, conjugate, n, conjugateCols);

            int trialBasisCols = CompactOrthonormalize(
                trial.AsSpan(0, n * trialCols),
                n,
                trialCols,
                trialBasis.AsSpan(0, n * trialCols));

            op.Apply(trialBasis.AsSpan(0, n * trialBasisCols), appliedTrialBasis.AsSpan(0, n * trialBasisCols), trialBasisCols);
            RotateToExtremalRitzBasis(trialBasis, appliedTrialBasis, n, trialBasisCols, k, theta, xWork, axWork, wantLargest);

            BuildDifferenceBlock(xWork, x, conjugate, n, k);
            ProjectBlockAgainst(conjugate, n, k, xWork, k);
            conjugateCols = CompactOrthonormalize(conjugate, n, k, conjugate);

            Swap(ref x, ref xWork);
            Swap(ref ax, ref axWork);
        }

        var outputList = new List<EigenPair>(k);
        for (int j = 0; j < k; j++)
        {
            double[] vector = new double[n];
            Array.Copy(x, j * n, vector, 0, n);
            outputList.Add(new EigenPair(theta[j], vector));
        }

        // Smallest ⇒ ascending; largest ⇒ descending.
        outputList.Sort((left, right) => wantLargest
            ? right.Lambda.CompareTo(left.Lambda)
            : left.Lambda.CompareTo(right.Lambda));

        return new Result(outputList, completedIterations, lastResidual, converged);
    }

    /// <summary>
    /// The <c>k</c> smallest eigenpairs of a dense symmetric matrix, sorted
    /// ascending — the iterative peer of <see cref="SpectralMath.BottomK(double[,], int, DenseEigenOptions)"/>.
    /// No null-space deflation (that is a graph-Laplacian concern handled by
    /// <c>Graphs.Spectral</c>); this is the plain "k smallest of an arbitrary
    /// symmetric matrix" entry. Use <see cref="Solve"/> directly if you need the
    /// iteration/residual/convergence telemetry.
    /// </summary>
    public static IReadOnlyList<EigenPair> BottomK(double[,] matrix, int k, Options? options = null)
    {
        if (k <= 0) return Array.Empty<EigenPair>();
        var op = new DenseSymmetricOperator(matrix);
        return Solve(op, k, (options ?? new Options()) with { WantLargest = false }).Eigenpairs;
    }

    /// <summary>
    /// Flat column-major overload of <see cref="BottomK(double[,], int, Options)"/>;
    /// <paramref name="n"/> is the dimension and the span length must equal <c>n × n</c>.
    /// </summary>
    public static IReadOnlyList<EigenPair> BottomK(ReadOnlySpan<double> flatColumnMajorMatrix, int n, int k, Options? options = null)
    {
        if (k <= 0) return Array.Empty<EigenPair>();
        var op = new DenseSymmetricOperator(flatColumnMajorMatrix.ToArray(), n);
        return Solve(op, k, (options ?? new Options()) with { WantLargest = false }).Eigenpairs;
    }

    /// <summary>
    /// The <c>k</c> largest eigenpairs of a dense symmetric matrix, sorted
    /// descending. Same contract as <see cref="BottomK(double[,], int, Options)"/>
    /// with <see cref="Options.WantLargest"/> flipped.
    /// </summary>
    public static IReadOnlyList<EigenPair> TopK(double[,] matrix, int k, Options? options = null)
    {
        if (k <= 0) return Array.Empty<EigenPair>();
        var op = new DenseSymmetricOperator(matrix);
        return Solve(op, k, (options ?? new Options()) with { WantLargest = true }).Eigenpairs;
    }

    private static void ApplyPreconditioner(
        IPreconditioner? preconditioner, double[] residual, double[] preconditioned, int n, int k)
    {
        if (preconditioner is null)
        {
            Array.Copy(residual, preconditioned, n * k);
            return;
        }

        preconditioner.Apply(residual.AsSpan(0, n * k), preconditioned.AsSpan(0, n * k), k);
    }

    private static void RotateToExtremalRitzBasis(
        double[] vectors,
        double[] appliedVectors,
        int n,
        int sourceCols,
        int takeCols,
        double[] eigenvalues,
        double[] outputVectors,
        double[] outputAppliedVectors,
        bool wantLargest)
    {
        var projected = BuildProjectedMatrix(vectors, appliedVectors, n, sourceCols);
        var eig = Eigen.DecomposeSymmetric(projected); // eigenvalues sorted descending
        var coeff = BuildExtremalEigenvectorMatrix(eig, takeCols, wantLargest);

        MultiplyBlockByCoefficients(vectors, n, sourceCols, coeff, takeCols, outputVectors);
        MultiplyBlockByCoefficients(appliedVectors, n, sourceCols, coeff, takeCols, outputAppliedVectors);

        int lastIndex = eig.Eigenvalues.Length - 1;
        for (int col = 0; col < takeCols; col++)
            eigenvalues[col] = wantLargest ? eig.Eigenvalues[col] : eig.Eigenvalues[lastIndex - col];
    }

    private static double[,] BuildProjectedMatrix(double[] left, double[] right, int n, int cols)
    {
        var projected = new double[cols, cols];
        for (int i = 0; i < cols; i++)
        {
            ReadOnlySpan<double> leftCol = left.AsSpan(i * n, n);
            for (int j = i; j < cols; j++)
            {
                double value = MatrixOps.Dot(leftCol, right.AsSpan(j * n, n));
                projected[i, j] = value;
                projected[j, i] = value;
            }
        }

        return projected;
    }

    private static double[,] BuildExtremalEigenvectorMatrix(EigenResult eig, int takeCols, bool wantLargest)
    {
        int sourceCols = eig.Eigenvalues.Length;
        var coeff = new double[sourceCols, takeCols];
        int lastIndex = sourceCols - 1;

        for (int col = 0; col < takeCols; col++)
        {
            // Descending order: index 0 is the largest eigenpair, lastIndex the smallest.
            double[] eigenvector = wantLargest ? eig.Eigenvectors[col] : eig.Eigenvectors[lastIndex - col];
            for (int row = 0; row < sourceCols; row++)
                coeff[row, col] = eigenvector[row];
        }

        return coeff;
    }

    private static void MultiplyBlockByCoefficients(
        double[] source,
        int n,
        int sourceCols,
        double[,] coeff,
        int targetCols,
        double[] destination)
    {
        Array.Clear(destination, 0, n * targetCols);
        for (int targetCol = 0; targetCol < targetCols; targetCol++)
        {
            Span<double> dst = destination.AsSpan(targetCol * n, n);
            for (int sourceCol = 0; sourceCol < sourceCols; sourceCol++)
            {
                double weight = coeff[sourceCol, targetCol];
                if (weight == 0.0) continue;

                ReadOnlySpan<double> src = source.AsSpan(sourceCol * n, n);
                for (int row = 0; row < n; row++)
                    dst[row] += weight * src[row];
            }
        }
    }

    private static void ComputeResidualBlock(double[] appliedVectors, double[] vectors, double[] eigenvalues, double[] residual, int n, int cols)
    {
        for (int col = 0; col < cols; col++)
        {
            double lambda = eigenvalues[col];
            ReadOnlySpan<double> axCol = appliedVectors.AsSpan(col * n, n);
            ReadOnlySpan<double> xCol = vectors.AsSpan(col * n, n);
            Span<double> rCol = residual.AsSpan(col * n, n);

            for (int row = 0; row < n; row++)
                rCol[row] = axCol[row] - lambda * xCol[row];
        }
    }

    private static double ComputeMaxColumnNorm(double[] block, int n, int cols)
    {
        double max = 0.0;
        for (int col = 0; col < cols; col++)
        {
            ReadOnlySpan<double> vector = block.AsSpan(col * n, n);
            double norm = Math.Sqrt(MatrixOps.Dot(vector, vector));
            if (norm > max) max = norm;
        }

        return max;
    }

    private static int AppendBlock(double[] destination, int destinationCols, double[] source, int n, int sourceCols)
    {
        for (int col = 0; col < sourceCols; col++)
            source.AsSpan(col * n, n).CopyTo(destination.AsSpan((destinationCols + col) * n, n));
        return destinationCols + sourceCols;
    }

    private static int CompactOrthonormalize(ReadOnlySpan<double> source, int n, int sourceCols, Span<double> destination, double tolerance = 1e-12)
    {
        int acceptedCols = 0;

        for (int sourceCol = 0; sourceCol < sourceCols; sourceCol++)
        {
            Span<double> candidate = destination.Slice(acceptedCols * n, n);
            source.Slice(sourceCol * n, n).CopyTo(candidate);

            for (int pass = 0; pass < 2; pass++)
            {
                for (int basisCol = 0; basisCol < acceptedCols; basisCol++)
                {
                    ReadOnlySpan<double> basis = destination.Slice(basisCol * n, n);
                    double projection = MatrixOps.Dot(candidate, basis);
                    if (projection != 0.0)
                        MatrixOps.ScaledSubtract(candidate, basis, projection, n);
                }
            }

            double norm = Math.Sqrt(MatrixOps.Dot(candidate, candidate));
            if (norm < tolerance)
            {
                candidate.Clear();
                continue;
            }

            double invNorm = 1.0 / norm;
            for (int row = 0; row < n; row++)
                candidate[row] *= invNorm;

            acceptedCols++;
        }

        return acceptedCols;
    }

    private static void ProjectBlockAgainst(double[] block, int n, int blockCols, double[] basis, int basisCols)
    {
        if (basisCols <= 0) return;

        for (int blockCol = 0; blockCol < blockCols; blockCol++)
        {
            Span<double> vector = block.AsSpan(blockCol * n, n);
            for (int pass = 0; pass < 2; pass++)
            {
                for (int basisCol = 0; basisCol < basisCols; basisCol++)
                {
                    ReadOnlySpan<double> basisVector = basis.AsSpan(basisCol * n, n);
                    double projection = MatrixOps.Dot(vector, basisVector);
                    if (projection != 0.0)
                        MatrixOps.ScaledSubtract(vector, basisVector, projection, n);
                }
            }
        }
    }

    private static void BuildDifferenceBlock(double[] left, double[] right, double[] destination, int n, int cols)
    {
        for (int col = 0; col < cols; col++)
        {
            ReadOnlySpan<double> leftCol = left.AsSpan(col * n, n);
            ReadOnlySpan<double> rightCol = right.AsSpan(col * n, n);
            Span<double> dstCol = destination.AsSpan(col * n, n);

            for (int row = 0; row < n; row++)
                dstCol[row] = leftCol[row] - rightCol[row];
        }
    }

    private static void Swap(ref double[] left, ref double[] right)
    {
        double[] tmp = left;
        left = right;
        right = tmp;
    }
}
