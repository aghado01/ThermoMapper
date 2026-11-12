using System;
using System.Linq;
using Graphs.Primitives;
using Maths.LinAlg;
using Xunit;

namespace VizCore.Tests;

/// <summary>
/// EigenFast's contract, read from the decomposition itself (spectral-solver brief, A2): eigenvalues descending, each pair's
/// residual ‖Av − λv‖ and the eigenvectors' orthogonality at rounding level, known spectra recovered, and the convergence it
/// reports describing the diagonalization it returns. One parity check against the reference Jacobi stays, as a qualification.
/// </summary>
public sealed class EigenFastTests
{
    private const double PivotTolerance = 1e-12;   // the solver's default

    /// <summary>A 2 × 2 block with equal diagonals, θ = 0 in the rotation formula, beside a decoupled entry: spectrum {5, 3, 1}.</summary>
    [Fact]
    public void EqualDiagonalBlock_RecoversItsSpectrum()
    {
        var matrix = new double[,] { { 2.0, 1.0, 0.0 }, { 1.0, 2.0, 0.0 }, { 0.0, 0.0, 5.0 } };
        EigenResult result = EigenFast.DecomposeSymmetric(matrix);

        AssertDecomposition(matrix, result);
        AssertSpectrum([5.0, 3.0, 1.0], result);
    }

    [Fact]
    public void NonSquare_Throws() => Assert.Throws<ArgumentException>(() => EigenFast.DecomposeSymmetric(new double[2, 3]));

    [Fact]
    public void KnownSpectra_Recovered()
    {
        (double[] Spectrum, double[,] Basis)[] fixtures =
        [
            ([5.0, 3.0, 1.0], PlanarRotation3(Math.PI / 6.0)),
            ([7.0, 4.0, 2.0, -1.0], NormalizedHadamard(4)),
            ([12.0, 9.0, 7.0, 5.0, 4.0, 3.0, 2.0, 0.5], NormalizedHadamard(8)),
            ([6.0, 6.0, 2.0, 1.0], NormalizedHadamard(4)),   // a repeated eigenvalue: its eigenspace, not its vectors, is determined
        ];
        foreach ((double[] spectrum, double[,] basis) in fixtures)
        {
            double[,] matrix = WithSpectrum(spectrum, basis);
            EigenResult result = EigenFast.DecomposeSymmetric(matrix);
            AssertDecomposition(matrix, result);
            AssertSpectrum(spectrum, result);
        }
    }

    /// <summary>A graph with two components has a Laplacian with two null modes.</summary>
    [Fact]
    public void DisconnectedLaplacian_HasTwoNullModes()
    {
        DisconnectedControl.Fixture control = DisconnectedControl.Generate();
        CsrGraph graph = FixtureGoldenHelpers.BuildDraftGraph(control.Points, k: 4, ensureConnected: false);
        double[,] laplacian = CombinatorialLaplacian(graph);
        EigenResult result = EigenFast.DecomposeSymmetric(laplacian);

        AssertDecomposition(laplacian, result);
        Assert.InRange(Math.Abs(result.Eigenvalues[^1]), 0.0, 1e-10);
        Assert.InRange(Math.Abs(result.Eigenvalues[^2]), 0.0, 1e-10);
    }

    /// <summary>A dense random matrix needs rotations: the solve converges in a handful of sweeps, every pivot of the returned
    /// diagonalization below the tolerance.</summary>
    [Fact]
    public void RandomSymmetric_ConvergesWithRotations()
    {
        double[,] matrix = RandomSymmetric(size: 50, seed: 1729);
        EigenResult result = EigenFast.DecomposeSymmetric(matrix);
        DenseEigenConvergence c = result.Convergence;

        AssertDecomposition(matrix, result);
        Assert.True(c.Converged);
        Assert.True(c.Rotations > 0);
        Assert.InRange(c.Sweeps, 2, 30);
        Assert.InRange(c.LargestOffDiagonal, 0.0, PivotTolerance);
        Assert.InRange(c.OffDiagonalNorm, 0.0, Math.Sqrt(50.0 * 49.0) * PivotTolerance);
        AssertReportedOffDiagonal(matrix, result);
    }

    /// <summary>
    /// The brief's witness: diagonal 1 to 32 and every off-diagonal 5·10⁻¹³, below the pivot tolerance each, while off(A) is about
    /// 1.6·10⁻¹¹, above it. Every pivot is skipped; the solve used to run to its cap of sweeps, rotating nothing. A sweep that
    /// rotates nothing now ends it: one sweep, no rotation, the same eigenpairs as under a cap of one.
    /// </summary>
    [Fact]
    public void SweepWithNoRotations_EndsTheSolve()
    {
        const int n = 32;
        var matrix = new double[n, n];
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++) matrix[i, j] = i == j ? i + 1.0 : 5e-13;
        EigenResult result = EigenFast.DecomposeSymmetric(matrix);
        EigenResult capped = EigenFast.DecomposeSymmetric(matrix, maxSweeps: 1);
        DenseEigenConvergence c = result.Convergence;

        Assert.Equal(new DenseEigenConvergence(1, 0, c.OffDiagonalNorm, 5e-13, true), c);
        Assert.Equal(Math.Sqrt(n * (n - 1.0)) * 5e-13, c.OffDiagonalNorm, 1e-24);
        Assert.True(c.OffDiagonalNorm > PivotTolerance);
        Assert.Equal(capped.Eigenvalues, result.Eigenvalues);
        Assert.Equal(capped.Convergence, c);
        for (int i = 0; i < n; i++)
        {
            Assert.Equal((double)(n - i), result.Eigenvalues[i]);
            Assert.Equal(capped.Eigenvectors[i], result.Eigenvectors[i]);
            Assert.InRange(Residual(matrix, result.Eigenvalues[i], result.Eigenvectors[i]), 0.0, Math.Sqrt(n - 1.0) * 5e-13 * 1.01);
        }
    }

    /// <summary>Capped at one sweep, a dense random matrix has not converged, and the reported off-diagonal is the returned
    /// diagonalization's: VᵀAV, recomputed from the result, agrees with it.</summary>
    [Fact]
    public void SweepCap_ReportsTheReturnedDiagonalization()
    {
        double[,] matrix = RandomSymmetric(size: 50, seed: 1729);
        EigenResult result = EigenFast.DecomposeSymmetric(matrix, maxSweeps: 1);
        DenseEigenConvergence c = result.Convergence;

        Assert.Equal(1, c.Sweeps);
        Assert.Equal(50L * 49 / 2, c.Rotations);   // every pivot of a dense random matrix rotates on the first sweep
        Assert.False(c.Converged);
        Assert.True(c.LargestOffDiagonal >= PivotTolerance);
        AssertReportedOffDiagonal(matrix, result);
    }

    /// <summary>Qualification: the reference Jacobi gives the same spectrum on a dense random matrix.</summary>
    [Fact]
    public void RandomSymmetric_AgreesWithReference()
    {
        double[,] matrix = RandomSymmetric(size: 50, seed: 1729);
        EigenResult reference = Eigen.DecomposeSymmetric(matrix);
        EigenResult fast = EigenFast.DecomposeSymmetric(matrix);

        Assert.True(reference.Convergence.Converged);
        for (int i = 0; i < 50; i++) Assert.InRange(Math.Abs(reference.Eigenvalues[i] - fast.Eigenvalues[i]), 0.0, 1e-12);
    }

    /// <summary>
    /// The caller-owned entries (spectral-solver brief, B2) give the allocating overload's bits: eigenvectors as flat columns, or as
    /// the columns of a double[,]; and eigenvalues alone, which skip the eigenvector rotations, are the same eigenvalues and the same
    /// convergence report.
    /// </summary>
    [Fact]
    public void CallerOwnedEntries_MatchTheAllocatingOverload()
    {
        const int n = 40;
        double[,] matrix = RandomSymmetric(size: n, seed: 271);
        double[] flat = new double[n * n];
        for (int c = 0; c < n; c++)
            for (int r = 0; r < n; r++) flat[c * n + r] = matrix[r, c];
        EigenResult allocated = EigenFast.DecomposeSymmetric(matrix);

        var values = new double[n];
        var vectors = new double[n * n];
        DenseEigenConvergence withVectors = EigenFast.DecomposeSymmetricInPlace(flat, n, values, vectors);
        var alone = new double[n];
        DenseEigenConvergence valuesOnly = EigenFast.DecomposeSymmetricInPlace(flat, n, alone, Span<double>.Empty);

        Assert.Equal(allocated.Eigenvalues, values);
        Assert.Equal(allocated.Eigenvalues, alone);
        Assert.Equal(allocated.Convergence, withVectors);
        Assert.Equal(withVectors, valuesOnly);
        for (int j = 0; j < n; j++) Assert.Equal(allocated.Eigenvectors[j], vectors.AsSpan(j * n, n).ToArray());

        EigenResult dispatched = DenseEigen.DecomposeSymmetric(matrix);
        var columnValues = new double[n];
        var columns = new double[n, n];
        DenseEigen.DecomposeSymmetricInto(matrix, columnValues, columns);
        Assert.Equal(dispatched.Eigenvalues, columnValues);
        for (int j = 0; j < n; j++)
            for (int i = 0; i < n; i++) Assert.Equal(dispatched.Eigenvectors[j][i], columns[i, j]);
    }

    /// <summary>A repeated eigenvalue's vectors keep the order the exchange sort always gave: diag(1, 1, 2) rotates nothing, and the
    /// exchanges leave the columns as e₂, e₁, e₀.</summary>
    [Fact]
    public void RepeatedEigenvalues_KeepTheExchangeOrder()
    {
        EigenResult result = EigenFast.DecomposeSymmetric(new[,] { { 1.0, 0.0, 0.0 }, { 0.0, 1.0, 0.0 }, { 0.0, 0.0, 2.0 } });

        Assert.Equal([2.0, 1.0, 1.0], result.Eigenvalues);
        Assert.Equal([0.0, 0.0, 1.0], result.Eigenvectors[0]);
        Assert.Equal([0.0, 1.0, 0.0], result.Eigenvectors[1]);
        Assert.Equal([1.0, 0.0, 0.0], result.Eigenvectors[2]);
    }

    /// <summary>Into caller-owned buffers a solve allocates nothing once the pools are warm: the working copy and the sort's scratch
    /// are pooled, and the convergence report is read without a closure.</summary>
    [Fact]
    public void CallerOwnedSolve_AllocatesNothingWhenWarm()
    {
        const int n = 24;
        double[,] matrix = RandomSymmetric(size: n, seed: 314);
        double[] flat = new double[n * n];
        for (int c = 0; c < n; c++)
            for (int r = 0; r < n; r++) flat[c * n + r] = matrix[r, c];
        var values = new double[n];
        var vectors = new double[n * n];
        for (int warm = 0; warm < 3; warm++) EigenFast.DecomposeSymmetricInPlace(flat, n, values, vectors);

        long least = long.MaxValue;
        for (int trial = 0; trial < 5; trial++)
        {
            long before = GC.GetAllocatedBytesForCurrentThread();
            EigenFast.DecomposeSymmetricInPlace(flat, n, values, vectors);
            least = Math.Min(least, GC.GetAllocatedBytesForCurrentThread() - before);
        }
        Assert.Equal(0L, least);
    }

    // Descending eigenvalues, residuals ‖Av − λv‖ at 10⁻¹⁰·max(1, ‖A‖_F) and orthogonality |VᵀV − I| at 10⁻¹².
    private static void AssertDecomposition(double[,] matrix, EigenResult result)
    {
        int n = matrix.GetLength(0);
        Assert.Equal(n, result.Eigenvalues.Length);
        Assert.Equal(n, result.Eigenvectors.Length);
        double scale = Math.Max(1.0, Math.Sqrt(matrix.Cast<double>().Sum(x => x * x)));
        for (int i = 0; i < n; i++)
        {
            if (i > 0) Assert.True(result.Eigenvalues[i - 1] >= result.Eigenvalues[i], "Eigenvalues are sorted descending.");
            Assert.InRange(Residual(matrix, result.Eigenvalues[i], result.Eigenvectors[i]), 0.0, 1e-10 * scale);
            for (int j = 0; j <= i; j++)
            {
                double dot = result.Eigenvectors[i].Zip(result.Eigenvectors[j], (a, b) => a * b).Sum();
                Assert.InRange(Math.Abs(dot - (i == j ? 1.0 : 0.0)), 0.0, 1e-12);
            }
        }
    }

    private static void AssertSpectrum(double[] spectrum, EigenResult result)
    {
        double[] expected = [.. spectrum.OrderByDescending(x => x)];
        double scale = expected.Max(Math.Abs);
        for (int i = 0; i < expected.Length; i++) Assert.InRange(Math.Abs(expected[i] - result.Eigenvalues[i]), 0.0, 1e-12 * scale);
    }

    // off(VᵀAV) from the result, against the reported off(A) of the solver's diagonalization: the two differ by rounding in
    // forming VᵀAV and in the solver's rotations, both about n·ε·‖A‖ per entry.
    private static void AssertReportedOffDiagonal(double[,] matrix, EigenResult result)
    {
        int n = matrix.GetLength(0);
        double squares = 0.0;
        for (int i = 0; i < n; i++)
        {
            double[] av = Multiply(matrix, result.Eigenvectors[i]);
            for (int j = 0; j < n; j++)
            {
                if (j == i) continue;
                double entry = result.Eigenvectors[j].Zip(av, (a, b) => a * b).Sum();
                squares += entry * entry;
            }
        }
        double recomputed = Math.Sqrt(squares);
        double frobenius = Math.Sqrt(matrix.Cast<double>().Sum(x => x * x));
        Assert.InRange(Math.Abs(recomputed - result.Convergence.OffDiagonalNorm), 0.0, 1e-11 * frobenius);
    }

    private static double Residual(double[,] matrix, double eigenvalue, double[] eigenvector)
    {
        double[] av = Multiply(matrix, eigenvector);
        double squares = 0.0;
        for (int row = 0; row < av.Length; row++)
        {
            double r = av[row] - eigenvalue * eigenvector[row];
            squares += r * r;
        }
        return Math.Sqrt(squares);
    }

    private static double[] Multiply(double[,] matrix, double[] vector)
    {
        int n = vector.Length;
        var product = new double[n];
        for (int row = 0; row < n; row++)
            for (int col = 0; col < n; col++) product[row] += matrix[row, col] * vector[col];
        return product;
    }

    private static double[,] RandomSymmetric(int size, int seed)
    {
        var random = new Random(seed);
        var matrix = new double[size, size];
        for (int i = 0; i < size; i++)
            for (int j = i; j < size; j++) matrix[i, j] = matrix[j, i] = random.NextDouble() * 2.0 - 1.0;
        return matrix;
    }

    // B·diag(λ)·Bᵀ for an orthogonal B.
    private static double[,] WithSpectrum(double[] eigenvalues, double[,] basis)
    {
        int n = eigenvalues.Length;
        var matrix = new double[n, n];
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                for (int k = 0; k < n; k++) matrix[i, j] += basis[i, k] * eigenvalues[k] * basis[j, k];
        return matrix;
    }

    private static double[,] PlanarRotation3(double theta)
    {
        double c = Math.Cos(theta), s = Math.Sin(theta);
        return new[,] { { c, -s, 0.0 }, { s, c, 0.0 }, { 0.0, 0.0, 1.0 } };
    }

    private static double[,] NormalizedHadamard(int size)
    {
        var h = new double[size, size];
        h[0, 0] = 1.0;
        for (int half = 1; half < size; half *= 2)
            for (int i = 0; i < half; i++)
                for (int j = 0; j < half; j++)
                {
                    h[i, j + half] = h[i, j];
                    h[i + half, j] = h[i, j];
                    h[i + half, j + half] = -h[i, j];
                }
        double scale = 1.0 / Math.Sqrt(size);
        for (int i = 0; i < size; i++)
            for (int j = 0; j < size; j++) h[i, j] *= scale;
        return h;
    }

    private static double[,] CombinatorialLaplacian(CsrGraph graph)
    {
        int n = graph.NodeCount;
        var laplacian = new double[n, n];
        for (int i = 0; i < n; i++)
        {
            double degree = 0.0;
            for (int edge = graph.RowPointers[i]; edge < graph.RowPointers[i + 1]; edge++)
            {
                degree += graph.Weights[edge];
                laplacian[i, graph.Targets[edge]] -= graph.Weights[edge];
            }
            laplacian[i, i] = degree;
        }
        return laplacian;
    }
}
