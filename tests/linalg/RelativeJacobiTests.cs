using System;
using Maths.LinAlg;
using Xunit;

namespace Maths.LinAlg.Tests;

public sealed class RelativeJacobiTests
{

    [Fact]
    public void Dispatch_UsesRequestedBuildKernelAndRelativeRule()
    {
        using var log = DecompositionLog.Begin();
        var a = new[,] { { 1.0, 5e-21 }, { 5e-21, 1e-40 } };
        var e = DenseEigen.DecomposeSymmetric(a);
        var values = new double[2];
        DenseEigen.DecomposeSymmetricInto(a, values, null);
        Assert.Equal(e.Eigenvalues, values);
        Assert.InRange(Math.Abs(values[1]/7.5e-41-1),0,1e-14);
#if EIGEN_REFERENCE
        Assert.Equal(2, log.Count(DecompositionKernel.ReferenceJacobi));
        Assert.Equal(0, log.Count(DecompositionKernel.FastJacobi));
#else
        Assert.Equal(0, log.Count(DecompositionKernel.ReferenceJacobi));
        Assert.Equal(2, log.Count(DecompositionKernel.FastJacobi));
#endif
    }

    [Fact]
    public void DecidedWitnesses_KeepWeakEigenvaluesAcrossScales()
    {
        foreach (int exponent in new[] { -240, -120, -20, 0, 20, 120, 240 })
        {
            double scale = Math.Pow(10, exponent);
            Check(new[,] { { scale, 0, 0 }, { 0, scale * 1e-20, scale * 0.5e-20 }, { 0, scale * 0.5e-20, scale * 1e-20 } },
                [scale, scale * 1.5e-20, scale * 0.5e-20]);
            Check(new[,] { { 0.0, scale }, { scale, 0.0 } }, [scale, -scale]);
            Check(new[,] { { scale, scale * 5e-21 }, { scale * 5e-21, scale * 1e-40 } }, [scale, scale * 7.5e-41]);
        }
    }

    [Fact]
    public void DensePositiveAndIndefinite_ResidualsInOriginalUnits()
    {
        foreach (int exponent in new[] { -280, -200, -100, 0, 100, 200, 280 })
        {
            double s = Math.Pow(10, exponent);
            Check(new[,] { { 2 * s, s, 0 }, { s, 2 * s, 0 }, { 0, 0, 5 * s } }, [5 * s, 3 * s, s]);
            // Hadamard congruence of diag(7,4,2,-1).
            Check(new[,] { { 3*s, 1.5*s, 2.5*s, 0.0 }, { 1.5*s, 3*s, 0.0, 2.5*s },
                { 2.5*s, 0.0, 3*s, 1.5*s }, { 0.0, 2.5*s, 1.5*s, 3*s } }, [7*s, 4*s, 2*s, -s]);
        }
    }

    [Fact]
    public void DenseEightPairs_ExerciseVectorizedRotationsAndFinalReports()
    {
        double[] spectrum = [8, 7, 6, 5, 4, 3, 2, -1];
        var h = new double[8,8];
        h[0,0] = 1;
        for (int half = 1; half < 8; half *= 2)
            for (int i = 0; i < half; i++)
                for (int j = 0; j < half; j++)
                {
                    h[i,j+half] = h[i,j];
                    h[i+half,j] = h[i,j];
                    h[i+half,j+half] = -h[i,j];
                }
        var a = new double[8,8];
        for (int i = 0; i < 8; i++)
            for (int j = 0; j < 8; j++)
                for (int k = 0; k < 8; k++) a[i,j] += h[i,k]*spectrum[k]*h[j,k]/8;
        Check(a,spectrum);
        EigenResult capped = DenseEigen.DecomposeSymmetric(a, maxSweeps: 1, tol: 1e-14);
        var flat = new double[64];
        for (int j = 0; j < 8; j++)
            for (int i = 0; i < 8; i++) flat[j*8+i] = a[i,j];
        var values = new double[8];
        Assert.Equal(capped.Convergence, DenseEigen.DecomposeSymmetric(flat,8,values,[],maxSweeps: 1,tol: 1e-14));
        double squares = 0, largest = 0;
        for (int p = 0; p < 8; p++)
            for (int q = p+1; q < 8; q++)
            {
                double pivot = 0;
                for (int i = 0; i < 8; i++)
                    for (int j = 0; j < 8; j++) pivot += capped.Eigenvectors[p][i]*a[i,j]*capped.Eigenvectors[q][j];
                largest = Math.Max(largest,Math.Abs(pivot));
                squares += 2*pivot*pivot;
            }
        Assert.InRange(Math.Abs(largest-capped.Convergence.LargestOffDiagonal),0,1e-13);
        Assert.InRange(Math.Abs(Math.Sqrt(squares)-capped.Convergence.OffDiagonalNorm),0,1e-13);
        Assert.False(capped.Convergence.Converged);
        Assert.False(capped.Convergence.Stalled);
    }

    [Fact]
    public void FailedAccuracyTest_WithUnrepresentableRotation_IsNotConvergence()
    {
        var a = new[,] { { 1e300, 1e-300 }, { 1e-300, 0.0 } };
        foreach (EigenResult e in SolveAll(a))
        {
            Assert.False(e.Convergence.Converged);
            Assert.True(e.Convergence.Stalled);
            Assert.Equal(1, e.Convergence.Sweeps);
            Assert.Equal(0L, e.Convergence.Rotations);
            Assert.Equal(1e-300, e.Convergence.LargestOffDiagonal);
        }
        double[] flat = [1e300, 1e-300, 1e-300, 0];
        foreach (bool vectors in new[] { false, true })
        {
            var values = new double[2];
            DenseEigenConvergence c = DenseEigen.DecomposeSymmetric(flat, 2, values, vectors ? new double[4] : []);
            Assert.False(c.Converged);
            Assert.True(c.Stalled);
            Assert.Equal(1e-300, c.LargestOffDiagonal);
        }
    }

    [Fact]
    public void CapAndZeroTolerance_ReportFinalPivotTests()
    {
        var a = new[,] { { 4.0, 1, 2 }, { 1, 3, 0.5 }, { 2, 0.5, 2 } };
        foreach (var solve in Solvers())
        {
            Assert.False(solve(a, 0, 1e-14).Convergence.Converged);
            Assert.False(solve(a, 1, 1e-14).Convergence.Converged);
            Assert.True(solve(new double[3,3], 1, 0).Convergence.Converged);
        }
        // Equality belongs to the relative pass test, including a nonzero pivot exactly on its threshold.
        foreach (var e in SolveAll(new[,] { { 1.0, 1e-14 }, { 1e-14, 1.0 } }))
        {
            Assert.True(e.Convergence.Converged);
            Assert.Equal(0L, e.Convergence.Rotations);
        }
    }

    private static Func<double[,], int, double, EigenResult>[] Solvers() =>
    [
        (a, cap, tol) => EigenFast.DecomposeSymmetric(a, cap, tol),
        (a, cap, tol) => Eigen.DecomposeSymmetric(a, cap, tol),
        (a, cap, tol) => DenseEigen.DecomposeSymmetric(a, cap, tol),
    ];

    private static EigenResult[] SolveAll(double[,] a)
    {
        var solvers = Solvers();
        return [solvers[0](a, 256, 1e-14), solvers[1](a, 256, 1e-14), solvers[2](a, 256, 1e-14)];
    }

    private static void Check(double[,] a, double[] expected)
    {
        foreach (EigenResult e in SolveAll(a)) Verify(a, expected, e);
        int n = expected.Length;
        var flat = new double[n*n];
        for (int j = 0; j < n; j++)
            for (int i = 0; i < n; i++) flat[j*n+i] = a[i,j];
        EigenResult allocated = DenseEigen.DecomposeSymmetric(flat, n, tol: 1e-14);
        Verify(a, expected, allocated);
        var values = new double[n];
        var columns = new double[n,n];
        DenseEigenConvergence c = DenseEigen.DecomposeSymmetricInto(a, values, columns, tol: 1e-14);
        Assert.Equal(allocated.Eigenvalues, values);
        Assert.Equal(allocated.Convergence, c);
        for (int j = 0; j < n; j++)
            for (int i = 0; i < n; i++) Assert.Equal(allocated.Eigenvectors[j][i], columns[i,j]);
        var alone = new double[n];
        Assert.Equal(c, DenseEigen.DecomposeSymmetric(flat, n, alone, [], tol: 1e-14));
        Assert.Equal(values, alone);
    }

    private static void Verify(double[,] a, double[] expected, EigenResult e)
    {
        Assert.True(e.Convergence.Converged);
        int n = expected.Length;
        double scale = Math.Abs(expected[0]);
        for (int j = 0; j < n; j++)
        {
            Assert.InRange(Math.Abs(e.Eigenvalues[j] / expected[j] - 1), 0, 2e-13);
            for (int i = 0; i < n; i++)
            {
                double r = -e.Eigenvalues[j] / scale * e.Eigenvectors[j][i];
                for (int k = 0; k < n; k++) r += a[i,k] / scale * e.Eigenvectors[j][k];
                Assert.InRange(Math.Abs(r), 0, 2e-13); // residual bound in original units: 2e-13*scale
                double dot = 0;
                for (int k = 0; k < n; k++) dot += e.Eigenvectors[i][k] * e.Eigenvectors[j][k];
                Assert.InRange(Math.Abs(dot - (i == j ? 1 : 0)), 0, 2e-13);
            }
        }
    }
}
