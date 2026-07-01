using System;
using System.Diagnostics;
using Maths.LinAlg;
using Xunit;
using Xunit.Abstractions;

namespace Maths.LinAlg.Tests;

public sealed class GramMatrixBenchmarks : IDisposable
{
    private readonly ITestOutputHelper _output;
    private readonly Stopwatch _sw = new();

    public GramMatrixBenchmarks(ITestOutputHelper output) => _output = output;

    [Theory]
    [InlineData(100, 50)]   // Small: n=100, p=50
    [InlineData(500, 100)]  // Medium: n=500, p=100  
    [InlineData(1000, 200)] // Large: n=1000, p=200
    [InlineData(2000, 300)] // XL: n=2000, p=300
    public void GramMatrixVsNaive_HeadToHead(int n, int p)
    {
        // Generate random row-major data
        var rng = new Random(42);
        double[][] data = new double[n][];
        for (int i = 0; i < n; i++)
        {
            data[i] = new double[p];
            for (int j = 0; j < p; j++)
                data[i][j] = rng.NextDouble();
        }

        // Compute column means for centering
        double[] mean = new double[p];
        for (int i = 0; i < n; i++)
            for (int j = 0; j < p; j++)
                mean[j] += data[i][j];
        for (int j = 0; j < p; j++) mean[j] /= n;

        // --- NAIVE: Per-pair inner products (O(p^2 n) redundant) ---
        _sw.Restart();
        double[,] gramNaive = new double[p, p];
        for (int i = 0; i < p; i++)
        {
            for (int j = i; j < p; j++)
            {
                double sum = 0;
                for (int k = 0; k < n; k++)
                {
                    double vi = data[k][i] - mean[i];
                    double vj = data[k][j] - mean[j];
                    sum += vi * vj;
                }
                gramNaive[i, j] = gramNaive[j, i] = sum;
            }
        }
        _sw.Stop();
        long naiveMs = _sw.ElapsedMilliseconds;

        // --- SIMD: Column-major transpose + SIMD dot ---
        _sw.Restart();
        double[][] cols = new double[p][];
        for (int j = 0; j < p; j++)
        {
            cols[j] = new double[n];
            double mu = mean[j];
            for (int i = 0; i < n; i++)
                cols[j][i] = data[i][j] - mu;
        }
        double[,] gramSimd = MatrixOps.ColumnGramMatrix(cols, n, p);
        _sw.Stop();
        long simdMs = _sw.ElapsedMilliseconds;

        // Verify correctness
        for (int i = 0; i < p; i++)
            for (int j = 0; j < p; j++)
                Assert.True(Math.Abs(gramNaive[i, j] - gramSimd[i, j]) < 1e-6,
                    $"Mismatch at ({i},{j}): naive={gramNaive[i,j]}, simd={gramSimd[i,j]}");

        double speedup = naiveMs / (double)Math.Max(simdMs, 1);
        _output.WriteLine($"n={n}, p={p}: Naive={naiveMs}ms, SIMD={simdMs}ms, Speedup={speedup:F2}x");

        // Assert SIMD wins or is competitive (allow 20% tolerance for small cases)
        if (n * p > 10000) // Only assert for reasonably large problems
            Assert.True(simdMs < naiveMs * 1.2, $"SIMD should beat or match naive for n={n}, p={p}");
    }

    [Fact]
    public void ResidualVarFromGram_Correctness()
    {
        // Verify the residual variance formula: tau = (G[i,i] - G[i,j]^2/G[j,j]) / n
        int n = 100, p = 10;
        var rng = new Random(42);
        double[][] data = new double[n][];
        for (int i = 0; i < n; i++)
        {
            data[i] = new double[p];
            for (int j = 0; j < p; j++)
                data[i][j] = rng.NextDouble();
        }

        // Center
        double[] mean = new double[p];
        for (int i = 0; i < n; i++)
            for (int j = 0; j < p; j++)
                mean[j] += data[i][j];
        for (int j = 0; j < p; j++) mean[j] /= n;

        // Transpose + Gram
        double[][] cols = new double[p][];
        for (int j = 0; j < p; j++)
        {
            cols[j] = new double[n];
            double mu = mean[j];
            for (int i = 0; i < n; i++)
                cols[j][i] = data[i][j] - mu;
        }
        double[,] gram = MatrixOps.ColumnGramMatrix(cols, n, p);

        // Compare tau from Gram vs direct computation
        for (int i = 0; i < p; i++)
        {
            for (int j = 0; j < p; j++)
            {
                if (i == j) continue;

                // Direct
                double sii = 0, sjj = 0, sij = 0;
                for (int k = 0; k < n; k++)
                {
                    double vi = data[k][i] - mean[i];
                    double vj = data[k][j] - mean[j];
                    sii += vi * vi;
                    sjj += vj * vj;
                    sij += vi * vj;
                }
                double tauDirect = (sjj > 1e-12 ? sii - sij * sij / sjj : sii) / n;

                // From Gram
                double gii = gram[i, i];
                double gjj = gram[j, j];
                double gij = gram[i, j];
                double tauGram = (gjj > 1e-12 ? gii - gij * gij / gjj : gii) / n;

                Assert.True(Math.Abs(tauDirect - tauGram) < 1e-9,
                    $"tau mismatch at ({i},{j}): direct={tauDirect}, gram={tauGram}");
            }
        }
    }

    public void Dispose() { }
}
