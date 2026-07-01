using System;
using System.Collections.Generic;
using System.Linq;
using Graphs.Spectral;
using Maths.LinAlg;
using Xunit;

namespace VizCore.Tests;

/// <summary>
/// M1 gate for the magnetic Laplacian (Seam A of the complex-analytic TDA thread):
/// the matrix-free Hermitian SpMV must match its dense embedding; <c>q=0</c> must
/// recover the ordinary undirected Laplacian; and fractional <c>q</c> must lift the
/// harmonic zero-mode off zero (Aharonov–Bohm flux frustrating the directed cycle).
/// The closed-form directed-cycle spectrum <c>λ_m = 2 − 2cos(2π(m/n + q))</c> is the
/// independent oracle — derived outside the implementation, not shared with it.
/// </summary>
public sealed class MagneticLaplacianTests
{
    private const int N = 6; // directed cycle 0→1→…→5→0

    [Fact]
    public void Apply_MatchesDenseEmbedding()
    {
        var op = MagneticLaplacianOperator.FromDirectedEdges(N, DirectedCycle(N), charge: 0.137);
        double[,] dense = op.ToDenseEmbedding();
        int dim = op.Dimension; // 2n

        // Embedding must be symmetric (LOBPCG assumes a symmetric operator).
        for (int i = 0; i < dim; i++)
            for (int j = 0; j < dim; j++)
                Assert.InRange(Math.Abs(dense[i, j] - dense[j, i]), 0.0, 1e-15);

        // Random 2-column block; matrix-free Apply must equal the dense matvec.
        const int cols = 2;
        var rng = new Random(3);
        var block = new double[dim * cols];
        for (int t = 0; t < block.Length; t++) block[t] = rng.NextDouble() * 2.0 - 1.0;

        var got = new double[dim * cols];
        op.Apply(block, got, cols);

        for (int c = 0; c < cols; c++)
            for (int i = 0; i < dim; i++)
            {
                double expected = 0.0;
                for (int j = 0; j < dim; j++)
                    expected += dense[i, j] * block[c * dim + j];
                Assert.InRange(Math.Abs(expected - got[c * dim + i]), 0.0, 1e-12);
            }
    }

    [Fact]
    public void ZeroCharge_RecoversUndirectedCycleSpectrum()
    {
        var op = MagneticLaplacianOperator.FromDirectedEdges(N, DirectedCycle(N), charge: 0.0);

        // Dense full spectrum (each eigenvalue doubled by the embedding) → stride 2.
        double[] distinct = Stride2(SpectralMath.BottomK(op.ToDenseEmbedding(), 2 * N));
        double[] analytic = AnalyticCycleSpectrum(N, 0.0);

        Assert.Equal(analytic.Length, distinct.Length);
        for (int m = 0; m < analytic.Length; m++)
            Assert.InRange(Math.Abs(distinct[m] - analytic[m]), 0.0, 1e-9);

        // The undirected cycle has a harmonic zero-mode.
        Assert.InRange(distinct[0], -1e-9, 1e-9);

        // Engine path: LOBPCG over the operator finds the same ~0 smallest eigenvalue.
        double[] engine = MagneticSpectral.BottomKEigenvalues(op, 1, SolverOptions());
        Assert.InRange(engine[0], -1e-6, 1e-6);
    }

    [Fact]
    public void FractionalCharge_LiftsZeroMode_MatchesAnalytic()
    {
        const double q = 1.0 / 12.0; // fractional flux: smallest mode frustrated off zero
        var op = MagneticLaplacianOperator.FromDirectedEdges(N, DirectedCycle(N), q);

        double[] analytic = AnalyticCycleSpectrum(N, q);
        double analyticMin = analytic[0]; // = 2 − 2cos(π/6) ≈ 0.2679

        // Engine path: the smallest eigenvalue is lifted and matches the oracle.
        double[] engine = MagneticSpectral.BottomKEigenvalues(op, 1, SolverOptions());
        Assert.True(engine[0] > 0.1, $"zero-mode not lifted: λ_min = {engine[0]}");
        Assert.InRange(Math.Abs(engine[0] - analyticMin), 0.0, 1e-4);

        // Full deduplicated spectrum matches the closed form.
        double[] distinct = Stride2(SpectralMath.BottomK(op.ToDenseEmbedding(), 2 * N));
        Assert.Equal(analytic.Length, distinct.Length);
        for (int m = 0; m < analytic.Length; m++)
            Assert.InRange(Math.Abs(distinct[m] - analytic[m]), 0.0, 1e-9);
    }

    private static (int from, int to)[] DirectedCycle(int n)
    {
        var edges = new (int, int)[n];
        for (int i = 0; i < n; i++) edges[i] = (i, (i + 1) % n);
        return edges;
    }

    private static double[] AnalyticCycleSpectrum(int n, double q)
    {
        var values = new double[n];
        for (int m = 0; m < n; m++)
            values[m] = 2.0 - 2.0 * Math.Cos(2.0 * Math.PI * ((double)m / n + q));
        Array.Sort(values);
        return values;
    }

    // The embedding doubles every eigenvalue, so identical values are adjacent in the
    // ascending dense spectrum; striding by two recovers the distinct H-spectrum.
    private static double[] Stride2(IReadOnlyList<EigenPair> ascending)
    {
        var outv = new List<double>(ascending.Count / 2);
        for (int i = 0; i < ascending.Count; i += 2) outv.Add(ascending[i].Lambda);
        return outv.ToArray();
    }

    private static LOBPCG.Options SolverOptions() =>
        new() { MaxIterations = 2000, Tolerance = 1e-11 };
}
