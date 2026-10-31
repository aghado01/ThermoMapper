using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using Maths.LinAlg;
using Maths.Rng;
using Xunit;
using Xunit.Abstractions;

namespace Maths.LinAlg.Tests;

/// <summary>
/// The wide-band factor (N-X2) against the dense one on block-banded Grams, the tensor designs of
/// <see cref="KhatriRao.FaceSplittingGram"/>. Where the matrix is zero beyond the band the two factors do the same operations
/// on the same nonzero terms in the same order, so every comparison here is exact, not within a tolerance:
/// <list type="bullet">
/// <item>the factor's entries, its log-determinant, its solves, its quadratic forms, its draws from one z, its status and its
/// condition estimate;</item>
/// <item>at the band's edge: the declared bandwidth is the Gram's, a wider declaration changes nothing, a full band is the
/// dense factor, and an entry beyond the band is no part of the matrix the factor sees;</item>
/// <item>near singularity: a nearly dependent pair of columns reports the dense factor's tiny relative pivot, an exactly
/// dependent one or an indefinite one stops at the dense factor's row, with its status, and the factor refuses to solve.</item>
/// </list>
/// The selection rule is checked at its boundaries and for what it returns.
/// </summary>
public sealed class WideBandCholeskyTests
{
    private readonly ITestOutputHelper _out;
    public WideBandCholeskyTests(ITestOutputHelper output) => _out = output;

    // A hat basis on `width` nodes: two adjacent nonzeros a row.
    private static double[,] Hats(double[] x, int width)
    {
        var b = new double[x.Length, width];
        for (int i = 0; i < x.Length; i++)
        {
            double t = x[i] * (width - 1);
            int j = Math.Min((int)t, width - 2);
            b[i, j] = 1.0 - (t - j);
            b[i, j + 1] = t - j;
        }
        return b;
    }

    // The Gram of a two-axis tensor design on scattered points, its half-bandwidth: block-banded, 9 × 8 nodes.
    private static (double[,] Gram, int Bandwidth) TensorGram(int seed = 1401, int nu1 = 9, int nu2 = 8, int samples = 700)
    {
        var rng = new Xoshiro256PlusPlus(seed);
        IReadOnlyList<double[,]> axes = [Hats([.. Enumerable.Range(0, samples).Select(_ => rng.NextDouble())], nu1),
                                         Hats([.. Enumerable.Range(0, samples).Select(_ => rng.NextDouble())], nu2)];
        return (KhatriRao.FaceSplittingGram(axes), KhatriRao.FaceSplittingHalfBandwidth(axes));
    }

    private static double[] Noise(int n, int seed)
    {
        var rng = new Xoshiro256PlusPlus(seed);
        return [.. Enumerable.Range(0, n).Select(_ => rng.NextDouble() - 0.5)];
    }

    private static void AssertSameFactor(DenseCholesky dense, WideBandCholesky band, int seed)
    {
        int n = dense.Dimension;
        Assert.Equal(n, band.Dimension);
        Assert.Equal(dense.Status, band.Status);
        Assert.Equal(dense.LogDet, band.LogDet);
        Assert.Equal(dense.FailedAt, band.FailedAt);
        for (int i = 0; i < n; i++)
        {
            Assert.Equal(dense.RelativePivot(i), band.RelativePivot(i));
            for (int j = 0; j < n; j++) Assert.Equal(dense.Factor(i, j), band.Factor(i, j));
        }
        double[] b = Noise(n, seed);
        Assert.Equal(dense.Solve(b), band.Solve(b));
        Assert.Equal(dense.InverseQuadraticForm(b), band.InverseQuadraticForm(b));
        var drawDense = new double[n];
        var drawBand = new double[n];
        dense.SampleInnovation(b, drawDense);
        band.SampleInnovation(b, drawBand);
        Assert.Equal(drawDense, drawBand);
        Assert.Equal(dense.EquilibratedCondition(), band.EquilibratedCondition());
    }

    private static (DenseCholesky Dense, WideBandCholesky Band) Both(double[,] a, int bandwidth)
    {
        int n = a.GetLength(0);
        var dense = new DenseCholesky(n);
        dense.Decompose(a);
        var band = new WideBandCholesky(n, bandwidth);
        band.Decompose(a);
        return (dense, band);
    }

    [Fact]
    public void OnABlockBandedGram_EverythingIsTheDenseFactors_Exactly()
    {
        (double[,] gram, int bandwidth) = TensorGram();
        Assert.Equal(72, gram.GetLength(0));
        Assert.Equal(9, bandwidth);   // p·ν₂ + p at p = 1, ν₂ = 8
        (DenseCholesky dense, WideBandCholesky band) = Both(gram, bandwidth);
        Assert.True(band.Status.PositiveDefinite);
        AssertSameFactor(dense, band, seed: 1);

        // Against the matrix itself: L·Lᵀ reconstructs it, and A·(A⁻¹b) returns b.
        int n = gram.GetLength(0);
        for (int i = 0; i < n; i++)
            for (int j = 0; j <= i; j++)
            {
                double sum = 0.0;
                for (int k = 0; k <= j; k++) sum += band.Factor(i, k) * band.Factor(j, k);
                Assert.True(Math.Abs(sum - gram[i, j]) <= 1e-12 * (1.0 + Math.Abs(gram[i, j])), $"LLᵀ[{i},{j}]");
            }
        double[] b = Noise(n, 2);
        double[] x = band.Solve(b);
        double worst = 0.0;
        for (int i = 0; i < n; i++)
        {
            double ax = 0.0;
            for (int j = 0; j < n; j++) ax += gram[i, j] * x[j];
            worst = Math.Max(worst, Math.Abs(ax - b[i]));
        }
        Assert.True(worst <= 1e-10, $"residual {worst}");
    }

    [Theory]
    [InlineData(0)]
    [InlineData(5)]
    [InlineData(1000)]
    public void ADeclaredBandwidthAboveTheGrams_ChangesNothing_AndAFullBandIsTheDenseFactor(int extra)
    {
        (double[,] gram, int bandwidth) = TensorGram(seed: 1402);
        (DenseCholesky dense, WideBandCholesky band) = Both(gram, bandwidth + extra);   // 1000 clamps to n − 1: the full band
        Assert.Equal(Math.Min(bandwidth + extra, gram.GetLength(0) - 1), band.Bandwidth);
        AssertSameFactor(dense, band, seed: 3);
    }

    [Fact]
    public void ADiagonalMatrix_IsABandOfZero()
    {
        var a = new double[20, 20];
        for (int i = 0; i < 20; i++) a[i, i] = 1.0 + i;
        (DenseCholesky dense, WideBandCholesky band) = Both(a, 0);
        AssertSameFactor(dense, band, seed: 4);
        Assert.Equal(Enumerable.Range(0, 20).Sum(i => Math.Log(1.0 + i)), band.LogDet, 12);
    }

    [Fact]
    public void AnEntryBeyondTheBand_IsNoPartOfTheMatrixTheFactorSees()
    {
        (double[,] gram, int bandwidth) = TensorGram(seed: 1403);
        int n = gram.GetLength(0);
        var corrupted = (double[,])gram.Clone();
        corrupted[n - 1, 0] = corrupted[0, n - 1] = 0.25 * Math.Sqrt(gram[0, 0] * gram[n - 1, n - 1]);   // outside the band
        var dense = new DenseCholesky(n);
        dense.Decompose(corrupted);
        var band = new WideBandCholesky(n, bandwidth);
        band.Decompose(corrupted);
        var clean = new WideBandCholesky(n, bandwidth);
        clean.Decompose(gram);
        Assert.Equal(clean.LogDet, band.LogDet);        // the band reads what lies inside it
        Assert.NotEqual(dense.LogDet, band.LogDet);     // the dense factor reads the whole matrix: the band is the caller's declaration
    }

    [Fact]
    public void ANearlyDependentPair_ReportsTheDenseFactorsTinyPivot_AnExactOneItsFailure()
    {
        const int n = 80;
        BandedMatrix Make(double neighbour)
        {
            // A pentadiagonal SPD matrix whose rows 40 and 41 are as alike as `neighbour` makes them.
            var a = new double[n, n];
            for (int i = 0; i < n; i++) a[i, i] = 3.0;
            for (int i = 0; i + 1 < n; i++) a[i, i + 1] = a[i + 1, i] = 0.7;
            for (int i = 0; i + 2 < n; i++) a[i, i + 2] = a[i + 2, i] = 0.2;
            a[40, 40] = a[41, 41] = 1.0;
            a[40, 41] = a[41, 40] = neighbour;
            for (int j = 38; j <= 39; j++) a[40, j] = a[j, 40] = a[41, j] = a[j, 41] = 0.0;
            for (int j = 42; j <= 43; j++) a[40, j] = a[j, 40] = a[41, j] = a[j, 41] = 0.0;
            return new BandedMatrix(a, 2);
        }

        // Nearly dependent: the relative pivot of row 41 is 1 − neighbour² ≈ 2·10⁻⁸, and every figure is the dense factor's.
        BandedMatrix near = Make(1.0 - 1e-8);
        (DenseCholesky dense, WideBandCholesky band) = Both(near.A, 2);
        Assert.True(band.Status.PositiveDefinite);
        Assert.True(band.Status.MinRelativePivot < 3e-8);
        AssertSameFactor(dense, band, seed: 5);
        Assert.True(band.Status.PivotRatio < 1e-7);   // the smallest pivot over the largest: 1/κ at most

        // Exactly dependent, and indefinite: the dense factor's row, its status, no solve.
        foreach (double neighbour in new[] { 1.0, 1.5 })
        {
            BandedMatrix singular = Make(neighbour);
            (DenseCholesky d, WideBandCholesky w) = Both(singular.A, 2);
            Assert.Equal(41, w.FailedAt);
            Assert.Equal(d.FailedAt, w.FailedAt);
            Assert.Equal(d.Status, w.Status);
            Assert.False(w.Status.PositiveDefinite);
            Assert.True(double.IsNaN(w.LogDet));
            Assert.Throws<InvalidOperationException>(() => w.Solve(new double[n]));
            Assert.Throws<InvalidOperationException>(() => w.SampleInnovation(new double[n], new double[n]));
            Assert.Throws<InvalidOperationException>(() => w.EquilibratedCondition());
            Assert.Throws<InvalidOperationException>(() => w.Factor(0, 0));
        }
    }

    private sealed record BandedMatrix(double[,] A, int Band);

    [Fact]
    public void TheFactorRefusesWhatItCannotHold()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => new WideBandCholesky(0, 1));
        Assert.Throws<ArgumentOutOfRangeException>(() => new WideBandCholesky(5, -1));
        var band = new WideBandCholesky(5, 1);
        Assert.Throws<InvalidOperationException>(() => band.Solve(new double[5]));   // nothing factored
        Assert.Throws<ArgumentException>(() => band.Decompose(new double[4, 4]));
    }

    [Fact]
    public void TheSelectionRule_ChoosesTheBandWhenItSavesBothStorageAndArithmetic()
    {
        Assert.False(SpdFactors.UsesWideBand(63, 0));    // below the minimum order
        Assert.True(SpdFactors.UsesWideBand(64, 0));
        Assert.True(SpdFactors.UsesWideBand(64, 31));     // 2(p + 1) = n
        Assert.False(SpdFactors.UsesWideBand(64, 32));
        Assert.False(SpdFactors.UsesWideBand(1000, -1));
        Assert.True(SpdFactors.UsesWideBand(1000, 499));
        Assert.False(SpdFactors.UsesWideBand(1000, 500));

        (double[,] gram, int bandwidth) = TensorGram(seed: 1404);
        ISpdFactor chosen = SpdFactors.ForBand(gram, bandwidth);
        Assert.IsType<WideBandCholesky>(chosen);
        ISpdFactor dense = SpdFactors.ForBand(gram, gram.GetLength(0) / 2 + 1);   // a band too wide to save anything
        Assert.IsType<DenseCholesky>(dense);
        Assert.Equal(dense.LogDet, chosen.LogDet);
        Assert.Equal(dense.Status, chosen.Status);

        var small = new double[10, 10];
        for (int i = 0; i < 10; i++) small[i, i] = 2.0;
        Assert.IsType<DenseCholesky>(SpdFactors.ForBand(small, 0));       // below the minimum order
        Assert.Throws<ArgumentException>(() => SpdFactors.ForBand(new double[3, 4], 1));
    }

    [Fact]
    public void ABandCostsLessThanTheDenseFactor_OnAWideBandedGram()
    {
        // Not an assertion on a clock, a figure for the selection rule: n = 1600, half-bandwidth 80.
        const int n = 1600, p = 80;
        var rng = new Xoshiro256PlusPlus(1405);
        var a = new double[n, n];
        for (int i = 0; i < n; i++)
        {
            a[i, i] = 2.0 * p;
            for (int j = Math.Max(0, i - p); j < i; j++) a[i, j] = a[j, i] = rng.NextDouble() - 0.5;
        }
        var dense = new DenseCholesky(n);
        var band = new WideBandCholesky(n, p);
        var clock = Stopwatch.StartNew();
        dense.Decompose(a);
        double denseSeconds = clock.Elapsed.TotalSeconds;
        clock.Restart();
        band.Decompose(a);
        double bandSeconds = clock.Elapsed.TotalSeconds;
        _out.WriteLine($"n = {n}, p = {p}: dense {denseSeconds:F3} s, wide-band {bandSeconds:F3} s");
        Assert.Equal(dense.LogDet, band.LogDet);
        Assert.Equal(dense.Solve(Noise(n, 6)), band.Solve(Noise(n, 6)));
    }
}
