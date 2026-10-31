using System;
using System.Collections.Generic;
using System.Linq;
using Maths.LinAlg;
using Maths.Rng;
using Xunit;
using Xunit.Abstractions;

namespace Maths.LinAlg.Tests;

/// <summary>
/// The face-splitting contraction (N-X1): the Gram and cross products of a tensor design X = B₁ ∗ B₂ ∗ … ∗ B_d, whose
/// row i is the row-wise Kronecker product of the axes' basis rows, against the explicit product on small designs.
/// <list type="bullet">
/// <item>Scattered samples, two and three axes, piecewise-linear hat bases (two nonzeros a row per axis), against XᵀX and Xᵀy
/// of the materialized product; the half-bandwidth of the Gram is p·ν₂ + p at p = 1, with the cells all hit.</item>
/// <item>Gridded samples, where the Gram is the Kronecker product of the axes' Grams, against the materialized product and
/// against the scattered accumulation over the same grid.</item>
/// <item>The scattered Khatri–Rao accumulation the new one shares its blocks with is unchanged, bit for bit.</item>
/// </list>
/// </summary>
public sealed class FaceSplittingTests
{
    private readonly ITestOutputHelper _out;
    public FaceSplittingTests(ITestOutputHelper output) => _out = output;

    // A piecewise-linear hat basis on `width` equally spaced nodes over [0, 1], evaluated at x: two adjacent nonzeros.
    private static double[] Hat(double x, int width)
    {
        var row = new double[width];
        double t = x * (width - 1);
        int j = Math.Min((int)t, width - 2);
        row[j] = 1.0 - (t - j);
        row[j + 1] = t - j;
        return row;
    }

    private static double[,] Basis(double[] x, int width)
    {
        var b = new double[x.Length, width];
        for (int i = 0; i < x.Length; i++)
        {
            double[] row = Hat(x[i], width);
            for (int j = 0; j < width; j++) b[i, j] = row[j];
        }
        return b;
    }

    // The explicit face-splitting product: row i is b¹ᵢ ⊗ b²ᵢ ⊗ …, columns lexicographic with the first axis slowest.
    private static double[,] Explicit(IReadOnlyList<double[,]> axes)
    {
        int rows = axes[0].GetLength(0);
        int width = axes.Aggregate(1, (acc, b) => acc * b.GetLength(1));
        var x = new double[rows, width];
        for (int i = 0; i < rows; i++)
        {
            double[] row = [1.0];
            foreach (double[,] b in axes)
            {
                int w = b.GetLength(1);
                var next = new double[row.Length * w];
                for (int c = 0; c < row.Length; c++)
                    for (int j = 0; j < w; j++) next[c * w + j] = row[c] * b[i, j];
                row = next;
            }
            for (int c = 0; c < width; c++) x[i, c] = row[c];
        }
        return x;
    }

    private static double[,] Gram(double[,] x)
    {
        int rows = x.GetLength(0), width = x.GetLength(1);
        var g = new double[width, width];
        for (int i = 0; i < width; i++)
            for (int j = 0; j < width; j++)
            {
                double sum = 0.0;
                for (int r = 0; r < rows; r++) sum += x[r, i] * x[r, j];
                g[i, j] = sum;
            }
        return g;
    }

    private static void AssertMatrix(double[,] expected, double[,] actual, double tolerance)
    {
        Assert.Equal(expected.GetLength(0), actual.GetLength(0));
        for (int i = 0; i < expected.GetLength(0); i++)
            for (int j = 0; j < expected.GetLength(1); j++)
                Assert.True(Math.Abs(expected[i, j] - actual[i, j]) <= tolerance * (1.0 + Math.Abs(expected[i, j])), $"[{i},{j}]: {expected[i, j]} against {actual[i, j]}");
    }

    private static double[] Points(Xoshiro256PlusPlus rng, int n)
    {
        var x = new double[n];
        for (int i = 0; i < n; i++) x[i] = rng.NextDouble();
        return x;
    }

    [Theory]
    [InlineData(2, 700)]
    [InlineData(3, 300)]
    public void ScatteredGram_AndCrossProducts_MatchTheExplicitProduct(int axisCount, int samples)
    {
        var rng = new Xoshiro256PlusPlus(seed: 1301 + axisCount);
        int[] widths = axisCount == 2 ? [7, 6] : [4, 3, 5];
        double[][] x = [.. widths.Select(_ => Points(rng, samples))];
        IReadOnlyList<double[,]> axes = [.. widths.Select((w, a) => Basis(x[a], w))];
        double[,] explicitX = Explicit(axes);
        int width = explicitX.GetLength(1);

        AssertMatrix(Gram(explicitX), KhatriRao.FaceSplittingGram(axes), 1e-13);

        double[] y = [.. Enumerable.Range(0, samples).Select(_ => 100.0 + rng.NextDouble())];
        double[] cross = KhatriRao.FaceSplittingCrossProducts(axes, y);
        for (int c = 0; c < width; c++)
        {
            double expected = 0.0;
            for (int r = 0; r < samples; r++) expected += explicitX[r, c] * y[r];
            Assert.True(Math.Abs(cross[c] - expected) <= 1e-12 * (1.0 + Math.Abs(expected)), $"Xᵀy[{c}]: {cross[c]} against {expected}");
        }
    }

    [Fact]
    public void TheHalfBandwidth_IsTheSpanOfARowsNonzeros_WhichIsPNuTwoPlusP_ForALinearTensorBasis()
    {
        var rng = new Xoshiro256PlusPlus(seed: 1305);
        const int samples = 800, nu1 = 8, nu2 = 7;
        IReadOnlyList<double[,]> axes = [Basis(Points(rng, samples), nu1), Basis(Points(rng, samples), nu2)];
        int p = 1;
        Assert.Equal(p * nu2 + p, KhatriRao.FaceSplittingHalfBandwidth(axes));

        // It is the Gram's own: no entry beyond it is nonzero, and one at it is.
        double[,] g = KhatriRao.FaceSplittingGram(axes);
        int band = KhatriRao.FaceSplittingHalfBandwidth(axes), beyond = 0, atEdge = 0;
        for (int i = 0; i < g.GetLength(0); i++)
            for (int j = 0; j < i; j++)
            {
                if (i - j > band && g[i, j] != 0.0) beyond++;
                if (i - j == band && g[i, j] != 0.0) atEdge++;
            }
        Assert.Equal(0, beyond);
        Assert.True(atEdge > 0);
    }

    [Fact]
    public void KroneckerGram_OnAGrid_IsTheMaterializedGram_AndTheScatteredOneOverTheGrid()
    {
        var rng = new Xoshiro256PlusPlus(seed: 1306);
        double[] x1 = Points(rng, 9), x2 = Points(rng, 11);
        double[,] b1 = Basis(x1, 5), b2 = Basis(x2, 6);

        // The grid's rows, (i1, i2) with i1 outer: the Kronecker product of the bases.
        var gridded = new double[x1.Length * x2.Length, 5 * 6];
        var left = new double[x1.Length * x2.Length, 5];
        var right = new double[x1.Length * x2.Length, 6];
        for (int i1 = 0; i1 < x1.Length; i1++)
            for (int i2 = 0; i2 < x2.Length; i2++)
            {
                int row = i1 * x2.Length + i2;
                for (int j = 0; j < 5; j++) left[row, j] = b1[i1, j];
                for (int k = 0; k < 6; k++) right[row, k] = b2[i2, k];
                for (int j = 0; j < 5; j++)
                    for (int k = 0; k < 6; k++) gridded[row, j * 6 + k] = b1[i1, j] * b2[i2, k];
            }
        double[,] kronecker = KhatriRao.KroneckerGram([b1, b2]);
        AssertMatrix(Gram(gridded), kronecker, 1e-13);
        AssertMatrix(kronecker, KhatriRao.FaceSplittingGram([left, right]), 1e-13);
    }

    [Fact]
    public void TheContraction_ChecksItsAxes()
    {
        double[,] a = new double[4, 3], b = new double[5, 3];
        Assert.Throws<ArgumentException>(() => KhatriRao.FaceSplittingGram([]));
        Assert.Throws<ArgumentException>(() => KhatriRao.FaceSplittingGram([a, b]));   // rows of one design
        Assert.Throws<ArgumentException>(() => KhatriRao.FaceSplittingCrossProducts([a, a], new double[3]));
        Assert.Throws<ArgumentException>(() => KhatriRao.KroneckerGram([new double[4, 0]]));
        Assert.Equal(9, KhatriRao.KroneckerGram([a, b]).GetLength(0));   // grids need not share row counts
    }

    // The scattered accumulation (the Hadamard-of-rows one) the face-splitting shares its cascade with is not moved.
    [Fact]
    public void ScatteredGram_IsUnchangedByTheSharedAccumulation_BitForBit()
    {
        const int n = 40, m = 50, terms = 6;
        var rng = new Xoshiro256PlusPlus(seed: 1307);
        var a = new double[terms][];
        var h = new double[terms][];
        for (int j = 0; j < terms; j++)
        {
            a[j] = new double[n];
            h[j] = new double[m];
            for (int i = 0; i < n; i++) a[j][i] = j == 0 ? 1.0 : Math.Max(rng.NextDouble() - 0.4, 0.0);
            for (int l = 0; l < m; l++) h[j][l] = j == 1 ? 1.0 : Math.Max(rng.NextDouble() - 0.3, 0.0);
        }
        var cases = new int[n * m + 123];
        var locations = new int[cases.Length];
        for (int r = 0; r < cases.Length; r++) { cases[r] = rng.NextInt(n); locations[r] = rng.NextInt(m); }
        double[,] g = KhatriRao.ScatteredGram(a, h, cases, locations);
        ulong hash = 14695981039346656037UL;
        for (int i = 0; i < terms; i++)
            for (int j = 0; j < terms; j++) hash = (hash ^ (ulong)BitConverter.DoubleToInt64Bits(g[i, j])) * 1099511628211UL;
        _out.WriteLine($"ScatteredGram hash {hash:x16}");
        Assert.Equal(0x6e4c55f9a3d70b21UL, hash);   // taken from the accumulation before it was shared
    }
}
