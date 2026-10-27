using System;
using Maths.LinAlg;
using Maths.Rng;
using Xunit;

namespace Maths.LinAlg.Tests;

/// <summary>
/// The Khatri–Rao contractions against the materialized design: column j = a_j ⊗ h_j over n·m rows, case outer.
/// The Hadamard Gram, its centred three-term form, the borders a birth needs, the response cross-products (centred
/// ones on a response with a large level), and the scattered accumulation over the full grid and over a subset. The
/// terms include a design-only term (h = 1), a functional-only term (a = 1), hinges with zeros, and factors whose
/// means are large against their spread.
/// </summary>
public sealed class KhatriRaoTests
{
    [Fact]
    public void Gram_AndBorder_MatchTheMaterializedDesign()
    {
        (double[][] a, double[][] h) = Factors(7, 9, new Xoshiro256PlusPlus(seed: 1101), intercept: true);
        double[,] x = Materialize(a, h);
        AssertMatrix(Gram(x), KhatriRao.Gram(a, h), 1e-14);

        double[,] g = KhatriRao.Gram(a, h);
        int last = a.Length - 1;
        double[] u = KhatriRao.Border(a[..last], h[..last], a[last], h[last], out double diagonal);
        for (int k = 0; k < last; k++) Assert.Equal(g[last, k], u[k]);
        Assert.Equal(g[last, last], diagonal);
    }

    [Fact]
    public void CentredGram_AndBorder_MatchExplicitCentring()
    {
        (double[][] a, double[][] h) = Factors(7, 9, new Xoshiro256PlusPlus(seed: 1102), intercept: false);
        double[,] xc = CentreColumns(Materialize(a, h));
        double[,] gc = KhatriRao.CentredGram(a, h);
        AssertMatrix(Gram(xc), gc, 1e-14);

        int last = a.Length - 1;
        double[] u = KhatriRao.CentredBorder(a[..last], h[..last], a[last], h[last], out double diagonal);
        for (int k = 0; k < last; k++) Assert.Equal(gc[last, k], u[k]);
        Assert.Equal(gc[last, last], diagonal);
    }

    [Fact]
    public void CrossProducts_MatchTheMaterializedDesign()
    {
        const int n = 7, m = 9;
        var rng = new Xoshiro256PlusPlus(seed: 1103);
        (double[][] a, double[][] h) = Factors(n, m, rng, intercept: false);
        var y = new double[n, m];
        for (int i = 0; i < n; i++) for (int l = 0; l < m; l++) y[i, l] = 1e3 + rng.NextDouble() - 0.5;   // a large level
        var vec = new double[n * m];
        for (int i = 0; i < n; i++) for (int l = 0; l < m; l++) vec[i * m + l] = y[i, l];

        double[,] x = Materialize(a, h);
        AssertVector(Transpose(x, vec), KhatriRao.CrossProducts(a, h, y), 1e-14);

        // Centred columns are orthogonal to the constant, so the reference reads the centred response.
        double[,] xc = CentreColumns(x);
        double level = PairwiseSum.Of(vec) / vec.Length;
        var centred = new double[vec.Length];
        for (int r = 0; r < vec.Length; r++) centred[r] = vec[r] - level;
        AssertVector(Transpose(xc, centred), KhatriRao.CentredCrossProducts(a, h, y), 1e-12);
    }

    [Fact]
    public void Scattered_OverTheFullGrid_IsTheCartesianShortcut()
    {
        const int n = 40, m = 50;                    // 2,000 samples: fifteen blocks through the cascade
        var rng = new Xoshiro256PlusPlus(seed: 1104);
        (double[][] a, double[][] h) = Factors(n, m, rng, intercept: true);
        var cases = new int[n * m];
        var locations = new int[n * m];
        var y = new double[n, m];
        var flat = new double[n * m];
        for (int i = 0; i < n; i++)
            for (int l = 0; l < m; l++)
            {
                cases[i * m + l] = i;
                locations[i * m + l] = l;
                y[i, l] = flat[i * m + l] = rng.NextDouble() - 0.5;
            }
        AssertMatrix(KhatriRao.Gram(a, h), KhatriRao.ScatteredGram(a, h, cases, locations), 1e-13);
        AssertVector(KhatriRao.CrossProducts(a, h, y), KhatriRao.ScatteredCrossProducts(a, h, cases, locations, flat), 1e-13);
    }

    [Fact]
    public void Scattered_OverASubset_MatchesItsRows()
    {
        const int n = 7, m = 9, samples = 30;
        var rng = new Xoshiro256PlusPlus(seed: 1105);
        (double[][] a, double[][] h) = Factors(n, m, rng, intercept: true);
        double[,] x = Materialize(a, h);
        var cases = new int[samples];
        var locations = new int[samples];
        var y = new double[samples];
        var rows = new double[samples, a.Length];
        for (int r = 0; r < samples; r++)
        {
            cases[r] = rng.NextInt(n);
            locations[r] = rng.NextInt(m);
            y[r] = rng.NextDouble() - 0.5;
            for (int j = 0; j < a.Length; j++) rows[r, j] = x[cases[r] * m + locations[r], j];
        }
        AssertMatrix(Gram(rows), KhatriRao.ScatteredGram(a, h, cases, locations), 1e-14);
        AssertVector(Transpose(rows, y), KhatriRao.ScatteredCrossProducts(a, h, cases, locations, y), 1e-14);
    }

    // Five terms: two products of factors near 1 (means large against spread), a hinge product with zeros, a
    // design-only and a functional-only term; R's intercept column first when asked.
    private static (double[][] A, double[][] H) Factors(int n, int m, Xoshiro256PlusPlus rng, bool intercept)
    {
        int count = intercept ? 6 : 5, j = 0;
        var a = new double[count][];
        var h = new double[count][];
        double[] Fill(int length, Func<double> f)
        {
            var v = new double[length];
            for (int i = 0; i < length; i++) v[i] = f();
            return v;
        }
        if (intercept)
        {
            a[j] = Fill(n, () => 1.0);
            h[j++] = Fill(m, () => 1.0);
        }
        for (int t = 0; t < 2; t++)
        {
            a[j] = Fill(n, () => 0.9 + 0.1 * rng.NextDouble());
            h[j++] = Fill(m, () => 0.95 + 0.05 * rng.NextDouble());
        }
        a[j] = Fill(n, () => Math.Max(rng.NextDouble() - 0.5, 0.0));
        h[j++] = Fill(m, () => Math.Max(rng.NextDouble() - 0.3, 0.0));
        a[j] = Fill(n, () => rng.NextDouble());
        h[j++] = Fill(m, () => 1.0);
        a[j] = Fill(n, () => 1.0);
        h[j] = Fill(m, () => rng.NextDouble());
        return (a, h);
    }

    private static double[,] Materialize(double[][] a, double[][] h)
    {
        int n = a[0].Length, m = h[0].Length;
        var x = new double[n * m, a.Length];
        for (int j = 0; j < a.Length; j++)
            for (int i = 0; i < n; i++)
                for (int l = 0; l < m; l++) x[i * m + l, j] = a[j][i] * h[j][l];
        return x;
    }

    private static double[,] CentreColumns(double[,] x)
    {
        int rows = x.GetLength(0), cols = x.GetLength(1);
        var c = (double[,])x.Clone();
        var column = new double[rows];
        for (int j = 0; j < cols; j++)
            for (int pass = 0; pass < 2; pass++)
            {
                for (int r = 0; r < rows; r++) column[r] = c[r, j];
                double mean = PairwiseSum.Of(column) / rows;
                for (int r = 0; r < rows; r++) c[r, j] -= mean;
            }
        return c;
    }

    private static double[,] Gram(double[,] x)
    {
        int rows = x.GetLength(0), cols = x.GetLength(1);
        var g = new double[cols, cols];
        var ci = new double[rows];
        var cj = new double[rows];
        for (int i = 0; i < cols; i++)
            for (int j = 0; j <= i; j++)
            {
                for (int r = 0; r < rows; r++)
                {
                    ci[r] = x[r, i];
                    cj[r] = x[r, j];
                }
                g[i, j] = g[j, i] = PairwiseSum.Dot(ci, cj);
            }
        return g;
    }

    private static double[] Transpose(double[,] x, double[] y)
    {
        int rows = x.GetLength(0), cols = x.GetLength(1);
        var result = new double[cols];
        var column = new double[rows];
        for (int j = 0; j < cols; j++)
        {
            for (int r = 0; r < rows; r++) column[r] = x[r, j];
            result[j] = PairwiseSum.Dot(column, y);
        }
        return result;
    }

    private static void AssertMatrix(double[,] expected, double[,] actual, double tolerance)
    {
        double scale = 0.0;
        foreach (double v in expected) scale = Math.Max(scale, Math.Abs(v));
        for (int i = 0; i < expected.GetLength(0); i++)
            for (int j = 0; j < expected.GetLength(1); j++)
                Assert.True(Math.Abs(expected[i, j] - actual[i, j]) <= tolerance * scale, $"[{i},{j}] {actual[i, j]:R} vs {expected[i, j]:R}");
    }

    private static void AssertVector(double[] expected, double[] actual, double tolerance)
    {
        double scale = 0.0;
        foreach (double v in expected) scale = Math.Max(scale, Math.Abs(v));
        Assert.Equal(expected.Length, actual.Length);
        for (int i = 0; i < expected.Length; i++)
            Assert.True(Math.Abs(expected[i] - actual[i]) <= tolerance * scale, $"[{i}] {actual[i]:R} vs {expected[i]:R}");
    }
}
