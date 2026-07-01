using System;
using Maths.Distance;
using Xunit;

namespace Maths.Distance.Tests;

public sealed class CosineTests
{
    [Fact]
    public void Similarity_Identical_Orthogonal_Opposite()
    {
        double[] u = { 1.0, 2.0, 3.0 };
        Assert.InRange(Cosine.Similarity(u, u), 1.0 - 1e-12, 1.0 + 1e-12);

        double[] a = { 1.0, 0.0 };
        double[] b = { 0.0, 1.0 };
        Assert.InRange(Cosine.Similarity(a, b), -1e-12, 1e-12);

        double[] c = { -1.0, 0.0 };
        Assert.InRange(Cosine.Similarity(a, c), -1.0 - 1e-12, -1.0 + 1e-12);
    }

    [Fact]
    public void Distance_IsOneMinusSimilarity()
    {
        double[] u = { 1.0, 2.0, 3.0 };
        double[] v = { 2.0, 1.0, 0.5 };
        Assert.Equal(1.0 - Cosine.Similarity(u, v), Cosine.Distance(u, v), 12);
    }

    /// <summary>
    /// All embeddings share a dominant common-mode direction. After ablating it, every
    /// embedding must be (near) orthogonal to that direction — exercises the SIMD
    /// Dot/ScaledSubtract path that backs IsotropizeInPlace.
    /// </summary>
    [Fact]
    public void IsotropizeInPlace_RemovesDominantDirection()
    {
        double[] u = Normalize(new[] { 1.0, 1.0, 0.0, 0.0 });
        var rng = new Random(5);
        int n = 16, d = 4;
        double[][] emb = new double[n][];
        for (int i = 0; i < n; i++)
        {
            double load = 10.0 + i; // large shared loading onto u
            emb[i] = new double[d];
            for (int j = 0; j < d; j++)
                emb[i][j] = load * u[j] + 0.01 * (rng.NextDouble() - 0.5);
        }

        Cosine.IsotropizeInPlace(emb, new[] { u });

        for (int i = 0; i < n; i++)
        {
            double proj = 0.0;
            for (int j = 0; j < d; j++) proj += emb[i][j] * u[j];
            Assert.InRange(proj, -1e-9, 1e-9);
        }
    }

    private static double[] Normalize(double[] v)
    {
        double norm = 0.0;
        foreach (double x in v) norm += x * x;
        norm = Math.Sqrt(norm);
        for (int i = 0; i < v.Length; i++) v[i] /= norm;
        return v;
    }
}
