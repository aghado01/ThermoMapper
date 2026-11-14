using System;
using Maths.Rng;

namespace Clustering.Geometric.EMeans.Tests;

/// <summary>Planted Gaussian blobs: the ground truth is the blob each point was drawn from.</summary>
internal static class Blobs
{
    /// <summary>Blob centres on hypercube corners of half-side <paramref name="halfSide"/>: blob <c>g</c> takes the
    /// corner whose signs are the bits of <c>g</c>, so two centres are at least <c>2·halfSide</c> apart.</summary>
    public static double[][] Corners(int count, int dimension, double halfSide)
    {
        if (count > 1 << Math.Min(dimension, 30)) throw new ArgumentException("More blobs than corners.");
        var centres = new double[count][];
        for (int g = 0; g < count; g++)
        {
            centres[g] = new double[dimension];
            for (int c = 0; c < dimension; c++)
                centres[g][c] = c < 31 && ((g >> c) & 1) == 1 ? halfSide : -halfSide;
        }
        return centres;
    }

    /// <summary>Draws <paramref name="perBlob"/> points around each centre, blob by blob.</summary>
    public static (double[][] Data, int[] Labels) Draw(double[][] centres, int perBlob, double sigma, int seed)
    {
        var rng = new Xoshiro256PlusPlus(seed);
        int d = centres[0].Length;
        var data = new double[centres.Length * perBlob][];
        var labels = new int[data.Length];
        for (int g = 0, i = 0; g < centres.Length; g++)
        {
            for (int p = 0; p < perBlob; p++, i++)
            {
                var x = new double[d];
                for (int c = 0; c < d; c++) x[c] = centres[g][c] + sigma * Normal(rng);
                data[i] = x;
                labels[i] = g;
            }
        }
        return (data, labels);
    }

    /// <summary>True when the two labelings induce the same partition, whatever the label names.</summary>
    public static bool SamePartition(int[] a, int[] b)
    {
        if (a.Length != b.Length) return false;
        var forward = new System.Collections.Generic.Dictionary<int, int>();
        var backward = new System.Collections.Generic.Dictionary<int, int>();
        for (int i = 0; i < a.Length; i++)
        {
            if (forward.TryGetValue(a[i], out int mapped))
            {
                if (mapped != b[i]) return false;
            }
            else
            {
                if (backward.ContainsKey(b[i])) return false;
                forward[a[i]] = b[i];
                backward[b[i]] = a[i];
            }
        }
        return true;
    }

    private static double Normal(Xoshiro256PlusPlus rng)
    {
        double u = 1.0 - rng.NextDouble();
        double v = rng.NextDouble();
        return Math.Sqrt(-2.0 * Math.Log(u)) * Math.Cos(2.0 * Math.PI * v);
    }
}
