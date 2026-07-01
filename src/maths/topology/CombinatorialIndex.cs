#nullable enable
using System;

namespace Maths.Topology;

/// <summary>
/// Combinatorial number system for k-subsets of {0,…,n−1} in colex order:
/// index = Σᵢ C(vᵢ, i+1) for vertices sorted ascending.
/// </summary>
public static class CombinatorialIndex
{
    public static int Index(ReadOnlySpan<int> verticesAscending)
    {
        int acc = 0;
        for (int i = 0; i < verticesAscending.Length; i++)
            acc += Binomial(verticesAscending[i], i + 1);
        return acc;
    }

    public static int[] Vertices(int index, int dimension)
    {
        int k = dimension + 1;
        var verts = new int[k];
        int remaining = index;
        for (int i = k - 1; i >= 0; i--)
        {
            verts[i] = i;
            while (Binomial(verts[i] + 1, i + 1) <= remaining)
                verts[i]++;
            remaining -= Binomial(verts[i], i + 1);
        }

        return verts;
    }

    public static long PackKey(int dimension, ReadOnlySpan<int> verticesAscending) =>
        PackKey(dimension, Index(verticesAscending));

    public static long PackKey(int dimension, int combinatorialIndex) =>
        ((long)dimension << 32) | (uint)combinatorialIndex;

    public static int Binomial(int n, int k)
    {
        if (k < 0 || n < k)
            return 0;
        if (k == 0)
            return 1;

        long result = 1;
        for (int i = 0; i < k; i++)
            result = result * (n - i) / (i + 1);
        return (int)result;
    }
}
