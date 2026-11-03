#nullable enable
using System;
using System.Collections.Generic;

using Maths.Topology;
namespace TDA.Ph;

/// <summary>Shared H0 union-find clearing and cutoff helpers for Ripser-style engines.</summary>
internal static class PersistenceClearing
{
    public static bool PassesCutoff(double birth, double death, double cutoff)
    {
        double p = Persistence(birth, death);
        return double.IsPositiveInfinity(p) || p >= cutoff;
    }

    public static double Persistence(double birth, double death) =>
        double.IsPositiveInfinity(death) ? double.PositiveInfinity : death - birth;

    /// <summary>
    /// H₀ by union-find over the 1-cells in filtration order, which is the pairing a full reduction finds. A 1-cell reads
    /// through its faces, so any cell complex works: two faces join, one joins its vertex to a root older than every cell
    /// (the far end lies in a quotiented subcomplex), and none closes a cycle by itself. With <paramref name="pairs"/>, also
    /// the pairs as cell indices: a vertex and the 1-cell that killed it, or −1.
    /// </summary>
    /// <returns>The bars, and the 1-cells that killed nothing, whose degree is the next to reduce.</returns>
    public static (List<Bar> Bars, HashSet<int> ToReduce) ComputeH0(
        IFiltration filtration, int maxDimension, double cutoff, List<(int Birth, int Death)>? pairs = null)
    {
        var bars = new List<Bar>();
        var toReduce = new HashSet<int>();
        int n = filtration.Count;
        int quotient = n;

        var parent = new int[n + 1];
        for (int i = 0; i <= n; i++)
            parent[i] = i;

        int Find(int x)
        {
            int r = x;
            while (parent[r] != r) r = parent[r];
            while (parent[x] != r) { int p = parent[x]; parent[x] = r; x = p; }
            return r;
        }

        for (int e = 0; e < n; e++)
        {
            if (filtration.GetDimension(e) != 1) continue;

            int[] faces = filtration.GetBoundaryIndices(e);
            if (faces.Length > 2)
            {
                pairs?.Clear();
                return ReduceH0Columns(filtration, maxDimension, cutoff, pairs);
            }

            if (faces.Length == 0)
            {
                toReduce.Add(e);
                continue;
            }

            int ri = Find(faces[0]);
            int rj = faces.Length == 2 ? Find(faces[1]) : quotient;
            if (ri == rj)
            {
                toReduce.Add(e);
                continue;
            }

            // The elder is the earlier cell, as a reduction's pivot leaves it, and the quotient outlives them all.
            int younger = ri == quotient ? rj : rj == quotient ? ri : Math.Max(ri, rj);
            int elder = younger == ri ? rj : ri;
            AddFinite(bars, pairs, filtration, younger, e, maxDimension, cutoff);
            parent[younger] = elder;
        }

        for (int i = 0; i < n; i++)
        {
            if (filtration.GetDimension(i) != 0) continue;
            if (Find(i) == i && maxDimension >= 0)
            {
                bars.Add(new Bar(filtration.GetBirth(i), double.PositiveInfinity, 0));
                pairs?.Add((i, -1));
            }
        }

        return (bars, toReduce);
    }

    /// <summary>
    /// H₀ for a general chain complex, where a 1-cell can have more than two faces and no union-find reads it: the columns of
    /// ∂₁ reduced. The result is the same on a complex the union-find can read.
    /// </summary>
    static (List<Bar> Bars, HashSet<int> ToReduce) ReduceH0Columns(
        IFiltration filtration, int maxDimension, double cutoff, List<(int Birth, int Death)>? pairs)
    {
        var bars = new List<Bar>();
        var toReduce = new HashSet<int>();
        var reduced = new Dictionary<int, SortedSet<int>>();
        int n = filtration.Count;

        for (int e = 0; e < n; e++)
        {
            if (filtration.GetDimension(e) != 1) continue;

            var col = new SortedSet<int>(filtration.GetBoundaryIndices(e));
            while (col.Count > 0 && reduced.TryGetValue(col.Max, out SortedSet<int>? owner))
                col.SymmetricExceptWith(owner);

            if (col.Count == 0)
            {
                toReduce.Add(e);
                continue;
            }

            reduced[col.Max] = col;
            AddFinite(bars, pairs, filtration, col.Max, e, maxDimension, cutoff);
        }

        for (int i = 0; i < n; i++)
        {
            if (filtration.GetDimension(i) != 0 || reduced.ContainsKey(i) || maxDimension < 0) continue;
            bars.Add(new Bar(filtration.GetBirth(i), double.PositiveInfinity, 0));
            pairs?.Add((i, -1));
        }

        return (bars, toReduce);
    }

    static void AddFinite(
        List<Bar> bars, List<(int Birth, int Death)>? pairs, IFiltration filtration, int vertex, int killer, int maxDimension, double cutoff)
    {
        double birth = filtration.GetBirth(vertex);
        double death = filtration.GetBirth(killer);
        if (maxDimension >= 0 && PassesCutoff(birth, death, cutoff))
        {
            bars.Add(new Bar(birth, death, 0));
            pairs?.Add((vertex, killer));
        }
    }
}
