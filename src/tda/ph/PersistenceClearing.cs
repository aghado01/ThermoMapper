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

    public static (List<Bar> Bars, HashSet<int> ToReduce, HashSet<int> ToSkip) ComputeH0(
        IFiltration filtration, int maxDimension, double cutoff)
    {
        var bars = new List<Bar>();
        var toReduce = new HashSet<int>();
        var toSkip = new HashSet<int>();
        int n = filtration.Count;

        var vertexIndex = new Dictionary<int, int>();
        var birth = new double[n];
        var parent = new int[n];
        var birthId = new int[n];
        int nextId = 0;

        for (int i = 0; i < n; i++)
        {
            parent[i] = i;
            if (filtration.GetDimension(i) == 0)
            {
                int v = filtration.GetVertices(i)[0];
                vertexIndex[v] = i;
                birth[i] = filtration.GetBirth(i);
                birthId[i] = nextId++;
            }
        }

        var edges = new List<int>();
        for (int i = 0; i < n; i++)
            if (filtration.GetDimension(i) == 1)
                edges.Add(i);

        edges.Sort((a, b) => filtration.GetBirth(a).CompareTo(filtration.GetBirth(b)));

        int Find(int x)
        {
            int r = x;
            while (parent[r] != r) r = parent[r];
            while (parent[x] != r) { int p = parent[x]; parent[x] = r; x = p; }
            return r;
        }

        foreach (int e in edges)
        {
            ReadOnlySpan<int> verts = filtration.GetVertices(e);
            int u = verts[0];
            int v = verts[1];
            int i = vertexIndex[u];
            int j = vertexIndex[v];
            int ri = Find(i);
            int rj = Find(j);
            double t = filtration.GetBirth(e);

            if (ri != rj)
            {
                int elder = birth[ri] > birth[rj] || (birth[ri] == birth[rj] && birthId[ri] > birthId[rj]) ? ri : rj;
                int younger = elder == ri ? rj : ri;
                if (maxDimension >= 0 && PassesCutoff(birth[younger], t, cutoff))
                    bars.Add(new Bar(birth[younger], t, 0));

                parent[younger] = elder;
                toSkip.Add(e);
            }
            else
            {
                toReduce.Add(e);
            }
        }

        for (int i = 0; i < n; i++)
        {
            if (filtration.GetDimension(i) != 0) continue;
            if (Find(i) == i && maxDimension >= 0)
                bars.Add(new Bar(birth[i], double.PositiveInfinity, 0));
        }

        return (bars, toReduce, toSkip);
    }
}
