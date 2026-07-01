#nullable enable
using System;
using System.Collections.Generic;
using System.Linq;

using Maths.Topology;
namespace TDA.Ph;

/// <summary>
/// Ordered simplicial filtration with combinatorial face index and cofacet lists.
/// </summary>
public sealed class SimplicialFiltration : IFiltration
{
    public IReadOnlyList<Simplex> Simplices { get; }
    public string Label { get; }
    public bool EmergentPairs => false;

    public int Count => Simplices.Count;

    private readonly Dictionary<long, int> _index;
    private readonly int[][] _cofacets;

    public SimplicialFiltration(IEnumerable<Simplex> simplices, string label = "")
    {
        ArgumentNullException.ThrowIfNull(simplices);

        var sorted = simplices
            .OrderBy(s => s.FiltrationValue)
            .ThenBy(s => s.Dimension)
            .ThenBy(s => CombinatorialIndex.Index(s.Vertices))
            .ToList();

        Simplices = sorted;
        Label = label;

        _index = new Dictionary<long, int>(sorted.Count);
        for (int i = 0; i < sorted.Count; i++)
        {
            Simplex s = sorted[i];
            long key = CombinatorialIndex.PackKey(s.Dimension, s.Vertices);
            _index[key] = i;
        }

        _cofacets = new int[sorted.Count][];
        var builders = new List<int>[sorted.Count];
        for (int i = 0; i < sorted.Count; i++)
            builders[i] = new List<int>();

        for (int j = 0; j < sorted.Count; j++)
        {
            int[] boundary = GetBoundaryIndicesInternal(j);
            for (int b = 0; b < boundary.Length; b++)
                builders[boundary[b]].Add(j);
        }

        for (int i = 0; i < sorted.Count; i++)
        {
            List<int> cof = builders[i];
            if (cof.Count == 0)
            {
                _cofacets[i] = Array.Empty<int>();
                continue;

            }

            cof.Sort();
            cof.Reverse();
            _cofacets[i] = cof.ToArray();
        }

        for (int j = 0; j < sorted.Count; j++)
        {
            var v = sorted[j].Vertices;
            if (v.Length <= 1) continue;
            for (int k = 0; k < v.Length; k++)
            {
                long faceKey = CombinatorialIndex.PackKey(v.Length - 2, RemoveAt(v, k));
                if (!_index.ContainsKey(faceKey))
                    throw new ArgumentException(
                        $"Simplex [{FormatVerts(v)}] at index {j} is missing face [{FormatVerts(RemoveAt(v, k))}].");
                int fi = _index[faceKey];
                if (fi > j)
                    throw new ArgumentException(
                        $"Face at index {fi} appears after coface [{FormatVerts(v)}] (index {j}).");
            }
        }
    }

    public int GetDimension(int simplexIndex) => Simplices[simplexIndex].Dimension;

    public double GetBirth(int simplexIndex) => Simplices[simplexIndex].FiltrationValue;

    public ReadOnlySpan<int> GetVertices(int simplexIndex) => Simplices[simplexIndex].Vertices;

    public int[] GetBoundaryIndices(int simplexIndex) => GetBoundaryIndicesInternal(simplexIndex);

    public int[] GetCoboundaryIndices(int simplexIndex) => _cofacets[simplexIndex];

    int[] GetBoundaryIndicesInternal(int j)
    {
        var v = Simplices[j].Vertices;
        if (v.Length <= 1) return Array.Empty<int>();
        var result = new int[v.Length];
        for (int k = 0; k < v.Length; k++)
            result[k] = _index[CombinatorialIndex.PackKey(v.Length - 2, RemoveAt(v, k))];
        return result;
    }

    static int[] RemoveAt(int[] v, int k)
    {
        var face = new int[v.Length - 1];
        for (int i = 0, w = 0; i < v.Length; i++)
            if (i != k) face[w++] = v[i];
        return face;
    }

    static string FormatVerts(int[] v) => string.Join(",", v);
}
