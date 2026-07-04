#nullable enable
using System;
using System.Collections.Generic;
using System.Linq;
using Graphs.Primitives;
using Graphs.Spectral;
using Maths.Topology;

namespace TDA.Ph;

/// <summary>
/// Graph-backed Vietoris–Rips filtration without a materialized <see cref="SimplicialFiltration"/>.
/// Discovers triangles from CSR adjacency; simplex order matches global filtration sort.
/// </summary>
public sealed class LazyRipsFiltration : IFiltration
{
    public string Label { get; }
    public bool EmergentPairs => false;
    public int Count => _simplices.Count;

    private readonly List<Simplex> _simplices;
    private readonly Dictionary<long, int> _index;
    private int[][] _cofacets;

    public LazyRipsFiltration(CsrGraph graph, FiltrationWeights weights, string label = "LazyRips")
    {
        ArgumentNullException.ThrowIfNull(weights);

        Label = label;
        _edgeWeights = ResolveEdgeWeights(graph, weights);

        var discovered = new List<Simplex>();
        int n = graph.NodeCount;

        for (int i = 0; i < n; i++)
            discovered.Add(new Simplex(0.0, i));

        foreach (var ((lo, hi), w) in _edgeWeights)
            discovered.Add(new Simplex(w, lo, hi));

        var triangleKeys = new Dictionary<long, Simplex>();
        foreach (var ((lo, hi), edgeBirth) in _edgeWeights)
        {
            foreach (int w in CommonNeighbors(graph, lo, hi))
            {
                double wU = _edgeWeights.GetValueOrDefault(Key(lo, w), double.PositiveInfinity);
                double wV = _edgeWeights.GetValueOrDefault(Key(hi, w), double.PositiveInfinity);
                double birth = Math.Max(edgeBirth, Math.Max(wU, wV));
                var tri = new Simplex(birth, lo, hi, w);
                triangleKeys[CombinatorialIndex.PackKey(tri.Dimension, tri.Vertices)] = tri;
            }
        }

        discovered.AddRange(triangleKeys.Values);

        var sorted = discovered
            .OrderBy(s => s.FiltrationValue)
            .ThenBy(s => s.Dimension)
            .ThenBy(s => CombinatorialIndex.Index(s.Vertices))
            .ToList();

        _simplices = sorted;
        _index = new Dictionary<long, int>(sorted.Count);
        for (int i = 0; i < sorted.Count; i++)
            _index[CombinatorialIndex.PackKey(sorted[i].Dimension, sorted[i].Vertices)] = i;

        BuildCofacetLists();
    }

    private readonly Dictionary<(int Lo, int Hi), double> _edgeWeights;

    public int GetDimension(int simplexIndex) => _simplices[simplexIndex].Dimension;

    public double GetBirth(int simplexIndex) => _simplices[simplexIndex].FiltrationValue;

    public ReadOnlySpan<int> GetVertices(int simplexIndex) => _simplices[simplexIndex].Vertices;

    public int[] GetBoundaryIndices(int simplexIndex) => GetBoundaryIndicesInternal(simplexIndex);

    public int[] GetCoboundaryIndices(int simplexIndex) => _cofacets[simplexIndex];

    void BuildCofacetLists()
    {
        var builders = new List<int>[_simplices.Count];
        for (int i = 0; i < _simplices.Count; i++)
            builders[i] = new List<int>();

        for (int j = 0; j < _simplices.Count; j++)
        {
            int[] boundary = GetBoundaryIndicesInternal(j);
            for (int b = 0; b < boundary.Length; b++)
                builders[boundary[b]].Add(j);
        }

        _cofacets = new int[_simplices.Count][];
        for (int i = 0; i < _simplices.Count; i++)
        {
            List<int> cof = builders[i];
            if (cof.Count == 0)
                _cofacets[i] = Array.Empty<int>();
            else
            {
                cof.Sort();
                cof.Reverse();
                _cofacets[i] = cof.ToArray();
            }
        }
    }

    static IEnumerable<int> CommonNeighbors(CsrGraph graph, int u, int v)
    {
        int startU = graph.RowPointers[u];
        int endU = graph.RowPointers[u + 1];
        int startV = graph.RowPointers[v];
        int endV = graph.RowPointers[v + 1];

        int i = startU;
        int j = startV;
        while (i < endU && j < endV)
        {
            int nu = graph.Targets[i];
            int nv = graph.Targets[j];
            if (nu == nv)
            {
                yield return nu;
                i++;
                j++;
            }
            else if (nu < nv)
                i++;
            else
                j++;
        }
    }

    int[] GetBoundaryIndicesInternal(int j)
    {
        int[] v = _simplices[j].Vertices;
        if (v.Length <= 1) return Array.Empty<int>();
        var result = new int[v.Length];
        for (int k = 0; k < v.Length; k++)
            result[k] = _index[CombinatorialIndex.PackKey(v.Length - 2, RemoveAt(v, k))];
        return result;
    }

    static (int Lo, int Hi) Key(int u, int v) => u < v ? (u, v) : (v, u);

    static int[] RemoveAt(int[] v, int k)
    {
        var face = new int[v.Length - 1];
        for (int i = 0, w = 0; i < v.Length; i++)
            if (i != k) face[w++] = v[i];
        return face;
    }

    static Dictionary<(int Lo, int Hi), double> ResolveEdgeWeights(CsrGraph g, FiltrationWeights weights) =>
        weights switch
        {
            RawDistanceWeights => BuildRawDistanceWeights(g),
            EffectiveResistanceWeights er => EffectiveResistance.BuildEdgeWeights(g, er.Eigenpairs),
            _ => throw new ArgumentException($"Unknown {nameof(FiltrationWeights)}: {weights.GetType().Name}", nameof(weights)),
        };

    static Dictionary<(int Lo, int Hi), double> BuildRawDistanceWeights(CsrGraph g)
    {
        var map = new Dictionary<(int Lo, int Hi), double>();
        for (int u = 0; u < g.NodeCount; u++)
        {
            int rowEnd = g.RowPointers[u + 1];
            for (int e = g.RowPointers[u]; e < rowEnd; e++)
            {
                int v = g.Targets[e];
                if (v <= u) continue;

                double w = g.Weights[e];
                if (w < 0.0)
                    throw new ArgumentException("Rips filtration requires non-negative edge weights (distances).");

                map[(u, v)] = w;
            }
        }

        return map;
    }
}
